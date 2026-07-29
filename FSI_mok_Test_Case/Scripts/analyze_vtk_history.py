#!/usr/bin/env python3
"""Extract reproducible stability diagnostics from Kratos legacy ASCII VTK files."""

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
from pathlib import Path


STEP_PATTERN = re.compile(r"^(?P<series>.+)_(?P<step>\d+)\.vtk$")
VECTOR_FIELDS = (
    "DISPLACEMENT",
    "ROTATION",
    "POINT_LOAD",
    "POINT_MOMENT",
    "REACTION",
    "REACTION_MOMENT",
)


def read_legacy_vtk(path):
    lines = path.read_text().splitlines()
    if len(lines) < 5 or lines[2].strip().upper() != "ASCII":
        raise RuntimeError(f"{path} is not a legacy ASCII VTK file")
    if not any("DATASET UNSTRUCTURED_GRID" in line for line in lines[:8]):
        raise RuntimeError(f"{path} is not an UNSTRUCTURED_GRID VTK file")
    points = []
    fields = {}
    index = 0
    while index < len(lines):
        parts = lines[index].split()
        if parts and parts[0] == "POINTS":
            count = int(parts[1])
            values = []
            index += 1
            while len(values) < 3 * count:
                values.extend(float(value) for value in lines[index].split())
                index += 1
            points = [values[3 * i:3 * i + 3] for i in range(count)]
            continue
        if len(parts) >= 3 and parts[0] == "FIELD":
            number_of_arrays = int(parts[2])
            index += 1
            for _ in range(number_of_arrays):
                if index >= len(lines):
                    raise RuntimeError(f"Truncated FIELD section in {path}")
                header = lines[index].split()
                if len(header) != 4:
                    raise RuntimeError(
                        f"Invalid FIELD array header in {path}: {lines[index]}"
                    )
                name = header[0]
                components = int(header[1])
                count = int(header[2])
                values = []
                index += 1
                while len(values) < components * count:
                    if index >= len(lines):
                        raise RuntimeError(
                            f"Truncated {name} array in {path}"
                        )
                    values.extend(float(value) for value in lines[index].split())
                    index += 1
                if len(values) != components * count:
                    raise RuntimeError(
                        f"Wrong value count for {name} in {path}"
                    )
                fields[name] = [
                    values[components * i:components * i + components]
                    for i in range(count)
                ]
            continue
        index += 1
    if not points:
        raise RuntimeError(f"No POINTS section found in {path}")
    return points, fields


def read_point_history(path):
    values = {}
    if not path or not Path(path).exists():
        return values
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            time = float(parts[0])
            value = float(parts[1])
        except ValueError:
            continue
        if math.isfinite(time) and math.isfinite(value):
            values[round(time, 10)] = value
    return values


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_vtk_directory(
    vtk_dir,
    dt,
    required_fields=(),
    max_valid_time=None,
    point_history=None,
    point_position=None,
    point_component=0,
    point_rel_l2_tolerance=1.0e-6,
):
    """Return a machine-readable integrity and point-history manifest."""
    vtk_dir = Path(vtk_dir)
    grouped = {}
    for path in sorted(vtk_dir.glob("*.vtk")):
        match = STEP_PATTERN.match(path.name)
        if match:
            grouped.setdefault(match.group("series"), []).append(
                (int(match.group("step")), path)
            )
    if not grouped:
        raise RuntimeError(f"No numbered VTK files found in {vtk_dir}")

    manifest = {
        "vtk_dir": str(vtk_dir.resolve()),
        "dt": dt,
        "required_fields": list(required_fields),
        "series": {},
        "errors": [],
    }
    parsed = {}
    for series_name, files in sorted(grouped.items()):
        files.sort()
        series_record = {
            "number_of_files": len(files),
            "steps": [step for step, _ in files],
            "files": [],
        }
        expected_point_count = None
        for step, path in files:
            record = {
                "file": path.name,
                "step": step,
                "time": step * dt,
                "bytes": path.stat().st_size,
                "sha256": file_sha256(path),
            }
            try:
                points, fields = read_legacy_vtk(path)
                record["point_count"] = len(points)
                record["fields"] = sorted(fields)
                record["field_size_mismatches"] = sorted(
                    name
                    for name, values in fields.items()
                    if len(values) != len(points)
                )
                record["all_values_finite"] = all(
                    math.isfinite(value)
                    for point in points
                    for value in point
                ) and all(
                    math.isfinite(value)
                    for values in fields.values()
                    for row in values
                    for value in row
                )
                record["missing_required_fields"] = sorted(
                    set(required_fields) - set(fields)
                )
                if expected_point_count is None:
                    expected_point_count = len(points)
                elif len(points) != expected_point_count:
                    record["point_count_changed"] = True
                if (
                    max_valid_time is None
                    or record["time"] <= max_valid_time + 1.0e-10
                ):
                    record["trusted_complete_output"] = True
                    parsed[(series_name, step)] = (points, fields)
                else:
                    record["trusted_complete_output"] = False
            except Exception as error:
                record["parse_error"] = f"{type(error).__name__}: {error}"
            if (
                "parse_error" in record
                or not record.get("all_values_finite", False)
                or record.get("missing_required_fields")
                or record.get("field_size_mismatches")
                or record.get("point_count_changed", False)
            ):
                manifest["errors"].append(
                    f"{series_name}/{path.name} failed integrity checks"
                )
            series_record["files"].append(record)

        trusted = [
            record
            for record in series_record["files"]
            if record.get("trusted_complete_output")
            and "parse_error" not in record
        ]
        if trusted:
            series_record["last_trusted_step"] = trusted[-1]["step"]
            series_record["last_trusted_time"] = trusted[-1]["time"]
        manifest["series"][series_name] = series_record

    if max_valid_time is not None:
        trusted_times = [
            record["time"]
            for series in manifest["series"].values()
            for record in series["files"]
            if record.get("trusted_complete_output")
            and "parse_error" not in record
        ]
        manifest["max_valid_time"] = max_valid_time
        manifest["last_vtk_time"] = max(trusted_times) if trusted_times else None
        if not trusted_times:
            manifest["errors"].append("No complete VTK at or before last_valid_time")
        elif max_valid_time - max(trusted_times) > 0.51 * dt:
            manifest["errors"].append(
                "Last complete VTK is more than half a time step behind "
                "last_valid_time"
            )

    if point_history and point_position is not None:
        history = read_point_history(point_history)
        candidates = []
        for (series_name, step), (points, fields) in parsed.items():
            if "DISPLACEMENT" in fields:
                candidates.append((series_name, step, points, fields))
        if not candidates:
            manifest["errors"].append(
                "No trusted DISPLACEMENT VTK available for point-history cross-check"
            )
        else:
            candidates.sort(key=lambda item: (item[0], item[1]))
            primary_series = candidates[0][0]
            first = next(item for item in candidates if item[0] == primary_series)
            point_index = min(
                range(len(first[2])),
                key=lambda i: sum(
                    (first[2][i][j] - point_position[j]) ** 2
                    for j in range(3)
                ),
            )
            point_distance = math.sqrt(
                sum(
                    (first[2][point_index][j] - point_position[j]) ** 2
                    for j in range(3)
                )
            )
            vtk_values = {}
            for series_name, step, _, fields in candidates:
                if series_name != primary_series:
                    continue
                time = round(step * dt, 10)
                vtk_values[time] = fields["DISPLACEMENT"][point_index][
                    point_component
                ]
            common_times = sorted(set(history) & set(vtk_values))
            if max_valid_time is not None:
                common_times = [
                    time
                    for time in common_times
                    if time <= max_valid_time + 1.0e-10
                ]
            denominator = sum(history[time] ** 2 for time in common_times)
            if not common_times or denominator <= 0.0:
                manifest["errors"].append(
                    "No nonzero common samples for VTK/point-history cross-check"
                )
                rel_l2 = None
            else:
                numerator = sum(
                    (vtk_values[time] - history[time]) ** 2
                    for time in common_times
                )
                rel_l2 = math.sqrt(numerator / denominator)
                if rel_l2 > point_rel_l2_tolerance:
                    manifest["errors"].append(
                        f"VTK/point-history rel_l2={rel_l2} exceeds "
                        f"{point_rel_l2_tolerance}"
                    )
            manifest["point_history_cross_check"] = {
                "history_file": str(Path(point_history).resolve()),
                "series": primary_series,
                "reference_position": list(point_position),
                "matched_point_index": point_index,
                "matched_point_distance": point_distance,
                "component": point_component,
                "common_steps": len(common_times),
                "rel_l2": rel_l2,
                "tolerance": point_rel_l2_tolerance,
            }

    manifest["valid"] = not manifest["errors"]
    return manifest


def cross(left, right):
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def vector_metrics(name, values, previous):
    result = {}
    if not values:
        return result
    result[f"{name}_max_norm"] = max(math.sqrt(sum(x * x for x in value)) for value in values)
    for component, label in enumerate("xyz"):
        component_values = [value[component] for value in values]
        result[f"{name}_min_{label}"] = min(component_values)
        result[f"{name}_max_{label}"] = max(component_values)
        result[f"{name}_sum_{label}"] = sum(component_values)
        result[f"{name}_tip_{label}"] = component_values[0]
    if previous and len(previous) == len(values):
        result[f"{name}_step_jump_max_norm"] = max(
            math.sqrt(sum((value[j] - old[j]) ** 2 for j in range(3)))
            for value, old in zip(values, previous)
        )
    else:
        result[f"{name}_step_jump_max_norm"] = 0.0
    return result


def analyze(vtk_dir, dt, root):
    files = []
    for path in vtk_dir.glob("*.vtk"):
        match = STEP_PATTERN.match(path.name)
        if match:
            files.append((int(match.group("step")), path))
    files.sort()
    if not files:
        raise RuntimeError(f"No numbered VTK files found in {vtk_dir}")

    rows = []
    previous_fields = {}
    for step, path in files:
        points, fields = read_legacy_vtk(path)
        row = {"step": step, "time": step * dt, "vtk_file": path.name}
        finite = True
        for name in VECTOR_FIELDS:
            values = fields.get(name)
            if values is None:
                continue
            finite = finite and all(math.isfinite(x) for value in values for x in value)
            row.update(vector_metrics(name, values, previous_fields.get(name)))
            previous_fields[name] = values

        loads = fields.get("POINT_LOAD")
        moments = fields.get("POINT_MOMENT")
        if loads and points and len(loads) == len(points):
            total_moment = [0.0, 0.0, 0.0]
            if moments and len(moments) == len(points):
                for value in moments:
                    for i in range(3):
                        total_moment[i] += value[i]
            for point, force in zip(points, loads):
                lever = [point[i] - root[i] for i in range(3)]
                lever_moment = cross(lever, force)
                for i in range(3):
                    total_moment[i] += lever_moment[i]
            for i, label in enumerate("xyz"):
                row[f"total_load_moment_about_root_{label}"] = total_moment[i]
        row["all_fields_finite"] = int(finite)
        rows.append(row)
    return rows


def summarize(rows):
    summary = {
        "number_of_steps": len(rows),
        "first_step": rows[0]["step"],
        "last_step": rows[-1]["step"],
        "first_time": rows[0]["time"],
        "last_time": rows[-1]["time"],
        "all_steps_finite": all(row["all_fields_finite"] for row in rows),
    }
    first_nonfinite = next((row for row in rows if not row["all_fields_finite"]), None)
    if first_nonfinite:
        summary["first_nonfinite_step"] = first_nonfinite["step"]
        summary["first_nonfinite_time"] = first_nonfinite["time"]

    metric_names = []
    for field in VECTOR_FIELDS:
        metric_names.extend((
            f"{field}_max_norm",
            f"{field}_step_jump_max_norm",
        ))
    metric_names.append("total_load_moment_about_root_z")
    for metric in metric_names:
        available = [row for row in rows if metric in row and math.isfinite(row[metric])]
        if not available:
            continue
        peak = max(available, key=lambda row: abs(row[metric]))
        values_before_peak = [abs(row[metric]) for row in available if row["step"] < peak["step"]]
        summary[metric] = {
            "peak_abs_value": abs(peak[metric]),
            "signed_value": peak[metric],
            "step": peak["step"],
            "time": peak["time"],
        }
        if values_before_peak:
            median = statistics.median(values_before_peak)
            summary[metric]["peak_over_prior_median"] = (
                abs(peak[metric]) / median if median > 0.0 else None
            )
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vtk_dir", type=Path)
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary")
    parser.add_argument("--root", type=float, nargs=3, default=(0.499, 0.0, 0.0))
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--required-fields", nargs="*", default=())
    parser.add_argument("--max-valid-time", type=float)
    parser.add_argument("--point-history", type=Path)
    parser.add_argument("--point-position", type=float, nargs=3)
    parser.add_argument("--point-component", type=int, default=0)
    parser.add_argument("--point-rel-l2-tolerance", type=float, default=1.0e-6)
    args = parser.parse_args()

    rows = analyze(args.vtk_dir, args.dt, args.root)
    field_names = []
    for row in rows:
        for name in row:
            if name not in field_names:
                field_names.append(name)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=field_names)
        writer.writeheader()
        writer.writerows(rows)
    summary_path = Path(args.summary) if args.summary else args.output.with_suffix(".json")
    summary_path.write_text(
        json.dumps(summarize(rows), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    if args.manifest:
        manifest = validate_vtk_directory(
            args.vtk_dir,
            args.dt,
            required_fields=args.required_fields,
            max_valid_time=args.max_valid_time,
            point_history=args.point_history,
            point_position=args.point_position,
            point_component=args.point_component,
            point_rel_l2_tolerance=args.point_rel_l2_tolerance,
        )
        args.manifest.write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        if not manifest["valid"]:
            raise SystemExit(1)
    print(f"Wrote {len(rows)} steps to {args.output}")
    print(f"Wrote summary to {summary_path}")


if __name__ == "__main__":
    main()
