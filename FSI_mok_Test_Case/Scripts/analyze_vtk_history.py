#!/usr/bin/env python3
"""Extract reproducible stability diagnostics from Kratos legacy ASCII VTK files."""

import argparse
import csv
import json
import math
import re
import statistics
from pathlib import Path


STEP_PATTERN = re.compile(r"_(\d+)\.vtk$")
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
        if len(parts) == 4 and parts[0] in VECTOR_FIELDS:
            name = parts[0]
            components = int(parts[1])
            count = int(parts[2])
            values = []
            index += 1
            while len(values) < components * count:
                values.extend(float(value) for value in lines[index].split())
                index += 1
            fields[name] = [
                values[components * i:components * i + components]
                for i in range(count)
            ]
            continue
        index += 1
    return points, fields


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
        match = STEP_PATTERN.search(path.name)
        if match:
            files.append((int(match.group(1)), path))
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
    summary_path.write_text(json.dumps(summarize(rows), indent=2) + "\n")
    print(f"Wrote {len(rows)} steps to {args.output}")
    print(f"Wrote summary to {summary_path}")


if __name__ == "__main__":
    main()
