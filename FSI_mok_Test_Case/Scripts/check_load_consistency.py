#!/usr/bin/env python3
"""Check force and moment conservation between fluid reactions and beam loads."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

from analyze_vtk_history import cross, read_legacy_vtk


STEP_PATTERN = re.compile(r"_(\d+)\.vtk$")


def read_submodelpart_node_ids(mdpa: Path, name: str) -> list[int]:
    ids = []
    in_part = False
    in_nodes = False
    for line in mdpa.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped == f"Begin SubModelPart {name}":
            in_part = True
        elif in_part and stripped == "Begin SubModelPartNodes":
            in_nodes = True
        elif in_nodes and stripped == "End SubModelPartNodes":
            break
        elif in_nodes and stripped:
            ids.append(int(stripped))
    if not ids:
        raise RuntimeError(f"No nodes found for SubModelPart {name} in {mdpa}")
    return ids


def numbered_files(directory: Path) -> dict[int, Path]:
    result = {}
    for path in directory.glob("*.vtk"):
        match = STEP_PATTERN.search(path.name)
        if match:
            result[int(match.group(1))] = path
    return result


def sum_vectors(values):
    return [sum(value[i] for value in values) for i in range(3)]


def total_moment(points, forces, moments, root):
    result = sum_vectors(moments) if moments else [0.0, 0.0, 0.0]
    for point, force in zip(points, forces):
        lever = [point[i] - root[i] for i in range(3)]
        contribution = cross(lever, force)
        for i in range(3):
            result[i] += contribution[i]
    return result


def norm(values):
    return math.sqrt(sum(value * value for value in values))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fluid_vtk_dir", type=Path)
    parser.add_argument("structure_vtk_dir", type=Path)
    parser.add_argument("--fluid-mdpa", type=Path, required=True)
    parser.add_argument("--scale", type=float, required=True)
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--root", type=float, nargs=3, default=(0.499, 0.0, 0.0))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.output.exists() or args.output.with_suffix(".json").exists():
        raise FileExistsError(f"Refusing to overwrite load-consistency output: {args.output}")
    interface_ids = read_submodelpart_node_ids(args.fluid_mdpa, "interface")
    fluid_files = numbered_files(args.fluid_vtk_dir)
    structure_files = numbered_files(args.structure_vtk_dir)
    common_steps = sorted(set(fluid_files) & set(structure_files))
    if not common_steps:
        raise RuntimeError("No common fluid/structure VTK steps")

    rows = []
    for step in common_steps:
        fluid_points, fluid_fields = read_legacy_vtk(fluid_files[step])
        structure_points, structure_fields = read_legacy_vtk(structure_files[step])
        if "REACTION" not in fluid_fields:
            raise RuntimeError(f"Fluid REACTION is absent at step {step}")
        interface_indices = [node_id - 1 for node_id in interface_ids]
        fluid_reactions = [fluid_fields["REACTION"][index] for index in interface_indices]
        fluid_interface_points = [fluid_points[index] for index in interface_indices]
        expected_force = [-args.scale * value for value in sum_vectors(fluid_reactions)]
        mapped_force = sum_vectors(structure_fields["POINT_LOAD"])
        fluid_moment = total_moment(fluid_interface_points, fluid_reactions, None, args.root)
        expected_moment = [-args.scale * value for value in fluid_moment]
        mapped_moment = total_moment(
            structure_points,
            structure_fields["POINT_LOAD"],
            structure_fields.get("POINT_MOMENT"),
            args.root,
        )
        force_error = [mapped_force[i] - expected_force[i] for i in range(3)]
        moment_error = [mapped_moment[i] - expected_moment[i] for i in range(3)]
        row = {"step": step, "time": step * args.dt}
        for label, values in (
            ("expected_force", expected_force),
            ("mapped_force", mapped_force),
            ("force_error", force_error),
            ("expected_moment", expected_moment),
            ("mapped_moment", mapped_moment),
            ("moment_error", moment_error),
        ):
            for component, component_label in enumerate("xyz"):
                row[f"{label}_{component_label}"] = values[component]
        row["force_error_norm"] = norm(force_error)
        row["force_relative_error"] = norm(force_error) / max(norm(expected_force), 1.0e-30)
        row["moment_error_norm"] = norm(moment_error)
        row["moment_relative_error"] = norm(moment_error) / max(norm(expected_moment), 1.0e-30)
        rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "fluid_vtk_dir": str(args.fluid_vtk_dir.resolve()),
        "structure_vtk_dir": str(args.structure_vtk_dir.resolve()),
        "fluid_mdpa": str(args.fluid_mdpa.resolve()),
        "interface_node_count": len(interface_ids),
        "scale": args.scale,
        "sign_convention": "expected beam load = -scale * fluid interface REACTION (swap_sign)",
        "first_step": common_steps[0],
        "last_step": common_steps[-1],
        "number_of_common_steps": len(common_steps),
        "max_force_relative_error": max(row["force_relative_error"] for row in rows),
        "max_moment_relative_error": max(row["moment_relative_error"] for row in rows),
        "max_force_error_norm": max(row["force_error_norm"] for row in rows),
        "max_moment_error_norm": max(row["moment_error_norm"] for row in rows),
        "note": "Legacy VTK output preserves ascending node-ID order for Parts_Fluid; interface node ID maps to zero-based index ID-1.",
    }
    args.output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
