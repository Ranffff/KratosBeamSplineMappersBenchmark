#!/usr/bin/env python3
"""Versioned analytical benchmark for small/finite beam-spline recovery.

This file is intentionally separate from every historical benchmark.  Each
mapper result is compared only with the prescribed rigid-body displacement.
"""

from __future__ import annotations

import csv
import argparse
import json
import math
from pathlib import Path

import KratosMultiphysics as KM
import KratosMultiphysics.LinearSolversApplication  # noqa: F401
import KratosMultiphysics.MappingApplication  # noqa: F401
import KratosMultiphysics.StructuralMechanicsApplication  # noqa: F401


ROOT = Path(__file__).resolve().parent
OUTPUT_ROOT = ROOT / "TestCase_Output" / "Analytical_Finite_Rotation_Levels_v2"
LENGTH = 1.0
WIDTH = 0.2
BEAM_ELEMENTS = 32
XI_DIVISIONS = 16
ETA_DIVISIONS = 80
ANGLES = tuple(math.radians(value) for value in (1.0, 30.0, 90.0, 170.0))


def array3(values):
    return KM.Array3([float(value) for value in values])


def finite_rotation_displacement(x: float, y: float, theta: float):
    cosine = math.cos(theta)
    sine = math.sin(theta)
    return (cosine * x - sine * y - x, sine * x + cosine * y - y, 0.0)


def create_beam(theta: float):
    model = KM.Model()
    beam = model.CreateModelPart("beam")
    beam.ProcessInfo[KM.DOMAIN_SIZE] = 3
    beam.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    beam.AddNodalSolutionStepVariable(KM.ROTATION)
    properties = beam.CreateNewProperties(1)
    for i in range(BEAM_ELEMENTS + 1):
        y = LENGTH * i / BEAM_ELEMENTS
        node = beam.CreateNewNode(i + 1, 0.0, y, 0.0)
        node.SetSolutionStepValue(KM.DISPLACEMENT, array3(finite_rotation_displacement(0.0, y, theta)))
        node.SetSolutionStepValue(KM.ROTATION, array3((0.0, 0.0, theta)))
        # A mapper is not allowed to overwrite current coordinates during
        # reference-configuration projection/search.
        node.X = node.X0 + 0.0123
        node.Y = node.Y0 - 0.0045
        node.Z = node.Z0 + 0.0021
    for i in range(BEAM_ELEMENTS):
        beam.CreateNewElement("CrBeamElement3D2N", i + 1, [i + 1, i + 2], properties)
    return model, beam


def create_surface(model):
    surface = model.CreateModelPart("surface")
    surface.ProcessInfo[KM.DOMAIN_SIZE] = 3
    surface.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    properties = surface.CreateNewProperties(1)
    for j in range(ETA_DIVISIONS + 1):
        y = LENGTH * j / ETA_DIVISIONS
        for i in range(XI_DIVISIONS + 1):
            x = -0.5 * WIDTH + WIDTH * i / XI_DIVISIONS
            node_id = j * (XI_DIVISIONS + 1) + i + 1
            surface.CreateNewNode(node_id, x, y, 0.0)
    element_id = 1
    for j in range(ETA_DIVISIONS):
        row = j * (XI_DIVISIONS + 1)
        next_row = (j + 1) * (XI_DIVISIONS + 1)
        for i in range(XI_DIVISIONS):
            surface.CreateNewElement(
                "ShellThinElementCorotational3D4N",
                element_id,
                [row + i + 1, row + i + 2, next_row + i + 2, next_row + i + 1],
                properties,
            )
            element_id += 1
    return surface


def run_angle(theta: float, mapper_name: str):
    model, beam = create_beam(theta)
    surface = create_surface(model)
    coordinates_before = {
        node.Id: (node.X, node.Y, node.Z) for node in beam.Nodes
    }
    if mapper_name in ("recovery_small", "recovery_finite"):
        mode = "small" if mapper_name == "recovery_small" else "finite"
        settings = KM.Parameters(r'''{
            "mapper_type": "beam_spline_mapper_with_recovery_of_rotations",
            "search_settings": {"search_radius": 0.30, "max_num_search_iterations": 30},
            "local_coord_tolerance": 0.25,
            "kernel_type": "gaussian",
            "kernel_radius": 0.50,
            "polynomial_level": 4,
            "rotation_recovery_mode": "small",
            "regularization": 1.0e-8,
            "echo_level": 0
        }''')
        settings["rotation_recovery_mode"].SetString(mode)
    else:
        settings = KM.Parameters(r'''{
            "mapper_type": "beam_spline_mapper",
            "search_settings": {"search_radius": 0.30, "max_num_search_iterations": 30},
            "local_coord_tolerance": 0.25,
            "echo_level": 0
        }''')
    mapper = KM.MapperFactory.CreateMapper(beam, surface, settings)
    coordinate_mutation = max(
        math.sqrt(
            (node.X - coordinates_before[node.Id][0]) ** 2
            + (node.Y - coordinates_before[node.Id][1]) ** 2
            + (node.Z - coordinates_before[node.Id][2]) ** 2
        )
        for node in beam.Nodes
    )
    mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT)
    error_squared = 0.0
    reference_squared = 0.0
    max_error = 0.0
    for node in surface.Nodes:
        mapped = node.GetSolutionStepValue(KM.DISPLACEMENT)
        exact = finite_rotation_displacement(node.X0, node.Y0, theta)
        error = math.sqrt(sum((float(mapped[i]) - exact[i]) ** 2 for i in range(3)))
        reference_norm_squared = sum(value * value for value in exact)
        error_squared += error * error
        reference_squared += reference_norm_squared
        max_error = max(max_error, error)
    return {
        "theta_rad": theta,
        "sin_theta": math.sin(theta),
        "theta_minus_sin_theta": theta - math.sin(theta),
        "relative_l2_error": math.sqrt(error_squared / reference_squared),
        "max_displacement_error": max_error,
        "surface_nodes": surface.NumberOfNodes(),
        "coordinate_mutation": coordinate_mutation,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mapper",
        choices=("recovery_small", "recovery_finite", "plain"),
        default="recovery_finite",
    )
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir or (
        OUTPUT_ROOT / args.mapper
    )
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing diagnostic: {output_dir}")
    output_dir.mkdir(parents=True)
    results = [run_angle(theta, args.mapper) for theta in ANGLES]
    (output_dir / "summary.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    with (output_dir / "summary.csv").open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    (output_dir / "parameters.json").write_text(
        json.dumps(
            {
                "beam_elements": BEAM_ELEMENTS,
                "surface_divisions": {"xi": XI_DIVISIONS, "eta": ETA_DIVISIONS},
                "mapper": args.mapper,
                "reference": "prescribed analytical Rz rigid-body displacement; no mapper-to-mapper reference",
                "relative_l2_definition": "sqrt(sum(||u_map-u_exact||^2)/sum(||u_exact||^2))",
                "exact_span_tolerance": 1.0e-10,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    for result in results:
        print(
            f"theta={result['theta_rad']:.2f}: rel_L2={result['relative_l2_error']:.6e}, "
            f"max_error={result['max_displacement_error']:.6e}"
        )


if __name__ == "__main__":
    main()
