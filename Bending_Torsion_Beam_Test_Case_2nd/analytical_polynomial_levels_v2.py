"""Analytical-only verification of numeric polynomial levels 0--4.

Historical scripts and outputs are deliberately not modified.  This driver
reuses their mesh/field construction but supplies only the recovery mapper and
compares every result with the prescribed analytical displacement field.
"""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "beam_mapper_benchmark.py"
OUTPUT = ROOT / "TestCase_Output" / "Analytical_Polynomial_Levels_v4"


def load_source():
    spec = importlib.util.spec_from_file_location("legacy_field_definitions", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def mapper_definition(level):
    return {
        "key": f"recovery_level_{level}",
        "label": f"Recovery polynomial level {level}",
        "settings": {
            "mapper_type": "beam_spline_mapper_with_recovery_of_rotations",
            "search_settings": {
                "search_radius": 0.30,
                "max_num_search_iterations": 30,
            },
            "local_coord_tolerance": 0.25,
            "kernel_type": "gaussian",
            "kernel_radius": 0.50,
            "polynomial_level": level,
            "rotation_recovery_mode": "small",
            "regularization": 0.0,
            "echo_level": 1 if level == 0 else 0,
        },
    }


def main():
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUTPUT}")
    OUTPUT.mkdir(parents=True)
    module = load_source()
    module.OUTPUT_ROOT = OUTPUT
    module.PLOT_ROOT = OUTPUT / "Plots"
    module.VTK_ROOT = OUTPUT / "VTK"
    max_curl_error = module.validate_rotation_consistency()

    results = []
    divisions = {"xi": 8, "eta": 40}
    for level in range(5):
        definition = mapper_definition(level)
        for case in module.BENCHMARKS:
            try:
                result = module.run_single_case(
                    "analytical_polynomial_levels_v4",
                    case,
                    definition,
                    32,
                    divisions,
                )
                result["status"] = "completed"
            except RuntimeError as error:
                result = {
                    "benchmark_case": case["key"],
                    "status": "expected_unisolvency_rejection" if level == 1 else "unexpected_failure",
                    "error_message": str(error),
                    "relative_l2_error": None,
                }
            result["polynomial_level"] = level
            result["rotation_recovery_mode"] = "small"
            result["reference"] = "prescribed_analytical_displacement"
            result["exact_span_tolerance"] = 1.0e-9
            result["exact_span_pass"] = (
                result["relative_l2_error"] <= 1.0e-9
                if case["key"] == "rigid_body_rotation" and result["relative_l2_error"] is not None
                else None
            )
            results.append(result)
            print(
                f"level={level} case={case['key']} "
                f"status={result['status']} rel_L2={result['relative_l2_error']}"
            )

    (OUTPUT / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    with (OUTPUT / "summary.csv").open("w", newline="") as stream:
        fieldnames = sorted({key for result in results for key in result})
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    (OUTPUT / "parameters.json").write_text(
        json.dumps(
            {
                "source_field_definitions": str(SOURCE),
                "beam_elements": 32,
                "surface_divisions": divisions,
                "regularization": 0.0,
                "max_analytical_theta_minus_half_curl": max_curl_error,
                "error": "sqrt(sum(||u_map-u_exact||^2)/sum(||u_exact||^2))",
                "reference_policy": "analytical only; no other mapper is a reference",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
