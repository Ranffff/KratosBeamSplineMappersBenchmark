#!/usr/bin/env python3
"""Validate recovery-small/finite tangent transpose by discrete virtual work."""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import sys
from pathlib import Path

import KratosMultiphysics as KM


CASE_ROOT = Path(__file__).resolve().parent.parent
SOURCE = (
    CASE_ROOT.parent
    / "Bending_Torsion_Beam_Test_Case_2nd"
    / "analytical_finite_tangent_work_v2.py"
)
# Canonical result directories. Each selected mode replaces only its own directory.
OUTPUT_ROOT = CASE_ROOT / "TestCase_Output" / "MapperVerification"
OUTPUT_DIRECTORIES = {
    "small": OUTPUT_ROOT / "Recovery_Small_Adjoint",
    "finite": OUTPUT_ROOT / "Recovery_Finite_Adjoint",
}


def create_recovery_mapper(origin, destination, mode: str):
    settings = KM.Parameters(f"""{{
        "mapper_type" : "beam_spline_mapper_with_recovery_of_rotations",
        "search_settings" : {{
            "search_radius" : 3.0,
            "max_num_search_iterations" : 30
        }},
        "local_coord_tolerance" : 0.25,
        "kernel_type" : "gaussian",
        "kernel_radius" : 0.50,
        "polynomial_level" : 0,
        "rotation_recovery_mode" : "{mode}",
        "regularization" : 1.0e-8,
        "echo_level" : 0
    }}""")
    return KM.MapperFactory.CreateMapper(origin, destination, settings)


def run_mode(mode: str, mesh_preset: str, case_limit: int | None) -> None:
    output_dir = OUTPUT_DIRECTORIES[mode].resolve()
    if output_dir.parent != OUTPUT_ROOT.resolve():
        raise RuntimeError(f"Unsafe configured output directory: {output_dir}")
    if output_dir.exists():
        shutil.rmtree(output_dir)
    spec = importlib.util.spec_from_file_location(f"recovery_{mode}_tangent_work", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.OUTPUT_ROOT = output_dir
    module.create_beam_spline_mapper = (
        lambda origin, destination: create_recovery_mapper(origin, destination, mode)
    )

    forwarded = [str(SOURCE), "--mesh-preset", mesh_preset]
    if case_limit is not None:
        forwarded += ["--case-limit", str(case_limit)]
    original_argv = sys.argv
    try:
        sys.argv = forwarded
        module.main()
    finally:
        sys.argv = original_argv


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("small", "finite", "both"), default="both")
    parser.add_argument("--mesh-preset", choices=("small", "medium", "large"), default="small")
    parser.add_argument("--case-limit", type=int)
    args = parser.parse_args()
    modes = ("small", "finite") if args.mode == "both" else (args.mode,)
    for mode in modes:
        run_mode(mode, args.mesh_preset, args.case_limit)


if __name__ == "__main__":
    main()
