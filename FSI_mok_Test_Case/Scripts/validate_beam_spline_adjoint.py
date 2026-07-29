#!/usr/bin/env python3
"""Validate BeamSplineMapper's forward tangent and adjoint virtual work."""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import sys
from pathlib import Path


CASE_ROOT = Path(__file__).resolve().parent.parent
KRATOS_ROOT = CASE_ROOT.parent
SOURCE = (
    KRATOS_ROOT
    / "Bending_Torsion_Beam_Test_Case_2nd"
    / "analytical_plain_tangent_work_v2.py"
)
# Canonical result directory. Each run replaces only this exact directory.
OUTPUT_DIRECTORY = (
    CASE_ROOT
    / "TestCase_Output"
    / "MapperVerification"
    / "BeamSplineMapper_Adjoint"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-preset", choices=("small", "medium", "large"), default="small")
    parser.add_argument("--case-limit", type=int)
    args = parser.parse_args()
    output_directory = OUTPUT_DIRECTORY.resolve()
    expected_parent = (CASE_ROOT / "TestCase_Output" / "MapperVerification").resolve()
    if output_directory.parent != expected_parent:
        raise RuntimeError(f"Unsafe configured output directory: {output_directory}")
    if output_directory.exists():
        shutil.rmtree(output_directory)

    spec = importlib.util.spec_from_file_location("plain_tangent_work", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.OUTPUT_ROOT = output_directory

    forwarded = [str(SOURCE), "--mesh-preset", args.mesh_preset]
    if args.case_limit is not None:
        forwarded += ["--case-limit", str(args.case_limit)]
    original_argv = sys.argv
    try:
        sys.argv = forwarded
        module.main()
    finally:
        sys.argv = original_argv


if __name__ == "__main__":
    main()
