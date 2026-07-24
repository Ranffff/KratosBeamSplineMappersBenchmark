#!/usr/bin/env python3
"""Validate BeamSplineMapper's forward tangent and adjoint virtual work."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path


CASE_ROOT = Path(__file__).resolve().parent.parent
KRATOS_ROOT = CASE_ROOT.parent
SOURCE = (
    KRATOS_ROOT
    / "Bending_Torsion_Beam_Test_Case_2nd"
    / "analytical_plain_tangent_work_v2.py"
)
DEFAULT_OUTPUT = CASE_ROOT / "TestCase_Output" / "MapperVerification" / "BeamSplineMapper_Adjoint"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--mesh-preset", choices=("small", "medium", "large"), default="small")
    parser.add_argument("--case-limit", type=int)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite {args.output_dir}")

    spec = importlib.util.spec_from_file_location("plain_tangent_work", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.OUTPUT_ROOT = args.output_dir.resolve()

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
