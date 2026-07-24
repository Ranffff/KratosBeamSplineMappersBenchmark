#!/usr/bin/env python3
"""Validate BeamSplineMapper forward mapping against analytical rigid rotation."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


CASE_ROOT = Path(__file__).resolve().parent.parent
SOURCE = (
    CASE_ROOT.parent
    / "Bending_Torsion_Beam_Test_Case_2nd"
    / "analytical_finite_rotation_levels_v2.py"
)
DEFAULT_OUTPUT = (
    CASE_ROOT
    / "TestCase_Output"
    / "MapperVerification"
    / "BeamSplineMapper_Forward"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite {args.output_dir}")
    subprocess.run(
        [
            sys.executable,
            str(SOURCE),
            "--mapper",
            "plain",
            "--output-dir",
            str(args.output_dir.resolve()),
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
