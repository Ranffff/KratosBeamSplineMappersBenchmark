#!/usr/bin/env python3
"""Validate BeamSplineMapper forward mapping against analytical rigid rotation."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


CASE_ROOT = Path(__file__).resolve().parent.parent
SOURCE = (
    CASE_ROOT.parent
    / "Bending_Torsion_Beam_Test_Case_2nd"
    / "analytical_finite_rotation_levels_v2.py"
)
# Canonical result directory. Each run replaces only this exact directory.
OUTPUT_DIRECTORY = (
    CASE_ROOT
    / "TestCase_Output"
    / "MapperVerification"
    / "BeamSplineMapper_Forward"
)


def main() -> None:
    output_directory = OUTPUT_DIRECTORY.resolve()
    expected_parent = (CASE_ROOT / "TestCase_Output" / "MapperVerification").resolve()
    if output_directory.parent != expected_parent:
        raise RuntimeError(f"Unsafe configured output directory: {output_directory}")
    if output_directory.exists():
        shutil.rmtree(output_directory)
    subprocess.run(
        [
            sys.executable,
            str(SOURCE),
            "--mapper",
            "plain",
            "--output-dir",
            str(output_directory),
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
