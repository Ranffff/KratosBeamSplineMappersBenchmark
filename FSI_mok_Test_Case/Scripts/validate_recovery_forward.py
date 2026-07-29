#!/usr/bin/env python3
"""Validate recovery-small/finite forward maps against analytical rigid rotation."""

from __future__ import annotations

import argparse
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
# Canonical result directories. Each selected mode replaces only its own directory.
OUTPUT_ROOT = CASE_ROOT / "TestCase_Output" / "MapperVerification"
OUTPUT_DIRECTORIES = {
    "small": OUTPUT_ROOT / "Recovery_Small_Forward",
    "finite": OUTPUT_ROOT / "Recovery_Finite_Forward",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("small", "finite", "both"), default="both")
    args = parser.parse_args()
    modes = ("small", "finite") if args.mode == "both" else (args.mode,)
    for mode in modes:
        output_dir = OUTPUT_DIRECTORIES[mode].resolve()
        if output_dir.parent != OUTPUT_ROOT.resolve():
            raise RuntimeError(f"Unsafe configured output directory: {output_dir}")
        if output_dir.exists():
            shutil.rmtree(output_dir)
        subprocess.run(
            [
                sys.executable,
                str(SOURCE),
                "--mapper",
                f"recovery_{mode}",
                "--output-dir",
                str(output_dir),
            ],
            check=True,
        )


if __name__ == "__main__":
    main()
