#!/usr/bin/env python3
"""Validate recovery-small/finite forward maps against analytical rigid rotation."""

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
DEFAULT_OUTPUT = CASE_ROOT / "TestCase_Output" / "MapperVerification"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("small", "finite", "both"), default="both")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    modes = ("small", "finite") if args.mode == "both" else (args.mode,)
    for mode in modes:
        output_dir = args.output_root.resolve() / f"Recovery_{mode.capitalize()}_Forward"
        if output_dir.exists():
            raise FileExistsError(f"Refusing to overwrite {output_dir}")
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
