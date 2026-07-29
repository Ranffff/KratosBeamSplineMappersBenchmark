#!/usr/bin/env python3
"""Reduce echo-level structural Newton histories by FSI step/coupling iteration."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


STEP_RE = re.compile(r"Fluid Dynamics Analysis: STEP:\s+(\d+)")
TIME_RE = re.compile(r"Fluid Dynamics Analysis: TIME:\s+([-+0-9.eE]+)")
COUPLING_RE = re.compile(r"Coupling iteration:.*?(\d+)\s*/\s*(\d+)")
DISP_RE = re.compile(r"RESIDUAL DISPLACEMENT CRITERION :: Ratio =\s*([-+0-9.eE]+);\s*Norm =\s*([-+0-9.eE]+)")
ROT_RE = re.compile(r"RESIDUAL ROTATION CRITERION :: Ratio =\s*([-+0-9.eE]+);\s*Norm =\s*([-+0-9.eE]+)")
CONVERGED_RE = re.compile(r"Convergence achieved after\s+(\d+)\s*/\s*(\d+) iterations")
NONCONVERGED_RE = re.compile(r"Solver did not converge for step\s+(\d+)")


def new_record(step, time, coupling_iteration, coupling_limit, line_number):
    return {
        "step": step,
        "time": time,
        "coupling_iteration": coupling_iteration,
        "coupling_limit": coupling_limit,
        "start_line": line_number,
        "end_line": None,
        "newton_iterations_reported": 0,
        "newton_limit": None,
        "converged": None,
        "displacement_ratio_first": None,
        "displacement_ratio_last": None,
        "displacement_ratio_max": None,
        "displacement_norm_last": None,
        "rotation_ratio_first": None,
        "rotation_ratio_last": None,
        "rotation_ratio_max": None,
        "rotation_norm_last": None,
    }


def finalize(records, current, line_number):
    if current is not None:
        current["end_line"] = line_number
        records.append(current)
    return None


def analyze(log_path: Path):
    records = []
    current = None
    step = None
    time = None
    lines = log_path.read_text(errors="replace").splitlines()
    for line_number, line in enumerate(lines, 1):
        step_match = STEP_RE.search(line)
        if step_match:
            step = int(step_match.group(1))
        time_match = TIME_RE.search(line)
        if time_match:
            time = float(time_match.group(1))
        coupling_match = COUPLING_RE.search(line)
        if coupling_match:
            current = finalize(records, current, line_number - 1)
            current = new_record(
                step,
                time,
                int(coupling_match.group(1)),
                int(coupling_match.group(2)),
                line_number,
            )
            continue
        if current is None:
            continue
        for pattern, prefix in ((DISP_RE, "displacement"), (ROT_RE, "rotation")):
            match = pattern.search(line)
            if match:
                ratio = float(match.group(1))
                value_norm = float(match.group(2))
                first_key = f"{prefix}_ratio_first"
                current[first_key] = ratio if current[first_key] is None else current[first_key]
                current[f"{prefix}_ratio_last"] = ratio
                max_key = f"{prefix}_ratio_max"
                current[max_key] = ratio if current[max_key] is None else max(current[max_key], ratio)
                current[f"{prefix}_norm_last"] = value_norm
        converged_match = CONVERGED_RE.search(line)
        if converged_match:
            current["newton_iterations_reported"] = int(converged_match.group(1))
            current["newton_limit"] = int(converged_match.group(2))
            current["converged"] = 1
        if NONCONVERGED_RE.search(line):
            current["newton_iterations_reported"] = max(
                current["newton_iterations_reported"],
                150,
            )
            current["newton_limit"] = 150
            current["converged"] = 0
    finalize(records, current, len(lines))
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("log", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix(".json").exists():
        raise FileExistsError(f"Refusing to overwrite Newton diagnostic: {args.output}")
    records = analyze(args.log)
    if not records:
        raise RuntimeError("No coupling/Newton histories found")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    failed = [record for record in records if record["converged"] == 0]
    completed = [record for record in records if record["converged"] == 1]
    summary = {
        "log": str(args.log.resolve()),
        "number_of_coupling_solves": len(records),
        "number_of_converged_structural_solves": len(completed),
        "number_of_failed_structural_solves": len(failed),
        "first_failed_structural_solve": failed[0] if failed else None,
        "maximum_converged_newton_iterations": max(
            completed,
            key=lambda record: record["newton_iterations_reported"],
        ) if completed else None,
    }
    args.output.with_suffix(".json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
