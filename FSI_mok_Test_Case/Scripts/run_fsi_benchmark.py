#!/usr/bin/env python3
import argparse
import hashlib
import json
import math
import os
import re
import signal
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
POINT_FILE_NAMES = {
    "NearestNeighbor": "point_A_thickness50_scaled_load",
    "BeamMapper_CoRotation": "point_A_beam_corotation_t50_scaled_load",
    "BeamSplineMapper_WithRotationalRecovery": "point_A_beam_spline_recovery_t50_scaled_load",
    "BeamSplineMapper": "point_A_beam_spline_t50_scaled_load",
}

CASE_DIRS = {
    "NearestNeighbor": "CoSimulation_Cases/NearestNeighbor",
    "BeamMapper_CoRotation": "CoSimulation_Cases/BeamMapper_CoRotation",
    "BeamSplineMapper_WithRotationalRecovery": "CoSimulation_Cases/BeamSplineMapper_WithRotationalRecovery",
    "BeamSplineMapper": "CoSimulation_Cases/BeamSplineMapper",
}

BEAM_SECTION_AT_SCALE_50 = {
    "CROSS_AREA": 0.25,
    "I22": 0.020833333333333336,
    "I33": 5.208333333333334e-7,
    "TORSIONAL_INERTIA": 2.0833333333333334e-6,
}


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, data):
    path.write_text(json.dumps(data, indent=4) + "\n")


def set_end_time(data, end_time):
    data["problem_data"]["end_time"] = end_time


def set_time_step(data, dt):
    settings = data["solver_settings"]
    if "time_stepping" in settings:
        settings["time_stepping"]["time_step"] = dt
    if "fluid_solver_settings" in settings:
        settings["fluid_solver_settings"]["time_stepping"]["time_step"] = dt


def set_coupling(data, alpha, iterations):
    settings = data["solver_settings"]
    settings["num_coupling_iterations"] = iterations
    accelerators = settings.setdefault("convergence_accelerators", [])
    if not accelerators:
        accelerators.append({})
    accelerators[0].clear()
    accelerators[0].update({
        "type": "constant_relaxation",
        "solver": "fluid",
        "data_name": "disp",
        "alpha": alpha,
    })
    fluid_data = settings["solvers"]["fluid"]["data"]
    fluid_data.setdefault("velocity", {
        "model_part_name": "FluidModelPart.interface",
        "dimension": 2,
        "variable_name": "VELOCITY",
    })
    settings["predictors"] = [{
        "type": "linear_derivative_based",
        "solver": "fluid",
        "data_name": "disp",
        "derivative_data_name": "velocity",
    }]


def set_mapper_options(data, kernel_radius, regularization, polynomial_level, rotation_recovery_mode):
    if all(value is None for value in (
        kernel_radius, regularization, polynomial_level, rotation_recovery_mode
    )):
        return
    operators = data["solver_settings"]["data_transfer_operators"]
    mapper = operators.get("beam_mapper")
    if mapper is None:
        return
    mapper_settings = mapper["mapper_settings"]
    if kernel_radius is not None:
        mapper_settings["kernel_radius"] = kernel_radius
    if regularization is not None:
        mapper_settings["regularization"] = regularization
    if polynomial_level is not None:
        mapper_settings.pop("polynomial_basis", None)
        mapper_settings["polynomial_level"] = polynomial_level
    if rotation_recovery_mode is not None:
        mapper_settings["rotation_recovery_mode"] = rotation_recovery_mode


def set_structural_controls(data, strategy, max_iterations, echo_level):
    settings = data["solver_settings"]
    if strategy is not None:
        settings.pop("line_search", None)
        settings.setdefault("solving_strategy_settings", {})["type"] = strategy
    if max_iterations is not None:
        settings["max_iteration"] = max_iterations
    if echo_level is not None:
        settings["echo_level"] = echo_level
        data["problem_data"]["echo_level"] = echo_level


def set_load_scale(data, scale):
    operations = data["solver_settings"].get("coupling_operations", {})
    operation = operations.get("scale_fluid_reaction")
    if operation is not None:
        operation["scaling_factor"] = scale


def set_beam_section_scale(materials, mdpa_text, scale):
    scaled_values = {
        name: value * scale / 50.0
        for name, value in BEAM_SECTION_AT_SCALE_50.items()
    }

    for prop in materials.get("properties", []):
        variables = prop.get("Material", {}).get("Variables", {})
        for name, value in scaled_values.items():
            if name in variables:
                variables[name] = value

    for name, value in scaled_values.items():
        mdpa_text = re.sub(
            rf"(^\s*{name}\s+).*$",
            rf"\g<1>{value:.17g}",
            mdpa_text,
            flags=re.MULTILINE,
        )
    mdpa_text = re.sub(
        r"linear thickness scale s = \d+(?:\.\d+)?",
        f"linear thickness scale s = {scale:g}",
        mdpa_text,
    )
    return mdpa_text


def set_solid_thickness_scale(materials, mdpa_text, scale):
    for prop in materials.get("properties", []):
        variables = prop.get("Material", {}).get("Variables", {})
        if "THICKNESS" in variables:
            variables["THICKNESS"] = scale
    return re.sub(
        r"(^\s*THICKNESS\s+).*$",
        rf"\g<1>{scale:.17g}",
        mdpa_text,
        flags=re.MULTILINE,
    )


def set_output_paths(csm, cfd, output_dir, point_file_name):
    rel_output = "../../" + str(output_dir.relative_to(ROOT))
    processes = csm["processes"]["list_other_processes"]
    for process in processes:
        params = process.get("Parameters", {})
        output_settings = params.get("output_file_settings")
        if output_settings is not None:
            output_settings["output_path"] = rel_output
            output_settings["file_name"] = point_file_name
    for process in csm["output_processes"].get("vtk_output", []):
        params = process.get("Parameters", {})
        if "output_path" in params:
            params["output_path"] = f"{rel_output}/vtk_output_mok_fsi_csd"
    for process in cfd["output_processes"].get("vtk_output", []):
        params = process.get("Parameters", {})
        if "output_path" in params:
            params["output_path"] = f"{rel_output}/vtk_output_mok_fsi_cfd"
        output_variables = params.get("nodal_solution_step_data_variables")
        if output_variables is not None and "REACTION" not in output_variables:
            output_variables.append("REACTION")


def read_point_history(path):
    values = {}
    if not path.exists():
        return values
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            time = float(parts[0])
            value = float(parts[1])
        except ValueError:
            continue
        if math.isfinite(time) and math.isfinite(value):
            values[round(time, 10)] = value
    return values


def monitor_run(case_dir, log_path, point_path, max_repeated_nonconvergence):
    env = os.environ.copy()
    env["PYTHONPATH"] = str((ROOT / ".." / "bin" / "Release").resolve()) + os.pathsep + env.get("PYTHONPATH", "")
    nonconv_re = re.compile(r"Solver did not converge for step\s+(\d+)")
    step_re = re.compile(r"Fluid Dynamics Analysis: STEP:\s+(\d+)")
    time_re = re.compile(r"Fluid Dynamics Analysis: TIME:\s+([-+0-9.eE]+)")
    last_nonconv_step = None
    repeated_nonconv = 0
    status = "completed"
    failure_step = None
    failure_time = None
    current_step = None
    current_time = None
    interrupted_signal = None
    process = subprocess.Popen(
        [sys.executable, "MainKratos.py"],
        cwd=case_dir,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    def request_graceful_stop(signum, _frame):
        nonlocal interrupted_signal
        interrupted_signal = signum
        if process.poll() is None:
            process.terminate()

    previous_handlers = {
        signum: signal.signal(signum, request_graceful_stop)
        for signum in (signal.SIGINT, signal.SIGTERM)
    }
    try:
        with log_path.open("w") as log:
            assert process.stdout is not None
            for line in process.stdout:
                log.write(line)
                if interrupted_signal is not None:
                    status = f"interrupted_signal_{signal.Signals(interrupted_signal).name}"
                    failure_step = current_step
                    failure_time = current_time
                    break
                step_match = step_re.search(line)
                if step_match:
                    current_step = int(step_match.group(1))
                time_match = time_re.search(line)
                if time_match:
                    current_time = float(time_match.group(1))
                match = nonconv_re.search(line)
                if match:
                    step = int(match.group(1))
                    repeated_nonconv = repeated_nonconv + 1 if step == last_nonconv_step else 1
                    last_nonconv_step = step
                    if repeated_nonconv >= max_repeated_nonconvergence:
                        status = f"stopped_repeated_nonconvergence_step_{step}"
                        failure_step = step
                        failure_time = current_time
                        process.terminate()
                        break
                if "CONVERGENCE WAS NOT ACHIEVED" in line:
                    failure_step = current_step
                    failure_time = current_time
                    status = f"stopped_coupling_nonconvergence_step_{current_step}"
                    process.terminate()
                    break
                if point_path.exists():
                    tail = point_path.read_text().splitlines()[-5:]
                    for point_line in tail:
                        if "nan" in point_line.lower() or "inf" in point_line.lower():
                            status = "stopped_nonfinite_point_output"
                            failure_step = current_step
                            failure_time = current_time
                            process.terminate()
                            break
                if status != "completed":
                    break
            if interrupted_signal is not None and status == "completed":
                status = f"interrupted_signal_{signal.Signals(interrupted_signal).name}"
                failure_step = current_step
                failure_time = current_time
            try:
                returncode = process.wait(timeout=20)
            except subprocess.TimeoutExpired:
                process.kill()
                returncode = process.wait()
    finally:
        for signum, previous_handler in previous_handlers.items():
            signal.signal(signum, previous_handler)
    if returncode != 0 and status == "completed":
        status = f"failed_returncode_{returncode}"
    return status, returncode, failure_step, failure_time


def sha256_text(text):
    return hashlib.sha256(text.encode()).hexdigest()


def write_trial_snapshot(output_dir, args, paths):
    snapshot_dir = output_dir / "input_snapshot"
    snapshot_dir.mkdir(parents=True, exist_ok=False)
    for path in paths:
        if path.exists():
            shutil.copy2(path, snapshot_dir / path.name)
    (snapshot_dir / "trial_arguments.json").write_text(
        json.dumps(vars(args), indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Run one FSI_mok stability trial in a separated output directory.")
    parser.add_argument("case", choices=sorted(CASE_DIRS))
    parser.add_argument("--tag", required=True)
    parser.add_argument("--end-time", type=float, default=10.0)
    parser.add_argument("--dt", type=float)
    parser.add_argument("--alpha", type=float, default=0.03)
    parser.add_argument("--iterations", type=int, default=100)
    parser.add_argument("--kernel-radius", type=float)
    parser.add_argument("--regularization", type=float)
    parser.add_argument("--polynomial-level", type=int, choices=range(5))
    parser.add_argument("--rotation-recovery-mode", choices=("small", "finite"))
    parser.add_argument("--scale", type=float, default=50.0)
    parser.add_argument("--max-repeated-nonconvergence", type=int, default=5)
    parser.add_argument("--structural-strategy", choices=("newton_raphson", "line_search"))
    parser.add_argument("--structural-max-iterations", type=int)
    parser.add_argument("--structural-echo-level", type=int)
    args = parser.parse_args()

    case_dir = ROOT / CASE_DIRS[args.case]
    scale_label = f"s{args.scale:g}".replace(".", "p")
    output_dir = ROOT / "TestCase_Output" / f"StabilityTrials_{scale_label}" / args.tag / args.case
    log_dir = output_dir / "benchmark_logs"
    if output_dir.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing trial output: {output_dir}. Use a new --tag.")
    log_dir.mkdir(parents=True, exist_ok=True)

    cosim_path = case_dir / "ProjectParametersCoSim.json"
    csm_path = case_dir / "ProjectParametersCSM.json"
    cfd_path = case_dir / "ProjectParametersCFD.json"
    materials_path = case_dir / "StructuralMaterials.json"
    mdpa_path = case_dir / "Mok_CSM.mdpa"
    originals = {
        path: path.read_text()
        for path in (cosim_path, csm_path, cfd_path, materials_path, mdpa_path)
        if path.exists()
    }
    original_hashes = {str(path): sha256_text(text) for path, text in originals.items()}
    status = "failed_before_run"
    returncode = None
    failure_step = None
    failure_time = None
    caught_exception = None

    point_file_name = POINT_FILE_NAMES[args.case]
    point_path = output_dir / f"{point_file_name}.dat"
    try:
        cosim = read_json(cosim_path)
        csm = read_json(csm_path)
        cfd = read_json(cfd_path)
        materials = read_json(materials_path) if materials_path.exists() else None
        for data in (cosim, csm, cfd):
            set_end_time(data, args.end_time)
        if args.dt is not None:
            set_time_step(csm, args.dt)
            set_time_step(cfd, args.dt)
        set_coupling(cosim, args.alpha, args.iterations)
        set_mapper_options(
            cosim,
            args.kernel_radius,
            args.regularization,
            args.polynomial_level,
            args.rotation_recovery_mode,
        )
        set_structural_controls(
            csm,
            args.structural_strategy,
            args.structural_max_iterations,
            args.structural_echo_level,
        )
        set_load_scale(cosim, args.scale)
        set_output_paths(csm, cfd, output_dir, point_file_name)
        if materials is not None and mdpa_path.exists():
            if args.case == "NearestNeighbor":
                mdpa_text = set_solid_thickness_scale(materials, mdpa_path.read_text(), args.scale)
            else:
                mdpa_text = set_beam_section_scale(materials, mdpa_path.read_text(), args.scale)
            write_json(materials_path, materials)
            mdpa_path.write_text(mdpa_text)
        write_json(cosim_path, cosim)
        write_json(csm_path, csm)
        write_json(cfd_path, cfd)

        write_trial_snapshot(
            output_dir,
            args,
            (cosim_path, csm_path, cfd_path, materials_path, mdpa_path),
        )

        status, returncode, failure_step, failure_time = monitor_run(
            case_dir,
            log_dir / f"{args.case}.log",
            point_path,
            args.max_repeated_nonconvergence,
        )
    except Exception as exc:
        caught_exception = f"{type(exc).__name__}: {exc}"
        status = f"failed_exception_{type(exc).__name__}"
    finally:
        for path, text in originals.items():
            path.write_text(text)

    restored_hashes = {
        str(path): sha256_text(path.read_text())
        for path in originals
        if path.exists()
    }
    restoration_verified = restored_hashes == original_hashes

    effective_dt = args.dt
    if effective_dt is None:
        original_csm = json.loads(originals[csm_path])
        effective_dt = original_csm["solver_settings"]["time_stepping"]["time_step"]
    valid_end_time = args.end_time
    if failure_time is not None:
        valid_end_time = min(valid_end_time, max(0.0, failure_time - effective_dt))

    summary = {
        "case": args.case,
        "tag": args.tag,
        "status": status,
        "returncode": returncode,
        "exception": caught_exception,
        "end_time_target": args.end_time,
        "dt": args.dt,
        "alpha": args.alpha,
        "iterations": args.iterations,
        "kernel_radius": args.kernel_radius,
        "regularization": args.regularization,
        "polynomial_level": args.polynomial_level,
        "rotation_recovery_mode": args.rotation_recovery_mode,
        "accuracy_reference": "external analytical verification only",
        "scale": args.scale,
        "structural_strategy": args.structural_strategy,
        "structural_max_iterations": args.structural_max_iterations,
        "structural_echo_level": args.structural_echo_level,
        "effective_dt": effective_dt,
        "failure_step": failure_step,
        "failure_time": failure_time,
        "valid_end_time_limit": valid_end_time,
        "input_snapshot_dir": str(output_dir / "input_snapshot"),
        "restoration_verified": restoration_verified,
        "point_file": str(point_path),
        "output_dir": str(output_dir),
    }
    history = read_point_history(point_path)
    if history:
        summary["last_finite_time"] = max(history)
        summary["last_finite_value"] = history[max(history)]
        valid_times = [time for time in history if time <= valid_end_time + 1.0e-10]
        if valid_times:
            summary["last_valid_time"] = max(valid_times)
            summary["last_valid_value"] = history[max(valid_times)]
    summary["rel_l2_skipped_reason"] = (
        "FSI is a stability/failure-mechanism trial; mapper accuracy is verified "
        "against prescribed analytical fields by the dedicated scripts."
    )
    summary_path = output_dir.parent / f"{args.case}_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    conclusion_lines = [
        f"status: {status}",
        f"last_valid_time: {summary.get('last_valid_time')}",
        "accuracy_reference: external analytical verification only",
        f"restoration_verified: {restoration_verified}",
    ]
    (output_dir.parent / "conclusion.txt").write_text("\n".join(conclusion_lines) + "\n")
    print(json.dumps(summary, indent=2))
    return 0 if status == "completed" and restoration_verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
