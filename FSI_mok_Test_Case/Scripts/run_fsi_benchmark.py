#!/usr/bin/env python3
import argparse
import csv
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

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
from analyze_vtk_history import validate_vtk_directory


ROOT = SCRIPT_DIR.parent

# Canonical FSI accuracy reference.  This is a physically consistent
# NearestNeighbor run with s=50, dt=0.05, alpha=0.03 and 60 coupling
# iterations, sampled at the same point A as the beam cases.
NN_REFERENCE_TAG = "nearest_neighbor_s50_dt005_a003_i60_t045_20260724"
NN_REFERENCE_POINT_FILE = (
    ROOT
    / "TestCase_Output"
    / "StabilityTrials_s50"
    / NN_REFERENCE_TAG
    / "NearestNeighbor"
    / "point_A_thickness50_scaled_load.dat"
)
NN_COMPARISON_SUMMARY_JSON = (
    ROOT / "TestCase_Output" / "StabilityTrials_s50" / "nn_rel_l2_comparison.json"
)
NN_COMPARISON_SUMMARY_CSV = NN_COMPARISON_SUMMARY_JSON.with_suffix(".csv")

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

BEAM_CSM_REQUIRED_VTK_FIELDS = (
    "DISPLACEMENT",
    "ROTATION",
    "POINT_LOAD",
    "POINT_MOMENT",
    "REACTION",
    "REACTION_MOMENT",
)
NN_CSM_REQUIRED_VTK_FIELDS = ("DISPLACEMENT", "POINT_LOAD")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path, data):
    path.write_text(
        json.dumps(data, indent=4, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def set_end_time(data, end_time):
    data["problem_data"]["end_time"] = end_time


def set_time_step(data, dt):
    settings = data["solver_settings"]
    if "time_stepping" in settings:
        settings["time_stepping"]["time_step"] = dt
    if "fluid_solver_settings" in settings:
        settings["fluid_solver_settings"]["time_stepping"]["time_step"] = dt


def set_coupling(data, alpha, iterations, predictor, accelerator):
    settings = data["solver_settings"]
    settings["num_coupling_iterations"] = iterations
    accelerators = settings.setdefault("convergence_accelerators", [])
    if not accelerators:
        accelerators.append({})
    accelerators[0].clear()
    if accelerator == "aitken":
        accelerators[0].update({
            "type": "aitken",
            "solver": "fluid",
            "data_name": "disp",
            "init_alpha": alpha,
            "init_alpha_max": alpha,
            "alpha_min": -0.5,
            "alpha_max": 0.5,
        })
    else:
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
    if predictor == "linear_derivative_based":
        settings["predictors"] = [{
            "type": "linear_derivative_based",
            "solver": "fluid",
            "data_name": "disp",
            "derivative_data_name": "velocity",
        }]
    else:
        settings["predictors"] = []


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
    mapper_type = mapper_settings.get("mapper_type")
    if mapper_type != "beam_spline_mapper_with_recovery_of_rotations":
        # The ordinary beam_spline_mapper has a fixed cubic Hermite-RBF
        # formulation and does not expose the recovery mapper's kernel,
        # regularization, polynomial, or rotation-mode controls.
        return
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
    point_output_index = 0
    for process in processes:
        params = process.get("Parameters", {})
        output_settings = params.get("output_file_settings")
        if output_settings is not None:
            output_settings["output_path"] = rel_output
            # The NN solid case has separate point-A and point-B processes.
            # Only point A is the canonical accuracy signal; assigning the
            # same filename to both processes corrupts the stream.
            if point_output_index == 0:
                output_settings["file_name"] = point_file_name
            point_output_index += 1
    for process in csm["output_processes"].get("vtk_output", []):
        params = process.get("Parameters", {})
        if "output_path" in params:
            params["output_path"] = f"{rel_output}/vtk_output_mok_fsi_csd"
        params["output_control_type"] = "step"
        params["output_interval"] = 1
        params["file_format"] = "ascii"
        params["output_precision"] = 15
    for process in cfd["output_processes"].get("vtk_output", []):
        params = process.get("Parameters", {})
        if "output_path" in params:
            params["output_path"] = f"{rel_output}/vtk_output_mok_fsi_cfd"
        params["output_control_type"] = "step"
        params["output_interval"] = 1
        params["file_format"] = "ascii"
        params["output_precision"] = 15
        output_variables = params.get("nodal_solution_step_data_variables")
        if output_variables is not None and "REACTION" not in output_variables:
            output_variables.append("REACTION")


def validate_output_configuration(case, csm, cfd, output_dir):
    """Fail before launching FSI if VTK/point streams are ambiguous."""
    errors = []
    required_csm = (
        NN_CSM_REQUIRED_VTK_FIELDS
        if case == "NearestNeighbor"
        else BEAM_CSM_REQUIRED_VTK_FIELDS
    )
    vtk_specs = {}
    for label, data, folder, required in (
        ("csd", csm, "vtk_output_mok_fsi_csd", required_csm),
        ("cfd", cfd, "vtk_output_mok_fsi_cfd", ()),
    ):
        processes = data.get("output_processes", {}).get("vtk_output", [])
        if len(processes) != 1:
            errors.append(f"{label}: expected one VTK output process")
            continue
        params = processes[0].get("Parameters", {})
        variables = tuple(params.get("nodal_solution_step_data_variables", ()))
        missing = sorted(set(required) - set(variables))
        if missing:
            errors.append(f"{label}: missing required fields {missing}")
        if params.get("output_control_type") != "step":
            errors.append(f"{label}: output_control_type must be step")
        if int(params.get("output_interval", 0)) != 1:
            errors.append(f"{label}: output_interval must be 1")
        if params.get("file_format") != "ascii":
            errors.append(f"{label}: file_format must be ascii")
        expected_suffix = f"/{folder}"
        if not str(params.get("output_path", "")).endswith(expected_suffix):
            errors.append(f"{label}: unexpected output_path")
        vtk_specs[label] = {
            "directory": str((output_dir / folder).resolve()),
            "required_fields": list(required if label == "csd" else variables),
            "declared_fields": list(variables),
            "file_format": params.get("file_format"),
            "output_precision": params.get("output_precision"),
            "output_interval": params.get("output_interval"),
        }

    point_processes = csm.get("processes", {}).get("list_other_processes", [])
    point_specs = []
    point_names = set()
    for process in point_processes:
        params = process.get("Parameters", {})
        settings = params.get("output_file_settings")
        if settings is None:
            continue
        name = settings.get("file_name")
        if name in point_names:
            errors.append(f"Duplicate point-output filename: {name}")
        point_names.add(name)
        point_specs.append(
            {
                "file_name": name,
                "position": params.get("position"),
                "output_variables": params.get("output_variables", []),
            }
        )
    if not point_specs:
        errors.append("No point-output process configured")
    elif not point_specs[0]["output_variables"] or point_specs[0][
        "output_variables"
    ][0] != "DISPLACEMENT_X":
        errors.append("Canonical point stream must start with DISPLACEMENT_X")

    if errors:
        raise RuntimeError("VTK preflight failed: " + "; ".join(errors))
    return {
        "valid": True,
        "vtk": vtk_specs,
        "point_outputs": point_specs,
    }


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


def compute_relative_l2(reference_path, result_path, max_time):
    """Compute rel L2 on common, finite sample times up to max_time.

    e_relL2 = sqrt(sum_t (u_result(t)-u_reference(t))^2
                     / sum_t u_reference(t)^2)

    No interpolation is used.  NaN/Inf values are discarded by
    read_point_history, and max_time must be the last valid/converged time of
    the result being assessed.
    """
    reference = read_point_history(reference_path)
    result = read_point_history(result_path)
    common_times = sorted(set(reference) & set(result))
    common_times = [
        time for time in common_times
        if time <= max_time + 1.0e-10
    ]
    if not common_times:
        raise RuntimeError(
            f"No common finite samples between {reference_path} and "
            f"{result_path} up to t={max_time}."
        )
    denominator = sum(reference[time] ** 2 for time in common_times)
    if denominator <= 0.0:
        raise RuntimeError(
            f"NearestNeighbor reference norm is zero on "
            f"{common_times[0]}...{common_times[-1]}."
        )
    numerator = sum(
        (result[time] - reference[time]) ** 2
        for time in common_times
    )
    return {
        "accuracy_reference": "NearestNeighbor point-A DISPLACEMENT_X",
        "reference_point_file": str(Path(reference_path).resolve()),
        "common_steps": len(common_times),
        "start_time": common_times[0],
        "end_time": common_times[-1],
        "rel_l2": math.sqrt(numerator / denominator),
    }


def write_conclusion(summary, conclusion_path):
    lines = [
        f"status: {summary.get('status')}",
        f"last_valid_time: {summary.get('last_valid_time')}",
        f"accuracy_reference: {summary.get('accuracy_reference')}",
        f"reference_point_file: {summary.get('reference_point_file')}",
        f"rel_l2_interval: {summary.get('start_time')}...{summary.get('end_time')}",
        f"rel_l2_common_steps: {summary.get('common_steps')}",
        f"rel_l2: {summary.get('rel_l2')}",
        f"vtk_valid: {summary.get('vtk_valid')}",
        f"vtk_validation_manifest: {summary.get('vtk_validation_manifest')}",
        f"restoration_verified: {summary.get('restoration_verified')}",
    ]
    conclusion_path.write_text("\n".join(lines) + "\n")


def update_existing_comparison(path):
    """Recompute NN rel L2 and update one existing trial summary in place."""
    candidate = path.resolve()
    output_root = (ROOT / "TestCase_Output").resolve()
    try:
        candidate.relative_to(output_root)
    except ValueError as exc:
        raise RuntimeError(
            f"Existing comparison target must be below {output_root}: "
            f"{candidate}"
        ) from exc

    if candidate.is_file():
        summary_path = candidate
    else:
        summary_paths = sorted(candidate.glob("*_summary.json"))
        if len(summary_paths) != 1:
            raise RuntimeError(
                f"Expected exactly one *_summary.json in {candidate}, "
                f"found {len(summary_paths)}."
            )
        summary_path = summary_paths[0]

    summary = read_json(summary_path)
    if float(summary.get("scale", 50.0)) != 50.0:
        raise RuntimeError(
            f"The canonical NN reference has scale 50, but {summary_path} "
            f"records scale={summary.get('scale')}."
        )
    result_path = Path(summary["point_file"])
    max_time = summary.get("last_valid_time")
    if max_time is None:
        raise RuntimeError(f"No last_valid_time in {summary_path}.")
    summary.update(
        compute_relative_l2(
            NN_REFERENCE_POINT_FILE,
            result_path,
            float(max_time),
        )
    )
    summary.pop("rel_l2_skipped_reason", None)
    write_json(summary_path, summary)
    write_conclusion(summary, summary_path.parent / "conclusion.txt")
    return {
        "summary_path": str(summary_path),
        "case": summary["case"],
        "tag": summary.get("tag"),
        "rotation_recovery_mode": summary.get("rotation_recovery_mode"),
        "common_steps": summary["common_steps"],
        "start_time": summary["start_time"],
        "end_time": summary["end_time"],
        "rel_l2": summary["rel_l2"],
    }


def monitor_run(case_dir, log_path, point_path, max_repeated_nonconvergence):
    env = os.environ.copy()
    env["PYTHONPATH"] = str((ROOT / ".." / "bin" / "Release").resolve()) + os.pathsep + env.get("PYTHONPATH", "")
    current_mapping_library = str(
        (ROOT / ".." / "build" / "applications" / "MappingApplication").resolve()
    )
    release_libraries = str((ROOT / ".." / "bin" / "Release" / "libs").resolve())
    env["LD_LIBRARY_PATH"] = (
        current_mapping_library
        + os.pathsep
        + release_libraries
        + os.pathsep
        + env.get("LD_LIBRARY_PATH", "")
    )
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
        # Keep terminal SIGINT/SIGTERM from killing runner and solver as one
        # process group. The runner handles the signal, terminates the solver,
        # and reaches its input-restoration and manifest-writing finally path.
        start_new_session=True,
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
        json.dumps(vars(args), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def main():
    parser = argparse.ArgumentParser(description="Run one FSI_mok stability trial in a separated output directory.")
    parser.add_argument("case", nargs="?", choices=sorted(CASE_DIRS))
    parser.add_argument("--tag")
    parser.add_argument(
        "--compare-existing",
        nargs="+",
        type=Path,
        help=(
            "Do not run FSI. Recompute NN-reference rel L2 for one or more "
            "existing trial directories and update their summary/conclusion."
        ),
    )
    parser.add_argument("--end-time", type=float, default=10.0)
    parser.add_argument("--dt", type=float)
    parser.add_argument("--alpha", type=float, default=0.03)
    parser.add_argument("--iterations", type=int, default=100)
    parser.add_argument(
        "--predictor",
        choices=("linear_derivative_based", "none"),
        default="linear_derivative_based",
    )
    parser.add_argument(
        "--accelerator",
        choices=("constant_relaxation", "aitken"),
        default="constant_relaxation",
    )
    parser.add_argument("--kernel-radius", type=float)
    parser.add_argument("--regularization", type=float)
    parser.add_argument("--polynomial-level", type=int, choices=range(8))
    parser.add_argument("--rotation-recovery-mode", choices=("small", "finite"))
    parser.add_argument("--scale", type=float, default=50.0)
    parser.add_argument("--max-repeated-nonconvergence", type=int, default=5)
    parser.add_argument("--structural-strategy", choices=("newton_raphson", "line_search"))
    parser.add_argument("--structural-max-iterations", type=int)
    parser.add_argument("--structural-echo-level", type=int)
    args = parser.parse_args()

    if args.compare_existing:
        if not NN_REFERENCE_POINT_FILE.exists():
            raise FileNotFoundError(
                f"Canonical NearestNeighbor reference is missing: "
                f"{NN_REFERENCE_POINT_FILE}"
            )
        comparisons = [
            update_existing_comparison(path)
            for path in args.compare_existing
        ]
        shared_end_time = min(item["end_time"] for item in comparisons)
        for item in comparisons:
            summary = read_json(Path(item["summary_path"]))
            shared = compute_relative_l2(
                NN_REFERENCE_POINT_FILE,
                Path(summary["point_file"]),
                shared_end_time,
            )
            item["shared_common_steps"] = shared["common_steps"]
            item["shared_start_time"] = shared["start_time"]
            item["shared_end_time"] = shared["end_time"]
            item["shared_rel_l2"] = shared["rel_l2"]
        comparison_record = {
            "reference": "NearestNeighbor point-A DISPLACEMENT_X",
            "reference_point_file": str(NN_REFERENCE_POINT_FILE),
            "definition": (
                "sqrt(sum_common_valid((u_mapper-u_NN)^2)"
                "/sum_common_valid(u_NN^2)); no interpolation"
            ),
            "comparisons": comparisons,
        }
        write_json(NN_COMPARISON_SUMMARY_JSON, comparison_record)
        with NN_COMPARISON_SUMMARY_CSV.open("w", newline="") as output:
            writer = csv.DictWriter(
                output,
                fieldnames=(
                    "case",
                    "tag",
                    "rotation_recovery_mode",
                    "common_steps",
                    "start_time",
                    "end_time",
                    "rel_l2",
                    "shared_common_steps",
                    "shared_start_time",
                    "shared_end_time",
                    "shared_rel_l2",
                    "summary_path",
                ),
            )
            writer.writeheader()
            writer.writerows(comparisons)
        print(json.dumps(comparison_record, indent=2, ensure_ascii=False))
        return 0
    if args.case is None or args.tag is None:
        parser.error("case and --tag are required unless --compare-existing is used")

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
    vtk_preflight = None

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
        set_coupling(
            cosim, args.alpha, args.iterations, args.predictor, args.accelerator
        )
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
        vtk_preflight = validate_output_configuration(
            args.case, csm, cfd, output_dir
        )
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
        "predictor": args.predictor,
        "accelerator": args.accelerator,
        "kernel_radius": args.kernel_radius,
        "regularization": args.regularization,
        "polynomial_level": args.polynomial_level,
        "rotation_recovery_mode": args.rotation_recovery_mode,
        "accuracy_reference": "NearestNeighbor point-A DISPLACEMENT_X",
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
    if args.scale == 50.0 and NN_REFERENCE_POINT_FILE.exists():
        try:
            summary.update(
                compute_relative_l2(
                    NN_REFERENCE_POINT_FILE,
                    point_path,
                    float(summary.get("last_valid_time", valid_end_time)),
                )
            )
        except RuntimeError as exc:
            summary["rel_l2_skipped_reason"] = str(exc)
    else:
        summary["rel_l2_skipped_reason"] = (
            f"Canonical same-scale NN reference unavailable: "
            f"{NN_REFERENCE_POINT_FILE}"
            if args.scale == 50.0
            else "Canonical NN reference has scale 50; scale mismatch."
        )

    vtk_validation = {
        "valid": False,
        "preflight": vtk_preflight,
        "errors": [],
    }
    if vtk_preflight is None:
        vtk_validation["errors"].append("VTK preflight did not complete")
    else:
        try:
            last_valid_time = float(
                summary.get("last_valid_time", valid_end_time)
            )
            point_position = vtk_preflight["point_outputs"][0]["position"]
            csd_spec = vtk_preflight["vtk"]["csd"]
            cfd_spec = vtk_preflight["vtk"]["cfd"]
            csd_manifest = validate_vtk_directory(
                csd_spec["directory"],
                effective_dt,
                required_fields=csd_spec["required_fields"],
                max_valid_time=last_valid_time,
                point_history=point_path,
                point_position=point_position,
                point_component=0,
                # PointOutputProcess selects by the requested coordinate while
                # the legacy VTK mesh stores the corresponding node coordinate
                # with an O(1e-8 m) geometry discrepancy in the Mok mesh.
                point_rel_l2_tolerance=1.0e-7,
            )
            cfd_manifest = validate_vtk_directory(
                cfd_spec["directory"],
                effective_dt,
                required_fields=cfd_spec["required_fields"],
                max_valid_time=last_valid_time,
            )
            vtk_validation.update(
                {
                    "csd": csd_manifest,
                    "cfd": cfd_manifest,
                    "valid": csd_manifest["valid"] and cfd_manifest["valid"],
                }
            )
            vtk_validation["errors"].extend(csd_manifest["errors"])
            vtk_validation["errors"].extend(cfd_manifest["errors"])
        except Exception as exc:
            vtk_validation["errors"].append(
                f"{type(exc).__name__}: {exc}"
            )
    vtk_manifest_path = output_dir / "vtk_validation_manifest.json"
    write_json(vtk_manifest_path, vtk_validation)
    summary["vtk_validation_manifest"] = str(vtk_manifest_path)
    summary["vtk_valid"] = vtk_validation["valid"]
    if not vtk_validation["valid"] and status == "completed":
        status = "failed_vtk_validation"
        summary["status"] = status

    summary_path = output_dir.parent / f"{args.case}_summary.json"
    write_json(summary_path, summary)
    write_conclusion(summary, output_dir.parent / "conclusion.txt")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0 if (
        status == "completed"
        and restoration_verified
        and vtk_validation["valid"]
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
