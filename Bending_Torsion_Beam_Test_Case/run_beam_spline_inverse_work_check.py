from pathlib import Path
from contextlib import contextmanager
import argparse
import ctypes
import datetime
import math
import os
import sys
import time

import KratosMultiphysics as KM
import KratosMultiphysics.LinearSolversApplication
import KratosMultiphysics.MappingApplication
import KratosMultiphysics.StructuralMechanicsApplication


ROOT = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT / "Bending_Torsion_Beam_Test_Case" / "Bending_Torsion_Beam_Test_Case"
BEAM_MDPA_PATH = CASE_DIR / "beam_geometry" / "beam_geometry_coarse.mdpa"
SURFACE_MDPA_PATH = CASE_DIR / "fluid_surface_mesh" / "fluid_surface_mesh_1x5.mdpa"
OUTPUT_ROOT = ROOT / "Bending_Torsion_Beam_Test_Case" / "TestCase_Output" / "InverseMap_Work_Check"
BEAM_LENGTH = 10.0
FINITE_DIFFERENCE_EPSILON = 1.0e-7
RELATIVE_TOLERANCE = 2.0e-5
ABSOLUTE_TOLERANCE = 2.0e-8
LIBC = ctypes.CDLL(None)


TEST_CASES = [
    ("bending_y", 0.25 * math.pi, "constant_y", "bending_y"),
    ("bending_z", 0.20 * math.pi, "constant_z", "bending_z"),
    ("torsion", 0.30 * math.pi, "torsional_ring", "torsion"),
    ("displacement_y", 1.0, "mixed", "mixed"),
]

MESH_PRESETS = {
    "small": (
        CASE_DIR / "beam_geometry" / "beam_geometry_coarse.mdpa",
        CASE_DIR / "fluid_surface_mesh" / "fluid_surface_mesh_1x5.mdpa",
    ),
    "dense_surface": (
        CASE_DIR / "beam_geometry" / "beam_geometry_coarse.mdpa",
        CASE_DIR / "fluid_surface_mesh" / "fluid_surface_mesh_10x100.mdpa",
    ),
    "fine": (
        CASE_DIR / "beam_geometry" / "beam_geometry.mdpa",
        CASE_DIR / "fluid_surface_mesh" / "fluid_surface_mesh_10x100.mdpa",
    ),
}


def prescribed_beam_kinematics(mode, x, parameter_value):
    xi = x / BEAM_LENGTH

    if mode == "bending_y":
        u = 0.0
        v = 0.5 * parameter_value * BEAM_LENGTH * xi * xi
        w = 0.0
        theta_x = 0.0
        theta_y = 0.0
        theta_z = parameter_value * xi
        return (u, v, w), (theta_x, theta_y, theta_z)

    if mode == "bending_z":
        u = 0.0
        v = 0.0
        w = 0.5 * parameter_value * BEAM_LENGTH * xi * xi
        theta_x = 0.0
        theta_y = -parameter_value * xi
        theta_z = 0.0
        return (u, v, w), (theta_x, theta_y, theta_z)

    if mode == "torsion":
        return (0.0, 0.0, 0.0), (parameter_value * xi, 0.0, 0.0)

    if mode == "displacement_y":
        u = 0.0
        v = 0.5 * parameter_value * (3.0 * xi * xi - xi * xi * xi)
        w = 0.0
        theta_x = 0.0
        theta_y = 0.0
        theta_z = 1.5 * parameter_value * (2.0 * xi - xi * xi) / BEAM_LENGTH
        return (u, v, w), (theta_x, theta_y, theta_z)

    raise RuntimeError(f"Unsupported kinematics mode: {mode}")


def perturbation(mode, x):
    xi = x / BEAM_LENGTH

    if mode == "bending_y":
        return (0.0, 0.15 * xi * xi, 0.0), (0.0, 0.0, 0.03 * xi)

    if mode == "bending_z":
        return (0.0, 0.0, -0.12 * xi * (1.0 + xi)), (0.0, 0.025 * xi, 0.0)

    if mode == "torsion":
        return (0.0, 0.0, 0.0), (0.04 * xi, 0.0, 0.0)

    if mode == "mixed":
        displacement = (
            0.03 * xi,
            0.09 * xi * xi,
            -0.07 * xi * (1.0 - 0.25 * xi),
        )
        rotation = (
            0.02 * xi,
            0.015 * (1.0 - xi) * xi,
            0.025 * xi * xi,
        )
        return displacement, rotation

    raise RuntimeError(f"Unsupported perturbation mode: {mode}")


def surface_force(force_mode, node):
    x = node.X0
    y = node.Y0
    z = node.Z0

    if force_mode == "constant_y":
        return (0.0, 1.0, 0.0)

    if force_mode == "constant_z":
        return (0.0, 0.0, -0.8)

    if force_mode == "torsional_ring":
        return (0.0, -z, y)

    if force_mode == "mixed":
        return (
            0.05 + 0.01 * x,
            0.8 + 0.05 * y + 0.02 * x,
            -0.35 + 0.04 * z - 0.01 * x,
        )

    raise RuntimeError(f"Unsupported force mode: {force_mode}")


def create_structure_model_part(name, kinematics_mode, parameter_value):
    model = KM.Model()
    structure = model.CreateModelPart(name)
    structure.ProcessInfo[KM.DOMAIN_SIZE] = 3
    structure.ProcessInfo[KM.TIME] = 0.0
    structure.ProcessInfo[KM.DELTA_TIME] = 1.0
    structure.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    structure.AddNodalSolutionStepVariable(KM.ROTATION)
    structure.AddNodalSolutionStepVariable(KM.FORCE)
    structure.AddNodalSolutionStepVariable(KM.MOMENT)

    KM.ModelPartIO(str(BEAM_MDPA_PATH.with_suffix(""))).ReadModelPart(structure)
    beam = structure.GetSubModelPart("Parts_Beam_beam")
    set_structure_kinematics(beam, kinematics_mode, parameter_value)
    return beam


def create_surface_model_part(name):
    model = KM.Model()
    surface_root = model.CreateModelPart(name)
    surface_root.ProcessInfo[KM.DOMAIN_SIZE] = 3
    surface_root.ProcessInfo[KM.TIME] = 0.0
    surface_root.ProcessInfo[KM.DELTA_TIME] = 1.0
    surface_root.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    surface_root.AddNodalSolutionStepVariable(KM.REACTION)

    KM.ModelPartIO(str(SURFACE_MDPA_PATH.with_suffix(""))).ReadModelPart(surface_root)
    return surface_root.GetSubModelPart("Parts_Shell_wet_surface")


def set_structure_kinematics(beam, kinematics_mode, parameter_value):
    for node in beam.Nodes:
        displacement, rotation = prescribed_beam_kinematics(kinematics_mode, node.X0, parameter_value)
        node.SetSolutionStepValue(KM.DISPLACEMENT, list(displacement))
        node.SetSolutionStepValue(KM.ROTATION, list(rotation))


def apply_perturbation(beam, perturbation_mode, epsilon):
    for node in beam.Nodes:
        displacement = list(node.GetSolutionStepValue(KM.DISPLACEMENT))
        rotation = list(node.GetSolutionStepValue(KM.ROTATION))
        delta_displacement, delta_rotation = perturbation(perturbation_mode, node.X0)
        for i in range(3):
            displacement[i] += epsilon * delta_displacement[i]
            rotation[i] += epsilon * delta_rotation[i]
        node.SetSolutionStepValue(KM.DISPLACEMENT, displacement)
        node.SetSolutionStepValue(KM.ROTATION, rotation)


def assign_surface_forces(surface, force_mode):
    for node in surface.Nodes:
        node.SetSolutionStepValue(KM.REACTION, list(surface_force(force_mode, node)))


def create_beam_spline_mapper(origin, destination):
    settings = KM.Parameters("""{
        "mapper_type" : "beam_spline_mapper",
        "search_settings" : {
            "search_radius" : 3.0,
            "max_num_search_iterations" : 30
        },
        "local_coord_tolerance" : 0.25,
        "echo_level" : 0
    }""")
    return KM.MapperFactory.CreateMapper(origin, destination, settings)


def vector_to_tuple(value):
    return (float(value[0]), float(value[1]), float(value[2]))


def dot(lhs, rhs):
    return sum(lhs[i] * rhs[i] for i in range(3))


def capture_surface_displacements(surface):
    return {
        node.Id: vector_to_tuple(node.GetSolutionStepValue(KM.DISPLACEMENT))
        for node in surface.Nodes
    }


def surface_directional_work(surface, base_displacements, epsilon):
    work = 0.0
    for node in surface.Nodes:
        force = vector_to_tuple(node.GetSolutionStepValue(KM.REACTION))
        displacement = vector_to_tuple(node.GetSolutionStepValue(KM.DISPLACEMENT))
        base_displacement = base_displacements[node.Id]
        displacement_derivative = tuple((displacement[i] - base_displacement[i]) / epsilon for i in range(3))
        work += dot(force, displacement_derivative)
    return work


def beam_generalized_work(beam, perturbation_mode):
    work = 0.0
    for node in beam.Nodes:
        force = vector_to_tuple(node.GetSolutionStepValue(KM.FORCE))
        moment = vector_to_tuple(node.GetSolutionStepValue(KM.MOMENT))
        delta_displacement, delta_rotation = perturbation(perturbation_mode, node.X0)
        work += dot(force, delta_displacement) + dot(moment, delta_rotation)
    return work


def run_case(case_index, kinematics_mode, parameter_value, force_mode, perturbation_mode, output_path):
    case_start_time = time.perf_counter()

    print("BEAM SPLINE INVERSE MAP ADJOINT CHECK")
    print("case_index:", case_index)
    print("kinematics_mode:", kinematics_mode)
    print("parameter_value:", parameter_value)
    print("force_mode:", force_mode)
    print("perturbation_mode:", perturbation_mode)
    print("epsilon:", FINITE_DIFFERENCE_EPSILON)
    print("beam mdpa:", BEAM_MDPA_PATH)
    print("surface mdpa:", SURFACE_MDPA_PATH)
    print("output_path:", output_path)
    print("")

    setup_start_time = time.perf_counter()
    beam = create_structure_model_part("origin_beam_structure", kinematics_mode, parameter_value)
    surface = create_surface_model_part("beam_spline_surface")
    mapper = create_beam_spline_mapper(beam, surface)
    setup_time = time.perf_counter() - setup_start_time

    forward_map_start_time = time.perf_counter()
    mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
    forward_map_time = time.perf_counter() - forward_map_start_time
    base_surface_displacements = capture_surface_displacements(surface)

    assign_surface_forces(surface, force_mode)
    inverse_map_start_time = time.perf_counter()
    mapper.InverseMap(KM.FORCE, KM.MOMENT, KM.REACTION, KM.Flags())
    inverse_map_time = time.perf_counter() - inverse_map_start_time
    beam_work = beam_generalized_work(beam, perturbation_mode)

    apply_perturbation(beam, perturbation_mode, FINITE_DIFFERENCE_EPSILON)
    perturbed_forward_map_start_time = time.perf_counter()
    mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
    perturbed_forward_map_time = time.perf_counter() - perturbed_forward_map_start_time
    surface_work = surface_directional_work(surface, base_surface_displacements, FINITE_DIFFERENCE_EPSILON)

    absolute_error = abs(surface_work - beam_work)
    reference = max(abs(surface_work), abs(beam_work), 1.0)
    relative_error = absolute_error / reference
    passed = absolute_error <= ABSOLUTE_TOLERANCE or relative_error <= RELATIVE_TOLERANCE
    total_time = time.perf_counter() - case_start_time

    print("mesh:", f"beam={beam.NumberOfNodes()} nodes/{beam.NumberOfElements()} elems,",
          f"surface={surface.NumberOfNodes()} nodes/{surface.NumberOfElements()} elems/{surface.NumberOfConditions()} conditions")
    print("surface_directional_work:", f"{surface_work:.16e}")
    print("beam_generalized_work:", f"{beam_work:.16e}")
    print("absolute_error:", f"{absolute_error:.16e}")
    print("relative_error:", f"{relative_error:.16e}")
    print("timings_seconds:")
    print("  setup:", f"{setup_time:.6f}")
    print("  forward_map:", f"{forward_map_time:.6f}")
    print("  inverse_map:", f"{inverse_map_time:.6f}")
    print("  perturbed_forward_map:", f"{perturbed_forward_map_time:.6f}")
    print("  total:", f"{total_time:.6f}")
    print("tolerances:", f"absolute={ABSOLUTE_TOLERANCE:.3e}", f"relative={RELATIVE_TOLERANCE:.3e}")
    print("status:", "PASS" if passed else "FAIL")

    return {
        "case_index": case_index,
        "kinematics_mode": kinematics_mode,
        "parameter_value": parameter_value,
        "force_mode": force_mode,
        "perturbation_mode": perturbation_mode,
        "surface_work": surface_work,
        "beam_work": beam_work,
        "absolute_error": absolute_error,
        "relative_error": relative_error,
        "setup_time": setup_time,
        "forward_map_time": forward_map_time,
        "inverse_map_time": inverse_map_time,
        "perturbed_forward_map_time": perturbed_forward_map_time,
        "total_time": total_time,
        "beam_nodes": beam.NumberOfNodes(),
        "beam_elements": beam.NumberOfElements(),
        "surface_nodes": surface.NumberOfNodes(),
        "surface_elements": surface.NumberOfElements(),
        "surface_conditions": surface.NumberOfConditions(),
        "passed": passed,
    }


@contextmanager
def redirect_stdout_to_file(output_file):
    output_file.parent.mkdir(parents=True, exist_ok=True)
    saved_stdout = os.dup(1)
    try:
        with open(output_file, "w") as log_file:
            flush_stdout()
            os.dup2(log_file.fileno(), 1)
            yield
    finally:
        flush_stdout()
        os.dup2(saved_stdout, 1)
        os.close(saved_stdout)


def flush_stdout():
    sys.stdout.flush()
    LIBC.fflush(None)


def format_pi_multiple(value):
    coefficient = value / math.pi
    return f"{coefficient:g}pi"


def output_folder_name(kinematics_mode, parameter_value):
    if kinematics_mode == "bending_y":
        return f"theta_y={format_pi_multiple(parameter_value)}"
    if kinematics_mode == "bending_z":
        return f"theta_z={format_pi_multiple(parameter_value)}"
    if kinematics_mode == "torsion":
        return f"torsion={format_pi_multiple(parameter_value)}"
    if kinematics_mode == "displacement_y":
        return f"u_y={parameter_value:g}"
    raise RuntimeError(f"Unsupported kinematics mode: {kinematics_mode}")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Check BeamSplineMapper InverseMap accuracy and timings with an adjoint finite-difference test."
    )
    parser.add_argument(
        "--mesh-preset",
        choices=sorted(MESH_PRESETS.keys()),
        default="small",
        help="Mesh pair to use for the check. Default: small.",
    )
    parser.add_argument(
        "--beam-mdpa",
        type=Path,
        default=None,
        help="Optional beam mdpa path overriding --mesh-preset.",
    )
    parser.add_argument(
        "--surface-mdpa",
        type=Path,
        default=None,
        help="Optional surface mdpa path overriding --mesh-preset.",
    )
    parser.add_argument(
        "--case-limit",
        type=int,
        default=None,
        help="Run only the first N cases. Useful for quick timing probes on large meshes.",
    )
    return parser.parse_args()


def main():
    global BEAM_MDPA_PATH
    global SURFACE_MDPA_PATH

    args = parse_arguments()
    preset_beam_mdpa_path, preset_surface_mdpa_path = MESH_PRESETS[args.mesh_preset]
    BEAM_MDPA_PATH = args.beam_mdpa or preset_beam_mdpa_path
    SURFACE_MDPA_PATH = args.surface_mdpa or preset_surface_mdpa_path

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    results = []
    test_cases = TEST_CASES[:args.case_limit] if args.case_limit is not None else TEST_CASES
    for case_index, (kinematics_mode, parameter_value, force_mode, perturbation_mode) in enumerate(test_cases, start=1):
        case_name = output_folder_name(kinematics_mode, parameter_value)
        output_path = OUTPUT_ROOT / case_name
        with redirect_stdout_to_file(output_path / "console_log.txt"):
            result = run_case(
                case_index,
                kinematics_mode,
                parameter_value,
                force_mode,
                perturbation_mode,
                output_path,
            )
        results.append(result)
        print(
            "Case completed",
            case_name,
            "status=",
            "PASS" if result["passed"] else "FAIL",
            "relative_error=",
            f"{result['relative_error']:.6e}",
            "inverse_map_s=",
            f"{result['inverse_map_time']:.6f}",
            "total_s=",
            f"{result['total_time']:.6f}",
        )

    summary_path = OUTPUT_ROOT / "summary.txt"
    with open(summary_path, "w") as summary_file:
        summary_file.write(
            f"# run_started: {datetime.datetime.now().isoformat(timespec='seconds')}\n"
        )
        summary_file.write(f"# mesh_preset: {args.mesh_preset}\n")
        summary_file.write(f"# beam_mdpa: {BEAM_MDPA_PATH}\n")
        summary_file.write(f"# surface_mdpa: {SURFACE_MDPA_PATH}\n")
        summary_file.write(f"# epsilon: {FINITE_DIFFERENCE_EPSILON:.16e}\n")
        summary_file.write(f"# absolute_tolerance: {ABSOLUTE_TOLERANCE:.16e}\n")
        summary_file.write(f"# relative_tolerance: {RELATIVE_TOLERANCE:.16e}\n")
        summary_file.write(
            "case_index status kinematics_mode parameter_value force_mode perturbation_mode "
            "beam_nodes beam_elements surface_nodes surface_elements surface_conditions "
            "surface_directional_work beam_generalized_work absolute_error relative_error "
            "setup_s forward_map_s inverse_map_s perturbed_forward_map_s total_s\n"
        )
        for result in results:
            summary_file.write(
                f"{result['case_index']} "
                f"{'PASS' if result['passed'] else 'FAIL'} "
                f"{result['kinematics_mode']} "
                f"{result['parameter_value']:.16e} "
                f"{result['force_mode']} "
                f"{result['perturbation_mode']} "
                f"{result['beam_nodes']} "
                f"{result['beam_elements']} "
                f"{result['surface_nodes']} "
                f"{result['surface_elements']} "
                f"{result['surface_conditions']} "
                f"{result['surface_work']:.16e} "
                f"{result['beam_work']:.16e} "
                f"{result['absolute_error']:.16e} "
                f"{result['relative_error']:.16e} "
                f"{result['setup_time']:.8f} "
                f"{result['forward_map_time']:.8f} "
                f"{result['inverse_map_time']:.8f} "
                f"{result['perturbed_forward_map_time']:.8f} "
                f"{result['total_time']:.8f}\n"
            )

        if results:
            summary_file.write("aggregate_timings_seconds\n")
            summary_file.write(
                "setup_s forward_map_s inverse_map_s perturbed_forward_map_s total_s\n"
            )
            summary_file.write(
                f"{sum(result['setup_time'] for result in results):.8f} "
                f"{sum(result['forward_map_time'] for result in results):.8f} "
                f"{sum(result['inverse_map_time'] for result in results):.8f} "
                f"{sum(result['perturbed_forward_map_time'] for result in results):.8f} "
                f"{sum(result['total_time'] for result in results):.8f}\n"
            )

    all_passed = all(result["passed"] for result in results)
    print("summary:", summary_path)
    print("overall_status:", "PASS" if all_passed else "FAIL")

    if not all_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
