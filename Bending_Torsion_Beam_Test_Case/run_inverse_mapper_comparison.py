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
BEAM_MDPA_PATH = CASE_DIR / "beam_geometry" / "beam_geometry.mdpa"
SURFACE_MDPA_PATH = CASE_DIR / "fluid_surface_mesh" / "fluid_surface_mesh_10x100.mdpa"
OUTPUT_ROOT = ROOT / "Bending_Torsion_Beam_Test_Case" / "TestCase_Output" / "InverseMap_Mapper_Comparison"
BEAM_LENGTH = 10.0
FINITE_DIFFERENCE_EPSILON = 1.0e-7
REFERENCE_MAPPER = "beam_mapper_corotation"
MAPPER_TYPES = (
    "beam_spline_mapper",
    "beam_mapper_linear",
    "beam_mapper_corotation",
)
LIBC = ctypes.CDLL(None)


TEST_CASES = [
    ("rotation", 0.25 * math.pi, "constant_y", "bending_y"),
    ("rotation", 0.50 * math.pi, "constant_y", "bending_y"),
    ("displacement", 1.0, "mixed", "mixed"),
    ("torsion", 0.30 * math.pi, "torsional_ring", "torsion"),
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

    if mode == "rotation":
        displacement = (0.0, 0.5 * parameter_value * BEAM_LENGTH * xi * xi, 0.0)
        rotation = (0.0, 0.0, parameter_value * xi)
        return displacement, rotation

    if mode == "displacement":
        displacement = (
            0.0,
            0.5 * parameter_value * (3.0 * xi * xi - xi * xi * xi),
            0.0,
        )
        rotation = (
            0.0,
            0.0,
            1.5 * parameter_value * (2.0 * xi - xi * xi) / BEAM_LENGTH,
        )
        return displacement, rotation

    if mode == "torsion":
        return (0.0, 0.0, 0.0), (parameter_value * xi, 0.0, 0.0)

    raise RuntimeError(f"Unsupported kinematics mode: {mode}")


def perturbation(mode, x):
    xi = x / BEAM_LENGTH

    if mode == "bending_y":
        return (0.0, 0.15 * xi * xi, 0.0), (0.0, 0.0, 0.03 * xi)

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

    if mode == "torsion":
        return (0.0, 0.0, 0.0), (0.04 * xi, 0.0, 0.0)

    raise RuntimeError(f"Unsupported perturbation mode: {mode}")


def surface_force(force_mode, node):
    x = node.X0
    y = node.Y0
    z = node.Z0

    if force_mode == "constant_y":
        return (0.0, 1.0, 0.0)

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


def create_mapper(origin, destination, mapper_type):
    if mapper_type == "beam_spline_mapper":
        settings = KM.Parameters("""{
            "mapper_type" : "beam_spline_mapper",
            "search_settings" : {
                "search_radius" : 3.0,
                "max_num_search_iterations" : 30
            },
            "local_coord_tolerance" : 0.25,
            "echo_level" : 0
        }""")
    elif mapper_type == "beam_mapper_linear":
        settings = KM.Parameters("""{
            "mapper_type" : "beam_mapper",
            "use_corotation" : false,
            "search_settings" : {
                "max_num_search_iterations" : 30,
                "search_radius" : 3.0
            },
            "echo_level" : 0
        }""")
    elif mapper_type == "beam_mapper_corotation":
        settings = KM.Parameters("""{
            "mapper_type" : "beam_mapper",
            "use_corotation" : true,
            "search_settings" : {
                "max_num_search_iterations" : 30,
                "search_radius" : 3.0
            },
            "echo_level" : 0
        }""")
    else:
        raise RuntimeError(f"Unsupported mapper type: {mapper_type}")

    return KM.MapperFactory.CreateMapper(origin, destination, settings)


def vector_to_tuple(value):
    return (float(value[0]), float(value[1]), float(value[2]))


def dot(lhs, rhs):
    return sum(lhs[i] * rhs[i] for i in range(3))


def norm_3d(value):
    return math.sqrt(sum(component * component for component in value))


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


def capture_beam_loads(beam):
    return {
        node.Id: (
            vector_to_tuple(node.GetSolutionStepValue(KM.FORCE)),
            vector_to_tuple(node.GetSolutionStepValue(KM.MOMENT)),
        )
        for node in beam.Nodes
    }


def run_mapper_case(mapper_type, kinematics_mode, parameter_value, force_mode, perturbation_mode):
    start_time = time.perf_counter()
    beam = create_structure_model_part(f"{mapper_type}_beam", kinematics_mode, parameter_value)
    surface = create_surface_model_part(f"{mapper_type}_surface")
    mapper = create_mapper(beam, surface, mapper_type)
    setup_time = time.perf_counter() - start_time

    forward_start = time.perf_counter()
    mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
    forward_time = time.perf_counter() - forward_start
    base_surface_displacements = capture_surface_displacements(surface)

    assign_surface_forces(surface, force_mode)
    inverse_start = time.perf_counter()
    mapper.InverseMap(KM.FORCE, KM.MOMENT, KM.REACTION, KM.Flags())
    inverse_time = time.perf_counter() - inverse_start
    beam_work = beam_generalized_work(beam, perturbation_mode)
    beam_loads = capture_beam_loads(beam)

    apply_perturbation(beam, perturbation_mode, FINITE_DIFFERENCE_EPSILON)
    perturbed_forward_start = time.perf_counter()
    mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
    perturbed_forward_time = time.perf_counter() - perturbed_forward_start
    surface_work = surface_directional_work(surface, base_surface_displacements, FINITE_DIFFERENCE_EPSILON)

    absolute_work_error = abs(surface_work - beam_work)
    work_reference = max(abs(surface_work), abs(beam_work), 1.0)
    relative_work_error = absolute_work_error / work_reference
    total_time = time.perf_counter() - start_time

    return {
        "mapper_type": mapper_type,
        "beam_nodes": beam.NumberOfNodes(),
        "beam_elements": beam.NumberOfElements(),
        "surface_nodes": surface.NumberOfNodes(),
        "surface_elements": surface.NumberOfElements(),
        "surface_conditions": surface.NumberOfConditions(),
        "surface_work": surface_work,
        "beam_work": beam_work,
        "absolute_work_error": absolute_work_error,
        "relative_work_error": relative_work_error,
        "beam_loads": beam_loads,
        "setup_time": setup_time,
        "forward_time": forward_time,
        "inverse_time": inverse_time,
        "perturbed_forward_time": perturbed_forward_time,
        "total_time": total_time,
    }


def compare_loads_to_reference(loads, reference_loads):
    max_force_diff = 0.0
    max_moment_diff = 0.0
    max_combined_diff = 0.0
    max_force_reference = 0.0
    max_moment_reference = 0.0
    max_combined_reference = 0.0
    max_combined_node = None

    for node_id, (force, moment) in loads.items():
        reference_force, reference_moment = reference_loads[node_id]
        force_diff = tuple(force[i] - reference_force[i] for i in range(3))
        moment_diff = tuple(moment[i] - reference_moment[i] for i in range(3))
        force_diff_norm = norm_3d(force_diff)
        moment_diff_norm = norm_3d(moment_diff)
        combined_diff_norm = math.sqrt(force_diff_norm * force_diff_norm + moment_diff_norm * moment_diff_norm)
        force_reference_norm = norm_3d(reference_force)
        moment_reference_norm = norm_3d(reference_moment)
        combined_reference_norm = math.sqrt(force_reference_norm * force_reference_norm + moment_reference_norm * moment_reference_norm)

        max_force_diff = max(max_force_diff, force_diff_norm)
        max_moment_diff = max(max_moment_diff, moment_diff_norm)
        max_force_reference = max(max_force_reference, force_reference_norm)
        max_moment_reference = max(max_moment_reference, moment_reference_norm)
        max_combined_reference = max(max_combined_reference, combined_reference_norm)
        if combined_diff_norm > max_combined_diff:
            max_combined_diff = combined_diff_norm
            max_combined_node = node_id

    return {
        "max_force_diff_to_reference": max_force_diff,
        "max_moment_diff_to_reference": max_moment_diff,
        "max_combined_load_diff_to_reference": max_combined_diff,
        "max_combined_load_diff_node": max_combined_node,
        "max_force_reference": max_force_reference,
        "max_moment_reference": max_moment_reference,
        "max_combined_load_reference": max_combined_reference,
        "combined_load_diff_percent": percentage(max_combined_diff, max_combined_reference),
    }


def percentage(error_norm, reference_norm):
    if reference_norm <= 0.0:
        return 0.0
    return 100.0 * error_norm / reference_norm


def run_case(case_index, kinematics_mode, parameter_value, force_mode, perturbation_mode, output_path):
    print("INVERSE MAP BEAM MAPPER COMPARISON")
    print("case_index:", case_index)
    print("kinematics_mode:", kinematics_mode)
    print("parameter_value:", parameter_value)
    print("force_mode:", force_mode)
    print("perturbation_mode:", perturbation_mode)
    print("epsilon:", FINITE_DIFFERENCE_EPSILON)
    print("reference_mapper:", REFERENCE_MAPPER)
    print("beam mdpa:", BEAM_MDPA_PATH)
    print("surface mdpa:", SURFACE_MDPA_PATH)
    print("output_path:", output_path)
    print("")

    mapper_results = {}
    for mapper_type in MAPPER_TYPES:
        result = run_mapper_case(
            mapper_type,
            kinematics_mode,
            parameter_value,
            force_mode,
            perturbation_mode,
        )
        mapper_results[mapper_type] = result
        print(
            "mapper:",
            mapper_type,
            "surface_work:",
            f"{result['surface_work']:.16e}",
            "beam_work:",
            f"{result['beam_work']:.16e}",
            "relative_work_error:",
            f"{result['relative_work_error']:.16e}",
            "inverse_s:",
            f"{result['inverse_time']:.6f}",
        )

    reference_loads = mapper_results[REFERENCE_MAPPER]["beam_loads"]
    case_results = []
    print("")
    print("load_comparison_reference:", REFERENCE_MAPPER)
    for mapper_type, result in mapper_results.items():
        load_metrics = compare_loads_to_reference(result["beam_loads"], reference_loads)
        beam_work_diff = abs(result["beam_work"] - mapper_results[REFERENCE_MAPPER]["beam_work"])
        beam_work_diff_percent = percentage(
            beam_work_diff,
            abs(mapper_results[REFERENCE_MAPPER]["beam_work"]),
        )
        result.update(load_metrics)
        result["beam_work_diff_to_reference"] = beam_work_diff
        result["beam_work_diff_to_reference_percent"] = beam_work_diff_percent
        case_results.append(result)

        print(
            mapper_type,
            "combined_load_diff_percent:",
            f"{load_metrics['combined_load_diff_percent']:.8f}",
            "max_combined_node:",
            load_metrics["max_combined_load_diff_node"],
            "beam_work_diff_percent:",
            f"{beam_work_diff_percent:.8f}",
        )

    return case_results


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
    if kinematics_mode == "rotation":
        return f"theta={format_pi_multiple(parameter_value)}"
    if kinematics_mode == "displacement":
        return f"u={parameter_value:g}"
    if kinematics_mode == "torsion":
        return f"torsion={format_pi_multiple(parameter_value)}"
    raise RuntimeError(f"Unsupported kinematics mode: {kinematics_mode}")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Compare InverseMap accuracy for BeamSplineMapper and BeamMapper variants."
    )
    parser.add_argument(
        "--mesh-preset",
        choices=sorted(MESH_PRESETS.keys()),
        default="small",
        help="Mesh pair to use for the comparison. Default: small.",
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
        help="Run only the first N cases.",
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

    all_results = []
    test_cases = TEST_CASES[:args.case_limit] if args.case_limit is not None else TEST_CASES
    for case_index, (kinematics_mode, parameter_value, force_mode, perturbation_mode) in enumerate(test_cases, start=1):
        case_name = output_folder_name(kinematics_mode, parameter_value)
        output_path = OUTPUT_ROOT / case_name
        with redirect_stdout_to_file(output_path / "console_log.txt"):
            case_results = run_case(
                case_index,
                kinematics_mode,
                parameter_value,
                force_mode,
                perturbation_mode,
                output_path,
            )
        all_results.extend(
            {
                **result,
                "case_index": case_index,
                "kinematics_mode": kinematics_mode,
                "parameter_value": parameter_value,
                "force_mode": force_mode,
                "perturbation_mode": perturbation_mode,
            }
            for result in case_results
        )
        spline_result = next(result for result in case_results if result["mapper_type"] == "beam_spline_mapper")
        print(
            "Case completed",
            case_name,
            "spline_relative_work_error=",
            f"{spline_result['relative_work_error']:.6e}",
            "spline_load_diff_to_corotation_percent=",
            f"{spline_result['combined_load_diff_percent']:.6f}",
            "spline_inverse_s=",
            f"{spline_result['inverse_time']:.6f}",
        )

    summary_path = OUTPUT_ROOT / "summary.txt"
    with open(summary_path, "w") as summary_file:
        summary_file.write(
            f"# run_started: {datetime.datetime.now().isoformat(timespec='seconds')}\n"
        )
        summary_file.write(f"# mesh_preset: {args.mesh_preset}\n")
        summary_file.write(f"# beam_mdpa: {BEAM_MDPA_PATH}\n")
        summary_file.write(f"# surface_mdpa: {SURFACE_MDPA_PATH}\n")
        summary_file.write(f"# reference_mapper: {REFERENCE_MAPPER}\n")
        summary_file.write(f"# epsilon: {FINITE_DIFFERENCE_EPSILON:.16e}\n")
        summary_file.write(
            "case_index mapper_type kinematics_mode parameter_value force_mode perturbation_mode "
            "beam_nodes beam_elements surface_nodes surface_elements surface_conditions "
            "surface_work beam_work absolute_work_error relative_work_error "
            "beam_work_diff_to_reference beam_work_diff_to_reference_percent "
            "max_force_diff_to_reference max_moment_diff_to_reference "
            "max_combined_load_diff_to_reference max_combined_load_diff_node combined_load_diff_percent "
            "setup_s forward_map_s inverse_map_s perturbed_forward_map_s total_s\n"
        )
        for result in all_results:
            summary_file.write(
                f"{result['case_index']} "
                f"{result['mapper_type']} "
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
                f"{result['absolute_work_error']:.16e} "
                f"{result['relative_work_error']:.16e} "
                f"{result['beam_work_diff_to_reference']:.16e} "
                f"{result['beam_work_diff_to_reference_percent']:.8f} "
                f"{result['max_force_diff_to_reference']:.16e} "
                f"{result['max_moment_diff_to_reference']:.16e} "
                f"{result['max_combined_load_diff_to_reference']:.16e} "
                f"{result['max_combined_load_diff_node']} "
                f"{result['combined_load_diff_percent']:.8f} "
                f"{result['setup_time']:.8f} "
                f"{result['forward_time']:.8f} "
                f"{result['inverse_time']:.8f} "
                f"{result['perturbed_forward_time']:.8f} "
                f"{result['total_time']:.8f}\n"
            )

    print("summary:", summary_path)


if __name__ == "__main__":
    main()
