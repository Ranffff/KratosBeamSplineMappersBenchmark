from pathlib import Path
from contextlib import contextmanager
import argparse
import ctypes
import datetime
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication
import KratosMultiphysics.StructuralMechanicsApplication


ROOT = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT / "Bending_Torsion_Beam_Test_Case" / "Bending_Torsion_Beam_Test_Case"
BEAM_MDPA_PATH = CASE_DIR / "beam_geometry" / "beam_geometry_coarse.mdpa"
OUTPUT_ROOT = (
    ROOT
    / "Bending_Torsion_Beam_Test_Case"
    / "TestCase_Output"
    / "Small_Rotation_Analytical_Comparison_2x20"
)
THESIS_FIGURE_DIR = (
    ROOT.parent
    / "BeamSplineDocumentation"
    / "Thesis"
    / "MasterThesis"
    / "images"
    / "verification"
    / "linear_mapper_comparison"
)

BEAM_SPLINE_COLOR = "#0072B2"
LINEAR_BEAM_COLOR = "#E69F00"
BAR_EDGE_COLOR = "#202020"
GRID_COLOR = "#D9D9D9"

BEAM_LENGTH = 10.0
CROSS_SECTION_SIZE = 1.0
SURFACE_CROSS_SECTION_DIVISIONS = 2
SURFACE_AXIAL_DIVISIONS = 20
LIBC = ctypes.CDLL(None)

# All rotation amplitudes remain below 0.05 rad. The transverse displacement
# fields satisfy theta_y=-dw/dx and theta_z=dv/dx.
SMALL_ROTATION_CASES = [
    {
        "name": "displacement_y",
        "tip_displacement_y": 0.10,
    },
    {
        "name": "displacement_z",
        "tip_displacement_z": -0.08,
    },
    {
        "name": "rotation_y",
        "tip_bending_rotation_y": -0.020,
    },
    {
        "name": "rotation_z",
        "tip_bending_rotation_z": 0.025,
    },
    {
        "name": "torsion_x",
        "tip_torsion": 0.030,
    },
    {
        "name": "combined_y_torsion",
        "tip_displacement_y": 0.08,
        "tip_bending_rotation_z": 0.015,
        "tip_torsion": 0.025,
    },
    {
        "name": "combined_biaxial",
        "tip_displacement_y": 0.07,
        "tip_displacement_z": -0.06,
        "tip_bending_rotation_y": -0.012,
        "tip_bending_rotation_z": 0.018,
    },
    {
        "name": "combined_3d",
        "tip_axial_displacement": 0.02,
        "tip_displacement_y": 0.08,
        "tip_displacement_z": 0.05,
        "tip_bending_rotation_y": -0.015,
        "tip_bending_rotation_z": 0.020,
        "tip_torsion": -0.025,
    },
]


def prescribed_beam_kinematics(x, test_case):
    """Return small-rotation Euler-Bernoulli beam-axis kinematics."""
    xi = x / BEAM_LENGTH

    axial_shape = xi * xi * (3.0 - 2.0 * xi)
    displacement_shape = 0.5 * xi * xi * (3.0 - xi)
    displacement_shape_derivative = 1.5 * xi * (2.0 - xi) / BEAM_LENGTH

    u = test_case.get("tip_axial_displacement", 0.0) * axial_shape
    v = test_case.get("tip_displacement_y", 0.0) * displacement_shape
    w = test_case.get("tip_displacement_z", 0.0) * displacement_shape
    dv_dx = (
        test_case.get("tip_displacement_y", 0.0)
        * displacement_shape_derivative
    )
    dw_dx = (
        test_case.get("tip_displacement_z", 0.0)
        * displacement_shape_derivative
    )

    tip_rotation_y = test_case.get("tip_bending_rotation_y", 0.0)
    tip_rotation_z = test_case.get("tip_bending_rotation_z", 0.0)
    v += 0.5 * tip_rotation_z * BEAM_LENGTH * xi * xi
    w -= 0.5 * tip_rotation_y * BEAM_LENGTH * xi * xi
    dv_dx += tip_rotation_z * xi
    dw_dx -= tip_rotation_y * xi

    displacement = (u, v, w)
    rotation = (
        test_case.get("tip_torsion", 0.0) * xi,
        -dw_dx,
        dv_dx,
    )
    return displacement, rotation


def analytical_surface_displacement(node, test_case):
    """Linear rigid-cross-section displacement: u_s=u_b+theta cross r."""
    beam_displacement, rotation = prescribed_beam_kinematics(node.X0, test_case)
    theta_x, theta_y, theta_z = rotation
    y = node.Y0
    z = node.Z0
    rotational_offset = (
        theta_y * z - theta_z * y,
        -theta_x * z,
        theta_x * y,
    )
    return tuple(
        beam_displacement[i] + rotational_offset[i] for i in range(3)
    )


def create_structure_solution(name, test_case):
    model = KM.Model()
    structure = model.CreateModelPart(name)
    structure.ProcessInfo[KM.DOMAIN_SIZE] = 3
    structure.ProcessInfo[KM.TIME] = 0.0
    structure.ProcessInfo[KM.DELTA_TIME] = 1.0
    structure.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    structure.AddNodalSolutionStepVariable(KM.ROTATION)

    KM.ModelPartIO(str(BEAM_MDPA_PATH.with_suffix(""))).ReadModelPart(structure)
    beam = structure.GetSubModelPart("Parts_Beam_beam")
    for node in beam.Nodes:
        displacement, rotation = prescribed_beam_kinematics(node.X0, test_case)
        node.SetSolutionStepValue(KM.DISPLACEMENT, list(displacement))
        node.SetSolutionStepValue(KM.ROTATION, list(rotation))
    return beam


def create_surface_model_part(name):
    model = KM.Model()
    surface_root = model.CreateModelPart(name)
    surface_root.ProcessInfo[KM.DOMAIN_SIZE] = 3
    surface_root.ProcessInfo[KM.TIME] = 0.0
    surface_root.ProcessInfo[KM.DELTA_TIME] = 1.0
    surface_root.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    surface_root.AddNodalSolutionStepVariable(KM.ROTATION)
    properties = surface_root.CreateNewProperties(0)

    half_size = 0.5 * CROSS_SECTION_SIZE
    perimeter = []
    for i in range(SURFACE_CROSS_SECTION_DIVISIONS):
        fraction = i / SURFACE_CROSS_SECTION_DIVISIONS
        perimeter.append((-half_size + CROSS_SECTION_SIZE * fraction, half_size))
    for i in range(SURFACE_CROSS_SECTION_DIVISIONS):
        fraction = i / SURFACE_CROSS_SECTION_DIVISIONS
        perimeter.append((half_size, half_size - CROSS_SECTION_SIZE * fraction))
    for i in range(SURFACE_CROSS_SECTION_DIVISIONS):
        fraction = i / SURFACE_CROSS_SECTION_DIVISIONS
        perimeter.append((half_size - CROSS_SECTION_SIZE * fraction, -half_size))
    for i in range(SURFACE_CROSS_SECTION_DIVISIONS):
        fraction = i / SURFACE_CROSS_SECTION_DIVISIONS
        perimeter.append((-half_size, -half_size + CROSS_SECTION_SIZE * fraction))

    nodes_per_station = len(perimeter)
    node_ids = []
    for axial_index in range(SURFACE_AXIAL_DIVISIONS + 1):
        x = BEAM_LENGTH * axial_index / SURFACE_AXIAL_DIVISIONS
        for perimeter_index, (y, z) in enumerate(perimeter):
            node_id = axial_index * nodes_per_station + perimeter_index + 1
            surface_root.CreateNewNode(node_id, x, y, z)
            node_ids.append(node_id)

    element_ids = []
    element_id = 1
    for axial_index in range(SURFACE_AXIAL_DIVISIONS):
        station_start = axial_index * nodes_per_station
        next_station_start = (axial_index + 1) * nodes_per_station
        for perimeter_index in range(nodes_per_station):
            next_perimeter_index = (perimeter_index + 1) % nodes_per_station
            connectivity = [
                station_start + perimeter_index + 1,
                next_station_start + perimeter_index + 1,
                next_station_start + next_perimeter_index + 1,
                station_start + next_perimeter_index + 1,
            ]
            surface_root.CreateNewElement(
                "ShellThinElementCorotational3D4N",
                element_id,
                connectivity,
                properties,
            )
            element_ids.append(element_id)
            element_id += 1

    surface = surface_root.CreateSubModelPart("Parts_Shell_wet_surface")
    surface.AddNodes(node_ids)
    surface.AddElements(element_ids)
    return surface


def create_analytical_surface(name, test_case):
    surface = create_surface_model_part(name)
    for node in surface.Nodes:
        node.SetSolutionStepValue(
            KM.DISPLACEMENT,
            list(analytical_surface_displacement(node, test_case)),
        )
    return surface


def run_mapping(origin, destination, mapper_type):
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
    # Co-rotational mapper intentionally disabled for this small-rotation study.
    # elif mapper_type == "beam_mapper_corotation":
    #     settings = KM.Parameters("""{
    #         "mapper_type" : "beam_mapper",
    #         "use_corotation" : true,
    #         "search_settings" : {
    #             "max_num_search_iterations" : 30,
    #             "search_radius" : 3.0
    #         },
    #         "echo_level" : 0
    #     }""")
    else:
        raise RuntimeError(f"Unsupported mapper type: {mapper_type}")

    mapper = KM.MapperFactory.CreateMapper(origin, destination, settings)
    mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())


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


def write_vtk(model_part, out_name, output_path, nodal_solution_step_variables):
    vtk_params = KM.Parameters(r'''{
        "model_part_name": "PLEASE_SPECIFY_MODEL_PART_NAME",
        "custom_name_prefix": "",
        "custom_name_postfix": "",
        "output_control_type": "step",
        "output_interval": 1.0,
        "output_path": "VTK_Output",
        "save_output_files_in_folder": true,
        "entity_type": "automatic",
        "file_format": "ascii",
        "output_precision": 7,
        "output_sub_model_parts": false,
        "write_deformed_configuration": true,
        "write_ids": true,
        "nodal_solution_step_data_variables": [],
        "nodal_data_value_variables": [],
        "element_data_value_variables": [],
        "condition_data_value_variables": [],
        "gauss_point_variables_in_elements": [],
        "gauss_point_variables_extrapolated_to_nodes": []
    }''')

    output_path.mkdir(parents=True, exist_ok=True)
    vtk_params["model_part_name"].SetString(model_part.Name)
    vtk_params["custom_name_prefix"].SetString(f"{out_name}_")
    vtk_params["output_path"].SetString(str(output_path))
    for variable_name in nodal_solution_step_variables:
        vtk_params["nodal_solution_step_data_variables"].Append(variable_name)
    KM.VtkOutput(model_part, vtk_params).PrintOutput()


def vector_to_tuple(value):
    return (float(value[0]), float(value[1]), float(value[2]))


def calculate_error_metrics(mapped_surface, analytical_surface):
    error_squared = 0.0
    reference_squared = 0.0
    maximum_error = 0.0
    maximum_error_node = None

    for node in mapped_surface.Nodes:
        mapped = vector_to_tuple(node.GetSolutionStepValue(KM.DISPLACEMENT))
        reference = vector_to_tuple(
            analytical_surface.GetNode(node.Id).GetSolutionStepValue(KM.DISPLACEMENT)
        )
        nodal_error_squared = sum(
            (mapped[i] - reference[i]) ** 2 for i in range(3)
        )
        nodal_error = math.sqrt(nodal_error_squared)
        error_squared += nodal_error_squared
        reference_squared += sum(component * component for component in reference)
        if nodal_error > maximum_error:
            maximum_error = nodal_error
            maximum_error_node = node.Id

    number_of_components = 3 * mapped_surface.NumberOfNodes()
    rmse = math.sqrt(error_squared / number_of_components)
    reference_rmse = math.sqrt(reference_squared / number_of_components)
    normalized_rmse_percent = (
        0.0 if reference_rmse <= 0.0 else 100.0 * rmse / reference_rmse
    )
    return {
        "rmse": rmse,
        "reference_rmse": reference_rmse,
        "normalized_rmse_percent": normalized_rmse_percent,
        "maximum_error": maximum_error,
        "maximum_error_node": maximum_error_node,
    }


def case_parameter_rows(test_case):
    parameter_names = (
        "tip_axial_displacement",
        "tip_displacement_y",
        "tip_displacement_z",
        "tip_bending_rotation_y",
        "tip_bending_rotation_z",
        "tip_torsion",
    )
    return [
        (parameter_name, test_case.get(parameter_name, 0.0))
        for parameter_name in parameter_names
    ]


def print_section(title):
    print("")
    print(f"[{title}]")


def print_key_values(rows):
    width = max(len(name) for name, _ in rows)
    for name, value in rows:
        print(f"  {name:<{width}} : {value}")


def run_case(test_case):
    output_path = OUTPUT_ROOT / test_case["name"]
    with redirect_stdout_to_file(output_path / "console_log.txt"):
        print("=" * 80)
        print("SMALL-ROTATION ANALYTICAL MAPPER COMPARISON")
        print("=" * 80)
        print_section("CASE CONFIGURATION")
        print_key_values(
            [
                ("test_case", test_case["name"]),
                *case_parameter_rows(test_case),
                ("analytical_reference", "u_surface = u_beam + theta cross r"),
                ("beam_mdpa", BEAM_MDPA_PATH),
                (
                    "surface_mesh",
                    f"{SURFACE_CROSS_SECTION_DIVISIONS}x"
                    f"{SURFACE_AXIAL_DIVISIONS}, programmatic",
                ),
                ("vtk_output_path", output_path),
            ]
        )

        beam = create_structure_solution(f"{test_case['name']}_beam", test_case)
        analytical_surface = create_analytical_surface(
            f"{test_case['name']}_analytical_surface", test_case
        )
        spline_surface = create_surface_model_part(
            f"{test_case['name']}_beam_spline_surface"
        )
        linear_surface = create_surface_model_part(
            f"{test_case['name']}_linear_beam_mapper_surface"
        )

        run_mapping(beam, spline_surface, "beam_spline_mapper")
        run_mapping(beam, linear_surface, "beam_mapper_linear")
        spline_metrics = calculate_error_metrics(
            spline_surface, analytical_surface
        )
        linear_metrics = calculate_error_metrics(
            linear_surface, analytical_surface
        )

        write_vtk(beam, "origin_beam", output_path, ["DISPLACEMENT", "ROTATION"])
        write_vtk(
            analytical_surface,
            "analytical_reference_surface",
            output_path,
            ["DISPLACEMENT"],
        )
        write_vtk(
            spline_surface,
            "beam_spline_mapped_surface",
            output_path,
            ["DISPLACEMENT"],
        )
        write_vtk(
            linear_surface,
            "linear_beam_mapper_mapped_surface",
            output_path,
            ["DISPLACEMENT"],
        )

        print_section("MESH")
        print_key_values(
            [
                ("beam", f"{beam.NumberOfNodes()} nodes, {beam.NumberOfElements()} elements"),
                (
                    "surface",
                    f"{analytical_surface.NumberOfNodes()} nodes, "
                    f"{analytical_surface.NumberOfElements()} elements",
                ),
            ]
        )
        print_section("ERROR SUMMARY")
        print(f"  {'mapper':<24} {'RMSE':>14} {'normalized RMSE':>18} {'maximum error':>16}")
        print(f"  {'-' * 24} {'-' * 14} {'-' * 18} {'-' * 16}")
        for mapper_name, metrics in (
            ("beam_spline_mapper", spline_metrics),
            ("linear_beam_mapper", linear_metrics),
        ):
            print(
                f"  {mapper_name:<24} "
                f"{metrics['rmse']:>14.8e} "
                f"{metrics['normalized_rmse_percent']:>17.8f}% "
                f"{metrics['maximum_error']:>16.8e}"
            )

    return {
        "test_case": test_case["name"],
        "parameters": dict(case_parameter_rows(test_case)),
        "surface_nodes": analytical_surface.NumberOfNodes(),
        "surface_elements": analytical_surface.NumberOfElements(),
        "spline": spline_metrics,
        "linear": linear_metrics,
    }


def write_summary(results):
    summary_path = OUTPUT_ROOT / "small_rotation_comparison_summary.txt"
    with open(summary_path, "a") as summary_file:
        summary_file.write("\n")
        summary_file.write(
            f"# run_started: {datetime.datetime.now().isoformat(timespec='seconds')}\n"
        )
        summary_file.write("# reference: analytical small-rotation rigid cross-section\n")
        summary_file.write("# relation: theta_y=-dw/dx, theta_z=dv/dx\n")
        summary_file.write("# co_rotation: disabled\n")
        summary_file.write(
            "test_case tip_axial_displacement tip_displacement_y "
            "tip_displacement_z tip_bending_rotation_y "
            "tip_bending_rotation_z tip_torsion surface_nodes surface_elements "
            "spline_rmse spline_normalized_rmse_percent spline_maximum_error "
            "linear_rmse linear_normalized_rmse_percent linear_maximum_error\n"
        )
        for result in results:
            parameters = result["parameters"]
            summary_file.write(
                f"{result['test_case']} "
                f"{parameters['tip_axial_displacement']:.16e} "
                f"{parameters['tip_displacement_y']:.16e} "
                f"{parameters['tip_displacement_z']:.16e} "
                f"{parameters['tip_bending_rotation_y']:.16e} "
                f"{parameters['tip_bending_rotation_z']:.16e} "
                f"{parameters['tip_torsion']:.16e} "
                f"{result['surface_nodes']} {result['surface_elements']} "
                f"{result['spline']['rmse']:.16e} "
                f"{result['spline']['normalized_rmse_percent']:.8f} "
                f"{result['spline']['maximum_error']:.16e} "
                f"{result['linear']['rmse']:.16e} "
                f"{result['linear']['normalized_rmse_percent']:.8f} "
                f"{result['linear']['maximum_error']:.16e}\n"
            )
    return summary_path


def read_latest_summary_results():
    """Read the latest complete result block without rerunning the cases."""
    summary_path = OUTPUT_ROOT / "small_rotation_comparison_summary.txt"
    if not summary_path.is_file():
        raise FileNotFoundError(f"No existing summary found at {summary_path}")

    header_prefix = "test_case tip_axial_displacement tip_displacement_y "
    lines = summary_path.read_text().splitlines()
    header_indices = [
        index for index, line in enumerate(lines) if line.startswith(header_prefix)
    ]
    if not header_indices:
        raise RuntimeError(f"No result block found in {summary_path}")

    results = []
    parameter_names = (
        "tip_axial_displacement",
        "tip_displacement_y",
        "tip_displacement_z",
        "tip_bending_rotation_y",
        "tip_bending_rotation_z",
        "tip_torsion",
    )
    for line in lines[header_indices[-1] + 1:]:
        if not line.strip():
            if results:
                break
            continue
        if line.startswith("#"):
            break

        values = line.split()
        parameters = {
            name: float(value)
            for name, value in zip(parameter_names, values[1:7])
        }
        results.append(
            {
                "test_case": values[0],
                "parameters": parameters,
                "surface_nodes": int(values[7]),
                "surface_elements": int(values[8]),
                "spline": {
                    "rmse": float(values[9]),
                    "normalized_rmse_percent": float(values[10]),
                    "maximum_error": float(values[11]),
                },
                "linear": {
                    "rmse": float(values[12]),
                    "normalized_rmse_percent": float(values[13]),
                    "maximum_error": float(values[14]),
                },
            }
        )

    if len(results) != len(SMALL_ROTATION_CASES):
        raise RuntimeError(
            f"Latest summary block contains {len(results)} cases; "
            f"expected {len(SMALL_ROTATION_CASES)}"
        )
    return results


def save_comparison_figure(figure, plot_path, filename):
    png_path = plot_path / filename
    svg_path = png_path.with_suffix(".svg")
    figure.savefig(png_path, dpi=300, bbox_inches="tight")
    figure.savefig(svg_path, bbox_inches="tight")

    THESIS_FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    thesis_filenames = {
        "rmse_by_case.png": "small_rotation_2x20__spline_vs_linear__rmse.png",
        "normalized_rmse_by_case.png": (
            "small_rotation_2x20__spline_vs_linear__normalized_rmse.png"
        ),
        "maximum_error_by_case.png": (
            "small_rotation_2x20__spline_vs_linear__max_nodal_error.png"
        ),
    }
    thesis_png_path = THESIS_FIGURE_DIR / thesis_filenames[filename]
    figure.savefig(thesis_png_path, dpi=300, bbox_inches="tight")
    figure.savefig(thesis_png_path.with_suffix(".svg"), bbox_inches="tight")


def create_plots(results):
    plot_path = OUTPUT_ROOT / "Summary_Plots"
    plot_path.mkdir(parents=True, exist_ok=True)
    case_labels = [result["test_case"].replace("_", " ") for result in results]
    case_positions = list(range(len(results)))
    bar_width = 0.36

    for metric_name, ylabel, filename in (
        ("rmse", r"$e_{\mathrm{RMSE}}$", "rmse_by_case.png"),
        (
            "normalized_rmse_percent",
            r"$e_{\mathrm{norm}}$ [\%]",
            "normalized_rmse_by_case.png",
        ),
        ("maximum_error", "Maximum nodal error", "maximum_error_by_case.png"),
    ):
        spline_values = [result["spline"][metric_name] for result in results]
        linear_values = [result["linear"][metric_name] for result in results]
        figure, axis = plt.subplots(figsize=(9.2, 4.8))
        axis.bar(
            [position - bar_width / 2 for position in case_positions],
            spline_values,
            bar_width,
            label="Beam-spline mapper",
            color=BEAM_SPLINE_COLOR,
            edgecolor=BAR_EDGE_COLOR,
            linewidth=0.8,
            hatch="///",
        )
        axis.bar(
            [position + bar_width / 2 for position in case_positions],
            linear_values,
            bar_width,
            label="Linear beam mapper",
            color=LINEAR_BEAM_COLOR,
            edgecolor=BAR_EDGE_COLOR,
            linewidth=0.8,
            hatch="xxx",
        )
        axis.set_xticks(case_positions)
        axis.set_xticklabels(case_labels, rotation=24, ha="right")
        axis.set_ylabel(ylabel)
        axis.grid(True, axis="y", color=GRID_COLOR, linewidth=0.6)
        axis.set_axisbelow(True)
        axis.legend()
        figure.tight_layout()
        save_comparison_figure(figure, plot_path, filename)
        plt.close(figure)


def validate_test_cases():
    names = [test_case["name"] for test_case in SMALL_ROTATION_CASES]
    if len(names) != len(set(names)):
        raise RuntimeError("Small-rotation test-case names must be unique")

    for test_case in SMALL_ROTATION_CASES:
        for sample_index in range(101):
            x = BEAM_LENGTH * sample_index / 100.0
            _, rotation = prescribed_beam_kinematics(x, test_case)
            if max(abs(component) for component in rotation) > 0.05 + 1.0e-12:
                raise RuntimeError(
                    f"Case {test_case['name']} exceeds the 0.05 rad "
                    "small-rotation limit"
                )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Compare beam-spline and linear beam mappers with an analytical "
            "small-rotation Euler-Bernoulli surface solution."
        )
    )
    parser.add_argument(
        "--case",
        choices=[test_case["name"] for test_case in SMALL_ROTATION_CASES],
        help="Run only the selected test case.",
    )
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Regenerate plots from the latest saved summary without mapping.",
    )
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    validate_test_cases()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    if arguments.plots_only:
        create_plots(read_latest_summary_results())
        print("plots regenerated:", OUTPUT_ROOT / "Summary_Plots")
        return

    selected_cases = [
        test_case
        for test_case in SMALL_ROTATION_CASES
        if arguments.case is None or test_case["name"] == arguments.case
    ]

    results = []
    for test_case in selected_cases:
        result = run_case(test_case)
        results.append(result)
        print(
            "Case completed",
            test_case["name"],
            "spline_normalized_rmse_percent=",
            result["spline"]["normalized_rmse_percent"],
            "linear_normalized_rmse_percent=",
            result["linear"]["normalized_rmse_percent"],
        )

    summary_path = write_summary(results)
    create_plots(results)
    print("summary:", summary_path)
    print("plots:", OUTPUT_ROOT / "Summary_Plots")


if __name__ == "__main__":
    main()
