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
STUDY_BEAM_MDPA_PATH = CASE_DIR / "beam_geometry" / "beam_geometry_coarse.mdpa"
OUTPUT_ROOT = (
    ROOT
    / "Bending_Torsion_Beam_Test_Case"
    / "TestCase_Output"
    / "Fluid_Mesh_Refinement"
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

BEAM_LENGTH = 10.0
CROSS_SECTION_SIZE = 1.0

# Colour-blind-friendly thesis palette. The Beam-spline data keep the same
# visual identity in every plot, with orange markers used only as emphasis.
BEAM_SPLINE_COLOR = "#0072B2"
DATA_MARKER_COLOR = "#E69F00"
GRID_COLOR = "#D9D9D9"
REFERENCE_LINE_COLOR = "#666666"

# Euler-Bernoulli-compatible prescribed kinematics. Each case combines
# axial displacement, bending in both transverse directions, and torsion.
# The cross-section rotations are always theta_y=-dw/dx, theta_z=dv/dx,
# while the twist follows Wang Eq. (5.31), theta_x=beta*x/L.
KINEMATIC_TEST_CASES = [
    {
        "name": "biaxial_bending_twist",
        "tip_axial_displacement": 0.20,
        "cubic_tip_displacement_y": 0.80,
        "cubic_tip_displacement_z": -0.60,
        "tip_bending_rotation_y": -0.25 * math.pi,
        "tip_bending_rotation_z": 0.35 * math.pi,
        "tip_torsion": 1.00 * math.pi,
    },
    {
        "name": "asymmetric_large_torsion",
        "tip_axial_displacement": 0.40,
        "cubic_tip_displacement_y": 1.50,
        "cubic_tip_displacement_z": 0.90,
        "tip_bending_rotation_y": -0.30 * math.pi,
        "tip_bending_rotation_z": 0.50 * math.pi,
        "tip_torsion": 2.00 * math.pi,
    },
    {
        "name": "reverse_bending_large_torsion",
        "tip_axial_displacement": -0.25,
        "cubic_tip_displacement_y": -1.10,
        "cubic_tip_displacement_z": 1.30,
        "tip_bending_rotation_y": 0.40 * math.pi,
        "tip_bending_rotation_z": -0.45 * math.pi,
        "tip_torsion": -1.50 * math.pi,
    },
]

# Each refinement halves both the axial and cross-section element sizes.
FLUID_MESH_LEVELS = [
    {"cross_section_divisions": 1, "axial_divisions": 10},
    {"cross_section_divisions": 2, "axial_divisions": 20},
    {"cross_section_divisions": 4, "axial_divisions": 40},
    {"cross_section_divisions": 8, "axial_divisions": 80},
    {"cross_section_divisions": 16, "axial_divisions": 160},
]

REFERENCE_CROSS_SECTION_DIVISIONS = 16
REFERENCE_AXIAL_DIVISIONS = 160
REFERENCE_BEAM_DIVISIONS = 320
FORCE_REBUILD_REFERENCE = False
LIBC = ctypes.CDLL(None)


def prescribed_beam_kinematics(x, test_case):
    """Return Euler-Bernoulli-compatible 3D beam-axis kinematics."""
    xi = x / BEAM_LENGTH

    # Smooth axial extension with zero displacement and slope at the clamp.
    axial_shape = xi * xi * (3.0 - 2.0 * xi)
    u = test_case["tip_axial_displacement"] * axial_shape

    # Cubic cantilever displacement fields with zero value and slope at x=0.
    displacement_shape = 0.5 * xi * xi * (3.0 - xi)
    displacement_shape_derivative = 1.5 * xi * (2.0 - xi) / BEAM_LENGTH
    v = test_case["cubic_tip_displacement_y"] * displacement_shape
    w = test_case["cubic_tip_displacement_z"] * displacement_shape
    dv_dx = test_case["cubic_tip_displacement_y"] * displacement_shape_derivative
    dw_dx = test_case["cubic_tip_displacement_z"] * displacement_shape_derivative

    # Add constant-curvature bending contributions in both transverse planes.
    v += (
        0.5
        * test_case["tip_bending_rotation_z"]
        * BEAM_LENGTH
        * xi
        * xi
    )
    w += (
        -0.5
        * test_case["tip_bending_rotation_y"]
        * BEAM_LENGTH
        * xi
        * xi
    )
    dv_dx += test_case["tip_bending_rotation_z"] * xi
    dw_dx -= test_case["tip_bending_rotation_y"] * xi

    displacement = (u, v, w)
    rotation = (
        test_case["tip_torsion"] * xi,
        -dw_dx,
        dv_dx,
    )
    return displacement, rotation


def create_structure_solution(name, beam_mdpa_path, test_case):
    model = KM.Model()
    structure = model.CreateModelPart(name)
    structure.ProcessInfo[KM.DOMAIN_SIZE] = 3
    structure.ProcessInfo[KM.TIME] = 0.0
    structure.ProcessInfo[KM.DELTA_TIME] = 1.0
    structure.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    structure.AddNodalSolutionStepVariable(KM.ROTATION)

    KM.ModelPartIO(str(beam_mdpa_path.with_suffix(""))).ReadModelPart(structure)
    beam = structure.GetSubModelPart("Parts_Beam_beam")

    for node in beam.Nodes:
        displacement, rotation = prescribed_beam_kinematics(node.X0, test_case)
        node.SetSolutionStepValue(KM.DISPLACEMENT, list(displacement))
        node.SetSolutionStepValue(KM.ROTATION, list(rotation))

    return beam


def create_reference_structure_solution(name, test_case):
    model = KM.Model()
    structure = model.CreateModelPart(name)
    structure.ProcessInfo[KM.DOMAIN_SIZE] = 3
    structure.ProcessInfo[KM.TIME] = 0.0
    structure.ProcessInfo[KM.DELTA_TIME] = 1.0
    structure.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    structure.AddNodalSolutionStepVariable(KM.ROTATION)
    properties = structure.CreateNewProperties(0)

    node_ids = []
    for axial_index in range(REFERENCE_BEAM_DIVISIONS + 1):
        node_id = axial_index + 1
        x = BEAM_LENGTH * axial_index / REFERENCE_BEAM_DIVISIONS
        node = structure.CreateNewNode(node_id, x, 0.0, 0.0)
        displacement, rotation = prescribed_beam_kinematics(x, test_case)
        node.SetSolutionStepValue(KM.DISPLACEMENT, list(displacement))
        node.SetSolutionStepValue(KM.ROTATION, list(rotation))
        node_ids.append(node_id)

    element_ids = []
    for axial_index in range(REFERENCE_BEAM_DIVISIONS):
        element_id = axial_index + 1
        structure.CreateNewElement(
            "CrBeamElement3D2N",
            element_id,
            [axial_index + 1, axial_index + 2],
            properties,
        )
        element_ids.append(element_id)

    beam = structure.CreateSubModelPart("Parts_Beam_beam")
    beam.AddNodes(node_ids)
    beam.AddElements(element_ids)
    return beam


def create_reference_surface_model_part(name):
    return create_structured_surface_model_part(
        name,
        REFERENCE_CROSS_SECTION_DIVISIONS,
        REFERENCE_AXIAL_DIVISIONS,
    )


def perimeter_coordinates(cross_section_divisions):
    half_size = 0.5 * CROSS_SECTION_SIZE
    coordinates = []

    # Follow the ordering used by the existing square-tube surface meshes:
    # top, right, bottom, and left faces, without duplicate corner nodes.
    for i in range(cross_section_divisions):
        fraction = i / cross_section_divisions
        coordinates.append((-half_size + CROSS_SECTION_SIZE * fraction, half_size))
    for i in range(cross_section_divisions):
        fraction = i / cross_section_divisions
        coordinates.append((half_size, half_size - CROSS_SECTION_SIZE * fraction))
    for i in range(cross_section_divisions):
        fraction = i / cross_section_divisions
        coordinates.append((half_size - CROSS_SECTION_SIZE * fraction, -half_size))
    for i in range(cross_section_divisions):
        fraction = i / cross_section_divisions
        coordinates.append((-half_size, -half_size + CROSS_SECTION_SIZE * fraction))

    return coordinates


def create_structured_surface_model_part(name, cross_section_divisions, axial_divisions):
    model = KM.Model()
    surface_root = model.CreateModelPart(name)
    surface_root.ProcessInfo[KM.DOMAIN_SIZE] = 3
    surface_root.ProcessInfo[KM.TIME] = 0.0
    surface_root.ProcessInfo[KM.DELTA_TIME] = 1.0
    surface_root.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    surface_root.AddNodalSolutionStepVariable(KM.ROTATION)
    properties = surface_root.CreateNewProperties(0)

    perimeter = perimeter_coordinates(cross_section_divisions)
    nodes_per_station = len(perimeter)
    node_ids = []

    for axial_index in range(axial_divisions + 1):
        x = BEAM_LENGTH * axial_index / axial_divisions
        for perimeter_index, (y, z) in enumerate(perimeter):
            node_id = axial_index * nodes_per_station + perimeter_index + 1
            surface_root.CreateNewNode(node_id, x, y, z)
            node_ids.append(node_id)

    element_ids = []
    element_id = 1
    for axial_index in range(axial_divisions):
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


def coordinate_to_perimeter_parameter(y, z):
    half_size = 0.5 * CROSS_SECTION_SIZE
    tolerance = 1.0e-10

    if math.isclose(z, half_size, abs_tol=tolerance):
        return (y + half_size) / CROSS_SECTION_SIZE
    if math.isclose(y, half_size, abs_tol=tolerance):
        return 1.0 + (half_size - z) / CROSS_SECTION_SIZE
    if math.isclose(z, -half_size, abs_tol=tolerance):
        return 2.0 + (half_size - y) / CROSS_SECTION_SIZE
    if math.isclose(y, -half_size, abs_tol=tolerance):
        return 3.0 + (z + half_size) / CROSS_SECTION_SIZE

    raise RuntimeError(f"Node ({y}, {z}) is not on the square cross-section boundary")


def capture_reference_field(reference_surface):
    values = {}
    perimeter_node_count = 4 * REFERENCE_CROSS_SECTION_DIVISIONS

    for node in reference_surface.Nodes:
        axial_index = int(round(node.X0 * REFERENCE_AXIAL_DIVISIONS / BEAM_LENGTH))
        perimeter_parameter = coordinate_to_perimeter_parameter(node.Y0, node.Z0)
        perimeter_index = int(
            round(perimeter_parameter * REFERENCE_CROSS_SECTION_DIVISIONS)
        ) % perimeter_node_count
        key = (axial_index, perimeter_index)
        if key in values:
            raise RuntimeError(f"Duplicate reference-grid point detected: {key}")
        values[key] = vector_to_tuple(node.GetSolutionStepValue(KM.DISPLACEMENT))

    expected_size = (
        (REFERENCE_AXIAL_DIVISIONS + 1)
        * 4
        * REFERENCE_CROSS_SECTION_DIVISIONS
    )
    if len(values) != expected_size:
        raise RuntimeError(
            f"Expected {expected_size} reference nodes, captured {len(values)}"
        )
    return values


def reference_vtk_path(output_path):
    return (
        output_path
        / "reference_corotation_mapped_surface_"
        "Parts_Shell_wet_surface_0_0.vtk"
    )


def load_reference_field_from_vtk(vtk_path):
    if FORCE_REBUILD_REFERENCE or not vtk_path.is_file():
        return None

    try:
        with open(vtk_path, "r") as vtk_file:
            lines = iter(vtk_file)
            points_header = next(
                line.strip() for line in lines if line.startswith("POINTS ")
            )
            number_of_points = int(points_header.split()[1])
            coordinates = [
                tuple(float(component) for component in next(lines).split())
                for _ in range(number_of_points)
            ]
            displacement_header = next(
                line.strip()
                for line in lines
                if line.startswith("DISPLACEMENT ")
            )
            displacement_parts = displacement_header.split()
            number_of_components = int(displacement_parts[1])
            number_of_displacements = int(displacement_parts[2])
            displacements = [
                tuple(float(component) for component in next(lines).split())
                for _ in range(number_of_displacements)
            ]
    except (OSError, StopIteration, ValueError):
        return None

    expected_size = (
        (REFERENCE_AXIAL_DIVISIONS + 1)
        * 4
        * REFERENCE_CROSS_SECTION_DIVISIONS
    )
    if (
        number_of_points != expected_size
        or number_of_components != 3
        or number_of_displacements != expected_size
    ):
        return None

    reference_field = {}
    perimeter_node_count = 4 * REFERENCE_CROSS_SECTION_DIVISIONS
    for value_index, displacement in enumerate(displacements):
        axial_index = value_index // perimeter_node_count
        perimeter_index = value_index % perimeter_node_count
        reference_field[(axial_index, perimeter_index)] = displacement

    if len(reference_field) != expected_size:
        return None
    return reference_field


def interpolate_reference_displacement(reference_field, x, y, z):
    axial_coordinate = x * REFERENCE_AXIAL_DIVISIONS / BEAM_LENGTH
    axial_index_0 = min(int(math.floor(axial_coordinate)), REFERENCE_AXIAL_DIVISIONS)
    axial_index_1 = min(axial_index_0 + 1, REFERENCE_AXIAL_DIVISIONS)
    axial_fraction = axial_coordinate - axial_index_0

    perimeter_node_count = 4 * REFERENCE_CROSS_SECTION_DIVISIONS
    perimeter_coordinate = (
        coordinate_to_perimeter_parameter(y, z)
        * REFERENCE_CROSS_SECTION_DIVISIONS
    )
    perimeter_index_0_unwrapped = int(math.floor(perimeter_coordinate))
    perimeter_index_0 = perimeter_index_0_unwrapped % perimeter_node_count
    perimeter_index_1 = (perimeter_index_0 + 1) % perimeter_node_count
    perimeter_fraction = perimeter_coordinate - perimeter_index_0_unwrapped

    value_00 = reference_field[(axial_index_0, perimeter_index_0)]
    value_01 = reference_field[(axial_index_0, perimeter_index_1)]
    value_10 = reference_field[(axial_index_1, perimeter_index_0)]
    value_11 = reference_field[(axial_index_1, perimeter_index_1)]

    interpolated = []
    for component in range(3):
        value_at_x0 = (
            (1.0 - perimeter_fraction) * value_00[component]
            + perimeter_fraction * value_01[component]
        )
        value_at_x1 = (
            (1.0 - perimeter_fraction) * value_10[component]
            + perimeter_fraction * value_11[component]
        )
        interpolated.append(
            (1.0 - axial_fraction) * value_at_x0
            + axial_fraction * value_at_x1
        )
    return tuple(interpolated)


def calculate_error_metrics(surface, reference_field):
    error_squared = 0.0
    reference_squared = 0.0

    for node in surface.Nodes:
        mapped = vector_to_tuple(node.GetSolutionStepValue(KM.DISPLACEMENT))
        reference = interpolate_reference_displacement(
            reference_field, node.X0, node.Y0, node.Z0
        )
        for component in range(3):
            error_squared += (mapped[component] - reference[component]) ** 2
            reference_squared += reference[component] ** 2

    # Wang, Eq. (5.30): component-wise displacement RMSE over 3*n_s DOFs.
    number_of_components = 3 * surface.NumberOfNodes()
    rmse = math.sqrt(error_squared / number_of_components)
    reference_rmse = math.sqrt(reference_squared / number_of_components)
    normalized_rmse_percent = (
        0.0 if reference_rmse <= 0.0 else 100.0 * rmse / reference_rmse
    )
    return rmse, reference_rmse, normalized_rmse_percent


def observed_convergence_order(coarse_error, fine_error):
    if coarse_error <= 0.0 or fine_error <= 0.0:
        return None
    # h_fine = h_coarse/2 for every adjacent pair in FLUID_MESH_LEVELS.
    return math.log(coarse_error / fine_error) / math.log(2.0)


def mesh_label(cross_section_divisions, axial_divisions):
    return f"{cross_section_divisions}x{axial_divisions}"


def format_pi_multiple(value):
    return f"{value / math.pi:g}pi"


def case_output_root(test_case):
    return OUTPUT_ROOT / test_case["name"]


def kinematics_log_rows(test_case):
    return [
        ("test_case", test_case["name"]),
        ("tip_axial_displacement", test_case["tip_axial_displacement"]),
        ("cubic_tip_displacement_y", test_case["cubic_tip_displacement_y"]),
        ("cubic_tip_displacement_z", test_case["cubic_tip_displacement_z"]),
        (
            "tip_bending_rotation_y",
            format_pi_multiple(test_case["tip_bending_rotation_y"]),
        ),
        (
            "tip_bending_rotation_z",
            format_pi_multiple(test_case["tip_bending_rotation_z"]),
        ),
        ("tip_torsion", format_pi_multiple(test_case["tip_torsion"])),
    ]


def validate_refinement_levels():
    previous_cross_section_divisions = None
    previous_axial_divisions = None

    for level in FLUID_MESH_LEVELS:
        cross_section_divisions = level["cross_section_divisions"]
        axial_divisions = level["axial_divisions"]

        if previous_cross_section_divisions is not None:
            if cross_section_divisions != 2 * previous_cross_section_divisions:
                raise RuntimeError("Cross-section divisions must double at each level")
            if axial_divisions != 2 * previous_axial_divisions:
                raise RuntimeError("Axial divisions must double at each level")

        if REFERENCE_CROSS_SECTION_DIVISIONS % cross_section_divisions != 0:
            raise RuntimeError("Every study cross-section mesh must divide the reference mesh")
        if REFERENCE_AXIAL_DIVISIONS % axial_divisions != 0:
            raise RuntimeError("Every study axial mesh must divide the reference mesh")

        previous_cross_section_divisions = cross_section_divisions
        previous_axial_divisions = axial_divisions

    finest_level = FLUID_MESH_LEVELS[-1]
    if (
        finest_level["cross_section_divisions"]
        != REFERENCE_CROSS_SECTION_DIVISIONS
        or finest_level["axial_divisions"] != REFERENCE_AXIAL_DIVISIONS
    ):
        raise RuntimeError("The reference mesh must equal the finest study mesh")

    if REFERENCE_BEAM_DIVISIONS % REFERENCE_AXIAL_DIVISIONS != 0:
        raise RuntimeError("The reference beam mesh must contain every fluid station")

    case_names = [test_case["name"] for test_case in KINEMATIC_TEST_CASES]
    if len(case_names) != len(set(case_names)):
        raise RuntimeError("Kinematic test-case names must be unique")

    required_case_parameters = {
        "tip_axial_displacement",
        "cubic_tip_displacement_y",
        "cubic_tip_displacement_z",
        "tip_bending_rotation_y",
        "tip_bending_rotation_z",
        "tip_torsion",
    }
    for test_case in KINEMATIC_TEST_CASES:
        missing_parameters = required_case_parameters.difference(test_case)
        if missing_parameters:
            raise RuntimeError(
                f"Test case {test_case['name']} is missing: "
                + ", ".join(sorted(missing_parameters))
            )


def print_section(title):
    print("")
    print(f"[{title}]")


def print_key_values(rows):
    width = max(len(name) for name, _ in rows)
    for name, value in rows:
        print(f"  {name:<{width}} : {value}")


def create_reference_solution(test_case):
    reference_label = mesh_label(
        REFERENCE_CROSS_SECTION_DIVISIONS, REFERENCE_AXIAL_DIVISIONS
    )
    torsion_label = format_pi_multiple(test_case["tip_torsion"])
    output_path = (
        case_output_root(test_case)
        / f"reference_{reference_label}_beam{REFERENCE_BEAM_DIVISIONS}"
        f"_torsion={torsion_label}"
    )
    vtk_path = reference_vtk_path(output_path)
    cached_reference_field = load_reference_field_from_vtk(vtk_path)
    if cached_reference_field is not None:
        print("reference VTK:", vtk_path)
        return cached_reference_field

    with redirect_stdout_to_file(output_path / "console_log.txt"):
        print("=" * 80)
        print("CO-ROTATIONAL FINEST-MESH REFERENCE")
        print("=" * 80)
        print_section("CASE CONFIGURATION")
        print_key_values(
            [
                ("reference_mapper", "beam_mapper(use_corotation=true)"),
                ("beam_divisions", REFERENCE_BEAM_DIVISIONS),
                ("beam_generation", "programmatic uniform line"),
                ("fluid_mesh", reference_label),
                ("surface_generation", "programmatic structured square tube"),
                *kinematics_log_rows(test_case),
                ("reference_vtk", vtk_path),
                ("vtk_output_path", output_path),
            ]
        )

        beam = create_reference_structure_solution(
            "reference_beam_structure", test_case
        )
        surface = create_reference_surface_model_part("reference_corotation_surface")
        run_mapping(beam, surface, "beam_mapper_corotation")

        print_section("MESH")
        print_key_values(
            [
                ("beam", f"{beam.NumberOfNodes()} nodes, {beam.NumberOfElements()} elements"),
                (
                    "surface",
                    f"{surface.NumberOfNodes()} nodes, {surface.NumberOfElements()} elements",
                ),
            ]
        )

        write_vtk(beam, "reference_beam_structure", output_path, ["DISPLACEMENT", "ROTATION"])
        write_vtk(
            surface,
            "reference_corotation_mapped_surface",
            output_path,
            ["DISPLACEMENT"],
        )
        reference_field = capture_reference_field(surface)
        if not vtk_path.is_file():
            raise RuntimeError(f"Reference VTK was not written: {vtk_path}")
        return reference_field


def run_refinement_case(test_case, level, reference_field, previous_error):
    cross_section_divisions = level["cross_section_divisions"]
    axial_divisions = level["axial_divisions"]
    label = mesh_label(cross_section_divisions, axial_divisions)
    reference_label = mesh_label(
        REFERENCE_CROSS_SECTION_DIVISIONS, REFERENCE_AXIAL_DIVISIONS
    )
    output_path = case_output_root(test_case) / label

    with redirect_stdout_to_file(output_path / "console_log.txt"):
        print("=" * 80)
        print("BEAM SPLINE FLUID-MESH REFINEMENT")
        print("=" * 80)
        print_section("CASE CONFIGURATION")
        print_key_values(
            [
                ("mapper", "beam_spline_mapper"),
                (
                    "reference_mapper",
                    f"beam_mapper(use_corotation=true), {reference_label}",
                ),
                ("beam_mdpa", STUDY_BEAM_MDPA_PATH),
                ("fluid_mesh", label),
                *kinematics_log_rows(test_case),
                ("vtk_output_path", output_path),
            ]
        )

        beam = create_structure_solution(
            f"beam_spline_{label}_beam", STUDY_BEAM_MDPA_PATH, test_case
        )
        surface = create_structured_surface_model_part(
            f"beam_spline_{label}_surface",
            cross_section_divisions,
            axial_divisions,
        )
        run_mapping(beam, surface, "beam_spline_mapper")
        rmse, reference_rmse, normalized_rmse_percent = calculate_error_metrics(
            surface, reference_field
        )
        convergence_order = (
            None
            if previous_error is None
            else observed_convergence_order(previous_error, rmse)
        )

        axial_element_size = BEAM_LENGTH / axial_divisions
        cross_section_element_size = CROSS_SECTION_SIZE / cross_section_divisions
        characteristic_mesh_size = max(
            axial_element_size, cross_section_element_size
        )

        write_vtk(beam, f"beam_spline_{label}_beam", output_path, ["DISPLACEMENT", "ROTATION"])
        write_vtk(
            surface,
            f"beam_spline_{label}_mapped_surface",
            output_path,
            ["DISPLACEMENT"],
        )

        print_section("MESH")
        print_key_values(
            [
                ("beam", f"{beam.NumberOfNodes()} nodes, {beam.NumberOfElements()} elements"),
                (
                    "surface",
                    f"{surface.NumberOfNodes()} nodes, {surface.NumberOfElements()} elements",
                ),
                ("axial_element_size", axial_element_size),
                ("cross_section_element_size", cross_section_element_size),
                ("characteristic_mesh_size", characteristic_mesh_size),
            ]
        )
        print_section("ERROR SUMMARY")
        print_key_values(
            [
                ("rmse_eq_5_30", f"{rmse:.8e}"),
                ("reference_rmse", f"{reference_rmse:.8e}"),
                ("normalized_rmse", f"{normalized_rmse_percent:.8f}%"),
                (
                    "convergence_order",
                    "not_available"
                    if convergence_order is None
                    else f"{convergence_order:.8f}",
                ),
            ]
        )

    return {
        "test_case": test_case["name"],
        "mesh_label": label,
        "cross_section_divisions": cross_section_divisions,
        "axial_divisions": axial_divisions,
        "axial_element_size": axial_element_size,
        "cross_section_element_size": cross_section_element_size,
        "characteristic_mesh_size": characteristic_mesh_size,
        "surface_nodes": surface.NumberOfNodes(),
        "surface_elements": surface.NumberOfElements(),
        "rmse": rmse,
        "reference_rmse": reference_rmse,
        "normalized_rmse_percent": normalized_rmse_percent,
        "convergence_order": convergence_order,
    }


def write_summary(test_case, results):
    summary_path = case_output_root(test_case) / "fluid_mesh_refinement_summary.txt"
    reference_label = mesh_label(
        REFERENCE_CROSS_SECTION_DIVISIONS, REFERENCE_AXIAL_DIVISIONS
    )
    with open(summary_path, "a") as summary_file:
        summary_file.write("\n")
        summary_file.write(
            f"# run_started: {datetime.datetime.now().isoformat(timespec='seconds')}\n"
        )
        summary_file.write("# study_mapper: beam_spline_mapper\n")
        for name, value in kinematics_log_rows(test_case):
            summary_file.write(f"# {name}: {value}\n")
        summary_file.write(f"# study_beam_mdpa: {STUDY_BEAM_MDPA_PATH}\n")
        summary_file.write(
            "# reference_mapper: beam_mapper(use_corotation=true), "
            f"programmatic fluid mesh {reference_label}\n"
        )
        summary_file.write(
            f"# reference_beam: programmatic uniform mesh, "
            f"{REFERENCE_BEAM_DIVISIONS} divisions\n"
        )
        summary_file.write("# rmse_definition: Wang Eq. (5.30)\n")
        summary_file.write(
            "# convergence_order: log(rmse_h/rmse_h_over_2)/log(2)\n"
        )
        summary_file.write(
            "mesh_label cross_section_divisions axial_divisions "
            "cross_section_element_size axial_element_size characteristic_mesh_size "
            "surface_nodes surface_elements rmse reference_rmse "
            "normalized_rmse_percent convergence_order\n"
        )
        for result in results:
            convergence_order = result["convergence_order"] 
            summary_file.write(
                f"{result['mesh_label']} "
                f"{result['cross_section_divisions']} "
                f"{result['axial_divisions']} "
                f"{result['cross_section_element_size']:.16e} "
                f"{result['axial_element_size']:.16e} "
                f"{result['characteristic_mesh_size']:.16e} "
                f"{result['surface_nodes']} "
                f"{result['surface_elements']} "
                f"{result['rmse']:.16e} "
                f"{result['reference_rmse']:.16e} "
                f"{result['normalized_rmse_percent']:.8f} "
                f"{'nan' if convergence_order is None else f'{convergence_order:.8f}'}\n"
            )
    return summary_path


def read_latest_summary_results(test_case):
    """Read the latest complete result block without rerunning the mapping."""
    summary_path = case_output_root(test_case) / "fluid_mesh_refinement_summary.txt"
    if not summary_path.is_file():
        raise FileNotFoundError(f"No existing summary found at {summary_path}")

    header_prefix = "mesh_label cross_section_divisions axial_divisions "
    lines = summary_path.read_text().splitlines()
    header_indices = [
        index for index, line in enumerate(lines) if line.startswith(header_prefix)
    ]
    if not header_indices:
        raise RuntimeError(f"No result block found in {summary_path}")

    results = []
    for line in lines[header_indices[-1] + 1:]:
        if not line.strip():
            if results:
                break
            continue
        if line.startswith("#"):
            break

        values = line.split()
        convergence_order = None if values[11] == "nan" else float(values[11])
        results.append(
            {
                "mesh_label": values[0],
                "cross_section_divisions": int(values[1]),
                "axial_divisions": int(values[2]),
                "cross_section_element_size": float(values[3]),
                "axial_element_size": float(values[4]),
                "characteristic_mesh_size": float(values[5]),
                "surface_nodes": int(values[6]),
                "surface_elements": int(values[7]),
                "rmse": float(values[8]),
                "reference_rmse": float(values[9]),
                "normalized_rmse_percent": float(values[10]),
                "convergence_order": convergence_order,
            }
        )

    if len(results) != len(FLUID_MESH_LEVELS):
        raise RuntimeError(
            f"Latest summary block in {summary_path} contains {len(results)} "
            f"levels; expected {len(FLUID_MESH_LEVELS)}"
        )
    return results


def configure_error_mesh_axis(axis, results):
    mesh_levels = list(range(len(results)))
    axis.set_xticks(mesh_levels)
    axis.set_xticklabels(
        [
            rf"${result['cross_section_divisions']}\times"
            rf"{result['axial_divisions']}$"
            for result in results
        ]
    )
    axis.set_xlim(-0.25, len(results) - 0.75)
    axis.set_xlabel("Surface mesh")
    axis.grid(True, axis="y", color=GRID_COLOR, linewidth=0.6)


def set_error_limits(axis, values):
    value_min = min(values)
    value_max = max(values)
    value_span = max(value_max - value_min, 0.05 * max(abs(value_max), 1.0e-12))
    axis.set_ylim(
        max(0.0, value_min - 0.20 * value_span),
        value_max + 0.30 * value_span,
    )


def save_figure(figure, output_path, thesis_filename=None):
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    figure.savefig(output_path.with_suffix(".svg"), bbox_inches="tight")
    if thesis_filename is not None:
        THESIS_FIGURE_DIR.mkdir(parents=True, exist_ok=True)
        thesis_output_path = THESIS_FIGURE_DIR / thesis_filename
        figure.savefig(thesis_output_path, dpi=300, bbox_inches="tight")
        figure.savefig(thesis_output_path.with_suffix(".svg"), bbox_inches="tight")


def create_plots(test_case, results):
    plot_path = case_output_root(test_case) / "Summary_Plots"
    plot_path.mkdir(parents=True, exist_ok=True)
    rmse_results = sorted(results, key=lambda result: result["rmse"])
    normalized_rmse_results = sorted(
        results, key=lambda result: result["normalized_rmse_percent"]
    )
    mesh_levels = list(range(len(results)))
    rmse_values = [result["rmse"] for result in rmse_results]
    normalized_rmse_values = [
        result["normalized_rmse_percent"] for result in normalized_rmse_results
    ]

    figure, axis = plt.subplots(figsize=(6.6, 3.8))
    axis.plot(
        mesh_levels,
        rmse_values,
        "o-",
        color=BEAM_SPLINE_COLOR,
        markerfacecolor=DATA_MARKER_COLOR,
        markeredgecolor=BEAM_SPLINE_COLOR,
        linewidth=1.6,
        markersize=5.2,
    )
    configure_error_mesh_axis(axis, rmse_results)
    set_error_limits(axis, rmse_values)
    axis.set_ylabel(r"$e_{\mathrm{RMSE}}$")
    figure.tight_layout()
    save_figure(
        figure,
        plot_path / "rmse_vs_fluid_mesh_size.png",
        (
            "biaxial_bending_twist_rmse_vs_fluid_mesh_size.png"
            if test_case["name"] == "biaxial_bending_twist"
            else None
        ),
    )
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(6.6, 3.8))
    axis.plot(
        mesh_levels,
        normalized_rmse_values,
        "o-",
        color=BEAM_SPLINE_COLOR,
        markerfacecolor=DATA_MARKER_COLOR,
        markeredgecolor=BEAM_SPLINE_COLOR,
        linewidth=1.6,
        markersize=5.2,
    )
    configure_error_mesh_axis(axis, normalized_rmse_results)
    set_error_limits(axis, normalized_rmse_values)
    axis.set_ylabel(r"$e_{\mathrm{norm}}$ [\%]")
    figure.tight_layout()
    save_figure(
        figure,
        plot_path / "normalized_rmse_vs_fluid_mesh_size.png",
        (
            "biaxial_bending_twist_normalized_rmse_vs_fluid_mesh_size.png"
            if test_case["name"] == "biaxial_bending_twist"
            else None
        ),
    )
    plt.close(figure)

    order_results = [
        result for result in results if result["convergence_order"] is not None
    ]
    order_results = list(reversed(order_results))
    orders = [result["convergence_order"] for result in order_results]
    refinement_steps = list(range(1, len(order_results) + 1))
    transition_labels = [
        (
            rf"$\frac{{{results[index]['cross_section_divisions']}\times"
            rf"{results[index]['axial_divisions']}}}"
            rf"{{{results[index - 1]['cross_section_divisions']}\times"
            rf"{results[index - 1]['axial_divisions']}}}$"
        )
        for index in range(len(results) - 1, 0, -1)
    ]

    figure, axis = plt.subplots(figsize=(6.6, 3.8))
    axis.plot(
        refinement_steps,
        orders,
        "o-",
        color=BEAM_SPLINE_COLOR,
        markerfacecolor=DATA_MARKER_COLOR,
        markeredgecolor=BEAM_SPLINE_COLOR,
        linewidth=1.6,
        markersize=5.2,
    )
    axis.set_xticks(refinement_steps)
    axis.set_xticklabels(transition_labels)
    axis.set_xlim(0.65, len(refinement_steps) + 0.35)
    axis.axhline(
        0.0,
        color=REFERENCE_LINE_COLOR,
        linestyle="--",
        linewidth=1.0,
    )

    order_min = min(orders)
    order_max = max(orders)
    order_span = max(order_max - order_min, 0.05)
    axis.set_ylim(
        min(0.0, order_min - 0.20 * order_span),
        order_max + 0.35 * order_span,
    )

    axis.set_xlabel("Fine/coarse surface mesh pair")
    axis.set_ylabel(r"Observed order $p$")
    axis.grid(True, axis="y", color=GRID_COLOR, linewidth=0.6)
    figure.tight_layout()
    save_figure(
        figure,
        plot_path / "convergence_order_vs_fluid_mesh_size.png",
        (
            "biaxial_bending_twist_observed_order.png"
            if test_case["name"] == "biaxial_bending_twist"
            else None
        ),
    )
    plt.close(figure)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run beam-spline fluid-mesh refinement studies."
    )
    parser.add_argument(
        "--case",
        choices=[test_case["name"] for test_case in KINEMATIC_TEST_CASES],
        help="Run only the selected kinematic test case.",
    )
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Regenerate plots from the latest saved summary without mapping.",
    )
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    validate_refinement_levels()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    selected_test_cases = [
        test_case
        for test_case in KINEMATIC_TEST_CASES
        if arguments.case is None or test_case["name"] == arguments.case
    ]

    for test_case in selected_test_cases:
        case_root = case_output_root(test_case)
        case_root.mkdir(parents=True, exist_ok=True)
        if arguments.plots_only:
            results = read_latest_summary_results(test_case)
            create_plots(test_case, results)
            print("plots regenerated:", case_root / "Summary_Plots")
            continue

        reference_field = create_reference_solution(test_case)

        results = []
        previous_error = None
        for level in FLUID_MESH_LEVELS:
            result = run_refinement_case(
                test_case, level, reference_field, previous_error
            )
            results.append(result)
            previous_error = result["rmse"]
            print(
                "Case completed",
                test_case["name"],
                result["mesh_label"],
                "rmse=",
                result["rmse"],
                "normalized_rmse_percent=",
                result["normalized_rmse_percent"],
                "convergence_order=",
                result["convergence_order"],
            )

        summary_path = write_summary(test_case, results)
        create_plots(test_case, results)
        print("summary:", summary_path)
        print("plots:", case_root / "Summary_Plots")


if __name__ == "__main__":
    main()
