# Minimal Kratos + NumPy 3D Euler-Bernoulli beam solver.
#
# Run from /root/dev/Kratos.
# **************************************************************************************************************
# Command-line interface:
#   General syntax:
#     python3 Bending_Torsion_Beam_Test_Case/BernoulliSolver.py [ProjectParameters.json|model.mdpa] [options]
#
#   Positional ProjectParameters.json:
#     python3 Bending_Torsion_Beam_Test_Case/BernoulliSolver.py path/to/ProjectParameters.json
#     Reads mdpa path, materials path, JSON processes, and JSON output settings
#     from the ProjectParameters file.
#
#   Positional model.mdpa:
#     python3 Bending_Torsion_Beam_Test_Case/BernoulliSolver.py path/to/model.mdpa [options]
#     Reads geometry and mdpa Properties from the mdpa file directly. Loads and
#     boundary conditions can be supplied through command-line options such as
#     --load, --fix-node, and --prescribe.
#
#   --project ProjectParameters.json
#     Explicitly selects the Kratos ProjectParameters file.
#
#   --materials StructuralMaterials.json
#     Overrides the material json file. Material data are read with priority:
#       StructuralMaterials.json > mdpa Properties > solver defaults
#
#   --load node_id fx fy fz mx my mz
#     Adds a nodal load to node_id. Values are in the global coordinate system:
#       fx, fy, fz : point forces in global X, Y, Z
#       mx, my, mz : point moments about global X, Y, Z
#     This option can be repeated. Repeated loads are accumulated.
#
#   --fix-node node_id
#     Fully fixes all six degrees of freedom of node_id to zero:
#       ux = uy = uz = rx = ry = rz = 0
#     This option can be repeated for multiple fixed nodes.
#
#   --prescribe node_id dof value
#     Prescribes one nodal degree of freedom. Valid dof names are:
#       ux, uy, uz : global displacement components
#       rx, ry, rz : global rotation components
#     This option can be repeated to define pinned, roller, or mixed supports.
#
#   --output vtk_output_directory
#     Sets the VTK output directory. The main model file is written as:
#       vtk_output_directory/BeamModelPart_0_0.vtk
#
#   Input model directory for validation cases:
#     Bending_Torsion_Beam_Test_Case/TestCase_Input/
#
#   Material file used by the validation cases:
#     Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json
# **************************************************************************************************************
# General project run:
#   python3 Bending_Torsion_Beam_Test_Case/BernoulliSolver.py \
#     --project Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/ProjectParameters.json \
#     --output Bending_Torsion_Beam_Test_Case/BernoulliSolver_python_TestCases/project_run
# **************************************************************************************************************
# CLI validation testcase 01: cantilever beam, free-end point load Fy = -1000.
#   Input:
#     Bending_Torsion_Beam_Test_Case/TestCase_Input/01_cantilever_point_load/beam_geometry_coarse.mdpa
#   Output:
#     Bending_Torsion_Beam_Test_Case/BernoulliSolver_python_TestCases/01_cantilever_point_load
#   Command:
#     python3 Bending_Torsion_Beam_Test_Case/BernoulliSolver.py \
#       Bending_Torsion_Beam_Test_Case/TestCase_Input/01_cantilever_point_load/beam_geometry_coarse.mdpa \
#       --materials Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json \
#       --fix-node 19 \
#       --load 1 0 -1000 0 0 0 0 \
#       --output Bending_Torsion_Beam_Test_Case/BernoulliSolver_python_TestCases/01_cantilever_point_load
#
# CLI validation testcase 02: simply supported beam, center point load Fy = -1000.
#   Input:
#     Bending_Torsion_Beam_Test_Case/TestCase_Input/02_simply_supported_point_load/beam_geometry_coarse.mdpa
#   Output:
#     Bending_Torsion_Beam_Test_Case/BernoulliSolver_python_TestCases/02_simply_supported_point_load
#   Command:
#     python3 Bending_Torsion_Beam_Test_Case/BernoulliSolver.py \
#       Bending_Torsion_Beam_Test_Case/TestCase_Input/02_simply_supported_point_load/beam_geometry_coarse.mdpa \
#       --materials Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json \
#       --prescribe 19 ux 0 --prescribe 19 uy 0 --prescribe 19 uz 0 --prescribe 19 rx 0 \
#       --prescribe 1 uy 0 --prescribe 1 uz 0 --prescribe 1 rx 0 \
#       --load 10 0 -1000 0 0 0 0 \
#       --output Bending_Torsion_Beam_Test_Case/BernoulliSolver_python_TestCases/02_simply_supported_point_load
#
# CLI validation testcase 03: simply supported beam, uniform distributed load qy = -1000.
#   Input:
#     Bending_Torsion_Beam_Test_Case/TestCase_Input/03_simply_supported_udl/simply_supported_udl.mdpa
#   Output:
#     Bending_Torsion_Beam_Test_Case/BernoulliSolver_python_TestCases/03_simply_supported_udl
#   Command:
#     python3 Bending_Torsion_Beam_Test_Case/BernoulliSolver.py \
#       Bending_Torsion_Beam_Test_Case/TestCase_Input/03_simply_supported_udl/simply_supported_udl.mdpa \
#       --materials Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json \
#       --prescribe 1 ux 0 --prescribe 1 uy 0 --prescribe 1 uz 0 --prescribe 1 rx 0 \
#       --prescribe 11 uy 0 --prescribe 11 uz 0 --prescribe 11 rx 0 \
#       --load 1 0 -500 0 0 0 -83.3333333333 \
#       --load 2 0 -1000 0 0 0 0 \
#       --load 3 0 -1000 0 0 0 0 \
#       --load 4 0 -1000 0 0 0 0 \
#       --load 5 0 -1000 0 0 0 0 \
#       --load 6 0 -1000 0 0 0 0 \
#       --load 7 0 -1000 0 0 0 0 \
#       --load 8 0 -1000 0 0 0 0 \
#       --load 9 0 -1000 0 0 0 0 \
#       --load 10 0 -1000 0 0 0 0 \
#       --load 11 0 -500 0 0 0 83.3333333333 \
#       --output Bending_Torsion_Beam_Test_Case/BernoulliSolver_python_TestCases/03_simply_supported_udl
#
# CLI validation testcase 04: fixed-free torsion beam, free-end torque Mx = +1000.
#   Input:
#     Bending_Torsion_Beam_Test_Case/TestCase_Input/04_torsion/beam_geometry_coarse.mdpa
#   Output:
#     Bending_Torsion_Beam_Test_Case/BernoulliSolver_python_TestCases/04_torsion
#   Command:
#     python3 Bending_Torsion_Beam_Test_Case/BernoulliSolver.py \
#       Bending_Torsion_Beam_Test_Case/TestCase_Input/04_torsion/beam_geometry_coarse.mdpa \
#       --materials Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json \
#       --fix-node 19 \
#       --load 1 0 0 0 1000 0 0 \
#       --output Bending_Torsion_Beam_Test_Case/BernoulliSolver_python_TestCases/04_torsion
#
# CLI validation testcase 05: fixed-free combined loading, Mx = +1000 and Mz = +500.
#   Input:
#     Bending_Torsion_Beam_Test_Case/TestCase_Input/05_combined_bending_torsion/beam_geometry_coarse.mdpa
#   Output:
#     Bending_Torsion_Beam_Test_Case/BernoulliSolver_python_TestCases/05_combined_bending_torsion
#   Command:
#     python3 Bending_Torsion_Beam_Test_Case/BernoulliSolver.py \
#       Bending_Torsion_Beam_Test_Case/TestCase_Input/05_combined_bending_torsion/beam_geometry_coarse.mdpa \
#       --materials Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json \
#       --fix-node 19 \
#       --load 1 0 0 0 1000 0 500 \
#       --output Bending_Torsion_Beam_Test_Case/BernoulliSolver_python_TestCases/05_combined_bending_torsion
# **************************************************************************************************************

import argparse
import json
from pathlib import Path

import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np


DOFS_PER_NODE = 6

POINT_LOAD = getattr(SMA, "POINT_LOAD")
POINT_LOAD_X = getattr(SMA, "POINT_LOAD_X")
POINT_LOAD_Y = getattr(SMA, "POINT_LOAD_Y")
POINT_LOAD_Z = getattr(SMA, "POINT_LOAD_Z")
POINT_MOMENT = getattr(SMA, "POINT_MOMENT")
POINT_MOMENT_X = getattr(SMA, "POINT_MOMENT_X")
POINT_MOMENT_Y = getattr(SMA, "POINT_MOMENT_Y")
POINT_MOMENT_Z = getattr(SMA, "POINT_MOMENT_Z")
CROSS_AREA = getattr(SMA, "CROSS_AREA")
I22 = getattr(SMA, "I22")
I33 = getattr(SMA, "I33")
TORSIONAL_INERTIA = getattr(SMA, "TORSIONAL_INERTIA")


DISPLACEMENT_COMPONENTS = [
    KM.DISPLACEMENT_X,
    KM.DISPLACEMENT_Y,
    KM.DISPLACEMENT_Z,
]
ROTATION_COMPONENTS = [
    KM.ROTATION_X,
    KM.ROTATION_Y,
    KM.ROTATION_Z,
]
POINT_LOAD_COMPONENTS = [
    POINT_LOAD_X,
    POINT_LOAD_Y,
    POINT_LOAD_Z,
]
POINT_MOMENT_COMPONENTS = [
    POINT_MOMENT_X,
    POINT_MOMENT_Y,
    POINT_MOMENT_Z,
]


def local_stiffness_3d_beam(E, G, A, Iy, Iz, J, L):
    """
    Local DOF order per node:
    [ux, uy, uz, rx, ry, rz]

    Element DOF order:
    [ux1, uy1, uz1, rx1, ry1, rz1, ux2, uy2, uz2, rx2, ry2, rz2]
    """
    k = np.zeros((12, 12))

    EA_L = E * A / L
    k[0, 0] = EA_L
    k[0, 6] = -EA_L
    k[6, 0] = -EA_L
    k[6, 6] = EA_L

    GJ_L = G * J / L
    k[3, 3] = GJ_L
    k[3, 9] = -GJ_L
    k[9, 3] = -GJ_L
    k[9, 9] = GJ_L

    c1 = 12.0 * E * Iz / L**3
    c2 = 6.0 * E * Iz / L**2
    c3 = 4.0 * E * Iz / L
    c4 = 2.0 * E * Iz / L
    dofs = [1, 5, 7, 11]
    kb = np.array([
        [c1, c2, -c1, c2],
        [c2, c3, -c2, c4],
        [-c1, -c2, c1, -c2],
        [c2, c4, -c2, c3],
    ])
    for i in range(4):
        for j in range(4):
            k[dofs[i], dofs[j]] += kb[i, j]

    c1 = 12.0 * E * Iy / L**3
    c2 = 6.0 * E * Iy / L**2
    c3 = 4.0 * E * Iy / L
    c4 = 2.0 * E * Iy / L
    dofs = [2, 4, 8, 10]
    kb = np.array([
        [c1, -c2, -c1, -c2],
        [-c2, c3, c2, c4],
        [-c1, c2, c1, c2],
        [-c2, c4, c2, c3],
    ])
    for i in range(4):
        for j in range(4):
            k[dofs[i], dofs[j]] += kb[i, j]

    return k


def strip_json_comments(text):
    result = []
    in_string = False
    escaped = False
    i = 0
    while i < len(text):
        char = text[i]
        if in_string:
            result.append(char)
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                in_string = False
            i += 1
            continue

        if char == '"':
            in_string = True
            result.append(char)
            i += 1
            continue

        if char == "/" and i + 1 < len(text) and text[i + 1] == "/":
            while i < len(text) and text[i] != "\n":
                i += 1
            continue

        result.append(char)
        i += 1
    return "".join(result)


def read_json_file(path):
    with open(path, "r", encoding="utf-8") as json_file:
        return json.loads(strip_json_comments(json_file.read()))


def with_mdpa_extension(path):
    path = Path(path)
    if path.suffix != ".mdpa":
        path = path.with_suffix(".mdpa")
    return path


def without_mdpa_extension(path):
    path = Path(path)
    if path.suffix == ".mdpa":
        return path.with_suffix("")
    return path


def path_exists_with_optional_mdpa_extension(path):
    path = Path(path)
    return path.exists() or with_mdpa_extension(path).exists()


def resolve_project_input_path(project_directory, relative_path):
    path = Path(relative_path)
    if path.is_absolute():
        return path

    path_from_project_directory = project_directory / path
    if path_exists_with_optional_mdpa_extension(path_from_project_directory):
        return path_from_project_directory

    path_from_case_directory = project_directory.parent / path
    if path_exists_with_optional_mdpa_extension(path_from_case_directory):
        return path_from_case_directory

    return path_from_project_directory


def read_project_settings(project_parameters_file):
    project_parameters_file = Path(project_parameters_file)
    project_parameters = read_json_file(project_parameters_file)
    project_directory = project_parameters_file.parent

    solver_settings = project_parameters["solver_settings"]
    if (
        "model_import_settings" not in solver_settings
        and "solvers" in solver_settings
        and "beam_structure" in solver_settings["solvers"]
    ):
        nested_file = solver_settings["solvers"]["beam_structure"]["solver_wrapper_settings"]["input_file"]
        nested_file = Path(nested_file)
        if not nested_file.suffix:
            nested_file = nested_file.with_suffix(".json")
        return read_project_settings(project_directory / nested_file)

    mdpa_file = resolve_project_input_path(
        project_directory,
        solver_settings["model_import_settings"]["input_filename"],
    )
    materials_file = resolve_project_input_path(
        project_directory,
        solver_settings["material_import_settings"]["materials_filename"],
    )

    output_path = Path("Bending_Torsion_Beam_Test_Case/BernoulliBeamSolverOutput_vtk")
    vtk_processes = project_parameters.get("output_processes", {}).get("vtk_output", [])
    if vtk_processes:
        vtk_parameters = vtk_processes[0]["Parameters"]
        if "output_path" in vtk_parameters:
            output_path = project_directory.parent / vtk_parameters["output_path"]

    return project_parameters, mdpa_file, materials_file, output_path


def read_material(material_file, model_part):
    material = {
        "E": 206.9e9,
        "nu": 0.29,
        "A": 1.0,
        "Iy": 1.0,
        "Iz": 1.0,
        "J": 1.0,
    }

    def apply_properties(properties):
        variable_map = {
            "E": KM.YOUNG_MODULUS,
            "nu": KM.POISSON_RATIO,
            "A": CROSS_AREA,
            "Iy": I22,
            "Iz": I33,
            "J": TORSIONAL_INERTIA,
        }
        for key, variable in variable_map.items():
            if properties.Has(variable):
                material[key] = properties[variable]

    if model_part.NumberOfElements() > 0:
        apply_properties(next(iter(model_part.Elements)).Properties)
    elif model_part.NumberOfProperties() > 0:
        apply_properties(next(iter(model_part.Properties)))

    if material_file and Path(material_file).exists():
        material_settings = read_json_file(material_file)
        variables = material_settings["properties"][0]["Material"]["Variables"]
        material["E"] = variables.get("YOUNG_MODULUS", material["E"])
        material["nu"] = variables.get("POISSON_RATIO", material["nu"])
        material["A"] = variables.get("CROSS_AREA", material["A"])
        material["Iy"] = variables.get("I22", material["Iy"])
        material["Iz"] = variables.get("I33", material["Iz"])
        material["J"] = variables.get("TORSIONAL_INERTIA", material["J"])

    material["G"] = material["E"] / (2.0 * (1.0 + material["nu"]))
    return material


def node_coordinates(node):
    return np.array([node.X, node.Y, node.Z], dtype=float)


def rotation_matrix_local_to_global(node1, node2):
    ex = node_coordinates(node2) - node_coordinates(node1)
    length = np.linalg.norm(ex)
    if length <= np.finfo(float).eps:
        raise RuntimeError("Zero length beam element.")
    ex /= length

    reference = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(ex, reference)) > 0.95:
        reference = np.array([0.0, 1.0, 0.0])

    ey = np.cross(reference, ex)
    ey /= np.linalg.norm(ey)
    ez = np.cross(ex, ey)

    return np.vstack((ex, ey, ez))


def transformation_matrix_12x12(rotation_matrix):
    T = np.zeros((12, 12))
    for block in range(4):
        start = block * 3
        T[start:start + 3, start:start + 3] = rotation_matrix
    return T


def resolve_model_part(root_model_part, model_part_name):
    names = [name for name in model_part_name.split(".") if name]
    if names and names[0] in (root_model_part.Name, "Structure"):
        names = names[1:]

    model_part = root_model_part
    for name in names:
        if not model_part.HasSubModelPart(name):
            raise RuntimeError(f"SubModelPart '{name}' not found while resolving '{model_part_name}'.")
        model_part = model_part.GetSubModelPart(name)
    return model_part


def evaluate_scalar_setting_at_node(value, node, time):
    if isinstance(value, (int, float)):
        return float(value)
    function = KM.GenericFunctionUtility(str(value))
    return function.RotateAndCallFunction(node.X, node.Y, node.Z, time, node.X0, node.Y0, node.Z0)


def set_vector_component(node, variable_name, component, value, add=False):
    if variable_name == "DISPLACEMENT":
        variable = DISPLACEMENT_COMPONENTS[component]
    elif variable_name == "ROTATION":
        variable = ROTATION_COMPONENTS[component]
    elif variable_name == "POINT_LOAD":
        variable = POINT_LOAD_COMPONENTS[component]
    elif variable_name == "POINT_MOMENT":
        variable = POINT_MOMENT_COMPONENTS[component]
    else:
        raise RuntimeError(f"Unsupported vector variable: {variable_name}")

    if add:
        value += node.GetSolutionStepValue(variable)
    node.SetSolutionStepValue(variable, value)


def fix_vector_component(node, variable_name, component):
    if variable_name == "DISPLACEMENT":
        node.Fix(DISPLACEMENT_COMPONENTS[component])
    elif variable_name == "ROTATION":
        node.Fix(ROTATION_COMPONENTS[component])
    else:
        raise RuntimeError(f"Only DISPLACEMENT and ROTATION can be constrained. Got: {variable_name}")


def parse_dof_name(dof_name):
    mapping = {
        "ux": ("DISPLACEMENT", 0),
        "uy": ("DISPLACEMENT", 1),
        "uz": ("DISPLACEMENT", 2),
        "rx": ("ROTATION", 0),
        "ry": ("ROTATION", 1),
        "rz": ("ROTATION", 2),
    }
    if dof_name not in mapping:
        raise RuntimeError("Unsupported dof name. Use ux, uy, uz, rx, ry, or rz.")
    return mapping[dof_name]


def apply_assign_vector_variable_process(root_model_part, process_parameters, time, fix_constrained_components):
    target_model_part = resolve_model_part(root_model_part, process_parameters["model_part_name"])
    variable_name = process_parameters["variable_name"]
    values = process_parameters["value"]
    constrained = [False, False, False]
    if fix_constrained_components and "constrained" in process_parameters:
        constrained = process_parameters["constrained"]

    for node in target_model_part.Nodes:
        for component in range(3):
            value = evaluate_scalar_setting_at_node(values[component], node, time)
            set_vector_component(node, variable_name, component, value, add=variable_name.startswith("POINT_"))
            if constrained[component]:
                fix_vector_component(node, variable_name, component)


def apply_vector_by_direction_to_condition_process(root_model_part, process_parameters):
    target_model_part = resolve_model_part(root_model_part, process_parameters["model_part_name"])
    variable_name = process_parameters["variable_name"]
    modulus = float(process_parameters["modulus"])
    value = np.array(process_parameters["direction"], dtype=float) * modulus

    for condition in target_model_part.Conditions:
        geometry = condition.GetGeometry()
        share = 1.0 / len(geometry)
        for node in geometry:
            for component in range(3):
                set_vector_component(node, variable_name, component, share * value[component], add=True)


def apply_project_processes(model_part, project_parameters):
    if not project_parameters or "processes" not in project_parameters:
        return

    process_time = project_parameters.get("problem_data", {}).get("end_time", 0.0)
    processes = project_parameters["processes"]

    for process in processes.get("constraints_process_list", []):
        process_name = process["process_name"]
        if process_name == "AssignVectorVariableProcess":
            apply_assign_vector_variable_process(model_part, process["Parameters"], process_time, True)
        else:
            print(f"Ignoring unsupported constraint process: {process_name}")

    for process in processes.get("loads_process_list", []):
        process_name = process["process_name"]
        if process_name == "AssignVectorByDirectionToConditionProcess":
            apply_vector_by_direction_to_condition_process(model_part, process["Parameters"])
        elif process_name == "AssignVectorVariableProcess":
            apply_assign_vector_variable_process(model_part, process["Parameters"], process_time, False)
        else:
            print(f"Ignoring unsupported load process: {process_name}")


def apply_domain_size_constraints(model_part, project_parameters):
    solver_settings = project_parameters.get("solver_settings", {}) if project_parameters else {}
    if solver_settings.get("domain_size") == 2:
        for node in model_part.Nodes:
            node.Fix(KM.DISPLACEMENT_Z)
            node.Fix(KM.ROTATION_Y)


def fix_all_node_dofs(model_part, node_id):
    node = model_part.GetNode(node_id)
    for component in range(3):
        set_vector_component(node, "DISPLACEMENT", component, 0.0)
        set_vector_component(node, "ROTATION", component, 0.0)
        fix_vector_component(node, "DISPLACEMENT", component)
        fix_vector_component(node, "ROTATION", component)


def prescribe_node_dof(model_part, node_id, dof_name, value):
    variable_name, component = parse_dof_name(dof_name)
    node = model_part.GetNode(node_id)
    set_vector_component(node, variable_name, component, value)
    fix_vector_component(node, variable_name, component)


def add_command_line_load(model_part, node_id, load):
    node = model_part.GetNode(node_id)
    for component in range(3):
        set_vector_component(node, "POINT_LOAD", component, load[component], add=True)
        set_vector_component(node, "POINT_MOMENT", component, load[component + 3], add=True)


def is_fixed_dof(node, local_dof):
    if local_dof < 3:
        return node.IsFixed(DISPLACEMENT_COMPONENTS[local_dof])
    return node.IsFixed(ROTATION_COMPONENTS[local_dof - 3])


def get_dof_value(node, local_dof):
    if local_dof < 3:
        return node.GetSolutionStepValue(DISPLACEMENT_COMPONENTS[local_dof])
    return node.GetSolutionStepValue(ROTATION_COMPONENTS[local_dof - 3])


def get_dof_load(node, local_dof):
    if local_dof < 3:
        return node.GetSolutionStepValue(POINT_LOAD_COMPONENTS[local_dof])
    return node.GetSolutionStepValue(POINT_MOMENT_COMPONENTS[local_dof - 3])


def set_dof_value(node, local_dof, value):
    if local_dof < 3:
        node.SetSolutionStepValue(DISPLACEMENT_COMPONENTS[local_dof], value)
    else:
        node.SetSolutionStepValue(ROTATION_COMPONENTS[local_dof - 3], value)


def full_dof(node_equation_ids, node_id, local_dof):
    return node_equation_ids[node_id] * DOFS_PER_NODE + local_dof


def assemble_global_stiffness(model_part, material, node_equation_ids):
    total_dofs = model_part.NumberOfNodes() * DOFS_PER_NODE
    K = np.zeros((total_dofs, total_dofs))

    for element in model_part.Elements:
        geometry = element.GetGeometry()
        if len(geometry) != 2:
            raise RuntimeError("Only 2-node beam elements are supported.")

        node1 = geometry[0]
        node2 = geometry[1]
        L = np.linalg.norm(node_coordinates(node2) - node_coordinates(node1))
        element_material = material.copy()

        properties = element.Properties
        if properties.Has(KM.YOUNG_MODULUS):
            element_material["E"] = properties[KM.YOUNG_MODULUS]
        if properties.Has(KM.POISSON_RATIO):
            element_material["nu"] = properties[KM.POISSON_RATIO]
        if properties.Has(CROSS_AREA):
            element_material["A"] = properties[CROSS_AREA]
        if properties.Has(I22):
            element_material["Iy"] = properties[I22]
        if properties.Has(I33):
            element_material["Iz"] = properties[I33]
        if properties.Has(TORSIONAL_INERTIA):
            element_material["J"] = properties[TORSIONAL_INERTIA]
        element_material["G"] = element_material["E"] / (2.0 * (1.0 + element_material["nu"]))

        k_local = local_stiffness_3d_beam(
            element_material["E"],
            element_material["G"],
            element_material["A"],
            element_material["Iy"],
            element_material["Iz"],
            element_material["J"],
            L,
        )
        T = transformation_matrix_12x12(rotation_matrix_local_to_global(node1, node2))
        k_global = T.T @ k_local @ T

        dofs = []
        for node in (node1, node2):
            for local_dof in range(DOFS_PER_NODE):
                dofs.append(full_dof(node_equation_ids, node.Id, local_dof))

        for i in range(12):
            for j in range(12):
                K[dofs[i], dofs[j]] += k_global[i, j]

    return K


def build_force_vector(model_part, node_equation_ids):
    total_dofs = model_part.NumberOfNodes() * DOFS_PER_NODE
    F = np.zeros(total_dofs)
    for node in model_part.Nodes:
        for local_dof in range(DOFS_PER_NODE):
            F[full_dof(node_equation_ids, node.Id, local_dof)] = get_dof_load(node, local_dof)
    return F


def solve_reduced_system(model_part, K, F, node_equation_ids):
    total_dofs = model_part.NumberOfNodes() * DOFS_PER_NODE
    fixed_dofs = []
    free_dofs = []
    U = np.zeros(total_dofs)

    # TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    print_nodal_displacement_vectors_to_terminal(model_part, "INITIAL DISPLACEMENT AND ROTATION VECTORS")

    for node in model_part.Nodes:
        for local_dof in range(DOFS_PER_NODE):
            equation_id = full_dof(node_equation_ids, node.Id, local_dof)
            U[equation_id] = get_dof_value(node, local_dof)
            if is_fixed_dof(node, local_dof):
                fixed_dofs.append(equation_id)
            else:
                free_dofs.append(equation_id)

    if free_dofs:
        K_ff = K[np.ix_(free_dofs, free_dofs)]
        F_f = F[free_dofs]
        if fixed_dofs:
            K_fs = K[np.ix_(free_dofs, fixed_dofs)]
            F_f = F_f - K_fs @ U[fixed_dofs]
        U[free_dofs] = np.linalg.solve(K_ff, F_f)

    print(f"Solved reduced system with {len(free_dofs)} free DOFs.")
    return U


def write_solution_to_nodes(model_part, U, node_equation_ids):
    for node in model_part.Nodes:
        for local_dof in range(DOFS_PER_NODE):
            set_dof_value(node, local_dof, U[full_dof(node_equation_ids, node.Id, local_dof)])

    # TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
    print_nodal_displacement_vectors_to_terminal(model_part, "FINAL DISPLACEMENT AND ROTATION VECTORS")


# TEMP DEBUG OUTPUT - REMOVE AFTER VERIFICATION
def print_nodal_displacement_vectors_to_terminal(model_part, label):
    print(f"BERNOULLI SOLVER DEBUG {label}")
    print("node_id, x, y, z, displacement_vector, rotation_vector")
    for node in sorted(model_part.Nodes, key=lambda item: item.Id):
        displacement = node.GetSolutionStepValue(KM.DISPLACEMENT)
        rotation = node.GetSolutionStepValue(KM.ROTATION)
        print(
            node.Id,
            node.X0,
            node.Y0,
            node.Z0,
            (displacement[0], displacement[1], displacement[2]),
            (rotation[0], rotation[1], rotation[2]),
        )


def write_vtk(model_part, output_path):
    vtk_output_parameters = KM.Parameters(r'''{
        "model_part_name" : "BeamModelPart",
        "output_path" : "vtk_output",
        "file_format" : "ascii",
        "output_precision" : 7,
        "output_sub_model_parts" : true,
        "write_deformed_configuration" : true,
        "nodal_solution_step_data_variables" : ["DISPLACEMENT", "ROTATION"],
        "nodal_data_value_variables" : [],
        "element_data_value_variables" : [],
        "condition_data_value_variables" : []
    }''')
    vtk_output_parameters["output_path"].SetString(str(output_path))
    KM.VtkOutput(model_part, vtk_output_parameters).PrintOutput()


def parse_arguments():
    parser = argparse.ArgumentParser(description="NumPy/Kratos 3D Euler-Bernoulli beam solver.")
    parser.add_argument(
        "input",
        nargs="?",
        default="Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/ProjectParameters.json",
        help="ProjectParameters.json or model.mdpa. Defaults to the beam geometry project.",
    )
    parser.add_argument("--project", help="Explicit ProjectParameters.json file.")
    parser.add_argument("--materials", help="StructuralMaterials.json override.")
    parser.add_argument("--output", help="VTK output directory override.")
    parser.add_argument(
        "--load",
        nargs=7,
        action="append",
        metavar=("NODE", "FX", "FY", "FZ", "MX", "MY", "MZ"),
        help="Add a nodal load/moment in global coordinates. Can be repeated.",
    )
    parser.add_argument("--fix-node", type=int, action="append", default=[], help="Fix all 6 DOFs of a node.")
    parser.add_argument(
        "--prescribe",
        nargs=3,
        action="append",
        metavar=("NODE", "DOF", "VALUE"),
        help="Prescribe one DOF: ux, uy, uz, rx, ry, rz. Can be repeated.",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    input_path = Path(args.project or args.input)
    input_is_project = input_path.suffix == ".json" or args.project is not None

    project_parameters = {}
    materials_file = None
    output_path = Path("Bending_Torsion_Beam_Test_Case/BernoulliBeamSolverOutput_vtk")

    if input_is_project:
        project_parameters, mdpa_file, materials_file, output_path = read_project_settings(input_path)
    else:
        mdpa_file = Path(args.input)
        materials_file = Path("Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json")

    if args.materials:
        materials_file = Path(args.materials)
    if args.output:
        output_path = Path(args.output)

    model = KM.Model()
    mp = model.CreateModelPart("BeamModelPart")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.ROTATION)
    mp.AddNodalSolutionStepVariable(POINT_LOAD)
    mp.AddNodalSolutionStepVariable(POINT_MOMENT)

    KM.ModelPartIO(str(without_mdpa_extension(mdpa_file))).ReadModelPart(mp)

    apply_domain_size_constraints(mp, project_parameters)
    apply_project_processes(mp, project_parameters)

    for node_id in args.fix_node:
        fix_all_node_dofs(mp, node_id)

    for prescribed in args.prescribe or []:
        node_id, dof_name, value = prescribed
        prescribe_node_dof(mp, int(node_id), dof_name, float(value))

    for load in args.load or []:
        node_id = int(load[0])
        add_command_line_load(mp, node_id, [float(value) for value in load[1:]])

    node_equation_ids = {node.Id: index for index, node in enumerate(mp.Nodes)}
    material = read_material(materials_file, mp)

    K = assemble_global_stiffness(mp, material, node_equation_ids)
    F = build_force_vector(mp, node_equation_ids)
    U = solve_reduced_system(mp, K, F, node_equation_ids)
    write_solution_to_nodes(mp, U, node_equation_ids)
    write_vtk(mp, output_path)

    print(f"Read {mp.NumberOfNodes()} nodes and {mp.NumberOfElements()} elements.")
    print(f"Wrote VTK output to: {output_path}")


if __name__ == "__main__":
    main()
