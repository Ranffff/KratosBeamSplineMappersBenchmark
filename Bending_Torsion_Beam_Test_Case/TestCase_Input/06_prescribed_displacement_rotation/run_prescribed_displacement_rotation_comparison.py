from pathlib import Path
import importlib.util
import math

import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication
import KratosMultiphysics.StructuralMechanicsApplication


REPO_ROOT = Path(__file__).resolve().parents[3]
CASE_DIR = Path(__file__).resolve().parent
SOLVER_PATH = REPO_ROOT / "Bending_Torsion_Beam_Test_Case" / "BernoulliSolver.py"
MDPA_PATH = CASE_DIR / "beam_geometry_coarse.mdpa"
MATERIALS_PATH = (
    REPO_ROOT
    / "Bending_Torsion_Beam_Test_Case"
    / "Bending_Torsion_Beam_Test_Case"
    / "beam_geometry"
    / "StructuralMaterials.json"
)

FIXED_NODE_ID = 19
PRESCRIBED_NODE_ID = 1
PRESCRIBED_DOF_VALUES = {
    "ux": 0.10,
    "uy": -0.25,
    "uz": 0.15,
    "rx": 0.03,
    "ry": -0.04,
    "rz": 0.05,
}


def import_bernoulli_solver():
    spec = importlib.util.spec_from_file_location("BernoulliSolver", SOLVER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def vector_to_tuple(value):
    return (float(value[0]), float(value[1]), float(value[2]))


def create_reference_solution():
    solver = import_bernoulli_solver()

    model = KM.Model()
    reference = model.CreateModelPart("reference_bernoulli_solution")
    reference.ProcessInfo[KM.DOMAIN_SIZE] = 3
    reference.ProcessInfo[KM.TIME] = 0.0
    reference.ProcessInfo[KM.DELTA_TIME] = 1.0
    reference.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    reference.AddNodalSolutionStepVariable(KM.ROTATION)
    reference.AddNodalSolutionStepVariable(solver.POINT_LOAD)
    reference.AddNodalSolutionStepVariable(solver.POINT_MOMENT)

    KM.ModelPartIO(str(MDPA_PATH.with_suffix(""))).ReadModelPart(reference)

    solver.fix_all_node_dofs(reference, FIXED_NODE_ID)
    for dof_name, value in PRESCRIBED_DOF_VALUES.items():
        solver.prescribe_node_dof(reference, PRESCRIBED_NODE_ID, dof_name, value)

    node_equation_ids = {node.Id: index for index, node in enumerate(reference.Nodes)}
    material = solver.read_material(MATERIALS_PATH, reference)
    stiffness = solver.assemble_global_stiffness(reference, material, node_equation_ids)
    force = solver.build_force_vector(reference, node_equation_ids)
    solution = solver.solve_reduced_system(reference, stiffness, force, node_equation_ids)
    solver.write_solution_to_nodes(reference, solution, node_equation_ids)

    return reference


def create_mapping_origin(reference):
    model = KM.Model()
    origin = model.CreateModelPart("mapping_origin")
    origin.ProcessInfo[KM.DOMAIN_SIZE] = 3
    origin.ProcessInfo[KM.TIME] = 0.0
    origin.ProcessInfo[KM.DELTA_TIME] = 1.0
    origin.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    origin.AddNodalSolutionStepVariable(KM.ROTATION)

    properties = origin.CreateNewProperties(0)
    for node in reference.Nodes:
        new_node = origin.CreateNewNode(node.Id, node.X0, node.Y0, node.Z0)
        new_node.SetSolutionStepValue(KM.DISPLACEMENT, list(vector_to_tuple(node.GetSolutionStepValue(KM.DISPLACEMENT))))
        new_node.SetSolutionStepValue(KM.ROTATION, list(vector_to_tuple(node.GetSolutionStepValue(KM.ROTATION))))

    for element in reference.Elements:
        node_ids = [node.Id for node in element.GetNodes()]
        origin.CreateNewElement("CrBeamElement3D2N", element.Id, node_ids, properties)

    return origin


def create_destination_nodes(reference, name):
    model = KM.Model()
    destination = model.CreateModelPart(name)
    destination.ProcessInfo[KM.DOMAIN_SIZE] = 3
    destination.ProcessInfo[KM.TIME] = 0.0
    destination.ProcessInfo[KM.DELTA_TIME] = 1.0
    destination.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

    for node in reference.Nodes:
        destination.CreateNewNode(node.Id, node.X0, node.Y0, node.Z0)

    return destination


def run_mapping(origin, destination, mapper_type):
    if mapper_type == "beam_spline_mapper":
        settings = KM.Parameters("""{
            "mapper_type" : "beam_spline_mapper",
            "search_settings" : {
                "search_radius" : 2.0,
                "max_num_search_iterations" : 30
            },
            "local_coord_tolerance" : 0.25,
            "echo_level" : 0
        }""")
    elif mapper_type == "beam_mapper":
        settings = KM.Parameters("""{
            "mapper_type" : "beam_mapper",
            "search_iterations" : 30,
            "use_corotation" : false,
            "echo_level" : 0
        }""")
    else:
        raise RuntimeError(f"Unsupported mapper type: {mapper_type}")

    mapper = KM.MapperFactory.CreateMapper(origin, destination, settings)
    mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())


def norm_3d(lhs, rhs):
    return math.sqrt(sum((lhs[i] - rhs[i]) ** 2 for i in range(3)))


def main():
    print("PRESCRIBED DISPLACEMENT AND ROTATION COMPARISON")
    print("input mdpa:", MDPA_PATH)
    print("fixed node:", FIXED_NODE_ID, "(all displacement and rotation DOFs)")
    print("prescribed node:", PRESCRIBED_NODE_ID, PRESCRIBED_DOF_VALUES)
    print("")

    reference = create_reference_solution()
    mapping_origin = create_mapping_origin(reference)

    spline_destination = create_destination_nodes(reference, "beam_spline_destination")
    beam_mapper_destination = create_destination_nodes(reference, "beam_mapper_destination")

    run_mapping(mapping_origin, spline_destination, "beam_spline_mapper")
    run_mapping(mapping_origin, beam_mapper_destination, "beam_mapper")

    max_spline_error = 0.0
    max_beam_mapper_error = 0.0
    max_spline_error_node = None
    max_beam_mapper_error_node = None

    print("node_id, x")
    print("  BernoulliSolver displacement")
    print("  Beam_Spline_mapper displacement, error_norm")
    print("  beam_mapper(use_corotation=false) displacement, error_norm")

    for reference_node in sorted(reference.Nodes, key=lambda node: node.X0):
        node_id = reference_node.Id
        reference_displacement = vector_to_tuple(reference_node.GetSolutionStepValue(KM.DISPLACEMENT))
        spline_displacement = vector_to_tuple(
            spline_destination.GetNode(node_id).GetSolutionStepValue(KM.DISPLACEMENT)
        )
        beam_mapper_displacement = vector_to_tuple(
            beam_mapper_destination.GetNode(node_id).GetSolutionStepValue(KM.DISPLACEMENT)
        )

        spline_error = norm_3d(spline_displacement, reference_displacement)
        beam_mapper_error = norm_3d(beam_mapper_displacement, reference_displacement)

        if spline_error > max_spline_error:
            max_spline_error = spline_error
            max_spline_error_node = node_id
        if beam_mapper_error > max_beam_mapper_error:
            max_beam_mapper_error = beam_mapper_error
            max_beam_mapper_error_node = node_id

        print(node_id, reference_node.X0)
        print(" ", reference_displacement)
        print(" ", spline_displacement, spline_error)
        print(" ", beam_mapper_displacement, beam_mapper_error)

    print("")
    print("max_beam_spline_mapper_error_norm:", max_spline_error)
    print("max_beam_spline_mapper_error_node:", max_spline_error_node)
    print("max_beam_mapper_error_norm:", max_beam_mapper_error)
    print("max_beam_mapper_error_node:", max_beam_mapper_error_node)


if __name__ == "__main__":
    main()
