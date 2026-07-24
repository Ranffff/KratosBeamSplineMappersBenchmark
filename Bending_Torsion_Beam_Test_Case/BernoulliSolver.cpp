// Minimal Kratos-based 3D Euler-Bernoulli beam solver.
/*
Running examples in /root/dev/Kratos after compiling BernoulliSolver.
**************************************************************************************************************
Command-line interface:
  General syntax:
    ./Bending_Torsion_Beam_Test_Case/BernoulliSolver [ProjectParameters.json|model.mdpa] [options]

  Positional ProjectParameters.json:
    ./Bending_Torsion_Beam_Test_Case/BernoulliSolver path/to/ProjectParameters.json
    Reads mdpa path, materials path, JSON processes, and JSON output settings
    from the ProjectParameters file.

  Positional model.mdpa:
    ./Bending_Torsion_Beam_Test_Case/BernoulliSolver path/to/model.mdpa [options]
    Reads geometry and mdpa Properties from the mdpa file directly. In this
    mode, no JSON load/BC processes are applied unless --project is also given;
    loads and boundary conditions should be supplied through command-line
    options such as --load, --fix-node, and --prescribe.

  --project ProjectParameters.json
    Explicitly selects the Kratos ProjectParameters file. Equivalent to passing
    the json file positionally. If both --project and a positional mdpa are
    supplied, the mdpa path overrides the mdpa path from the project, while JSON
    processes from the project are still applied.

  --materials StructuralMaterials.json
    Overrides the material json file. Material data are read with priority:
      StructuralMaterials.json > mdpa Properties > solver defaults

  --load node_id fx fy fz mx my mz
    Adds a nodal load to node_id. Values are in the global coordinate system:
      fx, fy, fz : point forces in global X, Y, Z
      mx, my, mz : point moments about global X, Y, Z
    This option can be repeated. Repeated loads are accumulated.
    Example:
      --load 1 0 -1000 0 0 0 0
    applies a point force Fy = -1000 at node 1.

  --fix-node node_id
    Fully fixes all six degrees of freedom of node_id to zero:
      ux = uy = uz = rx = ry = rz = 0
    This option can be repeated for multiple fixed nodes.

  --prescribe node_id dof value
    Prescribes one nodal degree of freedom. Valid dof names are:
      ux, uy, uz : global displacement components
      rx, ry, rz : global rotation components
    This option can be repeated to define pinned, roller, or mixed supports.
    Example:
      --prescribe 19 uy 0
    fixes the global Y displacement of node 19 to zero.

  --output vtk_output_directory
    Sets the VTK output directory. The main model file is written as:
      vtk_output_directory/BeamModelPart_0_0.vtk

  --help or -h
    Prints the compact usage line and exits.
**************************************************************************************************************
General project run:
  ./Bending_Torsion_Beam_Test_Case/BernoulliSolver \
    --project Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/ProjectParameters.json \
    --output Bending_Torsion_Beam_Test_Case/kratos_vtk_coarse
**************************************************************************************************************
CLI validation testcase 01: cantilever beam, free-end point load Fy = -1000.
  ./Bending_Torsion_Beam_Test_Case/BernoulliSolver \
    Bending_Torsion_Beam_Test_Case/TestCase_Input/01_cantilever_point_load/beam_geometry_coarse.mdpa \
    --materials Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json \
    --fix-node 19 \
    --load 1 0 -1000 0 0 0 0 \
    --output Bending_Torsion_Beam_Test_Case/BernoulliSolver_validation/cli_01_cantilever_point_load

CLI validation testcase 02: simply supported beam, center point load Fy = -1000.
  ./Bending_Torsion_Beam_Test_Case/BernoulliSolver \
    Bending_Torsion_Beam_Test_Case/TestCase_Input/02_simply_supported_point_load/beam_geometry_coarse.mdpa \
    --materials Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json \
    --prescribe 19 ux 0 --prescribe 19 uy 0 --prescribe 19 uz 0 --prescribe 19 rx 0 \
    --prescribe 1 uy 0 --prescribe 1 uz 0 --prescribe 1 rx 0 \
    --load 10 0 -1000 0 0 0 0 \
    --output Bending_Torsion_Beam_Test_Case/BernoulliSolver_validation/cli_02_simply_supported_point_load

CLI validation testcase 03: simply supported beam, uniform distributed load qy = -1000.
  This uses the stored testcase mdpa only as a mesh; all supports and equivalent
  nodal loads are specified below on the command line.
  ./Bending_Torsion_Beam_Test_Case/BernoulliSolver \
    Bending_Torsion_Beam_Test_Case/TestCase_Input/03_simply_supported_udl/simply_supported_udl.mdpa \
    --materials Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json \
    --prescribe 1 ux 0 --prescribe 1 uy 0 --prescribe 1 uz 0 --prescribe 1 rx 0 \
    --prescribe 11 uy 0 --prescribe 11 uz 0 --prescribe 11 rx 0 \
    --load 1 0 -500 0 0 0 -83.3333333333 \
    --load 2 0 -1000 0 0 0 0 --load 3 0 -1000 0 0 0 0 \
    --load 4 0 -1000 0 0 0 0 --load 5 0 -1000 0 0 0 0 \
    --load 6 0 -1000 0 0 0 0 --load 7 0 -1000 0 0 0 0 \
    --load 8 0 -1000 0 0 0 0 --load 9 0 -1000 0 0 0 0 \
    --load 10 0 -1000 0 0 0 0 \
    --load 11 0 -500 0 0 0 83.3333333333 \
    --output Bending_Torsion_Beam_Test_Case/BernoulliSolver_validation/cli_03_simply_supported_udl

CLI validation testcase 04: fixed-free torsion beam, free-end torque Mx = +1000.
  ./Bending_Torsion_Beam_Test_Case/BernoulliSolver \
    Bending_Torsion_Beam_Test_Case/TestCase_Input/04_torsion/beam_geometry_coarse.mdpa \
    --materials Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json \
    --fix-node 19 \
    --load 1 0 0 0 1000 0 0 \
    --output Bending_Torsion_Beam_Test_Case/BernoulliSolver_validation/cli_04_torsion

CLI validation testcase 05: fixed-free combined loading, Mx = +1000 and Mz = +500.
  ./Bending_Torsion_Beam_Test_Case/BernoulliSolver \
    Bending_Torsion_Beam_Test_Case/TestCase_Input/05_combined_bending_torsion/beam_geometry_coarse.mdpa \
    --materials Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json \
    --fix-node 19 \
    --load 1 0 0 0 1000 0 500 \
    --output Bending_Torsion_Beam_Test_Case/BernoulliSolver_validation/cli_05_combined_bending_torsion
*/

/************ 
Build and run from /root/dev/Kratos.

1. Compile:

Paste the following command in the terminal from /root/dev/Kratos:

clang++ -std=c++20 -O0 -g -DNDEBUG -DKRATOS_DEBUG -DBOOST_UBLAS_MOVE_SEMANTICS -DKRATOS_CORE=IMPORT,API -DKRATOS_MAJOR_VERSION=10 -DKRATOS_MINOR_VERSION=4 -DKRATOS_PYTHON -DKRATOS_SMP_OPENMP -DSTRUCTURAL_MECHANICS_APPLICATION=IMPORT,API -fopenmp=libomp Bending_Torsion_Beam_Test_Case/BernoulliSolver.cpp -o Bending_Torsion_Beam_Test_Case/BernoulliSolver -Ibuild/kratos -Ikratos -Iexternal_libraries/delaunator-cpp/include -Iexternal_libraries/tinyexpr -isystem /root/dev/boost_1_85_0 -isystem external_libraries -Iapplications/StructuralMechanicsApplication -Lbuild/kratos -Lbuild/applications/StructuralMechanicsApplication -lKratosStructuralMechanicsCore -lKratosCore -Wl,-rpath,/root/dev/Kratos/build/kratos -Wl,-rpath,/root/dev/Kratos/build/applications/StructuralMechanicsApplication

Note:
Use clang++ with libomp to match this Kratos build. Do not add
-DKRATOS_BUILD_TESTING here; it can instantiate unrelated test-only Kratos
templates and break this standalone compile.

2. Run:

./Bending_Torsion_Beam_Test_Case/BernoulliSolver \
  --project Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/ProjectParameters.json \
  --output Bending_Torsion_Beam_Test_Case/kratos_vtk_coarse
****************/
/*
OutputFile:
    Bending_Torsion_Beam_Test_Case/BernoulliBeamSolverOutput_vtk/cantilever_downward_load.vtk
*/



#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "includes/kernel.h"
#include "containers/model.h"
#include "includes/model_part_io.h"
#include "includes/variables.h"
#include "input_output/vtk_output.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "spaces/ublas_space.h"
#include "utilities/function_parser_utility.h"
#include "utilities/math_utils.h"

#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
namespace
{

constexpr std::size_t DofsPerNode = 6;

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using DenseSpaceType = UblasSpace<double, Matrix, Vector>;
using SolverType = SkylineLUFactorizationSolver<SparseSpaceType, DenseSpaceType>;
using SparseRows = std::vector<std::unordered_map<std::size_t, double>>;

struct ProjectSettings {
    std::filesystem::path ProjectParametersFile;
    std::filesystem::path BaseDirectory;
    std::filesystem::path MdpaFile;
    std::filesystem::path MaterialsFile;
    std::filesystem::path VtkOutputPath = "Bending_Torsion_Beam_Test_Case/BernoulliBeamSolverOutput_vtk";
    Parameters ParametersObject = Parameters(R"({})");
};

struct PrescribedDof {
    int NodeId = 0;
    std::string VariableName;
    std::size_t Component = 0;
    double Value = 0.0;
};

struct MaterialData {
    double E = 206.9e9;
    double Nu = 0.29;
    double A = 1.0;
    double Iy = 1.0;
    double Iz = 1.0;
    double J = 1.0;

    double G() const
    {
        return E / (2.0 * (1.0 + Nu));
    }
};

double Norm3(const array_1d<double, 3>& rV)
{
    return std::sqrt(rV[0] * rV[0] + rV[1] * rV[1] + rV[2] * rV[2]);
}

array_1d<double, 3> Normalize3(const array_1d<double, 3>& rV)
{
    const double norm = Norm3(rV);
    KRATOS_ERROR_IF(norm <= std::numeric_limits<double>::epsilon()) << "Zero length vector." << std::endl;
    return rV / norm;
}

void ApplyPropertiesMaterialValues(const Properties& rProperties, MaterialData& rData)
{
    if (rProperties.Has(YOUNG_MODULUS)) rData.E = rProperties[YOUNG_MODULUS];
    if (rProperties.Has(POISSON_RATIO)) rData.Nu = rProperties[POISSON_RATIO];
    if (rProperties.Has(CROSS_AREA)) rData.A = rProperties[CROSS_AREA];
    if (rProperties.Has(I22)) rData.Iy = rProperties[I22];
    if (rProperties.Has(I33)) rData.Iz = rProperties[I33];
    if (rProperties.Has(TORSIONAL_INERTIA)) rData.J = rProperties[TORSIONAL_INERTIA];
}

MaterialData ReadMaterial(const std::filesystem::path& rMaterialFile, const ModelPart& rModelPart)
{
    MaterialData data;
    if (rModelPart.NumberOfElements() > 0) {
        ApplyPropertiesMaterialValues(rModelPart.ElementsBegin()->GetProperties(), data);
    } else if (rModelPart.NumberOfProperties() > 0) {
        ApplyPropertiesMaterialValues(*rModelPart.PropertiesBegin(), data);
    }

    std::ifstream file(rMaterialFile);
    if (!file) {
        KRATOS_WARNING("BernoulliSolver") << "Material file not found: " << rMaterialFile
                                          << ". Using mdpa/default values." << std::endl;
        return data;
    }

    const std::string json((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    const Parameters parameters(json);
    const Parameters variables = parameters["properties"][0]["Material"]["Variables"];
    data.E = variables.Has("YOUNG_MODULUS") ? variables["YOUNG_MODULUS"].GetDouble() : data.E;
    data.Nu = variables.Has("POISSON_RATIO") ? variables["POISSON_RATIO"].GetDouble() : data.Nu;
    data.A = variables.Has("CROSS_AREA") ? variables["CROSS_AREA"].GetDouble() : data.A;
    data.Iy = variables.Has("I22") ? variables["I22"].GetDouble() : data.Iy;
    data.Iz = variables.Has("I33") ? variables["I33"].GetDouble() : data.Iz;
    data.J = variables.Has("TORSIONAL_INERTIA") ? variables["TORSIONAL_INERTIA"].GetDouble() : data.J;
    return data;
}

std::string ReadFileToString(const std::filesystem::path& rFileName)
{
    std::ifstream file(rFileName);
    KRATOS_ERROR_IF_NOT(file) << "Could not open file: " << rFileName << std::endl;
    return std::string((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
}

std::filesystem::path MakeAbsolutePath(
    const std::filesystem::path& rBaseDirectory,
    const std::filesystem::path& rPath)
{
    return rPath.is_absolute() ? rPath : rBaseDirectory / rPath;
}

bool PathExistsWithOptionalMdpaExtension(const std::filesystem::path& rPath)
{
    std::filesystem::path mdpa_path = rPath;
    mdpa_path.replace_extension(".mdpa");
    return std::filesystem::exists(rPath) || std::filesystem::exists(mdpa_path);
}

std::filesystem::path ResolveProjectInputPath(
    const std::filesystem::path& rProjectFileDirectory,
    const std::filesystem::path& rRelativePath)
{
    if (rRelativePath.is_absolute()) {
        return rRelativePath;
    }

    const std::filesystem::path path_from_project_directory = rProjectFileDirectory / rRelativePath;
    if (PathExistsWithOptionalMdpaExtension(path_from_project_directory)) {
        return path_from_project_directory;
    }

    const std::filesystem::path path_from_case_directory = rProjectFileDirectory.parent_path() / rRelativePath;
    if (PathExistsWithOptionalMdpaExtension(path_from_case_directory)) {
        return path_from_case_directory;
    }

    return path_from_project_directory;
}

std::filesystem::path WithoutMdpaExtension(std::filesystem::path path)
{
    if (path.extension() == ".mdpa") {
        path.replace_extension("");
    }
    return path;
}

ProjectSettings ReadProjectSettings(const std::filesystem::path& rProjectParametersFile)
{
    ProjectSettings settings;
    settings.ProjectParametersFile = rProjectParametersFile;
    settings.BaseDirectory = rProjectParametersFile.parent_path();
    settings.ParametersObject = Parameters(ReadFileToString(rProjectParametersFile));

    const Parameters solver_settings = settings.ParametersObject["solver_settings"];
    if (!solver_settings.Has("model_import_settings") &&
        solver_settings.Has("solvers") &&
        solver_settings["solvers"].Has("beam_structure")) {
        std::filesystem::path nested_project_file =
            solver_settings["solvers"]["beam_structure"]["solver_wrapper_settings"]["input_file"].GetString();
        if (nested_project_file.extension().empty()) {
            nested_project_file.replace_extension(".json");
        }
        return ReadProjectSettings(MakeAbsolutePath(settings.BaseDirectory, nested_project_file));
    }

    const std::string input_file_name =
        solver_settings["model_import_settings"]["input_filename"].GetString();
    const std::string materials_file_name =
        solver_settings["material_import_settings"]["materials_filename"].GetString();

    settings.MdpaFile = ResolveProjectInputPath(settings.BaseDirectory, input_file_name);
    settings.MaterialsFile = ResolveProjectInputPath(settings.BaseDirectory, materials_file_name);

    if (settings.ParametersObject.Has("output_processes") &&
        settings.ParametersObject["output_processes"].Has("vtk_output") &&
        settings.ParametersObject["output_processes"]["vtk_output"].size() > 0) {
        const Parameters vtk_parameters =
            settings.ParametersObject["output_processes"]["vtk_output"][0]["Parameters"];
        if (vtk_parameters.Has("output_path")) {
            settings.VtkOutputPath =
                MakeAbsolutePath(settings.BaseDirectory.parent_path(), vtk_parameters["output_path"].GetString());
        }
    }
    return settings;
}

ModelPart& ResolveModelPart(ModelPart& rRootModelPart, const std::string& rModelPartName)
{
    std::vector<std::string> names;
    std::stringstream name_stream(rModelPartName);
    std::string name;
    while (std::getline(name_stream, name, '.')) {
        if (!name.empty()) {
            names.push_back(name);
        }
    }

    std::size_t first_sub_model_part_index = 0;
    if (!names.empty() && (names.front() == rRootModelPart.Name() || names.front() == "Structure")) {
        first_sub_model_part_index = 1;
    }
    if (first_sub_model_part_index >= names.size()) {
        return rRootModelPart;
    }

    ModelPart* p_model_part = &rRootModelPart;
    for (std::size_t i = first_sub_model_part_index; i < names.size(); ++i) {
        KRATOS_ERROR_IF_NOT(p_model_part->HasSubModelPart(names[i]))
            << "SubModelPart '" << names[i] << "' not found while resolving '"
            << rModelPartName << "'." << std::endl;
        p_model_part = &p_model_part->GetSubModelPart(names[i]);
    }
    return *p_model_part;
}

double EvaluateScalarSettingAtNode(const Parameters rValue, const Node& rNode, const double Time)
{
    if (rValue.IsNumber()) {
        return rValue.GetDouble();
    }

    KRATOS_ERROR_IF_NOT(rValue.IsString())
        << "Vector process values must be numbers or function strings." << std::endl;

    GenericFunctionUtility function(rValue.GetString());
    return function.RotateAndCallFunction(
        rNode.X(), rNode.Y(), rNode.Z(), Time,
        rNode.X0(), rNode.Y0(), rNode.Z0());
}

void SetVectorComponent(Node& rNode, const std::string& rVariableName, const std::size_t Component, const double Value)
{
    if (rVariableName == "DISPLACEMENT") {
        rNode.FastGetSolutionStepValue(DISPLACEMENT)[Component] = Value;
    } else if (rVariableName == "ROTATION") {
        rNode.FastGetSolutionStepValue(ROTATION)[Component] = Value;
    } else if (rVariableName == "POINT_LOAD") {
        rNode.FastGetSolutionStepValue(POINT_LOAD)[Component] += Value;
    } else if (rVariableName == "POINT_MOMENT") {
        rNode.FastGetSolutionStepValue(POINT_MOMENT)[Component] += Value;
    } else {
        KRATOS_ERROR << "Unsupported vector variable in BernoulliSolver process: "
                     << rVariableName << std::endl;
    }
}

void FixVectorComponent(Node& rNode, const std::string& rVariableName, const std::size_t Component)
{
    if (rVariableName == "DISPLACEMENT") {
        if (Component == 0) rNode.Fix(DISPLACEMENT_X);
        if (Component == 1) rNode.Fix(DISPLACEMENT_Y);
        if (Component == 2) rNode.Fix(DISPLACEMENT_Z);
    } else if (rVariableName == "ROTATION") {
        if (Component == 0) rNode.Fix(ROTATION_X);
        if (Component == 1) rNode.Fix(ROTATION_Y);
        if (Component == 2) rNode.Fix(ROTATION_Z);
    } else {
        KRATOS_ERROR << "Only DISPLACEMENT and ROTATION can be constrained. Got: "
                     << rVariableName << std::endl;
    }
}

std::pair<std::string, std::size_t> ParseDofName(const std::string& rDofName)
{
    if (rDofName == "ux") return {"DISPLACEMENT", 0};
    if (rDofName == "uy") return {"DISPLACEMENT", 1};
    if (rDofName == "uz") return {"DISPLACEMENT", 2};
    if (rDofName == "rx") return {"ROTATION", 0};
    if (rDofName == "ry") return {"ROTATION", 1};
    if (rDofName == "rz") return {"ROTATION", 2};
    KRATOS_ERROR << "Unsupported dof name '" << rDofName
                 << "'. Use ux, uy, uz, rx, ry, or rz." << std::endl;
    return {"DISPLACEMENT", 0};
}

void FixAllNodeDofs(ModelPart& rModelPart, const int NodeId)
{
    Node& r_node = rModelPart.GetNode(NodeId);
    for (std::size_t component = 0; component < 3; ++component) {
        SetVectorComponent(r_node, "DISPLACEMENT", component, 0.0);
        FixVectorComponent(r_node, "DISPLACEMENT", component);
        SetVectorComponent(r_node, "ROTATION", component, 0.0);
        FixVectorComponent(r_node, "ROTATION", component);
    }
}

void PrescribeNodeDof(ModelPart& rModelPart, const PrescribedDof& rPrescribedDof)
{
    Node& r_node = rModelPart.GetNode(rPrescribedDof.NodeId);
    SetVectorComponent(
        r_node,
        rPrescribedDof.VariableName,
        rPrescribedDof.Component,
        rPrescribedDof.Value);
    FixVectorComponent(r_node, rPrescribedDof.VariableName, rPrescribedDof.Component);
}

Matrix LocalStiffness3DBeam(const MaterialData& rMaterial, const double L)
{
    Matrix k = ZeroMatrix(12, 12);

    const double EA_L = rMaterial.E * rMaterial.A / L;
    k(0, 0) = EA_L;
    k(0, 6) = -EA_L;
    k(6, 0) = -EA_L;
    k(6, 6) = EA_L;

    const double GJ_L = rMaterial.G() * rMaterial.J / L;
    k(3, 3) = GJ_L;
    k(3, 9) = -GJ_L;
    k(9, 3) = -GJ_L;
    k(9, 9) = GJ_L;

    auto add_block = [&k](const std::array<std::size_t, 4>& dofs, const Matrix& rBlock) {
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = 0; j < 4; ++j) {
                k(dofs[i], dofs[j]) += rBlock(i, j);
            }
        }
    };

    Matrix kb = ZeroMatrix(4, 4);
    double c1 = 12.0 * rMaterial.E * rMaterial.Iz / std::pow(L, 3);
    double c2 = 6.0 * rMaterial.E * rMaterial.Iz / std::pow(L, 2);
    double c3 = 4.0 * rMaterial.E * rMaterial.Iz / L;
    double c4 = 2.0 * rMaterial.E * rMaterial.Iz / L;
    kb(0, 0) = c1;  kb(0, 1) = c2;  kb(0, 2) = -c1; kb(0, 3) = c2;
    kb(1, 0) = c2;  kb(1, 1) = c3;  kb(1, 2) = -c2; kb(1, 3) = c4;
    kb(2, 0) = -c1; kb(2, 1) = -c2; kb(2, 2) = c1;  kb(2, 3) = -c2;
    kb(3, 0) = c2;  kb(3, 1) = c4;  kb(3, 2) = -c2; kb(3, 3) = c3;
    add_block({1, 5, 7, 11}, kb);

    c1 = 12.0 * rMaterial.E * rMaterial.Iy / std::pow(L, 3);
    c2 = 6.0 * rMaterial.E * rMaterial.Iy / std::pow(L, 2);
    c3 = 4.0 * rMaterial.E * rMaterial.Iy / L;
    c4 = 2.0 * rMaterial.E * rMaterial.Iy / L;
    kb = ZeroMatrix(4, 4);
    kb(0, 0) = c1;  kb(0, 1) = -c2; kb(0, 2) = -c1; kb(0, 3) = -c2;
    kb(1, 0) = -c2; kb(1, 1) = c3;  kb(1, 2) = c2;  kb(1, 3) = c4;
    kb(2, 0) = -c1; kb(2, 1) = c2;  kb(2, 2) = c1;  kb(2, 3) = c2;
    kb(3, 0) = -c2; kb(3, 1) = c4;  kb(3, 2) = c2;  kb(3, 3) = c3;
    add_block({2, 4, 8, 10}, kb);

    return k;
}

Matrix RotationMatrixLocalToGlobal(const Node& rNode1, const Node& rNode2)
{
    array_1d<double, 3> ex = rNode2.Coordinates() - rNode1.Coordinates();
    ex = Normalize3(ex);

    array_1d<double, 3> reference;
    reference[0] = 0.0; reference[1] = 0.0; reference[2] = 1.0;
    if (std::abs(inner_prod(ex, reference)) > 0.95) {
        reference[1] = 1.0;
        reference[2] = 0.0;
    }

    array_1d<double, 3> ey;
    MathUtils<double>::CrossProduct(ey, reference, ex);
    ey = Normalize3(ey);

    array_1d<double, 3> ez;
    MathUtils<double>::CrossProduct(ez, ex, ey);

    Matrix R = ZeroMatrix(3, 3);
    for (std::size_t i = 0; i < 3; ++i) {
        R(0, i) = ex[i];
        R(1, i) = ey[i];
        R(2, i) = ez[i];
    }
    return R;
}

Matrix TransformToGlobal(const Matrix& rLocalK, const Matrix& rR)
{
    Matrix T = ZeroMatrix(12, 12);
    for (std::size_t block = 0; block < 4; ++block) {
        const std::size_t offset = block * 3;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                T(offset + i, offset + j) = rR(i, j);
            }
        }
    }
    return prod(trans(T), Matrix(prod(rLocalK, T)));
}

std::size_t FullDof(const std::unordered_map<IndexType, std::size_t>& rNodeEquationIds, const IndexType NodeId, const std::size_t LocalDof)
{
    return rNodeEquationIds.at(NodeId) * DofsPerNode + LocalDof;
}

void AssembleStiffness(
    ModelPart& rModelPart,
    const MaterialData& rMaterial,
    const std::unordered_map<IndexType, std::size_t>& rNodeEquationIds,
    SparseRows& rK)
{
    rK.assign(rModelPart.NumberOfNodes() * DofsPerNode, {});

    for (auto& r_element : rModelPart.Elements()) {
        auto& r_geom = r_element.GetGeometry();
        KRATOS_ERROR_IF(r_geom.size() != 2) << "Only 2-node beam elements are supported." << std::endl;

        const double L = r_geom.Length();
        const Matrix local_k = LocalStiffness3DBeam(rMaterial, L);
        const Matrix global_k = TransformToGlobal(local_k, RotationMatrixLocalToGlobal(r_geom[0], r_geom[1]));

        std::array<std::size_t, 12> dofs;
        for (std::size_t a = 0; a < 2; ++a) {
            for (std::size_t d = 0; d < DofsPerNode; ++d) {
                dofs[a * DofsPerNode + d] = FullDof(rNodeEquationIds, r_geom[a].Id(), d);
            }
        }

        for (std::size_t i = 0; i < 12; ++i) {
            for (std::size_t j = 0; j < 12; ++j) {
                rK[dofs[i]][dofs[j]] += global_k(i, j);
            }
        }
    }
}

void AddCommandLineLoad(ModelPart& rModelPart, const int NodeId, const std::array<double, 6>& rLoad)
{
    auto& r_node = rModelPart.GetNode(NodeId);
    auto& r_force = r_node.FastGetSolutionStepValue(POINT_LOAD);
    auto& r_moment = r_node.FastGetSolutionStepValue(POINT_MOMENT);
    for (std::size_t i = 0; i < 3; ++i) {
        r_force[i] += rLoad[i];
        r_moment[i] += rLoad[i + 3];
    }
}

void ApplyAssignVectorVariableProcess(
    ModelPart& rRootModelPart,
    const Parameters rProcessParameters,
    const double Time,
    const bool FixConstrainedComponents)
{
    const std::string model_part_name = rProcessParameters["model_part_name"].GetString();
    const std::string variable_name = rProcessParameters["variable_name"].GetString();
    ModelPart& r_target_model_part = ResolveModelPart(rRootModelPart, model_part_name);

    const Parameters values = rProcessParameters["value"];
    KRATOS_ERROR_IF(values.size() != 3)
        << "Expected a 3-component value for " << variable_name << " in "
        << model_part_name << "." << std::endl;

    std::array<bool, 3> constrained = {false, false, false};
    if (FixConstrainedComponents && rProcessParameters.Has("constrained")) {
        const Parameters constrained_settings = rProcessParameters["constrained"];
        KRATOS_ERROR_IF(constrained_settings.size() != 3)
            << "Expected a 3-component constrained array for " << variable_name << "." << std::endl;
        for (std::size_t component = 0; component < 3; ++component) {
            constrained[component] = constrained_settings[component].GetBool();
        }
    }

    for (auto& r_node : r_target_model_part.Nodes()) {
        for (std::size_t component = 0; component < 3; ++component) {
            const double value = EvaluateScalarSettingAtNode(values[component], r_node, Time);
            SetVectorComponent(r_node, variable_name, component, value);
            if (constrained[component]) {
                FixVectorComponent(r_node, variable_name, component);
            }
        }
    }
}

void AddVectorToNode(Node& rNode, const std::string& rVariableName, const array_1d<double, 3>& rValue)
{
    if (rVariableName == "POINT_LOAD") {
        noalias(rNode.FastGetSolutionStepValue(POINT_LOAD)) += rValue;
    } else if (rVariableName == "POINT_MOMENT") {
        noalias(rNode.FastGetSolutionStepValue(POINT_MOMENT)) += rValue;
    } else {
        KRATOS_ERROR << "Unsupported load variable in BernoulliSolver: "
                     << rVariableName << std::endl;
    }
}

void ApplyVectorByDirectionToConditionProcess(
    ModelPart& rRootModelPart,
    const Parameters rProcessParameters)
{
    const std::string model_part_name = rProcessParameters["model_part_name"].GetString();
    const std::string variable_name = rProcessParameters["variable_name"].GetString();
    ModelPart& r_target_model_part = ResolveModelPart(rRootModelPart, model_part_name);

    const double modulus = rProcessParameters["modulus"].GetDouble();
    const Parameters direction_settings = rProcessParameters["direction"];
    KRATOS_ERROR_IF(direction_settings.size() != 3)
        << "Expected a 3-component direction for " << variable_name << " in "
        << model_part_name << "." << std::endl;

    array_1d<double, 3> value;
    for (std::size_t i = 0; i < 3; ++i) {
        value[i] = modulus * direction_settings[i].GetDouble();
    }

    for (auto& r_condition : r_target_model_part.Conditions()) {
        auto& r_geometry = r_condition.GetGeometry();
        const double share = 1.0 / static_cast<double>(r_geometry.size());
        for (auto& r_node : r_geometry) {
            array_1d<double, 3> nodal_value = share * value;
            AddVectorToNode(r_node, variable_name, nodal_value);
        }
    }
}

void ApplyProjectProcesses(
    ModelPart& rModelPart,
    const Parameters rProjectParameters)
{
    if (!rProjectParameters.Has("processes")) {
        return;
    }

    const Parameters process_settings = rProjectParameters["processes"];
    double process_time = 0.0;
    if (rProjectParameters.Has("problem_data") &&
        rProjectParameters["problem_data"].Has("end_time")) {
        process_time = rProjectParameters["problem_data"]["end_time"].GetDouble();
    }

    if (process_settings.Has("constraints_process_list")) {
        for (const auto& r_process : process_settings["constraints_process_list"]) {
            const std::string process_name = r_process["process_name"].GetString();
            if (process_name == "AssignVectorVariableProcess") {
                ApplyAssignVectorVariableProcess(rModelPart, r_process["Parameters"], process_time, true);
            } else {
                KRATOS_WARNING("BernoulliSolver")
                    << "Ignoring unsupported constraint process: " << process_name << std::endl;
            }
        }
    }

    if (process_settings.Has("loads_process_list")) {
        for (const auto& r_process : process_settings["loads_process_list"]) {
            const std::string process_name = r_process["process_name"].GetString();
            if (process_name == "AssignVectorByDirectionToConditionProcess") {
                ApplyVectorByDirectionToConditionProcess(rModelPart, r_process["Parameters"]);
            } else if (process_name == "AssignVectorVariableProcess") {
                ApplyAssignVectorVariableProcess(rModelPart, r_process["Parameters"], process_time, false);
            } else {
                KRATOS_WARNING("BernoulliSolver")
                    << "Ignoring unsupported load process: " << process_name << std::endl;
            }
        }
    }
}

void ApplyDomainSizeConstraints(ModelPart& rModelPart, const Parameters rProjectParameters)
{
    if (!rProjectParameters.Has("solver_settings") ||
        !rProjectParameters["solver_settings"].Has("domain_size")) {
        return;
    }

    const int domain_size = rProjectParameters["solver_settings"]["domain_size"].GetInt();
    if (domain_size == 2) {
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.Fix(DISPLACEMENT_Z);
            r_node.Fix(ROTATION_Y);
        }
    }
}

bool IsFixedDof(const Node& rNode, const std::size_t LocalDof)
{
    if (LocalDof == 0) return rNode.IsFixed(DISPLACEMENT_X);
    if (LocalDof == 1) return rNode.IsFixed(DISPLACEMENT_Y);
    if (LocalDof == 2) return rNode.IsFixed(DISPLACEMENT_Z);
    if (LocalDof == 3) return rNode.IsFixed(ROTATION_X);
    if (LocalDof == 4) return rNode.IsFixed(ROTATION_Y);
    return rNode.IsFixed(ROTATION_Z);
}

double GetDofValue(const Node& rNode, const std::size_t LocalDof)
{
    if (LocalDof < 3) {
        return rNode.FastGetSolutionStepValue(DISPLACEMENT)[LocalDof];
    }
    return rNode.FastGetSolutionStepValue(ROTATION)[LocalDof - 3];
}

double GetDofLoad(const Node& rNode, const std::size_t LocalDof)
{
    if (LocalDof < 3) {
        return rNode.FastGetSolutionStepValue(POINT_LOAD)[LocalDof];
    }
    return rNode.FastGetSolutionStepValue(POINT_MOMENT)[LocalDof - 3];
}

void SetDofValue(Node& rNode, const std::size_t LocalDof, const double Value)
{
    if (LocalDof < 3) {
        rNode.FastGetSolutionStepValue(DISPLACEMENT)[LocalDof] = Value;
    } else {
        rNode.FastGetSolutionStepValue(ROTATION)[LocalDof - 3] = Value;
    }
}

Vector SolveReducedSystem(
    ModelPart& rModelPart,
    const SparseRows& rK,
    const std::unordered_map<IndexType, std::size_t>& rNodeEquationIds)
{
    const std::size_t full_size = rModelPart.NumberOfNodes() * DofsPerNode;
    std::vector<int> full_to_free(full_size, -1);
    Vector fixed_values(full_size, 0.0);
    Vector rhs_full(full_size, 0.0);
    std::size_t n_free = 0;

    for (auto& r_node : rModelPart.Nodes()) {
        for (std::size_t d = 0; d < DofsPerNode; ++d) {
            const std::size_t equation_id = FullDof(rNodeEquationIds, r_node.Id(), d);
            rhs_full[equation_id] = GetDofLoad(r_node, d);
            fixed_values[equation_id] = GetDofValue(r_node, d);
            if (!IsFixedDof(r_node, d)) {
                full_to_free[equation_id] = static_cast<int>(n_free++);
            }
        }
    }

    if (n_free == 0) {
        KRATOS_INFO("BernoulliSolver")
            << "All DOFs are constrained by the JSON settings; writing prescribed values." << std::endl;
        return fixed_values;
    }

    CompressedMatrix Kff(n_free, n_free);
    Vector Ff(n_free, 0.0);
    for (std::size_t i_full = 0; i_full < full_size; ++i_full) {
        const int i_free = full_to_free[i_full];
        if (i_free < 0) {
            continue;
        }

        Ff[i_free] = rhs_full[i_full];
        for (const auto& [j_full, value] : rK[i_full]) {
            const int j_free = full_to_free[j_full];
            if (j_free >= 0) {
                Kff.insert_element(i_free, j_free, value);
            } else {
                Ff[i_free] -= value * fixed_values[j_full];
            }
        }
    }

    Vector Uf(n_free, 0.0);
    SolverType solver;
    solver.Solve(Kff, Uf, Ff);

    Vector U(full_size, 0.0);
    for (std::size_t i_full = 0; i_full < full_size; ++i_full) {
        U[i_full] = full_to_free[i_full] >= 0 ? Uf[full_to_free[i_full]] : fixed_values[i_full];
    }

    KRATOS_INFO("BernoulliSolver") << "Solved reduced system with " << n_free << " free DOFs." << std::endl;
    return U;
}

void WriteSolutionToNodes(
    ModelPart& rModelPart,
    const Vector& rU,
    const std::unordered_map<IndexType, std::size_t>& rNodeEquationIds)
{
    for (auto& r_node : rModelPart.Nodes()) {
        for (std::size_t d = 0; d < DofsPerNode; ++d) {
            SetDofValue(r_node, d, rU[FullDof(rNodeEquationIds, r_node.Id(), d)]);
        }
    }
}

void WriteVtk(ModelPart& rModelPart, const std::filesystem::path& rOutputPath)
{
    Parameters vtk_parameters(R"({
        "model_part_name" : "BeamModelPart",
        "output_path" : "Bending_Torsion_Beam_Test_Case/BernoulliBeamSolverOutput_vtk",
        "file_format" : "ascii",
        "output_precision" : 7,
        "output_sub_model_parts" : true,
        "write_deformed_configuration" : true,
        "nodal_solution_step_data_variables" : ["DISPLACEMENT", "ROTATION"],
        "nodal_data_value_variables" : [],
        "element_data_value_variables" : [],
        "condition_data_value_variables" : []
    })");
    vtk_parameters["output_path"].SetString(rOutputPath.string());
    VtkOutput(rModelPart, vtk_parameters).PrintOutput();
}

} // namespace
} // namespace Kratos

int main(int argc, char** argv)
{
    using namespace Kratos;

    try {
        Kernel kernel;
        auto p_structural_mechanics_application = Kratos::make_shared<KratosStructuralMechanicsApplication>();
        kernel.ImportApplication(p_structural_mechanics_application);

        std::filesystem::path project_parameters_file =
            "Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/ProjectParameters.json";
        std::filesystem::path mdpa_file_override;
        std::filesystem::path materials_file_override;
        std::filesystem::path output_path_override;
        std::vector<std::pair<int, std::array<double, 6>>> command_line_loads;
        std::vector<int> fixed_node_overrides;
        std::vector<PrescribedDof> prescribed_dof_overrides;
        bool project_parameters_were_provided = false;

        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--project" && i + 1 < argc) {
                project_parameters_file = argv[++i];
                project_parameters_were_provided = true;
            } else if (arg == "--materials" && i + 1 < argc) {
                materials_file_override = argv[++i];
            } else if (arg == "--output" && i + 1 < argc) {
                output_path_override = argv[++i];
            } else if (arg == "--load" && i + 7 < argc) {
                const int node_id = std::stoi(argv[++i]);
                std::array<double, 6> load{};
                for (double& r_value : load) {
                    r_value = std::stod(argv[++i]);
                }
                command_line_loads.emplace_back(node_id, load);
            } else if (arg == "--fix-node" && i + 1 < argc) {
                fixed_node_overrides.push_back(std::stoi(argv[++i]));
            } else if (arg == "--prescribe" && i + 3 < argc) {
                PrescribedDof prescribed_dof;
                prescribed_dof.NodeId = std::stoi(argv[++i]);
                const auto [variable_name, component] = ParseDofName(argv[++i]);
                prescribed_dof.VariableName = variable_name;
                prescribed_dof.Component = component;
                prescribed_dof.Value = std::stod(argv[++i]);
                prescribed_dof_overrides.push_back(prescribed_dof);
            } else if (arg == "--help" || arg == "-h") {
                std::cout << "Usage: " << argv[0] << " [ProjectParameters.json|mesh.mdpa] "
                          << "[--project ProjectParameters.json] "
                          << "[--materials StructuralMaterials.json] "
                          << "[--load node fx fy fz mx my mz] "
                          << "[--fix-node node] "
                          << "[--prescribe node ux|uy|uz|rx|ry|rz value] "
                          << "[--output vtk_output]\n";
                return 0;
            } else if (arg.rfind("--", 0) == 0) {
                KRATOS_ERROR << "Unknown option: " << arg << std::endl;
            } else if (std::filesystem::path(arg).extension() == ".json") {
                project_parameters_file = arg;
                project_parameters_were_provided = true;
            } else {
                mdpa_file_override = arg;
            }
        }

        ProjectSettings project_settings;
        if (!mdpa_file_override.empty() && !project_parameters_were_provided) {
            project_settings.ProjectParametersFile = project_parameters_file;
            project_settings.BaseDirectory = std::filesystem::current_path();
            project_settings.MdpaFile = mdpa_file_override;
            project_settings.MaterialsFile =
                materials_file_override.empty()
                    ? std::filesystem::path("Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json")
                    : materials_file_override;
        } else {
            project_settings = ReadProjectSettings(project_parameters_file);
        }
        const std::filesystem::path mdpa_file =
            mdpa_file_override.empty() ? project_settings.MdpaFile : mdpa_file_override;
        const std::filesystem::path materials_file =
            materials_file_override.empty() ? project_settings.MaterialsFile : materials_file_override;
        const std::filesystem::path output_path =
            output_path_override.empty() ? project_settings.VtkOutputPath : output_path_override;

        Model model;
        ModelPart& r_model_part = model.CreateModelPart("BeamModelPart");
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(ROTATION);
        r_model_part.AddNodalSolutionStepVariable(POINT_LOAD);
        r_model_part.AddNodalSolutionStepVariable(POINT_MOMENT);

        ModelPartIO(WithoutMdpaExtension(mdpa_file).string()).ReadModelPart(r_model_part);
        ApplyDomainSizeConstraints(r_model_part, project_settings.ParametersObject);
        ApplyProjectProcesses(r_model_part, project_settings.ParametersObject);
        for (const int node_id : fixed_node_overrides) {
            FixAllNodeDofs(r_model_part, node_id);
        }
        for (const auto& r_prescribed_dof : prescribed_dof_overrides) {
            PrescribeNodeDof(r_model_part, r_prescribed_dof);
        }
        for (const auto& [node_id, load] : command_line_loads) {
            AddCommandLineLoad(r_model_part, node_id, load);
        }

        std::unordered_map<IndexType, std::size_t> node_equation_ids;
        std::size_t node_counter = 0;
        for (const auto& r_node : r_model_part.Nodes()) {
            node_equation_ids[r_node.Id()] = node_counter++;
        }

        SparseRows K;
        AssembleStiffness(r_model_part, ReadMaterial(materials_file, r_model_part), node_equation_ids, K);
        const Vector U = SolveReducedSystem(r_model_part, K, node_equation_ids);
        WriteSolutionToNodes(r_model_part, U, node_equation_ids);
        WriteVtk(r_model_part, output_path);

        KRATOS_INFO("BernoulliSolver") << "Read " << r_model_part.NumberOfNodes() << " nodes and "
                                       << r_model_part.NumberOfElements() << " elements." << std::endl;
    } catch (const std::exception& r_error) {
        std::cerr << r_error.what() << std::endl;
        return 1;
    }
    return 0;
}
