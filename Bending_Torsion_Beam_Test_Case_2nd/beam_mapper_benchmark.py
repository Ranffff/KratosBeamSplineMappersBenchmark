"""Curl-consistent benchmark for beam mappers with rotational recovery.

This script implements the revised benchmark described in
BeamSplineDocumentation/Benchmark/rotation_recovery_benchmark_note.pdf.

The benchmark is a pure mapping test.  An analytical three-dimensional
displacement field g(x,y,z) is prescribed first.  The beam displacement is
g(x_b), the beam rotation is 0.5 curl(g)(x_b), and the mapped fluid-interface
displacement is compared against g(x_f).
"""

from pathlib import Path
import csv
import json
import math
import os
import shutil

os.environ.setdefault("MPLCONFIGDIR", "/tmp/kratos_documented_beam_mapper_benchmark_mpl")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import KratosMultiphysics as KM
import KratosMultiphysics.LinearSolversApplication
import KratosMultiphysics.MappingApplication
import KratosMultiphysics.StructuralMechanicsApplication


ROOT = Path(__file__).resolve().parent
OUTPUT_ROOT = ROOT / "TestCase_Output" / "Consistent_Beam_Mapper_Fluid_Convergence"
PLOT_ROOT = OUTPUT_ROOT / "Plots"
VTK_ROOT = OUTPUT_ROOT / "VTK"
GRID_COLOR = "#D9D9D9"

LENGTH = 1.0
WIDTH = 0.20
TRANSLATION_AMPLITUDE = 1.0e-3
ROTATION_AMPLITUDE = 1.0e-2

FINE_BEAM_ELEMENTS = 128
FLUID_REFINEMENT_DIVISIONS = (
    {"xi": 2, "eta": 10},
    {"xi": 4, "eta": 20},
    {"xi": 8, "eta": 40},
    {"xi": 16, "eta": 80},
    {"xi": 32, "eta": 160},
)

MAPPERS = (
    {
        "key": "beam_spline_no_rotational_recovery",
        "label": "Beam spline, no recovery",
        "settings": {
            "mapper_type": "beam_spline_mapper",
            "search_settings": {
                "search_radius": 0.30,
                "max_num_search_iterations": 30,
            },
            "local_coord_tolerance": 0.25,
            "echo_level": 0,
        },
        "color": "#0072B2",
        "marker": "o",
        "linestyle": "--",
    },
    {
        "key": "beam_spline_with_recovery_of_rotations",
        "label": "Beam spline, with recovery",
        "settings": {
            "mapper_type": "beam_spline_mapper_with_recovery_of_rotations",
            "search_settings": {
                "search_radius": 0.30,
                "max_num_search_iterations": 30,
            },
            "local_coord_tolerance": 0.25,
            "kernel_type": "gaussian",
            "kernel_radius": 0.50,
            "polynomial_basis": "high_order_line_adapted",
            "regularization": 1.0e-10,
            "echo_level": 0,
        },
        "color": "#009E73",
        "marker": "s",
        "linestyle": "-",
    },
    {
        "key": "beam_mapper_corotation",
        "label": "Beam mapper, co-rotation",
        "settings": {
            "mapper_type": "beam_mapper",
            "use_corotation": True,
            "search_settings": {
                "search_radius": 0.30,
                "max_num_search_iterations": 30,
            },
            "echo_level": 0,
        },
        "color": "#D55E00",
        "marker": "^",
        "linestyle": "-.",
    },
)

BENCHMARKS = (
    {
        "key": "rigid_body_rotation",
        "label": "Rigid-body rotation",
        "description": "g = [omega_y z, 0, u0 - omega_y x]",
    },
    {
        "key": "polynomial_recovery_field",
        "label": "Polynomial recovery field",
        "description": "g_z = a y^2 + b x",
    },
    {
        "key": "curl_compatible_beam_field",
        "label": "Curl-compatible beam field",
        "description": "g = [theta(y) z, 0, w(y) - theta(y) x]",
    },
)


def array3(values):
    return KM.Array3([float(values[0]), float(values[1]), float(values[2])])


def vector_to_tuple(value):
    return (float(value[0]), float(value[1]), float(value[2]))


def analytical_displacement(case_key, x, y, z):
    if case_key == "rigid_body_rotation":
        omega = (0.0, ROTATION_AMPLITUDE, 0.0)
        translation = (0.0, 0.0, TRANSLATION_AMPLITUDE)
        return (
            translation[0] + omega[1] * z - omega[2] * y,
            translation[1] + omega[2] * x - omega[0] * z,
            translation[2] + omega[0] * y - omega[1] * x,
        )
    if case_key == "polynomial_recovery_field":
        a = TRANSLATION_AMPLITUDE / (LENGTH * LENGTH)
        b = ROTATION_AMPLITUDE * WIDTH
        return (0.0, 0.0, a * y * y + b * x)
    if case_key == "curl_compatible_beam_field":
        w = TRANSLATION_AMPLITUDE * math.sin(math.pi * y / LENGTH)
        theta = ROTATION_AMPLITUDE * math.sin(math.pi * y / LENGTH)
        return (theta * z, 0.0, w - theta * x)
    raise RuntimeError(f"Unknown benchmark case: {case_key}")


def analytical_rotation(case_key, x, y, z):
    if case_key == "rigid_body_rotation":
        return (0.0, ROTATION_AMPLITUDE, 0.0)
    if case_key == "polynomial_recovery_field":
        a = TRANSLATION_AMPLITUDE / (LENGTH * LENGTH)
        b = ROTATION_AMPLITUDE * WIDTH
        return (a * y, -0.5 * b, 0.0)
    if case_key == "curl_compatible_beam_field":
        w_prime = (
            TRANSLATION_AMPLITUDE
            * math.pi
            / LENGTH
            * math.cos(math.pi * y / LENGTH)
        )
        theta = ROTATION_AMPLITUDE * math.sin(math.pi * y / LENGTH)
        theta_prime = (
            ROTATION_AMPLITUDE
            * math.pi
            / LENGTH
            * math.cos(math.pi * y / LENGTH)
        )
        return (
            0.5 * (w_prime - theta_prime * x),
            theta,
            -0.5 * theta_prime * z,
        )
    raise RuntimeError(f"Unknown benchmark case: {case_key}")


def numerical_half_curl(case_key, x, y, z, step=1.0e-6):
    def component(component_index, x_value, y_value, z_value):
        return analytical_displacement(case_key, x_value, y_value, z_value)[component_index]

    duz_dy = (component(2, x, y + step, z) - component(2, x, y - step, z)) / (2.0 * step)
    duy_dz = (component(1, x, y, z + step) - component(1, x, y, z - step)) / (2.0 * step)
    dux_dz = (component(0, x, y, z + step) - component(0, x, y, z - step)) / (2.0 * step)
    duz_dx = (component(2, x + step, y, z) - component(2, x - step, y, z)) / (2.0 * step)
    duy_dx = (component(1, x + step, y, z) - component(1, x - step, y, z)) / (2.0 * step)
    dux_dy = (component(0, x, y + step, z) - component(0, x, y - step, z)) / (2.0 * step)

    return (
        0.5 * (duz_dy - duy_dz),
        0.5 * (dux_dz - duz_dx),
        0.5 * (duy_dx - dux_dy),
    )


def validate_rotation_consistency():
    sample_points = (
        (-0.5 * WIDTH, 0.25 * LENGTH, 0.0),
        (0.0, 0.50 * LENGTH, 0.0),
        (0.5 * WIDTH, 0.75 * LENGTH, 0.0),
        (0.25 * WIDTH, 0.40 * LENGTH, 0.15 * WIDTH),
    )
    max_difference = 0.0
    for case in BENCHMARKS:
        for x, y, z in sample_points:
            analytical = analytical_rotation(case["key"], x, y, z)
            numerical = numerical_half_curl(case["key"], x, y, z)
            difference = math.sqrt(sum((analytical[i] - numerical[i]) ** 2 for i in range(3)))
            max_difference = max(max_difference, difference)
            if difference > 1.0e-8:
                raise RuntimeError(
                    "Rotation field is not consistent with theta = 0.5 curl(g): "
                    f"case={case['key']}, point=({x}, {y}, {z}), "
                    f"theta={analytical}, 0.5curl={numerical}, diff={difference:.3e}"
                )
    return max_difference


def displacement_and_rotation(case_key, s):
    x = 0.0
    y = s
    z = 0.0
    displacement = analytical_displacement(case_key, x, y, z)
    rotation = analytical_rotation(case_key, x, y, z)
    return displacement, rotation


def exact_surface_displacement(case_key, xi, eta):
    return analytical_displacement(case_key, xi, eta, 0.0)


def create_beam_model_part(case_key, number_of_elements):
    model = KM.Model()
    beam = model.CreateModelPart(f"beam_{case_key}_{number_of_elements}")
    beam.ProcessInfo[KM.DOMAIN_SIZE] = 3
    beam.ProcessInfo[KM.TIME] = 0.0
    beam.ProcessInfo[KM.DELTA_TIME] = 1.0
    for variable in (KM.DISPLACEMENT, KM.ROTATION):
        beam.AddNodalSolutionStepVariable(variable)

    properties = beam.CreateNewProperties(1)
    for i in range(number_of_elements + 1):
        s = LENGTH * i / number_of_elements
        node = beam.CreateNewNode(i + 1, 0.0, s, 0.0)
        displacement, rotation = displacement_and_rotation(case_key, s)
        node.SetSolutionStepValue(KM.DISPLACEMENT, array3(displacement))
        node.SetSolutionStepValue(KM.ROTATION, array3(rotation))

    for i in range(number_of_elements):
        beam.CreateNewElement("CrBeamElement3D2N", i + 1, [i + 1, i + 2], properties)

    return beam


def create_surface_model_part(case_key, xi_divisions, eta_divisions):
    model = KM.Model()
    surface = model.CreateModelPart(f"surface_{case_key}_{xi_divisions}x{eta_divisions}")
    surface.ProcessInfo[KM.DOMAIN_SIZE] = 3
    surface.ProcessInfo[KM.TIME] = 0.0
    surface.ProcessInfo[KM.DELTA_TIME] = 1.0
    surface.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    properties = surface.CreateNewProperties(1)

    for j in range(eta_divisions + 1):
        eta = LENGTH * j / eta_divisions
        for i in range(xi_divisions + 1):
            xi = -0.5 * WIDTH + WIDTH * i / xi_divisions
            node_id = j * (xi_divisions + 1) + i + 1
            surface.CreateNewNode(node_id, xi, eta, 0.0)

    element_id = 1
    for j in range(eta_divisions):
        row = j * (xi_divisions + 1)
        next_row = (j + 1) * (xi_divisions + 1)
        for i in range(xi_divisions):
            surface.CreateNewElement(
                "ShellThinElementCorotational3D4N",
                element_id,
                [row + i + 1, row + i + 2, next_row + i + 2, next_row + i + 1],
                properties,
            )
            element_id += 1

    return surface


def create_mapper(beam, surface, mapper_definition):
    return KM.MapperFactory.CreateMapper(
        beam,
        surface,
        KM.Parameters(json.dumps(mapper_definition["settings"])),
    )


def relative_discrete_l2_error(surface, case_key):
    error_squared = 0.0
    reference_squared = 0.0
    max_error = 0.0
    for node in surface.Nodes:
        mapped = vector_to_tuple(node.GetSolutionStepValue(KM.DISPLACEMENT))
        exact = exact_surface_displacement(case_key, node.X0, node.Y0)
        difference = tuple(mapped[i] - exact[i] for i in range(3))
        error_norm_squared = sum(value * value for value in difference)
        exact_norm_squared = sum(value * value for value in exact)
        error_squared += error_norm_squared
        reference_squared += exact_norm_squared
        max_error = max(max_error, math.sqrt(error_norm_squared))

    if reference_squared <= 0.0:
        raise RuntimeError(f"Zero exact displacement norm for case {case_key}")
    number_of_nodes = surface.NumberOfNodes()
    relative_l2_error = math.sqrt(error_squared) / math.sqrt(reference_squared)
    rmse_error = math.sqrt(error_squared / number_of_nodes)
    return relative_l2_error, rmse_error, max_error


def characteristic_fluid_size(xi_divisions, eta_divisions):
    return max(WIDTH / xi_divisions, LENGTH / eta_divisions)


def write_surface_vtk(surface, case_key, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    nodes = sorted(surface.Nodes, key=lambda node: node.Id)
    elements = sorted(surface.Elements, key=lambda element: element.Id)
    node_index = {node.Id: index for index, node in enumerate(nodes)}

    with path.open("w") as output:
        output.write("# vtk DataFile Version 3.0\n")
        output.write("Documented beam mapper benchmark\n")
        output.write("ASCII\n")
        output.write("DATASET UNSTRUCTURED_GRID\n")
        output.write(f"POINTS {len(nodes)} float\n")
        for node in nodes:
            output.write(f"{node.X0:.16e} {node.Y0:.16e} {node.Z0:.16e}\n")

        cell_size = sum(len(element.GetNodes()) + 1 for element in elements)
        output.write(f"CELLS {len(elements)} {cell_size}\n")
        for element in elements:
            ids = [node_index[node.Id] for node in element.GetNodes()]
            output.write(str(len(ids)) + " " + " ".join(str(i) for i in ids) + "\n")

        output.write(f"CELL_TYPES {len(elements)}\n")
        for _ in elements:
            output.write("9\n")

        output.write(f"POINT_DATA {len(nodes)}\n")
        for name, getter in (
            ("mapped_displacement", lambda node: vector_to_tuple(node.GetSolutionStepValue(KM.DISPLACEMENT))),
            ("exact_displacement", lambda node: exact_surface_displacement(case_key, node.X0, node.Y0)),
            (
                "displacement_error",
                lambda node: tuple(
                    vector_to_tuple(node.GetSolutionStepValue(KM.DISPLACEMENT))[i]
                    - exact_surface_displacement(case_key, node.X0, node.Y0)[i]
                    for i in range(3)
                ),
            ),
        ):
            output.write(f"VECTORS {name} float\n")
            for node in nodes:
                values = getter(node)
                output.write(f"{values[0]:.16e} {values[1]:.16e} {values[2]:.16e}\n")


def run_single_case(study_key, case, mapper_definition, beam_elements, fluid_divisions):
    beam = create_beam_model_part(case["key"], beam_elements)
    surface = create_surface_model_part(
        case["key"],
        fluid_divisions["xi"],
        fluid_divisions["eta"],
    )
    mapper = create_mapper(beam, surface, mapper_definition)
    mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT)

    relative_l2_error, rmse_error, max_error = relative_discrete_l2_error(surface, case["key"])
    characteristic_h = characteristic_fluid_size(fluid_divisions["xi"], fluid_divisions["eta"])
    refinement_label = f"fluid_{fluid_divisions['xi']}x{fluid_divisions['eta']}"

    vtk_path = (
        VTK_ROOT
        / study_key
        / case["key"]
        / mapper_definition["key"]
        / f"{refinement_label}.vtk"
    )
    write_surface_vtk(surface, case["key"], vtk_path)

    return {
        "study": study_key,
        "case": case["key"],
        "case_label": case["label"],
        "case_description": case["description"],
        "mapper": mapper_definition["key"],
        "mapper_label": mapper_definition["label"],
        "beam_elements": beam_elements,
        "fluid_xi_divisions": fluid_divisions["xi"],
        "fluid_eta_divisions": fluid_divisions["eta"],
        "surface_nodes": surface.NumberOfNodes(),
        "surface_elements": surface.NumberOfElements(),
        "characteristic_h": characteristic_h,
        "relative_l2_error": relative_l2_error,
        "rmse_displacement_error": rmse_error,
        "max_displacement_error": max_error,
        "vtk_file": str(vtk_path),
    }


def run_benchmark():
    results = []
    for case in BENCHMARKS:
        for fluid_divisions in FLUID_REFINEMENT_DIVISIONS:
            for mapper_definition in MAPPERS:
                print(
                    "Running "
                    f"{case['key']} | {mapper_definition['key']} | "
                    f"beam={FINE_BEAM_ELEMENTS}, "
                    f"fluid={fluid_divisions['xi']}x{fluid_divisions['eta']}",
                    flush=True,
                )
                results.append(
                    run_single_case(
                        "fluid_refinement",
                        case,
                        mapper_definition,
                        FINE_BEAM_ELEMENTS,
                        fluid_divisions,
                    )
                )
    return results


def write_results(results):
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_ROOT / "beam_mapper_benchmark_results.csv"
    json_path = OUTPUT_ROOT / "beam_mapper_benchmark_results.json"
    fieldnames = (
        "study",
        "case",
        "case_label",
        "mapper",
        "mapper_label",
        "beam_elements",
        "fluid_xi_divisions",
        "fluid_eta_divisions",
        "surface_nodes",
        "surface_elements",
        "characteristic_h",
        "relative_l2_error",
        "rmse_displacement_error",
        "max_displacement_error",
        "vtk_file",
    )
    with csv_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow({name: result[name] for name in fieldnames})
    with json_path.open("w") as output:
        json.dump(results, output, indent=2)
    return csv_path, json_path


def configure_plot_style():
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "lines.linewidth": 1.7,
            "lines.markersize": 5.2,
            "figure.dpi": 120,
            "savefig.dpi": 300,
        }
    )


def save_figure(figure, output_path):
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    figure.savefig(output_path.with_suffix(".svg"), bbox_inches="tight")


def create_convergence_plot(results, study_key, basename, metric_key, ylabel):
    selected = [result for result in results if result["study"] == study_key]
    figure, axes = plt.subplots(1, len(BENCHMARKS), figsize=(14.0, 4.1), sharey=False)
    for axis, case in zip(axes, BENCHMARKS):
        case_results = [result for result in selected if result["case"] == case["key"]]
        for mapper_definition in MAPPERS:
            curve = [
                result
                for result in case_results
                if result["mapper"] == mapper_definition["key"]
            ]
            curve = sorted(curve, key=lambda result: result["characteristic_h"], reverse=True)
            axis.loglog(
                [result["characteristic_h"] for result in curve],
                [max(result[metric_key], 1.0e-18) for result in curve],
                color=mapper_definition["color"],
                marker=mapper_definition["marker"],
                linestyle=mapper_definition["linestyle"],
                label=mapper_definition["label"],
            )
        axis.set_title(case["label"])
        axis.set_xlabel("Characteristic mesh size h")
        axis.grid(True, which="both", color=GRID_COLOR, linewidth=0.6)
        axis.invert_xaxis()
    axes[0].set_ylabel(ylabel)
    axes[-1].legend(frameon=True, framealpha=0.95, loc="best")
    figure.tight_layout()
    PLOT_ROOT.mkdir(parents=True, exist_ok=True)
    save_figure(figure, PLOT_ROOT / f"{basename}.png")
    plt.close(figure)


def create_plots(results):
    configure_plot_style()
    create_convergence_plot(
        results,
        "fluid_refinement",
        "relative_l2_error_fluid_refinement",
        "relative_l2_error",
        "Relative discrete L2 error",
    )
    create_convergence_plot(
        results,
        "fluid_refinement",
        "rmse_displacement_error_fluid_refinement",
        "rmse_displacement_error",
        "RMSE displacement error",
    )
    create_convergence_plot(
        results,
        "fluid_refinement",
        "max_displacement_error_fluid_refinement",
        "max_displacement_error",
        "Maximum displacement error",
    )


def validate_outputs(results, csv_path, json_path):
    expected = len(BENCHMARKS) * len(MAPPERS) * (
        len(FLUID_REFINEMENT_DIVISIONS)
    )
    if len(results) != expected:
        raise RuntimeError(f"Expected {expected} results, got {len(results)}")
    for path in (csv_path, json_path):
        if not path.is_file():
            raise RuntimeError(f"Missing result file: {path}")
    for basename in (
        "relative_l2_error_fluid_refinement",
        "rmse_displacement_error_fluid_refinement",
        "max_displacement_error_fluid_refinement",
    ):
        for extension in ("png", "svg"):
            path = PLOT_ROOT / f"{basename}.{extension}"
            if not path.is_file():
                raise RuntimeError(f"Missing plot: {path}")
    missing_vtk = [result["vtk_file"] for result in results if not Path(result["vtk_file"]).is_file()]
    if missing_vtk:
        raise RuntimeError("Missing VTK files: " + ", ".join(missing_vtk[:5]))


def print_summary(results):
    print("")
    print("Documented beam mapper benchmark")
    print(f"L = {LENGTH:g} m, c = {WIDTH:g} m, A = {TRANSLATION_AMPLITUDE:g} m, alpha = {ROTATION_AMPLITUDE:g} rad")
    print(f"Fixed fine beam: {FINE_BEAM_ELEMENTS} elements, h_beam = {LENGTH / FINE_BEAM_ELEMENTS:.3e}")
    print("Relative L2 error = sqrt(sum_i ||u_map_i - u_ref_i||^2) / sqrt(sum_i ||u_ref_i||^2)")
    print(f"{'study':<18} {'case':<36} {'mapper':<38} {'h_f':>10} {'rel L2':>12} {'RMSE':>12} {'max err':>12}")
    print("-" * 146)
    for result in results:
        print(
            f"{result['study']:<18} "
            f"{result['case']:<36} "
            f"{result['mapper']:<38} "
            f"{result['characteristic_h']:>10.3e} "
            f"{result['relative_l2_error']:>12.4e} "
            f"{result['rmse_displacement_error']:>12.4e} "
            f"{result['max_displacement_error']:>12.4e}"
        )


def main():
    if OUTPUT_ROOT.exists():
        shutil.rmtree(OUTPUT_ROOT)
    OUTPUT_ROOT.mkdir(parents=True)
    max_curl_difference = validate_rotation_consistency()
    results = run_benchmark()
    csv_path, json_path = write_results(results)
    create_plots(results)
    validate_outputs(results, csv_path, json_path)
    print_summary(results)
    print("")
    print(f"CSV: {csv_path}")
    print(f"JSON: {json_path}")
    print(f"Plots: {PLOT_ROOT}")
    print(f"VTK: {VTK_ROOT}")
    print(f"Max theta - 0.5 curl(g) check: {max_curl_difference:.3e}")


if __name__ == "__main__":
    main()
