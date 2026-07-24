"""Prototype benchmark for enriched line-adapted rotational modes.

This script does not call or modify the Kratos beam spline mapper.  It is a
standalone Hermite RBF prototype used to test whether an enriched polynomial
basis can reproduce the basic beam kinematic modes before touching the C++
mapper implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import json
import math
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/kratos_enriched_line_adapted_mpl")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
OUTPUT_ROOT = ROOT / "TestCase_Output" / "Enriched_Line_Adapted_Prototype"
PLOT_ROOT = OUTPUT_ROOT / "Plots"
GRID_COLOR = "#D9D9D9"

LENGTH = 1.0
WIDTH = 0.20
TRANSLATION_AMPLITUDE = 1.0e-3
ROTATION_AMPLITUDE = 1.0e-2
GAUSSIAN_RADIUS = 0.25
REGULARIZATION = 1.0e-10

BEAM_REFINEMENT_ELEMENTS = (4, 8, 16, 32, 64)
SURFACE_DIVISIONS = {"xi": 16, "eta": 80}


@dataclass(frozen=True)
class BasisDefinition:
    key: str
    label: str
    color: str
    marker: str
    linestyle: str
    enriched: bool


BASES = (
    BasisDefinition(
        key="line_adapted",
        label="line-adapted",
        color="#0072B2",
        marker="o",
        linestyle="--",
        enriched=False,
    ),
    BasisDefinition(
        key="enriched_line_adapted",
        label="enriched line-adapted",
        color="#009E73",
        marker="s",
        linestyle="-", 
        enriched=True,
    ),
)


FUNCTIONALS = (
    ((0, -1, 1.0),),
    ((1, -1, 1.0),),
    ((2, -1, 1.0),),
    ((1, 2, -1.0), (2, 1, 1.0)),
    ((0, 2, 1.0), (2, 0, -1.0)),
    ((0, 1, -1.0), (1, 0, 1.0)),
)


def gaussian_kernel(difference: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    h2 = GAUSSIAN_RADIUS * GAUSSIAN_RADIUS
    rho2 = float(np.dot(difference, difference))
    phi = math.exp(-0.5 * rho2 / h2)
    gradient = -difference / h2 * phi
    hessian = (np.outer(difference, difference) / (h2 * h2) - np.eye(3) / h2) * phi
    return phi, gradient, hessian


def apply_functionals(row_functional: int, column_functional: int, x_row: np.ndarray, x_col: np.ndarray) -> float:
    phi, gradient, hessian = gaussian_kernel(x_row - x_col)
    value = 0.0
    for row_component, row_derivative, row_coefficient in FUNCTIONALS[row_functional]:
        for col_component, col_derivative, col_coefficient in FUNCTIONALS[column_functional]:
            if row_component != col_component:
                continue
            if row_derivative < 0 and col_derivative < 0:
                derivative_value = phi
            elif row_derivative >= 0 and col_derivative < 0:
                derivative_value = gradient[row_derivative]
            elif row_derivative < 0 and col_derivative >= 0:
                derivative_value = -gradient[col_derivative]
            else:
                derivative_value = -hessian[row_derivative, col_derivative]
            value += row_coefficient * col_coefficient * derivative_value
    return value


def cross_matrix(axis: np.ndarray, coordinates: np.ndarray) -> np.ndarray:
    return np.cross(axis, coordinates)


def line_adapted_modes(coordinates: np.ndarray, enriched: bool) -> tuple[list[np.ndarray], list[np.ndarray]]:
    x, y, z = coordinates
    s = y
    ex = np.array([1.0, 0.0, 0.0])
    ey = np.array([0.0, 1.0, 0.0])
    ez = np.array([0.0, 0.0, 1.0])
    t = ey

    values: list[np.ndarray] = []
    curls: list[np.ndarray] = []

    for axis in (ex, ey, ez):
        values.append(axis.copy())
        curls.append(np.zeros(3))

    for axis in (ex, ey, ez):
        values.append(cross_matrix(axis, coordinates))
        curls.append(2.0 * axis)

    for axis in (ex, ey, ez):
        values.append(s * axis)
        curls.append(np.cross(t, axis))

    if enriched:
        for axis in (ex, ey, ez):
            rigid_rotation = cross_matrix(axis, coordinates)
            values.append(s * rigid_rotation)
            curls.append(np.cross(t, rigid_rotation) + 2.0 * s * axis)

    return values, curls


def beam_nodes(number_of_elements: int) -> np.ndarray:
    return np.array([[0.0, LENGTH * i / number_of_elements, 0.0] for i in range(number_of_elements + 1)])


def surface_nodes(xi_divisions: int, eta_divisions: int) -> np.ndarray:
    nodes = []
    for j in range(eta_divisions + 1):
        y = LENGTH * j / eta_divisions
        for i in range(xi_divisions + 1):
            x = -0.5 * WIDTH + WIDTH * i / xi_divisions
            nodes.append([x, y, 0.0])
    return np.array(nodes)


def analytical_state(case_key: str, coordinates: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x, y, z = coordinates
    if case_key == "constant_translation":
        return np.array([0.0, 0.0, TRANSLATION_AMPLITUDE]), np.zeros(3)
    if case_key == "constant_rotation":
        omega = np.array([0.0, ROTATION_AMPLITUDE, 0.0])
        return np.cross(omega, coordinates), omega
    if case_key == "linear_centerline_displacement":
        slope = 2.0 * TRANSLATION_AMPLITUDE
        return np.array([0.0, 0.0, slope * y]), np.array([0.5 * slope, 0.0, 0.0])
    if case_key == "linear_varying_rotation":
        theta = ROTATION_AMPLITUDE * y / LENGTH
        theta_prime = ROTATION_AMPLITUDE / LENGTH
        displacement = np.array([theta * z, 0.0, -theta * x])
        rotation = np.array([-0.5 * theta_prime * x, theta, -0.5 * theta_prime * z])
        return displacement, rotation
    if case_key == "curl_compatible_beam_field":
        w = TRANSLATION_AMPLITUDE * math.sin(math.pi * y / LENGTH)
        w_prime = TRANSLATION_AMPLITUDE * math.pi / LENGTH * math.cos(math.pi * y / LENGTH)
        theta = ROTATION_AMPLITUDE * math.sin(math.pi * y / LENGTH)
        theta_prime = ROTATION_AMPLITUDE * math.pi / LENGTH * math.cos(math.pi * y / LENGTH)
        displacement = np.array([theta * z, 0.0, w - theta * x])
        rotation = np.array([0.5 * (w_prime - theta_prime * x), theta, -0.5 * theta_prime * z])
        return displacement, rotation
    raise RuntimeError(f"Unknown case: {case_key}")


def build_polynomial_matrix(nodes: np.ndarray, enriched: bool) -> np.ndarray:
    number_of_nodes = nodes.shape[0]
    number_of_modes = 12 if enriched else 9
    matrix = np.zeros((6 * number_of_nodes, number_of_modes))
    for i, coordinates in enumerate(nodes):
        values, curls = line_adapted_modes(coordinates, enriched)
        for mode_index, (value, curl) in enumerate(zip(values, curls)):
            matrix[3 * i : 3 * i + 3, mode_index] = value
            matrix[3 * number_of_nodes + 3 * i : 3 * number_of_nodes + 3 * i + 3, mode_index] = curl
    return matrix


def build_system(nodes: np.ndarray, enriched: bool) -> np.ndarray:
    number_of_nodes = nodes.shape[0]
    interpolation_size = 6 * number_of_nodes
    polynomial_matrix = build_polynomial_matrix(nodes, enriched)
    system_size = interpolation_size + polynomial_matrix.shape[1]
    system = np.zeros((system_size, system_size))

    for i, row_coordinates in enumerate(nodes):
        for j, column_coordinates in enumerate(nodes):
            for row_functional in range(6):
                row = 3 * i + row_functional if row_functional < 3 else 3 * number_of_nodes + 3 * i + row_functional - 3
                for column_functional in range(6):
                    column = 3 * j + column_functional if column_functional < 3 else 3 * number_of_nodes + 3 * j + column_functional - 3
                    system[row, column] = apply_functionals(
                        row_functional,
                        column_functional,
                        row_coordinates,
                        column_coordinates,
                    )

    system[:interpolation_size, :interpolation_size] += REGULARIZATION * np.eye(interpolation_size)
    system[:interpolation_size, interpolation_size:] = polynomial_matrix
    system[interpolation_size:, :interpolation_size] = polynomial_matrix.T
    return system


def build_rhs(nodes: np.ndarray, case_key: str, number_of_modes: int) -> np.ndarray:
    number_of_nodes = nodes.shape[0]
    rhs = np.zeros(6 * number_of_nodes + number_of_modes)
    for i, coordinates in enumerate(nodes):
        displacement, rotation = analytical_state(case_key, coordinates)
        rhs[3 * i : 3 * i + 3] = displacement
        rhs[3 * number_of_nodes + 3 * i : 3 * number_of_nodes + 3 * i + 3] = 2.0 * rotation
    return rhs


def solve_coefficients(nodes: np.ndarray, case_key: str, enriched: bool) -> np.ndarray:
    system = build_system(nodes, enriched)
    number_of_modes = 12 if enriched else 9
    rhs = build_rhs(nodes, case_key, number_of_modes)
    try:
        return np.linalg.solve(system, rhs)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(system, rhs, rcond=1.0e-12)[0]


def evaluation_row(nodes: np.ndarray, coordinates: np.ndarray, enriched: bool) -> np.ndarray:
    number_of_nodes = nodes.shape[0]
    polynomial_values, _ = line_adapted_modes(coordinates, enriched)
    row = np.zeros((3, 6 * number_of_nodes + len(polynomial_values)))
    for component in range(3):
        for j, node_coordinates in enumerate(nodes):
            for column_functional in range(6):
                column = 3 * j + column_functional if column_functional < 3 else 3 * number_of_nodes + 3 * j + column_functional - 3
                row[component, column] = apply_functionals(component, column_functional, coordinates, node_coordinates)
        for mode_index, value in enumerate(polynomial_values):
            row[component, 6 * number_of_nodes + mode_index] = value[component]
    return row


def evaluate_displacements(nodes: np.ndarray, evaluation_nodes: np.ndarray, coefficients: np.ndarray, enriched: bool) -> np.ndarray:
    mapped = np.zeros_like(evaluation_nodes)
    for i, coordinates in enumerate(evaluation_nodes):
        mapped[i, :] = evaluation_row(nodes, coordinates, enriched) @ coefficients
    return mapped


def compute_metrics(mapped: np.ndarray, exact: np.ndarray) -> tuple[float, float, float]:
    differences = mapped - exact
    error_norms = np.linalg.norm(differences, axis=1)
    exact_norm_squared = float(np.sum(exact * exact))
    rmse = float(math.sqrt(np.mean(error_norms * error_norms)))
    max_error = float(np.max(error_norms))
    relative_l2 = float(math.sqrt(np.sum(error_norms * error_norms) / exact_norm_squared)) if exact_norm_squared > 0.0 else 0.0
    return rmse, max_error, relative_l2


def run_case(case_key: str, basis: BasisDefinition, number_of_elements: int, evaluation_nodes: np.ndarray) -> dict[str, object]:
    support_nodes = beam_nodes(number_of_elements)
    coefficients = solve_coefficients(support_nodes, case_key, basis.enriched)
    mapped = evaluate_displacements(support_nodes, evaluation_nodes, coefficients, basis.enriched)
    exact = np.array([analytical_state(case_key, coordinates)[0] for coordinates in evaluation_nodes])
    rmse, max_error, relative_l2 = compute_metrics(mapped, exact)
    return {
        "case": case_key,
        "basis": basis.key,
        "basis_label": basis.label,
        "kernel": "gaussian",
        "kernel_radius": GAUSSIAN_RADIUS,
        "regularization": REGULARIZATION,
        "beam_elements": number_of_elements,
        "characteristic_h": LENGTH / number_of_elements,
        "evaluation_nodes": int(evaluation_nodes.shape[0]),
        "rmse": rmse,
        "max_displacement_error": max_error,
        "relative_l2_error": relative_l2,
    }


def run_sanity_checks() -> list[dict[str, object]]:
    cases = (
        "constant_translation",
        "constant_rotation",
        "linear_centerline_displacement",
        "linear_varying_rotation",
    )
    evaluation_nodes = surface_nodes(8, 40)
    results = []
    for case_key in cases:
        for basis in BASES:
            results.append(run_case(case_key, basis, 16, evaluation_nodes))
    return results


def run_refinement_benchmark() -> list[dict[str, object]]:
    evaluation_nodes = surface_nodes(SURFACE_DIVISIONS["xi"], SURFACE_DIVISIONS["eta"])
    results = []
    for number_of_elements in BEAM_REFINEMENT_ELEMENTS:
        for basis in BASES:
            results.append(run_case("curl_compatible_beam_field", basis, number_of_elements, evaluation_nodes))
    return results


def write_results(results: list[dict[str, object]]) -> tuple[Path, Path]:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_ROOT / "enriched_line_adapted_prototype_results.csv"
    json_path = OUTPUT_ROOT / "enriched_line_adapted_prototype_results.json"
    fieldnames = (
        "case",
        "basis",
        "basis_label",
        "kernel",
        "kernel_radius",
        "regularization",
        "beam_elements",
        "characteristic_h",
        "evaluation_nodes",
        "rmse",
        "max_displacement_error",
        "relative_l2_error",
    )
    with csv_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow({name: result[name] for name in fieldnames})
    with json_path.open("w") as output:
        json.dump(results, output, indent=2)
    return csv_path, json_path


def configure_plot_style() -> None:
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


def save_figure(figure, output_path: Path) -> None:
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    figure.savefig(output_path.with_suffix(".svg"), bbox_inches="tight")


def create_metric_plot(results: list[dict[str, object]], metric_key: str, ylabel: str, basename: str) -> None:
    PLOT_ROOT.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(5.4, 4.0))
    for basis in BASES:
        curve = [r for r in results if r["case"] == "curl_compatible_beam_field" and r["basis"] == basis.key]
        curve = sorted(curve, key=lambda r: float(r["characteristic_h"]), reverse=True)
        axis.loglog(
            [r["characteristic_h"] for r in curve],
            [max(float(r[metric_key]), 1.0e-18) for r in curve],
            color=basis.color,
            marker=basis.marker,
            linestyle=basis.linestyle,
            label=basis.label,
        )
    axis.set_xlabel("Characteristic beam size h")
    axis.set_ylabel(ylabel)
    axis.grid(True, which="both", color=GRID_COLOR, linewidth=0.6)
    axis.invert_xaxis()
    axis.legend(frameon=True, framealpha=0.95, loc="best")
    figure.tight_layout()
    save_figure(figure, PLOT_ROOT / f"{basename}.png")
    plt.close(figure)


def create_plots(results: list[dict[str, object]]) -> None:
    configure_plot_style()
    create_metric_plot(results, "rmse", "RMSE", "rmse_curl_compatible_refinement")
    create_metric_plot(
        results,
        "max_displacement_error",
        "Max displacement error",
        "max_displacement_error_curl_compatible_refinement",
    )


def print_summary(results: list[dict[str, object]]) -> None:
    print("\nEnriched line-adapted prototype")
    print(f"kernel = gaussian, radius = {GAUSSIAN_RADIUS:g}, regularization = {REGULARIZATION:g}")
    print(f"{'case':34s} {'basis':24s} {'elems':>5s} {'RMSE':>12s} {'max error':>12s} {'rel L2':>12s}")
    print("-" * 105)
    for result in results:
        if result["case"] != "curl_compatible_beam_field" or result["beam_elements"] in (16, 64):
            print(
                f"{str(result['case']):34s} "
                f"{str(result['basis']):24s} "
                f"{int(result['beam_elements']):5d} "
                f"{float(result['rmse']):12.4e} "
                f"{float(result['max_displacement_error']):12.4e} "
                f"{float(result['relative_l2_error']):12.4e}"
            )


def main() -> None:
    sanity_results = run_sanity_checks()
    refinement_results = run_refinement_benchmark()
    results = sanity_results + refinement_results
    csv_path, json_path = write_results(results)
    create_plots(results)
    print_summary(results)
    print(f"\nCSV: {csv_path}")
    print(f"JSON: {json_path}")
    print(f"Plots: {PLOT_ROOT}")


if __name__ == "__main__":
    main()
