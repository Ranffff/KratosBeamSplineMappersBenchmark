#!/usr/bin/env python3
"""Reproduce the straight-beam recovery matrix and quantify its conditioning."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MDPA = ROOT / "CoSimulation_Cases" / "BeamSplineMapper_WithRotationalRecovery" / "Mok_CSM.mdpa"


OPERATORS = (
    ((0, -1, 1.0),),
    ((1, -1, 1.0),),
    ((2, -1, 1.0),),
    ((1, 2, -1.0), (2, 1, 1.0)),
    ((0, 2, 1.0), (2, 0, -1.0)),
    ((0, 1, -1.0), (1, 0, 1.0)),
)


def read_mdpa_nodes(path: Path) -> np.ndarray:
    coordinates = []
    inside = False
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped == "Begin Nodes":
            inside = True
            continue
        if stripped == "End Nodes":
            break
        if inside and stripped and not stripped.startswith("//"):
            fields = stripped.split()
            coordinates.append([float(fields[1]), float(fields[2]), float(fields[3])])
    if len(coordinates) < 2:
        raise RuntimeError(f"Could not read a beam chain from {path}")
    return np.asarray(coordinates, dtype=float)


def gaussian_data(difference: np.ndarray, radius: float, corrected_origin: bool):
    rho = np.linalg.norm(difference)
    h2 = radius * radius
    phi = np.exp(-0.5 * np.dot(difference, difference) / h2)
    gradient = -(difference / h2) * phi
    if rho < np.finfo(float).eps:
        hessian = (-np.eye(3) / h2) if corrected_origin else np.zeros((3, 3))
    else:
        hessian = (np.outer(difference, difference) / (h2 * h2) - np.eye(3) / h2) * phi
    return phi, gradient, hessian


def functional_block(row: np.ndarray, column: np.ndarray, radius: float, corrected_origin: bool) -> np.ndarray:
    phi, gradient, hessian = gaussian_data(row - column, radius, corrected_origin)
    block = np.zeros((6, 6))
    for i, row_operators in enumerate(OPERATORS):
        for j, column_operators in enumerate(OPERATORS):
            value = 0.0
            for row_component, row_direction, row_coefficient in row_operators:
                for column_component, column_direction, column_coefficient in column_operators:
                    if row_component != column_component:
                        continue
                    if row_direction < 0 and column_direction < 0:
                        derivative = phi
                    elif row_direction >= 0 and column_direction < 0:
                        derivative = gradient[row_direction]
                    elif row_direction < 0 and column_direction >= 0:
                        derivative = -gradient[column_direction]
                    else:
                        derivative = -hessian[row_direction, column_direction]
                    value += row_coefficient * column_coefficient * derivative
            block[i, j] = value
    return block


def polynomial_modes(point: np.ndarray, reference: np.ndarray, tangent: np.ndarray):
    values = np.zeros((30, 3))
    curls = np.zeros((30, 3))
    relative = point - reference
    s = float(np.dot(relative, tangent))
    axes = np.eye(3)

    values[0:3] = axes
    values[3] = np.array([0.0, -relative[2], relative[1]])
    values[4] = np.array([relative[2], 0.0, -relative[0]])
    values[5] = np.array([-relative[1], relative[0], 0.0])
    curls[3:6] = 2.0 * axes
    for axis_index, axis in enumerate(axes):
        values[6 + axis_index] = s * axis
        curls[6 + axis_index] = np.cross(tangent, axis)
        rigid_rotation = np.cross(axis, relative)
        values[9 + axis_index] = s * rigid_rotation
        curls[9 + axis_index] = np.cross(tangent, rigid_rotation) + 2.0 * s * axis

    mode = 12
    for degree in (2, 3):
        for axis in axes:
            values[mode] = s**degree * axis
            curls[mode] = degree * s ** (degree - 1) * np.cross(tangent, axis)
            mode += 1
    for degree in (2, 3, 4, 5):
        for axis in axes:
            rigid_rotation = np.cross(axis, relative)
            values[mode] = s**degree * rigid_rotation
            curls[mode] = (
                degree * s ** (degree - 1) * np.cross(tangent, rigid_rotation)
                + 2.0 * s**degree * axis
            )
            mode += 1
    assert mode == 30
    return values, curls


def interpolation_index(node: int, functional: int, number_of_nodes: int) -> int:
    if functional < 3:
        return 3 * node + functional
    return 3 * number_of_nodes + 3 * node + functional - 3


def build_system(coordinates: np.ndarray, radius: float, regularization: float, corrected_origin: bool):
    number_of_nodes = coordinates.shape[0]
    interpolation_size = 6 * number_of_nodes
    system = np.zeros((interpolation_size + 30, interpolation_size + 30))
    for i, row in enumerate(coordinates):
        for j, column in enumerate(coordinates):
            block = functional_block(row, column, radius, corrected_origin)
            for local_row in range(6):
                global_row = interpolation_index(i, local_row, number_of_nodes)
                for local_column in range(6):
                    global_column = interpolation_index(j, local_column, number_of_nodes)
                    system[global_row, global_column] = block[local_row, local_column]

    reference = coordinates[0]
    tangent = coordinates[-1] - coordinates[0]
    tangent /= np.linalg.norm(tangent)
    polynomial = np.zeros((interpolation_size, 30))
    for i, point in enumerate(coordinates):
        values, curls = polynomial_modes(point, reference, tangent)
        polynomial[3 * i : 3 * i + 3, :] = values.T
        polynomial[3 * number_of_nodes + 3 * i : 3 * number_of_nodes + 3 * i + 3, :] = curls.T
    system[:interpolation_size, interpolation_size:] = polynomial
    system[interpolation_size:, :interpolation_size] = polynomial.T
    system[np.arange(interpolation_size), np.arange(interpolation_size)] += regularization
    return system


def summarize(matrix: np.ndarray, number_of_nodes: int) -> dict:
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    threshold = np.finfo(float).eps * max(matrix.shape) * singular_values[0]
    row_scale = np.maximum(np.max(np.abs(matrix), axis=1), np.finfo(float).tiny)
    scale = 1.0 / np.sqrt(row_scale)
    equilibrated = scale[:, None] * matrix * scale[None, :]
    equilibrated_singular_values = np.linalg.svd(equilibrated, compute_uv=False)
    rng = np.random.default_rng(20260722)
    rhs = rng.standard_normal(matrix.shape[0])
    solution = np.linalg.solve(matrix, rhs)
    relative_residual = np.linalg.norm(matrix @ solution - rhs) / np.linalg.norm(rhs)
    curl_start = 3 * number_of_nodes
    return {
        "size": list(matrix.shape),
        "symmetry_relative_frobenius": float(np.linalg.norm(matrix - matrix.T) / np.linalg.norm(matrix)),
        "largest_singular_value": float(singular_values[0]),
        "smallest_singular_value": float(singular_values[-1]),
        "condition_number_2": float(singular_values[0] / singular_values[-1]),
        "numerical_rank": int(np.count_nonzero(singular_values > threshold)),
        "rank_threshold": float(threshold),
        "equilibrated_condition_number_2": float(equilibrated_singular_values[0] / equilibrated_singular_values[-1]),
        "random_solve_relative_residual": float(relative_residual),
        "first_node_curl_diagonal": [float(matrix[curl_start + component, curl_start + component]) for component in range(3)],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mdpa", type=Path, default=DEFAULT_MDPA)
    parser.add_argument("--kernel-radius", type=float, default=0.5)
    parser.add_argument("--regularization", type=float, default=1.0e-8)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "TestCase_Output" / "Diagnostics_20260722" / "RecoveryMatrixConditioning.json",
    )
    args = parser.parse_args()
    coordinates = read_mdpa_nodes(args.mdpa)
    results = {
        "mdpa": str(args.mdpa.resolve()),
        "number_of_support_nodes": int(coordinates.shape[0]),
        "kernel": "gaussian",
        "kernel_radius": args.kernel_radius,
        "regularization": args.regularization,
        "polynomial_level": 4,
        "original_zero_hessian": summarize(
            build_system(coordinates, args.kernel_radius, args.regularization, False),
            coordinates.shape[0],
        ),
        "corrected_analytic_hessian": summarize(
            build_system(coordinates, args.kernel_radius, args.regularization, True),
            coordinates.shape[0],
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite existing diagnostic: {args.output}")
    args.output.write_text(
        json.dumps(results, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(results, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
