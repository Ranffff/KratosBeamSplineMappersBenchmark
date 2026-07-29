#!/usr/bin/env python3
"""Contract and tangent checks for the two beam-spline mappers.

The script uses prescribed nodal generalized coordinates only.  Finite
differences verify the production analytical tangent; they are never used by
the mapper implementation.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
from pathlib import Path

import KratosMultiphysics as KM
import KratosMultiphysics.LinearSolversApplication  # noqa: F401
import KratosMultiphysics.MappingApplication  # noqa: F401
import KratosMultiphysics.StructuralMechanicsApplication  # noqa: F401


ROOT = Path(__file__).resolve().parent.parent
OUTPUT_ROOT = ROOT / "TestCase_Output" / "MapperVerification"
EPS = 1.0e-6
TANGENT_TOLERANCE = 1.0e-7
CONTRACT_TOLERANCE = 1.0e-12


def rel_l2(actual, reference, floor=1.0e-30):
    numerator = math.sqrt(sum((a - b) ** 2 for a, b in zip(actual, reference)))
    denominator = max(math.sqrt(sum(b * b for b in reference)), floor)
    return numerator / denominator


def make_model(number_of_beam_nodes=6):
    model = KM.Model()
    origin = model.CreateModelPart("beam")
    destination = model.CreateModelPart("surface")
    for model_part in (origin, destination):
        model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.ROTATION)
        model_part.AddNodalSolutionStepVariable(KM.FORCE)
        model_part.AddNodalSolutionStepVariable(KM.MOMENT)
        model_part.AddNodalSolutionStepVariable(KM.REACTION)

    properties = origin.CreateNewProperties(1)
    for i in range(number_of_beam_nodes):
        origin.CreateNewNode(i + 1, float(i), 0.0, 0.0)
    for i in range(number_of_beam_nodes - 1):
        origin.CreateNewElement(
            "CrBeamElement3D2N", i + 1, [i + 1, i + 2], properties
        )

    destination.CreateNewNode(101, 1.25, 0.35, -0.20)
    destination.CreateNewNode(102, number_of_beam_nodes - 2.25, -0.25, 0.30)
    set_state(origin)
    return model, origin, destination


def set_state(origin):
    length = max(node.X0 for node in origin.Nodes)
    for node in origin.Nodes:
        xi = node.X0 / max(length, 1.0)
        node.SetSolutionStepValue(
            KM.DISPLACEMENT,
            [0.03 * xi, 0.18 * xi * xi, -0.11 * xi * (1.0 - 0.2 * xi)],
        )
        node.SetSolutionStepValue(
            KM.ROTATION,
            [0.22 * xi, -0.17 * xi * (1.0 - xi), 0.31 * xi * xi],
        )


def perturbation(origin):
    length = max(node.X0 for node in origin.Nodes)
    values = {}
    for node in origin.Nodes:
        xi = node.X0 / max(length, 1.0)
        values[node.Id] = (
            (0.013 * (1.0 + xi), -0.027 * xi * xi, 0.019 * xi),
            (0.021 * xi, -0.016 * xi * (1.0 - xi), 0.018 * xi * xi),
        )
    return values


def mapper_settings(mapper_type, mode=None, level=4, legacy_basis=None):
    settings = {
        "mapper_type": mapper_type,
        "search_settings": {
            "search_radius": 2.0,
            "max_num_search_iterations": 30,
        },
        "local_coord_tolerance": 0.25,
        "echo_level": 0,
    }
    if "with_recovery" in mapper_type:
        settings.update(
            {
                "kernel_type": "gaussian",
                "kernel_radius": 0.50,
                "regularization": 1.0e-8,
            }
        )
        if legacy_basis is None:
            settings["polynomial_level"] = level
        else:
            settings["polynomial_basis"] = legacy_basis
        if mode is not None:
            settings["rotation_recovery_mode"] = mode
    return KM.Parameters(json.dumps(settings))


def make_mapper(origin, destination, mapper_type, mode=None, level=4, legacy_basis=None):
    return KM.MapperFactory.CreateMapper(
        origin,
        destination,
        mapper_settings(mapper_type, mode, level, legacy_basis),
    )


def expect_error(label, callback, expected_fragment):
    try:
        callback()
    except RuntimeError as error:
        message = str(error)
        if expected_fragment not in message:
            raise AssertionError(
                f"{label}: expected {expected_fragment!r} in {message!r}"
            ) from error
        return {"name": label, "status": "PASS", "message": message.splitlines()[0]}
    raise AssertionError(f"{label}: expected RuntimeError")


def capture_surface(destination):
    return [
        float(value)
        for node in destination.Nodes
        for value in node.GetSolutionStepValue(KM.DISPLACEMENT)
    ]


def capture_generalized_load(origin):
    result = []
    for node in origin.Nodes:
        result.extend(float(v) for v in node.GetSolutionStepValue(KM.FORCE))
        result.extend(float(v) for v in node.GetSolutionStepValue(KM.MOMENT))
    return result


def apply_perturbation(origin, direction, scale):
    for node in origin.Nodes:
        displacement = list(node.GetSolutionStepValue(KM.DISPLACEMENT))
        rotation = list(node.GetSolutionStepValue(KM.ROTATION))
        delta_displacement, delta_rotation = direction[node.Id]
        for i in range(3):
            displacement[i] += scale * delta_displacement[i]
            rotation[i] += scale * delta_rotation[i]
        node.SetSolutionStepValue(KM.DISPLACEMENT, displacement)
        node.SetSolutionStepValue(KM.ROTATION, rotation)


def flatten_direction(origin, direction):
    result = []
    for node in origin.Nodes:
        result.extend(direction[node.Id][0])
        result.extend(direction[node.Id][1])
    return result


def set_surface_basis_force(destination, basis_index):
    for node_index, node in enumerate(destination.Nodes):
        force = [0.0, 0.0, 0.0]
        local_index = basis_index - 3 * node_index
        if 0 <= local_index < 3:
            force[local_index] = 1.0
        node.SetSolutionStepValue(KM.REACTION, force)


def analytical_directional_tangent(mapper, origin, destination, direction):
    generalized_direction = flatten_direction(origin, direction)
    surface_dofs = 3 * destination.NumberOfNodes()
    result = []
    for basis_index in range(surface_dofs):
        set_surface_basis_force(destination, basis_index)
        mapper.InverseMap(KM.FORCE, KM.MOMENT, KM.REACTION, KM.Flags())
        transpose_row = capture_generalized_load(origin)
        result.append(sum(a * b for a, b in zip(transpose_row, generalized_direction)))
    return result


def finite_difference_directional_tangent(mapper, origin, destination, direction):
    apply_perturbation(origin, direction, EPS)
    mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
    plus = capture_surface(destination)
    apply_perturbation(origin, direction, -2.0 * EPS)
    mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
    minus = capture_surface(destination)
    apply_perturbation(origin, direction, EPS)
    return [(a - b) / (2.0 * EPS) for a, b in zip(plus, minus)]


def coordinate_snapshot(origin, destination):
    return [
        (node.Id, node.X0, node.Y0, node.Z0, node.X, node.Y, node.Z)
        for model_part in (origin, destination)
        for node in model_part.Nodes
    ]


def run_contract_checks():
    checks = []

    _, short_origin, short_destination = make_model(2)
    checks.append(
        expect_error(
            "explicit_level_4_is_strict",
            lambda: make_mapper(
                short_origin,
                short_destination,
                "beam_spline_mapper_with_recovery_of_rotations",
                "finite",
                4,
            ),
            "requires at least 5 beam support nodes",
        )
    )
    auto_mapper = make_mapper(
        short_origin,
        short_destination,
        "beam_spline_mapper_with_recovery_of_rotations",
        "finite",
        0,
    )
    auto_mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
    checks.append({"name": "level_0_auto_short_chain", "status": "PASS"})

    _, origin, destination = make_model()
    legacy_mapper = make_mapper(
        origin,
        destination,
        "beam_spline_mapper_with_recovery_of_rotations",
        mode="small",
        legacy_basis="line_adapted",
    )
    legacy_mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
    checks.append({"name": "legacy_polynomial_basis", "status": "PASS"})

    default_small = make_mapper(
        origin,
        destination,
        "beam_spline_mapper_with_recovery_of_rotations",
        mode=None,
        level=4,
    )
    default_small.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
    default_values = capture_surface(destination)
    explicit_small = make_mapper(
        origin,
        destination,
        "beam_spline_mapper_with_recovery_of_rotations",
        mode="small",
        level=4,
    )
    explicit_small.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
    small_error = rel_l2(capture_surface(destination), default_values)
    if small_error > CONTRACT_TOLERANCE:
        raise AssertionError(f"default small-mode rel_l2={small_error}")
    checks.append(
        {"name": "default_mode_is_small", "status": "PASS", "rel_l2": small_error}
    )

    for mapper_type, mode in (
        ("beam_spline_mapper", None),
        ("beam_spline_mapper_with_recovery_of_rotations", "finite"),
    ):
        _, zero_origin, zero_destination = make_model()
        for node in zero_origin.Nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT, [0.0, 0.0, 0.0])
            node.SetSolutionStepValue(KM.ROTATION, [0.0, 0.0, 0.0])
        zero_mapper = make_mapper(
            zero_origin, zero_destination, mapper_type, mode
        )
        for node in zero_destination.Nodes:
            node.SetSolutionStepValue(KM.REACTION, [0.2, -0.1, 0.3])
        zero_mapper.InverseMap(
            KM.FORCE, KM.MOMENT, KM.REACTION, KM.Flags()
        )
        checks.append(
            {
                "name": f"{mapper_type}_zero_initial_inverse",
                "status": "PASS",
            }
        )

        _, local_origin, local_destination = make_model()
        mapper = make_mapper(local_origin, local_destination, mapper_type, mode)
        checks.append(
            expect_error(
                f"{mapper_type}_inverse_before_map",
                lambda m=mapper: m.InverseMap(
                    KM.FORCE, KM.MOMENT, KM.REACTION, KM.Flags()
                ),
                "requires a successful preceding Map call",
            )
        )
        mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
        before = coordinate_snapshot(local_origin, local_destination)
        mapper.UpdateInterface(KM.Mapper.REMESHED)
        checks.append(
            expect_error(
                f"{mapper_type}_update_invalidates_state",
                lambda m=mapper: m.InverseMap(
                    KM.FORCE, KM.MOMENT, KM.REACTION, KM.Flags()
                ),
                "requires a successful preceding Map call",
            )
        )
        mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
        after = coordinate_snapshot(local_origin, local_destination)
        if before != after:
            raise AssertionError(f"{mapper_type} modified current/reference coordinates")
        checks.append(
            {"name": f"{mapper_type}_coordinate_immutability", "status": "PASS"}
        )

    return checks


def run_finite_tangent_check():
    _, origin, destination = make_model()
    mapper = make_mapper(
        origin,
        destination,
        "beam_spline_mapper_with_recovery_of_rotations",
        "finite",
        4,
    )
    mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
    direction = perturbation(origin)
    analytical = analytical_directional_tangent(
        mapper, origin, destination, direction
    )
    finite_difference = finite_difference_directional_tangent(
        mapper, origin, destination, direction
    )
    tangent_error = rel_l2(analytical, finite_difference)
    if tangent_error > TANGENT_TOLERANCE:
        raise AssertionError(f"finite directional tangent rel_l2={tangent_error}")

    mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT, KM.Flags())
    for node in destination.Nodes:
        node.SetSolutionStepValue(KM.REACTION, [0.7, -0.3, 0.2])
    mapper.InverseMap(KM.FORCE, KM.MOMENT, KM.REACTION, KM.Flags())
    normal = capture_generalized_load(origin)
    mapper.InverseMap(KM.FORCE, KM.MOMENT, KM.REACTION, KM.Mapper.SWAP_SIGN)
    swapped = capture_generalized_load(origin)
    swap_error = rel_l2(swapped, [-value for value in normal])
    if swap_error > CONTRACT_TOLERANCE:
        raise AssertionError(f"SWAP_SIGN rel_l2={swap_error}")

    return {
        "name": "finite_analytical_directional_tangent",
        "status": "PASS",
        "epsilon": EPS,
        "rel_l2": tangent_error,
        "tolerance": TANGENT_TOLERANCE,
        "swap_sign_rel_l2": swap_error,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tag",
        default="Mapper_Contracts_" + dt.datetime.now().strftime("%Y%m%dT%H%M%S"),
    )
    args = parser.parse_args()
    output = (OUTPUT_ROOT / args.tag).resolve()
    if output.parent != OUTPUT_ROOT.resolve():
        raise RuntimeError(f"Unsafe output tag: {args.tag}")
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")

    results = run_contract_checks()
    results.append(run_finite_tangent_check())
    output.mkdir(parents=True)
    summary = {
        "status": "PASS",
        "error_metric": "rel_l2",
        "results": results,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"summary: {output / 'summary.json'}")


if __name__ == "__main__":
    main()
