#!/usr/bin/env python3
"""Quantify CrBeamElement3D2N versus total rotation-vector semantics."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
from pathlib import Path

import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA


ROOT = Path(__file__).resolve().parent.parent
OUTPUT_ROOT = ROOT / "TestCase_Output" / "MapperVerification"
COMMUTING_TOLERANCE = 1.0e-10


def skew(theta):
    x, y, z = theta
    return [[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]]


def matmul(a, b):
    return [
        [sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)]
        for i in range(3)
    ]


def rotation_matrix(theta):
    angle = math.sqrt(sum(value * value for value in theta))
    k = skew(theta)
    k2 = matmul(k, k)
    if angle < 1.0e-8:
        a = 1.0 - angle * angle / 6.0
        b = 0.5 - angle * angle / 24.0
    else:
        a = math.sin(angle) / angle
        b = (1.0 - math.cos(angle)) / (angle * angle)
    return [
        [
            (1.0 if i == j else 0.0) + a * k[i][j] + b * k2[i][j]
            for j in range(3)
        ]
        for i in range(3)
    ]


def rel_l2(actual, reference):
    actual_flat = [value for row in actual for value in row]
    reference_flat = [value for row in reference for value in row]
    numerator = math.sqrt(
        sum((a - b) ** 2 for a, b in zip(actual_flat, reference_flat))
    )
    denominator = max(
        math.sqrt(sum(value * value for value in reference_flat)), 1.0e-30
    )
    return numerator / denominator


def create_element():
    model = KM.Model()
    model_part = model.CreateModelPart("crbeam")
    model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
    model_part.ProcessInfo[KM.IS_RESTARTED] = False
    for variable in (
        KM.DISPLACEMENT,
        KM.REACTION,
        KM.ROTATION,
        KM.REACTION_MOMENT,
        SMA.POINT_LOAD,
        SMA.POINT_MOMENT,
        KM.VOLUME_ACCELERATION,
        KM.VELOCITY,
        KM.ANGULAR_VELOCITY,
        KM.ACCELERATION,
        KM.ANGULAR_ACCELERATION,
    ):
        model_part.AddNodalSolutionStepVariable(variable)
    for variable, reaction in (
        (KM.DISPLACEMENT_X, KM.REACTION_X),
        (KM.DISPLACEMENT_Y, KM.REACTION_Y),
        (KM.DISPLACEMENT_Z, KM.REACTION_Z),
        (KM.ROTATION_X, KM.REACTION_MOMENT_X),
        (KM.ROTATION_Y, KM.REACTION_MOMENT_Y),
        (KM.ROTATION_Z, KM.REACTION_MOMENT_Z),
    ):
        KM.VariableUtils().AddDof(variable, reaction, model_part)

    model_part.SetBufferSize(2)
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
    properties = model_part.CreateNewProperties(1)
    properties.SetValue(KM.YOUNG_MODULUS, 210.0e9)
    properties.SetValue(KM.DENSITY, 7850.0)
    properties.SetValue(SMA.CROSS_AREA, 0.01)
    properties.SetValue(KM.POISSON_RATIO, 0.30)
    properties.SetValue(SMA.TORSIONAL_INERTIA, 1.0e-5)
    properties.SetValue(SMA.I22, 1.0e-5)
    properties.SetValue(SMA.I33, 1.0e-5)
    properties.SetValue(KM.VOLUME_ACCELERATION, [0.0, 0.0, 0.0])
    properties.SetValue(KM.CONSTITUTIVE_LAW, SMA.LinearElastic3DLaw())
    element = model_part.CreateNewElement(
        "CrBeamElement3D2N", 1, [1, 2], properties
    )
    element.SetValue(KM.LOCAL_AXIS_2, [0.0, 1.0, 0.0])
    model_part.CloneTimeStep(0.0)
    element.Initialize(model_part.ProcessInfo)
    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    scheme.Initialize(model_part)
    return model_part, element, scheme


def set_total_state(model_part, theta):
    rotation = rotation_matrix(theta)
    model_part.Nodes[1].SetSolutionStepValue(KM.DISPLACEMENT, [0.0, 0.0, 0.0])
    model_part.Nodes[2].SetSolutionStepValue(
        KM.DISPLACEMENT,
        [rotation[0][0] - 1.0, rotation[1][0], rotation[2][0]],
    )
    for node in model_part.Nodes:
        node.SetSolutionStepValue(KM.ROTATION, list(theta))


def element_orientation(element, process_info):
    axes = [
        element.CalculateOnIntegrationPoints(variable, process_info)[1]
        for variable in (KM.LOCAL_AXIS_1, KM.LOCAL_AXIS_2, KM.LOCAL_AXIS_3)
    ]
    return [[float(axes[j][i]) for j in range(3)] for i in range(3)]


def commit_increment(model_part, element, scheme, time, theta):
    model_part.CloneTimeStep(time)
    set_total_state(model_part, theta)
    matrix = KM.CompressedMatrix()
    increment = KM.Vector(0)
    residual = KM.Vector(0)
    scheme.InitializeNonLinIteration(model_part, matrix, increment, residual)
    orientation = element_orientation(element, model_part.ProcessInfo)
    scheme.FinalizeNonLinIteration(model_part, matrix, increment, residual)
    return orientation


def run_commuting_case():
    model_part, element, scheme = create_element()
    commit_increment(model_part, element, scheme, 1.0, [1.0e-4, 0.0, 0.0])
    theta = [2.0e-4, 0.0, 0.0]
    actual = commit_increment(model_part, element, scheme, 2.0, theta)
    expected = rotation_matrix(theta)
    error = rel_l2(actual, expected)
    if error > COMMUTING_TOLERANCE:
        raise AssertionError(f"commuting CrBeam orientation rel_l2={error}")
    return {
        "case": "single_axis_commuting",
        "status": "PASS",
        "rel_l2": error,
        "tolerance": COMMUTING_TOLERANCE,
    }


def run_noncommuting_case():
    model_part, element, scheme = create_element()
    commit_increment(model_part, element, scheme, 1.0, [0.30, 0.0, 0.0])
    theta = [0.30, 0.40, 0.0]
    actual = commit_increment(model_part, element, scheme, 2.0, theta)
    expected = rotation_matrix(theta)
    return {
        "case": "noncommuting_increment_diagnostic",
        "status": "DIAGNOSTIC",
        "rel_l2": rel_l2(actual, expected),
        "acceptance_gate": None,
        "interpretation": (
            "CrBeam internally composes incremental quaternions; exp(total nodal "
            "ROTATION) is not the same orientation for non-commuting increments."
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tag",
        default="CrBeam_Rotation_Compatibility_"
        + dt.datetime.now().strftime("%Y%m%dT%H%M%S"),
    )
    args = parser.parse_args()
    output = (OUTPUT_ROOT / args.tag).resolve()
    if output.parent != OUTPUT_ROOT.resolve():
        raise RuntimeError(f"Unsafe output tag: {args.tag}")
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")

    summary = {
        "status": "PASS",
        "error_metric": "rel_l2",
        "mapper_policy": (
            "Use nodal ROTATION as total exponential coordinates; accept "
            "single-axis/commuting compatibility and report general "
            "non-commuting CrBeam motion as a formulation boundary."
        ),
        "results": [run_commuting_case(), run_noncommuting_case()],
    }
    output.mkdir(parents=True)
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"summary: {output / 'summary.json'}")


if __name__ == "__main__":
    main()
