# FSI Mok Beam-Spline Mapper Benchmark

This directory contains a cleaned, runnable Mok FSI benchmark set for the
Kratos MappingApplication beam-spline mapper work.  It keeps the simulation
inputs, maintained analysis scripts, theoretical notes, and lightweight result
summaries.  Large generated VTK fields, input snapshots, Python caches, LaTeX
temporary files, and long console logs are intentionally not versioned.

## Required Kratos Branches

Use a Kratos build that contains the current mapper implementations:

- `mapping/beam_splines` for `beam_spline_mapper`
- `mapping/beam_spline_with_recovery_of_rotations` for
  `beam_spline_mapper_with_recovery_of_rotations`

The recovery branch contains both mapper types and is the recommended branch
for running all cases in this folder.

## Directory Layout

- `CoSimulation_Cases/`
  - `NearestNeighbor`
  - `BeamMapper_CoRotation`
  - `BeamSplineMapper`
  - `BeamSplineMapper_WithRotationalRecovery`
- `Scripts/`
  - analytical forward/adjoint verification
  - FSI benchmark runner
  - post-processing and diagnosis scripts
- `Notes/`
  - `Beam_Abstraction_Short_Note.tex/.pdf`
  - `BeamSpline_SO3_Verification_Report.tex/.pdf`
- `TestCase_Output/`
  - lightweight retained summaries, CSV files, JSON files, conclusions, and
    selected point histories

## Running A Case Directly

From a Kratos build tree, export the Python path to the Kratos binaries:

```bash
export PYTHONPATH=/path/to/Kratos/bin/Release
```

Then run one case:

```bash
cd FSI_mok_Test_Case/CoSimulation_Cases/BeamSplineMapper_WithRotationalRecovery
python3 MainKratos.py
```

The case writes output paths defined in its `ProjectParameters*.json` files.
For clean benchmark runs, prefer the runner below because it snapshots inputs,
restores modified files, and writes a compact summary.

## Analytical Verification

Run from the `FSI_mok_Test_Case` directory:

```bash
export PYTHONPATH=/path/to/Kratos/bin/Release

python3 Scripts/validate_beam_spline_forward.py
python3 Scripts/validate_beam_spline_adjoint.py
python3 Scripts/validate_recovery_forward.py
python3 Scripts/validate_recovery_adjoint.py
python3 Scripts/validate_mapper_contracts.py
python3 Scripts/validate_crbeam_rotation_compatibility.py
```

The analytical checks use prescribed fields as references.  They do not use
another mapper as the accuracy reference.

## FSI Benchmark Runner

Example finite-recovery gate run:

```bash
python3 Scripts/run_fsi_benchmark.py \
  BeamSplineMapper_WithRotationalRecovery \
  --tag gate_4p5_recovery_auto7_finite \
  --end-time 4.5 \
  --dt 0.05 \
  --alpha 0.03 \
  --iterations 60 \
  --scale 50 \
  --kernel-radius 0.50 \
  --regularization 1e-8 \
  --polynomial-level 0 \
  --rotation-recovery-mode finite
```

Available runner cases are:

- `NearestNeighbor`
- `BeamMapper_CoRotation`
- `BeamSplineMapper`
- `BeamSplineMapper_WithRotationalRecovery`

Recommended gate sequence:

```text
0.5 s -> 4.5 s -> 10 s -> 25 s
```

Do not promote a mapper to a longer gate until the shorter gate is understood.

## Retained Results

The repository keeps compact evidence only:

- analytical `summary.txt`, `summary.csv`, `summary.json`, and
  `parameters.json`
- FSI `*_summary.json`
- FSI `conclusion.txt`
- aggregate `nn_rel_l2_comparison.csv/.json`
- selected point-history `.dat` files

The repository intentionally excludes:

- `vtk_output_*`
- `input_snapshot`
- `benchmark_logs`
- `console_log.txt`
- `__pycache__`
- LaTeX auxiliary files such as `.aux`, `.log`, `.out`, `.fls`, and
  `.fdb_latexmk`

Re-run the scripts if full VTK fields or detailed logs are needed.

## Current Interpretation

The analytical mapper checks are the main correctness evidence.  The Mok FSI
results are nonlinear coupled-trajectory evidence and should be interpreted
together with the diagnostic summaries.  The rotational-recovery finite mode
is an Ahrem-style rotational recovery plus a rigid-section finite-rotation
hybrid; it is not a finite-rotation theorem from the original small-rotation
curl-recovery formulation.
