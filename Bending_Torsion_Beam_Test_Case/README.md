# Beam Splines Test Cases

Personal test and validation workspace for the Kratos
`beam_spline_mapper`. This repository contains beam and surface meshes,
forward- and inverse-mapping checks, mapper comparisons, a small
Euler-Bernoulli beam solver, and previously generated VTK/plot/report output.

This is an independent personal repository. It is not part of the Kratos
source tree, even though it is most convenient to clone it next to the local
Kratos files.

## What This Repository Is For

The main uses are:

- compare `beam_spline_mapper` with the linear and co-rotational beam mappers;
- validate forward structure-to-surface displacement mapping;
- validate inverse surface-to-structure force and moment mapping;
- check the adjoint/work consistency of `Map` and `InverseMap`;
- run small Euler-Bernoulli beam reference problems;
- retain meshes, logs, VTK files, plots, and LaTeX reports from experiments.

The mapper implementation itself is **not** maintained here. The scripts use
the `MappingApplication` compiled in the active Kratos build.

## Important Path Convention

Several scripts currently calculate their paths with:

```python
ROOT = Path(__file__).resolve().parents[1]
```

and then look for:

```text
ROOT/Bending_Torsion_Beam_Test_Case/
```

Therefore, clone this repository with the local directory name
`Bending_Torsion_Beam_Test_Case`. The recommended layout is:

```text
/root/dev/Kratos/
├── applications/
├── bin/
├── build/
└── Bending_Torsion_Beam_Test_Case/   # this repository
```

Recommended clone command:

```bash
cd /root/dev/Kratos
git clone git@github.com:Ranffff/beam_splines_testcases.git \
  Bending_Torsion_Beam_Test_Case
```

If the repository is cloned as `beam_splines_testcases`, either rename that
directory or update the `ROOT`, `CASE_DIR`, and `OUTPUT_ROOT` constants in the
top-level scripts.

## Requirements

The Kratos build used for these tests must provide:

- `KratosMultiphysics`;
- `MappingApplication`;
- `StructuralMechanicsApplication`;
- `LinearSolversApplication` for the inverse-mapping checks;
- `CoSimulationApplication` for the original co-simulation case;
- a compiled and registered `beam_spline_mapper`.

Python-side tools:

- Python 3;
- NumPy for `BernoulliSolver.py`;
- Matplotlib for summary plots.

The commands below assume:

```text
KRATOS_ROOT=/root/dev/Kratos
KRATOS_BUILD=/root/dev/Kratos/bin/Release
```

Set up the shell:

```bash
cd /root/dev/Kratos
export PYTHONPATH="$PWD/bin/Release${PYTHONPATH:+:$PYTHONPATH}"
```

Verify the required applications and mapper:

```bash
python3 -c "
import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.LinearSolversApplication
print('beam_spline_mapper registered:',
      KM.MapperFactory.HasMapper('beam_spline_mapper'))
"
```

The final line must report:

```text
beam_spline_mapper registered: True
```

If the mapper source was changed, rebuild the MappingApplication before
running these cases. From the Kratos root, for the current build layout:

```bash
cmake --build build/applications/MappingApplication
```

## Quick Start

All commands in this README are intended to be run from the Kratos root:

```bash
cd /root/dev/Kratos
export PYTHONPATH="$PWD/bin/Release${PYTHONPATH:+:$PYTHONPATH}"
```

Run one small inverse-map work-consistency case:

```bash
python3 Bending_Torsion_Beam_Test_Case/run_beam_spline_inverse_work_check.py \
  --mesh-preset small \
  --case-limit 1
```

The command prints a short completion line. Detailed output is written to:

```text
Bending_Torsion_Beam_Test_Case/TestCase_Output/InverseMap_Work_Check/
```

Check:

```text
TestCase_Output/InverseMap_Work_Check/summary.txt
```

The expected final status is `PASS`. A nonzero process exit code means at
least one work-consistency case exceeded its tolerance.

## Main Workflows

### 1. Inverse Map Work/Adjoint Check

Script:

```text
run_beam_spline_inverse_work_check.py
```

This is the most direct validation of `beam_spline_mapper::InverseMap`. It:

1. prescribes beam displacement and rotation;
2. maps the beam motion to the surface;
3. applies a surface force field;
4. inverse-maps the force and moment to the beam;
5. perturbs the beam kinematics;
6. compares surface directional work with beam generalized work.

Fast check:

```bash
python3 Bending_Torsion_Beam_Test_Case/run_beam_spline_inverse_work_check.py \
  --mesh-preset small \
  --case-limit 1
```

All four built-in cases:

```bash
python3 Bending_Torsion_Beam_Test_Case/run_beam_spline_inverse_work_check.py \
  --mesh-preset small
```

Available mesh presets:

| Preset | Beam | Surface | Intended use |
|---|---|---|---|
| `small` | coarse | 1x5 | smoke test and debugging |
| `dense_surface` | coarse | 10x100 | surface-resolution check |
| `fine` | fine | 10x100 | expensive full-resolution check |

Custom meshes can override the preset:

```bash
python3 Bending_Torsion_Beam_Test_Case/run_beam_spline_inverse_work_check.py \
  --beam-mdpa /absolute/path/to/beam.mdpa \
  --surface-mdpa /absolute/path/to/surface.mdpa
```

Outputs:

```text
TestCase_Output/InverseMap_Work_Check/
├── summary.txt
└── <case-name>/
    └── console_log.txt
```

### 2. Compare Inverse Mapping Between Mappers

Script:

```text
run_inverse_mapper_comparison.py
```

This runs the same inverse-load idea with:

- `beam_spline_mapper`;
- linear `beam_mapper`;
- co-rotational `beam_mapper`, used as the comparison reference.

Quick run:

```bash
python3 Bending_Torsion_Beam_Test_Case/run_inverse_mapper_comparison.py \
  --mesh-preset small \
  --case-limit 1
```

Full built-in comparison:

```bash
python3 Bending_Torsion_Beam_Test_Case/run_inverse_mapper_comparison.py \
  --mesh-preset small
```

For the denser combinations, replace `small` with `dense_surface` or `fine`.

Results are written to:

```text
TestCase_Output/InverseMap_Mapper_Comparison/
```

The main table is:

```text
TestCase_Output/InverseMap_Mapper_Comparison/summary.txt
```

Useful columns include work error, force/moment differences relative to the
co-rotational mapper, and setup/forward/inverse timings.

### 3. Forward Structure-to-Surface Sweep

Script:

```text
run_structure_to_surface_mapper_comparison.py
```

This compares surface displacements produced by:

- `beam_spline_mapper`;
- linear `beam_mapper`;
- co-rotational `beam_mapper`, used as the reference.

This script does not currently have command-line options. Before running,
edit the constants near the top:

```python
KINEMATICS_MODE = "rotation"
ERROR_TYPE = "normalized_rmse"
```

Supported modes:

| Mode | Values used |
|---|---|
| `rotation` | `TIP_ROTATIONS` |
| `displacement` | `TIP_DISPLACEMENTS` |
| `torsion` | `TIP_TORSIONS` |
| `combined` | `COMBINED_KINEMATICS_CASES` |

Also check these mesh constants near the top:

```python
BEAM_MDPA_PATH = ...
SURFACE_MDPA_PATH = ...
```

Run:

```bash
python3 \
  Bending_Torsion_Beam_Test_Case/run_structure_to_surface_mapper_comparison.py
```

Depending on `KINEMATICS_MODE`, output goes to one of:

```text
TestCase_Output/Rotation_1x5/
TestCase_Output/Displacement/
TestCase_Output/Torsion/
TestCase_Output/Combined/
```

Each parameter value gets a subdirectory with mapper VTK files and
`console_log.txt`. The script appends a new block to `sweep_summary.txt`;
repeated runs therefore preserve previous summary blocks rather than
replacing the file.

### 4. Generate Summary Plots

After running the forward sweeps:

```bash
python3 \
  Bending_Torsion_Beam_Test_Case/TestCase_Output/Summary_Plots/plot_mapper_sweep_summary_png.py
```

To parse individual `console_log.txt` files and generate RMSE CSV/plots for a
specific result directory:

```bash
python3 \
  Bending_Torsion_Beam_Test_Case/TestCase_Output/Summary_Plots/plot_rmse_from_console_logs.py \
  Bending_Torsion_Beam_Test_Case/TestCase_Output/Rotation_1x5
```

Matplotlib uses a local cache under `TestCase_Output/Summary_Plots/.matplotlib`
so the scripts also work in restricted environments without a writable home
cache.

### 5. Run the Python Euler-Bernoulli Solver

`BernoulliSolver.py` is a small 3D Euler-Bernoulli reference solver built with
Kratos model I/O and NumPy linear algebra. It accepts either:

- a Kratos `ProjectParameters.json`; or
- an `.mdpa` model plus command-line loads and boundary conditions.

Show all options:

```bash
python3 Bending_Torsion_Beam_Test_Case/BernoulliSolver.py --help
```

Example: fixed-free beam with an end torque:

```bash
python3 Bending_Torsion_Beam_Test_Case/BernoulliSolver.py \
  Bending_Torsion_Beam_Test_Case/TestCase_Input/04_torsion/beam_geometry_coarse.mdpa \
  --materials \
  Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json \
  --fix-node 19 \
  --load 1 0 0 0 1000 0 0 \
  --output \
  Bending_Torsion_Beam_Test_Case/BernoulliSolver_python_TestCases/04_torsion
```

Load syntax:

```text
--load node_id fx fy fz mx my mz
```

Boundary-condition syntax:

```text
--fix-node node_id
--prescribe node_id dof value
```

Valid prescribed DOFs are `ux`, `uy`, `uz`, `rx`, `ry`, and `rz`. Both
`--load` and `--prescribe` may be repeated.

The full commands for validation cases 01-05 are documented at the top of
`BernoulliSolver.py`.

### 6. Prescribed Displacement/Rotation Comparison

This focused case solves a beam with prescribed six-DOF end motion and then
maps its result:

```bash
python3 \
  Bending_Torsion_Beam_Test_Case/TestCase_Input/06_prescribed_displacement_rotation/run_prescribed_displacement_rotation_comparison.py
```

Case details and prescribed values are documented in:

```text
TestCase_Input/06_prescribed_displacement_rotation/README.md
```

### 7. Original Co-Simulation Case

The original coupled case lives under:

```text
Bending_Torsion_Beam_Test_Case/
```

Run it from that directory because `MainKratos.py` opens
`ProjectParametersCoSim.json` with a relative path:

```bash
cd /root/dev/Kratos/Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case
python3 MainKratos.py
```

The active mapper is configured in `ProjectParametersCoSim.json` under:

```text
solver_settings
  -> data_transfer_operators
  -> mapper
  -> mapper_settings
```

Check `mapper_type`, `use_corotation`, search radius, model-part names, and
solver names before reusing this case.

## Directory Guide

| Path | Purpose |
|---|---|
| `Bending_Torsion_Beam_Test_Case/` | original beam/surface co-simulation input |
| `TestCase_Input/` | reusable MDPA inputs for validation cases |
| `TestCase_Output/` | mapper logs, summaries, VTK files, and plots |
| `Mapper_Comparison_TestCases/` | additional mapper-comparison results |
| `BernoulliSolver.py` | Python Euler-Bernoulli reference solver |
| `BernoulliSolver.cpp` | Kratos-linked C++ beam solver experiment |
| `BernoulliBeamSolver_noKratos.cpp` | standalone C++ reference solver |
| `BernoulliSolver_python_TestCases/` | Python solver output and reports |
| `BernoulliSolver_TestCases/` | C++ solver output and report |
| `Archived_BeamSplineMapper_MathUtilsSolver/` | archived mapper implementation snapshot |
| `Figures/` | explanatory LaTeX/PDF figures |

Many generated outputs are intentionally committed because this repository is
also a record of test results. Expect VTK, GiD, LaTeX intermediate, executable,
cache, and `Zone.Identifier` files in the history.

## How to Change a Test Safely

Before changing a test:

1. confirm the active Kratos build and mapper registration;
2. start with `small` meshes and one case;
3. save or commit the current `summary.txt`/`sweep_summary.txt`;
4. change one of kinematics, force field, mesh, or mapper settings at a time;
5. inspect both the numerical summary and VTK output;
6. only then run `dense_surface` or `fine`.

The main places to edit are:

- `TEST_CASES` for inverse checks;
- `MESH_PRESETS` for reusable beam/surface mesh pairs;
- `prescribed_beam_kinematics(...)` for imposed beam motion;
- `surface_force(...)` for surface loading;
- `create_beam_spline_mapper(...)` or `create_mapper(...)` for mapper settings;
- `TIP_*` and `COMBINED_KINEMATICS_CASES` for forward sweeps.

The beam is assumed to have length `10.0` in the current analytical
kinematics. If a different beam geometry is used, update `BEAM_LENGTH`.

## Troubleshooting

### `ModuleNotFoundError: No module named 'KratosMultiphysics'`

The Kratos Python package is not on `PYTHONPATH`:

```bash
cd /root/dev/Kratos
export PYTHONPATH="$PWD/bin/Release${PYTHONPATH:+:$PYTHONPATH}"
```

### `beam_spline_mapper registered: False`

The active MappingApplication does not contain or register the mapper. Check
the source branch, rebuild MappingApplication, and ensure Python is importing
the intended `bin/Release`.

### Input MDPA file cannot be found

Check the local clone name. It must currently be:

```text
Bending_Torsion_Beam_Test_Case
```

and the commands should be launched from the Kratos root as shown above.

### Mapper projection/search failure

Check:

- beam and surface coordinates;
- `search_radius`;
- `max_num_search_iterations`;
- `local_coord_tolerance`;
- the expected beam and surface submodel-part names.

For debugging, increase mapper `echo_level` and use the `small` mesh.

### A fine inverse case is unexpectedly slow

`fine` uses the fine beam and dense 10x100 surface. First run:

```bash
python3 Bending_Torsion_Beam_Test_Case/run_beam_spline_inverse_work_check.py \
  --mesh-preset small \
  --case-limit 1
```

Then test `dense_surface`, and only then `fine`.

### Plots are empty

Run the corresponding forward sweep first and confirm that
`sweep_summary.txt` or the per-case `console_log.txt` files exist. Also check
that Matplotlib is installed in the same Python environment.

## Saving New Work

This directory is its own Git repository. From anywhere:

```bash
cd /root/dev/Kratos/Bending_Torsion_Beam_Test_Case
git status
git add <files-to-save>
git commit -m "Describe the testcase change"
git push
```

Do not use the parent Kratos repository to commit these test cases. The remote
for this repository should be:

```text
git@github.com:Ranffff/beam_splines_testcases.git
```

## Repository

Private GitHub repository:

```text
https://github.com/Ranffff/beam_splines_testcases
```

