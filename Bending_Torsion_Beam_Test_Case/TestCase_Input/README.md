# Bernoulli Solver Test Case Inputs

This folder stores the mdpa input model used by each validation case.

| Case | Input mdpa | Notes |
|---|---|---|
| 01_cantilever_point_load | `01_cantilever_point_load/beam_geometry_coarse.mdpa` | Shared 19-node, 18-element coarse beam. |
| 02_simply_supported_point_load | `02_simply_supported_point_load/beam_geometry_coarse.mdpa` | Same coarse beam, different command-line supports and load. |
| 03_simply_supported_udl | `03_simply_supported_udl/simply_supported_udl.mdpa` | 11-node, 10-element straight beam reconstructed from the previous C++ UDL VTK result. |
| 04_torsion | `04_torsion/beam_geometry_coarse.mdpa` | Same coarse beam, fixed-free torsion setup. |
| 05_combined_bending_torsion | `05_combined_bending_torsion/beam_geometry_coarse.mdpa` | Same coarse beam, combined end torque and bending moment. |
| 06_prescribed_displacement_rotation | `06_prescribed_displacement_rotation/beam_geometry_coarse.mdpa` | Same coarse beam, one end fixed and the opposite end prescribed in all displacement and rotation DOFs. Includes a mapper comparison script. |

Material properties are read from:

```text
Bending_Torsion_Beam_Test_Case/Bending_Torsion_Beam_Test_Case/beam_geometry/StructuralMaterials.json
```

Loads and boundary conditions are supplied by the solver command-line options for the validation runs.
