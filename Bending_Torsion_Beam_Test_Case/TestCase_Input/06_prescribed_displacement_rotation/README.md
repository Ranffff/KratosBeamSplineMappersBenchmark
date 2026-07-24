# Prescribed Displacement And Rotation Case

This case uses the shared coarse straight beam model:

```text
beam_geometry_coarse.mdpa
```

Boundary conditions:

```text
node 19: fixed ux, uy, uz, rx, ry, rz
node 1 : prescribed ux =  0.10
node 1 : prescribed uy = -0.25
node 1 : prescribed uz =  0.15
node 1 : prescribed rx =  0.03
node 1 : prescribed ry = -0.04
node 1 : prescribed rz =  0.05
```

Run:

```text
python3 Bending_Torsion_Beam_Test_Case/TestCase_Input/06_prescribed_displacement_rotation/run_prescribed_displacement_rotation_comparison.py
```

The script first solves the beam with `BernoulliSolver.py`, then maps the solved
beam displacements with `beam_spline_mapper` and `beam_mapper`
(`use_corotation=false`) to matching centerline destination nodes and reports
the displacement error norms against the Bernoulli solution.
