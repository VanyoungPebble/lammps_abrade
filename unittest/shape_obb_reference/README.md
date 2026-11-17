# Standalone SHAPE OBB reference checks

This directory verifies the geometry kernel without running LAMMPS. The frozen
fixtures come from an independent Python translation of the local SHAPE
`minboundbox.m` with metric `v` and level `3`. The translation retains the
original Euler-angle reconstruction, face candidates, three cyclic edge frames,
and projected minimum-area rectangle calculation. SciPy supplies Qhull convex
hulls independently of the C++ implementation.

These are **not results executed in MATLAB or Octave**. No such runtime was
available when the fixtures were created. The original source hash is recorded
in `reference_cases.json`; direct MATLAB validation remains a separate check.

From the repository root:

```sh
c++ -std=c++11 -O2 -Isrc/RIGID unittest/shape_obb_reference/kernel_driver.cpp src/RIGID/shape_obb.cpp -o /tmp/shape_obb_driver
python3 unittest/shape_obb_reference/reference_minboundbox.py --kernel /tmp/shape_obb_driver
```

Checking frozen fixtures requires only Python's standard library. Regenerating
them additionally requires NumPy and SciPy:

```sh
python3 unittest/shape_obb_reference/reference_minboundbox.py --regenerate --shape-source /path/to/SHAPE/lib/minboundbox/minboundbox.m
```

The checks cover a box before/after a rigid transform, the tetrahedron in the
original SHAPE example, generic asymmetric clouds, interior and duplicate points,
a sampled ellipsoid, scales of 1e-9 and 1e9, and invalid or degenerate clouds.
The existing 386-surface-node cube molecule is also tested in nine orientations;
its densely sampled, nearly coplanar faces exercise the hull's numerical
handling. The original molecule file is reused rather than copied into fixtures.

The two `coordinate_sensitive_cloud` fixtures deliberately retain a limitation
of SHAPE's level-3 heuristic: its additional perpendicular edge directions are
coordinate dependent, so an arbitrary cloud can have slightly different
candidate minima after a rotation. Comparisons must use the same coordinate
frame. This test suite does not replace that algorithm with a different,
rotation-invariant OBB definition.

The translated MATLAB algorithm retains its upstream BSD notice in
`LICENSE.minboundbox`.
