# Shape metrics example

`in.shape_metrics` moves and rotates a triangulated 2 × 3 × 4 box. It writes
the SHAPE `minVolume`, level-3 OBB axes and form indices through the existing
`compute rigid/local_abrade` interface. The molecule has eight surface nodes
and a central owning atom; its sparse surface is intended for a geometry
demonstration, not for contact-resolution studies.

Run from this directory using a LAMMPS build with GRANULAR, MOLECULE and RIGID:

```sh
/absolute/path/to/lmp -in in.shape_metrics
```

`shape.dump` should report full axis lengths `4, 3, 2`, `obb_valid = 1`,
volume `24`, and surface area `52` within floating-point tolerance. Shape
indices are dimensionless. The example also writes existing body properties,
atomic coordinates and a single-column axis output for inspection.

The source example contains no `offset`; append `offset` to the
`fix abrade` line to exercise the inward contact-node convention. The new
shape output reconstructs the corresponding physical surface in a copy.

## Isolated integration checks

The standard-library Python harness generates fresh test inputs and writes
all outputs to an empty directory supplied by the caller:

```sh
python3 validate_shape_metrics.py \
  --lmp /absolute/path/to/new/lmp \
  --baseline /absolute/path/to/original/lmp \
  --output-dir /tmp/shape-validation
```

It checks analytic axes and all new scalar formulas, `run 0`, offset on/off,
single-column output, one/two MPI ranks, translation and rotation across
periodic boundaries, consecutive runs, restart, and actual abrasive impacts
with one and two MPI ranks. For each processor count it compares old body
properties and atomic output exactly with shape output enabled/disabled, and
optionally with the original binary. If the impact deletes surface nodes,
the report records successful shape output after that deletion. It rejects
stale output directories and times out failed MPI runs.

Use `--mpiexec /absolute/path/to/mpiexec` to choose an MPI launcher. Add
`--include-remesh` to exercise actual node deletion through initial
`equalise`. If that case fails in both old and new builds with the existing
MPICH integer-datatype error, the report records this as a baseline limitation
rather than claiming full remesh validation. A successful remesh run must
actually remove nodes and retain valid shape output.

`validation.json` records executed cases, return codes, elapsed times and any
known baseline limitations. A known limitation sets the overall `passed`
field to false, even when all regular shape checks pass. The standalone OBB reference tests are maintained
separately from these LAMMPS integration checks.
