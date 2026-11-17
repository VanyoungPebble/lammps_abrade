#!/usr/bin/env python3
"""Small, isolated integration checks for compute rigid/local_abrade shape fields.

Only Python's standard library is required. All simulation files are written to
--output-dir, never to an existing production simulation directory.
"""

import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
import time


HERE = Path(__file__).resolve().parent
OLD_FIELDS = "mol volume surface_area abraded_volume wear_energy inertiax inertiay inertiaz"
SHAPE_FIELDS = (
    "l_obb i_obb s_obb sphericity_wadell obb_valid "
    "zingg_s_over_i_obb zingg_i_over_l_obb sphericity_krumbein_obb "
    "flatness_potticary_obb elongation_potticary_obb "
    "flatness_kong_obb elongation_kong_obb "
    "flatness_angelidakis_obb elongation_angelidakis_obb"
)


def box_molecule(path):
    """A 2 x 3 x 4 box, eight surface nodes and a central owner node."""
    vertices = [(0, 0, 0), (-1, -1.5, -2), (1, -1.5, -2),
                (1, 1.5, -2), (-1, 1.5, -2), (-1, -1.5, 2),
                (1, -1.5, 2), (1, 1.5, 2), (-1, 1.5, 2)]
    faces = [(2, 4, 3), (2, 5, 4), (6, 7, 8), (6, 8, 9),
             (2, 3, 7), (2, 7, 6), (3, 4, 8), (3, 8, 7),
             (4, 5, 9), (4, 9, 8), (5, 2, 6), (5, 6, 9)]
    sections = {
        "Coords": [f"{i} {x} {y} {z}" for i, (x, y, z) in enumerate(vertices, 1)],
        "Types": [f"{i} 1" for i in range(1, 10)],
        "Molecules": [f"{i} 1" for i in range(1, 10)],
        "Diameters": [f"{i} 0.1" for i in range(1, 10)],
        "Masses": [f"{i} 1" for i in range(1, 10)],
        "Angles": [f"{i} 1 {a} {b} {c}" for i, (a, b, c) in enumerate(faces, 1)],
    }
    text = "Analytic box for shape regression\n\n9 atoms\n12 angles\n"
    for name, rows in sections.items():
        text += "\n" + name + "\n\n" + "\n".join(rows) + "\n"
    path.write_text(text)


def output_commands(shape, every):
    text = f"""compute old all rigid/local_abrade abrade {OLD_FIELDS}
dump old all local {every} old.dump c_old[*]
dump_modify old sort 1 format float %.17g
dump positions all custom {every} atoms.dump id type mol xu yu zu vx vy vz fx fy fz radius
dump_modify positions sort id format float %.17g
"""
    if shape:
        text += f"""compute shape all rigid/local_abrade abrade mol {SHAPE_FIELDS}
dump shape all local {every} shape.dump c_shape[*]
dump_modify shape sort 1 format float %.17g
compute single all rigid/local_abrade abrade l_obb
dump single all local {every} single.dump c_single
dump_modify single format float %.17g
"""
    return text


def box_input(shape, offset=False, moving=False, restart=False):
    start = "read_restart state.10\n" if restart else """region box block -10 10 -10 10 -10 10
create_box 4 box angle/types 1 extra/angle/per/atom 10
"""
    create = "" if restart else """molecule Box data.Box
create_atoms 0 single 0 0 0 mol Box 1729871 rotate 37 1 2 3 units box
"""
    motion = """
variable ux atom -50*y+12000
variable uy atom 50*x
velocity all set v_ux v_uy 0
""" if moving and not restart else ""
    extra = " offset" if offset else ""
    extra += " infile state.10.rigid" if restart else ""
    end = "run 10\n" if restart else (
        "restart 10 state\nrun 10\nrun 10\n" if moving else "run 0\n")
    return f"""# Deliberately small isolated regression input
units cgs
boundary p p p
atom_style hybrid molecular sphere
newton on
{start}angle_style zero
angle_coeff * 75.0
{create}neighbor 0.4 multi
comm_modify mode multi vel yes
pair_style granular
pair_coeff * * hertz 1e8 0.55 tangential linear_history 1e8 1.0 0.39 damping coeff_restitution limit_damping
group owning id 1
group owned subtract all owning
set group owning type 4
neigh_modify collection/type 2 1*3 4
comm_modify cutoff/multi 1 6.0 cutoff/multi 2 6.0
neigh_modify exclude molecule/intra all
neigh_modify exclude group owning owned
neigh_modify exclude molecule/inter owning
{motion}fix abrade all rigid/abrade 1e20 0.25 1.42 molecule normals{extra}
{output_commands(shape, 5 if moving or restart else 1)}timestep 1e-4
thermo 10
thermo_style custom step atoms
{end}"""


def impact_input(shape, equalise=False):
    extra = " equalise" if equalise else ""
    return f"""units cgs
boundary p p p
atom_style hybrid molecular sphere
newton on
region box block -5 5 -5 5 -5 5
create_box 4 box angle/types 1 extra/angle/per/atom 10
angle_style zero
angle_coeff * 75.0
molecule Cube data.Cube
neighbor 0.2 multi
comm_modify mode multi vel yes
pair_style granular
pair_coeff * * hertz 6.542e9 0.55 tangential linear_history 3.822e9 1.0 0.39 damping coeff_restitution limit_damping
create_atoms 0 single 0 0 0 mol Cube 1729871 units box
create_atoms 2 single 3 0 0 units box
group owning id 1
group cube molecule 1
group surface subtract cube owning
group projectile id 388
set group owning type 4
set group projectile diameter 2
velocity projectile set -5000 0 0
neigh_modify collection/type 2 1*3 4
comm_modify cutoff/multi 1 0.0 cutoff/multi 2 5.0
neigh_modify exclude molecule/intra all
neigh_modify exclude group owning surface
neigh_modify exclude molecule/inter owning
fix projectile projectile nve/sphere
fix abrade cube rigid/abrade 1e8 0.25 1.42 molecule remesh 0.95{extra} normals stationary
{output_commands(shape, 100)}timestep 5e-7
thermo 100
thermo_style custom step atoms
run 1000
"""


def frames(path):
    """Parse local and atom dumps, preserving timesteps and exact numeric rows."""
    lines = path.read_text().splitlines()
    result = {}
    i = 0
    while i < len(lines):
        if lines[i] != "ITEM: TIMESTEP":
            raise AssertionError(f"Bad dump frame in {path}: {lines[i]}")
        step = int(lines[i + 1])
        count = int(lines[i + 3])
        rows = [[float(x) for x in row.split()] for row in lines[i + 9:i + 9 + count]]
        if len(rows) != count:
            raise AssertionError(f"Incomplete frame in {path}")
        result[step] = rows
        i += 9 + count
    return result


def close(actual, expected, *, atol=1e-8, rtol=1e-8):
    if not math.isclose(actual, expected, abs_tol=atol, rel_tol=rtol):
        raise AssertionError(f"{actual:.17g} != expected {expected:.17g}")


def check_box(path):
    result = frames(path / "shape.dump")
    expected = [4., 3., 2., math.pi**(1/3) * 144**(2/3) / 52, 1.,
                2/3, 3/4, (3/8)**(1/3), 2/9, 1/9, 1/3, 1/4,
                9/17-1/3, 8/17-1/3]
    for step, rows in result.items():
        if len(rows) != 1 or len(rows[0]) != 15:
            raise AssertionError(f"Wrong shape output dimensions at step {step}")
        close(rows[0][0], 1)
        for actual, target in zip(rows[0][1:], expected):
            close(actual, target)
    for rows in frames(path / "single.dump").values():
        close(rows[0][0], 4.)
    return result


def compare_frames(a, b, atol=0., rtol=0.):
    left, right = frames(a), frames(b)
    if left.keys() != right.keys():
        raise AssertionError(f"Different output steps: {a}, {b}")
    for step in left:
        if len(left[step]) != len(right[step]):
            raise AssertionError(f"Different row counts at step {step}")
        for xrow, yrow in zip(left[step], right[step]):
            if len(xrow) != len(yrow):
                raise AssertionError(f"Different column counts at step {step}")
            for x, y in zip(xrow, yrow):
                close(x, y, atol=atol, rtol=rtol)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lmp", required=True, type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--include-remesh", action="store_true",
                        help="also exercise actual equalise/remesh; records known baseline failures")
    args = parser.parse_args()
    args.lmp = args.lmp.resolve(strict=True)
    if args.baseline:
        args.baseline = args.baseline.resolve(strict=True)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        parser.error("--output-dir must be empty to preserve previous results")
    checks = []
    known_failures = []
    remesh_validated = False
    node_deletion_output_validated = False

    def run(name, executable, text, ranks=1, restart_from=None, allow_failure=False):
        case = output / name
        case.mkdir()
        box_molecule(case / "data.Box")
        shutil.copyfile(HERE.parent / "Stationary_Sphere_Impact/data.Cube", case / "data.Cube")
        if restart_from:
            for filename in ("state.10", "state.10.rigid"):
                shutil.copyfile(restart_from / filename, case / filename)
        (case / "in.test").write_text(text)
        command = [str(executable), "-in", "in.test", "-log", "log.lammps"]
        if ranks > 1:
            command = [args.mpiexec, "-n", str(ranks)] + command
        begin = time.monotonic()
        completed = subprocess.run(command, cwd=case, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, timeout=120)
        (case / "screen.txt").write_text(completed.stdout)
        elapsed = time.monotonic() - begin
        if completed.returncode and not allow_failure:
            raise AssertionError(f"{name} failed ({completed.returncode}):\n{completed.stdout[-5000:]}")
        checks.append({"case": name, "ranks": ranks, "returncode": completed.returncode,
                       "wall_seconds": elapsed})
        return case, completed

    for offset in (False, True):
        name = "box_offset" if offset else "box"
        one, _ = run(name, args.lmp, box_input(True, offset=offset))
        two, _ = run(name + "_mpi2", args.lmp, box_input(True, offset=offset), ranks=2)
        check_box(one)
        check_box(two)
        compare_frames(one / "shape.dump", two / "shape.dump", atol=1e-8, rtol=1e-8)

    moving, _ = run("moving", args.lmp, box_input(True, moving=True))
    check_box(moving)
    moving2, _ = run("moving_mpi2", args.lmp, box_input(True, moving=True), ranks=2)
    check_box(moving2)
    compare_frames(moving / "shape.dump", moving2 / "shape.dump", atol=1e-7, rtol=1e-7)
    restarted, _ = run("restart", args.lmp, box_input(True, restart=True), restart_from=moving)
    check_box(restarted)
    full, resumed = frames(moving / "shape.dump"), frames(restarted / "shape.dump")
    for step in resumed:
        for left, right in zip(full[step][0], resumed[step][0]):
            close(left, right, atol=1e-7, rtol=1e-7)

    impact, _ = run("impact_shape", args.lmp, impact_input(True))
    history = frames(impact / "old.dump")
    if not history[1000][0][1] < history[0][0][1]:
        raise AssertionError("Impact did not produce measurable abrasion")
    for step, rows in frames(impact / "shape.dump").items():
        row = rows[0]
        if row[5] != 1 or not row[1] >= row[2] >= row[3] > 0:
            raise AssertionError(f"Invalid abrasion OBB at step {step}: {row}")
        volume, area = history[step][0][1:3]
        close(row[4], math.pi**(1/3) * (6*volume)**(2/3)/area, atol=1e-12, rtol=1e-12)
    impact_old, _ = run("impact_old", args.lmp, impact_input(False))
    compare_frames(impact / "old.dump", impact_old / "old.dump")
    compare_frames(impact / "atoms.dump", impact_old / "atoms.dump")
    if args.baseline:
        baseline, _ = run("impact_baseline", args.baseline, impact_input(False))
        compare_frames(impact_old / "old.dump", baseline / "old.dump")
        compare_frames(impact_old / "atoms.dump", baseline / "atoms.dump")

    # Dynamic trajectories can differ between processor counts in the original
    # integrator. Compare enabled/disabled/original at the SAME processor count.
    impact2, _ = run("impact_shape_mpi2", args.lmp, impact_input(True), ranks=2)
    impact2_old, _ = run("impact_old_mpi2", args.lmp, impact_input(False), ranks=2)
    compare_frames(impact2 / "old.dump", impact2_old / "old.dump")
    compare_frames(impact2 / "atoms.dump", impact2_old / "atoms.dump")
    if args.baseline:
        baseline2, _ = run("impact_baseline_mpi2", args.baseline, impact_input(False), ranks=2)
        compare_frames(impact2_old / "old.dump", baseline2 / "old.dump")
        compare_frames(impact2_old / "atoms.dump", baseline2 / "atoms.dump")
    impact2_counts = [len(rows) for rows in frames(impact2 / "atoms.dump").values()]
    for step, rows in frames(impact2 / "shape.dump").items():
        if rows[0][5] != 1 or not rows[0][1] >= rows[0][2] >= rows[0][3] > 0:
            raise AssertionError(f"Invalid distributed abrasion/remesh OBB at step {step}")
    node_deletion_output_validated = min(impact2_counts) < impact2_counts[0]

    if args.include_remesh:
        remeshed, new_result = run("remesh_shape", args.lmp, impact_input(True, equalise=True),
                                   allow_failure=True)
        if new_result.returncode:
            if not args.baseline:
                raise AssertionError("Remesh failed; supply --baseline to establish whether this is pre-existing")
            _, base_result = run("remesh_baseline", args.baseline, impact_input(False, equalise=True),
                                 allow_failure=True)
            marker = "memcpy argument memory ranges overlap"
            if not base_result.returncode or marker not in base_result.stdout or marker not in new_result.stdout:
                raise AssertionError("Remesh failure does not match the known original MPI datatype failure")
            known_failures.append("Original rigid/abrade equalise/remesh aborts in MPICH before shape output; "
                                  "MPI_Allreduce uses MPI_LMP_BIGINT on int counts. Full remesh validation remains blocked.")
        else:
            counts = [len(rows) for rows in frames(remeshed / "atoms.dump").values()]
            if min(counts) >= 388:
                raise AssertionError("Remesh case did not delete any surface nodes")
            for rows in frames(remeshed / "shape.dump").values():
                if rows[0][5] != 1:
                    raise AssertionError("Invalid OBB after remeshing")
            remesh_validated = True

    report = {"checks": checks, "known_baseline_failures": known_failures,
              "passed": not known_failures, "remesh_validated": remesh_validated,
              "node_deletion_output_validated": node_deletion_output_validated,
              "impact_mpi2_atom_count_initial": impact2_counts[0],
              "impact_mpi2_atom_count_minimum": min(impact2_counts),
              "baseline": str(args.baseline) if args.baseline else None,
              "lmp": str(args.lmp)}
    (output / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Completed {len(checks)} integration runs; regular checks passed; "
          f"{len(known_failures)} known baseline limitation(s).")
    for message in known_failures:
        print(message)
    print(output / "validation.json")


if __name__ == "__main__":
    main()
