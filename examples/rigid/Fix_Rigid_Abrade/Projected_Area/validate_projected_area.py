"""Isolated analytic/MPI checks; all generated runs go under --out."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import time


HERE = Path(__file__).resolve().parent


def read_dump(path):
    lines = path.read_text().splitlines()
    frames, cursor = {}, 0
    while cursor < len(lines):
        if lines[cursor] != 'ITEM: TIMESTEP':
            raise ValueError(f'Unexpected dump header in {path}')
        step, count = int(lines[cursor+1]), int(lines[cursor+3])
        if step in frames or cursor+9+count > len(lines):
            raise ValueError(f'Duplicate or incomplete frame in {path}')
        frames[step] = [list(map(float, row.split()))
                        for row in lines[cursor+9:cursor+9+count]]
        cursor += 9+count
    return frames


def close(actual, expected, label, tolerance=3e-7):
    if not math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance):
        raise AssertionError(f'{label}: {actual} != {expected}')


def box_reference(rows, direction):
    # Exact parallelepiped projection from its three actual edge vectors.
    # rigid/abrade reinitializes omega, so prescribed nodal velocities are not
    # a valid substitute for the resulting orientation.
    positions = {int(row[0]): row[2:5] for row in rows}
    edges = [[positions[tag][k]-positions[2][k] for k in range(3)] for tag in (3,5,6)]
    for edge, length in zip(edges, (2,3,4)):
        close(math.sqrt(sum(x*x for x in edge)),length,'cuboid edge length',tolerance=1e-9)
    total = 0.0
    for a,b in ((edges[0],edges[1]),(edges[1],edges[2]),(edges[2],edges[0])):
        close(sum(x*y for x,y in zip(a,b)),0,'orthogonal cuboid edges',tolerance=1e-9)
        normal = [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
        total += abs(sum(x*y for x,y in zip(normal,direction)))
    return total/math.sqrt(sum(x*x for x in direction))


def with_measurements(template, directions):
    lines = [line for line in template.splitlines()
             if not line.startswith(('compute ap ', 'dump projected ', 'dump_modify projected '))]
    block = []
    for name, direction in directions.items():
        block += [f'compute {name} owned projected/area ' + ' '.join(map(str, direction)),
                  f'dump d_{name} all local 5 {name}.dump c_{name}[1] c_{name}[2]',
                  f'dump_modify d_{name} sort 1 format float %.17g']
    return '\n'.join(lines).replace('timestep 1e-4', '\n'.join(block)+'\ntimestep 1e-4')+'\n'


def gravel_reference(path):
    coordinates, faces, section = {}, [], None
    for line in path.read_text().splitlines():
        line = line.split('#')[0].strip()
        if not line:
            continue
        if line in ('Coords', 'Types', 'Molecules', 'Diameters', 'Masses', 'Angles'):
            section = line
        elif section == 'Coords':
            row = line.split()
            coordinates[int(row[0])] = list(map(float, row[1:4]))
        elif section == 'Angles':
            faces.append(list(map(int, line.split()[2:5])))
    # Independent reference for this known-convex initial template: facet formula.
    total = 0.0
    for face in faces:
        a, b, c = [coordinates[tag] for tag in face]
        by, bz, cy, cz = b[1]-a[1], b[2]-a[2], c[1]-a[1], c[2]-a[2]
        total += abs(by*cz-bz*cy)/4
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--lmp', type=Path, required=True)
    parser.add_argument('--mpiexec', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--gravel', type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    report = {'passed': False, 'runs': [], 'checks': [],
              'lmp': str(args.lmp.resolve()), 'scope': 'prototype correctness, not a performance benchmark'}
    template = (HERE/'in.projected_area').read_text()
    directions = {'ap': (1,0,0), 'ay': (0,1,0)}
    dynamic = with_measurements(template, directions)
    static = with_measurements(template, {
        'ap': (1,0,0), 'ay': (0,1,0), 'az': (0,0,1),
        'diagonal': (1,1,1), 'negative': (-1,0,0), 'scaled': (7,0,0),
    }).replace('run 20', 'run 0')

    def execute(name, text, ranks, binary=None, data=None, success=True):
        directory = args.out/f'{name}_mpi{ranks}'
        directory.mkdir(parents=True, exist_ok=True)
        (directory/'in.test').write_text(text)
        shutil.copyfile(data or HERE/'data.Box', directory/'data.Box')
        command = [str(args.mpiexec.resolve()), '-n', str(ranks),
                   str((binary or args.lmp).resolve()), '-in', 'in.test']
        started = time.monotonic()
        with (directory/'screen.txt').open('w') as screen:
            process = subprocess.run(command, cwd=directory, stdout=screen,
                                     stderr=subprocess.STDOUT, timeout=60,
                                     env={**os.environ, 'OMP_NUM_THREADS': '1'})
        report['runs'].append({'case': name, 'ranks': ranks, 'exit_code': process.returncode,
                               'directory': str(directory), 'seconds': time.monotonic()-started})
        if (process.returncode == 0) != success:
            raise RuntimeError(f'{name} MPI{ranks}: unexpected exit {process.returncode}\n'
                               +(directory/'screen.txt').read_text()[-2200:])
        print(f'{name} MPI{ranks}: expected exit', flush=True)
        return directory

    try:
        trajectories = {}
        for ranks in (1,2):
            directory = execute('analytic_box', static, ranks)
            for name, expected in {'ap':12, 'ay':8, 'az':6, 'diagonal':26/math.sqrt(3),
                                   'negative':12, 'scaled':12}.items():
                rows = read_dump(directory/f'{name}.dump')[0]
                if len(rows) != 1 or rows[0][0] != 1:
                    raise AssertionError('One row per molecule was not produced')
                close(rows[0][1], expected, name)
            report['checks'].append(f'analytic box/direction normalization MPI{ranks}')

            directory = execute('spin_pbc', dynamic, ranks)
            trajectories[ranks] = directory
            positions = read_dump(directory/'atoms.dump')
            for name, direction in directions.items():
                for step, rows in read_dump(directory/f'{name}.dump').items():
                    close(rows[0][1],box_reference(positions[step],direction),
                          f'{name} actual cuboid orientation step {step}',tolerance=1e-10)
            close(positions[20][0][2], 24, 'unwrapped COM crossed the periodic boundary')
            report['checks'].append(f'rotation and PBC using both yz/xz projections MPI{ranks}')

            multi = static.replace(
                'create_atoms 0 single 0 0 0 mol Box 1729871 rotate 0 0 0 1 units box',
                'create_atoms 0 single -5 0 0 mol Box 1729871 rotate 0 0 0 1 units box\n'
                'create_atoms 0 single 5 0 0 mol Box 1729872 rotate 90 0 0 1 units box'
            ).replace('group owning id 1\n', 'group owning id 1 10\n')
            directory = execute('two_molecules', multi, ranks)
            rows = read_dump(directory/'ap.dump')[0]
            if [row[0] for row in rows] != [1,2]:
                raise AssertionError('Molecule labels were merged or duplicated')
            close(rows[0][1],12,'molecule 1'); close(rows[1][1],8,'molecule 2')
            report['checks'].append(f'separate molecule rows and distributed input MPI{ranks}')

            empty = static.replace('compute ap owned ', 'group empty id 999999\ncompute ap empty ')
            directory = execute('empty_group', empty, ranks)
            if read_dump(directory/'ap.dump')[0] != []:
                raise AssertionError('Empty group did not produce zero rows')
            report['checks'].append(f'empty group MPI{ranks}')

            point_input = '''units cgs
boundary p p p
atom_style molecular
region box block 0 4 0 4 -2 2
create_box 1 box
create_atoms 1 single 1 1 0
create_atoms 1 single 3 1 0
create_atoms 1 single 2 1 0
mass 1 1
set group all mol 1
pair_style zero 1
pair_coeff * *
compute ap all projected/area 0 0 1
dump ap all local 1 ap.dump c_ap[1] c_ap[2]
dump_modify ap sort 1 format float %.17g
thermo_style custom step atoms
run 0
'''
            directory = execute('collinear', point_input, ranks)
            if not math.isnan(read_dump(directory/'ap.dump')[0][0][1]):
                raise AssertionError('Collinear Qhull failure did not produce NaN')
            if 'writing NaN' not in (directory/'screen.txt').read_text():
                raise AssertionError('Missing diagnostic for degenerate points')
            report['checks'].append(f'Qhull failure/NaN without random perturbation MPI{ranks}')

            changing = point_input.replace('create_atoms 1 single 2 1 0',
                                          'create_atoms 1 single 3 3 0\ncreate_atoms 1 single 1 3 0')
            changing += '''group corner id 3
displace_atoms corner move -1 -1 0 units box
run 1
group removed id 4
delete_atoms group removed compress no
run 1
'''
            directory = execute('move_delete_points', changing, ranks)
            frames = read_dump(directory/'ap.dump')
            for step, expected in [(0,4), (1,2), (2,1)]:
                close(frames[step][0][1], expected, f'updated point set step {step}')
            report['checks'].append(f'fresh coordinates and atom deletion MPI{ranks}')

            if args.baseline:
                plain = with_measurements(template, {})
                old = execute('baseline_no_measurement', plain, ranks, binary=args.baseline)
                new = execute('new_no_measurement', plain, ranks)
                if not ((old/'atoms.dump').read_bytes() == (new/'atoms.dump').read_bytes()
                        == (trajectories[ranks]/'atoms.dump').read_bytes()):
                    raise AssertionError('Enabling projected area changed positions/velocities/forces')
                report['checks'].append(f'baseline/new/measurement trajectories byte-identical MPI{ranks}')

            if args.gravel:
                directory = execute('actual_gravel_template', static, ranks, data=args.gravel)
                area = read_dump(directory/'ap.dump')[0][0][1]
                expected = gravel_reference(args.gravel)
                close(area,expected,'known-convex gravel facet reference',tolerance=2e-9)
                report['checks'].append(f'initial gravel area={area:.15g} vs facet reference MPI{ranks}')

        for name in ('ap','ay'):
            one = read_dump(trajectories[1]/f'{name}.dump')
            two = read_dump(trajectories[2]/f'{name}.dump')
            for step in one:
                close(one[step][0][1],two[step][0][1],f'MPI consistency {name} step {step}',tolerance=1e-10)
        report['checks'].append('1/2 rank area agreement')

        for name, original, replacement, expected_message in [
            ('zero_direction', 'projected/area 0 0 1', 'projected/area 0 0 0', 'direction must be finite and nonzero'),
            ('invalid_molecule', 'set group all mol 1', 'set group all mol 0', 'requires positive molecule IDs'),
        ]:
            directory = execute(name, point_input.replace(original,replacement),2,success=False)
            if expected_message not in (directory/'screen.txt').read_text():
                raise AssertionError(f'{name} failed for an unexpected reason')
            report['checks'].append(f'collective error: {name}')
        report['passed'] = True
    except Exception as error:
        report['failure'] = str(error)
        raise
    finally:
        report['binary_sha256'] = hashlib.sha256(args.lmp.read_bytes()).hexdigest()
        (args.out/'validation_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':True,'runs':len(report['runs']),'checks':len(report['checks'])}),flush=True)


if __name__ == '__main__':
    main()
