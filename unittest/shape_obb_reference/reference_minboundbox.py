#!/usr/bin/env python3
"""Independent numerical reference for SHAPE minboundbox(metric='v', level=3).

This is a NumPy/SciPy translation of the shipped MATLAB algorithm, including
euler123 and its three cyclic edge frames. It is NOT a MATLAB/Octave execution.
The frozen JSON fixtures can also be checked against the standalone C++ driver.
Only fixture regeneration requires NumPy and SciPy.

Derived from minboundbox.m by Johannes Korsawe and John D'Errico. See the
accompanying LICENSE.minboundbox for the required copyright and BSD terms.
"""

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent


def reference(points):
    import numpy as np
    from scipy.spatial import ConvexHull

    xyz = np.asarray(points, dtype=float)
    hull = ConvexHull(xyz, qhull_options="Qt")
    triangles = hull.simplices
    best = math.inf
    lengths = None

    def euler123(v1, v2, nv):
        beta = math.asin(float(np.clip(nv[0], -1.0, 1.0)))
        cb = math.cos(beta)
        if nv[0] == 1:
            alpha = math.asin(float(np.clip(v2[2], -1.0, 1.0)))
            gamma = 0.0
        else:
            alpha = math.acos(float(np.clip(nv[2] / cb, -1.0, 1.0)))
            if np.sign(nv[1]) != np.sign(-math.sin(alpha) * cb):
                alpha = -alpha
            sg = float(np.clip(-v2[0] / cb, -1.0, 1.0))
            gamma = math.asin(sg) if v1[0] >= 0 else -math.pi - math.asin(sg)
        ca, sa = math.cos(alpha), math.sin(alpha)
        cg, sg = math.cos(gamma), math.sin(gamma)
        sb = math.sin(beta)
        return np.array([
            [cb * cg, -cb * sg, sb],
            [sa * sb * cg + ca * sg, -sa * sb * sg + ca * cg, -sa * cb],
            [-ca * sb * cg + sa * sg, ca * sb * sg + sa * cg, ca * cb],
        ])

    def checkbox(v1, v2, nv):
        nonlocal best, lengths
        rot = euler123(v1, v2, nv)
        projected = xyz @ rot
        boundary = ConvexHull(projected[:, :2]).vertices
        boundary = np.r_[boundary, boundary[:1]]
        delta = np.diff(projected[boundary, :2], axis=0)
        angles = np.unique(np.mod(np.arctan2(delta[:, 1], delta[:, 0]), math.pi / 2))
        best_area = math.inf
        rot2 = None
        for angle in angles:
            c, s = math.cos(-angle), math.sin(-angle)
            candidate = np.array([[c, s], [-s, c]])
            extent = np.ptp(projected[boundary, :2] @ candidate, axis=0)
            area = float(np.prod(extent))
            if area < best_area:
                best_area, rot2 = area, candidate
        block = np.eye(3)
        block[:2, :2] = rot2
        extent = np.ptp(xyz @ (rot @ block), axis=0)
        volume = float(np.prod(extent))
        if volume < best:
            best, lengths = volume, sorted(map(float, extent))

    # Face candidates: exactly the Gram-Schmidt construction used in SHAPE.
    for i, j, k in triangles:
        v1 = xyz[j] - xyz[i]
        v1 /= np.linalg.norm(v1)
        v2 = xyz[k] - xyz[i]
        v2 -= np.dot(v1, v2) * v1
        v2 /= np.linalg.norm(v2)
        checkbox(v1, v2, np.cross(v1, v2))

    edges = sorted({tuple(sorted((int(face[a]), int(face[b]))))
                    for face in triangles for a, b in ((0, 1), (1, 2), (0, 2))})
    # Edge candidates, including the coordinate-dependent perpendicular choice.
    for i, j in edges:
        va = xyz[i] - xyz[j]
        vb = np.array([va[1], -va[0], 0.0])
        if np.sum(np.abs(vb)) == 0:
            vb = np.array([va[2], 0.0, -va[0]])
        va /= np.linalg.norm(va)
        vb /= np.linalg.norm(vb)
        nv = np.cross(va, vb)
        checkbox(va, vb, nv)
        checkbox(vb, nv, va)
        checkbox(nv, va, vb)
    return {"length_s_i_l": lengths, "volume": best}


def build_cases():
    import numpy as np
    from scipy.spatial.transform import Rotation

    rng = np.random.default_rng(20261116)
    rotation = Rotation.from_rotvec([0.31, -0.57, 0.23]).as_matrix()
    box = np.array(list(itertools.product((-0.5, 0.5), repeat=3))) * [1.3, 2.1, 3.7]
    ellipsoid = rng.normal(size=(25, 3))
    ellipsoid /= np.linalg.norm(ellipsoid, axis=1)[:, None]
    ellipsoid *= [0.7, 1.4, 2.6]
    cloud = rng.normal(size=(19, 3)) * [0.8, 1.2, 1.6]
    # SHAPE's additional perpendicular edge frames use global coordinate
    # directions. Keep one case demonstrating the resulting heuristic frame
    # dependence, rather than inadvertently testing an idealized algorithm.
    frame_rng = np.random.default_rng(901)
    for _ in range(2):
        frame_cloud = frame_rng.normal(size=(12, 3)) * [1, 1.5, 2]
        frame_rotated = frame_cloud @ Rotation.from_rotvec(frame_rng.normal(size=3)).as_matrix().T
    candidates = {
        "axis_aligned_box": box,
        "rotated_translated_box": box @ rotation.T + [4.3, -2.7, 1.6],
        "shape_documented_tetrahedron": np.array([[1, 0, 0], [0, 1, 0],
                                                   [0.1, 0.1, 0.9], [1, 1, 1]]),
        "asymmetric_cloud": cloud,
        "asymmetric_cloud_interior_duplicates": np.vstack([cloud, cloud.mean(axis=0), cloud[:3]]),
        "sampled_ellipsoid": ellipsoid @ rotation.T,
        "coordinate_sensitive_cloud": frame_cloud,
        "coordinate_sensitive_cloud_rotated": frame_rotated,
    }
    return [{"name": name, "points": points.tolist(), **reference(points)}
            for name, points in candidates.items()]


def check_kernel(executable, cases):
    def invoke(points):
        data = str(len(points)) + "\n" + "\n".join(
            " ".join(format(value, ".17g") for value in point) for point in points) + "\n"
        result = subprocess.run([str(executable)], input=data, text=True,
                                capture_output=True, check=True, timeout=30)
        fields = result.stdout.strip().split()
        if len(fields) != 4:
            raise AssertionError(f"Unexpected driver output: {result.stdout!r}")
        return bool(int(fields[0])), list(map(float, fields[1:]))

    count = 0
    for case in cases:
        valid, actual = invoke(case["points"])
        expected = case["length_s_i_l"]
        if not valid or not all(math.isclose(a, e, rel_tol=2e-8, abs_tol=2e-10)
                                for a, e in zip(actual, expected)):
            raise AssertionError(f"{case['name']}: actual={actual}, expected={expected}, valid={valid}")
        print(f"PASS {case['name']}: S/I/L={actual}")
        count += 1
        # Exercise scale handling without large-translation loss of significance.
        for scale in (1e-9, 1e9):
            valid, actual = invoke([[value * scale for value in p] for p in case["points"]])
            if not valid or not all(math.isclose(a / scale, e, rel_tol=2e-8, abs_tol=2e-10)
                                    for a, e in zip(actual, expected)):
                raise AssertionError(f"{case['name']} scale={scale}: {actual}, expected={expected}")
            count += 1

    degenerate = [[], [[0, 0, 0]], [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]],
                  [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
                  [[0, 0, 0], [math.nan, 1, 0], [1, 0, 0], [0, 0, 1]]]
    for points in degenerate:
        valid, _ = invoke(points)
        if valid:
            raise AssertionError(f"Degenerate cloud accepted: {points}")
        count += 1

    # The repository's 386-node STL cube has densely sampled, nearly coplanar
    # faces (approximately 1e-7 coordinate perturbations). This is a much more
    # demanding hull regression than an eight-corner box. Exercise it after
    # rigid rotations without duplicating the existing molecule data file.
    molecule = HERE.parents[1] / "examples/rigid/Fix_Rigid_Abrade/Stationary_Sphere_Impact/data.Cube"
    cube = []
    in_coords = False
    for line in molecule.read_text().splitlines():
        row = line.split()
        if row == ["Coords"]:
            in_coords = True
            continue
        if not in_coords or not row or row[0].startswith("#"):
            continue
        if len(row) != 4 or not row[0].isdigit():
            break
        if row[0] != "1":
            cube.append(list(map(float, row[1:])))
    if len(cube) != 386:
        raise AssertionError(f"Unexpected dense cube fixture: {len(cube)} vertices")
    angles = [(0.0, 0.0, 0.0)] + [(.23 + .11 * i, -.17 - .07 * i, .31 + .13 * i)
                                      for i in range(8)]
    for i, (alpha, beta, gamma) in enumerate(angles):
        ca, sa, cb, sb, cg, sg = (math.cos(alpha), math.sin(alpha), math.cos(beta),
                                math.sin(beta), math.cos(gamma), math.sin(gamma))
        rotated = []
        for x, y, z in cube:
            y, z = ca * y - sa * z, sa * y + ca * z
            x, z = cb * x + sb * z, -sb * x + cb * z
            x, y = cg * x - sg * y, sg * x + cg * y
            rotated.append([x, y, z])
        valid, actual = invoke(rotated)
        if not valid or not all(math.isclose(value, 1.61, rel_tol=0, abs_tol=5e-7)
                                for value in actual):
            raise AssertionError(f"Dense nearcoplanar cube rotation {i}: {actual}, valid={valid}")
        count += 1
    print("PASS 9 dense nearcoplanar 386-vertex cube orientations")
    print(f"PASS {count} standalone OBB checks")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regenerate", action="store_true")
    parser.add_argument("--shape-source", type=Path)
    parser.add_argument("--kernel", type=Path)
    args = parser.parse_args()
    fixture_path = HERE / "reference_cases.json"
    if args.regenerate:
        payload = {"provenance": "Independent NumPy/SciPy translation; not MATLAB execution",
                   "algorithm": "SHAPE minboundbox(metric='v', level=3)", "cases": build_cases()}
        if args.shape_source:
            payload["source_sha256"] = hashlib.sha256(args.shape_source.read_bytes()).hexdigest()
        fixture_path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
        print(f"Wrote {len(payload['cases'])} reference cases to {fixture_path}")
    if args.kernel:
        check_kernel(args.kernel.resolve(), json.loads(fixture_path.read_text())["cases"])


if __name__ == "__main__":
    main()
