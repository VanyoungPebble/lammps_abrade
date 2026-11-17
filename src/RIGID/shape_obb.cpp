/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Face/edge search adapted from minboundbox('v',3), as used by SHAPE:
   Johannes Korsawe (2015), with the 2D rectangle method of John D'Errico.
   The incremental 3D hull and 2D monotone-chain hull below are independent
   implementations, so no MATLAB, Qhull or other runtime library is needed.

   Copyright (c) 2015, Johannes Korsawe
   Copyright (c) 2012, John D'Errico
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   * Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.
   * Neither the name of the Volkswagen AG nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
------------------------------------------------------------------------- */

#include "shape_obb.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <utility>

namespace LAMMPS_NS {
namespace ShapeOBB {
namespace {

// All geometric tests below operate on coordinates scaled to unit extent.
constexpr double TOL = 1.0e-12;
using Point2 = std::array<double, 2>;

Point sub(const Point &a, const Point &b)
{
  return {{a[0] - b[0], a[1] - b[1], a[2] - b[2]}};
}

double dot(const Point &a, const Point &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Point cross(const Point &a, const Point &b)
{
  return {{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
           a[0] * b[1] - a[1] * b[0]}};
}

bool normalize(Point &a)
{
  const double norm = std::sqrt(dot(a, a));
  if (!(norm > TOL) || !std::isfinite(norm)) return false;
  for (double &x : a) x /= norm;
  return true;
}

struct Face {
  int a, b, c;
  Point normal;
};

bool make_face(int a, int b, int c, const std::vector<Point> &p, const Point &inside,
               Face &face)
{
  Point normal = cross(sub(p[b], p[a]), sub(p[c], p[a]));
  // This vector has area units.  A distance cutoff here would reject small
  // but well-defined facets on otherwise three-dimensional dense surfaces.
  const double area_norm = std::hypot(std::hypot(normal[0], normal[1]), normal[2]);
  if (!(area_norm > TOL * TOL) || !std::isfinite(area_norm)) return false;
  for (double &value : normal) value /= area_norm;
  if (dot(normal, sub(inside, p[a])) > 0.0) {
    std::swap(b, c);
    for (double &x : normal) x = -x;
  }
  face = {a, b, c, normal};
  return true;
}

// Maintain outward triangles. The initial tetrahedron's centroid remains
// inside every subsequent hull, giving a fixed and unambiguous orientation.
bool convex_hull(const std::vector<Point> &p, std::vector<Face> &faces)
{
  const int n = static_cast<int>(p.size());
  int seed[4] = {0, -1, -1, -1};
  double best = 0.0;
  for (int i = 1; i < n; ++i) {
    const Point d = sub(p[i], p[0]);
    const double value = dot(d, d);
    if (value > best) best = value, seed[1] = i;
  }
  if (seed[1] < 0 || best <= TOL * TOL) return false;
  Point line = sub(p[seed[1]], p[0]);
  if (!normalize(line)) return false;
  best = 0.0;
  for (int i = 1; i < n; ++i) {
    const Point d = cross(line, sub(p[i], p[0]));
    const double value = dot(d, d);
    if (value > best) best = value, seed[2] = i;
  }
  if (seed[2] < 0 || best <= TOL * TOL) return false;
  Point normal = cross(line, sub(p[seed[2]], p[0]));
  if (!normalize(normal)) return false;
  best = 0.0;
  for (int i = 1; i < n; ++i) {
    const double value = std::abs(dot(normal, sub(p[i], p[0])));
    if (value > best) best = value, seed[3] = i;
  }
  if (seed[3] < 0 || best <= TOL) return false;

  Point inside = {{0.0, 0.0, 0.0}};
  for (int i : seed)
    for (int k = 0; k < 3; ++k) inside[k] += p[i][k] * 0.25;
  const int initial[4][3] = {{0, 1, 2}, {0, 3, 1}, {0, 2, 3}, {1, 3, 2}};
  for (const auto &s : initial) {
    Face face;
    if (!make_face(seed[s[0]], seed[s[1]], seed[s[2]], p, inside, face)) return false;
    faces.push_back(face);
  }

  // Insert distant vertices first.  Lexicographic insertion of dense grids on
  // almost planar faces builds very skinny intermediate triangles, whose
  // normals amplify coordinate round-off.  Establish the outer corners first
  // instead; every remaining point is still tested against the resulting hull.
  std::vector<int> order;
  order.reserve(n);
  for (int i = 0; i < n; ++i) order.push_back(i);
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    const Point da = sub(p[a], inside), db = sub(p[b], inside);
    const double na = dot(da, da), nb = dot(db, db);
    return na != nb ? na > nb : a < b;
  });
  for (int i : order) {
    if (std::find(seed, seed + 4, i) != seed + 4) continue;
    std::vector<Face> retained;
    std::map<std::pair<int, int>, int> edges;
    bool outside = false;
    for (const Face &face : faces) {
      if (dot(face.normal, sub(p[i], p[face.a])) <= TOL) {
        retained.push_back(face);
        continue;
      }
      outside = true;
      const int vertices[3] = {face.a, face.b, face.c};
      for (int k = 0; k < 3; ++k) {
        const int a = vertices[k], b = vertices[(k + 1) % 3];
        ++edges[std::minmax(a, b)];
      }
    }
    if (!outside) continue;
    int added = 0;
    for (const auto &entry : edges) {
      if (entry.second == 2) continue;
      if (entry.second != 1) return false;
      Face face;
      if (!make_face(entry.first.first, entry.first.second, i, p, inside, face)) return false;
      retained.push_back(face);
      ++added;
    }
    if (added < 3) return false;
    faces.swap(retained);
  }

  // Reject a numerical hull failure rather than silently returning a box of
  // an incomplete surface. This also catches open/nonmanifold triangulations.
  std::map<std::pair<int, int>, int> edges;
  for (const Face &face : faces) {
    const int vertices[3] = {face.a, face.b, face.c};
    for (int k = 0; k < 3; ++k) ++edges[std::minmax(vertices[k], vertices[(k + 1) % 3])];
    for (const Point &point : p)
      if (dot(face.normal, sub(point, p[face.a])) > 16.0 * TOL) return false;
  }
  for (const auto &entry : edges)
    if (entry.second != 2) return false;
  return true;
}

double cross2(const Point2 &a, const Point2 &b, const Point2 &c)
{
  return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

std::vector<Point2> hull2(std::vector<Point2> p)
{
  std::sort(p.begin(), p.end());
  p.erase(std::unique(p.begin(), p.end()), p.end());
  if (p.size() < 3) return p;
  std::vector<Point2> hull(2 * p.size());
  std::size_t k = 0;
  for (const Point2 &point : p) {
    while (k >= 2 && cross2(hull[k - 2], hull[k - 1], point) <= 0.0) --k;
    hull[k++] = point;
  }
  const std::size_t lower = k;
  for (std::size_t i = p.size() - 1; i > 0; --i) {
    const Point2 &point = p[i - 1];
    while (k > lower && cross2(hull[k - 2], hull[k - 1], point) <= 0.0) --k;
    hull[k++] = point;
  }
  hull.resize(k - 1);
  return hull;
}

// For one normal, search every projected hull edge for the minimum-area
// rectangle. Its product with the fixed normal extent minimizes box volume.
void check_frame(const std::vector<Point> &points, const Point &u, const Point &v,
                 const Point &w, double &best_volume, std::array<double, 3> &best_length)
{
  std::vector<Point2> projected;
  projected.reserve(points.size());
  double zmin = std::numeric_limits<double>::infinity();
  double zmax = -zmin;
  for (const Point &point : points) {
    projected.push_back({{dot(point, u), dot(point, v)}});
    const double z = dot(point, w);
    zmin = std::min(zmin, z);
    zmax = std::max(zmax, z);
  }
  const double depth = zmax - zmin;
  if (!(depth > TOL)) return;
  const std::vector<Point2> hull = hull2(std::move(projected));
  if (hull.size() < 3) return;
  double best_area = std::numeric_limits<double>::infinity();
  std::array<double, 2> best_width = {{0.0, 0.0}};
  for (std::size_t i = 0; i < hull.size(); ++i) {
    const Point2 &a = hull[i], &b = hull[(i + 1) % hull.size()];
    const double dx = b[0] - a[0], dy = b[1] - a[1];
    const double norm = std::hypot(dx, dy);
    if (!(norm > TOL)) continue;
    const double c = dx / norm, s = dy / norm;
    double xmin = std::numeric_limits<double>::infinity(), ymin = xmin;
    double xmax = -xmin, ymax = -xmin;
    for (const Point2 &point : hull) {
      const double x = point[0] * c + point[1] * s;
      const double y = -point[0] * s + point[1] * c;
      xmin = std::min(xmin, x);
      xmax = std::max(xmax, x);
      ymin = std::min(ymin, y);
      ymax = std::max(ymax, y);
    }
    const double xwidth = xmax - xmin, ywidth = ymax - ymin;
    const double area = xwidth * ywidth;
    if (xwidth > TOL && ywidth > TOL && area < best_area) {
      best_area = area;
      best_width = {{xwidth, ywidth}};
    }
  }
  const double volume = best_area * depth;
  if (!std::isfinite(volume)) return;
  std::array<double, 3> length = {{best_width[0], best_width[1], depth}};
  std::sort(length.begin(), length.end());
  if (volume < best_volume) {
    best_volume = volume;
    best_length = length;
  }
}

}    // namespace

Result min_volume(const std::vector<Point> &input)
{
  Result result;
  if (input.size() < 4 || input.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    return result;
  Point low = input[0], high = input[0];
  for (const Point &point : input)
    for (int k = 0; k < 3; ++k) {
      if (!std::isfinite(point[k])) return result;
      low[k] = std::min(low[k], point[k]);
      high[k] = std::max(high[k], point[k]);
    }
  double scale = 0.0;
  for (int k = 0; k < 3; ++k) scale = std::max(scale, high[k] - low[k]);
  if (!(scale > 0.0) || !std::isfinite(scale)) return result;

  std::vector<Point> points;
  points.reserve(input.size());
  for (const Point &point : input) {
    Point normalized;
    for (int k = 0; k < 3; ++k) normalized[k] = (point[k] - low[k]) / scale;
    points.push_back(normalized);
  }
  std::sort(points.begin(), points.end());
  points.erase(std::unique(points.begin(), points.end()), points.end());
  if (points.size() < 4) return result;
  std::vector<Face> faces;
  if (!convex_hull(points, faces)) return result;

  std::set<std::pair<int, int>> edges;
  for (const Face &face : faces) {
    edges.insert(std::minmax(face.a, face.b));
    edges.insert(std::minmax(face.b, face.c));
    edges.insert(std::minmax(face.c, face.a));
  }
  double best_volume = std::numeric_limits<double>::infinity();
  std::array<double, 3> best_length = {{0.0, 0.0, 0.0}};
  for (const Face &face : faces) {
    Point u = sub(points[face.b], points[face.a]);
    if (!normalize(u)) continue;
    Point v = cross(face.normal, u);
    if (!normalize(v)) continue;
    // Hull topology determines candidates, but include every input point in
    // the extents so near-hull points cannot be clipped by hull tolerances.
    check_frame(points, u, v, face.normal, best_volume, best_length);
  }
  for (const auto &edge : edges) {
    Point u = sub(points[edge.first], points[edge.second]);
    Point v = {{u[1], -u[0], 0.0}};
    if (v[0] == 0.0 && v[1] == 0.0) v = {{u[2], 0.0, -u[0]}};
    if (!normalize(u) || !normalize(v)) continue;
    Point w = cross(u, v);
    if (!normalize(w)) continue;
    check_frame(points, u, v, w, best_volume, best_length);
    check_frame(points, v, w, u, best_volume, best_length);
    check_frame(points, w, u, v, best_volume, best_length);
  }
  if (!std::isfinite(best_volume)) return result;
  for (int k = 0; k < 3; ++k) {
    result.length[k] = best_length[k] * scale;
    if (!(result.length[k] > 0.0) || !std::isfinite(result.length[k])) return Result{};
  }
  result.valid = true;
  return result;
}

}    // namespace ShapeOBB
}    // namespace LAMMPS_NS
