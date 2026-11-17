/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_SHAPE_OBB_H
#define LMP_SHAPE_OBB_H

#include <array>
#include <vector>

namespace LAMMPS_NS {
namespace ShapeOBB {

using Point = std::array<double, 3>;

struct Result {
  bool valid = false;
  double length[3] = {0.0, 0.0, 0.0};  // full dimensions, sorted S <= I <= L
};

// The face/edge candidate search used by SHAPE's minboundbox('v',3).
// Finite, non-coplanar points are required. No input coordinates are changed.
Result min_volume(const std::vector<Point> &points);

}    // namespace ShapeOBB
}    // namespace LAMMPS_NS

#endif
