/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Distributed under the GNU General Public License; see LICENSE.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(projected/area,ComputeProjectedArea);
// clang-format on
#else

#ifndef LMP_COMPUTE_PROJECTED_AREA_H
#define LMP_COMPUTE_PROJECTED_AREA_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeProjectedArea : public Compute {
 public:
  ComputeProjectedArea(class LAMMPS *, int, char **);
  ~ComputeProjectedArea() override;
  void init() override;
  void compute_local() override;
  double memory_usage() override;

 private:
  double u[3], v[3];
  int nmax, warned;
};

}    // namespace LAMMPS_NS

#endif
#endif
