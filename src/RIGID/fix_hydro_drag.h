/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Distributed under the GNU General Public License; see LICENSE.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(hydro/drag,FixHydroDrag);
// clang-format on
#else

#ifndef LMP_FIX_HYDRO_DRAG_H
#define LMP_FIX_HYDRO_DRAG_H

#include "fix.h"

#include <array>
#include <vector>

namespace LAMMPS_NS {

// Hydrodynamic drag (Holzer & Sommerfeld 2008) on fix rigid/abrade bodies.
// This fix holds the fluid parameters and computes the drag, but it does not act on its own:
// fix rigid/abrade calls add_body_forces() when it sums body forces, so the
// drag goes onto each body's COM force and never into the per-atom forces
// that drive abrasion.

class FixHydroDrag : public Fix {
 public:
  FixHydroDrag(class LAMMPS *, int, char **);
  ~FixHydroDrag() override;
  int setmask() override;
  void init() override;
  double compute_vector(int) override;

  // called by fix rigid/abrade in compute_forces_and_torques()
  void add_body_forces(class FixRigidAbrade *);

 protected:
  double rho;           // fluid mass density
  double nu;            // fluid kinematic viscosity
  double flow[3];       // uniform, steady fluid velocity
  int buoyancy_flag;    // 1 = also add buoyancy, 0 = drag only

  // per-body work data, rebuilt at every call of add_body_forces()

  struct BodyData {
    double urel[3];                              // fluid velocity relative to the body COM
    double umag;                                 // |urel|; 0 = no relative flow, body skipped
    double e1[3], e2[3];                         // plane perpendicular to the flow, body coords
    std::vector<std::array<double, 2>> points;   // surface atoms projected onto (e1, e2)
    double rmax;                                 // largest radius among those atoms
  };
  std::vector<BodyData> bodydata;    // one entry per body owned by this process

  // output
  // global vector: total drag force (3), total buoyancy force (3)
  // local array: one row per owned body, see add_body_forces()

  double ftotal[6], ftotal_all[6];    // this process / summed over all processes
  int force_flag;                     // 1 if ftotal_all is up to date
  int nmax;                           // rows allocated in array_local

  void collect_projected_points(class FixRigidAbrade *);
};

}    // namespace LAMMPS_NS

#endif
#endif
