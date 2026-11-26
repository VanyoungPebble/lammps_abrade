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

namespace LAMMPS_NS {

// Hydrodynamic drag (Holzer & Sommerfeld 2008) on fix rigid/abrade bodies.
// This fix only stores the fluid parameters.  It does not act on its own:
// fix rigid/abrade calls add_body_forces() when it sums body forces, so the
// drag goes onto each body's COM force and never into the per-atom forces
// that drive abrasion.

class FixHydroDrag : public Fix {
 public:
  FixHydroDrag(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;

  // called by fix rigid/abrade in compute_forces_and_torques()
  void add_body_forces(class FixRigidAbrade *);

 protected:
  double rho;           // fluid mass density
  double nu;            // fluid kinematic viscosity
  double flow[3];       // uniform, steady fluid velocity
  int buoyancy_flag;    // 1 = also add buoyancy, 0 = drag only
};

}    // namespace LAMMPS_NS

#endif
#endif
