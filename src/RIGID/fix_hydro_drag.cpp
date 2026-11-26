/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Distributed under the GNU General Public License; see LICENSE.
------------------------------------------------------------------------- */

#include "fix_hydro_drag.h"

#include "comm.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   fix ID group hydro/drag rho <rho> nu <nu> flow <ux> <uy> <uz>
                           [buoyancy yes/no]
------------------------------------------------------------------------- */

FixHydroDrag::FixHydroDrag(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), rho(0.0), nu(0.0), buoyancy_flag(0)
{
  flow[0] = flow[1] = flow[2] = 0.0;

  // rho, nu and flow are required; remember which ones were given

  int rho_set = 0, nu_set = 0, flow_set = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "rho") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix hydro/drag rho", error);
      rho = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      rho_set = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "nu") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix hydro/drag nu", error);
      nu = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      nu_set = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "flow") == 0) {
      if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix hydro/drag flow", error);
      flow[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      flow[1] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      flow[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      flow_set = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg], "buoyancy") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix hydro/drag buoyancy", error);
      buoyancy_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown fix hydro/drag keyword: {}", arg[iarg]);
    }
  }

  if (!rho_set) error->all(FLERR, "Fix hydro/drag requires the rho keyword");
  if (!nu_set) error->all(FLERR, "Fix hydro/drag requires the nu keyword");
  if (!flow_set) error->all(FLERR, "Fix hydro/drag requires the flow keyword");
  if (rho <= 0.0) error->all(FLERR, "Fix hydro/drag rho must be > 0");
  if (nu <= 0.0) error->all(FLERR, "Fix hydro/drag nu must be > 0");
}

/* ----------------------------------------------------------------------
   no mask: this fix is never called by LAMMPS directly, only by
   fix rigid/abrade through add_body_forces()
------------------------------------------------------------------------- */

int FixHydroDrag::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHydroDrag::init()
{
  // first version: the projected area needs all vertices of a body on one
  // process, so only serial runs are supported

  if (comm->nprocs > 1) error->all(FLERR, "Fix hydro/drag currently supports only 1 MPI process");
}

/* ----------------------------------------------------------------------
   add drag (and optionally buoyancy) to the COM force of each body
   skeleton only: no force is added yet
------------------------------------------------------------------------- */

void FixHydroDrag::add_body_forces(FixRigidAbrade * /*rigid*/) {}
