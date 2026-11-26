/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Distributed under the GNU General Public License; see LICENSE.
------------------------------------------------------------------------- */

#include "fix_hydro_drag.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix_rigid_abrade.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullError.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <sstream>

using namespace LAMMPS_NS;
using MathConst::MY_PI;

namespace {
using Point = std::array<double, 2>;

/* ----------------------------------------------------------------------
   2d convex hull of points -> enclosed area and perimeter
   copied from projected_hull_area() in compute_projected_area.cpp and
   changed to also return the perimeter and to sort the points before
   removing duplicates; returns false if the hull cannot be computed
   (fewer than 3 distinct points, all points collinear, Qhull error)
------------------------------------------------------------------------- */

bool hull_area_perimeter(std::vector<Point> &points, double &area, double &perimeter)
{
  // remove exact duplicates, e.g. front and back corners of a cube
  // that is face-on to the flow project onto the same point

  std::sort(points.begin(), points.end());
  points.erase(std::unique(points.begin(), points.end()), points.end());
  if (points.size() < 3) return false;

  // shift and scale the points into a unit box to keep Qhull's
  // tolerances independent of units and position

  Point low = points.front(), high = low;
  for (const auto &p : points)
    for (int k = 0; k < 2; k++) {
      low[k] = std::min(low[k], p[k]);
      high[k] = std::max(high[k], p[k]);
    }
  const double scale = std::max(high[0] - low[0], high[1] - low[1]);
  if (!std::isfinite(scale) || scale <= 0.0) return false;

  std::vector<double> coordinates;
  coordinates.reserve(2 * points.size());
  for (const auto &p : points) {
    coordinates.push_back((p[0] - low[0]) / scale);
    coordinates.push_back((p[1] - low[1]) / scale);
  }

  std::ostringstream diagnostics;
  try {
    orgQhull::Qhull hull;
    hull.setErrorStream(&diagnostics);
    hull.disableOutputStream();
    hull.runQhull("hydro/drag", 2, static_cast<int>(points.size()), coordinates.data(), "Qt");
    // in two dimensions Qhull volume() is the enclosed area and area() is the perimeter
    area = hull.volume() * scale * scale;
    perimeter = hull.area() * scale;
    return std::isfinite(area) && area > 0.0 && std::isfinite(perimeter);
  } catch (const orgQhull::QhullError &) {
    return false;
  }
}

/* ----------------------------------------------------------------------
   two unit vectors e1, e2 spanning the plane perpendicular to the unit
   vector dir; same construction as in compute projected/area: use the
   coordinate axis least parallel to dir as a helper
------------------------------------------------------------------------- */

void plane_basis(const double *dir, double *e1, double *e2)
{
  int axis = 0;
  for (int k = 1; k < 3; k++)
    if (std::abs(dir[k]) < std::abs(dir[axis])) axis = k;
  double helper[3] = {0.0, 0.0, 0.0};
  helper[axis] = 1.0;

  MathExtra::cross3(dir, helper, e1);
  MathExtra::norm3(e1);
  MathExtra::cross3(dir, e1, e2);
}
}    // namespace

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

  // the projected area correction uses the atom radii

  if (!atom->radius_flag) error->all(FLERR, "Fix hydro/drag requires atom attribute radius");
}

/* ----------------------------------------------------------------------
   add drag (and optionally buoyancy) to the COM force of each body
   drag only so far: buoyancy is parsed but not applied yet
------------------------------------------------------------------------- */

void FixHydroDrag::add_body_forces(FixRigidAbrade *rigid)
{
  const int nbody = rigid->nlocal_body;
  bodydata.resize(nbody);

  // 1. per body: relative velocity, and the plane (e1, e2) perpendicular
  //    to the flow, in body coordinates
  //    the vertices (displace) are stored in body coordinates, so it is
  //    cheaper to rotate the one flow direction than all the vertices

  for (int ibody = 0; ibody < nbody; ibody++) {
    FixRigidAbrade::Body &b = rigid->body[ibody];
    BodyData &bd = bodydata[ibody];
    bd.points.clear();
    bd.rmax = 0.0;

    bd.urel[0] = flow[0] - b.vcm[0];
    bd.urel[1] = flow[1] - b.vcm[1];
    bd.urel[2] = flow[2] - b.vcm[2];
    bd.umag = MathExtra::len3(bd.urel);

    // no relative flow: no drag, and no flow direction to project along

    if (bd.umag == 0.0) continue;

    double uhat[3], dir[3];
    uhat[0] = bd.urel[0] / bd.umag;
    uhat[1] = bd.urel[1] / bd.umag;
    uhat[2] = bd.urel[2] / bd.umag;
    MathExtra::transpose_matvec(b.ex_space, b.ey_space, b.ez_space, uhat, dir);
    plane_basis(dir, bd.e1, bd.e2);
  }

  // 2. one pass over the atoms: project the surface atoms of each body

  collect_projected_points(rigid);

  // 3. per body: shape factors, projected area and drag force

  for (int ibody = 0; ibody < nbody; ibody++) {
    FixRigidAbrade::Body &b = rigid->body[ibody];
    BodyData &bd = bodydata[ibody];
    if (bd.umag == 0.0) continue;

    // shape factors from the body volume and surface area,
    // both kept up to date by fix rigid/abrade after abrasion
    // d = diameter of the sphere with the same volume
    // phi = sphericity, clamped to <= 1: mesh round-off can push it slightly
    //   above 1, and then -log10(phi) < 0 would make the drag formula NaN

    const double d = cbrt(6.0 * b.volume / MY_PI);
    double phi = MY_PI * d * d / b.surface_area;
    if (phi > 1.0) phi = 1.0;

    // projected area: 2d convex hull of the projected atom centres

    double area, perimeter;
    if (!hull_area_perimeter(bd.points, area, perimeter))
      error->one(FLERR, "Fix hydro/drag cannot compute the projected area of rigid body {}",
                 atom->tag[b.ilocal]);

    // atom radius correction: the hull of the atom centres grown outwards
    // by the atom radius r has the exact area A + P*r + pi*r^2
    // (A = hull area, P = hull perimeter; assumes one radius per body)
    // NOTE: always applied, independent of the offset keyword of
    //   fix rigid/abrade; without offset the body volume and surface area
    //   are computed from the atom centres, so this may need to be linked
    //   to offset later

    const double r = bd.rmax;
    const double aperp = area + perimeter * r + MY_PI * r * r;

    // transverse sphericity

    const double phiperp = 0.25 * MY_PI * d * d / aperp;

    // drag force, Holzer & Sommerfeld (2008), simplified form with
    // phi_parallel replaced by phi_perp:
    //   F = 1/2 rho (pi d^2/4) C_D |u| urel
    // with Re = |u| d / nu substituted and multiplied out,
    //   C_D |u| = a1 + a2 sqrt|u| + a3 |u|
    // so nothing is divided by |u| or Re (no 0*inf = NaN as |u| -> 0)
    // divide by ftm2v to convert mass*length/time^2 to force units

    const double nud = nu / d;
    const double a1 = (8.0 / sqrt(phiperp) + 16.0 / sqrt(phi)) * nud;
    const double a2 = 3.0 / pow(phi, 0.75) * sqrt(nud);
    const double a3 = 0.42 * pow(10.0, 0.4 * pow(-log10(phi), 0.2)) / phiperp;

    const double coeff = 0.5 * rho * (0.25 * MY_PI * d * d) *
        (a1 + a2 * sqrt(bd.umag) + a3 * bd.umag) / force->ftm2v;

    // add to the COM force only, not to the atoms, so drag never drives
    // abrasion; no torque

    b.fcm[0] += coeff * bd.urel[0];
    b.fcm[1] += coeff * bd.urel[1];
    b.fcm[2] += coeff * bd.urel[2];
  }
}

/* ----------------------------------------------------------------------
   bucket method: a single pass over the local atoms; each surface atom of
   a body with relative flow is projected onto that body's plane (e1, e2)
   and appended to that body's point list, and the largest atom radius of
   each body is recorded
   uses the atom centres in body coordinates (displace), so periodic
   images do not matter; the COM atom of each body is skipped
------------------------------------------------------------------------- */

void FixHydroDrag::collect_projected_points(FixRigidAbrade *rigid)
{
  const int nlocal = atom->nlocal;
  const int nbody = static_cast<int>(bodydata.size());
  tagint *tag = atom->tag;
  double *radius = atom->radius;

  for (int i = 0; i < nlocal; i++) {
    const int ibody = rigid->atom2body[i];
    if (ibody < 0 || ibody >= nbody) continue;      // not in a body owned by this process
    if (rigid->bodytag[i] == tag[i]) continue;      // COM atom of the body

    BodyData &bd = bodydata[ibody];
    if (bd.umag == 0.0) continue;

    const double *x = rigid->displace[i];
    bd.points.push_back({MathExtra::dot3(x, bd.e1), MathExtra::dot3(x, bd.e2)});
    bd.rmax = std::max(bd.rmax, radius[i]);
  }
}
