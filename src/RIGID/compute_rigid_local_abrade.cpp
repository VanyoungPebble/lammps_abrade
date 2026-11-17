// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_rigid_local_abrade.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "modify.h"
#include "fix_rigid_abrade.h"
#include "math_const.h"
#include "memory.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <vector>

using namespace LAMMPS_NS;

static constexpr int DELTA = 10000;

enum{ID,MOL,MASS,X,Y,Z,XU,YU,ZU,VX,VY,VZ,FX,FY,FZ,IX,IY,IZ,
     TQX,TQY,TQZ,OMEGAX,OMEGAY,OMEGAZ,ANGMOMX,ANGMOMY,ANGMOMZ,
     QUATW,QUATI,QUATJ,QUATK,INERTIAX,INERTIAY,INERTIAZ, VOLUME, SURFACE_AREA, ABRADED_VOLUME, WEAR_ENERGY,
     SPHERICITY_WADELL,L_OBB,I_OBB,S_OBB,OBB_VALID,ZINGG_S_OVER_I_OBB,ZINGG_I_OVER_L_OBB,
     SPHERICITY_KRUMBEIN_OBB,FLATNESS_POTTICARY_OBB,ELONGATION_POTTICARY_OBB,
     FLATNESS_KONG_OBB,ELONGATION_KONG_OBB,FLATNESS_ANGELIDAKIS_OBB,ELONGATION_ANGELIDAKIS_OBB};


/* ---------------------------------------------------------------------- */

ComputeRigidLocalAbrade::ComputeRigidLocalAbrade(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  rstyle(nullptr), idrigid(nullptr), fixrigid(nullptr), vlocal(nullptr), alocal(nullptr)
{
  if (narg < 5) error->all(FLERR,"Illegal compute rigid/local_abrade command");

  local_flag = 1;
  nvalues = narg - 4;

  idrigid = utils::strdup(arg[3]);

  rstyle = new int[nvalues];

  nvalues = 0;
  obb_flag = 0;
  for (int iarg = 4; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"id") == 0) rstyle[nvalues++] = ID;
    else if (strcmp(arg[iarg],"mol") == 0) rstyle[nvalues++] = MOL;
    else if (strcmp(arg[iarg],"mass") == 0) rstyle[nvalues++] = MASS;
    else if (strcmp(arg[iarg],"x") == 0) rstyle[nvalues++] = X;
    else if (strcmp(arg[iarg],"y") == 0) rstyle[nvalues++] = Y;
    else if (strcmp(arg[iarg],"z") == 0) rstyle[nvalues++] = Z;
    else if (strcmp(arg[iarg],"xu") == 0) rstyle[nvalues++] = XU;
    else if (strcmp(arg[iarg],"yu") == 0) rstyle[nvalues++] = YU;
    else if (strcmp(arg[iarg],"zu") == 0) rstyle[nvalues++] = ZU;
    else if (strcmp(arg[iarg],"vx") == 0) rstyle[nvalues++] = VX;
    else if (strcmp(arg[iarg],"vy") == 0) rstyle[nvalues++] = VY;
    else if (strcmp(arg[iarg],"vz") == 0) rstyle[nvalues++] = VZ;
    else if (strcmp(arg[iarg],"fx") == 0) rstyle[nvalues++] = FX;
    else if (strcmp(arg[iarg],"fy") == 0) rstyle[nvalues++] = FY;
    else if (strcmp(arg[iarg],"fz") == 0) rstyle[nvalues++] = FZ;
    else if (strcmp(arg[iarg],"ix") == 0) rstyle[nvalues++] = IX;
    else if (strcmp(arg[iarg],"iy") == 0) rstyle[nvalues++] = IY;
    else if (strcmp(arg[iarg],"iz") == 0) rstyle[nvalues++] = IZ;
    else if (strcmp(arg[iarg],"tqx") == 0) rstyle[nvalues++] = TQX;
    else if (strcmp(arg[iarg],"tqy") == 0) rstyle[nvalues++] = TQY;
    else if (strcmp(arg[iarg],"tqz") == 0) rstyle[nvalues++] = TQZ;
    else if (strcmp(arg[iarg],"omegax") == 0) rstyle[nvalues++] = OMEGAX;
    else if (strcmp(arg[iarg],"omegay") == 0) rstyle[nvalues++] = OMEGAY;
    else if (strcmp(arg[iarg],"omegaz") == 0) rstyle[nvalues++] = OMEGAZ;
    else if (strcmp(arg[iarg],"angmomx") == 0) rstyle[nvalues++] = ANGMOMX;
    else if (strcmp(arg[iarg],"angmomy") == 0) rstyle[nvalues++] = ANGMOMY;
    else if (strcmp(arg[iarg],"angmomz") == 0) rstyle[nvalues++] = ANGMOMZ;
    else if (strcmp(arg[iarg],"quatw") == 0) rstyle[nvalues++] = QUATW;
    else if (strcmp(arg[iarg],"quati") == 0) rstyle[nvalues++] = QUATI;
    else if (strcmp(arg[iarg],"quatj") == 0) rstyle[nvalues++] = QUATJ;
    else if (strcmp(arg[iarg],"quatk") == 0) rstyle[nvalues++] = QUATK;
    else if (strcmp(arg[iarg],"inertiax") == 0) rstyle[nvalues++] = INERTIAX;
    else if (strcmp(arg[iarg],"inertiay") == 0) rstyle[nvalues++] = INERTIAY;
    else if (strcmp(arg[iarg],"inertiaz") == 0) rstyle[nvalues++] = INERTIAZ;
    else if (strcmp(arg[iarg],"volume") == 0) rstyle[nvalues++] = VOLUME;
    else if (strcmp(arg[iarg],"surface_area") == 0) rstyle[nvalues++] = SURFACE_AREA;
    else if (strcmp(arg[iarg],"abraded_volume") == 0) rstyle[nvalues++] = ABRADED_VOLUME;
    else if (strcmp(arg[iarg],"wear_energy") == 0) rstyle[nvalues++] = WEAR_ENERGY;
    else if (strcmp(arg[iarg],"sphericity_wadell") == 0) rstyle[nvalues++] = SPHERICITY_WADELL;
    else if (strcmp(arg[iarg],"l_obb") == 0) rstyle[nvalues++] = L_OBB;
    else if (strcmp(arg[iarg],"i_obb") == 0) rstyle[nvalues++] = I_OBB;
    else if (strcmp(arg[iarg],"s_obb") == 0) rstyle[nvalues++] = S_OBB;
    else if (strcmp(arg[iarg],"obb_valid") == 0) rstyle[nvalues++] = OBB_VALID;
    else if (strcmp(arg[iarg],"zingg_s_over_i_obb") == 0) rstyle[nvalues++] = ZINGG_S_OVER_I_OBB;
    else if (strcmp(arg[iarg],"zingg_i_over_l_obb") == 0) rstyle[nvalues++] = ZINGG_I_OVER_L_OBB;
    else if (strcmp(arg[iarg],"sphericity_krumbein_obb") == 0) rstyle[nvalues++] = SPHERICITY_KRUMBEIN_OBB;
    else if (strcmp(arg[iarg],"flatness_potticary_obb") == 0) rstyle[nvalues++] = FLATNESS_POTTICARY_OBB;
    else if (strcmp(arg[iarg],"elongation_potticary_obb") == 0) rstyle[nvalues++] = ELONGATION_POTTICARY_OBB;
    else if (strcmp(arg[iarg],"flatness_kong_obb") == 0) rstyle[nvalues++] = FLATNESS_KONG_OBB;
    else if (strcmp(arg[iarg],"elongation_kong_obb") == 0) rstyle[nvalues++] = ELONGATION_KONG_OBB;
    else if (strcmp(arg[iarg],"flatness_angelidakis_obb") == 0) rstyle[nvalues++] = FLATNESS_ANGELIDAKIS_OBB;
    else if (strcmp(arg[iarg],"elongation_angelidakis_obb") == 0) rstyle[nvalues++] = ELONGATION_ANGELIDAKIS_OBB;
    else error->all(FLERR,"Invalid keyword in compute rigid/local_abrade_abrade command");

    if (rstyle[nvalues-1] >= L_OBB) obb_flag = 1;
  }

  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  ncount = nmax = 0;
  vlocal = nullptr;
  alocal = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeRigidLocalAbrade::~ComputeRigidLocalAbrade()
{
  memory->destroy(vlocal);
  memory->destroy(alocal);
  delete[] idrigid;
  delete[] rstyle;
}

/* ---------------------------------------------------------------------- */

void ComputeRigidLocalAbrade::init()
{
  // set fixrigid

  auto *ifix = modify->get_fix_by_id(idrigid);
  if (!ifix) error->all(FLERR,"FixRigidAbrade ID {} for compute rigid/local_abrade does not exist", idrigid);
  fixrigid = dynamic_cast<FixRigidAbrade *>(ifix);
  if (!fixrigid)
    error->all(FLERR,"Fix ID {} for compute rigid/local_abrade does not point to fix rigid/abrade", idrigid);

  // do initial memory allocation so that memory_usage() is correct

  ncount = compute_rigid(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputeRigidLocalAbrade::compute_local()
{
  invoked_local = update->ntimestep;

  // All ranks participate, including ranks with no owned bodies.  Do this
  // once before the count/fill passes, and only for requested OBB attributes.
  if (obb_flag) compute_obb();

  // count local entries and compute bond info

  ncount = compute_rigid(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  ncount = compute_rigid(1);
}

/* ----------------------------------------------------------------------
   Read the current physical surface without changing the rigid integrator.
   Owned vertices are unique; body coordinates avoid periodic image jumps.
   The first version gathers vertices on output, without persistent caches.
------------------------------------------------------------------------- */

void ComputeRigidLocalAbrade::compute_obb()
{
  std::vector<tagint> local_ids;
  std::vector<double> local_xyz;
  const double nan = std::numeric_limits<double>::quiet_NaN();

  for (int i = 0; i < atom->nlocal; i++) {
    const int ibody = fixrigid->atom2body[i];
    if (ibody < 0 || fixrigid->bodytag[i] == atom->tag[i]) continue;

    double p[3] = {fixrigid->displace[i][0], fixrigid->displace[i][1],
                   fixrigid->displace[i][2]};
    if (fixrigid->offset_flag) {
      // Integrator output stores inward-offset contact sphere centers.
      // Restore the same outward surface used for volume/area, in a copy.
      const double r = std::hypot(std::hypot(p[0], p[1]), p[2]);
      const double scale = r > 0.0 ? 1.0 + atom->radius[i]/r : nan;
      for (double &value : p) value *= scale;
    }
    local_ids.push_back(fixrigid->bodytag[i]);
    local_ids.push_back(atom->tag[i]);
    local_xyz.insert(local_xyz.end(), p, p+3);
  }

  int nprocs;
  MPI_Comm_size(world, &nprocs);
  int overflow = local_xyz.size() > static_cast<size_t>(std::numeric_limits<int>::max());
  int any_overflow;
  MPI_Allreduce(&overflow, &any_overflow, 1, MPI_INT, MPI_MAX, world);
  if (any_overflow)
    error->all(FLERR,"Too many surface vertices for compute rigid/local_abrade OBB");

  int nlocal = static_cast<int>(local_ids.size()/2);
  std::vector<int> counts(nprocs), idcounts(nprocs), xyzcounts(nprocs);
  std::vector<int> idoffsets(nprocs), xyzoffsets(nprocs);
  MPI_Allgather(&nlocal, 1, MPI_INT, counts.data(), 1, MPI_INT, world);
  long long total = 0;
  for (int p = 0; p < nprocs; p++) {
    total += counts[p];
    if (total > std::numeric_limits<int>::max()/3)
      error->all(FLERR,"Too many surface vertices for compute rigid/local_abrade OBB");
    idcounts[p] = 2*counts[p];
    xyzcounts[p] = 3*counts[p];
    idoffsets[p] = 2*static_cast<int>(total-counts[p]);
    xyzoffsets[p] = 3*static_cast<int>(total-counts[p]);
  }

  std::vector<tagint> ids(2*total);
  std::vector<double> xyz(3*total);
  MPI_Allgatherv(local_ids.data(), 2*nlocal, MPI_LMP_TAGINT, ids.data(),
                 idcounts.data(), idoffsets.data(), MPI_LMP_TAGINT, world);
  MPI_Allgatherv(local_xyz.data(), 3*nlocal, MPI_DOUBLE, xyz.data(),
                 xyzcounts.data(), xyzoffsets.data(), MPI_DOUBLE, world);

  // Index only bodies whose owner is selected locally.  Stable atom tags
  // order and deduplicate each surface independently of MPI decomposition.
  std::map<tagint, int> owned;
  std::vector<std::map<tagint, ShapeOBB::Point>> vertices(fixrigid->nlocal_body);
  obb_results.assign(fixrigid->nlocal_body, ShapeOBB::Result{});
  for (int b = 0; b < fixrigid->nlocal_body; b++) {
    const int i = fixrigid->body[b].ilocal;
    if (atom->mask[i] & groupbit) owned.emplace(atom->tag[i], b);
  }
  for (long long j = 0; j < total; j++) {
    const auto body = owned.find(ids[2*j]);
    if (body == owned.end()) continue;
    vertices[body->second].emplace(ids[2*j+1],
        ShapeOBB::Point{{xyz[3*j], xyz[3*j+1], xyz[3*j+2]}});
  }
  for (const auto &entry : owned) {
    const int b = entry.second;
    std::vector<ShapeOBB::Point> points;
    points.reserve(vertices[b].size());
    for (const auto &vertex : vertices[b]) points.push_back(vertex.second);
    obb_results[b] = ShapeOBB::min_volume(points);
  }
}

/* ----------------------------------------------------------------------
   count rigid bodies and compute rigid info on this proc
   if flag is set, compute requested info about rigid body
   owning atom of rigid body must be in group
------------------------------------------------------------------------- */

int ComputeRigidLocalAbrade::compute_rigid(int flag)
{
  int i,m,n,ibody;
  double *ptr;
  FixRigidAbrade::Body *body;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  m = 0;
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    ibody = fixrigid->bodyown[i];
    if (ibody < 0) continue;
    body = &fixrigid->body[ibody];

    if (flag) {
      if (nvalues == 1) ptr = &vlocal[m];
      else ptr = alocal[m];

      const double nan = std::numeric_limits<double>::quiet_NaN();
      const bool valid = obb_flag && obb_results[ibody].valid;
      const double S = valid ? obb_results[ibody].length[0] : nan;
      const double I = valid ? obb_results[ibody].length[1] : nan;
      const double L = valid ? obb_results[ibody].length[2] : nan;
      const double si = valid ? S/L : nan;
      const double ii = valid ? I/L : nan;

      for (n = 0; n < nvalues; n++) {
        switch (rstyle[n]) {
        case ID:
          ptr[n] = tag[body->ilocal];
          break;
        case MOL:
          ptr[n] = molecule[body->ilocal];
          break;
        case MASS:
          ptr[n] = body->mass;
          break;
        case X:
          ptr[n] = body->xcm[0];
          break;
        case Y:
          ptr[n] = body->xcm[1];
          break;
        case Z:
          ptr[n] = body->xcm[2];
          break;
        case XU:
          ptr[n] = body->xcm[0] +
            ((body->image & IMGMASK) - IMGMAX) * xprd;
          break;
        case YU:
          ptr[n] = body->xcm[1] +
            ((body->image >> IMGBITS & IMGMASK) - IMGMAX) * yprd;
          break;
        case ZU:
          ptr[n] = body->xcm[2] +
            ((body->image >> IMG2BITS) - IMGMAX) * zprd;
          break;
        case VX:
          ptr[n] = body->vcm[0];
          break;
        case VY:
          ptr[n] = body->vcm[1];
          break;
        case VZ:
          ptr[n] = body->vcm[2];
          break;
        case FX:
          ptr[n] = body->fcm[0];
          break;
        case FY:
          ptr[n] = body->fcm[1];
          break;
        case FZ:
          ptr[n] = body->fcm[2];
          break;
        case IX:
          ptr[n] = (body->image & IMGMASK) - IMGMAX;
          break;
        case IY:
          ptr[n] = (body->image >> IMGBITS & IMGMASK) - IMGMAX;
          break;
        case IZ:
          ptr[n] = (body->image >> IMG2BITS) - IMGMAX;
          break;
        case TQX:
          ptr[n] = body->torque[0];
          break;
        case TQY:
          ptr[n] = body->torque[1];
          break;
        case TQZ:
          ptr[n] = body->torque[2];
          break;
        case OMEGAX:
          ptr[n] = body->omega[0];
          break;
        case OMEGAY:
          ptr[n] = body->omega[1];
          break;
        case OMEGAZ:
          ptr[n] = body->omega[2];
          break;
        case ANGMOMX:
          ptr[n] = body->angmom[0];
          break;
        case ANGMOMY:
          ptr[n] = body->angmom[1];
          break;
        case ANGMOMZ:
          ptr[n] = body->angmom[2];
          break;
        case QUATW:
          ptr[n] = body->quat[0];
          break;
        case QUATI:
          ptr[n] = body->quat[1];
          break;
        case QUATJ:
          ptr[n] = body->quat[2];
          break;
        case QUATK:
          ptr[n] = body->quat[3];
          break;
        case INERTIAX:
          ptr[n] = body->inertia[0];
          break;
        case INERTIAY:
          ptr[n] = body->inertia[1];
          break;
        case INERTIAZ:
          ptr[n] = body->inertia[2];
          break;
        case VOLUME:
          ptr[n] = body->volume;
          break;
        case SURFACE_AREA:
          ptr[n] = body->surface_area;
          break;
        case ABRADED_VOLUME:
          ptr[n] = body->abraded_volume;
          break;
        case WEAR_ENERGY:
          ptr[n] = body->wear_energy;
          break;
        case SPHERICITY_WADELL:
          if (std::isfinite(body->volume) && body->volume > 0.0 &&
              std::isfinite(body->surface_area) && body->surface_area > 0.0) {
            const double r = std::cbrt(body->volume)/std::sqrt(body->surface_area);
            ptr[n] = std::cbrt(36.0*MathConst::MY_PI)*r*r;
          } else ptr[n] = nan;
          break;
        case L_OBB: ptr[n] = L; break;
        case I_OBB: ptr[n] = I; break;
        case S_OBB: ptr[n] = S; break;
        case OBB_VALID: ptr[n] = valid ? 1.0 : 0.0; break;
        case ZINGG_S_OVER_I_OBB: ptr[n] = S/I; break;
        case ZINGG_I_OVER_L_OBB: ptr[n] = ii; break;
        case SPHERICITY_KRUMBEIN_OBB: ptr[n] = std::cbrt(ii*si); break;
        case FLATNESS_POTTICARY_OBB: ptr[n] = 2.0*(ii-si)/(1.0+ii+si); break;
        case ELONGATION_POTTICARY_OBB: ptr[n] = (1.0-ii)/(1.0+ii+si); break;
        case FLATNESS_KONG_OBB: ptr[n] = 1.0-S/I; break;
        case ELONGATION_KONG_OBB: ptr[n] = 1.0-ii; break;
        case FLATNESS_ANGELIDAKIS_OBB: ptr[n] = ii*ii/(si+ii*ii)-si/(1.0+si); break;
        case ELONGATION_ANGELIDAKIS_OBB: ptr[n] = si/(si+ii*ii)-si/(1.0+si); break;
        }
      }
    }

    m++;
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRigidLocalAbrade::reallocate(int n)
{
  // grow vector_local or array_local

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vlocal);
    memory->create(vlocal,nmax,"rigid/local_abrade:vector_local");
    vector_local = vlocal;
  } else {
    memory->destroy(alocal);
    memory->create(alocal,nmax,nvalues,"rigid/local_abrade:array_local");
    array_local = alocal;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeRigidLocalAbrade::memory_usage()
{
  double bytes = (double)nmax*nvalues * sizeof(double);
  bytes += obb_results.capacity()*sizeof(ShapeOBB::Result);
  return bytes;
}
