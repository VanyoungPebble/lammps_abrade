/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Distributed under the GNU General Public License; see LICENSE.
------------------------------------------------------------------------- */

#include "compute_projected_area.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "update.h"

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullError.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

using namespace LAMMPS_NS;

namespace {
using Point = std::array<double, 2>;

struct Record {
  tagint molecule;
  Point point;
};

// Qhull implements the hull and area calculation.  This wrapper only removes
// exact duplicates and applies a uniform scale to keep units/translation benign.
double projected_hull_area(std::vector<Point> &points)
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  points.erase(std::unique(points.begin(), points.end()), points.end());
  if (points.size() < 3) return nan;

  Point low = points.front(), high = low;
  for (const auto &p : points)
    for (int k = 0; k < 2; k++) {
      low[k] = std::min(low[k], p[k]);
      high[k] = std::max(high[k], p[k]);
    }
  const double scale = std::max(high[0] - low[0], high[1] - low[1]);
  if (!std::isfinite(scale) || scale <= 0.0) return nan;

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
    hull.runQhull("projected/area", 2, static_cast<int>(points.size()),
                  coordinates.data(), "Qt");
    // In two dimensions Qhull volume() is enclosed area; area() is perimeter.
    // No QJ/random perturbation and no alternate hull implementation are used.
    const double area = (hull.volume() * scale) * scale;
    return std::isfinite(area) && area > 0.0 ? area : nan;
  } catch (const orgQhull::QhullError &) {
    return nan;
  }
}
}    // namespace

ComputeProjectedArea::ComputeProjectedArea(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), nmax(0), warned(0)
{
  if (narg != 6) error->all(FLERR, "Illegal compute projected/area command");

  double direction[3];
  for (int k = 0; k < 3; k++) {
    direction[k] = utils::numeric(FLERR, arg[3 + k], false, lmp);
    if (!std::isfinite(direction[k]))
      error->all(FLERR, "Compute projected/area direction must be finite and nonzero");
  }
  const double scale = std::max({std::abs(direction[0]), std::abs(direction[1]),
                                 std::abs(direction[2])});
  if (scale == 0.0)
    error->all(FLERR, "Compute projected/area direction must be finite and nonzero");
  for (double &value : direction) value /= scale;
  const double norm = std::hypot(std::hypot(direction[0], direction[1]), direction[2]);
  for (double &value : direction) value /= norm;

  // Choose the least parallel coordinate axis to construct a stable plane basis.
  int axis = 0;
  for (int k = 1; k < 3; k++)
    if (std::abs(direction[k]) < std::abs(direction[axis])) axis = k;
  double helper[3] = {0.0, 0.0, 0.0};
  helper[axis] = 1.0;
  u[0] = direction[1] * helper[2] - direction[2] * helper[1];
  u[1] = direction[2] * helper[0] - direction[0] * helper[2];
  u[2] = direction[0] * helper[1] - direction[1] * helper[0];
  const double unorm = std::hypot(std::hypot(u[0], u[1]), u[2]);
  for (double &value : u) value /= unorm;
  v[0] = direction[1] * u[2] - direction[2] * u[1];
  v[1] = direction[2] * u[0] - direction[0] * u[2];
  v[2] = direction[0] * u[1] - direction[1] * u[0];

  local_flag = 1;
  size_local_cols = 2;
  size_local_rows = 0;
  array_local = nullptr;
}

ComputeProjectedArea::~ComputeProjectedArea()
{
  memory->destroy(array_local);
}

void ComputeProjectedArea::init()
{
  if (!atom->molecule_flag)
    error->all(FLERR, "Compute projected/area requires atom attribute molecule");
}

void ComputeProjectedArea::compute_local()
{
  invoked_local = update->ntimestep;
  size_local_rows = 0;
  std::vector<tagint> local_ids;
  std::vector<double> local_xy;
  int invalid[2] = {0, 0};

  // Read current positions directly, independently of fix rigid/abrade and OBB.
  // Owned atoms only: MPI ghost copies never enter the point cloud.
  for (int i = 0; i < atom->nlocal; i++) {
    if (!(atom->mask[i] & groupbit)) continue;
    if (atom->molecule[i] <= 0) invalid[0] = 1;
    double position[3], xy[2] = {0.0, 0.0};
    domain->unmap(atom->x[i], atom->image[i], position);
    for (int k = 0; k < 3; k++) {
      if (!std::isfinite(position[k])) invalid[1] = 1;
      xy[0] += position[k] * u[k];
      xy[1] += position[k] * v[k];
    }
    if (!std::isfinite(xy[0]) || !std::isfinite(xy[1])) invalid[1] = 1;
    local_ids.push_back(atom->molecule[i]);
    local_xy.push_back(xy[0]);
    local_xy.push_back(xy[1]);
  }

  int any_invalid[2];
  MPI_Allreduce(invalid, any_invalid, 2, MPI_INT, MPI_MAX, world);
  if (any_invalid[0])
    error->all(FLERR, "Compute projected/area requires positive molecule IDs in its group");
  if (any_invalid[1])
    error->all(FLERR, "Compute projected/area encountered a non-finite position or projection");

  const long long local_count = static_cast<long long>(local_ids.size());
  long long total_count;
  MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG_LONG, MPI_SUM, world);
  if (total_count > std::numeric_limits<int>::max() / 2)
    error->all(FLERR, "Too many points for compute projected/area MPI gather");
  const int nlocal = static_cast<int>(local_count);
  std::vector<int> counts(comm->nprocs), offsets(comm->nprocs);
  std::vector<int> xycounts(comm->nprocs), xyoffsets(comm->nprocs);
  MPI_Gather(&nlocal, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, world);
  if (comm->me == 0) {
    int offset = 0;
    for (int p = 0; p < comm->nprocs; p++) {
      offsets[p] = offset;
      xycounts[p] = 2 * counts[p];
      xyoffsets[p] = 2 * offset;
      offset += counts[p];
    }
  }

  const int total = comm->me == 0 ? static_cast<int>(total_count) : 0;
  std::vector<tagint> ids(total);
  std::vector<double> xy(2 * total);
  MPI_Gatherv(local_ids.data(), nlocal, MPI_LMP_TAGINT, ids.data(), counts.data(),
              offsets.data(), MPI_LMP_TAGINT, 0, world);
  MPI_Gatherv(local_xy.data(), 2 * nlocal, MPI_DOUBLE, xy.data(), xycounts.data(),
              xyoffsets.data(), MPI_DOUBLE, 0, world);
  if (comm->me != 0) return;

  // A prototype root-only local table keeps MPI ownership out of the integrator.
  std::vector<Record> records;
  records.reserve(total);
  for (int i = 0; i < total; i++) records.push_back({ids[i], {{xy[2 * i], xy[2 * i + 1]}}});
  std::sort(records.begin(), records.end(), [](const Record &a, const Record &b) {
    return a.molecule != b.molecule ? a.molecule < b.molecule : a.point < b.point;
  });

  std::vector<std::array<double, 2>> rows;
  for (size_t begin = 0; begin < records.size();) {
    size_t end = begin + 1;
    while (end < records.size() && records[end].molecule == records[begin].molecule) end++;
    std::vector<Point> points;
    points.reserve(end - begin);
    for (size_t i = begin; i < end; i++) points.push_back(records[i].point);
    const double area = projected_hull_area(points);
    if (!std::isfinite(area) && !warned) {
      error->warning(FLERR, "Compute projected/area cannot form a 2D hull for molecule {}; "
                     "writing NaN (further hull warnings suppressed)", records[begin].molecule);
      warned = 1;
    }
    rows.push_back({{static_cast<double>(records[begin].molecule), area}});
    begin = end;
  }
  size_local_rows = static_cast<int>(rows.size());
  if (size_local_rows > nmax) {
    nmax = size_local_rows;
    memory->grow(array_local, nmax, 2, "projected/area:array_local");
  }
  for (int i = 0; i < size_local_rows; i++) {
    array_local[i][0] = rows[i][0];
    array_local[i][1] = rows[i][1];
  }
}

double ComputeProjectedArea::memory_usage()
{
  return static_cast<double>(nmax) * 2 * sizeof(double);
}
