.. index:: compute projected/area

compute projected/area command
==============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID projected/area ex ey ez

* ID, group-ID are documented in :doc:`compute <compute>` command
* ex, ey, ez = constant components of the projection direction (dimensionless)

Examples
""""""""

.. code-block:: LAMMPS

   compute ap owned projected/area 1 0 0
   dump projected all local 100 projected_area.dump c_ap[1] c_ap[2]
   dump_modify projected sort 1 format float %.17g

Description
"""""""""""

Project the current unwrapped atom centers in the group onto a plane normal
to the specified direction.  Group the two-dimensional points by positive
molecule ID, and use the Qhull library to compute each point set's convex
hull area.  The direction is normalized internally and must be finite and
nonzero.  It is fixed in space, not attached to the rotating particle.

This prototype reads atom positions and image flags directly.  It does not
access rigid-body displacement arrays, OBB calculations, triangle topology,
or stored surface normals.  Atom radii are not included.  For rigid/abrade,
select all surface atoms in each desired molecule, for example using the
existing ``owned`` group.  With ``offset`` enabled, the result describes the
inward-offset atom centers, not the restored triangular surface.

For a convex triangulated particle whose vertices are the selected atom
centers, this is its orthogonal projected area.  For a non-convex particle,
this is the projected **convex hull** area and can exceed its silhouette
area.  Separate molecules are processed independently, so inter-particle
shadowing is not included.

Only owned atoms, never MPI ghosts, are collected.  Projected points are
gathered to MPI rank 0, sorted by molecule ID and coordinates, and exact
duplicate coordinates are removed per molecule before calling Qhull.
All result rows reside on rank 0; the other ranks have zero local rows.
This simple gather is intended for initial trials rather than very large
particle counts.  Image flags must describe a consistent unwrapped image
of each molecule, including particles crossing periodic boundaries.

The compute runs only when requested, such as on dump timesteps.  It has
no forces, no restart state, and no feedback into abrasion or integration.
Current atom coordinates are read each time, so changes caused by rotation,
abrasion, or atom deletion are reflected on the next request.  The direction
does not automatically follow fluid velocity or particle-relative velocity.

Output info
"""""""""""

This compute produces a local array with two columns: molecule ID and
projected convex hull area (distance units squared).  A degenerate point
set or Qhull failure produces NaN; the first such case per compute issues
a warning.  Qhull's random perturbation option is not enabled and no
alternate convex hull algorithm is substituted.  Invalid positions or
non-positive molecule IDs are errors.  An empty group produces zero rows.

Restrictions
""""""""""""

This prototype is in the RIGID package and requires an atom style providing
molecule IDs.  Use CMake and an installed Qhull CMake package with the
``Qhull::qhullcpp`` and ``Qhull::qhullstatic_r`` targets.  This branch adds
these dependencies when RIGID is enabled.  Traditional make builds have
not been integrated.  Tested dependency: Qhull 2020.2 (8.0.2), official
C++ wrapper over the reentrant library.

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`compute rigid/local_abrade <compute_rigid_local_abrade>`

Default
"""""""

none
