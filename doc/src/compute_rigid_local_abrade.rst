.. index:: compute rigid/local_abrade

compute rigid/local_abrade command
==================================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID rigid/local_abrade rigidID input1 input2 ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* rigid/local_abrade = style name of this compute command
* rigidID = ID of fix rigid/abrade command
* input = one or more rigid body attributes

  .. parsed-literal::

       possible attributes = id, mol, mass, volume,
                             surface_area, abraded_volume, wear_energy,
                             x, y, z, xu, yu, zu, ix, iy, iz,
                             vx, vy, vz, fx, fy, fz,
                             omegax, omegay, omegaz,
                             angmomx, angmomy, angmomz,
                             quatw, quati, quatj, quatk,
                             tqx, tqy, tqz,
                             inertiax, inertiay, inertiaz,
                             sphericity_wadell, l_obb, i_obb, s_obb, obb_valid,
                             zingg_s_over_i_obb, zingg_i_over_l_obb,
                             sphericity_krumbein_obb,
                             flatness_potticary_obb, elongation_potticary_obb,
                             flatness_kong_obb, elongation_kong_obb,
                             flatness_angelidakis_obb, elongation_angelidakis_obb

  .. parsed-literal::

           id = atom ID of atom within body which owns body properties
           mol = molecule ID used to define body in :doc:`fix rigid/small <fix_rigid>` command
           volume = total volume of body
           mass = total mass of body
           surface_area = total surface area of the body
           abraded_volume = cumulative volume removed from the particle 
           wear_energy = cumulative dissipated energy
           x,y,z = center of mass coords of body
           xu,yu,zu = unwrapped center of mass coords of body
           ix,iy,iz = box image that the center of mass is in
           vx,vy,vz = center of mass velocities
           fx,fy,fz = force of center of mass
           omegax,omegay,omegaz = angular velocity of body
           angmomx,angmomy,angmomz = angular momentum of body
           quatw,quati,quatj,quatk = quaternion components for body
           tqx,tqy,tqz = torque on body
           inertiax,inertiay,inertiaz = diagonalized moments of inertia of body

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all rigid/local_abrade abrade id volume surface_area abraded_volume wear_energy

.. code-block:: LAMMPS

   compute shape all rigid/local_abrade abrade mol volume surface_area l_obb i_obb s_obb sphericity_wadell obb_valid
   dump shape_out all local 40000 particle_shape.dump c_shape[*]
   dump_modify shape_out sort 1 format float %.17g

Description
"""""""""""


An extension of the :doc:`rigid/local <compute_rigid_local>` command, *compute rigid/local_abrade* grants 
access to abrasion-specific rigid body attributes for rigid bodies defined by the :doc:`rigid/abrade <fix_rigid_abrade>`
command. The data is stored as local data so it can be accessed by other :doc:`output commands <Howto_output>` that 
process local data, such as the :doc:`compute reduce <compute_reduce>` or :doc:`dump local <dump>` commands.

Here is an example of how to use this compute to dump rigid body info to a file:

.. code-block:: LAMMPS

   compute 1 all rigid/local_abrade abrade id mass volume surface_area abraded_volume wear_energy
   dump 2 all local 500 dump.abrade index c_3[1] c_3[2] c_3[3] c_3[4] c_3[5] c_3[6] 

----------

This section explains the abrasion-specific rigid body attributes that can be specified.

The *volume* attribute is the total volume of the body as calculated through the triangulation stored in the :doc:`angles <angles>` comprising the surface atoms.

The *mass* is the product of the body's volume and mass density. 

The *surface_area* is the sum of associated areas for all atoms in the body.

The *abraded_volume* is the total abraded volume removed from the body. This cumulative value can be preserved between simulations through the use of the optional *infile* keyword in the :doc:`rigid/abrade <fix_rigid_abrade>` command.

The *wear_energy* is the total energy dissipated to displace surface atoms inwards along their normals. The energy dissipated to displace a surface atom :math:`i` is adapted from :ref:`(Capozza) <Capozza_compute>` as,

.. math::

   E_{w,i} = \left| \mathbf{F}_{i,n} \right| \left| \mathbf{u}_i \right| \Delta t

where, :math:`\mathbf{F}_{i,n}` is the imposed force aligned along the atom's associate normal, :math:`\mathbf{u}_i` is the total abrasion velocity, and :math:`\Delta t` is the employed timestep.

----------

The *l_obb*, *i_obb*, and *s_obb* attributes are the full long, intermediate,
and short dimensions of an oriented bounding box (OBB), respectively, with
:math:`0 < S \le I \le L`.  These dimensions use the ``minVolume``, level-3
search used by SHAPE :ref:`(Angelidakis) <Angelidakis_shape>`.  The search
tests convex-hull face orientations and the three frame permutations for
each convex-hull edge, with a minimum-area rectangle in each projection.
It is a bounding-box search, without a guarantee of the global minimum.
SHAPE's level-3 edge-frame candidates depend on the coordinate axes; to
reproduce these results offline, apply the reference algorithm to the same
surface in the same body principal frame, not to untransformed global dump
coordinates.  Different triangulations of coplanar convex-hull facets or
equal-volume box solutions can also change the selected dimensions.
It is distinct from PCA dimensions, inertia-equivalent ellipsoid axes, and
surface-fitted ellipsoid axes.

Only owned surface vertices are gathered; the center owning atom is
excluded.  Coordinates are evaluated in each body's principal frame, so
periodic images do not split a particle.  With the *offset* option, an
outward surface is reconstructed in a temporary copy to match the surface
used for the body's volume and area.  The fix and its coordinates are not
modified.  Output after a time step uses the completed abrasion and
remeshing state.

The *obb_valid* attribute is 1 for a successful three-dimensional box and
0 for an invalid or numerically degenerate point set.  Invalid box
dimensions and their derived indices are NaN, not zero or a substitute
axis definition.  The algorithm centers and scales the vertices internally;
dimensions below its relative geometric tolerance can be classified as
degenerate.  Nearly symmetric particles can have nonunique box directions.
Only the sorted dimensions, not directions, are returned.

The *sphericity_wadell* attribute uses the current body volume :math:`V`
and surface area :math:`A`:

.. math::

   \psi_W = \frac{\pi^{1/3}(6V)^{2/3}}{A}.

It does not require an OBB calculation and is NaN for nonpositive or
nonfinite volume or area.  Values are not clipped to the interval [0,1].
As with the existing volume and area outputs, it describes the triangulated
surface rather than an ideal analytic curved surface.

The remaining shape attributes use the OBB dimensions and the same
formulas as SHAPE.  Their names distinguish the different definitions:

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Attribute
     - Formula
   * - ``zingg_s_over_i_obb``
     - :math:`S/I`
   * - ``zingg_i_over_l_obb``
     - :math:`I/L`
   * - ``sphericity_krumbein_obb``
     - :math:`(IS/L^2)^{1/3}`
   * - ``flatness_potticary_obb``
     - :math:`2(I-S)/(L+I+S)`
   * - ``elongation_potticary_obb``
     - :math:`(L-I)/(L+I+S)`
   * - ``flatness_kong_obb``
     - :math:`(I-S)/I`
   * - ``elongation_kong_obb``
     - :math:`(L-I)/L`
   * - ``flatness_angelidakis_obb``
     - :math:`I^2/(LS+I^2)-S/(L+S)`
   * - ``elongation_angelidakis_obb``
     - :math:`LS/(LS+I^2)-S/(L+S)`

OBB-related fields trigger a collective vertex gather once per invocation
of this compute, followed by a box search on each selected body's owner.
No persistent geometry cache is used.  The first implementation replicates
the gathered coordinates on every MPI rank, so its temporary memory and
communication scale with the total number of surface vertices in the fix,
even when this compute selects fewer bodies.  Choose the output interval
to suit the particle count and mesh resolution.  Using only existing
attributes, or adding only *sphericity_wadell*, does not gather vertices or
run the OBB search.  Use *mol* or *id* to identify particles because row order
can change when bodies migrate between processors.

----------

Output info
"""""""""""

This compute calculates a local vector or local array depending on the
number of keywords.  The length of the vector or number of rows in the
array is the number of rigid bodies.  If a single keyword is
specified, a local vector is produced.  If two or more keywords are
specified, a local array is produced where the number of columns = the
number of keywords.  The vector or array can be accessed by any
command that uses local values from a compute as input.  See the
:doc:`Howto output <Howto_output>` page for an overview of LAMMPS
output options.

The vector or array values will be in whatever :doc:`units <units>` the
corresponding attribute is in:

* id,mol = unitless
* mass = mass units
* volume = volume units
* surface_area = area units
* abraded_volume = volume units
* wear_energy = energy units
* x,y,z and xy,yu,zu = distance units
* vx,vy,vz = velocity units
* fx,fy,fz = force units
* omegax,omegay,omegaz = radians/time units
* angmomx,angmomy,angmomz = mass\*distance\ :math:`^2`\ /time units
* quatw,quati,quatj,quatk = unitless
* tqx,tqy,tqz = torque units
* inertiax,inertiay,inertiaz = mass\*distance\ :math:`^2` units
* l_obb,i_obb,s_obb = distance units
* sphericity_wadell,obb_valid and all OBB-derived shape indices = unitless

Restrictions
""""""""""""

This compute is part of the RIGID package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`rigid/local <compute_rigid_local>`, :doc:`dump local <dump>`, :doc:`compute reduce <compute_reduce>`

Default
"""""""

none

----------

.. _Capozza_compute:

**(Capozza)** R. Capozza, K. J. Hanley, A comprehensive model of plastic wear based on the discrete element method, Powder Technol, 410, 117864 (2022).

.. _Angelidakis_shape:

**(Angelidakis)** V. Angelidakis, S. Nadimi, S. Utili, SHape Analyser for
Particle Engineering (SHAPE): Seamless characterisation and simplification
of particle morphology from imaging data, Computer Physics Communications,
265, 107983 (2021).
