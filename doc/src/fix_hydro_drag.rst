.. index:: fix hydro/drag

fix hydro/drag command
======================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID hydro/drag rho rho_f nu nu_f flow ux uy uz keyword value

* ID, group-ID are documented in :doc:`fix <fix>` command (the group is not used)
* hydro/drag = style name of this fix command
* rho rho_f = fluid mass density (mass/volume units)
* nu nu_f = fluid kinematic viscosity (length^2/time units)
* flow ux uy uz = uniform, steady fluid velocity (velocity units)
* zero or one keyword/value pair may be appended

  .. parsed-literal::

       *buoyancy* value = *yes* or *no*
         yes = also add the buoyancy force of the fully submerged body
         no = drag only

Examples
""""""""

.. code-block:: LAMMPS

   fix water all hydro/drag rho 1.0 nu 0.01 flow 50.0 0.0 0.0
   fix abrade all rigid/abrade ${H} 0.25 ${density} molecule normals drag water

   fix grav none gravity 981.0 vector 0 0 -1
   fix water all hydro/drag rho 1.0 nu 0.01 flow 50.0 0.0 0.0 buoyancy yes
   fix abrade all rigid/abrade ${H} 0.25 ${density} molecule normals gravity grav drag water

Description
"""""""""""

Add a hydrodynamic drag force, and optionally a buoyancy force, to each
body of a :doc:`fix rigid/abrade <fix_rigid_abrade>` command.  This fix
only holds the fluid parameters and computes the forces.  It does not act
on its own: the *drag* keyword of fix rigid/abrade must name this fix.
Fix rigid/abrade then calls it every timestep, after it has summed the
atom forces and added Langevin and gravity forces, and the forces are added
to the center-of-mass force of each body.  They are never added to the
per-atom forces, so they do not drive abrasion.  No torque is applied.

The drag follows the correlation of Holzer and Sommerfeld
:ref:`(Holzer) <Holzer>` in its simplified form, with the lengthwise
sphericity replaced by the crosswise sphericity:

.. math::

   C_D = \frac{8}{Re\sqrt{\Phi_\perp}} + \frac{16}{Re\sqrt{\Phi}}
       + \frac{3}{\sqrt{Re}\,\Phi^{3/4}}
       + \frac{0.42 \cdot 10^{0.4(-\log_{10}\Phi)^{0.2}}}{\Phi_\perp}

with :math:`Re = |u| d/\nu`, the relative velocity
:math:`\mathbf{u} = \mathbf{u}_f - \mathbf{v}_{cm}` of the fluid with respect
to the body center of mass, the volume-equivalent sphere diameter
:math:`d = (6V/\pi)^{1/3}`, the sphericity :math:`\Phi = \pi d^2/S` and the
crosswise sphericity :math:`\Phi_\perp = (\pi d^2/4)/A_\perp`.  The volume
*V* and surface area *S* are those computed by fix rigid/abrade, and are
updated when a body abrades.  :math:`\Phi` is clamped to at most 1.  The force

.. math::

   \mathbf{F}_D = \tfrac{1}{2}\rho_f \frac{\pi d^2}{4} C_D |u|\, \mathbf{u}
               = \tfrac{1}{2}\rho_f \frac{\pi d^2}{4}
                 \left[a_1 + a_2\sqrt{|u|} + a_3|u|\right] \mathbf{u}

is evaluated in the second form, where *Re* has been multiplied out, so
nothing is divided by :math:`|u|` and the force goes to zero smoothly as the
body approaches the fluid velocity.  A body with no relative velocity gets
no drag.

The projected area :math:`A_\perp` normal to the relative velocity is
computed every timestep.  The surface atom centers of the body (excluding
its center-of-mass atom), in body coordinates, are projected onto the plane
normal to the relative velocity, and the Qhull library computes the area *A*
and perimeter *P* of their two-dimensional convex hull.  The atom radius
*r* is then included exactly for the hull grown outwards by *r*:
:math:`A_\perp = A + P r + \pi r^2`, with *r* the largest radius of the
body's surface atoms.  This correction is always applied, independent of
the *offset* keyword of fix rigid/abrade.  For a non-convex body the convex
hull overestimates the silhouette.

With *buoyancy* = *yes*, the force :math:`-\rho_f V \mathbf{g}` is also
added to each body, using the gravity vector that fix rigid/abrade obtains
through its *gravity* keyword, which is then required.  The body is assumed
to be fully submerged and of uniform density, so no torque results.
Buoyancy does not depend on the flow and is added even without relative
velocity.

Output info
"""""""""""

This fix computes a global vector of length 6, which can be accessed by
various :doc:`output commands <Howto_output>`: the total drag force on all
bodies (1-3) and the total buoyancy force on all bodies (4-6), in force
units.  The vector values are "extensive".

This fix also computes a local array with one row per body and 10 columns,
which can be output with :doc:`dump local <dump>`: 1 body ID (atom ID of its
center-of-mass atom), 2 :math:`|u|`, 3 *Re*, 4 :math:`C_D`, 5 :math:`\Phi`,
6 :math:`\Phi_\perp`, 7 :math:`A_\perp`, 8-10 drag force.  Without relative
velocity, columns 3, 4 and 6-10 are zero.  *Re* and :math:`C_D` are
diagnostic only.  The values are those of the last timestep.  Because fix
rigid/abrade first calls this fix in the first timestep, the local array
has no rows at step 0 and the global vector is zero there.  Use a separate
dump local for this fix: it cannot share one with computes such as
:doc:`compute rigid/local_abrade <compute_rigid_local_abrade>` that have
rows at step 0.

No information about this fix is written to binary restart files.  None of
the :doc:`fix_modify <fix_modify>` options are relevant to this fix.

Restrictions
""""""""""""

This fix is part of the RIGID package and requires the Qhull library.  It
is only enabled if LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

This fix only works through the *drag* keyword of
:doc:`fix rigid/abrade <fix_rigid_abrade>`.  It currently supports only a
single MPI process.  It requires atoms with a radius.  The fluid velocity is
uniform and steady; lift forces, torques and velocity profiles are not
included.

Related commands
""""""""""""""""

:doc:`fix rigid/abrade <fix_rigid_abrade>`, :doc:`fix gravity <fix_gravity>`,
:doc:`compute projected/area <compute_projected_area>`

Default
"""""""

The default is *buoyancy* = *no*.

----------

.. _Holzer:

**(Holzer)** A. Holzer, M. Sommerfeld, New simple correlation formula for the drag coefficient of non-spherical particles, Powder Technol, 184, 361-365 (2008).
