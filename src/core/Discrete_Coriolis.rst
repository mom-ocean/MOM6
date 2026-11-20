.. _Discrete_Coriolis:

Discrete Coriolis Term
======================

.. _Coriolis:

Coriolis Term
-------------

In general, the discrete equations are written as simple difference equations based on the Arakawa C-grid as described in section `Horizontal grids <Discrete_Grids.html#horizontal-grids>`__. One of the more interesting exceptions is the Coriolis term. It is computed in the form shown in `momentum equation <General_Coordinate.html#eq-h-horz-momentum>`__, or:

.. math::

   \frac{( f + \zeta )}{h} \, \hat{\mathbf{z}} \times h \, \mathbf{u}

This term needs to be evaluated at :math:`u` points for the :math:`v` equation and vice versa, plus we need to keep the thickness, :math:`h`, positive definite. MOM6 contains a number of options for how to compute this term.

- SADOURNY75_ENERGY Sadourny :cite:`sadourny1975` figured out how to conserve energy or enstrophy but not both. This option is energy conserving. The term in the :math:`u` equation becomes:

.. math::

   \frac{1}{4 dx} \left( q_{i,j} (vh_{i+1,j} + vh_{i,j}) +
       q_{i,j-1} (vh_{i+1,j-1} + vh_{i,j-1}) \right)

where :math:`q = \frac{f + \zeta}{h}` and :math:`h` is an area-weighted average of the four thicknesses surrounding the :math:`q` point, such that it is guaranteed to be positive definite.

There is a variant on this scheme with the CORIOLIS_EN_DIS option. If true, two estimates of the thickness fluxes :math:`vh` are used to estimate the Coriolis term, and the one that dissipates energy relative to the other one is used.

- SADOURNY75_ENSTRO Also from :cite:`sadourny1975`, this option is enstrophy conserving.

.. math::

   \frac{1}{8 dx} ( q_{i,j} + q_{i,j-1} ) ((vh_{i+1,j} + vh_{i,j}) +
       (vh_{i+1,j-1} + vh_{i,j-1}) )

- ARAKAWA_LAMB81 From :cite:`arakawa1981` is a scheme which is both energy and enstrophy conserving. Its weaknesses are a large stencil and differing thickness stencils in the numerator and denominator. This scheme and several others (with differing values of :math:`a, b, c, d` and :math:`ep`) are implemented as:

.. _eq:Coriolis_abcd:

.. math::

   \begin{aligned}
       \frac{1}{dx} (a_{i,j} vh_{i+1,j} &+ b_{i,j} vh_{i,j} +
       d_{i,j} vh_{i+1,j-1} + c_{i,j} vh_{i,j-1} \\
       &+ ep_{i,j}*uh_{i-1,j} -
       ep_{i+1,j}*uh_{i+1,j})
   \end{aligned}

with

.. math::

   \begin{aligned}
      a_{i,j} &= \frac{1}{24} (2.0*(q_{i+1,j} + q_{i,j-1}) + (q_{i,j} + q_{i+1,j-1})) \\
      b_{i,j} &= \frac{1}{24} ((q_{i,j} + q_{i-1,j-1}) + 2.0*(q_{i-1,j} + q_{i,j-1})) \\
      c_{i,j} &= \frac{1}{24} (2.0*(q_{i,j} + q_{i-1,j-1}) + (q_{i-1,j} + q_{i,j-1})) \\
      d_{i,j} &= \frac{1}{24} ((q_{i+1,j} + q_{i,j-1}) + 2.0*(q_{i,j} + q_{i+1,j-1})) \\
      ep_{i,j} &= \frac{1}{24}((q_{i,j} - q_{i-1,j-1}) + (q_{i-1,j} - q_{i,j-1}))
   \end{aligned}

- ARAKAWA_HSU90 From :cite:`arakawa1990` is a scheme which always conserves energy and conserves enstrophy in the limit of non-divergent flow. This one has a larger stencil than Sadourny’s energy scheme, but it’s much better behaved in terms of handling vanishing layers than Arakawa and Lamb. This scheme is implemented with:

.. math::

   \frac{1}{dx} (a_{i,j} vh_{i+1,j} + b_{i,j} vh_{i,j} +
       d_{i,j} vh_{i+1,j-1} + c_{i,j} vh_{i,j-1})

and

.. math::

   \begin{aligned}
      a_{i,j} &= \frac{1}{12} (q_{i,j} + (q_{i+1,j} + q_{i,j-1})) \\
      b_{i,j} &= \frac{1}{12} (q_{i,j} + (q_{i-1,j} + q_{i,j-1})) \\
      c_{i,j} &= \frac{1}{12} (q_{i,j} + (q_{i-1,j-1} + q_{i,j-1})) \\
      d_{i,j} &= \frac{1}{12} (q_{i,j} + (q_{i+1,j-1} + q_{i,j-1}))
   \end{aligned}

- ARAKAWA_LAMB_BLEND This is a blending of Arakawa and Lamb, Arakawa and Hsu, and the Sadourny Energy scheme. There are weights CORIOLIS_BLEND_WT_LIN and CORIOLIS_BLEND_F_EFF_MAX to control this scheme. The equation is the same as for Arakawa and Lamb `Coriolis discretization <Discrete_Coriolis.html#eq-coriolis-abcd>`__, but the values of :math:`a, b, c, d` and :math:`ep` differ when the pure Arakawa and Lamb scheme breaks down due to thickness variations.

- ROBUST_ENSTRO An enstrophy-conserving scheme which is robust to vanishing layers.

Some of these options also support the BOUND_CORIOLIS flag. If true, the Coriolis terms in the :math:`u` equation are bounded by the four estimates of :math:`\frac{(f+\zeta)}{h}vh` from the four neighboring :math:`v` points, and similarly in the :math:`v` equation. This option would have no effect on the SADOURNY75_ENERGY scheme if it were possible to use centered difference thickness fluxes.

Note, if BOUND_CORIOLIS is on, it will also turn on the BOUND_CORIOLIS_BIHARM option by default. This option uses a viscosity that increases with the square of the velocity shears, so that the resulting viscous drag is of comparable magnitude to the Coriolis terms when the velocity differences between adjacent grid points is 0.5*BOUND_CORIOLIS_VEL.

.. _Coriolis_BC:

Wall boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~

Two sets of boundary conditions have been coded in the definition of relative vorticity. These are written as:

NOSLIP defined (in spherical coordinates):

.. math::

   \begin{aligned}
     \mbox{relvort} &= dv/dx \mbox{ (east $\&$ west)}, \mbox{ with } v = 0. \\
     \mbox{relvort} &= -\sec(\phi) * d(u \cos(\phi))/dy \mbox{ (north $\&$
     south)}, \mbox{ with } u = 0.
   \end{aligned}

Free slip (NOSLIP not defined):

.. math::

   \mbox{relvort} = 0 \mbox{ (all boundaries)}

with :math:`\phi` defined as latitude. The free slip boundary condition is much more natural on a C-grid.
