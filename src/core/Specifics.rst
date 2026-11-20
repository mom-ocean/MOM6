.. _Specifics:

Specifics
=========

.. _section_Specifics:

Specifics of the Ocean Model Equations
--------------------------------------

We here provide more details of the terms appearing in the ocean model equations described in `Generalized vertical coordinate equations <General_Coordinate.html#general-coordinate>`__.

.. _Horiz_mom_eq:

Horizontal Momentum Equation
----------------------------

Equation `momentum equation <General_Coordinate.html#eq-h-horz-momentum>`__ is the horizontal momentum equation written in its vector-invariant advective form [#note1]_ with :math:`\mathbf{u} = \hat{\mathbf{x}} \, u + \hat{\mathbf{y}} \, v` the horizontal velocity, :math:`p` the hydrostatic pressure, :math:`f` the Coriolis parameter, and

.. math::

   \zeta^{(r)} = \hat{\bf z} \cdot (\nabla_{r}\times \mathbf{u})

the vertical component of the vorticity using :math:`\nabla_{r}` for the curl operator. The discretization of the Coriolis term is the enstrophy conserving scheme of :cite:`sadourny1975`. The geopotential coordinate, :math:`z`, has a value :math:`z=0` at a resting ocean surface, :math:`z=\eta(x,y,t)` at the ocean free surface, and :math:`z=-H(x,y)` at the ocean bottom. We use the Boussinesq approximation (volume conserving kinematics) with :math:`\rho_{0} = 1035~\mbox{kg}~\mbox{m}^{-3}` the reference density. [#note2]_ Time and horizontal derivatives are computed holding the generalized vertical coordinate fixed rather than the geopotential

.. math::

   \frac{\partial}{\partial t} = \left[ \frac{\partial}{\partial t} \right]_{r} \qquad \nabla_{r} = \hat{\mathbf{x}} \left[ \frac{\partial}{\partial x} \right]_{r} + \hat{\mathbf{y}} \left[ \frac{\partial}{\partial y} \right]_{r}.

The transport of seawater crossing surfaces of constant :math:`r` is measured by the dia-surface velocity component (see section 6.7 of :cite:`SMGbook`)

.. math::

   \frac{\partial z}{\partial r} \, \frac{\mathrm{D}r}{\mathrm{D}t} = z_{r} \, \dot{r},

with :math:`z_{r}` the specific thickness that is assumed one-signed throughout the ocean, and :math:`\mathrm{D}/\mathrm{D}t` the material time derivative operator. In the ocean interior where :math:`r` is aligned with isopycnals, the dia-surface velocity becomes the diapycnal velocity whose value is directly related to irreversible processes such as mixing that act on potential temperature and salinity. In the unstratified mixed layers, :math:`r =z^{*}` so that :math:`z_{r} \dot{r} = (\partial z/\partial z^{*}) \, \mathrm{D}z^{*}/\mathrm{D}t`, which is close to the familiar vertical velocity component :math:`\mathrm{D}z/\mathrm{D}t`.

Viscous dissipation (Laplacian and biharmonic friction following :cite:`griffies2000`) and mechanical boundary forces (winds, bottom stress) contribute to the divergence of the deviatoric (symmetric and trace-free) stress tensor, :math:`\boldsymbol{\mathcal{F}} = \nabla \cdot \mathbf{\tau}`. MOM6 and the real ocean have no vertical sidewalls, and MOM6 treats all solid-earth boundaries with bottom stress parameterized as a quadratic drag.

.. _hydrostatic_balance:

Hydrostatic balance
-------------------

Equation `hydrostatic equation <General_Coordinate.html#eq-h-hydrostatic-equation>`__ is the discrete version of the hydrostatic balance. The horizontal pressure gradient force is implemented as a contact force following the method of :cite:`adcroft2008`. These equations differ from :cite:`bleck2002` who uses the Montgomery potential to calculate pressure gradient accelerations.

.. _Thickness_and_tracer:

Thickness and tracer equations
------------------------------

Volume conservation appears in the form of a prognostic flux-form layer thickness equation `thickness equation <General_Coordinate.html#eq-h-thickness-equation>`__, with the non-negative layer thickness given by

.. math::

   h = \frac{\partial z}{\partial r} \, \mathrm{d}r,

where :math:`\mathrm{d}r` is the thickness of a layer in :math:`r`-space (e.g., the density difference between target density classes or the thickness between target depths). The layer thickness increases where horizontal thickness fluxes converge, :math:`\nabla_r \cdot \left( h_k \, \mathbf{u} \right) < 0`, and where dia-surface flow converges, :math:`\delta_r (z_{r} \, \dot{r} ) < 0`. The volume flux :math:`h_k \mathbf{u}` is computed using the quasi-third order PPM scheme (:cite:`colella1984`) using a positive-definite limiter rather than the monotonic limiter. This last choice avoids limiting of positive extrema and thus retains third-order accuracy everywhere except near vanishing layers.

Transport in the thickness equation is discretized compatibly with that in the flux-form potential temperature and salinity equations `potential temperature equation <General_Coordinate.html#eq-h-temperature-equation>`__ and `salinity equation <General_Coordinate.html#eq-h-salinity-equation>`__. Compatibility is required to maintain global and layer integrated conservation properties for volume, heat, and salt. Tracer reconstruction for transport uses PPM with monotonic limiters but using third order interpolation for edge values. This reduces the size of the stencil which helps the computational efficiency of the transport scheme. The flux convergences, :math:`\boldsymbol{\mathcal{N}}_\theta^\gamma` and :math:`\boldsymbol{\mathcal{N}}_S^\gamma`, provide subgrid scale neutral diffusion for the potential temperature and salinity, whereas :math:`\delta_{r}J_{\theta}^{(z)}` and :math:`\delta_{r}J_{S}^{(z)}` provide subgrid scale vertical diffusion as well as boundary fluxes. In the interior, both subgrid fluxes vanish when their respective tracers are spatially uniform, thus ensuring that the tracer equation reduces to the thickness equation when the tracer is uniform.

Parameterized subgrid scale advection from the submesoscale (:cite:`fox-kemper2011`) and mesoscale (:cite:`gent1995`) parameterizations are combined with the lateral advection of thickness and tracer, thus providing a residual mean advective transport for the scalar fields. Furthermore, we implement subgrid advective terms solely as lateral transports, thus interpreting them as layer bolus transport as appropriate for vertical Lagrangian models rather than a three-dimensional eddy-induced advection as appropriate for vertical Eulerian models (see :cite:`mcdougall2001` for details).

.. _EOS:

Equation of state
-----------------

The equation of state, `equation of state equation <General_Coordinate.html#eq-h-equation-of-state>`__, determines *in situ* density as a function of potential temperature, salinity, and pressure. We evaluate the pressure in the equation of state according to :math:`- g \, \rho_{0} \, z`. Doing so maintains energetic consistency for the Boussinesq fluid according to section 2.4.3 of :cite:`GVbook`. We make use of the :cite:`wright1997` equation of state so that :math:`\theta` is potential temperature and :math:`S` is the practical salinity. Although MOM6 has the more updated equation of state from :cite:`TEOS2010`, the required changes for thermodynamic variables were implemented only after the basic model configuration was developed. Time constraints on model development prompted us to retain usage of :cite:`wright1997` for OM4.

The freezing point of seawater is approximated as

.. math::

   T_f = -0.054 S - 7.75\times10^{-08} p,

where :math:`p` is in units of Pascals and :math:`S` is in units of :math:`1\times10^{-3}`. When the local temperature anywhere in the ocean column falls below the freezing point, the water-equivalent volume of ice is calculated and the fusion heat locally added back to the ocean to raise the liquid seawater temperature back to the freezing point. The frozen water and salt are sent to the sea-ice model.


.. [#note1] The vector-invariant advective form is commonly used in models such as MOM6 using an Arakawa C-grid (e.g., see section 10 of :cite:`griffies2000-2`).

.. [#note2] MOM6 has an option for compressible non-Boussinesq flow (mass conserving kinematics). We chose the Boussinesq option largely based on legacy.
