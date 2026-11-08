.. _EPBL:

Energetically-constrained Planetary Boundary Layer
==================================================

We here describe a scheme for modeling the ocean surface boundary layer (OSBL) suitable for use in global climate models. It builds on the ideas in `Bulk Surface Mixed Layer <BML.html#bml>`__, bringing in some of the ideas from `Shear-driven mixing in Jackson <Internal_Vert_Mixing.html#subsection-kappa-shear>`__, to make an energetically consistent boundary layer suitable for use with a generalized vertical coordinate. Unlike in `Bulk Surface Mixed Layer <BML.html#bml>`__, variables are allowed to have vertical structure within the boundary layer. The downward turbulent flux of buoyant water by OSBL turbulence converts mechanical energy into potential energy as it mixes with less buoyant water at the base of the OSBL. As described in :cite:`reichl2018`, we focus on OSBL parameterizations that constrain this integrated potential energy conversion due to turbulent mixing.

The leading-order mean OSBL equation for arbitrary scalar :math:`\phi` is:

.. math::

   \frac{\partial \overline{\phi}}{\partial t} = - \frac{\partial}{\partial z}
      \overline{w^\prime \phi^\prime} + \nu_\phi \frac{\partial^2 \overline{\phi}}{\partial z^2}

where the symbols are as follows:

.. list-table:: Symbols used in TKE equation
   :header-rows: 1

   * - Symbol
     - Meaning
   * - :math:`u_i`
     - horizontal components of the velocity
   * - :math:`\phi`
     - arbitrary scalar (tracer) quantity
   * - :math:`w`
     - vertical component of the velocity
   * - :math:`\overline{w}`
     - ensemble average :math:`w`
   * - :math:`w^\prime`
     - fluctuations from :math:`\overline{w}`
   * - :math:`k`
     - turbulent kinetic energy (TKE)
   * - :math:`K_M`
     - turbulent mixing coefficient for momentum
   * - :math:`K_\phi`
     - turbulent mixing coefficient for :math:`\phi`
   * - :math:`\sigma_k`
     - turbulent Schmidt number
   * - :math:`b`
     - buoyancy
   * - :math:`\epsilon`
     - buoyancy turbulent dissipation rate


This equation describes the evolution of mean quantity :math:`\overline{\phi}` due to vertical processes, including the often negligible molecular mixing. We would like to parameterize the vertical mixing since we won’t be resolving all the relevant time and space scales.

We use the Boussinesq hypothesis for turbulence closure. This approximates the Reynolds stress terms using an eddy viscosity (eddy diffusivity for turbulent scalar fluxes):

.. math::

   \overline{u_i^\prime w^\prime} = - K_M \frac{\partial \overline{u_i}}{\partial z} ,

Similarly, the eddy diffusivity is used to parameterize turbulent scalar fluxes as:

.. math::

   \overline{\phi^\prime w^\prime} = - K_\phi \frac{\partial \overline{\phi}}{\partial z} ,

The parameters needed to close the system of equations are then reduced to the turbulent mixing coefficients, :math:`K_\phi` and :math:`K_M`.

We start with an equation for the turbulent kinetic energy (TKE):

.. math::

   \frac{\partial k}{\partial t} = \frac{\partial}{\partial z} \left( \frac{K_M}{\sigma_k}
     \frac{\partial k}{\partial z} \right) - \overline{u_i^\prime w^\prime} \frac{\partial \overline{u_i}}
     {\partial z} + \overline{w^\prime b^\prime} - \epsilon

Terms in this equation represent TKE storage (LHS), TKE flux convergence, shear production, buoyancy production, and dissipation.

.. _section_WMBL:

Well-mixed Boundary Layers (WMBL)
---------------------------------

Assuming steady state and other parameterizations, integrating vertically over the surface boundary layer, :cite:`reichl2018` obtains the form:

.. math::

   \frac{1}{2} H_{bl} w_e \Delta b = m_\ast u_\ast^3 - n_\ast \frac{H_{bl}}{2}
      B(H_{bl}) ,

with the following variables:

.. list-table:: Symbols used in integrated TKE equation
   :header-rows: 1

   * - Symbol
     - Meaning
   * - :math:`H_{bl}`
     - boundary layer thickness
   * - :math:`w_e`
     - entrainment velocity
   * - :math:`\Delta b`
     - change in buoyancy at base of mixed layer
   * - :math:`m_\ast`
     - sum of mechanical coefficients
   * - :math:`u_\ast`
     - friction velocity (:math:`u_\ast = (|\tau| / \rho_0)^{1/2}`)
   * - :math:`\tau`
     - wind stress
   * - :math:`n_\ast`
     - convective proportionality coefficient
   * -
     - 1 for stabilizing surface buoyancy flux, less otherwise
   * - :math:`B(H_{bl})`
     - surface buoyancy flux


.. _section_ePBL:

Energetics-based Planetary Boundary Layer
-----------------------------------------

Once again, the goal is to formulate a surface mixing scheme to find the turbulent eddy diffusivity (and viscosity) in a way that is suitable for use in a global climate model, using long timesteps and large grid spacing. After evaluating a well-mixed boundary layer (WMBL), the shear mixing of :cite:`jackson2008` (JHL, `Shear-driven mixing in Jackson <Internal_Vert_Mixing.html#subsection-kappa-shear>`__), as well as a more complete boundary layer scheme, it was decided to combine a number of these ideas into a new scheme:

.. math::

   K(z) = F_x(K_{ePBL}(z), K_{JHL}(z), K_n(z))

where :math:`F_x` is some unknown function of a new :math:`K_{ePBL}`, :math:`K_{JHL}`, the diffusivity due to shear as determined by :cite:`jackson2008`, and :math:`K_n`, the diffusivity from other ideas. We start by specifying the form of :math:`K_{ePBL}` as being:

.. math::

   K_{ePBL}(z) = C_K w_t l ,

where :math:`w_t` is a turbulent velocity scale, :math:`C_K` is a coefficient, and :math:`l` is a length scale.

.. _subsection_lengthscale:

Turbulent length scale
~~~~~~~~~~~~~~~~~~~~~~

We propose a form for the length scale as follows:

.. math::

   l = (z_0 + |z|) \times \max \left[ \frac{l_b}{H_{bl}} , \left(
      \frac{H_{bl} - |z|}{H_{bl}} \right)^\gamma \, \right] ,

where we have the following variables:

.. list-table:: Symbols used in ePBL length scale
   :header-rows: 1

   * - Symbol
     - Meaning
   * - :math:`H_{bl}`
     - boundary layer thickness
   * - :math:`z_0`
     - roughness length
   * - :math:`\gamma`
     - coefficient, 2 is as in KPP, :cite:`large1994`
   * - :math:`l_b`
     - bottom length scale


.. _subsection_velocityscale:

Turbulent velocity scale
~~~~~~~~~~~~~~~~~~~~~~~~

We do not predict the TKE prognostically and therefore approximate the vertical TKE profile to estimate :math:`w_t`. An estimate for the mechanical contribution to the velocity scale follows the standard two-equation approach. In one and two-equation second-order :math:`K` parameterizations the boundary condition for the TKE is typically employed as a flux boundary condition.

.. math::

   K \left. \frac{\partial k}{\partial z} \right|_{z=0} = c_\mu^0 u_\ast^3 .

The profile of :math:`k` decays in the vertical from :math:`k \propto (c_\mu^0)^{2/3} u_\ast^2` toward the base of the OSBL. Here we assume a similar relationship to estimate the mechanical contribution to the TKE profile. The value of :math:`w_t` due to mechanical sources, :math:`v_\ast`, is estimate as :math:`v_\ast (z=0) \propto (c_\mu^0)^{1/3} u_\ast` at the surface. Since we only parameterize OSBL turbulent mixing due to surface forcing, the value of the velocity scale is assumed to decay moving away from the surface. For simplicity we employ a linear decay in depth:

.. math::

   v_\ast (z) = (c_\mu^0)^{1/3} u_\ast \left( 1 - a \cdot \min \left[ 1,
      \frac{|z|}{H_{bl}} \right] \right) ,

where :math:`1 > a > 0` has the effect of making :math:`v_\ast(z=H_{bl}) > 0`. Making the constant coefficient :math:`a` close to one has the effect of reducing the mixing rate near the base of the boundary layer, thus producing a more diffuse entrainment region. Making :math:`a` close to zero has the effect of increasing the mixing at the base of the boundary layer, producing a more ‘step-like’ entrainment region.

An estimate for the buoyancy contribution is found utilizing the convective velocity scale:

.. math::

   w_\ast (z) = C_{w_\ast} \left( \int_z^0 \overline{w^\prime b^\prime} dz \right)^{1/3} ,

where :math:`C_{w_\ast}` is a non-dimensional empirical coefficient. Convection in one and two-equation closure causes a TKE profile that peaks below the surface. The quantity :math:`\overline{w^\prime b^\prime}` is solved for in ePBL as :math:`KN^2`.

These choices for the convective and mechanical components of the velocity scale in the OSBL are then added together to get an estimate for the total turbulent velocity scale:

.. math::

   w_t (z) = w_\ast (z) + v_\ast (z) .

The value of :math:`a` is arbitrarily chosen to be 0.95 here.

.. _subsection_ePBL_summary:

Summarizing the ePBL implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ePBL mixing coefficient is found by multiplying a velocity scale (`Turbulent velocity scale <EPBL.html#subsection-velocityscale>`__) by a length scale (`Turbulent length scale <EPBL.html#subsection-lengthscale>`__). The precise value of the coefficient :math:`C_K` used does not significantly alter the prescribed potential energy change constraint. A reasonable value is :math:`C_K \approx 0.55` to be consistent with other approaches (e.g. :cite:`umlauf2005`).

The boundary layer thickness (:math:`H_{bl}`) within ePBL is based on the depth where the energy requirement for turbulent mixing of density exceeds the available energy (`Well-mixed Boundary Layers (WMBL) <EPBL.html#section-wmbl>`__). :math:`H_{bl}` is determined by the energetic constraint imposed using the value of :math:`m_\ast` and :math:`n_\ast`. An iterative solver is required because :math:`m_\ast` and the mixing length are dependent on :math:`H_{bl}`.

We use a constant value for convectively driven TKE of :math:`n_\ast = 0.066`. The parameterizations for :math:`m_\ast` are formulated specifically for the regimes where :math:`K_{JHL}` is sensitive to model numerics :math:`(|f| \Delta t \approx 1)` (:cite:`reichl2018`).

.. _subsection_ePBL_JHL:

Combining ePBL and JHL mixing coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We now address the combination of the ePBL mixing coefficient and the JHL mixing coefficient. The function :math:`F_x` above cannot be the linear sum of :math:`K_{ePBL}` and :math:`K_{JHL}`. One reason this sum is not valid is because the JHL mixing coefficient is determined by resolved current shear, including that driven by the surface wind. The wind-driven current is also included in the ePBL mixing coefficient formulation. An alternative approach is therefore needed to avoid double counting.

:math:`K_{ePBL}` is not used at the equator as scalings are only investigated when :math:`|f| > 0`. The solution we employ is to use the maximum mixing coefficient of the two contributions,

.. math::

   K (z) = \max (K_{ePBL} (z), K_{JHL} (z)),

where :math:`m_\ast` (and hence :math:`K_{ePBL}`) is constrained to be small as :math:`|f| \rightarrow 0`. This form uses the JHL mixing coefficient when the ePBL coefficient is small.

This approach is reasonable when the wind-driven mixing dominates, since both JHL and ePBL give a similar solution when deployed optimally. One weakness of this approach is the tropical region, where the shear-driven ePBL :math:`m_\ast` coefficient is not formulated. The JHL parameterization is skillful to simulate this mixing, but does not include the contribution of convection. The convective portion of :math:`K_{ePBL}` should be combined with :math:`K_{JHL}` in the equatorial region when shear and convection occur together. Future research is warranted.

Finally, one should note that the mixing coefficient here (:math:`K`) is used for both diffusivity and viscosity, implying a turbulent Prandtl number of 1.0.

.. _subsection_Langmuir:

Langmuir circulation
~~~~~~~~~~~~~~~~~~~~

While only briefly alluded to in :cite:`reichl2018`, the MOM6 code implementing ePBL does support the option to add a Langmuir parameterization. There are in fact two options, both adjusting :math:`m_\ast`.
