.. _Equation_of_State:

Equation of State
=================

Within MOM6, there is a wrapper for the equation of state, so that all calls look the same from the rest of the model. The equation of state code has to calculate not just in situ or potential density, but also the compressibility and various derivatives of the density. There is also code for computing specific volume and the freezing temperature, and for converting between potential and conservative temperatures and between practical and reference (or absolute) salinity.

.. _Linear_EOS:

Linear Equation of State
------------------------

Compute the required quantities with uniform values for :math:`\alpha = \frac{\partial \rho}{\partial T}` and :math:`\beta = \frac{\partial \rho}{\partial S}`, (DRHO_DT, DRHO_DS in MOM_input, also uses RHO_T0_S0).

.. _Wright_EOS:

Wright reduced range Equation of State
--------------------------------------

Compute the required quantities using the equation of state from :cite:`wright1997` as a function of potential temperature and practical salinity, with coefficients based on the reduced-range (salinity from 28 to 38 PSU, temperature from -2 to 30 degC and pressure up to 5000 dbar) fit to the UNESCO 1981 data. This equation of state is in the form:

.. math::

   \alpha(s, \theta, p) = A(s, \theta) + \frac{\lambda(s, \theta)}{P(s, \theta) + p}

where :math:`A, \lambda` and :math:`P` are functions only of :math:`s` and :math:`\theta` and :math:`\alpha = 1/ \rho` is the specific volume. This form is useful for the pressure gradient computation as discussed in `Pressure Gradient Term <Discrete_PG.html#section-pg>`__. This EoS is selected by setting EQN_OF_STATE = WRIGHT or WRIGHT_RED, which are mathematically equivalent, but the latter is refactored for consistent answers between compiler settings.

.. _Wright_full_EOS:

Wright full range Equation of State
-----------------------------------

Compute the required quantities using the equation of state from :cite:`wright1997` as a function of potential temperature and practical salinity, with coefficients based on a fit to the UNESCO 1981 data over the full range of validity of that data (salinity from 0 to 40 PSU, temperatures from -2 to 40 degC, and pressures up to 10000 dbar). The functional form of the WRIGHT_FULL equation of state is the same as for WRIGHT or WRIGHT_RED, but with different coefficients.

.. _Jackett06_EOS:

Jackett et al. (2006) Equation of State
---------------------------------------

Compute the required quantities using the equation of state from Jackett et al. (2006) as a function of potential temperature and practical salinity, with coefficients based on a fit to the updated data that were later used to define the TEOS-10 equation of state over the full range of validity of that data (salinity from 0 to 42 PSU, temperatures from the freezing point to 40 degC, and pressures up to 8500 dbar), but focused on the “oceanographic funnel” of thermodynamic properties observed in the ocean. This equation of state is commonly used in realistic Hycom simulations.

.. _UNESCO_EOS:

UNESCO Equation of State
------------------------

Compute the required quantities using the equation of state from :cite:`jackett1995`, which uses potential temperature and practical salinity as state variables and is a fit to the 1981 UNESCO equation of state with the same functional form but a replacement of the temperature variable (the original uses *in situ* temperature).

.. _ROQUET_RHO_EOS:

ROQUET_RHO Equation of State
----------------------------

Compute the required quantities using the equation of state from :cite:`roquet2015`, which uses a 75-member polynomial for density as a function of conservative temperature and absolute salinity, in a fit to the output from the full TEOS-10 equation of state.

.. _ROQUET_SPV_EOS:

ROQUET_SPV Equation of State
----------------------------

Compute the required quantities using the specific volume oriented equation of state from :cite:`roquet2015`, which uses a 75-member polynomial for specific volume as a function of conservative temperature and absolute salinity, in a fit to the output from the full TEOS-10 equation of state.

.. _TEOS-10_EOS:

TEOS-10 Equation of State
-------------------------

Compute the required quantities using the equation of state from `TEOS-10 <http://www.teos-10.org/>`__, with calls directly to the subroutines in that code package.

.. _section_TFREEZE:

Freezing Temperature of Sea Water
---------------------------------

There are four choices for computing the freezing point of sea water:

- Linear The freezing temperature is a linear function of the salinity and pressure:

.. math::

   T_{Fr} = (T_{Fr0} + a\,S) + b\,P

where :math:`T_{Fr0},a,b` are constants which can be set in MOM_input (TFREEZE_S0_P0, DTFREEZE_DS, DTFREEZE_DP).

- Millero The :cite:`millero1978` equation is used to calculate the freezing point from practical salinity and pressure, but modified so that returns a potential temperature rather than an *in situ* temperature:

.. math::

   T_{Fr} = S(a + (b \sqrt{\max(S,0.0)} + c\, S)) + d\,P

where :math:`a,b, c, d` are fixed constants.

- TEOS-10 The TEOS-10 package is used to compute the freezing conservative temperature [degC] from absolute salinity [g/kg], and pressure [Pa]. This one or TEOS_poly must be used if you are using the ROQUET_RHO, ROQUET_SPV or TEOS-10 equation of state.

- TEOS_poly A 23-term polynomial fit refactored from the TEOS-10 package is used to compute the freezing conservative temperature [degC] from absolute salinity [g/kg], and pressure [Pa].
