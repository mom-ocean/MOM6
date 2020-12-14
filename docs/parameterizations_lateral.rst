Lateral Parameterizations
=========================

The following sub-grid scale parameterizations generally yield fluxes that act in the lateral direction.

Lateral viscosity
-----------------

Laplacian and bi-harmonic viscosities with linear and Smagorinsky options are implemented in MOM_hor_visc.

Gent-McWilliams/TEM/isopycnal height diffusion
----------------------------------------------

Lagrangian mean eddy mass transport is parameterized following :cite:`gent1990`, in MOM_thickness_diffuse.

The diffusivity coefficients are calculated in MOM_lateral_mixing_coeffs
and MOM_thickness_diffuse and includes constants and the :cite:`visbeck1997`
scaling.

A model of sub-grid scale Mesoscale Eddy Kinetic Energy (MEKE) is implement in MOM_MEKE and the associated diffusivity added in MOM_thickness_diffuse.
See :cite:`jansen2015` and :cite:`marshall2010`.

   :ref:`namespacemom__meke_1section_MEKE`

Backscatter
-----------

A parameterization of the upscale unresolved cascade utilizes MOM_MEKE
and negative Laplacian viscosity in MOM_hor_visc.

Mixed layer restratification by sub-mesoscale eddies
----------------------------------------------------

Mixed layer restratification from :cite:`fox-kemper2008` and
:cite:`fox-kemper2008-2` is implemented in MOM_mixed_layer_restrat.

Lateral diffusion
-----------------

See :ref:`Horizontal_Diffusion`.

Tidal forcing
-------------

Astronomical tidal forcings and self-attraction and loading are implement in MOM_tidal_forcing.

