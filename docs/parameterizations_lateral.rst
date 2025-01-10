Lateral Parameterizations
=========================

The following sub-grid scale parameterizations generally yield fluxes that act in the lateral direction.

Lateral viscosity
-----------------

Laplacian and bi-harmonic viscosities with linear and Smagorinsky options are implemented in MOM_hor_visc.

    :ref:`namespacemom__hor__visc_1section_horizontal_viscosity`

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
:cite:`fox-kemper2008-2` is implemented in MOM_mixed_layer_restrat,
which now also contains the mixed layer restratication comes from :cite: Bodner2023.

    :ref:`namespacemom__mixed__layer__restrat_1section_mle`

Interface filtering
-------------------

For layer mode, one can filter the interface thicknesses:

    :ref:`namespacemom__interface__filter_1section_interface_filter`

Lateral diffusion
-----------------

See :ref:`Horizontal_Diffusion`.

See also :ref:`namespacemom__lateral__mixing__coeffs_1section_Resolution_Function`

Tidal forcing
-------------

Astronomical tidal forcings and self-attraction and loading are implement in

    :ref:`namespacetidal__forcing_1section_tides`

The Love numbers are stored internally in MOM_load_love_numbers:

    :ref:`namespacemom__load__love__numbers_1section_Love_numbers`

while the self attraction and loading is computed in MOM_self_attr_load:

    :ref:`namespaceself__attr__load_1section_SAL`

The self attraction and loading needs spherical harmonics, computed in MOM_spherical_harmonics:

    :ref:`namespacemom__spherical__harmonics_1section_spherical_harmonics`

Tides can also be added via an open boundary tidal specification,
see `OBC wiki page <https://github.com/NOAA-GFDL/MOM6-examples/wiki/Open-Boundary-Conditions>`_.
