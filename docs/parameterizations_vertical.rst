Vertical Parameterizations
==========================

The following sub-grid scale parameterizations generally yield fluxes that act in the vertical direction, with no lateral components resolved by the model.

Upper boundary
--------------

K-profile parameterization (KPP)
  Provided by module MOM_KPP, uses the CVmix implementation of KPP.

   :ref:`CVMix_KPP`

Energetic Planetary Boundary Layer (ePBL)
  A energetically constrained boundary layer scheme following :cite:`reichl2018`. Implemented in MOM_energetic_PBL.

Bulk mixed layer (BML)
  A 2-layer bulk mixed layer used in pure-isopycnal model. Implemented in MOM_bulk_mixed_layer.

Interior and bottom-driven mixing
---------------------------------

Kappa-shear
  MOM_kappa_shear implement the shear-driven mixing of :cite:`jackson2008`.

Internal-tide driven mixing
  The schemes of :cite:`st_laurent2002`, :cite:`polzin2009`, and :cite:`melet2012`, are all implemented through MOM_set_diffusivity and MOM_diabatic_driver.

Vertical friction
-----------------

Vertical viscosity is implemented in MOM_vert_frict and coefficient computed in MOM_set_viscosity, although contributions to viscosity from other parameterizations are calculated in those respective modules (e.g. MOM_kappa_shear, MOM_KPP, MOM_energetic_PBL).

Vertical diffusion
------------------

Vertical diffusion of scalars is implemented in MOM_diabatic_driver although contributions to diffusion from other parameterizations are calculated in those respective modules (e.g. MOM_kappa_shear, MOM_KPP, MOM_energetic_PBL).

Radiation
---------

Opacity
  Ocean color is prescribed or dynamically calculated in converted into optical properties in MOM_opacity.

Short-wave absorption
  Optical properties from MOM_opacity are used to calculate the convergence of shortwave radiation penetrating from the upper surface in MOM_shortwave_abs.

Geothermal heating
------------------

Geothermal heat fluxes are implemented in MOM_geothermal.

Isopycnal-mode entrainment and diapycnal diffusion
--------------------------------------------------

Diapycnal diffusion in a layered isopycnal mode following :cite:`hallberg2000`, is implemented in MOM_entrain_diffuse.
