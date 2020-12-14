!> This module implements boundary forcing for MOM6.
module MOM_forcing_type

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_array_transform, only : rotate_array, rotate_vector, rotate_array_pair
use MOM_debugging,     only : hchksum, uvchksum
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, register_scalar_field
use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_alloc, query_averaging_enabled
use MOM_diag_mediator, only : enable_averages, enable_averaging, disable_averaging
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_EOS,           only : calculate_density_derivs, EOS_domain
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_opacity,       only : sumSWoverBands, optics_type, extract_optics_slice, optics_nbands
use MOM_spatial_means, only : global_area_integral, global_area_mean
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : surface, thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

use coupler_types_mod, only : coupler_2d_bc_type, coupler_type_spawn
use coupler_types_mod, only : coupler_type_increment_data, coupler_type_initialized
use coupler_types_mod, only : coupler_type_copy_data, coupler_type_destructor

implicit none ; private

#include <MOM_memory.h>

public extractFluxes1d, extractFluxes2d, optics_type
public MOM_forcing_chksum, MOM_mech_forcing_chksum
public calculateBuoyancyFlux1d, calculateBuoyancyFlux2d
public forcing_accumulate, fluxes_accumulate
public forcing_SinglePointPrint, mech_forcing_diags, forcing_diagnostics
public register_forcing_type_diags, allocate_forcing_type, deallocate_forcing_type
public copy_common_forcing_fields, allocate_mech_forcing, deallocate_mech_forcing
public set_derived_forcing_fields, copy_back_forcing_fields
public set_net_mass_forcing, get_net_mass_forcing
public rotate_forcing, rotate_mech_forcing

!> Allocate the fields of a (flux) forcing type, based on either a set of input
!! flags for each group of fields, or a pre-allocated reference forcing.
interface allocate_forcing_type
  module procedure allocate_forcing_by_group
  module procedure allocate_forcing_by_ref
end interface allocate_forcing_type

!> Allocate the fields of a mechanical forcing type, based on either a set of
!! input flags for each group of fields, or a pre-allocated reference forcing.
interface allocate_mech_forcing
  module procedure allocate_mech_forcing_by_group
  module procedure allocate_mech_forcing_from_ref
end interface allocate_mech_forcing

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Structure that contains pointers to the boundary forcing used to drive the
!! liquid ocean simulated by MOM.
!!
!! Data in this type is allocated in the module MOM_surface_forcing.F90, of which there
!! are three: solo, coupled, and ice-shelf. Alternatively, they are allocated in
!! MESO_surface_forcing.F90, which is a special case of solo_driver/MOM_surface_forcing.F90.
type, public :: forcing

  ! surface stress components and turbulent velocity scale
  real, pointer, dimension(:,:) :: &
    ustar         => NULL(), & !< surface friction velocity scale [Z T-1 ~> m s-1].
    ustar_gustless => NULL()   !< surface friction velocity scale without any
                               !! any augmentation for gustiness [Z T-1 ~> m s-1].

  ! surface buoyancy force, used when temperature is not a state variable
  real, pointer, dimension(:,:) :: &
    buoy          => NULL()  !< buoyancy flux [L2 T-3 ~> m2 s-3]

  ! radiative heat fluxes into the ocean [Q R Z T-1 ~> W m-2]
  real, pointer, dimension(:,:) :: &
    sw         => NULL(), & !< shortwave [Q R Z T-1 ~> W m-2]
    sw_vis_dir => NULL(), & !< visible, direct shortwave [Q R Z T-1 ~> W m-2]
    sw_vis_dif => NULL(), & !< visible, diffuse shortwave [Q R Z T-1 ~> W m-2]
    sw_nir_dir => NULL(), & !< near-IR, direct shortwave [Q R Z T-1 ~> W m-2]
    sw_nir_dif => NULL(), & !< near-IR, diffuse shortwave [Q R Z T-1 ~> W m-2]
    lw         => NULL()    !< longwave [Q R Z T-1 ~> W m-2] (typically negative)

  ! turbulent heat fluxes into the ocean [Q R Z T-1 ~> W m-2]
  real, pointer, dimension(:,:) :: &
    latent           => NULL(), & !< latent [Q R Z T-1 ~> W m-2] (typically < 0)
    sens             => NULL(), & !< sensible [Q R Z T-1 ~> W m-2] (typically negative)
    seaice_melt_heat => NULL(), & !< sea ice and snow melt or formation [Q R Z T-1 ~> W m-2] (typically negative)
    heat_added       => NULL()    !< additional heat flux from SST restoring or flux adjustments [Q R Z T-1 ~> W m-2]

  ! components of latent heat fluxes used for diagnostic purposes
  real, pointer, dimension(:,:) :: &
    latent_evap_diag    => NULL(), & !< latent [Q R Z T-1 ~> W m-2] from evaporating liquid water (typically < 0)
    latent_fprec_diag   => NULL(), & !< latent [Q R Z T-1 ~> W m-2] from melting fprec  (typically < 0)
    latent_frunoff_diag => NULL()    !< latent [Q R Z T-1 ~> W m-2] from melting frunoff (calving) (typically < 0)

  ! water mass fluxes into the ocean [R Z T-1 ~> kg m-2 s-1]; these fluxes impact the ocean mass
  real, pointer, dimension(:,:) :: &
    evap        => NULL(), & !< (-1)*fresh water flux evaporated out of the ocean [R Z T-1 ~> kg m-2 s-1]
    lprec       => NULL(), & !< precipitating liquid water into the ocean [R Z T-1 ~> kg m-2 s-1]
    fprec       => NULL(), & !< precipitating frozen water into the ocean [R Z T-1 ~> kg m-2 s-1]
    vprec       => NULL(), & !< virtual liquid precip associated w/ SSS restoring [R Z T-1 ~> kg m-2 s-1]
    lrunoff     => NULL(), & !< liquid river runoff entering ocean [R Z T-1 ~> kg m-2 s-1]
    frunoff     => NULL(), & !< frozen river runoff (calving) entering ocean [R Z T-1 ~> kg m-2 s-1]
    seaice_melt => NULL(), & !< snow/seaice melt (positive) or formation (negative) [R Z T-1 ~> kg m-2 s-1]
    netMassIn   => NULL(), & !< Sum of water mass flux out of the ocean [kg m-2 s-1]
    netMassOut  => NULL(), & !< Net water mass flux into of the ocean [kg m-2 s-1]
    netSalt     => NULL()    !< Net salt entering the ocean [kgSalt m-2 s-1]

  ! heat associated with water crossing ocean surface
  real, pointer, dimension(:,:) :: &
    heat_content_cond    => NULL(), & !< heat content associated with condensating water [Q R Z T-1 ~> W m-2]
    heat_content_lprec   => NULL(), & !< heat content associated with liquid >0 precip   [Q R Z T-1 ~> W m-2]
    heat_content_icemelt => NULL(), & !< heat content associated with snow and seaice
                                      !! melt and formation [Q R Z T-1 ~> W m-2]
    heat_content_fprec   => NULL(), & !< heat content associated with frozen precip      [Q R Z T-1 ~> W m-2]
    heat_content_vprec   => NULL(), & !< heat content associated with virtual >0 precip  [Q R Z T-1 ~> W m-2]
    heat_content_lrunoff => NULL(), & !< heat content associated with liquid runoff      [Q R Z T-1 ~> W m-2]
    heat_content_frunoff => NULL(), & !< heat content associated with frozen runoff      [Q R Z T-1 ~> W m-2]
    heat_content_massout => NULL(), & !< heat content associated with mass leaving ocean [Q R Z T-1 ~> W m-2]
    heat_content_massin  => NULL()    !< heat content associated with mass entering ocean [Q R Z T-1 ~> W m-2]

  ! salt mass flux (contributes to ocean mass only if non-Bouss )
  real, pointer, dimension(:,:) :: &
    salt_flux       => NULL(), & !< net salt flux into the ocean [R Z T-1 ~> kgSalt m-2 s-1]
    salt_flux_in    => NULL(), & !< salt flux provided to the ocean from coupler [R Z T-1 ~> kgSalt m-2 s-1]
    salt_flux_added => NULL()    !< additional salt flux from restoring or flux adjustment before adjustment
                                 !! to net zero [R Z T-1 ~> kgSalt m-2 s-1]

  ! applied surface pressure from other component models (e.g., atmos, sea ice, land ice)
  real, pointer, dimension(:,:) :: p_surf_full => NULL()
                !< Pressure at the top ocean interface [R L2 T-2 ~> Pa].
                !! if there is sea-ice, then p_surf_flux is at ice-ocean interface
  real, pointer, dimension(:,:) :: p_surf => NULL()
                !< Pressure at the top ocean interface [R L2 T-2 ~> Pa] as used to drive the ocean model.
                !! If p_surf is limited, p_surf may be smaller than p_surf_full, otherwise they are the same.
  real, pointer, dimension(:,:) :: p_surf_SSH => NULL()
                !< Pressure at the top ocean interface [R L2 T-2 ~> Pa] that is used in corrections to the sea surface
                !! height field that is passed back to the calling routines.
                !! p_surf_SSH may point to p_surf or to p_surf_full.
  logical :: accumulate_p_surf = .false. !< If true, the surface pressure due to the atmosphere
                                 !! and various types of ice needs to be accumulated, and the
                                 !! surface pressure explicitly reset to zero at the driver level
                                 !! when appropriate.

  ! tide related inputs
  real, pointer, dimension(:,:) :: &
    TKE_tidal     => NULL(), & !< tidal energy source driving mixing in bottom boundary layer [R Z3 T-3 ~> W m-2]
    ustar_tidal   => NULL()    !< tidal contribution to bottom ustar [Z T-1 ~> m s-1]

  ! iceberg related inputs
  real, pointer, dimension(:,:) :: &
    ustar_berg => NULL(), &   !< iceberg contribution to top ustar [Z T-1 ~> m s-1].
    area_berg  => NULL(), &   !< area of ocean surface covered by icebergs [m2 m-2]
    mass_berg  => NULL()      !< mass of icebergs [R Z ~> kg m-2]

  ! land ice-shelf related inputs
  real, pointer, dimension(:,:) :: ustar_shelf => NULL()  !< Friction velocity under ice-shelves [Z T-1 ~> m s-1].
                                 !! as computed by the ocean at the previous time step.
  real, pointer, dimension(:,:) :: frac_shelf_h => NULL() !< Fractional ice shelf coverage of
                                 !! h-cells, nondimensional from 0 to 1. This is only
                                 !! associated if ice shelves are enabled, and are
                                 !! exactly 0 away from shelves or on land.
  real, pointer, dimension(:,:) :: iceshelf_melt => NULL() !< Ice shelf melt rate (positive)
                                 !! or freezing (negative) [R Z T-1 ~> kg m-2 s-1]

  ! Scalars set by surface forcing modules
  real :: vPrecGlobalAdj = 0.     !< adjustment to restoring vprec to zero out global net [kg m-2 s-1]
  real :: saltFluxGlobalAdj = 0.  !< adjustment to restoring salt flux to zero out global net [kgSalt m-2 s-1]
  real :: netFWGlobalAdj = 0.     !< adjustment to net fresh water to zero out global net [kg m-2 s-1]
  real :: vPrecGlobalScl = 0.     !< scaling of restoring vprec to zero out global net ( -1..1 ) [nondim]
  real :: saltFluxGlobalScl = 0.  !< scaling of restoring salt flux to zero out global net ( -1..1 ) [nondim]
  real :: netFWGlobalScl = 0.     !< scaling of net fresh water to zero out global net ( -1..1 ) [nondim]

  logical :: fluxes_used = .true. !< If true, all of the heat, salt, and mass
                                  !! fluxes have been applied to the ocean.
  real :: dt_buoy_accum = -1.0    !< The amount of time over which the buoyancy fluxes
                                  !! should be applied [T ~> s].  If negative, this forcing
                                  !! type variable has not yet been inialized.
  logical :: gustless_accum_bug = .true. !< If true, use an incorrect expression in the time
                                  !! average of the gustless wind stress.
  real :: C_p                !< heat capacity of seawater [Q degC-1 ~> J kg-1 degC-1].
                             !! C_p is is the same value as in thermovar_ptrs_type.

  ! passive tracer surface fluxes
  type(coupler_2d_bc_type) :: tr_fluxes !< This structure contains arrays of
     !! of named fields used for passive tracer fluxes.
     !! All arrays in tr_fluxes use the coupler indexing, which has no halos.
     !! This is not a convenient convention, but imposed on MOM6 by the coupler.

  ! For internal error tracking
  integer :: num_msg = 0 !< Number of messages issued about excessive SW penetration
  integer :: max_msg = 2 !< Maximum number of messages to issue about excessive SW penetration

end type forcing

!> Structure that contains pointers to the mechanical forcing at the surface
!! used to drive the liquid ocean simulated by MOM.
!! Data in this type is allocated in the module MOM_surface_forcing.F90,
!! of which there are three versions:  solo, coupled, and ice-shelf.
type, public :: mech_forcing
  ! surface stress components and turbulent velocity scale
  real, pointer, dimension(:,:) :: &
    taux  => NULL(), & !< zonal wind stress [R L Z T-2 ~> Pa]
    tauy  => NULL(), & !< meridional wind stress [R L Z T-2 ~> Pa]
    ustar => NULL(), & !< surface friction velocity scale [Z T-1 ~> m s-1].
    net_mass_src => NULL() !< The net mass source to the ocean [kg m-2 s-1].

  ! applied surface pressure from other component models (e.g., atmos, sea ice, land ice)
  real, pointer, dimension(:,:) :: p_surf_full => NULL()
                !< Pressure at the top ocean interface [R L2 T-2 ~> Pa].
                !! if there is sea-ice, then p_surf_flux is at ice-ocean interface
  real, pointer, dimension(:,:) :: p_surf => NULL()
                !< Pressure at the top ocean interface [R L2 T-2 ~> Pa] as used to drive the ocean model.
                !! If p_surf is limited, p_surf may be smaller than p_surf_full, otherwise they are the same.
  real, pointer, dimension(:,:) :: p_surf_SSH => NULL()
                !< Pressure at the top ocean interface [R L2 T-2 ~> Pa] that is used in corrections
                !! to the sea surface height field that is passed back to the calling routines.
                !! p_surf_SSH may point to p_surf or to p_surf_full.

  ! iceberg related inputs
  real, pointer, dimension(:,:) :: &
    area_berg  => NULL(), &    !< fractional area of ocean surface covered by icebergs [m2 m-2]
    mass_berg  => NULL()       !< mass of icebergs per unit ocean area [R Z ~> kg m-2]

  ! land ice-shelf related inputs
  real, pointer, dimension(:,:) :: frac_shelf_u  => NULL() !< Fractional ice shelf coverage of u-cells,
                !! nondimensional from 0 to 1 [nondim]. This is only associated if ice shelves are enabled,
                !! and is exactly 0 away from shelves or on land.
  real, pointer, dimension(:,:) :: frac_shelf_v  => NULL() !< Fractional ice shelf coverage of v-cells,
                !! nondimensional from 0 to 1 [nondim]. This is only associated if ice shelves are enabled,
                !! and is exactly 0 away from shelves or on land.
  real, pointer, dimension(:,:) :: &
    rigidity_ice_u => NULL(), & !< Depth-integrated lateral viscosity of ice shelves or sea ice at
                                !! u-points [L4 Z-1 T-1 ~> m3 s-1]
    rigidity_ice_v => NULL()    !< Depth-integrated lateral viscosity of ice shelves or sea ice at
                                !! v-points [L4 Z-1 T-1 ~> m3 s-1]
  real :: dt_force_accum = -1.0 !< The amount of time over which the mechanical forcing fluxes
                                !! have been averaged [s].
  logical :: net_mass_src_set = .false. !< If true, an estimate of net_mass_src has been provided.
  logical :: accumulate_p_surf = .false. !< If true, the surface pressure due to the atmosphere
                                !! and various types of ice needs to be accumulated, and the
                                !! surface pressure explicitly reset to zero at the driver level
                                !! when appropriate.
  logical :: accumulate_rigidity = .false. !< If true, the rigidity due to various types of
                                !! ice needs to be accumulated, and the rigidity explicitly
                                !! reset to zero at the driver level when appropriate.

  logical :: initialized = .false. !< This indicates whether the appropriate arrays have been initialized.
end type mech_forcing

!> Structure that defines the id handles for the forcing type
type, public :: forcing_diags

  !>@{ Forcing diagnostic handles
  ! mass flux diagnostic handles
  integer :: id_prcme        = -1, id_evap        = -1
  integer :: id_precip       = -1, id_vprec       = -1
  integer :: id_lprec        = -1, id_fprec       = -1
  integer :: id_lrunoff      = -1, id_frunoff     = -1
  integer :: id_net_massout  = -1, id_net_massin  = -1
  integer :: id_massout_flux = -1, id_massin_flux = -1
  integer :: id_seaice_melt  = -1

  ! global area integrated mass flux diagnostic handles
  integer :: id_total_prcme        = -1, id_total_evap        = -1
  integer :: id_total_precip       = -1, id_total_vprec       = -1
  integer :: id_total_lprec        = -1, id_total_fprec       = -1
  integer :: id_total_lrunoff      = -1, id_total_frunoff     = -1
  integer :: id_total_net_massout  = -1, id_total_net_massin  = -1
  integer :: id_total_seaice_melt  = -1

  ! global area averaged mass flux diagnostic handles
  integer :: id_prcme_ga  = -1, id_evap_ga = -1
  integer :: id_lprec_ga  = -1, id_fprec_ga= -1
  integer :: id_precip_ga = -1, id_vprec_ga= -1

  ! heat flux diagnostic handles
  integer :: id_net_heat_coupler    = -1, id_net_heat_surface      = -1
  integer :: id_sens                = -1, id_LwLatSens             = -1
  integer :: id_sw                  = -1, id_lw                    = -1
  integer :: id_sw_vis              = -1, id_sw_nir                = -1
  integer :: id_lat_evap            = -1, id_lat_frunoff           = -1
  integer :: id_lat                 = -1, id_lat_fprec             = -1
  integer :: id_heat_content_lrunoff= -1, id_heat_content_frunoff  = -1
  integer :: id_heat_content_lprec  = -1, id_heat_content_fprec    = -1
  integer :: id_heat_content_cond   = -1, id_heat_content_surfwater= -1
  integer :: id_heat_content_vprec  = -1, id_heat_content_massout  = -1
  integer :: id_heat_added          = -1, id_heat_content_massin   = -1
  integer :: id_hfrainds            = -1, id_hfrunoffds            = -1
  integer :: id_seaice_melt_heat    = -1, id_heat_content_icemelt  = -1

  ! global area integrated heat flux diagnostic handles
  integer :: id_total_net_heat_coupler    = -1, id_total_net_heat_surface      = -1
  integer :: id_total_sens                = -1, id_total_LwLatSens             = -1
  integer :: id_total_sw                  = -1, id_total_lw                    = -1
  integer :: id_total_lat_evap            = -1, id_total_lat_frunoff           = -1
  integer :: id_total_lat                 = -1, id_total_lat_fprec             = -1
  integer :: id_total_heat_content_lrunoff= -1, id_total_heat_content_frunoff  = -1
  integer :: id_total_heat_content_lprec  = -1, id_total_heat_content_fprec    = -1
  integer :: id_total_heat_content_cond   = -1, id_total_heat_content_surfwater= -1
  integer :: id_total_heat_content_vprec  = -1, id_total_heat_content_massout  = -1
  integer :: id_total_heat_added          = -1, id_total_heat_content_massin   = -1
  integer :: id_total_seaice_melt_heat    = -1, id_total_heat_content_icemelt  = -1

  ! global area averaged heat flux diagnostic handles
  integer :: id_net_heat_coupler_ga = -1, id_net_heat_surface_ga = -1
  integer :: id_sens_ga             = -1, id_LwLatSens_ga        = -1
  integer :: id_sw_ga               = -1, id_lw_ga               = -1
  integer :: id_lat_ga              = -1

  ! salt flux diagnostic handles
  integer :: id_saltflux          = -1
  integer :: id_saltFluxIn        = -1
  integer :: id_saltFluxAdded     = -1

  integer :: id_total_saltflux        = -1
  integer :: id_total_saltFluxIn      = -1
  integer :: id_total_saltFluxAdded   = -1

  integer :: id_vPrecGlobalAdj    = -1
  integer :: id_vPrecGlobalScl    = -1
  integer :: id_saltFluxGlobalAdj = -1
  integer :: id_saltFluxGlobalScl = -1
  integer :: id_netFWGlobalAdj    = -1
  integer :: id_netFWGlobalScl    = -1

  ! momentum flux and forcing diagnostic handles
  integer :: id_taux  = -1
  integer :: id_tauy  = -1
  integer :: id_ustar = -1

  integer :: id_psurf     = -1
  integer :: id_TKE_tidal = -1
  integer :: id_buoy      = -1

  ! iceberg diagnostic handles
  integer :: id_ustar_berg = -1
  integer :: id_area_berg = -1
  integer :: id_mass_berg = -1

  ! Iceberg + Ice shelf diagnostic handles
  integer :: id_ustar_ice_cover = -1
  integer :: id_frac_ice_cover = -1
  !>@}

  integer :: id_clock_forcing = -1 !< CPU clock id

end type forcing_diags

contains

!> This subroutine extracts fluxes from the surface fluxes type. It works on a j-row
!! for optimization purposes. The 2d (i,j) wrapper is the next subroutine below.
!! This routine multiplies fluxes by dt, so that the result is an accumulation of fluxes
!! over a time step.
subroutine extractFluxes1d(G, GV, US, fluxes, optics, nsw, j, dt, &
                  FluxRescaleDepth, useRiverHeatContent, useCalvingHeatContent, &
                  h, T, netMassInOut, netMassOut, net_heat, net_salt, pen_SW_bnd, tv, &
                  aggregate_FW, nonpenSW, netmassInOut_rate, net_Heat_Rate, &
                  net_salt_rate, pen_sw_bnd_Rate, skip_diags)

  type(ocean_grid_type),    intent(in)    :: G              !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV             !< ocean vertical grid structure
  type(unit_scale_type),    intent(in)    :: US             !< A dimensional unit scaling type
  type(forcing),            intent(inout) :: fluxes         !< structure containing pointers to possible
                                                            !! forcing fields. NULL unused fields.
  type(optics_type),        pointer       :: optics         !< pointer to optics
  integer,                  intent(in)    :: nsw            !< number of bands of penetrating SW
  integer,                  intent(in)    :: j              !< j-index to work on
  real,                     intent(in)    :: dt             !< The time step for these fluxes [T ~> s]
  real,                     intent(in)    :: FluxRescaleDepth !< min ocean depth before fluxes
                                                            !! are scaled away [H ~> m or kg m-2]
  logical,                  intent(in)    :: useRiverHeatContent   !< logical for river heat content
  logical,                  intent(in)    :: useCalvingHeatContent !< logical for calving heat content
  real, dimension(SZI_(G),SZK_(G)), &
                            intent(in)    :: h              !< layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZK_(G)), &
                            intent(in)    :: T              !< layer temperatures [degC]
  real, dimension(SZI_(G)), intent(out)   :: netMassInOut   !< net mass flux (non-Bouss) or volume flux
                                                            !! (if Bouss) of water in/out of ocean over
                                                            !! a time step [H ~> m or kg m-2]
  real, dimension(SZI_(G)), intent(out)   :: netMassOut     !< net mass flux (non-Bouss) or volume flux
                                                            !! (if Bouss) of water leaving ocean surface
                                                            !! over a time step [H ~> m or kg m-2].
                                                            !! netMassOut < 0 means mass leaves ocean.
  real, dimension(SZI_(G)), intent(out)   :: net_heat       !< net heat at the surface accumulated over a
                                                            !! time step for coupler + restoring.
                                                            !! Exclude two terms from net_heat:
                                                            !! (1) downwelling (penetrative) SW,
                                                            !! (2) evaporation heat content,
                                                            !! (since do not yet know evap temperature).
                                                            !! [degC H ~> degC m or degC kg m-2].
  real, dimension(SZI_(G)), intent(out)   :: net_salt       !< surface salt flux into the ocean
                                                            !! accumulated over a time step
                                                            !! [ppt H ~> ppt m or ppt kg m-2].
  real, dimension(max(1,nsw),G%isd:G%ied), intent(out) :: pen_SW_bnd !< penetrating SW flux, split into bands.
                                                            !! [degC H ~> degC m or degC kg m-2]
                                                            !! and array size nsw x SZI_(G), where
                                                            !! nsw=number of SW bands in pen_SW_bnd.
                                                            !! This heat flux is not part of net_heat.
  type(thermo_var_ptrs),    intent(inout) :: tv             !< structure containing pointers to available
                                                            !! thermodynamic fields. Used to keep
                                                            !! track of the heat flux associated with net
                                                            !! mass fluxes into the ocean.
  logical,                  intent(in)    :: aggregate_FW   !< For determining how to aggregate forcing.
  real, dimension(SZI_(G)), &
                  optional, intent(out)   :: nonpenSW       !< Non-penetrating SW used in net_heat
                                                            !! [degC H ~> degC m or degC kg m-2].
                                                            !! Summed over SW bands when diagnosing nonpenSW.
  real, dimension(SZI_(G)), &
                  optional, intent(out)   :: net_Heat_rate  !< Rate of net surface heating
                                                            !! [degC H T-1 ~> degC m s-1 or degC kg m-2 s-1].
  real, dimension(SZI_(G)), &
                  optional, intent(out)   :: net_salt_rate  !< Surface salt flux into the ocean
                                                            !! [ppt H T-1 ~> ppt m s-1 or ppt kg m-2 s-1].
  real, dimension(SZI_(G)), &
                  optional, intent(out)   :: netmassInOut_rate !< Rate of net mass flux into the ocean
                                                            !! [H T-1 ~> m s-1 or kg m-2 s-1].
  real, dimension(max(1,nsw),G%isd:G%ied), &
                  optional, intent(out)   :: pen_sw_bnd_rate !< Rate of penetrative shortwave heating
                                                             !! [degC H T-1 ~> degC m s-1 or degC kg m-2 s-1].
  logical,        optional, intent(in)    :: skip_diags      !< If present and true, skip calculating diagnostics

  ! local
  real :: htot(SZI_(G))       ! total ocean depth [H ~> m or kg m-2]
  real :: Pen_sw_tot(SZI_(G)) ! sum across all bands of Pen_SW [degC H ~> degC m or degC kg m-2].
  real :: pen_sw_tot_rate(SZI_(G)) ! Summed rate of shortwave heating across bands
                              ! [degC H T-1 ~> degC m s-1 or degC kg m-2 s-1]
  real :: Ih_limit            ! inverse depth at which surface fluxes start to be limited
                              ! or 0 for no limiting [H-1 ~> m-1 or m2 kg-1]
  real :: scale               ! scale scales away fluxes if depth < FluxRescaleDepth
  real :: I_Cp                ! 1.0 / C_p [degC Q-1 ~> kg degC J-1]
  real :: I_Cp_Hconvert       ! Unit conversion factors divided by the heat capacity
                              ! [degC H R-1 Z-1 Q-1 ~> degC m3 J-1 or kg degC J-1]
  logical :: calculate_diags  ! Indicate to calculate/update diagnostic arrays
  character(len=200) :: mesg
  integer            :: is, ie, nz, i, k, n

  logical :: do_NHR, do_NSR, do_NMIOR, do_PSWBR

  !BGR-Jul 5,2017{
  ! Initializes/sets logicals if 'rates' are requested
  ! These factors are required for legacy reasons
  !  and therefore computed only when optional outputs are requested
  do_NHR = .false.
  do_NSR = .false.
  do_NMIOR = .false.
  do_PSWBR = .false.
  if (present(net_heat_rate)) do_NHR = .true.
  if (present(net_salt_rate)) do_NSR = .true.
  if (present(netmassinout_rate)) do_NMIOR = .true.
  if (present(pen_sw_bnd_rate)) do_PSWBR = .true.
  !}BGR

  Ih_limit  = 0.0 ; if (FluxRescaleDepth > 0.0) Ih_limit  = 1.0 / FluxRescaleDepth
  I_Cp      = 1.0 / fluxes%C_p
  I_Cp_Hconvert = 1.0 / (GV%H_to_RZ * fluxes%C_p)

  is = G%isc ; ie = G%iec ; nz = G%ke

  calculate_diags = .true.
  if (present(skip_diags)) calculate_diags = .not. skip_diags

  ! error checking

  if (nsw > 0) then ; if (nsw /= optics_nbands(optics)) call MOM_error(WARNING, &
    "mismatch in the number of bands of shortwave radiation in MOM_forcing_type extract_fluxes.")
  endif

  if (.not.associated(fluxes%sw)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%sw is not associated.")

  if (.not.associated(fluxes%lw)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%lw is not associated.")

  if (.not.associated(fluxes%latent)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%latent is not associated.")

  if (.not.associated(fluxes%sens)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%sens is not associated.")

  if (.not.associated(fluxes%evap)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: No evaporation defined.")

  if (.not.associated(fluxes%vprec)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%vprec not defined.")

  if ((.not.associated(fluxes%lprec)) .or. &
      (.not.associated(fluxes%fprec))) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: No precipitation defined.")

  do i=is,ie ; htot(i) = h(i,1) ; enddo
  do k=2,nz ; do i=is,ie ; htot(i) = htot(i) + h(i,k) ; enddo ; enddo

  if (nsw >= 1) then
    call extract_optics_slice(optics, j, G, GV, penSW_top=Pen_SW_bnd)
    if (do_PSWBR) call extract_optics_slice(optics, j, G, GV, penSW_top=Pen_SW_bnd_rate)
  endif

  do i=is,ie

    scale = 1.0 ; if ((Ih_limit > 0.0) .and. (htot(i)*Ih_limit < 1.0)) scale = htot(i)*Ih_limit

    ! Convert the penetrating shortwave forcing to (K * H) and reduce fluxes for shallow depths.
    ! (H=m for Bouss, H=kg/m2 for non-Bouss)
    Pen_sw_tot(i) = 0.0
    if (nsw >= 1) then
      do n=1,nsw
        Pen_SW_bnd(n,i) = I_Cp_Hconvert*scale*dt * max(0.0, Pen_SW_bnd(n,i))
        Pen_sw_tot(i)   = Pen_sw_tot(i) + Pen_SW_bnd(n,i)
      enddo
    else
      Pen_SW_bnd(1,i) = 0.0
    endif

    if (do_PSWBR) then  ! Repeat the above code w/ dt=1s for legacy reasons
      pen_sw_tot_rate(i) = 0.0
      if (nsw >= 1) then
        do n=1,nsw
          Pen_SW_bnd_rate(n,i) = I_Cp_Hconvert*scale * max(0.0, Pen_SW_bnd_rate(n,i))
          pen_sw_tot_rate(i) = pen_sw_tot_rate(i) + pen_sw_bnd_rate(n,i)
        enddo
      else
        pen_sw_bnd_rate(1,i) = 0.0
      endif
    endif

    ! net volume/mass of liquid and solid passing through surface boundary fluxes
    netMassInOut(i) = dt * (scale * &
                                   (((((( fluxes%lprec(i,j)        &
                                        + fluxes%fprec(i,j)      )  &
                                        + fluxes%evap(i,j)       )  &
                                        + fluxes%lrunoff(i,j)    )  &
                                        + fluxes%vprec(i,j)      )  &
                                        + fluxes%seaice_melt(i,j))  &
                                        + fluxes%frunoff(i,j)    ))

    if (do_NMIOr) then  ! Repeat the above code without multiplying by a timestep for legacy reasons
      netMassInOut_rate(i) = (scale * &
                                   (((((( fluxes%lprec(i,j)      &
                                        + fluxes%fprec(i,j)      )  &
                                        + fluxes%evap(i,j)       )  &
                                        + fluxes%lrunoff(i,j)    )  &
                                        + fluxes%vprec(i,j)      )  &
                                        + fluxes%seaice_melt(i,j))  &
                                        + fluxes%frunoff(i,j)   ))
    endif

    ! smg:
    ! for non-Bouss, we add/remove salt mass to total ocean mass. to conserve
    ! total salt mass ocean+ice, the sea ice model must lose mass when salt mass
    ! is added to the ocean, which may still need to be coded.  Not that the units
    ! of netMassInOut are still kg_m2, so no conversion to H should occur yet.
    if (.not.GV%Boussinesq .and. associated(fluxes%salt_flux)) then
      netMassInOut(i) = netMassInOut(i) + dt * (scale * fluxes%salt_flux(i,j))
      if (do_NMIOr) netMassInOut_rate(i) = netMassInOut_rate(i) + &
                                               (scale * fluxes%salt_flux(i,j))
    endif

    ! net volume/mass of water leaving the ocean.
    ! check that fluxes are < 0, which means mass is indeed leaving.
    netMassOut(i) = 0.0

    ! evap > 0 means condensating water is added into ocean.
    ! evap < 0 means evaporation of water from the ocean, in
    ! which case heat_content_evap is computed in MOM_diabatic_driver.F90
    if (fluxes%evap(i,j) < 0.0) netMassOut(i) = netMassOut(i) + fluxes%evap(i,j)
  !   if (associated(fluxes%heat_content_cond)) fluxes%heat_content_cond(i,j) = 0.0 !??? --AJA

    ! lprec < 0 means sea ice formation taking water from the ocean.
    ! smg: we should split the ice melt/formation from the lprec
    if (fluxes%lprec(i,j) < 0.0) netMassOut(i) = netMassOut(i) + fluxes%lprec(i,j)

    ! seaice_melt < 0 means sea ice formation taking water from the ocean.
    if (fluxes%seaice_melt(i,j) < 0.0) netMassOut(i) = netMassOut(i) + fluxes%seaice_melt(i,j)

    ! vprec < 0 means virtual evaporation arising from surface salinity restoring,
    ! in which case heat_content_vprec is computed in MOM_diabatic_driver.F90.
    if (fluxes%vprec(i,j) < 0.0) netMassOut(i) = netMassOut(i) + fluxes%vprec(i,j)

    netMassOut(i) = dt * scale * netMassOut(i)

    ! convert to H units (Bouss=meter or non-Bouss=kg/m^2)
    netMassInOut(i) = GV%RZ_to_H * netMassInOut(i)
    if (do_NMIOr) netMassInOut_rate(i) = GV%RZ_to_H * netMassInOut_rate(i)
    netMassOut(i)   = GV%RZ_to_H * netMassOut(i)

    ! surface heat fluxes from radiation and turbulent fluxes (K * H)
    ! (H=m for Bouss, H=kg/m2 for non-Bouss)

    ! CIME provides heat flux from snow&ice melt (seaice_melt_heat), so this is added below
    if (associated(fluxes%seaice_melt_heat)) then
      net_heat(i) = scale * dt * I_Cp_Hconvert * &
                    ( fluxes%sw(i,j) + (((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)) + &
                      fluxes%seaice_melt_heat(i,j)) )
      !Repeats above code w/ dt=1. for legacy reason
      if (do_NHR)  net_heat_rate(i) = scale * I_Cp_Hconvert * &
           ( fluxes%sw(i,j) + (((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)) + &
             fluxes%seaice_melt_heat(i,j)))
    else
      net_heat(i) = scale * dt * I_Cp_Hconvert * &
                    ( fluxes%sw(i,j) + ((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)) )
      !Repeats above code w/ dt=1. for legacy reason
      if (do_NHR)  net_heat_rate(i) = scale * I_Cp_Hconvert * &
           ( fluxes%sw(i,j) + ((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)) )
    endif

    ! Add heat flux from surface damping (restoring) (K * H) or flux adjustments.
    if (associated(fluxes%heat_added)) then
      net_heat(i) = net_heat(i) + (scale * (dt * I_Cp_Hconvert)) * fluxes%heat_added(i,j)
      if (do_NHR) net_heat_rate(i) = net_heat_rate(i) + (scale * I_Cp_Hconvert) * fluxes%heat_added(i,j)
    endif

    ! Add explicit heat flux for runoff (which is part of the ice-ocean boundary
    ! flux type). Runoff is otherwise added with a temperature of SST.
    if (useRiverHeatContent) then
      ! remove lrunoff*SST here, to counteract its addition elsewhere
      net_heat(i) = (net_heat(i) + (scale*(dt * I_Cp_Hconvert)) * fluxes%heat_content_lrunoff(i,j)) - &
                     (GV%RZ_to_H * (scale * dt)) * fluxes%lrunoff(i,j) * T(i,1)
      !BGR-Jul 5, 2017{
      !Intentionally neglect the following contribution to rate for legacy reasons.
      !if (do_NHR) net_heat_rate(i) = (net_heat_rate(i) + (scale*I_Cp_Hconvert) * fluxes%heat_content_lrunoff(i,j)) - &
      !               (GV%RZ_to_H * (scale)) * fluxes%lrunoff(i,j) * T(i,1)
      !}BGR
      if (calculate_diags .and. associated(tv%TempxPmE)) then
        tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + (scale * dt) * &
            (I_Cp*fluxes%heat_content_lrunoff(i,j) - fluxes%lrunoff(i,j)*T(i,1))
      endif
    endif

    ! Add explicit heat flux for calving (which is part of the ice-ocean boundary
    ! flux type). Calving is otherwise added with a temperature of SST.
    if (useCalvingHeatContent) then
      ! remove frunoff*SST here, to counteract its addition elsewhere
      net_heat(i) = net_heat(i) + (scale*(dt * I_Cp_Hconvert)) * fluxes%heat_content_frunoff(i,j) - &
                    (GV%RZ_to_H * (scale * dt)) * fluxes%frunoff(i,j) * T(i,1)
      !BGR-Jul 5, 2017{
      !Intentionally neglect the following contribution to rate for legacy reasons.
!      if (do_NHR) net_heat_rate(i) = net_heat_rate(i) + (scale*I_Cp_Hconvert) * fluxes%heat_content_frunoff(i,j) - &
!                    (GV%RZ_to_H * scale) * fluxes%frunoff(i,j) * T(i,1)
      !}BGR
      if (calculate_diags .and. associated(tv%TempxPmE)) then
        tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + (scale * dt) * &
            (I_Cp*fluxes%heat_content_frunoff(i,j) - fluxes%frunoff(i,j)*T(i,1))
      endif
    endif

! smg: new code
    ! add heat from all terms that may add mass to the ocean (K * H).
    ! if evap, lprec, or vprec < 0, then compute their heat content
    ! inside MOM_diabatic_driver.F90 and fill in fluxes%heat_content_massout.
    ! we do so since we do not here know the temperature
    ! of water leaving the ocean, as it could be leaving from more than
    ! one layer of the upper ocean in the case of very thin layers.
    ! When evap, lprec, or vprec > 0, then we know their heat content here
    ! via settings from inside of the appropriate config_src driver files.
!    if (associated(fluxes%heat_content_lprec)) then
!      net_heat(i) = net_heat(i) + scale * dt * I_Cp_Hconvert * &
!     (fluxes%heat_content_lprec(i,j)    + (fluxes%heat_content_fprec(i,j)   + &
!     (fluxes%heat_content_lrunoff(i,j)  + (fluxes%heat_content_frunoff(i,j) + &
!     (fluxes%heat_content_cond(i,j)     +  fluxes%heat_content_vprec(i,j))))))
!    endif

    if (fluxes%num_msg < fluxes%max_msg) then
      if (Pen_SW_tot(i) > 1.000001 * I_Cp_Hconvert*scale*dt*fluxes%sw(i,j)) then
        fluxes%num_msg = fluxes%num_msg + 1
        write(mesg,'("Penetrating shortwave of ",1pe17.10, &
                    &" exceeds total shortwave of ",1pe17.10,&
                    &" at ",1pg11.4,"E, "1pg11.4,"N.")') &
               Pen_SW_tot(i), I_Cp_Hconvert*scale*dt * fluxes%sw(i,j), &
               G%geoLonT(i,j), G%geoLatT(i,j)
        call MOM_error(WARNING,mesg)
      endif
    endif

    ! remove penetrative portion of the SW that is NOT absorbed within a
    ! tiny layer at the top of the ocean.
    net_heat(i) = net_heat(i) - Pen_SW_tot(i)
    !Repeat above code for 'rate' term
    if (do_NHR) net_heat_rate(i) = net_heat_rate(i) - Pen_SW_tot_rate(i)

    ! diagnose non-downwelling SW
    if (present(nonPenSW)) then
      nonPenSW(i) = scale * dt * I_Cp_Hconvert * fluxes%sw(i,j) - Pen_SW_tot(i)
    endif

    ! Salt fluxes
    Net_salt(i) = 0.0
    if (do_NSR) Net_salt_rate(i) = 0.0
    ! Convert salt_flux from kg (salt)/(m^2 * s) to
    ! Boussinesq: (ppt * m)
    ! non-Bouss:  (g/m^2)
    if (associated(fluxes%salt_flux)) then
      Net_salt(i) = (scale * dt * (1000.0 * fluxes%salt_flux(i,j))) * GV%RZ_to_H
      !Repeat above code for 'rate' term
      if (do_NSR) Net_salt_rate(i) = (scale * 1. * (1000.0 * fluxes%salt_flux(i,j))) * GV%RZ_to_H
    endif

    ! Diagnostics follow...
    if (calculate_diags) then

      ! Store Net_salt for unknown reason?
      if (associated(fluxes%salt_flux)) then
        ! This seems like a bad idea to me. -RWH
        if (calculate_diags) fluxes%netSalt(i,j) = US%kg_m2s_to_RZ_T*Net_salt(i)
      endif

      ! Initialize heat_content_massin that is diagnosed in mixedlayer_convection or
      ! applyBoundaryFluxes such that the meaning is as the sum of all incoming components.
      if (associated(fluxes%heat_content_massin))  then
        if (aggregate_FW) then
          if (netMassInOut(i) > 0.0) then ! net is "in"
            fluxes%heat_content_massin(i,j) = -fluxes%C_p * netMassOut(i) * T(i,1) * GV%H_to_RZ / dt
          else ! net is "out"
            fluxes%heat_content_massin(i,j) = fluxes%C_p * ( netMassInout(i) - netMassOut(i) ) * &
                                               T(i,1) * GV%H_to_RZ / dt
          endif
        else
          fluxes%heat_content_massin(i,j) = 0.
        endif
      endif

      ! Initialize heat_content_massout that is diagnosed in mixedlayer_convection or
      ! applyBoundaryFluxes such that the meaning is as the sum of all outgoing components.
      if (associated(fluxes%heat_content_massout)) then
        if (aggregate_FW) then
          if (netMassInOut(i) > 0.0) then ! net is "in"
            fluxes%heat_content_massout(i,j) = fluxes%C_p * netMassOut(i) * T(i,1) * GV%H_to_RZ / dt
          else ! net is "out"
            fluxes%heat_content_massout(i,j) = -fluxes%C_p * ( netMassInout(i) - netMassOut(i) ) * &
                                                T(i,1) * GV%H_to_RZ / dt
          endif
        else
          fluxes%heat_content_massout(i,j) = 0.0
        endif
      endif

      ! smg: we should remove sea ice melt from lprec!!!
      ! fluxes%lprec > 0 means ocean gains mass via liquid precipitation and/or sea ice melt.
      ! When atmosphere does not provide heat of this precipitation, the ocean assumes
      ! it enters the ocean at the SST.
      ! fluxes%lprec < 0 means ocean loses mass via sea ice formation. As we do not yet know
      ! the layer at which this mass is removed, we cannot compute it heat content. We must
      ! wait until MOM_diabatic_driver.F90.
      if (associated(fluxes%heat_content_lprec)) then
        if (fluxes%lprec(i,j) > 0.0) then
          fluxes%heat_content_lprec(i,j) = fluxes%C_p*fluxes%lprec(i,j)*T(i,1)
        else
          fluxes%heat_content_lprec(i,j) = 0.0
        endif
      endif

      ! fprec SHOULD enter ocean at 0degC if atmos model does not provide fprec heat content.
      ! However, we need to adjust netHeat above to reflect the difference between 0decC and SST
      ! and until we do so fprec is treated like lprec and enters at SST. -AJA
      if (associated(fluxes%heat_content_fprec)) then
        if (fluxes%fprec(i,j) > 0.0) then
          fluxes%heat_content_fprec(i,j) = fluxes%C_p*fluxes%fprec(i,j)*T(i,1)
        else
          fluxes%heat_content_fprec(i,j) = 0.0
        endif
      endif

      ! Following lprec and fprec, water flux due to sea ice melt (seaice_melt) enters at SST - GMM
      if (associated(fluxes%heat_content_icemelt)) then
        if (fluxes%seaice_melt(i,j) > 0.0) then
          fluxes%heat_content_icemelt(i,j) = fluxes%C_p*fluxes%seaice_melt(i,j)*T(i,1)
        else
          fluxes%heat_content_icemelt(i,j) = 0.0
        endif
      endif

      ! virtual precip associated with salinity restoring
      ! vprec > 0 means add water to ocean, assumed to be at SST
      ! vprec < 0 means remove water from ocean; set heat_content_vprec in MOM_diabatic_driver.F90
      if (associated(fluxes%heat_content_vprec)) then
        if (fluxes%vprec(i,j) > 0.0) then
          fluxes%heat_content_vprec(i,j) = fluxes%C_p*fluxes%vprec(i,j)*T(i,1)
        else
          fluxes%heat_content_vprec(i,j) = 0.0
        endif
      endif

      ! fluxes%evap < 0 means ocean loses mass due to evaporation.
      ! Evaporation leaves ocean surface at a temperature that has yet to be determined,
      ! since we do not know the precise layer that the water evaporates.  We therefore
      ! compute fluxes%heat_content_massout at the relevant point inside MOM_diabatic_driver.F90.
      ! fluxes%evap > 0 means ocean gains moisture via condensation.
      ! Condensation is assumed to drop into the ocean at the SST, just like lprec.
      if (associated(fluxes%heat_content_cond)) then
        if (fluxes%evap(i,j) > 0.0) then
          fluxes%heat_content_cond(i,j) = fluxes%C_p*fluxes%evap(i,j)*T(i,1)
        else
          fluxes%heat_content_cond(i,j) = 0.0
        endif
      endif

      ! Liquid runoff enters ocean at SST if land model does not provide runoff heat content.
      if (.not. useRiverHeatContent) then
        if (associated(fluxes%lrunoff) .and. associated(fluxes%heat_content_lrunoff)) then
          fluxes%heat_content_lrunoff(i,j) = fluxes%C_p*fluxes%lrunoff(i,j)*T(i,1)
        endif
      endif

      ! Icebergs enter ocean at SST if land model does not provide calving heat content.
      if (.not. useCalvingHeatContent) then
        if (associated(fluxes%frunoff) .and. associated(fluxes%heat_content_frunoff)) then
          fluxes%heat_content_frunoff(i,j) = fluxes%C_p*fluxes%frunoff(i,j)*T(i,1)
        endif
      endif

    endif ! calculate_diags

  enddo ! i-loop

end subroutine extractFluxes1d


!> 2d wrapper for 1d extract fluxes from surface fluxes type.
!! This subroutine extracts fluxes from the surface fluxes type. It multiplies the
!! fluxes by dt, so that the result is an accumulation of the fluxes over a time step.
subroutine extractFluxes2d(G, GV, US, fluxes, optics, nsw, dt, FluxRescaleDepth, &
                           useRiverHeatContent, useCalvingHeatContent, h, T, &
                           netMassInOut, netMassOut, net_heat, Net_salt, Pen_SW_bnd, tv, &
                           aggregate_FW)

  type(ocean_grid_type),            intent(in)    :: G              !< ocean grid structure
  type(verticalGrid_type),          intent(in)    :: GV             !< ocean vertical grid structure
  type(unit_scale_type),            intent(in)    :: US             !< A dimensional unit scaling type
  type(forcing),                    intent(inout) :: fluxes         !< structure containing pointers to forcing.
  type(optics_type),                pointer       :: optics         !< pointer to optics
  integer,                          intent(in)    :: nsw            !< number of bands of penetrating SW
  real,                             intent(in)    :: dt             !< The time step for these fluxes [T ~> s]
  real,                             intent(in)    :: FluxRescaleDepth !< min ocean depth before fluxes
                                                                    !! are scaled away [H ~> m or kg m-2]
  logical,                          intent(in)    :: useRiverHeatContent   !< logical for river heat content
  logical,                          intent(in)    :: useCalvingHeatContent !< logical for calving heat content
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                                    intent(in)    :: h              !< layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                                    intent(in)    :: T              !< layer temperatures [degC]
  real, dimension(SZI_(G),SZJ_(G)), intent(out)   :: netMassInOut   !< net mass flux (non-Bouss) or volume flux
                                                                    !! (if Bouss) of water in/out of ocean over
                                                                    !! a time step [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out)   :: netMassOut     !< net mass flux (non-Bouss) or volume flux
                                                                    !! (if Bouss) of water leaving ocean surface
                                                                    !! over a time step [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)), intent(out)   :: net_heat       !< net heat at the surface accumulated over a
                                                                    !! time step associated with coupler + restore.
                                                                    !! Exclude two terms from net_heat:
                                                                    !! (1) downwelling (penetrative) SW,
                                                                    !! (2) evaporation heat content,
                                                                    !! (since do not yet know temperature of evap).
                                                                    !! [degC H ~> degC m or degC kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out)   :: net_salt       !< surface salt flux into the ocean accumulated
                                                                    !! over a time step [ppt H ~> ppt m or ppt kg m-2]
  real, dimension(max(1,nsw),G%isd:G%ied,G%jsd:G%jed), intent(out) :: pen_SW_bnd !< penetrating SW flux, by frequency
                                                                    !! band [degC H ~> degC m or degC kg m-2] with array
                                                                    !! size nsw x SZI_(G), where nsw=number of SW bands
                                                                    !! in pen_SW_bnd. This heat flux is not in net_heat.
  type(thermo_var_ptrs),            intent(inout) :: tv             !< structure containing pointers to available
                                                                    !! thermodynamic fields. Here it is used to keep
                                                                    !! track of the heat flux associated with net
                                                                    !! mass fluxes into the ocean.
  logical,                          intent(in)    :: aggregate_FW   !< For determining how to aggregate the forcing.

  integer :: j
  !$OMP parallel do default(shared)
  do j=G%jsc, G%jec
    call extractFluxes1d(G, GV, US, fluxes, optics, nsw, j, dt, &
            FluxRescaleDepth, useRiverHeatContent, useCalvingHeatContent,&
            h(:,j,:), T(:,j,:), netMassInOut(:,j), netMassOut(:,j),              &
            net_heat(:,j), net_salt(:,j), pen_SW_bnd(:,:,j), tv, aggregate_FW)
  enddo

end subroutine extractFluxes2d


!> This routine calculates surface buoyancy flux by adding up the heat, FW & salt fluxes.
!! These are actual fluxes, with units of stuff per time. Setting dt=1 in the call to
!! extractFluxes routine allows us to get "stuf per time" rather than the time integrated
!! fluxes needed in other routines that call extractFluxes.
subroutine calculateBuoyancyFlux1d(G, GV, US, fluxes, optics, nsw, h, Temp, Salt, tv, j, &
                                   buoyancyFlux, netHeatMinusSW, netSalt, skip_diags)
  type(ocean_grid_type),                    intent(in)    :: G              !< ocean grid
  type(verticalGrid_type),                  intent(in)    :: GV             !< ocean vertical grid structure
  type(unit_scale_type),                    intent(in)    :: US             !< A dimensional unit scaling type
  type(forcing),                            intent(inout) :: fluxes         !< surface fluxes
  type(optics_type),                        pointer       :: optics         !< penetrating SW optics
  integer,                                  intent(in)    :: nsw            !< The number of frequency bands of
                                                                            !! penetrating shortwave radiation
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h              !< layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: Temp           !< prognostic temp [degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: Salt           !< salinity [ppt]
  type(thermo_var_ptrs),                    intent(inout) :: tv             !< thermodynamics type
  integer,                                  intent(in)    :: j              !< j-row to work on
  real, dimension(SZI_(G),SZK_(G)+1),       intent(inout) :: buoyancyFlux   !< buoyancy fluxes [L2 T-3 ~> m2 s-3]
  real, dimension(SZI_(G)),                 intent(inout) :: netHeatMinusSW !< surf Heat flux
                                                                      !! [degC H s-1 ~> degC m s-1 or degC kg m-2 s-1]
  real, dimension(SZI_(G)),                 intent(inout) :: netSalt        !< surf salt flux
                                                                      !! [ppt H s-1 ~> ppt m s-1 or ppt kg m-2 s-1]
  logical,                        optional, intent(in)    :: skip_diags     !< If present and true, skip calculating
                                                                            !! diagnostics inside extractFluxes1d()
  ! local variables
  integer                               :: k
  real, parameter                       :: dt = 1.    ! to return a rate from extractFluxes1d
  real, dimension(SZI_(G))              :: netH       ! net FW flux [H s-1 ~> m s-1 or kg m-2 s-1]
  real, dimension(SZI_(G))              :: netEvap    ! net FW flux leaving ocean via evaporation
                                                      ! [H s-1 ~> m s-1 or kg m-2 s-1]
  real, dimension(SZI_(G))              :: netHeat    ! net temp flux [degC H s-1 ~> degC m s-2 or degC kg m-2 s-1]
  real, dimension(max(nsw,1), SZI_(G))  :: penSWbnd   ! penetrating SW radiation by band
                                                      ! [degC H ~> degC m or degC kg m-2]
  real, dimension(SZI_(G))              :: pressure   ! pressure at the surface [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G))              :: dRhodT     ! density partial derivative wrt temp [R degC-1 ~> kg m-3 degC-1]
  real, dimension(SZI_(G))              :: dRhodS     ! density partial derivative wrt saln [R ppt-1 ~> kg m-3 ppt-1]
  real, dimension(SZI_(G),SZK_(G)+1)    :: netPen     ! The net penetrating shortwave radiation at each level
                                                      ! [degC H ~> degC m or degC kg m-2]

  logical :: useRiverHeatContent
  logical :: useCalvingHeatContent
  real    :: depthBeforeScalingFluxes  ! A depth scale [H ~> m or kg m-2]
  real    :: GoRho ! The gravitational acceleration divided by mean density times some
                   ! unit conversion factors [L2 H-1 s R-1 T-3 ~> m4 kg-1 s-2 or m7 kg-2 s-2]
  real    :: H_limit_fluxes            ! Another depth scale [H ~> m or kg m-2]
  integer :: i

  !  smg: what do we do when have heat fluxes from calving and river?
  useRiverHeatContent   = .False.
  useCalvingHeatContent = .False.

  depthBeforeScalingFluxes = max( GV%Angstrom_H, 1.e-30*GV%m_to_H )
  pressure(:) = 0.
  if (associated(tv%p_surf)) then ; do i=G%isc,G%iec ; pressure(i) = tv%p_surf(i,j) ; enddo ; endif
  GoRho       = (GV%g_Earth * GV%H_to_Z*US%T_to_s) / GV%Rho0

  H_limit_fluxes = depthBeforeScalingFluxes

  ! The surface forcing is contained in the fluxes type.
  ! We aggregate the thermodynamic forcing for a time step into the following:
  ! netH       = water added/removed via surface fluxes [H s-1 ~> m s-1 or kg m-2 s-1]
  ! netHeat    = heat via surface fluxes [degC H s-1 ~> degC m s-1 or degC kg m-2 s-1]
  ! netSalt    = salt via surface fluxes [ppt H s-1 ~> ppt m s-1 or gSalt m-2 s-1]
  ! Note that unlike other calls to extractFLuxes1d() that return the time-integrated flux
  ! this call returns the rate because dt=1
  call extractFluxes1d(G, GV, US, fluxes, optics, nsw, j, dt*US%s_to_T,               &
                depthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
                h(:,j,:), Temp(:,j,:), netH, netEvap, netHeatMinusSW,                 &
                netSalt, penSWbnd, tv, .false., skip_diags=skip_diags)

  ! Sum over bands and attenuate as a function of depth
  ! netPen is the netSW as a function of depth
  call sumSWoverBands(G, GV, US, h(:,j,:), optics_nbands(optics), optics, j, dt*US%s_to_T, &
                      H_limit_fluxes, .true., penSWbnd, netPen)

  ! Density derivatives
  call calculate_density_derivs(Temp(:,j,1), Salt(:,j,1), pressure, dRhodT, dRhodS, &
                                tv%eqn_of_state, EOS_domain(G%HI))

  ! Adjust netSalt to reflect dilution effect of FW flux
  netSalt(G%isc:G%iec) = netSalt(G%isc:G%iec) - Salt(G%isc:G%iec,j,1) * netH(G%isc:G%iec) ! ppt H/s

  ! Add in the SW heating for purposes of calculating the net
  ! surface buoyancy flux affecting the top layer.
  !netHeat(:) = netHeatMinusSW(:) + sum( penSWbnd, dim=1 )
  netHeat(G%isc:G%iec) = netHeatMinusSW(G%isc:G%iec) + netPen(G%isc:G%iec,1) ! K H/s

  ! Convert to a buoyancy flux, excluding penetrating SW heating
  buoyancyFlux(G%isc:G%iec,1) = - GoRho * ( dRhodS(G%isc:G%iec) * netSalt(G%isc:G%iec) + &
                                             dRhodT(G%isc:G%iec) * netHeat(G%isc:G%iec) ) ! [L2 T-3 ~> m2 s-3]
  ! We also have a penetrative buoyancy flux associated with penetrative SW
  do k=2, G%ke+1
    buoyancyFlux(G%isc:G%iec,k) = - GoRho * ( dRhodT(G%isc:G%iec) * netPen(G%isc:G%iec,k) ) ! [L2 T-3 ~> m2 s-3]
  enddo

end subroutine calculateBuoyancyFlux1d


!> Calculates surface buoyancy flux by adding up the heat, FW and salt fluxes,
!! for 2d arrays.  This is a wrapper for calculateBuoyancyFlux1d.
subroutine calculateBuoyancyFlux2d(G, GV, US, fluxes, optics, h, Temp, Salt, tv, &
                                   buoyancyFlux, netHeatMinusSW, netSalt, skip_diags)
  type(ocean_grid_type),                      intent(in)    :: G      !< ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV     !< ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  type(forcing),                              intent(inout) :: fluxes !< surface fluxes
  type(optics_type),                          pointer       :: optics !< SW ocean optics
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: h      !< layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: Temp   !< temperature [degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: Salt   !< salinity [ppt]
  type(thermo_var_ptrs),                      intent(inout) :: tv     !< thermodynamics type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: buoyancyFlux   !< buoyancy fluxes [L2 T-3 ~> m2 s-3]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(inout) :: netHeatMinusSW !< surf temp flux
                                                                              !! [degC H ~> degC m or degC kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(inout) :: netSalt        !< surf salt flux
                                                                              !! [ppt H ~> ppt m or ppt kg m-2]
  logical, optional,                          intent(in)    :: skip_diags     !< If present and true, skip calculating
                                                                              !! diagnostics inside extractFluxes1d()
  ! local variables
  real, dimension( SZI_(G) ) :: netT ! net temperature flux [degC H s-1 ~> degC m s-2 or degC kg m-2 s-1]
  real, dimension( SZI_(G) ) :: netS ! net saln flux !! [ppt H s-1 ~> ppt m s-1 or ppt kg m-2 s-1]
  integer :: j

  netT(G%isc:G%iec) = 0. ; netS(G%isc:G%iec) = 0.

  !$OMP parallel do default(shared) firstprivate(netT,netS)
  do j=G%jsc,G%jec
    call calculateBuoyancyFlux1d(G, GV, US, fluxes, optics, optics_nbands(optics), h, Temp, Salt, &
                                 tv, j, buoyancyFlux(:,j,:), netT, netS, skip_diags=skip_diags)
    if (present(netHeatMinusSW)) netHeatMinusSW(G%isc:G%iec,j) = netT(G%isc:G%iec)
    if (present(netSalt)) netSalt(G%isc:G%iec,j) = netS(G%isc:G%iec)
  enddo

end subroutine calculateBuoyancyFlux2d


!> Write out chksums for thermodynamic fluxes.
subroutine MOM_forcing_chksum(mesg, fluxes, G, US, haloshift)
  character(len=*),        intent(in) :: mesg      !< message
  type(forcing),           intent(in) :: fluxes    !< A structure containing thermodynamic forcing fields
  type(ocean_grid_type),   intent(in) :: G         !< grid type
  type(unit_scale_type),   intent(in) :: US        !< A dimensional unit scaling type
  integer, optional,       intent(in) :: haloshift !< shift in halo

  integer :: is, ie, js, je, nz, hshift
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  hshift = 1 ; if (present(haloshift)) hshift = haloshift

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  if (associated(fluxes%ustar)) &
    call hchksum(fluxes%ustar, mesg//" fluxes%ustar", G%HI, haloshift=hshift, scale=US%Z_to_m*US%s_to_T)
  if (associated(fluxes%buoy)) &
    call hchksum(fluxes%buoy, mesg//" fluxes%buoy ", G%HI, haloshift=hshift, scale=US%L_to_m**2*US%s_to_T**3)
  if (associated(fluxes%sw)) &
    call hchksum(fluxes%sw, mesg//" fluxes%sw", G%HI, haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%sw_vis_dir)) &
    call hchksum(fluxes%sw_vis_dir, mesg//" fluxes%sw_vis_dir", G%HI, haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%sw_vis_dif)) &
    call hchksum(fluxes%sw_vis_dif, mesg//" fluxes%sw_vis_dif", G%HI, haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%sw_nir_dir)) &
    call hchksum(fluxes%sw_nir_dir, mesg//" fluxes%sw_nir_dir", G%HI, haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%sw_nir_dif)) &
    call hchksum(fluxes%sw_nir_dif, mesg//" fluxes%sw_nir_dif", G%HI, haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%lw)) &
    call hchksum(fluxes%lw, mesg//" fluxes%lw", G%HI, haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%latent)) &
    call hchksum(fluxes%latent, mesg//" fluxes%latent", G%HI, haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%latent_evap_diag)) &
    call hchksum(fluxes%latent_evap_diag, mesg//" fluxes%latent_evap_diag", G%HI, &
                 haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%latent_fprec_diag)) &
    call hchksum(fluxes%latent_fprec_diag, mesg//" fluxes%latent_fprec_diag", G%HI, &
                 haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%latent_frunoff_diag)) &
    call hchksum(fluxes%latent_frunoff_diag, mesg//" fluxes%latent_frunoff_diag", G%HI, &
                 haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%sens)) &
    call hchksum(fluxes%sens, mesg//" fluxes%sens", G%HI, haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%evap)) &
    call hchksum(fluxes%evap, mesg//" fluxes%evap", G%HI, haloshift=hshift, scale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%lprec)) &
    call hchksum(fluxes%lprec, mesg//" fluxes%lprec", G%HI, haloshift=hshift, scale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%fprec)) &
    call hchksum(fluxes%fprec, mesg//" fluxes%fprec", G%HI, haloshift=hshift, scale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%vprec)) &
    call hchksum(fluxes%vprec, mesg//" fluxes%vprec", G%HI, haloshift=hshift, scale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%seaice_melt)) &
    call hchksum(fluxes%seaice_melt, mesg//" fluxes%seaice_melt", G%HI, haloshift=hshift, scale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%seaice_melt_heat)) &
    call hchksum(fluxes%seaice_melt_heat, mesg//" fluxes%seaice_melt_heat", G%HI, &
                 haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%p_surf)) &
    call hchksum(fluxes%p_surf, mesg//" fluxes%p_surf", G%HI, haloshift=hshift , scale=US%RL2_T2_to_Pa)
  if (associated(fluxes%salt_flux)) &
    call hchksum(fluxes%salt_flux, mesg//" fluxes%salt_flux", G%HI, haloshift=hshift, scale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%TKE_tidal)) &
    call hchksum(fluxes%TKE_tidal, mesg//" fluxes%TKE_tidal", G%HI, haloshift=hshift, &
                 scale=US%RZ3_T3_to_W_m2)
  if (associated(fluxes%ustar_tidal)) &
    call hchksum(fluxes%ustar_tidal, mesg//" fluxes%ustar_tidal", G%HI, haloshift=hshift, scale=US%Z_to_m*US%s_to_T)
  if (associated(fluxes%lrunoff)) &
    call hchksum(fluxes%lrunoff, mesg//" fluxes%lrunoff", G%HI, haloshift=hshift, scale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%frunoff)) &
    call hchksum(fluxes%frunoff, mesg//" fluxes%frunoff", G%HI, haloshift=hshift, scale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%heat_content_lrunoff)) &
    call hchksum(fluxes%heat_content_lrunoff, mesg//" fluxes%heat_content_lrunoff", G%HI, &
                 haloshift=hshift, scale=US%RZ_T_to_kg_m2s)
  if (associated(fluxes%heat_content_frunoff)) &
    call hchksum(fluxes%heat_content_frunoff, mesg//" fluxes%heat_content_frunoff", G%HI, &
                 haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_lprec)) &
    call hchksum(fluxes%heat_content_lprec, mesg//" fluxes%heat_content_lprec", G%HI,  &
                 haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_fprec)) &
    call hchksum(fluxes%heat_content_fprec, mesg//" fluxes%heat_content_fprec", G%HI, &
                 haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_icemelt)) &
    call hchksum(fluxes%heat_content_icemelt, mesg//" fluxes%heat_content_icemelt", G%HI, &
                 haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_cond)) &
    call hchksum(fluxes%heat_content_cond, mesg//" fluxes%heat_content_cond", G%HI, &
                 haloshift=hshift, scale=US%QRZ_T_to_W_m2)
  if (associated(fluxes%heat_content_massout)) &
    call hchksum(fluxes%heat_content_massout, mesg//" fluxes%heat_content_massout", G%HI, &
                 haloshift=hshift, scale=US%QRZ_T_to_W_m2)
end subroutine MOM_forcing_chksum

!> Write out chksums for the driving mechanical forces.
subroutine MOM_mech_forcing_chksum(mesg, forces, G, US, haloshift)
  character(len=*),        intent(in) :: mesg      !< message
  type(mech_forcing),      intent(in) :: forces    !< A structure with the driving mechanical forces
  type(ocean_grid_type),   intent(in) :: G         !< grid type
  type(unit_scale_type),   intent(in) :: US        !< A dimensional unit scaling type
  integer, optional,       intent(in) :: haloshift !< shift in halo

  integer :: is, ie, js, je, nz, hshift
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  hshift=1; if (present(haloshift)) hshift=haloshift

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  if (associated(forces%taux) .and. associated(forces%tauy)) &
    call uvchksum(mesg//" forces%tau[xy]", forces%taux, forces%tauy, G%HI, &
                  haloshift=hshift, symmetric=.true., scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
  if (associated(forces%p_surf)) &
    call hchksum(forces%p_surf, mesg//" forces%p_surf", G%HI, haloshift=hshift, scale=US%RL2_T2_to_Pa)
  if (associated(forces%ustar)) &
    call hchksum(forces%ustar, mesg//" forces%ustar", G%HI, haloshift=hshift, scale=US%Z_to_m*US%s_to_T)
  if (associated(forces%rigidity_ice_u) .and. associated(forces%rigidity_ice_v)) &
    call uvchksum(mesg//" forces%rigidity_ice_[uv]", forces%rigidity_ice_u, &
        forces%rigidity_ice_v, G%HI, haloshift=hshift, symmetric=.true., &
        scale=US%L_to_m**3*US%L_to_Z*US%s_to_T, scalar_pair=.true.)

end subroutine MOM_mech_forcing_chksum

!> Write out values of the mechanical forcing arrays at the i,j location. This is a debugging tool.
subroutine mech_forcing_SinglePointPrint(forces, G, i, j, mesg)
  type(mech_forcing),    intent(in) :: forces !< A structure with the driving mechanical forces
  type(ocean_grid_type), intent(in) :: G      !< Grid type
  character(len=*),      intent(in) :: mesg   !< Message
  integer,               intent(in) :: i      !< i-index
  integer,               intent(in) :: j      !< j-index

  write(0,'(2a)') 'MOM_forcing_type, forcing_SinglePointPrint: Called from ',mesg
  write(0,'(a,2es15.3)') 'MOM_forcing_type, forcing_SinglePointPrint: lon,lat = ',G%geoLonT(i,j),G%geoLatT(i,j)
  call locMsg(forces%taux,'taux')
  call locMsg(forces%tauy,'tauy')

  contains
  !> Format and write a message depending on associated state of array
  subroutine locMsg(array,aname)
    real, dimension(:,:), pointer :: array !< Array to write element from
    character(len=*)              :: aname !< Name of array

    if (associated(array)) then
      write(0,'(3a,es15.3)') 'MOM_forcing_type, mech_forcing_SinglePointPrint: ',trim(aname),' = ',array(i,j)
    else
      write(0,'(4a)') 'MOM_forcing_type, mech_forcing_SinglePointPrint: ',trim(aname),' is not associated.'
    endif
  end subroutine locMsg

end subroutine mech_forcing_SinglePointPrint

!> Write out values of the fluxes arrays at the i,j location. This is a debugging tool.
subroutine forcing_SinglePointPrint(fluxes, G, i, j, mesg)
  type(forcing),         intent(in) :: fluxes !< A structure containing thermodynamic forcing fields
  type(ocean_grid_type), intent(in) :: G      !< Grid type
  character(len=*),      intent(in) :: mesg   !< Message
  integer,               intent(in) :: i      !< i-index
  integer,               intent(in) :: j      !< j-index

  write(0,'(2a)') 'MOM_forcing_type, forcing_SinglePointPrint: Called from ',mesg
  write(0,'(a,2es15.3)') 'MOM_forcing_type, forcing_SinglePointPrint: lon,lat = ',G%geoLonT(i,j),G%geoLatT(i,j)
  call locMsg(fluxes%ustar,'ustar')
  call locMsg(fluxes%buoy,'buoy')
  call locMsg(fluxes%sw,'sw')
  call locMsg(fluxes%sw_vis_dir,'sw_vis_dir')
  call locMsg(fluxes%sw_vis_dif,'sw_vis_dif')
  call locMsg(fluxes%sw_nir_dir,'sw_nir_dir')
  call locMsg(fluxes%sw_nir_dif,'sw_nir_dif')
  call locMsg(fluxes%lw,'lw')
  call locMsg(fluxes%latent,'latent')
  call locMsg(fluxes%latent_evap_diag,'latent_evap_diag')
  call locMsg(fluxes%latent_fprec_diag,'latent_fprec_diag')
  call locMsg(fluxes%latent_frunoff_diag,'latent_frunoff_diag')
  call locMsg(fluxes%sens,'sens')
  call locMsg(fluxes%evap,'evap')
  call locMsg(fluxes%lprec,'lprec')
  call locMsg(fluxes%fprec,'fprec')
  call locMsg(fluxes%vprec,'vprec')
  call locMsg(fluxes%seaice_melt,'seaice_melt')
  call locMsg(fluxes%seaice_melt_heat,'seaice_melt_heat')
  call locMsg(fluxes%p_surf,'p_surf')
  call locMsg(fluxes%salt_flux,'salt_flux')
  call locMsg(fluxes%TKE_tidal,'TKE_tidal')
  call locMsg(fluxes%ustar_tidal,'ustar_tidal')
  call locMsg(fluxes%lrunoff,'lrunoff')
  call locMsg(fluxes%frunoff,'frunoff')
  call locMsg(fluxes%heat_content_lrunoff,'heat_content_lrunoff')
  call locMsg(fluxes%heat_content_frunoff,'heat_content_frunoff')
  call locMsg(fluxes%heat_content_lprec,'heat_content_lprec')
  call locMsg(fluxes%heat_content_fprec,'heat_content_fprec')
  call locMsg(fluxes%heat_content_icemelt,'heat_content_icemelt')
  call locMsg(fluxes%heat_content_vprec,'heat_content_vprec')
  call locMsg(fluxes%heat_content_cond,'heat_content_cond')
  call locMsg(fluxes%heat_content_cond,'heat_content_massout')

  contains
  !> Format and write a message depending on associated state of array
  subroutine locMsg(array,aname)
    real, dimension(:,:), pointer :: array !< Array to write element from
    character(len=*)              :: aname !< Name of array

    if (associated(array)) then
      write(0,'(3a,es15.3)') 'MOM_forcing_type, forcing_SinglePointPrint: ',trim(aname),' = ',array(i,j)
    else
      write(0,'(4a)') 'MOM_forcing_type, forcing_SinglePointPrint: ',trim(aname),' is not associated.'
    endif
  end subroutine locMsg

end subroutine forcing_SinglePointPrint


!> Register members of the forcing type for diagnostics
subroutine register_forcing_type_diags(Time, diag, US, use_temperature, handles, use_berg_fluxes)
  type(time_type),     intent(in)    :: Time            !< time type
  type(diag_ctrl),     intent(inout) :: diag            !< diagnostic control type
  type(unit_scale_type), intent(in)  :: US              !< A dimensional unit scaling type
  logical,             intent(in)    :: use_temperature !< True if T/S are in use
  type(forcing_diags), intent(inout) :: handles         !< handles for diagnostics
  logical, optional,   intent(in)    :: use_berg_fluxes !< If true, allow iceberg flux diagnostics

  ! Clock for forcing diagnostics
  handles%id_clock_forcing=cpu_clock_id('(Ocean forcing diagnostics)', grain=CLOCK_ROUTINE)


  handles%id_taux = register_diag_field('ocean_model', 'taux', diag%axesCu1, Time,  &
        'Zonal surface stress from ocean interactions with atmos and ice', &
        'Pa', conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
        standard_name='surface_downward_x_stress', cmor_field_name='tauuo',         &
        cmor_units='N m-2', cmor_long_name='Surface Downward X Stress',             &
        cmor_standard_name='surface_downward_x_stress')

  handles%id_tauy = register_diag_field('ocean_model', 'tauy', diag%axesCv1, Time,  &
        'Meridional surface stress ocean interactions with atmos and ice', &
        'Pa',  conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
        standard_name='surface_downward_y_stress', cmor_field_name='tauvo',        &
        cmor_units='N m-2', cmor_long_name='Surface Downward Y Stress',            &
        cmor_standard_name='surface_downward_y_stress')

  handles%id_ustar = register_diag_field('ocean_model', 'ustar', diag%axesT1, Time, &
      'Surface friction velocity = [(gustiness + tau_magnitude)/rho0]^(1/2)', &
      'm s-1', conversion=US%Z_to_m*US%s_to_T)

  if (present(use_berg_fluxes)) then
    if (use_berg_fluxes) then
      handles%id_ustar_berg = register_diag_field('ocean_model', 'ustar_berg', diag%axesT1, Time, &
          'Friction velocity below iceberg ', 'm s-1', conversion=US%Z_to_m*US%s_to_T)

      handles%id_area_berg = register_diag_field('ocean_model', 'area_berg', diag%axesT1, Time, &
          'Area of grid cell covered by iceberg ', 'm2 m-2')

      handles%id_mass_berg = register_diag_field('ocean_model', 'mass_berg', diag%axesT1, Time, &
          'Mass of icebergs ', 'kg m-2', conversion=US%RZ_to_kg_m2)

      handles%id_ustar_ice_cover = register_diag_field('ocean_model', 'ustar_ice_cover', diag%axesT1, Time, &
          'Friction velocity below iceberg and ice shelf together', 'm s-1', conversion=US%Z_to_m*US%s_to_T)

      handles%id_frac_ice_cover = register_diag_field('ocean_model', 'frac_ice_cover', diag%axesT1, Time, &
          'Area of grid cell below iceberg and ice shelf together ', 'm2 m-2')
    endif
  endif

  handles%id_psurf = register_diag_field('ocean_model', 'p_surf', diag%axesT1, Time, &
        'Pressure at ice-ocean or atmosphere-ocean interface', &
        'Pa', conversion=US%RL2_T2_to_Pa, cmor_field_name='pso', &
        cmor_long_name='Sea Water Pressure at Sea Water Surface', &
        cmor_standard_name='sea_water_pressure_at_sea_water_surface')

  handles%id_TKE_tidal = register_diag_field('ocean_model', 'TKE_tidal', diag%axesT1, Time, &
        'Tidal source of BBL mixing', 'W m-2', conversion=US%RZ3_T3_to_W_m2)

  if (.not. use_temperature) then
    handles%id_buoy = register_diag_field('ocean_model', 'buoy', diag%axesT1, Time, &
          'Buoyancy forcing', 'm2 s-3', conversion=US%L_to_m**2*US%s_to_T**3)
    return
  endif


  !===============================================================
  ! surface mass flux maps

  handles%id_prcme = register_diag_field('ocean_model', 'PRCmE', diag%axesT1, Time,                  &
        'Net surface water flux (precip+melt+lrunoff+ice calving-evap)', 'kg m-2 s-1', &
        standard_name='water_flux_into_sea_water', cmor_field_name='wfo',                            &
        cmor_standard_name='water_flux_into_sea_water',cmor_long_name='Water Flux Into Sea Water')
        ! This diagnostic is rescaled to MKS units when combined.

  handles%id_evap = register_diag_field('ocean_model', 'evap', diag%axesT1, Time, &
       'Evaporation/condensation at ocean surface (evaporation is negative)', &
       'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, &
       standard_name='water_evaporation_flux', cmor_field_name='evs', &
       cmor_standard_name='water_evaporation_flux', &
       cmor_long_name='Water Evaporation Flux Where Ice Free Ocean over Sea')

  ! smg: seaice_melt field requires updates to the sea ice model
  handles%id_seaice_melt = register_diag_field('ocean_model', 'seaice_melt',       &
     diag%axesT1, Time, 'water flux to ocean from snow/sea ice melting(> 0) or formation(< 0)', &
     'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, &
      standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics',     &
      cmor_field_name='fsitherm',                                                  &
      cmor_standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics',&
      cmor_long_name='water flux to ocean from sea ice melt(> 0) or form(< 0)')

  handles%id_precip = register_diag_field('ocean_model', 'precip', diag%axesT1, Time, &
        'Liquid + frozen precipitation into ocean', 'kg m-2 s-1')
        ! This diagnostic is rescaled to MKS units when combined.

  handles%id_fprec = register_diag_field('ocean_model', 'fprec', diag%axesT1, Time,     &
        'Frozen precipitation into ocean', &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s,   &
        standard_name='snowfall_flux', cmor_field_name='prsn',                          &
        cmor_standard_name='snowfall_flux', cmor_long_name='Snowfall Flux where Ice Free Ocean over Sea')

  handles%id_lprec = register_diag_field('ocean_model', 'lprec', diag%axesT1, Time,       &
        'Liquid precipitation into ocean', &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s,                 &
        standard_name='rainfall_flux',                                                    &
        cmor_field_name='prlq', cmor_standard_name='rainfall_flux',                       &
        cmor_long_name='Rainfall Flux where Ice Free Ocean over Sea')

  handles%id_vprec = register_diag_field('ocean_model', 'vprec', diag%axesT1, Time, &
        'Virtual liquid precip into ocean due to SSS restoring', &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_frunoff = register_diag_field('ocean_model', 'frunoff', diag%axesT1, Time,    &
        'Frozen runoff (calving) and iceberg melt into ocean', &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, &
        standard_name='water_flux_into_sea_water_from_icebergs',                           &
        cmor_field_name='ficeberg',                                                        &
        cmor_standard_name='water_flux_into_sea_water_from_icebergs',                      &
        cmor_long_name='Water Flux into Seawater from Icebergs')

  handles%id_lrunoff = register_diag_field('ocean_model', 'lrunoff', diag%axesT1, Time, &
        'Liquid runoff (rivers) into ocean', &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, &
        standard_name='water_flux_into_sea_water_from_rivers', cmor_field_name='friver',      &
        cmor_standard_name='water_flux_into_sea_water_from_rivers',                           &
        cmor_long_name='Water Flux into Sea Water From Rivers')

  handles%id_net_massout = register_diag_field('ocean_model', 'net_massout', diag%axesT1, Time, &
        'Net mass leaving the ocean due to evaporation, seaice formation', 'kg m-2 s-1')
        ! This diagnostic is rescaled to MKS units when combined.

  handles%id_net_massin  = register_diag_field('ocean_model', 'net_massin', diag%axesT1, Time, &
        'Net mass entering ocean due to precip, runoff, ice melt', 'kg m-2 s-1')
        ! This diagnostic is rescaled to MKS units when combined.

  handles%id_massout_flux = register_diag_field('ocean_model', 'massout_flux', diag%axesT1, Time, &
        'Net mass flux of freshwater out of the ocean (used in the boundary flux calculation)', &
         'kg m-2', conversion=diag%GV%H_to_kg_m2)
        ! This diagnostic is calculated in MKS units.

  handles%id_massin_flux  = register_diag_field('ocean_model', 'massin_flux', diag%axesT1, Time, &
        'Net mass flux of freshwater into the ocean (used in boundary flux calculation)', 'kg m-2')
        ! This diagnostic is calculated in MKS units.

  !=========================================================================
  ! area integrated surface mass transport, all are rescaled to MKS units before area integration.

  handles%id_total_prcme = register_scalar_field('ocean_model', 'total_PRCmE', Time, diag,         &
      long_name='Area integrated net surface water flux (precip+melt+liq runoff+ice calving-evap)',&
      units='kg s-1', standard_name='water_flux_into_sea_water_area_integrated',                   &
      cmor_field_name='total_wfo',                                                                 &
      cmor_standard_name='water_flux_into_sea_water_area_integrated',                              &
      cmor_long_name='Water Transport Into Sea Water Area Integrated')

  handles%id_total_evap = register_scalar_field('ocean_model', 'total_evap', Time, diag,&
      long_name='Area integrated evap/condense at ocean surface',                       &
      units='kg s-1', standard_name='water_evaporation_flux_area_integrated',           &
      cmor_field_name='total_evs',                                                      &
      cmor_standard_name='water_evaporation_flux_area_integrated',                      &
      cmor_long_name='Evaporation Where Ice Free Ocean over Sea Area Integrated')

  ! seaice_melt field requires updates to the sea ice model
  handles%id_total_seaice_melt = register_scalar_field('ocean_model', 'total_icemelt', Time, diag, &
      long_name='Area integrated sea ice melt (>0) or form (<0)', units='kg s-1',                      &
      standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics_area_integrated',         &
      cmor_field_name='total_fsitherm',                                                                &
      cmor_standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics_area_integrated',    &
      cmor_long_name='Water Melt/Form from Sea Ice Area Integrated')

  handles%id_total_precip = register_scalar_field('ocean_model', 'total_precip', Time, diag, &
      long_name='Area integrated liquid+frozen precip into ocean', units='kg s-1')

  handles%id_total_fprec = register_scalar_field('ocean_model', 'total_fprec', Time, diag,&
      long_name='Area integrated frozen precip into ocean', units='kg s-1',               &
      standard_name='snowfall_flux_area_integrated',                                      &
      cmor_field_name='total_prsn',                                                       &
      cmor_standard_name='snowfall_flux_area_integrated',                                 &
      cmor_long_name='Snowfall Flux where Ice Free Ocean over Sea Area Integrated')

  handles%id_total_lprec = register_scalar_field('ocean_model', 'total_lprec', Time, diag,&
      long_name='Area integrated liquid precip into ocean', units='kg s-1',               &
      standard_name='rainfall_flux_area_integrated',                                      &
      cmor_field_name='total_pr',                                                         &
      cmor_standard_name='rainfall_flux_area_integrated',                                 &
      cmor_long_name='Rainfall Flux where Ice Free Ocean over Sea Area Integrated')

  handles%id_total_vprec = register_scalar_field('ocean_model', 'total_vprec', Time, diag, &
      long_name='Area integrated virtual liquid precip due to SSS restoring', units='kg s-1')

  handles%id_total_frunoff = register_scalar_field('ocean_model', 'total_frunoff', Time, diag,    &
      long_name='Area integrated frozen runoff (calving) & iceberg melt into ocean', units='kg s-1',&
      cmor_field_name='total_ficeberg',                                                           &
      cmor_standard_name='water_flux_into_sea_water_from_icebergs_area_integrated',               &
      cmor_long_name='Water Flux into Seawater from Icebergs Area Integrated')

  handles%id_total_lrunoff = register_scalar_field('ocean_model', 'total_lrunoff', Time, diag,&
      long_name='Area integrated liquid runoff into ocean', units='kg s-1',                   &
      cmor_field_name='total_friver',                                                         &
      cmor_standard_name='water_flux_into_sea_water_from_rivers_area_integrated',             &
      cmor_long_name='Water Flux into Sea Water From Rivers Area Integrated')

  handles%id_total_net_massout = register_scalar_field('ocean_model', 'total_net_massout', Time, diag, &
      long_name='Area integrated mass leaving ocean due to evap and seaice form', units='kg s-1')

  handles%id_total_net_massin = register_scalar_field('ocean_model', 'total_net_massin', Time, diag, &
      long_name='Area integrated mass entering ocean due to predip, runoff, ice melt', units='kg s-1')

  !=========================================================================
  ! area averaged surface mass transport

  handles%id_prcme_ga = register_scalar_field('ocean_model', 'PRCmE_ga', Time, diag,             &
      long_name='Area averaged net surface water flux (precip+melt+liq runoff+ice calving-evap)',&
      units='kg m-2 s-1', standard_name='water_flux_into_sea_water_area_averaged',               &
      cmor_field_name='ave_wfo',                                                                 &
      cmor_standard_name='rainfall_flux_area_averaged',                                          &
      cmor_long_name='Water Transport Into Sea Water Area Averaged')

  handles%id_evap_ga = register_scalar_field('ocean_model', 'evap_ga', Time, diag,&
      long_name='Area averaged evap/condense at ocean surface',                   &
      units='kg m-2 s-1', standard_name='water_evaporation_flux_area_averaged',   &
      cmor_field_name='ave_evs',                                                  &
      cmor_standard_name='water_evaporation_flux_area_averaged',                  &
      cmor_long_name='Evaporation Where Ice Free Ocean over Sea Area Averaged')

 handles%id_lprec_ga = register_scalar_field('ocean_model', 'lprec_ga', Time, diag,&
      long_name='Area integrated liquid precip into ocean', units='kg m-2 s-1',    &
      standard_name='rainfall_flux_area_averaged',                                 &
      cmor_field_name='ave_pr',                                                    &
      cmor_standard_name='rainfall_flux_area_averaged',                            &
      cmor_long_name='Rainfall Flux where Ice Free Ocean over Sea Area Averaged')

 handles%id_fprec_ga = register_scalar_field('ocean_model', 'fprec_ga', Time, diag,&
      long_name='Area integrated frozen precip into ocean', units='kg m-2 s-1',    &
      standard_name='snowfall_flux_area_averaged',                                 &
      cmor_field_name='ave_prsn',                                                  &
      cmor_standard_name='snowfall_flux_area_averaged',                            &
      cmor_long_name='Snowfall Flux where Ice Free Ocean over Sea Area Averaged')

  handles%id_precip_ga = register_scalar_field('ocean_model', 'precip_ga', Time, diag, &
      long_name='Area averaged liquid+frozen precip into ocean', units='kg m-2 s-1')

  handles%id_vprec_ga = register_scalar_field('ocean_model', 'vrec_ga', Time, diag, &
      long_name='Area averaged virtual liquid precip due to SSS restoring', units='kg m-2 s-1')

  !===============================================================
  ! surface heat flux maps

  handles%id_heat_content_frunoff = register_diag_field('ocean_model', 'heat_content_frunoff', &
        diag%axesT1, Time, 'Heat content (relative to 0C) of solid runoff into ocean',         &
        'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='temperature_flux_due_to_solid_runoff_expressed_as_heat_flux_into_sea_water')

  handles%id_heat_content_lrunoff = register_diag_field('ocean_model', 'heat_content_lrunoff', &
        diag%axesT1, Time, 'Heat content (relative to 0C) of liquid runoff into ocean',        &
        'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='temperature_flux_due_to_runoff_expressed_as_heat_flux_into_sea_water')

  handles%id_hfrunoffds = register_diag_field('ocean_model', 'hfrunoffds',                            &
        diag%axesT1, Time, 'Heat content (relative to 0C) of liquid+solid runoff into ocean', &
        'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='temperature_flux_due_to_runoff_expressed_as_heat_flux_into_sea_water')

  handles%id_heat_content_lprec = register_diag_field('ocean_model', 'heat_content_lprec',             &
        diag%axesT1,Time,'Heat content (relative to 0degC) of liquid precip entering ocean',           &
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_heat_content_fprec = register_diag_field('ocean_model', 'heat_content_fprec',&
        diag%axesT1,Time,'Heat content (relative to 0degC) of frozen prec entering ocean',&
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_heat_content_icemelt = register_diag_field('ocean_model', 'heat_content_icemelt',&
        diag%axesT1,Time,'Heat content (relative to 0degC) of water flux due to sea ice melting/freezing',&
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_heat_content_vprec = register_diag_field('ocean_model', 'heat_content_vprec',   &
        diag%axesT1,Time,'Heat content (relative to 0degC) of virtual precip entering ocean',&
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_heat_content_cond = register_diag_field('ocean_model', 'heat_content_cond',   &
        diag%axesT1,Time,'Heat content (relative to 0degC) of water condensing into ocean',&
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_hfrainds = register_diag_field('ocean_model', 'hfrainds',                                 &
        diag%axesT1,Time,'Heat content (relative to 0degC) of liquid+frozen precip entering ocean',    &
        'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='temperature_flux_due_to_rainfall_expressed_as_heat_flux_into_sea_water',&
        cmor_long_name='Heat Content (relative to 0degC) of Liquid + Frozen Precipitation')

  handles%id_heat_content_surfwater = register_diag_field('ocean_model', 'heat_content_surfwater',&
         diag%axesT1, Time,                                                                       &
        'Heat content (relative to 0degC) of net water crossing ocean surface (frozen+liquid)',   &
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_heat_content_massout = register_diag_field('ocean_model', 'heat_content_massout',                      &
         diag%axesT1, Time,'Heat content (relative to 0degC) of net mass leaving ocean ocean via evap and ice form',&
        'W m-2', conversion=US%QRZ_T_to_W_m2,                                                      &
        cmor_field_name='hfevapds',                                                                                 &
        cmor_standard_name='temperature_flux_due_to_evaporation_expressed_as_heat_flux_out_of_sea_water',           &
        cmor_long_name='Heat Content (relative to 0degC) of Water Leaving Ocean via Evaporation and Ice Formation')

  handles%id_heat_content_massin = register_diag_field('ocean_model', 'heat_content_massin',   &
         diag%axesT1, Time,'Heat content (relative to 0degC) of net mass entering ocean ocean',&
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_net_heat_coupler = register_diag_field('ocean_model', 'net_heat_coupler',          &
        diag%axesT1,Time,'Surface ocean heat flux from SW+LW+latent+sensible+seaice_melt_heat (via the coupler)',&
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_net_heat_surface = register_diag_field('ocean_model', 'net_heat_surface',diag%axesT1, Time,  &
        'Surface ocean heat flux from SW+LW+lat+sens+mass transfer+frazil+restore+seaice_melt_heat or '// &
        'flux adjustments', &
        'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='surface_downward_heat_flux_in_sea_water', cmor_field_name='hfds',            &
        cmor_standard_name='surface_downward_heat_flux_in_sea_water',           &
        cmor_long_name='Surface ocean heat flux from SW+LW+latent+sensible+masstransfer+frazil+seaice_melt_heat')

  handles%id_sw = register_diag_field('ocean_model', 'SW', diag%axesT1, Time,  &
        'Shortwave radiation flux into ocean', 'W m-2', conversion=US%QRZ_T_to_W_m2,                        &
        standard_name='net_downward_shortwave_flux_at_sea_water_surface',      &
        cmor_field_name='rsntds',                                              &
        cmor_standard_name='net_downward_shortwave_flux_at_sea_water_surface', &
        cmor_long_name='Net Downward Shortwave Radiation at Sea Water Surface')
  handles%id_sw_vis = register_diag_field('ocean_model', 'sw_vis', diag%axesT1, Time,     &
        'Shortwave radiation direct and diffuse flux into the ocean in the visible band', &
        'W m-2', conversion=US%QRZ_T_to_W_m2)
  handles%id_sw_nir = register_diag_field('ocean_model', 'sw_nir', diag%axesT1, Time,     &
        'Shortwave radiation direct and diffuse flux into the ocean in the near-infrared band', &
        'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_LwLatSens = register_diag_field('ocean_model', 'LwLatSens', diag%axesT1, Time, &
        'Combined longwave, latent, and sensible heating at ocean surface', 'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_lw = register_diag_field('ocean_model', 'LW', diag%axesT1, Time, &
        'Longwave radiation flux into ocean', 'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='surface_net_downward_longwave_flux',                   &
        cmor_field_name='rlntds',                                             &
        cmor_standard_name='surface_net_downward_longwave_flux',              &
        cmor_long_name='Surface Net Downward Longwave Radiation')

  handles%id_lat = register_diag_field('ocean_model', 'latent', diag%axesT1, Time,                    &
        'Latent heat flux into ocean due to fusion and evaporation (negative means ocean heat loss)', &
        'W m-2', conversion=US%QRZ_T_to_W_m2, cmor_field_name='hflso',                                &
        cmor_standard_name='surface_downward_latent_heat_flux',                                       &
        cmor_long_name='Surface Downward Latent Heat Flux due to Evap + Melt Snow/Ice')

  handles%id_lat_evap = register_diag_field('ocean_model', 'latent_evap', diag%axesT1, Time, &
        'Latent heat flux into ocean due to evaporation/condensation', 'W m-2', conversion=US%QRZ_T_to_W_m2)

  handles%id_lat_fprec = register_diag_field('ocean_model', 'latent_fprec_diag', diag%axesT1, Time,&
        'Latent heat flux into ocean due to melting of frozen precipitation', 'W m-2', conversion=US%QRZ_T_to_W_m2, &
        cmor_field_name='hfsnthermds',                                                             &
        cmor_standard_name='heat_flux_into_sea_water_due_to_snow_thermodynamics',                  &
        cmor_long_name='Latent Heat to Melt Frozen Precipitation')

  handles%id_lat_frunoff = register_diag_field('ocean_model', 'latent_frunoff', diag%axesT1, Time, &
        'Latent heat flux into ocean due to melting of icebergs', 'W m-2', conversion=US%QRZ_T_to_W_m2, &
        cmor_field_name='hfibthermds',                                                             &
        cmor_standard_name='heat_flux_into_sea_water_due_to_iceberg_thermodynamics',               &
        cmor_long_name='Latent Heat to Melt Frozen Runoff/Iceberg')

  handles%id_sens = register_diag_field('ocean_model', 'sensible', diag%axesT1, Time, &
        'Sensible heat flux into ocean', 'W m-2', conversion=US%QRZ_T_to_W_m2,        &
        standard_name='surface_downward_sensible_heat_flux',                         &
        cmor_field_name='hfsso',                                                     &
        cmor_standard_name='surface_downward_sensible_heat_flux',                    &
        cmor_long_name='Surface Downward Sensible Heat Flux')

  handles%id_seaice_melt_heat = register_diag_field('ocean_model', 'seaice_melt_heat', diag%axesT1, Time,&
        'Heat flux into ocean due to snow and sea ice melt/freeze', 'W m-2', conversion=US%QRZ_T_to_W_m2, &
        standard_name='snow_ice_melt_heat_flux',                         &
  !GMM TODO cmor_field_name='hfsso',                                                     &
        cmor_standard_name='snow_ice_melt_heat_flux',                    &
        cmor_long_name='Heat flux into ocean from snow and sea ice melt')

  handles%id_heat_added = register_diag_field('ocean_model', 'heat_added', diag%axesT1, Time, &
        'Flux Adjustment or restoring surface heat flux into ocean', 'W m-2', conversion=US%QRZ_T_to_W_m2)


  !===============================================================
  ! area integrated surface heat fluxes

  handles%id_total_heat_content_frunoff = register_scalar_field('ocean_model',                     &
      'total_heat_content_frunoff', Time, diag,                                                    &
      long_name='Area integrated heat content (relative to 0C) of solid runoff',                   &
      units='W', cmor_field_name='total_hfsolidrunoffds',                                          &
      cmor_standard_name=                                                                          &
      'temperature_flux_due_to_solid_runoff_expressed_as_heat_flux_into_sea_water_area_integrated',&
      cmor_long_name=                                                                              &
      'Temperature Flux due to Solid Runoff Expressed as Heat Flux into Sea Water Area Integrated')

  handles%id_total_heat_content_lrunoff = register_scalar_field('ocean_model',               &
      'total_heat_content_lrunoff', Time, diag,                                              &
      long_name='Area integrated heat content (relative to 0C) of liquid runoff',            &
      units='W', cmor_field_name='total_hfrunoffds',                                         &
      cmor_standard_name=                                                                    &
      'temperature_flux_due_to_runoff_expressed_as_heat_flux_into_sea_water_area_integrated',&
      cmor_long_name=                                                                        &
      'Temperature Flux due to Runoff Expressed as Heat Flux into Sea Water Area Integrated')

  handles%id_total_heat_content_lprec = register_scalar_field('ocean_model',                   &
      'total_heat_content_lprec', Time, diag,                                                  &
      long_name='Area integrated heat content (relative to 0C) of liquid precip',              &
      units='W', cmor_field_name='total_hfrainds',                                             &
      cmor_standard_name=                                                                      &
      'temperature_flux_due_to_rainfall_expressed_as_heat_flux_into_sea_water_area_integrated',&
      cmor_long_name=                                                                          &
      'Temperature Flux due to Rainfall Expressed as Heat Flux into Sea Water Area Integrated')

  handles%id_total_heat_content_fprec = register_scalar_field('ocean_model',     &
      'total_heat_content_fprec', Time, diag,                                    &
      long_name='Area integrated heat content (relative to 0C) of frozen precip',&
      units='W')

  handles%id_total_heat_content_icemelt = register_scalar_field('ocean_model',     &
      'total_heat_content_icemelt', Time, diag,long_name=                          &
      'Area integrated heat content (relative to 0C) of water flux due sea ice melting/freezing', &
      units='W')

  handles%id_total_heat_content_vprec = register_scalar_field('ocean_model',      &
      'total_heat_content_vprec', Time, diag,                                     &
      long_name='Area integrated heat content (relative to 0C) of virtual precip',&
      units='W')

  handles%id_total_heat_content_cond = register_scalar_field('ocean_model',   &
      'total_heat_content_cond', Time, diag,                                  &
      long_name='Area integrated heat content (relative to 0C) of condensate',&
      units='W')

  handles%id_total_heat_content_surfwater = register_scalar_field('ocean_model',          &
      'total_heat_content_surfwater', Time, diag,                                         &
      long_name='Area integrated heat content (relative to 0C) of water crossing surface',&
      units='W')

  handles%id_total_heat_content_massout = register_scalar_field('ocean_model',                      &
      'total_heat_content_massout', Time, diag,                                                     &
      long_name='Area integrated heat content (relative to 0C) of water leaving ocean',             &
      units='W',                                                                                    &
      cmor_field_name='total_hfevapds',                                                             &
      cmor_standard_name=                                                                           &
      'temperature_flux_due_to_evaporation_expressed_as_heat_flux_out_of_sea_water_area_integrated',&
      cmor_long_name='Heat Flux Out of Sea Water due to Evaporating Water Area Integrated')

  handles%id_total_heat_content_massin = register_scalar_field('ocean_model',           &
      'total_heat_content_massin', Time, diag,                                          &
      long_name='Area integrated heat content (relative to 0C) of water entering ocean',&
      units='W')

  handles%id_total_net_heat_coupler = register_scalar_field('ocean_model',                       &
      'total_net_heat_coupler', Time, diag,                                                      &
      long_name='Area integrated surface heat flux from SW+LW+latent+sensible+seaice_melt_heat (via the coupler)',&
      units='W')

  handles%id_total_net_heat_surface = register_scalar_field('ocean_model',                      &
      'total_net_heat_surface', Time, diag,                                                     &
      long_name='Area integrated surface heat flux from SW+LW+lat+sens+mass+frazil+restore or flux adjustments', &
      units='W',                                                                                &
      cmor_field_name='total_hfds',                                                             &
      cmor_standard_name='surface_downward_heat_flux_in_sea_water_area_integrated',             &
      cmor_long_name=                                                                           &
      'Surface Ocean Heat Flux from SW+LW+latent+sensible+mass transfer+frazil Area Integrated')

  handles%id_total_sw = register_scalar_field('ocean_model',                                &
      'total_sw', Time, diag,                                                               &
      long_name='Area integrated net downward shortwave at sea water surface',              &
      units='W',                                                                            &
      cmor_field_name='total_rsntds',                                                       &
      cmor_standard_name='net_downward_shortwave_flux_at_sea_water_surface_area_integrated',&
      cmor_long_name=                                                                       &
      'Net Downward Shortwave Radiation at Sea Water Surface Area Integrated')

  handles%id_total_LwLatSens = register_scalar_field('ocean_model',&
      'total_LwLatSens', Time, diag,                               &
      long_name='Area integrated longwave+latent+sensible heating',&
      units='W')

  handles%id_total_lw = register_scalar_field('ocean_model',                  &
      'total_lw', Time, diag,                                                 &
      long_name='Area integrated net downward longwave at sea water surface', &
      units='W',                                                              &
      cmor_field_name='total_rlntds',                                         &
      cmor_standard_name='surface_net_downward_longwave_flux_area_integrated',&
      cmor_long_name=                                                         &
      'Surface Net Downward Longwave Radiation Area Integrated')

  handles%id_total_lat = register_scalar_field('ocean_model',                &
      'total_lat', Time, diag,                                               &
      long_name='Area integrated surface downward latent heat flux',         &
      units='W',                                                             &
      cmor_field_name='total_hflso',                                         &
      cmor_standard_name='surface_downward_latent_heat_flux_area_integrated',&
      cmor_long_name=                                                        &
      'Surface Downward Latent Heat Flux Area Integrated')

  handles%id_total_lat_evap = register_scalar_field('ocean_model',      &
      'total_lat_evap', Time, diag,                                     &
      long_name='Area integrated latent heat flux due to evap/condense',&
      units='W')

  handles%id_total_lat_fprec = register_scalar_field('ocean_model',                            &
      'total_lat_fprec', Time, diag,                                                           &
      long_name='Area integrated latent heat flux due to melting frozen precip',               &
      units='W',                                                                               &
      cmor_field_name='total_hfsnthermds',                                                     &
      cmor_standard_name='heat_flux_into_sea_water_due_to_snow_thermodynamics_area_integrated',&
      cmor_long_name=                                                                          &
      'Latent Heat to Melt Frozen Precipitation Area Integrated')

  handles%id_total_lat_frunoff = register_scalar_field('ocean_model',                             &
      'total_lat_frunoff', Time, diag,                                                            &
      long_name='Area integrated latent heat flux due to melting icebergs',                       &
      units='W',                                                                                  &
      cmor_field_name='total_hfibthermds',                                                        &
      cmor_standard_name='heat_flux_into_sea_water_due_to_iceberg_thermodynamics_area_integrated',&
      cmor_long_name=                                                                             &
      'Heat Flux into Sea Water due to Iceberg Thermodynamics Area Integrated')

  handles%id_total_sens = register_scalar_field('ocean_model',                 &
      'total_sens', Time, diag,                                                &
      long_name='Area integrated downward sensible heat flux',                 &
      units='W',                                                               &
      cmor_field_name='total_hfsso',                                           &
      cmor_standard_name='surface_downward_sensible_heat_flux_area_integrated',&
      cmor_long_name=                                                          &
      'Surface Downward Sensible Heat Flux Area Integrated')

  handles%id_total_heat_added = register_scalar_field('ocean_model',&
      'total_heat_adjustment', Time, diag,                               &
      long_name='Area integrated surface heat flux from restoring and/or flux adjustment',   &
      units='W')

  handles%id_total_seaice_melt_heat = register_scalar_field('ocean_model',&
      'total_seaice_melt_heat', Time, diag,                               &
      long_name='Area integrated surface heat flux from snow and sea ice melt',   &
      units='W')

  !===============================================================
  ! area averaged surface heat fluxes

  handles%id_net_heat_coupler_ga = register_scalar_field('ocean_model',                       &
      'net_heat_coupler_ga', Time, diag,                                                      &
      long_name='Area averaged surface heat flux from SW+LW+latent+sensible+seaice_melt_heat (via the coupler)',&
      units='W m-2')

  handles%id_net_heat_surface_ga = register_scalar_field('ocean_model',                       &
      'net_heat_surface_ga', Time, diag, long_name=                                                     &
      'Area averaged surface heat flux from SW+LW+lat+sens+mass+frazil+restore+seaice_melt_heat or flux adjustments', &
      units='W m-2',                                                                          &
      cmor_field_name='ave_hfds',                                                             &
      cmor_standard_name='surface_downward_heat_flux_in_sea_water_area_averaged',             &
      cmor_long_name=                                                                         &
      'Surface Ocean Heat Flux from SW+LW+latent+sensible+mass transfer+frazil Area Averaged')

  handles%id_sw_ga = register_scalar_field('ocean_model',                                 &
      'sw_ga', Time, diag,                                                                &
      long_name='Area averaged net downward shortwave at sea water surface',              &
      units='W m-2',                                                                      &
      cmor_field_name='ave_rsntds',                                                       &
      cmor_standard_name='net_downward_shortwave_flux_at_sea_water_surface_area_averaged',&
      cmor_long_name=                                                                     &
      'Net Downward Shortwave Radiation at Sea Water Surface Area Averaged')

  handles%id_LwLatSens_ga = register_scalar_field('ocean_model',&
      'LwLatSens_ga', Time, diag,                               &
      long_name='Area averaged longwave+latent+sensible heating',&
      units='W m-2')

  handles%id_lw_ga = register_scalar_field('ocean_model',                   &
      'lw_ga', Time, diag,                                                  &
      long_name='Area averaged net downward longwave at sea water surface', &
      units='W m-2',                                                        &
      cmor_field_name='ave_rlntds',                                         &
      cmor_standard_name='surface_net_downward_longwave_flux_area_averaged',&
      cmor_long_name=                                                       &
      'Surface Net Downward Longwave Radiation Area Averaged')

  handles%id_lat_ga = register_scalar_field('ocean_model',                 &
      'lat_ga', Time, diag,                                                &
      long_name='Area averaged surface downward latent heat flux',         &
      units='W m-2',                                                       &
      cmor_field_name='ave_hflso',                                         &
      cmor_standard_name='surface_downward_latent_heat_flux_area_averaged',&
      cmor_long_name=                                                      &
      'Surface Downward Latent Heat Flux Area Averaged')

  handles%id_sens_ga = register_scalar_field('ocean_model',                  &
      'sens_ga', Time, diag,                                                 &
      long_name='Area averaged downward sensible heat flux',                 &
      units='W m-2',                                                         &
      cmor_field_name='ave_hfsso',                                           &
      cmor_standard_name='surface_downward_sensible_heat_flux_area_averaged',&
      cmor_long_name=                                                        &
      'Surface Downward Sensible Heat Flux Area Averaged')


  !===============================================================
  ! maps of surface salt fluxes, virtual precip fluxes, and adjustments

  handles%id_saltflux = register_diag_field('ocean_model', 'salt_flux', diag%axesT1, Time,&
        'Net salt flux into ocean at surface (restoring + sea-ice)',                      &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, &
        cmor_field_name='sfdsi', cmor_standard_name='downward_sea_ice_basal_salt_flux', &
        cmor_long_name='Downward Sea Ice Basal Salt Flux')

  handles%id_saltFluxIn = register_diag_field('ocean_model', 'salt_flux_in', diag%axesT1, Time, &
        'Salt flux into ocean at surface from coupler', &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_saltFluxAdded = register_diag_field('ocean_model', 'salt_flux_added', &
        diag%axesT1,Time,'Salt flux into ocean at surface due to restoring or flux adjustment', &
        units='kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)

  handles%id_saltFluxGlobalAdj = register_scalar_field('ocean_model',              &
        'salt_flux_global_restoring_adjustment', Time, diag,                       &
        'Adjustment needed to balance net global salt flux into ocean at surface', &
         units='kg m-2 s-1') !, conversion=US%RZ_T_to_kg_m2s)

  handles%id_vPrecGlobalAdj = register_scalar_field('ocean_model',  &
        'vprec_global_adjustment', Time, diag,                      &
        'Adjustment needed to adjust net vprec into ocean to zero', &
        'kg m-2 s-1')

  handles%id_netFWGlobalAdj = register_scalar_field('ocean_model',       &
        'net_fresh_water_global_adjustment', Time, diag,                 &
        'Adjustment needed to adjust net fresh water into ocean to zero',&
        'kg m-2 s-1')

  handles%id_saltFluxGlobalScl = register_scalar_field('ocean_model',            &
        'salt_flux_global_restoring_scaling', Time, diag,                        &
        'Scaling applied to balance net global salt flux into ocean at surface', &
        'nondim')

  handles%id_vPrecGlobalScl = register_scalar_field('ocean_model',&
        'vprec_global_scaling', Time, diag,                       &
        'Scaling applied to adjust net vprec into ocean to zero', &
        'nondim')

  handles%id_netFWGlobalScl = register_scalar_field('ocean_model',      &
        'net_fresh_water_global_scaling', Time, diag,                   &
        'Scaling applied to adjust net fresh water into ocean to zero', &
        'nondim')

  !===============================================================
  ! area integrals of surface salt fluxes

  handles%id_total_saltflux = register_scalar_field('ocean_model',          &
      'total_salt_flux', Time, diag,                                        &
      long_name='Area integrated surface salt flux', units='kg s-1',        &
      cmor_field_name='total_sfdsi',                                        &
      cmor_standard_name='downward_sea_ice_basal_salt_flux_area_integrated',&
      cmor_long_name='Downward Sea Ice Basal Salt Flux Area Integrated')

  handles%id_total_saltFluxIn = register_scalar_field('ocean_model', 'total_salt_Flux_In', &
      Time, diag, long_name='Area integrated surface salt flux at surface from coupler', units='kg s-1')

  handles%id_total_saltFluxAdded = register_scalar_field('ocean_model', 'total_salt_Flux_Added', &
      Time, diag, long_name='Area integrated surface salt flux due to restoring or flux adjustment', units='kg s-1')


end subroutine register_forcing_type_diags

!> Accumulate the forcing over time steps, taking input from a mechanical forcing type
!! and a temporary forcing-flux type.
subroutine forcing_accumulate(flux_tmp, forces, fluxes, G, wt2)
  type(forcing),         intent(in)    :: flux_tmp !< A temporary structure with current
                                                 !!thermodynamic forcing fields
  type(mech_forcing),    intent(in)    :: forces !< A structure with the driving mechanical forces
  type(forcing),         intent(inout) :: fluxes !< A structure containing time-averaged
                                                 !! thermodynamic forcing fields
  type(ocean_grid_type), intent(inout) :: G      !< The ocean's grid structure
  real,                  intent(out)   :: wt2    !< The relative weight of the new fluxes

  ! This subroutine copies mechancal forcing from flux_tmp to fluxes and
  ! stores the time-weighted averages of the various buoyancy fluxes in fluxes,
  ! and increments the amount of time over which the buoyancy forcing should be
  ! applied, all via a call to fluxes accumulate.

  call fluxes_accumulate(flux_tmp, fluxes, G, wt2, forces)

end subroutine forcing_accumulate

!> Accumulate the thermodynamic fluxes over time steps
subroutine fluxes_accumulate(flux_tmp, fluxes, G, wt2, forces)
  type(forcing),             intent(in)    :: flux_tmp !< A temporary structure with current
                                                     !! thermodynamic forcing fields
  type(forcing),             intent(inout) :: fluxes !< A structure containing time-averaged
                                                     !! thermodynamic forcing fields
  type(ocean_grid_type),     intent(inout) :: G      !< The ocean's grid structure
  real,                      intent(out)   :: wt2    !< The relative weight of the new fluxes
  type(mech_forcing), optional, intent(in) :: forces !< A structure with the driving mechanical forces

  ! This subroutine copies mechancal forcing from flux_tmp to fluxes and
  ! stores the time-weighted averages of the various buoyancy fluxes in fluxes,
  ! and increments the amount of time over which the buoyancy forcing in fluxes should be
  ! applied based on the time interval stored in flux_tmp.

  real :: wt1
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, i0, j0
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isr, ier, jsr, jer
  is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
  Isq  = G%IscB  ; Ieq  = G%IecB   ; Jsq  = G%JscB  ; Jeq  = G%JecB
  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB


  if (fluxes%dt_buoy_accum < 0) call MOM_error(FATAL, "fluxes_accumulate: "//&
     "fluxes must be initialzed before it can be augmented.")

  ! wt1 is the relative weight of the previous fluxes.
  wt1 = fluxes%dt_buoy_accum / (fluxes%dt_buoy_accum + flux_tmp%dt_buoy_accum)
  wt2 = 1.0 - wt1 ! = flux_tmp%dt_buoy_accum / (fluxes%dt_buoy_accum + flux_tmp%dt_buoy_accum)
  fluxes%dt_buoy_accum = fluxes%dt_buoy_accum + flux_tmp%dt_buoy_accum

  ! Copy over the pressure fields and accumulate averages of ustar, either from the forcing
  ! type or from the temporary fluxes type.
  if (present(forces)) then
    do j=js,je ; do i=is,ie
      fluxes%p_surf(i,j) = forces%p_surf(i,j)
      fluxes%p_surf_full(i,j) = forces%p_surf_full(i,j)

      fluxes%ustar(i,j) = wt1*fluxes%ustar(i,j) + wt2*forces%ustar(i,j)
    enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      fluxes%p_surf(i,j) = flux_tmp%p_surf(i,j)
      fluxes%p_surf_full(i,j) = flux_tmp%p_surf_full(i,j)

      fluxes%ustar(i,j) = wt1*fluxes%ustar(i,j) + wt2*flux_tmp%ustar(i,j)
    enddo ; enddo
  endif

  ! Average the water, heat, and salt fluxes, and ustar.
  do j=js,je ; do i=is,ie
    if (fluxes%gustless_accum_bug) then
      fluxes%ustar_gustless(i,j) = flux_tmp%ustar_gustless(i,j)
    else
      fluxes%ustar_gustless(i,j) = wt1*fluxes%ustar_gustless(i,j) + wt2*flux_tmp%ustar_gustless(i,j)
    endif

    fluxes%evap(i,j) = wt1*fluxes%evap(i,j) + wt2*flux_tmp%evap(i,j)
    fluxes%lprec(i,j) = wt1*fluxes%lprec(i,j) + wt2*flux_tmp%lprec(i,j)
    fluxes%fprec(i,j) = wt1*fluxes%fprec(i,j) + wt2*flux_tmp%fprec(i,j)
    fluxes%vprec(i,j) = wt1*fluxes%vprec(i,j) + wt2*flux_tmp%vprec(i,j)
    fluxes%lrunoff(i,j) = wt1*fluxes%lrunoff(i,j) + wt2*flux_tmp%lrunoff(i,j)
    fluxes%frunoff(i,j) = wt1*fluxes%frunoff(i,j) + wt2*flux_tmp%frunoff(i,j)
    fluxes%seaice_melt(i,j) = wt1*fluxes%seaice_melt(i,j) + wt2*flux_tmp%seaice_melt(i,j)
    fluxes%sw(i,j) = wt1*fluxes%sw(i,j) + wt2*flux_tmp%sw(i,j)
    fluxes%sw_vis_dir(i,j) = wt1*fluxes%sw_vis_dir(i,j) + wt2*flux_tmp%sw_vis_dir(i,j)
    fluxes%sw_vis_dif(i,j) = wt1*fluxes%sw_vis_dif(i,j) + wt2*flux_tmp%sw_vis_dif(i,j)
    fluxes%sw_nir_dir(i,j) = wt1*fluxes%sw_nir_dir(i,j) + wt2*flux_tmp%sw_nir_dir(i,j)
    fluxes%sw_nir_dif(i,j) = wt1*fluxes%sw_nir_dif(i,j) + wt2*flux_tmp%sw_nir_dif(i,j)
    fluxes%lw(i,j) = wt1*fluxes%lw(i,j) + wt2*flux_tmp%lw(i,j)
    fluxes%latent(i,j) = wt1*fluxes%latent(i,j) + wt2*flux_tmp%latent(i,j)
    fluxes%sens(i,j) = wt1*fluxes%sens(i,j) + wt2*flux_tmp%sens(i,j)

    fluxes%salt_flux(i,j) = wt1*fluxes%salt_flux(i,j) + wt2*flux_tmp%salt_flux(i,j)
  enddo ; enddo
  if (associated(fluxes%heat_added) .and. associated(flux_tmp%heat_added)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_added(i,j) = wt1*fluxes%heat_added(i,j) + wt2*flux_tmp%heat_added(i,j)
    enddo ; enddo
  endif
  ! These might always be associated, in which case they can be combined?
  if (associated(fluxes%heat_content_cond) .and. associated(flux_tmp%heat_content_cond)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_cond(i,j) = wt1*fluxes%heat_content_cond(i,j) + wt2*flux_tmp%heat_content_cond(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_lprec) .and. associated(flux_tmp%heat_content_lprec)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_lprec(i,j) = wt1*fluxes%heat_content_lprec(i,j) + wt2*flux_tmp%heat_content_lprec(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_fprec) .and. associated(flux_tmp%heat_content_fprec)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_fprec(i,j) = wt1*fluxes%heat_content_fprec(i,j) + wt2*flux_tmp%heat_content_fprec(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_icemelt) .and. associated(flux_tmp%heat_content_icemelt)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_icemelt(i,j) = wt1*fluxes%heat_content_icemelt(i,j) + wt2*flux_tmp%heat_content_icemelt(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_vprec) .and. associated(flux_tmp%heat_content_vprec)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_vprec(i,j) = wt1*fluxes%heat_content_vprec(i,j) + wt2*flux_tmp%heat_content_vprec(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_lrunoff) .and. associated(flux_tmp%heat_content_lrunoff)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_lrunoff(i,j) = wt1*fluxes%heat_content_lrunoff(i,j) + wt2*flux_tmp%heat_content_lrunoff(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_frunoff) .and. associated(flux_tmp%heat_content_frunoff)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_frunoff(i,j) = wt1*fluxes%heat_content_frunoff(i,j) + wt2*flux_tmp%heat_content_frunoff(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%heat_content_icemelt) .and. associated(flux_tmp%heat_content_icemelt)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_content_icemelt(i,j) = wt1*fluxes%heat_content_icemelt(i,j) + wt2*flux_tmp%heat_content_icemelt(i,j)
    enddo ; enddo
  endif

  if (associated(fluxes%ustar_shelf) .and. associated(flux_tmp%ustar_shelf)) then
    do i=isd,ied ; do j=jsd,jed
      fluxes%ustar_shelf(i,j)  = flux_tmp%ustar_shelf(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%iceshelf_melt) .and. associated(flux_tmp%iceshelf_melt)) then
    do i=isd,ied ; do j=jsd,jed
      fluxes%iceshelf_melt(i,j)  = flux_tmp%iceshelf_melt(i,j)
    enddo ; enddo
  endif
  if (associated(fluxes%frac_shelf_h) .and. associated(flux_tmp%frac_shelf_h)) then
    do i=isd,ied ; do j=jsd,jed
      fluxes%frac_shelf_h(i,j)  = flux_tmp%frac_shelf_h(i,j)
    enddo ; enddo
  endif

  if (coupler_type_initialized(fluxes%tr_fluxes) .and. &
      coupler_type_initialized(flux_tmp%tr_fluxes)) &
    call coupler_type_increment_data(flux_tmp%tr_fluxes, fluxes%tr_fluxes, &
                              scale_factor=wt2, scale_prev=wt1)

end subroutine fluxes_accumulate

!> This subroutine copies the computational domains of common forcing fields
!! from a mech_forcing type to a (thermodynamic) forcing type.
subroutine copy_common_forcing_fields(forces, fluxes, G, skip_pres)
  type(mech_forcing),      intent(in)    :: forces   !< A structure with the driving mechanical forces
  type(forcing),           intent(inout) :: fluxes   !< A structure containing thermodynamic forcing fields
  type(ocean_grid_type),   intent(in)    :: G        !< grid type
  logical,       optional, intent(in)    :: skip_pres !< If present and true, do not copy pressure fields.

  logical :: do_pres
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do_pres = .true. ; if (present(skip_pres)) do_pres = .not.skip_pres

  if (associated(forces%ustar) .and. associated(fluxes%ustar)) then
    do j=js,je ; do i=is,ie
      fluxes%ustar(i,j) = forces%ustar(i,j)
    enddo ; enddo
  endif

  if (do_pres) then
    if (associated(forces%p_surf) .and. associated(fluxes%p_surf)) then
      do j=js,je ; do i=is,ie
        fluxes%p_surf(i,j) = forces%p_surf(i,j)
      enddo ; enddo
    endif

    if (associated(forces%p_surf_full) .and. associated(fluxes%p_surf_full)) then
      do j=js,je ; do i=is,ie
        fluxes%p_surf_full(i,j) = forces%p_surf_full(i,j)
      enddo ; enddo
    endif

    if (associated(forces%p_surf_SSH, forces%p_surf_full)) then
      fluxes%p_surf_SSH => fluxes%p_surf_full
    elseif (associated(forces%p_surf_SSH, forces%p_surf)) then
      fluxes%p_surf_SSH => fluxes%p_surf
    endif
  endif

end subroutine copy_common_forcing_fields

!> This subroutine calculates certain derived forcing fields based on information
!! from a mech_forcing type and stores them in a (thermodynamic) forcing type.
subroutine set_derived_forcing_fields(forces, fluxes, G, US, Rho0)
  type(mech_forcing),      intent(in)    :: forces   !< A structure with the driving mechanical forces
  type(forcing),           intent(inout) :: fluxes   !< A structure containing thermodynamic forcing fields
  type(ocean_grid_type),   intent(in)    :: G        !< grid type
  type(unit_scale_type),   intent(in)    :: US       !< A dimensional unit scaling type
  real,                    intent(in)    :: Rho0     !< A reference density of seawater [R ~> kg m-3],
                                                     !! as used to calculate ustar.

  real :: taux2, tauy2 ! Squared wind stress components [R2 L2 Z2 T-4 ~> Pa2].
  real :: Irho0        ! Inverse of the mean density rescaled to [Z L-1 R-1 ~> m3 kg-1]
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  Irho0 = US%L_to_Z / Rho0

  if (associated(forces%taux) .and. associated(forces%tauy) .and. &
      associated(fluxes%ustar_gustless)) then
    do j=js,je ; do i=is,ie
      taux2 = 0.0
      if ((G%mask2dCu(I-1,j) + G%mask2dCu(I,j)) > 0) &
        taux2 = (G%mask2dCu(I-1,j) * forces%taux(I-1,j)**2 + &
                 G%mask2dCu(I,j) * forces%taux(I,j)**2) / &
                (G%mask2dCu(I-1,j) + G%mask2dCu(I,j))
      tauy2 = 0.0
      if ((G%mask2dCv(i,J-1) + G%mask2dCv(i,J)) > 0) &
        tauy2 = (G%mask2dCv(i,J-1) * forces%tauy(i,J-1)**2 + &
                 G%mask2dCv(i,J) * forces%tauy(i,J)**2) / &
                (G%mask2dCv(i,J-1) + G%mask2dCv(i,J))

      if (fluxes%gustless_accum_bug) then
        ! This change is just for computational efficiency, but it is wrapped with another change.
        fluxes%ustar_gustless(i,j) = sqrt(US%L_to_Z * sqrt(taux2 + tauy2) / Rho0)
      else
        fluxes%ustar_gustless(i,j) = sqrt(sqrt(taux2 + tauy2) * Irho0)
      endif
    enddo ; enddo
  endif

end subroutine set_derived_forcing_fields


!> This subroutine determines the net mass source to the ocean from
!! a (thermodynamic) forcing type and stores it in a mech_forcing type.
subroutine set_net_mass_forcing(fluxes, forces, G, US)
  type(forcing),           intent(in)    :: fluxes   !< A structure containing thermodynamic forcing fields
  type(mech_forcing),      intent(inout) :: forces   !< A structure with the driving mechanical forces
  type(unit_scale_type),   intent(in)    :: US       !< A dimensional unit scaling type
  type(ocean_grid_type),   intent(in)    :: G        !< The ocean grid type

  if (associated(forces%net_mass_src)) &
    call get_net_mass_forcing(fluxes, G, US, forces%net_mass_src)

end subroutine set_net_mass_forcing

!> This subroutine calculates determines the net mass source to the ocean from
!! a (thermodynamic) forcing type and stores it in a provided array.
subroutine get_net_mass_forcing(fluxes, G, US, net_mass_src)
  type(forcing),                    intent(in)  :: fluxes !< A structure containing thermodynamic forcing fields
  type(ocean_grid_type),            intent(in)  :: G      !< The ocean grid type
  type(unit_scale_type),            intent(in)  :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: net_mass_src !< The net mass flux of water into the ocean
                                                          !! [kg m-2 s-1].

  real :: RZ_T_conversion ! A combination of scaling factors for mass fluxes [kg T m-2 s-1 R-1 Z-1 ~> 1]
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  RZ_T_conversion = US%RZ_T_to_kg_m2s

  net_mass_src(:,:) = 0.0
  if (associated(fluxes%lprec)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + RZ_T_conversion*fluxes%lprec(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%fprec)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + RZ_T_conversion*fluxes%fprec(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%vprec)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + RZ_T_conversion*fluxes%vprec(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%lrunoff)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + RZ_T_conversion*fluxes%lrunoff(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%frunoff)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + RZ_T_conversion*fluxes%frunoff(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%evap)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + RZ_T_conversion*fluxes%evap(i,j)
  enddo ; enddo ; endif
  if (associated(fluxes%seaice_melt)) then ; do j=js,je ; do i=is,ie
    net_mass_src(i,j) = net_mass_src(i,j) + RZ_T_conversion*fluxes%seaice_melt(i,j)
  enddo ; enddo ; endif

end subroutine get_net_mass_forcing

!> This subroutine copies the computational domains of common forcing fields
!! from a mech_forcing type to a (thermodynamic) forcing type.
subroutine copy_back_forcing_fields(fluxes, forces, G)
  type(forcing),           intent(in)    :: fluxes   !< A structure containing thermodynamic forcing fields
  type(mech_forcing),      intent(inout) :: forces   !< A structure with the driving mechanical forces
  type(ocean_grid_type),   intent(in)    :: G        !< grid type

  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (associated(forces%ustar) .and. associated(fluxes%ustar)) then
    do j=js,je ; do i=is,ie
      forces%ustar(i,j) = fluxes%ustar(i,j)
    enddo ; enddo
  endif

end subroutine copy_back_forcing_fields

!> Offer mechanical forcing fields for diagnostics for those
!! fields registered as part of register_forcing_type_diags.
subroutine mech_forcing_diags(forces_in, dt, G, time_end, diag, handles)
  type(mech_forcing), target, intent(in) :: forces_in !< mechanical forcing input fields
  real,                  intent(in)    :: dt       !< time step for the forcing [s]
  type(ocean_grid_type), intent(in)    :: G        !< grid type
  type(time_type),       intent(in)    :: time_end !< The end time of the diagnostic interval.
  type(diag_ctrl),       intent(inout) :: diag     !< diagnostic type
  type(forcing_diags),   intent(inout) :: handles  !< diagnostic id for diag_manager

  integer :: i,j,is,ie,js,je

  type(mech_forcing), pointer :: forces
  integer :: turns

  call cpu_clock_begin(handles%id_clock_forcing)

  ! NOTE: post_data expects data to be on the rotated index map, so any
  !   rotations must be applied before saving the output.
  turns = diag%G%HI%turns
  if (turns /= 0) then
    allocate(forces)
    call allocate_mech_forcing(forces_in, diag%G, forces)
    call rotate_mech_forcing(forces_in, turns, forces)
  else
    forces => forces_in
  endif

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  call enable_averaging(dt, time_end, diag)
  ! if (query_averaging_enabled(diag)) then

    if ((handles%id_taux > 0) .and. associated(forces%taux)) &
      call post_data(handles%id_taux, forces%taux, diag)

    if ((handles%id_tauy > 0) .and. associated(forces%tauy)) &
      call post_data(handles%id_tauy, forces%tauy, diag)

    if ((handles%id_mass_berg > 0) .and. associated(forces%mass_berg)) &
      call post_data(handles%id_mass_berg, forces%mass_berg, diag)

    if ((handles%id_area_berg > 0) .and. associated(forces%area_berg)) &
      call post_data(handles%id_area_berg, forces%area_berg, diag)

  ! endif

  call disable_averaging(diag)

  if (turns /= 0) then
    call deallocate_mech_forcing(forces)
    deallocate(forces)
  endif

  call cpu_clock_end(handles%id_clock_forcing)
end subroutine mech_forcing_diags


!> Offer buoyancy forcing fields for diagnostics for those
!! fields registered as part of register_forcing_type_diags.
subroutine forcing_diagnostics(fluxes_in, sfc_state, G_in, US, time_end, diag, handles)
  type(forcing), target, intent(in)    :: fluxes_in !< A structure containing thermodynamic forcing fields
  type(surface),         intent(in)    :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(ocean_grid_type), target, intent(in) :: G_in !< Input grid type
  type(unit_scale_type), intent(in)    :: US        !< A dimensional unit scaling type
  type(time_type),       intent(in)    :: time_end  !< The end time of the diagnostic interval.
  type(diag_ctrl),       intent(inout) :: diag      !< diagnostic regulator
  type(forcing_diags),   intent(inout) :: handles   !< diagnostic ids

  ! local
  type(ocean_grid_type), pointer :: G   ! Grid metric on model index map
  type(forcing), pointer :: fluxes      ! Fluxes on the model index map
  real, dimension(SZI_(diag%G),SZJ_(diag%G)) :: res
  real :: total_transport ! for diagnosing integrated boundary transport
  real :: ave_flux        ! for diagnosing averaged   boundary flux
  real :: RZ_T_conversion ! A combination of scaling factors for mass fluxes [kg T m-2 s-1 R-1 Z-1 ~> 1]
  real :: I_dt            ! inverse time step [T-1 ~> s-1]
  real :: ppt2mks         ! conversion between ppt and mks
  integer :: turns        ! Number of index quarter turns
  integer :: i,j,is,ie,js,je

  call cpu_clock_begin(handles%id_clock_forcing)

  ! NOTE: post_data expects data to be on the rotated index map, so any
  !   rotations must be applied before saving the output.
  turns = diag%G%HI%turns
  if (turns /= 0) then
    G => diag%G
    allocate(fluxes)
    call allocate_forcing_type(fluxes_in, G, fluxes)
    call rotate_forcing(fluxes_in, fluxes, turns)
  else
    G => G_in
    fluxes => fluxes_in
  endif

  RZ_T_conversion = US%RZ_T_to_kg_m2s
  I_dt    = 1.0 / fluxes%dt_buoy_accum
  ppt2mks = 1e-3
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  call enable_averages(fluxes%dt_buoy_accum, time_end, diag)
  ! if (query_averaging_enabled(diag)) then

    ! post the diagnostics for surface mass fluxes ==================================

    if (handles%id_prcme > 0 .or. handles%id_total_prcme > 0 .or. handles%id_prcme_ga > 0) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0
        if (associated(fluxes%lprec))       res(i,j) = res(i,j) + RZ_T_conversion*fluxes%lprec(i,j)
        if (associated(fluxes%fprec))       res(i,j) = res(i,j) + RZ_T_conversion*fluxes%fprec(i,j)
        ! fluxes%cond is not needed because it is derived from %evap > 0
        if (associated(fluxes%evap))        res(i,j) = res(i,j) + RZ_T_conversion*fluxes%evap(i,j)
        if (associated(fluxes%lrunoff))     res(i,j) = res(i,j) + RZ_T_conversion*fluxes%lrunoff(i,j)
        if (associated(fluxes%frunoff))     res(i,j) = res(i,j) + RZ_T_conversion*fluxes%frunoff(i,j)
        if (associated(fluxes%vprec))       res(i,j) = res(i,j) + RZ_T_conversion*fluxes%vprec(i,j)
        if (associated(fluxes%seaice_melt)) res(i,j) = res(i,j) + RZ_T_conversion*fluxes%seaice_melt(i,j)
      enddo ; enddo
      if (handles%id_prcme > 0) call post_data(handles%id_prcme, res, diag)
      if (handles%id_total_prcme > 0) then
        total_transport = global_area_integral(res, G)
        call post_data(handles%id_total_prcme, total_transport, diag)
      endif
      if (handles%id_prcme_ga > 0) then
        ave_flux = global_area_mean(res, G)
        call post_data(handles%id_prcme_ga, ave_flux, diag)
      endif
    endif

    if (handles%id_net_massout > 0 .or. handles%id_total_net_massout > 0) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0
        if (associated(fluxes%lprec)) then
          if (fluxes%lprec(i,j) < 0.0) res(i,j) = res(i,j) + RZ_T_conversion*fluxes%lprec(i,j)
        endif
        if (associated(fluxes%vprec)) then
          if (fluxes%vprec(i,j) < 0.0) res(i,j) = res(i,j) + RZ_T_conversion*fluxes%vprec(i,j)
        endif
        if (associated(fluxes%evap)) then
          if (fluxes%evap(i,j) < 0.0) res(i,j) = res(i,j) + RZ_T_conversion*fluxes%evap(i,j)
        endif
        if (associated(fluxes%seaice_melt)) then
          if (fluxes%seaice_melt(i,j) < 0.0) &
            res(i,j) = res(i,j) + RZ_T_conversion*fluxes%seaice_melt(i,j)
        endif
      enddo ; enddo
      if (handles%id_net_massout > 0) call post_data(handles%id_net_massout, res, diag)
      if (handles%id_total_net_massout > 0) then
        total_transport = global_area_integral(res, G)
        call post_data(handles%id_total_net_massout, total_transport, diag)
      endif
    endif

    if (handles%id_massout_flux > 0 .and. associated(fluxes%netMassOut)) &
      call post_data(handles%id_massout_flux,fluxes%netMassOut,diag)

    if (handles%id_net_massin > 0 .or. handles%id_total_net_massin > 0) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0

        if (associated(fluxes%fprec)) &
          res(i,j) = res(i,j) + RZ_T_conversion*fluxes%fprec(i,j)
        if (associated(fluxes%lrunoff)) &
          res(i,j) = res(i,j) + RZ_T_conversion*fluxes%lrunoff(i,j)
        if (associated(fluxes%frunoff)) &
          res(i,j) = res(i,j) + RZ_T_conversion*fluxes%frunoff(i,j)

        if (associated(fluxes%lprec)) then
          if (fluxes%lprec(i,j) > 0.0) res(i,j) = res(i,j) + RZ_T_conversion*fluxes%lprec(i,j)
        endif
        if (associated(fluxes%vprec)) then
          if (fluxes%vprec(i,j) > 0.0) res(i,j) = res(i,j) + RZ_T_conversion*fluxes%vprec(i,j)
        endif
        ! fluxes%cond is not needed because it is derived from %evap > 0
        if (associated(fluxes%evap)) then
          if (fluxes%evap(i,j) > 0.0) res(i,j) = res(i,j) + RZ_T_conversion*fluxes%evap(i,j)
        endif
        if (associated(fluxes%seaice_melt)) then
          if (fluxes%seaice_melt(i,j) > 0.0) &
            res(i,j) = res(i,j) + RZ_T_conversion*fluxes%seaice_melt(i,j)
        endif
      enddo ; enddo
      if (handles%id_net_massin > 0) call post_data(handles%id_net_massin, res, diag)
      if (handles%id_total_net_massin > 0) then
        total_transport = global_area_integral(res, G)
        call post_data(handles%id_total_net_massin, total_transport, diag)
      endif
    endif

    if (handles%id_massin_flux > 0 .and. associated(fluxes%netMassIn)) &
      call post_data(handles%id_massin_flux,fluxes%netMassIn,diag)

    if ((handles%id_evap > 0) .and. associated(fluxes%evap)) &
      call post_data(handles%id_evap, fluxes%evap, diag)
    if ((handles%id_total_evap > 0) .and. associated(fluxes%evap)) then
      total_transport = global_area_integral(fluxes%evap, G, scale=US%RZ_T_to_kg_m2s)
      call post_data(handles%id_total_evap, total_transport, diag)
    endif
    if ((handles%id_evap_ga > 0) .and. associated(fluxes%evap)) then
      ave_flux = global_area_mean(fluxes%evap, G, scale=US%RZ_T_to_kg_m2s)
      call post_data(handles%id_evap_ga, ave_flux, diag)
    endif

    if (associated(fluxes%lprec) .and. associated(fluxes%fprec)) then
      do j=js,je ; do i=is,ie
        res(i,j) = RZ_T_conversion* (fluxes%lprec(i,j) + fluxes%fprec(i,j))
      enddo ; enddo
      if (handles%id_precip > 0) call post_data(handles%id_precip, res, diag)
      if (handles%id_total_precip > 0) then
        total_transport = global_area_integral(res, G)
        call post_data(handles%id_total_precip, total_transport, diag)
      endif
      if (handles%id_precip_ga > 0) then
        ave_flux = global_area_mean(res, G)
        call post_data(handles%id_precip_ga, ave_flux, diag)
      endif
    endif

    if (associated(fluxes%lprec)) then
      if (handles%id_lprec > 0) call post_data(handles%id_lprec, fluxes%lprec, diag)
      if (handles%id_total_lprec > 0) then
        total_transport = global_area_integral(fluxes%lprec, G, scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_lprec, total_transport, diag)
      endif
      if (handles%id_lprec_ga > 0) then
        ave_flux = global_area_mean(fluxes%lprec, G, scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_lprec_ga, ave_flux, diag)
      endif
    endif

    if (associated(fluxes%fprec)) then
      if (handles%id_fprec > 0) call post_data(handles%id_fprec, fluxes%fprec, diag)
      if (handles%id_total_fprec > 0) then
        total_transport = global_area_integral(fluxes%fprec, G, scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_fprec, total_transport, diag)
      endif
      if (handles%id_fprec_ga > 0) then
        ave_flux = global_area_mean(fluxes%fprec, G, scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_fprec_ga, ave_flux, diag)
      endif
    endif

    if (associated(fluxes%vprec)) then
      if (handles%id_vprec > 0) call post_data(handles%id_vprec, fluxes%vprec, diag)
      if (handles%id_total_vprec > 0) then
        total_transport = global_area_integral(fluxes%vprec, G, scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_vprec, total_transport, diag)
      endif
      if (handles%id_vprec_ga > 0) then
        ave_flux = global_area_mean(fluxes%vprec, G, scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_vprec_ga, ave_flux, diag)
      endif
    endif

    if (associated(fluxes%lrunoff)) then
    if (handles%id_lrunoff > 0) call post_data(handles%id_lrunoff, fluxes%lrunoff, diag)
      if (handles%id_total_lrunoff > 0) then
        total_transport = global_area_integral(fluxes%lrunoff, G, scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_lrunoff, total_transport, diag)
      endif
    endif

    if (associated(fluxes%frunoff)) then
      if (handles%id_frunoff > 0) call post_data(handles%id_frunoff, fluxes%frunoff, diag)
      if (handles%id_total_frunoff > 0) then
        total_transport = global_area_integral(fluxes%frunoff, G, scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_frunoff, total_transport, diag)
      endif
    endif

    if (associated(fluxes%seaice_melt)) then
      if (handles%id_seaice_melt > 0) call post_data(handles%id_seaice_melt, fluxes%seaice_melt, diag)
      if (handles%id_total_seaice_melt > 0) then
        total_transport = global_area_integral(fluxes%seaice_melt, G, scale=US%RZ_T_to_kg_m2s)
        call post_data(handles%id_total_seaice_melt, total_transport, diag)
      endif
    endif

    ! post diagnostics for boundary heat fluxes ====================================

    if ((handles%id_heat_content_lrunoff > 0) .and. associated(fluxes%heat_content_lrunoff))  &
      call post_data(handles%id_heat_content_lrunoff, fluxes%heat_content_lrunoff, diag)
    if ((handles%id_total_heat_content_lrunoff > 0) .and. associated(fluxes%heat_content_lrunoff)) then
      total_transport = global_area_integral(fluxes%heat_content_lrunoff, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_lrunoff, total_transport, diag)
    endif

    if ((handles%id_heat_content_frunoff > 0) .and. associated(fluxes%heat_content_frunoff))  &
      call post_data(handles%id_heat_content_frunoff, fluxes%heat_content_frunoff, diag)
    if ((handles%id_total_heat_content_frunoff > 0) .and. associated(fluxes%heat_content_frunoff)) then
      total_transport = global_area_integral(fluxes%heat_content_frunoff, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_frunoff, total_transport, diag)
    endif

    if ((handles%id_heat_content_lprec > 0) .and. associated(fluxes%heat_content_lprec))      &
      call post_data(handles%id_heat_content_lprec, fluxes%heat_content_lprec, diag)
    if ((handles%id_total_heat_content_lprec > 0) .and. associated(fluxes%heat_content_lprec)) then
      total_transport = global_area_integral(fluxes%heat_content_lprec, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_lprec, total_transport, diag)
    endif

    if ((handles%id_heat_content_fprec > 0) .and. associated(fluxes%heat_content_fprec))      &
      call post_data(handles%id_heat_content_fprec, fluxes%heat_content_fprec, diag)
    if ((handles%id_total_heat_content_fprec > 0) .and. associated(fluxes%heat_content_fprec)) then
      total_transport = global_area_integral(fluxes%heat_content_fprec, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_fprec, total_transport, diag)
    endif

    if ((handles%id_heat_content_icemelt > 0) .and. associated(fluxes%heat_content_icemelt))      &
      call post_data(handles%id_heat_content_icemelt, fluxes%heat_content_icemelt, diag)
    if ((handles%id_total_heat_content_icemelt > 0) .and. associated(fluxes%heat_content_icemelt)) then
      total_transport = global_area_integral(fluxes%heat_content_icemelt, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_icemelt, total_transport, diag)
    endif

    if ((handles%id_heat_content_vprec > 0) .and. associated(fluxes%heat_content_vprec))      &
      call post_data(handles%id_heat_content_vprec, fluxes%heat_content_vprec, diag)
    if ((handles%id_total_heat_content_vprec > 0) .and. associated(fluxes%heat_content_vprec)) then
      total_transport = global_area_integral(fluxes%heat_content_vprec, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_vprec, total_transport, diag)
    endif

    if ((handles%id_heat_content_cond > 0) .and. associated(fluxes%heat_content_cond))        &
      call post_data(handles%id_heat_content_cond, fluxes%heat_content_cond, diag)
    if ((handles%id_total_heat_content_cond > 0) .and. associated(fluxes%heat_content_cond)) then
      total_transport = global_area_integral(fluxes%heat_content_cond, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_cond, total_transport, diag)
    endif

    if ((handles%id_heat_content_massout > 0) .and. associated(fluxes%heat_content_massout))  &
      call post_data(handles%id_heat_content_massout, fluxes%heat_content_massout, diag)
    if ((handles%id_total_heat_content_massout > 0) .and. associated(fluxes%heat_content_massout)) then
      total_transport = global_area_integral(fluxes%heat_content_massout, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_massout, total_transport, diag)
    endif

    if ((handles%id_heat_content_massin > 0) .and. associated(fluxes%heat_content_massin))  &
      call post_data(handles%id_heat_content_massin, fluxes%heat_content_massin, diag)
    if ((handles%id_total_heat_content_massin > 0) .and. associated(fluxes%heat_content_massin)) then
      total_transport = global_area_integral(fluxes%heat_content_massin, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_content_massin, total_transport, diag)
    endif

    if (handles%id_net_heat_coupler > 0 .or. handles%id_total_net_heat_coupler > 0 .or. &
        handles%id_net_heat_coupler_ga > 0. ) then
      do j=js,je ; do i=is,ie
      res(i,j) = 0.0
      if (associated(fluxes%LW))               res(i,j) = res(i,j) + fluxes%lw(i,j)
      if (associated(fluxes%latent))           res(i,j) = res(i,j) + fluxes%latent(i,j)
      if (associated(fluxes%sens))             res(i,j) = res(i,j) + fluxes%sens(i,j)
      if (associated(fluxes%SW))               res(i,j) = res(i,j) + fluxes%sw(i,j)
      if (associated(fluxes%seaice_melt_heat)) res(i,j) = res(i,j) + fluxes%seaice_melt_heat(i,j)
      enddo ; enddo
      if (handles%id_net_heat_coupler > 0) call post_data(handles%id_net_heat_coupler, res, diag)
      if (handles%id_total_net_heat_coupler > 0) then
        total_transport = global_area_integral(res, G, scale=US%QRZ_T_to_W_m2)
        call post_data(handles%id_total_net_heat_coupler, total_transport, diag)
      endif
      if (handles%id_net_heat_coupler_ga > 0) then
        ave_flux = global_area_mean(res, G, scale=US%QRZ_T_to_W_m2)
        call post_data(handles%id_net_heat_coupler_ga, ave_flux, diag)
      endif
    endif

    if (handles%id_net_heat_surface > 0 .or. handles%id_total_net_heat_surface > 0 .or. &
        handles%id_net_heat_surface_ga > 0. ) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0
        if (associated(fluxes%LW))               res(i,j) = res(i,j) + fluxes%lw(i,j)
        if (associated(fluxes%latent))           res(i,j) = res(i,j) + fluxes%latent(i,j)
        if (associated(fluxes%sens))             res(i,j) = res(i,j) + fluxes%sens(i,j)
        if (associated(fluxes%SW))               res(i,j) = res(i,j) + fluxes%sw(i,j)
        if (associated(fluxes%seaice_melt_heat)) res(i,j) = res(i,j) + fluxes%seaice_melt_heat(i,j)
        if (allocated(sfc_state%frazil))         res(i,j) = res(i,j) + sfc_state%frazil(i,j) * I_dt
        !if (associated(sfc_state%TempXpme)) then
        !  res(i,j) = res(i,j) + sfc_state%TempXpme(i,j) * fluxes%C_p * I_dt
        !else
          if (associated(fluxes%heat_content_lrunoff)) &
            res(i,j) = res(i,j) + fluxes%heat_content_lrunoff(i,j)
          if (associated(fluxes%heat_content_frunoff)) &
            res(i,j) = res(i,j) + fluxes%heat_content_frunoff(i,j)
          if (associated(fluxes%heat_content_lprec)) &
            res(i,j) = res(i,j) + fluxes%heat_content_lprec(i,j)
          if (associated(fluxes%heat_content_fprec)) &
            res(i,j) = res(i,j) + fluxes%heat_content_fprec(i,j)
          if (associated(fluxes%heat_content_icemelt)) &
            res(i,j) = res(i,j) + fluxes%heat_content_icemelt(i,j)
          if (associated(fluxes%heat_content_vprec)) &
            res(i,j) = res(i,j) + fluxes%heat_content_vprec(i,j)
          if (associated(fluxes%heat_content_cond)) &
            res(i,j) = res(i,j) + fluxes%heat_content_cond(i,j)
          if (associated(fluxes%heat_content_massout)) &
            res(i,j) = res(i,j) + fluxes%heat_content_massout(i,j)
        !endif
        if (associated(fluxes%heat_added)) res(i,j) = res(i,j) + fluxes%heat_added(i,j)
      enddo ; enddo
      if (handles%id_net_heat_surface > 0) call post_data(handles%id_net_heat_surface, res, diag)

      if (handles%id_total_net_heat_surface > 0) then
        total_transport = global_area_integral(res, G, scale=US%QRZ_T_to_W_m2)
        call post_data(handles%id_total_net_heat_surface, total_transport, diag)
      endif
      if (handles%id_net_heat_surface_ga > 0) then
        ave_flux = global_area_mean(res, G, scale=US%QRZ_T_to_W_m2)
        call post_data(handles%id_net_heat_surface_ga, ave_flux, diag)
      endif
    endif

    if (handles%id_heat_content_surfwater > 0 .or. handles%id_total_heat_content_surfwater > 0) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0
      ! if (associated(sfc_state%TempXpme)) then
      !   res(i,j) = res(i,j) + sfc_state%TempXpme(i,j) * fluxes%C_p * I_dt
      ! else
          if (associated(fluxes%heat_content_lrunoff)) res(i,j) = res(i,j) + fluxes%heat_content_lrunoff(i,j)
          if (associated(fluxes%heat_content_frunoff)) res(i,j) = res(i,j) + fluxes%heat_content_frunoff(i,j)
          if (associated(fluxes%heat_content_lprec))   res(i,j) = res(i,j) + fluxes%heat_content_lprec(i,j)
          if (associated(fluxes%heat_content_icemelt)) res(i,j) = res(i,j) + fluxes%heat_content_icemelt(i,j)
          if (associated(fluxes%heat_content_fprec))   res(i,j) = res(i,j) + fluxes%heat_content_fprec(i,j)
          if (associated(fluxes%heat_content_vprec))   res(i,j) = res(i,j) + fluxes%heat_content_vprec(i,j)
          if (associated(fluxes%heat_content_cond))    res(i,j) = res(i,j) + fluxes%heat_content_cond(i,j)
          if (associated(fluxes%heat_content_massout)) res(i,j) = res(i,j) + fluxes%heat_content_massout(i,j)
      ! endif
      enddo ; enddo
      if (handles%id_heat_content_surfwater > 0) call post_data(handles%id_heat_content_surfwater, res, diag)
      if (handles%id_total_heat_content_surfwater > 0) then
        total_transport = global_area_integral(res, G, scale=US%QRZ_T_to_W_m2)
        call post_data(handles%id_total_heat_content_surfwater, total_transport, diag)
      endif
    endif

    ! for OMIP, hfrunoffds = heat content of liquid plus frozen runoff
    if (handles%id_hfrunoffds > 0) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0
        if (associated(fluxes%heat_content_lrunoff)) res(i,j) = res(i,j) + fluxes%heat_content_lrunoff(i,j)
        if (associated(fluxes%heat_content_frunoff)) res(i,j) = res(i,j) + fluxes%heat_content_frunoff(i,j)
      enddo ; enddo
      call post_data(handles%id_hfrunoffds, res, diag)
    endif

    ! for OMIP, hfrainds = heat content of lprec + fprec + cond
    if (handles%id_hfrainds > 0) then
      do j=js,je ; do i=is,ie
        res(i,j) = 0.0
        if (associated(fluxes%heat_content_lprec)) res(i,j) = res(i,j) + fluxes%heat_content_lprec(i,j)
        if (associated(fluxes%heat_content_fprec)) res(i,j) = res(i,j) + fluxes%heat_content_fprec(i,j)
        if (associated(fluxes%heat_content_cond)) res(i,j) = res(i,j) + fluxes%heat_content_cond(i,j)
      enddo ; enddo
      call post_data(handles%id_hfrainds, res, diag)
    endif

    if ((handles%id_LwLatSens > 0) .and. associated(fluxes%lw) .and. &
         associated(fluxes%latent) .and. associated(fluxes%sens)) then
      do j=js,je ; do i=is,ie
        res(i,j) = (fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)
      enddo ; enddo
      call post_data(handles%id_LwLatSens, res, diag)
    endif

    if ((handles%id_total_LwLatSens > 0) .and. associated(fluxes%lw) .and. &
         associated(fluxes%latent) .and. associated(fluxes%sens)) then
      do j=js,je ; do i=is,ie
        res(i,j) = (fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)
      enddo ; enddo
      total_transport = global_area_integral(res, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_LwLatSens, total_transport, diag)
    endif

    if ((handles%id_LwLatSens_ga > 0) .and. associated(fluxes%lw) .and. &
         associated(fluxes%latent) .and. associated(fluxes%sens)) then
      do j=js,je ; do i=is,ie
        res(i,j) = ((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j))
      enddo ; enddo
      ave_flux = global_area_mean(res, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_LwLatSens_ga, ave_flux, diag)
    endif

    if ((handles%id_sw > 0) .and. associated(fluxes%sw)) then
      call post_data(handles%id_sw, fluxes%sw, diag)
    endif
    if ((handles%id_sw_vis > 0) .and. associated(fluxes%sw_vis_dir) .and. &
        associated(fluxes%sw_vis_dif)) then
      call post_data(handles%id_sw_vis, fluxes%sw_vis_dir+fluxes%sw_vis_dif, diag)
    endif
    if ((handles%id_sw_nir > 0) .and. associated(fluxes%sw_nir_dir) .and. &
        associated(fluxes%sw_nir_dif)) then
      call post_data(handles%id_sw_nir, fluxes%sw_nir_dir+fluxes%sw_nir_dif, diag)
    endif
    if ((handles%id_total_sw > 0) .and. associated(fluxes%sw)) then
      total_transport = global_area_integral(fluxes%sw, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_sw, total_transport, diag)
    endif
    if ((handles%id_sw_ga > 0) .and. associated(fluxes%sw)) then
      ave_flux = global_area_mean(fluxes%sw, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_sw_ga, ave_flux, diag)
    endif

    if ((handles%id_lw > 0) .and. associated(fluxes%lw)) then
      call post_data(handles%id_lw, fluxes%lw, diag)
    endif
    if ((handles%id_total_lw > 0) .and. associated(fluxes%lw)) then
      total_transport = global_area_integral(fluxes%lw, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_lw, total_transport, diag)
    endif
    if ((handles%id_lw_ga > 0) .and. associated(fluxes%lw)) then
      ave_flux = global_area_mean(fluxes%lw, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_lw_ga, ave_flux, diag)
    endif

    if ((handles%id_lat > 0) .and. associated(fluxes%latent)) then
      call post_data(handles%id_lat, fluxes%latent, diag)
    endif
    if ((handles%id_total_lat > 0) .and. associated(fluxes%latent)) then
      total_transport = global_area_integral(fluxes%latent, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_lat, total_transport, diag)
    endif
    if ((handles%id_lat_ga > 0) .and. associated(fluxes%latent)) then
      ave_flux = global_area_mean(fluxes%latent, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_lat_ga, ave_flux, diag)
    endif

    if ((handles%id_lat_evap > 0) .and. associated(fluxes%latent_evap_diag)) then
      call post_data(handles%id_lat_evap, fluxes%latent_evap_diag, diag)
    endif
    if ((handles%id_total_lat_evap > 0) .and. associated(fluxes%latent_evap_diag)) then
      total_transport = global_area_integral(fluxes%latent_evap_diag, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_lat_evap, total_transport, diag)
    endif

    if ((handles%id_lat_fprec > 0) .and. associated(fluxes%latent_fprec_diag)) then
      call post_data(handles%id_lat_fprec, fluxes%latent_fprec_diag, diag)
    endif
    if ((handles%id_total_lat_fprec > 0) .and. associated(fluxes%latent_fprec_diag)) then
      total_transport = global_area_integral(fluxes%latent_fprec_diag, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_lat_fprec, total_transport, diag)
    endif

    if ((handles%id_lat_frunoff > 0) .and. associated(fluxes%latent_frunoff_diag)) then
      call post_data(handles%id_lat_frunoff, fluxes%latent_frunoff_diag, diag)
    endif
    if (handles%id_total_lat_frunoff > 0 .and. associated(fluxes%latent_frunoff_diag)) then
      total_transport = global_area_integral(fluxes%latent_frunoff_diag, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_lat_frunoff, total_transport, diag)
    endif

    if ((handles%id_sens > 0) .and. associated(fluxes%sens)) then
      call post_data(handles%id_sens, fluxes%sens, diag)
    endif

    if ((handles%id_seaice_melt_heat > 0) .and. associated(fluxes%seaice_melt_heat)) then
      call post_data(handles%id_seaice_melt_heat, fluxes%seaice_melt_heat, diag)
    endif

    if ((handles%id_total_seaice_melt_heat > 0) .and. associated(fluxes%seaice_melt_heat)) then
      total_transport = global_area_integral(fluxes%seaice_melt_heat, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_seaice_melt_heat, total_transport, diag)
    endif

    if ((handles%id_total_sens > 0) .and. associated(fluxes%sens)) then
      total_transport = global_area_integral(fluxes%sens, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_sens, total_transport, diag)
    endif
    if ((handles%id_sens_ga > 0) .and. associated(fluxes%sens)) then
      ave_flux = global_area_mean(fluxes%sens, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_sens_ga, ave_flux, diag)
    endif

    if ((handles%id_heat_added > 0) .and. associated(fluxes%heat_added)) then
      call post_data(handles%id_heat_added, fluxes%heat_added, diag)
    endif

    if ((handles%id_total_heat_added > 0) .and. associated(fluxes%heat_added)) then
      total_transport = global_area_integral(fluxes%heat_added, G, scale=US%QRZ_T_to_W_m2)
      call post_data(handles%id_total_heat_added, total_transport, diag)
    endif


    ! post the diagnostics for boundary salt fluxes ==========================

    if ((handles%id_saltflux > 0) .and. associated(fluxes%salt_flux)) &
      call post_data(handles%id_saltflux, fluxes%salt_flux, diag)
    if ((handles%id_total_saltflux > 0) .and. associated(fluxes%salt_flux)) then
      total_transport = ppt2mks*global_area_integral(fluxes%salt_flux, G, scale=US%RZ_T_to_kg_m2s)
      call post_data(handles%id_total_saltflux, total_transport, diag)
    endif

    if ((handles%id_saltFluxAdded > 0) .and. associated(fluxes%salt_flux_added)) &
      call post_data(handles%id_saltFluxAdded, fluxes%salt_flux_added, diag)
    if ((handles%id_total_saltFluxAdded > 0) .and. associated(fluxes%salt_flux_added)) then
      total_transport = ppt2mks*global_area_integral(fluxes%salt_flux_added, G, scale=US%RZ_T_to_kg_m2s)
      call post_data(handles%id_total_saltFluxAdded, total_transport, diag)
    endif

    if (handles%id_saltFluxIn > 0 .and. associated(fluxes%salt_flux_in)) &
      call post_data(handles%id_saltFluxIn, fluxes%salt_flux_in, diag)
    if ((handles%id_total_saltFluxIn > 0) .and. associated(fluxes%salt_flux_in)) then
      total_transport = ppt2mks*global_area_integral(fluxes%salt_flux_in, G, scale=US%RZ_T_to_kg_m2s)
      call post_data(handles%id_total_saltFluxIn, total_transport, diag)
    endif

    if (handles%id_saltFluxGlobalAdj > 0)                                            &
      call post_data(handles%id_saltFluxGlobalAdj, fluxes%saltFluxGlobalAdj, diag)
    if (handles%id_vPrecGlobalAdj > 0)                                               &
      call post_data(handles%id_vPrecGlobalAdj, fluxes%vPrecGlobalAdj, diag)
    if (handles%id_netFWGlobalAdj > 0)                                               &
      call post_data(handles%id_netFWGlobalAdj, fluxes%netFWGlobalAdj, diag)
    if (handles%id_saltFluxGlobalScl > 0)                                            &
      call post_data(handles%id_saltFluxGlobalScl, fluxes%saltFluxGlobalScl, diag)
    if (handles%id_vPrecGlobalScl > 0)                                               &
      call post_data(handles%id_vPrecGlobalScl, fluxes%vPrecGlobalScl, diag)
    if (handles%id_netFWGlobalScl > 0)                                               &
      call post_data(handles%id_netFWGlobalScl, fluxes%netFWGlobalScl, diag)


    ! remaining boundary terms ==================================================

    if ((handles%id_psurf > 0) .and. associated(fluxes%p_surf))                      &
      call post_data(handles%id_psurf, fluxes%p_surf, diag)

    if ((handles%id_TKE_tidal > 0) .and. associated(fluxes%TKE_tidal))               &
      call post_data(handles%id_TKE_tidal, fluxes%TKE_tidal, diag)

    if ((handles%id_buoy > 0) .and. associated(fluxes%buoy))                         &
      call post_data(handles%id_buoy, fluxes%buoy, diag)

    if ((handles%id_ustar > 0) .and. associated(fluxes%ustar)) &
      call post_data(handles%id_ustar, fluxes%ustar, diag)

    if ((handles%id_ustar_berg > 0) .and. associated(fluxes%ustar_berg)) &
      call post_data(handles%id_ustar_berg, fluxes%ustar_berg, diag)

    if ((handles%id_frac_ice_cover > 0) .and. associated(fluxes%frac_shelf_h)) &
      call post_data(handles%id_frac_ice_cover, fluxes%frac_shelf_h, diag)

    if ((handles%id_ustar_ice_cover > 0) .and. associated(fluxes%ustar_shelf)) &
      call post_data(handles%id_ustar_ice_cover, fluxes%ustar_shelf, diag)

  ! endif  ! query_averaging_enabled
  call disable_averaging(diag)

  if (turns /= 0) then
    call deallocate_forcing_type(fluxes)
    deallocate(fluxes)
  endif

  call cpu_clock_end(handles%id_clock_forcing)
end subroutine forcing_diagnostics


!> Conditionally allocate fields within the forcing type
subroutine allocate_forcing_by_group(G, fluxes, water, heat, ustar, press, &
                                     shelf, iceberg, salt, fix_accum_bug)
  type(ocean_grid_type), intent(in) :: G       !< Ocean grid structure
  type(forcing),      intent(inout) :: fluxes  !< A structure containing thermodynamic forcing fields
  logical, optional,     intent(in) :: water   !< If present and true, allocate water fluxes
  logical, optional,     intent(in) :: heat    !< If present and true, allocate heat fluxes
  logical, optional,     intent(in) :: ustar   !< If present and true, allocate ustar and related fields
  logical, optional,     intent(in) :: press   !< If present and true, allocate p_surf and related fields
  logical, optional,     intent(in) :: shelf   !< If present and true, allocate fluxes for ice-shelf
  logical, optional,     intent(in) :: iceberg !< If present and true, allocate fluxes for icebergs
  logical, optional,     intent(in) :: salt    !< If present and true, allocate salt fluxes
  logical, optional,     intent(in) :: fix_accum_bug !< If present and true, avoid using a bug in
                                               !! accumulation of ustar_gustless

  ! Local variables
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  logical :: heat_water

  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB

  call myAlloc(fluxes%ustar,isd,ied,jsd,jed, ustar)
  call myAlloc(fluxes%ustar_gustless,isd,ied,jsd,jed, ustar)

  call myAlloc(fluxes%evap,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%lprec,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%fprec,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%vprec,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%lrunoff,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%frunoff,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%seaice_melt,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%netMassOut,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%netMassIn,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%netSalt,isd,ied,jsd,jed, water)
  call myAlloc(fluxes%seaice_melt_heat,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%sw,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%lw,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%latent,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%sens,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%latent_evap_diag,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%latent_fprec_diag,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%latent_frunoff_diag,isd,ied,jsd,jed, heat)

  call myAlloc(fluxes%salt_flux,isd,ied,jsd,jed, salt)

  if (present(heat) .and. present(water)) then ; if (heat .and. water) then
    call myAlloc(fluxes%heat_content_cond,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_icemelt,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_lprec,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_fprec,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_vprec,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_lrunoff,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_frunoff,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_massout,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_massin,isd,ied,jsd,jed, .true.)
  endif ; endif

  call myAlloc(fluxes%p_surf,isd,ied,jsd,jed, press)

  call myAlloc(fluxes%frac_shelf_h,isd,ied,jsd,jed, shelf)
  call myAlloc(fluxes%ustar_shelf,isd,ied,jsd,jed, shelf)
  call myAlloc(fluxes%iceshelf_melt,isd,ied,jsd,jed, shelf)

  !These fields should only on allocated when iceberg area is being passed through the coupler.
  call myAlloc(fluxes%ustar_berg,isd,ied,jsd,jed, iceberg)
  call myAlloc(fluxes%area_berg,isd,ied,jsd,jed, iceberg)
  call myAlloc(fluxes%mass_berg,isd,ied,jsd,jed, iceberg)

  if (present(fix_accum_bug)) fluxes%gustless_accum_bug = .not.fix_accum_bug
end subroutine allocate_forcing_by_group


subroutine allocate_forcing_by_ref(fluxes_ref, G, fluxes)
  type(forcing), intent(in) :: fluxes_ref  !< Reference fluxes
  type(ocean_grid_type), intent(in) :: G        !< Grid metric of target fluxes
  type(forcing), intent(out) :: fluxes     !< Target fluxes

  logical :: do_ustar, do_water, do_heat, do_salt, do_press, do_shelf, &
      do_iceberg, do_heat_added, do_buoy

  call get_forcing_groups(fluxes_ref, do_water, do_heat, do_ustar, do_press, &
      do_shelf, do_iceberg, do_salt, do_heat_added, do_buoy)

  call allocate_forcing_type(G, fluxes, do_water, do_heat, do_ustar, &
      do_press, do_shelf, do_iceberg, do_salt)

  ! The following fluxes would typically be allocated by the driver
  call myAlloc(fluxes%sw_vis_dir, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%sw_vis_dir))
  call myAlloc(fluxes%sw_vis_dif, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%sw_vis_dif))
  call myAlloc(fluxes%sw_nir_dir, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%sw_nir_dir))
  call myAlloc(fluxes%sw_nir_dif, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%sw_nir_dif))

  call myAlloc(fluxes%salt_flux_in, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%salt_flux_in))
  call myAlloc(fluxes%salt_flux_added, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%salt_flux_added))

  call myAlloc(fluxes%p_surf_full, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%p_surf_full))

  call myAlloc(fluxes%heat_added, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%heat_added))
  call myAlloc(fluxes%buoy, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%buoy))

  call myAlloc(fluxes%TKE_tidal, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%TKE_tidal))
  call myAlloc(fluxes%ustar_tidal, G%isd, G%ied, G%jsd, G%jed, &
      associated(fluxes_ref%ustar_tidal))

  ! This flag would normally be set by a control flag in allocate_forcing_type.
  ! Here we copy the flag from the reference forcing.
  fluxes%gustless_accum_bug = fluxes_ref%gustless_accum_bug
end subroutine allocate_forcing_by_ref


!> Conditionally allocate fields within the mechanical forcing type using
!! control flags.
subroutine allocate_mech_forcing_by_group(G, forces, stress, ustar, shelf, &
                                          press, iceberg)
  type(ocean_grid_type), intent(in) :: G       !< Ocean grid structure
  type(mech_forcing), intent(inout) :: forces  !< Forcing fields structure

  logical, optional,     intent(in) :: stress  !< If present and true, allocate taux, tauy
  logical, optional,     intent(in) :: ustar   !< If present and true, allocate ustar and related fields
  logical, optional,     intent(in) :: shelf   !< If present and true, allocate forces for ice-shelf
  logical, optional,     intent(in) :: press   !< If present and true, allocate p_surf and related fields
  logical, optional,     intent(in) :: iceberg !< If present and true, allocate forces for icebergs

  ! Local variables
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  logical :: heat_water

  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB

  call myAlloc(forces%taux,IsdB,IedB,jsd,jed, stress)
  call myAlloc(forces%tauy,isd,ied,JsdB,JedB, stress)

  call myAlloc(forces%ustar,isd,ied,jsd,jed, ustar)

  call myAlloc(forces%p_surf,isd,ied,jsd,jed, press)
  call myAlloc(forces%p_surf_full,isd,ied,jsd,jed, press)
  call myAlloc(forces%net_mass_src,isd,ied,jsd,jed, press)

  call myAlloc(forces%rigidity_ice_u,IsdB,IedB,jsd,jed, shelf)
  call myAlloc(forces%rigidity_ice_v,isd,ied,JsdB,JedB, shelf)
  call myAlloc(forces%frac_shelf_u,IsdB,IedB,jsd,jed, shelf)
  call myAlloc(forces%frac_shelf_v,isd,ied,JsdB,JedB, shelf)

  !These fields should only on allocated when iceberg area is being passed through the coupler.
  call myAlloc(forces%area_berg,isd,ied,jsd,jed, iceberg)
  call myAlloc(forces%mass_berg,isd,ied,jsd,jed, iceberg)
end subroutine allocate_mech_forcing_by_group


!> Conditionally allocate fields within the mechanical forcing type based on a
!! reference forcing.
subroutine allocate_mech_forcing_from_ref(forces_ref, G, forces)
  type(mech_forcing), intent(in) :: forces_ref  !< Reference forcing fields
  type(ocean_grid_type), intent(in) :: G      !< Grid metric of target forcing
  type(mech_forcing), intent(out) :: forces   !< Mechanical forcing fields

  logical :: do_stress, do_ustar, do_shelf, do_press, do_iceberg

  ! Identify the active fields in the reference forcing
  call get_mech_forcing_groups(forces_ref, do_stress, do_ustar, do_shelf, &
                              do_press, do_iceberg)

  call allocate_mech_forcing(G, forces, do_stress, do_ustar, do_shelf, &
                             do_press, do_iceberg)
end subroutine allocate_mech_forcing_from_ref


!> Return flags indicating which groups of forcings are allocated
subroutine get_forcing_groups(fluxes, water, heat, ustar, press, shelf, &
                             iceberg, salt, heat_added, buoy)
  type(forcing), intent(in) :: fluxes  !< Reference flux fields
  logical, intent(out) :: water   !< True if fluxes contains water-based fluxes
  logical, intent(out) :: heat    !< True if fluxes contains heat-based fluxes
  logical, intent(out) :: ustar   !< True if fluxes contains ustar fluxes
  logical, intent(out) :: press   !< True if fluxes contains surface pressure
  logical, intent(out) :: shelf   !< True if fluxes contains ice shelf fields
  logical, intent(out) :: iceberg !< True if fluxes contains iceberg fluxes
  logical, intent(out) :: salt    !< True if fluxes contains salt flux
  logical, intent(out) :: heat_added !< True if fluxes contains explicit heat
  logical, intent(out) :: buoy    !< True if fluxes contains buoyancy fluxes

  ! NOTE: heat, salt, heat_added, and buoy would typically depend on each other
  !   to some degree.  But since this would be enforced at the driver level,
  !   we handle them here as independent flags.

  ustar = associated(fluxes%ustar) &
      .and. associated(fluxes%ustar_gustless)
  ! TODO: Check for all associated fields, but for now just check one as a marker
  water = associated(fluxes%evap)
  heat = associated(fluxes%seaice_melt_heat)
  salt = associated(fluxes%salt_flux)
  press = associated(fluxes%p_surf)
  shelf = associated(fluxes%frac_shelf_h)
  iceberg = associated(fluxes%ustar_berg)
  heat_added = associated(fluxes%heat_added)
  buoy = associated(fluxes%buoy)
end subroutine get_forcing_groups


!> Return flags indicating which groups of mechanical forcings are allocated
subroutine get_mech_forcing_groups(forces, stress, ustar, shelf, press, iceberg)
  type(mech_forcing), intent(in) :: forces  !< Reference forcing fields
  logical, intent(out) :: stress  !< True if forces contains wind stress fields
  logical, intent(out) :: ustar   !< True if forces contains ustar field
  logical, intent(out) :: shelf   !< True if forces contains ice shelf fields
  logical, intent(out) :: press   !< True if forces contains pressure fields
  logical, intent(out) :: iceberg !< True if forces contains iceberg fields

  stress = associated(forces%taux) &
      .and. associated(forces%tauy)
  ustar = associated(forces%ustar)
  shelf = associated(forces%rigidity_ice_u) &
      .and. associated(forces%rigidity_ice_v) &
      .and. associated(forces%frac_shelf_u) &
      .and. associated(forces%frac_shelf_v)
  press = associated(forces%p_surf) &
      .and. associated(forces%p_surf_full) &
      .and. associated(forces%net_mass_src)
  iceberg = associated(forces%area_berg) &
      .and. associated(forces%mass_berg)
end subroutine get_mech_forcing_groups


!> Allocates and zeroes-out array.
subroutine myAlloc(array, is, ie, js, je, flag)
  real, dimension(:,:), pointer :: array !< Array to be allocated
  integer,           intent(in) :: is !< Start i-index
  integer,           intent(in) :: ie !< End i-index
  integer,           intent(in) :: js !< Start j-index
  integer,           intent(in) :: je !< End j-index
  logical, optional, intent(in) :: flag !< Flag to indicate to allocate

  if (present(flag)) then ; if (flag) then ; if (.not.associated(array)) then
    allocate(array(is:ie,js:je)) ; array(is:ie,js:je) = 0.0
  endif ; endif ; endif
end subroutine myAlloc

!> Deallocate the forcing type
subroutine deallocate_forcing_type(fluxes)
  type(forcing), intent(inout) :: fluxes !< Forcing fields structure

  if (associated(fluxes%ustar))                deallocate(fluxes%ustar)
  if (associated(fluxes%ustar_gustless))       deallocate(fluxes%ustar_gustless)
  if (associated(fluxes%buoy))                 deallocate(fluxes%buoy)
  if (associated(fluxes%sw))                   deallocate(fluxes%sw)
  if (associated(fluxes%seaice_melt_heat))     deallocate(fluxes%seaice_melt_heat)
  if (associated(fluxes%sw_vis_dir))           deallocate(fluxes%sw_vis_dir)
  if (associated(fluxes%sw_vis_dif))           deallocate(fluxes%sw_vis_dif)
  if (associated(fluxes%sw_nir_dir))           deallocate(fluxes%sw_nir_dir)
  if (associated(fluxes%sw_nir_dif))           deallocate(fluxes%sw_nir_dif)
  if (associated(fluxes%lw))                   deallocate(fluxes%lw)
  if (associated(fluxes%latent))               deallocate(fluxes%latent)
  if (associated(fluxes%latent_evap_diag))     deallocate(fluxes%latent_evap_diag)
  if (associated(fluxes%latent_fprec_diag))    deallocate(fluxes%latent_fprec_diag)
  if (associated(fluxes%latent_frunoff_diag))  deallocate(fluxes%latent_frunoff_diag)
  if (associated(fluxes%sens))                 deallocate(fluxes%sens)
  if (associated(fluxes%heat_added))           deallocate(fluxes%heat_added)
  if (associated(fluxes%heat_content_lrunoff)) deallocate(fluxes%heat_content_lrunoff)
  if (associated(fluxes%heat_content_frunoff)) deallocate(fluxes%heat_content_frunoff)
  if (associated(fluxes%heat_content_icemelt)) deallocate(fluxes%heat_content_icemelt)
  if (associated(fluxes%heat_content_lprec))   deallocate(fluxes%heat_content_lprec)
  if (associated(fluxes%heat_content_fprec))   deallocate(fluxes%heat_content_fprec)
  if (associated(fluxes%heat_content_cond))    deallocate(fluxes%heat_content_cond)
  if (associated(fluxes%heat_content_massout)) deallocate(fluxes%heat_content_massout)
  if (associated(fluxes%heat_content_massin))  deallocate(fluxes%heat_content_massin)
  if (associated(fluxes%evap))                 deallocate(fluxes%evap)
  if (associated(fluxes%lprec))                deallocate(fluxes%lprec)
  if (associated(fluxes%fprec))                deallocate(fluxes%fprec)
  if (associated(fluxes%vprec))                deallocate(fluxes%vprec)
  if (associated(fluxes%lrunoff))              deallocate(fluxes%lrunoff)
  if (associated(fluxes%frunoff))              deallocate(fluxes%frunoff)
  if (associated(fluxes%seaice_melt))          deallocate(fluxes%seaice_melt)
  if (associated(fluxes%salt_flux))            deallocate(fluxes%salt_flux)
  if (associated(fluxes%p_surf_full))          deallocate(fluxes%p_surf_full)
  if (associated(fluxes%p_surf))               deallocate(fluxes%p_surf)
  if (associated(fluxes%TKE_tidal))            deallocate(fluxes%TKE_tidal)
  if (associated(fluxes%ustar_tidal))          deallocate(fluxes%ustar_tidal)
  if (associated(fluxes%ustar_shelf))          deallocate(fluxes%ustar_shelf)
  if (associated(fluxes%iceshelf_melt))        deallocate(fluxes%iceshelf_melt)
  if (associated(fluxes%frac_shelf_h))         deallocate(fluxes%frac_shelf_h)
  if (associated(fluxes%ustar_berg))           deallocate(fluxes%ustar_berg)
  if (associated(fluxes%area_berg))            deallocate(fluxes%area_berg)
  if (associated(fluxes%mass_berg))            deallocate(fluxes%mass_berg)

  call coupler_type_destructor(fluxes%tr_fluxes)

end subroutine deallocate_forcing_type


!> Deallocate the mechanical forcing type
subroutine deallocate_mech_forcing(forces)
  type(mech_forcing), intent(inout) :: forces  !< Forcing fields structure

  if (associated(forces%taux))  deallocate(forces%taux)
  if (associated(forces%tauy))  deallocate(forces%tauy)
  if (associated(forces%ustar)) deallocate(forces%ustar)
  if (associated(forces%p_surf))         deallocate(forces%p_surf)
  if (associated(forces%p_surf_full))    deallocate(forces%p_surf_full)
  if (associated(forces%net_mass_src))   deallocate(forces%net_mass_src)
  if (associated(forces%rigidity_ice_u)) deallocate(forces%rigidity_ice_u)
  if (associated(forces%rigidity_ice_v)) deallocate(forces%rigidity_ice_v)
  if (associated(forces%frac_shelf_u))   deallocate(forces%frac_shelf_u)
  if (associated(forces%frac_shelf_v))   deallocate(forces%frac_shelf_v)
  if (associated(forces%area_berg))      deallocate(forces%area_berg)
  if (associated(forces%mass_berg))      deallocate(forces%mass_berg)

end subroutine deallocate_mech_forcing


!< Rotate the fluxes by a set number of quarter turns
subroutine rotate_forcing(fluxes_in, fluxes, turns)
  type(forcing), intent(in)  :: fluxes_in     !< Input forcing struct
  type(forcing), intent(inout) :: fluxes      !< Rotated forcing struct
  integer, intent(in) :: turns                !< Number of quarter turns

  logical :: do_ustar, do_water, do_heat, do_salt, do_press, do_shelf, &
      do_iceberg, do_heat_added, do_buoy

  call get_forcing_groups(fluxes_in, do_water, do_heat, do_ustar, do_press, &
      do_shelf, do_iceberg, do_salt, do_heat_added, do_buoy)

  if (do_ustar) then
    call rotate_array(fluxes_in%ustar, turns, fluxes%ustar)
    call rotate_array(fluxes_in%ustar_gustless, turns, fluxes%ustar_gustless)
  endif

  if (do_water) then
    call rotate_array(fluxes_in%evap, turns, fluxes%evap)
    call rotate_array(fluxes_in%lprec, turns, fluxes%lprec)
    call rotate_array(fluxes_in%fprec, turns, fluxes%fprec)
    call rotate_array(fluxes_in%vprec, turns, fluxes%vprec)
    call rotate_array(fluxes_in%lrunoff, turns, fluxes%lrunoff)
    call rotate_array(fluxes_in%frunoff, turns, fluxes%frunoff)
    call rotate_array(fluxes_in%seaice_melt, turns, fluxes%seaice_melt)
    call rotate_array(fluxes_in%netMassOut, turns, fluxes%netMassOut)
    call rotate_array(fluxes_in%netMassIn, turns, fluxes%netMassIn)
    call rotate_array(fluxes_in%netSalt, turns, fluxes%netSalt)
  endif

  if (do_heat) then
    call rotate_array(fluxes_in%seaice_melt_heat, turns, fluxes%seaice_melt_heat)
    call rotate_array(fluxes_in%sw, turns, fluxes%sw)
    call rotate_array(fluxes_in%lw, turns, fluxes%lw)
    call rotate_array(fluxes_in%latent, turns, fluxes%latent)
    call rotate_array(fluxes_in%sens, turns, fluxes%sens)
    call rotate_array(fluxes_in%latent_evap_diag, turns, fluxes%latent_evap_diag)
    call rotate_array(fluxes_in%latent_fprec_diag, turns, fluxes%latent_fprec_diag)
    call rotate_array(fluxes_in%latent_frunoff_diag, turns, fluxes%latent_frunoff_diag)
  endif

  if (do_salt) then
    call rotate_array(fluxes_in%salt_flux, turns, fluxes%salt_flux)
  endif

  if (do_heat .and. do_water) then
    call rotate_array(fluxes_in%heat_content_cond, turns, fluxes%heat_content_cond)
    call rotate_array(fluxes_in%heat_content_icemelt, turns, fluxes%heat_content_icemelt)
    call rotate_array(fluxes_in%heat_content_lprec, turns, fluxes%heat_content_lprec)
    call rotate_array(fluxes_in%heat_content_fprec, turns, fluxes%heat_content_fprec)
    call rotate_array(fluxes_in%heat_content_vprec, turns, fluxes%heat_content_vprec)
    call rotate_array(fluxes_in%heat_content_lrunoff, turns, fluxes%heat_content_lrunoff)
    call rotate_array(fluxes_in%heat_content_frunoff, turns, fluxes%heat_content_frunoff)
    call rotate_array(fluxes_in%heat_content_massout, turns, fluxes%heat_content_massout)
    call rotate_array(fluxes_in%heat_content_massin, turns, fluxes%heat_content_massin)
  endif

  if (do_press) then
    call rotate_array(fluxes_in%p_surf, turns, fluxes%p_surf)
  endif

  if (do_shelf) then
    call rotate_array(fluxes_in%frac_shelf_h, turns, fluxes%frac_shelf_h)
    call rotate_array(fluxes_in%ustar_shelf, turns, fluxes%ustar_shelf)
    call rotate_array(fluxes_in%iceshelf_melt, turns, fluxes%iceshelf_melt)
  endif

  if (do_iceberg) then
    call rotate_array(fluxes_in%ustar_berg, turns, fluxes%ustar_berg)
    call rotate_array(fluxes_in%area_berg, turns, fluxes%area_berg)
    call rotate_array(fluxes_in%iceshelf_melt, turns, fluxes%iceshelf_melt)
  endif

  if (do_heat_added) then
    call rotate_array(fluxes_in%heat_added, turns, fluxes%heat_added)
  endif

  ! The following fields are handled by drivers rather than control flags.
  if (associated(fluxes_in%sw_vis_dir)) &
    call rotate_array(fluxes_in%sw_vis_dir, turns, fluxes%sw_vis_dir)
  if (associated(fluxes_in%sw_vis_dif)) &
    call rotate_array(fluxes_in%sw_vis_dif, turns, fluxes%sw_vis_dif)
  if (associated(fluxes_in%sw_nir_dir)) &
    call rotate_array(fluxes_in%sw_nir_dir, turns, fluxes%sw_nir_dir)
  if (associated(fluxes_in%sw_nir_dif)) &
    call rotate_array(fluxes_in%sw_nir_dif, turns, fluxes%sw_nir_dif)

  if (associated(fluxes_in%salt_flux_in)) &
    call rotate_array(fluxes_in%salt_flux_in, turns, fluxes%salt_flux_in)
  if (associated(fluxes_in%salt_flux_added)) &
    call rotate_array(fluxes_in%salt_flux_added, turns, fluxes%salt_flux_added)

  if (associated(fluxes_in%p_surf_full)) &
    call rotate_array(fluxes_in%p_surf_full, turns, fluxes%p_surf_full)

  if (associated(fluxes_in%buoy)) &
    call rotate_array(fluxes_in%buoy, turns, fluxes%buoy)

  if (associated(fluxes_in%TKE_tidal)) &
    call rotate_array(fluxes_in%TKE_tidal, turns, fluxes%TKE_tidal)
  if (associated(fluxes_in%ustar_tidal)) &
    call rotate_array(fluxes_in%ustar_tidal, turns, fluxes%ustar_tidal)

  ! TODO: tracer flux rotation
  if (coupler_type_initialized(fluxes%tr_fluxes)) &
    call MOM_error(FATAL, "Rotation of tracer BC fluxes not yet implemented.")

  ! Scalars and flags
  fluxes%accumulate_p_surf = fluxes_in%accumulate_p_surf

  fluxes%vPrecGlobalAdj = fluxes_in%vPrecGlobalAdj
  fluxes%saltFluxGlobalAdj = fluxes_in%saltFluxGlobalAdj
  fluxes%netFWGlobalAdj = fluxes_in%netFWGlobalAdj
  fluxes%vPrecGlobalScl = fluxes_in%vPrecGlobalScl
  fluxes%saltFluxGlobalScl = fluxes_in%saltFluxGlobalScl
  fluxes%netFWGlobalScl = fluxes_in%netFWGlobalScl

  fluxes%fluxes_used = fluxes_in%fluxes_used
  fluxes%dt_buoy_accum = fluxes_in%dt_buoy_accum
  fluxes%C_p = fluxes_in%C_p
  ! NOTE: gustless_accum_bug is set during allocation

  fluxes%num_msg = fluxes_in%num_msg
  fluxes%max_msg = fluxes_in%max_msg
end subroutine rotate_forcing

!< Rotate the forcing fields from the input domain
subroutine rotate_mech_forcing(forces_in, turns, forces)
  type(mech_forcing), intent(in)  :: forces_in  !< Forcing on the input domain
  integer, intent(in) :: turns                  !< Number of quarter-turns
  type(mech_forcing), intent(inout) :: forces   !< Forcing on the rotated domain

  logical :: do_stress, do_ustar, do_shelf, do_press, do_iceberg

  call get_mech_forcing_groups(forces_in, do_stress, do_ustar, do_shelf, &
                              do_press, do_iceberg)

  if (do_stress) &
    call rotate_vector(forces_in%taux, forces_in%tauy, turns, &
        forces%taux, forces%tauy)

  if (do_ustar) &
    call rotate_array(forces_in%ustar, turns, forces%ustar)

  if (do_shelf) then
    call rotate_array_pair( &
      forces_in%rigidity_ice_u, forces_in%rigidity_ice_v, turns, &
      forces%rigidity_ice_u, forces%rigidity_ice_v &
    )
    call rotate_array_pair( &
      forces_in%frac_shelf_u, forces_in%frac_shelf_v, turns, &
      forces%frac_shelf_u, forces%frac_shelf_v &
    )
  endif

  if (do_press) then
    ! NOTE: p_surf_SSH either points to p_surf or p_surf_full
    call rotate_array(forces_in%p_surf, turns, forces%p_surf)
    call rotate_array(forces_in%p_surf_full, turns, forces%p_surf_full)
    call rotate_array(forces_in%net_mass_src, turns, forces%net_mass_src)
  endif

  if (do_iceberg) then
    call rotate_array(forces_in%area_berg, turns, forces%area_berg)
    call rotate_array(forces_in%mass_berg, turns, forces%mass_berg)
  endif

  ! Copy fields
  forces%dt_force_accum = forces_in%dt_force_accum
  forces%net_mass_src_set = forces_in%net_mass_src_set
  forces%accumulate_p_surf = forces_in%accumulate_p_surf
  forces%accumulate_rigidity = forces_in%accumulate_rigidity
  forces%initialized = forces_in%initialized
end subroutine rotate_mech_forcing

!> \namespace mom_forcing_type
!!
!! \section section_fluxes Boundary fluxes
!!
!! The ocean is a forced-dissipative system. Forcing occurs at the
!! boundaries, and this module mediates the various forcing terms
!! from momentum, heat, salt, and mass.  Boundary fluxes from other
!! tracers are treated by coupling to biogeochemical models. We
!! here present elements of how MOM6 assumes boundary fluxes are
!! passed into the ocean.
!!
!! Note that all fluxes are positive into the ocean. For surface
!! boundary fluxes, that means fluxes are positive downward.
!! For example, a positive shortwave flux warms the ocean.
!!
!! \subsection subsection_momentum_fluxes Surface boundary momentum fluxes
!!
!! The ocean surface exchanges momentum with the overlying atmosphere,
!! sea ice, and land ice. The momentum is exchanged as a horizontal
!! stress (Newtons per squared meter: N/m2) imposed on the upper ocean
!! grid cell.
!!
!! \subsection subsection_mass_fluxes Surface boundary mass fluxes
!!
!! The ocean gains or loses mass through evaporation, precipitation,
!! sea ice melt/form, and and river runoff.  Positive mass fluxes
!! add mass to the liquid ocean. The boundary mass flux units are
!! (kilogram per square meter per sec: kg/(m2/sec)).
!!
!! * Evaporation field can in fact represent a
!!   mass loss (evaporation) or mass gain (condensation in foggy areas).
!! * sea ice formation leads to mass moving from the liquid ocean to the
!!   ice model, and melt adds liquid to the ocean.
!! * Precipitation can be liquid or frozen (snow). Furthermore, in
!!   some versions of the GFDL coupler, precipitation can be negative.
!!   The reason is that the ice model combines precipitation with
!!   ice melt and ice formation. This limitation of the ice model
!!   diagnostics should be overcome future versions.
!! * River runoff can be liquid or frozen.  Frozen runoff is often
!!   associated with calving land-ice and/or ice bergs.
!!
!! \subsection subsection_salt_fluxes Surface boundary salt fluxes
!!
!! Over most of the ocean, there is no exchange of salt with the
!! atmosphere. However, the liquid ocean exchanges salt with sea ice.
!! When ice forms, it extracts salt from ice pockets and discharges the
!! salt into the liquid ocean. The salt concentration of sea ice
!! is therefore much lower (around 5ppt) than liquid seawater
!! (around 30-35ppt in high latitudes).
!!
!! For ocean-ice models run with a prescribed atmosphere, such as
!! in the CORE/OMMIP simulations, it is necessary to employ a surface
!! restoring term to the k=1 salinity equation, thus imposing a salt
!! flux onto the ocean even outside of sea ice regimes.  This salt
!! flux is non-physical, and represents a limitation of the ocean-ice
!! models run without an interactive atmosphere.  Sometimes this salt
!! flux is converted to an implied fresh water flux.  However, doing
!! so generally leads to changes in the sea level, unless a global
!! normalization is provided to zero-out the net water flux.
!! As a complement, for models with a restoring salt flux, one may
!! choose to zero-out the net salt entering the ocean. There are
!! pros/cons of each approach.
!!
!!
!! \subsection subsection_heat_fluxes Surface boundary heat fluxes
!!
!!  There are many terms that contribute to boundary-related heating
!!  of the k=1 surface model grid cell. We here outline details of
!!  this heat, with each term having units W/m2.
!!
!!  The net flux of heat crossing ocean surface is stored in the diagnostic
!!  array "hfds".  This array is computed as
!! \f[
!!  \mbox{hfds = shortwave + longwave + latent + sensible + mass transfer + frazil + restore + flux adjustments}
!! \f]
!!
!!  * shortwave (SW)  = shortwave radiation (always warms ocean)
!!  * longwave (LW)   = longwave radiation (generally cools ocean)
!!  * latent (LAT)    = turbulent latent heat loss due to evaporation
!!                      (liquid to vapor) or melt (snow to liquid); generally
!!                      cools the ocean
!!  * sensible (SENS) = turbulent heat transfer due to differences in
!!                      air-sea or ice-sea temperature
!!  * mass transfer (MASS) = heat transfer due to heat content of mass (e.g., E-P+R)
!!                      transferred across ocean surface; computed relative
!!                      to 0 Celsius
!!  * frazil (FRAZ)   = heat transferred to form frazil sea ice
!!                      (positive heating of liquid ocean)
!!  * restore (RES)   = heat from surface damping sometimes imposed
!!                      in non-coupled model simulations .
!!  * restore (flux adjustments)   = heat from surface flux adjustment.
!!
!!  \subsubsection subsubsection_SW Treatment of shortwave
!!
!!  The shortwave field itself is split into two pieces:
!!
!!  * shortwave       = penetrative SW + non-penetrative SW
!!  * non-penetrative = non-downwelling shortwave; portion of SW
!!                      totally absorbed in the k=1 cell.
!!                      The non-penetrative SW is combined with
!!                      LW+LAT+SENS+seaice_melt_heat in net_heat inside routine
!!                      extractFluxes1d. Notably, for many cases,
!!                      non-penetrative SW = 0.
!!  * penetrative     = that portion of shortwave penetrating below
!!                      a tiny surface layer. This is the downwelling
!!                      shortwave. Penetrative SW participates in
!!                      the penetrative SW heating of k=1,nz cells,
!!                      with the amount of penetration dependent on
!!                      optical properties.
!!
!! \subsubsection subsubsection_bdy_heating Convergence of heat into the k=1 cell
!!
!! The convergence of boundary-related heat into surface grid cell is
!! given by the difference in the net heat entering the top of the k=1
!! cell and the penetrative SW leaving the bottom of the cell.
!! \f{eqnarray*}
!!  Q(k=1) &=& \mbox{hfds} - \mbox{pen}\_\mbox{SW(leaving bottom of k=1)}
!!   \\    &=& \mbox{nonpen}\_\mbox{SW} + (\mbox{pen}\_\mbox{SW(enter k=1)}-\mbox{pen}\_\mbox{SW(leave k=1)})
!!                              + \mbox{LW+LAT+SENS+MASS+FRAZ+RES}
!!   \\    &=& \mbox{nonpen}\_\mbox{SW}+ \mbox{LW+LAT+SENS+MASS+FRAZ+RES}
!!                + [\mbox{pen}\_\mbox{SW(enter k=1)} - \mbox{pen}\_\mbox{SW(leave k=1)}]
!!   \f}
!! The convergence of the penetrative shortwave flux is given by
!! \f$ \mbox{pen}\_\mbox{SW (enter k)}-\mbox{pen}\_\mbox{SW (leave k)}\f$.  This term
!! appears for all cells k=1,nz.  It is diagnosed as "rsdoabsorb" inside module
!! MOM6/src/parameterizations/vertical/MOM_diabatic_aux.F90
!!

end module MOM_forcing_type
