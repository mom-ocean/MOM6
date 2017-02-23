!> This module implements boundary forcing for MOM6.
module MOM_forcing_type

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,     only : hchksum, uchksum, vchksum
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, register_scalar_field
use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_alloc, query_averaging_enabled
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_EOS,           only : calculate_density_derivs
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_shortwave_abs, only : sumSWoverBands, optics_type
use MOM_spatial_means, only : global_area_integral, global_area_mean
use MOM_variables,     only : surface, thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

use coupler_types_mod, only : coupler_2d_bc_type

implicit none ; private

#include <MOM_memory.h>

public extractFluxes1d, extractFluxes2d, MOM_forcing_chksum, optics_type
public calculateBuoyancyFlux1d, calculateBuoyancyFlux2d, forcing_accumulate
public forcing_SinglePointPrint, mech_forcing_diags, forcing_diagnostics
public register_forcing_type_diags, allocate_forcing_type, deallocate_forcing_type

!> Structure that contains pointers to the boundary forcing
!! used to drive the liquid ocean simulated by MOM.
!! Data in this type is allocated in the module
!! MOM_surface_forcing.F90, of which there are three:
!! solo, coupled, and ice-shelf. Alternatively, they are
!! allocated in MESO_surface_forcing.F90, which is a
!! special case of solo_driver/MOM_surface_forcing.F90.
type, public :: forcing

  ! Pointers in this module should be initialized to NULL.

  ! surface stress components and turbulent velocity scale
  real, pointer, dimension(:,:) :: &
  taux          => NULL(), & !< zonal wind stress (Pa)
  tauy          => NULL(), & !< meridional wind stress (Pa)
  ustar         => NULL()    !< surface friction velocity scale (m/s)

  ! surface buoyancy force
  real, pointer, dimension(:,:) :: &
  buoy          => NULL()  !< buoyancy flux (m^2/s^3)

  ! radiative heat fluxes into the ocean (W/m^2)
  real, pointer, dimension(:,:) :: &
  sw            => NULL(), & !< shortwave (W/m^2)
  sw_vis_dir    => NULL(), & !< visible, direct shortwave (W/m^2)
  sw_vis_dif    => NULL(), & !< visible, diffuse shortwave (W/m^2)
  sw_nir_dir    => NULL(), & !< near-IR, direct shortwave (W/m^2)
  sw_nir_dif    => NULL(), & !< near-IR, diffuse shortwave (W/m^2)
  lw            => NULL()    !< longwave (W/m^2) (typically negative)

  ! turbulent heat fluxes into the ocean (W/m^2)
  real, pointer, dimension(:,:) :: &
  latent         => NULL(), & !< latent (W/m^2)   (typically < 0)
  sens           => NULL(), & !< sensible (W/m^2) (typically negative)
  heat_added     => NULL()    !< additional heat flux from SST restoring or flux adjustments (W/m^2)

  ! components of latent heat fluxes used for diagnostic purposes
  real, pointer, dimension(:,:) :: &
  latent_evap_diag    => NULL(), & !< latent (W/m^2) from evaporating liquid water (typically < 0)
  latent_fprec_diag   => NULL(), & !< latent (W/m^2) from melting fprec  (typically < 0)
  latent_frunoff_diag => NULL()    !< latent (W/m^2) from melting frunoff (calving) (typically < 0)

  ! water mass fluxes into the ocean ( kg/(m^2 s) ); these fluxes impact the ocean mass
  real, pointer, dimension(:,:) :: &
  evap          => NULL(), & !< (-1)*fresh water flux evaporated out of the ocean ( kg/(m^2 s) )
  lprec         => NULL(), & !< precipitating liquid water into the ocean ( kg/(m^2 s) )
  fprec         => NULL(), & !< precipitating frozen water into the ocean ( kg/(m^2 s) )
  vprec         => NULL(), & !< virtual liquid precip associated w/ SSS restoring ( kg/(m^2 s) )
  lrunoff       => NULL(), & !< liquid river runoff entering ocean ( kg/(m^2 s) )
  frunoff       => NULL(), & !< frozen river runoff (calving) entering ocean ( kg/(m^2 s) )
  seaice_melt   => NULL(), & !< seaice melt (positive) or formation (negative) ( kg/(m^2 s) )
  netMassIn     => NULL(), & !< Sum of water mass flux out of the ocean ( kg/(m^2 s) )
  netMassOut    => NULL(), & !< Net water mass flux into of the ocean ( kg/(m^2 s) )
  netSalt       => NULL()    !< Net salt entering the ocean

  ! heat associated with water crossing ocean surface
  real, pointer, dimension(:,:) :: &
  heat_content_cond    => NULL(), & !< heat content associated with condensating water (W/m^2)
  heat_content_lprec   => NULL(), & !< heat content associated with liquid >0 precip   (W/m^2) (diagnostic)
  heat_content_fprec   => NULL(), & !< heat content associated with frozen precip      (W/m^2)
  heat_content_vprec   => NULL(), & !< heat content associated with virtual >0 precip  (W/m^2)
  heat_content_lrunoff => NULL(), & !< heat content associated with liquid runoff      (W/m^2)
  heat_content_frunoff => NULL(), & !< heat content associated with frozen runoff      (W/m^2)
  heat_content_icemelt => NULL(), & !< heat content associated with liquid sea ice     (W/m^2)
  heat_content_massout => NULL(), & !< heat content associated with mass leaving ocean (W/m^2)
  heat_content_massin  => NULL()    !< heat content associated with mass entering ocean (W/m^2)

  ! salt mass flux (contributes to ocean mass only if non-Bouss )
  real, pointer, dimension(:,:) :: &
  salt_flux         => NULL(), & !< net salt flux into the ocean ( kg salt/(m^2 s) )
  salt_flux_in      => NULL(), & !< salt flux provided to the ocean from coupler ( kg salt/(m^2 s) )
  salt_flux_added => NULL()    !< additional salt flux from restoring or flux adjustment before adjustment
                                 !! to net zero ( kg salt/(m^2 s) )

  ! applied surface pressure from other component models (e.g., atmos, sea ice, land ice)
  real, pointer, dimension(:,:) :: &
  p_surf_full   => NULL(), & !< Pressure at the top ocean interface (Pa).
                             !! if there is sea-ice, then p_surf_flux is at ice-ocean interface
  p_surf        => NULL(), & !< Pressure at the top ocean interface (Pa) as used
                             !! to drive the ocean model. If p_surf is limited,
                             !! p_surf may be smaller than p_surf_full,
                             !! otherwise they are the same.
  p_surf_SSH    => NULL()    !< Pressure at the top ocean interface that is used
                             !! in corrections to the sea surface height field
                             !! that is passed back to the calling routines.
                             !! This may point to p_surf or to p_surf_full.

  ! tide related inputs
  real, pointer, dimension(:,:) :: &
  TKE_tidal     => NULL(), & !< tidal energy source driving mixing in bottom boundary layer (W/m^2)
  ustar_tidal   => NULL()    !< tidal contribution to bottom ustar (m/s)

  ! iceberg related inputs
  real, pointer, dimension(:,:) :: &
  ustar_berg   => NULL(),&    !< iceberg contribution to top ustar (m/s)
  area_berg   => NULL(),&     !< area of ocean surface covered by icebergs (m2/m2)
  mass_berg   => NULL()     !< mass of icebergs (kg/m2)

  ! land ice-shelf related inputs
  real, pointer, dimension(:,:) :: &
  ustar_shelf   => NULL(), &   !< friction velocity under ice-shelves (m/s)
                               !! as computed by the ocean at the previous time step.
  frac_shelf_h  => NULL(), &   !< Fractional ice shelf coverage of h-, u-, and v-
  frac_shelf_u  => NULL(), &   !< cells, nondimensional from 0 to 1. These are only
  frac_shelf_v  => NULL(), &   !< associated if ice shelves are enabled, and are
                               !! exactly 0 away from shelves or on land.
  iceshelf_melt   => NULL(), & !< ice shelf melt rate (positive) or freezing (negative) ( m/year )
  rigidity_ice_u => NULL(),&   !< Depth-integrated lateral viscosity of ice
  rigidity_ice_v => NULL()     !< shelves or sea ice at u- or v-points (m3/s)

  ! Scalars set by surface forcing modules
  real :: vPrecGlobalAdj     !< adjustment to restoring vprec to zero out global net ( kg/(m^2 s) )
  real :: saltFluxGlobalAdj  !< adjustment to restoring salt flux to zero out global net ( kg salt/(m^2 s) )
  real :: netFWGlobalAdj     !< adjustment to net fresh water to zero out global net ( kg/(m^2 s) )
  real :: vPrecGlobalScl     !< scaling of restoring vprec to zero out global net ( -1..1 )
  real :: saltFluxGlobalScl  !< scaling of restoring salt flux to zero out global net ( -1..1 )
  real :: netFWGlobalScl     !< scaling of net fresh water to zero out global net ( -1..1 )

  logical :: fluxes_used = .true. !< If true, all of the heat, salt, and mass
                                  !! fluxes have been applied to the ocean.
  real :: dt_buoy_accum  = -1.0   !< The amount of time over which the buoyancy fluxes
                                  !! should be applied, in s.  If negative, this forcing
                                  !! type variable has not yet been inialized.

  ! heat capacity
  real :: C_p                !< heat capacity of seawater ( J/(K kg) ).
                             !! C_p is is the same value as in thermovar_ptrs_type.

  ! passive tracer surface fluxes
  type(coupler_2d_bc_type), pointer :: tr_fluxes => NULL() !< This structure
     !! may contain an array of named fields used for passive tracer fluxes.
     !! All arrays in tr_fluxes use the coupler indexing, which has no halos.
     !! This is not a convenient convention, but imposed on MOM6 by the coupler.

  ! For internal error tracking
  integer :: num_msg = 0 !< Number of messages issues about excessive SW penetration
  integer :: max_msg = 2 !< Maximum number of messages to issue about excessive SW penetration

end type forcing

!> Structure that defines the id handles for the forcing type
type, public :: forcing_diags

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

  ! momentum flux diagnostic handls
  integer :: id_taux  = -1
  integer :: id_tauy  = -1
  integer :: id_ustar = -1

  integer :: id_psurf     = -1
  integer :: id_TKE_tidal = -1
  integer :: id_buoy      = -1

  ! clock id handle
  integer :: id_clock_forcing

  ! iceberg id handle
  integer :: id_ustar_berg = -1
  integer :: id_area_berg = -1
  integer :: id_mass_berg = -1

  !Iceberg + Ice shelf
  integer :: id_ustar_ice_cover = -1
  integer :: id_frac_ice_cover = -1

end type forcing_diags

contains

!> This subroutine extracts fluxes from the surface fluxes type. It works on a j-row
!! for optimization purposes. The 2d (i,j) wrapper is the next subroutine below.
!! This routine multiplies fluxes by dt, so that the result is an accumulation of fluxes
!! over a time step.
subroutine extractFluxes1d(G, GV, fluxes, optics, nsw, j, dt,                           &
                  DepthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
                  h, T, netMassInOut, netMassOut, net_heat, net_salt, pen_SW_bnd, tv,   &
                  aggregate_FW_forcing, nonpenSW)

  type(ocean_grid_type),             intent(in)    :: G                        !< ocean grid structure
  type(verticalGrid_type),           intent(in)    :: GV                       !< ocean vertical grid structure
  type(forcing),                     intent(inout) :: fluxes                   !< structure containing pointers to possible
                                                                               !! forcing fields. NULL unused fields.
  type(optics_type),                 pointer       :: optics                   !< pointer to optics
  integer,                           intent(in)    :: nsw                      !< number of bands of penetrating SW
  integer,                           intent(in)    :: j                        !< j-index to work on
  real,                              intent(in)    :: dt                       !< time step in seconds
  real,                              intent(in)    :: DepthBeforeScalingFluxes !< min ocean depth before scale away fluxes (H)
  logical,                           intent(in)    :: useRiverHeatContent      !< logical for river heat content
  logical,                           intent(in)    :: useCalvingHeatContent    !< logical for calving heat content
  real, dimension(SZI_(G),SZK_(G)),  intent(in)    :: h                        !< layer thickness (in H units)
  real, dimension(SZI_(G),SZK_(G)),  intent(in)    :: T                        !< layer temperatures (deg C)
  real, dimension(SZI_(G)),          intent(out)   :: netMassInOut             !<  net mass flux (non-Bouss) or volume flux
                                                                               !! (if Bouss) of water in/out of ocean over
                                                                               !! a time step (H units)
  real, dimension(SZI_(G)),          intent(out)   :: netMassOut               !< net mass flux (non-Bouss) or volume flux
                                                                               !! (if Bouss) of water leaving ocean surface
                                                                               !! over a time step (H units).
                                                                               !! netMassOut < 0 means mass leaves ocean.
  real, dimension(SZI_(G)),          intent(out)   :: net_heat                 !< net heat at the surface accumulated over a
                                                                               !! time step for coupler + restoring.
                                                                               !! Exclude two terms from net_heat:
                                                                               !! (1) downwelling (penetrative) SW,
                                                                               !! (2) evaporation heat content,
                                                                               !! (since do not yet know evap temperature).
                                                                               !! Units of net_heat are (K * H).
  real, dimension(SZI_(G)),          intent(out)   :: net_salt                 !< surface salt flux into the ocean accumulated
                                                                               !! over a time step (ppt * H)
  real, dimension(:,:),              intent(out)   :: pen_SW_bnd               !< penetrating SW flux, split into bands.
                                                                               !! Units are (deg K * H) and array size
                                                                               !! nsw x SZI_(G), where nsw=number of SW bands
                                                                               !! in pen_SW_bnd. This heat flux is not part
                                                                               !! of net_heat.
  type(thermo_var_ptrs),             intent(inout) :: tv                       !< structure containing pointers to available
                                                                               !! thermodynamic fields. Used to keep
                                                                               !! track of the heat flux associated with net
                                                                               !! mass fluxes into the ocean.
  logical,                           intent(in)    :: aggregate_FW_forcing     !< For determining how to aggregate forcing.
  real, dimension(SZI_(G)), optional, intent(out)  :: nonpenSW                 !< non-downwelling SW; use in net_heat.
                                                                               !! Sum over SW bands when diagnosing nonpenSW.
                                                                               !! Units are (K * H).

  ! local
  real :: htot(SZI_(G))       ! total ocean depth (m for Bouss or kg/m^2 for non-Bouss)
  real :: Pen_sw_tot(SZI_(G)) ! sum across all bands of Pen_SW (K * H)
  real :: Ih_limit            ! inverse depth at which surface fluxes start to be limited (1/H)
  real :: scale               ! scale scales away fluxes if depth < DepthBeforeScalingFluxes
  real :: J_m2_to_H           ! converts J/m^2 to H units (m for Bouss and kg/m^2 for non-Bouss)
  real :: Irho0               ! 1.0 / Rho0
  real :: I_Cp                ! 1.0 / C_p

  character(len=200) :: mesg
  integer            :: is, ie, nz, i, k, n
  Ih_limit  = 1.0 / DepthBeforeScalingFluxes
  Irho0     = 1.0 / GV%Rho0
  I_Cp      = 1.0 / fluxes%C_p
  J_m2_to_H = 1.0 / (GV%H_to_kg_m2 * fluxes%C_p)

  is = G%isc ; ie = G%iec ; nz = G%ke


  ! error checking

  if (nsw > 0) then ; if (nsw /= optics%nbands) call MOM_error(WARNING, &
    "mismatch in the number of bands of shortwave radiation in MOM_forcing_type extract_fluxes.")
  endif

  if (.not.ASSOCIATED(fluxes%sw)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%sw is not associated.")

  if (.not.ASSOCIATED(fluxes%lw)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%lw is not associated.")

  if (.not.ASSOCIATED(fluxes%latent)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%latent is not associated.")

  if (.not.ASSOCIATED(fluxes%sens)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%sens is not associated.")

  if (.not.ASSOCIATED(fluxes%evap)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: No evaporation defined.")

  if (.not.ASSOCIATED(fluxes%vprec)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%vprec not defined.")

  if ((.not.ASSOCIATED(fluxes%lprec)) .or. &
      (.not.ASSOCIATED(fluxes%fprec))) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: No precipitation defined.")

  do i=is,ie ; htot(i) = h(i,1) ; enddo
  do k=2,nz ; do i=is,ie ; htot(i) = htot(i) + h(i,k) ; enddo ; enddo


  do i=is,ie

    scale = 1.0
    if (htot(i)*Ih_limit < 1.0) scale = htot(i)*Ih_limit

    ! Convert the penetrating shortwave forcing to (K * H)
    ! (H=m for Bouss, H=kg/m2 for non-Bouss)
    Pen_sw_tot(i) = 0.0
    if (nsw >= 1) then
      do n=1,nsw
        Pen_SW_bnd(n,i) = J_m2_to_H*scale*dt * max(0.0, optics%sw_pen_band(n,i,j))
        Pen_sw_tot(i)   = Pen_sw_tot(i) + Pen_SW_bnd(n,i)
      enddo
    else
      Pen_SW_bnd(1,i) = 0.0
    endif

    ! net volume/mass of liquid and solid passing through surface boundary fluxes
    netMassInOut(i) = dt * (scale * ((((( fluxes%lprec(i,j)      &
                                        + fluxes%fprec(i,j)   )  &
                                        + fluxes%evap(i,j)    )  &
                                        + fluxes%lrunoff(i,j) )  &
                                        + fluxes%vprec(i,j)   )  &
                                        + fluxes%frunoff(i,j) ) )

    ! smg:
    ! for non-Bouss, we add/remove salt mass to total ocean mass. to conserve
    ! total salt mass ocean+ice, the sea ice model must lose mass when
    ! salt mass is added to the ocean, which may still need to be coded.
    if (.not.GV%Boussinesq .and. ASSOCIATED(fluxes%salt_flux)) then
      netMassInOut(i) = netMassInOut(i) + (dt * GV%kg_m2_to_H) * (scale * fluxes%salt_flux(i,j))
    endif

    ! net volume/mass of water leaving the ocean.
    ! check that fluxes are < 0, which means mass is indeed leaving.
    netMassOut(i) = 0.0

    ! evap > 0 means condensating water is added into ocean.
    ! evap < 0 means evaporation of water from the ocean, in
    ! which case heat_content_evap is computed in MOM_diabatic_driver.F90
    if(fluxes%evap(i,j) < 0.0) then
      netMassOut(i) = netMassOut(i) + fluxes%evap(i,j)
  !   if(ASSOCIATED(fluxes%heat_content_cond)) fluxes%heat_content_cond(i,j) = 0.0 !??? --AJA
    endif

    ! lprec < 0 means sea ice formation taking water from the ocean.
    ! smg: we should split the ice melt/formation from the lprec
    if(fluxes%lprec(i,j) < 0.0) then
      netMassOut(i) = netMassOut(i) + fluxes%lprec(i,j)
    endif

    ! vprec < 0 means virtual evaporation arising from surface salinity restoring,
    ! in which case heat_content_vprec is computed in MOM_diabatic_driver.F90.
    if(fluxes%vprec(i,j) < 0.0) then
      netMassOut(i) = netMassOut(i) + fluxes%vprec(i,j)
    endif
    netMassOut(i) = dt * scale * netMassOut(i)

    ! convert to H units (Bouss=meter or non-Bouss=kg/m^2)
    netMassInOut(i) = GV%kg_m2_to_H * netMassInOut(i)
    netMassOut(i)   = GV%kg_m2_to_H * netMassOut(i)

    ! surface heat fluxes from radiation and turbulent fluxes (K * H)
    ! (H=m for Bouss, H=kg/m2 for non-Bouss)
    net_heat(i) = scale * dt * J_m2_to_H * &
                  ( fluxes%sw(i,j) +  ((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)) )

    ! Add heat flux from surface damping (restoring) (K * H) or flux adjustments.
    if (ASSOCIATED(fluxes%heat_added)) then
       net_heat(i) = net_heat(i) + (scale * (dt * J_m2_to_H)) * fluxes%heat_added(i,j)
    endif

    ! Add explicit heat flux for runoff (which is part of the ice-ocean boundary
    ! flux type). Runoff is otherwise added with a temperature of SST.
    if (useRiverHeatContent) then
      ! remove lrunoff*SST here, to counteract its addition elsewhere
      net_heat(i) = (net_heat(i) + (scale*(dt*J_m2_to_H)) * fluxes%heat_content_lrunoff(i,j)) - &
                     (GV%kg_m2_to_H * (scale * dt)) * fluxes%lrunoff(i,j) * T(i,1)
      if (ASSOCIATED(tv%TempxPmE)) then
        tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + (scale * dt) * &
            (I_Cp*fluxes%heat_content_lrunoff(i,j) - fluxes%lrunoff(i,j)*T(i,1))
      endif
    endif

    ! Add explicit heat flux for calving (which is part of the ice-ocean boundary
    ! flux type). Calving is otherwise added with a temperature of SST.
    if (useCalvingHeatContent) then
      ! remove frunoff*SST here, to counteract its addition elsewhere
      net_heat(i) = net_heat(i) + (scale*(dt*J_m2_to_H)) * fluxes%heat_content_frunoff(i,j) - &
                    (GV%kg_m2_to_H * (scale * dt)) * fluxes%frunoff(i,j) * T(i,1)
      if (ASSOCIATED(tv%TempxPmE)) then
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
!    if (ASSOCIATED(fluxes%heat_content_lprec)) then
!      net_heat(i) = net_heat(i) + scale * dt * J_m2_to_H *                    &
!     (fluxes%heat_content_lprec(i,j)    + (fluxes%heat_content_fprec(i,j)   + &
!     (fluxes%heat_content_lrunoff(i,j)  + (fluxes%heat_content_frunoff(i,j) + &
!     (fluxes%heat_content_cond(i,j)     +  fluxes%heat_content_vprec(i,j))))))
!    endif

    if (fluxes%num_msg < fluxes%max_msg) then
      if (Pen_SW_tot(i) > 1.000001*J_m2_to_H*scale*dt*fluxes%sw(i,j)) then
        fluxes%num_msg = fluxes%num_msg + 1
        write(mesg,'("Penetrating shortwave of ",1pe17.10, &
                    &" exceeds total shortwave of ",1pe17.10,&
                    &" at ",1pg11.4,"E, "1pg11.4,"N.")') &
               Pen_SW_tot(i),J_m2_to_H*scale*dt*fluxes%sw(i,j),&
               G%geoLonT(i,j),G%geoLatT(i,j)
        call MOM_error(WARNING,mesg)
      endif
    endif

    ! remove penetrative portion of the SW that is NOT absorbed within a
    ! tiny layer at the top of the ocean.
    net_heat(i) = net_heat(i) - Pen_SW_tot(i)

    ! diagnose non-downwelling SW
    if (present(nonPenSW)) then
      nonPenSW(i) = scale * dt * J_m2_to_H * fluxes%sw(i,j) - Pen_SW_tot(i)
    endif

    ! Salt fluxes
    Net_salt(i) = 0.0
    ! Convert salt_flux from kg (salt)/(m^2 * s) to
    ! Boussinesq: (ppt * m)
    ! non-Bouss:  (g/m^2)
    if (ASSOCIATED(fluxes%salt_flux)) then
      Net_salt(i) = (scale * dt * (1000.0 * fluxes%salt_flux(i,j))) * GV%kg_m2_to_H
      fluxes%netSalt(i,j) = Net_salt(i)
    endif
    ! Diagnostics follow...

    ! Initialize heat_content_massin that is diagnosed in mixedlayer_convection or
    ! applyBoundaryFluxes such that the meaning is as the sum of all incoming components.
    if(ASSOCIATED(fluxes%heat_content_massin))  then
      if (aggregate_FW_forcing) then
        if (netMassInOut(i) > 0.0) then ! net is "in"
          fluxes%heat_content_massin(i,j) = -fluxes%C_p * netMassOut(i) * T(i,1) * GV%H_to_kg_m2 / dt
        else ! net is "out"
          fluxes%heat_content_massin(i,j) = fluxes%C_p * ( netMassInout(i) - netMassOut(i) ) * T(i,1) * GV%H_to_kg_m2 / dt
        endif
      else
        fluxes%heat_content_massin(i,j) = 0.
      endif
    endif

    ! Initialize heat_content_massout that is diagnosed in mixedlayer_convection or
    ! applyBoundaryFluxes such that the meaning is as the sum of all outgoing components.
    if(ASSOCIATED(fluxes%heat_content_massout)) then
      if (aggregate_FW_forcing) then
        if (netMassInOut(i) > 0.0) then ! net is "in"
          fluxes%heat_content_massout(i,j) = fluxes%C_p * netMassOut(i) * T(i,1) * GV%H_to_kg_m2 / dt
        else ! net is "out"
          fluxes%heat_content_massout(i,j) = -fluxes%C_p * ( netMassInout(i) - netMassOut(i) ) * T(i,1) * GV%H_to_kg_m2 / dt
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
    if(ASSOCIATED(fluxes%heat_content_lprec)) then
      if (fluxes%lprec(i,j) > 0.0) then
        fluxes%heat_content_lprec(i,j) = fluxes%C_p*fluxes%lprec(i,j)*T(i,1)
      else
        fluxes%heat_content_lprec(i,j) = 0.0
      endif
    endif

    ! fprec SHOULD enter ocean at 0degC if atmos model does not provide fprec heat content.
    ! However, we need to adjust netHeat above to reflect the difference between 0decC and SST
    ! and until we do so fprec is treated like lprec and enters at SST. -AJA
    if(ASSOCIATED(fluxes%heat_content_fprec)) then
      if (fluxes%fprec(i,j) > 0.0) then
        fluxes%heat_content_fprec(i,j) = fluxes%C_p*fluxes%fprec(i,j)*T(i,1)
      else
        fluxes%heat_content_fprec(i,j) = 0.0
      endif
    endif

    ! virtual precip associated with salinity restoring
    ! vprec > 0 means add water to ocean, assumed to be at SST
    ! vprec < 0 means remove water from ocean; set heat_content_vprec in MOM_diabatic_driver.F90
    if(ASSOCIATED(fluxes%heat_content_vprec)) then
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
    if(ASSOCIATED(fluxes%heat_content_cond)) then
      if (fluxes%evap(i,j) > 0.0) then
        fluxes%heat_content_cond(i,j) = fluxes%C_p*fluxes%evap(i,j)*T(i,1)
      else
        fluxes%heat_content_cond(i,j) = 0.0
      endif
    endif

    ! Liquid runoff enters ocean at SST if land model does not provide runoff heat content.
    if (.not. useRiverHeatContent) then
      if (ASSOCIATED(fluxes%lrunoff) .and. ASSOCIATED(fluxes%heat_content_lrunoff)) then
        fluxes%heat_content_lrunoff(i,j) = fluxes%C_p*fluxes%lrunoff(i,j)*T(i,1)
      endif
    endif

    ! Icebergs enter ocean at SST if land model does not provide calving heat content.
    if (.not. useCalvingHeatContent) then
      if (ASSOCIATED(fluxes%frunoff) .and. ASSOCIATED(fluxes%heat_content_frunoff)) then
        fluxes%heat_content_frunoff(i,j) = fluxes%C_p*fluxes%frunoff(i,j)*T(i,1)
      endif
    endif

  enddo ! i-loop

end subroutine extractFluxes1d


!> 2d wrapper for 1d extract fluxes from surface fluxes type.
!! This subroutine extracts fluxes from the surface fluxes type. It multiplies the
!! fluxes by dt, so that the result is an accumulation of the fluxes over a time step.
subroutine extractFluxes2d(G, GV, fluxes, optics, nsw, dt,                                  &
                  DepthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
                  h, T, netMassInOut, netMassOut, net_heat, Net_salt, Pen_SW_bnd, tv,   &
                  aggregate_FW_forcing)

  type(ocean_grid_type),                 intent(in)    :: G                        !< ocean grid structure
  type(verticalGrid_type),               intent(in)    :: GV                       !< ocean vertical grid structure
  type(forcing),                         intent(inout) :: fluxes                   !< structure containing pointers to forcing.
  type(optics_type),                     pointer       :: optics                   !< pointer to optics
  integer,                               intent(in)    :: nsw                      !< number of bands of penetrating SW
  real,                                  intent(in)    :: dt                       !< time step in seconds
  real,                                  intent(in)    :: DepthBeforeScalingFluxes !< min ocean depth before scale away fluxes (H)
  logical,                               intent(in)    :: useRiverHeatContent      !< logical for river heat content
  logical,                               intent(in)    :: useCalvingHeatContent    !< logical for calving heat content
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h                        !< layer thickness (in H units)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: T                        !< layer temperatures (deg C)
  real, dimension(SZI_(G),SZJ_(G)),        intent(out) :: netMassInOut             !<  net mass flux (non-Bouss) or volume flux
                                                                                   !! (if Bouss) of water in/out of ocean over
                                                                                   !! a time step (H units)
  real, dimension(SZI_(G),SZJ_(G)),        intent(out) :: netMassOut               !< net mass flux (non-Bouss) or volume flux
                                                                                   !! (if Bouss) of water leaving ocean surface
                                                                                   !! over a time step (H units).
  real, dimension(SZI_(G),SZJ_(G)),        intent(out) :: net_heat                 !< net heat at the surface accumulated over a
                                                                                   !! time step associated with coupler + restore.
                                                                                   !! Exclude two terms from net_heat:
                                                                                   !! (1) downwelling (penetrative) SW,
                                                                                   !! (2) evaporation heat content,
                                                                                   !! (since do not yet know temperature of evap).
                                                                                   !! Units of net_heat are (K * H).
  real, dimension(SZI_(G),SZJ_(G)),        intent(out) :: net_salt                 !< surface salt flux into the ocean accumulated
                                                                                   !! over a time step (ppt * H)
  real, dimension(:,:,:),                intent(out)   :: pen_SW_bnd               !< penetrating shortwave flux, split into bands.
                                                                                   !! Units (deg K * H) & array size nsw x SZI_(G),
                                                                                   !! where nsw=number of SW bands in pen_SW_bnd.
                                                                                   !! This heat flux is not in net_heat.
  type(thermo_var_ptrs),                 intent(inout) :: tv                       !< structure containing pointers to available
                                                                                   !! thermodynamic fields. Here it is used to keep
                                                                                   !! track of the heat flux associated with net
                                                                                   !! mass fluxes into the ocean.
  logical,                               intent(in)    :: aggregate_FW_forcing     !< For determining how to aggregate the forcing.


  integer :: j
!$OMP parallel do default(none) shared(G, GV, fluxes, optics, nsw,dt,DepthBeforeScalingFluxes, &
!$OMP                                  useRiverHeatContent, useCalvingHeatContent,             &
!$OMP                                  h,T,netMassInOut,netMassOut,Net_heat,Net_salt,Pen_SW_bnd,tv, &
!$OMP                                  aggregate_FW_forcing)
  do j=G%jsc, G%jec
    call extractFluxes1d(G, GV, fluxes, optics, nsw, j, dt,                      &
            DepthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent,&
            h(:,j,:), T(:,j,:), netMassInOut(:,j), netMassOut(:,j),              &
            net_heat(:,j), net_salt(:,j), pen_SW_bnd(:,:,j), tv, aggregate_FW_forcing)
  enddo

end subroutine extractFluxes2d


!> This routine calculates surface buoyancy flux by adding up the heat, FW & salt fluxes.
!! These are actual fluxes, with units of stuff per time. Setting dt=1 in the call to
!! extractFluxes routine allows us to get "stuf per time" rather than the time integrated
!! fluxes needed in other routines that call extractFluxes.
subroutine calculateBuoyancyFlux1d(G, GV, fluxes, optics, h, Temp, Salt, tv, j, &
                                   buoyancyFlux, netHeatMinusSW, netSalt )


  type(ocean_grid_type),                    intent(in)    :: G              !< ocean grid
  type(verticalGrid_type),                  intent(in)    :: GV             !< ocean vertical grid structure
  type(forcing),                            intent(inout) :: fluxes         !< surface fluxes
  type(optics_type),                        pointer       :: optics         !< penetrating SW optics
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h              !< layer thickness (H)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: Temp           !< prognostic temp(deg C)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: Salt           !< salinity (ppt)
  type(thermo_var_ptrs),                    intent(inout) :: tv             !< thermodynamics type
  integer,                                  intent(in)    :: j              !< j-row to work on
  real, dimension(SZI_(G),SZK_(G)+1),       intent(inout) :: buoyancyFlux   !< buoyancy flux (m^2/s^3)
  real, dimension(SZI_(G)),                 intent(inout) :: netHeatMinusSW !< surf Heat flux (K H/s)
  real, dimension(SZI_(G)),                 intent(inout) :: netSalt        !< surf salt flux (ppt H/s)

  ! local variables
  integer                                   :: nsw, start, npts, k
  real, parameter                           :: dt = 1.    ! to return a rate from extractFluxes1d
  real, dimension( SZI_(G) )                :: netH       ! net FW flux (m/s for Bouss)
  real, dimension( SZI_(G) )                :: netEvap    ! net FW flux leaving ocean via evaporation (m/s for Bouss)
  real, dimension( SZI_(G) )                :: netHeat    ! net temp flux (K m/s)
  real, dimension( optics%nbands, SZI_(G) ) :: penSWbnd   ! SW penetration bands
  real, dimension( SZI_(G) )                :: pressure   ! pressurea the surface (Pa)
  real, dimension( SZI_(G) )                :: dRhodT     ! density partial derivative wrt temp
  real, dimension( SZI_(G) )                :: dRhodS     ! density partial derivative wrt saln
  real, dimension(SZI_(G),SZK_(G)+1)        :: netPen

  logical :: useRiverHeatContent
  logical :: useCalvingHeatContent
  real    :: depthBeforeScalingFluxes, GoRho
  real    :: H_limit_fluxes

  nsw = optics%nbands

  !  smg: what do we do when have heat fluxes from calving and river?
  useRiverHeatContent   = .False.
  useCalvingHeatContent = .False.

  depthBeforeScalingFluxes = max( GV%Angstrom, 1.e-30*GV%m_to_H )
  pressure(:) = 0. ! Ignore atmospheric pressure
  GoRho       = GV%g_Earth / GV%Rho0
  start       = 1 + G%isc - G%isd
  npts        = 1 + G%iec - G%isc

  H_limit_fluxes = depthBeforeScalingFluxes

  ! The surface forcing is contained in the fluxes type.
  ! We aggregate the thermodynamic forcing for a time step into the following:
  ! netH       = water (H units/s) added/removed via surface fluxes
  ! netHeat    = heat (degC * H/s) via surface fluxes
  ! netSalt    = salt ( g(salt)/m2 for non-Bouss and ppt*m for Bouss /s) via surface fluxes
  ! Note that unlike other calls to extractFLuxes1d() that return the time-integrated flux
  ! this call returns the rate because dt=1
  call extractFluxes1d(G, GV, fluxes, optics, nsw, j, dt,                                 &
                depthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
                h(:,j,:), Temp(:,j,:), netH, netEvap, netHeatMinusSW,                 &
                netSalt, penSWbnd, tv, .false.)

  ! Sum over bands and attenuate as a function of depth
  ! netPen is the netSW as a function of depth
  call sumSWoverBands(G, GV, h(:,j,:), optics%opacity_band(:,:,j,:), nsw, j, dt, &
                      H_limit_fluxes, .true., penSWbnd, netPen)

  ! Density derivatives
  call calculate_density_derivs(Temp(:,j,1), Salt(:,j,1), pressure, &
                                dRhodT, dRhodS, start, npts, tv%eqn_of_state)

  ! Adjust netSalt to reflect dilution effect of FW flux
  netSalt(G%isc:G%iec) = netSalt(G%isc:G%iec) - Salt(G%isc:G%iec,j,1) * netH(G%isc:G%iec) * GV%H_to_m ! ppt H/s

  ! Add in the SW heating for purposes of calculating the net
  ! surface buoyancy flux affecting the top layer.
  !netHeat(:) = netHeatMinusSW(:) + sum( penSWbnd(:,:), dim=1 )
  netHeat(G%isc:G%iec) = netHeatMinusSW(G%isc:G%iec) + netPen(G%isc:G%iec,1) ! K H/s

  ! Convert to a buoyancy flux, excluding penetrating SW heating
  buoyancyFlux(G%isc:G%iec,1) = - GoRho * ( dRhodS(G%isc:G%iec) * netSalt(G%isc:G%iec) + &
                                             dRhodT(G%isc:G%iec) * netHeat(G%isc:G%iec) ) * GV%H_to_m ! m^2/s^3
  ! We also have a penetrative buoyancy flux associated with penetrative SW
  do k=2, G%ke+1
    buoyancyFlux(G%isc:G%iec,k) = - GoRho * ( dRhodT(G%isc:G%iec) * netPen(G%isc:G%iec,k) ) * GV%H_to_m ! m^2/s^3
  enddo

end subroutine calculateBuoyancyFlux1d


!> Calculates surface buoyancy flux by adding up the heat, FW and salt fluxes,
!! for 2d arrays.  This is a wrapper for calculateBuoyancyFlux1d.
subroutine calculateBuoyancyFlux2d(G, GV, fluxes, optics, h, Temp, Salt, tv, &
                                   buoyancyFlux, netHeatMinusSW, netSalt)
  type(ocean_grid_type),                      intent(in)    :: G              !< ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV             !< ocean vertical grid structure
  type(forcing),                              intent(inout) :: fluxes         !< surface fluxes
  type(optics_type),                          pointer       :: optics         !< SW ocean optics
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: h              !< layer thickness (H)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: Temp           !< temperature (deg C)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: Salt           !< salinity (ppt)
  type(thermo_var_ptrs),                      intent(inout) :: tv             !< thermodynamics type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: buoyancyFlux   !< buoy flux (m^2/s^3)
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(inout) :: netHeatMinusSW !< surf temp flux (K H)
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(inout) :: netSalt        !< surf salt flux (ppt H)

  ! local variables
  real, dimension( SZI_(G) ) :: netT ! net temperature flux (K m/s)
  real, dimension( SZI_(G) ) :: netS ! net saln flux (ppt m/s)
  integer :: j

  netT(G%isc:G%iec) = 0. ; netS(G%isc:G%iec) = 0.

!$OMP parallel do default(none) shared(G,GV,fluxes,optics,h,Temp,Salt,tv,buoyancyFlux,&
!$OMP                                  netHeatMinusSW,netSalt)                     &
!$OMP                     firstprivate(netT,netS)
  do j = G%jsc, G%jec
    call calculateBuoyancyFlux1d(G, GV, fluxes, optics, h, Temp, Salt, tv, j, buoyancyFlux(:,j,:), netT, netS )
    if (present(netHeatMinusSW)) netHeatMinusSW(G%isc:G%iec,j) = netT(G%isc:G%iec)
    if (present(netSalt)) netSalt(G%isc:G%iec,j) = netS(G%isc:G%iec)
  enddo ! j

end subroutine calculateBuoyancyFlux2d


!> Write out chksums for basic state variables.
subroutine MOM_forcing_chksum(mesg, fluxes, G, haloshift)
  character(len=*),        intent(in) :: mesg      !< message
  type(forcing),           intent(in) :: fluxes    !< fluxes type
  type(ocean_grid_type),   intent(in) :: G         !< grid type
  integer, optional,       intent(in) :: haloshift !< shift in halo

  integer :: is, ie, js, je, nz, hshift
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  hshift=1; if (present(haloshift)) hshift=haloshift

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  if (associated(fluxes%taux)) &
    call uchksum(fluxes%taux, mesg//" fluxes%taux",G%HI,haloshift=1)
  if (associated(fluxes%tauy)) &
    call vchksum(fluxes%tauy, mesg//" fluxes%tauy",G%HI,haloshift=1)
  if (associated(fluxes%ustar)) &
    call hchksum(fluxes%ustar, mesg//" fluxes%ustar",G%HI,haloshift=hshift)
  if (associated(fluxes%buoy)) &
    call hchksum(fluxes%buoy, mesg//" fluxes%buoy ",G%HI,haloshift=hshift)
  if (associated(fluxes%sw)) &
    call hchksum(fluxes%sw, mesg//" fluxes%sw",G%HI,haloshift=hshift)
  if (associated(fluxes%sw_vis_dir)) &
    call hchksum(fluxes%sw_vis_dir, mesg//" fluxes%sw_vis_dir",G%HI,haloshift=hshift)
  if (associated(fluxes%sw_vis_dif)) &
    call hchksum(fluxes%sw_vis_dif, mesg//" fluxes%sw_vis_dif",G%HI,haloshift=hshift)
  if (associated(fluxes%sw_nir_dir)) &
    call hchksum(fluxes%sw_nir_dir, mesg//" fluxes%sw_nir_dir",G%HI,haloshift=hshift)
  if (associated(fluxes%sw_nir_dif)) &
    call hchksum(fluxes%sw_nir_dif, mesg//" fluxes%sw_nir_dif",G%HI,haloshift=hshift)
  if (associated(fluxes%lw)) &
    call hchksum(fluxes%lw, mesg//" fluxes%lw",G%HI,haloshift=hshift)
  if (associated(fluxes%latent)) &
    call hchksum(fluxes%latent, mesg//" fluxes%latent",G%HI,haloshift=hshift)
  if (associated(fluxes%latent_evap_diag)) &
    call hchksum(fluxes%latent_evap_diag, mesg//" fluxes%latent_evap_diag",G%HI,haloshift=hshift)
  if (associated(fluxes%latent_fprec_diag)) &
    call hchksum(fluxes%latent_fprec_diag, mesg//" fluxes%latent_fprec_diag",G%HI,haloshift=hshift)
  if (associated(fluxes%latent_frunoff_diag)) &
    call hchksum(fluxes%latent_frunoff_diag, mesg//" fluxes%latent_frunoff_diag",G%HI,haloshift=hshift)
  if (associated(fluxes%sens)) &
    call hchksum(fluxes%sens, mesg//" fluxes%sens",G%HI,haloshift=hshift)
  if (associated(fluxes%evap)) &
    call hchksum(fluxes%evap, mesg//" fluxes%evap",G%HI,haloshift=hshift)
  if (associated(fluxes%lprec)) &
    call hchksum(fluxes%lprec, mesg//" fluxes%lprec",G%HI,haloshift=hshift)
  if (associated(fluxes%fprec)) &
    call hchksum(fluxes%fprec, mesg//" fluxes%fprec",G%HI,haloshift=hshift)
  if (associated(fluxes%vprec)) &
    call hchksum(fluxes%vprec, mesg//" fluxes%vprec",G%HI,haloshift=hshift)
  if (associated(fluxes%seaice_melt)) &
    call hchksum(fluxes%seaice_melt, mesg//" fluxes%seaice_melt",G%HI,haloshift=hshift)
  if (associated(fluxes%p_surf)) &
    call hchksum(fluxes%p_surf, mesg//" fluxes%p_surf",G%HI,haloshift=hshift)
  if (associated(fluxes%salt_flux)) &
    call hchksum(fluxes%salt_flux, mesg//" fluxes%salt_flux",G%HI,haloshift=hshift)
  if (associated(fluxes%TKE_tidal)) &
    call hchksum(fluxes%TKE_tidal, mesg//" fluxes%TKE_tidal",G%HI,haloshift=hshift)
  if (associated(fluxes%ustar_tidal)) &
    call hchksum(fluxes%ustar_tidal, mesg//" fluxes%ustar_tidal",G%HI,haloshift=hshift)
  if (associated(fluxes%lrunoff)) &
    call hchksum(fluxes%lrunoff, mesg//" fluxes%lrunoff",G%HI,haloshift=hshift)
  if (associated(fluxes%frunoff)) &
    call hchksum(fluxes%frunoff, mesg//" fluxes%frunoff",G%HI,haloshift=hshift)
  if (associated(fluxes%heat_content_lrunoff)) &
    call hchksum(fluxes%heat_content_lrunoff, mesg//" fluxes%heat_content_lrunoff",G%HI,haloshift=hshift)
  if (associated(fluxes%heat_content_frunoff)) &
    call hchksum(fluxes%heat_content_frunoff, mesg//" fluxes%heat_content_frunoff",G%HI,haloshift=hshift)
  if (associated(fluxes%heat_content_lprec)) &
    call hchksum(fluxes%heat_content_lprec, mesg//" fluxes%heat_content_lprec",G%HI,haloshift=hshift)
  if (associated(fluxes%heat_content_fprec)) &
    call hchksum(fluxes%heat_content_fprec, mesg//" fluxes%heat_content_fprec",G%HI,haloshift=hshift)
  if (associated(fluxes%heat_content_cond)) &
    call hchksum(fluxes%heat_content_cond, mesg//" fluxes%heat_content_cond",G%HI,haloshift=hshift)
  if (associated(fluxes%heat_content_massout)) &
    call hchksum(fluxes%heat_content_massout, mesg//" fluxes%heat_content_massout",G%HI,haloshift=hshift)
end subroutine MOM_forcing_chksum


!> Write out values of the fluxes arrays at the i,j location. This is a debugging tool.
subroutine forcing_SinglePointPrint(fluxes, G, i, j, mesg)
  type(forcing),         intent(in) :: fluxes !< Fluxes type
  type(ocean_grid_type), intent(in) :: G      !< Grid type
  character(len=*),      intent(in) :: mesg   !< Message
  integer,               intent(in) :: i      !< i-index
  integer,               intent(in) :: j      !< j-index

  write(0,'(2a)') 'MOM_forcing_type, forcing_SinglePointPrint: Called from ',mesg
  write(0,'(a,2es15.3)') 'MOM_forcing_type, forcing_SinglePointPrint: lon,lat = ',G%geoLonT(i,j),G%geoLatT(i,j)
  call locMsg(fluxes%taux,'taux')
  call locMsg(fluxes%tauy,'tauy')
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
subroutine register_forcing_type_diags(Time, diag, use_temperature, handles, use_berg_fluxes)
  type(time_type),     intent(in)    :: Time            !< time type
  type(diag_ctrl),     intent(inout) :: diag            !< diagnostic control type
  logical,             intent(in)    :: use_temperature !< True if T/S are in use
  type(forcing_diags), intent(inout) :: handles         !< handles for diagnostics
  logical, optional,   intent(in)    :: use_berg_fluxes !< If true, allow iceberg flux diagnostics

  ! Clock for forcing diagnostics
  handles%id_clock_forcing=cpu_clock_id('(Ocean forcing diagnostics)', grain=CLOCK_ROUTINE)


  handles%id_taux = register_diag_field('ocean_model', 'taux', diag%axesCu1, Time,  &
        'Zonal surface stress from ocean interactions with atmos and ice', 'Pascal',&
        standard_name='surface_downward_x_stress', cmor_field_name='tauuo',         &
        cmor_units='N m-2', cmor_long_name='Surface Downward X Stress',             &
        cmor_standard_name='surface_downward_x_stress')

  handles%id_tauy = register_diag_field('ocean_model', 'tauy', diag%axesCv1, Time,  &
        'Meridional surface stress ocean interactions with atmos and ice', 'Pascal',&
         standard_name='surface_downward_y_stress', cmor_field_name='tauvo',        &
         cmor_units='N m-2', cmor_long_name='Surface Downward Y Stress',            &
         cmor_standard_name='surface_downward_y_stress')

  handles%id_ustar = register_diag_field('ocean_model', 'ustar', diag%axesT1, Time, &
      'Surface friction velocity = [(gustiness + tau_magnitude)/rho0]^(1/2)', 'meter second-1')

  if (present(use_berg_fluxes)) then
    if (use_berg_fluxes) then
      handles%id_ustar_berg = register_diag_field('ocean_model', 'ustar_berg', diag%axesT1, Time, &
          'Friction velocity below iceberg ', 'meter second-1')

      handles%id_area_berg = register_diag_field('ocean_model', 'area_berg', diag%axesT1, Time, &
          'Area of grid cell covered by iceberg ', 'm2/m2')

      handles%id_mass_berg = register_diag_field('ocean_model', 'mass_berg', diag%axesT1, Time, &
          'Mass of icebergs ', 'kg/m2')

      handles%id_ustar_ice_cover = register_diag_field('ocean_model', 'ustar_ice_cover', diag%axesT1, Time, &
          'Friction velocity below iceberg and ice shelf together', 'meter second-1')

      handles%id_frac_ice_cover = register_diag_field('ocean_model', 'frac_ice_cover', diag%axesT1, Time, &
          'Area of grid cell below iceberg and ice shelf together ', 'm2/m2')
    endif
  endif

  handles%id_psurf = register_diag_field('ocean_model', 'p_surf', diag%axesT1, Time,           &
        'Pressure at ice-ocean or atmosphere-ocean interface', 'Pascal', cmor_field_name='pso',&
        cmor_long_name='Sea Water Pressure at Sea Water Surface', cmor_units='Pa',             &
        cmor_standard_name='sea_water_pressure_at_sea_water_surface')

  handles%id_TKE_tidal = register_diag_field('ocean_model', 'TKE_tidal', diag%axesT1, Time, &
        'Tidal source of BBL mixing', 'Watt/m^2')

  if (.not. use_temperature) then
    handles%id_buoy = register_diag_field('ocean_model', 'buoy', diag%axesT1, Time, &
          'Buoyancy forcing', 'meter^2/second^3')
    return
  endif


  !===============================================================
  ! surface mass flux maps

  handles%id_prcme = register_diag_field('ocean_model', 'PRCmE', diag%axesT1, Time,                  &
        'Net surface water flux (precip+melt+lrunoff+ice calving-evap)', 'kilogram meter-2 second-1',&
        standard_name='water_flux_into_sea_water', cmor_field_name='wfo', cmor_units='kg m-2 s-1',   &
        cmor_standard_name='water_flux_into_sea_water',cmor_long_name='Water Flux Into Sea Water')

  handles%id_evap = register_diag_field('ocean_model', 'evap', diag%axesT1, Time,                         &
       'Evaporation/condensation at ocean surface (evaporation is negative)', 'kilogram meter-2 second-1',&
       standard_name='water_evaporation_flux', cmor_field_name='evs', cmor_units='kg m-2 s-1',            &
       cmor_standard_name='water_evaporation_flux',                                                       &
       cmor_long_name='Water Evaporation Flux Where Ice Free Ocean over Sea')

  ! smg: seaice_melt field requires updates to the sea ice model
  !handles%id_seaice_melt = register_diag_field('ocean_model', 'seaice_melt',       &
  !   diag%axesT1, Time, 'water flux to ocean from sea ice melt(> 0) or form(< 0)', &
  !   'kilogram/(meter^2 * second)',                                                &
  !    standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics',     &
  !    cmor_field_name='fsitherm', cmor_units='kg m-2 s-1',                         &
  !    cmor_standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics',&
  !    cmor_long_name='water flux to ocean from sea ice melt(> 0) or form(< 0)')

  handles%id_precip = register_diag_field('ocean_model', 'precip', diag%axesT1, Time, &
        'Liquid + frozen precipitation into ocean', 'kilogram/(meter^2 * second)')

  handles%id_fprec = register_diag_field('ocean_model', 'fprec', diag%axesT1, Time,     &
        'Frozen precipitation into ocean', 'kilogram meter-2 second-1',                 &
        standard_name='snowfall_flux', cmor_field_name='prsn', cmor_units='kg m-2 s-1', &
        cmor_standard_name='snowfall_flux', cmor_long_name='Snowfall Flux where Ice Free Ocean over Sea')

  handles%id_lprec = register_diag_field('ocean_model', 'lprec', diag%axesT1, Time,       &
        'Liquid precipitation into ocean', 'kilogram/(meter^2 * second)',                 &
        standard_name='rainfall_flux',                                                    &
        cmor_field_name='prlq', cmor_units='kg m-2 s-1', cmor_standard_name='rainfall_flux',&
        cmor_long_name='Rainfall Flux where Ice Free Ocean over Sea')

  handles%id_vprec = register_diag_field('ocean_model', 'vprec', diag%axesT1, Time, &
        'Virtual liquid precip into ocean due to SSS restoring', 'kilogram/(meter^2 second)')

  handles%id_frunoff = register_diag_field('ocean_model', 'frunoff', diag%axesT1, Time,    &
        'Frozen runoff (calving) and iceberg melt into ocean', 'kilogram/(meter^2 second)',&
        standard_name='water_flux_into_sea_water_from_icebergs',                           &
        cmor_field_name='ficeberg', cmor_units='kg m-2 s-1',                               &
        cmor_standard_name='water_flux_into_sea_water_from_icebergs',                      &
        cmor_long_name='Water Flux into Seawater from Icebergs')

  handles%id_lrunoff = register_diag_field('ocean_model', 'lrunoff', diag%axesT1, Time, &
        'Liquid runoff (rivers) into ocean', 'kilogram meter-2 second-1',                     &
        standard_name='water_flux_into_sea_water_from_rivers', cmor_field_name='friver',      &
        cmor_units='kg m-2 s-1', cmor_standard_name='water_flux_into_sea_water_from_rivers',  &
        cmor_long_name='Water Flux into Sea Water From Rivers')

  handles%id_net_massout = register_diag_field('ocean_model', 'net_massout', diag%axesT1, Time, &
        'Net mass leaving the ocean due to evaporation, seaice formation', 'kilogram meter-2 second-1')

  handles%id_net_massin  = register_diag_field('ocean_model', 'net_massin', diag%axesT1, Time, &
        'Net mass entering ocean due to precip, runoff, ice melt', 'kilogram meter-2 second-1')

  handles%id_massout_flux = register_diag_field('ocean_model', 'massout_flux', diag%axesT1, Time, &
        'Net mass flux of freshwater out of the ocean (used in the boundary flux calculation)', &
         'kilogram meter-2')

  handles%id_massin_flux  = register_diag_field('ocean_model', 'massin_flux', diag%axesT1, Time, &
        'Net mass mass flux of freshwater into the ocean (used in boundary flux calculation)', 'kilogram meter-2')
  !=========================================================================
  ! area integrated surface mass transport

  handles%id_total_prcme = register_scalar_field('ocean_model', 'total_PRCmE', Time, diag,         &
      long_name='Area integrated net surface water flux (precip+melt+liq runoff+ice calving-evap)',&
      units='kg/s', standard_name='water_flux_into_sea_water_area_integrated',                     &
      cmor_field_name='total_wfo', cmor_units='kg s-1',                                            &
      cmor_standard_name='water_flux_into_sea_water_area_integrated',                              &
      cmor_long_name='Water Transport Into Sea Water Area Integrated')

  handles%id_total_evap = register_scalar_field('ocean_model', 'total_evap', Time, diag,&
      long_name='Area integrated evap/condense at ocean surface',                       &
      units='kg/s', standard_name='water_evaporation_flux_area_integrated',             &
      cmor_field_name='total_evs', cmor_units='kg s-1',                                 &
      cmor_standard_name='water_evaporation_flux_area_integrated',                      &
      cmor_long_name='Evaporation Where Ice Free Ocean over Sea Area Integrated')

  ! seaice_melt field requires updates to the sea ice model
  !handles%id_total_seaice_melt = register_scalar_field('ocean_model', 'total_seaice_melt', Time, diag, &
  !    long_name='Area integrated sea ice melt (>0) or form (<0)', units='kg/s',                        &
  !    standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics_area_integrated',         &
  !    cmor_field_name='total_fsitherm', cmor_units='kg s-1',                                           &
  !    cmor_standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics_area_integrated',    &
  !    cmor_long_name='Water Melt/Form from Sea Ice Area Integrated')

  handles%id_total_precip = register_scalar_field('ocean_model', 'total_precip', Time, diag, &
      long_name='Area integrated liquid+frozen precip into ocean', units='kg/s')

  handles%id_total_fprec = register_scalar_field('ocean_model', 'total_fprec', Time, diag,&
      long_name='Area integrated frozen precip into ocean', units='kg/s',                 &
      standard_name='snowfall_flux_area_integrated',                                      &
      cmor_field_name='total_prsn', cmor_units='kg s-1',                                  &
      cmor_standard_name='snowfall_flux_area_integrated',                                 &
      cmor_long_name='Snowfall Flux where Ice Free Ocean over Sea Area Integrated')

  handles%id_total_lprec = register_scalar_field('ocean_model', 'total_lprec', Time, diag,&
      long_name='Area integrated liquid precip into ocean', units='kg/s',                 &
      standard_name='rainfall_flux_area_integrated',                                      &
      cmor_field_name='total_pr', cmor_units='kg s-1',                                    &
      cmor_standard_name='rainfall_flux_area_integrated',                                 &
      cmor_long_name='Rainfall Flux where Ice Free Ocean over Sea Area Integrated')

  handles%id_total_vprec = register_scalar_field('ocean_model', 'total_vprec', Time, diag, &
      long_name='Area integrated virtual liquid precip due to SSS restoring', units='kg/s')

  handles%id_total_frunoff = register_scalar_field('ocean_model', 'total_frunoff', Time, diag,    &
      long_name='Area integrated frozen runoff (calving) & iceberg melt into ocean', units='kg/s',&
      cmor_field_name='total_ficeberg', cmor_units='kg s-1',                                      &
      cmor_standard_name='water_flux_into_sea_water_from_icebergs_area_integrated',               &
      cmor_long_name='Water Flux into Seawater from Icebergs Area Integrated')

  handles%id_total_lrunoff = register_scalar_field('ocean_model', 'total_lrunoff', Time, diag,&
      long_name='Area integrated liquid runoff into ocean', units='kg/s',                     &
      cmor_field_name='total_friver', cmor_units='kg s-1',                                    &
      cmor_standard_name='water_flux_into_sea_water_from_rivers_area_integrated',             &
      cmor_long_name='Water Flux into Sea Water From Rivers Area Integrated')

  handles%id_total_net_massout = register_scalar_field('ocean_model', 'total_net_massout', Time, diag, &
      long_name='Area integrated mass leaving ocean due to evap and seaice form', units='kg/s')

  handles%id_total_net_massin = register_scalar_field('ocean_model', 'total_net_massin', Time, diag, &
      long_name='Area integrated mass entering ocean due to predip, runoff, ice melt', units='kg/s')

  !=========================================================================
  ! area averaged surface mass transport

  handles%id_prcme_ga = register_scalar_field('ocean_model', 'PRCmE_ga', Time, diag,             &
      long_name='Area averaged net surface water flux (precip+melt+liq runoff+ice calving-evap)',&
      units='kg m-2 s-1', standard_name='water_flux_into_sea_water_area_averaged',               &
      cmor_field_name='ave_wfo', cmor_units='kg m-2 s-1',                                        &
      cmor_standard_name='rainfall_flux_area_averaged',                                          &
      cmor_long_name='Water Transport Into Sea Water Area Averaged')

  handles%id_evap_ga = register_scalar_field('ocean_model', 'evap_ga', Time, diag,&
      long_name='Area averaged evap/condense at ocean surface',                   &
      units='kg m-2 s-1', standard_name='water_evaporation_flux_area_averaged',   &
      cmor_field_name='ave_evs', cmor_units='kg m-2 s-1',                         &
      cmor_standard_name='water_evaporation_flux_area_averaged',                  &
      cmor_long_name='Evaporation Where Ice Free Ocean over Sea Area Averaged')

 handles%id_lprec_ga = register_scalar_field('ocean_model', 'lprec_ga', Time, diag,&
      long_name='Area integrated liquid precip into ocean', units='kg m-2 s-1',    &
      standard_name='rainfall_flux_area_averaged',                                 &
      cmor_field_name='ave_pr', cmor_units='kg m-2 s-1',                           &
      cmor_standard_name='rainfall_flux_area_averaged',                            &
      cmor_long_name='Rainfall Flux where Ice Free Ocean over Sea Area Averaged')

 handles%id_fprec_ga = register_scalar_field('ocean_model', 'fprec_ga', Time, diag,&
      long_name='Area integrated frozen precip into ocean', units='kg m-2 s-1',    &
      standard_name='snowfall_flux_area_averaged',                                 &
      cmor_field_name='ave_prsn', cmor_units='kg m-2 s-1',                         &
      cmor_standard_name='snowfall_flux_area_averaged',                            &
      cmor_long_name='Snowfall Flux where Ice Free Ocean over Sea Area Averaged')

  handles%id_precip_ga = register_scalar_field('ocean_model', 'precip_ga', Time, diag, &
      long_name='Area averaged liquid+frozen precip into ocean', units='kg m-2 s-1')

  handles%id_vprec_ga = register_scalar_field('ocean_model', 'vrec_ga', Time, diag, &
      long_name='Area averaged virtual liquid precip due to SSS restoring', units='kg m-2 s-1')

  !===============================================================
  ! surface heat flux maps

  handles%id_heat_content_frunoff = register_diag_field('ocean_model', 'heat_content_frunoff',        &
        diag%axesT1, Time, 'Heat content (relative to 0C) of solid runoff into ocean', 'Watt meter-2',&
        standard_name='temperature_flux_due_to_solid_runoff_expressed_as_heat_flux_into_sea_water')

  handles%id_heat_content_lrunoff = register_diag_field('ocean_model', 'heat_content_lrunoff',         &
        diag%axesT1, Time, 'Heat content (relative to 0C) of liquid runoff into ocean', 'Watt meter-2',&
        standard_name='temperature_flux_due_to_runoff_expressed_as_heat_flux_into_sea_water')

  handles%id_hfrunoffds = register_diag_field('ocean_model', 'hfrunoffds',                            &
        diag%axesT1, Time, 'Heat content (relative to 0C) of liquid+solid runoff into ocean', 'W m-2',&
        standard_name='temperature_flux_due_to_runoff_expressed_as_heat_flux_into_sea_water')

  handles%id_heat_content_lprec = register_diag_field('ocean_model', 'heat_content_lprec',             &
        diag%axesT1,Time,'Heat content (relative to 0degC) of liquid precip entering ocean',           &
        'W/m^2')

  handles%id_heat_content_fprec = register_diag_field('ocean_model', 'heat_content_fprec',&
        diag%axesT1,Time,'Heat content (relative to 0degC) of frozen prec entering ocean',&
        'W/m^2')

  handles%id_heat_content_vprec = register_diag_field('ocean_model', 'heat_content_vprec',   &
        diag%axesT1,Time,'Heat content (relative to 0degC) of virtual precip entering ocean',&
        'Watt/m^2')

  handles%id_heat_content_cond = register_diag_field('ocean_model', 'heat_content_cond',   &
        diag%axesT1,Time,'Heat content (relative to 0degC) of water condensing into ocean',&
        'Watt/m^2')

  handles%id_hfrainds = register_diag_field('ocean_model', 'hfrainds',                                 &
        diag%axesT1,Time,'Heat content (relative to 0degC) of liquid+frozen precip entering ocean',    &
        'W/m^2',standard_name='temperature_flux_due_to_rainfall_expressed_as_heat_flux_into_sea_water',&
        cmor_long_name='Heat Content (relative to 0degC) of Liquid + Frozen Precipitation')

  handles%id_heat_content_surfwater = register_diag_field('ocean_model', 'heat_content_surfwater',&
         diag%axesT1, Time,                                                                       &
        'Heat content (relative to 0degC) of net water crossing ocean surface (frozen+liquid)',   &
        'Watt/m^2')

  handles%id_heat_content_massout = register_diag_field('ocean_model', 'heat_content_massout',                      &
         diag%axesT1, Time,'Heat content (relative to 0degC) of net mass leaving ocean ocean via evap and ice form',&
        'Watt/m^2',                                                                                                 &
        cmor_field_name='hfevapds', cmor_units='W m-2',                                                             &
        cmor_standard_name='temperature_flux_due_to_evaporation_expressed_as_heat_flux_out_of_sea_water',           &
        cmor_long_name='Heat Content (relative to 0degC) of Water Leaving Ocean via Evaporation and Ice Formation')

  handles%id_heat_content_massin = register_diag_field('ocean_model', 'heat_content_massin',   &
         diag%axesT1, Time,'Heat content (relative to 0degC) of net mass entering ocean ocean',&
        'Watt/m^2')

  handles%id_net_heat_coupler = register_diag_field('ocean_model', 'net_heat_coupler',          &
        diag%axesT1,Time,'Surface ocean heat flux from SW+LW+latent+sensible (via the coupler)',&
        'Watt/m^2')

  handles%id_net_heat_surface = register_diag_field('ocean_model', 'net_heat_surface',diag%axesT1,  &
        Time,'Surface ocean heat flux from SW+LW+lat+sens+mass transfer+frazil+restore or flux adjustments', 'Watt/m^2',&
        standard_name='surface_downward_heat_flux_in_sea_water', cmor_field_name='hfds',            &
        cmor_units='W m-2', cmor_standard_name='surface_downward_heat_flux_in_sea_water',           &
        cmor_long_name='Surface ocean heat flux from SW+LW+latent+sensible+masstransfer+frazil')

  handles%id_sw = register_diag_field('ocean_model', 'SW', diag%axesT1, Time,  &
        'Shortwave radiation flux into ocean', 'Watt meter-2',                 &
        standard_name='net_downward_shortwave_flux_at_sea_water_surface',      &
        cmor_field_name='rsntds', cmor_units='W m-2',                          &
        cmor_standard_name='net_downward_shortwave_flux_at_sea_water_surface', &
        cmor_long_name='Net Downward Shortwave Radiation at Sea Water Surface')
  handles%id_sw_vis = register_diag_field('ocean_model', 'sw_vis', diag%axesT1, Time,     &
        'Shortwave radiation direct and diffuse flux into the ocean in the visible band', &
        'Watt/m^2')
  handles%id_sw_nir = register_diag_field('ocean_model', 'sw_nir', diag%axesT1, Time,     &
        'Shortwave radiation direct and diffuse flux into the ocean in the near-infrared band', &
        'Watt/m^2')

  handles%id_LwLatSens = register_diag_field('ocean_model', 'LwLatSens', diag%axesT1, Time, &
        'Combined longwave, latent, and sensible heating at ocean surface', 'Watt/m^2')

  handles%id_lw = register_diag_field('ocean_model', 'LW', diag%axesT1, Time, &
        'Longwave radiation flux into ocean', 'Watt meter-2',                 &
        standard_name='surface_net_downward_longwave_flux',                   &
        cmor_field_name='rlntds', cmor_units='W m-2',                         &
        cmor_standard_name='surface_net_downward_longwave_flux',              &
        cmor_long_name='Surface Net Downward Longwave Radiation')

  handles%id_lat = register_diag_field('ocean_model', 'latent', diag%axesT1, Time,                    &
        'Latent heat flux into ocean due to fusion and evaporation (negative means ocean heat loss)', &
        'Watt meter-2', cmor_field_name='hflso', cmor_units='W m-2',                                  &
        cmor_standard_name='surface_downward_latent_heat_flux',                                       &
        cmor_long_name='Surface Downward Latent Heat Flux due to Evap + Melt Snow/Ice')

  handles%id_lat_evap = register_diag_field('ocean_model', 'latent_evap', diag%axesT1, Time, &
        'Latent heat flux into ocean due to evaporation/condensation', 'Watt/m^2')

  handles%id_lat_fprec = register_diag_field('ocean_model', 'latent_fprec_diag', diag%axesT1, Time,&
        'Latent heat flux into ocean due to melting of frozen precipitation', 'Watt meter-2',      &
        cmor_field_name='hfsnthermds', cmor_units='W m-2',                                         &
        cmor_standard_name='heat_flux_into_sea_water_due_to_snow_thermodynamics',                  &
        cmor_long_name='Latent Heat to Melt Frozen Precipitation')

  handles%id_lat_frunoff = register_diag_field('ocean_model', 'latent_frunoff', diag%axesT1, Time, &
        'Latent heat flux into ocean due to melting of icebergs', 'Watt/m^2',                      &
        cmor_field_name='hfibthermds', cmor_units='W m-2',                                         &
        cmor_standard_name='heat_flux_into_sea_water_due_to_iceberg_thermodynamics',               &
        cmor_long_name='Latent Heat to Melt Frozen Runoff/Iceberg')

  handles%id_sens = register_diag_field('ocean_model', 'sensible', diag%axesT1, Time,&
        'Sensible heat flux into ocean', 'Watt meter-2',                             &
        standard_name='surface_downward_sensible_heat_flux',                         &
        cmor_field_name='hfsso', cmor_units='W m-2',                                 &
        cmor_standard_name='surface_downward_sensible_heat_flux',                    &
        cmor_long_name='Surface Downward Sensible Heat Flux')

  handles%id_heat_added = register_diag_field('ocean_model', 'heat_added', diag%axesT1, Time, &
        'Flux Adjustment or restoring surface heat flux into ocean', 'Watt/m^2')


  !===============================================================
  ! area integrated surface heat fluxes

  handles%id_total_heat_content_frunoff = register_scalar_field('ocean_model',                     &
      'total_heat_content_frunoff', Time, diag,                                                    &
      long_name='Area integrated heat content (relative to 0C) of solid runoff',                   &
      units='Watt', cmor_field_name='total_hfsolidrunoffds', cmor_units='W',                       &
      cmor_standard_name=                                                                          &
      'temperature_flux_due_to_solid_runoff_expressed_as_heat_flux_into_sea_water_area_integrated',&
      cmor_long_name=                                                                              &
      'Temperature Flux due to Solid Runoff Expressed as Heat Flux into Sea Water Area Integrated')

  handles%id_total_heat_content_lrunoff = register_scalar_field('ocean_model',               &
      'total_heat_content_lrunoff', Time, diag,                                              &
      long_name='Area integrated heat content (relative to 0C) of liquid runoff',            &
      units='Watt', cmor_field_name='total_hfrunoffds', cmor_units='W',                      &
      cmor_standard_name=                                                                    &
      'temperature_flux_due_to_runoff_expressed_as_heat_flux_into_sea_water_area_integrated',&
      cmor_long_name=                                                                        &
      'Temperature Flux due to Runoff Expressed as Heat Flux into Sea Water Area Integrated')

  handles%id_total_heat_content_lprec = register_scalar_field('ocean_model',                   &
      'total_heat_content_lprec', Time, diag,                                                  &
      long_name='Area integrated heat content (relative to 0C) of liquid precip',              &
      units='Watt', cmor_field_name='total_hfrainds', cmor_units='W',                          &
      cmor_standard_name=                                                                      &
      'temperature_flux_due_to_rainfall_expressed_as_heat_flux_into_sea_water_area_integrated',&
      cmor_long_name=                                                                          &
      'Temperature Flux due to Rainfall Expressed as Heat Flux into Sea Water Area Integrated')

  handles%id_total_heat_content_fprec = register_scalar_field('ocean_model',     &
      'total_heat_content_fprec', Time, diag,                                    &
      long_name='Area integrated heat content (relative to 0C) of frozen precip',&
      units='Watt')

  handles%id_total_heat_content_vprec = register_scalar_field('ocean_model',      &
      'total_heat_content_vprec', Time, diag,                                     &
      long_name='Area integrated heat content (relative to 0C) of virtual precip',&
      units='Watt')

  handles%id_total_heat_content_cond = register_scalar_field('ocean_model',   &
      'total_heat_content_cond', Time, diag,                                  &
      long_name='Area integrated heat content (relative to 0C) of condensate',&
      units='Watt')

  handles%id_total_heat_content_surfwater = register_scalar_field('ocean_model',          &
      'total_heat_content_surfwater', Time, diag,                                         &
      long_name='Area integrated heat content (relative to 0C) of water crossing surface',&
      units='Watt')

  handles%id_total_heat_content_massout = register_scalar_field('ocean_model',                      &
      'total_heat_content_massout', Time, diag,                                                     &
      long_name='Area integrated heat content (relative to 0C) of water leaving ocean',             &
      units='Watt',                                                                                 &
      cmor_field_name='total_hfevapds', cmor_units='W',                                             &
      cmor_standard_name=                                                                           &
      'temperature_flux_due_to_evaporation_expressed_as_heat_flux_out_of_sea_water_area_integrated',&
      cmor_long_name='Heat Flux Out of Sea Water due to Evaporating Water Area Integrated')

  handles%id_total_heat_content_massin = register_scalar_field('ocean_model',           &
      'total_heat_content_massin', Time, diag,                                          &
      long_name='Area integrated heat content (relative to 0C) of water entering ocean',&
      units='Watt')

  handles%id_total_net_heat_coupler = register_scalar_field('ocean_model',                       &
      'total_net_heat_coupler', Time, diag,                                                      &
      long_name='Area integrated surface heat flux from SW+LW+latent+sensible (via the coupler)',&
      units='Watt')

  handles%id_total_net_heat_surface = register_scalar_field('ocean_model',                      &
      'total_net_heat_surface', Time, diag,                                                     &
      long_name='Area integrated surface heat flux from SW+LW+lat+sens+mass+frazil+restore or flux adjustments',    &
      units='Watt',                                                                             &
      cmor_field_name='total_hfds', cmor_units='W',                                             &
      cmor_standard_name='surface_downward_heat_flux_in_sea_water_area_integrated',             &
      cmor_long_name=                                                                           &
      'Surface Ocean Heat Flux from SW+LW+latent+sensible+mass transfer+frazil Area Integrated')

  handles%id_total_sw = register_scalar_field('ocean_model',                                &
      'total_sw', Time, diag,                                                               &
      long_name='Area integrated net downward shortwave at sea water surface',              &
      units='Watt',                                                                         &
      cmor_field_name='total_rsntds', cmor_units='W',                                       &
      cmor_standard_name='net_downward_shortwave_flux_at_sea_water_surface_area_integrated',&
      cmor_long_name=                                                                       &
      'Net Downward Shortwave Radiation at Sea Water Surface Area Integrated')

  handles%id_total_LwLatSens = register_scalar_field('ocean_model',&
      'total_LwLatSens', Time, diag,                               &
      long_name='Area integrated longwave+latent+sensible heating',&
      units='Watt')

  handles%id_total_lw = register_scalar_field('ocean_model',                  &
      'total_lw', Time, diag,                                                 &
      long_name='Area integrated net downward longwave at sea water surface', &
      units='Watt',                                                           &
      cmor_field_name='total_rlntds', cmor_units='W',                         &
      cmor_standard_name='surface_net_downward_longwave_flux_area_integrated',&
      cmor_long_name=                                                         &
      'Surface Net Downward Longwave Radiation Area Integrated')

  handles%id_total_lat = register_scalar_field('ocean_model',                &
      'total_lat', Time, diag,                                               &
      long_name='Area integrated surface downward latent heat flux',         &
      units='Watt',                                                          &
      cmor_field_name='total_hflso', cmor_units='W',                         &
      cmor_standard_name='surface_downward_latent_heat_flux_area_integrated',&
      cmor_long_name=                                                        &
      'Surface Downward Latent Heat Flux Area Integrated')

  handles%id_total_lat_evap = register_scalar_field('ocean_model',      &
      'total_lat_evap', Time, diag,                                     &
      long_name='Area integrated latent heat flux due to evap/condense',&
      units='Watt')

  handles%id_total_lat_fprec = register_scalar_field('ocean_model',                            &
      'total_lat_fprec', Time, diag,                                                           &
      long_name='Area integrated latent heat flux due to melting frozen precip',               &
      units='Watt',                                                                            &
      cmor_field_name='total_hfsnthermds', cmor_units='W',                                     &
      cmor_standard_name='heat_flux_into_sea_water_due_to_snow_thermodynamics_area_integrated',&
      cmor_long_name=                                                                          &
      'Latent Heat to Melt Frozen Precipitation Area Integrated')

  handles%id_total_lat_frunoff = register_scalar_field('ocean_model',                             &
      'total_lat_frunoff', Time, diag,                                                            &
      long_name='Area integrated latent heat flux due to melting icebergs',                       &
      units='Watt',                                                                               &
      cmor_field_name='total_hfibthermds', cmor_units='W',                                        &
      cmor_standard_name='heat_flux_into_sea_water_due_to_iceberg_thermodynamics_area_integrated',&
      cmor_long_name=                                                                             &
      'Heat Flux into Sea Water due to Iceberg Thermodynamics Area Integrated')

  handles%id_total_sens = register_scalar_field('ocean_model',                 &
      'total_sens', Time, diag,                                                &
      long_name='Area integrated downward sensible heat flux',                 &
      units='Watt',                                                            &
      cmor_field_name='total_hfsso', cmor_units='W',                           &
      cmor_standard_name='surface_downward_sensible_heat_flux_area_integrated',&
      cmor_long_name=                                                          &
      'Surface Downward Sensible Heat Flux Area Integrated')

  handles%id_total_heat_added = register_scalar_field('ocean_model',&
      'total_heat_adjustment', Time, diag,                               &
      long_name='Area integrated surface heat flux from restoring and/or flux adjustment',   &
      units='Watt')


  !===============================================================
  ! area averaged surface heat fluxes

  handles%id_net_heat_coupler_ga = register_scalar_field('ocean_model',                       &
      'net_heat_coupler_ga', Time, diag,                                                      &
      long_name='Area averaged surface heat flux from SW+LW+latent+sensible (via the coupler)',&
      units='W m-2')

  handles%id_net_heat_surface_ga = register_scalar_field('ocean_model',                       &
      'net_heat_surface_ga', Time, diag,                                                      &
      long_name='Area averaged surface heat flux from SW+LW+lat+sens+mass+frazil+restore or flux adjustments',    &
      units='W m-2',                                                                          &
      cmor_field_name='ave_hfds', cmor_units='W m-2',                                         &
      cmor_standard_name='surface_downward_heat_flux_in_sea_water_area_averaged',             &
      cmor_long_name=                                                                         &
      'Surface Ocean Heat Flux from SW+LW+latent+sensible+mass transfer+frazil Area Averaged')

  handles%id_sw_ga = register_scalar_field('ocean_model',                                 &
      'sw_ga', Time, diag,                                                                &
      long_name='Area averaged net downward shortwave at sea water surface',              &
      units='W m-2',                                                                      &
      cmor_field_name='ave_rsntds', cmor_units='W m-2',                                   &
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
      cmor_field_name='ave_rlntds', cmor_units='W m-2',                     &
      cmor_standard_name='surface_net_downward_longwave_flux_area_averaged',&
      cmor_long_name=                                                       &
      'Surface Net Downward Longwave Radiation Area Averaged')

  handles%id_lat_ga = register_scalar_field('ocean_model',                 &
      'lat_ga', Time, diag,                                                &
      long_name='Area averaged surface downward latent heat flux',         &
      units='W m-2',                                                       &
      cmor_field_name='ave_hflso', cmor_units='W m-2',                     &
      cmor_standard_name='surface_downward_latent_heat_flux_area_averaged',&
      cmor_long_name=                                                      &
      'Surface Downward Latent Heat Flux Area Averaged')

  handles%id_sens_ga = register_scalar_field('ocean_model',                  &
      'sens_ga', Time, diag,                                                 &
      long_name='Area averaged downward sensible heat flux',                 &
      units='W m-2',                                                         &
      cmor_field_name='ave_hfsso', cmor_units='W m-2',                       &
      cmor_standard_name='surface_downward_sensible_heat_flux_area_averaged',&
      cmor_long_name=                                                        &
      'Surface Downward Sensible Heat Flux Area Averaged')


  !===============================================================
  ! maps of surface salt fluxes, virtual precip fluxes, and adjustments

  handles%id_saltflux = register_diag_field('ocean_model', 'salt_flux', diag%axesT1, Time,&
        'Net salt flux into ocean at surface (restoring + sea-ice)',                      &
        'kilogram meter-2 second-1',cmor_field_name='sfdsi', cmor_units='kg m-2 s-1',     &
        cmor_standard_name='downward_sea_ice_basal_salt_flux',                            &
        cmor_long_name='Downward Sea Ice Basal Salt Flux')

  handles%id_saltFluxIn = register_diag_field('ocean_model', 'salt_flux_in', diag%axesT1, Time, &
        'Salt flux into ocean at surface from coupler', 'kilogram/(meter^2 * second)')

  handles%id_saltFluxAdded = register_diag_field('ocean_model', 'salt_flux_added', &
        diag%axesT1,Time,'Salt flux into ocean at surface due to restoring or flux adjustment',           &
        'kilogram/(meter^2 * second)')

  handles%id_saltFluxGlobalAdj = register_scalar_field('ocean_model',              &
        'salt_flux_global_restoring_adjustment', Time, diag,                       &
        'Adjustment needed to balance net global salt flux into ocean at surface', &
        'kilogram/(meter^2 * second)')

  handles%id_vPrecGlobalAdj = register_scalar_field('ocean_model',  &
        'vprec_global_adjustment', Time, diag,                      &
        'Adjustment needed to adjust net vprec into ocean to zero', &
        'kilogram/(meter^2 * second)')

  handles%id_netFWGlobalAdj = register_scalar_field('ocean_model',       &
        'net_fresh_water_global_adjustment', Time, diag,                 &
        'Adjustment needed to adjust net fresh water into ocean to zero',&
        'kilogram/(meter^2 * second)')

  handles%id_saltFluxGlobalScl = register_scalar_field('ocean_model',            &
        'salt_flux_global_restoring_scaling', Time, diag,                        &
        'Scaling applied to balance net global salt flux into ocean at surface', &
        '(nondim)')

  handles%id_vPrecGlobalScl = register_scalar_field('ocean_model',&
        'vprec_global_scaling', Time, diag,                       &
        'Scaling applied to adjust net vprec into ocean to zero', &
        '(nondim)')

  handles%id_netFWGlobalScl = register_scalar_field('ocean_model',      &
        'net_fresh_water_global_scaling', Time, diag,                   &
        'Scaling applied to adjust net fresh water into ocean to zero', &
        '(nondim)')

  !===============================================================
  ! area integrals of surface salt fluxes

  handles%id_total_saltflux = register_scalar_field('ocean_model',          &
      'total_salt_flux', Time, diag,                                        &
      long_name='Area integrated surface salt flux', units='kg',            &
      cmor_field_name='total_sfdsi',                                        &
      cmor_units='kg s-1',                                                  &
      cmor_standard_name='downward_sea_ice_basal_salt_flux_area_integrated',&
      cmor_long_name='Downward Sea Ice Basal Salt Flux Area Integrated')

  handles%id_total_saltFluxIn = register_scalar_field('ocean_model', 'total_salt_Flux_In', &
      Time, diag, long_name='Area integrated surface salt flux at surface from coupler', units='kg')

  handles%id_total_saltFluxAdded = register_scalar_field('ocean_model', 'total_salt_Flux_Added', &
      Time, diag, long_name='Area integrated surface salt flux due to restoring or flux adjustment', units='kg')


end subroutine register_forcing_type_diags

!> Accumulate the forcing over time steps
subroutine forcing_accumulate(flux_tmp, fluxes, dt, G, wt2)
  type(forcing),         intent(in)    :: flux_tmp
  type(forcing),         intent(inout) :: fluxes
  real,                  intent(in)    :: dt
  type(ocean_grid_type), intent(inout) :: G
  real,                  intent(out)   :: wt2

  ! This subroutine copies mechancal forcing from flux_tmp to fluxes and
  ! stores the time-weighted averages of the various buoyancy fluxes in fluxes,
  ! and increments the amount of time over which the buoyancy forcing should be
  ! applied.

  real :: wt1
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, i0, j0
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isr, ier, jsr, jer
  is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
  Isq  = G%IscB  ; Ieq  = G%IecB   ; Jsq  = G%JscB  ; Jeq  = G%JecB
  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB


  if (fluxes%dt_buoy_accum < 0) call MOM_error(FATAL, "forcing_accumulate: "//&
     "fluxes must be initialzed before it can be augmented.")

  ! wt1 is the relative weight of the previous fluxes.
  wt1 = fluxes%dt_buoy_accum / (fluxes%dt_buoy_accum + dt)
  wt2 = 1.0 - wt1 ! = dt / (fluxes%dt_buoy_accum + dt)
  fluxes%dt_buoy_accum = fluxes%dt_buoy_accum + dt

  ! Copy over the pressure and momentum flux fields.
  do j=js,je ; do i=is,ie
    fluxes%p_surf(i,j) = flux_tmp%p_surf(i,j)
    fluxes%p_surf_full(i,j) = flux_tmp%p_surf_full(i,j)
  enddo ; enddo
  do j=js,je ; do I=Isq,Ieq
    fluxes%taux(I,j) = flux_tmp%taux(I,j)
  enddo ; enddo
  do J=Jsq,Jeq ; do i=is,ie
    fluxes%tauy(i,J) = flux_tmp%tauy(i,J)
  enddo ; enddo

  ! Average the water, heat, and salt fluxes, and ustar.
  do j=js,je ; do i=is,ie
    fluxes%ustar(i,j) = wt1*fluxes%ustar(i,j) + wt2*flux_tmp%ustar(i,j)

    fluxes%evap(i,j) = wt1*fluxes%evap(i,j) + wt2*flux_tmp%evap(i,j)
    fluxes%lprec(i,j) = wt1*fluxes%lprec(i,j) + wt2*flux_tmp%lprec(i,j)
    fluxes%fprec(i,j) = wt1*fluxes%fprec(i,j) + wt2*flux_tmp%fprec(i,j)
    fluxes%vprec(i,j) = wt1*fluxes%vprec(i,j) + wt2*flux_tmp%vprec(i,j)
    fluxes%lrunoff(i,j) = wt1*fluxes%lrunoff(i,j) + wt2*flux_tmp%lrunoff(i,j)
    fluxes%frunoff(i,j) = wt1*fluxes%frunoff(i,j) + wt2*flux_tmp%frunoff(i,j)
 ! ### ADD LATER fluxes%seaice_melt(i,j) = wt1*fluxes%seaice_melt(i,j) + wt2*flux_tmp%seaice_melt(i,j)

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
  if (associated(fluxes%frac_shelf_u) .and. associated(flux_tmp%frac_shelf_u)) then
    do I=IsdB,IedB ; do j=jsd,jed
      fluxes%frac_shelf_u(I,j)  = flux_tmp%frac_shelf_u(I,j)
    enddo ; enddo
  endif
  if (associated(fluxes%frac_shelf_v) .and. associated(flux_tmp%frac_shelf_v)) then
    do i=isd,ied ; do J=JsdB,JedB
      fluxes%frac_shelf_v(i,J)  = flux_tmp%frac_shelf_v(i,J)
    enddo ; enddo
  endif
  if (associated(fluxes%rigidity_ice_u) .and. associated(flux_tmp%rigidity_ice_u)) then
    do I=isd,ied-1 ; do j=jsd,jed
      fluxes%rigidity_ice_u(I,j)  = flux_tmp%rigidity_ice_u(I,j)
    enddo ; enddo
  endif
  if (associated(fluxes%rigidity_ice_v) .and. associated(flux_tmp%rigidity_ice_v)) then
    do i=isd,ied ; do J=jsd,jed-1
      fluxes%rigidity_ice_v(i,J)  = flux_tmp%rigidity_ice_v(i,J)
    enddo ; enddo
  endif

  !### This needs to be replaced with an appropriate copy and average.
  fluxes%tr_fluxes => flux_tmp%tr_fluxes

end subroutine forcing_accumulate


!> Offer mechanical forcing fields for diagnostics for those
!! fields registered as part of register_forcing_type_diags.
subroutine mech_forcing_diags(fluxes, dt, G, diag, handles)
  type(forcing),         intent(in)    :: fluxes   !< fluxes type
  real,                  intent(in)    :: dt       !< time step
  type(ocean_grid_type), intent(in)    :: G        !< grid type
  type(diag_ctrl),       intent(in)    :: diag     !< diagnostic type
  type(forcing_diags),   intent(inout) :: handles  !< diagnostic id for diag_manager

  real, dimension(SZI_(G),SZJ_(G)) :: sum
  integer :: i,j,is,ie,js,je

  call cpu_clock_begin(handles%id_clock_forcing)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  if (query_averaging_enabled(diag)) then

    if ((handles%id_taux > 0) .and. ASSOCIATED(fluxes%taux)) &
      call post_data(handles%id_taux, fluxes%taux, diag)
    if ((handles%id_tauy > 0) .and. ASSOCIATED(fluxes%tauy)) &
      call post_data(handles%id_tauy, fluxes%tauy, diag)
    if ((handles%id_ustar > 0) .and. ASSOCIATED(fluxes%ustar)) &
      call post_data(handles%id_ustar, fluxes%ustar, diag)
    if (handles%id_ustar_berg > 0) &
      call post_data(handles%id_ustar_berg, fluxes%ustar_berg, diag)
    if (handles%id_area_berg > 0) &
      call post_data(handles%id_area_berg, fluxes%area_berg, diag)
    if (handles%id_mass_berg > 0) &
      call post_data(handles%id_mass_berg, fluxes%mass_berg, diag)
    if (handles%id_frac_ice_cover > 0) &
      call post_data(handles%id_frac_ice_cover, fluxes%frac_shelf_h, diag)
    if (handles%id_ustar_ice_cover > 0) &
      call post_data(handles%id_ustar_ice_cover, fluxes%ustar_shelf, diag)

  endif

  call cpu_clock_end(handles%id_clock_forcing)
end subroutine mech_forcing_diags


!> Offer buoyancy forcing fields for diagnostics for those
!! fields registered as part of register_forcing_type_diags.
subroutine forcing_diagnostics(fluxes, state, dt, G, diag, handles)
  type(forcing),         intent(in)    :: fluxes    !< flux type
  type(surface),         intent(in)    :: state     !< ocean state
  real,                  intent(in)    :: dt        !< time step
  type(ocean_grid_type), intent(in)    :: G         !< grid type
  type(diag_ctrl),       intent(in)    :: diag      !< diagnostic regulator
  type(forcing_diags),   intent(inout) :: handles   !< diagnostic ids

  ! local
  real, dimension(SZI_(G),SZJ_(G)) :: sum
  real :: total_transport ! for diagnosing integrated boundary transport
  real :: ave_flux        ! for diagnosing averaged   boundary flux
  real :: C_p             ! seawater heat capacity (J/(deg K * kg))
  real :: I_dt            ! inverse time step
  real :: ppt2mks         ! conversion between ppt and mks
  integer :: i,j,is,ie,js,je

  call cpu_clock_begin(handles%id_clock_forcing)

  C_p     = fluxes%C_p
  I_dt    = 1.0/dt
  ppt2mks = 1e-3
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec


  if (query_averaging_enabled(diag)) then


    ! post the diagnostics for surface mass fluxes ==================================

    if (handles%id_prcme > 0 .or. handles%id_total_prcme > 0 .or. handles%id_prcme_ga > 0) then
      sum(:,:) = 0.0
      if (ASSOCIATED(fluxes%lprec))       sum(:,:) = sum(:,:)+fluxes%lprec(:,:)
      if (ASSOCIATED(fluxes%fprec))       sum(:,:) = sum(:,:)+fluxes%fprec(:,:)
      ! fluxes%cond is not needed because it is derived from %evap > 0
      if (ASSOCIATED(fluxes%evap))        sum(:,:) = sum(:,:)+fluxes%evap(:,:)
      if (ASSOCIATED(fluxes%lrunoff))     sum(:,:) = sum(:,:)+fluxes%lrunoff(:,:)
      if (ASSOCIATED(fluxes%frunoff))     sum(:,:) = sum(:,:)+fluxes%frunoff(:,:)
      if (ASSOCIATED(fluxes%vprec))       sum(:,:) = sum(:,:)+fluxes%vprec(:,:)
      call post_data(handles%id_prcme, sum, diag)
      if(handles%id_total_prcme > 0) then
        total_transport = global_area_integral(sum,G)
        call post_data(handles%id_total_prcme, total_transport, diag)
      endif
      if(handles%id_prcme_ga > 0) then
        ave_flux = global_area_mean(sum,G)
        call post_data(handles%id_prcme_ga, ave_flux, diag)
      endif
    endif

    if(handles%id_net_massout > 0 .or. handles%id_total_net_massout > 0) then
      sum(:,:) = 0.0
      do j=js,je ; do i=is,ie
        if(fluxes%lprec(i,j) < 0.0) sum(i,j) = sum(i,j) + fluxes%lprec(i,j)
        if(fluxes%vprec(i,j) < 0.0) sum(i,j) = sum(i,j) + fluxes%vprec(i,j)
        if(fluxes%evap(i,j)  < 0.0) sum(i,j) = sum(i,j) + fluxes%evap(i,j)
      enddo ; enddo
      call post_data(handles%id_net_massout, sum, diag)
      if(handles%id_total_net_massout > 0) then
        total_transport = global_area_integral(sum,G)
        call post_data(handles%id_total_net_massout, total_transport, diag)
      endif
    endif

    if(handles%id_massout_flux > 0) call post_data(handles%id_massout_flux,fluxes%netMassOut,diag)

    if(handles%id_net_massin > 0 .or. handles%id_total_net_massin > 0) then
      sum(:,:) = 0.0
      do j=js,je ; do i=is,ie
        sum(i,j) = sum(i,j) + fluxes%fprec(i,j) + fluxes%lrunoff(i,j) + fluxes%frunoff(i,j)
        if(fluxes%lprec(i,j) > 0.0) sum(i,j) = sum(i,j) + fluxes%lprec(i,j)
        if(fluxes%vprec(i,j) > 0.0) sum(i,j) = sum(i,j) + fluxes%vprec(i,j)
        ! fluxes%cond is not needed because it is derived from %evap > 0
        if(fluxes%evap(i,j)  > 0.0) sum(i,j) = sum(i,j) + fluxes%evap(i,j)
      enddo ; enddo
      call post_data(handles%id_net_massin, sum, diag)
      if(handles%id_total_net_massin > 0) then
        total_transport = global_area_integral(sum,G)
        call post_data(handles%id_total_net_massin, total_transport, diag)
      endif
    endif

    if(handles%id_massin_flux > 0) call post_data(handles%id_massin_flux,fluxes%netMassIn,diag)

    if ((handles%id_evap > 0) .and. ASSOCIATED(fluxes%evap)) &
      call post_data(handles%id_evap, fluxes%evap, diag)
    if ((handles%id_total_evap > 0) .and. ASSOCIATED(fluxes%evap)) then
      total_transport = global_area_integral(fluxes%evap(:,:),G)
      call post_data(handles%id_total_evap, total_transport, diag)
    endif
    if ((handles%id_evap_ga > 0) .and. ASSOCIATED(fluxes%evap)) then
      ave_flux = global_area_mean(fluxes%evap(:,:),G)
      call post_data(handles%id_evap_ga, ave_flux, diag)
    endif

    if ((handles%id_precip > 0) .and. ASSOCIATED(fluxes%lprec) .and. ASSOCIATED(fluxes%fprec)) then
      sum(:,:) = fluxes%lprec(:,:) + fluxes%fprec(:,:)
      call post_data(handles%id_precip, sum, diag)
    endif
    if ((handles%id_total_precip > 0) .and. ASSOCIATED(fluxes%lprec) .and. ASSOCIATED(fluxes%fprec)) then
      sum(:,:) = fluxes%lprec(:,:) + fluxes%fprec(:,:)
      total_transport = global_area_integral(sum,G)
      call post_data(handles%id_total_precip, total_transport, diag)
    endif
    if ((handles%id_precip_ga > 0) .and. ASSOCIATED(fluxes%lprec) .and. ASSOCIATED(fluxes%fprec)) then
      sum(:,:) = fluxes%lprec(:,:) + fluxes%fprec(:,:)
      ave_flux = global_area_mean(sum,G)
      call post_data(handles%id_precip_ga, ave_flux, diag)
    endif

    if ((handles%id_lprec > 0) .and. ASSOCIATED(fluxes%lprec)) &
      call post_data(handles%id_lprec, fluxes%lprec, diag)
    if ((handles%id_total_lprec > 0) .and. ASSOCIATED(fluxes%lprec)) then
      total_transport = global_area_integral(fluxes%lprec(:,:),G)
      call post_data(handles%id_total_lprec, total_transport, diag)
    endif
    if ((handles%id_lprec_ga > 0) .and. ASSOCIATED(fluxes%lprec)) then
      sum(:,:) = fluxes%lprec(:,:)
      ave_flux = global_area_mean(sum,G)
      call post_data(handles%id_lprec_ga, ave_flux, diag)
    endif

    if ((handles%id_fprec > 0) .and. ASSOCIATED(fluxes%fprec)) &
      call post_data(handles%id_fprec, fluxes%fprec, diag)
    if ((handles%id_total_fprec > 0) .and. ASSOCIATED(fluxes%fprec)) then
      total_transport = global_area_integral(fluxes%fprec(:,:),G)
      call post_data(handles%id_total_fprec, total_transport, diag)
    endif
    if ((handles%id_fprec_ga > 0) .and. ASSOCIATED(fluxes%fprec)) then
      sum(:,:) = fluxes%fprec(:,:)
      ave_flux = global_area_mean(sum,G)
      call post_data(handles%id_fprec_ga, ave_flux, diag)
    endif

    if ((handles%id_vprec > 0) .and. ASSOCIATED(fluxes%vprec)) &
      call post_data(handles%id_vprec, fluxes%vprec, diag)
    if ((handles%id_total_vprec > 0) .and. ASSOCIATED(fluxes%vprec)) then
      total_transport = global_area_integral(fluxes%vprec(:,:),G)
      call post_data(handles%id_total_vprec, total_transport, diag)
    endif
    if ((handles%id_vprec_ga > 0) .and. ASSOCIATED(fluxes%vprec)) then
      sum(:,:) = fluxes%vprec(:,:)
      ave_flux = global_area_mean(sum,G)
      call post_data(handles%id_vprec_ga, ave_flux, diag)
    endif

    if ((handles%id_lrunoff > 0) .and. ASSOCIATED(fluxes%lrunoff)) &
      call post_data(handles%id_lrunoff, fluxes%lrunoff, diag)
    if ((handles%id_total_lrunoff > 0) .and. ASSOCIATED(fluxes%lrunoff)) then
      total_transport = global_area_integral(fluxes%lrunoff(:,:),G)
      call post_data(handles%id_total_lrunoff, total_transport, diag)
    endif

    if ((handles%id_frunoff > 0) .and. ASSOCIATED(fluxes%frunoff)) &
      call post_data(handles%id_frunoff, fluxes%frunoff, diag)
    if ((handles%id_total_frunoff > 0) .and. ASSOCIATED(fluxes%frunoff)) then
      total_transport = global_area_integral(fluxes%frunoff(:,:),G)
      call post_data(handles%id_total_frunoff, total_transport, diag)
    endif


    ! post diagnostics for boundary heat fluxes ====================================

    if ((handles%id_heat_content_lrunoff > 0) .and. ASSOCIATED(fluxes%heat_content_lrunoff))  &
      call post_data(handles%id_heat_content_lrunoff, fluxes%heat_content_lrunoff, diag)
    if ((handles%id_total_heat_content_lrunoff > 0) .and. ASSOCIATED(fluxes%heat_content_lrunoff)) then
      total_transport = global_area_integral(fluxes%heat_content_lrunoff(:,:),G)
      call post_data(handles%id_total_heat_content_lrunoff, total_transport, diag)
    endif

    if ((handles%id_heat_content_frunoff > 0) .and. ASSOCIATED(fluxes%heat_content_frunoff))  &
      call post_data(handles%id_heat_content_frunoff, fluxes%heat_content_frunoff, diag)
    if ((handles%id_total_heat_content_frunoff > 0) .and. ASSOCIATED(fluxes%heat_content_frunoff)) then
      total_transport = global_area_integral(fluxes%heat_content_frunoff(:,:),G)
      call post_data(handles%id_total_heat_content_frunoff, total_transport, diag)
    endif

    if ((handles%id_heat_content_lprec > 0) .and. ASSOCIATED(fluxes%heat_content_lprec))      &
      call post_data(handles%id_heat_content_lprec, fluxes%heat_content_lprec, diag)
    if ((handles%id_total_heat_content_lprec > 0) .and. ASSOCIATED(fluxes%heat_content_lprec)) then
      total_transport = global_area_integral(fluxes%heat_content_lprec(:,:),G)
      call post_data(handles%id_total_heat_content_lprec, total_transport, diag)
    endif

    if ((handles%id_heat_content_fprec > 0) .and. ASSOCIATED(fluxes%heat_content_fprec))      &
      call post_data(handles%id_heat_content_fprec, fluxes%heat_content_fprec, diag)
    if ((handles%id_total_heat_content_fprec > 0) .and. ASSOCIATED(fluxes%heat_content_fprec)) then
      total_transport = global_area_integral(fluxes%heat_content_fprec(:,:),G)
      call post_data(handles%id_total_heat_content_fprec, total_transport, diag)
    endif

    if ((handles%id_heat_content_vprec > 0) .and. ASSOCIATED(fluxes%heat_content_vprec))      &
      call post_data(handles%id_heat_content_vprec, fluxes%heat_content_vprec, diag)
    if ((handles%id_total_heat_content_vprec > 0) .and. ASSOCIATED(fluxes%heat_content_vprec)) then
      total_transport = global_area_integral(fluxes%heat_content_vprec(:,:),G)
      call post_data(handles%id_total_heat_content_vprec, total_transport, diag)
    endif

    if ((handles%id_heat_content_cond > 0) .and. ASSOCIATED(fluxes%heat_content_cond))        &
      call post_data(handles%id_heat_content_cond, fluxes%heat_content_cond, diag)
    if ((handles%id_total_heat_content_cond > 0) .and. ASSOCIATED(fluxes%heat_content_cond)) then
      total_transport = global_area_integral(fluxes%heat_content_cond(:,:),G)
      call post_data(handles%id_total_heat_content_cond, total_transport, diag)
    endif

    if ((handles%id_heat_content_massout > 0) .and. ASSOCIATED(fluxes%heat_content_massout))  &
      call post_data(handles%id_heat_content_massout, fluxes%heat_content_massout, diag)
    if ((handles%id_total_heat_content_massout > 0) .and. ASSOCIATED(fluxes%heat_content_massout)) then
      total_transport = global_area_integral(fluxes%heat_content_massout(:,:),G)
      call post_data(handles%id_total_heat_content_massout, total_transport, diag)
    endif

    if ((handles%id_heat_content_massin > 0) .and. ASSOCIATED(fluxes%heat_content_massin))  &
      call post_data(handles%id_heat_content_massin, fluxes%heat_content_massin, diag)
    if ((handles%id_total_heat_content_massin > 0) .and. ASSOCIATED(fluxes%heat_content_massin)) then
      total_transport = global_area_integral(fluxes%heat_content_massin(:,:),G)
      call post_data(handles%id_total_heat_content_massin, total_transport, diag)
    endif

    if (handles%id_net_heat_coupler > 0 .or. handles%id_total_net_heat_coupler > 0 .or. handles%id_net_heat_coupler_ga > 0. ) then
      sum(:,:) = 0.0
      if (ASSOCIATED(fluxes%LW))         sum(:,:) = sum(:,:) + fluxes%LW(:,:)
      if (ASSOCIATED(fluxes%latent))     sum(:,:) = sum(:,:) + fluxes%latent(:,:)
      if (ASSOCIATED(fluxes%sens))       sum(:,:) = sum(:,:) + fluxes%sens(:,:)
      if (ASSOCIATED(fluxes%SW))         sum(:,:) = sum(:,:) + fluxes%SW(:,:)
      call post_data(handles%id_net_heat_coupler, sum, diag)
      if(handles%id_total_net_heat_coupler > 0) then
        total_transport = global_area_integral(sum,G)
        call post_data(handles%id_total_net_heat_coupler, total_transport, diag)
      endif
      if(handles%id_net_heat_coupler_ga > 0) then
        ave_flux = global_area_mean(sum,G)
        call post_data(handles%id_net_heat_coupler_ga, ave_flux, diag)
      endif
    endif

    if (handles%id_net_heat_surface > 0 .or. handles%id_total_net_heat_surface > 0 .or. handles%id_net_heat_surface_ga > 0. ) then
      sum(:,:) = 0.0
      if (ASSOCIATED(fluxes%LW))                   sum(:,:) = sum(:,:) + fluxes%LW(:,:)
      if (ASSOCIATED(fluxes%latent))               sum(:,:) = sum(:,:) + fluxes%latent(:,:)
      if (ASSOCIATED(fluxes%sens))                 sum(:,:) = sum(:,:) + fluxes%sens(:,:)
      if (ASSOCIATED(fluxes%SW))                   sum(:,:) = sum(:,:) + fluxes%SW(:,:)
      if (ASSOCIATED(state%frazil))                sum(:,:) = sum(:,:) + state%frazil(:,:) * I_dt
    ! if (ASSOCIATED(state%TempXpme)) then
    !    sum(:,:) = sum(:,:) + state%TempXpme(:,:) * fluxes%C_p * I_dt
    ! else
      if (ASSOCIATED(fluxes%heat_content_lrunoff)) sum(:,:) = sum(:,:) + fluxes%heat_content_lrunoff(:,:)
      if (ASSOCIATED(fluxes%heat_content_frunoff)) sum(:,:) = sum(:,:) + fluxes%heat_content_frunoff(:,:)
      if (ASSOCIATED(fluxes%heat_content_lprec))   sum(:,:) = sum(:,:) + fluxes%heat_content_lprec(:,:)
      if (ASSOCIATED(fluxes%heat_content_fprec))   sum(:,:) = sum(:,:) + fluxes%heat_content_fprec(:,:)
      if (ASSOCIATED(fluxes%heat_content_vprec))   sum(:,:) = sum(:,:) + fluxes%heat_content_vprec(:,:)
      if (ASSOCIATED(fluxes%heat_content_cond))    sum(:,:) = sum(:,:) + fluxes%heat_content_cond(:,:)
      if (ASSOCIATED(fluxes%heat_content_massout)) sum(:,:) = sum(:,:) + fluxes%heat_content_massout(:,:)
    ! endif
      if (ASSOCIATED(fluxes%heat_added))         sum(:,:) = sum(:,:) + fluxes%heat_added(:,:)
      call post_data(handles%id_net_heat_surface, sum, diag)

      if(handles%id_total_net_heat_surface > 0) then
        total_transport = global_area_integral(sum,G)
        call post_data(handles%id_total_net_heat_surface, total_transport, diag)
      endif
      if(handles%id_net_heat_surface_ga > 0) then
        ave_flux = global_area_mean(sum,G)
        call post_data(handles%id_net_heat_surface_ga, ave_flux, diag)
      endif
    endif

    if (handles%id_heat_content_surfwater > 0 .or. handles%id_total_heat_content_surfwater > 0) then
      sum(:,:) = 0.0
    ! if (ASSOCIATED(state%TempXpme)) then
    !   sum(:,:) = sum(:,:) + state%TempXpme(:,:) * fluxes%C_p * I_dt
    ! else
        if (ASSOCIATED(fluxes%heat_content_lrunoff)) sum(:,:) = sum(:,:) + fluxes%heat_content_lrunoff(:,:)
        if (ASSOCIATED(fluxes%heat_content_frunoff)) sum(:,:) = sum(:,:) + fluxes%heat_content_frunoff(:,:)
        if (ASSOCIATED(fluxes%heat_content_lprec))   sum(:,:) = sum(:,:) + fluxes%heat_content_lprec(:,:)
        if (ASSOCIATED(fluxes%heat_content_fprec))   sum(:,:) = sum(:,:) + fluxes%heat_content_fprec(:,:)
        if (ASSOCIATED(fluxes%heat_content_vprec))   sum(:,:) = sum(:,:) + fluxes%heat_content_vprec(:,:)
        if (ASSOCIATED(fluxes%heat_content_cond))    sum(:,:) = sum(:,:) + fluxes%heat_content_cond(:,:)
        if (ASSOCIATED(fluxes%heat_content_massout)) sum(:,:) = sum(:,:) + fluxes%heat_content_massout(:,:)
    ! endif
      call post_data(handles%id_heat_content_surfwater, sum, diag)
      if(handles%id_total_heat_content_surfwater > 0) then
        total_transport = global_area_integral(sum,G)
        call post_data(handles%id_total_heat_content_surfwater, total_transport, diag)
      endif
    endif

    ! for OMIP, hfrunoffds = heat content of liquid plus frozen runoff
    if (handles%id_hfrunoffds > 0) then
      sum(:,:) = 0.0
      if(ASSOCIATED(fluxes%heat_content_lrunoff)) then
        sum(:,:) = sum(:,:) + fluxes%heat_content_lrunoff(:,:)
      endif
      if(ASSOCIATED(fluxes%heat_content_frunoff)) then
        sum(:,:) = sum(:,:) + fluxes%heat_content_frunoff(:,:)
      endif
      call post_data(handles%id_hfrunoffds, sum, diag)
    endif

    ! for OMIP, hfrainds = heat content of lprec + fprec + cond
    if (handles%id_hfrainds > 0) then
      sum(:,:) = 0.0
      if(ASSOCIATED(fluxes%heat_content_lprec)) then
        sum(:,:) = sum(:,:) + fluxes%heat_content_lprec(:,:)
      endif
      if(ASSOCIATED(fluxes%heat_content_fprec)) then
        sum(:,:) = sum(:,:) + fluxes%heat_content_fprec(:,:)
      endif
      if(ASSOCIATED(fluxes%heat_content_cond)) then
        sum(:,:) = sum(:,:) + fluxes%heat_content_cond(:,:)
      endif
      call post_data(handles%id_hfrainds, sum, diag)
    endif

    if ((handles%id_LwLatSens > 0) .and. ASSOCIATED(fluxes%lw) .and. &
         ASSOCIATED(fluxes%latent) .and. ASSOCIATED(fluxes%sens)) then
      sum(:,:) = (fluxes%lw(:,:) + fluxes%latent(:,:)) + fluxes%sens(:,:)
      call post_data(handles%id_LwLatSens, sum, diag)
    endif

    if ((handles%id_total_LwLatSens > 0) .and. ASSOCIATED(fluxes%lw) .and. &
         ASSOCIATED(fluxes%latent) .and. ASSOCIATED(fluxes%sens)) then
      sum(:,:) = (fluxes%lw(:,:) + fluxes%latent(:,:)) + fluxes%sens(:,:)
      total_transport = global_area_integral(sum,G)
      call post_data(handles%id_total_LwLatSens, total_transport, diag)
    endif

    if ((handles%id_LwLatSens_ga > 0) .and. ASSOCIATED(fluxes%lw) .and. &
         ASSOCIATED(fluxes%latent) .and. ASSOCIATED(fluxes%sens)) then
      sum(:,:) = (fluxes%lw(:,:) + fluxes%latent(:,:)) + fluxes%sens(:,:)
      ave_flux = global_area_mean(sum,G)
      call post_data(handles%id_LwLatSens_ga, ave_flux, diag)
    endif

    if ((handles%id_sw > 0) .and. ASSOCIATED(fluxes%sw)) then
      call post_data(handles%id_sw, fluxes%sw, diag)
    endif
    if ((handles%id_sw_vis > 0) .and. ASSOCIATED(fluxes%sw_vis_dir) .and. &
        ASSOCIATED(fluxes%sw_vis_dif)) then
      call post_data(handles%id_sw_vis, fluxes%sw_vis_dir+fluxes%sw_vis_dif, diag)
    endif
    if ((handles%id_sw_nir > 0) .and. ASSOCIATED(fluxes%sw_nir_dir) .and. &
        ASSOCIATED(fluxes%sw_nir_dif)) then
      call post_data(handles%id_sw_nir, fluxes%sw_nir_dir+fluxes%sw_nir_dif, diag)
    endif
    if ((handles%id_total_sw > 0) .and. ASSOCIATED(fluxes%sw)) then
      total_transport = global_area_integral(fluxes%sw,G)
      call post_data(handles%id_total_sw, total_transport, diag)
    endif
    if ((handles%id_sw_ga > 0) .and. ASSOCIATED(fluxes%sw)) then
      ave_flux = global_area_mean(fluxes%sw,G)
      call post_data(handles%id_sw_ga, ave_flux, diag)
    endif

    if ((handles%id_lw > 0) .and. ASSOCIATED(fluxes%lw)) then
      call post_data(handles%id_lw, fluxes%lw, diag)
    endif
    if ((handles%id_total_lw > 0) .and. ASSOCIATED(fluxes%lw)) then
      total_transport = global_area_integral(fluxes%lw,G)
      call post_data(handles%id_total_lw, total_transport, diag)
    endif
    if ((handles%id_lw_ga > 0) .and. ASSOCIATED(fluxes%lw)) then
      ave_flux = global_area_mean(fluxes%lw,G)
      call post_data(handles%id_lw_ga, ave_flux, diag)
    endif

    if ((handles%id_lat > 0) .and. ASSOCIATED(fluxes%latent)) then
      call post_data(handles%id_lat, fluxes%latent, diag)
    endif
    if ((handles%id_total_lat > 0) .and. ASSOCIATED(fluxes%latent)) then
      total_transport = global_area_integral(fluxes%latent,G)
      call post_data(handles%id_total_lat, total_transport, diag)
    endif
    if ((handles%id_lat_ga > 0) .and. ASSOCIATED(fluxes%latent)) then
      ave_flux = global_area_mean(fluxes%latent,G)
      call post_data(handles%id_lat_ga, ave_flux, diag)
    endif

    if ((handles%id_lat_evap > 0) .and. ASSOCIATED(fluxes%latent_evap_diag)) then
      call post_data(handles%id_lat_evap, fluxes%latent_evap_diag, diag)
    endif
    if ((handles%id_total_lat_evap > 0) .and. ASSOCIATED(fluxes%latent_evap_diag)) then
      total_transport = global_area_integral(fluxes%latent_evap_diag,G)
      call post_data(handles%id_total_lat_evap, total_transport, diag)
    endif

    if ((handles%id_lat_fprec > 0) .and. ASSOCIATED(fluxes%latent_fprec_diag)) then
      call post_data(handles%id_lat_fprec, fluxes%latent_fprec_diag, diag)
    endif
    if ((handles%id_total_lat_fprec > 0) .and. ASSOCIATED(fluxes%latent_fprec_diag)) then
      total_transport = global_area_integral(fluxes%latent_fprec_diag,G)
      call post_data(handles%id_total_lat_fprec, total_transport, diag)
    endif

    if ((handles%id_lat_frunoff > 0) .and. ASSOCIATED(fluxes%latent_frunoff_diag)) then
      call post_data(handles%id_lat_frunoff, fluxes%latent_frunoff_diag, diag)
    endif
    if(handles%id_total_lat_frunoff > 0 .and. ASSOCIATED(fluxes%latent_frunoff_diag)) then
      total_transport = global_area_integral(fluxes%latent_frunoff_diag,G)
      call post_data(handles%id_total_lat_frunoff, total_transport, diag)
    endif

    if ((handles%id_sens > 0) .and. ASSOCIATED(fluxes%sens)) then
      call post_data(handles%id_sens, fluxes%sens, diag)
    endif
    if ((handles%id_total_sens > 0) .and. ASSOCIATED(fluxes%sens)) then
      total_transport = global_area_integral(fluxes%sens,G)
      call post_data(handles%id_total_sens, total_transport, diag)
    endif
    if ((handles%id_sens_ga > 0) .and. ASSOCIATED(fluxes%sens)) then
      ave_flux = global_area_mean(fluxes%sens,G)
      call post_data(handles%id_sens_ga, ave_flux, diag)
    endif

    if ((handles%id_heat_added > 0) .and. ASSOCIATED(fluxes%heat_added)) then
      call post_data(handles%id_heat_added, fluxes%heat_added, diag)
    endif

    if ((handles%id_total_heat_added > 0) .and. ASSOCIATED(fluxes%heat_added)) then
      total_transport = global_area_integral(fluxes%heat_added,G)
      call post_data(handles%id_total_heat_added, total_transport, diag)
    endif


    ! post the diagnostics for boundary salt fluxes ==========================

    if ((handles%id_saltflux > 0) .and. ASSOCIATED(fluxes%salt_flux)) &
      call post_data(handles%id_saltflux, fluxes%salt_flux, diag)
    if ((handles%id_total_saltflux > 0) .and. ASSOCIATED(fluxes%salt_flux)) then
      total_transport = ppt2mks*global_area_integral(fluxes%salt_flux,G)
      call post_data(handles%id_total_saltflux, total_transport, diag)
    endif

    if ((handles%id_saltFluxAdded > 0) .and. ASSOCIATED(fluxes%salt_flux_added)) &
      call post_data(handles%id_saltFluxAdded, fluxes%salt_flux_added, diag)
    if ((handles%id_total_saltFluxAdded > 0) .and. ASSOCIATED(fluxes%salt_flux_added)) then
      total_transport = ppt2mks*global_area_integral(fluxes%salt_flux_added,G)
      call post_data(handles%id_total_saltFluxAdded, total_transport, diag)
    endif

    if (handles%id_saltFluxIn > 0 .and. ASSOCIATED(fluxes%salt_flux_in)) &
      call post_data(handles%id_saltFluxIn, fluxes%salt_flux_in, diag)
    if ((handles%id_total_saltFluxIn > 0) .and. ASSOCIATED(fluxes%salt_flux_in)) then
      total_transport = ppt2mks*global_area_integral(fluxes%salt_flux_in,G)
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

    if ((handles%id_psurf > 0) .and. ASSOCIATED(fluxes%p_surf))                      &
      call post_data(handles%id_psurf, fluxes%p_surf, diag)

    if ((handles%id_TKE_tidal > 0) .and. ASSOCIATED(fluxes%TKE_tidal))               &
      call post_data(handles%id_TKE_tidal, fluxes%TKE_tidal, diag)

    if ((handles%id_buoy > 0) .and. ASSOCIATED(fluxes%buoy))                         &
      call post_data(handles%id_buoy, fluxes%buoy, diag)


  endif

  call cpu_clock_end(handles%id_clock_forcing)
end subroutine forcing_diagnostics


!> Conditionally allocate fields within the forcing type
subroutine allocate_forcing_type(G, fluxes, stress, ustar, water, heat, shelf, press, iceberg)
  type(ocean_grid_type), intent(in) :: G       !< Ocean grid structure
  type(forcing),      intent(inout) :: fluxes  !< Forcing fields structure
  logical, optional,     intent(in) :: stress  !< If present and true, allocate taux, tauy
  logical, optional,     intent(in) :: ustar   !< If present and true, allocate ustar
  logical, optional,     intent(in) :: water   !< If present and true, allocate water fluxes
  logical, optional,     intent(in) :: heat    !< If present and true, allocate heat fluxes
  logical, optional,     intent(in) :: shelf   !< If present and true, allocate fluxes for ice-shelf
  logical, optional,     intent(in) :: press   !< If present and true, allocate p_surf
  logical, optional,     intent(in) :: iceberg !< If present and true, allocate fluxes for icebergs

  ! Local variables
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  logical :: heat_water

  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB

  call myAlloc(fluxes%taux,IsdB,IedB,jsd,jed, stress)
  call myAlloc(fluxes%tauy,isd,ied,JsdB,JedB, stress)
  call myAlloc(fluxes%ustar,isd,ied,jsd,jed, ustar)

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

  call myAlloc(fluxes%sw,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%lw,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%latent,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%sens,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%latent_evap_diag,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%latent_fprec_diag,isd,ied,jsd,jed, heat)
  call myAlloc(fluxes%latent_frunoff_diag,isd,ied,jsd,jed, heat)

  if (present(heat) .and. present(water)) then ; if (heat .and. water) then
    call myAlloc(fluxes%heat_content_cond,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_lprec,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_fprec,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_vprec,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_lrunoff,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_frunoff,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_massout,isd,ied,jsd,jed, .true.)
    call myAlloc(fluxes%heat_content_massin,isd,ied,jsd,jed, .true.)
  endif ; endif

  call myAlloc(fluxes%frac_shelf_h,isd,ied,jsd,jed, shelf)
  call myAlloc(fluxes%frac_shelf_u,IsdB,IedB,jsd,jed, shelf)
  call myAlloc(fluxes%frac_shelf_v,isd,ied,JsdB,JedB, shelf)
  call myAlloc(fluxes%ustar_shelf,isd,ied,jsd,jed, shelf)
  call myAlloc(fluxes%iceshelf_melt,isd,ied,jsd,jed, shelf)
  call myAlloc(fluxes%rigidity_ice_u,IsdB,IedB,jsd,jed, shelf)
  call myAlloc(fluxes%rigidity_ice_v,isd,ied,JsdB,JedB, shelf)

  call myAlloc(fluxes%p_surf,isd,ied,jsd,jed, press)

  !These fields should only on allocated when iceberg area is being passed through the coupler.
  call myAlloc(fluxes%ustar_berg,isd,ied,jsd,jed, iceberg)
  call myAlloc(fluxes%area_berg,isd,ied,jsd,jed, iceberg)
  call myAlloc(fluxes%mass_berg,isd,ied,jsd,jed, iceberg)
  contains

  !> Allocates and zeroes-out array.
  subroutine myAlloc(array, is, ie, js, je, flag)
    real, dimension(:,:), pointer :: array !< Array to be allocated
    integer,           intent(in) :: is !< Start i-index
    integer,           intent(in) :: ie !< End i-index
    integer,           intent(in) :: js !< Start j-index
    integer,           intent(in) :: je !< End j-index
    logical, optional, intent(in) :: flag !< Flag to indicate to allocate

    if (present(flag)) then
      if (flag) then
        if (.not.associated(array)) then
          ALLOCATE(array(is:ie,js:je))
          array(is:ie,js:je) = 0.0
        endif
      endif
    endif
  end subroutine myAlloc

end subroutine allocate_forcing_type


!> Deallocate the forcing type
subroutine deallocate_forcing_type(fluxes)
  type(forcing), intent(inout) :: fluxes !< Forcing fields structure

  if (associated(fluxes%taux))                 deallocate(fluxes%taux)
  if (associated(fluxes%tauy))                 deallocate(fluxes%tauy)
  if (associated(fluxes%ustar))                deallocate(fluxes%ustar)
  if (associated(fluxes%buoy))                 deallocate(fluxes%buoy)
  if (associated(fluxes%sw))                   deallocate(fluxes%sw)
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
  if (associated(fluxes%frac_shelf_u))         deallocate(fluxes%frac_shelf_u)
  if (associated(fluxes%frac_shelf_v))         deallocate(fluxes%frac_shelf_v)
  if (associated(fluxes%rigidity_ice_u))       deallocate(fluxes%rigidity_ice_u)
  if (associated(fluxes%rigidity_ice_v))       deallocate(fluxes%rigidity_ice_v)
  if (associated(fluxes%tr_fluxes))            deallocate(fluxes%tr_fluxes)
  if (associated(fluxes%ustar_berg))           deallocate(fluxes%ustar_berg)
  if (associated(fluxes%area_berg))            deallocate(fluxes%area_berg)
  if (associated(fluxes%mass_berg))            deallocate(fluxes%mass_berg)
end subroutine deallocate_forcing_type


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
!!                      LW+LAT+SENS in net_heat inside routine
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
!! \f{eqnarray*}{
!!  Q(k=1) &=& \mbox{hfds} - \mbox{pen_SW(leaving bottom of k=1)}
!!   \\    &=& \mbox{nonpen_SW} + (\mbox{pen_SW(enter k=1)}-\mbox{pen_SW(leave k=1)})
!!                              + \mbox{LW+LAT+SENS+MASS+FRAZ+RES}
!!   \\    &=& \mbox{nonpen_SW}+ \mbox{LW+LAT+SENS+MASS+FRAZ+RES}
!!                + [\mbox{pen_SW(enter k=1)} - \mbox{pen_SW(leave k=1)}]
!!   \f}
!! The convergence of the penetrative shortwave flux is given by
!! \f$ \mbox{pen_SW (enter k)}-\mbox{pen_SW (leave k)}\f$.  This term
!! appears for all cells k=1,nz.  It is diagnosed as "rsdoabsorb" inside module
!! MOM6/src/parameterizations/vertical/MOM_diabatic_aux.F90
!!

end module MOM_forcing_type
