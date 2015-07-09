module MOM_forcing_type
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

use MOM_checksums,     only : hchksum, qchksum, uchksum, vchksum, is_NaN
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, register_scalar_field
use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_alloc, query_averaging_enabled
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_EOS,           only : calculate_density_derivs
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_shortwave_abs, only : absorbRemainingSW, sumSWoverBands, optics_type
use MOM_spatial_means, only : global_area_integral 
use MOM_variables,     only : surface, thermo_var_ptrs

use coupler_types_mod, only : coupler_2d_bc_type

implicit none ; private

#include <MOM_memory.h>

public extractFluxes1d, extractFluxes2d, MOM_forcing_chksum, optics_type
public calculateBuoyancyFlux1d, calculateBuoyancyFlux2d, forcing_accumulate
public forcing_SinglePointPrint, mech_forcing_diags, forcing_diagnostics
public register_forcing_type_diags, allocate_forcing_type, deallocate_forcing_type

integer :: num_msg = 0
integer :: max_msg = 2

type, public :: forcing
  ! This structure contains pointers to the boundary 
  ! forcing used to drive the liquid ocean as part of MOM.  
  ! All fluxes are positive into the ocean. For surface
  ! fluxes, that means fluxes are positive downward.
  !
  ! Pointers should be initialized to NULL.
  !
  ! The data in this type is allocated in the module
  ! MOM_surface_forcing.F90, of which there are three:
  ! solo, coupled and ice-shelf.  Alternatively, they are 
  ! allocated in MESO_surface_forcing.F90, which is a 
  ! special case of solo_driver/MOM_surface_forcing.F90.

  real, pointer, dimension(:,:) :: &

    ! surface stress components and turbulent velocity scale  
    taux          => NULL(), & ! zonal wind stress (Pa)
    tauy          => NULL(), & ! meridional wind stress (Pa)
    ustar         => NULL(), & ! surface friction velocity scale (m/s)

    ! surface buoyancy force 
    buoy          => NULL(), & ! buoyancy flux (m^2/s^3)

    ! radiative heat fluxes into the ocean (W/m^2) 
    sw            => NULL(), & ! shortwave (W/m^2)
    sw_vis_dir    => NULL(), & ! visible, direct shortwave (W/m^2)
    sw_vis_dif    => NULL(), & ! visible, diffuse shortwave (W/m^2)
    sw_nir_dir    => NULL(), & ! near-IR, direct shortwave (W/m^2)
    sw_nir_dif    => NULL(), & ! near-IR, diffuse shortwave (W/m^2)
    lw            => NULL(), & ! longwave (W/m^2) (typically negative)

    ! turbulent heat fluxes into the ocean (W/m^2) 
    latent         => NULL(), & ! latent (W/m^2)   (typically < 0)
    sens           => NULL(), & ! sensible (W/m^2) (typically negative)
    heat_restore   => NULL(), & ! heat flux from SST restoring (W/m^2) in idealized simulations

    ! components of latent heat fluxes used for diagnostic purposes 
    latent_evap_diag    => NULL(), & ! latent (W/m^2) from evaporating liquid water (typically < 0)
    latent_fprec_diag   => NULL(), & ! latent (W/m^2) from melting fprec  (typically < 0)
    latent_frunoff_diag => NULL(), & ! latent (W/m^2) from melting frunoff (calving) (typically < 0)

    ! water mass fluxes into the ocean ( kg/(m^2 s) )
    ! these mass fluxes impact the ocean mass
    evap          => NULL(), & ! (-1)*fresh water flux evaporated out of the ocean ( kg/(m^2 s) )
    lprec         => NULL(), & ! precipitating liquid water into the ocean ( kg/(m^2 s) )
    fprec         => NULL(), & ! precipitating frozen water into the ocean ( kg/(m^2 s) )
    vprec         => NULL(), & ! virtual liquid precip associated w/ SSS restoring ( kg/(m^2 s) )
    lrunoff       => NULL(), & ! liquid river runoff entering ocean ( kg/(m^2 s) )
    frunoff       => NULL(), & ! frozen river runoff (calving) entering ocean ( kg/(m^2 s) )
    seaice_melt   => NULL(), & ! seaice melt (positive) or formation (negative) ( kg/(m^2 s) )

    ! heat associated with water crossing ocean surface 
    heat_content_cond    => NULL(), & ! heat content associated with condensating water (W/m^2)
    heat_content_lprec   => NULL(), & ! heat content associated with liquid >0 precip   (W/m^2) (diagnostic)
    heat_content_fprec   => NULL(), & ! heat content associated with frozen precip      (W/m^2)
    heat_content_vprec   => NULL(), & ! heat content associated with virtual >0 precip  (W/m^2)
    heat_content_lrunoff => NULL(), & ! heat content associated with liquid runoff      (W/m^2)
    heat_content_frunoff => NULL(), & ! heat content associated with frozen runoff      (W/m^2)
    heat_content_icemelt => NULL(), & ! heat content associated with liquid sea ice     (W/m^2)
    heat_content_massout => NULL(), & ! heat content associated with mass leaving ocean (W/m^2)
    heat_content_massin  => NULL(), & ! heat content associated with mass entering ocean (W/m^2)

    ! salt mass flux (contributes to ocean mass only if non-Bouss )
    salt_flux         => NULL(), & ! net salt flux into the ocean ( kg salt/(m^2 s) )
    salt_flux_in      => NULL(), & ! salt flux provided to the ocean from coupler ( kg salt/(m^2 s) )
    salt_flux_restore => NULL(), & ! restoring piece of salt flux before adjustment 
                                   ! to net zero ( kg salt/(m^2 s) )

    ! applied surface pressure from other component models (e.g., atmos, sea ice, land ice)
    p_surf_full   => NULL(), & ! pressure at the top ocean interface (Pa).
                               ! if there is sea-ice, then p_surf_flux is at ice-ocean interface
    p_surf        => NULL(), & ! pressure at top ocean interface (Pa) as used to drive the ocean model. 
                               ! if p_surf is limited, then p_surf may be smaller than p_surf_full,
                               ! otherwise they are the same.

    ! tide related inputs 
    TKE_tidal     => NULL(), & ! tidal energy source driving mixing in bottom boundary layer (W/m^2)
    ustar_tidal   => NULL(), & ! tidal contribution to bottom ustar (m/s)

    ! land ice-shelf related inputs 
    ustar_shelf   => NULL(), & ! friction velocity under ice-shelves (m/s)
                               ! as computed by the ocean at the previous time step.
    frac_shelf_h  => NULL(), & ! Fractional ice shelf coverage of h-, u-, and v-
    frac_shelf_u  => NULL(), & ! cells, nondimensional from 0 to 1. These are only
    frac_shelf_v  => NULL(), & ! associated if ice shelves are enabled, and are
                               ! exactly 0 away from shelves or on land.
    rigidity_ice_u => NULL(),& ! Depth-integrated lateral viscosity of
    rigidity_ice_v => NULL()   ! ice shelves at u- or v-points (m3/s)

  ! Scalars set by surface forcing modules
  real :: vPrecGlobalAdj     ! adjustment to restoring vprec to zero out global net ( kg/(m^2 s) )
  real :: saltFluxGlobalAdj  ! adjustment to restoring salt flux to zero out global net ( kg salt/(m^2 s) )
  real :: netFWGlobalAdj     ! adjustment to net fresh water to zero out global net ( kg/(m^2 s) )
  real :: vPrecGlobalScl     ! scaling of restoring vprec to zero out global net ( -1..1 )
  real :: saltFluxGlobalScl  ! scaling of restoring salt flux to zero out global net ( -1..1 )
  real :: netFWGlobalScl     ! scaling of net fresh water to zero out global net ( -1..1 )

  logical :: fluxes_used = .true. ! If true, all of the heat, salt, and mass
                             ! fluxes have been applied to the ocean.
  real :: dt_buoy_accum  = -1.0 ! The amount of time over which the buoyancy fluxes
                             ! should be applied, in s.  If negative, this forcing
                             ! type variable has not yet been inialized.

  ! heat capacity
  real :: C_p                ! heat capacity of seawater ( J/(K kg) )
                             ! C_p is is the same value as in thermovar_ptrs_type.

  ! passive tracer surface fluxes 
  type(coupler_2d_bc_type), pointer :: tr_fluxes  => NULL()
     ! This structure may contain an array of named fields used for passive tracer fluxes.
     ! All arrays in tr_fluxes use the coupler indexing, which has no halos. This is not 
     ! a convenient convention, but imposed on MOM6 by the coupler.  

end type forcing

type, public :: forcing_diags
  ! id handles for the forcing type

  ! mass flux diagnostic handles 
  integer :: id_prcme        = -1, id_evap        = -1
  integer :: id_precip       = -1, id_vprec       = -1
  integer :: id_lprec        = -1, id_fprec       = -1
  integer :: id_lrunoff      = -1, id_frunoff     = -1
  integer :: id_net_massout  = -1, id_net_massin  = -1
  integer :: id_seaice_melt  = -1

  ! global area integrated mass flux diagnostic handles
  integer :: id_total_prcme        = -1, id_total_evap        = -1
  integer :: id_total_precip       = -1, id_total_vprec       = -1
  integer :: id_total_lprec        = -1, id_total_fprec       = -1
  integer :: id_total_lrunoff      = -1, id_total_frunoff     = -1
  integer :: id_total_net_massout  = -1, id_total_net_massin  = -1
  integer :: id_total_seaice_melt  = -1

  ! heat flux diagnostic handles 
  integer :: id_net_heat_coupler    = -1, id_net_heat_surface      = -1
  integer :: id_sens                = -1, id_LwLatSens             = -1
  integer :: id_sw                  = -1, id_lw                    = -1
  integer :: id_lat_evap            = -1, id_lat_frunoff           = -1
  integer :: id_lat                 = -1, id_lat_fprec             = -1
  integer :: id_heat_content_lrunoff= -1, id_heat_content_frunoff  = -1
  integer :: id_heat_content_lprec  = -1, id_heat_content_fprec    = -1
  integer :: id_heat_content_cond   = -1, id_heat_content_surfwater= -1
  integer :: id_heat_content_vprec  = -1, id_heat_content_massout  = -1
  integer :: id_heat_restore        = -1, id_heat_content_massin   = -1

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
  integer :: id_total_heat_restore        = -1, id_total_heat_content_massin   = -1


  ! salt flux diagnostic handles 
  integer :: id_saltflux          = -1
  integer :: id_saltFluxIn        = -1
  integer :: id_saltFluxRestore   = -1 

  integer :: id_total_saltflux        = -1 
  integer :: id_total_saltFluxIn      = -1 
  integer :: id_total_saltFluxRestore = -1 

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

end type forcing_diags

contains

!> Extract fluxes from surface fluxes type. 
subroutine extractFluxes1d(G, fluxes, optics, nsw, j, dt,                               &
                  DepthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
                  h, T, netMassInOut, netMassOut, net_heat, net_salt, pen_SW_bnd, tv,   &
                  aggregate_FW_forcing)

!  This subroutine extracts fluxes from the surface fluxes type. It works on a j-row
!  for optimization purposes.  The 2d (i,j) wrapper is the next subroutine below. 

  type(ocean_grid_type),          intent(in)    :: G
  type(forcing),                  intent(inout) :: fluxes 
  type(optics_type),              pointer       :: optics
  integer,                        intent(in)    :: nsw
  integer,                        intent(in)    :: j
  real,                           intent(in)    :: dt
  real,                           intent(in)    :: DepthBeforeScalingFluxes
  logical,                        intent(in)    :: useRiverHeatContent
  logical,                        intent(in)    :: useCalvingHeatContent
  real, dimension(NIMEM_,NKMEM_), intent(in)    :: h
  real, dimension(NIMEM_,NKMEM_), intent(in)    :: T
  real, dimension(NIMEM_),        intent(out)   :: netMassInOut
  real, dimension(NIMEM_),        intent(out)   :: netMassOut 
  real, dimension(NIMEM_),        intent(out)   :: net_heat
  real, dimension(NIMEM_),        intent(out)   :: net_salt
  real, dimension(:,:),           intent(out)   :: pen_SW_bnd
  type(thermo_var_ptrs),          intent(inout) :: tv
  logical,                        intent(in)    :: aggregate_FW_forcing

!  (in)      G                        = ocean grid structure
!  (in)      fluxes                   = structure containing pointers to possible
!                                       forcing fields.  Unused fields have NULL ptrs.
!  (in)      nsw                      = number of bands of penetrating shortwave radiation
!  (in)      j                        = j-index to work on
!  (in)      dt                       = time step in seconds 
!  (in)      DepthBeforeScalingFluxes =  minimum ocean thickness to allow before scaling away fluxes in H
!  (in)      h                        =  layer thickness, in m for Bouss or (kg/m^2) for non-Bouss
!  (in)      T                        =  layer temperatures, in deg C

!  (out)     netMassInOut     = net mass flux (if non-Boussinesq) or volume flux (if Boussinesq)
!                               of water in/out of ocean over a time step (H units)
!  (out)     netMassOut       = net mass flux (if non-Boussinesq) or volume flux (if Boussinesq)
!                               of water leaving ocean surface over a time step (H units).
!                               netMassOut < 0 means mass leaves ocean.  
!  (out)     net_heat         = net heat at the surface over a time step associated with coupler
!                               and restoring. We exclude two terms form net_heat: (1) heat that 
!                               can leave bottom of surface cell via penetrative SW, (2) 
!                               evaporation heat content, since do not yet know temp of evaporation. 
!                               Units are (K * H).
!  (out)     net_salt         = surface salt flux into the ocean over a time step (psu * H)
!  (out)     pen_SW_bnd       = penetrating shortwave heating at the sea surface
!                               in each penetrating band, in K H, size nsw x NIMEM_.
!  (inout)   tv               = structure containing pointers to any available
!                               thermodynamic fields. Here it is used to keep track of the
!                               heat flux associated with net mass fluxes into the ocean.

  real :: htot(SZI_(G))        ! total ocean depth (m for Bouss or kg/m^2 for non-Bouss)
  real :: Pen_sw_tot(SZI_(G))  ! sum across all bands of Pen_SW (K * H)
  real :: Ih_limit             ! inverse depth at which surface fluxes start to be limited (1/H)
  real :: scale                ! scale scales away fluxes if depth < DepthBeforeScalingFluxes
  real :: J_m2_to_H            ! converts J/m^2 to H units (m for Bouss and kg/m^2 for non-Bouss)
  real :: Irho0                !  1.0 / Rho0
  real :: I_Cp                 !  1.0 / C_p

  character(len=200) :: mesg
  integer            :: is, ie, nz, i, k, n
  Ih_limit    = 1.0 / DepthBeforeScalingFluxes
  Irho0       = 1.0 / G%Rho0 
  I_Cp        = 1.0 / fluxes%C_p
  J_m2_to_H   = 1.0 / (G%H_to_kg_m2 * fluxes%C_p)

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
    if (.not.G%Boussinesq .and. ASSOCIATED(fluxes%salt_flux)) then 
      netMassInOut(i) = netMassInOut(i) + (dt * G%kg_m2_to_H) * (scale * fluxes%salt_flux(i,j))
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
    netMassInOut(i) = G%kg_m2_to_H * netMassInOut(i)
    netMassOut(i)   = G%kg_m2_to_H * netMassOut(i)


    ! surface heat fluxes from radiation and turbulent fluxes (K * H)
    ! (H=m for Bouss, H=kg/m2 for non-Bouss)
    net_heat(i) = scale * dt * J_m2_to_H * &
                  ( fluxes%sw(i,j) +  ((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)) )

    ! Add heat flux from surface damping (restoring) (K * H).
    if (ASSOCIATED(fluxes%heat_restore)) then 
       net_heat(i) = net_heat(i) + (scale * (dt * J_m2_to_H)) * fluxes%heat_restore(i,j)
    endif 

! smg: old code 
    if (useRiverHeatContent) then
      ! remove lrunoff*SST here, to counteract its addition elsewhere
      net_heat(i) = (net_heat(i) + (scale*(dt*J_m2_to_H)) * fluxes%heat_content_lrunoff(i,j)) - &
                     (G%kg_m2_to_H * (scale * dt)) * fluxes%lrunoff(i,j) * T(i,1)
      if (ASSOCIATED(tv%TempxPmE)) then
        tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + (scale * dt) * &
            (I_Cp*fluxes%heat_content_lrunoff(i,j) - fluxes%lrunoff(i,j)*T(i,1))
      endif
    endif

! smg: old code 
    if (useCalvingHeatContent) then
      ! remove frunoff*SST here, to counteract its addition elsewhere
      net_heat(i) = net_heat(i) + (scale*(dt*J_m2_to_H)) * fluxes%heat_content_frunoff(i,j) - &
                    (G%kg_m2_to_H * (scale * dt)) * fluxes%frunoff(i,j) * T(i,1)
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

    if (num_msg < max_msg) then
      if (Pen_SW_tot(i) > 1.000001*J_m2_to_H*scale*dt*fluxes%sw(i,j)) then
        num_msg = num_msg + 1
        write(mesg,'("Penetrating shortwave of ",1pe17.10, &
                    &" exceeds total shortwave of ",1pe17.10,&
                    &" at ",1pg11.4,"E, "1pg11.4,"N.")') &
               Pen_SW_tot(i),J_m2_to_H*scale*dt*fluxes%sw(i,j),&
               G%geoLonT(i,j),G%geoLatT(i,j)
        call MOM_error(WARNING,mesg) 
      endif
    endif

    ! remove penetrative portion of the SW that exits through bottom of the top cell
    net_heat(i) = net_heat(i) - Pen_SW_tot(i)


    ! Salt fluxes
    Net_salt(i) = 0.0
    ! Convert salt_flux from kg (salt)/(m^2 * s) to 
    ! Boussinesq: (ppt * m)
    ! non-Bouss:  (g/m^2)
    if (ASSOCIATED(fluxes%salt_flux)) &
      Net_salt(i) = (scale * dt * (1000.0 * fluxes%salt_flux(i,j))) * G%kg_m2_to_H

    ! Diagnostics follow...

    ! Initialize heat_content_massin that is diagnosed in mixedlayer_convection or
    ! applyBoundaryFluxes such that the meaning is as the sum of all incoming components.
    if(ASSOCIATED(fluxes%heat_content_massin))  then
      if (aggregate_FW_forcing) then
        if (netMassInOut(i) > 0.0) then ! net is "in"
          fluxes%heat_content_massin(i,j) = -fluxes%C_p * netMassOut(i) * T(i,1) * G%H_to_kg_m2 / dt
        else ! net is "out"
          fluxes%heat_content_massin(i,j) = fluxes%C_p * ( netMassInout(i) - netMassOut(i) ) * T(i,1) * G%H_to_kg_m2 / dt
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
          fluxes%heat_content_massout(i,j) = fluxes%C_p * netMassOut(i) * T(i,1) * G%H_to_kg_m2 / dt
        else ! net is "out"
          fluxes%heat_content_massout(i,j) = -fluxes%C_p * ( netMassInout(i) - netMassOut(i) ) * T(i,1) * G%H_to_kg_m2 / dt
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

    ! Assume liquid runoff enters ocean at SST if land model does not provide runoff heat content.
    if (.not. useRiverHeatContent) then
      if (ASSOCIATED(fluxes%lrunoff) .and. ASSOCIATED(fluxes%heat_content_lrunoff)) then
        fluxes%heat_content_lrunoff(i,j) = fluxes%C_p*fluxes%lrunoff(i,j)*T(i,1)
      endif
    endif

    ! Assume solid runoff enters ocean at 0degC if land model does not provide calving heat content.
    if (.not. useCalvingHeatContent) then
      if (ASSOCIATED(fluxes%heat_content_lrunoff)) then
        fluxes%heat_content_frunoff(i,j) = 0.0
      endif
    endif

  enddo ! i-loop

end subroutine extractFluxes1d


!> 2d wrapper for 1d extract fluxes from surface fluxes type. 
subroutine extractFluxes2d(G, fluxes, optics, nsw, dt,                                  &
                  DepthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
                  h, T, netMassInOut, netMassOut, net_heat, Net_salt, Pen_SW_bnd, tv,   &
                  aggregate_FW_forcing)

  !  This subroutine extracts fluxes from the surface fluxes type. It is a 
  !  wrapper for the 1d routine extractFluxes1d.

  type(ocean_grid_type),                 intent(in)    :: G
  type(forcing),                         intent(inout) :: fluxes
  type(optics_type),                     pointer       :: optics
  integer,                               intent(in)    :: nsw
  real,                                  intent(in)    :: dt
  real,                                  intent(in)    :: DepthBeforeScalingFluxes
  logical,                               intent(in)    :: useRiverHeatContent
  logical,                               intent(in)    :: useCalvingHeatContent
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: T
  real, dimension(NIMEM_,NJMEM_),        intent(out)   :: netMassInOut
  real, dimension(NIMEM_,NJMEM_),        intent(out)   :: netMassOut
  real, dimension(NIMEM_,NJMEM_),        intent(out)   :: net_heat
  real, dimension(NIMEM_,NJMEM_),        intent(out)   :: net_salt
  real, dimension(:,:,:),                intent(out)   :: pen_SW_bnd
  type(thermo_var_ptrs),                 intent(inout) :: tv
  logical,                               intent(in)    :: aggregate_FW_forcing


!  (in)      G                        = ocean grid structure
!  (in)      fluxes                   = structure containing pointers to possible
!                                       forcing fields.  Unused fields have NULL ptrs.
!  (in)      nsw                      = number of bands of penetrating shortwave radiation
!  (in)      dt                       = time step in seconds
!  (in)      DepthBeforeScalingFluxes =  minimum ocean thickness to allow before scaling away fluxes in H
!  (in)      h                        =  layer thickness, in m for Bouss or (kg/m^2) for non-Bouss
!  (in)      T                        =  layer temperatures, in deg C

!  (out)     netMassInOut = net mass flux (if non-Boussinesq) or volume flux (if Boussinesq)
!                           of water in/out of ocean surface over a time step (H)
!  (out)     netMassOut   = net mass flux (if non-Boussinesq) or volume flux (if Boussinesq)
!                           of water leaving ocean surface over a time step (H)
!  (out)     net_heat     = net heating at the surface over a time step associated with coupler
!                           and restoring;  i.e., net_heat=SW+LW+Latent+Sensible+river (K * H).
!                           This term misses the heat from precip-evap.
!  (out)     net_salt     = surface salt flux into the ocean over a time step (psu * H)
!  (out)     pen_SW_bnd   = penetrating shortwave heating at the sea surface
!                           in each penetrating band, in K H, size nsw x NIMEM_.
!  (inout)   tv           = structure containing pointers to any available
!                           thermodynamic fields. Here it is used to keep track of the
!                           heat flux associated with net mass fluxes into the ocean.

  integer :: j
!$OMP parallel do default(none) shared(G,fluxes, optics, nsw,dt,DepthBeforeScalingFluxes, &
!$OMP                                  useRiverHeatContent, useCalvingHeatContent,        &
!$OMP                                  h,T,netMassInOut,netMassOut,Net_heat,Net_salt,Pen_SW_bnd,tv, &
!$OMP                                  aggregate_FW_forcing)
  do j=G%jsc, G%jec
    call extractFluxes1d(G, fluxes, optics, nsw, j, dt,                          &
            DepthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent,&
            h(:,j,:), T(:,j,:), netMassInOut(:,j), netMassOut(:,j),              &
            net_heat(:,j), net_salt(:,j), pen_SW_bnd(:,:,j), tv, aggregate_FW_forcing)
  enddo

end subroutine extractFluxes2d


!> Compute surface buoyancy fluxes. 
subroutine calculateBuoyancyFlux1d(G, fluxes, optics, h, Temp, Salt, tv, j, &
                                   buoyancyFlux, netHeatMinusSW, netSalt )

  ! This subtourine calculates the surface buoyancy flux by adding up the heat, 
  ! FW and salt fluxes and linearizing about the surface state.

  type(ocean_grid_type),                 intent(in)    :: G              ! ocean grid
  type(forcing),                         intent(inout) :: fluxes         ! surface fluxes
  type(optics_type),                     pointer       :: optics         ! penetrating SW optics 
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h              ! layer thickness (H)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: Temp           ! prognostic temp(deg C)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: Salt           ! salinity (ppt)
  type(thermo_var_ptrs),                 intent(inout) :: tv             ! thermodynamics type
  integer,                               intent(in)    :: j              ! j-row to work on 
  real, dimension(NIMEM_,NK_INTERFACE_), intent(inout) :: buoyancyFlux   ! buoyancy flux (m^2/s^3)
  real, dimension(NIMEM_),               intent(inout) :: netHeatMinusSW ! surf Heat flux (K H)
  real, dimension(NIMEM_),               intent(inout) :: netSalt        ! surf salt flux (ppt H)


  ! Local variables
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

  depthBeforeScalingFluxes = max( G%Angstrom, 1.e-30*G%m_to_H )
  pressure(:) = 0. ! Ignore atmospheric pressure
  GoRho       = G%g_Earth / G%Rho0
  start       = 1 + G%isc - G%isd
  npts        = 1 + G%iec - G%isc

  H_limit_fluxes = depthBeforeScalingFluxes

  ! The surface forcing is contained in the fluxes type.
  ! We aggregate the thermodynamic forcing for a time step into the following:
  ! netH       = water (H units) added/removed via surface fluxes
  ! netHeat    = heat (degC * H) via surface fluxes
  ! netSalt    = salt ( g(salt)/m2 for non-Bouss and ppt*m for Bouss ) via surface fluxes
  call extractFluxes1d(G, fluxes, optics, nsw, j, dt,                                 &
                depthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
                h(:,j,:), Temp(:,j,:), netH, netEvap, netHeatMinusSW,                 &
                netSalt, penSWbnd, tv, .false.)

  ! Sum over bands and attenuate as a function of depth
  ! netPen is the netSW as a function of depth
  call sumSWoverBands(G, h(:,j,:), optics%opacity_band(:,:,j,:), nsw, j, dt, &
                      H_limit_fluxes, .true., penSWbnd, netPen)

  ! Density derivatives
  call calculate_density_derivs(Temp(:,j,1), Salt(:,j,1), pressure, &
                                dRhodT, dRhodS, start, npts, tv%eqn_of_state)

  ! Adjust netSalt to reflect dilution effect of FW flux
  netSalt(G%isc:G%iec) = netSalt(G%isc:G%iec) - Salt(G%isc:G%iec,j,1) * netH(G%isc:G%iec) * G%H_to_m

  ! Add in the SW heating for purposes of calculating the net
  ! surface buoyancy flux affecting the top layer.
  !netHeat(:) = netHeatMinusSW(:) + sum( penSWbnd(:,:), dim=1 )
  netHeat(G%isc:G%iec) = netHeatMinusSW(G%isc:G%iec) + netPen(G%isc:G%iec,1)

  ! Convert to a buoyancy flux, excluding penetrating SW heating
  buoyancyFlux(G%isc:G%iec,1) = - GoRho * ( dRhodS(G%isc:G%iec) * netSalt(G%isc:G%iec) + &
                                             dRhodT(G%isc:G%iec) * netHeat(G%isc:G%iec) ) * G%H_to_m ! m^2/s^3
  ! We also have a penetrative buoyancy flux associated with penetrative SW
  do k=2, G%ke+1
    buoyancyFlux(G%isc:G%iec,k) = - GoRho * ( dRhodT(G%isc:G%iec) * netPen(G%isc:G%iec,k) ) * G%H_to_m ! m^2/s^3
  enddo

end subroutine calculateBuoyancyFlux1d


!> 2d wrapper to compute surface buoyancy fluxes. 
subroutine calculateBuoyancyFlux2d(G, fluxes, optics, h, Temp, Salt, tv, &
                                   buoyancyFlux, netHeatMinusSW, netSalt)

! This subtourine calculates the surface buoyancy flux by adding up the heat, 
! FW and salt fluxes and linearizing about the surface state.
! This routine is a wrapper for calculateBuoyancyFlux1d.

  type(ocean_grid_type),                       intent(in)    :: G              ! ocean grid
  type(forcing),                               intent(inout) :: fluxes         ! surface fluxes
  type(optics_type),                           pointer       :: optics         ! SW ocean optics
  real, dimension(NIMEM_,NJMEM_,NKMEM_),       intent(in)    :: h              ! layer thickness (H)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),       intent(in)    :: Temp           ! temperature (deg C)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),       intent(in)    :: Salt           ! salinity (ppt)
  type(thermo_var_ptrs),                       intent(inout) :: tv             ! Thermodynamics type
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_),intent(inout) :: buoyancyFlux   ! buoy flux (m^2/s^3)
  real, dimension(NIMEM_,NJMEM_),optional,     intent(inout) :: netHeatMinusSW ! surf temp flux (K H)
  real, dimension(NIMEM_,NJMEM_),optional,     intent(inout) :: netSalt        ! surf salt flux (ppt H)

  real, dimension( SZI_(G) ) :: netT ! net temperature flux (K m/s)
  real, dimension( SZI_(G) ) :: netS ! net saln flux (ppt m/s)
  integer :: j

  netT(G%isc:G%iec) = 0. ; netS(G%isc:G%iec) = 0.

!$OMP parallel do default(none) shared(G,fluxes,optics,h,Temp,Salt,tv,buoyancyFlux,&
!$OMP                                  netHeatMinusSW,netSalt)                     &
!$OMP                     firstprivate(netT,netS)
  do j = G%jsc, G%jec
    call calculateBuoyancyFlux1d(G, fluxes, optics, h, Temp, Salt, tv, j, buoyancyFlux(:,j,:), netT, netS )
    if (present(netHeatMinusSW)) netHeatMinusSW(G%isc:G%iec,j) = netT(G%isc:G%iec)
    if (present(netSalt)) netSalt(G%isc:G%iec,j) = netS(G%isc:G%iec)
  enddo ! j

end subroutine calculateBuoyancyFlux2d


!> Write out chksums for basic state variables.
subroutine MOM_forcing_chksum(mesg, fluxes, G, haloshift)


  character(len=*),                    intent(in) :: mesg
  type(forcing),                       intent(in) :: fluxes
  type(ocean_grid_type),               intent(in) :: G
  integer, optional,                   intent(in) :: haloshift

! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m.
!  (in)      uh - Volume flux through zonal faces = u*h*dy, m3 s-1.
!  (in)      vh - Volume flux through meridional faces = v*h*dx, in m3 s-1.
!  (in)      G - The ocean grid structure.
  integer :: is, ie, js, je, nz, hshift
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  hshift=1; if (present(haloshift)) hshift=haloshift

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  if (associated(fluxes%taux)) &
    call uchksum(fluxes%taux, mesg//" fluxes%taux",G,haloshift=1)
  if (associated(fluxes%tauy)) &
    call vchksum(fluxes%tauy, mesg//" fluxes%tauy",G,haloshift=1)
  if (associated(fluxes%ustar)) &
    call hchksum(fluxes%ustar, mesg//" fluxes%ustar",G,haloshift=1)
  if (associated(fluxes%buoy)) &
    call hchksum(fluxes%buoy, mesg//" fluxes%buoy ",G,haloshift=hshift)
  if (associated(fluxes%sw)) &
    call hchksum(fluxes%sw, mesg//" fluxes%sw",G,haloshift=hshift)
  if (associated(fluxes%sw_vis_dir)) &
    call hchksum(fluxes%sw_vis_dir, mesg//" fluxes%sw_vis_dir",G,haloshift=hshift)
  if (associated(fluxes%sw_vis_dif)) &
    call hchksum(fluxes%sw_vis_dif, mesg//" fluxes%sw_vis_dif",G,haloshift=hshift)
  if (associated(fluxes%sw_nir_dir)) &
    call hchksum(fluxes%sw_nir_dir, mesg//" fluxes%sw_nir_dir",G,haloshift=hshift)
  if (associated(fluxes%sw_nir_dif)) &
    call hchksum(fluxes%sw_nir_dif, mesg//" fluxes%sw_nir_dif",G,haloshift=hshift)
  if (associated(fluxes%lw)) &
    call hchksum(fluxes%lw, mesg//" fluxes%lw",G,haloshift=hshift)
  if (associated(fluxes%latent)) &
    call hchksum(fluxes%latent, mesg//" fluxes%latent",G,haloshift=hshift)
  if (associated(fluxes%latent_evap_diag)) &
    call hchksum(fluxes%latent_evap_diag, mesg//" fluxes%latent_evap_diag",G,haloshift=hshift)
  if (associated(fluxes%latent_fprec_diag)) &
    call hchksum(fluxes%latent_fprec_diag, mesg//" fluxes%latent_fprec_diag",G,haloshift=hshift)
  if (associated(fluxes%latent_frunoff_diag)) &
    call hchksum(fluxes%latent_frunoff_diag, mesg//" fluxes%latent_frunoff_diag",G,haloshift=hshift)
  if (associated(fluxes%sens)) &
    call hchksum(fluxes%sens, mesg//" fluxes%sens",G,haloshift=hshift)
  if (associated(fluxes%evap)) &
    call hchksum(fluxes%evap, mesg//" fluxes%evap",G,haloshift=hshift)
  if (associated(fluxes%lprec)) &
    call hchksum(fluxes%lprec, mesg//" fluxes%lprec",G,haloshift=hshift)
  if (associated(fluxes%fprec)) &
    call hchksum(fluxes%fprec, mesg//" fluxes%fprec",G,haloshift=hshift)
  if (associated(fluxes%vprec)) &
    call hchksum(fluxes%vprec, mesg//" fluxes%vprec",G,haloshift=hshift)
  if (associated(fluxes%seaice_melt)) &
    call hchksum(fluxes%seaice_melt, mesg//" fluxes%seaice_melt",G,haloshift=hshift)
  if (associated(fluxes%p_surf)) &
    call hchksum(fluxes%p_surf, mesg//" fluxes%p_surf",G,haloshift=hshift)
  if (associated(fluxes%salt_flux)) &
    call hchksum(fluxes%salt_flux, mesg//" fluxes%salt_flux",G,haloshift=hshift)
  if (associated(fluxes%TKE_tidal)) &
    call hchksum(fluxes%TKE_tidal, mesg//" fluxes%TKE_tidal",G,haloshift=hshift)
  if (associated(fluxes%ustar_tidal)) &
    call hchksum(fluxes%ustar_tidal, mesg//" fluxes%ustar_tidal",G,haloshift=hshift)
  if (associated(fluxes%lrunoff)) &
    call hchksum(fluxes%lrunoff, mesg//" fluxes%lrunoff",G,haloshift=hshift)
  if (associated(fluxes%frunoff)) &
    call hchksum(fluxes%frunoff, mesg//" fluxes%frunoff",G,haloshift=hshift)
  if (associated(fluxes%heat_content_lrunoff)) &
    call hchksum(fluxes%heat_content_lrunoff, mesg//" fluxes%heat_content_lrunoff",G,haloshift=hshift)
  if (associated(fluxes%heat_content_frunoff)) &
    call hchksum(fluxes%heat_content_frunoff, mesg//" fluxes%heat_content_frunoff",G,haloshift=hshift)
  if (associated(fluxes%heat_content_lprec)) &
    call hchksum(fluxes%heat_content_lprec, mesg//" fluxes%heat_content_lprec",G,haloshift=hshift)
  if (associated(fluxes%heat_content_fprec)) &
    call hchksum(fluxes%heat_content_fprec, mesg//" fluxes%heat_content_fprec",G,haloshift=hshift)
  if (associated(fluxes%heat_content_cond)) &
    call hchksum(fluxes%heat_content_cond, mesg//" fluxes%heat_content_cond",G,haloshift=hshift)
  if (associated(fluxes%heat_content_massout)) &
    call hchksum(fluxes%heat_content_massout, mesg//" fluxes%heat_content_massout",G,haloshift=hshift)
end subroutine MOM_forcing_chksum


!> Write out values of the fluxes arrays at the i,j location
subroutine forcing_SinglePointPrint(fluxes, G, i, j, mesg)

  type(forcing),                       intent(in) :: fluxes
  type(ocean_grid_type),               intent(in) :: G
  character(len=*),                    intent(in) :: mesg
  integer,                             intent(in) :: i, j

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

  subroutine locMsg(array,aname)
  real, dimension(:,:), pointer :: array
  character(len=*)              :: aname

  if (associated(array)) then
    write(0,'(3a,es15.3)') 'MOM_forcing_type, forcing_SinglePointPrint: ',trim(aname),' = ',array(i,j)
  else
    write(0,'(4a)') 'MOM_forcing_type, forcing_SinglePointPrint: ',trim(aname),' is not associated.'
  endif

  end subroutine locMsg

end subroutine forcing_SinglePointPrint


!> Register members of the forcing type for diagnostics
subroutine register_forcing_type_diags(Time, diag, use_temperature, handles)
  type(time_type),     intent(in)    :: Time
  type(diag_ctrl),     intent(inout) :: diag
  logical,             intent(in)    :: use_temperature !< True if T/S are in use
  type(forcing_diags), intent(inout) :: handles

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
  handles%id_seaice_melt = register_diag_field('ocean_model', 'seaice_melt',       &
     diag%axesT1, Time, 'water flux to ocean from sea ice melt(> 0) or form(< 0)', &
     'kilogram/(meter^2 * second)',                                                &
      standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics',     &
      cmor_field_name='fsitherm', cmor_units='kg m-2 s-1',                         &
      cmor_standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics',&
      cmor_long_name='water flux to ocean from sea ice melt(> 0) or form(< 0)')

  handles%id_precip = register_diag_field('ocean_model', 'precip', diag%axesT1, Time, &
        'Liquid + frozen precipitation into ocean', 'kilogram/(meter^2 * second)')

  handles%id_fprec = register_diag_field('ocean_model', 'fprec', diag%axesT1, Time,     &
        'Frozen precipitation into ocean', 'kilogram meter-2 second-1',                 &
        standard_name='snowfall_flux', cmor_field_name='prsn', cmor_units='kg m-2 s-1', &
        cmor_standard_name='snowfall_flux', cmor_long_name='Snowfall Flux where Ice Free Ocean over Sea')

  handles%id_lprec = register_diag_field('ocean_model', 'lprec', diag%axesT1, Time,       &
        'Liquid precipitation into ocean', 'kilogram/(meter^2 * second)',                 &
        standard_name='rainfall_flux',                                                    &
        cmor_field_name='pr', cmor_units='kg m-2 s-1', cmor_standard_name='rainfall_flux',&
        cmor_long_name='Rainfall Flux where Ice Free Ocean over Sea')

  handles%id_vprec = register_diag_field('ocean_model', 'vprec', diag%axesT1, Time, &
        'Virtual liquid precip into ocean due to SSS restoring', 'kilogram/(meter^2 second)')

  handles%id_frunoff = register_diag_field('ocean_model', 'frunoff', diag%axesT1, Time,    &
        'Frozen runoff (calving) and iceberg melt into ocean', 'kilogram/(meter^2 second)',&
        standard_name='water_flux_into_sea_water_from_icebergs',                           &
        cmor_field_name='ficeberg', cmor_units='kg m-2 s-1',                               &
        cmor_standard_name='water_flux_into_sea_water_from_icebergs',                      &
        cmor_long_name='Water Flux into seawater from icebergs')

  handles%id_lrunoff = register_diag_field('ocean_model', 'lrunoff', diag%axesT1, Time, &
        'Liquid runoff (rivers) into ocean', 'kilogram meter-2 second-1',                     &
        standard_name='water_flux_into_sea_water_from_rivers', cmor_field_name='friver',      &
        cmor_units='kg m-2 s-1', cmor_standard_name='water_flux_into_sea_water_from_rivers',  &
        cmor_long_name='Water Flux into Sea Water From Rivers')

  handles%id_net_massout = register_diag_field('ocean_model', 'net_massout', diag%axesT1, Time, &
        'Net mass leaving the ocean due to evaporation, seaice formation', 'kilogram meter-2 second-1')

  handles%id_net_massin  = register_diag_field('ocean_model', 'net_massin', diag%axesT1, Time, &
        'Net mass entering ocean due to precip, runoff, ice melt', 'kilogram meter-2 second-1')


  !=========================================================================
  ! area integrated surface mass transport 

  handles%id_total_prcme = register_scalar_field('ocean_model', 'total_PRCmE', Time, diag,         &
      long_name='Area integrated net surface water flux (precip+melt+liq runoff+ice calving-evap)',&
      units='kg/s', standard_name='water_flux_into_sea_water_area_integrated')

  handles%id_total_evap = register_scalar_field('ocean_model', 'total_evap', Time, diag,&
      long_name='Area integrated evap/condense at ocean surface',                       &
      units='kg/s', standard_name='water_evaporation_flux_area_integrated')

  handles%id_total_seaice_melt = register_scalar_field('ocean_model', 'total_seaice_melt', Time, diag, &
      long_name='Area integrated sea ice melt (>0) or form (<0)', units='kg/s',                        &
      standard_name='water_flux_into_sea_water_due_to_sea_ice_thermodynamics_area_integrated')

  handles%id_total_precip = register_scalar_field('ocean_model', 'total_precip', Time, diag, &
      long_name='Area integrated liquid+frozen precip into ocean', units='kg/s')

  handles%id_total_fprec = register_scalar_field('ocean_model', 'total_fprec', Time, diag,&
      long_name='Area integrated frozen precip into ocean', units='kg/s',                 &
      standard_name='snowfall_flux_area_integrated')

  handles%id_total_lprec = register_scalar_field('ocean_model', 'total_lprec', Time, diag,&
      long_name='Area integrated liquid precip into ocean', units='kg/s',                 &
      standard_name='rainfall_flux_area_integrated')

  handles%id_total_vprec = register_scalar_field('ocean_model', 'total_vprec', Time, diag, &
      long_name='Area integrated virtual liquid precip due to SSS restoring)', units='kg/s')

  handles%id_total_frunoff = register_scalar_field('ocean_model', 'total_frunoff', Time, diag, &
      long_name='Area integrated frozen runoff (calving) & iceberg melt into ocean', units='kg/s')

  handles%id_total_lrunoff = register_scalar_field('ocean_model', 'total_lrunoff', Time, diag, &
      long_name='Area integrated liquid runoff into ocean', units='kg/s')

  handles%id_total_net_massout = register_scalar_field('ocean_model', 'total_net_massout', Time, diag, &
      long_name='Area integrated mass leaving ocean due to evap and seaice form', units='kg/s')

  handles%id_total_net_massin = register_scalar_field('ocean_model', 'total_net_massin', Time, diag, &
      long_name='Area integrated mass entering ocean due to predip, runoff, ice melt', units='kg/s')


  !===============================================================
  ! surface heat flux maps 

  handles%id_heat_content_frunoff = register_diag_field('ocean_model', 'heat_content_frunoff',          &
        diag%axesT1, Time, 'Heat content (relative to 0C) of solid runoff into ocean', 'Watt meter-2',  &
        standard_name='temperature_flux_due_to_solid_runoff_expressed_as_heat_flux_into_sea_water',     &
        cmor_field_name='hfsolidrunoffds', cmor_units='W m-2',                                          &
        cmor_standard_name='temperature_flux_due_to_solid_runoff_expressed_as_heat_flux_into_sea_water',&
        cmor_long_name='Temperature Flux due to Solid Runoff Expressed as Heat Flux into Sea Water')

  handles%id_heat_content_lrunoff = register_diag_field('ocean_model', 'heat_content_lrunoff',         &
        diag%axesT1, Time, 'Heat content (relative to 0C) of liquid runoff into ocean', 'Watt meter-2',&
        standard_name='temperature_flux_due_to_runoff_expressed_as_heat_flux_into_sea_water',          &
        cmor_field_name='hfrunoffds', cmor_units='W m-2',                                              &
        cmor_standard_name='temperature_flux_due_to_runoff_expressed_as_heat_flux_into_sea_water',     &
        cmor_long_name='Temperature Flux due to Runoff Expressed as Heat Flux into Sea Water')

  handles%id_heat_content_lprec = register_diag_field('ocean_model', 'heat_content_lprec',             &
        diag%axesT1,Time,'Heat content (relative to 0degC) of liquid precip entering ocean',           &
        'W/m^2',standard_name='temperature_flux_due_to_rainfall_expressed_as_heat_flux_into_sea_water',&
        cmor_field_name='hfrainds', cmor_units='W m-2',                                                &
        cmor_standard_name='temperature_flux_due_to_rainfall_expressed_as_heat_flux_into_sea_water',   &
        cmor_long_name='Heat Flux into Sea Water due to Liquid Precipitation')

  handles%id_heat_content_fprec = register_diag_field('ocean_model', 'heat_content_fprec',&
        diag%axesT1,Time,'Heat content (relative to 0degC) of frozen prec entering ocean',&
        'Watt/m^2')

  handles%id_heat_content_vprec = register_diag_field('ocean_model', 'heat_content_vprec',   &
        diag%axesT1,Time,'Heat content (relative to 0degC) of virtual precip entering ocean',&
        'Watt/m^2')

  handles%id_heat_content_cond = register_diag_field('ocean_model', 'heat_content_cond',   &
        diag%axesT1,Time,'Heat content (relative to 0degC) of water condensing into ocean',&
        'Watt/m^2')

  handles%id_heat_content_surfwater = register_diag_field('ocean_model', 'heat_content_surfwater',&
         diag%axesT1, Time,                                                                       &
        'Heat content (relative to 0degC) of net water crossing ocean surface (frozen+liquid)',   &
        'Watt/m^2')

  handles%id_heat_content_massout = register_diag_field('ocean_model', 'heat_content_massout',           &
         diag%axesT1, Time,'Heat content (relative to 0degC)of net mass leaving ocean ocean',            &
        'Watt/m^2',                                                                                      &
        cmor_field_name='hfevapds', cmor_units='W m-2',                                                  &
        cmor_standard_name='temperature_flux_due_to_evaporation_expressed_as_heat_flux_out_of_sea_water',&
        cmor_long_name='Heat flux out of sea water due to evaporating water')


  handles%id_heat_content_massin = register_diag_field('ocean_model', 'heat_content_massin',&
         diag%axesT1, Time,'Heat content (relative to 0degC)of net mass entering ocean ocean', &
        'Watt/m^2')

  handles%id_net_heat_coupler = register_diag_field('ocean_model', 'net_heat_coupler',          &
        diag%axesT1,Time,'Surface ocean heat flux from SW+LW+latent+sensible (via the coupler)',&
        'Watt/m^2')

  handles%id_net_heat_surface = register_diag_field('ocean_model', 'net_heat_surface',diag%axesT1,  &
        Time,'Surface ocean heat flux from SW+LW+lat+sens+mass transfer+frazil+restore', 'Watt/m^2',&
        standard_name='surface_downward_heat_flux_in_sea_water', cmor_field_name='hfds',            &
        cmor_units='W m-2', cmor_standard_name='surface_downard_heat_flux_in_sea_water',            &
        cmor_long_name='Surface ocean heat flux from SW+LW+latent+sensible+mass transfer+frazil')

  handles%id_sw = register_diag_field('ocean_model', 'SW', diag%axesT1, Time,                              &
        'Shortwave radiation flux into ocean', 'Watt meter-2',                                             &
        standard_name='surface_net_downward_shortwave_flux', cmor_field_name='rsntds', cmor_units='W m-2', &
        cmor_standard_name='net_downward_shortwave_flux_at_sea_water_surface',                             &
        cmor_long_name='Net Downward Shortwave Radiation at Sea Water Surface')

  handles%id_LwLatSens = register_diag_field('ocean_model', 'LwLatSens', diag%axesT1, Time, &
        'Combined longwave, latent, and sensible heating at ocean surface', 'Watt/m^2')

  handles%id_lw = register_diag_field('ocean_model', 'LW', diag%axesT1, Time,                            &
        'Longwave radiation flux into ocean', 'Watt meter-2',                                            &
        standard_name='surface_net_downward_longwave_flux', cmor_field_name='rlds', cmor_units='W m-2',  &
        cmor_standard_name='surface_net_downward_longwave_flux',                                         &
        cmor_long_name='Surface Net Downward Longwave Radiation')

  handles%id_lat = register_diag_field('ocean_model', 'latent', diag%axesT1, Time,                      &
        'Latent heat flux into ocean due to fusion and evaporation (negative means ocean losses heat)', &
        'Watt meter-2', cmor_field_name='hfls', cmor_units='W m-2',                                     &
        cmor_standard_name='surface_downward_latent_heat_flux',                                         &
        cmor_long_name='Surface Downward Latent Heat Flux')

  handles%id_lat_evap = register_diag_field('ocean_model', 'latent_evap', diag%axesT1, Time, &
        'Latent heat flux into ocean due to evaporation/condensation', 'Watt/m^2')

  handles%id_lat_fprec = register_diag_field('ocean_model', 'latent_fprec_diag', diag%axesT1, Time,&
        'Latent heat flux into ocean due to melting of frozen precipitation', 'Watt meter-2',      &
        cmor_field_name='hfsnthermds', cmor_units='W m-2',                                         &
        cmor_standard_name='heat_flux_into_sea_water_due_to_snow_thermodynamics',                  &
        cmor_long_name='Latent heat to melt frozen precipitation')

  handles%id_lat_frunoff = register_diag_field('ocean_model', 'latent_frunoff', diag%axesT1, Time, &
        'Latent heat flux into ocean due to melting of icebergs', 'Watt/m^2',                      &
        cmor_field_name='hfibthermds', cmor_units='W m-2',                                         &
        cmor_standard_name='heat_flux_into_sea_water_due_to_iceberg_thermodynamics',               &
        cmor_long_name='Heat flux into sea water due to iceberg thermodynamics')

  handles%id_sens = register_diag_field('ocean_model', 'sensible', diag%axesT1, Time,&
        'Sensible heat flux into ocean', 'Watt meter-2',                             &
        standard_name='surface_downward_sensible_heat_flux',                         &
        cmor_field_name='hfss', cmor_units='W m-2',                                  &
        cmor_standard_name='surface_downward_sensible_heat_flux',                    &
        cmor_long_name='Surface downward sensible heat flux')

  handles%id_heat_restore = register_diag_field('ocean_model', 'heat_restore', diag%axesT1, Time, &
        'Restoring surface heat flux into ocean', 'Watt/m^2')


  !===============================================================
  ! area integrated surface heat fluxes  

  handles%id_total_heat_content_frunoff = register_scalar_field('ocean_model',  &
      'total_heat_content_frunoff', Time, diag,                                 &
      long_name='Area integrated heat content (relative to 0C) of solid runoff',&
      units='Watt')

  handles%id_total_heat_content_lrunoff = register_scalar_field('ocean_model',   &
      'total_heat_content_lrunoff', Time, diag,                                  &
      long_name='Area integrated heat content (relative to 0C) of liquid runoff',&
      units='Watt')

  handles%id_total_heat_content_lprec = register_scalar_field('ocean_model',     &
      'total_heat_content_lprec', Time, diag,                                    &
      long_name='Area integrated heat content (relative to 0C) of liquid precip',&
      units='Watt')

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

  handles%id_total_heat_content_massout = register_scalar_field('ocean_model',         &
      'total_heat_content_massout', Time, diag,                                        &
      long_name='Area integrated heat content (relative to 0C) of water leaving ocean',&
      units='Watt')

  handles%id_total_heat_content_massin = register_scalar_field('ocean_model',           &
      'total_heat_content_massin', Time, diag,                                          &
      long_name='Area integrated heat content (relative to 0C) of water entering ocean',&
      units='Watt')

  handles%id_total_net_heat_coupler = register_scalar_field('ocean_model',                       &
      'total_net_heat_coupler', Time, diag,                                                      &
      long_name='Area integrated surface heat flux from SW+LW+latent+sensible (via the coupler)',&
      units='Watt')

  handles%id_total_net_heat_surface = register_scalar_field('ocean_model',                  &
      'total_net_heat_surface', Time, diag,                                                 &
      long_name='Area integrated surface heat flux from SW+LW+lat+sens+mass+frazil+restore',&
      units='Watt')

  handles%id_total_sw = register_scalar_field('ocean_model',                  &
      'total_sw', Time, diag,                                                 &
      long_name='Area integrated net downward shortwave at sea water surface',&
      units='Watt')

  handles%id_total_LwLatSens = register_scalar_field('ocean_model',&
      'total_LwLatSens', Time, diag,                               &
      long_name='Area integrated longwave+latent+sensible heating',&
      units='Watt')

  handles%id_total_lw = register_scalar_field('ocean_model',                 &
      'total_lw', Time, diag,                                                &
      long_name='Area integrated net downward longwave at sea water surface',&
      units='Watt')

  handles%id_total_lat = register_scalar_field('ocean_model',       &
      'total_lat', Time, diag,                                      &
      long_name='Area integrated surface downward latent heat flux',&
      units='Watt')

  handles%id_total_lat_evap = register_scalar_field('ocean_model',      &
      'total_lat_evap', Time, diag,                                     &
      long_name='Area integrated latent heat flux due to evap/condense',&
      units='Watt')

  handles%id_total_lat_fprec = register_scalar_field('ocean_model',             &
      'total_lat_fprec', Time, diag,                                            &
      long_name='Area integrated latent heat flux due to melting frozen precip',&
      units='Watt')

  handles%id_total_lat_frunoff = register_scalar_field('ocean_model',      &
      'total_lat_frunoff', Time, diag,                                     &
      long_name='Area integrated latent heat flux due to melting icebergs',&
      units='Watt')

  handles%id_total_sens = register_scalar_field('ocean_model',&
      'total_sens', Time, diag,                               &
      long_name='Area integrated downward sensible heat flux',&
      units='Watt')

  handles%id_total_heat_restore = register_scalar_field('ocean_model',&
      'total_heat_restore', Time, diag,                               &
      long_name='Area integrated surface heat flux from restoring',   &
      units='Watt')


  !===============================================================
  ! maps of surface salt fluxes, virtual precip fluxes, and adjustments  

  handles%id_saltflux = register_diag_field('ocean_model', 'salt_flux', diag%axesT1, Time,        &
        'Salt flux into ocean at surface', 'kilogram meter-2 second-1', cmor_field_name='sfdsi',  &
        cmor_units='kg m-2 s-1', cmor_standard_name='downward_sea_ice_basal_salt_flux',           &
        cmor_long_name='Downward Sea Ice Basal Salt Flux')

  handles%id_saltFluxIn = register_diag_field('ocean_model', 'salt_flux_in', diag%axesT1, Time, &
        'Salt flux into ocean at surface from coupler', 'kilogram/(meter^2 * second)')

  handles%id_saltFluxRestore = register_diag_field('ocean_model', 'salt_flux_restore', &
        diag%axesT1,Time,'Salt flux into ocean at surface due to restoring',           &
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

  handles%id_total_saltflux = register_scalar_field('ocean_model', 'total_salt_flux', &
      Time, diag, long_name='Area integrated surface salt flux', units='kg')

  handles%id_total_saltFluxIn = register_scalar_field('ocean_model', 'total_salt_Flux_In', &
      Time, diag, long_name='Area integrated surface salt flux at surface from coupler', units='kg')

  handles%id_total_saltFluxRestore = register_scalar_field('ocean_model', 'total_salt_Flux_Restore', &
      Time, diag, long_name='Area integrated surface salt flux due to restoring', units='kg')


end subroutine register_forcing_type_diags


subroutine forcing_accumulate(flux_tmp, fluxes, dt, G, wt2)
  type(forcing),         intent(in)    :: flux_tmp
  type(forcing),         intent(inout) :: fluxes
  real,                  intent(in)    :: dt
  type(ocean_grid_type), intent(inout) :: G
  real,                  intent(out)   :: wt2
  !   This subroutine copies mechancal forcing from flux_tmp to fluxes and
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
  if (associated(fluxes%heat_restore) .and. associated(flux_tmp%heat_restore)) then
    do j=js,je ; do i=is,ie
      fluxes%heat_restore(i,j) = wt1*fluxes%heat_restore(i,j) + wt2*flux_tmp%heat_restore(i,j)
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


!> Offers mechanical forcing fields for diagnostics
subroutine mech_forcing_diags(fluxes, dt, G, diag, handles)

! This subroutine offers forcing fields for diagnostics. 
! These fields must be registered in register_forcing_type_diags.

  type(forcing),         intent(in)    :: fluxes
  real,                  intent(in)    :: dt
  type(ocean_grid_type), intent(in)    :: G
  type(diag_ctrl),       intent(in)    :: diag
  type(forcing_diags),   intent(inout) :: handles

!  fluxes  = A structure containing pointers to any possible
!            forcing fields.  Unused fields are unallocated.
!  dt      = time step 
!  G       = ocean grid structure
!  diag    = structure used to regulate diagnostic output 
!  handles = ids for diagnostic manager 

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

  endif

  call cpu_clock_end(handles%id_clock_forcing)
end subroutine mech_forcing_diags


!> Offers buoyancy forcing fields for diagnostics
subroutine forcing_diagnostics(fluxes, state, dt, G, diag, handles)

! This subroutine offers forcing fields for diagnostics. 
! These fields must be registered in register_forcing_type_diags.

  type(forcing),         intent(in)    :: fluxes
  type(surface),         intent(in)    :: state 
  real,                  intent(in)    :: dt
  type(ocean_grid_type), intent(in)    :: G
  type(diag_ctrl),       intent(in)    :: diag
  type(forcing_diags),   intent(inout) :: handles

!  fluxes  = A structure containing pointers to any possible
!            forcing fields.  Unused fields are unallocated.
!  dt      = time step 
!  G       = ocean grid structure
!  diag    = structure used to regulate diagnostic output 
!  handles = ids for diagnostic manager 

  real, dimension(SZI_(G),SZJ_(G)) :: sum
  real :: total_transport ! for diagnosing integrated boundary transport
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

    if (handles%id_prcme > 0 .or. handles%id_total_prcme > 0) then
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

    if ((handles%id_evap > 0) .and. ASSOCIATED(fluxes%evap)) &
      call post_data(handles%id_evap, fluxes%evap, diag)
    if ((handles%id_total_evap > 0) .and. ASSOCIATED(fluxes%evap)) then 
      total_transport = global_area_integral(fluxes%evap(:,:),G)   
      call post_data(handles%id_total_evap, total_transport, diag)
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

    if ((handles%id_lprec > 0) .and. ASSOCIATED(fluxes%lprec)) &
      call post_data(handles%id_lprec, fluxes%lprec, diag)
    if ((handles%id_total_lprec > 0) .and. ASSOCIATED(fluxes%lprec)) then
      total_transport = global_area_integral(fluxes%lprec(:,:),G)   
      call post_data(handles%id_total_lprec, total_transport, diag)
    endif

    if ((handles%id_fprec > 0) .and. ASSOCIATED(fluxes%fprec)) &
      call post_data(handles%id_fprec, fluxes%fprec, diag)
    if ((handles%id_total_fprec > 0) .and. ASSOCIATED(fluxes%fprec)) then
      total_transport = global_area_integral(fluxes%fprec(:,:),G)   
      call post_data(handles%id_total_fprec, total_transport, diag)
    endif

    if ((handles%id_vprec > 0) .and. ASSOCIATED(fluxes%vprec)) &
      call post_data(handles%id_vprec, fluxes%vprec, diag)
    if ((handles%id_total_vprec > 0) .and. ASSOCIATED(fluxes%vprec)) then
      total_transport = global_area_integral(fluxes%vprec(:,:),G)   
      call post_data(handles%id_total_vprec, total_transport, diag)
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


    if (handles%id_net_heat_coupler > 0 .or. handles%id_total_net_heat_coupler > 0) then
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
    endif

    if (handles%id_net_heat_surface > 0 .or. handles%id_total_net_heat_surface > 0) then
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
      if (ASSOCIATED(fluxes%heat_restore))         sum(:,:) = sum(:,:) + fluxes%heat_restore(:,:)
      call post_data(handles%id_net_heat_surface, sum, diag)

      if(handles%id_total_net_heat_surface > 0) then 
        total_transport = global_area_integral(sum,G)   
        call post_data(handles%id_total_net_heat_surface, total_transport, diag)
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

    if ((handles%id_sw > 0) .and. ASSOCIATED(fluxes%sw)) then
      call post_data(handles%id_sw, fluxes%sw, diag)
    endif 
    if ((handles%id_total_sw > 0) .and. ASSOCIATED(fluxes%sw)) then
      total_transport = global_area_integral(fluxes%sw,G)   
      call post_data(handles%id_total_sw, total_transport, diag)
    endif 

    if ((handles%id_lw > 0) .and. ASSOCIATED(fluxes%lw)) then
      call post_data(handles%id_lw, fluxes%lw, diag)
    endif 
    if ((handles%id_total_lw > 0) .and. ASSOCIATED(fluxes%lw)) then
      total_transport = global_area_integral(fluxes%lw,G)   
      call post_data(handles%id_total_lw, total_transport, diag)
    endif 

    if ((handles%id_lat > 0) .and. ASSOCIATED(fluxes%latent)) then
      call post_data(handles%id_lat, fluxes%latent, diag)
    endif 
    if ((handles%id_total_lat > 0) .and. ASSOCIATED(fluxes%latent)) then
      total_transport = global_area_integral(fluxes%latent,G)   
      call post_data(handles%id_total_lat, total_transport, diag)
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

    if ((handles%id_heat_restore > 0) .and. ASSOCIATED(fluxes%heat_restore)) then
      call post_data(handles%id_heat_restore, fluxes%heat_restore, diag)
    endif 
    if ((handles%id_total_heat_restore > 0) .and. ASSOCIATED(fluxes%heat_restore)) then
      total_transport = global_area_integral(fluxes%heat_restore,G)   
      call post_data(handles%id_total_heat_restore, total_transport, diag)
    endif 


    ! post the diagnostics for boundary salt fluxes ==========================

    if ((handles%id_saltflux > 0) .and. ASSOCIATED(fluxes%salt_flux)) &
      call post_data(handles%id_saltflux, fluxes%salt_flux, diag)
    if ((handles%id_total_saltflux > 0) .and. ASSOCIATED(fluxes%salt_flux)) then
      total_transport = ppt2mks*global_area_integral(fluxes%salt_flux,G)   
      call post_data(handles%id_total_saltflux, total_transport, diag)
    endif

    if ((handles%id_saltFluxRestore > 0) .and. ASSOCIATED(fluxes%salt_flux_restore)) &
      call post_data(handles%id_saltFluxRestore, fluxes%salt_flux_restore, diag)
    if ((handles%id_total_saltFluxRestore > 0) .and. ASSOCIATED(fluxes%salt_flux_restore)) then
      total_transport = ppt2mks*global_area_integral(fluxes%salt_flux_restore,G)   
      call post_data(handles%id_total_saltFluxRestore, total_transport, diag)
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


    ! remamining boundary terms ==================================================

    if ((handles%id_psurf > 0) .and. ASSOCIATED(fluxes%p_surf))                      &
      call post_data(handles%id_psurf, fluxes%p_surf, diag)

    if ((handles%id_TKE_tidal > 0) .and. ASSOCIATED(fluxes%TKE_tidal))               &
      call post_data(handles%id_TKE_tidal, fluxes%TKE_tidal, diag)

    if ((handles%id_buoy > 0) .and. ASSOCIATED(fluxes%buoy))                         &
      call post_data(handles%id_buoy, fluxes%buoy, diag)


  endif

  call cpu_clock_end(handles%id_clock_forcing)
end subroutine forcing_diagnostics


!> Conditionally allocates fields within the forcing type
subroutine allocate_forcing_type(G, fluxes, stress, ustar, water, heat)
  type(ocean_grid_type), intent(in) :: G !< Ocean grid structure
  type(forcing),      intent(inout) :: fluxes !< Forcing fields structure
  logical, optional,     intent(in) :: stress !< If present and true, allocate taux, tauy
  logical, optional,     intent(in) :: ustar !< If present and true, allocate ustar
  logical, optional,     intent(in) :: water !< If present and true, allocate water fluxes
  logical, optional,     intent(in) :: heat  !< If present and true, allocate heat fluxes
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
 
  contains

  subroutine myAlloc(array, is, ie, js, je, flag)
    real, dimension(:,:), pointer :: array
    integer,           intent(in) :: is, ie, js, je !< Bounds
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


!> Deallocates the forcing type
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
  if (associated(fluxes%heat_restore))         deallocate(fluxes%heat_restore)
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
  if (associated(fluxes%frac_shelf_h))         deallocate(fluxes%frac_shelf_h)
  if (associated(fluxes%frac_shelf_u))         deallocate(fluxes%frac_shelf_u)
  if (associated(fluxes%frac_shelf_v))         deallocate(fluxes%frac_shelf_v)
  if (associated(fluxes%rigidity_ice_u))       deallocate(fluxes%rigidity_ice_u)
  if (associated(fluxes%rigidity_ice_v))       deallocate(fluxes%rigidity_ice_v)
  if (associated(fluxes%tr_fluxes))            deallocate(fluxes%tr_fluxes)
end subroutine deallocate_forcing_type

end module MOM_forcing_type
