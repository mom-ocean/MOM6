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

use MOM_checksums,     only : hchksum, qchksum, uchksum, vchksum
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : time_type
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_EOS,           only : calculate_density_derivs
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_variables,     only : thermo_var_ptrs

use coupler_types_mod, only : coupler_2d_bc_type

implicit none ; private

#include <MOM_memory.h>

public extractFluxes1d
public extractFluxes2d
public MOM_forcing_chksum
public absorbRemainingSW
public calculateBuoyancyFlux1d
public calculateBuoyancyFlux2d
public forcing_SinglePointPrint
public deallocate_forcing_type

integer :: num_msg = 0
integer :: max_msg = 2


type, public :: forcing
   ! This structure contains pointers to the boundary 
   ! forcing used to drive the liquid ocean as part of MOM.  
   ! All fluxes are positive downward (into the ocean).  
   ! Pointers to unused fluxes should be set to NULL.

  real, pointer, dimension(:,:) :: &

    ! surface stress components and turbulent velocity scale  
    taux          => NULL(), & ! zonal wind stress (Pa)
    tauy          => NULL(), & ! meridional wind stress (Pa)
    ustar         => NULL(), & ! surface friction velocity (m/s)

    ! surface buoyancy force 
    buoy          => NULL(), & ! buoyancy flux (m2/s3)

    ! radiative heat fluxes into the ocean (W/m2) 
    sw            => NULL(), & ! shortwave heat flux (W/m2)
    sw_vis_dir    => NULL(), & ! visible, direct shortwave heat flux (W/m2)
    sw_vis_dif    => NULL(), & ! visible, diffuse shortwave heat flux (W/m2)
    sw_nir_dir    => NULL(), & ! near-IR, direct shortwave heat flux (W/m2)
    sw_nir_dif    => NULL(), & ! near-IR, diffuse shortwave heat flux (W/m2)
    lw            => NULL(), & ! longwave heat flux (W/m2) (typically negative)

    ! turbulent heat fluxes into the ocean (W/m2) 
    latent        => NULL(), & ! latent heat flux (W/m2) (typically negative)
    sens          => NULL(), & ! sensible heat flux (W/m2) (typically negative)
    heat_restore  => NULL(), & ! heat flux from SST restoring (W/m2) in idealized simulations

    ! sensible heat associated with runoff and calving 
    runoff_hflx   => NULL(), & ! heat flux associated with liq_runoff (W/m2)
    calving_hflx  => NULL(), & ! heat flux associated with froz_runoff (W/m2)

    ! water mass fluxes into the ocean ( kg/(m2 s) )
    ! these mass fluxes impact the ocean mass
    evap          => NULL(), & ! (-1)*fresh water flux evaporated out of the ocean ( kg/(m2 s) )
    liq_precip    => NULL(), & ! precipitating liquid water into the ocean ( kg/(m2 s) )
    froz_precip   => NULL(), & ! frozen water into the ocean ( kg/(m2 s) )
    virt_precip   => NULL(), & ! virtual water associated w/ SSS restoring ( kg/(m2 s) )
    liq_runoff    => NULL(), & ! liquid river runoff  ( kg/(m2 s) )
    froz_runoff   => NULL(), & ! calving land ice ( kg/(m2 s) )

    ! salt mass flux (does not contribute to ocean mass)
     salt_flux    => NULL(), & ! net salt flux into the ocean ( kg salt/(m2 s) )

    ! applied surface pressure from other component models (e.g., atmos, sea ice, land ice)
    p_surf_full   => NULL(), & ! pressure at the top ocean interface (Pa).
                               ! if there is sea-ice, then p_surf_flux is at ice-ocean interface
    p_surf        => NULL(), & ! pressure at top ocean interface (Pa) as used to drive the ocean model. 
                               ! if p_surf is limited, then p_surf may be smaller than p_surf_full,
                               ! otherwise they are the same.

    ! tide related inputs 
    TKE_tidal     => NULL(), & ! tidal energy source driving mixing in bottom boundary layer (W/m2)
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

    ! heat capacity 
    real :: C_p                ! heat capacity of seawater ( J/(K kg) )
                               ! C_p is is the same value as in thermovar_ptrs_type.

    ! passive tracer surface fluxes 
    type(coupler_2d_bc_type), pointer :: tr_fluxes  => NULL()
       ! This structure may contain an array of named fields used for passive tracer fluxes.
       ! All arrays in tr_fluxes use the coupler indexing, which has no halos. This is not 
       ! a convenient convention, but imposed on MOM6 by the coupler.  

end type forcing


type, public :: optics_type
  ! Type for ocean optical properties 

  integer :: nbands    ! number of penetrating bands of SW radiation

  real, pointer, dimension(:,:,:,:) :: &
    opacity_band => NULL()  ! SW optical depth per unit thickness (1/m)
                            ! Number of radiation bands is most rapidly varying (first) index.

  real, pointer, dimension(:,:,:) :: &
    SW_pen_band  => NULL()  ! shortwave radiation (W/m2) at the surface in each of
                            ! the nbands bands that penetrates beyond the surface.
                            ! The most rapidly varying dimension is the band.

  real, pointer, dimension(:) ::   &
    min_wavelength_band => NULL(), & ! The range of wavelengths in each band of
    max_wavelength_band => NULL()    ! penetrating shortwave radiation (nm)

end type optics_type



contains

subroutine extractFluxes1d(G, fluxes, optics, nsw, j, dt,                               &
                  DepthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
                  h, T, Net_H, Net_heat, Net_salt, Pen_SW_bnd, tv)

!  This subroutine extracts fluxes from the surface fluxes type. It works on a j-row
!  for optimization purposes.  The 2d i,j wrapper is the next subroutine below. 

  type(ocean_grid_type),          intent(in)    :: G
  type(forcing),                  intent(in)    :: fluxes
  type(optics_type),              pointer       :: optics
  integer,                        intent(in)    :: nsw
  integer,                        intent(in)    :: j
  real,                           intent(in)    :: dt
  real,                           intent(in)    :: DepthBeforeScalingFluxes
  logical,                        intent(in)    :: useRiverHeatContent
  logical,                        intent(in)    :: useCalvingHeatContent
  real, dimension(NIMEM_,NKMEM_), intent(in)    :: h
  real, dimension(NIMEM_,NKMEM_), intent(in)    :: T
  real, dimension(NIMEM_),        intent(out)   :: Net_H
  real, dimension(NIMEM_),        intent(out)   :: Net_heat
  real, dimension(NIMEM_),        intent(out)   :: Net_salt
  real, dimension(:,:),           intent(out)   :: Pen_SW_bnd
  type(thermo_var_ptrs),          intent(inout) :: tv

!  (in)      G                        = ocean grid structure
!  (in)      fluxes                   = structure containing pointers to possible
!                                       forcing fields.  Unused fields have NULL ptrs.
!  (in)      nsw                      = number of bands of penetrating shortwave radiation
!  (in)      j                        = j-index to work on
!  (in)      dt                       = time step in seconds 
!  (in)      DepthBeforeScalingFluxes =  minimum ocean thickness to allow before scaling away fluxes
!  (in)      useRiverHeatContent      =  if true, apply river heat content
!  (in)      useCalvingHeatContent    =  if true, apply calving heat content
!  (in)      h                        =  layer thickness, in m for Bouss or (kg/m2) for non-Bouss
!  (in)      T                        =  layer temperatures, in deg C

!  (out)     Net_H      = net mass flux (if non-Boussinesq) or volume flux (if Boussinesq)
!                         i.e. fresh water flux (P+R-E) into ocean over a time step (units of H)
!  (out)     Net_heat   = net heating at the surface over a time step,
!                         exclusive of heating that appears in Pen_SW (K * H)
!  (out)     Net_salt   = surface salt flux into the ocean over a time step (psu * H)
!  (out)     Pen_SW_bnd = penetrating shortwave heating at the sea surface
!                         in each penetrating band, in K H, size nsw x NIMEM_.
!  (inout)   tv         = structure containing pointers to any available
!                         thermodynamic fields. Here it is used to keep track of the
!                         heat flux associated with net mass fluxes into the ocean.

  real :: htot(SZI_(G))        ! total ocean depth (m for Bouss or kg/m2 for non-Bouss)
  real :: Pen_sw_tot(SZI_(G))  ! sum across all bands of Pen_SW (K * H)
  real :: Ih_limit             ! inverse of total depth at which surface fluxes start to be limited (1/H)
  real :: scale                ! scale scales away fluxes if total depth is thinner than DepthBeforeScalingFluxes
  real :: Irho_cp              !  1.0 / (rho_0 * C_p).
  real :: Irho0                !  1.0 / Rho0
  real :: I_Cp                 !  1.0 / C_p

  character(len=200) :: mesg
  integer            :: is, ie, nz, i, k, n
  Ih_limit    = 1.0 / (DepthBeforeScalingFluxes * G%m_to_H)
  Irho0       = 1.0 / G%Rho0 
  I_Cp        = 1.0 / fluxes%C_p
  Irho_cp     = 1.0 / (G%H_to_kg_m2 * fluxes%C_p)

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
    "MOM_forcing_type extractFluxes1d: No evaporation defined in mixedlayer.")

  if (.not.ASSOCIATED(fluxes%virt_precip)) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: fluxes%virt_precip not defined in mixedlayer.")

  if ((.not.ASSOCIATED(fluxes%liq_precip)) .or. &
      (.not.ASSOCIATED(fluxes%froz_precip))) call MOM_error(FATAL, &
    "MOM_forcing_type extractFluxes1d: No precipitation defined in mixedlayer.")

  if (useRiverHeatContent .and. &
      .not.ASSOCIATED(fluxes%runoff_hflx)) call MOM_error(FATAL, &
        "MOM_forcing_type extractFluxes1d: fluxes%runoff_hflx must be "//&
        "assocated if USE_RIVER_HEAT_CONTENT is true.")

  if (useCalvingHeatContent .and. &
      .not.ASSOCIATED(fluxes%calving_hflx)) call MOM_error(FATAL, &
        "MOM_forcing_type extractFluxes1d: fluxes%calving_hflx must be "//&
        "assocated if USE_CALVING_HEAT_CONTENT is true.")

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
        Pen_SW_bnd(n,i) = Irho_cp*scale*dt * max(0.0, optics%sw_pen_band(n,i,j))
        Pen_sw_tot(i)   = Pen_sw_tot(i) + Pen_SW_bnd(n,i)
      enddo
    else
      Pen_SW_bnd(1,i) = 0.0
    endif

    ! Volume/mass fluxes
    Net_H(i) = dt * (scale * ((((( fluxes%liq_precip(i,j)    &
                                 + fluxes%froz_precip(i,j) ) &
                                 + fluxes%evap(i,j)        ) &
                                 + fluxes%liq_runoff(i,j)  ) &
                                 + fluxes%virt_precip(i,j) ) &
                                 + fluxes%froz_runoff(i,j) ) )

    Net_H(i) = G%kg_m2_to_H * Net_H(i)

    ! for non-Bouss, we add salt mass to total ocean mass
    ! smg: is there salt mass loss from ice model?  
    if (.not.G%Boussinesq .and. ASSOCIATED(fluxes%salt_flux)) &
      Net_H(i) = Net_H(i) + (dt * G%kg_m2_to_H) * (scale * fluxes%salt_flux(i,j))

    ! Heat fluxes
    Net_heat(i) = scale * dt * Irho_cp * ( fluxes%sw(i,j) +  ((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)) )

    if (ASSOCIATED(fluxes%heat_restore)) Net_heat(i) = Net_heat(i) + (scale * (dt * Irho_cp)) * fluxes%heat_restore(i,j)

    if (useRiverHeatContent) then
      ! remove liq_runoff*SST here, to counteract its addition elsewhere
      ! smg: this is awkward; needs to be cleaned up.  
      Net_heat(i) = (Net_heat(i) + (scale*(dt*Irho_cp)) * fluxes%runoff_hflx(i,j)) - &
                     (G%kg_m2_to_H * (scale * dt)) * fluxes%liq_runoff(i,j) * T(i,1)
      if (ASSOCIATED(tv%TempxPmE)) then
        tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + (scale * dt) * &
            (I_Cp*fluxes%runoff_hflx(i,j) - fluxes%liq_runoff(i,j)*T(i,1))
      endif    
    endif

    if (useCalvingHeatContent) then
      ! remove liq_runoff*SST here, to counteract its addition elsewhere
      ! smg: this is awkward; needs to be cleaned up.  
      Net_heat(i) = Net_heat(i) + (scale*(dt*Irho_cp)) * fluxes%calving_hflx(i,j) - &
                    (G%kg_m2_to_H * (scale * dt)) * fluxes%froz_runoff(i,j) * T(i,1)
      if (ASSOCIATED(tv%TempxPmE)) then
        tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + (scale * dt) * &
            (I_Cp*fluxes%calving_hflx(i,j) - fluxes%froz_runoff(i,j)*T(i,1))
      endif    
    endif

    if (num_msg < max_msg) then
      if (Pen_SW_tot(i) > 1.000001*Irho_cp*scale*dt*fluxes%sw(i,j)) then
        num_msg = num_msg + 1
        write(mesg,'("Penetrating shortwave of ",1pe17.10, &
                    &" exceeds total shortwave of ",1pe17.10,&
                    &" at ",1pg11.4,"E, "1pg11.4,"N.")') &
               Pen_SW_tot(i),Irho_cp*scale*dt*fluxes%sw(i,j),&
               G%geoLonT(i,j),G%geoLatT(i,j)
        call MOM_error(WARNING,mesg) 
      endif
    endif

    Net_heat(i) = Net_heat(i) - Pen_SW_tot(i)

    ! Salt fluxes
    Net_salt(i) = 0.0
    ! Convert salt_flux from kg (salt) m-2 s-1 to PSU m s-1.
    if (ASSOCIATED(fluxes%salt_flux)) &
      Net_salt(i) = (scale * dt * (1000.0 * fluxes%salt_flux(i,j))) * G%kg_m2_to_H

  enddo ! i-loop

end subroutine extractFluxes1d


subroutine extractFluxes2d(G, fluxes, optics, nsw, dt,                                  &
                  DepthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
                  h, T, Net_H, Net_heat, Net_salt, Pen_SW_bnd, tv)

!  This subroutine extracts fluxes from the surface fluxes type. It is a 
!  wrapper for the 1d routine extractFluxes1d.

  type(ocean_grid_type),                 intent(in)    :: G
  type(forcing),                         intent(in)    :: fluxes
  type(optics_type),                     pointer       :: optics
  integer,                               intent(in)    :: nsw
  real,                                  intent(in)    :: dt
  real,                                  intent(in)    :: DepthBeforeScalingFluxes
  logical,                               intent(in)    :: useRiverHeatContent
  logical,                               intent(in)    :: useCalvingHeatContent
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: T
  real, dimension(NIMEM_,NJMEM_),        intent(out)   :: Net_H
  real, dimension(NIMEM_,NJMEM_),        intent(out)   :: Net_heat
  real, dimension(NIMEM_,NJMEM_),        intent(out)   :: Net_salt
  real, dimension(:,:,:),                intent(out)   :: Pen_SW_bnd
  type(thermo_var_ptrs),                 intent(inout) :: tv

!  (in)      G                        = ocean grid structure
!  (in)      fluxes                   = structure containing pointers to possible
!                                       forcing fields.  Unused fields have NULL ptrs.
!  (in)      nsw                      = number of bands of penetrating shortwave radiation
!  (in)      dt                       = time step in seconds
!  (in)      DepthBeforeScalingFluxes =  minimum ocean thickness to allow before scaling away fluxes
!  (in)      useRiverHeatContent      =  if true, apply river heat content
!  (in)      useCalvingHeatContent    =  if true, apply calving heat content
!  (in)      h                        =  layer thickness, in m for Bouss or (kg/m2) for non-Bouss
!  (in)      T                        =  layer temperatures, in deg C

!  (out)     Net_H      = net mass flux (if non-Boussinesq) or volume flux (if Boussinesq)
!                         i.e. fresh water flux (P+R-E) into ocean over a time step (units of H)
!  (out)     Net_heat   = net heating at the surface over a time step,
!                         exclusive of heating that appears in Pen_SW (K * H)
!  (out)     Net_salt   = surface salt flux into the ocean over a time step (psu * H)
!  (out)     Pen_SW_bnd = penetrating shortwave heating at the sea surface
!                         in each penetrating band, in K H, size nsw x NIMEM_.
!  (inout)   tv         = structure containing pointers to any available
!                         thermodynamic fields. Here it is used to keep track of the
!                         heat flux associated with net mass fluxes into the ocean.

  integer :: j
  do j=G%jsc, G%jec
    call extractFluxes1d(G, fluxes, optics, nsw, j, dt,                           &
            DepthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
            h(:,j,:), T(:,j,:), Net_H(:,j), Net_heat(:,j), Net_salt(:,j), Pen_SW_bnd(:,:,j), tv)
  enddo

end subroutine extractFluxes2d


subroutine calculateBuoyancyFlux1d(G, fluxes, optics, h, Temp, Salt, tv, j, buoyancyFlux, netHeatMinusSW, netSalt )

! This subtourine calculates the surface buoyancy flux by adding up the heat, 
! FW and salt fluxes and linearizing about the surface state.

  type(ocean_grid_type),                 intent(in)    :: G              ! ocean grid
  type(forcing),                         intent(in)    :: fluxes         ! surface fluxes/forcing type
  type(optics_type),                     pointer       :: optics         ! optics for penetrating SW
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h              ! layer/level thicknesses (units of H)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: Temp           ! potential or conservative temp (degrees C)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: Salt           ! salinity (ppt)
  type(thermo_var_ptrs),                 intent(inout) :: tv             ! thermodynamics type (out needed for tv%TempxPmE ????)
  integer,                               intent(in)    :: j              ! j-index of row to work on
  real, dimension(NIMEM_,NK_INTERFACE_), intent(inout) :: buoyancyFlux   ! buoyancy flux (m2/s3)
  real, dimension(NIMEM_),               intent(inout) :: netHeatMinusSW ! heat flux excluding SW (K m/s)
  real, dimension(NIMEM_),               intent(inout) :: netSalt        ! net salt flux (ppt m/s)

  ! Local variables
  integer                                   :: nsw, start, npts, k
  integer, dimension(SZI_(G),SZK_(G))       :: ksort
  real, parameter                           :: dt = 1.    ! set to unity to return a rate from extractFluxes1d
  real, dimension( SZI_(G) )                :: netH       ! net FW flux (m/s for Bouss)
  real, dimension( SZI_(G) )                :: netHeat    ! net temp flux (K m/s)
  real, dimension( optics%nbands, SZI_(G) ) :: penSWbnd   ! SW penetration bands
  real, dimension( SZI_(G) )                :: pressure   ! pressurea the surface (Pa)
  real, dimension( SZI_(G) )                :: dRhodT     ! partial derivatives of density wrt temp
  real, dimension( SZI_(G) )                :: dRhodS     ! partial derivatives of density wrt saln
  real, dimension(SZI_(G),SZK_(G)+1)        :: netPen 

  logical :: useRiverHeatContent
  logical :: useCalvingHeatContent
  real    :: depthBeforeScalingFluxes, GoRho
  real    :: H_limit_fluxes

  nsw = optics%nbands

!  smg: what do we do when have heat fluxes from calving and river?  
  useRiverHeatContent   = .False.    
  useCalvingHeatContent = .False.  

  depthBeforeScalingFluxes = max( G%Angstrom, 1.e-30 )
  pressure(:) = 0. ! Ignore atmospheric pressure
  GoRho = G%g_Earth / G%Rho0
  start = 1 + G%isc - G%isd
  npts = 1 + G%iec - G%isc

  do k=1, G%ke
    ksort(:,k) = k
  enddo
  H_limit_fluxes = depthBeforeScalingFluxes

  ! Fetch the fresh-water, heat and salt fluxes
  ! netH is the fresh-water flux
  ! netSalt is the salt flux (typically zero except under sea-ice)
  ! netHeat is the heat flux EXCEPT the penetrating SW
  ! penSWbnd is the surface SW for each band
  call extractFluxes1d(G, fluxes, optics, nsw, j, dt,                                 &
                depthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
                h(:,j,:), Temp(:,j,:), netH, netHeatMinusSW, netSalt, penSWbnd, tv)

  ! Sum over bands and attenuate as a function of depth
  ! netPen is the netSW as a function of depth
  call sumSWoverBands(G, h(:,j,:), 0.*h(:,j,:), 0.*h(:,j,1), optics%opacity_band(:,:,j,:), nsw, j, dt, &
                      H_limit_fluxes, .true., .true., ksort, penSWbnd, netPen)

  ! Density derivatives
  call calculate_density_derivs(Temp(:,j,1), Salt(:,j,1), pressure, &
                                dRhodT, dRhodS, start, npts, tv%eqn_of_state)

  ! Adjust netSalt to reflect dilution effect of FW flux
  netSalt(:) = netSalt(:) - Salt(:,j,1) * netH * G%H_to_m

  ! Add in the SW heating for purposes of calculating the net surface buoyancy flux
  !netHeat(:) = netHeatMinusSW(:) + sum( penSWbnd(:,:), dim=1 )
  netHeat(:) = netHeatMinusSW(:) + netPen(:,1)

  ! Convert to a buoyancy flux, excluding penetrating SW heating
  buoyancyFlux(:,1) = - GoRho * ( dRhodS(:) * netSalt(:) + dRhodT(:) * netHeat(:) ) ! m2/s3
  ! We also have a penetrative buoyancy flux associated with penetrative SW
  do k=2, G%ke+1
    buoyancyFlux(:,k) = - GoRho * ( dRhodT(:) * netPen(:,k) ) ! m2/s3
  enddo

end subroutine calculateBuoyancyFlux1d



subroutine calculateBuoyancyFlux2d(G, fluxes, optics, h, Temp, Salt, tv, buoyancyFlux, netHeatMinusSW, netSalt)

! This subtourine calculates the surface buoyancy flux by adding up the heat, 
! FW and salt fluxes and linearizing about the surface state.
! This routine is a wrapper for calculateBuoyancyFlux1d.

  type(ocean_grid_type),                       intent(in)    :: G              ! ocean grid
  type(forcing),                               intent(in)    :: fluxes         ! surface fluxes/forcing type
  type(optics_type),                           pointer       :: optics         ! optics for penetrating SW
  real, dimension(NIMEM_,NJMEM_,NKMEM_),       intent(in)    :: h              ! layer/level thicknesses (units of H)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),       intent(in)    :: Temp           ! potential or conservative temp (deg C)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),       intent(in)    :: Salt           ! salinity (ppt)
  type(thermo_var_ptrs),                       intent(inout) :: tv             ! Thermodynamics type (out needed for tv%TempxPmE ????)
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_),intent(inout) :: buoyancyFlux   ! buoyancy flux (m2/s3)
  real, dimension(NIMEM_,NJMEM_),optional,     intent(inout) :: netHeatMinusSW ! temp flux excluding SW (K m/s)
  real, dimension(NIMEM_,NJMEM_),optional,     intent(inout) :: netSalt        ! net salt flux (ppt m/s)

  real, dimension( SZI_(G) ) :: netT ! net temperature flux (K m/s)
  real, dimension( SZI_(G) ) :: netS ! net saln flux (ppt m/s)
  integer :: j

  netT(:) = 0. ; netS(:) = 0.

!$OMP parallel do default(shared) firstprivate(netT,netS)
  do j = G%jsc, G%jec
    call calculateBuoyancyFlux1d(G, fluxes, optics, h, Temp, Salt, tv, j, buoyancyFlux(:,j,:), netT, netS )
    if (present(netHeatMinusSW)) netHeatMinusSW(:,j) = netT(:)
    if (present(netSalt)) netSalt(:,j) = netS(:)
  enddo ! j

end subroutine calculateBuoyancyFlux2d


subroutine absorbRemainingSW(G, h, eps, htot, opacity_band, nsw, j, dt, &
                             H_limit_fluxes, correctAbsorption, absorbAllSW, &
                             ksort, T, Ttot, Pen_SW_bnd)

! This subroutine applies shortwave heating below the mixed layer.  In
! addition, it causes all of the remaining SW radiation to be absorbed,
! provided that the total water column thickness is greater than
! H_limit_fluxes.  For thinner water columns, the heating is scaled down
! proportionately, the assumption being that the remaining heating (which is
! left in Pen_SW) should go into an (absent for now) ocean bottom sediment layer.

  type(ocean_grid_type),             intent(in)    :: G
  real, dimension(NIMEM_,NKMEM_),    intent(in)    :: h, eps
  real, dimension(NIMEM_),           intent(in)    :: htot
  real, dimension(:,:,:),            intent(in)    :: opacity_band
  integer,                           intent(in)    :: nsw
  integer,                           intent(in)    :: j
  real,                              intent(in)    :: dt
  real,                              intent(in)    :: H_limit_fluxes
  logical,                           intent(in)    :: correctAbsorption
  logical,                           intent(in)    :: absorbAllSW
  integer, dimension(NIMEM_,NKMEM_), intent(in)    :: ksort
  real, dimension(NIMEM_,NKMEM_),    intent(inout) :: T
  real, dimension(NIMEM_),           intent(inout) :: Ttot
  real, dimension(:,:),              intent(inout) :: Pen_SW_bnd

! Arguments:
!  (in)      G            = ocean grid structure
!  (in)      h            = layer thickness, in m or kg m-2. 
!                           units of h are referred to as "H" below.
!  (in)      eps          = small thickness that must remain in each layer, and
!                           which will not be subject to heating (units of H)
!  (in)      htot         = total mixed layer thickness, in H
!  (in)      opacity_band = opacity in each band of penetrating shortwave
!                           radiation (1/H). The indicies are band, i, k.
!  (in)      nsw          = number of bands of penetrating shortwave radiation
!  (in)      j            = j-index to work on
!  (in)      dt           = time step (seconds)
!  (inout)   ksort        = density-sorted k-indicies
!  (inout)   T            = layer potential temperatures (deg C)
!  (inout)   Ttot         = depth integrated mixed layer temperature (units of K H)
!  (inout)   Pen_SW_bnd   = penetrating shortwave heating in each band that
!                           hits the bottom and will be redistributed through
!                           the water column (units of K H), size nsw x NIMEM_.

  real :: h_heat(SZI_(G))              ! thickness of the water column that receives
                                       ! the remaining shortwave radiation (H units).
  real :: T_chg_above(SZI_(G),SZK_(G)) ! temperature change of all the thick layers 
                                       ! above a given layer (K). The net change in the 
                                       ! of a layer is the sum of T_chg_above from all 
                                       ! the layers below, plus any contribution from 
                                       ! absorbing radiation that hits the bottom.
  real :: T_chg(SZI_(G))               ! temperature change of thick layers due to
                                       ! the remaining shortwave radiation (K)
  real :: Pen_SW_rem(SZI_(G))          ! sum across all wavelength bands of the
                                       ! penetrating shortwave heating that hits the bottom
                                       ! and will be redistributed through the water column
                                       ! (in units of K H)
  real :: SW_trans                     ! fraction of shortwave radiation that is not
                                       ! absorbed in a layer (nondimensional)
  real :: unabsorbed                   ! fraction of the shortwave radiation that
                                       ! is not absorbed because the layers are too thin
  real :: Ih_limit                     ! inverse of the total depth at which the
                                       ! surface fluxes start to be limited (1/H)
  real :: h_min_heat                   ! minimum thickness layer that should get
                                       ! heated (H)
  real :: opt_depth                    ! optical depth of a layer (non-dim)
  real :: exp_OD                       ! exp(-opt_depth) (non-dim)
  real :: heat_bnd                     ! heating due to absorption in the current
                                       ! layer by the current band, including any piece that
                                       ! is moved upward (K H units)
  real :: SWa                          ! fraction of the absorbed shortwave that is
                                       ! moved to layers above with correctAbsorption (non-dim)
  logical :: SW_Remains                ! If true, some column has shortwave radiation that
                                       ! was not entirely absorbed

  integer :: is, ie, nz, i, k, ks, n
  SW_Remains = .false.

  h_min_heat = 2.0*G%Angstrom + G%H_subroundoff
  is = G%isc ; ie = G%iec ; nz = G%ke

  do i=is,ie ; h_heat(i) = htot(i) ; enddo

  ! Apply penetrating SW radiation to remaining parts of layers.  Excessively thin
  ! layers are not heated.
  do ks=1,nz

    do i=is,ie ; if (ksort(i,ks) > 0) then
      k = ksort(i,ks)

      T_chg_above(i,k) = 0.0

      if (h(i,k) > 1.5*eps(i,k)) then
        do n=1,nsw ; if (Pen_SW_bnd(n,i) > 0.0) then
          ! SW_trans is the SW that is transmitted THROUGH the layer
          opt_depth = h(i,k) * opacity_band(n,i,k)
          exp_OD = exp(-opt_depth)
          SW_trans = exp_OD
          ! Heating at a rate of less than 10-4 W m-2 = 10-3 K m / Century,
          ! and of the layer in question less than 1 K / Century, can be
          ! absorbed without further penetration.
          if ((nsw*Pen_SW_bnd(n,i)*SW_trans < G%m_to_H*2.5e-11*dt) .and. &
              (nsw*Pen_SW_bnd(n,i)*SW_trans < h(i,k)*dt*2.5e-8)) &
            SW_trans = 0.0

          if (correctAbsorption .and. (h_heat(i) > 0.0)) then
            if (opt_depth > 1e-5) then
              SWa = ((opt_depth + (opt_depth + 2.0)*exp_OD) - 2.0) / &
                ((opt_depth + opacity_band(n,i,k) * h_heat(i)) * &
                 (1.0 - exp_OD))
            else
              ! Use a Taylor's series expansion of the expression above for a
              ! more accurate form with very small layer optical depths.
              SWa = h(i,k) * (opt_depth * (1.0 - opt_depth)) / &
                ((h_heat(i) + h(i,k)) * (6.0 - 3.0*opt_depth))
            endif
            Heat_bnd = Pen_SW_bnd(n,i) * (1.0 - SW_trans)
            T_chg_above(i,k) = T_chg_above(i,k) + (SWa * Heat_bnd) / h_heat(i)
            T(i,k) = T(i,k) + ((1.0 - SWa) * Heat_bnd) / h(i,k)
          else
            T(i,k) = T(i,k) + Pen_SW_bnd(n,i) * (1.0 - SW_trans) / h(i,k)
          endif

          Pen_SW_bnd(n,i) = Pen_SW_bnd(n,i) * SW_trans
        endif ; enddo
      endif

      ! Add to the accumulated thickness above that could be heated.
      ! Only layers greater than h_min_heat thick should get heated.
      if (h(i,k) >= 2.0*h_min_heat) then
        h_heat(i) = h_heat(i) + h(i,k)
      elseif (h(i,k) > h_min_heat) then
        h_heat(i) = h_heat(i) + (2.0*h(i,k) - 2.0*h_min_heat)
      endif
    endif ; enddo ! i loop
  enddo ! k loop

  ! Apply heating above the layers in which it should have occurred to get the
  ! correct mean depth of the shortwave that should be absorbed by each layer.
!    if (correctAbsorption) then
!    endif

! if (.not.absorbAllSW .and. .not.correctAbsorption) return
  if (absorbAllSW) then

    ! Unless modified, there is no temperature change due to fluxes from the
    ! bottom.
    do i=is,ie ; T_chg(i) = 0.0 ; enddo

    ! If there is still shortwave radiation at this point, it could go into
    ! the bottom (with a bottom mud model), or it could be redistributed back
    ! through the water column.
    do i=is,ie
      Pen_SW_rem(i) = Pen_SW_bnd(1,i)
      do n=2,nsw ; Pen_SW_rem(i) = Pen_SW_rem(i) + Pen_SW_bnd(n,i) ; enddo
    enddo
    do i=is,ie ; if (Pen_SW_rem(i) > 0.0) SW_Remains = .true. ; enddo

    Ih_limit = 1.0 / (H_limit_fluxes * G%m_to_H)
    do i=is,ie ; if ((Pen_SW_rem(i) > 0.0) .and. (h_heat(i) > 0.0)) then
      if (h_heat(i)*Ih_limit >= 1.0) then
        T_chg(i) = Pen_SW_rem(i) / h_heat(i) ; unabsorbed = 0.0
      else
        T_chg(i) = Pen_SW_rem(i) * Ih_limit
        unabsorbed = 1.0 - h_heat(i)*Ih_limit
      endif
      do n=1,nsw ; Pen_SW_bnd(n,i) = unabsorbed * Pen_SW_bnd(n,i) ; enddo
    endif ; enddo

    do ks=nz,1,-1 ; do i=is,ie ; if (ksort(i,ks) > 0) then
      k = ksort(i,ks)
      if (T_chg(i) > 0.0) then
        ! Only layers greater than h_min_heat thick should get heated.
        if (h(i,k) >= 2.0*h_min_heat) then ; T(i,k) = T(i,k) + T_chg(i)
        elseif (h(i,k) > h_min_heat) then
          T(i,k) = T(i,k) + T_chg(i) * (2.0 - 2.0*h_min_heat/h(i,k))
        endif
      endif
      ! Increase the heating for layers above.
      T_chg(i) = T_chg(i) + T_chg_above(i,k)
    endif ; enddo ; enddo
    do i=is,ie ; Ttot(i) = Ttot(i) + T_chg(i) * htot(i) ; enddo
  endif ! absorbAllSW

end subroutine absorbRemainingSW


subroutine sumSWoverBands(G, h, eps, htot, opacity_band, nsw, j, dt, &
                          H_limit_fluxes, correctAbsorption, absorbAllSW, &
                          ksort, iPen_SW_bnd, netPen)

! This subroutine applies shortwave heating below the mixed layer.  In
! addition, it causes all of the remaining SW radiation to be absorbed,
! provided that the total water column thickness is greater than
! H_limit_fluxes.  For thinner water columns, the heating is scaled down
! proportionately, the assumption being that the remaining heating (which is
! left in Pen_SW) should go into an (absent for now) ocean bottom sediment layer.

  type(ocean_grid_type),                 intent(in)    :: G
  real, dimension(NIMEM_,NKMEM_),        intent(in)    :: h, eps
  real, dimension(NIMEM_),               intent(in)    :: htot
  real, dimension(:,:,:),                intent(in)    :: opacity_band
  integer,                               intent(in)    :: nsw
  integer,                               intent(in)    :: j
  real,                                  intent(in)    :: dt
  real,                                  intent(in)    :: H_limit_fluxes
  logical,                               intent(in)    :: correctAbsorption
  logical,                               intent(in)    :: absorbAllSW
  integer, dimension(NIMEM_,NKMEM_),     intent(in)    :: ksort
  real, dimension(:,:),                  intent(in)    :: iPen_SW_bnd
  real, dimension(NIMEM_,NK_INTERFACE_), intent(inout) :: netPen ! Units of K m

! Arguments:
!  (in)      G             = ocean grid structure
!  (in)      h             = layer thickness (m or kg/m2)
!                            units of h are referred to as H below.
!  (in)      eps           = (small) thickness that must remain in each layer, and
!                            which will not be subject to heating (H units)
!  (in)      htot          = total mixed layer thickness, in H.
!  (in)      opacity_band  = opacity in each band of penetrating shortwave
!                            radiation, in H-1. The indicies are band, i, k.
!  (in)      nsw           =  number of bands of penetrating shortwave radiation.
!  (in)      j             = j-index to work on
!  (in)      dt            = time step (seconds)
!  (inout)   ksort         = density-sorted k-indices.
!  (inout)   Pen_SW_bnd    = penetrating shortwave heating in each band that
!                            hits the bottom and will be redistributed through
!                            the water column (K H units) size nsw x NIMEM_.
!  (out)     netPen        = attenuated flux at interfaces, summed over bands (K m units)

  real :: h_heat(SZI_(G))     !  thickness of the water column that receives
                              !  remaining shortwave radiation, in H.
  real :: Pen_SW_rem(SZI_(G)) ! sum across all wavelength bands of the
                              ! penetrating shortwave heating that hits the bottom
                              ! and will be redistributed through the water column
                              ! (K H units)

  real, dimension(size(iPen_SW_bnd,1),size(iPen_SW_bnd,2)) :: Pen_SW_bnd
  real :: SW_trans        ! fraction of shortwave radiation that is not
                          ! absorbed in a layer (nondimensional)
  real :: unabsorbed      ! fraction of the shortwave radiation that
                          ! is not absorbed because the layers are too thin.
  real :: Ih_limit        ! inverse of the total depth at which the
                          ! surface fluxes start to be limited (1/H units)
  real :: h_min_heat      ! minimum thickness layer that should get heated (H units)
  real :: opt_depth       ! optical depth of a layer (non-dim)
  real :: exp_OD          ! exp(-opt_depth) (non-dim)
  real :: heat_bnd        ! heating due to absorption in the current
                          ! layer by the current band, including any piece that
                          ! is moved upward (K H units)
  real :: SWa             ! fraction of the absorbed shortwave that is
                          ! moved to layers above with correctAbsorption (non-dim)
  logical :: SW_Remains   ! If true, some column has shortwave radiation that
                          ! was not entirely absorbed.

  integer :: is, ie, nz, i, k, ks, n
  SW_Remains = .false.

  h_min_heat = 2.0*G%Angstrom + G%H_subroundoff
  is = G%isc ; ie = G%iec ; nz = G%ke

  pen_SW_bnd(:,:) = iPen_SW_bnd(:,:)
  do i=is,ie ; h_heat(i) = htot(i) ; enddo
  netPen(:,1) = sum( pen_SW_bnd(:,:), dim=1 ) ! Surface interface

  ! Apply penetrating SW radiation to remaining parts of layers.  Excessively thin
  ! layers are not heated.
  do ks=1,nz

    do i=is,ie ; if (ksort(i,ks) > 0) then
      k = ksort(i,ks)
      netPen(i,k+1) = 0.

      if (h(i,k) > 1.5*eps(i,k)) then
        do n=1,nsw ; if (Pen_SW_bnd(n,i) > 0.0) then
          ! SW_trans is the SW that is transmitted THROUGH the layer
          opt_depth = h(i,k) * opacity_band(n,i,k)
          exp_OD = exp(-opt_depth)
          SW_trans = exp_OD
          ! Heating at a rate of less than 10-4 W m-2 = 10-3 K m / Century,
          ! and of the layer in question less than 1 K / Century, can be
          ! absorbed without further penetration.
          if ((nsw*Pen_SW_bnd(n,i)*SW_trans < G%m_to_H*2.5e-11*dt) .and. &
              (nsw*Pen_SW_bnd(n,i)*SW_trans < h(i,k)*dt*2.5e-8)) &
            SW_trans = 0.0

          if (correctAbsorption .and. (h_heat(i) > 0.0)) then
            if (opt_depth > 1e-5) then
              SWa = ((opt_depth + (opt_depth + 2.0)*exp_OD) - 2.0) / &
                ((opt_depth + opacity_band(n,i,k) * h_heat(i)) * &
                 (1.0 - exp_OD))
            else
              ! Use a Taylor's series expansion of the expression above for a
              ! more accurate form with very small layer optical depths.
              SWa = h(i,k) * (opt_depth * (1.0 - opt_depth)) / &
                ((h_heat(i) + h(i,k)) * (6.0 - 3.0*opt_depth))
            endif
            Heat_bnd = Pen_SW_bnd(n,i) * (1.0 - SW_trans)
          endif

          Pen_SW_bnd(n,i) = Pen_SW_bnd(n,i) * SW_trans
          netPen(i,k+1) = netPen(i,k+1) + Pen_SW_bnd(n,i)
        endif ; enddo
      endif ! h(i,k) > 1.5*eps(i,k)

      ! Add to the accumulated thickness above that could be heated.
      ! Only layers greater than h_min_heat thick should get heated.
      if (h(i,k) >= 2.0*h_min_heat) then
        h_heat(i) = h_heat(i) + h(i,k)
      elseif (h(i,k) > h_min_heat) then
        h_heat(i) = h_heat(i) + (2.0*h(i,k) - 2.0*h_min_heat)
      endif
    endif ; enddo ! i loop
  enddo ! k loop

  ! Apply heating above the layers in which it should have occurred to get the
  ! correct mean depth of the shortwave that should be absorbed by each layer.
!    if (correctAbsorption) then
!    endif

! if (.not.absorbAllSW .and. .not.correctAbsorption) return
  if (absorbAllSW) then

    ! If there is still shortwave radiation at this point, it could go into
    ! the bottom (with a bottom mud model), or it could be redistributed back
    ! through the water column.
    do i=is,ie
      Pen_SW_rem(i) = Pen_SW_bnd(1,i)
      do n=2,nsw ; Pen_SW_rem(i) = Pen_SW_rem(i) + Pen_SW_bnd(n,i) ; enddo
    enddo
    do i=is,ie ; if (Pen_SW_rem(i) > 0.0) SW_Remains = .true. ; enddo

    Ih_limit = 1.0 / (H_limit_fluxes * G%m_to_H)
    do i=is,ie ; if ((Pen_SW_rem(i) > 0.0) .and. (h_heat(i) > 0.0)) then
      if (h_heat(i)*Ih_limit < 1.0) then
        unabsorbed = 1.0 - h_heat(i)*Ih_limit
      else
        unabsorbed = 0.0
      endif
      do n=1,nsw ; Pen_SW_bnd(n,i) = unabsorbed * Pen_SW_bnd(n,i) ; enddo
    endif ; enddo

  endif ! absorbAllSW

end subroutine sumSWoverBands


subroutine MOM_forcing_chksum(mesg, fluxes, G, haloshift)
! this subroutine writes out chksums for the model's basic state variables.

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
!  (in)      G - The ocean's grid structure.
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
  if (associated(fluxes%sens)) &
    call hchksum(fluxes%sens, mesg//" fluxes%sens",G,haloshift=hshift)
  if (associated(fluxes%evap)) &
    call hchksum(fluxes%evap, mesg//" fluxes%evap",G,haloshift=hshift)
  if (associated(fluxes%liq_precip)) &
    call hchksum(fluxes%liq_precip, mesg//" fluxes%liq_precip",G,haloshift=hshift)
  if (associated(fluxes%froz_precip)) &
    call hchksum(fluxes%froz_precip, mesg//" fluxes%froz_precip",G,haloshift=hshift)
  if (associated(fluxes%virt_precip)) &
    call hchksum(fluxes%virt_precip, mesg//" fluxes%virt_precip",G,haloshift=hshift)
  if (associated(fluxes%p_surf)) &
    call hchksum(fluxes%p_surf, mesg//" fluxes%p_surf",G,haloshift=hshift)
  if (associated(fluxes%salt_flux)) &
    call hchksum(fluxes%salt_flux, mesg//" fluxes%salt_flux",G,haloshift=hshift)
  if (associated(fluxes%TKE_tidal)) &
    call hchksum(fluxes%TKE_tidal, mesg//" fluxes%TKE_tidal",G,haloshift=hshift)
  if (associated(fluxes%ustar_tidal)) &
    call hchksum(fluxes%ustar_tidal, mesg//" fluxes%ustar_tidal",G,haloshift=hshift)
  if (associated(fluxes%liq_runoff)) &
    call hchksum(fluxes%liq_runoff, mesg//" fluxes%liq_runoff",G,haloshift=hshift)
  if (associated(fluxes%froz_runoff)) &
    call hchksum(fluxes%froz_runoff, mesg//" fluxes%froz_runoff",G,haloshift=hshift)
  if (associated(fluxes%runoff_hflx)) &
    call hchksum(fluxes%runoff_hflx, mesg//" fluxes%runoff_hflx",G,haloshift=hshift)
  if (associated(fluxes%calving_hflx)) &
    call hchksum(fluxes%calving_hflx, mesg//" fluxes%calving_hflx",G,haloshift=hshift)
end subroutine MOM_forcing_chksum

subroutine forcing_SinglePointPrint(fluxes, G, i, j, mesg)
! This subroutine writes out values of the fluxes arrays at the i,j location

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
  call locMsg(fluxes%sens,'sens')
  call locMsg(fluxes%evap,'evap')
  call locMsg(fluxes%liq_precip,'liq_precip')
  call locMsg(fluxes%froz_precip,'froz_precip')
  call locMsg(fluxes%virt_precip,'virt_precip')
  call locMsg(fluxes%p_surf,'p_surf')
  call locMsg(fluxes%salt_flux,'salt_flux')
  call locMsg(fluxes%TKE_tidal,'TKE_tidal')
  call locMsg(fluxes%ustar_tidal,'ustar_tidal')
  call locMsg(fluxes%liq_runoff,'liq_runoff')
  call locMsg(fluxes%froz_runoff,'froz_runoff')
  call locMsg(fluxes%runoff_hflx,'runoff_hflx')
  call locMsg(fluxes%calving_hflx,'calving_hflx')

  contains

  subroutine locMsg(array,aname)
  real, dimension(:,:), pointer :: array
  character(len=*) :: aname
  if (associated(array)) then
    write(0,'(3a,es15.3)') 'MOM_forcing_type, forcing_SinglePointPrint: ',trim(aname),' = ',array(i,j)
  else
    write(0,'(4a)') 'MOM_forcing_type, forcing_SinglePointPrint: ',trim(aname),' is not associated.'
  endif
  end subroutine locMsg
end subroutine forcing_SinglePointPrint

!> Deallocates the forcing type
subroutine deallocate_forcing_type(fluxes)
  type(forcing), intent(inout) :: fluxes
  if (associated(fluxes%taux))           deallocate(fluxes%taux)
  if (associated(fluxes%tauy))           deallocate(fluxes%tauy)
  if (associated(fluxes%ustar))          deallocate(fluxes%ustar)
  if (associated(fluxes%buoy))           deallocate(fluxes%buoy)
  if (associated(fluxes%sw))             deallocate(fluxes%sw)
  if (associated(fluxes%sw_vis_dir))     deallocate(fluxes%sw_vis_dir)
  if (associated(fluxes%sw_vis_dif))     deallocate(fluxes%sw_vis_dif)
  if (associated(fluxes%sw_nir_dir))     deallocate(fluxes%sw_nir_dir)
  if (associated(fluxes%sw_nir_dif))     deallocate(fluxes%sw_nir_dif)
  if (associated(fluxes%lw))             deallocate(fluxes%lw)
  if (associated(fluxes%latent))         deallocate(fluxes%latent)
  if (associated(fluxes%sens))           deallocate(fluxes%sens)
  if (associated(fluxes%heat_restore))   deallocate(fluxes%heat_restore)
  if (associated(fluxes%runoff_hflx))    deallocate(fluxes%runoff_hflx)
  if (associated(fluxes%calving_hflx))   deallocate(fluxes%calving_hflx)
  if (associated(fluxes%evap))           deallocate(fluxes%evap)
  if (associated(fluxes%liq_precip))     deallocate(fluxes%liq_precip)
  if (associated(fluxes%froz_precip))    deallocate(fluxes%froz_precip)
  if (associated(fluxes%virt_precip))    deallocate(fluxes%virt_precip)
  if (associated(fluxes%liq_runoff))     deallocate(fluxes%liq_runoff)
  if (associated(fluxes%froz_runoff))    deallocate(fluxes%froz_runoff)
  if (associated(fluxes%salt_flux))      deallocate(fluxes%salt_flux)
  if (associated(fluxes%p_surf_full))    deallocate(fluxes%p_surf_full)
  if (associated(fluxes%p_surf))         deallocate(fluxes%p_surf)
  if (associated(fluxes%TKE_tidal))      deallocate(fluxes%TKE_tidal)
  if (associated(fluxes%ustar_tidal))    deallocate(fluxes%ustar_tidal)
  if (associated(fluxes%ustar_shelf))    deallocate(fluxes%ustar_shelf)
  if (associated(fluxes%frac_shelf_h))   deallocate(fluxes%frac_shelf_h)
  if (associated(fluxes%frac_shelf_u))   deallocate(fluxes%frac_shelf_u)
  if (associated(fluxes%frac_shelf_v))   deallocate(fluxes%frac_shelf_v)
  if (associated(fluxes%rigidity_ice_u)) deallocate(fluxes%rigidity_ice_u)
  if (associated(fluxes%rigidity_ice_v)) deallocate(fluxes%rigidity_ice_v)
  if (associated(fluxes%tr_fluxes))      deallocate(fluxes%tr_fluxes)
end subroutine deallocate_forcing_type

end module MOM_forcing_type
