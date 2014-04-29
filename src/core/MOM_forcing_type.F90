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

use MOM_checksums, only : hchksum, qchksum, uchksum, vchksum
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : time_type
use MOM_domains, only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_EOS, only : calculate_density_derivs
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs

use coupler_types_mod, only : coupler_2d_bc_type

implicit none ; private

#include <MOM_memory.h>

public extractFluxes1d, extractFluxes2d
public MOM_forcing_chksum, absorbRemainingSW
public calculateBuoyancyFlux1d, calculateBuoyancyFlux2d
public forcing_SinglePointPrint

integer :: num_msg = 0, max_msg = 2

!   The following structure contains pointers to the forcing fields
! which may be used to drive MOM.  All fluxes are positive downward.
! Pointers to unused fluxes should be set to NULL.
type, public :: forcing
  real, pointer, dimension(:,:) :: &
    taux => NULL(), &       ! The zonal wind stress, in Pa.
    tauy => NULL(), &       ! The meridional wind stress, in Pa.
    ustar => NULL(), &      ! The surface friction velocity, in units of m s-1.
    buoy => NULL(), &       ! The buoyancy flux into the ocean in m2 s-3.

    sw => NULL(), &         ! The shortwave heat flux into the ocean, in W m-2.
    sw_vis_dir => NULL(), & ! The visible, direct shortwave heat flux into the
                            ! ocean, in W m-2.
    sw_vis_dif => NULL(), & ! The visible, diffuse shortwave heat flux into the
                            ! ocean, in W m-2.
    sw_nir_dir => NULL(), & ! The near-IR, direct shortwave heat flux into the
                            ! ocean, in W m-2.
    sw_nir_dif => NULL(), & ! The near-IR, diffuse shortwave heat
                            ! flux into the ocean, in W m-2.
    lw => NULL(), &         ! The longwave heat flux into the ocean, in W m-2.
                            ! This field is typically negative.
    latent => NULL(), &     ! The latent heat flux into the ocean, in W m-2.
                            ! This field is typically negative.
    sens => NULL(), &       ! The sensible heat flux into the ocean, in W m-2.
                            ! This field is typically negative.
    heat_restore => NULL(), & ! The heat flux into the ocean from temperature
                            ! restoring, in W m-2.
                            ! This field is typically negative.

    evap => NULL(), &       ! The negative of the fresh water flux out of the
                            ! ocean, in kg m-2 s-1.
    liq_precip => NULL(), & ! The liquid water flux into the ocean, in
                            ! kg m-2 s-1.
    froz_precip => NULL(), &  ! The frozen water flux into the ocean,
                            ! in kg m-2 s-1.
    virt_precip => NULL(), & ! The virtual water flux into the ocean associated
                            ! with salinity restoring, in kg m-2 s-1.
    runoff_hflx => NULL(), & ! Heat flux associated with liq_runoff in W m-2.
    calving_hflx => NULL(), & ! Heat flux associated with froz_runoff in W m-2.

    p_surf_full => NULL(), & ! The pressure at the top ocean interface, in Pa.
                             ! If there is sea-ice, this is at the interface
                             ! between the ice and ocean.
    p_surf => NULL(), &      ! The pressure at the top ocean interface, in Pa,
                             ! as used to drive the ocean model.  If p_surf is
                             ! limited, this may be smaller than p_surf_full,
                             ! otherwise they are the same.
    salt_flux => NULL(), &   ! The net salt flux into the ocean in kg Salt m-2 s-1.
    TKE_tidal => NULL(), &   ! The tidal source of energy driving mixing in the
                             ! bottom boundary layer, in W m-2.
    ustar_tidal => NULL(), & ! The tidal contribution to bottom ustar, in m s-1.
    liq_runoff => NULL(), &  ! Mass of river runoff in units of kg m-2 s-1.
    froz_runoff => NULL(), & ! Mass of calving in units of kg m-2 s-1.

    ustar_shelf => NULL(), & ! The friction velocity under ice-shelves in m s-1.
                             ! This was calculated by the ocean the previous
                             ! time step.
    frac_shelf_h => NULL(), &! Fractional ice shelf coverage of h-, u-, and v-
    frac_shelf_u => NULL(), &! cells, nondimensional from 0 to 1. These are only
    frac_shelf_v => NULL(), &! associated if ice shelves are enabled, and are
                             ! exactly 0 away from shelves or on land.
    rigidity_ice_u => NULL(),& ! The depth-integrated lateral viscosity of
    rigidity_ice_v => NULL()   ! ice shelves at u- or v-points, in m3 s-1.
  real :: C_p            !   The heat capacity of seawater, in J K-1 kg-1.
                         ! This is always set to the same value as is found in
                         ! the thermovar_ptrs_type.
  type(coupler_2d_bc_type), pointer :: tr_fluxes  => NULL()
                                            ! A structure that may contain an
                                            ! array of named fields used for
                                            ! passive tracer fluxes.
       !!! NOTE: ALL OF THE ARRAYS IN TR_FLUXES USE THE COUPLER'S INDEXING
       !!!       CONVENTION AND HAVE NO HALOS!  THIS IS DONE TO CONFORM TO
       !!!       THE TREATMENT IN MOM4, BUT I DON'T LIKE IT!
end type forcing


type, public :: optics_type
  integer :: nbands    ! The number of penetrating bands of SW radiation.
  real, pointer, dimension(:,:,:,:) :: &
    opacity_band => NULL()  ! The SW optical depth per unit thickness, m-1.
                       ! The number of bands of radiation is the most rapidly
                       ! varying (first) index.
  real, pointer, dimension(:,:,:) :: &
    SW_pen_band => NULL()  ! The shortwave radiation at the surface in each of
                       ! the nbands bands that penetrates beyond the surface,
                       ! in W m-2. The most rapidly varying dimension is the
                       ! band.
  real, pointer, dimension(:) :: &
    min_wavelength_band => NULL(), & ! The range of wavelengths in each band of
    max_wavelength_band => NULL()    ! penetrating shortwave radiation, in nm.
end type optics_type



contains

subroutine extractFluxes1d(G, fluxes, optics, nsw, j, dt, &
                  DepthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
                  h, T, Net_H, Net_heat, Net_salt, Pen_SW_bnd, tv)
  type(ocean_grid_type),         intent(in)  :: G
  type(forcing),                 intent(in)  :: fluxes
  type(optics_type),             pointer     :: optics
  integer,                       intent(in)  :: nsw, j
  real,                          intent(in)  :: dt, DepthBeforeScalingFluxes
  logical,                       intent(in)  :: useRiverHeatContent, &
                                                useCalvingHeatContent
  real, dimension(NIMEM_,NKMEM_), intent(in) :: h, T
  real, dimension(NIMEM_),       intent(out) :: Net_H, Net_heat, Net_salt
  real, dimension(:,:),          intent(out) :: Pen_SW_bnd
  type(thermo_var_ptrs),       intent(inout) :: tv
!   This subroutine extracts the relevant buoyancy fluxes from the surface
! fluxes type.

!  (in)      G - The ocean's grid structure.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      nsw - The number of bands of penetrating shortwave radiation.
!  (in)      j - The j-index to work on.
!  (in)      dt - The time step in s.
!  (in)      DepthBeforeScalingFluxes - Minimum ocean thickness to allow before scaling away fluxes
!  (in)      useRiverHeatContent - If true, apply river heat content
!  (in)      useCalvingHeatContent& - If true, apply calving heat content
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      T - Layer temperatures, in deg C.
!  (out)     Net_H - The net mass flux (if non-Boussinsq) or volume flux (if
!                    Boussinesq - i.e. the fresh water flux (P+R-E)) into the
!                    ocean over a time step, in H.
!  (out)     Net_heat - The net heating at the surface over a time step,
!                       exclusive of heating that appears in Pen_SW, in K H.
!  (out)     Net_salt - The surface salt flux into the ocean over a time step, psu H.
!  (out)     Pen_SW_bnd - The penetrating shortwave heating at the sea surface
!                         in each penetrating band, in K H, size nsw x NIMEM_.
!  (inout)   tv - A structure containing pointers to any available
!                 thermodynamic fields. Here it is used to keep track of the
!                 actual heat flux associated with net mass fluxes into the
!                 ocean.

  real :: htot(SZI_(G))  ! The total depth of the ocean in m or kg m-2.
  real :: Pen_sw_tot(SZI_(G))  ! The sum across all bands of Pen_SW, K H.
  real :: Ih_limit  !   The inverse of the total depth at which the
                    ! surface fluxes start to be limited, in H-1.
  real :: scale     !   Scale scales away the fluxes when the total depth is
                    ! thinner than DepthBeforeScalingFluxes.
  real :: Irho_cp       !  1.0 / rho_0 * c_p.
  real :: Irho0, I_Cp
  character(len=200) :: mesg
  integer :: is, ie, nz, i, k, n
  Ih_limit = 1.0 / (DepthBeforeScalingFluxes * G%m_to_H)
  Irho0 = 1.0 / G%Rho0 ; I_Cp = 1.0 / fluxes%C_p
  Irho_cp = 1.0 / (G%H_to_kg_m2 * fluxes%C_p)
  is = G%isc ; ie = G%iec ; nz = G%ke

  ! A whole load of error checking
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

    ! Convert the penetrating shortwave forcing to K m.
    Pen_sw_tot(i) = 0.0
    if (nsw >= 1) then
      do n=1,nsw
        Pen_SW_bnd(n,i) = Irho_cp*scale*dt * &
            max(0.0, optics%sw_pen_band(n,i,j))
        Pen_sw_tot(i) = Pen_sw_tot(i) + Pen_SW_bnd(n,i)
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
    if (.not.G%Boussinesq .and. ASSOCIATED(fluxes%salt_flux)) &
      Net_H(i) = Net_H(i) + (dt * G%kg_m2_to_H) * &
                    (scale * fluxes%salt_flux(i,j))

    ! Heat fluxes
    Net_heat(i) = scale * dt * Irho_cp * ( fluxes%sw(i,j) + &
         ((fluxes%lw(i,j) + fluxes%latent(i,j)) + fluxes%sens(i,j)) )

    if (ASSOCIATED(fluxes%heat_restore)) Net_heat(i) = Net_heat(i) + &
         (scale * (dt * Irho_cp)) * fluxes%heat_restore(i,j)

    if (useRiverHeatContent) then
      Net_heat(i) = (Net_heat(i) + (scale*(dt*Irho_cp)) * fluxes%runoff_hflx(i,j)) - &
                     (G%kg_m2_to_H * (scale * dt)) * fluxes%liq_runoff(i,j) * T(i,1)
      if (ASSOCIATED(tv%TempxPmE)) then
        tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + (scale * dt) * &
            (I_Cp*fluxes%runoff_hflx(i,j) - fluxes%liq_runoff(i,j)*T(i,1))
      endif    
    endif

    if (useCalvingHeatContent) then
      Net_heat(i) = Net_heat(i) + (scale*(dt*Irho_cp)) * fluxes%calving_hflx(i,j) - &
                    (G%kg_m2_to_H * (scale * dt)) * fluxes%froz_runoff(i,j) * T(i,1)
      if (ASSOCIATED(tv%TempxPmE)) then
        tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + (scale * dt) * &
            (I_Cp*fluxes%calving_hflx(i,j) - fluxes%froz_runoff(i,j)*T(i,1))
      endif    
    endif

    if (num_msg<max_msg) then
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


subroutine extractFluxes2d(G, fluxes, optics, nsw, dt, &
                  DepthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
                  h, T, Net_H, Net_heat, Net_salt, Pen_SW_bnd, tv)
  type(ocean_grid_type),         intent(in)  :: G
  type(forcing),                 intent(in)  :: fluxes
  type(optics_type),             pointer     :: optics
  integer,                       intent(in)  :: nsw
  real,                          intent(in)  :: dt, DepthBeforeScalingFluxes
  logical,                       intent(in)  :: useRiverHeatContent, &
                                                useCalvingHeatContent
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: h, T
  real, dimension(NIMEM_,NJMEM_),       intent(out) :: Net_H, Net_heat, Net_salt
  real, dimension(:,:,:),        intent(out) :: Pen_SW_bnd
  type(thermo_var_ptrs),       intent(inout) :: tv
!   This subroutine extracts the relevant buoyancy fluxes from the surface
! fluxes type.

!  (in)      G - The ocean's grid structure.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      nsw - The number of bands of penetrating shortwave radiation.
!  (in)      j - The j-index to work on.
!  (in)      dt - The time step in s.
!  (in)      DepthBeforeScalingFluxes - Minimum ocean thickness to allow before scaling away fluxes
!  (in)      useRiverHeatContent - If true, apply river heat content
!  (in)      useCalvingHeatContent& - If true, apply calving heat content
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      T - Layer temperatures, in deg C.
!  (out)     Net_H - The net mass flux (if non-Boussinsq) or volume flux (if
!                    Boussinesq - i.e. the fresh water flux (P+R-E)) into the
!                    ocean over a time step, in H.
!  (out)     Net_heat - The net heating at the surface over a time step,
!                       exclusive of heating that appears in Pen_SW, in K H.
!  (out)     Net_salt - The surface salt flux into the ocean over a time step, psu H.
!  (out)     Pen_SW_bnd - The penetrating shortwave heating at the sea surface
!                         in each penetrating band, in K H, size nsw x NIMEM_.
!  (inout)   tv - A structure containing pointers to any available
!                 thermodynamic fields. Here it is used to keep track of the
!                 actual heat flux associated with net mass fluxes into the
!                 ocean.
  integer :: j

  do j=G%jsc, G%jec
    call extractFluxes1d(G, fluxes, optics, nsw, j, dt, &
            DepthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
            h(:,j,:), T(:,j,:), Net_H(:,j), Net_heat(:,j), Net_salt(:,j), Pen_SW_bnd(:,:,j), tv)
  enddo

end subroutine extractFluxes2d


subroutine calculateBuoyancyFlux1d(G, fluxes, optics, h, Temp, Salt, tv, j, buoyancyFlux, netHeatMinusSW, netSalt )
! Calculates buoyancy flux by adding up the heat, FW and salt fluxes and linearizing
! about the surface state.

! Arguments
  type(ocean_grid_type),                 intent(in)    :: G      ! Ocean grid
  type(forcing),                         intent(in)    :: fluxes ! Surface fluxes/forcing type
  type(optics_type),                     pointer       :: optics ! Optics for penetrating SW
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h      ! Layer/level thicknesses (units of H)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: Temp   ! Pot. temperature (degrees C)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: Salt   ! Salinity (ppt)
  type(thermo_var_ptrs),                 intent(inout) :: tv     ! Thermodynamics type (out needed for tv%TempxPmE ????)
  integer,                               intent(in)    :: j      ! j-index of row to work on
  real, dimension(NIMEM_,NK_INTERFACE_), intent(inout) :: buoyancyFlux ! Buoyancy flux (m2/s3)
  real, dimension(NIMEM_),               intent(inout) :: netHeatMinusSW ! Heat flux excluding SW (K m/s)
  real, dimension(NIMEM_),               intent(inout) :: netSalt ! Salt flux (ppt m/s)

! Local variables
  integer :: nsw, start, npts, k
  real, parameter :: dt = 1. ! This is set to unity to return a rate from extractFluxes1d
  real, dimension( SZI_(G) ) :: netH, netHeat ! FW, heat fluxes in (m/s, K m/s, ppt m/s)
  real, dimension( optics%nbands, SZI_(G) ) :: penSWbnd ! SW penetration bands
  real, dimension( SZI_(G) ) :: pressure ! Pressurea the surface ( Pa )
  real, dimension( SZI_(G) ) :: dRhodT, dRhodS ! Derivatives of density
  logical :: useRiverHeatContent, useCalvingHeatContent
  real :: depthBeforeScalingFluxes, GoRho
  integer, dimension(SZI_(G),SZK_(G)) :: ksort
  real, dimension(SZI_(G),SZK_(G)+1) :: netPen
  real    :: H_limit_fluxes

  nsw = optics%nbands
  useRiverHeatContent = .False.    !  ????????????????
  useCalvingHeatContent = .False.  !  ????????????????
  depthBeforeScalingFluxes = max( G%Angstrom, 1.e-30 )
  pressure(:) = 0. ! Ignore atmospheric pressure
  GoRho = G%g_Earth / G%Rho0
  start = 1 + G%isc - G%isd
  npts = 1 + G%iec - G%isc

  do k=1, G%ke
    ksort(:,k) = k
  enddo
  H_limit_fluxes = depthBeforeScalingFluxes ! ?????????????

  ! Fetch the fresh-water, heat and salt fluxes
  ! netH is the fresh-water flux
  ! netSalt is the salt flux (typically zero except under sea-ice)
  ! netHeat is the heat flux EXCEPT the penetrating SW
  ! penSWbnd is the surface SW for each band
  call extractFluxes1d(G, fluxes, optics, nsw, j, dt, &
                depthBeforeScalingFluxes, useRiverHeatContent, useCalvingHeatContent, &
                h(:,j,:), Temp(:,j,:), netH, netHeatMinusSW, netSalt, penSWbnd, tv)

  ! Sum over bands and attenuate as a function of depth
  ! netPen is the netSW as a function of depth
  call sumSWoverBands(G, h(:,j,:), 0.*h(:,j,:), 0.*h(:,j,1), optics%opacity_band(:,:,j,:), nsw, j, dt, &
                      H_limit_fluxes, .true., .true., &
                      ksort, penSWbnd, netPen)

  ! Density derivatives
  call calculate_density_derivs(Temp(:,j,1), Salt(:,j,1), pressure, &
                                dRhodT, dRhodS, start, npts, tv%eqn_of_state)

  ! Adjust netSalt to reflect dillution effect of FW flux
  netSalt(:) = netSalt(:) - Salt(:,j,1) * netH * G%H_to_m

  ! Add in the SW heating for purposes of calculating the net surface buoyancy flux
 !netHeat(:) = netHeatMinusSW(:) + sum( penSWbnd(:,:), dim=1 )
  netHeat(:) = netHeatMinusSW(:) + netPen(:,1)

  ! Convert to a buoyancy flux, excluding penetrating SW heating
  buoyancyFlux(:,1) = - GoRho * ( dRhodS(:) * netSalt(:) + dRhodT(:) * netHeat(:) ) ! m2/s3
  ! We also have a penetrative buoyancy flux associaed with penetrative SW
  do k=2, G%ke+1
    buoyancyFlux(:,k) = - GoRho * ( dRhodT(:) * netPen(:,k) ) ! m2/s3
  enddo

end subroutine calculateBuoyancyFlux1d


subroutine calculateBuoyancyFlux2d(G, fluxes, optics, h, Temp, Salt, tv, buoyancyFlux, netHeatMinusSW, netSalt )
! Calculates buoyancy flux by adding up the heat, FW and salt fluxes and linearizing
! about the surface state.

! Arguments
  type(ocean_grid_type),                       intent(in)    :: G      ! Ocean grid
  type(forcing),                               intent(in)    :: fluxes ! Surface fluxes/forcing type
  type(optics_type),                           pointer       :: optics ! Optics for penetrating SW
  real, dimension(NIMEM_,NJMEM_,NKMEM_),       intent(in)    :: h      ! Layer/level thicknesses (units of H)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),       intent(in)    :: Temp   ! Pot. temperature (degrees C)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),       intent(in)    :: Salt   ! Salinity (ppt)
  type(thermo_var_ptrs),                       intent(inout) :: tv     ! Thermodynamics type (out needed for tv%TempxPmE ????)
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_),intent(inout) :: buoyancyFlux ! Buoyancy flux (m2/s3)
  real, dimension(NIMEM_,NJMEM_),optional,     intent(inout) :: netHeatMinusSW ! Heat flux excluding SW (K m/s)
  real, dimension(NIMEM_,NJMEM_),optional,     intent(inout) :: netSalt ! Salt flux (ppt m/s)

! Local variables
  real, dimension( SZI_(G) ) :: netT, netS ! Fluxes in (K m/s, ppt m/s)
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
  type(ocean_grid_type),          intent(in)    :: G
  real, dimension(NIMEM_,NKMEM_), intent(in)    :: h, eps
  real, dimension(NIMEM_),        intent(in)    :: htot
  real, dimension(:,:,:),         intent(in)    :: opacity_band
  integer,                        intent(in)    :: nsw
  integer,                        intent(in)    :: j
  real,                           intent(in)    :: dt
  real,                           intent(in)    :: H_limit_fluxes
  logical,                        intent(in)    :: correctAbsorption
  logical,                        intent(in)    :: absorbAllSW
  integer, dimension(NIMEM_,NKMEM_), intent(in) :: ksort
  real, dimension(NIMEM_,NKMEM_), intent(inout) :: T
  real, dimension(NIMEM_),        intent(inout) :: Ttot
  real, dimension(:,:),           intent(inout) :: Pen_SW_bnd
!   This subroutine applies shortwave heating below the mixed layer.  In
! addition, it causes all of the remaining SW radiation to be absorbed,
! provided that the total water column thickness is greater than
! H_limit_fluxes.  For thinner water columns, the heating is scaled down
! proportionately, the assumption being that the remaining heating (which is
! left in Pen_SW) should go into an (unincluded) ocean bottom sediment layer.

! Arguments:
!  (in)      G - The ocean's grid structure.
!  (in)      h - Layer thickness, in m or kg m-2. (Intent in)  The units
!                of h are referred to as H below.
!  (in)      eps - The (small) thickness that must remain in each layer, and
!                  which will not be subject to heating, in H.
!  (in)      htot - The total mixed layer thickness, in H.
!  (in)      opacity_band - The opacity in each band of penetrating shortwave
!                           radiation, in H-1. The indicies are band, i, k.
!  (in)      nsw - The number of bands of penetrating shortwave radiation.
!  (in)      j - The j-index to work on.
!  (in)      dt - The time step in s.
!  (inout)   ksort - The density-sorted k-indicies.
!  (inout)   T - The layer potential temperatures, in deg C.
!  (inout)   Ttot - The depth integrated mixed layer temperature, in K H.
!  (inout)   Pen_SW_bnd - The penetrating shortwave heating in each band that
!                         hits the bottom and will be redistributed through
!                         the water column, in K H, size nsw x NIMEM_.

  real :: h_heat(SZI_(G)) !   The thickness of the water column that receives
                          ! the remaining shortwave radiation, in H.
  real :: T_chg_above(SZI_(G),SZK_(G)) !  The temperature change of all the
                          ! thick layers above a given layer, in K.  The net
                          ! change in the heat of a layer is the sum of
                          ! T_chg_above from all the layers below, plus any
                          ! contribution from absorbing radiation that hits
                          ! the bottom.
  real :: T_chg(SZI_(G))  !   The temperature change of thick layers due to
                          ! the remaining shortwave radiation, in K.
  real :: Pen_SW_rem(SZI_(G)) ! The sum across all wavelength bands of the
                          ! penetrating shortwave heating that hits the bottom
                          ! and will be redistributed through the water column
                          ! in units of K H.
  real :: SW_trans  !   The fraction of shortwave radiation that is not
                    ! absorbed in a layer, nondimensional.
  real :: unabsorbed      !   The fraction of the shortwave radiation that
                          ! is not absorbed because the layers are too thin.
  real :: Ih_limit        !   The inverse of the total depth at which the
                          ! surface fluxes start to be limited, in H-1.
  real :: h_min_heat      !   The minimum thickness layer that should get
                          ! heated, in H.
  real :: opt_depth       !   The optical depth of a layer, non-dim.
  real :: exp_OD          !   exp(-opt_depth), non-dim.
  real :: heat_bnd        !   The heating due to absorption in the current
                          ! layer by the current band, including any piece that
                          ! is moved upward, in K H.
  real :: SWa             !   The fraction of the absorbed shortwave that is
                          ! moved to layers above with correctAbsorption, ND.
  logical :: SW_Remains   ! If true, some column has shortwave radiation that
                          ! was not entirely absorbed.
  integer :: is, ie, nz, i, k, ks, n
  SW_Remains = .false.

  h_min_heat = 2.0*G%Angstrom + G%H_subroundoff
  is = G%isc ; ie = G%iec ; nz = G%ke

  do i=is,ie ; h_heat(i) = htot(i) ; enddo

  ! Apply penetrating SW radiation to remaining parts of layers.  Excessively thin
  ! isopycnal layers are not heated.
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
  type(ocean_grid_type),          intent(in)    :: G
  real, dimension(NIMEM_,NKMEM_), intent(in)    :: h, eps
  real, dimension(NIMEM_),        intent(in)    :: htot
  real, dimension(:,:,:),         intent(in)    :: opacity_band
  integer,                        intent(in)    :: nsw
  integer,                        intent(in)    :: j
  real,                           intent(in)    :: dt
  real,                           intent(in)    :: H_limit_fluxes
  logical,                        intent(in)    :: correctAbsorption
  logical,                        intent(in)    :: absorbAllSW
  integer, dimension(NIMEM_,NKMEM_), intent(in) :: ksort
  real, dimension(:,:),           intent(in)    :: iPen_SW_bnd
  real, dimension(NIMEM_,NK_INTERFACE_), intent(inout) :: netPen ! Units of K m
!   This subroutine applies shortwave heating below the mixed layer.  In
! addition, it causes all of the remaining SW radiation to be absorbed,
! provided that the total water column thickness is greater than
! H_limit_fluxes.  For thinner water columns, the heating is scaled down
! proportionately, the assumption being that the remaining heating (which is
! left in Pen_SW) should go into an (unincluded) ocean bottom sediment layer.

! Arguments:
!  (in)      G - The ocean's grid structure.
!  (in)      h - Layer thickness, in m or kg m-2. (Intent in)  The units
!                of h are referred to as H below.
!  (in)      eps - The (small) thickness that must remain in each layer, and
!                  which will not be subject to heating, in H.
!  (in)      htot - The total mixed layer thickness, in H.
!  (in)      opacity_band - The opacity in each band of penetrating shortwave
!                           radiation, in H-1. The indicies are band, i, k.
!  (in)      nsw - The number of bands of penetrating shortwave radiation.
!  (in)      j - The j-index to work on.
!  (in)      dt - The time step in s.
!  (inout)   ksort - The density-sorted k-indicies.
!  (inout)   Pen_SW_bnd - The penetrating shortwave heating in each band that
!                         hits the bottom and will be redistributed through
!                         the water column, in K H, size nsw x NIMEM_.
!  (out)     netPen is the attenuated flux at interfaces, summed over bands (units of K m)

  real :: h_heat(SZI_(G)) !   The thickness of the water column that receives
                          ! the remaining shortwave radiation, in H.
  real :: Pen_SW_rem(SZI_(G)) ! The sum across all wavelength bands of the
                          ! penetrating shortwave heating that hits the bottom
                          ! and will be redistributed through the water column
                          ! in units of K H.
  real, dimension(size(iPen_SW_bnd,1),size(iPen_SW_bnd,2)) :: Pen_SW_bnd
  real :: SW_trans  !   The fraction of shortwave radiation that is not
                    ! absorbed in a layer, nondimensional.
  real :: unabsorbed      !   The fraction of the shortwave radiation that
                          ! is not absorbed because the layers are too thin.
  real :: Ih_limit        !   The inverse of the total depth at which the
                          ! surface fluxes start to be limited, in H-1.
  real :: h_min_heat      !   The minimum thickness layer that should get
                          ! heated, in H.
  real :: opt_depth       !   The optical depth of a layer, non-dim.
  real :: exp_OD          !   exp(-opt_depth), non-dim.
  real :: heat_bnd        !   The heating due to absorption in the current
                          ! layer by the current band, including any piece that
                          ! is moved upward, in K H.
  real :: SWa             !   The fraction of the absorbed shortwave that is
                          ! moved to layers above with correctAbsorption, ND.
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
  ! isopycnal layers are not heated.
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
  character(len=*),                    intent(in) :: mesg
  type(forcing),                       intent(in) :: fluxes
  type(ocean_grid_type),               intent(in) :: G
  integer, optional,                   intent(in) :: haloshift
!   This subroutine writes out chksums for the model's basic state variables.
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
  type(forcing),                       intent(in) :: fluxes
  type(ocean_grid_type),               intent(in) :: G
  character(len=*),                    intent(in) :: mesg
  integer,                             intent(in) :: i, j
!   This subroutine writes out values of the fluxes arrays at
! the i,j location

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

end module MOM_forcing_type
