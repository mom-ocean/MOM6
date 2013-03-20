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
use MOM_diag_mediator, only : time_type, diag_ptrs
use MOM_domains, only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs

use coupler_types_mod, only : coupler_2d_bc_type

implicit none ; private

#include <MOM_memory.h>

public extractFluxes1d, MOM_forcing_chksum

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

!  Convert the penetrating shortwave forcing to K m.
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

    Net_H(i) = dt * (scale * ((((( fluxes%liq_precip(i,j)    &
                                 + fluxes%froz_precip(i,j) ) &
                                 + fluxes%evap(i,j)        ) &
                                 + fluxes%liq_runoff(i,j)  ) &
                                 + fluxes%virt_precip(i,j) ) &
                                 + fluxes%froz_runoff(i,j) ) )
    Net_H(i) = G%kg_m2_to_H * Net_H(i)
    if (.not.G%Boussinesq .and. associated(fluxes%salt_flux)) &
      Net_H(i) = Net_H(i) + (dt * G%kg_m2_to_H) * &
                    (scale * fluxes%salt_flux(i,j))

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

    Net_salt(i) = 0.0
    ! Convert salt_flux from kg (salt) m-2 s-1 to PSU m s-1.
    if (associated(fluxes%salt_flux)) &
      Net_salt(i) = (scale * dt * (1000.0 * fluxes%salt_flux(i,j))) * G%kg_m2_to_H

  enddo

end subroutine extractFluxes1d

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
end subroutine MOM_forcing_chksum

end module MOM_forcing_type
