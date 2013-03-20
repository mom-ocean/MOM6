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

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : time_type, diag_ptrs
use MOM_domains, only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : forcing, thermo_var_ptrs, optics_type

implicit none ; private

#include <MOM_memory.h>

public extractFluxes1d

integer :: num_msg = 0, max_msg = 2

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

end module MOM_forcing_type
