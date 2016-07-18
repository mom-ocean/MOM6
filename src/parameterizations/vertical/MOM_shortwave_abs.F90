module MOM_shortwave_abs
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

use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public absorbRemainingSW, sumSWoverBands

type, public :: optics_type
  ! ocean optical properties

  integer :: nbands    ! number of penetrating bands of SW radiation

  real, pointer, dimension(:,:,:,:) :: &
    opacity_band => NULL()  ! SW optical depth per unit thickness (1/m)
                            ! Number of radiation bands is most rapidly varying (first) index.

  real, pointer, dimension(:,:,:) :: &
    SW_pen_band  => NULL()  ! shortwave radiation (W/m^2) at the surface in each of
                            ! the nbands bands that penetrates beyond the surface.
                            ! The most rapidly varying dimension is the band.

  real, pointer, dimension(:) ::   &
    min_wavelength_band => NULL(), & ! The range of wavelengths in each band of
    max_wavelength_band => NULL()    ! penetrating shortwave radiation (nm)

end type optics_type


contains

!> Apply shortwave heating below surface boundary layer.
subroutine absorbRemainingSW(G, GV, h, opacity_band, nsw, j, dt, H_limit_fluxes, &
                             adjustAbsorptionProfile, absorbAllSW, T, Pen_SW_bnd, &
                             eps, ksort, htot, Ttot, TKE, dSV_dT)

! This subroutine applies shortwave heating below the boundary layer (when running
! with the bulk mixed layer from GOLD) or throughout the water column.  In
! addition, it causes all of the remaining SW radiation to be absorbed,
! provided that the total water column thickness is greater than
! H_limit_fluxes.  For thinner water columns, the heating is scaled down
! proportionately, the assumption being that the remaining heating (which is
! left in Pen_SW) should go into an (absent for now) ocean bottom sediment layer.

  type(ocean_grid_type),             intent(in)    :: G
  type(verticalGrid_type),           intent(in)    :: GV
  real, dimension(SZI_(G),SZK_(G)),  intent(in)    :: h
  real, dimension(:,:,:),            intent(in)    :: opacity_band
  integer,                           intent(in)    :: nsw
  integer,                           intent(in)    :: j
  real,                              intent(in)    :: dt
  real,                              intent(in)    :: H_limit_fluxes
  logical,                           intent(in)    :: adjustAbsorptionProfile
  logical,                           intent(in)    :: absorbAllSW
  real, dimension(SZI_(G),SZK_(G)),  intent(inout) :: T
  real, dimension(:,:),              intent(inout) :: Pen_SW_bnd
  real, dimension(SZI_(G),SZK_(G)),    optional, intent(in)    :: eps
  integer, dimension(SZI_(G),SZK_(G)), optional, intent(in)    :: ksort
  real, dimension(SZI_(G)),            optional, intent(in)    :: htot
  real, dimension(SZI_(G)),            optional, intent(inout) :: Ttot
  real, dimension(SZI_(G),SZK_(G)),    optional, intent(in)    :: dSV_dT
  real, dimension(SZI_(G),SZK_(G)),    optional, intent(inout) :: TKE

! Arguments:
!  (in)    G            = the ocean grid structure.
!  (in)    GV           = The ocean's vertical grid structure.
!  (in)    h            = the layer thicknesses, in m or kg m-2.
!                         units of h are referred to as "H" below.
!  (in)    opacity_band = opacity in each band of penetrating shortwave
!                         radiation (1/H). The indicies are band, i, k.
!  (in)    nsw          = number of bands of penetrating shortwave radiation
!  (in)    j            = j-index to work on
!  (in)    dt           = time step (seconds)
!  (in)    H_limit_fluxes = if the total ocean depth is less than this, they
!                         are scaled away to avoid numerical instabilities. (H)
!                         This would not be necessary if a finite heat
!                         capacity mud-layer were added.
!  (in)    adjustAbsorptionProfile = if true, apply heating above the layers
!                         in which it should have occurred to get the correct
!                         mean depth (and potential energy change) of the
!                         shortwave that should be absorbed by each layer.
!  (in)    absorbAllSW  = if true, any shortwave radiation that hits the
!                         bottom is absorbed uniformly over the water column.
!  (inout) T            = layer potential/conservative temperatures (deg C)
!  (inout) Pen_SW_bnd   = penetrating shortwave heating in each band that
!                         hits the bottom and will be redistributed through
!                         the water column (units of K*H), size nsw x SZI_(G).

! These optional arguments apply when the bulk mixed layer is used
! but are unnecessary with other schemes.
!  (in,opt)    eps      = small thickness that must remain in each layer, and
!                         which will not be subject to heating (units of H)
!  (inout,opt) ksort    = density-sorted k-indicies
!  (in,opt)    htot     = total mixed layer thickness, in H
!  (inout,opt) Ttot     = depth integrated mixed layer temperature (units of K H)
!  (in,opt)    dSV_dT   = the partial derivative of specific volume with temperature, in m3 kg-1 K-1.
!  (inout,opt) TKE      = the TKE sink from mixing the heating throughout a layer, in J m-2.

  real, dimension(SZI_(G),SZK_(G)) :: &
    T_chg_above    ! A temperature change that will be applied to all the thick
                   ! layers above a given layer, in K.  This is only nonzero if
                   ! adjustAbsorptionProfile is true, in which case the net
                   ! change in the temperature of a layer is the sum of the
                   ! direct heating of that layer plus T_chg_above from all of
                   ! the layers below, plus any contribution from absorbing
                   ! radiation that hits the bottom.
  real, dimension(SZI_(G)) :: &
    h_heat, &      ! The thickness of the water column that will be heated by
                   ! any remaining shortwave radiation (H units).
    T_chg, &       ! The temperature change of thick layers due to the remaining
                   ! shortwave radiation and contributions from T_chg_above, in K.
    Pen_SW_rem     ! The sum across all wavelength bands of the penetrating shortwave
                   ! heating that hits the bottom and will be redistributed through
                   ! the water column (in units of K H)
  real :: SW_trans          ! fraction of shortwave radiation that is not
                            ! absorbed in a layer (nondimensional)
  real :: unabsorbed        ! fraction of the shortwave radiation that
                            ! is not absorbed because the layers are too thin
  real :: Ih_limit          ! inverse of the total depth at which the
                            ! surface fluxes start to be limited (1/H)
  real :: h_min_heat        ! minimum thickness layer that should get heated (H)
  real :: opt_depth         ! optical depth of a layer (non-dim)
  real :: exp_OD            ! exp(-opt_depth) (non-dim)
  real :: heat_bnd          ! heating due to absorption in the current
                            ! layer by the current band, including any piece that
                            ! is moved upward (K H units)
  real :: SWa               ! fraction of the absorbed shortwave that is
                            ! moved to layers above with adjustAbsorptionProfile (non-dim)
  real :: coSWa_frac        ! The fraction of SWa that is actually moved upward.
  real :: min_SW_heating    ! A minimum remaining shortwave heating rate that will be
                            ! simply absorbed in the next layer for computational
                            ! efficiency, instead of continuing to penetrate, in units
                            ! of K H s-1.  The default, 2.5e-11, is about 0.08 K m / century.
  real :: epsilon           ! A small thickness that must remain in each
                            ! layer, and which will not be subject to heating (units of H)
  real :: I_G_Earth
  real :: g_Hconv2
  logical :: SW_Remains     ! If true, some column has shortwave radiation that
                            ! was not entirely absorbed.
  logical :: TKE_calc       ! If true, calculate the implications to the
                            ! TKE budget of the shortwave heating.
  real :: C1_6, C1_60
  integer :: is, ie, nz, i, k, ks, n
  SW_Remains = .false.

  min_SW_heating = 2.5e-11

  h_min_heat = 2.0*GV%Angstrom + GV%H_subroundoff
  is = G%isc ; ie = G%iec ; nz = G%ke
  C1_6 = 1.0 / 6.0 ; C1_60 = 1.0 / 60.0

  TKE_calc = (present(TKE) .and. present(dSV_dT))
  g_Hconv2 = GV%g_Earth * GV%H_to_kg_m2**2

  h_heat(:) = 0.0
  if (present(htot)) then ; do i=is,ie ; h_heat(i) = htot(i) ; enddo ; endif

  ! Apply penetrating SW radiation to remaining parts of layers.
  ! Excessively thin layers are not heated to avoid runaway temps.
  do ks=1,nz ; do i=is,ie
    k = ks
    if (present(ksort)) then
      if (ksort(i,ks) <= 0) cycle
      k = ksort(i,ks)
    endif
    epsilon = 0.0 ; if (present(eps)) epsilon = eps(i,k)

    T_chg_above(i,k) = 0.0

    if (h(i,k) > 1.5*epsilon) then
      do n=1,nsw ; if (Pen_SW_bnd(n,i) > 0.0) then
        ! SW_trans is the SW that is transmitted THROUGH the layer
        opt_depth = h(i,k) * opacity_band(n,i,k)
        exp_OD = exp(-opt_depth)
        SW_trans = exp_OD

        ! Heating at a rate of less than 10-4 W m-2 = 10-3 K m / Century,
        ! and of the layer in question less than 1 K / Century, can be
        ! absorbed without further penetration.
        ! ###Make these numbers into parameters!
        if (nsw*Pen_SW_bnd(n,i)*SW_trans < &
            dt*min_SW_heating*min(GV%m_to_H,1e3*h(i,k)) ) SW_trans = 0.0

        Heat_bnd = Pen_SW_bnd(n,i) * (1.0 - SW_trans)
        if (adjustAbsorptionProfile .and. (h_heat(i) > 0.0)) then
          !   In this case, a fraction of the heating is applied to the
          ! overlying water so that the mean pressure at which the shortwave
          ! heating occurs is exactly what it would have been with a careful
          ! pressure-weighted averaging of the exponential heating profile,
          ! hence there should be no TKE budget requirements due to this
          ! layer.  Very clever, but this is also limited so that the
          ! water above is not heated at a faster rate than the layer
          ! actually being heated, i.e., SWA <= h_heat / (h_heat + h(i,k))
          ! and takes the energetics of the rest of the heating into account.
          ! (-RWH, ~7 years later.)
          if (opt_depth > 1e-5) then
            SWa = ((opt_depth + (opt_depth + 2.0)*exp_OD) - 2.0) / &
              ((opt_depth + opacity_band(n,i,k) * h_heat(i)) * &
               (1.0 - exp_OD))
          else
            ! Use Taylor series expansion of the expression above for a
            ! more accurate form with very small layer optical depths.
            SWa = h(i,k) * (opt_depth * (1.0 - opt_depth)) / &
              ((h_heat(i) + h(i,k)) * (6.0 - 3.0*opt_depth))
          endif
          coSWa_frac = 0.0
          if (SWa*(h_heat(i) + h(i,k)) > h_heat(i)) then
            coSWa_frac = (SWa*(h_heat(i) + h(i,k)) - h_heat(i) ) / &
                         (SWa*(h_heat(i) + h(i,k)))
            SWa = h_heat(i) / (h_heat(i) + h(i,k))
          endif

          T_chg_above(i,k) = T_chg_above(i,k) + (SWa * Heat_bnd) / h_heat(i)
          T(i,k) = T(i,k) + ((1.0 - SWa) * Heat_bnd) / h(i,k)
        else
          coSWa_frac = 1.0
          T(i,k) = T(i,k) + Pen_SW_bnd(n,i) * (1.0 - SW_trans) / h(i,k)
        endif

        if (TKE_calc) then
          if (opt_depth > 1e-2) then
            TKE(i,k) = TKE(i,k) - coSWa_frac*Heat_bnd*dSV_dT(i,k)* &
               (0.5*h(i,k)*g_Hconv2) * &
               (opt_depth*(1.0+exp_OD) - 2.0*(1.0-exp_OD)) / (opt_depth*(1.0-exp_OD))
          else
            ! Use Taylor series-derived approximation to the above expression
            ! that is well behaved and more accurate when opt_depth is small.
            TKE(i,k) = TKE(i,k) - coSWa_frac*Heat_bnd*dSV_dT(i,k)* &
               (0.5*h(i,k)*g_Hconv2) * &
               (C1_6*opt_depth * (1.0 - C1_60*opt_depth**2))
          endif
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
  enddo ; enddo ! i & k loops


! if (.not.absorbAllSW .and. .not.adjustAbsorptionProfile) return

  ! Unless modified, there is no temperature change due to fluxes from the bottom.
  do i=is,ie ; T_chg(i) = 0.0 ; enddo

  if (absorbAllSW) then
    ! If there is still shortwave radiation at this point, it could go into
    ! the bottom (with a bottom mud model), or it could be redistributed back
    ! through the water column.
    do i=is,ie
      Pen_SW_rem(i) = Pen_SW_bnd(1,i)
      do n=2,nsw ; Pen_SW_rem(i) = Pen_SW_rem(i) + Pen_SW_bnd(n,i) ; enddo
    enddo
    do i=is,ie ; if (Pen_SW_rem(i) > 0.0) SW_Remains = .true. ; enddo

    Ih_limit = 1.0 / H_limit_fluxes
    do i=is,ie ; if ((Pen_SW_rem(i) > 0.0) .and. (h_heat(i) > 0.0)) then
      if (h_heat(i)*Ih_limit >= 1.0) then
        T_chg(i) = Pen_SW_rem(i) / h_heat(i) ; unabsorbed = 0.0
      else
        T_chg(i) = Pen_SW_rem(i) * Ih_limit
        unabsorbed = 1.0 - h_heat(i)*Ih_limit
      endif
      do n=1,nsw ; Pen_SW_bnd(n,i) = unabsorbed * Pen_SW_bnd(n,i) ; enddo
    endif ; enddo
  endif ! absorbAllSW

  if (absorbAllSW .or. adjustAbsorptionProfile) then
    do ks=nz,1,-1 ; do i=is,ie
      k = ks
      if (present(ksort)) then
        if (ksort(i,ks) <= 0) cycle
        k = ksort(i,ks)
      endif

      if (T_chg(i) > 0.0) then
        ! Only layers greater than h_min_heat thick should get heated.
        if (h(i,k) >= 2.0*h_min_heat) then ; T(i,k) = T(i,k) + T_chg(i)
        elseif (h(i,k) > h_min_heat) then
          T(i,k) = T(i,k) + T_chg(i) * (2.0 - 2.0*h_min_heat/h(i,k))
        endif
      endif
      ! Increase the heating for layers above.
      T_chg(i) = T_chg(i) + T_chg_above(i,k)
    enddo ; enddo
    if (present(htot) .and. present(Ttot)) then
      do i=is,ie ; Ttot(i) = Ttot(i) + T_chg(i) * htot(i) ; enddo
    endif
  endif ! absorbAllSW .or. adjustAbsorptionProfile

end subroutine absorbRemainingSW


subroutine sumSWoverBands(G, GV, h, opacity_band, nsw, j, dt, &
                          H_limit_fluxes, absorbAllSW, iPen_SW_bnd, netPen)
! This subroutine calculates the total shortwave heat flux integrated over
! bands as a function of depth.  This routine is only called for computing
! buoyancy fluxes for use in KPP. This routine does not update the state.
  type(ocean_grid_type),                 intent(in)    :: G
  type(verticalGrid_type),               intent(in)    :: GV
  real, dimension(SZI_(G),SZK_(G)),      intent(in)    :: h
  real, dimension(:,:,:),                intent(in)    :: opacity_band
  integer,                               intent(in)    :: nsw
  integer,                               intent(in)    :: j
  real,                                  intent(in)    :: dt
  real,                                  intent(in)    :: H_limit_fluxes
  logical,                               intent(in)    :: absorbAllSW
  real, dimension(:,:),                  intent(in)    :: iPen_SW_bnd
  real, dimension(SZI_(G),SZK_(G)+1), intent(inout) :: netPen ! Units of K H

! Arguments:
!  (in)      G             = ocean grid structure
!  (in)      GV            = The ocean's vertical grid structure.
!  (in)      h             = layer thickness (units of m or kg/m^2);
!                            units of h are referred to as H below.
!  (in)      opacity_band  = opacity in each band of penetrating shortwave
!                            radiation, in m-1. The indicies are band, i, k.
!  (in)      nsw           = number of bands of penetrating shortwave radiation
!  (in)      j             = j-index to work on
!  (in)      dt            = time step (seconds)
!  (inout)   Pen_SW_bnd    = penetrating shortwave heating in each band that
!                            hits the bottom and will be redistributed through
!                            the water column (K H units); size nsw x SZI_(G).
!  (out)     netPen        = attenuated flux at interfaces, summed over bands (K H units)

  real :: h_heat(SZI_(G))     ! thickness of the water column that receives
                              ! remaining shortwave radiation, in H.
  real :: Pen_SW_rem(SZI_(G)) ! sum across all wavelength bands of the
                              ! penetrating shortwave heating that hits the bottom
                              ! and will be redistributed through the water column
                              ! (K H units)

  real, dimension(size(iPen_SW_bnd,1),size(iPen_SW_bnd,2)) :: Pen_SW_bnd
  real :: SW_trans        ! fraction of shortwave radiation not
                          ! absorbed in a layer (nondimensional)
  real :: unabsorbed      ! fraction of the shortwave radiation
                          ! not absorbed because the layers are too thin.
  real :: Ih_limit        ! inverse of the total depth at which the
                          ! surface fluxes start to be limited (1/H units)
  real :: h_min_heat      ! minimum thickness layer that should get heated (H units)
  real :: opt_depth       ! optical depth of a layer (non-dim)
  real :: exp_OD          ! exp(-opt_depth) (non-dim)
  logical :: SW_Remains   ! If true, some column has shortwave radiation that
                          ! was not entirely absorbed.

  integer :: is, ie, nz, i, k, ks, n
  SW_Remains = .false.

  h_min_heat = 2.0*GV%Angstrom + GV%H_subroundoff
  is = G%isc ; ie = G%iec ; nz = G%ke

  pen_SW_bnd(:,:) = iPen_SW_bnd(:,:)
  do i=is,ie ; h_heat(i) = 0.0 ; enddo
  netPen(:,1) = sum( pen_SW_bnd(:,:), dim=1 ) ! Surface interface

  ! Apply penetrating SW radiation to remaining parts of layers.
  ! Excessively thin layers are not heated to avoid runaway temps.
  do k=1,nz

    do i=is,ie
      netPen(i,k+1) = 0.

      if (h(i,k) > 0.0) then
        do n=1,nsw ; if (Pen_SW_bnd(n,i) > 0.0) then
          ! SW_trans is the SW that is transmitted THROUGH the layer
          opt_depth = h(i,k)*GV%H_to_m * opacity_band(n,i,k)
          exp_OD = exp(-opt_depth)
          SW_trans = exp_OD

          ! Heating at a rate of less than 10-4 W m-2 = 10-3 K m / Century,
          ! and of the layer in question less than 1 K / Century, can be
          ! absorbed without further penetration.
          if ((nsw*Pen_SW_bnd(n,i)*SW_trans < GV%m_to_H*2.5e-11*dt) .and. &
              (nsw*Pen_SW_bnd(n,i)*SW_trans < h(i,k)*dt*2.5e-8)) &
            SW_trans = 0.0

          Pen_SW_bnd(n,i) = Pen_SW_bnd(n,i) * SW_trans
          netPen(i,k+1)   = netPen(i,k+1) + Pen_SW_bnd(n,i)
        endif ; enddo
      endif ! h(i,k) > 0.0

      ! Add to the accumulated thickness above that could be heated.
      ! Only layers greater than h_min_heat thick should get heated.
      if (h(i,k) >= 2.0*h_min_heat) then
        h_heat(i) = h_heat(i) + h(i,k)
      elseif (h(i,k) > h_min_heat) then
        h_heat(i) = h_heat(i) + (2.0*h(i,k) - 2.0*h_min_heat)
      endif
    enddo ! i loop
  enddo ! k loop

  if (absorbAllSW) then

    ! If there is still shortwave radiation at this point, it could go into
    ! the bottom (with a bottom mud model), or it could be redistributed back
    ! through the water column.
    do i=is,ie
      Pen_SW_rem(i) = Pen_SW_bnd(1,i)
      do n=2,nsw ; Pen_SW_rem(i) = Pen_SW_rem(i) + Pen_SW_bnd(n,i) ; enddo
    enddo
    do i=is,ie ; if (Pen_SW_rem(i) > 0.0) SW_Remains = .true. ; enddo

    Ih_limit = 1.0 / H_limit_fluxes
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

end module MOM_shortwave_abs
