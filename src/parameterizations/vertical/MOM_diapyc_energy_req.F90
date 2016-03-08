module MOM_diapyc_energy_req
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

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!* By Robert Hallberg, May 2015                                        *
!*                                                                     *
!*   This module calculates the energy requirements of mixing.         *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_diag_mediator, only : diag_ctrl
use MOM_checksums, only : hchksum
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_EOS, only : calculate_specific_vol_derivs

implicit none ; private

#include <MOM_memory.h>

public diapyc_energy_req_init, diapyc_energy_req_calc, diapyc_energy_req_test, diapyc_energy_req_end

type, public :: diapyc_energy_req_CS ; private
  logical :: initialized = .false. ! A variable that is here because empty
                                   ! structures are not permitted by some compilers.
end type diapyc_energy_req_CS

contains

subroutine diapyc_energy_req_test(h_3d, dt, tv, G)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h_3d
  type(thermo_var_ptrs),                 intent(inout) :: tv
  real,                                  intent(in)    :: dt
  type(ocean_grid_type),                 intent(in)    :: G
! Arguments: h_3d -  Layer thickness before entrainment, in m or kg m-2.
!  (in/out)  tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.

  real, dimension(SZK_(G)) :: &
    T0, S0, &   ! T0 & S0 are columns of initial temperatures and salinities, in degC and g/kg.
    h_col       ! h_col is a column of thicknesses h at tracer points, in H (m or kg m-2).
  real, dimension(SZK_(G)+1) :: &
    Kd, &       ! A column of diapycnal diffusivities at interfaces, in m2 s-1.
    h_top, h_bot ! Distances from the top or bottom, in H.
  real :: ustar, absf, htot
  real :: energy_Kd ! The energy used by diapycnal mixing in W m-2.
  real :: tmp1  ! A temporary array.
  integer :: i, j, k, is, ie, js, je, nz, itt
  logical :: surface_BL, bottom_BL
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

!$OMP do
  do j=js,je ; do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
    htot = 0.0 ; h_top(1) = 0.0
    do k=1,nz
      T0(k) = tv%T(i,j,k) ; S0(k) = tv%S(i,j,k)
      h_col(k) = h_3d(i,j,k)
      h_top(K+1) = h_top(K) + h_col(k)
    enddo
    htot = h_top(nz+1)
    h_bot(nz+1) = 0.0
    do k=nz,1,-1
      h_bot(K) = h_bot(K+1) + h_col(k)
    enddo

    ustar = 0.01 ! Change this to being an input parameter?
    absf = 0.25*((abs(G%CoriolisBu(I-1,J-1)) + abs(G%CoriolisBu(I,J))) + &
                 (abs(G%CoriolisBu(I-1,J-1)) + abs(G%CoriolisBu(I,J))))
    Kd(1) = 0.0 ; Kd(nz+1) = 0.0
    do K=2,nz
      tmp1 = h_top(K) * h_bot(K) * G%GV%H_to_m
      Kd(k) = ustar * 0.41 * (tmp1*ustar) / (absf*tmp1 + htot*ustar)
    enddo
    call diapyc_energy_req_calc(h_col, T0, S0, Kd, energy_Kd, dt, tv, G)
  endif ; enddo ; enddo

end subroutine diapyc_energy_req_test

subroutine diapyc_energy_req_calc(h_in, T_in, S_in, Kd, energy_Kd, dt, tv, G)
  type(ocean_grid_type),                 intent(in)    :: G
  real, dimension(SZK_(G)),              intent(in)    :: h_in, T_in, S_in
  real, dimension(SZK_(G)+1),            intent(in)    :: Kd
  real,                                  intent(in)    :: dt
  real,                                  intent(out)   :: energy_Kd
  type(thermo_var_ptrs),                 intent(inout) :: tv
! Arguments: h_in -  Layer thickness before entrainment, in m or kg m-2.
!  (in)      T_in - The layer temperatures, in degC.
!  (in)      S_in - The layer salinities, in g kg-1.
!  (in)      Kd - The interfaces diapycnal diffusivities, in m2 s-1.
!  (in/out)  tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      dt - The amount of time covered by this call, in s.
!  (out)     energy_Kd - The column-integrated rate of energy consumption
!                 by diapycnal diffusion, in W m-2.
!  (in)      G - The ocean's grid structure.

!   This subroutine uses a substantially refactored tridiagonal equation for
! diapycnal mixing of temperature and salinity to estimate the potential energy
! change due to diapycnal mixing in a column of water.  It does this estimate
! 4 different ways, all of which should be equivalent, but reports only one.
! The various estimates are taken because they will later be used as templates
! for other bits of code.

  real, dimension(SZK_(G)) :: &
    p_lay, &    ! Average pressure of a layer, in Pa.
    dSV_dT, &   ! Partial derivative of specific volume with temperature, in m3 kg-1 K-1.
    dSV_dS, &   ! Partial derivative of specific volume with salinity, in m3 kg-1 / (g kg-1).
    T0, S0, &
    Te, Se, &
    dTe, dSe, & ! Running (1-way) estimates of temperature and salinity change.
    dT_to_dPE, & ! Partial derivative of column potential energy with the temperature
    dS_to_dPE, & ! and salinity changes within a layer, in J m-2 K-1 and J m-2 / (g kg-1).
    b_denom_1, & ! The first term in the denominator of b1, in m or kg m-2.
    c1, &       ! c1 is used by the tridiagonal solver, ND.
    h_tr        ! h_tr is h at tracer points with a h_neglect added to
                ! ensure positive definiteness, in m or kg m-2.
  real, dimension(SZK_(G)+1) :: &
    pres, &     ! Interface pressures in Pa.
    dT_to_dPE_a, &
    dS_to_dPE_a, &
    dT_to_dPE_b, &
    dS_to_dPE_b, &
    PE_chg_k, & ! The integrated potential energy change within a timestep due
                ! to the diffusivity at interface K, in J m-2.
    PEchg, &
    Kddt_h      ! The diapycnal diffusivity times a timestep divided by the
                ! average thicknesses around a layer, in m or kg m-2.
  real :: &
    b1              ! b1 is used by the tridiagonal solver, in m-1 or m2 kg-1.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in m.
  real :: dTe_term, dSe_term
  real :: Kddt_h_guess
  real :: htot
  real :: dT_k, dT_km1, dS_k, dS_km1  ! Temporary arrays
  real :: b1Kd                        ! Temporary arrays
  real :: Kd_rat, Kdr_denom, I_Kdr_denom  ! Temporary arrays
  real :: dTe_t2, dSe_t2, dT_km1_t2, dS_km1_t2, dT_k_t2, dS_k_t2

  ! The following are a bunch of diagnostic arrays for debugging purposes.
  real, dimension(SZK_(G)) :: &
    Ta, Sa
  real, dimension(SZK_(G)+1) :: &
    dPEa_dKd, dPEa_dKd_est, dPEa_dKd_err, dPEa_dKd_trunc, dPEa_dKd_err_norm, &
    dPEb_dKd, dPEb_dKd_est, dPEb_dKd_err, dPEb_dKd_trunc, dPEb_dKd_err_norm
  real :: PE_chg_tot1A, PE_chg_tot2A, T_chg_totA
  real :: PE_chg_tot1B, PE_chg_tot2B, T_chg_totB
  real, dimension(SZK_(G)+1)  :: dPEchg_dKd
  real :: PE_chg(6)
  real, dimension(6) :: dT_k_itt, dS_k_itt, dT_km1_itt, dS_km1_itt

  real :: I_G_Earth
  real :: dKd_rat_dKd, ddT_k_dKd, ddS_k_dKd, ddT_km1_dKd, ddS_km1_dKd
  integer :: k, nz, itt, max_itt
  logical :: surface_BL, bottom_BL, debug
  nz = G%ke
  h_neglect = G%GV%H_subroundoff

  I_G_Earth = 1.0 / G%G_earth
  surface_BL = .true. ; bottom_BL = .true. ; debug = .true.

  dPEa_dKd(:) = 0.0 ; dPEa_dKd_est(:) = 0.0 ; dPEa_dKd_err(:) = 0.0 ; dPEa_dKd_err_norm(:) = 0.0 ; dPEa_dKd_trunc(:) = 0.0
  dPEb_dKd(:) = 0.0 ; dPEb_dKd_est(:) = 0.0 ; dPEb_dKd_err(:) = 0.0 ; dPEb_dKd_err_norm(:) = 0.0 ; dPEb_dKd_trunc(:) = 0.0

  htot = 0.0 ; pres(1) = 0.0
  do k=1,nz
    T0(k) = T_in(k) ; S0(k) = S_in(k)
    h_tr(k) = h_in(k)
    htot = htot + h_tr(k)
    pres(K+1) = pres(K) + G%g_Earth * G%GV%H_to_kg_m2 * h_tr(k)
    p_lay(k) = 0.5*(pres(K) + pres(K+1))
  enddo
  do k=1,nz
    h_tr(k) = max(h_tr(k), 1e-15*htot)
  enddo

  ! Introduce a diffusive flux variable, Kddt_h(K) = ea(k) = eb(k-1)

  Kddt_h(1) = 0.0 ; Kddt_h(nz+1) = 0.0
  do K=2,nz
    Kddt_h(K) = min((G%GV%m_to_H**2*dt)*Kd(k) / (0.5*(h_tr(k-1) + h_tr(k))),1e3*htot)
  enddo

  ! Solve the tridiagonal equations for new temperatures.

  call calculate_specific_vol_derivs(T0, S0, p_lay, dSV_dT, dSV_dS, 1, nz, tv%eqn_of_state)

  do k=1,nz
    dT_to_dPE(k) = 0.5 * I_G_Earth *(pres(K+1)**2 - pres(K)**2) * dSV_dT(k)
    dS_to_dPE(k) = 0.5 * I_G_Earth *(pres(K+1)**2 - pres(K)**2) * dSV_dS(k)
  enddo

  PE_chg_k(1) = 0.0 ; PE_chg_k(nz+1) = 0.0
  PEchg(:) = 0.0 ; dPEchg_dKd(:) = 0.0

  if (surface_BL) then  ! This version is appropriate for a surface boundary layer.

    b_denom_1(1) = h_tr(1)
    dT_to_dPE_a(1) = dT_to_dPE(1)
    dS_to_dPE_a(1) = dS_to_dPE(1)
    do K=2,nz
      ! At this point, Kddt_h(K) will be unknown because its value may depend
      ! on how much energy is available.

      ! Precalculate some temporary expressions that are independent of Kddt_h_guess.
      if (K==2) then
        dT_km1_t2 = (T0(k)-T0(k-1))
        dS_km1_t2 = (S0(k)-S0(k-1))
        dTe_t2 = 0.0 ; dSe_t2 = 0.0
      else
        dTe_t2 = Kddt_h(K-1) * ((T0(k-2) - T0(k-1)) + dTe(k-2))
        dSe_t2 = Kddt_h(K-1) * ((S0(k-2) - S0(k-1)) + dSe(k-2))
        dT_km1_t2 = (T0(k)-T0(k-1)) - &
              (Kddt_h(K-1) / b_denom_1(k-1)) * ((T0(k-2) - T0(k-1)) + dTe(k-2))
        dS_km1_t2 = (S0(k)-S0(k-1)) - &
              (Kddt_h(K-1) / b_denom_1(k-1)) * ((S0(k-2) - S0(k-1)) + dSe(k-2))
      endif
      dTe_term = dTe_t2 + b_denom_1(k-1) * (T0(k-1)-T0(k))
      dSe_term = dSe_t2 + b_denom_1(k-1) * (S0(k-1)-S0(k))


      ! Find the energy change due to a guess at the strength of diffusion at interface K.

      Kddt_h_guess = Kddt_h(K)
      call find_PE_chg(Kddt_h_guess, h_tr(k), b_denom_1(k-1), &
                       dTe_term, dSe_term, dT_km1_t2, dS_km1_t2, &
                       dT_to_dPE(k), dS_to_dPE(k), dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), &
                       PE_chg_k(k), dPEa_dKd(k))

      if (debug) then
        do itt=1,5
          Kddt_h_guess = (1.0+0.01*(itt-3))*Kddt_h(K)

          call find_PE_chg(Kddt_h_guess, h_tr(k), b_denom_1(k-1), &
                           dTe_term, dSe_term, dT_km1_t2, dS_km1_t2, &
                           dT_to_dPE(k), dS_to_dPE(k), dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), &
                           PE_chg=PE_chg(itt))
        enddo
        ! Compare with a 4th-order finite difference estimate.
        dPEa_dKd_est(k) = (4.0*(PE_chg(4)-Pe_chg(2))/(0.02*Kddt_h(K)) - &
                               (PE_chg(5)-Pe_chg(1))/(0.04*Kddt_h(K))) / 3.0
        dPEa_dKd_trunc(k) = (PE_chg(4)-Pe_chg(2))/(0.02*Kddt_h(K)) - &
                            (PE_chg(5)-Pe_chg(1))/(0.04*Kddt_h(K))
        dPEa_dKd_err(k) = (dPEa_dKd_est(k) - dPEa_dKd(k))
        dPEa_dKd_err_norm(k) = (dPEa_dKd_est(k) - dPEa_dKd(k)) / &
                              (abs(dPEa_dKd_est(k)) + abs(dPEa_dKd(k)) + 1e-100)
      endif

      !   At this point, the final value of Kddt_h(K) is known, so the estimated
      ! properties for layer k-1 can be calculated.

      b1 = 1.0 / (b_denom_1(k-1) + Kddt_h(K))
      c1(K) = Kddt_h(K) * b1
      if (k==2) then
        Te(1) = b1*(h_tr(1)*T0(1))
        Se(1) = b1*(h_tr(1)*S0(1))
      else
        Te(k-1) = b1 * (h_tr(k-1) * T0(k-1) + Kddt_h(K-1) * Te(k-2))
        Se(k-1) = b1 * (h_tr(k-1) * S0(k-1) + Kddt_h(K-1) * Se(k-2))
      endif
      dTe(k-1) = b1 * ( Kddt_h(K)*(T0(k)-T0(k-1)) + dTe_t2 )
      dSe(k-1) = b1 * ( Kddt_h(K)*(S0(k)-S0(k-1)) + dSe_t2 )

      b_denom_1(k) = h_tr(k) + (b_denom_1(k-1) * b1) * Kddt_h(K)
      dT_to_dPE_a(k) = dT_to_dPE(k) + c1(K)*dT_to_dPE_a(k-1)
      dS_to_dPE_a(k) = dS_to_dPE(k) + c1(K)*dS_to_dPE_a(k-1)

    enddo

    b1 = 1.0 / (b_denom_1(nz))
    Te(nz) = b1 * (h_tr(nz) * T0(nz) + Kddt_h(nz) * Te(nz-1))
    Se(nz) = b1 * (h_tr(nz) * S0(nz) + Kddt_h(nz) * Se(nz-1))
    dTe(nz) = b1 * Kddt_h(nz) * ((T0(nz-1)-T0(nz)) + dTe(nz-1))
    dSe(nz) = b1 * Kddt_h(nz) * ((S0(nz-1)-S0(nz)) + dSe(nz-1))

    do k=nz-1,1,-1
      Te(k) = Te(k) + c1(K+1)*Te(k+1)
      Se(k) = Se(k) + c1(K+1)*Se(k+1)
    enddo

    if (debug) then
      do k=1,nz ; Ta(k) = Te(k) ; Sa(k) = Se(k) ; enddo
      PE_chg_tot1A = 0.0 ; PE_chg_tot2A = 0.0 ; T_chg_totA = 0.0
      do k=1,nz
        PE_chg_tot1A = PE_chg_tot1A + (dT_to_dPE(k) * (Te(k) - T0(k)) + &
                                      dS_to_dPE(k) * (Se(k) - S0(k)))
        T_chg_totA = T_chg_totA + h_tr(k) * (Te(k) - T0(k))
      enddo
      do K=2,nz
        PE_chg_tot2A = PE_chg_tot2A + PE_chg_k(K)
      enddo
    endif

  endif

  if (bottom_BL) then  ! This version is appropriate for a bottom boundary layer.

    b_denom_1(nz) = h_tr(nz)
    dT_to_dPE_b(nz) = dT_to_dPE(nz)
    dS_to_dPE_b(nz) = dS_to_dPE(nz)
    do K=nz,2,-1  ! Loop over interior interfaces.
      ! At this point, Kddt_h(K) will be unknown because its value may depend
      ! on how much energy is available.

      ! Precalculate some temporary expressions that are independent of Kddt_h_guess.
      if (K==nz) then
        dT_k_t2 = (T0(k-1)-T0(k))
        dS_k_t2 = (S0(k-1)-S0(k))
        dTe_t2 = 0.0 ; dSe_t2 = 0.0
      else
        dTe_t2 = Kddt_h(K+1) * ((T0(k+1) - T0(k)) + dTe(k+1))
        dSe_t2 = Kddt_h(K+1) * ((S0(k+1) - S0(k)) + dSe(k+1))
        dT_k_t2 = (T0(k-1)-T0(k)) - &
                (Kddt_h(k+1)/ b_denom_1(k)) * ((T0(k+1) - T0(k)) + dTe(k+1))
        dS_k_t2 = (S0(k-1)-S0(k)) - &
                (Kddt_h(k+1)/ b_denom_1(k)) * ((S0(k+1) - S0(k)) + dSe(k+1))
      endif
      dTe_term = dTe_t2 + b_denom_1(k) * (T0(k)-T0(k-1))
      dSe_term = dSe_t2 + b_denom_1(k) * (S0(k)-S0(k-1))

      ! Find the energy change due to a guess at the strength of diffusion at interface K.

      Kddt_h_guess = Kddt_h(K)
      call find_PE_chg(Kddt_h_guess, h_tr(k-1), b_denom_1(k), &
                       dTe_term, dSe_term, dT_k_t2, dS_k_t2, &
                       dT_to_dPE(k-1), dS_to_dPE(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                       PE_chg_k(K), dPEb_dKd(K))


      if (debug) then
        ! Compare with a 4th-order finite difference estimate.
        do itt=1,5
          Kddt_h_guess = (1.0+0.01*(itt-3))*Kddt_h(K)

          call find_PE_chg(Kddt_h_guess, h_tr(k-1), b_denom_1(k), &
                           dTe_term, dSe_term, dT_k_t2, dS_k_t2, &
                           dT_to_dPE(k-1), dS_to_dPE(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                           PE_chg=PE_chg(itt))
        enddo

        dPEb_dKd_est(k) = (4.0*(PE_chg(4)-Pe_chg(2))/(0.02*Kddt_h(K)) - &
                               (PE_chg(5)-Pe_chg(1))/(0.04*Kddt_h(K))) / 3.0
        dPEb_dKd_trunc(k) = (PE_chg(4)-Pe_chg(2))/(0.02*Kddt_h(K)) - &
                            (PE_chg(5)-Pe_chg(1))/(0.04*Kddt_h(K))
        dPEb_dKd_err(k) = (dPEb_dKd_est(k) - dPEb_dKd(k))
        dPEb_dKd_err_norm(k) = (dPEb_dKd_est(k) - dPEb_dKd(k)) / &
                              (abs(dPEb_dKd_est(k)) + abs(dPEb_dKd(k)) + 1e-100)
      endif

      !   At this point, the final value of Kddt_h(K) is known, so the estimated
      ! properties for layer k can be calculated.

      b1 = 1.0 / (b_denom_1(k) + Kddt_h(K))
      c1(K) = Kddt_h(K) * b1
      if (k==nz) then
        Te(nz) = b1* (h_tr(nz)*T0(nz))
        Se(nz) = b1* (h_tr(nz)*S0(nz))
      else
        Te(k) = b1 * (h_tr(k) * T0(k) + Kddt_h(K+1) * Te(k+1))
        Se(k) = b1 * (h_tr(k) * S0(k) + Kddt_h(k+1) * Se(k+1))
      endif
      dTe(k) = b1 * ( Kddt_h(K)*(T0(k-1)-T0(k)) + dTe_t2 )
      dSe(k) = b1 * ( Kddt_h(K)*(S0(k-1)-S0(k)) + dSe_t2 )

      b_denom_1(k-1) = h_tr(k-1) + (b_denom_1(k) * b1) * Kddt_h(K)
      dT_to_dPE_b(k-1) = dT_to_dPE(k-1) + c1(K)*dT_to_dPE_b(k)
      dS_to_dPE_b(k-1) = dS_to_dPE(k-1) + c1(K)*dS_to_dPE_b(k)

    enddo

    b1 = 1.0 / (b_denom_1(1))
    Te(1) = b1 * (h_tr(1) * T0(1) + Kddt_h(2) * Te(2))
    Se(1) = b1 * (h_tr(1) * S0(1) + Kddt_h(2) * Se(2))
    dTe(1) = b1 * Kddt_h(2) * ((T0(2)-T0(1)) + dTe(2))
    dSe(1) = b1 * Kddt_h(2) * ((S0(2)-S0(1)) + dSe(2))

    do k=2,nz
      Te(k) = Te(k) + c1(K)*Te(k-1)
      Se(k) = Se(k) + c1(K)*Se(k-1)
    enddo

    if (debug) then
      PE_chg_tot1B = 0.0 ; PE_chg_tot2B = 0.0 ; T_chg_totB = 0.0
      do k=1,nz
        PE_chg_tot1B = PE_chg_tot1B + (dT_to_dPE(k) * (Te(k) - T0(k)) + &
                                       dS_to_dPE(k) * (Se(k) - S0(k)))
        T_chg_totB = T_chg_totB + h_tr(k) * (Te(k) - T0(k))
      enddo
      do K=2,nz
        PE_chg_tot2B = PE_chg_tot2B + PE_chg_k(K)
      enddo
    endif

  endif

  energy_Kd = 0.0 ; do K=2,nz ; energy_Kd = energy_Kd + PE_chg_k(K) ; enddo
  energy_Kd = energy_Kd / dt
  K=nz

end subroutine diapyc_energy_req_calc

subroutine find_PE_chg(Kddt_h, h_k, b_den_1, dTe_term, dSe_term, &
                       dT_km1_t2, dS_km1_t2, dT_to_dPE_k, dS_to_dPE_k, &
                       dT_to_dPEa, dS_to_dPEa, PE_chg, dPE_chg_dKd)
  real, intent(in)  :: Kddt_h, h_k, b_den_1, dTe_term, dSe_term
  real, intent(in)  :: dT_km1_t2, dS_km1_t2, dT_to_dPE_k, dS_to_dPE_k
  real, intent(in)  :: dT_to_dPEa, dS_to_dPEa
  real, optional, intent(out) :: PE_chg
  real, optional, intent(out) :: dPE_chg_dKd

!   This subroutine determines the total potential energy change due to mixing
! at an interface, including all of the implicit effects of the prescribed
! mixing at interfaces above.  Everything here is derived by careful manipulation
! of the robust tridiagonal solvers used for tracers by MOM6.
!   The comments describing these arguments are for a downward mixing pass, but
! this routine can also be used for an upward pass with the sense of direction
! reversed.
!
! Arguments: Kddt_h - The diffusivity at an interface times the time step and
!                     divided by the average of the thicknesses around the
!                     interface, in units of H (m or kg-2).
!  (in)      h_k - The thickness of the layer below the interface, in H.
!  (in)      b_den_1 - The first first term in the denominator of the pivot
!                  for the tridiagonal solver, given by h_k plus a term that
!                  is a fraction (determined from the tridiagonal solver) of
!                  Kddt_h for the interface above, in H.
!  (in)      dTe_term - A diffusivity-independent term related to the
!                  temperature change in the layer below the interface, in K H.
!  (in)      dSe_term - A diffusivity-independent term related to the
!                  salinity change in the layer below the interface, in ppt H.
!  (in)      dT_km1_t2 - A diffusivity-independent term related to the
!                  temperature change in the layer above the interface, in K.
!  (in)      dS_km1_t2 - A diffusivity-independent term related to the
!                  salinity change in the layer above the interface, in ppt.
!  (in)      dT_to_dPE_k - A factor (1/2*pres*mass_lay*dSpec_vol/dT) relating
!                  a layer's temperature change to the change in column
!                  potential energy, in J m-2 K-1.
!  (in)      dS_to_dPE_k - A factor (1/2*pres*mass_lay*dSpec_vol/dS) relating
!                  a layer's salinity change to the change in column
!                  potential energy, in J m-2 ppt-1.
!  (in)      dT_to_dPEa - A factor (1/2*pres*mass_lay*dSpec_vol/dT) relating
!                  a layer's temperature change to the change in column
!                  potential energy, including all implicit diffusive changes
!                  in the temperatures of all the layers above, in J m-2 K-1.
!  (in)      dS_to_dPEa - A factor (1/2*pres*mass_lay*dSpec_vol/dS) relating
!                  a layer's salinity change to the change in column
!                  potential energy, including all implicit diffusive changes
!                  in the salinities of all the layers above, in J m-2 ppt-1.
!  (out,opt) PE_chg - The change in column potential energy from applying
!                  Kddt_h at the present interface, in J m-2.
!  (out,opt) dPE_chg_dKd - The partial derivative of PE_chg with Kddt_h,
!                  in units of J m-2 H-1.

  ! b_den_1 - The first term in the denominator of b1, in m or kg m-2.
  real :: b1            ! b1 is used by the tridiagonal solver, in H-1.
  real :: b1Kd          ! Temporary array (nondim.)
  real :: dT_k, dT_km1  ! Temporary arrays in K.
  real :: dS_k, dS_km1  ! Temporary arrays in ppt.
  real :: I_Kr_denom, dKr_dKd   ! Temporary arrays in H-2 and nondim.
  real :: ddT_k_dKd, ddT_km1_dKd ! Temporary arrays in K H-1.
  real :: ddS_k_dKd, ddS_km1_dKd ! Temporary arrays in ppt H-1.

  b1 = 1.0 / (b_den_1 + Kddt_h)
  b1Kd = Kddt_h*b1

  ! Start with the temperature change in layer k-1 due to the diffusivity at
  ! interface K without considering the effects of changes in layer k.

  ! Calculate the change in PE due to the diffusion at interface K
  ! if Kddt_h(K+1) = 0.
  I_Kr_denom = 1.0 / (h_k*b_den_1 + (b_den_1 + h_k)*Kddt_h)

  dT_k = (Kddt_h*I_Kr_denom) * dTe_term
  dS_k = (Kddt_h*I_Kr_denom) * dSe_term

  if (present(PE_chg)) then
    ! Find the change in energy due to diffusion with strength Kddt_h at this interface.
    ! Increment the temperature changes in layer k-1 due the changes in layer k.
    dT_km1 = b1Kd * ( dT_k + dT_km1_t2 )
    dS_km1 = b1Kd * ( dS_k + dS_km1_t2 )

    PE_chg = (dT_to_dPE_k * dT_k + dT_to_dPEa * dT_km1) + &
             (dS_to_dPE_k * dS_k + dS_to_dPEa * dS_km1)
  endif

  if (present(dPE_chg_dKd)) then
    ! Find the derivatives of the temperature and salinity changes with Kddt_h.
    dKr_dKd = (h_k*b_den_1) * I_Kr_denom**2

    ddT_k_dKd = dKr_dKd * dTe_term
    ddS_k_dKd = dKr_dKd * dSe_term
    ddT_km1_dKd = (b1**2 * b_den_1) * ( dT_k + dT_km1_t2 ) + b1Kd * ddT_k_dKd
    ddS_km1_dKd = (b1**2 * b_den_1) * ( dS_k + dS_km1_t2 ) + b1Kd * ddS_k_dKd

    ! Calculate the partial derivative of Pe_chg with Kddt_h.
    dPE_chg_dKd = (dT_to_dPE_k * ddT_k_dKd + dT_to_dPEa * ddT_km1_dKd) + &
                  (dS_to_dPE_k * ddS_k_dKd + dS_to_dPEa * ddS_km1_dKd)
  endif

end subroutine find_PE_chg

subroutine diapyc_energy_req_init(G, param_file, CS)
  type(ocean_grid_type),        intent(in) :: G
  type(param_file_type),        intent(in) :: param_file
  type(diapyc_energy_req_CS),   pointer    :: CS
! Arguments: param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  Reg - A pointer that is set to point to the tracer registry.
  integer, save :: init_calls = 0
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_diapyc_energy_req" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (.not.associated(CS)) then ; allocate(CS)
  else ; return ; endif

  CS%initialized = .true.

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")

end subroutine diapyc_energy_req_init

subroutine diapyc_energy_req_end(CS)
  type(diapyc_energy_req_CS), pointer :: CS
  if (associated(CS)) deallocate(CS)
end subroutine diapyc_energy_req_end

end module MOM_diapyc_energy_req
