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

use MOM_diag_mediator, only : diag_ctrl, Time_type, post_data_1d_k, register_diag_field
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_specific_vol_derivs, calculate_density

implicit none ; private

public diapyc_energy_req_init, diapyc_energy_req_calc, diapyc_energy_req_test, diapyc_energy_req_end

type, public :: diapyc_energy_req_CS ; private
  logical :: initialized = .false. ! A variable that is here because empty
                               ! structures are not permitted by some compilers.
  real :: test_Kh_scaling      ! A scaling factor for the diapycnal diffusivity.
  real :: ColHt_scaling        ! A scaling factor for the column height change
                               ! correction term.
  logical :: use_test_Kh_profile   ! If true, use the internal test diffusivity
                               ! profile in place of any that might be passed
                               ! in as an argument.
  type(diag_ctrl), pointer :: diag ! Structure used to regulate timing of diagnostic output
  integer :: id_ERt=-1, id_ERb=-1, id_ERc=-1, id_ERh=-1, id_Kddt=-1, id_Kd=-1
  integer :: id_CHCt=-1, id_CHCb=-1, id_CHCc=-1, id_CHCh=-1
  integer :: id_T0=-1, id_Tf=-1, id_S0=-1, id_Sf=-1, id_N2_0=-1, id_N2_f=-1
  integer :: id_h=-1, id_zInt=-1
end type diapyc_energy_req_CS

contains

subroutine diapyc_energy_req_test(h_3d, dt, tv, G, GV, CS, Kd_int)
  type(ocean_grid_type),          intent(in)    :: G
  type(verticalGrid_type),        intent(in)    :: GV
  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke), &
                                  intent(in)    :: h_3d
  type(thermo_var_ptrs),          intent(inout) :: tv
  real,                           intent(in)    :: dt
  type(diapyc_energy_req_CS),     pointer       :: CS
  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke+1), &
                        optional, intent(in)    :: Kd_int

! Arguments: h_3d -  Layer thickness before entrainment, in m or kg m-2.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - This module's control structure
!  (in,opt)  Kd_int -  Interface diffusivities.

  real, dimension(GV%ke) :: &
    T0, S0, &   ! T0 & S0 are columns of initial temperatures and salinities, in degC and g/kg.
    h_col       ! h_col is a column of thicknesses h at tracer points, in H (m or kg m-2).
  real, dimension(GV%ke+1) :: &
    Kd, &       ! A column of diapycnal diffusivities at interfaces, in m2 s-1.
    h_top, h_bot ! Distances from the top or bottom, in H.
  real :: ustar, absf, htot
  real :: energy_Kd ! The energy used by diapycnal mixing in W m-2.
  real :: tmp1  ! A temporary array.
  integer :: i, j, k, is, ie, js, je, nz, itt
  logical :: may_print
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "diapyc_energy_req_test: "// &
         "Module must be initialized before it is used.")

!$OMP do
  do j=js,je ; do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
    if (present(Kd_int) .and. .not.CS%use_test_Kh_profile) then
      do k=1,nz+1 ; Kd(K) = CS%test_Kh_scaling*Kd_int(i,j,K) ; enddo
    else
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
        tmp1 = h_top(K) * h_bot(K) * GV%H_to_m
        Kd(K) = CS%test_Kh_scaling * &
                ustar * 0.41 * (tmp1*ustar) / (absf*tmp1 + htot*ustar)
      enddo
    endif
    may_print = is_root_PE() .and. (i==ie) .and. (j==je)
    call diapyc_energy_req_calc(h_col, T0, S0, Kd, energy_Kd, dt, tv, G, GV, &
                                may_print=may_print, CS=CS)
  endif ; enddo ; enddo

end subroutine diapyc_energy_req_test

subroutine diapyc_energy_req_calc(h_in, T_in, S_in, Kd, energy_Kd, dt, tv, &
                                  G, GV, may_print, CS)
  type(ocean_grid_type),    intent(in)    :: G
  type(verticalGrid_type),  intent(in)    :: GV
  real, dimension(GV%ke),   intent(in)    :: h_in, T_in, S_in
  real, dimension(GV%ke+1), intent(in)    :: Kd
  real,                     intent(in)    :: dt
  real,                     intent(out)   :: energy_Kd
  type(thermo_var_ptrs),    intent(inout) :: tv
  logical,        optional, intent(in)    :: may_print
  type(diapyc_energy_req_CS), optional, pointer :: CS
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
!  (in)      GV - The ocean's vertical grid structure.
!  (in,opt)  may_print - If present and true, write out diagnostics of energy use.
!  (in,opt)  CS - This module's control structure

!   This subroutine uses a substantially refactored tridiagonal equation for
! diapycnal mixing of temperature and salinity to estimate the potential energy
! change due to diapycnal mixing in a column of water.  It does this estimate
! 4 different ways, all of which should be equivalent, but reports only one.
! The various estimates are taken because they will later be used as templates
! for other bits of code.

  real, dimension(GV%ke) :: &
    p_lay, &    ! Average pressure of a layer, in Pa.
    dSV_dT, &   ! Partial derivative of specific volume with temperature, in m3 kg-1 K-1.
    dSV_dS, &   ! Partial derivative of specific volume with salinity, in m3 kg-1 / (g kg-1).
    T0, S0, &   ! Initial temperatures and salinities.
    Te, Se, &   ! Running incomplete estimates of the new temperatures and salinities.
    Te_a, Se_a, &   ! Running incomplete estimates of the new temperatures and salinities.
    Te_b, Se_b, &   ! Running incomplete estimates of the new temperatures and salinities.
    Tf, Sf, &   ! New final values of the temperatures and salinities.
    dTe, dSe, & ! Running (1-way) estimates of temperature and salinity change.
    dTe_a, dSe_a, & ! Running (1-way) estimates of temperature and salinity change.
    dTe_b, dSe_b, & ! Running (1-way) estimates of temperature and salinity change.
    Th_a, &     !  An effective temperature times a thickness in the layer above,
                !  including implicit mixing effects with other yet higher layers, in K H.
    Sh_a, &     !  An effective salinity times a thickness in the layer above,
                !  including implicit mixing effects with other yet higher layers, in K H.
    Th_b, &     !  An effective temperature times a thickness in the layer below,
                !  including implicit mixing effects with other yet lower layers, in K H.
    Sh_b, &     !  An effective salinity times a thickness in the layer below,
                !  including implicit mixing effects with other yet lower layers, in K H.
    dT_to_dPE, & ! Partial derivative of column potential energy with the temperature
    dS_to_dPE, & ! and salinity changes within a layer, in J m-2 K-1 and J m-2 / (g kg-1).
    dT_to_dColHt, & ! Partial derivatives of the total column height with the temperature
    dS_to_dColHt, & ! and salinity changes within a layer, in m K-1 and m ppt-1.
    dT_to_dColHt_a, & ! Partial derivatives of the total column height with the temperature
    dS_to_dColHt_a, & ! and salinity changes within a layer, including the implicit effects
                    ! of mixing with layers higher in the water colun, in m K-1 and m ppt-1.
    dT_to_dColHt_b, & ! Partial derivatives of the total column height with the temperature
    dS_to_dColHt_b, & ! and salinity changes within a layer, including the implicit effects
                    ! of mixing with layers lower in the water colun, in m K-1 and m ppt-1.
    dT_to_dPE_a, &  ! Partial derivatives of column potential energy with the temperature
    dS_to_dPE_a, &  ! and salinity changes within a layer, including the implicit effects
                    ! of mixing with layers higher in the water column, in
                    ! units of J m-2 K-1 and J m-2 ppt-1.
    dT_to_dPE_b, &  ! Partial derivatives of column potential energy with the temperature
    dS_to_dPE_b, &  ! and salinity changes within a layer, including the implicit effects
                    ! of mixing with layers lower in the water column, in
                    ! units of J m-2 K-1 and J m-2 ppt-1.
    hp_a, &     ! An effective pivot thickness of the layer including the effects
                ! of coupling with layers above, in H.  This is the first term
                ! in the denominator of b1 in a downward-oriented tridiagonal solver.
    hp_b, &     ! An effective pivot thickness of the layer including the effects
                ! of coupling with layers below, in H.  This is the first term
                ! in the denominator of b1 in an upward-oriented tridiagonal solver.
    c1_a, &     ! c1_a is used by a downward-oriented tridiagonal solver, ND.
    c1_b, &     ! c1_b is used by an upward-oriented tridiagonal solver, ND.
    h_tr        ! h_tr is h at tracer points with a h_neglect added to
                ! ensure positive definiteness, in m or kg m-2.
  real, dimension(GV%ke+1) :: &
    pres, &     ! Interface pressures in Pa.
    z_Int, &    ! Interface heights relative to the surface, in m.
    N2, &       ! An estimate of the buoyancy frequency in s-2.
    Kddt_h_a, & !
    Kddt_h_b, & !
    Kddt_h, &   ! The diapycnal diffusivity times a timestep divided by the
                ! average thicknesses around a layer, in m or kg m-2.
    Kd_so_far   ! The value of Kddt_h that has been applied already in
                ! calculating the energy changes, in m or kg m-2.
  real, dimension(GV%ke+1,4) :: &
    PE_chg_k, &     ! The integrated potential energy change within a timestep due
                    ! to the diffusivity at interface K for 4 different orders of
                    ! accumulating the diffusivities, in J m-2.
    ColHt_cor_k     ! The correction to the potential energy change due to
                    ! changes in the net column height, in J m-2.
  real :: &
    b1              ! b1 is used by the tridiagonal solver, in m-1 or m2 kg-1.
  real :: &
    I_b1            ! The inverse of b1, in m or kg m-2.
  real :: Kd0, dKd
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in m.
  real :: dTe_term, dSe_term
  real :: Kddt_h_guess
  real :: dMass     ! The mass per unit area within a layer, in kg m-2.
  real :: dPres     ! The hydrostatic pressure change across a layer, in Pa.
  real :: dMKE_max  ! The maximum amount of mean kinetic energy that could be
                    ! converted to turbulent kinetic energy if the velocity in
                    ! the layer below an interface were homogenized with all of
                    ! the water above the interface, in J m-2 = kg s-2.
  real :: rho_here  ! The in-situ density, in kg m-3.
  real :: PE_change, ColHt_cor
  real :: htot
  real :: dT_k, dT_km1, dS_k, dS_km1  ! Temporary arrays
  real :: b1Kd                        ! Temporary arrays
  real :: Kd_rat, Kdr_denom, I_Kdr_denom  ! Temporary arrays
  real :: dTe_t2, dSe_t2, dT_km1_t2, dS_km1_t2, dT_k_t2, dS_k_t2
  logical :: do_print

  ! The following are a bunch of diagnostic arrays for debugging purposes.
  real, dimension(GV%ke) :: &
    Ta, Sa, Tb, Sb
  real, dimension(GV%ke+1) :: &
    dPEa_dKd, dPEa_dKd_est, dPEa_dKd_err, dPEa_dKd_trunc, dPEa_dKd_err_norm, &
    dPEb_dKd, dPEb_dKd_est, dPEb_dKd_err, dPEb_dKd_trunc, dPEb_dKd_err_norm
  real :: PE_chg_tot1A, PE_chg_tot2A, T_chg_totA
  real :: PE_chg_tot1B, PE_chg_tot2B, T_chg_totB
  real :: PE_chg_tot1C, PE_chg_tot2C, T_chg_totC
  real :: PE_chg_tot1D, PE_chg_tot2D, T_chg_totD
  real, dimension(GV%ke+1)  :: dPEchg_dKd
  real :: PE_chg(6)
  real, dimension(6) :: dT_k_itt, dS_k_itt, dT_km1_itt, dS_km1_itt

  real :: I_G_Earth
  real :: dKd_rat_dKd, ddT_k_dKd, ddS_k_dKd, ddT_km1_dKd, ddS_km1_dKd
  integer :: k, nz, itt, max_itt, k_cent
  logical :: surface_BL, bottom_BL, central, halves, debug
  logical :: old_PE_calc
  nz = G%ke
  h_neglect = GV%H_subroundoff

  I_G_Earth = 1.0 / GV%g_Earth
  debug = .true.

  surface_BL = .true. ; bottom_BL = .true. ; halves = .true.
  central = .true. ; K_cent = nz/2

  do_print = .false. ; if (present(may_print) .and. present(CS)) do_print = may_print

  dPEa_dKd(:) = 0.0 ; dPEa_dKd_est(:) = 0.0 ; dPEa_dKd_err(:) = 0.0 ; dPEa_dKd_err_norm(:) = 0.0 ; dPEa_dKd_trunc(:) = 0.0
  dPEb_dKd(:) = 0.0 ; dPEb_dKd_est(:) = 0.0 ; dPEb_dKd_err(:) = 0.0 ; dPEb_dKd_err_norm(:) = 0.0 ; dPEb_dKd_trunc(:) = 0.0

  htot = 0.0 ; pres(1) = 0.0 ; Z_int(1) = 0.0
  do k=1,nz
    T0(k) = T_in(k) ; S0(k) = S_in(k)
    h_tr(k) = h_in(k)
    htot = htot + h_tr(k)
    pres(K+1) = pres(K) + GV%g_Earth * GV%H_to_kg_m2 * h_tr(k)
    p_lay(k) = 0.5*(pres(K) + pres(K+1))
    Z_int(K+1) = Z_int(K) - GV%H_to_m * h_tr(k)
  enddo
  do k=1,nz
    h_tr(k) = max(h_tr(k), 1e-15*htot)
  enddo

  ! Introduce a diffusive flux variable, Kddt_h(K) = ea(k) = eb(k-1)

  Kddt_h(1) = 0.0 ; Kddt_h(nz+1) = 0.0
  do K=2,nz
    Kddt_h(K) = min((GV%m_to_H**2*dt)*Kd(k) / (0.5*(h_tr(k-1) + h_tr(k))),1e3*htot)
  enddo

  ! Solve the tridiagonal equations for new temperatures.

  call calculate_specific_vol_derivs(T0, S0, p_lay, dSV_dT, dSV_dS, 1, nz, tv%eqn_of_state)

  do k=1,nz
    dMass = GV%H_to_kg_m2 * h_tr(k)
    dPres = GV%g_Earth * dMass
    dT_to_dPE(k) = (dMass * (pres(K) + 0.5*dPres)) * dSV_dT(k)
    dS_to_dPE(k) = (dMass * (pres(K) + 0.5*dPres)) * dSV_dS(k)
    dT_to_dColHt(k) = dMass * dSV_dT(k) * CS%ColHt_scaling
    dS_to_dColHt(k) = dMass * dSV_dS(k) * CS%ColHt_scaling
  enddo

!  PE_chg_k(1) = 0.0 ; PE_chg_k(nz+1) = 0.0
  ! PEchg(:) = 0.0
  PE_chg_k(:,:) = 0.0 ; ColHt_cor_k(:,:) = 0.0
  dPEchg_dKd(:) = 0.0

  if (surface_BL) then  ! This version is appropriate for a surface boundary layer.
    old_PE_calc = .false.

    ! Set up values appropriate for no diffusivity.
    do k=1,nz
      hp_a(k) = h_tr(k) ; hp_b(k) = h_tr(k)
      dT_to_dPE_a(k) = dT_to_dPE(k) ; dS_to_dPE_a(k) = dS_to_dPE(k)
      dT_to_dPE_b(k) = dT_to_dPE(k) ; dS_to_dPE_b(k) = dS_to_dPE(k)
      dT_to_dColHt_a(k) = dT_to_dColHt(k) ; dS_to_dColHt_a(k) = dS_to_dColHt(k)
      dT_to_dColHt_b(k) = dT_to_dColHt(k) ; dS_to_dColHt_b(k) = dS_to_dColHt(k)
    enddo

    do K=2,nz
      ! At this point, Kddt_h(K) will be unknown because its value may depend
      ! on how much energy is available.

      ! Precalculate some temporary expressions that are independent of Kddt_h_guess.
      if (old_PE_calc) then
        if (K==2) then
          dT_km1_t2 = (T0(k)-T0(k-1))
          dS_km1_t2 = (S0(k)-S0(k-1))
          dTe_t2 = 0.0 ; dSe_t2 = 0.0
        else
          dTe_t2 = Kddt_h(K-1) * ((T0(k-2) - T0(k-1)) + dTe(k-2))
          dSe_t2 = Kddt_h(K-1) * ((S0(k-2) - S0(k-1)) + dSe(k-2))
          dT_km1_t2 = (T0(k)-T0(k-1)) - &
                (Kddt_h(K-1) / hp_a(k-1)) * ((T0(k-2) - T0(k-1)) + dTe(k-2))
          dS_km1_t2 = (S0(k)-S0(k-1)) - &
                (Kddt_h(K-1) / hp_a(k-1)) * ((S0(k-2) - S0(k-1)) + dSe(k-2))
        endif
        dTe_term = dTe_t2 + hp_a(k-1) * (T0(k-1)-T0(k))
        dSe_term = dSe_t2 + hp_a(k-1) * (S0(k-1)-S0(k))
      else
        if (K<=2) then
          Th_a(k-1) = h_tr(k-1) * T0(k-1) ; Sh_a(k-1) = h_tr(k-1) * S0(k-1)
        else
          Th_a(k-1) = h_tr(k-1) * T0(k-1) + Kddt_h(K-1) * Te(k-2)
          Sh_a(k-1) = h_tr(k-1) * S0(k-1) + Kddt_h(K-1) * Se(k-2)
        endif
        Th_b(k) = h_tr(k) * T0(k) ; Sh_b(k) = h_tr(k) * S0(k)
      endif


      ! Find the energy change due to a guess at the strength of diffusion at interface K.

      Kddt_h_guess = Kddt_h(K)
      if (old_PE_calc) then
        call find_PE_chg_orig(Kddt_h_guess, h_tr(k), hp_a(k-1), &
                         dTe_term, dSe_term, dT_km1_t2, dS_km1_t2, &
                         dT_to_dPE(k), dS_to_dPE(k), dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), &
                         pres(K), dT_to_dColHt(k), dS_to_dColHt(k), &
                         dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                         PE_chg_k(k,1), dPEa_dKd(k))
      else
        call find_PE_chg(0.0, Kddt_h_guess, hp_a(k-1), hp_b(k), &
                         Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                         dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                         pres(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                         dT_to_dColHt_b(k), dS_to_dColHt_b(k), &
                         PE_chg=PE_chg_k(K,1), dPEc_dKd=dPEa_dKd(K), &
                         ColHt_cor=ColHt_cor_k(K,1))
      endif

      if (debug) then
        do itt=1,5
          Kddt_h_guess = (1.0+0.01*(itt-3))*Kddt_h(K)

          if (old_PE_calc) then
            call find_PE_chg_orig(Kddt_h_guess, h_tr(k), hp_a(k-1), &
                             dTe_term, dSe_term, dT_km1_t2, dS_km1_t2, &
                             dT_to_dPE(k), dS_to_dPE(k), dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), &
                             pres(K), dT_to_dColHt(k), dS_to_dColHt(k), &
                             dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                             PE_chg=PE_chg(itt))
          else
            call find_PE_chg(0.0, Kddt_h_guess, hp_a(k-1), hp_b(k), &
                             Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                             dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                             pres(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                             dT_to_dColHt_b(k), dS_to_dColHt_b(k), &
                             PE_chg=PE_chg(itt))
          endif
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

      b1 = 1.0 / (hp_a(k-1) + Kddt_h(K))
      c1_a(K) = Kddt_h(K) * b1
      if (k==2) then
        Te(1) = b1*(h_tr(1)*T0(1))
        Se(1) = b1*(h_tr(1)*S0(1))
      else
        Te(k-1) = b1 * (h_tr(k-1) * T0(k-1) + Kddt_h(K-1) * Te(k-2))
        Se(k-1) = b1 * (h_tr(k-1) * S0(k-1) + Kddt_h(K-1) * Se(k-2))
      endif
      if (old_PE_calc) then
        dTe(k-1) = b1 * ( Kddt_h(K)*(T0(k)-T0(k-1)) + dTe_t2 )
        dSe(k-1) = b1 * ( Kddt_h(K)*(S0(k)-S0(k-1)) + dSe_t2 )
      endif

      hp_a(k) = h_tr(k) + (hp_a(k-1) * b1) * Kddt_h(K)
      dT_to_dPE_a(k) = dT_to_dPE(k) + c1_a(K)*dT_to_dPE_a(k-1)
      dS_to_dPE_a(k) = dS_to_dPE(k) + c1_a(K)*dS_to_dPE_a(k-1)
      dT_to_dColHt_a(k) = dT_to_dColHt(k) + c1_a(K)*dT_to_dColHt_a(k-1)
      dS_to_dColHt_a(k) = dS_to_dColHt(k) + c1_a(K)*dS_to_dColHt_a(k-1)

    enddo

    b1 = 1.0 / (hp_a(nz))
    Tf(nz) = b1 * (h_tr(nz) * T0(nz) + Kddt_h(nz) * Te(nz-1))
    Sf(nz) = b1 * (h_tr(nz) * S0(nz) + Kddt_h(nz) * Se(nz-1))
    if (old_PE_calc) then
      dTe(nz) = b1 * Kddt_h(nz) * ((T0(nz-1)-T0(nz)) + dTe(nz-1))
      dSe(nz) = b1 * Kddt_h(nz) * ((S0(nz-1)-S0(nz)) + dSe(nz-1))
    endif

    do k=nz-1,1,-1
      Tf(k) = Te(k) + c1_a(K+1)*Tf(k+1)
      Sf(k) = Se(k) + c1_a(K+1)*Sf(k+1)
    enddo

    if (debug) then
      do k=1,nz ; Ta(k) = Tf(k) ; Sa(k) = Sf(k) ; enddo
      PE_chg_tot1A = 0.0 ; PE_chg_tot2A = 0.0 ; T_chg_totA = 0.0
      do k=1,nz
        PE_chg_tot1A = PE_chg_tot1A + (dT_to_dPE(k) * (Tf(k) - T0(k)) + &
                                       dS_to_dPE(k) * (Sf(k) - S0(k)))
        T_chg_totA = T_chg_totA + h_tr(k) * (Tf(k) - T0(k))
      enddo
      do K=2,nz
        PE_chg_tot2A = PE_chg_tot2A + (PE_chg_k(K,1) - ColHt_cor_k(K,1))
      enddo
    endif

  endif

  if (bottom_BL) then  ! This version is appropriate for a bottom boundary layer.
    old_PE_calc = .false.

    ! Set up values appropriate for no diffusivity.
    do k=1,nz
      hp_a(k) = h_tr(k) ; hp_b(k) = h_tr(k)
      dT_to_dPE_a(k) = dT_to_dPE(k) ; dS_to_dPE_a(k) = dS_to_dPE(k)
      dT_to_dPE_b(k) = dT_to_dPE(k) ; dS_to_dPE_b(k) = dS_to_dPE(k)
      dT_to_dColHt_a(k) = dT_to_dColHt(k) ; dS_to_dColHt_a(k) = dS_to_dColHt(k)
      dT_to_dColHt_b(k) = dT_to_dColHt(k) ; dS_to_dColHt_b(k) = dS_to_dColHt(k)
    enddo

    do K=nz,2,-1  ! Loop over interior interfaces.
      ! At this point, Kddt_h(K) will be unknown because its value may depend
      ! on how much energy is available.

      ! Precalculate some temporary expressions that are independent of Kddt_h_guess.
      if (old_PE_calc) then
        if (K==nz) then
          dT_k_t2 = (T0(k-1)-T0(k))
          dS_k_t2 = (S0(k-1)-S0(k))
          dTe_t2 = 0.0 ; dSe_t2 = 0.0
        else
          dTe_t2 = Kddt_h(K+1) * ((T0(k+1) - T0(k)) + dTe(k+1))
          dSe_t2 = Kddt_h(K+1) * ((S0(k+1) - S0(k)) + dSe(k+1))
          dT_k_t2 = (T0(k-1)-T0(k)) - &
                  (Kddt_h(k+1)/ hp_b(k)) * ((T0(k+1) - T0(k)) + dTe(k+1))
          dS_k_t2 = (S0(k-1)-S0(k)) - &
                  (Kddt_h(k+1)/ hp_b(k)) * ((S0(k+1) - S0(k)) + dSe(k+1))
        endif
        dTe_term = dTe_t2 + hp_b(k) * (T0(k)-T0(k-1))
        dSe_term = dSe_t2 + hp_b(k) * (S0(k)-S0(k-1))
      else
        Th_a(k-1) = h_tr(k-1) * T0(k-1) ; Sh_a(k-1) = h_tr(k-1) * S0(k-1)
        if (K>=nz) then
          Th_b(k) = h_tr(k) * T0(k) ; Sh_b(k) = h_tr(k) * S0(k)
        else
          Th_b(k) = h_tr(k) * T0(k) + Kddt_h(K+1) * Te(k+1)
          Sh_b(k) = h_tr(k) * S0(k) + Kddt_h(k+1) * Se(k+1)
        endif
      endif

      ! Find the energy change due to a guess at the strength of diffusion at interface K.
      Kddt_h_guess = Kddt_h(K)

      if (old_PE_calc) then
        call find_PE_chg_orig(Kddt_h_guess, h_tr(k-1), hp_b(k), &
                         dTe_term, dSe_term, dT_k_t2, dS_k_t2, &
                         dT_to_dPE(k-1), dS_to_dPE(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                         pres(K), dT_to_dColHt(k-1), dS_to_dColHt(k-1), &
                         dT_to_dColHt_b(k), dS_to_dColHt_b(k), &
                         PE_chg=PE_chg_k(K,2), dPEc_dKd=dPEb_dKd(K))
      else
        call find_PE_chg(0.0, Kddt_h_guess, hp_a(k-1), hp_b(k), &
                         Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                         dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                         pres(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                         dT_to_dColHt_b(k), dS_to_dColHt_b(k), &
                         PE_chg=PE_chg_k(K,2), dPEc_dKd=dPEb_dKd(K), &
                         ColHt_cor=ColHt_cor_k(K,2))
      endif

      if (debug) then
        ! Compare with a 4th-order finite difference estimate.
        do itt=1,5
          Kddt_h_guess = (1.0+0.01*(itt-3))*Kddt_h(K)

          if (old_PE_calc) then
            call find_PE_chg_orig(Kddt_h_guess, h_tr(k-1), hp_b(k), &
                           dTe_term, dSe_term, dT_k_t2, dS_k_t2, &
                           dT_to_dPE(k-1), dS_to_dPE(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                           pres(K), dT_to_dColHt(k-1), dS_to_dColHt(k-1), &
                           dT_to_dColHt_b(k), dS_to_dColHt_b(k), &
                           PE_chg=PE_chg(itt))
          else
            call find_PE_chg(0.0, Kddt_h_guess, hp_a(k-1), hp_b(k), &
                             Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                             dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                             pres(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                             dT_to_dColHt_b(k), dS_to_dColHt_b(k), &
                             PE_chg=PE_chg(itt))
          endif
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

      b1 = 1.0 / (hp_b(k) + Kddt_h(K))
      c1_b(K) = Kddt_h(K) * b1
      if (k==nz) then
        Te(nz) = b1* (h_tr(nz)*T0(nz))
        Se(nz) = b1* (h_tr(nz)*S0(nz))
      else
        Te(k) = b1 * (h_tr(k) * T0(k) + Kddt_h(K+1) * Te(k+1))
        Se(k) = b1 * (h_tr(k) * S0(k) + Kddt_h(k+1) * Se(k+1))
      endif
      if (old_PE_calc) then
        dTe(k) = b1 * ( Kddt_h(K)*(T0(k-1)-T0(k)) + dTe_t2 )
        dSe(k) = b1 * ( Kddt_h(K)*(S0(k-1)-S0(k)) + dSe_t2 )
      endif

      hp_b(k-1) = h_tr(k-1) + (hp_b(k) * b1) * Kddt_h(K)
      dT_to_dPE_b(k-1) = dT_to_dPE(k-1) + c1_b(K)*dT_to_dPE_b(k)
      dS_to_dPE_b(k-1) = dS_to_dPE(k-1) + c1_b(K)*dS_to_dPE_b(k)
      dT_to_dColHt_b(k-1) = dT_to_dColHt(k-1) + c1_b(K)*dT_to_dColHt_b(k)
      dS_to_dColHt_b(k-1) = dS_to_dColHt(k-1) + c1_b(K)*dS_to_dColHt_b(k)

    enddo

    b1 = 1.0 / (hp_b(1))
    Tf(1) = b1 * (h_tr(1) * T0(1) + Kddt_h(2) * Te(2))
    Sf(1) = b1 * (h_tr(1) * S0(1) + Kddt_h(2) * Se(2))
    if (old_PE_calc) then
      dTe(1) = b1 * Kddt_h(2) * ((T0(2)-T0(1)) + dTe(2))
      dSe(1) = b1 * Kddt_h(2) * ((S0(2)-S0(1)) + dSe(2))
    endif

    do k=2,nz
      Tf(k) = Te(k) + c1_b(K)*Tf(k-1)
      Sf(k) = Se(k) + c1_b(K)*Sf(k-1)
    enddo

    if (debug) then
      do k=1,nz ; Tb(k) = Tf(k) ; Sb(k) = Sf(k) ; enddo
      PE_chg_tot1B = 0.0 ; PE_chg_tot2B = 0.0 ; T_chg_totB = 0.0
      do k=1,nz
        PE_chg_tot1B = PE_chg_tot1B + (dT_to_dPE(k) * (Tf(k) - T0(k)) + &
                                       dS_to_dPE(k) * (Sf(k) - S0(k)))
        T_chg_totB = T_chg_totB + h_tr(k) * (Tf(k) - T0(k))
      enddo
      do K=2,nz
        PE_chg_tot2B = PE_chg_tot2B + (PE_chg_k(K,2) - ColHt_cor_k(K,2))
      enddo
    endif

  endif

  if (central) then

    ! Set up values appropriate for no diffusivity.
    do k=1,nz
      hp_a(k) = h_tr(k) ; hp_b(k) = h_tr(k)
      dT_to_dPE_a(k) = dT_to_dPE(k) ; dS_to_dPE_a(k) = dS_to_dPE(k)
      dT_to_dPE_b(k) = dT_to_dPE(k) ; dS_to_dPE_b(k) = dS_to_dPE(k)
      dT_to_dColHt_a(k) = dT_to_dColHt(k) ; dS_to_dColHt_a(k) = dS_to_dColHt(k)
      dT_to_dColHt_b(k) = dT_to_dColHt(k) ; dS_to_dColHt_b(k) = dS_to_dColHt(k)
    enddo

    ! Calculate the dependencies on layers above.
    Kddt_h_a(1) = 0.0
    do K=2,nz  ! Loop over interior interfaces.
      ! First calculate some terms that are independent of the change in Kddt_h(K).
      Kd0 = 0.0  ! This might need to be changed - it is the already applied value of Kddt_h(K).
      if (K<=2) then
        Th_a(k-1) = h_tr(k-1) * T0(k-1) ; Sh_a(k-1) = h_tr(k-1) * S0(k-1)
      else
        Th_a(k-1) = h_tr(k-1) * T0(k-1) + Kddt_h(K-1) * Te_a(k-2)
        Sh_a(k-1) = h_tr(k-1) * S0(k-1) + Kddt_h(K-1) * Se_a(k-2)
      endif
      Th_b(k) = h_tr(k) * T0(k) ; Sh_b(k) = h_tr(k) * S0(k)

      Kddt_h_a(K) = 0.0
      if (K < K_cent) Kddt_h_a(K) = Kddt_h(K)

      dKd = Kddt_h_a(K)

      call find_PE_chg(Kd0, dKd, hp_a(k-1), hp_b(k), &
                       Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                       dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                       pres(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                       dT_to_dColHt_b(k), dS_to_dColHt_b(k), &
                       PE_chg=PE_change, ColHt_cor=ColHt_cor)
      PE_chg_k(K,3) = PE_change
      ColHt_cor_k(K,3) = ColHt_cor

      b1 = 1.0 / (hp_a(k-1) + Kddt_h_a(K))
      c1_a(K) = Kddt_h_a(K) * b1
      if (k==2) then
        Te_a(1) = b1*(h_tr(1)*T0(1))
        Se_a(1) = b1*(h_tr(1)*S0(1))
      else
        Te_a(k-1) = b1 * (h_tr(k-1) * T0(k-1) + Kddt_h_a(K-1) * Te_a(k-2))
        Se_a(k-1) = b1 * (h_tr(k-1) * S0(k-1) + Kddt_h_a(K-1) * Se_a(k-2))
      endif

      hp_a(k) = h_tr(k) + (hp_a(k-1) * b1) * Kddt_h_a(K)
      dT_to_dPE_a(k) = dT_to_dPE(k) + c1_a(K)*dT_to_dPE_a(k-1)
      dS_to_dPE_a(k) = dS_to_dPE(k) + c1_a(K)*dS_to_dPE_a(k-1)
      dT_to_dColHt_a(k) = dT_to_dColHt(k) + c1_a(K)*dT_to_dColHt_a(k-1)
      dS_to_dColHt_a(k) = dS_to_dColHt(k) + c1_a(K)*dS_to_dColHt_a(k-1)
    enddo

    ! Calculate the dependencies on layers below.
    Kddt_h_b(nz+1) = 0.0
    do K=nz,2,-1  ! Loop over interior interfaces.
      ! First calculate some terms that are independent of the change in Kddt_h(K).
      Kd0 = 0.0  ! This might need to be changed - it is the already applied value of Kddt_h(K).
!     if (K<=2) then
        Th_a(k-1) = h_tr(k-1) * T0(k-1) ; Sh_a(k-1) = h_tr(k-1) * S0(k-1)
!     else
!       Th_a(k-1) = h_tr(k-1) * T0(k-1) + Kddt_h(K-1) * Te_a(k-2)
!       Sh_a(k-1) = h_tr(k-1) * S0(k-1) + Kddt_h(K-1) * Se_a(k-2)
!     endif
      if (K>=nz) then
        Th_b(k) = h_tr(k) * T0(k) ; Sh_b(k) = h_tr(k) * S0(k)
      else
        Th_b(k) = h_tr(k) * T0(k) + Kddt_h(K+1) * Te_b(k+1)
        Sh_b(k) = h_tr(k) * S0(k) + Kddt_h(k+1) * Se_b(k+1)
      endif

      Kddt_h_b(K) = 0.0 ; if (K > K_cent) Kddt_h_b(K) = Kddt_h(K)
      dKd = Kddt_h_b(K)

      call find_PE_chg(Kd0, dKd, hp_a(k-1), hp_b(k), &
                       Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                       dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                       pres(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                       dT_to_dColHt_b(k), dS_to_dColHt_b(k), &
                       PE_chg=PE_change, ColHt_cor=ColHt_cor)
      PE_chg_k(K,3) = PE_chg_k(K,3) + PE_change
      ColHt_cor_k(K,3) = ColHt_cor_k(K,3) + ColHt_cor

      b1 = 1.0 / (hp_b(k) + Kddt_h_b(K))
      c1_b(K) = Kddt_h_b(K) * b1
      if (k==nz) then
        Te_b(k) = b1 * (h_tr(k)*T0(k))
        Se_b(k) = b1 * (h_tr(k)*S0(k))
      else
        Te_b(k) = b1 * (h_tr(k) * T0(k) + Kddt_h_b(K+1) * Te_b(k+1))
        Se_b(k) = b1 * (h_tr(k) * S0(k) + Kddt_h_b(k+1) * Se_b(k+1))
      endif

      hp_b(k-1) = h_tr(k-1) + (hp_b(k) * b1) * Kddt_h_b(K)
      dT_to_dPE_b(k-1) = dT_to_dPE(k-1) + c1_b(K)*dT_to_dPE_b(k)
      dS_to_dPE_b(k-1) = dS_to_dPE(k-1) + c1_b(K)*dS_to_dPE_b(k)
      dT_to_dColHt_b(k-1) = dT_to_dColHt(k-1) + c1_b(K)*dT_to_dColHt_b(k)
      dS_to_dColHt_b(k-1) = dS_to_dColHt(k-1) + c1_b(K)*dS_to_dColHt_b(k)

    enddo

    ! Calculate the final solution for the layers surrounding interface K_cent
    K = K_cent

    ! First calculate some terms that are independent of the change in Kddt_h(K).
    Kd0 = 0.0  ! This might need to be changed - it is the already applied value of Kddt_h(K).
    if (K<=2) then
      Th_a(k-1) = h_tr(k-1) * T0(k-1) ; Sh_a(k-1) = h_tr(k-1) * S0(k-1)
    else
      Th_a(k-1) = h_tr(k-1) * T0(k-1) + Kddt_h(K-1) * Te_a(k-2)
      Sh_a(k-1) = h_tr(k-1) * S0(k-1) + Kddt_h(K-1) * Se_a(k-2)
    endif
    if (K>=nz) then
      Th_b(k) = h_tr(k) * T0(k) ; Sh_b(k) = h_tr(k) * S0(k)
    else
      Th_b(k) = h_tr(k) * T0(k) + Kddt_h(K+1) * Te_b(k+1)
      Sh_b(k) = h_tr(k) * S0(k) + Kddt_h(k+1) * Se_b(k+1)
    endif

    dKd = Kddt_h(K)

    call find_PE_chg(Kd0, dKd, hp_a(k-1), hp_b(k), &
                     Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                     dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                     pres(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                     dT_to_dColHt_b(k), dS_to_dColHt_b(k), &
                     PE_chg=PE_change, ColHt_cor=ColHt_cor)
    PE_chg_k(K,3) = PE_chg_k(K,3) + PE_change
    ColHt_cor_k(K,3) = ColHt_cor_k(K,3) + ColHt_cor


    !   To derive the following, first to a partial update for the estimated
    ! temperatures and salinities in the layers around this interface:
    !   b1_a = 1.0 / (hp_a(k-1) + Kddt_h(K))
    !   b1_b = 1.0 / (hp_b(k) + Kddt_h(K))
    !   Te_up(k) = Th_b(k) * b1_b ; Se_up(k) = Sh_b(k) * b1_b
    !   Te_up(k-1) = Th_a(k-1) * b1_a ; Se_up(k-1) = Sh_a(k-1) * b1_a
    ! Find the final values of T & S for both layers around K_cent, using that
    !   c1_a(K) = Kddt_h(K) * b1_a ; c1_b(K) = Kddt_h(K) * b1_b
    !   Tf(K_cent-1) = Te_up(K_cent-1) + c1_a(K_cent)*Tf(K_cent)
    !   Tf(K_cent)   = Te_up(K_cent) + c1_b(K_cent)*Tf(K_cent-1)
    ! but further exploiting the expressions for c1_a and c1_b to avoid
    ! subtraction in the denominator, and use only a single division.
    b1 = 1.0 / (hp_a(k-1)*hp_b(k) + Kddt_h(K)*(hp_a(k-1) + hp_b(k)))
    Tf(k-1) = ((hp_b(k) + Kddt_h(K)) * Th_a(k-1) + Kddt_h(K) * Th_b(k) ) * b1
    Sf(k-1) = ((hp_b(k) + Kddt_h(K)) * Sh_a(k-1) + Kddt_h(K) * Sh_b(k) ) * b1
    Tf(k) = (Kddt_h(K) * Th_a(k-1) + (hp_a(k-1) + Kddt_h(K)) * Th_b(k) ) * b1
    Sf(k) = (Kddt_h(K) * Sh_a(k-1) + (hp_a(k-1) + Kddt_h(K)) * Sh_b(k) ) * b1

    c1_a(K) = Kddt_h(K) / (hp_a(k-1) + Kddt_h(K))
    c1_b(K) = Kddt_h(K) / (hp_b(k) + Kddt_h(K))

    ! Now update the other layer working outward from k_cent to determine the final
    ! temperatures and salinities.
    do k=K_cent-2,1,-1
      Tf(k) = Te_a(k) + c1_a(K+1)*Tf(k+1)
      Sf(k) = Se_a(k) + c1_a(K+1)*Sf(k+1)
    enddo
    do k=K_cent+1,nz
      Tf(k) = Te_b(k) + c1_b(K)*Tf(k-1)
      Sf(k) = Se_b(k) + c1_b(K)*Sf(k-1)
    enddo

    if (debug) then
      PE_chg_tot1C = 0.0 ; PE_chg_tot2C = 0.0 ; T_chg_totC = 0.0
      do k=1,nz
        PE_chg_tot1C = PE_chg_tot1C + (dT_to_dPE(k) * (Tf(k) - T0(k)) + &
                                       dS_to_dPE(k) * (Sf(k) - S0(k)))
        T_chg_totC = T_chg_totC + h_tr(k) * (Tf(k) - T0(k))
      enddo
      do K=2,nz
        PE_chg_tot2C = PE_chg_tot2C + (PE_chg_k(K,3) - ColHt_cor_k(K,3))
      enddo
    endif

  endif

  if (halves) then

    ! Set up values appropriate for no diffusivity.
    do k=1,nz
      hp_a(k) = h_tr(k) ; hp_b(k) = h_tr(k)
      dT_to_dPE_a(k) = dT_to_dPE(k) ; dS_to_dPE_a(k) = dS_to_dPE(k)
      dT_to_dPE_b(k) = dT_to_dPE(k) ; dS_to_dPE_b(k) = dS_to_dPE(k)
      dT_to_dColHt_a(k) = dT_to_dColHt(k) ; dS_to_dColHt_a(k) = dS_to_dColHt(k)
      dT_to_dColHt_b(k) = dT_to_dColHt(k) ; dS_to_dColHt_b(k) = dS_to_dColHt(k)
    enddo
    do K=1,nz+1
      Kd_so_far(K) = 0.0
    enddo

    ! Calculate the dependencies on layers above.
    Kddt_h_a(1) = 0.0
    do K=2,nz  ! Loop over interior interfaces.
      ! First calculate some terms that are independent of the change in Kddt_h(K).
      Kd0 = Kd_so_far(K)
      if (K<=2) then
        Th_a(k-1) = h_tr(k-1) * T0(k-1) ; Sh_a(k-1) = h_tr(k-1) * S0(k-1)
      else
        Th_a(k-1) = h_tr(k-1) * T0(k-1) + Kd_so_far(K-1) * Te(k-2)
        Sh_a(k-1) = h_tr(k-1) * S0(k-1) + Kd_so_far(K-1) * Se(k-2)
      endif
      Th_b(k) = h_tr(k) * T0(k) ; Sh_b(k) = h_tr(k) * S0(k)

      dKd = 0.5 * Kddt_h(K) - Kd_so_far(K)

      call find_PE_chg(Kd0, dKd, hp_a(k-1), hp_b(k), &
                       Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                       dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                       pres(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                       dT_to_dColHt_b(k), dS_to_dColHt_b(k), &
                       PE_chg=PE_change, ColHt_cor=ColHt_cor)

      PE_chg_k(K,4) = PE_change
      ColHt_cor_k(K,4) = ColHt_cor

      Kd_so_far(K) = Kd_so_far(K) + dKd

      b1 = 1.0 / (hp_a(k-1) + Kd_so_far(K))
      c1_a(K) = Kd_so_far(K) * b1
      if (k==2) then
        Te(1) = b1*(h_tr(1)*T0(1))
        Se(1) = b1*(h_tr(1)*S0(1))
      else
        Te(k-1) = b1 * (h_tr(k-1) * T0(k-1) + Kd_so_far(K-1) * Te(k-2))
        Se(k-1) = b1 * (h_tr(k-1) * S0(k-1) + Kd_so_far(K-1) * Se(k-2))
      endif

      hp_a(k) = h_tr(k) + (hp_a(k-1) * b1) * Kd_so_far(K)
      dT_to_dPE_a(k) = dT_to_dPE(k) + c1_a(K)*dT_to_dPE_a(k-1)
      dS_to_dPE_a(k) = dS_to_dPE(k) + c1_a(K)*dS_to_dPE_a(k-1)
      dT_to_dColHt_a(k) = dT_to_dColHt(k) + c1_a(K)*dT_to_dColHt_a(k-1)
      dS_to_dColHt_a(k) = dS_to_dColHt(k) + c1_a(K)*dS_to_dColHt_a(k-1)
    enddo

    ! Calculate the dependencies on layers below.
    do K=nz,2,-1  ! Loop over interior interfaces.
      ! First calculate some terms that are independent of the change in Kddt_h(K).
      Kd0 = Kd_so_far(K)
      if (K<=2) then
        Th_a(k-1) = h_tr(k-1) * T0(k-1) ; Sh_a(k-1) = h_tr(k-1) * S0(k-1)
      else
        Th_a(k-1) = h_tr(k-1) * T0(k-1) + Kd_so_far(K-1) * Te(k-2)
        Sh_a(k-1) = h_tr(k-1) * S0(k-1) + Kd_so_far(K-1) * Se(k-2)
      endif
      if (K>=nz) then
        Th_b(k) = h_tr(k) * T0(k) ; Sh_b(k) = h_tr(k) * S0(k)
      else
        Th_b(k) = h_tr(k) * T0(k) + Kd_so_far(K+1) * Te(k+1)
        Sh_b(k) = h_tr(k) * S0(k) + Kd_so_far(k+1) * Se(k+1)
      endif

      dKd = Kddt_h(K) - Kd_so_far(K)

      call find_PE_chg(Kd0, dKd, hp_a(k-1), hp_b(k), &
                       Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                       dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                       pres(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                       dT_to_dColHt_b(k), dS_to_dColHt_b(k), &
                       PE_chg=PE_change, ColHt_cor=ColHt_cor)

      PE_chg_k(K,4) = PE_chg_k(K,4) + PE_change
      ColHt_cor_k(K,4) = ColHt_cor_k(K,4) + ColHt_cor


      Kd_so_far(K) = Kd_so_far(K) + dKd

      b1 = 1.0 / (hp_b(k) + Kd_so_far(K))
      c1_b(K) = Kd_so_far(K) * b1
      if (k==nz) then
        Te(k) = b1 * (h_tr(k)*T0(k))
        Se(k) = b1 * (h_tr(k)*S0(k))
      else
        Te(k) = b1 * (h_tr(k) * T0(k) + Kd_so_far(K+1) * Te(k+1))
        Se(k) = b1 * (h_tr(k) * S0(k) + Kd_so_far(k+1) * Se(k+1))
      endif

      hp_b(k-1) = h_tr(k-1) + (hp_b(k) * b1) * Kd_so_far(K)
      dT_to_dPE_b(k-1) = dT_to_dPE(k-1) + c1_b(K)*dT_to_dPE_b(k)
      dS_to_dPE_b(k-1) = dS_to_dPE(k-1) + c1_b(K)*dS_to_dPE_b(k)
      dT_to_dColHt_b(k-1) = dT_to_dColHt(k-1) + c1_b(K)*dT_to_dColHt_b(k)
      dS_to_dColHt_b(k-1) = dS_to_dColHt(k-1) + c1_b(K)*dS_to_dColHt_b(k)

    enddo

    ! Now update the other layer working down from the top to determine the
    ! final temperatures and salinities.
    b1 = 1.0 / (hp_b(1))
    Tf(1) = b1 * (h_tr(1) * T0(1) + Kddt_h(2) * Te(2))
    Sf(1) = b1 * (h_tr(1) * S0(1) + Kddt_h(2) * Se(2))
    do k=2,nz
      Tf(k) = Te(k) + c1_b(K)*Tf(k-1)
      Sf(k) = Se(k) + c1_b(K)*Sf(k-1)
    enddo

    if (debug) then
      PE_chg_tot1D = 0.0 ; PE_chg_tot2D = 0.0 ; T_chg_totD = 0.0
      do k=1,nz
        PE_chg_tot1D = PE_chg_tot1D + (dT_to_dPE(k) * (Tf(k) - T0(k)) + &
                                       dS_to_dPE(k) * (Sf(k) - S0(k)))
        T_chg_totD = T_chg_totD + h_tr(k) * (Tf(k) - T0(k))
      enddo
      do K=2,nz
        PE_chg_tot2D = PE_chg_tot2D + (PE_chg_k(K,4) - ColHt_cor_k(K,4))
      enddo
    endif

  endif

  energy_Kd = 0.0 ; do K=2,nz ; energy_Kd = energy_Kd + PE_chg_k(K,1) ; enddo
  energy_Kd = energy_Kd / dt
  K=nz

  if (do_print) then
    if (CS%id_ERt>0) call post_data_1d_k(CS%id_ERt, PE_chg_k(:,1), CS%diag)
    if (CS%id_ERb>0) call post_data_1d_k(CS%id_ERb, PE_chg_k(:,2), CS%diag)
    if (CS%id_ERc>0) call post_data_1d_k(CS%id_ERc, PE_chg_k(:,3), CS%diag)
    if (CS%id_ERh>0) call post_data_1d_k(CS%id_ERh, PE_chg_k(:,4), CS%diag)
    if (CS%id_Kddt>0) call post_data_1d_k(CS%id_Kddt, GV%H_to_m*Kddt_h, CS%diag)
    if (CS%id_Kd>0)   call post_data_1d_k(CS%id_Kd, Kd, CS%diag)
    if (CS%id_h>0)    call post_data_1d_k(CS%id_h, GV%H_to_m*h_tr, CS%diag)
    if (CS%id_zInt>0) call post_data_1d_k(CS%id_zInt, Z_int, CS%diag)
    if (CS%id_CHCt>0) call post_data_1d_k(CS%id_CHCt, ColHt_cor_k(:,1), CS%diag)
    if (CS%id_CHCb>0) call post_data_1d_k(CS%id_CHCb, ColHt_cor_k(:,2), CS%diag)
    if (CS%id_CHCc>0) call post_data_1d_k(CS%id_CHCc, ColHt_cor_k(:,3), CS%diag)
    if (CS%id_CHCh>0) call post_data_1d_k(CS%id_CHCh, ColHt_cor_k(:,4), CS%diag)
    if (CS%id_T0>0) call post_data_1d_k(CS%id_T0, T0, CS%diag)
    if (CS%id_Tf>0) call post_data_1d_k(CS%id_Tf, Tf, CS%diag)
    if (CS%id_S0>0) call post_data_1d_k(CS%id_S0, S0, CS%diag)
    if (CS%id_Sf>0) call post_data_1d_k(CS%id_Sf, Sf, CS%diag)
    if (CS%id_N2_0>0) then
      N2(1) = 0.0 ; N2(nz+1) = 0.0
      do K=2,nz
        call calculate_density(0.5*(T0(k-1) + T0(k)), 0.5*(S0(k-1) + S0(k)), &
                               pres(K), rho_here, tv%eqn_of_state)
        N2(K) = (GV%g_Earth * rho_here / (0.5*GV%H_to_m*(h_tr(k-1) + h_tr(k)))) * &
                ( 0.5*(dSV_dT(k-1) + dSV_dT(k)) * (T0(k-1) - T0(k)) + &
                  0.5*(dSV_dS(k-1) + dSV_dS(k)) * (S0(k-1) - S0(k)) )
      enddo
      call post_data_1d_k(CS%id_N2_0, N2, CS%diag)
    endif
    if (CS%id_N2_f>0) then
      N2(1) = 0.0 ; N2(nz+1) = 0.0
      do K=2,nz
        call calculate_density(0.5*(Tf(k-1) + Tf(k)), 0.5*(Sf(k-1) + Sf(k)), &
                               pres(K), rho_here, tv%eqn_of_state)
        N2(K) = (GV%g_Earth * rho_here / (0.5*GV%H_to_m*(h_tr(k-1) + h_tr(k)))) * &
                ( 0.5*(dSV_dT(k-1) + dSV_dT(k)) * (Tf(k-1) - Tf(k)) + &
                  0.5*(dSV_dS(k-1) + dSV_dS(k)) * (Sf(k-1) - Sf(k)) )
      enddo
      call post_data_1d_k(CS%id_N2_f, N2, CS%diag)
    endif
  endif

end subroutine diapyc_energy_req_calc

!> This subroutine calculates the change in potential energy and or derivatives
!! for several changes in an interfaces's diapycnal diffusivity times a timestep.
subroutine find_PE_chg(Kddt_h0, dKddt_h, hp_a, hp_b, Th_a, Sh_a, Th_b, Sh_b, &
                       dT_to_dPE_a, dS_to_dPE_a, dT_to_dPE_b, dS_to_dPE_b, &
                       pres, dT_to_dColHt_a, dS_to_dColHt_a, dT_to_dColHt_b, dS_to_dColHt_b, &
                       PE_chg, dPEc_dKd, dPE_max, dPEc_dKd_0, ColHt_cor)
  real, intent(in)  :: Kddt_h0  !< The previously used diffusivity at an interface times
                                !! the time step and  divided by the average of the
                                !! thicknesses around the interface, in units of H (m or kg-2).
  real, intent(in)  :: dKddt_h  !< The trial change in the diffusivity at an interface times
                                !! the time step and  divided by the average of the
                                !! thicknesses around the interface, in units of H (m or kg-2).
  real, intent(in)  :: hp_a     !< The effective pivot thickness of the layer above the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above, in H.
  real, intent(in)  :: hp_b     !< The effective pivot thickness of the layer below the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above, in H.
  real, intent(in)  :: Th_a     !< An effective temperature times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers, in K H.
  real, intent(in)  :: Sh_a     !< An effective salinity times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers, in K H.
  real, intent(in)  :: Th_b     !< An effective temperature times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers, in K H.
  real, intent(in)  :: Sh_b     !< An effective salinity times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers, in K H.
  real, intent(in)  :: dT_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! potential energy, including all implicit diffusive changes
                                !! in the temperatures of all the layers above, in J m-2 K-1.
  real, intent(in)  :: dS_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! potential energy, including all implicit diffusive changes
                                !! in the salinities of all the layers above, in J m-2 ppt-1.
  real, intent(in)  :: dT_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! potential energy, including all implicit diffusive changes
                                !! in the temperatures of all the layers below, in J m-2 K-1.
  real, intent(in)  :: dS_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! potential energy, including all implicit diffusive changes
                                !! in the salinities of all the layers below, in J m-2 ppt-1.
  real, intent(in)  :: pres     !< The hydrostatic interface pressure, which is used to relate
                                !! the changes in column thickness to the energy that is radiated
                                !! as gravity waves and unavailable to drive mixing, in Pa.
  real, intent(in)  :: dT_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers above, in m K-1.
  real, intent(in)  :: dS_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers above, in m ppt-1.
  real, intent(in)  :: dT_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers below, in m K-1.
  real, intent(in)  :: dS_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers below, in m ppt-1.

  real, optional, intent(out) :: PE_chg   !< The change in column potential energy from applying
                                          !! Kddt_h at the present interface, in J m-2.
  real, optional, intent(out) :: dPEc_dKd !< The partial derivative of PE_chg with Kddt_h,
                                          !! in units of J m-2 H-1.
  real, optional, intent(out) :: dPE_max  !< The maximum change in column potential energy that could
                                          !! be realizedd by applying a huge value of Kddt_h at the
                                          !! present interface, in J m-2.
  real, optional, intent(out) :: dPEc_dKd_0 !< The partial derivative of PE_chg with Kddt_h in the
                                            !! limit where Kddt_h = 0, in J m-2 H-1.
  real, optional, intent(out) :: ColHt_cor  !< The correction to PE_chg that is made due to a net
                                            !! change in the column height, in J m-2.

  real :: hps ! The sum of the two effective pivot thicknesses, in H.
  real :: bdt1 ! A product of the two pivot thicknesses plus a diffusive term, in H2.
  real :: dT_c ! The core term in the expressions for the temperature changes, in K H2.
  real :: dS_c ! The core term in the expressions for the salinity changes, in psu H2.
  real :: PEc_core ! The diffusivity-independent core term in the expressions
                   ! for the potential energy changes, J m-3.
  real :: ColHt_core ! The diffusivity-independent core term in the expressions
                     ! for the column height changes, J m-3.
  real :: ColHt_chg  ! The change in the column height, in m.
  real :: y1   ! A local temporary term, in units of H-3 or H-4 in various contexts.

  !   The expression for the change in potential energy used here is derived
  ! from the expression for the final estimates of the changes in temperature
  ! and salinities, and then extensively manipulated to get it into its most
  ! succint form. The derivation is not necessarily obvious, but it demonstrably
  ! works by comparison with separate calculations of the energy changes after
  ! the tridiagonal solver for the final changes in temperature and salinity are
  ! applied.

  hps = hp_a + hp_b
  bdt1 = hp_a * hp_b + Kddt_h0 * hps
  dT_c = hp_a * Th_b - hp_b * Th_a
  dS_c = hp_a * Sh_b - hp_b * Sh_a
  PEc_core = hp_b * (dT_to_dPE_a * dT_c + dS_to_dPE_a * dS_c) - &
             hp_a * (dT_to_dPE_b * dT_c + dS_to_dPE_b * dS_c)
  ColHt_core = hp_b * (dT_to_dColHt_a * dT_c + dS_to_dColHt_a * dS_c) - &
               hp_a * (dT_to_dColHt_b * dT_c + dS_to_dColHt_b * dS_c)

  if (present(PE_chg)) then
    ! Find the change in column potential energy due to the change in the
    ! diffusivity at this interface by dKddt_h.
    y1 = dKddt_h / (bdt1 * (bdt1 + dKddt_h * hps))
    PE_chg = PEc_core * y1
    ColHt_chg = ColHt_core * y1
    if (ColHt_chg < 0.0) PE_chg = PE_chg - pres * ColHt_chg
    if (present(ColHt_cor)) ColHt_cor = -pres * min(ColHt_chg, 0.0)
  else if (present(ColHt_cor)) then
    y1 = dKddt_h / (bdt1 * (bdt1 + dKddt_h * hps))
    ColHt_cor = -pres * min(ColHt_core * y1, 0.0)
  endif

  if (present(dPEc_dKd)) then
    ! Find the derivative of the potential energy change with dKddt_h.
    y1 = 1.0 / (bdt1 + dKddt_h * hps)**2
    dPEc_dKd = PEc_core * y1
    ColHt_chg = ColHt_core * y1
    if (ColHt_chg < 0.0) dPEc_dKd = dPEc_dKd - pres * ColHt_chg
  endif

  if (present(dPE_max)) then
    ! This expression is the limit of PE_chg for infinite dKddt_h.
    y1 = 1.0 / (bdt1 * hps)
    dPE_max = PEc_core * y1
    ColHt_chg = ColHt_core * y1
    if (ColHt_chg < 0.0) dPE_max = dPE_max - pres * ColHt_chg
  endif

  if (present(dPEc_dKd_0)) then
    ! This expression is the limit of dPEc_dKd for dKddt_h = 0.
    y1 = 1.0 / bdt1**2
    dPEc_dKd_0 = PEc_core * y1
    ColHt_chg = ColHt_core * y1
    if (ColHt_chg < 0.0) dPEc_dKd_0 = dPEc_dKd_0 - pres * ColHt_chg
  endif

end subroutine find_PE_chg


!> This subroutine calculates the change in potential energy and or derivatives
!! for several changes in an interfaces's diapycnal diffusivity times a timestep
!! using the original form used in the first version of ePBL.
subroutine find_PE_chg_orig(Kddt_h, h_k, b_den_1, dTe_term, dSe_term, &
                            dT_km1_t2, dS_km1_t2, dT_to_dPE_k, dS_to_dPE_k, &
                            dT_to_dPEa, dS_to_dPEa, pres, dT_to_dColHt_k, &
                            dS_to_dColHt_k, dT_to_dColHta, dS_to_dColHta, &
                            PE_chg, dPEc_dKd, dPE_max, dPEc_dKd_0)
  real, intent(in)  :: Kddt_h   !< The diffusivity at an interface times the time step and
                                !! divided by the average of the thicknesses around the
                                !! interface, in units of H (m or kg-2).
  real, intent(in)  :: h_k      !< The thickness of the layer below the interface, in H.
  real, intent(in)  :: b_den_1  !< The first term in the denominator of the pivot
                                !! for the tridiagonal solver, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above, in H.
  real, intent(in)  :: dTe_term !< A diffusivity-independent term related to the
                                !! temperature change in the layer below the interface, in K H.
  real, intent(in)  :: dSe_term !< A diffusivity-independent term related to the
                                !! salinity change in the layer below the interface, in ppt H.
  real, intent(in)  :: dT_km1_t2 !< A diffusivity-independent term related to the
                                 !! temperature change in the layer above the interface, in K.
  real, intent(in)  :: dS_km1_t2 !< A diffusivity-independent term related to the
                                 !! salinity change in the layer above the interface, in ppt.
  real, intent(in)  :: pres      !< The hydrostatic interface pressure, which is used to relate
                                 !! the changes in column thickness to the energy that is radiated
                                 !! as gravity waves and unavailable to drive mixing, in Pa.
  real, intent(in)  :: dT_to_dPE_k !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! potential energy, including all implicit diffusive changes
                                 !! in the temperatures of all the layers below, in J m-2 K-1.
  real, intent(in)  :: dS_to_dPE_k !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! potential energy, including all implicit diffusive changes
                                 !! in the salinities of all the layers below, in J m-2 ppt-1.
  real, intent(in)  :: dT_to_dPEa !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! potential energy, including all implicit diffusive changes
                                 !! in the temperatures of all the layers above, in J m-2 K-1.
  real, intent(in)  :: dS_to_dPEa !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! potential energy, including all implicit diffusive changes
                                 !! in the salinities of all the layers above, in J m-2 ppt-1.
  real, intent(in)  :: dT_to_dColHt_k !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the temperatures of all the layers below, in m K-1.
  real, intent(in)  :: dS_to_dColHt_k !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the salinities of all the layers below, in m ppt-1.
  real, intent(in)  :: dT_to_dColHta !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the temperatures of all the layers above, in m K-1.
  real, intent(in)  :: dS_to_dColHta !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the salinities of all the layers above, in m ppt-1.

  real, optional, intent(out) :: PE_chg   !< The change in column potential energy from applying
                                          !! Kddt_h at the present interface, in J m-2.
  real, optional, intent(out) :: dPEc_dKd !< The partial derivative of PE_chg with Kddt_h,
                                          !! in units of J m-2 H-1.
  real, optional, intent(out) :: dPE_max  !< The maximum change in column potential energy that could
                                          !! be realizedd by applying a huge value of Kddt_h at the
                                          !! present interface, in J m-2.
  real, optional, intent(out) :: dPEc_dKd_0 !< The partial derivative of PE_chg with Kddt_h in the
                                            !! limit where Kddt_h = 0, in J m-2 H-1.

!   This subroutine determines the total potential energy change due to mixing
! at an interface, including all of the implicit effects of the prescribed
! mixing at interfaces above.  Everything here is derived by careful manipulation
! of the robust tridiagonal solvers used for tracers by MOM6.  The results are
! positive for mixing in a stably stratified environment.
!   The comments describing these arguments are for a downward mixing pass, but
! this routine can also be used for an upward pass with the sense of direction
! reversed.

  real :: b1            ! b1 is used by the tridiagonal solver, in H-1.
  real :: b1Kd          ! Temporary array (nondim.)
  real :: ColHt_chg     ! The change in column thickness in m.
  real :: dColHt_max    ! The change in column thickess for infinite diffusivity, in m.
  real :: dColHt_dKd    ! The partial derivative of column thickess with diffusivity, in s m-1.
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
    ColHt_chg = (dT_to_dColHt_k * dT_k + dT_to_dColHta * dT_km1) + &
                (dS_to_dColHt_k * dS_k + dS_to_dColHta * dS_km1)
    if (ColHt_chg < 0.0) PE_chg = PE_chg - pres * ColHt_chg
  endif

  if (present(dPEc_dKd)) then
    ! Find the derivatives of the temperature and salinity changes with Kddt_h.
    dKr_dKd = (h_k*b_den_1) * I_Kr_denom**2

    ddT_k_dKd = dKr_dKd * dTe_term
    ddS_k_dKd = dKr_dKd * dSe_term
    ddT_km1_dKd = (b1**2 * b_den_1) * ( dT_k + dT_km1_t2 ) + b1Kd * ddT_k_dKd
    ddS_km1_dKd = (b1**2 * b_den_1) * ( dS_k + dS_km1_t2 ) + b1Kd * ddS_k_dKd

    ! Calculate the partial derivative of Pe_chg with Kddt_h.
    dPEc_dKd = (dT_to_dPE_k * ddT_k_dKd + dT_to_dPEa * ddT_km1_dKd) + &
               (dS_to_dPE_k * ddS_k_dKd + dS_to_dPEa * ddS_km1_dKd)
    dColHt_dKd = (dT_to_dColHt_k * ddT_k_dKd + dT_to_dColHta * ddT_km1_dKd) + &
                 (dS_to_dColHt_k * ddS_k_dKd + dS_to_dColHta * ddS_km1_dKd)
    if (dColHt_dKd < 0.0) dPEc_dKd = dPEc_dKd - pres * dColHt_dKd
  endif

  if (present(dPE_max)) then
    ! This expression is the limit of PE_chg for infinite Kddt_h.
    dPE_max = (dT_to_dPEa * dT_km1_t2 + dS_to_dPEa * dS_km1_t2) + &
              ((dT_to_dPE_k + dT_to_dPEa) * dTe_term + &
               (dS_to_dPE_k + dS_to_dPEa) * dSe_term) / (b_den_1 + h_k)
    dColHt_max = (dT_to_dColHta * dT_km1_t2 + dS_to_dColHta * dS_km1_t2) + &
              ((dT_to_dColHt_k + dT_to_dColHta) * dTe_term + &
               (dS_to_dColHt_k + dS_to_dColHta) * dSe_term) / (b_den_1 + h_k)
    if (dColHt_max < 0.0) dPE_max = dPE_max - pres*dColHt_max
  endif

  if (present(dPEc_dKd_0)) then
    ! This expression is the limit of dPEc_dKd for Kddt_h = 0.
    dPEc_dKd_0 = (dT_to_dPEa * dT_km1_t2 + dS_to_dPEa * dS_km1_t2) / (b_den_1) + &
                 (dT_to_dPE_k * dTe_term + dS_to_dPE_k * dSe_term) / (h_k*b_den_1)
    dColHt_dKd = (dT_to_dColHta * dT_km1_t2 + dS_to_dColHta * dS_km1_t2) / (b_den_1) + &
                 (dT_to_dColHt_k * dTe_term + dS_to_dColHt_k * dSe_term) / (h_k*b_den_1)
    if (dColHt_dKd < 0.0) dPEc_dKd_0 = dPEc_dKd_0 - pres*dColHt_dKd
  endif

end subroutine find_PE_chg_orig

subroutine diapyc_energy_req_init(Time, G, param_file, diag, CS)
  type(time_type),            intent(in)    :: Time        !< model time
  type(ocean_grid_type),      intent(in)    :: G           !< model grid structure
  type(param_file_type),      intent(in)    :: param_file  !< file to parse for parameter values
  type(diag_ctrl),    target, intent(inout) :: diag        !< structure to regulate diagnostic output
  type(diapyc_energy_req_CS), pointer       :: CS          !< module control structure
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
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "ENERGY_REQ_KH_SCALING", CS%test_Kh_scaling, &
                 "A scaling factor for the diapycnal diffusivity used in \n"//&
                 "testing the energy requirements.", default=1.0, units="nondim")
  call get_param(param_file, mod, "ENERGY_REQ_COL_HT_SCALING", CS%ColHt_scaling, &
                 "A scaling factor for the column height change correction \n"//&
                 "used in testing the energy requirements.", default=1.0, units="nondim")
  call get_param(param_file, mod, "ENERGY_REQ_USE_TEST_PROFILE", &
                 CS%use_test_Kh_profile, &
                 "If true, use the internal test diffusivity profile in \n"//&
                 "place of any that might be passed in as an argument.", default=.false.)

  CS%id_ERt = register_diag_field('ocean_model', 'EnReqTest_ERt', diag%axesZi, Time, &
                 "Diffusivity Energy Requirements, top-down", "J m-2")
  CS%id_ERb = register_diag_field('ocean_model', 'EnReqTest_ERb', diag%axesZi, Time, &
                 "Diffusivity Energy Requirements, bottom-up", "J m-2")
  CS%id_ERc = register_diag_field('ocean_model', 'EnReqTest_ERc', diag%axesZi, Time, &
                 "Diffusivity Energy Requirements, center-last", "J m-2")
  CS%id_ERh = register_diag_field('ocean_model', 'EnReqTest_ERh', diag%axesZi, Time, &
                 "Diffusivity Energy Requirements, halves", "J m-2")
  CS%id_Kddt = register_diag_field('ocean_model', 'EnReqTest_Kddt', diag%axesZi, Time, &
                 "Implicit diffusive coupling coefficient", "m")
  CS%id_Kd = register_diag_field('ocean_model', 'EnReqTest_Kd', diag%axesZi, Time, &
                 "Diffusivity in test", "m2 s-1")
  CS%id_h   = register_diag_field('ocean_model', 'EnReqTest_h_lay', diag%axesZL, Time, &
                 "Test column layer thicknesses", "m")
  CS%id_zInt   = register_diag_field('ocean_model', 'EnReqTest_z_int', diag%axesZi, Time, &
                 "Test column layer interface heights", "m")
  CS%id_CHCt = register_diag_field('ocean_model', 'EnReqTest_CHCt', diag%axesZi, Time, &
                 "Column Height Correction to Energy Requirements, top-down", "J m-2")
  CS%id_CHCb = register_diag_field('ocean_model', 'EnReqTest_CHCb', diag%axesZi, Time, &
                 "Column Height Correction to Energy Requirements, bottom-up", "J m-2")
  CS%id_CHCc = register_diag_field('ocean_model', 'EnReqTest_CHCc', diag%axesZi, Time, &
                 "Column Height Correction to Energy Requirements, center-last", "J m-2")
  CS%id_CHCh = register_diag_field('ocean_model', 'EnReqTest_CHCh', diag%axesZi, Time, &
                 "Column Height Correction to Energy Requirements, halves", "J m-2")
  CS%id_T0 = register_diag_field('ocean_model', 'EnReqTest_T0', diag%axesZL, Time, &
                 "Temperature before mixing", "deg C")
  CS%id_Tf = register_diag_field('ocean_model', 'EnReqTest_Tf', diag%axesZL, Time, &
                 "Temperature after mixing", "deg C")
  CS%id_S0 = register_diag_field('ocean_model', 'EnReqTest_S0', diag%axesZL, Time, &
                 "Salinity before mixing", "g kg-1")
  CS%id_Sf = register_diag_field('ocean_model', 'EnReqTest_Sf', diag%axesZL, Time, &
                 "Salinity after mixing", "g kg-1")
  CS%id_N2_0 = register_diag_field('ocean_model', 'EnReqTest_N2_0', diag%axesZi, Time, &
                 "Squared buoyancy frequency before mixing", "second-2")
  CS%id_N2_f = register_diag_field('ocean_model', 'EnReqTest_N2_f', diag%axesZi, Time, &
                 "Squared buoyancy frequency after mixing", "second-2")

end subroutine diapyc_energy_req_init

subroutine diapyc_energy_req_end(CS)
  type(diapyc_energy_req_CS), pointer :: CS
  if (associated(CS)) deallocate(CS)
end subroutine diapyc_energy_req_end

end module MOM_diapyc_energy_req
