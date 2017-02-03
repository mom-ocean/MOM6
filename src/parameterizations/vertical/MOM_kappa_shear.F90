module MOM_kappa_shear
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
!*  By Laura Jackson and Robert Hallberg, 2006-2008                    *
!*                                                                     *
!*    This file contains the subroutines that determine the diapycnal  *
!*  diffusivity driven by resolved shears, as specified by the         *
!*  parameterizations described in Jackson and Hallberg (JPO, 2008).   *
!*                                                                     *
!*    The technique by which the 6 equations (for kappa, TKE, u, v, T, *
!*  and S) are solved simultaneously has been dramatically revised     *
!*  from the previous version. The previous version was not converging *
!*  in some cases, especially near the surface mixed layer, while the  *
!*  revised version does.  The revised version solves for kappa and    *
!*  TKE with shear and stratification fixed, then marches the density  *
!*  and velocities forward with an adaptive (and aggressive) time step *
!*  in a predictor-corrector-corrector emulation of a trapezoidal      *
!*  scheme.  Run-time-settable parameters determine the tolerence to   *
!*  which the kappa and TKE equations are solved and the minimum time  *
!*  step that can be taken.                                            *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_debugging, only : hchksum
use MOM_error_handler, only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>
#ifdef use_netCDF
#include <netcdf.inc>
#endif

public Calculate_kappa_shear, kappa_shear_init, kappa_shear_is_used

type, public :: Kappa_shear_CS ! ; private
  real    :: RiNo_crit       ! The critical shear Richardson number for
                             ! shear-entrainment. The theoretical value is 0.25.
                             ! The values found by Jackson et al. are 0.25-0.35.
  real    :: Shearmix_rate   ! A nondimensional rate scale for shear-driven
                             ! entrainment.  The value given by Jackson et al.
                             ! is 0.085-0.089.
  real    :: FRi_curvature   !   A constant giving the curvature of the function
                             ! of the Richardson number that relates shear to
                             ! sources in the kappa equation, Nondim.
                             ! The values found by Jackson et al. are -0.97 - -0.89.
  real    :: C_N             !   The coefficient for the decay of TKE due to
                             ! stratification (i.e. proportional to N*tke), ND.
                             ! The values found by Jackson et al. are 0.24-0.28.
  real    :: C_S             !   The coefficient for the decay of TKE due to
                             ! shear (i.e. proportional to |S|*tke), ND.
                             ! The values found by Jackson et al. are 0.14-0.12.
  real    :: lambda          !   The coefficient for the buoyancy length scale
                             ! in the kappa equation.  Nondimensional.
                             ! The values found by Jackson et al. are 0.82-0.81.
  real    :: lambda2_N_S     !   The square of the ratio of the coefficients of
                             ! the buoyancy and shear scales in the diffusivity
                             ! equation, 0 to eliminate the shear scale. Nondim.
  real    :: TKE_bg          !   The background level of TKE, in m2 s-2.
  real    :: kappa_0         !   The background diapycnal diffusivity, in m2 s-1.
  real    :: kappa_tol_err   !   The fractional error in kappa that is tolerated.
  real    :: Prandtl_turb    ! Prandtl number used to convert Kd_turb into viscosity.
  integer :: nkml            !   The number of layers in the mixed layer, as
                             ! treated in this routine.  If the pieces of the
                             ! mixed layer are not to be treated collectively,
                             ! nkml is set to 1.
  integer :: max_RiNo_it     ! The maximum number of iterations that may be used
                             ! to estimate the instantaneous shear-driven mixing.
  integer :: max_KS_it       ! The maximum number of iterations that may be used
                             ! to estimate the time-averaged diffusivity.
  logical :: eliminate_massless ! If true, massless layers are merged with neighboring
                             ! massive layers in this calculation.  I can think of
                             ! no good reason why this should be false.
  logical :: layer_stagger = .false.
  logical :: debug = .false.
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                            ! timing of diagnostic output.
  integer :: id_Kd_shear = -1, id_TKE = -1
  integer :: id_ILd2 = -1, id_dz_Int = -1
end type Kappa_shear_CS

! integer :: id_clock_project, id_clock_KQ, id_clock_avg, id_clock_setup
  character(len=40)  :: mod = "MOM_kappa_shear"  ! This module's name.

#undef  DEBUG
#undef  ADD_DIAGNOSTICS

contains

subroutine Calculate_kappa_shear(u_in, v_in, h, tv, p_surf, kappa_io, tke_io, &
                                 kv_io, dt, G, GV, CS, initialize_all)
  type(ocean_grid_type),                      intent(in)    :: G
  type(verticalGrid_type),                    intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: u_in
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: v_in
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: h
  type(thermo_var_ptrs),                      intent(in)    :: tv
  real, dimension(:,:),                       pointer       :: p_surf
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: kappa_io
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: tke_io
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: kv_io ! really intent(out)
  real,                                       intent(in)    :: dt
  type(Kappa_shear_CS),                       pointer       :: CS
  logical,                          optional, intent(in)    :: initialize_all
!
! ----------------------------------------------
! Subroutine for calculating diffusivity and TKE
! ----------------------------------------------
! Arguments: u_in - Initial zonal velocity, in m s-1. (Intent in)
!  (in)      v_in - Initial meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      p_surf - The pressure at the ocean surface in Pa (or NULL).
!  (in/out)  kappa_io - The diapycnal diffusivity at each interface
!                       (not layer!) in m2 s-1.  Initially this is the value
!                       from the previous timestep, which may accelerate the
!                       iteration toward convergence.
!  (in/out)  tke_io - The turbulent kinetic energy per unit mass at each
!                     interface (not layer!) in m2 s-2.  Initially this is the
!                     value from the previous timestep, which may accelerate
!                     the iteration toward convergence.
!  (in/out)  kv_io - The vertical viscosity at each interface
!                    (not layer!) in m2 s-1. This discards any previous value
!                    i.e. intent(out) and simply sets Kv = Prandtl * Kd_turb
!  (in)      dt - Time increment, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 kappa_shear_init.
!  (in,opt)  initialize_all - If present and false, the previous value of
!                             kappa is used to start the iterations.
  real, dimension(SZI_(G),SZK_(G)) :: &
    h_2d, &                         ! A 2-D version of h, but converted to m.
    u_2d, v_2d, T_2d, S_2d, rho_2d  ! 2-D versions of u_in, v_in, T, S, and rho.
  real, dimension(SZI_(G),SZK_(G)+1) :: &
    kappa_2d, tke_2d              ! 2-D versions of various kappa_io and tke_io.
  real, dimension(SZK_(G)) :: &
    u, &        ! The zonal velocity after a timestep of mixing, in m s-1.
    v, &        ! The meridional velocity after a timestep of mixing, in m s-1.
    Idz, &      ! The inverse of the distance between TKE points, in m.
    T, &        ! The potential temperature after a timestep of mixing, in C.
    Sal, &      ! The salinity after a timestep of mixing, in psu.
    dz, &       ! The layer thickness, in m.
    u0xdz, &    ! The initial zonal velocity times dz, in m2 s-1.
    v0xdz, &    ! The initial meridional velocity times dz, in m2 s-1.
    T0xdz, &    ! The initial temperature times dz, in C m.
    S0xdz, &    ! The initial salinity times dz, in PSU m.
    u_test, v_test, T_test, S_test
  real, dimension(SZK_(G)+1) :: &
    N2, &       ! The squared buoyancy frequency at an interface, in s-2.
    dz_Int, &   ! The extent of a finite-volume space surrounding an interface,
                ! as used in calculating kappa and TKE, in m.
    I_dz_int, & ! The inverse of the distance between velocity & density points
                ! above and below an interface, in m-1.  This is used to
                ! calculate N2, shear, and fluxes, and it might differ from
                ! 1/dz_Int, as they have different uses.
    S2, &       ! The squared shear at an interface, in s-2.
    a1, &       ! a1 is the coupling between adjacent interfaces in the TKE,
                ! velocity, and density equations, in m s-1 or m.
    c1, &       ! c1 is used in the tridiagonal (and similar) solvers.
    k_src, &    ! The shear-dependent source term in the kappa equation, in s-1.
    kappa_src, & ! The shear-dependent source term in the kappa equation in s-1.
    kappa, &    ! The shear-driven diapycnal diffusivity at an interface, in
                ! units of m2 s-1.
    tke, &      ! The Turbulent Kinetic Energy per unit mass at an interface,
                ! in units of m2 s-2.
    kappa_avg, & ! The time-weighted average of kappa, in m2 s-1.
    tke_avg, &  ! The time-weighted average of TKE, in m2 s-2.
    kappa_out, & ! The kappa that results from the kappa equation, in m2 s-1.
    kappa_mid, & ! The average of the initial and predictor estimates of kappa,
                ! in units of m2 s-1.
    tke_pred, & ! The value of TKE from a predictor step, in m2 s-2.
    kappa_pred, & ! The value of kappa from a predictor step, in m2 s-1.
    pressure, & ! The pressure at an interface, in Pa.
    T_int, &    ! The temperature interpolated to an interface, in C.
    Sal_int, &  ! The salinity interpolated to an interface, in psu.
    dbuoy_dT, & ! The partial derivatives of buoyancy with changes in
    dbuoy_dS, & ! temperature and salinity, in m s-2 K-1 and m s-2 psu-1.
    I_L2_bdry, &   ! The inverse of the square of twice the harmonic mean
                   ! distance to the top and bottom boundaries, in m-2.
    K_Q, &         ! Diffusivity divided by TKE, in s.
    K_Q_tmp, &     ! Diffusivity divided by TKE, in s.
    local_src_avg, & ! The time-integral of the local source, nondim.
    tol_min, & ! Minimum tolerated ksrc for the corrector step, in s-1.
    tol_max, & ! Maximum tolerated ksrc for the corrector step, in s-1.
    tol_chg, & ! The tolerated change integrated in time, nondim.
    dist_from_top, &  ! The distance from the top surface, in m.
    local_src     ! The sum of all sources of kappa, including kappa_src and
                  ! sources from the elliptic term, in s-1.
  real :: f2   ! The squared Coriolis parameter of each column, in s-2.
  real :: dist_from_bot ! The distance from the bottom surface, in m.

  real :: b1            ! The inverse of the pivot in the tridiagonal equations.
  real :: bd1           ! A term in the denominator of b1.
  real :: d1            ! 1 - c1 in the tridiagonal equations.
  real :: gR0           ! Rho_0 times g in kg m-2 s-2.
  real :: g_R0          ! g_R0 is g/Rho in m4 kg-1 s-2.
  real :: Norm          ! A factor that normalizes two weights to 1, in m-2.
  real :: tol_dksrc, tol2  ! ### Tolerances that need to be set better later.
  real :: tol_dksrc_low ! The tolerance for the fractional decrease in ksrc
                        ! within an iteration.  0 < tol_dksrc_low < 1.
  real :: Ri_crit       !   The critical shear Richardson number for shear-
                        ! driven mixing. The theoretical value is 0.25.
  real :: dt_rem        !   The remaining time to advance the solution, in s.
  real :: dt_now        !   The time step used in the current iteration, in s.
  real :: dt_wt         !   The fractional weight of the current iteration, ND.
  real :: dt_test       !   A time-step that is being tested for whether it
                        ! gives acceptably small changes in k_src, in s.
  real :: Idtt          !   Idtt = 1 / dt_test, in s-1.
  real :: dt_inc        !   An increment to dt_test that is being tested, in s.

  real :: dz_in_lay     !   The running sum of the thickness in a layer, in m.
  real :: k0dt          ! The background diffusivity times the timestep, in m2.
  real :: dz_massless   ! A layer thickness that is considered massless, in m.
  logical :: valid_dt   ! If true, all levels so far exhibit acceptably small
                        ! changes in k_src.
  logical :: use_temperature  !  If true, temperature and salinity have been
                        ! allocated and are being used as state variables.
  logical :: new_kappa = .true. ! If true, ignore the value of kappa from the
                        ! last call to this subroutine.

  integer, dimension(SZK_(G)+1) :: kc ! The index map between the original
                        ! interfaces and the interfaces with massless layers
                        ! merged into nearby massive layers.
  real, dimension(SZK_(G)+1) :: kf ! The fractional weight of interface kc+1 for
                        ! interpolating back to the original index space, ND.
  integer :: ks_kappa, ke_kappa  ! The k-range with nonzero kappas.
  integer :: dt_halvings   ! The number of times that the time-step is halved
                           ! in seeking an acceptable timestep.  If none is
                           ! found, dt_rem*0.5^dt_halvings is used.
  integer :: dt_refinements ! The number of 2-fold refinements that will be used
                           ! to estimate the maximum permitted time step.  I.e.,
                           ! the resolution is 1/2^dt_refinements.
  integer :: is, ie, js, je, i, j, k, nz, nzc, itt, itt_dt

  ! Arrays that are used in find_kappa_tke if it is included in this subroutine.
  real, dimension(SZK_(G)) :: &
    aQ, &       ! aQ is the coupling between adjacent interfaces in the TKE
                ! equations, in m s-1.
    dQdz        ! Half the partial derivative of TKE with depth, m s-2.
  real, dimension(SZK_(G)+1) :: &
    dK, &         ! The change in kappa, in m2 s-1.
    dQ, &         ! The change in TKE, in m2 s-1.
    cQ, cK, &     ! cQ and cK are the upward influences in the tridiagonal and
                  ! hexadiagonal solvers for the TKE and kappa equations, ND.
    I_Ld2, &      ! 1/Ld^2, where Ld is the effective decay length scale
                  ! for kappa, in units of m-2.
    TKE_decay, &  ! The local TKE decay rate in s-1.
    dQmdK, &      ! With Newton's method the change in dQ(K-1) due to dK(K), s.
    dKdQ, &       ! With Newton's method the change in dK(K) due to dQ(K), s-1.
    e1            ! The fractional change in a layer TKE due to a change in the
                  ! TKE of the layer above when all the kappas below are 0.
                  ! e1 is nondimensional, and 0 < e1 < 1.


  ! Diagnostics that should be deleted?
#ifdef ADD_DIAGNOSTICS
  real, dimension(SZK_(G)+1) :: &  ! Additional diagnostics.
    I_Ld2_1d
  real, dimension(SZI_(G),SZK_(G)+1) :: & ! 2-D versions of diagnostics.
    I_Ld2_2d, dz_Int_2d
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: & ! 3-D versions of diagnostics.
    I_Ld2_3d, dz_Int_3d
#endif
#ifdef DEBUG
  integer :: max_debug_itt ; parameter(max_debug_itt=20)
  real :: wt(SZK_(G)+1), wt_tot, I_wt_tot, wt_itt
  real, dimension(SZK_(G)+1) :: &
    Ri_k, tke_prev, dtke, dkappa, dtke_norm, &
    ksrc_av    ! The average through the iterations of k_src, in s-1.
  real, dimension(SZK_(G)+1,0:max_debug_itt) :: &
    tke_it1, N2_it1, Sh2_it1, ksrc_it1, kappa_it1, kprev_it1
  real, dimension(SZK_(G)+1,1:max_debug_itt) :: &
    dkappa_it1, wt_it1, K_Q_it1, d_dkappa_it1, dkappa_norm
  real, dimension(SZK_(G),0:max_debug_itt) :: &
    u_it1, v_it1, rho_it1, T_it1, S_it1
  real, dimension(0:max_debug_itt) :: &
    dk_wt_it1, dkpos_wt_it1, dkneg_wt_it1, k_mag
  real, dimension(max_debug_itt) ::  dt_it1
#endif
  is = G%isc ; ie = G%iec; js = G%jsc ; je = G%jec ; nz = G%ke

  ! These are hard-coded for now.  Perhaps these could be made dynamic later?
  ! tol_dksrc = 0.5*tol_ksrc_chg ; tol_dksrc_low = 1.0 - 1.0/tol_ksrc_chg ?
  tol_dksrc = 10.0 ; tol_dksrc_low = 0.95 ; tol2 = 2.0*CS%kappa_tol_err
  dt_refinements = 5 ! Selected so that 1/2^dt_refinements < 1-tol_dksrc_low

  use_temperature = .false. ; if (associated(tv%T)) use_temperature = .true.
  new_kappa = .true. ; if (present(initialize_all)) new_kappa = initialize_all

  Ri_crit = CS%Rino_crit
  gR0 = GV%Rho0*GV%g_Earth ; g_R0 = GV%g_Earth/GV%Rho0

  k0dt = dt*CS%kappa_0
  dz_massless = 0.1*sqrt(k0dt)

!$OMP parallel do default(none) shared(js,je,is,ie,nz,h,u_in,v_in,use_temperature,new_kappa, &
!$OMP                                  tv,G,GV,CS,kappa_io,dz_massless,k0dt,p_surf,gR0,g_R0, &
!$OMP                                  dt,tol_dksrc,tol_dksrc_low,tol2,dt_refinements,       &
#ifdef ADD_DIAGNOSTICS
!$OMP                                  I_Ld2_3d,dz_Int_3d,                                   &
#endif
!$OMP                                  Ri_crit,tke_io,kv_io)                                 &
!$OMP                          private(h_2d,u_2d,v_2d,T_2d,S_2d,rho_2d,kappa_2d,nzc,dz,      &
!$OMP                                  u0xdz,v0xdz,T0xdz,S0xdz,kc,Idz,kf,I_dz_int,dz_in_lay, &
!$OMP                                  dist_from_top,a1,b1,u,v,T,Sal,c1,d1,bd1,dz_Int,Norm,  &
!$OMP                                  dist_from_bot,I_L2_bdry,f2,pressure,T_int,Sal_int,    &
!$OMP                                  dbuoy_dT,dbuoy_dS,kappa,K_Q,N2,S2,dt_rem,kappa_avg,   &
!$OMP                                  tke_avg,local_src_avg,tke,kappa_out,kappa_src,        &
!$OMP                                  local_src,ks_kappa,ke_kappa,dt_now,dt_test,tol_max,   &
!$OMP                                  tol_min,tol_chg,u_test, v_test, T_test, S_test,       &
!$OMP                                  valid_dt,Idtt,k_src,dt_inc,dt_wt,kappa_mid,K_Q_tmp,   &
#ifdef ADD_DIAGNOSTICS
!$OMP                                  I_Ld2_1d,I_Ld2_2d, dz_Int_2d,                         &
#endif
!$OMP                                  tke_pred,kappa_pred,tke_2d)

  do j=js,je
    do k=1,nz ; do i=is,ie
      h_2d(i,k) = h(i,j,k)*GV%H_to_m
      u_2d(i,k) = u_in(i,j,k) ; v_2d(i,k) = v_in(i,j,k)
    enddo ; enddo
    if (use_temperature) then ; do k=1,nz ; do i=is,ie
      T_2d(i,k) = tv%T(i,j,k) ; S_2d(i,k) = tv%S(i,j,k)
    enddo ; enddo ; else ; do k=1,nz ; do i=is,ie
      rho_2d(i,k) = GV%Rlay(k) ! Could be tv%Rho(i,j,k) ?
    enddo ; enddo ; endif
    if (.not.new_kappa) then ; do K=1,nz+1 ; do i=is,ie
      kappa_2d(i,K) = kappa_io(i,j,K)
    enddo ; enddo ; endif

!---------------------------------------
! Work on each column.
!---------------------------------------
    do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
    ! call cpu_clock_begin(id_clock_setup)
      ! Store a transposed version of the initial arrays.
      ! Any elimination of massless layers would occur here.
      if (CS%eliminate_massless) then
        nzc = 1
        do k=1,nz
          ! Zero out the thicknesses of all layers, even if they are unused.
          dz(k) = 0.0 ; u0xdz(k) = 0.0 ; v0xdz(k) = 0.0
          T0xdz(k) = 0.0 ; S0xdz(k) = 0.0

          ! Add a new layer if this one has mass.
!          if ((dz(nzc) > 0.0) .and. (h_2d(i,k) > dz_massless)) nzc = nzc+1
          if ((k>CS%nkml) .and. (dz(nzc) > 0.0) .and. &
              (h_2d(i,k) > dz_massless)) nzc = nzc+1

          ! Only merge clusters of massless layers.
!         if ((dz(nzc) > dz_massless) .or. &
!             ((dz(nzc) > 0.0) .and. (h_2d(i,k) > dz_massless))) nzc = nzc+1

          kc(k) = nzc
          dz(nzc) = dz(nzc) + h_2d(i,k)
          u0xdz(nzc) = u0xdz(nzc) + u_2d(i,k)*h_2d(i,k)
          v0xdz(nzc) = v0xdz(nzc) + v_2d(i,k)*h_2d(i,k)
          if (use_temperature) then
            T0xdz(nzc) = T0xdz(nzc) + T_2d(i,k)*h_2d(i,k)
            S0xdz(nzc) = S0xdz(nzc) + S_2d(i,k)*h_2d(i,k)
          else
            T0xdz(nzc) = T0xdz(nzc) + rho_2d(i,k)*h_2d(i,k)
            S0xdz(nzc) = S0xdz(nzc) + rho_2d(i,k)*h_2d(i,k)
          endif
        enddo
        kc(nz+1) = nzc+1
        ! Set up Idz as the inverse of layer thicknesses.
        do k=1,nzc ; Idz(k) = 1.0 / dz(k) ; enddo

        !   Now determine kf, the fractional weight of interface kc when
        ! interpolating between interfaces kc and kc+1.
        kf(1) = 0.0 ; dz_in_lay = h_2d(i,1)
        do k=2,nz
          if (kc(k) > kc(k-1)) then
            kf(k) = 0.0 ; dz_in_lay = h_2d(i,k)
          else
            kf(k) = dz_in_lay*Idz(kc(k)) ; dz_in_lay = dz_in_lay + h_2d(i,k)
          endif
        enddo
        kf(nz+1) = 0.0
      else
        do k=1,nz
          dz(k) = h_2d(i,k)
          u0xdz(k) = u_2d(i,k)*dz(k) ; v0xdz(k) = v_2d(i,k)*dz(k)
        enddo
        if (use_temperature) then
          do k=1,nz
            T0xdz(k) = T_2d(i,k)*dz(k) ; S0xdz(k) = S_2d(i,k)*dz(k)
          enddo
        else
          do k=1,nz
            T0xdz(k) = rho_2d(i,k)*dz(k) ; S0xdz(k) = rho_2d(i,k)*dz(k)
          enddo
        endif
        nzc = nz
        do k=1,nzc+1 ; kc(k) = k ; kf(k) = 0.0 ; enddo
        ! Set up Idz as the inverse of layer thicknesses.
        do k=1,nzc ; Idz(k) = 1.0 / dz(k) ; enddo
      endif

      !   Set up I_dz_int as the inverse of the distance between
      ! adjacent layer centers.
      I_dz_int(1) = 2.0 / dz(1)
      dist_from_top(1) = 0.0
      do K=2,nzc
        I_dz_int(K) = 2.0 / (dz(k-1) + dz(k))
        dist_from_top(K) = dist_from_top(K-1) + dz(k-1)
      enddo
      I_dz_int(nzc+1) = 2.0 / dz(nzc)

      !   Determine the velocities and thicknesses after eliminating massless
      ! layers and applying a time-step of background diffusion.
      if (nzc > 1) then
        a1(2) = k0dt*I_dz_int(2)
        b1 = 1.0 / (dz(1)+a1(2))
        u(1) = b1 * u0xdz(1) ; v(1) = b1 * v0xdz(1)
        T(1) = b1 * T0xdz(1) ; Sal(1) = b1 * S0xdz(1)
        c1(2) = a1(2) * b1 ; d1 = dz(1) * b1 ! = 1 - c1
        do k=2,nzc-1
          bd1 = dz(k) + d1*a1(k)
          a1(k+1) = k0dt*I_dz_int(k+1)
          b1 = 1.0 / (bd1 + a1(k+1))
          u(k) = b1 * (u0xdz(k) + a1(k)*u(k-1))
          v(k) = b1 * (v0xdz(k) + a1(k)*v(k-1))
          T(k) = b1 * (T0xdz(k) + a1(k)*T(k-1))
          Sal(k) = b1 * (S0xdz(k) + a1(k)*Sal(k-1))
          c1(k+1) = a1(k+1) * b1 ; d1 = bd1 * b1 ! d1 = 1 - c1
        enddo
        ! rho or T and S have insulating boundary conditions, u & v use no-slip
        ! bottom boundary conditions (if kappa0 > 0).
        ! For no-slip bottom boundary conditions
        b1 = 1.0 / ((dz(nzc) + d1*a1(nzc)) + k0dt*I_dz_int(nzc+1))
        u(nzc) = b1 * (u0xdz(nzc) + a1(nzc)*u(nzc-1))
        v(nzc) = b1 * (v0xdz(nzc) + a1(nzc)*v(nzc-1))
        ! For insulating boundary conditions
        b1 = 1.0 / (dz(nzc) + d1*a1(nzc))
        T(nzc) = b1 * (T0xdz(nzc) + a1(nzc)*T(nzc-1))
        Sal(nzc) = b1 * (S0xdz(nzc) + a1(nzc)*Sal(nzc-1))
        do k=nzc-1,1,-1
          u(k) = u(k) + c1(k+1)*u(k+1) ; v(k) = v(k) + c1(k+1)*v(k+1)
          T(k) = T(k) + c1(k+1)*T(k+1) ; Sal(k) = Sal(k) + c1(k+1)*Sal(k+1)
        enddo
      else
        ! This is correct, but probably unnecessary.
        b1 = 1.0 / (dz(1) + k0dt*I_dz_int(2))
        u(1) = b1 * u0xdz(1) ; v(1) = b1 * v0xdz(1)
        b1 = 1.0 / dz(1)
        T(1) = b1 * T0xdz(1) ; Sal(1) = b1 * S0xdz(1)
      endif

      ! This uses half the harmonic mean of thicknesses to provide two estimates
      ! of the boundary between cells, and the inverse of the harmonic mean to
      ! weight the two estimates.  The net effect is that interfaces around thin
      ! layers have thin cells, and the total thickness adds up properly.
      ! The top- and bottom- interfaces have zero thickness, consistent with
      ! adding additional zero thickness layers.
      dz_Int(1) = 0.0 ; dz_Int(2) = dz(1)
      do K=2,nzc-1
        Norm = 1.0 / (dz(k)*(dz(k-1)+dz(k+1)) + 2.0*dz(k-1)*dz(k+1))
        dz_Int(K) = dz_Int(K) + dz(k) * ( ((dz(k)+dz(k+1)) * dz(k-1)) * Norm)
        dz_Int(K+1) = dz(k) * ( ((dz(k-1)+dz(k)) * dz(k+1)) * Norm)
      enddo
      dz_Int(nzc) = dz_Int(nzc) + dz(nzc) ; dz_Int(nzc+1) = 0.0

#ifdef ADD_DIAGNOSTICS
      do K=1,nzc+1 ; I_Ld2_1d(K) = 0.0 ; enddo
#endif

      dist_from_bot = 0.0
      do K=nzc,2,-1
        dist_from_bot = dist_from_bot + dz(k)
        I_L2_bdry(K) = (dist_from_top(K) + dist_from_bot)**2 / &
                         (dist_from_top(K) * dist_from_bot)**2
      enddo
      f2 = 0.25*((G%CoriolisBu(I,j)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                 (G%CoriolisBu(I,J-1)**2 + G%CoriolisBu(I-1,J)**2))

      ! Calculate thermodynamic coefficients and an initial estimate of N2.
      if (use_temperature) then
        pressure(1) = 0.0
        if (associated(p_surf)) pressure(1) = p_surf(i,j)
        do K=2,nzc
          pressure(K) = pressure(K-1) + gR0*dz(k-1)
          T_int(K) = 0.5*(T(k-1) + T(k))
          Sal_int(K) = 0.5*(Sal(k-1) + Sal(k))
        enddo
        call calculate_density_derivs(T_int, Sal_int, pressure, dbuoy_dT, &
                                      dbuoy_dS, 2, nzc-1, tv%eqn_of_state)
        do K=2,nzc
          dbuoy_dT(K) = -G_R0*dbuoy_dT(K)
          dbuoy_dS(K) = -G_R0*dbuoy_dS(K)
        enddo
      else
        do K=1,nzc+1 ; dbuoy_dT(K) = -G_R0 ; dbuoy_dS(K) = 0.0 ; enddo
      endif

    ! ----------------------------------------------------
    ! Calculate kappa, here defined at interfaces.
    ! ----------------------------------------------------
      if (new_kappa) then
        do K=1,nzc+1 ; kappa(K) = 1.0 ; K_Q(K) = 0.0 ; enddo
      else
        do K=1,nzc+1 ; kappa(K) = kappa_2d(i,K) ; K_Q(K) = 0.0 ; enddo
      endif

#ifdef DEBUG
      N2(1) = 0.0 ; N2(nzc+1) = 0.0
      do K=2,nzc
        N2(K) = max((dbuoy_dT(K) * (T0xdz(k-1)*Idz(k-1) - T0xdz(k)*Idz(k)) + &
                     dbuoy_dS(K) * (S0xdz(k-1)*Idz(k-1) - S0xdz(k)*Idz(k))) * &
                     I_dz_int(K), 0.0)
      enddo
      do k=1,nzc
        u_it1(k,0) = u0xdz(k)*Idz(k) ; v_it1(k,0) = v0xdz(k)*Idz(k)
        T_it1(k,0) = T0xdz(k)*Idz(k) ; S_it1(k,0) = S0xdz(k)*Idz(k)
      enddo
      do K=1,nzc+1
        kprev_it1(K,0) = kappa(K) ; kappa_it1(K,0) = kappa(K)
        tke_it1(K,0) = tke(K)
        N2_it1(K,0) = N2(K) ; Sh2_it1(K,0) = S2(K) ; ksrc_it1(K,0) = k_src(K)
      enddo
      do k=nzc+1,nz
        u_it1(k,0) = 0.0 ; v_it1(k,0) = 0.0
        T_it1(k,0) = 0.0 ; S_it1(k,0) = 0.0
        kprev_it1(K+1,0) = 0.0 ; kappa_it1(K+1,0) = 0.0 ; tke_it1(K+1,0) = 0.0
        N2_it1(K+1,0) = 0.0 ; Sh2_it1(K+1,0) = 0.0 ; ksrc_it1(K+1,0) = 0.0
      enddo
      do itt=1,max_debug_itt
        dt_it1(itt) = 0.0
        do k=1,nz
          u_it1(k,itt) = 0.0 ; v_it1(k,itt) = 0.0
          T_it1(k,itt) = 0.0 ; S_it1(k,itt) = 0.0
          rho_it1(k,itt) = 0.0
        enddo
        do K=1,nz+1
          kprev_it1(K,itt) = 0.0 ; kappa_it1(K,itt) = 0.0 ; tke_it1(K,itt) = 0.0
          N2_it1(K,itt) = 0.0 ; Sh2_it1(K,itt) = 0.0
          ksrc_it1(K,itt) = 0.0
          dkappa_it1(K,itt) = 0.0 ; wt_it1(K,itt) = 0.0
          K_Q_it1(K,itt) = 0.0 ; d_dkappa_it1(K,itt) = 0.0
        enddo
      enddo
      do K=1,nz+1 ; ksrc_av(K) = 0.0 ; enddo
#endif

      ! This call just calculates N2 and S2.
      call calculate_projected_state(kappa, u, v, T, Sal, 0.0, nzc, &
                                     dz, I_dz_int, dbuoy_dT, dbuoy_dS, &
                                     u, v, T, Sal, N2=N2, S2=S2)
    ! ----------------------------------------------------
    ! Iterate
    ! ----------------------------------------------------
      dt_rem = dt
      do K=1,nzc+1
        kappa_avg(K) = 0.0 ; tke_avg(K) = 0.0
        local_src_avg(K) = 0.0
        ! Use the grid spacings to scale errors in the source.
        if ( dz_Int(K) > 0.0 ) &
          local_src_avg(K) = 0.1 * k0dt * I_dz_int(K) / dz_Int(K)
      enddo

    ! call cpu_clock_end(id_clock_setup)

!       do itt=1,CS%max_RiNo_it
      do itt=1,CS%max_KS_it

    ! ----------------------------------------------------
    ! Calculate new values of u, v, rho, N^2 and S.
    ! ----------------------------------------------------
#ifdef DEBUG
        do K=1,nzc+1
          Ri_k(K) = 1e3 ; if (S2(K) > 1e-3*N2(K)) Ri_k(K) = N2(K) / S2(K)
          tke_prev(K) = tke(K)
        enddo
#endif

      ! call cpu_clock_begin(id_clock_KQ)
        call find_kappa_tke(N2, S2, kappa, Idz, dz_Int, I_L2_bdry, f2, &
                            nzc, CS, K_Q, tke, kappa_out, kappa_src, local_src)
      ! call cpu_clock_end(id_clock_KQ)

      ! call cpu_clock_begin(id_clock_avg)
        ! Determine the range of non-zero values of kappa_out.
        ks_kappa = nz+1 ; ke_kappa = 0
        do K=2,nzc ; if (kappa_out(K) > 0.0) then
          ks_kappa = K ; exit
        endif ; enddo
        do k=nzc,ks_kappa,-1 ; if (kappa_out(K) > 0.0) then
          ke_kappa = K ; exit
        endif ; enddo
        if (ke_kappa == nzc) kappa_out(nzc+1) = 0.0
      ! call cpu_clock_end(id_clock_avg)

        ! Determine how long to use this value of kappa (dt_now).

      ! call cpu_clock_begin(id_clock_project)
        if ((ke_kappa < ks_kappa) .or. (itt==CS%max_RiNo_it)) then
          dt_now = dt_rem
        else
          ! Limit dt_now so that |k_src(k)-kappa_src(k)| < tol * local_src(k)
          dt_test = dt_rem
          do K=2,nzc
            tol_max(K) = kappa_src(K) + tol_dksrc * local_src(K)
            tol_min(K) = kappa_src(K) - tol_dksrc_low * local_src(K)
            tol_chg(K) = tol2 * local_src_avg(K)
          enddo

          do itt_dt=1,(CS%max_KS_it+1-itt)/2
            !   The maximum number of times that the time-step is halved in
            ! seeking an acceptable timestep is reduced with each iteration,
            ! so that as the maximum number of iterations is approached, the
            ! whole remaining timestep is used.  Typically, an acceptable
            ! timestep is found long before the minimum is reached, so the
            ! value of max_KS_it may be unimportant, especially if it is large
            ! enough.
            call calculate_projected_state(kappa_out, u, v, T, Sal, 0.5*dt_test, nzc, &
                                           dz, I_dz_int, dbuoy_dT, dbuoy_dS, &
                                           u_test, v_test, T_test, S_test, N2, S2, &
                                           ks_int = ks_kappa, ke_int = ke_kappa)
            valid_dt = .true.
            Idtt = 1.0 / dt_test
            do K=max(ks_kappa-1,2),min(ke_kappa+1,nzc)
              if (N2(K) < Ri_crit * S2(K)) then ! Equivalent to Ri < Ri_crit.
                k_src(K) = (2.0 * CS%Shearmix_rate * sqrt(S2(K))) * &
                           ((Ri_crit*S2(K) - N2(K)) / (Ri_crit*S2(K) + CS%FRi_curvature*N2(K)))
                 if ((k_src(K) > max(tol_max(K), kappa_src(K) + Idtt*tol_chg(K))) .or. &
                     (k_src(K) < min(tol_min(K), kappa_src(K) - Idtt*tol_chg(K)))) then
                  valid_dt = .false. ; exit
                endif
              else
                if (0.0 < min(tol_min(K), kappa_src(K) - Idtt*tol_chg(K))) then
                  valid_dt = .false. ; k_src(K) = 0.0 ; exit
                endif
              endif
            enddo

            if (valid_dt) exit
            dt_test = 0.5*dt_test
          enddo
          if ((dt_test < dt_rem) .and. valid_dt) then
            dt_inc = 0.5*dt_test
            do itt_dt=1,dt_refinements
              call calculate_projected_state(kappa_out, u, v, T, Sal, &
                       0.5*(dt_test+dt_inc), nzc, dz, I_dz_int, dbuoy_dT, &
                       dbuoy_dS, u_test, v_test, T_test, S_test, N2, S2, &
                       ks_int = ks_kappa, ke_int = ke_kappa)
              valid_dt = .true.
              Idtt = 1.0 / (dt_test+dt_inc)
              do K=max(ks_kappa-1,2),min(ke_kappa+1,nzc)
                if (N2(K) < Ri_crit * S2(K)) then ! Equivalent to Ri < Ri_crit.
                  k_src(K) = (2.0 * CS%Shearmix_rate * sqrt(S2(K))) * &
                             ((Ri_crit*S2(K) - N2(K)) / &
                              (Ri_crit*S2(K) + CS%FRi_curvature*N2(K)))
                   if ((k_src(K) > max(tol_max(K), kappa_src(K) + Idtt*tol_chg(K))) .or. &
                       (k_src(K) < min(tol_min(K), kappa_src(K) - Idtt*tol_chg(K)))) then
                    valid_dt = .false. ; exit
                  endif
                else
                  if (0.0 < min(tol_min(K), kappa_src(K) - Idtt*tol_chg(K))) then
                    valid_dt = .false. ; k_src(K) = 0.0 ; exit
                  endif
                endif
              enddo

              if (valid_dt) dt_test = dt_test + dt_inc
              dt_inc = 0.5*dt_inc
            enddo
          else
            dt_inc = 0.0
          endif

           dt_now = min(dt_test*(1.0+CS%kappa_tol_err)+dt_inc,dt_rem)
          do K=2,nzc
            local_src_avg(K) = local_src_avg(K) + dt_now * local_src(K)
          enddo
        endif  ! Are all the values of kappa_out 0?
      ! call cpu_clock_end(id_clock_project)

        ! The state has already been projected forward. Now find new values of kappa.

        if (ke_kappa < ks_kappa) then
          ! There is no mixing now, and will not be again.
        ! call cpu_clock_begin(id_clock_avg)
          dt_wt = dt_rem / dt ; dt_rem = 0.0
          do K=1,nzc+1
            kappa_mid(K) = 0.0
            ! This would be here but does nothing.
            ! kappa_avg(K) = kappa_avg(K) + kappa_mid(K)*dt_wt
            tke_avg(K) = tke_avg(K) + dt_wt*tke(K)
#ifdef DEBUG
            tke_pred(K) = tke(K) ; kappa_pred(K) = 0.0 ; kappa(K) = 0.0
#endif
          enddo
        ! call cpu_clock_end(id_clock_avg)
        else
        ! call cpu_clock_begin(id_clock_project)
          call calculate_projected_state(kappa_out, u, v, T, Sal, dt_now, nzc, &
                                         dz, I_dz_int, dbuoy_dT, dbuoy_dS, &
                                         u_test, v_test, T_test, S_test, N2=N2, S2=S2, &
                                         ks_int = ks_kappa, ke_int = ke_kappa)
        ! call cpu_clock_end(id_clock_project)

        ! call cpu_clock_begin(id_clock_KQ)
          do K=1,nzc+1 ; K_Q_tmp(K) = K_Q(K) ; enddo
          call find_kappa_tke(N2, S2, kappa_out, Idz, dz_Int, I_L2_bdry, f2, &
                              nzc, CS, K_Q_tmp, tke_pred, kappa_pred)
        ! call cpu_clock_end(id_clock_KQ)

          ks_kappa = nz+1 ; ke_kappa = 0
          do K=1,nzc+1
            kappa_mid(K) = 0.5*(kappa_out(K) + kappa_pred(K))
            if ((kappa_mid(K) > 0.0) .and. (K<ks_kappa)) ks_kappa = K
            if (kappa_mid(K) > 0.0) ke_kappa = K
          enddo

        ! call cpu_clock_begin(id_clock_project)
          call calculate_projected_state(kappa_mid, u, v, T, Sal, dt_now, nzc, &
                                         dz, I_dz_int, dbuoy_dT, dbuoy_dS, &
                                         u_test, v_test, T_test, S_test, N2=N2, S2=S2, &
                                         ks_int = ks_kappa, ke_int = ke_kappa)
        ! call cpu_clock_end(id_clock_project)

        ! call cpu_clock_begin(id_clock_KQ)
          call find_kappa_tke(N2, S2, kappa_out, Idz, dz_Int, I_L2_bdry, f2, &
                              nzc, CS, K_Q, tke_pred, kappa_pred)
        ! call cpu_clock_end(id_clock_KQ)

        ! call cpu_clock_begin(id_clock_avg)
          dt_wt = dt_now / dt ; dt_rem = dt_rem - dt_now
          do K=1,nzc+1
            kappa_mid(K) = 0.5*(kappa_out(K) + kappa_pred(K))
            kappa_avg(K) = kappa_avg(K) + kappa_mid(K)*dt_wt
            tke_avg(K) = tke_avg(K) + dt_wt*0.5*(tke_pred(K) + tke(K))
            kappa(K) = kappa_pred(K) ! First guess for the next iteration.
          enddo
        ! call cpu_clock_end(id_clock_avg)
        endif

        if (dt_rem > 0.0) then
          ! Update the values of u, v, T, Sal, N2, and S2 for the next iteration.
        ! call cpu_clock_begin(id_clock_project)
          call calculate_projected_state(kappa_mid, u, v, T, Sal, dt_now, nzc, &
                                         dz, I_dz_int, dbuoy_dT, dbuoy_dS, &
                                         u, v, T, Sal, N2, S2)
        ! call cpu_clock_end(id_clock_project)
        endif

#ifdef DEBUG
        if (itt <= max_debug_itt) then
          dt_it1(itt) = dt_now
          dk_wt_it1(itt) = 0.0 ; dkpos_wt_it1(itt) = 0.0 ;  dkneg_wt_it1(itt) = 0.0
          k_mag(itt) = 0.0
          wt_itt = 1.0/real(itt) ; wt_tot = 0.0
          do K=1,nzc+1
            ksrc_av(K) = (1.0-wt_itt)*ksrc_av(K) + wt_itt*k_src(K)
            wt_tot = wt_tot + dz_Int(K) * ksrc_av(K)
          enddo
          ! Use Adcroft's 1/0=0 convention.
          I_wt_tot = 0.0 ; if (wt_tot > 0.0) I_wt_tot = 1.0/wt_tot

          do K=1,nzc+1
            wt(K) = (dz_Int(K)*ksrc_av(K)) * I_wt_tot
            k_mag(itt) = k_mag(itt) + wt(K)*kappa_mid(K)
            dkappa_it1(K,itt) = kappa_pred(K) - kappa_out(K)
            dk_wt_it1(itt) = dk_wt_it1(itt) + wt(K)*dkappa_it1(K,itt)
            if (dk > 0.0) then
              dkpos_wt_it1(itt) = dkpos_wt_it1(itt) + wt(K)*dkappa_it1(K,itt)
            else
              dkneg_wt_it1(itt) = dkneg_wt_it1(itt) + wt(K)*dkappa_it1(K,itt)
            endif
            wt_it1(K,itt) = wt(K)
          enddo
        endif
        do K=1,nzc+1
          Ri_k(K) = 1e3 ; if (N2(K) < 1e3 * S2(K)) Ri_k(K) = N2(K) / S2(K)
          dtke(K) = tke_pred(K) - tke(K)
          dtke_norm(K) = dtke(K) / (0.5*(tke(K) + tke_pred(K)))
          dkappa(K) = kappa_pred(K) - kappa_out(K)
        enddo
        if (itt <= max_debug_itt) then
          do k=1,nzc
            u_it1(k,itt) = u(k) ; v_it1(k,itt) = v(k)
            T_it1(k,itt) = T(k) ; S_it1(k,itt) = Sal(k)
          enddo
          do K=1,nzc+1
            kprev_it1(K,itt)=kappa_out(K)
            kappa_it1(K,itt)=kappa_mid(K) ; tke_it1(K,itt) = 0.5*(tke(K)+tke_pred(K))
            N2_it1(K,itt)=N2(K) ; Sh2_it1(K,itt)=S2(K)
            ksrc_it1(K,itt) = kappa_src(K)
            K_Q_it1(K,itt) = kappa_out(K) / (TKE(K))
            if (itt > 1) then
              if (abs(dkappa_it1(K,itt-1)) > 1e-20) &
                d_dkappa_it1(K,itt) = dkappa_it1(K,itt) / dkappa_it1(K,itt-1)
            endif
            dkappa_norm(K,itt) = dkappa(K) / max(0.5*(Kappa_pred(K) + kappa_out(K)), 1e-100)
          enddo
        endif
#endif

        if (dt_rem <= 0.0) exit

      enddo ! end itt loop

    ! call cpu_clock_begin(id_clock_setup)
    ! Extrapolate from the vertically reduced grid back to the original layers.
      if (nz == nzc) then
        do K=1,nz+1
          kappa_2d(i,K) = kappa_avg(K)
          tke_2d(i,K) = tke(K)
        enddo
      else
        do K=1,nz+1
          if (kf(K) == 0.0) then
            kappa_2d(i,K) = kappa_avg(kc(K))
            tke_2d(i,K) = tke_avg(kc(K))
          else
            kappa_2d(i,K) = (1.0-kf(K)) * kappa_avg(kc(K)) + &
                             kf(K) * kappa_avg(kc(K)+1)
            tke_2d(i,K) = (1.0-kf(K)) * tke_avg(kc(K)) + &
                           kf(K) * tke_avg(kc(K)+1)
          endif
        enddo
      endif
#ifdef ADD_DIAGNOSTICS
      I_Ld2_2d(i,1) = 0.0 ; dz_Int_2d(i,1) = dz_Int(1)
      do K=2,nzc
        I_Ld2_2d(i,K) = (N2(K) / CS%lambda**2 + f2) / &
                         max(TKE(K),1e-30) + I_L2_bdry(K)
        dz_Int_2d(i,K) = dz_Int(K)
      enddo
      I_Ld2_2d(i,nzc+1) = 0.0 ; dz_Int_2d(i,nzc+1) = dz_Int(nzc+1)
      do K=nzc+2,nz+1
        I_Ld2_2d(i,K) = 0.0 ; dz_Int_2d(i,K) = 0.0
      enddo
#endif
    ! call cpu_clock_end(id_clock_setup)
    else  ! Land points, still inside the i-loop.
      do K=1,nz+1
        kappa_2d(i,K) = 0.0 ; tke_2d(i,K) = 0.0
#ifdef ADD_DIAGNOSTICS
        I_Ld2_2d(i,K) = 0.0
        dz_Int_2d(i,K) = dz_Int(K)
#endif
      enddo
    endif ; enddo ! i-loop

    do K=1,nz+1 ; do i=is,ie
      kappa_io(i,j,K) = G%mask2dT(i,j) * kappa_2d(i,K)
      tke_io(i,j,K) = G%mask2dT(i,j) * tke_2d(i,K)
      kv_io(i,j,K) = ( G%mask2dT(i,j) * kappa_2d(i,K) ) * CS%Prandtl_turb
#ifdef ADD_DIAGNOSTICS
      I_Ld2_3d(i,j,K) = I_Ld2_2d(i,K)
      dz_Int_3d(i,j,K) = dz_Int_2d(i,K)
#endif
    enddo ; enddo

  enddo ! end of j-loop

  if (CS%debug) then
    call hchksum(kappa_io,"kappa",G%HI)
    call hchksum(tke_io,"tke",G%HI)
  endif

  if (CS%id_Kd_shear > 0) call post_data(CS%id_Kd_shear, kappa_io, CS%diag)
  if (CS%id_TKE > 0) call post_data(CS%id_TKE, tke_io, CS%diag)
#ifdef ADD_DIAGNOSTICS
  if (CS%id_ILd2 > 0) call post_data(CS%id_ILd2, I_Ld2_3d, CS%diag)
  if (CS%id_dz_Int > 0) call post_data(CS%id_dz_Int, dz_Int_3d, CS%diag)
#endif

  return

end subroutine Calculate_kappa_shear

subroutine calculate_projected_state(kappa, u0, v0, T0, S0, dt, nz, &
                                     dz, I_dz_int, dbuoy_dT, dbuoy_dS, &
                                     u, v, T, Sal, N2, S2, ks_int, ke_int)
!   This subroutine calculates the velocities, temperature and salinity that
! the water column will have after mixing for dt with diffusivities kappa.  It
! may also calculate the projected buoyancy frequency and shear.
  integer,               intent(in)  :: nz
  real, dimension(nz+1), intent(in)  :: kappa
  real, dimension(nz),   intent(in)  :: u0, v0, T0, S0, dz
  real, dimension(nz+1), intent(in)  :: I_dz_int, dbuoy_dT, dbuoy_dS
  real,                  intent(in)  :: dt
  real, dimension(nz),   intent(inout) :: u, v, T, Sal
  real, dimension(nz+1), optional, intent(inout) :: N2, S2
  integer, optional,     intent(in)  :: ks_int, ke_int
  ! Arguments: kappa - The diapycnal diffusivity at interfaces, in m2 s-1.
  !  (in)      Sh - The shear at interfaces, in s-1.
  !  (in)      u0 - The initial zonal velocity, in m s-1.
  !  (in)      v0 - The initial meridional velocity, in m s-1.
  !  (in)      T0 - The initial temperature, in C.
  !  (in)      S0 - The initial salinity, in PSU.
  !  (in)      nz - The number of layers (after eliminating massless layers?).
  !  (in)      dz - The grid spacing of layers, in m.
  !  (in)      I_dz_int - The inverse of the layer's thicknesses, in m-1.
  !  (in)      dbuoy_dT - The partial derivative of buoyancy with temperature,
  !                       in m s-2 C-1.
  !  (in)      dbuoy_dS - The partial derivative of buoyancy with salinity,
  !                       in m s-2 PSU-1.
  !  (in)      dt - The time step in s.
  !  (in)      nz - The number of layers to work on.
  !  (out)     u - The zonal velocity after dt, in m s-1.
  !  (out)     v - The meridional velocity after dt, in m s-1.
  !  (in)      T - The temperature after dt, in C.
  !  (in)      Sal - The salinity after dt, in PSU.
  !  (out)     N2 - The buoyancy frequency squared at interfaces, in s-2.
  !  (out)     S2 - The squared shear at interfaces, in s-2.
  !  (in,opt)  ks_int - The topmost k-index with a non-zero diffusivity.
  !  (in,opt)  ke_int - The bottommost k-index with a non-zero diffusivity.

 ! UNCOMMENT THE FOLLOWING IF NOT CONTAINED IN THE OUTER SUBROUTINE.
  real, dimension(nz+1) :: c1
  real :: a_a, a_b, b1, d1, bd1, b1nz_0
  integer :: k, ks, ke

  ks = 1 ; ke = nz
  if (present(ks_int)) ks = max(ks_int-1,1)
  if (present(ke_int)) ke = min(ke_int,nz)

  if (ks > ke) return

  if (dt > 0.0) then
    a_b = dt*(kappa(ks+1)*I_dz_int(ks+1))
    b1 = 1.0 / (dz(ks) + a_b)
    c1(ks+1) = a_b * b1 ; d1 = dz(ks) * b1 ! = 1 - c1

    u(ks) = (b1 * dz(ks))*u0(ks) ; v(ks) = (b1 * dz(ks))*v0(ks)
    T(ks) = (b1 * dz(ks))*T0(ks) ; Sal(ks) = (b1 * dz(ks))*S0(ks)
    do K=ks+1,ke-1
      a_a = a_b
      a_b = dt*(kappa(K+1)*I_dz_int(K+1))
      bd1 = dz(k) + d1*a_a
      b1 = 1.0 / (bd1 + a_b)
      c1(K+1) = a_b * b1 ; d1 = bd1 * b1 ! d1 = 1 - c1

      u(k) = b1 * (dz(k)*u0(k) + a_a*u(k-1))
      v(k) = b1 * (dz(k)*v0(k) + a_a*v(k-1))
      T(k) = b1 * (dz(k)*T0(k) + a_a*T(k-1))
      Sal(k) = b1 * (dz(k)*S0(k) + a_a*Sal(k-1))
    enddo
    !   T and S have insulating boundary conditions, u & v use no-slip
    ! bottom boundary conditions at the solid bottom.

    ! For insulating boundary conditions or mixing simply stopping, use...
    a_a = a_b
    b1 = 1.0 / (dz(ke) + d1*a_a)
    T(ke) = b1 * (dz(ke)*T0(ke) + a_a*T(ke-1))
    Sal(ke) = b1 * (dz(ke)*S0(ke) + a_a*Sal(ke-1))

    !   There is no distinction between the effective boundary conditions for
    ! tracers and velocities if the mixing is separated from the bottom, but if
    ! the mixing goes all the way to the bottom, use no-slip BCs for velocities.
    if (ke == nz) then
      a_b = dt*(kappa(nz+1)*I_dz_int(nz+1))
      b1nz_0 = 1.0 / ((dz(nz) + d1*a_a) + a_b)
    else
      b1nz_0 = b1
    endif
    u(ke) = b1nz_0 * (dz(ke)*u0(ke) + a_a*u(ke-1))
    v(ke) = b1nz_0 * (dz(ke)*v0(ke) + a_a*v(ke-1))

    do k=ke-1,ks,-1
      u(k) = u(k) + c1(k+1)*u(k+1)
      v(k) = v(k) + c1(k+1)*v(k+1)
      T(k) = T(k) + c1(k+1)*T(k+1)
      Sal(k) = Sal(k) + c1(k+1)*Sal(k+1)
    enddo
  else ! dt <= 0.0
    do k=1,nz
      u(k) = u0(k) ; v(k) = v0(k) ; T(k) = T0(k) ; Sal(k) = S0(k)
    enddo
  endif

  if (present(S2)) then
    S2(1) = 0.0 ; S2(nz+1) = 0.0
    if (ks > 1) &
      S2(ks) = ((u(ks)-u0(ks-1))**2 + (v(ks)-v0(ks-1))**2) * I_dz_int(ks)**2
    do K=ks+1,ke
      S2(K) = ((u(k)-u(k-1))**2 + (v(k)-v(k-1))**2) * I_dz_int(K)**2
    enddo
    if (ke<nz) &
      S2(ke+1) = ((u0(ke+1)-u(ke))**2 + (v0(ke+1)-v(ke))**2) * I_dz_int(ke+1)**2
  endif

  if (present(N2)) then
    N2(1) = 0.0 ; N2(nz+1) = 0.0
    if (ks > 1) &
      N2(ks) = max(0.0, I_dz_int(ks) * &
        (dbuoy_dT(ks) * (T0(ks-1)-T(ks)) + dbuoy_dS(ks) * (S0(ks-1)-Sal(ks))))
    do K=ks+1,ke
      N2(K) = max(0.0, I_dz_int(K) * &
        (dbuoy_dT(K) * (T(k-1)-T(k)) + dbuoy_dS(K) * (Sal(k-1)-Sal(k))))
    enddo
    if (ke<nz) &
      N2(ke+1) = max(0.0, I_dz_int(ke+1) * &
        (dbuoy_dT(ke+1) * (T(ke)-T0(ke+1)) + dbuoy_dS(ke+1) * (Sal(ke)-S0(ke+1))))
  endif

end subroutine calculate_projected_state

subroutine find_kappa_tke(N2, S2, kappa_in, Idz, dz_Int, I_L2_bdry, f2, &
                          nz, CS, K_Q, tke, kappa, kappa_src, local_src)
  integer,               intent(in)  :: nz
  real, dimension(nz+1), intent(in)  :: N2, S2, kappa_in, dz_Int, I_L2_bdry
  real, dimension(nz),   intent(in)  :: Idz
  real,                  intent(in)  :: f2
  type(Kappa_shear_CS),  pointer     :: CS
  real, dimension(nz+1), intent(inout) :: K_Q
  real, dimension(nz+1), intent(out) :: tke, kappa
  real, dimension(nz+1), optional, intent(out) :: kappa_src
  real, dimension(nz+1), optional, intent(out) :: local_src
!   This subroutine calculates new, consistent estimates of TKE and kappa.

! Arguments: N2 - The buoyancy frequency squared at interfaces, in s-2.
!  (in)      S2 - The squared shear at interfaces, in s-2.
!  (in)      kappa_in - The initial guess at the diffusivity, in m2 s-1.
!  (in)      Idz - The inverse grid spacing of layers, in m-1.
!  (in)      dz_Int - The thicknesses associated with interfaces, in m.
!  (in)      I_L2_bdry - The inverse of the squared distance to boundaries, m2.
!  (in)      f2 - The squared Coriolis parameter, in s-2.
!  (in)      nz - The number of layers to work on.
!  (in)      CS - A pointer to this module's control structure.
!  (inout)   K_Q - The shear-driven diapycnal diffusivity divided by the
!                  turbulent kinetic energy per unit mass at interfaces, in s.
!  (out)     tke - The turbulent kinetic energy per unit mass at interfaces,
!                  in units of m2 s-2.
!  (out)     kappa - The diapycnal diffusivity at interfaces, in m2 s-1.
!  (out,opt) kappa_src - The source term for kappa, in s-1.
!  (out,opt) local_src - The sum of all local sources for kappa, in s-1.

! UNCOMMENT THE FOLLOWING IF NOT CONTAINED IN Calculate_kappa_shear
  real, dimension(nz) :: &
    aQ, &       ! aQ is the coupling between adjacent interfaces in the TKE
                ! equations, in m s-1.
    dQdz        ! Half the partial derivative of TKE with depth, m s-2.
  real, dimension(nz+1) :: &
    dK, &         ! The change in kappa, in m2 s-1.
    dQ, &         ! The change in TKE, in m2 s-1.
    cQ, cK, &     ! cQ and cK are the upward influences in the tridiagonal and
                  ! hexadiagonal solvers for the TKE and kappa equations, ND.
    I_Ld2, &      ! 1/Ld^2, where Ld is the effective decay length scale
                  ! for kappa, in units of m-2.
    TKE_decay, &  ! The local TKE decay rate in s-1.
    k_src, &      ! The source term in the kappa equation, in s-1.
    dQmdK, &      ! With Newton's method the change in dQ(k-1) due to dK(k), s.
    dKdQ, &       ! With Newton's method the change in dK(k) due to dQ(k), s-1.
    e1            ! The fractional change in a layer TKE due to a change in the
                  ! TKE of the layer above when all the kappas below are 0.
                  ! e1 is nondimensional, and 0 < e1 < 1.
  real :: tke_src       ! The net source of TKE due to mixing against the shear
                        ! and stratification, in m2 s-3.  (For convenience,
                        ! a term involving the non-dissipation of q0 is also
                        ! included here.)
  real :: bQ, bK        ! The inverse of the pivot in the tridiagonal equations.
  real :: bd1           ! A term in the denominator of bQ or bK.
  real :: cQcomp, cKcomp ! 1 - cQ or 1 - cK in the tridiagonal equations.
  real :: c_s2          !   The coefficient for the decay of TKE due to
                        ! shear (i.e. proportional to |S|*tke), nondimensional.
  real :: c_n2          !   The coefficient for the decay of TKE due to
                        ! stratification (i.e. proportional to N*tke), nondim.
  real :: Ri_crit       !   The critical shear Richardson number for shear-
                        ! driven mixing. The theoretical value is 0.25.
  real :: q0            !   The background level of TKE, in m2 s-2.
  real :: Ilambda2      ! 1.0 / CS%lambda**2.
  real :: TKE_min       !   The minimum value of shear-driven TKE that can be
                        ! solved for, in m2 s-2.
  real :: kappa0        ! The background diapycnal diffusivity, in m2 s-1.
  real :: max_err       ! The maximum value of norm_err in a column, nondim.
  real :: kappa_trunc   ! Diffusivities smaller than this are rounded to 0, m2 s-1.

  real :: eden1, eden2, I_eden, ome  ! Variables used in calculating e1.
  real :: diffusive_src ! The diffusive source in the kappa equation, in m s-1.
  real :: chg_by_k0     ! The value of k_src that leads to an increase of
                        ! kappa_0 if only the diffusive term is a sink, in s-1.

  real :: kappa_mean    ! A mean value of kappa, in m2 s-1.
  real :: Newton_test   ! The value of relative error that will cause the next
                        ! iteration to use Newton's method.
  ! Temporary variables used in the Newton's method iterations.
  real :: decay_term, I_Q, kap_src, v1, v2

  real :: tol_err        ! The tolerance for max_err that determines when to
                         ! stop iterating.
  real :: Newton_err     ! The tolerance for max_err that determines when to
                         ! start using Newton's method.  Empirically, an initial
                         ! value of about 0.2 seems to be most efficient.
  real, parameter :: roundoff = 1.0e-16 ! A negligible fractional change in TKE.
                         ! This could be larger but performance gains are small.

  logical :: tke_noflux_bottom_BC = .false. ! Specify the boundary conditions
  logical :: tke_noflux_top_BC = .false.    ! that are applied to the TKE eqns.
  logical :: do_Newton    ! If .true., use Newton's method for the next iteration.
  logical :: abort_Newton ! If .true., an Newton's method has encountered a 0
                          ! pivot, and should not have been used.
  logical :: was_Newton   ! The value of do_Newton before checking convergence.
  logical :: within_tolerance ! If .true., all points are within tolerance to
                          ! enable this subroutine to return.
  integer :: ks_src, ke_src ! The range indices that have nonzero k_src.
  integer :: ks_kappa, ke_kappa, ke_tke   ! The ranges of k-indices that are or
  integer :: ks_kappa_prev, ke_kappa_prev ! were being worked on.
  integer :: itt, k, k2
#ifdef DEBUG
  integer :: max_debug_itt ; parameter(max_debug_itt=20)
  real :: K_err_lin, Q_err_lin
  real, dimension(nz+1) :: &
    kappa_prev, & ! The value of kappa at the start of the current iteration, in m2 s-1.
    TKE_prev   ! The value of TKE at the start of the current iteration, in m2 s-2.
  real, dimension(nz+1,1:max_debug_itt) :: &
    tke_it1, kappa_it1, kprev_it1, &  ! Various values from each iteration.
    dkappa_it1, K_Q_it1, d_dkappa_it1, dkappa_norm_it1
  real :: norm_err      ! The absolute change in kappa between iterations,
                        ! normalized by the value of kappa, nondim.
  real :: max_TKE_err, min_TKE_err, TKE_err(nz)  ! Various normalized TKE changes.
  integer :: it2
#endif

  c_N2 = CS%C_N**2 ; c_S2 = CS%C_S**2
  q0 = CS%TKE_bg ; kappa0 = CS%kappa_0 ; TKE_min = max(CS%TKE_bg,1.0E-20)
  Ri_crit = CS%Rino_crit
  Ilambda2 = 1.0 / CS%lambda**2
  kappa_trunc = 0.01*kappa0  ! ### CHANGE THIS HARD-WIRING LATER?
  do_Newton = .false. ; abort_Newton = .false.
  tol_err = CS%kappa_tol_err
  Newton_err = 0.2     ! This initial value may be automatically reduced later.

  ks_kappa = 2 ; ke_kappa = nz ; ks_kappa_prev = 2 ; ke_kappa_prev = nz

  ke_src = 0 ; ks_src = nz+1
  do K=2,nz
    if (N2(K) < Ri_crit * S2(K)) then ! Equivalent to Ri < Ri_crit.
!       Ri = N2(K) / S2(K)
!       k_src(K) = (2.0 * CS%Shearmix_rate * sqrt(S2(K))) * &
!                  ((Ri_crit - Ri) / (Ri_crit + CS%FRi_curvature*Ri))
      k_src(K) = (2.0 * CS%Shearmix_rate * sqrt(S2(K))) * &
                 ((Ri_crit*S2(K) - N2(K)) / (Ri_crit*S2(K) + CS%FRi_curvature*N2(K)))
      ke_src = K
      if (ks_src > k) ks_src = K
    else
      k_src(K) = 0.0
    endif
  enddo

  ! If there is no source anywhere, return kappa(K) = 0.
  if (ks_src > ke_src) then
    do K=1,nz+1
      kappa(K) = 0.0 ; K_Q(K) = 0.0 ; tke(K) = TKE_min
    enddo
    if (present(kappa_src)) then ; do K=1,nz+1 ; kappa_src(K) = 0.0 ; enddo ; endif
    if (present(local_src)) then ; do K=1,nz+1 ; local_src(K) = 0.0 ; enddo ; endif
    return
  endif

  do K=1,nz+1
    kappa(K) = kappa_in(K)
!     TKE_decay(K) = c_n*sqrt(N2(K)) + c_s*sqrt(S2(K)) ! The expression in JHL.
    TKE_decay(K) = sqrt(c_n2*N2(K) + c_s2*S2(K))
    if ((kappa(K) > 0.0) .and. (K_Q(K) > 0.0)) then
      TKE(K) = kappa(K) / K_Q(K)
    else
      TKE(K) = TKE_min
    endif
  enddo
  ! Apply boundary conditions to kappa.
  kappa(1) = 0.0 ; kappa(nz+1) = 0.0

  ! Calculate the term (e1) that allows changes in TKE to be calculated quickly
  ! below the deepest nonzero value of kappa.  If kappa = 0, below interface
  ! k-1, the final changes in TKE are related by dQ(K+1) = e1(K+1)*dQ(K).
  eden2 = kappa0 * Idz(nz)
  if (tke_noflux_bottom_BC) then
    eden1 = dz_Int(nz+1)*TKE_decay(nz+1)
    I_eden = 1.0 / (eden2 + eden1)
    e1(nz+1) = eden2 * I_eden ; ome = eden1 * I_eden
  else
    e1(nz+1) = 0.0 ; ome = 1.0
  endif
  do k=nz,2,-1
    eden1 = dz_Int(K)*TKE_decay(K) + ome * eden2
    eden2 = kappa0 * Idz(k-1)
    I_eden = 1.0 / (eden2 + eden1)
    e1(K) = eden2 * I_eden ; ome = eden1 * I_eden ! = 1-e1
  enddo
  e1(1) = 0.0


  ! Iterate here to convergence to within some tolerance of order tol_err.
  do itt=1,CS%max_RiNo_it

  ! ----------------------------------------------------
  ! Calculate TKE
  ! ----------------------------------------------------

#ifdef DEBUG
    do K=1,nz+1 ; kappa_prev(K) = kappa(K) ; TKE_prev(K) = TKE(K) ; enddo
#endif

    if (.not.do_Newton) then
      !   Use separate steps of the TKE and kappa equations, that are
      ! explicit in the nonlinear source terms, implicit in a linearized
      ! version of the nonlinear sink terms, and implicit in the linear
      ! terms.

      ke_tke = max(ke_kappa,ke_kappa_prev)+1
      ! aQ is the coupling between adjacent interfaces in m s-1.
      do k=1,min(ke_tke,nz)
        aQ(k) = (0.5*(kappa(K)+kappa(K+1))+kappa0) * Idz(k)
      enddo
      dQ(1) = -TKE(1)
      if (tke_noflux_top_BC) then
        tke_src = kappa0*S2(1) + q0 * TKE_decay(1) ! Uses that kappa(1) = 0
        bd1 = dz_Int(1) * TKE_decay(1)
        bQ = 1.0 / (bd1 +  aQ(1))
        tke(1) = bQ * (dz_Int(1)*tke_src)
        cQ(2) = aQ(1) * bQ ; cQcomp = bd1 * bQ ! = 1 - cQ
      else
        tke(1) = q0 ; cQ(2) = 0.0 ; cQcomp = 1.0
      endif
      do K=2,ke_tke-1
        dQ(K) = -TKE(K)
        tke_src = (kappa(K) + kappa0)*S2(K) + q0*TKE_decay(K)
        bd1 = dz_Int(K)*(TKE_decay(K) + N2(K)*K_Q(K)) + cQcomp*aQ(k-1)
        bQ = 1.0 / (bd1 + aQ(k))
        tke(K) = bQ * (dz_Int(K)*tke_src + aQ(k-1)*tke(K-1))
        cQ(K+1) = aQ(k) * bQ ; cQcomp = bd1 * bQ ! = 1 - cQ
      enddo
      if ((ke_tke == nz+1) .and. .not.(tke_noflux_bottom_BC)) then
        tke(nz+1) = TKE_min
        dQ(nz+1) = 0.0
      else
        k = ke_tke
        tke_src = kappa0*S2(K) + q0*TKE_decay(K) ! Uses that kappa(ke_tke) = 0
        if (K == nz+1) then
          dQ(K) = -TKE(K)
          bQ = 1.0 / (dz_Int(K)*TKE_decay(K) + cQcomp*aQ(k-1))
          tke(K) = max(TKE_min, bQ * (dz_Int(K)*tke_src + aQ(k-1)*tke(K-1)))
          dQ(K) = tke(K) + dQ(K)
        else
          bQ = 1.0 / ((dz_Int(K)*TKE_decay(K) + cQcomp*aQ(k-1)) + aQ(k))
          cQ(K+1) = aQ(k) * bQ
          ! Account for all changes deeper in the water column.
          dQ(K) = -TKE(K)
          tke(K) = max((bQ * (dz_Int(K)*tke_src + aQ(k-1)*tke(K-1)) + &
                        cQ(K+1)*(tke(K+1) - e1(K+1)*tke(K))) / &
                       (1.0 - cQ(K+1)*e1(K+1)), TKE_min)
          dQ(K) = tke(K) + dQ(K)

          ! Adjust TKE deeper in the water column in case ke_tke increases.
          ! This might not be strictly necessary?
          do K=ke_tke+1,nz+1
            dQ(K) = e1(K)*dQ(K-1)
            tke(K) = max(tke(K) + dQ(K), TKE_min)
            if (abs(dQ(K)) < roundoff*tke(K)) exit
          enddo
          do K2=K+1,nz
            if (dQ(K2) == 0.0) exit
            dQ(K2) = 0.0
          enddo
        endif
      endif
      do K=ke_tke-1,1,-1
        tke(K) = max(tke(K) + cQ(K+1)*tke(K+1), TKE_min)
        dQ(K) = tke(K) + dQ(K)
      enddo

  ! ----------------------------------------------------
  ! Calculate kappa, here defined at interfaces.
  ! ----------------------------------------------------

      ke_kappa_prev = ke_kappa ; ks_kappa_prev = ks_kappa

      dK(1) = 0.0 ! kappa takes boundary values of 0.
      cK(2) = 0.0 ; cKcomp = 1.0
      if (itt == 1) then ; dO K=2,nz
        I_Ld2(K) = (N2(K)*Ilambda2 + f2) / tke(K) + I_L2_bdry(K)
      enddo ; endif
      do K=2,nz
        dK(K) = -kappa(K)
        if (itt>1) &
          I_Ld2(K) = (N2(K)*Ilambda2 + f2) / tke(K) + I_L2_bdry(K)
        bd1 = dz_Int(K)*I_Ld2(K) + cKcomp*Idz(k-1)
        bK = 1.0 / (bd1 + Idz(k))

        kappa(K) = bK * (Idz(k-1)*kappa(K-1) + dz_Int(K) * k_src(K))
        cK(K+1) = Idz(k) * bK ; cKcomp = bd1 * bK ! = 1 - cK(K+1)

        ! Neglect values that are smaller than kappa_trunc.
        if (kappa(K) < cKcomp*kappa_trunc) then
          kappa(K) = 0.0
          if (K > ke_src) then ; ke_kappa = k-1 ; K_Q(K) = 0.0 ; exit ; endif
        elseif (kappa(K) < 2.0*cKcomp*kappa_trunc) then
          kappa(K) = 2.0 * (kappa(K) - cKcomp*kappa_trunc)
        endif
      enddo
      K_Q(ke_kappa) = kappa(ke_kappa) / tke(ke_kappa)
      dK(ke_kappa) = dK(ke_kappa) + kappa(ke_kappa)
      do K=ke_kappa+2,ke_kappa_prev
        dK(K) = -kappa(K) ; kappa(K) = 0.0 ; K_Q(K) = 0.0
      enddo
      do K=ke_kappa-1,2,-1
        kappa(K) = kappa(K) + cK(K+1)*kappa(K+1)
        ! Neglect values that are smaller than kappa_trunc.
        if (kappa(K) <= kappa_trunc) then
          kappa(K) = 0.0
          if (K < ks_src) then ; ks_kappa = k+1 ; K_Q(K) = 0.0 ; exit ; endif
        elseif (kappa(K) < 2.0*kappa_trunc) then
          kappa(K) = 2.0 * (kappa(K) - kappa_trunc)
        endif

        dK(K) = dK(K) + kappa(K)
        K_Q(K) = kappa(K) / tke(K)
      enddo
      do K=ks_kappa_prev,ks_kappa-2 ; kappa(K) = 0.0 ; K_Q(K) = 0.0 ; enddo

    else ! do_Newton is .true.
!   Once the solutions are close enough, use a Newton's method solver of the
!  whole system to accelerate convergence.
      ks_kappa_prev = ks_kappa ; ke_kappa_prev = ke_kappa ; ke_kappa = nz
      ks_kappa = 2
      dK(1) = 0.0 ; cK(2) = 0.0 ; cKcomp = 1.0 ; dKdQ(1) = 0.0
      aQ(1) = (0.5*(kappa(1)+kappa(2))+kappa0) * Idz(1)
      dQdz(1) = 0.5*(TKE(1) - TKE(2))*Idz(1)
      if (tke_noflux_top_BC) then
        tke_src = dz_Int(1) * (kappa0*S2(1) - (TKE(1) - q0)*TKE_decay(1)) - &
                  aQ(1) * (TKE(1) - TKE(2))

        bQ = 1.0 / (aQ(1) + dz_Int(1)*TKE_decay(1))
        cQ(2) = aQ(1) * bQ
        cQcomp = (dz_Int(1)*TKE_decay(1)) * bQ ! = 1 - cQ(2)
        dQmdK(2) = -dQdz(1) * bQ
        dQ(1) = bQ * tke_src
      else
        dQ(1) = 0.0 ; cQ(2) = 0.0 ; cQcomp = 1.0 ; dQmdK(2) = 0.0
      endif
      do K=2,nz
        I_Q = 1.0 / TKE(K)
        I_Ld2(K) = (N2(K)*Ilambda2 + f2) * I_Q + I_L2_bdry(K)

        kap_src = dz_Int(K) * (k_src(K) - I_Ld2(K)*kappa(K)) + &
                  Idz(k-1)*(kappa(K-1)-kappa(K)) - Idz(k)*(kappa(K)-kappa(K+1))

        ! Ensure that the pivot is always positive, and that 0 <= cK <= 1.
        ! Otherwise do not use Newton's method.
        decay_term = -Idz(k-1)*dQmdK(K)*dKdQ(K-1) + dz_Int(K)*I_Ld2(K)
        if (decay_term < 0.0) then ; abort_Newton = .true. ; exit ; endif
        bK = 1.0 / (Idz(k) + Idz(k-1)*cKcomp + decay_term)

        cK(K+1) = bK * Idz(k)
        cKcomp = bK * (Idz(k-1)*cKcomp + decay_term) ! = 1-cK(K+1)
        dKdQ(K) = bK * (Idz(k-1)*dKdQ(K-1)*cQ(K) + &
                        (N2(K)*Ilambda2 + f2)*I_Q**2*kappa(K))
        dK(K) = bK * (kap_src + Idz(k-1)*dK(K-1) + Idz(k-1)*dKdQ(K-1)*dQ(K-1))

        ! Truncate away negligibly small values of kappa.
        if (dK(K) <= cKcomp*(kappa_trunc - kappa(K))) then
          dK(K) = -cKcomp*kappa(K)
!         if (K > ke_src) then ; ke_kappa = k-1 ; K_Q(K) = 0.0 ; exit ; endif
        elseif (dK(K) < cKcomp*(2.0*kappa_trunc - kappa(K))) then
          dK(K) = 2.0 * dK(K) - cKcomp*(2.0*kappa_trunc - kappa(K))
        endif

        ! Solve for dQ(K)...
        aQ(k) = (0.5*(kappa(K)+kappa(K+1))+kappa0) * Idz(k)
        dQdz(k) = 0.5*(TKE(K) - TKE(K+1))*Idz(k)
        tke_src = dz_Int(K) * ((kappa(k) + kappa0)*S2(k) - kappa(k)*N2(k) - &
                               (TKE(k) - q0)*TKE_decay(k)) - &
                  (aQ(k) * (TKE(K) - TKE(K+1)) - aQ(k-1) * (TKE(K-1) - TKE(K)))
        v1 = aQ(k-1) + dQdz(k-1)*dKdQ(K-1)
        v2 = (v1*dQmdK(K) + dQdz(k-1)*cK(K)) + &
             ((dQdz(k-1) - dQdz(k)) + dz_Int(K)*(S2(K) - N2(K)))

        ! Ensure that the pivot is always positive, and that 0 <= cQ <= 1.
        ! Otherwise do not use Newton's method.
        decay_term = dz_Int(K)*TKE_decay(K) - dQdz(k-1)*dKdQ(K-1)*cQ(K) - v2*dKdQ(K)
        if (decay_term < 0.0) then ; abort_Newton = .true. ; exit ; endif
        bQ = 1.0 / (aQ(k) + (cQcomp*aQ(k-1) + decay_term))

        cQ(K+1) = aQ(k) * bQ
        cQcomp = (cQcomp*aQ(k-1) + decay_term) * bQ
        dQmdK(K+1) = (v2 * cK(K+1) - dQdz(k)) * bQ

        ! Ensure that TKE+dQ will not drop below 0.5*TKE.
        dQ(K) = max(bQ * ((v1 * dQ(K-1) + dQdz(k-1)*dK(k-1)) + &
                          (v2 * dK(K) + tke_src)), cQcomp*(-0.5*TKE(K)))

        ! Check whether the next layer will be affected by any nonzero kappas.
        if ((itt > 1) .and. (K > ke_src) .and. (dK(K) == 0.0) .and. &
            ((kappa(K) + kappa(K+1)) == 0.0)) then
        ! Could also do  .and. (bQ*abs(tke_src) < roundoff*TKE(K)) then
          ke_kappa = k-1 ; exit
        endif
      enddo
      if ((ke_kappa == nz) .and. (.not. abort_Newton)) then
        dK(nz+1) = 0.0 ; dKdQ(nz+1) = 0.0
        if (tke_noflux_bottom_BC) then
          K = nz+1
          tke_src = dz_Int(K) * (kappa0*S2(K) - (TKE(K) - q0)*TKE_decay(K)) + &
                    aQ(k-1) * (TKE(K-1) - TKE(K))

          v1 = aQ(k-1) + dQdz(k-1)*dKdQ(K-1)
          decay_term = max(0.0, dz_Int(K)*TKE_decay(K) - dQdz(k-1)*dKdQ(K-1)*cQ(K))
          if (decay_term < 0.0) then
            abort_Newton = .true.
          else
            bQ = 1.0 / (aQ(k) + (cQcomp*aQ(k-1) + decay_term))
          ! Ensure that TKE+dQ will not drop below 0.5*TKE.
            dQ(K) = max(bQ * ((v1 * dQ(K-1) + dQdz(k-1)*dK(K-1)) + tke_src), &
                        -0.5*TKE(K))
            TKE(K) = max(TKE(K) + dQ(K), TKE_min)
          endif
        else
          dQ(nz+1) = 0.0
        endif
      elseif (.not. abort_Newton) then
        ! Alter the first-guess determination of dQ(K).
        dQ(ke_kappa+1) = dQ(ke_kappa+1) / (1.0 - cQ(ke_kappa+2)*e1(ke_kappa+2))
        TKE(ke_kappa+1) = max(TKE(ke_kappa+1) + dQ(ke_kappa+1), TKE_min)
        do k=ke_kappa+2,nz+1
#ifdef DEBUG
          if (K < nz+1) then
          ! Ignore this source?
            aQ(k) = (0.5*(kappa(K)+kappa(K+1))+kappa0) * Idz(k)
            tke_src = (dz_Int(K) * (kappa0*S2(K) - (TKE(K)-q0)*TKE_decay(K)) - &
                       (aQ(k) * (TKE(K) - TKE(K+1)) - &
                        aQ(k-1) * (TKE(K-1) - TKE(K))) ) / &
                       (aQ(k) + (aQ(k-1) + dz_Int(K)*TKE_decay(K)))
          endif
#endif
          dK(K) = 0.0
        ! Ensure that TKE+dQ will not drop below 0.5*TKE.
          dQ(K) = max(e1(K)*dQ(K-1),-0.5*TKE(K))
          TKE(K) = max(TKE(K) + dQ(K), TKE_min)
          if (abs(dQ(K)) < roundoff*TKE(K)) exit
        enddo
#ifdef DEBUG
        do K2=K+1,ke_kappa_prev+1 ; dQ(K2) = 0.0 ; dK(K2) = 0.0 ; enddo
        do K=K2,nz+1 ; if (dQ(K) == 0.0) exit ; dQ(K) = 0.0 ; dK(K) = 0.0 ; enddo
#endif
      endif
      if (.not. abort_Newton) then
        do K=ke_kappa,2,-1
          ! Ensure that TKE+dQ will not drop below 0.5*TKE.
          dQ(K) = max(dQ(K) + (cQ(K+1)*dQ(K+1) + dQmdK(K+1) * dK(K+1)), &
                      -0.5*TKE(K))
          TKE(K) = max(TKE(K) + dQ(K), TKE_min)
          dK(K) = dK(K) + (cK(K+1)*dK(K+1) + dKdQ(K) * dQ(K))
          ! Truncate away negligibly small values of kappa.
          if (dK(K) <= kappa_trunc - kappa(K)) then
            dK(K) = -kappa(K)
            kappa(K) = 0.0
            if ((K < ks_src) .and. (K+1 > ks_kappa)) ks_kappa = K+1
          elseif (dK(K) < 2.0*kappa_trunc - kappa(K)) then
            dK(K) =  2.0*dK(K) - (2.0*kappa_trunc - kappa(K))
            kappa(K) = max(kappa(K) + dK(K), 0.0) ! The max is for paranoia.
            if (K<=ks_kappa) ks_kappa = 2
          else
            kappa(K) = kappa(K) + dK(K)
            if (K<=ks_kappa) ks_kappa = 2
          endif
        enddo
        dQ(1) = max(dQ(1) + cQ(2)*dQ(2) + dQmdK(2) * dK(2), TKE_min - TKE(1))
        TKE(1) = max(TKE(1) + dQ(1), TKE_min)
        dK(1) = 0.0
      endif

#ifdef DEBUG
  ! Check these solutions for consistency.
      do K=2,nz
        ! In these equations, K_err_lin and Q_err_lin should be at round-off levels
        ! compared with the dominant terms, perhaps, dz_Int*I_Ld2*kappa and
        ! dz_Int*TKE_decay*TKE.  The exception is where, either 1) the decay term has been
        ! been increased to ensure a positive pivot, or 2) negative TKEs have been
        ! truncated, or 3) small or negative kappas have been rounded toward 0.
        I_Q = 1.0 / TKE(K)
        I_Ld2(K) = (N2(K)*Ilambda2 + f2) * I_Q + I_L2_bdry(K)

        kap_src = dz_Int(K) * (k_src(K) - I_Ld2(K)*kappa_prev(K)) + &
                  (Idz(k-1)*(kappa_prev(k-1)-kappa_prev(k)) - &
                   Idz(k)*(kappa_prev(k)-kappa_prev(k+1)))
        K_err_lin =  -Idz(k-1)*(dK(K-1)-dK(K)) + Idz(k)*(dK(K)-dK(K+1)) + &
                     dz_Int(K)*I_Ld2(K)*dK(K) - kap_src - &
                     (N2(K)*Ilambda2 + f2)*I_Q**2*kappa_prev(K) * dQ(K)

        tke_src = dz_Int(K) * ((kappa_prev(K) + kappa0)*S2(K) - &
                     kappa_prev(K)*N2(K) - (TKE_prev(K) - q0)*TKE_decay(K)) - &
                  (aQ(k) * (TKE_prev(K) - TKE_prev(K+1)) - &
                   aQ(k-1) * (TKE_prev(K-1) - TKE_prev(K)))
        Q_err_lin = (aQ(k-1) * (dQ(K-1)-dQ(K)) - aQ(k) * (dQ(k)-dQ(k+1))) - &
          0.5*(TKE_prev(K)-TKE_prev(K+1))*Idz(k)  * (dK(K) + dK(K+1)) - &
          0.5*(TKE_prev(K)-TKE_prev(K-1))*Idz(k-1)* (dK(K-1) + dK(K)) + &
          dz_Int(K) * (dK(K) * (S2(K) - N2(K)) - dQ(K)*TKE_decay(K)) + tke_src
      enddo
#endif
    endif  ! End of the Newton's method solver.

    ! Test kappa for convergence...
#ifdef DEBUG
    max_err = 0.0 ; max_TKE_err = 0.0 ; min_TKE_err = 0.0
    do K=min(ks_kappa,ks_kappa_prev),max(ke_kappa,ke_kappa_prev)
      norm_err = abs(kappa(K) - kappa_prev(K)) / &
                    (kappa0 + 0.5*(kappa(K) + kappa_prev(K)))
      if (max_err < norm_err) max_err = norm_err

      TKE_err(K) = dQ(K) / (tke(K) - 0.5*dQ(K))
      if (TKE_err(K) > max_TKE_err) max_TKE_err = TKE_err(K)
      if (TKE_err(K) < min_TKE_err) min_TKE_err = TKE_err(K)
    enddo
    if (do_Newton) then
      if (max(max_err,max_TKE_err,-min_TKE_err) >= 2.0*Newton_err) then
        do_Newton = .false. ; abort_Newton = .true.
      endif
    else
      if (max(max_err,max_TKE_err,-min_TKE_err) < Newton_err) do_Newton = .true.
    endif
    within_tolerance = (max_err < tol_err)
#else
 !   max_err = 0.0
    if ((tol_err < Newton_err) .and. (.not.abort_Newton)) then
      !   A lower tolerance is used to switch to Newton's method than to
      ! switch back.
      Newton_test = Newton_err ; if (do_Newton) Newton_test = 2.0*Newton_err
      was_Newton = do_Newton
      within_tolerance = .true. ; do_Newton = .true.
      do K=min(ks_kappa,ks_kappa_prev),max(ke_kappa,ke_kappa_prev)
        kappa_mean = kappa0 + (kappa(K) - 0.5*dK(K))
        if (abs(dK(K)) > Newton_test * kappa_mean) then
          if (do_Newton) abort_Newton = .true.
          within_tolerance = .false. ; do_Newton = .false. ; exit
        elseif (abs(dK(K)) > tol_err * kappa_mean) then
          within_tolerance = .false. ; if (.not.do_Newton) exit
        endif
        if (abs(dQ(K)) > Newton_test*(tke(K) - 0.5*dQ(K))) then
          if (do_Newton) abort_Newton = .true.
          do_Newton = .false. ; if (.not.within_tolerance) exit
        endif
      enddo

    else  ! Newton's method will not be used again, so no need to check.
      within_tolerance = .true.
      do K=min(ks_kappa,ks_kappa_prev),max(ke_kappa,ke_kappa_prev)
        if (abs(dK(K)) > tol_err * (kappa0 + (kappa(K) - 0.5*dK(K)))) then
          within_tolerance = .false. ;  exit
        endif
      enddo
    endif
#endif

    if (abort_Newton) then
      do_Newton = .false. ; abort_Newton = .false.
      ! We went to Newton too quickly last time, so restrict the tolerance.
      Newton_err = 0.5*Newton_err
      ke_kappa_prev = nz
      do K=2,nz ; K_Q(K) = kappa(K) / max(TKE(K), TKE_min) ; enddo
    endif

#ifdef DEBUG
    if (itt <= max_debug_itt) then
      do K=1,nz+1
        kprev_it1(K,itt)=kappa_prev(K)
        kappa_it1(K,itt)=kappa(K) ; tke_it1(K,itt) = tke(K)
        dkappa_it1(K,itt) = kappa(K) - kappa_prev(K)
        dkappa_norm_it1(K,itt) = (kappa(K) - kappa_prev(K)) / &
            (kappa0 + 0.5*(kappa(K) + kappa_prev(K)))
        K_Q_it1(K,itt) = kappa(K) / max(TKE(K),TKE_min)
        d_dkappa_it1(K,itt) = 0.0
        if (itt > 1) then ; if (abs(dkappa_it1(K,itt-1)) > 1e-20) &
            d_dkappa_it1(K,itt) = dkappa_it1(K,itt) / dkappa_it1(K,itt-1)
        endif
      enddo
    endif
#endif

    if (within_tolerance) exit

  enddo

#ifdef DEBUG
  do it2=itt+1,max_debug_itt ; do K=1,nz+1
    kprev_it1(K,it2) = 0.0 ; kappa_it1(K,it2) = 0.0 ; tke_it1(K,it2) = 0.0
    dkappa_it1(K,it2) = 0.0 ; K_Q_it1(K,it2) = 0.0 ; d_dkappa_it1(K,it2) = 0.0
  enddo ; enddo
#endif

  if (do_Newton) then  ! K_Q needs to be calculated.
    do K=1,ks_kappa-1 ;  K_Q(K) = 0.0 ; enddo
    do K=ks_kappa,ke_kappa ; K_Q(K) = kappa(K) / TKE(K) ; enddo
    do K=ke_kappa+1,nz+1 ; K_Q(K) = 0.0 ; enddo
  endif

  if (present(local_src)) then
    local_src(1) = 0.0 ; local_src(nz+1) = 0.0
    do K=2,nz
      diffusive_src = Idz(k-1)*(kappa(K-1)-kappa(K)) + &
                  Idz(k)*(kappa(K+1)-kappa(K))
      chg_by_k0 = kappa0 * ((Idz(k-1)+Idz(k)) / dz_Int(K) + I_Ld2(K))
      if (diffusive_src <= 0.0) then
        local_src(K) = k_src(K) + chg_by_k0
      else
        local_src(K) = (k_src(K) + chg_by_k0) + diffusive_src / dz_Int(K)
      endif
    enddo
  endif
  if (present(kappa_src)) then
    kappa_src(1) = 0.0 ; kappa_src(nz+1) = 0.0
    do K=2,nz
      kappa_src(K) = k_src(K)
    enddo
  endif

end subroutine find_kappa_tke


logical function kappa_shear_init(Time, G, GV, param_file, diag, CS)
  type(time_type),         intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(verticalGrid_type), intent(in)    :: GV
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ctrl), target, intent(inout) :: diag
  type(Kappa_shear_CS),    pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (returns) kappa_shear_init - True if module is to be used, False otherwise
  logical :: merge_mixedlayer
! This include declares and sets the variable "version".
#include "version_variable.h"
  real :: KD_normal ! The KD of the main model, read here only as a parameter
                    ! for setting the default of KD_SMOOTH
  if (associated(CS)) then
    call MOM_error(WARNING, "kappa_shear_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  !   The Jackson-Hallberg-Legg shear mixing parameterization uses the following
  ! 6 nondimensional coefficients.  That paper gives 3 best fit parameter sets.
  !    Ri_Crit  Rate    FRi_Curv  K_buoy  TKE_N  TKE_Shear
  ! p1: 0.25    0.089    -0.97     0.82    0.24    0.14
  ! p2: 0.30    0.085    -0.94     0.86    0.26    0.13
  ! p3: 0.35    0.088    -0.89     0.81    0.28    0.12
  !   Future research will reveal how these should be modified to take
  ! subgridscale inhomogeneity into account.

! Set default, read and log parameters
  call log_version(param_file, mod, version, &
    "Parameterization of shear-driven turbulence following Jackson, Hallberg and Legg, JPO 2008")
  call get_param(param_file, mod, "USE_JACKSON_PARAM", kappa_shear_init, &
                 "If true, use the Jackson-Hallberg-Legg (JPO 2008) \n"//&
                 "shear mixing parameterization.", default=.false.)
  call get_param(param_file, mod, "RINO_CRIT", CS%RiNo_crit, &
                 "The critical Richardson number for shear mixing.", &
                 units="nondim", default=0.25)
  call get_param(param_file, mod, "SHEARMIX_RATE", CS%Shearmix_rate, &
                 "A nondimensional rate scale for shear-driven entrainment.\n"//&
                 "Jackson et al find values in the range of 0.085-0.089.", &
                 units="nondim", default=0.089)
  call get_param(param_file, mod, "MAX_RINO_IT", CS%max_RiNo_it, &
                 "The maximum number of iterations that may be used to \n"//&
                 "estimate the Richardson number driven mixing.", &
                 units="nondim", default=50)
  call get_param(param_file, mod, "KD", KD_normal, default=1.0e-7, do_not_log=.true.)
  call get_param(param_file, mod, "KD_KAPPA_SHEAR_0", CS%kappa_0, &
                 "The background diffusivity that is used to smooth the \n"//&
                 "density and shear profiles before solving for the \n"//&
                 "diffusivities. Defaults to value of KD.", units="m2 s-1", default=KD_normal)
  call get_param(param_file, mod, "FRI_CURVATURE", CS%FRi_curvature, &
                 "The nondimensional curvature of the function of the \n"//&
                 "Richardson number in the kappa source term in the \n"//&
                 "Jackson et al. scheme.", units="nondim", default=-0.97)
  call get_param(param_file, mod, "TKE_N_DECAY_CONST", CS%C_N, &
                 "The coefficient for the decay of TKE due to \n"//&
                 "stratification (i.e. proportional to N*tke). \n"//&
                 "The values found by Jackson et al. are 0.24-0.28.", &
                 units="nondim", default=0.24)
!  call get_param(param_file, mod, "LAYER_KAPPA_STAGGER", CS%layer_stagger, &
!                 default=.false.)
  call get_param(param_file, mod, "TKE_SHEAR_DECAY_CONST", CS%C_S, &
                 "The coefficient for the decay of TKE due to shear (i.e. \n"//&
                 "proportional to |S|*tke). The values found by Jackson \n"//&
                 "et al. are 0.14-0.12.", units="nondim", default=0.14)
  call get_param(param_file, mod, "KAPPA_BUOY_SCALE_COEF", CS%lambda, &
                 "The coefficient for the buoyancy length scale in the \n"//&
                 "kappa equation.  The values found by Jackson et al. are \n"//&
                 "in the range of 0.81-0.86.", units="nondim", default=0.82)
  call get_param(param_file, mod, "KAPPA_N_OVER_S_SCALE_COEF2", CS%lambda2_N_S, &
                 "The square of the ratio of the coefficients of the \n"//&
                 "buoyancy and shear scales in the diffusivity equation, \n"//&
                 "Set this to 0 (the default) to eliminate the shear scale. \n"//&
                 "This is only used if USE_JACKSON_PARAM is true.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mod, "KAPPA_SHEAR_TOL_ERR", CS%kappa_tol_err, &
                 "The fractional error in kappa that is tolerated. \n"//&
                 "Iteration stops when changes between subsequent \n"//&
                 "iterations are smaller than this everywhere in a \n"//&
                 "column.  The peak diffusivities usually converge most \n"//&
                 "rapidly, and have much smaller errors than this.", &
                 units="nondim", default=0.1)
  call get_param(param_file, mod, "TKE_BACKGROUND", CS%TKE_bg, &
                 "A background level of TKE used in the first iteration \n"//&
                 "of the kappa equation.  TKE_BACKGROUND could be 0.", &
                 units="m2 s-2", default=0.0)
  call get_param(param_file, mod, "KAPPA_SHEAR_ELIM_MASSLESS", CS%eliminate_massless, &
                 "If true, massless layers are merged with neighboring \n"//&
                 "massive layers in this calculation.  The default is \n"//&
                 "true and I can think of no good reason why it should \n"//&
                 "be false. This is only used if USE_JACKSON_PARAM is true.", &
                 default=.true.)
  call get_param(param_file, mod, "MAX_KAPPA_SHEAR_IT", CS%max_KS_it, &
                 "The maximum number of iterations that may be used to \n"//&
                 "estimate the time-averaged diffusivity.", units="nondim", &
                 default=13)
  call get_param(param_file, mod, "PRANDTL_TURB", CS%Prandtl_turb, &
                 "The turbulent Prandtl number applied to shear \n"//&
                 "instability.", units="nondim", default=1.0, do_not_log=.true.)
  call get_param(param_file, mod, "DEBUG_KAPPA_SHEAR", CS%debug, &
                 "If true, write debugging data for the kappa-shear code. \n"//&
                 "Caution: this option is _very_ verbose and should only \n"//&
                 "be used in single-column mode!", default=.false.)

!    id_clock_KQ = cpu_clock_id('Ocean KS kappa_shear',grain=CLOCK_ROUTINE)
!    id_clock_avg = cpu_clock_id('Ocean KS avg',grain=CLOCK_ROUTINE)
!    id_clock_project = cpu_clock_id('Ocean KS project',grain=CLOCK_ROUTINE)
!    id_clock_setup = cpu_clock_id('Ocean KS setup',grain=CLOCK_ROUTINE)

  CS%nkml = 1
  if (GV%nkml>0) then
    call get_param(param_file, mod, "KAPPA_SHEAR_MERGE_ML",merge_mixedlayer, &
                 "If true, combine the mixed layers together before \n"//&
                 "solving the kappa-shear equations.", default=.true.)
    if (merge_mixedlayer) CS%nkml = GV%nkml
  endif

! Forego remainder of initialization if not using this scheme
  if (.not. kappa_shear_init) return

  CS%diag => diag

  CS%id_Kd_shear = register_diag_field('ocean_model','Kd_shear',diag%axesTi,Time, &
      'Shear-driven Diapycnal Diffusivity', 'meter2 second-1')
  CS%id_TKE = register_diag_field('ocean_model','TKE_shear',diag%axesTi,Time, &
      'Shear-driven Turbulent Kinetic Energy', 'meter2 second-2')
#ifdef ADD_DIAGNOSTICS
  CS%id_ILd2 = register_diag_field('ocean_model','ILd2_shear',diag%axesTi,Time, &
      'Inverse kappa decay scale at interfaces', 'meter-2')
  CS%id_dz_Int = register_diag_field('ocean_model','dz_Int_shear',diag%axesTi,Time, &
      'Finite volume thickness of interfaces', 'meter')
#endif

end function kappa_shear_init

logical function kappa_shear_is_used(param_file)
! Reads the parameter "USE_JACKSON_PARAM" and returns state.
!   This function allows other modules to know whether this parameterization will
! be used without needing to duplicate the log entry.
  type(param_file_type), intent(in) :: param_file
  call get_param(param_file, mod, "USE_JACKSON_PARAM", kappa_shear_is_used, &
                 default=.false., do_not_log = .true.)
end function kappa_shear_is_used

end module MOM_kappa_shear
