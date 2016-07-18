!> Calculations of isoneutral slopes and stratification.
module MOM_isopycnal_slopes

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : int_specific_vol_dp, calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>

public calc_isoneutral_slopes

contains

subroutine calc_isoneutral_slopes(G, GV, h, e, tv, dt_kappa_smooth, &
                                  slope_x, slope_y, N2_u, N2_v, halo)
  type(ocean_grid_type),                       intent(in)    :: G
  type(verticalGrid_type),                     intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),    intent(in)    :: h
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1),  intent(in)    :: e
  type(thermo_var_ptrs),                       intent(in)    :: tv
  real,                                        intent(in)    :: dt_kappa_smooth
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: slope_x
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1), intent(inout) :: slope_y
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: N2_u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1), intent(inout) :: N2_v
  optional                                                     :: N2_u, N2_v
  integer, optional,                             intent(in)    :: halo
  ! Local variables
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)) :: &
    T, &          ! The temperature (or density) in C, with the values in
                  ! in massless layers filled vertically by diffusion.
    S, &          ! The filled salinity, in PSU, with the values in
                  ! in massless layers filled vertically by diffusion.
    Rho           ! Density itself, when a nonlinear equation of state is
                  ! not in use.
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)+1) :: &
    pres          ! The pressure at an interface, in Pa.
  real, dimension(SZIB_(G)) :: &
    drho_dT_u, &  ! The derivatives of density with temperature and
    drho_dS_u     ! salinity at u points, in kg m-3 K-1 and kg m-3 psu-1.
  real, dimension(SZI_(G)) :: &
    drho_dT_v, &  ! The derivatives of density with temperature and
    drho_dS_v     ! salinity at v points, in kg m-3 K-1 and kg m-3 psu-1.
  real, dimension(SZIB_(G)) :: &
    T_u, S_u, &   ! Temperature, salinity, and pressure on the interface at
    pres_u        ! the u-point in the horizontal.
  real, dimension(SZI_(G)) :: &
    T_v, S_v, &   ! Temperature, salinity, and pressure on the interface at
    pres_v        ! the v-point in the horizontal.
  real :: drdiA, drdiB  ! Along layer zonal- and meridional- potential density
  real :: drdjA, drdjB  ! gradients in the layers above (A) and below(B) the
                        ! interface times the grid spacing, in kg m-3.
  real :: drdkL, drdkR  ! Vertical density differences across an interface,
                        ! in kg m-3.
  real :: hg2A, hg2B, hg2L, hg2R
  real :: haA, haB, haL, haR
  real :: dzaL, dzaR
  real :: wtA, wtB, wtL, wtR
  real :: drdx, drdy, drdz  ! Zonal, meridional, and vertical density gradients,
                            ! in units of kg m-4.
  real :: Slope         ! The slope of density surfaces, calculated in a way
                        ! that is always between -1 and 1.
  real :: mag_grad2     ! The squared magnitude of the 3-d density gradient, in kg2 m-8.
  real :: slope2_Ratio  ! The ratio of the slope squared to slope_max squared.
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected, in H.
  real :: h_neglect2    ! h_neglect^2, in H2.
  real :: dz_neglect    ! A thickness in m that is so small it is usually lost
                        ! in roundoff and can be neglected, in m.
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  real :: G_Rho0, N2, dzN2,  H_x(SZIB_(G)), H_y(SZI_(G))

  logical :: present_N2_u, present_N2_v
  integer :: is, ie, js, je, nz, IsdB
  integer :: i, j, k

  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  else
    is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  endif
  nz = G%ke ; IsdB = G%IsdB

  h_neglect = GV%H_subroundoff ; h_neglect2 = h_neglect**2
  dz_neglect = GV%H_subroundoff*GV%H_to_m

  use_EOS = associated(tv%eqn_of_state)

  present_N2_u = PRESENT(N2_u)
  present_N2_v = PRESENT(N2_v)
  G_Rho0 = GV%g_Earth / GV%Rho0
  if (present_N2_u) then
    do j=js,je ; do I=is-1,ie
      N2_u(I,j,1) = 0.
      N2_u(I,j,nz+1) = 0.
    enddo ; enddo
  endif
  if (present_N2_v) then
    do J=js-1,je ; do i=is,ie
      N2_v(i,J,1) = 0.
      N2_v(i,J,nz+1) = 0.
    enddo ; enddo
  endif

  if (use_EOS) then
    if (present(halo)) then
      call vert_fill_TS(h, tv%T, tv%S, dt_kappa_smooth, 1.0, T, S, G, GV, halo+1)
    else
      call vert_fill_TS(h, tv%T, tv%S, dt_kappa_smooth, 1.0, T, S, G, GV, 1)
    endif
  endif

  ! Find the maximum and minimum permitted streamfunction.
!$OMP parallel default(none) shared(is,ie,js,je,pres,GV,h,nz)
!$OMP do
  do j=js-1,je+1 ; do i=is-1,ie+1
    pres(i,j,1) = 0.0  ! ### This should be atmospheric pressure.
    pres(i,j,2) = pres(i,j,1) + GV%H_to_Pa*h(i,j,1)
  enddo ; enddo
!$OMP do
  do j=js-1,je+1
    do k=2,nz ; do i=is-1,ie+1
      pres(i,j,K+1) = pres(i,j,K) + GV%H_to_Pa*h(i,j,k)
    enddo ; enddo
  enddo
!$OMP end parallel

!$OMP parallel do default(none) shared(nz,is,ie,js,je,use_EOS,G,GV,pres,T,S, &
!$OMP                                  IsdB,tv,h,h_neglect,e,dz_neglect,  &
!$OMP                                  h_neglect2,present_N2_u,G_Rho0,N2_u,slope_x) &
!$OMP                          private(drdiA,drdiB,drdkL,drdkR,pres_u,T_u,S_u,      &
!$OMP                                  drho_dT_u,drho_dS_u,hg2A,hg2B,hg2L,hg2R,haA, &
!$OMP                                  haB,haL,haR,dzaL,dzaR,wtA,wtB,wtL,wtR,drdz,  &
!$OMP                                  drdx,mag_grad2,Slope,slope2_Ratio)
  do j=js,je ; do K=nz,2,-1
    if (.not.(use_EOS)) then
      drdiA = 0.0 ; drdiB = 0.0
      drdkL = GV%Rlay(k)-GV%Rlay(k-1) ; drdkR = GV%Rlay(k)-GV%Rlay(k-1)
    endif

    ! Calculate the zonal isopycnal slope.
    if (use_EOS) then
      do I=is-1,ie
        pres_u(I) = 0.5*(pres(i,j,K) + pres(i+1,j,K))
        T_u(I) = 0.25*((T(i,j,k) + T(i+1,j,k)) + (T(i,j,k-1) + T(i+1,j,k-1)))
        S_u(I) = 0.25*((S(i,j,k) + S(i+1,j,k)) + (S(i,j,k-1) + S(i+1,j,k-1)))
      enddo
      call calculate_density_derivs(T_u, S_u, pres_u, drho_dT_u, &
                   drho_dS_u, (is-IsdB+1)-1, ie-is+2, tv%eqn_of_state)
    endif

    do I=is-1,ie
      if (use_EOS) then
        ! Estimate the horizontal density gradients along layers.
        drdiA = drho_dT_u(I) * (T(i+1,j,k-1)-T(i,j,k-1)) + &
                drho_dS_u(I) * (S(i+1,j,k-1)-S(i,j,k-1))
        drdiB = drho_dT_u(I) * (T(i+1,j,k)-T(i,j,k)) + &
                drho_dS_u(I) * (S(i+1,j,k)-S(i,j,k))

        ! Estimate the vertical density gradients times the grid spacing.
        drdkL = (drho_dT_u(I) * (T(i,j,k)-T(i,j,k-1)) + &
                 drho_dS_u(I) * (S(i,j,k)-S(i,j,k-1)))
        drdkR = (drho_dT_u(I) * (T(i+1,j,k)-T(i+1,j,k-1)) + &
                 drho_dS_u(I) * (S(i+1,j,k)-S(i+1,j,k-1)))
      endif


      if (use_EOS) then
        hg2A = h(i,j,k-1)*h(i+1,j,k-1) + h_neglect2
        hg2B = h(i,j,k)*h(i+1,j,k) + h_neglect2
        hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
        hg2R = h(i+1,j,k-1)*h(i+1,j,k) + h_neglect2
        haA = 0.5*(h(i,j,k-1) + h(i+1,j,k-1))
        haB = 0.5*(h(i,j,k) + h(i+1,j,k)) + h_neglect
        haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
        haR = 0.5*(h(i+1,j,k-1) + h(i+1,j,k)) + h_neglect
        if (GV%Boussinesq) then
          dzaL = haL * GV%H_to_m ; dzaR = haR * GV%H_to_m
        else
          dzaL = 0.5*(e(i,j,K-1) - e(i,j,K+1)) + dz_neglect
          dzaR = 0.5*(e(i+1,j,K-1) - e(i+1,j,K+1)) + dz_neglect
        endif
        ! Use the harmonic mean thicknesses to weight the horizontal gradients.
        ! These unnormalized weights have been rearranged to minimize divisions.
        wtA = hg2A*haB ; wtB = hg2B*haA
        wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

        drdz = (wtL * drdkL + wtR * drdkR) / (dzaL*wtL + dzaR*wtR)
        ! The expression for drdz above is mathematically equivalent to:
        !   drdz = ((hg2L/haL) * drdkL/dzaL + (hg2R/haR) * drdkR/dzaR) / &
        !          ((hg2L/haL) + (hg2R/haR))
        ! This is the gradient of density along geopotentials.
        drdx = ((wtA * drdiA + wtB * drdiB) / (wtA + wtB) - &
                drdz * (e(i,j,K)-e(i+1,j,K))) * G%IdxCu(I,j)

        ! This estimate of slope is accurate for small slopes, but bounded
        ! to be between -1 and 1.
        mag_grad2 = drdx**2 + drdz**2
        if (mag_grad2 > 0.0) then
          slope_x(I,j,K) = drdx / sqrt(mag_grad2)
        else ! Just in case mag_grad2 = 0 ever.
          slope_x(I,j,K) = 0.0
        endif

        if (present_N2_u) N2_u(I,j,k) = G_Rho0 * drdz ! Square of Brunt-Vaisala frequency (s-2)

      else ! With .not.use_EOS, the layers are constant density.
        slope_x(I,j,K) = ((e(i,j,K)-e(i+1,j,K))*G%IdxCu(I,j)) * GV%m_to_H
      endif

    enddo ! I
  enddo ; enddo ! end of j-loop

    ! Calculate the meridional isopycnal slope.
!$OMP parallel do default(none) shared(nz,is,ie,js,je,use_EOS,G,GV,pres,T,S, &
!$OMP                                  IsdB,tv,h,h_neglect,e,dz_neglect,  &
!$OMP                                  h_neglect2,present_N2_v,G_Rho0,N2_v,slope_y) &
!$OMP                          private(drdjA,drdjB,drdkL,drdkR,pres_v,T_v,S_v,      &
!$OMP                                  drho_dT_v,drho_dS_v,hg2A,hg2B,hg2L,hg2R,haA, &
!$OMP                                  haB,haL,haR,dzaL,dzaR,wtA,wtB,wtL,wtR,drdz,  &
!$OMP                                  drdy,mag_grad2,Slope,slope2_Ratio)
  do j=js-1,je ; do K=nz,2,-1
    if (.not.(use_EOS)) then
      drdjA = 0.0 ; drdjB = 0.0
      drdkL = GV%Rlay(k)-GV%Rlay(k-1) ; drdkR = GV%Rlay(k)-GV%Rlay(k-1)
    endif

    if (use_EOS) then
      do i=is,ie
        pres_v(i) = 0.5*(pres(i,j,K) + pres(i,j+1,K))
        T_v(i) = 0.25*((T(i,j,k) + T(i,j+1,k)) + (T(i,j,k-1) + T(i,j+1,k-1)))
        S_v(i) = 0.25*((S(i,j,k) + S(i,j+1,k)) + (S(i,j,k-1) + S(i,j+1,k-1)))
      enddo
      call calculate_density_derivs(T_v, S_v, pres_v, drho_dT_v, &
                   drho_dS_v, is, ie-is+1, tv%eqn_of_state)
    endif
    do i=is,ie
      if (use_EOS) then
        ! Estimate the horizontal density gradients along layers.
        drdjA = drho_dT_v(i) * (T(i,j+1,k-1)-T(i,j,k-1)) + &
                drho_dS_v(i) * (S(i,j+1,k-1)-S(i,j,k-1))
        drdjB = drho_dT_v(i) * (T(i,j+1,k)-T(i,j,k)) + &
                drho_dS_v(i) * (S(i,j+1,k)-S(i,j,k))

        ! Estimate the vertical density gradients times the grid spacing.
        drdkL = (drho_dT_v(i) * (T(i,j,k)-T(i,j,k-1)) + &
                 drho_dS_v(i) * (S(i,j,k)-S(i,j,k-1)))
        drdkR = (drho_dT_v(i) * (T(i,j+1,k)-T(i,j+1,k-1)) + &
                 drho_dS_v(i) * (S(i,j+1,k)-S(i,j+1,k-1)))
      endif

      if (use_EOS) then
        hg2A = h(i,j,k-1)*h(i,j+1,k-1) + h_neglect2
        hg2B = h(i,j,k)*h(i,j+1,k) + h_neglect2
        hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
        hg2R = h(i,j+1,k-1)*h(i,j+1,k) + h_neglect2
        haA = 0.5*(h(i,j,k-1) + h(i,j+1,k-1)) + h_neglect
        haB = 0.5*(h(i,j,k) + h(i,j+1,k)) + h_neglect
        haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
        haR = 0.5*(h(i,j+1,k-1) + h(i,j+1,k)) + h_neglect
        if (GV%Boussinesq) then
          dzaL = haL * GV%H_to_m ; dzaR = haR * GV%H_to_m
        else
          dzaL = 0.5*(e(i,j,K-1) - e(i,j,K+1)) + dz_neglect
          dzaR = 0.5*(e(i,j+1,K-1) - e(i,j+1,K+1)) + dz_neglect
        endif
        ! Use the harmonic mean thicknesses to weight the horizontal gradients.
        ! These unnormalized weights have been rearranged to minimize divisions.
        wtA = hg2A*haB ; wtB = hg2B*haA
        wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

        drdz = (wtL * drdkL + wtR * drdkR) / (dzaL*wtL + dzaR*wtR)
        ! The expression for drdz above is mathematically equivalent to:
        !   drdz = ((hg2L/haL) * drdkL/dzaL + (hg2R/haR) * drdkR/dzaR) / &
        !          ((hg2L/haL) + (hg2R/haR))
        ! This is the gradient of density along geopotentials.
        drdy = ((wtA * drdjA + wtB * drdjB) / (wtA + wtB) - &
                drdz * (e(i,j,K)-e(i,j+1,K))) * G%IdyCv(i,J)

        ! This estimate of slope is accurate for small slopes, but bounded
        ! to be between -1 and 1.
        mag_grad2 = drdy**2 + drdz**2
        if (mag_grad2 > 0.0) then
          slope_y(i,J,K) = drdy / sqrt(mag_grad2)
        else ! Just in case mag_grad2 = 0 ever.
          slope_y(i,J,K) = 0.0
        endif

        if (present_N2_v) N2_v(i,J,k) = G_Rho0 * drdz ! Square of Brunt-Vaisala frequency (s-2)

      else ! With .not.use_EOS, the layers are constant density.
        slope_y(i,J,K) = ((e(i,j,K)-e(i,j+1,K))*G%IdyCv(i,J)) * GV%m_to_H
      endif

    enddo ! i
  enddo ; enddo ! end of j-loop

end subroutine calc_isoneutral_slopes

subroutine vert_fill_TS(h, T_in, S_in, kappa, dt, T_f, S_f, G, GV, halo_here)
  type(ocean_grid_type),                    intent(in)    :: G
  type(verticalGrid_type),                  intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: T_in
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: S_in
  real,                                     intent(in)    :: kappa
  real,                                     intent(in)    :: dt
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out)   :: T_f
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out)   :: S_f
  integer,                        optional, intent(in)    :: halo_here
!    This subroutine fills massless layers with sensible values of two
!*  tracer arrays (nominally temperature and salinity) by diffusing
!*  vertically with a (small?) constant diffusivity.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (in)      T_in - The input temperature, in K.
!  (in)      S_in - The input salinity, in psu.
!  (in)      kappa - The diapycnal diffusivity, in m2 s-1.
!  (in)      dt - Time increment in s.
!  (out)     T_f - The filled temperature, in K.
!  (out)     S_f - The filled salinity, in psu.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in,opt)  halo_here - the number of halo points to work on, 0 by default.

  real :: ent(SZI_(G),SZK_(G)+1)   ! The diffusive entrainment (kappa*dt)/dz
                                   ! between layers in a timestep in m or kg m-2.
  real :: b1(SZI_(G)), d1(SZI_(G)) ! b1, c1, and d1 are variables used by the
  real :: c1(SZI_(G),SZK_(G))      ! tridiagonal solver.
  real :: kap_dt_x2                ! The product of 2*kappa*dt, converted to
                                   ! the same units as h, in m2 or kg2 m-4.
  real :: h_neglect                ! A negligible thickness, in m or kg m-2, to
                                   ! allow for zero thicknesses.
  integer :: i, j, k, is, ie, js, je, nz, halo

  halo=0 ; if (present(halo_here)) halo = max(halo_here,0)

  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  nz = G%ke

  kap_dt_x2 = (2.0*kappa*dt)*GV%m_to_H**2
  h_neglect = GV%H_subroundoff

  if (kap_dt_x2 <= 0.0) then
!$OMP parallel do default(none) shared(is,ie,js,je,nz,T_f,T_in,S_f,S_in)
    do k=1,nz ; do j=js,je ; do i=is,ie
      T_f(i,j,k) = T_in(i,j,k) ; S_f(i,j,k) = S_in(i,j,k)
    enddo ; enddo ; enddo
  else
!$OMP parallel do default(none) private(ent,b1,d1,c1)   &
!$OMP                            shared(is,ie,js,je,nz,kap_dt_x2,h,h_neglect,T_f,S_f,T_in,S_in)
    do j=js,je
      do i=is,ie
        ent(i,2) = kap_dt_x2 / ((h(i,j,1)+h(i,j,2)) + h_neglect)
        b1(i) = 1.0 / (h(i,j,1)+ent(i,2))
        d1(i) = b1(i) * h(i,j,1)
        T_f(i,j,1) = (b1(i)*h(i,j,1))*T_in(i,j,1)
        S_f(i,j,1) = (b1(i)*h(i,j,1))*S_in(i,j,1)
      enddo
      do k=2,nz-1 ; do i=is,ie
        ent(i,K+1) = kap_dt_x2 / ((h(i,j,k)+h(i,j,k+1)) + h_neglect)
        c1(i,k) = ent(i,K) * b1(i)
        b1(i) = 1.0 / ((h(i,j,k) + d1(i)*ent(i,K)) + ent(i,K+1))
        d1(i) = b1(i) * (h(i,j,k) + d1(i)*ent(i,K))
        T_f(i,j,k) = b1(i) * (h(i,j,k)*T_in(i,j,k) + ent(i,K)*T_f(i,j,k-1))
        S_f(i,j,k) = b1(i) * (h(i,j,k)*S_in(i,j,k) + ent(i,K)*S_f(i,j,k-1))
      enddo ; enddo
      do i=is,ie
        c1(i,nz) = ent(i,nz) * b1(i)
        b1(i) = 1.0 / (h(i,j,nz) + d1(i)*ent(i,nz) + h_neglect)
        T_f(i,j,nz) = b1(i) * (h(i,j,nz)*T_in(i,j,nz) + ent(i,nz)*T_f(i,j,nz-1))
        S_f(i,j,nz) = b1(i) * (h(i,j,nz)*S_in(i,j,nz) + ent(i,nz)*S_f(i,j,nz-1))
      enddo
      do k=nz-1,1,-1 ; do i=is,ie
        T_f(i,j,k) = T_f(i,j,k) + c1(i,k+1)*T_f(i,j,k+1)
        S_f(i,j,k) = S_f(i,j,k) + c1(i,k+1)*S_f(i,j,k+1)
      enddo ; enddo
    enddo
  endif

end subroutine vert_fill_TS

end module MOM_isopycnal_slopes
