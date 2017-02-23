!> Accelerations due to the Coriolis force and momentum advection
module MOM_CoriolisAdv

! This file is part of MOM6. See LICENSE.md for the license.

!> \author Robert Hallberg, April 1994 - June 2002

use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type
use MOM_string_functions, only : uppercase
use MOM_variables,     only : accel_diag_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

public CorAdCalc, CoriolisAdv_init, CoriolisAdv_end

#include <MOM_memory.h>

!> Control structure for mom_coriolisadv
type, public :: CoriolisAdv_CS ; private
  integer :: Coriolis_Scheme !< Selects the discretization for the Coriolis terms.
                             !! Valid values are:
                             !! - SADOURNY75_ENERGY - Sadourny, 1975
                             !! - ARAKAWA_HSU90     - Arakawa & Hsu, 1990, Energy & non-div. Enstrophy
                             !! - ROBUST_ENSTRO     - Pseudo-enstrophy scheme
                             !! - SADOURNY75_ENSTRO - Sadourny, JAS 1975, Enstrophy
                             !! - ARAKAWA_LAMB81    - Arakawa & Lamb, MWR 1981, Energy & Enstrophy
                             !! - ARAKAWA_LAMB_BLEND - A blend of Arakawa & Lamb with Arakawa & Hsu and Sadourny energy.
                             !! The default, SADOURNY75_ENERGY, is the safest choice then the
                             !! deformation radius is poorly resolved.
  integer :: KE_Scheme       !< KE_SCHEME selects the discretization for
                             !! the kinetic energy. Valid values are:
                             !!  KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV
  integer :: PV_Adv_Scheme   !< PV_ADV_SCHEME selects the discretization for PV advection
                             !! Valid values are:
                             !! - PV_ADV_CENTERED - centered (aka Sadourny, 75)
                             !! - PV_ADV_UPWIND1  - upwind, first order
  real    :: F_eff_max_blend !< The factor by which the maximum effective Coriolis
                             !! acceleration from any point can be increased when
                             !! blending different discretizations with the
                             !! ARAKAWA_LAMB_BLEND Coriolis scheme.  This must be
                             !! greater than 2.0, and is 4.0 by default.
  real    :: wt_lin_blend    !< A weighting value beyond which the blending between
                             !! Sadourny and Arakawa & Hsu goes linearly to 0.
                             !! This must be between 1 and 1e-15, often 1/8.
  logical :: no_slip         !< If true, no slip boundary conditions are used.
                             !! Otherwise free slip boundary conditions are assumed.
                             !! The implementation of the free slip boundary
                             !! conditions on a C-grid is much cleaner than the
                             !! no slip boundary conditions. The use of free slip
                             !! b.c.s is strongly encouraged. The no slip b.c.s
                             !! are not implemented with the biharmonic viscosity.
  logical :: bound_Coriolis  !< If true, the Coriolis terms at u points are
                             !! bounded by the four estimates of (f+rv)v from the
                             !! four neighboring v points, and similarly at v
                             !! points.  This option would have no effect on the
                             !! SADOURNY75_ENERGY scheme if it were possible to
                             !! use centered difference thickness fluxes.
  logical :: Coriolis_En_Dis !< If CORIOLIS_EN_DIS is defined, two estimates of
                             !! the thickness fluxes are used to estimate the
                             !! Coriolis term, and the one that dissipates energy
                             !! relative to the other one is used.  This is only
                             !! available at present if Coriolis scheme is
                             !! SADOURNY75_ENERGY.
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the timing of diagnostic output.
  integer :: id_rv = -1, id_PV = -1, id_gKEu = -1, id_gKEv = -1
  integer :: id_rvxu = -1, id_rvxv = -1
end type CoriolisAdv_CS

! Enumeration values for Coriolis_Scheme
integer, parameter :: SADOURNY75_ENERGY = 1
integer, parameter :: ARAKAWA_HSU90     = 2
integer, parameter :: ROBUST_ENSTRO     = 3
integer, parameter :: SADOURNY75_ENSTRO = 4
integer, parameter :: ARAKAWA_LAMB81    = 5
integer, parameter :: AL_BLEND          = 6
character*(20), parameter :: SADOURNY75_ENERGY_STRING = "SADOURNY75_ENERGY"
character*(20), parameter :: ARAKAWA_HSU_STRING = "ARAKAWA_HSU90"
character*(20), parameter :: ROBUST_ENSTRO_STRING = "ROBUST_ENSTRO"
character*(20), parameter :: SADOURNY75_ENSTRO_STRING = "SADOURNY75_ENSTRO"
character*(20), parameter :: ARAKAWA_LAMB_STRING = "ARAKAWA_LAMB81"
character*(20), parameter :: AL_BLEND_STRING = "ARAKAWA_LAMB_BLEND"
! Enumeration values for KE_Scheme
integer, parameter :: KE_ARAKAWA        = 10
integer, parameter :: KE_SIMPLE_GUDONOV = 11
integer, parameter :: KE_GUDONOV        = 12
character*(20), parameter :: KE_ARAKAWA_STRING = "KE_ARAKAWA"
character*(20), parameter :: KE_SIMPLE_GUDONOV_STRING = "KE_SIMPLE_GUDONOV"
character*(20), parameter :: KE_GUDONOV_STRING = "KE_GUDONOV"
! Enumeration values for PV_Adv_Scheme
integer, parameter :: PV_ADV_CENTERED   = 21
integer, parameter :: PV_ADV_UPWIND1    = 22
character*(20), parameter :: PV_ADV_CENTERED_STRING = "PV_ADV_CENTERED"
character*(20), parameter :: PV_ADV_UPWIND1_STRING = "PV_ADV_UPWIND1"

contains

!> Calculates the Coriolis and momentum advection contributions to the acceleration.
subroutine CorAdCalc(u, v, h, uh, vh, CAu, CAv, OBC, AD, G, GV, CS)
  type(ocean_grid_type),                     intent(in)    :: G !< Ocen grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u !< Zonal velocity (m/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v !< Meridional velocity (m/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h !< Layer thickness (m or kg/m2)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: uh !< Zonal transport u*h*dy (m3/s or kg/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: vh !< Meridional transport v*h*dx (m3/s or kg/s)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out)   :: CAu !< Zonal acceleration due to Coriolis
                                                                  !! and momentum advection, in m/s2.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out)   :: CAv !< Meridional acceleration due to Coriolis
                                                                  !! and momentum advection, in m/s2.
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary control structure
  type(accel_diag_ptrs),                     intent(inout) :: AD !< Storage for acceleration diagnostics
  type(CoriolisAdv_CS),                      pointer       :: CS !< Control structure for MOM_CoriolisAdv
  ! Local variables

  real, dimension(SZIB_(G),SZJB_(G)) :: &
    q, &        ! Layer potential vorticity, in m-1 s-1.
    Ih_q, &     ! The inverse of thickness interpolated to q pointes, in
                ! units of m-1 or m2 kg-1.
    Area_q      ! The sum of the ocean areas at the 4 adjacent thickness
                ! points, in m2.

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    a, b, c, d  ! a, b, c, & d are combinations of the potential vorticities
                ! surrounding an h grid point.  At small scales, a = q/4,
                ! b = q/4, etc.  All are in units of m-1 s-1 or m2 kg-1 s-1,
                ! and use the indexing of the corresponding u point.

  real, dimension(SZI_(G),SZJ_(G)) :: &
    Area_h, &   ! The ocean area at h points, in m2.  Area_h is used to find the
                ! average thickness in the denominator of q.  0 for land points.
    KE          ! Kinetic energy per unit mass, KE = (u^2 + v^2)/2, in m2 s-2.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    KEx         ! The zonal gradient of Kinetic energy per unit mass,
                ! KEx = d/dx KE, in m s-2.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    KEy         ! The meridonal gradient of Kinetic energy per unit mass,
                ! KEy = d/dy KE, in m s-2.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    uh_min, uh_max, &   ! The smallest and largest estimates of the volume
    vh_min, vh_max, &   ! fluxes through the faces (i.e. u*h*dy & v*h*dx),
                        ! in m3 s-1 or kg s-1.
    ep_u, ep_v  ! Additional pseudo-Coriolis terms in the Arakawa and Lamb
                ! discretization, in m-1 s-1 or m2 kg-1 s-1.
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    dvdx,dudy, &! Contributions to the circulation around q-points (m2 s-1)
    abs_vort, & ! Absolute vorticity at q-points, in s-1.
    q2, &       ! Relative vorticity over thickness.
    max_fvq, &  ! The maximum or minimum of the
    min_fvq, &  ! adjacent values of (-u) or v times
    max_fuq, &  ! the absolute vorticity, in m s-2.
    min_fuq     ! All are defined at q points.
  real, dimension(SZIB_(G),SZJB_(G),SZK_(G)) :: &
    PV, &       ! A diagnostic array of the potential vorticities, in m-1 s-1.
    RV          ! A diagnostic array of the relative vorticities, in s-1.
  real :: fv1, fv2, fu1, fu2   ! (f+rv)*v or (f+rv)*u in m s-2.
  real :: max_fv, max_fu       ! The maximum or minimum of the neighbor-
  real :: min_fv, min_fu       ! max(min)_fu(v)q, in m s-2.

  real, parameter :: C1_12=1.0/12.0 ! C1_12 = 1/12
  real, parameter :: C1_24=1.0/24.0 ! C1_24 = 1/24
  real :: absolute_vorticity     ! Absolute vorticity, in s-1.
  real :: relative_vorticity     ! Relative vorticity, in s-1.
  real :: Ih                     ! Inverse of thickness, m-1 or m2 kg-1.
  real :: max_Ihq, min_Ihq       ! The maximum and minimum of the nearby Ihq.
  real :: hArea_q                ! The sum of area times thickness of the cells
                                 ! surrounding a q point, in m3 or kg.
  real :: h_neglect              ! A thickness that is so small it is usually
                                 ! lost in roundoff and can be neglected, in m.
  real :: temp1, temp2           ! Temporary variables, in m2 s-2.
  real, parameter :: eps_vel=1.0e-10 ! A tiny, positive velocity, in m s-1.

  real :: uhc, vhc               ! Centered estimates of uh and vh in m3 s-1 or kg s-1.
  real :: uhm, vhm               ! The input estimates of uh and vh in m3 s-1 or kg s-1.
  real :: c1, c2, c3, slope      ! Nondimensional parameters for the Coriolis
                                 ! limiter scheme.

  real :: Fe_m2         ! Nondimensional temporary variables asssociated with
  real :: rat_lin       ! the ARAKAWA_LAMB_BLEND scheme.
  real :: rat_m1        ! The ratio of the maximum neighboring inverse thickness
                        ! to the minimum inverse thickness minus 1. rat_m1 >= 0.
  real :: AL_wt         ! The relative weight of the Arakawa & Lamb scheme to the
                        ! Arakawa & Hsu scheme, nondimensional between 0 and 1.
  real :: Sad_wt        ! The relative weight of the Sadourny energy scheme to
                        ! the other two with the ARAKAWA_LAMB_BLEND scheme,
                        ! nondimensional between 0 and 1.

  real :: Heff1, Heff2  ! Temporary effective H at U or V points in m or kg m-2.
  real :: Heff3, Heff4  ! Temporary effective H at U or V points in m or kg m-2.
  real :: h_tiny        ! A very small thickness, in m or kg m-2.
  real :: UHeff, VHeff  ! More temporary variables, in m3 s-1 or kg s-1.
  real :: QUHeff,QVHeff ! More temporary variables, in m3 s-2 or kg s-2.
  integer :: i, j, k, n, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz

! To work, the following fields must be set outside of the usual
! is to ie range before this subroutine is called:
!  v[is-1,ie+1,ie+2], u[is-1,ie+1], vh[ie+1], uh[is-1], and
!  h[is-1,ie+1,ie+2].
! In the y-direction, the following fields must be set:
!  v[js-1,je+1], u[js-1,je+1,je+2], vh[js-1], uh[je+1], and
!  h[js-1,je+1,je+2].

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "MOM_CoriolisAdv: Module must be initialized before it is used.")
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = G%ke
  h_neglect = GV%H_subroundoff
  h_tiny = GV%Angstrom  ! Perhaps this should be set to h_neglect instead.

!$OMP parallel default(none) shared(u,v,h,uh,vh,CAu,CAv,G,CS,AD,Area_h,Area_q,nz,RV,PV, &
!$OMP                               is,ie,js,je,Isq,Ieq,Jsq,Jeq,h_neglect,h_tiny,OBC)   &
!$OMP                       private(relative_vorticity,absolute_vorticity,Ih,hArea_q,q, &
!$OMP                               abs_vort,Ih_q,fv1,fv2,fu1,fu2,max_fvq,max_fuq,      &
!$OMP                               min_fvq,min_fuq,q2,a,b,c,d,ep_u,ep_v,Fe_m2,rat_lin, &
!$OMP                               min_Ihq,max_Ihq,rat_m1,AL_wt,Sad_wt,c1,c2,c3,slope, &
!$OMP                               uhc,uhm,uh_min,uh_max,vhc,vhm,vh_min,vh_max,KE,KEx, &
!$OMP                               KEy,temp1,temp2,Heff1,Heff2,Heff3,Heff4,VHeff,      &
!$OMP                               QVHeff,max_fv,min_fv,UHeff,QUHeff,max_fu,min_fu)
!$OMP do
  do j=Jsq-1,Jeq+2 ; do I=Isq-1,Ieq+2
    Area_h(i,j) = G%mask2dT(i,j) * G%areaT(i,j)
  enddo ; enddo
!$OMP do
  do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
    Area_q(i,j) = (Area_h(i,j) + Area_h(i+1,j+1)) + &
                  (Area_h(i+1,j) + Area_h(i,j+1))
  enddo ; enddo

!$OMP do
  do k=1,nz
    ! Here the second order accurate layer potential vorticities, q,
    ! are calculated.  hq is  second order accurate in space.  Relative
    ! vorticity is second order accurate everywhere with free slip b.c.s,
    ! but only first order accurate at boundaries with no slip b.c.s.
    ! First calculate the contributions to the circulation around the q-point.
    do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
      dvdx(I,J) = v(i+1,J,k)*G%dyCv(i+1,J) - v(i,J,k)*G%dyCv(i,J)
      dudy(I,J) = u(I,j+1,k)*G%dxCu(I,j+1) - u(I,j,k)*G%dxCu(I,j)
    enddo ; enddo
    ! Adjust circulation components to relative vorticity on open boundaries.
    if (associated(OBC)) then ; if (OBC%zero_vorticity .or. OBC%freeslip_vorticity) then
      do n=1,OBC%number_of_segments
        if (OBC%segment(n)%is_N_or_S) then
          J = OBC%segment(n)%HI%JsdB
          do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
            if (OBC%zero_vorticity) then
              dvdx(I,J) = 0.
              dudy(I,J) = 0.
            elseif (OBC%freeslip_vorticity) then
              dudy(I,J) = 0.
            endif
          enddo
        elseif (OBC%segment(n)%is_E_or_W) then
          I = OBC%segment(n)%HI%IsdB
          do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
            if (OBC%zero_vorticity) then
              dvdx(I,J) = 0.
              dudy(I,J) = 0.
            elseif (OBC%freeslip_vorticity) then
              dvdx(I,J) = 0.
            endif
          enddo
        endif
      enddo
    endif ; endif
    do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
      if (CS%no_slip ) then
        relative_vorticity = (2.0-G%mask2dBu(I,J)) * (dvdx(I,J) - dudy(I,J)) * &
            (G%IdxBu(I,J) * G%IdyBu(I,J)) ! ### Using G%IareaBu(I,J) changes answers.
      else
        relative_vorticity = G%mask2dBu(I,J) * (dvdx(I,J) - dudy(I,J)) * &
            (G%IdxBu(I,J) * G%IdyBu(I,J)) ! ### Using G%IareaBu(I,J) changes answers.
      endif
      absolute_vorticity = G%CoriolisBu(I,J) + relative_vorticity
      Ih = 0.0
      if (Area_q(i,j) > 0.0) then
        hArea_q = ((Area_h(i,j) * h(i,j,k) + Area_h(i+1,j+1) * h(i+1,j+1,k)) + &
                   (Area_h(i,j+1) * h(i,j+1,k) + Area_h(i+1,j) * h(i+1,j,k)))
        Ih = Area_q(i,j) / (hArea_q + h_neglect*Area_q(i,j))
      endif
      q(I,J) = absolute_vorticity * Ih
      abs_vort(I,J) = absolute_vorticity
      Ih_q(I,J) = Ih

      if (CS%bound_Coriolis) then
        fv1 = absolute_vorticity*v(i+1,J,k)
        fv2 = absolute_vorticity*v(i,J,k)
        fu1 = -absolute_vorticity*u(I,j+1,k)
        fu2 = -absolute_vorticity*u(I,j,k)
        if (fv1 > fv2) then
          max_fvq(I,J) = fv1 ; min_fvq(I,J) = fv2
        else
          max_fvq(I,J) = fv2 ; min_fvq(I,J) = fv1
        endif
        if (fu1 > fu2) then
          max_fuq(I,J) = fu1 ; min_fuq(I,J) = fu2
        else
          max_fuq(I,J) = fu2 ; min_fuq(I,J) = fu1
        endif
      endif

      if (CS%id_rv > 0) RV(I,J,k) = relative_vorticity
      if (CS%id_PV > 0) PV(I,J,k) = q(I,J)
      if (ASSOCIATED(AD%rv_x_v) .or. ASSOCIATED(AD%rv_x_u)) &
        q2(I,J) = relative_vorticity * Ih
    enddo ; enddo

    !   a, b, c, and d are combinations of neighboring potential
    ! vorticities which form the Arakawa and Hsu vorticity advection
    ! scheme.  All are defined at u grid points.

    if (CS%Coriolis_Scheme == ARAKAWA_HSU90) then
      do j=Jsq,Jeq+1
        do I=is-1,Ieq
          a(I,j) = (q(I,J) + (q(I+1,J) + q(I,J-1))) * C1_12
          d(I,j) = ((q(I,J) + q(I+1,J-1)) + q(I,J-1)) * C1_12
        enddo
        do I=Isq,Ieq
          b(I,j) = (q(I,J) + (q(I-1,J) + q(I,J-1))) * C1_12
          c(I,j) = ((q(I,J) + q(I-1,J-1)) + q(I,J-1)) * C1_12
        enddo
      enddo
    elseif (CS%Coriolis_Scheme == ARAKAWA_LAMB81) then
      do j=Jsq,Jeq+1 ; do I=Isq,Ieq+1
        a(I-1,j) = (2.0*(q(I,J) + q(I-1,J-1)) + (q(I-1,J) + q(I,J-1))) * C1_24
        d(I-1,j) = ((q(I,j) + q(I-1,J-1)) + 2.0*(q(I-1,J) + q(I,J-1))) * C1_24
        b(I,j) =   ((q(I,J) + q(I-1,J-1)) + 2.0*(q(I-1,J) + q(I,J-1))) * C1_24
        c(I,j) =   (2.0*(q(I,J) + q(I-1,J-1)) + (q(I-1,J) + q(I,J-1))) * C1_24
        ep_u(i,j) = ((q(I,J) - q(I-1,J-1)) + (q(I-1,J) - q(I,J-1))) * C1_24
        ep_v(i,j) = (-(q(I,J) - q(I-1,J-1)) + (q(I-1,J) - q(I,J-1))) * C1_24
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == AL_BLEND) then
      Fe_m2 = CS%F_eff_max_blend - 2.0
      rat_lin = 1.5 * Fe_m2 / max(CS%wt_lin_blend, 1.0e-16)

      ! This allows the code to always give Sadourny Energy
      if (CS%F_eff_max_blend <= 2.0) then ; Fe_m2 = -1. ; rat_lin = -1.0 ; endif

      do j=Jsq,Jeq+1 ; do I=Isq,Ieq+1
        min_Ihq = MIN(Ih_q(I-1,J-1), Ih_q(I,J-1), Ih_q(I-1,J), Ih_q(I,J))
        max_Ihq = MAX(Ih_q(I-1,J-1), Ih_q(I,J-1), Ih_q(I-1,J), Ih_q(I,J))
        rat_m1 = 1.0e15
        if (max_Ihq < 1.0e15*min_Ihq) rat_m1 = max_Ihq / min_Ihq - 1.0
        ! The weights used here are designed to keep the effective Coriolis
        ! acceleration from any one point on its neighbors within a factor
        ! of F_eff_max.  The minimum permitted value is 2 (the factor for
        ! Sadourny's energy conserving scheme).

        ! Determine the relative weights of Arakawa & Lamb vs. Arakawa and Hsu.
        if (rat_m1 <= Fe_m2) then ; AL_wt = 1.0
        elseif (rat_m1 < 1.5*Fe_m2) then ; AL_wt = 3.0*Fe_m2 / rat_m1 - 2.0
        else ; AL_wt = 0.0 ; endif

        ! Determine the relative weights of Sadourny Energy vs. the other two.
        if (rat_m1 <= 1.5*Fe_m2) then ; Sad_wt = 0.0
        elseif (rat_m1 <= rat_lin) then
          Sad_wt = 1.0 - (1.5*Fe_m2) / rat_m1
        elseif (rat_m1 < 2.0*rat_lin) then
          Sad_wt = 1.0 - (CS%wt_lin_blend / rat_lin) * (rat_m1 - 2.0*rat_lin)
        else ; Sad_wt = 1.0 ; endif

        a(I-1,j) = Sad_wt * 0.25 * q(I-1,J) + (1.0 - Sad_wt) * &
                   ( ((2.0-AL_wt)* q(I-1,J) + AL_wt*q(I,J-1)) + &
                      2.0 * (q(I,J) + q(I-1,J-1)) ) * C1_24
        d(I-1,j) = Sad_wt * 0.25 * q(I-1,J-1) + (1.0 - Sad_wt) * &
                   ( ((2.0-AL_wt)* q(I-1,J-1) + AL_wt*q(I,J)) + &
                      2.0 * (q(I-1,J) + q(I,J-1)) ) * C1_24
        b(I,j) =   Sad_wt * 0.25 * q(I,J) + (1.0 - Sad_wt) * &
                   ( ((2.0-AL_wt)* q(I,J) + AL_wt*q(I-1,J-1)) + &
                      2.0 * (q(I-1,J) + q(I,J-1)) ) * C1_24
        c(I,j) =   Sad_wt * 0.25 * q(I,J-1) + (1.0 - Sad_wt) * &
                   ( ((2.0-AL_wt)* q(I,J-1) + AL_wt*q(I-1,J)) + &
                      2.0 * (q(I,J) + q(I-1,J-1)) ) * C1_24
        ep_u(i,j) = AL_wt  * ((q(I,J) - q(I-1,J-1)) + (q(I-1,J) - q(I,J-1))) * C1_24
        ep_v(i,j) = AL_wt * (-(q(I,J) - q(I-1,J-1)) + (q(I-1,J) - q(I,J-1))) * C1_24
      enddo ; enddo
    endif

    if (CS%Coriolis_En_Dis) then
    !  c1 = 1.0-1.5*RANGE ; c2 = 1.0-RANGE ; c3 = 2.0 ; slope = 0.5
      c1 = 1.0-1.5*0.5 ; c2 = 1.0-0.5 ; c3 = 2.0 ; slope = 0.5

      do J=Jsq,Jeq+1 ; do i=is-1,ie
        uhc = 0.5 * (G%dy_Cu(I,j) * u(I,j,k)) * (h(i,j,k) + h(i+1,j,k))
        uhm = uh(I,j,k)
        ! This sometimes matters with some types of open boundary conditions.
        if (G%dy_Cu(I,j) == 0.0) uhc = uhm

        if (abs(uhc) < 0.1*abs(uhm)) then
          uhm = 10.0*uhc
        elseif (abs(uhc) > c1*abs(uhm)) then
          if (abs(uhc) < c2*abs(uhm)) then ; uhc = (3.0*uhc+(1.0-c2*3.0)*uhm)
          elseif (abs(uhc) <= c3*abs(uhm)) then ; uhc = uhm
          else ; uhc = slope*uhc+(1.0-c3*slope)*uhm
          endif
        endif

        if (uhc > uhm) then
          uh_min(i,j) = uhm ; uh_max(i,j) = uhc
        else
          uh_max(i,j) = uhm ; uh_min(i,j) = uhc
        endif
      enddo ; enddo
      do j=js-1,je ; do I=Isq,Ieq+1
        vhc = 0.5 * (G%dx_Cv(i,J) * v(i,J,k)) * (h(i,j,k) + h(i,j+1,k))
        vhm = vh(i,J,k)
        ! This sometimes matters with some types of open boundary conditions.
        if (G%dx_Cv(i,J) == 0.0) vhc = vhm

        if (abs(vhc) < 0.1*abs(vhm)) then
          vhm = 10.0*vhc
        elseif (abs(vhc) > c1*abs(vhm)) then
          if (abs(vhc) < c2*abs(vhm)) then ; vhc = (3.0*vhc+(1.0-c2*3.0)*vhm)
          else if (abs(vhc) <= c3*abs(vhm)) then ; vhc = vhm
          else ; vhc = slope*vhc+(1.0-c3*slope)*vhm
          endif
        endif

        if (vhc > vhm) then
          vh_min(i,j) = vhm ; vh_max(i,j) = vhc
        else
          vh_max(i,j) = vhm ; vh_min(i,j) = vhc
        endif
      enddo ; enddo
    endif

    ! Calculate KE and the gradient of KE
    call gradKE(u, v, h, uh, vh, KE, KEx, KEy, k, OBC, G, CS)

    ! Calculate the tendencies of zonal velocity due to the Coriolis
    ! force and momentum advection.  On a Cartesian grid, this is
    !     CAu =  q * vh - d(KE)/dx.
    if (CS%Coriolis_Scheme == SADOURNY75_ENERGY) then
      if (CS%Coriolis_En_Dis) then
        ! Energy dissipating biased scheme, Hallberg 200x
        do j=js,je ; do I=Isq,Ieq
          if (q(I,J)*u(I,j,k) == 0.0) then
            temp1 = q(I,J) * ( (vh_max(i,j)+vh_max(i+1,j)) &
                             + (vh_min(i,j)+vh_min(i+1,j)) )*0.5
          elseif (q(I,J)*u(I,j,k) < 0.0) then
            temp1 = q(I,J) * (vh_max(i,j)+vh_max(i+1,j))
          else
            temp1 = q(I,J) * (vh_min(i,j)+vh_min(i+1,j))
          endif
          if (q(I,J-1)*u(I,j,k) == 0.0) then
            temp2 = q(I,J-1) * ( (vh_max(i,j-1)+vh_max(i+1,j-1)) &
                               + (vh_min(i,j-1)+vh_min(i+1,j-1)) )*0.5
          elseif (q(I,J-1)*u(I,j,k) < 0.0) then
            temp2 = q(I,J-1) * (vh_max(i,j-1)+vh_max(i+1,j-1))
          else
            temp2 = q(I,J-1) * (vh_min(i,j-1)+vh_min(i+1,j-1))
          endif
          CAu(I,j,k) = 0.25 * G%IdxCu(I,j) * (temp1 + temp2)
        enddo ; enddo
      else
        ! Energy conserving scheme, Sadourny 1975
        do j=js,je ; do I=Isq,Ieq
          CAu(I,j,k) = 0.25 * &
            (q(I,J) * (vh(i+1,J,k) + vh(i,J,k)) + &
             q(I,J-1) * (vh(i,J-1,k) + vh(i+1,J-1,k))) * G%IdxCu(I,j)
        enddo ; enddo
      endif
    elseif (CS%Coriolis_Scheme == SADOURNY75_ENSTRO) then
      do j=js,je ; do I=Isq,Ieq
        CAu(I,j,k) = 0.125 * (G%IdxCu(I,j) * (q(I,J) + q(I,J-1))) * &
                     ((vh(i+1,J,k) + vh(i,J,k)) + (vh(i,J-1,k) + vh(i+1,J-1,k)))
      enddo ; enddo
    elseif ((CS%Coriolis_Scheme == ARAKAWA_HSU90) .or. &
            (CS%Coriolis_Scheme == ARAKAWA_LAMB81) .or. &
            (CS%Coriolis_Scheme == AL_BLEND)) then
      ! (Global) Energy and (Local) Enstrophy conserving, Arakawa & Hsu 1990
      do j=js,je ; do I=Isq,Ieq
        CAu(I,j,k) = ((a(I,j) * vh(i+1,J,k) + &
                       c(I,j) * vh(i,J-1,k))  &
                    + (b(I,j) * vh(i,J,k) +   &
                       d(I,j) * vh(i+1,J-1,k))) * G%IdxCu(I,j)
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == ROBUST_ENSTRO) then
      ! An enstrophy conserving scheme robust to vanishing layers
      ! Note: Heffs are in lieu of h_at_v that should be returned by the
      !       continuity solver. AJA
      do j=js,je ; do I=Isq,Ieq
        Heff1 = abs(vh(i,J,k)*G%IdxCv(i,J))/(eps_vel+abs(v(i,J,k)))
        Heff1 = max(Heff1,min(h(i,j,k),h(i,j+1,k)))
        Heff1 = min(Heff1,max(h(i,j,k),h(i,j+1,k)))
        Heff2 = abs(vh(i,J-1,k)*G%IdxCv(i,J-1))/(eps_vel+abs(v(i,J-1,k)))
        Heff2 = max(Heff2,min(h(i,j-1,k),h(i,j,k)))
        Heff2 = min(Heff2,max(h(i,j-1,k),h(i,j,k)))
        Heff3 = abs(vh(i+1,J,k)*G%IdxCv(i+1,J))/(eps_vel+abs(v(i+1,J,k)))
        Heff3 = max(Heff3,min(h(i+1,j,k),h(i+1,j+1,k)))
        Heff3 = min(Heff3,max(h(i+1,j,k),h(i+1,j+1,k)))
        Heff4 = abs(vh(i+1,J-1,k)*G%IdxCv(i+1,J-1))/(eps_vel+abs(v(i+1,J-1,k)))
        Heff4 = max(Heff4,min(h(i+1,j-1,k),h(i+1,j,k)))
        Heff4 = min(Heff4,max(h(i+1,j-1,k),h(i+1,j,k)))
        if (CS%PV_Adv_Scheme == PV_ADV_CENTERED) then
          CAu(I,j,k) = 0.5*(abs_vort(I,J)+abs_vort(I,J-1)) * &
                       ((vh(i  ,J  ,k)+vh(i+1,J-1,k)) +      &
                        (vh(i  ,J-1,k)+vh(i+1,J  ,k)) ) /    &
                       (h_tiny +((Heff1+Heff4) +(Heff2+Heff3)) ) * G%IdxCu(I,j)
        elseif (CS%PV_Adv_Scheme == PV_ADV_UPWIND1) then
          VHeff = ((vh(i  ,J  ,k)+vh(i+1,J-1,k)) +      &
                   (vh(i  ,J-1,k)+vh(i+1,J  ,k)) )
          QVHeff = 0.5*( (abs_vort(I,J)+abs_vort(I,J-1))*VHeff &
                        -(abs_vort(I,J)-abs_vort(I,J-1))*abs(VHeff) )
          CAu(I,j,k) = QVHeff / &
                     (h_tiny +((Heff1+Heff4) +(Heff2+Heff3)) ) * G%IdxCu(I,j)
        endif
      enddo ; enddo
    endif
    ! Add in the additonal terms with Arakawa & Lamb.
    if ((CS%Coriolis_Scheme == ARAKAWA_LAMB81) .or. &
        (CS%Coriolis_Scheme == AL_BLEND)) then ; do j=js,je ; do I=Isq,Ieq
      CAu(I,j,k) = CAu(I,j,k) + &
            (ep_u(i,j)*uh(I-1,j,k) - ep_u(i+1,j)*uh(I+1,j,k)) * G%IdxCu(I,j)
    enddo ; enddo ; endif


    if (CS%bound_Coriolis) then
      do j=js,je ; do I=Isq,Ieq
        max_fv = MAX(max_fvq(I,J),max_fvq(I,J-1))
        min_fv = MIN(min_fvq(I,J),min_fvq(I,J-1))
       ! CAu(I,j,k) = min( CAu(I,j,k), max_fv )
       ! CAu(I,j,k) = max( CAu(I,j,k), min_fv )
        if (CAu(I,j,k) > max_fv) then
            CAu(I,j,k) = max_fv
        else
          if (CAu(I,j,k) < min_fv) CAu(I,j,k) = min_fv
        endif
      enddo ; enddo
    endif

    ! Term - d(KE)/dx.
    do j=js,je ; do I=Isq,Ieq
      CAu(I,j,k) = CAu(I,j,k) - KEx(I,j)
      if (ASSOCIATED(AD%gradKEu)) AD%gradKEu(I,j,k) = -KEx(I,j)
    enddo ; enddo


    ! Calculate the tendencies of meridional velocity due to the Coriolis
    ! force and momentum advection.  On a Cartesian grid, this is
    !     CAv = - q * uh - d(KE)/dy.
    if (CS%Coriolis_Scheme == SADOURNY75_ENERGY) then
      if (CS%Coriolis_En_Dis) then
        ! Energy dissipating biased scheme, Hallberg 200x
        do J=Jsq,Jeq ; do i=is,ie
          if (q(I-1,J)*v(i,J,k) == 0.0) then
            temp1 = q(I-1,J) * ( (uh_max(i-1,j)+uh_max(i-1,j+1)) &
                               + (uh_min(i-1,j)+uh_min(i-1,j+1)) )*0.5
          elseif (q(I-1,J)*v(i,J,k) > 0.0) then
            temp1 = q(I-1,J) * (uh_max(i-1,j)+uh_max(i-1,j+1))
          else
            temp1 = q(I-1,J) * (uh_min(i-1,j)+uh_min(i-1,j+1))
          endif
          if (q(I,J)*v(i,J,k) == 0.0) then
            temp2 = q(I,J) * ( (uh_max(i,j)+uh_max(i,j+1)) &
                             + (uh_min(i,j)+uh_min(i,j+1)) )*0.5
          elseif (q(I,J)*v(i,J,k) > 0.0) then
            temp2 = q(I,J) * (uh_max(i,j)+uh_max(i,j+1))
          else
            temp2 = q(I,J) * (uh_min(i,j)+uh_min(i,j+1))
          endif
          CAv(i,J,k) = - 0.25 * G%IdyCv(i,J) * (temp1 + temp2)
        enddo ; enddo
      else
        ! Energy conserving scheme, Sadourny 1975
        do J=Jsq,Jeq ; do i=is,ie
          CAv(i,J,k) = - 0.25* &
              (q(I-1,J)*(uh(I-1,j,k) + uh(I-1,j+1,k)) + &
               q(I,J)*(uh(I,j,k) + uh(I,j+1,k))) * G%IdyCv(i,J)
        enddo ; enddo
      endif
    elseif (CS%Coriolis_Scheme == SADOURNY75_ENSTRO) then
      do J=Jsq,Jeq ; do i=is,ie
        CAv(i,J,k) = -0.125 * (G%IdyCv(i,J) * (q(I-1,J) + q(I,J))) * &
                     ((uh(I-1,j,k) + uh(I-1,j+1,k)) + (uh(I,j,k) + uh(I,j+1,k)))
      enddo ; enddo
    elseif ((CS%Coriolis_Scheme == ARAKAWA_HSU90) .or. &
            (CS%Coriolis_Scheme == ARAKAWA_LAMB81) .or. &
            (CS%Coriolis_Scheme == AL_BLEND)) then
      ! (Global) Energy and (Local) Enstrophy conserving, Arakawa & Hsu 1990
      do J=Jsq,Jeq ; do i=is,ie
        CAv(i,J,k) = - ((a(I-1,j)   * uh(I-1,j,k) + &
                         c(I,j+1)   * uh(I,j+1,k))  &
                      + (b(I,j)     * uh(I,j,k) +   &
                         d(I-1,j+1) * uh(I-1,j+1,k))) * G%IdyCv(i,J)
      enddo ; enddo
    elseif (CS%Coriolis_Scheme == ROBUST_ENSTRO) then
      ! An enstrophy conserving scheme robust to vanishing layers
      ! Note: Heffs are in lieu of h_at_u that should be returned by the
      !       continuity solver. AJA
      do J=Jsq,Jeq ; do i=is,ie
        Heff1 = abs(uh(I,j,k)*G%IdyCu(I,j))/(eps_vel+abs(u(I,j,k)))
        Heff1 = max(Heff1,min(h(i,j,k),h(i+1,j,k)))
        Heff1 = min(Heff1,max(h(i,j,k),h(i+1,j,k)))
        Heff2 = abs(uh(I-1,j,k)*G%IdyCu(I-1,j))/(eps_vel+abs(u(I-1,j,k)))
        Heff2 = max(Heff2,min(h(i-1,j,k),h(i,j,k)))
        Heff2 = min(Heff2,max(h(i-1,j,k),h(i,j,k)))
        Heff3 = abs(uh(I,j+1,k)*G%IdyCu(I,j+1))/(eps_vel+abs(u(I,j+1,k)))
        Heff3 = max(Heff3,min(h(i,j+1,k),h(i+1,j+1,k)))
        Heff3 = min(Heff3,max(h(i,j+1,k),h(i+1,j+1,k)))
        Heff4 = abs(uh(I-1,j+1,k)*G%IdyCu(I-1,j+1))/(eps_vel+abs(u(I-1,j+1,k)))
        Heff4 = max(Heff4,min(h(i-1,j+1,k),h(i,j+1,k)))
        Heff4 = min(Heff4,max(h(i-1,j+1,k),h(i,j+1,k)))
        if (CS%PV_Adv_Scheme == PV_ADV_CENTERED) then
          CAv(i,J,k) = - 0.5*(abs_vort(I,J)+abs_vort(I-1,J)) * &
                         ((uh(I  ,j  ,k)+uh(I-1,j+1,k)) +      &
                          (uh(I-1,j  ,k)+uh(I  ,j+1,k)) ) /    &
                      (h_tiny + ((Heff1+Heff4) +(Heff2+Heff3)) ) * G%IdyCv(i,J)
        elseif (CS%PV_Adv_Scheme == PV_ADV_UPWIND1) then
          UHeff = ((uh(I  ,j  ,k)+uh(I-1,j+1,k)) +      &
                   (uh(I-1,j  ,k)+uh(I  ,j+1,k)) )
          QUHeff = 0.5*( (abs_vort(I,J)+abs_vort(I-1,J))*UHeff &
                        -(abs_vort(I,J)-abs_vort(I-1,J))*abs(UHeff) )
          CAv(i,J,k) = - QUHeff / &
                       (h_tiny + ((Heff1+Heff4) +(Heff2+Heff3)) ) * G%IdyCv(i,J)
        endif
      enddo ; enddo
    endif
    ! Add in the additonal terms with Arakawa & Lamb.
    if ((CS%Coriolis_Scheme == ARAKAWA_LAMB81) .or. &
        (CS%Coriolis_Scheme == AL_BLEND)) then ; do J=Jsq,Jeq ; do i=is,ie
      CAv(i,J,k) = CAv(i,J,k) + &
            (ep_v(i,j)*vh(i,J-1,k) - ep_v(i,j+1)*vh(i,J+1,k)) * G%IdyCv(i,J)
    enddo ; enddo ; endif

    if (CS%bound_Coriolis) then
      do J=Jsq,Jeq ; do i=is,ie
        max_fu = MAX(max_fuq(I,J),max_fuq(I-1,J))
        min_fu = MIN(min_fuq(I,J),min_fuq(I-1,J))
        if (CAv(i,J,k) > max_fu) then
            CAv(i,J,k) = max_fu
        else
          if (CAv(i,J,k) < min_fu) CAv(i,J,k) = min_fu
        endif
      enddo ; enddo
    endif

    ! Term - d(KE)/dy.
    do J=Jsq,Jeq ; do i=is,ie
      CAv(i,J,k) = CAv(i,J,k) - KEy(i,J)
      if (ASSOCIATED(AD%gradKEv)) AD%gradKEv(i,J,k) = -KEy(i,J)
    enddo ; enddo

    if (ASSOCIATED(AD%rv_x_u) .or. ASSOCIATED(AD%rv_x_v)) then
      ! Calculate the Coriolis-like acceleration due to relative vorticity.
      if (CS%Coriolis_Scheme == SADOURNY75_ENERGY) then
        if (ASSOCIATED(AD%rv_x_u)) then
          do J=Jsq,Jeq ; do i=is,ie
            AD%rv_x_u(i,J,k) = - 0.25* &
              (q2(I-1,j)*(uh(I-1,j,k) + uh(I-1,j+1,k)) + &
               q2(I,j)*(uh(I,j,k) + uh(I,j+1,k))) * G%IdyCv(i,J)
          enddo ; enddo
        endif

        if (ASSOCIATED(AD%rv_x_v)) then
          do j=js,je ; do I=Isq,Ieq
            AD%rv_x_v(I,j,k) = 0.25 * &
              (q2(I,j) * (vh(i+1,J,k) + vh(i,J,k)) + &
               q2(I,j-1) * (vh(i,J-1,k) + vh(i+1,J-1,k))) * G%IdxCu(I,j)
          enddo ; enddo
        endif
      else
        if (ASSOCIATED(AD%rv_x_u)) then
          do J=Jsq,Jeq ; do i=is,ie
            AD%rv_x_u(i,J,k) = -G%IdyCv(i,J) * C1_12 * &
              ((q2(I,J) + q2(I-1,J) + q2(I-1,J-1)) * uh(I-1,j,k) + &
               (q2(I-1,J) + q2(I,J) + q2(I,J-1)) * uh(I,j,k) + &
               (q2(I-1,J) + q2(I,J+1) + q2(I,J)) * uh(I,j+1,k) + &
               (q2(I,J) + q2(I-1,J+1) + q2(I-1,J)) * uh(I-1,j+1,k))
          enddo ; enddo
        endif

        if (ASSOCIATED(AD%rv_x_v)) then
          do j=js,je ; do I=Isq,Ieq
            AD%rv_x_v(I,j,k) = G%IdxCu(I,j) * C1_12 * &
              ((q2(I+1,J) + q2(I,J) + q2(I,J-1)) * vh(i+1,J,k) + &
               (q2(I-1,J) + q2(I,J) + q2(I,J-1)) * vh(i,J,k) + &
               (q2(I-1,J-1) + q2(I,J) + q2(I,J-1)) * vh(i,J-1,k) + &
               (q2(I+1,J-1) + q2(I,J) + q2(I,J-1)) * vh(i+1,J-1,k))
          enddo ; enddo
        endif
      endif
    endif

  enddo ! k-loop.
!$OMP end parallel

  ! Here the various Coriolis-related derived quantities are offered for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_rv > 0) call post_data(CS%id_rv, RV, CS%diag)
    if (CS%id_PV > 0) call post_data(CS%id_PV, PV, CS%diag)
    if (CS%id_gKEu>0) call post_data(CS%id_gKEu, AD%gradKEu, CS%diag)
    if (CS%id_gKEv>0) call post_data(CS%id_gKEv, AD%gradKEv, CS%diag)
    if (CS%id_rvxu > 0) call post_data(CS%id_rvxu, AD%rv_x_u, CS%diag)
    if (CS%id_rvxv > 0) call post_data(CS%id_rvxv, AD%rv_x_v, CS%diag)
  endif

end subroutine CorAdCalc


!> Calculates the acceleration due to the gradient of kinetic energy.
subroutine gradKE(u, v, h, uh, vh, KE, KEx, KEy, k, OBC, G, CS)
  type(ocean_grid_type),                      intent(in)  :: G !< Ocen grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)),  intent(in)  :: u !< Zonal velocity (m/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)),  intent(in)  :: v !< Meridional velocity (m/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)  :: h !< Layer thickness (m or kg/m2)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)),  intent(in)  :: uh !< Zonal transport u*h*dy (m3/s or kg/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)),  intent(in)  :: vh !< Meridional transport v*h*dx (m3/s or kg/s)
  real, dimension(SZI_(G) ,SZJ_(G) ),         intent(out) :: KE !< Kinetic energy (m2/s2)
  real, dimension(SZIB_(G),SZJ_(G) ),         intent(out) :: KEx !< Zonal acceleration due to kinetic
                                                                 !! energy gradient (m/s2)
  real, dimension(SZI_(G) ,SZJB_(G)),         intent(out) :: KEy !< Meridional acceleration due to kinetic
                                                                 !! energy gradient (m/s2)
  integer,                                    intent(in)  :: k !< Layer number to calculate for
  type(ocean_OBC_type),                       pointer     :: OBC !< Open boundary control structure
  type(CoriolisAdv_CS),                       pointer     :: CS !< Control structure for MOM_CoriolisAdv
  ! Local variables
  real :: um, up, vm, vp         ! Temporary variables with units of m s-1.
  real :: um2, up2, vm2, vp2     ! Temporary variables with units of m2 s-2.
  real :: um2a, up2a, vm2a, vp2a ! Temporary variables with units of m4 s-2.
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, n

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB


  ! Calculate KE (Kinetic energy for use in the -grad(KE) acceleration term).
  if (CS%KE_Scheme.eq.KE_ARAKAWA) then
    ! The following calculation of Kinetic energy includes the metric terms
    ! identified in Arakawa & Lamb 1982 as important for KE conservation.  It
    ! also includes the possibility of partially-blocked tracer cell faces.
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      KE(i,j) = ( ( G%areaCu( I ,j)*(u( I ,j,k)*u( I ,j,k))   &
                   +G%areaCu(I-1,j)*(u(I-1,j,k)*u(I-1,j,k)) ) &
                 +( G%areaCv(i, J )*(v(i, J ,k)*v(i, J ,k))   &
                   +G%areaCv(i,J-1)*(v(i,J-1,k)*v(i,J-1,k)) ) &
                )*0.25*G%IareaT(i,j)
    enddo ; enddo
  elseif (CS%KE_Scheme.eq.KE_SIMPLE_GUDONOV) then
    ! The following discretization of KE is based on the one-dimensinal Gudonov
    ! scheme which does not take into account any geometric factors
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      up = 0.5*( u(I-1,j,k) + ABS( u(I-1,j,k) ) ) ; up2 = up*up
      um = 0.5*( u( I ,j,k) - ABS( u( I ,j,k) ) ) ; um2 = um*um
      vp = 0.5*( v(i,J-1,k) + ABS( v(i,J-1,k) ) ) ; vp2 = vp*vp
      vm = 0.5*( v(i, J ,k) - ABS( v(i, J ,k) ) ) ; vm2 = vm*vm
      KE(i,j) = ( max(up2,um2) + max(vp2,vm2) ) *0.5
    enddo ; enddo
  elseif (CS%KE_Scheme.eq.KE_GUDONOV) then
    ! The following discretization of KE is based on the one-dimensinal Gudonov
    ! scheme but has been adapted to take horizontal grid factors into account
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      up = 0.5*( u(I-1,j,k) + ABS( u(I-1,j,k) ) ) ; up2a = up*up*G%areaCu(I-1,j)
      um = 0.5*( u( I ,j,k) - ABS( u( I ,j,k) ) ) ; um2a = um*um*G%areaCu( I ,j)
      vp = 0.5*( v(i,J-1,k) + ABS( v(i,J-1,k) ) ) ; vp2a = vp*vp*G%areaCv(i,J-1)
      vm = 0.5*( v(i, J ,k) - ABS( v(i, J ,k) ) ) ; vm2a = vm*vm*G%areaCv(i, J )
      KE(i,j) = ( max(um2a,up2a) + max(vm2a,vp2a) )*0.5*G%IareaT(i,j)
    enddo ; enddo
  endif

  ! Term - d(KE)/dx.
  do j=js,je ; do I=Isq,Ieq
    KEx(I,j) = (KE(i+1,j) - KE(i,j)) * G%IdxCu(I,j)
  enddo ; enddo

  ! Term - d(KE)/dy.
  do J=Jsq,Jeq ; do i=is,ie
    KEy(i,J) = (KE(i,j+1) - KE(i,j)) * G%IdyCv(i,J)
  enddo ; enddo

  if (associated(OBC)) then
    do n=1,OBC%number_of_segments
      if (OBC%segment(n)%is_N_or_S) then
        do i=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
          KEy(i,OBC%segment(n)%HI%JsdB) = 0.
        enddo
      elseif (OBC%segment(n)%is_E_or_W) then
        do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
          KEx(OBC%segment(n)%HI%IsdB,j) = 0.
        enddo
      endif
    enddo
  endif

end subroutine gradKE

!> Initializes the control structure for coriolisadv_cs
subroutine CoriolisAdv_init(Time, G, param_file, diag, AD, CS)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G !< Ocean grid structure
  type(param_file_type),   intent(in)    :: param_file !< Runtime parameter handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(accel_diag_ptrs),   target, intent(inout) :: AD !< Strorage for acceleration diagnostics
  type(CoriolisAdv_CS),    pointer       :: CS !< Control structure fro MOM_CoriolisAdv
  ! Local variables
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_CoriolisAdv" ! This module's name.
  character(len=20)  :: tmpstr
  character(len=400) :: mesg
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "CoriolisAdv_init called with associated control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag ; CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "NOSLIP", CS%no_slip, &
                 "If true, no slip boundary conditions are used; otherwise \n"//&
                 "free slip boundary conditions are assumed. The \n"//&
                 "implementation of the free slip BCs on a C-grid is much \n"//&
                 "cleaner than the no slip BCs. The use of free slip BCs \n"//&
                 "is strongly encouraged, and no slip BCs are not used with \n"//&
                 "the biharmonic viscosity.", default=.false.)

  call get_param(param_file, mod, "CORIOLIS_EN_DIS", CS%Coriolis_En_Dis, &
                 "If true, two estimates of the thickness fluxes are used \n"//&
                 "to estimate the Coriolis term, and the one that \n"//&
                 "dissipates energy relative to the other one is used.", &
                 default=.false.)

  ! Set %Coriolis_Scheme
  ! (Select the baseline discretization for the Coriolis term)
  call get_param(param_file, mod, "CORIOLIS_SCHEME", tmpstr, &
                 "CORIOLIS_SCHEME selects the discretization for the \n"//&
                 "Coriolis terms. Valid values are: \n"//&
                 "\t SADOURNY75_ENERGY - Sadourny, 1975; energy cons. \n"//&
                 "\t ARAKAWA_HSU90     - Arakawa & Hsu, 1990 \n"//&
                 "\t SADOURNY75_ENSTRO - Sadourny, 1975; enstrophy cons. \n"//&
                 "\t ARAKAWA_LAMB81    - Arakawa & Lamb, 1981; En. + Enst.\n"//&
                 "\t ARAKAWA_LAMB_BLEND - A blend of Arakawa & Lamb with \n"//&
                 "\t                      Arakawa & Hsu and Sadourny energy", &
                 default=SADOURNY75_ENERGY_STRING)
  tmpstr = uppercase(tmpstr)
  select case (tmpstr)
    case (SADOURNY75_ENERGY_STRING)
      CS%Coriolis_Scheme = SADOURNY75_ENERGY
    case (ARAKAWA_HSU_STRING)
      CS%Coriolis_Scheme = ARAKAWA_HSU90
    case (SADOURNY75_ENSTRO_STRING)
      CS%Coriolis_Scheme = SADOURNY75_ENSTRO
    case (ARAKAWA_LAMB_STRING)
      CS%Coriolis_Scheme = ARAKAWA_LAMB81
    case (AL_BLEND_STRING)
      CS%Coriolis_Scheme = AL_BLEND
    case (ROBUST_ENSTRO_STRING)
      CS%Coriolis_Scheme = ROBUST_ENSTRO
      CS%Coriolis_En_Dis = .false.
    case default
      call MOM_mesg('CoriolisAdv_init: Coriolis_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: Unrecognized setting "// &
            "#define CORIOLIS_SCHEME "//trim(tmpstr)//" found in input file.")
  end select
  if (CS%Coriolis_Scheme == AL_BLEND) then
    call get_param(param_file, mod, "CORIOLIS_BLEND_WT_LIN", CS%wt_lin_blend, &
                 "A weighting value for the ratio of inverse thicknesses, \n"//&
                 "beyond which the blending between Sadourny Energy and \n"//&
                 "Arakawa & Hsu goes linearly to 0 when CORIOLIS_SCHEME \n"//&
                 "is ARAWAKA_LAMB_BLEND. This must be between 1 and 1e-16.", &
                 units="nondim", default=0.125)
    call get_param(param_file, mod, "CORIOLIS_BLEND_F_EFF_MAX", CS%F_eff_max_blend, &
                 "The factor by which the maximum effective Coriolis \n"//&
                 "acceleration from any point can be increased when \n"//&
                 "blending different discretizations with the \n"//&
                 "ARAKAWA_LAMB_BLEND Coriolis scheme.  This must be \n"//&
                 "greater than 2.0 (the max value for Sadourny energy).", &
                 units="nondim", default=4.0)
    CS%wt_lin_blend = min(1.0, max(CS%wt_lin_blend,1e-16))
    if (CS%F_eff_max_blend < 2.0) call MOM_error(WARNING, "CoriolisAdv_init: "//&
           "CORIOLIS_BLEND_F_EFF_MAX should be at least 2.")
  endif

  mesg = "If true, the Coriolis terms at u-points are bounded by \n"//&
         "the four estimates of (f+rv)v from the four neighboring \n"//&
         "v-points, and similarly at v-points."
  if (CS%Coriolis_En_Dis .and. (CS%Coriolis_Scheme == SADOURNY75_ENERGY)) then
    mesg = trim(mesg)//"  This option is \n"//&
                 "always effectively false with CORIOLIS_EN_DIS defined and \n"//&
                 "CORIOLIS_SCHEME set to "//trim(SADOURNY75_ENERGY_STRING)//"."
  else
    mesg = trim(mesg)//"  This option would \n"//&
                 "have no effect on the SADOURNY Coriolis scheme if it \n"//&
                 "were possible to use centered difference thickness fluxes."
  endif
  call get_param(param_file, mod, "BOUND_CORIOLIS", CS%bound_Coriolis, mesg, &
                 default=.false.)
  if ((CS%Coriolis_En_Dis .and. (CS%Coriolis_Scheme == SADOURNY75_ENERGY)) .or. &
      (CS%Coriolis_Scheme == ROBUST_ENSTRO)) CS%bound_Coriolis = .false.

  ! Set KE_Scheme (selects discretization of KE)
  call get_param(param_file, mod, "KE_SCHEME", tmpstr, &
                 "KE_SCHEME selects the discretization for acceleration \n"//&
                 "due to the kinetic energy gradient. Valid values are: \n"//&
                 "\t KE_ARAKAWA, KE_SIMPLE_GUDONOV, KE_GUDONOV", &
                 default=KE_ARAKAWA_STRING)
  tmpstr = uppercase(tmpstr)
  select case (tmpstr)
    case (KE_ARAKAWA_STRING); CS%KE_Scheme = KE_ARAKAWA
    case (KE_SIMPLE_GUDONOV_STRING); CS%KE_Scheme = KE_SIMPLE_GUDONOV
    case (KE_GUDONOV_STRING); CS%KE_Scheme = KE_GUDONOV
    case default
      call MOM_mesg('CoriolisAdv_init: KE_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: "// &
               "#define KE_SCHEME "//trim(tmpstr)//" in input file is invalid.")
  end select

  ! Set PV_Adv_Scheme (selects discretization of PV advection)
  call get_param(param_file, mod, "PV_ADV_SCHEME", tmpstr, &
                 "PV_ADV_SCHEME selects the discretization for PV \n"//&
                 "advection. Valid values are: \n"//&
                 "\t PV_ADV_CENTERED - centered (aka Sadourny, 75) \n"//&
                 "\t PV_ADV_UPWIND1  - upwind, first order", &
                 default=PV_ADV_CENTERED_STRING)
  select case (uppercase(tmpstr))
    case (PV_ADV_CENTERED_STRING)
      CS%PV_Adv_Scheme = PV_ADV_CENTERED
    case (PV_ADV_UPWIND1_STRING)
      CS%PV_Adv_Scheme = PV_ADV_UPWIND1
    case default
      call MOM_mesg('CoriolisAdv_init: PV_Adv_Scheme ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "CoriolisAdv_init: "// &
                     "#DEFINE PV_ADV_SCHEME in input file is invalid.")
  end select

  CS%id_rv = register_diag_field('ocean_model', 'RV', diag%axesBL, Time, &
     'Relative Vorticity', 'second-1')

  CS%id_PV = register_diag_field('ocean_model', 'PV', diag%axesBL, Time, &
     'Potential Vorticity', 'meter-1 second-1')

  CS%id_gKEu = register_diag_field('ocean_model', 'gKEu', diag%axesCuL, Time, &
     'Zonal Acceleration from Grad. Kinetic Energy', 'meter-1 second-2')
  if (CS%id_gKEu > 0) call safe_alloc_ptr(AD%gradKEu,IsdB,IedB,jsd,jed,nz)

  CS%id_gKEv = register_diag_field('ocean_model', 'gKEv', diag%axesCvL, Time, &
     'Meridional Acceleration from Grad. Kinetic Energy', 'meter-1 second-2')
  if (CS%id_gKEv > 0) call safe_alloc_ptr(AD%gradKEv,isd,ied,JsdB,JedB,nz)

  CS%id_rvxu = register_diag_field('ocean_model', 'rvxu', diag%axesCvL, Time, &
     'Meridional Acceleration from Relative Vorticity', 'meter-1 second-2')
  if (CS%id_rvxu > 0) call safe_alloc_ptr(AD%rv_x_u,isd,ied,JsdB,JedB,nz)

  CS%id_rvxv = register_diag_field('ocean_model', 'rvxv', diag%axesCuL, Time, &
     'Zonal Acceleration from Relative Vorticity', 'meter-1 second-2')
  if (CS%id_rvxv > 0) call safe_alloc_ptr(AD%rv_x_v,IsdB,IedB,jsd,jed,nz)

end subroutine CoriolisAdv_init

!> Destructor for coriolisadv_cs
subroutine CoriolisAdv_end(CS)
  type(CoriolisAdv_CS), pointer :: CS !< Control structure fro MOM_CoriolisAdv
  deallocate(CS)
end subroutine CoriolisAdv_end

!> \namespace mom_coriolisadv
!!
!! This file contains the subroutine that calculates the time
!! derivatives of the velocities due to Coriolis acceleration and
!! momentum advection.  This subroutine uses either a vorticity
!! advection scheme from Arakawa and Hsu, Mon. Wea. Rev. 1990, or
!! Sadourny's (JAS 1975) energy conserving scheme.  Both have been
!! modified to use general orthogonal coordinates as described in
!! Arakawa and Lamb, Mon. Wea. Rev. 1981.  Both schemes are second
!! order accurate, and allow for vanishingly small layer thicknesses.
!! The Arakawa and Hsu scheme globally conserves both total energy
!! and potential enstrophy in the limit of nondivergent flow.
!! Sadourny's energy conserving scheme conserves energy if the flow
!! is nondivergent or centered difference thickness fluxes are used.
!!
!! Two sets of boundary conditions have been coded in the
!! definition of relative vorticity.  These are written as:
!! NOSLIP defined (in spherical coordinates):
!!   relvort = dv/dx (east & west), with v = 0.
!!   relvort = -sec(Q) * d(u cos(Q))/dy (north & south), with u = 0.
!!
!! NOSLIP not defined (free slip):
!!   relvort = 0 (all boundaries)
!!
!! with Q temporarily defined as latitude.  The free slip boundary
!! condition is much more natural on a C-grid.
!!
!! A small fragment of the grid is shown below:
!! \verbatim
!!
!!    j+1  x ^ x ^ x   At x:  q, CoriolisBu
!!    j+1  > o > o >   At ^:  v, CAv, vh
!!    j    x ^ x ^ x   At >:  u, CAu, uh, a, b, c, d
!!    j    > o > o >   At o:  h, KE
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!! \endverbatim
!!
!! The boundaries always run through q grid points (x).

end module MOM_CoriolisAdv
