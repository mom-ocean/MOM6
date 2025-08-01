!> Calculations of isoneutral slopes and stratification.
module MOM_isopycnal_slopes

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,     only : hchksum, uvchksum
use MOM_error_handler, only : MOM_error, FATAL
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_EOS,           only : calculate_density_derivs, calculate_density_second_derivs, EOS_domain
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S

implicit none ; private

#include <MOM_memory.h>

public calc_isoneutral_slopes, vert_fill_TS

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> Calculate isopycnal slopes, and optionally return other stratification dependent functions such as N^2
!! and dz*S^2*g-prime used, or calculable from factors used, during the calculation.
subroutine calc_isoneutral_slopes(G, GV, US, h, e, tv, dt_kappa_smooth, use_stanley, slope_x, slope_y, &
                                  N2_u, N2_v, dzu, dzv, dzSxN, dzSyN, halo, OBC, OBC_N2)
  type(ocean_grid_type),                       intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                     intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                       intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: e    !< Interface heights [Z ~> m]
  type(thermo_var_ptrs),                       intent(in)    :: tv   !< A structure pointing to various
                                                                     !! thermodynamic variables
  real,                                        intent(in)    :: dt_kappa_smooth !< A smoothing vertical
                                                                     !! diffusivity times a smoothing
                                                                     !! timescale [H Z ~> m2 or kg m-1]
  logical,                                     intent(in)    :: use_stanley !< turn on stanley param in slope
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: slope_x !< Isopycnal slope in i-dir [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: slope_y !< Isopycnal slope in j-dir [Z L-1 ~> nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: N2_u !< Brunt-Vaisala frequency squared at
                                                                     !! interfaces between u-points [L2 Z-2 T-2 ~> s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: N2_v !< Brunt-Vaisala frequency squared at
                                                                     !! interfaces between v-points [L2 Z-2 T-2 ~> s-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: dzu  !< Z-thickness at u-points [Z ~> m]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: dzv  !< Z-thickness at v-points [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: dzSxN !< Z-thickness times zonal slope contribution to
                                                                     !! Eady growth rate at u-points. [Z T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: dzSyN !< Z-thickness times meridional slope contrib. to
                                                                     !! Eady growth rate at v-points. [Z T-1 ~> m s-1]
  integer,                           optional, intent(in)    :: halo !< Halo width over which to compute
  type(ocean_OBC_type),              optional, pointer       :: OBC  !< Open boundaries control structure.
  logical,                           optional, intent(in)    :: OBC_N2 !< If present and true, use interior data
                                                                     !! to calculate stratification at open boundary
                                                                     !! condition faces.

  ! Local variables
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)) :: &
    T, &          ! The temperature [C ~> degC], with the values in
                  ! in massless layers filled vertically by diffusion.
    S             ! The filled salinity [S ~> ppt], with the values in
                  ! in massless layers filled vertically by diffusion.
  real, dimension(SZI_(G), SZJ_(G),SZK_(GV)+1) :: &
    pres          ! The pressure at an interface [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G)) :: scrap ! An array to pass to calculate_density_second_derivs() that is
                  ! set there but will be ignored, it is used simultaneously with four different
                  ! inconsistent units of [R S-1 C-1 ~> kg m-3 degC-1 ppt-1], [R S-2 ~> kg m-3 ppt-2],
                  ! [T2 S-1 L-2 ~> kg m-3 ppt-1 Pa-1] and [T2 C-1 L-2 ~> kg m-3 degC-1 Pa-1].
  real, dimension(SZIB_(G)) :: &
    drho_dT_u, &  ! The derivative of density with temperature at u points [R C-1 ~> kg m-3 degC-1].
    drho_dS_u     ! The derivative of density with salinity at u points [R S-1 ~> kg m-3 ppt-1].
  real, dimension(SZI_(G)) :: &
    drho_dT_v, &  ! The derivative of density with temperature at v points [R C-1 ~> kg m-3 degC-1].
    drho_dS_v, &  ! The derivative of density with salinity at v points [R S-1 ~> kg m-3 ppt-1].
    drho_dT_dT_h, & ! The second derivative of density with temperature at h points [R C-2 ~> kg m-3 degC-2]
    drho_dT_dT_hr ! The second derivative of density with temperature at h (+1) points [R C-2 ~> kg m-3 degC-2]
  real, dimension(SZIB_(G)) :: &
    T_u, &        ! Temperature on the interface at the u-point [C ~> degC].
    S_u, &        ! Salinity on the interface at the u-point [S ~> ppt].
    GxSpV_u, &    ! Gravitiational acceleration times the specific volume at an interface
                  ! at the u-points [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1]
    pres_u        ! Pressure on the interface at the u-point [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G)) :: &
    T_v, &        ! Temperature on the interface at the v-point [C ~> degC].
    S_v, &        ! Salinity on the interface at the v-point [S ~> ppt].
    GxSpV_v, &    ! Gravitiational acceleration times the specific volume at an interface
                  ! at the v-points [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1]
    pres_v        ! Pressure on the interface at the v-point [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G)) :: &
    T_h, &        ! Temperature on the interface at the h-point [C ~> degC].
    S_h, &        ! Salinity on the interface at the h-point [S ~> ppt]
    pres_h, &     ! Pressure on the interface at the h-point [R L2 T-2 ~> Pa].
    T_hr, &       ! Temperature on the interface at the h (+1) point [C ~> degC].
    S_hr, &       ! Salinity on the interface at the h (+1) point [S ~> ppt]
    pres_hr       ! Pressure on the interface at the h (+1) point [R L2 T-2 ~> Pa].
  real :: drdiA, drdiB  ! Along layer zonal potential density  gradients in the layers above (A)
                        ! and below (B) the interface times the grid spacing [R ~> kg m-3].
  real :: drdjA, drdjB  ! Along layer meridional potential density  gradients in the layers above (A)
                        ! and below (B) the interface times the grid spacing [R ~> kg m-3].
  real :: drdkL, drdkR  ! Vertical density differences across an interface [R ~> kg m-3].
  real :: hg2A, hg2B    ! Squares of geometric mean thicknesses [H2 ~> m2 or kg2 m-4].
  real :: hg2L, hg2R    ! Squares of geometric mean thicknesses [H2 ~> m2 or kg2 m-4].
  real :: haA, haB, haL, haR  ! Arithmetic mean thicknesses [H ~> m or kg m-2].
  real :: dzaL, dzaR    ! Temporary thicknesses in eta units [Z ~> m].
  real :: wtA, wtB      ! Unnormalized weights of the slopes above and below [H3 ~> m3 or kg3 m-6]
  real :: wtL, wtR      ! Unnormalized weights of the slopes to the left and right [H3 Z ~> m4 or kg3 m-5]
  real :: drdx, drdy    ! Zonal and meridional density gradients [R L-1 ~> kg m-4].
  real :: drdz          ! Vertical density gradient [R Z-1 ~> kg m-4].
  real :: slope         ! The slope of density surfaces, calculated in a way
                        ! that is always between -1 and 1. [Z L-1 ~> nondim]
  real :: mag_grad2     ! The squared magnitude of the 3-d density gradient [R2 Z-2 ~> kg2 m-8].
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: h_neglect2    ! h_neglect^2 [H2 ~> m2 or kg2 m-4].
  real :: dz_neglect    ! A change in interface heights that is so small it is usually lost
                        ! in roundoff and can be neglected [Z ~> m].
  logical :: use_EOS    ! If true, density is calculated from T & S using an equation of state.
  real :: G_Rho0        ! The gravitational acceleration divided by density [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1]

  logical :: present_N2_u, present_N2_v
  logical :: local_open_u_BC, local_open_v_BC ! True if u- or v-face OBCs exist anywhere in the global domain.
  logical :: OBC_friendly  ! If true, open boundary conditions are in use and only interior data should
                        ! be used to calculate N2 at OBC faces.
  integer, dimension(2) :: EOSdom_u  ! The shifted I-computational domain to use for equation of
                                     ! state calculations at u-points.
  integer, dimension(2) :: EOSdom_v  ! The shifted i-computational domain to use for equation of
                                     ! state calculations at v-points.
  integer, dimension(2) :: EOSdom_h1 ! The shifted i-computational domain to use for equation of
                                     ! state calculations at h points with 1 extra halo point
  integer :: is, ie, js, je, nz, IsdB
  integer :: i, j, k

  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
    EOSdom_h1(:) = EOS_domain(G%HI, halo=halo+1)
  else
    is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
    EOSdom_h1(:) = EOS_domain(G%HI, halo=1)
  endif
  EOSdom_u(1) = is-1 - (G%IsdB-1) ; EOSdom_u(2) = ie - (G%IsdB-1)
  EOSdom_v(:) = EOS_domain(G%HI, halo=halo)

  nz = GV%ke ; IsdB = G%IsdB


  h_neglect = GV%H_subroundoff ; h_neglect2 = h_neglect**2
  dz_neglect = GV%dZ_subroundoff

  local_open_u_BC = .false.
  local_open_v_BC = .false.
  OBC_friendly = .false.
  if (present(OBC)) then ; if (associated(OBC)) then
    local_open_u_BC = OBC%open_u_BCs_exist_globally
    local_open_v_BC = OBC%open_v_BCs_exist_globally
    if (present(OBC_N2)) OBC_friendly = OBC_N2
  endif ; endif

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
  if (present(dzu)) then
    do j=js,je ; do I=is-1,ie
      dzu(I,j,1) = 0.
      dzu(I,j,nz+1) = 0.
    enddo ; enddo
  endif
  if (present(dzv)) then
    do J=js-1,je ; do i=is,ie
      dzv(i,J,1) = 0.
      dzv(i,J,nz+1) = 0.
    enddo ; enddo
  endif
  if (present(dzSxN)) then
    do j=js,je ; do I=is-1,ie
      dzSxN(I,j,1) = 0.
      dzSxN(I,j,nz+1) = 0.
    enddo ; enddo
  endif
  if (present(dzSyN)) then
    do J=js-1,je ; do i=is,ie
      dzSyN(i,J,1) = 0.
      dzSyN(i,J,nz+1) = 0.
    enddo ; enddo
  endif

  if (use_EOS) then
    if (present(halo)) then
      call vert_fill_TS(h, tv%T, tv%S, dt_kappa_smooth, T, S, G, GV, US, halo+1)
    else
      call vert_fill_TS(h, tv%T, tv%S, dt_kappa_smooth, T, S, G, GV, US, 1)
    endif
  endif

  if ((use_EOS .and. allocated(tv%SpV_avg) .and. (tv%valid_SpV_halo < 1)) .and. &
      (present_N2_u .or. present(dzSxN) .or. present_N2_v .or. present(dzSyN))) then
    if (tv%valid_SpV_halo < 0) then
      call MOM_error(FATAL, "calc_isoneutral_slopes called in fully non-Boussinesq mode "//&
                            "with invalid values of SpV_avg.")
    else
      call MOM_error(FATAL, "calc_isoneutral_slopes called in fully non-Boussinesq mode "//&
                            "with insufficiently large SpV_avg halos of width 0 but 1 is needed.")
    endif
  endif

  ! Find the maximum and minimum permitted streamfunction.
  if (associated(tv%p_surf)) then
    !$OMP parallel do default(shared)
    do j=js-1,je+1 ; do i=is-1,ie+1
      pres(i,j,1) = tv%p_surf(i,j)
    enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do j=js-1,je+1 ; do i=is-1,ie+1
      pres(i,j,1) = 0.0
    enddo ; enddo
  endif
  !$OMP parallel do default(shared)
  do j=js-1,je+1
    do k=1,nz ; do i=is-1,ie+1
      pres(i,j,K+1) = pres(i,j,K) + GV%g_Earth * GV%H_to_RZ * h(i,j,k)
    enddo ; enddo
  enddo

  do I=is-1,ie
    GxSpV_u(I) = G_Rho0  ! This will be changed if both use_EOS and allocated(tv%SpV_avg) are true
  enddo
  !$OMP parallel do default(none) shared(nz,is,ie,js,je,IsdB,use_EOS,G,GV,US,pres,T,S,tv,h,e, &
  !$OMP                                  h_neglect,dz_neglect,h_neglect2, &
  !$OMP                                  present_N2_u,G_Rho0,N2_u,slope_x,dzSxN,EOSdom_u,EOSdom_h1, &
  !$OMP                                  local_open_u_BC,dzu,OBC,use_stanley,OBC_friendly) &
  !$OMP                          private(drdiA,drdiB,drdkL,drdkR,pres_u,T_u,S_u,      &
  !$OMP                                  drho_dT_u,drho_dS_u,hg2A,hg2B,hg2L,hg2R,haA, &
  !$OMP                                  drho_dT_dT_h,scrap,pres_h,T_h,S_h, &
  !$OMP                                  haB,haL,haR,dzaL,dzaR,wtA,wtB,wtL,wtR,drdz,  &
  !$OMP                                  drdx,mag_grad2,slope) &
  !$OMP                          firstprivate(GxSpV_u)
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
      if (OBC_friendly) then
        if (OBC%u_E_OBCs_on_PE .and. (j>=OBC%js_u_E_obc) .and. (j<=OBC%je_u_E_obc)) then
          do I = max(is-1, OBC%Is_u_E_obc), min(ie, OBC%Ie_u_E_obc)
            if (OBC%segnum_u(I,j) > 0) then !  OBC_DIRECTION_E
              pres_u(I) = pres(i,j,K)
              T_u(I) = 0.5*(T(i,j,k) + T(i,j,k-1))
              S_u(I) = 0.5*(S(i,j,k) + S(i,j,k-1))
            endif
          enddo
        endif
        if (OBC%u_W_OBCs_on_PE .and. (j>=OBC%js_u_W_obc) .and. (j<=OBC%je_u_W_obc)) then
          do I = max(is-1, OBC%Is_u_W_obc), min(ie, OBC%Ie_u_W_obc)
            if (OBC%segnum_u(I,j) < 0) then !  OBC_DIRECTION_W
              pres_u(I) = pres(i+1,j,K)
              T_u(I) = 0.5*(T(i+1,j,k) + T(i+1,j,k-1))
              S_u(I) = 0.5*(S(i+1,j,k) + S(i+1,j,k-1))
            endif
          enddo
        endif
      endif
      call calculate_density_derivs(T_u, S_u, pres_u, drho_dT_u, drho_dS_u, &
                                    tv%eqn_of_state, EOSdom_u)
      if (present_N2_u .or. (present(dzSxN))) then
        if (allocated(tv%SpV_avg)) then
          do I=is-1,ie
            GxSpV_u(I) = GV%g_Earth *  0.25* ((tv%SpV_avg(i,j,k) + tv%SpV_avg(i+1,j,k)) + &
                                              (tv%SpV_avg(i,j,k-1) + tv%SpV_avg(i+1,j,k-1)))
          enddo
        endif
      endif
    endif

    if (use_stanley) then
      do i=is-1,ie+1
        pres_h(i) = pres(i,j,K)
        T_h(i) = 0.5*(T(i,j,k) + T(i,j,k-1))
        S_h(i) = 0.5*(S(i,j,k) + S(i,j,k-1))
      enddo
      ! The second line below would correspond to arguments
      !            drho_dS_dS, drho_dS_dT, drho_dT_dT, drho_dS_dP, drho_dT_dP, &
      call calculate_density_second_derivs(T_h, S_h, pres_h, &
                   scrap, scrap, drho_dT_dT_h, scrap, scrap, &
                   tv%eqn_of_state, dom=EOSdom_h1)
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
      if (use_stanley) then
        ! Correction to the horizontal density gradient due to nonlinearity in
        ! the EOS rectifying SGS temperature anomalies
        drdiA = drdiA + 0.5 * ((drho_dT_dT_h(i+1) * tv%varT(i+1,j,k-1)) - &
                              (drho_dT_dT_h(i) * tv%varT(i,j,k-1)) )
        drdiB = drdiB + 0.5 * ((drho_dT_dT_h(i+1) * tv%varT(i+1,j,k)) - &
                              (drho_dT_dT_h(i) * tv%varT(i,j,k)) )
      endif

      hg2A = h(i,j,k-1)*h(i+1,j,k-1) + h_neglect2
      hg2B = h(i,j,k)*h(i+1,j,k) + h_neglect2
      hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
      hg2R = h(i+1,j,k-1)*h(i+1,j,k) + h_neglect2
      haA = 0.5*(h(i,j,k-1) + h(i+1,j,k-1)) + h_neglect
      haB = 0.5*(h(i,j,k) + h(i+1,j,k)) + h_neglect
      haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
      haR = 0.5*(h(i+1,j,k-1) + h(i+1,j,k)) + h_neglect
      if (GV%Boussinesq) then
        dzaL = haL * GV%H_to_Z ; dzaR = haR * GV%H_to_Z
      else
        dzaL = 0.5*(e(i,j,K-1) - e(i,j,K+1)) + dz_neglect
        dzaR = 0.5*(e(i+1,j,K-1) - e(i+1,j,K+1)) + dz_neglect
      endif
      if (present(dzu)) dzu(I,j,K) = 0.5*( dzaL + dzaR )
      ! Use the harmonic mean thicknesses to weight the horizontal gradients.
      ! These unnormalized weights have been rearranged to minimize divisions.
      wtA = hg2A*haB ; wtB = hg2B*haA
      wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

      drdz = ((wtL * drdkL) + (wtR * drdkR)) / ((dzaL*wtL) + (dzaR*wtR))
      ! The expression for drdz above is mathematically equivalent to:
      !   drdz = ((hg2L/haL) * drdkL/dzaL + (hg2R/haR) * drdkR/dzaR) / &
      !          ((hg2L/haL) + (hg2R/haR))
      ! which is an estimate of the gradient of density across geopotentials.
      if (present_N2_u) then
        if (OBC_friendly) then ; if (OBC%segnum_u(I,j) /= 0) then
          if (OBC%segnum_u(I,j) > 0) then !  OBC_DIRECTION_E
            drdz = drdkL / dzaL  ! Note that drdz is not used for slopes at OBC faces.
            if (use_EOS .and. allocated(tv%SpV_avg)) &
              GxSpV_u(I) = GV%g_Earth * 0.5 * (tv%SpV_avg(i,j,k) + tv%SpV_avg(i,j,k-1))
          elseif (OBC%segnum_u(I,j) < 0) then !  OBC_DIRECTION_W
            drdz = drdkR / dzaR
            if (use_EOS .and. allocated(tv%SpV_avg)) &
              GxSpV_u(I) = GV%g_Earth * 0.5 * (tv%SpV_avg(i+1,j,k) + tv%SpV_avg(i+1,j,k-1))
          endif
        endif ; endif

        N2_u(I,j,K) = GxSpV_u(I) * drdz * G%mask2dCu(I,j) ! Square of buoyancy freq. [L2 Z-2 T-2 ~> s-2]
      endif

      if (use_EOS) then
        drdx = ((wtA * drdiA + wtB * drdiB) / (wtA + wtB) - &
                drdz * (e(i,j,K)-e(i+1,j,K))) * G%IdxCu(I,j)

        ! This estimate of slope is accurate for small slopes, but bounded
        ! to be between -1 and 1.
        mag_grad2 = (US%Z_to_L*drdx)**2 + drdz**2
        if (mag_grad2 > 0.0) then
          slope = drdx / sqrt(mag_grad2)
        else ! Just in case mag_grad2 = 0 ever.
          slope = 0.0
        endif
      else ! With .not.use_EOS, the layers are constant density.
        slope = (e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j)
      endif

      if (local_open_u_BC) then
        if (OBC%segnum_u(I,j) /= 0) then
          if (OBC%segment(abs(OBC%segnum_u(I,j)))%open) then
            slope = 0.
            ! This and/or the masking code below is to make slopes match inside
            ! land mask. Might not be necessary except for DEBUG output.
!           if (OBC%segnum_u(I,j) > 0) then !  OBC_DIRECTION_E
!             slope_x(I+1,j,K) = 0.
!           else
!             slope_x(I-1,j,K) = 0.
!           endif
          endif
        endif
        slope = slope * max(G%mask2dT(i,j), G%mask2dT(i+1,j))
      endif

      slope_x(I,j,K) = slope
      if (present(dzSxN)) &
        dzSxN(I,j,K) = sqrt( GxSpV_u(I) * max(0., (wtL * ( dzaL * drdkL )) &
                                                + (wtR * ( dzaR * drdkR ))) / (wtL + wtR) ) & ! dz * N
                       * abs(slope) * G%mask2dCu(I,j) ! x-direction contribution to S^2

    enddo ! I
  enddo ; enddo ! end of j-loop

  do i=is,ie
    GxSpV_v(i) = G_Rho0  !This will be changed if both use_EOS and allocated(tv%SpV_avg) are true
  enddo
  ! Calculate the meridional isopycnal slope.
  !$OMP parallel do default(none) shared(nz,is,ie,js,je,IsdB,use_EOS,G,GV,US,pres,T,S,tv, &
  !$OMP                                  h,h_neglect,e,dz_neglect, &
  !$OMP                                  h_neglect2,present_N2_v,G_Rho0,N2_v,slope_y,dzSyN,EOSdom_v, &
  !$OMP                                  dzv,local_open_v_BC,OBC,use_stanley,OBC_friendly) &
  !$OMP                          private(drdjA,drdjB,drdkL,drdkR,pres_v,T_v,S_v,      &
  !$OMP                                  drho_dT_v,drho_dS_v,hg2A,hg2B,hg2L,hg2R,haA, &
  !$OMP                                  drho_dT_dT_h,scrap,pres_h,T_h,S_h, &
  !$OMP                                  drho_dT_dT_hr,pres_hr,T_hr,S_hr,             &
  !$OMP                                  haB,haL,haR,dzaL,dzaR,wtA,wtB,wtL,wtR,drdz,  &
  !$OMP                                  drdy,mag_grad2,slope) &
  !$OMP                          firstprivate(GxSpV_v)
  do J=js-1,je ; do K=nz,2,-1
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
      if (OBC_friendly) then
        if (OBC%v_N_OBCs_on_PE .and. (J>=OBC%Js_v_N_obc) .and. (J<=OBC%Je_v_N_obc)) then
          do i = max(is, OBC%is_v_N_obc), min(ie, OBC%ie_v_N_obc)
            if (OBC%segnum_v(i,J) > 0) then !  OBC_DIRECTION_N
              pres_v(i) = pres(i,j,K)
              T_v(i) = 0.5*(T(i,j,k) + T(i,j,k-1))
              S_v(i) = 0.5*(S(i,j,k) + S(i,j,k-1))
            endif
          enddo
        endif
        if (OBC%v_S_OBCs_on_PE .and. (J>=OBC%Js_v_S_obc) .and. (J<=OBC%Je_v_S_obc)) then
          do i = max(is, OBC%is_v_S_obc), min(ie, OBC%ie_v_S_obc)
            if (OBC%segnum_v(i,J) < 0) then !  OBC_DIRECTION_S
              pres_v(i) = pres(i,j+1,K)
              T_v(i) = 0.5*(T(i,j+1,k) + T(i,j+1,k-1))
              S_v(i) = 0.5*(S(i,j+1,k) + S(i,j+1,k-1))
            endif
          enddo
        endif
      endif
      call calculate_density_derivs(T_v, S_v, pres_v, drho_dT_v, drho_dS_v, &
                                    tv%eqn_of_state, EOSdom_v)

      if ((present_N2_v) .or. (present(dzSyN))) then
        if (allocated(tv%SpV_avg)) then
          do i=is,ie
            GxSpV_v(i) = GV%g_Earth *  0.25* ((tv%SpV_avg(i,j,k) + tv%SpV_avg(i,j+1,k)) + &
                                              (tv%SpV_avg(i,j,k-1) + tv%SpV_avg(i,j+1,k-1)))
          enddo
        endif
      endif
    endif

    if (use_stanley) then
      do i=is,ie
        pres_h(i) = pres(i,j,K)
        T_h(i) = 0.5*(T(i,j,k) + T(i,j,k-1))
        S_h(i) = 0.5*(S(i,j,k) + S(i,j,k-1))

        pres_hr(i) = pres(i,j+1,K)
        T_hr(i) = 0.5*(T(i,j+1,k) + T(i,j+1,k-1))
        S_hr(i) = 0.5*(S(i,j+1,k) + S(i,j+1,k-1))
      enddo
      ! The second line below would correspond to arguments
      !            drho_dS_dS, drho_dS_dT, drho_dT_dT, drho_dS_dP, drho_dT_dP, &
      call calculate_density_second_derivs(T_h, S_h, pres_h, &
                   scrap, scrap, drho_dT_dT_h, scrap, scrap, &
                   tv%eqn_of_state, dom=EOSdom_v)
      call calculate_density_second_derivs(T_hr, S_hr, pres_hr, &
                   scrap, scrap, drho_dT_dT_hr, scrap, scrap, &
                   tv%eqn_of_state, dom=EOSdom_v)
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
      if (use_stanley) then
        ! Correction to the horizontal density gradient due to nonlinearity in
        ! the EOS rectifying SGS temperature anomalies
        drdjA = drdjA + 0.5 * ((drho_dT_dT_hr(i) * tv%varT(i,j+1,k-1)) - &
                              (drho_dT_dT_h(i) * tv%varT(i,j,k-1)) )
        drdjB = drdjB + 0.5 * ((drho_dT_dT_hr(i) * tv%varT(i,j+1,k)) - &
                              (drho_dT_dT_h(i) * tv%varT(i,j,k)) )
      endif

      hg2A = h(i,j,k-1)*h(i,j+1,k-1) + h_neglect2
      hg2B = h(i,j,k)*h(i,j+1,k) + h_neglect2
      hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
      hg2R = h(i,j+1,k-1)*h(i,j+1,k) + h_neglect2
      haA = 0.5*(h(i,j,k-1) + h(i,j+1,k-1)) + h_neglect
      haB = 0.5*(h(i,j,k) + h(i,j+1,k)) + h_neglect
      haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
      haR = 0.5*(h(i,j+1,k-1) + h(i,j+1,k)) + h_neglect
      if (GV%Boussinesq) then
        dzaL = haL * GV%H_to_Z ; dzaR = haR * GV%H_to_Z
      else
        dzaL = 0.5*(e(i,j,K-1) - e(i,j,K+1)) + dz_neglect
        dzaR = 0.5*(e(i,j+1,K-1) - e(i,j+1,K+1)) + dz_neglect
      endif
      if (present(dzv)) dzv(i,J,K) = 0.5*( dzaL + dzaR )
      ! Use the harmonic mean thicknesses to weight the horizontal gradients.
      ! These unnormalized weights have been rearranged to minimize divisions.
      wtA = hg2A*haB ; wtB = hg2B*haA
      wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

      drdz = ((wtL * drdkL) + (wtR * drdkR)) / ((dzaL*wtL) + (dzaR*wtR))
      ! The expression for drdz above is mathematically equivalent to:
      !   drdz = ((hg2L/haL) * drdkL/dzaL + (hg2R/haR) * drdkR/dzaR) / &
      !          ((hg2L/haL) + (hg2R/haR))
      ! which is an estimate of the gradient of density across geopotentials.
      if (present_N2_v) then
        if (OBC_friendly) then ; if (OBC%segnum_v(i,J) /= 0) then
          if (OBC%segnum_v(i,J) > 0) then !  OBC_DIRECTION_N
            drdz = drdkL / dzaL  ! Note that drdz is not used for slopes at OBC faces.
            if (use_EOS .and. allocated(tv%SpV_avg)) &
              GxSpV_v(i) = GV%g_Earth * 0.5 * (tv%SpV_avg(i,j,k) + tv%SpV_avg(i,j,k-1))
          elseif (OBC%segnum_v(i,J) < 0) then !  OBC_DIRECTION_S
            drdz = drdkL / dzaL
            if (use_EOS .and. allocated(tv%SpV_avg)) &
              GxSpV_v(i) = GV%g_Earth * 0.5 * (tv%SpV_avg(i,j+1,k) + tv%SpV_avg(i,j+1,k-1))
          endif
        endif ; endif

        N2_v(i,J,K) = GxSpV_v(i) * drdz * G%mask2dCv(i,J) ! Square of buoyancy freq. [L2 Z-2 T-2 ~> s-2]
      endif

      if (use_EOS) then
        drdy = ((wtA * drdjA + wtB * drdjB) / (wtA + wtB) - &
                drdz * (e(i,j,K)-e(i,j+1,K))) * G%IdyCv(i,J)

        ! This estimate of slope is accurate for small slopes, but bounded
        ! to be between -1 and 1.
        mag_grad2 = (US%Z_to_L*drdy)**2 + drdz**2
        if (mag_grad2 > 0.0) then
          slope = drdy / sqrt(mag_grad2)
        else ! Just in case mag_grad2 = 0 ever.
          slope = 0.0
        endif
      else ! With .not.use_EOS, the layers are constant density.
        slope = (e(i,j+1,K)-e(i,j,K)) * G%IdyCv(i,J)
      endif

      if (local_open_v_BC) then
        if (OBC%segnum_v(i,J) /= 0) then
          if (OBC%segment(abs(OBC%segnum_v(i,J)))%open) then
            slope = 0.
            ! This and/or the masking code below is to make slopes match inside
            ! land mask. Might not be necessary except for DEBUG output.
!           if (OBC%segnum_v(i,J)) > 0) then ! OBC_DIRECTION_N
!             slope_y(i,J+1,K) = 0.
!           else
!             slope_y(i,J-1,K) = 0.
!           endif
          endif
        endif
        slope = slope * max(G%mask2dT(i,j), G%mask2dT(i,j+1))
      endif
      slope_y(i,J,K) = slope
      if (present(dzSyN)) &
        dzSyN(i,J,K) = sqrt( GxSpV_v(i) * max(0., (wtL * ( dzaL * drdkL )) &
                                                + (wtR * ( dzaR * drdkR ))) / (wtL + wtR) ) & ! dz * N
                        * abs(slope) * G%mask2dCv(i,J) ! x-direction contribution to S^2

    enddo ! i
  enddo ; enddo ! end of j-loop

end subroutine calc_isoneutral_slopes

!> Returns tracer arrays (nominally T and S) with massless layers filled with
!! sensible values, by diffusing vertically with a small but constant diffusivity.
subroutine vert_fill_TS(h, T_in, S_in, kappa_dt, T_f, S_f, G, GV, US, halo_here, larger_h_denom)
  type(ocean_grid_type),                     intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                     intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: T_in !< Input temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: S_in !< Input salinity [S ~> ppt]
  real,                                      intent(in)  :: kappa_dt !< A vertical diffusivity to use for smoothing
                                                                 !! times a smoothing timescale [H Z ~> m2 or kg m-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T_f  !< Filled temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S_f  !< Filled salinity [S ~> ppt]
  integer,                         optional, intent(in)  :: halo_here !< Number of halo points to work on,
                                                                 !! 0 by default
  logical,                         optional, intent(in)  :: larger_h_denom !< Present and true, add a large
                                                                 !! enough minimal thickness in the denominator of
                                                                 !! the flux calculations so that the fluxes are
                                                                 !! never so large as eliminate the transmission
                                                                 !! of information across groups of massless layers.
  ! Local variables
  real :: ent(SZI_(G),SZK_(GV)+1)  ! The diffusive entrainment (kappa*dt)/dz
                                   ! between layers in a timestep [H ~> m or kg m-2].
  real :: b1(SZI_(G))              ! A variable used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1]
  real :: d1(SZI_(G))              ! A variable used by the tridiagonal solver [nondim], d1 = 1 - c1.
  real :: c1(SZI_(G),SZK_(GV))     ! A variable used by the tridiagonal solver [nondim].
  real :: kap_dt_x2                ! The 2*kappa_dt converted to H units [H2 ~> m2 or kg2 m-4].
  real :: h_neglect                ! A negligible thickness [H ~> m or kg m-2], to allow for zero thicknesses.
  real :: h0                       ! A negligible thickness to allow for zero thickness layers without
                                   ! completely decoupling groups of layers [H ~> m or kg m-2].
                                   ! Often 0 < h_neglect << h0.
  real :: h_tr                     ! h_tr is h at tracer points with a tiny thickness
                                   ! added to ensure positive definiteness [H ~> m or kg m-2].
  integer :: i, j, k, is, ie, js, je, nz, halo

  halo=0 ; if (present(halo_here)) halo = max(halo_here,0)

  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo ; nz = GV%ke

  h_neglect = GV%H_subroundoff
  ! The use of the fixed rescaling factor in the next line avoids an extra call to thickness_to_dz()
  ! and the use of an extra 3-d array of vertical distnaces across layers (dz).  This would be more
  ! physically consistent, but it would also be more expensive, and given that this routine applies
  ! a small (but arbitrary) amount of mixing to clean up the properties of nearly massless layers,
  ! the added expense is hard to justify.
  kap_dt_x2 = (2.0*kappa_dt) * (US%Z_to_m*GV%m_to_H) ! Usually the latter term is GV%Z_to_H.
  h0 = h_neglect
  if (present(larger_h_denom)) then
    if (larger_h_denom) h0 = 1.0e-16*sqrt(0.5*kap_dt_x2)
  endif

  if (kap_dt_x2 <= 0.0) then
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=js,je ; do i=is,ie
      T_f(i,j,k) = T_in(i,j,k) ; S_f(i,j,k) = S_in(i,j,k)
    enddo ; enddo ; enddo
  else
    !$OMP parallel do default(shared) private(ent,b1,d1,c1,h_tr)
    do j=js,je
      do i=is,ie
        ent(i,2) = kap_dt_x2 / ((h(i,j,1)+h(i,j,2)) + h0)
        h_tr = h(i,j,1) + h_neglect
        b1(i) = 1.0 / (h_tr + ent(i,2))
        d1(i) = b1(i) * h_tr
        T_f(i,j,1) = (b1(i)*h_tr)*T_in(i,j,1)
        S_f(i,j,1) = (b1(i)*h_tr)*S_in(i,j,1)
      enddo
      do k=2,nz-1 ; do i=is,ie
        ent(i,K+1) = kap_dt_x2 / ((h(i,j,k)+h(i,j,k+1)) + h0)
        h_tr = h(i,j,k) + h_neglect
        c1(i,k) = ent(i,K) * b1(i)
        b1(i) = 1.0 / ((h_tr + d1(i)*ent(i,K)) + ent(i,K+1))
        d1(i) = b1(i) * (h_tr + d1(i)*ent(i,K))
        T_f(i,j,k) = b1(i) * (h_tr*T_in(i,j,k) + ent(i,K)*T_f(i,j,k-1))
        S_f(i,j,k) = b1(i) * (h_tr*S_in(i,j,k) + ent(i,K)*S_f(i,j,k-1))
      enddo ; enddo
      do i=is,ie
        c1(i,nz) = ent(i,nz) * b1(i)
        h_tr = h(i,j,nz) + h_neglect
        b1(i) = 1.0 / (h_tr + d1(i)*ent(i,nz))
        T_f(i,j,nz) = b1(i) * (h_tr*T_in(i,j,nz) + ent(i,nz)*T_f(i,j,nz-1))
        S_f(i,j,nz) = b1(i) * (h_tr*S_in(i,j,nz) + ent(i,nz)*S_f(i,j,nz-1))
      enddo
      do k=nz-1,1,-1 ; do i=is,ie
        T_f(i,j,k) = T_f(i,j,k) + c1(i,k+1)*T_f(i,j,k+1)
        S_f(i,j,k) = S_f(i,j,k) + c1(i,k+1)*S_f(i,j,k+1)
      enddo ; enddo
    enddo
  endif

end subroutine vert_fill_TS

end module MOM_isopycnal_slopes
