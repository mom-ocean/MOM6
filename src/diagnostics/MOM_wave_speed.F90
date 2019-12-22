!> Routines for calculating baroclinic wave speeds
module MOM_wave_speed

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : log_version
use MOM_grid, only : ocean_grid_type
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>

public wave_speed, wave_speeds, wave_speed_init, wave_speed_set_param

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for MOM_wave_speed
type, public :: wave_speed_CS ; private
  logical :: use_ebt_mode = .false.    !< If true, calculate the equivalent barotropic wave speed instead
                                       !! of the first baroclinic wave speed.
                                       !! This parameter controls the default behavior of wave_speed() which
                                       !! can be overridden by optional arguments.
  real :: mono_N2_column_fraction = 0. !< The lower fraction of water column over which N2 is limited as
                                       !! monotonic for the purposes of calculating the equivalent barotropic
                                       !! wave speed. This parameter controls the default behavior of
                                       !! wave_speed() which can be overridden by optional arguments.
  real :: mono_N2_depth = -1.          !< The depth below which N2 is limited as monotonic for the purposes of
                                       !! calculating the equivalent barotropic wave speed [Z ~> m].
                                       !! This parameter controls the default behavior of wave_speed() which
                                       !! can be overridden by optional arguments.
  type(remapping_CS) :: remapping_CS   !< Used for vertical remapping when calculating equivalent barotropic
                                       !! mode structure.
  logical :: remap_answers_2018 = .true.  !< If true, use the order of arithmetic and expressions that
                                       !! recover the remapping answers from 2018.  If false, use more
                                       !! robust forms of the same remapping expressions.
  type(diag_ctrl), pointer :: diag     !< Diagnostics control structure
end type wave_speed_CS

contains

!> Calculates the wave speed of the first baroclinic mode.
subroutine wave_speed(h, tv, G, GV, US, cg1, CS, full_halos, use_ebt_mode, &
                      mono_N2_column_fraction, mono_N2_depth, modal_structure)
  type(ocean_grid_type),            intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),          intent(in)  :: GV !< Vertical grid structure
  type(unit_scale_type),            intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                                    intent(in)  :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),            intent(in)  :: tv !< Thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: cg1 !< First mode internal wave speed [L T-1 ~> m s-1]
  type(wave_speed_CS),              pointer     :: CS !< Control structure for MOM_wave_speed
  logical, optional,                intent(in)  :: full_halos !< If true, do the calculation
                                          !! over the entire computational domain.
  logical, optional,                intent(in)  :: use_ebt_mode !< If true, use the equivalent
                                          !! barotropic mode instead of the first baroclinic mode.
  real, optional,                   intent(in)  :: mono_N2_column_fraction !< The lower fraction
                                          !! of water column over which N2 is limited as monotonic
                                          !! for the purposes of calculating vertical modal structure.
  real, optional,                   intent(in)  :: mono_N2_depth !< A depth below which N2 is limited as
                                          !! monotonic for the purposes of calculating vertical
                                          !! modal structure [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
        optional,                   intent(out) :: modal_structure !< Normalized model structure [nondim]

  ! Local variables
  real, dimension(SZK_(G)+1) :: &
    dRho_dT, &    ! Partial derivative of density with temperature [R degC-1 ~> kg m-3 degC-1]
    dRho_dS, &    ! Partial derivative of density with salinity [R ppt-1 ~> kg m-3 ppt-1]
    pres, &       ! Interface pressure [Pa]
    T_int, &      ! Temperature interpolated to interfaces [degC]
    S_int, &      ! Salinity interpolated to interfaces [ppt]
    gprime        ! The reduced gravity across each interface [L2 Z-1 T-2 ~> m s-2].
  real, dimension(SZK_(G)) :: &
    Igl, Igu, Igd ! The inverse of the reduced gravity across an interface times
                  ! the thickness of the layer below (Igl) or above (Igu) it.
                  ! Their sum, Igd, is provided for the tridiagonal solver.  [T2 L-2 ~> s2 m-2]
  real, dimension(SZK_(G),SZI_(G)) :: &
    Hf, &         ! Layer thicknesses after very thin layers are combined [Z ~> m]
    Tf, &         ! Layer temperatures after very thin layers are combined [degC]
    Sf, &         ! Layer salinities after very thin layers are combined [ppt]
    Rf            ! Layer densities after very thin layers are combined [R ~> kg m-3]
  real, dimension(SZK_(G)) :: &
    Hc, &         ! A column of layer thicknesses after convective istabilities are removed [Z ~> m]
    Tc, &         ! A column of layer temperatures after convective istabilities are removed [degC]
    Sc, &         ! A column of layer salinites after convective istabilities are removed [ppt]
    Rc, &         ! A column of layer densities after convective istabilities are removed [R ~> kg m-3]
    Hc_H          ! Hc(:) rescaled from Z to thickness units [H ~> m or kg m-2]
  real :: det, ddet, detKm1, detKm2, ddetKm1, ddetKm2
  real :: lam     ! The eigenvalue [T2 L-2 ~> s m-1]
  real :: dlam    ! The change in estimates of the eigenvalue [T2 L-2 ~> s m-1]
  real :: lam0    ! The first guess of the eigenvalue [T2 L-2 ~> s m-1]
  real :: min_h_frac ! [nondim]
  real :: Z_to_Pa  ! A conversion factor from thicknesses (in Z) to pressure (in Pa)
  real, dimension(SZI_(G)) :: &
    htot, hmin, &  ! Thicknesses [Z ~> m]
    H_here, &      ! A thickness [Z ~> m]
    HxT_here, &    ! A layer integrated temperature [degC Z ~> degC m]
    HxS_here, &    ! A layer integrated salinity [ppt Z ~> ppt m]
    HxR_here       ! A layer integrated density [R Z ~> kg m-2]
  real :: speed2_tot ! overestimate of the mode-1 speed squared [L2 T-2 ~> m2 s-2]
  real :: I_Hnew   ! The inverse of a new layer thickness [Z-1 ~> m-1]
  real :: drxh_sum ! The sum of density diffrences across interfaces times thicknesses [R Z ~> kg m-2]
  real :: L2_to_Z2 ! A scaling factor squared from units of lateral distances to depths [Z2 L-2 ~> 1].
  real, parameter :: tol1  = 0.0001, tol2 = 0.001
  real, pointer, dimension(:,:,:) :: T => NULL(), S => NULL()
  real :: g_Rho0  ! G_Earth/Rho0 [L2 T-2 Z-1 R-1 ~> m4 s-2 kg-1].
  real :: c2_scale ! A scaling factor for wave speeds to help control the growth of the determinant
                   ! and its derivative with lam between rows of the Thomas algorithm solver.  The
                   ! exact value should not matter for the final result if it is an even power of 2.
  real :: rescale, I_rescale
  integer :: kf(SZI_(G))
  integer, parameter :: max_itt = 10
  real :: lam_it(max_itt), det_it(max_itt), ddet_it(max_itt)
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  integer :: kc
  integer :: i, j, k, k2, itt, is, ie, js, je, nz
  real :: hw, sum_hc
  real :: gp      ! A limited local copy of gprime [L2 Z-1 T-2 ~> m s-2]
  real :: N2min   ! A minimum buoyancy frequency [T-2 ~> s-2]
  logical :: l_use_ebt_mode, calc_modal_structure
  real :: l_mono_N2_column_fraction, l_mono_N2_depth
  real :: mode_struct(SZK_(G)), ms_min, ms_max, ms_sq

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_wave_speed: "// &
           "Module must be initialized before it is used.")
  if (present(full_halos)) then ; if (full_halos) then
    is = G%isd ; ie = G%ied ; js = G%jsd ; je = G%jed
  endif ; endif

  L2_to_Z2 = US%L_to_Z**2

  l_use_ebt_mode = CS%use_ebt_mode
  if (present(use_ebt_mode)) l_use_ebt_mode = use_ebt_mode
  l_mono_N2_column_fraction = CS%mono_N2_column_fraction
  if (present(mono_N2_column_fraction)) l_mono_N2_column_fraction = mono_N2_column_fraction
  l_mono_N2_depth = CS%mono_N2_depth
  if (present(mono_N2_depth)) l_mono_N2_depth = mono_N2_depth
  calc_modal_structure = l_use_ebt_mode
  if (present(modal_structure)) calc_modal_structure = .true.
  if (calc_modal_structure) then
    do k=1,nz; do j=js,je; do i=is,ie
      modal_structure(i,j,k) = 0.0
    enddo ; enddo ; enddo
  endif

  S => tv%S ; T => tv%T
  g_Rho0 = GV%g_Earth / GV%Rho0
  Z_to_Pa = GV%Z_to_H * GV%H_to_Pa
  use_EOS = associated(tv%eqn_of_state)

  rescale = 1024.0**4 ; I_rescale = 1.0/rescale
  ! The following two lines give identical results:
  ! c2_scale = 16.0 * US%m_s_to_L_T**2
  c2_scale = US%m_s_to_L_T**2

  min_h_frac = tol1 / real(nz)
!$OMP parallel do default(none) shared(is,ie,js,je,nz,h,G,GV,US,min_h_frac,use_EOS,T,S,tv,&
!$OMP                                  calc_modal_structure,l_use_ebt_mode,modal_structure, &
!$OMP                                  l_mono_N2_column_fraction,l_mono_N2_depth,CS,   &
!$OMP                                  Z_to_Pa,cg1,g_Rho0,rescale,I_rescale,L2_to_Z2,c2_scale)  &
!$OMP                          private(htot,hmin,kf,H_here,HxT_here,HxS_here,HxR_here, &
!$OMP                                  Hf,Tf,Sf,Rf,pres,T_int,S_int,drho_dT,           &
!$OMP                                  drho_dS,drxh_sum,kc,Hc,Hc_H,Tc,Sc,I_Hnew,gprime,&
!$OMP                                  Rc,speed2_tot,Igl,Igu,Igd,lam0,lam,lam_it,dlam, &
!$OMP                                  mode_struct,sum_hc,N2min,gp,hw,                 &
!$OMP                                  ms_min,ms_max,ms_sq,                            &
!$OMP                                  det,ddet,detKm1,ddetKm1,detKm2,ddetKm2,det_it,ddet_it)
  do j=js,je
    !   First merge very thin layers with the one above (or below if they are
    ! at the top).  This also transposes the row order so that columns can
    ! be worked upon one at a time.
    do i=is,ie ; htot(i) = 0.0 ; enddo
    do k=1,nz ; do i=is,ie ; htot(i) = htot(i) + h(i,j,k)*GV%H_to_Z ; enddo ; enddo

    do i=is,ie
      hmin(i) = htot(i)*min_h_frac ; kf(i) = 1 ; H_here(i) = 0.0
      HxT_here(i) = 0.0 ; HxS_here(i) = 0.0 ; HxR_here(i) = 0.0
    enddo
    if (use_EOS) then
      do k=1,nz ; do i=is,ie
        if ((H_here(i) > hmin(i)) .and. (h(i,j,k)*GV%H_to_Z > hmin(i))) then
          Hf(kf(i),i) = H_here(i)
          Tf(kf(i),i) = HxT_here(i) / H_here(i)
          Sf(kf(i),i) = HxS_here(i) / H_here(i)
          kf(i) = kf(i) + 1

          ! Start a new layer
          H_here(i) = h(i,j,k)*GV%H_to_Z
          HxT_here(i) = (h(i,j,k)*GV%H_to_Z)*T(i,j,k)
          HxS_here(i) = (h(i,j,k)*GV%H_to_Z)*S(i,j,k)
        else
          H_here(i) = H_here(i) + h(i,j,k)*GV%H_to_Z
          HxT_here(i) = HxT_here(i) + (h(i,j,k)*GV%H_to_Z)*T(i,j,k)
          HxS_here(i) = HxS_here(i) + (h(i,j,k)*GV%H_to_Z)*S(i,j,k)
        endif
      enddo ; enddo
      do i=is,ie ; if (H_here(i) > 0.0) then
        Hf(kf(i),i) = H_here(i)
        Tf(kf(i),i) = HxT_here(i) / H_here(i)
        Sf(kf(i),i) = HxS_here(i) / H_here(i)
      endif ; enddo
    else
      do k=1,nz ; do i=is,ie
        if ((H_here(i) > hmin(i)) .and. (h(i,j,k)*GV%H_to_Z > hmin(i))) then
          Hf(kf(i),i) = H_here(i) ; Rf(kf(i),i) = HxR_here(i) / H_here(i)
          kf(i) = kf(i) + 1

          ! Start a new layer
          H_here(i) = h(i,j,k)*GV%H_to_Z
          HxR_here(i) = (h(i,j,k)*GV%H_to_Z)*GV%Rlay(k)
        else
          H_here(i) = H_here(i) + h(i,j,k)*GV%H_to_Z
          HxR_here(i) = HxR_here(i) + (h(i,j,k)*GV%H_to_Z)*GV%Rlay(k)
        endif
      enddo ; enddo
      do i=is,ie ; if (H_here(i) > 0.0) then
        Hf(kf(i),i) = H_here(i) ; Rf(kf(i),i) = HxR_here(i) / H_here(i)
      endif ; enddo
    endif

    ! From this point, we can work on individual columns without causing memory
    ! to have page faults.
    do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
      if (use_EOS) then
        pres(1) = 0.0
        do k=2,kf(i)
          pres(k) = pres(k-1) + Z_to_Pa*Hf(k-1,i)
          T_int(k) = 0.5*(Tf(k,i)+Tf(k-1,i))
          S_int(k) = 0.5*(Sf(k,i)+Sf(k-1,i))
        enddo
        call calculate_density_derivs(T_int, S_int, pres, drho_dT, drho_dS, 2, &
                                      kf(i)-1, tv%eqn_of_state, scale=US%kg_m3_to_R)

        ! Sum the reduced gravities to find out how small a density difference
        ! is negligibly small.
        drxh_sum = 0.0
        do k=2,kf(i)
          drxh_sum = drxh_sum + 0.5*(Hf(k-1,i)+Hf(k,i)) * &
              max(0.0,drho_dT(k)*(Tf(k,i)-Tf(k-1,i)) + &
                      drho_dS(k)*(Sf(k,i)-Sf(k-1,i)))
        enddo
      else
        drxh_sum = 0.0
        do k=2,kf(i)
          drxh_sum = drxh_sum + 0.5*(Hf(k-1,i)+Hf(k,i)) * &
                            max(0.0,Rf(k,i)-Rf(k-1,i))
        enddo
      endif

      if (calc_modal_structure) then
        mode_struct(:) = 0.
      endif

  !   Find gprime across each internal interface, taking care of convective
  ! instabilities by merging layers.
      if (drxh_sum <= 0.0) then
        cg1(i,j) = 0.0
      else
        ! Merge layers to eliminate convective instabilities or exceedingly
        ! small reduced gravities.
        if (use_EOS) then
          kc = 1
          Hc(1) = Hf(1,i) ; Tc(1) = Tf(1,i) ; Sc(1) = Sf(1,i)
          do k=2,kf(i)
            if ((drho_dT(k)*(Tf(k,i)-Tc(kc)) + drho_dS(k)*(Sf(k,i)-Sc(kc))) * &
                (Hc(kc) + Hf(k,i)) < 2.0 * tol2*drxh_sum) then
              ! Merge this layer with the one above and backtrack.
              I_Hnew = 1.0 / (Hc(kc) + Hf(k,i))
              Tc(kc) = (Hc(kc)*Tc(kc) + Hf(k,i)*Tf(k,i)) * I_Hnew
              Sc(kc) = (Hc(kc)*Sc(kc) + Hf(k,i)*Sf(k,i)) * I_Hnew
              Hc(kc) = (Hc(kc) + Hf(k,i))
              ! Backtrack to remove any convective instabilities above...  Note
              ! that the tolerance is a factor of two larger, to avoid limit how
              ! far back we go.
              do k2=kc,2,-1
                if ((drho_dT(k2)*(Tc(k2)-Tc(k2-1)) + drho_dS(k2)*(Sc(k2)-Sc(k2-1))) * &
                    (Hc(k2) + Hc(k2-1)) < tol2*drxh_sum) then
                  ! Merge the two bottommost layers.  At this point kc = k2.
                  I_Hnew = 1.0 / (Hc(kc) + Hc(kc-1))
                  Tc(kc-1) = (Hc(kc)*Tc(kc) + Hc(kc-1)*Tc(kc-1)) * I_Hnew
                  Sc(kc-1) = (Hc(kc)*Sc(kc) + Hc(kc-1)*Sc(kc-1)) * I_Hnew
                  Hc(kc-1) = (Hc(kc) + Hc(kc-1))
                  kc = kc - 1
                else ; exit ; endif
              enddo
            else
              ! Add a new layer to the column.
              kc = kc + 1
              drho_dS(kc) = drho_dS(k) ; drho_dT(kc) = drho_dT(k)
              Tc(kc) = Tf(k,i) ; Sc(kc) = Sf(k,i) ; Hc(kc) = Hf(k,i)
            endif
          enddo
          ! At this point there are kc layers and the gprimes should be positive.
          do k=2,kc ! Revisit this if non-Boussinesq.
            gprime(k) = g_Rho0 * (drho_dT(k)*(Tc(k)-Tc(k-1)) + &
                                  drho_dS(k)*(Sc(k)-Sc(k-1)))
          enddo
        else  ! .not.use_EOS
          ! Do the same with density directly...
          kc = 1
          Hc(1) = Hf(1,i) ; Rc(1) = Rf(1,i)
          do k=2,kf(i)
            if ((Rf(k,i) - Rc(kc)) * (Hc(kc) + Hf(k,i)) < 2.0*tol2*drxh_sum) then
              ! Merge this layer with the one above and backtrack.
              Rc(kc) = (Hc(kc)*Rc(kc) + Hf(k,i)*Rf(k,i)) / (Hc(kc) + Hf(k,i))
              Hc(kc) = (Hc(kc) + Hf(k,i))
              ! Backtrack to remove any convective instabilities above...  Note
              ! that the tolerance is a factor of two larger, to avoid limit how
              ! far back we go.
              do k2=kc,2,-1
                if ((Rc(k2)-Rc(k2-1)) * (Hc(k2)+Hc(k2-1)) < tol2*drxh_sum) then
                  ! Merge the two bottommost layers.  At this point kc = k2.
                  Rc(kc-1) = (Hc(kc)*Rc(kc) + Hc(kc-1)*Rc(kc-1)) / (Hc(kc) + Hc(kc-1))
                  Hc(kc-1) = (Hc(kc) + Hc(kc-1))
                  kc = kc - 1
                else ; exit ; endif
              enddo
            else
              ! Add a new layer to the column.
              kc = kc + 1
              Rc(kc) = Rf(k,i) ; Hc(kc) = Hf(k,i)
            endif
          enddo
          ! At this point there are kc layers and the gprimes should be positive.
          do k=2,kc ! Revisit this if non-Boussinesq.
            gprime(k) = g_Rho0 * (Rc(k)-Rc(k-1))
          enddo
        endif  ! use_EOS

        ! Sum the contributions from all of the interfaces to give an over-estimate
        ! of the first-mode wave speed.  Also populate Igl and Igu which are the
        ! non-leading diagonals of the tridiagonal matrix.
        if (kc >= 2) then
          speed2_tot = 0.0
          if (l_use_ebt_mode) then
            Igu(1) = 0. ! Neumann condition for pressure modes
            sum_hc = Hc(1)
            N2min = L2_to_Z2*gprime(2)/Hc(1)
            do k=2,kc
              hw = 0.5*(Hc(k-1)+Hc(k))
              gp = gprime(K)
              if (l_mono_N2_column_fraction>0. .or. l_mono_N2_depth>=0.) then
                if ( ((G%bathyT(i,j)-sum_hc < l_mono_N2_column_fraction*G%bathyT(i,j)) .or. &
                      ((l_mono_N2_depth >= 0.) .and. (sum_hc > l_mono_N2_depth))) .and. &
                     (L2_to_Z2*gp > N2min*hw) ) then
                  ! Filters out regions where N2 increases with depth but only in a lower fraction
                  ! of the water column or below a certain depth.
                  gp = US%Z_to_L**2 * (N2min*hw)
                else
                  N2min = L2_to_Z2 * gp/hw
                endif
              endif
              Igu(k) = 1.0/(gp*Hc(k))
              Igl(k-1) = 1.0/(gp*Hc(k-1))
              speed2_tot = speed2_tot + gprime(k)*(Hc(k-1)+Hc(k))*0.707
              sum_hc = sum_hc + Hc(k)
            enddo
           !Igl(kc) = 0. ! Neumann condition for pressure modes
            Igl(kc) = 2.*Igu(kc) ! Dirichlet condition for pressure modes
          else ! .not. l_use_ebt_mode
            do K=2,kc
              Igl(K) = 1.0/(gprime(K)*Hc(k)) ; Igu(K) = 1.0/(gprime(K)*Hc(k-1))
              speed2_tot = speed2_tot + gprime(k)*(Hc(k-1)+Hc(k))
            enddo
          endif

          if (calc_modal_structure) then
            mode_struct(1:kc) = 1. ! Uniform flow, first guess
          endif

          ! Overestimate the speed to start with.
          if (calc_modal_structure) then
            lam0 = 0.5 / speed2_tot ; lam = lam0
          else
            lam0 = 1.0 / speed2_tot ; lam = lam0
          endif
          ! Find the determinant and its derivative with lam.
          do itt=1,max_itt
            lam_it(itt) = lam
            if (l_use_ebt_mode) then
              ! This initialization of det,ddet imply Neumann boundary conditions so that first 3 rows
              ! of the matrix are
              !    /   b(1)-lam  igl(1)      0        0     0  ...  \
              !    |  igu(2)    b(2)-lam   igl(2)     0     0  ...  |
              !    |    0        igu(3)   b(3)-lam  igl(3)  0  ...  |
              ! which is consistent if the eigenvalue problem is for horizontal velocity or pressure modes.
             !detKm1 = c2_scale*(Igl(1)-lam) ; ddetKm1 = -1.0*c2_scale
             !det = (Igu(2)+Igl(2)-lam)*detKm1 - (Igu(2)*Igl(1)) ; ddet = (Igu(2)+Igl(2)-lam)*ddetKm1 - detKm1
              detKm1 = 1.0 ; ddetKm1 = 0.0
              det = (Igl(1)-lam) ; ddet = -1.0
              if (kc>1) then
                ! Shift variables and rescale rows to avoid over- or underflow.
                detKm2 = c2_scale*detKm1 ; ddetKm2 = c2_scale*ddetKm1
                detKm1 = c2_scale*det    ; ddetKm1 = c2_scale*ddet
                det = (Igu(2)+Igl(2)-lam)*detKm1 - (Igu(2)*Igl(1))*detKm2
                ddet = (Igu(2)+Igl(2)-lam)*ddetKm1 - (Igu(2)*Igl(1))*ddetKm2 - detKm1
              endif
              ! The last two rows of the pressure equation matrix are
              !    |    ...  0  igu(kc-1)  b(kc-1)-lam  igl(kc-1)  |
              !    \    ...  0     0        igu(kc)     b(kc)-lam  /
            else
              ! This initialization of det,ddet imply Dirichlet boundary conditions so that first 3 rows
              ! of the matrix are
              !    /  b(2)-lam  igl(2)      0       0     0  ...  |
              !    |  igu(3)  b(3)-lam   igl(3)     0     0  ...  |
              !    |    0       igu43)  b(4)-lam  igl(4)  0  ...  |
              ! which is consistent if the eigenvalue problem is for vertical velocity modes.
              detKm1 = 1.0 ; ddetKm1 = 0.0
              det = (Igu(2) + Igl(2) - lam) ; ddet = -1.0
              ! The last three rows of the w equation matrix are
              !    |    ...   0  igu(kc-1)  b(kc-1)-lam  igl(kc-1)     0       |
              !    |    ...   0     0        igu(kc-1)  b(kc-1)-lam  igl(kc-1) |
              !    \    ...   0     0           0        igu(kc)    b(kc)-lam  /
            endif
            do k=3,kc
              ! Shift variables and rescale rows to avoid over- or underflow.
              detKm2 = c2_scale*detKm1 ; ddetKm2 = c2_scale*ddetKm1
              detKm1 = c2_scale*det    ; ddetKm1 = c2_scale*ddet

              det = (Igu(k)+Igl(k)-lam)*detKm1 - (Igu(k)*Igl(k-1))*detKm2
              ddet = (Igu(k)+Igl(k)-lam)*ddetKm1 - (Igu(k)*Igl(k-1))*ddetKm2 - detKm1

              ! Rescale det & ddet if det is getting too large or too small.
              if (abs(det) > rescale) then
                det = I_rescale*det ; detKm1 = I_rescale*detKm1
                ddet = I_rescale*ddet ; ddetKm1 = I_rescale*ddetKm1
              elseif (abs(det) < I_rescale) then
                det = rescale*det ; detKm1 = rescale*detKm1
                ddet = rescale*ddet ; ddetKm1 = rescale*ddetKm1
              endif
            enddo
            ! Use Newton's method iteration to find a new estimate of lam.
            det_it(itt) = det ; ddet_it(itt) = ddet

            if ((ddet >= 0.0) .or. (-det > -0.5*lam*ddet)) then
              ! lam was not an under-estimate, as intended, so Newton's method
              ! may not be reliable; lam must be reduced, but not by more
              ! than half.
              lam = 0.5 * lam
              dlam = -lam
            else  ! Newton's method is OK.
              dlam = - det / ddet
              lam = lam + dlam
            endif

            if (calc_modal_structure) then
              do k = 1,kc
                Igd(k) = Igu(k) + Igl(k)
              enddo
              call tdma6(kc, -Igu, Igd, -Igl, lam, mode_struct)
              ms_min = mode_struct(1)
              ms_max = mode_struct(1)
              ms_sq = mode_struct(1)**2
              do k = 2,kc
                ms_min = min(ms_min, mode_struct(k))
                ms_max = max(ms_max, mode_struct(k))
                ms_sq = ms_sq + mode_struct(k)**2
              enddo
              if (ms_min<0. .and. ms_max>0.) then ! Any zero crossings => lam is too high
                lam = 0.5 * ( lam - dlam )
                dlam = -lam
                mode_struct(1:kc) = abs(mode_struct(1:kc)) / sqrt( ms_sq )
              else
                mode_struct(1:kc) = mode_struct(1:kc) / sqrt( ms_sq )
              endif
            endif

            if (abs(dlam) < tol2*lam) exit
          enddo

          cg1(i,j) = 0.0
          if (lam > 0.0) cg1(i,j) = 1.0 / sqrt(lam)

          if (present(modal_structure)) then
            if (mode_struct(1)/=0.) then ! Normalize
              mode_struct(1:kc) = mode_struct(1:kc) / mode_struct(1)
            else
              mode_struct(1:kc)=0.
            endif
            ! Note that remapping_core_h requires that the same units be used
            ! for both the source and target grid thicknesses, here [H ~> m or kg m-2].
            do k = 1,kc
              Hc_H(k) = GV%Z_to_H * Hc(k)
            enddo
            if (CS%remap_answers_2018) then
              call remapping_core_h(CS%remapping_CS, kc, Hc_H(:), mode_struct, &
                                    nz, h(i,j,:), modal_structure(i,j,:), &
                                    1.0e-30*GV%m_to_H, 1.0e-10*GV%m_to_H)
            else
              call remapping_core_h(CS%remapping_CS, kc, Hc_H(:), mode_struct, &
                                    nz, h(i,j,:), modal_structure(i,j,:), &
                                    GV%H_subroundoff, GV%H_subroundoff)
            endif
          endif
        else
          cg1(i,j) = 0.0
          if (present(modal_structure)) modal_structure(i,j,:) = 0.
        endif
      endif ! cg1 /= 0.0
    else
      cg1(i,j) = 0.0 ! This is a land point.
      if (present(modal_structure)) modal_structure(i,j,:) = 0.
    endif ; enddo ! i-loop
  enddo ! j-loop

end subroutine wave_speed

!> Solve a non-symmetric tridiagonal problem with a scalar contribution to the leading diagonal.
!! This uses the Thomas algorithm rather than the Hallberg algorithm since the matrix is not symmetric.
subroutine tdma6(n, a, b, c, lam, y)
  integer,            intent(in)    :: n !< Number of rows of matrix
  real, dimension(n), intent(in)    :: a !< Lower diagonal   [T2 L-2 ~> s2 m-2]
  real, dimension(n), intent(in)    :: b !< Leading diagonal [T2 L-2 ~> s2 m-2]
  real, dimension(n), intent(in)    :: c !< Upper diagonal   [T2 L-2 ~> s2 m-2]
  real,               intent(in)    :: lam !< Scalar subtracted from leading diagonal [T2 L-2 ~> s2 m-2]
  real, dimension(n), intent(inout) :: y !< RHS on entry, result on exit
  ! Local variables
  integer :: k, l
  real :: beta(n), lambda  ! Temporary variables in [T2 L-2 ~> s2 m-2]
  real :: I_beta(n)        ! Temporary variables in [L2 T-2 ~> m2 s-2]
  real :: yy(n)            ! A temporary variable with the same units as y on entry.

  lambda = lam
  beta(1) = b(1) - lambda
  if (beta(1)==0.) then ! lam was chosen too perfectly
    ! Change lambda and redo this first row
    lambda = (1. + 1.e-5) * lambda
    beta(1) = b(1) - lambda
  endif
  I_beta(1) = 1. / beta(1)
  yy(1) = y(1)
  do k = 2, n
    beta(k) = ( b(k) - lambda ) - a(k) * c(k-1) * I_beta(k-1)
    ! Perhaps the following 0 needs to become a tolerance to handle underflow?
    if (beta(k)==0.) then ! lam was chosen too perfectly
      ! Change lambda and redo everything up to row k
      lambda = (1. + 1.e-5) * lambda
      I_beta(1) = 1. / ( b(1) - lambda )
      do l = 2, k
        I_beta(l) = 1. / ( ( b(l) - lambda ) - a(l) * c(l-1) * I_beta(l-1) )
        yy(l) = y(l) - a(l) * yy(l-1) * I_beta(l-1)
      enddo
    else
      I_beta(k) = 1. / beta(k)
    endif
    yy(k) = y(k) - a(k) * yy(k-1) * I_beta(k-1)
  enddo
  ! The units of y change by a factor of [L2 T-2] in the following lines.
  y(n) = yy(n) * I_beta(n)
  do k = n-1, 1, -1
    y(k) = ( yy(k) - c(k) * y(k+1) ) * I_beta(k)
  enddo
end subroutine tdma6

!> Calculates the wave speeds for the first few barolinic modes.
subroutine wave_speeds(h, tv, G, GV, US, nmodes, cn, CS, full_halos)
  type(ocean_grid_type),                    intent(in)  :: G !< Ocean grid structure
  type(verticalGrid_type),                  intent(in)  :: GV !< Vertical grid structure
  type(unit_scale_type),                    intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)  :: h !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                    intent(in)  :: tv !< Thermodynamic variables
  integer,                                  intent(in)  :: nmodes !< Number of modes
  real, dimension(G%isd:G%ied,G%jsd:G%jed,nmodes), intent(out) :: cn !< Waves speeds [L T-1 ~> m s-1]
  type(wave_speed_CS), optional,            pointer     :: CS !< Control structure for MOM_wave_speed
  logical,             optional,            intent(in)  :: full_halos !< If true, do the calculation
                                                                      !! over the entire computational domain.
  ! Local variables
  real, dimension(SZK_(G)+1) :: &
    dRho_dT, &    ! Partial derivative of density with temperature [R degC-1 ~> kg m-3 degC-1]
    dRho_dS, &    ! Partial derivative of density with salinity [R ppt-1 ~> kg m-3 ppt-1]
    pres, &       ! Interface pressure [Pa]
    T_int, &      ! Temperature interpolated to interfaces [degC]
    S_int, &      ! Salinity interpolated to interfaces [ppt]
    gprime        ! The reduced gravity across each interface [L2 Z-1 T-2 ~> m s-2].
  real, dimension(SZK_(G)) :: &
    Igl, Igu      ! The inverse of the reduced gravity across an interface times
                  ! the thickness of the layer below (Igl) or above (Igu) it, in [T2 L-2 ~> s2 m-2].
  real, dimension(SZK_(G)-1) :: &
    a_diag, b_diag, c_diag
                  ! diagonals of tridiagonal matrix; one value for each
                  ! interface (excluding surface and bottom) [T2 L-2 ~> s2 m-2]
  real, dimension(SZK_(G),SZI_(G)) :: &
    Hf, &         ! Layer thicknesses after very thin layers are combined [Z ~> m]
    Tf, &         ! Layer temperatures after very thin layers are combined [degC]
    Sf, &         ! Layer salinities after very thin layers are combined [ppt]
    Rf            ! Layer densities after very thin layers are combined [R ~> kg m-3]
  real, dimension(SZK_(G)) :: &
    Hc, &         ! A column of layer thicknesses after convective istabilities are removed [Z ~> m]
    Tc, &         ! A column of layer temperatures after convective istabilities are removed [degC]
    Sc, &         ! A column of layer salinites after convective istabilities are removed [ppt]
    Rc            ! A column of layer densities after convective istabilities are removed [R ~> kg m-3]
  real :: c1_thresh  ! if c1 is below this value, don't bother calculating
                     ! cn values for higher modes [L T-1 ~> m s-1]
  real :: det, ddet       ! determinant & its derivative of eigen system
  real :: lam_1           ! approximate mode-1 eigenvalue [T2 L-2 ~> s2 m-2]
  real :: lam_n           ! approximate mode-n eigenvalue [T2 L-2 ~> s2 m-2]
  real :: dlam            ! increment in lam for Newton's method [T2 L-2 ~> s2 m-2]
  real :: lamMin          ! minimum lam value for root searching range [T2 L-2 ~> s2 m-2]
  real :: lamMax          ! maximum lam value for root searching range [T2 L-2 ~> s2 m-2]
  real :: lamInc          ! width of moving window for root searching [T2 L-2 ~> s2 m-2]
  real :: det_l,det_r     ! determinant value at left and right of window
  real :: ddet_l,ddet_r   ! derivative of determinant at left and right of window
  real :: det_sub,ddet_sub! derivative of determinant at subinterval endpoint
  real :: xl,xr           ! lam guesses at left and right of window [T2 L-2 ~> s2 m-2]
  real :: xl_sub          ! lam guess at left of subinterval window [T2 L-2 ~> s2 m-2]
  real,dimension(nmodes) :: &
          xbl,xbr         ! lam guesses bracketing a zero-crossing (root) [T2 L-2 ~> s2 m-2]
  integer :: numint       ! number of widows (intervals) in root searching range
  integer :: nrootsfound  ! number of extra roots found (not including 1st root)
  real :: min_h_frac
  real :: Z_to_Pa  ! A conversion factor from thicknesses (in Z) to pressure (in Pa)
  real, dimension(SZI_(G)) :: &
    htot, hmin, &  ! Thicknesses [Z ~> m]
    H_here, &      ! A thickness [Z ~> m]
    HxT_here, &    ! A layer integrated temperature [degC Z ~> degC m]
    HxS_here, &    ! A layer integrated salinity [ppt Z ~> ppt m]
    HxR_here       ! A layer integrated density [R Z ~> kg m-2]
  real :: speed2_tot ! overestimate of the mode-1 speed squared [L2 T-2 ~> m2 s-2]
  real :: speed2_min ! minimum mode speed (squared) to consider in root searching [L2 T-2 ~> m2 s-2]
  real, parameter :: reduct_factor = 0.5
                     ! factor used in setting speed2_min [nondim]
  real :: I_Hnew   ! The inverse of a new layer thickness [Z-1 ~> m-1]
  real :: drxh_sum ! The sum of density diffrences across interfaces times thicknesses [R Z ~> kg m-2]
  real, parameter :: tol1  = 0.0001, tol2 = 0.001
  real, pointer, dimension(:,:,:) :: T => NULL(), S => NULL()
  real :: g_Rho0  ! G_Earth/Rho0 [L2 T-2 Z-1 R-1 ~> m4 s-2 kg-1].
  integer :: kf(SZI_(G))
  integer, parameter :: max_itt = 10
  logical :: use_EOS    ! If true, density is calculated from T & S using the equation of state.
  real, dimension(SZK_(G)+1) :: z_int
  ! real, dimension(SZK_(G)+1) :: N2  ! The local squared buoyancy frequency [T-2 ~> s-2]
  integer :: nsub       ! number of subintervals used for root finding
  integer, parameter :: sub_it_max = 4
                        ! maximum number of times to subdivide interval
                        ! for root finding (# intervals = 2**sub_it_max)
  logical :: sub_rootfound ! if true, subdivision has located root
  integer :: kc, nrows
  integer :: sub, sub_it
  integer :: i, j, k, k2, itt, is, ie, js, je, nz, row, iint, m, ig, jg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (present(CS)) then
    if (.not. associated(CS)) call MOM_error(FATAL, "MOM_wave_speed: "// &
           "Module must be initialized before it is used.")
  endif

  if (present(full_halos)) then ; if (full_halos) then
    is = G%isd ; ie = G%ied ; js = G%jsd ; je = G%jed
  endif ; endif

  S => tv%S ; T => tv%T
  g_Rho0 = GV%g_Earth / GV%Rho0
  use_EOS = associated(tv%eqn_of_state)
  Z_to_Pa = GV%Z_to_H * GV%H_to_Pa
  c1_thresh = 0.01*US%m_s_to_L_T

  min_h_frac = tol1 / real(nz)
  !$OMP parallel do default(private) shared(is,ie,js,je,nz,h,G,GV,US,min_h_frac,use_EOS,T,S, &
  !$OMP                                     Z_to_Pa,tv,cn,g_Rho0,nmodes)
  do j=js,je
    !   First merge very thin layers with the one above (or below if they are
    ! at the top).  This also transposes the row order so that columns can
    ! be worked upon one at a time.
    do i=is,ie ; htot(i) = 0.0 ; enddo
    do k=1,nz ; do i=is,ie ; htot(i) = htot(i) + h(i,j,k)*GV%H_to_Z ; enddo ; enddo

    do i=is,ie
      hmin(i) = htot(i)*min_h_frac ; kf(i) = 1 ; H_here(i) = 0.0
      HxT_here(i) = 0.0 ; HxS_here(i) = 0.0 ; HxR_here(i) = 0.0
    enddo
    if (use_EOS) then
      do k=1,nz ; do i=is,ie
        if ((H_here(i) > hmin(i)) .and. (h(i,j,k)*GV%H_to_Z > hmin(i))) then
          Hf(kf(i),i) = H_here(i)
          Tf(kf(i),i) = HxT_here(i) / H_here(i)
          Sf(kf(i),i) = HxS_here(i) / H_here(i)
          kf(i) = kf(i) + 1

          ! Start a new layer
          H_here(i) = h(i,j,k)*GV%H_to_Z
          HxT_here(i) = (h(i,j,k)*GV%H_to_Z)*T(i,j,k)
          HxS_here(i) = (h(i,j,k)*GV%H_to_Z)*S(i,j,k)
        else
          H_here(i) = H_here(i) + h(i,j,k)*GV%H_to_Z
          HxT_here(i) = HxT_here(i) + (h(i,j,k)*GV%H_to_Z)*T(i,j,k)
          HxS_here(i) = HxS_here(i) + (h(i,j,k)*GV%H_to_Z)*S(i,j,k)
        endif
      enddo ; enddo
      do i=is,ie ; if (H_here(i) > 0.0) then
        Hf(kf(i),i) = H_here(i)
        Tf(kf(i),i) = HxT_here(i) / H_here(i)
        Sf(kf(i),i) = HxS_here(i) / H_here(i)
      endif ; enddo
    else
      do k=1,nz ; do i=is,ie
        if ((H_here(i) > hmin(i)) .and. (h(i,j,k)*GV%H_to_Z > hmin(i))) then
          Hf(kf(i),i) = H_here(i) ; Rf(kf(i),i) = HxR_here(i) / H_here(i)
          kf(i) = kf(i) + 1

          ! Start a new layer
          H_here(i) = h(i,j,k)*GV%H_to_Z
          HxR_here(i) = (h(i,j,k)*GV%H_to_Z)*GV%Rlay(k)
        else
          H_here(i) = H_here(i) + h(i,j,k)*GV%H_to_Z
          HxR_here(i) = HxR_here(i) + (h(i,j,k)*GV%H_to_Z)*GV%Rlay(k)
        endif
      enddo ; enddo
      do i=is,ie ; if (H_here(i) > 0.0) then
        Hf(kf(i),i) = H_here(i) ; Rf(kf(i),i) = HxR_here(i) / H_here(i)
      endif ; enddo
    endif

    ! From this point, we can work on individual columns without causing memory
    ! to have page faults.
    do i=is,ie
      if (G%mask2dT(i,j) > 0.5) then
        if (use_EOS) then
          pres(1) = 0.0
          do k=2,kf(i)
            pres(k) = pres(k-1) + Z_to_Pa*Hf(k-1,i)
            T_int(k) = 0.5*(Tf(k,i)+Tf(k-1,i))
            S_int(k) = 0.5*(Sf(k,i)+Sf(k-1,i))
          enddo
          call calculate_density_derivs(T_int, S_int, pres, drho_dT, drho_dS, 2, &
                                        kf(i)-1, tv%eqn_of_state, scale=US%kg_m3_to_R)

          ! Sum the reduced gravities to find out how small a density difference
          ! is negligibly small.
          drxh_sum = 0.0
          do k=2,kf(i)
            drxh_sum = drxh_sum + 0.5*(Hf(k-1,i)+Hf(k,i)) * &
                max(0.0,drho_dT(k)*(Tf(k,i)-Tf(k-1,i)) + &
                        drho_dS(k)*(Sf(k,i)-Sf(k-1,i)))
          enddo
        else
          drxh_sum = 0.0
          do k=2,kf(i)
            drxh_sum = drxh_sum + 0.5*(Hf(k-1,i)+Hf(k,i)) * &
                              max(0.0,Rf(k,i)-Rf(k-1,i))
          enddo
        endif
    !   Find gprime across each internal interface, taking care of convective
    ! instabilities by merging layers.
        if (drxh_sum <= 0.0) then
          cn(i,j,:) = 0.0
        else
          ! Merge layers to eliminate convective instabilities or exceedingly
          ! small reduced gravities.
          if (use_EOS) then
            kc = 1
            Hc(1) = Hf(1,i) ; Tc(1) = Tf(1,i) ; Sc(1) = Sf(1,i)
            do k=2,kf(i)
              if ((drho_dT(k)*(Tf(k,i)-Tc(kc)) + drho_dS(k)*(Sf(k,i)-Sc(kc))) * &
                  (Hc(kc) + Hf(k,i)) < 2.0 * tol2*drxh_sum) then
                ! Merge this layer with the one above and backtrack.
                I_Hnew = 1.0 / (Hc(kc) + Hf(k,i))
                Tc(kc) = (Hc(kc)*Tc(kc) + Hf(k,i)*Tf(k,i)) * I_Hnew
                Sc(kc) = (Hc(kc)*Sc(kc) + Hf(k,i)*Sf(k,i)) * I_Hnew
                Hc(kc) = (Hc(kc) + Hf(k,i))
                ! Backtrack to remove any convective instabilities above...  Note
                ! that the tolerance is a factor of two larger, to avoid limit how
                ! far back we go.
                do k2=kc,2,-1
                  if ((drho_dT(k2)*(Tc(k2)-Tc(k2-1)) + drho_dS(k2)*(Sc(k2)-Sc(k2-1))) * &
                      (Hc(k2) + Hc(k2-1)) < tol2*drxh_sum) then
                    ! Merge the two bottommost layers.  At this point kc = k2.
                    I_Hnew = 1.0 / (Hc(kc) + Hc(kc-1))
                    Tc(kc-1) = (Hc(kc)*Tc(kc) + Hc(kc-1)*Tc(kc-1)) * I_Hnew
                    Sc(kc-1) = (Hc(kc)*Sc(kc) + Hc(kc-1)*Sc(kc-1)) * I_Hnew
                    Hc(kc-1) = (Hc(kc) + Hc(kc-1))
                    kc = kc - 1
                  else ; exit ; endif
                enddo
              else
                ! Add a new layer to the column.
                kc = kc + 1
                drho_dS(kc) = drho_dS(k) ; drho_dT(kc) = drho_dT(k)
                Tc(kc) = Tf(k,i) ; Sc(kc) = Sf(k,i) ; Hc(kc) = Hf(k,i)
              endif
            enddo
            ! At this point there are kc layers and the gprimes should be positive.
            do k=2,kc ! Revisit this if non-Boussinesq.
              gprime(k) = g_Rho0 * (drho_dT(k)*(Tc(k)-Tc(k-1)) + &
                                    drho_dS(k)*(Sc(k)-Sc(k-1)))
            enddo
          else  ! .not.use_EOS
            ! Do the same with density directly...
            kc = 1
            Hc(1) = Hf(1,i) ; Rc(1) = Rf(1,i)
            do k=2,kf(i)
              if ((Rf(k,i) - Rc(kc)) * (Hc(kc) + Hf(k,i)) < 2.0*tol2*drxh_sum) then
                ! Merge this layer with the one above and backtrack.
                Rc(kc) = (Hc(kc)*Rc(kc) + Hf(k,i)*Rf(k,i)) / (Hc(kc) + Hf(k,i))
                Hc(kc) = (Hc(kc) + Hf(k,i))
                ! Backtrack to remove any convective instabilities above...  Note
                ! that the tolerance is a factor of two larger, to avoid limit how
                ! far back we go.
                do k2=kc,2,-1
                  if ((Rc(k2)-Rc(k2-1)) * (Hc(k2)+Hc(k2-1)) < tol2*drxh_sum) then
                    ! Merge the two bottommost layers.  At this point kc = k2.
                    Rc(kc-1) = (Hc(kc)*Rc(kc) + Hc(kc-1)*Rc(kc-1)) / (Hc(kc) + Hc(kc-1))
                    Hc(kc-1) = (Hc(kc) + Hc(kc-1))
                    kc = kc - 1
                  else ; exit ; endif
                enddo
              else
                ! Add a new layer to the column.
                kc = kc + 1
                Rc(kc) = Rf(k,i) ; Hc(kc) = Hf(k,i)
              endif
            enddo
            ! At this point there are kc layers and the gprimes should be positive.
            do k=2,kc ! Revisit this if non-Boussinesq.
              gprime(k) = g_Rho0 * (Rc(k)-Rc(k-1))
            enddo
          endif  ! use_EOS

          !-----------------NOW FIND WAVE SPEEDS---------------------------------------
          ig = i + G%idg_offset ; jg = j + G%jdg_offset
          !   Sum the contributions from all of the interfaces to give an over-estimate
          ! of the first-mode wave speed.
          if (kc >= 2) then
            ! Set depth at surface
            z_int(1) = 0.0
            ! initialize speed2_tot
            speed2_tot = 0.0
            ! Calculate Igu, Igl, depth, and N2 at each interior interface
            ! [excludes surface (K=1) and bottom (K=kc+1)]
            do K=2,kc
              Igl(K) = 1.0/(gprime(K)*Hc(k)) ; Igu(K) = 1.0/(gprime(K)*Hc(k-1))
              z_int(K) = z_int(K-1) + Hc(k-1)
              ! N2(K) = US%L_to_Z**2*gprime(K)/(0.5*(Hc(k)+Hc(k-1)))
              speed2_tot = speed2_tot + gprime(K)*(Hc(k-1)+Hc(k))
            enddo
            ! Set stratification for surface and bottom (setting equal to nearest interface for now)
            ! N2(1) = N2(2) ; N2(kc+1) = N2(kc)
            ! Calculate depth at bottom
            z_int(kc+1) = z_int(kc)+Hc(kc)
            ! check that thicknesses sum to total depth
            if (abs(z_int(kc+1)-htot(i)) > 1.e-12*htot(i)) then
              call MOM_error(FATAL, "wave_structure: mismatch in total depths")
            endif

            ! Define the diagonals of the tridiagonal matrix
            ! First, populate interior rows
            do K=3,kc-1
              row = K-1 ! indexing for TD matrix rows
              a_diag(row) = -Igu(K)
              b_diag(row) = Igu(K)+Igl(K)
              c_diag(row) = -Igl(K)
            enddo
            ! Populate top row of tridiagonal matrix
            K=2 ; row = K-1
            a_diag(row) = 0.0
            b_diag(row) = Igu(K)+Igl(K)
            c_diag(row) = -Igl(K)
            ! Populate bottom row of tridiagonal matrix
            K=kc ; row = K-1
            a_diag(row) = -Igu(K)
            b_diag(row) = Igu(K)+Igl(K)
            c_diag(row) = 0.0
            ! Total number of rows in the matrix = number of interior interfaces
            nrows = kc-1

            ! Under estimate the first eigenvalue to start with.
            lam_1 = 1.0 / speed2_tot

            ! Find the first eigen value
            do itt=1,max_itt
              ! calculate the determinant of (A-lam_1*I)
              call tridiag_det(a_diag(1:nrows),b_diag(1:nrows),c_diag(1:nrows), &
                                      nrows,lam_1,det,ddet, row_scale=US%m_s_to_L_T**2)
              ! Use Newton's method iteration to find a new estimate of lam_1
              !det = det_it(itt) ; ddet = ddet_it(itt)
              if ((ddet >= 0.0) .or. (-det > -0.5*lam_1*ddet)) then
                ! lam_1 was not an under-estimate, as intended, so Newton's method
                ! may not be reliable; lam_1 must be reduced, but not by more
                ! than half.
                lam_1 = 0.5 * lam_1
              else  ! Newton's method is OK.
                dlam = - det / ddet
                lam_1 = lam_1 + dlam
                if (abs(dlam) < tol2*lam_1) then
                  ! calculate 1st mode speed
                  if (lam_1 > 0.0) cn(i,j,1) = 1.0 / sqrt(lam_1)
                  exit
                endif
              endif
            enddo

            ! Find other eigen values if c1 is of significant magnitude, > cn_thresh
            nrootsfound = 0    ! number of extra roots found (not including 1st root)
            if (nmodes>1 .and. kc>=nmodes+1 .and. cn(i,j,1)>c1_thresh) then
              ! Set the the range to look for the other desired eigen values
              ! set min value just greater than the 1st root (found above)
              lamMin = lam_1*(1.0 + tol2)
              ! set max value based on a low guess at wavespeed for highest mode
              speed2_min = (reduct_factor*cn(i,j,1)/real(nmodes))**2
              lamMax = 1.0 / speed2_min
              ! set width of interval (not sure about this - BDM)
              lamInc = 0.5*lam_1
              ! set number of intervals within search range
              numint = nint((lamMax - lamMin)/lamInc)

              !   Find intervals containing zero-crossings (roots) of the determinant
              ! that are beyond the first root

              ! find det_l of first interval (det at left endpoint)
              call tridiag_det(a_diag(1:nrows),b_diag(1:nrows),c_diag(1:nrows), &
                               nrows,lamMin,det_l,ddet_l, row_scale=US%m_s_to_L_T**2)
              ! move interval window looking for zero-crossings************************
              do iint=1,numint
                xr = lamMin + lamInc * iint
                xl = xr - lamInc
                call tridiag_det(a_diag(1:nrows),b_diag(1:nrows),c_diag(1:nrows), &
                                 nrows,xr,det_r,ddet_r, row_scale=US%m_s_to_L_T**2)
                if (det_l*det_r < 0.0) then  ! if function changes sign
                  if (det_l*ddet_l < 0.0) then ! if function at left is headed to zero
                    nrootsfound = nrootsfound + 1
                    xbl(nrootsfound) = xl
                    xbr(nrootsfound) = xr
                  else
                    !   function changes sign but has a local max/min in interval,
                    ! try subdividing interval as many times as necessary (or sub_it_max).
                    ! loop that increases number of subintervals:
                    !call MOM_error(WARNING, "determinant changes sign"// &
                    !            "but has a local max/min in interval;"//&
                    !            " reduce increment in lam.")
                    ! begin subdivision loop -------------------------------------------
                    sub_rootfound = .false. ! initialize
                    do sub_it=1,sub_it_max
                      nsub = 2**sub_it ! number of subintervals; nsub=2,4,8,...
                      ! loop over each subinterval:
                      do sub=1,nsub-1,2 ! only check odds; sub = 1; 1,3; 1,3,5,7;...
                        xl_sub = xl + lamInc/(nsub)*sub
                        call tridiag_det(a_diag(1:nrows),b_diag(1:nrows),c_diag(1:nrows), &
                                 nrows,xl_sub,det_sub,ddet_sub, row_scale=US%m_s_to_L_T**2)
                        if (det_sub*det_r < 0.0) then  ! if function changes sign
                          if (det_sub*ddet_sub < 0.0) then ! if function at left is headed to zero
                            sub_rootfound = .true.
                            nrootsfound = nrootsfound + 1
                            xbl(nrootsfound) = xl_sub
                            xbr(nrootsfound) = xr
                            exit ! exit sub loop
                          endif ! headed toward zero
                        endif ! sign change
                      enddo ! sub-loop
                      if (sub_rootfound) exit ! root has been found, exit sub_it loop
                      !   Otherwise, function changes sign but has a local max/min in one of the
                      ! sub intervals, try subdividing again unless sub_it_max has been reached.
                      if (sub_it == sub_it_max) then
                        call MOM_error(WARNING, "wave_speed: root not found "// &
                                       " after sub_it_max subdivisions of original"// &
                                       " interval.")
                      endif ! sub_it == sub_it_max
                    enddo ! sub_it-loop-------------------------------------------------
                  endif ! det_l*ddet_l < 0.0
                endif ! det_l*det_r < 0.0
                ! exit iint-loop if all desired roots have been found
                if (nrootsfound >= nmodes-1) then
                  ! exit if all additional roots found
                  exit
                elseif (iint == numint) then
                  ! oops, lamMax not large enough - could add code to increase (BDM)
                  ! set unfound modes to zero for now (BDM)
                  cn(i,j,nrootsfound+2:nmodes) = 0.0
                else
                  ! else shift interval and keep looking until nmodes or numint is reached
                  det_l = det_r
                  ddet_l = ddet_r
                endif
              enddo ! iint-loop

              ! Use Newton's method to find the roots within the identified windows
              do m=1,nrootsfound ! loop over the root-containing widows (excluding 1st mode)
                lam_n = xbl(m) ! first guess is left edge of window
                do itt=1,max_itt
                  ! calculate the determinant of (A-lam_n*I)
                  call tridiag_det(a_diag(1:nrows),b_diag(1:nrows),c_diag(1:nrows), &
                                   nrows,lam_n,det,ddet, row_scale=US%m_s_to_L_T**2)
                  ! Use Newton's method to find a new estimate of lam_n
                  dlam = - det / ddet
                  lam_n = lam_n + dlam
                  if (abs(dlam) < tol2*lam_1) then
                    ! calculate nth mode speed
                    if (lam_n > 0.0) cn(i,j,m+1) = 1.0 / sqrt(lam_n)
                    exit
                  endif ! within tol
                enddo ! itt-loop
              enddo ! n-loop
            else
              cn(i,j,2:nmodes) = 0.0 ! else too small to worry about
            endif ! if nmodes>1 .and. kc>nmodes .and. c1>c1_thresh
          else
            cn(i,j,:) = 0.0
          endif ! if more than 2 layers
        endif ! if drxh_sum < 0
      else
        cn(i,j,:) = 0.0 ! This is a land point.
      endif ! if not land
    enddo ! i-loop
  enddo ! j-loop

end subroutine wave_speeds

!> Calculate the determinant of a tridiagonal matrix with diagonals a,b-lam,c and its derivative
!! with lam, where lam is constant across rows.  Only the ratio of det to its derivative and their
!! signs are typically used, so internal rescaling by consistent factors are used to avoid
!! over- or underflow.
subroutine tridiag_det(a, b, c, nrows, lam, det_out, ddet_out, row_scale)
  real, dimension(:), intent(in) :: a !< Lower diagonal of matrix (first entry = 0)
  real, dimension(:), intent(in) :: b !< Leading diagonal of matrix (excluding lam)
  real, dimension(:), intent(in) :: c !< Upper diagonal of matrix (last entry = 0)
  integer,            intent(in) :: nrows !< Size of matrix
  real,               intent(in) :: lam !< Value subtracted from b
  real,               intent(out):: det_out !< Determinant
  real,               intent(out):: ddet_out !< Derivative of determinant w.r.t. lam
  real,     optional, intent(in) :: row_scale !< A scaling factor of the rows of the
                                      !! matrix to limit the growth of the determinant
  ! Local variables
  real, dimension(nrows) :: det ! value of recursion function
  real, dimension(nrows) :: ddet ! value of recursion function for derivative
  real, parameter:: rescale = 1024.0**4 ! max value of determinant allowed before rescaling
  real :: rscl
  real :: I_rescale ! inverse of rescale
  integer :: n      ! row (layer interface) index

  if (size(b) /= nrows) call MOM_error(WARNING, "Diagonal b must be same length as nrows.")
  if (size(a) /= nrows) call MOM_error(WARNING, "Diagonal a must be same length as nrows.")
  if (size(c) /= nrows) call MOM_error(WARNING, "Diagonal c must be same length as nrows.")

  I_rescale = 1.0/rescale
  rscl = 1.0 ; if (present(row_scale)) rscl = row_scale

  det(1) = 1.0      ; ddet(1) = 0.0
  det(2) = b(2)-lam ; ddet(2) = -1.0
  do n=3,nrows
    det(n)  = rscl*(b(n)-lam)*det(n-1)  - rscl*(a(n)*c(n-1))*det(n-2)
    ddet(n) = rscl*(b(n)-lam)*ddet(n-1) - rscl*(a(n)*c(n-1))*ddet(n-2) - det(n-1)
    ! Rescale det & ddet if det is getting too large or too small to avoid overflow or underflow.
    if (abs(det(n)) > rescale) then
      det(n)  = I_rescale*det(n)  ; det(n-1)  = I_rescale*det(n-1)
      ddet(n) = I_rescale*ddet(n) ; ddet(n-1) = I_rescale*ddet(n-1)
    elseif (abs(det(n)) < I_rescale) then
      det(n)  = rescale*det(n)  ; det(n-1)  = rescale*det(n-1)
      ddet(n) = rescale*ddet(n) ; ddet(n-1) = rescale*ddet(n-1)
    endif
  enddo
  det_out = det(nrows)
  ddet_out = ddet(nrows) / rscl

end subroutine tridiag_det

!> Initialize control structure for MOM_wave_speed
subroutine wave_speed_init(CS, use_ebt_mode, mono_N2_column_fraction, mono_N2_depth, remap_answers_2018)
  type(wave_speed_CS), pointer :: CS !< Control structure for MOM_wave_speed
  logical, optional, intent(in) :: use_ebt_mode  !< If true, use the equivalent
                                     !! barotropic mode instead of the first baroclinic mode.
  real,    optional, intent(in) :: mono_N2_column_fraction !< The lower fraction of water column over
                                     !! which N2 is limited as monotonic for the purposes of
                                     !! calculating the vertical modal structure.
  real,    optional, intent(in) :: mono_N2_depth !< The depth below which N2 is limited
                                     !! as monotonic for the purposes of calculating the
                                     !! vertical modal structure [Z ~> m].
  logical, optional, intent(in) :: remap_answers_2018 !< If true, use the order of arithmetic and expressions
                                      !! that recover the remapping answers from 2018.  Otherwise
                                      !! use more robust but mathematically equivalent expressions.


  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_wave_speed"  ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "wave_speed_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  ! Write all relevant parameters to the model log.
  call log_version(mdl, version)

  call wave_speed_set_param(CS, use_ebt_mode=use_ebt_mode, mono_N2_column_fraction=mono_N2_column_fraction)

  call initialize_remapping(CS%remapping_CS, 'PLM', boundary_extrapolation=.false., &
                            answers_2018=CS%remap_answers_2018)

end subroutine wave_speed_init

!> Sets internal parameters for MOM_wave_speed
subroutine wave_speed_set_param(CS, use_ebt_mode, mono_N2_column_fraction, mono_N2_depth, remap_answers_2018)
  type(wave_speed_CS), pointer  :: CS !< Control structure for MOM_wave_speed
  logical, optional, intent(in) :: use_ebt_mode  !< If true, use the equivalent
                                      !! barotropic mode instead of the first baroclinic mode.
  real,    optional, intent(in) :: mono_N2_column_fraction !< The lower fraction of water column over
                                      !! which N2 is limited as monotonic for the purposes of
                                      !! calculating the vertical modal structure.
  real,    optional, intent(in) :: mono_N2_depth !< The depth below which N2 is limited
                                      !! as monotonic for the purposes of calculating the
                                      !! vertical modal structure [Z ~> m].
  logical, optional, intent(in) :: remap_answers_2018 !< If true, use the order of arithmetic and expressions
                                      !! that recover the remapping answers from 2018.  Otherwise
                                      !! use more robust but mathematically equivalent expressions.

  if (.not.associated(CS)) call MOM_error(FATAL, &
     "wave_speed_set_param called with an associated control structure.")

  if (present(use_ebt_mode)) CS%use_ebt_mode = use_ebt_mode
  if (present(mono_N2_column_fraction)) CS%mono_N2_column_fraction = mono_N2_column_fraction
  if (present(mono_N2_depth)) CS%mono_N2_depth = mono_N2_depth
  if (present(remap_answers_2018)) CS%remap_answers_2018 = remap_answers_2018

end subroutine wave_speed_set_param

!> \namespace mom_wave_speed

!!
!! Subroutine wave_speed() solves for the first baroclinic mode wave speed.  (It could
!! solve for all the wave speeds, but the iterative approach taken here means
!! that this is not particularly efficient.)
!!
!! If `e(k)` is the perturbation interface height, this means solving for the
!! smallest eigenvalue (`lam` = 1/c^2) of the system
!!
!! \verbatim
!!   -Igu(k)*e(k-1) + (Igu(k)+Igl(k)-lam)*e(k) - Igl(k)*e(k+1) = 0.0
!! \endverbatim
!!
!! with rigid lid boundary conditions e(1) = e(nz+1) = 0.0 giving
!!
!! \verbatim
!!   (Igu(2)+Igl(2)-lam)*e(2) - Igl(2)*e(3) = 0.0
!!   -Igu(nz)*e(nz-1) + (Igu(nz)+Igl(nz)-lam)*e(nz) = 0.0
!! \endverbatim
!!
!! Here
!! \verbatim
!!   Igl(k) = 1.0/(gprime(k)*h(k)) ; Igu(k) = 1.0/(gprime(k)*h(k-1))
!! \endverbatim
!!
!! Alternately, these same eigenvalues can be found from the second smallest
!! eigenvalue of the Montgomery potential (M(k)) calculation:
!!
!! \verbatim
!!   -Igl(k)*M(k-1) + (Igl(k)+Igu(k+1)-lam)*M(k) - Igu(k+1)*M(k+1) = 0.0
!! \endverbatim
!!
!! with rigid lid and flat bottom boundary conditions
!!
!! \verbatim
!!   (Igu(2)-lam)*M(1) - Igu(2)*M(2) = 0.0
!!   -Igl(nz)*M(nz-1) + (Igl(nz)-lam)*M(nz) = 0.0
!! \endverbatim
!!
!! Note that the barotropic mode has been eliminated from the rigid lid
!! interface height equations, hence the matrix is one row smaller.  Without
!! the rigid lid, the top boundary condition is simpler to implement with
!! the M equations.

end module MOM_wave_speed
