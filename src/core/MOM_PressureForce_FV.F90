!> Finite volume pressure gradient (integrated by quadrature or analytically)
module MOM_PressureForce_FV

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, register_diag_field
use MOM_diag_mediator, only : safe_alloc_ptr, diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_PressureForce_Mont, only : set_pbce_Bouss, set_pbce_nonBouss
use MOM_self_attr_load, only : calc_SAL, SAL_CS
use MOM_tidal_forcing, only : calc_tidal_forcing, tidal_forcing_CS
use MOM_tidal_forcing, only : calc_tidal_forcing_legacy
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_domain
use MOM_density_integrals, only : int_density_dz, int_specific_vol_dp
use MOM_density_integrals, only : int_density_dz_generic_plm, int_density_dz_generic_ppm
use MOM_density_integrals, only : int_spec_vol_dp_generic_plm
use MOM_density_integrals, only : int_density_dz_generic_pcm, int_spec_vol_dp_generic_pcm
use MOM_density_integrals, only : diagnose_mass_weight_Z, diagnose_mass_weight_p
use MOM_ALE, only : TS_PLM_edge_values, TS_PPM_edge_values, ALE_CS

implicit none ; private

#include <MOM_memory.h>

public PressureForce_FV_init
public PressureForce_FV_Bouss, PressureForce_FV_nonBouss

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Finite volume pressure gradient control structure
type, public :: PressureForce_FV_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: calculate_SAL  !< If true, calculate self-attraction and loading.
  logical :: tides          !< If true, apply tidal momentum forcing.
  real    :: Rho0           !< The density used in the Boussinesq
                            !! approximation [R ~> kg m-3].
  real    :: GFS_scale      !< A scaling of the surface pressure gradients to
                            !! allow the use of a reduced gravity model [nondim].
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                            !! timing of diagnostic output.
  integer :: MassWghtInterp !< A flag indicating whether and how to use mass weighting in T/S interpolation
  logical :: use_inaccurate_pgf_rho_anom !< If true, uses the older and less accurate
                            !! method to calculate density anomalies, as used prior to
                            !! March 2018.
  logical :: boundary_extrap !< Indicate whether high-order boundary
                            !! extrapolation should be used within boundary cells

  logical :: reconstruct    !< If true, polynomial profiles of T & S will be
                            !! reconstructed and used in the integrals for the
                            !! finite volume pressure gradient calculation.
                            !! The default depends on whether regridding is being used.

  integer :: Recon_Scheme   !< Order of the polynomial of the reconstruction of T & S
                            !! for the finite volume pressure gradient calculation.
                            !! By the default (1) is for a piecewise linear method

  logical :: use_SSH_in_Z0p !< If true, adjust the height at which the pressure used in the
                            !! equation of state is 0 to account for the displacement of the sea
                            !! surface including adjustments for atmospheric or sea-ice pressure.
  logical :: use_stanley_pgf  !< If true, turn on Stanley parameterization in the PGF
  integer :: tides_answer_date !< Recover old answers with tides in Boussinesq mode
  integer :: id_e_tide = -1 !< Diagnostic identifier
  integer :: id_e_tide_eq = -1 !< Diagnostic identifier
  integer :: id_e_tide_sal = -1 !< Diagnostic identifier
  integer :: id_e_sal = -1 !< Diagnostic identifier
  integer :: id_rho_pgf = -1 !< Diagnostic identifier
  integer :: id_rho_stanley_pgf = -1 !< Diagnostic identifier
  integer :: id_p_stanley = -1 !< Diagnostic identifier
  integer :: id_MassWt_u = -1 !< Diagnostic identifier
  integer :: id_MassWt_v = -1 !< Diagnostic identifier
  type(SAL_CS), pointer :: SAL_CSp => NULL() !< SAL control structure
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL() !< Tides control structure
end type PressureForce_FV_CS

contains

!> \brief Non-Boussinesq analytically-integrated finite volume form of pressure gradient
!!
!! Determines the acceleration due to hydrostatic pressure forces, using
!! the analytic finite volume form of the Pressure gradient, and does not
!! make the Boussinesq approximation.
!!
!! To work, the following fields must be set outside of the usual (is:ie,js:je)
!! range before this subroutine is called:
!!   h(isB:ie+1,jsB:je+1), T(isB:ie+1,jsB:je+1), and S(isB:ie+1,jsB:je+1).
subroutine PressureForce_FV_nonBouss(h, tv, PFu, PFv, G, GV, US, CS, ALE_CSp, p_atm, pbce, eta)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thickness [H ~> kg m-2]
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: PFu !< Zonal acceleration [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: PFv !< Meridional acceleration [L T-2 ~> m s-2]
  type(PressureForce_FV_CS),                  intent(in)  :: CS  !< Finite volume PGF control structure
  type(ALE_CS),                               pointer     :: ALE_CSp !< ALE control structure
  real, dimension(:,:),                       pointer     :: p_atm !< The pressure at the ice-ocean
                                                           !! or atmosphere-ocean interface [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), optional, intent(out) :: pbce !< The baroclinic pressure
                                                           !! anomaly in each layer due to eta anomalies
                                                           !! [L2 T-2 H-1 ~> m4 s-2 kg-1].
  real, dimension(SZI_(G),SZJ_(G)),          optional, intent(out) :: eta !< The total column mass used to
                                                           !! calculate PFu and PFv [H ~> kg m-2].
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: p ! Interface pressure [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target :: &
    T_tmp, &    ! Temporary array of temperatures where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [C ~> degC].
    S_tmp       ! Temporary array of salinities where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [S ~> ppt].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    S_t, S_b, & ! Top and bottom edge values for linear reconstructions
                ! of salinity within each layer [S ~> ppt].
    T_t, T_b    ! Top and bottom edge values for linear reconstructions
                ! of temperature within each layer [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    dza, &      ! The change in geopotential anomaly between the top and bottom
                ! of a layer [L2 T-2 ~> m2 s-2].
    intp_dza    ! The vertical integral in depth of the pressure anomaly less
                ! the pressure anomaly at the top of the layer [R L4 T-4 ~> Pa m2 s-2].
  real, dimension(SZI_(G),SZJ_(G))  :: &
    dp, &       ! The (positive) change in pressure across a layer [R L2 T-2 ~> Pa].
    SSH, &      ! The sea surface height anomaly, in depth units [Z ~> m].
    e_sal, &    ! The bottom geopotential anomaly due to self-attraction and loading [Z ~> m].
    e_tide_eq,  & ! The bottom geopotential anomaly due to tidal forces from astronomical sources [Z ~> m].
    e_tide_sal, & ! The bottom geopotential anomaly due to harmonic self-attraction and loading
                  ! specific to tides [Z ~> m].
    e_sal_tide, & ! The summation of self-attraction and loading and tidal forcing [Z ~> m].
    dM          ! The barotropic adjustment to the Montgomery potential to
                ! account for a reduced gravity model [L2 T-2 ~> m2 s-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    za          ! The geopotential anomaly (i.e. g*e + alpha_0*pressure) at the
                ! interfaces [L2 T-2 ~> m2 s-2].

  real, dimension(SZI_(G)) :: Rho_cv_BL !  The coordinate potential density in the deepest variable
                ! density near-surface layer [R ~> kg m-3].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: &
    intx_za     ! The zonal integral of the geopotential anomaly along the
                ! interfaces, divided by the grid spacing [L2 T-2 ~> m2 s-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    intx_dza    ! The change in intx_za through a layer [L2 T-2 ~> m2 s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: &
    inty_za     ! The meridional integral of the geopotential anomaly along the
                ! interfaces, divided by the grid spacing [L2 T-2 ~> m2 s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    inty_dza    ! The change in inty_za through a layer [L2 T-2 ~> m2 s-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    MassWt_u    ! The fractional mass weighting at a u-point [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    MassWt_v    ! The fractional mass weighting at a v-point [nondim].
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density, [R L2 T-2 ~> Pa] (usually 2e7 Pa = 2000 dbar).

  real :: dp_neglect         ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [R L2 T-2 ~> Pa].
  real :: I_gEarth           ! The inverse of GV%g_Earth [T2 Z L-2 ~> s2 m-1]
  real :: alpha_anom         ! The in-situ specific volume, averaged over a
                             ! layer, less alpha_ref [R-1 ~> m3 kg-1].
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_ALE         ! If true, use an ALE pressure reconstruction.
  logical :: use_EOS         ! If true, density is calculated from T & S using an equation of state.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.

  real :: alpha_ref     ! A reference specific volume [R-1 ~> m3 kg-1] that is used
                        ! to reduce the impact of truncation errors.
  real :: rho_in_situ(SZI_(G)) ! The in situ density [R ~> kg m-3].
  real :: Pa_to_H       ! A factor to convert from Pa to the thickness units (H)
                        ! [H T2 R-1 L-2 ~> m Pa-1 or kg m-2 Pa-1].
  real :: H_to_RL2_T2   ! A factor to convert from thickness units (H) to pressure
                        ! units [R L2 T-2 H-1 ~> Pa m-1 or Pa m2 kg-1].
!  real :: oneatm       ! 1 standard atmosphere of pressure in [R L2 T-2 ~> Pa]
  real, parameter :: C1_6 = 1.0/6.0  ! [nondim]
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  nkmb=GV%nk_rho_varies
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  EOSdom(1) = Isq - (G%isd-1) ;  EOSdom(2) = G%iec+1 - (G%isd-1)

  if (.not.CS%initialized) call MOM_error(FATAL, &
       "MOM_PressureForce_FV_nonBouss: Module must be initialized before it is used.")

  if (CS%use_stanley_pgf) call MOM_error(FATAL, &
       "MOM_PressureForce_FV_nonBouss: The Stanley parameterization is not yet"//&
       "implemented in non-Boussinesq mode.")

  use_p_atm = associated(p_atm)
  use_EOS = associated(tv%eqn_of_state)
  use_ALE = .false.
  if (associated(ALE_CSp)) use_ALE = CS%reconstruct .and. use_EOS

  H_to_RL2_T2 = GV%g_Earth*GV%H_to_RZ
  dp_neglect = GV%g_Earth*GV%H_to_RZ * GV%H_subroundoff
  alpha_ref = 1.0 / CS%Rho0
  I_gEarth = 1.0 / GV%g_Earth

  if ((CS%id_MassWt_u > 0) .or. (CS%id_MassWt_v > 0)) then
    MassWt_u(:,:,:) = 0.0 ; MassWt_v(:,:,:) = 0.0
  endif

  if (use_p_atm) then
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      p(i,j,1) = p_atm(i,j)
    enddo ; enddo
  else
    ! oneatm = 101325.0 * US%Pa_to_RL2_T2 ! 1 atm scaled to [R L2 T-2 ~> Pa]
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      p(i,j,1) = 0.0 ! or oneatm
    enddo ; enddo
  endif
  !$OMP parallel do default(shared)
  do j=Jsq,Jeq+1 ; do k=2,nz+1 ; do i=Isq,Ieq+1
    p(i,j,K) = p(i,j,K-1) + H_to_RL2_T2 * h(i,j,k-1)
  enddo ; enddo ; enddo

  if (use_EOS) then
  !   With a bulk mixed layer, replace the T & S of any layers that are
  ! lighter than the buffer layer with the properties of the buffer
  ! layer.  These layers will be massless anyway, and it avoids any
  ! formal calculations with hydrostatically unstable profiles.
    if (nkmb>0) then
      tv_tmp%T => T_tmp ; tv_tmp%S => S_tmp
      tv_tmp%eqn_of_state => tv%eqn_of_state
      do i=Isq,Ieq+1 ; p_ref(i) = tv%P_Ref ; enddo
      !$OMP parallel do default(shared) private(Rho_cv_BL)
      do j=Jsq,Jeq+1
        do k=1,nkmb ; do i=Isq,Ieq+1
          tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
        enddo ; enddo
        call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, Rho_cv_BL(:), &
                               tv%eqn_of_state, EOSdom)
        do k=nkmb+1,nz ; do i=Isq,Ieq+1
          if (GV%Rlay(k) < Rho_cv_BL(i)) then
            tv_tmp%T(i,j,k) = tv%T(i,j,nkmb) ; tv_tmp%S(i,j,k) = tv%S(i,j,nkmb)
          else
            tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
          endif
        enddo ; enddo
      enddo
    else
      tv_tmp%T => tv%T ; tv_tmp%S => tv%S
      tv_tmp%eqn_of_state => tv%eqn_of_state
    endif
  endif

  ! If regridding is activated, do a linear reconstruction of salinity
  ! and temperature across each layer. The subscripts 't' and 'b' refer
  ! to top and bottom values within each layer (these are the only degrees
  ! of freedom needed to know the linear profile).
  if ( use_ALE ) then
    if ( CS%Recon_Scheme == 1 ) then
      call TS_PLM_edge_values(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, CS%boundary_extrap)
    elseif ( CS%Recon_Scheme == 2) then
      call TS_PPM_edge_values(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, CS%boundary_extrap)
    endif
  endif

  !$OMP parallel do default(shared) private(alpha_anom,dp)
  do k=1,nz
    ! Calculate 4 integrals through the layer that are required in the
    ! subsequent calculation.
    if (use_EOS) then
      if ( use_ALE .and. CS%Recon_Scheme > 0 ) then
        if ( CS%Recon_Scheme == 1 ) then
          call int_spec_vol_dp_generic_plm( T_t(:,:,k), T_b(:,:,k), S_t(:,:,k), S_b(:,:,k), &
                    p(:,:,K), p(:,:,K+1), alpha_ref, dp_neglect, p(:,:,nz+1), G%HI, &
                    tv%eqn_of_state, US, dza(:,:,k), intp_dza(:,:,k), intx_dza(:,:,k), inty_dza(:,:,k), &
                    MassWghtInterp=CS%MassWghtInterp)
        elseif ( CS%Recon_Scheme == 2 ) then
          call MOM_error(FATAL, "PressureForce_FV_nonBouss: "//&
                         "int_spec_vol_dp_generic_ppm does not exist yet.")
        !  call int_spec_vol_dp_generic_ppm ( tv%T(:,:,k), T_t(:,:,k), T_b(:,:,k), &
        !            tv%S(:,:,k), S_t(:,:,k), S_b(:,:,k), p(:,:,K), p(:,:,K+1), &
        !            alpha_ref, G%HI, tv%eqn_of_state, dza(:,:,k), intp_dza(:,:,k), &
        !            intx_dza(:,:,k), inty_dza(:,:,k))
        endif
      else
        call int_specific_vol_dp(tv_tmp%T(:,:,k), tv_tmp%S(:,:,k), p(:,:,K), &
                               p(:,:,K+1), alpha_ref, G%HI, tv%eqn_of_state, &
                               US, dza(:,:,k), intp_dza(:,:,k), intx_dza(:,:,k), &
                               inty_dza(:,:,k), bathyP=p(:,:,nz+1), dP_tiny=dp_neglect, &
                               MassWghtInterp=CS%MassWghtInterp)
      endif
      if ((CS%id_MassWt_u > 0) .or. (CS%id_MassWt_v > 0)) &
        call diagnose_mass_weight_p(p(:,:,K), p(:,:,K+1), dp_neglect, p(:,:,nz+1), G%HI, &
                                    MassWt_u(:,:,k), MassWt_v(:,:,k))
    else
      alpha_anom = 1.0 / GV%Rlay(k) - alpha_ref
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dp(i,j) = H_to_RL2_T2 * h(i,j,k)
        dza(i,j,k) = alpha_anom * dp(i,j)
        intp_dza(i,j,k) = 0.5 * alpha_anom * dp(i,j)**2
      enddo ; enddo
      do j=js,je ; do I=Isq,Ieq
        intx_dza(i,j,k) = 0.5 * alpha_anom * (dp(i,j)+dp(i+1,j))
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        inty_dza(i,j,k) = 0.5 * alpha_anom * (dp(i,j)+dp(i,j+1))
      enddo ; enddo
    endif
  enddo

  !   The bottom geopotential anomaly is calculated first so that the increments
  ! to the geopotential anomalies can be reused.  Alternately, the surface
  ! geopotential could be calculated directly with separate calls to
  ! int_specific_vol_dp with alpha_ref=0, and the anomalies used going
  ! downward, which would relieve the need for dza, intp_dza, intx_dza, and
  ! inty_dza to be 3-D arrays.

  ! Sum vertically to determine the surface geopotential anomaly.
  !$OMP parallel do default(shared)
  do j=Jsq,Jeq+1
    do i=Isq,Ieq+1
      za(i,j,nz+1) = alpha_ref*p(i,j,nz+1) - GV%g_Earth*G%bathyT(i,j)
    enddo
    do k=nz,1,-1 ; do i=Isq,Ieq+1
      za(i,j,K) = za(i,j,K+1) + dza(i,j,k)
    enddo ; enddo
  enddo

  ! Calculate and add the self-attraction and loading geopotential anomaly.
  if (CS%calculate_SAL) then
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      SSH(i,j) = (za(i,j,1) - alpha_ref*p(i,j,1)) * I_gEarth - G%Z_ref &
                 - max(-G%bathyT(i,j)-G%Z_ref, 0.0)
    enddo ; enddo
    call calc_SAL(SSH, e_sal, G, CS%SAL_CSp, tmp_scale=US%Z_to_m)

    if ((CS%tides_answer_date>20230630) .or. (.not.GV%semi_Boussinesq) .or. (.not.CS%tides)) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        za(i,j,1) = za(i,j,1) - GV%g_Earth * e_sal(i,j)
      enddo ; enddo
    endif
  endif

  ! Calculate and add the tidal geopotential anomaly.
  if (CS%tides) then
    if ((CS%tides_answer_date>20230630) .or. (.not.GV%semi_Boussinesq)) then
      call calc_tidal_forcing(CS%Time, e_tide_eq, e_tide_sal, G, US, CS%tides_CSp)
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        za(i,j,1) = za(i,j,1) - GV%g_Earth * (e_tide_eq(i,j) + e_tide_sal(i,j))
      enddo ; enddo
    else  ! This block recreates older answers with tides.
      if (.not.CS%calculate_SAL) e_sal(:,:) = 0.0
      call calc_tidal_forcing_legacy(CS%Time, e_sal, e_sal_tide, e_tide_eq, e_tide_sal, &
                                     G, US, CS%tides_CSp)
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        za(i,j,1) = za(i,j,1) - GV%g_Earth * e_sal_tide(i,j)
      enddo ; enddo
    endif
  endif

  ! Find the height anomalies at the interfaces.  If there are no tides and no SAL,
  ! there is no need to correct za, but omitting this changes answers at roundoff.
  do k=1,nz
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      za(i,j,K+1) = za(i,j,K) - dza(i,j,k)
    enddo ; enddo
  enddo

  !   This order of integrating upward and then downward again is necessary with
  ! a nonlinear equation of state, so that the surface geopotentials will go
  ! linearly between the values at thickness points, but the bottom geopotentials
  ! will not now be linear at the sub-grid-scale.  Doing this ensures no motion
  ! with flat isopycnals, even with a nonlinear equation of state.
  !   With an ice-shelf or icebergs, this linearity condition might need to be applied
  ! to a sub-surface interface.
  !$OMP parallel do default(shared)
  do j=js,je ; do I=Isq,Ieq
    intx_za(I,j,1) = 0.5*(za(i,j,1) + za(i+1,j,1))
  enddo ; enddo
  do k=1,nz
    !$OMP parallel do default(shared)
    do j=js,je ; do I=Isq,Ieq
      intx_za(I,j,K+1) = intx_za(I,j,K) - intx_dza(I,j,k)
    enddo ; enddo
  enddo

  !$OMP parallel do default(shared)
  do J=Jsq,Jeq ; do i=is,ie
    inty_za(i,J,1) = 0.5*(za(i,j,1) + za(i,j+1,1))
  enddo ; enddo
  do k=1,nz
    !$OMP parallel do default(shared)
    do J=Jsq,Jeq ; do i=is,ie
      inty_za(i,J,K+1) = inty_za(i,J,K) - inty_dza(i,J,k)
    enddo ; enddo
  enddo

  !$OMP parallel do default(shared) private(dp)
  do k=1,nz
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      dp(i,j) = H_to_RL2_T2 * h(i,j,k)
    enddo ; enddo

    ! Find the horizontal pressure gradient accelerations.
    ! These expressions for the accelerations have been carefully checked in
    ! a set of idealized cases, and should be bug-free.
    do j=js,je ; do I=Isq,Ieq
      PFu(I,j,k) = ( ((za(i,j,K+1)*dp(i,j) + intp_dza(i,j,k)) - &
                      (za(i+1,j,K+1)*dp(i+1,j) + intp_dza(i+1,j,k))) + &
                     ((dp(i+1,j) - dp(i,j)) * intx_za(I,j,K+1) - &
                      (p(i+1,j,K) - p(i,j,K)) * intx_dza(I,j,k)) ) * &
                   (2.0*G%IdxCu(I,j) / ((dp(i,j) + dp(i+1,j)) + dp_neglect))
    enddo ; enddo

    do J=Jsq,Jeq ; do i=is,ie
      PFv(i,J,k) = (((za(i,j,K+1)*dp(i,j) + intp_dza(i,j,k)) - &
                     (za(i,j+1,K+1)*dp(i,j+1) + intp_dza(i,j+1,k))) + &
                    ((dp(i,j+1) - dp(i,j)) * inty_za(i,J,K+1) - &
                     (p(i,j+1,K) - p(i,j,K)) * inty_dza(i,J,k))) * &
                    (2.0*G%IdyCv(i,J) / ((dp(i,j) + dp(i,j+1)) + dp_neglect))
    enddo ; enddo
  enddo

  if (CS%GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    if (use_EOS) then
      !$OMP parallel do default(shared) private(rho_in_situ)
      do j=Jsq,Jeq+1
        call calculate_density(tv_tmp%T(:,j,1), tv_tmp%S(:,j,1), p(:,j,1), rho_in_situ, &
                               tv%eqn_of_state, EOSdom)

        do i=Isq,Ieq+1
          dM(i,j) = (CS%GFS_scale - 1.0) * (p(i,j,1)*(1.0/rho_in_situ(i) - alpha_ref) + za(i,j,1))
        enddo
      enddo
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dM(i,j) = (CS%GFS_scale - 1.0) * (p(i,j,1)*(1.0/GV%Rlay(1) - alpha_ref) + za(i,j,1))
      enddo ; enddo
    endif

    !$OMP parallel do default(shared)
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        PFu(I,j,k) = PFu(I,j,k) - (dM(i+1,j) - dM(i,j)) * G%IdxCu(I,j)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv(i,J,k) = PFv(i,J,k) - (dM(i,j+1) - dM(i,j)) * G%IdyCv(i,J)
      enddo ; enddo
    enddo
  endif

  if (present(pbce)) then
    call set_pbce_nonBouss(p, tv_tmp, G, GV, US, CS%GFS_scale, pbce)
  endif

  if (present(eta)) then
    Pa_to_H = 1.0 / (GV%g_Earth * GV%H_to_RZ)
    if (use_p_atm) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = (p(i,j,nz+1) - p_atm(i,j))*Pa_to_H ! eta has the same units as h.
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = p(i,j,nz+1)*Pa_to_H ! eta has the same units as h.
      enddo ; enddo
    endif
  endif

  ! To be consistent with old runs, tidal forcing diagnostic also includes total SAL.
  ! New diagnostics are given for each individual field.
  if (CS%id_e_tide>0) call post_data(CS%id_e_tide, e_sal+e_tide_eq+e_tide_sal, CS%diag)
  if (CS%id_e_sal>0) call post_data(CS%id_e_sal, e_sal, CS%diag)
  if (CS%id_e_tide_eq>0) call post_data(CS%id_e_tide_eq, e_tide_eq, CS%diag)
  if (CS%id_e_tide_sal>0) call post_data(CS%id_e_tide_sal, e_tide_sal, CS%diag)
  if (CS%id_MassWt_u>0) call post_data(CS%id_MassWt_u, MassWt_u, CS%diag)
  if (CS%id_MassWt_v>0) call post_data(CS%id_MassWt_v, MassWt_v, CS%diag)

end subroutine PressureForce_FV_nonBouss

!> \brief Boussinesq analytically-integrated finite volume form of pressure gradient
!!
!! Determines the acceleration due to hydrostatic pressure forces, using
!! the finite volume form of the terms and analytic integrals in depth.
!!
!! To work, the following fields must be set outside of the usual (is:ie,js:je)
!! range before this subroutine is called:
!!   h(isB:ie+1,jsB:je+1), T(isB:ie+1,jsB:je+1), and S(isB:ie+1,jsB:je+1).
subroutine PressureForce_FV_Bouss(h, tv, PFu, PFv, G, GV, US, CS, ALE_CSp, p_atm, pbce, eta)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thickness [H ~> m]
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: PFu !< Zonal acceleration [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: PFv !< Meridional acceleration [L T-2 ~> m s-2]
  type(PressureForce_FV_CS),                  intent(in)  :: CS  !< Finite volume PGF control structure
  type(ALE_CS),                               pointer     :: ALE_CSp !< ALE control structure
  real, dimension(:,:),                       pointer     :: p_atm !< The pressure at the ice-ocean
                                                         !! or atmosphere-ocean interface [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), optional, intent(out) :: pbce !< The baroclinic pressure
                                                         !! anomaly in each layer due to eta anomalies
                                                         !! [L2 T-2 H-1 ~> m s-2].
  real, dimension(SZI_(G),SZJ_(G)),          optional, intent(out) :: eta !< The sea-surface height used to
                                                         !! calculate PFu and PFv [H ~> m], with any
                                                         !! tidal contributions.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: e ! Interface height in depth units [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G))  :: &
    e_sal_tide, & ! The summation of self-attraction and loading and tidal forcing [Z ~> m].
    e_sal, &      ! The bottom geopotential anomaly due to self-attraction and loading [Z ~> m].
    e_tide_eq,  & ! The bottom geopotential anomaly due to tidal forces from astronomical sources
                  ! [Z ~> m].
    e_tide_sal, & ! The bottom geopotential anomaly due to harmonic self-attraction and loading
                  ! specific to tides [Z ~> m].
    Z_0p, &       ! The height at which the pressure used in the equation of state is 0 [Z ~> m]
    SSH, &      ! The sea surface height anomaly, in depth units [Z ~> m].
    dM          ! The barotropic adjustment to the Montgomery potential to
                ! account for a reduced gravity model [L2 T-2 ~> m2 s-2].
  real, dimension(SZI_(G)) :: &
    Rho_cv_BL   ! The coordinate potential density in the deepest variable
                ! density near-surface layer [R ~> kg m-3].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    dz_geo      ! The change in geopotential thickness through a layer [L2 T-2 ~> m2 s-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    pa          ! The pressure anomaly (i.e. pressure + g*RHO_0*e) at the
                ! the interface atop a layer [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    dpa, &      ! The change in pressure anomaly between the top and bottom
                ! of a layer [R L2 T-2 ~> Pa].
    intz_dpa    ! The vertical integral in depth of the pressure anomaly less the
                ! pressure anomaly at the top of the layer [H R L2 T-2 ~> m Pa].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: &
    intx_pa     ! The zonal integral of the pressure anomaly along the interface
                ! atop a layer, divided by the grid spacing [R L2 T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    intx_dpa    ! The change in intx_pa through a layer [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: &
    inty_pa     ! The meridional integral of the pressure anomaly along the
                ! interface atop a layer, divided by the grid spacing [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    inty_dpa    ! The change in inty_pa through a layer [R L2 T-2 ~> Pa].

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target :: &
    T_tmp, &    ! Temporary array of temperatures where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [C ~> degC].
    S_tmp       ! Temporary array of salinities where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [S ~> ppt].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    S_t, S_b, & ! Top and bottom edge values for linear reconstructions
                ! of salinity within each layer [S ~> ppt].
    T_t, T_b    ! Top and bottom edge values for linear reconstructions
                ! of temperature within each layer [C ~> degC].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    MassWt_u    ! The fractional mass weighting at a u-point [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    MassWt_v    ! The fractional mass weighting at a v-point [nondim].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    rho_pgf, rho_stanley_pgf ! Density [R ~> kg m-3] from EOS with and without SGS T variance
                             ! in Stanley parameterization.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    p_stanley   ! Pressure [R L2 T-2 ~> Pa] estimated with Rho_0
  real :: zeros(SZI_(G))     ! An array of zero values that can be used as an argument [various]
  real :: rho_in_situ(SZI_(G)) ! The in situ density [R ~> kg m-3].
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density, [R L2 T-2 ~> Pa] (usually 2e7 Pa = 2000 dbar).
  real :: p0(SZI_(G))        ! An array of zeros to use for pressure [R L2 T-2 ~> Pa].
  real :: h_neglect          ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [H ~> m].
  real :: I_Rho0             ! The inverse of the Boussinesq reference density [R-1 ~> m3 kg-1].
  real :: G_Rho0             ! G_Earth / Rho0 in [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1].
  real :: I_g_rho            ! The inverse of the density times the gravitational acceleration [Z T2 L-2 R-1 ~> m Pa-1]
  real :: rho_ref            ! The reference density [R ~> kg m-3].
  real :: dz_neglect         ! A minimal thickness [Z ~> m], like e.
  real :: H_to_RL2_T2        ! A factor to convert from thickness units (H) to pressure
                             ! units [R L2 T-2 H-1 ~> Pa m-1 or Pa m2 kg-1].
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_ALE         ! If true, use an ALE pressure reconstruction.
  logical :: use_EOS         ! If true, density is calculated from T & S using an equation of state.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.
  real, parameter :: C1_6 = 1.0/6.0 ! [nondim]
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer, dimension(2) :: EOSdom_h ! The i-computational domain for the equation of state at tracer points
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer :: i, j, k

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  nkmb=GV%nk_rho_varies
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  EOSdom(1) = Isq - (G%isd-1) ;  EOSdom(2) = G%iec+1 - (G%isd-1)

  if (.not.CS%initialized) call MOM_error(FATAL, &
       "MOM_PressureForce_FV_Bouss: Module must be initialized before it is used.")

  use_p_atm = associated(p_atm)
  use_EOS = associated(tv%eqn_of_state)
  do i=Isq,Ieq+1 ; p0(i) = 0.0 ; enddo
  use_ALE = .false.
  if (associated(ALE_CSp)) use_ALE = CS%reconstruct .and. use_EOS

  h_neglect = GV%H_subroundoff
  dz_neglect = GV%dZ_subroundoff
  I_Rho0 = 1.0 / GV%Rho0
  G_Rho0 = GV%g_Earth / GV%Rho0
  rho_ref = CS%Rho0

  if ((CS%id_MassWt_u > 0) .or. (CS%id_MassWt_v > 0)) then
    MassWt_u(:,:,:) = 0.0 ; MassWt_v(:,:,:) = 0.0
  endif

  if (CS%tides_answer_date>20230630) then
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e(i,j,nz+1) = -G%bathyT(i,j)
    enddo ; enddo

    ! Calculate and add the self-attraction and loading geopotential anomaly.
    if (CS%calculate_SAL) then
      !   Determine the surface height anomaly for calculating self attraction
      ! and loading.  This should really be based on bottom pressure anomalies,
      ! but that is not yet implemented, and the current form is correct for
      ! barotropic tides.
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1
        do i=Isq,Ieq+1
          SSH(i,j) = min(-G%bathyT(i,j) - G%Z_ref, 0.0)
        enddo
        do k=1,nz ; do i=Isq,Ieq+1
          SSH(i,j) = SSH(i,j) + h(i,j,k)*GV%H_to_Z
        enddo ; enddo
      enddo
      call calc_SAL(SSH, e_sal, G, CS%SAL_CSp, tmp_scale=US%Z_to_m)
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        e(i,j,nz+1) = e(i,j,nz+1) - e_sal(i,j)
      enddo ; enddo
    endif

    ! Calculate and add the tidal geopotential anomaly.
    if (CS%tides) then
      call calc_tidal_forcing(CS%Time, e_tide_eq, e_tide_sal, G, US, CS%tides_CSp)
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        e(i,j,nz+1) = e(i,j,nz+1) - (e_tide_eq(i,j) + e_tide_sal(i,j))
      enddo ; enddo
    endif
  else  ! Old answers
    ! Calculate and add the self-attraction and loading geopotential anomaly.
    if (CS%calculate_SAL) then
      !   Determine the surface height anomaly for calculating self attraction
      ! and loading.  This should really be based on bottom pressure anomalies,
      ! but that is not yet implemented, and the current form is correct for
      ! barotropic tides.
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1
        do i=Isq,Ieq+1
          SSH(i,j) = min(-G%bathyT(i,j) - G%Z_ref, 0.0)
        enddo
        do k=1,nz ; do i=Isq,Ieq+1
          SSH(i,j) = SSH(i,j) + h(i,j,k)*GV%H_to_Z
        enddo ; enddo
      enddo
      call calc_SAL(SSH, e_sal, G, CS%SAL_CSp, tmp_scale=US%Z_to_m)
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        e_sal(i,j) = 0.0
      enddo ; enddo
    endif

    ! Calculate and add the tidal geopotential anomaly.
    if (CS%tides) then
      call calc_tidal_forcing_legacy(CS%Time, e_sal, e_sal_tide, e_tide_eq, e_tide_sal, &
                                     G, US, CS%tides_CSp)
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        e(i,j,nz+1) = -(G%bathyT(i,j) + e_sal_tide(i,j))
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        e(i,j,nz+1) = -(G%bathyT(i,j) + e_sal(i,j))
      enddo ; enddo
    endif
  endif

  !$OMP parallel do default(shared)
  do j=Jsq,Jeq+1 ; do k=nz,1,-1 ; do i=Isq,Ieq+1
    e(i,j,K) = e(i,j,K+1) + h(i,j,k)*GV%H_to_Z
  enddo ; enddo ; enddo

  if (use_EOS) then
    if (nkmb>0) then
      ! With a bulk mixed layer, replace the T & S of any layers that are lighter than the buffer
      ! layer with the properties of the buffer layer.  These layers will be massless anyway, and
      ! it avoids any formal calculations with hydrostatically unstable profiles.
      tv_tmp%T => T_tmp ; tv_tmp%S => S_tmp
      tv_tmp%eqn_of_state => tv%eqn_of_state

      do i=Isq,Ieq+1 ; p_ref(i) = tv%P_Ref ; enddo
      !$OMP parallel do default(shared) private(Rho_cv_BL)
      do j=Jsq,Jeq+1
        do k=1,nkmb ; do i=Isq,Ieq+1
          tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
        enddo ; enddo
        call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, Rho_cv_BL(:), &
                               tv%eqn_of_state, EOSdom)

        do k=nkmb+1,nz ; do i=Isq,Ieq+1
          if (GV%Rlay(k) < Rho_cv_BL(i)) then
            tv_tmp%T(i,j,k) = tv%T(i,j,nkmb) ; tv_tmp%S(i,j,k) = tv%S(i,j,nkmb)
          else
            tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
          endif
        enddo ; enddo
      enddo
    else
      tv_tmp%T => tv%T ; tv_tmp%S => tv%S
      tv_tmp%eqn_of_state => tv%eqn_of_state
    endif
  endif

  ! If regridding is activated, do a linear reconstruction of salinity
  ! and temperature across each layer. The subscripts 't' and 'b' refer
  ! to top and bottom values within each layer (these are the only degrees
  ! of freedom needed to know the linear profile).
  if ( use_ALE ) then
    if ( CS%Recon_Scheme == 1 ) then
      call TS_PLM_edge_values(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, CS%boundary_extrap)
    elseif ( CS%Recon_Scheme == 2 ) then
      call TS_PPM_edge_values(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, CS%boundary_extrap)
    endif
  endif

  ! Set the surface boundary conditions on pressure anomaly and its horizontal
  ! integrals, assuming that the surface pressure anomaly varies linearly
  ! in x and y.
  if (use_p_atm) then
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      pa(i,j,1) = (rho_ref*GV%g_Earth)*(e(i,j,1) - G%Z_ref) + p_atm(i,j)
    enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      pa(i,j,1) = (rho_ref*GV%g_Earth)*(e(i,j,1) - G%Z_ref)
    enddo ; enddo
  endif

  if (CS%use_SSH_in_Z0p .and. use_p_atm) then
    I_g_rho = 1.0 / (CS%rho0*GV%g_Earth)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      Z_0p(i,j) = e(i,j,1) + p_atm(i,j) * I_g_rho
    enddo ; enddo
  elseif (CS%use_SSH_in_Z0p) then
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      Z_0p(i,j) = e(i,j,1)
    enddo ; enddo
  else
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      Z_0p(i,j) = G%Z_ref
    enddo ; enddo
  endif

  do k=1,nz
    ! Calculate 4 integrals through the layer that are required in the
    ! subsequent calculation.
    if (use_EOS) then
      ! The following routine computes the integrals that are needed to
      ! calculate the pressure gradient force. Linear profiles for T and S are
      ! assumed when regridding is activated. Otherwise, the previous version
      ! is used, whereby densities within each layer are constant no matter
      ! where the layers are located.
      if ( use_ALE .and. CS%Recon_Scheme > 0 ) then
        if ( CS%Recon_Scheme == 1 ) then
          call int_density_dz_generic_plm(k, tv, T_t, T_b, S_t, S_b, e, &
                    rho_ref, CS%Rho0, GV%g_Earth, dz_neglect, G%bathyT, &
                    G%HI, GV, tv%eqn_of_state, US, CS%use_stanley_pgf, dpa(:,:,k), intz_dpa(:,:,k), &
                    intx_dpa(:,:,k), inty_dpa(:,:,k), &
                    MassWghtInterp=CS%MassWghtInterp, &
                    use_inaccurate_form=CS%use_inaccurate_pgf_rho_anom, Z_0p=Z_0p)
        elseif ( CS%Recon_Scheme == 2 ) then
          call int_density_dz_generic_ppm(k, tv, T_t, T_b, S_t, S_b, e, &
                    rho_ref, CS%Rho0, GV%g_Earth, dz_neglect, G%bathyT, &
                    G%HI, GV, tv%eqn_of_state, US, CS%use_stanley_pgf, dpa(:,:,k), intz_dpa(:,:,k), &
                    intx_dpa(:,:,k), inty_dpa(:,:,k), &
                    MassWghtInterp=CS%MassWghtInterp, Z_0p=Z_0p)
        endif
      else
        call int_density_dz(tv_tmp%T(:,:,k), tv_tmp%S(:,:,k), e(:,:,K), e(:,:,K+1), &
                  rho_ref, CS%Rho0, GV%g_Earth, G%HI, tv%eqn_of_state, US, dpa(:,:,k), &
                  intz_dpa(:,:,k), intx_dpa(:,:,k), inty_dpa(:,:,k), G%bathyT, dz_neglect, &
                  CS%MassWghtInterp, Z_0p=Z_0p)
      endif
      if (GV%Z_to_H /= 1.0) then
        !$OMP parallel do default(shared)
        do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
          intz_dpa(i,j,k) = intz_dpa(i,j,k)*GV%Z_to_H
        enddo ; enddo
      endif
      if ((CS%id_MassWt_u > 0) .or. (CS%id_MassWt_v > 0)) &
        call diagnose_mass_weight_Z(e(:,:,K), e(:,:,K+1), dz_neglect, G%bathyT, G%HI, &
                                    MassWt_u(:,:,k), MassWt_v(:,:,k))
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dz_geo(i,j) = GV%g_Earth * GV%H_to_Z*h(i,j,k)
        dpa(i,j,k) = (GV%Rlay(k) - rho_ref) * dz_geo(i,j)
        intz_dpa(i,j,k) = 0.5*(GV%Rlay(k) - rho_ref) * dz_geo(i,j)*h(i,j,k)
      enddo ; enddo
      !$OMP parallel do default(shared)
      do j=js,je ; do I=Isq,Ieq
        intx_dpa(I,j,k) = 0.5*(GV%Rlay(k) - rho_ref) * (dz_geo(i,j) + dz_geo(i+1,j))
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=Jsq,Jeq ; do i=is,ie
        inty_dpa(i,J,k) = 0.5*(GV%Rlay(k) - rho_ref) * (dz_geo(i,j) + dz_geo(i,j+1))
      enddo ; enddo
    endif
  enddo

  ! Set the pressure anomalies at the interfaces.
  do k=1,nz
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      pa(i,j,K+1) = pa(i,j,K) + dpa(i,j,k)
    enddo ; enddo
  enddo

  ! Set the surface boundary conditions on the horizontally integrated pressure anomaly,
  ! assuming that the surface pressure anomaly varies linearly in x and y.
  ! If there is an ice-shelf or icebergs, this linear variation would need to be applied
  ! to an interior interface.
  !$OMP parallel do default(shared)
  do j=js,je ; do I=Isq,Ieq
    intx_pa(I,j,1) = 0.5*(pa(i,j,1) + pa(i+1,j,1))
  enddo ; enddo
  do k=1,nz
    !$OMP parallel do default(shared)
    do j=js,je ; do I=Isq,Ieq
      intx_pa(I,j,K+1) = intx_pa(I,j,K) + intx_dpa(I,j,k)
    enddo ; enddo
  enddo

  !$OMP parallel do default(shared)
  do J=Jsq,Jeq ; do i=is,ie
    inty_pa(i,J,1) = 0.5*(pa(i,j,1) + pa(i,j+1,1))
  enddo ; enddo
  do k=1,nz
    !$OMP parallel do default(shared)
    do J=Jsq,Jeq ; do i=is,ie
      inty_pa(i,J,K+1) = inty_pa(i,J,K) + inty_dpa(i,J,k)
    enddo ; enddo
  enddo

  ! Compute pressure gradient in x direction
  !$OMP parallel do default(shared)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    PFu(I,j,k) = (((pa(i,j,K)*h(i,j,k) + intz_dpa(i,j,k)) - &
                   (pa(i+1,j,K)*h(i+1,j,k) + intz_dpa(i+1,j,k))) + &
                  ((h(i+1,j,k) - h(i,j,k)) * intx_pa(I,j,K) - &
                   (e(i+1,j,K+1) - e(i,j,K+1)) * intx_dpa(I,j,k) * GV%Z_to_H)) * &
                 ((2.0*I_Rho0*G%IdxCu(I,j)) / &
                  ((h(i,j,k) + h(i+1,j,k)) + h_neglect))
  enddo ; enddo ; enddo

  ! Compute pressure gradient in y direction
  !$OMP parallel do default(shared)
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    PFv(i,J,k) = (((pa(i,j,K)*h(i,j,k) + intz_dpa(i,j,k)) - &
                   (pa(i,j+1,K)*h(i,j+1,k) + intz_dpa(i,j+1,k))) + &
                  ((h(i,j+1,k) - h(i,j,k)) * inty_pa(i,J,K) - &
                   (e(i,j+1,K+1) - e(i,j,K+1)) * inty_dpa(i,J,k) * GV%Z_to_H)) * &
                 ((2.0*I_Rho0*G%IdyCv(i,J)) / &
                  ((h(i,j,k) + h(i,j+1,k)) + h_neglect))
  enddo ; enddo ; enddo

  if (CS%GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    if (use_EOS) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1
        if (use_p_atm) then
          call calculate_density(tv_tmp%T(:,j,1), tv_tmp%S(:,j,1), p_atm(:,j), rho_in_situ, &
                                 tv%eqn_of_state, EOSdom)
        else
          call calculate_density(tv_tmp%T(:,j,1), tv_tmp%S(:,j,1), p0, rho_in_situ, &
                                 tv%eqn_of_state, EOSdom)
        endif
        do i=Isq,Ieq+1
          dM(i,j) = (CS%GFS_scale - 1.0) * (G_Rho0 * rho_in_situ(i)) * (e(i,j,1) - G%Z_ref)
        enddo
      enddo
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dM(i,j) = (CS%GFS_scale - 1.0) * (G_Rho0 * GV%Rlay(1)) * (e(i,j,1) - G%Z_ref)
      enddo ; enddo
    endif

    !$OMP parallel do default(shared)
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        PFu(I,j,k) = PFu(I,j,k) - (dM(i+1,j) - dM(i,j)) * G%IdxCu(I,j)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv(i,J,k) = PFv(i,J,k) - (dM(i,j+1) - dM(i,j)) * G%IdyCv(i,J)
      enddo ; enddo
    enddo
  endif

  if (present(pbce)) then
    call set_pbce_Bouss(e, tv_tmp, G, GV, US, CS%Rho0, CS%GFS_scale, pbce)
  endif

  if (present(eta)) then
    ! eta is the sea surface height relative to a time-invariant geoid, for comparison with
    ! what is used for eta in btstep.  See how e was calculated about 200 lines above.
    if (CS%tides_answer_date>20230630) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = e(i,j,1)*GV%Z_to_H
      enddo ; enddo
      if (CS%tides) then
        !$OMP parallel do default(shared)
        do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
          eta(i,j) = eta(i,j) + (e_tide_eq(i,j)+e_tide_sal(i,j))*GV%Z_to_H
        enddo ; enddo
      endif
      if (CS%calculate_SAL) then
        !$OMP parallel do default(shared)
        do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
          eta(i,j) = eta(i,j) + e_sal(i,j)*GV%Z_to_H
        enddo ; enddo
      endif
    else ! Old answers
      if (CS%tides) then
        !$OMP parallel do default(shared)
        do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
          eta(i,j) = e(i,j,1)*GV%Z_to_H + (e_sal_tide(i,j))*GV%Z_to_H
        enddo ; enddo
      else
        !$OMP parallel do default(shared)
        do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
          eta(i,j) = (e(i,j,1) + e_sal(i,j))*GV%Z_to_H
        enddo ; enddo
      endif
    endif
  endif

  if (CS%use_stanley_pgf) then
    ! Calculated diagnostics related to the Stanley parameterization
    zeros(:) = 0.0
    EOSdom_h(:) = EOS_domain(G%HI)
    if ((CS%id_p_stanley>0) .or. (CS%id_rho_pgf>0) .or. (CS%id_rho_stanley_pgf>0)) then
      ! Find the pressure at the mid-point of each layer.
      H_to_RL2_T2 = GV%g_Earth*GV%H_to_RZ
      if (use_p_atm) then
        do j=js,je ; do i=is,ie
          p_stanley(i,j,1) = 0.5*h(i,j,1) * H_to_RL2_T2 + p_atm(i,j)
        enddo ; enddo
      else
        do j=js,je ; do i=is,ie
          p_stanley(i,j,1) = 0.5*h(i,j,1) * H_to_RL2_T2
        enddo ; enddo
      endif
      do k=2,nz ; do j=js,je ; do i=is,ie
        p_stanley(i,j,k) = p_stanley(i,j,k-1) + 0.5*(h(i,j,k-1) + h(i,j,k)) * H_to_RL2_T2
      enddo ; enddo ; enddo
    endif
    if (CS%id_p_stanley>0) call post_data(CS%id_p_stanley, p_stanley, CS%diag)
    if (CS%id_rho_pgf>0) then
      do k=1,nz ; do j=js,je
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_stanley(:,j,k), zeros, &
                               zeros, zeros, rho_pgf(:,j,k), tv%eqn_of_state, EOSdom_h)
      enddo ; enddo
      call post_data(CS%id_rho_pgf, rho_pgf, CS%diag)
    endif
    if (CS%id_rho_stanley_pgf>0) then
      do k=1,nz ; do j=js,je
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_stanley(:,j,k), tv%varT(:,j,k), &
                               zeros, zeros, rho_stanley_pgf(:,j,k), tv%eqn_of_state, EOSdom_h)
      enddo ; enddo
      call post_data(CS%id_rho_stanley_pgf, rho_stanley_pgf, CS%diag)
    endif
  endif

  ! To be consistent with old runs, tidal forcing diagnostic also includes total SAL.
  ! New diagnostics are given for each individual field.
  if (CS%id_e_tide>0) call post_data(CS%id_e_tide, e_sal_tide, CS%diag)
  if (CS%id_e_sal>0) call post_data(CS%id_e_sal, e_sal, CS%diag)
  if (CS%id_e_tide_eq>0) call post_data(CS%id_e_tide_eq, e_tide_eq, CS%diag)
  if (CS%id_e_tide_sal>0) call post_data(CS%id_e_tide_sal, e_tide_sal, CS%diag)
  if (CS%id_MassWt_u>0) call post_data(CS%id_MassWt_u, MassWt_u, CS%diag)
  if (CS%id_MassWt_v>0) call post_data(CS%id_MassWt_v, MassWt_v, CS%diag)

  if (CS%id_rho_pgf>0) call post_data(CS%id_rho_pgf, rho_pgf, CS%diag)
  if (CS%id_rho_stanley_pgf>0) call post_data(CS%id_rho_stanley_pgf, rho_stanley_pgf, CS%diag)
  if (CS%id_p_stanley>0) call post_data(CS%id_p_stanley, p_stanley, CS%diag)

end subroutine PressureForce_FV_Bouss

!> Initializes the finite volume pressure gradient control structure
subroutine PressureForce_FV_init(Time, G, GV, US, param_file, diag, CS, SAL_CSp, tides_CSp)
  type(time_type), target,    intent(in)    :: Time !< Current model time
  type(ocean_grid_type),      intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),      intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),      intent(in)    :: param_file !< Parameter file handles
  type(diag_ctrl), target,    intent(inout) :: diag !< Diagnostics control structure
  type(PressureForce_FV_CS),  intent(inout) :: CS !< Finite volume PGF control structure
  type(SAL_CS),           intent(in), target, optional :: SAL_CSp !< SAL control structure
  type(tidal_forcing_CS), intent(in), target, optional :: tides_CSp !< Tides control structure

  ! Local variables
  real :: Stanley_coeff    ! Coefficient relating the temperature gradient and sub-gridscale
                           ! temperature variance [nondim]
  integer :: default_answer_date ! Global answer date
  logical :: useMassWghtInterp ! If true, use near-bottom mass weighting for T and S
  logical :: MassWghtInterp_NonBous_bug ! If true, use a buggy mass weighting when non-Boussinesq
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl  ! This module's name.
  logical :: use_ALE       ! If true, use the Vertical Lagrangian Remap algorithm

  CS%initialized = .true.
  CS%diag => diag ; CS%Time => Time
  if (present(tides_CSp)) &
    CS%tides_CSp => tides_CSp
  if (present(SAL_CSp)) &
    CS%SAL_CSp => SAL_CSp

  mdl = "MOM_PressureForce_FV"
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "RHO_PGF_REF", CS%Rho0, &
                 "The reference density that is subtracted off when calculating pressure "//&
                 "gradient forces.  Its inverse is subtracted off of specific volumes when "//&
                 "in non-Boussinesq mode.  The default is RHO_0.", &
                 units="kg m-3", default=GV%Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "TIDES", CS%tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  if (CS%tides) then
    call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
      "This sets the default value for the various _ANSWER_DATE parameters.", &
      default=99991231)
    call get_param(param_file, mdl, "TIDES_ANSWER_DATE", CS%tides_answer_date, &
      "The vintage of self-attraction and loading (SAL) and tidal forcing calculations in "//&
      "Boussinesq mode. Values below 20230701 recover the old answers in which the SAL is "//&
      "part of the tidal forcing calculation.  The change is due to a reordered summation "//&
      "and the difference is only at bit level.", default=20230630)
  endif
  call get_param(param_file, mdl, "CALCULATE_SAL", CS%calculate_SAL, &
                 "If true, calculate self-attraction and loading.", default=CS%tides)
  call get_param(param_file, mdl, "SSH_IN_EOS_PRESSURE_FOR_PGF", CS%use_SSH_in_Z0p, &
                 "If true, include contributions from the sea surface height in the height-based "//&
                 "pressure used in the equation of state calculations for the Boussinesq pressure "//&
                 "gradient forces, including adjustments for atmospheric or sea-ice pressure.", &
                 default=.false., do_not_log=.not.GV%Boussinesq)

  call get_param(param_file, "MOM", "USE_REGRIDDING", use_ALE, &
                 "If True, use the ALE algorithm (regridding/remapping). "//&
                 "If False, use the layered isopycnal algorithm.", default=.false. )
  call get_param(param_file, mdl, "MASS_WEIGHT_IN_PRESSURE_GRADIENT", useMassWghtInterp, &
                 "If true, use mass weighting when interpolating T/S for "//&
                 "integrals near the bathymetry in FV pressure gradient "//&
                 "calculations.", default=.false.)
  call get_param(param_file, mdl, "MASS_WEIGHT_IN_PGF_NONBOUS_BUG", MassWghtInterp_NonBous_bug, &
                 "If true, use a masking bug in non-Boussinesq calculations with mass weighting "//&
                 "when interpolating T/S for integrals near the bathymetry in FV pressure "//&
                 "gradient calculations.", &
                 default=.false., do_not_log=(GV%Boussinesq .or. (.not.useMassWghtInterp)))
  CS%MassWghtInterp = 0
  if (useMassWghtInterp) &
    CS%MassWghtInterp = ibset(CS%MassWghtInterp, 0) ! Same as CS%MassWghtInterp + 1
  if ((.not.GV%Boussinesq) .and. MassWghtInterp_NonBous_bug) &
    CS%MassWghtInterp = ibset(CS%MassWghtInterp, 3) ! Same as CS%MassWghtInterp + 8
  call get_param(param_file, mdl, "USE_INACCURATE_PGF_RHO_ANOM", CS%use_inaccurate_pgf_rho_anom, &
                 "If true, use a form of the PGF that uses the reference density "//&
                 "in an inaccurate way. This is not recommended.", default=.false.)
  call get_param(param_file, mdl, "RECONSTRUCT_FOR_PRESSURE", CS%reconstruct, &
                 "If True, use vertical reconstruction of T & S within "//&
                 "the integrals of the FV pressure gradient calculation. "//&
                 "If False, use the constant-by-layer algorithm. "//&
                 "The default is set by USE_REGRIDDING.", &
                 default=use_ALE )
  call get_param(param_file, mdl, "PRESSURE_RECONSTRUCTION_SCHEME", CS%Recon_Scheme, &
                 "Order of vertical reconstruction of T/S to use in the "//&
                 "integrals within the FV pressure gradient calculation.\n"//&
                 " 0: PCM or no reconstruction.\n"//&
                 " 1: PLM reconstruction.\n"//&
                 " 2: PPM reconstruction.", default=1)
  call get_param(param_file, mdl, "BOUNDARY_EXTRAPOLATION_PRESSURE", CS%boundary_extrap, &
                 "If true, the reconstruction of T & S for pressure in "//&
                 "boundary cells is extrapolated, rather than using PCM "//&
                 "in these cells. If true, the same order polynomial is "//&
                 "used as is used for the interior cells.", default=.true.)
  call get_param(param_file, mdl, "USE_STANLEY_PGF", CS%use_stanley_pgf, &
                 "If true, turn on Stanley SGS T variance parameterization "// &
                 "in PGF code.", default=.false.)
  if (CS%use_stanley_pgf) then
    call get_param(param_file, mdl, "STANLEY_COEFF", Stanley_coeff, &
                 "Coefficient correlating the temperature gradient and SGS T variance.", &
                 units="nondim", default=-1.0, do_not_log=.true.)
    if (Stanley_coeff < 0.0) call MOM_error(FATAL, &
                 "STANLEY_COEFF must be set >= 0 if USE_STANLEY_PGF is true.")

    CS%id_rho_pgf = register_diag_field('ocean_model', 'rho_pgf', diag%axesTL, &
        Time, 'rho in PGF', 'kg m-3', conversion=US%R_to_kg_m3)
    CS%id_rho_stanley_pgf = register_diag_field('ocean_model', 'rho_stanley_pgf', diag%axesTL, &
        Time, 'rho in PGF with Stanley correction', 'kg m-3', conversion=US%R_to_kg_m3)
    CS%id_p_stanley = register_diag_field('ocean_model', 'p_stanley', diag%axesTL, &
        Time, 'p in PGF with Stanley correction', 'Pa', conversion=US%RL2_T2_to_Pa)
  endif
  if (CS%calculate_SAL) then
    CS%id_e_sal = register_diag_field('ocean_model', 'e_sal', diag%axesT1, &
        Time, 'Self-attraction and loading height anomaly', 'meter', conversion=US%Z_to_m)
  endif
  if (CS%tides) then
    CS%id_e_tide = register_diag_field('ocean_model', 'e_tidal', diag%axesT1, Time, &
        'Tidal Forcing Astronomical and SAL Height Anomaly', 'meter', conversion=US%Z_to_m)
    CS%id_e_tide_eq  = register_diag_field('ocean_model', 'e_tide_eq', diag%axesT1, Time, &
        'Equilibrium tides height anomaly', 'meter', conversion=US%Z_to_m)
    CS%id_e_tide_sal = register_diag_field('ocean_model', 'e_tide_sal', diag%axesT1, Time, &
        'Read-in tidal self-attraction and loading height anomaly', 'meter', conversion=US%Z_to_m)
  endif

  CS%id_MassWt_u = register_diag_field('ocean_model', 'MassWt_u', diag%axesCuL, Time, &
        'The fractional mass weighting at u-point PGF calculations', 'nondim')
  CS%id_MassWt_v = register_diag_field('ocean_model', 'MassWt_v', diag%axesCvL, Time, &
        'The fractional mass weighting at v-point PGF calculations', 'nondim')

  CS%GFS_scale = 1.0
  if (GV%g_prime(1) /= GV%g_Earth) CS%GFS_scale = GV%g_prime(1) / GV%g_Earth

  call log_param(param_file, mdl, "GFS / G_EARTH", CS%GFS_scale, units="nondim")

end subroutine PressureForce_FV_init

!> \namespace mom_pressureforce_fv
!!
!! Provides the Boussinesq and non-Boussinesq forms of horizontal accelerations
!! due to pressure gradients using a vertically integrated finite volume form,
!! as described by Adcroft et al., 2008. Integration in the vertical is made
!! either by quadrature or analytically.
!!
!! This form eliminates the thermobaric instabilities that had been a problem with
!! previous forms of the pressure gradient force calculation, as described by
!! Hallberg, 2005.
!!
!! Adcroft, A., R. Hallberg, and M. Harrison, 2008: A finite volume discretization
!! of the pressure gradient force using analytic integration. Ocean Modelling, 22,
!! 106-113. http://doi.org/10.1016/j.ocemod.2008.02.001
!!
!! Hallberg, 2005: A thermobaric instability of Lagrangian vertical coordinate
!! ocean models. Ocean Modelling, 8, 279-300.
!! http://dx.doi.org/10.1016/j.ocemod.2004.01.001

end module MOM_PressureForce_FV
