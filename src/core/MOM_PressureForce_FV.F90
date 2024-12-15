!> Finite volume pressure gradient (integrated by quadrature or analytically)
module MOM_PressureForce_FV

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging, only : hchksum, uvchksum
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
use MOM_variables, only : thermo_var_ptrs, accel_diag_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_spec_vol, EOS_domain
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
  logical :: calculate_SAL = .false. !< If true, calculate self-attraction and loading.
  logical :: sal_use_bpa = .false. !< If true, use bottom pressure anomaly instead of SSH
                                   !! to calculate SAL.
  logical :: tides = .false.       !< If true, apply tidal momentum forcing.
  real    :: rho_ref        !< The reference density that is subtracted off when calculating pressure
                            !! gradient forces [R ~> kg m-3].
  logical :: rho_ref_bug    !< If true, recover a bug that mixes GV%Rho0 and CS%rho_ref in Boussinesq mode.
  real    :: GFS_scale      !< A scaling of the surface pressure gradients to
                            !! allow the use of a reduced gravity model [nondim].
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                            !! timing of diagnostic output.
  integer :: MassWghtInterp !< A flag indicating whether and how to use mass weighting in T/S interpolation
  logical :: correction_intxpa !< If true, apply a correction to the value of intxpa at a selected
                            !! interface under ice, using matching at the end values along with a
                            !! 5-point quadrature integral of the hydrostatic pressure or height
                            !! changes along that interface.  The selected interface is either at the
                            !! ocean's surface or in the interior, depending on reset_intxpa_integral.
  logical :: reset_intxpa_integral !< If true and the surface displacement between adjacent cells
                            !! exceeds the vertical grid spacing, reset intxpa at the interface below
                            !! a trusted interior cell.  (This often applies in ice shelf cavities.)
  logical :: MassWghtInterpVanOnly !< If true, don't do mass weighting of T/S interpolation unless vanished
  logical :: reset_intxpa_flattest !< If true, use flattest interface rather than top for reset integral
                                   !! in cases where no best nonvanished interface
  real    :: h_nonvanished  !< A minimal layer thickness that indicates that a layer is thick enough
                            !! to usefully reestimate the pressure integral across the interface
                            !! below it [H ~> m or kg m-2]
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

  logical :: debug          !< If true, write verbose checksums for debugging purposes.
  logical :: use_SSH_in_Z0p !< If true, adjust the height at which the pressure used in the
                            !! equation of state is 0 to account for the displacement of the sea
                            !! surface including adjustments for atmospheric or sea-ice pressure.
  logical :: use_stanley_pgf  !< If true, turn on Stanley parameterization in the PGF
  logical :: bq_sal_tides = .false. !< If true, use an alternative method for SAL and tides
                                    !! in Boussinesq mode
  integer :: tides_answer_date = 99991231 !< Recover old answers with tides
  integer :: id_e_tide = -1 !< Diagnostic identifier
  integer :: id_e_tidal_eq = -1 !< Diagnostic identifier
  integer :: id_e_tidal_sal = -1 !< Diagnostic identifier
  integer :: id_e_sal = -1 !< Diagnostic identifier
  integer :: id_rho_pgf = -1 !< Diagnostic identifier
  integer :: id_rho_stanley_pgf = -1 !< Diagnostic identifier
  integer :: id_p_stanley = -1 !< Diagnostic identifier
  integer :: id_MassWt_u = -1 !< Diagnostic identifier
  integer :: id_MassWt_v = -1 !< Diagnostic identifier
  integer :: id_sal_u = -1 !< Diagnostic identifier
  integer :: id_sal_v = -1 !< Diagnostic identifier
  integer :: id_tides_u = -1 !< Diagnostic identifier
  integer :: id_tides_v = -1 !< Diagnostic identifier
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
subroutine PressureForce_FV_nonBouss(h, tv, PFu, PFv, G, GV, US, CS, ALE_CSp, ADp, p_atm, pbce, eta)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thickness [H ~> kg m-2]
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: PFu !< Zonal acceleration [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: PFv !< Meridional acceleration [L T-2 ~> m s-2]
  type(PressureForce_FV_CS),                  intent(in)  :: CS  !< Finite volume PGF control structure
  type(ALE_CS),                               pointer     :: ALE_CSp !< ALE control structure
  type(accel_diag_ptrs),                      pointer     :: ADp !< Acceleration diagnostic pointers
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
    SSH, &      ! Sea surfae height anomaly for self-attraction and loading. Used if
                ! CALCULATE_SAL is True and SAL_USE_BPA is False [Z ~> m].
    pbot, &     ! Total bottom pressure for self-attraction and loading. Used if
                ! CALCULATE_SAL is True and SAL_USE_BPA is True [R L2 T-2 ~> Pa].
    e_sal, &    ! The bottom geopotential anomaly due to self-attraction and loading [Z ~> m].
    e_tidal_eq,  & ! The bottom geopotential anomaly due to tidal forces from astronomical sources [Z ~> m].
    e_tidal_sal, & ! The bottom geopotential anomaly due to harmonic self-attraction and loading
                  ! specific to tides [Z ~> m].
    e_sal_and_tide, & ! The summation of self-attraction and loading and tidal forcing, used for recovering
                  ! old answers only [Z ~> m].
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
  real, dimension(SZI_(G),SZJ_(G)) :: &
    T_top, &    ! Temperature of top layer used with correction_intxpa [C ~> degC]
    S_top, &    ! Salinity of top layer used with correction_intxpa [S ~> ppt]
    SpV_top     ! Specific volume anomaly of top layer used with correction_intxpa [R-1 ~> m3 kg-1]
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    intx_za_cor ! Correction for curvature in intx_za [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    inty_za_cor ! Correction for curvature in inty_za [L2 T-2 ~> m2 s-2]

  ! These variables are used with reset_intxpa_integral.  The values are taken from different
  ! interfaces as a function of position.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    T_int_W, T_int_E, & ! Temperatures on the reference interface to the east and west of a u-point [C ~> degC]
    S_int_W, S_int_E, & ! Salinities on the reference interface to the east and west of a u-point [S ~> ppt]
    p_int_W, p_int_E, & ! Pressures on the reference interface to the east and west of a u-point [R L2 T-2 ~> Pa]
    SpV_x_W, SpV_x_E, & ! Specific volume anomalies on the reference interface to the east and west
                        ! of a u-point [R-1 ~> m3 kg-1]
    intx_za_nonlin, &   ! Deviations in the previous version of intx_pa for the reference interface
                        ! from the value that would be obtained from assuming that pressure varies
                        ! linearly with depth along that interface [R L2 T-2 ~> Pa].
    dp_int_x, &         ! The change in x in pressure along the reference interface [R L2 T-2 ~> Pa]
    intx_za_cor_ri      ! The correction to intx_za based on the reference interface calculations [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    T_int_S, T_int_N, & ! Temperatures on the reference interface to the north and south of a v-point [C ~> degC]
    S_int_S, S_int_N, & ! Salinities on the reference interface to the north and south of a v-point [S ~> ppt]
    p_int_S, p_int_N, & ! Pressures on the reference interface to the north and south of a v-point [R L2 T-2 ~> Pa]
    SpV_y_S, SpV_y_N, & ! Specific volume anomalies on the reference interface to the north and south
                        ! of a v-point [R L2 T-2 ~> Pa]
    inty_za_nonlin, &   ! Deviations in the previous version of intx_pa for the reference interface
                        ! from the value that would be obtained from assuming that pressure varies
                        ! linearly with depth along that interface [L2 T-2 ~> m2 s-2].
    dp_int_y, &         ! The change in y in geopotenial height along the reference interface [R L2 T-2 ~> Pa]
    inty_za_cor_ri      ! The correction to inty_za based on the reference interface calculations [L2 T-2 ~> m2 s-2]
  logical, dimension(SZIB_(G),SZJ_(G)) :: &
    seek_x_cor          ! If true, try to find a u-point interface that would provide a better estimate
                        ! of the curvature terms in the intx_pa.
  logical, dimension(SZI_(G),SZJB_(G)) :: &
    seek_y_cor          ! If true, try to find a v-point interface that would provide a better estimate
                        ! of the curvature terms in the inty_pa.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    delta_p_x           ! If using flattest interface for reset integral, store x interface
                        ! differences [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    delta_p_y           ! If using flattest interface for reset integral, store y interface
                        ! differences [R L2 T-2 ~> Pa]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    MassWt_u    ! The fractional mass weighting at a u-point [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    MassWt_v    ! The fractional mass weighting at a v-point [nondim].
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density, [R L2 T-2 ~> Pa] (usually 2e7 Pa = 2000 dbar).
  real :: dp_sfc             ! The change in surface pressure between adjacent cells [R L2 T-2 ~> Pa]

  real :: dp_neglect         ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [R L2 T-2 ~> Pa].
  real :: p_nonvanished      ! nonvanshed pressure [R L2 T-2 ~> Pa]
  real :: I_gEarth           ! The inverse of GV%g_Earth [T2 Z L-2 ~> s2 m-1]
  real :: alpha_anom         ! The in-situ specific volume, averaged over a
                             ! layer, less alpha_ref [R-1 ~> m3 kg-1].
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_ALE         ! If true, use an ALE pressure reconstruction.
  logical :: use_EOS         ! If true, density is calculated from T & S using an equation of state.
  logical :: do_more_k       ! If true, there are still points where a flatter interface remains to be found.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.

  real :: alpha_ref     ! A reference specific volume [R-1 ~> m3 kg-1] that is used
                        ! to reduce the impact of truncation errors.
  real :: rho_in_situ(SZI_(G)) ! The in situ density [R ~> kg m-3].
  real :: Pa_to_H       ! A factor to convert from Pa to the thickness units (H)
                        ! [H T2 R-1 L-2 ~> m Pa-1 or kg m-2 Pa-1].
  real :: H_to_RL2_T2   ! A factor to convert from thickness units (H) to pressure
                        ! units [R L2 T-2 H-1 ~> Pa m-1 or Pa m2 kg-1].
  real :: T5(5)         ! Temperatures and salinities at five quadrature points [C ~> degC]
  real :: S5(5)         ! Salinities at five quadrature points [S ~> ppt]
  real :: p5(5)         ! Pressures at five quadrature points for use with the equation of state [R L2 T-2 ~> Pa]
  real :: SpV5(5)       ! Specific volume anomalies at five quadrature points [R-1 ~> m3 kg-1]
  real :: wt_R          ! A weighting factor [nondim]

  !  real :: oneatm       ! 1 standard atmosphere of pressure in [R L2 T-2 ~> Pa]
  real, parameter :: C1_6 = 1.0/6.0  ! [nondim]
  real, parameter :: C1_90 = 1.0/90.0  ! A rational constant [nondim]
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer, dimension(2) :: EOSdom_u ! The i-computational domain for the equation of state at u-velocity points
  integer, dimension(2) :: EOSdom_v ! The i-computational domain for the equation of state at v-velocity points
  integer :: i, j, k, m

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  nkmb=GV%nk_rho_varies
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  EOSdom(1) = Isq - (G%isd-1) ;  EOSdom(2) = G%iec+1 - (G%isd-1)
  EOSdom_u(1) = Isq - (G%IsdB-1) ; EOSdom_u(2) = Ieq - (G%IsdB-1)
  EOSdom_v(1) = is - (G%isd-1)   ; EOSdom_v(2) = ie - (G%isd-1)

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
  alpha_ref = 1.0 / CS%rho_ref
  I_gEarth = 1.0 / GV%g_Earth
  p_nonvanished = GV%g_Earth*GV%H_to_RZ*CS%h_nonvanished

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
  if ( use_ALE .and. (CS%Recon_Scheme == 1) ) then
      call TS_PLM_edge_values(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, CS%boundary_extrap)
  elseif ( use_ALE .and. (CS%Recon_Scheme == 2) ) then
      call TS_PPM_edge_values(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, CS%boundary_extrap)
  elseif (CS%reset_intxpa_integral) then
    do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      T_b(i,j,k) = tv%T(i,j,k) ; S_b(i,j,k) = tv%S(i,j,k)
    enddo ; enddo ; enddo
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
                    P_surf=p(:,:,1), MassWghtInterp=CS%MassWghtInterp, &
                    MassWghtInterpVanOnly=CS%MassWghtInterpVanOnly, p_nv=p_nonvanished)
        elseif ( CS%Recon_Scheme == 2 ) then
          call MOM_error(FATAL, "PressureForce_FV_nonBouss: "//&
                         "int_spec_vol_dp_generic_ppm does not exist yet.")
        !  call int_spec_vol_dp_generic_ppm ( tv%T(:,:,k), T_t(:,:,k), T_b(:,:,k), &
        !            tv%S(:,:,k), S_t(:,:,k), S_b(:,:,k), p(:,:,K), p(:,:,K+1), &
        !            alpha_ref, G%HI, tv%eqn_of_state, dza(:,:,k), intp_dza(:,:,k), &
        !            intx_dza(:,:,k), inty_dza(:,:,k), P_surf=p(:,:,1), MassWghtInterp=CS%MassWghtInterp)
        endif
      else
        call int_specific_vol_dp(tv_tmp%T(:,:,k), tv_tmp%S(:,:,k), p(:,:,K), &
                               p(:,:,K+1), alpha_ref, G%HI, tv%eqn_of_state, &
                               US, dza(:,:,k), intp_dza(:,:,k), intx_dza(:,:,k), &
                               inty_dza(:,:,k), bathyP=p(:,:,nz+1), P_surf=p(:,:,1), dP_tiny=dp_neglect, &
                               MassWghtInterp=CS%MassWghtInterp, &
                               MassWghtInterpVanOnly=CS%MassWghtInterpVanOnly, p_nv=p_nonvanished)
      endif
      if ((CS%id_MassWt_u > 0) .or. (CS%id_MassWt_v > 0)) &
        call diagnose_mass_weight_p(p(:,:,K), p(:,:,K+1), p(:,:,nz+1), p(:,:,1), dp_neglect, CS%MassWghtInterp, &
                                    G%HI, MassWt_u(:,:,k), MassWt_v(:,:,k), &
                                    MassWghtInterpVanOnly=CS%MassWghtInterpVanOnly, p_nv=p_nonvanished)
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

  ! Calculate and add self-attraction and loading (SAL) geopotential height anomaly to interface height.
  if (CS%calculate_SAL) then
    if (CS%sal_use_bpa) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        pbot(i,j) = p(i,j,nz+1)
      enddo ; enddo
      call calc_SAL(pbot, e_sal, G, CS%SAL_CSp, tmp_scale=US%Z_to_m)
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        SSH(i,j) = (za(i,j,1) - alpha_ref*p(i,j,1)) * I_gEarth - G%Z_ref &
                  - max(-G%bathyT(i,j)-G%Z_ref, 0.0)
      enddo ; enddo
      call calc_SAL(SSH, e_sal, G, CS%SAL_CSp, tmp_scale=US%Z_to_m)
    endif

    ! This gives new answers after the change of separating SAL from tidal forcing module.
    if ((CS%tides_answer_date>20230630) .or. (.not.GV%semi_Boussinesq) .or. (.not.CS%tides)) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        za(i,j,1) = za(i,j,1) - GV%g_Earth * e_sal(i,j)
      enddo ; enddo
    endif
  endif

  ! Calculate and add tidal geopotential height anomaly to interface height.
  if (CS%tides) then
    if ((CS%tides_answer_date>20230630) .or. (.not.GV%semi_Boussinesq)) then
      call calc_tidal_forcing(CS%Time, e_tidal_eq, e_tidal_sal, G, US, CS%tides_CSp)
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        za(i,j,1) = za(i,j,1) - GV%g_Earth * (e_tidal_eq(i,j) + e_tidal_sal(i,j))
      enddo ; enddo
    else  ! This block recreates older answers with tides.
      if (.not.CS%calculate_SAL) e_sal(:,:) = 0.0
      call calc_tidal_forcing_legacy(CS%Time, e_sal, e_sal_and_tide, e_tidal_eq, e_tidal_sal, &
                                     G, US, CS%tides_CSp)
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        za(i,j,1) = za(i,j,1) - GV%g_Earth * e_sal_and_tide(i,j)
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

  if (CS%debug) then
    call hchksum(za, "Pre-correction za", G%HI, haloshift=1, unscale=US%L_T_to_m_s**2)
    call hchksum(p, "Pre-correction p", G%HI, haloshift=1, unscale=US%RL2_T2_to_Pa)
  endif

  !   With an ice-shelf or icebergs, this linearity condition might need to be applied
  ! to a sub-surface interface.
  if (CS%correction_intxpa .or. CS%reset_intxpa_integral) then
    ! Determine surface temperature and salinity for use in the pressure gradient corrections
    if (use_ALE .and. (CS%Recon_Scheme > 0)) then
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        T_top(i,j) = T_t(i,j,1) ; S_top(i,j) = S_t(i,j,1)
      enddo ; enddo
    else
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        T_top(i,j) = tv%T(i,j,1) ; S_top(i,j) = tv%S(i,j,1)
      enddo ; enddo
    endif
  endif

  if (CS%correction_intxpa) then
    ! This version makes a 5 point quadrature correction for hydrostatic variations in surface
    ! pressure under ice.
    !$OMP parallel do default(shared) private(dp_sfc,T5,S5,p5,wt_R,SpV5)
    do j=js,je ; do I=Isq,Ieq
      intx_za_cor(I,j) = 0.0
      dp_sfc = (p(i+1,j,1) - p(i,j,1))
      ! If the changes in pressure and height anomaly were explicable by just a hydrostatic balance,
      ! the implied specific volume would be   SpV_implied = alpha_ref - (dza_x / dp_x)
      if (dp_sfc * (alpha_ref*dp_sfc - (za(i+1,j,1)-za(i,j,1))) > 0.0) then
        T5(1) = T_top(i,j) ; T5(5) = T_top(i+1,j)
        S5(1) = S_top(i,j) ; S5(5) = S_top(i+1,j)
        p5(1) = p(i,j,1)   ; p5(5) = p(i+1,j,1)
        do m=2,4
          wt_R =  0.25*real(m-1)
          T5(m) = T5(1) + (T5(5)-T5(1))*wt_R
          S5(m) = S5(1) + (S5(5)-S5(1))*wt_R
          p5(m) = p5(1) + (p5(5)-p5(1))*wt_R
        enddo !m
        call calculate_spec_vol(T5, S5, p5, SpV5, tv%eqn_of_state, spv_ref=alpha_ref)
        ! See the Boussinesq calculation of inty_pa_cor for the derivation of the following expression.
        intx_za_cor(I,j) = C1_90 * (4.75*(SpV5(5)-SpV5(1)) + 5.5*(SpV5(4)-SpV5(2))) * dp_sfc
        ! Note the consistency with the linear form below because (4.75 + 5.5/2) / 90 = 1/12
      endif
      intx_za(I,j,1) = 0.5*(za(i,j,1) + za(i+1,j,1)) + intx_za_cor(I,j)
    enddo ; enddo
    !$OMP parallel do default(shared) private(dp_sfc,T5,S5,p5,wt_R,SpV5)
    do J=Jsq,Jeq ; do i=is,ie
      inty_za_cor(i,J) = 0.0
      dp_sfc = (p(i,j+1,1) - p(i,j,1))
      if (dp_sfc * (alpha_ref*dp_sfc - (za(i,j+1,1)-za(i,j,1))) > 0.0) then
        ! The pressure/depth relationship has a positive implied specific volume.
        T5(1) = T_top(i,j) ; T5(5) = T_top(i,j+1)
        S5(1) = S_top(i,j) ; S5(5) = S_top(i,j+1)
        p5(1) = p(i,j,1)   ; p5(5) = p(i,j+1,1)
        do m=2,4
          wt_R =  0.25*real(m-1)
          T5(m) = T5(1) + (T5(5)-T5(1))*wt_R
          S5(m) = S5(1) + (S5(5)-S5(1))*wt_R
          p5(m) = p5(1) + (p5(5)-p5(1))*wt_R
        enddo !m
        call calculate_spec_vol(T5, S5, p5, SpV5, tv%eqn_of_state, spv_ref=alpha_ref)
        ! See the Boussinesq calculation of inty_pa_cor for the derivation of the following expression.
        inty_za_cor(i,J) = C1_90 * (4.75*(SpV5(5)-SpV5(1)) + 5.5*(SpV5(4)-SpV5(2))) * dp_sfc
      endif
      inty_za(i,J,1) = 0.5*(za(i,j,1) + za(i,j+1,1)) + inty_za_cor(i,J)
    enddo ; enddo
  else
    !   This order of integrating upward and then downward again is necessary with
    ! a nonlinear equation of state, so that the surface geopotentials will go
    ! linearly between the values at thickness points, but the bottom geopotentials
    ! will not now be linear at the sub-grid-scale.  Doing this ensures no motion
    ! with flat isopycnals, even with a nonlinear equation of state.
    !$OMP parallel do default(shared)
    do j=js,je ; do I=Isq,Ieq
      intx_za(I,j,1) = 0.5*(za(i,j,1) + za(i+1,j,1))
    enddo ; enddo
    !$OMP parallel do default(shared)
    do J=Jsq,Jeq ; do i=is,ie
      inty_za(i,J,1) = 0.5*(za(i,j,1) + za(i,j+1,1))
    enddo ; enddo
  endif

  do k=1,nz
    !$OMP parallel do default(shared)
    do j=js,je ; do I=Isq,Ieq
      intx_za(I,j,K+1) = intx_za(I,j,K) - intx_dza(I,j,k)
    enddo ; enddo
  enddo
  do k=1,nz
    !$OMP parallel do default(shared)
    do J=Jsq,Jeq ; do i=is,ie
      inty_za(i,J,K+1) = inty_za(i,J,K) - inty_dza(i,J,k)
    enddo ; enddo
  enddo

  if (CS%debug) then
    call uvchksum("Prelim int[xy]_za", intx_za, inty_za, G%HI, haloshift=0, &
                  symmetric=G%Domain%symmetric, scalar_pair=.true., unscale=US%L_T_to_m_s**2)
    call uvchksum("Prelim int[xy]_dza", intx_dza, inty_dza, G%HI, haloshift=0, &
                  symmetric=G%Domain%symmetric, scalar_pair=.true., unscale=US%L_T_to_m_s**2)
  endif

  if (CS%reset_intxpa_integral) then
    ! Having stored the pressure gradient info, we can work out where the first nonvanished layers is
    ! reset intx_za there, then adjust intx_za throughout the water column.

    ! Zero out the 2-d arrays that will be set from various reference interfaces.
    T_int_W(:,:) = 0.0 ; S_int_W(:,:) = 0.0 ; p_int_W(:,:) = 0.0
    T_int_E(:,:) = 0.0 ; S_int_E(:,:) = 0.0 ; p_int_E(:,:) = 0.0
    intx_za_nonlin(:,:) = 0.0 ; intx_za_cor_ri(:,:) = 0.0 ; dp_int_x(:,:) = 0.0
    do j=js,je ; do I=Isq,Ieq
      seek_x_cor(I,j) = (G%mask2dCu(I,j) > 0.)
      delta_p_x(I,j)  = 0.0
    enddo ; enddo

    do j=js,je ; do I=Isq,Ieq ; if (seek_x_cor(I,j)) then
      if ((p(i+1,j,2) >= p(i,j,1)) .and. (p(i,j,2) >= p(i+1,j,1))) then
        ! This is the typical case in the open ocean, so use the topmost interface.
        T_int_W(I,j) = T_top(i,j) ; T_int_E(I,j) = T_top(i+1,j)
        S_int_W(I,j) = S_top(i,j) ; S_int_E(I,j) = S_top(i+1,j)
        p_int_W(I,j) = p(i,j,1) ; p_int_E(I,j) = p(i+1,j,1)
        intx_za_nonlin(I,j) = intx_za(I,j,1) - 0.5*(za(i,j,1) + za(i+1,j,1))
        dp_int_x(I,j) = p(i+1,j,1)-p(i,j,1)
        seek_x_cor(I,j) = .false.
      endif
    endif ; enddo ; enddo

    do k=1,nz
      do_more_k = .false.
      do j=js,je ; do I=Isq,Ieq ; if (seek_x_cor(I,j)) then
        ! Find the topmost layer for which both sides are nonvanished and mass-weighting is not
        ! activated in the subgrid interpolation.
        if (((h(i,j,k) > CS%h_nonvanished) .and. (h(i+1,j,k) > CS%h_nonvanished)) .and. &
            (max(0., p(i,j,1)-p(i+1,j,K+1), p(i+1,j,1)-p(i,j,K+1)) <= 0.0)) then
          ! Store properties at the bottom of this cell to get a "good estimate" for intxpa at
          ! the interface below this cell (it might have quadratic pressure dependence if sloped)
          T_int_W(I,j) = T_b(i,j,k) ; T_int_E(I,j) = T_b(i+1,j,k)
          S_int_W(I,j) = S_b(i,j,k) ; S_int_E(I,j) = S_b(i+1,j,k)
          p_int_W(I,j) = p(i,j,K+1) ; p_int_E(I,j) = p(i+1,j,K+1)

          intx_za_nonlin(I,j) = intx_za(I,j,K+1) - 0.5*(za(i,j,K+1) + za(i+1,j,K+1))
          dp_int_x(I,j) = p(i+1,j,K+1)-p(i,j,K+1)
          seek_x_cor(I,j) = .false.
        else
          do_more_k = .true.
        endif
      endif ; enddo ; enddo
      if (.not.do_more_k) exit  ! All reference interfaces have been found, so stop working downward.
    enddo

    if (do_more_k) then
      if (CS%reset_intxpa_flattest) then
      ! There are still points where a correction is needed, so use flattest interface
        do j=js,je ; do I=Isq,Ieq ; if (seek_x_cor(I,j)) then
          ! choose top layer first
          T_int_W(I,j) = T_top(i,j) ; T_int_E(I,j) = T_top(i+1,j)
          S_int_W(I,j) = S_top(i,j) ; S_int_E(I,j) = S_top(i+1,j)
          p_int_W(I,j) = p(i,j,1) ; p_int_E(I,j) = p(i+1,j,1)
          intx_za_nonlin(I,j) = intx_za(I,j,1) - 0.5*(za(i,j,1) + za(i+1,j,1))
          dp_int_x(I,j) = p(i+1,j,1)-p(i,j,1)
          delta_p_x(I,j) = abs(p(i+1,j,1)-p(i,j,1))
          do k=1,nz
            if (abs(p(i+1,j,k+1)-p(i,j,k+1)) < delta_p_x(I,j)) then
              ! bottom of layer is less sloped than top. Use this layer
              delta_p_x(I,j) = abs(p(i+1,j,k+1)-p(i,j,k+1))
              T_int_W(I,j) = T_b(i,j,k) ; T_int_E(I,j) = T_b(i+1,j,k)
              S_int_W(I,j) = S_b(i,j,k) ; S_int_E(I,j) = S_b(i+1,j,k)
              p_int_W(I,j) = p(i,j,K+1) ; p_int_E(I,j) = p(i+1,j,K+1)
              intx_za_nonlin(I,j) = intx_za(I,j,K+1) - 0.5*(za(i,j,K+1) + za(i+1,j,K+1))
              dp_int_x(I,j) = p(i+1,j,K+1)-p(i,j,K+1)
            endif
          enddo
        seek_x_cor(I,j) = .false.
        endif; enddo; enddo;
      else
        ! There are still points where a correction is needed, so use the top interface.
        do j=js,je ; do I=Isq,Ieq ; if (seek_x_cor(I,j)) then
          T_int_W(I,j) = T_top(i,j) ; T_int_E(I,j) = T_top(i+1,j)
          S_int_W(I,j) = S_top(i,j) ; S_int_E(I,j) = S_top(i+1,j)
          p_int_W(I,j) = p(i,j,1) ; p_int_E(I,j) = p(i+1,j,1)
          intx_za_nonlin(I,j) = intx_za(I,j,1) - 0.5*(za(i,j,1) + za(i+1,j,1))
          dp_int_x(I,j) = p(i+1,j,1)-p(i,j,1)
          seek_x_cor(I,j) = .false.
        endif ; enddo ; enddo
      endif
    endif

    do j=js,je
      do I=Isq,Ieq
        ! This expression assumes that temperature and salinity vary linearly with pressure
        ! between the corners of the reference interfaces found above to get a correction to
        ! intx_pa that takes nonlinearities in the equation of state into account.
        ! It is derived from a 5 point quadrature estimate of the integral with a large-scale
        ! linear correction so that the pressures and heights match at the end-points.  It turns
        ! out that this linear correction cancels out the mid-point specific volume.
        ! This can be used without masking because dp_int_x and intx_za_nonlin are 0 over land.
        T5(1) = T_Int_W(I,j) ; S5(1) = S_Int_W(I,j) ; p5(1) = p_Int_W(I,j)
        T5(5) = T_Int_E(I,j) ; S5(5) = S_Int_E(I,j) ; p5(5) = p_Int_E(I,j)
        T5(2) = 0.25*(3.0*T5(1) + T5(5)) ; T5(4) = 0.25*(3.0*T5(5) + T5(1)) ; T5(3) = 0.5*(T5(5) + T5(1))
        S5(2) = 0.25*(3.0*S5(1) + S5(5)) ; S5(4) = 0.25*(3.0*S5(5) + S5(1)) ; S5(3) = 0.5*(S5(5) + S5(1))
        p5(2) = 0.25*(3.0*p5(1) + p5(5)) ; p5(4) = 0.25*(3.0*p5(5) + p5(1)) ; p5(3) = 0.5*(p5(5) + p5(1))
        call calculate_spec_vol(T5, S5, p5, SpV5, tv%eqn_of_state, spv_ref=alpha_ref)

        ! Note the consistency with the linear form below because (4.75 + 5.5/2) / 90 = 1/12
        intx_za_cor_ri(I,j) = C1_90 * (4.75*(SpV5(5)-SpV5(1)) + 5.5*(SpV5(4)-SpV5(2))) * &
                                      dp_int_x(I,j) - intx_za_nonlin(I,j)
      enddo
    enddo

    ! Repeat the calculations above for v-velocity points.
    T_int_S(:,:) = 0.0 ; S_int_S(:,:) = 0.0 ; p_int_S(:,:) = 0.0
    T_int_N(:,:) = 0.0 ; S_int_N(:,:) = 0.0 ; p_int_N(:,:) = 0.0
    inty_za_nonlin(:,:) = 0.0 ; inty_za_cor_ri(:,:) = 0.0 ; dp_int_y(:,:) = 0.0
    do J=Jsq,Jeq ; do i=is,ie
      seek_y_cor(i,J) = (G%mask2dCv(i,J) > 0.)
      delta_p_y(i,J) = 0.0
    enddo ; enddo

    do J=Jsq,Jeq ; do i=is,ie ; if (seek_y_cor(i,J)) then
      if ((p(i,j+1,2) >= p(i,j,1)) .and. (p(i,j,2) >= p(i,j+1,1))) then
        ! This is the typical case in the open ocean, so use the topmost interface.
        T_int_S(i,J) = T_top(i,j) ; T_int_N(i,J) = T_top(i,j+1)
        S_int_S(i,J) = S_top(i,j) ; S_int_N(i,J) = S_top(i,j+1)
        p_int_S(i,J) = p(i,j,1) ; p_int_N(i,J) = p(i,j+1,1)
        inty_za_nonlin(i,J) = inty_za(i,J,1) - 0.5*(za(i,j,1) + za(i,j+1,1))
        dp_int_y(i,J) = p(i,j+1,1) - p(i,j,1)
        seek_y_cor(i,J) = .false.
      endif
    endif ; enddo ; enddo

    do k=1,nz
      do_more_k = .false.
      do J=Jsq,Jeq ; do i=is,ie ; if (seek_y_cor(i,J)) then
        ! Find the topmost layer for which both sides are nonvanished and mass-weighting is not
        ! activated in the subgrid interpolation.
        if (((h(i,j,k) > CS%h_nonvanished) .and. (h(i,j+1,k) > CS%h_nonvanished)) .and. &
            (max(0., p(i,j,1)-p(i,j+1,K+1), p(i,j+1,1)-p(i,j,K+1)) <= 0.0)) then
          ! Store properties at the bottom of this cell to get a "good estimate" for intypa at
          ! the interface below this cell (it might have quadratic pressure dependence if sloped)
          T_int_S(i,J) = T_b(i,j,k) ; T_int_N(i,J) = T_b(i,j+1,k)
          S_int_S(i,J) = S_b(i,j,k) ; S_int_N(i,J) = S_b(i,j+1,k)
          p_int_S(i,J) = p(i,j,K+1) ; p_int_N(i,J) = p(i,j+1,K+1)
          inty_za_nonlin(i,J) = inty_za(i,J,K+1) - 0.5*(za(i,j,K+1) + za(i,j+1,K+1))
          dp_int_y(i,J) = p(i,j+1,K+1) - p(i,j,K+1)
          seek_y_cor(i,J) = .false.
        else
          do_more_k = .true.
        endif
      endif ; enddo ; enddo
      if (.not.do_more_k) exit  ! All reference interfaces have been found, so stop working downward.
    enddo

    if (do_more_k) then
      if (CS%reset_intxpa_flattest) then
        ! There are still points where a correction is needed, so use flattest interface.
        do J=Jsq,Jeq ; do i=is,ie ; if (seek_y_cor(i,J)) then
          ! choose top interface first
          T_int_S(i,J) = T_top(i,j) ; T_int_N(i,J) = T_top(i,j+1)
          S_int_S(i,J) = S_top(i,j) ; S_int_N(i,J) = S_top(i,j+1)
          p_int_S(i,J) = p(i,j,1) ; p_int_N(i,J) = p(i,j+1,1)
          inty_za_nonlin(i,J) = inty_za(i,J,1) - 0.5*(za(i,j,1) + za(i,j+1,1))
          dp_int_y(i,J) = p(i,j+1,1) - p(i,j,1)
          delta_p_y(i,J) = abs(p(i,j+1,1)-p(i,j,1))
          do k=1,nz
            if (abs(p(i,j+1,k+1)-p(i,j,k+1)) < delta_p_y(i,J)) then
              ! bottom of layer is less sloped than top. Use this layer
              delta_p_y(i,J) = abs(p(i,j+1,k+1)-p(i,j,k+1))
              T_int_S(i,J) = T_b(i,j,k) ; T_int_N(i,J) = T_b(i,j+1,k)
              S_int_S(i,J) = S_b(i,j,k) ; S_int_N(i,J) = S_b(i,j+1,k)
              p_int_S(i,J) = p(i,j,K+1) ; p_int_N(i,J) = p(i,j+1,K+1)
              inty_za_nonlin(i,J) = inty_za(i,J,K+1) - 0.5*(za(i,j,K+1) + za(i,j+1,K+1))
              dp_int_y(i,J) = p(i,j+1,K+1) - p(i,j,K+1)
            endif
          enddo
          seek_y_cor(i,J) = .false.
        endif ; enddo ; enddo
      else
        ! There are still points where a correction is needed, so use the top interface.
        do J=Jsq,Jeq ; do i=is,ie ; if (seek_y_cor(i,J)) then
          T_int_S(i,J) = T_top(i,j) ; T_int_N(i,J) = T_top(i,j+1)
          S_int_S(i,J) = S_top(i,j) ; S_int_N(i,J) = S_top(i,j+1)
          p_int_S(i,J) = p(i,j,1) ; p_int_N(i,J) = p(i,j+1,1)
          inty_za_nonlin(i,J) = inty_za(i,J,1) - 0.5*(za(i,j,1) + za(i,j+1,1))
          dp_int_y(i,J) = p(i,j+1,1) - p(i,j,1)
          seek_y_cor(i,J) = .false.
        endif ; enddo ; enddo
      endif
    endif

    do J=Jsq,Jeq
      do i=is,ie
        ! This expression assumes that temperature and salinity vary linearly with pressure
        ! between the corners of the reference interfaces found above to get a correction to
        ! intx_pa that takes nonlinearities in the equation of state into account.
        ! It is derived from a 5 point quadrature estimate of the integral with a large-scale
        ! linear correction so that the pressures and heights match at the end-points.  It turns
        ! out that this linear correction cancels out the mid-point specific volume.
        ! This can be used without masking because dp_int_x and intx_za_nonlin are 0 over land.
        T5(1) = T_Int_S(i,J) ; S5(1) = S_Int_S(i,J) ; p5(1) = p_Int_S(i,J)
        T5(5) = T_Int_N(i,J) ; S5(5) = S_Int_N(i,J) ; p5(5) = p_Int_N(i,J)
        T5(2) = 0.25*(3.0*T5(1) + T5(5)) ; T5(4) = 0.25*(3.0*T5(5) + T5(1)) ; T5(3) = 0.5*(T5(5) + T5(1))
        S5(2) = 0.25*(3.0*S5(1) + S5(5)) ; S5(4) = 0.25*(3.0*S5(5) + S5(1)) ; S5(3) = 0.5*(S5(5) + S5(1))
        p5(2) = 0.25*(3.0*p5(1) + p5(5)) ; p5(4) = 0.25*(3.0*p5(5) + p5(1)) ; p5(3) = 0.5*(p5(5) + p5(1))
        call calculate_spec_vol(T5, S5, p5, SpV5, tv%eqn_of_state, spv_ref=alpha_ref)

        ! Note the consistency with the linear form below because (4.75 + 5.5/2) / 90 = 1/12
        inty_za_cor_ri(i,J) = C1_90 * (4.75*(SpV5(5)-SpV5(1)) + 5.5*(SpV5(4)-SpV5(2))) * &
                                      dp_int_y(i,J) - inty_za_nonlin(i,J)
      enddo
    enddo

    if (CS%debug) then
      call uvchksum("Pre-reset int[xy]_za", intx_za, inty_za, G%HI, haloshift=0, &
                    symmetric=G%Domain%symmetric, scalar_pair=.true., unscale=US%L_T_to_m_s**2)
      call uvchksum("int[xy]_za_cor", intx_za_cor_ri, inty_za_cor_ri, G%HI, haloshift=0, &
                    symmetric=G%Domain%symmetric, scalar_pair=.true., unscale=US%L_T_to_m_s**2)
      call uvchksum("int[xy]_za_nonlin", intx_za_nonlin, inty_za_nonlin, G%HI, haloshift=0, &
                    symmetric=G%Domain%symmetric, scalar_pair=.true., unscale=US%L_T_to_m_s**2)
      call uvchksum("dp_int_[xy]", dp_int_x, dp_int_y, G%HI, haloshift=0, &
                    symmetric=G%Domain%symmetric, unscale=US%RL2_T2_to_Pa)
    endif

    ! Correct intx_pa and inty_pa at each interface using vertically constant corrections.
    do K=1,nz+1 ; do j=js,je ; do I=Isq,Ieq
      intx_za(I,j,K) = intx_za(I,j,K) + intx_za_cor_ri(I,j)
    enddo ; enddo ; enddo

    do K=1,nz+1 ; do J=Jsq,Jeq ; do i=is,ie
      inty_za(i,J,K) = inty_za(i,J,K) + inty_za_cor_ri(i,J)
    enddo ; enddo ; enddo

    if (CS%debug) then
      call uvchksum("Post-reset int[xy]_za", intx_za, inty_za, G%HI, haloshift=0, &
                    symmetric=G%Domain%symmetric, scalar_pair=.true., unscale=US%L_T_to_m_s**2)
    endif

  endif ! intx_za and inty_za have now been reset to reflect the properties of an unimpeded interface.

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

  if (CS%id_MassWt_u>0) call post_data(CS%id_MassWt_u, MassWt_u, CS%diag)
  if (CS%id_MassWt_v>0) call post_data(CS%id_MassWt_v, MassWt_v, CS%diag)

  ! Diagnostics for tidal forcing and SAL height anomaly
  if (CS%id_e_tide>0) then
    ! To be consistent with old runs, tidal forcing diagnostic also includes total SAL.
    ! New diagnostics are given for each individual field.
    if (CS%tides_answer_date>20230630) then ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e_sal_and_tide(i,j) = e_sal(i,j) + e_tidal_eq(i,j) + e_tidal_sal(i,j)
    enddo ; enddo ; endif
    call post_data(CS%id_e_tide, e_sal_and_tide, CS%diag)
  endif
  if (CS%id_e_sal>0) call post_data(CS%id_e_sal, e_sal, CS%diag)
  if (CS%id_e_tidal_eq>0) call post_data(CS%id_e_tidal_eq, e_tidal_eq, CS%diag)
  if (CS%id_e_tidal_sal>0) call post_data(CS%id_e_tidal_sal, e_tidal_sal, CS%diag)

  ! Diagnostics for tidal forcing and SAL horizontal gradients
  if (CS%calculate_SAL .and. (associated(ADp%sal_u) .or. associated(ADp%sal_v))) then
    if (CS%tides) then ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e_sal(i,j) = e_sal(i,j) + e_tidal_sal(i,j)
    enddo ; enddo ; endif
    if (associated(ADp%sal_u)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      ADp%sal_u(I,j,k) = (e_sal(i+1,j) - e_sal(i,j)) * GV%g_Earth * G%IdxCu(I,j)
    enddo ; enddo ; enddo ; endif
    if (associated(ADp%sal_v)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      ADp%sal_v(i,J,k) = (e_sal(i,j+1) - e_sal(i,j)) * GV%g_Earth * G%IdyCv(i,J)
    enddo ; enddo ; enddo ; endif
    if (CS%id_sal_u>0) call post_data(CS%id_sal_u, ADp%sal_u, CS%diag)
    if (CS%id_sal_v>0) call post_data(CS%id_sal_v, ADp%sal_v, CS%diag)
  endif

  if (CS%tides .and. (associated(ADp%tides_u) .or. associated(ADp%tides_v))) then
    if (associated(ADp%tides_u)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      ADp%tides_u(I,j,k) = (e_tidal_eq(i+1,j) - e_tidal_eq(i,j)) * GV%g_Earth * G%IdxCu(I,j)
    enddo ; enddo ; enddo ; endif
    if (associated(ADp%tides_v)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      ADp%tides_v(i,J,k) = (e_tidal_eq(i,j+1) - e_tidal_eq(i,j)) * GV%g_Earth * G%IdyCv(i,J)
    enddo ; enddo ; enddo ; endif
    if (CS%id_tides_u>0) call post_data(CS%id_tides_u, ADp%tides_u, CS%diag)
    if (CS%id_tides_v>0) call post_data(CS%id_tides_v, ADp%tides_v, CS%diag)
  endif
end subroutine PressureForce_FV_nonBouss

!> \brief Boussinesq analytically-integrated finite volume form of pressure gradient
!!
!! Determines the acceleration due to hydrostatic pressure forces, using
!! the finite volume form of the terms and analytic integrals in depth.
!!
!! To work, the following fields must be set outside of the usual (is:ie,js:je)
!! range before this subroutine is called:
!!   h(isB:ie+1,jsB:je+1), T(isB:ie+1,jsB:je+1), and S(isB:ie+1,jsB:je+1).
subroutine PressureForce_FV_Bouss(h, tv, PFu, PFv, G, GV, US, CS, ALE_CSp, ADp, p_atm, pbce, eta)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thickness [H ~> m]
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: PFu !< Zonal acceleration [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: PFv !< Meridional acceleration [L T-2 ~> m s-2]
  type(PressureForce_FV_CS),                  intent(in)  :: CS  !< Finite volume PGF control structure
  type(ALE_CS),                               pointer     :: ALE_CSp !< ALE control structure
  type(accel_diag_ptrs),                      pointer     :: ADp !< Acceleration diagnostic pointers
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
    e_sal_and_tide, & ! The summation of self-attraction and loading and tidal forcing [Z ~> m].
    e_sal, &      ! The bottom geopotential anomaly due to self-attraction and loading [Z ~> m].
    e_tidal_eq,  & ! The bottom geopotential anomaly due to tidal forces from astronomical sources
                  ! [Z ~> m].
    e_tidal_sal, & ! The bottom geopotential anomaly due to harmonic self-attraction and loading
                  ! specific to tides [Z ~> m].
    Z_0p, &       ! The height at which the pressure used in the equation of state is 0 [Z ~> m]
    SSH, &      ! Sea surfae height anomaly for self-attraction and loading. Used if
                ! CALCULATE_SAL is True and SAL_USE_BPA is False [Z ~> m].
    pbot, &     ! Total bottom pressure for self-attraction and loading. Used if
                ! CALCULATE_SAL is True and SAL_USE_BPA is True [R L2 T-2 ~> Pa].
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
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    intx_pa_cor ! Correction for curvature in intx_pa [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    inty_pa_cor ! Correction for curvature in inty_pa [R L2 T-2 ~> Pa]

  ! These variables are used with reset_intxpa_integral.  The values are taken from different
  ! interfaces as a function of position.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    T_int_W, T_int_E, & ! Temperatures on the reference interface to the east and west of a u-point [C ~> degC]
    S_int_W, S_int_E, & ! Salinities on the reference interface to the east and west of a u-point [S ~> ppt]
    p_int_W, p_int_E, & ! Pressures on the reference interface to the east and west of a u-point [R L2 T-2 ~> Pa]
    rho_x_W, rho_x_E, & ! Density anomalies on the reference interface to the east and west
                        ! of a u-point [R ~> kg m-3]
    intx_pa_nonlin, &   ! Deviations in the previous version of intx_pa for the reference interface
                        ! from the value that would be obtained from assuming that pressure varies
                        ! linearly with depth along that interface [R L2 T-2 ~> Pa].
    dgeo_x, &           ! The change in x in geopotenial height along the reference interface [L2 T-2 ~> m2 s-2]
    intx_pa_cor_ri      ! The correction to intx_pa based on the reference interface calculations [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    T_int_S, T_int_N, & ! Temperatures on the reference interface to the north and south of a v-point [C ~> degC]
    S_int_S, S_int_N, & ! Salinities on the reference interface to the north and south of a v-point [S ~> ppt]
    p_int_S, p_int_N, & ! Pressures on the reference interface to the north and south of a v-point [R L2 T-2 ~> Pa]
    rho_y_S, rho_y_N, & ! Density anomalies on the reference interface to the north and south
                        ! of a v-point [R ~> kg m-3]
    inty_pa_nonlin, &   ! Deviations in the previous version of intx_pa for the reference interface
                        ! from the value that would be obtained from assuming that pressure varies
                        ! linearly with depth along that interface [R L2 T-2 ~> Pa].
    dgeo_y, &           ! The change in y in geopotenial height along the reference interface [L2 T-2 ~> m2 s-2]
    inty_pa_cor_ri      ! The correction to inty_pa based on the reference interface calculations [R L2 T-2 ~> Pa]
  logical, dimension(SZIB_(G),SZJ_(G)) :: &
    seek_x_cor          ! If true, try to find a u-point interface that would provide a better estimate
                        ! of the curvature terms in the intx_pa.
  logical, dimension(SZI_(G),SZJB_(G)) :: &
    seek_y_cor          ! If true, try to find a v-point interface that would provide a better estimate
                        ! of the curvature terms in the inty_pa.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    delta_z_x           ! If using flattest interface for reset integral, store x interface differences [Z ~> m]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    delta_z_y           ! If using flattest interface for reset integral, store y interface differences [Z ~> m]

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
  real, dimension(SZI_(G),SZJ_(G)) :: &
    T_top, &    ! Temperature of top layer used with correction_intxpa [C ~> degC]
    S_top, &    ! Salinity of top layer used with correction_intxpa [S ~> ppt]
    rho_top     ! Density anomaly of top layer used in calculating intx_pa_cor and inty_pa_cor
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    rho_pgf, rho_stanley_pgf ! Density [R ~> kg m-3] from EOS with and without SGS T variance
                             ! in Stanley parameterization.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    p_stanley   ! Pressure [R L2 T-2 ~> Pa] estimated with Rho_0
  real :: zeros(SZI_(G))     ! An array of zero values that can be used as an argument [various]
  real :: rho_in_situ(SZI_(G)) ! The in situ density [R ~> kg m-3].
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density, [R L2 T-2 ~> Pa] (usually 2e7 Pa = 2000 dbar).
  real :: p_surf_EOS(SZI_(G))  ! The pressure at the ocean surface determined from the surface height,
                             ! consistent with what is used in the density integral routines [R L2 T-2 ~> Pa]
  real :: p0(SZI_(G))        ! An array of zeros to use for pressure [R L2 T-2 ~> Pa].
  real :: dz_geo_sfc         ! The change in surface geopotential height between adjacent cells [L2 T-2 ~> m2 s-2]
  real :: GxRho0             ! The gravitational acceleration times mean ocean density [R L2 Z-1 T-2 ~> Pa m-1]
  real :: GxRho_ref          ! The gravitational acceleration times reference density [R L2 Z-1 T-2 ~> Pa m-1]
  real :: rho0_int_density   ! Rho0 used in int_density_dz_* subroutines [R ~> kg m-3]
  real :: rho0_set_pbce      ! Rho0 used in set_pbce_Bouss subroutine [R ~> kg m-3]
  real :: h_neglect          ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [H ~> m].
  real :: I_Rho0             ! The inverse of the Boussinesq reference density [R-1 ~> m3 kg-1].
  real :: G_Rho0             ! G_Earth / Rho_0 in [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1].
  real :: I_g_rho            ! The inverse of the density times the gravitational acceleration [Z T2 L-2 R-1 ~> m Pa-1]
  real :: rho_ref            ! The reference density [R ~> kg m-3].
  real :: dz_neglect         ! A minimal thickness [Z ~> m], like e.
  real :: dz_nonvanished     ! A small thickness considered to be vanished for mass weighting [Z ~> m]
  real :: H_to_RL2_T2        ! A factor to convert from thickness units (H) to pressure
                             ! units [R L2 T-2 H-1 ~> Pa m-1 or Pa m2 kg-1].
  real :: T5(5)         ! Temperatures and salinities at five quadrature points [C ~> degC]
  real :: S5(5)         ! Salinities at five quadrature points [S ~> ppt]
  real :: p5(5)         ! Full pressures at five quadrature points for use with the equation of state [R L2 T-2 ~> Pa]
  real :: pa5(5)        ! The pressure anomaly (i.e. pressure + g*RHO_0*e) at five quadrature points [R L2 T-2 ~> Pa].
  real :: r5(5)         ! Densities at five quadrature points [R ~> kg m-3]
  real :: wt_R          ! A weighting factor [nondim]
  real, parameter :: C1_6 = 1.0/6.0    ! A rational constant [nondim]
  real, parameter :: C1_90 = 1.0/90.0  ! A rational constant [nondim]
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_ALE         ! If true, use an ALE pressure reconstruction.
  logical :: use_EOS         ! If true, density is calculated from T & S using an equation of state.
  logical :: do_more_k       ! If true, there are still points where a flatter interface remains to be found.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer, dimension(2) :: EOSdom_h ! The i-computational domain for the equation of state at tracer points
  integer, dimension(2) :: EOSdom_u ! The i-computational domain for the equation of state at u-velocity points
  integer, dimension(2) :: EOSdom_v ! The i-computational domain for the equation of state at v-velocity points
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer :: i, j, k, m, k2

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  nkmb=GV%nk_rho_varies
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  EOSdom(1) = Isq - (G%isd-1) ;  EOSdom(2) = G%iec+1 - (G%isd-1)
  EOSdom_u(1) = Isq - (G%IsdB-1) ; EOSdom_u(2) = Ieq - (G%IsdB-1)
  EOSdom_v(1) = is - (G%isd-1)   ; EOSdom_v(2) = ie - (G%isd-1)

  if (.not.CS%initialized) call MOM_error(FATAL, &
       "MOM_PressureForce_FV_Bouss: Module must be initialized before it is used.")

  use_p_atm = associated(p_atm)
  use_EOS = associated(tv%eqn_of_state)
  do i=Isq,Ieq+1 ; p0(i) = 0.0 ; enddo
  use_ALE = .false.
  if (associated(ALE_CSp)) use_ALE = CS%reconstruct .and. use_EOS

  h_neglect = GV%H_subroundoff
  dz_neglect = GV%dZ_subroundoff
  dz_nonvanished = GV%H_to_Z*CS%h_nonvanished
  I_Rho0 = 1.0 / GV%Rho0
  G_Rho0 = GV%g_Earth / GV%Rho0
  GxRho0 = GV%g_Earth * GV%Rho0
  rho_ref = CS%rho_ref

  if (CS%rho_ref_bug) then
    rho0_int_density = rho_ref
    rho0_set_pbce = rho_ref
    GxRho_ref = GxRho0
    I_g_rho = 1.0 / (rho_ref * GV%g_Earth)
  else
    rho0_int_density = GV%Rho0
    rho0_set_pbce = GV%Rho0
    GxRho_ref = GV%g_Earth * rho_ref
    I_g_rho = 1.0 / (GV%rho0 * GV%g_Earth)
  endif

  if ((CS%id_MassWt_u > 0) .or. (CS%id_MassWt_v > 0)) then
    MassWt_u(:,:,:) = 0.0 ; MassWt_v(:,:,:) = 0.0
  endif

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    e(i,j,nz+1) = -G%bathyT(i,j)
  enddo ; enddo

  ! The following two if-blocks are used to recover old answers for self-attraction and loading
  ! (SAL) and tides only. The old algorithm moves interface heights before density calculations,
  ! and therefore is incorrect without SSH_IN_EOS_PRESSURE_FOR_PGF=True (added in August 2024).
  ! See the code right after Pa calculation loop for the new algorithm.

  ! Calculate and add SAL geopotential anomaly to interface height (old answers)
  if (CS%calculate_SAL .and. CS%tides_answer_date<=20250131) then
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

    if (CS%tides_answer_date>20230630) then ! answers_date between [20230701, 20250131]
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        e(i,j,nz+1) = e(i,j,nz+1) - e_sal(i,j)
      enddo ; enddo
    endif
  endif

  ! Calculate and add tidal geopotential anomaly to interface height (old answers)
  if (CS%tides .and. CS%tides_answer_date<=20250131) then
    if (CS%tides_answer_date>20230630) then ! answers_date between [20230701, 20250131]
      call calc_tidal_forcing(CS%Time, e_tidal_eq, e_tidal_sal, G, US, CS%tides_CSp)
     !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        e(i,j,nz+1) = e(i,j,nz+1) - (e_tidal_eq(i,j) + e_tidal_sal(i,j))
      enddo ; enddo
    else  ! answers_date before 20230701
      if (.not.CS%calculate_SAL) e_sal(:,:) = 0.0
      call calc_tidal_forcing_legacy(CS%Time, e_sal, e_sal_and_tide, e_tidal_eq, e_tidal_sal, &
                                     G, US, CS%tides_CSp)
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        e(i,j,nz+1) = e(i,j,nz+1) - e_sal_and_tide(i,j)
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
  if ( use_ALE .and. (CS%Recon_Scheme == 1) ) then
    call TS_PLM_edge_values(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, CS%boundary_extrap)
  elseif ( use_ALE .and. (CS%Recon_Scheme == 2) ) then
    call TS_PPM_edge_values(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, CS%boundary_extrap)
  elseif (CS%reset_intxpa_integral) then
    do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      T_b(i,j,k) = tv%T(i,j,k) ; S_b(i,j,k) = tv%S(i,j,k)
    enddo ; enddo ; enddo
  endif

  ! Set the surface boundary conditions on pressure anomaly and its horizontal
  ! integrals, assuming that the surface pressure anomaly varies linearly
  ! in x and y.
  if (use_p_atm) then
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      pa(i,j,1) = GxRho_ref * (e(i,j,1) - G%Z_ref) + p_atm(i,j)
    enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      pa(i,j,1) = GxRho_ref * (e(i,j,1) - G%Z_ref)
    enddo ; enddo
  endif

  if (CS%use_SSH_in_Z0p .and. use_p_atm) then
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
                    rho_ref, rho0_int_density, GV%g_Earth, dz_neglect, G%bathyT, &
                    G%HI, GV, tv%eqn_of_state, US, CS%use_stanley_pgf, dpa(:,:,k), intz_dpa(:,:,k), &
                    intx_dpa(:,:,k), inty_dpa(:,:,k), &
                    MassWghtInterp=CS%MassWghtInterp, &
                    use_inaccurate_form=CS%use_inaccurate_pgf_rho_anom, Z_0p=Z_0p, &
                    MassWghtInterpVanOnly=CS%MassWghtInterpVanOnly, h_nv=dz_nonvanished)
        elseif ( CS%Recon_Scheme == 2 ) then
          call int_density_dz_generic_ppm(k, tv, T_t, T_b, S_t, S_b, e, &
                    rho_ref, rho0_int_density, GV%g_Earth, dz_neglect, G%bathyT, &
                    G%HI, GV, tv%eqn_of_state, US, CS%use_stanley_pgf, dpa(:,:,k), intz_dpa(:,:,k), &
                    intx_dpa(:,:,k), inty_dpa(:,:,k), &
                    MassWghtInterp=CS%MassWghtInterp, Z_0p=Z_0p, &
                    MassWghtInterpVanOnly=CS%MassWghtInterpVanOnly, h_nv=dz_nonvanished)
        endif
      else
        call int_density_dz(tv_tmp%T(:,:,k), tv_tmp%S(:,:,k), e(:,:,K), e(:,:,K+1), &
                  rho_ref, rho0_int_density, GV%g_Earth, G%HI, tv%eqn_of_state, US, dpa(:,:,k), &
                  intz_dpa(:,:,k), intx_dpa(:,:,k), inty_dpa(:,:,k), G%bathyT, e(:,:,1), dz_neglect, &
                  CS%MassWghtInterp, Z_0p=Z_0p, &
                  MassWghtInterpVanOnly=CS%MassWghtInterpVanOnly, h_nv=dz_nonvanished)
      endif
      if (GV%Z_to_H /= 1.0) then
        !$OMP parallel do default(shared)
        do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
          intz_dpa(i,j,k) = intz_dpa(i,j,k)*GV%Z_to_H
        enddo ; enddo
      endif
      if ((CS%id_MassWt_u > 0) .or. (CS%id_MassWt_v > 0)) &
        call diagnose_mass_weight_Z(e(:,:,K), e(:,:,K+1), G%bathyT, e(:,:,1), dz_neglect, CS%MassWghtInterp, &
                                    G%HI, MassWt_u(:,:,k), MassWt_v(:,:,k), &
                                    MassWghtInterpVanOnly=CS%MassWghtInterpVanOnly, h_nv=CS%h_nonvanished)
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

  ! Calculate and add SAL geopotential anomaly to interface height (new answers)
  if (CS%calculate_SAL .and. CS%tides_answer_date>20250131) then
    if (CS%sal_use_bpa) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        pbot(i,j) = pa(i,j,nz+1) - GxRho_ref * (e(i,j,nz+1) - G%Z_ref)
      enddo ; enddo
      call calc_SAL(pbot, e_sal, G, CS%SAL_CSp, tmp_scale=US%Z_to_m)
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        SSH(i,j) = e(i,j,1) - max(-G%bathyT(i,j) - G%Z_ref, 0.0) ! Remove topography above sea level
      enddo ; enddo
      call calc_SAL(SSH, e_sal, G, CS%SAL_CSp, tmp_scale=US%Z_to_m)
    endif
    if (.not.CS%bq_sal_tides) then ; do K=1,nz+1
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        e(i,j,K) = e(i,j,K) - e_sal(i,j)
        pa(i,j,K) = pa(i,j,K) - GxRho_ref * e_sal(i,j)
      enddo ; enddo
    enddo ; endif
  endif

  ! Calculate and add tidal geopotential anomaly to interface height (new answers)
  if (CS%tides .and. CS%tides_answer_date>20250131) then
    call calc_tidal_forcing(CS%Time, e_tidal_eq, e_tidal_sal, G, US, CS%tides_CSp)
    if (.not.CS%bq_sal_tides) then ; do K=1,nz+1
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        e(i,j,K) = e(i,j,K) - (e_tidal_eq(i,j) + e_tidal_sal(i,j))
        pa(i,j,K) = pa(i,j,K) - GxRho_ref * (e_tidal_eq(i,j) + e_tidal_sal(i,j))
      enddo ; enddo
    enddo ; endif
  endif

  if (CS%correction_intxpa .or. CS%reset_intxpa_integral) then
    ! Determine surface temperature and salinity for use in the pressure gradient corrections
    if (use_ALE .and. (CS%Recon_Scheme > 0)) then
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        T_top(i,j) = T_t(i,j,1) ; S_top(i,j) = S_t(i,j,1)
      enddo ; enddo
    else
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        T_top(i,j) = tv%T(i,j,1) ; S_top(i,j) = tv%S(i,j,1)
      enddo ; enddo
    endif
  endif

  if (CS%correction_intxpa) then
    ! Determine surface density for use in the pressure gradient corrections
    !$OMP parallel do default(shared) private(p_surf_EOS)
    do j=Jsq,Jeq+1
      ! P_surf_EOS here is consistent with the pressure that is used in the int_density_dz routines.
      do i=Isq,Ieq+1 ; p_surf_EOS(i) = -GxRho0*(e(i,j,1) - Z_0p(i,j)) ; enddo
      call calculate_density(T_top(:,j), S_top(:,j), p_surf_EOS, rho_top(:,j), &
                             tv%eqn_of_state, EOSdom, rho_ref=rho_ref)
    enddo

    if (CS%debug) then
      call hchksum(rho_top, "intx_pa rho_top", G%HI, haloshift=1, unscale=US%R_to_kg_m3)
      call hchksum(e(:,:,1), "intx_pa e(1)", G%HI, haloshift=1, unscale=US%Z_to_m)
      call hchksum(pa(:,:,1), "intx_pa pa(1)", G%HI, haloshift=1, unscale=US%RL2_T2_to_Pa)
    endif

    ! This version attempts to correct for hydrostatic variations in surface pressure under ice.
    !$OMP parallel do default(shared) private(dz_geo_sfc)
    do j=js,je ; do I=Isq,Ieq
      intx_pa_cor(I,j) = 0.0
      dz_geo_sfc = GV%g_Earth * (e(i+1,j,1)-e(i,j,1))
      if ((dz_geo_sfc * rho_ref - (pa(i+1,j,1)-pa(i,j,1)))*dz_geo_sfc > 0.0) then
        ! The pressure/depth relationship has a positive implied density given by
        !   rho_implied = rho_ref - (pa(i+1,j,1)-pa(i,j,1)) / dz_geo_sfc
        if (-dz_geo_sfc * (pa(i+1,j,1)-pa(i,j,1)) > &
            0.25*((rho_top(i+1,j)+rho_top(i,j))-2.0*rho_ref) * dz_geo_sfc**2) then
          ! The pressure difference is at least half the size of the difference expected by hydrostatic
          ! balance.  This test gets rid of pressure differences that are small, e.g. open ocean.
          ! Use 5 point quadrature to calculate intxpa
          T5(1) = T_top(i,j) ; T5(5) = T_top(i+1,j)
          S5(1) = S_top(i,j) ; S5(5) = S_top(i+1,j)
          pa5(1) = pa(i,j,1) ; pa5(5) = pa(i+1,j,1)
          ! Pressure input to density EOS is consistent with the pressure used in the int_density_dz routines.
          p5(1) = -GxRho0*(e(i,j,1) - Z_0p(i,j))
          p5(5) = -GxRho0*(e(i+1,j,1) - Z_0p(i,j))
          do m=2,4
            wt_R =  0.25*real(m-1)
            T5(m) = T5(1) + (T5(5)-T5(1))*wt_R
            S5(m) = S5(1) + (S5(5)-S5(1))*wt_R
            p5(m) = p5(1) + (p5(5)-p5(1))*wt_R
          enddo !m
          call calculate_density(T5, S5, p5, r5, tv%eqn_of_state, rho_ref=rho_ref)

          ! Use a trapezoidal rule integral of the hydrostatic equation to determine the pressure
          ! anomalies at 5 equally spaced points along the interface, and then use Boole's rule
          ! quadrature to find the integrated correction to the integral of pressure along the interface.
          ! The derivation for this expression is shown below in the y-direction version.
          intx_pa_cor(I,j) = C1_90 * (4.75*(r5(5)-r5(1)) + 5.5*(r5(4)-r5(2))) * dz_geo_sfc
          ! Note that (4.75 + 5.5/2) / 90 = 1/12, so this is consistent with the linear result below.
        endif
      endif
      intx_pa(I,j,1) = 0.5*(pa(i,j,1) + pa(i+1,j,1)) + intx_pa_cor(I,j)
    enddo ; enddo
    !$OMP parallel do default(shared) private(dz_geo_sfc)
    do J=Jsq,Jeq ; do i=is,ie
      inty_pa_cor(i,J) = 0.0
      dz_geo_sfc = GV%g_Earth * (e(i,j+1,1)-e(i,j,1))
      if ((dz_geo_sfc * rho_ref - (pa(i,j+1,1)-pa(i,j,1)))*dz_geo_sfc > 0.0) then
        ! The pressure/depth relationship has a positive implied density
        if (-dz_geo_sfc * (pa(i,j+1,1)-pa(i,j,1)) > &
            0.25*((rho_top(i,j+1)+rho_top(i,j))-2.0*rho_ref) * dz_geo_sfc**2) then
          ! The pressure difference is at least half the size of the difference expected by hydrostatic
          ! balance.  This test gets rid of pressure differences that are small, e.g. open ocean.
          ! Use 5 point quadrature to calculate intypa
          T5(1) = T_top(i,j) ; T5(5) = T_top(i,j+1)
          S5(1) = S_top(i,j) ; S5(5) = S_top(i,j+1)
          pa5(1) = pa(i,j,1) ; pa5(5) = pa(i,j+1,1)
          ! Pressure input to density EOS is consistent with the pressure used in the int_density_dz routines.
          p5(1) = -GxRho0*(e(i,j,1) - Z_0p(i,j))
          p5(5) = -GxRho0*(e(i,j+1,1) - Z_0p(i,j))

          do m=2,4
            wt_R =  0.25*real(m-1)
            T5(m) = T5(1) + (T5(5)-T5(1))*wt_R
            S5(m) = S5(1) + (S5(5)-S5(1))*wt_R
            p5(m) = p5(1) + (p5(5)-p5(1))*wt_R
          enddo !m
          call calculate_density(T5, S5, p5, r5, tv%eqn_of_state, rho_ref=rho_ref)

          ! Use a trapezoidal rule integral of the hydrostatic equation to determine the pressure
          ! anomalies at 5 equally spaced points along the interface, and then use Boole's rule
          ! quadrature to find the integrated correction to the integral of pressure along the interface.
          inty_pa_cor(i,J) = C1_90 * (4.75*(r5(5)-r5(1)) + 5.5*(r5(4)-r5(2))) * dz_geo_sfc

          ! The derivation of this correction follows:

          ! Make pressure curvature a difference from the linear fit of pressure between the two points
          ! (which is equivalent to taking 4 trapezoidal rule integrals of the hydrostatic equation on
          ! sub-segments), with a constant slope that is chosen so that the pressure anomalies at the
          ! two ends of the segment agree with their known values.
          ! d_geo_8 = 0.125*dz_geo_sfc
          ! dpa_subseg = 0.25*(pa5(5)-pa5(1)) + &
          !              0.25*d_geo_8 * ((r5(5)+r5(1)) + 2.0*((r5(4)+r5(2)) + r5(3)))
          ! do m=2,4
          !   pa5(m) = pa5(m-1) + dpa_subseg - d_geo_8*(r5(m)+r5(m-1)))
          ! enddo

          ! Explicitly finding expressions for the incremental pressures from the recursion relation above:
          ! pa5(2) = 0.25*(3.*pa5(1) + pa5(5)) + 0.25*d_geo_8 * ( (r5(5)-3.*r5(1)) + 2.0*((r5(4)-r5(2)) + r5(3)) )
          ! ! pa5(3) = 0.5*(pa5(1) + pa5(5)) + 0.25*d_geo_8 * &
          ! !   ( (r5(5)+r5(1)) + 2.0*((r5(4)+r5(2)) + r5(3)) + &
          ! !     (r5(5)-3.*r5(1)) + 2.0*((r5(4)-r5(2)) + r5(3)) - 4.*(r5(3)+r5(2)) )
          ! pa5(3) = 0.5*(pa5(1) + pa5(5)) + d_geo_8 * (0.5*(r5(5)-r5(1)) + (r5(4)-r5(2)) )
          ! ! pa5(4) = 0.25*(pa5(1) + 3.0*pa5(5)) + 0.25*d_geo_8 * &
          ! !   (2.0*(r5(5)-r5(1)) + 4.0*(r5(4)-r5(2)) + (r5(5)+r5(1)) + &
          ! !    2.0*(r5(4)+r5(2)) + 2.0*r5(3) - 4.*(r5(4)+r5(3)))
          ! pa5(4) = 0.25*(pa5(1) + 3.0*pa5(5)) + 0.25*d_geo_8 * ( (3.*r5(5)-r5(1)) + 2.0*((r5(4)-r5(2)) - r5(3)) )
          ! ! pa5(5) = pa5(5) + 0.25*d_geo_8 * &
          ! !     ( (3.*r5(5)-r5(1)) + 2.0*((r5(4)-r5(2)) - r5(3)) + &
          ! !      ((r5(5)+r5(1)) + 2.0*((r5(4)+r5(2)) + r5(3))) - 4.*(r5(5)+r5(4)) )
          ! pa5(5) = pa5(5)  ! As it should.

          ! From these:
          ! pa5(2) + pa5(4) = (pa5(1) + pa5(5)) + 0.25*d_geo_8 * &
          !     ( (r5(5)-3.*r5(1)) + 2.0*((r5(4)-r5(2)) + r5(3)) + (3.*r5(5)-r5(1)) + 2.0*((r5(4)-r5(2)) - r5(3))
          ! pa5(2) + pa5(4) = (pa5(1) + pa5(5)) + d_geo_8 * ( (r5(5)-r5(1)) + (r5(4)-r5(2)) )

          ! Get the correction from the difference between the 5-point quadrature integral of pa5 and
          ! its trapezoidal rule integral as:
          ! inty_pa_cor(i,J) = C1_90*(7.0*(pa5(1)+pa5(5)) + 32.0*(pa5(2)+pa5(4)) + 12.0*pa5(3)) - 0.5*(pa5(1)+pa5(5)))
          ! inty_pa_cor(i,J) = C1_90*((32.0*(pa5(2)+pa5(4)) + 12.0*pa5(3)) - 38.0*(pa5(1)+pa5(5)))
          ! inty_pa_cor(i,J) = C1_90*d_geo_8 * ((32.0*( (r5(5)-r5(1)) + (r5(4)-r5(2)) ) + &
          !                                      (6.*(r5(5)-r5(1)) + 12.0*(r5(4)-r5(2)) ))
          ! inty_pa_cor(i,J) = C1_90*d_geo_8 * ( 38.0*(r5(5)-r5(1)) + 44.0*(r5(4)-r5(2)) )
        endif
      endif
      inty_pa(i,J,1) = 0.5*(pa(i,j,1) + pa(i,j+1,1)) + inty_pa_cor(i,J)
    enddo ; enddo

    if (CS%debug) then
      call uvchksum("int[xy]_pa_cor", intx_pa_cor, inty_pa_cor, G%HI, haloshift=0, &
                    symmetric=G%Domain%symmetric, scalar_pair=.true., unscale=US%RL2_T2_to_Pa)
      call uvchksum("int[xy]_pa(1)", intx_pa(:,:,1), inty_pa(:,:,1), G%HI, haloshift=0, &
                    symmetric=G%Domain%symmetric, scalar_pair=.true., unscale=US%RL2_T2_to_Pa)
    endif

  else
    ! Set the surface boundary conditions on the horizontally integrated pressure anomaly,
    ! assuming that the surface pressure anomaly varies linearly in x and y.
    ! If there is an ice-shelf or icebergs, this linear variation would need to be applied
    ! to an interior interface.
    !$OMP parallel do default(shared)
    do j=js,je ; do I=Isq,Ieq
      intx_pa(I,j,1) = 0.5*(pa(i,j,1) + pa(i+1,j,1))
    enddo ; enddo
    !$OMP parallel do default(shared)
    do J=Jsq,Jeq ; do i=is,ie
      inty_pa(i,J,1) = 0.5*(pa(i,j,1) + pa(i,j+1,1))
    enddo ; enddo
  endif

  do k=1,nz
    !$OMP parallel do default(shared)
    do j=js,je ; do I=Isq,Ieq
      intx_pa(I,j,K+1) = intx_pa(I,j,K) + intx_dpa(I,j,k)
    enddo ; enddo
  enddo
  do k=1,nz
    !$OMP parallel do default(shared)
    do J=Jsq,Jeq ; do i=is,ie
      inty_pa(i,J,K+1) = inty_pa(i,J,K) + inty_dpa(i,J,k)
    enddo ; enddo
  enddo

  if (CS%reset_intxpa_integral) then
    ! Having stored the pressure gradient info, we can work out where the first nonvanished layers is
    ! reset intxpa there, then adjust intxpa throughout the water column.

    ! Zero out the 2-d arrays that will be set from various reference interfaces.
    T_int_W(:,:) = 0.0 ; S_int_W(:,:) = 0.0 ; p_int_W(:,:) = 0.0
    T_int_E(:,:) = 0.0 ; S_int_E(:,:) = 0.0 ; p_int_E(:,:) = 0.0
    intx_pa_nonlin(:,:) = 0.0 ; dgeo_x(:,:) = 0.0 ; intx_pa_cor_ri(:,:) = 0.0
    do j=js,je ; do I=Isq,Ieq
      seek_x_cor(I,j) = (G%mask2dCu(I,j) > 0.)
      delta_z_x(I,j)  = 0.0
    enddo ; enddo

    do j=js,je ; do I=Isq,Ieq ; if (seek_x_cor(I,j)) then
      if ((e(i+1,j,2) <= e(i,j,1)) .and. (e(i,j,2) <= e(i+1,j,1))) then
        ! This is a typical case in the open ocean, so use the topmost interface.
        T_int_W(I,j) = T_top(i,j) ; T_int_E(I,j) = T_top(i+1,j)
        S_int_W(I,j) = S_top(i,j) ; S_int_E(I,j) = S_top(i+1,j)
        p_int_W(I,j) = -GxRho0*(e(i,j,1) - Z_0p(i,j))
        p_int_E(I,j) = -GxRho0*(e(i+1,j,1) - Z_0p(i,j))
        intx_pa_nonlin(I,j) = intx_pa(I,j,1) - 0.5*(pa(i,j,1) + pa(i+1,j,1))
        dgeo_x(I,j) = GV%g_Earth * (e(i+1,j,1)-e(i,j,1))
        seek_x_cor(I,j) = .false.
      endif
    endif ; enddo ; enddo

    do k=1,nz
      do_more_k = .false.
      do j=js,je ; do I=Isq,Ieq ; if (seek_x_cor(I,j)) then
        ! Find the topmost layer for which both sides are nonvanished and mass-weighting is not
        ! activated in the subgrid interpolation.
        if (((h(i,j,k) > CS%h_nonvanished) .and. (h(i+1,j,k) > CS%h_nonvanished)) .and. &
            (max(0., e(i+1,j,K+1)-e(i,j,1), e(i,j,K+1)-e(i+1,j,1)) <= 0.0)) then
          ! Store properties at the bottom of this cell to get a "good estimate" for intxpa at
          ! the interface below this cell (it might have quadratic pressure dependence if sloped)
          T_int_W(I,j) = T_b(i,j,k) ; T_int_E(I,j) = T_b(i+1,j,k)
          S_int_W(I,j) = S_b(i,j,k) ; S_int_E(I,j) = S_b(i+1,j,k)
          ! These pressures are only used for the equation of state, and are only a function of
          ! height, consistent with the expressions in the int_density_dz routines.
          p_int_W(I,j) = -GxRho0*(e(i,j,K+1) - Z_0p(i,j))
          p_int_E(I,j) = -GxRho0*(e(i+1,j,K+1) - Z_0p(i,j))

          intx_pa_nonlin(I,j) = intx_pa(I,j,K+1) - 0.5*(pa(i,j,K+1) + pa(i+1,j,K+1))
          dgeo_x(I,j) = GV%g_Earth * (e(i+1,j,K+1)-e(i,j,K+1))
          seek_x_cor(I,j) = .false.
        else
          do_more_k = .true.
        endif
      endif ; enddo ; enddo
      if (.not.do_more_k) exit  ! All reference interfaces have been found, so stop working downward.
    enddo

    if (do_more_k) then
      if (CS%reset_intxpa_flattest) then
      ! There are still points where a correction is needed, so use flattest interface
        do j=js,je ; do I=Isq,Ieq ; if (seek_x_cor(I,j)) then
          ! choose top layer first
          T_int_W(I,j) = T_top(i,j) ; T_int_E(I,j) = T_top(i+1,j)
          S_int_W(I,j) = S_top(i,j) ; S_int_E(I,j) = S_top(i+1,j)
          p_int_W(I,j) = -GxRho0*(e(i,j,1) - Z_0p(i,j))
          p_int_E(I,j) = -GxRho0*(e(i+1,j,1) - Z_0p(i,j))
          intx_pa_nonlin(I,j) = intx_pa(I,j,1) - 0.5*(pa(i,j,1) + pa(i+1,j,1))
          dgeo_x(I,j) = GV%g_Earth * (e(i+1,j,1)-e(i,j,1))
          delta_z_x(I,j) = abs(e(i+1,j,1)-e(i,j,1))
          do k=1,nz
            if (abs(e(i+1,j,k+1)-e(i,j,k+1)) < delta_z_x(I,j)) then
              ! bottom of layer is less sloped than top. Use this layer
              delta_z_x(I,j) = abs(e(i+1,j,k+1)-e(i,j,k+1))
              T_int_W(I,j) = T_b(i,j,k) ; T_int_E(I,j) = T_b(i+1,j,k)
              S_int_W(I,j) = S_b(i,j,k) ; S_int_E(I,j) = S_b(i+1,j,k)
              p_int_W(I,j) = -GxRho0*(e(i,j,K+1) - Z_0p(i,j))
              p_int_E(I,j) = -GxRho0*(e(i+1,j,K+1) - Z_0p(i,j))
              intx_pa_nonlin(I,j) = intx_pa(I,j,K+1) - 0.5*(pa(i,j,K+1) + pa(i+1,j,K+1))
              dgeo_x(I,j) = GV%g_Earth * (e(i+1,j,K+1)-e(i,j,K+1))
            endif
          enddo
          seek_x_cor(I,j) = .false.
        endif ; enddo ; enddo
      else
        ! There are still points where a correction is needed, so use the top interface for lack of a better idea?
        do j=js,je ; do I=Isq,Ieq ; if (seek_x_cor(I,j)) then
          T_int_W(I,j) = T_top(i,j) ; T_int_E(I,j) = T_top(i+1,j)
          S_int_W(I,j) = S_top(i,j) ; S_int_E(I,j) = S_top(i+1,j)
          p_int_W(I,j) = -GxRho0*(e(i,j,1) - Z_0p(i,j))
          p_int_E(I,j) = -GxRho0*(e(i+1,j,1) - Z_0p(i,j))
          intx_pa_nonlin(I,j) = intx_pa(I,j,1) - 0.5*(pa(i,j,1) + pa(i+1,j,1))
          dgeo_x(I,j) = GV%g_Earth * (e(i+1,j,1)-e(i,j,1))
          seek_x_cor(I,j) = .false.
        endif ; enddo ; enddo
      endif
    endif

    do j=js,je
      do I=Isq,Ieq
        ! This expression assumes that temperature and salinity vary linearly with hieght
        ! between the corners of the reference interfaces found above to get a correction to
        ! intx_pa that takes nonlinearities in the equation of state into account.
        ! It is derived from a 5 point quadrature estimate of the integral with a large-scale
        ! linear correction so that the pressures and heights match at the end-points.  It turns
        ! out that this linear correction cancels out the mid-point density anomaly.
        ! This can be used without masking because dgeo_x and intx_pa_nonlin are 0 over land.
        T5(1) = T_Int_W(I,j) ; S5(1) = S_Int_W(I,j) ; p5(1) = p_Int_W(I,j)
        T5(5) = T_Int_E(I,j) ; S5(5) = S_Int_E(I,j) ; p5(5) = p_Int_E(I,j)
        T5(2) = 0.25*(3.0*T5(1) + T5(5)) ; T5(4) = 0.25*(3.0*T5(5) + T5(1)) ; T5(3) = 0.5*(T5(5) + T5(1))
        S5(2) = 0.25*(3.0*S5(1) + S5(5)) ; S5(4) = 0.25*(3.0*S5(5) + S5(1)) ; S5(3) = 0.5*(S5(5) + S5(1))
        p5(2) = 0.25*(3.0*p5(1) + p5(5)) ; p5(4) = 0.25*(3.0*p5(5) + p5(1)) ; p5(3) = 0.5*(p5(5) + p5(1))
        call calculate_density(T5, S5, p5, r5, tv%eqn_of_state, rho_ref=rho_ref)

        ! Note the consistency with the linear form below because (4.75 + 5.5/2) / 90 = 1/12
        intx_pa_cor_ri(I,j) = C1_90 * (4.75*(r5(5)-r5(1)) + 5.5*(r5(4)-r5(2))) * dgeo_x(I,j) - &
                              intx_pa_nonlin(I,j)
      enddo
    enddo

    ! Repeat the calculations above for v-velocity points.
    T_int_S(:,:) = 0.0 ; S_int_S(:,:) = 0.0 ; p_int_S(:,:) = 0.0
    T_int_N(:,:) = 0.0 ; S_int_N(:,:) = 0.0 ; p_int_N(:,:) = 0.0
    inty_pa_nonlin(:,:) = 0.0 ; dgeo_y(:,:) = 0.0 ; inty_pa_cor_ri(:,:) = 0.0
    do J=Jsq,Jeq ; do i=is,ie
      seek_y_cor(i,J) = (G%mask2dCv(i,J) > 0.)
      delta_z_y(i,J)  = 0.0
    enddo ; enddo

    do J=Jsq,Jeq ; do i=is,ie ; if (seek_y_cor(i,J)) then
      if ((e(i,j+1,2) <= e(i,j,1)) .and. (e(i,j,2) <= e(i,j+1,1))) then
        ! This is a typical case in the open ocean, so use the topmost interface.
        T_int_S(i,J) = T_top(i,j) ; T_int_N(i,J) = T_top(i,j+1)
        S_int_S(i,J) = S_top(i,j) ; S_int_N(i,J) = S_top(i,j+1)
        p_int_S(i,J) = -GxRho0*(e(i,j,1) - Z_0p(i,j))
        p_int_N(i,J) = -GxRho0*(e(i,j+1,1) - Z_0p(i,j))
        inty_pa_nonlin(i,J) = inty_pa(i,J,1) - 0.5*(pa(i,j,1) + pa(i,j+1,1))
        dgeo_y(i,J) = GV%g_Earth * (e(i,j+1,1)-e(i,j,1))
        seek_y_cor(i,J) = .false.
      endif
    endif ; enddo ; enddo

    do k=1,nz
      do_more_k = .false.
      do J=Jsq,Jeq ; do i=is,ie ; if (seek_y_cor(i,J)) then
        ! Find the topmost layer for which both sides are nonvanished and mass-weighting is not
        ! activated in the subgrid interpolation.
        if (((h(i,j,k) > CS%h_nonvanished) .and. (h(i,j+1,k) > CS%h_nonvanished)) .and. &
            (max(0., e(i,j+1,K+1)-e(i,j,1), e(i,j,K+1)-e(i,j+1,1)) <= 0.0)) then
          ! Store properties at the bottom of this cell to get a "good estimate" for intypa at
          ! the interface below this cell (it might have quadratic pressure dependence if sloped)
          T_int_S(i,J) = T_b(i,j,k) ; T_int_N(i,J) = T_b(i,j+1,k)
          S_int_S(i,J) = S_b(i,j,k) ; S_int_N(i,J) = S_b(i,j+1,k)
          ! These pressures are only used for the equation of state, and are only a function of
          ! height, consistent with the expressions in the int_density_dz routines.
          p_int_S(i,J) = -GxRho0*(e(i,j,K+1) - Z_0p(i,j))
          p_int_N(i,J) = -GxRho0*(e(i,j+1,K+1) - Z_0p(i,j))
          inty_pa_nonlin(i,J) = inty_pa(i,J,K+1) - 0.5*(pa(i,j,K+1) + pa(i,j+1,K+1))
          dgeo_y(i,J) = GV%g_Earth * (e(i,j+1,K+1)-e(i,j,K+1))
          seek_y_cor(i,J) = .false.
        else
          do_more_k = .true.
        endif
      endif ; enddo ; enddo
      if (.not.do_more_k) exit  ! All reference interfaces have been found, so stop working downward.
    enddo

    if (do_more_k) then
      if (CS%reset_intxpa_flattest) then
        ! There are still points where a correction is needed, so use flattest interface.
        do J=Jsq,Jeq ; do i=is,ie ; if (seek_y_cor(i,J)) then
          ! choose top interface first
          T_int_S(i,J) = T_top(i,j) ; T_int_N(i,J) = T_top(i,j+1)
          S_int_S(i,J) = S_top(i,j) ; S_int_N(i,J) = S_top(i,j+1)
          p_int_S(i,J) = -GxRho0*(e(i,j,1) - Z_0p(i,j))
          p_int_N(i,J) = -GxRho0*(e(i,j+1,1) - Z_0p(i,j))
          inty_pa_nonlin(i,J) = inty_pa(i,J,1) - 0.5*(pa(i,j,1) + pa(i,j+1,1))
          dgeo_y(i,J) = GV%g_Earth * (e(i,j+1,1)-e(i,j,1))
          delta_z_y(i,J) = abs(e(i,j+1,1)-e(i,j,1))
          do k=1,nz
            if (abs(e(i,j+1,k+1)-e(i,j,k+1)) < delta_z_y(i,J)) then
              ! bottom of layer is less sloped than top. Use this layer
              delta_z_y(i,J) = abs(e(i,j+1,k+1)-e(i,j,k+1))
              T_int_S(i,J) = T_b(i,j,k) ; T_int_N(i,J) = T_b(i,j+1,k)
              S_int_S(i,J) = S_b(i,j,k) ; S_int_N(i,J) = S_b(i,j+1,k)
              p_int_S(i,J) = -GxRho0*(e(i,j,k+1) - Z_0p(i,j))
              p_int_N(i,J) = -GxRho0*(e(i,j+1,k+1) - Z_0p(i,j))
              inty_pa_nonlin(i,J) = inty_pa(i,J,k+1) - 0.5*(pa(i,j,k+1) + pa(i,j+1,k+1))
              dgeo_y(i,J) = GV%g_Earth * (e(i,j+1,k+1)-e(i,j,k+1))
            endif
          enddo
          seek_y_cor(i,J) = .false.
        endif ; enddo ; enddo
      else
        ! There are still points where a correction is needed, so use the top interface for lack of a better idea?
        do J=Jsq,Jeq ; do i=is,ie ; if (seek_y_cor(i,J)) then
          T_int_S(i,J) = T_top(i,j) ; T_int_N(i,J) = T_top(i,j+1)
          S_int_S(i,J) = S_top(i,j) ; S_int_N(i,J) = S_top(i,j+1)
          p_int_S(i,J) = -GxRho0*(e(i,j,1) - Z_0p(i,j))
          p_int_N(i,J) = -GxRho0*(e(i,j+1,1) - Z_0p(i,j))
          inty_pa_nonlin(i,J) = inty_pa(i,J,1) - 0.5*(pa(i,j,1) + pa(i,j+1,1))
          dgeo_y(i,J) = GV%g_Earth * (e(i,j+1,1)-e(i,j,1))
          seek_y_cor(i,J) = .false.
        endif ; enddo ; enddo
      endif
    endif

    do J=Jsq,Jeq
      do i=is,ie
        ! This expression assumes that temperature and salinity vary linearly with hieght
        ! between the corners of the reference interfaces found above to get a correction to
        ! intx_pa that takes nonlinearities in the equation of state into account.
        ! It is derived from a 5 point quadrature estimate of the integral with a large-scale
        ! linear correction so that the pressures and heights match at the end-points.  It turns
        ! out that this linear correction cancels out the mid-point density anomaly.
        ! This can be used without masking because dgeo_y and inty_pa_nonlin are 0 over land.
        T5(1) = T_Int_S(i,J) ; S5(1) = S_Int_S(i,J) ; p5(1) = p_Int_S(i,J)
        T5(5) = T_Int_N(i,J) ; S5(5) = S_Int_N(i,J) ; p5(5) = p_Int_N(i,J)
        T5(2) = 0.25*(3.0*T5(1) + T5(5)) ; T5(4) = 0.25*(3.0*T5(5) + T5(1)) ; T5(3) = 0.5*(T5(5) + T5(1))
        S5(2) = 0.25*(3.0*S5(1) + S5(5)) ; S5(4) = 0.25*(3.0*S5(5) + S5(1)) ; S5(3) = 0.5*(S5(5) + S5(1))
        p5(2) = 0.25*(3.0*p5(1) + p5(5)) ; p5(4) = 0.25*(3.0*p5(5) + p5(1)) ; p5(3) = 0.5*(p5(5) + p5(1))
        call calculate_density(T5, S5, p5, r5, tv%eqn_of_state, rho_ref=rho_ref)

        ! Note the consistency with the linear form below because (4.75 + 5.5/2) / 90 = 1/12
        inty_pa_cor_ri(i,J) = C1_90 * (4.75*(r5(5)-r5(1)) + 5.5*(r5(4)-r5(2))) * dgeo_y(i,J) - &
                              inty_pa_nonlin(i,J)
      enddo
    enddo

    ! Correct intx_pa and inty_pa at each interface using vertically constant corrections.
    do K=1,nz+1 ; do j=js,je ; do I=Isq,Ieq
      intx_pa(I,j,K) = intx_pa(I,j,K) + intx_pa_cor_ri(I,j)
    enddo ; enddo ; enddo

    do K=1,nz+1 ; do J=Jsq,Jeq ; do i=is,ie
      inty_pa(i,J,K) = inty_pa(i,J,K) + inty_pa_cor_ri(i,J)
    enddo ; enddo ; enddo
  endif ! intx_pa and inty_pa have now been reset to reflect the properties of an unimpeded interface.

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

  ! Calculate SAL geopotential anomaly and add its gradient to pressure gradient force
  if (CS%calculate_SAL .and. CS%tides_answer_date>20230630 .and. CS%bq_sal_tides) then
    !$OMP parallel do default(shared)
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        PFu(I,j,k) = PFu(I,j,k) + (e_sal(i+1,j) - e_sal(i,j)) * GV%g_Earth * G%IdxCu(I,j)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv(i,J,k) = PFv(i,J,k) + (e_sal(i,j+1) - e_sal(i,j)) * GV%g_Earth * G%IdyCv(i,J)
      enddo ; enddo
    enddo
  endif

  ! Calculate tidal geopotential anomaly and add its gradient to pressure gradient force
  if (CS%tides .and. CS%tides_answer_date>20230630 .and. CS%bq_sal_tides) then
    !$OMP parallel do default(shared)
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        PFu(I,j,k) = PFu(I,j,k) + ((e_tidal_eq(i+1,j) + e_tidal_sal(i+1,j)) &
          - (e_tidal_eq(i,j) + e_tidal_sal(i,j))) * GV%g_Earth * G%IdxCu(I,j)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv(i,J,k) = PFv(i,J,k) + ((e_tidal_eq(i,j+1) + e_tidal_sal(i,j+1)) &
          - (e_tidal_eq(i,j) + e_tidal_sal(i,j))) * GV%g_Earth * G%IdyCv(i,J)
      enddo ; enddo
    enddo
  endif

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
    call set_pbce_Bouss(e, tv_tmp, G, GV, US, rho0_set_pbce, CS%GFS_scale, pbce)
  endif

  if (present(eta)) then
    ! eta is the sea surface height relative to a time-invariant geoid, for comparison with
    ! what is used for eta in btstep.  See how e was calculated about 200 lines above.
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta(i,j) = e(i,j,1)*GV%Z_to_H
    enddo ; enddo
    if (CS%tides .and. (.not.CS%bq_sal_tides)) then
      if (CS%tides_answer_date>20230630) then
        !$OMP parallel do default(shared)
        do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
          eta(i,j) = eta(i,j) + (e_tidal_eq(i,j)+e_tidal_sal(i,j))*GV%Z_to_H
        enddo ; enddo
      else
        !$OMP parallel do default(shared)
        do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
          eta(i,j) = eta(i,j) + e_sal_and_tide(i,j)*GV%Z_to_H
        enddo ; enddo
      endif
    endif
    if (CS%calculate_SAL .and. (CS%tides_answer_date>20230630) .and. (.not.CS%bq_sal_tides)) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = eta(i,j) + e_sal(i,j)*GV%Z_to_H
      enddo ; enddo
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

  if (CS%id_MassWt_u>0) call post_data(CS%id_MassWt_u, MassWt_u, CS%diag)
  if (CS%id_MassWt_v>0) call post_data(CS%id_MassWt_v, MassWt_v, CS%diag)

  if (CS%id_rho_pgf>0) call post_data(CS%id_rho_pgf, rho_pgf, CS%diag)
  if (CS%id_rho_stanley_pgf>0) call post_data(CS%id_rho_stanley_pgf, rho_stanley_pgf, CS%diag)
  if (CS%id_p_stanley>0) call post_data(CS%id_p_stanley, p_stanley, CS%diag)

  ! Diagnostics for tidal forcing and SAL height anomaly
  if (CS%id_e_tide>0) then
    ! To be consistent with old runs, tidal forcing diagnostic also includes total SAL.
    ! New diagnostics are given for each individual field.
    if (CS%tides_answer_date>20230630) then ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e_sal_and_tide(i,j) = e_sal(i,j) + e_tidal_eq(i,j) + e_tidal_sal(i,j)
    enddo ; enddo ; endif
    call post_data(CS%id_e_tide, e_sal_and_tide, CS%diag)
  endif
  if (CS%id_e_sal>0) call post_data(CS%id_e_sal, e_sal, CS%diag)
  if (CS%id_e_tidal_eq>0) call post_data(CS%id_e_tidal_eq, e_tidal_eq, CS%diag)
  if (CS%id_e_tidal_sal>0) call post_data(CS%id_e_tidal_sal, e_tidal_sal, CS%diag)

  ! Diagnostics for tidal forcing and SAL horizontal gradients
  if (CS%calculate_SAL .and. ((associated(ADp%sal_u) .or. associated(ADp%sal_v)))) then
    if (CS%tides) then ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e_sal(i,j) = e_sal(i,j) + e_tidal_sal(i,j)
    enddo ; enddo ; endif
    if (CS%bq_sal_tides) then
      ! sal_u = ( e(i+1) - e(i) ) * g / dx
      if (associated(ADp%sal_u)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
        ADp%sal_u(I,j,k) = (e_sal(i+1,j) - e_sal(i,j)) * GV%g_Earth * G%IdxCu(I,j)
      enddo ; enddo ; enddo ; endif
      if (associated(ADp%sal_v)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
        ADp%sal_v(i,J,k) = (e_sal(i,j+1) - e_sal(i,j)) * GV%g_Earth * G%IdyCv(i,J)
      enddo ; enddo ; enddo ; endif
    else
      ! sal_u = ( e(i+1) - e(i) ) * g / dx * (rho(k) / rho0)
      if (associated(ADp%sal_u)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
        ADp%sal_u(I,j,k) = (e_sal(i+1,j) - e_sal(i,j)) * G%IdxCu(I,j) * I_Rho0 * &
          (2.0 * intx_dpa(I,j,k) * GV%Z_to_H / ((h(i,j,k) + h(i+1,j,k)) + h_neglect) + GxRho_ref)
      enddo ; enddo ; enddo ; endif
      if (associated(ADp%sal_v)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
        ADp%sal_v(i,J,k) = (e_sal(i,j+1) - e_sal(i,j)) * G%IdyCv(i,J) * I_Rho0 * &
          (2.0 * inty_dpa(i,J,k) * GV%Z_to_H / ((h(i,j,k) + h(i,j+1,k)) + h_neglect) + GxRho_ref)
      enddo ; enddo ; enddo ; endif
    endif
    if (CS%id_sal_u>0) call post_data(CS%id_sal_u, ADp%sal_u, CS%diag)
    if (CS%id_sal_v>0) call post_data(CS%id_sal_v, ADp%sal_v, CS%diag)
  endif

  if (CS%tides .and. ((associated(ADp%tides_u) .or. associated(ADp%tides_v)))) then
    if (CS%bq_sal_tides) then
      ! tides_u = ( e(i+1) - e(i) ) * g / dx
      if (associated(ADp%tides_u)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
        ADp%tides_u(I,j,k) = (e_tidal_eq(i+1,j) - e_tidal_eq(i,j)) * GV%g_Earth * G%IdxCu(I,j)
      enddo ; enddo ; enddo ; endif
      if (associated(ADp%tides_v)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
        ADp%tides_v(i,J,k) = (e_tidal_eq(i,j+1) - e_tidal_eq(i,j)) * GV%g_Earth * G%IdyCv(i,J)
      enddo ; enddo ; enddo ; endif
    else
      ! tides_u = ( e(i+1) - e(i) ) * g / dx * (rho(k) / rho0)
      if (associated(ADp%tides_u)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
        ADp%tides_u(I,j,k) = (e_tidal_eq(i+1,j) - e_tidal_eq(i,j)) * G%IdxCu(I,j) * I_Rho0 * &
          (2.0 * intx_dpa(I,j,k) * GV%Z_to_H / ((h(i,j,k) + h(i+1,j,k)) + h_neglect) + GxRho_ref)
      enddo ; enddo ; enddo ; endif
      if (associated(ADp%tides_v)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
        ADp%tides_v(i,J,k) = (e_tidal_eq(i,j+1) - e_tidal_eq(i,j)) * G%IdyCv(i,J) * I_Rho0 * &
          (2.0 * inty_dpa(i,J,k) * GV%Z_to_H / ((h(i,j,k) + h(i,j+1,k)) + h_neglect) + GxRho_ref)
      enddo ; enddo ; enddo ; endif
    endif
    if (CS%id_tides_u>0) call post_data(CS%id_tides_u, ADp%tides_u, CS%diag)
    if (CS%id_tides_v>0) call post_data(CS%id_tides_v, ADp%tides_v, CS%diag)
  endif
end subroutine PressureForce_FV_Bouss

!> Initializes the finite volume pressure gradient control structure
subroutine PressureForce_FV_init(Time, G, GV, US, param_file, diag, CS, ADp, SAL_CSp, tides_CSp)
  type(time_type), target,    intent(in)    :: Time !< Current model time
  type(ocean_grid_type),      intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),      intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),      intent(in)    :: param_file !< Parameter file handles
  type(diag_ctrl), target,    intent(inout) :: diag !< Diagnostics control structure
  type(PressureForce_FV_CS),  intent(inout) :: CS !< Finite volume PGF control structure
  type(accel_diag_ptrs),      pointer       :: ADp !< Acceleration diagnostic pointers
  type(SAL_CS),           intent(in), target, optional :: SAL_CSp !< SAL control structure
  type(tidal_forcing_CS), intent(in), target, optional :: tides_CSp !< Tides control structure

  ! Local variables
  real :: Stanley_coeff    ! Coefficient relating the temperature gradient and sub-gridscale
                           ! temperature variance [nondim]
  integer :: default_answer_date ! Global answer date
  logical :: use_temperature   ! If true, temperature and salinity are used as state variables.
  logical :: use_EOS           ! If true, density calculated from T & S using an equation of state.
  logical :: useMassWghtInterp ! If true, use near-bottom mass weighting for T and S
  logical :: MassWghtInterpTop ! If true, use near-surface mass weighting for T and S under ice shelves
  logical :: MassWghtInterp_NonBous_bug ! If true, use a buggy mass weighting when non-Boussinesq
  logical :: MassWghtInterpVanOnly ! If true, turn of mass weighting unless one side is vanished
  logical :: enable_bugs  ! If true, the defaults for recently added bug-fix flags are set to
                          ! recreate the bugs, or if false bugs are only used if actively selected.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl  ! This module's name.
  logical :: use_ALE       ! If true, use the Vertical Lagrangian Remap algorithm
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz

  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%initialized = .true.
  CS%diag => diag ; CS%Time => Time
  if (present(tides_CSp)) &
    CS%tides_CSp => tides_CSp
  if (present(SAL_CSp)) &
    CS%SAL_CSp => SAL_CSp

  mdl = "MOM_PressureForce_FV"
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true., do_not_log=.true.)
  call get_param(param_file, mdl, "RHO_PGF_REF", CS%rho_ref, &
                 "The reference density that is subtracted off when calculating pressure "//&
                 "gradient forces.  Its inverse is subtracted off of specific volumes when "//&
                 "in non-Boussinesq mode.  The default is RHO_0.", &
                 units="kg m-3", default=GV%Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
  call get_param(param_file, mdl, "RHO_PGF_REF_BUG", CS%rho_ref_bug, &
                 "If true, recover a bug that RHO_0 (the mean seawater density in Boussinesq mode) "//&
                 "and RHO_PGF_REF (the subtracted reference density in finite volume pressure "//&
                 "gradient forces) are incorrectly interchanged in several instances in Boussinesq mode.", &
                 default=enable_bugs)
  call get_param(param_file, mdl, "TIDES", CS%tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  if (CS%tides) then
    call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
    call get_param(param_file, mdl, "TIDES_ANSWER_DATE", CS%tides_answer_date, "The vintage of "//&
                  "self-attraction and loading (SAL) and tidal forcing calculations.  Setting "//&
                  "dates before 20230701 recovers old answers (Boussinesq and non-Boussinesq "//&
                  "modes) when SAL is part of the tidal forcing calculation.  The answer "//&
                  "difference is only at bit level and due to a reordered summation.  Setting "//&
                  "dates before 20250201 recovers answers (Boussinesq mode) that interface "//&
                  "heights are modified before pressure force integrals are calculated.", &
                  default=default_answer_date, do_not_log=(.not.CS%tides))
  endif
  call get_param(param_file, mdl, "CALCULATE_SAL", CS%calculate_SAL, &
                 "If true, calculate self-attraction and loading.", default=CS%tides)
  if (CS%calculate_SAL) &
    call get_param(param_file, '', "SAL_USE_BPA", CS%sal_use_bpa, default=.false., &
                   do_not_log=.true.)
  if ((CS%tides .or. CS%calculate_SAL) .and. GV%Boussinesq) &
    call get_param(param_file, mdl, "BOUSSINESQ_SAL_TIDES", CS%bq_sal_tides, "If true, "//&
                   "in Boussinesq mode, use an alternative method to include self-attraction "//&
                   "and loading (SAL) and tidal forcings in pressure gradient, in which their "//&
                   "gradients are calculated separately, instead of adding geopotential "//&
                   "anomalies as corrections to the interface height.  This alternative method "//&
                   "elimates a baroclinic component of the SAL and tidal forcings.", &
                   default=.false.)
  call get_param(param_file, "MOM", "ENABLE_THERMODYNAMICS", use_temperature, &
                 "If true, Temperature and salinity are used as state variables.", &
                 default=.true., do_not_log=.true.)
  call get_param(param_file, "MOM", "USE_EOS", use_EOS, &
                 "If true,  density is calculated from temperature and "//&
                 "salinity with an equation of state.  If USE_EOS is "//&
                 "true, ENABLE_THERMODYNAMICS must be true as well.", &
                 default=use_temperature, do_not_log=.true.)

  call get_param(param_file, mdl, "SSH_IN_EOS_PRESSURE_FOR_PGF", CS%use_SSH_in_Z0p, &
                 "If true, include contributions from the sea surface height in the height-based "//&
                 "pressure used in the equation of state calculations for the Boussinesq pressure "//&
                 "gradient forces, including adjustments for atmospheric or sea-ice pressure.", &
                 default=.false., do_not_log=.not.GV%Boussinesq)
  if (CS%tides .and. CS%tides_answer_date<=20250131 .and. CS%use_SSH_in_Z0p) &
    call MOM_error(FATAL, trim(mdl) // ", PressureForce_FV_init: SSH_IN_EOS_PRESSURE_FOR_PGF "//&
                   "needs to be FALSE to recover tide answers before 20250131.")

  call get_param(param_file, "MOM", "USE_REGRIDDING", use_ALE, &
                 "If True, use the ALE algorithm (regridding/remapping). "//&
                 "If False, use the layered isopycnal algorithm.", default=.false. )
  call get_param(param_file, mdl, "MASS_WEIGHT_IN_PRESSURE_GRADIENT", useMassWghtInterp, &
                 "If true, use mass weighting when interpolating T/S for integrals "//&
                 "near the bathymetry in FV pressure gradient calculations.", &
                 default=.false.)
  call get_param(param_file, mdl, "MASS_WEIGHT_IN_PRESSURE_GRADIENT_TOP", MassWghtInterpTop, &
                 "If true and MASS_WEIGHT_IN_PRESSURE_GRADIENT is true, use mass weighting when "//&
                 "interpolating T/S for integrals near the top of the water column in FV "//&
                 "pressure gradient calculations. ", &
                 default=.false.) !### Change Default to MASS_WEIGHT_IN_PRESSURE_GRADIENT?
  call get_param(param_file, mdl, "MASS_WEIGHT_IN_PGF_NONBOUS_BUG", MassWghtInterp_NonBous_bug, &
                 "If true, use a masking bug in non-Boussinesq calculations with mass weighting "//&
                 "when interpolating T/S for integrals near the bathymetry in FV pressure "//&
                 "gradient calculations.", &
                 default=.false., do_not_log=(GV%Boussinesq .or. (.not.useMassWghtInterp)))
  call get_param(param_file, mdl, "MASS_WEIGHT_IN_PGF_VANISHED_ONLY", CS%MassWghtInterpVanOnly, &
                 "If true, use mass weighting when interpolating T/S for integrals "//&
                 "only if one side is vanished according to RESET_INTXPA_H_NONVANISHED. ", &
                 default=.false.)

  CS%MassWghtInterp = 0
  if (useMassWghtInterp) &
    CS%MassWghtInterp = ibset(CS%MassWghtInterp, 0) ! Same as CS%MassWghtInterp + 1
  if (MassWghtInterpTop) &
    CS%MassWghtInterp = ibset(CS%MassWghtInterp, 1) ! Same as CS%MassWghtInterp + 2
  if ((.not.GV%Boussinesq) .and. MassWghtInterp_NonBous_bug) &
    CS%MassWghtInterp = ibset(CS%MassWghtInterp, 3) ! Same as CS%MassWghtInterp + 8

  call get_param(param_file, mdl, "CORRECTION_INTXPA", CS%correction_intxpa, &
                 "If true, use a correction for surface pressure curvature in intx_pa.", &
                 default=.false., do_not_log=.not.use_EOS)
  call get_param(param_file, mdl, "RESET_INTXPA_INTEGRAL", CS%reset_intxpa_integral, &
                 "If true, reset INTXPA to match pressures at first nonvanished cell. "//&
                 "Includes pressure correction.", default=.false., do_not_log=.not.use_EOS)
  call get_param(param_file, mdl, "RESET_INTXPA_INTEGRAL_FLATTEST", CS%reset_intxpa_flattest, &
                 "If true, use flattest interface as reference interface where there is no "//&
                 "better choice for RESET_INTXPA_INTEGRAL. Otherwise, use surface interface.", &
                 default=.false., do_not_log=.not.use_EOS)
  if (.not.use_EOS) then  ! These options do nothing without an equation of state.
    CS%correction_intxpa = .false.
    CS%reset_intxpa_integral = .false.
    CS%reset_intxpa_flattest = .false.
  endif
  call get_param(param_file, mdl, "RESET_INTXPA_H_NONVANISHED", CS%h_nonvanished, &
                 "A minimal layer thickness that indicates that a layer is thick enough to usefully "//&
                 "reestimate the pressure integral across the interface below.", &
                 default=1.0e-6, units="m", scale=GV%m_to_H, do_not_log=.not.CS%reset_intxpa_integral)
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
    CS%id_e_sal = register_diag_field('ocean_model', 'e_sal', diag%axesT1, Time, &
        'Self-attraction and loading height anomaly', 'meter', conversion=US%Z_to_m)
    CS%id_sal_u = register_diag_field('ocean_model', 'SAL_u', diag%axesCuL, Time, &
        'Zonal Acceleration due to self-attraction and loading', 'm s-2', conversion=US%L_T2_to_m_s2)
    CS%id_sal_v = register_diag_field('ocean_model', 'SAL_v', diag%axesCvL, Time, &
        'Meridional Acceleration due to self-attraction and loading', 'm s-2', conversion=US%L_T2_to_m_s2)
    if (CS%id_sal_u > 0) &
      call safe_alloc_ptr(ADp%sal_u, IsdB, IedB, jsd, jed, nz)
    if (CS%id_sal_v > 0) &
      call safe_alloc_ptr(ADp%sal_v, isd, ied, JsdB, JedB, nz)
  endif
  if (CS%tides) then
    CS%id_e_tide = register_diag_field('ocean_model', 'e_tidal', diag%axesT1, Time, &
        'Tidal Forcing Astronomical and SAL Height Anomaly', 'meter', conversion=US%Z_to_m)
    CS%id_e_tidal_eq  = register_diag_field('ocean_model', 'e_tide_eq', diag%axesT1, Time, &
        'Equilibrium tides height anomaly', 'meter', conversion=US%Z_to_m)
    CS%id_e_tidal_sal = register_diag_field('ocean_model', 'e_tide_sal', diag%axesT1, Time, &
        'Read-in tidal self-attraction and loading height anomaly', 'meter', conversion=US%Z_to_m)
    CS%id_tides_u = register_diag_field('ocean_model', 'tides_u', diag%axesCuL, Time, &
        'Zonal Acceleration due to tidal forcing', 'm s-2', conversion=US%L_T2_to_m_s2)
    CS%id_tides_v = register_diag_field('ocean_model', 'tides_v', diag%axesCvL, Time, &
        'Meridional Acceleration due to tidal forcing', 'm s-2', conversion=US%L_T2_to_m_s2)
    if (CS%id_tides_u > 0) &
      call safe_alloc_ptr(ADp%tides_u, IsdB, IedB, jsd, jed, nz)
    if (CS%id_tides_v > 0) &
      call safe_alloc_ptr(ADp%tides_v, isd, ied, JsdB, JedB, nz)
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
