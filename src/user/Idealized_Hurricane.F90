!> Forcing for the idealized hurricane and SCM_idealized_hurricane examples.
module Idealized_hurricane

! This file is part of MOM6. See LICENSE.md for the license.

! History
!--------
! November 2014: Origination.
! October 2018: Renamed module from SCM_idealized_hurricane to idealized_hurricane
!               This module is no longer exclusively for use in SCM mode.
!               Legacy code that can be deleted is at the bottom (currently maintained
!               only to preserve exact answers in SCM mode).
!               The T/S initializations have been removed since they are redundant
!               w/ T/S initializations in CVMix_tests (which should be moved
!               into the main state_initialization to their utility
!               for multiple example cases).
! December 2024: Removed the legacy subroutine SCM_idealized_hurricane_wind_forcing

use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, mech_forcing
use MOM_forcing_type, only : allocate_mech_forcing
use MOM_grid, only : ocean_grid_type
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/), time_type_to_real
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs, surface
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public idealized_hurricane_wind_init !Public interface to initialize the idealized
                                     ! hurricane wind profile.
public idealized_hurricane_wind_forcing !Public interface to update the idealized
                                        ! hurricane wind profile.

!> Container for parameters describing idealized wind structure
type, public :: idealized_hurricane_CS ; private

  ! Parameters used to compute Holland radial wind profile
  real    :: rho_a                !< Mean air density [R ~> kg m-3]
  real    :: pressure_ambient     !< Pressure at surface of ambient air [R L2 T-2 ~> Pa]
  real    :: pressure_central     !< Pressure at surface at hurricane center [R L2 T-2 ~> Pa]
  real    :: rad_max_wind         !< Radius of maximum winds [L ~> m]
  real    :: rad_edge             !< Radius of the edge of the hurricane, normalized by
                                  !! the radius of maximum winds [nondim]
  real    :: rad_ambient          !< Radius at which the winds are at their ambient background values,
                                  !! normalized by the radius of maximum winds [nondim]
  real    :: max_windspeed        !< Maximum wind speeds [L T-1 ~> m s-1]
  real    :: hurr_translation_spd !< Hurricane translation speed [L T-1 ~> m s-1]
  real    :: hurr_translation_dir !< Hurricane translation direction [radians]
  real    :: gustiness            !< Gustiness (used in u*) [R Z2 T-2 ~> Pa]
  real    :: Rho0                 !< A reference ocean density [R ~> kg m-3]
  real    :: Hurr_cen_Y0          !< The initial y position of the hurricane
                                  !!  This experiment is conducted in a Cartesian
                                  !!  grid and this is assumed to be in meters [L ~> m]
  real    :: Hurr_cen_X0          !< The initial x position of the hurricane
                                  !!  This experiment is conducted in a Cartesian
                                  !!  grid and this is assumed to be in meters [L ~> m]
  real    :: Holland_B            !< Parameter 'B' from the Holland formula [nondim]
  logical :: relative_tau         !< A logical to take difference between wind
                                  !! and surface currents to compute the stress
  integer :: answer_date          !< The vintage of the expressions in the idealized hurricane
                                  !! test case.  Values below 20190101 recover the answers
                                  !! from the end of 2018, while higher values use expressions
                                  !! that are rescalable and respect rotational symmetry.
  ! Parameters used in a simple wind-speed dependent expression for C_drag
  real :: Cd_calm       !< The drag coefficient with weak relative winds [nondim]
  real :: calm_speed    !< The relative wind speed below which the drag coefficient takes its
                        !! calm value [L T-1 ~> m s-1]
  real :: Cd_windy      !< The drag coefficient with strong relative winds [nondim]
  real :: windy_speed   !< The relative wind speed below which the drag coefficient takes its
                        !! windy value [L T-1 ~> m s-1]
  real :: dCd_dU10      !< The partial derivative of the drag coefficient times 1000 with the 10 m
                        !! wind speed for intermediate wind speeds [T L-1 ~> s m-1]
  real :: Cd_intercept  !< The zero-wind intercept times 1000 of the linear fit for the drag
                        !! coefficient for the intermediate speeds where there is a linear
                        !! dependence on the 10 m wind speed [nondim]

  ! Parameters used to set the inflow angle as a function of radius and maximum wind speed
  real :: A0_0          !< The zero-radius, zero-speed intercept of the axisymmetric inflow angle [degrees]
  real :: A0_Rnorm      !< The normalized radius dependence of the axisymmetric inflow angle [degrees]
  real :: A0_speed      !< The maximum wind speed dependence of the axisymmetric inflow angle
                        !! [degrees T L-1 ~> degrees s m-1]
  real :: A1_0          !< The zero-radius, zero-speed intercept of the normalized inflow angle
                        !! asymmetry [degrees]
  real :: A1_Rnorm      !< The normalized radius dependence of the normalized inflow angle asymmetry [degrees]
  real :: A1_speed      !< The translation speed dependence of the normalized inflow angle asymmetry
                        !! [degrees T L-1 ~> degrees s m-1]
  real :: P1_0          !< The zero-radius, zero-speed intercept of the angle difference between the
                        !! translation direction and the inflow direction [degrees]
  real :: P1_Rnorm      !< The normalized radius dependence of the angle difference between the
                        !! translation direction and the inflow direction [degrees]
  real :: P1_speed      !< The translation speed dependence of the angle difference between the
                        !! translation direction and the inflow direction [degrees T L-1 ~> degrees s m-1]

  ! Parameters used if in SCM (single column model) mode
  logical :: SCM_mode   !< If true this being used in Single Column Model mode
  logical :: edge_taper_bug !< If true and SCM_mode is true, use a bug that does all of the tapering
                        !! and inflow angle calculations for radii between RAD_EDGE and RAD_AMBIENT
                        !! as though they were at RAD_EDGE.
  real :: f_column      !< Coriolis parameter used in the single column mode idealized
                        !! hurricane wind profile [T-1 ~> s-1]
  logical :: BR_Bench   !< A "benchmark" configuration (which is meant to
                        !! provide identical wind to reproduce a previous
                        !! experiment, where that wind formula contained an error)
  real    :: dy_from_center  !< (Fixed) distance in y from storm center path [L ~> m]

  real :: pi      !< The circumference of a circle divided by its diameter [nondim]
  real :: Deg2Rad !< The conversion factor from degrees to radians [radian degree-1]

end type

character(len=40)  :: mdl = "idealized_hurricane" !< This module's name.

contains

!> Initializes wind profile for the SCM idealized hurricane example
subroutine idealized_hurricane_wind_init(Time, G, US, param_file, CS)
  type(time_type),               intent(in) :: Time   !< Model time
  type(ocean_grid_type),         intent(in) :: G      !< Grid structure
  type(unit_scale_type),         intent(in) :: US     !< A dimensional unit scaling type
  type(param_file_type),         intent(in) :: param_file !< Input parameter structure
  type(idealized_hurricane_CS),  pointer    :: CS     !< Parameter container for this module

  ! Local variables
  real :: dP  ! The pressure difference across the hurricane [R L2 T-2 ~> Pa]
  real :: C   ! A temporary variable in units of the square root of a specific volume [sqrt(m3 kg-1)]
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  logical :: continuous_Cd  ! If true, use a continuous form for the simple drag coefficient as a
                 ! function of wind speed with the idealized hurricane.  When this is false, the
                 ! linear shape for the mid-range wind speeds is specified separately.

  ! This include declares and sets the variable "version".
# include "version_variable.h"

  if (associated(CS)) then
    call MOM_error(FATAL, "idealized_hurricane_wind_init called "// &
                          "with an associated control structure.")
    return
  endif

  allocate(CS)

  CS%pi = 4.0*atan(1.0)
  CS%Deg2Rad = CS%pi/180.

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  ! Parameters for computing a wind profile
  call get_param(param_file, mdl, "IDL_HURR_RHO_AIR", CS%rho_a, &
                 "Air density used to compute the idealized hurricane wind profile.", &
                 units='kg/m3', default=1.2, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "IDL_HURR_AMBIENT_PRESSURE", CS%pressure_ambient, &
                 "Ambient pressure used in the idealized hurricane wind profile.", &
                 units='Pa', default=101200., scale=US%Pa_to_RL2_T2)
  call get_param(param_file, mdl, "IDL_HURR_CENTRAL_PRESSURE", CS%pressure_central, &
                 "Central pressure used in the idealized hurricane wind profile.", &
                 units='Pa', default=96800., scale=US%Pa_to_RL2_T2)
  call get_param(param_file, mdl, "IDL_HURR_RAD_MAX_WIND", CS%rad_max_wind, &
                 "Radius of maximum winds used in the idealized hurricane wind profile.", &
                 units='m', default=50.e3, scale=US%m_to_L)
  call get_param(param_file, mdl, "IDL_HURR_RAD_EDGE", CS%rad_edge, &
                 "Radius of the edge of the hurricane, normalized by the radius of maximum winds.", &
                 units='nondim', default=10.0)
  call get_param(param_file, mdl, "IDL_HURR_RAD_AMBIENT", CS%rad_ambient, &
                 "Radius at which the winds are at their ambient background values, "//&
                 "normalized by the radius of maximum winds.", &
                 units='nondim', default=CS%rad_edge+2.0)
  call get_param(param_file, mdl, "IDL_HURR_MAX_WIND", CS%max_windspeed, &
                 "Maximum wind speed used in the idealized hurricane wind profile.", &
                 units='m/s', default=65., scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "IDL_HURR_TRAN_SPEED", CS%hurr_translation_spd, &
                 "Translation speed of hurricane used in the idealized hurricane wind profile.", &
                 units='m/s', default=5.0, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "IDL_HURR_TRAN_DIR", CS%hurr_translation_dir, &
                 "Translation direction (towards) of hurricane used in the "//&
                 "idealized hurricane wind profile.", &
                 units='degrees', default=180.0, scale=CS%Deg2Rad)
  call get_param(param_file, mdl, "IDL_HURR_X0", CS%Hurr_cen_X0, &
                 "Idealized Hurricane initial X position", &
                 units='m', default=0., scale=US%m_to_L)
  call get_param(param_file, mdl, "IDL_HURR_Y0", CS%Hurr_cen_Y0, &
                 "Idealized Hurricane initial Y position", &
                 units='m', default=0., scale=US%m_to_L)
  call get_param(param_file, mdl, "IDL_HURR_TAU_CURR_REL", CS%relative_tau, &
                 "Current relative stress switch used in the idealized hurricane wind profile.", &
                 default=.false.)

  call get_param(param_file, mdl, "IDL_HURR_AXI_INFLOW_0", CS%A0_0, &
                 "The zero-radius asymmetry, zero-speed intercept of the axisymmetric inflow "//&
                 "angle for the parametric idealized hurricane.", &
                 default=-14.33, units="degrees")
  call get_param(param_file, mdl, "IDL_HURR_AXI_INFLOW_RNORM", CS%A0_Rnorm, &
                 "The normalized radius dependence of the axisymmetric inflow angle "//&
                 "for the parametric idealized hurricane.", &
                 default=-0.9, units="degrees")
  call get_param(param_file, mdl, "IDL_HURR_AXI_INFLOW_MAX_SPEED", CS%A0_speed, &
                 "The maximum wind speed dependence of the axisymmetric inflow angle "//&
                 "for the parametric idealized hurricane.", &
                 default=-0.09, units="degrees s m-1", scale=US%L_T_to_m_s)
  call get_param(param_file, mdl, "IDL_HURR_ASYM_INFLOW_0", CS%A1_0, &
                 "The zero-radius, zero-speed intercept of the normalized inflow angle asymmetry "//&
                 "for the parametric idealized hurricane.", &
                 default=0.14, units="degrees")
  call get_param(param_file, mdl, "IDL_HURR_ASYM_INFLOW_RNORM", CS%A1_Rnorm, &
                 "The normalized radius dependence of the normalized inflow angle asymmetry "//&
                 "for the parametric idealized hurricane.", &
                 default=0.04, units="degrees")
  call get_param(param_file, mdl, "IDL_HURR_ASYM_INFLOW_TR_SPEED", CS%A1_speed, &
                 "The translation speed dependence of the normalized inflow angle asymmetry "//&
                 "for the parametric idealized hurricane.", &
                 default=0.05, units="degrees s m-1", scale=US%L_T_to_m_s)
  call get_param(param_file, mdl, "IDL_HURR_INFLOW_DANGLE_0", CS%P1_0, &
                 "The zero-radius, zero-speed intercept of the angle difference between the "//&
                 "translation direction and the inflow direction "//&
                 "for the parametric idealized hurricane.", &
                 default=85.31, units="degrees")
  call get_param(param_file, mdl, "IDL_HURR_INFLOW_DANGLE_RNORM", CS%P1_Rnorm, &
                 "The normalized radius dependence of the angle difference between the "//&
                 "translation direction and the inflow direction "//&
                 "for the parametric idealized hurricane.", &
                 default=6.88, units="degrees")
  call get_param(param_file, mdl, "IDL_HURR_INFLOW_DANGLE_TR_SPEED", CS%P1_speed, &
                 "The translation speed dependence of the angle difference between the "//&
                 "translation direction and the inflow direction"//&
                 "for the parametric idealized hurricane.", &
                 default=-9.60, units="degrees s m-1", scale=US%L_T_to_m_s)

  ! Parameters for SCM mode
  call get_param(param_file, mdl, "IDL_HURR_SCM_BR_BENCH", CS%BR_Bench, &
                 "Single column mode benchmark case switch, which is "// &
                 "invoking a modification (bug) in the wind profile meant to "//&
                 "reproduce a previous implementation.", default=.false.)
  call get_param(param_file, mdl, "IDL_HURR_SCM", CS%SCM_mode, &
                 "Single Column mode switch used in the SCM idealized hurricane wind profile.", &
                 default=.false.)
  call get_param(param_file, mdl, "IDL_HURR_SCM_EDGE_TAPER_BUG", CS%edge_taper_bug, &
                 "If true and IDL_HURR_SCM is true, use a bug that does all of the tapering and "//&
                 "inflow angle calculations for radii between RAD_EDGE and RAD_AMBIENT as though "//&
                 "they were at RAD_EDGE.", &
                 default=.false., do_not_log=.not.CS%SCM_mode)
  if (.not.CS%SCM_mode) CS%edge_taper_bug = .false.
  call get_param(param_file, mdl, "IDL_HURR_SCM_LOCY", CS%dy_from_center, &
                 "Y distance of station used in the SCM idealized hurricane wind profile.", &
                 units='m', default=50.e3, scale=US%m_to_L)
  call get_param(param_file, mdl, "IDL_HURR_SCM_CORIOLIS", CS%f_column, &
                 "Coriolis parameter used in the single column mode idealized hurricane wind profile.", &
                 units='s-1', default=5.5659e-05, scale=US%T_to_s, do_not_log=.not.CS%BR_Bench) ! (CS%SCM_mode)

  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "IDL_HURR_ANSWER_DATE", CS%answer_date, &
                 "The vintage of the expressions in the idealized hurricane test case.  "//&
                 "Values below 20190101 recover the answers from the end of 2018, while higher "//&
                 "values use expressions that are rescalable and respect rotational symmetry.", &
                 default=default_answer_date)

  ! Parameters for the simple Cdrag expression
  call get_param(param_file, mdl, "IDL_HURR_CD_CALM", CS%Cd_calm, &
                 "The drag coefficient with weak relative winds "//&
                 "for the simple drag coefficient expression used with the idealized hurricane.", &
                 units='nondim', default=1.2e-3)
  call get_param(param_file, mdl, "IDL_HURR_CD_CALM_SPEED", CS%calm_speed, &
                 "The relative wind speed below which the drag coefficient takes its calm value "//&
                 "for the simple drag coefficient expression used with the idealized hurricane.", &
                 units='m s-1', default=11.0, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "IDL_HURR_CD_WINDY", CS%Cd_windy, &
                 "The drag coefficient with strong relative winds "//&
                 "for the simple drag coefficient expression used with the idealized hurricane.", &
                 units='nondim', default=1.8e-3)
  call get_param(param_file, mdl, "IDL_HURR_CD_WINDY_SPEED", CS%windy_speed, &
                 "The relative wind speed below which the drag coefficient takes its windy value "//&
                 "for the simple drag coefficient expression used with the idealized hurricane.", &
                 units='m s-1', default=20.0, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "IDL_HURR_CD_CONTINUOUS", continuous_Cd, &
                 "If true, use a continuous form for the simple drag coefficient as a function of "//&
                 "wind speed with the idealized hurricane.  When this is false, the linear shape "//&
                 "for the mid-range wind speeds is specified separately.", &
                 default=.false.)
  call get_param(param_file, mdl, "IDL_HURR_CD_DCD_DU10", CS%dCd_dU10, &
                 "The partial derivative of the drag coefficient times 1000 with the 10 m wind speed "//&
                 "for the simple drag coefficient expression used with the idealized hurricane.", &
                 units="s m-1", default=0.065, scale=US%L_T_to_m_s, do_not_log=continuous_Cd)
  call get_param(param_file, mdl, "IDL_HURR_CD_INTERCEPT", CS%Cd_intercept, &
                 "The zero-wind intercept times 1000 of the linear fit for the drag coefficient "//&
                 "for the intermediate speeds where there is a linear dependence on the 10 m wind speed "//&
                 "for the simple drag coefficient expression used with the idealized hurricane.", &
                 units="nondim", default=0.49, do_not_log=continuous_Cd)
  if (continuous_Cd) then
    if (CS%windy_speed > CS%calm_speed) then
      CS%dCd_dU10 = (CS%Cd_windy - CS%Cd_calm) / (CS%windy_speed - CS%calm_speed)
      CS%Cd_intercept = CS%Cd_calm - CS%dCd_dU10 * CS%calm_speed
    else
      CS%dCd_dU10 = 0.0
      CS%Cd_intercept = CS%Cd_windy
    endif
  endif


  ! The following parameters are model run-time parameters which are used
  ! and logged elsewhere and so should not be logged here. The default
  ! value should be consistent with the rest of the model.
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINSEQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R, do_not_log=.true.)
  call get_param(param_file, mdl, "GUST_CONST", CS%gustiness, &
                 "The background gustiness in the winds.", &
                 units="Pa", default=0.0, scale=US%Pa_to_RLZ_T2*US%L_to_Z, do_not_log=.true.)

  if (CS%rad_edge >= CS%rad_ambient) call MOM_error(FATAL, &
    "idealized_hurricane_wind_init: IDL_HURR_RAD_AMBIENT must be larger than IDL_HURR_RAD_EDGE.")

  dP = CS%pressure_ambient - CS%pressure_central
  if (CS%answer_date < 20190101) then
    C = CS%max_windspeed / sqrt( US%R_to_kg_m3 * dP )
    CS%Holland_B = C**2 * US%R_to_kg_m3*CS%rho_a * exp(1.0)
  else
    CS%Holland_B = CS%max_windspeed**2 * CS%rho_a * exp(1.0) / dP
  endif

end subroutine idealized_hurricane_wind_init

!> Computes the surface wind for the idealized hurricane test cases
subroutine idealized_hurricane_wind_forcing(sfc_state, forces, day, G, US, CS)
  type(surface),                intent(in)    :: sfc_state  !< Surface state structure
  type(mech_forcing),           intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),              intent(in)    :: day    !< Time in days
  type(ocean_grid_type),        intent(inout) :: G      !< Grid structure
  type(unit_scale_type),        intent(in)    :: US     !< A dimensional unit scaling type
  type(idealized_hurricane_CS), pointer       :: CS     !< Container for idealized hurricane parameters

  ! Local variables
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  real :: TX, TY      !< wind stress components [R L Z T-2 ~> Pa]
  real :: Uocn, Vocn  !< Surface ocean velocity components [L T-1 ~> m s-1]
  real :: YY, XX      !< storm relative position [L ~> m]
  real :: XC, YC      !< Storm center location [L ~> m]
  real :: f_local     !< Local Coriolis parameter [T-1 ~> s-1]
  real :: fbench      !< The benchmark 'f' value [T-1 ~> s-1]
  real :: fbench_fac  !< A factor that is set to 0 to use the
                      !!  benchmark 'f' value [nondim]
  real :: rel_tau_fac !< A factor that is set to 0 to disable
                      !!  current relative stress calculation [nondim]

  ! Bounds for loops and memory allocation
  is = G%isc    ; ie = G%iec    ; js = G%jsc    ; je = G%jec
  Isq = G%IscB  ; Ieq = G%IecB  ; Jsq = G%JscB  ; Jeq = G%JecB
  isd = G%isd   ; ied = G%ied   ; jsd = G%jsd   ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if ((G%grid_unit_to_L <= 0.) .and. (.not.CS%SCM_mode)) call MOM_error(FATAL, "Idealized_Hurricane.F90: " //&
          "idealized_hurricane_wind_forcing is only set to work with Cartesian axis units.")

  ! Allocate the forcing arrays, if necessary.
  call allocate_mech_forcing(G, forces, stress=.true., ustar=.true., tau_mag=.true.)

  if (CS%relative_tau) then
    REL_TAU_FAC = 1.
  else
    REL_TAU_FAC = 0. !Multiplied to 0 surface current
  endif

  !> Compute storm center location
  XC = CS%Hurr_cen_X0 + (time_type_to_real(day)*US%s_to_T * CS%hurr_translation_spd * &
       cos(CS%hurr_translation_dir))
  YC = CS%Hurr_cen_Y0 + (time_type_to_real(day)*US%s_to_T * CS%hurr_translation_spd * &
       sin(CS%hurr_translation_dir))

  if (CS%BR_Bench) then
    ! f reset to value used in generated wind for benchmark test
    fbench = CS%f_column
    fbench_fac = 0.0
  else
    fbench = 0.0
    fbench_fac = 1.0
  endif

  !> Computes taux
  do j=js,je
    do I=is-1,Ieq
      Uocn = sfc_state%u(I,j) * REL_TAU_FAC
      if (CS%answer_date < 20190101) then
        Vocn = 0.25*(sfc_state%v(i,J)+sfc_state%v(i+1,J-1)&
                    +sfc_state%v(i+1,J)+sfc_state%v(i,J-1))*REL_TAU_FAC
      else
        Vocn = 0.25*((sfc_state%v(i,J)+sfc_state%v(i+1,J-1)) +&
                     (sfc_state%v(i+1,J)+sfc_state%v(i,J-1))) * REL_TAU_FAC
      endif
      f_local = abs(0.5*(G%CoriolisBu(I,J)+G%CoriolisBu(I,J-1)))*fbench_fac + fbench
      ! Calculate position as a function of time.
      if (CS%SCM_mode) then
        YY = CS%dy_from_center - YC
        XX = -XC
      else
        YY = G%geoLatCu(I,j) * G%grid_unit_to_L - YC
        XX = G%geoLonCu(I,j) * G%grid_unit_to_L - XC
      endif
      call idealized_hurricane_wind_profile(CS, US, f_local, YY, XX, Uocn, Vocn, TX, TY)
      forces%taux(I,j) = G%mask2dCu(I,j) * TX
    enddo
  enddo
  !> Computes tauy
  do J=js-1,Jeq
    do i=is,ie
      if (CS%answer_date < 20190101) then
        Uocn = 0.25*(sfc_state%u(I,j)+sfc_state%u(I-1,j+1) + &
                     sfc_state%u(I-1,j)+sfc_state%u(I,j+1))*REL_TAU_FAC
      else
        Uocn = 0.25*((sfc_state%u(I,j)+sfc_state%u(I-1,j+1)) + &
                     (sfc_state%u(I-1,j)+sfc_state%u(I,j+1))) * REL_TAU_FAC
      endif
      Vocn = sfc_state%v(i,J) * REL_TAU_FAC
      f_local = abs(0.5*(G%CoriolisBu(I-1,J)+G%CoriolisBu(I,J)))*fbench_fac + fbench
      ! Calculate position as a function of time.
      if (CS%SCM_mode) then
        YY = CS%dy_from_center - YC
        XX = -XC
      else
        YY = G%geoLatCv(i,J) * G%grid_unit_to_L - YC
        XX = G%geoLonCv(i,J) * G%grid_unit_to_L - XC
      endif
      call idealized_hurricane_wind_profile(CS, US, f_local, YY, XX, Uocn, Vocn, TX, TY)
      forces%tauy(i,J) = G%mask2dCv(i,J) * TY
    enddo
  enddo

  !> Get Ustar
  if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
    !  This expression can be changed if desired, but need not be.
    forces%ustar(i,j) = G%mask2dT(i,j) * sqrt(CS%gustiness/CS%Rho0 + &
            US%L_to_Z * sqrt(0.5*((forces%taux(I-1,j)**2) + (forces%taux(I,j)**2)) + &
                             0.5*((forces%tauy(i,J-1)**2) + (forces%tauy(i,J)**2)))/CS%Rho0)
  enddo ; enddo ; endif

  !> Get tau_mag [R Z2 T-2 ~> Pa]
  if (associated(forces%tau_mag)) then ; do j=js,je ; do i=is,ie
    forces%tau_mag(i,j) = G%mask2dT(i,j) * (CS%gustiness + &
            US%L_to_Z * sqrt(0.5*((forces%taux(I-1,j)**2) + (forces%taux(I,j)**2)) + &
                             0.5*((forces%tauy(i,J-1)**2) + (forces%tauy(i,J)**2))))
  enddo ; enddo ; endif

end subroutine idealized_hurricane_wind_forcing

!> Calculate the wind speed at a location as a function of time.
subroutine idealized_hurricane_wind_profile(CS, US, absf, YY, XX, UOCN, VOCN, Tx, Ty)
  ! Author: Brandon Reichl
  ! Date: Nov-20-2014
  !       Aug-14-2018 Generalized for non-SCM configuration

  ! Input parameters
  type(idealized_hurricane_CS), pointer       :: CS   !< Container for idealized hurricane parameters
  type(unit_scale_type),        intent(in)    :: US     !< A dimensional unit scaling type
  real, intent(in)  :: absf !< Input Coriolis magnitude [T-1 ~> s-1]
  real, intent(in)  :: YY   !< Location in m relative to center y [L ~> m]
  real, intent(in)  :: XX   !< Location in m relative to center x [L ~> m]
  real, intent(in)  :: UOCN !< X surface current [L T-1 ~> m s-1]
  real, intent(in)  :: VOCN !< Y surface current [L T-1 ~> m s-1]
  real, intent(out) :: Tx   !< X stress [R L Z T-2 ~> Pa]
  real, intent(out) :: Ty   !< Y stress [R L Z T-2 ~> Pa]

  ! Local variables

  ! Wind profile terms
  real :: U10  ! The 10 m wind speed [L T-1 ~> m s-1]
  real :: radius    ! The distance from the hurricane center [L ~> m]
  real :: radius10  ! The distance from the hurricane center to its edge [L ~> m]
  real :: radius_km ! The distance from the hurricane center, perhaps in km [L ~> m] or [1000 L ~> km]
  real :: du10 ! The magnitude of the difference between the 10 m wind and the ocean flow [L T-1 ~> m s-1]
  real :: du   ! The difference between the zonal 10 m wind and the zonal ocean flow [L T-1 ~> m s-1]
  real :: dv   ! The difference between the meridional 10 m wind and the zonal ocean flow [L T-1 ~> m s-1]
  real :: Cd   ! The drag coefficient [nondim]
  ! These variables with weird units are only used with pre-20240501 expressions
  real :: radiusB   ! A rescaled radius in m raised to the variable power CS%Holland_B [m^B]
  real :: Holland_A ! Parameter 'A' from the Holland formula, in units of m raised to Holland_B [m^B]
  real :: Holland_AxBxDP ! 'A' x 'B' x (Pressure Ambient-Pressure central)
                         ! for the Holland profile calculation [m^B R L2 T-2 ~> m^B Pa]
  real :: tmp  ! A temporary variable [m^B R L T-1 ~> m^B kg m-2 s-1]
  ! These variables are used with expressions from 20240501 or later
  real :: dP    ! The pressure difference across the hurricane [R L2 T-2 ~> Pa]
  real :: tmpA  ! A temporary variable [R L2 T-2 ~> Pa]
  real :: tmpB  ! A temporary variable [R L T-1 ~> kg m-2 s-1]
  real :: rad_max_rad_B  ! The radius of maximum wind divided by the distance from the center raised
                ! to the power of Holland_B [nondim]
  real :: rad_rad_max    ! The radius normalized by the radius of maximum winds [nondim]

  !Wind angle variables
  real :: Alph ! The wind inflow angle (positive outward) [radians]
  real :: Rstr ! A function of the position normalized by the radius of maximum winds [nondim]
  real :: A0   ! The axisymmetric inflow angle [degrees]
  real :: A1   ! The inflow angle asymmetry [degrees]
  real :: P1   ! The angle difference between the translation direction and the inflow direction [radians]
  real :: Adir ! The angle of the direction from the center to a point [radians]
  real :: V_TS ! Meridional hurricane translation speed [L T-1 ~> m s-1]
  real :: U_TS ! Zonal hurricane translation speed [L T-1 ~> m s-1]

  ! Implementing Holland (1980) parametric wind profile

  radius = SQRT((XX**2) + (YY**2))
  rad_rad_max = radius / CS%rad_max_wind

  ! rkm - r converted to km for Holland prof.
  !       used in km due to error, correct implementation should
  !       not need rkm, but to match winds w/ experiment this must
  !       be maintained.  Causes winds far from storm center to be a
  !       couple of m/s higher than the correct Holland prof.
  if (CS%BR_Bench) then
    radius_km = radius/1000.
  else
    ! if not comparing to benchmark, then use correct Holland prof.
    radius_km = radius
  endif

  !/
  ! Calculate U10 in the interior (inside of the hurricane edge radius),
  ! while adjusting U10 to 0 outside of the ambient wind radius.
  if (CS%answer_date < 20190101) then
    radiusB = (US%L_to_m*radius)**CS%Holland_B
    Holland_A = (US%L_to_m*CS%rad_max_wind)**CS%Holland_B
    if ( (radius > 0.001*CS%rad_max_wind) .and. (radius < CS%rad_edge*CS%rad_max_wind) ) then
      Holland_AxBxDP = Holland_A*CS%Holland_B*(CS%pressure_ambient - CS%pressure_central)
      U10 = sqrt(Holland_AxBxDP*exp(-Holland_A/radiusB) / (CS%rho_a*radiusB) + &
                 0.25*(radius_km*absf)**2) - 0.5*radius_km*absf
    elseif ( (radius > CS%rad_edge*CS%rad_max_wind) .and. (radius < CS%rad_ambient*CS%rad_max_wind) ) then
      if (CS%edge_taper_bug) then  ! This recreates a bug that was in SCM_idealized_hurricane_wind_forcing.
        radius = CS%rad_edge * CS%rad_max_wind
        rad_rad_max = CS%rad_edge
      endif

      radius10 = CS%rad_max_wind*CS%rad_edge
      if (CS%BR_Bench) then
        radius_km = radius10/1000.
      else
        radius_km = radius10
      endif
      radiusB = (US%L_to_m*radius10)**CS%Holland_B

      Holland_AxBxDP = Holland_A*CS%Holland_B*(CS%pressure_ambient - CS%pressure_central)
      U10 = (sqrt(Holland_AxBxDp*exp(-Holland_A/radiusB) / (CS%rho_a*radiusB) + &
                  0.25*(radius_km*absf)**2) - 0.5*radius_km*absf) &
             * (CS%rad_ambient - radius/CS%rad_max_wind) / (CS%rad_ambient - CS%rad_edge)
    else
      U10 = 0.
    endif
  elseif (CS%answer_date < 20240501) then
    ! This is mathematically equivalent to that is above but more accurate.
    radiusB = (US%L_to_m*radius)**CS%Holland_B
    Holland_A = (US%L_to_m*CS%rad_max_wind)**CS%Holland_B
    if ( (radius > 0.001*CS%rad_max_wind) .and. (radius < CS%rad_edge*CS%rad_max_wind) ) then
      Holland_AxBxDP = Holland_A*CS%Holland_B*(CS%pressure_ambient - CS%pressure_central)
      tmp = ( 0.5*radius_km*absf) * (CS%rho_a*radiusB)
      U10 = (Holland_AxBxDP * exp(-Holland_A/radiusB)) / &
            ( tmp + sqrt(Holland_AxBxDP*exp(-Holland_A/radiusB) * (CS%rho_a*radiusB) + tmp**2) )
    elseif ( (radius > CS%rad_edge*CS%rad_max_wind) .and. (radius < CS%rad_ambient*CS%rad_max_wind) ) then
      if (CS%edge_taper_bug) then  ! This recreates a bug that was in SCM_idealized_hurricane_wind_forcing.
        radius = CS%rad_edge * CS%rad_max_wind
        rad_rad_max = CS%rad_edge
      endif

      radius_km = CS%rad_edge * CS%rad_max_wind
      if (CS%BR_Bench) radius_km = radius_km/1000.
      radiusB = (CS%rad_edge*US%L_to_m*CS%rad_max_wind)**CS%Holland_B
      tmp = ( 0.5*radius_km*absf) * (CS%rho_a*radiusB)
      Holland_AxBxDP = Holland_A*CS%Holland_B*(CS%pressure_ambient - CS%pressure_central)
      U10 = ((CS%rad_ambient/(CS%rad_ambient - CS%rad_edge)) - &
             radius/((CS%rad_ambient - CS%rad_edge)*CS%rad_max_wind)) * &
            (Holland_AxBxDp*exp(-Holland_A/radiusB) ) / &
            ( tmp + sqrt(Holland_AxBxDp*exp(-Holland_A/radiusB) * (CS%rho_a*radiusB) + tmp**2) )
    else
      U10 = 0.0
    endif
  else
    ! This is mathematically equivalent to the expressions above, but allows for full
    ! dimensional consistency testing.
    dP = CS%pressure_ambient - CS%pressure_central
    if ( (rad_rad_max > 0.001) .and. (rad_rad_max <= CS%rad_edge) ) then
      rad_max_rad_B = (rad_rad_max)**(-CS%Holland_B)
      tmpA = (rad_max_rad_B*CS%Holland_B) * dp
      tmpB = (0.5*radius_km*absf) * CS%rho_a
      U10 = ( tmpA * exp(-rad_max_rad_B) ) / &
            ( tmpB + sqrt( (tmpA * CS%rho_a) * exp(-rad_max_rad_B) + tmpB**2) )
    elseif ( (rad_rad_max > CS%rad_edge) .and. (rad_rad_max < CS%rad_ambient) ) then
      if (CS%edge_taper_bug) then  ! This recreates a bug that was in SCM_idealized_hurricane_wind_forcing.
        radius = CS%rad_edge * CS%rad_max_wind
        rad_rad_max = CS%rad_edge
      endif

      radius_km = CS%rad_edge * CS%rad_max_wind
      if (CS%BR_Bench) radius_km = radius_km * 0.001
      rad_max_rad_B = CS%rad_edge**(-CS%Holland_B)
      tmpA = (rad_max_rad_B*CS%Holland_B) * dp
      tmpB = (0.5*radius_km*absf) * CS%rho_a
      U10 = ((CS%rad_ambient - rad_rad_max) * ( tmpA * exp(-rad_max_rad_B) )) / &
            ((CS%rad_ambient - CS%rad_edge) * &
             ( tmpB + sqrt((tmpA * CS%rho_a) * exp(-rad_max_rad_B) + tmpB**2) ) )
    else
      U10 = 0.0
    endif
  endif

  Adir = atan2(YY,XX)

  !\

  ! Wind angle model following Zhang and Ulhorn (2012)
  ! ALPH is inflow angle positive outward.
  RSTR = min(CS%rad_edge, rad_rad_max)
  if (CS%answer_date < 20240501) then
    A0 = CS%A0_Rnorm*RSTR + CS%A0_speed*CS%max_windspeed + CS%A0_0
    A1 = -A0*(CS%A1_Rnorm*RSTR + CS%A1_speed*CS%hurr_translation_spd + CS%A1_0)
    P1 = (CS%P1_Rnorm*RSTR + CS%P1_speed*CS%hurr_translation_spd + CS%P1_0) * CS%Deg2Rad
    ALPH = A0 - A1*cos(CS%hurr_translation_dir-Adir-P1)
    if ( (radius > CS%rad_edge*CS%rad_max_wind) .and. (radius < CS%rad_ambient*CS%rad_max_wind) ) then
      ALPH = ALPH*(CS%rad_ambient - rad_rad_max) / (CS%rad_ambient - CS%rad_edge)
    elseif (radius > CS%rad_ambient*CS%rad_max_wind) then  ! This should be >= to avoid a jump at CS%rad_ambient
      ALPH = 0.0
    endif
    ALPH = ALPH * CS%Deg2Rad
  else
    A0 = (CS%A0_Rnorm*RSTR + CS%A0_speed*CS%max_windspeed) + CS%A0_0
    A1 = -A0*((CS%A1_Rnorm*RSTR + CS%A1_speed*CS%hurr_translation_spd) + CS%A1_0)
    P1 = ((CS%P1_Rnorm*RSTR + CS%P1_speed*CS%hurr_translation_spd) + CS%P1_0) * CS%Deg2Rad
    ALPH = (A0 - A1*cos((CS%hurr_translation_dir- Adir) - P1) ) * CS%Deg2Rad
    if (rad_rad_max > CS%rad_edge) &
      ALPH = ALPH * (max(CS%rad_ambient - rad_rad_max, 0.0) / (CS%rad_ambient - CS%rad_edge))
  endif

  ! Calculate translation speed components
  U_TS = CS%hurr_translation_spd * 0.5*cos(CS%hurr_translation_dir)
  V_TS = CS%hurr_translation_spd * 0.5*sin(CS%hurr_translation_dir)

  ! Set output (relative) winds
  dU = U10*sin(Adir-CS%pi-Alph) - Uocn + U_TS
  dV = U10*cos(Adir-Alph) - Vocn + V_TS

  !  Use a simple drag coefficient as a function of U10 (from Sullivan et al., 2010)
  du10 = sqrt((du**2) + (dv**2))
  Cd = simple_wind_scaled_Cd(u10, du10, CS)

  ! Compute stress vector
  TX = US%L_to_Z * CS%rho_a * Cd * du10 * dU
  TY = US%L_to_Z * CS%rho_a * Cd * du10 * dV
end subroutine idealized_hurricane_wind_profile

!> This function returns the air-sea drag coefficient using a simple function of the air-sea velocity difference.
function simple_wind_scaled_Cd(u10, du10, CS) result(Cd)
  real,                      intent(in) :: U10  !< The 10 m wind speed [L T-1 ~> m s-1]
  real,                      intent(in) :: du10 !< The magnitude of the difference between the 10 m wind
                                                !! and the ocean flow [L T-1 ~> m s-1]
  type(idealized_hurricane_CS), pointer :: CS   !< Container for SCM parameters
  real :: Cd  ! Air-sea drag coefficient [nondim]

  ! Note that these expressions are discontinuous at dU10 = 11 and 20 m s-1.
  if (dU10 < CS%calm_speed) then
    Cd = CS%Cd_calm
  elseif (dU10 < CS%windy_speed) then
    if (CS%answer_date < 20190101) then
      Cd = (CS%Cd_intercept + CS%dCd_dU10 * U10 )*0.001
    else
      Cd = (CS%Cd_intercept + CS%dCd_dU10 * dU10 )*0.001
    endif
  else
    Cd = CS%Cd_windy
  endif

end function simple_wind_scaled_Cd

end module idealized_hurricane
