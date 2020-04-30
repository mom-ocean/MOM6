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
!               for multiple example cases)..
! To do
! 1. Remove the legacy SCM_idealized_hurricane_wind_forcing code
! 2. Make the hurricane-to-background wind transition a runtime parameter
!

use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, mech_forcing
use MOM_forcing_type, only : allocate_forcing_type, allocate_mech_forcing
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
public SCM_idealized_hurricane_wind_forcing !Public interface to the legacy idealized
                                        ! hurricane wind profile for SCM.

!> Container for parameters describing idealized wind structure
type, public :: idealized_hurricane_CS ; private

  ! Parameters used to compute Holland radial wind profile
  real    :: rho_a                !< Mean air density [R ~> kg m-3]
  real    :: pressure_ambient     !< Pressure at surface of ambient air [R L2 T-2 ~> Pa]
  real    :: pressure_central     !< Pressure at surface at hurricane center [R L2 T-2 ~> Pa]
  real    :: rad_max_wind         !< Radius of maximum winds [L ~> m]
  real    :: max_windspeed        !< Maximum wind speeds [L T-1 ~> m s-1]
  real    :: hurr_translation_spd !< Hurricane translation speed [L T-1 ~> m s-1]
  real    :: hurr_translation_dir !< Hurricane translation direction [radians]
  real    :: gustiness            !< Gustiness (optional, used in u*) [R L Z T-2 ~> Pa]
  real    :: Rho0                 !< A reference ocean density [R ~> kg m-3]
  real    :: Hurr_cen_Y0          !< The initial y position of the hurricane
                                  !!  This experiment is conducted in a Cartesian
                                  !!  grid and this is assumed to be in meters [L ~> m]
  real    :: Hurr_cen_X0          !< The initial x position of the hurricane
                                  !!  This experiment is conducted in a Cartesian
                                  !!  grid and this is assumed to be in meters [L ~> m]
  real    :: Holland_A            !< Parameter 'A' from the Holland formula [nondim]
  real    :: Holland_B            !< Parameter 'B' from the Holland formula [nondim]
  real    :: Holland_AxBxDP       !< 'A' x 'B' x (Pressure Ambient-Pressure central)
                                  !! for the Holland prorfile calculation [R L2 T-2 ~> Pa]
  logical :: relative_tau         !< A logical to take difference between wind
                                  !! and surface currents to compute the stress
  logical :: answers_2018         !< If true, use expressions driving the idealized hurricane test
                                  !! case that recover the answers from the end of 2018.  Otherwise use
                                  !! expressions that are rescalable and respect rotational symmetry.

  ! Parameters used if in SCM (single column model) mode
  logical :: SCM_mode        !< If true this being used in Single Column Model mode
  logical :: BR_BENCH        !< A "benchmark" configuration (which is meant to
                             !!  provide identical wind to reproduce a previous
                             !!  experiment, where that wind formula contained
                             !!  an error)
  real    :: dy_from_center  !< (Fixed) distance in y from storm center path [L ~> m]

  ! Par
  real :: PI      !< Mathematical constant
  real :: Deg2Rad !< Mathematical constant

end type

! This include declares and sets the variable "version".
#include "version_variable.h"

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
  real :: C
  logical :: default_2018_answers ! The default setting for the various 2018_ANSWERS flags.

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
                 units='Pa', default=101200., scale=US%m_s_to_L_T**2*US%kg_m3_to_R)
  call get_param(param_file, mdl, "IDL_HURR_CENTRAL_PRESSURE", CS%pressure_central, &
                 "Central pressure used in the idealized hurricane wind profile.", &
                 units='Pa', default=96800., scale=US%m_s_to_L_T**2*US%kg_m3_to_R)
  call get_param(param_file, mdl, "IDL_HURR_RAD_MAX_WIND", &
                 CS%rad_max_wind, "Radius of maximum winds used in the "//&
                 "idealized hurricane wind profile.", units='m', &
                 default=50.e3, scale=US%m_to_L)
  call get_param(param_file, mdl, "IDL_HURR_MAX_WIND", CS%max_windspeed, &
                 "Maximum wind speed used in the idealized hurricane"// &
                 "wind profile.", units='m/s', default=65., scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "IDL_HURR_TRAN_SPEED", CS%hurr_translation_spd, &
                 "Translation speed of hurricane used in the idealized "//&
                 "hurricane wind profile.", units='m/s', default=5.0, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "IDL_HURR_TRAN_DIR", CS%hurr_translation_dir, &
                 "Translation direction (towards) of hurricane used in the "//&
                 "idealized hurricane wind profile.", units='degrees', &
                 default=180.0, scale=CS%Deg2Rad)
  call get_param(param_file, mdl, "IDL_HURR_X0", CS%Hurr_cen_X0, &
                 "Idealized Hurricane initial X position", &
                 units='m', default=0., scale=US%m_to_L)
  call get_param(param_file, mdl, "IDL_HURR_Y0", CS%Hurr_cen_Y0, &
                 "Idealized Hurricane initial Y position", &
                 units='m', default=0., scale=US%m_to_L)
  call get_param(param_file, mdl, "IDL_HURR_TAU_CURR_REL", CS%relative_tau, &
                 "Current relative stress switch "//&
                 "used in the idealized hurricane wind profile.", &
                 units='', default=.false.)

  ! Parameters for SCM mode
  call get_param(param_file, mdl, "IDL_HURR_SCM_BR_BENCH", CS%BR_BENCH, &
                 "Single column mode benchmark case switch, which is "// &
                 "invoking a modification (bug) in the wind profile meant to "//&
                 "reproduce a previous implementation.", units='', default=.false.)
  call get_param(param_file, mdl, "IDL_HURR_SCM", CS%SCM_MODE, &
                 "Single Column mode switch "//&
                 "used in the SCM idealized hurricane wind profile.", &
                 units='', default=.false.)
  call get_param(param_file, mdl, "IDL_HURR_SCM_LOCY", CS%dy_from_center, &
                 "Y distance of station used in the SCM idealized hurricane "//&
                 "wind profile.", units='m', default=50.e3, scale=US%m_to_L)
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.true.)
  call get_param(param_file, mdl, "IDL_HURR_2018_ANSWERS", CS%answers_2018, &
                 "If true, use expressions driving the idealized hurricane test case that recover "//&
                 "the answers from the end of 2018.  Otherwise use expressions that are rescalable "//&
                 "and respect rotational symmetry.", default=default_2018_answers)

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
                 "The background gustiness in the winds.", units="Pa", &
                 default=0.0, scale=US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z, do_not_log=.true.)


  if (CS%BR_BENCH) then
    CS%rho_a = 1.2*US%kg_m3_to_R
  endif
  dP = CS%pressure_ambient - CS%pressure_central
  if (CS%answers_2018) then
    C = CS%max_windspeed / sqrt( US%R_to_kg_m3 * dP )
    CS%Holland_B = C**2 * US%R_to_kg_m3*CS%rho_a * exp(1.0)
  else
    CS%Holland_B = CS%max_windspeed**2 * CS%rho_a * exp(1.0) / dP
  endif
  CS%Holland_A = (US%L_to_m*CS%rad_max_wind)**CS%Holland_B
  CS%Holland_AxBxDP = CS%Holland_A*CS%Holland_B*dP

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

  ! Allocate the forcing arrays, if necessary.
  call allocate_mech_forcing(G, forces, stress=.true., ustar=.true.)

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
    fbench = 5.5659e-05 * US%T_to_s
    fbench_fac = 0.0
  else
    fbench = 0.0
    fbench_fac = 1.0
  endif

  !> Computes taux
  do j=js,je
    do I=is-1,Ieq
      Uocn = sfc_state%u(I,j) * REL_TAU_FAC
      if (CS%answers_2018) then
        Vocn = 0.25*(sfc_state%v(i,J)+sfc_state%v(i+1,J-1)&
                    +sfc_state%v(i+1,J)+sfc_state%v(i,J-1))*REL_TAU_FAC
      else
        Vocn =0.25*((sfc_state%v(i,J)+sfc_state%v(i+1,J-1)) +&
                    (sfc_state%v(i+1,J)+sfc_state%v(i,J-1))) * REL_TAU_FAC
      endif
      f_local = abs(0.5*(G%CoriolisBu(I,J)+G%CoriolisBu(I,J-1)))*fbench_fac + fbench
      ! Calculate position as a function of time.
      if (CS%SCM_mode) then
        YY = YC + CS%dy_from_center
        XX = XC
      else
        YY = G%geoLatCu(I,j)*1000.*US%m_to_L - YC
        XX = G%geoLonCu(I,j)*1000.*US%m_to_L - XC
      endif
      call idealized_hurricane_wind_profile(CS, US, f_local, YY, XX, Uocn, Vocn, TX, TY)
      forces%taux(I,j) = G%mask2dCu(I,j) * TX
    enddo
  enddo
  !> Computes tauy
  do J=js-1,Jeq
    do i=is,ie
      if (CS%answers_2018) then
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
        YY = YC + CS%dy_from_center
        XX = XC
      else
        YY = G%geoLatCv(i,J)*1000.*US%m_to_L - YC
        XX = G%geoLonCv(i,J)*1000.*US%m_to_L - XC
      endif
      call idealized_hurricane_wind_profile(CS, US, f_local, YY, XX, Uocn, Vocn, TX, TY)
      forces%tauy(i,J) = G%mask2dCv(i,J) * TY
    enddo
  enddo

  !> Get Ustar
  do j=js,je
    do i=is,ie
      !  This expression can be changed if desired, but need not be.
      forces%ustar(i,j) = G%mask2dT(i,j) * sqrt(US%L_to_Z * (CS%gustiness/CS%Rho0 + &
              sqrt(0.5*(forces%taux(I-1,j)**2 + forces%taux(I,j)**2) + &
                   0.5*(forces%tauy(i,J-1)**2 + forces%tauy(i,J)**2))/CS%Rho0))
    enddo
  enddo

  return
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
  real :: radius10  ! 10 times the distance from the hurricane center [L ~> m]
  real :: radius_km ! The distance from the hurricane center, perhaps in km [L ~> m] or [1000 L ~> km]
  real :: radiusB
  real :: tmp  ! A temporary variable [R L T-1 ~> kg m-2 s-1]
  real :: du10 ! The magnitude of the difference between the 10 m wind and the ocean flow [L T-1 ~> m s-1]
  real :: du   ! The difference between the zonal 10 m wind and the zonal ocean flow [L T-1 ~> m s-1]
  real :: dv   ! The difference between the meridional 10 m wind and the zonal ocean flow [L T-1 ~> m s-1]
  real :: CD

  !Wind angle variables
  real :: Alph !< The resulting inflow angle (positive outward)
  real :: Rstr
  real :: A0
  real :: A1
  real :: P1
  real :: Adir
  real :: V_TS ! Meridional hurricane translation speed [L T-1 ~> m s-1]
  real :: U_TS ! Zonal hurricane translation speed [L T-1 ~> m s-1]

  ! Implementing Holland (1980) parameteric wind profile

  radius = SQRT(XX**2 + YY**2)

  !/ BGR
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
  radiusB = (US%L_to_m*radius)**CS%Holland_B

  !/
  ! Calculate U10 in the interior (inside of 10x radius of maximum wind),
  ! while adjusting U10 to 0 outside of 12x radius of maximum wind.
  if (CS%answers_2018) then
    if ( (radius > 0.001*CS%rad_max_wind) .and. (radius < 10.*CS%rad_max_wind) ) then
      U10 = sqrt(CS%Holland_AxBxDP*exp(-CS%Holland_A/radiusB) / (CS%rho_a*radiusB) + &
                 0.25*(radius_km*absf)**2) - 0.5*radius_km*absf
    elseif ( (radius > 10.*CS%rad_max_wind) .and. (radius < 15.*CS%rad_max_wind) ) then
      radius10 = CS%rad_max_wind*10.
      if (CS%BR_Bench) then
        radius_km = radius10/1000.
      else
        radius_km = radius10
      endif
      radiusB = (US%L_to_m*radius10)**CS%Holland_B

      U10 = (sqrt(CS%Holland_AxBxDp*exp(-CS%Holland_A/radiusB) / (CS%rho_a*radiusB) + &
                  0.25*(radius_km*absf)**2) - 0.5*radius_km*absf) &
             * (15. - radius/CS%rad_max_wind)/5.
    else
      U10 = 0.
    endif
  else  ! This is mathematically equivalent to that is above but more accurate.
    if ( (radius > 0.001*CS%rad_max_wind) .and. (radius < 10.*CS%rad_max_wind) ) then
      tmp = ( 0.5*radius_km*absf) * (CS%rho_a*radiusB)
      U10 = (CS%Holland_AxBxDP * exp(-CS%Holland_A/radiusB)) / &
            ( tmp + sqrt(CS%Holland_AxBxDP*exp(-CS%Holland_A/radiusB) * (CS%rho_a*radiusB) + tmp**2) )
    elseif ( (radius > 10.*CS%rad_max_wind) .and. (radius < 15.*CS%rad_max_wind) ) then
      radius_km = 10.0 * CS%rad_max_wind
      if (CS%BR_Bench) radius_km = radius_km/1000.
      radiusB = (10.0*US%L_to_m*CS%rad_max_wind)**CS%Holland_B
      tmp = ( 0.5*radius_km*absf) * (CS%rho_a*radiusB)
      U10 = (3.0 - radius/(5.0*CS%rad_max_wind)) * (CS%Holland_AxBxDp*exp(-CS%Holland_A/radiusB) ) / &
            ( tmp + sqrt(CS%Holland_AxBxDp*exp(-CS%Holland_A/radiusB) * (CS%rho_a*radiusB) + tmp**2) )
    else
      U10 = 0.0
    endif
  endif

  Adir = atan2(YY,XX)

  !\

  ! Wind angle model following Zhang and Ulhorn (2012)
  ! ALPH is inflow angle positive outward.
  RSTR = min(10., radius / CS%rad_max_wind)
  A0 = -0.9*RSTR - 0.09*US%L_T_to_m_s*CS%max_windspeed - 14.33
  A1 = -A0*(0.04*RSTR + 0.05*US%L_T_to_m_s*CS%hurr_translation_spd + 0.14)
  P1 = (6.88*RSTR - 9.60*US%L_T_to_m_s*CS%hurr_translation_spd + 85.31) * CS%Deg2Rad
  ALPH = A0 - A1*cos(CS%hurr_translation_dir-Adir-P1)
  if ( (radius > 10.*CS%rad_max_wind) .and.&
       (radius < 15.*CS%rad_max_wind) ) then
     ALPH = ALPH*(15.0 - radius/CS%rad_max_wind)/5.
  elseif (radius > 15.*CS%rad_max_wind) then
     ALPH = 0.0
  endif
  ALPH = ALPH * CS%Deg2Rad

  ! Calculate translation speed components
  U_TS = CS%hurr_translation_spd * 0.5*cos(CS%hurr_translation_dir)
  V_TS = CS%hurr_translation_spd * 0.5*sin(CS%hurr_translation_dir)

  ! Set output (relative) winds
  dU = U10*sin(Adir-CS%Pi-Alph) - Uocn + U_TS
  dV = U10*cos(Adir-Alph) - Vocn + V_TS

  !  Use a simple drag coefficient as a function of U10 (from Sullivan et al., 2010)
  du10 = sqrt(du**2+dv**2)
  if (dU10 < 11.0*US%m_s_to_L_T) then
    Cd = 1.2e-3
  elseif (dU10 < 20.0*US%m_s_to_L_T) then
    if (CS%answers_2018) then
      Cd = (0.49 + 0.065*US%L_T_to_m_s*U10)*1.e-3
    else
      Cd = (0.49 + 0.065*US%L_T_to_m_s*dU10)*1.e-3
    endif
  else
    Cd = 1.8e-3
  endif

  ! Compute stress vector
  TX = US%L_to_Z * CS%rho_a * Cd * sqrt(dU**2 + dV**2) * dU
  TY = US%L_to_Z * CS%rho_a * Cd * sqrt(dU**2 + dV**2) * dV

end subroutine idealized_hurricane_wind_profile

!> This subroutine is primarily needed as a legacy for reproducing answers.
!! It is included as an additional subroutine rather than padded into the previous
!! routine with flags to ease its eventual removal.  Its functionality is replaced
!! with the new routines and it can be deleted when answer changes are acceptable.
subroutine SCM_idealized_hurricane_wind_forcing(sfc_state, forces, day, G, US, CS)
  type(surface),                intent(in)    :: sfc_state  !< Surface state structure
  type(mech_forcing),           intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),              intent(in)    :: day    !< Time in days
  type(ocean_grid_type),        intent(inout) :: G      !< Grid structure
  type(unit_scale_type),        intent(in)    :: US     !< A dimensional unit scaling type
  type(idealized_hurricane_CS), pointer       :: CS     !< Container for SCM parameters
  ! Local variables
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real :: pie, Deg2Rad
  real :: du10 ! The magnitude of the difference between the 10 m wind and the ocean flow [L T-1 ~> m s-1]
  real :: U10  ! The 10 m wind speed [L T-1 ~> m s-1]
  real :: A, B, C ! For wind profile expression
  real :: rad  ! The distance from the hurricane center [L ~> m]
  real :: rkm  ! The distance from the hurricane center, sometimes scaled to km [L ~> m] or [1000 L ~> km]
  real :: f_local  ! The local Coriolis parameter [T-1 ~> s-1]
  real :: xx  ! x-position [L ~> m]
  real :: t0 !for location
  real :: dP  ! The pressure difference across the hurricane [R L2 T-2 ~> Pa]
  real :: rB
  real :: Cd ! Air-sea drag coefficient
  real :: Uocn, Vocn ! Surface ocean velocity components [L T-1 ~> m s-1]
  real :: dU, dV ! Air-sea differential motion [L T-1 ~> m s-1]
  !Wind angle variables
  real :: Alph,Rstr, A0, A1, P1, Adir, transdir
  real :: V_TS, U_TS ! Components of the translation speed [L T-1 ~> m s-1]
  logical :: BR_Bench
  ! Bounds for loops and memory allocation
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! Allocate the forcing arrays, if necessary.

  call allocate_mech_forcing(G, forces, stress=.true., ustar=.true.)
  pie = 4.0*atan(1.0) ; Deg2Rad = pie/180.
  !/ BR
  ! Implementing Holland (1980) parameteric wind profile
  !------------------------------------------------------|
  BR_Bench = .true.   !true if comparing to LES runs     |
  t0 = 129600.        !TC 'eye' crosses (0,0) at 36 hours|
  transdir = pie      !translation direction (-x)        |
  !------------------------------------------------------|
  dP = CS%pressure_ambient - CS%pressure_central
  if (CS%answers_2018) then
    C = CS%max_windspeed / sqrt( US%R_to_kg_m3*dP )
    B = C**2 * US%R_to_kg_m3*CS%rho_a * exp(1.0)
    if (BR_Bench) then ! rho_a reset to value used in generated wind for benchmark test
       B = C**2 * 1.2 * exp(1.0)
    endif
  elseif (BR_Bench) then ! rho_a reset to value used in generated wind for benchmark test
    B = (CS%max_windspeed**2 / dP ) * 1.2*US%kg_m3_to_R * exp(1.0)
  else
    B = (CS%max_windspeed**2 /dP ) * CS%rho_a * exp(1.0)
  endif

  A = (US%L_to_m*CS%rad_max_wind / 1000.)**B
  f_local = G%CoriolisBu(is,js) ! f=f(x,y) but in the SCM is constant
  if (BR_Bench) then
    ! f reset to value used in generated wind for benchmark test
    f_local = 5.5659e-05*US%T_to_s
  endif
  !/ BR
  ! Calculate x position as a function of time.
  xx = US%s_to_T*( t0 - time_type_to_real(day)) * CS%hurr_translation_spd * cos(transdir)
  rad = sqrt(xx**2 + CS%dy_from_center**2)
  !/ BR
  ! rkm - rad converted to km for Holland prof.
  !       used in km due to error, correct implementation should
  !       not need rkm, but to match winds w/ experiment this must
  !       be maintained.  Causes winds far from storm center to be a
  !       couple of m/s higher than the correct Holland prof.
  if (BR_Bench) then
     rkm = rad/1000.
     rB = (US%L_to_m*rkm)**B
  else
     ! if not comparing to benchmark, then use correct Holland prof.
     rkm = rad
     rB = (US%L_to_m*rad)**B
  endif
  !/ BR
  ! Calculate U10 in the interior (inside of 10x radius of maximum wind),
  ! while adjusting U10 to 0 outside of 12x radius of maximum wind.
  ! Note that rho_a is set to 1.2 following generated wind for experiment
  if (rad > 0.001*CS%rad_max_wind .AND. rad < 10.*CS%rad_max_wind) then
    U10 = sqrt( A*B*dP*exp(-A/rB)/(1.2*US%kg_m3_to_R*rB) + 0.25*(rkm*f_local)**2 ) - 0.5*rkm*f_local
  elseif (rad > 10.*CS%rad_max_wind .AND. rad < 12.*CS%rad_max_wind) then
    rad=(CS%rad_max_wind)*10.
    if (BR_Bench) then
       rkm = rad/1000.
       rB = (US%L_to_m*rkm)**B
    else
       rkm = rad
       rB = (US%L_to_m*rad)**B
    endif
    U10 = ( sqrt( A*B*dP*exp(-A/rB)/(1.2*US%kg_m3_to_R*rB) + 0.25*(rkm*f_local)**2 ) - 0.5*rkm*f_local) &
          * (12. - rad/CS%rad_max_wind)/2.
  else
    U10 = 0.
  endif
  Adir = atan2(CS%dy_from_center,xx)

  !/ BR
  ! Wind angle model following Zhang and Ulhorn (2012)
  ! ALPH is inflow angle positive outward.
  RSTR = min(10., rad / CS%rad_max_wind)
  A0 = -0.9*RSTR - 0.09*US%L_T_to_m_s*CS%max_windspeed - 14.33
  A1 = -A0 *(0.04*RSTR + 0.05*US%L_T_to_m_s*CS%hurr_translation_spd + 0.14)
  P1 = (6.88*RSTR - 9.60*US%L_T_to_m_s*CS%hurr_translation_spd + 85.31)*pie/180.
  ALPH = A0 - A1*cos( (TRANSDIR - ADIR ) - P1)
  if (rad > 10.*CS%rad_max_wind .AND. rad < 12.*CS%rad_max_wind) then
    ALPH = ALPH* (12. - rad/CS%rad_max_wind)/2.
  elseif (rad > 12.*CS%rad_max_wind) then
    ALPH = 0.0
  endif
  ALPH = ALPH * Deg2Rad
 !/BR
  ! Prepare for wind calculation
  ! X_TS is component of translation speed added to wind vector
  ! due to background steering wind.
  U_TS = CS%hurr_translation_spd*0.5*cos(transdir)
  V_TS = CS%hurr_translation_spd*0.5*sin(transdir)

  ! Set the surface wind stresses, in [Pa]. A positive taux
  ! accelerates the ocean to the (pseudo-)east.
  !   The i-loop extends to is-1 so that taux can be used later in the
  ! calculation of ustar - otherwise the lower bound would be Isq.
  do j=js,je ; do I=is-1,Ieq
    !/BR
    ! Turn off surface current for stress calculation to be
    ! consistent with test case.
    Uocn = 0. ! sfc_state%u(I,j)
    Vocn = 0. ! 0.25*( (sfc_state%v(i,J) + sfc_state%v(i+1,J-1)) + &
              !        (sfc_state%v(i+1,J) + sfc_state%v(i,J-1)) )
    !/BR
    ! Wind vector calculated from location/direction (sin/cos flipped b/c
    ! cyclonic wind is 90 deg. phase shifted from position angle).
    dU = U10*sin(Adir-pie-Alph) - Uocn + U_TS
    dV = U10*cos(Adir-Alph) - Vocn + V_TS
    !/----------------------------------------------------|
    !BR
    !  Add a simple drag coefficient as a function of U10 |
    !/----------------------------------------------------|
    du10 = sqrt(du**2+dv**2)
    if (dU10 < 11.0*US%m_s_to_L_T) then
      Cd = 1.2e-3
    elseif (dU10 < 20.0*US%m_s_to_L_T) then
      if (CS%answers_2018) then
        Cd = (0.49 + 0.065 * US%L_T_to_m_s*U10 )*0.001
      else
        Cd = (0.49 + 0.065 * US%L_T_to_m_s*dU10 )*0.001
      endif
    else
      Cd = 0.0018
    endif
    forces%taux(I,j) = CS%rho_a * US%L_to_Z * G%mask2dCu(I,j) * Cd*du10*dU
  enddo ; enddo
  !/BR
  ! See notes above
  do J=js-1,Jeq ; do i=is,ie
    Uocn = 0. ! 0.25*( (sfc_state%u(I,j) + sfc_state%u(I-1,j+1)) + &
              !        (sfc_state%u(I-1,j) + sfc_state%u(I,j+1)) )
    Vocn = 0. ! sfc_state%v(i,J)
    dU = U10*sin(Adir-pie-Alph) - Uocn + U_TS
    dV = U10*cos(Adir-Alph) - Vocn + V_TS
    du10=sqrt(du**2+dv**2)
    if (dU10 < 11.0*US%m_s_to_L_T) then
      Cd = 1.2e-3
    elseif (dU10 < 20.0*US%m_s_to_L_T) then
      if (CS%answers_2018) then
        Cd = (0.49 + 0.065 * US%L_T_to_m_s*U10 )*0.001
      else
        Cd = (0.49 + 0.065 * US%L_T_to_m_s*dU10 )*0.001
      endif
    else
      Cd = 0.0018
    endif
    forces%tauy(I,j) = CS%rho_a * US%L_to_Z * G%mask2dCv(I,j) * Cd*dU10*dV
  enddo ; enddo
  ! Set the surface friction velocity [Z T-1 ~> m s-1]. ustar is always positive.
  do j=js,je ; do i=is,ie
    !  This expression can be changed if desired, but need not be.
    forces%ustar(i,j) = G%mask2dT(i,j) * sqrt(US%L_to_Z * (CS%gustiness/CS%Rho0 + &
            sqrt(0.5*(forces%taux(I-1,j)**2 + forces%taux(I,j)**2) + &
                 0.5*(forces%tauy(i,J-1)**2 + forces%tauy(i,J)**2))/CS%Rho0))
  enddo ; enddo

end subroutine SCM_idealized_hurricane_wind_forcing

end module idealized_hurricane
