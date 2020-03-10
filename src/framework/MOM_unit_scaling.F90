!> Provides a transparent unit rescaling type to facilitate dimensional consistency testing
module MOM_unit_scaling

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type

implicit none ; private

public unit_scaling_init, unit_scaling_end, fix_restart_unit_scaling

!> Describes various unit conversion factors
type, public :: unit_scale_type
  real :: m_to_Z !< A constant that translates distances in meters to the units of depth.
  real :: Z_to_m !< A constant that translates distances in the units of depth to meters.
  real :: m_to_L !< A constant that translates lengths in meters to the units of horizontal lengths.
  real :: L_to_m !< A constant that translates lengths in the units of horizontal lengths to meters.
  real :: s_to_T !< A constant that translates time intervals in seconds to the units of time.
  real :: T_to_s !< A constant that translates the units of time to seconds.
  real :: R_to_kg_m3 !< A constant that translates the units of density to kilograms per meter cubed.
  real :: kg_m3_to_R !< A constant that translates kilograms per meter cubed to the units of density.
  real :: Q_to_J_kg  !< A constant that translates the units of enthalpy to Joules per kilogram.
  real :: J_kg_to_Q  !< A constant that translates Joules per kilogram to the units of enthalpy.

  ! These are useful combinations of the fundamental scale conversion factors above.
  real :: Z_to_L !< Convert vertical distances to lateral lengths
  real :: L_to_Z !< Convert vertical distances to lateral lengths
  real :: L_T_to_m_s !< Convert lateral velocities from L T-1 to m s-1.
  real :: m_s_to_L_T !< Convert lateral velocities from m s-1 to L T-1.
  real :: L_T2_to_m_s2 !< Convert lateral accelerations from L T-2 to m s-2.
  real :: Z2_T_to_m2_s !< Convert vertical diffusivities from Z2 T-1 to m2 s-1.
  real :: m2_s_to_Z2_T !< Convert vertical diffusivities from m2 s-1 to Z2 T-1.
  real :: W_m2_to_QRZ_T !< Convert heat fluxes from W m-2 to Q R Z T-1.
  real :: QRZ_T_to_W_m2 !< Convert heat fluxes from Q R Z T-1 to W m-2.
  ! Not used enough:  real :: kg_m2_to_RZ   !< Convert mass loads from kg m-2 to R Z.
  real :: RZ_to_kg_m2   !< Convert mass loads from R Z to kg m-2.
  real :: kg_m2s_to_RZ_T  !< Convert mass fluxes from kg m-2 s-1 to R Z T-1.
  real :: RZ_T_to_kg_m2s  !< Convert mass fluxes from R Z T-1 to kg m-2 s-1.
  real :: RZ3_T3_to_W_m2  !< Convert turbulent kinetic energy fluxes from R Z3 T-3 to W m-2.

  ! These are used for changing scaling across restarts.
  real :: m_to_Z_restart = 0.0 !< A copy of the m_to_Z that is used in restart files.
  real :: m_to_L_restart = 0.0 !< A copy of the m_to_L that is used in restart files.
  real :: s_to_T_restart = 0.0 !< A copy of the s_to_T that is used in restart files.
  real :: kg_m3_to_R_restart = 0.0 !< A copy of the kg_m3_to_R that is used in restart files.
  real :: J_kg_to_Q_restart = 0.0 !< A copy of the J_kg_to_Q that is used in restart files.
end type unit_scale_type

contains

!> Allocates and initializes the ocean model unit scaling type
subroutine unit_scaling_init( param_file, US )
  type(param_file_type), intent(in) :: param_file !< Parameter file handle/type
  type(unit_scale_type), pointer    :: US         !< A dimensional unit scaling type

  ! This routine initializes a unit_scale_type structure (US).

  ! Local variables
  integer :: Z_power, L_power, T_power, R_power, Q_power
  real    :: Z_rescale_factor, L_rescale_factor, T_rescale_factor, R_rescale_factor, Q_rescale_factor
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=16) :: mdl = "MOM_unit_scaling"

  if (associated(US)) call MOM_error(FATAL, &
     'unit_scaling_init: called with an associated US pointer.')
  allocate(US)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, &
                   "Parameters for doing unit scaling of variables.")
  call get_param(param_file, mdl, "Z_RESCALE_POWER", Z_power, &
                 "An integer power of 2 that is used to rescale the model's "//&
                 "intenal units of depths and heights.  Valid values range from -300 to 300.", &
                 units="nondim", default=0, debuggingParam=.true.)
  call get_param(param_file, mdl, "L_RESCALE_POWER", L_power, &
                 "An integer power of 2 that is used to rescale the model's "//&
                 "intenal units of lateral distances.  Valid values range from -300 to 300.", &
                 units="nondim", default=0, debuggingParam=.true.)
  call get_param(param_file, mdl, "T_RESCALE_POWER", T_power, &
                 "An integer power of 2 that is used to rescale the model's "//&
                 "intenal units of time.  Valid values range from -300 to 300.", &
                 units="nondim", default=0, debuggingParam=.true.)
  call get_param(param_file, mdl, "R_RESCALE_POWER", R_power, &
                 "An integer power of 2 that is used to rescale the model's "//&
                 "intenal units of density.  Valid values range from -300 to 300.", &
                 units="nondim", default=0, debuggingParam=.true.)
  call get_param(param_file, mdl, "Q_RESCALE_POWER", Q_power, &
                 "An integer power of 2 that is used to rescale the model's "//&
                 "intenal units of heat content.  Valid values range from -300 to 300.", &
                 units="nondim", default=0, debuggingParam=.true.)
  if (abs(Z_power) > 300) call MOM_error(FATAL, "unit_scaling_init: "//&
                 "Z_RESCALE_POWER is outside of the valid range of -300 to 300.")
  if (abs(L_power) > 300) call MOM_error(FATAL, "unit_scaling_init: "//&
                 "L_RESCALE_POWER is outside of the valid range of -300 to 300.")
  if (abs(T_power) > 300) call MOM_error(FATAL, "unit_scaling_init: "//&
                 "T_RESCALE_POWER is outside of the valid range of -300 to 300.")
  if (abs(R_power) > 300) call MOM_error(FATAL, "unit_scaling_init: "//&
                 "R_RESCALE_POWER is outside of the valid range of -300 to 300.")
  if (abs(Q_power) > 300) call MOM_error(FATAL, "unit_scaling_init: "//&
                 "Q_RESCALE_POWER is outside of the valid range of -300 to 300.")

  Z_rescale_factor = 1.0
  if (Z_power /= 0) Z_rescale_factor = 2.0**Z_power
  US%Z_to_m = 1.0 * Z_rescale_factor
  US%m_to_Z = 1.0 / Z_rescale_factor

  L_rescale_factor = 1.0
  if (L_power /= 0) L_rescale_factor = 2.0**L_power
  US%L_to_m = 1.0 * L_rescale_factor
  US%m_to_L = 1.0 / L_rescale_factor

  T_rescale_factor = 1.0
  if (T_power /= 0) T_rescale_factor = 2.0**T_power
  US%T_to_s = 1.0 * T_rescale_factor
  US%s_to_T = 1.0 / T_rescale_factor

  R_rescale_factor = 1.0
  if (R_power /= 0) R_rescale_factor = 2.0**R_power
  US%R_to_kg_m3 = 1.0 * R_rescale_factor
  US%kg_m3_to_R = 1.0 / R_rescale_factor

  Q_Rescale_factor = 1.0
  if (Q_power /= 0) Q_Rescale_factor = 2.0**Q_power
  US%Q_to_J_kg = 1.0 * Q_Rescale_factor
  US%J_kg_to_Q = 1.0 / Q_Rescale_factor

  ! These are useful combinations of the fundamental scale conversion factors set above.
  US%Z_to_L = US%Z_to_m * US%m_to_L
  US%L_to_Z = US%L_to_m * US%m_to_Z
  US%L_T_to_m_s = US%L_to_m * US%s_to_T
  US%m_s_to_L_T = US%m_to_L * US%T_to_s
  US%L_T2_to_m_s2 = US%L_to_m * US%s_to_T**2
  ! It does not look like US%m_s2_to_L_T2 would be used, so it does not exist.
  US%Z2_T_to_m2_s = US%Z_to_m**2 * US%s_to_T
  US%m2_s_to_Z2_T = US%m_to_Z**2 * US%T_to_s
  ! It does not seem like US%kg_m2_to_RZ would be used enough in MOM6 to justify its existence.
  US%RZ_to_kg_m2  = US%R_to_kg_m3 * US%Z_to_m
  US%kg_m2s_to_RZ_T = US%kg_m3_to_R * US%m_to_Z * US%T_to_s
  US%RZ_T_to_kg_m2s = US%R_to_kg_m3 * US%Z_to_m * US%s_to_T
  US%RZ3_T3_to_W_m2 = US%R_to_kg_m3 * US%Z_to_m**3 * US%s_to_T**3
  US%W_m2_to_QRZ_T = US%J_kg_to_Q * US%kg_m3_to_R * US%m_to_Z * US%T_to_s
  US%QRZ_T_to_W_m2 = US%Q_to_J_kg * US%R_to_kg_m3 * US%Z_to_m * US%s_to_T

end subroutine unit_scaling_init

!> Set the unit scaling factors for output to restart files to the unit scaling
!! factors for this run.
subroutine fix_restart_unit_scaling(US)
  type(unit_scale_type), intent(inout) :: US !< A dimensional unit scaling type

  US%m_to_Z_restart = US%m_to_Z
  US%m_to_L_restart = US%m_to_L
  US%s_to_T_restart = US%s_to_T
  US%kg_m3_to_R_restart = US%kg_m3_to_R
  US%J_kg_to_Q_restart = US%J_kg_to_Q

end subroutine fix_restart_unit_scaling

!> Deallocates a unit scaling structure.
subroutine unit_scaling_end( US )
  type(unit_scale_type), pointer :: US !< A dimensional unit scaling type

  deallocate( US )

end subroutine unit_scaling_end

end module MOM_unit_scaling
