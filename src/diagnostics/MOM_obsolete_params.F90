!> Methods for testing for, and list of, obsolete run-time parameters.
module MOM_obsolete_params

! This file is part of MOM6. See LICENSE.md for the license.
! This module was first conceived and written by Robert Hallberg, July 2010.

use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : read_param, log_version, param_file_type

implicit none ; private

#include <MOM_memory.h>

public find_obsolete_params
public obsolete_logical, obsolete_int, obsolete_real, obsolete_char

contains

!> Scans input parameter file for list obsolete parameters.
subroutine find_obsolete_params(param_file)
  type(param_file_type), intent(in) :: param_file !< Structure containing parameter file data.
  ! Local variables
  character(len=40)  :: mdl = "find_obsolete_params" ! This module's name.
! This include declares and sets the variable "version".
#include "version_variable.h"
  integer :: l_seg, nseg
  logical :: test_logic, split
  character(len=40)  :: temp_string

  if (.not.is_root_pe()) return

  call obsolete_logical(param_file, "BLOCKED_ANALYTIC_FV_PGF", &
       hint="BLOCKED_ANALYTIC_FV_PGF is no longer available.")

  call obsolete_logical(param_file, "ADD_KV_SLOW", &
       hint="This option is no longer needed, nor supported.")

  call obsolete_char(param_file, "OBC_CONFIG", &
       hint="Instead use OBC_USER_CONFIG and use the new segments protocol.")
  call obsolete_char(param_file, "READ_OBC_ETA", &
       hint="Instead use OBC_SEGMENT_XXX_DATA.")
  call obsolete_char(param_file, "READ_OBC_UV", &
       hint="Instead use OBC_SEGMENT_XXX_DATA.")
  call obsolete_char(param_file, "READ_OBC_TS", &
       hint="Instead use OBC_SEGMENT_XXX_DATA.")
  call obsolete_char(param_file, "EXTEND_OBC_SEGMENTS", &
       hint="This option is no longer needed, nor supported.")
  call obsolete_char(param_file, "MEKE_VISCOSITY_COEFF", &
       hint="This option has been replaced by MEKE_VISCOSITY_COEFF_KU and \n" //&
            " MEKE_VISCOSITY_COEFF_AU. Please set these parameters instead.")
  nseg = 0
  call read_param(param_file, "OBC_NUMBER_OF_SEGMENTS", nseg)
  do l_seg = 1,nseg
    write(temp_string(1:22),"('OBC_SEGMENT_',i3.3,'_TNUDGE')") l_seg
    call obsolete_real(param_file, temp_string, &
         hint="Instead use OBC_SEGMENT_xxx_VELOCITY_NUDGING_TIMESCALES.")
  enddo

  call obsolete_logical(param_file, "CONVERT_THICKNESS_UNITS", .true.)
  call obsolete_logical(param_file, "MASK_MASSLESS_TRACERS", .false.)

  call obsolete_logical(param_file, "SALT_REJECT_BELOW_ML", .false.)
  call obsolete_logical(param_file, "MLE_USE_MLD_AVE_BUG", .false.)
  call obsolete_logical(param_file, "CORRECT_DENSITY", .true.)
  call obsolete_char(param_file, "WINDSTRESS_STAGGER", warning_val="C", &
                     hint="Use WIND_STAGGER instead.")

  call obsolete_char(param_file, "DIAG_REMAP_Z_GRID_DEF", &
                     hint="Use NUM_DIAG_COORDS, DIAG_COORDS and DIAG_COORD_DEF_Z")

  call obsolete_real(param_file, "VSTAR_SCALE_FACTOR", hint="Use EPBL_VEL_SCALE_FACTOR instead.")

  call obsolete_real(param_file, "VSTAR_SCALE_COEF")
  call obsolete_real(param_file, "ZSTAR_RIGID_SURFACE_THRESHOLD")
  call obsolete_logical(param_file, "HENYEY_IGW_BACKGROUND_NEW")

  call obsolete_real(param_file, "SLIGHT_DZ_SURFACE")
  call obsolete_int(param_file, "SLIGHT_NZ_SURFACE_FIXED")
  call obsolete_real(param_file, "SLIGHT_SURFACE_AVG_DEPTH")
  call obsolete_real(param_file, "SLIGHT_NLAY_TO_INTERIOR")
  call obsolete_logical(param_file, "SLIGHT_FIX_HALOCLINES")
  call obsolete_real(param_file, "HALOCLINE_FILTER_LENGTH")
  call obsolete_real(param_file, "HALOCLINE_STRAT_TOL")

  ! Test for inconsistent parameter settings.
  split = .true. ; test_logic = .false.
  call read_param(param_file,"SPLIT",split)
  call read_param(param_file,"DYNAMIC_SURFACE_PRESSURE",test_logic)
  if (test_logic .and. .not.split) call MOM_ERROR(FATAL, &
    "find_obsolete_params: #define DYNAMIC_SURFACE_PRESSURE is not yet "//&
    "implemented without #define SPLIT.")
  call obsolete_char(param_file, "CONTINUITY_SCHEME", warning_val="PPM", &
                     hint="Only one continuity scheme is available so this need not be specified.")
  call obsolete_real(param_file, "ETA_TOLERANCE_AUX", only_warn=.true.)
  call obsolete_real(param_file, "BT_MASS_SOURCE_LIMIT", 0.0)
  call obsolete_real(param_file, "FIRST_GUESS_SURFACE_LAYER_DEPTH")
  call obsolete_logical(param_file, "CORRECT_SURFACE_LAYER_AVERAGE")
  call obsolete_int(param_file, "SEAMOUNT_LENGTH_SCALE", hint="Use SEAMOUNT_X_LENGTH_SCALE instead.")
  call obsolete_int(param_file, "USE_LATERAL_BOUNDARY_DIFFUSION", &
                    hint="Use USE_HORIZONTAL_BOUNDARY_DIFFUSION instead.")

  call obsolete_logical(param_file, "MSTAR_FIXED", hint="Instead use MSTAR_MODE.")
  call obsolete_logical(param_file, "USE_VISBECK_SLOPE_BUG", .false.)
  call obsolete_logical(param_file, "Use_PP81", hint="get_param is case sensitive so use USE_PP81.")

  call obsolete_logical(param_file, "ALLOW_CLOCKS_IN_OMP_LOOPS", .true.)
  call obsolete_logical(param_file, "LARGE_FILE_SUPPORT", .true.)
  call obsolete_real(param_file, "MIN_Z_DIAG_INTERVAL")
  call obsolete_char(param_file, "Z_OUTPUT_GRID_FILE")

  call obsolete_logical(param_file, "CFL_BASED_TRUNCATIONS", .true.)
  call obsolete_logical(param_file, "KD_BACKGROUND_VIA_KDML_BUG", .false.)
  call obsolete_logical(param_file, "USE_DIABATIC_TIME_BUG", .false.)

  call read_param(param_file, "INTERPOLATE_SPONGE_TIME_SPACE", test_logic)
  call obsolete_logical(param_file, "NEW_SPONGES", warning_val=test_logic, &
                        hint="Use INTERPOLATE_SPONGE_TIME_SPACE instead.")

  test_logic = .true. ; call read_param(param_file, "BOUND_KH", test_logic)
  call obsolete_logical(param_file, "BETTER_BOUND_KH", warning_val=test_logic, hint="Use BOUND_KH alone.")
  test_logic = .true. ; call read_param(param_file, "BOUND_AH", test_logic)
  call obsolete_logical(param_file, "BETTER_BOUND_AH", warning_val=test_logic, hint="Use BOUND_AH alone.")

  test_logic = .false. ; call read_param(param_file, "UNSPLIT_DT_VISC_BUG", test_logic)
  call obsolete_logical(param_file, "FIX_UNSPLIT_DT_VISC_BUG", warning_val=(.not.test_logic), &
                        hint="Use UNSPLIT_DT_VISC_BUG instead, but with the reversed meaning.")

  call obsolete_logical(param_file, "SMOOTH_RI", hint="Instead use N_SMOOTH_RI.")

  call obsolete_logical(param_file, "INTERNAL_TIDE_CORNER_ADVECT", .false.)
  call obsolete_logical(param_file, "TIDE_USE_SAL_SCALAR", hint="Use SAL_SCALAR_APPROX instead.")
  call obsolete_logical(param_file, "TIDAL_SAL_SHT", hint="Use SAL_HARMONICS instead.")
  call obsolete_int(param_file, "TIDAL_SAL_SHT_DEGREE", hint="Use SAL_HARMONICS_DEGREE instead.")
  call obsolete_real(param_file, "RHO_E", hint="Use RHO_SOLID_EARTH instead.")
  call obsolete_logical(param_file, "DEFAULT_2018_ANSWERS", hint="Instead use DEFAULT_ANSWER_DATE.")

  call obsolete_logical(param_file, "SURFACE_FORCING_2018_ANSWERS", &
                        hint="Instead use SURFACE_FORCING_ANSWER_DATE.")
  call obsolete_logical(param_file, "WIND_GYRES_2018_ANSWERS", &
                        hint="Instead use WIND_GYRES_ANSWER_DATE.")

  call obsolete_logical(param_file, "BAROTROPIC_2018_ANSWERS", &
                        hint="Instead use BAROTROPIC_ANSWER_DATE.")
  call obsolete_logical(param_file, "EPBL_2018_ANSWERS", hint="Instead use EPBL_ANSWER_DATE.")
  call obsolete_logical(param_file, "HOR_REGRID_2018_ANSWERS", &
                        hint="Instead use HOR_REGRID_ANSWER_DATE.")
  call obsolete_logical(param_file, "HOR_VISC_2018_ANSWERS", &
                        hint="Instead use HOR_VISC_ANSWER_DATE.")
  call obsolete_logical(param_file, "IDL_HURR_2018_ANSWERS", &
                        hint="Instead use IDL_HURR_ANSWER_DATE.")
  call obsolete_logical(param_file, "MEKE_GEOMETRIC_2018_ANSWERS", &
                        hint="Instead use MEKE_GEOMETRIC_ANSWER_DATE.")
  call obsolete_logical(param_file, "ODA_2018_ANSWERS", hint="Instead use ODA_ANSWER_DATE.")
  call obsolete_logical(param_file, "OPTICS_2018_ANSWERS", hint="Instead use OPTICS_ANSWER_DATE.")
  call obsolete_logical(param_file, "REGULARIZE_LAYERS_2018_ANSWERS", &
                        hint="Instead use REGULARIZE_LAYERS_ANSWER_DATE.")
  call obsolete_logical(param_file, "REMAPPING_2018_ANSWERS", &
                        hint="Instead use REMAPPING_ANSWER_DATE.")
  call obsolete_logical(param_file, "SET_DIFF_2018_ANSWERS", &
                        hint="Instead use SET_DIFF_ANSWER_DATE.")
  call obsolete_logical(param_file, "SET_VISC_2018_ANSWERS", &
                        hint="Instead use SET_VISC_ANSWER_DATE.")
  call obsolete_logical(param_file, "SURFACE_2018_ANSWERS", hint="Instead use SURFACE_ANSWER_DATE.")
  call obsolete_logical(param_file, "TIDAL_MIXING_2018_ANSWERS", &
                        hint="Instead use TIDAL_MIXING_ANSWER_DATE.")
  call obsolete_logical(param_file, "VERT_FRICTION_2018_ANSWERS", &
                        hint="Instead use VERT_FRICTION_ANSWER_DATE.")

  call obsolete_logical(param_file, "USE_GRID_SPACE_DIAGNOSTIC_AXES", &
                        hint="Instead use USE_INDEX_DIAGNOSTIC_AXIS.")

  ! Write the file version number to the model log.
  call log_version(param_file, mdl, version)

end subroutine find_obsolete_params

!> Test for presence of obsolete LOGICAL in parameter file.
subroutine obsolete_logical(param_file, varname, warning_val, hint)
  type(param_file_type), intent(in) :: param_file  !< Structure containing parameter file data.
  character(len=*),      intent(in) :: varname     !< Name of obsolete LOGICAL parameter.
  logical,     optional, intent(in) :: warning_val !< An allowed value that causes a warning instead of an error.
  character(len=*), optional, intent(in) :: hint   !< A hint to the user about what to do.
  ! Local variables
  logical :: test_logic, fatal_err
  logical :: var_is_set  ! True if this value was read by read_param.
  character(len=128) :: hint_msg

  test_logic = .false. ; call read_param(param_file, varname, test_logic, set=var_is_set)
  fatal_err = .true.
  if (var_is_set .and. present(warning_val)) fatal_err = (warning_val .neqv. test_logic)
  hint_msg = " " ; if (present(hint)) hint_msg = hint

  if (var_is_set) then
    if (fatal_err) then
      call MOM_ERROR(FATAL, "MOM_obsolete_params: "//trim(varname)//   &
           " is an obsolete run-time flag, and should not be used. "// &
           trim(hint_msg))
    else
      call MOM_ERROR(WARNING, "MOM_obsolete_params: "//trim(varname)// &
           " is an obsolete run-time flag. "//trim(hint_msg))
    endif
  endif

end subroutine obsolete_logical

!> Test for presence of obsolete STRING in parameter file.
subroutine obsolete_char(param_file, varname, warning_val, hint)
  type(param_file_type), intent(in) :: param_file !< Structure containing parameter file data.
  character(len=*),      intent(in) :: varname    !< Name of obsolete STRING parameter.
  character(len=*), optional, intent(in) :: warning_val !< An allowed value that causes a warning instead of an error.
  character(len=*), optional, intent(in) :: hint  !< A hint to the user about what to do.
  ! Local variables
  character(len=200) :: test_string, hint_msg
  logical :: var_is_set  ! True if this value was read by read_param.
  logical :: only_warn

  test_string = ''; call read_param(param_file, varname, test_string, set=var_is_set)
  hint_msg = " " ; if (present(hint)) hint_msg = hint

  if (var_is_set) then
    only_warn = .false.
    if (present(warning_val)) then ! Check if test_string and warning_val are the same.
      if (len_trim(warning_val) == len_trim(test_string)) then
        if (index(trim(test_string), trim(warning_val)) == 1) only_warn = .true.
      endif
    endif

    if (only_warn) then
      call MOM_ERROR(WARNING, &
             "MOM_obsolete_params: "//trim(varname)// &
             " is an obsolete run-time flag. "//trim(hint_msg))
    else
      call MOM_ERROR(FATAL, &
             "MOM_obsolete_params: "//trim(varname)// &
             " is an obsolete run-time flag, and should not be used. "//trim(hint_msg))
    endif
  endif
end subroutine obsolete_char

!> Test for presence of obsolete REAL in parameter file.
subroutine obsolete_real(param_file, varname, warning_val, hint, only_warn)
  type(param_file_type), intent(in) :: param_file  !< Structure containing parameter file data.
  character(len=*),      intent(in) :: varname     !< Name of obsolete REAL parameter.
  real,        optional, intent(in) :: warning_val !< An allowed value that causes a warning instead of an error.
  character(len=*), optional, intent(in) :: hint   !< A hint to the user about what to do.
  logical,     optional, intent(in) :: only_warn   !< If present and true, issue warnings instead of fatal errors.

  ! Local variables
  real :: test_val, warn_val
  logical :: var_is_set  ! True if this value was read by read_param.
  logical :: issue_warning
  character(len=128) :: hint_msg

  test_val = -9e35; call read_param(param_file, varname, test_val, set=var_is_set)
  warn_val = -9e35; if (present(warning_val)) warn_val = warning_val
  hint_msg = " " ; if (present(hint)) hint_msg = hint
  issue_warning = .false. ; if (present(only_warn)) issue_warning = only_warn

  if (var_is_set) then
    if ((test_val == warn_val) .or. issue_warning) then
      call MOM_ERROR(WARNING, "MOM_obsolete_params: "//trim(varname)// &
         " is an obsolete run-time flag. "//trim(hint_msg))
    else
      call MOM_ERROR(FATAL, "MOM_obsolete_params: "//trim(varname)//   &
           " is an obsolete run-time flag, and should not be used. "// &
           trim(hint_msg))
    endif
  endif
end subroutine obsolete_real

!> Test for presence of obsolete INTEGER in parameter file.
subroutine obsolete_int(param_file, varname, warning_val, hint)
  type(param_file_type), intent(in) :: param_file  !< Structure containing parameter file data.
  character(len=*),      intent(in) :: varname     !< Name of obsolete INTEGER parameter.
  integer,     optional, intent(in) :: warning_val !< An allowed value that causes a warning instead of an error.
  character(len=*), optional, intent(in) :: hint   !< A hint to the user about what to do.
  ! Local variables
  logical :: var_is_set  ! True if this value was read by read_param.
  integer :: test_val, warn_val
  character(len=128) :: hint_msg

  test_val = -123456788; call read_param(param_file, varname, test_val, set=var_is_set)
  warn_val = -123456788; if (present(warning_val)) warn_val = warning_val
  hint_msg = " " ; if (present(hint)) hint_msg = hint

  if (var_is_set) then
    if (test_val == warn_val) then
      call MOM_ERROR(WARNING, "MOM_obsolete_params: "//trim(varname)// &
         " is an obsolete run-time flag. "//trim(hint_msg))
    else
      call MOM_ERROR(FATAL, "MOM_obsolete_params: "//trim(varname)//   &
           " is an obsolete run-time flag, and should not be used. "// &
           trim(hint_msg))
    endif
  endif
end subroutine obsolete_int

end module MOM_obsolete_params
