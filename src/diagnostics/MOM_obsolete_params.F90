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
  integer :: test_int, l_seg, nseg
  logical :: test_logic, test_logic2, test_logic3, split
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

  call obsolete_logical(param_file, "MASK_MASSLESS_TRACERS", .false.)

  call obsolete_logical(param_file, "SALT_REJECT_BELOW_ML", .false.)
  call obsolete_logical(param_file, "MLE_USE_MLD_AVE_BUG", .false.)
  call obsolete_logical(param_file, "KG_BG_2D_BUG", .false.)
  call obsolete_logical(param_file, "CORRECT_DENSITY", .true.)
  call obsolete_char(param_file, "WINDSTRESS_STAGGER", warning_val="C", &
                     hint="Use WIND_STAGGER instead.")

  call obsolete_char(param_file, "DIAG_REMAP_Z_GRID_DEF", &
                     hint="Use NUM_DIAG_COORDS, DIAG_COORDS and DIAG_COORD_DEF_Z")

  call obsolete_real(param_file, "VSTAR_SCALE_FACTOR", hint="Use EPBL_VEL_SCALE_FACTOR instead.")
  call obsolete_logical(param_file, "ORIG_MLD_ITERATION", .false.)

  call obsolete_real(param_file, "VSTAR_SCALE_COEF")
  call obsolete_real(param_file, "ZSTAR_RIGID_SURFACE_THRESHOLD")

  ! Test for inconsistent parameter settings.
  split = .true. ; test_logic = .false.
  call read_param(param_file,"SPLIT",split)
  call read_param(param_file,"DYNAMIC_SURFACE_PRESSURE",test_logic)
  if (test_logic .and. .not.split) call MOM_ERROR(FATAL, &
    "find_obsolete_params: #define DYNAMIC_SURFACE_PRESSURE is not yet "//&
    "implemented without #define SPLIT.")

  call obsolete_real(param_file, "BT_MASS_SOURCE_LIMIT", 0.0)

  call obsolete_int(param_file, "SEAMOUNT_LENGTH_SCALE", hint="Use SEAMOUNT_X_LENGTH_SCALE instead.")

  call obsolete_logical(param_file, "MSTAR_FIXED", hint="Instead use MSTAR_MODE.")
  call obsolete_logical(param_file, "USE_VISBECK_SLOPE_BUG", .false.)

  call obsolete_real(param_file, "MIN_Z_DIAG_INTERVAL")
  call obsolete_char(param_file, "Z_OUTPUT_GRID_FILE")

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
  character(len=128) :: hint_msg

  test_logic = .false. ; call read_param(param_file, varname,test_logic)
  fatal_err = .true.
  if (present(warning_val)) fatal_err = (warning_val .neqv. .true.)
  hint_msg = " " ; if (present(hint)) hint_msg = hint

  if (test_logic) then
    if (fatal_err) then
      call MOM_ERROR(FATAL, "MOM_obsolete_params: "//trim(varname)//   &
           " is an obsolete run-time flag, and should not be used. "// &
           trim(hint_msg))
    else
      call MOM_ERROR(WARNING, "MOM_obsolete_params: "//trim(varname)// &
           " is an obsolete run-time flag. "//trim(hint_msg))
    endif
  endif

  test_logic = .true. ; call read_param(param_file, varname, test_logic)
  fatal_err = .true.
  if (present(warning_val)) fatal_err = (warning_val .neqv. .false.)

  if (.not.test_logic) then
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
  logical :: only_warn

  test_string = ''; call read_param(param_file, varname, test_string)
  hint_msg = " " ; if (present(hint)) hint_msg = hint

  if (len_trim(test_string) > 0) then
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
subroutine obsolete_real(param_file, varname, warning_val, hint)
  type(param_file_type), intent(in) :: param_file  !< Structure containing parameter file data.
  character(len=*),      intent(in) :: varname     !< Name of obsolete REAL parameter.
  real,        optional, intent(in) :: warning_val !< An allowed value that causes a warning instead of an error.
  character(len=*), optional, intent(in) :: hint   !< A hint to the user about what to do.
  ! Local variables
  real :: test_val, warn_val
  character(len=128) :: hint_msg

  test_val = -9e35; call read_param(param_file, varname, test_val)
  warn_val = -9e35; if (present(warning_val)) warn_val = warning_val
  hint_msg = " " ; if (present(hint)) hint_msg = hint

  if (test_val /= -9e35) then
    if (test_val == warn_val) then
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
  integer :: test_val, warn_val
  character(len=128) :: hint_msg

  test_val = -123456788; call read_param(param_file, varname, test_val)
  warn_val = -123456788; if (present(warning_val)) warn_val = warning_val
  hint_msg = " " ; if (present(hint)) hint_msg = hint

  if (test_val /= -123456788) then
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
