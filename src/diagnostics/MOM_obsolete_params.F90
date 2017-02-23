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
  character(len=40)  :: mod = "find_obsolete_params" ! This module's name.
! This include declares and sets the variable "version".
#include "version_variable.h"
  integer :: test_int
  logical :: test_logic, test_logic2, test_logic3, split

  if (.not.is_root_pe()) return

  call obsolete_int(param_file, "NTSTEP", &
       hint="Instead use DT_THERM to set the thermodynamic time-step.")

  call obsolete_logical(param_file, "JACOBIAN_PGF", .false., &
       hint="Instead use ANALYTIC_FV_PGF.")

  call obsolete_logical(param_file, "SADOURNY", &
       hint="Instead use CORIOLIS_SCHEME='SADOURNY'.")

  call obsolete_logical(param_file, "ARITHMETIC_BT_THICK", &
       hint="Instead use BT_THICK_SCHEME='ARITHMETIC'.")

  call obsolete_logical(param_file, "HYBRID_BT_THICK", &
       hint="Instead use BT_THICK_SCHEME='HYBRID'.")

  call obsolete_logical(param_file, "BT_CONT_BT_THICK", &
       hint="Instead use BT_THICK_SCHEME='FROM_BT_CONT'.")

  call obsolete_logical(param_file, "APPLY_OBC_U", &
       hint="Instead use OBC_NUMBER_SEGMENTS>0 and use the new segments protocol.")
  call obsolete_logical(param_file, "APPLY_OBC_V", &
       hint="Instead use OBC_NUMBER_SEGMENTS>0 and use the new segments protocol.")
  call obsolete_logical(param_file, "APPLY_OBC_V_FLATHER_NORTH", &
       hint="Instead use OBC_NUMBER_SEGMENTS>0 and use the new segments protocol.")
  call obsolete_logical(param_file, "APPLY_OBC_V_FLATHER_SOUTH", &
       hint="Instead use OBC_NUMBER_SEGMENTS>0 and use the new segments protocol.")
  call obsolete_logical(param_file, "APPLY_OBC_U_FLATHER_EAST", &
       hint="Instead use OBC_NUMBER_SEGMENTS>0 and use the new segments protocol.")
  call obsolete_logical(param_file, "APPLY_OBC_U_FLATHER_WEST", &
       hint="Instead use OBC_NUMBER_SEGMENTS>0 and use the new segments protocol.")
  call obsolete_char(param_file, "OBC_CONFIG", &
       hint="Instead use OBC_USER_CONFIG and use the new segments protocol.")

  test_logic3 = .true. ; call read_param(param_file,"ENABLE_THERMODYNAMICS",test_logic3)
  test_logic = .true. ; call read_param(param_file,"TEMPERATURE",test_logic)
  test_logic2 = .false. ; call read_param(param_file,"TEMPERATURE",test_logic2)
  if (test_logic .eqv. test_logic2) then ; if (test_logic .eqv. test_logic3) then
    call MOM_ERROR(WARNING, "find_obsolete_params: "// &
         "TEMPERATURE is an obsolete run-time flag, but is set consistently with \n"//&
         "  ENABLE_THERMODYNAMICS.")
  else
    call MOM_ERROR(FATAL, "find_obsolete_params: "// &
         "TEMPERATURE is an obsolete run-time flag.  Use ENABLE_THERMODYNAMICS instead.")
  endif ; endif

  test_logic = test_logic3 ; call read_param(param_file,"NONLINEAR_EOS",test_logic)
  if (test_logic .neqv. test_logic3) then
    call MOM_error(WARNING, "find_obsolete_params: "// &
          "NONLINEAR_EOS is an obsolete option.  Instead define " // &
          "USE_EOS to use an equation of state to calculate density.")
  endif

! test_logic = .true. ; call read_param(param_file,"USE_RIVER_HEAT_CONTENT",test_logic)
! test_logic2 = .false. ; call read_param(param_file,"USE_RIVER_HEAT_CONTENT",test_logic2)
! if (test_logic .eqv. test_logic2) call MOM_ERROR(FATAL, "find_obsolete_params: "// &
!        "USE_RIVER_HEAT_CONTENT, is an obsolete run-time flag.")

! test_logic = .true. ; call read_param(param_file,"USE_CALVING_HEAT_CONTENT",test_logic)
! test_logic2 = .false. ; call read_param(param_file,"USE_CALVING_HEAT_CONTENT",test_logic2)
! if (test_logic .eqv. test_logic2) call MOM_ERROR(FATAL, "find_obsolete_params: "// &
!        "USE_CALVING_HEAT_CONTENT, is an obsolete run-time flag.")

  call obsolete_int(param_file, "NXTOT")
  call obsolete_int(param_file, "NYTOT")
  call obsolete_int(param_file, "NZ")
  call obsolete_int(param_file, "NXPROC")
  call obsolete_int(param_file, "NYPROC")
  call obsolete_int(param_file, "NXPROC_IO")
  call obsolete_int(param_file, "NYPROC_IO")
  call obsolete_int(param_file, "NXHALO")
  call obsolete_int(param_file, "NYHALO")
  call obsolete_int(param_file, "ML_PRESORT_NZ_CONV_ADJ")

  call obsolete_int(param_file, "NIPROC_IO", hint="Use IO_LAYOUT=#,# instead.")
  call obsolete_int(param_file, "NJPROC_IO", hint="Use IO_LAYOUT=#,# instead.")

  call obsolete_real(param_file, "BT_COR_SLOW_RATE", 0.0)
  call obsolete_real(param_file, "BT_COR_FRAC", 1.0)

  call obsolete_logical(param_file, "BT_INCLUDE_UDHDT", .false.)

  call obsolete_logical(param_file, "RIGA_SET_DIFFUSIVITY", .false.)
  call obsolete_logical(param_file, "RIGA_ITIDE_BUGS", .false.)
  call obsolete_logical(param_file, "RIGA_ENTRAINMENT_FOIBLES", .false.)
  call obsolete_logical(param_file, "RIGA_TRACER_DIFFUSE_BUGS", .false.)
  call obsolete_logical(param_file, "RIGA_KAPPA_SHEAR_BUGS1", .false.)
  call obsolete_logical(param_file, "RIGA_KAPPA_SHEAR_BUGS2", .false.)
  call obsolete_logical(param_file, "CONT_PPM_RIGA_BUGS", .false.)
  call obsolete_logical(param_file, "USE_REPRODUCING_SUM", .true.)
  call obsolete_logical(param_file, "SLOW_BITWISE_GLOBAL_FORCING_SUMS", .false.)
  call obsolete_logical(param_file, "ALWAYS_WRITE_GEOM")
  call obsolete_real(param_file, "I_ZETA")

  call obsolete_logical(param_file, "REF_COMPRESS_3D")
  call obsolete_char(param_file, "COMPRESS_FILE")
  call obsolete_char(param_file, "REF_COMPRESS_FILE_TEMP")
  call obsolete_char(param_file, "REF_COMPRESS_FILE_SALT")
  call obsolete_char(param_file, "REF_COMPRESS_FILE_DEPTH")
  call obsolete_char(param_file, "DIAG_REMAP_Z_GRID_DEF", "Use NUM_DIAG_COORDS, DIAG_COORDS and DIAG_COORD_DEF_Z")

  call obsolete_logical(param_file, "OLD_RESTRAT_PARAM", .false.)
  call obsolete_real(param_file, "ML_RESTRAT_COEF", 0.0)
  call obsolete_logical(param_file, "FULL_THICKNESSDIFFUSE", .true.)
  call obsolete_logical(param_file, "DIFFUSE_ISOPYCNALS", .true.)

  call obsolete_logical(param_file, "MOREL_PEN_SW")
  call obsolete_logical(param_file, "MANIZZA_PEN_SW")

  call obsolete_logical(param_file, "USE_H2000_SHEAR_MIXING", .false.)
  call obsolete_real(param_file, "SHEARMIX_LAT_EQ", 0.0)
  call obsolete_real(param_file, "RINO_CRIT_EQ")
  call obsolete_real(param_file, "SHEARMIX_RATE_EQ")

  call obsolete_logical(param_file, "CONTINUITY_PPM", .true.)

  call obsolete_logical(param_file, "USE_LOCAL_PREF", .true.)
  call obsolete_logical(param_file, "USE_LOCAL_PREF_CORRECT", .true.)
  test_logic = .false. ; call read_param(param_file, "USE_JACKSON_PARAM", test_logic)
  call obsolete_logical(param_file, "RINOMIX", test_logic)
  call obsolete_logical(param_file, "NORMALIZED_SUM_OUT", .true.)

  call obsolete_real(param_file, "RLAY_RANGE")
  call obsolete_real(param_file, "RLAY_REF")

  call obsolete_real(param_file, "HMIX")
  call obsolete_real(param_file, "VSTAR_SCALE_COEF")
  call obsolete_real(param_file, "ZSTAR_RIGID_SURFACE_THRESHOLD")

  test_int = -1 ; call read_param(param_file,"ML_RADIATION_CODING",test_int)
  if (test_int == 1) call MOM_ERROR(FATAL, "find_obsolete_params: "// &
    "ML_RADIATION_CODING is an obsolete option and the code previously "//&
    "used by setting it to 1 has been eliminated.")
  if (test_int /= -1) call MOM_ERROR(WARNING, "find_obsolete_params: "// &
    "ML_RADIATION_CODING is an obsolete option.")

  ! Test for inconsistent parameter settings.
  split = .true. ; test_logic = .false.
  call read_param(param_file,"SPLIT",split)
  call read_param(param_file,"DYNAMIC_SURFACE_PRESSURE",test_logic)
  if (test_logic .and. .not.split) call MOM_ERROR(FATAL, &
    "find_obsolete_params: #define DYNAMIC_SURFACE_PRESSURE is not yet "//&
    "implemented without #define SPLIT.")

  call read_param(param_file,"USE_LEGACY_SPLIT",test_logic)
  if (.not.(split .and. test_logic)) then
    call obsolete_logical(param_file, "FLUX_BT_COUPLING", .false.)
    call obsolete_logical(param_file, "READJUST_BT_TRANS", .false.)
    call obsolete_logical(param_file, "RESCALE_BE_FACE_AREAS", .false.)
    call obsolete_logical(param_file, "APPLY_BT_DRAG", .true.)
  endif

  call obsolete_int(param_file, "SEAMOUNT_LENGTH_SCALE", hint="Use SEAMOUNT_X_LENGTH_SCALE instead.")

  ! Write the file version number to the model log.
  call log_version(param_file, mod, version)

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
subroutine obsolete_char(param_file, varname, hint)
  type(param_file_type), intent(in) :: param_file !< Structure containing parameter file data.
  character(len=*),      intent(in) :: varname    !< Name of obsolete STRING parameter.
  character(len=*), optional, intent(in) :: hint  !< A hint to the user about what to do.
  ! Local variables
  character(len=200) :: test_string, hint_msg

  test_string = ''; call read_param(param_file, varname, test_string)
  hint_msg = " " ; if (present(hint)) hint_msg = hint

  if (len_trim(test_string) > 0) call MOM_ERROR(FATAL,                 &
           "MOM_obsolete_params: "//trim(varname)//                    &
           " is an obsolete run-time flag, and should not be used. "// &
           trim(hint_msg))

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
