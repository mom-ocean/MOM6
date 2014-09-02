module MOM_obsolete_params
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
!
!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, July 2010                                      *
!*                                                                     *
!*    The subroutine in this module checks whether obsolete parameters *
!*  are being used, and stops the model or issues a warning as         *
!*  appropriate.                                                       *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : read_param, log_version, param_file_type

implicit none ; private

#include <MOM_memory.h>

public find_obsolete_params

contains

subroutine find_obsolete_params(param_file)
  type(param_file_type),       intent(in)    :: param_file
! Argument: param_file - A structure indicating the open file to parse for
!                        model parameter values.
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "find_obsolete_params"  ! This module's name.
  integer :: test_int
  logical :: test_logic, test_logic2, test_logic3, split

  if (.not.is_root_pe()) return

  test_int = -1 ; call read_param(param_file,"NTSTEP",test_int)
  if (test_int /= -1) call MOM_ERROR(FATAL, "find_obsolete_params: "// &
      "NTSTEP is an obsolete option.  Instead #define DT_THERM \n" // &
      "  to set the thermodynamic time step.")

  test_logic = .false. ; call read_param(param_file,"JACOBIAN_PGF",test_logic)
  if (test_logic) call MOM_ERROR(FATAL, "find_obsolete_params: "// &
         "JACOBIAN_PGF is an obsolete run-time flag.  Instead use \n"// &
         "  #define ANALYTIC_FV_PGF.")
  test_logic = .true. ; call read_param(param_file,"JACOBIAN_PGF",test_logic)
  if (.not.test_logic) call MOM_ERROR(WARNING, "find_obsolete_params: "// &
         "JACOBIAN_PGF is an obsolete run-time flag.")

  test_logic = .true. ; call read_param(param_file,"SADOURNY",test_logic)
  test_logic2 = .false. ; call read_param(param_file,"SADOURNY",test_logic2)
  if (test_logic .eqv. test_logic2) call MOM_ERROR(FATAL, "find_obsolete_params: "// &
         "SADOURNY is an obsolete run-time flag.  Use CORIOLIS_SCHEME instead.")

  test_logic = .true. ; call read_param(param_file,"ARITHMETIC_BT_THICK",test_logic)
  test_logic2 = .false. ; call read_param(param_file,"ARITHMETIC_BT_THICK",test_logic2)
  if (test_logic .eqv. test_logic2) call MOM_ERROR(FATAL, "find_obsolete_params: "// &
         "ARITHMETIC_BT_THICK is an obsolete run-time flag.  Use \n"//&
         "  #define BT_THICK_SCHEME ARITHMETIC instead.")

  test_logic = .true. ; call read_param(param_file,"HYBRID_BT_THICK",test_logic)
  test_logic2 = .false. ; call read_param(param_file,"HYBRID_BT_THICK",test_logic2)
  if (test_logic .eqv. test_logic2) call MOM_ERROR(WARNING, "find_obsolete_params: "// &
         "  HYBRID_BT_THICK is an obsolete run-time flag.  Use \n"//&
         "#define BT_THICK_SCHEME HYBRID instead.")

  test_logic = .true. ; call read_param(param_file,"BT_CONT_BT_THICK",test_logic)
  test_logic2 = .false. ; call read_param(param_file,"BT_CONT_BT_THICK",test_logic2)
  if (test_logic .eqv. test_logic2) call MOM_ERROR(WARNING, "find_obsolete_params: "// &
         "BT_CONT_BT_THICK is an obsolete run-time flag.  Use \n"//&
         "  #define BT_THICK_SCHEME FROM_BT_CONT instead.")

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

  call test_obsolete_int(param_file, "NXTOT")
  call test_obsolete_int(param_file, "NYTOT")
  call test_obsolete_int(param_file, "NZ")
  call test_obsolete_int(param_file, "NXPROC")
  call test_obsolete_int(param_file, "NYPROC")
  call test_obsolete_int(param_file, "NXPROC_IO")
  call test_obsolete_int(param_file, "NYPROC_IO")
  call test_obsolete_int(param_file, "NXHALO")
  call test_obsolete_int(param_file, "NYHALO")
  call test_obsolete_int(param_file, "ML_PRESORT_NZ_CONV_ADJ")

  call test_obsolete_real(param_file, "BT_COR_SLOW_RATE", 0.0)
  call test_obsolete_real(param_file, "BT_COR_FRAC", 1.0)

  call test_obsolete_logical(param_file, "BT_INCLUDE_UDHDT", .false.)

  call test_obsolete_logical(param_file, "RIGA_SET_DIFFUSIVITY", .false.)
  call test_obsolete_logical(param_file, "RIGA_ITIDE_BUGS", .false.)
  call test_obsolete_logical(param_file, "RIGA_ENTRAINMENT_FOIBLES", .false.)
  call test_obsolete_logical(param_file, "RIGA_TRACER_DIFFUSE_BUGS", .false.)
  call test_obsolete_logical(param_file, "RIGA_KAPPA_SHEAR_BUGS1", .false.)
  call test_obsolete_logical(param_file, "RIGA_KAPPA_SHEAR_BUGS2", .false.)
  call test_obsolete_logical(param_file, "CONT_PPM_RIGA_BUGS", .false.)
  call test_obsolete_logical(param_file, "USE_REPRODUCING_SUM", .true.)
  call test_obsolete_logical(param_file, "SLOW_BITWISE_GLOBAL_FORCING_SUMS", .false.)
  call test_obsolete_logical(param_file, "ALWAYS_WRITE_GEOM")
  call test_obsolete_real(param_file, "I_ZETA")

  call test_obsolete_logical(param_file, "REF_COMPRESS_3D")
  call test_obsolete_char(param_file, "COMPRESS_FILE")
  call test_obsolete_char(param_file, "REF_COMPRESS_FILE_TEMP")
  call test_obsolete_char(param_file, "REF_COMPRESS_FILE_SALT")
  call test_obsolete_char(param_file, "REF_COMPRESS_FILE_DEPTH")

  call test_obsolete_logical(param_file, "OLD_RESTRAT_PARAM", .false.)
  call test_obsolete_real(param_file, "ML_RESTRAT_COEF", 0.0)
  call test_obsolete_logical(param_file, "FULL_THICKNESSDIFFUSE", .true.)
  call test_obsolete_logical(param_file, "DIFFUSE_ISOPYCNALS", .true.)

  call test_obsolete_logical(param_file, "MOREL_PEN_SW")
  call test_obsolete_logical(param_file, "MANIZZA_PEN_SW")

  call test_obsolete_logical(param_file, "USE_H2000_SHEAR_MIXING", .false.)
  call test_obsolete_real(param_file, "SHEARMIX_LAT_EQ", 0.0)
  call test_obsolete_real(param_file, "RINO_CRIT_EQ")
  call test_obsolete_real(param_file, "SHEARMIX_RATE_EQ")

  call test_obsolete_logical(param_file, "CONTINUITY_PPM", .true.)

  call test_obsolete_logical(param_file, "USE_LOCAL_PREF", .true.)
  call test_obsolete_logical(param_file, "USE_LOCAL_PREF_CORRECT", .true.)
  test_logic = .false. ; call read_param(param_file, "USE_JACKSON_PARAM", test_logic)
  call test_obsolete_logical(param_file, "RINOMIX", test_logic)
  call test_obsolete_logical(param_file, "NORMALIZED_SUM_OUT", .true.)

  call test_obsolete_real(param_file, "RLAY_RANGE")
  call test_obsolete_real(param_file, "RLAY_REF")

  call test_obsolete_real(param_file, "HMIX")

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
    call test_obsolete_logical(param_file, "FLUX_BT_COUPLING", .false.)
    call test_obsolete_logical(param_file, "READJUST_BT_TRANS", .false.)
    call test_obsolete_logical(param_file, "RESCALE_BE_FACE_AREAS", .false.)
    call test_obsolete_logical(param_file, "APPLY_BT_DRAG", .true.)
  endif

  ! Write the file version number to the model log.
  call log_version(param_file, mod, version)

end subroutine find_obsolete_params

subroutine test_obsolete_logical(param_file, varname, warning_val)
  type(param_file_type), intent(in) :: param_file
  character(len=*),      intent(in) :: varname
  logical,     optional, intent(in) :: warning_val

  logical :: test_logic, fatal_err

  test_logic = .false. ; call read_param(param_file, varname,test_logic)
  fatal_err = .true.
  if (present(warning_val)) fatal_err = (warning_val .neqv. .true.)

  if (test_logic) then
    if (fatal_err) then
      call MOM_ERROR(FATAL, "find_obsolete_params: "//trim(varname)// &
           " is an obsolete run-time flag, and should not be used.")
    else
      call MOM_ERROR(WARNING, "find_obsolete_params: "//trim(varname)// &
         " is an obsolete run-time flag.")
    endif
  endif

  test_logic = .true. ; call read_param(param_file, varname, test_logic)
  fatal_err = .true.
  if (present(warning_val)) fatal_err = (warning_val .neqv. .false.)

  if (.not.test_logic) then
    if (fatal_err) then
      call MOM_ERROR(FATAL, "find_obsolete_params: "//trim(varname)// &
           " is an obsolete run-time flag, and should not be used.")
    else
      call MOM_ERROR(WARNING, "find_obsolete_params: "//trim(varname)// &
         " is an obsolete run-time flag.")
    endif
  endif
  
end subroutine test_obsolete_logical

subroutine test_obsolete_char(param_file, varname)
  type(param_file_type), intent(in) :: param_file
  character(len=*),      intent(in) :: varname
  character(len=200) :: test_string
  
  test_string = ""; call read_param(param_file, varname, test_string)
  
  if (len_trim(test_string) > 0) call MOM_ERROR(FATAL, &
           "find_obsolete_params: "//trim(varname)// &
           " is an obsolete run-time flag, and should not be used.")

end subroutine test_obsolete_char

subroutine test_obsolete_real(param_file, varname, warning_val)
  type(param_file_type), intent(in) :: param_file
  character(len=*),      intent(in) :: varname
  real,        optional, intent(in) :: warning_val

  real :: test_val, warn_val
  
  test_val = -9e35; call read_param(param_file, varname, test_val)
  warn_val = -9e35; if (present(warning_val)) warn_val = warning_val
  
  if (test_val /= -9e35) then
    if (test_val == warn_val) then
      call MOM_ERROR(WARNING, "find_obsolete_params: "//trim(varname)// &
         " is an obsolete run-time flag.")
    else
      call MOM_ERROR(FATAL, "find_obsolete_params: "//trim(varname)// &
           " is an obsolete run-time flag, and should not be used.")
    endif
  endif
end subroutine test_obsolete_real

subroutine test_obsolete_int(param_file, varname, warning_val)
  type(param_file_type), intent(in) :: param_file
  character(len=*),      intent(in) :: varname
  integer,     optional, intent(in) :: warning_val

  integer :: test_val, warn_val
  
  test_val = -123456788; call read_param(param_file, varname, test_val)
  warn_val = -123456788; if (present(warning_val)) warn_val = warning_val
  
  if (test_val /= -123456788) then
    if (test_val == warn_val) then
      call MOM_ERROR(WARNING, "find_obsolete_params: "//trim(varname)// &
         " is an obsolete run-time flag.")
    else
      call MOM_ERROR(FATAL, "find_obsolete_params: "//trim(varname)// &
           " is an obsolete run-time flag, and should not be used.")
    endif
  endif
end subroutine test_obsolete_int

end module MOM_obsolete_params
