!> Routines for error handling and I/O management
module MOM_error_infra

! This file is part of MOM6. See LICENSE.md for the license.

use mpp_mod, only : mpp_error, mpp_pe, mpp_root_pe, mpp_stdlog=>stdlog, mpp_stdout=>stdout
use mpp_mod, only : NOTE, WARNING, FATAL

implicit none ; private

public :: MOM_err, is_root_pe, stdlog, stdout
!> Integer parameters encoding the severity of an error message
public :: NOTE, WARNING, FATAL

contains

!> MOM_err writes an error message, and may cause the run to stop depending on the
!! severity of the error.
subroutine MOM_err(severity, message)
  integer,           intent(in) :: severity !< The severity level of this error
  character(len=*),  intent(in) :: message  !< A message to write out

  call mpp_error(severity, message)
end subroutine MOM_err

!> stdout returns the standard Fortran unit number for output
integer function stdout()
  stdout = mpp_stdout()
end function stdout

!> stdlog returns the standard Fortran unit number to use to log messages
integer function stdlog()
  stdlog = mpp_stdlog()
end function stdlog

!> is_root_pe returns .true. if the current PE is the root PE.
logical function is_root_pe()
  is_root_pe = .false.
  if (mpp_pe() == mpp_root_pe()) is_root_pe = .true.
end function is_root_pe

end module MOM_error_infra
