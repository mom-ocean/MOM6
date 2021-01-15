!> Routines for error handling and I/O management
module MOM_error_infra

! This file is part of MOM6. See LICENSE.md for the license.

use mpp_mod, only : MOM_err => mpp_error, NOTE, WARNING, FATAL
use mpp_mod, only : mpp_pe, mpp_root_pe, stdlog, stdout

implicit none ; private

public MOM_err, NOTE, WARNING, FATAL, is_root_pe, stdlog, stdout

contains

! MOM_err writes an error message, and may stop the run depending on the
! severity of the error.

!> is_root_pe returns .true. if the current PE is the root PE.
function is_root_pe()
  logical :: is_root_pe
  is_root_pe = .false.
  if (mpp_pe() == mpp_root_pe()) is_root_pe = .true.
end function is_root_pe

end module MOM_error_infra
