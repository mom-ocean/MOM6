module MOM_error_handler
!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*    This file is a part of MOM.  See MOM.F90 for licensing.          *
!*                                                                     *
!*  By R. Hallberg, 2005-2012.                                         *
!*                                                                     *
!*    This module wraps the mpp_mod error handling code and the        *
!*  mpp functions stdlog() and stdout() that return open unit numbers. *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use mpp_mod, only : mpp_error, NOTE, WARNING, FATAL
use mpp_mod, only : mpp_pe, mpp_root_pe, stdlog, stdout

implicit none ; private

public MOM_error, MOM_mesg, NOTE, WARNING, FATAL, is_root_pe, stdlog, stdout
public MOM_set_verbosity, MOM_get_verbosity, MOM_verbose_enough

! Verbosity level:
!  0 - FATAL messages only
!  1 - FATAL + WARNING messages only
!  2 - FATAL + WARNING + NOTE messages only [default]
!  3 - above + informational
!  4 -
!  5 -
!  6 - above + call tree
!  7 -
!  8 -
!  9 - anything and everything (also set with #define DEBUG)
integer :: verbosity = 2
! Note that this is a module variable rather than contained in a
! type passed by argument (preferred for most data) for convenience
! and to reduce obfiscation of code

contains

function is_root_pe()
  ! This returns .true. if the current PE is the root PE.
  logical :: is_root_pe
  is_root_pe = .false.
  if (mpp_pe() == mpp_root_pe()) is_root_pe = .true.
  return
end function is_root_pe

subroutine MOM_mesg(message, verb, all_print)
  character(len=*), intent(in)  :: message
  integer, optional, intent(in) :: verb
  logical, optional, intent(in) :: all_print
  ! This provides a convenient interface for writing an informative comment.
  integer :: verb_msg
  logical :: write_msg
  
  write_msg = is_root_pe()
  if (present(all_print)) write_msg = write_msg .or. all_print

  verb_msg = 2 ; if (present(verb)) verb_msg = verb
  if (write_msg .and. (verbosity >= verb_msg)) call mpp_error(NOTE, message)

end subroutine MOM_mesg

subroutine MOM_error(level, message, all_print)
  integer,           intent(in) :: level
  character(len=*),  intent(in) :: message
  logical, optional, intent(in) :: all_print
  ! This provides a convenient interface for writing an mpp_error message
  ! with run-time filter based on a verbosity.
  logical :: write_msg
  
  write_msg = is_root_pe()
  if (present(all_print)) write_msg = write_msg .or. all_print
  
  select case (level)
    case (NOTE)
      if (write_msg.and.verbosity>=2) call mpp_error(NOTE, message)
    case (WARNING)
      if (write_msg.and.verbosity>=1) call mpp_error(WARNING, message)
    case (FATAL)
      if (verbosity>=0) call mpp_error(FATAL, message)
    case default
      call mpp_error(level, message)
  end select
end subroutine MOM_error

subroutine MOM_set_verbosity(verb)
  integer, intent(in) :: verb
  character(len=80) :: msg
  if (verb>0 .and. verb<10) then
    verbosity=verb
  else
    write(msg(1:80),'("Attempt to set verbosity outside of range (0-9). verb=",I0)') verb
    call MOM_error(FATAL,msg)
  endif
end subroutine MOM_set_verbosity

function MOM_get_verbosity()
  integer :: MOM_get_verbosity
  MOM_get_verbosity = verbosity
end function MOM_get_verbosity  

function MOM_verbose_enough(verb)
  integer, intent(in) :: verb
  logical :: MOM_verbose_enough

  MOM_verbose_enough = (verbosity >= verb)
end function MOM_verbose_enough

end module MOM_error_handler
