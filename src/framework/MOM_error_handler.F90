!> Routines for error handling and I/O management
module MOM_error_handler

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_infra, only : MOM_err, is_root_pe, stdlog, stdout, NOTE, WARNING, FATAL

implicit none ; private

! These routines are found in this module.
public :: MOM_error, MOM_mesg, assert
public :: MOM_set_verbosity, MOM_get_verbosity, MOM_verbose_enough
public :: callTree_showQuery, callTree_enter, callTree_leave, callTree_waypoint
! These routines are simply passed-through from MOM_error_infra
public :: is_root_pe, stdlog, stdout
!> Integer parameters encoding the severity of an error message
public :: NOTE, WARNING, FATAL

integer :: verbosity = 6
!< Verbosity level:
!!  0 - FATAL messages only
!!  1 - FATAL + WARNING messages only
!!  2 - FATAL + WARNING + NOTE messages only [default]
!!  3 - above + informational
!!  4 -
!!  5 -
!!  6 - above + call tree
!!  7 -
!!  8 -
!!  9 - anything and everything (also set with DEBUG=True)

!   Note that this module default will only hold until the
! VERBOSITY parameter is parsed and the given default imposed.
! We set it to 6 here so that the call tree will print before
! the parser has been initialized
!   Also note that this is a module variable rather than contained in
! a type passed by argument (preferred for most data) for convenience
! and to reduce obfuscation of code

integer :: callTreeIndentLevel = 0
!< The level of calling within the call tree

contains

!> This provides a convenient interface for writing an informative comment, depending
!! on the model's current verbosity setting and the verbosity level for this message.
subroutine MOM_mesg(message, verb, all_print)
  character(len=*), intent(in)  :: message !< A message to write out
  integer, optional, intent(in) :: verb !< A level of verbosity for this message
  logical, optional, intent(in) :: all_print !< If present and true, any PEs are
                                             !! able to write this message.
  ! This provides a convenient interface for writing an informative comment.
  integer :: verb_msg
  logical :: write_msg

  write_msg = is_root_pe()
  if (present(all_print)) write_msg = write_msg .or. all_print

  verb_msg = 2 ; if (present(verb)) verb_msg = verb
  if (write_msg .and. (verbosity >= verb_msg)) call MOM_err(NOTE, message)

end subroutine MOM_mesg

!> This provides a convenient interface for writing an error message
!! with run-time filter based on a verbosity and the severity of the error.
subroutine MOM_error(level, message, all_print)
  integer,           intent(in) :: level !< The severity level of this message
  character(len=*),  intent(in) :: message !< A message to write out
  logical, optional, intent(in) :: all_print !< If present and true, any PEs are
                                             !! able to write this message.
  ! This provides a convenient interface for writing an error message
  ! with run-time filter based on a verbosity.
  logical :: write_msg

  write_msg = is_root_pe()
  if (present(all_print)) write_msg = write_msg .or. all_print

  select case (level)
    case (NOTE)
      if (write_msg.and.verbosity>=2) call MOM_err(NOTE, message)
    case (WARNING)
      if (write_msg.and.verbosity>=1) call MOM_err(WARNING, message)
    case (FATAL)
      if (verbosity>=0) call MOM_err(FATAL, message)
    case default
      call MOM_err(level, message)
  end select
end subroutine MOM_error

!> This subroutine sets the level of verbosity filtering MOM error messages
subroutine MOM_set_verbosity(verb)
  integer, intent(in) :: verb !< A level of verbosity to set
  character(len=80) :: msg
  if (verb>0 .and. verb<10) then
    verbosity=verb
  else
    write(msg(1:80),'("Attempt to set verbosity outside of range (0-9). verb=",I0)') verb
    call MOM_error(FATAL,msg)
  endif
end subroutine MOM_set_verbosity

!> This subroutine gets the level of verbosity filtering MOM error messages
function MOM_get_verbosity()
  integer :: MOM_get_verbosity
  MOM_get_verbosity = verbosity
end function MOM_get_verbosity

!> This tests whether the level of verbosity filtering MOM error messages is
!! sufficient to write a message of verbosity level verb
function MOM_verbose_enough(verb)
  integer, intent(in) :: verb !< A level of verbosity to test
  logical :: MOM_verbose_enough
  MOM_verbose_enough = (verbosity >= verb)
end function MOM_verbose_enough

!> Returns True, if the verbosity>=6 indicating to show the call tree
function callTree_showQuery()
  ! Local variables
  logical :: callTree_showQuery
  callTree_showQuery = (verbosity >= 6)
end function callTree_showQuery

!> Writes a message about entering a subroutine if call tree reporting is active
subroutine callTree_enter(mesg,n)
  character(len=*),  intent(in) :: mesg !< Message to write
  integer, optional, intent(in) :: n !< An optional integer to write at end of message
  ! Local variables
  character(len=8) :: nAsString
  callTreeIndentLevel = callTreeIndentLevel + 1
  if (verbosity<6) return
  if (is_root_pe()) then
    nAsString = ''
    if (present(n)) then
      write(nAsString(1:8),'(i8)') n
      call MOM_err(NOTE, 'callTree: '// &
        repeat('   ',callTreeIndentLevel-1)//'loop '//trim(mesg)//trim(nAsString))
    else
      call MOM_err(NOTE, 'callTree: '// &
        repeat('   ',callTreeIndentLevel-1)//'---> '//trim(mesg))
    endif
  endif
end subroutine callTree_enter

!> Writes a message about leaving a subroutine if call tree reporting is active
subroutine callTree_leave(mesg)
  character(len=*) :: mesg !< Message to write
  if (callTreeIndentLevel<1) write(0,*) 'callTree_leave: error callTreeIndentLevel=',callTreeIndentLevel,trim(mesg)
  callTreeIndentLevel = callTreeIndentLevel - 1
  if (verbosity<6) return
  if (is_root_pe()) call MOM_err(NOTE, 'callTree: '// &
        repeat('   ',callTreeIndentLevel)//'<--- '//trim(mesg))
end subroutine callTree_leave

!> Writes a message about reaching a milestone if call tree reporting is active
subroutine callTree_waypoint(mesg,n)
  character(len=*),  intent(in) :: mesg !< Message to write
  integer, optional, intent(in) :: n !< An optional integer to write at end of message
  ! Local variables
  character(len=8) :: nAsString
  if (callTreeIndentLevel<0) write(0,*) 'callTree_waypoint: error callTreeIndentLevel=',callTreeIndentLevel,trim(mesg)
  if (verbosity<6) return
  if (is_root_pe()) then
    nAsString = ''
    if (present(n)) then
      write(nAsString(1:8),'(i8)') n
      call MOM_err(NOTE, 'callTree: '// &
        repeat('   ',callTreeIndentLevel)//'loop '//trim(mesg)//trim(nAsString))
    else
      call MOM_err(NOTE, 'callTree: '// &
        repeat('   ',callTreeIndentLevel)//'o '//trim(mesg))
    endif
  endif
end subroutine callTree_waypoint

!> Issues a FATAL error if the assertion fails, i.e. the first argument is false.
subroutine assert(logical_arg, msg)
  logical, intent(in) :: logical_arg !< If false causes a FATAL error
  character(len=*), intent(in) :: msg !< Message to issue in case of failed assertion

  if (.not. logical_arg) then
    call MOM_error(FATAL, msg)
  endif

end subroutine assert

end module MOM_error_handler
