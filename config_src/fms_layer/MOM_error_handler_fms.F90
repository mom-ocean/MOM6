!> FMS-based implementation of MOM_error_handler subprograms
submodule (MOM_error_handler) MOM_error_handler_fms

  ! NOTE: These are currently imported from the MOM_error_handler module.
  !   If uncommented, this line would cause a namespace conflict.
  !   Once all references to mpp_error and NOTE are removed from
  !   MOM_error_module, this line should be re-enabled.
  !
  !use mpp_mod, only : mpp_error, NOTE

  implicit none
contains
  !> This provides a convenient interface for writing an informative comment.
  module subroutine MOM_mesg(message, verb, all_print)
  !module subroutine MOM_mesg(message, verb, all_print)
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
    if (write_msg .and. (verbosity >= verb_msg)) call mpp_error(NOTE, message)
  end subroutine MOM_mesg
end submodule MOM_error_handler_fms
