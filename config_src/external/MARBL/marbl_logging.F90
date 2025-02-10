!> A non-functioning template of the MARBL logging module
module marbl_logging

    implicit none
    private

    !> A non-functioning template of the marbl status log type
    type, public :: marbl_status_log_entry_type
        integer :: ElementInd  !< dummy index
        logical :: lonly_master_writes  !< dummy flag
        character(len=0) :: LogMessage  !< dummy message
        type(marbl_status_log_entry_type), pointer :: next  !< dummy pointer
    end type marbl_status_log_entry_type

    !> A non-functioning template of the marbl status log type
    type, public :: marbl_log_type
        logical, public  :: labort_marbl  !< dummy flag
        type(marbl_status_log_entry_type), pointer :: FullLog  !< dummy pointer
    contains
        procedure, public :: log_error_trace  !< dummy trace routine
        procedure, public :: erase  !< dummy erase routine
    end type marbl_log_type

contains

    !> dummy trace routine
    subroutine log_error_trace(self, RoutineName, CodeLoc, ElemInd)
        class(marbl_log_type), intent(inout) :: self
        character(len=*),      intent(in)    :: RoutineName, CodeLoc
        integer, optional,     intent(in)    :: ElemInd
    end subroutine log_error_trace

    !> dummy erase routine
    subroutine erase(self)
        class(marbl_log_type), intent(inout) :: self
    end subroutine erase

end module marbl_logging