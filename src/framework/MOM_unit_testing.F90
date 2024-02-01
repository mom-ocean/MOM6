module MOM_unit_testing

use posix, only : chmod
use posix, only : sigsetjmp
use posix, only : sigjmp_buf

use MOM_coms, only : num_PEs, sync_PEs
use MOM_error_handler, only : is_root_pe
use MOM_error_handler, only : disable_fatal_errors
use MOM_error_handler, only : enable_fatal_errors

implicit none ; private

public :: string
public :: create_test_file
public :: delete_test_file
public :: TestSuite


!> String container type
type :: string
  character(len=:), allocatable :: s
    !< Internal character array of string
end type string


!> String constructor
interface string
  module procedure init_string_char
  module procedure init_string_int
end interface string


!> A generalized instance of a unit test function
type :: UnitTest
  private
  procedure(), nopass, pointer :: proc => null()
    !< Unit test function/subroutine
  procedure(), nopass, pointer :: cleanup => null()
    !< Cleanup function to be run after proc
  character(len=:), allocatable :: name
    !< Unit test name (usually set to name of proc)
  logical :: is_fatal
    !< True if proc() is expected to fail
contains
  procedure :: run => run_unit_test
    !< Run the unit test function, proc
end type UnitTest


!> Unit test constructor
interface UnitTest
  module procedure create_unit_test_basic
  module procedure create_unit_test_full
end interface UnitTest


!> Collection of unit tests
type :: TestSuite
  private
  type(UnitTestNode), pointer :: head => null()
    !< Head of the unit test linked list
  type(UnitTestNode), pointer :: tail => null()
    !< Tail of the unit test linked list (pre-allocated and unconfigured)

  ! Public API
  procedure(), nopass, pointer, public :: cleanup => null()
    !< Default cleanup function for unit tests in suite
contains
  private
  procedure :: add_basic => add_unit_test_basic
    !< Add a unit test without a cleanup function
  procedure :: add_full => add_unit_test_full
    !< Add a unit test with an explicit cleanup function
  generic, public :: add => add_basic, add_full
    !< Add a unit test to the test suite
  procedure, public :: run => run_test_suite
    !< Run all unit tests in the suite
end type TestSuite


!> TestSuite constructor
interface TestSuite
  module procedure create_test_suite
end interface TestSuite


!> UnitTest node of TestSuite's linked list
type :: UnitTestNode
  private
  type(UnitTest), pointer :: test => null()
    !< Node contents
  type(UnitTestNode), pointer :: next => null()
    !< Pointer to next node in list
end type UnitTestNode

contains

!> Return a new unit test without a cleanup function
function create_unit_test_basic(proc, name, fatal) result(test)
  procedure() :: proc
    !< Subroutine which defines the unit test
  character(len=*), intent(in) :: name
    !< Name of the unit test
  logical, intent(in), optional :: fatal
    !< True if the test is expected to raise a FATAL error
  type(UnitTest) :: test

  procedure(), pointer :: cleanup
  cleanup => null()

  test = create_unit_test_full(proc, name, fatal, cleanup)
end function create_unit_test_basic


!> Return a new unit test with an explicit cleanup function
function create_unit_test_full(proc, name, fatal, cleanup) result(test)
  procedure() :: proc
    !< Subroutine which defines the unit test
  character(len=*), intent(in) :: name
    !< Name of the unit test
  logical, optional :: fatal
    !< True if the test is expected to raise a FATAL error
  procedure() :: cleanup
    !< Cleanup subroutine, called after test
  type(UnitTest) :: test

  test%proc => proc
  test%name = name
  test%is_fatal = .false.
  if (present(fatal)) test%is_fatal = fatal
  test%cleanup => cleanup
end function create_unit_test_full


!> Launch a unit test with a custom cleanup procedure
subroutine run_unit_test(test)
  class(UnitTest), intent(in) :: test

  type(sigjmp_buf) :: env
  integer :: rc

  call sync_PEs

  ! FIXME: Some FATAL tests under MPI are unable to recover after jumpback, so
  !   we disable these tests for now.
  if (test%is_fatal .and. num_PEs() > 1) return

  if (test%is_fatal) then
    rc = sigsetjmp(env, 1)
    if (rc == 0) then
      call disable_fatal_errors(env)
      call test%proc
    endif
    call enable_fatal_errors
  else
    call test%proc
  endif

  if (associated(test%cleanup)) call test%cleanup
end subroutine run_unit_test


!> Return a new test suite
function create_test_suite() result(suite)
  type(TestSuite) :: suite

  ! Setup the head node, but do not populate it
  allocate(suite%head)
  suite%tail => suite%head
end function create_test_suite


subroutine add_unit_test_basic(suite, test, name, fatal)
  class(TestSuite), intent(inout) :: suite
  procedure() :: test
  character(len=*), intent(in) :: name
  logical, intent(in), optional :: fatal

  procedure(), pointer :: cleanup

  cleanup => null()
  if (associated(suite%cleanup)) cleanup => suite%cleanup

  call add_unit_test_full(suite, test, name, fatal, cleanup)
end subroutine add_unit_test_basic


subroutine add_unit_test_full(suite, test, name, fatal, cleanup)
  class(TestSuite), intent(inout) :: suite
  procedure() :: test
  character(len=*), intent(in) :: name
  procedure() :: cleanup
  logical, intent(in), optional :: fatal

  type(UnitTest), pointer :: utest
  type(UnitTestNode), pointer :: node

  ! Populate the current tail
  allocate(utest)
  utest = UnitTest(test, name, fatal, cleanup)
  suite%tail%test => utest

  ! Create and append the new (empty) node, and update the tail
  allocate(node)
  suite%tail%next => node
  suite%tail => node
end subroutine add_unit_test_full


subroutine run_test_suite(suite)
  class(TestSuite), intent(in) :: suite

  type(UnitTestNode), pointer :: node

  node => suite%head
  do while(associated(node%test))
    ! TODO: Capture FMS stdout/stderr
    print '(/a)', "=== "//node%test%name

    call node%test%run
    if (associated(node%test%cleanup)) call node%test%cleanup

    node => node%next
  enddo
end subroutine run_test_suite


!> Initialize string with a character array.
function init_string_char(c) result(str)
  character(len=*), dimension(:), intent(in) :: c
    !< List of character arrays
  type(string), dimension(size(c)) :: str
    !< String output

  integer :: i

  do i = 1, size(c)
    str(i)%s = c(i)
  enddo
end function init_string_char


!> Convert an integer to a string
function init_string_int(n) result(str)
  integer, intent(in) :: n
    !< Integer input
  type(string) :: str
    !< String output

  ! TODO: Estimate this with integer arithmetic
  character(1 + floor(log10(real(abs(n)))) + (1 - sign(1, n))/2) :: chr

  write(chr, '(i0)') n
  str = string(chr)
end function init_string_int


!> Create a text file for unit testing
subroutine create_test_file(filename, lines, mode)
  character(len=*), intent(in) :: filename
    !< Name of file to be created
  type(string), intent(in), optional :: lines(:)
    !< list of strings to write to file
  integer, optional, intent(in) :: mode
    !< Permissions of new file

  integer :: param_unit
  integer :: i
  integer :: rc
  logical :: sync

  if (is_root_PE()) then
    open(newunit=param_unit, file=filename, status='replace')
    if (present(lines)) then
      do i = 1, size(lines)
        write(param_unit, '(a)') lines(i)%s
      enddo
    endif
    close(param_unit)
    if (present(mode)) rc = chmod(filename, mode)
  endif
  call sync_PEs
end subroutine create_test_file


!> Delete a file created during testing
subroutine delete_test_file(filename)
  character(len=*), intent(in) :: filename
    !< Name of file to be deleted

  logical :: is_file, is_open
  integer :: io_unit

  if (is_root_PE()) then
    inquire(file=filename, exist=is_file, opened=is_open, number=io_unit)

    if (is_file) then
      if (.not. is_open) open(newunit=io_unit, file=filename)
      close(io_unit, status='delete')
    endif
  endif
  call sync_PEs
end subroutine delete_test_file

end module MOM_unit_testing
