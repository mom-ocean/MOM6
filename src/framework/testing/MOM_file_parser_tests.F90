module MOM_file_parser_tests

use posix, only : chmod

use MOM_file_parser, only : param_file_type
use MOM_file_parser, only : open_param_file
use MOM_file_parser, only : close_param_file
use MOM_file_parser, only : read_param
use MOM_file_parser, only : log_param
use MOM_file_parser, only : get_param
use MOM_file_parser, only : log_version
use MOM_file_parser, only : clearParameterBlock
use MOM_file_parser, only : openParameterBlock
use MOM_file_parser, only : closeParameterBlock

use MOM_time_manager, only : time_type
use MOM_time_manager, only : set_date
use MOM_time_manager, only : set_ticks_per_second
use MOM_time_manager, only : set_calendar_type
use MOM_time_manager, only : NOLEAP, NO_CALENDAR

use MOM_error_handler, only : assert
use MOM_error_handler, only : MOM_error
use MOM_error_handler, only : FATAL

use MOM_unit_testing, only : TestSuite
use MOM_unit_testing, only : string
use MOM_unit_testing, only : create_test_file
use MOM_unit_testing, only : delete_test_file

implicit none ; private

public :: run_file_parser_tests

character(len=*), parameter :: param_filename = 'TEST_input'
character(len=*), parameter :: missing_param_filename = 'MISSING_input'
character(len=*), parameter :: netcdf_param_filename = 'TEST_input.nc'

character(len=*), parameter :: sample_param_name = 'SAMPLE_PARAMETER'
character(len=*), parameter :: missing_param_name = 'MISSING_PARAMETER'

character(len=*), parameter :: module_name = "SAMPLE_module"
character(len=*), parameter :: module_version = "SAMPLE_version"
character(len=*), parameter :: module_desc = "Description here"

character(len=9), parameter :: param_docfiles(4) = [ &
  "all      ", &
  "debugging", &
  "layout   ", &
  "short    " &
]

contains

subroutine test_open_param_file
  type(param_file_type) :: param

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call close_param_file(param)
end subroutine test_open_param_file


subroutine test_close_param_file_quiet
  type(param_file_type) :: param

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call close_param_file(param, quiet_close=.true.)
end subroutine test_close_param_file_quiet


subroutine test_open_param_file_component
  type(param_file_type) :: param
  integer :: i

  call create_test_file(param_filename)

  call open_param_file(param_filename, param, component="TEST")
  call close_param_file(param, component="TEST")
end subroutine test_open_param_file_component


subroutine cleanup_open_param_file_component
  integer :: i

  call delete_test_file(param_filename)
  do i = 1, 4
    call delete_test_file("TEST_parameter_doc."//param_docfiles(i))
  enddo
end subroutine cleanup_open_param_file_component


subroutine test_open_param_file_docdir
  ! TODO: Make a new directory...?
  type(param_file_type) :: param

  call create_test_file(param_filename)

  call open_param_file(param_filename, param, doc_file_dir='./')
  call close_param_file(param)
end subroutine test_open_param_file_docdir


subroutine test_open_param_file_empty_filename
  type(param_file_type) :: param

  call open_param_file('', param)
  ! FATAL; return to program
end subroutine test_open_param_file_empty_filename


subroutine test_open_param_file_long_name
  !> Store filename in a variable longer than FILENAME_LENGTH
  type(param_file_type) :: param
  character(len=250) :: long_filename

  long_filename = param_filename

  call create_test_file(long_filename)

  call open_param_file(long_filename, param)
  call close_param_file(param)
end subroutine test_open_param_file_long_name


subroutine test_missing_param_file
  type(param_file_type) :: param
  logical :: file_exists

  inquire(file=missing_param_filename, exist=file_exists)
  if (file_exists) call MOM_error(FATAL, "Missing file already exists!")

  call open_param_file(missing_param_filename, param)
  ! FATAL; return to program
end subroutine test_missing_param_file


subroutine test_open_param_file_ioerr
  type(param_file_type) :: param
  ! NOTE: Induce an I/O error in open() by making the file unreadable

  call create_test_file(param_filename, mode=int(o'000'))

  call open_param_file(param_filename, param)
  ! FATAL; return to program
end subroutine test_open_param_file_ioerr


subroutine cleanup_open_param_file_ioerr
  integer :: rc

  rc = chmod(param_filename, int(o'700'))
  call cleanup_file_parser()
end subroutine cleanup_open_param_file_ioerr


subroutine test_open_param_file_netcdf
  type(param_file_type) :: param

  call create_test_file(netcdf_param_filename)

  call open_param_file(netcdf_param_filename, param)
  ! FATAL; return to program
end subroutine test_open_param_file_netcdf


subroutine cleanup_open_param_file_netcdf
  integer :: param_unit
  logical :: is_open

  call delete_test_file(netcdf_param_filename)
end subroutine cleanup_open_param_file_netcdf


subroutine test_open_param_file_checkable
  type(param_file_type) :: param

  call create_test_file(param_filename)

  call open_param_file(param_filename, param, checkable=.false.)
  call close_param_file(param)
end subroutine test_open_param_file_checkable


subroutine test_reopen_param_file
  type(param_file_type) :: param

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call open_param_file(param_filename, param)
  call close_param_file(param)
end subroutine test_reopen_param_file


subroutine test_open_param_file_no_doc
  type(param_file_type) :: param
  type(string) :: lines(1)

  lines(1) = string('DOCUMENT_FILE = ""')
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call close_param_file(param)
end subroutine test_open_param_file_no_doc


subroutine test_read_param_int
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(1)
  character(len=*), parameter :: sample_input = '123'
  integer, parameter :: sample_result = 123

  lines = string(sample_param_name // ' = ' // sample_input)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(sample == sample_result, 'Incorrect value')
end subroutine test_read_param_int


subroutine test_read_param_int_missing
  type(param_file_type) :: param
  integer :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call read_param(param, missing_param_name, sample, fail_if_missing=.true.)
  ! FATAL; return to program
end subroutine test_read_param_int_missing


subroutine test_read_param_int_undefined
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(1)

  lines = string('#undef ' // sample_param_name)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample, fail_if_missing=.true.)
  ! FATAL; return to program
end subroutine test_read_param_int_undefined


subroutine test_read_param_int_type_err
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(1)

  lines = string(sample_param_name // ' = not_an_integer')
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_int_type_err


subroutine test_read_param_int_array
  type(param_file_type) :: param
  integer :: sample(3)
  type(string) :: lines(1)
  character(len=*), parameter :: sample_input = '1, 2, 3'
  integer, parameter :: sample_result(3) = [1, 2, 3]

  lines = string(sample_param_name // ' = ' // sample_input)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(all(sample == sample_result), 'Incorrect value')
end subroutine test_read_param_int_array


subroutine test_read_param_int_array_missing
  type(param_file_type) :: param
  integer :: sample(3)

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call read_param(param, missing_param_name, sample, fail_if_missing=.true.)
  ! FATAL; return to program
end subroutine test_read_param_int_array_missing


subroutine test_read_param_int_array_undefined
  type(param_file_type) :: param
  integer :: sample(3)
  type(string) :: lines(1)

  lines = string('#undef ' // sample_param_name)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample, fail_if_missing=.true.)
  ! FATAL; return to program
end subroutine test_read_param_int_array_undefined


subroutine test_read_param_int_array_type_err
  type(param_file_type) :: param
  integer :: sample(3)
  type(string) :: lines(1)

  lines = string(sample_param_name // ' = not_an_int_array')
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_int_array_type_err


subroutine test_read_param_real
  type(param_file_type) :: param
  real :: sample
  type(string) :: lines(1)
  character(len=*), parameter :: sample_input = '3.14'
  real, parameter :: sample_result = 3.14

  lines = string(sample_param_name // ' = ' // sample_input)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(sample == sample_result, 'Incorrect value')
end subroutine test_read_param_real


subroutine test_read_param_real_missing
  type(param_file_type) :: param
  real :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call read_param(param, missing_param_name, sample, fail_if_missing=.true.)
  ! FATAL; return to program
end subroutine test_read_param_real_missing


subroutine test_read_param_real_undefined
  type(param_file_type) :: param
  real :: sample
  type(string) :: lines(1)

  lines = string('#undef ' // sample_param_name)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample, fail_if_missing=.true.)
  ! FATAL; return to program
end subroutine test_read_param_real_undefined


subroutine test_read_param_real_type_err
  type(param_file_type) :: param
  real :: sample
  type(string) :: lines(1)

  lines = string(sample_param_name // ' = not_a_real')
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_real_type_err


subroutine test_read_param_real_array
  type(param_file_type) :: param
  real :: sample(3)
  type(string) :: lines(1)
  character(len=*), parameter :: sample_input = '1., 2., 3.'
  real, parameter :: sample_result(3) = [1., 2., 3.]

  lines = string(sample_param_name // ' = ' // sample_input)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(all(sample == sample_result), 'Incorrect value')
end subroutine test_read_param_real_array


subroutine test_read_param_real_array_missing
  type(param_file_type) :: param
  real :: sample(3)

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call read_param(param, missing_param_name, sample, fail_if_missing=.true.)
  ! FATAL; return to program
end subroutine test_read_param_real_array_missing


subroutine test_read_param_real_array_undefined
  type(param_file_type) :: param
  real :: sample(3)
  type(string) :: lines(1)

  lines = string('#undef ' // sample_param_name)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample, fail_if_missing=.true.)
  ! FATAL; return to program
end subroutine test_read_param_real_array_undefined


subroutine test_read_param_real_array_type_err
  type(param_file_type) :: param
  real :: sample(3)
  type(string) :: lines(1)

  lines = string(sample_param_name // ' = not_a_real_array')
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_real_array_type_err


subroutine test_read_param_logical
  type(param_file_type) :: param
  logical :: sample
  type(string) :: lines(1)
  character(len=*), parameter :: sample_input = 'True'
  logical, parameter :: sample_result = .true.

  lines = string(sample_param_name // ' = ' // sample_input)

  !lines = string(sample_param_name // ' = True')
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(sample .eqv. sample_result, 'Incorrect value')
end subroutine test_read_param_logical


subroutine test_read_param_logical_missing
  type(param_file_type) :: param
  logical :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call read_param(param, missing_param_name, sample, fail_if_missing=.true.)
  ! FATAL; return to program
end subroutine test_read_param_logical_missing


subroutine test_read_param_char_no_delim
  type(param_file_type) :: param
  character(len=8) :: sample
  type(string) :: lines(1)
  character(len=*), parameter :: sample_input = "abcdefgh"
  character(len=*), parameter :: sample_result = "abcdefgh"

  lines = string(sample_param_name // ' = ' // sample_input)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(sample == sample_result, 'Incorrect value')
end subroutine test_read_param_char_no_delim


subroutine test_read_param_char_quote_delim
  type(param_file_type) :: param
  character(len=8) :: sample
  type(string) :: lines(1)
  character(len=*), parameter :: sample_input = '"abcdefgh"'
  character(len=*), parameter :: sample_result = "abcdefgh"

  lines = string(sample_param_name // ' = ' // sample_input)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(sample == sample_result, 'Incorrect value')
end subroutine test_read_param_char_quote_delim


subroutine test_read_param_char_apostrophe_delim
  type(param_file_type) :: param
  character(len=8) :: sample
  type(string) :: lines(1)
  character(len=*), parameter :: sample_input = "'abcdefgh'"
  character(len=*), parameter :: sample_result = "abcdefgh"

  lines = string(sample_param_name // " = " // sample_input)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(sample == sample_result, 'Incorrect value')
end subroutine test_read_param_char_apostrophe_delim


subroutine test_read_param_char_missing
  type(param_file_type) :: param
  character(len=8) :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call read_param(param, missing_param_name, sample, fail_if_missing=.true.)
  ! FATAL; return to program
end subroutine test_read_param_char_missing


subroutine test_read_param_char_array
  type(param_file_type) :: param
  character(len=3) :: sample(3)
  type(string) :: lines(1)
  character(len=*), parameter :: sample_input = '"abc", "def", "ghi"'
  character(len=*), parameter :: sample_result(3) = ["abc", "def", "ghi"]

  lines = string(sample_param_name // ' = ' // sample_input)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(all(sample == sample_result), 'Incorrect value')
end subroutine test_read_param_char_array


subroutine test_read_param_char_array_missing
  type(param_file_type) :: param
  character(len=8) :: sample(3)

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call read_param(param, missing_param_name, sample, fail_if_missing=.true.)
  ! FATAL; return to program
end subroutine test_read_param_char_array_missing


subroutine test_read_param_time_date
  type(param_file_type) :: param
  type(time_type) :: sample
  type(string) :: lines(1)

  lines = string(sample_param_name // ' = 1980-01-01 00:00:00')
  call create_test_file(param_filename, lines)

  call set_calendar_type(NOLEAP)
  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)
end subroutine test_read_param_time_date


subroutine test_read_param_time_date_bad_format
  type(param_file_type) :: param
  type(time_type) :: sample
  type(string) :: lines(1)

  lines = string(sample_param_name // ' = 1980--01--01 00::00::00')
  call create_test_file(param_filename, lines)

  call set_calendar_type(NOLEAP)
  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_time_date_bad_format


subroutine test_read_param_time_tuple
  type(param_file_type) :: param
  type(time_type) :: sample
  type(string) :: lines(1)

  lines = string(sample_param_name // ' = 1980,1,1,0,0,0')
  call create_test_file(param_filename, lines)

  call set_calendar_type(NOLEAP)
  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)
end subroutine test_read_param_time_tuple


subroutine test_read_param_time_bad_tuple
  type(param_file_type) :: param
  type(time_type) :: sample
  type(string) :: lines(1)

  lines = string(sample_param_name // ' = 1980, 1')
  call create_test_file(param_filename, lines)

  call set_calendar_type(NOLEAP)
  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_time_bad_tuple


subroutine test_read_param_time_bad_tuple_values
  type(param_file_type) :: param
  type(time_type) :: sample
  type(string) :: lines(1)

  lines = string(sample_param_name // ' = 0, 0, 0, 0, 0, 0')
  call create_test_file(param_filename, lines)

  call set_calendar_type(NOLEAP)
  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_time_bad_tuple_values


subroutine test_read_param_time_unit
  type(param_file_type) :: param
  type(time_type) :: sample
  type(string) :: lines(1)

  lines = string(sample_param_name // ' = 0.5')
  call create_test_file(param_filename, lines)

  call set_calendar_type(NOLEAP)
  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample, timeunit=86400.)
  call close_param_file(param)
end subroutine test_read_param_time_unit


subroutine test_read_param_time_missing
  type(param_file_type) :: param
  type(time_type) :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call read_param(param, missing_param_name, sample, fail_if_missing=.true.)
  ! FATAL; return to program
end subroutine test_read_param_time_missing


subroutine test_read_param_time_undefined
  type(param_file_type) :: param
  type(time_type) :: sample
  type(string) :: lines(1)

  lines = string('#undef ' // sample_param_name)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample, fail_if_missing=.true.)
  ! FATAL; return to program
end subroutine test_read_param_time_undefined


subroutine test_read_param_time_type_err
  type(param_file_type) :: param
  type(time_type) :: sample
  type(string) :: lines(1)

  lines = string(sample_param_name // ' = 1., 2., 3., 4., 5., 6.')
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_time_type_err

! Generic parameter tests

subroutine test_read_param_unused_fatal
  type(param_file_type) :: param
  type(string) :: lines(2)

  lines = [ &
      string('FATAL_UNUSED_PARAMS = True'), &
      string(sample_param_name // ' = 1') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call close_param_file(param)
  ! FATAL; return to program
end subroutine test_read_param_unused_fatal


subroutine test_read_param_replace_tabs
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(1)
  character(len=*), parameter :: sample_input = "1"
  integer, parameter :: sample_result = 1
  character, parameter :: tab = achar(9)

  lines = string(sample_param_name // tab // '=' // tab // sample_input)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(sample == sample_result, 'Incorrect value')
end subroutine test_read_param_replace_tabs


subroutine test_read_param_pad_equals
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(1)
  character(len=*), parameter :: sample_input = "1"
  integer, parameter :: sample_result = 1

  lines = string(sample_param_name // '=' // sample_input)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(sample == sample_result, 'Incorrect value')
end subroutine test_read_param_pad_equals


subroutine test_read_param_multiline_param
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(2)
  integer, parameter :: sample_result = 1
  character, parameter :: backslash = achar(92)

  lines = [ &
      string(sample_param_name // ' = ' // backslash), &
      string('  1') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(sample == sample_result, 'Incorrect result')
end subroutine test_read_param_multiline_param


subroutine test_read_param_multiline_param_unclosed
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(1)
  character, parameter :: backslash = achar(92)

  lines = string(sample_param_name // ' = ' // backslash)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  ! FATAL; return to program
end subroutine test_read_param_multiline_param_unclosed


subroutine test_read_param_multiline_comment
  type(param_file_type) :: param
  integer :: sample

  type(string) :: lines(6)

  lines = [ &
      string('/* First C comment line'), &
      string('   Second C comment line */'), &
      string('// First C++ comment line'), &
      string('// Second C++ comment line'), &
      string('! First Fortran comment line'), &
      string('! Second Fortran comment line') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call close_param_file(param)
end subroutine test_read_param_multiline_comment


subroutine test_read_param_multiline_comment_unclosed
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(1)

  lines = string('/* Unclosed C comment')
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  ! FATAL; return to program
end subroutine test_read_param_multiline_comment_unclosed


subroutine test_read_param_misplaced_quote
  type(param_file_type) :: param
  character(len=20) :: sample
  type(string) :: lines(1)

  lines = string(sample_param_name // ' = "abc')
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  ! FATAL; return to program
end subroutine test_read_param_misplaced_quote


subroutine test_read_param_define
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(1)
  integer, parameter :: sample_result = 2

  lines = string('#define ' // sample_param_name // ' 2')
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(sample == sample_result, 'Incorrect value')
end subroutine test_read_param_define


subroutine test_read_param_define_as_flag
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(1)

  lines = string('#define ' // sample_param_name)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)
end subroutine test_read_param_define_as_flag


subroutine test_read_param_override
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(2)
  integer, parameter :: sample_result = 2

  lines = [ &
      string(sample_param_name // ' = 1'), &
      string('#override ' // sample_param_name // ' = 2') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)

  call assert(sample == sample_result, 'Incorrect value')
end subroutine test_read_param_override


subroutine test_read_param_override_misplaced
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(1)

  lines(1) = string('#define #override ' // sample_param_name // ' = 1')
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  ! FATAL; return to program
end subroutine test_read_param_override_misplaced


subroutine test_read_param_override_twice
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(3)

  lines = [ &
      string(sample_param_name // ' = 1'), &
      string('#override ' // sample_param_name // ' = 2'), &
      string('#override ' // sample_param_name // ' = 3') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_override_twice


subroutine test_read_param_override_repeat
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(3)

  lines = [ &
      string(sample_param_name // ' = 1'), &
      string('#override ' // sample_param_name // ' = 2'), &
      string('#override ' // sample_param_name // ' = 2') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_override_repeat


subroutine test_read_param_override_warn_chain
  type(param_file_type) :: param
  integer :: sample
  character(len=*), parameter :: other_param_name = 'OTHER_PARAMETER'
  type(string) :: lines(4)

  lines = [ &
      string(other_param_name // ' = 1'), &
      string(sample_param_name // ' = 2'), &
      string('#override ' // other_param_name // ' = 3'), &
      string('#override ' // sample_param_name // ' = 4') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  ! First invoke the "other" override, adding it to the chain
  call read_param(param, other_param_name, sample)
  ! Now invoke the "sample" override, with "other" in the chain
  call read_param(param, sample_param_name, sample)
  ! Finally, re-invoke the "other" override, having already been issued.
  call read_param(param, other_param_name, sample)
  call close_param_file(param)
end subroutine test_read_param_override_warn_chain


subroutine test_read_param_assign_after_override
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(2)

  lines = [ &
      string('#override ' // sample_param_name // ' = 2'), &
      string(sample_param_name // ' = 3') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)
end subroutine test_read_param_assign_after_override


subroutine test_read_param_override_no_def
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(1)

  lines(1) = string('#override ' // sample_param_name)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_override_no_def


subroutine test_read_param_assign_twice
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(2)

  lines = [ &
      string(sample_param_name // ' = 1'), &
      string(sample_param_name // ' = 2') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_assign_twice


subroutine test_read_param_assign_repeat
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(2)

  lines = [ &
      string(sample_param_name // ' = 1'), &
      string(sample_param_name // ' = 1') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call close_param_file(param)
end subroutine test_read_param_assign_repeat


subroutine test_read_param_null_stmt
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(1)

  lines(1) = string(sample_param_name)
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_null_stmt


subroutine test_read_param_assign_in_define
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(1)

  lines = string('#define ' // sample_param_name // ' = 1')
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  ! FATAL; return to program
end subroutine test_read_param_assign_in_define

!-- Blocks

subroutine test_read_param_block
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(3)
  integer, parameter :: sample_result = 123

  lines = [ &
      string('ABC%'), &
      string('ABC%' // sample_param_name // ' = 123'), &
      string('%ABC') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call openParameterBlock(param, 'ABC')
  call read_param(param, sample_param_name, sample)
  call closeParameterBlock(param)
  call clearParameterBlock(param)
  call close_param_file(param)

  call assert(sample == sample_result, 'Incorrect value')
end subroutine test_read_param_block


! TODO: This test fails due to an implementation issue.
subroutine test_read_param_block_stack
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(5)

  lines = [ &
      string('ABC%'), &
      string('DEF%'), &
      string(sample_param_name // ' = 123'), &
      string('DEF%'), &
      string('%ABC') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call openParameterBlock(param, 'ABC')
  call openParameterBlock(param, 'DEF')
  call read_param(param, sample_param_name, sample)
  call closeParameterBlock(param)
  call clearParameterBlock(param)
  call close_param_file(param)
end subroutine test_read_param_block_stack


! NOTE: This is a simpler version of the block_stack test which works
subroutine test_read_param_block_inline_stack
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(3)

  lines = [ &
      string('ABC%'), &
      string('DEF%' // sample_param_name // ' = 123'), &
      string('%ABC') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call openParameterBlock(param, 'ABC')
  call openParameterBlock(param, 'DEF')
  call read_param(param, sample_param_name, sample)
  call closeParameterBlock(param)
  call clearParameterBlock(param)
  call close_param_file(param)
end subroutine test_read_param_block_inline_stack


subroutine test_read_param_block_empty_pop
  type(param_file_type) :: param

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call openParameterBlock(param, '%')
  call openParameterBlock(param, 'ABC')
  call closeParameterBlock(param)
  call closeParameterBlock(param)
  ! FATAL; return to program
end subroutine test_read_param_block_empty_pop


subroutine test_read_param_block_close_unnamed
  type(param_file_type) :: param
  type(string) :: lines(2)

  lines = [ &
      string('ABC%'), &
      string('%ABC') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call openParameterBlock(param, 'ABC')
  call closeParameterBlock(param)
  call closeParameterBlock(param)
  ! FATAL; return to program
end subroutine test_read_param_block_close_unnamed


subroutine test_read_param_block_close_unopened
  type(param_file_type) :: param
  type(string) :: lines(1)

  lines = string('%CBA')
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  ! FATAL; return to program
end subroutine test_read_param_block_close_unopened


subroutine test_read_param_block_unmatched
  type(param_file_type) :: param
  type(string) :: lines(2)

  lines = [ &
      string('ABC%'), &
      string('%CBA') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  ! FATAL; return to program
end subroutine test_read_param_block_unmatched


subroutine test_open_unallocated_block
  type(param_file_type) :: param
  character(len=*), parameter :: block_name = "ABC"

  call openParameterBlock(param, block_name)
  ! FATAL; return to program
end subroutine test_open_unallocated_block


subroutine test_close_unallocated_block
  type(param_file_type) :: param

  call closeParameterBlock(param)
  ! FATAL; return to program
end subroutine test_close_unallocated_block


subroutine test_clear_unallocated_block
  type(param_file_type) :: param

  call clearParameterBlock(param)
  ! FATAL; return to program
end subroutine test_clear_unallocated_block


subroutine test_read_param_block_outside_block
  type(param_file_type) :: param
  integer :: sample
  type(string) :: lines(3)

  lines = [ &
      string('ABC%'), &
      string(sample_param_name // ' = 1'), &
      string('%ABC') &
  ]
  call create_test_file(param_filename, lines)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
end subroutine test_read_param_block_outside_block

!---

subroutine test_log_version_cs
  type(param_file_type) :: param

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call log_version(param, module_name, module_version, desc=module_desc)
  call close_param_file(param)
end subroutine test_log_version_cs


subroutine test_log_version_plain
  call log_version(module_name, module_version)
end subroutine test_log_version_plain


subroutine test_log_param_int
  type(param_file_type) :: param
  integer, parameter :: sample = 1
  character(len=*), parameter :: desc = "Parameter description"

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call log_param(param, module_name, sample_param_name, sample, desc=desc)
  call close_param_file(param)
end subroutine test_log_param_int


subroutine test_log_param_int_array
  type(param_file_type) :: param
  integer, parameter :: sample(3) = [1, 2, 3]
  character(len=*), parameter :: desc = "Parameter description"

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call log_param(param, module_name, sample_param_name, sample, desc=desc)
  call close_param_file(param)
end subroutine test_log_param_int_array


subroutine test_log_param_real
  type(param_file_type) :: param
  real, parameter :: sample = 1.
  character(len=*), parameter :: desc = "Parameter description"

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call log_param(param, module_name, sample_param_name, sample, desc=desc)
  call close_param_file(param)
end subroutine test_log_param_real


subroutine test_log_param_real_array
  type(param_file_type) :: param
  real, parameter :: sample(3) = [1., 2., 3.]
  character(len=*), parameter :: desc = "Parameter description"

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call log_param(param, module_name, sample_param_name, sample, desc=desc)
  call close_param_file(param)
end subroutine test_log_param_real_array


subroutine test_log_param_time
  type(param_file_type) :: param
  type(time_type) :: sample
  character(len=*), parameter :: desc = "Parameter description"
  type(string) :: lines(1)

  lines = string(sample_param_name // ' = 1980,1,1,0,0,0')

  call set_calendar_type(NOLEAP)
  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call read_param(param, sample_param_name, sample)
  call log_param(param, module_name, sample_param_name, sample, desc=desc)
  call close_param_file(param)
end subroutine test_log_param_time


subroutine test_log_param_time_as_date
  type(param_file_type) :: param
  type(time_type) :: sample
  character(len=*), parameter :: desc = "Parameter description"

  call set_calendar_type(NOLEAP)
  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  sample = set_date(1980, 1, 1, 0, 0, 0)
  call log_param(param, module_name, sample_param_name, sample, desc=desc, &
      log_date=.true.)
  call close_param_file(param)
end subroutine test_log_param_time_as_date


subroutine test_log_param_time_as_date_default
  type(param_file_type) :: param
  type(time_type) :: sample
  type(time_type) :: default_date
  character(len=*), parameter :: desc = "Parameter description"

  call set_calendar_type(NOLEAP)
  call create_test_file(param_filename)

  call open_param_file(param_filename, param)

  call set_ticks_per_second(60)
  default_date = set_date(1980, 1, 1, 0, 0, 0, 30)
  call log_param(param, module_name, sample_param_name, sample, desc=desc, &
      log_date=.true., default=default_date)

  call set_ticks_per_second(300)
  default_date = set_date(1980, 1, 1, 0, 0, 0, 150)
  call log_param(param, module_name, sample_param_name, sample, desc=desc, &
      log_date=.true., default=default_date)

  call close_param_file(param)
end subroutine test_log_param_time_as_date_default


subroutine test_log_param_time_as_date_tick
  type(param_file_type) :: param
  type(time_type) :: sample
  character(len=*), parameter :: desc = "Parameter description"

  call set_calendar_type(NOLEAP)
  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call log_param(param, module_name, sample_param_name, sample, desc=desc, &
      log_date=.true.)
  call close_param_file(param)
end subroutine test_log_param_time_as_date_tick


subroutine test_log_param_time_with_unit
  type(param_file_type) :: param
  type(time_type) :: sample
  type(time_type) :: default_date
  character(len=*), parameter :: desc = "Parameter description"
  character(len=*), parameter :: sample_units = "days since whatever"

  call set_calendar_type(NOLEAP)
  call create_test_file(param_filename)

  call set_ticks_per_second(60)
  sample = set_date(1980, 1, 1, 0, 0, 0, 30)

  default_date = set_date(1980, 1, 1, 0, 0, 0, 30)

  call open_param_file(param_filename, param)
  call log_param(param, module_name, sample_param_name, sample, desc=desc, &
      units=sample_units, timeunit=86400., default=default_date)
  call close_param_file(param)
end subroutine test_log_param_time_with_unit


subroutine test_log_param_time_with_timeunit
  type(param_file_type) :: param
  type(time_type) :: sample
  integer :: i
  character(len=*), parameter :: desc = "Parameter description"
  real, parameter :: timeunits(5) = [1., 3600., 86400., 3.1e7, 1e8]

  call set_calendar_type(NOLEAP)
  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  do i = 1,5
    call log_param(param, module_name, sample_param_name, sample, desc=desc, &
        timeunit=timeunits(i))
  enddo
  call close_param_file(param)
end subroutine test_log_param_time_with_timeunit

!----

subroutine test_get_param_int
  type(param_file_type) :: param
  integer :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample)
  call close_param_file(param)
end subroutine test_get_param_int


subroutine test_get_param_int_no_read_no_log
  type(param_file_type) :: param
  integer :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample, &
      do_not_read=.true., do_not_log=.true.)
  call close_param_file(param)
end subroutine test_get_param_int_no_read_no_log


subroutine test_get_param_int_array
  type(param_file_type) :: param
  integer :: sample(3)

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample)
  call close_param_file(param)
end subroutine test_get_param_int_array


subroutine test_get_param_int_array_no_read_no_log
  type(param_file_type) :: param
  integer :: sample(3)

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample, &
      do_not_read=.true., do_not_log=.true.)
  call close_param_file(param)
end subroutine test_get_param_int_array_no_read_no_log


subroutine test_get_param_real
  type(param_file_type) :: param
  real :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample)
  call close_param_file(param)
end subroutine test_get_param_real


subroutine test_get_param_real_no_read_no_log
  type(param_file_type) :: param
  real :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample, &
      do_not_read=.true., do_not_log=.true.)
  call close_param_file(param)
end subroutine test_get_param_real_no_read_no_log


subroutine test_get_param_real_array
  type(param_file_type) :: param
  real :: sample(3)

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample)
  call close_param_file(param)
end subroutine test_get_param_real_array


subroutine test_get_param_real_array_no_read_no_log
  type(param_file_type) :: param
  real :: sample(3)

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample, &
      do_not_read=.true., do_not_log=.true.)
  call close_param_file(param)
end subroutine test_get_param_real_array_no_read_no_log


subroutine test_get_param_char
  type(param_file_type) :: param
  character(len=8) :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample)
  call close_param_file(param)
end subroutine test_get_param_char


subroutine test_get_param_char_no_read_no_log
  type(param_file_type) :: param
  character(len=8) :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample, &
      do_not_read=.true., do_not_log=.true.)
  call close_param_file(param)
end subroutine test_get_param_char_no_read_no_log


subroutine test_get_param_char_array
  type(param_file_type) :: param
  character(len=8) :: sample(3)

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample)
  call close_param_file(param)
end subroutine test_get_param_char_array


subroutine test_get_param_logical
  type(param_file_type) :: param
  logical :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample)
  call close_param_file(param)
end subroutine test_get_param_logical


subroutine test_get_param_logical_no_read_no_log
  type(param_file_type) :: param
  logical :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample, &
      do_not_read=.true., do_not_log=.true.)
  call close_param_file(param)
end subroutine test_get_param_logical_no_read_no_log


subroutine test_get_param_logical_default
  type(param_file_type) :: param
  logical :: sample
  logical, parameter :: default_value = .false.

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample, &
      default=default_value)
  call close_param_file(param)
end subroutine test_get_param_logical_default


subroutine test_get_param_time
  type(param_file_type) :: param
  type(time_type) :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample)
  call close_param_file(param)
end subroutine test_get_param_time


subroutine test_get_param_time_no_read_no_log
  type(param_file_type) :: param
  type(time_type) :: sample

  call create_test_file(param_filename)

  call open_param_file(param_filename, param)
  call get_param(param, module_name, sample_param_name, sample, &
      do_not_read=.true., do_not_log=.true.)
  call close_param_file(param)
end subroutine test_get_param_time_no_read_no_log


! Utility functions
! TODO: Move to a generic testing module

subroutine cleanup_file_parser
  integer :: i

  call delete_test_file(param_filename)
  do i = 1, 4
    call delete_test_file("MOM_parameter_doc."//param_docfiles(i))
  enddo

  call set_calendar_type(NO_CALENDAR)
end subroutine cleanup_file_parser


subroutine run_file_parser_tests
  ! testing...
  type(TestSuite) :: suite

  ! Delete any pre-existing test parameter files
  call cleanup_file_parser

  ! Build the test suite
  suite = TestSuite()
  suite%cleanup => cleanup_file_parser

  call suite%add(test_open_param_file, "test_open_param_file")

  call suite%add(test_close_param_file_quiet, "test_close_param_file_quiet")

  call suite%add(test_open_param_file_component, "test_open_param_file_component", &
      cleanup=cleanup_open_param_file_component)

  call suite%add(test_open_param_file_docdir, "test_open_param_file_docdir")

  call suite%add(test_open_param_file_empty_filename, &
      "test_open_param_file_empty_filename", fatal=.true.)

  call suite%add(test_open_param_file_long_name, &
      "test_open_param_file_longname")

  call suite%add(test_missing_param_file, "test_missing_param_file", &
      fatal=.true.)

  call suite%add(test_open_param_file_ioerr, "test_open_param_file_ioerr", &
      fatal=.true., cleanup=cleanup_open_param_file_ioerr)

  call suite%add(test_open_param_file_checkable, &
      "test_open_param_file_checkable")

  call suite%add(test_reopen_param_file, "test_reopen_param_file")

  call suite%add(test_open_param_file_netcdf, "test_open_param_file_netcdf", &
      fatal=.true., cleanup=cleanup_open_param_file_netcdf)

  call suite%add(test_open_param_file_no_doc, "test_open_param_file_no_doc")

  call suite%add(test_read_param_int, "test_read_param_int")

  call suite%add(test_read_param_int_missing, "test_read_param_int_missing", &
      fatal=.true.)

  call suite%add(test_read_param_int_undefined, &
      "test_read_param_int_undefined", fatal=.true.)

  call suite%add(test_read_param_int_type_err, &
      "test_read_param_int_type_err", fatal=.true.)

  call suite%add(test_read_param_int_array, "test_read_param_int_array")

  call suite%add(test_read_param_int_array_missing, &
      "test_read_param_int_array_missing", fatal=.true.)

  call suite%add(test_read_param_int_array_undefined, &
    "test_read_param_int_array_undefined", fatal=.true.)

  call suite%add(test_read_param_int_array_type_err, &
      "test_read_param_int_array_type_err", fatal=.true.)

  call suite%add(test_read_param_real, "test_read_param_real")

  call suite%add(test_read_param_real_missing, &
      "test_read_param_real_missing", fatal=.true.)

  call suite%add(test_read_param_real_undefined, &
      "test_read_param_real_undefined", fatal=.true.)

  call suite%add(test_read_param_real_type_err, &
      "test_read_param_real_type_err", fatal=.true.)

  call suite%add(test_read_param_real_array, "test_read_param_real_array")

  call suite%add(test_read_param_real_array_missing, &
      "test_read_param_real_array_missing", fatal=.true.)

  call suite%add(test_read_param_real_array_undefined, &
      "test_read_param_real_array_undefined", fatal=.true.)

  call suite%add(test_read_param_real_array_type_err, &
      "test_read_param_real_array_type_err", fatal=.true.)

  call suite%add(test_read_param_logical, "test_read_param_logical")

  call suite%add(test_read_param_logical_missing, &
      "test_read_param_logical_missing", fatal=.true.)

  call suite%add(test_read_param_char_no_delim, &
      "test_read_param_char_no_delim")

  call suite%add(test_read_param_char_quote_delim, &
      "test_read_param_char_quote_delim")

  call suite%add(test_read_param_char_apostrophe_delim, &
      "test_read_param_char_apostrophe_delim")

  call suite%add(test_read_param_char_missing, &
      "test_read_param_char_missing", fatal=.true.)

  call suite%add(test_read_param_char_array, "test_read_param_char_array")

  call suite%add(test_read_param_char_array_missing, &
      "test_read_param_char_array_missing", fatal=.true.)

  call suite%add(test_read_param_time_date, "test_read_param_time_date")

  call suite%add(test_read_param_time_date_bad_format, &
      "test_read_param_time_date_bad_format", fatal=.true.)

  call suite%add(test_read_param_time_tuple, "test_read_param_time_tuple")

  call suite%add(test_read_param_time_bad_tuple, &
      "test_read_param_time_bad_tuple", fatal=.true.)

  call suite%add(test_read_param_time_bad_tuple_values, &
      "test_read_param_time_bad_tuple_values", fatal=.true.)

  call suite%add(test_read_param_time_missing, &
      "test_read_param_time_missing", fatal=.true.)

  call suite%add(test_read_param_time_undefined, &
      "test_read_param_time_undefined", fatal=.true.)

  call suite%add(test_read_param_time_type_err, &
      "test_read_param_time_type_err", fatal=.true.)

  call suite%add(test_read_param_time_unit, "test_read_param_time_unit")

  call suite%add(test_read_param_unused_fatal, &
    "test_read_param_unused_fatal", fatal=.true.)

  call suite%add(test_read_param_multiline_comment, &
      "test_read_param_multiline_comment")

  call suite%add(test_read_param_multiline_comment_unclosed, &
      "test_read_param_multiline_comment_unclosed", fatal=.true.)

  call suite%add(test_read_param_multiline_param, &
      "test_read_param_multiline_param")

  call suite%add(test_read_param_multiline_param_unclosed, &
      "test_read_param_multiline_param_unclosed", fatal=.true.)

  call suite%add(test_read_param_replace_tabs, "test_read_param_replace_tabs")

  call suite%add(test_read_param_pad_equals, "test_read_param_pad_equals")

  call suite%add(test_read_param_misplaced_quote, &
      "test_read_param_misplaced_quote", fatal=.true.)

  call suite%add(test_read_param_define, "test_read_param_define")

  call suite%add(test_read_param_define_as_flag, &
      "test_read_param_define_as_flag")

  call suite%add(test_read_param_override, "test_read_param_override")

  call suite%add(test_read_param_override_misplaced, &
      "test_read_param_override_misplaced", fatal=.true.)

  call suite%add(test_read_param_override_twice, &
      "test_read_param_override_twice", fatal=.true.)

  call suite%add(test_read_param_override_repeat, &
      "test_read_param_override_repeat", fatal=.true.)

  call suite%add(test_read_param_override_warn_chain, &
      "test_read_param_override_warn_chain")

  call suite%add(test_read_param_override_no_def, &
      "test_read_param_override_no_def", fatal=.true.)

  call suite%add(test_read_param_assign_after_override, &
      "test_read_param_assign_after_override")

  call suite%add(test_read_param_assign_twice, &
      "test_read_param_assign_twice", fatal=.true.)

  call suite%add(test_read_param_assign_repeat, &
      "test_read_param_assign_repeat")

  call suite%add(test_read_param_null_stmt, "test_read_param_null_stmt", &
      fatal=.true.)

  call suite%add(test_read_param_assign_in_define, &
      "test_read_param_assign_in_define", fatal=.true.)

  call suite%add(test_read_param_block, "test_read_param_block")

  ! FIXME: Test does not pass
  !call suite%add(test_read_param_block_stack, "test_read_param_block_stack")

  call suite%add(test_read_param_block_inline_stack, &
      "test_read_param_block_inline_stack")

  call suite%add(test_read_param_block_empty_pop, &
      "test_read_param_block_empty_pop", fatal=.true.)

  call suite%add(test_read_param_block_close_unopened, &
      "test_read_param_block_close_unopened", fatal=.true.)

  call suite%add(test_read_param_block_close_unnamed, &
      "test_read_param_block_close_unnamed", fatal=.true.)

  call suite%add(test_read_param_block_unmatched, &
      "test_read_param_block_unmatched", fatal=.true.)

  call suite%add(test_read_param_block_outside_block, &
      "test_read_param_block_outside_block")

  call suite%add(test_open_unallocated_block, "test_open_unallocated_block", &
      fatal=.true.)

  call suite%add(test_close_unallocated_block, &
      "test_close_unallocated_block", fatal=.true.)

  call suite%add(test_clear_unallocated_block, &
      "test_clear_unallocated_block", fatal=.true.)

  call suite%add(test_log_version_cs, "test_log_version_cs")

  call suite%add(test_log_version_plain, "test_log_version_plain")

  call suite%add(test_log_param_int, "test_log_param_int")

  call suite%add(test_log_param_int_array, "test_log_param_int_array")

  call suite%add(test_log_param_real, "test_log_param_real")

  call suite%add(test_log_param_real_array, "test_log_param_real_array")

  call suite%add(test_log_param_time, "test_log_param_time")

  call suite%add(test_log_param_time_as_date, "test_log_param_time_as_date")

  call suite%add(test_log_param_time_as_date_default, &
      "test_log_param_time_as_date_default")

  call suite%add(test_log_param_time_as_date_tick, &
      "test_log_param_time_as_date_tick")

  call suite%add(test_log_param_time_with_unit, &
      "test_log_param_time_with_unit")

  call suite%add(test_log_param_time_with_timeunit, &
      "test_log_param_time_with_timeunit")

  call suite%add(test_get_param_int, "test_get_param_int")

  call suite%add(test_get_param_int_no_read_no_log, &
      "test_get_param_int_no_read_no_log")

  call suite%add(test_get_param_int_array, "test_get_param_int_array")

  call suite%add(test_get_param_int_array_no_read_no_log, &
      "test_get_param_int_array_no_read_no_log")

  call suite%add(test_get_param_real, "test_get_param_real")

  call suite%add(test_get_param_real_no_read_no_log, &
      "test_get_param_real_n_read_no_log")

  call suite%add(test_get_param_real_array, "test_get_param_real_array")

  call suite%add(test_get_param_real_array_no_read_no_log, &
      "test_get_param_real_array_no_read_no_log")

  call suite%add(test_get_param_char, "test_get_param_char")

  call suite%add(test_get_param_char_no_read_no_log, &
      "test_get_param_char_no_read_no_log")

  call suite%add(test_get_param_char_array, "test_get_param_char_array")

  call suite%add(test_get_param_logical, "test_get_param_logical")

  call suite%add(test_get_param_logical_default, &
      "test_get_param_logical_default")

  call suite%add(test_get_param_logical_no_read_no_log, &
      "test_get_param_logical_no_read_no_log")

  call suite%add(test_get_param_time, "test_get_param_time")

  call suite%add(test_get_param_time_no_read_no_log, &
      "test_get_param_time_np_read_no_log")

  call suite%run()
end subroutine run_file_parser_tests

end module MOM_file_parser_tests
