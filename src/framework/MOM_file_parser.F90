module MOM_file_parser
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

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg and Alistair Adcroft, updated 9/2013.           *
!*                                                                     *
!*    The subroutines here parse a set of input files for the value    *
!*  a named parameter and sets that parameter at run time.  Currently  *
!*  these files use use one of several formats:                        *
!*    #define VAR       ! To set the logical VAR to true.              *
!*    VAR = True        ! To set the logical VAR to true.              *
!*    #undef VAR        ! To set the logical VAR to false.             *
!*    VAR = False       ! To set the logical VAR to false.             *
!*    #define VAR 999   ! To set the real or integer VAR to 999.       *
!*    VAR = 999         ! To set the real or integer VAR to 999.       *
!*    #override VAR = 888 ! To override a previously set value.        *
!*    VAR = 1.1, 2.2, 3.3 ! To set an array of real values.            *
!*                                                                     *
!*  In addition, when set by the get_param interface, the values of    *
!*  parameters are automatically logged, along with defaults, units,   *
!*  and a description.  It is an error for a variable to be overridden *
!*  more than once, and MOM6 has a facility to check for unused lines  *
!*  to set variables, which may indicate miss-spelled or archaic       *
!*  parameters.  Parameter names are case-specific, and lines may use  *
!*  a F90 or C++ style comment, starting with ! or //.                 *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_coms, only : root_PE, broadcast
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_error_handler, only : is_root_pe, stdlog, stdout
use MOM_time_manager, only : set_time, get_time, time_type, get_ticks_per_second
use MOM_document, only : doc_param, doc_module, doc_init, doc_end, doc_type
use MOM_document, only : doc_openBlock, doc_closeBlock
use MOM_string_functions, only : left_int, left_ints, slasher
use MOM_string_functions, only : left_real, left_reals

implicit none ; private

integer, parameter, public :: MAX_PARAM_FILES = 5 ! Maximum number of parameter files.
integer, parameter :: INPUT_STR_LENGTH = 200 ! Maximum linelength in parameter file.
integer, parameter :: FILENAME_LENGTH = 200  ! Maximum number of characters in
                                             ! file names.

! The all_PEs_read option should be eliminated with post-riga shared code.
logical :: all_PEs_read = .false.

! Defaults
logical, parameter :: report_unused_default = .false.
logical, parameter :: unused_params_fatal_default = .false.
logical, parameter :: log_to_stdout_default = .false.
logical, parameter :: complete_doc_default = .true.
logical, parameter :: minimal_doc_default = .true.

type, private :: file_data_type ; private
  integer :: num_lines = 0
  character(len=INPUT_STR_LENGTH), pointer, dimension(:) :: line => NULL()
  logical,                         pointer, dimension(:) :: line_used => NULL()
end type file_data_type

type :: link_parameter ; private
  type(link_parameter), pointer :: next => NULL() ! Facilitates linked list
  character(len=80) :: name                       ! Parameter name
  logical :: hasIssuedOverrideWarning = .false.   ! Has a default value
end type link_parameter

type :: parameter_block ; private
  character(len=240) :: name = ''                 ! Parameter name
end type parameter_block

type, public :: param_file_type ; private
  integer  :: nfiles = 0            ! The number of open files.
  integer  :: iounit(MAX_PARAM_FILES)   ! The unit number of an open file.
  character(len=FILENAME_LENGTH) :: filename(MAX_PARAM_FILES) ! The names of the open files.
  logical  :: NetCDF_file(MAX_PARAM_FILES)! If true, the input file is in NetCDF.
                                    ! This is not yet implemented.
  type(file_data_type) :: param_data(MAX_PARAM_FILES) ! Structures that contain
                                    ! the valid data lines from the parameter
                                    ! files, enabling all subsequent reads of
                                    ! parameter data to occur internally.
  logical  :: report_unused = report_unused_default ! If true, report any
                                    ! parameter lines that are not used in the run.
  logical  :: unused_params_fatal = unused_params_fatal_default  ! If true, kill
                                    ! the run if there are any unused parameters.
  logical  :: log_to_stdout = log_to_stdout_default ! If true, all log
                                    ! messages are also sent to stdout.
  logical  :: log_open = .false.    ! True if the log file has been opened.
  integer  :: stdout, stdlog        ! The units from stdout() and stdlog().
  character(len=240) :: doc_file    ! A file where all run-time parameters, their
                                    ! settings and defaults are documented.
  logical  :: complete_doc = complete_doc_default ! If true, document all
                                    ! run-time parameters.
  logical  :: minimal_doc = minimal_doc_default ! If true, document only those
                                    ! run-time parameters that differ from defaults.
  type(doc_type), pointer :: doc => NULL() ! A structure that contains information
                                    ! related to parameter documentation.
  type(link_parameter), pointer :: chain => NULL() ! Facilitates linked list
  type(parameter_block), pointer :: blockName => NULL() ! Name of active parameter block
end type param_file_type

public read_param, open_param_file, close_param_file, log_param, log_version
public doc_param, get_param
public clearParameterBlock, openParameterBlock, closeParameterBlock

interface read_param
  module procedure read_param_int, read_param_real, read_param_logical, &
                   read_param_char, read_param_char_array, read_param_time, &
                   read_param_int_array, read_param_real_array
end interface
interface log_param
  module procedure log_param_int, log_param_real, log_param_logical, &
                   log_param_char, log_param_time, &
                   log_param_int_array, log_param_real_array
end interface
interface get_param
  module procedure get_param_int, get_param_real, get_param_logical, &
                   get_param_char, get_param_char_array, get_param_time, &
                   get_param_int_array, get_param_real_array
end interface
interface log_version
  module procedure log_version_cs, log_version_plain
end interface

contains

subroutine open_param_file(filename, CS, checkable, component, doc_file_dir)
  character(len=*),           intent(in) :: filename
  type(param_file_type),   intent(inout) :: CS
  logical,          optional, intent(in) :: checkable
  character(len=*), optional, intent(in) :: component
  character(len=*), optional, intent(in) :: doc_file_dir

  logical :: file_exists, unit_in_use, Netcdf_file, may_check
  integer :: ios, iounit, strlen, i
  character(len=240) :: doc_path
  type(parameter_block), pointer :: block

  may_check = .true. ; if (present(checkable)) may_check = checkable

  ! Check for non-blank filename
  strlen = len_trim(filename)
  if (strlen == 0) then
    call MOM_error(FATAL, "open_param_file: Input file has not been specified.")
  endif

  ! Check that this file has not already been opened
  if (CS%nfiles > 0) then
    inquire(file=trim(filename), number=iounit)
    if (iounit /= -1) then
      do i = 1, CS%nfiles
        if (CS%iounit(i) == iounit) then
          if (trim(CS%filename(1)) /= trim(filename)) then
            call MOM_error(FATAL, &
              "open_param_file: internal inconsistency! "//trim(filename)// &
              " is registered as open but has the wrong unit number!")
          else
            call MOM_error(WARNING, &
              "open_param_file: file "//trim(filename)// &
              " has already been opened. This should NOT happen!"// &
              " Did you specify the same file twice in a namelist?")
            return
          endif ! filenames
        endif ! unit numbers
      enddo ! i
    endif
  endif

  ! Check that the file exists to readstdlog
  inquire(file=trim(filename), exist=file_exists)
  if (.not.file_exists) call MOM_error(FATAL, &
      "open_param_file: Input file "// trim(filename)//" does not exist.")

  Netcdf_file = .false.
  if (strlen > 3) then
    if (filename(strlen-2:strlen) == ".nc") Netcdf_file = .true.
  endif

  if (Netcdf_file) &
    call MOM_error(FATAL,"open_param_file: NetCDF files are not yet supported.")

  if (all_PEs_read .or. is_root_pe()) then
    ! Find an unused unit number.
    do iounit=10,512
      INQUIRE(iounit,OPENED=unit_in_use) ; if (.not.unit_in_use) exit
    enddo
    if (iounit >= 512) call MOM_error(FATAL, &
        "open_param_file: No unused file unit could be found.")

    ! Open the parameter file.
    open(iounit, file=trim(filename), access='SEQUENTIAL', &
         form='FORMATTED', action='READ', position='REWIND', iostat=ios)
    if (ios /= 0) call MOM_error(FATAL, "open_param_file: Error opening "// &
                                       trim(filename))
  else
    iounit = 1
  endif

  ! Store/register the unit and details
  i = CS%nfiles + 1
  CS%nfiles = i
  CS%iounit(i) = iounit
  CS%filename(i) = filename
  CS%NetCDF_file(i) = Netcdf_file
  allocate(block) ; block%name = '' ; CS%blockName => block

  call MOM_mesg("open_param_file: "// trim(filename)// &
                 " has been opened successfully.", 5)

  call populate_param_data(iounit, filename, CS%param_data(i))

  call read_param(CS,"SEND_LOG_TO_STDOUT",CS%log_to_stdout)
  call read_param(CS,"REPORT_UNUSED_PARAMS",CS%report_unused)
  call read_param(CS,"FATAL_UNUSED_PARAMS",CS%unused_params_fatal)
  CS%doc_file = "MOM_parameter_doc"
  if (present(component)) CS%doc_file = trim(component)//"_parameter_doc"
  call read_param(CS,"DOCUMENT_FILE", CS%doc_file)
  if (.not.may_check) then
    CS%report_unused = .false.
    CS%unused_params_fatal = .false.
  endif

  ! Open the log file.
  CS%stdlog = stdlog() ; CS%stdout = stdout()
  CS%log_open = (stdlog() > 0)

  doc_path = CS%doc_file
  if (len_trim(CS%doc_file) > 0) then
    CS%complete_doc = complete_doc_default
    call read_param(CS, "COMPLETE_DOCUMENTATION", CS%complete_doc)
    CS%minimal_doc = minimal_doc_default
    call read_param(CS, "MINIMAL_DOCUMENTATION", CS%minimal_doc)
    if (present(doc_file_dir)) then ; if (len_trim(doc_file_dir) > 0) then
      doc_path = trim(slasher(doc_file_dir))//trim(CS%doc_file)
    endif ; endif
  else
    CS%complete_doc = .false.
    CS%minimal_doc = .false.
  endif
  call doc_init(doc_path, CS%doc, CS%minimal_doc, CS%complete_doc)

end subroutine open_param_file

subroutine close_param_file(CS, quiet_close, component)
  type(param_file_type),   intent(inout) :: CS
  logical,          optional, intent(in) :: quiet_close
  character(len=*), optional, intent(in) :: component
! Arguments: CS - the param_file_type to close
!  (in,opt)  quiet_close - if present and true, do not do any logging with this
!                          call.
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=128) :: docfile_default
  character(len=40)  :: mod   ! This module's name.
  integer :: i, n, num_unused

  if (present(quiet_close)) then ; if (quiet_close) then
    do i = 1, CS%nfiles
      CS%iounit(i) = -1
      CS%filename(i) = ''
      CS%NetCDF_file(i) = .false.
      deallocate (CS%param_data(i)%line)
      deallocate (CS%param_data(i)%line_used)
    enddo
    CS%log_open = .false.
    call doc_end(CS%doc)
    return
  endif ; endif

  ! Log the parameters for the parser.
  mod = "MOM_file_parser"
  call log_version(CS, mod, version, "")
  call log_param(CS, mod, "SEND_LOG_TO_STDOUT", &
                        CS%log_to_stdout, &
                 "If true, all log messages are also sent to stdout.", &
                 default=log_to_stdout_default)
  call log_param(CS, mod, "REPORT_UNUSED_PARAMS", &
                        CS%report_unused, &
                 "If true, report any parameter lines that are not used \n"//&
                 "in the run.", default=report_unused_default)
  call log_param(CS, mod, "FATAL_UNUSED_PARAMS", &
                        CS%unused_params_fatal, &
                 "If true, kill the run if there are any unused \n"//&
                 "parameters.", default=unused_params_fatal_default)
  docfile_default = "MOM_parameter_doc"
  if (present(component)) docfile_default = trim(component)//"_parameter_doc"
  call log_param(CS, mod, "DOCUMENT_FILE", CS%doc_file, &
                 "The basename for files where run-time parameters, their\n"//&
                 "settings, units and defaults are documented. Blank will\n"//&
                 "disable all parameter documentation.", default=docfile_default)
  if (len_trim(CS%doc_file) > 0) then
    call log_param(CS, mod, "COMPLETE_DOCUMENTATION", &
                   CS%complete_doc, &
                  "If true, all run-time parameters are\n"//&
                  "documented in "//trim(CS%doc_file)//&
                  ".all .", default=complete_doc_default)
    call log_param(CS, mod, "MINIMAL_DOCUMENTATION", &
                   CS%minimal_doc, &
                  "If true, non-default run-time parameters are\n"//&
                  "documented in "//trim(CS%doc_file)//&
                  ".short .", default=minimal_doc_default)
  endif

  num_unused = 0
  do i = 1, CS%nfiles
    ! only root pe has the file open
    if (all_PEs_read .or. is_root_pe()) close(CS%iounit(i))
    call MOM_mesg("close_param_file: "// trim(CS%filename(i))// &
                 " has been closed successfully.", 5)

    ! Check for unused lines.
    if (is_root_pe() .and. (CS%report_unused .or. &
                            CS%unused_params_fatal)) then
      do n=1,CS%param_data(i)%num_lines
        if (.not.CS%param_data(i)%line_used(n)) then
          num_unused = num_unused + 1
          if (CS%report_unused) &
            call MOM_error(WARNING, "Unused line in "//trim(CS%filename(i))// &
                            " : "//trim(CS%param_data(i)%line(n)))
        endif
      enddo
    endif

    CS%iounit(i) = -1
    CS%filename(i) = ''
    CS%NetCDF_file(i) = .false.
    deallocate (CS%param_data(i)%line)
    deallocate (CS%param_data(i)%line_used)
  enddo

  if (is_root_pe() .and. (num_unused>0) .and. CS%unused_params_fatal) &
    call MOM_error(FATAL, "Run stopped because of unused parameter lines.")

  CS%log_open = .false.
  call doc_end(CS%doc)

end subroutine close_param_file

subroutine populate_param_data(iounit, filename, param_data)
  integer,                 intent(in) :: iounit
  character(len=*),        intent(in) :: filename
  type(file_data_type), intent(inout) :: param_data

  character(len=INPUT_STR_LENGTH) :: line
  integer :: num_lines
  logical :: inMultiLineComment

! Find the number of keyword lines in a parameter file
! Allocate the space to hold the lines in param_data%line
! Populate param_data%line with the keyword lines from parameter file

  if (iounit <= 0) return

  if (all_PEs_read .or. is_root_pe()) then
    ! rewind the parameter file
    rewind(iounit)

    ! count the number of valid entries in the parameter file
    num_lines = 0
    inMultiLineComment = .false.
    do while(.true.)
      read(iounit, '(a)', end=8, err=9) line
      line = replaceTabs(line)
      if (inMultiLineComment) then
        if (closeMultiLineComment(line)) inMultiLineComment=.false.
      else
        if (lastNonCommentNonBlank(line)>0) num_lines = num_lines + 1
        if (openMultiLineComment(line)) inMultiLineComment=.true.
      endif
    enddo ! while (.true.)
 8  continue ! get here when read() reaches EOF

    if (inMultiLineComment .and. is_root_pe()) &
      call MOM_error(FATAL, 'MOM_file_parser : A C-style multi-line comment '// &
                      '(/* ... */) was not closed before the end of '//trim(filename))

    ! allocate space to hold contents of the parameter file
    param_data%num_lines = num_lines
  endif  ! (is_root_pe())

  ! Broadcast the number of valid entries in parameter file
  if (.not. all_PEs_read) then
    call broadcast(param_data%num_lines, root_pe())
  endif

  ! Set up the space for storing the actual lines.
  num_lines = param_data%num_lines
  allocate (param_data%line(num_lines))
  allocate (param_data%line_used(num_lines))
  param_data%line(:) = ' '
  param_data%line_used(:) = .false.

  ! Read the actual lines.
  if (all_PEs_read .or. is_root_pe()) then
    ! rewind the parameter file
    rewind(iounit)

    ! Populate param_data%line
    num_lines = 0
    do while(.true.)
      read(iounit, '(a)', end=18, err=9) line
      line = replaceTabs(line)
      if (inMultiLineComment) then
        if (closeMultiLineComment(line)) inMultiLineComment=.false.
      else
        if (lastNonCommentNonBlank(line)>0) then
          line = removeComments(line)
          line = simplifyWhiteSpace(line(:len_trim(line)))
          num_lines = num_lines + 1
          param_data%line(num_lines) = line
        endif
        if (openMultiLineComment(line)) inMultiLineComment=.true.
      endif
    enddo ! while (.true.)
18  continue ! get here when read() reaches EOF

    if (num_lines /= param_data%num_lines) &
      call MOM_error(FATAL, 'MOM_file_parser : Found different number of '// &
                      'valid lines on second reading of '//trim(filename))
  endif  ! (is_root_pe())

  ! Broadcast the populated array param_data%line
  if (.not. all_PEs_read) then
    call broadcast(param_data%line, INPUT_STR_LENGTH, root_pe())
  endif

  return

9 call MOM_error(FATAL, "MOM_file_parser : "//&
                  "Error while reading file "//trim(filename))

end subroutine populate_param_data

function openMultiLineComment(string)
  character(len=*), intent(in) :: string
  logical                      :: openMultiLineComment
! True if a /* appears on this line without a closing */
  integer :: icom, last
  openMultiLineComment = .false.
  last = lastNonCommentIndex(string)+1
  icom = index(string(last:), "/*")
  if (icom > 0) then
    openMultiLineComment=.true.
    last = last+icom+1
  endif
  icom = index(string(last:), "*/") ; if (icom > 0) openMultiLineComment=.false.
end function openMultiLineComment

function closeMultiLineComment(string)
  character(len=*), intent(in) :: string
  logical                      :: closeMultiLineComment
! True if a */ appears on this line
  closeMultiLineComment = .false.
  if (index(string, "*/")>0) closeMultiLineComment=.true.
end function closeMultiLineComment

function lastNonCommentIndex(string)
  character(len=*), intent(in) :: string
  integer                      :: lastNonCommentIndex
! Find position of last character before any comments
! This s/r is the only place where a comment needs to be defined
  integer :: icom, last
  last = len_trim(string)
  icom = index(string(:last), "!") ; if (icom > 0) last = icom-1 ! F90 style
  icom = index(string(:last), "//") ; if (icom > 0) last = icom-1 ! C+ style
  icom = index(string(:last), "/*") ; if (icom > 0) last = icom-1 ! C style
  lastNonCommentIndex = last
end function lastNonCommentIndex

function lastNonCommentNonBlank(string)
  character(len=*), intent(in) :: string
  integer                      :: lastNonCommentNonBlank
! Find position of last non-blank character before any comments
  lastNonCommentNonBlank = len_trim(string(:lastNonCommentIndex(string))) ! Ignore remaining trailing blanks
end function lastNonCommentNonBlank

function replaceTabs(string)
  character(len=*), intent(in) :: string
  character(len=len(string))   :: replaceTabs
! Returns string with tabs replaced by a ablank
  integer :: i
  do i=1, len(string)
    if (string(i:i)==achar(9)) then
      replaceTabs(i:i)=" "
    else
      replaceTabs(i:i)=string(i:i)
    endif
  enddo
end function replaceTabs

function removeComments(string)
  character(len=*), intent(in) :: string
  character(len=len(string))   :: removeComments
! Trims comments and leading blanks from string
  integer :: last
  removeComments=repeat(" ",len(string))
  last = lastNonCommentNonBlank(string)
  removeComments(:last)=adjustl(string(:last)) ! Copy only the non-comment part of string
end function removeComments

function simplifyWhiteSpace(string)
  character(len=*), intent(in) :: string
  character(len=len(string)+16)   :: simplifyWhiteSpace
! Constructs a string with all repeated whitespace replaced with single blanks
! and insert white space where it helps delineate tokens (e.g. around =)
  integer :: i,j
  logical :: nonBlank = .false., insideString = .false.
  character(len=1) :: quoteChar=" "
  nonBlank  = .false.; insideString = .false. ! NOTE: For some reason this line is needed??
  i=0
  simplifyWhiteSpace=repeat(" ",len(string)+16)
  do j=1,len_trim(string)
    if (insideString) then ! Do not change formatting inside strings
      i=i+1
      simplifyWhiteSpace(i:i)=string(j:j)
      if (string(j:j)==quoteChar) insideString=.false. ! End of string
    else ! The following is outside of string delimiters
      if (string(j:j)==" " .or. string(j:j)==achar(9)) then ! Space or tab
        if (nonBlank) then ! Only copy a blank if the preceeding character was non-blank
          i=i+1
          simplifyWhiteSpace(i:i)=" " ! Not string(j:j) so that tabs are replace by blanks
          nonBlank=.false.
        endif
      elseif (string(j:j)=='"' .or. string(j:j)=="'") then ! Start a sting
        i=i+1
        simplifyWhiteSpace(i:i)=string(j:j)
        insideString=.true.
        quoteChar=string(j:j) ! Keep copy of starting quote
        nonBlank=.true.       ! For exit from string
      elseif (string(j:j)=='=') then
        ! Insert spaces if this character is "=" so that line contains " = "
        if (nonBlank) then
          i=i+1
          simplifyWhiteSpace(i:i)=" "
        endif
        i=i+2
        simplifyWhiteSpace(i-1:i)=string(j:j)//" "
        nonBlank=.false.
      else ! All other characters
        i=i+1
        simplifyWhiteSpace(i:i)=string(j:j)
        nonBlank=.true.
      endif
    endif ! if (insideString)
  enddo ! j
  if (insideString) then ! A missing close quote should be flagged
    if (is_root_pe()) call MOM_error(FATAL, &
      "There is a mismatched quote in the parameter file line: "// &
      trim(string))
  endif
end function simplifyWhiteSpace

subroutine read_param_int(CS, varname, value, fail_if_missing)
  type(param_file_type),  intent(in) :: CS
  character(len=*),       intent(in) :: varname
  integer,             intent(inout) :: value
  logical,      optional, intent(in) :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical            :: found, defined

  call get_variable_line(CS, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read(value_string(1),*,err = 1001) value
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call MOM_error(FATAL,'read_param_int: Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call MOM_error(FATAL,'read_param_int: Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
  endif
  return
 1001 call MOM_error(FATAL,'read_param_int: read error for integer variable '//trim(varname)// &
                             ' parsing "'//trim(value_string(1))//'"')
end subroutine read_param_int

subroutine read_param_int_array(CS, varname, value, fail_if_missing)
  type(param_file_type),  intent(in) :: CS
  character(len=*),       intent(in) :: varname
  integer,             intent(inout) :: value(:)
  logical,      optional, intent(in) :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical            :: found, defined

  call get_variable_line(CS, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read(value_string(1),*,end=991,err=1002) value
 991 return
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call MOM_error(FATAL,'read_param_int_array: Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call MOM_error(FATAL,'read_param_int_array: Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
  endif
  return
 1002 call MOM_error(FATAL,'read_param_int_array: read error for integer array '//trim(varname)// &
                             ' parsing "'//trim(value_string(1))//'"')
end subroutine read_param_int_array

subroutine read_param_real(CS, varname, value, fail_if_missing)
  type(param_file_type),  intent(in) :: CS
  character(len=*),       intent(in) :: varname
  real,                intent(inout) :: value
  logical,      optional, intent(in) :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical            :: found, defined

  call get_variable_line(CS, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read(value_string(1),*,err=1003) value
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call MOM_error(FATAL,'read_param_real: Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call MOM_error(FATAL,'read_param_real: Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
  endif
  return
 1003 call MOM_error(FATAL,'read_param_real: read error for real variable '//trim(varname)// &
                             ' parsing "'//trim(value_string(1))//'"')
end subroutine read_param_real

subroutine read_param_real_array(CS, varname, value, fail_if_missing)
  type(param_file_type),  intent(in) :: CS
  character(len=*),       intent(in) :: varname
  real,                intent(inout) :: value(:)
  logical,      optional, intent(in) :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical                         :: found, defined

  call get_variable_line(CS, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read(value_string(1),*,end=991,err=1004) value
 991 return
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call MOM_error(FATAL,'read_param_real_array: Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call MOM_error(FATAL,'read_param_real_array: Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
  endif
  return
 1004 call MOM_error(FATAL,'read_param_real_array: read error for real array '//trim(varname)// &
                             ' parsing "'//trim(value_string(1))//'"')
end subroutine read_param_real_array

subroutine read_param_char(CS, varname, value, fail_if_missing)
  type(param_file_type),  intent(in) :: CS
  character(len=*),       intent(in) :: varname
  character(len=*),    intent(inout) :: value
  logical,      optional, intent(in) :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical            :: found, defined

  call get_variable_line(CS, varname, found, defined, value_string)
  if (found) then
    value = trim(strip_quotes(value_string(1)))
  elseif (present(fail_if_missing)) then ; if (fail_if_missing) then
    call MOM_error(FATAL,'Unable to find variable '//trim(varname)// &
                         ' in any input files.')
  endif ; endif

end subroutine read_param_char

subroutine read_param_char_array(CS, varname, value, fail_if_missing)
  type(param_file_type),  intent(in) :: CS
  character(len=*),       intent(in) :: varname
  character(len=*),    intent(inout) :: value(:)
  logical,      optional, intent(in) :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=INPUT_STR_LENGTH) :: value_string(1), loc_string
  logical            :: found, defined
  integer            :: i, i_out

  call get_variable_line(CS, varname, found, defined, value_string)
  if (found) then
    loc_string = trim(value_string(1))
    i = index(loc_string,",")
    i_out = 1
    do while(i>0)
      value(i_out) = trim(strip_quotes(loc_string(:i-1)))
      i_out = i_out+1
      loc_string = trim(adjustl(loc_string(i+1:)))
      i = index(loc_string,",")
    enddo
    if (len_trim(loc_string)>0) then
      value(i_out) = trim(strip_quotes(adjustl(loc_string)))
      i_out = i_out+1
    endif
    do i=i_out,SIZE(value) ; value(i) = " " ; enddo
  elseif (present(fail_if_missing)) then ; if (fail_if_missing) then
    call MOM_error(FATAL,'Unable to find variable '//trim(varname)// &
                         ' in any input files.')
  endif ; endif

end subroutine read_param_char_array

subroutine read_param_logical(CS, varname, value, fail_if_missing)
  type(param_file_type),  intent(in) :: CS
  character(len=*),       intent(in) :: varname
  logical,             intent(inout) :: value
  logical,      optional, intent(in) :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical            :: found, defined

  call get_variable_line(CS, varname, found, defined, value_string, paramIsLogical=.true.)
  if (found) then
    value = defined
  elseif (present(fail_if_missing)) then ; if (fail_if_missing) then
    call MOM_error(FATAL,'Unable to find variable '//trim(varname)// &
                         ' in any input files.')
  endif ; endif
end subroutine read_param_logical


subroutine read_param_time(CS, varname, value, timeunit, fail_if_missing)
  type(param_file_type),  intent(in) :: CS
  character(len=*),       intent(in) :: varname
  type(time_type),     intent(inout) :: value
  real,         optional, intent(in) :: timeunit
  logical,      optional, intent(in) :: fail_if_missing
! This subroutine determines the value of an integer model parameter
! from a parameter file. The arguments are the unit of the open file
! which is to be read, the (case-sensitive) variable name, the variable
! where the value is to be stored, and (optionally) a flag indicating
! whether to fail if this parameter can not be found.  The unique argument
! to read time is the number of seconds to use as the unit of time being read.
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical            :: found, defined
  real               :: real_time, time_unit
  integer            :: days, secs

  call get_variable_line(CS, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    time_unit = 1.0 ; if (present(timeunit)) time_unit = timeunit
    read( value_string(1), *) real_time
    days = int(real_time*(time_unit/86400.0))
    secs = int(floor((real_time*(time_unit/86400.0)-days)*86400.0 + 0.5))
    value = set_time(secs, days)
  else
    if (present(fail_if_missing)) then ; if (fail_if_missing) then
      if (.not.found) then
        call MOM_error(FATAL,'Unable to find variable '//trim(varname)// &
                             ' in any input files.')
      else
        call MOM_error(FATAL,'Variable '//trim(varname)// &
                             ' found but not set in input files.')
      endif
    endif ; endif
  endif
end subroutine read_param_time

function strip_quotes(val_str)
  character(len=*) :: val_str
  character(len=INPUT_STR_LENGTH) :: strip_quotes
  ! Local variables
  integer :: i
  strip_quotes = val_str
  i = index(strip_quotes,ACHAR(34)) ! Double quote
  do while (i>0)
    if (i > 1) then ; strip_quotes = strip_quotes(:i-1)//strip_quotes(i+1:)
    else ; strip_quotes = strip_quotes(2:) ; endif
    i = index(strip_quotes,ACHAR(34)) ! Double quote
  enddo
  i = index(strip_quotes,ACHAR(39)) ! Single quote
  do while (i>0)
    if (i > 1) then ; strip_quotes = strip_quotes(:i-1)//strip_quotes(i+1:)
    else ; strip_quotes = strip_quotes(2:) ; endif
    i = index(strip_quotes,ACHAR(39)) ! Single quote
  enddo
end function strip_quotes

subroutine get_variable_line(CS, varname, found, defined, value_string, paramIsLogical)
  type(param_file_type),  intent(in) :: CS
  character(len=*),       intent(in) :: varname
  logical,               intent(out) :: found, defined
  character(len=*),      intent(out) :: value_string(:)
  logical, optional,      intent(in) :: paramIsLogical

  character(len=INPUT_STR_LENGTH) :: val_str, lname, origLine
  character(len=INPUT_STR_LENGTH) :: line, continuationBuffer, blockName
  character(len=FILENAME_LENGTH)  :: filename
  integer            :: is, id, isd, isu, ise, iso, verbose, ipf
  integer            :: last, last1, ival, oval, max_vals, count, contBufSize
  character(len=52)  :: set
  logical            :: found_override, found_equals
  logical            :: found_define, found_undef
  logical            :: force_cycle, defined_in_line, continuedLine
  logical            :: variableKindIsLogical, valueIsSame
  logical            :: inWrongBlock, fullPathParameter
  logical, parameter :: requireNamedClose = .false.
  set = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
  continuationBuffer = repeat(" ",INPUT_STR_LENGTH)
  contBufSize = 0
  verbose = 1

  variableKindIsLogical=.false.
  if (present(paramIsLogical)) variableKindIsLogical = paramIsLogical

  ! Find the first instance (if any) where the named variable is found, and
  ! return variables indicating whether this variable is defined and the string
  ! that contains the value of this variable.
  found = .false.
  oval = 0; ival = 0;
  max_vals = SIZE(value_string)
  do is=1,max_vals ; value_string(is) = " " ; enddo

  paramfile_loop: do ipf = 1, CS%nfiles
    filename = CS%filename(ipf)
    continuedLine = .false.
    blockName = ''

    ! Scan through each line of the file
    do count = 1, CS%param_data(ipf)%num_lines
      line = CS%param_data(ipf)%line(count)
      last = len_trim(line)

      last1 = max(1,last)
      ! Check if line ends in continuation character (either & or \)
      ! Note achar(92) is a backslash
      if (line(last1:last1) == achar(92).or.line(last1:last1) == "&") then
        continuationBuffer(contBufSize+1:contBufSize+len_trim(line))=line(:last-1)
        contBufSize=contBufSize + len_trim(line)-1
        continuedLine = .true.
        if (count==CS%param_data(ipf)%num_lines .and. is_root_pe()) &
           call MOM_error(FATAL, "MOM_file_parser : the last line"// &
                 " of the file ends in a continuation character but"// &
                 " there are no more lines to read. "// &
                 " Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
        cycle ! cycle inorder to append the next line of the file
      elseif (continuedLine) then
        ! If we reached this point then this is the end of line continuation
        continuationBuffer(contBufSize+1:contBufSize+len_trim(line))=line(:last)
        line = continuationBuffer
        continuationBuffer=repeat(" ",INPUT_STR_LENGTH) ! Clear for next use
        contBufSize = 0
        continuedLine = .false.
        last = len_trim(line)
      endif

      origLine = trim(line) ! Keep original for error messages

      ! Check for '#override' at start of line
      found_override = .false.; found_define = .false.; found_undef = .false.
      iso = index(line(:last), "#override " )!; if (is > 0) found_override = .true.
      if (iso>1) call MOM_error(FATAL, "MOM_file_parser : #override was found "// &
                 " but was not the first keyword."// &
                 " Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
      if (iso==1) then
        found_override = .true.
        if (index(line(:last), "#override define ")==1) found_define = .true.
        if (index(line(:last), "#override undef ")==1) found_undef = .true.
        line = trim(adjustl(line(iso+10:last))); last = len_trim(line)
      endif

      ! Check for start of fortran namelist, ie. '&namelist'
      if (index(line(:last),'&')==1) then
        iso=index(line(:last),' ')
        if (iso>0) then ! possibly simething else on this line
          blockName = pushBlockLevel(blockName,line(2:iso-1))
          line=trim(adjustl(line(iso:last)))
          last=len_trim(line)
          if (last==0) cycle ! nothing else on this line
        else ! just the namelist on this line
          if (len_trim(blockName)>0) then
            blockName = trim(blockName) // '%' //trim(line(2:last))
          else
            blockName = trim(line(2:last))
          endif
          call flag_line_as_read(CS%param_data(ipf)%line_used,count)
          cycle
        endif
      endif

      ! Newer form of parameter block, block%, %block or block%param or
      iso=index(line(:last),'%')
      fullPathParameter = .false.
      if (iso==1) then ! % is first character means this is a close
        if (len_trim(blockName)==0 .and. is_root_pe()) call MOM_error(FATAL, &
            'get_variable_line: An extra close block was encountered. Line="'// &
            trim(line(:last))//'"' )
        if (last>1 .and. trim(blockName)/=trim(line(2:last)) .and. is_root_pe()) &
            call MOM_error(FATAL, 'get_variable_line: A named close for a parameter'// &
            ' block did not match the open block. Line="'//trim(line(:last))//'"' )
        if (last==1 .and. requireNamedClose) & ! line = '%' is a generic (unnamed) close
            call MOM_error(FATAL, 'get_variable_line: A named close for a parameter'// &
            ' block is required but found "%". Block="'//trim(blockName)//'"' )
        blockName = popBlockLevel(blockName)
        call flag_line_as_read(CS%param_data(ipf)%line_used,count)
      elseif (iso==last) then ! This is a new block if % is last character
        blockName = pushBlockLevel(blockName, line(:iso-1))
        call flag_line_as_read(CS%param_data(ipf)%line_used,count)
      else ! This is of the form block%parameter = ... (full path parameter)
        iso=index(line(:last),'%',.true.)
        ! Check that the parameter block names on the line matches the state set by the caller
        if (iso>0 .and. trim(CS%blockName%name)==trim(line(:iso-1))) then
          fullPathParameter = .true.
          line = trim(line(iso+1:last)) ! Strip away the block name for subsequent processing
          last = len_trim(line)
        endif
      endif

      ! We should only interpret this line if this block is the active block
      inWrongBlock = .false.
      if (len_trim(blockName)>0) then ! In a namelist block in file
        if (trim(CS%blockName%name)/=trim(blockName)) inWrongBlock = .true. ! Not in the required block
      endif
      if (len_trim(CS%blockName%name)>0) then ! In a namelist block in the model
        if (trim(CS%blockName%name)/=trim(blockName)) inWrongBlock = .true. ! Not in the required block
      endif

      ! Check for termination of a fortran namelist (with a '/')
      if (line(last:last)=='/') then
        if (len_trim(blockName)==0 .and. is_root_pe()) call MOM_error(FATAL, &
            'get_variable_line: An extra namelist/block end was encountered. Line="'// &
            trim(line(:last))//'"' )
        blockName = popBlockLevel(blockName)
        last = last - 1 ! Ignore the termination character from here on
      endif
      if (inWrongBlock .and. .not. fullPathParameter) then
        if (index(" "//line(:last+1), " "//trim(varname)//" ")>0) &
          call MOM_error(WARNING,"MOM_file_parser : "//trim(varname)// &
               ' found outside of block '//trim(CS%blockName%name)//'%. Ignoring.')
        cycle
      endif

      ! Determine whether this line mentions the named parameter or not
      if (index(" "//line(:last)//" ", " "//trim(varname)//" ") == 0) cycle

      ! Detect keywords
      found_equals = .false.
      isd = index(line(:last), "define" )!; if (isd > 0) found_define = .true.
      isu = index(line(:last), "undef" )!; if (isu > 0) found_undef = .true.
      ise = index(line(:last), " = " ); if (ise > 1) found_equals = .true.
      if (index(line(:last), "#define ")==1) found_define = .true.
      if (index(line(:last), "#undef ")==1) found_undef = .true.

      ! Check for missing, mutually exclusive or incomplete keywords
      if (is_root_pe()) then
        if (.not. (found_define .or. found_undef .or. found_equals)) &
               call MOM_error(FATAL, "MOM_file_parser : the parameter name '"// &
                 trim(varname)//"' was found without define or undef."// &
                 " Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
        if (found_define .and. found_undef) call MOM_error(FATAL, &
                 "MOM_file_parser : Both 'undef' and 'define' occur."// &
                 " Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
        if (found_equals .and. (found_define .or. found_undef)) &
               call MOM_error(FATAL, &
                 "MOM_file_parser : Both 'a=b' and 'undef/define' syntax occur."// &
                 " Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
        if (found_override .and. .not. (found_define .or. found_undef .or. found_equals)) &
               call MOM_error(FATAL, "MOM_file_parser : override was found "// &
                 " without a define or undef."// &
                 " Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
      endif

      ! Interpret the line and collect values, if any
      if (found_define) then
        ! Move starting pointer to first letter of defined name.
        is = isd + 5 + scan(line(isd+6:last), set)

        id = scan(line(is:last), ' ')  ! Find space between name and value
        if ( id == 0 ) then
          ! There is no space so the name is simply being defined.
          lname = trim(line(is:last))
          if (trim(lname) /= trim(varname)) cycle
          val_str = " "
        else
          ! There is a string or number after the name.
          lname = trim(line(is:is+id-1))
          if (trim(lname) /= trim(varname)) cycle
          val_str = trim(adjustl(line(is+id:last)))
        endif
        found = .true. ; defined_in_line = .true.
      elseif (found_undef) then
        ! Move starting pointer to first letter of undefined name.
        is = isu + 4 + scan(line(isu+5:last), set)

        id = scan(line(is:last), ' ')  ! Find the first space after the name.
        if (id > 0) last = is + id - 1
        lname = trim(line(is:last))
        if (trim(lname) /= trim(varname)) cycle
        val_str = " "
        found = .true. ; defined_in_line = .false.
      elseif (found_equals) then
        ! Move starting pointer to first letter of defined name.
        is = scan(line(1:ise), set)
        lname = trim(line(is:ise-1))
        if (trim(lname) /= trim(varname)) cycle
        val_str = trim(adjustl(line(ise+3:last)))
        if (variableKindIsLogical) then ! Special handling for logicals
          read(val_str(:len_trim(val_str)),*) defined_in_line
        else
          defined_in_line = .true.
        endif
        found = .true.
      else
        call MOM_error(FATAL, "MOM_file_parser (non-root PE?): the parameter name '"// &
           trim(varname)//"' was found without an assignment, define or undef."// &
           " Line: '"//trim(line(:last))//"'"//" in file "//trim(filename)//".")
      endif

      ! This line has now been used.
      call flag_line_as_read(CS%param_data(ipf)%line_used,count)

      ! Detect inconsistencies
      force_cycle = .false.
      valueIsSame = (trim(val_str) == trim(value_string(max_vals)))
      if (found_override .and. (oval >= max_vals)) then
        if (is_root_pe()) then
          if ((defined_in_line .neqv. defined) .or. .not. valueIsSame) then
            call MOM_error(FATAL,"MOM_file_parser : "//trim(varname)// &
                     " found with multiple inconsistent overrides."// &
                     " Line A: '"//trim(value_string(max_vals))//"'"//&
                     " Line B: '"//trim(line(:last))//"'"//&
                     " in file "//trim(filename)//" caused the model failure.")
          else
            call MOM_error(WARNING,"MOM_file_parser : "//trim(varname)// &
                     " over-ridden more times than is permitted."// &
                     " Line: '"//trim(line(:last))//"'"//&
                     " in file "//trim(filename)//" is being ignored.")
          endif
        endif
        force_cycle = .true.
      endif
      if (.not.found_override .and. (oval > 0)) then
        if (is_root_pe()) &
          call MOM_error(WARNING,"MOM_file_parser : "//trim(varname)// &
                   " has already been over-ridden."// &
                   " Line: '"//trim(line(:last))//"'"//&
                   " in file "//trim(filename)//" is being ignored.")
        force_cycle = .true.
      endif
      if (.not.found_override .and. (ival >= max_vals)) then
        if (is_root_pe()) then
          if ((defined_in_line .neqv. defined) .or. .not. valueIsSame) then
            call MOM_error(FATAL,"MOM_file_parser : "//trim(varname)// &
                     " found with multiple inconsistent definitions."// &
                     " Line A: '"//trim(value_string(max_vals))//"'"//&
                     " Line B: '"//trim(line(:last))//"'"//&
                     " in file "//trim(filename)//" caused the model failure.")
          else
            call MOM_error(WARNING,"MOM_file_parser : "//trim(varname)// &
                     " occurs more times than is permitted."// &
                     " Line: '"//trim(line(:last))//"'"//&
                     " in file "//trim(filename)//" is being ignored.")
          endif
        endif
        force_cycle = .true.
      endif
      if (force_cycle) cycle

      ! Store new values
      if (found_override) then
        oval = oval + 1
        value_string(oval) = trim(val_str)
        defined = defined_in_line
        if (verbose > 0 .and. ival > 0 .and. is_root_pe() .and. &
            .not. overrideWarningHasBeenIssued(CS%chain, trim(varname)) ) &
          call MOM_error(WARNING,"MOM_file_parser : "//trim(varname)// &
                 " over-ridden.  Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
      else ! (.not. found_overide)
        ival = ival + 1
        value_string(ival) = trim(val_str)
        defined = defined_in_line
        if (verbose > 1 .and. is_root_pe()) &
          call MOM_error(WARNING,"MOM_file_parser : "//trim(varname)// &
                 " set.  Line: '"//trim(line(:last))//"'"//&
                 " in file "//trim(filename)//".")
      endif

    enddo ! CS%param_data(ipf)%num_lines

    if (len_trim(blockName)>0 .and. is_root_pe()) call MOM_error(FATAL, &
      'A namelist/parameter block was not closed. Last open block appears '// &
      'to be "'//trim(blockName)//'".')

  enddo paramfile_loop

end subroutine get_variable_line

subroutine flag_line_as_read(line_used,count)
  logical, dimension(:), pointer :: line_used
  integer,            intent(in) :: count
  line_used(count) = .true.
end subroutine flag_line_as_read

function overrideWarningHasBeenIssued(chain, varName)
  type(link_parameter), pointer    :: chain
  character(len=*),     intent(in) :: varName
  logical                          :: overrideWarningHasBeenIssued
! Returns true if an override warning has been issued for the variable varName
  type(link_parameter), pointer :: newLink, this
  overrideWarningHasBeenIssued = .false.
  this => chain
  do while( associated(this) )
    if (trim(varName) == trim(this%name)) then
      overrideWarningHasBeenIssued = .true.
      return
    endif
    this => this%next
  enddo
  allocate(newLink)
  newLink%name = trim(varName)
  newLink%hasIssuedOverrideWarning = .true.
  newLink%next => chain
  chain => newLink
end function overrideWarningHasBeenIssued

! The following subroutines write out to a log file.

!> Log the version of a module to a log file and/or stdout, and/or to the
!! parameter documentation file.
subroutine log_version_cs(CS, modulename, version, desc)
  type(param_file_type),      intent(in) :: CS         !< File parser type
  character(len=*),           intent(in) :: modulename !< Name of calling module
  character(len=*),           intent(in) :: version    !< Version string of module
  character(len=*), optional, intent(in) :: desc       !< Module description
  ! Local variables
  character(len=240) :: mesg

  mesg = trim(modulename)//": "//trim(version)
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  if (present(desc)) call doc_module(CS%doc, modulename, desc)

end subroutine log_version_cs

!> Log the version of a module to a log file and/or stdout.
subroutine log_version_plain(modulename, version)
  character(len=*),           intent(in) :: modulename !< Name of calling module
  character(len=*),           intent(in) :: version    !< Version string of module
  ! Local variables
  character(len=240) :: mesg

  mesg = trim(modulename)//": "//trim(version)
  if (is_root_pe()) then
    write(stdlog(),'(a)') trim(mesg)
  endif

end subroutine log_version_plain

subroutine log_param_int(CS, modulename, varname, value, desc, units, &
                         default, layoutParam)
  type(param_file_type),      intent(in) :: CS
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  integer,                    intent(in) :: value
  character(len=*), optional, intent(in) :: desc, units
  integer,          optional, intent(in) :: default
  logical,          optional, intent(in) :: layoutParam
! This subroutine writes the value of an integer parameter to a log file,
! along with its name and the module it came from.
  character(len=240) :: mesg, myunits

  write(mesg, '("  ",a," ",a,": ",a)') trim(modulename), trim(varname), trim(left_int(value))
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits=" "; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default, &
                   layoutParam=layoutParam)

end subroutine log_param_int

subroutine log_param_int_array(CS, modulename, varname, value, desc, &
                               units, default, layoutParam)
  type(param_file_type),      intent(in) :: CS
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  integer,                    intent(in) :: value(:)
  character(len=*), optional, intent(in) :: desc, units
  integer,          optional, intent(in) :: default
  logical,          optional, intent(in) :: layoutParam
! This subroutine writes the value of an integer parameter to a log file,
! along with its name and the module it came from.
  character(len=1320) :: mesg
  character(len=240) :: myunits

  write(mesg, '("  ",a," ",a,": ",A)') trim(modulename), trim(varname), trim(left_ints(value))
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits=" "; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default, &
                   layoutParam=layoutParam)

end subroutine log_param_int_array

subroutine log_param_real(CS, modulename, varname, value, desc, units, &
                          default)
  type(param_file_type),      intent(in) :: CS
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  real,                       intent(in) :: value
  character(len=*), optional, intent(in) :: desc, units
  real,             optional, intent(in) :: default
! This subroutine writes the value of a real parameter to a log file,
! along with its name and the module it came from.
  character(len=240) :: mesg, myunits

  write(mesg, '("  ",a," ",a,": ",a)') &
    trim(modulename), trim(varname), trim(left_real(value))
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits="not defined"; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default)

end subroutine log_param_real

subroutine log_param_real_array(CS, modulename, varname, value, desc, &
                                units, default)
  type(param_file_type),      intent(in) :: CS
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  real,                       intent(in) :: value(:)
  character(len=*), optional, intent(in) :: desc, units
  real,             optional, intent(in) :: default
! This subroutine writes the value of a real parameter to a log file,
! along with its name and the module it came from.
  character(len=1320) :: mesg
  character(len=240) :: myunits

 !write(mesg, '("  ",a," ",a,": ",ES19.12,99(",",ES19.12))') &
 !write(mesg, '("  ",a," ",a,": ",G,99(",",G))') &
 !  trim(modulename), trim(varname), value
  write(mesg, '("  ",a," ",a,": ",a)') &
    trim(modulename), trim(varname), trim(left_reals(value))
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits="not defined"; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default)

end subroutine log_param_real_array

subroutine log_param_logical(CS, modulename, varname, value, desc, &
                             units, default, layoutParam)
  type(param_file_type),      intent(in) :: CS
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  logical,                    intent(in) :: value
  character(len=*), optional, intent(in) :: desc, units
  logical,          optional, intent(in) :: default
  logical,          optional, intent(in) :: layoutParam
! This subroutine writes the value of a logical parameter to a log file,
! along with its name and the module it came from.
  character(len=240) :: mesg, myunits

  if (value) then
    write(mesg, '("  ",a," ",a,": True")') trim(modulename), trim(varname)
  else
    write(mesg, '("  ",a," ",a,": False")') trim(modulename), trim(varname)
  endif
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits="Boolean"; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default, &
                   layoutParam=layoutParam)

end subroutine log_param_logical

subroutine log_param_char(CS, modulename, varname, value, desc, units, &
                          default, layoutParam)
  type(param_file_type),      intent(in) :: CS
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  character(len=*),           intent(in) :: value
  character(len=*), optional, intent(in) :: desc, units
  character(len=*), optional, intent(in) :: default
  logical,          optional, intent(in) :: layoutParam
! This subroutine writes the value of a character string parameter to a log
! file, along with its name and the module it came from.
  character(len=240) :: mesg, myunits

  write(mesg, '("  ",a," ",a,": ",a)') &
    trim(modulename), trim(varname), trim(value)
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits=" "; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default, &
                   layoutParam=layoutParam)

end subroutine log_param_char

subroutine log_param_time(CS, modulename, varname, value, desc, units, &
                          default, timeunit)
  type(param_file_type),      intent(in) :: CS
  character(len=*),           intent(in) :: modulename
  character(len=*),           intent(in) :: varname
  type(time_type),            intent(in) :: value
  character(len=*), optional, intent(in) :: desc, units
  type(time_type),  optional, intent(in) :: default
  real,             optional, intent(in) :: timeunit
! This subroutine writes the value of a time-type parameter to a log file,
! along with its name and the module it came from.
  real :: real_time, real_default
  logical :: use_timeunit = .false.
  character(len=240) :: mesg, myunits
  integer :: days, secs, ticks

  call get_time(value, secs, days, ticks)

  if (ticks == 0) then
    write(mesg, '("  ",a," ",a," (Time): ",i0,":",i0)') trim(modulename), &
       trim(varname), days, secs
  else
    write(mesg, '("  ",a," ",a," (Time): ",i0,":",i0,":",i0)') trim(modulename), &
       trim(varname), days, secs, ticks
  endif
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  if (present(desc)) then
    if (present(timeunit)) use_timeunit = (timeunit > 0.0)
    if (use_timeunit) then
      if (present(units)) then
        write(myunits(1:240),'(A)') trim(units)
      else
        if (abs(timeunit-1.0) < 0.01) then ; myunits = "seconds"
        elseif (abs(timeunit-3600.0) < 1.0) then ; myunits = "hours"
        elseif (abs(timeunit-86400.0) < 1.0) then ; myunits = "days"
        elseif (abs(timeunit-3.1e7) < 1.0e6) then ; myunits = "years"
        else ; write(myunits,'(es8.2," sec")') timeunit ; endif
      endif
      real_time = (86400.0/timeunit)*days + secs/timeunit
      if (ticks > 0) real_time = real_time + &
                           real(ticks) / (timeunit*get_ticks_per_second())
      if (present(default)) then
        call get_time(default, secs, days, ticks)
        real_default = (86400.0/timeunit)*days + secs/timeunit
        if (ticks > 0) real_default = real_default + &
                           real(ticks) / (timeunit*get_ticks_per_second())
        call doc_param(CS%doc, varname, desc, myunits, real_time, real_default)
      else
        call doc_param(CS%doc, varname, desc, myunits, real_time)
      endif
    else
      myunits='not defined'; if (present(units)) write(myunits(1:240),'(A)') trim(units)
      call doc_param(CS%doc, varname, desc, myunits, value, default)
    endif
  endif

end subroutine log_param_time


subroutine get_param_int(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, &
               static_value, layoutParam)
  type(param_file_type),      intent(in)    :: CS
  character(len=*),           intent(in)    :: modulename
  character(len=*),           intent(in)    :: varname
  integer,                    intent(inout) :: value
  character(len=*), optional, intent(in)    :: desc, units
  integer,          optional, intent(in)    :: default, static_value
  logical,          optional, intent(in)    :: fail_if_missing
  logical,          optional, intent(in)    :: do_not_read, do_not_log
  logical,          optional, intent(in)    :: layoutParam
! This subroutine writes the value of a real parameter to a log file,
! along with its name and the module it came from.
  logical :: do_read, do_log

  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (do_read) then
    if (present(default)) value = default
    if (present(static_value)) value = static_value
    call read_param_int(CS, varname, value, fail_if_missing)
  endif

  if (do_log) then
    call log_param_int(CS, modulename, varname, value, desc, units, &
                       default, layoutParam)
  endif

end subroutine get_param_int

subroutine get_param_int_array(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, &
               static_value, layoutParam)
  type(param_file_type),      intent(in)    :: CS
  character(len=*),           intent(in)    :: modulename
  character(len=*),           intent(in)    :: varname
  integer,                    intent(inout) :: value(:)
  character(len=*), optional, intent(in)    :: desc, units
  integer,          optional, intent(in)    :: default, static_value
  logical,          optional, intent(in)    :: fail_if_missing
  logical,          optional, intent(in)    :: do_not_read, do_not_log
  logical,          optional, intent(in)    :: layoutParam
! This subroutine writes the value of a real parameter to a log file,
! along with its name and the module it came from.
  logical :: do_read, do_log

  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (do_read) then
    if (present(default)) then ; value(:) = default ; endif
    if (present(static_value)) then ; value(:) = static_value ; endif
    call read_param_int_array(CS, varname, value, fail_if_missing)
  endif

  if (do_log) then
    call log_param_int_array(CS, modulename, varname, value, desc, &
                             units, default, layoutParam)
  endif

end subroutine get_param_int_array

subroutine get_param_real(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, static_value)
  type(param_file_type),      intent(in)    :: CS
  character(len=*),           intent(in)    :: modulename
  character(len=*),           intent(in)    :: varname
  real,                       intent(inout) :: value
  character(len=*), optional, intent(in)    :: desc, units
  real,             optional, intent(in)    :: default, static_value
  logical,          optional, intent(in)    :: fail_if_missing
  logical,          optional, intent(in)    :: do_not_read, do_not_log
! This subroutine writes the value of a real parameter to a log file,
! along with its name and the module it came from.
  logical :: do_read, do_log

  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (do_read) then
    if (present(default)) value = default
    if (present(static_value)) value = static_value
    call read_param_real(CS, varname, value, fail_if_missing)
  endif

  if (do_log) then
    call log_param_real(CS, modulename, varname, value, desc, units, &
                        default)
  endif

end subroutine get_param_real

subroutine get_param_real_array(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, static_value)
  type(param_file_type),      intent(in)    :: CS
  character(len=*),           intent(in)    :: modulename
  character(len=*),           intent(in)    :: varname
  real,                       intent(inout) :: value(:)
  character(len=*), optional, intent(in)    :: desc, units
  real,             optional, intent(in)    :: default, static_value
  logical,          optional, intent(in)    :: fail_if_missing
  logical,          optional, intent(in)    :: do_not_read, do_not_log
! This subroutine writes the value of a real parameter to a log file,
! along with its name and the module it came from.
  logical :: do_read, do_log

  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (do_read) then
    if (present(default)) then ; value(:) = default ; endif
    if (present(static_value)) then ; value(:) = static_value ; endif
    call read_param_real_array(CS, varname, value, fail_if_missing)
  endif

  if (do_log) then
    call log_param_real_array(CS, modulename, varname, value, desc, &
                              units, default)
  endif

end subroutine get_param_real_array

subroutine get_param_char(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, &
               static_value, layoutParam)
  type(param_file_type),      intent(in)    :: CS
  character(len=*),           intent(in)    :: modulename
  character(len=*),           intent(in)    :: varname
  character(len=*),           intent(inout) :: value
  character(len=*), optional, intent(in)    :: desc, units
  character(len=*), optional, intent(in)    :: default, static_value
  logical,          optional, intent(in)    :: fail_if_missing
  logical,          optional, intent(in)    :: do_not_read, do_not_log
  logical,          optional, intent(in)    :: layoutParam
! This subroutine writes the value of a real parameter to a log file,
! along with its name and the module it came from.
  logical :: do_read, do_log

  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (do_read) then
    if (present(default)) value = default
    if (present(static_value)) value = static_value
    call read_param_char(CS, varname, value, fail_if_missing)
  endif

  if (do_log) then
    call log_param_char(CS, modulename, varname, value, desc, units, &
                        default, layoutParam)
  endif

end subroutine get_param_char

subroutine get_param_char_array(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, static_value)
  type(param_file_type),      intent(in)    :: CS
  character(len=*),           intent(in)    :: modulename
  character(len=*),           intent(in)    :: varname
  character(len=*),           intent(inout) :: value(:)
  character(len=*), optional, intent(in)    :: desc, units
  character(len=*), optional, intent(in)    :: default, static_value
  logical,          optional, intent(in)    :: fail_if_missing
  logical,          optional, intent(in)    :: do_not_read, do_not_log
! This subroutine writes the value of a real parameter to a log file,
! along with its name and the module it came from.
  logical :: do_read, do_log
  integer :: i, len_tot, len_val
  character(len=240) :: cat_val

  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (do_read) then
    if (present(default)) then ; value(:) = default ; endif
    if (present(static_value)) then ; value(:) = static_value ; endif
    call read_param_char_array(CS, varname, value, fail_if_missing)
  endif

  if (do_log) then
    cat_val = trim(value(1)); len_tot = len_trim(value(1))
    do i=2,size(value)
      len_val = len_trim(value(i))
      if ((len_val > 0) .and. (len_tot + len_val + 2 < 240)) then
        cat_val = trim(cat_val)//ACHAR(34)// ", "//ACHAR(34)//trim(value(i))
        len_tot = len_tot + len_val
      endif
    enddo
    call log_param_char(CS, modulename, varname, cat_val, desc, &
                        units, default)
  endif

end subroutine get_param_char_array

subroutine get_param_logical(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, &
               static_value, layoutParam)
  type(param_file_type),      intent(in)    :: CS
  character(len=*),           intent(in)    :: modulename
  character(len=*),           intent(in)    :: varname
  logical,                    intent(inout) :: value
  character(len=*), optional, intent(in)    :: desc, units
  logical,          optional, intent(in)    :: default, static_value
  logical,          optional, intent(in)    :: fail_if_missing
  logical,          optional, intent(in)    :: do_not_read, do_not_log
  logical,          optional, intent(in)    :: layoutParam
! This subroutine writes the value of a real parameter to a log file,
! along with its name and the module it came from.
  logical :: do_read, do_log

  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (do_read) then
    if (present(default)) value = default
    if (present(static_value)) value = static_value
    call read_param_logical(CS, varname, value, fail_if_missing)
  endif

  if (do_log) then
    call log_param_logical(CS, modulename, varname, value, desc, &
                           units, default, layoutParam)
  endif

end subroutine get_param_logical

subroutine get_param_time(CS, modulename, varname, value, desc, units, &
                          default, fail_if_missing, do_not_read, do_not_log, &
                          timeunit, static_value)
  type(param_file_type),      intent(in)    :: CS
  character(len=*),           intent(in)    :: modulename
  character(len=*),           intent(in)    :: varname
  type(time_type),            intent(inout) :: value
  character(len=*), optional, intent(in)    :: desc, units
  type(time_type),  optional, intent(in)    :: default, static_value
  logical,          optional, intent(in)    :: fail_if_missing
  logical,          optional, intent(in)    :: do_not_read, do_not_log
  real,             optional, intent(in)    :: timeunit
! This subroutine writes the value of a real parameter to a log file,
! along with its name and the module it came from.
  logical :: do_read, do_log

  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log

  if (do_read) then
    if (present(default)) value = default
    if (present(static_value)) value = static_value
    call read_param_time(CS, varname, value, timeunit, fail_if_missing)
  endif

  if (do_log) then
    call log_param_time(CS, modulename, varname, value, desc, &
                           units, default, timeunit)
  endif

end subroutine get_param_time

! -----------------------------------------------------------------------------

subroutine clearParameterBlock(CS)
  type(param_file_type), intent(in) :: CS
! Resets the parameter block name to blank
  type(parameter_block), pointer :: block
  if (associated(CS%blockName)) then
    block => CS%blockName
    block%name = ''
  else
    if (is_root_pe()) call MOM_error(FATAL, &
      'clearParameterBlock: A clear was attempted before allocation.')
  endif
end subroutine clearParameterBlock

subroutine openParameterBlock(CS,blockName,desc)
  type(param_file_type),      intent(in) :: CS
  character(len=*),           intent(in) :: blockName
  character(len=*), optional, intent(in) :: desc
! Tags blockName onto the end of the active parameter block name
  type(parameter_block), pointer :: block
  if (associated(CS%blockName)) then
    block => CS%blockName
    block%name = pushBlockLevel(block%name,blockName)
    call doc_openBlock(CS%doc,block%name,desc)
  else
    if (is_root_pe()) call MOM_error(FATAL, &
      'openParameterBlock: A push was attempted before allocation.')
  endif
end subroutine openParameterBlock

subroutine closeParameterBlock(CS)
  type(param_file_type), intent(in) :: CS
! Remove the lowest level of recursion from the active block name
  type(parameter_block), pointer :: block

  if (associated(CS%blockName)) then
    block => CS%blockName
    if (is_root_pe().and.len_trim(block%name)==0) call MOM_error(FATAL, &
      'closeParameterBlock: A pop was attempted on an empty stack. ("'//&
      trim(block%name)//'")')
    call doc_closeBlock(CS%doc,block%name)
  else
    if (is_root_pe()) call MOM_error(FATAL, &
      'closeParameterBlock: A pop was attempted before allocation.')
  endif
  block%name = popBlockLevel(block%name)
end subroutine closeParameterBlock

function pushBlockLevel(oldblockName,newBlockName)
  character(len=*),        intent(in) :: oldBlockName, newBlockName
  character(len=len(oldBlockName)+40) :: pushBlockLevel
! Extends block name (deeper level of parameter block)
  if (len_trim(oldBlockName)>0) then
    pushBlockLevel=trim(oldBlockName)//'%'//trim(newBlockName)
  else
    pushBlockLevel=trim(newBlockName)
  endif
end function pushBlockLevel

function popBlockLevel(oldblockName)
  character(len=*),        intent(in) :: oldBlockName
  character(len=len(oldBlockName)+40) :: popBlockLevel
! Truncates block name (shallower level of parameter block)
  integer :: i
  i = index(trim(oldBlockName), '%', .true.)
  if (i>1) then
    popBlockLevel = trim(oldBlockName(1:i-1))
  elseif (i==0) then
    popBlockLevel = ''
  else ! i==1
    if (is_root_pe()) call MOM_error(FATAL, &
      'popBlockLevel: A pop was attempted leaving an empty block name.')
  endif
end function popBlockLevel

end module MOM_file_parser
