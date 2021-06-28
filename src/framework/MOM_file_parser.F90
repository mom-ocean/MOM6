!> The MOM6 facility to parse input files for runtime parameters
module MOM_file_parser

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms, only : root_PE, broadcast
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_error_handler, only : is_root_pe, stdlog, stdout
use MOM_time_manager, only : get_time, time_type, get_ticks_per_second
use MOM_time_manager, only : set_date, get_date, real_to_time, operator(-), set_time
use MOM_document, only : doc_param, doc_module, doc_init, doc_end, doc_type
use MOM_document, only : doc_openBlock, doc_closeBlock
use MOM_string_functions, only : left_int, left_ints, slasher
use MOM_string_functions, only : left_real, left_reals

implicit none ; private

integer, parameter, public :: MAX_PARAM_FILES = 5 !< Maximum number of parameter files.
integer, parameter :: INPUT_STR_LENGTH = 320 !< Maximum line length in parameter file.
integer, parameter :: FILENAME_LENGTH = 200  !< Maximum number of characters in file names.

! The all_PEs_read option should be eliminated with post-riga shared code.
logical :: all_PEs_read = .false. !< If true, all PEs read the input files
                                  !! TODO: Eliminate this parameter

!>@{ Default values for parameters
logical, parameter :: report_unused_default = .true.
logical, parameter :: unused_params_fatal_default = .false.
logical, parameter :: log_to_stdout_default = .false.
logical, parameter :: complete_doc_default = .true.
logical, parameter :: minimal_doc_default = .true.
!>@}

!> The valid lines extracted from an input parameter file without comments
type, private :: file_data_type ; private
  integer :: num_lines = 0 !< The number of lines in this type
  character(len=INPUT_STR_LENGTH), pointer, dimension(:) :: line => NULL() !< The line content
  logical,                         pointer, dimension(:) :: line_used => NULL() !< If true, the line has been read
end type file_data_type

!> A link in the list of variables that have already had override warnings issued
type :: link_parameter ; private
  type(link_parameter), pointer :: next => NULL() !< Facilitates linked list
  character(len=80) :: name                       !< Parameter name
  logical :: hasIssuedOverrideWarning = .false.   !< Has a default value
end type link_parameter

!> Specify the active parameter block
type :: parameter_block ; private
  character(len=240) :: name = ''   !< The active parameter block name
end type parameter_block

!> A structure that can be parsed to read and document run-time parameters.
type, public :: param_file_type ; private
  integer  :: nfiles = 0            !< The number of open files.
  integer  :: iounit(MAX_PARAM_FILES) !< The unit numbers of open files.
  character(len=FILENAME_LENGTH)  :: filename(MAX_PARAM_FILES) !< The names of the open files.
  logical  :: NetCDF_file(MAX_PARAM_FILES) !< If true, the input file is in NetCDF.
                                    ! This is not yet implemented.
  type(file_data_type) :: param_data(MAX_PARAM_FILES) !< Structures that contain
                                    !! the valid data lines from the parameter
                                    !! files, enabling all subsequent reads of
                                    !! parameter data to occur internally.
  logical  :: report_unused = report_unused_default !< If true, report any
                                    !! parameter lines that are not used in the run.
  logical  :: unused_params_fatal = unused_params_fatal_default  !< If true, kill
                                    !! the run if there are any unused parameters.
  logical  :: log_to_stdout = log_to_stdout_default !< If true, all log
                                    !! messages are also sent to stdout.
  logical  :: log_open = .false.    !< True if the log file has been opened.
  integer  :: stdout                !< The unit number from stdout().
  integer  :: stdlog                !< The unit number from stdlog().
  character(len=240) :: doc_file    !< A file where all run-time parameters, their
                                    !! settings and defaults are documented.
  logical  :: complete_doc = complete_doc_default !< If true, document all
                                    !! run-time parameters.
  logical  :: minimal_doc = minimal_doc_default !< If true, document only those
                                    !! run-time parameters that differ from defaults.
  type(doc_type), pointer :: doc => NULL() !< A structure that contains information
                                    !! related to parameter documentation.
  type(link_parameter), pointer :: chain => NULL() !< Facilitates linked list
  type(parameter_block), pointer :: blockName => NULL() !< Name of active parameter block
end type param_file_type

public read_param, open_param_file, close_param_file, log_param, log_version
public doc_param, get_param
public clearParameterBlock, openParameterBlock, closeParameterBlock

!> An overloaded interface to read various types of parameters
interface read_param
  module procedure read_param_int, read_param_real, read_param_logical, &
                   read_param_char, read_param_char_array, read_param_time, &
                   read_param_int_array, read_param_real_array
end interface
!> An overloaded interface to log the values of various types of parameters
interface log_param
  module procedure log_param_int, log_param_real, log_param_logical, &
                   log_param_char, log_param_time, &
                   log_param_int_array, log_param_real_array
end interface
!> An overloaded interface to read and log the values of various types of parameters
interface get_param
  module procedure get_param_int, get_param_real, get_param_logical, &
                   get_param_char, get_param_char_array, get_param_time, &
                   get_param_int_array, get_param_real_array
end interface

!> An overloaded interface to log version information about modules
interface log_version
  module procedure log_version_cs, log_version_plain
end interface

contains

!> Make the contents of a parameter input file availalble in a param_file_type
subroutine open_param_file(filename, CS, checkable, component, doc_file_dir)
  character(len=*),           intent(in) :: filename !< An input file name, optionally with the full path
  type(param_file_type),   intent(inout) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  logical,          optional, intent(in) :: checkable   !< If this is false, it disables checks of this
                                         !! file for unused parameters.  The default is True.
  character(len=*), optional, intent(in) :: component   !< If present, this component name is used
                                         !! to generate parameter documentation file names; the default is"MOM"
  character(len=*), optional, intent(in) :: doc_file_dir !< An optional directory in which to write out
                                         !! the documentation files.  The default is effectively './'.

  ! Local variables
  logical :: file_exists, unit_in_use, Netcdf_file, may_check
  integer :: ios, iounit, strlen, i
  character(len=240) :: doc_path
  type(parameter_block), pointer :: block => NULL()

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

  if (associated(CS%blockName)) deallocate(CS%blockName)
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
  call doc_init(doc_path, CS%doc, minimal=CS%minimal_doc, complete=CS%complete_doc, &
                layout=CS%complete_doc, debugging=CS%complete_doc)

end subroutine open_param_file

!> Close any open input files and deallocate memory associated with this param_file_type.
!! To use this type again, open_param_file would have to be called again.
subroutine close_param_file(CS, quiet_close, component)
  type(param_file_type),   intent(inout) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  logical,          optional, intent(in) :: quiet_close !< if present and true, do not do any
                                         !! logging with this call.
  character(len=*), optional, intent(in) :: component   !< If present, this component name is used
                                         !! to generate parameter documentation file names
  ! Local variables
  logical :: all_default
  character(len=128) :: docfile_default
  character(len=40)  :: mdl   ! This module's name.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  integer :: i, n, num_unused

  if (present(quiet_close)) then ; if (quiet_close) then
    do i = 1, CS%nfiles
      if (all_PEs_read .or. is_root_pe()) close(CS%iounit(i))
      call MOM_mesg("close_param_file: "// trim(CS%filename(i))// &
                    " has been closed successfully.", 5)
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
  docfile_default = "MOM_parameter_doc"
  if (present(component)) docfile_default = trim(component)//"_parameter_doc"

  all_default = (CS%log_to_stdout .eqv. log_to_stdout_default)
  all_default = all_default .and. (trim(CS%doc_file) == trim(docfile_default))
  if (len_trim(CS%doc_file) > 0) then
    all_default = all_default .and. (CS%complete_doc .eqv. complete_doc_default)
    all_default = all_default .and. (CS%minimal_doc .eqv. minimal_doc_default)
  endif

  mdl = "MOM_file_parser"
  call log_version(CS, mdl, version, "", debugging=.true., log_to_all=.true., all_default=all_default)
  call log_param(CS, mdl, "SEND_LOG_TO_STDOUT", CS%log_to_stdout, &
                 "If true, all log messages are also sent to stdout.", &
                 default=log_to_stdout_default)
  call log_param(CS, mdl, "REPORT_UNUSED_PARAMS", CS%report_unused, &
                 "If true, report any parameter lines that are not used "//&
                 "in the run.", default=report_unused_default, &
                 debuggingParam=.true.)
  call log_param(CS, mdl, "FATAL_UNUSED_PARAMS", CS%unused_params_fatal, &
                 "If true, kill the run if there are any unused "//&
                 "parameters.", default=unused_params_fatal_default, &
                 debuggingParam=.true.)
  call log_param(CS, mdl, "DOCUMENT_FILE", CS%doc_file, &
                 "The basename for files where run-time parameters, their "//&
                 "settings, units and defaults are documented. Blank will "//&
                 "disable all parameter documentation.", default=docfile_default)
  if (len_trim(CS%doc_file) > 0) then
    call log_param(CS, mdl, "COMPLETE_DOCUMENTATION",  CS%complete_doc, &
                  "If true, all run-time parameters are "//&
                  "documented in "//trim(CS%doc_file)//&
                  ".all .", default=complete_doc_default)
    call log_param(CS, mdl, "MINIMAL_DOCUMENTATION", CS%minimal_doc, &
                  "If true, non-default run-time parameters are "//&
                  "documented in "//trim(CS%doc_file)//&
                  ".short .", default=minimal_doc_default)
  endif

  num_unused = 0
  do i = 1, CS%nfiles
    if (is_root_pe() .and. (CS%report_unused .or. &
                            CS%unused_params_fatal)) then
      ! Check for unused lines.
      do n=1,CS%param_data(i)%num_lines
        if (.not.CS%param_data(i)%line_used(n)) then
          num_unused = num_unused + 1
          if (CS%report_unused) &
            call MOM_error(WARNING, "Unused line in "//trim(CS%filename(i))// &
                            " : "//trim(CS%param_data(i)%line(n)))
        endif
      enddo
    endif

    if (all_PEs_read .or. is_root_pe()) close(CS%iounit(i))
    call MOM_mesg("close_param_file: "// trim(CS%filename(i))// &
                  " has been closed successfully.", 5)
    CS%iounit(i) = -1
    CS%filename(i) = ''
    CS%NetCDF_file(i) = .false.
    deallocate (CS%param_data(i)%line)
    deallocate (CS%param_data(i)%line_used)
  enddo
  deallocate(CS%blockName)

  if (is_root_pe() .and. (num_unused>0) .and. CS%unused_params_fatal) &
    call MOM_error(FATAL, "Run stopped because of unused parameter lines.")

  CS%log_open = .false.
  call doc_end(CS%doc)

end subroutine close_param_file

!> Read the contents of a parameter input file, and store the contents in a
!! file_data_type after removing comments and simplifying white space
subroutine populate_param_data(iounit, filename, param_data)
  integer,                 intent(in) :: iounit !< The IO unit number that is open for filename
  character(len=*),        intent(in) :: filename !< An input file name, optionally with the full path
  type(file_data_type), intent(inout) :: param_data !< A list of the input lines that set parameters
                                                !! after comments have been stripped out.

  ! Local variables
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


!> Return True if a /* appears on this line without a closing */
function openMultiLineComment(string)
  character(len=*), intent(in) :: string  !< The input string to process
  logical                      :: openMultiLineComment

  ! Local variables
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

!> Return True if a */ appears on this line
function closeMultiLineComment(string)
  character(len=*), intent(in) :: string  !< The input string to process
  logical                      :: closeMultiLineComment
! True if a */ appears on this line
  closeMultiLineComment = .false.
  if (index(string, "*/")>0) closeMultiLineComment=.true.
end function closeMultiLineComment

!> Find position of last character before any comments, As marked by "!", "//", or "/*"
!! following F90, C++, or C syntax
function lastNonCommentIndex(string)
  character(len=*), intent(in) :: string  !< The input string to process
  integer                      :: lastNonCommentIndex

  ! Local variables
  integer :: icom, last

  ! This subroutine is the only place where a comment needs to be defined
  last = len_trim(string)
  icom = index(string(:last), "!") ; if (icom > 0) last = icom-1 ! F90 style
  icom = index(string(:last), "//") ; if (icom > 0) last = icom-1 ! C++ style
  icom = index(string(:last), "/*") ; if (icom > 0) last = icom-1 ! C style
  lastNonCommentIndex = last
end function lastNonCommentIndex

!> Find position of last non-blank character before any comments
function lastNonCommentNonBlank(string)
  character(len=*), intent(in) :: string  !< The input string to process
  integer                      :: lastNonCommentNonBlank

  lastNonCommentNonBlank = len_trim(string(:lastNonCommentIndex(string))) ! Ignore remaining trailing blanks
end function lastNonCommentNonBlank

!> Returns a string with tabs replaced by a blank
function replaceTabs(string)
  character(len=*), intent(in) :: string  !< The input string to process
  character(len=len(string))   :: replaceTabs

  integer :: i

  do i=1, len(string)
    if (string(i:i)==achar(9)) then
      replaceTabs(i:i)=" "
    else
      replaceTabs(i:i)=string(i:i)
    endif
  enddo
end function replaceTabs

!> Trims comments and leading blanks from string
function removeComments(string)
  character(len=*), intent(in) :: string  !< The input string to process
  character(len=len(string))   :: removeComments

  integer :: last

  removeComments=repeat(" ",len(string))
  last = lastNonCommentNonBlank(string)
  removeComments(:last)=adjustl(string(:last)) ! Copy only the non-comment part of string
end function removeComments

!> Constructs a string with all repeated whitespace replaced with single blanks
!! and insert white space where it helps delineate tokens (e.g. around =)
function simplifyWhiteSpace(string)
  character(len=*), intent(in) :: string !< A string to modify to simpify white space
  character(len=len(string)+16)   :: simplifyWhiteSpace

  ! Local variables
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

!> This subroutine reads the value of an integer model parameter from a parameter file.
subroutine read_param_int(CS, varname, value, fail_if_missing)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  integer,             intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,      optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  ! Local variables
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

!> This subroutine reads the values of an array of integer model parameters from a parameter file.
subroutine read_param_int_array(CS, varname, value, fail_if_missing)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  integer, dimension(:),  intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,      optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  ! Local variables
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

!> This subroutine reads the value of a real model parameter from a parameter file.
subroutine read_param_real(CS, varname, value, fail_if_missing, scale)
  type(param_file_type), intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),      intent(in) :: varname !< The case-sensitive name of the parameter to read
  real,               intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,     optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  real,        optional, intent(in) :: scale   !< A scaling factor that the parameter is multiplied
                                         !! by before it is returned.

  ! Local variables
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical            :: found, defined

  call get_variable_line(CS, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read(value_string(1),*,err=1003) value
    if (present(scale)) value = scale*value
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

!> This subroutine reads the values of an array of real model parameters from a parameter file.
subroutine read_param_real_array(CS, varname, value, fail_if_missing, scale)
  type(param_file_type), intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),      intent(in) :: varname !< The case-sensitive name of the parameter to read
  real, dimension(:), intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,     optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  real,        optional, intent(in) :: scale   !< A scaling factor that the parameter is multiplied
                                         !! by before it is returned.

  ! Local variables
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  logical                         :: found, defined

  call get_variable_line(CS, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    read(value_string(1),*,end=991,err=1004) value
991 continue
    if (present(scale)) value(:) = scale*value(:)
    return
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

!> This subroutine reads the value of a character string model parameter from a parameter file.
subroutine read_param_char(CS, varname, value, fail_if_missing)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  character(len=*),    intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,      optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  ! Local variables
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

!> This subroutine reads the values of an array of character string model parameters from a parameter file.
subroutine read_param_char_array(CS, varname, value, fail_if_missing)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  character(len=*), dimension(:), intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,      optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file

  ! Local variables
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

!> This subroutine reads the value of a logical model parameter from a parameter file.
subroutine read_param_logical(CS, varname, value, fail_if_missing)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  logical,             intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  logical,      optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file

  ! Local variables
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

!> This subroutine reads the value of a time_type model parameter from a parameter file.
subroutine read_param_time(CS, varname, value, timeunit, fail_if_missing, date_format)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  type(time_type),     intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file
  real,         optional, intent(in) :: timeunit !< The number of seconds in a time unit for real-number input.
  logical,      optional, intent(in) :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,     optional, intent(out) :: date_format !< If present, this indicates whether this
                                         !! parameter was read in a date format, so that it can
                                         !! later be logged in the same format.

  ! Local variables
  character(len=INPUT_STR_LENGTH) :: value_string(1)
  character(len=240) :: err_msg
  logical            :: found, defined
  real               :: real_time, time_unit
  integer            :: vals(7)

  if (present(date_format)) date_format = .false.

  call get_variable_line(CS, varname, found, defined, value_string)
  if (found .and. defined .and. (LEN_TRIM(value_string(1)) > 0)) then
    ! Determine whether value string should be parsed for a real number
    ! or a date, in either a string format or a comma-delimited list of values.
    if ((INDEX(value_string(1),'-') > 0) .and. &
        (INDEX(value_string(1),'-',back=.true.) > INDEX(value_string(1),'-'))) then
      ! There are two dashes, so this must be a date format.
      value = set_date(value_string(1), err_msg=err_msg)
      if (LEN_TRIM(err_msg) > 0) call MOM_error(FATAL,'read_param_time: '//&
          trim(err_msg)//' in integer list read error for time-type variable '//&
          trim(varname)// ' parsing "'//trim(value_string(1))//'"')
      if (present(date_format)) date_format = .true.
    elseif (INDEX(value_string(1),',') > 0) then
      ! Initialize vals with an invalid date.
      vals(:) = (/ -999, -999, -999, 0, 0, 0, 0 /)
      read(value_string(1), *, end=995, err=1005) vals
      995 continue
      if ((vals(1) < 0) .or. (vals(2) < 0) .or. (vals(3) < 0)) &
        call MOM_error(FATAL,'read_param_time: integer list read error for time-type variable '//&
                       trim(varname)// ' parsing "'//trim(value_string(1))//'"')
      value = set_date(vals(1), vals(2), vals(3), vals(4), vals(5), vals(6), &
                       vals(7), err_msg=err_msg)
      if (LEN_TRIM(err_msg) > 0) call MOM_error(FATAL,'read_param_time: '//&
          trim(err_msg)//' in integer list read error for time-type variable '//&
          trim(varname)// ' parsing "'//trim(value_string(1))//'"')
      if (present(date_format)) date_format = .true.
    else
      time_unit = 1.0 ; if (present(timeunit)) time_unit = timeunit
      read( value_string(1), *) real_time
      value = real_to_time(real_time*time_unit)
    endif
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
  return

  1005 call MOM_error(FATAL, 'read_param_time: read error for time-type variable '//&
                             trim(varname)// ' parsing "'//trim(value_string(1))//'"')
end subroutine read_param_time

!> This function removes single and double quotes from a character string
function strip_quotes(val_str)
  character(len=*) :: val_str !< The character string to work on
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

!> This subtoutine extracts the contents of lines in the param_file_type that refer to
!! a named parameter.  The value_string that is returned must be interepreted in a way
!! that depends on the type of this variable.
subroutine get_variable_line(CS, varname, found, defined, value_string, paramIsLogical)
  type(param_file_type),  intent(in) :: CS      !< The control structure for the file_parser module,
                                                !! it is also a structure to parse for run-time parameters
  character(len=*),       intent(in) :: varname !< The case-sensitive name of the parameter to read
  logical,               intent(out) :: found   !< If true, this parameter has been found in CS
  logical,               intent(out) :: defined !< If true, this parameter is set (or true) in the CS
  character(len=*),      intent(out) :: value_string(:) !< A string that encodes the new value
  logical, optional,      intent(in) :: paramIsLogical  !< If true, this is a logical parameter
                                                !! that can be simply defined without parsing a value_string.

  ! Local variables
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
  oval = 0; ival = 0
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

!> Record that a line has been used to set a parameter
subroutine flag_line_as_read(line_used, count)
  logical, dimension(:), pointer :: line_used !< A structure indicating which lines have been read
  integer,            intent(in) :: count !< The parameter on this line number has been read
  line_used(count) = .true.
end subroutine flag_line_as_read

!> Returns true if an override warning has been issued for the variable varName
function overrideWarningHasBeenIssued(chain, varName)
  type(link_parameter), pointer    :: chain   !< The linked list of variables that have already had
                                              !! override warnings issued
  character(len=*),     intent(in) :: varName !< The name of the variable being queried for warnings
  logical                          :: overrideWarningHasBeenIssued
  ! Local variables
  type(link_parameter), pointer :: newLink => NULL(), this => NULL()

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
subroutine log_version_cs(CS, modulename, version, desc, log_to_all, all_default, layout, debugging)
  type(param_file_type),      intent(in) :: CS         !< File parser type
  character(len=*),           intent(in) :: modulename !< Name of calling module
  character(len=*),           intent(in) :: version    !< Version string of module
  character(len=*), optional, intent(in) :: desc       !< Module description
  logical,          optional, intent(in) :: log_to_all !< If present and true, log this parameter to the
                                                       !! ..._doc.all files, even if this module also has layout
                                                       !! or debugging parameters.
  logical,          optional, intent(in) :: all_default !< If true, all parameters take their default values.
  logical,          optional, intent(in) :: layout     !< If present and true, this module has layout parameters.
  logical,          optional, intent(in) :: debugging  !< If present and true, this module has debugging parameters.
  ! Local variables
  character(len=240) :: mesg

  mesg = trim(modulename)//": "//trim(version)
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  if (present(desc)) call doc_module(CS%doc, modulename, desc, log_to_all, all_default, layout, debugging)

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

!> Log the name and value of an integer model parameter in documentation files.
subroutine log_param_int(CS, modulename, varname, value, desc, units, &
                         default, layoutParam, debuggingParam, like_default)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the module using this parameter
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  integer,                    intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in) :: units   !< The units of this parameter
  integer,          optional, intent(in) :: default !< The default value of the parameter
  logical,          optional, intent(in) :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.

  character(len=240) :: mesg, myunits

  write(mesg, '("  ",a," ",a,": ",a)') trim(modulename), trim(varname), trim(left_int(value))
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits=" "; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default, &
                   layoutParam=layoutParam, debuggingParam=debuggingParam, like_default=like_default)

end subroutine log_param_int

!> Log the name and values of an array of integer model parameter in documentation files.
subroutine log_param_int_array(CS, modulename, varname, value, desc, &
                               units, default, layoutParam, debuggingParam, like_default)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the module using this parameter
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  integer, dimension(:),      intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in) :: units   !< The units of this parameter
  integer,          optional, intent(in) :: default !< The default value of the parameter
  logical,          optional, intent(in) :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.

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
                   layoutParam=layoutParam, debuggingParam=debuggingParam, like_default=like_default)

end subroutine log_param_int_array

!> Log the name and value of a real model parameter in documentation files.
subroutine log_param_real(CS, modulename, varname, value, desc, units, &
                          default, debuggingParam, like_default)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the calling module
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  real,                       intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in) :: units   !< The units of this parameter
  real,             optional, intent(in) :: default !< The default value of the parameter
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.

  character(len=240) :: mesg, myunits

  write(mesg, '("  ",a," ",a,": ",a)') &
    trim(modulename), trim(varname), trim(left_real(value))
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits="not defined"; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default, &
                   debuggingParam=debuggingParam, like_default=like_default)

end subroutine log_param_real

!> Log the name and values of an array of real model parameter in documentation files.
subroutine log_param_real_array(CS, modulename, varname, value, desc, &
                                units, default, debuggingParam, like_default)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the calling module
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  real, dimension(:),         intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                             !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in) :: units   !< The units of this parameter
  real,             optional, intent(in) :: default !< The default value of the parameter
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.

  character(len=:), allocatable :: mesg
  character(len=240) :: myunits

 !write(mesg, '("  ",a," ",a,": ",ES19.12,99(",",ES19.12))') &
 !write(mesg, '("  ",a," ",a,": ",G,99(",",G))') &
 !  trim(modulename), trim(varname), value
  mesg = "  " // trim(modulename) // " " // trim(varname) // ": " // trim(left_reals(value))
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits="not defined"; if (present(units)) write(myunits(1:240),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default, &
                   debuggingParam=debuggingParam, like_default=like_default)

end subroutine log_param_real_array

!> Log the name and value of a logical model parameter in documentation files.
subroutine log_param_logical(CS, modulename, varname, value, desc, &
                             units, default, layoutParam, debuggingParam, like_default)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the calling module
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  logical,                    intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in) :: units   !< The units of this parameter
  logical,          optional, intent(in) :: default !< The default value of the parameter
  logical,          optional, intent(in) :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.

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
                   layoutParam=layoutParam, debuggingParam=debuggingParam, like_default=like_default)

end subroutine log_param_logical

!> Log the name and value of a character string model parameter in documentation files.
subroutine log_param_char(CS, modulename, varname, value, desc, units, &
                          default, layoutParam, debuggingParam, like_default)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the calling module
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  character(len=*),           intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in) :: units   !< The units of this parameter
  character(len=*), optional, intent(in) :: default !< The default value of the parameter
  logical,          optional, intent(in) :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.

  character(len=1024) :: mesg, myunits

  write(mesg, '("  ",a," ",a,": ",a)') &
    trim(modulename), trim(varname), trim(value)
  if (is_root_pe()) then
    if (CS%log_open) write(CS%stdlog,'(a)') trim(mesg)
    if (CS%log_to_stdout) write(CS%stdout,'(a)') trim(mesg)
  endif

  myunits=" "; if (present(units)) write(myunits(1:1024),'(A)') trim(units)
  if (present(desc)) &
    call doc_param(CS%doc, varname, desc, myunits, value, default, &
                   layoutParam=layoutParam, debuggingParam=debuggingParam, like_default=like_default)

end subroutine log_param_char

!> This subroutine writes the value of a time-type parameter to a log file,
!! along with its name and the module it came from.
subroutine log_param_time(CS, modulename, varname, value, desc, units, &
                          default, timeunit, layoutParam, debuggingParam, log_date, like_default)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: modulename !< The name of the calling module
  character(len=*),           intent(in) :: varname !< The name of the parameter to log
  type(time_type),            intent(in) :: value   !< The value of the parameter to log
  character(len=*), optional, intent(in) :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in) :: units   !< The units of this parameter
  type(time_type),  optional, intent(in) :: default !< The default value of the parameter
  real,             optional, intent(in) :: timeunit !< The number of seconds in a time unit for
                                         !! real-number output.
  logical,          optional, intent(in) :: log_date   !< If true, log the time_type in date format.
                                         !! If missing the default is false.
  logical,          optional, intent(in) :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in) :: like_default !< If present and true, log this parameter as
                                         !! though it has the default value, even if there is no default.

  ! Local variables
  real :: real_time, real_default
  logical :: use_timeunit, date_format
  character(len=240) :: mesg, myunits
  character(len=80) :: date_string, default_string
  integer :: days, secs, ticks, ticks_per_sec

  use_timeunit = .false.
  date_format = .false. ; if (present(log_date)) date_format = log_date

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
    if (date_format) then
      myunits='[date]'

      date_string = convert_date_to_string(value)
      if (present(default)) then
        default_string = convert_date_to_string(default)
        call doc_param(CS%doc, varname, desc, myunits, date_string, &
                       default=default_string, layoutParam=layoutParam, &
                       debuggingParam=debuggingParam, like_default=like_default)
      else
        call doc_param(CS%doc, varname, desc, myunits, date_string, &
                       layoutParam=layoutParam, debuggingParam=debuggingParam, like_default=like_default)
      endif
    elseif (use_timeunit) then
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
        call doc_param(CS%doc, varname, desc, myunits, real_time, real_default, like_default=like_default)
      else
        call doc_param(CS%doc, varname, desc, myunits, real_time, like_default=like_default)
      endif
    else
      call doc_param(CS%doc, varname, desc, value, default, units=units, like_default=like_default)
    endif
  endif

end subroutine log_param_time

!> This function converts a date into a string, valid with ticks and for dates up to year 99,999,999
function convert_date_to_string(date) result(date_string)
  type(time_type), intent(in) :: date !< The date to be translated into a string.
  character(len=40) :: date_string    !< A date string in a format like YYYY-MM-DD HH:MM:SS.sss

  ! Local variables
  character(len=40) :: sub_string
  real    :: real_secs
  integer :: yrs, mons, days, hours, mins, secs, ticks, ticks_per_sec

  call get_date(date, yrs, mons, days, hours, mins, secs, ticks)
  write (date_string, '(i8.4)') yrs
  write (sub_string, '("-", i2.2, "-", I2.2, " ", i2.2, ":", i2.2, ":")') &
         mons, days, hours, mins
  date_string = trim(adjustl(date_string)) // trim(sub_string)
  if (ticks > 0) then
    ticks_per_sec = get_ticks_per_second()
    real_secs = secs + ticks/ticks_per_sec
    if (ticks_per_sec <= 100) then
      write (sub_string, '(F7.3)') real_secs
    else
      write (sub_string, '(F10.6)') real_secs
    endif
  else
    write (sub_string, '(i2.2)') secs
  endif
  date_string = trim(date_string) // trim(adjustl(sub_string))

end function convert_date_to_string

!> This subroutine reads the value of an integer model parameter from a parameter file
!! and logs it in documentation files.
subroutine get_param_int(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, &
               static_value, layoutParam, debuggingParam)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  integer,                    intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  integer,          optional, intent(in)    :: default !< The default value of the parameter
  integer,          optional, intent(in)    :: static_value !< If this parameter is static, it takes
                                         !! this value, which can be compared for consistency with
                                         !! what is in the parameter file.
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  logical,          optional, intent(in)    :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file

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
                       default, layoutParam, debuggingParam)
  endif

end subroutine get_param_int

!> This subroutine reads the values of an array of integer model parameters from a parameter file
!! and logs them in documentation files.
subroutine get_param_int_array(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, &
               static_value, layoutParam, debuggingParam)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  integer, dimension(:),      intent(inout) :: value   !< The value of the parameter that may be reset
                                         !! from the parameter file
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  integer,          optional, intent(in)    :: default !< The default value of the parameter
  integer,          optional, intent(in)    :: static_value !< If this parameter is static, it takes
                                         !! this value, which can be compared for consistency with
                                         !! what is in the parameter file.
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  logical,          optional, intent(in)    :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file

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
                             units, default, layoutParam, debuggingParam)
  endif

end subroutine get_param_int_array

!> This subroutine reads the value of a real model parameter from a parameter file
!! and logs it in documentation files.
subroutine get_param_real(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, &
               static_value, debuggingParam, scale, unscaled)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  real,                       intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  real,             optional, intent(in)    :: default !< The default value of the parameter
  real,             optional, intent(in)    :: static_value !< If this parameter is static, it takes
                                         !! this value, which can be compared for consistency with
                                         !! what is in the parameter file.
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  real,             optional, intent(in)    :: scale   !< A scaling factor that the parameter is
                                         !! multiplied by before it is returned.
  real,             optional, intent(out)   :: unscaled !< The value of the parameter that would be
                                         !! returned without any multiplication by a scaling factor.

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
                        default, debuggingParam)
  endif

  if (present(unscaled)) unscaled = value
  if (present(scale)) value = scale*value

end subroutine get_param_real

!> This subroutine reads the values of an array of real model parameters from a parameter file
!! and logs them in documentation files.
subroutine get_param_real_array(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, debuggingParam, &
               static_value, scale, unscaled)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  real, dimension(:),         intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  real,             optional, intent(in)    :: default !< The default value of the parameter
  real,             optional, intent(in)    :: static_value !< If this parameter is static, it takes
                                         !! this value, which can be compared for consistency with
                                         !! what is in the parameter file.
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  real,             optional, intent(in)    :: scale   !< A scaling factor that the parameter is
                                         !! multiplied by before it is returned.
  real, dimension(:), optional, intent(out) :: unscaled !< The value of the parameter that would be
                                         !! returned without any multiplication by a scaling factor.

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
                              units, default, debuggingParam)
  endif

  if (present(unscaled)) unscaled(:) = value(:)
  if (present(scale)) value(:) = scale*value(:)

end subroutine get_param_real_array

!> This subroutine reads the value of a character string model parameter from a parameter file
!! and logs it in documentation files.
subroutine get_param_char(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, &
               static_value, layoutParam, debuggingParam)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  character(len=*),           intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  character(len=*), optional, intent(in)    :: default !< The default value of the parameter
  character(len=*), optional, intent(in)    :: static_value !< If this parameter is static, it takes
                                         !! this value, which can be compared for consistency with
                                         !! what is in the parameter file.
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  logical,          optional, intent(in)    :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file

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
                        default, layoutParam, debuggingParam)
  endif

end subroutine get_param_char

!> This subroutine reads the values of an array of character string model parameters
!! from a parameter file and logs them in documentation files.
subroutine get_param_char_array(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, static_value)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  character(len=*), dimension(:), intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  character(len=*), optional, intent(in)    :: default !< The default value of the parameter
  character(len=*), optional, intent(in)    :: static_value !< If this parameter is static, it takes
                                         !! this value, which can be compared for consistency with
                                         !! what is in the parameter file.
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files

  ! Local variables
  logical :: do_read, do_log
  integer :: i, len_tot, len_val
  character(len=1024) :: cat_val

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

!> This subroutine reads the value of a logical model parameter from a parameter file
!! and logs it in documentation files.
subroutine get_param_logical(CS, modulename, varname, value, desc, units, &
               default, fail_if_missing, do_not_read, do_not_log, &
               static_value, layoutParam, debuggingParam)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  logical,                    intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  logical,          optional, intent(in)    :: default !< The default value of the parameter
  logical,          optional, intent(in)    :: static_value !< If this parameter is static, it takes
                                         !! this value, which can be compared for consistency with
                                         !! what is in the parameter file.
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  logical,          optional, intent(in)    :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file

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
                           units, default, layoutParam, debuggingParam)
  endif

end subroutine get_param_logical

!> This subroutine reads the value of a time-type model parameter from a parameter file
!! and logs it in documentation files.
subroutine get_param_time(CS, modulename, varname, value, desc, units, &
                          default, fail_if_missing, do_not_read, do_not_log, &
                          timeunit, static_value, layoutParam, debuggingParam, &
                          log_as_date)
  type(param_file_type),      intent(in)    :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in)    :: modulename !< The name of the calling module
  character(len=*),           intent(in)    :: varname !< The case-sensitive name of the parameter to read
  type(time_type),            intent(inout) :: value   !< The value of the parameter that may be
                                         !! read from the parameter file and logged
  character(len=*), optional, intent(in)    :: desc    !< A description of this variable; if not
                                         !! present, this parameter is not written to a doc file
  character(len=*), optional, intent(in)    :: units   !< The units of this parameter
  type(time_type),  optional, intent(in)    :: default !< The default value of the parameter
  type(time_type),  optional, intent(in)    :: static_value !< If this parameter is static, it takes
                                         !! this value, which can be compared for consistency with
                                         !! what is in the parameter file.
  logical,          optional, intent(in)    :: fail_if_missing !< If present and true, a fatal error occurs
                                         !! if this variable is not found in the parameter file
  logical,          optional, intent(in)    :: do_not_read  !< If present and true, do not read a
                                         !! value for this parameter, although it might be logged.
  logical,          optional, intent(in)    :: do_not_log !< If present and true, do not log this
                                         !! parameter to the documentation files
  real,             optional, intent(in)    :: timeunit !< The number of seconds in a time unit for
                                         !! real-number input to be translated to a time.
  logical,          optional, intent(in)    :: layoutParam !< If present and true, this parameter is
                                         !! logged in the layout parameter file
  logical,          optional, intent(in)    :: debuggingParam !< If present and true, this parameter is
                                         !! logged in the debugging parameter file
  logical,          optional, intent(in)    :: log_as_date  !< If true, log the time_type in date
                                         !! format. The default is false.

  logical :: do_read, do_log, date_format, log_date

  do_read = .true. ; if (present(do_not_read)) do_read = .not.do_not_read
  do_log  = .true. ; if (present(do_not_log))  do_log  = .not.do_not_log
  log_date = .false.

  if (do_read) then
    if (present(default)) value = default
    if (present(static_value)) value = static_value
    call read_param_time(CS, varname, value, timeunit, fail_if_missing, date_format=log_date)
  endif

  if (do_log) then
    if (present(log_as_date)) log_date = log_as_date
    call log_param_time(CS, modulename, varname, value, desc, units, default, &
                        timeunit, layoutParam=layoutParam, &
                        debuggingParam=debuggingParam, log_date=log_date)
  endif

end subroutine get_param_time

! -----------------------------------------------------------------------------

!> Resets the parameter block name to blank
subroutine clearParameterBlock(CS)
  type(param_file_type), intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters

  type(parameter_block), pointer :: block => NULL()
  if (associated(CS%blockName)) then
    block => CS%blockName
    block%name = ''
  else
    if (is_root_pe()) call MOM_error(FATAL, &
      'clearParameterBlock: A clear was attempted before allocation.')
  endif
end subroutine clearParameterBlock

!> Tags blockName onto the end of the active parameter block name
subroutine openParameterBlock(CS,blockName,desc)
  type(param_file_type),      intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters
  character(len=*),           intent(in) :: blockName !< The name of a parameter block being added
  character(len=*), optional, intent(in) :: desc    !< A description of the parameter block being added

  type(parameter_block), pointer :: block => NULL()
  if (associated(CS%blockName)) then
    block => CS%blockName
    block%name = pushBlockLevel(block%name,blockName)
    call doc_openBlock(CS%doc,block%name,desc)
  else
    if (is_root_pe()) call MOM_error(FATAL, &
      'openParameterBlock: A push was attempted before allocation.')
  endif
end subroutine openParameterBlock

!> Remove the lowest level of recursion from the active block name
subroutine closeParameterBlock(CS)
  type(param_file_type), intent(in) :: CS      !< The control structure for the file_parser module,
                                         !! it is also a structure to parse for run-time parameters

  type(parameter_block), pointer :: block => NULL()

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

!> Extends block name (deeper level of parameter block)
function pushBlockLevel(oldblockName,newBlockName)
  character(len=*),        intent(in) :: oldBlockName  !< A sequence of hierarchical parameter block names
  character(len=*),        intent(in) :: newBlockName  !< A new block name to add to the end of the sequence
  character(len=len(oldBlockName)+40) :: pushBlockLevel

  if (len_trim(oldBlockName)>0) then
    pushBlockLevel=trim(oldBlockName)//'%'//trim(newBlockName)
  else
    pushBlockLevel=trim(newBlockName)
  endif
end function pushBlockLevel

!> Truncates block name (shallower level of parameter block)
function popBlockLevel(oldblockName)
  character(len=*),        intent(in) :: oldBlockName !< A sequence of hierarchical parameter block names
  character(len=len(oldBlockName)+40) :: popBlockLevel

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

!> \namespace mom_file_parser
!!
!!  By Robert Hallberg and Alistair Adcroft, updated 9/2013.
!!
!!    The subroutines here parse a set of input files for the value
!!  a named parameter and sets that parameter at run time.  Currently
!!  these files use use one of several formats:
!!    \#define VAR      ! To set the logical VAR to true.
!!    VAR = True        ! To set the logical VAR to true.
!!    \#undef VAR       ! To set the logical VAR to false.
!!    VAR = False       ! To set the logical VAR to false.
!!    \#define VAR 999  ! To set the real or integer VAR to 999.
!!    VAR = 999         ! To set the real or integer VAR to 999.
!!    \#override VAR = 888 ! To override a previously set value.
!!    VAR = 1.1, 2.2, 3.3  ! To set an array of real values.
  ! Note that in the comments above, dOxygen translates \# to # .
!!
!!  In addition, when set by the get_param interface, the values of
!!  parameters are automatically logged, along with defaults, units,
!!  and a description.  It is an error for a variable to be overridden
!!  more than once, and MOM6 has a facility to check for unused lines
!!  to set variables, which may indicate miss-spelled or archaic
!!  parameters.  Parameter names are case-specific, and lines may use
!!  a F90 or C++ style comment, starting with ! or //.

end module MOM_file_parser
