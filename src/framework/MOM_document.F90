!> The subroutines here provide hooks for document generation functions at
!! various levels of granularity.
module MOM_document

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_time_manager,  only : time_type, operator(==), get_time, get_ticks_per_second
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe

implicit none ; private

public doc_param, doc_subroutine, doc_function, doc_module, doc_init, doc_end
public doc_openBlock, doc_closeBlock

!> Document parameter values
interface doc_param
  module procedure doc_param_none, &
                   doc_param_logical, doc_param_logical_array, &
                   doc_param_int,     doc_param_int_array, &
                   doc_param_real,    doc_param_real_array, &
                   doc_param_char, &
                   doc_param_time
end interface

integer, parameter :: mLen = 1240 !< Length of interface/message strings

!> A structure that controls where the documentation occurs, its veborsity and formatting.
type, public :: doc_type ; private
  integer :: unitAll = -1           !< The open unit number for docFileBase + .all.
  integer :: unitShort = -1         !< The open unit number for docFileBase + .short.
  integer :: unitLayout = -1        !< The open unit number for docFileBase + .layout.
  integer :: unitDebugging  = -1    !< The open unit number for docFileBase + .debugging.
  logical :: filesAreOpen = .false. !< True if any files were successfully opened.
  character(len=mLen) :: docFileBase = '' !< The basename of the files where run-time
                                    !! parameters, settings and defaults are documented.
  logical :: complete = .true.      !< If true, document all parameters.
  logical :: minimal = .true.       !< If true, document non-default parameters.
  logical :: layout = .true.        !< If true, document layout parameters.
  logical :: debugging = .true.     !< If true, document debugging parameters.
  logical :: defineSyntax = .false. !< If true, use '\#def' syntax instead of a=b syntax
  logical :: warnOnConflicts = .false. !< Cause a WARNING error if defaults differ.
  integer :: commentColumn = 32     !< Number of spaces before the comment marker.
  integer :: max_line_len = 112     !< The maximum length of message lines.
  type(link_msg), pointer :: chain_msg => NULL() !< Database of messages
  character(len=240) :: blockPrefix = '' !< The full name of the current block.
end type doc_type

!> A linked list of the parameter documentation messages that have been issued so far.
type :: link_msg ; private
  type(link_msg), pointer :: next => NULL()  !< Facilitates linked list
  character(len=80) :: name                  !< Parameter name
  character(len=620) :: msg                  !< Parameter value and default
end type link_msg

character(len=4), parameter :: STRING_TRUE  = 'True'  !< A string for true logicals
character(len=5), parameter :: STRING_FALSE = 'False' !< A string for false logicals

contains

! ----------------------------------------------------------------------

!> This subroutine handles parameter documentation with no value.
subroutine doc_param_none(doc, varname, desc, units)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: varname !< The name of the parameter being documented
  character(len=*), intent(in) :: desc    !< A description of the parameter being documented
  character(len=*), intent(in) :: units   !< The units of the parameter being documented
! This subroutine handles parameter documentation with no value.
  integer :: numspc
  character(len=mLen) :: mesg

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    numspc = max(1,doc%commentColumn-8-len_trim(varname))
    mesg = "#define "//trim(varname)//repeat(" ",numspc)//"!"
    if (len_trim(units) > 0) mesg = trim(mesg)//"   ["//trim(units)//"]"

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc)
  endif
end subroutine doc_param_none

!> This subroutine handles parameter documentation for logicals.
subroutine doc_param_logical(doc, varname, desc, units, val, default, &
                             layoutParam, debuggingParam)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  logical,           intent(in) :: val     !< The value of this parameter
  logical, optional, intent(in) :: default !< The default value of this parameter
  logical, optional, intent(in) :: layoutParam !< If present and true, this is a layout parameter.
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
! This subroutine handles parameter documentation for logicals.
  character(len=mLen) :: mesg
  logical :: equalsDefault

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    if (val) then
      mesg = define_string(doc, varname, STRING_TRUE, units)
    else
      mesg = undef_string(doc, varname, units)
    endif

    equalsDefault = .false.
    if (present(default)) then
      if (val .eqv. default) equalsDefault = .true.
      if (default) then
        mesg = trim(mesg)//" default = "//STRING_TRUE
      else
        mesg = trim(mesg)//" default = "//STRING_FALSE
      endif
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, &
                             layoutParam=layoutParam, debuggingParam=debuggingParam)
  endif
end subroutine doc_param_logical

!> This subroutine handles parameter documentation for arrays of logicals.
subroutine doc_param_logical_array(doc, varname, desc, units, vals, default, &
                                   layoutParam, debuggingParam)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  logical,           intent(in) :: vals(:) !< The array of values to record
  logical, optional, intent(in) :: default !< The default value of this parameter
  logical, optional, intent(in) :: layoutParam !< If present and true, this is a layout parameter.
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
! This subroutine handles parameter documentation for arrays of logicals.
  integer :: i
  character(len=mLen) :: mesg
  character(len=mLen) :: valstring
  logical :: equalsDefault

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    if (vals(1)) then ; valstring = STRING_TRUE ; else ; valstring = STRING_FALSE ; endif
    do i=2,min(size(vals),128)
      if (vals(i)) then
        valstring = trim(valstring)//", "//STRING_TRUE
      else
        valstring = trim(valstring)//", "//STRING_FALSE
      endif
    enddo

    mesg = define_string(doc, varname, valstring, units)

  equalsDefault = .false.
    if (present(default)) then
      equalsDefault = .true.
      do i=1,size(vals) ; if (vals(i) .neqv. default) equalsDefault = .false. ; enddo
      if (default) then
        mesg = trim(mesg)//" default = "//STRING_TRUE
      else
        mesg = trim(mesg)//" default = "//STRING_FALSE
      endif
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, &
                             layoutParam=layoutParam, debuggingParam=debuggingParam)
  endif
end subroutine doc_param_logical_array

!> This subroutine handles parameter documentation for integers.
subroutine doc_param_int(doc, varname, desc, units, val, default, &
                         layoutParam, debuggingParam)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  integer,           intent(in) :: val     !< The value of this parameter
  integer, optional, intent(in) :: default !< The default value of this parameter
  logical, optional, intent(in) :: layoutParam !< If present and true, this is a layout parameter.
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
! This subroutine handles parameter documentation for integers.
  character(len=mLen) :: mesg
  character(len=doc%commentColumn)  :: valstring
  logical :: equalsDefault

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    valstring = int_string(val)
    mesg = define_string(doc, varname, valstring, units)

    equalsDefault = .false.
    if (present(default)) then
      if (val == default) equalsDefault = .true.
      mesg = trim(mesg)//" default = "//(trim(int_string(default)))
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, &
                             layoutParam=layoutParam, debuggingParam=debuggingParam)
  endif
end subroutine doc_param_int

!> This subroutine handles parameter documentation for arrays of integers.
subroutine doc_param_int_array(doc, varname, desc, units, vals, default, &
                               layoutParam, debuggingParam)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  integer,           intent(in) :: vals(:) !< The array of values to record
  integer, optional, intent(in) :: default !< The default value of this parameter
  logical, optional, intent(in) :: layoutParam !< If present and true, this is a layout parameter.
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
! This subroutine handles parameter documentation for arrays of integers.
  integer :: i
  character(len=mLen) :: mesg
  character(len=mLen)  :: valstring
  logical :: equalsDefault

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    valstring = int_string(vals(1))
    do i=2,min(size(vals),128)
      valstring = trim(valstring)//", "//trim(int_string(vals(i)))
    enddo

    mesg = define_string(doc, varname, valstring, units)

    equalsDefault = .false.
    if (present(default)) then
      equalsDefault = .true.
      do i=1,size(vals) ; if (vals(i) /= default) equalsDefault = .false. ; enddo
      mesg = trim(mesg)//" default = "//(trim(int_string(default)))
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, &
                             layoutParam=layoutParam, debuggingParam=debuggingParam)
  endif

end subroutine doc_param_int_array

!> This subroutine handles parameter documentation for reals.
subroutine doc_param_real(doc, varname, desc, units, val, default, debuggingParam)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  real,              intent(in) :: val     !< The value of this parameter
  real,    optional, intent(in) :: default !< The default value of this parameter
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
! This subroutine handles parameter documentation for reals.
  character(len=mLen) :: mesg
  character(len=doc%commentColumn)  :: valstring
  logical :: equalsDefault

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    valstring = real_string(val)
    mesg = define_string(doc, varname, valstring, units)

    equalsDefault = .false.
    if (present(default)) then
      if (val == default) equalsDefault = .true.
      mesg = trim(mesg)//" default = "//trim(real_string(default))
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, debuggingParam=debuggingParam)
  endif
end subroutine doc_param_real

!> This subroutine handles parameter documentation for arrays of reals.
subroutine doc_param_real_array(doc, varname, desc, units, vals, default, debuggingParam)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  real,              intent(in) :: vals(:) !< The array of values to record
  real,    optional, intent(in) :: default !< The default value of this parameter
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
! This subroutine handles parameter documentation for arrays of reals.
  integer :: i
  character(len=mLen) :: mesg
  character(len=mLen) :: valstring
  logical :: equalsDefault

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    valstring = trim(real_array_string(vals(:)))

    mesg = define_string(doc, varname, valstring, units)

    equalsDefault = .false.
    if (present(default)) then
      equalsDefault = .true.
      do i=1,size(vals) ; if (vals(i) /= default) equalsDefault = .false. ; enddo
      mesg = trim(mesg)//" default = "//trim(real_string(default))
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, debuggingParam=debuggingParam)
  endif

end subroutine doc_param_real_array

!> This subroutine handles parameter documentation for character strings.
subroutine doc_param_char(doc, varname, desc, units, val, default, &
                          layoutParam, debuggingParam)
  type(doc_type),    pointer    :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: varname !< The name of the parameter being documented
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  character(len=*),  intent(in) :: units   !< The units of the parameter being documented
  character(len=*),  intent(in) :: val     !< The value of the parameter
  character(len=*), &
           optional, intent(in) :: default !< The default value of this parameter
  logical, optional, intent(in) :: layoutParam !< If present and true, this is a layout parameter.
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.
! This subroutine handles parameter documentation for character strings.
  character(len=mLen) :: mesg
  logical :: equalsDefault

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    mesg = define_string(doc, varname, '"'//trim(val)//'"', units)

    equalsDefault = .false.
    if (present(default)) then
      if (trim(val) == trim(default)) equalsDefault = .true.
      mesg = trim(mesg)//' default = "'//trim(adjustl(default))//'"'
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, &
                             layoutParam=layoutParam, debuggingParam=debuggingParam)
  endif

end subroutine doc_param_char

!> This subroutine handles documentation for opening a parameter block.
subroutine doc_openBlock(doc, blockName, desc)
  type(doc_type),   pointer    :: doc       !< A pointer to a structure that controls where the
                                            !! documentation occurs and its formatting
  character(len=*), intent(in) :: blockName !< The name of the parameter block being opened
  character(len=*), optional, intent(in) :: desc !< A description of the parameter block being opened
! This subroutine handles documentation for opening a parameter block.
  character(len=mLen) :: mesg
  character(len=doc%commentColumn) :: valstring

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    mesg = trim(blockName)//'%'

    if (present(desc)) then
      call writeMessageAndDesc(doc, mesg, desc)
    else
      call writeMessageAndDesc(doc, mesg, '')
    endif
  endif
  doc%blockPrefix = trim(doc%blockPrefix)//trim(blockName)//'%'
end subroutine doc_openBlock

!> This subroutine handles documentation for closing a parameter block.
subroutine doc_closeBlock(doc, blockName)
  type(doc_type),   pointer    :: doc !< A pointer to a structure that controls where the
                                      !! documentation occurs and its formatting
  character(len=*), intent(in) :: blockName !< The name of the parameter block being closed
! This subroutine handles documentation for closing a parameter block.
  character(len=mLen) :: mesg
  character(len=doc%commentColumn) :: valstring
  integer :: i

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    mesg = '%'//trim(blockName)

    call writeMessageAndDesc(doc, mesg, '')
  endif
  i = index(trim(doc%blockPrefix), trim(blockName)//'%', .true.)
  if (i>1) then
    doc%blockPrefix = trim(doc%blockPrefix(1:i-1))
  else
    doc%blockPrefix = ''
  endif
end subroutine doc_closeBlock

!> This subroutine handles parameter documentation for time-type variables.
subroutine doc_param_time(doc, varname, desc, val, default, units, debuggingParam)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: varname !< The name of the parameter being documented
  character(len=*), intent(in) :: desc    !< A description of the parameter being documented
  type(time_type),  intent(in) :: val     !< The value of the parameter
  type(time_type),  optional, intent(in) :: default !< The default value of this parameter
  character(len=*), optional, intent(in) :: units   !< The units of the parameter being documented
  logical,          optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.

  ! Local varables
  character(len=mLen)              :: mesg          ! The output message
  character(len=doc%commentColumn) :: valstring     ! A string with the formatted value.
  logical                          :: equalsDefault ! True if val = default.

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    valstring = time_string(val)
    if (present(units)) then
      mesg = define_string(doc, varname, valstring, units)
    else
      mesg = define_string(doc, varname, valstring, "[days : seconds]")
    endif

    equalsDefault = .false.
    if (present(default)) then
      if (val == default) equalsDefault = .true.
      mesg = trim(mesg)//" default = "//trim(time_string(default))
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, debuggingParam=debuggingParam)
  endif

end subroutine doc_param_time

!> This subroutine writes out the message and description to the documetation files.
subroutine writeMessageAndDesc(doc, vmesg, desc, valueWasDefault, indent, &
                               layoutParam, debuggingParam)
  type(doc_type),    intent(in) :: doc     !< A pointer to a structure that controls where the
                                           !! documentation occurs and its formatting
  character(len=*),  intent(in) :: vmesg   !< A message with the parameter name, units, and default value.
  character(len=*),  intent(in) :: desc    !< A description of the parameter being documented
  logical, optional, intent(in) :: valueWasDefault !< If true, this parameter has its default value
  integer, optional, intent(in) :: indent      !< An amount by which to indent this message
  logical, optional, intent(in) :: layoutParam !< If present and true, this is a layout parameter.
  logical, optional, intent(in) :: debuggingParam !< If present and true, this is a debugging parameter.

  ! Local variables
  character(len=mLen) :: mesg          ! A full line of a message including indents.
  character(len=mLen) :: mesg_text     ! A line of message text without preliminary indents.
  integer :: start_ind = 1             ! The starting index in the description for the next line.
  integer :: nl_ind, tab_ind, end_ind  ! The indices of new-lines, tabs, and the end of a line.
  integer :: len_text, len_tab, len_nl ! The lengths of the text string, tabs and new-lines.
  integer :: len_cor                   ! The permitted length corrected for tab sizes in a line.
  integer :: len_desc                  ! The non-whitespace length of the description.
  integer :: substr_start              ! The starting index of a substring to search for tabs.
  integer :: indnt, msg_pad            ! Space counts used to format a message.
  logical :: msg_done, reset_msg_pad   ! Logicals used to format messages.
  logical :: all, short, layout, debug ! Flags indicating which files to write into.

  layout = .false. ; if (present(layoutParam)) layout = layoutParam
  debug = .false. ; if (present(debuggingParam)) debug = debuggingParam
  all = doc%complete .and. (doc%unitAll > 0) .and. .not. (layout .or. debug)
  short = doc%minimal .and. (doc%unitShort > 0) .and. .not. (layout .or. debug)
  if (present(valueWasDefault)) short = short .and. (.not. valueWasDefault)

  if (all) write(doc%unitAll, '(a)') trim(vmesg)
  if (short) write(doc%unitShort, '(a)') trim(vmesg)
  if (layout) write(doc%unitLayout, '(a)') trim(vmesg)
  if (debug) write(doc%unitDebugging, '(a)') trim(vmesg)

  if (len_trim(desc) == 0) return

  len_tab = len_trim("_\t_") - 2
  len_nl = len_trim("_\n_") - 2

  indnt = doc%commentColumn ; if (present(indent)) indnt = indent
  len_text = doc%max_line_len - (indnt + 2)
  start_ind = 1 ; msg_pad = 0 ; msg_done = .false.
  do
    if (len_trim(desc(start_ind:)) < 1) exit

    len_cor = len_text - msg_pad

    substr_start = start_ind
    len_desc = len_trim(desc)
    do ! Adjust the available line length for anomalies in the size of tabs, counting \t as 2 spaces.
      if (substr_start >= start_ind+len_cor) exit
      tab_ind = index(desc(substr_start:min(len_desc,start_ind+len_cor)), "\t")
      if (tab_ind == 0) exit
      substr_start = substr_start + tab_ind
      len_cor = len_cor + (len_tab - 2)
    enddo

    nl_ind = index(desc(start_ind:), "\n")
    end_ind = 0
    if ((nl_ind > 0) .and. (len_trim(desc(start_ind:start_ind+nl_ind-2)) > len_cor)) then
      ! This line is too long despite the new-line character.  Look for an earlier space to break.
      end_ind = scan(desc(start_ind:start_ind+len_cor), " ", back=.true.) - 1
      if (end_ind > 0) nl_ind = 0
    elseif ((nl_ind == 0) .and. (len_trim(desc(start_ind:)) > len_cor)) then
      ! This line is too long and does not have a new-line character.  Look for a space to break.
      end_ind = scan(desc(start_ind:start_ind+len_cor), " ", back=.true.) - 1
    endif

    reset_msg_pad = .false.
    if (nl_ind > 0) then
      mesg_text = trim(desc(start_ind:start_ind+nl_ind-2))
      start_ind = start_ind + nl_ind + len_nl - 1
      reset_msg_pad = .true.
    elseif (end_ind > 0) then
      mesg_text = trim(desc(start_ind:start_ind+end_ind))
      start_ind = start_ind + end_ind + 1
      ! Adjust the starting point to move past leading spaces.
      start_ind = start_ind + (len_trim(desc(start_ind:)) - len_trim(adjustl(desc(start_ind:))))
    else
      mesg_text = trim(desc(start_ind:))
      msg_done = .true.
    endif

    do ; tab_ind = index(mesg_text, "\t") ! Replace \t with 2 spaces.
      if (tab_ind == 0) exit
      mesg_text(tab_ind:) = "  "//trim(mesg_text(tab_ind+len_tab:))
    enddo

    mesg = repeat(" ",indnt)//"! "//repeat(" ",msg_pad)//trim(mesg_text)

    if (reset_msg_pad) then
      msg_pad = 0
    elseif (msg_pad == 0) then ! Indent continuation lines.
      msg_pad = len_trim(mesg_text) - len_trim(adjustl(mesg_text))
      ! If already indented, indent an additional 2 spaces.
      if (msg_pad >= 2) msg_pad = msg_pad + 2
    endif

    if (all) write(doc%unitAll, '(a)') trim(mesg)
    if (short) write(doc%unitShort, '(a)') trim(mesg)
    if (layout) write(doc%unitLayout, '(a)') trim(mesg)
    if (debug) write(doc%unitDebugging, '(a)') trim(mesg)

    if (msg_done) exit
  enddo

end subroutine writeMessageAndDesc

! ----------------------------------------------------------------------

!> This function returns a string with a time type formatted as seconds (perhaps including a
!! fractional number of seconds) and days
function time_string(time)
  type(time_type), intent(in) :: time !< The time type being translated
  character(len=40) :: time_string

  ! Local variables
  integer :: secs, days, ticks, ticks_per_sec

  call get_time(Time, secs, days, ticks)

  time_string = trim(adjustl(int_string(days))) // ":" // trim(adjustl(int_string(secs)))
  if (ticks /= 0) then
    ticks_per_sec = get_ticks_per_second()
    time_string = trim(time_string) // ":" // &
                  trim(adjustl(int_string(ticks)))//"/"//trim(adjustl(int_string(ticks_per_sec)))
  endif

end function time_string

!> This function returns a string with a real formatted like '(G)'
function real_string(val)
  real, intent(in)  :: val !< The value being written into a string
  character(len=32) :: real_string
! This function returns a string with a real formatted like '(G)'
  integer :: len, ind

  if ((abs(val) < 1.0e4) .and. (abs(val) >= 1.0e-3)) then
    write(real_string, '(F30.11)') val
    if (.not.testFormattedFloatIsReal(real_string,val)) then
      write(real_string, '(F30.12)') val
      if (.not.testFormattedFloatIsReal(real_string,val)) then
        write(real_string, '(F30.13)') val
        if (.not.testFormattedFloatIsReal(real_string,val)) then
          write(real_string, '(F30.14)') val
          if (.not.testFormattedFloatIsReal(real_string,val)) then
            write(real_string, '(F30.15)') val
            if (.not.testFormattedFloatIsReal(real_string,val)) then
              write(real_string, '(F30.16)') val
            endif
          endif
        endif
      endif
    endif
    do
      len = len_trim(real_string)
      if ((len<2) .or. (real_string(len-1:len) == ".0") .or. &
          (real_string(len:len) /= "0")) exit
      real_string(len:len) = " "
    enddo
  elseif (val == 0.) then
    real_string = "0.0"
  else
    if ((abs(val) <= 1.0e-100) .or. (abs(val) >= 1.0e100)) then
      write(real_string(1:32), '(ES24.14E3)') val
      if (.not.testFormattedFloatIsReal(real_string,val)) &
        write(real_string(1:32), '(ES24.15E3)') val
    else
      write(real_string(1:32), '(ES23.14)') val
      if (.not.testFormattedFloatIsReal(real_string,val)) &
        write(real_string(1:32), '(ES23.15)') val
    endif
    do
      ind = index(real_string,"0E")
      if (ind == 0) exit
      if (real_string(ind-1:ind-1) == ".") exit
      real_string = real_string(1:ind-1)//real_string(ind+1:)
    enddo
  endif
  real_string = adjustl(real_string)
end function real_string

!> Returns a character string of a comma-separated, compact formatted, reals
!> e.g. "1., 2., 5*3., 5.E2", that give the list of values.
function real_array_string(vals, sep)
  character(len=1320)    :: real_array_string !< The output string listing vals
  real,      intent(in)  :: vals(:) !< The array of values to record
  character(len=*), &
    optional, intent(in) :: sep     !< The separator between successive values,
                                    !! by default it is ', '.
! Returns a character string of a comma-separated, compact formatted, reals
! e.g. "1., 2., 5*3., 5.E2"
  ! Local variables
  integer :: j, n, b, ns
  logical :: doWrite
  character(len=10) :: separator
  n=1 ; doWrite=.true. ; real_array_string='' ; b=1
  if (present(sep)) then
    separator=sep ; ns=len(sep)
  else
    separator=', ' ; ns=2
  endif
  do j=1,size(vals)
    doWrite=.true.
    if (j<size(vals)) then
      if (vals(j)==vals(j+1)) then
        n=n+1
        doWrite=.false.
      endif
    endif
    if (doWrite) then
      if (b>1) then ! Write separator if a number has already been written
        write(real_array_string(b:),'(A)') separator
        b=b+ns
      endif
      if (n>1) then
        write(real_array_string(b:),'(A,"*",A)') trim(int_string(n)),trim(real_string(vals(j)))
      else
        write(real_array_string(b:),'(A)') trim(real_string(vals(j)))
      endif
      n=1 ; b=len_trim(real_array_string)+1
    endif
  enddo
end function real_array_string

!> This function tests whether a real value is encoded in a string.
function testFormattedFloatIsReal(str, val)
  character(len=*), intent(in) :: str !< The string that match val
  real,             intent(in) :: val !< The value being tested
  logical                      :: testFormattedFloatIsReal
  ! Local variables
  real :: scannedVal

  read(str(1:),*) scannedVal
  if (scannedVal == val) then
    testFormattedFloatIsReal=.true.
  else
    testFormattedFloatIsReal=.false.
  endif
end function testFormattedFloatIsReal

!> This function returns a string with an integer formatted like '(I)'
function int_string(val)
  integer, intent(in)  :: val !< The value being written into a string
  character(len=24)    :: int_string
! This function returns a string with an integer formatted like '(I)'
  write(int_string, '(i24)') val
  int_string = adjustl(int_string)
end function int_string

!> This function returns a string with an logical formatted like '(L)'
function logical_string(val)
  logical, intent(in)  :: val !< The value being written into a string
  character(len=24)    :: logical_string
! This function returns a string with an logical formatted like '(L)'
  write(logical_string, '(l24)') val
  logical_string = adjustl(logical_string)
end function logical_string

!> This function returns a string for formatted parameter assignment
function define_string(doc, varName, valString, units)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: varName !< The name of the parameter being documented
  character(len=*), intent(in) :: valString !< A string containing the value of the parameter
  character(len=*), intent(in) :: units   !< The units of the parameter being documented
  character(len=mLen) :: define_string
! This function returns a string for formatted parameter assignment
  integer :: numSpaces
  define_string = repeat(" ",mLen) ! Blank everything for safety
  if (doc%defineSyntax) then
    define_string = "#define "//trim(varName)//" "//valString
  else
    define_string = trim(varName)//" = "//valString
  endif
  numSpaces = max(1, doc%commentColumn - len_trim(define_string) )
  define_string = trim(define_string)//repeat(" ",numSpaces)//"!"
  if (len_trim(units) > 0) define_string = trim(define_string)//"   ["//trim(units)//"]"
end function define_string

!> This function returns a string for formatted false logicals
function undef_string(doc, varName, units)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: varName !< The name of the parameter being documented
  character(len=*), intent(in) :: units   !< The units of the parameter being documented
  character(len=mLen) :: undef_string
! This function returns a string for formatted false logicals
  integer :: numSpaces
  undef_string = repeat(" ",240) ! Blank everything for safety
  undef_string = "#undef "//trim(varName)
  if (doc%defineSyntax) then
    undef_string = "#undef "//trim(varName)
  else
    undef_string = trim(varName)//" = "//STRING_FALSE
  endif
  numSpaces = max(1, doc%commentColumn - len_trim(undef_string) )
  undef_string = trim(undef_string)//repeat(" ",numSpaces)//"!"
  if (len_trim(units) > 0) undef_string = trim(undef_string)//"   ["//trim(units)//"]"
end function undef_string

! ----------------------------------------------------------------------

!> This subroutine handles the module documentation
subroutine doc_module(doc, modname, desc)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: modname !< The name of the module being documented
  character(len=*), intent(in) :: desc    !< A description of the module being documented
! This subroutine handles the module documentation
  character(len=mLen) :: mesg

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    call writeMessageAndDesc(doc, '', '') ! Blank line for delineation
    mesg = "! === module "//trim(modname)//" ==="
    call writeMessageAndDesc(doc, mesg, desc, indent=0)
  endif
end subroutine doc_module

!> This subroutine handles the subroutine documentation
subroutine doc_subroutine(doc, modname, subname, desc)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: modname !< The name of the module being documented
  character(len=*), intent(in) :: subname !< The name of the subroutine being documented
  character(len=*), intent(in) :: desc    !< A description of the subroutine being documented
! This subroutine handles the subroutine documentation
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

end subroutine doc_subroutine

!> This subroutine handles the function documentation
subroutine doc_function(doc, modname, fnname, desc)
  type(doc_type),   pointer    :: doc     !< A pointer to a structure that controls where the
                                          !! documentation occurs and its formatting
  character(len=*), intent(in) :: modname !< The name of the module being documented
  character(len=*), intent(in) :: fnname  !< The name of the function being documented
  character(len=*), intent(in) :: desc    !< A description of the function being documented
! This subroutine handles the function documentation
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

end subroutine doc_function

! ----------------------------------------------------------------------

!> Initialize the parameter documentation
subroutine doc_init(docFileBase, doc, minimal, complete, layout, debugging)
  character(len=*),  intent(in)  :: docFileBase !< The base file name for this set of parameters,
                                             !! for example MOM_parameter_doc
  type(doc_type),    pointer     :: doc      !< A pointer to a structure that controls where the
                                             !! documentation occurs and its formatting
  logical, optional, intent(in)  :: minimal  !< If present and true, write out the files (.short) documenting
                                             !! those parameters that do not take on their default values.
  logical, optional, intent(in)  :: complete !< If present and true, write out the (.all) files documenting all
                                             !! parameters
  logical, optional, intent(in)  :: layout   !< If present and true, write out the (.layout) files documenting
                                             !! the layout parameters
  logical, optional, intent(in)  :: debugging !< If present and true, write out the (.debugging) files documenting
                                             !! the debugging parameters

  if (.not. associated(doc)) then
    allocate(doc)
  endif

  doc%docFileBase = docFileBase
  if (present(minimal)) doc%minimal = minimal
  if (present(complete)) doc%complete = complete
  if (present(layout)) doc%layout = layout
  if (present(debugging)) doc%debugging = debugging

end subroutine doc_init

!> This subroutine allocates and populates a structure that controls where the
!! documentation occurs and its formatting, and opens up the files controlled
!! by this structure
subroutine open_doc_file(doc)
  type(doc_type), pointer :: doc !< A pointer to a structure that controls where the
                                 !! documentation occurs and its formatting

  logical :: opened, new_file
  integer :: ios
  character(len=240) :: fileName

  if (.not. (is_root_pe() .and. associated(doc))) return

  if ((len_trim(doc%docFileBase) > 0) .and. doc%complete .and. (doc%unitAll<0)) then
    new_file = .true. ; if (doc%unitAll /= -1) new_file = .false.
    doc%unitAll = find_unused_unit_number()

    write(fileName(1:240),'(a)') trim(doc%docFileBase)//'.all'
    if (new_file) then
      open(doc%unitAll, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='REPLACE', iostat=ios)
      write(doc%unitAll, '(a)') &
       '! This file was written by the model and records all non-layout '//&
       'or debugging parameters used at run-time.'
    else ! This file is being reopened, and should be appended.
      open(doc%unitAll, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='OLD', position='APPEND', iostat=ios)
    endif
    inquire(doc%unitAll, opened=opened)
    if ((.not.opened) .or. (ios /= 0)) then
      call MOM_error(FATAL, "Failed to open doc file "//trim(fileName)//".")
    endif
    doc%filesAreOpen = .true.
  endif

  if ((len_trim(doc%docFileBase) > 0) .and. doc%minimal .and. (doc%unitShort<0)) then
    new_file = .true. ; if (doc%unitShort /= -1) new_file = .false.
    doc%unitShort = find_unused_unit_number()

    write(fileName(1:240),'(a)') trim(doc%docFileBase)//'.short'
    if (new_file) then
      open(doc%unitShort, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='REPLACE', iostat=ios)
      write(doc%unitShort, '(a)') &
       '! This file was written by the model and records the non-default parameters used at run-time.'
    else ! This file is being reopened, and should be appended.
      open(doc%unitShort, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='OLD', position='APPEND', iostat=ios)
    endif
    inquire(doc%unitShort, opened=opened)
    if ((.not.opened) .or. (ios /= 0)) then
      call MOM_error(FATAL, "Failed to open doc file "//trim(fileName)//".")
    endif
    doc%filesAreOpen = .true.
  endif

  if ((len_trim(doc%docFileBase) > 0) .and. doc%layout .and. (doc%unitLayout<0)) then
    new_file = .true. ; if (doc%unitLayout /= -1) new_file = .false.
    doc%unitLayout = find_unused_unit_number()

    write(fileName(1:240),'(a)') trim(doc%docFileBase)//'.layout'
    if (new_file) then
      open(doc%unitLayout, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='REPLACE', iostat=ios)
      write(doc%unitLayout, '(a)') &
       '! This file was written by the model and records the layout parameters used at run-time.'
    else ! This file is being reopened, and should be appended.
      open(doc%unitLayout, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='OLD', position='APPEND', iostat=ios)
    endif
    inquire(doc%unitLayout, opened=opened)
    if ((.not.opened) .or. (ios /= 0)) then
      call MOM_error(FATAL, "Failed to open doc file "//trim(fileName)//".")
    endif
    doc%filesAreOpen = .true.
  endif

  if ((len_trim(doc%docFileBase) > 0) .and. doc%debugging .and. (doc%unitDebugging<0)) then
    new_file = .true. ; if (doc%unitDebugging /= -1) new_file = .false.
    doc%unitDebugging = find_unused_unit_number()

    write(fileName(1:240),'(a)') trim(doc%docFileBase)//'.debugging'
    if (new_file) then
      open(doc%unitDebugging, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='REPLACE', iostat=ios)
      write(doc%unitDebugging, '(a)') &
       '! This file was written by the model and records the debugging parameters used at run-time.'
    else ! This file is being reopened, and should be appended.
      open(doc%unitDebugging, file=trim(fileName), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='OLD', position='APPEND', iostat=ios)
    endif
    inquire(doc%unitDebugging, opened=opened)
    if ((.not.opened) .or. (ios /= 0)) then
      call MOM_error(FATAL, "Failed to open doc file "//trim(fileName)//".")
    endif
    doc%filesAreOpen = .true.
  endif

end subroutine open_doc_file

!> Find an unused unit number, returning >0 if found, and triggering a FATAL error if not.
function find_unused_unit_number()
! Find an unused unit number.
! Returns >0 if found. FATAL if not.
  integer :: find_unused_unit_number
  logical :: opened
  do find_unused_unit_number=512,42,-1
    inquire( find_unused_unit_number, opened=opened)
    if (.not.opened) exit
  enddo
  if (opened) call MOM_error(FATAL, &
    "doc_init failed to find an unused unit number.")
end function find_unused_unit_number

!> This subroutine closes the the files controlled by doc, and sets flags in
!! doc to indicate that parameterization is no longer permitted.
subroutine doc_end(doc)
  type(doc_type), pointer :: doc !< A pointer to a structure that controls where the
                                 !! documentation occurs and its formatting
  type(link_msg), pointer :: this => NULL(), next => NULL()

  if (.not.associated(doc)) return

  if (doc%unitAll > 0) then
    close(doc%unitAll)
    doc%unitAll = -2
  endif

  if (doc%unitShort > 0) then
    close(doc%unitShort)
    doc%unitShort = -2
  endif

  if (doc%unitLayout > 0) then
    close(doc%unitLayout)
    doc%unitLayout = -2
  endif

  if (doc%unitDebugging > 0) then
    close(doc%unitDebugging)
    doc%unitDebugging = -2
  endif

  doc%filesAreOpen = .false.

  this => doc%chain_msg
  do while( associated(this) )
    next => this%next
    deallocate(this)
    this => next
  enddo
end subroutine doc_end

! -----------------------------------------------------------------------------

!> Returns true if documentation has already been written
function mesgHasBeenDocumented(doc,varName,mesg)
  type(doc_type),   pointer     :: doc  !< A pointer to a structure that controls where the
                                        !! documentation occurs and its formatting
  character(len=*), intent(in)  :: varName !< The name of the parameter being documented
  character(len=*), intent(in)  :: mesg !< A message with parameter values, defaults, and descriptions
                                        !! to compare with the message that was written previously
  logical                       :: mesgHasBeenDocumented
! Returns true if documentation has already been written
  type(link_msg), pointer :: newLink => NULL(), this => NULL(), last => NULL()

  mesgHasBeenDocumented = .false.

!!if (mesg(1:1) == '!') return ! Ignore commented parameters

  ! Search through list for this parameter
  last => NULL()
  this => doc%chain_msg
  do while( associated(this) )
    if (trim(doc%blockPrefix)//trim(varName) == trim(this%name)) then
      mesgHasBeenDocumented = .true.
      if (trim(mesg) == trim(this%msg)) return
      ! If we fail the above test then cause an error
      if (mesg(1:1) == '!') return ! Do not cause error for commented parameters
      call MOM_error(WARNING, "Previous msg:"//trim(this%msg))
      call MOM_error(WARNING, "New message :"//trim(mesg))
      call MOM_error(WARNING, "Encountered inconsistent documentation line for parameter "&
                     //trim(varName)//"!")
    endif
    last => this
    this => this%next
  enddo

  ! Allocate a new link
  allocate(newLink)
  newLink%name = trim(doc%blockPrefix)//trim(varName)
  newLink%msg = trim(mesg)
  newLink%next => NULL()
  if (.not. associated(doc%chain_msg)) then
    doc%chain_msg => newLink
  else
    if (.not. associated(last)) call MOM_error(FATAL, &
         "Unassociated LINK in mesgHasBeenDocumented: "//trim(mesg))
    last%next => newLink
  endif
end function mesgHasBeenDocumented

end module MOM_document
