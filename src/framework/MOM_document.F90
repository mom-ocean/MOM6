module MOM_document
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
!*    The subroutines here provide hooks for document generation       *
!*  functions at various levels of granularity.                        *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_time_manager, only : time_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe

implicit none ; private

public doc_param, doc_subroutine, doc_function, doc_module, doc_init, doc_end
public doc_openBlock, doc_closeBlock

interface doc_param
  module procedure doc_param_none, &
                   doc_param_logical, doc_param_logical_array, &
                   doc_param_int,     doc_param_int_array, &
                   doc_param_real,    doc_param_real_array, &
                   doc_param_char, &
                   doc_param_time
end interface

integer, parameter :: mLen = 1240 ! Length of interface/message strings

type, public :: doc_type ; private
  integer :: unitAll = -1           ! The open unit number for docFileBase + .all.
  integer :: unitShort = -1         ! The open unit number for docFileBase + .short.
  integer :: unitLayout = -1        ! The open unit number for docFileBase + .layout.
  logical :: filesAreOpen = .false. ! True if any files were successfully opened.
  character(len=mLen) :: docFileBase = '' ! The basename of the files where run-time
                                    ! parameters, settings and defaults are documented.
  logical :: complete = .true.      ! If true, document all parameters.
  logical :: minimal = .true.       ! If true, document non-default parameters.
  logical :: layout = .true.        ! If true, document layout parameters.
  logical :: defineSyntax = .false. ! If true, use #def syntax instead of a=b syntax
  logical :: warnOnConflicts = .false. ! Cause a WARNING error if defaults differ.
  integer :: commentColumn = 32     ! Number of spaces before the comment marker.
  type(link_msg), pointer :: chain_msg => NULL() ! Db of messages
  character(len=240) :: blockPrefix = '' ! The full name of the current block.
end type doc_type

type :: link_msg ; private
  type(link_msg), pointer :: next => NULL()  ! Facilitates linked list
  character(len=80) :: name                  ! Parameter name
  character(len=620) :: msg                  ! Parameter value and default
end type link_msg

character(len=4), parameter :: STRING_TRUE = 'True'
character(len=5), parameter :: STRING_FALSE = 'False'

contains

! ----------------------------------------------------------------------

subroutine doc_param_none(doc, varname, desc, units)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
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

subroutine doc_param_logical(doc, varname, desc, units, val, default, layoutParam)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  logical,          intent(in) :: val
  logical,          optional, intent(in) :: default
  logical,          optional, intent(in) :: layoutParam
! This subroutine handles parameter documentation for logicals.
  character(len=mLen) :: mesg
  logical :: equalsDefault

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    if (val) then
      mesg = define_string(doc,varname,STRING_TRUE,units)
    else
      mesg = undef_string(doc,varname,units)
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
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, layoutParam=layoutParam)
  endif
end subroutine doc_param_logical

subroutine doc_param_logical_array(doc, varname, desc, units, vals, default, layoutParam)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  logical,          intent(in) :: vals(:)
  logical,          optional, intent(in) :: default
  logical,          optional, intent(in) :: layoutParam
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

    mesg = define_string(doc,varname,valstring,units)

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
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, layoutParam=layoutParam)
  endif
end subroutine doc_param_logical_array

subroutine doc_param_int(doc, varname, desc, units, val, default, layoutParam)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  integer,          intent(in) :: val
  integer,          optional, intent(in) :: default
  logical,          optional, intent(in) :: layoutParam
! This subroutine handles parameter documentation for integers.
  character(len=mLen) :: mesg
  character(len=doc%commentColumn)  :: valstring
  logical :: equalsDefault

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    valstring = int_string(val)
    mesg = define_string(doc,varname,valstring,units)

    equalsDefault = .false.
    if (present(default)) then
      if (val == default) equalsDefault = .true.
      mesg = trim(mesg)//" default = "//(trim(int_string(default)))
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, layoutParam=layoutParam)
  endif
end subroutine doc_param_int

subroutine doc_param_int_array(doc, varname, desc, units, vals, default, layoutParam)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  integer,          intent(in) :: vals(:)
  integer,          optional, intent(in) :: default
  logical,          optional, intent(in) :: layoutParam
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

    mesg = define_string(doc,varname,valstring,units)

    equalsDefault = .false.
    if (present(default)) then
      equalsDefault = .true.
      do i=1,size(vals) ; if (vals(i) /= default) equalsDefault = .false. ; enddo
      mesg = trim(mesg)//" default = "//(trim(int_string(default)))
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, layoutParam=layoutParam)
  endif

end subroutine doc_param_int_array

subroutine doc_param_real(doc, varname, desc, units, val, default)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  real,             intent(in) :: val
  real,             optional, intent(in) :: default
! This subroutine handles parameter documentation for reals.
  character(len=mLen) :: mesg
  character(len=doc%commentColumn)  :: valstring
  logical :: equalsDefault

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    valstring = real_string(val)
    mesg = define_string(doc,varname,valstring,units)

    equalsDefault = .false.
    if (present(default)) then
      if (val == default) equalsDefault = .true.
      mesg = trim(mesg)//" default = "//trim(real_string(default))
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault)
  endif
end subroutine doc_param_real

subroutine doc_param_real_array(doc, varname, desc, units, vals, default)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  real,             intent(in) :: vals(:)
  real,             optional, intent(in) :: default
! This subroutine handles parameter documentation for arrays of reals.
  integer :: i
  character(len=mLen) :: mesg
  character(len=mLen) :: valstring
  logical :: equalsDefault

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    valstring = trim(real_array_string(vals(:)))

    mesg = define_string(doc,varname,valstring,units)

    equalsDefault = .false.
    if (present(default)) then
      equalsDefault = .true.
      do i=1,size(vals) ; if (vals(i) /= default) equalsDefault = .false. ; enddo
      mesg = trim(mesg)//" default = "//trim(real_string(default))
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault)
  endif

end subroutine doc_param_real_array

subroutine doc_param_char(doc, varname, desc, units, val, default, layoutParam)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  character(len=*), intent(in) :: val
  character(len=*), optional, intent(in) :: default
  logical,          optional, intent(in) :: layoutParam
! This subroutine handles parameter documentation for character strings.
  character(len=mLen) :: mesg
  logical :: equalsDefault

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (doc%filesAreOpen) then
    mesg = define_string(doc,varname,'"'//trim(val)//'"',units)

    equalsDefault = .false.
    if (present(default)) then
      if (trim(val) == trim(default)) equalsDefault = .true.
      mesg = trim(mesg)//' default = "'//trim(adjustl(default))//'"'
    endif

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault, layoutParam=layoutParam)
  endif

end subroutine doc_param_char

subroutine doc_openBlock(doc, blockName, desc)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: blockName
  character(len=*), optional, intent(in) :: desc
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

subroutine doc_closeBlock(doc, blockName)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: blockName
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

subroutine doc_param_time(doc, varname, desc, units, val, default)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  type(time_type),  intent(in) :: val
  type(time_type),  optional, intent(in) :: default
! This subroutine handles parameter documentation for time-type variables.
!  ### This needs to be written properly!
  integer :: numspc
  character(len=mLen) :: mesg
  logical :: equalsDefault

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  equalsDefault = .false.
  if (doc%filesAreOpen) then
    numspc = max(1,doc%commentColumn-18-len_trim(varname))
    mesg = "#define "//trim(varname)//" Time-type"//repeat(" ",numspc)//"!"
    if (len_trim(units) > 0) mesg = trim(mesg)//"   ["//trim(units)//"]"

    if (mesgHasBeenDocumented(doc, varName, mesg)) return ! Avoid duplicates
    call writeMessageAndDesc(doc, mesg, desc, equalsDefault)
  endif

end subroutine doc_param_time

subroutine writeMessageAndDesc(doc, vmesg, desc, valueWasDefault, indent, layoutParam)
  type(doc_type),             intent(in) :: doc
  character(len=*),           intent(in) :: vmesg, desc
  logical,          optional, intent(in) :: valueWasDefault
  integer,          optional, intent(in) :: indent
  logical,          optional, intent(in) :: layoutParam
  character(len=mLen) :: mesg
  integer :: start_ind = 1, end_ind, indnt, tab, len_tab, len_nl
  logical :: all, short, layout

  layout = .false.
  if (present(layoutParam)) layout = layoutParam
  all = doc%complete .and. (doc%unitAll > 0) .and. .not. layout
  short = doc%minimal .and. (doc%unitShort > 0) .and. .not. layout
  if (present(valueWasDefault)) short = short .and. (.not. valueWasDefault)

  if (all) write(doc%unitAll, '(a)') trim(vmesg)
  if (short) write(doc%unitShort, '(a)') trim(vmesg)
  if (layout) write(doc%unitLayout, '(a)') trim(vmesg)

  if (len_trim(desc) == 0) return

  len_tab = len_trim("_\t_") - 2
  len_nl = len_trim("_\n_") -2

  indnt = doc%commentColumn ; if (present(indent)) indnt = indent
  start_ind = 1
  do
    if (len_trim(desc(start_ind:)) < 1) exit

    end_ind = index(desc(start_ind:), "\n")

    if (end_ind > 0) then
      mesg = repeat(" ",indnt)//"! "//trim(desc(start_ind:start_ind+end_ind-2))
      start_ind = start_ind + end_ind - 1 + len_nl

      do ; tab = index(mesg, "\t")
        if (tab == 0) exit
        mesg(tab:) = "  "//trim(mesg(tab+len_tab:))
      enddo
      if (all) write(doc%unitAll, '(a)') trim(mesg)
      if (short) write(doc%unitShort, '(a)') trim(mesg)
      if (layout) write(doc%unitLayout, '(a)') trim(mesg)
    else
      mesg = repeat(" ",indnt)//"! "//trim(desc(start_ind:))
      do ; tab = index(mesg, "\t")
        if (tab == 0) exit
        mesg(tab:) = "  "//trim(mesg(tab+len_tab:))
      enddo
      if (all) write(doc%unitAll, '(a)') trim(mesg)
      if (short) write(doc%unitShort, '(a)') trim(mesg)
      if (layout) write(doc%unitLayout, '(a)') trim(mesg)
      exit
    endif

  enddo
end subroutine writeMessageAndDesc

! ----------------------------------------------------------------------

function real_string(val)
  real, intent(in)  :: val
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
    write(real_string(1:32), '(ES23.14)') val
    if (.not.testFormattedFloatIsReal(real_string,val)) then
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

function real_array_string(vals,sep)
  character(len=1320) :: real_array_string
  real, intent(in)  :: vals(:)
  character(len=*), optional :: sep
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

function testFormattedFloatIsReal(str, val)
  character(len=*), intent(in) :: str
  real,             intent(in) :: val
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

function int_string(val)
  integer, intent(in)  :: val
  character(len=24)    :: int_string
! This function returns a string with an integer formatted like '(I)'
  write(int_string, '(i24)') val
  int_string = adjustl(int_string)
end function int_string

function logical_string(val)
  logical, intent(in)  :: val
  character(len=24)    :: logical_string
! This function returns a string with an logical formatted like '(L)'
  write(logical_string, '(l24)') val
  logical_string = adjustl(logical_string)
end function logical_string

function define_string(doc,varName,valString,units)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varName, valString, units
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

function undef_string(doc,varName,units)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varName, units
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

subroutine doc_module(doc, modname, desc)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: modname, desc
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

subroutine doc_subroutine(doc, modname, subname, desc)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: modname, subname, desc
! This subroutine handles the subroutine documentation
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

end subroutine doc_subroutine

subroutine doc_function(doc, modname, fnname, desc)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: modname, fnname, desc
! This subroutine handles the function documentation
  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

end subroutine doc_function

! ----------------------------------------------------------------------

subroutine doc_init(docFileBase, doc, minimal, complete)
  character(len=*),  intent(in)  :: docFileBase
  type(doc_type),    pointer     :: doc
  logical, optional, intent(in)  :: minimal, complete
! Arguments: docFileBase - The name of the doc file.
!  (inout)   doc - The doc_type to populate.

  if (.not. associated(doc)) then
    allocate(doc)
  endif

  doc%docFileBase = docFileBase
  if (present(minimal)) doc%minimal = minimal
  if (present(minimal)) doc%complete = complete
end subroutine doc_init

subroutine open_doc_file(doc)
  type(doc_type), pointer    :: doc

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
       '! This file was written by the model and records all non-layout parameters used at run-time.'
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

end subroutine open_doc_file

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

subroutine doc_end(doc)
  type(doc_type), pointer :: doc
  type(link_msg), pointer :: this, next

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

  doc%filesAreOpen = .false.

  this => doc%chain_msg
  do while( associated(this) )
    next => this%next
    deallocate(this)
    this => next
  enddo
end subroutine doc_end

! -----------------------------------------------------------------------------

function mesgHasBeenDocumented(doc,varName,mesg)
  type(doc_type),   pointer     :: doc
  character(len=*), intent(in)  :: varName, mesg
  logical                       :: mesgHasBeenDocumented
! Returns true if documentation has already been written
  type(link_msg), pointer :: newLink, this, last

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
