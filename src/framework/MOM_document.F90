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

interface doc_param
  module procedure doc_param_none, &
                   doc_param_logical, doc_param_logical_array, &
                   doc_param_int,     doc_param_int_array, &
                   doc_param_real,    doc_param_real_array, &
                   doc_param_char, &
                   doc_param_time
end interface

type :: link_logical ; private
  type(link_logical), pointer :: next => NULL() ! Facilitates linked list
  character(len=80) :: name                        ! Parameter name
  logical :: default                            ! Parameter value and default
  logical :: hasAdefaultValue = .false.         ! Has a default value
  logical :: hasBeenDocumented = .false.        ! Has been written to documentation
end type link_logical

type :: link_int ; private
  type(link_int), pointer :: next => NULL()  ! Facilitates linked list
  character(len=80) :: name                     ! Parameter name
  integer :: default                         ! Parameter value and default
  logical :: hasAdefaultValue = .false.      ! Has a default value
  logical :: hasBeenDocumented = .false.     ! Has been written to documentation
end type link_int

type :: link_real ; private
  type(link_real), pointer :: next => NULL() ! Facilitates linked list
  character(len=80) :: name                     ! Parameter name
  real :: default                            ! Parameter value and default
  logical :: hasAdefaultValue = .false.      ! Has a default value
  logical :: hasBeenDocumented = .false.     ! Has been written to documentation
end type link_real

type :: link_char ; private
  type(link_char), pointer :: next => NULL() ! Facilitates linked list
  character(len=80) :: name                  ! Parameter name
  character(len=120) :: default              ! Parameter value and default
  logical :: hasAdefaultValue = .false.      ! Has a default value
  logical :: hasBeenDocumented = .false.     ! Has been written to documentation
end type link_char

type, public :: doc_type ; private
  integer :: unit = -1              ! The open unit number for docfile.
  character(len=240) :: docfile = '' ! A file where all run-time parameters,
                                    ! their settings and defaults are documented.
  logical :: minimal = .false.      ! If true, only document non-default values.
  logical :: defineSyntax = .true.  ! If true, use the a=b syntax
  logical :: warnOnConflicts = .false. ! Cause a WARNING error if defaults differ.
  integer :: commentColumn = 32     ! Number of spaces before the comment marker.
  type(link_logical), pointer :: chain_logicals => NULL() ! Db of logicals
  type(link_int),     pointer :: chain_ints => NULL()     ! Db of integers
  type(link_real),    pointer :: chain_reals => NULL()    ! Db of reals
  type(link_char),    pointer :: chain_chars => NULL()    ! Db of characters
end type doc_type

character(len=4), parameter :: STRING_TRUE = 'True'
character(len=5), parameter :: STRING_FALSE = 'False'

contains

! ----------------------------------------------------------------------

subroutine doc_param_none(doc, varname, desc, units)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
! This subroutine handles parameter documentation with no value.
  integer :: numspc
  character(len=240) :: mesg

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (is_root_pe() .and. (doc%unit > 0)) then
    numspc = max(1,doc%commentColumn-8-len_trim(varname))
    mesg = "#define "//trim(varname)//repeat(" ",numspc)//"!"
    if (len_trim(units) > 0) mesg = trim(mesg)//"   ["//trim(units)//"]"

    write(doc%unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(doc, desc)
  endif 
end subroutine doc_param_none

subroutine doc_param_logical(doc, varname, desc, units, val, default)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  logical,          intent(in) :: val
  logical,          optional, intent(in) :: default
! This subroutine handles parameter documentation for logicals.
  character(len=240) :: mesg

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (is_root_pe() .and. (doc%unit > 0)) then
    if (present(default) .and. doc%minimal) then
      if (val .eqv. default) return
    endif

    if (val) then
      mesg = define_string(doc,varname,STRING_TRUE,units)
    else
      mesg = undef_string(doc,varname,units)
    endif

    if (present(default)) then
      if (default) then
        mesg = trim(mesg)//" default = "//STRING_TRUE
      else
        mesg = trim(mesg)//" default = "//STRING_FALSE
      endif
    endif

    write(doc%unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(doc, desc)
  endif
end subroutine doc_param_logical

subroutine doc_param_logical_array(doc, varname, desc, units, vals, default)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  logical,          intent(in) :: vals(:)
  logical,          optional, intent(in) :: default
! This subroutine handles parameter documentation for arrays of logicals.
  logical :: go_back
  integer :: i
  character(len=1280) :: mesg
  character(len=1240)  :: valstring

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (is_root_pe() .and. (doc%unit > 0)) then
    ! If documentation has already been written and has the
    ! same defaults then skip re-writing it.
    if (logicalHasBeenDocumented(doc,varname,default)) return

    if (present(default) .and. doc%minimal) then
      go_back = .true.
      do i=1,size(vals) ; if (vals(i) .neqv. default) go_back = .false. ; enddo
      if (go_back) return
    endif

    if (vals(1)) then ; valstring = STRING_TRUE ; else ; valstring = STRING_FALSE ; endif
    do i=2,min(size(vals),128)
      if (vals(i)) then
        valstring = trim(valstring)//", "//STRING_TRUE
      else
        valstring = trim(valstring)//", "//STRING_FALSE
      endif
    enddo

    mesg = define_string(doc,varname,valstring,units)

    if (present(default)) then
      if (default) then
        mesg = trim(mesg)//" default = "//STRING_TRUE
      else
        mesg = trim(mesg)//" default = "//STRING_FALSE
      endif
    endif

    write(doc%unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(doc, desc)
  endif
end subroutine doc_param_logical_array

subroutine doc_param_int(doc, varname, desc, units, val, default)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  integer,          intent(in) :: val
  integer,          optional, intent(in) :: default
! This subroutine handles parameter documentation for integers.
  character(len=240) :: mesg
  character(len=doc%commentColumn)  :: valstring

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (is_root_pe() .and. (doc%unit > 0)) then
    ! If documentation has already been written and has the
    ! same defaults then skip re-writing it.
    if (intHasBeenDocumented(doc,varname,default)) return

    if (present(default) .and. doc%minimal) then
      if (val == default) return
    endif

    valstring = int_string(val)
    mesg = define_string(doc,varname,valstring,units)

    if (present(default)) &
      mesg = trim(mesg)//" default = "//(trim(int_string(default)))

    write(doc%unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(doc, desc)
  endif
end subroutine doc_param_int

subroutine doc_param_int_array(doc, varname, desc, units, vals, default)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  integer,          intent(in) :: vals(:)
  integer,          optional, intent(in) :: default
! This subroutine handles parameter documentation for arrays of integers.
  logical :: go_back
  integer :: i
  character(len=1280) :: mesg
  character(len=1240)  :: valstring

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (is_root_pe() .and. (doc%unit > 0)) then
    if (present(default) .and. doc%minimal) then
      go_back = .true.
      do i=1,size(vals) ; if (vals(i) /= default) go_back = .false. ; enddo
      if (go_back) return
    endif

    valstring = int_string(vals(1))
    do i=2,min(size(vals),128)
      valstring = trim(valstring)//", "//trim(int_string(vals(i)))
    enddo

    mesg = define_string(doc,varname,valstring,units)

    if (present(default)) &
      mesg = trim(mesg)//" default = "//(trim(int_string(default)))

    write(doc%unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(doc, desc)
  endif

end subroutine doc_param_int_array

subroutine doc_param_real(doc, varname, desc, units, val, default)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  real,             intent(in) :: val
  real,             optional, intent(in) :: default
! This subroutine handles parameter documentation for reals.
  character(len=240) :: mesg
  character(len=doc%commentColumn)  :: valstring

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (is_root_pe() .and. (doc%unit > 0)) then
    ! If documentation has already been written and has the
    ! same defaults then skip re-writing it.
    if (realHasBeenDocumented(doc,varname,default)) return
    
    if (present(default) .and. doc%minimal) then
      if (val == default) return
    endif

    valstring = real_string(val)
    mesg = define_string(doc,varname,valstring,units)

    if (present(default)) &
      mesg = trim(mesg)//" default = "//trim(real_string(default))

    write(doc%unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(doc, desc)
  endif
end subroutine doc_param_real

subroutine doc_param_real_array(doc, varname, desc, units, vals, default)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  real,             intent(in) :: vals(:)
  real,             optional, intent(in) :: default
! This subroutine handles parameter documentation for arrays of reals.
  logical :: go_back
  integer :: i
  character(len=1280) :: mesg
  character(len=1240)  :: valstring

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (is_root_pe() .and. (doc%unit > 0)) then
    if (present(default) .and. doc%minimal) then
      go_back = .true.
      do i=1,size(vals) ; if (vals(i) /= default) go_back = .false. ; enddo
      if (go_back) return
    endif

    valstring = real_string(vals(1))
    do i=2,min(size(vals),128)
      valstring = trim(valstring)//", "//trim(real_string(vals(i)))
    enddo

    mesg = define_string(doc,varname,valstring,units)

    if (present(default)) &
      mesg = trim(mesg)//" default = "//trim(real_string(default))

    write(doc%unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(doc, desc)
  endif

end subroutine doc_param_real_array

subroutine doc_param_char(doc, varname, desc, units, val, default)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  character(len=*), intent(in) :: val
  character(len=*), optional, intent(in) :: default
! This subroutine handles parameter documentation for character strings.
  character(len=240) :: mesg

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (is_root_pe() .and. (doc%unit > 0)) then
    ! If documentation has already been written and has the
    ! same defaults then skip re-writing it.
    if (charHasBeenDocumented(doc,varname,default)) return
    
    if (present(default) .and. doc%minimal) then
      if (trim(val) == trim(default)) return
    endif

    mesg = define_string(doc,varname,'"'//trim(val)//'"',units)

    if (present(default)) &
      mesg = trim(mesg)//' default = "'//trim(adjustl(default))//'"'

    write(doc%unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(doc, desc)
  endif

end subroutine doc_param_char

subroutine doc_param_time(doc, varname, desc, units, val, default)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varname, desc, units
  type(time_type),  intent(in) :: val
  type(time_type),  optional, intent(in) :: default
! This subroutine handles parameter documentation for time-type variables.
!  ### This needs to be written properly!
  integer :: numspc
  character(len=240) :: mesg

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (is_root_pe() .and. (doc%unit > 0)) then
    numspc = max(1,doc%commentColumn-18-len_trim(varname))
    mesg = "#define "//trim(varname)//" Time-type"//repeat(" ",numspc)//"!"
    if (len_trim(units) > 0) mesg = trim(mesg)//"   ["//trim(units)//"]"

    write(doc%unit, '(a)') trim(mesg)
    
    if (len_trim(desc) > 0) call write_desc(doc, desc)
  endif 

end subroutine doc_param_time

subroutine write_desc(doc, desc, indent)
  type(doc_type),             intent(in) :: doc
  character(len=*),           intent(in) :: desc
  integer,          optional, intent(in) :: indent
  character(len=240) :: mesg
  integer :: start_ind = 1, end_ind, indnt, tab, len_tab, len_nl

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
      write(doc%unit, '(a)') trim(mesg)
    else
      mesg = repeat(" ",indnt)//"! "//trim(desc(start_ind:))
      do ; tab = index(mesg, "\t")
        if (tab == 0) exit
        mesg(tab:) = "  "//trim(mesg(tab+len_tab:))
      enddo
      write(doc%unit, '(a)') trim(mesg)
      exit
    endif

  enddo
end subroutine write_desc

! ----------------------------------------------------------------------

function real_string(val)
  real, intent(in)  :: val
  character(len=32) :: real_string
! This function returns a string with a real formatted like '(G)'
  integer :: len, ind

  if ((abs(val) < 1.0e4) .and. (abs(val) >= 1.0e-3)) then
    write(real_string, '(F30.12)') val
    do
      len = len_trim(real_string)
      if ((len<2) .or. (real_string(len-1:len) == ".0") .or. &
          (real_string(len:len) /= "0")) exit
      real_string(len:len) = " "
    enddo
  elseif (val == 0) then
    real_string = "0.0"
  else
    write(real_string(1:32), '(ES23.15)') val
    do
      ind = index(real_string,"0E")
      if (ind == 0) exit
      if (real_string(ind-1:ind-1) == ".") exit
      real_string = real_string(1:ind-1)//real_string(ind+1:)
    enddo
  endif
  real_string = adjustl(real_string)
  
end function real_string

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
  character(len=240) :: define_string
! This function returns a string for formatted parameter assignment
  integer :: numSpaces
  define_string = repeat(" ",240) ! Blank everything for safety
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
  character(len=240) :: undef_string
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
  character(len=240) :: mesg

  if (.not. (is_root_pe() .and. associated(doc))) return
  call open_doc_file(doc)

  if (is_root_pe() .and. (doc%unit > 0)) then
    mesg = "    !  Parameters of module "//trim(modname)
    write(doc%unit, '(a)') trim(mesg)
    if (len_trim(desc) > 0) call write_desc(doc, desc, 8)
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
  
subroutine doc_init(docfile, doc, minimal)
  character(len=*),  intent(in)  :: docfile
  type(doc_type),    pointer     :: doc
  logical, optional, intent(in)  :: minimal
! Arguments: docfile - The name of the doc file.
!  (inout)   doc - The doc_type to populate.

  if (.not. associated(doc)) then
    allocate(doc)
  endif

  doc%docfile = docfile
  if (present(minimal)) doc%minimal = minimal
end subroutine doc_init

subroutine open_doc_file(doc)
  type(doc_type), pointer    :: doc

  logical :: opened, new_file
  integer :: ios, unit

  if (.not. (is_root_pe() .and. associated(doc))) return

  if (is_root_pe() .and. (len_trim(doc%docfile) > 0) .and. (doc%unit<0)) then
    new_file = .true. ; if (doc%unit /= -1) new_file = .false.
    ! Find an unused unit number.
    do unit=512,42,-1
      inquire( unit, opened=opened)
      if (.not.opened) exit
    enddo

    if (opened) call MOM_error(FATAL, &
        "doc_init failed to find an unused unit number.")
    doc%unit = unit

    if (new_file) then
      open(doc%unit, file=trim(doc%docfile), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='REPLACE', iostat=ios)
    else ! This file is being reopened, and should be appended.
      open(doc%unit, file=trim(doc%docfile), access='SEQUENTIAL', form='FORMATTED', &
           action='WRITE', status='OLD', position='APPEND', iostat=ios)
    endif
    inquire(doc%unit, opened=opened)
    if ((.not.opened) .or. (ios /= 0)) then
      call MOM_error(FATAL, "Failed to open doc file "//trim(doc%docfile)//".")
    endif
  endif

end subroutine open_doc_file

subroutine doc_end(doc)
  type(doc_type), pointer :: doc

  if (.not.associated(doc)) return

  if (doc%unit > 0) then
    close(doc%unit)
    doc%unit = -2
  endif

end subroutine doc_end

! -----------------------------------------------------------------------------

function logicalHasBeenDocumented(doc,varName,defaultValue)
  type(doc_type),   pointer     :: doc
  character(len=*), intent(in)  :: varName
  logical, optional, intent(in) :: defaultValue
  logical                       :: logicalHasBeenDocumented
! Returns true is documentation has been written and current arguments are
! consistent with previous documentation
  type(link_logical), pointer :: newRealLink, this

  allocate(newRealLink)
  newRealLink%name = varName
  if (present(defaultValue)) then
    newRealLink%default = defaultValue
    newRealLink%hasAdefaultValue = .true.
  else
    newRealLink%hasAdefaultValue = .false.
  endif
  newRealLink%next => NULL()

  ! Now search through list for this parameter
  logicalHasBeenDocumented = .false.
  if (.not. associated(doc%chain_logicals)) then
    newRealLink%hasBeenDocumented = .true.
    doc%chain_logicals => newRealLink
  else
    this => doc%chain_logicals
    do while( associated(this) )
      if (trim(varName) == trim(this%name)) then
        deallocate(newRealLink)
        if (this%hasAdefaultValue) then
          if (present(defaultValue)) then
            if (defaultValue .eqv. this%default) then
             logicalHasBeenDocumented = .true.
            else
             if (doc%warnOnConflicts) call MOM_error(WARNING,   &
               'Conflicting default: '//trim(varName)//' '//     &
               trim(logical_string(defaultValue))//' '//trim(logical_string(this%default)))
            endif
          endif
        endif
        return
      endif
      this => this%next
    enddo
    newRealLink%next => doc%chain_logicals ! If not encountered in the linked list, add it
    doc%chain_logicals => newRealLink
  endif
end function logicalHasBeenDocumented

function intHasBeenDocumented(doc,varName,defaultValue)
  type(doc_type),   pointer       :: doc
  character(len=*), intent(in)    :: varName
  integer, optional,   intent(in) :: defaultValue
  logical                         :: intHasBeenDocumented
! Returns true is documentation has been written and current arguments are
! consistent with previous documentation
  type(link_int), pointer :: newRealLink, this

  allocate(newRealLink)
  newRealLink%name = varName
  if (present(defaultValue)) then
    newRealLink%default = defaultValue
    newRealLink%hasAdefaultValue = .true.
  else
    newRealLink%hasAdefaultValue = .false.
  endif
  newRealLink%next => NULL()

  ! Now search through list for this parameter
  intHasBeenDocumented = .false.
  if (.not. associated(doc%chain_ints)) then
    newRealLink%hasBeenDocumented = .true.
    doc%chain_ints => newRealLink
  else
    this => doc%chain_ints
    do while( associated(this) )
      if (trim(varName) == trim(this%name)) then
        deallocate(newRealLink)
        if (this%hasAdefaultValue) then
          if (present(defaultValue)) then
            if (defaultValue==this%default) then
             intHasBeenDocumented = .true.
            else
             if (doc%warnOnConflicts) call MOM_error(WARNING,   &
               'Conflicting default: '//trim(varName)//' '//     &
               trim(int_string(defaultValue))//' '//trim(int_string(this%default)))
            endif
          endif
        endif
        return
      endif
      this => this%next
    enddo
    newRealLink%next => doc%chain_ints ! If not encountered in the linked list, add it
    doc%chain_ints => newRealLink
  endif
end function intHasBeenDocumented

function realHasBeenDocumented(doc,varName,defaultValue)
  type(doc_type),   pointer    :: doc
  character(len=*), intent(in) :: varName
  real, optional,   intent(in) :: defaultValue
  logical                      :: realHasBeenDocumented
! Returns true is documentation has been written and current arguments are
! consistent with previous documentation
  type(link_real), pointer :: newRealLink, this

  allocate(newRealLink)
  newRealLink%name = varName
  if (present(defaultValue)) then
    newRealLink%default = defaultValue
    newRealLink%hasAdefaultValue = .true.
  else
    newRealLink%hasAdefaultValue = .false.
  endif
  newRealLink%next => NULL()

  ! Now search through list for this parameter
  realHasBeenDocumented = .false.
  if (.not. associated(doc%chain_reals)) then
    newRealLink%hasBeenDocumented = .true.
    doc%chain_reals => newRealLink
  else
    this => doc%chain_reals
    do while( associated(this) )
      if (trim(varName) == trim(this%name)) then
        deallocate(newRealLink)
        if (this%hasAdefaultValue) then
          if (present(defaultValue)) then
            if (defaultValue==this%default) then
             realHasBeenDocumented = .true.
            else
             if (doc%warnOnConflicts) call MOM_error(WARNING,   &
               'Conflicting default: '//trim(varName)//' '//     &
               trim(real_string(defaultValue))//' '//trim(real_string(this%default)))
            endif
          endif
        endif
        return
      endif
      this => this%next
    enddo
    newRealLink%next => doc%chain_reals ! If not encountered in the linked list, add it
    doc%chain_reals => newRealLink
  endif
end function realHasBeenDocumented

function charHasBeenDocumented(doc,varName,defaultValue)
  type(doc_type),             pointer    :: doc
  character(len=*),           intent(in) :: varName
  character(len=*), optional, intent(in) :: defaultValue
  logical                                :: charHasBeenDocumented
! Returns true is documentation has been written and current arguments are
! consistent with previous documentation
  type(link_char), pointer :: newRealLink, this

  allocate(newRealLink)
  newRealLink%name = varName
  if (present(defaultValue)) then
    newRealLink%default = defaultValue
    newRealLink%hasAdefaultValue = .true.
  else
    newRealLink%hasAdefaultValue = .false.
  endif
  newRealLink%next => NULL()

  ! Now search through list for this parameter
  charHasBeenDocumented = .false.
  if (.not. associated(doc%chain_chars)) then
    newRealLink%hasBeenDocumented = .true.
    doc%chain_chars => newRealLink
  else
    this => doc%chain_chars
    do while( associated(this) )
      if (trim(varName) == trim(this%name)) then
        deallocate(newRealLink)
        if (this%hasAdefaultValue) then
          if (present(defaultValue)) then
            if (trim(defaultValue)==trim(this%default)) then
             charHasBeenDocumented = .true.
            else
             if (doc%warnOnConflicts) call MOM_error(WARNING,    &
               'Conflicting default: '//trim(varName)//' "'//     &
               trim(defaultValue)//'" "'//trim(this%default)//'"')
            endif
          endif
        endif
        return
      endif
      this => this%next
    enddo
    newRealLink%next => doc%chain_chars ! If not encountered in the linked list, add it
    doc%chain_chars => newRealLink
  endif
end function charHasBeenDocumented

end module MOM_document
