module MOM_string_functions
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
!*  By Robert Hallberg, June 2005.                                     *
!*                                                                     *
!*    The subroutines here parse a set of input files for the value    *
!*  a named parameter and sets that parameter at run time.  Currently  *
!*  these files use the same format as the header file MOM_memory.h.   *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

implicit none ; private

public lowercase, uppercase
public left_int, left_ints
public left_real, left_reals

contains

function lowercase(input_string)
!   This function returns a string in which all uppercase letters have been
! replaced by their lowercase counterparts.  It is loosely based on the
! lowercase function in mpp_util.F90.
  ! Arguments
  character(len=*),     intent(in) :: input_string
  character(len=len(input_string)) :: lowercase
  ! Local variables
  integer, parameter :: co=iachar('a')-iachar('A') ! case offset
  integer :: k

  lowercase = input_string
  do k=1, len_trim(input_string)
    if (lowercase(k:k) >= 'A' .and. lowercase(k:k) <= 'Z') &
        lowercase(k:k) = achar(ichar(lowercase(k:k))+co)
  end do
end function lowercase

function uppercase(input_string)
  character(len=*),     intent(in) :: input_string
  character(len=len(input_string)) :: uppercase
!   This function returns a string in which all lowercase letters have been
! replaced by their uppercase counterparts.  It is loosely based on the
! uppercase function in mpp_util.F90.
  ! Arguments
  integer, parameter :: co=iachar('A')-iachar('a') ! case offset
  integer :: k

  uppercase = input_string
  do k=1, len_trim(input_string)
    if (uppercase(k:k) >= 'a' .and. uppercase(k:k) <= 'z') &
        uppercase(k:k) = achar(ichar(uppercase(k:k))+co)
  end do
end function uppercase

function left_int(i)
! Returns a character string of a left-formatted integer
! e.g. "123       "  (assumes 19 digit maximum)
  ! Arguments
  character(len=19) :: left_int
  integer, intent(in) :: i
  ! Local variables
  character(len=19) :: tmp
  write(tmp(1:19),'(I19)') i
  write(left_int(1:19),'(A)') adjustl(tmp)
end function left_int

function left_ints(i)
! Returns a character string of a comma-separated, compact formatted,
! integers  e.g. "1, 2, 3, 4"
  ! Arguments
  character(len=1320) :: left_ints
  integer, intent(in) :: i(:)
  ! Local variables
  character(len=1320) :: tmp
  integer :: j
  write(left_ints(1:1320),'(A)') trim(left_int(i(1)))
  if (size(i)>1) then
    do j=2,size(i)
      tmp=left_ints
      write(left_ints(1:1320),'(A,", ",A)') trim(tmp),trim(left_int(i(j)))
    enddo
  endif
end function left_ints

function left_real(r)
! Returns a character string of a left-formatted real using either
! F or E format with trailing 0s removed e.g. "12.345    "
  ! Arguments
  character(len=30) :: left_real
  real, intent(in) :: r
  ! Local variables
  character(len=30) :: tmp,expon
  integer :: i

  write(tmp(1:30),*) r," "
  i=index(tmp,'E')
  ! Save exponent
  if (i>0) then
    write(expon(1:30),'(A)') trim(tmp(i:len_trim(tmp)))
    tmp(i:len_trim(tmp))=' '
  else
    expon=' '
  endif
  ! Left justify
  write(left_real(1:30),'(A)') adjustl(tmp)
  ! Strip trailing zeros
  do i=len_trim(left_real),2,-1
    if (left_real(i:i).ne.'0') exit
    left_real(i:i)=' '
  enddo
  ! Re-attach exponent
  write(left_real(i+1:len(left_real)),'(A)') trim(expon)
end function left_real

function left_reals(r)
! Returns a character string of a comma-separated, compact formatted, reals
! e.g. "1., 2., 3., 5.E2"
  ! Arguments
  character(len=1320) :: left_reals
  real, intent(in) :: r(:)
  ! Local variables
  character(len=1320) :: tmp
  integer :: j
  write(left_reals(1:1320),'(A)') trim(left_real(r(1)))
  if (size(r)>1) then
    do j=2,size(r)
      tmp=left_reals
      write(left_reals(1:1320),'(A,", ",A)') trim(tmp),trim(left_real(r(j)))
    enddo
  endif
end function left_reals

end module MOM_string_functions
