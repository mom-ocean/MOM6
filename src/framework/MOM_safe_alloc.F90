module MOM_safe_alloc

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
!*    The subroutines here provide a convenient way to safely allocate *
!*  memory without accidentally reallocating a pointer and causing a   *
!*  memory leak.                                                       *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

implicit none ; private

public safe_alloc_ptr, safe_alloc_alloc

interface safe_alloc_ptr
  module procedure safe_alloc_ptr_3d_2arg, safe_alloc_ptr_2d_2arg
  module procedure safe_alloc_ptr_3d, safe_alloc_ptr_2d, safe_alloc_ptr_1d
end interface safe_alloc_ptr

interface safe_alloc_alloc
  module procedure safe_alloc_allocatable_3d, safe_alloc_allocatable_2d
end interface safe_alloc_alloc

!   This combined interface might work with a later version of Fortran, but
! it fails with the gnu F90 compiler.
!
! interface safe_alloc
!   module procedure safe_alloc_ptr_3d_2arg, safe_alloc_ptr_2d_2arg
!   module procedure safe_alloc_ptr_3d, safe_alloc_ptr_2d, safe_alloc_ptr_1d
!   module procedure safe_alloc_allocatable_3d, safe_alloc_allocatable_2d
! end interface safe_alloc

contains

subroutine safe_alloc_ptr_1d(ptr, i1, i2)
  real, pointer :: ptr(:)
  integer, intent(in) :: i1
  integer, optional, intent(in) :: i2
  if (.not.ASSOCIATED(ptr)) then
    if (present(i2)) then
      allocate(ptr(i1:i2))
    else
      allocate(ptr(i1))
    endif
    ptr(:) = 0.0
  endif
end subroutine safe_alloc_ptr_1d

subroutine safe_alloc_ptr_2d_2arg(ptr, ni, nj)
  real, pointer :: ptr(:,:)
  integer, intent(in) :: ni, nj
  if (.not.ASSOCIATED(ptr)) then
    allocate(ptr(ni,nj))
    ptr(:,:) = 0.0
  endif
end subroutine safe_alloc_ptr_2d_2arg

subroutine safe_alloc_ptr_3d_2arg(ptr, ni, nj, nk)
  real, pointer :: ptr(:,:,:)
  integer, intent(in) :: ni, nj, nk
  if (.not.ASSOCIATED(ptr)) then
    allocate(ptr(ni,nj,nk))
    ptr(:,:,:) = 0.0
  endif
end subroutine safe_alloc_ptr_3d_2arg

subroutine safe_alloc_ptr_2d(ptr, is, ie, js, je)
  real, pointer :: ptr(:,:)
  integer, intent(in) :: is, ie, js, je
  if (.not.ASSOCIATED(ptr)) then
    allocate(ptr(is:ie,js:je))
    ptr(:,:) = 0.0
  endif
end subroutine safe_alloc_ptr_2d

subroutine safe_alloc_ptr_3d(ptr, is, ie, js, je, nk)
  real, pointer :: ptr(:,:,:)
  integer, intent(in) :: is, ie, js, je, nk
  if (.not.ASSOCIATED(ptr)) then
    allocate(ptr(is:ie,js:je,nk))
    ptr(:,:,:) = 0.0
  endif
end subroutine safe_alloc_ptr_3d

subroutine safe_alloc_allocatable_2d(ptr, is, ie, js, je)
  real, allocatable :: ptr(:,:)
  integer, intent(in) :: is, ie, js, je
  if (.not.ALLOCATED(ptr)) then
    allocate(ptr(is:ie,js:je))
    ptr(:,:) = 0.0
  endif
end subroutine safe_alloc_allocatable_2d

subroutine safe_alloc_allocatable_3d(ptr, is, ie, js, je, nk)
  real, allocatable :: ptr(:,:,:)
  integer, intent(in) :: is, ie, js, je, nk
  if (.not.ALLOCATED(ptr)) then
    allocate(ptr(is:ie,js:je,nk))
    ptr(:,:,:) = 0.0
  endif
end subroutine safe_alloc_allocatable_3d

end module MOM_safe_alloc
