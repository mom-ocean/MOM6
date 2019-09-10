!> Convenience functions for safely allocating memory without
!! accidentally reallocating pointer and causing memory leaks.
module MOM_safe_alloc

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public safe_alloc_ptr, safe_alloc_alloc

!> Allocate a pointer to a 1-d, 2-d or 3-d array
interface safe_alloc_ptr
  module procedure safe_alloc_ptr_3d_3arg,  safe_alloc_ptr_3d_6arg, safe_alloc_ptr_2d_2arg
  module procedure safe_alloc_ptr_3d, safe_alloc_ptr_2d, safe_alloc_ptr_1d
end interface safe_alloc_ptr

!> Allocate a 2-d or 3-d allocatable array
interface safe_alloc_alloc
  module procedure safe_alloc_allocatable_3d, safe_alloc_allocatable_2d
  module procedure safe_alloc_allocatable_3d_6arg
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

!> Allocate a pointer to a 1-d array
subroutine safe_alloc_ptr_1d(ptr, i1, i2)
  real, dimension(:), pointer :: ptr !< A pointer to allocate
  integer,            intent(in) :: i1 !< The size of the array, or its starting index if i2 is present
  integer, optional,  intent(in) :: i2 !< The ending index of the array
  if (.not.associated(ptr)) then
    if (present(i2)) then
      allocate(ptr(i1:i2))
    else
      allocate(ptr(i1))
    endif
    ptr(:) = 0.0
  endif
end subroutine safe_alloc_ptr_1d

!> Allocate a pointer to a 2-d array based on its dimension sizes
subroutine safe_alloc_ptr_2d_2arg(ptr, ni, nj)
  real, dimension(:,:), pointer :: ptr !< A pointer to allocate
  integer, intent(in) :: ni !< The size of the 1st dimension of the array
  integer, intent(in) :: nj !< The size of the 2nd dimension of the array
  if (.not.associated(ptr)) then
    allocate(ptr(ni,nj))
    ptr(:,:) = 0.0
  endif
end subroutine safe_alloc_ptr_2d_2arg

!> Allocate a pointer to a 3-d array based on its dimension sizes
subroutine safe_alloc_ptr_3d_3arg(ptr, ni, nj, nk)
  real, dimension(:,:,:), pointer :: ptr !< A pointer to allocate
  integer, intent(in) :: ni !< The size of the 1st dimension of the array
  integer, intent(in) :: nj !< The size of the 2nd dimension of the array
  integer, intent(in) :: nk !< The size of the 3rd dimension of the array
  if (.not.associated(ptr)) then
    allocate(ptr(ni,nj,nk))
    ptr(:,:,:) = 0.0
  endif
end subroutine safe_alloc_ptr_3d_3arg

!> Allocate a pointer to a 2-d array based on its index starting and ending values
subroutine safe_alloc_ptr_2d(ptr, is, ie, js, je)
  real, dimension(:,:), pointer :: ptr !< A pointer to allocate
  integer, intent(in) :: is !< The start index to allocate for the 1st dimension
  integer, intent(in) :: ie !< The end index to allocate for the 1st dimension
  integer, intent(in) :: js !< The start index to allocate for the 2nd dimension
  integer, intent(in) :: je !< The end index to allocate for the 2nd dimension
  if (.not.associated(ptr)) then
    allocate(ptr(is:ie,js:je))
    ptr(:,:) = 0.0
  endif
end subroutine safe_alloc_ptr_2d

!> Allocate a pointer to a 3-d array based on its index starting and ending values
subroutine safe_alloc_ptr_3d(ptr, is, ie, js, je, nk)
  real, dimension(:,:,:), pointer :: ptr !< A pointer to allocate
  integer, intent(in) :: is !< The start index to allocate for the 1st dimension
  integer, intent(in) :: ie !< The end index to allocate for the 1st dimension
  integer, intent(in) :: js !< The start index to allocate for the 2nd dimension
  integer, intent(in) :: je !< The end index to allocate for the 2nd dimension
  integer, intent(in) :: nk !< The size to allocate for the 3rd dimension
  if (.not.associated(ptr)) then
    allocate(ptr(is:ie,js:je,nk))
    ptr(:,:,:) = 0.0
  endif
end subroutine safe_alloc_ptr_3d

!> Allocate a pointer to a 3-d array based on its index starting and ending values
subroutine safe_alloc_ptr_3d_6arg(ptr, is, ie, js, je, ks, ke)
  real, dimension(:,:,:), pointer :: ptr !< A pointer to allocate
  integer, intent(in) :: is !< The start index to allocate for the 1st dimension
  integer, intent(in) :: ie !< The end index to allocate for the 1st dimension
  integer, intent(in) :: js !< The start index to allocate for the 2nd dimension
  integer, intent(in) :: je !< The end index to allocate for the 2nd dimension
  integer, intent(in) :: ks !< The start index to allocate for the 3rd dimension
  integer, intent(in) :: ke !< The end index to allocate for the 3rd dimension
  if (.not.associated(ptr)) then
    allocate(ptr(is:ie,js:je,ks:ke))
    ptr(:,:,:) = 0.0
  endif
end subroutine safe_alloc_ptr_3d_6arg


!> Allocate a 2-d allocatable array based on its index starting and ending values
subroutine safe_alloc_allocatable_2d(ptr, is, ie, js, je)
  real, dimension(:,:), allocatable :: ptr !< An allocatable array to allocate
  integer, intent(in) :: is !< The start index to allocate for the 1st dimension
  integer, intent(in) :: ie !< The end index to allocate for the 1st dimension
  integer, intent(in) :: js !< The start index to allocate for the 2nd dimension
  integer, intent(in) :: je !< The end index to allocate for the 2nd dimension
  if (.not.allocated(ptr)) then
    allocate(ptr(is:ie,js:je))
    ptr(:,:) = 0.0
  endif
end subroutine safe_alloc_allocatable_2d

!> Allocate a 3-d allocatable array based on its index starting and ending values
!! and k-index size
subroutine safe_alloc_allocatable_3d(ptr, is, ie, js, je, nk)
  real, dimension(:,:,:), allocatable :: ptr !< An allocatable array to allocate
  integer, intent(in) :: is !< The start index to allocate for the 1st dimension
  integer, intent(in) :: ie !< The end index to allocate for the 1st dimension
  integer, intent(in) :: js !< The start index to allocate for the 2nd dimension
  integer, intent(in) :: je !< The end index to allocate for the 2nd dimension
  integer, intent(in) :: nk !< The size to allocate for the 3rd dimension
  if (.not.allocated(ptr)) then
    allocate(ptr(is:ie,js:je,nk))
    ptr(:,:,:) = 0.0
  endif
end subroutine safe_alloc_allocatable_3d

!> Allocate a 3-d allocatable array based on its 6 index starting and ending values
subroutine safe_alloc_allocatable_3d_6arg(ptr, is, ie, js, je, ks, ke)
  real, dimension(:,:,:), allocatable :: ptr !< An allocatable array to allocate
  integer, intent(in) :: is !< The start index to allocate for the 1st dimension
  integer, intent(in) :: ie !< The end index to allocate for the 1st dimension
  integer, intent(in) :: js !< The start index to allocate for the 2nd dimension
  integer, intent(in) :: je !< The end index to allocate for the 2nd dimension
  integer, intent(in) :: ks !< The start index to allocate for the 3rd dimension
  integer, intent(in) :: ke !< The end index to allocate for the 3rd dimension
  if (.not.allocated(ptr)) then
    allocate(ptr(is:ie,js:je,ks:ke))
    ptr(:,:,:) = 0.0
  endif
end subroutine safe_alloc_allocatable_3d_6arg

end module MOM_safe_alloc
