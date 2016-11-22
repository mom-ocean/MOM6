module MOM_transform_test

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

use MOM_coms, only : PE_here, root_PE
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : log_version, get_param, param_file_type
use MOM_error_handler,  only : callTree_enter, callTree_leave

use mpp_domains_mod, only : mpp_get_compute_domain, mpp_set_compute_domain, mpp_copy_domain, domain2d
use mpp_mod, only : mpp_gather, mpp_max
use ensemble_manager_mod, only : get_ensemble_size, get_ensemble_id, get_ensemble_pelist

implicit none ; private

public :: MOM_transform_test_init, transform_test_started
public :: transform, transform_and_swap
public :: do_transform_test, do_transform_on_this_pe
public :: transform_compare, undo_transform
public :: transform_pointer, undo_transform_pointer
public :: transform_allocatable, undo_transform_allocatable
public :: transform_domain
public :: swap_pointer

interface swap_pointer
  module procedure swap_pointer_2d, swap_pointer_3d
end interface

interface transform
  module procedure transform_2d, transform_3d
end interface

interface transform_allocatable
  module procedure transform_allocatable_2d, transform_allocatable_3d, &
                   transform_allocatable_4d
end interface

interface undo_transform_allocatable
  module procedure transform_allocatable_2d
end interface

interface transform_pointer
  module procedure transform_pointer_2d_ptr_log
  module procedure transform_pointer_2d_ptr, transform_pointer_3d_ptr
end interface

interface undo_transform_pointer
  module procedure transform_pointer_2d_ptr_log
  module procedure transform_pointer_2d_ptr, transform_pointer_3d_ptr
end interface

interface undo_transform
  module procedure undo_transform_2d, undo_transform_3d
end interface

interface transform_compare
  module procedure transform_compare_scalar
  module procedure transform_compare_1d, transform_compare_2d, transform_compare_3d
end interface

interface transform_and_swap
  module procedure transform_and_swap_2d, transform_and_swap_3d
end interface

integer, parameter :: TRANSFORM_TRANSPOSE = 1
integer, parameter :: TRANSFORM_ROT90 = 2

!> Whether or not we're in a transform test run
logical :: transform_test = .false.
!> Whether the transform being done on this PE?
logical :: transform_on_this_pe = .false.
!> Whether the test has started. No comparisons are done before this
! flag is set.
logical :: test_started = .false.

integer :: transform_type = TRANSFORM_TRANSPOSE

contains

! =====================================================================

!> MOM_transform_test_init initializes the module.
subroutine MOM_transform_test_init(param_file)
  type(param_file_type),   intent(in)    :: param_file
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_transform_test" ! This module's name.

  integer, dimension(6) :: ensemble_size
  character(len=9) :: trans_type_str
  integer, dimension(:, :), allocatable :: ensemble_pelist

  if (test_started) then
    return
  endif

  call log_version(param_file, mod, version)

  call get_param(param_file, mod, "TRANSFORM_TEST", &
                 trans_type_str, &
                 "Whether or not to run a transformation \n"//&
                 "test. This involves transposing or rotating all \n"//&
                 "model inputs. This is a testing feature that can be \n"//&
                 "used to help find horizontal indexing errors. \n"//&
                 "This can be either 'TRANSPOSE', ROT90 or NONE' \n", &
                 default="NONE")

  if (trim(trans_type_str) == 'TRANSPOSE') then
    transform_type = TRANSFORM_TRANSPOSE
    transform_test = .true.
  elseif (trim(trans_type_str) == 'ROT90') then
    transform_type = TRANSFORM_ROT90
    transform_test = .true.
  elseif (trim(trans_type_str) == 'NONE') then
    transform_test = .false.
  else
    call MOM_error(FATAL, &
         "TRANSFORM_TEST: TRANSFORM_TYPE must be 'TRANSPOSE', 'ROT90' or 'NONE'")
  endif

  ! Check that we're running as an ensemble. And that each ensemble member uses 1 PE.
  if (transform_test) then
    ensemble_size = get_ensemble_size()
    if (ensemble_size(1) == 1) then
      call MOM_error(FATAL, &
                     "TRANSFORM_TEST: Model not within ensemble")
    endif
    if (ensemble_size(2) /= 1) then
      call MOM_error(FATAL, &
                     "TRANSFORM_TEST: must have 1 PE per ensemble")
    endif

    allocate(ensemble_pelist(ensemble_size(1), ensemble_size(2)))
    call get_ensemble_pelist(ensemble_pelist)

    ! For this test the root PE will be the transformed run and the other
    ! will be the vanilla run.
    if (PE_here() == ensemble_pelist(1, 1)) then
      transform_on_this_pe = .true.
    endif

    deallocate(ensemble_pelist)

  endif

  test_started = transform_test

end subroutine MOM_transform_test_init

function transform_test_started()
    logical :: transform_test_started

    transform_test_started = test_started

end function transform_test_started

function do_transform_on_this_pe()
    logical :: do_transform_on_this_pe

    do_transform_on_this_pe = transform_on_this_pe

end function do_transform_on_this_pe

function do_transform_test()
    logical :: do_transform_test

    do_transform_test = transform_test

end function do_transform_test

subroutine swap_pointer_2d(arrayA, arrayB)
  real, dimension(:,:), pointer, intent(inout) :: arrayA
  real, dimension(:,:), pointer, intent(inout) :: arrayB

  real, dimension(:,:), allocatable :: tmp

  allocate(tmp(size(arrayA, 1), size(arrayA, 2)))
  tmp(:,:) = arrayA(:,:)
  arrayA(:,:) = arrayB(:,:)
  arrayB(:,:) = tmp(:,:)
  deallocate(tmp)

end subroutine swap_pointer_2d

subroutine swap_pointer_3d(arrayA, arrayB)
  real, dimension(:,:,:), pointer, intent(inout) :: arrayA
  real, dimension(:,:,:), pointer, intent(inout) :: arrayB

  real, dimension(:,:,:), allocatable :: tmp

  allocate(tmp(size(arrayA, 1), size(arrayA, 2), size(arrayA, 3)))
  tmp(:,:,:) = arrayA(:,:,:)
  arrayA(:,:,:) = arrayB(:,:,:)
  arrayB(:,:,:) = tmp(:,:,:)
  deallocate(tmp)

end subroutine swap_pointer_3d

subroutine transform_domain(Domain)
  type(domain2d), intent(inout) :: Domain

  integer :: xbegin, xend, ybegin, yend

  call mpp_get_compute_domain(Domain, xbegin, xend, ybegin, yend)
  call mpp_set_compute_domain(Domain, ybegin, yend, xbegin, xend)

end subroutine transform_domain

!< 
subroutine transform_2d(arrayIn, arrayOut)
  real, dimension(:,:), intent(in) :: arrayIn !< Input array
  real, dimension(:,:), intent(out) :: arrayOut !< Transformed input

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_2d: should not be called on this PE.')
  endif

  if (transform_type == TRANSFORM_TRANSPOSE) then
      call transpose_2d(arrayIn, arrayOut)
  else
      ! 90 degree rotation
      call rot90_2d(arrayIn, arrayOut, 1)
  endif

end subroutine transform_2d

subroutine transform_3d(arrayIn, arrayOut)
  real, dimension(:,:,:), intent(in) :: arrayIn !< Input array
  real, dimension(:,:,:), intent(out) :: arrayOut !< Transformed input

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_3d: should not be called on this PE.')
  endif

  if (transform_type == TRANSFORM_TRANSPOSE) then
      call transpose_3d(arrayIn, arrayOut)
  else
      call rot90_3d(arrayIn, arrayOut, 1)
  endif

end subroutine transform_3d

!< Transform an allocatable array. After this call input may have
! a different shape.
subroutine transform_allocatable_2d(array)
  real, dimension(:,:), allocatable, intent(inout) :: array

  real, allocatable, dimension(:,:) :: tmp
  integer :: isz, jsz

  if (.not. allocated(array)) then
    call MOM_error(FATAL, 'transform_2d: array not allocated.')
  endif

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_2d: should not be called on this PE.')
  endif

  isz = size(array, 1)
  jsz = size(array, 2)

  allocate(tmp(isz, jsz))
  tmp(:, :) = array(:, :)
  deallocate(array)
  allocate(array(jsz, isz))

  if (transform_type == TRANSFORM_TRANSPOSE) then
    call transpose_2d(tmp, array)
  else
    call rot90_2d(tmp, array, 1)
  endif

  deallocate(tmp)

end subroutine transform_allocatable_2d

subroutine transform_pointer_2d_ptr(array)
  real, dimension(:,:), pointer, contiguous, intent(inout) :: array

  real, allocatable, dimension(:,:) :: tmp
  integer :: is, ie, js, je

  if (.not. associated(array)) then
    call MOM_error(FATAL, 'transform_2d_ptr: array not associated.')
  endif

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_2d: should not be called on this PE.')
  endif

  is = lbound(array, 1)
  ie = ubound(array, 1)
  js = lbound(array, 2)
  je = ubound(array, 2)

  allocate(tmp(js:je, is:ie))
  if (transform_type == TRANSFORM_TRANSPOSE) then
    call transpose_2d(array, tmp)
  else
    call rot90_2d(array, tmp, 1)
  endif

  array(js:je, is:ie) => array
  array(:, :) = tmp(:, :)

  deallocate(tmp)

end subroutine transform_pointer_2d_ptr

subroutine transform_pointer_2d_ptr_log(array)
  logical, dimension(:,:), pointer, contiguous, intent(inout) :: array

  logical, allocatable, dimension(:,:) :: tmp
  integer :: is, ie, js, je

  if (.not. associated(array)) then
    call MOM_error(FATAL, 'transform_2d_ptr: array not associated.')
  endif

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_2d: should not be called on this PE.')
  endif

  is = lbound(array, 1)
  ie = ubound(array, 1)
  js = lbound(array, 2)
  je = ubound(array, 2)

  allocate(tmp(js:je, is:ie))
  if (transform_type == TRANSFORM_TRANSPOSE) then
    tmp = transpose(array)
  else
    call MOM_error(FATAL, 'rot90_2d: not supported for logical arguments.')
  endif

  array(js:je, is:ie) => array
  array(:, :) = tmp(:, :)

  deallocate(tmp)

end subroutine transform_pointer_2d_ptr_log

subroutine transform_pointer_3d_ptr(array)
  real, dimension(:,:, :), pointer, contiguous, intent(inout) :: array

  real, allocatable, dimension(:,:,:) :: tmp
  integer :: is, ie, js, je, ks, ke

  if (.not. associated(array)) then
    call MOM_error(FATAL, 'transform_2d_ptr: array not associated.')
  endif

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_2d: should not be called on this PE.')
  endif

  is = lbound(array, 1)
  ie = ubound(array, 1)
  js = lbound(array, 2)
  je = ubound(array, 2)
  ks = lbound(array, 3)
  ke = ubound(array, 3)

  allocate(tmp(js:je, is:ie, ks:ke))
  if (transform_type == TRANSFORM_TRANSPOSE) then
    call transpose_3d(array, tmp)
  else
    call rot90_3d(array, tmp, 1)
  endif
  array(js:je, is:ie, ks:ke) => array

  array(:, :, :) = tmp(:, :, :)

  deallocate(tmp)

end subroutine transform_pointer_3d_ptr

!< Transform an allocatable array. After this call input may have
! a different shape.
subroutine transform_allocatable_3d(array)
  real, dimension(:, :, :), allocatable, intent(inout) :: array

  real, allocatable, dimension(:, :, :) :: tmp
  integer :: isz, jsz, ksz

  if (.not. allocated(array)) then
    call MOM_error(FATAL, 'transform_3d: array not associated.')
  endif

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_3d: should not be called on this PE.')
  endif

  isz = size(array, 1)
  jsz = size(array, 2)
  ksz = size(array, 3)

  allocate(tmp(isz, jsz, ksz))
  tmp(:, :, :) = array(:, :, :)
  deallocate(array)
  allocate(array(jsz, isz, ksz))

  if (transform_type == TRANSFORM_TRANSPOSE) then
    call transpose_3d(tmp, array)
  else
    call rot90_3d(tmp, array, 1)
  endif

  deallocate(tmp)

end subroutine transform_allocatable_3d

!< Transform an allocatable array. After this call input may have
! a different shape.
subroutine transform_allocatable_3d_ptr(array)
  real, dimension(:, :, :), pointer, intent(inout) :: array

  real, allocatable, dimension(:, :, :) :: tmp
  integer :: isz, jsz, ksz

  if (.not. associated(array)) then
    call MOM_error(FATAL, 'transform_3d_ptr: array not associated.')
  endif

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_3d: should not be called on this PE.')
  endif

  isz = size(array, 1)
  jsz = size(array, 2)
  ksz = size(array, 3)

  allocate(tmp(isz, jsz, ksz))
  tmp(:, :, :) = array(:, :, :)
  deallocate(array)
  allocate(array(jsz, isz, ksz))

  if (transform_type == TRANSFORM_TRANSPOSE) then
    call transpose_3d(tmp, array)
  else
    call rot90_3d(tmp, array, 1)
  endif

  deallocate(tmp)

end subroutine transform_allocatable_3d_ptr


!< Transform an allocatable array. After this call input may have
! a different shape.
subroutine transform_allocatable_4d(array)
  real, dimension(:, :, :, :), allocatable, intent(inout) :: array

  real, allocatable, dimension(:, :, :, :) :: tmp
  integer :: isz, jsz, ksz, lsz

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_4d: should not be called on this PE.')
  endif

  isz = size(array, 1)
  jsz = size(array, 2)
  ksz = size(array, 3)
  lsz = size(array, 4)

  allocate(tmp(isz, jsz, ksz, lsz))
  tmp(:, :, :, :) = array(:, :, :, :)
  deallocate(array)
  allocate(array(jsz, isz, ksz, lsz))

  if (transform_type == TRANSFORM_TRANSPOSE) then
    call transpose_4d(tmp, array)
  else
    call MOM_error(FATAL, 'transform_4d: NotImplemented.')
  endif

  deallocate(tmp)

end subroutine transform_allocatable_4d

subroutine undo_transform_2d(original, undone)
  real, dimension(:,:), intent(in) :: original  !< The transformed array
  real, dimension(:,:), intent(out) :: undone !< The un-transformed array

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'undo_transform_2d: should not be called on this PE.')
  endif

  if (transform_type == TRANSFORM_TRANSPOSE) then
      call transpose_2d(original, undone)
  else
      call rot90_2d(original, undone, 3)
  endif

end subroutine undo_transform_2d


subroutine undo_transform_3d(original, reversed)
  real, dimension(:,:,:), intent(in) :: original  !< The transformed array
  real, dimension(:,:,:), intent(out) :: reversed !< The un-transformed array

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'undo_transform_3d: should not be called on this PE.')
  endif

  if (transform_type == TRANSFORM_TRANSPOSE) then
      call transpose_3d(original, reversed)
  else
      call rot90_3d(original, reversed, 3)
  endif

end subroutine undo_transform_3d


subroutine transpose_2d(arrayIn, arrayOut)
  real, dimension(:,:), intent(in) :: arrayIn !< Array to be transposed
  real, dimension(:,:), intent(out) :: arrayOut !< Transposed array

  if (size(arrayIn, 1) /= size(arrayOut, 2) &
      .or. size(arrayIn, 2) /= size(arrayOut, 1)) then
    call MOM_error(FATAL, 'transform_2d: array shapes incompatible.')
  endif

  arrayOut(:, :) = transpose(arrayIn(:, :))

end subroutine transpose_2d

subroutine transpose_3d(arrayIn, arrayOut)
  real, dimension(:,:,:), intent(in) :: arrayIn !< The array to be transposed
  real, dimension(:,:,:), intent(inout) :: arrayOut !< Transposed array

  integer :: k

  do k=lbound(arrayIn, 3), ubound(arrayIn, 3)
     call transpose_2d(arrayIn(:, :, k), arrayOut(:, :, k))
  enddo

end subroutine transpose_3d

subroutine transpose_4d(arrayIn, arrayOut)
  real, dimension(:,:,:,:), intent(in) :: arrayIn !< The array to be transposed
  real, dimension(:,:,:,:), intent(inout) :: arrayOut !< Transposed array

  integer :: l

  do l=lbound(arrayIn, 4), ubound(arrayIn, 4)
     call transpose_3d(arrayIn(:, :, : , l), arrayOut(:, :, :, l))
  enddo

end subroutine transpose_4d

subroutine transform_and_swap_2d(arrayA, arrayB)
  real, intent(inout), dimension(:,:) :: arrayA, arrayB

  real, allocatable, dimension(:,:) :: tmp

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_and_swap_2d: should not be called on this PE.')
  endif

  if (size(arrayA, 1) /= size(arrayB, 1) .or. size(arrayA, 2) /= size(arrayB, 2)) then
    call MOM_error(FATAL, 'transform_and_swap_2d: arrays shapes incompatible.')
  endif

  allocate(tmp(size(arrayA, 1), size(arrayA, 2)))

  tmp(:, :) = arrayA(:, :)

  if (transform_type == TRANSFORM_TRANSPOSE) then
    call transpose_2d(arrayB, arrayA)
    call transpose_2d(tmp, arrayB)
  else
    call rot90_2d(arrayB, arrayA, 1)
    call rot90_2d(tmp, arrayB, 1)
  endif

  deallocate(tmp)

end subroutine transform_and_swap_2d

subroutine transform_and_swap_3d(arrayA, arrayB)
  real, intent(inout), dimension(:,:,:) :: arrayA, arrayB

  real, allocatable, dimension(:,:,:) :: tmp

  if (.not. transform_on_this_pe) then
    call MOM_error(FATAL, 'transform_and_swap_3d: should not be called on this PE.')
  endif

  if (size(arrayA, 1) /= size(arrayB, 2) .or. &
          size(arrayA, 2) /= size(arrayB, 1) .or. &
          size(arrayA, 3) /= size(arrayB, 3)) then
    call MOM_error(FATAL, 'trans_and_swap_3d: array shapes incompatible.')
  endif

  allocate(tmp(size(arrayA, 1), size(arrayA, 2), size(arrayA, 3)))

  tmp(:, :, :) = arrayA(:, :, :)

  if (transform_type == TRANSFORM_TRANSPOSE) then
    call transpose_3d(arrayB, arrayA)
    call transpose_3d(tmp, arrayB)
  else
    call rot90_3d(arrayB, arrayA, 1)
    call rot90_3d(tmp, arrayB, 1)
  endif

  deallocate(tmp)

end subroutine transform_and_swap_3d

subroutine transform_compare_2d(arrayA, arrayB, ret)
  real, dimension(:, :), intent(in) :: arrayA, arrayB
  integer, intent(out) :: ret
  integer :: retA, retB

  real, allocatable, dimension(:,:) :: tmp

  ret = 1
  if (.not. test_started) then
    ret = 0
    return
  endif

  if (transform_on_this_pe) then
    allocate(tmp(size(arrayA, 2), size(arrayA, 1)))
    call undo_transform_2d(arrayA, tmp)
    call ensemble_compare_1d(reshape(tmp, (/ size(tmp) /)), retA)

    if (retA /= 0) then
      call write_to_netcdf_2d(tmp, 'transform_test_A_trans.nc')
    endif
    deallocate(tmp)

    allocate(tmp(size(arrayB, 2), size(arrayB, 1)))
    call undo_transform_2d(arrayB, tmp)
    call ensemble_compare_1d(reshape(tmp, (/ size(tmp) /)), retB)

    if (retB /= 0) then
      call write_to_netcdf_2d(tmp, 'transform_test_B_trans.nc')
    endif
    deallocate(tmp)
  else
    call ensemble_compare_1d(reshape(arrayB, (/ size(arrayB) /)), retA)

    if (retA /= 0) then
      call write_to_netcdf_2d(arrayB, 'transform_test_B_vanilla.nc')
    endif

    call ensemble_compare_1d(reshape(arrayA, (/ size(arrayA) /)), retB)

    if (ret /= 0) then
      call write_to_netcdf_2d(arrayB, 'transform_test_A_vanilla.nc')
    endif
  endif

  ret = retA + retB

end subroutine transform_compare_2d

subroutine transform_compare_3d(arrayA, arrayB, ret)
  real, dimension(:, :, :), intent(in) :: arrayA, arrayB
  integer, intent(out) :: ret

  real, allocatable, dimension(:, :, :) :: tmp
  integer :: retA, retB

  if (.not. test_started) then
    ret = 0
    return
  endif

  if (transform_on_this_pe) then
    allocate(tmp(size(arrayA, 2), size(arrayA, 1), size(arrayA, 3)))
    call undo_transform_3d(arrayA, tmp)
    call ensemble_compare_1d(reshape(tmp, (/ size(tmp) /)), retA)

    if (retA /= 0) then
      call write_to_netcdf_3d(tmp, 'transform_test_A_trans.nc')
    endif
    deallocate(tmp)

    allocate(tmp(size(arrayB, 2), size(arrayB, 1), size(arrayB, 3)))
    call undo_transform_3d(arrayB, tmp)
    call ensemble_compare_1d(reshape(tmp, (/ size(tmp) /)), retB)

    if (retB /= 0) then
      call write_to_netcdf_3d(tmp, 'transform_test_B_trans.nc')
    endif
    deallocate(tmp)

  else

    call ensemble_compare_1d(reshape(arrayB, (/ size(arrayB) /)), retA)

    if (retA /= 0) then
      call write_to_netcdf_3d(arrayB, 'transform_test_B_vanilla.nc')
    endif

    call ensemble_compare_1d(reshape(arrayA, (/ size(arrayA) /)), retB)

    if (retB /= 0) then
      call write_to_netcdf_3d(arrayB, 'transform_test_A_vanilla.nc')
    endif
  endif

  ret = retA + retB

end subroutine transform_compare_3d

subroutine transform_compare_1d(arrayA, arrayB, ret)
  real, intent(in), dimension(:) :: arrayA, arrayB
  integer, intent(out) :: ret

  if (.not. test_started) then
    ret = 0
    return
  endif

  if (transform_on_this_pe) then
    call ensemble_compare_1d(arrayA, ret)
  else
    call ensemble_compare_1d(arrayB, ret)
  endif

end subroutine transform_compare_1d

subroutine transform_compare_scalar(a, b, ret)
  integer, intent(in) :: a, b
  integer, intent(out) :: ret

  real, dimension(1) :: tmp

  if (.not. test_started) then
    ret = 0
    return
  endif

  if (transform_on_this_pe) then
    tmp(1) = real(a)
  else
    tmp(1) = real(b)
  endif

  call ensemble_compare_1d(tmp, ret)

end subroutine transform_compare_scalar

subroutine ensemble_compare_1d(sbuf, ret)
  real, intent(in), dimension(:) :: sbuf
  integer, intent(out) :: ret

  integer, dimension(:, :), allocatable :: ensemble_pelist
  integer, dimension(6) :: ensemble_size
  real, dimension(:), allocatable :: rbuf
  integer :: sbuf_size, e, i
  real :: a, b

  ret = 0

  ! If we are running in an ensemble then communicate with other member.
  ensemble_size = get_ensemble_size()
  if (ensemble_size(1) == 1) then
    return
  endif

  ! ensemble_size(1) is the number of ensemble members
  ! ensemble_size(2) is the number of pes per ensemble
  sbuf_size = size(sbuf)
  allocate(rbuf(ensemble_size(1) * sbuf_size))

  allocate(ensemble_pelist(ensemble_size(1), ensemble_size(2)))
  call get_ensemble_pelist(ensemble_pelist)

  ! Gather to root of every ensemble member.
  call mpp_gather(sbuf, rbuf, ensemble_pelist(:, 1))

  ! Check that sbuf is the same on all pes.
  if (transform_on_this_pe) then
    do e=0,ensemble_size(1)-1
      do i=1,sbuf_size
        a = sbuf(i)
        b = rbuf(e*sbuf_size + i)
        if ((a /= b) .and. (abs(a - b) /= 0.0)) then
          print*, a - b
          ret = i
          exit
        endif
      enddo
    enddo
  endif

  deallocate(rbuf)

  call mpp_max(ret, ensemble_pelist(:, 1))

end subroutine ensemble_compare_1d


subroutine rot90_2d(arrayIn, arrayOut, nrot90)
  real, dimension(:,:), intent(in) :: arrayIn !< Array to be rotated
  real, dimension(:,:), intent(out) :: arrayOut !< Rotated array
  integer, intent(in) :: nrot90 !< Number of 90 degree rotations to perform

  if (.not. nrot90 < 4) then
    call MOM_error(FATAL, 'rot90_2d: nrot should be < 4')
  endif

  if (modulo(nrot90, 2) == 2) then
    if (size(arrayIn, 1) /= size(arrayOut, 1) .or. size(arrayIn, 2) /= size(arrayOut, 2)) then
      call MOM_error(FATAL, 'rot90_2d: 180 deg rotation bad array shapes.')
    endif
  else
    if (size(arrayIn, 1) /= size(arrayOut, 2) .or. size(arrayIn, 2) /= size(arrayOut, 1)) then
      call MOM_error(FATAL, 'rot90_2d: 90 deg rotation bad array shapes.')
    endif
  endif

  if (nrot90 == 1) then
    ! transpose, reverse rows
    arrayOut(:, :) = transpose(arrayIn(:, :))
    arrayOut(:, :) = arrayOut(:, ubound(arrayOut, 2):lbound(arrayOut, 2):-1)
  elseif (nrot90 == 2) then
    ! reverse both rows and cols
    arrayOut(:, :) = arrayIn(ubound(arrayIn, 1):lbound(arrayIn, 1):-1, &
                             ubound(arrayIn, 2):lbound(arrayIn, 2):-1)
  elseif (nrot90 == 3) then
    ! transpose, reverse cols
    arrayOut(:,:) = transpose(arrayIn(:, :))
    arrayOut(:, :) = arrayOut(ubound(arrayOut, 1):lbound(arrayOut, 1):-1, :)
  endif

end subroutine rot90_2d

subroutine rot90_3d(arrayIn, arrayOut, nrot90)
  real, dimension(:,:,:), intent(in) :: arrayIn !< The array to be rotated
  real, dimension(:,:,:), intent(inout) :: arrayOut !< Rotated array
  integer, intent(in) :: nrot90 !< Number of 90 degree rotations to perform

  integer :: k

  if (.not. nrot90 < 4) then
    call MOM_error(FATAL, 'rot90_3d: nrot should be < 4')
  endif

  do k=lbound(arrayIn, 3), ubound(arrayIn, 3)
     call rot90_2d(arrayIn(:, :, k), arrayOut(:, :, k), nrot90)
  enddo

end subroutine rot90_3d

subroutine write_to_netcdf_3d(array, file_name)
  use netcdf
  implicit none

  real, intent(in), dimension(:,:,:) :: array
  character(len=*), intent(in) :: file_name

  integer :: file_id, xdim_id, ydim_id, zdim_id
  integer :: array_id
  integer, dimension(3) :: arrdims
  character(len=*), parameter :: arrunit = 'nondim'

  integer :: i, j, k
  integer :: ierr

  i = size(array,1)
  j = size(array,2)
  k = size(array,3)

  ! create the file
  ierr = nf90_create(path=trim(file_name), cmode=NF90_CLOBBER, ncid=file_id)

  ! define the dimensions
  ierr = nf90_def_dim(file_id, 'X', i, xdim_id)
  ierr = nf90_def_dim(file_id, 'Y', j, ydim_id)
  ierr = nf90_def_dim(file_id, 'Z', k, zdim_id)

  ! now that the dimensions are defined, we can define variables on them,...
  arrdims = (/ xdim_id, ydim_id, zdim_id /)
  ierr = nf90_def_var(file_id, 'Array',  NF90_DOUBLE, arrdims, array_id)

  ! ...and assign units to them as an attribute
  ierr = nf90_put_att(file_id, array_id, "units", arrunit)

  ! done defining
  ierr = nf90_enddef(file_id)

  ! Write out the values
  ierr = nf90_put_var(file_id, array_id, array)

  ! close; done
  ierr = nf90_close(file_id)
end subroutine write_to_netcdf_3d


subroutine write_to_netcdf_2d(array, file_name)
  use netcdf
  implicit none

  real, intent(in), dimension(:,:) :: array
  character(len=*), intent(in) :: file_name

  integer :: file_id, xdim_id, ydim_id
  integer :: array_id
  integer, dimension(2) :: arrdims
  character(len=*), parameter :: arrunit = 'nondim'

  integer :: i, j
  integer :: ierr

  i = size(array,1)
  j = size(array,2)

  ! create the file
  ierr = nf90_create(path=trim(file_name), cmode=NF90_CLOBBER, ncid=file_id)

  ! define the dimensions
  ierr = nf90_def_dim(file_id, 'X', i, xdim_id)
  ierr = nf90_def_dim(file_id, 'Y', j, ydim_id)

  ! now that the dimensions are defined, we can define variables on them,...
  arrdims = (/ xdim_id, ydim_id /)
  ierr = nf90_def_var(file_id, 'Array',  NF90_DOUBLE, arrdims, array_id)

  ! ...and assign units to them as an attribute
  ierr = nf90_put_att(file_id, array_id, "units", arrunit)

  ! done defining
  ierr = nf90_enddef(file_id)

  ! Write out the values
  ierr = nf90_put_var(file_id, array_id, array)

  ! close; done
  ierr = nf90_close(file_id)
end subroutine write_to_netcdf_2d

end module MOM_transform_test

