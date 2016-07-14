module MOM_spatial_means

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

use MOM_coms, only : EFP_type, operator(+), operator(-), assignment(=)
use MOM_coms, only : EFP_to_real, real_to_EFP, EFP_list_sum_across_PEs
use MOM_coms, only : reproducing_sum
use MOM_coms, only : query_EFP_overflow_error, reset_EFP_overflow_error
use MOM_error_handler, only : MOM_error, NOTE, WARNING, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public :: global_i_mean, global_j_mean
public :: global_area_mean, global_layer_mean
public :: global_area_integral
public :: global_volume_mean
public :: adjust_area_mean_to_zero

contains

function global_area_mean(var,G)
  type(ocean_grid_type),             intent(in)  :: G
  real, dimension(SZI_(G), SZJ_(G)), intent(in)  :: var
  real, dimension(SZI_(G), SZJ_(G))              :: tmpForSumming
  real :: global_area_mean

  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  tmpForSumming(:,:) = 0.
  do j=js,je ; do i=is, ie
    tmpForSumming(i,j) = ( var(i,j) * (G%areaT(i,j) * G%mask2dT(i,j)) )
  enddo ; enddo
  global_area_mean = reproducing_sum( tmpForSumming ) * G%IareaT_global

end function global_area_mean

function global_area_integral(var,G)
  type(ocean_grid_type),             intent(in)  :: G
  real, dimension(SZI_(G), SZJ_(G)), intent(in)  :: var
  real, dimension(SZI_(G), SZJ_(G))              :: tmpForSumming
  real :: global_area_integral

  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  tmpForSumming(:,:) = 0.
  do j=js,je ; do i=is, ie
    tmpForSumming(i,j) = ( var(i,j) * (G%areaT(i,j) * G%mask2dT(i,j)) )
  enddo ; enddo
  global_area_integral = reproducing_sum( tmpForSumming )

end function global_area_integral

function global_layer_mean(var, h, G, GV)
  type(ocean_grid_type),                     intent(in)  :: G
  type(verticalGrid_type),                   intent(in)  :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: var
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h
  real, dimension(SZK_(GV))                   :: global_layer_mean

  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)) :: tmpForSumming, weight
  real, dimension(SZK_(GV)) :: scalarij, weightij
  real, dimension(SZK_(GV)) :: global_temp_scalar, global_weight_scalar
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  tmpForSumming(:,:,:) = 0. ; weight(:,:,:) = 0.

  do k=1,nz ; do j=js,je ; do i=is,ie
    weight(i,j,k)  =  h(i,j,k) * (G%areaT(i,j) * G%mask2dT(i,j))
    tmpForSumming(i,j,k) =  var(i,j,k) * weight(i,j,k)
  enddo ; enddo ; enddo

  global_temp_scalar   = reproducing_sum(tmpForSumming,sums=scalarij)
  global_weight_scalar = reproducing_sum(weight,sums=weightij)

  do k=1, nz
    global_layer_mean(k) = scalarij(k) / weightij(k)
  enddo

end function global_layer_mean

function global_volume_mean(var, h, G, GV)
  type(ocean_grid_type),                     intent(in)  :: G
  type(verticalGrid_type),                   intent(in)  :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: var
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h
  real :: global_volume_mean

  real :: weight_here
  real, dimension(SZI_(G), SZJ_(G)) :: tmpForSumming, sum_weight
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  tmpForSumming(:,:) = 0. ; sum_weight(:,:) = 0.

  do k=1,nz ; do j=js,je ; do i=is,ie
    weight_here  =  h(i,j,k) * (G%areaT(i,j) * G%mask2dT(i,j))
    tmpForSumming(i,j) = tmpForSumming(i,j) + var(i,j,k) * weight_here
    sum_weight(i,j) = sum_weight(i,j) + weight_here
  enddo ; enddo ; enddo
  global_volume_mean = (reproducing_sum(tmpForSumming)) / &
                       (reproducing_sum(sum_weight))

end function global_volume_mean


subroutine global_i_mean(array, i_mean, G, mask)
  type(ocean_grid_type),            intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: array
  real, dimension(SZJ_(G)),         intent(out)   :: i_mean
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: mask

!    This subroutine determines the global mean of a field along rows of
!  constant i, returning it in a 1-d array using the local indexing.

! Arguments: array - The 2-d array whose i-mean is to be taken.
!  (out)     i_mean - Global mean of array along its i-axis.
!  (in)      G - The ocean's grid structure.
!  (in)      mask - An array used for weighting the i-mean.

  type(EFP_type), allocatable, dimension(:) :: asum, mask_sum
  real :: mask_sum_r
  integer :: is, ie, js, je, idg_off, jdg_off
  integer :: i, j

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  idg_off = G%idg_offset ; jdg_off = G%jdg_offset

  call reset_EFP_overflow_error()

  allocate(asum(G%jsg:G%jeg))
  if (present(mask)) then
    allocate(mask_sum(G%jsg:G%jeg))

    do j=G%jsg,G%jeg
      asum(j) = real_to_EFP(0.0) ; mask_sum(j) = real_to_EFP(0.0)
    enddo

    do i=is,ie ; do j=js,je
      asum(j+jdg_off) = asum(j+jdg_off) + real_to_EFP(array(i,j)*mask(i,j))
      mask_sum(j+jdg_off) = mask_sum(j+jdg_off) + real_to_EFP(mask(i,j))
    enddo ; enddo

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_i_mean overflow error occurred before sums across PEs.")

    call EFP_list_sum_across_PEs(asum(G%jsg:G%jeg), G%jeg-G%jsg+1)
    call EFP_list_sum_across_PEs(mask_sum(G%jsg:G%jeg), G%jeg-G%jsg+1)

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_i_mean overflow error occurred during sums across PEs.")

    do j=js,je
      mask_sum_r = EFP_to_real(mask_sum(j+jdg_off))
      if (mask_sum_r == 0.0 ) then ; i_mean(j) = 0.0 ; else
        i_mean(j) = EFP_to_real(asum(j+jdg_off)) / mask_sum_r
      endif
    enddo

    deallocate(mask_sum)
  else
    do j=G%jsg,G%jeg ; asum(j) = real_to_EFP(0.0) ; enddo

    do i=is,ie ; do j=js,je
      asum(j+jdg_off) = asum(j+jdg_off) + real_to_EFP(array(i,j))
    enddo ; enddo

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_i_mean overflow error occurred before sum across PEs.")

    call EFP_list_sum_across_PEs(asum(G%jsg:G%jeg), G%jeg-G%jsg+1)

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_i_mean overflow error occurred during sum across PEs.")

    do j=js,je
      i_mean(j) = EFP_to_real(asum(j+jdg_off)) / real(G%ieg-G%isg+1)
    enddo
  endif

  deallocate(asum)

end subroutine global_i_mean

subroutine global_j_mean(array, j_mean, G, mask)
  type(ocean_grid_type),            intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: array
  real, dimension(SZI_(G)),         intent(out)   :: j_mean
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: mask

!    This subroutine determines the global mean of a field along rows of
!  constant j, returning it in a 1-d array using the local indexing.

! Arguments: array - The 2-d array whose j-mean is to be taken.
!  (out)     j_mean - Global mean of array along its j-axis.
!  (in)      G - The ocean's grid structure.
!  (in)      mask - An array used for weighting the j-mean.

  type(EFP_type), allocatable, dimension(:) :: asum, mask_sum
  real :: mask_sum_r
  integer :: is, ie, js, je, idg_off, jdg_off
  integer :: i, j

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  idg_off = G%idg_offset ; jdg_off = G%jdg_offset

  call reset_EFP_overflow_error()

  allocate(asum(G%isg:G%ieg))
  if (present(mask)) then
    allocate (mask_sum(G%isg:G%ieg))

    do i=G%isg,G%ieg
      asum(i) = real_to_EFP(0.0) ; mask_sum(i) = real_to_EFP(0.0)
    enddo

    do i=is,ie ; do j=js,je
      asum(i+idg_off) = asum(i+idg_off) + real_to_EFP(array(i,j)*mask(i,j))
      mask_sum(i+idg_off) = mask_sum(i+idg_off) + real_to_EFP(mask(i,j))
    enddo ; enddo

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_j_mean overflow error occurred before sums across PEs.")

    call EFP_list_sum_across_PEs(asum(G%isg:G%ieg), G%ieg-G%isg+1)
    call EFP_list_sum_across_PEs(mask_sum(G%isg:G%ieg), G%ieg-G%isg+1)

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_j_mean overflow error occurred during sums across PEs.")

    do i=is,ie
      mask_sum_r = EFP_to_real(mask_sum(i+idg_off))
      if (mask_sum_r == 0.0 ) then ; j_mean(i) = 0.0 ; else
        j_mean(i) = EFP_to_real(asum(i+idg_off)) / mask_sum_r
      endif
    enddo

    deallocate(mask_sum)
  else
    do i=G%isg,G%ieg ; asum(i) = real_to_EFP(0.0) ; enddo

    do i=is,ie ; do j=js,je
      asum(i+idg_off) = asum(i+idg_off) + real_to_EFP(array(i,j))
    enddo ; enddo

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_j_mean overflow error occurred before sum across PEs.")

    call EFP_list_sum_across_PEs(asum(G%isg:G%ieg), G%ieg-G%isg+1)

    if (query_EFP_overflow_error()) call MOM_error(FATAL, &
      "global_j_mean overflow error occurred during sum across PEs.")

    do i=is,ie
      j_mean(i) = EFP_to_real(asum(i+idg_off)) / real(G%jeg-G%jsg+1)
    enddo
  endif

  deallocate(asum)

end subroutine global_j_mean

!> Adjust 2d array such that area mean is zero without moving the zero contour
subroutine adjust_area_mean_to_zero(array, G, scaling)
  type(ocean_grid_type),            intent(inout) :: G       !< Grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: array   !< 2D array to be adjusted
  real, optional,                   intent(out)   :: scaling !< The scaling factor used
  ! Local variables
  real, dimension(SZI_(G), SZJ_(G)) :: posVals, negVals, areaXposVals, areaXnegVals
  integer :: i,j
  real :: areaIntPosVals, areaIntNegVals, posScale, negScale

  areaXposVals(:,:) = 0.
  areaXnegVals(:,:) = 0.

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    posVals(i,j) = max(0., array(i,j))
    areaXposVals(i,j) = G%areaT(i,j) * posVals(i,j)
    negVals(i,j) = min(0., array(i,j))
    areaXnegVals(i,j) = G%areaT(i,j) * negVals(i,j)
  enddo ; enddo

  areaIntPosVals = reproducing_sum( areaXposVals )
  areaIntNegVals = reproducing_sum( areaXnegVals )

  posScale = 0.0 ; negScale = 0.0
  if ((areaIntPosVals>0.).and.(areaIntNegVals<0.)) then ! Only adjust if possible
    if (areaIntPosVals>-areaIntNegVals) then ! Scale down positive values
      posScale = - areaIntNegVals / areaIntPosVals
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        array(i,j) = (posScale * posVals(i,j)) + negVals(i,j)
      enddo ; enddo
    elseif (areaIntPosVals<-areaIntNegVals) then ! Scale down negative values
      negScale = - areaIntPosVals / areaIntNegVals
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        array(i,j) = posVals(i,j) + (negScale * negVals(i,j))
      enddo ; enddo
    endif
  endif
  if (present(scaling)) scaling = posScale - negScale

end subroutine adjust_area_mean_to_zero

end module MOM_spatial_means
