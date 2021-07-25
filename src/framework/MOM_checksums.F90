!> Routines to calculate checksums of various array and vector types
module MOM_checksums

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_array_transform, only : rotate_array, rotate_array_pair, rotate_vector
use MOM_array_transform, only : allocate_rotated_array
use MOM_coms,            only : PE_here, root_PE, num_PEs, sum_across_PEs
use MOM_coms,            only : min_across_PEs, max_across_PEs
use MOM_coms,            only : reproducing_sum, field_chksum
use MOM_error_handler,   only : MOM_error, FATAL, is_root_pe
use MOM_file_parser,     only : log_version, param_file_type
use MOM_hor_index,       only : hor_index_type, rotate_hor_index

use iso_fortran_env,     only : error_unit, int32, int64

implicit none ; private

public :: chksum0, zchksum, rotated_field_chksum
public :: hchksum, Bchksum, uchksum, vchksum, qchksum, is_NaN, chksum
public :: hchksum_pair, uvchksum, Bchksum_pair
public :: MOM_checksums_init

!> Checksums a pair of arrays (2d or 3d) staggered at tracer points
interface hchksum_pair
  module procedure chksum_pair_h_2d, chksum_pair_h_3d
end interface

!> Checksums a pair velocity arrays (2d or 3d) staggered at C-grid locations
interface uvchksum
  module procedure chksum_uv_2d, chksum_uv_3d
end interface

!> Checksums an array (2d or 3d) staggered at C-grid u points.
interface uchksum
  module procedure chksum_u_2d, chksum_u_3d
end interface

!> Checksums an array (2d or 3d) staggered at C-grid v points.
interface vchksum
  module procedure chksum_v_2d, chksum_v_3d
end interface

!> Checksums a pair of arrays (2d or 3d) staggered at corner points
interface Bchksum_pair
  module procedure chksum_pair_B_2d, chksum_pair_B_3d
end interface

!> Checksums an array (2d or 3d) staggered at tracer points.
interface hchksum
  module procedure chksum_h_2d, chksum_h_3d
end interface

!> Checksums an array (2d or 3d) staggered at corner points.
interface Bchksum
  module procedure chksum_B_2d, chksum_B_3d
end interface

!> This is an older interface that has been renamed Bchksum
interface qchksum
  module procedure chksum_B_2d, chksum_B_3d
end interface

!> This is an older interface for 1-, 2-, or 3-D checksums
interface chksum
  module procedure chksum1d, chksum2d, chksum3d
end interface

!> Write a message with either checksums or numerical statistics of arrays
interface chk_sum_msg
  module procedure chk_sum_msg1, chk_sum_msg2, chk_sum_msg3, chk_sum_msg5
end interface

!> Returns .true. if any element of x is a NaN, and .false. otherwise.
interface is_NaN
  module procedure is_NaN_0d, is_NaN_1d, is_NaN_2d, is_NaN_3d
end interface

!> Rotate and compute the checksum of a field
interface rotated_field_chksum
  module procedure rotated_field_chksum_real_0d
  module procedure rotated_field_chksum_real_1d
  module procedure rotated_field_chksum_real_2d
  module procedure rotated_field_chksum_real_3d
  module procedure rotated_field_chksum_real_4d
end interface rotated_field_chksum

integer, parameter :: bc_modulus = 1000000000 !< Modulus of checksum bitcount
integer, parameter :: default_shift=0 !< The default array shift
logical :: calculateStatistics=.true. !< If true, report min, max and mean.
logical :: writeChksums=.true. !< If true, report the bitcount checksum
logical :: checkForNaNs=.true. !< If true, checks array for NaNs and cause
                               !! FATAL error is any are found

contains

!> Checksum a scalar field (consistent with array checksums)
subroutine chksum0(scalar, mesg, scale, logunit)
  real, intent(in) :: scalar                !< The array to be checksummed
  character(len=*), intent(in) :: mesg     !< An identifying message
  real, optional, intent(in) :: scale      !< A scaling factor for this array.
  integer, optional, intent(in) :: logunit !< IO unit for checksum logging

  real :: scaling   !< Explicit rescaling factor
  integer :: iounit !< Log IO unit
  real :: rs        !< Rescaled scalar
  integer :: bc     !< Scalar bitcount

  if (checkForNaNs .and. is_NaN(scalar)) &
    call chksum_error(FATAL, 'NaN detected: '//trim(mesg))

  scaling = 1.0 ; if (present(scale)) scaling = scale
  iounit = error_unit; if(present(logunit)) iounit = logunit

  if (calculateStatistics) then
    rs = scaling * scalar
    if (is_root_pe()) &
      call chk_sum_msg(" scalar:", rs, rs, rs, mesg, iounit)
  endif

  if (.not. writeChksums) return

  bc = mod(bitcount(abs(scaling * scalar)), bc_modulus)
  if (is_root_pe()) &
    call chk_sum_msg(" scalar:", bc, mesg, iounit)

end subroutine chksum0


!> Checksum a 1d array (typically a column).
subroutine zchksum(array, mesg, scale, logunit)
  real, dimension(:), intent(in) :: array  !< The array to be checksummed
  character(len=*), intent(in) :: mesg     !< An identifying message
  real, optional, intent(in) :: scale      !< A scaling factor for this array.
  integer, optional, intent(in) :: logunit !< IO unit for checksum logging

  real, allocatable, dimension(:) :: rescaled_array
  real :: scaling
  integer :: iounit !< Log IO unit
  integer :: k
  real :: aMean, aMin, aMax
  integer :: bc0

  if (checkForNaNs) then
    if (is_NaN(array(:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
  endif

  scaling = 1.0 ; if (present(scale)) scaling = scale
  iounit = error_unit; if(present(logunit)) iounit = logunit

  if (calculateStatistics) then
    if (present(scale)) then
      allocate(rescaled_array(LBOUND(array,1):UBOUND(array,1)))
      rescaled_array(:) = 0.0
      do k=1, size(array, 1)
        rescaled_array(k) = scale * array(k)
      enddo

      call subStats(rescaled_array, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(array, aMean, aMin, aMax)
    endif

    if (is_root_pe()) &
      call chk_sum_msg(" column:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not. writeChksums) return

  bc0 = subchk(array, scaling)
  if (is_root_pe()) call chk_sum_msg(" column:", bc0, mesg, iounit)

  contains

  integer function subchk(array, scale)
    real, dimension(:), intent(in) :: array !< The array to be checksummed
    real, intent(in) :: scale !< A scaling factor for this array.
    integer :: k, bc
    subchk = 0
    do k=LBOUND(array, 1), UBOUND(array, 1)
      bc = bitcount(abs(scale * array(k)))
      subchk = subchk + bc
    enddo
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(array, aMean, aMin, aMax)
    real, dimension(:), intent(in) :: array !< The array to be checksummed
    real, intent(out) :: aMean !< Array mean
    real, intent(out) :: aMin !< Array minimum
    real, intent(out) :: aMax !< Array maximum

    integer :: k, n

    aMin = array(1)
    aMax = array(1)
    n = 0
    do k=LBOUND(array,1), UBOUND(array,1)
      aMin = min(aMin, array(k))
      aMax = max(aMax, array(k))
      n = n + 1
    enddo
    aMean = sum(array(:)) / real(n)
  end subroutine subStats
end subroutine zchksum

!> Checksums on a pair of 2d arrays staggered at tracer points.
subroutine chksum_pair_h_2d(mesg, arrayA, arrayB, HI, haloshift, omit_corners, &
                            scale, logunit, scalar_pair)
  character(len=*),                 intent(in) :: mesg !< Identifying messages
  type(hor_index_type),   target,   intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:), target, intent(in) :: arrayA !< The first array to be checksummed
  real, dimension(HI%isd:,HI%jsd:), target, intent(in) :: arrayB !< The second array to be checksummed
  integer,                optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                   optional, intent(in) :: scale     !< A scaling factor for this array.
  integer,                optional, intent(in) :: logunit !< IO unit for checksum logging
  logical,                optional, intent(in) :: scalar_pair !< If true, then the arrays describe
                                                              !! a scalar, rather than vector
  logical :: vector_pair
  integer :: turns
  type(hor_index_type), pointer :: HI_in
  real, dimension(:,:), pointer :: arrayA_in, arrayB_in

  vector_pair = .true.
  if (present(scalar_pair)) vector_pair = .not. scalar_pair

  turns = HI%turns
  if (modulo(turns, 4) /= 0) then
    ! Rotate field back to the input grid
    allocate(HI_in)
    call rotate_hor_index(HI, -turns, HI_in)
    allocate(arrayA_in(HI_in%isd:HI_in%ied, HI_in%jsd:HI_in%jed))
    allocate(arrayB_in(HI_in%isd:HI_in%ied, HI_in%jsd:HI_in%jed))

    if (vector_pair) then
      call rotate_vector(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    else
      call rotate_array_pair(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    endif
  else
    HI_in => HI
    arrayA_in => arrayA
    arrayB_in => arrayB
  endif

  if (present(haloshift)) then
    call chksum_h_2d(arrayA_in, 'x '//mesg, HI_in, haloshift, omit_corners, &
                     scale=scale, logunit=logunit)
    call chksum_h_2d(arrayB_in, 'y '//mesg, HI_in, haloshift, omit_corners, &
                     scale=scale, logunit=logunit)
  else
    call chksum_h_2d(arrayA_in, 'x '//mesg, HI_in, scale=scale, logunit=logunit)
    call chksum_h_2d(arrayB_in, 'y '//mesg, HI_in, scale=scale, logunit=logunit)
  endif
end subroutine chksum_pair_h_2d

!> Checksums on a pair of 3d arrays staggered at tracer points.
subroutine chksum_pair_h_3d(mesg, arrayA, arrayB, HI, haloshift, omit_corners, &
                            scale, logunit, scalar_pair)
  character(len=*),                    intent(in) :: mesg !< Identifying messages
  type(hor_index_type),      target,   intent(in) :: HI   !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:, :), target, intent(in) :: arrayA !< The first array to be checksummed
  real, dimension(HI%isd:,HI%jsd:, :), target, intent(in) :: arrayB !< The second array to be checksummed
  integer,                   optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                   optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                      optional, intent(in) :: scale     !< A scaling factor for this array.
  integer,                   optional, intent(in) :: logunit !< IO unit for checksum logging

  logical,                optional, intent(in) :: scalar_pair !< If true, then the arrays describe
                                                              !! a scalar, rather than vector
  logical :: vector_pair
  integer :: turns
  type(hor_index_type), pointer :: HI_in
  real, dimension(:,:,:), pointer :: arrayA_in, arrayB_in

  vector_pair = .true.
  if (present(scalar_pair)) vector_pair = .not. scalar_pair

  turns = HI%turns
  if (modulo(turns, 4) /= 0) then
    ! Rotate field back to the input grid
    allocate(HI_in)
    call rotate_hor_index(HI, -turns, HI_in)
    allocate(arrayA_in(HI_in%isd:HI_in%ied, HI_in%jsd:HI_in%jed, size(arrayA, 3)))
    allocate(arrayB_in(HI_in%isd:HI_in%ied, HI_in%jsd:HI_in%jed, size(arrayB, 3)))

    if (vector_pair) then
      call rotate_vector(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    else
      call rotate_array_pair(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    endif
  else
    HI_in => HI
    arrayA_in => arrayA
    arrayB_in => arrayB
  endif

  if (present(haloshift)) then
    call chksum_h_3d(arrayA_in, 'x '//mesg, HI_in, haloshift, omit_corners, &
                     scale=scale, logunit=logunit)
    call chksum_h_3d(arrayB_in, 'y '//mesg, HI_in, haloshift, omit_corners, &
                     scale=scale, logunit=logunit)
  else
    call chksum_h_3d(arrayA_in, 'x '//mesg, HI_in, scale=scale, logunit=logunit)
    call chksum_h_3d(arrayB_in, 'y '//mesg, HI_in, scale=scale, logunit=logunit)
  endif

  ! NOTE: automatic deallocation of array[AB]_in
end subroutine chksum_pair_h_3d

!> Checksums a 2d array staggered at tracer points.
subroutine chksum_h_2d(array_m, mesg, HI_m, haloshift, omit_corners, scale, logunit)
  type(hor_index_type), target, intent(in) :: HI_m    !< Horizontal index bounds of the model grid
  real, dimension(HI_m%isd:,HI_m%jsd:), target, intent(in) :: array_m !< Field array on the model grid
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,               optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                  optional, intent(in) :: scale     !< A scaling factor for this array.
  integer, optional, intent(in) :: logunit !< IO unit for checksum logging

  real, pointer :: array(:,:)           ! Field array on the input grid
  real, allocatable, dimension(:,:) :: rescaled_array
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling
  integer :: iounit !< Log IO unit
  integer :: i, j
  real :: aMean, aMin, aMax
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners
  integer :: turns                      ! Quarter turns from input to model grid

  ! Rotate array to the input grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    allocate(array(HI%isd:HI%ied, HI%jsd:HI%jed))
    call rotate_array(array_m, -turns, array)
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%jsc:HI%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0 ; if (present(scale)) scaling = scale
  iounit = error_unit; if(present(logunit)) iounit = logunit

  if (calculateStatistics) then
    if (present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2)) )
      rescaled_array(:,:) = 0.0
      do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
        rescaled_array(i,j) = scale*array(i,j)
      enddo ; enddo
      call subStats(HI, rescaled_array, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, aMean, aMin, aMax)
    endif

    if (is_root_pe()) &
      call chk_sum_msg("h-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_h_2d: haloshift =',hshift
    write(0,*) 'chksum_h_2d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_h_2d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_h_2d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  if (hshift==0) then
    if (is_root_pe()) call chk_sum_msg("h-point:", bc0, mesg, iounit)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (do_corners) then
    bcSW = subchk(array, HI, -hshift, -hshift, scaling)
    bcSE = subchk(array, HI, hshift, -hshift, scaling)
    bcNW = subchk(array, HI, -hshift, hshift, scaling)
    bcNE = subchk(array, HI, hshift, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg("h-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
  else
    bcS = subchk(array, HI, 0, -hshift, scaling)
    bcE = subchk(array, HI, hshift, 0, scaling)
    bcW = subchk(array, HI, -hshift, 0, scaling)
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg_NSEW("h-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
  endif

  contains
  integer function subchk(array, HI, di, dj, scale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%jsd:), intent(in) :: array !< The array to be checksummed
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: scale !< A scaling factor for this array.
    integer :: i, j, bc
    subchk = 0
    do j=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(scale*array(i,j)))
      subchk = subchk + bc
    enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%jsd:), intent(in) :: array !< The array to be checksummed
    real, intent(out) :: aMean !< Array mean
    real, intent(out) :: aMin !< Array minimum
    real, intent(out) :: aMax !< Array maximum

    integer :: i, j, n

    aMin = array(HI%isc,HI%jsc)
    aMax = array(HI%isc,HI%jsc)
    n = 0
    do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,j))
      aMax = max(aMax, array(i,j))
      n = n + 1
    enddo ; enddo
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec))
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end subroutine chksum_h_2d

!> Checksums on a pair of 2d arrays staggered at q-points.
subroutine chksum_pair_B_2d(mesg, arrayA, arrayB, HI, haloshift, symmetric, &
                            omit_corners, scale, logunit, scalar_pair)
  character(len=*),                 intent(in) :: mesg   !< Identifying messages
  type(hor_index_type),   target,   intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:), target, intent(in) :: arrayA !< The first array to be checksummed
  real, dimension(HI%isd:,HI%jsd:), target, intent(in) :: arrayB !< The second array to be checksummed
  logical,                optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                            !! symmetric computational domain.
  integer,                optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                   optional, intent(in) :: scale     !< A scaling factor for this array.
  integer,                optional, intent(in) :: logunit !< IO unit for checksum logging
  logical,                optional, intent(in) :: scalar_pair !< If true, then the arrays describe
                                                              !! a scalar, rather than vector

  logical :: sym
  logical :: vector_pair
  integer :: turns
  type(hor_index_type), pointer :: HI_in
  real, dimension(:,:), pointer :: arrayA_in, arrayB_in

  vector_pair = .true.
  if (present(scalar_pair)) vector_pair = .not. scalar_pair

  turns = HI%turns
  if (modulo(turns, 4) /= 0) then
    ! Rotate field back to the input grid
    allocate(HI_in)
    call rotate_hor_index(HI, -turns, HI_in)
    allocate(arrayA_in(HI_in%IsdB:HI_in%IedB, HI_in%JsdB:HI_in%JedB))
    allocate(arrayB_in(HI_in%IsdB:HI_in%IedB, HI_in%JsdB:HI_in%JedB))

    if (vector_pair) then
      call rotate_vector(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    else
      call rotate_array_pair(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    endif
  else
    HI_in => HI
    arrayA_in => arrayA
    arrayB_in => arrayB
  endif

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if (present(haloshift)) then
    call chksum_B_2d(arrayA_in, 'x '//mesg, HI_in, haloshift, symmetric=sym, &
                     omit_corners=omit_corners, scale=scale, logunit=logunit)
    call chksum_B_2d(arrayB_in, 'y '//mesg, HI_in, haloshift, symmetric=sym, &
                     omit_corners=omit_corners, scale=scale, logunit=logunit)
  else
    call chksum_B_2d(arrayA_in, 'x '//mesg, HI_in, symmetric=sym, scale=scale, &
                     logunit=logunit)
    call chksum_B_2d(arrayB_in, 'y '//mesg, HI_in, symmetric=sym, scale=scale, &
                     logunit=logunit)
  endif

end subroutine chksum_pair_B_2d

!> Checksums on a pair of 3d arrays staggered at q-points.
subroutine chksum_pair_B_3d(mesg, arrayA, arrayB, HI, haloshift, symmetric, &
                            omit_corners, scale, logunit, scalar_pair)
  character(len=*),                    intent(in) :: mesg !< Identifying messages
  type(hor_index_type),      target,   intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%JsdB:, :), target, intent(in) :: arrayA !< The first array to be checksummed
  real, dimension(HI%IsdB:,HI%JsdB:, :), target, intent(in) :: arrayB !< The second array to be checksummed
  integer,                   optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                   optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                               !! symmetric computational domain.
  logical,                   optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                      optional, intent(in) :: scale     !< A scaling factor for this array.
  integer,                   optional, intent(in) :: logunit !< IO unit for checksum logging
  logical,                   optional, intent(in) :: scalar_pair !< If true, then the arrays describe
                                                              !! a scalar, rather than vector

  logical :: sym
  logical :: vector_pair
  integer :: turns
  type(hor_index_type), pointer :: HI_in
  real, dimension(:,:,:), pointer :: arrayA_in, arrayB_in

  vector_pair = .true.
  if (present(scalar_pair)) vector_pair = .not. scalar_pair

  turns = HI%turns
  if (modulo(turns, 4) /= 0) then
    ! Rotate field back to the input grid
    allocate(HI_in)
    call rotate_hor_index(HI, -turns, HI_in)
    allocate(arrayA_in(HI_in%IsdB:HI_in%IedB, HI_in%JsdB:HI_in%JedB, size(arrayA, 3)))
    allocate(arrayB_in(HI_in%IsdB:HI_in%IedB, HI_in%JsdB:HI_in%JedB, size(arrayB, 3)))

    if (vector_pair) then
      call rotate_vector(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    else
      call rotate_array_pair(arrayA, arrayB, -turns, arrayA_in, arrayB_in)
    endif
  else
    HI_in => HI
    arrayA_in => arrayA
    arrayB_in => arrayB
  endif

  if (present(haloshift)) then
    call chksum_B_3d(arrayA_in, 'x '//mesg, HI_in, haloshift, symmetric, &
                     omit_corners, scale=scale, logunit=logunit)
    call chksum_B_3d(arrayB_in, 'y '//mesg, HI_in, haloshift, symmetric, &
                     omit_corners, scale=scale, logunit=logunit)
  else
    call chksum_B_3d(arrayA_in, 'x '//mesg, HI_in, symmetric=symmetric, scale=scale, &
                     logunit=logunit)
    call chksum_B_3d(arrayB_in, 'y '//mesg, HI_in, symmetric=symmetric, scale=scale, &
                     logunit=logunit)
  endif
end subroutine chksum_pair_B_3d

!> Checksums a 2d array staggered at corner points.
subroutine chksum_B_2d(array_m, mesg, HI_m, haloshift, symmetric, omit_corners, &
                       scale, logunit)
  type(hor_index_type), target, intent(in) :: HI_m     !< A horizontal index type
  real, dimension(HI_m%IsdB:,HI_m%JsdB:), &
                        target, intent(in) :: array_m !< The array to be checksummed
  character(len=*),     intent(in) :: mesg  !< An identifying message
  integer,    optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,    optional, intent(in) :: symmetric !< If true, do the checksums on the
                                                !! full symmetric computational domain.
  logical,    optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,       optional, intent(in) :: scale     !< A scaling factor for this array.
  integer, optional, intent(in) :: logunit !< IO unit for checksum logging

  real, pointer :: array(:,:)           ! Field array on the input grid
  real, allocatable, dimension(:,:) :: rescaled_array
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling
  integer :: iounit !< Log IO unit
  integer :: i, j, Is, Js
  real :: aMean, aMin, aMax
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats
  integer :: turns                      ! Quarter turns from input to model grid

  ! Rotate array to the input grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    allocate(array(HI%IsdB:HI%IedB, HI%JsdB:HI%JedB))
    call rotate_array(array_m, -turns, array)
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%JscB:HI%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0 ; if (present(scale)) scaling = scale
  iounit = error_unit; if(present(logunit)) iounit = logunit
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then
    if (present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2)) )
      rescaled_array(:,:) = 0.0
      Is = HI%isc ; if (sym_stats) Is = HI%isc-1
      Js = HI%jsc ; if (sym_stats) Js = HI%jsc-1
      do J=Js,HI%JecB ; do I=Is,HI%IecB
        rescaled_array(I,J) = scale*array(I,J)
      enddo ; enddo
      call subStats(HI, rescaled_array, sym_stats, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, sym_stats, aMean, aMin, aMax)
    endif
    if (is_root_pe()) &
      call chk_sum_msg("B-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%iscB-hshift<HI%isdB .or. HI%iecB+hshift>HI%iedB .or. &
       HI%jscB-hshift<HI%jsdB .or. HI%jecB+hshift>HI%jedB ) then
    write(0,*) 'chksum_B_2d: haloshift =',hshift
    write(0,*) 'chksum_B_2d: isd,isc,iec,ied=',HI%isdB,HI%iscB,HI%iecB,HI%iedB
    write(0,*) 'chksum_B_2d: jsd,jsc,jec,jed=',HI%jsdB,HI%jscB,HI%jecB,HI%jedB
    call chksum_error(FATAL,'Error in chksum_B_2d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if ((hshift==0) .and. .not.sym) then
    if (is_root_pe()) call chk_sum_msg("B-point:", bc0, mesg, iounit)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (do_corners) then
    if (sym) then
      bcSW = subchk(array, HI, -hshift-1, -hshift-1, scaling)
      bcSE = subchk(array, HI, hshift, -hshift-1, scaling)
      bcNW = subchk(array, HI, -hshift-1, hshift, scaling)
    else
      bcSW = subchk(array, HI, -hshift, -hshift, scaling)
      bcSE = subchk(array, HI, hshift, -hshift, scaling)
      bcNW = subchk(array, HI, -hshift, hshift, scaling)
    endif
    bcNE = subchk(array, HI, hshift, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg("B-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
  else
    bcS = subchk(array, HI, 0, -hshift, scaling)
    bcE = subchk(array, HI, hshift, 0, scaling)
    bcW = subchk(array, HI, -hshift, 0, scaling)
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg_NSEW("B-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
  endif

  contains

  integer function subchk(array, HI, di, dj, scale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%JsdB:), intent(in) :: array !< The array to be checksummed
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: scale !< A scaling factor for this array.
    integer :: i, j, bc
    subchk = 0
    ! This line deliberately uses the h-point computational domain.
    do J=HI%jsc+dj,HI%jec+dj ; do I=HI%isc+di,HI%iec+di
      bc = bitcount(abs(scale*array(I,J)))
      subchk = subchk + bc
    enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, sym_stats, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%JsdB:), intent(in) :: array !< The array to be checksummed
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.
    real, intent(out) :: aMean !< Array mean
    real, intent(out) :: aMin !< Array minimum
    real, intent(out) :: aMax !< Array maximum

    integer :: i, j, n, IsB, JsB

    IsB = HI%isc ; if (sym_stats) IsB = HI%isc-1
    JsB = HI%jsc ; if (sym_stats) JsB = HI%jsc-1

    aMin = array(HI%isc,HI%jsc) ; aMax = aMin
    do J=JsB,HI%JecB ; do I=IsB,HI%IecB
      aMin = min(aMin, array(I,J))
      aMax = max(aMax, array(I,J))
    enddo ; enddo
    ! This line deliberately uses the h-point computational domain.
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec))
    n = (1 + HI%jec - HI%jsc) * (1 + HI%iec - HI%isc)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end subroutine chksum_B_2d

!> Checksums a pair of 2d velocity arrays staggered at C-grid locations
subroutine chksum_uv_2d(mesg, arrayU, arrayV, HI, haloshift, symmetric, &
                        omit_corners, scale, logunit, scalar_pair)
  character(len=*),                  intent(in) :: mesg   !< Identifying messages
  type(hor_index_type),    target,   intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%jsd:), target, intent(in) :: arrayU !< The u-component array to be checksummed
  real, dimension(HI%isd:,HI%JsdB:), target, intent(in) :: arrayV !< The v-component array to be checksummed
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                 optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                             !! symmetric computational domain.
  logical,                 optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A scaling factor for these arrays.
  integer,                 optional, intent(in) :: logunit !< IO unit for checksum logging
  logical,                 optional, intent(in) :: scalar_pair !< If true, then the arrays describe a
                                                               !! a scalar, rather than vector
  logical :: vector_pair
  integer :: turns
  type(hor_index_type), pointer :: HI_in
  real, dimension(:,:), pointer :: arrayU_in, arrayV_in

  vector_pair = .true.
  if (present(scalar_pair)) vector_pair = .not. scalar_pair

  turns = HI%turns
  if (modulo(turns, 4) /= 0) then
    ! Rotate field back to the input grid
    allocate(HI_in)
    call rotate_hor_index(HI, -turns, HI_in)
    allocate(arrayU_in(HI_in%IsdB:HI_in%IedB, HI_in%jsd:HI_in%jed))
    allocate(arrayV_in(HI_in%isd:HI_in%ied, HI_in%JsdB:HI_in%JedB))

    if (vector_pair) then
      call rotate_vector(arrayU, arrayV, -turns, arrayU_in, arrayV_in)
    else
      call rotate_array_pair(arrayU, arrayV, -turns, arrayU_in, arrayV_in)
    endif
  else
    HI_in => HI
    arrayU_in => arrayU
    arrayV_in => arrayV
  endif

  if (present(haloshift)) then
    call chksum_u_2d(arrayU_in, 'u '//mesg, HI_in, haloshift, symmetric, &
                     omit_corners, scale=scale, logunit=logunit)
    call chksum_v_2d(arrayV_in, 'v '//mesg, HI_in, haloshift, symmetric, &
                     omit_corners, scale=scale, logunit=logunit)
  else
    call chksum_u_2d(arrayU_in, 'u '//mesg, HI_in, symmetric=symmetric, &
                     scale=scale, logunit=logunit)
    call chksum_v_2d(arrayV_in, 'v '//mesg, HI_in, symmetric=symmetric, &
                     scale=scale, logunit=logunit)
  endif
end subroutine chksum_uv_2d

!> Checksums a pair of 3d velocity arrays staggered at C-grid locations
subroutine chksum_uv_3d(mesg, arrayU, arrayV, HI, haloshift, symmetric, &
                        omit_corners, scale, logunit, scalar_pair)
  character(len=*),                    intent(in) :: mesg   !< Identifying messages
  type(hor_index_type),      target,   intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%jsd:,:), target, intent(in) :: arrayU !< The u-component array to be checksummed
  real, dimension(HI%isd:,HI%JsdB:,:), target, intent(in) :: arrayV !< The v-component array to be checksummed
  integer,                   optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                   optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                               !! symmetric computational domain.
  logical,                   optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                      optional, intent(in) :: scale     !< A scaling factor for these arrays.
  integer,                   optional, intent(in) :: logunit !< IO unit for checksum logging
  logical,                 optional, intent(in) :: scalar_pair !< If true, then the arrays describe a
                                                               !! a scalar, rather than vector
  logical :: vector_pair
  integer :: turns
  type(hor_index_type), pointer :: HI_in
  real, dimension(:,:,:), pointer :: arrayU_in, arrayV_in

  vector_pair = .true.
  if (present(scalar_pair)) vector_pair = .not. scalar_pair

  turns = HI%turns
  if (modulo(turns, 4) /= 0) then
    ! Rotate field back to the input grid
    allocate(HI_in)
    call rotate_hor_index(HI, -turns, HI_in)
    allocate(arrayU_in(HI_in%IsdB:HI_in%IedB, HI_in%jsd:HI_in%jed, size(arrayU, 3)))
    allocate(arrayV_in(HI_in%isd:HI_in%ied, HI_in%JsdB:HI_in%JedB, size(arrayV, 3)))

    if (vector_pair) then
      call rotate_vector(arrayU, arrayV, -turns, arrayU_in, arrayV_in)
    else
      call rotate_array_pair(arrayU, arrayV, -turns, arrayU_in, arrayV_in)
    endif
  else
    HI_in => HI
    arrayU_in => arrayU
    arrayV_in => arrayV
  endif

  if (present(haloshift)) then
    call chksum_u_3d(arrayU_in, 'u '//mesg, HI_in, haloshift, symmetric, &
                     omit_corners, scale=scale, logunit=logunit)
    call chksum_v_3d(arrayV_in, 'v '//mesg, HI_in, haloshift, symmetric, &
                     omit_corners, scale=scale, logunit=logunit)
  else
    call chksum_u_3d(arrayU_in, 'u '//mesg, HI_in, symmetric=symmetric, &
                     scale=scale, logunit=logunit)
    call chksum_v_3d(arrayV_in, 'v '//mesg, HI_in, symmetric=symmetric, &
                     scale=scale, logunit=logunit)
  endif
end subroutine chksum_uv_3d

!> Checksums a 2d array staggered at C-grid u points.
subroutine chksum_u_2d(array_m, mesg, HI_m, haloshift, symmetric, omit_corners, &
                       scale, logunit)
  type(hor_index_type),  target,   intent(in) :: HI_m     !< A horizontal index type
  real, dimension(HI_m%IsdB:,HI_m%jsd:), target, intent(in) :: array_m !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,               optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                           !! symmetric computational domain.
  logical,               optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A scaling factor for this array.
  integer, optional, intent(in) :: logunit !< IO unit for checksum logging

  real, pointer :: array(:,:)           ! Field array on the input grid
  real, allocatable, dimension(:,:) :: rescaled_array
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling
  integer :: iounit !< Log IO unit
  integer :: i, j, Is
  real :: aMean, aMin, aMax
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats
  integer :: turns                      ! Quarter turns from input to model grid

  ! Rotate array to the input grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    if (modulo(turns, 2) /= 0) then
      ! Arrays originating from v-points must be handled by vchksum
      allocate(array(HI%isd:HI%ied, HI%JsdB:HI%JedB))
      call rotate_array(array_m, -turns, array)
      call vchksum(array, mesg, HI, haloshift, symmetric, omit_corners, scale, logunit)
      return
    else
      allocate(array(HI%IsdB:HI%IedB, HI%jsd:HI%jed))
      call rotate_array(array_m, -turns, array)
    endif
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%jsc:HI%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0 ; if (present(scale)) scaling = scale
  iounit = error_unit; if(present(logunit)) iounit = logunit
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then
    if (present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2)) )
      rescaled_array(:,:) = 0.0
      Is = HI%isc ; if (sym_stats) Is = HI%isc-1
      do j=HI%jsc,HI%jec ; do I=Is,HI%IecB
        rescaled_array(I,j) = scale*array(I,j)
      enddo ; enddo
      call subStats(HI, rescaled_array, sym_stats, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, sym_stats, aMean, aMin, aMax)
    endif

    if (is_root_pe()) &
      call chk_sum_msg("u-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%iedB-HI%iecB

  if ( HI%iscB-hshift<HI%isdB .or. HI%iecB+hshift>HI%iedB .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_u_2d: haloshift =',hshift
    write(0,*) 'chksum_u_2d: isd,isc,iec,ied=',HI%isdB,HI%iscB,HI%iecB,HI%iedB
    write(0,*) 'chksum_u_2d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_u_2d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if ((hshift==0) .and. .not.sym) then
    if (is_root_pe()) call chk_sum_msg("u-point:", bc0, mesg, iounit)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (hshift==0) then
    bcW = subchk(array, HI, -hshift-1, 0, scaling)
    if (is_root_pe()) call chk_sum_msg_W("u-point:", bc0, bcW, mesg, iounit)
  elseif (do_corners) then
    if (sym) then
      bcSW = subchk(array, HI, -hshift-1, -hshift, scaling)
      bcNW = subchk(array, HI, -hshift-1, hshift, scaling)
    else
      bcSW = subchk(array, HI, -hshift, -hshift, scaling)
      bcNW = subchk(array, HI, -hshift, hshift, scaling)
    endif
    bcSE = subchk(array, HI, hshift, -hshift, scaling)
    bcNE = subchk(array, HI, hshift, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg("u-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
  else
    bcS = subchk(array, HI, 0, -hshift, scaling)
    bcE = subchk(array, HI, hshift, 0, scaling)
    if (sym) then
      bcW = subchk(array, HI, -hshift-1, 0, scaling)
    else
      bcW = subchk(array, HI, -hshift, 0, scaling)
    endif
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg_NSEW("u-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
  endif

  contains

  integer function subchk(array, HI, di, dj, scale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: array !< The array to be checksummed
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: scale !< A scaling factor for this array.
    integer :: i, j, bc
    subchk = 0
    ! This line deliberately uses the h-point computational domain.
    do j=HI%jsc+dj,HI%jec+dj ; do I=HI%isc+di,HI%iec+di
      bc = bitcount(abs(scale*array(I,j)))
      subchk = subchk + bc
    enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, sym_stats, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: array !< The array to be checksummed
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.
    real, intent(out) :: aMean !< Array mean
    real, intent(out) :: aMin !< Array minimum
    real, intent(out) :: aMax !< Array maximum

    integer :: i, j, n, IsB

    IsB = HI%isc ; if (sym_stats) IsB = HI%isc-1

    aMin = array(HI%isc,HI%jsc) ; aMax = aMin
    do j=HI%jsc,HI%jec ; do I=IsB,HI%IecB
      aMin = min(aMin, array(I,j))
      aMax = max(aMax, array(I,j))
    enddo ; enddo
    ! This line deliberately uses the h-point computational domain.
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec))
    n = (1 + HI%jec - HI%jsc) * (1 + HI%iec - HI%isc)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end subroutine chksum_u_2d

!> Checksums a 2d array staggered at C-grid v points.
subroutine chksum_v_2d(array_m, mesg, HI_m, haloshift, symmetric, omit_corners, &
                       scale, logunit)
  type(hor_index_type),  target,   intent(in) :: HI_m      !< A horizontal index type
  real, dimension(HI_m%isd:,HI_m%JsdB:), target, intent(in) :: array_m !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,               optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                           !! symmetric computational domain.
  logical,               optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                  optional, intent(in) :: scale     !< A scaling factor for this array.
  integer, optional, intent(in) :: logunit !< IO unit for checksum logging

  real, pointer :: array(:,:)           ! Field array on the input grid
  real, allocatable, dimension(:,:) :: rescaled_array
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling
  integer :: iounit !< Log IO unit
  integer :: i, j, Js
  real :: aMean, aMin, aMax
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats
  integer :: turns                      ! Quarter turns from input to model grid

  ! Rotate array to the input grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    if (modulo(turns, 2) /= 0) then
      ! Arrays originating from u-points must be handled by uchksum
      allocate(array(HI%IsdB:HI%IedB, HI%jsd:HI%jed))
      call rotate_array(array_m, -turns, array)
      call uchksum(array, mesg, HI, haloshift, symmetric, omit_corners, scale, logunit)
      return
    else
      allocate(array(HI%isd:HI%ied, HI%JsdB:HI%JedB))
      call rotate_array(array_m, -turns, array)
    endif
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%JscB:HI%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0 ; if (present(scale)) scaling = scale
  iounit = error_unit; if(present(logunit)) iounit = logunit
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then
    if (present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2)) )
      rescaled_array(:,:) = 0.0
      Js = HI%jsc ; if (sym_stats) Js = HI%jsc-1
      do J=Js,HI%JecB ; do i=HI%isc,HI%iec
        rescaled_array(i,J) = scale*array(i,J)
      enddo ; enddo
      call subStats(HI, rescaled_array, sym_stats, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, sym_stats, aMean, aMin, aMax)
    endif

    if (is_root_pe()) &
      call chk_sum_msg("v-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jscB-hshift<HI%jsdB .or. HI%jecB+hshift>HI%jedB ) then
    write(0,*) 'chksum_v_2d: haloshift =',hshift
    write(0,*) 'chksum_v_2d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_v_2d: jsd,jsc,jec,jed=',HI%jsdB,HI%jscB,HI%jecB,HI%jedB
    call chksum_error(FATAL,'Error in chksum_v_2d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if ((hshift==0) .and. .not.sym) then
    if (is_root_pe()) call chk_sum_msg("v-point:", bc0, mesg, iounit)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (hshift==0) then
    bcS = subchk(array, HI, 0, -hshift-1, scaling)
    if (is_root_pe()) call chk_sum_msg_S("v-point:", bc0, bcS, mesg, iounit)
  elseif (do_corners) then
    if (sym) then
      bcSW = subchk(array, HI, -hshift, -hshift-1, scaling)
      bcSE = subchk(array, HI, hshift, -hshift-1, scaling)
    else
      bcSW = subchk(array, HI, -hshift, -hshift, scaling)
      bcSE = subchk(array, HI, hshift, -hshift, scaling)
    endif
    bcNW = subchk(array, HI, -hshift, hshift, scaling)
    bcNE = subchk(array, HI, hshift, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg("v-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
  else
    if (sym) then
      bcS = subchk(array, HI, 0, -hshift-1, scaling)
    else
      bcS = subchk(array, HI, 0, -hshift, scaling)
    endif
    bcE = subchk(array, HI, hshift, 0, scaling)
    bcW = subchk(array, HI, -hshift, 0, scaling)
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg_NSEW("v-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
  endif

  contains

  integer function subchk(array, HI, di, dj, scale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%JsdB:), intent(in) :: array !< The array to be checksummed
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: scale !< A scaling factor for this array.
    integer :: i, j, bc
    subchk = 0
    ! This line deliberately uses the h-point computational domain.
    do J=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(scale*array(i,J)))
      subchk = subchk + bc
    enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, sym_stats, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%JsdB:), intent(in) :: array !< The array to be checksummed
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.
    real, intent(out) :: aMean !< Array mean
    real, intent(out) :: aMin !< Array minimum
    real, intent(out) :: aMax !< Array maximum

    integer :: i, j, n, JsB

    JsB = HI%jsc ; if (sym_stats) JsB = HI%jsc-1

    aMin = array(HI%isc,HI%jsc) ; aMax = aMin
    do J=JsB,HI%JecB ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,J))
      aMax = max(aMax, array(i,J))
    enddo ; enddo
    ! This line deliberately uses the h-computational domain.
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec))
    n = (1 + HI%jec - HI%jsc) * (1 + HI%iec - HI%isc)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end subroutine chksum_v_2d

!> Checksums a 3d array staggered at tracer points.
subroutine chksum_h_3d(array_m, mesg, HI_m, haloshift, omit_corners, scale, logunit)
  type(hor_index_type),    target,   intent(in) :: HI_m !< A horizontal index type
  real, dimension(HI_m%isd:,HI_m%jsd:,:), target, intent(in) :: array_m !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                 optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A scaling factor for this array.
  integer, optional, intent(in) :: logunit !< IO unit for checksum logging

  real, pointer :: array(:,:,:)         ! Field array on the input grid
  real, allocatable, dimension(:,:,:) :: rescaled_array
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling
  integer :: iounit !< Log IO unit
  integer :: i, j, k
  real :: aMean, aMin, aMax
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners
  integer :: turns                      ! Quarter turns from input to model grid

  ! Rotate array to the input grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    allocate(array(HI%isd:HI%ied, HI%jsd:HI%jed, size(array_m, 3)))
    call rotate_array(array_m, -turns, array)
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0 ; if (present(scale)) scaling = scale
  iounit = error_unit; if(present(logunit)) iounit = logunit

  if (calculateStatistics) then
    if (present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2), &
                               LBOUND(array,3):UBOUND(array,3)) )
      rescaled_array(:,:,:) = 0.0
      do k=1,size(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
        rescaled_array(i,j,k) = scale*array(i,j,k)
      enddo ; enddo ; enddo

      call subStats(HI, rescaled_array, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, aMean, aMin, aMax)
    endif

    if (is_root_pe()) &
      call chk_sum_msg("h-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_h_3d: haloshift =',hshift
    write(0,*) 'chksum_h_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_h_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_h_3d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  if (hshift==0) then
    if (is_root_pe()) call chk_sum_msg("h-point:", bc0, mesg, iounit)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (do_corners) then
    bcSW = subchk(array, HI, -hshift, -hshift, scaling)
    bcSE = subchk(array, HI, hshift, -hshift, scaling)
    bcNW = subchk(array, HI, -hshift, hshift, scaling)
    bcNE = subchk(array, HI, hshift, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg("h-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
  else
    bcS = subchk(array, HI, 0, -hshift, scaling)
    bcE = subchk(array, HI, hshift, 0, scaling)
    bcW = subchk(array, HI, -hshift, 0, scaling)
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg_NSEW("h-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
  endif

  contains

  integer function subchk(array, HI, di, dj, scale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: array !< The array to be checksummed
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: scale !< A scaling factor for this array.
    integer :: i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(scale*array(i,j,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: array !< The array to be checksummed
    real, intent(out) :: aMean !<  Array mean
    real, intent(out) :: aMin !< Array minimum
    real, intent(out) :: aMax !< Array maximum

    integer :: i, j, k, n

    aMin = array(HI%isc,HI%jsc,1)
    aMax = array(HI%isc,HI%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,j,k))
      aMax = max(aMax, array(i,j,k))
      n = n + 1
    enddo ; enddo ; enddo
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end subroutine chksum_h_3d

!> Checksums a 3d array staggered at corner points.
subroutine chksum_B_3d(array_m, mesg, HI_m, haloshift, symmetric, omit_corners, &
                       scale, logunit)
  type(hor_index_type),     target,   intent(in) :: HI_m !< A horizontal index type
  real, dimension(HI_m%IsdB:,HI_m%JsdB:,:), target, intent(in) :: array_m !< The array to be checksummed
  character(len=*),                   intent(in) :: mesg  !< An identifying message
  integer,                  optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                  optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                              !! symmetric computational domain.
  logical,                  optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                     optional, intent(in) :: scale     !< A scaling factor for this array.
  integer, optional, intent(in) :: logunit !< IO unit for checksum logging

  real, pointer :: array(:,:,:)         ! Field array on the input grid
  real, allocatable, dimension(:,:,:) :: rescaled_array
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling
  integer :: iounit !< Log IO unit
  integer :: i, j, k, Is, Js
  real :: aMean, aMin, aMax
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats
  integer :: turns                      ! Quarter turns from input to model grid

  ! Rotate array to the input grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    allocate(array(HI%IsdB:HI%IedB, HI%JsdB:HI%JedB, size(array_m, 3)))
    call rotate_array(array_m, -turns, array)
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%JscB:HI%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0 ; if (present(scale)) scaling = scale
  iounit = error_unit; if(present(logunit)) iounit = logunit
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then
    if (present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2), &
                               LBOUND(array,3):UBOUND(array,3)) )
      rescaled_array(:,:,:) = 0.0
      Is = HI%isc ; if (sym_stats) Is = HI%isc-1
      Js = HI%jsc ; if (sym_stats) Js = HI%jsc-1
      do k=1,size(array,3) ; do J=Js,HI%JecB ; do I=Is,HI%IecB
        rescaled_array(I,J,k) = scale*array(I,J,k)
      enddo ; enddo ; enddo
      call subStats(HI, rescaled_array, sym_stats, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, sym_stats, aMean, aMin, aMax)
    endif

    if (is_root_pe()) &
      call chk_sum_msg("B-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_B_3d: haloshift =',hshift
    write(0,*) 'chksum_B_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_B_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_B_3d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if ((hshift==0) .and. .not.sym) then
    if (is_root_pe()) call chk_sum_msg("B-point:", bc0, mesg, iounit)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (do_corners) then
    if (sym) then
      bcSW = subchk(array, HI, -hshift-1, -hshift-1, scaling)
      bcSE = subchk(array, HI, hshift, -hshift-1, scaling)
      bcNW = subchk(array, HI, -hshift-1, hshift, scaling)
    else
      bcSW = subchk(array, HI, -hshift-1, -hshift-1, scaling)
      bcSE = subchk(array, HI, hshift, -hshift-1, scaling)
      bcNW = subchk(array, HI, -hshift-1, hshift, scaling)
    endif
    bcNE = subchk(array, HI, hshift, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg("B-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
  else
    if (sym) then
      bcS = subchk(array, HI, 0, -hshift-1, scaling)
      bcW = subchk(array, HI, -hshift-1, 0, scaling)
    else
      bcS = subchk(array, HI, 0, -hshift, scaling)
      bcW = subchk(array, HI, -hshift, 0, scaling)
    endif
    bcE = subchk(array, HI, hshift, 0, scaling)
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg_NSEW("B-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
  endif

  contains

  integer function subchk(array, HI, di, dj, scale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: scale !< A scaling factor for this array.
    integer :: i, j, k, bc
    subchk = 0
    ! This line deliberately uses the h-point computational domain.
    do k=LBOUND(array,3),UBOUND(array,3) ; do J=HI%jsc+dj,HI%jec+dj ; do I=HI%isc+di,HI%iec+di
      bc = bitcount(abs(scale*array(I,J,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, sym_stats, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.
    real, intent(out) :: aMean !< Array mean
    real, intent(out) :: aMin !< Array minimum
    real, intent(out) :: aMax !< Array maximum

    integer :: i, j, k, n, IsB, JsB

    IsB = HI%isc ; if (sym_stats) IsB = HI%isc-1
    JsB = HI%jsc ; if (sym_stats) JsB = HI%jsc-1

    aMin = array(HI%isc,HI%jsc,1) ; aMax = aMin
    do k=LBOUND(array,3),UBOUND(array,3) ; do J=JsB,HI%JecB ; do I=IsB,HI%IecB
      aMin = min(aMin, array(I,J,k))
      aMax = max(aMax, array(I,J,k))
    enddo ; enddo ; enddo
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))
    n = (1 + HI%jec - HI%jsc) * (1 + HI%iec - HI%isc) * size(array,3)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end subroutine chksum_B_3d

!> Checksums a 3d array staggered at C-grid u points.
subroutine chksum_u_3d(array_m, mesg, HI_m, haloshift, symmetric, omit_corners, &
                       scale, logunit)
  type(hor_index_type),    target,   intent(in) :: HI_m !< A horizontal index type
  real, dimension(HI_m%isdB:,HI_m%Jsd:,:), target, intent(in) :: array_m !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                 optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                             !! symmetric computational domain.
  logical,                 optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A scaling factor for this array.
  integer, optional, intent(in) :: logunit !< IO unit for checksum logging

  real, pointer :: array(:,:,:)         ! Field array on the input grid
  real, allocatable, dimension(:,:,:) :: rescaled_array
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling
  integer :: iounit !< Log IO unit
  integer :: i, j, k, Is
  real :: aMean, aMin, aMax
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats
  integer :: turns                      ! Quarter turns from input to model grid

  ! Rotate array to the input grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    if (modulo(turns, 2) /= 0) then
      ! Arrays originating from v-points must be handled by vchksum
      allocate(array(HI%isd:HI%ied, HI%JsdB:HI%JedB, size(array_m, 3)))
      call rotate_array(array_m, -turns, array)
      call vchksum(array, mesg, HI, haloshift, symmetric, omit_corners, scale, logunit)
      return
    else
      allocate(array(HI%IsdB:HI%IedB, HI%jsd:HI%jed, size(array_m, 3)))
      call rotate_array(array_m, -turns, array)
    endif
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%jsc:HI%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0 ; if (present(scale)) scaling = scale
  iounit = error_unit; if(present(logunit)) iounit = logunit
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then
    if (present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2), &
                               LBOUND(array,3):UBOUND(array,3)) )
      rescaled_array(:,:,:) = 0.0
      Is = HI%isc ; if (sym_stats) Is = HI%isc-1
      do k=1,size(array,3) ; do j=HI%jsc,HI%jec ; do I=Is,HI%IecB
        rescaled_array(I,j,k) = scale*array(I,j,k)
      enddo ; enddo ; enddo
      call subStats(HI, rescaled_array, sym_stats, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, sym_stats, aMean, aMin, aMax)
    endif
    if (is_root_pe()) &
      call chk_sum_msg("u-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_u_3d: haloshift =',hshift
    write(0,*) 'chksum_u_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_u_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_u_3d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if ((hshift==0) .and. .not.sym) then
    if (is_root_pe()) call chk_sum_msg("u-point:", bc0, mesg, iounit)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (hshift==0) then
    bcW = subchk(array, HI, -hshift-1, 0, scaling)
    if (is_root_pe()) call chk_sum_msg_W("u-point:", bc0, bcW, mesg, iounit)
  elseif (do_corners) then
    if (sym) then
      bcSW = subchk(array, HI, -hshift-1, -hshift, scaling)
      bcNW = subchk(array, HI, -hshift-1, hshift, scaling)
    else
      bcSW = subchk(array, HI, -hshift, -hshift, scaling)
      bcNW = subchk(array, HI, -hshift, hshift, scaling)
    endif
    bcSE = subchk(array, HI, hshift, -hshift, scaling)
    bcNE = subchk(array, HI, hshift, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg("u-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
  else
    bcS = subchk(array, HI, 0, -hshift, scaling)
    bcE = subchk(array, HI, hshift, 0, scaling)
    if (sym) then
      bcW = subchk(array, HI, -hshift-1, 0, scaling)
    else
      bcW = subchk(array, HI, -hshift, 0, scaling)
    endif
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg_NSEW("u-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
  endif

  contains

  integer function subchk(array, HI, di, dj, scale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%jsd:,:), intent(in) :: array !< The array to be checksummed
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: scale !< A scaling factor for this array.
    integer :: i, j, k, bc
    subchk = 0
    ! This line deliberately uses the h-point computational domain.
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc+dj,HI%jec+dj ; do I=HI%isc+di,HI%iec+di
      bc = bitcount(abs(scale*array(I,j,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  subroutine subStats(HI, array, sym_stats, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%jsd:,:), intent(in) :: array !< The array to be checksummed
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.
    real, intent(out) :: aMean !< Array mean
    real, intent(out) :: aMin !< Array minimum
    real, intent(out) :: aMax !< Array maximum

    integer :: i, j, k, n, IsB

    IsB = HI%isc ; if (sym_stats) IsB = HI%isc-1

    aMin = array(HI%isc,HI%jsc,1) ; aMax = aMin
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do I=IsB,HI%IecB
      aMin = min(aMin, array(I,j,k))
      aMax = max(aMax, array(I,j,k))
    enddo ; enddo ; enddo
    ! This line deliberately uses the h-point computational domain.
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))
    n = (1 + HI%jec - HI%jsc) * (1 + HI%iec - HI%isc) * size(array,3)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end subroutine chksum_u_3d

!> Checksums a 3d array staggered at C-grid v points.
subroutine chksum_v_3d(array_m, mesg, HI_m, haloshift, symmetric, omit_corners, &
                       scale, logunit)
  type(hor_index_type),    target,   intent(in) :: HI_m !< A horizontal index type
  real, dimension(HI_m%isd:,HI_m%JsdB:,:), target, intent(in) :: array_m !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                 optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                             !! symmetric computational domain.
  logical,                 optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A scaling factor for this array.
  integer, optional, intent(in) :: logunit !< IO unit for checksum logging

  real, pointer :: array(:,:,:)         ! Field array on the input grid
  real, allocatable, dimension(:,:,:) :: rescaled_array
  type(hor_index_type), pointer :: HI   ! Horizontal index bounds of the input grid
  real :: scaling
  integer :: iounit !< Log IO unit
  integer :: i, j, k, Js
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  real :: aMean, aMin, aMax
  logical :: do_corners, sym, sym_stats
  integer :: turns                      ! Quarter turns from input to model grid

  ! Rotate array to the input grid
  turns = HI_m%turns
  if (modulo(turns, 4) /= 0) then
    allocate(HI)
    call rotate_hor_index(HI_m, -turns, HI)
    if (modulo(turns, 2) /= 0) then
      ! Arrays originating from u-points must be handled by uchksum
      allocate(array(HI%IsdB:HI%IedB, HI%jsd:HI%jed, size(array_m, 3)))
      call rotate_array(array_m, -turns, array)
      call uchksum(array, mesg, HI, haloshift, symmetric, omit_corners, scale, logunit)
      return
    else
      allocate(array(HI%isd:HI%ied, HI%JsdB:HI%JedB, size(array_m, 3)))
      call rotate_array(array_m, -turns, array)
    endif
  else
    HI => HI_m
    array => array_m
  endif

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%JscB:HI%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  scaling = 1.0 ; if (present(scale)) scaling = scale
  iounit = error_unit; if(present(logunit)) iounit = logunit
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then
    if (present(scale)) then
      allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                               LBOUND(array,2):UBOUND(array,2), &
                               LBOUND(array,3):UBOUND(array,3)) )
      rescaled_array(:,:,:) = 0.0
      Js = HI%jsc ; if (sym_stats) Js = HI%jsc-1
      do k=1,size(array,3) ; do J=Js,HI%JecB ; do i=HI%isc,HI%iec
        rescaled_array(i,J,k) = scale*array(i,J,k)
      enddo ; enddo ; enddo
      call subStats(HI, rescaled_array, sym_stats, aMean, aMin, aMax)
      deallocate(rescaled_array)
    else
      call subStats(HI, array, sym_stats, aMean, aMin, aMax)
    endif
    if (is_root_pe()) &
      call chk_sum_msg("v-point:", aMean, aMin, aMax, mesg, iounit)
  endif

  if (.not.writeChksums) return

  hshift = default_shift
  if (present(haloshift)) hshift = haloshift
  if (hshift<0) hshift = HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_v_3d: haloshift =',hshift
    write(0,*) 'chksum_v_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_v_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_v_3d '//trim(mesg))
  endif

  bc0 = subchk(array, HI, 0, 0, scaling)

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if ((hshift==0) .and. .not.sym) then
    if (is_root_pe()) call chk_sum_msg("v-point:", bc0, mesg, iounit)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (hshift==0) then
    bcS = subchk(array, HI, 0, -hshift-1, scaling)
    if (is_root_pe()) call chk_sum_msg_S("v-point:", bc0, bcS, mesg, iounit)
  elseif (do_corners) then
    if (sym) then
      bcSW = subchk(array, HI, -hshift, -hshift-1, scaling)
      bcSE = subchk(array, HI, hshift, -hshift-1, scaling)
    else
      bcSW = subchk(array, HI, -hshift, -hshift, scaling)
      bcSE = subchk(array, HI, hshift, -hshift, scaling)
    endif
    bcNW = subchk(array, HI, -hshift, hshift, scaling)
    bcNE = subchk(array, HI, hshift, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg("v-point:", bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
  else
    if (sym) then
      bcS = subchk(array, HI, 0, -hshift-1, scaling)
    else
      bcS = subchk(array, HI, 0, -hshift, scaling)
    endif
    bcE = subchk(array, HI, hshift, 0, scaling)
    bcW = subchk(array, HI, -hshift, 0, scaling)
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) &
      call chk_sum_msg_NSEW("v-point:", bc0, bcN, bcS, bcE, bcW, mesg, iounit)
  endif

  contains

  integer function subchk(array, HI, di, dj, scale)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
    integer, intent(in) :: di    !< i- direction array shift for this checksum
    integer, intent(in) :: dj    !< j- direction array shift for this checksum
    real, intent(in)    :: scale !< A scaling factor for this array.
    integer :: i, j, k, bc
    subchk = 0
    ! This line deliberately uses the h-point computational domain.
    do k=LBOUND(array,3),UBOUND(array,3) ; do J=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(scale*array(i,J,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk, bc_modulus)
  end function subchk

  !subroutine subStats(HI, array, mesg, sym_stats)
  subroutine subStats(HI, array, sym_stats, aMean, aMin, aMax)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.
    real, intent(out) :: aMean   !< Mean of array over domain
    real, intent(out) :: aMin    !< Minimum of array over domain
    real, intent(out) :: aMax    !< Maximum of array over domain

    integer :: i, j, k, n, JsB

    JsB = HI%jsc ; if (sym_stats) JsB = HI%jsc-1

    aMin = array(HI%isc,HI%jsc,1) ; aMax = aMin
    do k=LBOUND(array,3),UBOUND(array,3) ; do J=JsB,HI%JecB ; do i=HI%isc,HI%iec
      aMin = min(aMin, array(i,J,k))
      aMax = max(aMax, array(i,J,k))
    enddo ; enddo ; enddo
    ! This line deliberately uses the h-point computational domain.
    aMean = reproducing_sum(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))
    n = (1 + HI%jec - HI%jsc) * (1 + HI%iec - HI%isc) * size(array,3)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
  end subroutine subStats

end subroutine chksum_v_3d

!   These are the older version of chksum that do not take the grid staggering
! into account.

!> chksum1d does a checksum of a 1-dimensional array.
subroutine chksum1d(array, mesg, start_i, end_i, compare_PEs)
  real, dimension(:), intent(in) :: array   !< The array to be summed (index starts at 1).
  character(len=*),   intent(in) :: mesg    !< An identifying message.
  integer, optional,  intent(in) :: start_i !< The starting index for the sum (default 1)
  integer, optional,  intent(in) :: end_i   !< The ending index for the sum (default all)
  logical, optional,  intent(in) :: compare_PEs !< If true, compare across PEs instead of summing
                                                !! and list the root_PE value (default true)

  integer :: is, ie, i, bc, sum1, sum_bc
  real :: sum
  real, allocatable :: sum_here(:)
  logical :: compare
  integer :: pe_num   ! pe number of the data
  integer :: nPEs     ! Total number of processsors

  is = LBOUND(array,1) ; ie = UBOUND(array,1)
  if (present(start_i)) is = start_i
  if (present(end_i)) ie = end_i
  compare = .true. ; if (present(compare_PEs)) compare = compare_PEs

  sum = 0.0 ; sum_bc = 0
  do i=is,ie
    sum = sum + array(i)
    bc = bitcount(ABS(array(i)))
    sum_bc = sum_bc + bc
  enddo

  pe_num = pe_here() + 1 - root_pe() ; nPEs = num_pes()
  allocate(sum_here(nPEs)) ; sum_here(:) = 0.0 ; sum_here(pe_num) = sum
  call sum_across_PEs(sum_here,nPEs)

  sum1 = sum_bc
  call sum_across_PEs(sum1)

  if (.not.compare) then
    sum = 0.0
    do i=1,nPEs ; sum = sum + sum_here(i) ; enddo
    sum_bc = sum1
  elseif (is_root_pe()) then
    if (sum1 /= nPEs*sum_bc) &
      write(0, '(A40," bitcounts do not match across PEs: ",I12,1X,I12)') &
            mesg, sum1, nPEs*sum_bc
    do i=1,nPEs ; if (sum /= sum_here(i)) then
      write(0, '(A40," PE ",i4," sum mismatches root_PE: ",3(ES22.13,1X))') &
            mesg, i, sum_here(i), sum, sum_here(i)-sum
    endif ; enddo
  endif
  deallocate(sum_here)

  if (is_root_pe()) &
    write(0,'(A50,1X,ES25.16,1X,I12)') mesg, sum, sum_bc

end subroutine chksum1d

!   These are the older version of chksum that do not take the grid staggering
! into account.

!> chksum2d does a checksum of all data in a 2-d array.
subroutine chksum2d(array, mesg)

  real, dimension(:,:) :: array !< The array to be checksummed
  character(len=*) :: mesg  !< An identifying message

  integer :: xs,xe,ys,ye,i,j,sum1,bc
  real :: sum

  xs = LBOUND(array,1) ; xe = UBOUND(array,1)
  ys = LBOUND(array,2) ; ye = UBOUND(array,2)

  sum = 0.0 ; sum1 = 0
  do i=xs,xe ; do j=ys,ye
    bc = bitcount(abs(array(i,j)))
    sum1 = sum1 + bc
  enddo ; enddo
  call sum_across_PEs(sum1)

  sum = reproducing_sum(array(:,:))

  if (is_root_pe()) &
    write(0,'(A50,1X,ES25.16,1X,I12)') mesg, sum, sum1
!    write(0,'(A40,1X,Z16.16,1X,Z16.16,1X,ES25.16,1X,I12)') &
!      mesg, sum, sum1, sum, sum1

end subroutine chksum2d

!> chksum3d does a checksum of all data in a 2-d array.
subroutine chksum3d(array, mesg)

  real, dimension(:,:,:) :: array !< The array to be checksummed
  character(len=*) :: mesg  !< An identifying message

  integer :: xs,xe,ys,ye,zs,ze,i,j,k, bc,sum1
  real :: sum

  xs = LBOUND(array,1) ; xe = UBOUND(array,1)
  ys = LBOUND(array,2) ; ye = UBOUND(array,2)
  zs = LBOUND(array,3) ; ze = UBOUND(array,3)

  sum = 0.0 ; sum1 = 0
  do i=xs,xe ; do j=ys,ye ; do k=zs,ze
    bc = bitcount(ABS(array(i,j,k)))
    sum1 = sum1 + bc
  enddo ; enddo ; enddo

  call sum_across_PEs(sum1)
  sum = reproducing_sum(array(:,:,:))

  if (is_root_pe()) &
    write(0,'(A50,1X,ES25.16,1X,I12)') mesg, sum, sum1
!    write(0,'(A40,1X,Z16.16,1X,Z16.16,1X,ES25.16,1X,I12)') &
!      mesg, sum, sum1, sum, sum1

end subroutine chksum3d

!> This function returns .true. if x is a NaN, and .false. otherwise.
function is_NaN_0d(x)
  real, intent(in) :: x !< The value to be checked for NaNs.
  logical :: is_NaN_0d

 !is_NaN_0d = (((x < 0.0) .and. (x >= 0.0)) .or. &
 !          (.not.(x < 0.0) .and. .not.(x >= 0.0)))
  if (((x < 0.0) .and. (x >= 0.0)) .or. &
            (.not.(x < 0.0) .and. .not.(x >= 0.0))) then
    is_NaN_0d = .true.
  else
    is_NaN_0d = .false.
  endif

end function is_NaN_0d

!> Returns .true. if any element of x is a NaN, and .false. otherwise.
function is_NaN_1d(x, skip_mpp)
  real, dimension(:), intent(in) :: x !< The array to be checked for NaNs.
  logical,  optional, intent(in) :: skip_mpp  !< If true, only check this array only
                                              !! on the local PE (default false).
  logical :: is_NaN_1d

  integer :: i, n
  logical :: global_check

  n = 0
  do i = LBOUND(x,1), UBOUND(x,1)
    if (is_NaN_0d(x(i))) n = n + 1
  enddo
  global_check = .true.
  if (present(skip_mpp)) global_check = .not.skip_mpp

  if (global_check) call sum_across_PEs(n)
  is_NaN_1d = .false.
  if (n>0) is_NaN_1d = .true.

end function is_NaN_1d

!> Returns .true. if any element of x is a NaN, and .false. otherwise.
function is_NaN_2d(x)
  real, dimension(:,:), intent(in) :: x !< The array to be checked for NaNs.
  logical :: is_NaN_2d

  integer :: i, j, n

  n = 0
  do j = LBOUND(x,2), UBOUND(x,2) ; do i = LBOUND(x,1), UBOUND(x,1)
    if (is_NaN_0d(x(i,j))) n = n + 1
  enddo ; enddo
  call sum_across_PEs(n)
  is_NaN_2d = .false.
  if (n>0) is_NaN_2d = .true.

end function is_NaN_2d

!> Returns .true. if any element of x is a NaN, and .false. otherwise.
function is_NaN_3d(x)
  real, dimension(:,:,:), intent(in) :: x !< The array to be checked for NaNs.
  logical :: is_NaN_3d

  integer :: i, j, k, n

  n = 0
  do k = LBOUND(x,3), UBOUND(x,3)
    do j = LBOUND(x,2), UBOUND(x,2) ; do i = LBOUND(x,1), UBOUND(x,1)
      if (is_NaN_0d(x(i,j,k))) n = n + 1
    enddo ; enddo
  enddo
  call sum_across_PEs(n)
  is_NaN_3d = .false.
  if (n>0) is_NaN_3d = .true.

end function is_NaN_3d

! The following set of routines do a checksum across the computational domain of
! a field, with the potential for rotation of this field and masking.

!> Compute the field checksum of a scalar.
function rotated_field_chksum_real_0d(field, pelist, mask_val, turns) &
    result(chksum)
  real,              intent(in) :: field      !< Input scalar
  integer, optional, intent(in) :: pelist(:)  !< PE list of ranks to checksum
  real,    optional, intent(in) :: mask_val   !< FMS mask value
  integer, optional, intent(in) :: turns      !< Number of quarter turns
  integer(kind=int64) :: chksum               !< checksum of scalar

  if (present(turns)) call MOM_error(FATAL, "Rotation not supported for 0d fields.")

  chksum = field_chksum(field, pelist=pelist, mask_val=mask_val)
end function rotated_field_chksum_real_0d


!> Compute the field checksum of a 1d field.
function rotated_field_chksum_real_1d(field, pelist, mask_val, turns) &
    result(chksum)
  real, dimension(:), intent(in) :: field     !< Input array
  integer,  optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,     optional, intent(in) :: mask_val  !< FMS mask value
  integer,  optional, intent(in) :: turns     !< Number of quarter turns
  integer(kind=int64) :: chksum               !< checksum of array

  if (present(turns)) call MOM_error(FATAL, "Rotation not supported for 1d fields.")

  chksum = field_chksum(field, pelist=pelist, mask_val=mask_val)
end function rotated_field_chksum_real_1d


!> Compute the field checksum of a rotated 2d field.
function rotated_field_chksum_real_2d(field, pelist, mask_val, turns) &
    result(chksum)
  real, dimension(:,:),     intent(in) :: field     !< Unrotated input field
  integer,        optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,           optional, intent(in) :: mask_val  !< FMS mask value
  integer,        optional, intent(in) :: turns     !< Number of quarter turns
  integer(kind=int64) :: chksum                     !< checksum of array

  ! Local variables
  real, allocatable :: field_rot(:,:)  ! A rotated version of field, with the same units
  integer :: qturns ! The number of quarter turns through which to rotate field

  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)

  if (qturns == 0) then
    chksum = field_chksum(field, pelist=pelist, mask_val=mask_val)
  else
    call allocate_rotated_array(field, [1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    chksum = field_chksum(field_rot, pelist=pelist, mask_val=mask_val)
    deallocate(field_rot)
  endif
end function rotated_field_chksum_real_2d

!> Compute the field checksum of a rotated 3d field.
function rotated_field_chksum_real_3d(field, pelist, mask_val, turns) &
    result(chksum)
  real, dimension(:,:,:),   intent(in) :: field     !< Unrotated input field
  integer,        optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,           optional, intent(in) :: mask_val  !< FMS mask value
  integer,        optional, intent(in) :: turns     !< Number of quarter turns
  integer(kind=int64) :: chksum                     !< checksum of array

  ! Local variables
  real, allocatable :: field_rot(:,:,:)  ! A rotated version of field, with the same units
  integer :: qturns ! The number of quarter turns through which to rotate field

  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)

  if (qturns == 0) then
    chksum = field_chksum(field, pelist=pelist, mask_val=mask_val)
  else
    call allocate_rotated_array(field, [1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    chksum = field_chksum(field_rot, pelist=pelist, mask_val=mask_val)
    deallocate(field_rot)
  endif
end function rotated_field_chksum_real_3d

!> Compute the field checksum of a rotated 4d field.
function rotated_field_chksum_real_4d(field, pelist, mask_val, turns) &
    result(chksum)
  real, dimension(:,:,:,:), intent(in) :: field     !< Unrotated input field
  integer,        optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,           optional, intent(in) :: mask_val  !< FMS mask value
  integer,        optional, intent(in) :: turns     !< Number of quarter turns
  integer(kind=int64) :: chksum                     !< checksum of array

  ! Local variables
  real, allocatable :: field_rot(:,:,:,:)  ! A rotated version of field, with the same units
  integer :: qturns ! The number of quarter turns through which to rotate field

  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)

  if (qturns == 0) then
    chksum = field_chksum(field, pelist=pelist, mask_val=mask_val)
  else
    call allocate_rotated_array(field, [1,1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    chksum = field_chksum(field_rot, pelist=pelist, mask_val=mask_val)
    deallocate(field_rot)
  endif
end function rotated_field_chksum_real_4d


!> Write a message including the checksum of the non-shifted array
subroutine chk_sum_msg1(fmsg, bc0, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: iounit !< Checksum logger IO unit

  if (is_root_pe()) &
    write(iounit, '(A,1(A,I10,X),A)') fmsg, " c=", bc0, trim(mesg)
end subroutine chk_sum_msg1

!> Write a message including checksums of non-shifted and diagonally shifted arrays
subroutine chk_sum_msg5(fmsg, bc0, bcSW, bcSE, bcNW, bcNE, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcSW !< The bitcount for SW shifted array
  integer,          intent(in) :: bcSE !< The bitcount for SE shifted array
  integer,          intent(in) :: bcNW !< The bitcount for NW shifted array
  integer,          intent(in) :: bcNE !< The bitcount for NE shifted array
  integer,          intent(in) :: iounit !< Checksum logger IO unit

  if (is_root_pe()) write(iounit, '(A,5(A,I10,1X),A)') &
    fmsg, " c=", bc0, "sw=", bcSW, "se=", bcSE, "nw=", bcNW, "ne=", bcNE, trim(mesg)
end subroutine chk_sum_msg5

!> Write a message including checksums of non-shifted and laterally shifted arrays
subroutine chk_sum_msg_NSEW(fmsg, bc0, bcN, bcS, bcE, bcW, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcN !< The bitcount for N shifted array
  integer,          intent(in) :: bcS !< The bitcount for S shifted array
  integer,          intent(in) :: bcE !< The bitcount for E shifted array
  integer,          intent(in) :: bcW !< The bitcount for W shifted array
  integer,          intent(in) :: iounit !< Checksum logger IO unit

  if (is_root_pe()) write(iounit, '(A,5(A,I10,1X),A)') &
    fmsg, " c=", bc0, "N=", bcN, "S=", bcS, "E=", bcE, "W=", bcW, trim(mesg)
end subroutine chk_sum_msg_NSEW

!> Write a message including checksums of non-shifted and southward shifted arrays
subroutine chk_sum_msg_S(fmsg, bc0, bcS, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcS  !< The bitcount of the south-shifted array
  integer,          intent(in) :: iounit !< Checksum logger IO unit

  if (is_root_pe()) write(iounit, '(A,2(A,I10,1X),A)') &
    fmsg, " c=", bc0, "S=", bcS, trim(mesg)
end subroutine chk_sum_msg_S

!> Write a message including checksums of non-shifted and westward shifted arrays
subroutine chk_sum_msg_W(fmsg, bc0, bcW, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcW  !< The bitcount of the west-shifted array
  integer,          intent(in) :: iounit !< Checksum logger IO unit

  if (is_root_pe()) write(iounit, '(A,2(A,I10,1X),A)') &
    fmsg, " c=", bc0, "W=", bcW, trim(mesg)
end subroutine chk_sum_msg_W

!> Write a message including checksums of non-shifted and southwestward shifted arrays
subroutine chk_sum_msg2(fmsg, bc0, bcSW, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcSW !< The bitcount of the southwest-shifted array
  integer,          intent(in) :: iounit !< Checksum logger IO unit

  if (is_root_pe()) write(iounit, '(A,2(A,I9,1X),A)') &
    fmsg, " c=", bc0, "s/w=", bcSW, trim(mesg)
end subroutine chk_sum_msg2

!> Write a message including the global mean, maximum and minimum of an array
subroutine chk_sum_msg3(fmsg, aMean, aMin, aMax, mesg, iounit)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  real,             intent(in) :: aMean !< The mean value of the array
  real,             intent(in) :: aMin !< The minimum value of the array
  real,             intent(in) :: aMax !< The maximum value of the array
  integer,          intent(in) :: iounit !< Checksum logger IO unit

  ! NOTE: We add zero to aMin and aMax to remove any negative zeros.
  ! This is due to inconsistencies of signed zero in local vs MPI calculations.

  if (is_root_pe()) write(iounit, '(A,3(A,ES25.16,1X),A)') &
    fmsg, " mean=", aMean, "min=", (0. + aMin), "max=", (0. + aMax), trim(mesg)
end subroutine chk_sum_msg3

!> MOM_checksums_init initializes the MOM_checksums module. As it happens, the
!! only thing that it does is to log the version of this module.
subroutine MOM_checksums_init(param_file)
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_checksums" ! This module's name.

  call log_version(param_file, mdl, version)

end subroutine MOM_checksums_init

!> A wrapper for MOM_error used in the checksum code
subroutine chksum_error(signal, message)
  ! Wrapper for MOM_error to help place specific break points in debuggers
  integer, intent(in) :: signal !< An error severity level, such as FATAL or WARNING
  character(len=*), intent(in) :: message !< An error message
  call MOM_error(signal, message)
end subroutine chksum_error

!> Does a bitcount of a number by first casting to an integer and then using BTEST
!! to check bit by bit
integer function bitcount(x)
  real, intent(in) :: x !< Number to be bitcount

  integer, parameter :: xk = kind(x)  !< Kind type of x

  ! NOTE: Assumes that reals and integers of kind=xk are the same size
  bitcount = popcnt(transfer(x, 1_xk))
end function bitcount

end module MOM_checksums
