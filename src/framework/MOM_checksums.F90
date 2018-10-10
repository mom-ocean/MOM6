!> Routines to calculate checksums of various array and vector types
module MOM_checksums

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms, only : PE_here, root_PE, num_PEs, sum_across_PEs
use MOM_coms, only : min_across_PEs, max_across_PEs
use MOM_coms, only : reproducing_sum
use MOM_error_handler, only : MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : log_version, param_file_type
use MOM_hor_index, only : hor_index_type

implicit none ; private

public :: hchksum, Bchksum, uchksum, vchksum, qchksum, is_NaN, chksum
public :: hchksum_pair, uvchksum, Bchksum_pair
public :: chksum_general
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

!> Return the bitcount of an array
interface chksum_general
  module procedure chksum_general_1d, chksum_general_2d, chksum_general_3d
end interface

integer, parameter :: default_shift=0 !< The default array shift
logical :: calculateStatistics=.true. !< If true, report min, max and mean.
logical :: writeChksums=.true. !< If true, report the bitcount checksum
logical :: checkForNaNs=.true. !< If true, checks array for NaNs and cause
                               !! FATAL error is any are found

contains

!> Checksums on a pair of 2d arrays staggered at tracer points.
subroutine chksum_pair_h_2d(mesg, arrayA, arrayB, HI, haloshift, omit_corners, scale)
  character(len=*),                 intent(in) :: mesg !< Identifying messages
  type(hor_index_type),             intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:), intent(in) :: arrayA !< The first array to be checksummed
  real, dimension(HI%isd:,HI%jsd:), intent(in) :: arrayB !< The second array to be checksummed
  integer,                optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                   optional, intent(in) :: scale     !< A scaling factor for this array.

  if (present(haloshift)) then
    call chksum_h_2d(arrayA, 'x '//mesg, HI, haloshift, omit_corners, scale=scale)
    call chksum_h_2d(arrayB, 'y '//mesg, HI, haloshift, omit_corners, scale=scale)
  else
    call chksum_h_2d(arrayA, 'x '//mesg, HI, scale=scale)
    call chksum_h_2d(arrayB, 'y '//mesg, HI, scale=scale)
  endif

end subroutine chksum_pair_h_2d

!> Checksums on a pair of 3d arrays staggered at tracer points.
subroutine chksum_pair_h_3d(mesg, arrayA, arrayB, HI, haloshift, omit_corners, scale)
  character(len=*),                    intent(in) :: mesg !< Identifying messages
  type(hor_index_type),                intent(in) :: HI   !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:, :), intent(in) :: arrayA !< The first array to be checksummed
  real, dimension(HI%isd:,HI%jsd:, :), intent(in) :: arrayB !< The second array to be checksummed
  integer,                   optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                   optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                      optional, intent(in) :: scale     !< A scaling factor for this array.

  if (present(haloshift)) then
    call chksum_h_3d(arrayA, 'x '//mesg, HI, haloshift, omit_corners, scale=scale)
    call chksum_h_3d(arrayB, 'y '//mesg, HI, haloshift, omit_corners, scale=scale)
  else
    call chksum_h_3d(arrayA, 'x '//mesg, HI, scale=scale)
    call chksum_h_3d(arrayB, 'y '//mesg, HI, scale=scale)
  endif

end subroutine chksum_pair_h_3d

!> Checksums a 2d array staggered at tracer points.
subroutine chksum_h_2d(array, mesg, HI, haloshift, omit_corners, scale)
  type(hor_index_type),            intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:), intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,               optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                  optional, intent(in) :: scale     !< A scaling factor for this array.

  real, allocatable, dimension(:,:) :: rescaled_array
  real :: scaling
  integer :: i, j
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%jsc:HI%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif
  scaling = 1.0 ; if (present(scale)) scaling = scale

  if (calculateStatistics) then ; if (present(scale)) then
    allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                             LBOUND(array,2):UBOUND(array,2)) )
    rescaled_array(:,:) = 0.0
    do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      rescaled_array(i,j) = scale*array(i,j)
    enddo ; enddo
    call subStats(HI, rescaled_array, mesg)
    deallocate(rescaled_array)
  else
    call subStats(HI, array, mesg)
  endif ; endif

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
    if (is_root_pe()) call chk_sum_msg("h-point:",bc0,mesg)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (do_corners) then
    bcSW = subchk(array, HI, -hshift, -hshift, scaling)
    bcSE = subchk(array, HI, hshift, -hshift, scaling)
    bcNW = subchk(array, HI, -hshift, hshift, scaling)
    bcNE = subchk(array, HI, hshift, hshift, scaling)

    if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  else
    bcS = subchk(array, HI, 0, -hshift, scaling)
    bcE = subchk(array, HI, hshift, 0, scaling)
    bcW = subchk(array, HI, -hshift, 0, scaling)
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) call chk_sum_msg_NSEW("h-point:",bc0,bcN,bcS,bcE,bcW,mesg)
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
    do j=HI%jsc+dj,HI%jec+dj; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(scale*array(i,j)))
      subchk = subchk + bc
    enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%jsd:), intent(in) :: array !< The array to be checksummed
    character(len=*), intent(in) :: mesg  !< An identifying message

    integer :: i, j, n
    real :: aMean, aMin, aMax

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
    if (is_root_pe()) call chk_sum_msg("h-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_h_2d

!> Checksums on a pair of 2d arrays staggered at q-points.
subroutine chksum_pair_B_2d(mesg, arrayA, arrayB, HI, haloshift, symmetric, omit_corners, scale)
  character(len=*),                 intent(in) :: mesg   !< Identifying messages
  type(hor_index_type),             intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:), intent(in) :: arrayA !< The first array to be checksummed
  real, dimension(HI%isd:,HI%jsd:), intent(in) :: arrayB !< The second array to be checksummed
  logical,                optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                            !! symmetric computational domain.
  integer,                optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                   optional, intent(in) :: scale     !< A scaling factor for this array.

  logical :: sym

  sym = .false. ; if (present(symmetric)) sym = symmetric

  if (present(haloshift)) then
    call chksum_B_2d(arrayA, 'x '//mesg, HI, haloshift, symmetric=sym, &
                     omit_corners=omit_corners, scale=scale)
    call chksum_B_2d(arrayB, 'y '//mesg, HI, haloshift, symmetric=sym, &
                     omit_corners=omit_corners, scale=scale)
  else
    call chksum_B_2d(arrayA, 'x '//mesg, HI, symmetric=sym, scale=scale)
    call chksum_B_2d(arrayB, 'y '//mesg, HI, symmetric=sym, scale=scale)
  endif

end subroutine chksum_pair_B_2d

!> Checksums on a pair of 3d arrays staggered at q-points.
subroutine chksum_pair_B_3d(mesg, arrayA, arrayB, HI, haloshift, symmetric, omit_corners, scale)
  character(len=*),                    intent(in) :: mesg !< Identifying messages
  type(hor_index_type),                intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%JsdB:, :), intent(in) :: arrayA !< The first array to be checksummed
  real, dimension(HI%IsdB:,HI%JsdB:, :), intent(in) :: arrayB !< The second array to be checksummed
  integer,                   optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                   optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                               !! symmetric computational domain.
  logical,                   optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                      optional, intent(in) :: scale     !< A scaling factor for this array.

  logical :: sym

  if (present(haloshift)) then
    call chksum_B_3d(arrayA, 'x '//mesg, HI, haloshift, symmetric, &
                     omit_corners, scale=scale)
    call chksum_B_3d(arrayB, 'y '//mesg, HI, haloshift, symmetric, &
                     omit_corners, scale=scale)
  else
    call chksum_B_3d(arrayA, 'x '//mesg, HI, symmetric=symmetric, scale=scale)
    call chksum_B_3d(arrayB, 'y '//mesg, HI, symmetric=symmetric, scale=scale)
  endif

end subroutine chksum_pair_B_3d

!> Checksums a 2d array staggered at corner points.
subroutine chksum_B_2d(array, mesg, HI, haloshift, symmetric, omit_corners, scale)
  type(hor_index_type), intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%JsdB:), &
                        intent(in) :: array !< The array to be checksummed
  character(len=*),     intent(in) :: mesg  !< An identifying message
  integer,    optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,    optional, intent(in) :: symmetric !< If true, do the checksums on the
                                                !! full symmetric computational domain.
  logical,    optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,       optional, intent(in) :: scale     !< A scaling factor for this array.

  real, allocatable, dimension(:,:) :: rescaled_array
  real :: scaling
  integer :: i, j, Is, Js
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%JscB:HI%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif
  scaling = 1.0 ; if (present(scale)) scaling = scale
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then ; if (present(scale)) then
    allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                             LBOUND(array,2):UBOUND(array,2)) )
    rescaled_array(:,:) = 0.0
    Is = HI%isc ; if (sym_stats) Is = HI%isc-1
    Js = HI%jsc ; if (sym_stats) Js = HI%jsc-1
    do J=Js,HI%JecB ; do I=Is,HI%IecB
      rescaled_array(I,J) = scale*array(I,J)
    enddo ; enddo
    call subStats(HI, rescaled_array, mesg, sym_stats)
    deallocate(rescaled_array)
  else
    call subStats(HI, array, mesg, sym_stats)
  endif ; endif

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
    if (is_root_pe()) call chk_sum_msg("B-point:",bc0,mesg)
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

    if (is_root_pe()) call chk_sum_msg("B-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  else
    bcS = subchk(array, HI, 0, -hshift, scaling)
    bcE = subchk(array, HI, hshift, 0, scaling)
    bcW = subchk(array, HI, -hshift, 0, scaling)
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) call chk_sum_msg_NSEW("B-point:",bc0,bcN,bcS,bcE,bcW,mesg)
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
    do J=HI%jsc+dj,HI%jec+dj; do I=HI%isc+di,HI%iec+di
      bc = bitcount(abs(scale*array(I,J)))
      subchk = subchk + bc
    enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg, sym_stats)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%JsdB:), intent(in) :: array !< The array to be checksummed
    character(len=*), intent(in) :: mesg      !< An identifying message
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.

    integer :: i, j, n, IsB, JsB
    real :: aMean, aMin, aMax

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
    if (is_root_pe()) call chk_sum_msg("B-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_B_2d

!> Checksums a pair of 2d velocity arrays staggered at C-grid locations
subroutine chksum_uv_2d(mesg, arrayU, arrayV, HI, haloshift, symmetric, omit_corners, scale)
  character(len=*),                  intent(in) :: mesg   !< Identifying messages
  type(hor_index_type),              intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: arrayU !< The u-component array to be checksummed
  real, dimension(HI%isd:,HI%JsdB:), intent(in) :: arrayV !< The v-component array to be checksummed
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                 optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                             !! symmetric computational domain.
  logical,                 optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A scaling factor for these arrays.

  if (present(haloshift)) then
    call chksum_u_2d(arrayU, 'u '//mesg, HI, haloshift, symmetric, omit_corners, scale)
    call chksum_v_2d(arrayV, 'v '//mesg, HI, haloshift, symmetric, omit_corners, scale)
  else
    call chksum_u_2d(arrayU, 'u '//mesg, HI, symmetric=symmetric)
    call chksum_v_2d(arrayV, 'v '//mesg, HI, symmetric=symmetric)
  endif

end subroutine chksum_uv_2d

!> Checksums a pair of 3d velocity arrays staggered at C-grid locations
subroutine chksum_uv_3d(mesg, arrayU, arrayV, HI, haloshift, symmetric, omit_corners, scale)
  character(len=*),                    intent(in) :: mesg   !< Identifying messages
  type(hor_index_type),                intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%jsd:,:), intent(in) :: arrayU !< The u-component array to be checksummed
  real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: arrayV !< The v-component array to be checksummed
  integer,                   optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                   optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                               !! symmetric computational domain.
  logical,                   optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                      optional, intent(in) :: scale     !< A scaling factor for these arrays.

  if (present(haloshift)) then
    call chksum_u_3d(arrayU, 'u '//mesg, HI, haloshift, symmetric, omit_corners, scale)
    call chksum_v_3d(arrayV, 'v '//mesg, HI, haloshift, symmetric, omit_corners, scale)
  else
    call chksum_u_3d(arrayU, 'u '//mesg, HI, symmetric=symmetric)
    call chksum_v_3d(arrayV, 'v '//mesg, HI, symmetric=symmetric)
  endif

end subroutine chksum_uv_3d

!> Checksums a 2d array staggered at C-grid u points.
subroutine chksum_u_2d(array, mesg, HI, haloshift, symmetric, omit_corners, scale)
  type(hor_index_type),           intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,               optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                           !! symmetric computational domain.
  logical,               optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A scaling factor for this array.

  real, allocatable, dimension(:,:) :: rescaled_array
  real :: scaling
  integer :: i, j, Is
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%jsc:HI%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif
  scaling = 1.0 ; if (present(scale)) scaling = scale

  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then ; if (present(scale)) then
    allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                             LBOUND(array,2):UBOUND(array,2)) )
    rescaled_array(:,:) = 0.0
    Is = HI%isc ; if (sym_stats) Is = HI%isc-1
    do j=HI%jsc,HI%jec ; do I=Is,HI%IecB
      rescaled_array(I,j) = scale*array(I,j)
    enddo ; enddo
    call subStats(HI, rescaled_array, mesg, sym_stats)
    deallocate(rescaled_array)
  else
    call subStats(HI, array, mesg, sym_stats)
  endif ; endif

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
    if (is_root_pe()) call chk_sum_msg("u-point:",bc0,mesg)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (hshift==0) then
    bcW = subchk(array, HI, -hshift-1, 0, scaling)
    if (is_root_pe()) call chk_sum_msg_W("u-point:",bc0,bcW,mesg)
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

    if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  else
    bcS = subchk(array, HI, 0, -hshift, scaling)
    bcE = subchk(array, HI, hshift, 0, scaling)
    if (sym) then
      bcW = subchk(array, HI, -hshift-1, 0, scaling)
    else
      bcW = subchk(array, HI, -hshift, 0, scaling)
    endif
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) call chk_sum_msg_NSEW("u-point:",bc0,bcN,bcS,bcE,bcW,mesg)
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
    do j=HI%jsc+dj,HI%jec+dj; do I=HI%isc+di,HI%iec+di
      bc = bitcount(abs(scale*array(I,j)))
      subchk = subchk + bc
    enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg, sym_stats)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: array !< The array to be checksummed
    character(len=*), intent(in) :: mesg      !< An identifying message
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.

    integer :: i, j, n, IsB
    real :: aMean, aMin, aMax

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
    if (is_root_pe()) call chk_sum_msg("u-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_u_2d

!> Checksums a 2d array staggered at C-grid v points.
subroutine chksum_v_2d(array, mesg, HI, haloshift, symmetric, omit_corners, scale)
  type(hor_index_type),           intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%JsdB:), intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,               optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                           !! symmetric computational domain.
  logical,               optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                  optional, intent(in) :: scale     !< A scaling factor for this array.

  real, allocatable, dimension(:,:) :: rescaled_array
  real :: scaling
  integer :: i, j, Js
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%JscB:HI%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif
  scaling = 1.0 ; if (present(scale)) scaling = scale

  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then ; if (present(scale)) then
    allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                             LBOUND(array,2):UBOUND(array,2)) )
    rescaled_array(:,:) = 0.0
    Js = HI%jsc ; if (sym_stats) Js = HI%jsc-1
    do J=Js,HI%JecB ; do i=HI%isc,HI%iec
      rescaled_array(i,J) = scale*array(i,J)
    enddo ; enddo
    call subStats(HI, rescaled_array, mesg, sym_stats)
    deallocate(rescaled_array)
  else
    call subStats(HI, array, mesg, sym_stats)
  endif ; endif

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
    if (is_root_pe()) call chk_sum_msg("v-point:",bc0,mesg)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (hshift==0) then
    bcS = subchk(array, HI, 0, -hshift-1, scaling)
    if (is_root_pe()) call chk_sum_msg_S("v-point:",bc0,bcS,mesg)
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

    if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  else
    if (sym) then
      bcS = subchk(array, HI, 0, -hshift-1, scaling)
    else
      bcS = subchk(array, HI, 0, -hshift, scaling)
    endif
    bcE = subchk(array, HI, hshift, 0, scaling)
    bcW = subchk(array, HI, -hshift, 0, scaling)
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) call chk_sum_msg_NSEW("v-point:",bc0,bcN,bcS,bcE,bcW,mesg)
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
    do J=HI%jsc+dj,HI%jec+dj; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(scale*array(i,J)))
      subchk = subchk + bc
    enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg, sym_stats)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%JsdB:), intent(in) :: array !< The array to be checksummed
    character(len=*), intent(in) :: mesg      !< An identifying message
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.

    integer :: i, j, n, JsB
    real :: aMean, aMin, aMax

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
    if (is_root_pe()) call chk_sum_msg("v-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_v_2d

!> Checksums a 3d array staggered at tracer points.
subroutine chksum_h_3d(array, mesg, HI, haloshift, omit_corners, scale)
  type(hor_index_type),             intent(in) :: HI !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:,:),  intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                 optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A scaling factor for this array.

  real, allocatable, dimension(:,:,:) :: rescaled_array
  real :: scaling
  integer :: i, j, k
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif
  scaling = 1.0 ; if (present(scale)) scaling = scale

  if (calculateStatistics) then ; if (present(scale)) then
    allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                             LBOUND(array,2):UBOUND(array,2), &
                             LBOUND(array,3):UBOUND(array,3)) )
    rescaled_array(:,:,:) = 0.0
    do k=1,size(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      rescaled_array(i,j,k) = scale*array(i,j,k)
    enddo ; enddo ; enddo

    call subStats(HI, rescaled_array, mesg)
    deallocate(rescaled_array)
  else
    call subStats(HI, array, mesg)
  endif ; endif

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
    if (is_root_pe()) call chk_sum_msg("h-point:",bc0,mesg)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (do_corners) then
    bcSW = subchk(array, HI, -hshift, -hshift, scaling)
    bcSE = subchk(array, HI, hshift, -hshift, scaling)
    bcNW = subchk(array, HI, -hshift, hshift, scaling)
    bcNE = subchk(array, HI, hshift, hshift, scaling)

    if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  else
    bcS = subchk(array, HI, 0, -hshift, scaling)
    bcE = subchk(array, HI, hshift, 0, scaling)
    bcW = subchk(array, HI, -hshift, 0, scaling)
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) call chk_sum_msg_NSEW("h-point:",bc0,bcN,bcS,bcE,bcW,mesg)
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
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: array !< The array to be checksummed
    character(len=*), intent(in) :: mesg  !< An identifying message

    integer :: i, j, k, n
    real :: aMean, aMin, aMax

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
    if (is_root_pe()) call chk_sum_msg("h-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_h_3d

!> Checksums a 3d array staggered at corner points.
subroutine chksum_B_3d(array, mesg, HI, haloshift, symmetric, omit_corners, scale)
  type(hor_index_type),              intent(in) :: HI !< A horizontal index type
  real, dimension(HI%IsdB:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                   intent(in) :: mesg  !< An identifying message
  integer,                  optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                  optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                              !! symmetric computational domain.
  logical,                  optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                     optional, intent(in) :: scale     !< A scaling factor for this array.

  real, allocatable, dimension(:,:,:) :: rescaled_array
  real :: scaling
  integer :: i, j, k, Is, Js
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%JscB:HI%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif
  scaling = 1.0 ; if (present(scale)) scaling = scale
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then ; if (present(scale)) then
    allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                             LBOUND(array,2):UBOUND(array,2), &
                             LBOUND(array,3):UBOUND(array,3)) )
    rescaled_array(:,:,:) = 0.0
    Is = HI%isc ; if (sym_stats) Is = HI%isc-1
    Js = HI%jsc ; if (sym_stats) Js = HI%jsc-1
    do k=1,size(array,3) ; do J=Js,HI%JecB ; do I=Is,HI%IecB
      rescaled_array(I,J,k) = scale*array(I,J,k)
    enddo ; enddo ; enddo
    call subStats(HI, rescaled_array, mesg, sym_stats)
    deallocate(rescaled_array)
  else
    call subStats(HI, array, mesg, sym_stats)
  endif ; endif

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
    if (is_root_pe()) call chk_sum_msg("B-point:",bc0,mesg)
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

    if (is_root_pe()) call chk_sum_msg("B-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
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

    if (is_root_pe()) call chk_sum_msg_NSEW("B-point:",bc0,bcN,bcS,bcE,bcW,mesg)
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
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg, sym_stats)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
    character(len=*), intent(in) :: mesg      !< An identifying message
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.

    integer :: i, j, k, n, IsB, JsB
    real :: aMean, aMin, aMax

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
    if (is_root_pe()) call chk_sum_msg("B-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_B_3d

!> Checksums a 3d array staggered at C-grid u points.
subroutine chksum_u_3d(array, mesg, HI, haloshift, symmetric, omit_corners, scale)
  type(hor_index_type),             intent(in) :: HI !< A horizontal index type
  real, dimension(HI%isdB:,HI%Jsd:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                 optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                             !! symmetric computational domain.
  logical,                 optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A scaling factor for this array.

  real, allocatable, dimension(:,:,:) :: rescaled_array
  real :: scaling
  integer :: i, j, k, Is
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%jsc:HI%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif
  scaling = 1.0 ; if (present(scale)) scaling = scale
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then ; if (present(scale)) then
    allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                             LBOUND(array,2):UBOUND(array,2), &
                             LBOUND(array,3):UBOUND(array,3)) )
    rescaled_array(:,:,:) = 0.0
    Is = HI%isc ; if (sym_stats) Is = HI%isc-1
    do k=1,size(array,3) ; do j=HI%jsc,HI%jec ; do I=Is,HI%IecB
      rescaled_array(I,j,k) = scale*array(I,j,k)
    enddo ; enddo ; enddo
    call subStats(HI, rescaled_array, mesg, sym_stats)
    deallocate(rescaled_array)
  else
    call subStats(HI, array, mesg, sym_stats)
  endif ; endif

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
    if (is_root_pe()) call chk_sum_msg("u-point:",bc0,mesg)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (hshift==0) then
    bcW = subchk(array, HI, -hshift-1, 0, scaling)
    if (is_root_pe()) call chk_sum_msg_W("u-point:",bc0,bcW,mesg)
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

    if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  else
    bcS = subchk(array, HI, 0, -hshift, scaling)
    bcE = subchk(array, HI, hshift, 0, scaling)
    if (sym) then
      bcW = subchk(array, HI, -hshift-1, 0, scaling)
    else
      bcW = subchk(array, HI, -hshift, 0, scaling)
    endif
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) call chk_sum_msg_NSEW("u-point:",bc0,bcN,bcS,bcE,bcW,mesg)
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
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg, sym_stats)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%IsdB:,HI%jsd:,:), intent(in) :: array !< The array to be checksummed
    character(len=*), intent(in) :: mesg      !< An identifying message
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.

    integer :: i, j, k, n, IsB
    real :: aMean, aMin, aMax

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
    if (is_root_pe()) call chk_sum_msg("u-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_u_3d

!> Return the bitcount of an arbitrarily sized 3d array
integer function chksum_general_3d( array, scale_factor, istart, iend, jstart, jend, kstart, kend ) &
                            result(subchk)
  real, dimension(:,:,:), intent(in) :: array    !< Array to be checksummed
  real,    optional, intent(in) :: scale_factor  !< Factor to scale array by before checksum
  integer, optional, intent(in) :: istart        !< Starting index in the i-direction
  integer, optional, intent(in) :: iend          !< Ending index in the i-direction
  integer, optional, intent(in) :: jstart        !< Starting index in the j-direction
  integer, optional, intent(in) :: jend          !< Ending index in the j-direction
  integer, optional, intent(in) :: kstart        !< Starting index in the k-direction
  integer, optional, intent(in) :: kend          !< Ending index in the k-direction
  integer :: i, j, k, bc, is, ie, js, je, ks, ke
  real :: scale

  ! By default do not scale
  scale = 1.
  if (present(scale_factor)) scale = scale_factor

  ! Set the loop indices based on full array
  is = LBOUND(array,1) ; ie = UBOUND(array,1)
  js = LBOUND(array,2) ; je = UBOUND(array,2)
  ks = LBOUND(array,3) ; ke = UBOUND(array,3)

  ! Override indices if subdomain requested
  if (present(istart)) is = istart ; if (present(iend)) ie = iend
  if (present(jstart)) js = jstart ; if (present(jend)) je = jend
  if (present(kstart)) ks = kstart ; if (present(kend)) ke = kend

  subchk = 0
  do k=ks,ke ; do j=js,je ; do i=is,ie
    bc = bitcount(abs(scale*array(i,j,k)))
    subchk = subchk + bc
  enddo ; enddo ; enddo
  call sum_across_PEs(subchk)
  subchk=mod(subchk,1000000000)
end function chksum_general_3d

!> Return the bitcount of an arbitrarily sized 2d array by promotion to a 3d array
integer function chksum_general_2d( array_2d, scale_factor, istart, iend, jstart, jend )
  real, dimension(:,:), intent(in) :: array_2d   !< Array to be checksummed
  real,    optional, intent(in) :: scale_factor  !< Factor to scale array by before checksum
  integer, optional, intent(in) :: istart        !< Starting index in the i-direction
  integer, optional, intent(in) :: iend          !< Ending index in the i-direction
  integer, optional, intent(in) :: jstart        !< Starting index in the j-direction
  integer, optional, intent(in) :: jend          !< Ending index in the j-direction
  integer :: is, ie, js, je
  real, dimension(:,:,:), allocatable :: array_3d !< Promotion from 2d to 3d array

  is = LBOUND(array_2d,1) ; ie = UBOUND(array_2d,1)
  js = LBOUND(array_2d,2) ; je = UBOUND(array_2d,2)
  allocate(array_3d(is:ie, js:je,1))
  array_3d(:,:,1) = array_2d(:,:)
  chksum_general_2d = chksum_general_3d( array_3d, scale_factor, istart, iend, jstart, jend )
  deallocate(array_3d)
end function chksum_general_2d

!> Return the bitcount of an arbitrarily sized 1d array by promotion to a 3d array
integer function chksum_general_1d( array_1d, scale_factor, istart, iend )
  real, dimension(:), intent(in) :: array_1d      !< Array to be checksummed
  real,    optional,  intent(in) :: scale_factor  !< Factor to scale array by before checksum
  integer, optional,  intent(in) :: istart        !< Starting index in the i-direction
  integer, optional,  intent(in) :: iend          !< Ending index in the i-direction
  integer :: is, ie
  real, dimension(:,:,:), allocatable :: array_3d !< Promotion from 2d to 3d array

  is = LBOUND(array_1d,1) ; ie = UBOUND(array_1d,1)
  allocate(array_3d(is:ie, 1,1))
  array_3d(:,1,1) = array_1d(:)
  chksum_general_1d = chksum_general_3d( array_3d, scale_factor, istart, iend )
  deallocate(array_3d)
end function chksum_general_1d

!> Checksums a 3d array staggered at C-grid v points.
subroutine chksum_v_3d(array, mesg, HI, haloshift, symmetric, omit_corners, scale)
  type(hor_index_type),             intent(in) :: HI !< A horizontal index type
  real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)
  logical,                 optional, intent(in) :: symmetric !< If true, do the checksums on the full
                                                             !! symmetric computational domain.
  logical,                 optional, intent(in) :: omit_corners !< If true, avoid checking diagonal shifts
  real,                    optional, intent(in) :: scale     !< A scaling factor for this array.

  real, allocatable, dimension(:,:,:) :: rescaled_array
  real :: scaling
  integer :: i, j, k, Js
  integer :: bc0, bcSW, bcSE, bcNW, bcNE, hshift
  integer :: bcN, bcS, bcE, bcW
  logical :: do_corners, sym, sym_stats

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%JscB:HI%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif
  scaling = 1.0 ; if (present(scale)) scaling = scale
  sym_stats = .false. ; if (present(symmetric)) sym_stats = symmetric
  if (present(haloshift)) then ; if (haloshift > 0) sym_stats = .true. ; endif

  if (calculateStatistics) then ; if (present(scale)) then
    allocate( rescaled_array(LBOUND(array,1):UBOUND(array,1), &
                             LBOUND(array,2):UBOUND(array,2), &
                             LBOUND(array,3):UBOUND(array,3)) )
    rescaled_array(:,:,:) = 0.0
    Js = HI%jsc ; if (sym_stats) Js = HI%jsc-1
    do k=1,size(array,3) ; do J=Js,HI%JecB ; do i=HI%isc,HI%iec
      rescaled_array(i,J,k) = scale*array(i,J,k)
    enddo ; enddo ; enddo
    call subStats(HI, rescaled_array, mesg, sym_stats)
    deallocate(rescaled_array)
  else
    call subStats(HI, array, mesg, sym_stats)
  endif ; endif

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
    if (is_root_pe()) call chk_sum_msg("v-point:",bc0,mesg)
    return
  endif

  do_corners = .true. ; if (present(omit_corners)) do_corners = .not.omit_corners

  if (hshift==0) then
    bcS = subchk(array, HI, 0, -hshift-1, scaling)
    if (is_root_pe()) call chk_sum_msg_S("v-point:",bc0,bcS,mesg)
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

    if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  else
    if (sym) then
      bcS = subchk(array, HI, 0, -hshift-1, scaling)
    else
      bcS = subchk(array, HI, 0, -hshift, scaling)
    endif
    bcE = subchk(array, HI, hshift, 0, scaling)
    bcW = subchk(array, HI, -hshift, 0, scaling)
    bcN = subchk(array, HI, 0, hshift, scaling)

    if (is_root_pe()) call chk_sum_msg_NSEW("v-point:",bc0,bcN,bcS,bcE,bcW,mesg)
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
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg, sym_stats)
    type(hor_index_type), intent(in) ::  HI     !< A horizontal index type
    real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
    character(len=*), intent(in) :: mesg      !< An identifying message
    logical,          intent(in) :: sym_stats !< If true, evaluate the statistics on the
                                              !! full symmetric computational domain.

    integer :: i, j, k, n, JsB
    real :: aMean, aMin, aMax

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
    if (is_root_pe()) call chk_sum_msg("v-point:",aMean,aMin,aMax,mesg)
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
  logical :: call_mpp

  n = 0
  do i = LBOUND(x,1), UBOUND(x,1)
    if (is_NaN_0d(x(i))) n = n + 1
  enddo
  call_mpp = .true.
  if (present(skip_mpp)) call_mpp = .not.skip_mpp

  if (call_mpp) call sum_across_PEs(n)
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

!> Write a message including the checksum of the non-shifted array
subroutine chk_sum_msg1(fmsg,bc0,mesg)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  if (is_root_pe()) write(0,'(A,1(A,I10,X),A)') fmsg," c=",bc0,trim(mesg)
end subroutine chk_sum_msg1

!> Write a message including checksums of non-shifted and diagonally shifted arrays
subroutine chk_sum_msg5(fmsg,bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcSW !< The bitcount for SW shifted array
  integer,          intent(in) :: bcSE !< The bitcount for SE shifted array
  integer,          intent(in) :: bcNW !< The bitcount for NW shifted array
  integer,          intent(in) :: bcNE !< The bitcount for NE shifted array
  if (is_root_pe()) write(0,'(A,5(A,I10,1X),A)') &
     fmsg," c=",bc0,"sw=",bcSW,"se=",bcSE,"nw=",bcNW,"ne=",bcNE,trim(mesg)
end subroutine chk_sum_msg5

!> Write a message including checksums of non-shifted and laterally shifted arrays
subroutine chk_sum_msg_NSEW(fmsg,bc0,bcN,bcS,bcE,bcW,mesg)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcN !< The bitcount for N shifted array
  integer,          intent(in) :: bcS !< The bitcount for S shifted array
  integer,          intent(in) :: bcE !< The bitcount for E shifted array
  integer,          intent(in) :: bcW !< The bitcount for W shifted array
  if (is_root_pe()) write(0,'(A,5(A,I10,1X),A)') &
     fmsg," c=",bc0,"N=",bcN,"S=",bcS,"E=",bcE,"W=",bcW,trim(mesg)
end subroutine chk_sum_msg_NSEW

!> Write a message including checksums of non-shifted and southward shifted arrays
subroutine chk_sum_msg_S(fmsg,bc0,bcS,mesg)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcS  !< The bitcount of the south-shifted array
  if (is_root_pe()) write(0,'(A,2(A,I10,1X),A)') &
     fmsg," c=",bc0,"S=",bcS,trim(mesg)
end subroutine chk_sum_msg_S

!> Write a message including checksums of non-shifted and westward shifted arrays
subroutine chk_sum_msg_W(fmsg,bc0,bcW,mesg)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcW  !< The bitcount of the west-shifted array
  if (is_root_pe()) write(0,'(A,2(A,I10,1X),A)') &
     fmsg," c=",bc0,"W=",bcW,trim(mesg)
end subroutine chk_sum_msg_W

!> Write a message including checksums of non-shifted and southwestward shifted arrays
subroutine chk_sum_msg2(fmsg,bc0,bcSW,mesg)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  integer,          intent(in) :: bc0  !< The bitcount of the non-shifted array
  integer,          intent(in) :: bcSW !< The bitcount of the southwest-shifted array
  if (is_root_pe()) write(0,'(A,2(A,I9,1X),A)') &
     fmsg," c=",bc0,"s/w=",bcSW,trim(mesg)
end subroutine chk_sum_msg2

!> Write a message including the global mean, maximum and minimum of an array
subroutine chk_sum_msg3(fmsg,aMean,aMin,aMax,mesg)
  character(len=*), intent(in) :: fmsg !< A checksum code-location specific preamble
  character(len=*), intent(in) :: mesg !< An identifying message supplied by top-level caller
  real,             intent(in) :: aMean !< The mean value of the array
  real,             intent(in) :: aMin !< The minimum value of the array
  real,             intent(in) :: aMax !< The maximum value of the array
  if (is_root_pe()) write(0,'(A,3(A,ES25.16,1X),A)') &
     fmsg," mean=",aMean,"min=",aMin,"max=",aMax,trim(mesg)
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
integer function bitcount( x )
  real :: x !< Number to be bitcount

  ! Local variables
  integer(kind(x)) :: y !< Store the integer representation of the memory used by x
  integer :: bit

  bitcount = 0
  y = transfer(x,y)

  ! Fortran standard says that bit indexing start at 0
  do bit = 0, bit_size(y)-1
    if (BTEST(y,bit)) bitcount = bitcount+1
  enddo

end function bitcount

end module MOM_checksums
