module MOM_checksums

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

use MOM_coms, only : PE_here, root_PE, num_PEs, sum_across_PEs
use MOM_coms, only : min_across_PEs, max_across_PEs
use MOM_error_handler, only : MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type

implicit none ; private

public :: hchksum, qchksum, uchksum, vchksum, chksum, is_NaN
public :: totalStuff, totalTandS
public :: MOM_checksums_init

interface hchksum
  module procedure chksum_h_2d
  module procedure chksum_h_3d
  module procedure chksum_h_2d_G
  module procedure chksum_h_3d_G
end interface

interface qchksum
  module procedure chksum_q_2d
  module procedure chksum_q_3d
  module procedure chksum_q_2d_G
  module procedure chksum_q_3d_G
end interface

interface uchksum
  module procedure chksum_u_2d
  module procedure chksum_u_3d
  module procedure chksum_u_2d_G
  module procedure chksum_u_3d_G
end interface

interface vchksum
  module procedure chksum_v_2d
  module procedure chksum_v_3d
  module procedure chksum_v_2d_G
  module procedure chksum_v_3d_G
end interface

interface chksum
  module procedure chksum1d
!  module procedure chksum2d
!  module procedure chksum3d
end interface

interface chk_sum_msg
  module procedure chk_sum_msg1
  module procedure chk_sum_msg3
  module procedure chk_sum_msg5
end interface

interface is_NaN
  module procedure is_NaN_0d
  module procedure is_NaN_1d
  module procedure is_NaN_2d
  module procedure is_NaN_3d
end interface

integer, parameter :: default_shift=0
logical :: calculateStatistics=.false. ! If true, report min, max and mean
                               ! instead of the bitcount checksum
logical :: checkForNaNs=.true. ! If true, checks array for NaNs and cause
                               ! FATAL error is any are found

contains

! =====================================================================

!> chksum_h_2d_G performs checksums on a 2d array staggered at tracer points.
subroutine chksum_h_2d_G(array, mesg, G, haloshift)
  type(ocean_grid_type),           intent(in) :: G     !< The ocean grid type
  real, dimension(G%isd:,G%jsd:),  intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  call chksum_h_2d(array, mesg, G%HI, haloshift)
end subroutine chksum_h_2d_G

!> chksum_q_2d_G performs checksums on a 2d array staggered at corner points.
subroutine chksum_q_2d_G(array, mesg, G, haloshift)
  type(ocean_grid_type),            intent(in) :: G     !< The ocean grid type
  real, dimension(G%IsdB:,G%JsdB:), intent(in) :: array !< The array to be checksummed
  character(len=*),                 intent(in) :: mesg  !< An identifying message
  integer,                optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  call chksum_q_2d(array, mesg, G%HI, haloshift)
end subroutine chksum_q_2d_G

!> chksum_u_2d_G performs checksums on a 2d array staggered at C-grid u points.
subroutine chksum_u_2d_G(array, mesg, G, haloshift)
  type(ocean_grid_type),           intent(in) :: G     !< The ocean grid type
  real, dimension(G%IsdB:,G%jsd:), intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  call chksum_u_2d(array, mesg, G%HI, haloshift)
end subroutine chksum_u_2d_G

!> chksum_v_2d_G performs checksums on a 2d array staggered at C-grid v points.
subroutine chksum_v_2d_G(array, mesg, G, haloshift)
  type(ocean_grid_type),           intent(in) :: G     !< The ocean grid type
  real, dimension(G%isd:,G%JsdB:), intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  call chksum_v_2d(array, mesg, G%HI, haloshift)
end subroutine chksum_v_2d_G

!> chksum_h_3d performs checksums on a 3d array staggered at tracer points.
subroutine chksum_h_3d_G(array, mesg, G, haloshift)
  type(ocean_grid_type),             intent(in) :: G !< The ocean grid type
  real, dimension(G%isd:,G%jsd:,:),  intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  call chksum_h_3d(array, mesg, G%HI, haloshift)
end subroutine chksum_h_3d_G

!> chksum_q_3d_G performs checksums on a 3d array staggered at corner points.
subroutine chksum_q_3d_G(array, mesg, G, haloshift)
  type(ocean_grid_type),              intent(in) :: G !< The ocean grid type
  real, dimension(G%IsdB:,G%JsdB:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                   intent(in) :: mesg  !< An identifying message
  integer,                  optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  call chksum_q_3d(array, mesg, G%HI, haloshift)
end subroutine chksum_q_3d_G

!> chksum_u_3d_G performs checksums on a 3d array staggered at C-grid u points.
subroutine chksum_u_3d_G(array, mesg, G, haloshift)
  type(ocean_grid_type),             intent(in) :: G !< The ocean grid type
  real, dimension(G%isdB:,G%Jsd:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  call chksum_u_3d(array, mesg, G%HI, haloshift)
end subroutine chksum_u_3d_G


!> chksum_v_3d_G performs checksums on a 3d array staggered at C-grid v points.
subroutine chksum_v_3d_G(array, mesg, G, haloshift)
  type(ocean_grid_type),             intent(in) :: G !< The ocean grid type
  real, dimension(G%isd:,G%JsdB:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  call chksum_v_3d(array, mesg, G%HI, haloshift)
end subroutine chksum_v_3d_G

! The following versions work with a horizontal index type instead of the
! ocean grid type.

!> chksum_h_2d performs checksums on a 2d array staggered at tracer points.
subroutine chksum_h_2d(array, mesg, HI, haloshift)
  type(hor_index_type),           intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:),  intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%jsc:HI%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(HI, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. &
       HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. &
       HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_h_2d: haloshift =',hshift
    write(0,*) 'chksum_h_2d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_h_2d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_h_2d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("h-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%jsd:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=HI%jsc+dj,HI%jec+dj; do i=HI%isc+di,HI%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%jsd:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(HI%isc,HI%jsc)
    aMax = array(HI%isc,HI%jsc)
    n = 0
    do j=HI%jsc,HI%jec; do i=HI%isc,HI%iec
        aMean = aMean + array(i,j)
        aMin = min(aMin, array(i,j))
        aMax = max(aMax, array(i,j))
        n = n + 1
    enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("h-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_h_2d

! =====================================================================

!> chksum_q_2d performs checksums on a 2d array staggered at corner points.
subroutine chksum_q_2d(array, mesg, HI, haloshift)
  type(hor_index_type),            intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%JsdB:), intent(in) :: array !< The array to be checksummed
  character(len=*),                 intent(in) :: mesg  !< An identifying message
  integer,                optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%JscB:HI%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(HI, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. &
       HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. &
       HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_q_2d: haloshift =',hshift
    write(0,*) 'chksum_q_2d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_q_2d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_q_2d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("q-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("q-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%JsdB:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=HI%jsc+dj,HI%jec+dj; do i=HI%isc+di,HI%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%JsdB:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(HI%isc,HI%jsc)
    aMax = array(HI%isc,HI%jsc)
    n = 0
    do j=HI%jsc,HI%jec; do i=HI%isc,HI%iec
        aMean = aMean + array(i,j)
        aMin = min(aMin, array(i,j))
        aMax = max(aMax, array(i,j))
        n = n + 1
    enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("q-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_q_2d

! =====================================================================

!> chksum_u_2d performs checksums on a 2d array staggered at C-grid u points.
subroutine chksum_u_2d(array, mesg, HI, haloshift)
  type(hor_index_type),           intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%jsc:HI%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(HI, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. &
       HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. &
       HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_u_2d: haloshift =',hshift
    write(0,*) 'chksum_u_2d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_u_2d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_u_2d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("u-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=HI%jsc+dj,HI%jec+dj; do i=HI%isc+di,HI%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%jsd:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(HI%isc,HI%jsc)
    aMax = array(HI%isc,HI%jsc)
    n = 0
    do j=HI%jsc,HI%jec; do i=HI%isc,HI%iec
        aMean = aMean + array(i,j)
        aMin = min(aMin, array(i,j))
        aMax = max(aMax, array(i,j))
        n = n + 1
    enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("u-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_u_2d

! =====================================================================

!> chksum_v_2d performs checksums on a 2d array staggered at C-grid v points.
subroutine chksum_v_2d(array, mesg, HI, haloshift)
  type(hor_index_type),           intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%JsdB:), intent(in) :: array !< The array to be checksummed
  character(len=*),                intent(in) :: mesg  !< An identifying message
  integer,               optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%JscB:HI%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(HI, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. &
       HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. &
       HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_v_2d: haloshift =',hshift
    write(0,*) 'chksum_v_2d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_v_2d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_v_2d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("v-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%JsdB:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=HI%jsc+dj,HI%jec+dj; do i=HI%isc+di,HI%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%JsdB:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(HI%isc,HI%jsc)
    aMax = array(HI%isc,HI%jsc)
    n = 0
    do j=HI%jsc,HI%jec; do i=HI%isc,HI%iec
        aMean = aMean + array(i,j)
        aMin = min(aMin, array(i,j))
        aMax = max(aMax, array(i,j))
        n = n + 1
    enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("v-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_v_2d

! =====================================================================

!> chksum_h_3d performs checksums on a 3d array staggered at tracer points.
subroutine chksum_h_3d(array, mesg, HI, haloshift)
  type(hor_index_type),             intent(in) :: HI !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:,:),  intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%jsc:HI%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(HI, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. &
       HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. &
       HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_h_3d: haloshift =',hshift
    write(0,*) 'chksum_h_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_h_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_h_3d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("h-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(array(i,j,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(HI%isc,HI%jsc,1)
    aMax = array(HI%isc,HI%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      aMean = aMean + array(i,j,k)
      aMin = min(aMin, array(i,j,k))
      aMax = max(aMax, array(i,j,k))
      n = n + 1
    enddo ; enddo ; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("h-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_h_3d

! =====================================================================

!> chksum_q_3d performs checksums on a 3d array staggered at corner points.
subroutine chksum_q_3d(array, mesg, HI, haloshift)
  type(hor_index_type),              intent(in) :: HI !< A horizontal index type
  real, dimension(HI%IsdB:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                   intent(in) :: mesg  !< An identifying message
  integer,                  optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%JscB:HI%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(HI, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. &
       HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. &
       HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_q_3d: haloshift =',hshift
    write(0,*) 'chksum_q_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_q_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_q_3d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
    if (is_root_pe()) call chk_sum_msg("q-point:",bc0,mesg)
    return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("q-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%JsdB:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(array(i,j,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%JsdB:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(HI%isc,HI%jsc,1)
    aMax = array(HI%isc,HI%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      aMean = aMean + array(i,j,k)
      aMin = min(aMin, array(i,j,k))
      aMax = max(aMax, array(i,j,k))
      n = n + 1
    enddo ; enddo ; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("q-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_q_3d

! =====================================================================

!> chksum_u_3d performs checksums on a 3d array staggered at C-grid u points.
subroutine chksum_u_3d(array, mesg, HI, haloshift)
  type(hor_index_type),             intent(in) :: HI !< A horizontal index type
  real, dimension(HI%isdB:,HI%Jsd:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(HI%IscB:HI%IecB,HI%jsc:HI%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(HI, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. &
       HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. &
       HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_u_3d: haloshift =',hshift
    write(0,*) 'chksum_u_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_u_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_u_3d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
    if (is_root_pe()) call chk_sum_msg("u-point:",bc0,mesg)
    return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%jsd:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(array(i,j,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%IsdB:,HI%jsd:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(HI%isc,HI%jsc,1)
    aMax = array(HI%isc,HI%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      aMean = aMean + array(i,j,k)
      aMin = min(aMin, array(i,j,k))
      aMax = max(aMax, array(i,j,k))
      n = n + 1
    enddo ; enddo ; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("u-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_u_3d

! =====================================================================

!> chksum_v_3d performs checksums on a 3d array staggered at C-grid v points.
subroutine chksum_v_3d(array, mesg, HI, haloshift)
  type(hor_index_type),             intent(in) :: HI !< A horizontal index type
  real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: array !< The array to be checksummed
  character(len=*),                  intent(in) :: mesg  !< An identifying message
  integer,                 optional, intent(in) :: haloshift !< The width of halos to check (default 0)

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(HI%isc:HI%iec,HI%JscB:HI%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(HI, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=HI%ied-HI%iec

  if ( HI%isc-hshift<HI%isd .or. &
       HI%iec+hshift>HI%ied .or. &
       HI%jsc-hshift<HI%jsd .or. &
       HI%jec+hshift>HI%jed ) then
    write(0,*) 'chksum_v_3d: haloshift =',hshift
    write(0,*) 'chksum_v_3d: isd,isc,iec,ied=',HI%isd,HI%isc,HI%iec,HI%ied
    write(0,*) 'chksum_v_3d: jsd,jsc,jec,jed=',HI%jsd,HI%jsc,HI%jec,HI%jed
    call chksum_error(FATAL,'Error in chksum_v_3d '//trim(mesg))
  endif

  bc0=subchk(array, HI, 0, 0)

  if (hshift==0) then
    if (is_root_pe()) call chk_sum_msg("v-point:",bc0,mesg)
    return
  endif

  bcSW=subchk(array, HI, -hshift, -hshift)
  bcSE=subchk(array, HI, hshift, -hshift)
  bcNW=subchk(array, HI, -hshift, hshift)
  bcNE=subchk(array, HI, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, HI, di, dj)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc+dj,HI%jec+dj ; do i=HI%isc+di,HI%iec+di
      bc = bitcount(abs(array(i,j,k)))
      subchk = subchk + bc
    enddo ; enddo ; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(HI, array, mesg)
    type(hor_index_type), intent(in) :: HI
    real, dimension(HI%isd:,HI%JsdB:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(HI%isc,HI%jsc,1)
    aMax = array(HI%isc,HI%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3) ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      aMean = aMean + array(i,j,k)
      aMin = min(aMin, array(i,j,k))
      aMax = max(aMax, array(i,j,k))
      n = n + 1
    enddo ; enddo ; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("v-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_v_3d


! =====================================================================

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
  integer :: bitcount
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
    write(0,'(A40,1X,ES22.13,1X,I12)') mesg, sum, sum_bc

end subroutine chksum1d

! =====================================================================
!   These are the older version of chksum that do not take the grid staggering
! into account.

subroutine chksum2d(array, mesg, start_x, end_x, start_y, end_y)

  real, dimension(:,:) :: array
  character(len=*) :: mesg
  integer, optional :: start_x, end_x, start_y, end_y

  integer :: bitcount
  integer :: xs,xe,ys,ye,i,j,sum1,bc
  real :: sum
  real, allocatable :: sum_here(:)
  integer :: pe_num   ! pe number of the data

  xs = LBOUND(array,1) ; xe = UBOUND(array,1)
  ys = LBOUND(array,2) ; ye = UBOUND(array,2)
  if (present(start_x)) xs = start_x
  if (present(end_x  )) xe = end_x
  if (present(start_y)) ys = start_y
  if (present(end_y  )) ye = end_y

  sum = 0.0 ; sum1 = 0
  do i=xs,xe ; do j=ys,ye
    sum = sum + array(i,j)
    bc = bitcount(abs(array(i,j)))
    sum1 = sum1 + bc
  enddo ; enddo

  pe_num = pe_here() + 1 - root_pe()
  allocate(sum_here(num_pes())) ; sum_here(:) = 0.0
  sum_here(pe_num) = sum
  call sum_across_PEs(sum_here,num_pes())
  sum = 0.0
  do i=1,num_pes() ; sum = sum + sum_here(i) ; enddo
!   call sum_across_PEs(sum)
  call sum_across_PEs(sum1)
  deallocate(sum_here)

  if (is_root_pe()) &
    write(0,'(A40,1X,ES22.13,1X,I12)') &
      mesg, sum, sum1
!    write(0,'(A40,1X,Z16.16,1X,Z16.16,1X,ES25.16,1X,I12)') &
!      mesg, sum, sum1, sum, sum1

end subroutine chksum2d

subroutine chksum3d(array, mesg, start_x, end_x, start_y, end_y, start_z, end_z)

  real, dimension(:,:,:) :: array
  character(len=*) :: mesg
  integer, optional :: start_x, end_x, start_y, end_y, start_z, end_z

  integer :: bitcount
  integer :: xs,xe,ys,ye,zs,ze,i,j,k, bc,sum1
  real :: sum
  real, allocatable :: sum_here(:)
  integer :: pe_num   ! pe number of the data

  xs = LBOUND(array,1) ; xe = UBOUND(array,1)
  ys = LBOUND(array,2) ; ye = UBOUND(array,2)
  zs = LBOUND(array,3) ; ze = UBOUND(array,3)
  if (present(start_x)) xs = start_x
  if (present(end_x  )) xe = end_x
  if (present(start_y)) ys = start_y
  if (present(end_y  )) ye = end_y
  if (present(start_z)) zs = start_z
  if (present(end_z  )) ze = end_z

  sum = 0.0 ; sum1 = 0
  do i=xs,xe ; do j=ys,ye ; do k=zs,ze
    sum = sum + array(i,j,k)
    bc = bitcount(ABS(array(i,j,k)))
    sum1 = sum1 + bc
  enddo ; enddo ; enddo

  pe_num = pe_here() + 1 - root_pe()
  allocate(sum_here(num_pes())) ; sum_here(:) = 0.0
  sum_here(pe_num) = sum
  call sum_across_PEs(sum_here,num_pes())
  sum = 0.0
  do i=1,num_pes() ; sum = sum + sum_here(i) ; enddo
!   call sum_across_PEs(sum)
  call sum_across_PEs(sum1)
  deallocate(sum_here)

  if (is_root_pe()) &
    write(0,'(A40,1X,ES22.13,1X,I12)') &
      mesg, sum, sum1
!    write(0,'(A40,1X,Z16.16,1X,Z16.16,1X,ES25.16,1X,I12)') &
!      mesg, sum, sum1, sum, sum1

end subroutine chksum3d

! =====================================================================

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

! =====================================================================

!> This function returns .true. if any element of x is a NaN, and .false. otherwise.
function is_NaN_1d(x, skip_mpp)
  real, dimension(:), intent(in) :: x !< The array to be checked for NaNs.
  logical :: is_NaN_1d
  logical, optional :: skip_mpp  !< If true, only check this array only on the local PE (default false).

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

! =====================================================================

!> This function returns .true. if any element of x is a NaN, and .false. otherwise.
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

! =====================================================================

!> This function returns .true. if any element of x is a NaN, and .false. otherwise.
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

! =====================================================================

!> This function returns the sum over computational domain of all
!! processors of hThick*stuff, where stuff is a 3-d array at tracer points.
function totalStuff(G, hThick, stuff)
  type(ocean_grid_type),            intent(in) :: G      !< The ocean grid type
  real, dimension(G%isd:,G%jsd:,:), intent(in) :: hThick !< The array of thicknesses to use as weights
  real, dimension(G%isd:,G%jsd:,:), intent(in) :: stuff  !< The array of stuff to be summed 
  real                                         :: totalStuff
  integer :: i, j, k

  totalStuff = 0.
  do k = 1, G%ke ; do j = G%jsc, G%jec ; do i = G%isc, G%iec
    totalStuff = totalStuff + hThick(i,j,k) * stuff(i,j,k) * G%areaT(i,j)
  enddo ; enddo ; enddo
  call sum_across_PEs(totalStuff)

end function totalStuff

! =====================================================================

!> This subroutine display the total thickness, temperature and salinity
!! as well as the change since the last call.
!! NOTE: This subroutine uses "save" data which is not thread safe and is purely
!! for extreme debugging without a proper debugger.
subroutine totalTandS(G, hThick, temperature, salinity, mesg)
  type(ocean_grid_type),            intent(in) :: G      !< The ocean grid type
  real, dimension(G%isd:,G%jsd:,:), intent(in) :: hThick !< The array of thicknesses to use as weights
  real, dimension(G%isd:,G%jsd:,:), intent(in) :: temperature !< The temperature field to sum
  real, dimension(G%isd:,G%jsd:,:), intent(in) :: salinity    !< The salinity field to sum
  character(len=*),                 intent(in) :: mesg        !< An identifying message

  ! NOTE: This subroutine uses "save" data which is not thread safe and is purely for
  ! extreme debugging without a proper debugger.
  real, save :: totalH = 0., totalT = 0., totalS = 0.

  logical, save :: firstCall = .true.
  real :: thisH, thisT, thisS, delH, delT, delS
  integer :: i, j, k

  thisH = 0.
  do k = 1, G%ke ; do j = G%jsc, G%jec ; do i = G%isc, G%iec
    thisH = thisH + hThick(i,j,k) * G%areaT(i,j)
  enddo ; enddo ; enddo
  call sum_across_PEs(thisH)
  thisT = totalStuff(G, hThick, temperature)
  thisS = totalStuff(G, hThick, salinity)

  if (is_root_pe()) then
    if (firstCall) then
      totalH = thisH ; totalT = thisT ; totalS = thisS
      write(0,*) 'Totals H,T,S:',thisH,thisT,thisS,' ',mesg
      firstCall = .false.
    else
      delH = thisH - totalH
      delT = thisT - totalT
      delS = thisS - totalS
      totalH = thisH ; totalT = thisT ; totalS = thisS
      write(0,*) 'Tot/del H,T,S:',thisH,thisT,thisS,delH,delT,delS,' ',mesg
    endif
  endif

end subroutine totalTandS

! =====================================================================

subroutine chk_sum_msg1(fmsg,bc0,mesg)
  character(len=*), intent(in) :: fmsg, mesg
  integer,          intent(in) :: bc0
  if (is_root_pe()) write(0,'(A,1(A,I10,X),A)') fmsg," c=",bc0,mesg
end subroutine chk_sum_msg1

! =====================================================================

subroutine chk_sum_msg5(fmsg,bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  character(len=*), intent(in) :: fmsg, mesg
  integer,          intent(in) :: bc0,bcSW,bcSE,bcNW,bcNE
  if (is_root_pe()) write(0,'(A,5(A,I10,1X),A)') &
     fmsg," c=",bc0,"sw=",bcSW,"se=",bcSE,"nw=",bcNW,"ne=",bcNE,mesg
end subroutine chk_sum_msg5

! =====================================================================

subroutine chk_sum_msg3(fmsg,aMean,aMin,aMax,mesg)
  character(len=*), intent(in) :: fmsg, mesg
  real,             intent(in) :: aMean,aMin,aMax
  if (is_root_pe()) write(0,'(A,3(A,ES12.4,1X),A)') &
     fmsg," mean=",aMean,"min=",aMin,"max=",aMax,mesg
end subroutine chk_sum_msg3

! =====================================================================

!> MOM_checksums_init initializes the MOM_checksums module. As it happens, the
!! only thing that it does is to log the version of this module.
subroutine MOM_checksums_init(param_file)
  type(param_file_type),   intent(in)    :: param_file
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_checksums" ! This module's name.

  call log_version(param_file, mod, version)

end subroutine MOM_checksums_init

! =====================================================================

subroutine chksum_error(signal, message)
  ! Wrapper for MOM_error to help place specific break points in
  ! debuggers
  integer, intent(in) :: signal
  character(len=*), intent(in) :: message
  call MOM_error(signal, message)
end subroutine chksum_error

! =====================================================================

end module MOM_checksums
