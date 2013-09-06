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

implicit none ; private

public :: hchksum, qchksum, uchksum, vchksum, chksum, is_NaN
public :: totalStuff, totalTandS
public :: MOM_checksums_init

interface hchksum
! module procedure hchksum2d
! module procedure hchksum3d
  module procedure chksum_h_2d
  module procedure chksum_h_3d
end interface

interface qchksum
! module procedure hchksum2d
! module procedure hchksum3d
  module procedure chksum_q_2d
  module procedure chksum_q_3d
end interface

interface uchksum
! module procedure hchksum2d
! module procedure hchksum3d
  module procedure chksum_u_2d
  module procedure chksum_u_3d
end interface

interface vchksum
! module procedure hchksum2d
! module procedure hchksum3d
  module procedure chksum_v_2d
  module procedure chksum_v_3d
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

subroutine chksum_h_2d(array, mesg, G, haloshift)
  type(ocean_grid_type), intent(in) :: G
  real, dimension(G%isd:,G%jsd:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%isc:G%iec,G%jsc:G%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_h_2d: haloshift =',hshift
    write(0,*) 'chksum_h_2d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_h_2d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_h_2d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("h-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%jsd:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%jsd:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc)
    aMax = array(G%isc,G%jsc)
    n = 0
    do j=G%jsc,G%jec; do i=G%isc,G%iec
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

subroutine chksum_q_2d(array, mesg, G, haloshift)
  type(ocean_grid_type), intent(in) :: G
  real, dimension(G%IsdB:,G%JsdB:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%IscB:G%IecB,G%JscB:G%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_q_2d: haloshift =',hshift
    write(0,*) 'chksum_q_2d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_q_2d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_q_2d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("q-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("q-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%JsdB:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%JsdB:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc)
    aMax = array(G%isc,G%jsc)
    n = 0
    do j=G%jsc,G%jec; do i=G%isc,G%iec
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

subroutine chksum_u_2d(array, mesg, G, haloshift)
  type(ocean_grid_type), intent(in) :: G
  real, dimension(G%IsdB:,G%jsd:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%IscB:G%IecB,G%jsc:G%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_u_2d: haloshift =',hshift
    write(0,*) 'chksum_u_2d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_u_2d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_u_2d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("u-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%jsd:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%jsd:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc)
    aMax = array(G%isc,G%jsc)
    n = 0
    do j=G%jsc,G%jec; do i=G%isc,G%iec
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

subroutine chksum_v_2d(array, mesg, G, haloshift)
  type(ocean_grid_type), intent(in) :: G
  real, dimension(G%isd:,G%JsdB:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%isc:G%iec,G%JscB:G%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_v_2d: haloshift =',hshift
    write(0,*) 'chksum_v_2d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_v_2d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_v_2d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("v-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%JsdB:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%JsdB:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc)
    aMax = array(G%isc,G%jsc)
    n = 0
    do j=G%jsc,G%jec; do i=G%isc,G%iec
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

subroutine chksum_h_3d(array, mesg, G, haloshift)
  type(ocean_grid_type), intent(in) :: G
  real, dimension(G%isd:,G%jsd:,:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%isc:G%iec,G%jsc:G%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_h_3d: haloshift =',hshift
    write(0,*) 'chksum_h_3d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_h_3d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_h_3d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("h-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%jsd:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j,k)))
        subchk = subchk + bc
    enddo; enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%jsd:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc,1)
    aMax = array(G%isc,G%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc,G%jec; do i=G%isc,G%iec
        aMean = aMean + array(i,j,k)
        aMin = min(aMin, array(i,j,k))
        aMax = max(aMax, array(i,j,k))
        n = n + 1
    enddo; enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("h-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_h_3d

! =====================================================================

subroutine chksum_q_3d(array, mesg, G, haloshift)
  type(ocean_grid_type), intent(in) :: G
  real, dimension(G%IsdB:,G%JsdB:,:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%IscB:G%IecB,G%JscB:G%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_q_3d: haloshift =',hshift
    write(0,*) 'chksum_q_3d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_q_3d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_q_3d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("q-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("q-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%JsdB:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j,k)))
        subchk = subchk + bc
    enddo; enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%JsdB:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc,1)
    aMax = array(G%isc,G%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc,G%jec; do i=G%isc,G%iec
        aMean = aMean + array(i,j,k)
        aMin = min(aMin, array(i,j,k))
        aMax = max(aMax, array(i,j,k))
        n = n + 1
    enddo; enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("q-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_q_3d

! =====================================================================

subroutine chksum_u_3d(array, mesg, G, haloshift)
  type(ocean_grid_type), intent(in) :: G
  real, dimension(G%IsdB:,G%jsd:,:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%IscB:G%IecB,G%jsc:G%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_u_3d: haloshift =',hshift
    write(0,*) 'chksum_u_3d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_u_3d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_u_3d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("u-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%jsd:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j,k)))
        subchk = subchk + bc
    enddo; enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%jsd:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc,1)
    aMax = array(G%isc,G%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc,G%jec; do i=G%isc,G%iec
        aMean = aMean + array(i,j,k)
        aMin = min(aMin, array(i,j,k))
        aMax = max(aMax, array(i,j,k))
        n = n + 1
    enddo; enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("u-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_u_3d

! =====================================================================

subroutine chksum_v_3d(array, mesg, G, haloshift)
  type(ocean_grid_type), intent(in) :: G
  real, dimension(G%isd:,G%JsdB:,:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%isc:G%iec,G%JscB:G%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_v_3d: haloshift =',hshift
    write(0,*) 'chksum_v_3d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_v_3d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_v_3d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("v-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%JsdB:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j,k)))
        subchk = subchk + bc
    enddo; enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(ocean_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%JsdB:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc,1)
    aMax = array(G%isc,G%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc,G%jec; do i=G%isc,G%iec
        aMean = aMean + array(i,j,k)
        aMin = min(aMin, array(i,j,k))
        aMax = max(aMax, array(i,j,k))
        n = n + 1
    enddo; enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("v-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_v_3d

! =====================================================================

subroutine hchksum2d(array, mesg, start_x, end_x, start_y, end_y, haloshift)

  real, dimension(:,:) :: array
  character(len=*) :: mesg
  integer :: start_x, end_x, start_y, end_y
  integer, optional :: haloshift

  integer :: xs,xe,ys,ye,bc0,bcSW,bcSE,bcNW,bcNE,hshift

  xs = LBOUND(array,1) ; xe = UBOUND(array,1)
  ys = LBOUND(array,2) ; ye = UBOUND(array,2)
  if ( xs>start_x-1 .or. xe<end_x+1 .or. ys>start_y-1 .or. ye<end_y+1) then
    write(0,*) 'hchksum2d: must pass full array with haloes!'
    write(0,*) 'hchksum2d: xs,xe,ys,ye=',xs,xe,ys,ye
    write(0,*) 'hchksum2d: start_x,end_x,start_y,end_y=',start_x,end_x,start_y,end_y
    call chksum_error(FATAL,'Error in hchksum2d '//trim(mesg))
  endif

  bc0=hsubsum2d(array, start_x, end_x, start_y, end_y, 0, 0)
  if (present(haloshift)) then
    if (haloshift==0) then
      if (is_root_pe()) write(0,'(A40,5(1X,A,I12))') mesg,"c=",bc0
      return
    endif
    hshift=haloshift
  else
    hshift=1
  endif

  bcSW=hsubsum2d(array, start_x, end_x, start_y, end_y, -hshift, -hshift)
  bcSE=hsubsum2d(array, start_x, end_x, start_y, end_y, hshift, -hshift)
  bcNW=hsubsum2d(array, start_x, end_x, start_y, end_y, -hshift, hshift)
  bcNE=hsubsum2d(array, start_x, end_x, start_y, end_y, hshift, hshift)

  if (is_root_pe()) &
  write(0,'(A40,5(1X,A,I12))') mesg,"c=",bc0,"sw=",bcSW,"se=",bcSE,"nw=",bcNW,"ne=",bcNE

  contains

  integer function hsubsum2d(array, start_x, end_x, start_y, end_y, di, dj)
  real, dimension(:,:) :: array
  integer :: start_x, end_x, start_y, end_y, di, dj
  integer :: bitcount
  integer :: i, j, bc
  hsubsum2d = 0
  do j=start_y+dj,end_y+dj
    do i=start_x+di,end_x+di
      bc = bitcount(abs(array(i,j)))
      hsubsum2d = hsubsum2d + bc
    enddo
  enddo
  call sum_across_PEs(hsubsum2d)
  end function hsubsum2d
end subroutine hchksum2d

! =====================================================================

subroutine hchksum3d(array, mesg, start_x, end_x, start_y, end_y, haloshift)

  real, dimension(:,:,:) :: array
  character(len=*) :: mesg
  integer :: start_x, end_x, start_y, end_y
  integer, optional :: haloshift

  integer :: xs,xe,ys,ye,bc0,bcSW,bcSE,bcNW,bcNE,hshift

  xs = LBOUND(array,1) ; xe = UBOUND(array,1)
  ys = LBOUND(array,2) ; ye = UBOUND(array,2)
  if ( xs>start_x-1 .or. xe<end_x+1 .or. ys>start_y-1 .or. ye<end_y+1) then
    write(0,*) 'hchksum2d: must pass full array with haloes!'
    write(0,*) 'hchksum2d: xs,xe,ys,ye=',xs,xe,ys,ye
    write(0,*) 'hchksum2d: start_x,end_x,start_y,end_y=',start_x,end_x,start_y,end_y
    call chksum_error(FATAL,'Error in hchksum3d '//trim(mesg))
  endif

  bc0=hsubsum3d(array, start_x, end_x, start_y, end_y, 0, 0)
  if (present(haloshift)) then
    if (haloshift==0) then
      if (is_root_pe()) write(0,'(A40,5(1X,A,I12))') mesg,"c=",bc0
      return
    endif
    hshift=haloshift
  else
    hshift=1
  endif

  bcSW=hsubsum3d(array, start_x, end_x, start_y, end_y, -hshift, -hshift)
  bcSE=hsubsum3d(array, start_x, end_x, start_y, end_y, hshift, -hshift)
  bcNW=hsubsum3d(array, start_x, end_x, start_y, end_y, -hshift, hshift)
  bcNE=hsubsum3d(array, start_x, end_x, start_y, end_y, hshift, hshift)

  if (is_root_pe()) &
  write(0,'(A40,5(1X,A,I12))') mesg,"c=",bc0,"sw=",bcSW,"se=",bcSE,"nw=",bcNW,"ne=",bcNE

  contains

  integer function hsubsum3d(array, start_x, end_x, start_y, end_y, di, dj)
  real, dimension(:,:,:) :: array
  integer :: start_x, end_x, start_y, end_y, di, dj
  integer :: bitcount
  integer :: i, j, k, bc
  hsubsum3d = 0
  do k=LBOUND(array,3),UBOUND(array,3)
    do j=start_y+dj,end_y+dj
      do i=start_x+di,end_x+di
        bc = bitcount(abs(array(i,j,k)))
        hsubsum3d = hsubsum3d + bc
      enddo
    enddo
  enddo
  call sum_across_PEs(hsubsum3d)
  end function hsubsum3d
end subroutine hchksum3d

! =====================================================================

!   These are the older version of chksum that do not take the grid staggering
! into account.
subroutine chksum1d(array, mesg, start_x, end_x)

  real, dimension(:) :: array
  character(len=*) :: mesg
  integer, optional :: start_x, end_x

  integer :: xs,xe,i,bc,sum1
  integer :: bitcount
  real :: sum
  real, allocatable :: sum_here(:)
  integer :: pe_num   ! pe number of the data

  xs = LBOUND(array,1) ; xe = UBOUND(array,1)
  if (present(start_x)) xs = start_x
  if (present(end_x)) xe = end_x

  sum = 0.0 ; sum1 = 0
  do i=xs,xe
    sum = sum + array(i)
    bc = bitcount(ABS(array(i)))
    sum1 = sum1 + bc
  enddo

  pe_num = pe_here() + 1 - root_pe()
  allocate(sum_here(num_pes())) ; sum_here(:) = 0.0
  sum_here(pe_num) = sum
  call sum_across_PEs(sum_here,num_pes())
  sum = 0.0
  do i=1,num_pes() ; sum = sum + sum_here(i) ; enddo
!   call sum_across_PEs(sum)
  call sum_across_PEs(sum1)

  if (is_root_pe()) &
    write(*,'(A40,1X,ES22.13,1X,I12)') mesg, sum, sum1

end subroutine chksum1d

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
  do i=xs,xe
    do j=ys,ye
      sum = sum + array(i,j)
      bc = bitcount(abs(array(i,j)))
      sum1 = sum1 + bc
    enddo
  enddo

  pe_num = pe_here() + 1 - root_pe()
  allocate(sum_here(num_pes())) ; sum_here(:) = 0.0
  sum_here(pe_num) = sum
  call sum_across_PEs(sum_here,num_pes())
  sum = 0.0
  do i=1,num_pes() ; sum = sum + sum_here(i) ; enddo
!   call sum_across_PEs(sum)
  call sum_across_PEs(sum1)

  if (is_root_pe()) &
    write(*,'(A40,1X,ES22.13,1X,I12)') &
      mesg, sum, sum1
!    write(*,'(A40,1X,Z16.16,1X,Z16.16,1X,ES25.16,1X,I12)') &
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
  do i=xs,xe
    do j=ys,ye
      do k=zs,ze
        sum = sum + array(i,j,k)
        bc = bitcount(ABS(array(i,j,k)))
        sum1 = sum1 + bc
      enddo
    enddo
  enddo

  pe_num = pe_here() + 1 - root_pe()
  allocate(sum_here(num_pes())) ; sum_here(:) = 0.0
  sum_here(pe_num) = sum
  call sum_across_PEs(sum_here,num_pes())
  sum = 0.0
  do i=1,num_pes() ; sum = sum + sum_here(i) ; enddo
!   call sum_across_PEs(sum)
  call sum_across_PEs(sum1)

  if (is_root_pe()) &
    write(*,'(A40,1X,ES22.13,1X,I12)') &
      mesg, sum, sum1
!    write(*,'(A40,1X,Z16.16,1X,Z16.16,1X,ES25.16,1X,I12)') &
!      mesg, sum, sum1, sum, sum1

end subroutine chksum3d

! =====================================================================

function is_NaN_0d(x)
  real, intent(in) :: x
  logical :: is_NaN_0d
! This subroutine returns .true. if x is a NaN, and .false. otherwise.

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

function is_NaN_1d(x,skip_mpp)
  real, dimension(:), intent(in) :: x
  logical :: is_NaN_1d
  logical, optional :: skip_mpp
! This subroutine returns .true. if any x is a NaN, and .false. otherwise.
  integer :: i, n
  logical :: call_mpp

  n = 0
  do i = LBOUND(x,1), UBOUND(x,1)
    if (is_NaN_0d(x(i))) n = n + 1
  enddo
  call_mpp = .true.
  if (present(skip_mpp)) then
    if (skip_mpp) call_mpp = .false.
  endif
  if (call_mpp) call sum_across_PEs(n)
  is_NaN_1d = .false.
  if (n>0) is_NaN_1d = .true.

end function is_NaN_1d

! =====================================================================

function is_NaN_2d(x)
  real, dimension(:,:), intent(in) :: x
  logical :: is_NaN_2d
! This subroutine returns .true. if any x is a NaN, and .false. otherwise.
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

function is_NaN_3d(x)
  real, dimension(:,:,:), intent(in) :: x
  logical :: is_NaN_3d
! This subroutine returns .true. if any x is a NaN, and .false. otherwise.
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

function totalStuff(G, hThick, stuff)
  type(ocean_grid_type),            intent(in) :: G
  real, dimension(G%isd:,G%jsd:,:), intent(in) :: hThick, stuff
  real                                         :: totalStuff
! This subroutine returns sum over computational domain of hThick*Stuff
  integer :: i, j, k

  totalStuff = 0.
  do k = 1, G%ke ; do j = G%jsc, G%jec ; do i= G%isc, G%iec
    totalStuff = totalStuff + hThick(i,j,k) * stuff(i,j,k) * G%areaT(i,j)
  enddo ; enddo ; enddo
  call sum_across_PEs(totalStuff)

end function totalStuff

! =====================================================================

subroutine totalTandS(G, hThick, temperature, salinity, mesg)
  type(ocean_grid_type),            intent(in) :: G
  real, dimension(G%isd:,G%jsd:,:), intent(in) :: hThick, temperature, salinity
  character(len=*),                 intent(in) :: mesg
! This subroutine display the total thickness, temperature and salinity
! as well as the change since the last call.
! NOTE: This uses "save" data which is not thread safe and is purely for
! extreme debugging without a proper debugger.
  real, save :: totalH = 0., totalT = 0., totalS = 0.
  logical, save :: firstCall = .true.
  real :: thisH, thisT, thisS, delH, delT, delS
  integer :: i, j, k

  thisH = 0.
  do k = 1, G%ke ; do j = G%jsc, G%jec ; do i= G%isc, G%iec
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
  if (is_root_pe()) write(0,'(A,1(A,I9,X),A)') fmsg," c=",bc0,mesg
end subroutine chk_sum_msg1

! =====================================================================

subroutine chk_sum_msg5(fmsg,bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  character(len=*), intent(in) :: fmsg, mesg
  integer,          intent(in) :: bc0,bcSW,bcSE,bcNW,bcNE
  if (is_root_pe()) write(0,'(A,5(A,I9,1X),A)') &
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
