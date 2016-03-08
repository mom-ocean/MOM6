module MOM_checksum_packages

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

!   This module provdes a several routines that do check-sums of groups
! of variables in the various dynamic solver routines.

use MOM_checksums, only : hchksum, uchksum, vchksum
use MOM_checksums, only : totalStuff
use MOM_domains, only : sum_across_PEs, min_across_PEs, max_across_PEs
use MOM_error_handler, only : MOM_mesg, is_root_pe
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs

implicit none ; private

public MOM_state_chksum, MOM_thermo_chksum, MOM_accel_chksum
public MOM_state_stats

interface MOM_state_chksum
  module procedure MOM_state_chksum_5arg
  module procedure MOM_state_chksum_3arg
end interface

#include <MOM_memory.h>

type :: stats
  private
  real :: minimum = 1.E34, maximum = -1.E34, average = 0.
end type stats

contains

! =============================================================================

subroutine MOM_state_chksum_5arg(mesg, u, v, h, uh, vh, G, haloshift)
  character(len=*),                       intent(in) :: mesg
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in) :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in) :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in) :: h
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in) :: uh
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in) :: vh
  type(ocean_grid_type),                  intent(in) :: G
  integer, optional,                      intent(in) :: haloshift
!   This subroutine writes out chksums for the model's basic state variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m.
!  (in)      uh - Volume flux through zonal faces = u*h*dy, m3 s-1.
!  (in)      vh - Volume flux through meridional faces = v*h*dx, in m3 s-1.
!  (in)      G - The ocean's grid structure.
  integer :: is, ie, js, je, nz, hs
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  hs=1; if (present(haloshift)) hs=haloshift
  call uchksum(u, mesg//" u",G,haloshift=hs)
  call vchksum(v, mesg//" v",G,haloshift=hs)
  call hchksum(G%GV%H_to_kg_m2*h, mesg//" h",G,haloshift=hs)
  call uchksum(G%GV%H_to_kg_m2*uh, mesg//" uh",G,haloshift=hs)
  call vchksum(G%GV%H_to_kg_m2*vh, mesg//" vh",G,haloshift=hs)
end subroutine MOM_state_chksum_5arg

! =============================================================================

subroutine MOM_state_chksum_3arg(mesg, u, v, h, G, haloshift)
  character(len=*),                       intent(in) :: mesg
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in) :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in) :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in) :: h
  type(ocean_grid_type),                  intent(in) :: G
  integer, optional,                      intent(in) :: haloshift
!   This subroutine writes out chksums for the model's basic state variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m.
!  (in)      uh - Volume flux through zonal faces = u*h*dy, m3 s-1.
!  (in)      vh - Volume flux through meridional faces = v*h*dx, in m3 s-1.
!  (in)      G - The ocean's grid structure.
  integer :: is, ie, js, je, nz, hs
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  hs=1; if (present(haloshift)) hs=haloshift
  call uchksum(u, mesg//" u",G,haloshift=hs)
  call vchksum(v, mesg//" v",G,haloshift=hs)
  call hchksum(G%GV%H_to_kg_m2*h, mesg//" h",G,haloshift=hs)
end subroutine MOM_state_chksum_3arg

! =============================================================================

subroutine MOM_thermo_chksum(mesg, tv, G, haloshift)
  character(len=*),         intent(in) :: mesg
  type(thermo_var_ptrs),    intent(in) :: tv
  type(ocean_grid_type),    intent(in) :: G
  integer, optional,        intent(in) :: haloshift
!   This subroutine writes out chksums for the model's thermodynamic state
! variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      tv - A structure containing pointers to any thermodynamic
!                 fields that are in use.
!  (in)      G - The ocean's grid structure.
  integer :: is, ie, js, je, nz, hs
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  hs=1; if (present(haloshift)) hs=haloshift

  if (associated(tv%T)) call hchksum(tv%T, mesg//" T",G,haloshift=hs)
  if (associated(tv%S)) call hchksum(tv%S, mesg//" S",G,haloshift=hs)
  if (associated(tv%frazil)) call hchksum(tv%frazil, mesg//" frazil",G,haloshift=hs)
  if (associated(tv%salt_deficit)) call hchksum(tv%salt_deficit, mesg//" salt deficit",G,haloshift=hs)

end subroutine MOM_thermo_chksum

! =============================================================================

subroutine MOM_accel_chksum(mesg, CAu, CAv, PFu, PFv, diffu, diffv, G, pbce, &
                            u_accel_bt, v_accel_bt)
  character(len=*),                       intent(in) :: mesg
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in) :: CAu
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in) :: CAv
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in) :: PFu
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in) :: PFv
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in) :: diffu
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in) :: diffv
  type(ocean_grid_type),                  intent(in) :: G
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  optional, intent(in) :: pbce
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), optional, intent(in) :: u_accel_bt
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), optional, intent(in) :: v_accel_bt
!   This subroutine writes out chksums for the model's accelerations.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      CAu - Zonal acceleration due to Coriolis and momentum
!                  advection terms, in m s-2.
!  (in)      CAv - Meridional acceleration due to Coriolis and
!                  momentum advection terms, in m s-2.
!  (in)      PFu - Zonal acceleration due to pressure gradients
!                  (equal to -dM/dx) in m s-2.
!  (in)      PFv - Meridional acceleration due to pressure
!                  gradients (equal to -dM/dy) in m s-2.
!  (in)      diffu - Zonal acceleration due to convergence of the
!                    along-isopycnal stress tensor, in m s-2.
!  (in)      diffv - Meridional acceleration due to convergence of
!                    the along-isopycnal stress tensor, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      pbce - the baroclinic pressure anomaly in each layer
!                   due to free surface height anomalies, in m s-2.
!                   pbce points to a space with nz layers or NULL.
!  (in)      u_accel_bt - The zonal acceleration from terms in the barotropic
!                         solver, in m s-2.
!  (in)      v_accel_bt - The meridional acceleration from terms in the
!                         barotropic solver, in m s-2.
  integer :: is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  call uchksum(CAu, mesg//" CAu",G,haloshift=0)
  call vchksum(CAv, mesg//" CAv",G,haloshift=0)
  call uchksum(PFu, mesg//" PFu",G,haloshift=0)
  call vchksum(PFv, mesg//" PFv",G,haloshift=0)
  call uchksum(diffu, mesg//" diffu",G,haloshift=0)
  call vchksum(diffv, mesg//" diffv",G,haloshift=0)
  if (present(pbce)) &
    call hchksum(G%GV%kg_m2_to_H*pbce, mesg//" pbce",G,haloshift=0)
  if (present(u_accel_bt)) &
    call uchksum(u_accel_bt, mesg//" u_accel_bt",G,haloshift=0)
  if (present(v_accel_bt)) &
    call vchksum(v_accel_bt, mesg//" v_accel_bt",G,haloshift=0)
end subroutine MOM_accel_chksum

! =============================================================================

subroutine MOM_state_stats(mesg, u, v, h, Temp, Salt, G, allowChange, permitDiminishing)
  character(len=*),                       intent(in) :: mesg
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in) :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in) :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in) :: h
  real, pointer, dimension(:,:,:),        intent(in) :: Temp, Salt
  type(ocean_grid_type),                  intent(in) :: G
  logical, optional,                      intent(in) :: allowChange, permitDiminishing
!   This subroutine monitors statistics for the model's state variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in) u - Zonal velocity, in m s-1.
!  (in) v - Meridional velocity, in m s-1.
!  (in) h - Layer thickness, in m.
!  (in) T - Temperature, in degree C.
!  (in) S - Salinity, in ppt.
!  (in) G - The ocean's grid structure.
!  (in) allowChange - do not flag an error if the statistics change
!  (in) permitDiminishing - do not flag an error if the extrema are diminishing
  integer :: is, ie, js, je, nz, i, j, k
  real :: Vol, dV, Area, h_minimum
  type(stats) :: T, S, delT, delS
  type(stats), save :: oldT, oldS     ! NOTE: save data is not normally allowed but
  logical, save :: firstCall = .true. ! we use it for debugging purposes here on the
  logical :: do_TS
  real, save :: oldVol                ! assumption we will not turn this on with threads
  character(len=80) :: lMsg
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  do_TS = associated(Temp) .and. associated(Salt)

  ! First collect local stats
  Area = 0. ; Vol = 0.
  do j = js, je ; do i = is, ie
    Area = Area + G%areaT(i,j)
  enddo ; enddo
  T%minimum = 1.E34 ; T%maximum = -1.E34 ; T%average = 0.
  S%minimum = 1.E34 ; S%maximum = -1.E34 ; S%average = 0.
  h_minimum = 1.E34
  do k = 1, nz ; do j = js, je ; do i = is, ie
    if (G%mask2dT(i,j)>0.) then
      dV = G%areaT(i,j)*h(i,j,k) ; Vol = Vol + dV
      if (do_TS .and. h(i,j,k)>0.) then
        T%minimum = min( T%minimum, Temp(i,j,k) ) ; T%maximum = max( T%maximum, Temp(i,j,k) )
        T%average = T%average + dV*Temp(i,j,k)
        S%minimum = min( S%minimum, Salt(i,j,k) ) ; S%maximum = max( S%maximum, Salt(i,j,k) )
        S%average = S%average + dV*Salt(i,j,k)
      endif
      if (h_minimum > h(i,j,k)) h_minimum = h(i,j,k)
    endif
  enddo ; enddo ; enddo
  call sum_across_PEs( Area ) ; call sum_across_PEs( Vol )
  if (do_TS) then
    call min_across_PEs( T%minimum ) ; call max_across_PEs( T%maximum ) ; call sum_across_PEs( T%average )
    call min_across_PEs( S%minimum ) ; call max_across_PEs( S%maximum ) ; call sum_across_PEs( S%average )
    T%average = T%average / Vol ; S%average = S%average / Vol
  endif
  if (is_root_pe()) then
    if (.not.firstCall) then
      dV = Vol - oldVol
      delT%minimum = T%minimum - oldT%minimum ; delT%maximum = T%maximum - oldT%maximum
      delT%average = T%average - oldT%average
      delS%minimum = S%minimum - oldS%minimum ; delS%maximum = S%maximum - oldS%maximum
      delS%average = S%average - oldS%average
      write(lMsg(1:80),'(2(a,es12.4))') 'Mean thickness =',Vol/Area,' frac. delta=',dV/Vol
      call MOM_mesg(lMsg//trim(mesg))
      if (do_TS) then
        write(lMsg(1:80),'(a,3es12.4)') 'Temp min/mean/max =',T%minimum,T%average,T%maximum
        call MOM_mesg(lMsg//trim(mesg))
        write(lMsg(1:80),'(a,3es12.4)') 'delT min/mean/max =',delT%minimum,delT%average,delT%maximum
        call MOM_mesg(lMsg//trim(mesg))
        write(lMsg(1:80),'(a,3es12.4)') 'Salt min/mean/max =',S%minimum,S%average,S%maximum
        call MOM_mesg(lMsg//trim(mesg))
        write(lMsg(1:80),'(a,3es12.4)') 'delS min/mean/max =',delS%minimum,delS%average,delS%maximum
        call MOM_mesg(lMsg//trim(mesg))
      endif
    else
      write(lMsg(1:80),'(a,es12.4)') 'Mean thickness =',Vol/Area
      call MOM_mesg(lMsg//trim(mesg))
      if (do_TS) then
        write(lMsg(1:80),'(a,3es12.4)') 'Temp min/mean/max =',T%minimum,T%average,T%maximum
        call MOM_mesg(lMsg//trim(mesg))
        write(lMsg(1:80),'(a,3es12.4)') 'Salt min/mean/max =',S%minimum,S%average,S%maximum
        call MOM_mesg(lMsg//trim(mesg))
      endif
    endif
  endif
  firstCall = .false. ; oldVol = Vol
  oldT%minimum = T%minimum ; oldT%maximum = T%maximum ; oldT%average = T%average
  oldS%minimum = S%minimum ; oldS%maximum = S%maximum ; oldS%average = S%average

  if (do_TS .and. T%minimum<-5.0) then
    do j = js, je ; do i = is, ie
      if (minval(Temp(i,j,:)) == T%minimum) then
        write(0,'(a,2f12.5)') 'x,y=',G%geoLonT(i,j),G%geoLatT(i,j)
        write(0,'(a3,3a12)') 'k','h','Temp','Salt'
        do k = 1, nz
          write(0,'(i3,3es12.4)') k,h(i,j,k),Temp(i,j,k),Salt(i,j,k)
        enddo
        stop 'Extremum detected'
      endif
    enddo ; enddo
  endif

  if (h_minimum<0.0) then
    do j = js, je ; do i = is, ie
      if (minval(h(i,j,:)) == h_minimum) then
        write(0,'(a,2f12.5)') 'x,y=',G%geoLonT(i,j),G%geoLatT(i,j)
        write(0,'(a3,3a12)') 'k','h','Temp','Salt'
        do k = 1, nz
          write(0,'(i3,3es12.4)') k,h(i,j,k),Temp(i,j,k),Salt(i,j,k)
        enddo
        stop 'Negative thickness detected'
      endif
    enddo ; enddo
  endif

end subroutine MOM_state_stats

end module MOM_checksum_packages
