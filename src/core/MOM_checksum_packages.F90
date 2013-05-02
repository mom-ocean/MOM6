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
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs

implicit none ; private

public MOM_state_chksum, MOM_thermo_chksum, MOM_accel_chksum

#include <MOM_memory.h>

contains

! =============================================================================

subroutine MOM_state_chksum(mesg, u, v, h, uh, vh, G, haloshift)
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
  call hchksum(G%H_to_kg_m2*h, mesg//" h",G,haloshift=hs)
  call uchksum(G%H_to_kg_m2*uh, mesg//" uh",G,haloshift=hs)
  call vchksum(G%H_to_kg_m2*vh, mesg//" vh",G,haloshift=hs)
end subroutine MOM_state_chksum

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
    call hchksum(G%kg_m2_to_H*pbce, mesg//" pbce",G,haloshift=0)
  if (present(u_accel_bt)) &
    call uchksum(u_accel_bt, mesg//" u_accel_bt",G,haloshift=0)
  if (present(v_accel_bt)) &
    call vchksum(v_accel_bt, mesg//" v_accel_bt",G,haloshift=0)
end subroutine MOM_accel_chksum

end module MOM_checksum_packages
