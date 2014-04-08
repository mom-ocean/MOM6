module MOM_tracer_registry
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
!* By Robert Hallberg, May 2013                                        *
!*                                                                     *
!*   This module contains the tracer_registry_type and the subroutines *
!* that handle the registration of tracers and related subroutines.    *
!* The primary subroutine, register_tracer, is called to indicate the  *
!* tracers that will be advected and diffused.                         *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_diag_mediator, only : diag_ctrl
use MOM_checksums, only : hchksum
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type

implicit none ; private

#include <MOM_memory.h>

public register_tracer, tracer_registry_init, MOM_tracer_chksum
public add_tracer_diagnostics, add_tracer_OBC_values
public tracer_vertdiff, tracer_registry_end

type, public :: tracer_type
  real, dimension(:,:,:), pointer :: t => NULL()
                     ! The array containing the tracer concentration.
  real :: OBC_inflow_conc = 0.0  ! A tracer concentration for generic inflows.
  real, dimension(:,:,:), pointer :: OBC_in_u => NULL(), OBC_in_v => NULL()
             ! These arrays contain structured values for flow into the domain
             ! that are specified in open boundary conditions through u- and
             ! v- faces of the tracer cell.
  real, dimension(:,:,:), pointer :: ad_x => NULL(), ad_y => NULL()
             ! The arrays in which x- & y- advective fluxes are stored.
  real, dimension(:,:,:), pointer :: df_x => NULL(), df_y => NULL()
             ! The arrays in which x- & y- diffusive fluxes are stored.
  real, dimension(:,:), pointer :: ad2d_x => NULL(), ad2d_y => NULL()
             ! The arrays in which vertically summed x- & y- advective fluxes
             ! are stored in units of CONC m3 s-1..
  real, dimension(:,:), pointer :: df2d_x => NULL(), df2d_y => NULL()
             ! The arrays in which vertically summed x- & y- diffusive fluxes
             ! are stored in units of CONC m3 s-1..
  character(len=32) :: name  ! A tracer name for error messages.
end type tracer_type

type, public :: tracer_registry_type
  integer :: ntr = 0        ! The number of registered tracers.
  type(tracer_type) :: Tr(MAX_FIELDS_)  ! The array of registered tracers.
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
end type tracer_registry_type

contains

subroutine register_tracer(tr1, name, param_file, Reg, ad_x, ad_y, &
                           df_x, df_y, OBC_inflow, OBC_in_u, OBC_in_v, &
                           ad_2d_x, ad_2d_y, df_2d_x, df_2d_y)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), target :: tr1
  character(len=*), intent(in)               :: name
  type(param_file_type), intent(in)          :: param_file
  type(tracer_registry_type), pointer        :: Reg
  real, pointer, dimension(:,:,:), optional  :: ad_x, ad_y, df_x, df_y
  real, intent(in), optional                 :: OBC_inflow
  real, pointer, dimension(:,:,:), optional  :: OBC_in_u, OBC_in_v
  real, dimension(:,:),   pointer, optional  :: ad_2d_x, ad_2d_y, df_2d_x, df_2d_y
! This subroutine registers a tracer to be advected and horizontally
! diffused.

! Arguments: tr1 - The pointer to the tracer, in arbitrary concentration units (CONC).
!  (in)      name - The name to be used in messages about the tracer.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  Reg - A pointer to the tracer registry.
!  (in)      ad_x - An array in which zonal advective fluxes are stored in
!                   units of CONC m3 s-1.
!  (in)      ad_y - An array in which meridional advective fluxes are stored
!                   in units of CONC m3 s-1.
!  (in)      df_x - An array in which zonal diffusive fluxes are stored in
!                   units of CONC m3 s-1.
!  (in)      df_y - An array in which meridional diffusive fluxes are stored
!                   in units of CONC m3 s-1.
!  (in)      OBC_inflow - The value of the tracer for all inflows via the open
!                         boundary conditions for which OBC_in_u or OBC_in_v are
!                         not specified, in the same units as tr (CONC).
!  (in)      OBC_in_u - The value of the tracer at inflows through u-faces of
!                       tracer cells, in the same units as tr (CONC).
!  (in)      OBC_in_v - The value of the tracer at inflows through v-faces of
!                       tracer cells, in the same units as tr (CONC).
!  (in,opt)  ad_2d_x - An array in which the vertically summed zonal advective
!                   fluxes are stored in units of CONC m3 s-1.
!  (in,opt)  ad_2d_y - An array in which the vertically summed meridional advective
!                   fluxes are stored in units of CONC m3 s-1.
!  (in,opt)  df_2d_x - An array in which the vertically summed zonal diffusive
!                   fluxes are stored in units of CONC m3 s-1.
!  (in,opt)  df_2d_y - An array in which the vertically summed meridional diffusive
!                   fluxes are stored in units of CONC m3 s-1.

  integer :: ntr
  type(tracer_type) :: temp
  character(len=256) :: mesg    ! Message for error messages.

  if (.not. associated(Reg)) call tracer_registry_init(param_file, Reg)

  if (Reg%ntr>=MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ in MOM_memory.h to at least ",I3," to allow for &
        &all the tracers being registered via register_tracer.")') Reg%ntr+1
    call MOM_error(FATAL,"MOM register_tracer: "//mesg)
  endif
  Reg%ntr = Reg%ntr + 1
  ntr = Reg%ntr

  Reg%Tr(ntr)%name = trim(name)
  Reg%Tr(ntr)%t => tr1

  if (present(ad_x)) then ; if (associated(ad_x)) Reg%Tr(ntr)%ad_x => ad_x ; endif
  if (present(ad_y)) then ; if (associated(ad_y)) Reg%Tr(ntr)%ad_y => ad_y ; endif
  if (present(df_x)) then ; if (associated(df_x)) Reg%Tr(ntr)%df_x => df_x ; endif
  if (present(df_y)) then ; if (associated(df_y)) Reg%Tr(ntr)%df_y => df_y ; endif
  if (present(OBC_inflow)) Reg%Tr(ntr)%OBC_inflow_conc = OBC_inflow
  if (present(OBC_in_u)) then ; if (associated(OBC_in_u)) &
                                    Reg%Tr(ntr)%OBC_in_u => OBC_in_u ; endif
  if (present(OBC_in_v)) then ; if (associated(OBC_in_v)) &
                                    Reg%Tr(ntr)%OBC_in_v => OBC_in_v ; endif
  if (present(ad_2d_x)) then ; if (associated(ad_2d_x)) Reg%Tr(ntr)%ad2d_x => ad_2d_x ; endif
  if (present(ad_2d_y)) then ; if (associated(ad_2d_y)) Reg%Tr(ntr)%ad2d_y => ad_2d_y ; endif
  if (present(df_2d_x)) then ; if (associated(df_2d_x)) Reg%Tr(ntr)%df2d_x => df_2d_x ; endif
  if (present(df_2d_y)) then ; if (associated(df_2d_y)) Reg%Tr(ntr)%df2d_y => df_2d_y ; endif

end subroutine register_tracer

subroutine add_tracer_OBC_values(name, Reg, OBC_inflow, OBC_in_u, OBC_in_v)
  character(len=*), intent(in)               :: name
  type(tracer_registry_type), pointer        :: Reg
  real, intent(in), optional                 :: OBC_inflow
  real, pointer, dimension(:,:,:), optional  :: OBC_in_u, OBC_in_v
! This subroutine adds open boundary condition concentrations for a tracer that
! has previously been registered by a call to register_tracer.

! Arguments: name - The name of the tracer for which the diagnostic pointers.
!  (in/out)  Reg - A pointer to the tracer registry.
!  (in)      OBC_inflow - The value of the tracer for all inflows via the open
!                         boundary conditions for which OBC_in_u or OBC_in_v are
!                         not specified, in the same units as tr (CONC).
!  (in)      OBC_in_u - The value of the tracer at inflows through u-faces of
!                       tracer cells, in the same units as tr (CONC).
!  (in)      OBC_in_v - The value of the tracer at inflows through v-faces of
  integer :: m

  if (.not. associated(Reg)) call MOM_error(FATAL, "add_tracer_OBC_values :"// &
       "register_tracer must be called before add_tracer_OBC_values")

  do m=1,Reg%ntr ; if (Reg%Tr(m)%name == trim(name)) exit ; enddo

  if (m <= Reg%ntr) then
    if (present(OBC_inflow)) Reg%Tr(m)%OBC_inflow_conc = OBC_inflow
    if (present(OBC_in_u)) then ; if (associated(OBC_in_u)) &
                                      Reg%Tr(m)%OBC_in_u => OBC_in_u ; endif
    if (present(OBC_in_v)) then ; if (associated(OBC_in_v)) &
                                      Reg%Tr(m)%OBC_in_v => OBC_in_v ; endif
  else
    call MOM_error(FATAL, "MOM_tracer: register_tracer must be called for "//&
             trim(name)//" before add_tracer_OBC_values is called for it.")
  endif

end subroutine add_tracer_OBC_values

subroutine add_tracer_diagnostics(name, Reg, ad_x, ad_y, df_x, df_y, &
                                  ad_2d_x, ad_2d_y, df_2d_x, df_2d_y)
  character(len=*), intent(in)              :: name
  type(tracer_registry_type), pointer       :: Reg
  real, dimension(:,:,:), pointer, optional :: ad_x, ad_y, df_x, df_y
  real, dimension(:,:),   pointer, optional :: ad_2d_x, ad_2d_y, df_2d_x, df_2d_y
! This subroutine adds diagnostic arrays for a tracer that has previously been
! registered by a call to register_tracer.

! Arguments: name - The name of the tracer for which the diagnostic pointers.
!  (in/out)  Reg - A pointer to the tracer registry.
!  (in,opt)  ad_x - An array in which zonal advective fluxes are stored in
!                   units of CONC m3 s-1.
!  (in,opt)  ad_y - An array in which meridional advective fluxes are stored
!                   in units of CONC m3 s-1.
!  (in,opt)  df_x - An array in which zonal diffusive fluxes are stored in
!                   units of CONC m3 s-1.
!  (in,opt)  df_y - An array in which meridional diffusive fluxes are stored
!                   in units of CONC m3 s-1.
!  (in,opt)  ad_2d_x - An array in which the vertically summed zonal advective
!                   fluxes are stored in units of CONC m3 s-1.
!  (in,opt)  ad_2d_y - An array in which the vertically summed meridional advective
!                   fluxes are stored in units of CONC m3 s-1.
!  (in,opt)  df_2d_x - An array in which the vertically summed zonal diffusive
!                   fluxes are stored in units of CONC m3 s-1.
!  (in,opt)  df_2d_y - An array in which the vertically summed meridional diffusive
!                   fluxes are stored in units of CONC m3 s-1.
  integer :: m

  if (.not. associated(Reg)) call MOM_error(FATAL, "add_tracer_diagnostics: "// &
       "register_tracer must be called before add_tracer_diagnostics")

  do m=1,Reg%ntr ; if (Reg%Tr(m)%name == trim(name)) exit ; enddo

  if (m <= Reg%ntr) then
    if (present(ad_x)) then ; if (associated(ad_x)) Reg%Tr(m)%ad_x => ad_x ; endif
    if (present(ad_y)) then ; if (associated(ad_y)) Reg%Tr(m)%ad_y => ad_y ; endif
    if (present(df_x)) then ; if (associated(df_x)) Reg%Tr(m)%df_x => df_x ; endif
    if (present(df_y)) then ; if (associated(df_y)) Reg%Tr(m)%df_y => df_y ; endif
    if (present(ad_2d_x)) then ; if (associated(ad_2d_x)) Reg%Tr(m)%ad2d_x => ad_2d_x ; endif
    if (present(ad_2d_y)) then ; if (associated(ad_2d_y)) Reg%Tr(m)%ad2d_y => ad_2d_y ; endif
    if (present(df_2d_x)) then ; if (associated(df_2d_x)) Reg%Tr(m)%df2d_x => df_2d_x ; endif
    if (present(df_2d_y)) then ; if (associated(df_2d_y)) Reg%Tr(m)%df2d_y => df_2d_y ; endif
  else
    call MOM_error(FATAL, "MOM_tracer: register_tracer must be called for "//&
             trim(name)//" before add_tracer_diagnostics is called for it.")
  endif

end subroutine add_tracer_diagnostics

subroutine tracer_vertdiff(h_old, ea, eb, dt, tr, G, &
                           sfc_flux, btm_flux, btm_reservoir, sink_rate)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h_old, ea, eb
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: tr
  real,                                  intent(in)    :: dt
  type(ocean_grid_type),                 intent(in)    :: G
  real, dimension(NIMEM_,NJMEM_), optional, intent(in) :: sfc_flux
  real, dimension(NIMEM_,NJMEM_), optional, intent(in) :: btm_flux
  real, dimension(NIMEM_,NJMEM_), optional, intent(inout) :: btm_reservoir
  real,                           optional, intent(in) :: sink_rate
! Arguments: h_old -  Layer thickness before entrainment, in m or kg m-2.
!  (in)      ea - The amount of fluid entrained from the layer above, in the
!                 same units as h_old, i.e. m or kg m-2.
!  (in)      eb - The amount of fluid entrained from the layer below, in the
!                 same units as h_old, i.e. m or kg m-2
!  (inout)   tr - The tracer concentration, in concentration units (CU).
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!  (in,opt)  sfc_flux - The surface flux of the tracer, in CU kg m-2 s-1.
!  (in,opt)  btm_flux - The (negative upward) bottom flux of the tracer,
!                       in units of CU kg m-2 s-1.
!  (inout,opt) btm_reservoir - The amount of tracer in a bottom reservoir, in
!                              units of CU kg m-2. (was CU m)
!  (in,opt)  sink_rate - The rate at which the tracer sinks, in m s-1.

!   This subroutine solves a tridiagonal equation for the final tracer
! concentrations after the dual-entrainments, and possibly sinking or surface
! and bottom sources, are applied.  The sinking is implemented with an
! fully implicit upwind advection scheme.
 
  real :: sink_dist ! The distance the tracer sinks in a time step, in m or kg m-2.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    sfc_src, &      ! The time-integrated surface source of the tracer, in
                    ! units of m or kg m-2 times a concentration.
    btm_src         ! The time-integrated bottom source of the tracer, in
                    ! units of m or kg m-2  times a concentration.
  real, dimension(SZI_(G)) :: &
    b1, &           ! b1 is used by the tridiagonal solver, in m-1 or m2 kg-1.
    d1              ! d1=1-c1 is used by the tridiagonal solver, nondimensional.
  real :: c1(SZI_(G),SZK_(G))     ! c1 is used by the tridiagonal solver, ND.
  real :: h_minus_dsink(SZI_(G),SZK_(G))  ! The layer thickness minus the
                    ! difference in sinking rates across the layer, in m or kg m-2.
                    ! By construction, 0 <= h_minus_dsink < h_old.
  real :: sink(SZI_(G),SZK_(G)+1) ! The tracer's sinking distances at the
                    ! interfaces, limited to prevent characteristics from
                    ! crossing within a single timestep, in m or kg m-2.
  real :: b_denom_1 ! The first term in the denominator of b1, in m or kg m-2.
  real :: h_tr      ! h_tr is h at tracer points with a h_neglect added to
                    ! ensure positive definiteness, in m or kg m-2.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in m.
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  h_neglect = G%H_subroundoff

  sink_dist = 0.0
  if (present(sink_rate)) sink_dist = (dt*sink_rate) * G%m_to_H
  do j=js,je; do i=is,ie ; sfc_src(i,j) = 0.0 ; btm_src(i,j) = 0.0 ; enddo; enddo
  if (present(sfc_flux)) then
     do j = js, je; do i = is,ie
        sfc_src(i,j) = (sfc_flux(i,j)*dt) * G%kg_m2_to_H
     enddo; enddo
  endif
  if (present(btm_flux)) then
     do j = js, je; do i = is,ie
        btm_src(i,j) = (btm_flux(i,j)*dt) * G%kg_m2_to_H
     enddo; enddo
  endif


  if (present(sink_rate)) then
!$OMP parallel do default(shared) private(sink,h_minus_dsink,b_denom_1,b1,d1,h_tr,c1)
    do j=js,je
      ! Find the sinking rates at all interfaces, limiting them if necesary
      ! so that the characteristics do not cross within a timestep.
      !   If a non-constant sinking rate were used, that would be incorprated
      ! here.
      if (present(btm_reservoir)) then
        do i=is,ie ; sink(i,nz+1) = sink_dist ; enddo
        do k=2,nz ; do i=is,ie
          sink(i,K) = sink_dist ; h_minus_dsink(i,k) = h_old(i,j,k)
        enddo ; enddo
      else
        do i=is,ie ; sink(i,nz+1) = 0.0 ; enddo
        ! Find the limited sinking distance at the interfaces.
        do k=nz,2,-1 ; do i=is,ie
          if (sink(i,K+1) >= sink_dist) then
            sink(i,K) = sink_dist
            h_minus_dsink(i,k) = h_old(i,j,k) + (sink(i,K+1) - sink(i,K))
          elseif (sink(i,K+1) + h_old(i,j,k) < sink_dist) then
            sink(i,K) = sink(i,K+1) + h_old(i,j,k)
            h_minus_dsink(i,k) = 0.0
          else
            sink(i,K) = sink_dist
            h_minus_dsink(i,k) = (h_old(i,j,k) + sink(i,K+1)) - sink(i,K)
          endif
        enddo ; enddo
      endif
      do i=is,ie
        sink(i,1) = 0.0 ; h_minus_dsink(i,1) = (h_old(i,j,1) + sink(i,2))
      enddo

      ! Now solve the tridiagonal equation for the tracer concentrations.
      do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
        b_denom_1 = h_minus_dsink(i,1) + ea(i,j,1) + h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,1))
        d1(i) = b_denom_1 * b1(i)
        h_tr = h_old(i,j,1) + h_neglect
        tr(i,j,1) = b1(i)*(h_tr*tr(i,j,1) + sfc_src(i,j))
      endif ; enddo
      do k=2,nz-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
        c1(i,k) = eb(i,j,k-1) * b1(i)
        b_denom_1 = h_minus_dsink(i,k) + d1(i) * (ea(i,j,k) + sink(i,K)) + &
                    h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
        d1(i) = b_denom_1 * b1(i)
        h_tr = h_old(i,j,k) + h_neglect
        tr(i,j,k) = b1(i) * (h_tr * tr(i,j,k) + &
                             (ea(i,j,k) + sink(i,K)) * tr(i,j,k-1))
      endif ; enddo ; enddo
      do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
        c1(i,nz) = eb(i,j,nz-1) * b1(i)
        b_denom_1 = h_minus_dsink(i,nz) + d1(i) * (ea(i,j,nz) + sink(i,nz)) + &
                    h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,nz))
        h_tr = h_old(i,j,nz) + h_neglect
        tr(i,j,nz) = b1(i) * ((h_tr * tr(i,j,nz) + btm_src(i,j)) + &
                              (ea(i,j,nz) + sink(i,nz)) * tr(i,j,nz-1))
      endif ; enddo
      if (present(btm_reservoir)) then ; do i=is,ie ; if (G%mask2dT(i,j)>0.5) then
        btm_reservoir(i,j) = btm_reservoir(i,j) + &
                             (sink(i,nz+1)*tr(i,j,nz)) * G%H_to_kg_m2
      endif ; enddo ; endif

      do k=nz-1,1,-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
        tr(i,j,k) = tr(i,j,k) + c1(i,k+1)*tr(i,j,k+1)
      endif ; enddo ; enddo
    enddo
  else
!$OMP parallel do default(shared) private(h_tr,b_denom_1,b1,d1,c1)
    do j=js,je
      do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
        h_tr = h_old(i,j,1) + h_neglect
        b_denom_1 = h_tr + ea(i,j,1)
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,1))
        d1(i) = b_denom_1 * b1(i)
        tr(i,j,1) = b1(i)*(h_tr*tr(i,j,1) + sfc_src(i,j))
       endif 
      enddo
      do k=2,nz-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
        c1(i,k) = eb(i,j,k-1) * b1(i)
        h_tr = h_old(i,j,k) + h_neglect
        b_denom_1 = h_tr + d1(i) * ea(i,j,k)
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
        d1(i) = b_denom_1 * b1(i)
        tr(i,j,k) = b1(i) * (h_tr * tr(i,j,k) + ea(i,j,k) * tr(i,j,k-1))
      endif ; enddo ; enddo
      do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
        c1(i,nz) = eb(i,j,nz-1) * b1(i)
        h_tr = h_old(i,j,nz) + h_neglect
        b1(i) = 1.0 / (h_tr + d1(i) * ea(i,j,nz) + eb(i,j,nz))
        tr(i,j,nz) = b1(i) * ((h_tr * tr(i,j,nz) + btm_src(i,j)) + &
                              ea(i,j,nz) * tr(i,j,nz-1))
      endif ; enddo
      do k=nz-1,1,-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
        tr(i,j,k) = tr(i,j,k) + c1(i,k+1)*tr(i,j,k+1)
      endif ; enddo ; enddo
    enddo
  endif

end subroutine tracer_vertdiff

subroutine MOM_tracer_chksum(mesg, Tr, ntr, G)
  character(len=*),         intent(in) :: mesg
  type(tracer_type),        intent(in) :: Tr(:)
  integer,                  intent(in) :: ntr
  type(ocean_grid_type),    intent(in) :: G
!   This subroutine writes out chksums for the model's thermodynamic state
! variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      Tr - An array of all of the registered tracers.
!  (in)      ntr - The number of registered tracers.
!  (in)      G - The ocean's grid structure.
  integer :: is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  do m=1,ntr
    call hchksum(Tr(m)%t, mesg//trim(Tr(m)%name), G)
  enddo
end subroutine MOM_tracer_chksum

subroutine tracer_registry_init(param_file, Reg)
  type(param_file_type),      intent(in) :: param_file
  type(tracer_registry_type), pointer    :: Reg
! Arguments: param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  Reg - A pointer that is set to point to the tracer registry.
  integer, save :: init_calls = 0
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_tracer_registry" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (.not.associated(Reg)) then ; allocate(Reg)
  else ; return ; endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")

  init_calls = init_calls + 1
  if (init_calls > 1) then
    write(mesg,'("tracer_registry_init called ",I3, &
      &" times with different registry pointers.")') init_calls
    if (is_root_pe()) call MOM_error(WARNING,"MOM_tracer"//mesg)
  endif

end subroutine tracer_registry_init

subroutine tracer_registry_end(Reg)
  type(tracer_registry_type), pointer :: Reg
  if (associated(Reg)) deallocate(Reg)
end subroutine tracer_registry_end

end module MOM_tracer_registry
