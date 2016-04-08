module MOM_continuity_PPM
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
!*  By Robert Hallberg and Alistair Adcroft, September 2006 - .        *
!*                                                                     *
!*    This program contains the subroutine that advects layer          *
!*  thickness.  The scheme here uses a Piecewise-Parabolic method with *
!*  a positive definite limiter.                                       *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, vh                                    *
!*    j    x ^ x ^ x   At >:  u, uh                                    *
!*    j    > o > o >   At o:  h, hin                                   *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : ocean_OBC_type, BT_cont_type, OBC_SIMPLE
use MOM_variables, only : OBC_FLATHER_E, OBC_FLATHER_W, OBC_FLATHER_N, OBC_FLATHER_S
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public continuity_PPM, continuity_PPM_init, continuity_PPM_end

integer :: id_clock_update, id_clock_correct

type, public :: continuity_PPM_CS ; private
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  logical :: upwind_1st      ! If true, use a first-order upwind scheme.
  logical :: monotonic       ! If true, use the Colella & Woodward monotonic
                             ! limiter; otherwise use a simple positive
                             ! definite limiter.
  logical :: simple_2nd      ! If true, use a simple second order (arithmetic
                             ! mean) interpolation of the edge values instead
                             ! of the higher order interpolation.
  real :: tol_eta            ! The tolerance for free-surface height
                             ! discrepancies between the barotropic solution and
                             ! the sum of the layer thicknesses, in m.
  real :: tol_vel            ! The tolerance for barotropic velocity
                             ! discrepancies between the barotropic solution and
                             ! the sum of the layer thicknesses, in m s-1.
  real :: tol_eta_aux        ! The tolerance for free-surface height
                             ! discrepancies between the barotropic solution and
                             ! the sum of the layer thicknesses when calculating
                             ! the auxiliary corrected velocities, in m.
  real :: CFL_limit_adjust   ! The maximum CFL of the adjusted velocities, ND.
  logical :: aggress_adjust  ! If true, allow the adjusted velocities to have a
                             ! relative CFL change up to 0.5.  False by default.
  logical :: vol_CFL         ! If true, use the ratio of the open face lengths
                             ! to the tracer cell areas when estimating CFL
                             ! numbers.  Without aggress_adjust, the default is
                             ! false; it is always true with.
  logical :: better_iter     ! If true, stop corrective iterations using a
                             ! velocity-based criterion and only stop if the
                             ! iteration is better than all predecessors.
  logical :: use_visc_rem_max ! If true, use more appropriate limiting bounds
                             ! for corrections in strongly viscous columns.
  logical :: marginal_faces  ! If true, use the marginal face areas from the
                             ! continuity solver for use as the weights in the
                             ! barotropic solver.  Otherwise use the transport
                             ! averaged areas.
end type continuity_PPM_CS

type :: loop_bounds_type ; private
  integer :: ish, ieh, jsh, jeh
end type loop_bounds_type

contains

subroutine continuity_PPM(u, v, hin, h, uh, vh, dt, G, GV, CS, uhbt, vhbt, OBC, &
                          visc_rem_u, visc_rem_v, u_cor, v_cor, &
                          uhbt_aux, vhbt_aux, u_cor_aux, v_cor_aux, BT_cont)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in)    :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in)    :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: hin
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(inout) :: h
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(out)   :: uh
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(out)   :: vh
  real,                                   intent(in)    :: dt
  type(ocean_grid_type),                  intent(inout) :: G
  type(continuity_PPM_CS),                pointer       :: CS
  type(verticalGrid_type),                intent(in)    :: GV
  real, dimension(NIMEMB_,NJMEM_),        intent(in),  optional :: uhbt
  real, dimension(NIMEM_,NJMEMB_),        intent(in),  optional :: vhbt
  type(ocean_OBC_type),                   pointer,     optional :: OBC
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in),  optional :: visc_rem_u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in),  optional :: visc_rem_v
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(out), optional :: u_cor
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(out), optional :: v_cor
  real, dimension(NIMEMB_,NJMEM_),        intent(in),  optional :: uhbt_aux
  real, dimension(NIMEM_,NJMEMB_),        intent(in),  optional :: vhbt_aux
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(out), optional :: u_cor_aux
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(out), optional :: v_cor_aux
  type(BT_cont_type),                     pointer,     optional :: BT_cont
!    This subroutine time steps the layer thicknesses, using a monotonically
!  limit, directionally split PPM scheme, based on Lin (1994).  In the following
!  documentation, H is used for the units of thickness (usually m or kg m-2.)

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      hin - Initial layer thickness, in H.
!  (out)     h - Final layer thickness, in H.
!  (out)     uh - Volume flux through zonal faces = u*h*dy, H m2 s-1.
!  (out)     vh - Volume flux through meridional faces = v*h*dx,
!                  in H m2 s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 continuity_PPM_init.
!  (in, opt) uhbt - The summed volume flux through zonal faces, H m2 s-1.
!  (in, opt) vhbt - The summed volume flux through meridional faces, H m2 s-1.
!  (in, opt) OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in, opt) visc_rem_u - Both the fraction of the momentum originally in a
!  (in, opt) visc_rem_v - layer that remains after a time-step of viscosity,
!                         and the fraction of a time-step's worth of a
!                         barotropic acceleration that a layer experiences
!                         after viscosity is applied, in the zonal (_u) and
!                         meridional (_v) directions.  Nondimensional between
!                         0 (at the bottom) and 1 (far above the bottom).
!  (out, opt) u_cor - The zonal velocities that give uhbt as the depth-
!                     integrated transport, in m s-1.
!  (out, opt) v_cor - The meridional velocities that give vhbt as the
!                     depth-integrated transport, in m s-1.
!  (in, opt) uhbt_aux - A second set of summed volume fluxes through zonal
!  (in, opt) vhbt_aux - and meridional faces, both in H m2 s-1.
!  (out, opt) u_cor_aux - The zonal and meridional velocities that give uhbt_aux
!  (out, opt) v_cor_aux - and vhbt_aux as the depth-integrated transports,
!                         both in m s-1.
!  (out, opt) BT_cont - A structure with elements that describe the effective
!                       open face areas as a function of barotropic flow.

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    h_input      ! Left and right face thicknesses, in H.
  real :: h_min
  type(loop_bounds_type) :: LB
  integer :: is, ie, js, je, nz, stensil
  integer :: i, j, k

  logical :: apply_OBC_u, apply_OBC_v, x_first
  logical :: apply_OBC_u_flather_east, apply_OBC_u_flather_west
  logical :: apply_OBC_v_flather_north, apply_OBC_v_flather_south
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  h_min = GV%Angstrom

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "MOM_continuity_PPM: Module must be initialized before it is used.")
  x_first = (MOD(G%first_direction,2) == 0)
  if (present(OBC)) then ; if (associated(OBC)) then
    apply_OBC_u = OBC%apply_OBC_u ; apply_OBC_v = OBC%apply_OBC_v
  endif ; endif

  apply_OBC_u_flather_east = .false. ; apply_OBC_u_flather_west = .false.
  apply_OBC_v_flather_north = .false. ; apply_OBC_v_flather_south = .false.
  if (present(OBC)) then ; if (associated(OBC)) then
    apply_OBC_u_flather_east = OBC%apply_OBC_u_flather_east
    apply_OBC_u_flather_west = OBC%apply_OBC_u_flather_west
    apply_OBC_v_flather_north = OBC%apply_OBC_v_flather_north
    apply_OBC_v_flather_south = OBC%apply_OBC_v_flather_south
    !   If an OBC is being applied, copy the input thicknesses so that the
    ! OBC code works even if hin == h.
    if (apply_OBC_u_flather_east .or. apply_OBC_u_flather_west .or. &
        apply_OBC_v_flather_north .or. apply_OBC_v_flather_south) &
      h_input(:,:,:) = hin(:,:,:)
  endif ; endif

  if (present(visc_rem_u) .neqv. present(visc_rem_v)) call MOM_error(FATAL, &
      "MOM_continuity_PPM: Either both visc_rem_u and visc_rem_v or neither"// &
      " one must be present in call to continuity_PPM.")

  stensil = 3 ; if (CS%simple_2nd) stensil = 2 ; if (CS%upwind_1st) stensil = 1

  if (x_first) then
  !    First, advect zonally.
    LB%ish = G%isc ; LB%ieh = G%iec
    LB%jsh = G%jsc-stensil ; LB%jeh = G%jec+stensil
    call zonal_mass_flux(u, hin, uh, dt, G, GV, CS, LB, uhbt, OBC, visc_rem_u, &
                         u_cor, uhbt_aux, u_cor_aux, BT_cont)

    call cpu_clock_begin(id_clock_update)
!$OMP parallel do default(none) shared(LB,nz,G,uh,hin,dt,h)
    do k=1,nz ; do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j,k) = hin(i,j,k) - dt* G%IareaT(i,j) * (uh(I,j,k) - uh(I-1,j,k))
  !   Uncomment this line to prevent underflow.
  !   if (h(i,j,k) < h_min) h(i,j,k) = h_min
    enddo ; enddo ; enddo
    call cpu_clock_end(id_clock_update)

    if (apply_OBC_u_flather_east .or. apply_OBC_u_flather_west) then
      do k=1,nz ; do j=LB%jsh,LB%jeh
        do I=LB%ish,LB%ieh+1
          if (OBC%OBC_mask_u(I-1,j) .and. (OBC%OBC_kind_u(I-1,j) == OBC_FLATHER_E)) &
            h(i,j,k) = h_input(i-1,j,k)
        enddo
        do i=LB%ish-1,LB%ieh
          if (OBC%OBC_mask_u(I,j) .and. (OBC%OBC_kind_u(I,j) == OBC_FLATHER_W)) &
            h(i,j,k) = h_input(i+1,j,k)
        enddo
      enddo ; enddo
    endif
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec

  !    Now advect meridionally, using the updated thicknesses to determine
  !  the fluxes.
    call meridional_mass_flux(v, h, vh, dt, G, GV, CS, LB, vhbt, OBC, visc_rem_v, &
                              v_cor, vhbt_aux, v_cor_aux, BT_cont)

    call cpu_clock_begin(id_clock_update)
!$OMP parallel do default(none) private(h_min) shared(nz,LB,h,dt,G,vh)
    do k=1,nz ; do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j,k) = h(i,j,k) - dt*G%IareaT(i,j) * (vh(i,J,k) - vh(i,J-1,k))
  !   This line prevents underflow.
      if (h(i,j,k) < h_min) h(i,j,k) = h_min
    enddo ; enddo ; enddo
    call cpu_clock_end(id_clock_update)

    if (apply_OBC_v_flather_north .or. apply_OBC_v_flather_south) then
      do k=1,nz
        do J=LB%jsh,LB%jeh+1 ; do i=LB%ish-1,LB%ieh+1
          if (OBC%OBC_mask_v(i,J-1) .and. (OBC%OBC_kind_v(i,J-1) == OBC_FLATHER_N)) &
            h(i,j,k) = h_input(i,j-1,k)
        enddo ; enddo
        do J=LB%jsh-1,LB%jeh ; do i=LB%ish-1,LB%ieh+1
          if (OBC%OBC_mask_v(i,J) .and. (OBC%OBC_kind_v(i,J) == OBC_FLATHER_S)) &
            h(i,j,k) = h_input(i,j+1,k)
        enddo ; enddo
      enddo
    endif
  else  ! .not. x_first
  !    First, advect meridionally, so set the loop bounds accordingly.
    LB%ish = G%isc-stensil ; LB%ieh = G%iec+stensil
    LB%jsh = G%jsc ; LB%jeh = G%jec

    call meridional_mass_flux(v, hin, vh, dt, G, GV, CS, LB, vhbt, OBC, visc_rem_v, &
                              v_cor, vhbt_aux, v_cor_aux, BT_cont)

    call cpu_clock_begin(id_clock_update)
!$OMP parallel do default(none) shared(nz,LB,h,hin,dt,G,vh)
    do k=1,nz ; do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j,k) = hin(i,j,k) - dt*G%IareaT(i,j) * (vh(i,J,k) - vh(i,J-1,k))
    enddo ; enddo ; enddo
    call cpu_clock_end(id_clock_update)

    if (apply_OBC_v_flather_north .or. apply_OBC_v_flather_south) then
      do k=1,nz
        do J=LB%jsh,LB%jeh+1 ; do i=LB%ish-1,LB%ieh+1
          if (OBC%OBC_mask_v(i,J-1) .and. (OBC%OBC_kind_v(i,J-1) == OBC_FLATHER_N)) &
            h(i,j,k) = h_input(i,j-1,k)
        enddo ; enddo
        do J=LB%jsh-1,LB%jeh ; do i=LB%ish-1,LB%ieh+1
          if (OBC%OBC_mask_v(i,J) .and. (OBC%OBC_kind_v(i,J) == OBC_FLATHER_S)) &
            h(i,j,k) = h_input(i,j+1,k)
        enddo ; enddo
      enddo
    endif

  !    Now advect zonally, using the updated thicknesses to determine
  !  the fluxes.
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
    call zonal_mass_flux(u, h, uh, dt, G, GV, CS, LB, uhbt, OBC, visc_rem_u, &
                         u_cor, uhbt_aux, u_cor_aux, BT_cont)

    call cpu_clock_begin(id_clock_update)
!$OMP parallel do default(none) private(h_min) shared(nz,LB,h,dt,G,uh)
    do k=1,nz ; do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j,k) = h(i,j,k) - dt* G%IareaT(i,j) * (uh(I,j,k) - uh(I-1,j,k))
  !   This line prevents underflow.
      if (h(i,j,k) < h_min) h(i,j,k) = h_min
    enddo ; enddo ; enddo
    call cpu_clock_end(id_clock_update)

    if (apply_OBC_u_flather_east .or. apply_OBC_u_flather_west) then
      do k=1,nz ; do j=LB%jsh,LB%jeh
        do I=LB%ish,LB%ieh+1
          if (OBC%OBC_mask_u(I-1,j) .and. (OBC%OBC_kind_u(I-1,j) == OBC_FLATHER_E)) &
            h(i,j,k) = h_input(i-1,j,k)
        enddo
        do i=LB%ish-1,LB%ieh
          if (OBC%OBC_mask_u(I,j) .and. (OBC%OBC_kind_u(I,j) == OBC_FLATHER_W)) &
            h(i,j,k) = h_input(i+1,j,k)
        enddo
      enddo ; enddo
    endif
  endif

end subroutine continuity_PPM

subroutine zonal_mass_flux(u, h_in, uh, dt, G, GV, CS, LB, uhbt, OBC, &
                           visc_rem_u, u_cor, uhbt_aux, u_cor_aux, BT_cont)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in)    :: u
  real,  dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h_in
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(out)   :: uh
  real,                                   intent(in)    :: dt
  type(ocean_grid_type),                  intent(inout) :: G
  type(verticalGrid_type),                intent(in)    :: GV
  type(continuity_PPM_CS),                pointer       :: CS
  type(loop_bounds_type),                 intent(in)    :: LB
  type(ocean_OBC_type),                   pointer,     optional :: OBC
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in),  optional :: visc_rem_u
  real, dimension(NIMEMB_,NJMEM_),        intent(in),  optional :: uhbt, uhbt_aux
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(out), optional :: u_cor, u_cor_aux
  type(BT_cont_type),                     pointer,     optional :: BT_cont
!   This subroutine calculates the mass or volume fluxes through the zonal
! faces, and other related quantities.
! Arguments: u - Zonal velocity, in m s-1.
!  (in)      h_in - Layer thickness used to calculate the fluxes, in H.
!  (out)     uh - Volume flux through zonal faces = u*h*dy, H m2 s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 continuity_PPM_init.
!  (in)      LB - A structure with the active loop bounds.
!  (in, opt) uhbt - The summed volume flux through zonal faces, H m2 s-1.
!  (in, opt) OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in, opt) visc_rem_u - Both the fraction of the momentum originally in a
!                         layer that remains after a time-step of viscosity,
!                         and the fraction of a time-step's worth of a
!                         barotropic acceleration that a layer experiences
!                         after viscosity is applied.  Nondimensional between
!                         0 (at the bottom) and 1 (far above the bottom).
!  (out, opt) u_cor - The zonal velocitiess (u with a barotropic correction)
!                     that give uhbt as the depth-integrated transport, m s-1.
!  (in, opt) uhbt_aux - A second set of summed volume fluxes through zonal
!                       faces, in H m2 s-1.
!  (out, opt) u_cor_aux - The zonal velocities (u with a barotropic correction)
!                         that give uhbt_aux as the depth-integrated transports,
!                         in m s-1.
!  (out, opt) BT_cont - A structure with elements that describe the effective
!                       open face areas as a function of barotropic flow.

  real, dimension(SZIB_(G),SZK_(G)) :: &
    duhdu      ! Partial derivative of uh with u, in H m.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    hl, hr      ! Left and right face thicknesses, in H.
  real, dimension(SZIB_(G)) :: &
    du, &      ! Corrective barotropic change in the velocity, in m s-1.
    du_min_CFL, & ! Min/max limits on du correction
    du_max_CFL, & ! to avoid CFL violations
    duhdu_tot_0, & ! Summed partial derivative of uh with u, in H m.
    uh_tot_0, & ! Summed transport with no barotropic correction in H m2 s-1.
    visc_rem_max  ! The column maximum of visc_rem.
  logical, dimension(SZIB_(G)) :: do_i
  real, dimension(SZIB_(G),SZK_(G)) :: &
    visc_rem      ! A 2-D copy of visc_rem_u or an array of 1's.
  real :: I_vrm   ! 1.0 / visc_rem_max, nondim.
  real :: CFL_dt  ! The maximum CFL ratio of the adjusted velocities divided by
                  ! the time step, in s-1.
  real :: I_dt    ! 1.0 / dt, in s-1.
  real :: du_lim  ! The velocity change that give a relative CFL of 1, in m s-1.
  real :: dx_E, dx_W ! Effective x-grid spacings to the east and west, in m.
  integer :: i, j, k, ish, ieh, jsh, jeh, nz
  logical :: do_aux, apply_OBC_u, use_visc_rem, set_BT_cont, any_simple_OBC

  do_aux = (present(uhbt_aux) .and. present(u_cor_aux))
  use_visc_rem = present(visc_rem_u)
  apply_OBC_u = .false. ; set_BT_cont = .false.
  if (present(BT_cont)) set_BT_cont = (associated(BT_cont))
  if (present(OBC)) then ; if (associated(OBC)) then
    apply_OBC_u = OBC%apply_OBC_u
  endif ; endif
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = G%ke

  CFL_dt = CS%CFL_limit_adjust / dt
  I_dt = 1.0 / dt
  if (CS%aggress_adjust) CFL_dt = I_dt

  call cpu_clock_begin(id_clock_update)
!$OMP parallel do default(none) shared(ish,ieh,jsh,jeh,nz,CS,hl,h_in,hr,G,GV,LB,visc_rem)
  do k=1,nz
    ! This sets hl and hr.
    if (CS%upwind_1st) then
      do j=jsh,jeh ; do i=ish-1,ieh+1
        hl(i,j,k) = h_in(i,j,k) ; hr(i,j,k) = h_in(i,j,k)
      enddo ; enddo
    else
      call PPM_reconstruction_x(h_in(:,:,k), hl(:,:,k), hr(:,:,k), G, LB, &
                                2.0*GV%Angstrom, CS%monotonic, simple_2nd=CS%simple_2nd)
    endif
    do I=ish-1,ieh ; visc_rem(I,k) = 1.0 ; enddo
  enddo
  call cpu_clock_end(id_clock_update)

  call cpu_clock_begin(id_clock_correct)
!$OMP parallel do default(none) shared(ish,ieh,jsh,jeh,nz,u,h_in,hL,hR,use_visc_rem,visc_rem_u,  &
!$OMP                                  uh,dt,G,GV,CS,apply_OBC_u,OBC,uhbt,do_aux,set_BT_cont,    &
!$OMP                                  CFL_dt,I_dt,u_cor,uhbt_aux,u_cor_aux,BT_cont) &
!$OMP                          private(do_i,duhdu,du,du_max_CFL,du_min_CFL,uh_tot_0,duhdu_tot_0, &
!$OMP                                  visc_rem_max, I_vrm, du_lim, dx_E, dx_W, any_simple_OBC ) &
!$OMP      firstprivate(visc_rem)
  do j=jsh,jeh
    do I=ish-1,ieh ; do_i(I) = .true. ; visc_rem_max(I) = 0.0 ; enddo
    ! Set uh and duhdu.
    do k=1,nz
      if (use_visc_rem) then ; do I=ish-1,ieh
        visc_rem(I,k) = visc_rem_u(I,j,k)
        visc_rem_max(I) = max(visc_rem_max(I), visc_rem(I,k))
      enddo ; endif
      call zonal_flux_layer(u(:,j,k), h_in(:,j,k), hL(:,j,k), hR(:,j,k), &
                            uh(:,j,k), duhdu(:,k), visc_rem(:,k), &
                            dt, G, j, ish, ieh, do_i, CS%vol_CFL)
      if (apply_OBC_u) then ; do I=ish-1,ieh
        if (OBC%OBC_mask_u(I,j) .and. (OBC%OBC_kind_u(I,j) == OBC_SIMPLE)) &
          uh(I,j,k) = OBC%uh(I,j,k)
      enddo ; endif
    enddo

    if ((.not.use_visc_rem).or.(.not.CS%use_visc_rem_max)) then ; do I=ish-1,ieh
      visc_rem_max(I) = 1.0
    enddo ; endif

    if (present(uhbt) .or. do_aux .or. set_BT_cont) then
      !   Set limits on du that will keep the CFL number between -1 and 1.
      ! This should be adequate to keep the root bracketed in all cases.
      do I=ish-1,ieh
        I_vrm = 0.0
        if (visc_rem_max(I) > 0.0) I_vrm = 1.0 / visc_rem_max(I)
        if (CS%vol_CFL) then
          dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
          dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
        else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif
        du_max_CFL(I) = 2.0* (CFL_dt * dx_W) * I_vrm
        du_min_CFL(I) = -2.0 * (CFL_dt * dx_E) * I_vrm
        uh_tot_0(I) = 0.0 ; duhdu_tot_0(I) = 0.0
      enddo
      do k=1,nz ; do I=ish-1,ieh
        duhdu_tot_0(I) = duhdu_tot_0(I) + duhdu(I,k)
        uh_tot_0(I) = uh_tot_0(I) + uh(I,j,k)
      enddo ; enddo
      if (use_visc_rem) then
        if (CS%aggress_adjust) then
          do k=1,nz ; do I=ish-1,ieh
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            du_lim = 0.499*((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k)))
            if (du_max_CFL(I) * visc_rem(I,k) > du_lim) &
              du_max_CFL(I) = du_lim / visc_rem(I,k)

            du_lim = 0.499*((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k)))
            if (du_min_CFL(I) * visc_rem(I,k) < du_lim) &
              du_min_CFL(I) = du_lim / visc_rem(I,k)
          enddo ; enddo
        else
          do k=1,nz ; do I=ish-1,ieh
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            if (du_max_CFL(I) * visc_rem(I,k) > dx_W*CFL_dt - u(I,j,k)) &
              du_max_CFL(I) = (dx_W*CFL_dt - u(I,j,k)) / visc_rem(I,k)
            if (du_min_CFL(I) * visc_rem(I,k) < -dx_E*CFL_dt - u(I,j,k)) &
              du_min_CFL(I) = -(dx_E*CFL_dt + u(I,j,k)) / visc_rem(I,k)
          enddo ; enddo
        endif
      else
        if (CS%aggress_adjust) then
          do k=1,nz ; do I=ish-1,ieh
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            du_max_CFL(I) = MIN(du_max_CFL(I), 0.499 * &
                        ((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k))) )
            du_min_CFL(I) = MAX(du_min_CFL(I), 0.499 * &
                        ((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k))) )
          enddo ; enddo
        else
          do k=1,nz ; do I=ish-1,ieh
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            du_max_CFL(I) = MIN(du_max_CFL(I), dx_W*CFL_dt - u(I,j,k))
            du_min_CFL(I) = MAX(du_min_CFL(I), -(dx_E*CFL_dt + u(I,j,k)))
          enddo ; enddo
        endif
      endif
      do I=ish-1,ieh
        du_max_CFL(I) = max(du_max_CFL(I),0.0)
        du_min_CFL(I) = min(du_min_CFL(I),0.0)
      enddo

      ! Up to this point, everything is shared between uhbt and uhbt_aux.

      any_simple_OBC = .false.
      if (present(uhbt) .or. do_aux .or. set_BT_cont) then
        if (.not.apply_OBC_u) then ; do I=ish-1,ieh
          do_i(I) = .true.
        enddo ; else ; do I=ish-1,ieh
          do_i(I) = .not.(OBC%OBC_mask_u(I,j) .and. &
                          (OBC%OBC_kind_u(I,j) == OBC_SIMPLE))
          if (.not.do_i(I)) any_simple_OBC = .true.
        enddo ; endif
      endif

      if (present(uhbt)) then
        call zonal_flux_adjust(u, h_in, hL, hR, uhbt(:,j), uh_tot_0, &
                               duhdu_tot_0, du, du_max_CFL, du_min_CFL, dt, G, &
                               CS, visc_rem, j, ish, ieh, do_I, .true., uh)

        if (present(u_cor)) then ; do k=1,nz
          do I=ish-1,ieh ; u_cor(I,j,k) = u(I,j,k) + du(I) * visc_rem(I,k) ; enddo
          if (apply_OBC_u) then ; do I=ish-1,ieh
            if (OBC%OBC_mask_u(I,j) .and. (OBC%OBC_kind_u(I,j) == OBC_SIMPLE)) &
              u_cor(I,j,k) = OBC%u(I,j,k)
          enddo ; endif
        enddo ; endif ! u-corrected

      endif

      if (do_aux) then
        call zonal_flux_adjust(u, h_in, hL, hR, uhbt_aux(:,j), uh_tot_0, &
                               duhdu_tot_0, du, du_max_CFL, du_min_CFL, dt, G, &
                               CS, visc_rem, j, ish, ieh, do_I, .false.)

        do k=1,nz
          do I=ish-1,ieh ; u_cor_aux(I,j,k) = u(I,j,k) + du(I) * visc_rem(I,k) ; enddo
          if (apply_OBC_u) then ; do I=ish-1,ieh
            if (OBC%OBC_mask_u(I,j) .and. (OBC%OBC_kind_u(I,j) == OBC_SIMPLE)) &
              u_cor_aux(I,j,k) = OBC%u(I,j,k)
          enddo ; endif
        enddo
      endif ! do_aux

      if (set_BT_cont) then
        call set_zonal_BT_cont(u, h_in, hL, hR, BT_cont, uh_tot_0, duhdu_tot_0,&
                               du_max_CFL, du_min_CFL, dt, G, CS, visc_rem, &
                               visc_rem_max, j, ish, ieh, do_i)
        if (any_simple_OBC) then
          do I=ish-1,ieh
            do_i(I) = (OBC%OBC_mask_u(I,j) .and. &
                       (OBC%OBC_kind_u(I,j) == OBC_SIMPLE))
            if (do_i(I)) BT_cont%Fa_u_W0(I,j) = GV%H_subroundoff*G%dy_Cu(I,j)
          enddo
          do k=1,nz ; do I=ish-1,ieh ; if (do_i(I)) then
            if (abs(OBC%u(I,j,k)) > 0.0) &
              BT_cont%Fa_u_W0(I,j) = BT_cont%Fa_u_W0(I,j) + &
                   OBC%uh(I,j,k) / OBC%u(I,j,k)
          endif ; enddo ; enddo
          do I=ish-1,ieh ; if (do_i(I)) then
            BT_cont%Fa_u_E0(I,j) = BT_cont%Fa_u_W0(I,j)
            BT_cont%Fa_u_WW(I,j) = BT_cont%Fa_u_W0(I,j)
            BT_cont%Fa_u_EE(I,j) = BT_cont%Fa_u_W0(I,j)
          endif ; enddo
        endif
      endif ! set_BT_cont

    endif ! present(uhbt) or do_aux or set_BT_cont
  enddo ! j-loop
  call cpu_clock_end(id_clock_correct)

  if  (set_BT_cont) then ; if (associated(BT_cont%h_u)) then
    if (present(u_cor)) then
      call zonal_face_thickness(u_cor, h_in, hL, hR, BT_cont%h_u, dt, G, LB, &
                                CS%vol_CFL, CS%marginal_faces, visc_rem_u)
    else
      call zonal_face_thickness(u, h_in, hL, hR, BT_cont%h_u, dt, G, LB, &
                                CS%vol_CFL, CS%marginal_faces, visc_rem_u)
    endif
  endif ; endif

end subroutine zonal_mass_flux

subroutine zonal_flux_layer(u, h, hL, hR, uh, duhdu, visc_rem, dt, G, j, &
                            ish, ieh, do_i, vol_CFL)
  real, dimension(NIMEMB_),    intent(in)    :: u, visc_rem
  real, dimension(NIMEM_),     intent(in)    :: h, hL, hR
  real, dimension(NIMEMB_),    intent(inout) :: uh, duhdu
  real,                        intent(in)    :: dt
  type(ocean_grid_type),       intent(inout) :: G
  integer,                     intent(in)    :: j, ish, ieh
  logical, dimension(NIMEMB_), intent(in)    :: do_i
  logical,                     intent(in)    :: vol_CFL
!   This subroutines evaluates the zonal mass or volume fluxes in a layer.
!
! Arguments: u - Zonal velocity, in m s-1.
!  (in)      h - Layer thickness used to calculate the fluxes, in H.
!  (in)      hL, hR - Left- and right- thicknesses in the reconstruction, in H.
!  (out)     uh - The zonal mass or volume transport, in H m2 s-1.
!  (out)     duhdu - The partial derivative of uh with u, in H m.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      visc_rem - Both the fraction of the momentum originally in a
!                       layer that remains after a time-step of viscosity,
!                       and the fraction of a time-step's worth of a
!                       barotropic acceleration that a layer experiences
!                       after viscosity is applied.  Nondimensional between
!                       0 (at the bottom) and 1 (far above the bottom).
!  (in)      j, ish, ieh - The index range to work on.
!  (in)      do_i - A logical flag indiciating which I values to work on.
!  (in)      vol_CFL - If true, rescale the ratio of face areas to the cell
!                      areas when estimating the CFL number.
  real :: CFL  ! The CFL number based on the local velocity and grid spacing, ND.
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h_in.
  real :: h_marg ! The marginal thickness of a flux, in H.
  integer :: i

  do I=ish-1,ieh ; if (do_i(I)) then
    ! Set new values of uh and duhdu.
    if (u(I) > 0.0) then
      if (vol_CFL) then ; CFL = (u(I) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
      else ; CFL = u(I) * dt * G%IdxT(i,j) ; endif
      curv_3 = hL(i) + hR(i) - 2.0*h(i)
      uh(I) = G%dy_Cu(I,j) * u(I) * &
          (hR(i) + CFL * (0.5*(hL(i) - hR(i)) + curv_3*(CFL - 1.5)))
      h_marg = hR(i) + CFL * ((hL(i) - hR(i)) + 3.0*curv_3*(CFL - 1.0))
    elseif (u(I) < 0.0) then
      if (vol_CFL) then ; CFL = (-u(I) * dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
      else ; CFL = -u(I) * dt * G%IdxT(i+1,j) ; endif
      curv_3 = hL(i+1) + hR(i+1) - 2.0*h(i+1)
      uh(I) = G%dy_Cu(I,j) * u(I) * &
          (hL(i+1) + CFL * (0.5*(hR(i+1)-hL(i+1)) + curv_3*(CFL - 1.5)))
      h_marg = hL(i+1) + CFL * ((hR(i+1)-hL(i+1)) + 3.0*curv_3*(CFL - 1.0))
    else
      uh(I) = 0.0
      h_marg = 0.5 * (hl(i+1) + hr(i))
    endif
    duhdu(I) = G%dy_Cu(I,j) * h_marg * visc_rem(I)
  endif ; enddo

end subroutine zonal_flux_layer

subroutine zonal_face_thickness(u, h, hL, hR, h_u, dt, G, LB, vol_CFL, &
                                marginal, visc_rem_u)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in)    :: u
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: h, hL, hR
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(inout) :: h_u
  real,                                   intent(in)    :: dt
  type(ocean_grid_type),                  intent(inout) :: G
  type(loop_bounds_type),                 intent(in)    :: LB
  logical,                                intent(in)    :: vol_CFL
  logical,                                intent(in)    :: marginal
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in),  optional :: visc_rem_u
!   This subroutines sets the effective interface thickness at each zonal
! velocity point.
!
! Arguments: u - Zonal velocity, in m s-1.
!  (in)      h - Layer thickness used to calculate the fluxes, in H.
!  (in)      hL, hR - Left- and right- thicknesses in the reconstruction, in H.
!  (out)     h_u - The thickness at zonal faces, in H.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      LB - A structure with the active thickness loop bounds.
!  (in)      vol_CFL - If true, rescale the ratio of face areas to the cell
!                      areas when estimating the CFL number.
!  (in)      marginal - If true, report the marginal face thicknesses; otherwise
!                       report transport-averaged thicknesses.
!  (in, opt) visc_rem_u - Both the fraction of the momentum originally in a
!                         layer that remains after a time-step of viscosity,
!                         and the fraction of a time-step's worth of a
!                         barotropic acceleration that a layer experiences
!                         after viscosity is applied.  Nondimensional between
!                         0 (at the bottom) and 1 (far above the bottom).
  real :: CFL  ! The CFL number based on the local velocity and grid spacing, ND.
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h_in.
  real :: h_avg  ! The average thickness of a flux, in H.
  real :: h_marg ! The marginal thickness of a flux, in H.
  integer :: i, j, k, ish, ieh, jsh, jeh, nz
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = G%ke

!$OMP parallel default(none) shared(ish,ieh,jsh,jeh,nz,u,vol_CFL,dt,G, &
!$OMP                               hL,hR,h,h_u,visc_rem_u,marginal) &
!$OMP                       private(CFL,curv_3,h_marg,h_avg)
!$OMP do
  do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
    if (u(I,j,k) > 0.0) then
      if (vol_CFL) then ; CFL = (u(I,j,k) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
      else ; CFL = u(I,j,k) * dt * G%IdxT(i,j) ; endif
      curv_3 = hL(i,j,k) + hR(i,j,k) - 2.0*h(i,j,k)
      h_avg = hR(i,j,k) + CFL * (0.5*(hL(i,j,k) - hR(i,j,k)) + curv_3*(CFL - 1.5))
      h_marg = hR(i,j,k) + CFL * ((hL(i,j,k) - hR(i,j,k)) + 3.0*curv_3*(CFL - 1.0))
    elseif (u(I,j,k) < 0.0) then
      if (vol_CFL) then ; CFL = (-u(I,j,k)*dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
      else ; CFL = -u(I,j,k) * dt * G%IdxT(i+1,j) ; endif
      curv_3 = hL(i+1,j,k) + hR(i+1,j,k) - 2.0*h(i+1,j,k)
      h_avg = hL(i+1,j,k) + CFL * (0.5*(hR(i+1,j,k)-hL(i+1,j,k)) + curv_3*(CFL - 1.5))
      h_marg = hL(i+1,j,k) + CFL * ((hR(i+1,j,k)-hL(i+1,j,k)) + &
                                    3.0*curv_3*(CFL - 1.0))
    else
      h_avg = 0.5 * (hl(i+1,j,k) + hr(i,j,k))
      !   The choice to use the arithmetic mean here is somewhat arbitrariy, but
      ! it should be noted that hl(i+1,j,k) and hr(i,j,k) are usually the same.
      h_marg = 0.5 * (hl(i+1,j,k) + hr(i,j,k))
 !    h_marg = (2.0 * hl(i+1,j,k) * hr(i,j,k)) / &
 !             (hl(i+1,j,k) + hr(i,j,k) + GV%H_subroundoff)
    endif

    if (marginal) then ; h_u(I,j,k) = h_marg
    else ; h_u(I,j,k) = h_avg ; endif
  enddo; enddo ; enddo
  if (present(visc_rem_u)) then
!$OMP do
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
      h_u(I,j,k) = h_u(I,j,k) * visc_rem_u(I,j,k)
    enddo ; enddo ; enddo
  endif
!$OMP end parallel

end subroutine zonal_face_thickness

subroutine zonal_flux_adjust(u, h_in, hL, hR, uhbt, uh_tot_0, duhdu_tot_0, &
                             du, du_max_CFL, du_min_CFL, dt, G, CS, visc_rem, &
                             j, ish, ieh, do_i_in, full_precision, uh_3d)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in) :: u
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in) :: h_in, hL, hR
  real, dimension(NIMEMB_,NKMEM_),        intent(in) :: visc_rem
  real, dimension(NIMEMB_),            intent(in)    :: uhbt
  real, dimension(NIMEMB_),            intent(in)    :: du_max_CFL, du_min_CFL
  real, dimension(NIMEMB_),            intent(in)    :: uh_tot_0, duhdu_tot_0
  real, dimension(NIMEMB_),            intent(out)   :: du
  real,                                intent(in)    :: dt
  type(ocean_grid_type),               intent(inout) :: G
  type(continuity_PPM_CS),             pointer       :: CS
  integer,                             intent(in)    :: j, ish, ieh
  logical, dimension(NIMEMB_),         intent(in)    :: do_i_in
  logical,                             intent(in),    optional :: full_precision
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(inout), optional :: uh_3d

!   This subroutine returns the barotropic velocity adjustment that gives the
! desired barotropic (layer-summed) transport.
!
! Arguments: u - Zonal velocity, in m s-1.
!  (in)      h_in - Layer thickness used to calculate the fluxes, in H.
!  (in)      hL, hR - Left- and right- thicknesses in the reconstruction, in H.
!  (in)      uhbt - Summed volume flux through zonal faces, in H m2 s-1.
!  (in)      uh_tot_0 - The summed transport with 0 adjustment, in H m2 s-1.
!  (in)      duhdu_tot_0 - The partial derivative of du_err with du at 0
!                          adjustment, in H m.
!  (out)     du - The barotropic velocity adjustment, in m s-1.
!  (in)      du_max_CFL, du_min_CFL - The acceptable range of du, in m s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 continuity_PPM_init.
!  (in)      visc_rem - Both the fraction of the momentum originally in a
!                       layer that remains after a time-step of viscosity,
!                       and the fraction of a time-step's worth of a
!                       barotropic acceleration that a layer experiences
!                       after viscosity is applied.  Nondimensional between
!                       0 (at the bottom) and 1 (far above the bottom).
!  (in)      j, ish, ieh - The index range to work on.
!  (in)      do_I_in - A logical flag indiciating which I values to work on.
!  (in, opt) full_precision - A flag indicating how carefully to iterate.  The
!                             default is .true. (more accurate).
!  (out,opt) uh_3d - Volume flux through zonal faces = u*h*dy, H m2 s-1.

  real, dimension(SZIB_(G),SZK_(G)) :: &
    uh_aux, &  ! An auxiliary zonal volume flux, in H m s-1.
    duhdu      ! Partial derivative of uh with u, in H m.
  real, dimension(SZIB_(G)) :: &
    uh_err, &  ! Difference between uhbt and the summed uh, in H m2 s-1.
    uh_err_best, & ! The smallest value of uh_err found so far, in H m2 s-1.
    u_new, &   ! The velocity with the correction added, in m s-1.
    duhdu_tot,&! Summed partial derivative of uh with u, in H m.
    du_min, &  ! Min/max limits on du correction based on CFL limits
    du_max     ! and previous iterations, in m s-1.
  real :: du_prev ! The previous value of du, in m s-1.
  real :: ddu    ! The change in du from the previous iteration, in m s-1.
  real :: tol_eta ! The tolerance for the current iteration, in m.
  real :: tol_vel ! The tolerance for velocity in the current iteration, m s-1.
  integer :: i, k, nz, itt, max_itts = 20
  logical :: full_prec, domore, do_i(SZIB_(G))

  nz = G%ke
  full_prec = .true. ; if (present(full_precision)) full_prec = full_precision

  uh_aux(:,:) = 0.0 ; duhdu(:,:) = 0.0

  if (present(uh_3d)) then ; do k=1,nz ; do I=ish-1,ieh
    uh_aux(i,k) = uh_3d(I,j,k)
  enddo ; enddo ; endif

  do I=ish-1,ieh
    du(I) = 0.0 ; do_i(I) = do_i_in(I)
    du_max(I) = du_max_CFL(I) ; du_min(I) = du_min_CFL(I)
    uh_err(I) = uh_tot_0(I) - uhbt(I) ; duhdu_tot(I) = duhdu_tot_0(I)
    uh_err_best(I) = abs(uh_err(I))
  enddo

  do itt=1,max_itts
    if (full_prec) then
      select case (itt)
        case (:1) ; tol_eta = 1e-6 * CS%tol_eta
        case (2)  ; tol_eta = 1e-4 * CS%tol_eta
        case (3)  ; tol_eta = 1e-2 * CS%tol_eta
        case default ; tol_eta = CS%tol_eta
      end select
    else
      tol_eta = CS%tol_eta_aux ; if (itt<=1) tol_eta = 1e-6 * CS%tol_eta_aux
    endif
    tol_vel = CS%tol_vel

    do I=ish-1,ieh
      if (uh_err(I) > 0.0) then ; du_max(I) = du(I)
      elseif (uh_err(I) < 0.0) then ; du_min(I) = du(I)
      else ; do_I(I) = .false. ; endif
    enddo
    domore = .false.
    do I=ish-1,ieh ; if (do_i(I)) then
      if ((dt*min(G%IareaT(i,j),G%IareaT(i+1,j))*abs(uh_err(I)) > tol_eta) .or.&
          (CS%better_iter .and. ((abs(uh_err(I)) > tol_vel * duhdu_tot(I)) .or.&
                                 (abs(uh_err(I)) > uh_err_best(I))) )) then
        !   Use Newton's method, provided it stays bounded.  Otherwise bisect
        ! the value with the appropriate bound.
        if (full_prec) then
          ddu = -uh_err(I) / duhdu_tot(I)
          du_prev = du(I)
          du(I) = du(I) + ddu
          if (abs(ddu) < 1.0e-15*abs(du(I))) then
            do_i(I) = .false. ! ddu is small enough to quit.
          elseif (ddu > 0.0) then
            if (du(I) >= du_max(I)) then
              du(I) = 0.5*(du_prev + du_max(I))
              if (du_max(I) - du_prev < 1.0e-15*abs(du(I))) do_i(I) = .false.
            endif
          else ! ddu < 0.0
            if (du(I) <= du_min(I)) then
              du(I) = 0.5*(du_prev + du_min(I))
              if (du_prev - du_min(I) < 1.0e-15*abs(du(I))) do_i(I) = .false.
            endif
          endif
        else
          !   Use Newton's method, provided it stays bounded, just like above.
          du(I) = du(I) - uh_err(I) / duhdu_tot(I)
          if ((du(I) >= du_max(I)) .or. (du(I) <= du_min(I))) &
            du(I) = 0.5*(du_max(I) + du_min(I))
        endif
        if (do_i(I)) domore = .true.
      else
        do_i(I) = .false.
      endif
    endif ; enddo
    if (.not.domore) exit

    if ((itt < max_itts) .or. present(uh_3d)) then ; do k=1,nz
      do I=ish-1,ieh ; u_new(I) = u(I,j,k) + du(I) * visc_rem(I,k) ; enddo
      call zonal_flux_layer(u_new, h_in(:,j,k), hL(:,j,k), hR(:,j,k), &
                            uh_aux(:,k), duhdu(:,k), visc_rem(:,k), &
                            dt, G, j, ish, ieh, do_i, CS%vol_CFL)
    enddo ; endif

    if (itt < max_itts) then
      do I=ish-1,ieh
        uh_err(I) = -uhbt(I) ; duhdu_tot(I) = 0.0
      enddo
      do k=1,nz ; do I=ish-1,ieh
        uh_err(I) = uh_err(I) + uh_aux(I,k)
        duhdu_tot(I) = duhdu_tot(I) + duhdu(I,k)
      enddo ; enddo
      do I=ish-1,ieh
        uh_err_best(I) = min(uh_err_best(I), abs(uh_err(I)))
      enddo
    endif
  enddo ! itt-loop
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.

  if (present(uh_3d)) then ; do k=1,nz ; do I=ish-1,ieh
    uh_3d(I,j,k) = uh_aux(I,k)
  enddo ; enddo ; endif

end subroutine zonal_flux_adjust

subroutine set_zonal_BT_cont(u, h_in, hL, hR, BT_cont, uh_tot_0, duhdu_tot_0, &
                             du_max_CFL, du_min_CFL, dt, G, CS, visc_rem, &
                             visc_rem_max, j, ish, ieh, do_i)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in) :: u
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in) :: h_in, hL, hR
  type(BT_cont_type),                  intent(inout) :: BT_cont
  real, dimension(NIMEMB_),            intent(in)    :: uh_tot_0, duhdu_tot_0
  real, dimension(NIMEMB_),            intent(in)    :: du_max_CFL, du_min_CFL
  real,                                intent(in)    :: dt
  type(ocean_grid_type),               intent(inout) :: G
  type(continuity_PPM_CS),             pointer       :: CS
  real, dimension(NIMEMB_,NKMEM_),     intent(in)    :: visc_rem
  real, dimension(NIMEMB_),            intent(in)    :: visc_rem_max
  integer,                             intent(in)    :: j, ish, ieh
  logical, dimension(NIMEMB_),         intent(in)    :: do_i
!   This subroutine sets of a structure that describes the zonal barotropic
! volume or mass fluxes as a function of barotropic flow to agree closely with
! the sum of the layer's transports.
!
! Arguments: u - Zonal velocity, in m s-1.
!  (in)      h_in - Layer thickness used to calculate the fluxes, in H.
!  (in)      hL, hR - Left- and right- thicknesses in the reconstruction, in H.
!  (inout)   BT_cont - A structure with elements that describe the effective
!                      open face areas as a function of barotropic flow.
!  (in)      uh_tot_0 - The summed transport with 0 adjustment, in H m2 s-1.
!  (in)      duhdu_tot_0 - The partial derivative of du_err with du at 0
!                          adjustment, in H m.
!  (out)     du - The barotropic velocity adjustment, in m s-1.
!  (in)      du_max_CFL, du_min_CFL - The acceptable range of du, in m s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 continuity_PPM_init.
!  (in)      visc_rem - Both the fraction of the momentum originally in a
!                       layer that remains after a time-step of viscosity,
!                       and the fraction of a time-step's worth of a
!                       barotropic acceleration that a layer experiences
!                       after viscosity is applied.  Nondimensional between
!                       0 (at the bottom) and 1 (far above the bottom).
!  (in)      j, ish, ieh - The index range to work on.
!  (in)      do_I - A logical flag indiciating which I values to work on.
  real, dimension(SZIB_(G)) :: &
    du0, &        ! The barotropic velocity increment that gives 0 transport, m s-1.
    duL, duR, &   ! The barotropic velocity increments that give the westerly
                  ! (duL) and easterly (duR) test velocities.
    zeros, &      ! An array of full of 0's.
    du_CFL, &     ! The velocity increment that corresponds to CFL_min, in m s-1.
    u_L, u_R, &   ! The westerly (u_L), easterly (u_R), and zero-barotropic
    u_0, &        ! transport (u_0) layer test velocities, in m s-1.
    FA_marg_L, &  ! The effective layer marginal face areas with the westerly
    FA_marg_R, &  ! (_L), easterly (_R), and zero-barotropic (_0) test
    FA_marg_0, &  ! velocities, in H m.
    uh_L, uh_R, & ! The layer transports with the westerly (_L), easterly (_R),
    uh_0, &       ! and zero-barotropic (_0) test velocities, in H m2 s-1.
    FAmt_L, FAmt_R, & ! The summed effective marginal face areas for the 3
    FAmt_0, &     ! test velocities, in H m.
    uhtot_L, &    ! The summed transport with the westerly (uhtot_L) and
    uhtot_R       ! and easterly (uhtot_R) test velocities, in H m2 s-1.
  real :: FA_0    ! The effective face area with 0 barotropic transport, in m H.
  real :: FA_avg  ! The average effective face area, in m H, nominally given by
                  ! the realized transport divided by the barotropic velocity.
  real :: visc_rem_lim ! The larger of visc_rem and min_visc_rem, ND.  This
                       ! limiting is necessary to keep the inverse of visc_rem
                       ! from leading to large CFL numbers.
  real :: min_visc_rem ! The smallest permitted value for visc_rem that is used
                       ! in finding the barotropic velocity that changes the
                       ! flow direction.  This is necessary to keep the inverse
                       ! of visc_rem from leading to large CFL numbers.
  real :: CFL_min ! A minimal increment in the CFL to try to ensure that the
                  ! flow is truly upwind, ND.
  real :: Idt     ! The inverse of the time step, in s-1.
  logical :: domore
  integer :: i, k, nz

  nz = G%ke ; Idt = 1.0/dt
  min_visc_rem = 0.1 ; CFL_min = 1e-6

 ! Diagnose the zero-transport correction, du0.
  do I=ish-1,ieh ; zeros(I) = 0.0 ; enddo
  call zonal_flux_adjust(u, h_in, hL, hR, zeros, uh_tot_0, &
                         duhdu_tot_0, du0, du_max_CFL, du_min_CFL, dt, G, &
                         CS, visc_rem, j, ish, ieh, do_I, .true.)

  ! Determine the westerly- and easterly- fluxes.  Choose a sufficiently
  ! negative velocity correction for the easterly-flux, and a sufficiently
  ! positive correction for the westerly-flux.
  domore = .false.
  do I=ish-1,ieh
    if (do_i(I)) domore = .true.
    du_CFL(I) = (CFL_min * Idt) * G%dxCu(I,j)
    duR(I) = min(0.0,du0(I) - du_CFL(I))
    duL(I) = max(0.0,du0(I) + du_CFL(I))
    FAmt_L(I) = 0.0 ; FAmt_R(I) = 0.0 ; FAmt_0(I) = 0.0
    uhtot_L(I) = 0.0 ; uhtot_R(I) = 0.0
  enddo

  if (.not.domore) then
    do k=1,nz ; do I=ish-1,ieh
      BT_cont%FA_u_W0(I,j) = 0.0 ; BT_cont%FA_u_WW(I,j) = 0.0
      BT_cont%uBT_WW(I,j) = 0.0
      BT_cont%FA_u_E0(I,j) = 0.0 ; BT_cont%FA_u_EE(I,j) = 0.0
      BT_cont%uBT_EE(I,j) = 0.0
    enddo ; enddo
    return
  endif

  do k=1,nz ; do I=ish-1,ieh ; if (do_i(I)) then
    visc_rem_lim = max(visc_rem(I,k), min_visc_rem*visc_rem_max(I))
    if (u(I,j,k) + duR(I)*visc_rem_lim > -du_CFL(I)*visc_rem(I,k)) &
      duR(I) = -(u(I,j,k) + du_CFL(I)*visc_rem(I,k)) / visc_rem_lim
    if (u(I,j,k) + duL(I)*visc_rem_lim < du_CFL(I)*visc_rem(I,k)) &
      duL(I) = -(u(I,j,k) - du_CFL(I)*visc_rem(I,k)) / visc_rem_lim
  endif ; enddo ; enddo

  do k=1,nz
    do I=ish-1,ieh ; if (do_i(I)) then
      u_L(I) = u(I,j,k) + duL(I) * visc_rem(I,k)
      u_R(I) = u(I,j,k) + duR(I) * visc_rem(I,k)
      u_0(I) = u(I,j,k) + du0(I) * visc_rem(I,k)
    endif ; enddo
    call zonal_flux_layer(u_0, h_in(:,j,k), hL(:,j,k), hR(:,j,k), uh_0, &
                          FA_marg_0, visc_rem(:,k), dt, G, j, ish, ieh, do_i, &
                          CS%vol_CFL)
    call zonal_flux_layer(u_L, h_in(:,j,k), hL(:,j,k), hR(:,j,k), uh_L, &
                          FA_marg_L, visc_rem(:,k), dt, G, j, ish, ieh, do_i, &
                          CS%vol_CFL)
    call zonal_flux_layer(u_R, h_in(:,j,k), hL(:,j,k), hR(:,j,k), uh_R, &
                          FA_marg_R, visc_rem(:,k), dt, G, j, ish, ieh, do_i, &
                          CS%vol_CFL)
    do I=ish-1,ieh ; if (do_i(I)) then
      FAmt_0(I) = FAmt_0(I) + FA_marg_0(I)
      FAmt_L(I) = FAmt_L(I) + FA_marg_L(I)
      FAmt_R(I) = FAmt_R(I) + FA_marg_R(I)
      uhtot_L(I) = uhtot_L(I) + uh_L(I)
      uhtot_R(I) = uhtot_R(I) + uh_R(I)
    endif ; enddo
  enddo
  do I=ish-1,ieh ; if (do_i(I)) then
    FA_0 = FAmt_0(I) ; FA_avg = FAmt_0(I)
    if ((duL(I) - du0(I)) /= 0.0) &
      FA_avg = uhtot_L(I) / (duL(I) - du0(I))
    if (FA_avg > max(FA_0, FAmt_L(I))) then ; FA_avg = max(FA_0, FAmt_L(I))
    elseif (FA_avg < min(FA_0, FAmt_L(I))) then ; FA_0 = FA_avg ; endif

    BT_cont%FA_u_W0(I,j) = FA_0 ; BT_cont%FA_u_WW(I,j) = FAmt_L(I)
    if (abs(FA_0-FAmt_L(I)) <= 1e-12*FA_0) then ; BT_cont%uBT_WW(I,j) = 0.0 ; else
      BT_cont%uBT_WW(I,j) = (1.5 * (duL(I) - du0(I))) * &
                            ((FAmt_L(I) - FA_avg) / (FAmt_L(I) - FA_0))
    endif

    FA_0 = FAmt_0(I) ; FA_avg = FAmt_0(I)
    if ((duR(I) - du0(I)) /= 0.0) &
      FA_avg = uhtot_R(I) / (duR(I) - du0(I))
    if (FA_avg > max(FA_0, FAmt_R(I))) then ; FA_avg = max(FA_0, FAmt_R(I))
    elseif (FA_avg < min(FA_0, FAmt_R(I))) then ; FA_0 = FA_avg ; endif

    BT_cont%FA_u_E0(I,j) = FA_0 ; BT_cont%FA_u_EE(I,j) = FAmt_R(I)
    if (abs(FAmt_R(I) - FA_0) <= 1e-12*FA_0) then ; BT_cont%uBT_EE(I,j) = 0.0 ; else
      BT_cont%uBT_EE(I,j) = (1.5 * (duR(I) - du0(I))) * &
                            ((FAmt_R(I) - FA_avg) / (FAmt_R(I) - FA_0))
    endif
  else
    BT_cont%FA_u_W0(I,j) = 0.0 ; BT_cont%FA_u_WW(I,j) = 0.0
    BT_cont%uBT_WW(I,j) = 0.0
    BT_cont%FA_u_E0(I,j) = 0.0 ; BT_cont%FA_u_EE(I,j) = 0.0
    BT_cont%uBT_EE(I,j) = 0.0
  endif ; enddo

end subroutine set_zonal_BT_cont

subroutine meridional_mass_flux(v, h_in, vh, dt, G, GV, CS, LB, vhbt, OBC, &
                                visc_rem_v, v_cor, vhbt_aux, v_cor_aux, BT_cont)
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in)    :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: h_in
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(out)   :: vh
  real,                                   intent(in)    :: dt
  type(ocean_grid_type),                  intent(inout) :: G
  type(verticalGrid_type),                intent(in)    :: GV
  type(continuity_PPM_CS),                pointer       :: CS
  type(loop_bounds_type),                 intent(in)    :: LB
  type(ocean_OBC_type),                   pointer,     optional :: OBC
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in),  optional :: visc_rem_v
  real, dimension(NIMEM_,NJMEMB_),        intent(in),  optional :: vhbt, vhbt_aux
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(out), optional :: v_cor, v_cor_aux
  type(BT_cont_type),                     pointer,     optional :: BT_cont
!   This subroutine calculates the mass or volume fluxes through the meridional
! faces, and other related quantities.
! Arguments: v - Meridional velocity, in m s-1.
!  (in)      h_in - Layer thickness used to calculate the fluxes, in H.
!  (out)     vh - Volume flux through meridional faces = v*h*dy, H m2 s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 continuity_PPM_init.
!  (in)      LB - A structure with the active loop bounds.
!  (in, opt) vhbt - The summed volume flux through meridional faces, H m2 s-1.
!  (in, opt) OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in, opt) visc_rem_v - Both the fraction of the momentum originally in a
!                         layer that remains after a time-step of viscosity,
!                         and the fraction of a time-step's worth of a
!                         barotropic acceleration that a layer experiences
!                         after viscosity is applied.  Nondimensional between
!                         0 (at the bottom) and 1 (far above the bottom).
!  (out, opt) v_cor - The meridional velocities (v with a barotropic corection)
!                     that give vhbt as the depth-integrated transport, m s-1.
!  (in, opt) vhbt_aux - A second set of summed volume fluxes through meridional
!                       faces, in H m2 s-1.
!  (out, opt) v_cor_aux - The meridional velocities (v with a barotropic
!                         correction) that give vhbt_aux as the depth-integrated
!                         transports, in m s-1.
!  (out, opt) BT_cont - A structure with elements that describe the effective
!                       open face areas as a function of barotropic flow.
  real, dimension(SZI_(G),SZK_(G)) :: &
    dvhdv      ! Partial derivative of vh with v, in m2.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    hl, hr      ! Left and right face thicknesses, in m.
  real, dimension(SZI_(G)) :: &
    dv, &      ! Corrective barotropic change in the velocity, in m s-1.
    dv_min_CFL, & ! Min/max limits on dv correction
    dv_max_CFL, & ! to avoid CFL violations
    dvhdv_tot_0, & ! Summed partial derivative of vh with v, in H m.
    vh_tot_0, &   ! Summed transport with no barotropic correction in H m2 s-1.
    visc_rem_max  ! The column maximum of visc_rem.
  logical, dimension(SZI_(G)) :: do_i
  real, dimension(SZI_(G),SZK_(G)) :: &
    visc_rem      ! A 2-D copy of visc_rem_v or an array of 1's.
  real :: I_vrm   ! 1.0 / visc_rem_max, nondim.
  real :: CFL_dt  ! The maximum CFL ratio of the adjusted velocities divided by
                  ! the time step, in s-1.
  real :: I_dt    ! 1.0 / dt, in s-1.
  real :: dv_lim  ! The velocity change that give a relative CFL of 1, in m s-1.
  real :: dy_N, dy_S ! Effective y-grid spacings to the north and south, in m.
  integer :: i, j, k, ish, ieh, jsh, jeh, nz
  logical :: do_aux, apply_OBC_v, use_visc_rem, set_BT_cont, any_simple_OBC

  do_aux = (present(vhbt_aux) .and. present(v_cor_aux))
  use_visc_rem = present(visc_rem_v)
  apply_OBC_v = .false. ; set_BT_cont = .false.
  if (present(BT_cont)) set_BT_cont = (associated(BT_cont))
  if (present(OBC)) then ; if (associated(OBC)) then
    apply_OBC_v = OBC%apply_OBC_v
  endif ; endif
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = G%ke

  CFL_dt = CS%CFL_limit_adjust / dt
  I_dt = 1.0 / dt
  if (CS%aggress_adjust) CFL_dt = I_dt

  call cpu_clock_begin(id_clock_update)
!$OMP parallel do default(none) shared(nz,ish,ieh,jsh,jeh,h_in,hl,hr,G,GV,LB,CS,visc_rem)
  do k=1,nz
    ! This sets hl and hr.
    if (CS%upwind_1st) then
      do j=jsh-1,jeh+1 ; do i=ish,ieh
        hl(i,j,k) = h_in(i,j,k) ; hr(i,j,k) = h_in(i,j,k)
      enddo ; enddo
    else
      call PPM_reconstruction_y(h_in(:,:,k), hl(:,:,k), hr(:,:,k), G, LB, &
                                2.0*GV%Angstrom, CS%monotonic, simple_2nd=CS%simple_2nd)
    endif
    do i=ish,ieh ; visc_rem(i,k) = 1.0 ; enddo
  enddo
  call cpu_clock_end(id_clock_update)

  call cpu_clock_begin(id_clock_correct)
!$OMP parallel do default(none) shared(ish,ieh,jsh,jeh,nz,v,h_in,hL,hR,vh,use_visc_rem, &
!$OMP                                  visc_rem_v,dt,G,GV,CS,apply_OBC_v,OBC,vhbt,do_aux, &
!$OMP                                  set_BT_cont,CFL_dt,I_dt,v_cor,vhbt_aux,          &
!$OMP                                  v_cor_aux,BT_cont )                              &
!$OMP                          private(do_i,dvhdv,dv,dv_max_CFL,dv_min_CFL,vh_tot_0,    &
!$OMP                                  dvhdv_tot_0,visc_rem_max,I_vrm,dv_lim,dy_N,      &
!$OMP                                  dy_S,any_simple_OBC ) &
!$OMP                     firstprivate(visc_rem)
  do J=jsh-1,jeh
    do i=ish,ieh ; do_i(i) = .true. ; visc_rem_max(I) = 0.0 ; enddo
    ! This sets vh and dvhdv.
    do k=1,nz
      if (use_visc_rem) then ; do i=ish,ieh
        visc_rem(i,k) = visc_rem_v(i,J,k)
        visc_rem_max(i) = max(visc_rem_max(i), visc_rem(i,k))
      enddo ; endif
      call merid_flux_layer(v(:,J,k), h_in(:,:,k), hL(:,:,k), hR(:,:,k), &
                            vh(:,J,k), dvhdv(:,k), visc_rem(:,k), &
                            dt, G, J, ish, ieh, do_i, CS%vol_CFL)
      if (apply_OBC_v) then ; do i=ish,ieh
        if (OBC%OBC_mask_v(i,J) .and. (OBC%OBC_kind_v(i,J) == OBC_SIMPLE)) &
          vh(i,J,k) = OBC%vh(i,J,k)
      enddo ; endif
    enddo ! k-loop
    if ((.not.use_visc_rem) .or. (.not.CS%use_visc_rem_max)) then ; do i=ish,ieh
      visc_rem_max(i) = 1.0
    enddo ; endif

    if (present(vhbt) .or. do_aux .or. set_BT_cont) then
      !   Set limits on dv that will keep the CFL number between -1 and 1.
      ! This should be adequate to keep the root bracketed in all cases.
      do i=ish,ieh
        I_vrm = 0.0
        if (visc_rem_max(i) > 0.0) I_vrm = 1.0 / visc_rem_max(i)
        if (CS%vol_CFL) then
          dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(i,J), 1000.0*G%dyT(i,j))
          dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(i,J), 1000.0*G%dyT(i,j+1))
        else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
        dv_max_CFL(i) = 2.0 * (CFL_dt * dy_S) * I_vrm
        dv_min_CFL(i) = -2.0 * (CFL_dt * dy_N) * I_vrm
        vh_tot_0(i) = 0.0 ; dvhdv_tot_0(i) = 0.0
      enddo
      do k=1,nz ; do i=ish,ieh
        dvhdv_tot_0(i) = dvhdv_tot_0(i) + dvhdv(i,k)
        vh_tot_0(i) = vh_tot_0(i) + vh(i,J,k)
      enddo ; enddo

      if (use_visc_rem) then
        if (CS%aggress_adjust) then
          do k=1,nz ; do i=ish,ieh
            if (CS%vol_CFL) then
              dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(I,j), 1000.0*G%dyT(i,j))
              dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(I,j), 1000.0*G%dyT(i,j+1))
            else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
            dv_lim = 0.499*((dy_S*I_dt - v(i,J,k)) + MIN(0.0,v(i,J-1,k)))
            if (dv_max_CFL(i) * visc_rem(i,k) > dv_lim) &
              dv_max_CFL(i) = dv_lim / visc_rem(i,k)

            dv_lim = 0.499*((-dy_N*CFL_dt - v(i,J,k)) + MAX(0.0,v(i,J+1,k)))
            if (dv_min_CFL(i) * visc_rem(i,k) < dv_lim) &
              dv_min_CFL(i) = dv_lim / visc_rem(i,k)
          enddo ; enddo
        else
          do k=1,nz ; do i=ish,ieh
            if (CS%vol_CFL) then
              dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(I,j), 1000.0*G%dyT(i,j))
              dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(I,j), 1000.0*G%dyT(i,j+1))
            else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
            if (dv_max_CFL(i) * visc_rem(i,k) > dy_S*CFL_dt - v(i,J,k)) &
              dv_max_CFL(i) = (dy_S*CFL_dt - v(i,J,k)) / visc_rem(i,k)
            if (dv_min_CFL(i) * visc_rem(i,k) < -dy_N*CFL_dt - v(i,J,k)) &
              dv_min_CFL(i) = -(dy_N*CFL_dt + v(i,J,k)) / visc_rem(i,k)
          enddo ; enddo
        endif
      else
        if (CS%aggress_adjust) then
          do k=1,nz ; do i=ish,ieh
            if (CS%vol_CFL) then
              dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(I,j), 1000.0*G%dyT(i,j))
              dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(I,j), 1000.0*G%dyT(i,j+1))
            else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
            dv_max_CFL(i) = min(dv_max_CFL(i), 0.499 * &
                        ((dy_S*I_dt - v(i,J,k)) + MIN(0.0,v(i,J-1,k))) )
            dv_min_CFL(i) = max(dv_min_CFL(i), 0.499 * &
                        ((-dy_N*I_dt - v(i,J,k)) + MAX(0.0,v(i,J+1,k))) )
          enddo ; enddo
        else
          do k=1,nz ; do i=ish,ieh
            if (CS%vol_CFL) then
              dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(I,j), 1000.0*G%dyT(i,j))
              dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(I,j), 1000.0*G%dyT(i,j+1))
            else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
            dv_max_CFL(i) = min(dv_max_CFL(i), dy_S*CFL_dt - v(i,J,k))
            dv_min_CFL(i) = max(dv_min_CFL(i), -(dy_N*CFL_dt + v(i,J,k)))
          enddo ; enddo
        endif
      endif
      do i=ish,ieh
        dv_max_CFL(i) = max(dv_max_CFL(i),0.0)
        dv_min_CFL(i) = min(dv_min_CFL(i),0.0)
      enddo

      ! Up to this point, everything is shared between vhbt and vhbt_aux.

      any_simple_OBC = .false.
      if (present(vhbt) .or. do_aux .or. set_BT_cont) then
        if (.not.apply_OBC_v) then ; do i=ish,ieh
          do_i(i) = .true.
        enddo ; else ; do i=ish,ieh
          do_i(i) = .not.(OBC%OBC_mask_v(i,J) .and. &
                          (OBC%OBC_kind_v(i,J) == OBC_SIMPLE))
          if (.not.do_i(I)) any_simple_OBC = .true.
        enddo ; endif
      endif

      if (present(vhbt)) then
        call meridional_flux_adjust(v, h_in, hL, hR, vhbt(:,J), vh_tot_0, &
                               dvhdv_tot_0, dv, dv_max_CFL, dv_min_CFL, dt, G, &
                               CS, visc_rem, j, ish, ieh, do_i, .true., vh)

        if (present(v_cor)) then ; do k=1,nz
          do i=ish,ieh ; v_cor(i,J,k) = v(i,J,k) + dv(i) * visc_rem(i,k) ; enddo
          if (apply_OBC_v) then ; do i=ish,ieh
            if (OBC%OBC_mask_v(i,J) .and. (OBC%OBC_kind_v(i,J) == OBC_SIMPLE)) &
              v_cor(i,J,k) = OBC%v(i,J,k)
          enddo ; endif
        enddo ; endif ! v-corrected
      endif

      if (do_aux) then
        call meridional_flux_adjust(v, h_in, hL, hR, vhbt_aux(:,J), vh_tot_0, &
                               dvhdv_tot_0, dv, dv_max_CFL, dv_min_CFL, dt, G, &
                               CS, visc_rem, j, ish, ieh, do_i, .false.)

        do k=1,nz
          do i=ish,ieh ; v_cor_aux(i,J,k) = v(i,J,k) + dv(i) * visc_rem(i,k) ; enddo
          if (apply_OBC_v) then ; do i=ish,ieh
            if (OBC%OBC_mask_v(i,J) .and. (OBC%OBC_kind_v(i,J) == OBC_SIMPLE)) &
              v_cor_aux(i,J,k) = OBC%v(i,J,k)
          enddo ; endif
        enddo
      endif ! do_aux

      if (set_BT_cont) then
        call set_merid_BT_cont(v, h_in, hL, hR, BT_cont, vh_tot_0, dvhdv_tot_0,&
                               dv_max_CFL, dv_min_CFL, dt, G, CS, visc_rem, &
                               visc_rem_max, J, ish, ieh, do_i)
        if (any_simple_OBC) then
          do i=ish,ieh
            do_i(i) = (OBC%OBC_mask_v(i,J) .and. &
                       (OBC%OBC_kind_v(i,J) == OBC_SIMPLE))
            if (do_i(i)) BT_cont%Fa_v_S0(i,J) = GV%H_subroundoff*G%dx_Cv(I,j)
          enddo
          do k=1,nz ; do i=ish,ieh ; if (do_i(i)) then
            if (abs(OBC%v(i,J,k)) > 0.0) &
              BT_cont%Fa_v_S0(i,J) = BT_cont%Fa_v_S0(i,J) + &
                   OBC%vh(i,J,k) / OBC%v(i,J,k)
          endif ; enddo ; enddo
          do i=ish,ieh ; if (do_i(i)) then
            BT_cont%Fa_v_N0(i,J) = BT_cont%Fa_v_S0(i,J)
            BT_cont%Fa_v_SS(i,J) = BT_cont%Fa_v_S0(i,J)
            BT_cont%Fa_v_NN(i,J) = BT_cont%Fa_v_S0(i,J)
          endif ; enddo
        endif
      endif ! set_BT_cont

    endif ! present(vhbt) or do_aux or set_BT_cont
  enddo ! j-loop
  call cpu_clock_end(id_clock_correct)

  if (set_BT_cont) then ; if (associated(BT_cont%h_v)) then
    if (present(v_cor)) then
      call merid_face_thickness(v_cor, h_in, hL, hR, BT_cont%h_v, dt, G, LB, &
                                CS%vol_CFL, CS%marginal_faces, visc_rem_v)
    else
      call merid_face_thickness(v, h_in, hL, hR, BT_cont%h_v, dt, G, LB, &
                                CS%vol_CFL, CS%marginal_faces, visc_rem_v)
    endif
  endif ; endif

end subroutine meridional_mass_flux

subroutine merid_flux_layer(v, h, hL, hR, vh, dvhdv, visc_rem, dt, G, J, &
                            ish, ieh, do_i, vol_CFL)
  real, dimension(NIMEM_),        intent(in)    :: v, visc_rem
  real, dimension(NIMEM_,NJMEM_), intent(in)    :: h, hL, hR
  real, dimension(NIMEM_),        intent(inout) :: vh, dvhdv
  real,                           intent(in)    :: dt
  type(ocean_grid_type),          intent(inout) :: G
  integer,                        intent(in)    :: J, ish, ieh
  logical, dimension(NIMEM_),     intent(in)    :: do_i
  logical,                        intent(in)    :: vol_CFL
!   This subroutines evaluates the meridional mass or volume fluxes in a layer.
!
! Arguments: v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness used to calculate the fluxes, in H.
!  (in)      hL, hR - Left- and right- thicknesses in the reconstruction, in H.
!  (out)     vh - The meridional mass or volume transport, in H m2 s-1.
!  (out)     dvhdv - The partial derivative of vh with v, in H m.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      visc_rem - Both the fraction of the momentum originally in a
!                       layer that remains after a time-step of viscosity,
!                       and the fraction of a time-step's worth of a
!                       barotropic acceleration that a layer experiences
!                       after viscosity is applied.  Nondimensional between
!                       0 (at the bottom) and 1 (far above the bottom).
!  (in)      J, ish, ieh - The index range to work on.
!  (in)      do_i - A logical flag indiciating which i values to work on.
!  (in)      vol_CFL - If true, rescale the ratio of face areas to the cell
!                       areas when estimating the CFL number.
  real :: CFL ! The CFL number based on the local velocity and grid spacing, ND.
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h_in.
  real :: h_marg ! The marginal thickness of a flux, in m.
  integer :: i

  do i=ish,ieh ; if (do_i(i)) then
    if (v(i) > 0.0) then
      if (vol_CFL) then ; CFL = (v(i) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
      else ; CFL = v(i) * dt * G%IdyT(i,j) ; endif
      curv_3 = hL(i,j) + hR(i,j) - 2.0*h(i,j)
      vh(i) = G%dx_Cv(i,J) * v(i) * ( hR(i,j) + CFL * &
          (0.5*(hL(i,j) - hR(i,j)) + curv_3*(CFL - 1.5)) )
      h_marg = hR(i,j) + CFL * ((hL(i,j) - hR(i,j)) + &
                                  3.0*curv_3*(CFL - 1.0))
    elseif (v(i) < 0.0) then
      if (vol_CFL) then ; CFL = (-v(i) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
      else ; CFL = -v(i) * dt * G%IdyT(i,j+1) ; endif
      curv_3 = hL(i,j+1) + hR(i,j+1) - 2.0*h(i,j+1)
      vh(i) = G%dx_Cv(i,J) * v(i) * ( hL(i,j+1) + CFL * &
          (0.5*(hR(i,j+1)-hL(i,j+1)) + curv_3*(CFL - 1.5)) )
      h_marg = hL(i,j+1) + CFL * ((hR(i,j+1)-hL(i,j+1)) + &
                                    3.0*curv_3*(CFL - 1.0))
    else
      vh(i) = 0.0
      h_marg = 0.5 * (hl(i,j+1) + hr(i,j))
    endif
    dvhdv(i) = G%dx_Cv(i,J) * h_marg * visc_rem(i)
  endif ; enddo

end subroutine merid_flux_layer

subroutine merid_face_thickness(v, h, hL, hR, h_v, dt, G, LB, vol_CFL, &
                                marginal, visc_rem_v)
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in)    :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: h, hL, hR
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(inout) :: h_v
  real,                                   intent(in)    :: dt
  type(ocean_grid_type),                  intent(inout) :: G
  type(loop_bounds_type),                 intent(in)    :: LB
  logical,                                intent(in)    :: vol_CFL
  logical,                                intent(in)    :: marginal
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in),  optional :: visc_rem_v
!   This subroutines sets the effective interface thickness at each meridional
! velocity point.
!
! Arguments: v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness used to calculate the fluxes, in H.
!  (in)      hL, hR - Left- and right- thicknesses in the reconstruction, in H.
!  (out)     h_v - The thickness at meridional faces, in H.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      LB - A structure with the active thickness loop bounds.
!  (in)      vol_CFL - If true, rescale the ratio of face areas to the cell
!                       areas when estimating the CFL number.
!  (in)      marginal - If true, report the marginal face thicknesses; otherwise
!                       report transport-averaged thicknesses.
!  (in, opt) visc_rem_v - Both the fraction of the momentum originally in a
!                         layer that remains after a time-step of viscosity,
!                         and the fraction of a time-step's worth of a
!                         barotropic acceleration that a layer experiences
!                         after viscosity is applied.  Nondimensional between
!                         0 (at the bottom) and 1 (far above the bottom).
  real :: CFL ! The CFL number based on the local velocity and grid spacing, ND.
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h_in.
  real :: h_avg  ! The average thickness of a flux, in H.
  real :: h_marg ! The marginal thickness of a flux, in H.
  integer :: i, j, k, ish, ieh, jsh, jeh, nz
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = G%ke

!$OMP parallel default(none) shared(ish,ieh,jsh,jeh,nz,v,vol_CFL,dt,G, &
!$OMP                               hL,hR,h,h_v,visc_rem_v,marginal) &
!$OMP                       private(CFL,curv_3,h_marg,h_avg)
!$OMP do
  do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
    if (v(i,J,k) > 0.0) then
      if (vol_CFL) then ; CFL = (v(i,J,k) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
      else ; CFL = v(i,J,k) * dt * G%IdyT(i,j) ; endif
      curv_3 = hL(i,j,k) + hR(i,j,k) - 2.0*h(i,j,k)
      h_avg = hR(i,j,k) + CFL * (0.5*(hL(i,j,k) - hR(i,j,k)) + curv_3*(CFL - 1.5))
      h_marg = hR(i,j,k) + CFL * ((hL(i,j,k) - hR(i,j,k)) + &
                                3.0*curv_3*(CFL - 1.0))
    elseif (v(i,J,k) < 0.0) then
      if (vol_CFL) then ; CFL = (-v(i,J,k)*dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
      else ; CFL = -v(i,J,k) * dt * G%IdyT(i,j+1) ; endif
      curv_3 = hL(i,j+1,k) + hR(i,j+1,k) - 2.0*h(i,j+1,k)
      h_avg = hL(i,j+1,k) + CFL * (0.5*(hR(i,j+1,k)-hL(i,j+1,k)) + curv_3*(CFL - 1.5))
      h_marg = hL(i,j+1,k) + CFL * ((hR(i,j+1,k)-hL(i,j+1,k)) + &
                                    3.0*curv_3*(CFL - 1.0))
    else
      h_avg = 0.5 * (hl(i,j+1,k) + hr(i,j,k))
      !   The choice to use the arithmetic mean here is somewhat arbitrariy, but
      ! it should be noted that hl(i+1,j,k) and hr(i,j,k) are usually the same.
      h_marg = 0.5 * (hl(i,j+1,k) + hr(i,j,k))
 !    h_marg = (2.0 * hl(i,j+1,k) * hr(i,j,k)) / &
 !             (hl(i,j+1,k) + hr(i,j,k) + GV%H_subroundoff)
    endif

    if (marginal) then ; h_v(i,J,k) = h_marg
    else ; h_v(i,J,k) = h_avg ; endif
  enddo ; enddo ; enddo

  if (present(visc_rem_v)) then
!$OMP do
    do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
      h_v(i,J,k) = h_v(i,J,k) * visc_rem_v(i,J,k)
    enddo ; enddo ; enddo
  endif
!$OMP end parallel

end subroutine merid_face_thickness

subroutine meridional_flux_adjust(v, h_in, hL, hR, vhbt, vh_tot_0, dvhdv_tot_0, &
                             dv, dv_max_CFL, dv_min_CFL, dt, G, CS, visc_rem, &
                             j, ish, ieh, do_i_in, full_precision, vh_3d)
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in) :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in) :: h_in, hL, hR
  real, dimension(NIMEM_,NKMEM_),         intent(in) :: visc_rem
  real, dimension(NIMEM_),             intent(in)    :: vhbt, dv_max_CFL, dv_min_CFL
  real, dimension(NIMEM_),             intent(in)    :: vh_tot_0, dvhdv_tot_0
  real, dimension(NIMEM_),             intent(out)   :: dv
  real,                                intent(in)    :: dt
  type(ocean_grid_type),               intent(inout) :: G
  type(continuity_PPM_CS),             pointer       :: CS
  integer,                             intent(in)    :: j, ish, ieh
  logical, dimension(NIMEM_),          intent(in)    :: do_i_in
  logical,                             intent(in),    optional :: full_precision
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(inout), optional :: vh_3d

!   This subroutine returns the barotropic velocity adjustment that gives the
! desired barotropic (layer-summed) transport.
! Arguments: v - Meridional velocity, in m s-1.
!  (in)      h_in - Layer thickness used to calculate the fluxes, in H.
!  (in)      hL, hR - Left- and right- thicknesses in the reconstruction, in H.
!  (in)      vhbt - Summed volume flux through meridional faces, in H m2 s-1.
!  (in)      vh_tot_0 - The summed transport with 0 adjustment, in H m2 s-1.
!  (in)      dvhdv_tot_0 - The partial derivative of dv_err with dv at 0
!                          adjustment, in H m.
!  (out)     dv - The barotropic velocity adjustment, in m s-1.
!  (in)      dv_max_CFL, du_min_CFL - The acceptable range of dv, in m s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 continuity_PPM_init.
!  (in)      visc_rem - Both the fraction of the momentum originally in a
!                       layer that remains after a time-step of viscosity,
!                       and the fraction of a time-step's worth of a
!                       barotropic acceleration that a layer experiences
!                       after viscosity is applied.  Nondimensional between
!                       0 (at the bottom) and 1 (far above the bottom).
!  (in)      j, ish, ieh - The index range to work on.
!  (in)      do_i_in - A logical flag indiciating which i values to work on.
!  (in, opt) full_precision - A flag indicating how carefully to iterate.  The
!                             default is .true. (more accurate).
!  (out,opt) v_3d - Volume flux through meridional faces = v*h*dx, H m2 s-1.

  real, dimension(SZI_(G),SZK_(G)) :: &
    vh_aux, &  ! An auxiliary meridonal volume flux, in H m s-1.
    dvhdv      ! Partial derivative of vh with v, in H m.
  real, dimension(SZI_(G)) :: &
    vh_err, &  ! Difference between vhbt and the summed vh, in H m2 s-1.
    vh_err_best, & ! The smallest value of vh_err found so far, in H m2 s-1.
    v_new, &   ! The velocity with the correction added, in m s-1.
    dvhdv_tot,&! Summed partial derivative of vh with u, in H m.
    dv_min, &  ! Min/max limits on dv correction based on CFL limits
    dv_max     ! and previous iterations, in m s-1.
  real :: dv_prev ! The previous value of dv, in m s-1.
  real :: ddv    ! The change in dv from the previous iteration, in m s-1.
  real :: tol_eta ! The tolerance for the current iteration, in m.
  real :: tol_vel ! The tolerance for velocity in the current iteration, m s-1.
  integer :: i, k, nz, itt, max_itts = 20
  logical :: full_prec, domore, do_i(SZI_(G))

  nz = G%ke
  full_prec = .true. ; if (present(full_precision)) full_prec = full_precision

  vh_aux(:,:) = 0.0 ; dvhdv(:,:) = 0.0

  if (present(vh_3d)) then ; do k=1,nz ; do i=ish,ieh
    vh_aux(i,k) = vh_3d(i,J,k)
  enddo ; enddo ; endif

  do i=ish,ieh
    dv(i) = 0.0 ; do_i(i) = do_i_in(i)
    dv_max(i) = dv_max_CFL(i) ; dv_min(i) = dv_min_CFL(i)
    vh_err(i) = vh_tot_0(i) - vhbt(i) ; dvhdv_tot(i) = dvhdv_tot_0(i)
    vh_err_best(i) = abs(vh_err(i))
  enddo

  do itt=1,max_itts
    if (full_prec) then
      select case (itt)
        case (:1) ; tol_eta = 1e-6 * CS%tol_eta
        case (2)  ; tol_eta = 1e-4 * CS%tol_eta
        case (3)  ; tol_eta = 1e-2 * CS%tol_eta
        case default ; tol_eta = CS%tol_eta
      end select
    else
      tol_eta = CS%tol_eta_aux ; if (itt<=1) tol_eta = 1e-6 * CS%tol_eta_aux
    endif
    tol_vel = CS%tol_vel

    do i=ish,ieh
      if (vh_err(i) > 0.0) then ; dv_max(i) = dv(i)
      elseif (vh_err(i) < 0.0) then ; dv_min(i) = dv(i)
      else ; do_i(i) = .false. ; endif
    enddo
    domore = .false.
    do i=ish,ieh ; if (do_i(i)) then
      if ((dt*min(G%IareaT(i,j),G%IareaT(i,j+1))*abs(vh_err(i)) > tol_eta) .or.&
          (CS%better_iter .and. ((abs(vh_err(i)) > tol_vel * dvhdv_tot(i)) .or.&
                                 (abs(vh_err(i)) > vh_err_best(i))) )) then
        !   Use Newton's method, provided it stays bounded.  Otherwise bisect
        ! the value with the appropriate bound.
        if (full_prec) then
          ddv = -vh_err(i) / dvhdv_tot(i)
          dv_prev = dv(i)
          dv(i) = dv(i) + ddv
          if (abs(ddv) < 1.0e-15*abs(dv(i))) then
            do_i(i) = .false. ! ddv is small enough to quit.
          elseif (ddv > 0.0) then
            if (dv(i) >= dv_max(i)) then
              dv(i) = 0.5*(dv_prev + dv_max(i))
              if (dv_max(i) - dv_prev < 1.0e-15*abs(dv(i))) do_i(i) = .false.
            endif
          else ! dvv(i) < 0.0
            if (dv(i) <= dv_min(i)) then
              dv(i) = 0.5*(dv_prev + dv_min(i))
              if (dv_prev - dv_min(i) < 1.0e-15*abs(dv(i))) do_i(i) = .false.
            endif
          endif
        else
          !   Use Newton's method, provided it stays bounded, just like above.
          dv(i) = dv(i) - vh_err(i) / dvhdv_tot(i)
          if ((dv(i) >= dv_max(i)) .or. (dv(i) <= dv_min(i))) &
            dv(i) = 0.5*(dv_max(i) + dv_min(i))
        endif
        if (do_i(i)) domore = .true.
      else
        do_i(i) = .false.
      endif
    endif ; enddo
    if (.not.domore) exit

    if ((itt < max_itts) .or. present(vh_3d)) then ; do k=1,nz
      do i=ish,ieh ; v_new(i) = v(i,J,k) + dv(i) * visc_rem(i,k) ; enddo
      call merid_flux_layer(v_new, h_in(:,:,k), hL(:,:,k), hR(:,:,k), &
                            vh_aux(:,k), dvhdv(:,k), visc_rem(:,k), &
                            dt, G, J, ish, ieh, do_i, CS%vol_CFL)
    enddo ; endif

    if (itt < max_itts) then
      do i=ish,ieh
        vh_err(i) = -vhbt(i) ; dvhdv_tot(i) = 0.0
      enddo
      do k=1,nz ; do i=ish,ieh
        vh_err(i) = vh_err(i) + vh_aux(i,k)
        dvhdv_tot(i) = dvhdv_tot(i) + dvhdv(i,k)
      enddo ; enddo
      do i=ish,ieh
        vh_err_best(i) = min(vh_err_best(i), abs(vh_err(i)))
      enddo
    endif
  enddo ! itt-loop
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.

  if (present(vh_3d)) then ; do k=1,nz ; do i=ish,ieh
    vh_3d(i,J,k) = vh_aux(i,k)
  enddo ; enddo ; endif

end subroutine meridional_flux_adjust

subroutine set_merid_BT_cont(v, h_in, hL, hR, BT_cont, vh_tot_0, dvhdv_tot_0, &
                             dv_max_CFL, dv_min_CFL, dt, G, CS, visc_rem, &
                             visc_rem_max, j, ish, ieh, do_i)
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in) :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in) :: h_in, hL, hR
  type(BT_cont_type),                  intent(inout) :: BT_cont
  real, dimension(NIMEM_),             intent(in)    :: vh_tot_0, dvhdv_tot_0
  real, dimension(NIMEM_),             intent(in)    :: dv_max_CFL, dv_min_CFL
  real,                                intent(in)    :: dt
  type(ocean_grid_type),               intent(inout) :: G
  type(continuity_PPM_CS),             pointer       :: CS
  real, dimension(NIMEM_,NKMEM_),      intent(in)    :: visc_rem
  real, dimension(NIMEM_),             intent(in)    :: visc_rem_max
  integer,                             intent(in)    :: j, ish, ieh
  logical, dimension(NIMEM_),          intent(in)    :: do_i
!   This subroutine sets of a structure that describes the meridional
! barotropic volume or mass fluxes as a function of barotropic flow to agree
! closely with the sum of the layer's transports.
!
! Arguments: v - Meridional velocity, in m s-1.
!  (in)      h_in - Layer thickness used to calculate the fluxes, in H.
!  (in)      hL, hR - Left- and right- thicknesses in the reconstruction, in H.
!  (inout)   BT_cont - A structure with elements that describe the effective
!                      open face areas as a function of barotropic flow.
!  (in)      vh_tot_0 - The summed transport with 0 adjustment, in H m2 s-1.
!  (in)      dvhdv_tot_0 - The partial derivative of dv_err with dv at 0
!                          adjustment, in H m.
!  (in)      dv_max_CFL, du_min_CFL - The acceptable range of dv, in m s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 continuity_PPM_init.
!  (in)      visc_rem - Both the fraction of the momentum originally in a
!                       layer that remains after a time-step of viscosity,
!                       and the fraction of a time-step's worth of a
!                       barotropic acceleration that a layer experiences
!                       after viscosity is applied.  Nondimensional between
!                       0 (at the bottom) and 1 (far above the bottom).
!  (in)      J, ish, ieh - The index range to work on.
!  (in)      do_i - A logical flag indiciating which i values to work on.

  real, dimension(SZI_(G)) :: &
    dv0, &        ! The barotropic velocity increment that gives 0 transport, m s-1.
    dvL, dvR, &   ! The barotropic velocity increments that give the southerly
                  ! (dvL) and northerly (dvR) test velocities.
    zeros, &      ! An array of full of 0's.
    dv_CFL, &     ! The velocity increment that corresponds to CFL_min, in m s-1.
    v_L, v_R, &   ! The southerly (v_L), northerly (v_R), and zero-barotropic
    v_0, &        ! transport (v_0) layer test velocities, in m s-1.
    FA_marg_L, &  ! The effective layer marginal face areas with the southerly
    FA_marg_R, &  ! (_L), northerly (_R), and zero-barotropic (_0) test
    FA_marg_0, &  ! velocities, in H m.
    vh_L, vh_R, & ! The layer transports with the southerly (_L), northerly (_R)
    vh_0, &       ! and zero-barotropic (_0) test velocities, in H m2 s-1.
    FAmt_L, FAmt_R, & ! The summed effective marginal face areas for the 3
    FAmt_0, &     ! test velocities, in H m.
    vhtot_L, &    ! The summed transport with the southerly (vhtot_L) and
    vhtot_R       ! and northerly (vhtot_R) test velocities, in H m2 s-1.
  real :: FA_0    ! The effective face area with 0 barotropic transport, in m H.
  real :: FA_avg  ! The average effective face area, in m H, nominally given by
                  ! the realized transport divided by the barotropic velocity.
  real :: visc_rem_lim ! The larger of visc_rem and min_visc_rem, ND.  This
                       ! limiting is necessary to keep the inverse of visc_rem
                       ! from leading to large CFL numbers.
  real :: min_visc_rem ! The smallest permitted value for visc_rem that is used
                       ! in finding the barotropic velocity that changes the
                       ! flow direction.  This is necessary to keep the inverse
                       ! of visc_rem from leading to large CFL numbers.
  real :: CFL_min ! A minimal increment in the CFL to try to ensure that the
                  ! flow is truly upwind, ND.
  real :: Idt     ! The inverse of the time step, in s-1.
  logical :: domore
  integer :: i, k, nz

  nz = G%ke ; Idt = 1.0/dt
  min_visc_rem = 0.1 ; CFL_min = 1e-6

 ! Diagnose the zero-transport correction, dv0.
  do i=ish,ieh ; zeros(i) = 0.0 ; enddo
  call meridional_flux_adjust(v, h_in, hL, hR, zeros, vh_tot_0, &
                         dvhdv_tot_0, dv0, dv_max_CFL, dv_min_CFL, dt, G, &
                         CS, visc_rem, j, ish, ieh, do_i, .true.)

  !   Determine the southerly- and northerly- fluxes.  Choose a sufficiently
  ! negative velocity correction for the northerly-flux, and a sufficiently
  ! positive correction for the southerly-flux.
  domore = .false.
  do i=ish,ieh ; if (do_i(i)) then
    domore = .true.
    dv_CFL(i) = (CFL_min * Idt) * G%dyCv(i,J)
    dvR(i) = min(0.0,dv0(i) - dv_CFL(i))
    dvL(i) = max(0.0,dv0(i) + dv_CFL(i))
    FAmt_L(i) = 0.0 ; FAmt_R(i) = 0.0 ; FAmt_0(i) = 0.0
    vhtot_L(i) = 0.0 ; vhtot_R(i) = 0.0
  endif ; enddo

  if (.not.domore) then
    do k=1,nz ; do i=ish,ieh
      BT_cont%FA_v_S0(i,J) = 0.0 ; BT_cont%FA_v_SS(i,J) = 0.0
      BT_cont%vBT_SS(i,J) = 0.0
      BT_cont%FA_v_N0(i,J) = 0.0 ; BT_cont%FA_v_NN(i,J) = 0.0
      BT_cont%vBT_NN(i,J) = 0.0
    enddo ; enddo
    return
  endif

  do k=1,nz ; do i=ish,ieh ; if (do_i(i)) then
    visc_rem_lim = max(visc_rem(i,k), min_visc_rem*visc_rem_max(i))
    if (v(i,J,k) + dvR(i)*visc_rem_lim > -dv_CFL(i)*visc_rem(i,k)) &
      dvR(i) = -(v(i,J,k) + dv_CFL(i)*visc_rem(i,k)) / visc_rem_lim
    if (v(i,J,k) + dvL(i)*visc_rem_lim < dv_CFL(i)*visc_rem(i,k)) &
      dvL(i) = -(v(i,J,k) - dv_CFL(i)*visc_rem(i,k)) / visc_rem_lim
  endif ; enddo ; enddo
  do k=1,nz
    do i=ish,ieh ; if (do_i(i)) then
      v_L(i) = v(I,j,k) + dvL(i) * visc_rem(i,k)
      v_R(i) = v(I,j,k) + dvR(i) * visc_rem(i,k)
      v_0(i) = v(I,j,k) + dv0(i) * visc_rem(i,k)
    endif ; enddo
    call merid_flux_layer(v_0, h_in(:,:,k), hL(:,:,k), hR(:,:,k), vh_0, &
                          FA_marg_0, visc_rem(:,k), dt, G, J, ish, ieh, do_i, &
                          CS%vol_CFL)
    call merid_flux_layer(v_L, h_in(:,:,k), hL(:,:,k), hR(:,:,k), vh_L, &
                          FA_marg_L, visc_rem(:,k), dt, G, J, ish, ieh, do_i, &
                          CS%vol_CFL)
    call merid_flux_layer(v_R, h_in(:,:,k), hL(:,:,k), hR(:,:,k), vh_R, &
                          FA_marg_R, visc_rem(:,k), dt, G, J, ish, ieh, do_i, &
                          CS%vol_CFL)
    do i=ish,ieh ; if (do_i(i)) then
      FAmt_0(i) = FAmt_0(i) + FA_marg_0(i)
      FAmt_L(i) = FAmt_L(i) + FA_marg_L(i)
      FAmt_R(i) = FAmt_R(i) + FA_marg_R(i)
      vhtot_L(i) = vhtot_L(i) + vh_L(i)
      vhtot_R(i) = vhtot_R(i) + vh_R(i)
    endif ; enddo
  enddo
  do i=ish,ieh ; if (do_i(i)) then
    FA_0 = FAmt_0(i) ; FA_avg = FAmt_0(i)
    if ((dvL(i) - dv0(i)) /= 0.0) &
      FA_avg = vhtot_L(i) / (dvL(i) - dv0(i))
    if (FA_avg > max(FA_0, FAmt_L(i))) then ; FA_avg = max(FA_0, FAmt_L(i))
    elseif (FA_avg < min(FA_0, FAmt_L(i))) then ; FA_0 = FA_avg ; endif
    BT_cont%FA_v_S0(i,J) = FA_0 ; BT_cont%FA_v_SS(i,J) = FAmt_L(i)
    if (abs(FA_0-FAmt_L(i)) <= 1e-12*FA_0) then ; BT_cont%vBT_SS(i,J) = 0.0 ; else
      BT_cont%vBT_SS(i,J) = (1.5 * (dvL(i) - dv0(i))) * &
                   ((FAmt_L(i) - FA_avg) / (FAmt_L(i) - FA_0))
    endif

    FA_0 = FAmt_0(i) ; FA_avg = FAmt_0(i)
    if ((dvR(i) - dv0(i)) /= 0.0) &
      FA_avg = vhtot_R(i) / (dvR(i) - dv0(i))
    if (FA_avg > max(FA_0, FAmt_R(i))) then ; FA_avg = max(FA_0, FAmt_R(i))
    elseif (FA_avg < min(FA_0, FAmt_R(i))) then ; FA_0 = FA_avg ; endif
    BT_cont%FA_v_N0(i,J) = FA_0 ; BT_cont%FA_v_NN(i,J) = FAmt_R(i)
    if (abs(FAmt_R(i) - FA_0) <= 1e-12*FA_0) then ; BT_cont%vBT_NN(i,J) = 0.0 ; else
      BT_cont%vBT_NN(i,J) = (1.5 * (dvR(i) - dv0(i))) * &
                   ((FAmt_R(i) - FA_avg) / (FAmt_R(i) - FA_0))
    endif
  else
    BT_cont%FA_v_S0(i,J) = 0.0 ; BT_cont%FA_v_SS(i,J) = 0.0
    BT_cont%vBT_SS(i,J) = 0.0
    BT_cont%FA_v_N0(i,J) = 0.0 ; BT_cont%FA_v_NN(i,J) = 0.0
    BT_cont%vBT_NN(i,J) = 0.0
  endif ; enddo

end subroutine set_merid_BT_cont

subroutine PPM_reconstruction_x(h_in, h_l, h_r, G, LB, h_min, monotonic, simple_2nd)
  real, dimension(NIMEM_,NJMEM_), intent(in)  :: h_in
  real, dimension(NIMEM_,NJMEM_), intent(out) :: h_l, h_r
  type(ocean_grid_type),          intent(in)  :: G
  type(loop_bounds_type),         intent(in)  :: LB
  real,                           intent(in)  :: h_min
  logical, optional,              intent(in)  :: monotonic
  logical, optional,              intent(in)  :: simple_2nd
! This subroutine calculates left/right edge valus for PPM reconstruction.
! Arguments: h_in    - thickness of layer (2D)
!  (out)     h_l,h_r - left/right edge value of reconstruction (2D)
!  (in)      G - The ocean's grid structure.
!  (in)      LB - A structure with the active loop bounds.
!  (in)      h_min   - The minimum thickness that can be obtained by a
!                      concave parabolic fit.
!  (in, opt) monotonic - If true, use the Colella & Woodward monotonic limiter.
!                        Otherwise use a simple positive-definite limiter.
!  (in, opt) simple_2nd - If true, use the arithmetic mean thicknesses as the
!                         default edge values for a simple 2nd order scheme.

! Local variables with useful mnemonic names.
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slopes.
  real, parameter :: oneSixth = 1./6.
  real :: h_ip1, h_im1
  real :: dMx, dMn
  logical :: use_CW84, use_2nd
  character(len=256) :: mesg
  integer :: i, j, isl, iel, jsl, jel, stensil

  use_CW84 = .false. ; if (present(monotonic)) use_CW84 = monotonic
  use_2nd = .false. ; if (present(simple_2nd)) use_2nd = simple_2nd
  isl = LB%ish-1 ; iel = LB%ieh+1 ; jsl = LB%jsh ; jel = LB%jeh

  ! This is the stensil of the reconstruction, not the scheme overall.
  stensil = 2 ; if (use_2nd) stensil = 1

  if ((isl-stensil < G%isd) .or. (iel+stensil > G%ied)) then
    write(mesg,'("In MOM_continuity_PPM, PPM_reconstruction_x called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               stensil + max(G%isd-isl,iel-G%ied)
    call MOM_error(FATAL,mesg)
  endif
  if ((jsl < G%jsd) .or. (jel > G%jed)) then
    write(mesg,'("In MOM_continuity_PPM, PPM_reconstruction_x called with a ", &
               & "y-halo that needs to be increased by ",i2,".")') &
               max(G%jsd-jsl,jel-G%jed)
    call MOM_error(FATAL,mesg)
  endif

  if (use_2nd) then
    do j=jsl,jel ; do i=isl,iel
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j)
      h_l(i,j) = 0.5*( h_im1 + h_in(i,j) )
      h_r(i,j) = 0.5*( h_ip1 + h_in(i,j) )
    enddo ; enddo
  else
    do j=jsl,jel ; do i=isl-1,iel+1
      if ((G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j)) == 0.0) then
        slp(i,j) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j) = 0.5 * (h_in(i+1,j) - h_in(i-1,j))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i+1,j), h_in(i-1,j), h_in(i,j)) - h_in(i,j)
        dMn = h_in(i,j) - min(h_in(i+1,j), h_in(i-1,j), h_in(i,j))
        slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dMx, dMn))
                ! * (G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j))
      endif
    enddo; enddo

    do j=jsl,jel ; do i=isl,iel
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
    ! h_im1 = h_in(i-1,j,k) ; if (G%mask2dT(i-1,j) < 0.5) h_im1 = h_in(i,j)
    ! h_ip1 = h_in(i+1,j,k) ; if (G%mask2dT(i+1,j) < 0.5) h_ip1 = h_in(i,j)
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_l(i,j) = 0.5*( h_im1 + h_in(i,j) ) + oneSixth*( slp(i-1,j) - slp(i,j) )
      h_r(i,j) = 0.5*( h_ip1 + h_in(i,j) ) + oneSixth*( slp(i,j) - slp(i+1,j) )
    enddo; enddo
  endif

  if (use_CW84) then
    call PPM_limit_CW84(h_in, h_l, h_r, isl, iel, jsl, jel)
  else
    call PPM_limit_pos(h_in, h_l, h_r, h_min, isl, iel, jsl, jel)
  endif

  return
end subroutine PPM_reconstruction_x

subroutine PPM_reconstruction_y(h_in, h_l, h_r, G, LB, h_min, monotonic, simple_2nd)
  real, dimension(NIMEM_,NJMEM_), intent(in)  :: h_in
  real, dimension(NIMEM_,NJMEM_), intent(out) :: h_l, h_r
  type(ocean_grid_type),          intent(in)  :: G
  type(loop_bounds_type),         intent(in)  :: LB
  real,                           intent(in)  :: h_min
  logical, optional,              intent(in)  :: monotonic
  logical, optional,              intent(in)  :: simple_2nd
! This subroutine calculates left/right edge valus for PPM reconstruction.
! Arguments: h_in    - thickness of layer (2D)
!  (out)     h_l,h_r - left/right edge value of reconstruction (2D)
!  (in)      G - The ocean's grid structure.
!  (in)      LB - A structure with the active loop bounds.
!  (in)      h_min   - The minimum thickness that can be obtained by a
!                      concave parabolic fit.
!  (in, opt) monotonic - If true, use the Colella & Woodward monotonic limiter.
!                        Otherwise use a simple positive-definite limiter.
!  (in, opt) simple_2nd - If true, use the arithmetic mean thicknesses as the
!                         default edge values for a simple 2nd order scheme.

! Local variables with useful mnemonic names.
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slopes.
  real, parameter :: oneSixth = 1./6.
  real :: h_jp1, h_jm1
  real :: dMx, dMn
  logical :: use_CW84, use_2nd
  character(len=256) :: mesg
  integer :: i, j, isl, iel, jsl, jel, stensil

  use_CW84 = .false. ; if (present(monotonic)) use_CW84 = monotonic
  use_2nd = .false. ; if (present(simple_2nd)) use_2nd = simple_2nd
  isl = LB%ish ; iel = LB%ieh ; jsl = LB%jsh-1 ; jel = LB%jeh+1

  ! This is the stensil of the reconstruction, not the scheme overall.
  stensil = 2 ; if (use_2nd) stensil = 1

  if ((isl < G%isd) .or. (iel > G%ied)) then
    write(mesg,'("In MOM_continuity_PPM, PPM_reconstruction_y called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               max(G%isd-isl,iel-G%ied)
    call MOM_error(FATAL,mesg)
  endif
  if ((jsl-stensil < G%jsd) .or. (jel+stensil > G%jed)) then
    write(mesg,'("In MOM_continuity_PPM, PPM_reconstruction_y called with a ", &
                 & "y-halo that needs to be increased by ",i2,".")') &
                 stensil + max(G%jsd-jsl,jel-G%jed)
    call MOM_error(FATAL,mesg)
  endif

  if (use_2nd) then
    do j=jsl,jel ; do i=isl,iel
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j)
      h_l(i,j) = 0.5*( h_jm1 + h_in(i,j) )
      h_r(i,j) = 0.5*( h_jp1 + h_in(i,j) )
    enddo ; enddo
  else
    do j=jsl-1,jel+1 ; do i=isl,iel
      if ((G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1)) == 0.0) then
        slp(i,j) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j) = 0.5 * (h_in(i,j+1) - h_in(i,j-1))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i,j+1), h_in(i,j-1), h_in(i,j)) - h_in(i,j)
        dMn = h_in(i,j) - min(h_in(i,j+1), h_in(i,j-1), h_in(i,j))
        slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dMx, dMn))
                ! * (G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1))
      endif
    enddo ; enddo

    do j=jsl,jel ; do i=isl,iel
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_l(i,j) = 0.5*( h_jm1 + h_in(i,j) ) + oneSixth*( slp(i,j-1) - slp(i,j) )
      h_r(i,j) = 0.5*( h_jp1 + h_in(i,j) ) + oneSixth*( slp(i,j) - slp(i,j+1) )
    enddo ; enddo
  endif

  if (use_CW84) then
    call PPM_limit_CW84(h_in, h_l, h_r, isl, iel, jsl, jel)
  else
    call PPM_limit_pos(h_in, h_l, h_r, h_min, isl, iel, jsl, jel)
  endif

  return
end subroutine PPM_reconstruction_y

subroutine PPM_limit_pos(h_in, h_L, h_R, h_min, iis, iie, jis, jie)
  real, dimension(NIMEM_,NJMEM_), intent(in)     :: h_in
  real, dimension(NIMEM_,NJMEM_), intent(inout)  :: h_L, h_R
  real,                           intent(in)     :: h_min
  integer,                        intent(in)     :: iis, iie, jis, jie
! This subroutine limits the left/right edge values of the PPM reconstruction
! to give a reconstruction that is positive-definite.  Here this is
! reinterpreted as giving a constant thickness if the mean thickness is less
! than h_min, with a minimum of h_min otherwise.
! Arguments: h_in    - thickness of layer (2D)
!  (inout)   h_L     - left edge value (2D)
!  (inout)   h_R     - right edge value (2D)
!  (in)      h_min   - The minimum thickness that can be obtained by a
!                      concave parabolic fit.
!  (in)      iis, iie, jis, jie - Index range for computation.

! Local variables
  real    :: curv, dh, scale
  character(len=256) :: mesg
  integer :: i,j

  do j=jis,jie ; do i=iis,iie
    ! This limiter prevents undershooting minima within the domain with
    ! values less than h_min.
    curv = 3.0*(h_L(i,j) + h_R(i,j) - 2.0*h_in(i,j))
    if (curv > 0.0) then ! Only minima are limited.
      dh = h_R(i,j) - h_L(i,j)
      if (abs(dh) < curv) then ! The parabola's minimum is within the cell.
        if (h_in(i,j) <= h_min) then
          h_L(i,j) = h_in(i,j) ; h_R(i,j) = h_in(i,j)
        elseif (12.0*curv*(h_in(i,j) - h_min) < (curv**2 + 3.0*dh**2)) then
          ! The minimum value is h_in - (curv^2 + 3*dh^2)/(12*curv), and must
          ! be limited in this case.  0 < scale < 1.
          scale = 12.0*curv*(h_in(i,j) - h_min) / (curv**2 + 3.0*dh**2)
          h_L(i,j) = h_in(i,j) + scale*(h_L(i,j) - h_in(i,j))
          h_R(i,j) = h_in(i,j) + scale*(h_R(i,j) - h_in(i,j))
        endif
      endif
    endif
  enddo ; enddo

end subroutine PPM_limit_pos

subroutine PPM_limit_CW84(h_in, h_l, h_r, iis, iie, jis, jie)
  real, dimension(NIMEM_,NJMEM_), intent(in)     :: h_in
  real, dimension(NIMEM_,NJMEM_), intent(inout)  :: h_l, h_r
  integer,                        intent(in)     :: iis, iie, jis, jie
! This subroutine limits the left/right edge values of the PPM reconstruction
! according to the monotonic prescription of Colella and Woodward, 1984.
! Arguments: h_in    - thickness of layer (2D)
!  (inout)   h_l     - left edge value (2D)
!  (inout)   h_r     - right edge value (2D)
!  (in)      iis, iie, jis, jie - Index range for computation.

! Local variables
  real    :: h_i, RLdiff, RLdiff2, RLmean, FunFac
  character(len=256) :: mesg
  integer :: i,j

  do j=jis,jie ; do i=iis,iie
    ! This limiter monotonizes the parabola following
    ! Colella and Woodward, 1984, Eq. 1.10
    h_i = h_in(i,j)
    if ( ( h_r(i,j) - h_i ) * ( h_i - h_l(i,j) ) <= 0. ) then
      h_l(i,j) = h_i ; h_r(i,j) = h_i
    else
      RLdiff = h_r(i,j) - h_l(i,j)            ! Difference of edge values
      RLmean = 0.5 * ( h_r(i,j) + h_l(i,j) )  ! Mean of edge values
      FunFac = 6. * RLdiff * ( h_i - RLmean ) ! Some funny factor
      RLdiff2 = RLdiff * RLdiff               ! Square of difference
      if ( FunFac >  RLdiff2 ) h_l(i,j) = 3. * h_i - 2. * h_r(i,j)
      if ( FunFac < -RLdiff2 ) h_r(i,j) = 3. * h_i - 2. * h_l(i,j)
    endif
  enddo ; enddo

  return
end subroutine PPM_limit_CW84

function ratio_max(a, b, maxrat) result(ratio)
  real, intent(in) :: a, b, maxrat
  real :: ratio

  if (abs(a) > abs(maxrat*b)) then
    ratio = maxrat
  else
    ratio = a / b
  endif
end function ratio_max

subroutine continuity_PPM_init(Time, G, GV, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(verticalGrid_type), intent(in)    :: GV
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ctrl), target, intent(inout) :: diag
  type(continuity_PPM_CS), pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_continuity_PPM" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "continuity_PPM_init called with associated control structure.")
    return
  endif
  allocate(CS)

! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "MONOTONIC_CONTINUITY", CS%monotonic, &
                 "If true, CONTINUITY_PPM uses the Colella and Woodward \n"//&
                 "monotonic limiter.  The default (false) is to use a \n"//&
                 "simple positive definite limiter.", default=.false.)
  call get_param(param_file, mod, "SIMPLE_2ND_PPM_CONTINUITY", CS%simple_2nd, &
                 "If true, CONTINUITY_PPM uses a simple 2nd order \n"//&
                 "(arithmetic mean) interpolation of the edge values. \n"//&
                 "This may give better PV conservation propterties. While \n"//&
                 "it formally reduces the accuracy of the continuity \n"//&
                 "solver itself in the strongly advective limit, it does \n"//&
                 "not reduce the overall order of accuracy of the dynamic \n"//&
                 "core.", default=.false.)
  call get_param(param_file, mod, "UPWIND_1ST_CONTINUITY", CS%upwind_1st, &
                 "If true, CONTINUITY_PPM becomes a 1st-order upwind \n"//&
                 "continuity solver.  This scheme is highly diffusive \n"//&
                 "but may be useful for debugging or in single-column \n"//&
                 "mode where its minimal stensil is useful.", default=.false.)
  call get_param(param_file, mod, "ETA_TOLERANCE", CS%tol_eta, &
                 "The tolerance for the differences between the \n"//&
                 "barotropic and baroclinic estimates of the sea surface \n"//&
                 "height due to the fluxes through each face.  The total \n"//&
                 "tolerance for SSH is 4 times this value.  The default \n"//&
                 "is 0.5*NK*ANGSTROM, and this should not be set less x\n"//&
                 "than about 10^-15*MAXIMUM_DEPTH.", units="m", &
                 default=0.5*G%ke*GV%Angstrom_z)

  call get_param(param_file, mod, "ETA_TOLERANCE_AUX", CS%tol_eta_aux, &
                 "The tolerance for free-surface height discrepancies \n"//&
                 "between the barotropic solution and the sum of the \n"//&
                 "layer thicknesses when calculating the auxiliary \n"//&
                 "corrected velocities. By default, this is the same as \n"//&
                 "ETA_TOLERANCE, but can be made larger for efficiency.", &
                 units="m", default=CS%tol_eta)
  call get_param(param_file, mod, "VELOCITY_TOLERANCE", CS%tol_vel, &
                 "The tolerance for barotropic velocity discrepancies \n"//&
                 "between the barotropic solution and  the sum of the \n"//&
                 "layer thicknesses.", units="m s-1", default=3.0e8) ! The speed of light is the default.

  call get_param(param_file, mod, "CONT_PPM_AGGRESS_ADJUST", CS%aggress_adjust,&
                 "If true, allow the adjusted velocities to have a \n"//&
                 "relative CFL change up to 0.5.", default=.false.)
  CS%vol_CFL = CS%aggress_adjust
  call get_param(param_file, mod, "CONT_PPM_VOLUME_BASED_CFL", CS%vol_CFL, &
                 "If true, use the ratio of the open face lengths to the \n"//&
                 "tracer cell areas when estimating CFL numbers.  The \n"//&
                 "default is set by CONT_PPM_AGGRESS_ADJUST.", &
                 default=CS%aggress_adjust, do_not_read=CS%aggress_adjust)
  call get_param(param_file, mod, "CONTINUITY_CFL_LIMIT", CS%CFL_limit_adjust, &
                 "The maximum CFL of the adjusted velocities.", units="nondim", &
                 default=0.5)
  call get_param(param_file, mod, "CONT_PPM_BETTER_ITER", CS%better_iter, &
                 "If true, stop corrective iterations using a velocity \n"//&
                 "based criterion and only stop if the iteration is \n"//&
                 "better than all predecessors.", default=.true.)
  call get_param(param_file, mod, "CONT_PPM_USE_VISC_REM_MAX", &
                                 CS%use_visc_rem_max, &
                 "If true, use more appropriate limiting bounds for \n"//&
                 "corrections in strongly viscous columns.", default=.true.)
  call get_param(param_file, mod, "CONT_PPM_MARGINAL_FACE_AREAS", CS%marginal_faces, &
                 "If true, use the marginal face areas from the continuity \n"//&
                 "solver for use as the weights in the barotropic solver. \n"//&
                 "Otherwise use the transport averaged areas.", default=.true.)

  CS%diag => diag

  id_clock_update = cpu_clock_id('(Ocean continuity update)', grain=CLOCK_ROUTINE)
  id_clock_correct = cpu_clock_id('(Ocean continuity correction)', grain=CLOCK_ROUTINE)

  CS%tol_eta = CS%tol_eta * GV%m_to_H
  CS%tol_eta_aux = CS%tol_eta_aux * GV%m_to_H

end subroutine continuity_PPM_init

subroutine continuity_PPM_end(CS)
  type(continuity_PPM_CS), pointer :: CS
  deallocate(CS)
end subroutine continuity_PPM_end

end module MOM_continuity_PPM
