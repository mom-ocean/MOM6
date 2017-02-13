module MOM_dynamics_unsplit

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
!*  By Robert Hallberg, 1993-2012                                      *
!*                                                                     *
!*    This file contains code that does the time-stepping of the       *
!*  adiabatic dynamic core, in this case with an unsplit third-order   *
!*  Runge-Kutta time stepping scheme for the momentum and a forward-   *
!*  backward coupling between the momentum and continuity equations.   *
!*  This was the orignal unsplit time stepping scheme used in early    *
!*  versions of HIM and its precuror.  While it is very simple and     *
!*  accurate, it is much less efficient that the split time stepping   *
!*  scheme for realistic oceanographic applications.  It has been      *
!*  retained for all of these years primarily to verify that the split *
!*  scheme is giving the right answers, and to debug the failings of   *
!*  the split scheme when it is not.  The split time stepping scheme   *
!*  is now sufficiently robust that it should be first choice for      *
!*  almost any conceivable application, except perhaps from cases      *
!*  with just a few layers for which the exact timing of the high-     *
!*  frequency barotropic gravity waves is of paramount importance.     *
!*  This scheme is slightly more efficient than the other unsplit      *
!*  scheme that can be found in MOM_dynamics_unsplit_RK2.F90.          *
!*                                                                     *
!*    The subroutine step_MOM_dyn_unsplit actually does the time       *
!*  stepping, while register_restarts_dyn_unsplit  sets the fields     *
!*  that are found in a full restart file with this scheme, and        *
!*  initialize_dyn_unsplit  initializes the cpu clocks that are        *                                      *
!*  used in this module.  For largely historical reasons, this module  *
!*  does not have its own control structure, but shares the same       *
!*  control structure with MOM.F90 and the other MOM_dynamics_...      *
!*  modules.                                                           *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, CoriolisBu                            *
!*    j+1  > o > o >   At ^:  v, PFv, CAv, vh, diffv, tauy, vbt, vhtr  *
!*    j    x ^ x ^ x   At >:  u, PFu, CAu, uh, diffu, taux, ubt, uhtr  *
!*    j    > o > o >   At o:  h, bathyT, eta, T, S, tr                 *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1                                                  *
!*           i  i+1                                                    *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**


use MOM_variables, only : vertvisc_type, thermo_var_ptrs
use MOM_variables, only : accel_diag_ptrs, ocean_internal_state, cont_diag_ptrs
use MOM_forcing_type, only : forcing
use MOM_checksum_packages, only : MOM_thermo_chksum, MOM_state_chksum, MOM_accel_chksum
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
use MOM_cpu_clock, only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : diag_mediator_init, enable_averaging
use MOM_diag_mediator, only : disable_averaging, post_data, safe_alloc_ptr
use MOM_diag_mediator, only : register_diag_field, register_static_field
use MOM_diag_mediator, only : set_diag_mediator_grid, diag_ctrl, diag_update_remap_grids
use MOM_domains, only : MOM_domains_init, pass_var, pass_vector
use MOM_domains, only : pass_var_start, pass_var_complete
use MOM_domains, only : pass_vector_start, pass_vector_complete
use MOM_domains, only : To_South, To_West, To_All, CGRID_NE, SCALAR_PAIR
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_io, only : MOM_io_init, vardesc
use MOM_restart, only : register_restart_field, query_initialized, save_restart
use MOM_restart, only : restart_init, MOM_restart_CS
use MOM_time_manager, only : time_type, set_time, time_type_to_real, operator(+)
use MOM_time_manager, only : operator(-), operator(>), operator(*), operator(/)

use MOM_ALE, only : ALE_CS
use MOM_continuity, only : continuity, continuity_init, continuity_CS
use MOM_CoriolisAdv, only : CorAdCalc, CoriolisAdv_init, CoriolisAdv_CS
use MOM_debugging, only : check_redundant
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_hor_visc, only : horizontal_viscosity, hor_visc_init, hor_visc_CS
use MOM_interface_heights, only : find_eta
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE_types, only : MEKE_type
use MOM_open_boundary, only : ocean_OBC_type
use MOM_open_boundary, only : radiation_open_bdry_conds
use MOM_boundary_update, only : update_OBC_data
use MOM_PressureForce, only : PressureForce, PressureForce_init, PressureForce_CS
use MOM_set_visc, only : set_viscous_BBL, set_viscous_ML, set_visc_CS
use MOM_tidal_forcing, only : tidal_forcing_init, tidal_forcing_CS
use MOM_vert_friction, only : vertvisc, vertvisc_coef
use MOM_vert_friction, only : vertvisc_limit_vel, vertvisc_init, vertvisc_CS
use MOM_verticalGrid, only : verticalGrid_type, get_thickness_units
use MOM_verticalGrid, only : get_flux_units, get_tr_flux_units

implicit none ; private

#include <MOM_memory.h>
type, public :: MOM_dyn_unsplit_CS ; private
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    CAu, &    ! CAu = f*v - u.grad(u) in m s-2.
    PFu, &    ! PFu = -dM/dx, in m s-2.
    diffu     ! Zonal acceleration due to convergence of the along-isopycnal
              ! stress tensor, in m s-2.

  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    CAv, &    ! CAv = -f*u - u.grad(v) in m s-2.
    PFv, &    ! PFv = -dM/dy, in m s-2.
    diffv     ! Meridional acceleration due to convergence of the
              ! along-isopycnal stress tensor, in m s-2.

  real, pointer, dimension(:,:) :: taux_bot => NULL(), tauy_bot => NULL()
    ! The frictional bottom stresses from the ocean to the seafloor, in Pa.

  logical :: debug           ! If true, write verbose checksums for debugging purposes.

  logical :: module_is_initialized = .false.

  integer :: id_uh = -1, id_vh = -1
  integer :: id_PFu = -1, id_PFv = -1, id_CAu = -1, id_CAv = -1

  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                                   ! timing of diagnostic output.
  type(accel_diag_ptrs), pointer :: ADp ! A structure pointing to the various
                                   ! accelerations in the momentum equations,
                                   ! which can later be used to calculate
                                   ! derived diagnostics like energy budgets.
  type(cont_diag_ptrs), pointer :: CDp ! A structure with pointers to various
                                   ! terms in the continuity equations,
                                   ! which can later be used to calculate
                                   ! derived diagnostics like energy budgets.
! The remainder of the structure is pointers to child subroutines' control strings.
  type(hor_visc_CS), pointer :: hor_visc_CSp => NULL()
  type(continuity_CS), pointer :: continuity_CSp => NULL()
  type(CoriolisAdv_CS), pointer :: CoriolisAdv_CSp => NULL()
  type(PressureForce_CS), pointer :: PressureForce_CSp => NULL()
  type(vertvisc_CS), pointer :: vertvisc_CSp => NULL()
  type(set_visc_CS), pointer :: set_visc_CSp => NULL()
  type(ocean_OBC_type), pointer :: OBC => NULL() ! A pointer to an open boundary
     ! condition type that specifies whether, where, and  what open boundary
     ! conditions are used.  If no open BCs are used, this pointer stays
     ! nullified.  Flather OBCs use open boundary_CS as well.
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL()

! This is a copy of the pointer in the top-level control structure.
  type(ALE_CS), pointer :: ALE_CSp => NULL()

end type MOM_dyn_unsplit_CS

public step_MOM_dyn_unsplit, register_restarts_dyn_unsplit
public initialize_dyn_unsplit, end_dyn_unsplit

integer :: id_clock_Cor, id_clock_pres, id_clock_vertvisc
integer :: id_clock_continuity, id_clock_horvisc, id_clock_mom_update
integer :: id_clock_pass, id_clock_pass_init

contains

! =============================================================================

subroutine step_MOM_dyn_unsplit(u, v, h, tv, visc, Time_local, dt, fluxes, &
                  p_surf_begin, p_surf_end, uh, vh, uhtr, vhtr, eta_av, G, GV, CS, &
                  VarMix, MEKE)
  type(ocean_grid_type),                     intent(inout) :: G
  type(verticalGrid_type),                   intent(in)    :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h
  type(thermo_var_ptrs),                     intent(in)    :: tv
  type(vertvisc_type),                       intent(inout) :: visc
  type(time_type),                           intent(in)    :: Time_local
  real,                                      intent(in)    :: dt
  type(forcing),                             intent(in)    :: fluxes
  real, dimension(:,:),                      pointer       :: p_surf_begin, p_surf_end
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uh
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vh
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uhtr
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vhtr
  real, dimension(SZI_(G),SZJ_(G)),          intent(out)   :: eta_av
  type(MOM_dyn_unsplit_CS),                  pointer       :: CS
  type(VarMix_CS),                           pointer       :: VarMix
  type(MEKE_type),                           pointer       :: MEKE
! Arguments: u - The input and output zonal velocity, in m s-1.
!  (inout)   v - The input and output meridional velocity, in m s-1.
!  (inout)   h - The input and output layer thicknesses, in m or kg m-2,
!                depending on whether the Boussinesq approximation is made.
!  (in)      tv - a structure pointing to various thermodynamic variables.
!  (inout)   visc - A structure containing vertical viscosities, bottom drag
!                   viscosities, and related fields.
!  (in)      Time_local - The model time at the end of the time step.
!  (in)      dt - The time step in s.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      p_surf_begin - A pointer (perhaps NULL) to the surface pressure
!                     at the beginning of this dynamic step, in Pa.
!  (in)      p_surf_end - A pointer (perhaps NULL) to the surface pressure
!                     at the end of this dynamic step, in Pa.
!  (inout)   uh - The zonal volume or mass transport, in m3 s-1 or kg s-1.
!  (inout)   vh - The meridional volume or mass transport, in m3 s-1 or kg s-1.
!  (inout)   uhtr - The accumulated zonal volume or mass transport since the last
!                   tracer advection, in m3 or kg.
!  (inout)   vhtr - The accumulated meridional volume or mass transport since the last
!                   tracer advection, in m3 or kg.
!  (out)     eta_av - The time-mean free surface height or column mass, in m or
!                     kg m-2.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure set up by initialize_dyn_unsplit.
!  (in)      VarMix - A pointer to a structure with fields that specify the
!                     spatially variable viscosities.
!  (inout)   MEKE - A pointer to a structure containing fields related to
!                   the Mesoscale Eddy Kinetic Energy.

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_av, hp
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: up, upp
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: vp, vpp
  real, dimension(:,:), pointer :: p_surf
  real :: dt_pred   ! The time step for the predictor part of the baroclinic
                    ! time stepping.
  logical :: dyn_p_surf
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  dt_pred = dt / 3.0

  h_av(:,:,:) = 0; hp(:,:,:) = 0
  up(:,:,:) = 0; upp(:,:,:) = 0
  vp(:,:,:) = 0; vpp(:,:,:) = 0

  dyn_p_surf = associated(p_surf_begin) .and. associated(p_surf_end)
  if (dyn_p_surf) then
    call safe_alloc_ptr(p_surf,G%isd,G%ied,G%jsd,G%jed) ; p_surf(:,:) = 0.0
  else
    p_surf => fluxes%p_surf
  endif

! Matsuno's third order accurate three step scheme is used to step
! all of the fields except h.  h is stepped separately.

  if (CS%debug) then
    call MOM_state_chksum("Start First Predictor ", u, v, h, uh, vh, G, GV)
  endif

! diffu = horizontal viscosity terms (u,h)
  call enable_averaging(dt,Time_local, CS%diag)
  call cpu_clock_begin(id_clock_horvisc)
  call horizontal_viscosity(u, v, h, CS%diffu, CS%diffv, MEKE, Varmix, &
                            G, GV, CS%hor_visc_CSp)
  call cpu_clock_end(id_clock_horvisc)
  call disable_averaging(CS%diag)

! uh = u*h
! hp = h + dt/2 div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(u, v, h, hp, uh, vh, dt*0.5, G, GV, CS%continuity_CSp, &
                  OBC=CS%OBC)
  call cpu_clock_end(id_clock_continuity)
  call cpu_clock_begin(id_clock_pass)
  call pass_var(hp, G%Domain)
  call pass_vector(uh, vh, G%Domain)
  call cpu_clock_end(id_clock_pass)

  call enable_averaging(0.5*dt,Time_local-set_time(int(0.5*dt)), CS%diag)
!   Here the first half of the thickness fluxes are offered for averaging.
  if (CS%id_uh > 0) call post_data(CS%id_uh, uh, CS%diag)
  if (CS%id_vh > 0) call post_data(CS%id_vh, vh, CS%diag)
  call disable_averaging(CS%diag)

! h_av = (h + hp)/2
! u = u + dt diffu
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      h_av(i,j,k) = (h(i,j,k) + hp(i,j,k)) * 0.5
    enddo ; enddo
    do j=js,je ; do I=Isq,Ieq
      u(I,j,k) = u(I,j,k) + dt * CS%diffu(I,j,k) * G%mask2dCu(I,j)
    enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie
      v(i,J,k) = v(i,J,k) + dt * CS%diffv(i,J,k) * G%mask2dCv(i,J)
    enddo ; enddo
    do j=js-2,je+2 ; do I=Isq-2,Ieq+2
      uhtr(i,j,k) = uhtr(i,j,k) + 0.5*dt*uh(i,j,k)
    enddo ; enddo
    do J=Jsq-2,Jeq+2 ; do i=is-2,ie+2
      vhtr(i,j,k) = vhtr(i,j,k) + 0.5*dt*vh(i,j,k)
    enddo ; enddo
  enddo
  call cpu_clock_end(id_clock_mom_update)
  call cpu_clock_begin(id_clock_pass)
  call pass_vector(u, v, G%Domain)
  call cpu_clock_end(id_clock_pass)

! CAu = -(f+zeta)/h_av vh + d/dx KE
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(u, v, h_av, uh, vh, CS%CAu, CS%CAv, CS%OBC, CS%ADp, &
                 G, GV, CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! PFu = d/dx M(h_av,T,S)
  call cpu_clock_begin(id_clock_pres)
  if (dyn_p_surf) then ; do j=js-2,je+2 ; do i=is-2,ie+2
    p_surf(i,j) = 0.75*p_surf_begin(i,j) + 0.25*p_surf_end(i,j)
  enddo ; enddo ; endif
  call PressureForce(h_av, tv, CS%PFu, CS%PFv, G, GV, &
                     CS%PressureForce_CSp, CS%ALE_CSp, p_surf)
  call cpu_clock_end(id_clock_pres)

  if (associated(CS%OBC)) then; if (CS%OBC%update_OBC) then
    call update_OBC_data(CS%OBC, G, GV, tv, h, Time_local)
  endif; endif

! up = u + dt_pred * (PFu + CAu)
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    up(I,j,k) = G%mask2dCu(I,j) * (u(I,j,k) + dt_pred * &
                               (CS%PFu(I,j,k) + CS%CAu(I,j,k)))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    vp(i,J,k) = G%mask2dCv(i,J) * (v(i,J,k) + dt_pred * &
                               (CS%PFv(i,J,k) + CS%CAv(i,J,k)))
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) then
    call MOM_state_chksum("Predictor 1", up, vp, h_av, uh, vh, G, GV)
    call MOM_accel_chksum("Predictor 1 accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv,&
                          CS%diffu, CS%diffv, G, GV)
  endif

! visc contains viscosity and BBL thickness (u_in,h_in)
  if (visc%calc_bbl) then
    call enable_averaging(visc%bbl_calc_time_interval, &
              Time_local+set_time(int(visc%bbl_calc_time_interval-dt)), CS%diag)
    call set_viscous_BBL(u, v, h_av, tv, visc, G, GV, CS%set_visc_CSp)
    call disable_averaging(CS%diag)
    call cpu_clock_begin(id_clock_pass)
    if (associated(visc%Ray_u) .and. associated(visc%Ray_v)) &
      call pass_vector(visc%Ray_u, visc%Ray_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE)
    if (associated(visc%kv_bbl_u) .and. associated(visc%kv_bbl_v)) then
      call pass_vector(visc%bbl_thick_u, visc%bbl_thick_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE, complete=.false.)
      call pass_vector(visc%kv_bbl_u, visc%kv_bbl_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE)
    endif
    call cpu_clock_end(id_clock_pass)
    visc%calc_bbl = .false.
  endif

 ! up <- up + dt/2 d/dz visc d/dz up
  call cpu_clock_begin(id_clock_vertvisc)
  call enable_averaging(dt, Time_local, CS%diag)
  call set_viscous_ML(u, v, h_av, tv, fluxes, visc, dt*0.5, G, GV, &
                      CS%set_visc_CSp)
  call disable_averaging(CS%diag)
  call vertvisc_coef(up, vp, h_av, fluxes, visc, dt*0.5, G, GV, CS%vertvisc_CSp)
  call vertvisc(up, vp, h_av, fluxes, visc, dt*0.5, CS%OBC, CS%ADp, CS%CDp, &
                G, GV, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)
  call cpu_clock_begin(id_clock_pass)
  call pass_vector(up, vp, G%Domain)
  call cpu_clock_end(id_clock_pass)

! uh = up * hp
! h_av = hp + dt/2 div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(up, vp, hp, h_av, uh, vh, &
                  (0.5*dt), G, GV, CS%continuity_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_continuity)
  call cpu_clock_begin(id_clock_pass)
  call pass_var(h_av, G%Domain)
  call pass_vector(uh, vh, G%Domain)
  call cpu_clock_end(id_clock_pass)

! h_av <- (hp + h_av)/2
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    h_av(i,j,k) = (hp(i,j,k) + h_av(i,j,k)) * 0.5
  enddo ; enddo ; enddo

! CAu = -(f+zeta(up))/h_av vh + d/dx KE(up)
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(up, vp, h_av, uh, vh, CS%CAu, CS%CAv, CS%OBC, CS%ADp, &
                 G, GV, CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! PFu = d/dx M(h_av,T,S)
  call cpu_clock_begin(id_clock_pres)
  if (dyn_p_surf) then ; do j=js-2,je+2 ; do i=is-2,ie+2
    p_surf(i,j) = 0.25*p_surf_begin(i,j) + 0.75*p_surf_end(i,j)
  enddo ; enddo ; endif
  call PressureForce(h_av, tv, CS%PFu, CS%PFv, G, GV, &
                     CS%PressureForce_CSp, CS%ALE_CSp, p_surf)
  call cpu_clock_end(id_clock_pres)

  if (associated(CS%OBC)) then; if (CS%OBC%update_OBC) then
    call update_OBC_data(CS%OBC, G, GV, tv, h, Time_local)
  endif; endif

! upp = u + dt/2 * ( PFu + CAu )
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    upp(I,j,k) = G%mask2dCu(I,j) * (u(I,j,k) + dt * 0.5 * &
            (CS%PFu(I,j,k) + CS%CAu(I,j,k)))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    vpp(i,J,k) = G%mask2dCv(i,J) * (v(i,J,k) + dt * 0.5 * &
            (CS%PFv(i,J,k) + CS%CAv(i,J,k)))
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) then
    call MOM_state_chksum("Predictor 2", upp, vpp, h_av, uh, vh, G, GV)
    call MOM_accel_chksum("Predictor 2 accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv,&
                          CS%diffu, CS%diffv, G, GV)
  endif

! upp <- upp + dt/2 d/dz visc d/dz upp
  call cpu_clock_begin(id_clock_vertvisc)
  call vertvisc_coef(upp, vpp, hp, fluxes, visc, dt*0.5, G, GV, CS%vertvisc_CSp)
  call vertvisc(upp, vpp, hp, fluxes, visc, dt*0.5, CS%OBC, CS%ADp, CS%CDp, &
                G, GV, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)
  call cpu_clock_begin(id_clock_pass)
  call pass_vector(upp, vpp, G%Domain)
  call cpu_clock_end(id_clock_pass)

! uh = upp * hp
! h = hp + dt/2 div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(upp, vpp, hp, h, uh, vh, &
                  (dt*0.5), G, GV, CS%continuity_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_continuity)
  call cpu_clock_begin(id_clock_pass)
  call pass_var(h, G%Domain)
  call pass_vector(uh, vh, G%Domain)
  call cpu_clock_end(id_clock_pass)
  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  call diag_update_remap_grids(CS%diag)

  call enable_averaging(0.5*dt, Time_local, CS%diag)
!   Here the second half of the thickness fluxes are offered for averaging.
  if (CS%id_uh > 0) call post_data(CS%id_uh, uh, CS%diag)
  if (CS%id_vh > 0) call post_data(CS%id_vh, vh, CS%diag)
  call disable_averaging(CS%diag)
  call enable_averaging(dt, Time_local, CS%diag)

! h_av = (h + hp)/2
  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      h_av(i,j,k) = 0.5*(h(i,j,k) + hp(i,j,k))
    enddo ; enddo
    do j=js-2,je+2 ; do I=Isq-2,Ieq+2
      uhtr(i,j,k) = uhtr(i,j,k) + 0.5*dt*uh(i,j,k)
    enddo ; enddo
    do J=Jsq-2,Jeq+2 ; do i=is-2,ie+2
      vhtr(i,j,k) = vhtr(i,j,k) + 0.5*dt*vh(i,j,k)
    enddo ; enddo
  enddo

! CAu = -(f+zeta(upp))/h_av vh + d/dx KE(upp)
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(upp, vpp, h_av, uh, vh, CS%CAu, CS%CAv, CS%OBC, CS%ADp, &
                 G, GV, CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! PFu = d/dx M(h_av,T,S)
  call cpu_clock_begin(id_clock_pres)
  call PressureForce(h_av, tv, CS%PFu, CS%PFv, G, GV, &
                     CS%PressureForce_CSp, CS%ALE_CSp, p_surf)
  call cpu_clock_end(id_clock_pres)

  if (associated(CS%OBC)) then; if (CS%OBC%update_OBC) then
    call update_OBC_data(CS%OBC, G, GV, tv, h, Time_local)
  endif; endif

! u = u + dt * ( PFu + CAu )
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u(I,j,k) = G%mask2dCu(I,j) * (u(I,j,k) + dt * &
            (CS%PFu(I,j,k) + CS%CAu(I,j,k)))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v(i,J,k) = G%mask2dCv(i,J) * (v(i,J,k) + dt * &
            (CS%PFv(i,J,k) + CS%CAv(i,J,k)))
  enddo ; enddo ; enddo

! u <- u + dt d/dz visc d/dz u
  call cpu_clock_begin(id_clock_vertvisc)
  call vertvisc_coef(u, v, h_av, fluxes, visc, dt, G, GV, CS%vertvisc_CSp)
  call vertvisc(u, v, h_av, fluxes, visc, dt, CS%OBC, CS%ADp, CS%CDp, &
                G, GV, CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot)
  call cpu_clock_end(id_clock_vertvisc)
  call cpu_clock_begin(id_clock_pass)
  call pass_vector(u, v, G%Domain)
  call cpu_clock_end(id_clock_pass)

  if (CS%debug) then
    call MOM_state_chksum("Corrector", u, v, h, uh, vh, G, GV)
    call MOM_accel_chksum("Corrector accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv, &
                          CS%diffu, CS%diffv, G, GV)
  endif

  if (GV%Boussinesq) then
    do j=js,je ; do i=is,ie ; eta_av(i,j) = -G%bathyT(i,j) ; enddo ; enddo
  else
    do j=js,je ; do i=is,ie ; eta_av(i,j) = 0.0 ; enddo ; enddo
  endif
  do k=1,nz ; do j=js,je ; do i=is,ie
    eta_av(i,j) = eta_av(i,j) + h_av(i,j,k)
  enddo ; enddo ; enddo

  if (dyn_p_surf) deallocate(p_surf)

!   Here various terms used in to update the momentum equations are
! offered for averaging.
  if (CS%id_PFu > 0) call post_data(CS%id_PFu, CS%PFu, CS%diag)
  if (CS%id_PFv > 0) call post_data(CS%id_PFv, CS%PFv, CS%diag)
  if (CS%id_CAu > 0) call post_data(CS%id_CAu, CS%CAu, CS%diag)
  if (CS%id_CAv > 0) call post_data(CS%id_CAv, CS%CAv, CS%diag)

end subroutine step_MOM_dyn_unsplit

! =============================================================================

subroutine register_restarts_dyn_unsplit(HI, GV, param_file, CS, restart_CS)
  type(hor_index_type),         intent(in)    :: HI
  type(verticalGrid_type),      intent(in)    :: GV
  type(param_file_type),        intent(in)    :: param_file
  type(MOM_dyn_unsplit_CS),     pointer       :: CS
  type(MOM_restart_CS),         pointer       :: restart_CS
!   This subroutine sets up any auxiliary restart variables that are specific
! to the unsplit time stepping scheme.  All variables registered here should
! have the ability to be recreated if they are not present in a restart file.

! Arguments: HI - A horizontal index type structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (inout)   CS - The control structure set up by initialize_dyn_unsplit.
!  (inout)   restart_CS - A pointer to the restart control structure.

  type(vardesc) :: vd
  character(len=40)  :: mod = "MOM_dynamics_unsplit" ! This module's name.
  character(len=48) :: thickness_units, flux_units
  integer :: isd, ied, jsd, jed, nz, IsdB, IedB, JsdB, JedB
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke
  IsdB = HI%IsdB ; IedB = HI%IedB ; JsdB = HI%JsdB ; JedB = HI%JedB

! This is where a control structure that is specific to this module would be allocated.
  if (associated(CS)) then
    call MOM_error(WARNING, "register_restarts_dyn_unsplit called with an associated "// &
                             "control structure.")
    return
  endif
  allocate(CS)

  ALLOC_(CS%diffu(IsdB:IedB,jsd:jed,nz)) ; CS%diffu(:,:,:) = 0.0
  ALLOC_(CS%diffv(isd:ied,JsdB:JedB,nz)) ; CS%diffv(:,:,:) = 0.0
  ALLOC_(CS%CAu(IsdB:IedB,jsd:jed,nz)) ; CS%CAu(:,:,:) = 0.0
  ALLOC_(CS%CAv(isd:ied,JsdB:JedB,nz)) ; CS%CAv(:,:,:) = 0.0
  ALLOC_(CS%PFu(IsdB:IedB,jsd:jed,nz)) ; CS%PFu(:,:,:) = 0.0
  ALLOC_(CS%PFv(isd:ied,JsdB:JedB,nz)) ; CS%PFv(:,:,:) = 0.0

  thickness_units = get_thickness_units(GV)
  flux_units = get_flux_units(GV)

!  No extra restart fields are needed with this time stepping scheme.

end subroutine register_restarts_dyn_unsplit

subroutine initialize_dyn_unsplit(u, v, h, Time, G, GV, param_file, diag, CS, &
                                  restart_CS, Accel_diag, Cont_diag, MIS, &
                                  OBC, ALE_CSp, setVisc_CSp, visc, dirs, ntrunc)
  type(ocean_grid_type),                     intent(inout) :: G
  type(verticalGrid_type),                   intent(in)    :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout) :: h
  type(time_type),                   target, intent(in)    :: Time
  type(param_file_type),                     intent(in)    :: param_file
  type(diag_ctrl),                   target, intent(inout) :: diag
  type(MOM_dyn_unsplit_CS),                  pointer       :: CS
  type(MOM_restart_CS),                      pointer       :: restart_CS
  type(accel_diag_ptrs),             target, intent(inout) :: Accel_diag
  type(cont_diag_ptrs),              target, intent(inout) :: Cont_diag
  type(ocean_internal_state),                intent(inout) :: MIS
  type(ocean_OBC_type),                      pointer       :: OBC
  type(ALE_CS),                              pointer       :: ALE_CSp
  type(set_visc_CS),                         pointer       :: setVisc_CSp
  type(vertvisc_type),                       intent(inout) :: visc
  type(directories),                         intent(in)    :: dirs
  integer, target,                           intent(inout) :: ntrunc

! Arguments: u - The zonal velocity, in m s-1.
!  (inout)   v - The meridional velocity, in m s-1.
!  (inout)   h - The layer thicknesses, in m or kg m-2, depending on whether
!                the Boussinesq approximation is made.
!  (in)      Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (inout)   CS - The control structure set up by initialize_dyn_unsplit.
!  (in)      restart_CS - A pointer to the restart control structure.
!  (inout)   Accel_diag - A set of pointers to the various accelerations in
!                  the momentum equations, which can be used for later derived
!                  diagnostics, like energy budgets.
!  (inout)   Cont_diag - A structure with pointers to various terms in the
!                   continuity equations.
!  (inout)   MIS - The "MOM6 Internal State" structure, used to pass around
!                  pointers to various arrays for diagnostic purposes.
!  (in)      OBC - If open boundary conditions are used, this points to the
!                  ocean_OBC_type that was set up in MOM_initialization.
!  (in)      ALE_CS - This points to the ALE control structure.
!  (in)      setVisc_CSp - This points to the set_visc control structure.
!  (inout)   visc - A structure containing vertical viscosities, bottom drag
!                   viscosities, and related fields.
!  (in)      dirs - A structure containing several relevant directory paths.
!  (in)      ntrunc - A target for the variable that records the number of times
!                     the velocity is truncated (this should be 0).

  !   This subroutine initializes all of the variables that are used by this
  ! dynamic core, including diagnostics and the cpu clocks.
  character(len=40) :: mod = "MOM_dynamics_unsplit" ! This module's name.
  character(len=48) :: thickness_units, flux_units
  logical :: use_tides
  integer :: isd, ied, jsd, jed, nz, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(CS)) call MOM_error(FATAL, &
      "initialize_dyn_unsplit called with an unassociated control structure.")
  if (CS%module_is_initialized) then
    call MOM_error(WARNING, "initialize_dyn_unsplit called with a control "// &
                            "structure that has already been initialized.")
    return
  endif
  CS%module_is_initialized = .true.

  CS%diag => diag

  call get_param(param_file, mod, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", default=.false.)
  call get_param(param_file, mod, "TIDES", use_tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)

  allocate(CS%taux_bot(IsdB:IedB,jsd:jed)) ; CS%taux_bot(:,:) = 0.0
  allocate(CS%tauy_bot(isd:ied,JsdB:JedB)) ; CS%tauy_bot(:,:) = 0.0

  MIS%diffu => CS%diffu ; MIS%diffv => CS%diffv
  MIS%PFu => CS%PFu ; MIS%PFv => CS%PFv
  MIS%CAu => CS%CAu ; MIS%CAv => CS%CAv

  CS%ADp => Accel_diag ; CS%CDp => Cont_diag
  Accel_diag%diffu => CS%diffu ; Accel_diag%diffv => CS%diffv
  Accel_diag%PFu => CS%PFu ; Accel_diag%PFv => CS%PFv
  Accel_diag%CAu => CS%CAu ; Accel_diag%CAv => CS%CAv

  call continuity_init(Time, G, GV, param_file, diag, CS%continuity_CSp)
  call CoriolisAdv_init(Time, G, param_file, diag, CS%ADp, CS%CoriolisAdv_CSp)
  if (use_tides) call tidal_forcing_init(Time, G, param_file, CS%tides_CSp)
  call PressureForce_init(Time, G, GV, param_file, diag, CS%PressureForce_CSp, &
                          CS%tides_CSp)
  call hor_visc_init(Time, G, param_file, diag, CS%hor_visc_CSp)
  call vertvisc_init(MIS, Time, G, GV, param_file, diag, CS%ADp, dirs, &
                     ntrunc, CS%vertvisc_CSp)
  if (.not.associated(setVisc_CSp)) call MOM_error(FATAL, &
    "initialize_dyn_unsplit called with setVisc_CSp unassociated.")
  CS%set_visc_CSp => setVisc_CSp

  if (associated(ALE_CSp)) CS%ALE_CSp => ALE_CSp
  if (associated(OBC)) CS%OBC => OBC

  flux_units = get_flux_units(GV)
  CS%id_uh = register_diag_field('ocean_model', 'uh', diag%axesCuL, Time, &
      'Zonal Thickness Flux', flux_units, y_cell_method='sum', v_extensive=.true.)
  CS%id_vh = register_diag_field('ocean_model', 'vh', diag%axesCvL, Time, &
      'Meridional Thickness Flux', flux_units, x_cell_method='sum', v_extensive=.true.)
  CS%id_CAu = register_diag_field('ocean_model', 'CAu', diag%axesCuL, Time, &
      'Zonal Coriolis and Advective Acceleration', 'meter second-2')
  CS%id_CAv = register_diag_field('ocean_model', 'CAv', diag%axesCvL, Time, &
      'Meridional Coriolis and Advective Acceleration', 'meter second-2')
  CS%id_PFu = register_diag_field('ocean_model', 'PFu', diag%axesCuL, Time, &
      'Zonal Pressure Force Acceleration', 'meter second-2')
  CS%id_PFv = register_diag_field('ocean_model', 'PFv', diag%axesCvL, Time, &
      'Meridional Pressure Force Acceleration', 'meter second-2')

  id_clock_Cor = cpu_clock_id('(Ocean Coriolis & mom advection)', grain=CLOCK_MODULE)
  id_clock_continuity = cpu_clock_id('(Ocean continuity equation)', grain=CLOCK_MODULE)
  id_clock_pres = cpu_clock_id('(Ocean pressure force)', grain=CLOCK_MODULE)
  id_clock_vertvisc = cpu_clock_id('(Ocean vertical viscosity)', grain=CLOCK_MODULE)
  id_clock_horvisc = cpu_clock_id('(Ocean horizontal viscosity)', grain=CLOCK_MODULE)
  id_clock_mom_update = cpu_clock_id('(Ocean momentum increments)', grain=CLOCK_MODULE)
  id_clock_pass = cpu_clock_id('(Ocean message passing)', grain=CLOCK_MODULE)
  id_clock_pass_init = cpu_clock_id('(Ocean init message passing)', grain=CLOCK_ROUTINE)

end subroutine initialize_dyn_unsplit

subroutine end_dyn_unsplit(CS)
  type(MOM_dyn_unsplit_CS), pointer :: CS
!  (inout)    CS - The control structure set up by initialize_dyn_unsplit.

  DEALLOC_(CS%diffu) ; DEALLOC_(CS%diffv)
  DEALLOC_(CS%CAu)   ; DEALLOC_(CS%CAv)
  DEALLOC_(CS%PFu)   ; DEALLOC_(CS%PFv)

  deallocate(CS)
end subroutine end_dyn_unsplit

end module MOM_dynamics_unsplit
