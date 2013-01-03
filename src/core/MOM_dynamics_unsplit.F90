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
!*    j+1  x ^ x ^ x   At x:  q, f                                     *
!*    j+1  > o > o >   At ^:  v, PFv, CAv, vh, diffv, tauy, vbt, vhtr  *
!*    j    x ^ x ^ x   At >:  u, PFu, CAu, uh, diffu, taux, ubt, uhtr  *
!*    j    > o > o >   At o:  h, D, eta, T, S, tr                      *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1                                                  *
!*           i  i+1                                                    *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**


use MOM_variables, only : directories, vertvisc_type, ocean_OBC_type
use MOM_variables, only : BT_cont_type, alloc_bt_cont_type, dealloc_bt_cont_type
use MOM_variables, only : &
  forcing, &      ! A structure containing pointers to the forcing fields
                  ! which may be used to drive MOM.  All fluxes are
                  ! positive downward.
  thermo_var_ptrs ! A structure containing pointers to an assortment of
                  ! thermodynamic fields that may be available, including
                  ! potential temperature, salinity and mixed layer density.

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
use MOM_cpu_clock, only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : diag_mediator_init, enable_averaging
use MOM_diag_mediator, only : disable_averaging, post_data, safe_alloc_ptr
use MOM_diag_mediator, only : register_diag_field, register_static_field
use MOM_diag_mediator, only : set_diag_mediator_grid, diag_ptrs
use MOM_domains, only : MOM_domains_init, pass_var, pass_vector
use MOM_domains, only : pass_var_start, pass_var_complete
use MOM_domains, only : pass_vector_start, pass_vector_complete
use MOM_domains, only : To_South, To_West, To_All, CGRID_NE, SCALAR_PAIR
use MOM_checksums, only : MOM_checksums_init, hchksum, uchksum, vchksum
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : MOM_set_verbosity
use MOM_file_parser, only : read_param, log_param, log_version, param_file_type
use MOM_io, only : MOM_io_init, vardesc
use MOM_obsolete_params, only : find_obsolete_params
use MOM_restart, only : register_restart_field, query_initialized, save_restart
use MOM_restart, only : restart_init, MOM_restart_CS
use MOM_time_manager, only : time_type, set_time, time_type_to_real, operator(+)
use MOM_time_manager, only : operator(-), operator(>), operator(*), operator(/)
use MOM_initialization, only : MOM_initialize, Get_MOM_Input
use MOM_initialization, only : MOM_initialization_struct

use MOM_continuity, only : continuity, continuity_init, continuity_CS
use MOM_CoriolisAdv, only : CorAdCalc, CoriolisAdv_init, CoriolisAdv_CS
use MOM_diabatic_driver, only : diabatic, diabatic_driver_init, diabatic_CS
use MOM_diagnostics, only : calculate_diagnostic_fields, MOM_diagnostics_init
use MOM_diagnostics, only : diagnostics_CS
use MOM_diag_to_Z, only : calculate_Z_diag_fields, calculate_Z_transport
use MOM_diag_to_Z, only : MOM_diag_to_Z_init, register_Z_tracer, diag_to_Z_CS
use MOM_diag_to_Z, only : MOM_diag_to_Z_end
use MOM_EOS, only : select_eqn_of_state
use MOM_error_checking, only : check_redundant
use MOM_grid, only : MOM_grid_init, ocean_grid_type, get_thickness_units
use MOM_grid, only : get_flux_units, get_tr_flux_units
use MOM_hor_visc, only : horizontal_viscosity, hor_visc_init, hor_visc_CS
use MOM_lateral_mixing_coeffs, only : calc_slope_function, VarMix_init
use MOM_lateral_mixing_coeffs, only : calc_resoln_function, VarMix_CS
use MOM_interface_heights, only : find_eta
use MOM_MEKE, only : MEKE_init, MEKE_alloc_register_restart, step_forward_MEKE, MEKE_CS
use MOM_MEKE_types, only : MEKE_type
use MOM_mixed_layer_restrat, only : mixedlayer_restrat, mixedlayer_restrat_init, mixedlayer_restrat_CS
use MOM_open_boundary, only : Radiation_Open_Bdry_Conds, open_boundary_init
use MOM_open_boundary, only : open_boundary_CS
use MOM_PressureForce, only : PressureForce, PressureForce_init, PressureForce_CS
use MOM_thickness_diffuse, only : thickness_diffuse, thickness_diffuse_init, thickness_diffuse_CS
use MOM_tidal_forcing, only : tidal_forcing_init, tidal_forcing_CS
use MOM_tracer, only : advect_tracer, register_tracer, add_tracer_diagnostics
use MOM_tracer, only : add_tracer_2d_diagnostics, tracer_hordiff
use MOM_tracer, only : advect_tracer_init, advect_tracer_diag_init, advect_tracer_CS
use MOM_tracer_flow_control, only : call_tracer_register, tracer_flow_control_CS
use MOM_tracer_flow_control, only : tracer_flow_control_init, call_tracer_surface_state
use MOM_vert_friction, only : vertvisc, vertvisc_coef, vertvisc_remnant
use MOM_vert_friction, only : vertvisc_limit_vel, vertvisc_init, vertvisc_CS
use MOM_set_visc, only : set_viscous_BBL, set_viscous_ML, set_visc_init, set_visc_CS
use MOM_CS_type, only : MOM_control_struct
use MOM_CS_type, only : MOM_state_chksum, MOM_thermo_chksum, MOM_accel_chksum

implicit none ; private

#include <MOM_memory.h>

public step_MOM_dyn_unsplit, register_restarts_dyn_unsplit
public initialize_dyn_unsplit

integer :: id_clock_Cor, id_clock_pres, id_clock_vertvisc
integer :: id_clock_horvisc, id_clock_mom_update
integer :: id_clock_continuity, id_clock_thick_diff, id_clock_ml_restrat
integer :: id_clock_pass, id_clock_pass_init

contains

! =============================================================================

subroutine step_MOM_dyn_unsplit(u, v, h, Time_local, dt, fluxes, &
                  p_surf_begin, p_surf_end, uh, vh, eta_av, G, CS)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(inout) :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(inout) :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(inout) :: h
  type(time_type),                        intent(in)    :: Time_local
  real,                                   intent(in)    :: dt
  type(forcing),                          intent(in)    :: fluxes
  real, dimension(:,:),                   pointer       :: p_surf_begin, p_surf_end
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(inout) :: uh
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(inout) :: vh
  real, dimension(NIMEM_,NJMEM_),         intent(out)   :: eta_av
  type(ocean_grid_type),                  intent(inout) :: G
  type(MOM_control_struct),               pointer      :: CS
! Arguments: u - The input and output zonal velocity, in m s-1.
!  (in)      v - The input and output meridional velocity, in m s-1.
!  (in)      h - The input and output layer thicknesses, in m or kg m-2,
!                depending on whether the Boussinesq approximation is made.
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
!  (out)     eta_av - The time-mean free surface height or column mass, in m or
!                     kg m-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure set up by initialize_MOM.

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_av, hp
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: up, upp
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: vp, vpp
  real, dimension(:,:), pointer :: p_surf
  real :: dt_pred   ! The time step for the predictor part of the baroclinic
                    ! time stepping.
  logical :: dyn_p_surf
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  dt_pred = dt / 3.0

  h_av(:,:,:) = 0; hp(:,:,:) = 0
  up(:,:,:) = 0; upp(:,:,:) = 0
  vp(:,:,:) = 0; vpp(:,:,:) = 0

  dyn_p_surf = CS%interp_p_surf .and. associated(p_surf_begin) .and. &
               associated(p_surf_end)
  if (dyn_p_surf) then
    call safe_alloc_ptr(p_surf,G%isd,G%ied,G%jsd,G%jed) ; p_surf(:,:) = 0.0
  else
    p_surf => fluxes%p_surf
  endif

! Matsuno's third order accurate three step scheme is used to step
! all of the fields except h.  h is stepped separately.

  if (CS%debug) then
    call MOM_state_chksum("Start First Predictor ", u, v, h, uh, vh, G)
  endif

! diffu = horizontal viscosity terms (u,h)
  call enable_averaging(dt,Time_local, CS%diag)
  call cpu_clock_begin(id_clock_horvisc)
  call horizontal_viscosity(u, v, h, CS%diffu, CS%diffv, CS%MEKE, CS%Varmix, &
                            G, CS%hor_visc_CSp)
  call cpu_clock_end(id_clock_horvisc)
  call disable_averaging(CS%diag)

! uh = u*h
! hp = h + dt/2 div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(u, v, h, hp, uh, vh, dt*0.5, G, CS%continuity_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_continuity)
  call cpu_clock_begin(id_clock_pass)
  call pass_var(hp, G%Domain)
  call pass_vector(uh, vh, G%Domain)
  call cpu_clock_end(id_clock_pass)

  call enable_averaging(0.5*dt,Time_local-set_time(int(0.5*dt)), CS%diag)
!   Here the thickness fluxes are offered for averaging.
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
      u(i,j,k) = u(i,j,k) + dt * CS%diffu(i,j,k) * G%umask(i,j)
    enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie
      v(i,j,k) = v(i,j,k) + dt * CS%diffv(i,j,k) * G%vmask(i,j)
    enddo ; enddo
    do j=js-2,je+2 ; do I=Isq-2,Ieq+2
      CS%uhtr(i,j,k) = CS%uhtr(i,j,k) + 0.5*dt*uh(i,j,k)
    enddo ; enddo
    do J=Jsq-2,Jeq+2 ; do i=is-2,ie+2
      CS%vhtr(i,j,k) = CS%vhtr(i,j,k) + 0.5*dt*vh(i,j,k)
    enddo ; enddo
  enddo
  call cpu_clock_end(id_clock_mom_update)
  call cpu_clock_begin(id_clock_pass)
  call pass_vector(u, v, G%Domain)
  call cpu_clock_end(id_clock_pass)

! CAu = -(f+zeta)/h_av vh + d/dx KE
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(u, v, h_av, uh, vh, CS%CAu, CS%CAv, G, CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! PFu = d/dx M(h_av,T,S)
  call cpu_clock_begin(id_clock_pres)
  if (dyn_p_surf) then ; do j=js-2,je+2 ; do i=is-2,ie+2
    p_surf(i,j) = 0.75*p_surf_begin(i,j) + 0.25*p_surf_end(i,j)
  enddo ; enddo ; endif
  call PressureForce(h_av, CS%tv, CS%PFu, CS%PFv, G, &
                     CS%PressureForce_CSp, CS%regridding_opts, p_surf)
  call cpu_clock_end(id_clock_pres)

! up = u + dt_pred * (PFu + CAu)
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    up(i,j,k) = G%umask(i,j) * (u(i,j,k) + dt_pred * &
                               (CS%PFu(i,j,k) + CS%CAu(i,j,k)))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    vp(i,j,k) = G%vmask(i,j) * (v(i,j,k) + dt_pred * &
                               (CS%PFv(i,j,k) + CS%CAv(i,j,k)))
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) then
    call MOM_state_chksum("Predictor 1", up, vp, h_av, uh, vh, G)
    call MOM_accel_chksum("Predictor 1 accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv,&
                          CS%diffu, CS%diffv, G)
  endif

! CS%visc contains viscosity and BBL thickness (u_in,h_in)
  if (CS%calc_bbl) then
    call enable_averaging(CS%bbl_calc_time_interval, &
                          Time_local-set_time(int(dt)), CS%diag)
    call set_viscous_BBL(u, v, h_av, CS%tv, CS%visc, G, CS%set_visc_CSp)
    call cpu_clock_begin(id_clock_pass)
    if (associated(CS%visc%Ray_u) .and. associated(CS%visc%Ray_v)) &
      call pass_vector(CS%visc%Ray_u, CS%visc%Ray_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE)
    if (associated(CS%visc%kv_bbl_u) .and. associated(CS%visc%kv_bbl_v)) then
      call pass_vector(CS%visc%bbl_thick_u, CS%visc%bbl_thick_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE, complete=.false.)
      call pass_vector(CS%visc%kv_bbl_u, CS%visc%kv_bbl_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE)
    endif
    call cpu_clock_end(id_clock_pass)
    call disable_averaging(CS%diag)
    CS%calc_bbl = .false.
  endif

 ! up <- up + dt/2 d/dz visc d/dz up
  call cpu_clock_begin(id_clock_vertvisc)
  call enable_averaging(dt, Time_local, CS%diag)
  call set_viscous_ML(u, v, h_av, CS%tv, fluxes, CS%visc, dt*0.5, G, &
                      CS%set_visc_CSp)
  call disable_averaging(CS%diag)
  call vertvisc_coef(up, vp, h_av, fluxes, CS%visc, dt*0.5, G, CS%vertvisc_CSp)
  call vertvisc(up, vp, h_av, fluxes, CS%visc, dt*0.5, CS%OBC, G, &
                CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)
  call cpu_clock_begin(id_clock_pass)
  call pass_vector(up, vp, G%Domain)
  call cpu_clock_end(id_clock_pass)

! uh = up * hp
! h_av = hp + dt/2 div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(up, vp, hp, h_av, uh, vh, &
                  (0.5*dt), G, CS%continuity_CSp, OBC=CS%OBC)
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
  call CorAdCalc(up, vp, h_av, uh, vh, CS%CAu, CS%CAv, &
                 G, CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! PFu = d/dx M(h_av,T,S)
  call cpu_clock_begin(id_clock_pres)
  if (dyn_p_surf) then ; do j=js-2,je+2 ; do i=is-2,ie+2
    p_surf(i,j) = 0.25*p_surf_begin(i,j) + 0.75*p_surf_end(i,j)
  enddo ; enddo ; endif
  call PressureForce(h_av, CS%tv, CS%PFu, CS%PFv, G, &
                     CS%PressureForce_CSp, CS%regridding_opts, p_surf)
  call cpu_clock_end(id_clock_pres)

! upp = u + dt/2 * ( PFu + CAu )
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    upp(i,j,k) = G%umask(i,j) * (u(i,j,k) + dt * 0.5 * &
            (CS%PFu(i,j,k) + CS%CAu(i,j,k)))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    vpp(i,j,k) = G%vmask(i,j) * (v(i,j,k) + dt * 0.5 * &
            (CS%PFv(i,j,k) + CS%CAv(i,j,k)))
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) then
    call MOM_state_chksum("Predictor 2", upp, vpp, h_av, uh, vh, G)
    call MOM_accel_chksum("Predictor 2 accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv,&
                          CS%diffu, CS%diffv, G)
  endif

! upp <- upp + dt/2 d/dz visc d/dz upp
  call cpu_clock_begin(id_clock_vertvisc)
  call vertvisc_coef(upp, vpp, hp, fluxes, CS%visc, dt*0.5, G, CS%vertvisc_CSp)
  call vertvisc(upp, vpp, hp, fluxes, CS%visc, dt*0.5, CS%OBC, G, &
                CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)
  call cpu_clock_begin(id_clock_pass)
  call pass_vector(upp, vpp, G%Domain)
  call cpu_clock_end(id_clock_pass)

! uh = upp * hp
! h = hp + dt/2 div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(upp, vpp, hp, h, uh, vh, &
                  (dt*0.5), G, CS%continuity_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_continuity)
  call cpu_clock_begin(id_clock_pass)
  call pass_var(h, G%Domain)
  call pass_vector(uh, vh, G%Domain)
  call cpu_clock_end(id_clock_pass)

  call enable_averaging(0.5*dt, Time_local, CS%diag)
!   Here the thickness fluxes are offered for averaging.
  if (CS%id_uh > 0) call post_data(CS%id_uh, uh, CS%diag)
  if (CS%id_vh > 0) call post_data(CS%id_vh, vh, CS%diag)
  call disable_averaging(CS%diag)

! h_av = (h + hp)/2
  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      h_av(i,j,k) = 0.5*(h(i,j,k) + hp(i,j,k))
    enddo ; enddo
    do j=js-2,je+2 ; do I=Isq-2,Ieq+2
      CS%uhtr(i,j,k) = CS%uhtr(i,j,k) + 0.5*dt*uh(i,j,k)
    enddo ; enddo
    do J=Jsq-2,Jeq+2 ; do i=is-2,ie+2
      CS%vhtr(i,j,k) = CS%vhtr(i,j,k) + 0.5*dt*vh(i,j,k)
    enddo ; enddo
  enddo

  call enable_averaging(dt,Time_local, CS%diag)

! CAu = -(f+zeta(upp))/h_av vh + d/dx KE(upp)
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(upp, vpp, h_av, uh, vh, CS%CAu, CS%CAv, &
                 G, CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! PFu = d/dx M(h_av,T,S)
  call cpu_clock_begin(id_clock_pres)
  call PressureForce(h_av, CS%tv, CS%PFu, CS%PFv, G, &
                     CS%PressureForce_CSp, CS%regridding_opts, p_surf)
  call cpu_clock_end(id_clock_pres)

! u = u + dt * ( PFu + CAu )
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u(i,j,k) = G%umask(i,j) * (u(i,j,k) + dt * &
            (CS%PFu(i,j,k) + CS%CAu(i,j,k)))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v(i,j,k) = G%vmask(i,j) * (v(i,j,k) + dt * &
            (CS%PFv(i,j,k) + CS%CAv(i,j,k)))
  enddo ; enddo ; enddo

! u <- u + dt d/dz visc d/dz u
  call cpu_clock_begin(id_clock_vertvisc)
  call vertvisc_coef(u, v, h_av, fluxes, CS%visc, dt, G, CS%vertvisc_CSp)
  call vertvisc(u, v, h_av, fluxes, CS%visc, dt, CS%OBC, G, &
                CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot)
  call cpu_clock_end(id_clock_vertvisc)
  call cpu_clock_begin(id_clock_pass)
  call pass_vector(u, v, G%Domain)
  call cpu_clock_end(id_clock_pass)

  if (CS%thickness_diffuse .and. .not.CS%thickness_diffuse_first) then
    call cpu_clock_begin(id_clock_thick_diff)
    if (associated(CS%VarMix)) &
      call calc_slope_function(h, CS%tv, G, CS%VarMix)
    call thickness_diffuse(h, CS%uhtr, CS%vhtr, CS%tv, dt, G, &
                           CS%MEKE, CS%VarMix, CS%thickness_diffuse_CSp)
    call cpu_clock_end(id_clock_thick_diff)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(h, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  if (CS%mixedlayer_restrat) then
    call cpu_clock_begin(id_clock_ml_restrat)
    call mixedlayer_restrat(h, CS%uhtr ,CS%vhtr, CS%tv, fluxes, dt, &
                            G, CS%mixedlayer_restrat_CSp)
    call cpu_clock_end(id_clock_ml_restrat)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(h, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  if (associated(CS%MEKE)) then
    call step_forward_MEKE(CS%MEKE, h, CS%visc, dt, G, CS%MEKE_CSp)
  endif

  if (CS%debug) then
    call MOM_state_chksum("Corrector", u, v, h, uh, vh, G)
    call MOM_accel_chksum("Corrector accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv, &
                          CS%diffu, CS%diffv, G)
  endif

  if (G%Boussinesq) then
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

subroutine register_restarts_dyn_unsplit(G, param_file, CS, restart_CS)
  type(ocean_grid_type),     intent(in) :: G
  type(param_file_type),     intent(in) :: param_file
  type(MOM_control_struct), intent(in) :: CS
  type(MOM_restart_CS),     pointer    :: restart_CS
!   This subroutine sets up any auxiliary restart variables that are specific
! to the unsplit time stepping scheme.  All variables registered here should
! have the ability to be recreated if they are not present in a restart file.

! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      CS - The control structure set up by initialize_MOM.
!  (in)      restart_CS - A pointer to the restart control structure.

  type(vardesc) :: vd
  character(len=48) :: thickness_units, flux_units

! This is where a control structure that is specific to this module would be allocated.
! if (associated(CS_unsplit)) then
!   call MOM_error(WARNING, "register_restarts_unsplit called with an associated "// &
!                            "control structure.")
!   return
! endif
! allocate(CS_unsplit)

  thickness_units = get_thickness_units(G)
  flux_units = get_flux_units(G)

!  No extra restart fields are needed with this time stepping scheme.

end subroutine register_restarts_dyn_unsplit

subroutine initialize_dyn_unsplit(u, v, h, Time, G, param_file, diag, CS, restart_CS)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(inout) :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(inout) :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_) , intent(inout) :: h
  type(time_type),                target, intent(in)    :: Time
  type(ocean_grid_type),                  intent(inout) :: G
  type(param_file_type),                  intent(in)    :: param_file
  type(diag_ptrs),                target, intent(inout) :: diag
  type(MOM_control_struct),               intent(inout) :: CS
  type(MOM_restart_CS),                   pointer       :: restart_CS

  !   This subroutine initializes any variables that are specific to this time
  ! stepping scheme, including the cpu clocks.

  id_clock_Cor = cpu_clock_id('(Ocean Coriolis & mom advection)', grain=CLOCK_MODULE)
  id_clock_continuity = cpu_clock_id('(Ocean continuity equation)', grain=CLOCK_MODULE)
  id_clock_pres = cpu_clock_id('(Ocean pressure force)', grain=CLOCK_MODULE)
  id_clock_vertvisc = cpu_clock_id('(Ocean vertical viscosity)', grain=CLOCK_MODULE)
  id_clock_horvisc = cpu_clock_id('(Ocean horizontal viscosity)', grain=CLOCK_MODULE)
  id_clock_mom_update = cpu_clock_id('(Ocean momentum increments)', grain=CLOCK_MODULE)
  id_clock_pass = cpu_clock_id('(Ocean message passing)', grain=CLOCK_MODULE)
  id_clock_pass_init = cpu_clock_id('(Ocean init message passing)', grain=CLOCK_ROUTINE)
  if (CS%thickness_diffuse) &
    id_clock_thick_diff = cpu_clock_id('(Ocean thickness diffusion)', grain=CLOCK_MODULE)
  if (CS%mixedlayer_restrat) &
    id_clock_ml_restrat = cpu_clock_id('(Ocean mixed layer restrat)', grain=CLOCK_MODULE)

end subroutine initialize_dyn_unsplit

end module MOM_dynamics_unsplit
