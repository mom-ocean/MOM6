module MOM_dynamics_split_RK2

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
!*  By Robert Hallberg and Alistair Adcroft, 1996-2012                 *
!*                                                                     *
!*    This file contains code that does the time-stepping of the       *
!*  adiabatic dynamic core, in this case with mode-splitting between   *
!*  the baroclinic and barotropic modes and a pseudo-second order      *
!*  Runge-Kutta time stepping scheme for the baroclinic momentum       *
!*  equation and a forward-backward coupling between the baroclinic    *
!*  momentum and continuity equations.  This split time-stepping       *
!*  scheme is described in detail in Hallberg (JCP, 1997), with the    *
!*  primary issues related to exact tracer conservation and how to     *
!*  ensure consistency between the barotropic and layered estimates    *
!*  of the free surface height described carefully in Hallberg and     *
!*  Adcroft (Ocean Modelling, 2009).  This was the time stepping code  *
!*  that is used for most GOLD applications, including GFDL's ESM2G    *
!*  Earth system model, and all of the examples provided with the      *
!*  MOM code (although several of these solutions are routinely        *
!*  verified by comparison with the slower unsplit schemes).           *
!*                                                                     *
!*    The subroutine step_MOM_dyn_split_RK2 actually does the time     *
!*  stepping, while register_restarts_dyn_split_RK2 sets the fields    *
!*  that are found in a full restart file with this scheme, and        *
!*  initialize_dyn_split_RK2 initializes the cpu clocks that are       *                                      *
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


use MOM_variables, only : vertvisc_type, ocean_OBC_type, thermo_var_ptrs
use MOM_variables, only : BT_cont_type, alloc_bt_cont_type, dealloc_bt_cont_type
use MOM_variables, only : ocean_internal_state
use MOM_forcing_type, only : forcing

use MOM_checksum_packages, only : MOM_thermo_chksum, MOM_state_chksum, MOM_accel_chksum
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
use MOM_file_parser, only : read_param, get_param, log_version, param_file_type
use MOM_io, only : MOM_io_init, vardesc
use MOM_restart, only : register_restart_field, query_initialized, save_restart
use MOM_restart, only : restart_init, MOM_restart_CS
use MOM_time_manager, only : time_type, set_time, time_type_to_real, operator(+)
use MOM_time_manager, only : operator(-), operator(>), operator(*), operator(/)

use MOM_barotropic, only : barotropic_init, btstep, btcalc, bt_mass_source
use MOM_barotropic, only : register_barotropic_restarts, set_dtbt, barotropic_CS
use MOM_continuity, only : continuity, continuity_init, continuity_CS
use MOM_CoriolisAdv, only : CorAdCalc, CoriolisAdv_init, CoriolisAdv_CS
use MOM_diabatic_driver, only : diabatic, diabatic_driver_init, diabatic_CS
use MOM_EOS, only : select_eqn_of_state
use MOM_error_checking, only : check_redundant
use MOM_grid, only : MOM_grid_init, ocean_grid_type, get_thickness_units
use MOM_grid, only : get_flux_units, get_tr_flux_units
use MOM_hor_visc, only : horizontal_viscosity, hor_visc_init, hor_visc_CS
use MOM_interface_heights, only : find_eta
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE_types, only : MEKE_type
use MOM_open_boundary, only : Radiation_Open_Bdry_Conds, open_boundary_init
use MOM_open_boundary, only : open_boundary_CS
use MOM_PressureForce, only : PressureForce, PressureForce_init, PressureForce_CS
use MOM_tidal_forcing, only : tidal_forcing_init, tidal_forcing_CS
use MOM_vert_friction, only : vertvisc, vertvisc_coef, vertvisc_remnant
use MOM_vert_friction, only : vertvisc_limit_vel, vertvisc_init, vertvisc_CS
use MOM_set_visc, only : set_viscous_BBL, set_viscous_ML, set_visc_init, set_visc_CS
use MOM_CS_type, only : MOM_control_struct, MOM_dyn_control_struct

implicit none ; private

#include <MOM_memory.h>

public step_MOM_dyn_split_RK2, register_restarts_dyn_split_RK2
public adjustments_dyn_split_RK2
public initialize_dyn_split_RK2, end_dyn_split_RK2

integer :: id_clock_Cor, id_clock_pres, id_clock_vertvisc
integer :: id_clock_horvisc, id_clock_mom_update
integer :: id_clock_continuity, id_clock_thick_diff
integer :: id_clock_btstep, id_clock_btcalc, id_clock_btforce
integer :: id_clock_pass, id_clock_pass_init

contains

! =============================================================================

subroutine step_MOM_dyn_split_RK2(u, v, h, tv, visc, &
                 Time_local, dt, fluxes, p_surf_begin, p_surf_end, &
                 dt_since_flux, dt_therm, uh, vh, uhtr, vhtr, eta_av, &
                 G, CS, calc_dtbt, VarMix, MEKE)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), target, intent(inout) :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), target, intent(inout) :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(inout) :: h
  type(thermo_var_ptrs),                  intent(in)    :: tv
  type(vertvisc_type),                    intent(inout) :: visc
  type(time_type),                        intent(in)    :: Time_local
  real,                                   intent(in)    :: dt
  type(forcing),                          intent(in)    :: fluxes
  real, dimension(:,:),                   pointer       :: p_surf_begin, p_surf_end
  real,                                   intent(in)    :: dt_since_flux, dt_therm
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), target, intent(inout) :: uh
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), target, intent(inout) :: vh
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(inout) :: uhtr
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(inout) :: vhtr
  real, dimension(NIMEM_,NJMEM_),         intent(out)   :: eta_av
  type(ocean_grid_type),                  intent(inout) :: G
  type(MOM_dyn_control_struct),           pointer       :: CS
  logical,                                intent(in)    :: calc_dtbt
  type(VarMix_CS),                        pointer       :: VarMix
  type(MEKE_type),                        pointer       :: MEKE
! Arguments: u - The zonal velocity, in m s-1.
!  (inout)   v - The meridional velocity, in m s-1.
!  (inout)   h - The layer thicknesses, in m or kg m-2, depending on
!                whether the Boussinesq approximation is made.
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
!  (in)      dt_since_flux - The elapsed time since fluxes were applied, in s.
!  (in)      dt_therm - The thermodynamic time step, in s.
!  (inout)   uh - The zonal volume or mass transport, in m3 s-1 or kg s-1.
!  (inout)   vh - The meridional volume or mass transport, in m3 s-1 or kg s-1.
!  (inout)   uhtr - The accumulated zonal volume or mass transport since the last
!                   tracer advection, in m3 or kg.
!  (inout)   vhtr - The accumulated meridional volume or mass transport since the last
!                   tracer advection, in m3 or kg.
!  (out)     eta_av - The free surface height or column mass time-averaged
!                     over a time step, in m or kg m-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure set up by initialize_MOM.
!  (in)      calc_dtbt - If true, recalculate the barotropic time step.
!  (in)      VarMix - A pointer to a structure with fields that specify the
!                     spatially variable viscosities.
!  (inout)   MEKE - A pointer to a structure containing fields related to
!                   the Mesoscale Eddy Kinetic Energy.

  real :: dt_pred   ! The time step for the predictor part of the baroclinic
                    ! time stepping.

  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: &
    up   ! Predicted zonal velocitiy in m s-1.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: &
    vp   ! Predicted meridional velocitiy in m s-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G))  :: &
    hp   ! Predicted thickness in m or kg m-2 (H).

  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: u_bc_accel
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: v_bc_accel
    ! u_bc_accel and v_bc_accel are the summed baroclinic accelerations of each
    ! layer calculated by the non-barotropic part of the model, both in m s-2.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), target :: uh_in
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), target :: vh_in
    ! uh_in and vh_in are the zonal or meridional mass transports that would be
    ! obtained using the initial velocities, both in m3 s-1 or kg s-1.
  real, dimension(SZIB_(G),SZJ_(G)) :: uhbt_out
  real, dimension(SZI_(G),SZJB_(G)) :: vhbt_out
    ! uhbt_out and vhbt_out are the vertically summed transports from the
    ! barotropic solver based on its final velocities, both in m3 s-1 or kg s-1.
  real, dimension(SZI_(G),SZJ_(G)) :: eta_pred
    ! eta_pred is the predictor value of the free surface height or column mass,
    ! in m or kg m-2.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), target :: u_adj
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), target :: v_adj
    ! u_adj and v_adj are the zonal or meridional velocities after u and v
    ! have been barotropically adjusted so the resulting transports match
    ! uhbt_out and vhbt_out, both in m s-1.

  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: u_old_rad_OBC
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: v_old_rad_OBC
  real, dimension(SZI_(G),SZJ_(G),SZK_(G))  :: h_old_rad_OBC
    ! u_old_rad_OBC and v_old_rad_OBC are the starting velocities, which are
    ! saved for use in the Flather open boundary condition code, both in m s-1.
  
  real :: Pa_to_eta ! A factor that converts pressures to the units of eta.
  real, pointer, dimension(:,:) :: &
    p_surf => NULL(), eta_PF_start => NULL(), &
    taux_bot => NULL(), tauy_bot => NULL(), &
    uhbt_in, vhbt_in, eta
  real, pointer, dimension(:,:,:) :: &
    uh_ptr => NULL(), u_ptr => NULL(),  vh_ptr => NULL(), v_ptr => NULL(), &
    u_init => NULL(), v_init => NULL(), & ! Pointers to u and v or u_adj and v_adj.
    u_av, & ! The zonal velocity time-averaged over a time step, in m s-1.
    v_av, & ! The meridional velocity time-averaged over a time step, in m s-1.
    h_av    ! The layer thickness time-averaged over a time step, in m or
            ! kg m-2.
  real :: Idt
  logical :: dyn_p_surf
  logical :: BT_cont_BT_thick ! If true, use the BT_cont_type to estimate the
                              ! relative weightings of the layers in calculating
                              ! the barotropic accelerations.
  integer :: pid_Ray, pid_bbl_h, pid_kv_bbl, pid_eta_PF, pid_eta, pid_visc
  integer :: pid_h, pid_u, pid_u_av, pid_uh, pid_uhbt_in
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  u_av => CS%u_av ; v_av => CS%v_av ; h_av => CS%h_av
  eta => CS%eta ; uhbt_in => CS%uhbt_in ; vhbt_in => CS%vhbt_in
  Idt = 1.0 / dt

  up(:,:,:) = 0.0 ; vp(:,:,:) = 0.0 ; hp(:,:,:) = h(:,:,:)

  if (CS%debug) then
    call MOM_state_chksum("Start predictor ", u, v, h, uh, vh, G)
    call check_redundant("Start predictor u ", u, v, G)
    call check_redundant("Start predictor uh ", uh, vh, G)
  endif

  dyn_p_surf = CS%interp_p_surf .and. associated(p_surf_begin) .and. &
               associated(p_surf_end)
  if (dyn_p_surf) then
    p_surf => p_surf_end
    call safe_alloc_ptr(eta_PF_start,G%isd,G%ied,G%jsd,G%jed)
    eta_PF_start(:,:) = 0.0
  else
    p_surf => fluxes%p_surf
  endif

  if (associated(CS%OBC)) then
    do k=1,nz ; do j=js,je ; do I=is-2,ie+1
      u_old_rad_OBC(I,j,k) = u(I,j,k)
    enddo ; enddo ; enddo
    do k=1,nz ; do J=js-2,je+1 ; do i=is,ie
      v_old_rad_OBC(i,J,k) = v(i,J,k)
    enddo ; enddo ; enddo
    do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      h_old_rad_OBC(i,j,k) = h(i,j,k)
    enddo ; enddo ; enddo
  endif

  BT_cont_BT_thick = .false.
  if (associated(CS%BT_cont)) BT_cont_BT_thick = &
    (associated(CS%BT_cont%h_u) .and. associated(CS%BT_cont%h_v))

  if (CS%split_bottom_stress) then
    taux_bot => CS%taux_bot ; tauy_bot => CS%tauy_bot
  endif

  if (visc%calc_bbl) then
    ! Calculate the BBL properties and store them inside visc (u,h).
    call cpu_clock_begin(id_clock_vertvisc)
    call enable_averaging(visc%bbl_calc_time_interval, &
                          Time_local-set_time(int(dt)), CS%diag)
    call set_viscous_BBL(u, v, h, tv, visc, G, CS%set_visc_CSp)
    call disable_averaging(CS%diag)
    call cpu_clock_end(id_clock_vertvisc)

    call cpu_clock_begin(id_clock_pass)
    if (G%nonblocking_updates) then   
      if (associated(visc%Ray_u) .and. associated(visc%Ray_v)) &
        pid_Ray = pass_vector_start(visc%Ray_u, visc%Ray_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(visc%bbl_thick_u) .and. associated(visc%bbl_thick_v)) &
        pid_bbl_h = pass_vector_start(visc%bbl_thick_u, visc%bbl_thick_v, &
                       G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(visc%kv_bbl_u) .and. associated(visc%kv_bbl_v)) &
        pid_kv_bbl = pass_vector_start(visc%kv_bbl_u, visc%kv_bbl_v, &
                       G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
      ! visc%calc_bbl will be set to .false. when the message passing is complete.
    else
      if (associated(visc%Ray_u) .and. associated(visc%Ray_v)) &
        call pass_vector(visc%Ray_u, visc%Ray_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(visc%kv_bbl_u) .and. associated(visc%kv_bbl_v)) then
        call pass_vector(visc%bbl_thick_u, visc%bbl_thick_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE, complete=.false.)
        call pass_vector(visc%kv_bbl_u, visc%kv_bbl_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE)
      endif
      visc%calc_bbl = .false.
    endif
    call cpu_clock_end(id_clock_pass)
  endif

! PFu = d/dx M(h,T,S)
! pbce = dM/deta
  if (CS%begw == 0.0) call enable_averaging(dt, Time_local, CS%diag)
  call cpu_clock_begin(id_clock_pres)
  call PressureForce(h, tv, CS%PFu, CS%PFv, G, CS%PressureForce_CSp, &
                     CS%ALE_CSp, p_surf, CS%pbce, CS%eta_PF)
  if (dyn_p_surf) then
    if (G%Boussinesq) then
      Pa_to_eta = 1.0 / (G%Rho0*G%g_Earth)
    else
      Pa_to_eta = 1.0 / G%H_to_Pa
    endif
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_PF_start(i,j) = CS%eta_PF(i,j) - Pa_to_eta * &
                          (p_surf_begin(i,j) - p_surf_end(i,j))
    enddo ; enddo
  endif
  call cpu_clock_end(id_clock_pres)
  call disable_averaging(CS%diag)

  if (G%nonblocking_updates) then
    call cpu_clock_begin(id_clock_pass)
    pid_eta_PF = pass_var_start(CS%eta_PF, G%Domain)
    pid_eta = pass_var_start(eta, G%Domain)
    if (CS%readjust_velocity) &
      pid_uhbt_in = pass_vector_start(uhbt_in, vhbt_in, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

! CAu = -(f+zeta_av)/h_av vh + d/dx KE_av
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(u_av, v_av, h_av, uh, vh, CS%CAu, CS%CAv, G,CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! u_bc_accel = CAu + PFu + diffu(u[n-1])
  call cpu_clock_begin(id_clock_btforce)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u_bc_accel(I,j,k) = (CS%Cau(I,j,k) + CS%PFu(I,j,k)) + CS%diffu(I,j,k)
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v_bc_accel(i,J,k) = (CS%Cav(i,J,k) + CS%PFv(i,J,k)) + CS%diffv(i,J,k)
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_btforce)

  if (CS%debug) then
    call check_redundant("pre-btstep CS%Ca ", CS%Cau, CS%Cav, G)
    call check_redundant("pre-btstep CS%PF ", CS%PFu, CS%PFv, G)
    call check_redundant("pre-btstep CS%diff ", CS%diffu, CS%diffv, G)
    call check_redundant("pre-btstep u_bc_accel ", u_bc_accel, v_bc_accel, G)
  endif

  if (G%nonblocking_updates) then
    call cpu_clock_begin(id_clock_pass)
    if (visc%calc_bbl) then
      if (associated(visc%Ray_u) .and. associated(visc%Ray_v)) &
        call pass_vector_complete(pid_Ray, visc%Ray_u, visc%Ray_v, G%Domain, &
                         To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(visc%bbl_thick_u) .and. associated(visc%bbl_thick_v)) &
        call pass_vector_complete(pid_bbl_h, visc%bbl_thick_u, visc%bbl_thick_v, &
                       G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(visc%kv_bbl_u) .and. associated(visc%kv_bbl_v)) &
        call pass_vector_complete(pid_kv_bbl, visc%kv_bbl_u, visc%kv_bbl_v, &
                       G%Domain, To_All+SCALAR_PAIR, CGRID_NE)

      ! visc%calc_bbl is set to .false. now that the message passing is completed.
      visc%calc_bbl = .false.
    endif
    call pass_var_complete(pid_eta_PF, CS%eta_PF, G%Domain)
    call pass_var_complete(pid_eta, eta, G%Domain)
    if (CS%readjust_velocity) &
      call pass_vector_complete(pid_uhbt_in, uhbt_in, vhbt_in, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  call cpu_clock_begin(id_clock_vertvisc)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    up(i,j,k) = G%mask2dCu(i,j) * (u(i,j,k) + dt * u_bc_accel(I,j,k))
  enddo ; enddo ;  enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    vp(i,j,k) = G%mask2dCv(i,j) * (v(i,j,k) + dt * v_bc_accel(i,J,k))
  enddo ; enddo ;  enddo
  call enable_averaging(dt, Time_local, CS%diag)
  call set_viscous_ML(u, v, h, tv, fluxes, visc, dt, G, &
                      CS%set_visc_CSp)
  call disable_averaging(CS%diag)

  call vertvisc_coef(up, vp, h, fluxes, visc, dt, G, CS%vertvisc_CSp)
  call vertvisc_remnant(visc, CS%visc_rem_u, CS%visc_rem_v, dt, G, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)

  call cpu_clock_begin(id_clock_pass)
  if (G%nonblocking_updates) then
    pid_visc = pass_vector_start(CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                                 To_All+SCALAR_PAIR, CGRID_NE)
  else
    call pass_var(CS%eta_PF, G%Domain, complete=.false.)
    call pass_var(eta, G%Domain)
    if (CS%readjust_velocity) call pass_vector(uhbt_in, vhbt_in, G%Domain)
    call pass_vector(CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE)
  endif
  call cpu_clock_end(id_clock_pass)

  call cpu_clock_begin(id_clock_btcalc)
  ! Calculate the relative layer weights for determining barotropic quantities.
  if (.not.BT_cont_BT_thick) &
    call btcalc(h, G, CS%barotropic_CSp)
  call bt_mass_source(h, eta, fluxes, .true., dt_therm, dt_since_flux, &
                      G, CS%barotropic_CSp)
  call cpu_clock_end(id_clock_btcalc)

  if (G%nonblocking_updates) then
    call cpu_clock_begin(id_clock_pass)
    call pass_vector_complete(pid_visc, CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE)
    call cpu_clock_end(id_clock_pass)
  endif

! u_accel_bt = layer accelerations due to barotropic solver
  if (CS%flux_BT_coupling) then
    call cpu_clock_begin(id_clock_continuity)
    if (CS%readjust_velocity) then
      ! Adjust the input velocites so that their transports match uhbt_out & vhbt_out.
      call continuity(u, v, h, hp, uh_in, vh_in, dt, G, &
                      CS%continuity_CSp, uhbt_in, vhbt_in, CS%OBC, &
                      CS%visc_rem_u, CS%visc_rem_v, u_adj, v_adj, &
                      BT_cont=CS%BT_cont)
      u_init => u_adj ; v_init => v_adj
      if (ASSOCIATED(CS%diag%du_adj2)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
        CS%diag%du_adj2(I,j,k) = u_adj(I,j,k) - u(I,j,k)
      enddo ; enddo ; enddo ; endif
      if (ASSOCIATED(CS%diag%dv_adj2)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
        CS%diag%dv_adj2(i,J,k) = v_adj(i,J,k) - v(i,J,k)
      enddo ; enddo ; enddo ; endif
      CS%readjust_velocity = .false.
    else
      call continuity(u, v, h, hp, uh_in, vh_in, dt, G, &
                      CS%continuity_CSp, OBC=CS%OBC, BT_cont=CS%BT_cont)
!###   call continuity(u, v, h, hp, uh_in, vh_in, dt, G, &
!###                   CS%continuity_CSp, OBC=CS%OBC, visc_rem_u=CS%visc_rem_u, &
!###                      visc_rem_v=CS%visc_rem_v, BT_cont=CS%BT_cont)
      u_init => u ; v_init => v
    endif
    call cpu_clock_end(id_clock_continuity)

    if (BT_cont_BT_thick) then
      call cpu_clock_begin(id_clock_pass)
      call pass_vector(CS%BT_cont%h_u, CS%BT_cont%h_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE)
      call cpu_clock_end(id_clock_pass)
      call btcalc(h, G, CS%barotropic_CSp, CS%BT_cont%h_u, CS%BT_cont%h_v)
    endif
    call cpu_clock_begin(id_clock_btstep)
    if (calc_dtbt) call set_dtbt(G, CS%barotropic_CSp, eta, CS%pbce, CS%BT_cont)
    call btstep(.true., uh_in, vh_in, eta, dt, u_bc_accel, v_bc_accel, &
                fluxes, CS%pbce, CS%eta_PF, uh, vh, CS%u_accel_bt, &
                CS%v_accel_bt, eta_pred, CS%uhbt, CS%vhbt, G, &
                CS%barotropic_CSp, CS%visc_rem_u, CS%visc_rem_v, &
                uhbt_out = uhbt_out, vhbt_out = vhbt_out, OBC = CS%OBC, &
                BT_cont = CS%BT_cont, eta_PF_start = eta_PF_start, &
                taux_bot=taux_bot, tauy_bot=tauy_bot)
    call cpu_clock_end(id_clock_btstep)
  else
    
    if (associated(CS%BT_cont) .or. CS%BT_use_layer_fluxes) then
      call cpu_clock_begin(id_clock_continuity)
      call continuity(u, v, h, hp, uh_in, vh_in, dt, G, &
                      CS%continuity_CSp, OBC=CS%OBC, visc_rem_u=CS%visc_rem_u, &
                      visc_rem_v=CS%visc_rem_v, BT_cont=CS%BT_cont)
      call cpu_clock_end(id_clock_continuity)
      if (BT_cont_BT_thick) then
        call cpu_clock_begin(id_clock_pass)
        call pass_vector(CS%BT_cont%h_u, CS%BT_cont%h_v, G%Domain, &
                         To_All+SCALAR_PAIR, CGRID_NE)
        call cpu_clock_end(id_clock_pass)
        call btcalc(h, G, CS%barotropic_CSp, CS%BT_cont%h_u, CS%BT_cont%h_v)
      endif
    endif
    
    if (CS%BT_use_layer_fluxes) then
      uh_ptr => uh_in; vh_ptr => vh_in; u_ptr => u; v_ptr => v
    endif

    u_init => u ; v_init => v
    call cpu_clock_begin(id_clock_btstep)
    if (calc_dtbt) call set_dtbt(G, CS%barotropic_CSp, eta, CS%pbce)
    call btstep(.false., u, v, eta, dt, u_bc_accel, v_bc_accel, &
                fluxes, CS%pbce, CS%eta_PF, u_av, v_av, CS%u_accel_bt, &
                CS%v_accel_bt, eta_pred, CS%uhbt, CS%vhbt, G, CS%barotropic_CSp,&
                CS%visc_rem_u, CS%visc_rem_v, OBC=CS%OBC, &
                BT_cont = CS%BT_cont, eta_PF_start=eta_PF_start, &
                taux_bot=taux_bot, tauy_bot=tauy_bot, &
                uh0=uh_ptr, vh0=vh_ptr, u_uh0=u_ptr, v_vh0=v_ptr)
    call cpu_clock_end(id_clock_btstep)
  endif

! up = u + dt_pred*( u_bc_accel + u_accel_bt )
  dt_pred = dt * CS%be
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    vp(i,J,k) = G%mask2dCv(i,J) * (v_init(i,J,k) + dt_pred * &
                    (v_bc_accel(i,J,k) + CS%v_accel_bt(i,J,k)))
  enddo ; enddo ;  enddo

  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    up(i,j,k) = G%mask2dCu(i,j) * (u_init(i,j,k) + dt_pred  * &
                    (u_bc_accel(I,j,k) + CS%u_accel_bt(I,j,k)))
  enddo ; enddo ;  enddo
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) then
    call uchksum(up,"Predictor 1 u",G,haloshift=0)
    call vchksum(vp,"Predictor 1 v",G,haloshift=0)
    call hchksum(G%H_to_kg_m2*h,"Predictor 1 h",G,haloshift=1)
    call uchksum(G%H_to_kg_m2*uh,"Predictor 1 uh",G,haloshift=2)
    call vchksum(G%H_to_kg_m2*vh,"Predictor 1 vh",G,haloshift=2)
!   call MOM_state_chksum("Predictor 1", up, vp, h, uh, vh, G, haloshift=1)
    call MOM_accel_chksum("Predictor accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv, &
             CS%diffu, CS%diffv, G, CS%pbce, CS%u_accel_bt, CS%v_accel_bt)
    call MOM_state_chksum("Predictor 1 init", u_init, v_init, h, uh, vh, G, haloshift=2)
    call check_redundant("Predictor 1 up", up, vp, G)
    call check_redundant("Predictor 1 uh", uh, vh, G)
  endif

! up <- up + dt_pred d/dz visc d/dz up
! u_av  <- u_av  + dt_pred d/dz visc d/dz u_av
  call cpu_clock_begin(id_clock_vertvisc)
  call vertvisc_coef(up, vp, h, fluxes, visc, dt_pred, G, CS%vertvisc_CSp)
  call vertvisc(up, vp, h, fluxes, visc, dt_pred, CS%OBC, G, &
                CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot)
  if (G%nonblocking_updates) then
    call cpu_clock_end(id_clock_vertvisc) ; call cpu_clock_begin(id_clock_pass)
    pid_u = pass_vector_start(up, vp, G%Domain)
    call cpu_clock_end(id_clock_pass) ; call cpu_clock_begin(id_clock_vertvisc)
  endif
  call vertvisc_remnant(visc, CS%visc_rem_u, CS%visc_rem_v, dt_pred, G, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)

  call cpu_clock_begin(id_clock_pass)
  call pass_vector(CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                   To_All+SCALAR_PAIR, CGRID_NE)
  if (G%nonblocking_updates) then
    call pass_vector_complete(pid_u, up, vp, G%Domain)
  else
    call pass_vector(up, vp, G%Domain)
  endif
  call cpu_clock_end(id_clock_pass)

! uh = u_av * h
! hp = h + dt * div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(up, vp, h, hp, uh, vh, dt, G, CS%continuity_CSp, &
                  CS%uhbt, CS%vhbt, CS%OBC, CS%visc_rem_u, CS%visc_rem_v, &
                  u_av, v_av, BT_cont=CS%BT_cont)
  call cpu_clock_end(id_clock_continuity)

  call cpu_clock_begin(id_clock_pass)
  call pass_var(hp, G%Domain)
  if (G%nonblocking_updates) then
    pid_u_av = pass_vector_start(u_av, v_av, G%Domain)
    pid_uh = pass_vector_start(uh(:,:,:), vh(:,:,:), G%Domain)
  else
    call pass_vector(u_av, v_av, G%Domain, complete=.false.)
    call pass_vector(uh(:,:,:), vh(:,:,:), G%Domain)
  endif
  call cpu_clock_end(id_clock_pass)

  if (associated(CS%OBC)) then
    call Radiation_Open_Bdry_Conds(CS%OBC, u_av, u_old_rad_OBC, v_av, &
             v_old_rad_OBC, hp, h_old_rad_OBC, G, CS%open_boundary_CSp)
  endif

! h_av = (h + hp)/2
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    h_av(i,j,k) = 0.5*(h(i,j,k) + hp(i,j,k))
  enddo ; enddo ; enddo

! The correction phase of the time step starts here.
  call enable_averaging(dt, Time_local, CS%diag)

!   Calculate a revised estimate of the free-surface height correction to be
! used in the next call to btstep.  This call is at this point so that
! hp can be changed if CS%begw /= 0.
! eta_cor = ...                 (hidden inside CS%barotropic_CSp)
  call cpu_clock_begin(id_clock_btcalc)
  call bt_mass_source(hp, eta_pred, fluxes, .false., dt_therm, &
                      dt_since_flux+dt, G, CS%barotropic_CSp)
  call cpu_clock_end(id_clock_btcalc)

  if (CS%begw /= 0.0) then
    ! hp <- (1-begw)*h_in + begw*hp
    ! Back up hp to the value it would have had after a time-step of
    ! begw*dt.  hp is not used again until recalculated by continuity.
    do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      hp(i,j,k) = (1.0-CS%begw)*h(i,j,k) + CS%begw*hp(i,j,k)
    enddo ; enddo ; enddo

! PFu = d/dx M(hp,T,S)
! pbce = dM/deta
    call cpu_clock_begin(id_clock_pres)
    call PressureForce(hp, tv, CS%PFu, CS%PFv, G, &
                       CS%PressureForce_CSp, CS%ALE_CSp, &
                       p_surf, CS%pbce, CS%eta_PF)
    call cpu_clock_end(id_clock_pres)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(CS%eta_PF, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  if (G%nonblocking_updates) then
    call cpu_clock_begin(id_clock_pass)
    call pass_vector_complete(pid_u_av, u_av, v_av, G%Domain)
    call pass_vector_complete(pid_uh, uh(:,:,:), vh(:,:,:), G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  if (BT_cont_BT_thick) then
    call cpu_clock_begin(id_clock_pass)
    call pass_vector(CS%BT_cont%h_u, CS%BT_cont%h_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE)
    call cpu_clock_end(id_clock_pass)
    call btcalc(h, G, CS%barotropic_CSp, CS%BT_cont%h_u, CS%BT_cont%h_v)
  endif

  if (CS%debug) then
    call MOM_state_chksum("Predictor ", up, vp, hp, uh, vh, G)
    call uchksum(u_av,"Predictor avg u",G,haloshift=1)
    call vchksum(v_av,"Predictor avg v",G,haloshift=1)
    call hchksum(G%H_to_kg_m2*h_av,"Predictor avg h",G,haloshift=0)
  ! call MOM_state_chksum("Predictor avg ", u_av, v_av,  h_av,uh, vh, G)
    call check_redundant("Predictor up ", up, vp, G)
    call check_redundant("Predictor uh ", uh, vh, G)
  endif

! diffu = horizontal viscosity terms (u_av)
  call cpu_clock_begin(id_clock_horvisc)
  call horizontal_viscosity(u_av, v_av, h_av, CS%diffu, CS%diffv, &
                            MEKE, Varmix, G, CS%hor_visc_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_horvisc)

! CAu = -(f+zeta_av)/h_av vh + d/dx KE_av
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(u_av, v_av, h_av, uh, vh, CS%CAu, CS%CAv, G,CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! Calculate the momentum forcing terms for the barotropic equations.

! u_bc_accel = CAu + PFu + diffu(u[n-1])
  call cpu_clock_begin(id_clock_btforce)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u_bc_accel(I,j,k) = (CS%Cau(I,j,k) + CS%PFu(I,j,k)) + CS%diffu(I,j,k)
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v_bc_accel(i,J,k) = (CS%Cav(i,J,k) + CS%PFv(i,J,k)) + CS%diffv(i,J,k)
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_btforce)

  if (CS%debug) then
    call check_redundant("corr pre-btstep CS%Ca ", CS%Cau, CS%Cav, G)
    call check_redundant("corr pre-btstep CS%PF ", CS%PFu, CS%PFv, G)
    call check_redundant("corr pre-btstep CS%diff ", CS%diffu, CS%diffv, G)
    call check_redundant("corr pre-btstep u_bc_accel ", u_bc_accel, v_bc_accel, G)
  endif

! u_accel_bt = layer accelerations due to barotropic solver
! pbce = dM/deta
  call cpu_clock_begin(id_clock_btstep)
  if (CS%flux_BT_coupling) then
    call btstep(.true., uh_in, vh_in, eta, dt, u_bc_accel, v_bc_accel, &
                fluxes, CS%pbce, CS%eta_PF, uh, vh, CS%u_accel_bt, &
                CS%v_accel_bt, eta, CS%uhbt, CS%vhbt, G, &
                CS%barotropic_CSp, CS%visc_rem_u, CS%visc_rem_v, etaav=eta_av, &
                uhbt_out = uhbt_out, vhbt_out = vhbt_out, OBC=CS%OBC, &
                BT_cont = CS%BT_cont, eta_PF_start = eta_PF_start, &
                taux_bot=taux_bot, tauy_bot=tauy_bot)
  else
    if (CS%BT_use_layer_fluxes) then
      uh_ptr => uh ; vh_ptr => vh ; u_ptr => u_av ; v_ptr => v_av
    endif

    call btstep(.false., u, v, eta, dt, u_bc_accel, v_bc_accel, &
                fluxes, CS%pbce, CS%eta_PF, u_av, v_av, CS%u_accel_bt, &
                CS%v_accel_bt, eta, CS%uhbt, CS%vhbt, G, &
                CS%barotropic_CSp, CS%visc_rem_u, CS%visc_rem_v, &
                etaav=eta_av, OBC=CS%OBC, &
                BT_cont = CS%BT_cont, eta_PF_start=eta_PF_start, &
                taux_bot=taux_bot, tauy_bot=tauy_bot, &
                uh0=uh_ptr, vh0=vh_ptr, u_uh0=u_ptr, v_vh0=v_ptr)
  endif
  call cpu_clock_end(id_clock_btstep)

  if (CS%debug) then
    call check_redundant("u_accel_bt ", CS%u_accel_bt, CS%v_accel_bt, G)
  endif

! u = u + dt*( u_bc_accel + u_accel_bt )
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u(i,j,k) = G%mask2dCu(i,j) * (u_init(i,j,k) + dt * &
                    (u_bc_accel(I,j,k) + CS%u_accel_bt(I,j,k)))
  enddo ; enddo ; enddo
  if (ASSOCIATED(CS%diag%PFu_tot)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      CS%diag%PFu_tot(i,j,k) = CS%PFu(i,j,k)
  enddo ; enddo ; enddo ; endif
  if (ASSOCIATED(CS%diag%CAu_tot)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      CS%diag%CAu_tot(i,j,k) = CS%CAu(i,j,k) !+ CS%u_accel_bt(i,j) - CS%diag%PFu_bt(i,j)
  enddo ; enddo ; enddo ; endif

  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v(i,j,k) = G%mask2dCv(i,j) * (v_init(i,j,k) + dt * &
                    (v_bc_accel(i,J,k) + CS%v_accel_bt(i,J,k)))
  enddo ; enddo ; enddo
  if (ASSOCIATED(CS%diag%PFv_tot)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    CS%diag%PFv_tot(i,j,k) = CS%PFv(i,j,k)
  enddo ; enddo ; enddo ; endif
  if (ASSOCIATED(CS%diag%CAv_tot)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    CS%diag%CAv_tot(i,j,k) = CS%CAv(i,j,k) !+ CS%v_accel_bt(i,j) - CS%diag%PFv_bt(i,j)
  enddo ; enddo ; enddo ; endif
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) then
    call uchksum(u,"Corrector 1 u",G,haloshift=0)
    call vchksum(v,"Corrector 1 v",G,haloshift=0)
    call hchksum(G%H_to_kg_m2*h,"Corrector 1 h",G,haloshift=2)
    call uchksum(G%H_to_kg_m2*uh,"Corrector 1 uh",G,haloshift=2)
    call vchksum(G%H_to_kg_m2*vh,"Corrector 1 vh",G,haloshift=2)
  ! call MOM_state_chksum("Corrector 1", u, v, h, uh, vh, G, haloshift=1)
    call MOM_accel_chksum("Corrector accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv, &
             CS%diffu, CS%diffv, G, CS%pbce, CS%u_accel_bt, CS%v_accel_bt)
  endif

! u <- u + dt d/dz visc d/dz u
! u_av <- u_av + dt d/dz visc d/dz u_av
  call cpu_clock_begin(id_clock_vertvisc)
  call vertvisc_coef(u, v, h, fluxes, visc, dt, G, CS%vertvisc_CSp)
  call vertvisc(u, v, h, fluxes, visc, dt, CS%OBC, G, &
                CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot)
  if (G%nonblocking_updates) then
    call cpu_clock_end(id_clock_vertvisc) ; call cpu_clock_begin(id_clock_pass)
    pid_u = pass_vector_start(u, v, G%Domain)
    call cpu_clock_end(id_clock_pass) ; call cpu_clock_begin(id_clock_vertvisc)
  endif
  call vertvisc_remnant(visc, CS%visc_rem_u, CS%visc_rem_v, dt, G, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)

! Later, h_av = (h_in + h_out)/2, but for now use h_av to store h_in.
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    h_av(i,j,k) = h(i,j,k)
  enddo ; enddo ; enddo

  call cpu_clock_begin(id_clock_pass)
  call pass_vector(CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                   To_All+SCALAR_PAIR, CGRID_NE)
  if (G%nonblocking_updates) then
    call pass_vector_complete(pid_u, u, v, G%Domain)
  else
    call pass_vector(u, v, G%Domain)
  endif
  call cpu_clock_end(id_clock_pass)

! uh = u_av * h
! h = h + dt * div . uh
  if (CS%flux_BT_coupling) then
    ! u_av and v_av adjusted so their mass transports match uhbt and vhbt.
    ! Also, determine the values of u and v so that their transports
    ! that agree with uhbt_out and vhbt_out.
    if (ASSOCIATED(CS%diag%du_adj)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      CS%diag%du_adj(I,j,k) = u(I,j,k)
    enddo ; enddo ; enddo ; endif
    if (ASSOCIATED(CS%diag%dv_adj)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      CS%diag%dv_adj(i,J,k) = v(i,J,k)
    enddo ; enddo ; enddo ; endif
    call cpu_clock_begin(id_clock_continuity)
    call continuity(u, v, h, h, uh, vh, dt, G, &
                    CS%continuity_CSp, CS%uhbt, CS%vhbt, CS%OBC, &
                    CS%visc_rem_u, CS%visc_rem_v, u_av, v_av, &
                    uhbt_out, vhbt_out, u, v)
    call cpu_clock_end(id_clock_continuity)
    if (G%nonblocking_updates) then
      call cpu_clock_begin(id_clock_pass)
      pid_h = pass_var_start(h, G%Domain)
      call cpu_clock_end(id_clock_pass)
    endif
    if (ASSOCIATED(CS%diag%du_adj)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      CS%diag%du_adj(I,j,k) = u(I,j,k) - CS%diag%du_adj(I,j,k)
    enddo ; enddo ; enddo ; endif
    if (ASSOCIATED(CS%diag%dv_adj)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      CS%diag%dv_adj(i,J,k) = v(i,J,k) - CS%diag%dv_adj(i,J,k)
    enddo ; enddo ; enddo ; endif

    call cpu_clock_begin(id_clock_vertvisc)
    call vertvisc_limit_vel(u, v, h_av, fluxes, visc, dt, G, CS%vertvisc_CSp)
    if (G%nonblocking_updates) then
      call cpu_clock_end(id_clock_vertvisc) ; call cpu_clock_begin(id_clock_pass)
      pid_u = pass_vector_start(u, v, G%Domain)
      call cpu_clock_end(id_clock_pass) ; call cpu_clock_begin(id_clock_vertvisc)
    endif
    call vertvisc_limit_vel(u_av, v_av, h_av, fluxes, visc, dt, G, CS%vertvisc_CSp)
    call cpu_clock_end(id_clock_vertvisc)

    call cpu_clock_begin(id_clock_pass)
    if (G%nonblocking_updates) then
      call pass_var_complete(pid_h, h, G%Domain)
      call pass_vector_complete(pid_u, u, v, G%Domain)
    else
      call pass_var(h, G%Domain)
      call pass_vector(u, v, G%Domain, complete=.false.)
    endif
    call cpu_clock_end(id_clock_pass)
  else
    ! u_av and v_av adjusted so their mass transports match uhbt and vhbt.
    call cpu_clock_begin(id_clock_continuity)
    call continuity(u, v, h, h, uh, vh, dt, G, &
                    CS%continuity_CSp, CS%uhbt, CS%vhbt, CS%OBC, &
                    CS%visc_rem_u, CS%visc_rem_v, u_av, v_av)
    call cpu_clock_end(id_clock_continuity)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(h, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  call cpu_clock_begin(id_clock_pass)
  if (G%nonblocking_updates) then
    pid_uh = pass_vector_start(uh(:,:,:), vh(:,:,:), G%Domain)
    pid_u_av = pass_vector_start(u_av, v_av, G%Domain)
  else
    call pass_vector(u_av, v_av, G%Domain, complete=.false.)
    call pass_vector(uh(:,:,:), vh(:,:,:), G%Domain)
  endif
  call cpu_clock_end(id_clock_pass)

  if (associated(CS%OBC)) then
    call Radiation_Open_Bdry_Conds(CS%OBC, u, u_old_rad_OBC, v, &
             v_old_rad_OBC, h, h_old_rad_OBC, G, CS%open_boundary_CSp)
  endif

! h_av = (h_in + h_out)/2 . Going in to this line, h_av = h_in.
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    h_av(i,j,k) = 0.5*(h_av(i,j,k) + h(i,j,k))
  enddo ; enddo ; enddo

  if (G%nonblocking_updates) then
    call cpu_clock_begin(id_clock_pass)
    call pass_vector_complete(pid_uh, uh(:,:,:), vh(:,:,:), G%Domain)
    call pass_vector_complete(pid_u_av, u_av, v_av, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  do k=1,nz ; do j=js-2,je+2 ; do I=Isq-2,Ieq+2
    uhtr(I,j,k) = uhtr(I,j,k) + uh(I,j,k)*dt
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq-2,Jeq+2 ; do i=is-2,ie+2
    vhtr(i,J,k) = vhtr(i,J,k) + vh(i,J,k)*dt
  enddo ; enddo ; enddo

!   The time-averaged free surface height has already been set by the last
!  call to btstep.

!   Here various terms used in to update the momentum equations are
! offered for averaging.
  if (CS%id_PFu > 0) call post_data(CS%id_PFu, CS%diag%PFu_tot, CS%diag)
  if (CS%id_PFv > 0) call post_data(CS%id_PFv, CS%diag%PFv_tot, CS%diag)
  if (CS%id_CAu > 0) call post_data(CS%id_CAu, CS%diag%CAu_tot, CS%diag)
  if (CS%id_CAv > 0) call post_data(CS%id_CAv, CS%diag%CAv_tot, CS%diag)

!   Here the thickness fluxes are offered for averaging.
  if (CS%id_uh > 0) call post_data(CS%id_uh, uh, CS%diag)
  if (CS%id_vh > 0) call post_data(CS%id_vh, vh, CS%diag)
  if (CS%id_uav > 0) call post_data(CS%id_uav, u_av, CS%diag)
  if (CS%id_vav > 0) call post_data(CS%id_vav, v_av, CS%diag)
  if (CS%id_u_BT_accel > 0) call post_data(CS%id_u_BT_accel, CS%u_accel_bt, CS%diag)
  if (CS%id_v_BT_accel > 0) call post_data(CS%id_v_BT_accel, CS%v_accel_bt, CS%diag)
  if (CS%id_du_adj > 0) call post_data(CS%id_du_adj, CS%diag%du_adj, CS%diag)
  if (CS%id_dv_adj > 0) call post_data(CS%id_dv_adj, CS%diag%dv_adj, CS%diag)
  if (CS%id_du_adj2 > 0) call post_data(CS%id_du_adj2, CS%diag%du_adj2, CS%diag)
  if (CS%id_dv_adj2 > 0) call post_data(CS%id_dv_adj2, CS%diag%dv_adj2, CS%diag)
  if (CS%debug) then
    call MOM_state_chksum("Corrector ", u, v, h, uh, vh, G)
    call uchksum(u_av,"Corrector avg u",G,haloshift=1)
    call vchksum(v_av,"Corrector avg v",G,haloshift=1)
    call hchksum(G%H_to_kg_m2*h_av,"Corrector avg h",G,haloshift=1)
 !  call MOM_state_chksum("Corrector avg ", u_av, v_av, h_av, uh, vh, G)
  endif

end subroutine step_MOM_dyn_split_RK2

! =============================================================================

subroutine adjustments_dyn_split_RK2(u, v, h, dt, G, CS)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in)    :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in)    :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: h
  real,                                   intent(in)    :: dt
  type(ocean_grid_type),                  intent(inout) :: G
  type(MOM_dyn_control_struct),           pointer       :: CS
 
! Arguments: u - The zonal velocity, in m s-1.
!  (in)      v - The meridional velocity, in m s-1.
!  (in)      h - The layer thicknesses, in m or kg m-2, depending on
!                whether the Boussinesq approximation is made.
!  (in)      dt - The time step in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure set up by initialize_MOM.

  ! Temporary arrays to contain layer thickness fluxes in m3 s-1 or kg s-1.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: uh_temp 
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: vh_temp
  ! A temporary array to contain layer projected thicknesses in m or kg m-2.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G))  :: h_temp
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (CS%readjust_BT_trans) then
    call cpu_clock_begin(id_clock_continuity)
    call continuity(u, v, h, h_temp, uh_temp, vh_temp, dt, G, &
                    CS%continuity_CSp, OBC=CS%OBC)
    call cpu_clock_end(id_clock_continuity)

    do j=js,je ; do I=is-1,ie ; CS%uhbt_in(I,j) = uh_temp(I,j,1) ; enddo ; enddo
    do k=2,nz ; do j=js,je ; do I=is-1,ie
      CS%uhbt_in(I,j) = CS%uhbt_in(I,j) + uh_temp(I,j,k)
    enddo ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; CS%vhbt_in(i,J) = vh_temp(i,J,1) ; enddo ; enddo
    do k=2,nz ; do J=js-1,je ; do i=is,ie
      CS%vhbt_in(i,J) = CS%vhbt_in(i,J) + vh_temp(i,J,k)
    enddo ; enddo ; enddo
    CS%readjust_velocity = .true.
  endif

end subroutine adjustments_dyn_split_RK2

! =============================================================================

subroutine register_restarts_dyn_split_RK2(G, param_file, CS, restart_CS, uh, vh)
  type(ocean_grid_type),         intent(in)    :: G
  type(param_file_type),         intent(in)    :: param_file
  type(MOM_dyn_control_struct),  intent(inout) :: CS
  type(MOM_restart_CS),          pointer       :: restart_CS
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), target, intent(inout) :: uh
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), target, intent(inout) :: vh
!   This subroutine sets up any auxiliary restart variables that are specific
! to the unsplit time stepping scheme.  All variables registered here should
! have the ability to be recreated if they are not present in a restart file.

! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      CS - The control structure set up by initialize_MOM.
!  (in)      restart_CS - A pointer to the restart control structure.
!  (inout)   uh - The zonal volume or mass transport, in m3 s-1 or kg s-1.
!  (inout)   vh - The meridional volume or mass transport, in m3 s-1 or kg s-1.

  type(vardesc) :: vd
  character(len=48) :: thickness_units, flux_units
  logical :: adiabatic
  integer :: isd, ied, jsd, jed, nz, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

! This is where a control structure that is specific to this module would be allocated.
! if (associated(CS_split)) then
!   call MOM_error(WARNING, "register_restarts_split_RK2 called with an associated "// &
!                            "control structure.")
!   return
! endif
! allocate(CS_split)

  call get_param(param_file, "MOM", "BE", CS%be, &
                 "If SPLIT is true, BE determines the relative weighting \n"//&
                 "of a  2nd-order Runga-Kutta baroclinic time stepping \n"//&
                 "scheme (0.5) and a backward Euler scheme (1) that is \n"//&
                 "used for the Coriolis and inertial terms.  BE may be \n"//&
                 "from 0.5 to 1, but instability may occur near 0.5. \n"//&
                 "BE is also applicable if SPLIT is false and USE_RK2 \n"//&
                 "is true.", units="nondim", default=0.6)
  call get_param(param_file, "MOM", "BEGW", CS%begw, &
                 "If SPILT is true, BEGW is a number from 0 to 1 that \n"//&
                 "controls the extent to which the treatment of gravity \n"//&
                 "waves is forward-backward (0) or simulated backward \n"//&
                 "Euler (1).  0 is almost always used.\n"//&
                 "If SPLIT is false and USE_RK2 is true, BEGW can be \n"//&
                 "between 0 and 0.5 to damp gravity waves.", &
                 units="nondim", default=0.0)

  call get_param(param_file, "MOM", "FLUX_BT_COUPLING", CS%flux_BT_coupling, &
                 "If true, use mass fluxes to ensure consistency between \n"//&
                 "the baroclinic and barotropic modes. This is only used \n"//&
                 "if SPLIT is true.", default=.false.)
  call get_param(param_file, "MOM", "READJUST_BT_TRANS", CS%readjust_BT_trans, &
                 "If true, make a barotropic adjustment to the layer \n"//&
                 "velocities after the thermodynamic part of the step \n"//&
                 "to ensure that the interaction between the thermodynamics \n"//&
                 "and the continuity solver do not change the barotropic \n"//&
                 "transport.  This is only used if FLUX_BT_COUPLING and \n"//&
                 "SPLIT are true.", default=.false.)
  call get_param(param_file, "MOM", "SPLIT_BOTTOM_STRESS", CS%split_bottom_stress, &
                 "If true, provide the bottom stress calculated by the \n"//&
                 "vertical viscosity to the barotropic solver.", default=.false.)
  call get_param(param_file, "MOM", "BT_USE_LAYER_FLUXES", CS%BT_use_layer_fluxes, &
                 "If true, use the summed layered fluxes plus an \n"//&
                 "adjustment due to the change in the barotropic velocity \n"//&
                 "in the barotropic continuity equation.", default=.true.)
  adiabatic=.false. ; call read_param(param_file, "ADIABATIC", adiabatic)
  if (.not.CS%flux_BT_coupling .or. adiabatic) CS%readjust_BT_trans = .false.

  ALLOC_(CS%diffu(IsdB:IedB,jsd:jed,nz)) ; CS%diffu(:,:,:) = 0.0
  ALLOC_(CS%diffv(isd:ied,JsdB:JedB,nz)) ; CS%diffv(:,:,:) = 0.0
  ALLOC_(CS%CAu(IsdB:IedB,jsd:jed,nz)) ; CS%CAu(:,:,:) = 0.0
  ALLOC_(CS%CAv(isd:ied,JsdB:JedB,nz)) ; CS%CAv(:,:,:) = 0.0
  ALLOC_(CS%PFu(IsdB:IedB,jsd:jed,nz)) ; CS%PFu(:,:,:) = 0.0
  ALLOC_(CS%PFv(isd:ied,JsdB:JedB,nz)) ; CS%PFv(:,:,:) = 0.0

  ALLOC_(CS%eta(isd:ied,jsd:jed))       ; CS%eta(:,:) = 0.0
  ALLOC_(CS%u_av(IsdB:IedB,jsd:jed,nz)) ; CS%u_av(:,:,:) = 0.0
  ALLOC_(CS%v_av(isd:ied,JsdB:JedB,nz)) ; CS%v_av(:,:,:) = 0.0
  ALLOC_(CS%h_av(isd:ied,jsd:jed,nz))   ; CS%h_av(:,:,:) = G%Angstrom
  ALLOC_(CS%uhbt_in(IsdB:IedB,jsd:jed)) ; CS%uhbt_in(:,:) = 0.0
  ALLOC_(CS%vhbt_in(isd:ied,JsdB:JedB)) ; CS%vhbt_in(:,:) = 0.0

  thickness_units = get_thickness_units(G)
  flux_units = get_flux_units(G)

 ! if (G%Boussinesq) then
    vd = vardesc("sfc","Free surface Height",'h','1','s',thickness_units)
 ! else
 !   vd(1) = vardesc("ocean_mass?","Ocean column mass",'h','1','s',"kg meter-2")
 ! endif
  call register_restart_field(CS%eta, vd, .false., restart_CS)

  vd = vardesc("u2","Auxiliary Zonal velocity",'u','L','s',"meter second-1")
  call register_restart_field(CS%u_av, vd, .false., restart_CS)

  vd = vardesc("v2","Auxiliary Meridional velocity",'v','L','s',"meter second-1")
  call register_restart_field(CS%v_av, vd, .false., restart_CS)

  vd = vardesc("h2","Auxiliary Layer Thickness",'h','L','s',thickness_units)
  call register_restart_field(CS%h_av, vd, .false., restart_CS)

  vd = vardesc("uh","Zonal thickness flux",'u','L','s',flux_units)
  call register_restart_field(uh, vd, .false., restart_CS)

  vd = vardesc("vh","Meridional thickness flux",'v','L','s',flux_units)
  call register_restart_field(vh, vd, .false., restart_CS)

  vd = vardesc("diffu","Zonal horizontal viscous acceleration",'u','L','s', &
               "meter second-2")
  call register_restart_field(CS%diffu, vd, .false., restart_CS)

  vd = vardesc("diffv","Meridional horizontal viscous acceleration",'v','L','s',&
               "meter second-2")
  call register_restart_field(CS%diffv, vd, .false., restart_CS)

  call register_barotropic_restarts(G, param_file, CS%barotropic_CSp, &
                                    restart_CS)

  if (CS%readjust_bt_trans) then
    vd = vardesc("uhbt_in","Final instantaneous barotropic zonal thickness flux",&
                 'u','1','s',flux_units)
    call register_restart_field(CS%uhbt_in, vd, .false., restart_CS)

    vd = vardesc("vhbt_in","Final instantaneous barotropic meridional thickness flux",&
                 'v','1','s',flux_units)
    call register_restart_field(CS%vhbt_in, vd, .false., restart_CS)
  endif

end subroutine register_restarts_dyn_split_RK2

subroutine initialize_dyn_split_RK2(u, v, h, uh, vh, Time, G, param_file, &
                                    diag, CS, restart_CS, dt, MIS, VarMix, MEKE)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(inout) :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(inout) :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_) , intent(inout) :: h
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), target, intent(inout) :: uh
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), target, intent(inout) :: vh
  type(time_type),                target, intent(in)    :: Time
  type(ocean_grid_type),                  intent(inout) :: G
  type(param_file_type),                  intent(in)    :: param_file
  type(diag_ptrs),                target, intent(inout) :: diag
  type(MOM_dyn_control_struct),   target, intent(inout) :: CS
  type(MOM_restart_CS),                   pointer       :: restart_CS
  real,                                   intent(in)    :: dt
  type(ocean_internal_state),             intent(inout) :: MIS
  type(VarMix_CS),                        pointer       :: VarMix
  type(MEKE_type),                        pointer       :: MEKE
  

  !   This subroutine initializes any variables that are specific to this time
  ! stepping scheme, including the cpu clocks.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_tmp
  character(len=48) :: thickness_units, flux_units
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  allocate(CS%taux_bot(IsdB:IedB,jsd:jed)) ; CS%taux_bot(:,:) = 0.0
  allocate(CS%tauy_bot(isd:ied,JsdB:JedB)) ; CS%tauy_bot(:,:) = 0.0

  ALLOC_(CS%uhbt(IsdB:IedB,jsd:jed))   ; CS%uhbt(:,:) = 0.0
  ALLOC_(CS%vhbt(isd:ied,JsdB:JedB))   ; CS%vhbt(:,:) = 0.0
  ALLOC_(CS%visc_rem_u(IsdB:IedB,jsd:jed,nz)) ; CS%visc_rem_u(:,:,:) = 0.0
  ALLOC_(CS%visc_rem_v(isd:ied,JsdB:JedB,nz)) ; CS%visc_rem_v(:,:,:) = 0.0
  ALLOC_(CS%eta_PF(isd:ied,jsd:jed))   ; CS%eta_PF(:,:) = 0.0
  ALLOC_(CS%pbce(isd:ied,jsd:jed,nz))  ; CS%pbce(:,:,:) = 0.0

  ALLOC_(CS%u_accel_bt(IsdB:IedB,jsd:jed,nz)) ; CS%u_accel_bt(:,:,:) = 0.0
  ALLOC_(CS%v_accel_bt(isd:ied,JsdB:JedB,nz)) ; CS%v_accel_bt(:,:,:) = 0.0

  MIS%diffu => CS%diffu ; MIS%diffv => CS%diffv
  MIS%PFu => CS%PFu ; MIS%PFv => CS%PFv
  MIS%CAu => CS%CAu ; MIS%CAv => CS%CAv
  MIS%pbce => CS%pbce
  MIS%u_accel_bt => CS%u_accel_bt ; MIS%v_accel_bt => CS%v_accel_bt
  MIS%u_av => CS%u_av ; MIS%v_av => CS%v_av

  call continuity_init(Time, G, param_file, diag, CS%continuity_CSp)
  call CoriolisAdv_init(Time, G, param_file, diag, CS%CoriolisAdv_CSp)
  call PressureForce_init(Time, G, param_file, diag, CS%PressureForce_CSp, &
                          CS%tides_CSp)
  call hor_visc_init(Time, G, param_file, diag, CS%hor_visc_CSp)

  if (.not. query_initialized(CS%eta,"sfc",restart_CS))  then
    ! Estimate eta based on the layer thicknesses - h.  With the Boussinesq
    ! approximation, eta is the free surface height anomaly, while without it
    ! eta is the mass of ocean per unit area.  eta always has the same
    ! dimensions as h, either m or kg m-3.  
    !   CS%eta(:,:) = 0.0 already from initialization.
    if (G%Boussinesq) then
      do j=js,je ; do i=is,ie ; CS%eta(i,j) = -G%bathyT(i,j) ; enddo ; enddo
    endif
    do k=1,nz ; do j=js,je ; do i=is,ie
       CS%eta(i,j) = CS%eta(i,j) + h(i,j,k)
    enddo ; enddo ; enddo
  endif  

  call barotropic_init(u, v, h, CS%eta, Time, G, &
                       param_file, diag, CS%barotropic_CSp, restart_CS, &
                       CS%BT_cont, CS%tides_CSp)

  if (.not. query_initialized(CS%diffu,"diffu",restart_CS) .or. &
      .not. query_initialized(CS%diffv,"diffv",restart_CS)) &
    call horizontal_viscosity(u, v, h, CS%diffu, CS%diffv, MEKE, VarMix, &
                              G, CS%hor_visc_CSp)
  if (.not. query_initialized(CS%u_av,"u2", restart_CS) .or. &
      .not. query_initialized(CS%u_av,"v2", restart_CS)) then
    CS%u_av(:,:,:) = u(:,:,:)
    CS%v_av(:,:,:) = v(:,:,:)
  endif
! This call is just here to initialize uh and vh.
  if (.not. query_initialized(uh,"uh",restart_CS) .or. &
      .not. query_initialized(vh,"vh",restart_CS)) then
    h_tmp(:,:,:) = h(:,:,:)
    call continuity(u, v, h, h_tmp, uh, vh, dt, G, CS%continuity_CSp, OBC=CS%OBC)
    call cpu_clock_begin(id_clock_pass_init)
    call pass_var(h_tmp, G%Domain)
    call cpu_clock_end(id_clock_pass_init)
    CS%h_av(:,:,:) = 0.5*(h(:,:,:) + h_tmp(:,:,:))
  else
    if (.not. query_initialized(CS%h_av,"h2",restart_CS)) &
      CS%h_av(:,:,:) = h(:,:,:)
  endif

  !   Determine whether there is a barotropic transport that is to be used
  ! to adjust the layers' velocities.
  CS%readjust_velocity = .false.
  if (CS%readjust_BT_trans) then
    if (query_initialized(CS%uhbt_in,"uhbt_in",restart_CS) .and. &
        query_initialized(CS%vhbt_in,"vhbt_in",restart_CS)) then
      CS%readjust_velocity = .true.
      call cpu_clock_begin(id_clock_pass_init)
      call pass_vector(CS%uhbt_in, CS%vhbt_in, G%Domain)
      call cpu_clock_end(id_clock_pass_init)
    endif
  endif

  call cpu_clock_begin(id_clock_pass_init)
  call pass_vector(CS%u_av,CS%v_av, G%Domain)
  call pass_var(CS%h_av, G%Domain)
  call pass_vector(uh, vh, G%Domain)
  call cpu_clock_end(id_clock_pass_init)

  flux_units = get_flux_units(G)
  CS%id_uh = register_diag_field('ocean_model', 'uh', G%axesCuL, Time, &
      'Zonal Thickness Flux', flux_units)
  CS%id_vh = register_diag_field('ocean_model', 'vh', G%axesCvL, Time, &
      'Meridional Thickness Flux', flux_units)
  CS%id_CAu = register_diag_field('ocean_model', 'CAu', G%axesCuL, Time, &
      'Zonal Coriolis and Advective Acceleration', 'meter second-2')
  CS%id_CAv = register_diag_field('ocean_model', 'CAv', G%axesCvL, Time, &
      'Meridional Coriolis and Advective Acceleration', 'meter second-2')
  CS%id_PFu = register_diag_field('ocean_model', 'PFu', G%axesCuL, Time, &
      'Zonal Pressure Force Acceleration', 'meter second-2')
  CS%id_PFv = register_diag_field('ocean_model', 'PFv', G%axesCvL, Time, &
      'Meridional Pressure Force Acceleration', 'meter second-2')
  if (CS%id_PFu > 0) call safe_alloc_ptr(diag%PFu_tot,IsdB,IedB,jsd,jed,nz)
  if (CS%id_PFv > 0) call safe_alloc_ptr(diag%PFv_tot,isd,ied,JsdB,JedB,nz)
  if (CS%id_CAu > 0) call safe_alloc_ptr(diag%CAu_tot,IsdB,IedB,jsd,jed,nz)
  if (CS%id_CAv > 0) call safe_alloc_ptr(diag%CAv_tot,isd,ied,JsdB,JedB,nz)

  CS%id_uav = register_diag_field('ocean_model', 'uav', G%axesCuL, Time, &
      'Barotropic-step Averaged Zonal Velocity', 'meter second-1')
  CS%id_vav = register_diag_field('ocean_model', 'vav', G%axesCvL, Time, &
      'Barotropic-step Averaged Meridional Velocity', 'meter second-1')


  if (CS%flux_BT_coupling) then
    CS%id_du_adj = register_diag_field('ocean_model', 'du_adj', G%axesCuL, Time, &
        'Zonal velocity Adjustment 1', 'meter second-1')
    CS%id_dv_adj = register_diag_field('ocean_model', 'dv_adj', G%axesCvL, Time, &
        'Meridional velocity Adjustment 1', 'meter second-1')
    if (CS%id_du_adj > 0) call safe_alloc_ptr(CS%diag%du_adj,IsdB,IedB,jsd,jed,nz)
    if (CS%id_dv_adj > 0) call safe_alloc_ptr(CS%diag%dv_adj,isd,ied,JsdB,JedB,nz)
    if (CS%readjust_BT_trans) then
      CS%id_du_adj2 = register_diag_field('ocean_model', 'du_adj2', G%axesCuL, Time, &
          'Zonal velocity Adjustment 2', 'meter second-1')
      CS%id_dv_adj2 = register_diag_field('ocean_model', 'dv_adj2', G%axesCvL, Time, &
          'Meridional velocity Adjustment 2', 'meter second-1')
      if (CS%id_du_adj2 > 0) call safe_alloc_ptr(CS%diag%du_adj2,IsdB,IedB,jsd,jed,nz)
      if (CS%id_dv_adj2 > 0) call safe_alloc_ptr(CS%diag%dv_adj2,isd,ied,JsdB,JedB,nz)
    endif
  endif
  CS%id_u_BT_accel = register_diag_field('ocean_model', 'u_BT_accel', G%axesCuL, Time, &
    'Barotropic Anomaly Zonal Acceleration', 'meter second-1')
  CS%id_v_BT_accel = register_diag_field('ocean_model', 'v_BT_accel', G%axesCvL, Time, &
    'Barotropic Anomaly Meridional Acceleration', 'meter second-1')

  if (CS%debug_truncations) then
    if (CS%flux_BT_coupling) then
      call safe_alloc_ptr(diag%du_adj,IsdB,IedB,jsd,jed,nz)
      call safe_alloc_ptr(diag%dv_adj,isd,ied,JsdB,JedB,nz)
    endif
    if (CS%flux_BT_coupling .and. CS%readjust_BT_trans) then
      call safe_alloc_ptr(diag%du_adj2,IsdB,IedB,jsd,jed,nz)
      call safe_alloc_ptr(diag%dv_adj2,isd,ied,JsdB,JedB,nz)
    endif
  endif

  id_clock_Cor = cpu_clock_id('(Ocean Coriolis & mom advection)', grain=CLOCK_MODULE)
  id_clock_continuity = cpu_clock_id('(Ocean continuity equation)', grain=CLOCK_MODULE)
  id_clock_pres = cpu_clock_id('(Ocean pressure force)', grain=CLOCK_MODULE)
  id_clock_vertvisc = cpu_clock_id('(Ocean vertical viscosity)', grain=CLOCK_MODULE)
  id_clock_horvisc = cpu_clock_id('(Ocean horizontal viscosity)', grain=CLOCK_MODULE)
  id_clock_mom_update = cpu_clock_id('(Ocean momentum increments)', grain=CLOCK_MODULE)
  id_clock_pass = cpu_clock_id('(Ocean message passing)', grain=CLOCK_MODULE)
  id_clock_pass_init = cpu_clock_id('(Ocean init message passing)', grain=CLOCK_ROUTINE)
  id_clock_btcalc = cpu_clock_id('(Ocean barotropic mode calc)', grain=CLOCK_MODULE)
  id_clock_btstep = cpu_clock_id('(Ocean barotropic mode stepping)', grain=CLOCK_MODULE)
  id_clock_btforce = cpu_clock_id('(Ocean barotropic forcing calc)', grain=CLOCK_MODULE)


end subroutine initialize_dyn_split_RK2

subroutine end_dyn_split_RK2(CS)
  type(MOM_dyn_control_struct), pointer :: CS

  DEALLOC_(CS%diffu) ; DEALLOC_(CS%diffv)
  DEALLOC_(CS%CAu)   ; DEALLOC_(CS%CAv)
  DEALLOC_(CS%PFu)   ; DEALLOC_(CS%PFv)
  
  DEALLOC_(CS%uhbt) ; DEALLOC_(CS%vhbt)
  DEALLOC_(CS%u_accel_bt) ; DEALLOC_(CS%v_accel_bt)
  DEALLOC_(CS%visc_rem_u) ; DEALLOC_(CS%visc_rem_v)

  DEALLOC_(CS%eta) ; DEALLOC_(CS%eta_PF) ; DEALLOC_(CS%pbce)
  DEALLOC_(CS%h_av) ; DEALLOC_(CS%u_av) ; DEALLOC_(CS%v_av)
  DEALLOC_(CS%uhbt_in) ; DEALLOC_(CS%vhbt_in)

  call dealloc_BT_cont_type(CS%BT_cont)

  deallocate(CS)
end subroutine end_dyn_split_RK2

end module MOM_dynamics_split_RK2
