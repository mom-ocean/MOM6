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

use MOM_barotropic, only : barotropic_init, btstep, btcalc, bt_mass_source
use MOM_barotropic, only : register_barotropic_restarts, set_dtbt, barotropic_CS
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

public step_MOM_dyn_split_RK2, register_restarts_dyn_split_RK2
public initialize_dyn_split_RK2

integer :: id_clock_Cor, id_clock_pres, id_clock_vertvisc
integer :: id_clock_horvisc, id_clock_mom_update
integer :: id_clock_continuity, id_clock_thick_diff, id_clock_ml_restrat
integer :: id_clock_btstep, id_clock_btcalc, id_clock_btforce
integer :: id_clock_pass, id_clock_pass_init

contains

! =============================================================================

subroutine step_MOM_dyn_split_RK2(u_in, v_in, h_in, eta_in, uhbt_in, vhbt_in, &
                 Time_local, dt, fluxes, p_surf_begin, p_surf_end, &
                 dt_since_flux, dt_therm, uh, vh, u_av, v_av, h_av, eta_av, &
                 u_out, v_out, h_out, eta_out, G, CS)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), target, intent(in) :: u_in
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), target, intent(in) :: v_in
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: h_in
  real, dimension(NIMEM_,NJMEM_),         intent(inout) :: eta_in
  real, dimension(NIMEMB_,NJMEM_),        intent(inout) :: uhbt_in
  real, dimension(NIMEM_,NJMEMB_),        intent(inout) :: vhbt_in
  type(time_type),                        intent(in)    :: Time_local
  real,                                   intent(in)    :: dt
  type(forcing),                          intent(in)    :: fluxes
  real, dimension(:,:),                   pointer       :: p_surf_begin, p_surf_end
  real,                                   intent(in)    :: dt_since_flux, dt_therm
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), target, intent(inout) :: uh
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), target, intent(inout) :: vh
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), target, intent(inout) :: u_av
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), target, intent(inout) :: v_av
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(inout) :: h_av
  real, dimension(NIMEM_,NJMEM_),         intent(out)   :: eta_av
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(out)   :: u_out
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(out)   :: v_out
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(out)   :: h_out
  real, dimension(NIMEM_,NJMEM_),         intent(out)   :: eta_out
  type(ocean_grid_type),                  intent(inout) :: G
  type(MOM_control_struct),               pointer       :: CS
! Arguments: u_in - The input zonal velocity, in m s-1.
!  (in)      v_in - The input meridional velocity, in m s-1.
!  (in)      h_in - The input layer thicknesses, in m or kg m-2, depending on
!                   whether the Boussinesq approximation is made.
!  (in)      eta_in - The input free surface height or column mass, in m or
!                     kg m-2.  (Intent inout for halo updates.)
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
!  (inout)   u_av - The zonal velocity time-averaged over a time step, in m s-1.
!  (inout)   v_av - The meridional velocity time-averaged over a time step, in m s-1.
!  (inout)   h_av - The layer thickness time-averaged over a time step, in m or
!                   kg m-2.
!  (out)     eta_av - The free surface height or column mass time-averaged
!                     over a time step, in m or kg m-2.
!  (out)     u_out - The output zonal velocity, in m s-1.
!  (out)     v_out - The output meridional velocity, in m s-1.
!  (out)     h_out - The output layer thicknesses, in m or kg m-2.
!  (out)     eta_out - The output free surface height or column mass, in m or
!                      kg m-2.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure set up by initialize_MOM.

  real :: dt_pred   ! The time step for the predictor part of the baroclinic
                    ! time stepping.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_tmp
    ! A temporary estimated thickness, in m.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: uh_tmp
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: vh_tmp
    ! Temporary transports, in m3 s-1 or kg s-1.
  real, dimension(SZIB_(G),SZJ_(G)) :: u_dhdt
  real, dimension(SZI_(G),SZJB_(G)) :: v_dhdt
    ! Vertically summed transport tendencies in m3 s-2 or kg s-2.

  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: u_bc_accel
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: v_bc_accel
    ! u_bc_accel and v_bc_accel are the summed baroclinic accelerations of each
    ! layer calculated by the non-barotropic part of the model, both in m s-2.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), target :: uh_in
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), target :: vh_in
    ! uh_in and vh_in are the zonal or meridional mass transports that would be
    ! obtained using the velocities u_in and v_in, both in m3 s-1 or kg s-1.
  real, dimension(SZIB_(G),SZJ_(G)) :: uhbt_out
  real, dimension(SZI_(G),SZJB_(G)) :: vhbt_out
    ! uhbt_out and vhbt_out are the vertically summed transports from the
    ! barotropic solver based on its final velocities, both in m3 s-1 or kg s-1.
  real, dimension(SZI_(G),SZJ_(G)) :: eta_pred
    ! eta_pred is the predictor value of the free surface height or column mass,
    ! in m or kg m-2.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), target :: u_adj
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), target :: v_adj
    ! u_adj and v_adj are the zonal or meridional velocities after u_in and v_in
    ! have been barotropically adjusted so the resulting transports match
    ! uhbt_out and vhbt_out, both in m s-1.
  real :: Pa_to_eta ! A factor that converts pressures to the units of eta.
  real, pointer, dimension(:,:,:) :: u_init, v_init  ! Pointers to u_in and v_in
                                                     ! or u_adj and v_adj.
  real, pointer, dimension(:,:)   :: p_surf => NULL(), eta_PF_start => NULL()
  real, pointer, dimension(:,:)   :: taux_bot => NULL(), tauy_bot => NULL()
  real, pointer, dimension(:,:,:) :: uh_ptr => NULL(), u_ptr => NULL()
  real, pointer, dimension(:,:,:) :: vh_ptr => NULL(), v_ptr => NULL()
  real :: Idt
  logical :: dyn_p_surf
  logical :: BT_cont_BT_thick ! If true, use the BT_cont_type to estimate the
                              ! relative weightings of the layers in calculating
                              ! the barotropic accelerations.
  integer :: pid_Ray, pid_bbl_h, pid_kv_bbl, pid_eta_PF, pid_eta_in, pid_visc
  integer :: pid_h, pid_u, pid_u_av, pid_uh, pid_uhbt_in
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  Idt = 1.0 / dt

  if (CS%debug) then
    call MOM_state_chksum("Start predictor ", u_in, v_in, h_in, uh, vh, G)
    call check_redundant("Start predictor u_in ", u_in, v_in, G)
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

  BT_cont_BT_thick = .false.
  if (associated(CS%BT_cont)) BT_cont_BT_thick = &
    (associated(CS%BT_cont%h_u) .and. associated(CS%BT_cont%h_v))

  if (CS%split_bottom_stress) then
    taux_bot => CS%taux_bot ; tauy_bot => CS%tauy_bot
  endif

  if (CS%calc_bbl) then
    ! Calculate the BBL properties and store them inside CS%visc (u_in,h_in).
    call cpu_clock_begin(id_clock_vertvisc)
    call enable_averaging(CS%bbl_calc_time_interval, &
                          Time_local-set_time(int(dt)), CS%diag)
    call set_viscous_BBL(u_in, v_in, h_in, CS%tv, CS%visc, G, CS%set_visc_CSp)
    call disable_averaging(CS%diag)
    call cpu_clock_end(id_clock_vertvisc)

    call cpu_clock_begin(id_clock_pass)
    if (G%nonblocking_updates) then   
      if (associated(CS%visc%Ray_u) .and. associated(CS%visc%Ray_v)) &
        pid_Ray = pass_vector_start(CS%visc%Ray_u, CS%visc%Ray_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(CS%visc%bbl_thick_u) .and. associated(CS%visc%bbl_thick_v)) &
        pid_bbl_h = pass_vector_start(CS%visc%bbl_thick_u, CS%visc%bbl_thick_v, &
                       G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(CS%visc%kv_bbl_u) .and. associated(CS%visc%kv_bbl_v)) &
        pid_kv_bbl = pass_vector_start(CS%visc%kv_bbl_u, CS%visc%kv_bbl_v, &
                       G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
      ! CS%calc_bbl will be set to .false. when the message passing is complete.
    else
      if (associated(CS%visc%Ray_u) .and. associated(CS%visc%Ray_v)) &
        call pass_vector(CS%visc%Ray_u, CS%visc%Ray_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(CS%visc%kv_bbl_u) .and. associated(CS%visc%kv_bbl_v)) then
        call pass_vector(CS%visc%bbl_thick_u, CS%visc%bbl_thick_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE, complete=.false.)
        call pass_vector(CS%visc%kv_bbl_u, CS%visc%kv_bbl_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE)
      endif
      CS%calc_bbl = .false.
    endif
    call cpu_clock_end(id_clock_pass)
  endif

! PFu = d/dx M(h_in,T,S)
! pbce = dM/deta
  if (CS%begw == 0.0) call enable_averaging(dt, Time_local, CS%diag)
  call cpu_clock_begin(id_clock_pres)
  call PressureForce(h_in, CS%tv, CS%PFu, CS%PFv, G, CS%PressureForce_CSp, &
                     CS%regridding_opts, p_surf, CS%pbce, CS%eta_PF)
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
    pid_eta_PF = pass_var_start(CS%eta_PF(:,:), G%Domain)
    pid_eta_in = pass_var_start(eta_in(:,:), G%Domain)
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
    if (CS%calc_bbl) then
      if (associated(CS%visc%Ray_u) .and. associated(CS%visc%Ray_v)) &
        call pass_vector_complete(pid_Ray, CS%visc%Ray_u, CS%visc%Ray_v, G%Domain, &
                         To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(CS%visc%bbl_thick_u) .and. associated(CS%visc%bbl_thick_v)) &
        call pass_vector_complete(pid_bbl_h, CS%visc%bbl_thick_u, CS%visc%bbl_thick_v, &
                       G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
      if (associated(CS%visc%kv_bbl_u) .and. associated(CS%visc%kv_bbl_v)) &
        call pass_vector_complete(pid_kv_bbl, CS%visc%kv_bbl_u, CS%visc%kv_bbl_v, &
                       G%Domain, To_All+SCALAR_PAIR, CGRID_NE)

      ! CS%calc_bbl is set to .false. now that the message passing is completed.
      CS%calc_bbl = .false.
    endif
    call pass_var_complete(pid_eta_PF, CS%eta_PF(:,:), G%Domain)
    call pass_var_complete(pid_eta_in, eta_in(:,:), G%Domain)
    if (CS%readjust_velocity) &
      call pass_vector_complete(pid_uhbt_in, uhbt_in, vhbt_in, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  call cpu_clock_begin(id_clock_vertvisc)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u_out(i,j,k) = G%umask(i,j) * (u_in(i,j,k) + dt * u_bc_accel(I,j,k))
  enddo ; enddo ;  enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v_out(i,j,k) = G%vmask(i,j) * (v_in(i,j,k) + dt * v_bc_accel(i,J,k))
  enddo ; enddo ;  enddo
  call enable_averaging(dt, Time_local, CS%diag)
  call set_viscous_ML(u_in, v_in, h_in, CS%tv, fluxes, CS%visc, dt, G, &
                      CS%set_visc_CSp)
  call disable_averaging(CS%diag)

  call vertvisc_coef(u_out, v_out, h_in, fluxes, CS%visc, dt, G, CS%vertvisc_CSp)
  call vertvisc_remnant(CS%visc, CS%visc_rem_u, CS%visc_rem_v, dt, G, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)

  call cpu_clock_begin(id_clock_pass)
  if (G%nonblocking_updates) then
    pid_visc = pass_vector_start(CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                                 To_All+SCALAR_PAIR, CGRID_NE)
  else
    call pass_var(CS%eta_PF(:,:), G%Domain, complete=.false.)
    call pass_var(eta_in(:,:), G%Domain)
    if (CS%readjust_velocity) call pass_vector(uhbt_in, vhbt_in, G%Domain)
    call pass_vector(CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                     To_All+SCALAR_PAIR, CGRID_NE)
  endif
  call cpu_clock_end(id_clock_pass)

  call cpu_clock_begin(id_clock_btcalc)
  ! Calculate the relative layer weights for determining barotropic quantities.
  if (.not.BT_cont_BT_thick) &
    call btcalc(h_in, G, CS%barotropic_CSp)
  call bt_mass_source(h_in, eta_in, fluxes, .true., dt_therm, dt_since_flux, &
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
      call continuity(u_in, v_in, h_in, h_out, uh_in, vh_in, dt, G, &
                      CS%continuity_CSp, uhbt_in, vhbt_in, CS%OBC, &
                      CS%visc_rem_u, CS%visc_rem_v, u_adj, v_adj, &
                      BT_cont=CS%BT_cont)
      u_init => u_adj ; v_init => v_adj
      if (ASSOCIATED(CS%diag%du_adj2)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
        CS%diag%du_adj2(I,j,k) = u_adj(I,j,k) - u_in(I,j,k)
      enddo ; enddo ; enddo ; endif
      if (ASSOCIATED(CS%diag%dv_adj2)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
        CS%diag%dv_adj2(i,J,k) = v_adj(i,J,k) - v_in(i,J,k)
      enddo ; enddo ; enddo ; endif
    else
      call continuity(u_in, v_in, h_in, h_out, uh_in, vh_in, dt, G, &
                      CS%continuity_CSp, OBC=CS%OBC, BT_cont=CS%BT_cont)
!###   call continuity(u_in, v_in, h_in, h_out, uh_in, vh_in, dt, G, &
!###                   CS%continuity_CSp, OBC=CS%OBC, visc_rem_u=CS%visc_rem_u, &
!###                      visc_rem_v=CS%visc_rem_v, BT_cont=CS%BT_cont)
      u_init => u_in ; v_init => v_in
    endif
    call cpu_clock_end(id_clock_continuity)

    if (BT_cont_BT_thick) then
      call cpu_clock_begin(id_clock_pass)
      call pass_vector(CS%BT_cont%h_u, CS%BT_cont%h_v, G%Domain, &
                       To_All+SCALAR_PAIR, CGRID_NE)
      call cpu_clock_end(id_clock_pass)
      call btcalc(h_in, G, CS%barotropic_CSp, CS%BT_cont%h_u, CS%BT_cont%h_v)
    endif
    call cpu_clock_begin(id_clock_btstep)
    if (CS%calc_dtbt) call set_dtbt(G, CS%barotropic_CSp, eta_in, CS%pbce, &
                                    CS%BT_cont)
    call btstep(.true., uh_in, vh_in, eta_in, dt, u_bc_accel, v_bc_accel, &
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
      call continuity(u_in, v_in, h_in, h_out, uh_in, vh_in, dt, G, &
                      CS%continuity_CSp, OBC=CS%OBC, visc_rem_u=CS%visc_rem_u, &
                      visc_rem_v=CS%visc_rem_v, BT_cont=CS%BT_cont)
      call cpu_clock_end(id_clock_continuity)
      if (BT_cont_BT_thick) then
        call cpu_clock_begin(id_clock_pass)
        call pass_vector(CS%BT_cont%h_u, CS%BT_cont%h_v, G%Domain, &
                         To_All+SCALAR_PAIR, CGRID_NE)
        call cpu_clock_end(id_clock_pass)
        call btcalc(h_in, G, CS%barotropic_CSp, CS%BT_cont%h_u, CS%BT_cont%h_v)
      endif
    endif
    
    if (CS%BT_use_layer_fluxes) then
      uh_ptr => uh_in; vh_ptr => vh_in; u_ptr => u_in; v_ptr => v_in
    endif

    u_init => u_in ; v_init => v_in
    call cpu_clock_begin(id_clock_btstep)
    if (CS%calc_dtbt) call set_dtbt(G, CS%barotropic_CSp, eta_in, CS%pbce)
    call btstep(.false., u_in, v_in, eta_in, dt, u_bc_accel, v_bc_accel, &
                fluxes, CS%pbce, CS%eta_PF, u_av, v_av, CS%u_accel_bt, &
                CS%v_accel_bt, eta_pred, CS%uhbt, CS%vhbt, G, CS%barotropic_CSp,&
                CS%visc_rem_u, CS%visc_rem_v, OBC=CS%OBC, &
                BT_cont = CS%BT_cont, eta_PF_start=eta_PF_start, &
                taux_bot=taux_bot, tauy_bot=tauy_bot, &
                uh0=uh_ptr, vh0=vh_ptr, u_uh0=u_ptr, v_vh0=v_ptr)
    call cpu_clock_end(id_clock_btstep)
  endif

! u_out = u_in + dt_pred*( u_bc_accel + u_accel_bt )
  dt_pred = dt * CS%be
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v_out(i,J,k) = G%vmask(i,J) * (v_init(i,J,k) + dt_pred * &
                    (v_bc_accel(i,J,k) + CS%v_accel_bt(i,J,k)))
  enddo ; enddo ;  enddo

  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u_out(i,j,k) = G%umask(i,j) * (u_init(i,j,k) + dt_pred  * &
                    (u_bc_accel(I,j,k) + CS%u_accel_bt(I,j,k)))
  enddo ; enddo ;  enddo
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) then
    call uchksum(u_out,"Predictor 1 u",G,haloshift=0)
    call vchksum(v_out,"Predictor 1 v",G,haloshift=0)
    call hchksum(G%H_to_kg_m2*h_in,"Predictor 1 h",G,haloshift=1)
    call uchksum(G%H_to_kg_m2*uh,"Predictor 1 uh",G,haloshift=2)
    call vchksum(G%H_to_kg_m2*vh,"Predictor 1 vh",G,haloshift=2)
!   call MOM_state_chksum("Predictor 1", u_out, v_out, h_in, uh, vh, G, haloshift=1)
    call MOM_accel_chksum("Predictor accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv, &
             CS%diffu, CS%diffv, G, CS%pbce, CS%u_accel_bt, CS%v_accel_bt)
    call MOM_state_chksum("Predictor 1 init", u_init, v_init, h_in, uh, vh, G, haloshift=2)
    call check_redundant("Predictor 1 u_out", u_out, v_out, G)
    call check_redundant("Predictor 1 uh", uh, vh, G)
  endif

! u_out <- u_out + dt_pred d/dz visc d/dz u_out
! u_av  <- u_av  + dt_pred d/dz visc d/dz u_av
  call cpu_clock_begin(id_clock_vertvisc)
  call vertvisc_coef(u_out, v_out, h_in, fluxes, CS%visc, dt_pred, G, CS%vertvisc_CSp)
  call vertvisc(u_out, v_out, h_in, fluxes, CS%visc, dt_pred, CS%OBC, G, &
                CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot)
  if (G%nonblocking_updates) then
    call cpu_clock_end(id_clock_vertvisc) ; call cpu_clock_begin(id_clock_pass)
    pid_u = pass_vector_start(u_out, v_out, G%Domain)
    call cpu_clock_end(id_clock_pass) ; call cpu_clock_begin(id_clock_vertvisc)
  endif
  call vertvisc_remnant(CS%visc, CS%visc_rem_u, CS%visc_rem_v, dt_pred, G, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)

  call cpu_clock_begin(id_clock_pass)
  call pass_vector(CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                   To_All+SCALAR_PAIR, CGRID_NE)
  if (G%nonblocking_updates) then
    call pass_vector_complete(pid_u, u_out, v_out, G%Domain)
  else
    call pass_vector(u_out, v_out, G%Domain)
  endif
  call cpu_clock_end(id_clock_pass)

! uh = u_av * h_in
! h_out = h_in + dt * div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(u_out, v_out, h_in, h_out, uh, vh, dt, G, CS%continuity_CSp, &
                  CS%uhbt, CS%vhbt, CS%OBC, CS%visc_rem_u, CS%visc_rem_v, &
                  u_av, v_av, BT_cont=CS%BT_cont)
  call cpu_clock_end(id_clock_continuity)

  call cpu_clock_begin(id_clock_pass)
  call pass_var(h_out, G%Domain)
  if (G%nonblocking_updates) then
    pid_u_av = pass_vector_start(u_av, v_av, G%Domain)
    pid_uh = pass_vector_start(uh(:,:,:), vh(:,:,:), G%Domain)
  else
    call pass_vector(u_av, v_av, G%Domain, complete=.false.)
    call pass_vector(uh(:,:,:), vh(:,:,:), G%Domain)
  endif
  call cpu_clock_end(id_clock_pass)

  if (associated(CS%OBC)) then
    call Radiation_Open_Bdry_Conds(CS%OBC, u_av, u_in, v_av, v_in, &
                                   h_out, h_in, G, CS%open_boundary_CSp)
  endif

  if (CS%BT_include_udhdt) then
    call cpu_clock_begin(id_clock_continuity)
    ! Estimate u dh_dt for driving the barotropic solver.
    call continuity(u_av, v_av, h_out, h_tmp, uh_tmp, vh_tmp, dt, G, &
                    CS%continuity_CSp, OBC=CS%OBC)
    do j=js,je ; do I=Isq,Ieq ; u_dhdt(I,j) = 0.0 ; enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie ; v_dhdt(i,J) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js,je ; do I=is-1,ie
      u_dhdt(I,j) = u_dhdt(I,j) + (uh_tmp(I,j,k) - uh(I,j,k)) * Idt
    enddo ; enddo ; enddo
    do k=1,nz ; do J=js-1,je ; do i=is,ie
      v_dhdt(i,J) = v_dhdt(i,J) + (vh_tmp(i,J,k) - vh(i,J,k)) * Idt
    enddo ; enddo ; enddo
    call cpu_clock_end(id_clock_continuity)
  endif

! h_av = (h_in + h_out)/2
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    h_av(i,j,k) = 0.5*(h_in(i,j,k) + h_out(i,j,k))
  enddo ; enddo ; enddo

! The correction phase of the time step starts here.
  call enable_averaging(dt, Time_local, CS%diag)

!   Calculate a revised estimate of the free-surface height correction to be
! used in the next call to btstep.  This call is at this point so that
! h_out can be changed if CS%begw /= 0.
! eta_cor = ...                 (hidden inside CS%barotropic_CSp)
  call cpu_clock_begin(id_clock_btcalc)
  call bt_mass_source(h_out, eta_pred, fluxes, .false., dt_therm, &
                      dt_since_flux+dt, G, CS%barotropic_CSp)
  call cpu_clock_end(id_clock_btcalc)

  if (CS%begw /= 0.0) then
    ! h_out <- (1-begw)*h_in + begw*h_out
    ! Back up h_out to the value it would have had after a time-step of
    ! begw*dt.  h_out is not used again until recalculated by continuity.
    do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      h_out(i,j,k) = (1.0-CS%begw)*h_in(i,j,k) + CS%begw*h_out(i,j,k)
    enddo ; enddo ; enddo

! PFu = d/dx M(h_out,T,S)
! pbce = dM/deta
    call cpu_clock_begin(id_clock_pres)
    call PressureForce(h_out, CS%tv, CS%PFu, CS%PFv, G, &
                       CS%PressureForce_CSp, CS%regridding_opts, &
                       p_surf, CS%pbce, CS%eta_PF)
    call cpu_clock_end(id_clock_pres)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(CS%eta_PF(:,:), G%Domain)
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
    call btcalc(h_in, G, CS%barotropic_CSp, CS%BT_cont%h_u, CS%BT_cont%h_v)
  endif

  if (CS%debug) then
    call MOM_state_chksum("Predictor ", u_out, v_out, h_out, uh, vh, G)
    call uchksum(u_av,"Predictor avg u",G,haloshift=1)
    call vchksum(v_av,"Predictor avg v",G,haloshift=1)
    call hchksum(G%H_to_kg_m2*h_av,"Predictor avg h",G,haloshift=0)
  ! call MOM_state_chksum("Predictor avg ", u_av, v_av,  h_av,uh, vh, G)
    call check_redundant("Predictor u_out ", u_out, v_out, G)
    call check_redundant("Predictor u_out ", uh, vh, G)
  endif

! diffu = horizontal viscosity terms (u_av)
  call cpu_clock_begin(id_clock_horvisc)
  call horizontal_viscosity(u_av, v_av, h_av, CS%diffu, CS%diffv, &
                            CS%MEKE, CS%Varmix, G, CS%hor_visc_CSp, OBC=CS%OBC)
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
  if (CS%flux_BT_coupling .and. CS%BT_include_udhdt) then
    call btstep(.true., uh_in, vh_in, eta_in, dt, u_bc_accel, v_bc_accel, &
                fluxes, CS%pbce, CS%eta_PF, uh, vh, CS%u_accel_bt, &
                CS%v_accel_bt, eta_out, CS%uhbt, CS%vhbt, G, &
                CS%barotropic_CSp, CS%visc_rem_u, CS%visc_rem_v, etaav=eta_av, &
                uhbt_out = uhbt_out, vhbt_out = vhbt_out, OBC=CS%OBC, &
                BT_cont = CS%BT_cont, eta_PF_start = eta_PF_start, &
                sum_u_dhdt=u_dhdt, sum_v_dhdt=v_dhdt, &
                taux_bot=taux_bot, tauy_bot=tauy_bot)
  elseif (CS%flux_BT_coupling) then
    call btstep(.true., uh_in, vh_in, eta_in, dt, u_bc_accel, v_bc_accel, &
                fluxes, CS%pbce, CS%eta_PF, uh, vh, CS%u_accel_bt, &
                CS%v_accel_bt, eta_out, CS%uhbt, CS%vhbt, G, &
                CS%barotropic_CSp, CS%visc_rem_u, CS%visc_rem_v, etaav=eta_av, &
                uhbt_out = uhbt_out, vhbt_out = vhbt_out, OBC=CS%OBC, &
                BT_cont = CS%BT_cont, eta_PF_start = eta_PF_start, &
                taux_bot=taux_bot, tauy_bot=tauy_bot)
  else
    if (CS%BT_use_layer_fluxes) then
      uh_ptr => uh ; vh_ptr => vh ; u_ptr => u_av ; v_ptr => v_av
    endif

    call btstep(.false., u_in, v_in, eta_in, dt, u_bc_accel, v_bc_accel, &
                fluxes, CS%pbce, CS%eta_PF, u_av, v_av, CS%u_accel_bt, &
                CS%v_accel_bt, eta_out, CS%uhbt, CS%vhbt, G, &
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

! u_out = u_in + dt*( u_bc_accel + u_accel_bt )
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u_out(i,j,k) = G%umask(i,j) * (u_init(i,j,k) + dt * &
                    (u_bc_accel(I,j,k) + CS%u_accel_bt(I,j,k)))
  enddo ; enddo ; enddo
  if (ASSOCIATED(CS%diag%PFu_tot)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      CS%diag%PFu_tot(i,j,k) = CS%PFu(i,j,k)
  enddo ; enddo ; enddo ; endif
  if (ASSOCIATED(CS%diag%CAu_tot)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      CS%diag%CAu_tot(i,j,k) = CS%CAu(i,j,k) !+ CS%u_accel_bt(i,j) - CS%diag%PFu_bt(i,j)
  enddo ; enddo ; enddo ; endif

  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v_out(i,j,k) = G%vmask(i,j) * (v_init(i,j,k) + dt * &
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
    call uchksum(u_out,"Corrector 1 u",G,haloshift=0)
    call vchksum(v_out,"Corrector 1 v",G,haloshift=0)
    call hchksum(G%H_to_kg_m2*h_in,"Corrector 1 h",G,haloshift=2)
    call uchksum(G%H_to_kg_m2*uh,"Corrector 1 uh",G,haloshift=2)
    call vchksum(G%H_to_kg_m2*vh,"Corrector 1 vh",G,haloshift=2)
  ! call MOM_state_chksum("Corrector 1", u_out, v_out, h_in, uh, vh, G, haloshift=1)
    call MOM_accel_chksum("Corrector accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv, &
             CS%diffu, CS%diffv, G, CS%pbce, CS%u_accel_bt, CS%v_accel_bt)
  endif

! u_out <- u_out + dt d/dz visc d/dz u_out
! u_av  <- u_av  + dt d/dz visc d/dz u_av
  call cpu_clock_begin(id_clock_vertvisc)
  call vertvisc_coef(u_out, v_out, h_in, fluxes, CS%visc, dt, G, CS%vertvisc_CSp)
  call vertvisc(u_out, v_out, h_in, fluxes, CS%visc, dt, CS%OBC, G, &
                CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot)
  if (G%nonblocking_updates) then
    call cpu_clock_end(id_clock_vertvisc) ; call cpu_clock_begin(id_clock_pass)
    pid_u = pass_vector_start(u_out, v_out, G%Domain)
    call cpu_clock_end(id_clock_pass) ; call cpu_clock_begin(id_clock_vertvisc)
  endif
  call vertvisc_remnant(CS%visc, CS%visc_rem_u, CS%visc_rem_v, dt, G, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)

  call cpu_clock_begin(id_clock_pass)
  call pass_vector(CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                   To_All+SCALAR_PAIR, CGRID_NE)
  if (G%nonblocking_updates) then
    call pass_vector_complete(pid_u, u_out, v_out, G%Domain)
  else
    call pass_vector(u_out, v_out, G%Domain)
  endif
  call cpu_clock_end(id_clock_pass)

! uh = u_av * h_in
! h_out = h_in + dt * div . uh
  if (CS%flux_BT_coupling) then
    ! u_av and v_av adjusted so their mass transports match uhbt and vhbt.
    ! Also, determine the values of u_out and v_out so that their transports
    ! that agree with uhbt_out and vhbt_out.
    if (ASSOCIATED(CS%diag%du_adj)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      CS%diag%du_adj(I,j,k) = u_out(I,j,k)
    enddo ; enddo ; enddo ; endif
    if (ASSOCIATED(CS%diag%dv_adj)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      CS%diag%dv_adj(i,J,k) = v_out(i,J,k)
    enddo ; enddo ; enddo ; endif
    if (CS%BT_include_udhdt) then
      do j=js,je ; do I=is-1,ie ; uhbt_in(I,j) = uhbt_out(I,j) ; enddo ; enddo
      do J=js-1,je ; do i=is,ie ; vhbt_in(i,J) = vhbt_out(i,J) ; enddo ; enddo
      CS%readjust_velocity = .true.
    endif
    call cpu_clock_begin(id_clock_continuity)
    call continuity(u_out, v_out, h_in, h_out, uh, vh, dt, G, &
                    CS%continuity_CSp, CS%uhbt, CS%vhbt, CS%OBC, &
                    CS%visc_rem_u, CS%visc_rem_v, u_av, v_av, &
                    uhbt_out, vhbt_out, u_out, v_out)
    call cpu_clock_end(id_clock_continuity)
    if (G%nonblocking_updates) then
      call cpu_clock_begin(id_clock_pass)
      pid_h = pass_var_start(h_out, G%Domain)
      call cpu_clock_end(id_clock_pass)
    endif
    if (ASSOCIATED(CS%diag%du_adj)) then ; do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      CS%diag%du_adj(I,j,k) = u_out(I,j,k) - CS%diag%du_adj(I,j,k)
    enddo ; enddo ; enddo ; endif
    if (ASSOCIATED(CS%diag%dv_adj)) then ; do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      CS%diag%dv_adj(i,J,k) = v_out(i,J,k) - CS%diag%dv_adj(i,J,k)
    enddo ; enddo ; enddo ; endif

    call cpu_clock_begin(id_clock_vertvisc)
    call vertvisc_limit_vel(u_out, v_out, h_in, fluxes, CS%visc, dt, G, CS%vertvisc_CSp)
    if (G%nonblocking_updates) then
      call cpu_clock_end(id_clock_vertvisc) ; call cpu_clock_begin(id_clock_pass)
      pid_u = pass_vector_start(u_out, v_out, G%Domain)
      call cpu_clock_end(id_clock_pass) ; call cpu_clock_begin(id_clock_vertvisc)
    endif
    call vertvisc_limit_vel(u_av, v_av, h_in, fluxes, CS%visc, dt, G, CS%vertvisc_CSp)
    call cpu_clock_end(id_clock_vertvisc)

    call cpu_clock_begin(id_clock_pass)
    if (G%nonblocking_updates) then
      call pass_var_complete(pid_h, h_out, G%Domain)
      call pass_vector_complete(pid_u, u_out, v_out, G%Domain)
    else
      call pass_var(h_out, G%Domain)
      call pass_vector(u_out, v_out, G%Domain, complete=.false.)
    endif
    call cpu_clock_end(id_clock_pass)
  else
    ! u_av and v_av adjusted so their mass transports match uhbt and vhbt.
    call cpu_clock_begin(id_clock_continuity)
    call continuity(u_out, v_out, h_in, h_out, uh, vh, dt, G, &
                    CS%continuity_CSp, CS%uhbt, CS%vhbt, CS%OBC, &
                    CS%visc_rem_u, CS%visc_rem_v, u_av, v_av)
    call cpu_clock_end(id_clock_continuity)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(h_out, G%Domain)
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
    call Radiation_Open_Bdry_Conds(CS%OBC, u_out, u_in, v_out, v_in, &
                                   h_out, h_in, G, CS%open_boundary_CSp)
  endif

! h_av = (h_in + h_out)/2
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    h_av(i,j,k) = 0.5*(h_in(i,j,k) + h_out(i,j,k))
  enddo ; enddo ; enddo

  if (G%nonblocking_updates) then
    call cpu_clock_begin(id_clock_pass)
    call pass_vector_complete(pid_uh, uh(:,:,:), vh(:,:,:), G%Domain)
    call pass_vector_complete(pid_u_av, u_av, v_av, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  do k=1,nz ; do j=js-2,je+2 ; do I=Isq-2,Ieq+2
    CS%uhtr(I,j,k) = CS%uhtr(I,j,k) + uh(I,j,k)*dt
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq-2,Jeq+2 ; do i=is-2,ie+2
    CS%vhtr(i,J,k) = CS%vhtr(i,J,k) + vh(i,J,k)*dt
  enddo ; enddo ; enddo

  if (CS%thickness_diffuse .and. .not.CS%thickness_diffuse_first) then
    call cpu_clock_begin(id_clock_thick_diff)
    if (associated(CS%VarMix)) &
      call calc_slope_function(h_out, CS%tv, G, CS%VarMix)
    call thickness_diffuse(h_out, CS%uhtr, CS%vhtr, CS%tv, dt, G, &
                           CS%MEKE, CS%VarMix, CS%thickness_diffuse_CSp)
    call cpu_clock_end(id_clock_thick_diff)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(h_out, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  if (CS%mixedlayer_restrat) then
    call cpu_clock_begin(id_clock_ml_restrat)
    call mixedlayer_restrat(h_out, CS%uhtr,CS%vhtr,CS%tv, fluxes, dt, &
                            G, CS%mixedlayer_restrat_CSp)
    call cpu_clock_end(id_clock_ml_restrat)
    call cpu_clock_begin(id_clock_pass)
    call pass_var(h_out, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  if (associated(CS%MEKE)) then
    call step_forward_MEKE(CS%MEKE, h_out, CS%visc, dt, G, CS%MEKE_CSp)
  endif

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
  if (CS%BT_include_udhdt) then
    if (CS%id_h_dudt > 0) call post_data(CS%id_h_dudt, u_dhdt, CS%diag)
    if (CS%id_h_dvdt > 0) call post_data(CS%id_h_dvdt, v_dhdt, CS%diag)
  endif
  if (CS%debug) then
    call MOM_state_chksum("Corrector ", u_out, v_out, h_out, uh, vh, G)
    call uchksum(u_av,"Corrector avg u",G,haloshift=1)
    call vchksum(v_av,"Corrector avg v",G,haloshift=1)
    call hchksum(G%H_to_kg_m2*h_av,"Corrector avg h",G,haloshift=1)
 !  call MOM_state_chksum("Corrector avg ", u_av, v_av, h_av, uh, vh, G)
  endif

end subroutine step_MOM_dyn_split_RK2

! =============================================================================

subroutine register_restarts_dyn_split_RK2(G, param_file, CS, restart_CS)
  type(ocean_grid_type),     intent(in)    :: G
  type(param_file_type),     intent(in)    :: param_file
  type(MOM_control_struct), intent(inout) :: CS
  type(MOM_restart_CS),     pointer       :: restart_CS
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
! if (associated(CS_split)) then
!   call MOM_error(WARNING, "register_restarts_split_RK2 called with an associated "// &
!                            "control structure.")
!   return
! endif
! allocate(CS_split)

  thickness_units = get_thickness_units(G)
  flux_units = get_flux_units(G)

 ! if (G%Boussinesq) then
    vd = vardesc("sfc","Free surface Height",'h','1','s',thickness_units, 'd')
 ! else
 !   vd(1) = vardesc("ocean_mass?","Ocean column mass",'h','1','s',"kg meter-2", 'd')
 ! endif
  call register_restart_field(CS%eta, CS%eta, vd, .false., CS%restart_CSp)

  vd = vardesc("u2","Auxiliary Zonal velocity",'u','L','s',"meter second-1", 'd')
  call register_restart_field(CS%u_av, CS%u_av, vd, .false., CS%restart_CSp)

  vd = vardesc("v2","Auxiliary Meridional velocity",'v','L','s',"meter second-1", 'd')
  call register_restart_field(CS%v_av, CS%v_av, vd, .false., CS%restart_CSp)

  vd = vardesc("h2","Auxiliary Layer Thickness",'h','L','s',thickness_units, 'd')
  call register_restart_field(CS%h_av, CS%h_av, vd, .false., CS%restart_CSp)

  vd = vardesc("uh","Zonal thickness flux",'u','L','s',flux_units, 'd')
  call register_restart_field(CS%uh, CS%uh, vd, .false., CS%restart_CSp)

  vd = vardesc("vh","Meridional thickness flux",'v','L','s',flux_units, 'd')
  call register_restart_field(CS%vh, CS%vh, vd, .false., CS%restart_CSp)

  vd = vardesc("diffu","Zonal horizontal viscous acceleration",'u','L','s', &
               "meter second-2", 'd')
  call register_restart_field(CS%diffu, CS%diffu, vd, .false., CS%restart_CSp)

  vd = vardesc("diffv","Meridional horizontal viscous acceleration",'v','L','s',&
               "meter second-2", 'd')
  call register_restart_field(CS%diffv, CS%diffv, vd, .false., CS%restart_CSp)

  call register_barotropic_restarts(G, param_file, CS%barotropic_CSp, &
                                    CS%restart_CSp)

  if (CS%readjust_bt_trans) then
    vd = vardesc("uhbt_in","Final instantaneous barotropic zonal thickness flux",&
                 'u','1','s',flux_units, 'd')
    call register_restart_field(CS%uhbt_in, CS%uhbt_in, vd, .false., CS%restart_CSp)

    vd = vardesc("vhbt_in","Final instantaneous barotropic meridional thickness flux",&
                 'v','1','s',flux_units, 'd')
    call register_restart_field(CS%vhbt_in, CS%vhbt_in, vd, .false., CS%restart_CSp)
  endif

end subroutine register_restarts_dyn_split_RK2

subroutine initialize_dyn_split_RK2(u, v, h, Time, G, param_file, diag, CS, restart_CS)
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
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_tmp
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

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

  if (CS%thickness_diffuse) &
    id_clock_thick_diff = cpu_clock_id('(Ocean thickness diffusion)', grain=CLOCK_MODULE)
  if (CS%mixedlayer_restrat) &
    id_clock_ml_restrat = cpu_clock_id('(Ocean mixed layer restrat)', grain=CLOCK_MODULE)

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

  call barotropic_init(u, v, h, CS%eta(:,:), Time, G, &
                       param_file, diag, CS%barotropic_CSp, restart_CS, &
                       CS%BT_cont, CS%tides_CSp)

  if (.not. query_initialized(CS%diffu,"diffu",restart_CS) .or. &
      .not. query_initialized(CS%diffv,"diffv",restart_CS)) &
    call horizontal_viscosity(u, v, h, CS%diffu, CS%diffv, CS%MEKE, CS%VarMix, &
                              G, CS%hor_visc_CSp)
  if (.not. query_initialized(CS%u_av,"u2", restart_CS) .or. &
      .not. query_initialized(CS%u_av,"v2", restart_CS)) then
    CS%u_av(:,:,:) = u(:,:,:)
    CS%v_av(:,:,:) = v(:,:,:)
  endif
! This call is just here to initialize uh and vh.
  if (.not. query_initialized(CS%uh,"uh",restart_CS) .or. &
      .not. query_initialized(CS%vh,"vh",restart_CS)) then
    h_tmp(:,:,:) = h(:,:,:)
    call continuity(u, v, h, h_tmp, CS%uh, CS%vh, &
                    CS%dt, G, CS%continuity_CSp, OBC=CS%OBC)
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
  if (CS%readjust_BT_trans .or. CS%BT_include_udhdt) then
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
  call pass_vector(CS%uh, CS%vh, G%Domain)
  call cpu_clock_end(id_clock_pass_init)

end subroutine initialize_dyn_split_RK2

end module MOM_dynamics_split_RK2
