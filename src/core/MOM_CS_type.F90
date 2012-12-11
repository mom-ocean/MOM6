module MOM_CS_type

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

!   This module provdes a control stucture (collection of types and pointers)
! that is used throughout the main high-level algorithm (i.e., the code in
! MOM.F90) and is shared between the variants on the algorithm (e.g. split RK2
! and unsplit RK3).

use MOM_variables, only : directories, vertvisc_type, ocean_OBC_type
use MOM_variables, only : BT_cont_type
use MOM_variables, only : &
  forcing, &      ! A structure containing pointers to the forcing fields
                  ! which may be used to drive MOM.  All fluxes are
                  ! positive downward.
  surface, &      ! A structure containing pointers to various fields which
                  ! may be used describe the surface state of MOM, and
                  ! which will be returned to the calling program
  thermo_var_ptrs, & ! A structure containing pointers to an assortment of
                  ! thermodynamic fields that may be available, including
                  ! potential temperature, salinity and mixed layer density.
  ocean_internal_state  ! A structure containing pointers to most of the above.

use MOM_checksums, only : hchksum, uchksum, vchksum
use MOM_diag_mediator, only : diag_ptrs
use MOM_restart, only : MOM_restart_CS
use MOM_time_manager, only : time_type

use MOM_barotropic, only : barotropic_CS
use MOM_continuity, only : continuity_CS
use MOM_CoriolisAdv, only : CoriolisAdv_CS
use MOM_diabatic_driver, only : diabatic_CS
use MOM_diagnostics, only : diagnostics_CS
use MOM_diag_to_Z, only : diag_to_Z_CS
use MOM_grid, only : ocean_grid_type
use MOM_hor_visc, only : hor_visc_CS
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE, only : MEKE_CS
use MOM_MEKE_types, only : MEKE_type
use MOM_mixed_layer_restrat, only : mixedlayer_restrat_CS
use MOM_open_boundary, only : open_boundary_CS
use MOM_PressureForce, only : PressureForce_CS
use MOM_thickness_diffuse, only : thickness_diffuse_CS
use MOM_tidal_forcing, only : tidal_forcing_CS
use MOM_tracer, only : advect_tracer_CS
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_vert_friction, only : vertvisc_CS
use MOM_set_visc, only : set_visc_CS
use regrid_defs, only: regridding_opts_t

implicit none ; private

public MOM_state_chksum, MOM_thermo_chksum, MOM_accel_chksum

#include <MOM_memory.h>

type, public :: MOM_control_struct
  real PTR_, dimension(NXMEMQP_,NYMEM_,NKMEM_,C2_) :: &
    u         ! Zonal velocity, in m s-1.
  real PTR_, dimension(NXMEM_,NYMEMQP_,NKMEM_,C2_) :: &
    v         ! Meridional velocity, in m s-1.
  real PTR_, dimension(NXMEM_,NYMEM_,NKMEM_,C2_) :: &
    h         ! Layer thickness, in m or kg m-2 (H).
  real PTR_, dimension(NXMEM_,NYMEM_) :: &
    eta       ! Instantaneous free surface height, in m.
  real PTR_, dimension(NXMEM_,NYMEM_,NKMEM_) :: &
    T, &      ! Potential temperature in C.
    S         ! Salinity in PSU.
  real PTR_, dimension(NXMEM_,NYMEM_,NKMEM_) :: &
    h_aux     ! Work array for remapping (same units as h).
  real PTR_, dimension(NXMEMQP_,NYMEM_,NKMEM_) :: &
    uh, &     ! uh = u * h * dy at u grid points in m3 s-1.
    CAu, &    ! CAu = f*v - u.grad(u) in m s-2.
    PFu, &    ! PFu = -dM/dx, in m s-2.
    diffu, &  ! Zonal acceleration due to convergence of the along-isopycnal
              ! stress tensor, in m s-2.
    visc_rem_u, & ! Both the fraction of the zonal momentum originally in a
              ! layer that remains after a time-step of viscosity, and the
              ! fraction of a time-step's worth of a barotropic acceleration
              ! that a layer experiences after viscosity is applied.
              ! Nondimensional between 0 (at the bottom) and 1 (far above).
    uhtr      ! Accumlated zonal thickness fluxes used to advect tracers, in m3.
  real PTR_, dimension(NXMEM_,NYMEMQP_,NKMEM_) :: &
    vh, &     ! vh = v * h * dx at v grid points in m3 s-1.
    CAv, &    ! CAv = -f*u - u.grad(v) in m s-2.
    PFv, &    ! PFv = -dM/dy, in m s-2.
    diffv, &  ! Meridional acceleration due to convergence of the
              ! along-isopycnal stress tensor, in m s-2.
    visc_rem_v, & ! Both the fraction of the meridional momentum originally in
              ! a layer that remains after a time-step of viscosity, and the
              ! fraction of a time-step's worth of a barotropic acceleration
              ! that a layer experiences after viscosity is applied.
              ! Nondimensional between 0 (at the bottom) and 1 (far above).
    vhtr      ! Accumlated meridional thickness fluxes used to advect tracers, in m3.
  real PTR_, dimension(NXMEM_,NYMEM_) :: &
    ave_ssh, &! The time-averaged sea surface height in m.
    eta_PF    ! The instantaneous SSH used in calculating PFu and PFv, in m.
  real PTR_, dimension(NXMEMQP_,NYMEM_) :: uhbt
  real PTR_, dimension(NXMEM_,NYMEMQP_) :: vhbt
    ! uhbt and vhbt are the average volume or mass fluxes determined by the
    ! barotropic solver in m3 s-1 or kg s-1.  uhbt and vhbt should (roughly?) 
    ! equal the verticals sum of uh and vh, respectively.
  real PTR_, dimension(NXMEMQP_,NYMEM_) :: uhbt_in
  real PTR_, dimension(NXMEM_,NYMEMQP_) :: vhbt_in
    ! uhbt_in and vhbt_in are the vertically summed transports from based on
    ! the final thicknessses and velocities from the previous dynamics time
    ! step, both in units of m3 s-1 or kg s-1.
! The following 6 variables are only used with the split time stepping scheme.
  real PTR_, dimension(NXMEM_,NYMEM_,NKMEM_) :: pbce
      ! pbce times eta gives the baroclinic pressure anomaly in each layer due
      ! to free surface height anomalies.  pbce has units of m2 H-1 s-2.
  real PTR_, dimension(NXMEMQP_,NYMEM_,NKMEM_) :: u_accel_bt
  real PTR_, dimension(NXMEM_,NYMEMQP_,NKMEM_) :: v_accel_bt
    ! u_accel_bt and v_accel_bt are layers' accelerations due to
    ! the difference between the accelerations from the barotropic calculation
    ! baroclinic accelerations that were fed into the barotropic
    ! calculation, in m s-2.
  real PTR_, dimension(NXMEMQP_,NYMEM_,NKMEM_) :: u_av
  real PTR_, dimension(NXMEM_,NYMEMQP_,NKMEM_) :: v_av
    ! u_av and v_av are the layer velocities with the vertical mean replaced by
    ! the time-mean barotropic velocity over a baroclinic timestep, in m s-1.
  real PTR_, dimension(NXMEM_,NYMEM_,NKMEM_)  :: h_av
    ! The arithmetic mean of two successive layer thicknesses, in m or kg m-2.
  real, pointer, dimension(:,:) :: taux_bot => NULL(), tauy_bot => NULL()
    ! The frictional bottom stresses from the ocean to the seafloor, in Pa.
  real, pointer, dimension(:,:,:) :: &
    u_prev => NULL(), &  ! The previous values of u and v, stored for
    v_prev => NULL()     ! diagnostic purposes.

  type(ocean_grid_type) :: grid ! A structure containing metrics and grid info.
  type(thermo_var_ptrs) :: tv ! A structure containing pointers to an assortment
                              ! of thermodynamic fields that may be available.
  type(diag_ptrs) :: diag     ! A structure containing pointers to
                              ! diagnostic fields that might be calculated
                              ! and shared between modules.
  type(vertvisc_type) :: visc ! A structure containing vertical viscosities,
                              ! bottom drag viscosities, and related fields.
  type(BT_cont_type), pointer :: BT_cont => NULL()
                              ! A structure with elements that describe the
                              ! effective summed open face areas as a function
                              ! of barotropic flow.
  type(MEKE_type), pointer :: MEKE => NULL()  ! A structure containing fields
                              ! related to the Mesoscale Eddy Kinetic Energy.

  logical :: split           ! If true, use the split time stepping scheme.
  logical :: use_RK2         ! If true, use RK2 instead of RK3 for time-
                             ! stepping in unsplit mode.
  logical :: split_bottom_stress  ! If true, provide the bottom stress
                             ! calculated by the vertical viscosity to the
                             ! barotropic solver.
  logical :: adiabatic       ! If true, there are no diapycnal mass fluxes, and
                             ! the subroutine calls to calculate and apply such
                             ! diapycnal fluxes are eliminated.
  logical :: use_temperature ! If true, temperature and salinity are used as
                             ! state variables.
  logical :: use_frazil      ! If true, water freezes if it gets too cold, and
                             ! the accumulated heat deficit is returned in the
                             ! surface state.
  logical :: bound_salinity  ! If true, salt is added to keep the salinity above
                             ! a minimum value, and the deficit is reported.
  logical :: bulkmixedlayer  ! If true, a refined bulk mixed layer is used with
                             ! nkml sublayers and nkbl buffer layer.
  logical :: thickness_diffuse ! If true, interfaces are diffused with a
                             ! coefficient of KHTH.
  logical :: thickness_diffuse_first ! If true, diffuse thickness before dynamics.
  logical :: mixedlayer_restrat ! If true, a density-gradient dependent
                             ! restratifying flow is imposed in the mixed layer.
  logical :: debug           ! If true, write verbose checksums for debugging purposes.
  logical :: debug_truncations  ! If true, make sure that all diagnostics that
                             ! could be useful for debugging any truncations are
                             ! calculated.

  real    :: dt              ! The (baroclinic) dynamics time step, in s.
  real    :: dt_therm        ! The thermodynamics time step, in s.
  real    :: be              ! A nondimensional number from 0.5 to 1 that controls
                             ! the backward weighting of the time stepping scheme.
  real    :: begw            ! A nondimensional number from 0 to 1 that controls
                             ! the extent to which the treatment of gravity waves
                             ! is forward-backward (0) or simulated backward
                             ! Euler (1).  0 is almost always used.
  type(time_type) :: Z_diag_interval  !   The amount of time between calls to
                             ! calculate Z-space diagnostics.
  type(time_type) :: Z_diag_time  ! The next time at which Z-space diagnostics
                             ! should be calculated.
  real    :: Hmix            ! The diagnostic mixed layer thickness in m when
                             ! the bulk mixed layer is not used.
  logical :: calc_bbl ! If true, the BBL viscosity and thickness need to be
                      ! calculated. This only applies with BOTTOMDRAGLAW true.
  real :: bbl_calc_time_interval ! The amount of time to use in diagnostics of
                                 ! the BBL properties.
  real :: missing=-1.0e34    ! The missing data value for masked fields.

  integer :: ntrunc          ! The number of times the velocity has been
                             ! truncated since the last call to write_energy.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(ocean_OBC_type), pointer :: OBC => NULL() ! A pointer to an open boundary
                             ! condition type that specifies whether, where, and
                             ! what open boundary conditions are used.  If no
                             ! open BCs are used, this pointer stays nullified.
  real :: rel_time = 0.0     ! Relative time in s since the start
                             ! of the current execution.

  ! This is to allow the previous, velocity-based coupling with between the
  ! baroclinic and barotropic modes.
  logical :: interp_p_surf     ! If true, linearly interpolate the surface
                               ! pressure over the coupling time step, using 
                               ! the specified value at the end of the coupling
                               ! step. False by default.
  real    :: smooth_ssh_passes ! If greater than 0, apply this number of 
                               ! spatial smoothing passes to the sea surface
                               ! heights that are reported back to the calling
                               ! program by MOM.  The default is 0.
  logical :: p_surf_prev_set   ! If true, p_surf_prev has been properly set from
                               ! a previous time-step or the ocean restart file.
                               ! This is only valid when interp_p_surf is true.
  logical :: flux_BT_coupling  ! If true, use volume fluxes, not velocities,
                               ! to couple the baroclinic and barotropic modes.
  logical :: readjust_BT_trans ! If true, readjust the barotropic transport of
                               ! the input velocities to agree with CS%uhbt_in
                               ! and CS%vhbt_in after the diabatic step.
  logical :: BT_include_udhdt  ! If true, estimate the sum of u dh/dt and v dh/dt
                               ! and provide them to the barotropic solver.
  logical :: BT_use_layer_fluxes ! If true, use the summed layered fluxes plus
                               ! an adjustment due to a changed barotropic
                               ! velocity in the barotropic continuity equation.
  real    :: dtbt_reset_period ! The time interval in seconds between dynamic
                               ! recalculation of the barotropic time step.  If
                               ! this is negative, it is never calculated, and
                               ! if it is 0, it is calculated every step.
  logical :: calc_dtbt         ! If true, calculate the barotropic time-step
                               ! dynamically.
  logical :: readjust_velocity ! A flag that varies with time that determines
                               ! whether the velocities currently need to be
                               ! readjusted to agree with CS%uhbt_in and
                               ! CS%vhbt_in.  This is only used if 
                               ! CS%readjust_BT_trans or BT_include_udhdt are true.
  logical :: check_bad_surface_vals ! If true, scans the surface state for
                                    ! ridiculous values
  real    :: bad_val_ssh_max   ! Maximum SSH before triggering bad value message
  real    :: bad_val_sst_max   ! Maximum SST before triggering bad value message
  real    :: bad_val_sst_min   ! Minimum SST before triggering bad value message
  real    :: bad_val_sss_max   ! Maximum SSS before triggering bad value message

  real, pointer, dimension(:,:) :: &
    p_surf_prev, &  ! The value of the surface pressure at the end of the
                    ! previous call to step_MOM, in Pa.
    p_surf_begin, & ! The values of the surface pressure at the beginning and
    p_surf_end      ! end of a call to step_MOM_dyn_..., in Pa.

  ! Arrays that can be used to store advective and diffusive tracer fluxes.
  real, pointer, dimension(:,:,:) :: &
    T_adx => NULL(), T_ady => NULL(), T_diffx => NULL(), T_diffy => NULL(), &
    S_adx => NULL(), S_ady => NULL(), S_diffx => NULL(), S_diffy => NULL()
  ! Arrays that can be used to store vertically integrated advective and
  ! diffusive tracer fluxes.
  real, pointer, dimension(:,:) :: &
    T_adx_2d => NULL(), T_ady_2d => NULL(), T_diffx_2d => NULL(), T_diffy_2d => NULL(), &
    S_adx_2d => NULL(), S_ady_2d => NULL(), S_diffx_2d => NULL(), S_diffy_2d => NULL(), &
    SST_sq => NULL()

! The following are the ids of various diagnostics.
  integer :: id_u = -1, id_v = -1, id_h = -1, id_uh = -1, id_vh = -1
  integer :: id_uav = -1, id_vav = -1
  integer :: id_T = -1, id_S = -1, id_ssh = -1, id_fraz = -1
  integer :: id_salt_deficit = -1, id_Heat_PmE = -1, id_intern_heat = -1
  integer :: id_du_adj = -1, id_dv_adj = -1, id_du_adj2 = -1, id_dv_adj2 = -1
  integer :: id_h_dudt = -1, id_h_dvdt = -1
  integer :: id_sst = -1, id_sst_sq = -1, id_sss = -1, id_ssu = -1, id_ssv = -1
  integer :: id_speed = -1, id_ssh_inst = -1

  integer :: id_PFu = -1, id_PFv = -1, id_CAu = -1, id_CAv = -1
  integer :: id_u_BT_accel = -1, id_v_BT_accel = -1
  integer :: id_Tadx = -1, id_Tady = -1, id_Tdiffx = -1, id_Tdiffy = -1
  integer :: id_Sadx = -1, id_Sady = -1, id_Sdiffx = -1, id_Sdiffy = -1
  integer :: id_Tadx_2d = -1, id_Tady_2d = -1, id_Tdiffx_2d = -1, id_Tdiffy_2d = -1
  integer :: id_Sadx_2d = -1, id_Sady_2d = -1, id_Sdiffx_2d = -1, id_Sdiffy_2d = -1
  integer :: id_u_predia = -1, id_v_predia = -1, id_h_predia = -1
  integer :: id_T_predia = -1, id_S_predia = -1, id_e_predia = -1
! The remainder of the structure is pointers to child subroutines' control strings.
  type(hor_visc_CS), pointer :: hor_visc_CSp => NULL()
  type(continuity_CS), pointer :: continuity_CSp => NULL()
  type(CoriolisAdv_CS), pointer :: CoriolisAdv_CSp => NULL()
  type(PressureForce_CS), pointer :: PressureForce_CSp => NULL()
  type(barotropic_CS), pointer :: barotropic_CSp => NULL()
  type(vertvisc_CS), pointer :: vertvisc_CSp => NULL()
  type(set_visc_CS), pointer :: set_visc_CSp => NULL()
  type(diabatic_CS), pointer :: diabatic_CSp => NULL()
  type(thickness_diffuse_CS), pointer :: thickness_diffuse_CSp => NULL()
  type(mixedlayer_restrat_CS), pointer :: mixedlayer_restrat_CSp => NULL()
  type(open_boundary_CS), pointer :: open_boundary_CSp => NULL()
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL()
  type(MEKE_CS),  pointer :: MEKE_CSp => NULL()
  type(VarMix_CS),  pointer :: VarMix => NULL()
  type(advect_tracer_CS), pointer :: tracer_CSp => NULL()
  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp => NULL()
  type(diagnostics_CS), pointer :: diagnostics_CSp => NULL()
  type(diag_to_Z_CS), pointer :: diag_to_Z_CSp => NULL()
  type(MOM_restart_CS),  pointer :: restart_CSp => NULL()
  type(regridding_opts_t) :: regridding_opts ! Why is this not a pointer? - AJA
end type MOM_control_struct

contains

! =============================================================================

subroutine MOM_state_chksum(mesg, u, v, h, uh, vh, G, haloshift)
  character(len=*),                       intent(in) :: mesg
  real, dimension(NXMEMQ_,NYMEM_,NKMEM_), intent(in) :: u
  real, dimension(NXMEM_,NYMEMQ_,NKMEM_), intent(in) :: v
  real, dimension(NXMEM_,NYMEM_,NKMEM_),  intent(in) :: h
  real, dimension(NXMEMQ_,NYMEM_,NKMEM_), intent(in) :: uh
  real, dimension(NXMEM_,NYMEMQ_,NKMEM_), intent(in) :: vh
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
  real, dimension(NXMEMQ_,NYMEM_,NKMEM_), intent(in) :: CAu
  real, dimension(NXMEM_,NYMEMQ_,NKMEM_), intent(in) :: CAv
  real, dimension(NXMEMQ_,NYMEM_,NKMEM_), intent(in) :: PFu
  real, dimension(NXMEM_,NYMEMQ_,NKMEM_), intent(in) :: PFv
  real, dimension(NXMEMQ_,NYMEM_,NKMEM_), intent(in) :: diffu
  real, dimension(NXMEM_,NYMEMQ_,NKMEM_), intent(in) :: diffv
  type(ocean_grid_type),                  intent(in) :: G
  real, dimension(NXMEM_,NYMEM_,NKMEM_),  optional, intent(in) :: pbce
  real, dimension(NXMEMQ_,NYMEM_,NKMEM_), optional, intent(in) :: u_accel_bt
  real, dimension(NXMEM_,NYMEMQ_,NKMEM_), optional, intent(in) :: v_accel_bt
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

end module MOM_CS_type
