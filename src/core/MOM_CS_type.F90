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

use MOM_variables, only : vertvisc_type, ocean_OBC_type
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
use MOM_ALE, only : ALE_CS

implicit none ; private

#include <MOM_memory.h>

type, public :: MOM_control_struct
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: &
    h, &      ! Layer thickness, in m or kg m-2 (H).
    T, &      ! Potential temperature in C.
    S, &      ! Salinity in PSU.
    h_aux     ! Work array for remapping (same units as h).
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    u, &      ! Zonal velocity, in m s-1.
    uh, &     ! uh = u * h * dy at u grid points in m3 s-1.
    uhtr      ! Accumlated zonal thickness fluxes used to advect tracers, in m3.
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    v, &      ! Meridional velocity, in m s-1.
    vh, &     ! vh = v * h * dx at v grid points in m3 s-1.
    vhtr      ! Accumlated meridional thickness fluxes used to advect tracers, in m3.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    ave_ssh   ! The time-averaged sea surface height in m.

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
  type(MEKE_type), pointer :: MEKE => NULL()  ! A structure containing fields
                              ! related to the Mesoscale Eddy Kinetic Energy.

  logical :: split           ! If true, use the split time stepping scheme.
  logical :: use_RK2         ! If true, use RK2 instead of RK3 for time-
                             ! stepping in unsplit mode.
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
  logical :: useALEalgorithm ! If true, use the ALE algorithm rather than layered
                             ! isopycnal/stacked shallow water mode. This logical is
                             ! set by calling the function useRegridding() from the
                             ! MOM_regridding module.

  real    :: dt              ! The (baroclinic) dynamics time step, in s.
  real    :: dt_therm        ! The thermodynamics time step, in s.
  type(time_type) :: Z_diag_interval  !   The amount of time between calls to
                             ! calculate Z-space diagnostics.
  type(time_type) :: Z_diag_time  ! The next time at which Z-space diagnostics
                             ! should be calculated.
  real    :: Hmix            ! The diagnostic mixed layer thickness in m when
                             ! the bulk mixed layer is not used.
  real :: missing=-1.0e34    ! The missing data value for masked fields.

  integer :: ntrunc          ! The number of times the velocity has been
                             ! truncated since the last call to write_energy.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  real :: rel_time = 0.0     ! Relative time in s since the start
                             ! of the current execution.

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

  real    :: dtbt_reset_period ! The time interval in seconds between dynamic
                               ! recalculation of the barotropic time step.  If
                               ! this is negative, it is never calculated, and
                               ! if it is 0, it is calculated every step.

  logical :: check_bad_surface_vals ! If true, scans the surface state for
                                    ! ridiculous values
  real    :: bad_val_ssh_max   ! Maximum SSH before triggering bad value message
  real    :: bad_val_sst_max   ! Maximum SST before triggering bad value message
  real    :: bad_val_sst_min   ! Minimum SST before triggering bad value message
  real    :: bad_val_sss_max   ! Maximum SSS before triggering bad value message

  real, pointer, dimension(:,:) :: &
    p_surf_prev => NULL(), &  ! The value of the surface pressure at the end of
                              ! the previous call to step_MOM, in Pa.
    p_surf_begin => NULL(), & ! The values of the surface pressure at the start
    p_surf_end => NULL()      ! and end of a call to step_MOM_dyn_..., in Pa.

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
  integer :: id_u = -1, id_v = -1, id_h = -1
  integer :: id_T = -1, id_S = -1, id_ssh = -1, id_fraz = -1
  integer :: id_salt_deficit = -1, id_Heat_PmE = -1, id_intern_heat = -1
  integer :: id_sst = -1, id_sst_sq = -1, id_sss = -1, id_ssu = -1, id_ssv = -1
  integer :: id_speed = -1, id_ssh_inst = -1

  integer :: id_Tadx = -1, id_Tady = -1, id_Tdiffx = -1, id_Tdiffy = -1
  integer :: id_Sadx = -1, id_Sady = -1, id_Sdiffx = -1, id_Sdiffy = -1
  integer :: id_Tadx_2d = -1, id_Tady_2d = -1, id_Tdiffx_2d = -1, id_Tdiffy_2d = -1
  integer :: id_Sadx_2d = -1, id_Sady_2d = -1, id_Sdiffx_2d = -1, id_Sdiffy_2d = -1
  integer :: id_u_predia = -1, id_v_predia = -1, id_h_predia = -1
  integer :: id_T_predia = -1, id_S_predia = -1, id_e_predia = -1
! The remainder of the structure is pointers to child subroutines' control strings.
  type(MOM_dyn_control_struct), pointer :: dyn_CSp => NULL()
 

  type(diabatic_CS), pointer :: diabatic_CSp => NULL()
  type(thickness_diffuse_CS), pointer :: thickness_diffuse_CSp => NULL()
  type(mixedlayer_restrat_CS), pointer :: mixedlayer_restrat_CSp => NULL()
  type(MEKE_CS),  pointer :: MEKE_CSp => NULL()
  type(VarMix_CS),  pointer :: VarMix => NULL()
  type(advect_tracer_CS), pointer :: tracer_CSp => NULL()
  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp => NULL()
  type(diagnostics_CS), pointer :: diagnostics_CSp => NULL()
  type(diag_to_Z_CS), pointer :: diag_to_Z_CSp => NULL()
  type(MOM_restart_CS),  pointer :: restart_CSp => NULL()
  type(ocean_OBC_type), pointer :: OBC => NULL()
  type(ALE_CS), pointer :: ALE_CSp => NULL()
end type MOM_control_struct

type, public :: MOM_dyn_control_struct
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    CAu, &    ! CAu = f*v - u.grad(u) in m s-2.
    PFu, &    ! PFu = -dM/dx, in m s-2.
    diffu, &  ! Zonal acceleration due to convergence of the along-isopycnal
              ! stress tensor, in m s-2.
    visc_rem_u, & ! Both the fraction of the zonal momentum originally in a
              ! layer that remains after a time-step of viscosity, and the
              ! fraction of a time-step's worth of a barotropic acceleration
              ! that a layer experiences after viscosity is applied.
              ! Nondimensional between 0 (at the bottom) and 1 (far above).
    u_accel_bt ! The layers' zonal accelerations due to the difference between
              ! the barotropic accelerations and the baroclinic accelerations
              ! that were fed into the barotopic calculation, in m s-2.
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    CAv, &    ! CAv = -f*u - u.grad(v) in m s-2.
    PFv, &    ! PFv = -dM/dy, in m s-2.
    diffv, &  ! Meridional acceleration due to convergence of the
              ! along-isopycnal stress tensor, in m s-2.
    visc_rem_v, & ! Both the fraction of the meridional momentum originally in
              ! a layer that remains after a time-step of viscosity, and the
              ! fraction of a time-step's worth of a barotropic acceleration
              ! that a layer experiences after viscosity is applied.
              ! Nondimensional between 0 (at the bottom) and 1 (far above).
    v_accel_bt ! The layers' meridional accelerations due to the difference between
              ! the barotropic accelerations and the baroclinic accelerations
              ! that were fed into the barotopic calculation, in m s-2.

! The following variables are only used with the split time stepping scheme.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    eta       ! Instantaneous free surface height, in m.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: u_av
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: v_av
    ! u_av and v_av are the layer velocities with the vertical mean replaced by
    ! the time-mean barotropic velocity over a baroclinic timestep, in m s-1.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_)  :: h_av
    ! The arithmetic mean of two successive layer thicknesses, in m or kg m-2.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    eta_PF    ! The instantaneous SSH used in calculating PFu and PFv, in m.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: uhbt
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: vhbt
    ! uhbt and vhbt are the average volume or mass fluxes determined by the
    ! barotropic solver in m3 s-1 or kg s-1.  uhbt and vhbt should (roughly?) 
    ! equal the verticals sum of uh and vh, respectively.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: uhbt_in
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: vhbt_in
    ! uhbt_in and vhbt_in are the vertically summed transports from based on
    ! the final thicknessses and velocities from the previous dynamics time
    ! step, both in units of m3 s-1 or kg s-1.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: pbce
      ! pbce times eta gives the baroclinic pressure anomaly in each layer due
      ! to free surface height anomalies.  pbce has units of m2 H-1 s-2.

  real, pointer, dimension(:,:) :: taux_bot => NULL(), tauy_bot => NULL()
    ! The frictional bottom stresses from the ocean to the seafloor, in Pa.
  type(BT_cont_type), pointer :: BT_cont => NULL()
                              ! A structure with elements that describe the
                              ! effective summed open face areas as a function
                              ! of barotropic flow.

  ! This is to allow the previous, velocity-based coupling with between the
  ! baroclinic and barotropic modes.
  logical :: flux_BT_coupling  ! If true, use volume fluxes, not velocities,
                               ! to couple the baroclinic and barotropic modes.
  logical :: BT_use_layer_fluxes ! If true, use the summed layered fluxes plus
                               ! an adjustment due to a changed barotropic
                               ! velocity in the barotropic continuity equation.
  logical :: split_bottom_stress  ! If true, provide the bottom stress
                               ! calculated by the vertical viscosity to the
                               ! barotropic solver.
  logical :: readjust_BT_trans ! If true, readjust the barotropic transport of
                               ! the input velocities to agree with CS%uhbt_in
                               ! and CS%vhbt_in after the diabatic step.
  logical :: readjust_velocity ! A flag that varies with time that determines
                               ! whether the velocities currently need to be
                               ! readjusted to agree with CS%uhbt_in and
                               ! CS%vhbt_in.  This is only used if 
                               ! CS%readjust_BT_trans is true.
  logical :: calc_dtbt         ! If true, calculate the barotropic time-step
                               ! dynamically.

  real    :: be              ! A nondimensional number from 0.5 to 1 that controls
                             ! the backward weighting of the time stepping scheme.
  real    :: begw            ! A nondimensional number from 0 to 1 that controls
                             ! the extent to which the treatment of gravity waves
                             ! is forward-backward (0) or simulated backward
                             ! Euler (1).  0 is almost always used.
  logical :: debug           ! If true, write verbose checksums for debugging purposes.

  logical :: module_is_initialized = .false.

  integer :: id_uh = -1, id_vh = -1
  integer :: id_PFu = -1, id_PFv = -1, id_CAu = -1, id_CAv = -1

! Split scheme only.
  integer :: id_uav = -1, id_vav = -1
  integer :: id_du_adj = -1, id_dv_adj = -1, id_du_adj2 = -1, id_dv_adj2 = -1
  integer :: id_u_BT_accel = -1, id_v_BT_accel = -1

  type(diag_ptrs), pointer :: diag ! A structure containing pointers to
                                   ! diagnostic fields that might be calculated
                                   ! and shared between modules.
! The remainder of the structure is pointers to child subroutines' control strings.
  type(hor_visc_CS), pointer :: hor_visc_CSp => NULL()
  type(continuity_CS), pointer :: continuity_CSp => NULL()
  type(CoriolisAdv_CS), pointer :: CoriolisAdv_CSp => NULL()
  type(PressureForce_CS), pointer :: PressureForce_CSp => NULL()
  type(barotropic_CS), pointer :: barotropic_CSp => NULL()
  type(vertvisc_CS), pointer :: vertvisc_CSp => NULL()
  type(set_visc_CS), pointer :: set_visc_CSp => NULL()
  type(open_boundary_CS), pointer :: open_boundary_CSp => NULL()
  type(ocean_OBC_type), pointer :: OBC => NULL() ! A pointer to an open boundary
     ! condition type that specifies whether, where, and  what open boundary
     ! conditions are used.  If no open BCs are used, this pointer stays
     ! nullified.  Flather OBCs use open boundary_CS as well.
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL()

! This is a copy of the pointer in the top-level control structure.
  type(ALE_CS), pointer :: ALE_CSp => NULL()

end type MOM_dyn_control_struct

contains

end module MOM_CS_type
