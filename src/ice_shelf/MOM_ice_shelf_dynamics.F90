!> Implements a crude placeholder for a later implementation of full
!! ice shelf dynamics.
module MOM_ice_shelf_dynamics

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_ROUTINE
use MOM_IS_diag_mediator, only : post_data=>post_IS_data
use MOM_IS_diag_mediator, only : register_diag_field=>register_MOM_IS_diag_field, safe_alloc_ptr
!use MOM_IS_diag_mediator, only : MOM_IS_diag_mediator_init, set_IS_diag_mediator_grid
use MOM_IS_diag_mediator, only : diag_ctrl, time_type, enable_averages, disable_averaging
use MOM_domains, only : MOM_domains_init, clone_MOM_domain
use MOM_domains, only : pass_var, pass_vector, TO_ALL, CGRID_NE, BGRID_NE, CORNER, CENTER
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_grid, only : MOM_grid_init, ocean_grid_type
use MOM_io, only : file_exists, slasher, MOM_read_data
use MOM_io, only : open_ASCII_file, get_filename_appendix
use MOM_io, only : APPEND_FILE, WRITEONLY_FILE
use MOM_restart, only : register_restart_field, MOM_restart_CS
use MOM_time_manager, only : time_type, get_time, set_time, time_type_to_real, operator(>)
use MOM_time_manager,  only : operator(+), operator(-), operator(*), operator(/)
use MOM_time_manager,  only : operator(/=), operator(<=), operator(>=), operator(<)
use MOM_unit_scaling, only : unit_scale_type, unit_scaling_init
!MJH use MOM_ice_shelf_initialize, only : initialize_ice_shelf_boundary
use MOM_ice_shelf_state, only : ice_shelf_state
use MOM_coms, only : reproducing_sum, max_across_PEs, min_across_PEs
use MOM_checksums, only : hchksum, qchksum
use MOM_ice_shelf_initialize, only : initialize_ice_shelf_boundary_channel,initialize_ice_flow_from_file
use MOM_ice_shelf_initialize, only : initialize_ice_shelf_boundary_from_file,initialize_ice_C_basal_friction
use MOM_ice_shelf_initialize, only : initialize_ice_AGlen
implicit none ; private

#include <MOM_memory.h>

public register_ice_shelf_dyn_restarts, initialize_ice_shelf_dyn, update_ice_shelf, IS_dynamics_post_data
public ice_time_step_CFL, ice_shelf_dyn_end, change_in_draft, write_ice_shelf_energy
public shelf_advance_front, ice_shelf_min_thickness_calve, calve_to_mask, volume_above_floatation
public masked_var_grounded

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure for the ice shelf dynamics.
type, public :: ice_shelf_dyn_CS ; private
  real, pointer, dimension(:,:) :: u_shelf => NULL() !< the zonal velocity of the ice shelf/sheet
                                       !! on q-points (B grid) [L T-1 ~> m s-1]
  real, pointer, dimension(:,:) :: v_shelf => NULL() !< the meridional velocity of the ice shelf/sheet
                                       !! on q-points (B grid) [L T-1 ~> m s-1]
  real, pointer, dimension(:,:) :: taudx_shelf => NULL() !< the zonal driving stress of the ice shelf/sheet
                                       !! on q-points (C grid) [R L2 T-2 ~> Pa]
  real, pointer, dimension(:,:) :: taudy_shelf => NULL() !< the meridional driving stress of the ice shelf/sheet
                                       !! on q-points (C grid) [R L2 T-2 ~> Pa]
  real, pointer, dimension(:,:) :: sx_shelf => NULL() !< the zonal surface slope of the ice shelf/sheet
                                       !! on q-points (B grid) [nondim]
  real, pointer, dimension(:,:) :: sy_shelf => NULL() !< the meridional surface slope of the ice shelf/sheet
                                       !! on q-points (B grid) [nondim]
  real, pointer, dimension(:,:) :: u_face_mask => NULL() !< mask for velocity boundary conditions on the C-grid
                                       !! u-face - this is because the FEM cares about FACES THAT GET INTEGRATED OVER,
                                       !! not vertices. Will represent boundary conditions on computational boundary
                                       !! (or permanent boundary between fast-moving and near-stagnant ice
                                       !! FOR NOW: 1=interior bdry, 0=no-flow boundary, 2=stress bdry condition,
                                       !! 3=inhomogeneous Dirichlet boundary for u and v, 4=flux boundary: at these
                                       !! faces a flux will be specified which will override velocities; a homogeneous
                                       !! velocity condition will be specified (this seems to give the solver less
                                       !! difficulty)  5=inhomogenous Dirichlet boundary for u only. 6=inhomogenous
                                       !! Dirichlet boundary for v only
  real, pointer, dimension(:,:) :: v_face_mask => NULL()  !< A mask for velocity boundary conditions on the C-grid
                                       !! v-face, with valued defined similarly to u_face_mask, but 5 is Dirichlet for v
                                       !! and 6 is Dirichlet for u
  real, pointer, dimension(:,:) :: u_face_mask_bdry => NULL() !< A duplicate copy of u_face_mask?
  real, pointer, dimension(:,:) :: v_face_mask_bdry => NULL() !< A duplicate copy of v_face_mask?
  real, pointer, dimension(:,:) :: u_flux_bdry_val => NULL() !< The ice volume flux per unit face length into the cell
                                       !! through open boundary u-faces (where u_face_mask=4) [Z L T-1 ~> m2 s-1]
  real, pointer, dimension(:,:) :: v_flux_bdry_val => NULL() !< The ice volume flux per unit face length into the cell
                                       !! through open boundary v-faces (where v_face_mask=4) [Z L T-1 ~> m2 s-1]??
   ! needed where u_face_mask is equal to 4, similarly for v_face_mask
  real, pointer, dimension(:,:) :: umask => NULL()      !< u-mask on the actual degrees of freedom (B grid)
                                       !! 1=normal node, 3=inhomogeneous boundary node,
                                       !!  0 - no flow node (will also get ice-free nodes)
  real, pointer, dimension(:,:) :: vmask => NULL()      !< v-mask on the actual degrees of freedom (B grid)
                                       !! 1=normal node, 3=inhomogeneous boundary node,
                                       !!  0 - no flow node (will also get ice-free nodes)
  real, pointer, dimension(:,:) :: calve_mask => NULL() !< a mask to prevent the ice shelf front from
                                          !! advancing past its initial position (but it may retreat)
  real, pointer, dimension(:,:) :: t_shelf => NULL() !< Vertically integrated temperature in the ice shelf/stream,
                                                     !! on corner-points (B grid) [C ~> degC]
  real, pointer, dimension(:,:) :: tmask => NULL()   !< A mask on tracer points that is 1 where there is ice.
  real, pointer, dimension(:,:,:) :: ice_visc => NULL() !< Area and depth-integrated Glen's law ice viscosity
                                                        !!  (Pa m3 s) in [R L4 Z T-1 ~> kg m2 s-1].
                                                        !!  at either 1 (cell-centered) or 4 quadrature points per cell
  real, pointer, dimension(:,:) :: AGlen_visc => NULL() !< Ice-stiffness parameter in Glen's law ice viscosity,
                                                      !! often in [Pa-3 s-1] if n_Glen is 3.
  real, pointer, dimension(:,:) :: u_bdry_val => NULL() !< The zonal ice velocity at inflowing boundaries
                                       !! [L yr-1 ~> m yr-1]
  real, pointer, dimension(:,:) :: v_bdry_val => NULL() !< The meridional ice velocity at inflowing boundaries
                                       !! [L yr-1 ~> m yr-1]
  real, pointer, dimension(:,:) :: h_bdry_val => NULL() !< The ice thickness at inflowing boundaries [Z ~> m].
  real, pointer, dimension(:,:) :: t_bdry_val => NULL() !< The ice temperature at inflowing boundaries [C ~> degC].

  real, pointer, dimension(:,:) :: bed_elev => NULL()  !< The bed elevation used for ice dynamics [Z ~> m],
                                                       !! relative to mean sea-level.  This is
                                                       !! the same as G%bathyT+Z_ref, when below sea-level.
                                                       !! Sign convention: positive below sea-level, negative above.

  real, pointer, dimension(:,:) :: basal_traction => NULL() !< The area-integrated taub_beta field
                                                            !! (m2 Pa s m-1, or kg s-1) related to the nonlinear part
                                                            !! of "linearized" basal stress (Pa) [R L3 T-1 ~> kg s-1]
                !!  The exact form depends on basal law exponent and/or whether flow is "hybridized" a la Goldberg 2011
  real, pointer, dimension(:,:) :: C_basal_friction => NULL()!< Coefficient in sliding law tau_b = C u^(n_basal_fric),
                                                            !!  units= Pa (s m-1)^(n_basal_fric)
  real, pointer, dimension(:,:) :: OD_rt => NULL()         !< A running total for calculating OD_av [Z ~> m].
  real, pointer, dimension(:,:) :: ground_frac_rt => NULL() !< A running total for calculating ground_frac.
  real, pointer, dimension(:,:) :: OD_av => NULL()         !< The time average open ocean depth [Z ~> m].
  real, pointer, dimension(:,:) :: ground_frac => NULL()   !< Fraction of the time a cell is "exposed", i.e. the column
                               !! thickness is below a threshold and interacting with the rock [nondim].  When this
                               !! is 1, the ice-shelf is grounded
  real, pointer, dimension(:,:) :: float_cond => NULL()   !< If GL_regularize=true, indicates cells containing
                                                !! the grounding line (float_cond=1) or not (float_cond=0)
  real, pointer, dimension(:,:,:,:) :: Phi => NULL() !< The gradients of bilinear basis elements at Gaussian
                                                !! 4 quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real, pointer, dimension(:,:,:) :: PhiC => NULL()  !< The gradients of bilinear basis elements at 1 cell-centered
                                                !! quadrature point per cell [L-1 ~> m-1].
  real, pointer, dimension(:,:,:,:,:,:) :: Phisub => NULL() !< Quadrature structure weights at subgridscale
                                                !!  locations for finite element calculations [nondim]
  integer :: OD_rt_counter = 0 !< A counter of the number of contributions to OD_rt.

  real :: velocity_update_time_step !< The time interval over which to update the ice shelf velocity
                    !! using the nonlinear elliptic equation, or 0 to update every timestep [T ~> s].
                    ! DNGoldberg thinks this should be done no more often than about once a day
                    ! (maybe longer) because it will depend on ocean values  that are averaged over
                    ! this time interval, and solving for the equilibrated flow will begin to lose
                    ! meaning if it is done too frequently.
  real :: elapsed_velocity_time  !< The elapsed time since the ice velocities were last updated [T ~> s].

  real :: g_Earth      !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real :: density_ice  !< A typical density of ice [R ~> kg m-3].
  real :: Cp_ice       !< The heat capacity of fresh ice [Q C-1 ~> J kg-1 degC-1].

  logical :: advect_shelf !< If true (default), advect ice shelf and evolve thickness
  logical :: reentrant_x !< If true, the domain is zonally reentrant
  logical :: reentrant_y !< If true, the domain is meridionally reentrant
  logical :: alternate_first_direction_IS !< If true, alternate whether the x- or y-direction
                                          !! updates occur first in directionally split parts of the calculation.
  integer :: first_direction_IS !< An integer that indicates which direction is
                                !! to be updated first in directionally split
                                !! parts of the ice sheet calculation (e.g. advection).
  real    :: first_dir_restart_IS = -1.0 !< A real copy of CS%first_direction_IS for use in restart files
  integer :: visc_qps !< The number of quadrature points per cell (1 or 4) on which to calculate ice viscosity.
  character(len=40) :: ice_viscosity_compute !< Specifies whether the ice viscosity is computed internally
                                   !! according to Glen's flow law; is constant (for debugging purposes)
                                   !! or using observed strain rates and read from a file
  logical :: GL_regularize  !< Specifies whether to regularize the floatation condition
                            !! at the grounding line as in Goldberg Holland Schoof 2009
  integer :: n_sub_regularize
                            !< partition of cell over which to integrate for
                            !! interpolated grounding line the (rectangular) is
                            !! divided into nxn equally-sized rectangles, over which
                            !!  basal contribution is integrated (iterative quadrature)
  logical :: GL_couple      !< whether to let the floatation condition be
                            !! determined by ocean column thickness means update_OD_ffrac
                            !! will be called (note: GL_regularize and GL_couple
                            !! should be exclusive)

  real    :: CFL_factor     !< A factor used to limit subcycled advective timestep in uncoupled runs
                            !! i.e. dt <= CFL_factor * min(dx / u) [nondim]

  real :: min_h_shelf !< The minimum ice thickness used during ice dynamics [L ~> m].
  real :: min_basal_traction !< The minimum basal traction for grounded ice (Pa m-1 s) [R L T-1 ~> kg m-2 s-1]
  real :: max_surface_slope !< The maximum allowed ice-sheet surface slope (to ignore, set to zero) [nondim]
  real :: min_ice_visc !< The minimum allowed Glen's law ice viscosity (Pa s), in [R L2 T-1 ~> kg m-1 s-1].

  real :: n_glen            !< Nonlinearity exponent in Glen's Law [nondim]
  real :: eps_glen_min      !< Min. strain rate to avoid infinite Glen's law viscosity, [T-1 ~> s-1].
  real :: n_basal_fric      !< Exponent in sliding law tau_b = C u^(m_slide) [nondim]
  logical :: CoulombFriction !< Use Coulomb friction law (Schoof 2005, Gagliardini et al 2007)
  real :: CF_MinN           !< Minimum Coulomb friction effective pressure [Pa]
  real :: CF_PostPeak       !< Coulomb friction post peak exponent [nondim]
  real :: CF_Max            !< Coulomb friction maximum coefficient [nondim]
  real :: density_ocean_avg !< A typical ocean density [R ~> kg m-3].  This does not affect ocean
                            !! circulation or thermodynamics.  It is used to estimate the
                            !! gravitational driving force at the shelf front (until we think of
                            !! a better way to do it, but any difference will be negligible).
  real :: thresh_float_col_depth !< The water column depth over which the shelf if considered to be floating
  logical :: moving_shelf_front  !< Specify whether to advance shelf front (and calve).
  logical :: calve_to_mask       !< If true, calve off the ice shelf when it passes the edge of a mask.
  real :: min_thickness_simple_calve !< min. ice shelf thickness criteria for calving [Z ~> m].
  real :: T_shelf_missing   !< An ice shelf temperature to use where there is no ice shelf [C ~> degC]
  real :: cg_tolerance !< The tolerance in the CG solver, relative to initial residual, that
                       !! determines when to stop the conjugate gradient iterations [nondim].
  real :: nonlinear_tolerance !< The fractional nonlinear tolerance, relative to the initial error,
                              !! that sets when to stop the iterative velocity solver [nondim]
  integer :: cg_max_iterations !< The maximum number of iterations that can be used in the CG solver
  integer :: nonlin_solve_err_mode  !< 1: exit vel solve based on nonlin residual
                    !! 2: exit based on "fixed point" metric (|u - u_last| / |u| < tol) where | | is infty-norm
                    !! 3: exit based on change of norm

  ! for write_ice_shelf_energy
  type(time_type) :: energysavedays            !< The interval between writing the energies
                                               !! and other integral quantities of the run.
  type(time_type) :: energysavedays_geometric  !< The starting interval for computing a geometric
                                               !! progression of time deltas between calls to
                                               !! write_energy. This interval will increase by a factor of 2.
                                               !! after each call to write_energy.
  logical         :: energysave_geometric      !< Logical to control whether calls to write_energy should
                                               !! follow a geometric progression
  type(time_type) :: write_energy_time         !< The next time to write to the energy file.
  type(time_type) :: geometric_end_time        !< Time at which to stop the geometric progression
                                               !! of calls to write_energy and revert to the standard
                                               !! energysavedays interval
  real    :: timeunit           !< The length of the units for the time axis and certain input parameters
                                !! including ENERGYSAVEDAYS [s].
  type(time_type) :: Start_time !< The start time of the simulation.
                                ! Start_time is set in MOM_initialization.F90
  integer :: prev_IS_energy_calls = 0 !< The number of times write_ice_shelf_energy has been called.
  integer :: IS_fileenergy_ascii   !< The unit number of the ascii version of the energy file.
  character(len=200) :: IS_energyfile  !< The name of the ice sheet energy file with path.

  ! ids for outputting intermediate thickness in advection subroutine (debugging)
  !integer :: id_h_after_uflux = -1, id_h_after_vflux = -1, id_h_after_adv = -1

  logical :: debug                !< If true, write verbose checksums for debugging purposes
                                  !! and use reproducible sums
  logical :: module_is_initialized = .false. !< True if this module has been initialized.

  !>@{ Diagnostic handles
  integer :: id_u_shelf = -1, id_v_shelf = -1, id_shelf_speed, id_t_shelf = -1, &
             id_taudx_shelf = -1, id_taudy_shelf = -1, id_taud_shelf = -1, id_bed_elev = -1, &
             id_ground_frac = -1, id_col_thick = -1, id_OD_av = -1, id_float_cond = -1, &
             id_u_mask = -1, id_v_mask = -1, id_ufb_mask =-1, id_vfb_mask = -1, id_t_mask = -1, &
             id_sx_shelf = -1, id_sy_shelf = -1, id_surf_slope_mag_shelf, &
             id_duHdx = -1, id_dvHdy = -1, id_fluxdiv = -1, &
             id_strainrate_xx = -1, id_strainrate_yy = -1, id_strainrate_xy = -1, &
             id_pstrainrate_1 = -1, id_pstrainrate_2, &
             id_devstress_xx = -1, id_devstress_yy = -1, id_devstress_xy = -1, &
             id_pdevstress_1 = -1, id_pdevstress_2 = -1

  !>@}
  ! ids for outputting intermediate thickness in advection subroutine (debugging)
  !>@{ Diagnostic handles for debugging
  integer :: id_h_after_uflux = -1, id_h_after_vflux = -1, id_h_after_adv = -1, &
             id_visc_shelf = -1, id_taub = -1
  !>@}
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to control diagnostic output.

end type ice_shelf_dyn_CS

!> A container for loop bounds
type :: loop_bounds_type ; private
  !>@{ Loop bounds
  integer :: ish, ieh, jsh, jeh
  !>@}
end type loop_bounds_type

contains

!> used for flux limiting in advective subroutines Van Leer limiter (source: Wikipedia)
!! The return value is between 0 and 2 [nondim].
function slope_limiter(num, denom)
  real, intent(in)    :: num   !< The numerator of the ratio used in the Van Leer slope limiter
  real, intent(in)    :: denom !< The denominator of the ratio used in the Van Leer slope limiter
  real :: slope_limiter ! The slope limiter value, between 0 and 2 [nondim].
  real :: r  ! The ratio of num/denom [nondim]

  if (denom == 0) then
    slope_limiter = 0
  elseif (num*denom <= 0) then
    slope_limiter = 0
  else
    r = num/denom
    slope_limiter = (r+abs(r))/(1+abs(r))
  endif

end function slope_limiter

!> Calculate area of quadrilateral.
function quad_area (X, Y)
  real, dimension(4), intent(in) :: X !< The x-positions of the vertices of the quadrilateral [L ~> m].
  real, dimension(4), intent(in) :: Y !< The y-positions of the vertices of the quadrilateral [L ~> m].
  real :: quad_area ! Computed area [L2 ~> m2]
  real :: p2, q2, a2, c2, b2, d2

! X and Y must be passed in the form
    !  3 - 4
    !  |   |
    !  1 - 2

  p2 = ( ((X(4)-X(1))**2) + ((Y(4)-Y(1))**2) ) ; q2 = ( ((X(3)-X(2))**2) + ((Y(3)-Y(2))**2) )
  a2 = ( ((X(3)-X(4))**2) + ((Y(3)-Y(4))**2) ) ; c2 = ( ((X(1)-X(2))**2) + ((Y(1)-Y(2))**2) )
  b2 = ( ((X(2)-X(4))**2) + ((Y(2)-Y(4))**2) ) ; d2 = ( ((X(3)-X(1))**2) + ((Y(3)-Y(1))**2) )
  quad_area = .25 * sqrt(4*P2*Q2-(B2+D2-A2-C2)**2)

end function quad_area

!> This subroutine is used to register any fields related to the ice shelf
!! dynamics that should be written to or read from the restart file.
subroutine register_ice_shelf_dyn_restarts(G, US, param_file, CS, restart_CS)
  type(ocean_grid_type),  intent(inout) :: G    !< The grid type describing the ice shelf grid.
  type(unit_scale_type),  intent(in)    :: US   !< A structure containing unit conversion factors
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(ice_shelf_dyn_CS), pointer       :: CS !< A pointer to the ice shelf dynamics control structure
  type(MOM_restart_CS),   intent(inout) :: restart_CS !< MOM restart control struct

  ! Local variables
  real :: T_shelf_missing ! An ice shelf temperature to use where there is no ice shelf [C ~> degC]
  logical :: shelf_mass_is_dynamic, override_shelf_movement, active_shelf_dynamics
  character(len=40)  :: mdl = "MOM_ice_shelf_dyn"  ! This module's name.
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(FATAL, "MOM_ice_shelf_dyn.F90, register_ice_shelf_dyn_restarts: "// &
                          "called with an associated control structure.")
    return
  endif
  allocate(CS)

  override_shelf_movement = .false. ; active_shelf_dynamics = .false.
  call get_param(param_file, mdl, "DYNAMIC_SHELF_MASS", shelf_mass_is_dynamic, &
                 "If true, the ice sheet mass can evolve with time.", &
                 default=.false., do_not_log=.true.)
  if (shelf_mass_is_dynamic) then
    call get_param(param_file, mdl, "OVERRIDE_SHELF_MOVEMENT", override_shelf_movement, &
                 "If true, user provided code specifies the ice-shelf "//&
                 "movement instead of the dynamic ice model.", default=.false., do_not_log=.true.)
    active_shelf_dynamics = .not.override_shelf_movement
  endif

  if (active_shelf_dynamics) then
    call get_param(param_file, mdl, "MISSING_SHELF_TEMPERATURE", T_shelf_missing, &
                 "An ice shelf temperature to use where there is no ice shelf.",&
                 units="degC", default=-10.0, scale=US%degC_to_C, do_not_log=.true.)

    call get_param(param_file, mdl, "NUMBER_OF_ICE_VISCOSITY_QUADRATURE_POINTS", CS%visc_qps, &
                 "Number of ice viscosity quadrature points. Either 1 (cell-centered) for 4", &
                  units="none", default=1)
    if (CS%visc_qps/=1 .and. CS%visc_qps/=4) call MOM_error (FATAL, &
      "NUMBER OF ICE_VISCOSITY_QUADRATURE_POINTS must be 1 or 4")

    call get_param(param_file, mdl, "FIRST_DIRECTION_IS", CS%first_direction_IS, &
                 "An integer that indicates which direction goes first "//&
                 "in parts of the code that use directionally split "//&
                 "updates (e.g. advection), with even numbers (or 0) used for x- first "//&
                 "and odd numbers used for y-first.", default=0)
    call get_param(param_file, mdl, "ALTERNATE_FIRST_DIRECTION_IS", CS%alternate_first_direction_IS, &
                 "If true, after every advection call, alternate whether the x- or y- "//&
                 "direction advection updates occur first. "//&
                 "If this is true, FIRST_DIRECTION applies at the start of a new run or if "//&
                 "the next first direction can not be found in the restart file.", default=.false.)

    allocate(CS%u_shelf(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%v_shelf(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%t_shelf(isd:ied,jsd:jed), source=T_shelf_missing) ! [C ~> degC]
    allocate(CS%ice_visc(isd:ied,jsd:jed,CS%visc_qps), source=0.0)
    allocate(CS%AGlen_visc(isd:ied,jsd:jed), source=2.261e-25) ! [Pa-3 s-1]
    allocate(CS%basal_traction(isd:ied,jsd:jed), source=0.0)   ! [R L3 T-1 ~> kg s-1]
    allocate(CS%C_basal_friction(isd:ied,jsd:jed), source=5.0e10) ! [Pa (s m-1)^n_sliding]
    allocate(CS%OD_av(isd:ied,jsd:jed), source=0.0)
    allocate(CS%ground_frac(isd:ied,jsd:jed), source=0.0)
    allocate(CS%taudx_shelf(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%taudy_shelf(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%sx_shelf(isd:ied,jsd:jed), source=0.0)
    allocate(CS%sy_shelf(isd:ied,jsd:jed), source=0.0)
    allocate(CS%bed_elev(isd:ied,jsd:jed), source=0.0)
    allocate(CS%u_bdry_val(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%v_bdry_val(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%u_face_mask_bdry(IsdB:IedB,JsdB:JedB), source=-2.0)
    allocate(CS%v_face_mask_bdry(IsdB:iedB,JsdB:JedB), source=-2.0)
    allocate(CS%h_bdry_val(isd:ied,jsd:jed), source=0.0)

   ! additional restarts for ice shelf state
    call register_restart_field(CS%u_shelf, "u_shelf", .false., restart_CS, &
                                "ice sheet/shelf u-velocity", &
                                units="m s-1", conversion=US%L_T_to_m_s, hor_grid='Bu')
    call register_restart_field(CS%v_shelf, "v_shelf", .false., restart_CS, &
                                "ice sheet/shelf v-velocity", &
                                units="m s-1", conversion=US%L_T_to_m_s, hor_grid='Bu')
    call register_restart_field(CS%u_bdry_val, "u_bdry_val", .false., restart_CS, &
                                "ice sheet/shelf boundary u-velocity", &
                                units="m s-1", conversion=US%L_T_to_m_s, hor_grid='Bu')
    call register_restart_field(CS%v_bdry_val, "v_bdry_val", .false., restart_CS, &
                                "ice sheet/shelf boundary v-velocity", &
                                units="m s-1", conversion=US%L_T_to_m_s, hor_grid='Bu')
    call register_restart_field(CS%u_face_mask_bdry, "u_face_mask_bdry", .false., restart_CS, &
                                "ice sheet/shelf boundary u-mask", "nondim", hor_grid='Bu')
    call register_restart_field(CS%v_face_mask_bdry, "v_face_mask_bdry", .false., restart_CS, &
                                "ice sheet/shelf boundary v-mask", "nondim", hor_grid='Bu')

    call register_restart_field(CS%OD_av, "OD_av", .true., restart_CS, &
                                "Average open ocean depth in a cell", "m", conversion=US%Z_to_m)
    call register_restart_field(CS%ground_frac, "ground_frac", .true., restart_CS, &
                                "fractional degree of grounding", "nondim")
    call register_restart_field(CS%C_basal_friction, "C_basal_friction", .true., restart_CS, &
                                "basal sliding coefficients", "Pa (s m-1)^n_sliding")
    call register_restart_field(CS%AGlen_visc, "AGlen_visc", .true., restart_CS, &
                                "ice-stiffness parameter", "Pa-3 s-1")
    call register_restart_field(CS%h_bdry_val, "h_bdry_val", .false., restart_CS, &
                                "ice thickness at the boundary", "m", conversion=US%Z_to_m)
    call register_restart_field(CS%bed_elev, "bed elevation", .true., restart_CS, &
                                "bed elevation", "m", conversion=US%Z_to_m)
    call register_restart_field(CS%first_dir_restart_IS, "first_direction_IS", .false., restart_CS, &
                                "Indicator of the first direction in split ice shelf calculations.", "nondim")
  endif

end subroutine register_ice_shelf_dyn_restarts

!> Initializes shelf model data, parameters and diagnostics
subroutine initialize_ice_shelf_dyn(param_file, Time, ISS, CS, G, US, diag, new_sim, Cp_ice, &
                                    Input_start_time, directory, solo_ice_sheet_in)
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(time_type),         intent(inout) :: Time !< The clock that that will indicate the model time
  type(ice_shelf_state),   intent(in)    :: ISS  !< A structure with elements that describe
                                                 !! the ice-shelf state
  type(ice_shelf_dyn_CS),  pointer       :: CS   !< A pointer to the ice shelf dynamics control structure
  type(ocean_grid_type),   intent(inout) :: G    !< The grid type describing the ice shelf grid.
  type(unit_scale_type),   intent(in)    :: US   !< A structure containing unit conversion factors
  type(diag_ctrl), target, intent(in)    :: diag !< A structure that is used to regulate the diagnostic output.
  logical,                 intent(in)    :: new_sim !< If true this is a new simulation, otherwise
                                                 !! has been started from a restart file.
  real,                    intent(in)    :: Cp_ice !< Heat capacity of ice [Q C-1 ~> J kg-1 degC-1]
  type(time_type),         intent(in)    :: Input_start_time !< The start time of the simulation.
  character(len=*),        intent(in)    :: directory  !< The directory where the ice sheet energy file goes.
  logical,       optional, intent(in)    :: solo_ice_sheet_in !< If present, this indicates whether
                                                 !! a solo ice-sheet driver.

  ! Local variables
  real    :: T_shelf_bdry ! A default ice shelf temperature to use for ice flowing
                          ! in through open boundaries [C ~> degC]
  !This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=200) :: IC_file,filename,inputdir
  character(len=40)  :: var_name
  character(len=40)  :: mdl = "MOM_ice_shelf_dyn"  ! This module's name.
  logical :: shelf_mass_is_dynamic, override_shelf_movement, active_shelf_dynamics
  logical :: debug
  integer :: i, j, isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq, iters
  character(len=200) :: IS_energyfile  ! The name of the energy file.
  character(len=32) :: filename_appendix = '' ! FMS appendix to filename for ensemble runs

  Isdq = G%isdB ; Iedq = G%iedB ; Jsdq = G%jsdB ; Jedq = G%jedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) then
    call MOM_error(FATAL, "MOM_ice_shelf_dyn.F90, initialize_ice_shelf_dyn: "// &
                          "called with an associated control structure.")
    return
  endif
  if (CS%module_is_initialized) then
    call MOM_error(WARNING, "MOM_ice_shelf_dyn.F90, initialize_ice_shelf_dyn was "//&
             "called with a control structure that has already been initialized.")
  endif
  CS%module_is_initialized = .true.

  CS%diag => diag ! ; CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DEBUG", debug, default=.false.)
  call get_param(param_file, mdl, "DEBUG_IS", CS%debug, &
                 "If true, write verbose debugging messages for the ice shelf.", &
                 default=debug)
  call get_param(param_file, mdl, "DYNAMIC_SHELF_MASS", shelf_mass_is_dynamic, &
                 "If true, the ice sheet mass can evolve with time.", &
                 default=.false.)
  override_shelf_movement = .false. ; active_shelf_dynamics = .false.
  if (shelf_mass_is_dynamic) then
    call get_param(param_file, mdl, "OVERRIDE_SHELF_MOVEMENT", override_shelf_movement, &
                 "If true, user provided code specifies the ice-shelf "//&
                 "movement instead of the dynamic ice model.", default=.false., do_not_log=.true.)
    active_shelf_dynamics = .not.override_shelf_movement

    call get_param(param_file, mdl, "GROUNDING_LINE_INTERPOLATE", CS%GL_regularize, &
                 "If true, regularize the floatation condition at the "//&
                 "grounding line as in Goldberg Holland Schoof 2009.", default=.false.)
    call get_param(param_file, mdl, "GROUNDING_LINE_INTERP_SUBGRID_N", CS%n_sub_regularize, &
                 "The number of sub-partitions of each cell over which to "//&
                 "integrate for the interpolated grounding line. Each cell "//&
                 "is divided into NxN equally-sized rectangles, over which the "//&
                 "basal contribution is integrated by iterative quadrature.", &
                 default=0)
    call get_param(param_file, mdl, "GROUNDING_LINE_COUPLE", CS%GL_couple, &
                 "If true, let the floatation condition be determined by "//&
                 "ocean column thickness. This means that update_OD_ffrac "//&
                 "will be called.  GL_REGULARIZE and GL_COUPLE are exclusive.", &
                 default=.false., do_not_log=CS%GL_regularize)
    if (CS%GL_regularize) CS%GL_couple = .false.
    if (present(solo_ice_sheet_in)) then
      if (solo_ice_sheet_in) CS%GL_couple = .false.
    endif
    if (CS%GL_regularize .and. (CS%n_sub_regularize == 0)) call MOM_error (FATAL, &
      "GROUNDING_LINE_INTERP_SUBGRID_N must be a positive integer if GL regularization is used")
    call get_param(param_file, mdl, "ICE_SHELF_CFL_FACTOR", CS%CFL_factor, &
                 "A factor used to limit timestep as CFL_FACTOR * min (\Delta x / u). "//&
                 "This is only used with an ice-only model.", units="nondim", default=0.25)
  endif
  call get_param(param_file, mdl, "RHO_0", CS%density_ocean_avg, &
                 "avg ocean density used in floatation cond", &
                 units="kg m-3", default=1035., scale=US%kg_m3_to_R)
  if (active_shelf_dynamics) then
    call get_param(param_file, mdl, "ICE_VELOCITY_TIMESTEP", CS%velocity_update_time_step, &
                 "seconds between ice velocity calcs", units="s", scale=US%s_to_T, &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, "G_EARTH", CS%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default=9.80, scale=US%m_s_to_L_T**2*US%Z_to_m)

    call get_param(param_file, mdl, "MIN_H_SHELF", CS%min_h_shelf, &
                 "min. ice thickness used during ice dynamics", &
                  units="m", default=0.,scale=US%m_to_L)
    call get_param(param_file, mdl, "MIN_BASAL_TRACTION", CS%min_basal_traction, &
                 "min. allowed basal traction. Input is in [Pa m-1 yr], but is converted when read in to [Pa m-1 s]", &
                 units="Pa m-1 yr", default=0., scale=365.0*86400.0*US%Pa_to_RLZ_T2*US%L_T_to_m_s)
    call get_param(param_file, mdl, "MAX_SURFACE_SLOPE", CS%max_surface_slope, &
                 "max. allowed ice-sheet surface slope. To ignore, set to zero.", &
                 units="none", default=0., scale=US%m_to_Z/US%m_to_L)
    call get_param(param_file, mdl, "MIN_ICE_VISC", CS%min_ice_visc, &
                 "min. allowed Glen's law ice viscosity", &
                 units="Pa s", default=0., scale=US%Pa_to_RL2_T2*US%s_to_T)

    call get_param(param_file, mdl, "GLEN_EXPONENT", CS%n_glen, &
                 "nonlinearity exponent in Glen's Law", &
                  units="none", default=3.)
    call get_param(param_file, mdl, "MIN_STRAIN_RATE_GLEN", CS%eps_glen_min, &
                 "min. strain rate to avoid infinite Glen's law viscosity", &
                 units="s-1", default=1.e-19, scale=US%T_to_s)
    call get_param(param_file, mdl, "BASAL_FRICTION_EXP", CS%n_basal_fric, &
                 "Exponent in sliding law \tau_b = C u^(n_basal_fric)", &
                 units="none", fail_if_missing=.true.)
    call get_param(param_file, mdl, "USE_COULOMB_FRICTION", CS%CoulombFriction, &
                 "Use Coulomb Friction Law", &
                 units="none", default=.false., fail_if_missing=.false.)
    call get_param(param_file, mdl, "CF_MinN", CS%CF_MinN, &
                 "Minimum Coulomb friction effective pressure", &
                 units="Pa", default=1.0, fail_if_missing=.false.)
    call get_param(param_file, mdl, "CF_PostPeak", CS%CF_PostPeak, &
                 "Coulomb friction post peak exponent", &
                 units="none", default=1.0, fail_if_missing=.false.)
    call get_param(param_file, mdl, "CF_Max", CS%CF_Max, &
                 "Coulomb friction maximum coefficient", &
                 units="none", default=0.5, fail_if_missing=.false.)

    call get_param(param_file, mdl, "DENSITY_ICE", CS%density_ice, &
                 "A typical density of ice.", units="kg m-3", default=917.0, scale=US%kg_m3_to_R)
    call get_param(param_file, mdl, "CONJUGATE_GRADIENT_TOLERANCE", CS%cg_tolerance, &
                "tolerance in CG solver, relative to initial residual", units="nondim", default=1.e-6)
    call get_param(param_file, mdl, "ICE_NONLINEAR_TOLERANCE", CS%nonlinear_tolerance, &
                "nonlin tolerance in iterative velocity solve", units="nondim", default=1.e-6)
    call get_param(param_file, mdl, "CONJUGATE_GRADIENT_MAXIT", CS%cg_max_iterations, &
                "max iteratiions in CG solver", default=2000)
    call get_param(param_file, mdl, "THRESH_FLOAT_COL_DEPTH", CS%thresh_float_col_depth, &
                "min ocean thickness to consider ice *floating*; "//&
                "will only be important with use of tides", &
                units="m", default=1.e-3, scale=US%m_to_Z)
    call get_param(param_file, mdl, "NONLIN_SOLVE_ERR_MODE", CS%nonlin_solve_err_mode, &
                "Choose whether nonlin error in vel solve is based on nonlinear "//&
                "residual (1), relative change since last iteration (2), or change in norm (3)", default=3)

    call get_param(param_file, mdl, "SHELF_MOVING_FRONT", CS%moving_shelf_front, &
                 "Specify whether to advance shelf front (and calve).", &
                 default=.false.)
    call get_param(param_file, mdl, "CALVE_TO_MASK", CS%calve_to_mask, &
                 "If true, do not allow an ice shelf where prohibited by a mask.", &
                 default=.false.)
    call get_param(param_file, mdl, "ADVECT_SHELF", CS%advect_shelf, &
                 "If true, advect ice shelf and evolve thickness", &
                 default=.true.)
    call get_param(param_file, mdl, "REENTRANT_X", CS%reentrant_x, &
                 " If true, the domain is zonally reentrant.", &
                 default=.false.)
    call get_param(param_file, mdl, "REENTRANT_Y", CS%reentrant_y, &
                 " If true, the domain is meridionally reentrant.", &
                 default=.false.)
    call get_param(param_file, mdl, "ICE_VISCOSITY_COMPUTE", CS%ice_viscosity_compute, &
                 "If MODEL, compute ice viscosity internally using 1 or 4 quadrature points,"//&
                 "if OBS read from a file,"//&
                 "if CONSTANT a constant value (for debugging).", &
                 default="MODEL")

    if ((CS%visc_qps/=1) .and. (trim(CS%ice_viscosity_compute) /= "MODEL")) then
      call MOM_error(FATAL, "NUMBER_OF_ICE_VISCOSITY_QUADRATURE_POINTS must be 1 unless ICE_VISCOSITY_COMPUTE==MODEL.")
    endif
    call get_param(param_file, mdl, "INFLOW_SHELF_TEMPERATURE", T_shelf_bdry, &
                 "A default ice shelf temperature to use for ice flowing in through "//&
                 "open boundaries.", units="degC", default=-15.0, scale=US%degC_to_C)
  endif
  call get_param(param_file, mdl, "MISSING_SHELF_TEMPERATURE", CS%T_shelf_missing, &
                 "An ice shelf temperature to use where there is no ice shelf.",&
                 units="degC", default=-10.0, scale=US%degC_to_C)
  call get_param(param_file, mdl, "MIN_THICKNESS_SIMPLE_CALVE", CS%min_thickness_simple_calve, &
                 "Min thickness rule for the VERY simple calving law",&
                 units="m", default=0.0, scale=US%m_to_Z)
  CS%Cp_ice = Cp_ice !Heat capacity of ice (J kg-1 K-1), needed for heat flux of any bergs calved from
                     !the ice shelf and for ice sheet temperature solver
  !for write_ice_shelf_energy
      ! Note that the units of CS%Timeunit are the MKS units of [s].
  call get_param(param_file, mdl, "TIMEUNIT", CS%Timeunit, &
    "The time unit in seconds a number of input fields", &
    units="s", default=86400.0)
  if (CS%Timeunit < 0.0) CS%Timeunit = 86400.0
  call get_param(param_file, mdl, "ENERGYSAVEDAYS",CS%energysavedays, &
    "The interval in units of TIMEUNIT between saves of the "//&
    "energies of the run and other globally summed diagnostics.",&
    default=set_time(0,days=1), timeunit=CS%Timeunit)
  call get_param(param_file, mdl, "ENERGYSAVEDAYS_GEOMETRIC",CS%energysavedays_geometric, &
    "The starting interval in units of TIMEUNIT for the first call "//&
    "to save the energies of the run and other globally summed diagnostics. "//&
    "The interval increases by a factor of 2. after each call to write_ice_shelf_energy.",&
    default=set_time(seconds=0), timeunit=CS%Timeunit)
  if ((time_type_to_real(CS%energysavedays_geometric) > 0.) .and. &
    (CS%energysavedays_geometric < CS%energysavedays)) then
    CS%energysave_geometric = .true.
  else
    CS%energysave_geometric = .false.
  endif
  CS%Start_time = Input_start_time
  call get_param(param_file, mdl, "ICE_SHELF_ENERGYFILE", IS_energyfile, &
                 "The file to use to write the energies and globally "//&
                 "summed diagnostics.", default="ice_shelf.stats")
  !query fms_io if there is a filename_appendix (for ensemble runs)
  call get_filename_appendix(filename_appendix)
  if (len_trim(filename_appendix) > 0) then
    IS_energyfile = trim(IS_energyfile) //'.'//trim(filename_appendix)
  endif

  CS%IS_energyfile = trim(slasher(directory))//trim(IS_energyfile)
  call log_param(param_file, mdl, "output_path/ENERGYFILE", CS%IS_energyfile)
#ifdef STATSLABEL
  CS%IS_energyfile = trim(CS%IS_energyfile)//"."//trim(adjustl(STATSLABEL))
#endif

  ! Allocate memory in the ice shelf dynamics control structure that was not
  ! previously allocated for registration for restarts.

  if (active_shelf_dynamics) then
    allocate( CS%t_bdry_val(isd:ied,jsd:jed), source=T_shelf_bdry) ! [C ~> degC]
    allocate( CS%u_face_mask(Isdq:Iedq,Jsdq:Jedq), source=0.0)
    allocate( CS%v_face_mask(Isdq:Iedq,Jsdq:Jedq), source=0.0)
    allocate( CS%u_flux_bdry_val(Isdq:Iedq,jsd:jed), source=0.0)
    allocate( CS%v_flux_bdry_val(isd:ied,Jsdq:Jedq), source=0.0)
    allocate( CS%umask(Isdq:Iedq,Jsdq:Jedq), source=-1.0)
    allocate( CS%vmask(Isdq:Iedq,Jsdq:Jedq), source=-1.0)
    allocate( CS%tmask(Isdq:Iedq,Jsdq:Jedq), source=-1.0)
    allocate( CS%float_cond(isd:ied,jsd:jed))

    CS%OD_rt_counter = 0
    allocate( CS%OD_rt(isd:ied,jsd:jed), source=0.0)
    allocate( CS%ground_frac_rt(isd:ied,jsd:jed), source=0.0)

    if (CS%calve_to_mask) then
      allocate( CS%calve_mask(isd:ied,jsd:jed), source=0.0)
    endif

    allocate(CS%Phi(1:8,1:4,isd:ied,jsd:jed), source=0.0)
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      call bilinear_shape_fn_grid(G, i, j, CS%Phi(:,:,i,j))
    enddo; enddo

    if (CS%GL_regularize) then
      allocate(CS%Phisub(2,2,CS%n_sub_regularize,CS%n_sub_regularize,2,2), source=0.0)
      call bilinear_shape_functions_subgrid(CS%Phisub, CS%n_sub_regularize)
    endif

    if ((trim(CS%ice_viscosity_compute) == "MODEL") .and. CS%visc_qps==1) then
      !for calculating viscosity and 1 cell-centered quadrature point per cell
      allocate(CS%PhiC(1:8,G%isc:G%iec,G%jsc:G%jec), source=0.0)
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        call bilinear_shape_fn_grid_1qp(G, i, j, CS%PhiC(:,i,j))
      enddo; enddo
    endif

    CS%elapsed_velocity_time = 0.0

    call update_velocity_masks(CS, G, ISS%hmask, CS%umask, CS%vmask, CS%u_face_mask, CS%v_face_mask)
  endif

  ! Take additional initialization steps, for example of dependent variables.
  if (active_shelf_dynamics .and. .not.new_sim) then

    call pass_var(CS%OD_av,G%domain, complete=.false.)
    call pass_var(CS%ground_frac, G%domain, complete=.false.)
    call pass_var(CS%basal_traction, G%domain, complete=.false.)
    call pass_var(CS%AGlen_visc, G%domain, complete=.false.)
    call pass_var(CS%bed_elev, G%domain, complete=.false.)
    call pass_var(CS%C_basal_friction, G%domain, complete=.false.)
    call pass_var(CS%h_bdry_val, G%domain, complete=.true.)
    call pass_var(CS%ice_visc, G%domain)

    call pass_vector(CS%u_bdry_val, CS%v_bdry_val, G%domain, TO_ALL, BGRID_NE, complete=.false.)
    call pass_vector(CS%u_face_mask_bdry, CS%v_face_mask_bdry, G%domain, TO_ALL, BGRID_NE, complete=.true.)
    call update_velocity_masks(CS, G, ISS%hmask, CS%umask, CS%vmask, CS%u_face_mask, CS%v_face_mask)

    ! This is unfortunately necessary (?); if grid is not symmetric the boundary values
    ! of u and v are otherwise not set till the end of the first linear solve, and so
    ! viscosity is not calculated correctly.
    ! This has to occur after init_boundary_values or some of the arrays on the
    ! right hand side have not been set up yet.
    if (.not. G%symmetric) then
      do j=G%jsd,G%jed ; do i=G%isd,G%ied
        if ((i+G%idg_offset) == (G%domain%nihalo+1)) then
          if (CS%u_face_mask(I-1,j) == 3) then
            CS%u_shelf(I-1,J-1) = CS%u_bdry_val(I-1,J-1)
            CS%u_shelf(I-1,J) = CS%u_bdry_val(I-1,J)
            CS%v_shelf(I-1,J-1) = CS%v_bdry_val(I-1,J-1)
            CS%v_shelf(I-1,J) = CS%v_bdry_val(I-1,J)
          elseif (CS%u_face_mask(I-1,j) == 5) then
            CS%u_shelf(I-1,J-1) = CS%u_bdry_val(I-1,J-1)
            CS%u_shelf(I-1,J) = CS%u_bdry_val(I-1,J)
          elseif (CS%u_face_mask(I-1,j) == 6) then
            CS%v_shelf(I-1,J-1) = CS%v_bdry_val(I-1,J-1)
            CS%v_shelf(I-1,J) = CS%v_bdry_val(I-1,J)
          endif
        endif
        if ((j+G%jdg_offset) == (G%domain%njhalo+1)) then
          if (CS%v_face_mask(i,J-1) == 3) then
            CS%v_shelf(I-1,J-1) = CS%v_bdry_val(I-1,J-1)
            CS%v_shelf(I,J-1) = CS%v_bdry_val(I,J-1)
            CS%u_shelf(I-1,J-1) = CS%u_bdry_val(I-1,J-1)
            CS%u_shelf(I,J-1) = CS%u_bdry_val(I,J-1)
          elseif (CS%v_face_mask(i,J-1) == 5) then
            CS%v_shelf(I-1,J-1) = CS%v_bdry_val(I-1,J-1)
            CS%v_shelf(I,J-1) = CS%v_bdry_val(I,J-1)
          elseif (CS%v_face_mask(i,J-1) == 6) then
            CS%u_shelf(I-1,J-1) = CS%u_bdry_val(I-1,J-1)
            CS%u_shelf(I,J-1) = CS%u_bdry_val(I,J-1)
          endif
        endif
      enddo ; enddo
    endif
    call pass_vector(CS%u_shelf, CS%v_shelf, G%domain, TO_ALL, BGRID_NE)
  endif

  if (active_shelf_dynamics) then
    if (CS%first_dir_restart_IS > -1.0) then
      CS%first_direction_IS = modulo(NINT(CS%first_dir_restart_IS), 2)
    else
      CS%first_dir_restart_IS = real(modulo(CS%first_direction_IS, 2))
    endif

    ! If we are calving to a mask, i.e. if a mask exists where a shelf cannot, read the mask from a file.
    if (CS%calve_to_mask) then
      call MOM_mesg("  MOM_ice_shelf.F90, initialize_ice_shelf: reading calving_mask")

      call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
      inputdir = slasher(inputdir)
      call get_param(param_file, mdl, "CALVING_MASK_FILE", IC_file, &
                   "The file with a mask for where calving might occur.", &
                   default="ice_shelf_h.nc")
      call get_param(param_file, mdl, "CALVING_MASK_VARNAME", var_name, &
                   "The variable to use in masking calving.", &
                   default="area_shelf_h")

      filename = trim(inputdir)//trim(IC_file)
      call log_param(param_file, mdl, "INPUTDIR/CALVING_MASK_FILE", filename)
      if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
         " calving mask file: Unable to open "//trim(filename))

      call MOM_read_data(filename,trim(var_name),CS%calve_mask,G%Domain)
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        if (CS%calve_mask(i,j) > 0.0) CS%calve_mask(i,j) = 1.0
      enddo ; enddo
      call pass_var(CS%calve_mask,G%domain)
    endif

    ! initialize basal friction coefficients
    if (new_sim) then
      call initialize_ice_C_basal_friction(CS%C_basal_friction, G, US, param_file)
      call pass_var(CS%C_basal_friction, G%domain, complete=.false.)

      ! initialize ice-stiffness AGlen
      call initialize_ice_AGlen(CS%AGlen_visc, CS%ice_viscosity_compute, G, US, param_file)
      call pass_var(CS%AGlen_visc, G%domain, complete=.false.)

      !initialize boundary conditions
      call initialize_ice_shelf_boundary_from_file(CS%u_face_mask_bdry, CS%v_face_mask_bdry, &
                  CS%u_bdry_val, CS%v_bdry_val, CS%umask, CS%vmask, CS%h_bdry_val, &
                  ISS%hmask,  ISS%h_shelf, G, US, param_file )
      call pass_var(ISS%hmask, G%domain, complete=.false.)
      call pass_var(CS%h_bdry_val, G%domain, complete=.true.)
      call pass_vector(CS%u_bdry_val, CS%v_bdry_val, G%domain, TO_ALL, BGRID_NE, complete=.false.)
      call pass_vector(CS%u_face_mask_bdry, CS%v_face_mask_bdry, G%domain, TO_ALL, BGRID_NE, complete=.false.)

      !initialize ice flow characteristic (velocities, bed elevation under the grounded part, etc) from file
      call initialize_ice_flow_from_file(CS%bed_elev,CS%u_shelf, CS%v_shelf, CS%ground_frac, &
                  G, US, param_file)
      call pass_vector(CS%u_shelf, CS%v_shelf, G%domain, TO_ALL, BGRID_NE, complete=.true.)
      call pass_var(CS%ground_frac, G%domain, complete=.false.)
      call pass_var(CS%bed_elev, G%domain, complete=.true.)
      call update_velocity_masks(CS, G, ISS%hmask, CS%umask, CS%vmask, CS%u_face_mask, CS%v_face_mask)

      do J=Jsdq,Jedq ; do I=Isdq,Iedq
        if (CS%umask(I,J) == 3) then
          CS%u_shelf(I,J) = CS%u_bdry_val(I,J)
        elseif (CS%umask(I,J) == 0) then
          CS%u_shelf(I,J) = 0
        endif
        if (CS%vmask(I,J) == 3) then
          CS%v_shelf(I,J) = CS%v_bdry_val(I,J)
        elseif (CS%vmask(I,J) == 0) then
          CS%v_shelf(I,J) = 0
        endif
      enddo ; enddo
    endif

  ! Register diagnostics.
    CS%id_u_shelf = register_diag_field('ice_shelf_model','u_shelf',CS%diag%axesB1, Time, &
       'x-velocity of ice', 'm yr-1', conversion=365.0*86400.0*US%L_T_to_m_s)
    CS%id_v_shelf = register_diag_field('ice_shelf_model','v_shelf',CS%diag%axesB1, Time, &
       'y-velocity of ice', 'm yr-1', conversion=365.0*86400.0*US%L_T_to_m_s)
    CS%id_shelf_speed = register_diag_field('ice_shelf_model','shelf_speed',CS%diag%axesB1, Time, &
       'speed of of ice shelf', 'm yr-1', conversion=365.0*86400.0*US%L_T_to_m_s)
    CS%id_taudx_shelf = register_diag_field('ice_shelf_model','taudx_shelf',CS%diag%axesB1, Time, &
       'x-driving stress of ice', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_taudy_shelf = register_diag_field('ice_shelf_model','taudy_shelf',CS%diag%axesB1, Time, &
       'y-driving stress of ice', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_taud_shelf = register_diag_field('ice_shelf_model','taud_shelf',CS%diag%axesB1, Time, &
       'magnitude of driving stress of ice', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_sx_shelf = register_diag_field('ice_shelf_model','sx_shelf',CS%diag%axesB1, Time, &
       'x-surface slope of ice', 'none')
    CS%id_sy_shelf = register_diag_field('ice_shelf_model','sy_shelf',CS%diag%axesB1, Time, &
       'y-surface slope of ice', 'none')
    CS%id_surf_slope_mag_shelf = register_diag_field('ice_shelf_model','surf_slope_mag_shelf', CS%diag%axesB1, Time, &
       'magnitude of surface slope of ice', 'none')
    CS%id_u_mask = register_diag_field('ice_shelf_model','u_mask',CS%diag%axesB1, Time, &
       'mask for u-nodes', 'none')
    CS%id_v_mask = register_diag_field('ice_shelf_model','v_mask',CS%diag%axesB1, Time, &
       'mask for v-nodes', 'none')
    CS%id_ground_frac = register_diag_field('ice_shelf_model','ice_ground_frac',CS%diag%axesT1, Time, &
       'fraction of cell that is grounded', 'none')
    CS%id_float_cond = register_diag_field('ice_shelf_model','float_cond',CS%diag%axesT1, Time, &
       'sub-cell grounding cells', 'none')
    CS%id_col_thick = register_diag_field('ice_shelf_model','col_thick',CS%diag%axesT1, Time, &
       'ocean column thickness passed to ice model', 'm', conversion=US%Z_to_m)
    CS%id_visc_shelf = register_diag_field('ice_shelf_model','ice_visc',CS%diag%axesT1, Time, &
       'vi-viscosity', 'Pa m s', conversion=US%RL2_T2_to_Pa*US%Z_to_m*US%T_to_s) !vertically integrated viscosity
    CS%id_taub = register_diag_field('ice_shelf_model','taub_beta',CS%diag%axesT1, Time, &
       'taub', 'MPa s m-1', conversion=1e-6*US%RL2_T2_to_Pa/(365.0*86400.0*US%L_T_to_m_s))
    CS%id_OD_av = register_diag_field('ice_shelf_model','OD_av',CS%diag%axesT1, Time, &
       'intermediate ocean column thickness passed to ice model', 'm', conversion=US%Z_to_m)

    CS%id_duHdx = register_diag_field('ice_shelf_model','duHdx',CS%diag%axesT1, Time, &
       'x-component of ice-sheet flux divergence', 'm yr-1', conversion=365.0*86400.0*US%Z_to_m*US%s_to_T)
    CS%id_dvHdy = register_diag_field('ice_shelf_model','dvHdy',CS%diag%axesT1, Time, &
       'y-component of ice-sheet flux divergence', 'm yr-1', conversion=365.0*86400.0*US%Z_to_m*US%s_to_T)
    CS%id_fluxdiv = register_diag_field('ice_shelf_model','fluxdiv',CS%diag%axesT1, Time, &
       'ice-sheet flux divergence', 'm yr-1', conversion=365.0*86400.0*US%Z_to_m*US%s_to_T)
    CS%id_strainrate_xx = register_diag_field('ice_shelf_model','strainrate_xx',CS%diag%axesT1, Time, &
       'x-component of ice-shelf strain-rate', 'yr-1', conversion=365.0*86400.0*US%s_to_T)
    CS%id_strainrate_yy = register_diag_field('ice_shelf_model','strainrate_yy',CS%diag%axesT1, Time, &
       'y-component of ice-shelf strain-rate', 'yr-1', conversion=365.0*86400.0*US%s_to_T)
    CS%id_strainrate_xy = register_diag_field('ice_shelf_model','strainrate_xy',CS%diag%axesT1, Time, &
       'xy-component of ice-shelf strain-rate', 'yr-1', conversion=365.0*86400.0*US%s_to_T)
    CS%id_pstrainrate_1 = register_diag_field('ice_shelf_model','pstrainrate_1',CS%diag%axesT1, Time, &
       'max principal horizontal ice-shelf strain-rate', 'yr-1', conversion=365.0*86400.0*US%s_to_T)
    CS%id_pstrainrate_2 = register_diag_field('ice_shelf_model','pstrainrate_2',CS%diag%axesT1, Time, &
       'min principal horizontal ice-shelf strain-rate', 'yr-1', conversion=365.0*86400.0*US%s_to_T)
    CS%id_devstress_xx = register_diag_field('ice_shelf_model','devstress_xx',CS%diag%axesT1, Time, &
       'x-component of ice-shelf deviatoric stress', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_devstress_yy = register_diag_field('ice_shelf_model','devstress_yy',CS%diag%axesT1, Time, &
       'y-component of ice-shelf deviatoric stress', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_devstress_xy = register_diag_field('ice_shelf_model','devstress_xy',CS%diag%axesT1, Time, &
       'xy-component of ice-shelf deviatoric stress', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_pdevstress_1 = register_diag_field('ice_shelf_model','pdevstress_1',CS%diag%axesT1, Time, &
       'max principal horizontal ice-shelf deviatoric stress', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_pdevstress_2 = register_diag_field('ice_shelf_model','pdevstress_2',CS%diag%axesT1, Time, &
       'min principal ice-shelf deviatoric stress', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)

    !Update these variables so that they are nonzero in case
    !IS_dynamics_post_data is called before update_ice_shelf
    if (CS%id_taudx_shelf>0 .or. CS%id_taudy_shelf>0) &
      call calc_shelf_driving_stress(CS, ISS, G, US, CS%taudx_shelf, CS%taudy_shelf, CS%OD_av)
    if (CS%id_taub>0) &
      call calc_shelf_taub(CS, ISS, G, US, CS%u_shelf, CS%v_shelf)
    if (CS%id_visc_shelf>0) &
      call calc_shelf_visc(CS, ISS, G, US, CS%u_shelf, CS%v_shelf)
  endif

  if (new_sim) then
    call update_OD_ffrac_uncoupled(CS, G, ISS%h_shelf(:,:))
  endif

end subroutine initialize_ice_shelf_dyn


subroutine initialize_diagnostic_fields(CS, ISS, G, US, Time)
  type(ice_shelf_dyn_CS), intent(inout) :: CS  !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G   !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US   !< A structure containing unit conversion factors
  type(time_type),        intent(in)    :: Time !< The current model time

  integer         :: i, j, iters, isd, ied, jsd, jed
  real            :: rhoi_rhow
  real            :: OD  ! Depth of open water below the ice shelf [Z ~> m]
  type(time_type) :: dummy_time
!
  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  dummy_time = set_time(0,0)
  isd=G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j=jsd,jed
    do i=isd,ied
      OD = CS%bed_elev(i,j) - rhoi_rhow * max(ISS%h_shelf(i,j),CS%min_h_shelf)
      if (OD >= 0) then
    ! ice thickness does not take up whole ocean column -> floating
        CS%OD_av(i,j) = OD
        CS%ground_frac(i,j) = 0.
      else
        CS%OD_av(i,j) = 0.
        CS%ground_frac(i,j) = 1.
      endif
    enddo
  enddo

  call ice_shelf_solve_outer(CS, ISS, G, US, CS%u_shelf, CS%v_shelf,CS%taudx_shelf,CS%taudy_shelf, iters, Time)
end subroutine initialize_diagnostic_fields

!> This function returns the global maximum advective timestep that can be taken based on the current
!! ice velocities.  Because it involves finding a global minimum, it can be surprisingly expensive.
function ice_time_step_CFL(CS, ISS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS  !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G   !< The grid structure used by the ice shelf.
  real :: ice_time_step_CFL !< The maximum permitted timestep based on the ice velocities [T ~> s].

  real :: dt_local, min_dt ! These should be the minimum stable timesteps at a CFL of 1 [T ~> s]
  real :: min_vel          ! A minimal velocity for estimating a timestep [L T-1 ~> m s-1]
  integer :: i, j

  min_dt = 5.0e17*G%US%s_to_T ! The starting maximum is roughly the lifetime of the universe.
  min_vel = (1.0e-12/(365.0*86400.0)) * G%US%m_s_to_L_T
  do j=G%jsc,G%jec ; do i=G%isc,G%iec ; if (ISS%hmask(i,j) == 1.0 .or. ISS%hmask(i,j)==3) then
    dt_local = 2.0*G%areaT(i,j) / &
       (((G%dyCu(I,j)  * max(abs(CS%u_shelf(I,J)  + CS%u_shelf(I,j-1)), min_vel)) + &
         (G%dyCu(I-1,j)* max(abs(CS%u_shelf(I-1,J)+ CS%u_shelf(I-1,j-1)), min_vel))) + &
        ((G%dxCv(i,J)  * max(abs(CS%v_shelf(i,J)  + CS%v_shelf(i-1,J)), min_vel)) + &
         (G%dxCv(i,J-1)* max(abs(CS%v_shelf(i,J-1)+ CS%v_shelf(i-1,J-1)), min_vel))))

    min_dt = min(min_dt, dt_local)
  endif ; enddo ; enddo ! i- and j- loops

  call min_across_PEs(min_dt)

  ice_time_step_CFL = CS%CFL_factor * min_dt

end function ice_time_step_CFL

!> This subroutine updates the ice shelf velocities, mass, stresses and properties due to the
!! ice shelf dynamics.
subroutine update_ice_shelf(CS, ISS, G, US, time_step, Time, calve_ice_shelf_bergs, &
                            ocean_mass, coupled_grounding, must_update_vel)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                              !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real,                   intent(in)    :: time_step !< time step [T ~> s]
  type(time_type),        intent(in)    :: Time !< The current model time
  logical,                intent(in)    :: calve_ice_shelf_bergs !< To convert ice flux through front
                                                                 !! to bergs
  real, dimension(SZDI_(G),SZDJ_(G)), &
                optional, intent(in)    :: ocean_mass !< If present this is the mass per unit area
                                              !! of the ocean [R Z ~> kg m-2].
  logical,      optional, intent(in)    :: coupled_grounding !< If true, the grounding line is
                                              !! determined by coupled ice-ocean dynamics
  logical,      optional, intent(in)    :: must_update_vel !< Always update the ice velocities if true.
  integer :: iters
  logical :: update_ice_vel, coupled_GL

  update_ice_vel = .false.
  if (present(must_update_vel)) update_ice_vel = must_update_vel

  coupled_GL = .false.
  if (present(ocean_mass) .and. present(coupled_grounding)) coupled_GL = coupled_grounding
!
  if (CS%advect_shelf) then
    call ice_shelf_advect(CS, ISS, G, time_step, Time, calve_ice_shelf_bergs)
    if (CS%alternate_first_direction_IS) then
      CS%first_direction_IS = modulo(CS%first_direction_IS+1,2)
      CS%first_dir_restart_IS = real(CS%first_direction_IS)
    endif
  endif
  CS%elapsed_velocity_time = CS%elapsed_velocity_time + time_step
  if (CS%elapsed_velocity_time >= CS%velocity_update_time_step) update_ice_vel = .true.

  if (coupled_GL) then
    call update_OD_ffrac(CS, G, US, ocean_mass, update_ice_vel)
  elseif (update_ice_vel) then
    call update_OD_ffrac_uncoupled(CS, G, ISS%h_shelf(:,:))
    CS%GL_couple=.false.
  endif

  if (update_ice_vel) then
    call ice_shelf_solve_outer(CS, ISS, G, US, CS%u_shelf, CS%v_shelf,CS%taudx_shelf,CS%taudy_shelf, iters, Time)
    CS%elapsed_velocity_time = 0.0
  endif

! call ice_shelf_temp(CS, ISS, G, US, time_step, ISS%water_flux, Time)

end subroutine update_ice_shelf

subroutine volume_above_floatation(CS, G, ISS, vaf, hemisphere)
  type(ice_shelf_dyn_CS), intent(in) :: CS !< The ice shelf dynamics control structure
  type(ocean_grid_type),  intent(in) :: G  !< The grid structure used by the ice shelf.
  type(ice_shelf_state),  intent(in) :: ISS !< A structure with elements that describe
                                            !! the ice-shelf state
  real, intent(out) :: vaf !< area integrated volume above floatation [m3]
  integer, optional, intent(in) :: hemisphere !< 0 for Antarctica only, 1 for Greenland only. Otherwise, all ice sheets
  integer :: IS_ID ! local copy of hemisphere
  real, dimension(SZI_(G),SZJ_(G))  :: vaf_cell !< cell-wise volume above floatation [m3]
  integer, dimension(SZI_(G),SZJ_(G))  :: mask ! a mask for active cells depending on hemisphere indicated
  integer :: is,ie,js,je,i,j
  real :: rhoi_rhow, rhow_rhoi

  if (CS%GL_couple) &
    call MOM_error(FATAL, "MOM_ice_shelf_dyn, volume above floatation calculation assumes GL_couple=.FALSE..")

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  rhow_rhoi = CS%density_ocean_avg / CS%density_ice
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (present(hemisphere)) then
    IS_ID=hemisphere
  else
    IS_ID=-1
  endif

  mask(:,:)=0
  if (IS_ID==0) then     !Antarctica (S. Hemisphere) only
    do j = js,je; do i = is,ie
      if (ISS%hmask(i,j)>0 .and. G%geoLatT(i,j)<=0.0) mask(i,j)=1
    enddo; enddo
  elseif (IS_ID==1) then !Greenland (N. Hemisphere) only
    do j = js,je; do i = is,ie
      if (ISS%hmask(i,j)>0 .and. G%geoLatT(i,j)>0.0)  mask(i,j)=1
    enddo; enddo
  else                   !All ice sheets
    mask(is:ie,js:je)=ISS%hmask(is:ie,js:je)
  endif

  vaf_cell(:,:)=0.0
  do j = js,je; do i = is,ie
    if (mask(i,j)>0) then
      if (CS%bed_elev(i,j) <= 0) then
        !grounded above sea level
        vaf_cell(i,j)= (ISS%h_shelf(i,j) * G%US%Z_to_m) * (ISS%area_shelf_h(i,j) * G%US%L_to_m**2)
      else
        !grounded if vaf_cell(i,j) > 0
        vaf_cell(i,j) = (max(ISS%h_shelf(i,j) - rhow_rhoi * CS%bed_elev(i,j), 0.0) * G%US%Z_to_m) * &
                        (ISS%area_shelf_h(i,j) * G%US%L_to_m**2)
      endif
    endif
  enddo; enddo

  vaf = reproducing_sum(vaf_cell)
end subroutine volume_above_floatation

!> multiplies a variable with the ice sheet grounding fraction
subroutine masked_var_grounded(G,CS,var,varout)
  type(ocean_grid_type), intent(in) :: G !< The grid structure used by the ice shelf.
  type(ice_shelf_dyn_CS), intent(in) :: CS !< The ice shelf dynamics control structure
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: var !< variable in
  real, dimension(SZI_(G),SZJ_(G)), intent(out)  :: varout !<variable out
  integer :: i,j
  do j = G%jsc,G%jec; do i = G%isc,G%iec
      varout(i,j) = var(i,j) * CS%ground_frac(i,j)
  enddo; enddo
end subroutine masked_var_grounded

!> Ice shelf dynamics post_data calls
subroutine IS_dynamics_post_data(time_step, Time, CS, ISS, G)
  real :: time_step !< Length of time for post data averaging [T ~> s].
  type(time_type),        intent(in)    :: Time !< The current model time
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDIB_(G),SZDJB_(G))  :: taud_x, taud_y, taud  ! area-averaged driving stress [R L2 T-2 ~> Pa]
  real, dimension(SZDI_(G),SZDJ_(G))  :: ice_visc ! area-averaged vertically integrated ice viscosity
                                                  !! [R L2 Z T-1 ~> Pa s m]
  real, dimension(SZDI_(G),SZDJ_(G))  :: basal_tr ! area-averaged taub_beta field related to basal traction,
                                                  !! [R L1 T-1 ~> Pa s m-1]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: surf_slope ! the surface slope of the ice shelf/sheet [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: ice_speed ! ice sheet flow speed [L T-1 ~> m s-1]

  integer :: i,j
    call enable_averages(time_step, Time, CS%diag)
    if (CS%id_col_thick > 0) call post_data(CS%id_col_thick, CS%OD_av, CS%diag)
    if (CS%id_u_shelf > 0) call post_data(CS%id_u_shelf, CS%u_shelf, CS%diag)
    if (CS%id_v_shelf > 0) call post_data(CS%id_v_shelf, CS%v_shelf, CS%diag)
    if (CS%id_shelf_speed > 0) then
      do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
        ice_speed(I,J) = sqrt((CS%u_shelf(I,J)**2) + (CS%v_shelf(I,J)**2))
      enddo ; enddo
      call post_data(CS%id_shelf_speed, ice_speed, CS%diag)
    endif
!   if (CS%id_t_shelf > 0) call post_data(CS%id_t_shelf, CS%t_shelf, CS%diag)
    if (CS%id_taudx_shelf > 0) then
      do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
        taud_x(I,J) = CS%taudx_shelf(I,J)*G%IareaBu(I,J)
      enddo ; enddo
      call post_data(CS%id_taudx_shelf, taud_x, CS%diag)
    endif
    if (CS%id_taudy_shelf > 0) then
      do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
        taud_y(I,J) = CS%taudy_shelf(I,J)*G%IareaBu(I,J)
      enddo ; enddo
      call post_data(CS%id_taudy_shelf, taud_y, CS%diag)
    endif
    if (CS%id_taud_shelf > 0) then
      do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
        taud(I,J) = sqrt((CS%taudx_shelf(I,J)**2)+(CS%taudy_shelf(I,J)**2))*G%IareaBu(I,J)
      enddo ; enddo
      call post_data(CS%id_taud_shelf, taud, CS%diag)
    endif
    if (CS%id_sx_shelf > 0) call post_data(CS%id_sx_shelf, CS%sx_shelf, CS%diag)
    if (CS%id_sy_shelf > 0) call post_data(CS%id_sy_shelf, CS%sy_shelf, CS%diag)
    if (CS%id_surf_slope_mag_shelf > 0) then
      do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
        surf_slope(I,J) = sqrt((CS%sx_shelf(I,J)**2)+(CS%sy_shelf(I,J)**2))
      enddo ; enddo
      call post_data(CS%id_surf_slope_mag_shelf, surf_slope, CS%diag)
    endif
    if (CS%id_ground_frac > 0) call post_data(CS%id_ground_frac, CS%ground_frac, CS%diag)
    if (CS%id_float_cond > 0) call post_data(CS%id_float_cond, CS%float_cond, CS%diag)
    if (CS%id_OD_av >0) call post_data(CS%id_OD_av, CS%OD_av,CS%diag)
    if (CS%id_visc_shelf > 0) then
      call ice_visc_diag(CS,G,ice_visc)
      call post_data(CS%id_visc_shelf, ice_visc, CS%diag)
    endif
    if (CS%id_taub > 0) then
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        basal_tr(i,j) = CS%basal_traction(i,j)*G%IareaT(i,j)
      enddo ; enddo
      call post_data(CS%id_taub, basal_tr, CS%diag)
    endif
    if (CS%id_u_mask > 0) call post_data(CS%id_u_mask, CS%umask, CS%diag)
    if (CS%id_v_mask > 0) call post_data(CS%id_v_mask, CS%vmask, CS%diag)
    if (CS%id_ufb_mask > 0) call post_data(CS%id_ufb_mask, CS%u_face_mask_bdry, CS%diag)
    if (CS%id_vfb_mask > 0) call post_data(CS%id_vfb_mask, CS%v_face_mask_bdry, CS%diag)
!   if (CS%id_t_mask > 0) call post_data(CS%id_t_mask, CS%tmask, CS%diag)

    if (CS%id_duHdx > 0         .or. CS%id_dvHdy > 0         .or. CS%id_fluxdiv > 0       .or. &
        CS%id_devstress_xx > 0  .or. CS%id_devstress_yy > 0  .or. CS%id_devstress_xy > 0  .or. &
        CS%id_strainrate_xx > 0 .or. CS%id_strainrate_yy > 0 .or. CS%id_strainrate_xy > 0 .or. &
        CS%id_pdevstress_1 > 0  .or. CS%id_pdevstress_2 > 0  .or. &
        CS%id_pstrainrate_1 > 0 .or. CS%id_pstrainrate_2 > 0) then
      call IS_dynamics_post_data_2(CS, ISS, G)
    endif

    call disable_averaging(CS%diag)
end subroutine IS_dynamics_post_data

!> Calculate cell-centered, area-averaged, vertically integrated ice viscosity for diagnostics
subroutine ice_visc_diag(CS,G,ice_visc)
  type(ice_shelf_dyn_CS), intent(in) :: CS !< The ice shelf dynamics control structure
  type(ocean_grid_type),  intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), intent(out)  :: ice_visc !< area-averaged vertically integrated ice viscosity
                                                               !! [R L2 Z T-1 ~> Pa s m]
  integer :: i,j

  ice_visc(:,:)=0.0
  if (CS%visc_qps==4) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      ice_visc(i,j) = (0.25 * G%IareaT(i,j)) * &
        ((CS%ice_visc(i,j,1) + CS%ice_visc(i,j,4)) + (CS%ice_visc(i,j,2) + CS%ice_visc(i,j,3)))
    enddo ; enddo
  else
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      ice_visc(i,j) = CS%ice_visc(i,j,1)*G%IareaT(i,j)
    enddo ; enddo
  endif
end subroutine ice_visc_diag

!>  Writes the total ice shelf kinetic energy and mass to an ascii file
subroutine write_ice_shelf_energy(CS, G, US, mass, area, day, time_step)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: mass !< The mass per unit area of the ice shelf
                                                !! or sheet [R Z ~> kg m-2]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(in)    :: area !< The ice shelf or ice sheet area [L2 ~> m2]
  type(time_type),         intent(in)    :: day !< The current model time.
  type(time_type),  optional, intent(in) :: time_step !< The current time step
  ! Local variables
  type(time_type) :: dt ! A time_type version of the timestep.
  real, dimension(SZDI_(G),SZDJ_(G)) :: tmp1 ! A temporary array used in reproducing sums [various]
  real :: KE_tot, mass_tot, KE_scale_factor, mass_scale_factor
  integer :: is, ie, js, je, isr, ier, jsr, jer, i, j
  character(len=32)  :: mesg_intro, time_units, day_str, n_str, date_str
  integer :: start_of_day, num_days
  real    :: reday  ! Time in units given by CS%Timeunit, but often [days]

  ! write_energy_time is the next integral multiple of energysavedays.
  if (present(time_step)) then
    dt = time_step
  else
    dt = set_time(seconds=2)
  endif

   !CS%prev_IS_energy_calls tracks the ice sheet step, which is outputted in the energy file.
  if (CS%prev_IS_energy_calls == 0) then
    if (CS%energysave_geometric) then
      if (CS%energysavedays_geometric < CS%energysavedays) then
        CS%write_energy_time = day + CS%energysavedays_geometric
        CS%geometric_end_time = CS%Start_time + CS%energysavedays * &
          (1 + (day - CS%Start_time) / CS%energysavedays)
      else
        CS%write_energy_time = CS%Start_time + CS%energysavedays * &
          (1 + (day - CS%Start_time) / CS%energysavedays)
      endif
    else
      CS%write_energy_time = CS%Start_time + CS%energysavedays * &
        (1 + (day - CS%Start_time) / CS%energysavedays)
    endif
  elseif (day + (dt/2) <= CS%write_energy_time) then
    CS%prev_IS_energy_calls = CS%prev_IS_energy_calls + 1
    return  ! Do not write this step
  else ! Determine the next write time before proceeding
    if (CS%energysave_geometric) then
      if (CS%write_energy_time + CS%energysavedays_geometric >= &
          CS%geometric_end_time) then
        CS%write_energy_time = CS%geometric_end_time
        CS%energysave_geometric = .false.  ! stop geometric progression
      else
        CS%write_energy_time = CS%write_energy_time + CS%energysavedays_geometric
      endif
      CS%energysavedays_geometric = CS%energysavedays_geometric*2
    else
      CS%write_energy_time = CS%write_energy_time + CS%energysavedays
    endif
  endif

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isr = is - (G%isd-1) ; ier = ie - (G%isd-1) ; jsr = js - (G%jsd-1) ; jer = je - (G%jsd-1)

  !calculate KE using cell-centered ice shelf velocity
  tmp1(:,:)=0.0
  KE_scale_factor = US%L_to_m**2 * (US%RZ_to_kg_m2 * US%L_T_to_m_s**2)
  do j=js,je ; do i=is,ie
    tmp1(i,j) = (KE_scale_factor * 0.03125) * (mass(i,j) * area(i,j)) * &
      ((((CS%u_shelf(I-1,J-1)+CS%u_shelf(I,J))+(CS%u_shelf(I,J-1)+CS%u_shelf(I-1,J)))**2) + &
       (((CS%v_shelf(I-1,J-1)+CS%v_shelf(I,J))+(CS%v_shelf(I,J-1)+CS%v_shelf(I-1,J)))**2))
  enddo; enddo

  KE_tot = reproducing_sum(tmp1, isr, ier, jsr, jer)

  !calculate mass
  tmp1(:,:)=0.0
  mass_scale_factor = US%L_to_m**2 * US%RZ_to_kg_m2
  do j=js,je ; do i=is,ie
    tmp1(i,j) =  mass_scale_factor * (mass(i,j) * area(i,j))
  enddo; enddo

  mass_tot = reproducing_sum(tmp1, isr, ier, jsr, jer)

  if (is_root_pe()) then  ! Only the root PE actually writes anything.
    if (day > CS%Start_time) then
      call open_ASCII_file(CS%IS_fileenergy_ascii, trim(CS%IS_energyfile), action=APPEND_FILE)
    else
      call open_ASCII_file(CS%IS_fileenergy_ascii, trim(CS%IS_energyfile), action=WRITEONLY_FILE)
      if (abs(CS%timeunit - 86400.0) < 1.0) then
        write(CS%IS_fileenergy_ascii,'("  Step,",7x,"Day,",8x,"Energy/Mass,",13x,"Total Mass")')
        write(CS%IS_fileenergy_ascii,'(12x,"[days]",10x,"[m2 s-2]",17x,"[kg]")')
      else
        if ((CS%timeunit >= 0.99) .and. (CS%timeunit < 1.01)) then
          time_units = "           [seconds]     "
        elseif ((CS%timeunit >= 3599.0) .and. (CS%timeunit < 3601.0)) then
          time_units = "            [hours]      "
        elseif ((CS%timeunit >= 86399.0) .and. (CS%timeunit < 86401.0)) then
          time_units = "             [days]      "
        elseif ((CS%timeunit >= 3.0e7) .and. (CS%timeunit < 3.2e7)) then
          time_units = "            [years]      "
        else
          write(time_units,'(9x,"[",es8.2," s]    ")') CS%timeunit
        endif

        write(CS%IS_fileenergy_ascii,'("  Step,",7x,"Time,",7x,"Energy/Mass,",13x,"Total Mass")')
        write(CS%IS_fileenergy_ascii,'(A25,3x,"[m2 s-2]",17x,"[kg]")') time_units
      endif
    endif

    call get_time(day, start_of_day, num_days)

    if (abs(CS%timeunit - 86400.0) < 1.0) then
      reday = REAL(num_days)+ (REAL(start_of_day)/86400.0)
    else
      reday = REAL(num_days)*(86400.0/CS%timeunit) + REAL(start_of_day)/abs(CS%timeunit)
    endif

    if (reday < 1.0e8) then ;      write(day_str, '(F12.3)') reday
    elseif (reday < 1.0e11) then ; write(day_str, '(F15.3)') reday
    else ;                         write(day_str, '(ES15.9)') reday ; endif

    if     (CS%prev_IS_energy_calls < 1000000)   then ; write(n_str, '(I6)') CS%prev_IS_energy_calls
    elseif (CS%prev_IS_energy_calls < 10000000)  then ; write(n_str, '(I7)') CS%prev_IS_energy_calls
    elseif (CS%prev_IS_energy_calls < 100000000) then ; write(n_str, '(I8)') CS%prev_IS_energy_calls
    else                        ; write(n_str, '(I10)') CS%prev_IS_energy_calls ; endif

    write(CS%IS_fileenergy_ascii,'(A,",",A,", En ",ES22.16,", M ",ES11.5)') &
      trim(n_str), trim(day_str), KE_tot/mass_tot, mass_tot
  endif

  CS%prev_IS_energy_calls = CS%prev_IS_energy_calls + 1
end subroutine write_ice_shelf_energy

!> This subroutine takes the velocity (on the Bgrid) and timesteps h_t = - div (uh) once.
!! Additionally, it will update the volume of ice in partially-filled cells, and update
!! hmask accordingly
subroutine ice_shelf_advect(CS, ISS, G, time_step, Time, calve_ice_shelf_bergs)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  real,                   intent(in)    :: time_step !< time step [T ~> s]
  type(time_type),        intent(in)    :: Time !< The current model time
  logical,                intent(in)    :: calve_ice_shelf_bergs !< If true, track ice shelf flux through a
                                               !! static ice shelf, so that it can be converted into icebergs

! 3/8/11 DNG
!
!    This subroutine takes the velocity (on the Bgrid) and timesteps h_t = - div (uh) once.
!    ADDITIONALLY, it will update the volume of ice in partially-filled cells, and update
!        hmask accordingly
!
!    The flux overflows are included here. That is because they will be used to advect 3D scalars
!    into partial cells

  real, dimension(SZDI_(G),SZDJ_(G))   :: h_after_flux1, h_after_flux2 ! Ice thicknesses [Z ~> m].
  real, dimension(SZDIB_(G),SZDJ_(G))  :: uh_ice  ! The accumulated zonal ice volume flux [Z L2 ~> m3]
  real, dimension(SZDI_(G),SZDJB_(G))  :: vh_ice  ! The accumulated meridional ice volume flux [Z L2 ~> m3]
  type(loop_bounds_type) :: LB
  integer                           :: isd, ied, jsd, jed, i, j, isc, iec, jsc, jec, stencil

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  uh_ice(:,:) = 0.0
  vh_ice(:,:) = 0.0

  h_after_flux1(:,:) = 0.0
  h_after_flux2(:,:) = 0.0
  ! call MOM_mesg("MOM_ice_shelf.F90: ice_shelf_advect called")

  do j=jsd,jed ; do i=isd,ied ; if (CS%h_bdry_val(i,j) /= 0.0) then
    ISS%h_shelf(i,j) = CS%h_bdry_val(i,j)
  endif ; enddo ; enddo

  stencil = 2
  if (modulo(CS%first_direction_IS,2)==0) then
    !x first
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc-stencil ; LB%jeh = G%jec+stencil
    if (LB%jsh < jsd) call MOM_error(FATAL, &
      "ice_shelf_advect:  Halo is too small for the ice thickness advection stencil.")
    call ice_shelf_advect_thickness_x(CS, G, LB, time_step, ISS%hmask, ISS%h_shelf, h_after_flux1, uh_ice)
    call pass_var(h_after_flux1, G%domain)
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
    call ice_shelf_advect_thickness_y(CS, G, LB, time_step, ISS%hmask, h_after_flux1, h_after_flux2, vh_ice)
  else
    ! y first
    LB%ish = G%isc-stencil ; LB%ieh = G%iec+stencil ; LB%jsh = G%jsc ; LB%jeh = G%jec
    if (LB%ish < isd) call MOM_error(FATAL, &
      "ice_shelf_advect:  Halo is too small for the ice thickness advection stencil.")
    call ice_shelf_advect_thickness_y(CS, G, LB, time_step, ISS%hmask, ISS%h_shelf, h_after_flux1, vh_ice)
    call pass_var(h_after_flux1, G%domain)
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
    call ice_shelf_advect_thickness_x(CS, G, LB, time_step, ISS%hmask, h_after_flux1, h_after_flux2, uh_ice)
  endif
  call pass_var(h_after_flux2, G%domain)

  do j=jsd,jed
    do i=isd,ied
      if (ISS%hmask(i,j) == 1) ISS%h_shelf(i,j) = h_after_flux2(i,j)
    enddo
  enddo

  if (CS%moving_shelf_front) then
    call shelf_advance_front(CS, ISS, G, ISS%hmask, uh_ice, vh_ice)
    if (CS%min_thickness_simple_calve > 0.0) then
      call ice_shelf_min_thickness_calve(G, ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, &
                                         CS%min_thickness_simple_calve)
    endif
    if (CS%calve_to_mask) then
      call calve_to_mask(G, ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, CS%calve_mask)
    endif
  elseif (calve_ice_shelf_bergs) then
    !advect the front to create partially-filled cells
    call shelf_advance_front(CS, ISS, G, ISS%hmask, uh_ice, vh_ice)
    !add mass of the partially-filled cells to calving field, which is used to initialize icebergs
    !Then, remove the partially-filled cells from the ice shelf
    ISS%calving(:,:)=0.0
    ISS%calving_hflx(:,:)=0.0
    do j=jsc,jec; do i=isc,iec
      if (ISS%hmask(i,j)==2) then
        ISS%calving(i,j) = (ISS%h_shelf(i,j) * CS%density_ice) * &
                           (ISS%area_shelf_h(i,j) * G%IareaT(i,j)) / time_step
        ISS%calving_hflx(i,j) = (CS%Cp_ice * CS%t_shelf(i,j)) * &
                                ((ISS%h_shelf(i,j) * CS%density_ice) * &
                                (ISS%area_shelf_h(i,j) * G%IareaT(i,j)))
        ISS%h_shelf(i,j) = 0.0; ISS%area_shelf_h(i,j) = 0.0; ISS%hmask(i,j) = 0.0
      endif
    enddo; enddo
  endif

  do j=jsc,jec; do i=isc,iec
    ISS%mass_shelf(i,j) = ISS%h_shelf(i,j) * CS%density_ice
  enddo; enddo

  call pass_var(ISS%mass_shelf, G%domain, complete=.false.)
  call pass_var(ISS%h_shelf, G%domain, complete=.false.)
  call pass_var(ISS%area_shelf_h, G%domain, complete=.false.)
  call pass_var(ISS%hmask, G%domain, complete=.true.)

  call update_velocity_masks(CS, G, ISS%hmask, CS%umask, CS%vmask, CS%u_face_mask, CS%v_face_mask)

end subroutine ice_shelf_advect

!>This subroutine computes u- and v-velocities of the ice shelf iterating on non-linear ice viscosity
!subroutine ice_shelf_solve_outer(CS, ISS, G, US, u_shlf, v_shlf, iters, time)
subroutine ice_shelf_solve_outer(CS, ISS, G, US, u_shlf, v_shlf, taudx, taudy, iters, Time)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_shlf  !< The zonal ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_shlf  !< The meridional ice shelf velocity at vertices [L T-1 ~> m s-1]
  integer,                intent(out)   :: iters !< The number of iterations used in the solver.
  type(time_type),        intent(in)    :: Time !< The current model time

  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: taudx !< Driving x-stress at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: taudy !< Driving y-stress at q-points [R L3 Z T-2 ~> kg m s-2]
  !real, dimension(SZDIB_(G),SZDJB_(G)) :: u_bdry_cont ! Boundary u-stress contribution [R L3 Z T-2 ~> kg m s-2]
  !real, dimension(SZDIB_(G),SZDJB_(G)) :: v_bdry_cont ! Boundary v-stress contribution [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: Au, Av ! The retarding lateral stress contributions [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: u_last, v_last ! Previous velocities [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: H_node ! Ice shelf thickness at corners [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)) :: float_cond ! If GL_regularize=true, indicates cells containing
                                                ! the grounding line (float_cond=1) or not (float_cond=0)
  real, dimension(SZDIB_(G),SZDJB_(G)) :: Normvec  ! Used for convergence
  character(len=160) :: mesg  ! The text of an error message
  integer :: conv_flag, i, j, k,l, iter, nodefloat
  integer :: Isdq, Iedq, Jsdq, Jedq, isd, ied, jsd, jed
  integer :: Iscq, Iecq, Jscq, Jecq, isc, iec, jsc, jec
  real    :: err_max, err_tempu, err_tempv, err_init, max_vel, tempu, tempv, Norm, PrevNorm
  real    :: rhoi_rhow ! The density of ice divided by a typical water density [nondim]
  integer :: Is_sum, Js_sum, Ie_sum, Je_sum ! Loop bounds for global sums or arrays starting at 1.
  integer :: Iscq_sv, Jscq_sv ! Starting loop bound for sum_vec

  Isdq = G%IsdB ; Iedq = G%IedB ; Jsdq = G%JsdB ; Jedq = G%JedB
  Iscq = G%IscB ; Iecq = G%IecB ; Jscq = G%JscB ; Jecq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  rhoi_rhow = CS%density_ice / CS%density_ocean_avg

  taudx(:,:) = 0.0 ; taudy(:,:) = 0.0
  Au(:,:) = 0.0 ; Av(:,:) = 0.0

  ! need to make these conditional on GL interpolation
  CS%float_cond(:,:) = 0.0 ; H_node(:,:) = 0.0
  !CS%ground_frac(:,:) = 0.0

  if (.not. CS%GL_couple) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (rhoi_rhow * max(ISS%h_shelf(i,j),CS%min_h_shelf) - CS%bed_elev(i,j) > 0) then
        CS%ground_frac(i,j) = 1.0
        CS%OD_av(i,j) =0.0
      endif
    enddo ; enddo
  endif

  call calc_shelf_driving_stress(CS, ISS, G, US, taudx, taudy, CS%OD_av)
  call pass_vector(taudx, taudy, G%domain, TO_ALL, BGRID_NE)
  ! This is to determine which cells contain the grounding line, the criterion being that the cell
  ! is ice-covered, with some nodes floating and some grounded flotation condition is estimated by
  ! assuming topography is cellwise constant and H is bilinear in a cell; floating where
  ! rho_i/rho_w * H_node - D is negative

  ! need to make this conditional on GL interp

  if (CS%GL_regularize) then

    call interpolate_H_to_B(G, ISS%h_shelf, ISS%hmask, H_node, CS%min_h_shelf)

    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      nodefloat = 0

      do l=0,1 ; do k=0,1
        if ((ISS%hmask(i,j) == 1 .or. ISS%hmask(i,j)==3) .and. &
            (rhoi_rhow * H_node(i-1+k,j-1+l) - CS%bed_elev(i,j) <= 0)) then
          nodefloat = nodefloat + 1
        endif
      enddo ; enddo
      if ((nodefloat > 0) .and. (nodefloat < 4)) then
        CS%float_cond(i,j) = 1.0
        CS%ground_frac(i,j) = 1.0
      endif
    enddo ; enddo

    call pass_var(CS%float_cond, G%Domain, complete=.false.)
    call pass_var(CS%ground_frac, G%domain, complete=.false.)

  endif

  call calc_shelf_taub(CS, ISS, G, US, u_shlf, v_shlf)
  call pass_var(CS%basal_traction, G%domain, complete=.true.)
  call calc_shelf_visc(CS, ISS, G, US, u_shlf, v_shlf)
  call pass_var(CS%ice_visc, G%domain)

  ! This makes sure basal stress is only applied when it is supposed to be
  if (CS%GL_regularize) then
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      if (CS%ground_frac(i,j)/=1.0) CS%basal_traction(i,j) = 0.0
    enddo ; enddo
  else
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      CS%basal_traction(i,j) = CS%basal_traction(i,j) * CS%ground_frac(i,j)
    enddo ; enddo
  endif

  if (CS%nonlin_solve_err_mode == 1) then

    Au(:,:) = 0.0 ; Av(:,:) = 0.0

    call CG_action(CS, Au, Av, u_shlf, v_shlf, CS%Phi, CS%Phisub, CS%umask, CS%vmask, ISS%hmask, H_node, &
                   CS%ice_visc, CS%float_cond, CS%bed_elev, CS%basal_traction, &
                   G, US, G%isc-1, G%iec+1, G%jsc-1, G%jec+1, rhoi_rhow)
    call pass_vector(Au, Av, G%domain, TO_ALL, BGRID_NE)

    err_init = 0 ; err_tempu = 0 ; err_tempv = 0
    do J=G%JscB,G%JecB ; do I=G%IscB,G%IecB
      if (CS%umask(I,J) == 1) then
        err_tempu = ABS(Au(I,J) - taudx(I,J))
        if (err_tempu >= err_init) err_init = err_tempu
      endif
      if (CS%vmask(I,J) == 1) then
        err_tempv = ABS(Av(I,J) - taudy(I,J))
        if (err_tempv >= err_init) err_init = err_tempv
      endif
    enddo ; enddo

    call max_across_PEs(err_init)
  elseif (CS%nonlin_solve_err_mode == 3) then
    Normvec=0.0

    ! Determine the loop limits for sums, bearing in mind that the arrays will be starting at 1.
    ! Includes the edge of the tile is at the western/southern bdry (if symmetric)
    if ((isc+G%idg_offset==G%isg) .and. (.not. CS%reentrant_x)) then
      Is_sum = Iscq + (1-Isdq) ; Iscq_sv = Iscq
    else
      Is_sum = isc  + (1-Isdq) ; Iscq_sv = isc
    endif
    if ((jsc+G%jdg_offset==G%jsg) .and. (.not. CS%reentrant_y)) then
      Js_sum = Jscq + (1-Jsdq) ; Jscq_sv = Jscq
    else
      Js_sum = jsc + (1-Jsdq) ; Jscq_sv = jsc
    endif
    Ie_sum = Iecq + (1-Isdq) ; Je_sum = Jecq + (1-Jsdq)

    do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
      if (CS%umask(I,J) == 1) Normvec(I,J) = Normvec(I,J) + (u_shlf(I,J)**2 * US%L_T_to_m_s**2)
      if (CS%vmask(I,J) == 1) Normvec(I,J) = Normvec(I,J) + (v_shlf(I,J)**2 * US%L_T_to_m_s**2)
    enddo ; enddo
    Norm = reproducing_sum( Normvec, Is_sum, Ie_sum, Js_sum, Je_sum )
    Norm = sqrt(Norm)
  endif

  u_last(:,:) = u_shlf(:,:) ; v_last(:,:) = v_shlf(:,:)

  !! begin loop

  do iter=1,50

    call ice_shelf_solve_inner(CS, ISS, G, US, u_shlf, v_shlf, taudx, taudy, H_node, CS%float_cond, &
                               ISS%hmask, conv_flag, iters, time, CS%Phi, CS%Phisub)

    if (CS%debug) then
      call qchksum(u_shlf, "u shelf", G%HI, haloshift=2, unscale=US%L_T_to_m_s)
      call qchksum(v_shlf, "v shelf", G%HI, haloshift=2, unscale=US%L_T_to_m_s)
    endif

    write(mesg,*) "ice_shelf_solve_outer: linear solve done in ",iters," iterations"
    call MOM_mesg(mesg, 5)

    call calc_shelf_taub(CS, ISS, G, US, u_shlf, v_shlf)
    call pass_var(CS%basal_traction, G%domain, complete=.true.)
    call calc_shelf_visc(CS, ISS, G, US, u_shlf, v_shlf)
    call pass_var(CS%ice_visc, G%domain)

    ! makes sure basal stress is only applied when it is supposed to be
    if (CS%GL_regularize) then
      do j=G%jsd,G%jed ; do i=G%isd,G%ied
        if (CS%ground_frac(i,j)/=1.0) CS%basal_traction(i,j) = 0.0
      enddo ; enddo
    else
      do j=G%jsd,G%jed ; do i=G%isd,G%ied
        CS%basal_traction(i,j) = CS%basal_traction(i,j) * CS%ground_frac(i,j)
      enddo ; enddo
    endif

    if (CS%nonlin_solve_err_mode == 1) then

      Au(:,:) = 0 ; Av(:,:) = 0

      call CG_action(CS, Au, Av, u_shlf, v_shlf, CS%Phi, CS%Phisub, CS%umask, CS%vmask, ISS%hmask, H_node, &
                     CS%ice_visc, CS%float_cond, CS%bed_elev, CS%basal_traction, &
                     G, US, G%isc-1, G%iec+1, G%jsc-1, G%jec+1, rhoi_rhow)

      call pass_vector(Au, Av, G%domain, TO_ALL, BGRID_NE)

      err_max = 0

      do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
        if (CS%umask(I,J) == 1) then
          err_tempu = ABS(Au(I,J) - taudx(I,J))
          if (err_tempu >= err_max) err_max = err_tempu
        endif
        if (CS%vmask(I,J) == 1) then
          err_tempv = ABS(Av(I,J) - taudy(I,J))
          if (err_tempv >= err_max) err_max = err_tempv
        endif
      enddo ; enddo

      call max_across_PEs(err_max)

    elseif (CS%nonlin_solve_err_mode == 2) then

      err_max=0. ;  max_vel = 0 ; tempu = 0 ; tempv = 0 ; err_tempu = 0
      do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
        if (CS%umask(I,J) == 1) then
          err_tempu = ABS(u_last(I,J)-u_shlf(I,J))
          if (err_tempu >= err_max) err_max = err_tempu
          tempu = u_shlf(I,J)
        else
          tempu = 0.0
        endif
        if (CS%vmask(I,J) == 1) then
          err_tempv = MAX(ABS(v_last(I,J)-v_shlf(I,J)), err_tempu)
          if (err_tempv >= err_max) err_max = err_tempv
          tempv = SQRT((v_shlf(I,J)**2) + (tempu**2))
        endif
        if (tempv >= max_vel) max_vel = tempv
      enddo ; enddo

      u_last(:,:) = u_shlf(:,:)
      v_last(:,:) = v_shlf(:,:)

      call max_across_PEs(max_vel)
      call max_across_PEs(err_max)
      err_init = max_vel

    elseif (CS%nonlin_solve_err_mode == 3) then
      PrevNorm=Norm; Norm=0.0; Normvec=0.0
      do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
        if (CS%umask(I,J) == 1) Normvec(I,J) = Normvec(I,J) + (u_shlf(I,J)**2 * US%L_T_to_m_s**2)
        if (CS%vmask(I,J) == 1) Normvec(I,J) = Normvec(I,J) + (v_shlf(I,J)**2 * US%L_T_to_m_s**2)
      enddo; enddo
      Norm = reproducing_sum( Normvec, Is_sum, Ie_sum, Js_sum, Je_sum )
      Norm = sqrt(Norm)
      err_max=2.*abs(Norm-PrevNorm); err_init=Norm+PrevNorm
    endif

    write(mesg,*) "ice_shelf_solve_outer: nonlinear fractional residual = ", err_max/err_init
    call MOM_mesg(mesg, 5)

    if (err_max <= CS%nonlinear_tolerance * err_init) then
      exit
    endif

  enddo

  write(mesg,*) "ice_shelf_solve_outer: nonlinear fractional residual = ", err_max/err_init
  call MOM_mesg(mesg)
  write(mesg,*) "ice_shelf_solve_outer: exiting nonlinear solve after ",iter," iterations"
  call MOM_mesg(mesg)

end subroutine ice_shelf_solve_outer

subroutine ice_shelf_solve_inner(CS, ISS, G, US, u_shlf, v_shlf, taudx, taudy, H_node, float_cond, &
                                 hmask, conv_flag, iters, time, Phi, Phisub)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                           !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_shlf  !< The zonal ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_shlf  !< The meridional ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: taudx !< The x-direction driving stress [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: taudy  !< The y-direction driving stress [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal (corner)
                                             !! points [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: float_cond !< If GL_regularize=true, indicates cells containing
                                                !! the grounding line (float_cond=1) or not (float_cond=0)
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  integer,                intent(out)   :: conv_flag !< A flag indicating whether (1) or not (0) the
                                           !! iterations have converged to the specified tolerance
  integer,                intent(out)   :: iters !< The number of iterations used in the solver.
  type(time_type),        intent(in)    :: Time !< The current model time
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: Phi !< The gradients of bilinear basis elements at Gaussian
                                             !! quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real, dimension(:,:,:,:,:,:), &
                          intent(in)    :: Phisub !< Quadrature structure weights at subgridscale
                                            !! locations for finite element calculations [nondim]
! one linear solve (nonlinear iteration) of the solution for velocity

! in this subroutine:
!    RHS = taud
!    diagonal of matrix is found (for Jacobi precondition)
!    CG iteration is carried out for max. iterations or until convergence

! assumed - u, v, taud, visc, basal_traction are valid on the halo

  real, dimension(SZDIB_(G),SZDJB_(G)) ::  &
                        Ru, Rv, &     ! Residuals in the stress calculations [R L3 Z T-2 ~> m kg s-2]
                        Ru_old, Rv_old, & ! Previous values of Ru and Rv [R L3 Z T-2 ~> m kg s-2]
                        Zu, Zv, & ! Contributions to velocity changes [L T-1 ~> m s-1]
                        Zu_old, Zv_old, & ! Previous values of Zu and Zv [L T-1 ~> m s-1]
                        DIAGu, DIAGv, & ! Diagonals with units like Ru/Zu [R L2 Z T-1 ~> kg s-1]
                        RHSu, RHSv, & ! Right hand side of the stress balance [R L3 Z T-2 ~> m kg s-2]
                        Au, Av, & ! The retarding lateral stress contributions [R L3 Z T-2 ~> kg m s-2]
                        Du, Dv, & ! Velocity changes [L T-1 ~> m s-1]
                        sum_vec, sum_vec_2, sum_vec_3 !, &
                        !ubd, vbd   ! Boundary stress contributions [R L3 Z T-2 ~> kg m s-2]
  real    :: beta_k, dot_p1, resid0tol2, cg_halo, max_cg_halo
  real    :: alpha_k     ! A scaling factor for iterative corrections [nondim]
  real    :: resid_scale ! A scaling factor for redimensionalizing the global residuals [m2 L-2 ~> 1]
                         ! [m2 L-2 ~> 1] [R L3 Z T-2 ~> m kg s-2]
  real    :: resid2_scale ! A scaling factor for redimensionalizing the global squared residuals
                         ! [m2 L-2 ~> 1] [R L3 Z T-2 ~> m kg s-2]
  real    :: rhoi_rhow  ! The density of ice divided by a typical water density [nondim]
  integer :: iter, i, j, isd, ied, jsd, jed, isc, iec, jsc, jec, is, js, ie, je
  integer :: Is_sum, Js_sum, Ie_sum, Je_sum ! Loop bounds for global sums or arrays starting at 1.
  integer :: Isdq, Iedq, Jsdq, Jedq, Iscq, Iecq, Jscq, Jecq, nx_halo, ny_halo
  integer :: Iscq_sv, Jscq_sv ! Starting loop bound for sum_vec

  Isdq = G%IsdB ; Iedq = G%IedB ; Jsdq = G%JsdB ; Jedq = G%JedB
  Iscq = G%IscB ; Iecq = G%IecB ; Jscq = G%JscB ; Jecq = G%JecB
  ny_halo = G%domain%njhalo ; nx_halo = G%domain%nihalo
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg

  Zu(:,:) = 0 ; Zv(:,:) = 0 ; DIAGu(:,:) = 0 ; DIAGv(:,:) = 0
  Ru(:,:) = 0 ; Rv(:,:) = 0 ; Au(:,:) = 0 ; Av(:,:) = 0 ; RHSu(:,:) = 0 ; RHSv(:,:) = 0
  Du(:,:) = 0 ; Dv(:,:) = 0
  dot_p1 = 0

  ! Determine the loop limits for sums, bearing in mind that the arrays will be starting at 1.
  ! Includes the edge of the tile is at the western/southern bdry (if symmetric)
  if ((isc+G%idg_offset==G%isg) .and. (.not. CS%reentrant_x)) then
    Is_sum = Iscq + (1-Isdq) ; Iscq_sv = Iscq
  else
    Is_sum = isc  + (1-Isdq) ; Iscq_sv = isc
  endif
  if ((jsc+G%jdg_offset==G%jsg) .and. (.not. CS%reentrant_y)) then
    Js_sum = Jscq + (1-Jsdq) ; Jscq_sv = Jscq
  else
    Js_sum = jsc + (1-Jsdq) ; Jscq_sv = jsc
  endif
  Ie_sum = Iecq + (1-Isdq) ; Je_sum = Jecq + (1-Jsdq)

  RHSu(:,:) = taudx(:,:) ; RHSv(:,:) = taudy(:,:)

  call pass_vector(RHSu, RHSv, G%domain, TO_ALL, BGRID_NE, complete=.false.)

  call matrix_diagonal(CS, G, US, float_cond, H_node, CS%ice_visc, CS%basal_traction, &
                       hmask, rhoi_rhow, Phi, Phisub, DIAGu, DIAGv)

  call pass_vector(DIAGu, DIAGv, G%domain, TO_ALL, BGRID_NE, complete=.false.)

  call CG_action(CS, Au, Av, u_shlf, v_shlf, Phi, Phisub, CS%umask, CS%vmask, hmask, &
                 H_node, CS%ice_visc, float_cond, CS%bed_elev, CS%basal_traction, &
                 G, US, isc-1, iec+1, jsc-1, jec+1, rhoi_rhow)

  call pass_vector(Au, Av, G%domain, TO_ALL, BGRID_NE, complete=.true.)

  Ru(:,:) = (RHSu(:,:) - Au(:,:)) ; Rv(:,:) = (RHSv(:,:) - Av(:,:))
  resid_scale = (US%L_to_m**2*US%s_to_T)*(US%RZ_to_kg_m2*US%L_T_to_m_s**2)
  resid2_scale = ((US%RZ_to_kg_m2*US%L_to_m)*US%L_T_to_m_s**2)**2

  sum_vec(:,:) = 0.0
  do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
    if (CS%umask(I,J) == 1) sum_vec(I,J) = resid2_scale*Ru(I,J)**2
    if (CS%vmask(I,J) == 1) sum_vec(I,J) = sum_vec(I,J) + resid2_scale*Rv(I,J)**2
  enddo ; enddo

  !resid0 = sqrt(reproducing_sum( sum_vec, Is_sum, Ie_sum, Js_sum, Je_sum ))
  resid0tol2 = CS%cg_tolerance**2 * reproducing_sum( sum_vec, Is_sum, Ie_sum, Js_sum, Je_sum )

  do J=Jsdq,Jedq ; do I=Isdq,Iedq
    if (CS%umask(I,J) == 1 .AND.(DIAGu(I,J)/=0)) Zu(I,J) = Ru(I,J) / DIAGu(I,J)
    if (CS%vmask(I,J) == 1 .AND.(DIAGv(I,J)/=0)) Zv(I,J) = Rv(I,J) / DIAGv(I,J)
  enddo ; enddo

  Du(:,:) = Zu(:,:) ; Dv(:,:) = Zv(:,:)

  if (G%symmetric) then
    max_cg_halo=min(nx_halo,ny_halo)
  else
    max_cg_halo=min(nx_halo,ny_halo)-1
  endif
  cg_halo = max_cg_halo
  conv_flag = 0

  !!!!!!!!!!!!!!!!!!
  !!              !!
  !! MAIN CG LOOP !!
  !!              !!
  !!!!!!!!!!!!!!!!!!

  ! initially, c-grid data is valid up to 3 halo nodes out

  do iter = 1,CS%cg_max_iterations

    ! we can never assume that any arrays are legit more than 3 vertices past
    ! the computational domain - this is their state in the initial iteration

    is = isc - cg_halo ; ie = Iecq + cg_halo
    js = jsc - cg_halo ; je = Jecq + cg_halo

    Au(:,:) = 0 ; Av(:,:) = 0

    call CG_action(CS, Au, Av, Du, Dv, Phi, Phisub, CS%umask, CS%vmask, hmask, &
                   H_node, CS%ice_visc, float_cond, CS%bed_elev, CS%basal_traction, &
                   G, US, is, ie, js, je, rhoi_rhow)

    ! Au, Av valid region moves in by 1

    call pass_vector(Au,Av,G%domain, TO_ALL, BGRID_NE)

    sum_vec(:,:) = 0.0 ; sum_vec_2(:,:) = 0.0

    do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
      if (CS%umask(I,J) == 1) then
        sum_vec(I,J)   = resid_scale * (Zu(I,J) * Ru(I,J))
        sum_vec_2(I,J) = resid_scale * (Du(I,J) * Au(I,J))
        Ru_old(I,J) = Ru(I,J) ; Zu_old(I,J) = Zu(I,J)
      endif
      if (CS%vmask(I,J) == 1) then
        sum_vec(I,J)   = sum_vec(I,J)   + resid_scale * (Zv(I,J) * Rv(I,J))
        sum_vec_2(I,J) = sum_vec_2(I,J) + resid_scale * (Dv(I,J) * Av(I,J))
        Rv_old(I,J) = Rv(I,J) ; Zv_old(I,J) = Zv(I,J)
      endif
    enddo ; enddo

    alpha_k = reproducing_sum( sum_vec, Is_sum, Ie_sum, Js_sum, Je_sum ) / &
              reproducing_sum( sum_vec_2, Is_sum, Ie_sum, Js_sum, Je_sum )

    do J=js,je-1 ; do I=is,ie-1
      if (CS%umask(I,J) == 1) then
        u_shlf(I,J) = u_shlf(I,J) + alpha_k * Du(I,J)
        Ru(I,J) = Ru(I,J) - alpha_k * Au(I,J)
        if (DIAGu(I,J)/=0) Zu(I,J) = Ru(I,J) / DIAGu(I,J)
      endif
      if (CS%vmask(I,J) == 1) then
        v_shlf(I,J) = v_shlf(I,J) + alpha_k * Dv(I,J)
        Rv(I,J) = Rv(I,J) - alpha_k * Av(I,J)
        if (DIAGv(I,J)/=0) Zv(I,J) = Rv(I,J) / DIAGv(I,J)
      endif
    enddo; enddo


    ! R,u,v,Z valid region moves in by 1

    ! beta_k = (Z \dot R) / (Zold \dot Rold)
    sum_vec(:,:) = 0.0 ; sum_vec_2(:,:) = 0.0 ; sum_vec_3(:,:) = 0.0

    do J=jscq_sv,jecq ; do i=iscq_sv,iecq
      if (CS%umask(I,J) == 1) then
        sum_vec(I,J)   = resid_scale  * (Zu(I,J) * Ru(I,J))
        sum_vec_2(I,J) = resid_scale  * (Zu_old(I,J) * Ru_old(I,J))
        sum_vec_3(I,J) = resid2_scale * Ru(I,J)**2
      endif
      if (CS%vmask(I,J) == 1) then
        sum_vec(I,J)   = sum_vec(I,J)   + resid_scale  * (Zv(I,J) * Rv(I,J))
        sum_vec_2(I,J) = sum_vec_2(I,J) + resid_scale  * (Zv_old(I,J) * Rv_old(I,J))
        sum_vec_3(I,J) = sum_vec_3(I,J) + resid2_scale * Rv(I,J)**2
      endif
    enddo ; enddo

    beta_k = reproducing_sum(sum_vec, Is_sum, Ie_sum, Js_sum, Je_sum ) / &
             reproducing_sum(sum_vec_2, Is_sum, Ie_sum, Js_sum, Je_sum )

    do J=js,je-1 ; do I=is,ie-1
        if (CS%umask(I,J) == 1) Du(I,J) = Zu(I,J) + beta_k * Du(I,J)
        if (CS%vmask(I,J) == 1) Dv(I,J) = Zv(I,J) + beta_k * Dv(I,J)
    enddo ; enddo

   ! D valid region moves in by 1

    dot_p1 = reproducing_sum( sum_vec_3, Is_sum, Ie_sum, Js_sum, Je_sum )

    !if sqrt(dot_p1) <= (CS%cg_tolerance * resid0)
    if (dot_p1 <= resid0tol2) then
      iters = iter
      conv_flag = 1
      exit
    endif

    cg_halo = cg_halo - 1

    if (cg_halo == 0) then
     ! pass vectors
      call pass_vector(Du, Dv, G%domain, TO_ALL, BGRID_NE, complete=.false.)
      call pass_vector(u_shlf, v_shlf, G%domain, TO_ALL, BGRID_NE, complete=.false.)
      call pass_vector(Ru, Rv, G%domain, TO_ALL, BGRID_NE, complete=.true.)
      cg_halo = max_cg_halo
    endif

  enddo ! end of CG loop

  do J=Jsdq,Jedq ; do I=Isdq,Iedq
      if (CS%umask(I,J) == 3) then
        u_shlf(I,J) = CS%u_bdry_val(I,J)
      elseif (CS%umask(I,J) == 0) then
        u_shlf(I,J) = 0
      endif

      if (CS%vmask(I,J) == 3) then
        v_shlf(I,J) = CS%v_bdry_val(I,J)
      elseif (CS%vmask(I,J) == 0) then
        v_shlf(I,J) = 0
      endif
  enddo ; enddo

  call pass_vector(u_shlf, v_shlf, G%domain, TO_ALL, BGRID_NE)
  if (conv_flag == 0) then
    iters = CS%cg_max_iterations
  endif

end subroutine ice_shelf_solve_inner

subroutine ice_shelf_advect_thickness_x(CS, G, LB, time_step, hmask, h0, h_after_uflux, uh_ice)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(loop_bounds_type), intent(in)    :: LB   !< Loop bounds structure.
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h0 !< The initial ice shelf thicknesses [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_uflux !< The ice shelf thicknesses after
                                              !! the zonal mass fluxes [Z ~> m].
  real, dimension(SZDIB_(G),SZDJ_(G)), &
                          intent(inout) :: uh_ice !< The accumulated zonal ice volume flux [Z L2 ~> m3]

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization


  integer :: i, j
  integer :: ish, ieh, jsh, jeh
  real :: u_face     ! Zonal velocity at a face [L T-1 ~> m s-1]
  real :: h_face     ! Thickness at a face for transport [Z ~> m]
  real :: slope_lim  ! The value of the slope limiter, in the range of 0 to 2 [nondim]

!  is = G%isc-2 ; ie = G%iec+2 ; js = G%jsc ; je = G%jec
!  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh

  ! hmask coded values: 1) fully covered; 2) partly covered - no export; 3) Specified boundary condition
  ! relevant u_face_mask coded values: 1) Normal interior point; 4) Specified flux BC

  do j=jsh,jeh ; do I=ish-1,ieh
    if (CS%u_face_mask(I,j) == 4.) then ! The flux itself is a specified boundary condition.
      uh_ice(I,j) = (time_step * G%dyCu(I,j)) * CS%u_flux_bdry_val(I,j)
    elseif ((hmask(i,j) == 1 .or. hmask(i,j) == 3) .or. (hmask(i+1,j) == 1 .or. hmask(i+1,j) == 3)) then
      u_face = 0.5 * (CS%u_shelf(I,J-1) + CS%u_shelf(I,J))
      h_face = 0.0 ! This will apply when the source cell is iceless or not fully ice covered.

      if (u_face > 0) then
        if (hmask(i,j) == 3) then ! This is a open boundary inflow from the west
          h_face = CS%h_bdry_val(i,j)
        elseif (hmask(i,j) == 1) then ! There can be eastward flow through this face.
          if ((hmask(i-1,j) == 1 .or. hmask(i-1,j) == 3) .and. &
            (hmask(i+1,j) == 1 .or. hmask(i+1,j) == 3)) then
            slope_lim = slope_limiter(h0(i,j)-h0(i-1,j), h0(i+1,j)-h0(i,j))
            ! This is a 2nd-order centered scheme with a slope limiter.  We could try PPM here.
            h_face = h0(i,j) - slope_lim * (0.5 * (h0(i,j)-h0(i+1,j)))
          else
            h_face = h0(i,j)
          endif
        endif
      else
        if (hmask(i+1,j) == 3) then ! This is a open boundary inflow from the east
          h_face = CS%h_bdry_val(i+1,j)
        elseif (hmask(i+1,j) == 1) then
          if ((hmask(i,j) == 1 .or. hmask(i,j) == 3) .and. &
            (hmask(i+2,j) == 1 .or. hmask(i+2,j) == 3)) then
            slope_lim = slope_limiter(h0(i+1,j)-h0(i,j), h0(i+2,j)-h0(i+1,j))
            h_face = h0(i+1,j) - slope_lim * (0.5 * (h0(i+2,j)-h0(i+1,j)))
          else
            h_face = h0(i+1,j)
          endif
        endif
      endif

      uh_ice(I,j) = (time_step * G%dyCu(I,j)) * (u_face * h_face)
    else
      uh_ice(I,j) = 0.0
    endif
  enddo ; enddo

  do j=jsh,jeh ; do i=ish,ieh
    if (hmask(i,j) /= 3) &
      h_after_uflux(i,j) = h0(i,j) + (uh_ice(I-1,j) - uh_ice(I,j)) * G%IareaT(i,j)

     ! Update the masks of cells that have gone from no ice to partial ice.
    if ((hmask(i,j) == 0) .and. ((uh_ice(I-1,j) > 0.0) .or. (uh_ice(I,j) < 0.0))) hmask(i,j) = 2
  enddo ; enddo

end subroutine ice_shelf_advect_thickness_x

subroutine ice_shelf_advect_thickness_y(CS, G, LB, time_step, hmask, h0, h_after_vflux, vh_ice)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(loop_bounds_type), intent(in)    :: LB !< Loop bounds structure.
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: hmask !< A mask indicating which tracer points are
                                              !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h0 !< The initial ice shelf thicknesses [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_vflux !< The ice shelf thicknesses after
                                              !! the meridional mass fluxes [Z ~> m].
  real, dimension(SZDI_(G),SZDJB_(G)), &
                          intent(inout) :: vh_ice !< The accumulated meridional ice volume flux [Z L2 ~> m3]

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization


  integer :: i, j
  integer :: ish, ieh, jsh, jeh
  real :: v_face     ! Pseudo-meridional velocity at a face [L T-1 ~> m s-1]
  real :: h_face     ! Thickness at a face for transport [Z ~> m]
  real :: slope_lim  ! The value of the slope limiter, in the range of 0 to 2 [nondim]

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh

  ! hmask coded values: 1) fully covered; 2) partly covered - no export; 3) Specified boundary condition
  ! relevant u_face_mask coded values: 1) Normal interior point; 4) Specified flux BC

  do J=jsh-1,jeh ; do i=ish,ieh
    if (CS%v_face_mask(i,J) == 4.) then ! The flux itself is a specified boundary condition.
      vh_ice(i,J) = (time_step * G%dxCv(i,J)) * CS%v_flux_bdry_val(i,J)
    elseif ((hmask(i,j) == 1 .or. hmask(i,j) == 3) .or. (hmask(i,j+1) == 1 .or. hmask(i,j+1) == 3)) then
      v_face = 0.5 * (CS%v_shelf(I-1,J) + CS%v_shelf(I,J))
      h_face = 0.0 ! This will apply when the source cell is iceless or not fully ice covered.

      if (v_face > 0) then
        if (hmask(i,j) == 3) then ! This is a open boundary inflow from the south
          h_face = CS%h_bdry_val(i,j)
        elseif (hmask(i,j) == 1) then ! There can be northward flow through this face.
          if ((hmask(i,j-1) == 1 .or. hmask(i,j-1) == 3) .and. &
            (hmask(i,j+1) == 1 .or. hmask(i,j+1) == 3)) then
            slope_lim = slope_limiter(h0(i,j)-h0(i,j-1), h0(i,j+1)-h0(i,j))
            ! This is a 2nd-order centered scheme with a slope limiter.  We could try PPM here.
            h_face = h0(i,j) - slope_lim * (0.5 * (h0(i,j)-h0(i,j+1)))
          else
            h_face = h0(i,j)
          endif
        endif
      else
        if (hmask(i,j+1) == 3) then ! This is a open boundary inflow from the north
          h_face = CS%h_bdry_val(i,j+1)
        elseif (hmask(i,j+1) == 1) then
          if ((hmask(i,j) == 1 .or. hmask(i,j) == 3) .and. &
            (hmask(i,j+2) == 1 .or. hmask(i,j+2) == 3)) then
            slope_lim = slope_limiter(h0(i,j+1)-h0(i,j), h0(i,j+2)-h0(i,j+1))
            h_face = h0(i,j+1) - slope_lim * (0.5 * (h0(i,j+2)-h0(i,j+1)))
          else
            h_face = h0(i,j+1)
          endif
        endif
      endif

      vh_ice(i,J) = (time_step * G%dxCv(i,J)) * (v_face * h_face)
    else
      vh_ice(i,J) = 0.0
    endif
  enddo ; enddo

  do j=jsh,jeh ; do i=ish,ieh
    if (hmask(i,j) /= 3) &
      h_after_vflux(i,j) = h0(i,j) + (vh_ice(i,J-1) - vh_ice(i,J)) * G%IareaT(i,j)

    ! Update the masks of cells that have gone from no ice to partial ice.
    if ((hmask(i,j) == 0) .and. ((vh_ice(i,J-1) > 0.0) .or. (vh_ice(i,J) < 0.0))) hmask(i,j) = 2
  enddo ; enddo

end subroutine ice_shelf_advect_thickness_y

subroutine shelf_advance_front(CS, ISS, G, hmask, uh_ice, vh_ice)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                           !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: hmask !< A mask indicating which tracer points are
                                              !! partly or fully covered by an ice-shelf
  real, dimension(SZDIB_(G),SZDJ_(G)), &
                          intent(inout) :: uh_ice !< The accumulated zonal ice volume flux [Z L2 ~> m3]
  real, dimension(SZDI_(G),SZDJB_(G)), &
                          intent(inout) :: vh_ice !< The accumulated meridional ice volume flux [Z L2 ~> m3]

  ! in this subroutine we go through the computational cells only and, if they are empty or partial cells,
  ! we find the reference thickness and update the shelf mass and partial area fraction and the hmask if necessary

  ! if any cells go from partial to complete, we then must set the thickness, update hmask accordingly,
  ! and divide the overflow across the adjacent EMPTY (not partly-covered) cells.
  ! (it is highly unlikely there will not be any; in which case this will need to be rethought.)

  ! most likely there will only be one "overflow". If not, though, a pass_var of all relevant variables
  ! is done; there will therefore be a loop which, in practice, will hopefully not have to go through
  ! many iterations

  ! when 3d advected scalars are introduced, they will be impacted by what is done here

  ! flux_enter(isd:ied,jsd:jed,1:4): if cell is not ice-covered, gives flux of ice into cell from kth boundary
  !
  !   from eastern neighbor:  flux_enter(:,:,1)
  !   from western neighbor:  flux_enter(:,:,2)
  !   from southern neighbor: flux_enter(:,:,3)
  !   from northern neighbor: flux_enter(:,:,4)
  !
  !        o--- (4) ---o
  !        |           |
  !       (1)         (2)
  !        |           |
  !        o--- (3) ---o
  !

  integer :: i, j, isc, iec, jsc, jec, n_flux, k, iter_count
  integer :: i_off, j_off
  integer :: iter_flag

  real :: h_reference ! A reference thicknesss based on neighboring cells [Z ~> m]
  real :: h_reference_ew !contribution to reference thickness from east + west cells [Z ~> m]
  real :: h_reference_ns !contribution to reference thickness from north + south cells [Z ~> m]
  real :: tot_flux    ! The total ice mass flux [Z L2 ~> m3]
  real :: tot_flux_ew ! The contribution to total ice mass flux from east + west cells [Z L2 ~> m3]
  real :: tot_flux_ns ! The contribution to total ice mass flux from north + south cells [Z L2 ~> m3]
  real :: partial_vol ! The volume covered by ice shelf [Z L2 ~> m3]
  real :: dxdyh       ! Cell area [L2 ~> m2]
  character(len=160) :: mesg  ! The text of an error message
  integer, dimension(4) :: mapi, mapj, new_partial
  real, dimension(SZDI_(G),SZDJ_(G),4) :: flux_enter  ! The ice volume flux into the
                                              ! cell through the 4 cell boundaries [Z L2 ~> m3].
  real, dimension(SZDI_(G),SZDJ_(G),4) :: flux_enter_replace ! An updated ice volume flux into the
                                              ! cell through the 4 cell boundaries [Z L2 ~> m3].

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  i_off = G%idg_offset ; j_off = G%jdg_offset
  iter_count = 0 ; iter_flag = 1

  flux_enter(:,:,:) = 0.0
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    if ((hmask(i,j) == 0) .or. (hmask(i,j) == 2)) then
      flux_enter(i,j,1) = max(uh_ice(I-1,j), 0.0)
      flux_enter(i,j,2) = max(-uh_ice(I,j), 0.0)
      flux_enter(i,j,3) = max(vh_ice(i,J-1), 0.0)
      flux_enter(i,j,4) = max(-vh_ice(i,J), 0.0)
    endif
  enddo ; enddo

  mapi(1) = -1 ; mapi(2) = 1 ; mapi(3:4) = 0
  mapj(3) = -1 ; mapj(4) = 1 ; mapj(1:2) = 0

  do while (iter_flag == 1)

    iter_flag = 0

    if (iter_count > 0) then
      flux_enter(:,:,:) = flux_enter_replace(:,:,:)
    endif
    flux_enter_replace(:,:,:) = 0.0

    iter_count = iter_count + 1

    ! if iter_count >= 3 then some halo updates need to be done...
    if (iter_count==3) then
      call MOM_error(FATAL, "MOM_ice_shelf_dyn.F90, shelf_advance_front iter >=3.")
    endif

    do j=jsc-1,jec+1

      if (CS%reentrant_y .OR. (((j+j_off) <= G%domain%njglobal) .AND. &
          ((j+j_off) >= 1))) then

        do i=isc-1,iec+1

          if (CS%reentrant_x .OR. (((i+i_off) <= G%domain%niglobal) .AND. &
              ((i+i_off) >= 1))) then
            ! first get reference thickness by averaging over cells that are fluxing into this cell
            n_flux = 0
            h_reference_ew = 0.0
            h_reference_ns = 0.0
            tot_flux_ew = 0.0
            tot_flux_ns = 0.0

            do k=1,2
              if (flux_enter(i,j,k) > 0) then
                n_flux = n_flux + 1
                h_reference_ew = h_reference_ew + flux_enter(i,j,k) * ISS%h_shelf(i+2*k-3,j)
                !h_reference = h_reference + ISS%h_shelf(i+2*k-3,j)
                tot_flux_ew = tot_flux_ew + flux_enter(i,j,k)
                flux_enter(i,j,k) = 0.0
              endif
            enddo

            do k=1,2
              if (flux_enter(i,j,k+2) > 0) then
                n_flux = n_flux + 1
                h_reference_ns = h_reference_ns + flux_enter(i,j,k+2) * ISS%h_shelf(i,j+2*k-3)
                !h_reference = h_reference + ISS%h_shelf(i,j+2*k-3)
                tot_flux_ns = tot_flux_ns + flux_enter(i,j,k+2)
                flux_enter(i,j,k+2) = 0.0
              endif
            enddo

            h_reference = h_reference_ew + h_reference_ns
            tot_flux = tot_flux_ew + tot_flux_ns

            if (n_flux > 0) then
              dxdyh = G%areaT(i,j)
              h_reference = h_reference / tot_flux
              !h_reference = h_reference / real(n_flux)
              partial_vol = ISS%h_shelf(i,j) * ISS%area_shelf_h(i,j) + tot_flux

              if ((partial_vol / G%areaT(i,j)) == h_reference) then ! cell is exactly covered, no overflow
                if (ISS%hmask(i,j)/=3) ISS%hmask(i,j) = 1
                ISS%h_shelf(i,j) = h_reference
                ISS%area_shelf_h(i,j) = G%areaT(i,j)
              elseif ((partial_vol / G%areaT(i,j)) < h_reference) then
                ISS%hmask(i,j) = 2
               !  ISS%mass_shelf(i,j) = partial_vol * CS%density_ice
                ISS%area_shelf_h(i,j) = partial_vol / h_reference
                ISS%h_shelf(i,j) = h_reference
              else

                if (ISS%hmask(i,j)/=3) ISS%hmask(i,j) = 1
                ISS%area_shelf_h(i,j) = G%areaT(i,j)
                !h_temp(i,j) = h_reference
                partial_vol = partial_vol - h_reference * G%areaT(i,j)

                iter_flag  = 1

                n_flux = 0 ; new_partial(:) = 0

                do k=1,2
                  if (CS%u_face_mask(I-2+k,j) == 2) then
                    n_flux = n_flux + 1
                  elseif (ISS%hmask(i+2*k-3,j) == 0) then
                    n_flux = n_flux + 1
                    new_partial(k) = 1
                  endif
                  if (CS%v_face_mask(i,J-2+k) == 2) then
                    n_flux = n_flux + 1
                  elseif (ISS%hmask(i,j+2*k-3) == 0) then
                    n_flux = n_flux + 1
                    new_partial(k+2) = 1
                  endif
                enddo

                if (n_flux == 0) then ! there is nowhere to put the extra ice!
                  ISS%h_shelf(i,j) = h_reference + partial_vol / G%areaT(i,j)
                else
                  ISS%h_shelf(i,j) = h_reference

                  do k=1,2
                    if (new_partial(k) == 1) &
                      flux_enter_replace(i+2*k-3,j,3-k) = partial_vol / real(n_flux)
                    if (new_partial(k+2) == 1) &
                      flux_enter_replace(i,j+2*k-3,5-k) = partial_vol / real(n_flux)
                  enddo
                endif

              endif ! Parital_vol test.
            endif ! n_flux gt 0 test.

          endif
        enddo ! j-loop
      endif
    enddo

  !  call max_across_PEs(iter_flag)

  enddo ! End of do while(iter_flag) loop

  call max_across_PEs(iter_count)

  if (is_root_pe() .and. (iter_count > 1)) then
    write(mesg,*) "shelf_advance_front: ", iter_count, " max iterations"
    call MOM_mesg(mesg, 5)
  endif

end subroutine shelf_advance_front

!> Apply a very simple calving law using a minimum thickness rule
subroutine ice_shelf_min_thickness_calve(G, h_shelf, area_shelf_h, hmask, thickness_calve, halo)
  type(ocean_grid_type), intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: area_shelf_h !< The area per cell covered by
                                             !! the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real,                  intent(in)    :: thickness_calve !< The thickness at which to trigger calving [Z ~> m].
  integer,     optional, intent(in)    :: halo  !< The number of halo points to use.  If not present,
                                                !! work on the entire data domain.
  integer :: i, j, is, ie, js, je

  if (present(halo)) then
    is = G%isc - halo ; ie = G%iec + halo ; js = G%jsc - halo ; je = G%jec + halo
  else
    is = G%isd ; ie = G%ied ; js = G%jsd ; je = G%jed
  endif

  do j=js,je ; do i=is,ie
!    if ((h_shelf(i,j) < CS%thickness_calve) .and. (hmask(i,j) == 1) .and. &
!        (CS%ground_frac(i,j) == 0.0)) then
    if ((h_shelf(i,j) < thickness_calve) .and. (area_shelf_h(i,j) > 0.)) then
      h_shelf(i,j) = 0.0
      area_shelf_h(i,j) = 0.0
      hmask(i,j) = 0.0
    endif
  enddo ; enddo

end subroutine ice_shelf_min_thickness_calve

subroutine calve_to_mask(G, h_shelf, area_shelf_h, hmask, calve_mask)
  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: area_shelf_h !< The area per cell covered by
                                                             !! the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: hmask !< A mask indicating which tracer points are
                                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), intent(in)    :: calve_mask !< A mask that indicates where the ice
                                                             !! shelf can exist, and where it will calve.

  integer                        :: i,j

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((calve_mask(i,j) == 0.0) .and. (hmask(i,j) /= 0.0)) then
      h_shelf(i,j) = 0.0
      area_shelf_h(i,j) = 0.0
      hmask(i,j) = 0.0
    endif
  enddo ; enddo

end subroutine calve_to_mask

!> Calculate driving stress using cell-centered bed elevation and ice thickness
subroutine calc_shelf_driving_stress(CS, ISS, G, US, taudx, taudy, OD)
  type(ice_shelf_dyn_CS), intent(in)   :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state), intent(in)    :: ISS !< A structure with elements that describe
                                             !! the ice-shelf state
  type(ocean_grid_type), intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: OD  !< ocean floor depth at tracer points [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: taudx  !< X-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: taudy  !< Y-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]


! driving stress!

! ! taudx and taudy will hold driving stress in the x- and y- directions when done.
!    they will sit on the BGrid, and so their size depends on whether the grid is symmetric
!
! Since this is a finite element solve, they will actually have the form \int \Phi_i rho g h \nabla s
!
! OD -this is important and we do not yet know where (in MOM) it will come from. It represents
!     "average" ocean depth -- and is needed to find surface elevation
!    (it is assumed that base_ice = bed + OD)

  real, dimension(SIZE(OD,1),SIZE(OD,2))  :: S     ! surface elevation [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)) :: sx_e, sy_e !element contributions to driving stress
  real    :: rho, rhow, rhoi_rhow ! Ice and ocean densities [R ~> kg m-3]
  real    :: sx, sy    ! Ice shelf top slopes [Z L-1 ~> nondim]
  real    :: neumann_val ! [R Z L2 T-2 ~> kg s-2]
  real    :: dxh, dyh,Dx,Dy  ! Local grid spacing [L ~> m]
  real    :: grav      ! The gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real    :: scale     ! Scaling factor used to ensure surface slope magnitude does not exceed CS%max_surface_slope
  integer :: i, j, iscq, iecq, jscq, jecq, isd, jsd, ied, jed, is, js, iegq, jegq
  integer :: giec, gjec, gisc, gjsc, cnt, isc, jsc, iec, jec
  integer :: i_off, j_off

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
!  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
!  iegq = G%iegB ; jegq = G%jegB
!  gisc = G%domain%nihalo+1 ; gjsc = G%domain%njhalo+1
  gisc = 1 ; gjsc = 1
!  giec = G%domain%niglobal+G%domain%nihalo ; gjec = G%domain%njglobal+G%domain%njhalo
  giec = G%domain%niglobal ; gjec = G%domain%njglobal
!  is = iscq - 1; js = jscq - 1
  i_off = G%idg_offset ; j_off = G%jdg_offset


  rho =  CS%density_ice
  rhow = CS%density_ocean_avg
  grav = CS%g_Earth
  rhoi_rhow = rho/rhow
  ! prelim - go through and calculate S

  if (CS%GL_couple) then
    do j=jsc-G%domain%njhalo,jec+G%domain%njhalo
      do i=isc-G%domain%nihalo,iec+G%domain%nihalo
        S(i,j) = -CS%bed_elev(i,j) + (OD(i,j) + max(ISS%h_shelf(i,j),CS%min_h_shelf))
      enddo
    enddo
  else
    ! check whether the ice is floating or grounded
    do j=jsc-G%domain%njhalo,jec+G%domain%njhalo
      do i=isc-G%domain%nihalo,iec+G%domain%nihalo
        if (rhoi_rhow * max(ISS%h_shelf(i,j),CS%min_h_shelf) - CS%bed_elev(i,j) <= 0) then
          S(i,j) = (1 - rhoi_rhow)*max(ISS%h_shelf(i,j),CS%min_h_shelf)
        else
          S(i,j) = max(ISS%h_shelf(i,j),CS%min_h_shelf)-CS%bed_elev(i,j)
        endif
      enddo
    enddo
  endif

  call pass_var(S, G%domain)

  sx_e(:,:)=0.0; sy_e(:,:)=0.0

  do j=jsc-1,jec+1
    do i=isc-1,iec+1
      cnt = 0
      sx = 0
      sy = 0
      dxh = G%dxT(i,j)
      dyh = G%dyT(i,j)
      Dx=dxh
      Dy=dyh
      if (ISS%hmask(i,j) == 1 .or. ISS%hmask(i,j) == 3) then
        ! we are inside the global computational bdry, at an ice-filled cell

        ! calculate sx
        if (((i+i_off) == gisc) .and. (.not. CS%reentrant_x)) then ! at west computational bdry
         if (ISS%hmask(i+1,j) == 1 .or. ISS%hmask(i+1,j) == 3) then
            sx = (S(i+1,j)-S(i,j))/dxh
          else
            sx = 0
          endif
        elseif (((i+i_off) == giec) .and. (.not. CS%reentrant_x)) then ! at east computational bdry
          if (ISS%hmask(i-1,j) == 1 .or. ISS%hmask(i-1,j) == 3) then
            sx = (S(i,j)-S(i-1,j))/dxh
          else
            sx = 0
          endif
        else ! interior
          if (ISS%hmask(i+1,j) == 1 .or. ISS%hmask(i+1,j) == 3) then
            cnt = cnt+1
            Dx = dxh + G%dxT(i+1,j)
            sx = S(i+1,j)
          else
            sx = S(i,j)
          endif
          if (ISS%hmask(i-1,j) == 1 .or. ISS%hmask(i-1,j) == 3) then
            cnt = cnt+1
            Dx = dxh + G%dxT(i-1,j)
            sx = sx - S(i-1,j)
          else
            sx = sx - S(i,j)
          endif
          if (cnt == 0) then
            sx = 0
          else
            sx = sx / Dx
          endif
        endif

        cnt = 0

        ! calculate sy, similarly
        if (((j+j_off) == gjsc) .and. (.not. CS%reentrant_y)) then ! at south computational bdry
          if (ISS%hmask(i,j+1) == 1 .or. ISS%hmask(i,j+1) == 3) then
            sy = (S(i,j+1)-S(i,j))/dyh
          else
            sy = 0
          endif
        elseif (((j+j_off) == gjec) .and. (.not. CS%reentrant_y)) then ! at north computational bdry
          if (ISS%hmask(i,j-1) == 1 .or. ISS%hmask(i,j-1) == 3) then
            sy = (S(i,j)-S(i,j-1))/dyh
          else
            sy = 0
          endif
        else ! interior
          if (ISS%hmask(i,j+1) == 1 .or. ISS%hmask(i,j+1) == 3) then
            cnt = cnt+1
            Dy = dyh + G%dyT(i,j+1)
            sy = S(i,j+1)
          else
            sy = S(i,j)
          endif
          if (ISS%hmask(i,j-1) == 1 .or. ISS%hmask(i,j-1) == 3) then
            cnt = cnt+1
            Dy = dyh + G%dyT(i,j-1)
            sy = sy - S(i,j-1)
          else
            sy = sy - S(i,j)
          endif
          if (cnt == 0) then
            sy = 0
          else
            sy = sy / Dy
          endif
        endif

        if (CS%max_surface_slope>0) then
          scale = min(CS%max_surface_slope/sqrt((sx**2)+(sy**2)),1.0)
          sx = scale*sx; sy = scale*sy
        endif

        sx_e(i,j) = (-.25 * G%areaT(i,j)) * ((rho * grav) * (max(ISS%h_shelf(i,j),CS%min_h_shelf) * sx))
        sy_e(i,j) = (-.25 * G%areaT(i,j)) * ((rho * grav) * (max(ISS%h_shelf(i,j),CS%min_h_shelf) * sy))

        CS%sx_shelf(i,j) = sx ; CS%sy_shelf(i,j) = sy

        !Stress (Neumann) boundary conditions
        if (CS%ground_frac(i,j) == 1) then
          neumann_val = ((.5 * grav) * (rho * max(ISS%h_shelf(i,j),CS%min_h_shelf)**2 - rhow * CS%bed_elev(i,j)**2))
        else
          neumann_val = (.5 * grav) * ((1-rho/rhow) * (rho * max(ISS%h_shelf(i,j),CS%min_h_shelf)**2))
        endif
        if ((CS%u_face_mask_bdry(I-1,j) == 2) .OR. &
          ((ISS%hmask(i-1,j) == 0 .OR. ISS%hmask(i-1,j) == 2) .AND. (CS%reentrant_x .OR. (i+i_off /= gisc)))) then
          ! left face of the cell is at a stress boundary
          ! the depth-integrated longitudinal stress is equal to the difference of depth-integrated
          ! pressure on either side of the face
          ! on the ice side, it is rho g h^2 / 2
          ! on the ocean side, it is rhow g (delta OD)^2 / 2
          ! OD can be zero under the ice; but it is ASSUMED on the ice-free side of the face, topography elevation
          !     is not above the base of the ice in the current cell

          ! Note the negative sign due to the direction of the normal vector
          taudx(I-1,J-1) = taudx(I-1,J-1) - .5 * dyh * neumann_val
          taudx(I-1,J) = taudx(I-1,J) - .5 * dyh * neumann_val
        endif

        if ((CS%u_face_mask_bdry(I,j) == 2) .OR. &
          ((ISS%hmask(i+1,j) == 0 .OR. ISS%hmask(i+1,j) == 2) .and. (CS%reentrant_x .OR. (i+i_off /= giec)))) then
          ! east face of the cell is at a stress boundary
          taudx(I,J-1) = taudx(I,J-1) + .5 * dyh * neumann_val
          taudx(I,J) = taudx(I,J) + .5 * dyh * neumann_val
        endif

        if ((CS%v_face_mask_bdry(i,J-1) == 2) .OR. &
          ((ISS%hmask(i,j-1) == 0 .OR. ISS%hmask(i,j-1) == 2) .and. (CS%reentrant_y .OR. (j+j_off /= gjsc)))) then
          ! south face of the cell is at a stress boundary
          taudy(I-1,J-1) = taudy(I-1,J-1) - .5 * dxh * neumann_val
          taudy(I,J-1) = taudy(I,J-1) - .5 * dxh * neumann_val
        endif

        if ((CS%v_face_mask_bdry(i,J) == 2) .OR. &
          ((ISS%hmask(i,j+1) == 0 .OR. ISS%hmask(i,j+1) == 2) .and. (CS%reentrant_y .OR. (j+j_off /= gjec)))) then
          ! north face of the cell is at a stress boundary
          taudy(I-1,J) = taudy(I-1,J) + .5 * dxh * neumann_val
          taudy(I,J) = taudy(I,J) + .5 * dxh * neumann_val
        endif
      endif
    enddo
  enddo

  do J=jsc-2,jec+1; do I=isc-2,iec+1
    taudx(I,J) = taudx(I,J) + ((sx_e(i,j)+sx_e(i+1,j+1)) + (sx_e(i+1,j)+sx_e(i,j+1)))
    taudy(I,J) = taudy(I,J) + ((sy_e(i,j)+sy_e(i+1,j+1)) + (sy_e(i+1,j)+sy_e(i,j+1)))
  enddo; enddo
end subroutine calc_shelf_driving_stress

subroutine CG_action(CS, uret, vret, u_shlf, v_shlf, Phi, Phisub, umask, vmask, hmask, H_node, &
                     ice_visc, float_cond, bathyT, basal_trac, G, US, is, ie, js, je, dens_ratio)

  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                         intent(inout) :: uret !< The retarding stresses working at u-points [R L3 Z T-2 ~> kg m s-2].
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                         intent(inout) :: vret !< The retarding stresses working at v-points [R L3 Z T-2 ~> kg m s-2].
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                         intent(in)   :: Phi !< The gradients of bilinear basis elements at Gaussian
                                             !! quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real, dimension(:,:,:,:,:,:), &
                         intent(in)    :: Phisub !< Quadrature structure weights at subgridscale
                                            !! locations for finite element calculations [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: u_shlf  !< The zonal ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: v_shlf  !< The meridional ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: umask !< A coded mask indicating the nature of the
                                             !! zonal flow at the corner point
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: vmask !< A coded mask indicating the nature of the
                                             !! meridional flow at the corner point
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: H_node !< The ice shelf thickness at nodal (corner)
                                             !! points [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G),CS%visc_qps), &
                         intent(in)    :: ice_visc !< A field related to the ice viscosity from Glen's
                                               !! flow law [R L4 Z T-1 ~> kg m2 s-1].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: float_cond !< If GL_regularize=true, an array indicating where the ice
                                                !! shelf is floating: 0 if floating, 1 if not
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: bathyT !< The depth of ocean bathymetry at tracer points
                                                 !! relative to sea-level [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: basal_trac  !< Area-integrated taub_beta field related to the nonlinear
                                                !! part of the "linearized" basal stress [R L3 T-1 ~> kg s-1].

  real,                  intent(in)    :: dens_ratio !< The density of ice divided by the density
                                                     !! of seawater, nondimensional
  type(unit_scale_type), intent(in)    :: US  !< A structure containing unit conversion factors
  integer,               intent(in)    :: is  !< The starting i-index to work on
  integer,               intent(in)    :: ie  !< The ending i-index to work on
  integer,               intent(in)    :: js  !< The starting j-index to work on
  integer,               intent(in)    :: je  !< The ending j-index to work on

! the linear action of the matrix on (u,v) with bilinear finite elements
! as of now everything is passed in so no grid pointers or anything of the sort have to be dereferenced,
! but this may change pursuant to conversations with others
!
! is & ie are the cells over which the iteration is done; this may change between calls to this subroutine
!     in order to make less frequent halo updates

! the linear action of the matrix on (u,v) with bilinear finite elements
! Phi has the form
! Phi(k,q,i,j) - applies to cell i,j

    !  3 - 4
    !  |   |
    !  1 - 2

! Phi(2*k-1,q,i,j) gives d(Phi_k)/dx at quadrature point q
! Phi(2*k,q,i,j) gives d(Phi_k)/dy at quadrature point q
! Phi_k is equal to 1 at vertex k, and 0 at vertex l /= k, and bilinear

  real :: ux, uy, vx, vy ! Components of velocity shears or divergence [T-1 ~> s-1]
  real :: uq, vq  ! Interpolated velocities [L T-1 ~> m s-1]
  integer :: iq, jq, iphi, jphi, i, j, ilq, jlq, Itgt, Jtgt, qp, qpv
  logical :: visc_qp4
  real, dimension(2) :: xquad
  real, dimension(2,2) :: Ucell, Vcell, Hcell, Usub, Vsub
  real, dimension(2,2,4) :: uret_qp, vret_qp
  real, dimension(SZDIB_(G),SZDJB_(G),4) :: uret_b, vret_b

  xquad(1) = .5 * (1-sqrt(1./3)) ; xquad(2) = .5 * (1+sqrt(1./3))

  if (CS%visc_qps == 4) then
    visc_qp4=.true.
  else
    visc_qp4=.false.
    qpv = 1
  endif

  uret(:,:) = 0.0; vret(:,:)=0.0
  uret_b(:,:,:)=0.0 ; vret_b(:,:,:)=0.0

  do j=js,je ; do i=is,ie ; if (hmask(i,j) == 1 .or. hmask(i,j)==3) then

    uret_qp(:,:,:)=0.0; vret_qp(:,:,:)=0.0

      do iq=1,2 ; do jq=1,2

        qp = 2*(jq-1)+iq !current quad point

        uq = ((u_shlf(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
              (u_shlf(I,J) * (xquad(iq) * xquad(jq)))) + &
             ((u_shlf(I,J-1) * (xquad(iq) * xquad(3-jq))) + &
              (u_shlf(I-1,J) * (xquad(3-iq) * xquad(jq))))

        vq = ((v_shlf(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
              (v_shlf(I,J) * (xquad(iq) * xquad(jq)))) + &
             ((v_shlf(I,J-1) * (xquad(iq) * xquad(3-jq))) + &
              (v_shlf(I-1,J) * (xquad(3-iq) * xquad(jq))))

        ux = ((u_shlf(I-1,J-1) * Phi(1,qp,i,j)) + &
              (u_shlf(I,J) * Phi(7,qp,i,j))) + &
             ((u_shlf(I,J-1) * Phi(3,qp,i,j)) + &
              (u_shlf(I-1,J) * Phi(5,qp,i,j)))

        vx = ((v_shlf(I-1,J-1) * Phi(1,qp,i,j)) + &
              (v_shlf(I,J) * Phi(7,qp,i,j))) + &
             ((v_shlf(I,J-1) * Phi(3,qp,i,j)) + &
              (v_shlf(I-1,J) * Phi(5,qp,i,j)))

        uy = ((u_shlf(I-1,J-1) * Phi(2,qp,i,j)) + &
              (u_shlf(I,J) * Phi(8,qp,i,j))) + &
             ((u_shlf(I,J-1) * Phi(4,qp,i,j)) + &
              (u_shlf(I-1,J) * Phi(6,qp,i,j)))

        vy = ((v_shlf(I-1,J-1) * Phi(2,qp,i,j)) + &
              (v_shlf(I,J) * Phi(8,qp,i,j))) + &
             ((v_shlf(I,J-1) * Phi(4,qp,i,j)) + &
              (v_shlf(I-1,J) * Phi(6,qp,i,j)))

        if (visc_qp4) qpv = qp !current quad point for viscosity

        do jphi=1,2 ; Jtgt = J-2+jphi ; do iphi=1,2 ; Itgt = I-2+iphi
          if (umask(Itgt,Jtgt) == 1) uret_qp(iphi,jphi,qp) = ice_visc(i,j,qpv) * &
            (((4*ux+2*vy) * Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)) + &
            ((uy+vx) * Phi(2*(2*(jphi-1)+iphi),qp,i,j)))
          if (vmask(Itgt,Jtgt) == 1) vret_qp(iphi,jphi,qp) = ice_visc(i,j,qpv) * &
            (((uy+vx) * Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)) + &
            ((4*vy+2*ux) * Phi(2*(2*(jphi-1)+iphi),qp,i,j)))

          if (float_cond(i,j) == 0) then
            ilq = 1 ; if (iq == iphi) ilq = 2
            jlq = 1 ; if (jq == jphi) jlq = 2
            if (umask(Itgt,Jtgt) == 1) uret_qp(iphi,jphi,qp) = uret_qp(iphi,jphi,qp) +  &
              ((basal_trac(i,j) * uq) * (xquad(ilq) * xquad(jlq)))
            if (vmask(Itgt,Jtgt) == 1) vret_qp(iphi,jphi,qp) = vret_qp(iphi,jphi,qp) +  &
              ((basal_trac(i,j) * vq) * (xquad(ilq) * xquad(jlq)))
          endif
        enddo ; enddo
      enddo ; enddo

      !element contribution to SW node (node 1, which sees the current element as element 4)
      uret_b(I-1,J-1,4) = 0.25*((uret_qp(1,1,1)+uret_qp(1,1,4))+(uret_qp(1,1,2)+uret_qp(1,1,3)))
      vret_b(I-1,J-1,4) = 0.25*((vret_qp(1,1,1)+vret_qp(1,1,4))+(vret_qp(1,1,2)+vret_qp(1,1,3)))

      !element contribution to NW node (node 3, which sees the current element as element 2)
      uret_b(I-1,J  ,2) = 0.25*((uret_qp(1,2,1)+uret_qp(1,2,4))+(uret_qp(1,2,2)+uret_qp(1,2,3)))
      vret_b(I-1,J  ,2) = 0.25*((vret_qp(1,2,1)+vret_qp(1,2,4))+(vret_qp(1,2,2)+vret_qp(1,2,3)))

      !element contribution to SE node (node 2, which sees the current element as element 3)
      uret_b(I  ,J-1,3) = 0.25*((uret_qp(2,1,1)+uret_qp(2,1,4))+(uret_qp(2,1,2)+uret_qp(2,1,3)))
      vret_b(I  ,J-1,3) = 0.25*((vret_qp(2,1,1)+vret_qp(2,1,4))+(vret_qp(2,1,2)+vret_qp(2,1,3)))

      !element contribution to NE node (node 4, which sees the current element as element 1)
      uret_b(I  ,J  ,1) = 0.25*((uret_qp(2,2,1)+uret_qp(2,2,4))+(uret_qp(2,2,2)+uret_qp(2,2,3)))
      vret_b(I  ,J  ,1) = 0.25*((vret_qp(2,2,1)+vret_qp(2,2,4))+(vret_qp(2,2,2)+vret_qp(2,2,3)))

      if (float_cond(i,j) == 1) then
        Ucell(:,:) = u_shlf(I-1:I,J-1:J) ; Vcell(:,:) = v_shlf(I-1:I,J-1:J)
        Hcell(:,:) = H_node(I-1:I,J-1:J)

        call CG_action_subgrid_basal(Phisub, Hcell, Ucell, Vcell, &
                                     bathyT(i,j), dens_ratio, Usub, Vsub)

        if (umask(I-1,J-1) == 1) uret_b(I-1,J-1,4) = uret_b(I-1,J-1,4) + (Usub(1,1) * basal_trac(i,j))
        if (umask(I-1,J  ) == 1) uret_b(I-1,J  ,2) = uret_b(I-1,J  ,2) + (Usub(1,2) * basal_trac(i,j))
        if (umask(I  ,J-1) == 1) uret_b(I  ,J-1,3) = uret_b(I  ,J-1,3) + (Usub(2,1) * basal_trac(i,j))
        if (umask(I  ,J  ) == 1) uret_b(I  ,J  ,1) = uret_b(I  ,J  ,1) + (Usub(2,2) * basal_trac(i,j))

        if (vmask(I-1,J-1) == 1) vret_b(I-1,J-1,4) = vret_b(I-1,J-1,4) + (Vsub(1,1) * basal_trac(i,j))
        if (vmask(I-1,J  ) == 1) vret_b(I-1,J  ,2) = vret_b(I-1,J  ,2) + (Vsub(1,2) * basal_trac(i,j))
        if (vmask(I  ,J-1) == 1) vret_b(I  ,J-1,3) = vret_b(I  ,J-1,3) + (Vsub(2,1) * basal_trac(i,j))
        if (vmask(I  ,J  ) == 1) vret_b(I  ,J  ,1) = vret_b(I  ,J  ,1) + (Vsub(2,2) * basal_trac(i,j))
      endif
  endif ; enddo ; enddo

  do J=js-1,je ; do I=is-1,ie
    uret(I,J) = (uret_b(I,J,1)+uret_b(I,J,4)) + (uret_b(I,J,2)+uret_b(I,J,3))
    vret(I,J) = (vret_b(I,J,1)+vret_b(I,J,4)) + (vret_b(I,J,2)+vret_b(I,J,3))
  enddo; enddo

end subroutine CG_action

subroutine CG_action_subgrid_basal(Phisub, H, U, V, bathyT, dens_ratio, Ucontr, Vcontr)
  real, dimension(:,:,:,:,:,:), &
                        intent(in)    :: Phisub !< Quadrature structure weights at subgridscale
                                            !! locations for finite element calculations [nondim]
  real, dimension(2,2), intent(in)    :: H  !< The ice shelf thickness at nodal (corner) points [Z ~> m].
  real, dimension(2,2), intent(in)    :: U  !< The zonal ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(2,2), intent(in)    :: V  !< The meridional ice shelf velocity at vertices [L T-1 ~> m s-1]
  real,                 intent(in)    :: bathyT !< The depth of ocean bathymetry at tracer points
                                            !! relative to sea-level [Z ~> m].
  real,                 intent(in)    :: dens_ratio !< The density of ice divided by the density
                                            !! of seawater [nondim]
  real, dimension(2,2), intent(out)   :: Ucontr !< The areal average of u-velocities where the ice shelf
                                            !! is grounded, or 0 where it is floating [L T-1 ~> m s-1].
  real, dimension(2,2), intent(out)   :: Vcontr !< The areal average of v-velocities where the ice shelf
                                            !! is grounded, or 0 where it is floating [L T-1 ~> m s-1].

  real, dimension(SIZE(Phisub,3),SIZE(Phisub,3),2,2) :: Ucontr_sub, Vcontr_sub ! The contributions to Ucontr and Vcontr
                                                                               !! at each sub-cell
  real, dimension(2,2,SIZE(Phisub,3),SIZE(Phisub,3)) :: uloc_arr !The local sub-cell u-velocity [L T-1 ~> m s-1]
  real, dimension(2,2,SIZE(Phisub,3),SIZE(Phisub,3)) :: vloc_arr !The local sub-cell v-velocity [L T-1 ~> m s-1]
  real, dimension(2,2) :: Ucontr_q, Vcontr_q !Contributions to a node from each quadrature point in a sub-grid cell
  real    :: subarea ! The fractional sub-cell area [nondim]
  real    :: hloc    ! The local sub-cell ice thickness [Z ~> m]
  integer :: nsub, i, j, qx, qy, m, n

  nsub = size(Phisub,3)
  subarea = 1.0 / (nsub**2)

  uloc_arr(:,:,:,:) = 0.0; vloc_arr(:,:,:,:)=0.0

  do j=1,nsub ; do i=1,nsub;  do qy=1,2 ; do qx=1,2
    hloc = ((Phisub(qx,qy,i,j,1,1)*H(1,1)) + (Phisub(qx,qy,i,j,2,2)*H(2,2))) + &
           ((Phisub(qx,qy,i,j,1,2)*H(1,2)) + (Phisub(qx,qy,i,j,2,1)*H(2,1)))
    if (dens_ratio * hloc - bathyT > 0) then
      uloc_arr(qx,qy,i,j) = (((Phisub(qx,qy,i,j,1,1) * U(1,1)) + (Phisub(qx,qy,i,j,2,2) * U(2,2))) + &
                             ((Phisub(qx,qy,i,j,1,2) * U(1,2)) + (Phisub(qx,qy,i,j,2,1) * U(2,1))))
      vloc_arr(qx,qy,i,j) = (((Phisub(qx,qy,i,j,1,1) * V(1,1)) + (Phisub(qx,qy,i,j,2,2) * V(2,2))) + &
                             ((Phisub(qx,qy,i,j,1,2) * V(1,2)) + (Phisub(qx,qy,i,j,2,1) * V(2,1))))
    endif
  enddo; enddo ; enddo ; enddo

  do n=1,2 ; do m=1,2 ; do j=1,nsub ; do i=1,nsub
    do qy=1,2 ; do qx=1,2
      !calculate quadrature point contributions for the sub-cell, to each node
        Ucontr_q(qx,qy) = Phisub(qx,qy,i,j,m,n) * uloc_arr(qx,qy,i,j)
        Vcontr_q(qx,qy) = Phisub(qx,qy,i,j,m,n) * vloc_arr(qx,qy,i,j)
    enddo; enddo

    !calculate sub-cell contribution to each node by summing up quadrature point contributions from the sub-cell
    Ucontr_sub(i,j,m,n) = (subarea * 0.25) * ((Ucontr_q(1,1) + Ucontr_q(2,2)) + (Ucontr_q(1,2)+Ucontr_q(2,1)))
    Vcontr_sub(i,j,m,n) = (subarea * 0.25) * ((Vcontr_q(1,1) + Vcontr_q(2,2)) + (Vcontr_q(1,2)+Vcontr_q(2,1)))
  enddo; enddo ; enddo ; enddo

  !sum up the sub-cell contributions to each node
  do n=1,2 ; do m=1,2
    call sum_square_matrix(Ucontr(m,n),Ucontr_sub(:,:,m,n),nsub)
    call sum_square_matrix(Vcontr(m,n),Vcontr_sub(:,:,m,n),nsub)
  enddo ; enddo

end subroutine CG_action_subgrid_basal


!! Returns the sum of the elements in a square matrix. This sum is bitwise identical even if the matrices are rotated.
subroutine sum_square_matrix(sum_out, mat_in, n)
  integer, intent(in) :: n !< The length and width of each matrix in mat_in
  real, dimension(n,n), intent(in) :: mat_in !< The n x n matrix whose elements will be summed
  real, intent(out) :: sum_out !< The sum of the elements of matrix mat_in
  integer :: s0,e0,s1,e1

  sum_out=0.0

  s0=1; e0=n

  !start by summing elements on outer edges of matrix
  do while (s0<e0)

    !corners
    sum_out = sum_out + ( (mat_in(s0,s0) + mat_in(e0,e0)) + (mat_in(e0,s0) + mat_in(s0,e0)) )

    s1=s0+1; e1=e0-1

    do while (s1<e1) !non-corners

      sum_out = sum_out + &
                ( ( (mat_in(s0,s1) + mat_in(s1,s0)) + (mat_in(e0,e1) + mat_in(e1,e0)) ) + &
                  ( (mat_in(e1,s0) + mat_in(e0,s1)) + (mat_in(s1,e0) + mat_in(s0,e1)) ) )

      s1=s1+1 ; e1=e1-1
    enddo

    !center element of an edge
    if (s1==e1) sum_out = sum_out + ( (mat_in(s1,s0) + mat_in(e1,e0)) + (mat_in(e0,e1) + mat_in(s0,s1)) )

    s0=s0+1 ; e0=e0-1 !next loop iteration using new edges that are one element inward of the current edges
  enddo

  !center element of entire matrix
  if (s0==e0) sum_out = sum_out + mat_in(s0,e0)

end subroutine sum_square_matrix

!> returns the diagonal entries of the matrix for a Jacobi preconditioning
subroutine matrix_diagonal(CS, G, US, float_cond, H_node, ice_visc, basal_trac, hmask, dens_ratio, &
                           Phi, Phisub, u_diagonal, v_diagonal)

  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: float_cond !< If GL_regularize=true, indicates cells containing
                                                !! the grounding line (float_cond=1) or not (float_cond=0)
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal
                                                 !! (corner) points [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G),CS%visc_qps), &
                          intent(in)    :: ice_visc !< A field related to the ice viscosity from Glen's
                                                !! flow law [R L4 Z T-1 ~> kg m2 s-1].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: basal_trac !< Area-integrated taub_beta field related to the nonlinear
                                                !! part of the "linearized" basal stress [R L3 T-1 ~> kg s-1].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real,                   intent(in)    :: dens_ratio !< The density of ice divided by the density
                                                     !! of seawater [nondim]
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: Phi !< The gradients of bilinear basis elements at Gaussian
                                             !! quadrature points surrounding the cell vertices [L-1 ~> m-1]
  real, dimension(:,:,:,:,:,:), intent(in) :: Phisub !< Quadrature structure weights at subgridscale
                                            !! locations for finite element calculations [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_diagonal !< The diagonal elements of the u-velocity
                                            !! matrix from the left-hand side of the solver [R L2 Z T-1 ~> kg s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_diagonal  !< The diagonal elements of the v-velocity
                                            !! matrix from the left-hand side of the solver [R L2 Z T-1 ~> kg s-1]


! returns the diagonal entries of the matrix for a Jacobi preconditioning

  real :: ux, uy, vx, vy ! Interpolated weight gradients [L-1 ~> m-1]
  real :: uq, vq
  real, dimension(2)   :: xquad
  real, dimension(2,2) :: Hcell, sub_ground
  real, dimension(2,2,4) :: u_diag_qp, v_diag_qp
  real, dimension(SZDIB_(G),SZDJB_(G),4) :: u_diag_b, v_diag_b
  logical :: visc_qp4
  integer :: i, j, isc, jsc, iec, jec, iphi, jphi, iq, jq, ilq, jlq, Itgt, Jtgt, qp, qpv

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  xquad(1) = .5 * (1-sqrt(1./3)) ; xquad(2) = .5 * (1+sqrt(1./3))

  if (CS%visc_qps == 4) then
    visc_qp4=.true.
  else
    visc_qp4=.false.
    qpv = 1
  endif

  u_diag_b(:,:,:)=0.0
  v_diag_b(:,:,:)=0.0

  do j=jsc-1,jec+1 ; do i=isc-1,iec+1 ; if (hmask(i,j) == 1 .or. hmask(i,j)==3) then

    ! Phi(2*i-1,j) gives d(Phi_i)/dx at quadrature point j
    ! Phi(2*i,j) gives d(Phi_i)/dy at quadrature point j

    u_diag_qp(:,:,:)=0.0; v_diag_qp(:,:,:)=0.0

    do iq=1,2 ; do jq=1,2

      qp = 2*(jq-1)+iq !current quad point
      if (visc_qp4) qpv = qp !current quad point for viscosity

      do jphi=1,2 ; Jtgt = J-2+jphi ; do iphi=1,2 ; Itgt = I-2+iphi

        ilq = 1 ; if (iq == iphi) ilq = 2
        jlq = 1 ; if (jq == jphi) jlq = 2

        if (CS%umask(Itgt,Jtgt) == 1) then

          ux = Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)
          uy = Phi(2*(2*(jphi-1)+iphi),qp,i,j)
          vx = 0.
          vy = 0.

          u_diag_qp(iphi,jphi,qp) = &
            ice_visc(i,j,qpv) * (((4*ux+2*vy) * Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)) + &
            ((uy+vx) * Phi(2*(2*(jphi-1)+iphi),qp,i,j)))

          if (float_cond(i,j) == 0) then
            uq = xquad(ilq) * xquad(jlq)
            u_diag_qp(iphi,jphi,qp) = u_diag_qp(iphi,jphi,qp) + &
              (basal_trac(i,j) * uq) * (xquad(ilq) * xquad(jlq))
          endif
        endif

        if (CS%vmask(Itgt,Jtgt) == 1) then

          vx = Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)
          vy = Phi(2*(2*(jphi-1)+iphi),qp,i,j)
          ux = 0.
          uy = 0.

          v_diag_qp(iphi,jphi,qp) = &
            ice_visc(i,j,qpv) * (((uy+vx) * Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)) + &
            ((4*vy+2*ux) * Phi(2*(2*(jphi-1)+iphi),qp,i,j)))

          if (float_cond(i,j) == 0) then
            vq = xquad(ilq) * xquad(jlq)
            v_diag_qp(iphi,jphi,qp) = v_diag_qp(iphi,jphi,qp) + &
              (basal_trac(i,j) * vq) * (xquad(ilq) * xquad(jlq))
          endif
        endif
      enddo ; enddo
    enddo ; enddo

    !element contribution to SW node (node 1, which sees the current element as element 4)
    u_diag_b(I-1,J-1,4) = 0.25*((u_diag_qp(1,1,1)+u_diag_qp(1,1,4))+(u_diag_qp(1,1,2)+u_diag_qp(1,1,3)))
    v_diag_b(I-1,J-1,4) = 0.25*((v_diag_qp(1,1,1)+v_diag_qp(1,1,4))+(v_diag_qp(1,1,2)+v_diag_qp(1,1,3)))

    !element contribution to NW node (node 3, which sees the current element as element 2)
    u_diag_b(I-1,J  ,2) = 0.25*((u_diag_qp(1,2,1)+u_diag_qp(1,2,4))+(u_diag_qp(1,2,2)+u_diag_qp(1,2,3)))
    v_diag_b(I-1,J  ,2) = 0.25*((v_diag_qp(1,2,1)+v_diag_qp(1,2,4))+(v_diag_qp(1,2,2)+v_diag_qp(1,2,3)))

    !element contribution to SE node (node 2, which sees the current element as element 3)
    u_diag_b(I  ,J-1,3) = 0.25*((u_diag_qp(2,1,1)+u_diag_qp(2,1,4))+(u_diag_qp(2,1,2)+u_diag_qp(2,1,3)))
    v_diag_b(I  ,J-1,3) = 0.25*((v_diag_qp(2,1,1)+v_diag_qp(2,1,4))+(v_diag_qp(2,1,2)+v_diag_qp(2,1,3)))

    !element contribution to NE node (node 4, which sees the current element as element 1)
    u_diag_b(I  ,J  ,1) = 0.25*((u_diag_qp(2,2,1)+u_diag_qp(2,2,4))+(u_diag_qp(2,2,2)+u_diag_qp(2,2,3)))
    v_diag_b(I  ,J  ,1) = 0.25*((v_diag_qp(2,2,1)+v_diag_qp(2,2,4))+(v_diag_qp(2,2,2)+v_diag_qp(2,2,3)))

    if (float_cond(i,j) == 1) then
      Hcell(:,:) = H_node(i-1:i,j-1:j)
      call CG_diagonal_subgrid_basal(Phisub, Hcell, CS%bed_elev(i,j), dens_ratio, sub_ground)

        if (CS%umask(I-1,J-1) == 1) u_diag_b(I-1,J-1,4) = u_diag_b(I-1,J-1,4) + (sub_ground(1,1) * basal_trac(i,j))
        if (CS%umask(I-1,J  ) == 1) u_diag_b(I-1,J  ,2) = u_diag_b(I-1,J  ,2) + (sub_ground(1,2) * basal_trac(i,j))
        if (CS%umask(I  ,J-1) == 1) u_diag_b(I  ,J-1,3) = u_diag_b(I  ,J-1,3) + (sub_ground(2,1) * basal_trac(i,j))
        if (CS%umask(I  ,J  ) == 1) u_diag_b(I  ,J  ,1) = u_diag_b(I  ,J  ,1) + (sub_ground(2,2) * basal_trac(i,j))

        if (CS%vmask(I-1,J-1) == 1) v_diag_b(I-1,J-1,4) = v_diag_b(I-1,J-1,4) + (sub_ground(1,1) * basal_trac(i,j))
        if (CS%vmask(I-1,J  ) == 1) v_diag_b(I-1,J  ,2) = v_diag_b(I-1,J  ,2) + (sub_ground(1,2) * basal_trac(i,j))
        if (CS%vmask(I  ,J-1) == 1) v_diag_b(I  ,J-1,3) = v_diag_b(I  ,J-1,3) + (sub_ground(2,1) * basal_trac(i,j))
        if (CS%vmask(I  ,J  ) == 1) v_diag_b(I  ,J  ,1) = v_diag_b(I  ,J  ,1) + (sub_ground(2,2) * basal_trac(i,j))
    endif
  endif ; enddo ; enddo

  do J=jsc-2,jec+1 ; do I=isc-2,iec+1
    u_diagonal(I,J) = (u_diag_b(I,J,1)+u_diag_b(I,J,4)) + (u_diag_b(I,J,2)+u_diag_b(I,J,3))
    v_diagonal(I,J) = (v_diag_b(I,J,1)+v_diag_b(I,J,4)) + (v_diag_b(I,J,2)+v_diag_b(I,J,3))
  enddo ; enddo

end subroutine matrix_diagonal

subroutine CG_diagonal_subgrid_basal (Phisub, H_node, bathyT, dens_ratio, f_grnd)
  real, dimension(:,:,:,:,:,:), &
                        intent(in) :: Phisub !< Quadrature structure weights at subgridscale
                                             !! locations for finite element calculations [nondim]
  real, dimension(2,2), intent(in) :: H_node !< The ice shelf thickness at nodal (corner)
                                             !! points [Z ~> m].
  real,              intent(in)    :: bathyT !< The depth of ocean bathymetry at tracer points [Z ~> m].
  real,              intent(in)    :: dens_ratio !< The density of ice divided by the density
                                                 !! of seawater [nondim]
  real, dimension(2,2), intent(out) :: f_grnd !< The weighted fraction of the sub-cell where the ice shelf
                                              !! is grounded [nondim]

  real, dimension(SIZE(Phisub,3),SIZE(Phisub,3),2,2) :: f_grnd_sub ! The contributions to nodal f_grnd
                                                                   !! from each sub-cell
  integer, dimension(2,2,SIZE(Phisub,3),SIZE(Phisub,3)) :: grnd_stat !0 at floating quad points, 1 at grounded
  real, dimension(2,2) :: f_grnd_q  !Contributions to a node from each quadrature point in a sub-grid cell
  real    :: subarea ! The fractional sub-cell area [nondim]
  real    :: hloc    ! The local sub-region thickness [Z ~> m]
  integer :: nsub, i, j, qx, qy, m, n

  nsub = size(Phisub,3)
  subarea = 1.0 / (nsub**2)

  grnd_stat(:,:,:,:)=0

  do j=1,nsub ; do i=1,nsub;  do qy=1,2 ; do qx=1,2
    hloc = ((Phisub(qx,qy,i,j,1,1)*H_node(1,1)) + (Phisub(qx,qy,i,j,2,2)*H_node(2,2))) + &
           ((Phisub(qx,qy,i,j,1,2)*H_node(1,2)) + (Phisub(qx,qy,i,j,2,1)*H_node(2,1)))
    if (dens_ratio * hloc - bathyT > 0) grnd_stat(qx,qy,i,j) = 1
  enddo; enddo ; enddo ; enddo

  do n=1,2 ; do m=1,2 ; do j=1,nsub ; do i=1,nsub
    do qy=1,2 ; do qx = 1,2
        f_grnd_q(qx,qy) = grnd_stat(qx,qy,i,j) * Phisub(qx,qy,i,j,m,n)**2
    enddo ; enddo
    !calculate sub-cell contribution to each node by summing up quadrature point contributions from the sub-cell
    f_grnd_sub(i,j,m,n) = (subarea * 0.25) * ((f_grnd_q(1,1) + f_grnd_q(2,2)) + (f_grnd_q(1,2)+f_grnd_q(2,1)))
  enddo ; enddo ; enddo ; enddo

  !sum up the sub-cell contributions to each node
  do n=1,2 ; do m=1,2
   call sum_square_matrix(f_grnd(m,n),f_grnd_sub(:,:,m,n),nsub)
  enddo ; enddo

end subroutine CG_diagonal_subgrid_basal

!> Post_data calls related to ice-sheet flux divergence, strain-rate, and deviatoric stress
subroutine IS_dynamics_post_data_2(CS, ISS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDIB_(G),SZDJB_(G)) :: H_node ! Ice shelf thickness at corners [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)) :: Hu  ! Ice shelf u_flux at corners [Z L T-1 ~> m2 s-1].
  real, dimension(SZDIB_(G),SZDJB_(G)) :: Hv  ! Ice shelf v_flux at corners [Z L T-1 ~> m2 s-1].
  real, dimension(SZDI_(G),SZDJ_(G)) :: Hux  ! Ice shelf d(u_flux)/dx at cell centers [Z T-1 ~> m s-1].
  real, dimension(SZDI_(G),SZDJ_(G)) :: Hvy  ! Ice shelf d(v_flux)/dy at cell centers [Z T-1 ~> m s-1].
  real, dimension(SZDI_(G),SZDJ_(G)) :: flux_div ! horizontal flux divergence div(uH) [Z T-1 ~> m s-1].
  real, dimension(SZDI_(G),SZDJ_(G),3) :: strain_rate ! strain-rate components xx,yy, and xy [T-1 ~> s-1]
  real, dimension(SZDI_(G),SZDJ_(G),2) :: p_strain_rate ! horizontal principal strain-rates [T-1 ~> s-1]
  real, dimension(SZDI_(G),SZDJ_(G),3) :: dev_stress ! deviatoric stress components xx,yy, and xy [R L Z T-2 ~> Pa]
  real, dimension(SZDI_(G),SZDJ_(G),2) :: p_dev_stress ! horizontal principal deviatoric stress [R L Z T-2 ~> Pa]
  real, dimension(SZDI_(G),SZDJ_(G))  :: ice_visc ! area-averaged ice viscosity [R L2 T-1 ~> Pa s]
  real :: p1,p2 ! Used to calculate strain-rate principal components [T-1 ~> s-1]
  integer :: i, j

  !Allocate the gradient basis functions for 1 cell-centered quadrature point per cell
  if (.not. associated(CS%PhiC)) then
    allocate(CS%PhiC(1:8,G%isc:G%iec,G%jsc:G%jec), source=0.0)
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      call bilinear_shape_fn_grid_1qp(G, i, j, CS%PhiC(:,i,j))
    enddo; enddo
  endif

  !Calculate flux divergence and its components
  if (CS%id_duHdx > 0 .or. CS%id_dvHdy > 0 .or. CS%id_fluxdiv > 0) then
    call interpolate_H_to_B(G, ISS%h_shelf, ISS%hmask, H_node, CS%min_h_shelf)

    Hu(:,:) = 0.0; Hv(:,:) = 0.0; Hux(:,:) = 0.0 ; Hvy(:,:) = 0.0 ; flux_div(:,:) = 0.0
    do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
      if (CS%umask(I,J) > 0) then
        Hu(I,J) = (H_node(I,J) * CS%u_shelf(I,J))
      endif
      if (CS%vmask(I,J) > 0) then
        Hv(I,J) = (H_node(I,J) * CS%v_shelf(I,J))
      endif
    enddo; enddo

    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 3)) then
        !components of flux divergence at cell centers
        Hux(i,j) = (((Hu(I-1,J-1) * CS%PhiC(1,i,j)) + (Hu(I,J  ) * CS%PhiC(7,i,j))) + &
                    ((Hu(I-1,J  ) * CS%PhiC(5,i,j)) + (Hu(I,J-1) * CS%PhiC(3,i,j))))

        Hvy(i,j) = (((Hv(I-1,J-1) * CS%PhiC(2,i,j)) + (Hv(I,J  ) * CS%PhiC(8,i,j))) + &
                    ((Hv(I-1,J  ) * CS%PhiC(6,i,j)) + (Hv(I,J-1) * CS%PhiC(4,i,j))))
        flux_div(i,j) = Hux(i,j) + Hvy(i,j)
      endif
    enddo ; enddo

    if (CS%id_duHdx > 0)   call post_data(CS%id_duHdx, Hux, CS%diag)
    if (CS%id_dvHdy > 0)   call post_data(CS%id_dvHdy, Hvy, CS%diag)
    if (CS%id_fluxdiv > 0) call post_data(CS%id_fluxdiv, flux_div, CS%diag)
  endif

  if (CS%id_devstress_xx > 0  .or. CS%id_devstress_yy > 0  .or. CS%id_devstress_xy > 0  .or. &
      CS%id_strainrate_xx > 0 .or. CS%id_strainrate_yy > 0 .or. CS%id_strainrate_xy > 0 .or. &
      CS%id_pdevstress_1 > 0  .or. CS%id_pdevstress_2 > 0  .or. &
      CS%id_pstrainrate_1 > 0 .or. CS%id_pstrainrate_2 > 0) then

    strain_rate(:,:,:) = 0.0
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      !strain-rates at cell centers
      if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 3)) then
        !strain_rate(:,:,1) = strain_rate_xx(:,:) = ux(:,:)
        strain_rate(i,j,1) = (((CS%u_shelf(I-1,J-1) * CS%PhiC(1,i,j)) + (CS%u_shelf(I,J  ) * CS%PhiC(7,i,j))) + &
                              ((CS%u_shelf(I-1,J  ) * CS%PhiC(5,i,j)) + (CS%u_shelf(I,J-1) * CS%PhiC(3,i,j))))
        !strain_rate(:,:,2) = strain_rate_yy(:,:) = uy(:,:)
        strain_rate(i,j,2) = (((CS%v_shelf(I-1,J-1) * CS%PhiC(2,i,j)) + (CS%v_shelf(I,J  ) * CS%PhiC(8,i,j))) + &
                              ((CS%v_shelf(I-1,J  ) * CS%PhiC(6,i,j)) + (CS%v_shelf(I,J-1) * CS%PhiC(4,i,j))))
        !strain_rate(:,:,3) = strain_rate_xy(:,:) = 0.5 * (uy(:,:) + vy(:,:))
        strain_rate(i,j,3) = 0.5 * ((((CS%u_shelf(I-1,J-1) * CS%PhiC(2,i,j)) + (CS%u_shelf(I,J  ) * CS%PhiC(8,i,j))) + &
                                     ((CS%u_shelf(I-1,J  ) * CS%PhiC(6,i,j)) + (CS%u_shelf(I,J-1) * CS%PhiC(4,i,j))))+ &
                                    (((CS%v_shelf(I-1,J-1) * CS%PhiC(1,i,j)) + (CS%v_shelf(I,J  ) * CS%PhiC(7,i,j))) + &
                                     ((CS%v_shelf(I-1,J  ) * CS%PhiC(5,i,j)) + (CS%v_shelf(I,J-1) * CS%PhiC(3,i,j)))))
      endif
    enddo ; enddo


    if (CS%id_strainrate_xx > 0) call post_data(CS%id_strainrate_xx, strain_rate(:,:,1), CS%diag)
    if (CS%id_strainrate_yy > 0) call post_data(CS%id_strainrate_yy, strain_rate(:,:,2), CS%diag)
    if (CS%id_strainrate_xy > 0) call post_data(CS%id_strainrate_xy, strain_rate(:,:,3), CS%diag)

    if (CS%id_pstrainrate_1 > 0 .or. CS%id_pstrainrate_2 > 0 .or. &
        CS%id_pdevstress_1  > 0 .or. CS%id_pdevstress_2  > 0) then
      p_strain_rate(:,:,:) = 0.0
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        p1 = 0.5*( strain_rate(i,j,1) + strain_rate(i,j,2))
        p2 = sqrt( (( 0.5 * (strain_rate(i,j,1) - strain_rate(i,j,2)) )**2) + (strain_rate(i,j,3)**2) )
        p_strain_rate(i,j,1) = p1+p2 !Max horizontal principal strain-rate
        p_strain_rate(i,j,2) = p1-p2 !Min horizontal principal strain-rate
      enddo ; enddo

      if (CS%id_pstrainrate_1 > 0) call post_data(CS%id_pstrainrate_1, p_strain_rate(:,:,1), CS%diag)
      if (CS%id_pstrainrate_2 > 0) call post_data(CS%id_pstrainrate_2, p_strain_rate(:,:,2), CS%diag)
    endif

    if (CS%id_devstress_xx > 0 .or. CS%id_devstress_yy > 0 .or. CS%id_devstress_xy > 0 .or. &
        CS%id_pdevstress_1 > 0 .or. CS%id_pdevstress_2 > 0) then

      call ice_visc_diag(CS,G,ice_visc)

      if (CS%id_devstress_xx > 0 .or. CS%id_devstress_yy > 0 .or. CS%id_devstress_xy > 0) then
        dev_stress(:,:,:)=0.0
        do j=G%jsc,G%jec ; do i=G%isc,G%iec
          if (ISS%h_shelf(i,j)>0) then
            dev_stress(i,j,1) = 2*ice_visc(i,j)*strain_rate(i,j,1)/ISS%h_shelf(i,j) !deviatoric stress xx
            dev_stress(i,j,2) = 2*ice_visc(i,j)*strain_rate(i,j,2)/ISS%h_shelf(i,j) !deviatoric stress yy
            dev_stress(i,j,3) = 2*ice_visc(i,j)*strain_rate(i,j,3)/ISS%h_shelf(i,j) !deviatoric stress xy
          endif
        enddo; enddo
        if (CS%id_devstress_xx > 0) call post_data(CS%id_devstress_xx, dev_stress(:,:,1), CS%diag)
        if (CS%id_devstress_yy > 0) call post_data(CS%id_devstress_yy, dev_stress(:,:,2), CS%diag)
        if (CS%id_devstress_xy > 0) call post_data(CS%id_devstress_xy, dev_stress(:,:,3), CS%diag)
      endif

      if (CS%id_pdevstress_1 > 0 .or. CS%id_pdevstress_2 > 0) then
        p_dev_stress(:,:,:)=0.0
        do j=G%jsc,G%jec ; do i=G%isc,G%iec
          if (ISS%h_shelf(i,j)>0) then
            p_dev_stress(i,j,1) = 2*ice_visc(i,j)*p_strain_rate(i,j,1)/ISS%h_shelf(i,j) !max horiz principal dev stress
            p_dev_stress(i,j,2) = 2*ice_visc(i,j)*p_strain_rate(i,j,2)/ISS%h_shelf(i,j) !min horiz principal dev stress
          endif
        enddo; enddo
        if (CS%id_pdevstress_1 > 0) call post_data(CS%id_pdevstress_1, p_dev_stress(:,:,1), CS%diag)
        if (CS%id_pdevstress_2 > 0) call post_data(CS%id_pdevstress_2, p_dev_stress(:,:,2), CS%diag)
      endif
    endif
  endif
end subroutine IS_dynamics_post_data_2

!> Update depth integrated viscosity, based on horizontal strain rates
subroutine calc_shelf_visc(CS, ISS, G, US, u_shlf, v_shlf)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                          intent(inout) :: u_shlf !< The zonal ice shelf velocity [L T-1 ~> m s-1].
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                          intent(inout) :: v_shlf !< The meridional ice shelf velocity [L T-1 ~> m s-1].

! update DEPTH_INTEGRATED viscosity, based on horizontal strain rates - this is for bilinear FEM solve


! this may be subject to change later... to make it "hybrid"
!  real, dimension(SZDIB_(G),SZDJB_(G)) ::  eII, ux, uy, vx, vy
  integer :: i, j, iscq, iecq, jscq, jecq, isd, jsd, ied, jed, iegq, jegq, iq, jq
  integer :: giec, gjec, gisc, gjsc, isc, jsc, iec, jec, is, js
  real :: Visc_coef, n_g
  real :: ux, uy, vx, vy
  real :: eps_min   ! Velocity shears [T-1 ~> s-1]
  logical :: model_qp1, model_qp4

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
  iegq = G%iegB ; jegq = G%jegB
  gisc = G%domain%nihalo+1 ; gjsc = G%domain%njhalo+1
  giec = G%domain%niglobal+gisc ; gjec = G%domain%njglobal+gjsc
  is = iscq - 1; js = jscq - 1

  if (trim(CS%ice_viscosity_compute) == "MODEL") then
    if (CS%visc_qps==1) then
      model_qp1=.true.
      model_qp4=.false.
    else
      model_qp1=.false.
      model_qp4=.true.
    endif
  endif

  n_g = CS%n_glen; eps_min = CS%eps_glen_min

  do j=jsc,jec ; do i=isc,iec

    if ((ISS%hmask(i,j) == 1) .OR. (ISS%hmask(i,j) == 3)) then

      if (trim(CS%ice_viscosity_compute) == "CONSTANT") then
        CS%ice_visc(i,j,1) = 1e15 * (US%kg_m3_to_R*US%m_to_L*US%m_s_to_L_T) * &
                             (G%areaT(i,j) * max(ISS%h_shelf(i,j),CS%min_h_shelf))
        ! constant viscocity for debugging
      elseif (trim(CS%ice_viscosity_compute) == "OBS") then
        if (CS%AGlen_visc(i,j) >0) then
          CS%ice_visc(i,j,1) = (G%areaT(i,j) * max(ISS%h_shelf(i,j),CS%min_h_shelf)) * &
                               max(CS%AGlen_visc(i,j) ,CS%min_ice_visc)
        endif
        ! Here CS%Aglen_visc(i,j) is the ice viscosity [Pa s ~> R L2 T-1] computed from obs and read from a file
      elseif (model_qp1) then
        !calculate viscosity at 1 cell-centered quadrature point per cell

        Visc_coef = (CS%AGlen_visc(i,j))**(-1./n_g)
        ! Units of Aglen_visc [Pa-(n_g) s-1]

        ux = ((u_shlf(I-1,J-1) * CS%PhiC(1,i,j)) + &
              (u_shlf(I,J) * CS%PhiC(7,i,j))) + &
             ((u_shlf(I-1,J) * CS%PhiC(5,i,j)) + &
              (u_shlf(I,J-1) * CS%PhiC(3,i,j)))

        vx = ((v_shlf(I-1,J-1) * CS%PhiC(1,i,j)) + &
              (v_shlf(I,J) * CS%PhiC(7,i,j))) + &
             ((v_shlf(I-1,J) * CS%PhiC(5,i,j)) + &
              (v_shlf(I,J-1) * CS%PhiC(3,i,j)))

        uy = ((u_shlf(I-1,J-1) * CS%PhiC(2,i,j)) + &
              (u_shlf(I,J) * CS%PhiC(8,i,j))) + &
             ((u_shlf(I-1,J) * CS%PhiC(6,i,j)) + &
              (u_shlf(I,J-1) * CS%PhiC(4,i,j)))

        vy = ((v_shlf(I-1,J-1) * CS%PhiC(2,i,j)) + &
              (v_shlf(I,J) * CS%PhiC(8,i,j))) + &
             ((v_shlf(I-1,J) * CS%PhiC(6,i,j)) + &
              (v_shlf(I,J-1) * CS%PhiC(4,i,j)))

        CS%ice_visc(i,j,1) = (G%areaT(i,j) * max(ISS%h_shelf(i,j),CS%min_h_shelf)) * &
          max(0.5 * Visc_coef * &
          (US%s_to_T**2 * (((ux**2) + (vy**2)) + ((ux*vy) + 0.25*((uy+vx)**2)) + eps_min**2))**((1.-n_g)/(2.*n_g)) * &
          (US%Pa_to_RL2_T2*US%s_to_T),CS%min_ice_visc)
      elseif (model_qp4) then
        !calculate viscosity at 4 quadrature points per cell

        Visc_coef = (CS%AGlen_visc(i,j))**(-1./n_g)

        do iq=1,2 ; do jq=1,2

          ux = ((u_shlf(I-1,J-1) * CS%Phi(1,2*(jq-1)+iq,i,j)) + &
                (u_shlf(I,J) * CS%Phi(7,2*(jq-1)+iq,i,j))) + &
               ((u_shlf(I,J-1) * CS%Phi(3,2*(jq-1)+iq,i,j)) + &
                (u_shlf(I-1,J) * CS%Phi(5,2*(jq-1)+iq,i,j)))

          vx = ((v_shlf(I-1,J-1) * CS%Phi(1,2*(jq-1)+iq,i,j)) + &
                (v_shlf(I,J) * CS%Phi(7,2*(jq-1)+iq,i,j))) + &
               ((v_shlf(I,J-1) * CS%Phi(3,2*(jq-1)+iq,i,j)) + &
                (v_shlf(I-1,J) * CS%Phi(5,2*(jq-1)+iq,i,j)))

          uy = ((u_shlf(I-1,J-1) * CS%Phi(2,2*(jq-1)+iq,i,j)) + &
                (u_shlf(I,J) * CS%Phi(8,2*(jq-1)+iq,i,j))) + &
               ((u_shlf(I,J-1) * CS%Phi(4,2*(jq-1)+iq,i,j)) + &
                (u_shlf(I-1,J) * CS%Phi(6,2*(jq-1)+iq,i,j)))

          vy = ((v_shlf(I-1,J-1) * CS%Phi(2,2*(jq-1)+iq,i,j)) + &
                (v_shlf(I,J) * CS%Phi(8,2*(jq-1)+iq,i,j))) + &
               ((v_shlf(I,J-1) * CS%Phi(4,2*(jq-1)+iq,i,j)) + &
                (v_shlf(I-1,J) * CS%Phi(6,2*(jq-1)+iq,i,j)))

          CS%ice_visc(i,j,2*(jq-1)+iq) = (G%areaT(i,j) * max(ISS%h_shelf(i,j),CS%min_h_shelf)) * &
            max(0.5 * Visc_coef * &
            (US%s_to_T**2 * (((ux**2) + (vy**2)) + ((ux*vy) + 0.25*((uy+vx)**2)) + eps_min**2))**((1.-n_g)/(2.*n_g)) * &
            (US%Pa_to_RL2_T2*US%s_to_T),CS%min_ice_visc)
        enddo; enddo
      endif
    endif
  enddo ; enddo

end subroutine calc_shelf_visc


!> Update basal shear
subroutine calc_shelf_taub(CS, ISS, G, US, u_shlf, v_shlf)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                          intent(inout) :: u_shlf !< The zonal ice shelf velocity [L T-1 ~> m s-1].
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                          intent(inout) :: v_shlf !< The meridional ice shelf velocity [L T-1 ~> m s-1].

! also this subroutine updates the nonlinear part of the basal traction

! this may be subject to change later... to make it "hybrid"

  integer :: i, j, iscq, iecq, jscq, jecq, isd, jsd, ied, jed, iegq, jegq
  integer :: giec, gjec, gisc, gjsc, isc, jsc, iec, jec, is, js
  real :: umid, vmid, unorm, eps_min ! Velocities [L T-1 ~> m s-1]
  real :: alpha !Coulomb coefficient [nondim]
  real :: Hf !"floatation thickness" for Coulomb friction [Z ~> m]
  real :: fN !Effective pressure (ice pressure - ocean pressure) for Coulomb friction [Pa]
  real :: fB !for Coulomb Friction [(T L-1)^CS%CF_PostPeak ~> (s m-1)^CS%CF_PostPeak]
  real :: fN_scale !To convert effective pressure to mks units during Coulomb friction [Pa T2 R-1 L-2 ~> 1]

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
  iegq = G%iegB ; jegq = G%jegB
  gisc = G%domain%nihalo+1 ; gjsc = G%domain%njhalo+1
  giec = G%domain%niglobal+gisc ; gjec = G%domain%njglobal+gjsc
  is = iscq - 1; js = jscq - 1

  eps_min = CS%eps_glen_min

  if (CS%CoulombFriction) then
    if (CS%CF_PostPeak/=1.0) THEN
      alpha = (CS%CF_PostPeak-1.0)**(CS%CF_PostPeak-1.0) / CS%CF_PostPeak**CS%CF_PostPeak ![nondim]
    else
      alpha = 1.0
    endif
    fN_scale = US%R_to_kg_m3 * US%L_T_to_m_s**2
  endif

  do j=jsd+1,jed
    do i=isd+1,ied
      if ((ISS%hmask(i,j) == 1) .OR. (ISS%hmask(i,j) == 3)) then
        umid = ((u_shlf(I,J) + u_shlf(I-1,J-1)) + (u_shlf(I,J-1) + u_shlf(I-1,J))) * 0.25
        vmid = ((v_shlf(I,J) + v_shlf(I-1,J-1)) + (v_shlf(I,J-1) + v_shlf(I-1,J))) * 0.25
        unorm = US%L_T_to_m_s * sqrt( ((umid**2) + (vmid**2)) + (eps_min**2 * (G%dxT(i,j)**2 + G%dyT(i,j)**2)) )

        !Coulomb friction (Schoof 2005, Gagliardini et al 2007)
        if (CS%CoulombFriction) then
          !Effective pressure
          Hf = max((CS%density_ocean_avg/CS%density_ice) * CS%bed_elev(i,j), 0.0)
          fN = max(fN_scale*((CS%density_ice * CS%g_Earth) * (max(ISS%h_shelf(i,j),CS%min_h_shelf) - Hf)),CS%CF_MinN)
          fB = alpha * (CS%C_basal_friction(i,j) / (CS%CF_Max * fN))**(CS%CF_PostPeak/CS%n_basal_fric)

          CS%basal_traction(i,j) = ((G%areaT(i,j) * CS%C_basal_friction(i,j)) * &
            (unorm**(CS%n_basal_fric-1.0) / (1.0 + fB * unorm**CS%CF_PostPeak)**(CS%n_basal_fric))) * &
            (US%Pa_to_RLZ_T2*US%L_T_to_m_s)
        else
          !linear (CS%n_basal_fric=1) or "Weertman"/power-law (CS%n_basal_fric /= 1)
          CS%basal_traction(i,j) = ((G%areaT(i,j) * CS%C_basal_friction(i,j)) * (unorm**(CS%n_basal_fric-1))) * &
                                   (US%Pa_to_RLZ_T2*US%L_T_to_m_s)
        endif

        CS%basal_traction(i,j)=max(CS%basal_traction(i,j), CS%min_basal_traction * G%areaT(i,j))
      endif
    enddo
  enddo

end subroutine calc_shelf_taub

subroutine update_OD_ffrac(CS, G, US, ocean_mass, find_avg)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type), intent(in)     :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: ocean_mass !< The mass per unit area of the ocean [R Z ~> kg m-2].
  logical,                intent(in)    :: find_avg !< If true, find the average of OD and ffrac, and
                                              !! reset the underlying running sums to 0.

  integer :: isc, iec, jsc, jec, i, j
  real    :: I_rho_ocean ! A typical specific volume of the ocean [R-1 ~> m3 kg-1]
  real    :: I_counter

  I_rho_ocean = 1.0 / CS%density_ocean_avg

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  do j=jsc,jec ; do i=isc,iec
    CS%OD_rt(i,j) = CS%OD_rt(i,j) + ocean_mass(i,j)*I_rho_ocean
    if (ocean_mass(i,j)*I_rho_ocean > CS%thresh_float_col_depth) then
      CS%ground_frac_rt(i,j) = CS%ground_frac_rt(i,j) + 1.0
    endif
  enddo ; enddo
  CS%OD_rt_counter = CS%OD_rt_counter + 1

  if (find_avg) then
    I_counter = 1.0 / real(CS%OD_rt_counter)
    do j=jsc,jec ; do i=isc,iec
      CS%ground_frac(i,j) = 1.0 - (CS%ground_frac_rt(i,j) * I_counter)
      CS%OD_av(i,j) = CS%OD_rt(i,j) * I_counter

      CS%OD_rt(i,j) = 0.0 ; CS%ground_frac_rt(i,j) = 0.0; CS%OD_rt_counter = 0
    enddo ; enddo

    call pass_var(CS%ground_frac, G%domain, complete=.false.)
    call pass_var(CS%OD_av, G%domain, complete=.true.)
  endif

end subroutine update_OD_ffrac

subroutine update_OD_ffrac_uncoupled(CS, G, h_shelf)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_shelf !< the thickness of the ice shelf [Z ~> m].

  integer :: i, j, isd, ied, jsd, jed
  real    :: rhoi_rhow, OD

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j=jsd,jed
    do i=isd,ied
      OD = CS%bed_elev(i,j) - rhoi_rhow * max(h_shelf(i,j),CS%min_h_shelf)
      if (OD >= 0) then
    ! ice thickness does not take up whole ocean column -> floating
        CS%OD_av(i,j) = OD
        CS%ground_frac(i,j) = 0.
      else
        CS%OD_av(i,j) = 0.
        CS%ground_frac(i,j) = 1.
      endif
    enddo
  enddo

end subroutine update_OD_ffrac_uncoupled

subroutine change_in_draft(CS, G, h_shelf0, h_shelf1, ddraft)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_shelf0 !< the previous thickness of the ice shelf [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_shelf1 !< the current thickness of the ice shelf [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout)    :: ddraft !< the change in shelf draft thickness
  real :: b0,b1
  integer :: i, j, isc, iec, jsc, jec
  real    :: rhoi_rhow, OD

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ddraft = 0.0

  do j=jsc,jec
    do i=isc,iec

      b0=0.0; b1=0.0

      if (h_shelf0(i,j)>0.0) then
        OD = CS%bed_elev(i,j) - rhoi_rhow * h_shelf0(i,j)
        if (OD >= 0) then
          !floating
          b0 = rhoi_rhow * h_shelf0(i,j)
        else
          b0 = CS%bed_elev(i,j)
        endif
      endif

      if (h_shelf1(i,j)>0.0) then
        OD = CS%bed_elev(i,j) - rhoi_rhow * h_shelf1(i,j)
        if (OD >= 0) then
          !floating
          b1 = rhoi_rhow * h_shelf1(i,j)
        else
          b1 = CS%bed_elev(i,j)
        endif
      endif

      ddraft(i,j) = b1-b0
    enddo
  enddo
end subroutine change_in_draft

!> This subroutine calculates the gradients of bilinear basis elements that
!! that are centered at the vertices of the cell.  Values are calculated at
!! points of gaussian quadrature.
subroutine bilinear_shape_functions (X, Y, Phi, area)
  real, dimension(4),   intent(in)    :: X   !< The x-positions of the vertices of the quadrilateral [L ~> m].
  real, dimension(4),   intent(in)    :: Y   !< The y-positions of the vertices of the quadrilateral [L ~> m].
  real, dimension(8,4), intent(inout) :: Phi !< The gradients of bilinear basis elements at Gaussian
                                             !! quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real,                 intent(out)   :: area !< The quadrilateral cell area [L2 ~> m2].

! X and Y must be passed in the form
    !  3 - 4
    !  |   |
    !  1 - 2

! this subroutine calculates the gradients of bilinear basis elements that
! that are centered at the vertices of the cell. values are calculated at
! points of gaussian quadrature. (in 1D: .5 * (1 +/- sqrt(1/3)) for [0,1])
!     (ordered in same way as vertices)
!
! Phi(2*i-1,j) gives d(Phi_i)/dx at quadrature point j
! Phi(2*i,j) gives d(Phi_i)/dy at quadrature point j
! Phi_i is equal to 1 at vertex i, and 0 at vertex k /= i, and bilinear
!
! This should be a one-off; once per nonlinear solve? once per lifetime?
! ... will all cells have the same shape and dimension?

  real, dimension(4) :: xquad, yquad ! [nondim]
  real :: a,b,c,d  ! Various lengths [L ~> m]
  real :: xexp, yexp ! [nondim]
  integer :: node, qpoint, xnode, ynode

  xquad(1:3:2) = .5 * (1-sqrt(1./3)) ; yquad(1:2) = .5 * (1-sqrt(1./3))
  xquad(2:4:2) = .5 * (1+sqrt(1./3)) ; yquad(3:4) = .5 * (1+sqrt(1./3))

  do qpoint=1,4

    a = ((-X(1)*(1-yquad(qpoint)))+(X(4)*yquad(qpoint))) + ((X(2)*(1-yquad(qpoint)))-(X(3)*yquad(qpoint))) !d(x)/d(x*)
    b = ((-Y(1)*(1-yquad(qpoint)))+(Y(4)*yquad(qpoint))) + ((Y(2)*(1-yquad(qpoint)))-(Y(3)*yquad(qpoint))) !d(y)/d(x*)
    c = ((-X(1)*(1-xquad(qpoint)))+(X(4)*xquad(qpoint))) + ((-X(2)*xquad(qpoint))+(X(3)*(1-xquad(qpoint))))!d(x)/d(y*)
    d = ((-Y(1)*(1-xquad(qpoint)))+(Y(4)*xquad(qpoint))) + ((-Y(2)*xquad(qpoint))+(Y(3)*(1-xquad(qpoint))))!d(y)/d(y*)

    do node=1,4

      xnode = 2-mod(node,2) ; ynode = ceiling(REAL(node)/2)

      if (ynode == 1) then
        yexp = 1-yquad(qpoint)
      else
        yexp = yquad(qpoint)
      endif

      if (1 == xnode) then
        xexp = 1-xquad(qpoint)
      else
        xexp = xquad(qpoint)
      endif

      Phi(2*node-1,qpoint) = ( d * (2 * xnode - 3) * yexp - b * (2 * ynode - 3) * xexp) / ((a*d)-(b*c))
      Phi(2*node,qpoint)   = (-c * (2 * xnode - 3) * yexp + a * (2 * ynode - 3) * xexp) / ((a*d)-(b*c))

    enddo
  enddo

  area = quad_area(X, Y)

end subroutine bilinear_shape_functions

!> This subroutine calculates the gradients of bilinear basis elements that are centered at the
!! vertices of the cell using a locally orthogoal MOM6 grid.  Values are calculated at
!! points of gaussian quadrature.
subroutine bilinear_shape_fn_grid(G, i, j, Phi)
  type(ocean_grid_type), intent(in)    :: G  !< The grid structure used by the ice shelf.
  integer,               intent(in)    :: i   !< The i-index in the grid to work on.
  integer,               intent(in)    :: j   !< The j-index in the grid to work on.
  real, dimension(8,4),  intent(inout) :: Phi !< The gradients of bilinear basis elements at Gaussian
                                              !! quadrature points surrounding the cell vertices [L-1 ~> m-1].

! This subroutine calculates the gradients of bilinear basis elements that
! that are centered at the vertices of the cell.  The values are calculated at
! points of gaussian quadrature. (in 1D: .5 * (1 +/- sqrt(1/3)) for [0,1])
!     (ordered in same way as vertices)
!
! Phi(2*i-1,j) gives d(Phi_i)/dx at quadrature point j
! Phi(2*i,j) gives d(Phi_i)/dy at quadrature point j
! Phi_i is equal to 1 at vertex i, and 0 at vertex k /= i, and bilinear
!
! This should be a one-off; once per nonlinear solve? once per lifetime?

  real, dimension(4) :: xquad, yquad ! [nondim]
  real :: a, d       ! Interpolated grid spacings [L ~> m]
  real :: xexp, yexp ! [nondim]
  integer :: node, qpoint, xnode, ynode

  xquad(1:3:2) = .5 * (1-sqrt(1./3)) ; yquad(1:2) = .5 * (1-sqrt(1./3))
  xquad(2:4:2) = .5 * (1+sqrt(1./3)) ; yquad(3:4) = .5 * (1+sqrt(1./3))

  do qpoint=1,4
    if (J>1) then
      a = (G%dxCv(i,J-1) * (1-yquad(qpoint))) + (G%dxCv(i,J) * yquad(qpoint)) ! d(x)/d(x*)
    else
      a = G%dxCv(i,J) !* yquad(qpoint) ! d(x)/d(x*)
    endif
    if (I>1) then
      d = (G%dyCu(I-1,j) * (1-xquad(qpoint))) + (G%dyCu(I,j) * xquad(qpoint)) ! d(y)/d(y*)
    else
      d = G%dyCu(I,j) !* xquad(qpoint)
    endif
!    a = G%dxCv(i,J-1) * (1-yquad(qpoint)) + G%dxCv(i,J) * yquad(qpoint) ! d(x)/d(x*)
!    d = G%dyCu(I-1,j) * (1-xquad(qpoint)) + G%dyCu(I,j) * xquad(qpoint) ! d(y)/d(y*)

    do node=1,4
      xnode = 2-mod(node,2) ; ynode = ceiling(REAL(node)/2)

      if (ynode == 1) then
        yexp = 1-yquad(qpoint)
      else
        yexp = yquad(qpoint)
      endif

      if (1 == xnode) then
        xexp = 1-xquad(qpoint)
      else
        xexp = xquad(qpoint)
      endif

      Phi(2*node-1,qpoint) = ( (d * (2 * xnode - 3)) * yexp ) / (a*d)
      Phi(2*node,qpoint)   = ( (a * (2 * ynode - 3)) * xexp ) / (a*d)

    enddo
  enddo

end subroutine bilinear_shape_fn_grid

!> This subroutine calculates the gradients of bilinear basis elements that are centered at the
!! vertices of the cell using a locally orthogoal MOM6 grid.  Values are calculated at
!! a sinlge cell-centered quadrature point, which should match the grid cell h-point
subroutine bilinear_shape_fn_grid_1qp(G, i, j, Phi)
  type(ocean_grid_type), intent(in)    :: G  !< The grid structure used by the ice shelf.
  integer,               intent(in)    :: i   !< The i-index in the grid to work on.
  integer,               intent(in)    :: j   !< The j-index in the grid to work on.
  real, dimension(8),    intent(inout) :: Phi !< The gradients of bilinear basis elements at Gaussian
                                              !! quadrature points surrounding the cell vertices [L-1 ~> m-1].

! This subroutine calculates the gradients of bilinear basis elements that
! that are centered at the vertices of the cell.  The values are calculated at
! a cell-cented point of gaussian quadrature. (in 1D: .5 for [0,1])
!     (ordered in same way as vertices)
!
! Phi(2*i-1) gives d(Phi_i)/dx at the quadrature point
! Phi(2*i) gives d(Phi_i)/dy at the quadrature point
! Phi_i is equal to 1 at vertex i, and 0 at vertex k /= i, and bilinear

  real :: a, d       ! Interpolated grid spacings [L ~> m]
  real :: xexp=0.5, yexp=0.5 ! [nondim]
  integer :: node, qpoint, xnode, ynode

    ! d(x)/d(x*)
    if (J>1) then
      a = 0.5 * (G%dxCv(i,J-1) + G%dxCv(i,J))
    else
      a = G%dxCv(i,J)
    endif

    ! d(y)/d(y*)
    if (I>1) then
      d = 0.5 * (G%dyCu(I-1,j) + G%dyCu(I,j))
    else
      d = G%dyCu(I,j)
    endif

    do node=1,4
      xnode = 2-mod(node,2) ; ynode = ceiling(REAL(node)/2)
      Phi(2*node-1) = ( (d * (2 * xnode - 3)) * yexp ) / (a*d)
      Phi(2*node)   = ( (a * (2 * ynode - 3)) * xexp ) / (a*d)
    enddo
end subroutine bilinear_shape_fn_grid_1qp


subroutine bilinear_shape_functions_subgrid(Phisub, nsub)
  integer, intent(in)    :: nsub   !< The number of subgridscale quadrature locations in each direction
  real, dimension(2,2,nsub,nsub,2,2), &
           intent(inout) :: Phisub !< Quadrature structure weights at subgridscale
                                   !! locations for finite element calculations [nondim]

  ! this subroutine is a helper for interpolation of floatation condition
  ! for the purposes of evaluating the terms \int (u,v) \phi_i dx dy in a cell that is
  !     in partial floatation
  ! the array Phisub contains the values of \phi_i (where i is a node of the cell)
  !     at quad point j
  ! i think this general approach may not work for nonrectangular elements...
  !

  ! Phisub(q1,q2,i,j,k,l)
  !  q1: quad point x-index
  !  q2: quad point y-index
  !  i: subgrid index in x-direction
  !  j: subgrid index in y-direction
  !  k: basis function x-index
  !  l: basis function y-index

  ! e.g. k=1,l=1 => node 1
  !      q1=2,q2=1 => quad point 2

    !  3 - 4
    !  |   |
    !  1 - 2

  integer :: i, j, qx, qy
  real,dimension(2)    :: xquad
  real                 :: x0, y0, x, y, fracx

  xquad(1) = .5 * (1-sqrt(1./3)) ; xquad(2) = .5 * (1+sqrt(1./3))
  fracx = 1.0/real(nsub)

  do j=1,nsub ; do i=1,nsub
    x0 = (i-1) * fracx ; y0 = (j-1) * fracx
    do qy=1,2 ; do qx=1,2
      x = x0 + fracx*xquad(qx)
      y = y0 + fracx*xquad(qy)
      Phisub(qx,qy,i,j,1,1) = (1.0-x) * (1.0-y)
      Phisub(qx,qy,i,j,1,2) = (1.0-x) * y
      Phisub(qx,qy,i,j,2,1) = x * (1.0-y)
      Phisub(qx,qy,i,j,2,2) = x * y
    enddo ; enddo
  enddo ; enddo

end subroutine bilinear_shape_functions_subgrid


subroutine update_velocity_masks(CS, G, hmask, umask, vmask, u_face_mask, v_face_mask)
  type(ice_shelf_dyn_CS),intent(in)    :: CS !< A pointer to the ice shelf dynamics control structure
  type(ocean_grid_type), intent(inout) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: umask !< A coded mask indicating the nature of the
                                             !! zonal flow at the corner point
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: vmask !< A coded mask indicating the nature of the
                                             !! meridional flow at the corner point
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: u_face_mask !< A coded mask for velocities at the C-grid u-face
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: v_face_mask !< A coded mask for velocities at the C-grid v-face
  ! sets masks for velocity solve
  ! ignores the fact that their might be ice-free cells - this only considers the computational boundary

  ! !!!IMPORTANT!!! relies on thickness mask - assumed that this is called after hmask has been updated & halo-updated

  integer :: i, j, k, iscq, iecq, jscq, jecq, isd, jsd, is, js, iegq, jegq
  integer :: giec, gjec, gisc, gjsc, isc, jsc, iec, jec

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  isd = G%isd ; jsd = G%jsd
  iegq = G%iegB ; jegq = G%jegB
  gisc = G%Domain%nihalo ; gjsc = G%Domain%njhalo
  giec = G%Domain%niglobal+gisc ; gjec = G%Domain%njglobal+gjsc

  umask(:,:) = 0 ; vmask(:,:) = 0
  u_face_mask(:,:) = 0 ; v_face_mask(:,:) = 0

  if (G%symmetric) then
    is = isd ; js = jsd
  else
    is = isd+1 ; js = jsd+1
  endif

  do j=js,G%jed; do i=is,G%ied
    if (hmask(i,j) == 1 .or. hmask(i,j)==3) then
      umask(I-1:I,J-1:J)=1
      vmask(I-1:I,J-1:J)=1
    endif
  enddo; enddo

  do j=js,G%jed
    do i=is,G%ied

      if ((hmask(i,j) == 1) .OR. (hmask(i,j) == 3)) then

        do k=0,1

          select case (int(CS%u_face_mask_bdry(I-1+k,j)))
            case (5)
              umask(I-1+k,J-1:J) = 3.
              u_face_mask(I-1+k,j) = 5.
            case (3)
              umask(I-1+k,J-1:J) = 3.
              vmask(I-1+k,J-1:J) = 3.
              u_face_mask(I-1+k,j) = 3.
            case (6)
              vmask(I-1+k,J-1:J) = 3.
              u_face_mask(I-1+k,j) = 6.
            case (2)
              u_face_mask(I-1+k,j) = 2.
            case (4)
              umask(I-1+k,J-1:J) = 0.
              u_face_mask(I-1+k,j) = 4.
            case (0)
              umask(I-1+k,J-1:J) = 0.
              u_face_mask(I-1+k,j) = 0.
            case (1)  ! stress free x-boundary
              umask(I-1+k,J-1:J) = 0.
            case default
              umask(I-1+k,J-1) = max(1. , umask(I-1+k,J-1))
              umask(I-1+k,J)   = max(1. , umask(I-1+k,J))
          end select
        enddo

        do k=0,1

          select case (int(CS%v_face_mask_bdry(i,J-1+k)))
            case (5)
              vmask(I-1:I,J-1+k) = 3.
              v_face_mask(i,J-1+k) = 5.
            case (3)
              vmask(I-1:I,J-1+k) = 3.
              umask(I-1:I,J-1+k) = 3.
              v_face_mask(i,J-1+k) = 3.
            case (6)
              umask(I-1:I,J-1+k) = 3.
              v_face_mask(i,J-1+k) = 6.
            case (2)
              v_face_mask(i,J-1+k) = 2.
            case (4)
              vmask(I-1:I,J-1+k) = 0.
              v_face_mask(i,J-1+k) = 4.
            case (0)
              vmask(I-1:I,J-1+k) = 0.
              v_face_mask(i,J-1+k) = 0.
            case (1) ! stress free y-boundary
              vmask(I-1:I,J-1+k) = 0.
            case default
              vmask(I-1,J-1+k) = max(1. , vmask(I-1,J-1+k))
              vmask(I,J-1+k)   = max(1. , vmask(I,J-1+k))
          end select
        enddo


        if (i < G%ied) then
          if ((hmask(i+1,j) == 0) .OR. (hmask(i+1,j) == 2)) then
            ! east boundary or adjacent to unfilled cell
            u_face_mask(I,j) = 2.
          endif
        endif

        if (i > G%isd) then
          if ((hmask(i-1,j) == 0) .OR. (hmask(i-1,j) == 2)) then
            !adjacent to unfilled cell
            u_face_mask(I-1,j) = 2.
          endif
        endif

        if (j > G%jsd) then
          if ((hmask(i,j-1) == 0) .OR. (hmask(i,j-1) == 2)) then
            !adjacent to unfilled cell
            v_face_mask(i,J-1) = 2.
          endif
        endif

        if (j < G%jed) then
          if ((hmask(i,j+1) == 0) .OR. (hmask(i,j+1) == 2)) then
            !adjacent to unfilled cell
            v_face_mask(i,j) = 2.
          endif
        endif


      endif

    enddo
  enddo

  ! note: if the grid is nonsymmetric, there is a part that will not be transferred with a halo update
  ! so this subroutine must update its own symmetric part of the halo

  call pass_vector(u_face_mask, v_face_mask, G%domain, TO_ALL, CGRID_NE)
  call pass_vector(umask, vmask, G%domain, TO_ALL, BGRID_NE)

end subroutine update_velocity_masks

!> Interpolate the ice shelf thickness from tracer point to nodal points,
!! subject to a mask.
subroutine interpolate_H_to_B(G, h_shelf, hmask, H_node, min_h_shelf)
  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: h_shelf !< The ice shelf thickness at tracer points [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: H_node !< The ice shelf thickness at nodal (corner)
                                             !! points [Z ~> m].
  real, intent(in) :: min_h_shelf !< The minimum ice thickness used during ice dynamics [L ~> m].

  integer :: i, j, isc, iec, jsc, jec, num_h, k, l, ic, jc
  real    :: h_arr(2,2)

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  H_node(:,:) = 0.0

  ! H_node is node-centered; average over all cells that share that node
  ! if no (active) cells share the node then its value there is irrelevant

  do j=jsc-1,jec
    do i=isc-1,iec
      num_h = 0
      do l=1,2; jc=j-1+l; do k=1,2; ic=i-1+k
        if (hmask(ic,jc) == 1.0 .or. hmask(ic,jc) == 3.0) then
          h_arr(k,l)=max(h_shelf(ic,jc),min_h_shelf)
          num_h = num_h + 1
        else
          h_arr(k,l)=0.0
        endif
        if (num_h > 0) then
          H_node(i,j) = ((h_arr(1,1)+h_arr(2,2))+(h_arr(1,2)+h_arr(2,1))) / num_h
        endif
      enddo; enddo
    enddo
  enddo

  call pass_var(H_node, G%domain,position=CORNER)

end subroutine interpolate_H_to_B

!> Deallocates all memory associated with the ice shelf dynamics module
subroutine ice_shelf_dyn_end(CS)
  type(ice_shelf_dyn_CS), pointer   :: CS !< A pointer to the ice shelf dynamics control structure

  if (.not.associated(CS)) return

  deallocate(CS%u_shelf, CS%v_shelf)
  deallocate(CS%taudx_shelf, CS%taudy_shelf)
  deallocate(CS%t_shelf, CS%tmask)
  deallocate(CS%u_bdry_val, CS%v_bdry_val)
  deallocate(CS%u_face_mask, CS%v_face_mask)
  deallocate(CS%umask, CS%vmask)
  deallocate(CS%u_face_mask_bdry, CS%v_face_mask_bdry)
  deallocate(CS%h_bdry_val)
  deallocate(CS%float_cond)

  deallocate(CS%ice_visc, CS%AGlen_visc)
  deallocate(CS%basal_traction,CS%C_basal_friction)
  deallocate(CS%OD_rt, CS%OD_av)
  deallocate(CS%t_bdry_val, CS%bed_elev)
  deallocate(CS%ground_frac, CS%ground_frac_rt)

  deallocate(CS)

end subroutine ice_shelf_dyn_end


!> This subroutine updates the vertically averaged ice shelf temperature.
subroutine ice_shelf_temp(CS, ISS, G, US, time_step, melt_rate, Time)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: melt_rate !< basal melt rate [R Z T-1 ~> kg m-2 s-1]
  type(time_type),        intent(in)    :: Time !< The current model time

!    This subroutine takes the velocity (on the Bgrid) and timesteps
!      (HT)_t = - div (uHT) + (adot Tsurf -bdot Tbot) once and then calculates T=HT/H
!
!    The flux overflows are included here. That is because they will be used to advect 3D scalars
!    into partial cells

  real, dimension(SZDI_(G),SZDJ_(G))   :: th_after_uflux, th_after_vflux, TH ! Integrated temperatures [C Z ~> degC m]
  integer                           :: isd, ied, jsd, jed, i, j, isc, iec, jsc, jec
  real :: Tsurf ! Surface air temperature [C ~> degC].  This is hard coded but should be an input argument.
  real :: adot  ! A surface heat exchange coefficient [R Z T-1 ~> kg m-2 s-1].


  ! For now adot and Tsurf are defined here adot=surf acc 0.1m/yr, Tsurf=-20oC, vary them later
  adot = (0.1/(365.0*86400.0))*US%m_to_Z*US%T_to_s * CS%density_ice
  Tsurf = -20.0*US%degC_to_C

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  th_after_uflux(:,:) = 0.0
  th_after_vflux(:,:) = 0.0

  do j=jsd,jed ; do i=isd,ied
!    if (ISS%hmask(i,j) > 1) then
    if ((ISS%hmask(i,j) == 3) .or. (ISS%hmask(i,j) == -2)) then
      CS%t_shelf(i,j) = CS%t_bdry_val(i,j)
    endif
  enddo ; enddo

  do j=jsd,jed ; do i=isd,ied
    ! Convert the averge temperature to a depth integrated temperature.
    TH(i,j) = CS%t_shelf(i,j)*ISS%h_shelf(i,j)
  enddo ; enddo


  call ice_shelf_advect_temp_x(CS, G, time_step, ISS%hmask, TH, th_after_uflux)
  call ice_shelf_advect_temp_y(CS, G, time_step, ISS%hmask, th_after_uflux, th_after_vflux)

  do j=jsc,jec ; do i=isc,iec
    ! Convert the integrated temperature back to the average temperature.
!   if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 2)) then
    if (ISS%h_shelf(i,j) > 0.0) then
      CS%t_shelf(i,j) = th_after_vflux(i,j) / ISS%h_shelf(i,j)
    else
      CS%t_shelf(i,j) = CS%T_shelf_missing
    endif
!   endif

    if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 2)) then
      if (ISS%h_shelf(i,j) > 0.0) then
        CS%t_shelf(i,j) = CS%t_shelf(i,j) + &
            time_step*(adot*Tsurf - melt_rate(i,j)*ISS%tfreeze(i,j))/(CS%density_ice*ISS%h_shelf(i,j))
      else
        ! the ice is about to melt away in this case set thickness, area, and mask to zero
        ! NOTE: not mass conservative, should maybe scale salt & heat flux for this cell
        CS%t_shelf(i,j) = CS%T_shelf_missing
        CS%tmask(i,j) = 0.0
      endif
    elseif (ISS%hmask(i,j) == 0) then
      CS%t_shelf(i,j) = CS%T_shelf_missing
    elseif ((ISS%hmask(i,j) == 3) .or. (ISS%hmask(i,j) == -2)) then
      CS%t_shelf(i,j) = CS%t_bdry_val(i,j)
    endif
  enddo ; enddo

  call pass_var(CS%t_shelf, G%domain, complete=.false.)
  call pass_var(CS%tmask, G%domain, complete=.true.)

  if (CS%debug) then
    call hchksum(CS%t_shelf, "temp after front", G%HI, haloshift=3, unscale=US%C_to_degC)
  endif

end subroutine ice_shelf_temp


subroutine ice_shelf_advect_temp_x(CS, G, time_step, hmask, h0, h_after_uflux)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h0 !< The initial ice shelf thicknesses times temperature [C Z ~> degC m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_uflux !< The ice shelf thicknesses times temperature after
                                              !! the zonal mass fluxes [C Z ~> degC m]

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization

  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  integer :: i_off, j_off
  logical :: at_east_bdry, at_west_bdry
  real, dimension(-2:2) :: stencil ! A copy of the neighboring thicknesses times temperatures [C Z ~> degC m]
  real :: u_face     ! Zonal velocity at a face, positive if out [L T-1 ~> m s-1]
  real :: flux_diff  ! The difference in fluxes [C Z ~> degC m]
  real :: phi        ! A limiting ratio [nondim]

  is = G%isc-2 ; ie = G%iec+2 ; js = G%jsc ; je = G%jec ; isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset

  do j=jsd+1,jed-1
    if (((j+j_off) <= G%domain%njglobal+G%domain%njhalo) .AND. &
        ((j+j_off) >= G%domain%njhalo+1)) then ! based on mehmet's code - only if btw north & south boundaries

      stencil(:) = 0.0 ! This is probably unnecessary, as the code is written
!     if (i+i_off == G%domain%nihalo+G%domain%nihalo)
      do i=is,ie

        if (((i+i_off) <= G%domain%niglobal+G%domain%nihalo) .AND. &
             ((i+i_off) >= G%domain%nihalo+1)) then

          if (i+i_off == G%domain%nihalo+1) then
            at_west_bdry=.true.
          else
            at_west_bdry=.false.
          endif

          if (i+i_off == G%domain%niglobal+G%domain%nihalo) then
            at_east_bdry=.true.
          else
            at_east_bdry=.false.
          endif

          if (hmask(i,j) == 1) then

            h_after_uflux(i,j) = h0(i,j)

            stencil(:) = h0(i-2:i+2,j)  ! fine as long has nx_halo >= 2

            flux_diff = 0

            ! 1ST DO LEFT FACE

            if (CS%u_face_mask(I-1,j) == 4.) then

              flux_diff = flux_diff + G%dyCu(I-1,j) * time_step * CS%u_flux_bdry_val(I-1,j) * &
                               CS%t_bdry_val(i-1,j) / G%areaT(i,j)
            else

              ! get u-velocity at center of left face
              u_face = 0.5 * (CS%u_shelf(I-1,J-1) + CS%u_shelf(I-1,J))

              if (u_face > 0) then !flux is into cell - we need info from h(i-2), h(i-1) if available

              ! i may not cover all the cases.. but i cover the realistic ones

                if (at_west_bdry .AND. (hmask(i-1,j) == 3)) then ! at western bdry but there is a
                              ! thickness bdry condition, and the stencil contains it
                  flux_diff = flux_diff + ABS(u_face) * G%dyCu(I-1,j) * time_step * stencil(-1) / G%areaT(i,j)

                elseif (hmask(i-1,j) * hmask(i-2,j) == 1) then  ! h(i-2) and h(i-1) are valid
                  phi = slope_limiter(stencil(-1)-stencil(-2), stencil(0)-stencil(-1))
                  flux_diff = flux_diff + ((ABS(u_face) * G%dyCu(I-1,j)* time_step / G%areaT(i,j)) * &
                           (stencil(-1) - (phi * (stencil(-1)-stencil(0))/2)))

                else                            ! h(i-1) is valid
                                    ! (o.w. flux would most likely be out of cell)
                                    !  but h(i-2) is not

                  flux_diff = flux_diff + ABS(u_face) * G%dyCu(I-1,j) * time_step / G%areaT(i,j) * stencil(-1)

                endif

              elseif (u_face < 0) then !flux is out of cell - we need info from h(i-1), h(i+1) if available
                if (hmask(i-1,j) * hmask(i+1,j) == 1) then         ! h(i-1) and h(i+1) are both valid
                  phi = slope_limiter(stencil(0)-stencil(1), stencil(-1)-stencil(0))
                  flux_diff = flux_diff - ((ABS(u_face) * G%dyCu(I-1,j) * time_step / G%areaT(i,j)) * &
                             (stencil(0) - (phi * (stencil(0)-stencil(-1))/2)))

                else
                  flux_diff = flux_diff - ABS(u_face) * G%dyCu(I-1,j) * time_step / G%areaT(i,j) * stencil(0)
                endif
              endif
            endif

            ! NEXT DO RIGHT FACE

            ! get u-velocity at center of eastern face

            if (CS%u_face_mask(I,j) == 4.) then

              flux_diff = flux_diff + G%dyCu(I,j) * time_step * CS%u_flux_bdry_val(I,j) *&
                               CS%t_bdry_val(i+1,j) / G%areaT(i,j)
            else

              u_face = 0.5 * (CS%u_shelf(I,J-1) + CS%u_shelf(I,J))

              if (u_face < 0) then !flux is into cell - we need info from h(i+2), h(i+1) if available

                if (at_east_bdry .AND. (hmask(i+1,j) == 3)) then ! at eastern bdry but there is a
                                            ! thickness bdry condition, and the stencil contains it

                  flux_diff = flux_diff + ABS(u_face) * G%dyCu(I,j) * time_step * stencil(1) / G%areaT(i,j)

                elseif (hmask(i+1,j) * hmask(i+2,j) == 1) then  ! h(i+2) and h(i+1) are valid

                  phi = slope_limiter(stencil(1)-stencil(2), stencil(0)-stencil(1))
                  flux_diff = flux_diff + ((ABS(u_face) * G%dyCu(I,j) * time_step / G%areaT(i,j)) * &
                      (stencil(1) - (phi * (stencil(1)-stencil(0))/2)))

                else                            ! h(i+1) is valid
                                            ! (o.w. flux would most likely be out of cell)
                                            !  but h(i+2) is not

                  flux_diff = flux_diff + ABS(u_face) * G%dyCu(I,j) * time_step / G%areaT(i,j) * stencil(1)

                endif

              elseif (u_face > 0) then !flux is out of cell - we need info from h(i-1), h(i+1) if available

                if (hmask(i-1,j) * hmask(i+1,j) == 1) then         ! h(i-1) and h(i+1) are both valid

                  phi = slope_limiter(stencil(0)-stencil(-1), stencil(1)-stencil(0))
                  flux_diff = flux_diff - ((ABS(u_face) * G%dyCu(I,j) * time_step / G%areaT(i,j)) * &
                      (stencil(0) - (phi * (stencil(0)-stencil(1))/2)))

                else  ! h(i+1) is valid (o.w. flux would most likely be out of cell) but h(i+2) is not

                  flux_diff = flux_diff - ABS(u_face) * G%dyCu(I,j) * time_step / G%areaT(i,j) * stencil(0)

                endif

              endif

              h_after_uflux(i,j) = h_after_uflux(i,j) + flux_diff

            endif

          endif

        endif

      enddo ! i loop

    endif

  enddo ! j loop

end subroutine ice_shelf_advect_temp_x

subroutine ice_shelf_advect_temp_y(CS, G, time_step, hmask, h_after_uflux, h_after_vflux)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_after_uflux !< The ice shelf thicknesses times temperature after
                                              !! the zonal mass fluxes [C Z ~> degC m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_vflux !< The ice shelf thicknesses times temperature after
                                              !! the meridional mass fluxes [C Z ~> degC m]

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization

  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  integer :: i_off, j_off
  logical :: at_north_bdry, at_south_bdry
  real, dimension(-2:2) :: stencil ! A copy of the neighboring thicknesses times temperatures [C Z ~> degC m]
  real :: v_face     ! Pseudo-meridional velocity at a cell face, positive if out [L T-1 ~> m s-1]
  real :: flux_diff  ! The difference in fluxes [C Z ~> degC m]
  real :: phi

  is = G%isc ; ie = G%iec ; js = G%jsc-1 ; je = G%jec+1 ; isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset

  do i=isd+2,ied-2
    if (((i+i_off) <= G%domain%niglobal+G%domain%nihalo) .AND. &
       ((i+i_off) >= G%domain%nihalo+1)) then  ! based on mehmet's code - only if btw east & west boundaries

      stencil(:) = 0.0 ! This is probably unnecessary, as the code is written

      do j=js,je

        if (((j+j_off) <= G%domain%njglobal+G%domain%njhalo) .AND. &
             ((j+j_off) >= G%domain%njhalo+1)) then

          if (j+j_off == G%domain%njhalo+1) then
            at_south_bdry=.true.
          else
            at_south_bdry=.false.
          endif
          if (j+j_off == G%domain%njglobal+G%domain%njhalo) then
            at_north_bdry=.true.
          else
            at_north_bdry=.false.
          endif

          if (hmask(i,j) == 1) then
            h_after_vflux(i,j) = h_after_uflux(i,j)

            stencil(:) = h_after_uflux(i,j-2:j+2)  ! fine as long has ny_halo >= 2
            flux_diff = 0

            ! 1ST DO south FACE

            if (CS%v_face_mask(i,J-1) == 4.) then

              flux_diff = flux_diff + G%dxCv(i,J-1) * time_step * CS%v_flux_bdry_val(i,J-1) * &
                                 CS%t_bdry_val(i,j-1)/ G%areaT(i,j)
            else

              ! get u-velocity at center of west face
              v_face = 0.5 * (CS%v_shelf(I-1,J-1) + CS%v_shelf(I,J-1))

              if (v_face > 0) then !flux is into cell - we need info from h(j-2), h(j-1) if available

                ! i may not cover all the cases.. but i cover the realistic ones

                if (at_south_bdry .AND. (hmask(i,j-1) == 3)) then ! at western bdry but there is a
                                            ! thickness bdry condition, and the stencil contains it
                  flux_diff = flux_diff + ABS(v_face) * G%dxCv(i,J-1) * time_step * stencil(-1) / G%areaT(i,j)

                elseif (hmask(i,j-1) * hmask(i,j-2) == 1) then  ! h(j-2) and h(j-1) are valid

                  phi = slope_limiter(stencil(-1)-stencil(-2), stencil(0)-stencil(-1))
                  flux_diff = flux_diff + ((ABS(v_face) * G%dxCv(i,J-1) * time_step / G%areaT(i,j)) * &
                      (stencil(-1) - (phi * (stencil(-1)-stencil(0))/2)))

                else     ! h(j-1) is valid
                         ! (o.w. flux would most likely be out of cell)
                         !  but h(j-2) is not
                  flux_diff = flux_diff + ABS(v_face) * G%dxCv(i,J-1) * time_step / G%areaT(i,j) * stencil(-1)
                endif

              elseif (v_face < 0) then !flux is out of cell - we need info from h(j-1), h(j+1) if available

                if (hmask(i,j-1) * hmask(i,j+1) == 1) then  ! h(j-1) and h(j+1) are both valid
                  phi = slope_limiter(stencil(0)-stencil(1), stencil(-1)-stencil(0))
                  flux_diff = flux_diff - ((ABS(v_face) * G%dxCv(i,J-1) * time_step / G%areaT(i,j)) * &
                      (stencil(0) - (phi * (stencil(0)-stencil(-1))/2)))
                else
                  flux_diff = flux_diff - ABS(v_face) * G%dxCv(i,J-1) * time_step / G%areaT(i,j) * stencil(0)
                endif

              endif

            endif

            ! NEXT DO north FACE

            if (CS%v_face_mask(i,J) == 4.) then
              flux_diff = flux_diff + G%dxCv(i,J) * time_step * CS%v_flux_bdry_val(i,J) *&
                               CS%t_bdry_val(i,j+1)/ G%areaT(i,j)
            else

            ! get u-velocity at center of east face
              v_face = 0.5 * (CS%v_shelf(I-1,J) + CS%v_shelf(I,J))

              if (v_face < 0) then !flux is into cell - we need info from h(j+2), h(j+1) if available

                if (at_north_bdry .AND. (hmask(i,j+1) == 3)) then ! at eastern bdry but there is a
                                            ! thickness bdry condition, and the stencil contains it
                  flux_diff = flux_diff + ABS(v_face) * G%dxCv(i,J) * time_step * stencil(1) / G%areaT(i,j)
                elseif (hmask(i,j+1) * hmask(i,j+2) == 1) then  ! h(j+2) and h(j+1) are valid
                  phi = slope_limiter (stencil(1)-stencil(2), stencil(0)-stencil(1))
                  flux_diff = flux_diff + ((ABS(v_face) * G%dxCv(i,J) * time_step / G%areaT(i,j)) * &
                      (stencil(1) - (phi * (stencil(1)-stencil(0))/2)))
                else     ! h(j+1) is valid
                         ! (o.w. flux would most likely be out of cell)
                         !  but h(j+2) is not
                  flux_diff = flux_diff + ABS(v_face) * G%dxCv(i,J) * time_step / G%areaT(i,j) * stencil(1)
                endif

              elseif (v_face > 0) then !flux is out of cell - we need info from h(j-1), h(j+1) if available

                if (hmask(i,j-1) * hmask(i,j+1) == 1) then         ! h(j-1) and h(j+1) are both valid
                  phi = slope_limiter (stencil(0)-stencil(-1), stencil(1)-stencil(0))
                  flux_diff = flux_diff - ((ABS(v_face) * G%dxCv(i,J) * time_step / G%areaT(i,j)) * &
                      (stencil(0) - (phi * (stencil(0)-stencil(1))/2)))
                else   ! h(j+1) is valid
                       ! (o.w. flux would most likely be out of cell)
                       !  but h(j+2) is not
                  flux_diff = flux_diff - ABS(v_face) * G%dxCv(i,J) * time_step / G%areaT(i,j) * stencil(0)
                endif

              endif

            endif

            h_after_vflux(i,j) = h_after_vflux(i,j) + flux_diff
          endif
        endif
      enddo ! j loop
    endif
  enddo ! i loop

end subroutine ice_shelf_advect_temp_y

end module MOM_ice_shelf_dynamics
