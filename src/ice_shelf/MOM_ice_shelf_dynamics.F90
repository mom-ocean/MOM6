!> Implements a crude placeholder for a later implementation of full
!! ice shelf dynamics.
module MOM_ice_shelf_dynamics

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_mediator_init, set_diag_mediator_grid
use MOM_diag_mediator, only : diag_ctrl, time_type, enable_averages, disable_averaging
use MOM_domains, only : MOM_domains_init, clone_MOM_domain
use MOM_domains, only : pass_var, pass_vector, TO_ALL, CGRID_NE, BGRID_NE, CORNER
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_grid, only : MOM_grid_init, ocean_grid_type
use MOM_io, only : file_exists, slasher, MOM_read_data
use MOM_restart, only : register_restart_field, query_initialized
use MOM_restart, only : MOM_restart_CS
use MOM_time_manager, only : time_type, set_time
use MOM_unit_scaling, only : unit_scale_type, unit_scaling_init
!MJH use MOM_ice_shelf_initialize, only : initialize_ice_shelf_boundary
use MOM_ice_shelf_state, only : ice_shelf_state
use MOM_coms, only : reproducing_sum, sum_across_PEs, max_across_PEs, min_across_PEs
use MOM_checksums, only : hchksum, qchksum

implicit none ; private

#include <MOM_memory.h>

public register_ice_shelf_dyn_restarts, initialize_ice_shelf_dyn, update_ice_shelf
public ice_time_step_CFL, ice_shelf_dyn_end
public shelf_advance_front, ice_shelf_min_thickness_calve, calve_to_mask

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

  real, pointer, dimension(:,:) :: u_face_mask => NULL() !< mask for velocity boundary conditions on the C-grid
                                       !! u-face - this is because the FEM cares about FACES THAT GET INTEGRATED OVER,
                                       !! not vertices. Will represent boundary conditions on computational boundary
                                       !! (or permanent boundary between fast-moving and near-stagnant ice
                                       !! FOR NOW: 1=interior bdry, 0=no-flow boundary, 2=stress bdry condition,
                                       !! 3=inhomogeneous Dirichlet boundary, 4=flux boundary: at these faces a flux
                                       !! will be specified which will override velocities; a homogeneous velocity
                                       !! condition will be specified (this seems to give the solver less difficulty)
  real, pointer, dimension(:,:) :: v_face_mask => NULL()  !< A mask for velocity boundary conditions on the C-grid
                                       !! v-face, with valued defined similarly to u_face_mask.
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
                                                     !! on corner-points (B grid) [degC]
  real, pointer, dimension(:,:) :: tmask => NULL()   !< A mask on tracer points that is 1 where there is ice.
  real, pointer, dimension(:,:) :: ice_visc => NULL()   !< Glen's law ice viscosity, often in [R L4 Z T-1 ~> kg m2 s-1].
  real, pointer, dimension(:,:) :: thickness_bdry_val => NULL() !< The ice thickness at an inflowing boundary [Z ~> m].
  real, pointer, dimension(:,:) :: u_bdry_val => NULL() !< The zonal ice velocity at inflowing boundaries
                                       !! [L yr-1 ~> m yr-1]
  real, pointer, dimension(:,:) :: v_bdry_val => NULL() !< The meridional ice velocity at inflowing boundaries
                                       !! [L yr-1 ~> m yr-1]
  real, pointer, dimension(:,:) :: h_bdry_val => NULL() !< The ice thickness at inflowing boundaries [m].
  real, pointer, dimension(:,:) :: t_bdry_val => NULL() !< The ice temperature at inflowing boundaries [degC].

  real, pointer, dimension(:,:) :: basal_traction => NULL() !< The area integrated nonlinear part of "linearized"
                                                            !! basal stress [R Z L2 T-1 ~> kg s-1].
                !!  The exact form depends on basal law exponent and/or whether flow is "hybridized" a la Goldberg 2011

  real, pointer, dimension(:,:) :: OD_rt => NULL()         !< A running total for calculating OD_av.
  real, pointer, dimension(:,:) :: ground_frac_rt => NULL() !< A running total for calculating ground_frac.
  real, pointer, dimension(:,:) :: OD_av => NULL()         !< The time average open ocean depth [Z ~> m].
  real, pointer, dimension(:,:) :: ground_frac => NULL()   !< Fraction of the time a cell is "exposed", i.e. the column
                               !! thickness is below a threshold and interacting with the rock [nondim].  When this
                               !! is 1, the ice-shelf is grounded
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
                            !! i.e. dt <= CFL_factor * min(dx / u)

  real :: A_glen_isothermal !< Ice viscosity parameter in Glen's Law, [Pa-3 s-1].
  real :: n_glen            !< Nonlinearity exponent in Glen's Law
  real :: eps_glen_min      !< Min. strain rate to avoid infinite Glen's law viscosity, [year-1].
  real :: C_basal_friction  !< Coefficient in sliding law tau_b = C u^(n_basal_fric), in
                            !!  units= Pa (m yr-1)-(n_basal_fric)
  real :: n_basal_fric      !< Exponent in sliding law tau_b = C u^(m_slide)
  real :: density_ocean_avg !< A typical ocean density [R ~> kg m-3].  This does not affect ocean
                            !! circulation or thermodynamics.  It is used to estimate the
                            !! gravitational driving force at the shelf front (until we think of
                            !! a better way to do it, but any difference will be negligible).
  real :: thresh_float_col_depth !< The water column depth over which the shelf if considered to be floating
  logical :: moving_shelf_front  !< Specify whether to advance shelf front (and calve).
  logical :: calve_to_mask       !< If true, calve off the ice shelf when it passes the edge of a mask.
  real :: min_thickness_simple_calve !< min. ice shelf thickness criteria for calving [Z ~> m].

  real :: cg_tolerance !< The tolerance in the CG solver, relative to initial residual, that
                       !! determines when to stop the conjugate gradient iterations.
  real :: nonlinear_tolerance !< The fractional nonlinear tolerance, relative to the initial error,
                              !! that sets when to stop the iterative velocity solver
  integer :: cg_max_iterations !< The maximum number of iterations that can be used in the CG solver
  integer :: nonlin_solve_err_mode  !< 1: exit vel solve based on nonlin residual
                    !! 2: exit based on "fixed point" metric (|u - u_last| / |u| < tol) where | | is infty-norm

  ! ids for outputting intermediate thickness in advection subroutine (debugging)
  !integer :: id_h_after_uflux = -1, id_h_after_vflux = -1, id_h_after_adv = -1

  logical :: debug                !< If true, write verbose checksums for debugging purposes
                                  !! and use reproducible sums
  logical :: module_is_initialized = .false. !< True if this module has been initialized.

  !>@{ Diagnostic handles
  integer :: id_u_shelf = -1, id_v_shelf = -1, id_t_shelf = -1, &
             id_ground_frac = -1, id_col_thick = -1, id_OD_av = -1, &
             id_u_mask = -1, id_v_mask = -1, id_t_mask = -1
  !>@}
  ! ids for outputting intermediate thickness in advection subroutine (debugging)
  !integer :: id_h_after_uflux = -1, id_h_after_vflux = -1, id_h_after_adv = -1

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

  p2 = (X(4)-X(1))**2 + (Y(4)-Y(1))**2 ; q2 = (X(3)-X(2))**2 + (Y(3)-Y(2))**2
  a2 = (X(3)-X(4))**2 + (Y(3)-Y(4))**2 ; c2 = (X(1)-X(2))**2 + (Y(1)-Y(2))**2
  b2 = (X(2)-X(4))**2 + (Y(2)-Y(4))**2 ; d2 = (X(3)-X(1))**2 + (Y(3)-Y(1))**2
  quad_area = .25 * sqrt(4*P2*Q2-(B2+D2-A2-C2)**2)

end function quad_area

!> This subroutine is used to register any fields related to the ice shelf
!! dynamics that should be written to or read from the restart file.
subroutine register_ice_shelf_dyn_restarts(G, param_file, CS, restart_CS)
  type(ocean_grid_type),  intent(inout) :: G    !< The grid type describing the ice shelf grid.
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(ice_shelf_dyn_CS), pointer       :: CS !< A pointer to the ice shelf dynamics control structure
  type(MOM_restart_CS),   pointer       :: restart_CS !< A pointer to the restart control structure.

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
    allocate( CS%u_shelf(IsdB:IedB,JsdB:JedB) ) ; CS%u_shelf(:,:) = 0.0
    allocate( CS%v_shelf(IsdB:IedB,JsdB:JedB) ) ; CS%v_shelf(:,:) = 0.0
    allocate( CS%t_shelf(isd:ied,jsd:jed) )   ; CS%t_shelf(:,:) = -10.0
    allocate( CS%ice_visc(isd:ied,jsd:jed) )    ; CS%ice_visc(:,:) = 0.0
    allocate( CS%basal_traction(isd:ied,jsd:jed) ) ; CS%basal_traction(:,:) = 0.0
    allocate( CS%OD_av(isd:ied,jsd:jed) )       ; CS%OD_av(:,:) = 0.0
    allocate( CS%ground_frac(isd:ied,jsd:jed) )  ; CS%ground_frac(:,:) = 0.0

    ! additional restarts for ice shelf state
    call register_restart_field(CS%u_shelf, "u_shelf", .false., restart_CS, &
                                "ice sheet/shelf u-velocity", "m s-1", hor_grid='Bu')
    call register_restart_field(CS%v_shelf, "v_shelf", .false., restart_CS, &
                                "ice sheet/shelf v-velocity", "m s-1", hor_grid='Bu')
    call register_restart_field(CS%t_shelf, "t_shelf", .true., restart_CS, &
                                "ice sheet/shelf vertically averaged temperature", "deg C")
    call register_restart_field(CS%OD_av, "OD_av", .true., restart_CS, &
                                "Average open ocean depth in a cell","m")
    call register_restart_field(CS%ground_frac, "ground_frac", .true., restart_CS, &
                                "fractional degree of grounding", "nondim")
    call register_restart_field(CS%ice_visc, "viscosity", .true., restart_CS, &
                                "Volume integrated Glens law ice viscosity", "kg m2 s-1")
    call register_restart_field(CS%basal_traction, "tau_b_beta", .true., restart_CS, &
                                "The area integrated basal traction coefficient", "kg s-1")
  endif

end subroutine register_ice_shelf_dyn_restarts

!> Initializes shelf model data, parameters and diagnostics
subroutine initialize_ice_shelf_dyn(param_file, Time, ISS, CS, G, US, diag, new_sim, solo_ice_sheet_in)
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(ocean_grid_type),   pointer       :: ocn_grid   !< The calling ocean model's horizontal grid structure
  type(time_type),         intent(inout) :: Time !< The clock that that will indicate the model time
  type(ice_shelf_state),   intent(in)    :: ISS  !< A structure with elements that describe
                                                 !! the ice-shelf state
  type(ice_shelf_dyn_CS),  pointer       :: CS   !< A pointer to the ice shelf dynamics control structure
  type(ocean_grid_type),   intent(inout) :: G    !< The grid type describing the ice shelf grid.
  type(unit_scale_type),   intent(in)    :: US   !< A structure containing unit conversion factors
  type(diag_ctrl), target, intent(in)    :: diag !< A structure that is used to regulate the diagnostic output.
  logical,                 intent(in)    :: new_sim !< If true this is a new simulation, otherwise
                                                 !! has been started from a restart file.
  logical,       optional, intent(in)    :: solo_ice_sheet_in !< If present, this indicates whether
                                                 !! a solo ice-sheet driver.

  ! Local variables
  real    :: Z_rescale  ! A rescaling factor for heights from the representation in
                        ! a restart file to the internal representation in this run.
  real    :: vel_rescale ! A rescaling factor for horizontal velocities from the representation
                        ! in a restart file to the internal representation in this run.
  !This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=200) :: config
  character(len=200) :: IC_file,filename,inputdir
  character(len=40)  :: var_name
  character(len=40)  :: mdl = "MOM_ice_shelf_dyn"  ! This module's name.
  logical :: shelf_mass_is_dynamic, override_shelf_movement, active_shelf_dynamics
  logical :: debug
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq, iters

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
    if (CS%GL_regularize .and. (CS%n_sub_regularize == 0)) call MOM_error (FATAL, &
      "GROUNDING_LINE_INTERP_SUBGRID_N must be a positive integer if GL regularization is used")
    call get_param(param_file, mdl, "ICE_SHELF_CFL_FACTOR", CS%CFL_factor, &
                 "A factor used to limit timestep as CFL_FACTOR * min (\Delta x / u). "//&
                 "This is only used with an ice-only model.", default=0.25)
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
                 units="m s-2", default = 9.80, scale=US%m_s_to_L_T**2*US%Z_to_m)

    call get_param(param_file, mdl, "A_GLEN_ISOTHERM", CS%A_glen_isothermal, &
                 "Ice viscosity parameter in Glen's Law", &
                 units="Pa-3 yr-1", default=9.461e-18, scale=1.0/(365.0*86400.0))
                 ! This default is equivalent to 3.0001e-25 Pa-3 s-1, appropriate at about -10 C.
    call get_param(param_file, mdl, "GLEN_EXPONENT", CS%n_glen, &
                 "nonlinearity exponent in Glen's Law", &
                  units="none", default=3.)
    call get_param(param_file, mdl, "MIN_STRAIN_RATE_GLEN", CS%eps_glen_min, &
                 "min. strain rate to avoid infinite Glen's law viscosity", &
                 units="a-1", default=1.e-12, scale=US%T_to_s/(365.0*86400.0))
    call get_param(param_file, mdl, "BASAL_FRICTION_EXP", CS%n_basal_fric, &
                 "Exponent in sliding law \tau_b = C u^(n_basal_fric)", &
                 units="none", fail_if_missing=.true.)
    call get_param(param_file, mdl, "BASAL_FRICTION_COEFF", CS%C_basal_friction, &
                 "Coefficient in sliding law \tau_b = C u^(n_basal_fric)", &
                 units="Pa (m yr-1)-(n_basal_fric)", scale=US%kg_m2s_to_RZ_T*((365.0*86400.0)**CS%n_basal_fric), &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, "DENSITY_ICE", CS%density_ice, &
                 "A typical density of ice.", units="kg m-3", default=917.0, scale=US%kg_m3_to_R)
    call get_param(param_file, mdl, "CONJUGATE_GRADIENT_TOLERANCE", CS%cg_tolerance, &
                "tolerance in CG solver, relative to initial residual", default=1.e-6)
    call get_param(param_file, mdl, "ICE_NONLINEAR_TOLERANCE", CS%nonlinear_tolerance, &
                "nonlin tolerance in iterative velocity solve",default=1.e-6)
    call get_param(param_file, mdl, "CONJUGATE_GRADIENT_MAXIT", CS%cg_max_iterations, &
                "max iteratiions in CG solver", default=2000)
    call get_param(param_file, mdl, "THRESH_FLOAT_COL_DEPTH", CS%thresh_float_col_depth, &
                "min ocean thickness to consider ice *floating*; "//&
                "will only be important with use of tides", &
                units="m", default=1.e-3, scale=US%m_to_Z)
    call get_param(param_file, mdl, "NONLIN_SOLVE_ERR_MODE", CS%nonlin_solve_err_mode, &
                "Choose whether nonlin error in vel solve is based on nonlinear "//&
                "residual (1) or relative change since last iteration (2)", default=1)

    call get_param(param_file, mdl, "SHELF_MOVING_FRONT", CS%moving_shelf_front, &
                 "Specify whether to advance shelf front (and calve).", &
                 default=.true.)
    call get_param(param_file, mdl, "CALVE_TO_MASK", CS%calve_to_mask, &
                 "If true, do not allow an ice shelf where prohibited by a mask.", &
                 default=.false.)
  endif
  call get_param(param_file, mdl, "MIN_THICKNESS_SIMPLE_CALVE", CS%min_thickness_simple_calve, &
                 "Min thickness rule for the VERY simple calving law",&
                 units="m", default=0.0, scale=US%m_to_Z)

  ! Allocate memory in the ice shelf dynamics control structure that was not
  ! previously allocated for registration for restarts.
  ! OVS vertically integrated Temperature

  if (active_shelf_dynamics) then
    ! DNG
    allocate( CS%u_bdry_val(Isdq:Iedq,Jsdq:Jedq) ) ; CS%u_bdry_val(:,:) = 0.0
    allocate( CS%v_bdry_val(Isdq:Iedq,Jsdq:Jedq) ) ; CS%v_bdry_val(:,:) = 0.0
    allocate( CS%t_bdry_val(isd:ied,jsd:jed) )   ; CS%t_bdry_val(:,:) = -15.0
    allocate( CS%h_bdry_val(isd:ied,jsd:jed) ) ; CS%h_bdry_val(:,:) = 0.0
    allocate( CS%thickness_bdry_val(isd:ied,jsd:jed) ) ; CS%thickness_bdry_val(:,:) = 0.0
    allocate( CS%u_face_mask(Isdq:Iedq,jsd:jed) ) ; CS%u_face_mask(:,:) = 0.0
    allocate( CS%v_face_mask(isd:ied,Jsdq:Jedq) ) ; CS%v_face_mask(:,:) = 0.0
    allocate( CS%u_face_mask_bdry(Isdq:Iedq,jsd:jed) ) ; CS%u_face_mask_bdry(:,:) = -2.0
    allocate( CS%v_face_mask_bdry(isd:ied,Jsdq:Jedq) ) ; CS%v_face_mask_bdry(:,:) = -2.0
    allocate( CS%u_flux_bdry_val(Isdq:Iedq,jsd:jed) ) ; CS%u_flux_bdry_val(:,:) = 0.0
    allocate( CS%v_flux_bdry_val(isd:ied,Jsdq:Jedq) ) ; CS%v_flux_bdry_val(:,:) = 0.0
    allocate( CS%umask(Isdq:Iedq,Jsdq:Jedq) ) ; CS%umask(:,:) = -1.0
    allocate( CS%vmask(Isdq:Iedq,Jsdq:Jedq) ) ; CS%vmask(:,:) = -1.0
    allocate( CS%tmask(Isdq:Iedq,Jsdq:Jedq) ) ; CS%tmask(:,:) = -1.0

    CS%OD_rt_counter = 0
    allocate( CS%OD_rt(isd:ied,jsd:jed) ) ; CS%OD_rt(:,:) = 0.0
    allocate( CS%ground_frac_rt(isd:ied,jsd:jed) ) ; CS%ground_frac_rt(:,:) = 0.0

    if (CS%calve_to_mask) then
      allocate( CS%calve_mask(isd:ied,jsd:jed) ) ; CS%calve_mask(:,:) = 0.0
    endif

    CS%elapsed_velocity_time = 0.0

    call update_velocity_masks(CS, G, ISS%hmask, CS%umask, CS%vmask, CS%u_face_mask, CS%v_face_mask)
  endif

  ! Take additional initialization steps, for example of dependent variables.
  if (active_shelf_dynamics .and. .not.new_sim) then
    if ((US%m_to_Z_restart /= 0.0) .and. (US%m_to_Z_restart /= US%m_to_Z)) then
      Z_rescale = US%m_to_Z / US%m_to_Z_restart
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        CS%OD_av(i,j) = Z_rescale * CS%OD_av(i,j)
      enddo ; enddo
    endif

    if ((US%m_to_L_restart*US%s_to_T_restart /= 0.0) .and. &
        (US%m_to_L_restart /= US%m_s_to_L_T*US%s_to_T_restart)) then
      vel_rescale = US%m_s_to_L_T*US%s_to_T_restart / US%m_to_L_restart
      do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
        CS%u_shelf(I,J) = vel_rescale * CS%u_shelf(I,J)
        CS%v_shelf(I,J) = vel_rescale * CS%v_shelf(I,J)
      enddo ; enddo
    endif

    ! this is unfortunately necessary; if grid is not symmetric the boundary values
    !  of u and v are otherwise not set till the end of the first linear solve, and so
    !  viscosity is not calculated correctly.
    ! This has to occur after init_boundary_values or some of the arrays on the
    ! right hand side have not been set up yet.
    if (.not. G%symmetric) then
      do j=G%jsd,G%jed ; do i=G%isd,G%ied
        if (((i+G%idg_offset) == (G%domain%nihalo+1)).and.(CS%u_face_mask(I-1,j) == 3)) then
          CS%u_shelf(I-1,J-1) = CS%u_bdry_val(I-1,J-1)
          CS%u_shelf(I-1,J) = CS%u_bdry_val(I-1,J)
          CS%v_shelf(I-1,J-1) = CS%v_bdry_val(I-1,J-1)
          CS%v_shelf(I-1,J) = CS%v_bdry_val(I-1,J)
        endif
        if (((j+G%jdg_offset) == (G%domain%njhalo+1)).and.(CS%v_face_mask(i,J-1) == 3)) then
          CS%u_shelf(I-1,J-1) = CS%u_bdry_val(I-1,J-1)
          CS%u_shelf(I,J-1) = CS%u_bdry_val(I,J-1)
          CS%v_shelf(I-1,J-1) = CS%v_bdry_val(I-1,J-1)
          CS%v_shelf(I,J-1) = CS%v_bdry_val(I,J-1)
        endif
      enddo ; enddo
    endif

    call pass_var(CS%OD_av,G%domain)
    call pass_var(CS%ground_frac,G%domain)
    call pass_var(CS%ice_visc,G%domain)
    call pass_var(CS%basal_traction, G%domain)
    call pass_vector(CS%u_shelf, CS%v_shelf, G%domain, TO_ALL, BGRID_NE)
  endif

  if (active_shelf_dynamics) then
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

!    call init_boundary_values(CS, G, time, ISS%hmask, CS%input_flux, CS%input_thickness, new_sim)

    if (new_sim) then
      call MOM_mesg("MOM_ice_shelf.F90, initialize_ice_shelf: initialize ice velocity.")
      call update_OD_ffrac_uncoupled(CS, G, ISS%h_shelf(:,:))
      call ice_shelf_solve_outer(CS, ISS, G, US, CS%u_shelf, CS%v_shelf, iters, Time)

      if (CS%id_u_shelf > 0) call post_data(CS%id_u_shelf, CS%u_shelf, CS%diag)
      if (CS%id_v_shelf > 0) call post_data(CS%id_v_shelf, CS%v_shelf,CS%diag)
    endif

  ! Register diagnostics.
    CS%id_u_shelf = register_diag_field('ocean_model','u_shelf',CS%diag%axesCu1, Time, &
       'x-velocity of ice', 'm yr-1', conversion=365.0*86400.0*US%L_T_to_m_s)
    CS%id_v_shelf = register_diag_field('ocean_model','v_shelf',CS%diag%axesCv1, Time, &
       'y-velocity of ice', 'm yr-1', conversion=365.0*86400.0*US%L_T_to_m_s)
    CS%id_u_mask = register_diag_field('ocean_model','u_mask',CS%diag%axesCu1, Time, &
       'mask for u-nodes', 'none')
    CS%id_v_mask = register_diag_field('ocean_model','v_mask',CS%diag%axesCv1, Time, &
       'mask for v-nodes', 'none')
!    CS%id_surf_elev = register_diag_field('ocean_model','ice_surf',CS%diag%axesT1, Time, &
!       'ice surf elev', 'm')
    CS%id_ground_frac = register_diag_field('ocean_model','ice_ground_frac',CS%diag%axesT1, Time, &
       'fraction of cell that is grounded', 'none')
    CS%id_col_thick = register_diag_field('ocean_model','col_thick',CS%diag%axesT1, Time, &
       'ocean column thickness passed to ice model', 'm', conversion=US%Z_to_m)
    CS%id_OD_av = register_diag_field('ocean_model','OD_av',CS%diag%axesT1, Time, &
       'intermediate ocean column thickness passed to ice model', 'm', conversion=US%Z_to_m)
    !CS%id_h_after_uflux = register_diag_field('ocean_model','h_after_uflux',CS%diag%axesh1, Time, &
    !   'thickness after u flux ', 'none')
    !CS%id_h_after_vflux = register_diag_field('ocean_model','h_after_vflux',CS%diag%axesh1, Time, &
    !   'thickness after v flux ', 'none')
    !CS%id_h_after_adv = register_diag_field('ocean_model','h_after_adv',CS%diag%axesh1, Time, &
    !   'thickness after front adv ', 'none')

!!! OVS vertically integrated temperature
    CS%id_t_shelf = register_diag_field('ocean_model','t_shelf',CS%diag%axesT1, Time, &
       'T of ice', 'oC')
    CS%id_t_mask = register_diag_field('ocean_model','tmask',CS%diag%axesT1, Time, &
       'mask for T-nodes', 'none')
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

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  dummy_time = set_time(0,0)
  isd=G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j=jsd,jed
    do i=isd,ied
      OD = G%bathyT(i,j) - rhoi_rhow * ISS%h_shelf(i,j)
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

  call ice_shelf_solve_outer(CS, ISS, G, US, CS%u_shelf, CS%v_shelf, iters, dummy_time)

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
  do j=G%jsc,G%jec ; do i=G%isc,G%iec ; if (ISS%hmask(i,j) == 1.0) then
    dt_local = 2.0*G%areaT(i,j) / &
       ((G%dyCu(I,j)  * max(abs(CS%u_shelf(I,J)  + CS%u_shelf(I,j-1)), min_vel) + &
         G%dyCu(I-1,j)* max(abs(CS%u_shelf(I-1,J)+ CS%u_shelf(I-1,j-1)), min_vel)) + &
        (G%dxCv(i,J)  * max(abs(CS%v_shelf(i,J)  + CS%v_shelf(i-1,J)), min_vel) + &
         G%dxCv(i,J-1)* max(abs(CS%v_shelf(i,J-1)+ CS%v_shelf(i-1,J-1)), min_vel)))

    min_dt = min(min_dt, dt_local)
  endif ; enddo ; enddo ! i- and j- loops

  call min_across_PEs(min_dt)

  ice_time_step_CFL = CS%CFL_factor * min_dt

end function ice_time_step_CFL

!> This subroutine updates the ice shelf velocities, mass, stresses and properties due to the
!! ice shelf dynamics.
subroutine update_ice_shelf(CS, ISS, G, US, time_step, Time, ocean_mass, coupled_grounding, must_update_vel)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                              !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real,                   intent(in)    :: time_step !< time step [T ~> s]
  type(time_type),        intent(in)    :: Time !< The current model time
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

  call ice_shelf_advect(CS, ISS, G, time_step, Time)
  CS%elapsed_velocity_time = CS%elapsed_velocity_time + time_step
  if (CS%elapsed_velocity_time >= CS%velocity_update_time_step) update_ice_vel = .true.

  if (coupled_GL) then
    call update_OD_ffrac(CS, G, US, ocean_mass, update_ice_vel)
  elseif (update_ice_vel) then
    call update_OD_ffrac_uncoupled(CS, G, ISS%h_shelf(:,:))
  endif

  if (update_ice_vel) then
    call ice_shelf_solve_outer(CS, ISS, G, US, CS%u_shelf, CS%v_shelf, iters, Time)
  endif

  call ice_shelf_temp(CS, ISS, G, US, time_step, ISS%water_flux, Time)

  if (update_ice_vel) then
    call enable_averages(CS%elapsed_velocity_time, Time, CS%diag)
    if (CS%id_col_thick > 0) call post_data(CS%id_col_thick, CS%OD_av, CS%diag)
    if (CS%id_u_shelf > 0) call post_data(CS%id_u_shelf, CS%u_shelf, CS%diag)
    if (CS%id_v_shelf > 0) call post_data(CS%id_v_shelf, CS%v_shelf, CS%diag)
    if (CS%id_t_shelf > 0) call post_data(CS%id_t_shelf,CS%t_shelf,CS%diag)
    if (CS%id_ground_frac > 0) call post_data(CS%id_ground_frac, CS%ground_frac,CS%diag)
    if (CS%id_OD_av >0) call post_data(CS%id_OD_av, CS%OD_av,CS%diag)

    if (CS%id_u_mask > 0) call post_data(CS%id_u_mask,CS%umask,CS%diag)
    if (CS%id_v_mask > 0) call post_data(CS%id_v_mask,CS%vmask,CS%diag)
    if (CS%id_t_mask > 0) call post_data(CS%id_t_mask,CS%tmask,CS%diag)

    call disable_averaging(CS%diag)

    CS%elapsed_velocity_time = 0.0
  endif

end subroutine update_ice_shelf

!> This subroutine takes the velocity (on the Bgrid) and timesteps h_t = - div (uh) once.
!! Additionally, it will update the volume of ice in partially-filled cells, and update
!! hmask accordingly
subroutine ice_shelf_advect(CS, ISS, G, time_step, Time)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  real,                   intent(in)    :: time_step !< time step [T ~> s]
  type(time_type),        intent(in)    :: Time !< The current model time


! 3/8/11 DNG
!
!    This subroutine takes the velocity (on the Bgrid) and timesteps h_t = - div (uh) once.
!    ADDITIONALLY, it will update the volume of ice in partially-filled cells, and update
!        hmask accordingly
!
!    The flux overflows are included here. That is because they will be used to advect 3D scalars
!    into partial cells

  real, dimension(SZDI_(G),SZDJ_(G))   :: h_after_uflux, h_after_vflux ! Ice thicknesses [Z ~> m].
  real, dimension(SZDIB_(G),SZDJ_(G))  :: uh_ice  ! The accumulated zonal ice volume flux [Z L2 ~> m3]
  real, dimension(SZDI_(G),SZDJB_(G))  :: vh_ice  ! The accumulated meridional ice volume flux [Z L2 ~> m3]
  type(loop_bounds_type) :: LB
  integer                           :: isd, ied, jsd, jed, i, j, isc, iec, jsc, jec, stencil

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  uh_ice(:,:) = 0.0
  vh_ice(:,:) = 0.0

  h_after_uflux(:,:) = 0.0
  h_after_vflux(:,:) = 0.0
  ! call MOM_mesg("MOM_ice_shelf.F90: ice_shelf_advect called")

  do j=jsd,jed ; do i=isd,ied ; if (CS%thickness_bdry_val(i,j) /= 0.0) then
    ISS%h_shelf(i,j) = CS%thickness_bdry_val(i,j)
  endif ; enddo ; enddo

  stencil = 2
  LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc-stencil ; LB%jeh = G%jec+stencil
  if (LB%jsh < jsd) call MOM_error(FATAL, &
    "ice_shelf_advect:  Halo is too small for the ice thickness advection stencil.")

  call ice_shelf_advect_thickness_x(CS, G, LB, time_step, ISS%hmask, ISS%h_shelf, h_after_uflux, uh_ice)

!  call enable_averages(time_step, Time, CS%diag)
!  call pass_var(h_after_uflux, G%domain)
!  if (CS%id_h_after_uflux > 0) call post_data(CS%id_h_after_uflux, h_after_uflux, CS%diag)
!  call disable_averaging(CS%diag)

  LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
  call ice_shelf_advect_thickness_y(CS, G, LB, time_step, ISS%hmask, h_after_uflux, h_after_vflux, vh_ice)

!  call enable_averages(time_step, Time, CS%diag)
!  call pass_var(h_after_vflux, G%domain)
!  if (CS%id_h_after_vflux > 0) call post_data(CS%id_h_after_vflux, h_after_vflux, CS%diag)
!  call disable_averaging(CS%diag)

  do j=jsd,jed
    do i=isd,ied
      if (ISS%hmask(i,j) == 1) ISS%h_shelf(i,j) = h_after_vflux(i,j)
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
  endif

  !call enable_averages(time_step, Time, CS%diag)
  !if (CS%id_h_after_adv > 0) call post_data(CS%id_h_after_adv, ISS%h_shelf, CS%diag)
  !call disable_averaging(CS%diag)

  !call change_thickness_using_melt(ISS, G, time_step, fluxes, CS%density_ice)

  call update_velocity_masks(CS, G, ISS%hmask, CS%umask, CS%vmask, CS%u_face_mask, CS%v_face_mask)

end subroutine ice_shelf_advect

subroutine ice_shelf_solve_outer(CS, ISS, G, US, u_shlf, v_shlf, iters, time)
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

  real, dimension(SZDIB_(G),SZDJB_(G)) :: taudx, taudy ! Driving stresses at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: u_bdry_cont ! Boundary u-stress contribution [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: v_bdry_cont ! Boundary v-stress contribution [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: Au, Av ! The retarding lateral stress contributions [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: err_u, err_v
  real, dimension(SZDIB_(G),SZDJB_(G)) :: u_last, v_last ! Previous velocities [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: H_node ! Ice shelf thickness at corners [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)) :: float_cond ! An array indicating where the ice
                                                ! shelf is floating: 0 if floating, 1 if not.
  character(len=160) :: mesg  ! The text of an error message
  integer :: conv_flag, i, j, k,l, iter
  integer :: isdq, iedq, jsdq, jedq, isd, ied, jsd, jed, nodefloat, nsub
  real    :: err_max, err_tempu, err_tempv, err_init, area, max_vel, tempu, tempv
  real    :: rhoi_rhow ! The density of ice divided by a typical water density [nondim]
  real, pointer, dimension(:,:,:,:) :: Phi => NULL() ! The gradients of bilinear basis elements at Gaussian
                                                ! quadrature points surrounding the cell vertices [m-1].
  real, pointer, dimension(:,:,:,:,:,:) :: Phisub => NULL() ! Quadrature structure weights at subgridscale
                                                !  locations for finite element calculations [nondim]
  character(2)                :: iternum
  character(2)                :: numproc

  ! for GL interpolation
  nsub = CS%n_sub_regularize

  isdq = G%isdB ; iedq = G%iedB ; jsdq = G%jsdB ; jedq = G%jedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  rhoi_rhow = CS%density_ice / CS%density_ocean_avg

  taudx(:,:) = 0.0 ; taudy(:,:) = 0.0
  u_bdry_cont(:,:) = 0.0 ; v_bdry_cont(:,:) = 0.0
  Au(:,:) = 0.0 ; Av(:,:) = 0.0

  ! need to make these conditional on GL interpolation
  float_cond(:,:) = 0.0 ; H_node(:,:) = 0.0
  allocate(Phisub(nsub,nsub,2,2,2,2)) ; Phisub(:,:,:,:,:,:) = 0.0

  call calc_shelf_driving_stress(CS, ISS, G, US, taudx, taudy, CS%OD_av)

  ! This is to determine which cells contain the grounding line, the criterion being that the cell
  ! is ice-covered, with some nodes floating and some grounded flotation condition is estimated by
  ! assuming topography is cellwise constant and H is bilinear in a cell; floating where
  ! rho_i/rho_w * H_node - D is negative

  ! need to make this conditional on GL interp

  if (CS%GL_regularize) then

    call interpolate_H_to_B(G, ISS%h_shelf, ISS%hmask, H_node)

    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      nodefloat = 0

      do l=0,1 ; do k=0,1
        if ((ISS%hmask(i,j) == 1) .and. &
            (rhoi_rhow * H_node(i-1+k,j-1+l) - G%bathyT(i,j) <= 0)) then
          nodefloat = nodefloat + 1
        endif
      enddo ; enddo
      if ((nodefloat > 0) .and. (nodefloat < 4)) then
        float_cond(i,j) = 1.0
        CS%ground_frac(i,j) = 1.0
      endif
    enddo ; enddo

    call pass_var(float_cond, G%Domain)

    call bilinear_shape_functions_subgrid(Phisub, nsub)

  endif

  ! must prepare Phi
  allocate(Phi(1:8,1:4,isd:ied,jsd:jed)) ; Phi(:,:,:,:) = 0.0

  do j=jsd,jed ; do i=isd,ied
    call bilinear_shape_fn_grid(G, i, j, Phi(:,:,i,j))
  enddo ; enddo

  call calc_shelf_visc(CS, ISS, G, US, u_shlf, v_shlf)

  call pass_var(CS%ice_visc, G%domain)
  call pass_var(CS%basal_traction, G%domain)

  ! This makes sure basal stress is only applied when it is supposed to be
  do j=G%jsd,G%jed ; do i=G%isd,G%ied
    CS%basal_traction(i,j) = CS%basal_traction(i,j) * CS%ground_frac(i,j)
  enddo ; enddo

  call apply_boundary_values(CS, ISS, G, US, time, Phisub, H_node, CS%ice_visc, &
                             CS%basal_traction, float_cond, rhoi_rhow, u_bdry_cont, v_bdry_cont)

  Au(:,:) = 0.0 ; Av(:,:) = 0.0

  call CG_action(Au, Av, u_shlf, v_shlf, Phi, Phisub, CS%umask, CS%vmask, ISS%hmask, H_node, &
                 CS%ice_visc, float_cond, G%bathyT, CS%basal_traction, &
                 G, US, G%isc-1, G%iec+1, G%jsc-1, G%jec+1, rhoi_rhow)

  if (CS%nonlin_solve_err_mode == 1) then
    err_init = 0 ; err_tempu = 0 ; err_tempv = 0
    do J=G%IscB,G%JecB ; do I=G%IscB,G%IecB
      if (CS%umask(I,J) == 1) then
        err_tempu = ABS(Au(I,J) + u_bdry_cont(I,J) - taudx(I,J))
        if (err_tempu >= err_init) err_init = err_tempu
      endif
      if (CS%vmask(I,J) == 1) then
        err_tempv = ABS(Av(I,J) + v_bdry_cont(I,J) - taudy(I,J))
        if (err_tempv >= err_init) err_init = err_tempv
      endif
    enddo ; enddo

    call max_across_PEs(err_init)
  endif

  u_last(:,:) = u_shlf(:,:) ; v_last(:,:) = v_shlf(:,:)

  !! begin loop

  do iter=1,100

    call ice_shelf_solve_inner(CS, ISS, G, US, u_shlf, v_shlf, taudx, taudy, H_node, float_cond, &
                               ISS%hmask, conv_flag, iters, time, Phi, Phisub)

    if (CS%debug) then
      call qchksum(u_shlf, "u shelf", G%HI, haloshift=2, scale=US%L_T_to_m_s)
      call qchksum(v_shlf, "v shelf", G%HI, haloshift=2, scale=US%L_T_to_m_s)
    endif

    write(mesg,*) "ice_shelf_solve_outer: linear solve done in ",iters," iterations"
    call MOM_mesg(mesg, 5)

    call calc_shelf_visc(CS, ISS, G, US, u_shlf, v_shlf)
    call pass_var(CS%ice_visc, G%domain)
    call pass_var(CS%basal_traction, G%domain)

    ! makes sure basal stress is only applied when it is supposed to be
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      CS%basal_traction(i,j) = CS%basal_traction(i,j) * CS%ground_frac(i,j)
    enddo ; enddo

    u_bdry_cont(:,:) = 0 ; v_bdry_cont(:,:) = 0

    call apply_boundary_values(CS, ISS, G, US, time, Phisub, H_node, CS%ice_visc, &
                               CS%basal_traction, float_cond, rhoi_rhow, u_bdry_cont, v_bdry_cont)

    Au(:,:) = 0 ; Av(:,:) = 0

    call CG_action(Au, Av, u_shlf, v_shlf, Phi, Phisub, CS%umask, CS%vmask, ISS%hmask, H_node, &
                   CS%ice_visc, float_cond, G%bathyT, CS%basal_traction, &
                   G, US, G%isc-1, G%iec+1, G%jsc-1, G%jec+1, rhoi_rhow)

    err_max = 0

    if (CS%nonlin_solve_err_mode == 1) then

      do J=G%jscB,G%jecB ; do I=G%jscB,G%iecB
        if (CS%umask(I,J) == 1) then
          err_tempu = ABS(Au(I,J) + u_bdry_cont(I,J) - taudx(I,J))
          if (err_tempu >= err_max) err_max = err_tempu
        endif
        if (CS%vmask(I,J) == 1) then
          err_tempv = ABS(Av(I,J) + v_bdry_cont(I,J) - taudy(I,J))
          if (err_tempv >= err_max) err_max = err_tempv
        endif
      enddo ; enddo

      call max_across_PEs(err_max)

    elseif (CS%nonlin_solve_err_mode == 2) then

      max_vel = 0 ; tempu = 0 ; tempv = 0
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
          tempv = SQRT(v_shlf(I,J)**2 + tempu**2)
        endif
        if (tempv >= max_vel) max_vel = tempv
      enddo ; enddo

      u_last(:,:) = u_shlf(:,:)
      v_last(:,:) = v_shlf(:,:)

      call max_across_PEs(max_vel)
      call max_across_PEs(err_max)
      err_init = max_vel
    endif

    write(mesg,*) "ice_shelf_solve_outer: nonlinear fractional residual = ", err_max/err_init
    call MOM_mesg(mesg, 5)

    if (err_max <= CS%nonlinear_tolerance * err_init) then
      write(mesg,*) "ice_shelf_solve_outer: exiting nonlinear solve after ",iter," iterations"
      call MOM_mesg(mesg, 5)
      exit
    endif

  enddo

  deallocate(Phi)
  deallocate(Phisub)

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
                          intent(in)    :: float_cond !< An array indicating where the ice
                                                !! shelf is floating: 0 if floating, 1 if not.
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
!    boundary contributions are added to taud to get the RHS
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
                        ubd, vbd, &   ! Boundary stress contributions [R L3 Z T-2 ~> kg m s-2]
                        Au, Av, & ! The retarding lateral stress contributions [R L3 Z T-2 ~> kg m s-2]
                        Du, Dv, & ! Velocity changes [L T-1 ~> m s-1]
                        sum_vec, sum_vec_2
  real    :: tol, beta_k, area, dot_p1, resid0, cg_halo
  real    :: num, denom
  real    :: alpha_k     ! A scaling factor for iterative corrections [nondim]
  real    :: resid_scale ! A scaling factor for redimensionalizing the global residuals [m2 L-2 ~> 1]
                         ! [m2 L-2 ~> 1] [R L3 Z T-2 ~> m kg s-2]
  real    :: resid2_scale ! A scaling factor for redimensionalizing the global squared residuals
                         ! [m2 L-2 ~> 1] [R L3 Z T-2 ~> m kg s-2]
  real    :: rhoi_rhow  ! The density of ice divided by a typical water density [nondim]
  integer :: iter, i, j, isd, ied, jsd, jed, isc, iec, jsc, jec, is, js, ie, je
  integer :: Is_sum, Js_sum, Ie_sum, Je_sum ! Loop bounds for global sums or arrays starting at 1.
  integer :: isdq, iedq, jsdq, jedq, iscq, iecq, jscq, jecq, nx_halo, ny_halo

  isdq = G%isdB ; iedq = G%iedB ; jsdq = G%jsdB ; jedq = G%jedB
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  ny_halo = G%domain%njhalo ; nx_halo = G%domain%nihalo
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg

  Zu(:,:) = 0 ; Zv(:,:) = 0 ; DIAGu(:,:) = 0 ; DIAGv(:,:) = 0
  Ru(:,:) = 0 ; Rv(:,:) = 0 ; Au(:,:) = 0 ; Av(:,:) = 0
  Du(:,:) = 0 ; Dv(:,:) = 0 ; ubd(:,:) = 0 ; vbd(:,:) = 0
  dot_p1 = 0

  ! Determine the loop limits for sums, bearing in mind that the arrays will be starting at 1.
  Is_sum = G%isc + (1-G%IsdB)
  Ie_sum = G%iecB + (1-G%IsdB)
  ! Include the edge if tile is at the western bdry;  Should add a test to avoid this if reentrant.
  if (G%isc+G%idg_offset==G%isg) Is_sum = G%IscB + (1-G%IsdB)

  Js_sum = G%jsc + (1-G%JsdB)
  Je_sum = G%jecB + (1-G%JsdB)
  ! Include the edge if tile is at the southern bdry;  Should add a test to avoid this if reentrant.
  if (G%jsc+G%jdg_offset==G%jsg) Js_sum = G%JscB + (1-G%JsdB)

  call apply_boundary_values(CS, ISS, G, US, time, Phisub, H_node, CS%ice_visc, &
                             CS%basal_traction, float_cond, rhoi_rhow, ubd, vbd)

  RHSu(:,:) = taudx(:,:) - ubd(:,:)
  RHSv(:,:) = taudy(:,:) - vbd(:,:)

  call pass_vector(RHSu, RHSv, G%domain, TO_ALL, BGRID_NE)

  call matrix_diagonal(CS, G, US, float_cond, H_node, CS%ice_visc, CS%basal_traction, &
                       hmask, rhoi_rhow, Phisub, DIAGu, DIAGv)

  call pass_vector(DIAGu, DIAGv, G%domain, TO_ALL, BGRID_NE)

  call CG_action(Au, Av, u_shlf, v_shlf, Phi, Phisub, CS%umask, CS%vmask, hmask, &
                 H_node, CS%ice_visc, float_cond, G%bathyT, CS%basal_traction, &
                 G, US, isc-1, iec+1, jsc-1, jec+1, rhoi_rhow)

  call pass_vector(Au, Av, G%domain, TO_ALL, BGRID_NE)

  Ru(:,:) = (RHSu(:,:) - Au(:,:))
  Rv(:,:) = (RHSv(:,:) - Av(:,:))

  resid_scale = US%L_to_m**2*US%s_to_T*US%RZ_to_kg_m2*US%L_T_to_m_s**2
  resid2_scale = (US%RZ_to_kg_m2*US%L_to_m*US%L_T_to_m_s**2)**2

  sum_vec(:,:) = 0.0
  do j=jscq,jecq ; do i=iscq,iecq
    if (CS%umask(I,J) == 1) sum_vec(I,J) = resid2_scale*Ru(I,J)**2
    if (CS%vmask(I,J) == 1) sum_vec(I,J) = sum_vec(I,J) + resid2_scale*Rv(I,J)**2
  enddo ; enddo

  dot_p1 = reproducing_sum( sum_vec, Js_sum, Ie_sum, Js_sum, Je_sum )

  resid0 = sqrt(dot_p1)

  do j=jsdq,jedq
    do i=isdq,iedq
      if (CS%umask(I,J) == 1) Zu(I,J) = Ru(I,J) / DIAGu(I,J)
      if (CS%vmask(I,J) == 1) Zv(I,J) = Rv(I,J) / DIAGv(I,J)
    enddo
  enddo

  Du(:,:) = Zu(:,:) ; Dv(:,:) = Zv(:,:)

  cg_halo = 3
  conv_flag = 0

  !!!!!!!!!!!!!!!!!!
  !!              !!
  !! MAIN CG LOOP !!
  !!              !!
  !!!!!!!!!!!!!!!!!!

  ! initially, c-grid data is valid up to 3 halo nodes out

  do iter = 1,CS%cg_max_iterations

    ! assume asymmetry
    ! thus we can never assume that any arrays are legit more than 3 vertices past
    ! the computational domain - this is their state in the initial iteration


    is = isc - cg_halo ; ie = iecq + cg_halo
    js = jscq - cg_halo ; je = jecq + cg_halo

    Au(:,:) = 0 ; Av(:,:) = 0

    call CG_action(Au, Av, Du, Dv, Phi, Phisub, CS%umask, CS%vmask, hmask, &
                   H_node, CS%ice_visc, float_cond, G%bathyT, CS%basal_traction, &
                   G, US, is, ie, js, je, rhoi_rhow)

    ! Au, Av valid region moves in by 1


    sum_vec(:,:) = 0.0 ; sum_vec_2(:,:) = 0.0

    do j=jscq,jecq ; do i=iscq,iecq
      if (CS%umask(I,J) == 1) then
        sum_vec(I,J) = resid_scale * Zu(I,J) * Ru(I,J)
        sum_vec_2(I,J) = resid_scale * Du(I,J) * Au(I,J)
      endif
      if (CS%vmask(I,J) == 1) then
        sum_vec(I,J) = sum_vec(I,J) + resid_scale * Zv(I,J) * Rv(I,J)
        sum_vec_2(I,J) = sum_vec_2(I,J) + resid_scale * Dv(I,J) * Av(I,J)
      endif
    enddo ; enddo

    alpha_k = reproducing_sum( sum_vec, Is_sum, Ie_sum, Js_sum, Je_sum ) / &
              reproducing_sum( sum_vec_2, Is_sum, Ie_sum, Js_sum, Je_sum )


    do j=jsd,jed ; do i=isd,ied
      if (CS%umask(I,J) == 1) u_shlf(I,J) = u_shlf(I,J) + alpha_k * Du(I,J)
      if (CS%vmask(I,J) == 1) v_shlf(I,J) = v_shlf(I,J) + alpha_k * Dv(I,J)
    enddo ; enddo

    do j=jsd,jed ; do i=isd,ied
      if (CS%umask(I,J) == 1) then
        Ru_old(I,J) = Ru(I,J) ; Zu_old(I,J) = Zu(I,J)
      endif
      if (CS%vmask(I,J) == 1) then
        Rv_old(I,J) = Rv(I,J) ; Zv_old(I,J) = Zv(I,J)
      endif
    enddo ; enddo

!    Ru(:,:) = Ru(:,:) - alpha_k * Au(:,:)
!    Rv(:,:) = Rv(:,:) - alpha_k * Av(:,:)

    do j=jsd,jed
      do i=isd,ied
        if (CS%umask(I,J) == 1) Ru(I,J) = Ru(I,J) - alpha_k * Au(I,J)
        if (CS%vmask(I,J) == 1) Rv(I,J) = Rv(I,J) - alpha_k * Av(I,J)
      enddo
    enddo

    do j=jsdq,jedq
      do i=isdq,iedq
        if (CS%umask(I,J) == 1) then
          Zu(I,J) = Ru(I,J) / DIAGu(I,J)
        endif
        if (CS%vmask(I,J) == 1) then
          Zv(I,J) = Rv(I,J) / DIAGv(I,J)
        endif
      enddo
    enddo

    ! R,u,v,Z valid region moves in by 1

    ! beta_k = (Z \dot R) / (Zold \dot Rold}
    sum_vec(:,:) = 0.0 ; sum_vec_2(:,:) = 0.0

    do j=jscq,jecq ; do i=iscq,iecq
      if (CS%umask(I,J) == 1) then
        sum_vec(I,J) = resid_scale * Zu(I,J) * Ru(I,J)
        sum_vec_2(I,J) = resid_scale * Zu_old(I,J) * Ru_old(I,J)
      endif
      if (CS%vmask(I,J) == 1) then
        sum_vec(I,J) = sum_vec(I,J) + resid_scale * Zv(I,J) * Rv(I,J)
        sum_vec_2(I,J) = sum_vec_2(I,J) + resid_scale * Zv_old(I,J) * Rv_old(I,J)
      endif
    enddo ; enddo

    beta_k = reproducing_sum(sum_vec, Is_sum, Ie_sum, Js_sum, Je_sum ) / &
             reproducing_sum(sum_vec_2, Is_sum, Ie_sum, Js_sum, Je_sum )

!    Du(:,:) = Zu(:,:) + beta_k * Du(:,:)
!    Dv(:,:) = Zv(:,:) + beta_k * Dv(:,:)

    do j=jsd,jed
      do i=isd,ied
        if (CS%umask(I,J) == 1) Du(I,J) = Zu(I,J) + beta_k * Du(I,J)
        if (CS%vmask(I,J) == 1) Dv(I,J) = Zv(I,J) + beta_k * Dv(I,J)
      enddo
    enddo

   ! D valid region moves in by 1

    sum_vec(:,:) = 0.0
    do j=jscq,jecq ; do i=iscq,iecq
      if (CS%umask(I,J) == 1) sum_vec(I,J) = resid2_scale*Ru(I,J)**2
      if (CS%vmask(I,J) == 1) sum_vec(I,J) = sum_vec(I,J) + resid2_scale*Rv(I,J)**2
    enddo ; enddo

    dot_p1 = reproducing_sum( sum_vec, Is_sum, Ie_sum, Js_sum, Je_sum )
    dot_p1 = sqrt(dot_p1)

    if (dot_p1 <= CS%cg_tolerance * resid0) then
      iters = iter
      conv_flag = 1
      exit
    endif

    cg_halo = cg_halo - 1

    if (cg_halo == 0) then
      ! pass vectors
      call pass_vector(Du, Dv, G%domain, TO_ALL, BGRID_NE)
      call pass_vector(u_shlf, v_shlf, G%domain, TO_ALL, BGRID_NE)
      call pass_vector(Ru, Rv, G%domain, TO_ALL, BGRID_NE)
      cg_halo = 3
    endif

  enddo ! end of CG loop

  do j=jsdq,jedq
    do i=isdq,iedq
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
    enddo
  enddo

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
  real :: u_face     ! Zonal velocity at a face [L Z-1 ~> m s-1]
  real :: h_face     ! Thickness at a face for transport [Z ~> m]
  real :: slope_lim  ! The value of the slope limiter, in the range of 0 to 2 [nondim]

!  is = G%isc-2 ; ie = G%iec+2 ; js = G%jsc ; je = G%jec
!  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
!  i_off = G%idg_offset ; j_off = G%jdg_offset

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh

  ! hmask coded values: 1) fully covered; 2) partly covered - no export; 3) Specified boundary condition
  ! relevant u_face_mask coded values: 1) Normal interior point; 4) Specified flux BC

  do j=jsh,jeh ; do I=ish-1,ieh
    if (CS%u_face_mask(I,j) == 4.) then ! The flux itself is a specified boundary condition.
      uh_ice(I,j) = time_step * G%dyCu(I,j) * CS%u_flux_bdry_val(I,j)
    elseif ((hmask(i,j)==1) .or. (hmask(i+1,j) == 1)) then
      u_face = 0.5 * (CS%u_shelf(I,J-1) + CS%u_shelf(I,J))
      h_face = 0.0 ! This will apply when the source cell is iceless or not fully ice covered.

      if (u_face > 0) then
        if (hmask(i,j) == 3) then ! This is a open boundary inflow from the west
          h_face = CS%thickness_bdry_val(i,j)
        elseif (hmask(i,j) == 1) then ! There can be eastward flow through this face.
          if ((hmask(i-1,j) == 1) .and. (hmask(i+1,j) == 1)) then
            slope_lim = slope_limiter(h0(i,j)-h0(i-1,j), h0(i+1,j)-h0(i,j))
            ! This is a 2nd-order centered scheme with a slope limiter.  We could try PPM here.
            h_face = h0(i,j) - slope_lim * 0.5 * (h0(i,j)-h0(i+1,j))
          else
            h_face = h0(i,j)
          endif
        endif
      else
        if (hmask(i+1,j) == 3) then ! This is a open boundary inflow from the east
          h_face = CS%thickness_bdry_val(i+1,j)
        elseif (hmask(i+1,j) == 1) then
          if ((hmask(i,j) == 1) .and. (hmask(i+2,j) == 1)) then
            slope_lim = slope_limiter(h0(i+1,j)-h0(i,j), h0(i+2,j)-h0(i+1,j))
            h_face = h0(i+1,j) - slope_lim * 0.5 * (h0(i+1,j)-h0(i,j))
          else
            h_face = h0(i+1,j)
          endif
        endif
      endif

      uh_ice(I,j) = time_step * G%dyCu(I,j) * u_face * h_face
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
  real :: v_face     ! Pseudo-meridional velocity at a face [L Z-1 ~> m s-1]
  real :: h_face     ! Thickness at a face for transport [Z ~> m]
  real :: slope_lim  ! The value of the slope limiter, in the range of 0 to 2 [nondim]

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh

  ! hmask coded values: 1) fully covered; 2) partly covered - no export; 3) Specified boundary condition
  ! relevant u_face_mask coded values: 1) Normal interior point; 4) Specified flux BC

  do J=jsh-1,jeh ; do i=ish,ieh
    if (CS%v_face_mask(i,J) == 4.) then ! The flux itself is a specified boundary condition.
      vh_ice(i,J) = time_step * G%dxCv(i,J) * CS%v_flux_bdry_val(i,J)
    elseif ((hmask(i,j)==1) .or. (hmask(i,j+1) == 1)) then

      v_face = 0.5 * (CS%v_shelf(I-1,J) + CS%v_shelf(I,J))
      h_face = 0.0 ! This will apply when the source cell is iceless or not fully ice covered.

      if (v_face > 0) then
        if (hmask(i,j) == 3) then ! This is a open boundary inflow from the south
          h_face = CS%thickness_bdry_val(i,j)
        elseif (hmask(i,j) == 1) then ! There can be northtward flow through this face.
          if ((hmask(i,j-1) == 1) .and. (hmask(i,j+1) == 1)) then
            slope_lim = slope_limiter(h0(i,j)-h0(i,j-1), h0(i,j+1)-h0(i,j))
            ! This is a 2nd-order centered scheme with a slope limiter.  We could try PPM here.
            h_face = h0(i,j) - slope_lim * 0.5 * (h0(i,j)-h0(i,j+1))
          else
            h_face = h0(i,j)
          endif
        endif
      else
        if (hmask(i,j+1) == 3) then ! This is a open boundary inflow from the north
          h_face = CS%thickness_bdry_val(i,j+1)
        elseif (hmask(i,j+1) == 1) then
          if ((hmask(i,j) == 1) .and. (hmask(i,j+2) == 1)) then
            slope_lim = slope_limiter(h0(i,j+1)-h0(i,j), h0(i,j+2)-h0(i,j+1))
            h_face = h0(i,j+1) - slope_lim * 0.5 * (h0(i,j+1)-h0(i,j))
          else
            h_face = h0(i,j+1)
          endif
        endif
      endif

      vh_ice(i,J) = time_step * G%dxCv(i,J) * v_face * h_face
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

  integer :: i, j, isc, iec, jsc, jec, n_flux, k, l, iter_count
  integer :: i_off, j_off
  integer :: iter_flag

  real :: h_reference ! A reference thicknesss based on neighboring cells [Z ~> m]
  real :: tot_flux    ! The total ice mass flux [Z L2 ~> m3]
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

    do j=jsc-1,jec+1

      if (((j+j_off) <= G%domain%njglobal+G%domain%njhalo) .AND. &
         ((j+j_off) >= G%domain%njhalo+1)) then

        do i=isc-1,iec+1

          if (((i+i_off) <= G%domain%niglobal+G%domain%nihalo) .AND. &
              ((i+i_off) >= G%domain%nihalo+1)) then
        ! first get reference thickness by averaging over cells that are fluxing into this cell
            n_flux = 0
            h_reference = 0.0
            tot_flux = 0.0

            do k=1,2
              if (flux_enter(i,j,k) > 0) then
                n_flux = n_flux + 1
                h_reference = h_reference + ISS%h_shelf(i+2*k-3,j)
                tot_flux = tot_flux + flux_enter(i,j,k)
                flux_enter(i,j,k) = 0.0
              endif
            enddo

            do k=1,2
              if (flux_enter(i,j,k+2) > 0) then
                n_flux = n_flux + 1
                h_reference = h_reference + ISS%h_shelf(i,j+2*k-3)
                tot_flux = tot_flux + flux_enter(i,j,k+2)
                flux_enter(i,j,k+2) = 0.0
              endif
            enddo

            if (n_flux > 0) then
              dxdyh = G%areaT(i,j)
              h_reference = h_reference / real(n_flux)
              partial_vol = ISS%h_shelf(i,j) * ISS%area_shelf_h(i,j) + tot_flux

              if ((partial_vol / G%areaT(i,j)) == h_reference) then ! cell is exactly covered, no overflow
                ISS%hmask(i,j) = 1
                ISS%h_shelf(i,j) = h_reference
                ISS%area_shelf_h(i,j) = G%areaT(i,j)
              elseif ((partial_vol / G%areaT(i,j)) < h_reference) then
                ISS%hmask(i,j) = 2
               !  ISS%mass_shelf(i,j) = partial_vol * CS%density_ice
                ISS%area_shelf_h(i,j) = partial_vol / h_reference
                ISS%h_shelf(i,j) = h_reference
              else

                ISS%hmask(i,j) = 1
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

subroutine calc_shelf_driving_stress(CS, ISS, G, US, taudx, taudy, OD)
  type(ice_shelf_dyn_CS), intent(in)   :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state), intent(in)    :: ISS !< A structure with elements that describe
                                             !! the ice-shelf state
  type(ocean_grid_type), intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: OD  !< ocean floor depth at tracer points [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: taudx  !< X-direction driving stress at q-points [kg L s-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: taudy  !< Y-direction driving stress at q-points [kg L s-2 ~> kg m s-2]
                                                  ! This will become [R L3 Z T-2 ~> kg m s-2]

! driving stress!

! ! taudx and taudy will hold driving stress in the x- and y- directions when done.
!    they will sit on the BGrid, and so their size depends on whether the grid is symmetric
!
! Since this is a finite element solve, they will actually have the form \int \Phi_i rho g h \nabla s
!
! OD -this is important and we do not yet know where (in MOM) it will come from. It represents
!     "average" ocean depth -- and is needed to find surface elevation
!    (it is assumed that base_ice = bed + OD)

  real, dimension(SIZE(OD,1),SIZE(OD,2))  :: S, &     ! surface elevation [Z ~> m].
                            BASE     ! basal elevation of shelf/stream [Z ~> m].


  real    :: rho, rhow ! Ice and ocean densities [R ~> kg m-3]
  real    :: sx, sy    ! Ice shelf top slopes [Z L-1 ~> m s-1]
  real    :: neumann_val ! [R Z L2 T-2 ~> kg s-2]
  real    :: dxh, dyh  ! Local grid spacing [L ~> m]
  real    :: grav      ! The gravitational acceleration [L2 Z-1 T-2 ~> m s-2]

  integer :: i, j, iscq, iecq, jscq, jecq, isd, jsd, is, js, iegq, jegq
  integer :: giec, gjec, gisc, gjsc, cnt, isc, jsc, iec, jec
  integer :: i_off, j_off

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  isd = G%isd ; jsd = G%jsd
  iegq = G%iegB ; jegq = G%jegB
  gisc = G%domain%nihalo+1 ; gjsc = G%domain%njhalo+1
  giec = G%domain%niglobal+G%domain%nihalo ; gjec = G%domain%njglobal+G%domain%njhalo
  is = iscq - 1; js = jscq - 1
  i_off = G%idg_offset ; j_off = G%jdg_offset

  rho =  CS%density_ice
  rhow = CS%density_ocean_avg
  grav = CS%g_Earth

  ! prelim - go through and calculate S

  ! or is this faster?
  BASE(:,:) = -G%bathyT(:,:) + OD(:,:)
  S(:,:) = BASE(:,:) + ISS%h_shelf(:,:)

  do j=jsc-1,jec+1
    do i=isc-1,iec+1
      cnt = 0
      sx = 0
      sy = 0
      dxh = G%dxT(i,j)
      dyh = G%dyT(i,j)

      if (ISS%hmask(i,j) == 1) then ! we are inside the global computational bdry, at an ice-filled cell

        ! calculate sx
        if ((i+i_off) == gisc) then ! at left computational bdry
          if (ISS%hmask(i+1,j) == 1) then
            sx = (S(i+1,j)-S(i,j))/dxh
          else
            sx = 0
          endif
        elseif ((i+i_off) == giec) then ! at east computational bdry
          if (ISS%hmask(i-1,j) == 1) then
            sx = (S(i,j)-S(i-1,j))/dxh
          else
            sx = 0
          endif
        else ! interior
          if (ISS%hmask(i+1,j) == 1) then
            cnt = cnt+1
            sx = S(i+1,j)
          else
            sx = S(i,j)
          endif
          if (ISS%hmask(i-1,j) == 1) then
            cnt = cnt+1
            sx = sx - S(i-1,j)
          else
            sx = sx - S(i,j)
          endif
          if (cnt == 0) then
            sx = 0
          else
            sx = sx / (cnt * dxh)
          endif
        endif

        cnt = 0

        ! calculate sy, similarly
        if ((j+j_off) == gjsc) then ! at south computational bdry
          if (ISS%hmask(i,j+1) == 1) then
            sy = (S(i,j+1)-S(i,j))/dyh
          else
            sy = 0
          endif
        elseif ((j+j_off) == gjec) then ! at nprth computational bdry
          if (ISS%hmask(i,j-1) == 1) then
            sy = (S(i,j)-S(i,j-1))/dyh
          else
            sy = 0
          endif
        else ! interior
          if (ISS%hmask(i,j+1) == 1) then
            cnt = cnt+1
            sy = S(i,j+1)
          else
            sy = S(i,j)
          endif
          if (ISS%hmask(i,j-1) == 1) then
            cnt = cnt+1
            sy = sy - S(i,j-1)
          else
            sy = sy - S(i,j)
          endif
          if (cnt == 0) then
            sy = 0
          else
            sy = sy / (cnt * dyh)
          endif
        endif

        ! SW vertex
        taudx(I-1,J-1) = taudx(I-1,J-1) - .25 * rho * grav * ISS%h_shelf(i,j) * sx * G%areaT(i,j)
        taudy(I-1,J-1) = taudy(I-1,J-1) - .25 * rho * grav * ISS%h_shelf(i,j) * sy * G%areaT(i,j)

        ! SE vertex
        taudx(I,J-1) = taudx(I,J-1) - .25 * rho * grav * ISS%h_shelf(i,j) * sx * G%areaT(i,j)
        taudy(I,J-1) = taudy(I,J-1) - .25 * rho * grav * ISS%h_shelf(i,j) * sy * G%areaT(i,j)

        ! NW vertex
        taudx(I-1,J) = taudx(I-1,J) - .25 * rho * grav * ISS%h_shelf(i,j) * sx * G%areaT(i,j)
        taudy(I-1,J) = taudy(I-1,J) - .25 * rho * grav * ISS%h_shelf(i,j) * sy * G%areaT(i,j)

        ! NE vertex
        taudx(I,J) = taudx(I,J) - .25 * rho * grav * ISS%h_shelf(i,j) * sx * G%areaT(i,j)
        taudy(I,J) = taudy(I,J) - .25 * rho * grav * ISS%h_shelf(i,j) * sy * G%areaT(i,j)

        if (CS%ground_frac(i,j) == 1) then
          neumann_val = .5 * grav * (rho * ISS%h_shelf(i,j)**2 - rhow * G%bathyT(i,j)**2)
        else
          neumann_val = .5 * grav * (1-rho/rhow) * rho * ISS%h_shelf(i,j)**2
        endif

        if ((CS%u_face_mask(I-1,j) == 2) .OR. (ISS%hmask(i-1,j) == 0) .OR. (ISS%hmask(i-1,j) == 2) ) then
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

        if ((CS%u_face_mask(I,j) == 2) .OR. (ISS%hmask(i+1,j) == 0) .OR. (ISS%hmask(i+1,j) == 2) ) then
          ! east face of the cell is at a stress boundary
          taudx(I,J-1) = taudx(I,J-1) + .5 * dyh * neumann_val
          taudx(I,J) = taudx(I,J) + .5 * dyh * neumann_val
        endif

        if ((CS%v_face_mask(i,J-1) == 2) .OR. (ISS%hmask(i,j-1) == 0) .OR. (ISS%hmask(i,j-1) == 2) ) then
          ! south face of the cell is at a stress boundary
          taudy(I-1,J-1) = taudy(I-1,J-1) - .5 * dxh * neumann_val
          taudy(I,J-1) = taudy(I,J-1) - .5 * dxh * neumann_val
        endif

        if ((CS%v_face_mask(i,J) == 2) .OR. (ISS%hmask(i,j+1) == 0) .OR. (ISS%hmask(i,j+1) == 2) ) then
          ! north face of the cell is at a stress boundary
          taudy(I-1,J) = taudy(I-1,J) + .5 * dxh * neumann_val
          taudy(I,J) = taudy(I,J) + .5 * dxh * neumann_val
        endif

      endif
    enddo
  enddo

end subroutine calc_shelf_driving_stress

subroutine init_boundary_values(CS, G, time, hmask, input_flux, input_thick, new_sim)
  type(ice_shelf_dyn_CS),intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type), intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(time_type),       intent(in)    :: Time !< The current model time
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real,                  intent(in)    :: input_flux !< The integrated inward ice thickness flux per
                                             !! unit face length [Z L T-1 ~> m2 s-1]
  real,                  intent(in)    :: input_thick !< The ice thickness at boundaries [Z ~> m].
  logical,     optional, intent(in)    :: new_sim !< If present and false, this run is being restarted

! this will be a per-setup function. the boundary values of thickness and velocity
! (and possibly other variables) will be updated in this function

! FOR RESTARTING PURPOSES: if grid is not symmetric and the model is restarted, we will
!               need to update those velocity points not *technically* in any
!               computational domain -- if this function gets moves to another module,
!               DO NOT TAKE THE RESTARTING BIT WITH IT
  integer :: i, j , isd, jsd, ied, jed
  integer :: gjec, gisc, gjsc, cnt, isc, jsc, iec, jec
  integer :: i_off, j_off

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset

  ! this loop results in some values being set twice but... eh.

  do j=jsd,jed
    do i=isd,ied

      if (hmask(i,j) == 3) then
        CS%thickness_bdry_val(i,j) = input_thick
      endif

      if ((hmask(i,j) == 0) .or. (hmask(i,j) == 1) .or. (hmask(i,j) == 2)) then
        if ((i <= iec).and.(i >= isc)) then
          if (CS%u_face_mask(I-1,j) == 3) then
            CS%u_bdry_val(I-1,J-1) = (1 - ((G%geoLatBu(I-1,J-1) - 0.5*G%len_lat)*2./G%len_lat)**2) * &
                  1.5 * input_flux / input_thick
            CS%u_bdry_val(I-1,J) = (1 - ((G%geoLatBu(I-1,J) - 0.5*G%len_lat)*2./G%len_lat)**2) * &
                  1.5 * input_flux / input_thick
          endif
        endif
      endif

      if (.not.(new_sim)) then
        if (.not. G%symmetric) then
          if (((i+i_off) == (G%domain%nihalo+1)).and.(CS%u_face_mask(I-1,j) == 3)) then
            CS%u_shelf(I-1,J-1) = CS%u_bdry_val(I-1,J-1)
            CS%u_shelf(I-1,J) = CS%u_bdry_val(I-1,J)
            CS%v_shelf(I-1,J-1) = CS%v_bdry_val(I-1,J-1)
            CS%v_shelf(I-1,J) = CS%v_bdry_val(I-1,J)
          endif
          if (((j+j_off) == (G%domain%njhalo+1)).and.(CS%v_face_mask(i,J-1) == 3)) then
            CS%u_shelf(I-1,J-1) = CS%u_bdry_val(I-1,J-1)
            CS%u_shelf(I,J-1) = CS%u_bdry_val(I,J-1)
            CS%v_shelf(I-1,J-1) = CS%v_bdry_val(I-1,J-1)
            CS%v_shelf(I,J-1) = CS%v_bdry_val(I,J-1)
          endif
        endif
      endif
    enddo
  enddo

end subroutine init_boundary_values


subroutine CG_action(uret, vret, u_shlf, v_shlf, Phi, Phisub, umask, vmask, hmask, H_node, &
                     ice_visc, float_cond, bathyT, basal_trac, G, US, is, ie, js, je, dens_ratio)

  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                         intent(inout) :: uret !< The retarding stresses working at u-points [R L3 Z T-2 ~> kg m s-2].
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                         intent(inout) :: vret !< The retarding stresses working at v-points [R L3 Z T-2 ~> kg m s-2].
  real, dimension(SZDI_(G),SZDJ_(G),8,4), &
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
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: ice_visc !< A field related to the ice viscosity from Glen's
                                               !! flow law [R L4 Z T-1 ~> kg m2 s-1]. The exact form
                                               !!  and units depend on the basal law exponent.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: float_cond !< An array indicating where the ice
                                                !! shelf is floating: 0 if floating, 1 if not.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: bathyT !< The depth of ocean bathymetry at tracer points [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: basal_trac  !< A field related to the nonlinear part of the
                                                !! "linearized" basal stress [R Z T-1 ~> kg m-2 s-1].
                ! and/or whether flow is "hybridized"
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
  integer :: iq, jq, iphi, jphi, i, j, ilq, jlq, Itgt, Jtgt
  real, dimension(2) :: xquad
  real, dimension(2,2) :: Ucell, Vcell, Hcell, Usub, Vsub

  xquad(1) = .5 * (1-sqrt(1./3)) ; xquad(2) = .5 * (1+sqrt(1./3))

  do j=js,je ; do i=is,ie ; if (hmask(i,j) == 1) then

      do iq=1,2 ; do jq=1,2

        uq = u_shlf(I-1,J-1) * xquad(3-iq) * xquad(3-jq) + &
             u_shlf(I,J-1) * xquad(iq) * xquad(3-jq) + &
             u_shlf(I-1,J) * xquad(3-iq) * xquad(jq) + &
             u_shlf(I,J) * xquad(iq) * xquad(jq)

        vq = v_shlf(I-1,J-1) * xquad(3-iq) * xquad(3-jq) + &
             v_shlf(I,J-1) * xquad(iq) * xquad(3-jq) + &
             v_shlf(I-1,J) * xquad(3-iq) * xquad(jq) + &
             v_shlf(I,J) * xquad(iq) * xquad(jq)

        ux = u_shlf(I-1,J-1) * Phi(1,2*(jq-1)+iq,i,j) + &
             u_shlf(I,J-1) * Phi(3,2*(jq-1)+iq,i,j) + &
             u_shlf(I-1,J) * Phi(5,2*(jq-1)+iq,i,j) + &
             u_shlf(I,J) * Phi(7,2*(jq-1)+iq,i,j)

        vx = v_shlf(I-1,J-1) * Phi(1,2*(jq-1)+iq,i,j) + &
             v_shlf(I,J-1) * Phi(3,2*(jq-1)+iq,i,j) + &
             v_shlf(I-1,J) * Phi(5,2*(jq-1)+iq,i,j) + &
             v_shlf(I,J) * Phi(7,2*(jq-1)+iq,i,j)

        uy = u_shlf(I-1,J-1) * Phi(2,2*(jq-1)+iq,i,j) + &
             u_shlf(I,J-1) * Phi(4,2*(jq-1)+iq,i,j) + &
             u_shlf(I-1,J) * Phi(6,2*(jq-1)+iq,i,j) + &
             u_shlf(I,J) * Phi(8,2*(jq-1)+iq,i,j)

        vy = v_shlf(I-1,j-1) * Phi(2,2*(jq-1)+iq,i,j) + &
             v_shlf(I,J-1) * Phi(4,2*(jq-1)+iq,i,j) + &
             v_shlf(I-1,J) * Phi(6,2*(jq-1)+iq,i,j) + &
             v_shlf(I,J) * Phi(8,2*(jq-1)+iq,i,j)

        do iphi=1,2 ; do jphi=1,2 ; Itgt = I-2+iphi ; Jtgt = J-2-jphi
          if (umask(Itgt,Jtgt) == 1) uret(Itgt,Jtgt) = uret(Itgt,Jtgt) + 0.25 * ice_visc(i,j) * &
                               ((4*ux+2*vy) * Phi(2*(2*(jphi-1)+iphi)-1,2*(jq-1)+iq,i,j) + &
                                    (uy+vx) * Phi(2*(2*(jphi-1)+iphi),2*(jq-1)+iq,i,j))
          if (vmask(Itgt,Jtgt) == 1) vret(Itgt,Jtgt) = vret(Itgt,Jtgt) + 0.25 * ice_visc(i,j) * &
                                   ((uy+vx) * Phi(2*(2*(jphi-1)+iphi)-1,2*(jq-1)+iq,i,j) + &
                                (4*vy+2*ux) * Phi(2*(2*(jphi-1)+iphi),2*(jq-1)+iq,i,j))

          if (float_cond(i,j) == 0) then
            ilq = 1 ; if (iq == iphi) ilq = 2
            jlq = 1 ; if (jq == jphi) jlq = 2
            if (umask(Itgt,Jtgt) == 1) uret(Itgt,Jtgt) = uret(Itgt,Jtgt) +  &
                  0.25 * basal_trac(i,j) * uq * xquad(ilq) * xquad(jlq)
            if (vmask(Itgt,Jtgt) == 1) vret(Itgt,Jtgt) = vret(Itgt,Jtgt) +  &
                  0.25 * basal_trac(i,j) * vq * xquad(ilq) * xquad(jlq)
          endif
        enddo ; enddo
      enddo ; enddo

      if (float_cond(i,j) == 1) then
        Ucell(:,:) = u_shlf(I-1:I,J-1:J) ; Vcell(:,:) = v_shlf(I-1:I,J-1:J)
        Hcell(:,:) = H_node(i-1:i,j-1:j)
        call CG_action_subgrid_basal(Phisub, Hcell, Ucell, Vcell, bathyT(i,j), dens_ratio, Usub, Vsub)

        if (umask(I-1,J-1)==1) uret(I-1,J-1) = uret(I-1,J-1) + Usub(1,1) * basal_trac(i,j)
        if (umask(I-1,J) == 1) uret(I-1,J) = uret(I-1,J) + Usub(1,2) * basal_trac(i,j)
        if (umask(I,J-1) == 1) uret(I,J-1) = uret(I,J-1) + Usub(2,1) * basal_trac(i,j)
        if (umask(I,J) == 1)   uret(I,J)   = uret(I,J) + Usub(2,2) * basal_trac(i,j)

        if (vmask(I-1,J-1)==1) vret(I-1,J-1) = vret(I-1,J-1) + Vsub(1,1) * basal_trac(i,j)
        if (vmask(I-1,J) == 1) vret(I-1,J) = vret(I-1,J) + Vsub(1,2) * basal_trac(i,j)
        if (vmask(I,J-1) == 1) vret(I,J-1) = vret(I,J-1) + Vsub(2,1) * basal_trac(i,j)
        if (vmask(I,J) == 1)   vret(I,J)   = vret(I,J) + Vsub(2,2) * basal_trac(i,j)
      endif

  endif ; enddo ; enddo

end subroutine CG_action

subroutine CG_action_subgrid_basal(Phisub, H, U, V, bathyT, dens_ratio, Ucontr, Vcontr)
  real, dimension(:,:,:,:,:,:), &
                        intent(in)    :: Phisub !< Quadrature structure weights at subgridscale
                                            !! locations for finite element calculations [nondim]
  real, dimension(2,2), intent(in)    :: H  !< The ice shelf thickness at nodal (corner) points [Z ~> m].
  real, dimension(2,2), intent(in)    :: U  !< The zonal ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(2,2), intent(in)    :: V  !< The meridional ice shelf velocity at vertices [L T-1 ~> m s-1]
  real,                 intent(in)    :: bathyT !< The depth of ocean bathymetry at tracer points [Z ~> m].
  real,                 intent(in)    :: dens_ratio !< The density of ice divided by the density
                                            !! of seawater [nondim]
  real, dimension(2,2), intent(out)   :: Ucontr !< The areal average of u-velocities where the ice shelf
                                            !! is grounded, or 0 where it is floating [L T-1 ~> m s-1].
  real, dimension(2,2), intent(out)   :: Vcontr !< The areal average of v-velocities where the ice shelf
                                            !! is grounded, or 0 where it is floating [L T-1 ~> m s-1].

  real    :: subarea ! The fractional sub-cell area [nondim]
  real    :: hloc    ! The local sub-cell ice thickness [Z ~> m]
  integer :: nsub, i, j, qx, qy, m, n

  nsub = size(Phisub,1)
  subarea = 1.0 / (nsub**2)

  do n=1,2 ; do m=1,2
    Ucontr(m,n) = 0.0 ; Vcontr(m,n) = 0.0
    do qy=1,2 ; do qx=1,2 ; do j=1,nsub ; do i=1,nsub
      hloc = (Phisub(i,j,1,1,qx,qy)*H(1,1) + Phisub(i,j,2,2,qx,qy)*H(2,2)) + &
             (Phisub(i,j,1,2,qx,qy)*H(1,2) + Phisub(i,j,2,1,qx,qy)*H(2,1))
      if (dens_ratio * hloc - bathyT > 0) then
        Ucontr(m,n) = Ucontr(m,n) + subarea * 0.25 * Phisub(i,j,m,n,qx,qy) * &
             ((Phisub(i,j,1,1,qx,qy) * U(1,1) + Phisub(i,j,2,2,qx,qy) * U(2,2)) + &
              (Phisub(i,j,1,2,qx,qy) * U(1,2) + Phisub(i,j,2,1,qx,qy) * U(2,1)))
        Vcontr(m,n) = Vcontr(m,n) + subarea * 0.25 * Phisub(i,j,m,n,qx,qy) * &
             ((Phisub(i,j,1,1,qx,qy) * V(1,1) + Phisub(i,j,2,2,qx,qy) * V(2,2)) + &
              (Phisub(i,j,1,2,qx,qy) * V(1,2) + Phisub(i,j,2,1,qx,qy) * V(2,1)))
      endif
    enddo ; enddo ; enddo ; enddo
  enddo ; enddo

end subroutine CG_action_subgrid_basal

!> returns the diagonal entries of the matrix for a Jacobi preconditioning
subroutine matrix_diagonal(CS, G, US, float_cond, H_node, ice_visc, basal_trac, hmask, dens_ratio, &
                           Phisub, u_diagonal, v_diagonal)

  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: float_cond !< An array indicating where the ice
                                                !! shelf is floating: 0 if floating, 1 if not.
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal
                                                 !! (corner) points [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: ice_visc !< A field related to the ice viscosity from Glen's
                                                !! flow law [R L4 Z T-1 ~> kg m2 s-1]. The exact form
                                                !!  and units depend on the basal law exponent.
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: basal_trac !< A field related to the nonlinear part of the
                                                !! "linearized" basal stress [R Z T-1 ~> kg m-2 s-1].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real,                   intent(in)    :: dens_ratio !< The density of ice divided by the density
                                                     !! of seawater [nondim]
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
  real, dimension(8,4) :: Phi ! Weight gradients [L-1 ~> m-1]
  real, dimension(2)   :: xquad
  real, dimension(2,2) :: Hcell, sub_ground
  integer :: i, j, is, js, cnt, isc, jsc, iec, jec, iphi, jphi, iq, jq, ilq, jlq, Itgt, Jtgt

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  xquad(1) = .5 * (1-sqrt(1./3)) ; xquad(2) = .5 * (1+sqrt(1./3))

  do j=jsc-1,jec+1 ; do i=isc-1,iec+1 ; if (hmask(i,j) == 1) then

    call bilinear_shape_fn_grid(G, i, j, Phi)

    ! Phi(2*i-1,j) gives d(Phi_i)/dx at quadrature point j
    ! Phi(2*i,j) gives d(Phi_i)/dy at quadrature point j

    do iq=1,2 ; do jq=1,2 ; do iphi=1,2 ; do jphi=1,2 ; Itgt = I-2+iphi ; Jtgt = J-2-jphi
      ilq = 1 ; if (iq == iphi) ilq = 2
      jlq = 1 ; if (jq == jphi) jlq = 2

      if (CS%umask(Itgt,Jtgt) == 1) then

        ux = Phi(2*(2*(jphi-1)+iphi)-1, 2*(jq-1)+iq)
        uy = Phi(2*(2*(jphi-1)+iphi), 2*(jq-1)+iq)
        vx = 0.
        vy = 0.

        u_diagonal(Itgt,Jtgt) = u_diagonal(Itgt,Jtgt) + &
              0.25 * ice_visc(i,j) * ((4*ux+2*vy) * Phi(2*(2*(jphi-1)+iphi)-1,2*(jq-1)+iq) + &
                                      (uy+vy) * Phi(2*(2*(jphi-1)+iphi),2*(jq-1)+iq))

        if (float_cond(i,j) == 0) then
          uq = xquad(ilq) * xquad(jlq)
          u_diagonal(Itgt,Jtgt) = u_diagonal(Itgt,Jtgt) + &
              0.25 * basal_trac(i,j) * uq * xquad(ilq) * xquad(jlq)
        endif
      endif

      if (CS%vmask(Itgt,Jtgt) == 1) then

        vx = Phi(2*(2*(jphi-1)+iphi)-1, 2*(jq-1)+iq)
        vy = Phi(2*(2*(jphi-1)+iphi), 2*(jq-1)+iq)
        ux = 0.
        uy = 0.

        v_diagonal(Itgt,Jtgt) = v_diagonal(Itgt,Jtgt) + &
              0.25 * ice_visc(i,j) * ((uy+vx) * Phi(2*(2*(jphi-1)+iphi)-1,2*(jq-1)+iq) + &
                                  (4*vy+2*ux) * Phi(2*(2*(jphi-1)+iphi),2*(jq-1)+iq))

        if (float_cond(i,j) == 0) then
          vq = xquad(ilq) * xquad(jlq)
          v_diagonal(Itgt,Jtgt) = v_diagonal(Itgt,Jtgt) + &
                0.25 * basal_trac(i,j) * vq * xquad(ilq) * xquad(jlq)
        endif
      endif
    enddo ; enddo ; enddo ; enddo

    if (float_cond(i,j) == 1) then
      Hcell(:,:) = H_node(i-1:i,j-1:j)
      call CG_diagonal_subgrid_basal(Phisub, Hcell, G%bathyT(i,j), dens_ratio, sub_ground)
      do iphi=1,2 ; do jphi=1,2 ; Itgt = I-2+iphi ; Jtgt = J-2-jphi
        if (CS%umask(Itgt,Jtgt) == 1) then
          u_diagonal(Itgt,Jtgt) = u_diagonal(Itgt,Jtgt) + sub_ground(iphi,jphi) * basal_trac(i,j)
          v_diagonal(Itgt,Jtgt) = v_diagonal(Itgt,Jtgt) + sub_ground(iphi,jphi) * basal_trac(i,j)
        endif
      enddo ; enddo
    endif
  endif ; enddo ; enddo

end subroutine matrix_diagonal

subroutine CG_diagonal_subgrid_basal (Phisub, H_node, bathyT, dens_ratio, sub_grnd)
  real, dimension(:,:,:,:,:,:), &
                        intent(in) :: Phisub !< Quadrature structure weights at subgridscale
                                             !! locations for finite element calculations [nondim]
  real, dimension(2,2), intent(in) :: H_node !< The ice shelf thickness at nodal (corner)
                                             !! points [Z ~> m].
  real,              intent(in)    :: bathyT !< The depth of ocean bathymetry at tracer points [Z ~> m].
  real,              intent(in)    :: dens_ratio !< The density of ice divided by the density
                                                 !! of seawater [nondim]
  real, dimension(2,2), intent(out) :: sub_grnd !< The weighted fraction of the sub-cell where the ice shelf
                                                !! is grounded [nondim]

  ! bathyT = cellwise-constant bed elevation

  real    :: subarea ! The fractional sub-cell area [nondim]
  real    :: hloc    ! The local sub-region thickness [Z ~> m]
  integer :: nsub, i, j, k, l, qx, qy, m, n

  nsub = size(Phisub,1)
  subarea = 1.0 / (nsub**2)

  sub_grnd(:,:) = 0.0
  do m=1,2 ; do n=1,2 ; do j=1,nsub ; do i=1,nsub ; do qx=1,2 ; do qy = 1,2

    hloc = (Phisub(i,j,1,1,qx,qy)*H_node(1,1) + Phisub(i,j,2,2,qx,qy)*H_node(2,2)) + &
           (Phisub(i,j,1,2,qx,qy)*H_node(1,2) + Phisub(i,j,2,1,qx,qy)*H_node(2,1))

    if (dens_ratio * hloc - bathyT > 0) then
      sub_grnd(m,n) = sub_grnd(m,n) + subarea * 0.25 * Phisub(i,j,m,n,qx,qy)**2
    endif

  enddo ; enddo ; enddo ; enddo ; enddo ; enddo

end subroutine CG_diagonal_subgrid_basal


subroutine apply_boundary_values(CS, ISS, G, US, time, Phisub, H_node, ice_visc, basal_trac, float_cond, &
                                 dens_ratio, u_bdry_contr, v_bdry_contr)

  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  type(time_type),        intent(in)    :: Time !< The current model time
  real, dimension(:,:,:,:,:,:), &
                          intent(in)    :: Phisub !< Quadrature structure weights at subgridscale
                                            !! locations for finite element calculations [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal
                                                 !! (corner) points [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: ice_visc !< A field related to the ice viscosity from Glen's
                                                !! flow law. The exact form and units depend on the
                                                !! basal law exponent.  [R L4 Z T-1 ~> kg m2 s-1].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: basal_trac !< A field related to the nonlinear part of the
                                                !! "linearized" basal stress [R Z T-1 ~> kg m-2 s-1].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: float_cond !< An array indicating where the ice
                                                !! shelf is floating: 0 if floating, 1 if not.
  real,                   intent(in)    :: dens_ratio !< The density of ice divided by the density
                                                     !! of seawater, nondimensional
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_bdry_contr !< Zonal force contributions due to the
                                                        !! open boundaries [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_bdry_contr !< Meridional force contributions due to the
                                                        !! open boundaries [R L3 Z T-2 ~> kg m s-2]

! this will be a per-setup function. the boundary values of thickness and velocity
! (and possibly other variables) will be updated in this function

  real, dimension(8,4)  :: Phi
  real, dimension(2) :: xquad
  real :: ux, uy, vx, vy ! Components of velocity shears or divergence [T-1 ~> s-1]
  real :: uq, vq  ! Interpolated velocities [L T-1 ~> m s-1]
  real :: area
  real, dimension(2,2) :: Ucell,Vcell,Hcell,Usubcontr,Vsubcontr
  integer :: i, j, isc, jsc, iec, jec, iq, jq, iphi, jphi, ilq, jlq, Itgt, Jtgt

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  xquad(1) = .5 * (1-sqrt(1./3)) ; xquad(2) = .5 * (1+sqrt(1./3))

  do j=jsc-1,jec+1 ; do i=isc-1,iec+1 ; if (ISS%hmask(i,j) == 1) then

    ! process this cell if any corners have umask set to non-dirichlet bdry.
    ! NOTE: vmask not considered, probably should be

    if ((CS%umask(I-1,J-1) == 3) .OR. (CS%umask(I,J-1) == 3) .OR. &
        (CS%umask(I-1,J) == 3) .OR. (CS%umask(I,J) == 3)) then

      call bilinear_shape_fn_grid(G, i, j, Phi)

      ! Phi(2*i-1,j) gives d(Phi_i)/dx at quadrature point j
      ! Phi(2*i,j) gives d(Phi_i)/dy at quadrature point j

      do iq=1,2 ; do jq=1,2

        uq = CS%u_bdry_val(I-1,J-1) * xquad(3-iq) * xquad(3-jq) + &
             CS%u_bdry_val(I,J-1) * xquad(iq) * xquad(3-jq) + &
             CS%u_bdry_val(I-1,J) * xquad(3-iq) * xquad(jq) + &
             CS%u_bdry_val(I,J) * xquad(iq) * xquad(jq)

        vq = CS%v_bdry_val(I-1,J-1) * xquad(3-iq) * xquad(3-jq) + &
             CS%v_bdry_val(I,J-1) * xquad(iq) * xquad(3-jq) + &
             CS%v_bdry_val(I-1,J) * xquad(3-iq) * xquad(jq) + &
             CS%v_bdry_val(I,J) * xquad(iq) * xquad(jq)

        ux = CS%u_bdry_val(I-1,J-1) * Phi(1,2*(jq-1)+iq) + &
             CS%u_bdry_val(I,J-1) * Phi(3,2*(jq-1)+iq) + &
             CS%u_bdry_val(I-1,J) * Phi(5,2*(jq-1)+iq) + &
             CS%u_bdry_val(I,J) * Phi(7,2*(jq-1)+iq)

        vx = CS%v_bdry_val(I-1,J-1) * Phi(1,2*(jq-1)+iq) + &
             CS%v_bdry_val(I,J-1) * Phi(3,2*(jq-1)+iq) + &
             CS%v_bdry_val(I-1,J) * Phi(5,2*(jq-1)+iq) + &
             CS%v_bdry_val(I,J) * Phi(7,2*(jq-1)+iq)

        uy = CS%u_bdry_val(I-1,J-1) * Phi(2,2*(jq-1)+iq) + &
             CS%u_bdry_val(I,J-1) * Phi(4,2*(jq-1)+iq) + &
             CS%u_bdry_val(I-1,J) * Phi(6,2*(jq-1)+iq) + &
             CS%u_bdry_val(I,J) * Phi(8,2*(jq-1)+iq)

        vy = CS%v_bdry_val(I-1,J-1) * Phi(2,2*(jq-1)+iq) + &
             CS%v_bdry_val(I,J-1) * Phi(4,2*(jq-1)+iq) + &
             CS%v_bdry_val(I-1,J) * Phi(6,2*(jq-1)+iq) + &
             CS%v_bdry_val(I,J) * Phi(8,2*(jq-1)+iq)

        do iphi=1,2 ; do jphi=1,2 ; Itgt = I-2+iphi ; Jtgt = J-2-jphi
          ilq = 1 ; if (iq == iphi) ilq = 2
          jlq = 1 ; if (jq == jphi) jlq = 2

          if (CS%umask(Itgt,Jtgt) == 1) then
            u_bdry_contr(Itgt,Jtgt) = u_bdry_contr(Itgt,Jtgt) + &
               0.25 * ice_visc(i,j) * ( (4*ux+2*vy) * Phi(2*(2*(jphi-1)+iphi)-1,2*(jq-1)+iq) + &
                                            (uy+vx) * Phi(2*(2*(jphi-1)+iphi),2*(jq-1)+iq) )

            if (float_cond(i,j) == 0) then
              u_bdry_contr(Itgt,Jtgt) = u_bdry_contr(Itgt,Jtgt) + &
                0.25 * basal_trac(i,j) * uq * xquad(ilq) * xquad(jlq)
            endif
          endif

          if (CS%vmask(Itgt,Jtgt) == 1) then
            v_bdry_contr(Itgt,Jtgt) = v_bdry_contr(Itgt,Jtgt) + &
                0.25 *  ice_visc(i,j) * ( (uy+vx) * Phi(2*(2*(jphi-1)+iphi)-1,2*(jq-1)+iq) + &
                                      (4*vy+2*ux) * Phi(2*(2*(jphi-1)+iphi),2*(jq-1)+iq) )

            if (float_cond(i,j) == 0) then
              v_bdry_contr(Itgt,Jtgt) = v_bdry_contr(Itgt,Jtgt) + &
                  0.25 * basal_trac(i,j) * vq * xquad(ilq) * xquad(jlq)
            endif
          endif
        enddo ; enddo
      enddo ; enddo

      if (float_cond(i,j) == 1) then
        Ucell(:,:) = CS%u_bdry_val(i-1:i,j-1:j) ; Vcell(:,:) = CS%v_bdry_val(i-1:i,j-1:j)
        Hcell(:,:) = H_node(i-1:i,j-1:j)
        call CG_action_subgrid_basal(Phisub, Hcell, Ucell, Vcell, G%bathyT(i,j), &
                                     dens_ratio, Usubcontr, Vsubcontr)

        if (CS%umask(I-1,J-1)==1) u_bdry_contr(I-1,J-1) = u_bdry_contr(I-1,J-1) + Usubcontr(1,1) * basal_trac(i,j)
        if (CS%umask(I-1,J) == 1) u_bdry_contr(I-1,J) = u_bdry_contr(I-1,J) + Usubcontr(1,2) * basal_trac(i,j)
        if (CS%umask(I,J-1) == 1) u_bdry_contr(I,J-1) = u_bdry_contr(I,J-1) + Usubcontr(2,1) * basal_trac(i,j)
        if (CS%umask(I,J) == 1)   u_bdry_contr(I,J)   = u_bdry_contr(I,J) + Usubcontr(2,2) * basal_trac(i,j)

        if (CS%vmask(I-1,J-1)==1) v_bdry_contr(I-1,J-1) = v_bdry_contr(I-1,J-1) + Vsubcontr(1,1) * basal_trac(i,j)
        if (CS%vmask(I-1,J) == 1) v_bdry_contr(I-1,J) = v_bdry_contr(I-1,J) + Vsubcontr(1,2) * basal_trac(i,j)
        if (CS%vmask(I,J-1) == 1) v_bdry_contr(I,J-1) = v_bdry_contr(I,J-1) + Vsubcontr(2,1) * basal_trac(i,j)
        if (CS%vmask(I,J) == 1)   v_bdry_contr(I,J)   = v_bdry_contr(I,J) + Vsubcontr(2,2) * basal_trac(i,j)
      endif
    endif
  endif ; enddo ; enddo

end subroutine apply_boundary_values

!> Update depth integrated viscosity, based on horizontal strain rates, and also update the
!! nonlinear part of the basal traction.
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
! so there is an "upper" and "lower" bilinear viscosity

! also this subroutine updates the nonlinear part of the basal traction

! this may be subject to change later... to make it "hybrid"

  integer :: i, j, iscq, iecq, jscq, jecq, isd, jsd, ied, jed, iegq, jegq
  integer :: giec, gjec, gisc, gjsc, cnt, isc, jsc, iec, jec, is, js
  real :: Visc_coef, n_g
  real :: ux, uy, vx, vy, eps_min ! Velocity shears [T-1 ~> s-1]
  real :: umid, vmid, unorm ! Velocities [L T-1 ~> m s-1]

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
  iegq = G%iegB ; jegq = G%jegB
  gisc = G%domain%nihalo+1 ; gjsc = G%domain%njhalo+1
  giec = G%domain%niglobal+gisc ; gjec = G%domain%njglobal+gjsc
  is = iscq - 1; js = jscq - 1

  n_g = CS%n_glen; eps_min = CS%eps_glen_min

  Visc_coef = US%kg_m2s_to_RZ_T*US%m_to_L*US%Z_to_L*(CS%A_glen_isothermal)**(1./CS%n_glen)

  do j=jsd+1,jed-1
    do i=isd+1,ied-1

      if (ISS%hmask(i,j) == 1) then
        ux = ((u_shlf(I,J) + u_shlf(I,J-1)) - (u_shlf(I-1,J) + u_shlf(I-1,J-1))) / (2*G%dxT(i,j))
        vx = ((v_shlf(I,J) + v_shlf(I,J-1)) - (v_shlf(I-1,J) + v_shlf(I-1,J-1))) / (2*G%dxT(i,j))
        uy = ((u_shlf(I,J) + u_shlf(I-1,J)) - (u_shlf(I,J-1) + u_shlf(I-1,J-1))) / (2*G%dyT(i,j))
        vy = ((v_shlf(I,J) + v_shlf(I-1,J)) - (v_shlf(I,J-1) + v_shlf(I-1,J-1))) / (2*G%dyT(i,j))
        CS%ice_visc(i,j) = 0.5 * Visc_coef * (G%areaT(i,j) * ISS%h_shelf(i,j)) * &
             (US%s_to_T**2 * (ux**2 + vy**2 + ux*vy + 0.25*(uy+vx)**2 + eps_min**2))**((1.-n_g)/(2.*n_g))

        umid = ((u_shlf(I,J) + u_shlf(I-1,J-1)) + (u_shlf(I,J-1) + u_shlf(I-1,J))) * 0.25
        vmid = ((v_shlf(I,J) + v_shlf(I-1,J-1)) + (v_shlf(I,J-1) + v_shlf(I-1,J))) * 0.25
        unorm = sqrt(umid**2 + vmid**2 + eps_min**2*(G%dxT(i,j)**2 + G%dyT(i,j)**2))
        CS%basal_traction(i,j) = G%areaT(i,j) * CS%C_basal_friction * (US%L_T_to_m_s*unorm)**(CS%n_basal_fric-1)
      endif
    enddo
  enddo

end subroutine calc_shelf_visc

subroutine update_OD_ffrac(CS, G, US, ocean_mass, find_avg)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type), intent(in)     :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: ocean_mass !< The mass per unit area of the ocean [kg m-2].
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

      CS%OD_rt(i,j) = 0.0 ; CS%ground_frac_rt(i,j) = 0.0
    enddo ; enddo

    call pass_var(CS%ground_frac, G%domain)
    call pass_var(CS%OD_av, G%domain)
  endif

end subroutine update_OD_ffrac

subroutine update_OD_ffrac_uncoupled(CS, G, h_shelf)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_shelf !< the thickness of the ice shelf [Z ~> m].

  integer :: i, j, iters, isd, ied, jsd, jed
  real    :: rhoi_rhow, OD

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j=jsd,jed
    do i=isd,ied
      OD = G%bathyT(i,j) - rhoi_rhow * h_shelf(i,j)
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
  integer :: node, qpoint, xnode, xq, ynode, yq

  xquad(1:3:2) = .5 * (1-sqrt(1./3)) ; yquad(1:2) = .5 * (1-sqrt(1./3))
  xquad(2:4:2) = .5 * (1+sqrt(1./3)) ; yquad(3:4) = .5 * (1+sqrt(1./3))

  do qpoint=1,4

    a = -X(1)*(1-yquad(qpoint)) + X(2)*(1-yquad(qpoint)) - X(3)*yquad(qpoint) + X(4)*yquad(qpoint) ! d(x)/d(x*)
    b = -Y(1)*(1-yquad(qpoint)) + Y(2)*(1-yquad(qpoint)) - Y(3)*yquad(qpoint) + Y(4)*yquad(qpoint) ! d(y)/d(x*)
    c = -X(1)*(1-xquad(qpoint)) - X(2)*(xquad(qpoint)) + X(3)*(1-xquad(qpoint)) + X(4)*(xquad(qpoint)) ! d(x)/d(y*)
    d = -Y(1)*(1-xquad(qpoint)) - Y(2)*(xquad(qpoint)) + Y(3)*(1-xquad(qpoint)) + Y(4)*(xquad(qpoint)) ! d(y)/d(y*)

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

      Phi(2*node-1,qpoint) = ( d * (2 * xnode - 3) * yexp - b * (2 * ynode - 3) * xexp) / (a*d-b*c)
      Phi(2*node,qpoint)   = (-c * (2 * xnode - 3) * yexp + a * (2 * ynode - 3) * xexp) / (a*d-b*c)

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
  integer :: node, qpoint, xnode, xq, ynode, yq

  xquad(1:3:2) = .5 * (1-sqrt(1./3)) ; yquad(1:2) = .5 * (1-sqrt(1./3))
  xquad(2:4:2) = .5 * (1+sqrt(1./3)) ; yquad(3:4) = .5 * (1+sqrt(1./3))

  do qpoint=1,4
    a = G%dxCv(i,J-1) * (1-yquad(qpoint)) + G%dxCv(i,J) * yquad(qpoint) ! d(x)/d(x*)
    d = G%dyCu(I-1,j) * (1-xquad(qpoint)) + G%dyCu(I,j) * xquad(qpoint) ! d(y)/d(y*)

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

      Phi(2*node-1,qpoint) = ( d * (2 * xnode - 3) * yexp ) / (a*d)
      Phi(2*node,qpoint)   = ( a * (2 * ynode - 3) * xexp ) / (a*d)

    enddo
  enddo

end subroutine bilinear_shape_fn_grid


subroutine bilinear_shape_functions_subgrid(Phisub, nsub)
  real, dimension(nsub,nsub,2,2,2,2), &
           intent(inout) :: Phisub !< Quadrature structure weights at subgridscale
                                   !! locations for finite element calculations [nondim]
  integer, intent(in)    :: nsub   !< The number of subgridscale quadrature locations in each direction

  ! this subroutine is a helper for interpolation of floatation condition
  ! for the purposes of evaluating the terms \int (u,v) \phi_i dx dy in a cell that is
  !     in partial floatation
  ! the array Phisub contains the values of \phi_i (where i is a node of the cell)
  !     at quad point j
  ! i think this general approach may not work for nonrectangular elements...
  !

  ! Phisub(i,j,k,l,q1,q2)
  !  i: subgrid index in x-direction
  !  j: subgrid index in y-direction
  !  k: basis function x-index
  !  l: basis function y-index
  !  q1: quad point x-index
  !  q2: quad point y-index

  ! e.g. k=1,l=1 => node 1
  !      q1=2,q2=1 => quad point 2

    !  3 - 4
    !  |   |
    !  1 - 2

  integer :: i, j, k, l, qx, qy, indx, indy
  real,dimension(2)    :: xquad
  real                 :: x0, y0, x, y, val, fracx

  xquad(1) = .5 * (1-sqrt(1./3)) ; xquad(2) = .5 * (1+sqrt(1./3))
  fracx = 1.0/real(nsub)

  do j=1,nsub ; do i=1,nsub
    x0 = (i-1) * fracx ; y0 = (j-1) * fracx
    do qy=1,2 ; do qx=1,2
      x = x0 + fracx*xquad(qx)
      y = y0 + fracx*xquad(qy)
      Phisub(i,j,1,1,qx,qy) = (1.0-x) * (1.0-y)
      Phisub(i,j,1,2,qx,qy) = (1.0-x) * y
      Phisub(i,j,2,1,qx,qy) = x * (1.0-y)
      Phisub(i,j,2,2,qx,qy) = x * y
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
  real, dimension(SZDIB_(G),SZDJ_(G)), &
                         intent(out)   :: u_face_mask !< A coded mask for velocities at the C-grid u-face
  real, dimension(SZDI_(G),SZDJB_(G)), &
                         intent(out)   :: v_face_mask !< A coded mask for velocities at the C-grid v-face
  ! sets masks for velocity solve
  ! ignores the fact that their might be ice-free cells - this only considers the computational boundary

  ! !!!IMPORTANT!!! relies on thickness mask - assumed that this is called after hmask has been updated & halo-updated

  integer :: i, j, k, iscq, iecq, jscq, jecq, isd, jsd, is, js, iegq, jegq
  integer :: giec, gjec, gisc, gjsc, isc, jsc, iec, jec
  integer :: i_off, j_off

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  i_off = G%idg_offset ; j_off = G%jdg_offset
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

  do j=js,G%jed
    do i=is,G%ied

      if (hmask(i,j) == 1) then

        umask(I-1:I,j-1:j) = 1.
        vmask(I-1:I,j-1:j) = 1.

        do k=0,1

          select case (int(CS%u_face_mask_bdry(I-1+k,j)))
            case (3)
              umask(I-1+k,J-1:J)=3.
              vmask(I-1+k,J-1:J)=0.
              u_face_mask(I-1+k,j)=3.
            case (2)
              u_face_mask(I-1+k,j)=2.
            case (4)
              umask(I-1+k,J-1:J)=0.
              vmask(I-1+k,J-1:J)=0.
              u_face_mask(I-1+k,j)=4.
            case (0)
              umask(I-1+k,J-1:J)=0.
              vmask(I-1+k,J-1:J)=0.
              u_face_mask(I-1+k,j)=0.
            case (1)  ! stress free x-boundary
              umask(I-1+k,J-1:J)=0.
            case default
          end select
        enddo

        do k=0,1

          select case (int(CS%v_face_mask_bdry(i,J-1+k)))
            case (3)
              vmask(I-1:I,J-1+k)=3.
              umask(I-1:I,J-1+k)=0.
              v_face_mask(i,J-1+k)=3.
            case (2)
              v_face_mask(i,J-1+k)=2.
            case (4)
              umask(I-1:I,J-1+k)=0.
              vmask(I-1:I,J-1+k)=0.
              v_face_mask(i,J-1+k)=4.
            case (0)
              umask(I-1:I,J-1+k)=0.
              vmask(I-1:I,J-1+k)=0.
              v_face_mask(i,J-1+k)=0.
            case (1) ! stress free y-boundary
              vmask(I-1:I,J-1+k)=0.
            case default
          end select
        enddo

        !if (CS%u_face_mask_bdry(I-1,j) >= 0) then ! Western boundary
        !  u_face_mask(I-1,j) = CS%u_face_mask_bdry(I-1,j)
        !  umask(I-1,J-1:J) = 3.
        !  vmask(I-1,J-1:J) = 0.
        !endif

        !if (j_off+j == gjsc+1) then ! SoutherN boundary
        !  v_face_mask(i,J-1) = 0.
        !  umask(I-1:I,J-1) = 0.
        !  vmask(I-1:I,J-1) = 0.
        !elseif (j_off+j == gjec) then ! Northern boundary
        !  v_face_mask(i,J) = 0.
        !  umask(I-1:I,J) = 0.
        !  vmask(I-1:I,J) = 0.
        !endif

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
subroutine interpolate_H_to_B(G, h_shelf, hmask, H_node)
  type(ocean_grid_type), intent(inout) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: h_shelf !< The ice shelf thickness at tracer points [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: H_node !< The ice shelf thickness at nodal (corner)
                                             !! points [Z ~> m].

  integer :: i, j, isc, iec, jsc, jec, num_h, k, l
  real    :: summ

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  H_node(:,:) = 0.0

  ! H_node is node-centered; average over all cells that share that node
  ! if no (active) cells share the node then its value there is irrelevant

  do j=jsc-1,jec
    do i=isc-1,iec
      summ = 0.0
      num_h = 0
      do k=0,1
        do l=0,1
          if (hmask(i+k,j+l) == 1.0) then
            summ = summ + h_shelf(i+k,j+l)
            num_h = num_h + 1
          endif
        enddo
      enddo
      if (num_h > 0) then
        H_node(i,j) = summ / num_h
      endif
    enddo
  enddo

  call pass_var(H_node, G%domain, position=CORNER)

end subroutine interpolate_H_to_B

!> Deallocates all memory associated with the ice shelf dynamics module
subroutine ice_shelf_dyn_end(CS)
  type(ice_shelf_dyn_CS), pointer   :: CS !< A pointer to the ice shelf dynamics control structure

  if (.not.associated(CS)) return

  deallocate(CS%u_shelf, CS%v_shelf)
  deallocate(CS%t_shelf, CS%tmask)
  deallocate(CS%u_bdry_val, CS%v_bdry_val, CS%t_bdry_val)
  deallocate(CS%u_face_mask, CS%v_face_mask)
  deallocate(CS%umask, CS%vmask)

  deallocate(CS%ice_visc, CS%basal_traction)
  deallocate(CS%OD_rt, CS%OD_av)
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

! 5/23/12 OVS
!    This subroutine takes the velocity (on the Bgrid) and timesteps
!      (HT)_t = - div (uHT) + (adot Tsurf -bdot Tbot) once and then calculates T=HT/H
!
!    The flux overflows are included here. That is because they will be used to advect 3D scalars
!    into partial cells

  real, dimension(SZDI_(G),SZDJ_(G))   :: th_after_uflux, th_after_vflux, TH
  integer                           :: isd, ied, jsd, jed, i, j, isc, iec, jsc, jec
  real :: Tsurf ! Surface air temperature.  This is hard coded but should be an input argument.
  real :: adot  ! A surface heat exchange coefficient divided by the heat capacity of
                ! ice [R Z T-1 degC-1 ~> kg m-2 s-1 degC-1].


  ! For now adot and Tsurf are defined here adot=surf acc 0.1m/yr, Tsurf=-20oC, vary them later
  adot = (0.1/(365.0*86400.0))*US%m_to_Z*US%T_to_s * CS%density_ice
  Tsurf = -20.0

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

!  call enable_averages(time_step, Time, CS%diag)
!  call pass_var(h_after_uflux, G%domain)
!  call pass_var(h_after_vflux, G%domain)
!  if (CS%id_h_after_uflux > 0) call post_data(CS%id_h_after_uflux, h_after_uflux, CS%diag)
!  if (CS%id_h_after_vflux > 0) call post_data(CS%id_h_after_vflux, h_after_vflux, CS%diag)
!  call disable_averaging(CS%diag)

  call ice_shelf_advect_temp_x(CS, G, time_step, ISS%hmask, TH, th_after_uflux)
  call ice_shelf_advect_temp_y(CS, G, time_step, ISS%hmask, th_after_uflux, th_after_vflux)

  do j=jsc,jec ; do i=isc,iec
    ! Convert the integrated temperature back to the average temperature.
!   if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 2)) then
    if (ISS%h_shelf(i,j) > 0.0) then
      CS%t_shelf(i,j) = th_after_vflux(i,j) / ISS%h_shelf(i,j)
    else
      CS%t_shelf(i,j) = -10.0
    endif
!   endif

    if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 2)) then
      if (ISS%h_shelf(i,j) > 0.0) then
        CS%t_shelf(i,j) = CS%t_shelf(i,j) + &
            time_step*(adot*Tsurf - melt_rate(i,j)*ISS%tfreeze(i,j))/(CS%density_ice*ISS%h_shelf(i,j))
      else
        ! the ice is about to melt away in this case set thickness, area, and mask to zero
        ! NOTE: not mass conservative, should maybe scale salt & heat flux for this cell
        CS%t_shelf(i,j) = -10.0
        CS%tmask(i,j) = 0.0
      endif
    elseif (ISS%hmask(i,j) == 0) then
      CS%t_shelf(i,j) = -10.0
    elseif ((ISS%hmask(i,j) == 3) .or. (ISS%hmask(i,j) == -2)) then
      CS%t_shelf(i,j) = CS%t_bdry_val(i,j)
    endif
  enddo ; enddo

  call pass_var(CS%t_shelf, G%domain)
  call pass_var(CS%tmask, G%domain)

  if (CS%debug) then
    call hchksum(CS%t_shelf, "temp after front", G%HI, haloshift=3)
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
                          intent(in)    :: h0 !< The initial ice shelf thicknesses [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_uflux !< The ice shelf thicknesses after
                                              !! the zonal mass fluxes [Z ~> m].

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization

  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, gjed, gied
  integer :: i_off, j_off
  logical :: at_east_bdry, at_west_bdry, one_off_west_bdry, one_off_east_bdry
  real, dimension(-2:2) :: stencil
  real :: u_face     ! Zonal velocity at a face, positive if out {L T-1 ~> m s-1]
  real :: flux_diff, phi

  character (len=1)        :: debug_str


  is = G%isc-2 ; ie = G%iec+2 ; js = G%jsc ; je = G%jec ; isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset

  do j=jsd+1,jed-1
    if (((j+j_off) <= G%domain%njglobal+G%domain%njhalo) .AND. &
        ((j+j_off) >= G%domain%njhalo+1)) then ! based on mehmet's code - only if btw north & south boundaries

      stencil(:) = -1
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
                  flux_diff = flux_diff + ABS(u_face) * G%dyCu(I-1,j)* time_step / G%areaT(i,j) * &
                           (stencil(-1) - phi * (stencil(-1)-stencil(0))/2)

                else                            ! h(i-1) is valid
                                    ! (o.w. flux would most likely be out of cell)
                                    !  but h(i-2) is not

                  flux_diff = flux_diff + ABS(u_face) * G%dyCu(I-1,j) * time_step / G%areaT(i,j) * stencil(-1)

                endif

              elseif (u_face < 0) then !flux is out of cell - we need info from h(i-1), h(i+1) if available
                if (hmask(i-1,j) * hmask(i+1,j) == 1) then         ! h(i-1) and h(i+1) are both valid
                  phi = slope_limiter(stencil(0)-stencil(1), stencil(-1)-stencil(0))
                  flux_diff = flux_diff - ABS(u_face) * G%dyCu(I-1,j) * time_step / G%areaT(i,j) * &
                             (stencil(0) - phi * (stencil(0)-stencil(-1))/2)

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
                  flux_diff = flux_diff + ABS(u_face) * G%dyCu(I,j) * time_step / G%areaT(i,j) * &
                      (stencil(1) - phi * (stencil(1)-stencil(0))/2)

                else                            ! h(i+1) is valid
                                            ! (o.w. flux would most likely be out of cell)
                                            !  but h(i+2) is not

                  flux_diff = flux_diff + ABS(u_face) * G%dyCu(I,j) * time_step / G%areaT(i,j) * stencil(1)

                endif

              elseif (u_face > 0) then !flux is out of cell - we need info from h(i-1), h(i+1) if available

                if (hmask(i-1,j) * hmask(i+1,j) == 1) then         ! h(i-1) and h(i+1) are both valid

                  phi = slope_limiter(stencil(0)-stencil(-1), stencil(1)-stencil(0))
                  flux_diff = flux_diff - ABS(u_face) * G%dyCu(I,j) * time_step / G%areaT(i,j) * &
                      (stencil(0) - phi * (stencil(0)-stencil(1))/2)

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
                          intent(in)    :: h_after_uflux !< The ice shelf thicknesses after
                                              !! the zonal mass fluxes [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_vflux !< The ice shelf thicknesses after
                                              !! the meridional mass fluxes [Z ~> m].

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization

  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, gjed, gied
  integer :: i_off, j_off
  logical :: at_north_bdry, at_south_bdry, one_off_west_bdry, one_off_east_bdry
  real, dimension(-2:2) :: stencil
  real :: v_face     ! Pseudo-meridional velocity at a cell face, positive if out {L T-1 ~> m s-1]
  real :: flux_diff, phi
  character(len=1)        :: debug_str

  is = G%isc ; ie = G%iec ; js = G%jsc-1 ; je = G%jec+1 ; isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset

  do i=isd+2,ied-2
    if (((i+i_off) <= G%domain%niglobal+G%domain%nihalo) .AND. &
       ((i+i_off) >= G%domain%nihalo+1)) then  ! based on mehmet's code - only if btw east & west boundaries

      stencil(:) = -1

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
                  flux_diff = flux_diff + ABS(v_face) * G%dxCv(i,J-1) * time_step / G%areaT(i,j) * &
                      (stencil(-1) - phi * (stencil(-1)-stencil(0))/2)

                else     ! h(j-1) is valid
                         ! (o.w. flux would most likely be out of cell)
                         !  but h(j-2) is not
                  flux_diff = flux_diff + ABS(v_face) * G%dxCv(i,J-1) * time_step / G%areaT(i,j) * stencil(-1)
                endif

              elseif (v_face < 0) then !flux is out of cell - we need info from h(j-1), h(j+1) if available

                if (hmask(i,j-1) * hmask(i,j+1) == 1) then  ! h(j-1) and h(j+1) are both valid
                  phi = slope_limiter(stencil(0)-stencil(1), stencil(-1)-stencil(0))
                  flux_diff = flux_diff - ABS(v_face) * G%dxCv(i,J-1) * time_step / G%areaT(i,j) * &
                      (stencil(0) - phi * (stencil(0)-stencil(-1))/2)
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
                  flux_diff = flux_diff + ABS(v_face) * G%dxCv(i,J) * time_step / G%areaT(i,j) * &
                      (stencil(1) - phi * (stencil(1)-stencil(0))/2)
                else     ! h(j+1) is valid
                         ! (o.w. flux would most likely be out of cell)
                         !  but h(j+2) is not
                  flux_diff = flux_diff + ABS(v_face) * G%dxCv(i,J) * time_step / G%areaT(i,j) * stencil(1)
                endif

              elseif (v_face > 0) then !flux is out of cell - we need info from h(j-1), h(j+1) if available

                if (hmask(i,j-1) * hmask(i,j+1) == 1) then         ! h(j-1) and h(j+1) are both valid
                  phi = slope_limiter (stencil(0)-stencil(-1), stencil(1)-stencil(0))
                  flux_diff = flux_diff - ABS(v_face) * G%dxCv(i,J) * time_step / G%areaT(i,j) * &
                      (stencil(0) - phi * (stencil(0)-stencil(1))/2)
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
