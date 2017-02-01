module MOM_ice_shelf
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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
! Derived from code by Chris Little, early 2010.                               !
!                                                                              !
! This file implements the thermodynamic aspects of ocean / ice-shelf inter-   !
! actions, along with a crude placeholder for a later implementation of full   !
! ice shelf dynamics, all using the MOM framework and coding style.            !
!
! NOTE: THERE ARE A NUMBER OF SUBROUTINES WITH "TRIANGLE" IN THE NAME; THESE
! HAVE NOT BEEN TESTED AND SHOULD PROBABLY BE PHASED OUT

!   The ice-sheet dynamics subroutines do the following:
!  initialize_shelf_mass - Initializes the ice shelf mass distribution.
!      - Initializes h_shelf, h_mask, area_shelf_h
!      - CURRENTLY: initializes mass_shelf as well, but this is unnecessary, as mass_shelf is initialized based on
!             h_shelf and density_ice immediately afterwards. Possibly subroutine should be renamed
!  update_shelf_mass - Does nothing for now, but in some versions calls
!                      USER_update_shelf_mass.
!  ice_shelf_solve_outer - Orchestrates the calls to calculate the shelf
!      - outer loop calls ice_shelf_solve_inner
!                          stresses and checks for error tolerances.
!         Max iteration count for outer loop currently fixed at 100 iteration
!      - tolerance (and error evaluation) can be set through input file
!      - updates u_shelf, v_shelf, ice_visc_bilinear, taub_beta_eff_bilinear
!  ice_shelf_solve_inner - Conjugate Gradient solve of matrix solve for ice_shelf_solve_outer
!      - Jacobi Preconditioner - basically diagonal of matrix (not sure if it is effective at all)
!      - modifies u_shelf and v_shelf only
!      - max iteration count can be set through input file
!      - tolerance (and error evaluation) can be set through input file
!                  (ISSUE:  Too many mpp_sum calls?)
!    calc_shelf_driving_stress - Determine the driving stresses using h_shelf, (water) column thickness, bathymetry
!            - does not modify any permanent arrays
!    init_boundary_values -
!    bilinear_shape_functions - shape function for FEM solve using (convex) quadrilateral elements and bilinear nodal basis
!    calc_shelf_visc_bilinear - Glen's law viscosity and nonlinear sliding law (called by ice_shelf_solve_outer)
!    calc_shelf_visc_triangular - LET'S TAKE THIS OUT
!    apply_boundary_values_bilinear - same as CG_action_bilinear, but input is zero except for dirichlet bdry conds
!    apply_boundary_values_triangle - LET'S TAKE THIS OUT
!    CG_action_bilinear - Effect of matrix (that is never explicitly constructed)
!        on vector space of Degrees of Freedom (DoFs) in velocity solve
!    CG_action_triangular -LET'S TAKE THIS OUT
!      matrix_diagonal_bilinear - Returns the diagonal entries of a matrix for preconditioning.
!                  (ISSUE:  No need to use control structure - add arguments.
!      matrix_diagonal_triangle - LET'S TAKE THIS OUT
!  ice_shelf_advect - Given the melt rate and velocities, it advects the ice shelf THICKNESS
!        - modified h_shelf, area_shelf_h, hmask
!        (maybe should updater mass_shelf as well ???)
!    ice_shelf_advect_thickness_x, ice_shelf_advect_thickness_y - These
!        subroutines determine the mass fluxes through the faces.
!                  (ISSUE: duplicative flux calls for shared faces?)
!    ice_shelf_advance_front - Iteratively determine the ice-shelf front location.
!           - IF ice_shelf_advect_thickness_x,y are modified to avoid
!       dupe face processing, THIS NEEDS TO BE MODIFIED TOO
!       as it depends on arrays modified in those functions
!       (if in doubt consult DNG)
!    update_velocity_masks - Controls which elements of u_shelf and v_shelf are considered DoFs in linear solve
!    solo_time_step - called only in ice-only mode.
!    shelf_calc_flux - after melt rate & fluxes are calculated, ice dynamics are done. currently mass_shelf is
! updated immediately after ice_shelf_advect.
!
!
!   NOTES: be aware that hmask(:,:) has a number of functions; it is used for front advancement,
! for subroutines in the velocity solve, and for thickness boundary conditions (this last one may be removed).
! in other words, interfering with its updates will have implications you might not expect.
!
!  Overall issues: Many variables need better documentation and units and the
!                  subgrid on which they are discretized.

! DNG 4/09/11 : due to a misunderstanding (i confused a SYMMETRIC GRID
!      a SOUTHWEST GRID there is a variable called "isym" that appears
!      throughout in array loops. i am leaving it in for now,
!      though uniformly setting it to zero
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_mediator_init, set_diag_mediator_grid
use MOM_diag_mediator, only : diag_ctrl, time_type, enable_averaging, disable_averaging
use MOM_domains, only : MOM_domains_init, clone_MOM_domain
use MOM_domains, only : pass_var, pass_vector, TO_ALL, CGRID_NE, BGRID_NE
use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_grid, only : MOM_grid_init, ocean_grid_type
use MOM_grid_initialize, only : set_grid_metrics
use MOM_fixed_initialization, only : MOM_initialize_topography
use MOM_fixed_initialization, only : MOM_initialize_rotation
use user_initialization, only : user_initialize_topography
use MOM_io, only : field_exists, file_exists, read_data, write_version_number
use MOM_io, only : slasher, vardesc, var_desc, fieldtype
use MOM_io, only : write_field, close_file, SINGLE_FILE, MULTIPLE
use MOM_restart, only : register_restart_field, query_initialized, save_restart
use MOM_restart, only : restart_init, restore_state, MOM_restart_CS
use MOM_time_manager, only : time_type, set_time, time_type_to_real
use MOM_transcribe_grid, only : copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid
!use MOM_variables, only : forcing, surface
use MOM_variables, only : surface
use MOM_forcing_type, only : forcing, allocate_forcing_type
use MOM_get_input, only : directories, Get_MOM_input
use MOM_EOS, only : calculate_density, calculate_density_derivs, calculate_TFreeze
use MOM_EOS, only : EOS_type, EOS_init
!MJHuse MOM_ice_shelf_initialize, only : initialize_ice_shelf_boundary, initialize_ice_thickness
use MOM_ice_shelf_initialize, only : initialize_ice_thickness
use user_shelf_init, only : USER_initialize_shelf_mass, USER_update_shelf_mass
use user_shelf_init, only : user_ice_shelf_CS
use constants_mod,      only: GRAV
use mpp_mod, only : mpp_sum, mpp_max, mpp_min, mpp_pe, mpp_npes, mpp_sync
use MOM_coms, only : reproducing_sum
use MOM_debugging, only : hchksum, qchksum, chksum, uchksum, vchksum

implicit none ; private

#include <MOM_memory.h>
#ifdef SYMMETRIC_LAND_ICE
#  define GRID_SYM_ .true.
#  define NILIMB_SYM_ NIMEMB_SYM_
#  define NJLIMB_SYM_ NJMEMB_SYM_
#  define ISUMSTART_INT_ CS%grid%iscB+1
#  define JSUMSTART_INT_ CS%grid%jscB+1
#else
#  define GRID_SYM_ .false.
#  define NILIMB_SYM_ NIMEMB_
#  define NJLIMB_SYM_ NJMEMB_
#  define ISUMSTART_INT_ CS%grid%iscB
#  define JSUMSTART_INT_ CS%grid%jscB
#endif

public shelf_calc_flux, add_shelf_flux, initialize_ice_shelf, ice_shelf_end
public ice_shelf_save_restart, solo_time_step

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
type, public :: ice_shelf_CS ; private
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()
  type(ocean_grid_type) :: grid !< Grid for the ice-shelf model
!  type(dyn_horgrid_type), pointer :: dG  !< Dynamic grid for the ice-shelf model
  type(ocean_grid_type), pointer :: ocn_grid => NULL() !< A pointer to the ocean model grid
  ! The rest is private
  real ::   flux_factor = 1.0
  character(len=128) :: restart_output_dir = ' '
  real, pointer, dimension(:,:) :: &
    mass_shelf => NULL(), &   ! The mass per unit area of the ice shelf or sheet, in kg m-2.
    area_shelf_h => NULL(), & ! The area per cell covered by the ice shelf, in m2.

    t_flux => NULL(), &   ! The UPWARD sensible ocean heat flux at the ocean-ice
                    ! interface, in W m-2.
    salt_flux => NULL(), & ! The downward salt flux at the ocean-ice interface, in kg m-2 s-1.
    lprec => NULL(), &    ! The downward liquid water flux at the ocean-ice interface,
                    ! in kg m-2 s-1.
    ! Perhaps these diagnostics should only be kept with the call?
    exch_vel_t => NULL(), &
    exch_vel_s => NULL(), &
    utide   => NULL(), &
    tfreeze => NULL(), &  ! The freezing point potential temperature an the ice-ocean
                    ! interface, in deg C.
    tflux_shelf => NULL(), & ! The UPWARD diffusive heat flux in the ice shelf at the
                    ! ice-ocean interface, in W m-2.
!!! DNG !!!
    u_shelf => NULL(), & ! the zonal (?) velocity of the ice shelf/sheet... in meters per second???
          ! on q-points (B grid)
    v_shelf => NULL(), & ! the meridional velocity of the ice shelf/sheet... m/s ??
          ! on q-points (B grid)
    h_shelf => NULL(), & ! the thickness of the shelf in m... redundant with mass
          ! but may make code more readable
    hmask => NULL(),&    ! used to indicate ice-covered cells, as well as partially-covered
          ! 1: fully covered, solve for velocity here
          !   (for now all ice-covered cells are treated the same, this may change)
          ! 2: partially covered, do not solve for velocity
          ! 0: no ice in cell.
          ! 3: bdry condition on thickness set - not in computational domain
          ! -2 : default (out of computational boundary, and not = 3

          ! NOTE: hmask will change over time and NEEDS TO BE MAINTAINED
          !    otherwise the wrong nodes will be included in velocity calcs.
    u_face_mask => NULL(), v_face_mask => NULL(), &
   ! masks for velocity boundary conditions - on *C GRID* - this is because the FEM solution
   !      cares about FACES THAT GET INTEGRATED OVER, not vertices
   ! Will represent boundary conditions on computational boundary (or permanent boundary
   !      between fast-moving and near-stagnant ice
   ! FOR NOW: 1=interior bdry, 0=no-flow boundary, 2=stress bdry condition, 3=inhomogeneous dirichlet boundary,
   !          4=flux boundary: at these faces a flux will be specified which will override velocities;
   !                           a homogeneous velocity condition will be specified (this seems to give the solver less difficulty)
    u_face_mask_boundary => NULL(), v_face_mask_boundary => NULL(), &
    u_flux_boundary_values => NULL(), v_flux_boundary_values => NULL(), &
   ! needed where u_face_mask is equal to 4, similary for v_face_mask
    umask => NULL(), vmask => NULL(), &
   ! masks on the actual degrees of freedom (B grid) -
   !   1=normal node, 3=inhomogeneous boundary node, 0 - no flow node (will also get ice-free nodes)
    calve_mask => NULL(), & ! a mask to prevent the ice shelf front from advancing past its initial position (but it may retreat)

!!! OVS !!!
    t_shelf => NULL(), & ! veritcally integrated temperature the ice shelf/stream... oC
          ! on q-points (B grid)
     tmask => NULL(), &
  ! masks for temperature boundary conditions ???
    ice_visc_bilinear => NULL(), &
    ice_visc_lower_tri => NULL(), &
    ice_visc_upper_tri => NULL(), &
    thickness_boundary_values => NULL(), &
    u_boundary_values => NULL(), &
    v_boundary_values => NULL(), &
    h_boundary_values => NULL(), &
!!! OVS !!!
    t_boundary_values => NULL(), &

    taub_beta_eff_bilinear => NULL(), & ! nonlinear part of "linearized" basal stress - exact form depends on basal law exponent
                ! and/or whether flow is "hybridized" a la Goldberg 2011
    taub_beta_eff_lower_tri => NULL(), &
    taub_beta_eff_upper_tri => NULL(), &

    OD_rt => NULL(), float_frac_rt => NULL(), &
    OD_av => NULL(), float_frac => NULL()  !! two arrays that represent averages of ocean values that are maintained
                  !! within the ice shelf module and updated based on the "ocean state".
                  !! OD_av is ocean depth, and float_frac is the average amount of time
                  !! a cell is "exposed", i.e. the column thickness is below a threshold.
                  !! both are averaged over the time of a diagnostic (ice velocity)

                       !! [if float_frac = 1 ==> grounded; obv. counterintuitive; might fix]

  real :: ustar_bg     ! A minimum value for ustar under ice shelves, in m s-1.
  real :: cdrag        ! drag coefficient under ice shelves , non-dimensional.
  real :: g_Earth      ! The gravitational acceleration in m s-2.
  real :: Cp           ! The heat capacity of sea water, in J kg-1 K-1.
  real :: Rho0         ! A reference ocean density in kg/m3.
  real :: Cp_ice       ! The heat capacity of fresh ice, in J kg-1 K-1.
  real :: gamma_t      !   The (fixed) turbulent exchange velocity in the
                       ! 2-equation formulation, in m s-1.
  real :: Salin_ice    ! The salinity of shelf ice, in PSU.
  real :: Temp_ice     ! The core temperature of shelf ice, in C.
  real :: kv_ice       ! The viscosity of ice, in m2 s-1.
  real :: density_ice  ! A typical density of ice, in kg m-3.
  real :: kv_molec     ! The molecular kinematic viscosity of sea water, m2 s-1.
  real :: kd_molec_salt  ! The molecular diffusivity of salt, in m2 s-1.
  real :: kd_molec_temp  ! The molecular diffusivity of heat, in m2 s-1.
  real :: Lat_fusion   ! The latent heat of fusion, in J kg-1.
  real :: Gamma_T_3EQ  !  Nondimensional heat-transfer coefficient, used in the 3Eq. formulation
                       !  This number should be specified by the user.
  real :: col_thick_melt_threshold ! if the mixed layer is below this threshold, melt rate
                                  ! is not calculated for the cell

!!!! PHYSICAL AND NUMERICAL PARAMETERS FOR ICE DYNAMICS !!!!!!

  real :: time_step    ! this is the shortest timestep that the ice shelf sees, and
            ! is equal to the forcing timestep (it is passed in when the shelf
            ! is initialized - so need to reorganize MOM driver.
            ! it will be the prognistic timestep ... maybe.

!!! all need to be initialized

  logical :: solo_ice_sheet ! whether the ice model is running without being coupled to the ocean
  logical :: GL_regularize ! whether to regularize the floatation condition at the grounding line
                           !   a la Goldberg Holland Schoof 2009
  integer :: n_sub_regularize
                           ! partition of cell over which to integrate for interpolated grounding line
                           !  the (rectangular) is divided into nxn equally-sized rectangles, over which
                           !  basal contribution is integrated (iterative quadrature)
  logical :: GL_couple     ! whether to let the floatation condition be determined by ocean column thickness
                           !   means update_OD_ffrac will be called
                           !   (note: GL_regularize and GL_couple should be exclusive)

  real :: A_glen_isothermal
  real :: n_glen
  real :: eps_glen_min
  real :: C_basal_friction
  real :: n_basal_friction
  real :: density_ocean_avg    ! this does not affect ocean circulation OR thermodynamics
                ! it is to estimate the gravitational driving force at the shelf front
                ! (until we think of a better way to do it- but any difference will be negligible)
  real :: thresh_float_col_depth ! the water column depth over which the shelf if considered to be floating
  logical :: moving_shelf_front
  logical :: calve_to_mask
  real :: min_thickness_simple_calve
  real :: input_flux
  real :: input_thickness

  real :: len_lat ! this really should be a Grid or Domain field


  real :: velocity_update_time_step ! the time to update the velocity through the nonlinear
                    ! elliptic equation. i think this should be done no more often than
                    ! ~ once a day (maybe longer) because it will depend on ocean values
                    ! that are averaged over this time interval, and the solve will begin
                    ! to lose meaning if it is done too frequently
  integer :: velocity_update_sub_counter ! there is no outer loop for the velocity solve; the counter will have to be stored
  integer :: velocity_update_counter ! the "outer" timestep number
  integer :: nstep_velocity        ! ~ (velocity_update_time_step / time_step)

  real :: cg_tolerance, nonlinear_tolerance
  integer :: cg_max_iterations
  integer :: nonlin_solve_err_mode  ! 1: exit vel solve based on nonlin residual
                    ! 2: exit based on "fixed point" metric (|u - u_last| / |u| < tol where | | is infty-norm
  real    :: CFL_factor            ! in uncoupled run, how to limit subcycled advective timestep
                      ! i.e. dt = CFL_factor * min (dx / u)
  logical :: use_reproducing_sums ! use new reproducing sums of Bob & Alistair for global sums
                                  ! NOTE: for this to work all tiles must have the same & of
                                  !       elements. this means thatif a symmetric grid is being
                                  !       used, the southwest nodes of the southwest tiles will not
                                  !       be included in the


  logical :: switch_var ! for debdugging - a switch to ensure some event happens only once

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(time_type) :: Time ! The component's time.
  type(EOS_type), pointer :: eqn_of_state => NULL() ! Type that indicates the
                                        ! equation of state to use.
  logical :: shelf_mass_is_dynamic ! True if the ice shelf mass changes with
                       ! time.
  logical :: override_shelf_movement ! If true, user code specifies the shelf
                       ! movement instead of using the dynamic ice-shelf mode.
  logical :: isthermo  ! True if the ice shelf can exchange heat and mass with
                       ! the underlying ocean.
  logical :: threeeq   ! If true, the 3 equation consistency equations are
                       ! used to calculate the flux at the ocean-ice interface.
  logical :: insulator ! If true, ice shelf is a perfect insulator
  logical :: const_gamma    ! If true, gamma_T is specified by the user.
  logical :: find_salt_root ! If true, if true find Sbdry using a quadratic eq.
  real    :: cutoff_depth ! depth above which melt is set to zero (>= 0).
  real    :: lambda1, lambda2, lambda3 ! liquidus coeffs. Needed if
                                       ! find_salt_root = true
  integer :: id_melt = -1, id_exch_vel_s = -1, id_exch_vel_t = -1, &
             id_tfreeze = -1, id_tfl_shelf = -1, &
             id_thermal_driving = -1, id_haline_driving = -1, &
             id_u_ml = -1, id_v_ml = -1, id_sbdry = -1, &
             id_u_shelf = -1, id_v_shelf = -1, id_h_shelf = -1, id_h_mask = -1, &
             id_u_mask = -1, id_v_mask = -1, id_t_shelf = -1, id_t_mask = -1, &
             id_surf_elev = -1, id_bathym = -1, id_float_frac = -1, id_col_thick = -1, &
             id_area_shelf_h = -1, id_OD_av = -1, id_float_frac_rt = -1,&
             id_ustar_shelf = -1, id_shelf_mass = -1, id_mass_flux = -1

  ! ids for outputting intermediate thickness in advection subroutine (debugging)
  !integer :: id_h_after_uflux = -1, id_h_after_vflux = -1, id_h_after_adv = -1

  type(diag_ctrl), pointer :: diag     ! A structure that is used to control diagnostic
                              ! output.
  type(user_ice_shelf_CS), pointer :: user_CS => NULL()

  logical :: write_output_to_file ! this is for seeing arrays w/out netcdf capability
  logical :: debug                ! If true, write verbose checksums for debugging purposes
                                  ! and use reproducible sums
end type ice_shelf_CS

integer :: id_clock_shelf, id_clock_pass

contains

! ------DNG-------------------------------

! FUNCTION CALLS

function slope_limiter (num, denom)
  real, intent(in)     :: num
  real, intent(in)    :: denom
  real :: slope_limiter
  real :: r

! used for flux limiting in advective subroutines
! Van Leer limiter (source: Wikipedia!)

  if (denom .eq. 0) then
    slope_limiter = 0
  elseif (num*denom .le. 0) then
    slope_limiter = 0
  else
    r = num/denom
    slope_limiter = (r+abs(r))/(1+abs(r))
  endif

end function slope_limiter

function quad_area (X, Y)
  real, dimension(4), intent(in) :: X
  real, dimension(4), intent(in) :: Y
  real :: quad_area, p2, q2, a2, c2, b2, d2

! area of quadrilateral
! X and Y must be passed in the form
    !  3 - 4
    !  |   |
    !  1 - 2

  p2 = (X(4)-X(1))**2 + (Y(4)-Y(1))**2 ; q2 = (X(3)-X(2))**2 + (Y(3)-Y(2))**2
  a2 = (X(3)-X(4))**2 + (Y(3)-Y(4))**2 ; c2 = (X(1)-X(2))**2 + (Y(1)-Y(2))**2
  b2 = (X(2)-X(4))**2 + (Y(2)-Y(4))**2 ; d2 = (X(3)-X(1))**2 + (Y(3)-Y(1))**2
  quad_area = .25 * sqrt(4*P2*Q2-(B2+D2-A2-C2)**2)

end function quad_area

!-----------------------

subroutine shelf_calc_flux(state, fluxes, Time, time_step, CS)
  type(surface),         intent(inout) :: state
  type(forcing),         intent(inout) :: fluxes
  type(time_type),       intent(in)    :: Time
  real,                  intent(in)    :: time_step
  type(ice_shelf_CS),    pointer       :: CS

! This sets up the fluxes between the ocean and an ice-shelf.
!
! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (inout)   fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      Time - Start time of the fluxes.
!  (in)      time_step - Length of time over which these fluxes
!                        will be applied, in s.
!  (in)      CS - A pointer to the control structure returned by a previous
!                 call to initialize_ice_shelf.

  real, dimension(SZI_(CS%grid)) :: &
    Rhoml, &   ! Ocean mixed layer density in kg m-3.
    dR0_dT, &  ! Partial derivative of the mixed layer density
               ! with temperature, in units of kg m-3 K-1.
    dR0_dS, &  ! Partial derivative of the mixed layer density
               ! with salinity, in units of kg m-3 psu-1.
    p_int   ! The pressure at the ice-ocean interface, in Pa.

  real, dimension(:,:), allocatable :: mass_flux  ! total mass flux of freshwater across
  real, dimension(:,:), allocatable :: haline_driving  ! (SSS - S_boundary) ice-ocean interface, positive for melting and negative for freezing. This is computed as part of the ISOMIP diagnostics.
  real, parameter :: VK    = 0.40     ! Von Karman's constant - dimensionless
  real :: ZETA_N = 0.052   ! The fraction of the boundary layer over which the
                           ! viscosity is linearly increasing. (Was 1/8. Why?)
  real, parameter :: RC    = 0.20     ! critical flux Richardson number.
  real :: I_ZETA_N  ! The inverse of ZETA_N.
  real :: LF, I_LF  ! Latent Heat of fusion (J kg-1) and its inverse.
  real :: I_VK      ! The inverse of VK.
  real :: PR, SC    ! The Prandtl number and Schmidt number, nondim.

  ! 3 equation formulation variables
  real, dimension(:,:), allocatable :: Sbdry !Salinities in the ocean at the interface with the ice shelf, in PSU.
  real :: Sbdry_it
  real :: Sbdry1, Sbdry2, S_a, S_b, S_c  ! use to find salt roots
  real :: dS_it     ! The interface salinity change during an iteration, in PSU.
  real :: hBL_neut  ! The neutral boundary layer thickness, in m.
  real :: hBL_neut_h_molec ! The ratio of the neutral boundary layer thickness
                           ! to the molecular boundary layer thickness, ND.
  real :: wT_flux ! The vertical fluxes of heat and buoyancy just inside the
  real :: wB_flux ! ocean, in C m s-1 and m2 s-3, ###CURRENTLY POSITIVE UPWARD.
  real :: dB_dS  ! The derivative of buoyancy with salinity, in m s-2 PSU-1.
  real :: dB_dT  ! The derivative of buoyancy with temperature, in m s-2 C-1.
  real :: I_n_star, n_star_term, absf
  real :: dIns_dwB  ! The partial derivative of I_n_star with wB_flux, in ???.
  real :: dT_ustar, dS_ustar
  real :: ustar_h
  real :: Gam_turb
  real :: Gam_mol_t, Gam_mol_s
  real :: RhoCp
  real :: I_RhoLF
  real :: ln_neut
  real :: mass_exch
  real :: Sb_min, Sb_max
  real :: dS_min, dS_max
  ! Variables used in iterating for wB_flux.
  real :: wB_flux_new, DwB, dDwB_dwB_in
  real :: I_Gam_T, I_Gam_S, dG_dwB, iDens
  real :: u_at_h, v_at_h, Isqrt2
  logical :: Sb_min_set, Sb_max_set
  character(4) :: stepnum
  character(2) :: procnum

  type(ocean_grid_type), pointer :: G
  real, parameter :: c2_3 = 2.0/3.0
  integer :: i, j, is, ie, js, je, ied, jed, it1, it3, iters_vel_solve
  real, parameter :: rho_fw = 1000.0 ! fresh water density
  if (.not. associated(CS)) call MOM_error(FATAL, "shelf_calc_flux: "// &
       "initialize_ice_shelf must be called before shelf_calc_flux.")
  call cpu_clock_begin(id_clock_shelf)

  G => CS%grid
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; ied = G%ied ; jed = G%jed
  I_ZETA_N = 1.0 / ZETA_N
  LF = CS%Lat_fusion
  I_RhoLF = 1.0/(CS%Rho0*LF)
  I_LF = 1.0 / LF
  SC = CS%kv_molec/CS%kd_molec_salt
  PR = CS%kv_molec/CS%kd_molec_temp
  I_VK = 1.0/VK
  RhoCp = CS%Rho0 * CS%Cp
  Isqrt2 = 1.0/sqrt(2.0)

!first calculate molecular component
  Gam_mol_t = 12.5 * (PR**c2_3) - 6
  Gam_mol_s = 12.5 * (SC**c2_3) - 6

  iDens = 1.0/CS%density_ocean_avg

  ! GM, zero some fields of the ice shelf structure (ice_shelf_CS)
  ! these fields are already set to zero during initialization
  ! However, they seem to be changed somewhere and, for diagnostic
  ! reasons, it is better to set them to zero again.
  CS%tflux_shelf(:,:) = 0.0; CS%exch_vel_t(:,:) = 0.0
  CS%lprec(:,:) = 0.0; CS%exch_vel_s(:,:) = 0.0
  CS%salt_flux(:,:) = 0.0; CS%t_flux(:,:) = 0.0
  CS%tfreeze(:,:) = 0.0
  ! define Sbdry to avoid Run-Time Check Failure, when melt is not computed.
  ALLOCATE ( haline_driving(G%ied,G%jed) ); haline_driving(:,:) = 0.0
  ALLOCATE ( Sbdry(G%ied,G%jed) ); Sbdry(:,:) = state%sss(:,:)

  if (CS%shelf_mass_is_dynamic .and. CS%override_shelf_movement) &
                                  call update_shelf_mass(CS, Time)

  do j=js,je
    ! Find the pressure at the ice-ocean interface, averaged only over the
    ! part of the cell covered by ice shelf.
    do i=is,ie ; p_int(i) = CS%g_Earth * CS%mass_shelf(i,j) ; enddo

    ! Calculate insitu densities and expansion coefficients
    call calculate_density(state%sst(:,j),state%sss(:,j), p_int, &
             Rhoml(:), is, ie-is+1, CS%eqn_of_state)
    call calculate_density_derivs(state%sst(:,j), state%sss(:,j), p_int, &
             dR0_dT, dR0_dS, is, ie-is+1, CS%eqn_of_state)

    do i=is,ie

      ! DNG - to allow this everywhere Hml>0.0 allows for melting under grounded cells
      !       propose instead to allow where Hml > [some threshold]

      if ((iDens*state%ocean_mass(i,j) > CS%col_thick_melt_threshold) .and. &
          (CS%area_shelf_h(i,j) > 0.0) .and. &
          (CS%isthermo) .and. (state%Hml(i,j) > 0.0) ) then

        if (CS%threeeq) then
          !   Iteratively determine a self-consistent set of fluxes, with the ocean
          ! salinity just below the ice-shelf as the variable that is being
          ! iterated for.
          ! ### SHOULD I SET USTAR_SHELF YET?

          u_at_h = state%u(i,j)
          v_at_h = state%v(i,j)

          fluxes%ustar_shelf(i,j)= sqrt(CS%cdrag*((u_at_h**2.0 + v_at_h**2.0) +&
                                                    CS%utide(i,j)**2))

          ustar_h = MAX(CS%ustar_bg, fluxes%ustar_shelf(i,j))

          fluxes%ustar_shelf(i,j) = ustar_h

          if (associated(state%taux_shelf) .and. associated(state%tauy_shelf)) then
            state%taux_shelf(i,j) = ustar_h*ustar_h*CS%Rho0*Isqrt2
            state%tauy_shelf(i,j) = state%taux_shelf(i,j)
          endif

          ! Estimate the neutral ocean boundary layer thickness as the minimum of the
          ! reported ocean mixed layer thickness and the neutral Ekman depth.
          absf = 0.25*((abs(G%CoriolisBu(I,J)) + abs(G%CoriolisBu(I-1,J-1))) + &
                       (abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I-1,J))))
          if (absf*state%Hml(i,j) <= VK*ustar_h) then ; hBL_neut = state%Hml(i,j)
          else ; hBL_neut = (VK*ustar_h) / absf ; endif
          hBL_neut_h_molec = ZETA_N * ((hBL_neut * ustar_h) / (5.0 * CS%Kv_molec))

          ! Determine the mixed layer buoyancy flux, wB_flux.
          dB_dS = (CS%g_Earth / Rhoml(i)) * dR0_dS(i)
          dB_dT = (CS%g_Earth / Rhoml(i)) * dR0_dT(i)
          ln_neut = 0.0 ; if (hBL_neut_h_molec > 1.0) ln_neut = log(hBL_neut_h_molec)

          if (CS%find_salt_root) then
            ! read liquidus parameters

            S_a = CS%lambda1 * CS%Gamma_T_3EQ * CS%Cp
!            S_b = -CS%Gamma_T_3EQ*(CS%lambda2-CS%lambda3*p_int(i)-state%sst(i,j)) &
!               -LF*CS%Gamma_T_3EQ/35.0

            S_b = CS%Gamma_T_3EQ*CS%Cp*(CS%lambda2+CS%lambda3*p_int(i)- &
                  state%sst(i,j))-LF*CS%Gamma_T_3EQ/35.0
            S_c = LF*(CS%Gamma_T_3EQ/35.0)*state%sss(i,j)

            Sbdry1 = (-S_b + SQRT(S_b*S_b-4*S_a*S_c))/(2*S_a)
            Sbdry2 = (-S_b - SQRT(S_b*S_b-4*S_a*S_c))/(2*S_a)
            Sbdry(i,j) = MAX(Sbdry1, Sbdry2)
            ! Safety check
            if (Sbdry(i,j) < 0.) then
               write(*,*)'state%sss(i,j)',state%sss(i,j)
               write(*,*)'S_a, S_b, S_c',S_a, S_b, S_c
               write(*,*)'I,J,Sbdry1,Sbdry2',i,j,Sbdry1,Sbdry2
               call MOM_error(FATAL, &
                  "shelf_calc_flux: Negative salinity (Sbdry).")
            endif
          else
            ! Guess sss as the iteration starting point for the boundary salinity.
            Sbdry(i,j) = state%sss(i,j) ; Sb_max_set = .false.
            Sb_min_set = .false.
          endif !find_salt_root

          do it1 = 1,20
            ! Determine the potential temperature at the ice-ocean interface.
            call calculate_TFreeze(Sbdry(i,j), p_int(i), CS%tfreeze(i,j), CS%eqn_of_state)

            dT_ustar = (state%sst(i,j) - CS%tfreeze(i,j)) * ustar_h
            dS_ustar = (state%sss(i,j) - Sbdry(i,j)) * ustar_h

            ! First, determine the buoyancy flux assuming no effects of stability
            ! on the turbulence.  Following H & J '99, this limit also applies
            ! when the buoyancy flux is destabilizing.

            if (CS%const_gamma) then ! if using a constant gamma_T
               ! note the different form, here I_Gam_T is NOT 1/Gam_T!
               I_Gam_T = CS%Gamma_T_3EQ
               I_Gam_S = CS%Gamma_T_3EQ/35.
            else
               Gam_turb = I_VK * (ln_neut + (0.5 * I_ZETA_N - 1.0))
               I_Gam_T = 1.0 / (Gam_mol_t + Gam_turb)
               I_Gam_S = 1.0 / (Gam_mol_s + Gam_turb)
            endif

            wT_flux = dT_ustar * I_Gam_T
            wB_flux = dB_dS * (dS_ustar * I_Gam_S) + dB_dT * wT_flux

            if (wB_flux > 0.0) then
              ! The buoyancy flux is stabilizing and will reduce the tubulent
              ! fluxes, and iteration is required.
              n_star_term = (ZETA_N/RC) * (hBL_neut * VK) / ustar_h**3
              do it3 = 1,30
               ! n_star <= 1.0 is the ratio of working boundary layer thickness
               ! to the neutral thickness.
               ! hBL = n_star*hBL_neut ; hSub = 1/8*n_star*hBL

                I_n_star = sqrt(1.0 + n_star_term * wB_flux)
                dIns_dwB = 0.5 * n_star_term / I_n_star
                if (hBL_neut_h_molec > I_n_star**2) then
                  Gam_turb = I_VK * ((ln_neut - 2.0*log(I_n_star)) + &
                                    (0.5*I_ZETA_N*I_n_star - 1.0))
                  dG_dwB =  I_VK * ( -2.0 / I_n_star + (0.5 * I_ZETA_N)) * dIns_dwB
                else
                  !   The layer dominated by molecular viscosity is smaller than
                  ! the assumed boundary layer.  This should be rare!
                  Gam_turb = I_VK * (0.5 * I_ZETA_N*I_n_star - 1.0)
                  dG_dwB = I_VK * (0.5 * I_ZETA_N) * dIns_dwB
                endif

                if (CS%const_gamma) then ! if using a constant gamma_T
                   ! note the different form, here I_Gam_T is NOT 1/Gam_T!
                   I_Gam_T = CS%Gamma_T_3EQ
                   I_Gam_S = CS%Gamma_T_3EQ/35.
                else
                  I_Gam_T = 1.0 / (Gam_mol_t + Gam_turb)
                  I_Gam_S = 1.0 / (Gam_mol_s + Gam_turb)
                endif

                wT_flux = dT_ustar * I_Gam_T
                wB_flux_new = dB_dS * (dS_ustar * I_Gam_S) + dB_dT * wT_flux

                ! Find the root where dwB = 0.0
                DwB = wB_flux_new - wB_flux
                if (abs(wB_flux_new - wB_flux) < &
                    1e-4*(abs(wB_flux_new) + abs(wB_flux))) exit

                dDwB_dwB_in = -dG_dwB * (dB_dS * (dS_ustar * I_Gam_S**2) + &
                                         dB_dT * (dT_ustar * I_Gam_T**2)) - 1.0
                ! This is Newton's method without any bounds.
                ! ### SHOULD BOUNDS BE NEEDED?
                wB_flux_new = wB_flux - DwB / dDwB_dwB_in
              enddo !it3
            endif

            CS%t_flux(i,j)  = RhoCp * wT_flux
            CS%exch_vel_t(i,j) = ustar_h * I_Gam_T
            CS%exch_vel_s(i,j) = ustar_h * I_Gam_S

    !Calculate the heat flux inside the ice shelf.

    !vertical adv/diff as in H+J 1999, eqns (26) & approx from (31).
    ! Q_ice = rho_ice * CS%CP_Ice * K_ice * dT/dz (at interface)
    !vertical adv/diff as in H+J 199, eqs (31) & (26)...
    !  dT/dz ~= min( (lprec/(rho_ice*K_ice))*(CS%Temp_Ice-T_freeze) , 0.0 )
    !If this approximation is not made, iterations are required... See H+J Fig 3.

            if (CS%t_flux(i,j) <= 0.0) then  ! Freezing occurs, so zero ice heat flux.
              CS%lprec(i,j) = I_LF * CS%t_flux(i,j)
              CS%tflux_shelf(i,j) = 0.0
            else
              if (CS%insulator) then
                 !no conduction/perfect insulator
                 CS%tflux_shelf(i,j) = 0.0
                 CS%lprec(i,j) = I_LF * (- CS%tflux_shelf(i,j) + CS%t_flux(i,j))

              else
                 ! With melting, from H&J 1999, eqs (31) & (26)...
                 !   Q_ice ~= cp_ice * (CS%Temp_Ice-T_freeze) * lprec
                 !   RhoLF*lprec = Q_ice + CS%t_flux(i,j)
                 !   lprec = (CS%t_flux(i,j)) / (LF + cp_ice * (T_freeze-CS%Temp_Ice))
                 CS%lprec(i,j) = CS%t_flux(i,j) / &
                       (LF + CS%CP_Ice * (CS%Tfreeze(i,j) - CS%Temp_Ice))

                CS%tflux_shelf(i,j) = CS%t_flux(i,j) - LF*CS%lprec(i,j)
              endif

            endif
            !other options: dTi/dz linear through shelf
            !    dTi_dz = (CS%Temp_Ice - CS%tfreeze(i,j))/G%draft(i,j)
            !    CS%tflux_shelf(i,j) = - Rho_Ice * CS%CP_Ice * KTI * dTi_dz


            if (CS%find_salt_root) then
              exit ! no need to do interaction, so exit loop
            else

              mass_exch = CS%exch_vel_s(i,j) * CS%Rho0
              Sbdry_it = (state%sss(i,j) * mass_exch + CS%Salin_ice * &
                          CS%lprec(i,j)) / (mass_exch + CS%lprec(i,j))
              dS_it = Sbdry_it - Sbdry(i,j)
              if (abs(dS_it) < 1e-4*(0.5*(state%sss(i,j) + Sbdry(i,j) + 1.e-10))) exit


              if (dS_it < 0.0) then ! Sbdry is now the upper bound.
                if (Sb_max_set .and. (Sbdry(i,j) > Sb_max)) &
                call MOM_error(FATAL,"shelf_calc_flux: Irregular iteration for Sbdry (max).")
                Sb_max = Sbdry(i,j) ; dS_max = dS_it ; Sb_max_set = .true.
              else ! Sbdry is now the lower bound.
                if (Sb_min_set .and. (Sbdry(i,j) < Sb_min)) &
                   call MOM_error(FATAL, &
                   "shelf_calc_flux: Irregular iteration for Sbdry (min).")
                   Sb_min = Sbdry(i,j) ; dS_min = dS_it ; Sb_min_set = .true.
              endif ! dS_it < 0.0

              if (Sb_min_set .and. Sb_max_set) then
                 ! Use the false position method for the next iteration.
                 Sbdry(i,j) = Sb_min + (Sb_max-Sb_min) * &
                             (dS_min / (dS_min - dS_max))
              else
                 Sbdry(i,j) = Sbdry_it
              endif ! Sb_min_set

              Sbdry(i,j) = Sbdry_it
            endif ! CS%find_salt_root

          enddo !it1
  ! Check for non-convergence and/or non-boundedness?

        else
          !   In the 2-equation form, the mixed layer turbulent exchange velocity
          ! is specified and large enough that the ocean salinity at the interface
          ! is about the same as the boundary layer salinity.

          call calculate_TFreeze(state%sss(i,j), p_int(i), CS%tfreeze(i,j), CS%eqn_of_state)

          CS%exch_vel_t(i,j) = CS%gamma_t
          CS%t_flux(i,j) = RhoCp * CS%exch_vel_t(i,j) * (state%sst(i,j) - CS%tfreeze(i,j))
          CS%tflux_shelf(i,j) = 0.0
          CS%lprec(i,j) = I_LF * CS%t_flux(i,j)
          Sbdry(i,j) = 0.0
        endif
      else !not shelf
        CS%t_flux(i,j) = 0.0
      endif

!      haline_driving(:,:) = state%sss(i,j) - Sbdry(i,j)

    enddo ! i-loop
  enddo ! j-loop

  ! CS%lprec = precipitating liquid water into the ocean ( kg/(m^2 s) )
  ! We want melt in m/year
  if (CS%const_gamma) then ! use ISOMIP+ eq. with rho_fw
    fluxes%iceshelf_melt = CS%lprec  * (86400.0*365.0/rho_fw)
  else ! use original eq.
    fluxes%iceshelf_melt = CS%lprec  * (86400.0*365.0/CS%density_ice)
  endif

  do j=js,je
    do i=is,ie
      if ((iDens*state%ocean_mass(i,j) > CS%col_thick_melt_threshold) .and. &
          (CS%area_shelf_h(i,j) > 0.0) .and. &
          (CS%isthermo) .and. (state%Hml(i,j) > 0.0) ) then

         ! Set melt to zero above a cutoff pressure
         ! (CS%Rho0*CS%cutoff_depth*CS%g_Earth) this is needed for the isomip
         ! test case.
         if ((CS%g_Earth * CS%mass_shelf(i,j)) < CS%Rho0*CS%cutoff_depth* &
            CS%g_Earth) then
              CS%lprec(i,j) = 0.0
              fluxes%iceshelf_melt(i,j) = 0.0
         endif
         ! Compute haline driving, which is one of the diags. used in ISOMIP
         haline_driving(i,j) = (CS%lprec(i,j) * Sbdry(i,j)) / &
                               (CS%Rho0 * CS%exch_vel_s(i,j))

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!Safety checks !!!!!!!!!!!!!!!!!!!!!!!!!
         !1)Check if haline_driving computed above is consistent with
         ! haline_driving = state%sss - Sbdry
         !if (fluxes%iceshelf_melt(i,j) /= 0.0) then
         !   if (haline_driving(i,j) /= (state%sss(i,j) - Sbdry(i,j))) then
         !      write(*,*)'Something is wrong at i,j',i,j
         !      write(*,*)'haline_driving, sss-Sbdry',haline_driving(i,j), &
         !                (state%sss(i,j) - Sbdry(i,j))
         !     call MOM_error(FATAL, &
         !            "shelf_calc_flux: Inconsistency in melt and haline_driving")
         !   endif
         !endif

         ! 2) check if |melt| > 0 when star_shelf = 0.
         ! this should never happen
         if (abs(fluxes%iceshelf_melt(i,j))>0.0) then
             if (fluxes%ustar_shelf(i,j) == 0.0) then
                write(*,*)'Something is wrong at i,j',i,j
                call MOM_error(FATAL, &
                     "shelf_calc_flux: |melt| > 0 and star_shelf = 0.")
             endif
          endif
      endif ! area_shelf_h
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!End of safety checks !!!!!!!!!!!!!!!!!!!
     enddo ! i-loop
   enddo ! j-loop

  ! mass flux (kg/s), part of ISOMIP diags.
  ALLOCATE ( mass_flux(G%ied,G%jed) ); mass_flux(:,:) = 0.0
  mass_flux = (CS%lprec) * CS%area_shelf_h

  if (CS%DEBUG) then
   call hchksum (fluxes%iceshelf_melt, "melt rate", G%HI, haloshift=0)
   call hchksum (fluxes%ustar_shelf, "ustar_shelf calc", G%HI, haloshift=0)
   call hchksum (state%Hml, "Hml", G%HI, haloshift=0)
  endif

  if (CS%shelf_mass_is_dynamic) then
    call cpu_clock_begin(id_clock_pass)
    call pass_var(CS%area_shelf_h, G%domain, complete=.false.)
    call pass_var(CS%mass_shelf, G%domain)
    call cpu_clock_end(id_clock_pass)
  endif

  call add_shelf_flux(G, CS, state, fluxes)

  ! now the thermodynamic data is passed on... time to update the ice dynamic quantities

  if (CS%shelf_mass_is_dynamic .and. .not.CS%override_shelf_movement) then

    ! advect the ice shelf, and advance the front. Calving will be in here somewhere as well..
    ! when we decide on how to do it

    ! note time_step is [s] and lprec is [kg / m^2 / s]

    call ice_shelf_advect (CS, time_step, CS%lprec, Time)


    do j=G%jsd,G%jed
      do i=G%isd,G%ied
        if ((CS%hmask(i,j) .eq. 1) .or. (CS%hmask(i,j) .eq. 2)) then
          CS%mass_shelf(i,j) = CS%h_shelf(i,j)*CS%density_ice
        endif
      enddo
    enddo


    CS%velocity_update_sub_counter = CS%velocity_update_sub_counter+1

    if (CS%GL_couple .and. .not. CS%solo_ice_sheet) then
      call update_OD_ffrac (CS, state%ocean_mass, CS%velocity_update_sub_counter, CS%nstep_velocity, CS%time_step, CS%velocity_update_time_step)
    else
      call update_OD_ffrac_uncoupled (CS)
    endif

    if (CS%velocity_update_sub_counter .eq. CS%nstep_velocity) then

      if (is_root_pe()) write(*,*) "ABOUT TO CALL VELOCITY SOLVER"

      call ice_shelf_solve_outer (CS, CS%u_shelf, CS%v_shelf, 1, iters_vel_solve, Time)

      CS%velocity_update_sub_counter = 0

    endif
  endif

  call enable_averaging(time_step,Time,CS%diag)
   if (CS%id_shelf_mass > 0) call post_data(CS%id_shelf_mass, CS%mass_shelf, CS%diag)
   if (CS%id_area_shelf_h > 0) call post_data(CS%id_area_shelf_h, CS%area_shelf_h, CS%diag)
   if (CS%id_ustar_shelf > 0) call post_data(CS%id_ustar_shelf, fluxes%ustar_shelf, CS%diag)
   if (CS%id_melt > 0) call post_data(CS%id_melt, fluxes%iceshelf_melt, CS%diag)
   if (CS%id_thermal_driving > 0) call post_data(CS%id_thermal_driving, (state%sst-CS%tfreeze), CS%diag)
   if (CS%id_Sbdry > 0) call post_data(CS%id_Sbdry, Sbdry, CS%diag)
   if (CS%id_haline_driving > 0) call post_data(CS%id_haline_driving, haline_driving, CS%diag)
   if (CS%id_mass_flux > 0) call post_data(CS%id_mass_flux, mass_flux, CS%diag)
   if (CS%id_u_ml > 0) call post_data(CS%id_u_ml,state%u,CS%diag)
   if (CS%id_v_ml > 0) call post_data(CS%id_v_ml,state%v,CS%diag)
   if (CS%id_tfreeze > 0) call post_data(CS%id_tfreeze, CS%tfreeze, CS%diag)
   if (CS%id_tfl_shelf > 0) call post_data(CS%id_tfl_shelf, CS%tflux_shelf, CS%diag)
   if (CS%id_exch_vel_t > 0) call post_data(CS%id_exch_vel_t, CS%exch_vel_t, CS%diag)
   if (CS%id_exch_vel_s > 0) call post_data(CS%id_exch_vel_s, CS%exch_vel_s, CS%diag)
   if (CS%id_col_thick > 0) call post_data(CS%id_col_thick, CS%OD_av, CS%diag)
   if (CS%id_h_shelf > 0) call post_data(CS%id_h_shelf,CS%h_shelf,CS%diag)
   if (CS%id_h_mask > 0) call post_data(CS%id_h_mask,CS%hmask,CS%diag)
   if (CS%id_u_shelf > 0) call post_data(CS%id_u_shelf,CS%u_shelf,CS%diag)
   if (CS%id_v_shelf > 0) call post_data(CS%id_v_shelf,CS%v_shelf,CS%diag)
   if (CS%id_float_frac > 0) call post_data(CS%id_float_frac,CS%float_frac,CS%diag)
   if (CS%id_OD_av >0) call post_data(CS%id_OD_av,CS%OD_av,CS%diag)
   if (CS%id_float_frac_rt>0) call post_data(CS%id_float_frac_rt,CS%float_frac_rt,CS%diag)
  call disable_averaging(CS%diag)

  call cpu_clock_end(id_clock_shelf)

end subroutine shelf_calc_flux

subroutine add_shelf_flux(G, CS, state, fluxes)
  type(ocean_grid_type),              intent(inout)    :: G
  type(ice_shelf_CS),                 intent(inout)    :: CS
  type(surface),                      intent(inout)    :: state
  type(forcing),                      intent(inout) :: fluxes
! Arguments:
!  (in)      fluxes - A structure of surface fluxes that may be used.
!  (in)      visc - A structure containing vertical viscosities, bottom boundary
!                   layer properies, and related fields.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - This module's control structure.
 !need to use visc variables
 !time step therm v. dynamic?
  real :: Irho0         ! The inverse of the mean density in m3 kg-1.
  real :: frac_area     ! The fractional area covered by the ice shelf, nondim.
  real :: taux2, tauy2  ! The squared surface stresses, in Pa.
  real :: asu1, asu2    ! Ocean areas covered by ice shelves at neighboring u-
  real :: asv1, asv2    ! and v-points, in m2.
  real :: fraz          ! refreezing rate in kg m-2 s-1
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed

  Irho0 = 1.0 / CS%Rho0
  ! Determine ustar and the square magnitude of the velocity in the
  ! bottom boundary layer. Together these give the TKE source and
  ! vertical decay scale.
  if (CS%shelf_mass_is_dynamic) then
    do j=jsd,jed ; do i=isd,ied
      if (G%areaT(i,j) > 0.0) &
        fluxes%frac_shelf_h(i,j) = CS%area_shelf_h(i,j) / G%areaT(i,j)
    enddo ; enddo
    !do I=isd,ied-1 ; do j=isd,jed
    do j=jsd,jed ; do i=isd,ied-1 ! ### changed stride order; i->ied-1?
      fluxes%frac_shelf_u(I,j) = 0.0
      if ((G%areaT(i,j) + G%areaT(i+1,j) > 0.0)) & ! .and. (G%dxdy_u(I,j) > 0.0)) &
        fluxes%frac_shelf_u(I,j) = ((CS%area_shelf_h(i,j) + CS%area_shelf_h(i+1,j)) / &
                                    (G%areaT(i,j) + G%areaT(i+1,j)))
      fluxes%rigidity_ice_u(I,j) = (CS%kv_ice / CS%density_ice) * &
                                    min(CS%mass_shelf(i,j), CS%mass_shelf(i+1,j))
    enddo ; enddo
    do j=jsd,jed-1 ; do i=isd,ied ! ### change stride order; j->jed-1?
    !do i=isd,ied ; do J=isd,jed-1
      fluxes%frac_shelf_v(i,J) = 0.0
      if ((G%areaT(i,j) + G%areaT(i,j+1) > 0.0)) & ! .and. (G%dxdy_v(i,J) > 0.0)) &
        fluxes%frac_shelf_v(i,J) = ((CS%area_shelf_h(i,j) + CS%area_shelf_h(i,j+1)) / &
                                    (G%areaT(i,j) + G%areaT(i,j+1)))
      fluxes%rigidity_ice_v(i,J) = (CS%kv_ice / CS%density_ice) * &
                                    max(CS%mass_shelf(i,j), CS%mass_shelf(i,j+1))
    enddo ; enddo
    call pass_vector(fluxes%frac_shelf_u, fluxes%frac_shelf_v, G%domain, TO_ALL, CGRID_NE)
  else
    ! This is needed because rigidity is potentially modified in the coupler. Reset
    ! in the ice shelf cavity: MJH

    do j=jsd,jed ; do i=isd,ied-1 ! changed stride
      fluxes%rigidity_ice_u(I,j) = (CS%kv_ice / CS%density_ice) * &
                    min(CS%mass_shelf(i,j), CS%mass_shelf(i+1,j))
    enddo ; enddo

    do j=jsd,jed-1 ; do i=isd,ied ! changed stride
      fluxes%rigidity_ice_v(i,J) = (CS%kv_ice / CS%density_ice) * &
                    max(CS%mass_shelf(i,j), CS%mass_shelf(i,j+1))
    enddo ; enddo
  endif

  if (CS%debug) then
    if (associated(state%taux_shelf) .and. associated(state%tauy_shelf)) then
      call uchksum(state%taux_shelf, "taux_shelf", G%HI, haloshift=0)
      call vchksum(state%tauy_shelf, "tauy_shelf", G%HI, haloshift=0)
      call uchksum(fluxes%rigidity_ice_u, "rigidity_ice_u", G%HI, haloshift=0)
      call vchksum(fluxes%rigidity_ice_v, "rigidity_ice_v", G%HI, haloshift=0)
      call uchksum(fluxes%frac_shelf_u, "frac_shelf_u", G%HI, haloshift=0)
      call vchksum(fluxes%frac_shelf_v, "frac_shelf_v", G%HI, haloshift=0)
    endif
  endif

  if (associated(state%taux_shelf) .and. associated(state%tauy_shelf)) then
    call pass_vector(state%taux_shelf, state%tauy_shelf, G%domain, TO_ALL, CGRID_NE)
  endif

  if (associated(fluxes%sw_vis_dir)) fluxes%sw_vis_dir = 0.0
  if (associated(fluxes%sw_vis_dif)) fluxes%sw_vis_dif = 0.0
  if (associated(fluxes%sw_nir_dir)) fluxes%sw_nir_dir = 0.0
  if (associated(fluxes%sw_nir_dif)) fluxes%sw_nir_dif = 0.0

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    frac_area = fluxes%frac_shelf_h(i,j)
    if (frac_area > 0.0) then
      ! ### THIS SHOULD BE AN AREA WEIGHTED AVERAGE OF THE ustar_shelf POINTS.
      taux2 = 0.0 ; tauy2 = 0.0
      asu1 = fluxes%frac_shelf_u(i-1,j) * (G%areaT(i-1,j) + G%areaT(i,j)) ! G%dxdy_u(i-1,j)
      asu2 = fluxes%frac_shelf_u(i,j) * (G%areaT(i,j) + G%areaT(i+1,j)) ! G%dxdy_u(i,j)
      asv1 = fluxes%frac_shelf_v(i,j-1) * (G%areaT(i,j-1) + G%areaT(i,j)) ! G%dxdy_v(i,j-1)
      asv2 = fluxes%frac_shelf_v(i,j) * (G%areaT(i,j) + G%areaT(i,j+1)) ! G%dxdy_v(i,j)
      if ((asu1 + asu2 > 0.0) .and. associated(state%taux_shelf)) &
        taux2 = (asu1 * state%taux_shelf(i-1,j)**2 + &
                 asu2 * state%taux_shelf(i,j)**2  ) / (asu1 + asu2)
      if ((asv1 + asv2 > 0.0) .and. associated(state%tauy_shelf)) &
        tauy2 = (asv1 * state%tauy_shelf(i,j-1)**2 + &
                 asv2 * state%tauy_shelf(i,j)**2  ) / (asv1 + asv2)

      ! GM: melting is computed using ustar_shelf (and not ustar), which has already
      ! been passed, so believe we do not need to update fluxes%ustar.
      !fluxes%ustar(i,j) = MAX(CS%ustar_bg, sqrt(Irho0 * sqrt(taux2 + tauy2)))


      if (associated(fluxes%sw)) fluxes%sw(i,j) = 0.0
      if (associated(fluxes%lw)) fluxes%lw(i,j) = 0.0
      if (associated(fluxes%latent)) fluxes%latent(i,j) = 0.0
      if (associated(fluxes%evap)) fluxes%evap(i,j) = 0.0
      if (associated(fluxes%lprec)) then
        if (CS%lprec(i,j) > 0.0 ) then
          fluxes%lprec(i,j) =  frac_area*CS%lprec(i,j)*CS%flux_factor
        else
          fluxes%evap(i,j) = frac_area*CS%lprec(i,j)*CS%flux_factor
        endif
      endif


      ! Add frazil formation diagnosed by the ocean model (J m-2) in the
      ! form of surface layer evaporation (kg m-2 s-1). Update lprec in the
      ! control structure for diagnostic purposes.

      if (associated(state%frazil)) then
        fraz = state%frazil(i,j) / CS%time_step / CS%Lat_fusion
        if (associated(fluxes%evap)) fluxes%evap(i,j) = fluxes%evap(i,j) - fraz
        CS%lprec(i,j)=CS%lprec(i,j) - fraz
        state%frazil(i,j) = 0.0
      endif

      if (associated(fluxes%sens)) fluxes%sens(i,j) = -frac_area*CS%t_flux(i,j)*CS%flux_factor
      if (associated(fluxes%salt_flux)) fluxes%salt_flux(i,j) = frac_area * CS%salt_flux(i,j)*CS%flux_factor
      if (associated(fluxes%p_surf)) fluxes%p_surf(i,j) = frac_area * CS%g_Earth * CS%mass_shelf(i,j)
      ! Same for IOB%p
      if (associated(fluxes%p_surf_full) ) fluxes%p_surf_full(i,j) = &
           frac_area * CS%g_Earth * CS%mass_shelf(i,j)

    endif
  enddo ; enddo

  ! If the shelf mass is changing, the fluxes%rigidity_ice_[uv] needs to be
  ! updated here.

  if (CS%shelf_mass_is_dynamic) then
    do j=G%jsc,G%jec ; do i=G%isc-1,G%iec
      fluxes%rigidity_ice_u(I,j) = (CS%kv_ice / CS%density_ice) * &
                                    max(CS%mass_shelf(i,j), CS%mass_shelf(i+1,j))
    enddo ; enddo

    do j=G%jsc-1,G%jec ; do i=G%isc,G%iec
      fluxes%rigidity_ice_v(i,J) = (CS%kv_ice / CS%density_ice) * &
                                    max(CS%mass_shelf(i,j), CS%mass_shelf(i,j+1))
    enddo ; enddo
  endif
end subroutine add_shelf_flux


! subroutine add_shelf_flux_IOB(CS, state, fluxes)
! !  type(ice_ocean_boundary_type),              intent(inout)    :: IOB
!   type(ice_shelf_CS),                 intent(in)    :: CS
!   type(surface),                      intent(inout)    :: state
!   type(forcing),                      intent(inout) :: fluxes
! ! Arguments:
! !  (in)      fluxes - A structure of surface fluxes that may be used.
! !  (in)      visc - A structure containing vertical viscosities, bottom boundary
! !                   layer properies, and related fields.
! !  (in)      G - The ocean's grid structure.
! !  (in)      CS - This module's control structure.
!  !need to use visc variables
!  !time step therm v. dynamic?
!   real :: Irho0         ! The inverse of the mean density in m3 kg-1.
!   real :: frac_area     ! The fractional area covered by the ice shelf, nondim.
!   real :: taux2, tauy2  ! The squared surface stresses, in Pa.
!   real :: asu1, asu2    ! Ocean areas covered by ice shelves at neighboring u-
!   real :: asv1, asv2    ! and v-points, in m2.
!   integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
!   type(ocean_grid_type), pointer :: G

!   G=>CS%grid
!   is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
!   isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed

!   Irho0 = 1.0 / CS%Rho0
!   ! Determine ustar and the square magnitude of the velocity in the
!   ! bottom boundary layer. Together these give the TKE source and
!   ! vertical decay scale.
!   if (CS%shelf_mass_is_dynamic) then
!     do j=jsd,jed ; do i=isd,ied
!       if (G%areaT(i,j) > 0.0) &
!         fluxes%frac_shelf_h(i,j) = CS%area_shelf_h(i,j) / G%areaT(i,j)
!     enddo ; enddo
!     !do I=isd,ied-1 ; do j=isd,jed
!     do j=jsd,jed ; do i=isd,ied-1 ! ### changed stride order; i->ied-1?
!       fluxes%frac_shelf_u(I,j) = 0.0
!       if ((G%areaT(i,j) + G%areaT(i+1,j) > 0.0)) & ! .and. (G%dxdy_u(I,j) > 0.0)) &
!         fluxes%frac_shelf_u(I,j) = ((CS%area_shelf_h(i,j) + CS%area_shelf_h(i+1,j)) / &
!                                     (G%areaT(i,j) + G%areaT(i+1,j)))
!       fluxes%rigidity_ice_u(I,j) = (CS%kv_ice / CS%density_ice) * &
!                                     min(CS%mass_shelf(i,j), CS%mass_shelf(i+1,j))
!     enddo ; enddo
!     do j=jsd,jed-1 ; do i=isd,ied ! ### change stride order; j->jed-1?
!     !do i=isd,ied ; do J=isd,jed-1
!       fluxes%frac_shelf_v(i,J) = 0.0
!       if ((G%areaT(i,j) + G%areaT(i,j+1) > 0.0)) & ! .and. (G%dxdy_v(i,J) > 0.0)) &
!         fluxes%frac_shelf_v(i,J) = ((CS%area_shelf_h(i,j) + CS%area_shelf_h(i,j+1)) / &
!                                     (G%areaT(i,j) + G%areaT(i,j+1)))
!       fluxes%rigidity_ice_v(i,J) = (CS%kv_ice / CS%density_ice) * &
!                                     min(CS%mass_shelf(i,j), CS%mass_shelf(i,j+1))
!     enddo ; enddo
!     call pass_vector(fluxes%frac_shelf_u, fluxes%frac_shelf_v, G%domain, TO_ALL, CGRID_NE)
!   endif

!   if (CS%debug) then
!     if (associated(state%taux_shelf)) then
!       call uchksum(state%taux_shelf, "taux_shelf", G%HI, haloshift=0)
!     endif
!     if (associated(state%tauy_shelf)) then
!       call vchksum(state%tauy_shelf, "tauy_shelf", G%HI, haloshift=0)
!     endif
!   endif

!   if (associated(state%taux_shelf) .and. associated(state%tauy_shelf)) then
!     call pass_vector(state%taux_shelf, state%tauy_shelf, G%domain, TO_ALL, CGRID_NE)
!   endif

!   do j=G%jsc,G%jec ; do i=G%isc,G%iec
!     frac_area = fluxes%frac_shelf_h(i,j)
!     if (frac_area > 0.0) then
!       ! ### THIS SHOULD BE AN AREA WEIGHTED AVERAGE OF THE ustar_shelf POINTS.
!       taux2 = 0.0 ; tauy2 = 0.0
!       asu1 = fluxes%frac_shelf_u(i-1,j) * (G%areaT(i-1,j) + G%areaT(i,j)) ! G%dxdy_u(i-1,j)
!       asu2 = fluxes%frac_shelf_u(i,j) * (G%areaT(i,j) + G%areaT(i+1,j)) ! G%dxdy_u(i,j)
!       asv1 = fluxes%frac_shelf_v(i,j-1) * (G%areaT(i,j-1) + G%areaT(i,j)) ! G%dxdy_v(i,j-1)
!       asv2 = fluxes%frac_shelf_v(i,j) * (G%areaT(i,j) + G%areaT(i,j+1)) ! G%dxdy_v(i,j)
!       if ((asu1 + asu2 > 0.0) .and. associated(state%taux_shelf)) &
!         taux2 = (asu1 * state%taux_shelf(i-1,j)**2 + &
!                  asu2 * state%taux_shelf(i,j)**2  ) / (asu1 + asu2)
!       if ((asv1 + asv2 > 0.0) .and. associated(state%tauy_shelf)) &
!         tauy2 = (asv1 * state%tauy_shelf(i,j-1)**2 + &
!                  asv2 * state%tauy_shelf(i,j)**2  ) / (asv1 + asv2)
!       fluxes%ustar_shelf(i,j) = MAX(CS%ustar_bg, sqrt(Irho0 * sqrt(taux2 + tauy2)))

!       if (CS%lprec(i,j) > 0.0) then
!         fluxes%lprec(i,j) = fluxes%lprec(i,j) + frac_area*CS%lprec(i,j)
!         ! Same for IOB%lprec
!       else
!         fluxes%evap(i,j) = fluxes%evap(i,j) + frac_area*CS%lprec(i,j)
!         ! Same for -1*IOB%q_flux
!       endif
!       fluxes%sens(i,j) = fluxes%sens(i,j) - frac_area*CS%t_flux(i,j)
!       ! Same for -1*IOB%t_flux
!     ! fluxes%salt_flux(i,j) = fluxes%salt_flux(i,j) + frac_area * CS%salt_flux(i,j)
!     ! ! Same for IOB%salt_flux.
!       fluxes%p_surf(i,j) = fluxes%p_surf(i,j) + &
!                            frac_area * CS%g_Earth * CS%mass_shelf(i,j)
!       ! Same for IOB%p
!       if (associated(fluxes%p_surf_full)) fluxes%p_surf_full(i,j) = &
!            fluxes%p_surf_full(i,j) + frac_area * CS%g_Earth * CS%mass_shelf(i,j)
!     endif
!   enddo ; enddo

!   if (CS%debug) then
!     call hchksum(fluxes%ustar_shelf, "ustar_shelf", G%HI, haloshift=0)
!   endif

!   ! If the shelf mass is changing, the fluxes%rigidity_ice_[uv] needs to be
!   ! updated here.

!   if (CS%shelf_mass_is_dynamic) then
!     do j=G%jsc,G%jec ; do i=G%isc-1,G%iec
!       fluxes%rigidity_ice_u(I,j) = (CS%kv_ice / CS%density_ice) * &
!                                     min(CS%mass_shelf(i,j), CS%mass_shelf(i+1,j))
!     enddo ; enddo

!     do j=G%jsc-1,G%jec ; do i=G%isc,G%iec
!       fluxes%rigidity_ice_v(i,J) = (CS%kv_ice / CS%density_ice) * &
!                                     min(CS%mass_shelf(i,j), CS%mass_shelf(i,j+1))
!     enddo ; enddo
!   endif
! end subroutine add_shelf_flux_IOB


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! shelf_model_init - initializes shelf model data, parameters and diagnostics  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine initialize_ice_shelf(param_file, ocn_grid, Time, CS, diag, fluxes, Time_in, solo_ice_sheet_in)
  type(param_file_type), intent(in) :: param_file
  type(ocean_grid_type), pointer    :: ocn_grid
  type(time_type),    intent(inout)   :: Time
  type(ice_shelf_CS), pointer         :: CS
  type(diag_ctrl), target, intent(in) :: diag
  type(forcing), optional, intent(inout) :: fluxes
  type(time_type), optional, intent(in)  :: Time_in
  logical, optional,intent(in)         :: solo_ice_sheet_in

  type(ocean_grid_type), pointer :: G, OG ! Convenience pointers
  type(directories)  :: dirs
  type(vardesc) :: vd
  type(dyn_horgrid_type), pointer :: dG => NULL()
  real :: cdrag, drag_bg_vel
  logical :: new_sim, save_IC, var_force
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=200) :: config
  character(len=200) :: IC_file,filename,inputdir
  character(len=40)  :: var_name
  character(len=40)  :: mod = "MOM_ice_shelf"  ! This module's name.
  character(len=2)   :: procnum
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq, iters
  integer :: wd_halos(2)
  logical :: solo_ice_sheet, read_TideAmp
  character(len=240) :: Tideamp_file
  real    :: utide
  if (associated(CS)) then
    call MOM_error(FATAL, "MOM_ice_shelf.F90, initialize_ice_shelf: "// &
                          "called with an associated control structure.")
    return
  endif
  allocate(CS)

  !   Go through all of the infrastructure initialization calls, since this is
  ! being treated as an independent component that just happens to use the
  ! MOM's grid and infrastructure.
  call Get_MOM_Input(dirs=dirs)

  ! Set up the ice-shelf domain and grid
  wd_halos(:)=0
  call MOM_domains_init(CS%grid%domain, param_file, min_halo=wd_halos, symmetric=GRID_SYM_)
! call diag_mediator_init(CS%grid,param_file,CS%diag) ! this needs to be fixed - will probably break when not using coupled driver 0
  call MOM_grid_init(CS%grid, param_file)

  call create_dyn_horgrid(dG, CS%grid%HI)
  call clone_MOM_domain(CS%grid%Domain, dG%Domain)

  call set_grid_metrics(dG, param_file)
! call set_diag_mediator_grid(CS%grid, CS%diag)

  ! The ocean grid is possibly different
  if (associated(ocn_grid)) CS%ocn_grid => ocn_grid

  ! Convenience pointers
  G => CS%grid
  OG => CS%ocn_grid

  if (is_root_pe()) then
   write(0,*) 'OG: ', OG%isd, OG%isc, OG%iec, OG%ied, OG%jsd, OG%jsc, OG%jsd, OG%jed
   write(0,*) 'IG: ', G%isd, G%isc, G%iec, G%ied, G%jsd, G%jsc, G%jsd, G%jed
  endif

  CS%Time = Time ! ### This might not be in the right place?
  CS%diag => diag

  ! Are we being called from the solo ice-sheet driver? When called by the ocean
  ! model solo_ice_sheet_in is not preset.
  solo_ice_sheet = .false.
  if (present(solo_ice_sheet_in)) solo_ice_sheet = solo_ice_sheet_in
  CS%solo_ice_sheet = solo_ice_sheet

  if (present(Time_in)) Time = Time_in

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
  Isdq = G%IsdB ; Iedq = G%IedB ; Jsdq = G%JsdB ; Jedq = G%JedB

  CS%Lat_fusion = 3.34e5
  CS%override_shelf_movement = .false.

  CS%use_reproducing_sums = .false.
  CS%switch_var = .false.

  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "DEBUG_IS", CS%debug, default=.false.)
  call get_param(param_file, mod, "DYNAMIC_SHELF_MASS", CS%shelf_mass_is_dynamic, &
                 "If true, the ice sheet mass can evolve with time.", &
                 default=.false.)
  if (CS%shelf_mass_is_dynamic) then
    call get_param(param_file, mod, "OVERRIDE_SHELF_MOVEMENT", CS%override_shelf_movement, &
                 "If true, user provided code specifies the ice-shelf \n"//&
                 "movement instead of the dynamic ice model.", default=.false.)
    call get_param(param_file, mod, "GROUNDING_LINE_INTERPOLATE", CS%GL_regularize, &
                 "THIS PARAMETER NEEDS A DESCRIPTION.", default=.false.)
    call get_param(param_file, mod, "GROUNDING_LINE_INTERP_SUBGRID_N", CS%n_sub_regularize, &
                 "THIS PARAMETER NEEDS A DESCRIPTION.", default=0)
    call get_param(param_file, mod, "GROUNDING_LINE_COUPLE", CS%GL_couple, &
                 "THIS PARAMETER NEEDS A DESCRIPTION.", default=.false.)
    if (CS%GL_regularize) CS%GL_couple = .false.
    if (CS%GL_regularize .and. (CS%n_sub_regularize.eq.0)) call MOM_error (FATAL, &
      "GROUNDING_LINE_INTERP_SUBGRID_N must be a positive integer if GL regularization is used")
  endif
  call get_param(param_file, mod, "SHELF_THERMO", CS%isthermo, &
                 "If true, use a thermodynamically interactive ice shelf.", &
                 default=.false.)
  call get_param(param_file, mod, "SHELF_THREE_EQN", CS%threeeq, &
                 "If true, use the three equation expression of \n"//&
                 "consistency to calculate the fluxes at the ice-ocean \n"//&
                 "interface.", default=.true.)
  call get_param(param_file, mod, "SHELF_INSULATOR", CS%insulator, &
                 "If true, the ice shelf is a perfect insulatior \n"//&
                 "(no conduction).", default=.false.)
  call get_param(param_file, mod, "MELTING_CUTOFF_DEPTH", CS%cutoff_depth, &
                 "Depth above which the melt is set to zero (it must be >= 0) \n"//&
                 "Default value won't affect the solution.", default=0.0)
  if (CS%cutoff_depth < 0.) &
     call MOM_error(WARNING,"Initialize_ice_shelf: MELTING_CUTOFF_DEPTH must be >= 0.")

  call get_param(param_file, mod, "SHELF_3EQ_GAMMA", CS%const_gamma, &
                 "If true, user specifies a constant nondimensional heat-transfer coefficient \n"//&
                 "(GAMMA_T_3EQ), from which the salt-transfer coefficient is then computed \n"//&
                 " as GAMMA_T_3EQ/35. This is used with SHELF_THREE_EQN.", default=.false.)
  if (CS%const_gamma) call get_param(param_file, mod, "SHELF_3EQ_GAMMA_T", CS%Gamma_T_3EQ, &
                 "Nondimensional heat-transfer coefficient.",default=2.2E-2, &
                  units="nondim.", fail_if_missing=.true.)
  if (CS%threeeq) &
    call get_param(param_file, mod, "SHELF_S_ROOT", CS%find_salt_root, &
                 "If SHELF_S_ROOT = True, salinity at the ice/ocean interface (Sbdry) \n "//&
                 "is computed from a quadratic equation. Otherwise, the previous \n"//&
                 "interactive method to estimate Sbdry is used.", default=.false.)
  if (CS%find_salt_root) then ! read liquidus coeffs.
     call get_param(param_file, mod, "TFREEZE_S0_P0",CS%lambda1, &
                 "this is the freezing potential temperature at \n"//&
                 "S=0, P=0.", units="deg C", default=0.0, do_not_log=.true.)
    call get_param(param_file, mod, "DTFREEZE_DS",CS%lambda1, &
                 "this is the derivative of the freezing potential \n"//&
                 "temperature with salinity.", &
                 units="deg C PSU-1", default=-0.054, do_not_log=.true.)
    call get_param(param_file, mod, "DTFREEZE_DP",CS%lambda3, &
                 "this is the derivative of the freezing potential \n"//&
                 "temperature with pressure.", &
                 units="deg C Pa-1", default=0.0, do_not_log=.true.)

  endif

  if (.not.CS%threeeq) &
    call get_param(param_file, mod, "SHELF_2EQ_GAMMA_T", CS%gamma_t, &
                 "If SHELF_THREE_EQN is false, this the fixed turbulent \n"//&
                 "exchange velocity at the ice-ocean interface.", &
                 units="m s-1", fail_if_missing=.true.)

  call get_param(param_file, mod, "G_EARTH", CS%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)
  call get_param(param_file, mod, "C_P", CS%Cp, &
                 "The heat capacity of sea water.", units="J kg-1 K-1", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0) !### MAKE THIS A SEPARATE PARAMETER.
  call get_param(param_file, mod, "C_P_ICE", CS%Cp_ice, &
                 "The heat capacity of ice.", units="J kg-1 K-1", &
                 default=2.10e3)

  call get_param(param_file, mod, "ICE_SHELF_FLUX_FACTOR", CS%flux_factor, &
                 "Non-dimensional factor applied to shelf thermodynamic \n"//&
                 "fluxes.", units="none", default=1.0)

  call get_param(param_file, mod, "KV_ICE", CS%kv_ice, &
                 "The viscosity of the ice.", units="m2 s-1", default=1.0e10)
  call get_param(param_file, mod, "KV_MOLECULAR", CS%kv_molec, &
                 "The molecular kinimatic viscosity of sea water at the \n"//&
                 "freezing temperature.", units="m2 s-1", default=1.95e-6)
  call get_param(param_file, mod, "ICE_SHELF_SALINITY", CS%Salin_ice, &
                 "The salinity of the ice inside the ice shelf.", units="PSU", &
                 default=0.0)
  call get_param(param_file, mod, "ICE_SHELF_TEMPERATURE", CS%Temp_ice, &
                 "The temperature at the center of the ice shelf.", &
                 units = "degC", default=-15.0)
  call get_param(param_file, mod, "KD_SALT_MOLECULAR", CS%kd_molec_salt, &
                 "The molecular diffusivity of salt in sea water at the \n"//&
                 "freezing point.", units="m2 s-1", default=8.02e-10)
  call get_param(param_file, mod, "KD_TEMP_MOLECULAR", CS%kd_molec_temp, &
                 "The molecular diffusivity of heat in sea water at the \n"//&
                 "freezing point.", units="m2 s-1", default=1.41e-7)
  call get_param(param_file, mod, "RHO_0", CS%density_ocean_avg, &
                 "avg ocean density used in floatation cond", &
                 units="kg m-3", default=1035.)
  call get_param(param_file, mod, "DT_FORCING", CS%time_step, &
                 "The time step for changing forcing, coupling with other \n"//&
                 "components, or potentially writing certain diagnostics. \n"//&
                 "The default value is given by DT.", units="s", default=0.0)
  call get_param(param_file, mod, "SHELF_DIAG_TIMESTEP", CS%velocity_update_time_step, &
                 "A timestep to use for diagnostics of the shelf.", default=0.0)

  call get_param(param_file, mod, "COL_THICK_MELT_THRESHOLD", CS%col_thick_melt_threshold, &
                 "The minimum ML thickness where melting is allowed.", units="m", &
                 default=0.0)

  call get_param(param_file, mod, "READ_TIDEAMP", read_TIDEAMP, &
                 "If true, read a file (given by TIDEAMP_FILE) containing \n"//&
                 "the tidal amplitude with INT_TIDE_DISSIPATION.", default=.false.)

  call safe_alloc_ptr(CS%utide,isd,ied,jsd,jed)   ; CS%utide(:,:) = 0.0

  if (read_TIDEAMP) then
    call get_param(param_file, mod, "TIDEAMP_FILE", TideAmp_file, &
                 "The path to the file containing the spatially varying \n"//&
                 "tidal amplitudes.", &
                 default="tideamp.nc")
    call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    TideAmp_file = trim(inputdir) // trim(TideAmp_file)
    call read_data(TideAmp_file,'tideamp',CS%utide,domain=G%domain%mpp_domain,timelevel=1)
  else
    call get_param(param_file, mod, "UTIDE", utide, &
                 "The constant tidal amplitude used with INT_TIDE_DISSIPATION.", &
                 units="m s-1", default=0.0)
    CS%utide = utide
  endif

  call EOS_init(param_file, CS%eqn_of_state)

  !! new parameters that need to be in MOM_input

  if (CS%shelf_mass_is_dynamic .and. .not.CS%override_shelf_movement) then

    call get_param(param_file, mod, "A_GLEN_ISOTHERM", CS%A_glen_isothermal, &
                 "Ice viscosity parameter in Glen's Law", &
                 units="Pa -1/3 a", default=9.461e-18)
    call get_param(param_file, mod, "GLEN_EXPONENT", CS%n_glen, &
                 "nonlinearity exponent in Glen's Law", &
                  units="none", default=3.)
    call get_param(param_file, mod, "MIN_STRAIN_RATE_GLEN", CS%eps_glen_min, &
                 "min. strain rate to avoid infinite Glen's law viscosity", &
                 units="a-1", default=1.e-12)
    call get_param(param_file, mod, "BASAL_FRICTION_COEFF", CS%C_basal_friction, &
                 "ceofficient in sliding law \tau_b = C u^(n_basal_friction)", &
                 units="Pa (m-a)-(n_basal_friction)", fail_if_missing=.true.)
    call get_param(param_file, mod, "BASAL_FRICTION_EXP", CS%n_basal_friction, &
                 "exponent in sliding law \tau_b = C u^(m_slide)", &
                 units="none", fail_if_missing=.true.)
    call get_param(param_file, mod, "DENSITY_ICE", CS%density_ice, &
                 "A typical density of ice.", units="kg m-3", default=917.0)

    call get_param(param_file, mod, "INPUT_FLUX_ICE_SHELF", CS%input_flux, &
                 "volume flux at upstream boundary", &
                 units="m2 s-1", default=0.)
    call get_param(param_file, mod, "INPUT_THICK_ICE_SHELF", CS%input_thickness, &
                 "flux thickness at upstream boundary", &
                 units="m", default=1000.)
    call get_param(param_file, mod, "ICE_VELOCITY_TIMESTEP", CS%velocity_update_time_step, &
                 "seconds between ice velocity calcs", units="s", &
                 fail_if_missing=.true.)

    call get_param(param_file, mod, "CONJUGATE_GRADIENT_TOLERANCE", CS%cg_tolerance, &
        "tolerance in CG solver, relative to initial residual", default=1.e-6)
    call get_param(param_file, mod, "ICE_NONLINEAR_TOLERANCE", &
        CS%nonlinear_tolerance,"nonlin tolerance in iterative velocity solve",default=1.e-6)
    call get_param(param_file, mod, "CONJUGATE_GRADIENT_MAXIT", CS%cg_max_iterations, &
        "max iteratiions in CG solver", default=2000)
    call get_param(param_file, mod, "THRESH_FLOAT_COL_DEPTH", CS%thresh_float_col_depth, &
        "min ocean thickness to consider ice *floating*; \n"// &
        "will only be important with use of tides", &
        units="m",default=1.e-3)

    call get_param(param_file, mod, "SHELF_MOVING_FRONT", CS%moving_shelf_front, &
                 "whether or not to advance shelf front (and calve..)")
    call get_param(param_file, mod, "CALVE_TO_MASK", CS%calve_to_mask, &
                 "if true, do not allow an ice shelf where prohibited by a mask")
    call get_param(param_file, mod, "MIN_THICKNESS_SIMPLE_CALVE", CS%min_thickness_simple_calve, &
                 "min thickness rule for VERY simple calving law",&
                 units="m", default=0.0)

    call get_param(param_file, mod, "ICE_SHELF_CFL_FACTOR", CS%CFL_factor, &
        "limit timestep as a factor of min (\Delta x / u); \n"// &
        "only important for ice-only model", &
        default=0.25)
    call get_param(param_file, mod, "NONLIN_SOLVE_ERR_MODE", CS%nonlin_solve_err_mode, &
        "choose whether nonlin error in vel solve is based on nonlinear residual (1) \n"// &
        "or relative change since last iteration (2)", &
        default=1)


    if (CS%debug) CS%use_reproducing_sums = .true.

    CS%nstep_velocity = FLOOR (CS%velocity_update_time_step / CS%time_step)
    CS%velocity_update_counter = 0
    CS%velocity_update_sub_counter = 0
  else
    CS%nstep_velocity = 0
    ! This is here because of inconsistent defaults.  I don't know why.  RWH
    call get_param(param_file, mod, "DENSITY_ICE", CS%density_ice, &
                 "A typical density of ice.", units="kg m-3", default=900.0)
  endif
  call get_param(param_file, mod, "WRITE_OUTPUT_TO_FILE", &
        CS%write_output_to_file, "for debugging purposes",default=.false.)

  call get_param(param_file, mod, "USTAR_SHELF_BG", CS%ustar_bg, &
                 "The minimum value of ustar under ice sheves.", units="m s-1", &
                 default=0.0)
  call get_param(param_file, mod, "CDRAG_SHELF", cdrag, &
       "CDRAG is the drag coefficient relating the magnitude of \n"//&
       "the velocity field to the surface stress.", units="nondim", &
       default=0.003)
  CS%cdrag = cdrag
  if (CS%ustar_bg <= 0.0) then
    call get_param(param_file, mod, "DRAG_BG_VEL_SHELF", drag_bg_vel, &
                 "DRAG_BG_VEL is either the assumed bottom velocity (with \n"//&
                 "LINEAR_DRAG) or an unresolved  velocity that is \n"//&
                 "combined with the resolved velocity to estimate the \n"//&
                 "velocity magnitude.", units="m s-1", default=0.0)
    if (CS%cdrag*drag_bg_vel > 0.0) CS%ustar_bg = sqrt(CS%cdrag)*drag_bg_vel

  endif

  ! Allocate  and initialize variables
  allocate( CS%mass_shelf(isd:ied,jsd:jed) )   ; CS%mass_shelf(:,:) = 0.0
  allocate( CS%area_shelf_h(isd:ied,jsd:jed) ) ; CS%area_shelf_h(:,:) = 0.0
  allocate( CS%t_flux(isd:ied,jsd:jed) )       ; CS%t_flux(:,:) = 0.0
  allocate( CS%lprec(isd:ied,jsd:jed) )        ; CS%lprec(:,:) = 0.0
  allocate( CS%salt_flux(isd:ied,jsd:jed) )    ; CS%salt_flux(:,:) = 0.0

  allocate( CS%tflux_shelf(isd:ied,jsd:jed) ) ; CS%tflux_shelf(:,:) = 0.0
  allocate( CS%tfreeze(isd:ied,jsd:jed) )     ; CS%tfreeze(:,:) = 0.0
  allocate( CS%exch_vel_s(isd:ied,jsd:jed) )  ; CS%exch_vel_s(:,:) = 0.0
  allocate( CS%exch_vel_t(isd:ied,jsd:jed) )  ; CS%exch_vel_t(:,:) = 0.0

  allocate ( CS%h_shelf(isd:ied,jsd:jed) )   ; CS%h_shelf(:,:) = 0.0
  allocate ( CS%hmask(isd:ied,jsd:jed) )   ; CS%hmask(:,:) = -2.0


  ! OVS vertically integrated Temperature
  allocate ( CS%t_shelf(isd:ied,jsd:jed) )   ; CS%t_shelf(:,:) = -10.0
  allocate ( CS%t_boundary_values(isd:ied,jsd:jed) )   ; CS%t_boundary_values(:,:) = -15.0
  allocate ( CS%tmask(Isdq:Iedq,Jsdq:Jedq) ) ; CS%tmask(:,:) = -1.0

  if (CS%shelf_mass_is_dynamic .and. .not.CS%override_shelf_movement) then
    ! DNG
    allocate ( CS%u_shelf(Isdq:Iedq,Jsdq:Jedq) ) ; CS%u_shelf(:,:) = 0.0
    allocate ( CS%v_shelf(Isdq:Iedq,Jsdq:Jedq) ) ; CS%v_shelf(:,:) = 0.0
    allocate ( CS%u_boundary_values(Isdq:Iedq,Jsdq:Jedq) ) ; CS%u_boundary_values(:,:) = 0.0
    allocate ( CS%v_boundary_values(Isdq:Iedq,Jsdq:Jedq) ) ; CS%v_boundary_values(:,:) = 0.0
    allocate ( CS%h_boundary_values(isd:ied,jsd:jed) ) ; CS%h_boundary_values(:,:) = 0.0
    allocate ( CS%thickness_boundary_values(isd:ied,jsd:jed) ) ; CS%thickness_boundary_values(:,:) = 0.0
    allocate ( CS%ice_visc_bilinear(isd:ied,jsd:jed) ) ; CS%ice_visc_bilinear(:,:) = 0.0
    allocate ( CS%ice_visc_lower_tri(isd:ied,jsd:jed) ) ; CS%ice_visc_lower_tri = 0.0
    allocate ( CS%ice_visc_upper_tri(isd:ied,jsd:jed) ) ; CS%ice_visc_upper_tri = 0.0
    allocate ( CS%u_face_mask(Isdq:Iedq,jsd:jed) ) ; CS%u_face_mask(:,:) = 0.0
    allocate ( CS%v_face_mask(isd:ied,Jsdq:Jedq) ) ; CS%v_face_mask(:,:) = 0.0
    allocate ( CS%u_face_mask_boundary(Isdq:Iedq,jsd:jed) ) ; CS%u_face_mask_boundary(:,:) = -2.0
    allocate ( CS%v_face_mask_boundary(isd:ied,Jsdq:Jedq) ) ; CS%v_face_mask_boundary(:,:) = -2.0
    allocate ( CS%u_flux_boundary_values(Isdq:Iedq,jsd:jed) ) ; CS%u_flux_boundary_values(:,:) = 0.0
    allocate ( CS%v_flux_boundary_values(isd:ied,Jsdq:Jedq) ) ; CS%v_flux_boundary_values(:,:) = 0.0
    allocate ( CS%umask(Isdq:Iedq,Jsdq:Jedq) ) ; CS%umask(:,:) = -1.0
    allocate ( CS%vmask(Isdq:Iedq,Jsdq:Jedq) ) ; CS%vmask(:,:) = -1.0

    allocate ( CS%taub_beta_eff_bilinear(isd:ied,jsd:jed) ) ; CS%taub_beta_eff_bilinear(:,:) = 0.0
    allocate ( CS%taub_beta_eff_upper_tri(isd:ied,jsd:jed) ) ; CS%taub_beta_eff_upper_tri(:,:) = 0.0
    allocate ( CS%taub_beta_eff_lower_tri(isd:ied,jsd:jed) ) ; CS%taub_beta_eff_lower_tri(:,:) = 0.0
    allocate ( CS%OD_rt(isd:ied,jsd:jed) ) ; CS%OD_rt(:,:) = 0.0
    allocate ( CS%OD_av(isd:ied,jsd:jed) ) ; CS%OD_av(:,:) = 0.0
    allocate ( CS%float_frac(isd:ied,jsd:jed) ) ; CS%float_frac(:,:) = 0.0
    allocate ( CS%float_frac_rt(isd:ied,jsd:jed) ) ; CS%float_frac_rt(:,:) = 0.0

    if (CS%calve_to_mask) then
      allocate ( CS%calve_mask (isd:ied,jsd:jed) ) ; CS%calve_mask(:,:) = 0.0
    endif

  endif

  ! Allocate the arrays for passing ice-shelf data through the forcing type.
  if (.not. solo_ice_sheet) then
    if (is_root_pe())  print *,"initialize_ice_shelf: allocating fluxes"
       ! GM: the following assures that water/heat fluxes are just allocated
       ! when SHELF_THERMO = True. These fluxes are necessary if one wants to
       ! use either ENERGETICS_SFC_PBL (ALE mode) or BULKMIXEDLAYER (layer mode).
       call allocate_forcing_type(G, fluxes, ustar=.true., shelf=.true., &
                                 press=.true., water=CS%isthermo, heat=CS%isthermo)
  else
    if (is_root_pe())  print *,"allocating fluxes in solo mode"
    call allocate_forcing_type(G, fluxes, ustar=.true., shelf=.true., press=.true.)
  endif

! Set up the bottom depth, G%D either analytically or from file
  call MOM_initialize_topography(G%bathyT, G%max_depth, dG, param_file)
! Set up the Coriolis parameter, G%f, usually analytically.
  call MOM_initialize_rotation(G%CoriolisBu, dG, param_file)
  call copy_dyngrid_to_MOM_grid(dG, CS%grid)

  call destroy_dyn_horgrid(dG)

  ! Set up the restarts.
  call restart_init(param_file, CS%restart_CSp, "Shelf.res")
  vd = var_desc("shelf_mass","kg m-2","Ice shelf mass",z_grid='1')
  call register_restart_field(CS%mass_shelf, vd, .true., CS%restart_CSp)
  vd = var_desc("shelf_area","m2","Ice shelf area in cell",z_grid='1')
  call register_restart_field(CS%area_shelf_h, vd, .true., CS%restart_CSp)
  vd = var_desc("h_shelf","m","ice sheet/shelf thickness",z_grid='1')
  call register_restart_field(CS%h_shelf, vd, .true., CS%restart_CSp)

  if (CS%shelf_mass_is_dynamic .and. .not.CS%override_shelf_movement) then
    ! additional restarts for ice shelf state
    vd = var_desc("u_shelf","m s-1","ice sheet/shelf velocity",'q',z_grid='1')
    call register_restart_field(CS%u_shelf, vd, .true., CS%restart_CSp)
    vd = var_desc("v_shelf","m s-1","ice sheet/shelf velocity",'q',z_grid='1')
    call register_restart_field(CS%v_shelf, vd, .true., CS%restart_CSp)
    !vd = var_desc("h_shelf","m","ice sheet/shelf thickness",z_grid='1')
    !call register_restart_field(CS%h_shelf, vd, .true., CS%restart_CSp)

    vd = var_desc("h_mask","none","ice sheet/shelf thickness mask",z_grid='1')
    call register_restart_field(CS%hmask, vd, .true., CS%restart_CSp)

    ! OVS vertically integrated stream/shelf temperature
    vd = var_desc("t_shelf","deg C","ice sheet/shelf temperature",z_grid='1')
    call register_restart_field(CS%t_shelf, vd, .true., CS%restart_CSp)


  !  vd = var_desc("area_shelf_h","m-2","ice-covered area of a cell",z_grid='1')
  !  call register_restart_field(CS%area_shelf_h, CS%area_shelf_h, vd, .true., CS%restart_CSp)

    vd = var_desc("OD_av","m","avg ocean depth in a cell",z_grid='1')
    call register_restart_field(CS%OD_av, vd, .true., CS%restart_CSp)

  !  vd = var_desc("OD_av_rt","m","avg ocean depth in a cell, intermed",z_grid='1')
  !  call register_restart_field(CS%OD_av_rt, CS%OD_av_rt, vd, .true., CS%restart_CSp)

    vd = var_desc("float_frac","m","degree of grounding",z_grid='1')
    call register_restart_field(CS%float_frac, vd, .true., CS%restart_CSp)

  !  vd = var_desc("float_frac_rt","m","degree of grounding, intermed",z_grid='1')
  !  call register_restart_field(CS%float_frac_rt, CS%float_frac_rt, vd, .true., CS%restart_CSp)

    vd = var_desc("viscosity","m","glens law ice visc",z_grid='1')
    call register_restart_field(CS%ice_visc_bilinear, vd, .true., CS%restart_CSp)
    vd = var_desc("tau_b_beta","m","coefficient of basal traction",z_grid='1')
    call register_restart_field(CS%taub_beta_eff_bilinear, vd, .true., CS%restart_CSp)
  endif

! GM - I think we do not need to save ustar_shelf and iceshelf_melt in the restart file

!  if (.not. solo_ice_sheet) then
!    vd = var_desc("ustar_shelf","m s-1","Friction velocity under ice shelves",z_grid='1')
!    call register_restart_field(fluxes%ustar_shelf, vd, .true., CS%restart_CSp)
!    vd = var_desc("iceshelf_melt","m year-1","Ice Shelf Melt Rate",z_grid='1')
!    call register_restart_field(fluxes%iceshelf_melt, vd, .true., CS%restart_CSp)
!  endif

  CS%restart_output_dir = dirs%restart_output_dir

  new_sim = .false.
  if ((dirs%input_filename(1:1) == 'n') .and. &
      (LEN_TRIM(dirs%input_filename) == 1)) new_sim = .true.

  if (CS%override_shelf_movement) then
  ! This call is always made because parameters and types may be set
  ! inside, in addition to initializing the mass arrays.
    call initialize_shelf_mass(G, param_file, CS)
!  else if (CS%shelf_mass_is_dynamic) then
!    call initialize_ice_shelf_boundary ( CS%u_face_mask_boundary, CS%v_face_mask_boundary, &
!                                         CS%u_flux_boundary_values, CS%v_flux_boundary_values, &
!                                         CS%u_boundary_values, CS%v_boundary_values, CS%h_boundary_values, &
!                                         CS%hmask, G, param_file)
  end if

  if (CS%shelf_mass_is_dynamic .and. .not. CS%override_shelf_movement) then
    ! the only reason to initialize boundary conds is if the shelf is dynamic

    !MJHcall initialize_ice_shelf_boundary ( CS%u_face_mask_boundary, CS%v_face_mask_boundary, &
    !MJH                                     CS%u_flux_boundary_values, CS%v_flux_boundary_values, &
    !MJH                                     CS%u_boundary_values, CS%v_boundary_values, CS%h_boundary_values, &
    !MJH                                     CS%hmask, G, param_file)

  end if

  if (new_sim .and. .not. CS%override_shelf_movement) then

    ! This model is initialized internally or from a file.
    call initialize_ice_thickness (CS%h_shelf, CS%area_shelf_h, CS%hmask, G, param_file)

    ! next make sure mass is consistent with thickness
    do j=G%jsd,G%jed
      do i=G%isd,G%ied
        if ((CS%hmask(i,j) .eq. 1) .or. (CS%hmask(i,j) .eq. 2)) then
          CS%mass_shelf(i,j) = CS%h_shelf(i,j)*CS%density_ice
        endif
      enddo
    enddo

  ! else ! Previous block for new_sim=.T., this block restores the state.
  elseif (.not.new_sim) then
!    This line calls a subroutine that reads the initial conditions  !
!  from a previously generated (restart?) file.                      !
    call restore_state(dirs%input_filename, dirs%restart_input_dir, Time, &
                       G, CS%restart_CSp)

    ! i think this call isnt necessary - all it does is set hmask to 3 at
    ! the dirichlet boundary, and now this is done elsewhere
    !  call initialize_shelf_mass(G, param_file, CS, .false.)

    if (CS%shelf_mass_is_dynamic .and. .not.CS%override_shelf_movement) then

      ! this is unfortunately necessary; if grid is not symmetric the boundary values
      !  of u and v are otherwise not set till the end of the first linear solve, and so
      !  viscosity is not calculated correctly
      if (.not. G%symmetric) then
        do j=G%jsd,G%jed
          do i=G%isd,G%ied
            if (((i+G%idg_offset) .eq. (G%domain%nihalo+1)).and.(CS%u_face_mask(i-1,j).eq.3)) then
              CS%u_shelf (i-1,j-1) = CS%u_boundary_values (i-1,j-1)
              CS%u_shelf (i-1,j) = CS%u_boundary_values (i-1,j)
            endif
            if (((j+G%jdg_offset) .eq. (G%domain%njhalo+1)).and.(CS%v_face_mask(i,j-1).eq.3)) then
              CS%u_shelf (i-1,j-1) = CS%u_boundary_values (i-1,j-1)
              CS%u_shelf (i,j-1) = CS%u_boundary_values (i,j-1)
            endif
          enddo
        enddo
      endif

      call pass_var (CS%area_shelf_h,G%domain)
      call pass_var (CS%h_shelf,G%domain)
      call pass_var (CS%hmask,G%domain)
      call pass_var (CS%OD_av,G%domain)
      call pass_var (CS%float_frac,G%domain)
      call pass_var (CS%ice_visc_bilinear,G%domain)
      call pass_var (CS%taub_beta_eff_bilinear,G%domain)
      call pass_vector(CS%u_shelf, CS%v_shelf, G%domain, TO_ALL, BGRID_NE)

      if (is_root_pe()) PRINT *, "RESTORING ICE SHELF FROM FILE!!!!!!!!!!!!!"
    endif

  endif ! .not. new_sim

  CS%Time = Time

  ! Transfer the appropriate fields to the forcing type.
  if (CS%shelf_mass_is_dynamic .and. .not.CS%override_shelf_movement) then
    call cpu_clock_begin(id_clock_pass)
    call pass_var(CS%area_shelf_h, G%domain)
    call pass_var(G%bathyT, G%domain)
    call pass_var(CS%hmask, G%domain)

    call update_velocity_masks (CS)

    call pass_var(CS%h_shelf, G%domain)
    call pass_var(CS%mass_shelf, G%domain)
    call cpu_clock_end(id_clock_pass)
  endif
    call pass_var(CS%area_shelf_h, G%domain)

  do j=jsd,jed ; do i=isd,ied ! changed stride
    if (CS%area_shelf_h(i,j) > G%areaT(i,j)) then
      call MOM_error(WARNING,"Initialize_ice_shelf: area_shelf_h exceeds G%areaT.")
      CS%area_shelf_h(i,j) = G%areaT(i,j)
    endif
   !if (.not. solo_ice_sheet) then
    if (G%areaT(i,j) > 0.0) fluxes%frac_shelf_h(i,j) = CS%area_shelf_h(i,j) / G%areaT(i,j)
    if (associated(fluxes%p_surf)) &
      fluxes%p_surf(i,j) = fluxes%p_surf(i,j) + &
        fluxes%frac_shelf_h(i,j) * (CS%g_Earth * CS%mass_shelf(i,j))
    if (associated(fluxes%p_surf_full)) &
      fluxes%p_surf_full(i,j) = fluxes%p_surf_full(i,j) + &
        fluxes%frac_shelf_h(i,j) * (CS%g_Earth * CS%mass_shelf(i,j))
   !endif
  enddo ; enddo

  if (.not. solo_ice_sheet) then
    do j=jsd,jed ; do i=isd,ied-1 ! changed stride
    !do I=isd,ied-1 ; do j=isd,jed
    fluxes%frac_shelf_u(I,j) = 0.0
    if ((G%areaT(i,j) + G%areaT(i+1,j) > 0.0)) & ! .and. (G%dxdy_u(I,j) > 0.0)) &
    fluxes%frac_shelf_u(I,j) = ((CS%area_shelf_h(i,j) + CS%area_shelf_h(i+1,j)) / &
                    (G%areaT(i,j) + G%areaT(i+1,j)))
    fluxes%rigidity_ice_u(I,j) = (CS%kv_ice / CS%density_ice) * &
                    min(CS%mass_shelf(i,j), CS%mass_shelf(i+1,j))
    enddo ; enddo


    do j=jsd,jed-1 ; do i=isd,ied ! changed stride
    !do i=isd,ied ; do J=isd,jed-1
    fluxes%frac_shelf_v(i,J) = 0.0
    if ((G%areaT(i,j) + G%areaT(i,j+1) > 0.0)) & ! .and. (G%dxdy_v(i,J) > 0.0)) &
    fluxes%frac_shelf_v(i,J) = ((CS%area_shelf_h(i,j) + CS%area_shelf_h(i,j+1)) / &
                    (G%areaT(i,j) + G%areaT(i,j+1)))
    fluxes%rigidity_ice_v(i,J) = (CS%kv_ice / CS%density_ice) * &
                    min(CS%mass_shelf(i,j), CS%mass_shelf(i,j+1))
    enddo ; enddo
  endif

  !GM, is this needed?
  !write (procnum,'(I2)') mpp_pe()

  if (.not. solo_ice_sheet) then
  call pass_vector(fluxes%frac_shelf_u, fluxes%frac_shelf_v, G%domain, TO_ALL, CGRID_NE)
  endif
 ! call savearray2 ('frac_shelf_u'//procnum,fluxes%frac_shelf_u,CS%write_output_to_file)
 ! call savearray2 ('frac_shelf_v'//procnum,fluxes%frac_shelf_v,CS%write_output_to_file)
 ! call savearray2 ('frac_shelf_h'//procnum,fluxes%frac_shelf_h,CS%write_output_to_file)
 ! call savearray2 ('area_shelf_h'//procnum,CS%area_shelf_h,CS%write_output_to_file)

  ! if we are calving to a mask, i.e. if a mask exists where a shelf cannot, then we read
  ! the mask from a file

  if (CS%shelf_mass_is_dynamic .and. CS%calve_to_mask .and. &
           .not.CS%override_shelf_movement) then

    call MOM_mesg("  MOM_ice_shelf.F90, initialize_ice_shelf: reading calving_mask")

    call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    call get_param(param_file, mod, "CALVING_MASK_FILE", IC_file, &
                 "The file with a mask for where calving might occur.", &
                 default="ice_shelf_h.nc")
    call get_param(param_file, mod, "CALVING_MASK_VARNAME", var_name, &
                 "The variable to use in masking calving.", &
                 default="area_shelf_h")

    filename = trim(inputdir)//trim(IC_file)
    call log_param(param_file, mod, "INPUTDIR/CALVING_MASK_FILE", filename)
    if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " calving mask file: Unable to open "//trim(filename))

    call read_data(filename,trim(var_name),CS%calve_mask,domain=G%Domain%mpp_domain)
    do j=G%jsc,G%jec
      do i=G%isc,G%iec
        if (CS%calve_mask(i,j) > 0.0) CS%calve_mask(i,j) = 1.0
      enddo
    enddo

    call pass_var (CS%calve_mask,G%domain)
  endif

  if (CS%shelf_mass_is_dynamic .and. .not.CS%override_shelf_movement) then
!    call init_boundary_values (CS, time, CS%input_flux, CS%input_thickness, new_sim)

    if (.not. CS%isthermo) then
      CS%lprec(:,:) = 0.0
    endif


    if (new_sim) then
      if (is_root_pe()) print *,"NEW SIM: initialize velocity"
      call update_OD_ffrac_uncoupled (CS)
      call ice_shelf_solve_outer (CS, CS%u_shelf, CS%v_shelf, 1, iters, Time)

!      write (procnum,'(I2)') mpp_pe()

      if (CS%id_u_shelf > 0) call post_data(CS%id_u_shelf,CS%u_shelf,CS%diag)
      if (CS%id_v_shelf > 0) call post_data(CS%id_v_shelf,CS%v_shelf,CS%diag)
    endif
  endif

  call get_param(param_file, mod, "SAVE_INITIAL_CONDS", save_IC, &
                 "If true, save the ice shelf initial conditions.", &
                 default=.false.)
  if (save_IC) call get_param(param_file, mod, "SHELF_IC_OUTPUT_FILE", IC_file,&
                 "The name-root of the output file for the ice shelf \n"//&
                 "initial conditions.", default="MOM_Shelf_IC")

  if (save_IC .and. .not.((dirs%input_filename(1:1) == 'r') .and. &
                          (LEN_TRIM(dirs%input_filename) == 1))) then

    call save_restart(dirs%output_directory, CS%Time, G, &
                      CS%restart_CSp, filename=IC_file)
  endif


  CS%id_area_shelf_h = register_diag_field('ocean_model', 'area_shelf_h', CS%diag%axesT1, CS%Time, &
     'Ice Shelf Area in cell', 'meter-2')
  CS%id_shelf_mass = register_diag_field('ocean_model', 'shelf_mass', CS%diag%axesT1, CS%Time, &
     'mass of shelf', 'kg/m^2')
  CS%id_mass_flux = register_diag_field('ocean_model', 'mass_flux', CS%diag%axesT1,&
     CS%Time,'Total mass flux of freshwater across the ice-ocean interface.', 'kg/s')
  CS%id_melt = register_diag_field('ocean_model', 'melt', CS%diag%axesT1, CS%Time, &
     'Ice Shelf Melt Rate', 'meter year-1')
  CS%id_thermal_driving = register_diag_field('ocean_model', 'thermal_driving', CS%diag%axesT1, CS%Time, &
     'pot. temp. in the boundary layer minus freezing pot. temp. at the ice-ocean interface.', 'Celsius')
  CS%id_haline_driving = register_diag_field('ocean_model', 'haline_driving', CS%diag%axesT1, CS%Time, &
     'salinity in the boundary layer minus salinity at the ice-ocean interface.', 'PPT')
  CS%id_Sbdry = register_diag_field('ocean_model', 'sbdry', CS%diag%axesT1, CS%Time, &
     'salinity at the ice-ocean interface.', 'PPT')
  CS%id_u_ml = register_diag_field('ocean_model', 'u_ml', CS%diag%axesT1, CS%Time, &
     'Eastward vel. in the boundary layer (used to compute ustar)', 'meter second-1')
  CS%id_v_ml = register_diag_field('ocean_model', 'v_ml', CS%diag%axesT1, CS%Time, &
     'Northward vel. in the boundary layer (used to compute ustar)', 'meter second-1')
  CS%id_exch_vel_s = register_diag_field('ocean_model', 'exch_vel_s', CS%diag%axesT1, CS%Time, &
     'Sub-shelf salinity exchange velocity', 'meter second-1')
  CS%id_exch_vel_t = register_diag_field('ocean_model', 'exch_vel_t', CS%diag%axesT1, CS%Time, &
     'Sub-shelf thermal exchange velocity', 'meter second-1')
  CS%id_tfreeze = register_diag_field('ocean_model', 'tfreeze', CS%diag%axesT1, CS%Time, &
     'In Situ Freezing point at ice shelf interface', 'degrees Celsius')
  CS%id_tfl_shelf = register_diag_field('ocean_model', 'tflux_shelf', CS%diag%axesT1, CS%Time, &
     'Heat conduction into ice shelf', 'Watts meter-2')
  CS%id_ustar_shelf = register_diag_field('ocean_model', 'ustar_shelf', CS%diag%axesT1, CS%Time, &
     'Fric vel under shelf', 'm/s')

  if (CS%shelf_mass_is_dynamic .and. .not.CS%override_shelf_movement) then
    CS%id_u_shelf = register_diag_field('ocean_model','u_shelf',CS%diag%axesB1,CS%Time, &
       'x-velocity of ice', 'm year')
    CS%id_v_shelf = register_diag_field('ocean_model','v_shelf',CS%diag%axesB1,CS%Time, &
       'y-velocity of ice', 'm year')
    CS%id_u_mask = register_diag_field('ocean_model','u_mask',CS%diag%axesB1,CS%Time, &
       'mask for u-nodes', 'none')
    CS%id_v_mask = register_diag_field('ocean_model','v_mask',CS%diag%axesB1,CS%Time, &
       'mask for v-nodes', 'none')
    CS%id_h_mask = register_diag_field('ocean_model','h_mask',CS%diag%axesT1,CS%Time, &
       'ice shelf thickness', 'none')
    CS%id_surf_elev = register_diag_field('ocean_model','ice_surf',CS%diag%axesT1,CS%Time, &
       'ice surf elev', 'm')
    CS%id_float_frac = register_diag_field('ocean_model','ice_float_frac',CS%diag%axesT1,CS%Time, &
       'fraction of cell that is floating (sort of)', 'none')
    CS%id_col_thick = register_diag_field('ocean_model','col_thick',CS%diag%axesT1,CS%Time, &
       'ocean column thickness passed to ice model', 'm')
    CS%id_OD_av = register_diag_field('ocean_model','OD_av',CS%diag%axesT1,CS%Time, &
       'intermediate ocean column thickness passed to ice model', 'm')
    CS%id_float_frac_rt = register_diag_field('ocean_model','float_frac_rt',CS%diag%axesT1,CS%Time, &
       'timesteps where cell is floating ', 'none')
    !CS%id_h_after_uflux = register_diag_field('ocean_model','h_after_uflux',CS%diag%axesh1,CS%Time, &
    !   'thickness after u flux ', 'none')
    !CS%id_h_after_vflux = register_diag_field('ocean_model','h_after_vflux',CS%diag%axesh1,CS%Time, &
    !   'thickness after v flux ', 'none')
    !CS%id_h_after_adv = register_diag_field('ocean_model','h_after_adv',CS%diag%axesh1,CS%Time, &
    !   'thickness after front adv ', 'none')

!!! OVS vertically integrated temperature
    CS%id_t_shelf = register_diag_field('ocean_model','t_shelf',CS%diag%axesT1,CS%Time, &
       'T of ice', 'oC')
    CS%id_t_mask = register_diag_field('ocean_model','tmask',CS%diag%axesT1,CS%Time, &
       'mask for T-nodes', 'none')
  endif

  id_clock_shelf = cpu_clock_id('Ice shelf', grain=CLOCK_COMPONENT)
  id_clock_pass = cpu_clock_id(' Ice shelf halo updates', grain=CLOCK_ROUTINE)

end subroutine initialize_ice_shelf

subroutine initialize_shelf_mass(G, param_file, CS, new_sim)

  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  type(ice_shelf_CS),    pointer    :: CS
  logical, optional            :: new_sim

  integer :: i, j, is, ie, js, je
  logical :: read_shelf_area, new_sim_2
  character(len=240) :: config, inputdir, shelf_file, filename
  character(len=120) :: shelf_mass_var  ! The name of shelf mass in the file.
  character(len=120) :: shelf_area_var ! The name of shelf area in the file.
  character(len=40)  :: mod = "MOM_ice_shelf"
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.not. present(new_sim)) then
    new_sim_2 = .true.
  else
    new_sim_2 = .false.
  endif

  call get_param(param_file, mod, "ICE_SHELF_CONFIG", config, &
                 "A string that specifies how the ice shelf is \n"//&
                 "initialized. Valid options include:\n"//&
                 " \tfile\t Read from a file.\n"//&
                 " \tzero\t Set shelf mass to 0 everywhere.\n"//&
                 " \tUSER\t Call USER_initialize_shelf_mass.\n", &
                 fail_if_missing=.true.)

  ! ### THIS NEEDS TO BE WRITTEN MORE COMPLETELY, WITH ADDITIONAL CASES?

  select case ( trim(config) )
    case ("file")
      call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
      inputdir = slasher(inputdir)

      call get_param(param_file, mod, "SHELF_FILE", shelf_file, &
                 "The file from which to read the shelf mass.", &
                 default="shelf_mass.nc")
      call get_param(param_file, mod, "SHELF_MASS_VAR", shelf_mass_var, &
                 "The variable in SHELF_FILE with the shelf mass.", &
                 default="shelf_mass")
      call get_param(param_file, mod, "READ_SHELF_AREA", read_shelf_area, &
                 "If true, also read the area covered by ice-shelf from SHELF_FILE.", &
                 default=.true.)
      if (read_shelf_area) &
        call get_param(param_file, mod, "SHELF_AREA_VAR", shelf_area_var, &
                 "The variable in SHELF_FILE with the shelf area.", &
                 default="shelf_area")

      filename = trim(inputdir)//trim(shelf_file)
      call log_param(param_file, mod, "INPUTDIR/SHELF_FILE", filename)
      if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
           " initialize_shelf_mass: Unable to open "//trim(filename))

      call read_data(filename, trim(shelf_mass_var), CS%mass_shelf, &
                     domain=G%Domain%mpp_domain)
      if (read_shelf_area) then
        call read_data(filename, trim(shelf_area_var), CS%area_shelf_h, &
                       domain=G%Domain%mpp_domain)
      else
        do j=js,je ; do i=is,ie
          CS%area_shelf_h(i,j) = 0.0
          if (CS%mass_shelf(i,j) > 0.0) CS%area_shelf_h(i,j) = G%areaT(i,j)
        enddo ; enddo
      endif

    case ("zero")
      do j=js,je ; do i=is,ie
        CS%mass_shelf(i,j) = 0.0
        CS%area_shelf_h(i,j) = 0.0
      enddo ; enddo

    case ("USER")
      call USER_initialize_shelf_mass(CS%mass_shelf, CS%area_shelf_h, &
               CS%h_shelf, CS%hmask, G, CS%user_CS, param_file, new_sim_2)

    case default ;  call MOM_error(FATAL,"initialize_ice_shelf: "// &
      "Unrecognized ice shelf setup "//trim(config))
  end select

end subroutine initialize_shelf_mass

subroutine update_shelf_mass(CS, Time)
  type(ice_shelf_CS),         pointer    :: CS
  type(time_type),            intent(in) :: Time

  call USER_update_shelf_mass(CS%mass_shelf, CS%area_shelf_h, CS%h_shelf, &
                              CS%hmask, CS%grid, CS%user_CS, Time, .true.)

end subroutine update_shelf_mass

subroutine initialize_diagnostic_fields (CS, FE, Time)
  type(ice_shelf_CS), pointer    :: CS
  integer             :: FE
  type(time_type),            intent(in) :: Time

  type(ocean_grid_type), pointer :: G
  integer             :: i, j, iters, isd, ied, jsd, jed
  real                 :: rhoi, rhow, OD
  type(time_type)          :: dummy_time
  real,dimension(:,:),pointer     :: OD_av, float_frac, h_shelf

  G => CS%grid
  rhoi = CS%density_ice
  rhow = CS%density_ocean_avg
  dummy_time = set_time (0,0)
  OD_av => CS%OD_av
  h_shelf => CS%h_shelf
  float_frac => CS%float_frac
  isd=G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j=jsd,jed
    do i=isd,ied
      OD = G%bathyT(i,j) - rhoi/rhow * h_shelf (i,j)
      if (OD.ge.0) then
    ! ice thickness does not take up whole ocean column -> floating
        OD_av (i,j) = OD
        float_frac(i,j) = 0.
      else
        OD_av (i,j) = 0.
        float_frac(i,j) = 1.
      endif
    enddo
  enddo

  call ice_shelf_solve_outer (CS, CS%u_shelf, CS%v_shelf, FE, iters, dummy_time)

end subroutine initialize_diagnostic_fields

subroutine ice_shelf_save_restart(CS, Time, directory, time_stamped, filename_suffix)
  type(ice_shelf_CS),         pointer    :: CS
  type(time_type),            intent(in) :: Time
  character(len=*), optional, intent(in) :: directory
  logical,          optional, intent(in) :: time_stamped
  character(len=*), optional, intent(in) :: filename_suffix

! Arguments: CS - A structure containing the internal ocean state (in).
!  (in)      Time - The model time at this call.  This is needed for mpp_write calls.
!  (in, opt) directory - An optional directory into which to write these restart files.
!  (in, opt) time_stamped - If true, the restart file names include
!                           a unique time stamp.  The default is false.
!  (in, opt) filename_suffix - An optional suffix (e.g., a time-stamp) to append
!                              to the restart file names.
  character(len=200) :: restart_dir
  character(2) :: procnum

!  write (procnum,'(I2)') mpp_pe()

  !### THESE ARE ONLY HERE FOR DEBUGGING?
! call savearray2 ("U_before_"//"p"//trim(procnum),CS%u_shelf,CS%write_output_to_file)
! call savearray2 ("V_before_"//"p"//trim(procnum),CS%v_shelf,CS%write_output_to_file)
! call savearray2 ("H_before_"//"p"//trim(procnum),CS%h_shelf,CS%write_output_to_file)
! call savearray2 ("Hmask_before_"//"p"//trim(procnum),CS%hmask,CS%write_output_to_file)
! call savearray2 ("Harea_before_"//"p"//trim(procnum),CS%area_shelf_h,CS%write_output_to_file)
! call savearray2 ("Visc_before_"//"p"//trim(procnum),CS%ice_visc_bilinear,CS%write_output_to_file)
! call savearray2 ("taub_before_"//"p"//trim(procnum),CS%taub_beta_eff_bilinear,CS%write_output_to_file)
!  call savearray2 ("taub_before_"//"p"//trim(procnum),CS%taub_beta_eff_bilinear,CS%write_output_to_file)
  if (present(directory)) then ; restart_dir = directory
  else ; restart_dir = CS%restart_output_dir ; endif

  call save_restart(restart_dir, Time, CS%grid, CS%restart_CSp, time_stamped)

end subroutine ice_shelf_save_restart


subroutine ice_shelf_advect(CS, time_step, melt_rate, Time)
  type(ice_shelf_CS),         pointer    :: CS
  real,                       intent(in) :: time_step
  real,pointer,dimension(:,:),intent(in) :: melt_rate
  type(time_type)             :: Time

! time_step: time step in sec
! melt_rate: basal melt rate in kg/m^2/s

! 3/8/11 DNG
! Arguments:
! CS - A structure containing the ice shelf state - including current velocities
! h0 - an array containing the thickness at the beginning of the call
! h_after_uflux - an array containing the thickness after advection in u-direction
! h_after_vflux - similar
!
!    This subroutine takes the velocity (on the Bgrid) and timesteps h_t = - div (uh) once.
!    ADDITIONALLY, it will update the volume of ice in partially-filled cells, and update
!        hmask accordingly
!
!    The flux overflows are included here. That is because they will be used to advect 3D scalars
!    into partial cells

  !
  ! flux_enter: this is to capture flow into partially covered cells; it gives the mass flux into a given
  ! cell across its boundaries.
  ! ###Perhaps flux_enter should be changed into u-face and v-face
  ! ###fluxes, which can then be used in halo updates, etc.
  !
  !   from left neighbor:   flux_enter (:,:,1)
  !   from right neighbor:  flux_enter (:,:,2)
  !   from bottom neighbor: flux_enter (:,:,3)
  !   from top neighbor:    flux_enter (:,:,4)
  !
  !  THESE ARE NOT CONSISTENT ==> FIND OUT WHAT YOU IMPLEMENTED

  ! flux_enter(isd:ied,jsd:jed,1:4): if cell is not ice-covered, gives flux of ice into cell from kth boundary
  !
  !   o--- (4) ---o
  !   |           |
  !  (1)         (2)
  !   |           |
  !   o--- (3) ---o
  !

  type(ocean_grid_type), pointer :: G
  real, dimension(size(CS%h_shelf,1),size(CS%h_shelf,2))   :: h_after_uflux, h_after_vflux
  real, dimension(size(CS%h_shelf,1),size(CS%h_shelf,2),4) :: flux_enter
  integer                           :: isd, ied, jsd, jed, i, j, isc, iec, jsc, jec
  real                              :: rho, spy, thick_bd
  real, dimension(:,:), pointer     :: hmask
  character(len=2)                  :: procnum

  hmask => CS%hmask
  G => CS%grid
  rho = CS%density_ice
  spy = 365 * 86400 ! seconds per year; is there a global constant for this?  No - it is dependent upon a calendar.

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  flux_enter (:,:,:) = 0.0

  h_after_uflux (:,:) = 0.0
  h_after_vflux (:,:) = 0.0
!   if (is_root_pe()) write(*,*) "ice_shelf_advect called"

  do j=jsd,jed
    do i=isd,ied
      thick_bd = CS%thickness_boundary_values(i,j)
      if (thick_bd .ne. 0.0) then
          CS%h_shelf(i,j) = CS%thickness_boundary_values(i,j)
      endif
    enddo
  enddo

  call ice_shelf_advect_thickness_x (CS, time_step/spy, CS%h_shelf, h_after_uflux, flux_enter)

!  call enable_averaging(time_step,Time,CS%diag)
 ! call pass_var (h_after_uflux, G%domain)
!  if (CS%id_h_after_uflux > 0) call post_data(CS%id_h_after_uflux, h_after_uflux, CS%diag)
!  call disable_averaging(CS%diag)

  call ice_shelf_advect_thickness_y (CS, time_step/spy, h_after_uflux, h_after_vflux, flux_enter)

!  call enable_averaging(time_step,Time,CS%diag)
!  call pass_var (h_after_vflux, G%domain)
!  if (CS%id_h_after_vflux > 0) call post_data(CS%id_h_after_vflux, h_after_vflux, CS%diag)
!  call disable_averaging(CS%diag)

  do j=jsd,jed
    do i=isd,ied
      if (CS%hmask(i,j) .eq. 1) then
        CS%h_shelf (i,j) = h_after_vflux(i,j)
      endif
    enddo
  enddo

  if (CS%moving_shelf_front) then
    call shelf_advance_front (CS, flux_enter)
    if (CS%min_thickness_simple_calve > 0.0) then
      call ice_shelf_min_thickness_calve (CS, CS%h_shelf, CS%area_shelf_h, CS%hmask)
    endif
    if (CS%calve_to_mask) then
      call calve_to_mask (CS, CS%h_shelf, CS%area_shelf_h, CS%hmask, CS%calve_mask)
    endif
  endif

  !call enable_averaging(time_step,Time,CS%diag)
  !if (CS%id_h_after_adv > 0) call post_data(CS%id_h_after_adv, CS%h_shelf, CS%diag)
  !call disable_averaging(CS%diag)

  do j=jsc,jec
    do i=isc,iec
      if ((CS%hmask(i,j) .eq. 1) .or. (CS%hmask(i,j) .eq. 2)) then
        if (melt_rate (i,j) / rho * time_step .lt. CS%h_shelf (i,j)) then
          CS%h_shelf (i,j) = CS%h_shelf (i,j) - melt_rate (i,j) / rho * time_step
        else
          ! the ice is about to melt away
          ! in this case set thickness, area, and mask to zero
          ! NOTE: not mass conservative
          ! should maybe scale salt & heat flux for this cell

          CS%h_shelf(i,j) = 0.0
          CS%hmask(i,j) = 0.0
          CS%area_shelf_h(i,j) = 0.0
        endif
      endif
    enddo
  enddo

  call update_velocity_masks (CS)

  call pass_var(CS%area_shelf_h, G%domain)
  call pass_var(CS%h_shelf, G%domain)
  call pass_var(CS%hmask, G%domain)

  if (CS%DEBUG) then
    call hchksum (CS%h_shelf, "h after front", G%HI, haloshift=3)
    call hchksum (CS%h_shelf, "shelf area after front", G%HI, haloshift=3)
  endif


end subroutine ice_shelf_advect

subroutine ice_shelf_solve_outer (CS, u, v, FE, iters, time)
  type(ice_shelf_CS),                     pointer       :: CS
  real, dimension(NILIMB_SYM_,NJLIMB_SYM_), intent(inout) :: u, v
  integer,                                intent(in)    :: FE
  integer,                                intent(out)   :: iters
  type(time_type),                        intent(in)    :: time

  real, dimension(:,:), pointer :: TAUDX, TAUDY, u_prev_iterate, v_prev_iterate, &
                        u_bdry_cont, v_bdry_cont, Au, Av, err_u, err_v, &
                        geolonq, geolatq, u_last, v_last, float_cond, H_node
  type(ocean_grid_type), pointer      :: G
  integer                 :: conv_flag, i, j, k,l, iter, isym, &
                        isdq, iedq, jsdq, jedq, isd, ied, jsd, jed, isumstart, jsumstart, nodefloat, nsub
  real                     :: err_max, err_tempu, err_tempv, err_init, area, max_vel, tempu, tempv, rhoi, rhow
  real, pointer, dimension(:,:,:,:) :: Phi
  real, pointer, dimension(:,:,:,:,:,:) :: Phisub
  real, dimension (8,4)       :: Phi_temp
  real, dimension (2,2)       :: X,Y
  character(2)                :: iternum
  character(2)                :: procnum, numproc

  ! for GL interpolation - need to make this a readable parameter
  nsub = CS%n_sub_regularize

  G => CS%grid
  isdq = G%isdB ; iedq = G%iedB ; jsdq = G%jsdB ; jedq = G%jedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  rhoi = CS%density_ice
  rhow = CS%density_ocean_avg
  ALLOCATE (TAUDX (isdq:iedq,jsdq:jedq) ) ; TAUDX(:,:)=0
  ALLOCATE (TAUDY (isdq:iedq,jsdq:jedq) ) ; TAUDY(:,:)=0
  ALLOCATE (u_prev_iterate (isdq:iedq,jsdq:jedq) )
  ALLOCATE (v_prev_iterate (isdq:iedq,jsdq:jedq) )
  ALLOCATE (u_bdry_cont (isdq:iedq,jsdq:jedq) ) ; u_bdry_cont(:,:)=0
  ALLOCATE (v_bdry_cont (isdq:iedq,jsdq:jedq) ) ; v_bdry_cont(:,:)=0
  ALLOCATE (Au (isdq:iedq,jsdq:jedq) ) ; Au(:,:)=0
  ALLOCATE (Av (isdq:iedq,jsdq:jedq) ) ; Av(:,:)=0
  ALLOCATE (err_u (isdq:iedq,jsdq:jedq) )
  ALLOCATE (err_v (isdq:iedq,jsdq:jedq) )
  ALLOCATE (u_last (isdq:iedq,jsdq:jedq) )
  ALLOCATE (v_last (isdq:iedq,jsdq:jedq) )

  ! need to make these conditional on GL interpolation
  ALLOCATE (float_cond (G%isd:G%ied,G%jsd:G%jed)) ; float_cond(:,:)=0
    ALLOCATE (H_node (G%isdB:G%iedB,G%jsdB:G%jedB)) ; H_node(:,:)=0
    ALLOCATE (Phisub (nsub,nsub,2,2,2,2)) ; Phisub = 0.0

  geolonq => G%geoLonBu ; geolatq => G%geoLatBu

  if (G%isc+G%idg_offset==G%isg) then
  ! tile is at west bdry
    isumstart = G%iscB
  else
  ! tile is interior
    isumstart = ISUMSTART_INT_
  endif

  if (G%jsc+G%jdg_offset==G%jsg) then
  ! tile is at south bdry
    jsumstart = G%jscB
  else
  ! tile is interior
    jsumstart = JSUMSTART_INT_
  endif

  call calc_shelf_driving_stress (CS, TAUDX, TAUDY, CS%OD_av, FE)

  ! this is to determine which cells contain the grounding line,
  !  the criterion being that the cell is ice-covered, with some nodes
  !  floating and some grounded
  ! floatation condition is estimated by assuming topography is cellwise constant
  !  and H is bilinear in a cell; floating where rho_i/rho_w * H_node + D is nonpositive

  ! need to make this conditional on GL interp

  if (CS%GL_regularize) then

    call interpolate_H_to_B (CS, CS%h_shelf, CS%hmask, H_node)
    call savearray2 ("H_node",H_node,CS%write_output_to_file)

    do j=G%jsc,G%jec
      do i=G%isc,G%iec
        nodefloat = 0
        do k=0,1
          do l=0,1
            if ((CS%hmask(i,j) .eq. 1) .and. &
              (rhoi/rhow * H_node(i-1+k,j-1+l) - G%bathyT(i,j) .le. 0)) then
              nodefloat = nodefloat + 1
            endif
          enddo
        enddo
        if ((nodefloat .gt. 0) .and. (nodefloat .lt. 4)) then
          !print *,"nodefloat",nodefloat
          float_cond (i,j) = 1.0
          CS%float_frac (i,j) = 1.0
        endif
      enddo
    enddo
    call savearray2 ("float_cond",float_cond,CS%write_output_to_file)

    call pass_var (float_cond, G%Domain)

    call bilinear_shape_functions_subgrid (Phisub, nsub)

    call savearray2("Phisub1111",Phisub(:,:,1,1,1,1),CS%write_output_to_file)

  endif

  ! make above conditional

  u_prev_iterate (:,:) = u(:,:)
  v_prev_iterate (:,:) = v(:,:)

  isym=0

  ! must prepare phi
  if (FE .eq. 1) then
    allocate (Phi (isd:ied,jsd:jed,1:8,1:4)) ; Phi(:,:,:,:)=0

    do j=jsd,jed
      do i=isd,ied

        if (((i .gt. isd) .and. (j .gt. jsd)) .or. (isym .eq. 1)) then
          X(:,:) = geolonq (i-1:i,j-1:j)*1000
          Y(:,:) = geolatq (i-1:i,j-1:j)*1000
        else
          X(2,:) = geolonq(i,j)*1000
          X(1,:) = geolonq(i,j)*1000-G%dxT(i,j)
          Y(:,2) = geolatq(i,j)*1000
          Y(:,1) = geolatq(i,j)*1000-G%dyT(i,j)
        endif

        call bilinear_shape_functions (X, Y, Phi_temp, area)
        Phi (i,j,:,:) = Phi_temp

      enddo
    enddo
  endif

  if (FE .eq. 1) then
      call calc_shelf_visc_bilinear (CS, u, v)

      call pass_var (CS%ice_visc_bilinear, G%domain)
      call pass_var (CS%taub_beta_eff_bilinear, G%domain)
  else
      call calc_shelf_visc_triangular (CS,u,v)

      call pass_var (CS%ice_visc_upper_tri, G%domain)
      call pass_var (CS%taub_beta_eff_upper_tri, G%domain)
      call pass_var (CS%ice_visc_lower_tri, G%domain)
      call pass_var (CS%taub_beta_eff_lower_tri, G%domain)
  endif

  ! makes sure basal stress is only applied when it is supposed to be

  do j=G%jsd,G%jed
    do i=G%isd,G%ied
      if (FE .eq. 1) then
        CS%taub_beta_eff_bilinear (i,j) = CS%taub_beta_eff_bilinear (i,j) * CS%float_frac (i,j)
      else
        CS%taub_beta_eff_upper_tri (i,j) = CS%taub_beta_eff_upper_tri (i,j) * CS%float_frac (i,j)
        CS%taub_beta_eff_lower_tri (i,j) = CS%taub_beta_eff_lower_tri (i,j) * CS%float_frac (i,j)
      endif
    enddo
  enddo

  if (FE .eq. 1) then
    call apply_boundary_values_bilinear (CS, time, Phisub, H_node, float_cond, &
      rhoi/rhow, u_bdry_cont, v_bdry_cont)
  elseif (FE .eq. 2) then
    call apply_boundary_values_triangle (CS, time, u_bdry_cont, v_bdry_cont)
  endif

  Au(:,:) = 0.0 ; Av(:,:) = 0.0

  if (FE .eq. 1) then
    call CG_action_bilinear (Au, Av, u, v, Phi, Phisub, CS%umask, CS%vmask, CS%hmask, H_node, &
              CS%ice_visc_bilinear, float_cond, G%bathyT, CS%taub_beta_eff_bilinear, G%areaT, &
              G%isc-1, G%iec+1, G%jsc-1, G%jec+1, rhoi/rhow)
  elseif (FE .eq. 2) then
    call CG_action_triangular (Au, Av, u, v, CS%umask, CS%vmask, CS%hmask, CS%ice_visc_upper_tri, &
              CS%ice_visc_lower_tri, CS%taub_beta_eff_upper_tri, CS%taub_beta_eff_lower_tri, &
              G%dxT, G%dyT, G%areaT, G%isc-1, G%iec+1, G%jsc-1, G%jec+1, isym)
  endif

!  write (procnum,'(I2)') mpp_pe()


  err_init = 0 ; err_tempu = 0; err_tempv = 0
  do j=jsumstart,G%jecB
    do i=isumstart,G%iecB
      if (CS%umask(i,j) .eq. 1) then
        err_tempu = ABS (Au(i,j) + u_bdry_cont(i,j) - TAUDX(i,j))
      endif
      if (CS%vmask(i,j) .eq. 1) then
        err_tempv = MAX(ABS (Av(i,j) + v_bdry_cont(i,j) - TAUDY(i,j)), err_tempu)
      endif
      if (err_tempv .ge. err_init) then
        err_init = err_tempv
      endif
    enddo
  enddo

  call mpp_max (err_init)

  if (is_root_pe()) print *,"INITIAL nonlinear residual: ",err_init

  u_last(:,:) = u(:,:) ; v_last(:,:) = v(:,:)

  !! begin loop

  do iter=1,100


    call ice_shelf_solve_inner (CS, u, v, TAUDX, TAUDY, H_node, float_cond, &
                                FE, conv_flag, iters, time, Phi, Phisub)


    if (CS%DEBUG) then
      call qchksum (u, "u shelf", G%HI, haloshift=2)
      call qchksum (v, "v shelf", G%HI, haloshift=2)
    endif

    if (is_root_pe()) print *,"linear solve done",iters," iterations"

    if (FE .eq. 1) then
      call calc_shelf_visc_bilinear (CS,u,v)
      call pass_var (CS%ice_visc_bilinear, G%domain)
      call pass_var (CS%taub_beta_eff_bilinear, G%domain)
    else
      call calc_shelf_visc_triangular (CS,u,v)
      call pass_var (CS%ice_visc_upper_tri, G%domain)
      call pass_var (CS%taub_beta_eff_upper_tri, G%domain)
      call pass_var (CS%ice_visc_lower_tri, G%domain)
      call pass_var (CS%taub_beta_eff_lower_tri, G%domain)
    endif

    if (iter .eq. 1) then
!      call savearray2 ("visc1",CS%ice_visc_bilinear,CS%write_output_to_file)
    endif

    ! makes sure basal stress is only applied when it is supposed to be

    do j=G%jsd,G%jed
      do i=G%isd,G%ied
        if (FE .eq. 1) then
          CS%taub_beta_eff_bilinear (i,j) = CS%taub_beta_eff_bilinear (i,j) * CS%float_frac (i,j)
        else
          CS%taub_beta_eff_upper_tri (i,j) = CS%taub_beta_eff_upper_tri (i,j) * CS%float_frac (i,j)
          CS%taub_beta_eff_lower_tri (i,j) = CS%taub_beta_eff_lower_tri (i,j) * CS%float_frac (i,j)
        endif
      enddo
    enddo

    u_bdry_cont (:,:) = 0 ; v_bdry_cont (:,:) = 0

    if (FE .eq. 1) then
      call apply_boundary_values_bilinear (CS, time, Phisub, H_node, float_cond, &
        rhoi/rhow, u_bdry_cont, v_bdry_cont)
    elseif (FE .eq. 2) then
      call apply_boundary_values_triangle (CS, time, u_bdry_cont, v_bdry_cont)
    endif

    Au(:,:) = 0 ; Av(:,:) = 0

    if (FE .eq. 1) then
      call CG_action_bilinear (Au, Av, u, v, Phi, Phisub, CS%umask, CS%vmask, CS%hmask, H_node, &
              CS%ice_visc_bilinear, float_cond, G%bathyT, CS%taub_beta_eff_bilinear, G%areaT, G%isc-1, &
              G%iec+1, G%jsc-1, G%jec+1, rhoi/rhow)
    elseif (FE .eq. 2) then
      call CG_action_triangular (Au, Av, u, v, CS%umask, CS%vmask, CS%hmask, CS%ice_visc_upper_tri, &
              CS%ice_visc_lower_tri, CS%taub_beta_eff_upper_tri, CS%taub_beta_eff_lower_tri, &
              G%dxT, G%dyT, G%areaT, G%isc-1, G%iec+1, G%jsc-1, G%jec+1, isym)
    endif

    err_max = 0

      if (CS%nonlin_solve_err_mode .eq. 1) then

      do j=jsumstart,G%jecB
        do i=isumstart,G%iecB
          if (CS%umask(i,j) .eq. 1) then
            err_tempu = ABS (Au(i,j) + u_bdry_cont(i,j) - TAUDX(i,j))
          endif
          if (CS%vmask(i,j) .eq. 1) then
            err_tempv = MAX(ABS (Av(i,j) + v_bdry_cont(i,j) - TAUDY(i,j)), err_tempu)
          endif
          if (err_tempv .ge. err_max) then
            err_max = err_tempv
          endif
        enddo
      enddo

      call mpp_max (err_max)

    elseif (CS%nonlin_solve_err_mode .eq. 2) then

      max_vel = 0 ; tempu = 0 ; tempv = 0

      do j=jsumstart,G%jecB
        do i=isumstart,G%iecB
          if (CS%umask(i,j) .eq. 1) then
            err_tempu = ABS (u_last(i,j)-u(i,j))
            tempu = u(i,j)
          endif
          if (CS%vmask(i,j) .eq. 1) then
            err_tempv = MAX(ABS (v_last(i,j)- v(i,j)), err_tempu)
            tempv = SQRT(v(i,j)**2+tempu**2)
          endif
          if (err_tempv .ge. err_max) then
            err_max = err_tempv
          endif
          if  (tempv .ge. max_vel) then
            max_vel = tempv
          endif
        enddo
      enddo

      u_last (:,:) = u(:,:)
      v_last (:,:) = v(:,:)

      call mpp_max (max_vel)
      call mpp_max (err_max)
      err_init = max_vel

    endif

    if (is_root_pe()) print *,"nonlinear residual: ",err_max/err_init

    if (err_max .le. CS%nonlinear_tolerance * err_init) then
      if (is_root_pe()) &
        print *,"exiting nonlinear solve after ",iter," iterations"
      exit
    endif

  enddo

  !write (procnum,'(I1)') mpp_pe()
  !write (numproc,'(I1)') mpp_npes()

  DEALLOCATE (TAUDX)
  DEALLOCATE (TAUDY)
  DEALLOCATE (u_prev_iterate)
  DEALLOCATE (v_prev_iterate)
  DEALLOCATE (u_bdry_cont)
  DEALLOCATE (v_bdry_cont)
  DEALLOCATE (Au)
  DEALLOCATE (Av)
  DEALLOCATE (err_u)
  DEALLOCATE (err_v)
  DEALLOCATE (u_last)
  DEALLOCATE (v_last)
  DEALLOCATE (H_node)
  DEALLOCATE (float_cond)
  DEALLOCATE (Phisub)

end subroutine ice_shelf_solve_outer

subroutine ice_shelf_solve_inner (CS, u, v, taudx, taudy, H_node, float_cond, FE, conv_flag, iters, time, Phi, Phisub)
  type(ice_shelf_CS),         pointer    :: CS
  real, dimension(NILIMB_SYM_,NJLIMB_SYM_), intent(inout)  :: u, v
  real, dimension(NILIMB_SYM_,NJLIMB_SYM_), intent(in)     :: taudx, taudy, H_node
  real, dimension(:,:),intent(in)                        :: float_cond
  integer, intent(in)          :: FE
  integer, intent(out)         :: conv_flag, iters
  type(time_type)              :: time
  real, pointer, dimension(:,:,:,:)      :: Phi
  real, dimension (:,:,:,:,:,:),pointer :: Phisub

! one linear solve (nonlinear iteration) of the solution for velocity

! in this subroutine:
!    boundary contributions are added to taud to get the RHS
!    diagonal of matrix is found (for Jacobi precondition)
!    CG iteration is carried out for max. iterations or until convergence

! assumed - u, v, taud, visc, beta_eff are valid on the halo


  real, dimension(:,:), pointer :: hmask, umask, vmask, u_bdry, v_bdry, &
                        visc, visc_lo, beta, beta_lo, geolonq, geolatq
  real, dimension(LBOUND(u,1):UBOUND(u,1),LBOUND(u,2):UBOUND(u,2)) ::  &
                        Ru, Rv, Zu, Zv, DIAGu, DIAGv, RHSu, RHSv, &
                        ubd, vbd, Au, Av, Du, Dv, &
                        Zu_old, Zv_old, Ru_old, Rv_old, &
                        sum_vec, sum_vec_2
  integer                 :: iter, i, j, isym, isd, ied, jsd, jed, &
                        isc, iec, jsc, jec, is, js, ie, je, isumstart, jsumstart, &
                        isdq, iedq, jsdq, jedq, iscq, iecq, jscq, jecq, nx_halo, ny_halo
  real                     :: tol, beta_k, alpha_k, area, dot_p1, dot_p2, resid0, cg_halo, dot_p1a, dot_p2a
  type(ocean_grid_type), pointer     :: G
  character(1)                       :: procnum
  character(2)                       :: gridsize

  real, dimension (8,4)              :: Phi_temp
  real, dimension (2,2)              :: X,Y

  hmask => CS%hmask
  umask => CS%umask
  vmask => CS%vmask
  u_bdry => CS%u_boundary_values
  v_bdry => CS%v_boundary_values

  G => CS%grid
  geolonq => G%geoLonBu
  geolatq => G%geoLatBu
  hmask => CS%hmask
  isdq = G%isdB ; iedq = G%iedB ; jsdq = G%jsdB ; jedq = G%jedB
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  ny_halo = G%domain%njhalo ; nx_halo = G%domain%nihalo
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  Zu(:,:) = 0 ; Zv(:,:) = 0 ; DIAGu(:,:) = 0 ; DIAGv(:,:) = 0
  Ru(:,:) = 0 ; Rv (:,:) = 0 ; Au (:,:) = 0 ; Av (:,:) = 0
  Du(:,:) = 0 ; Dv (:,:) = 0 ; ubd(:,:) = 0 ; vbd(:,:) = 0
  dot_p1 = 0 ; dot_p2 = 0

!   if (G%symmetric) then
!     isym = 1
!   else
!     isym = 0
!   endif

  isym = 0

  if (G%isc+G%idg_offset==G%isg) then
  ! tile is at west bdry
    isumstart = G%iscB
  else
  ! tile is interior
    isumstart = ISUMSTART_INT_
  endif

  if (G%jsc+G%jdg_offset==G%jsg) then
  ! tile is at south bdry
    jsumstart = G%jscB
  else
  ! tile is interior
    jsumstart = JSUMSTART_INT_
  endif

  if (FE .eq. 1) then
    visc => CS%ice_visc_bilinear
    beta => CS%taub_beta_eff_bilinear
  elseif (FE .eq. 2) then
    visc => CS%ice_visc_upper_tri
    visc_lo => CS%ice_visc_lower_tri
    beta => CS%taub_beta_eff_upper_tri
    beta_lo => CS%taub_beta_eff_lower_tri
  endif

  if (FE .eq. 1) then
    call apply_boundary_values_bilinear (CS, time, Phisub, H_node, float_cond, &
      CS%density_ice/CS%density_ocean_avg, ubd, vbd)
  elseif (FE .eq. 2) then
    call apply_boundary_values_triangle (CS, time, ubd, vbd)
  endif

  RHSu(:,:) = taudx(:,:) - ubd(:,:)
  RHSv(:,:) = taudy(:,:) - vbd(:,:)


  call pass_vector(RHSu, RHSv, G%domain, TO_ALL, BGRID_NE)


  if (FE .eq. 1) then
    call matrix_diagonal_bilinear(CS, float_cond, H_node, &
      CS%density_ice/CS%density_ocean_avg, Phisub, DIAGu, DIAGv)
!    DIAGu(:,:) = 1 ; DIAGv(:,:) = 1
  elseif (FE .eq. 2) then
    call matrix_diagonal_triangle (CS, DIAGu, DIAGv)
    DIAGu(:,:) = 1 ; DIAGv(:,:) = 1
  endif

  call pass_vector(DIAGu, DIAGv, G%domain, TO_ALL, BGRID_NE)



  if (FE .eq. 1) then
    call CG_action_bilinear (Au, Av, u, v, Phi, Phisub, umask, vmask, hmask, &
            H_node, visc, float_cond, G%bathyT, beta, G%areaT, isc-1, iec+1, jsc-1, &
            jec+1, CS%density_ice/CS%density_ocean_avg)
  elseif (FE .eq. 2) then
    call CG_action_triangular (Au, Av, u, v, umask, vmask, hmask, visc, visc_lo, &
            beta, beta_lo, G%dxT, G%dyT, G%areaT, isc-1, iec+1, jsc-1, jec+1, isym)
  endif

  call pass_vector(Au, Av, G%domain, TO_ALL, BGRID_NE)

  Ru(:,:) = RHSu(:,:) - Au(:,:) ; Rv(:,:) = RHSv(:,:) - Av(:,:)

  if (.not. CS%use_reproducing_sums) then

    do j=jsumstart,jecq
      do i=isumstart,iecq
        if (umask(i,j) .eq. 1) dot_p1 = dot_p1 + Ru(i,j)**2
        if (vmask(i,j) .eq. 1) dot_p1 = dot_p1 + Rv(i,j)**2
      enddo
    enddo

    call mpp_sum (dot_p1)

  else

    sum_vec(:,:) = 0.0

    do j=JSUMSTART_INT_,jecq
      do i=ISUMSTART_INT_,iecq
        if (umask(i,j) .eq. 1) sum_vec(i,j) = Ru(i,j)**2
        if (vmask(i,j) .eq. 1) sum_vec(i,j) = sum_vec(i,j) + Rv(i,j)**2
      enddo
    enddo

    dot_p1 = reproducing_sum ( sum_vec, ISUMSTART_INT_, iecq, &
                                        JSUMSTART_INT_, jecq )

  endif

  resid0 = sqrt (dot_p1)

  do j=jsdq,jedq
    do i=isdq,iedq
      if (umask(i,j) .eq. 1) Zu(i,j) = Ru (i,j) / DIAGu (i,j)
      if (vmask(i,j) .eq. 1) Zv(i,j) = Rv (i,j) / DIAGv (i,j)
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

    if (FE .eq. 1) then

      call CG_action_bilinear (Au, Av, Du, Dv, Phi, Phisub, umask, vmask, hmask, &
            H_node, visc, float_cond, G%bathyT, beta, G%areaT, is, ie, js, &
            je, CS%density_ice/CS%density_ocean_avg)

    elseif (FE .eq. 2) then

      call CG_action_triangular (Au, Av, Du, Dv, umask, vmask, hmask, visc, visc_lo, &
                beta, beta_lo, G%dxT, G%dyT, G%areaT, is, ie, js, je, isym)
    endif


    ! Au, Av valid region moves in by 1

    if ( .not. CS%use_reproducing_sums) then


      ! alpha_k = (Z \dot R) / (D \dot AD}
      dot_p1 = 0 ; dot_p2 = 0
      do j=jsumstart,jecq
        do i=isumstart,iecq
          if (umask(i,j) .eq. 1) then
            dot_p1 = dot_p1 + Zu(i,j)*Ru(i,j)
            dot_p2 = dot_p2 + Du(i,j)*Au(i,j)
          endif
          if (vmask(i,j) .eq. 1) then
              dot_p1 = dot_p1 + Zv(i,j)*Rv(i,j)
              dot_p2 = dot_p2 + Dv(i,j)*Av(i,j)
          endif
        enddo
      enddo
      call mpp_sum (dot_p1) ; call mpp_sum (dot_p2)
    else

      sum_vec(:,:) = 0.0 ; sum_vec_2(:,:) = 0.0

      do j=jscq,jecq
        do i=iscq,iecq
          if (umask(i,j) .eq. 1) sum_vec(i,j) = Zu(i,j) * Ru(i,j)
          if (vmask(i,j) .eq. 1) sum_vec(i,j) = sum_vec(i,j) + &
                                                Zv(i,j) * Rv(i,j)

          if (umask(i,j) .eq. 1) sum_vec_2(i,j) = Du(i,j) * Au(i,j)
          if (vmask(i,j) .eq. 1) sum_vec_2(i,j) = sum_vec_2(i,j) + &
                                                Dv(i,j) * Av(i,j)
        enddo
      enddo

      dot_p1 = reproducing_sum ( sum_vec, iscq, iecq, &
                                          jscq, jecq )

      dot_p2 = reproducing_sum ( sum_vec_2, iscq, iecq, &
                                          jscq, jecq )

    endif

    alpha_k = dot_p1/dot_p2

    !### These should probably use explicit index notation so that they are
    !### not applied outside of the valid range.  - RWH

    ! u(:,:) = u(:,:) + alpha_k * Du(:,:)
    ! v(:,:) = v(:,:) + alpha_k * Dv(:,:)

    do j=jsd,jed
      do i=isd,ied
        if (umask(i,j) .eq. 1) u(i,j) = u(i,j) + alpha_k * Du(i,j)
        if (vmask(i,j) .eq. 1) v(i,j) = v(i,j) + alpha_k * Dv(i,j)
      enddo
    enddo

    do j=jsd,jed
      do i=isd,ied
        if (umask(i,j) .eq. 1) then
          Ru_old(i,j) = Ru(i,j) ; Zu_old(i,j) = Zu(i,j)
        endif
        if (vmask(i,j) .eq. 1) then
          Rv_old(i,j) = Rv(i,j) ; Zv_old(i,j) = Zv(i,j)
        endif
      enddo
    enddo

!    Ru(:,:) = Ru(:,:) - alpha_k * Au(:,:)
!    Rv(:,:) = Rv(:,:) - alpha_k * Av(:,:)

    do j=jsd,jed
      do i=isd,ied
        if (umask(i,j) .eq. 1) Ru(i,j) = Ru(i,j) - alpha_k * Au(i,j)
        if (vmask(i,j) .eq. 1) Rv(i,j) = Rv(i,j) - alpha_k * Av(i,j)
      enddo
    enddo


    do j=jsdq,jedq
      do i=isdq,iedq
        if (umask(i,j) .eq. 1) then
          Zu(i,j) = Ru (i,j) / DIAGu (i,j)
        endif
        if (vmask(i,j) .eq. 1) then
          Zv(i,j) = Rv (i,j) / DIAGv (i,j)
        endif
      enddo
    enddo

    ! R,u,v,Z valid region moves in by 1

    if (.not. CS%use_reproducing_sums) then

    ! beta_k = (Z \dot R) / (Zold \dot Rold}
      dot_p1 = 0 ; dot_p2 = 0
      do j=jsumstart,jecq
        do i=isumstart,iecq
          if (umask(i,j) .eq. 1) then
            dot_p1 = dot_p1 + Zu(i,j)*Ru(i,j)
            dot_p2 = dot_p2 + Zu_old(i,j)*Ru_old(i,j)
          endif
          if (vmask(i,j) .eq. 1) then
            dot_p1 = dot_p1 + Zv(i,j)*Rv(i,j)
            dot_p2 = dot_p2 + Zv_old(i,j)*Rv_old(i,j)
          endif
        enddo
      enddo
      call mpp_sum (dot_p1) ; call mpp_sum (dot_p2)


    else

      sum_vec(:,:) = 0.0 ; sum_vec_2(:,:) = 0.0

      do j=JSUMSTART_INT_,jecq
        do i=ISUMSTART_INT_,iecq
          if (umask(i,j) .eq. 1) sum_vec(i,j) = Zu(i,j) * Ru(i,j)
          if (vmask(i,j) .eq. 1) sum_vec(i,j) = sum_vec(i,j) + &
                                                Zv(i,j) * Rv(i,j)

          if (umask(i,j) .eq. 1) sum_vec_2(i,j) = Zu_old(i,j) * Ru_old(i,j)
          if (vmask(i,j) .eq. 1) sum_vec_2(i,j) = sum_vec_2(i,j) + &
                                                Zv_old(i,j) * Rv_old(i,j)
        enddo
      enddo


      dot_p1 = reproducing_sum ( sum_vec, ISUMSTART_INT_, iecq, &
                                          JSUMSTART_INT_, jecq )

      dot_p2 = reproducing_sum ( sum_vec_2, ISUMSTART_INT_, iecq, &
                                          JSUMSTART_INT_, jecq )

    endif

    beta_k = dot_p1/dot_p2


!    Du(:,:) = Zu(:,:) + beta_k * Du(:,:)
!    Dv(:,:) = Zv(:,:) + beta_k * Dv(:,:)

    do j=jsd,jed
      do i=isd,ied
        if (umask(i,j) .eq. 1) Du(i,j) = Zu(i,j) + beta_k * Du(i,j)
        if (vmask(i,j) .eq. 1) Dv(i,j) = Zv(i,j) + beta_k * Dv(i,j)
      enddo
    enddo

   ! D valid region moves in by 1

    dot_p1 = 0

    if (.not. CS%use_reproducing_sums) then

      do j=jsumstart,jecq
        do i=isumstart,iecq
          if (umask(i,j) .eq. 1) then
            dot_p1 = dot_p1 + Ru(i,j)**2
          endif
          if (vmask(i,j) .eq. 1) then
            dot_p1 = dot_p1 + Rv(i,j)**2
          endif
        enddo
      enddo
      call mpp_sum (dot_p1)

    else

      sum_vec(:,:) = 0.0

      do j=JSUMSTART_INT_,jecq
        do i=ISUMSTART_INT_,iecq
          if (umask(i,j) .eq. 1) sum_vec(i,j) = Ru(i,j)**2
          if (vmask(i,j) .eq. 1) sum_vec(i,j) = sum_vec(i,j) + Rv(i,j)**2
        enddo
      enddo

      dot_p1 = reproducing_sum ( sum_vec, ISUMSTART_INT_, iecq, &
                                        JSUMSTART_INT_, jecq )

!      if (is_root_pe()) print *, dot_p1
!      if (is_root_pe()) print *, dot_p1a

    endif

    dot_p1 = sqrt (dot_p1)

!    if (mpp_pe () == 0) then
!         print *,"|r|",dot_p1
!     endif

    if (dot_p1 .le. CS%cg_tolerance * resid0) then
      iters = iter
      conv_flag = 1
      exit
    endif

    cg_halo = cg_halo - 1

    if (cg_halo .eq. 0) then
      ! pass vectors
      call pass_vector(Du, Dv, G%domain, TO_ALL, BGRID_NE)
      call pass_vector(u, v, G%domain, TO_ALL, BGRID_NE)
      call pass_vector(Ru, Rv, G%domain, TO_ALL, BGRID_NE)
      cg_halo = 3
    endif

  enddo ! end of CG loop

  do j=jsdq,jedq
    do i=isdq,iedq
      if (umask(i,j) .eq. 3) then
        u(i,j) = u_bdry(i,j)
      elseif (umask(i,j) .eq. 0) then
        u(i,j) = 0
      endif

      if (vmask(i,j) .eq. 3) then
        v(i,j) = v_bdry(i,j)
      elseif (vmask(i,j) .eq. 0) then
        v(i,j) = 0
      endif
    enddo
  enddo

  call pass_vector (u,v, G%domain, TO_ALL, BGRID_NE)

  if (conv_flag .eq. 0) then
    iters = CS%cg_max_iterations
  endif

end subroutine ice_shelf_solve_inner

subroutine ice_shelf_advect_thickness_x (CS, time_step, h0, h_after_uflux, flux_enter)
  type(ice_shelf_CS),         pointer    :: CS
  real,                       intent(in) :: time_step
  real, dimension(:,:), intent(in) :: h0
  real, dimension(:,:), intent(inout) :: h_after_uflux
  real, dimension(:,:,:), intent(inout) :: flux_enter

  ! use will be made of CS%hmask here - its value at the boundary will be zero, just like uncovered cells

  ! if there is an input bdry condition, the thickness there will be set in initialization

  ! flux_enter(isd:ied,jsd:jed,1:4): if cell is not ice-covered, gives flux of ice into cell from kth boundary
  !
  !   from left neighbor:   flux_enter (:,:,1)
  !   from right neighbor:  flux_enter (:,:,2)
  !   from bottom neighbor: flux_enter (:,:,3)
  !   from top neighbor:    flux_enter (:,:,4)
  !
  !        o--- (4) ---o
  !        |           |
  !       (1)         (2)
  !        |           |
  !        o--- (3) ---o
  !

  integer :: isym, i, j, is, ie, js, je, isd, ied, jsd, jed, gjed, gied
  integer :: i_off, j_off
  logical :: at_east_bdry, at_west_bdry, one_off_west_bdry, one_off_east_bdry
  type(ocean_grid_type), pointer :: G
  real, dimension(-2:2) :: stencil
  real, dimension(:,:), pointer  :: hmask, u_face_mask, u_flux_boundary_values
  real :: u_face, &  ! positive if out
      flux_diff_cell, phi, dxh, dyh, dxdyh

  character (len=1)        :: debug_str, procnum

!   if (CS%grid%symmetric) then
!     isym = 1
!   else
!     isym = 0
!   endif

  isym = 0

  G => CS%grid
  hmask => CS%hmask
  u_face_mask => CS%u_face_mask
  u_flux_boundary_values => CS%u_flux_boundary_values
  is = G%isc-2 ; ie = G%iec+2 ; js = G%jsc ; je = G%jec ; isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset

  do j=jsd+1,jed-1
    if (((j+j_off) .le. G%domain%njglobal+G%domain%njhalo) .AND. &
        ((j+j_off) .ge. G%domain%njhalo+1)) then ! based on mehmet's code - only if btw north & south boundaries

      stencil(:) = -1
!     if (i+i_off .eq. G%domain%nihalo+G%domain%nihalo)
      do i=is,ie

        if (((i+i_off) .le. G%domain%niglobal+G%domain%nihalo) .AND. &
             ((i+i_off) .ge. G%domain%nihalo+1)) then

          if (i+i_off .eq. G%domain%nihalo+1) then
            at_west_bdry=.true.
          else
            at_west_bdry=.false.
          endif

          if (i+i_off .eq. G%domain%niglobal+G%domain%nihalo) then
            at_east_bdry=.true.
          else
            at_east_bdry=.false.
          endif

          if (hmask(i,j) .eq. 1) then

            dxh = G%dxT(i,j) ; dyh = G%dyT(i,j) ; dxdyh = G%areaT(i,j)

            h_after_uflux(i,j) = h0(i,j)

            stencil(:) = h0(i-2:i+2,j)  ! fine as long has nx_halo >= 2

            flux_diff_cell = 0

            ! 1ST DO LEFT FACE

            if (u_face_mask (i-1,j) .eq. 4.) then

              flux_diff_cell = flux_diff_cell + dyh * time_step * u_flux_boundary_values (i-1,j) / dxdyh

            else

              ! get u-velocity at center of left face
              u_face = 0.5 * (CS%u_shelf(i-1,j-1) + CS%u_shelf(i-1,j))

  !            if (at_west_bdry .and. (i .eq. G%isc)) then
  !                print *, j, u_face, stencil(-1)
  !            endif

              if (u_face .gt. 0) then !flux is into cell - we need info from h(i-2), h(i-1) if available

              ! i may not cover all the cases.. but i cover the realistic ones

                if (at_west_bdry .AND. (hmask(i-1,j).eq.3)) then ! at western bdry but there is a thickness bdry condition,
                              ! and the stencil contains it
                  stencil (-1) = CS%thickness_boundary_values(i-1,j)
                  flux_diff_cell = flux_diff_cell + ABS(u_face) * dyh * time_step * stencil(-1) / dxdyh

                elseif (hmask(i-1,j) * hmask(i-2,j) .eq. 1) then  ! h(i-2) and h(i-1) are valid
                  phi = slope_limiter (stencil(-1)-stencil(-2), stencil(0)-stencil(-1))
                  flux_diff_cell = flux_diff_cell + ABS(u_face) * dyh* time_step / dxdyh * &
                           (stencil(-1) - phi * (stencil(-1)-stencil(0))/2)

                else                            ! h(i-1) is valid
                                    ! (o.w. flux would most likely be out of cell)
                                    !  but h(i-2) is not

                  flux_diff_cell = flux_diff_cell + ABS(u_face) * dyh * time_step / dxdyh * stencil(-1)

                endif

              elseif (u_face .lt. 0) then !flux is out of cell - we need info from h(i-1), h(i+1) if available
                if (hmask(i-1,j) * hmask(i+1,j) .eq. 1) then         ! h(i-1) and h(i+1) are both valid
                  phi = slope_limiter (stencil(0)-stencil(1), stencil(-1)-stencil(0))
                  flux_diff_cell = flux_diff_cell - ABS(u_face) * dyh * time_step / dxdyh * &
                             (stencil(0) - phi * (stencil(0)-stencil(-1))/2)

                else
                  flux_diff_cell = flux_diff_cell - ABS(u_face) * dyh * time_step / dxdyh * stencil(0)

                  if ((hmask(i-1,j) .eq. 0) .OR. (hmask(i-1,j) .eq. 2)) then
                    flux_enter(i-1,j,2) = ABS(u_face) * dyh * time_step * stencil(0)
                  endif
                endif
              endif
            endif

            ! NEXT DO RIGHT FACE

            ! get u-velocity at center of right face

            if (u_face_mask (i+1,j) .eq. 4.) then

              flux_diff_cell = flux_diff_cell + dyh * time_step * u_flux_boundary_values (i+1,j) / dxdyh

            else

              u_face = 0.5 * (CS%u_shelf(i,j-1) + CS%u_shelf(i,j))

              if (u_face .lt. 0) then !flux is into cell - we need info from h(i+2), h(i+1) if available

                if (at_east_bdry .AND. (hmask(i+1,j).eq.3)) then ! at eastern bdry but there is a thickness bdry condition,
                                            ! and the stencil contains it

                  flux_diff_cell = flux_diff_cell + ABS(u_face) * dyh * time_step * stencil(1) / dxdyh

                elseif (hmask(i+1,j) * hmask(i+2,j) .eq. 1) then  ! h(i+2) and h(i+1) are valid

                  phi = slope_limiter (stencil(1)-stencil(2), stencil(0)-stencil(1))
                  flux_diff_cell = flux_diff_cell + ABS(u_face) * dyh * time_step / dxdyh * &
                      (stencil(1) - phi * (stencil(1)-stencil(0))/2)

                else                            ! h(i+1) is valid
                                            ! (o.w. flux would most likely be out of cell)
                                            !  but h(i+2) is not

                  flux_diff_cell = flux_diff_cell + ABS(u_face) * dyh * time_step / dxdyh * stencil(1)

                endif

              elseif (u_face .gt. 0) then !flux is out of cell - we need info from h(i-1), h(i+1) if available

                if (hmask(i-1,j) * hmask(i+1,j) .eq. 1) then         ! h(i-1) and h(i+1) are both valid

                  phi = slope_limiter (stencil(0)-stencil(-1), stencil(1)-stencil(0))
                  flux_diff_cell = flux_diff_cell - ABS(u_face) * dyh * time_step / dxdyh * &
                      (stencil(0) - phi * (stencil(0)-stencil(1))/2)

                else                            ! h(i+1) is valid
                                            ! (o.w. flux would most likely be out of cell)
                                            !  but h(i+2) is not

                  flux_diff_cell = flux_diff_cell - ABS(u_face) * dyh * time_step / dxdyh * stencil(0)

                  if ((hmask(i+1,j) .eq. 0) .OR. (hmask(i+1,j) .eq. 2)) then
                    flux_enter(i+1,j,1) = ABS(u_face) * dyh * time_step  * stencil(0)
                  endif

                endif

              endif

              h_after_uflux(i,j) = h_after_uflux(i,j) + flux_diff_cell

            endif

          elseif ((hmask(i,j) .eq. 0) .OR. (hmask(i,j) .eq. 2)) then

            if (at_west_bdry .AND. (hmask(i-1,j) .EQ. 3)) then
              u_face = 0.5 * (CS%u_shelf(i-1,j-1) + CS%u_shelf(i-1,j))
              flux_enter (i,j,1) = ABS(u_face) * G%dyT(i,j) * time_step * CS%thickness_boundary_values(i-1,j)
            elseif (u_face_mask (i-1,j) .eq. 4.) then
              flux_enter (i,j,1) = G%dyT(i,j) * time_step * u_flux_boundary_values (i-1,j)
            endif

            if (at_east_bdry .AND. (hmask(i+1,j) .EQ. 3)) then
              u_face = 0.5 * (CS%u_shelf(i,j-1) + CS%u_shelf(i,j))
              flux_enter(i,j,2) = ABS(u_face) * G%dyT(i,j) * time_step * CS%thickness_boundary_values(i+1,j)
            elseif (u_face_mask (i+1,j) .eq. 4.) then
              flux_enter (i,j,2) = G%dyT(i,j) * time_step * u_flux_boundary_values (i+1,j)
            endif

            if ((i .eq. is) .AND. (hmask(i,j) .eq. 0) .AND. (hmask(i-1,j) .eq. 1)) then
              ! this is solely for the purposes of keeping the mask consistent while advancing the front without having
              ! to call pass_var - if cell is empty and cell to left is ice-covered then this cell will become partly covered

              hmask(i,j) = 2
            elseif ((i .eq. ie) .AND. (hmask(i,j) .eq. 0) .AND. (hmask(i+1,j) .eq. 1)) then
              ! this is solely for the purposes of keeping the mask consistent while advancing the front without having
              ! to call pass_var - if cell is empty and cell to left is ice-covered then this cell will become partly covered

              hmask(i,j) = 2

            endif

          endif

        endif

      enddo ! i loop

    endif

  enddo ! j loop

!  write (procnum,'(I1)') mpp_pe()

end subroutine ice_shelf_advect_thickness_x

subroutine ice_shelf_advect_thickness_y (CS, time_step, h_after_uflux, h_after_vflux, flux_enter)
  type(ice_shelf_CS),         pointer    :: CS
  real,                       intent(in) :: time_step
  real, dimension(:,:), intent(in) :: h_after_uflux
  real, dimension(:,:), intent(inout) :: h_after_vflux
  real, dimension(:,:,:), intent(inout) :: flux_enter

  ! use will be made of CS%hmask here - its value at the boundary will be zero, just like uncovered cells

  ! if there is an input bdry condition, the thickness there will be set in initialization

  ! flux_enter(isd:ied,jsd:jed,1:4): if cell is not ice-covered, gives flux of ice into cell from kth boundary
  !
  !   from left neighbor:   flux_enter (:,:,1)
  !   from right neighbor:  flux_enter (:,:,2)
  !   from bottom neighbor: flux_enter (:,:,3)
  !   from top neighbor:    flux_enter (:,:,4)
  !
  !        o--- (4) ---o
  !        |           |
  !       (1)         (2)
  !        |           |
  !        o--- (3) ---o
  !

  integer :: isym, i, j, is, ie, js, je, isd, ied, jsd, jed, gjed, gied
  integer :: i_off, j_off
  logical :: at_north_bdry, at_south_bdry, one_off_west_bdry, one_off_east_bdry
  type(ocean_grid_type), pointer :: G
  real, dimension(-2:2) :: stencil
  real, dimension(:,:), pointer  :: hmask, v_face_mask, v_flux_boundary_values
  real :: v_face, &  ! positive if out
      flux_diff_cell, phi, dxh, dyh, dxdyh
  character(len=1)        :: debug_str, procnum

!   if (CS%grid%symmetric) then
!     isym = 1
!   else
!     isym = 0
!   endif

  isym = 0

  G => CS%grid
  hmask => CS%hmask
  v_face_mask => CS%v_face_mask
  v_flux_boundary_values => CS%v_flux_boundary_values
  is = G%isc ; ie = G%iec ; js = G%jsc-1 ; je = G%jec+1 ; isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset

  do i=isd+2,ied-2
    if (((i+i_off) .le. G%domain%niglobal+G%domain%nihalo) .AND. &
       ((i+i_off) .ge. G%domain%nihalo+1)) then  ! based on mehmet's code - only if btw east & west boundaries

      stencil(:) = -1

      do j=js,je

        if (((j+j_off) .le. G%domain%njglobal+G%domain%njhalo) .AND. &
             ((j+j_off) .ge. G%domain%njhalo+1)) then

          if (j+j_off .eq. G%domain%njhalo+1) then
            at_south_bdry=.true.
          else
            at_south_bdry=.false.
          endif

          if (j+j_off .eq. G%domain%njglobal+G%domain%njhalo) then
            at_north_bdry=.true.
          else
            at_north_bdry=.false.
          endif

          if (hmask(i,j) .eq. 1) then
            dxh = G%dxT(i,j) ; dyh = G%dyT(i,j) ; dxdyh = G%areaT(i,j)
            h_after_vflux (i,j) = h_after_uflux (i,j)

            stencil (:) = h_after_uflux (i,j-2:j+2)  ! fine as long has ny_halo >= 2
            flux_diff_cell = 0

            ! 1ST DO south FACE

            if (v_face_mask (i,j-1) .eq. 4.) then

              flux_diff_cell = flux_diff_cell + dxh * time_step * v_flux_boundary_values (i,j-1) / dxdyh

            else

              ! get u-velocity at center of left face
              v_face = 0.5 * (CS%v_shelf(i-1,j-1) + CS%v_shelf(i,j-1))

              if (v_face .gt. 0) then !flux is into cell - we need info from h(j-2), h(j-1) if available

                ! i may not cover all the cases.. but i cover the realistic ones

                if (at_south_bdry .AND. (hmask(i,j-1).eq.3)) then ! at western bdry but there is a thickness bdry condition,
                                            ! and the stencil contains it
                  flux_diff_cell = flux_diff_cell + ABS(v_face) * dxh * time_step * stencil(-1) / dxdyh

                elseif (hmask(i,j-1) * hmask(i,j-2) .eq. 1) then  ! h(j-2) and h(j-1) are valid

                  phi = slope_limiter (stencil(-1)-stencil(-2), stencil(0)-stencil(-1))
                  flux_diff_cell = flux_diff_cell + ABS(v_face) * dxh * time_step / dxdyh * &
                      (stencil(-1) - phi * (stencil(-1)-stencil(0))/2)

                else     ! h(j-1) is valid
                         ! (o.w. flux would most likely be out of cell)
                         !  but h(j-2) is not
                  flux_diff_cell = flux_diff_cell + ABS(v_face) * dxh * time_step / dxdyh * stencil(-1)
                endif

              elseif (v_face .lt. 0) then !flux is out of cell - we need info from h(j-1), h(j+1) if available

                if (hmask(i,j-1) * hmask(i,j+1) .eq. 1) then  ! h(j-1) and h(j+1) are both valid
                  phi = slope_limiter (stencil(0)-stencil(1), stencil(-1)-stencil(0))
                  flux_diff_cell = flux_diff_cell - ABS(v_face) * dxh * time_step / dxdyh * &
                      (stencil(0) - phi * (stencil(0)-stencil(-1))/2)
                else
                  flux_diff_cell = flux_diff_cell - ABS(v_face) * dxh * time_step / dxdyh * stencil(0)

                  if ((hmask(i,j-1) .eq. 0) .OR. (hmask(i,j-1) .eq. 2)) then
                    flux_enter(i,j-1,4) = ABS(v_face) * dyh * time_step * stencil(0)
                  endif

                endif

              endif

            endif

            ! NEXT DO north FACE

            if (v_face_mask(i,j+1) .eq. 4.) then

              flux_diff_cell = flux_diff_cell + dxh * time_step * v_flux_boundary_values (i,j+1) / dxdyh

            else

            ! get u-velocity at center of right face
              v_face = 0.5 * (CS%v_shelf(i-1,j) + CS%v_shelf(i,j))

              if (v_face .lt. 0) then !flux is into cell - we need info from h(j+2), h(j+1) if available

                if (at_north_bdry .AND. (hmask(i,j+1).eq.3)) then ! at eastern bdry but there is a thickness bdry condition,
                                            ! and the stencil contains it
                  flux_diff_cell = flux_diff_cell + ABS(v_face) * dxh * time_step * stencil(1) / dxdyh
                elseif (hmask(i,j+1) * hmask(i,j+2) .eq. 1) then  ! h(j+2) and h(j+1) are valid
                  phi = slope_limiter (stencil(1)-stencil(2), stencil(0)-stencil(1))
                  flux_diff_cell = flux_diff_cell + ABS(v_face) * dxh * time_step / dxdyh * &
                      (stencil(1) - phi * (stencil(1)-stencil(0))/2)
                else     ! h(j+1) is valid
                         ! (o.w. flux would most likely be out of cell)
                         !  but h(j+2) is not
                  flux_diff_cell = flux_diff_cell + ABS(v_face) * dxh * time_step / dxdyh * stencil(1)
                endif

              elseif (v_face .gt. 0) then !flux is out of cell - we need info from h(j-1), h(j+1) if available

                if (hmask(i,j-1) * hmask(i,j+1) .eq. 1) then         ! h(j-1) and h(j+1) are both valid
                  phi = slope_limiter (stencil(0)-stencil(-1), stencil(1)-stencil(0))
                  flux_diff_cell = flux_diff_cell - ABS(v_face) * dxh * time_step / dxdyh * &
                      (stencil(0) - phi * (stencil(0)-stencil(1))/2)
                else   ! h(j+1) is valid
                       ! (o.w. flux would most likely be out of cell)
                       !  but h(j+2) is not
                  flux_diff_cell = flux_diff_cell - ABS(v_face) * dxh * time_step / dxdyh * stencil(0)
                  if ((hmask(i,j+1) .eq. 0) .OR. (hmask(i,j+1) .eq. 2)) then
                    flux_enter(i,j+1,3) = ABS(v_face) * dxh * time_step * stencil(0)
                  endif
                endif

              endif

            endif

            h_after_vflux (i,j) = h_after_vflux (i,j) + flux_diff_cell

          elseif ((hmask(i,j) .eq. 0) .OR. (hmask(i,j) .eq. 2)) then

            if (at_south_bdry .AND. (hmask(i,j-1) .EQ. 3)) then
              v_face = 0.5 * (CS%u_shelf(i-1,j-1) + CS%u_shelf(i,j-1))
              flux_enter (i,j,3) = ABS(v_face) * G%dxT(i,j) * time_step * CS%thickness_boundary_values(i,j-1)
            elseif (v_face_mask(i,j-1) .eq. 4.) then
              flux_enter (i,j,3) = G%dxT(i,j) * time_step * v_flux_boundary_values (i,j-1)
            endif

            if (at_north_bdry .AND. (hmask(i,j+1) .EQ. 3)) then
              v_face = 0.5 * (CS%u_shelf(i-1,j) + CS%u_shelf(i,j))
              flux_enter (i,j,4) = ABS(v_face) * G%dxT(i,j) * time_step * CS%thickness_boundary_values(i,j+1)
            elseif (v_face_mask(i,j+1) .eq. 4.) then
              flux_enter (i,j,4) = G%dxT(i,j) * time_step * v_flux_boundary_values (i,j+1)
            endif

            if ((j .eq. js) .AND. (hmask(i,j) .eq. 0) .AND. (hmask(i,j-1) .eq. 1)) then
                ! this is solely for the purposes of keeping the mask consistent while advancing the front without having
                ! to call pass_var - if cell is empty and cell to left is ice-covered then this cell will become partly covered
              hmask (i,j) = 2
            elseif ((j .eq. je) .AND. (hmask(i,j) .eq. 0) .AND. (hmask(i,j+1) .eq. 1)) then
                ! this is solely for the purposes of keeping the mask consistent while advancing the front without having
                ! to call pass_var - if cell is empty and cell to left is ice-covered then this cell will become partly covered
              hmask (i,j) = 2
            endif

          endif
        endif
      enddo ! j loop
    endif
  enddo ! i loop

  !write (procnum,'(I1)') mpp_pe()

end subroutine ice_shelf_advect_thickness_y

subroutine shelf_advance_front (CS, flux_enter)
  type(ice_shelf_CS),         pointer    :: CS
  real, dimension(:,:,:), intent(inout)  :: flux_enter

  ! in this subroutine we go through the computational cells only and, if they are empty or partial cells,
  ! we find the reference thickness and update the shelf mass and partial area fraction and the hmask if necessary

  ! if any cells go from partial to complete, we then must set the thickness, update hmask accordingly,
  ! and divide the overflow across the adjacent EMPTY (not partly-covered) cells.
  ! (it is highly unlikely there will not be any; in which case this will need to be rethought.)

  ! most likely there will only be one "overflow". if not, though, a pass_var of all relevant variables
  ! is done; there will therefore be a loop which, in practice, will hopefully not have to go through
  ! many iterations

  ! when 3d advected scalars are introduced, they will be impacted by what is done here

  ! flux_enter(isd:ied,jsd:jed,1:4): if cell is not ice-covered, gives flux of ice into cell from kth boundary
  !
  !   from left neighbor:   flux_enter (:,:,1)
  !   from right neighbor:  flux_enter (:,:,2)
  !   from bottom neighbor: flux_enter (:,:,3)
  !   from top neighbor:    flux_enter (:,:,4)
  !
  !        o--- (4) ---o
  !        |           |
  !       (1)         (2)
  !        |           |
  !        o--- (3) ---o
  !

  integer :: i, j, isc, iec, jsc, jec, n_flux, k, l, iter_count, isym
  integer :: i_off, j_off
  integer :: iter_flag
  type(ocean_grid_type), pointer :: G
  real, dimension(:,:), pointer  :: hmask, mass_shelf, area_shelf_h, u_face_mask, v_face_mask, h_shelf
  real :: h_reference, dxh, dyh, dxdyh, rho, partial_vol, tot_flux
  integer, dimension(4) :: mapi, mapj, new_partial
!   real, dimension(size(flux_enter,1),size(flux_enter,2),size(flux_enter,2)) :: flux_enter_replace
  real, dimension (:,:,:), pointer :: flux_enter_replace => NULL()

  G => CS%grid
  h_shelf => CS%h_shelf
  hmask => CS%hmask
  mass_shelf => CS%mass_shelf
  area_shelf_h => CS%area_shelf_h
  u_face_mask => CS%u_face_mask
  v_face_mask => CS%v_face_mask
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  i_off = G%idg_offset ; j_off = G%jdg_offset
  rho = CS%density_ice
  iter_count = 0 ; iter_flag = 1

!   if (G%symmetric) then
!     isym = 1
!   else
!     isym = 0
!   endif

  isym = 0

  mapi(1) = -1 ; mapi(2) = 1 ; mapi(3:4) = 0
  mapj(3) = -1 ; mapj(4) = 1 ; mapj(1:2) = 0

  do while (iter_flag .eq. 1)

    iter_flag = 0

    if (iter_count .gt. 0) then
      flux_enter (:,:,:) = flux_enter_replace(:,:,:)
      flux_enter_replace (:,:,:) = 0.0
    endif

    iter_count = iter_count + 1

    ! if iter_count .ge. 3 then some halo updates need to be done...



    do j=jsc-1,jec+1

      if (((j+j_off) .le. G%domain%njglobal+G%domain%njhalo) .AND. &
         ((j+j_off) .ge. G%domain%njhalo+1)) then

      do i=isc-1,iec+1

         if (((i+i_off) .le. G%domain%niglobal+G%domain%nihalo) .AND. &
             ((i+i_off) .ge. G%domain%nihalo+1)) then
        ! first get reference thickness by averaging over cells that are fluxing into this cell
            n_flux = 0
            h_reference = 0.0
            tot_flux = 0.0

            do k=1,2
              if (flux_enter(i,j,k) .gt. 0) then
                n_flux = n_flux + 1
                h_reference = h_reference + h_shelf(i+2*k-3,j)
                tot_flux = tot_flux + flux_enter(i,j,k)
                flux_enter(i,j,k) = 0.0
              endif
            enddo

            do k=1,2
              if (flux_enter(i,j,k+2) .gt. 0) then
                n_flux = n_flux + 1
                h_reference = h_reference + h_shelf (i,j+2*k-3)
                tot_flux = tot_flux + flux_enter(i,j,k+2)
                flux_enter (i,j,k+2) = 0.0
              endif
            enddo

            if (n_flux .gt. 0) then
              dxdyh = G%areaT(i,j)
              h_reference = h_reference / real(n_flux)
              partial_vol = h_shelf (i,j) * area_shelf_h (i,j) + tot_flux

              if ((partial_vol / dxdyh) .eq. h_reference) then ! cell is exactly covered, no overflow
                hmask (i,j) = 1
                h_shelf (i,j) = h_reference
                area_shelf_h(i,j) = dxdyh
              elseif ((partial_vol / dxdyh) .lt. h_reference) then
                hmask (i,j) = 2
        !         mass_shelf (i,j) = partial_vol * rho
                area_shelf_h (i,j) = partial_vol / h_reference
                h_shelf (i,j) = h_reference
              else
                if (.not. associated (flux_enter_replace)) then
                  allocate ( flux_enter_replace (G%isd:G%ied,G%jsd:G%jed,1:4) )
                  flux_enter_replace (:,:,:) = 0.0
                endif

                hmask (i,j) = 1
                area_shelf_h(i,j) = dxdyh
                !h_temp (i,j) = h_reference
                partial_vol = partial_vol - h_reference * dxdyh

                iter_flag  = 1

                n_flux = 0 ; new_partial (:) = 0

                do k=1,2
                  if (u_face_mask (i-2+k,j) .eq. 2) then
                    n_flux = n_flux + 1
                  elseif (hmask (i+2*k-3,j) .eq. 0) then
                    n_flux = n_flux + 1
                    new_partial (k) = 1
                  endif
                enddo
                do k=1,2
                  if (v_face_mask (i,j-2+k) .eq. 2) then
                    n_flux = n_flux + 1
                  elseif (hmask (i,j+2*k-3) .eq. 0) then
                    n_flux = n_flux + 1
                    new_partial (k+2) = 1
                  endif
                enddo

                if (n_flux .eq. 0) then ! there is nowhere to put the extra ice!
                  h_shelf(i,j) = h_reference + partial_vol / dxdyh
                else
                  h_shelf(i,j) = h_reference

                  do k=1,2
                    if (new_partial(k) .eq. 1) &
                      flux_enter_replace (i+2*k-3,j,3-k) = partial_vol / real(n_flux)
                  enddo
                  do k=1,2 ! ### Combine these two loops?
                    if (new_partial(k+2) .eq. 1) &
                      flux_enter_replace(i,j+2*k-3,5-k) = partial_vol / real(n_flux)
                  enddo
                endif

              endif ! Parital_vol test.
            endif ! n_flux gt 0 test.

          endif
        enddo ! j-loop
      endif
    enddo

  !  call mpp_max(iter_flag)

  enddo ! End of do while(iter_flag) loop

  call mpp_max(iter_count)

  if(is_root_pe() .and. (iter_count.gt.1)) print *, iter_count, "MAX ITERATIONS,ADVANCE FRONT"

  if (associated(flux_enter_replace)) DEALLOCATE(flux_enter_replace)

end subroutine shelf_advance_front


subroutine ice_shelf_min_thickness_calve (CS, h_shelf, area_shelf_h,hmask)
  type(ice_shelf_CS), pointer       :: CS
  real, dimension(:,:), intent(inout)   :: h_shelf, area_shelf_h, hmask

  type(ocean_grid_type), pointer :: G
  integer                        :: i,j

  G => CS%grid

  do j=G%jsd,G%jed
    do i=G%isd,G%ied
      if ((h_shelf(i,j) .lt. CS%min_thickness_simple_calve) .and. (hmask(i,j).eq.1) .and. &
           (CS%float_frac(i,j) .eq. 0.0)) then
        h_shelf(i,j) = 0.0
        area_shelf_h(i,j) = 0.0
        hmask(i,j) = 0.0
      endif
    enddo
  enddo

end subroutine ice_shelf_min_thickness_calve

subroutine calve_to_mask (CS, h_shelf, area_shelf_h, hmask, calve_mask)
  type(ice_shelf_CS), pointer       :: CS
  real, dimension(:,:), intent(inout)   :: h_shelf, area_shelf_h, hmask, calve_mask

  type(ocean_grid_type), pointer :: G
  integer                        :: i,j

  G => CS%grid

  if (CS%calve_to_mask) then
    do j=G%jsc,G%jec
      do i=G%isc,G%iec
        if ((calve_mask(i,j) .eq. 0.0) .and. (hmask(i,j) .ne. 0.0)) then
          h_shelf(i,j) = 0.0
          area_shelf_h(i,j) = 0.0
          hmask(i,j) = 0.0
        endif
      enddo
    enddo
  endif

end subroutine calve_to_mask

subroutine calc_shelf_driving_stress (CS, TAUD_X, TAUD_Y, OD, FE)
  type(ice_shelf_CS),         pointer   :: CS
  real, dimension(:,:), intent(in)    :: OD
  real, dimension(NILIMB_SYM_,NJLIMB_SYM_), intent(inout)    :: TAUD_X, TAUD_Y
  integer, intent(in)            :: FE

! driving stress!

! ! TAUD_X and TAUD_Y will hold driving stress in the x- and y- directions when done.
!    they will sit on the BGrid, and so their size depends on whether the grid is symmetric
!
! Since this is a finite element solve, they will actually have the form \int \phi_i rho g h \nabla s
!
! OD -this is important and we do not yet know where (in MOM) it will come from. It represents
!     "average" ocean depth -- and is needed to find surface elevation
!    (it is assumed that base_ice = bed + OD)

! FE : 1 if bilinear, 2 if triangular linear FE

  real, dimension (:,:), pointer :: D, & ! ocean floor depth
                                    H, &  ! ice shelf thickness
                          hmask, u_face_mask, v_face_mask, float_frac
  real, dimension (SIZE(OD,1),SIZE(OD,2))  :: S, &     ! surface elevation
                            BASE     ! basal elevation of shelf/stream
  character(1)                   :: procnum


  real      :: rho, rhow, sx, sy, neumann_val, dxh, dyh, dxdyh

  type(ocean_grid_type), pointer :: G
  integer :: isym, i, j, iscq, iecq, jscq, jecq, isd, jsd, is, js, iegq, jegq, giec, gjec, gisc, gjsc, cnt, isc, jsc, iec, jec
  integer :: i_off, j_off

  G => CS%grid

  isym = 0
  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  isd = G%isd ; jsd = G%jsd
  iegq = G%iegB ; jegq = G%jegB
  gisc = G%domain%nihalo+1 ; gjsc = G%domain%njhalo+1
  giec = G%domain%niglobal+G%domain%nihalo ; gjec = G%domain%njglobal+G%domain%njhalo
  is = iscq - (1-isym); js = jscq - (1-isym)
  i_off = G%idg_offset ; j_off = G%jdg_offset

  D => G%bathyT
  H => CS%h_shelf
  float_frac => CS%float_frac
  hmask => CS%hmask
  u_face_mask => CS%u_face_mask
  v_face_mask => CS%v_face_mask
  rho = CS%density_ice
  rhow = CS%density_ocean_avg

  call savearray2 ("H",H,CS%write_output_to_file)
!  call savearray2 ("hmask",hmask,CS%write_output_to_file)
  call savearray2 ("u_face_mask", CS%u_face_mask_boundary,CS%write_output_to_file)
  call savearray2 ("umask", CS%umask,CS%write_output_to_file)
  call savearray2 ("v_face_mask", CS%v_face_mask_boundary,CS%write_output_to_file)
  call savearray2 ("vmask", CS%vmask,CS%write_output_to_file)

!   if (G%symmetric) then
!     isym=1
!   else
!     isym=0
!   endif

  isym = 0

  ! prelim - go through and calculate S

  ! or is this faster?
  BASE(:,:) = -D(:,:) + OD(:,:)
  S(:,:) = BASE(:,:) + H(:,:)

!  write (procnum,'(I1)') mpp_pe()

  do j=jsc-1,jec+1
    do i=isc-1,iec+1
      cnt = 0
      sx = 0
      sy = 0
      dxh = G%dxT(i,j)
      dyh = G%dyT(i,j)
      dxdyh = G%areaT(i,j)
!     print *,dxh," ",dyh," ",dxdyh

      if (hmask(i,j) .eq. 1) then ! we are inside the global computational bdry, at an ice-filled cell

        ! calculate sx
        if ((i+i_off) .eq. gisc) then ! at left computational bdry
          if (hmask(i+1,j) .eq. 1) then
            sx = (S(i+1,j)-S(i,j))/dxh
          else
            sx = 0
          endif
        elseif ((i+i_off) .eq. giec) then ! at right computational bdry
          if (hmask(i-1,j) .eq. 1) then
            sx = (S(i,j)-S(i-1,j))/dxh
          else
            sx=0
          endif
        else ! interior
          if (hmask(i+1,j) .eq. 1) then
            cnt = cnt+1
                sx = S(i+1,j)
          else
            sx = S(i,j)
              endif
              if (hmask(i-1,j) .eq. 1) then
            cnt = cnt+1
                sx = sx - S(i-1,j)
          else
            sx = sx - S(i,j)
              endif
          if (cnt .eq. 0) then
            sx=0
          else
            sx = sx / (cnt * dxh)
          endif
        endif

        cnt = 0

        ! calculate sy, similarly
        if ((j+j_off) .eq. gjsc) then ! at south computational bdry
          if (hmask(i,j+1) .eq. 1) then
            sy = (S(i,j+1)-S(i,j))/dyh
          else
            sy = 0
          endif
        elseif ((j+j_off) .eq. gjec) then ! at nprth computational bdry
          if (hmask(i,j-1) .eq. 1) then
            sy = (S(i,j)-S(i,j-1))/dyh
          else
            sy = 0
          endif
        else ! interior
          if (hmask(i,j+1) .eq. 1) then
            cnt = cnt+1
            sy = S(i,j+1)
          else
            sy = S(i,j)
          endif
          if (hmask(i,j-1) .eq. 1) then
            cnt = cnt+1
            sy = sy - S(i,j-1)
          else
            sy = sy - S(i,j)
          endif
          if (cnt .eq. 0) then
            sy=0
          else
            sy = sy / (cnt * dyh)
          endif
        endif


        if (FE .eq. 1) then

          ! SW vertex
          taud_x (i-1,j-1) = taud_x (i-1,j-1) - .25 * rho * grav * H(i,j) * sx * dxdyh
          taud_y(i-1,j-1) = taud_y(i-1,j-1) - .25 * rho * grav * H(i,j) * sy * dxdyh

          ! SE vertex
          taud_x(i,j-1) = taud_x(i,j-1) - .25 * rho * grav * H(i,j) * sx * dxdyh
          taud_y(i,j-1) = taud_y(i,j-1) - .25 * rho * grav * H(i,j) * sy * dxdyh

          ! NW vertex
          taud_x(i-1,j) = taud_x(i-1,j) - .25 * rho * grav * H(i,j) * sx * dxdyh
          taud_y(i-1,j) = taud_y(i-1,j) - .25 * rho * grav * H(i,j) * sy * dxdyh

          ! NE vertex
          taud_x(i,j) = taud_x(i,j) - .25 * rho * grav * H(i,j) * sx * dxdyh
          taud_y(i,j) = taud_y(i,j) - .25 * rho * grav * H(i,j) * sy * dxdyh


        else

          ! SW vertex
          taud_x(i-1,j-1) = taud_x(i-1,j-1) - (1./6) * rho * grav * H(i,j) * sx * dxdyh
          taud_y(i-1,j-1) = taud_y(i-1,j-1) - (1./6) * rho * grav * H(i,j) * sy * dxdyh

          ! SE vertex
          taud_x(i,j-1) = taud_x(i,j-1) - (1./3) * rho * grav * H(i,j) * sx * dxdyh
          taud_y(i,j-1) = taud_y(i,j-1) - (1./3) * rho * grav * H(i,j) * sy * dxdyh

          ! NW vertex
          taud_x(i-1,j) = taud_x(i-1,j) - (1./3) * rho * grav * H(i,j) * sx * dxdyh
          taud_y(i-1,j) = taud_y(i-1,j) - (1./3) * rho * grav * H(i,j) * sy * dxdyh

          ! NE vertex
          taud_x(i,j) = taud_x(i,j) - (1./6) * rho * grav * H(i,j) * sx * dxdyh
          taud_y(i,j) = taud_y(i,j) - (1./6) * rho * grav * H(i,j) * sy * dxdyh

        endif

        if (float_frac(i,j) .eq. 1) then
          neumann_val = .5 * grav * (rho * H (i,j) ** 2 - rhow * D(i,j) ** 2)
        else
          neumann_val = .5 * grav * (1-rho/rhow) * rho * H(i,j) ** 2
        endif


        if ((u_face_mask(i-1,j) .eq. 2) .OR. (hmask(i-1,j) .eq. 0) .OR. (hmask(i-1,j) .eq. 2) ) then ! left face of the cell is at a stress boundary
          ! the depth-integrated longitudinal stress is equal to the difference of depth-integrated pressure on either side of the face
          ! on the ice side, it is rho g h^2 / 2
          ! on the ocean side, it is rhow g (delta OD)^2 / 2
          ! OD can be zero under the ice; but it is ASSUMED on the ice-free side of the face, topography elevation is not above the base of the
          !     ice in the current cell
          taud_x(i-1,j-1) = taud_x(i-1,j-1) - .5 * dyh * neumann_val  ! note negative sign is due to direction of normal vector
          taud_x(i-1,j) = taud_x(i-1,j) - .5 * dyh * neumann_val
        endif

        if ((u_face_mask(i,j) .eq. 2) .OR. (hmask(i+1,j) .eq. 0) .OR. (hmask(i+1,j) .eq. 2) ) then ! right face of the cell is at a stress boundary
          taud_x(i,j-1) = taud_x(i,j-1) + .5 * dyh * neumann_val
          taud_x(i,j) = taud_x(i,j) + .5 * dyh * neumann_val
        endif

        if ((v_face_mask(i,j-1) .eq. 2) .OR. (hmask(i,j-1) .eq. 0) .OR. (hmask(i,j-1) .eq. 2) ) then ! south face of the cell is at a stress boundary
          taud_y(i-1,j-1) = taud_y(i-1,j-1) - .5 * dxh * neumann_val
          taud_y(i,j-1) = taud_y(i,j-1) - .5 * dxh * neumann_val
        endif

        if ((v_face_mask(i,j) .eq. 2) .OR. (hmask(i,j+1) .eq. 0) .OR. (hmask(i,j+1) .eq. 2) ) then ! north face of the cell is at a stress boundary
          taud_y(i-1,j) = taud_y(i-1,j) + .5 * dxh * neumann_val ! note negative sign is due to direction of normal vector
          taud_y(i,j) = taud_y(i,j) + .5 * dxh * neumann_val
        endif

      endif
    enddo
  enddo


!   call savearray2 ("Taux"//"p"//procnum,taud_x,CS%write_output_to_file)
!   call savearray2 ("Tauy"//"p"//procnum,taud_y,CS%write_output_to_file)

end subroutine calc_shelf_driving_stress

subroutine init_boundary_values (CS, time, input_flux, input_thick, new_sim)
  type(time_type),       intent(in)    :: Time
  type(ice_shelf_CS),    pointer       :: CS
  real, intent(in)               :: input_flux, input_thick
  logical, optional               :: new_sim

! this will be a per-setup function. the boundary values of thickness and velocity
! (and possibly other variables) will be updated in this function

! FOR RESTARTING PURPOSES: if grid is not symmetric and the model is restarted, we will
!               need to update those velocity points not *technically* in any
!               computational domain -- if this function gets moves to another module,
!               DO NOT TAKE THE RESTARTING BIT WITH IT

  real, dimension (:,:) , pointer      :: thickness_boundary_values, &
                          u_boundary_values, &
                          v_boundary_values, &
                          u_face_mask, v_face_mask, hmask
  type(ocean_grid_type), pointer :: G
  integer :: isym, i, j, iscq, iecq, jscq, jecq, isd, jsd, ied, jed, iegq, jegq, giec, gjec, gisc, gjsc, cnt, isc, jsc, iec, jec
  integer :: i_off, j_off
  real :: A, n, ux, uy, vx, vy, eps_min, domain_width

  G => CS%grid

!   if (G%symmetric) then
!     isym=1
!   else
!     isym=0
!   endif

  isym = 0

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
!   iscq = G%iscq ; iecq = G%iecq ; jscq = G%jscq ; jecq = G%jecq
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
!   iegq = G%iegq ; jegq = G%jegq
  i_off = G%idg_offset ; j_off = G%jdg_offset

  thickness_boundary_values => CS%thickness_boundary_values
  u_boundary_values => CS%u_boundary_values ; v_boundary_values => CS%v_boundary_values
  u_face_mask => CS%u_face_mask ; v_face_mask => CS%v_face_mask ; hmask => CS%hmask

  domain_width = CS%len_lat

  ! this loop results in some values being set twice but... eh.

  do j=jsd,jed
    do i=isd,ied

!      if ((i .eq. 4) .AND. ((mpp_pe() .eq. 0) .or. (mpp_pe() .eq. 6))) then
!    print *,hmask(i,j),i,j,mpp_pe()
!      endif

      if (hmask(i,j) .eq. 3) then
        thickness_boundary_values (i,j) = input_thick
      endif

      if ((hmask(i,j) .eq. 0) .or. (hmask(i,j) .eq. 1) .or. (hmask(i,j) .eq. 2)) then
        if ((i.le.iec).and.(i.ge.isc)) then
          if (u_face_mask (i-1,j) .eq. 3) then
            u_boundary_values (i-1,j-1) = (1 - ((G%geoLatBu(i-1,j-1) - 0.5*CS%len_lat)*2./CS%len_lat)**2) * &
                  1.5 * input_flux / input_thick
            u_boundary_values (i-1,j) = (1 - ((G%geoLatBu(i-1,j) - 0.5*CS%len_lat)*2./CS%len_lat)**2) * &
                  1.5 * input_flux / input_thick
          endif
        endif
      endif

      if (.not.(new_sim)) then
        if (.not. G%symmetric) then
          if (((i+i_off) .eq. (G%domain%nihalo+1)).and.(u_face_mask(i-1,j).eq.3)) then
            CS%u_shelf (i-1,j-1) = u_boundary_values (i-1,j-1)
            CS%u_shelf (i-1,j) = u_boundary_values (i-1,j)
!            print *, u_boundary_values (i-1,j)
          endif
          if (((j+j_off) .eq. (G%domain%njhalo+1)).and.(v_face_mask(i,j-1).eq.3)) then
            CS%u_shelf (i-1,j-1) = u_boundary_values (i-1,j-1)
            CS%u_shelf (i,j-1) = u_boundary_values (i,j-1)
          endif
        endif
      endif
    enddo
  enddo

end subroutine init_boundary_values

subroutine CG_action_triangular (uret, vret, u, v, umask, vmask, hmask, nu_upper, nu_lower, &
                beta_upper, beta_lower, dxh, dyh, dxdyh, is, ie, js, je, isym)

real, dimension (:,:), intent (inout)  :: uret, vret
real, dimension (:,:), intent (in)     :: u, v
real, dimension (:,:), intent (in)     :: umask, vmask
real, dimension (:,:), intent (in)     :: hmask, nu_upper, nu_lower, beta_upper, beta_lower
real, dimension (:,:), intent (in)     :: dxh, dyh, dxdyh
integer, intent(in)               :: is, ie, js, je, isym

! the linear action of the matrix on (u,v) with triangular finite elements
! as of now everything is passed in so no grid pointers or anything of the sort have to be dereferenced,
! but this may change pursuant to conversations with others
!
! is & ie are the cells over which the iteration is done; this may change between calls to this subroutine
!     in order to make less frequent halo updates
! isym = 1 if grid is symmetric, 0 o.w.

  real :: ux, uy, vx, vy
  integer :: i,j

  do i=is,ie
    do j=js,je

      if (hmask(i,j) .eq. 1) then ! this cell's vertices contain degrees of freedom

        ux = (u(i,j-1)-u(i-1,j-1))/dxh(i,j)
        vx = (v(i,j-1)-v(i-1,j-1))/dxh(i,j)
        uy = (u(i-1,j)-u(i-1,j-1))/dyh(i,j)
        vy = (v(i-1,j)-v(i-1,j-1))/dyh(i,j)

        if (umask(i,j-1) .eq. 1) then ! this (bot right) is a degree of freedom node

          uret(i,j-1) = uret(i,j-1) + &
              .5 * dxdyh(i,j) * nu_lower (i,j) * ((4*ux+2*vy) * (1./dxh(i,j)) + (uy+vy) * (0./dyh(i,j)))

          vret(i,j-1) = vret(i,j-1) + &
              .5 * dxdyh(i,j) * nu_lower (i,j) * ((uy+vx) * (1./dxh(i,j)) + (4*vy+2*ux) * (0./dyh(i,j)))

          uret(i,j-1) = uret(i,j-1) + &
              beta_lower(i,j) * dxdyh(i,j) * 1./24 * (u(i-1,j-1) + &
                                      u(i-1,j) + u(i,j-1))

          vret(i,j-1) = vret(i,j-1) + &
              beta_lower(i,j) * dxdyh(i,j) * 1./24 * (v(i-1,j-1) + &
                                      v(i-1,j) + v(i,j-1))
        endif

        if (umask(i-1,j) .eq. 1) then ! this (top left) is a degree of freedom node

          uret(i-1,j) = uret(i-1,j) + &
              .5 * dxdyh(i,j) * nu_lower (i,j) * ((4*ux+2*vy) * (0./dxh(i,j)) + (uy+vy) * (1./dyh(i,j)))

          vret(i-1,j) = vret(i-1,j) + &
              .5 * dxdyh(i,j) * nu_lower (i,j) * ((uy+vx) * (0./dxh(i,j)) + (4*vy+2*ux) * (1./dyh(i,j)))

          uret(i,j-1) = uret(i,j-1) + &
              beta_lower(i,j) * dxdyh(i,j) * 1./24 * (u(i-1,j-1) + &
                                      u(i-1,j) + u(i,j-1))

          vret(i,j-1) = vret(i,j-1) + &
              beta_lower(i,j) * dxdyh(i,j) * 1./24 * (v(i-1,j-1) + &
                                      v(i-1,j) + v(i,j-1))
        endif

        if (umask(i-1,j-1) .eq. 1) then ! this (bot left) is a degree of freedom node

          uret(i-1,j-1) = uret(i-1,j-1) + &
              .5 * dxdyh(i,j) * nu_upper (i,j) * ((4*ux+2*vy) * (-1./dxh(i,j)) + (uy+vy) * (-1./dyh(i,j)))

          vret(i-1,j-1) = vret(i-1,j-1) + &
              .5 * dxdyh(i,j) * nu_upper (i,j) * ((uy+vx) * (-1./dxh(i,j)) + (4*vy+2*ux) * (-1./dyh(i,j)))

          uret(i-1,j-1) = uret(i-1,j-1) + &
              beta_lower(i,j) * dxdyh(i,j) * 1./24 * (u(i-1,j-1) + &
                                      u(i-1,j) + u(i,j-1))

          vret(i-1,j-1) = vret(i-1,j-1) + &
              beta_lower(i,j) * dxdyh(i,j) * 1./24 * (v(i-1,j-1) + &
                                      v(i-1,j) + v(i,j-1))
        endif


        ux = (u(i,j)-u(i-1,j))/dxh(i,j)
        vx = (v(i,j)-v(i-1,j))/dxh(i,j)
        uy = (u(i,j)-u(i,j-1))/dyh(i,j)
        vy = (v(i,j)-v(i,j-1))/dyh(i,j)

        if (umask(i,j-1) .eq. 1) then ! this (bot right) is a degree of freedom node

          uret(i,j-1) = uret(i,j-1) + &
              .5 * dxdyh(i,j) * nu_upper (i,j) * ((4*ux+2*vy) * (0./dxh(i,j)) + (uy+vy) * (-1./dyh(i,j)))

          vret(i,j-1) = vret(i,j-1) + &
              .5 * dxdyh(i,j) * nu_upper (i,j) * ((uy+vx) * (0./dxh(i,j)) + (4*vy+2*ux) * (-1./dyh(i,j)))

          uret(i,j-1) = uret(i,j-1) + &
              beta_upper(i,j) * dxdyh(i,j) * 1./24 * (u(i,j) + &
                                      u(i-1,j) + u(i,j-1))

          vret(i,j-1) = vret(i,j-1) + &
              beta_upper(i,j) * dxdyh(i,j) * 1./24 * (u(i,j) + &
                                      u(i-1,j) + u(i,j-1))
        endif

        if (umask(i-1,j) .eq. 1) then ! this (top left) is a degree of freedom node

          uret(i-1,j) = uret(i-1,j) + &
              .5 * dxdyh(i,j) * nu_upper (i,j) * ((4*ux+2*vy) * (-1./dxh(i,j)) + (uy+vy) * (0./dyh(i,j)))

          vret(i-1,j) = vret(i-1,j) + &
              .5 * dxdyh(i,j) * nu_upper (i,j) * ((uy+vx) * (-1./dxh(i,j)) + (4*vy+2*ux) * (0./dyh(i,j)))

          uret(i,j-1) = uret(i,j-1) + &
              beta_upper(i,j) * dxdyh(i,j) * 1./24 * (u(i,j) + &
                                      u(i-1,j) + u(i,j-1))

          vret(i,j-1) = vret(i,j-1) + &
              beta_upper(i,j) * dxdyh(i,j) * 1./24 * (u(i,j) + &
                                      u(i-1,j) + u(i,j-1))
        endif

        if (umask(i,j) .eq. 1) then ! this (top right) is a degree of freedom node

          uret(i,j) = uret(i,j) + &
              .5 * dxdyh(i,j) * nu_upper (i,j) * ((4*ux+2*vy) * (1./dxh(i,j)) + (uy+vy) * (1./dyh(i,j)))

          vret(i,j) = vret(i,j) + &
              .5 * dxdyh(i,j) * nu_upper (i,j) * ((uy+vx) * (1./dxh(i,j)) + (4*vy+2*ux) * (1./dyh(i,j)))

          uret(i,j) = uret(i,j) + &
              beta_upper(i,j) * dxdyh(i,j) * 1./24 * (u(i,j) + &
                                      u(i-1,j) + u(i,j-1))

          vret(i,j) = vret(i,j) + &
              beta_upper(i,j) * dxdyh(i,j) * 1./24 * (u(i,j) + &
                                      u(i-1,j) + u(i,j-1))
        endif

      endif

    enddo
  enddo

end subroutine CG_action_triangular

subroutine CG_action_bilinear (uret, vret, u, v, Phi, Phisub, umask, vmask, hmask, H_node, &
                nu, float_cond, D, beta, dxdyh, is, ie, js, je, dens_ratio)

real, dimension (NILIMB_SYM_,NJLIMB_SYM_), intent (inout)  :: uret, vret
real, dimension (:,:,:,:), pointer :: Phi
real, dimension (:,:,:,:,:,:),pointer :: Phisub
real, dimension (NILIMB_SYM_,NJLIMB_SYM_), intent (in)     :: u, v
real, dimension (NILIMB_SYM_,NJLIMB_SYM_), intent (in)     :: umask, vmask, H_node
real, dimension (:,:), intent (in)     :: hmask, nu, float_cond, D, beta, dxdyh
real, intent(in)                       :: dens_ratio
integer, intent(in)               :: is, ie, js, je

! the linear action of the matrix on (u,v) with triangular finite elements
! as of now everything is passed in so no grid pointers or anything of the sort have to be dereferenced,
! but this may change pursuant to conversations with others
!
! is & ie are the cells over which the iteration is done; this may change between calls to this subroutine
!     in order to make less frequent halo updates
! isym = 1 if grid is symmetric, 0 o.w.

! the linear action of the matrix on (u,v) with triangular finite elements
! Phi has the form
! Phi (i,j,k,q) - applies to cell i,j

    !  3 - 4
    !  |   |
    !  1 - 2

! Phi (i,j,2*k-1,q) gives d(Phi_k)/dx at quadrature point q
! Phi (i,j,2*k,q) gives d(Phi_k)/dy at quadrature point q
! Phi_k is equal to 1 at vertex k, and 0 at vertex l .ne. k, and bilinear

  real :: ux, vx, uy, vy, uq, vq, area, basel
  integer :: iq, jq, iphi, jphi, i, j, ilq, jlq
  real, dimension(2) :: xquad
  real, dimension(2,2) :: Ucell,Vcell,Hcell,Usubcontr,Vsubcontr,Ucontr

  xquad(1) = .5 * (1-sqrt(1./3)) ; xquad(2) = .5 * (1+sqrt(1./3))

  do j=js,je
    do i=is,ie ; if (hmask(i,j) .eq. 1) then
!     dxh = G%dxh(i,j)
!     dyh = G%dyh(i,j)
!
!     X(:,:) = geolonq (i-1:i,j-1:j)
!     Y(:,:) = geolatq (i-1:i,j-1:j)
!
!     call bilinear_shape_functions (X, Y, Phi, area)

    ! X and Y must be passed in the form
        !  3 - 4
        !  |   |
        !  1 - 2
    ! Phi (2*i-1,j) gives d(Phi_i)/dx at quadrature point j
    ! Phi (2*i,j) gives d(Phi_i)/dy at quadrature point j

      area = dxdyh(i,j)

        Ucontr=0
      do iq=1,2 ; do jq=1,2


        if (iq .eq. 2) then
            ilq = 2
        else
            ilq = 1
        endif

        if (jq .eq. 2) then
            jlq = 2
        else
            jlq = 1
        endif

        uq = u(i-1,j-1) * xquad(3-iq) * xquad(3-jq) + &
        u(i,j-1) * xquad(iq) * xquad(3-jq) + &
        u(i-1,j) * xquad(3-iq) * xquad(jq) + &
        u(i,j) * xquad(iq) * xquad(jq)

        vq = v(i-1,j-1) * xquad(3-iq) * xquad(3-jq) + &
        v(i,j-1) * xquad(iq) * xquad(3-jq) + &
        v(i-1,j) * xquad(3-iq) * xquad(jq) + &
        v(i,j) * xquad(iq) * xquad(jq)

        ux = u(i-1,j-1) * Phi(i,j,1,2*(jq-1)+iq) + &
        u(i,j-1) * Phi(i,j,3,2*(jq-1)+iq) + &
        u(i-1,j) * Phi(i,j,5,2*(jq-1)+iq) + &
        u(i,j) * Phi(i,j,7,2*(jq-1)+iq)

        vx = v(i-1,j-1) * Phi(i,j,1,2*(jq-1)+iq) + &
        v(i,j-1) * Phi(i,j,3,2*(jq-1)+iq) + &
        v(i-1,j) * Phi(i,j,5,2*(jq-1)+iq) + &
        v(i,j) * Phi(i,j,7,2*(jq-1)+iq)

        uy = u(i-1,j-1) * Phi(i,j,2,2*(jq-1)+iq) + &
        u(i,j-1) * Phi(i,j,4,2*(jq-1)+iq) + &
        u(i-1,j) * Phi(i,j,6,2*(jq-1)+iq) + &
        u(i,j) * Phi(i,j,8,2*(jq-1)+iq)

        vy = v(i-1,j-1) * Phi(i,j,2,2*(jq-1)+iq) + &
        v(i,j-1) * Phi(i,j,4,2*(jq-1)+iq) + &
        v(i-1,j) * Phi(i,j,6,2*(jq-1)+iq) + &
        v(i,j) * Phi(i,j,8,2*(jq-1)+iq)

        do iphi=1,2 ; do jphi=1,2
          if (umask (i-2+iphi,j-2+jphi) .eq. 1) then

            uret (i-2+iphi,j-2+jphi) = uret (i-2+iphi,j-2+jphi) + &
                .25 * area * nu (i,j) * ((4*ux+2*vy) * Phi(i,j,2*(2*(jphi-1)+iphi)-1,2*(jq-1)+iq) + &
                                (uy+vx) * Phi(i,j,2*(2*(jphi-1)+iphi),2*(jq-1)+iq))
          endif
          if (vmask (i-2+iphi,j-2+jphi) .eq. 1) then

            vret (i-2+iphi,j-2+jphi) = vret (i-2+iphi,j-2+jphi) + &
                .25 * area * nu (i,j) * ((uy+vx) * Phi(i,j,2*(2*(jphi-1)+iphi)-1,2*(jq-1)+iq) + &
                                (4*vy+2*ux) * Phi(i,j,2*(2*(jphi-1)+iphi),2*(jq-1)+iq))
          endif

          if (iq .eq. iphi) then
            ilq = 2
          else
            ilq = 1
          endif

          if (jq .eq. jphi) then
            jlq = 2
          else
            jlq = 1
          endif

          if (float_cond(i,j) .eq. 0) then

            if (umask (i-2+iphi,j-2+jphi) .eq. 1) then

              uret (i-2+iphi,j-2+jphi) = uret (i-2+iphi,j-2+jphi) + &
                .25 * beta(i,j) * area * uq * xquad(ilq) * xquad(jlq)

            endif

            if (vmask (i-2+iphi,j-2+jphi) .eq. 1) then

              vret (i-2+iphi,j-2+jphi) = vret (i-2+iphi,j-2+jphi) + &
                .25 * beta(i,j) * area * vq * xquad(ilq) * xquad(jlq)

            endif

          endif
              Ucontr(iphi,jphi) = Ucontr(iphi,jphi) + .25 * area * uq * xquad(ilq) * xquad(jlq) * beta(i,j)
!              if((i.eq.27) .and. (j.eq.8) .and. (iphi.eq.1) .and. (jphi.eq.1))  print *, "grid", uq, .25 * area * uq * xquad(ilq) * xquad(jlq)

          !endif
        enddo ; enddo
      enddo ; enddo

      if (float_cond(i,j) .eq. 1) then
        Usubcontr = 0.0 ; Vsubcontr = 0.0 ; basel = D(i,j)
        Ucell(:,:) = u(i-1:i,j-1:j) ; Vcell(:,:) = v(i-1:i,j-1:j) ; Hcell(:,:) = H_node(i-1:i,j-1:j)
        call CG_action_subgrid_basal_bilinear &
            (Phisub, Hcell, Ucell, Vcell, area, basel, dens_ratio, Usubcontr, Vsubcontr, i, j)
        do iphi=1,2 ; do jphi=1,2
          if (umask (i-2+iphi,j-2+jphi) .eq. 1) then
            uret (i-2+iphi,j-2+jphi) = uret (i-2+iphi,j-2+jphi) + Usubcontr (iphi,jphi) * beta(i,j)
          endif
          if (vmask (i-2+iphi,j-2+jphi) .eq. 1) then
            vret (i-2+iphi,j-2+jphi) = vret (i-2+iphi,j-2+jphi) + Vsubcontr (iphi,jphi) * beta(i,j)
            !if ( (iphi.eq.1) .and. (jphi.eq.1)) print *,  i,j, Usubcontr (iphi,jphi) * beta(i,j), " ", Ucontr(iphi,jphi)
          endif
        enddo ; enddo
      endif

    endif
  enddo ; enddo

end subroutine CG_action_bilinear

subroutine CG_action_subgrid_basal_bilinear (Phisub, H, U, V, DXDYH, D, dens_ratio, Ucontr, Vcontr, iin, jin)
  real, pointer, dimension(:,:,:,:,:,:) :: Phisub
  real, dimension(2,2), intent(in) :: H,U,V
  real, intent(in)                 :: DXDYH, D, dens_ratio
  real, dimension(2,2), intent(inout) :: Ucontr, Vcontr
  integer, optional, intent(in)              :: iin, jin

  ! D = cellwise-constant bed elevation

  integer              :: nsub, i, j, k, l, qx, qy, m, n, i_m, j_m
  real                 :: subarea, hloc, uq, vq

  nsub = size(Phisub,1)
  subarea = DXDYH / (nsub**2)


  if (.not. present(iin)) then
   i_m = -1
  else
   i_m = iin
  endif

  if (.not. present(jin)) then
   j_m = -1
  else
   j_m = jin
  endif


  do m=1,2
    do n=1,2
      do j=1,nsub
        do i=1,nsub
          do qx=1,2
            do qy = 1,2

              hloc = Phisub(i,j,1,1,qx,qy)*H(1,1)+Phisub(i,j,1,2,qx,qy)*H(1,2)+&
                Phisub(i,j,2,1,qx,qy)*H(2,1)+Phisub(i,j,2,2,qx,qy)*H(2,2)

              if (dens_ratio * hloc - D .gt. 0) then
              !if (.true.) then
                uq = 0 ; vq = 0
                do k=1,2
                  do l=1,2
                    !Ucontr (m,n) = Ucontr (m,n) + subarea * 0.25 * Phisub(i,j,m,n,qx,qy) * Phisub(i,j,k,l,qx,qy) * U(k,l)
                    !Vcontr (m,n) = Vcontr (m,n) + subarea * 0.25 * Phisub(i,j,m,n,qx,qy) * Phisub(i,j,k,l,qx,qy) * V(k,l)
                    uq = uq + Phisub(i,j,k,l,qx,qy) * U(k,l) ; vq = vq + Phisub(i,j,k,l,qx,qy) * V(k,l)
                  enddo
                enddo

                Ucontr (m,n) = Ucontr (m,n) + subarea * 0.25 * Phisub(i,j,m,n,qx,qy) * uq
                Vcontr (m,n) = Vcontr (m,n) + subarea * 0.25 * Phisub(i,j,m,n,qx,qy) * vq

 !               if ((i_m .eq. 27) .and. (j_m .eq. 8) .and. (m.eq.1) .and. (n.eq.1)) print *, "in subgrid", uq,  Phisub(i,j,m,n,qx,qy)

              endif

            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine CG_action_subgrid_basal_bilinear

subroutine matrix_diagonal_triangle (CS, u_diagonal, v_diagonal)

  type(ice_shelf_CS),    pointer       :: CS
  real, dimension (:,:), intent(inout) :: u_diagonal, v_diagonal

! returns the diagonal entries of the matrix for a Jacobi preconditioning

  real, pointer, dimension (:,:)       :: umask, vmask, &
                          nu_lower, nu_upper, beta_lower, beta_upper, hmask
  type(ocean_grid_type), pointer :: G
  integer :: isym, i, j, is, js, cnt, isc, jsc, iec, jec
  real :: A, n, ux, uy, vx, vy, eps_min, domain_width, dxh, dyh, dxdyh

  G => CS%grid

!   if (G%symmetric) then
!     isym=1
!   else
!     isym=0
!   endif

  isym = 0

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  umask => CS%umask ; vmask => CS%vmask ; hmask => CS%hmask
  nu_lower => CS%ice_visc_lower_tri ; nu_upper => CS%ice_visc_upper_tri
  beta_lower => CS%taub_beta_eff_lower_tri ; beta_upper => CS%taub_beta_eff_upper_tri

  do i=isc-1,iec+1  ; do j=jsc-1,jec+1 ; if (hmask(i,j) .eq. 1) then
    dxh = G%dxT(i,j)
    dyh = G%dyT(i,j)
    dxdyh = G%areaT(i,j)

    if (umask (i,j-1) .eq. 1) then ! this (bot right) is a degree of freedom node

      ux = 1./dxh ; uy = 0./dyh
      vx = 0. ; vy = 0.

      u_diagonal (i,j-1) = u_diagonal (i,j-1) + &
          .5 * dxdyh * nu_lower (i,j) * ((4*ux+2*vy) * (1./dxh) + (uy+vy) * (0./dyh))

      u_diagonal (i,j-1) = u_diagonal (i,j-1) + &
          beta_lower(i,j) * dxdyh * 1./24

      ux = 0. ; uy = 0.
      vx = 1./dxh ; vy = 0./dyh

      v_diagonal (i,j-1) = v_diagonal (i,j-1) + &
          .5 * dxdyh * nu_lower (i,j) * ((uy+vx) * (1./dxh) + (4*vy+2*ux) * (0./dyh))

      v_diagonal (i,j-1) = v_diagonal (i,j-1) + &
          beta_lower(i,j) * dxdyh * 1./24

      ux = 0./dxh ; uy = -1./dyh
      vx = 0. ; vy = 0.

      u_diagonal (i,j-1) = u_diagonal (i,j-1) + &
          .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (0./dxh) + (uy+vy) * (-1./dyh))

      u_diagonal (i,j-1) = u_diagonal (i,j-1) + &
          beta_upper(i,j) * dxdyh * 1./24

      vx = 0./dxh ; vy = -1./dyh
      ux = 0. ; uy = 0.

      v_diagonal (i,j-1) = v_diagonal (i,j-1) + &
          .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (0./dxh) + (4*vy+2*ux) * (-1./dyh))

      v_diagonal (i,j-1) = v_diagonal (i,j-1) + &
          beta_upper(i,j) * dxdyh * 1./24

    endif

    if (umask (i-1,j) .eq. 1) then ! this (top left) is a degree of freedom node

      ux = 0./dxh ; uy = 1./dyh
      vx = 0. ; vy = 0.

      u_diagonal (i-1,j) = u_diagonal (i-1,j) + &
          .5 * dxdyh * nu_lower (i,j) * ((4*ux+2*vy) * (0./dxh) + (uy+vy) * (1./dyh))

      u_diagonal (i,j-1) = u_diagonal (i,j-1) + &
          beta_lower(i,j) * dxdyh * 1./24

      ux = 0. ; uy = 0.
      vx = 0./dxh ; vy = 1./dyh

      v_diagonal (i-1,j) = v_diagonal (i-1,j) + &
          .5 * dxdyh * nu_lower (i,j) * ((uy+vx) * (0./dxh) + (4*vy+2*ux) * (1./dyh))

      v_diagonal (i,j-1) = v_diagonal (i,j-1) + &
          beta_lower(i,j) * dxdyh * 1./24

      ux = -1./dxh ; uy = 0./dyh
      vx = 0. ; vy = 0.

      u_diagonal (i-1,j) = u_diagonal (i-1,j) + &
          .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (-1./dxh) + (uy+vy) * (0./dyh))

      u_diagonal (i,j-1) = u_diagonal (i,j-1) + &
          beta_upper(i,j) * dxdyh * 1./24

      vx = -1./dxh ; vy = 0./dyh
      ux = 0. ; uy = 0.

      v_diagonal (i-1,j) = v_diagonal (i-1,j) + &
          .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (-1./dxh) + (4*vy+2*ux) * (0./dyh))

      v_diagonal (i,j-1) = v_diagonal (i,j-1) + &
          beta_upper(i,j) * dxdyh * 1./24

    endif

    if (umask (i-1,j-1) .eq. 1) then ! this (bot left) is a degree of freedom node

      ux = -1./dxh ; uy = -1./dyh
      vx = 0. ; vy = 0.

      u_diagonal (i-1,j-1) = u_diagonal (i-1,j-1) + &
          .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (-1./dxh) + (uy+vy) * (-1./dyh))

      u_diagonal (i-1,j-1) = u_diagonal (i-1,j-1) + &
          beta_lower(i,j) * dxdyh * 1./24

      vx = -1./dxh ; vy = -1./dyh
      ux = 0. ; uy = 0.

      v_diagonal (i-1,j-1) = v_diagonal (i-1,j-1) + &
          .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (-1./dxh) + (4*vy+2*ux) * (-1./dyh))

      v_diagonal (i-1,j-1) = v_diagonal (i-1,j-1) + &
          beta_lower(i,j) * dxdyh * 1./24
    endif

    if (umask (i,j) .eq. 1) then ! this (top right) is a degree of freedom node

      ux = 1./ dxh ; uy = 1./dyh
      vx = 0. ; vy = 0.

      u_diagonal (i,j) = u_diagonal (i,j) + &
          .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (1./dxh) + (uy+vy) * (1./dyh))

      u_diagonal (i,j) = u_diagonal (i,j) + &
          beta_upper(i,j) * dxdyh * 1./24

      vx = 1./ dxh ; vy = 1./dyh
      ux = 0. ; uy = 0.

      v_diagonal (i,j) = v_diagonal (i,j) + &
          .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (1./dxh) + (4*vy+2*ux) * (1./dyh))

      v_diagonal (i,j) = v_diagonal (i,j) + &
          beta_upper(i,j) * dxdyh * 1./24

    endif
  endif ; enddo ; enddo

end subroutine matrix_diagonal_triangle

subroutine matrix_diagonal_bilinear(CS, float_cond, H_node, dens_ratio, Phisub, u_diagonal, v_diagonal)

  type(ice_shelf_CS),    pointer       :: CS
  real, dimension (NILIMB_SYM_,NJLIMB_SYM_), intent(in) :: H_node
  real                                :: dens_ratio
  real, dimension (:,:), intent(in) :: float_cond
  real, dimension (:,:,:,:,:,:),pointer :: Phisub
  real, dimension (NILIMB_SYM_,NJLIMB_SYM_), intent(inout) :: u_diagonal, v_diagonal


! returns the diagonal entries of the matrix for a Jacobi preconditioning

  real, dimension (:,:), pointer       :: umask, vmask, hmask, &
                          nu, beta
  type(ocean_grid_type), pointer :: G
  integer :: isym, i, j, is, js, cnt, isc, jsc, iec, jec, iphi, jphi, iq, jq, ilq, jlq
  real :: A, n, ux, uy, vx, vy, eps_min, domain_width, dxh, dyh, dxdyh, area, uq, vq, basel
  real, dimension(8,4)  :: Phi
  real, dimension(4) :: X, Y
  real, dimension(2) :: xquad
  real, dimension(2,2) :: Hcell,Usubcontr,Vsubcontr

  G => CS%grid

!   if (G%symmetric) then
!     isym=1
!   else
!     isym=0
!   endif

  isym = 0

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  umask => CS%umask ; vmask => CS%vmask ; hmask => CS%hmask
  nu => CS%ice_visc_bilinear
  beta => CS%taub_beta_eff_bilinear

  xquad(1) = .5 * (1-sqrt(1./3)) ; xquad(2) = .5 * (1+sqrt(1./3))

! X and Y must be passed in the form
    !  3 - 4
    !  |   |
    !  1 - 2
! Phi (2*i-1,j) gives d(Phi_i)/dx at quadrature point j
! Phi (2*i,j) gives d(Phi_i)/dy at quadrature point j

  do j=jsc-1,jec+1 ; do i=isc-1,iec+1 ; if (hmask(i,j) .eq. 1) then

    dxh = G%dxT(i,j)
    dyh = G%dyT(i,j)
    dxdyh = G%areaT(i,j)

    X(1:2) = G%geoLonBu (i-1:i,j-1)*1000
    X(3:4) = G%geoLonBu (i-1:i,j) *1000
    Y(1:2) = G%geoLatBu (i-1:i,j-1) *1000
    Y(3:4) = G%geoLatBu (i-1:i,j)*1000

    call bilinear_shape_functions (X, Y, Phi, area)

    ! X and Y must be passed in the form
        !  3 - 4
        !  |   |
        !  1 - 2
    ! Phi (2*i-1,j) gives d(Phi_i)/dx at quadrature point j
    ! Phi (2*i,j) gives d(Phi_i)/dy at quadrature point j

    do iq=1,2 ; do jq=1,2

      do iphi=1,2 ; do jphi=1,2

          if (iq .eq. iphi) then
            ilq = 2
          else
            ilq = 1
          endif

          if (jq .eq. jphi) then
            jlq = 2
          else
            jlq = 1
          endif

        if (umask (i-2+iphi,j-2+jphi) .eq. 1) then

          ux = Phi (2*(2*(jphi-1)+iphi)-1, 2*(jq-1)+iq)
          uy = Phi (2*(2*(jphi-1)+iphi), 2*(jq-1)+iq)
          vx = 0.
          vy = 0.

          u_diagonal (i-2+iphi,j-2+jphi) = u_diagonal (i-2+iphi,j-2+jphi) + &
              .25 * dxdyh * nu (i,j) * ((4*ux+2*vy) * Phi(2*(2*(jphi-1)+iphi)-1,2*(jq-1)+iq) + &
                              (uy+vy) * Phi(2*(2*(jphi-1)+iphi),2*(jq-1)+iq))

          uq = xquad(ilq) * xquad(jlq)

          if (float_cond(i,j) .eq. 0) then
            u_diagonal (i-2+iphi,j-2+jphi) = u_diagonal (i-2+iphi,j-2+jphi) + &
                .25 * beta(i,j) * dxdyh * uq * xquad(ilq) * xquad(jlq)
          endif

        endif

        if (vmask (i-2+iphi,j-2+jphi) .eq. 1) then

          vx = Phi (2*(2*(jphi-1)+iphi)-1, 2*(jq-1)+iq)
          vy = Phi (2*(2*(jphi-1)+iphi), 2*(jq-1)+iq)
          ux = 0.
          uy = 0.

          v_diagonal (i-2+iphi,j-2+jphi) = v_diagonal (i-2+iphi,j-2+jphi) + &
              .25 * dxdyh * nu (i,j) * ((uy+vx) * Phi(2*(2*(jphi-1)+iphi)-1,2*(jq-1)+iq) + &
                              (4*vy+2*ux) * Phi(2*(2*(jphi-1)+iphi),2*(jq-1)+iq))

          vq = xquad(ilq) * xquad(jlq)

          if (float_cond(i,j) .eq. 0) then
            v_diagonal (i-2+iphi,j-2+jphi) = v_diagonal (i-2+iphi,j-2+jphi) + &
                .25 * beta(i,j) * dxdyh * vq * xquad(ilq) * xquad(jlq)
          endif

        endif
      enddo ; enddo
    enddo ; enddo
    if (float_cond(i,j) .eq. 1) then
      Usubcontr = 0.0 ; Vsubcontr = 0.0 ; basel = G%bathyT(i,j)
      Hcell(:,:) = H_node(i-1:i,j-1:j)
      call CG_diagonal_subgrid_basal_bilinear &
          (Phisub, Hcell, dxdyh, basel, dens_ratio, Usubcontr, Vsubcontr)
      do iphi=1,2 ; do jphi=1,2
        if (umask (i-2+iphi,j-2+jphi) .eq. 1) then
          u_diagonal (i-2+iphi,j-2+jphi) = u_diagonal (i-2+iphi,j-2+jphi) + Usubcontr (iphi,jphi) * beta(i,j)
          v_diagonal (i-2+iphi,j-2+jphi) = v_diagonal (i-2+iphi,j-2+jphi) + Vsubcontr (iphi,jphi) * beta(i,j)
        endif
      enddo ; enddo
    endif
  endif ; enddo ; enddo

end subroutine matrix_diagonal_bilinear

subroutine CG_diagonal_subgrid_basal_bilinear (Phisub, H, DXDYH, D, dens_ratio, Ucontr, Vcontr)
  real, pointer, dimension(:,:,:,:,:,:) :: Phisub
  real, dimension(2,2), intent(in) :: H
  real, intent(in)                 :: DXDYH, D, dens_ratio
  real, dimension(2,2), intent(inout) :: Ucontr, Vcontr

  ! D = cellwise-constant bed elevation

  integer              :: nsub, i, j, k, l, qx, qy, m, n
  real                 :: subarea, hloc

  nsub = size(Phisub,1)
  subarea = DXDYH / (nsub**2)

  do m=1,2
    do n=1,2
      do j=1,nsub
        do i=1,nsub
          do qx=1,2
            do qy = 1,2

              hloc = Phisub(i,j,1,1,qx,qy)*H(1,1)+Phisub(i,j,1,2,qx,qy)*H(1,2)+&
                Phisub(i,j,2,1,qx,qy)*H(2,1)+Phisub(i,j,2,2,qx,qy)*H(2,2)

              if (dens_ratio * hloc - D .gt. 0) then
                Ucontr (m,n) = Ucontr (m,n) + subarea * 0.25 * Phisub(i,j,m,n,qx,qy)**2
                Vcontr (m,n) = Vcontr (m,n) + subarea * 0.25 * Phisub(i,j,m,n,qx,qy)**2
              endif


            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine CG_diagonal_subgrid_basal_bilinear


subroutine apply_boundary_values_triangle (CS, time, u_boundary_contr, v_boundary_contr)

  type(time_type),       intent(in)    :: Time
  type(ice_shelf_CS),    pointer       :: CS
  real, dimension (:,:), intent(inout) :: u_boundary_contr, v_boundary_contr

! this will be a per-setup function. the boundary values of thickness and velocity
! (and possibly other variables) will be updated in this function

  real, pointer, dimension (:,:)       :: u_boundary_values, &
                          v_boundary_values, &
                          umask, vmask, hmask, &
                          nu_lower, nu_upper, beta_lower, beta_upper
  type(ocean_grid_type), pointer :: G
  integer :: isym, i, j, cnt, isc, jsc, iec, jec
  real :: A, n, ux, uy, vx, vy, eps_min, domain_width, dxh, dyh, dxdyh

  G => CS%grid

!   if (G%symmetric) then
!     isym=1
!   else
!     isym=0
!   endif

  isym = 0

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  u_boundary_values => CS%u_boundary_values
  v_boundary_values => CS%v_boundary_values
  umask => CS%umask ; vmask => CS%vmask ; hmask => CS%hmask
  nu_lower => CS%ice_visc_lower_tri ; nu_upper => CS%ice_visc_upper_tri
  beta_lower => CS%taub_beta_eff_lower_tri ; beta_upper => CS%taub_beta_eff_upper_tri

  domain_width = CS%len_lat

  do i=isc-1,iec+1 ; do j=jsc-1,jec+1 ; if (hmask(i,j) .eq. 1) then

    if ((umask(i-1,j-1) .eq. 3) .OR. (umask(i,j-1) .eq. 3) .OR. (umask(i-1,j) .eq. 3)) then

      dxh = G%dxT(i,j)
      dyh = G%dyT(i,j)
      dxdyh = G%areaT(i,j)

      ux = (u_boundary_values(i,j-1)-u_boundary_values(i-1,j-1))/dxh
      vx = (v_boundary_values(i,j-1)-v_boundary_values(i-1,j-1))/dxh
      uy = (u_boundary_values(i-1,j)-u_boundary_values(i-1,j-1))/dyh
      vy = (v_boundary_values(i-1,j)-v_boundary_values(i-1,j-1))/dyh

      if (umask (i,j-1) .eq. 1) then ! this (bot right) is a degree of freedom node

        u_boundary_contr (i,j-1) = u_boundary_contr (i,j-1) + &
            .5 * dxdyh * nu_lower (i,j) * ((4*ux+2*vy) * (1./dxh) + (uy+vy) * (0./dyh))

        v_boundary_contr (i,j-1) = v_boundary_contr (i,j-1) + &
            .5 * dxdyh * nu_lower (i,j) * ((uy+vx) * (1./dxh) + (4*vy+2*ux) * (0./dyh))

        u_boundary_contr (i,j-1) = u_boundary_contr (i,j-1) + &
            beta_lower(i,j) * dxdyh * 1./24 * (u_boundary_values(i-1,j-1) + &
                          u_boundary_values(i-1,j) + u_boundary_values(i,j-1))

        v_boundary_contr (i,j-1) = v_boundary_contr (i,j-1) + &
            beta_lower(i,j) * dxdyh * 1./24 * (v_boundary_values(i-1,j-1) + &
                         v_boundary_values(i-1,j) + v_boundary_values(i,j-1))
      endif

      if (umask (i-1,j) .eq. 1) then ! this (top left) is a degree of freedom node

        u_boundary_contr (i-1,j) = u_boundary_contr (i-1,j) + &
            .5 * dxdyh * nu_lower (i,j) * ((4*ux+2*vy) * (0./dxh) + (uy+vy) * (1./dyh))

        v_boundary_contr (i-1,j) = v_boundary_contr (i-1,j) + &
            .5 * dxdyh * nu_lower (i,j) * ((uy+vx) * (0./dxh) + (4*vy+2*ux) * (1./dyh))

        u_boundary_contr (i,j-1) = u_boundary_contr (i,j-1) + &
            beta_lower(i,j) * dxdyh * 1./24 * (u_boundary_values(i-1,j-1) + &
                            u_boundary_values(i-1,j) + u_boundary_values(i,j-1))

        v_boundary_contr (i,j-1) = v_boundary_contr (i,j-1) + &
            beta_lower(i,j) * dxdyh * 1./24 * (v_boundary_values(i-1,j-1) + &
                            v_boundary_values(i-1,j) + v_boundary_values(i,j-1))
      endif

      if (umask (i-1,j-1) .eq. 1) then ! this (bot left) is a degree of freedom node

        u_boundary_contr (i-1,j-1) = u_boundary_contr (i-1,j-1) + &
            .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (-1./dxh) + (uy+vy) * (-1./dyh))

        v_boundary_contr (i-1,j-1) = v_boundary_contr (i-1,j-1) + &
            .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (-1./dxh) + (4*vy+2*ux) * (-1./dyh))

        u_boundary_contr (i-1,j-1) = u_boundary_contr (i-1,j-1) + &
            beta_lower(i,j) * dxdyh * 1./24 * (u_boundary_values(i-1,j-1) + &
                            u_boundary_values(i-1,j) + u_boundary_values(i,j-1))

        v_boundary_contr (i-1,j-1) = v_boundary_contr (i-1,j-1) + &
            beta_lower(i,j) * dxdyh * 1./24 * (v_boundary_values(i-1,j-1) + &
                            v_boundary_values(i-1,j) + v_boundary_values(i,j-1))
      endif

    endif

    if ((umask(i,j) .eq. 3) .OR. (umask(i,j-1) .eq. 3) .OR. (umask(i-1,j) .eq. 3)) then

      dxh = G%dxT(i,j)
      dyh = G%dyT(i,j)
      dxdyh = G%areaT(i,j)

      ux = (u_boundary_values(i,j)-u_boundary_values(i-1,j))/dxh
      vx = (v_boundary_values(i,j)-v_boundary_values(i-1,j))/dxh
      uy = (u_boundary_values(i,j)-u_boundary_values(i,j-1))/dyh
      vy = (v_boundary_values(i,j)-v_boundary_values(i,j-1))/dyh

      if (umask (i,j-1) .eq. 1) then ! this (bot right) is a degree of freedom node

          u_boundary_contr (i,j-1) = u_boundary_contr (i,j-1) + &
              .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (0./dxh) + (uy+vy) * (-1./dyh))

          v_boundary_contr (i,j-1) = v_boundary_contr (i,j-1) + &
              .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (0./dxh) + (4*vy+2*ux) * (-1./dyh))

          u_boundary_contr (i,j-1) = u_boundary_contr (i,j-1) + &
              beta_upper(i,j) * dxdyh * 1./24 * (u_boundary_values(i,j) + &
                                      u_boundary_values(i-1,j) +  &
                                u_boundary_values(i,j-1))

          v_boundary_contr (i,j-1) = v_boundary_contr (i,j-1) + &
              beta_upper(i,j) * dxdyh * 1./24 * (u_boundary_values(i,j) + &
                                      u_boundary_values(i-1,j) +  &
                                u_boundary_values(i,j-1))
      endif

      if (umask (i-1,j) .eq. 1) then ! this (top left) is a degree of freedom node

        u_boundary_contr (i-1,j) = u_boundary_contr (i-1,j) + &
            .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (-1./dxh) + (uy+vy) * (0./dyh))

        v_boundary_contr (i-1,j) = v_boundary_contr (i-1,j) + &
            .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (-1./dxh) + (4*vy+2*ux) * (0./dyh))

        u_boundary_contr (i,j-1) = u_boundary_contr (i,j-1) + &
            beta_upper(i,j) * dxdyh * 1./24 * (u_boundary_values(i,j) + &
                                    u_boundary_values(i-1,j) +  &
                              u_boundary_values(i,j-1))

        v_boundary_contr (i,j-1) = v_boundary_contr (i,j-1) + &
            beta_upper(i,j) * dxdyh * 1./24 * (u_boundary_values(i,j) + &
                                    u_boundary_values(i-1,j) +  &
                              u_boundary_values(i,j-1))
      endif

      if (umask (i,j) .eq. 1) then ! this (top right) is a degree of freedom node

        u_boundary_contr (i,j) = u_boundary_contr (i,j) + &
            .5 * dxdyh * nu_upper (i,j) * ((4*ux+2*vy) * (1./dxh) + (uy+vy) * (1./dyh))

        v_boundary_contr (i,j) = v_boundary_contr (i,j) + &
            .5 * dxdyh * nu_upper (i,j) * ((uy+vx) * (1./dxh) + (4*vy+2*ux) * (1./dyh))

        u_boundary_contr (i,j) = u_boundary_contr (i,j) + &
            beta_upper(i,j) * dxdyh * 1./24 * (u_boundary_values(i,j) + &
                                    u_boundary_values(i-1,j) +  &
                              u_boundary_values(i,j-1))

        v_boundary_contr (i,j) = v_boundary_contr (i,j) + &
            beta_upper(i,j) * dxdyh * 1./24 * (u_boundary_values(i,j) + &
                                    u_boundary_values(i-1,j) +  &
                              u_boundary_values(i,j-1))
      endif


    endif
  endif ; enddo ; enddo

end subroutine apply_boundary_values_triangle

subroutine apply_boundary_values_bilinear (CS, time, Phisub, H_node, float_cond, dens_ratio, u_boundary_contr, v_boundary_contr)

  type(time_type),       intent(in)    :: Time
  real, dimension (:,:,:,:,:,:),pointer:: Phisub
  type(ice_shelf_CS),    pointer       :: CS
  real, dimension (NILIMB_SYM_,NJLIMB_SYM_), intent (in)     :: H_node
  real, dimension (:,:), intent (in)   :: float_cond
  real                                 :: dens_ratio
  real, dimension (NILIMB_SYM_,NJLIMB_SYM_), intent(inout) :: u_boundary_contr, v_boundary_contr

! this will be a per-setup function. the boundary values of thickness and velocity
! (and possibly other variables) will be updated in this function

  real, pointer, dimension (:,:)       :: u_boundary_values, &
                          v_boundary_values, &
                          umask, vmask, &
                          nu, beta, hmask
  real, dimension(8,4)  :: Phi
  real, dimension(4) :: X, Y
  real, dimension(2) :: xquad
  type(ocean_grid_type), pointer :: G
  integer :: isym, i, j, isc, jsc, iec, jec, iq, jq, iphi, jphi, ilq, jlq
  real :: A, n, ux, uy, vx, vy, eps_min, domain_width, dxh, dyh, dxdyh, uq, vq, area, basel
  real, dimension(2,2) :: Ucell,Vcell,Hcell,Usubcontr,Vsubcontr

  G => CS%grid

!   if (G%symmetric) then
!     isym=1
!   else
!     isym=0
!   endif

  isym = 0

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  u_boundary_values => CS%u_boundary_values
  v_boundary_values => CS%v_boundary_values
  umask => CS%umask ; vmask => CS%vmask ; hmask => CS%hmask
  nu => CS%ice_visc_bilinear
  beta => CS%taub_beta_eff_bilinear

  xquad(1) = .5 * (1-sqrt(1./3)) ; xquad(2) = .5 * (1+sqrt(1./3))

! X and Y must be passed in the form
    !  3 - 4
    !  |   |
    !  1 - 2
! Phi (2*i-1,j) gives d(Phi_i)/dx at quadrature point j
! Phi (2*i,j) gives d(Phi_i)/dy at quadrature point j

  do j=jsc-1,jec+1 ; do i=isc-1,iec+1 ; if (hmask(i,j) .eq. 1) then

    ! process this cell if any corners have umask set to non-dirichlet bdry.
    ! NOTE: vmask not considered, probably should be

    if ((umask(i-1,j-1) .eq. 3) .OR. (umask(i,j-1) .eq. 3) .OR. &
        (umask(i-1,j) .eq. 3) .OR. (umask(i,j) .eq. 3)) then


      dxh = G%dxT(i,j)
      dyh = G%dyT(i,j)
      dxdyh = G%areaT(i,j)

      X(1:2) = G%geoLonBu (i-1:i,j-1)*1000
      X(3:4) = G%geoLonBu (i-1:i,j)*1000
      Y(1:2) = G%geoLatBu (i-1:i,j-1)*1000
      Y(3:4) = G%geoLatBu (i-1:i,j)*1000

      call bilinear_shape_functions (X, Y, Phi, area)

      ! X and Y must be passed in the form
          !  3 - 4
                 !  |   |
          !  1 - 2
      ! Phi (2*i-1,j) gives d(Phi_i)/dx at quadrature point j
      ! Phi (2*i,j) gives d(Phi_i)/dy at quadrature point j



      do iq=1,2 ; do jq=1,2

        uq = u_boundary_values(i-1,j-1) * xquad(3-iq) * xquad(3-jq) + &
              u_boundary_values(i,j-1) * xquad(iq) * xquad(3-jq) + &
             u_boundary_values(i-1,j) * xquad(3-iq) * xquad(jq) + &
             u_boundary_values(i,j) * xquad(iq) * xquad(jq)

        vq = v_boundary_values(i-1,j-1) * xquad(3-iq) * xquad(3-jq) + &
              v_boundary_values(i,j-1) * xquad(iq) * xquad(3-jq) + &
             v_boundary_values(i-1,j) * xquad(3-iq) * xquad(jq) + &
             v_boundary_values(i,j) * xquad(iq) * xquad(jq)

        ux = u_boundary_values(i-1,j-1) * Phi(1,2*(jq-1)+iq) + &
              u_boundary_values(i,j-1) * Phi(3,2*(jq-1)+iq) + &
             u_boundary_values(i-1,j) * Phi(5,2*(jq-1)+iq) + &
             u_boundary_values(i,j) * Phi(7,2*(jq-1)+iq)

        vx = v_boundary_values(i-1,j-1) * Phi(1,2*(jq-1)+iq) + &
              v_boundary_values(i,j-1) * Phi(3,2*(jq-1)+iq) + &
             v_boundary_values(i-1,j) * Phi(5,2*(jq-1)+iq) + &
             v_boundary_values(i,j) * Phi(7,2*(jq-1)+iq)

        uy = u_boundary_values(i-1,j-1) * Phi(2,2*(jq-1)+iq) + &
              u_boundary_values(i,j-1) * Phi(4,2*(jq-1)+iq) + &
             u_boundary_values(i-1,j) * Phi(6,2*(jq-1)+iq) + &
             u_boundary_values(i,j) * Phi(8,2*(jq-1)+iq)

        vy = v_boundary_values(i-1,j-1) * Phi(2,2*(jq-1)+iq) + &
              v_boundary_values(i,j-1) * Phi(4,2*(jq-1)+iq) + &
             v_boundary_values(i-1,j) * Phi(6,2*(jq-1)+iq) + &
             v_boundary_values(i,j) * Phi(8,2*(jq-1)+iq)

        do iphi=1,2 ; do jphi=1,2

          if (iq .eq. iphi) then
            ilq = 2
          else
            ilq = 1
          endif

          if (jq .eq. jphi) then
            jlq = 2
          else
            jlq = 1
          endif

          if (umask (i-2+iphi,j-2+jphi) .eq. 1) then


            u_boundary_contr (i-2+iphi,j-2+jphi) = u_boundary_contr (i-2+iphi,j-2+jphi) + &
            .25 * dxdyh * nu (i,j) * ( (4*ux+2*vy) * Phi(2*(2*(jphi-1)+iphi)-1,2*(jq-1)+iq) + &
                         (uy+vx) * Phi(2*(2*(jphi-1)+iphi),2*(jq-1)+iq) )

            if (float_cond(i,j) .eq. 0) then
              u_boundary_contr (i-2+iphi,j-2+jphi) = u_boundary_contr (i-2+iphi,j-2+jphi) + &
                .25 * beta(i,j) * dxdyh * uq * xquad(ilq) * xquad(jlq)
            endif

          endif

          if (vmask (i-2+iphi,j-2+jphi) .eq. 1) then


            v_boundary_contr (i-2+iphi,j-2+jphi) = v_boundary_contr (i-2+iphi,j-2+jphi) + &
              .25 * dxdyh * nu (i,j) * ( (uy+vx) * Phi(2*(2*(jphi-1)+iphi)-1,2*(jq-1)+iq) + &
                           (4*vy+2*ux) * Phi(2*(2*(jphi-1)+iphi),2*(jq-1)+iq))

            if (float_cond(i,j) .eq. 0) then
              v_boundary_contr (i-2+iphi,j-2+jphi) = v_boundary_contr (i-2+iphi,j-2+jphi) + &
                .25 * beta(i,j) * dxdyh * vq * xquad(ilq) * xquad(jlq)
            endif

          endif
        enddo ; enddo
      enddo ; enddo

      if (float_cond(i,j) .eq. 1) then
        Usubcontr = 0.0 ; Vsubcontr = 0.0 ; basel = G%bathyT(i,j)
        Ucell(:,:) = u_boundary_values(i-1:i,j-1:j) ; Vcell(:,:) = v_boundary_values(i-1:i,j-1:j)
        Hcell(:,:) = H_node(i-1:i,j-1:j)
        call CG_action_subgrid_basal_bilinear &
            (Phisub, Hcell, Ucell, Vcell, dxdyh, basel, dens_ratio, Usubcontr, Vsubcontr)
        do iphi=1,2 ; do jphi = 1,2
          if (umask (i-2+iphi,j-2+jphi) .eq. 1) then
            u_boundary_contr (i-2+iphi,j-2+jphi) = u_boundary_contr (i-2+iphi,j-2+jphi) + &
              Usubcontr(iphi,jphi) * beta (i,j)
          endif
          if (vmask (i-2+iphi,j-2+jphi) .eq. 1) then
            v_boundary_contr (i-2+iphi,j-2+jphi) = v_boundary_contr (i-2+iphi,j-2+jphi) + &
              Vsubcontr(iphi,jphi) * beta (i,j)
          endif
        enddo ; enddo
      endif
    endif
  endif ; enddo ; enddo

end subroutine apply_boundary_values_bilinear

subroutine calc_shelf_visc_triangular (CS,u,v)
  type(ice_shelf_CS),         pointer   :: CS
  real, dimension(:,:), intent(inout)    :: u, v

! update DEPTH_INTEGRATED viscosity, based on horizontal strain rates - this is for triangle FEM solve so there is
! an "upper" and "lower" triangular viscosity

! also this subroutine updates the nonlinear part of the basal traction

! this may be subject to change later... to make it "hybrid"

  real, pointer, dimension (:,:)    :: nu_lower , &
                         nu_upper, &
                       beta_eff_lower, &
                       beta_eff_upper
  real, pointer, dimension (:,:)    :: H,    &! thickness
                       hmask

  type(ocean_grid_type), pointer :: G
  integer :: isym, i, j, iscq, iecq, jscq, jecq, isd, jsd, ied, jed, iegq, jegq, giec, gjec, gisc, gjsc, cnt, isc, jsc, iec, jec, is, js
  real :: A, n, ux, uy, vx, vy, eps_min, umid, vmid, unorm, C_basal_friction, n_basal_friction, dxh, dyh, dxdyh

  G => CS%grid

  if (G%symmetric) then
     isym = 1
  else
     isym = 0
  endif

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  isd = G%isd ; jsd = G%jsd ; ied = G%isd ; jed = G%jsd
  iegq = G%iegB ; jegq = G%jegB
  gisc = G%domain%nihalo+1 ; gjsc = G%domain%njhalo+1
  giec = G%domain%niglobal+gisc ; gjec = G%domain%njglobal+gjsc
  is = iscq - (1-isym); js = jscq - (1-isym)

  A = CS%A_glen_isothermal ; n = CS%n_glen; eps_min = CS%eps_glen_min

  H => CS%h_shelf
  hmask => CS%hmask
  nu_upper => CS%ice_visc_upper_tri
  nu_lower => CS%ice_visc_lower_tri
  beta_eff_upper => CS%taub_beta_eff_upper_tri
  beta_eff_lower => CS%taub_beta_eff_lower_tri

  C_basal_friction = CS%C_basal_friction ; n_basal_friction = CS%n_basal_friction

  do i=isd,ied
    do j=jsd,jed

      dxh = G%dxT(i,j)
      dyh = G%dyT(i,j)
      dxdyh = G%areaT(i,j)

      if (hmask (i,j) .eq. 1) then
        ux = (u(i,j-1)-u(i-1,j-1)) / dxh
        vx = (v(i,j-1)-v(i-1,j-1)) / dxh
        uy = (u(i-1,j)-u(i-1,j-1)) / dyh
        vy = (v(i-1,j)-v(i-1,j-1)) / dyh

        nu_lower(i,j) = A**(-1/n) * (ux**2+vy**2+ux*vy+0.25*(uy+vx)**2+eps_min**2) ** ((1-n)/(2*n)) * H(i,j)
        umid = 1./3 * (u(i-1,j-1)+u(i-1,j)+u(i,j-1))
        vmid = 1./3 * (v(i-1,j-1)+v(i-1,j)+v(i,j-1))
        unorm = sqrt (umid**2+vmid**2+(eps_min*dxh)**2) ; beta_eff_lower (i,j) = C_basal_friction * unorm ** (n_basal_friction-1)

        ux = (u(i,j)-u(i-1,j)) / dxh
        vx = (v(i,j)-v(i-1,j)) / dxh
        uy = (u(i,j)-u(i,j-1)) / dyh
        vy = (u(i,j)-u(i,j-1)) / dyh

        nu_upper(i,j) = A**(-1/n) * (ux**2+vy**2+ux*vy+0.25*(uy+vx)**2+eps_min**2) ** ((1-n)/(2*n)) * H(i,j)
        umid = 1./3 * (u(i,j)+u(i-1,j)+u(i,j-1))
        vmid = 1./3 * (v(i,j)+v(i-1,j)+v(i,j-1))
        unorm = sqrt (umid**2+vmid**2+(eps_min*dxh)**2) ; beta_eff_upper (i,j) = C_basal_friction * unorm ** (n_basal_friction-1)

      endif
    enddo
  enddo

end subroutine calc_shelf_visc_triangular

subroutine calc_shelf_visc_bilinear (CS, u, v)
  type(ice_shelf_CS),         pointer   :: CS
  real, dimension(NILIMB_SYM_,NJLIMB_SYM_), intent(inout)    :: u, v

! update DEPTH_INTEGRATED viscosity, based on horizontal strain rates - this is for triangle FEM solve so there is
! an "upper" and "lower" triangular viscosity

! also this subroutine updates the nonlinear part of the basal traction

! this may be subject to change later... to make it "hybrid"

  real, pointer, dimension (:,:)    :: nu, &
                       beta
  real, pointer, dimension (:,:)    :: H,    &! thickness
                       hmask

  type(ocean_grid_type), pointer :: G
  integer :: isym, i, j, iscq, iecq, jscq, jecq, isd, jsd, ied, jed, iegq, jegq, giec, gjec, gisc, gjsc, cnt, isc, jsc, iec, jec, is, js
  real :: A, n, ux, uy, vx, vy, eps_min, umid, vmid, unorm, C_basal_friction, n_basal_friction, dxh, dyh, dxdyh

  G => CS%grid

  isym=0
  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
  iegq = G%iegB ; jegq = G%jegB
  gisc = G%domain%nihalo+1 ; gjsc = G%domain%njhalo+1
  giec = G%domain%niglobal+gisc ; gjec = G%domain%njglobal+gjsc
  is = iscq - (1-isym); js = jscq - (1-isym)

  A = CS%A_glen_isothermal ; n = CS%n_glen; eps_min = CS%eps_glen_min
  C_basal_friction = CS%C_basal_friction ; n_basal_friction = CS%n_basal_friction

  H => CS%h_shelf
  hmask => CS%hmask
  nu => CS%ice_visc_bilinear
  beta => CS%taub_beta_eff_bilinear

  do j=jsd+1,jed-1
    do i=isd+1,ied-1

      dxh = G%dxT(i,j)
      dyh = G%dyT(i,j)
      dxdyh = G%areaT(i,j)

      if (hmask (i,j) .eq. 1) then
        ux = (u(i,j) + u(i,j-1) - u(i-1,j) - u(i-1,j-1)) / (2*dxh)
        vx = (v(i,j) + v(i,j-1) - v(i-1,j) - v(i-1,j-1)) / (2*dxh)
        uy = (u(i,j) - u(i,j-1) + u(i-1,j) - u(i-1,j-1)) / (2*dyh)
        vy = (v(i,j) - v(i,j-1) + v(i-1,j) - v(i-1,j-1)) / (2*dyh)

        nu(i,j) = .5 * A**(-1/n) * (ux**2+vy**2+ux*vy+0.25*(uy+vx)**2+eps_min**2) ** ((1-n)/(2*n)) * H(i,j)

        umid = (u(i,j) + u(i,j-1) + u(i-1,j) + u(i-1,j-1))/4
        vmid = (v(i,j) + v(i,j-1) + v(i-1,j) + v(i-1,j-1))/4
        unorm = sqrt (umid**2+vmid**2+(eps_min*dxh)**2) ; beta (i,j) = C_basal_friction * unorm ** (n_basal_friction-1)
      endif
    enddo
  enddo

end subroutine calc_shelf_visc_bilinear

subroutine update_OD_ffrac (CS, ocean_mass, counter, nstep_velocity, time_step, velocity_update_time_step)
  type(ice_shelf_CS),         pointer   :: CS
  real, pointer, dimension(:,:)        :: ocean_mass
  integer,intent(in)            :: counter
  integer,intent(in)            :: nstep_velocity
  real,intent(in)            :: time_step
  real,intent(in)            :: velocity_update_time_step

  type(ocean_grid_type), pointer :: G
  integer :: isc, iec, jsc, jec, i, j
  real      :: threshold_col_depth, rho_ocean, inv_rho_ocean

  threshold_col_depth = CS%thresh_float_col_depth

  G=>CS%grid

  rho_ocean = CS%density_ocean_avg
  inv_rho_ocean = 1./rho_ocean

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  do j=jsc,jec
    do i=isc,iec
      CS%OD_rt(i,j) = CS%OD_rt(i,j) +  ocean_mass(i,j)*inv_rho_ocean
      if (ocean_mass(i,j) > threshold_col_depth*rho_ocean) then
        CS%float_frac_rt(i,j) = CS%float_frac_rt(i,j) + 1.0
      endif
    enddo
  enddo

  if (counter .eq. nstep_velocity) then

    do j=jsc,jec
      do i=isc,iec
        CS%float_frac(i,j) = 1.0 - (CS%float_frac_rt(i,j) / real(nstep_velocity))
!    if ((CS%float_frac(i,j) .gt. 0) .and. (CS%float_frac(i,j) .lt. 1)) then
!        print *,"PARTLY GROUNDED", CS%float_frac(i,j),i,j,mpp_pe()
!    endif
        CS%OD_av(i,j) = CS%OD_rt(i,j) / real(nstep_velocity)

        CS%OD_rt(i,j) = 0.0 ; CS%float_frac_rt(i,j) = 0.0
      enddo
    enddo

    call pass_var(CS%float_frac, G%domain)
    call pass_var(CS%OD_av, G%domain)

  endif

end subroutine update_OD_ffrac

subroutine update_OD_ffrac_uncoupled (CS)
  type(ice_shelf_CS), pointer    :: CS

  type(ocean_grid_type), pointer :: G
  integer             :: i, j, iters, isd, ied, jsd, jed
  real                 :: rhoi, rhow, OD
  type(time_type)          :: dummy_time
  real,dimension(:,:),pointer     :: OD_av, float_frac, h_shelf


  G => CS%grid
  rhoi = CS%density_ice
  rhow = CS%density_ocean_avg
  dummy_time = set_time (0,0)
  OD_av => CS%OD_av
  h_shelf => CS%h_shelf
  float_frac => CS%float_frac
  isd=G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

!   print *,"rhow",rhow,"rho",rhoi

  do j=jsd,jed
    do i=isd,ied
      OD = G%bathyT(i,j) - rhoi/rhow * h_shelf (i,j)
      if (OD.ge.0) then
    ! ice thickness does not take up whole ocean column -> floating
        OD_av (i,j) = OD
        float_frac(i,j) = 0.
      else
        OD_av (i,j) = 0.
        float_frac(i,j) = 1.
      endif
    enddo
  enddo


end subroutine update_OD_ffrac_uncoupled

subroutine bilinear_shape_functions (X, Y, Phi, area)
  real, dimension(4), intent(in) :: X, Y
  real, dimension(8,4), intent (inout) :: Phi
  real, intent (out) :: area

! X and Y must be passed in the form
    !  3 - 4
    !  |   |
    !  1 - 2

! this subroutine calculates the gradients of bilinear basis elements that
! that are centered at the vertices of the cell. values are calculated at
! points of gaussian quadrature. (in 1D: .5 * (1 +/- sqrt(1/3)) for [0,1])
!     (ordered in same way as vertices)
!
! Phi (2*i-1,j) gives d(Phi_i)/dx at quadrature point j
! Phi (2*i,j) gives d(Phi_i)/dy at quadrature point j
! Phi_i is equal to 1 at vertex i, and 0 at vertex k .ne. i, and bilinear
!
! This should be a one-off; once per nonlinear solve? once per lifetime?
! ... will all cells have the same shape and dimension?

  real, dimension (4) :: xquad, yquad
  integer :: node, qpoint, xnode, xq, ynode, yq
  real :: a,b,c,d,e,f,xexp,yexp

  xquad(1:3:2) = .5 * (1-sqrt(1./3)) ; yquad(1:2) = .5 * (1-sqrt(1./3))
  xquad(2:4:2) = .5 * (1+sqrt(1./3)) ; yquad(3:4) = .5 * (1+sqrt(1./3))

  do qpoint=1,4

    a = -X(1)*(1-yquad(qpoint)) + X(2)*(1-yquad(qpoint)) - X(3)*yquad(qpoint) + X(4)*yquad(qpoint) ! d(x)/d(x*)
    b = -Y(1)*(1-yquad(qpoint)) + Y(2)*(1-yquad(qpoint)) - Y(3)*yquad(qpoint) + Y(4)*yquad(qpoint) ! d(y)/d(x*)
    c = -X(1)*(1-xquad(qpoint)) - X(2)*(xquad(qpoint)) + X(3)*(1-xquad(qpoint)) + X(4)*(xquad(qpoint)) ! d(x)/d(y*)
    d = -Y(1)*(1-xquad(qpoint)) - Y(2)*(xquad(qpoint)) + Y(3)*(1-xquad(qpoint)) + Y(4)*(xquad(qpoint)) ! d(y)/d(y*)

    do node=1,4

      xnode = 2-mod(node,2) ; ynode = ceiling(REAL(node)/2)

      if (ynode .eq. 1) then
        yexp = 1-yquad(qpoint)
      else
        yexp = yquad(qpoint)
      endif

      if (1 .eq. xnode) then
        xexp = 1-xquad(qpoint)
      else
        xexp = xquad(qpoint)
      endif

      Phi (2*node-1,qpoint) = ( d * (2 * xnode - 3) * yexp - b * (2 * ynode - 3) * xexp) / (a*d-b*c)
      Phi (2*node,qpoint) = ( -c * (2 * xnode - 3) * yexp + a * (2 * ynode - 3) * xexp) / (a*d-b*c)

    enddo
  enddo

  area = quad_area (X,Y)

end subroutine bilinear_shape_functions


subroutine bilinear_shape_functions_subgrid (Phisub, nsub)
  real, dimension(nsub,nsub,2,2,2,2), intent(inout) :: Phisub
  integer                                           :: nsub

  ! this subroutine is a helper for interpolation of floatation condition
  ! for the purposes of evaluating the terms \int (u,v) \phi_i dx dy in a cell that is
  !     in partial floatation
  ! the array Phisub contains the values of \phi_i (where i is a node of the cell)
  !     at quad point j
  ! i think this general approach may not work for nonrectangular elements...
  !

  ! Phisub (i,j,k,l,q1,q2)
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

  do j=1,nsub
    do i=1,nsub
      x0 = (i-1) * fracx ; y0 = (j-1) * fracx
      do qx=1,2
        do qy=1,2
          x = x0 + fracx*xquad(qx)
          y = y0 + fracx*xquad(qy)
          do k=1,2
            do l=1,2
              val = 1.0
              if (k .eq. 1) then
                val = val * (1.0-x)
              else
                val = val * x
              endif
              if (l .eq. 1) then
                val = val * (1.0-y)
              else
                val = val * y
              endif
              Phisub (i,j,k,l,qx,qy) = val
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

!  print *, Phisub(1,1,2,2,1,1),Phisub(1,1,2,2,1,2),Phisub(1,1,2,2,2,1),Phisub(1,1,2,2,2,2)


end subroutine bilinear_shape_functions_subgrid


subroutine update_velocity_masks (CS)
  type(ice_shelf_CS),    pointer    :: CS

  ! sets masks for velocity solve
  ! ignores the fact that their might be ice-free cells - this only considers the computational boundary

  ! !!!!IMPORTANT!!!! relies on thickness mask - assumed that this is called after hmask has been updated (and halo-updated)

  integer :: isym, i, j, iscq, iecq, jscq, jecq, isd, jsd, is, js, iegq, jegq, giec, gjec, gisc, gjsc, isc, jsc, iec, jec, k
  integer :: i_off, j_off
  type(ocean_grid_type), pointer :: G
  real, dimension(:,:), pointer  :: umask, vmask, u_face_mask, v_face_mask, hmask, u_face_mask_boundary, v_face_mask_boundary

  G => CS%grid
  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  i_off = G%idg_offset ; j_off = G%jdg_offset
  isd = G%isd ; jsd = G%jsd
  iegq = G%iegB ; jegq = G%jegB
  gisc = G%Domain%nihalo ; gjsc = G%Domain%njhalo
  giec = G%Domain%niglobal+gisc ; gjec = G%Domain%njglobal+gjsc

  umask => CS%umask
  vmask => CS%vmask
  u_face_mask => CS%u_face_mask
  v_face_mask => CS%v_face_mask
  u_face_mask_boundary => CS%u_face_mask_boundary
  v_face_mask_boundary => CS%v_face_mask_boundary
  hmask => CS%hmask


!   if (G%symmetric) then
!     isym=1
!   else
!     isym=0
!   endif

  isym = 0

  umask (:,:) = 0 ; vmask (:,:) = 0
  u_face_mask (:,:) = 0 ; v_face_mask (:,:) = 0

  if (G%symmetric) then
   is = isd ; js = jsd
  else
   is = isd+1 ; js = jsd+1
  endif

  do j=js,G%jed
    do i=is,G%ied

      if (hmask(i,j) .eq. 1) then

        umask(i-1:i,j-1:j) = 1.
        vmask(i-1:i,j-1:j) = 1.

        do k=0,1

          select case (int(u_face_mask_boundary(i-1+k,j)))
            case (3)
              umask(i-1+k,j-1:j)=3.
              vmask(i-1+k,j-1:j)=0.
              u_face_mask(i-1+k,j)=3.
            case (2)
              u_face_mask(i-1+k,j)=2.
            case (4)
              umask(i-1+k,j-1:j)=0.
              vmask(i-1+k,j-1:j)=0.
              u_face_mask(i-1+k,j)=4.
            case (0)
              umask(i-1+k,j-1:j)=0.
              vmask(i-1+k,j-1:j)=0.
              u_face_mask(i-1+k,j)=0.
            case (1)  ! stress free x-boundary
              umask(i-1+k,j-1:j)=0.
            case default
          end select
        enddo

        do k=0,1

          select case (int(v_face_mask_boundary(i,j-1+k)))
            case (3)
              vmask(i-1:i,j-1+k)=3.
              umask(i-1:i,j-1+k)=0.
              v_face_mask(i,j-1+k)=3.
            case (2)
              v_face_mask(i,j-1+k)=2.
            case (4)
              umask(i-1:i,j-1+k)=0.
              vmask(i-1:i,j-1+k)=0.
              v_face_mask(i,j-1+k)=4.
            case (0)
              umask(i-1:i,j-1+k)=0.
              vmask(i-1:i,j-1+k)=0.
              u_face_mask(i,j-1+k)=0.
            case (1) ! stress free y-boundary
              vmask(i-1:i,j-1+k)=0.
            case default
          end select
        enddo

        !if (u_face_mask_boundary(i-1,j).geq.0) then !left boundary
        !  u_face_mask (i-1,j) = u_face_mask_boundary(i-1,j)
        !  umask (i-1,j-1:j) = 3.
        !  vmask (i-1,j-1:j) = 0.
        !endif

        !if (j_off+j .eq. gjsc+1) then !bot boundary
        !  v_face_mask (i,j-1) = 0.
        !  umask (i-1:i,j-1) = 0.
        !  vmask (i-1:i,j-1) = 0.
        !elseif (j_off+j .eq. gjec) then !top boundary
        !  v_face_mask (i,j) = 0.
        !  umask (i-1:i,j) = 0.
        !  vmask (i-1:i,j) = 0.
        !endif

        if (i .lt. G%ied) then
          if ((hmask(i+1,j) .eq. 0) &
              .OR. (hmask(i+1,j) .eq. 2)) then
            !right boundary or adjacent to unfilled cell
            u_face_mask (i,j) = 2.
          endif
        endif

        if (i .gt. G%isd) then
          if ((hmask(i-1,j) .eq. 0) .OR. (hmask(i-1,j) .eq. 2)) then
            !adjacent to unfilled cell
            u_face_mask (i-1,j) = 2.
          endif
        endif

        if (j .gt. G%jsd) then
          if ((hmask(i,j-1) .eq. 0) .OR. (hmask(i,j-1) .eq. 2)) then
            !adjacent to unfilled cell
            v_face_mask (i,j-1) = 2.
          endif
        endif

        if (j .lt. G%jed) then
          if ((hmask(i,j+1) .eq. 0) .OR. (hmask(i,j+1) .eq. 2)) then
            !adjacent to unfilled cell
            v_face_mask (i,j) = 2.
          endif
        endif


      endif

    enddo
  enddo

  ! note: if the grid is nonsymmetric, there is a part that will not be transferred with a halo update
  ! so this subroutine must update its own symmetric part of the halo

  call pass_vector(u_face_mask, v_face_mask, G%domain, TO_ALL, CGRID_NE)
  call pass_vector (umask,vmask,G%domain,TO_ALL,BGRID_NE)

end subroutine update_velocity_masks


subroutine interpolate_H_to_B (CS, h_shelf, hmask, H_node)
  type(ice_shelf_CS), pointer                            :: CS
  real, dimension (:,:), intent(in)                      :: h_shelf, hmask
  real, dimension (NILIMB_SYM_,NJLIMB_SYM_), intent(inout) :: H_node

  type(ocean_grid_type), pointer :: G
  integer                        :: i, j, isc, iec, jsc, jec, num_h, k, l
  real                           :: summ

  G => CS%grid
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
          if (hmask (i+k,j+l) .eq. 1.0) then
            summ = summ + h_shelf (i+k,j+l)
            num_h = num_h + 1
          endif
        enddo
      enddo
      if (num_h .gt. 0) then
        H_node(i,j) = summ / num_h
      endif
    enddo
  enddo

  call pass_var(H_node, G%domain)

end subroutine interpolate_H_to_B


subroutine ice_shelf_end(CS)
  type(ice_shelf_CS), pointer   :: CS
  ! This subroutine deallocates all memory associated with this module.

  if (.not.associated(CS)) return

  deallocate(CS%mass_shelf) ; deallocate(CS%area_shelf_h)
  deallocate(CS%t_flux) ; deallocate(CS%lprec)
  deallocate(CS%salt_flux)

  deallocate(CS%tflux_shelf) ; deallocate(CS%tfreeze);
  deallocate(CS%exch_vel_t) ; deallocate(CS%exch_vel_s)

  deallocate(CS%h_shelf) ; deallocate(CS%hmask)

  if (CS%shelf_mass_is_dynamic .and. .not.CS%override_shelf_movement) then
    deallocate(CS%u_shelf) ; deallocate(CS%v_shelf)
!!! OVS !!!
    deallocate(CS%t_shelf); deallocate(CS%tmask);
    deallocate(CS%t_boundary_values)
    deallocate(CS%u_boundary_values) ; deallocate(CS%v_boundary_values)
    deallocate(CS%ice_visc_bilinear)
    deallocate(CS%ice_visc_lower_tri) ; deallocate(CS%ice_visc_upper_tri)
    deallocate(CS%u_face_mask) ; deallocate(CS%v_face_mask)
    deallocate(CS%umask) ; deallocate(CS%vmask)

    deallocate(CS%taub_beta_eff_bilinear)
    deallocate(CS%taub_beta_eff_upper_tri)
    deallocate(CS%taub_beta_eff_lower_tri)
    deallocate(CS%OD_rt) ; deallocate(CS%OD_av)
    deallocate(CS%float_frac) ; deallocate(CS%float_frac_rt)
  endif

  deallocate(CS)

end subroutine ice_shelf_end

subroutine savearray2(fname,A,flag)

! print 2-D array to file

! this is here strictly for debug purposes

CHARACTER(*),intent(in)         :: fname
! This change is to allow the code to compile with the GNU compiler.
! DOUBLE PRECISION,DIMENSION(:,:),intent(in)      :: A
REAL, DIMENSION(:,:), intent(in)      :: A
LOGICAL               :: flag

INTEGER               :: M,N,i,j,iock,lh,FIN
CHARACTER(23000)         :: ln
CHARACTER(17)            :: sing
CHARACTER(9)            :: STR
CHARACTER(7)            :: FMT1

if (.NOT. flag) then
 return
endif

PRINT *,"WRITING ARRAY " // fname

FIN=7
M = size(A,1)
N = size(A,2)

OPEN(unit=fin,FILE=fname,STATUS='REPLACE',ACCESS='SEQUENTIAL',&
   ACTION='WRITE',IOSTAT=iock)

IF(M .gt. 1300) THEN
   WRITE(fin) 'SECOND DIMENSION TOO LARGE'
   CLOSE(fin)
   RETURN
END IF

DO i=1,M
  WRITE(ln,'(E17.9)') A(i,1)
  DO j=2,N
    WRITE(sing,'(E17.9)') A(i,j)
    ln = TRIM(ln) // ' ' // TRIM(sing)
  END DO


  IF(i.eq.1) THEN

   lh = LEN(TRIM(ln))

   FMT1 = '(A'

   SELECT CASE (lh)
     CASE(1:9)
    WRITE(FMT1(3:3),'(I1)') lh

     CASE(10:99)
       WRITE(FMT1(3:4),'(I2)') lh

     CASE(100:999)
       WRITE(FMT1(3:5),'(I3)') lh

     CASE(1000:9999)
       WRITE(FMT1(3:6),'(I4)') lh

   END SELECT

   FMT1 = TRIM(FMT1) // ')'

  END IF

  WRITE(UNIT=fin,IOSTAT=iock,FMT=TRIM(FMT1)) TRIM(ln)

  IF(iock .ne. 0) THEN
     PRINT*,iock
  END IF
END DO

CLOSE(FIN)

end subroutine savearray2


subroutine solo_time_step (CS, time_step, n, Time, min_time_step_in)
  type(ice_shelf_CS), pointer    :: CS
  real,intent(in)      :: time_step
  integer, intent(inout)      :: n
  type(time_type)      :: Time
  real,optional,intent(in)   :: min_time_step_in

  type(ocean_grid_type), pointer :: G
  integer          :: is, iec, js, jec, i, j, ki, kj, iters
  real             :: ratio, min_ratio, time_step_remain, local_u_max, &
               local_v_max, time_step_int, min_time_step,spy,dumtimeprint
  real, dimension(:,:), pointer  :: u_shelf, v_shelf, hmask, umask, vmask
  logical         :: flag
  type (time_type)      :: dummy
  character(2)         :: procnum
  character(4)         :: stepnum

  CS%velocity_update_sub_counter = CS%velocity_update_sub_counter + 1
  spy = 365 * 86400
  G => CS%grid
  u_shelf => CS%u_shelf
  v_shelf => CS%v_shelf
  hmask => CS%hmask
  umask => CS%umask
  vmask => CS%vmask
  time_step_remain = time_step
  if (.not. (present (min_time_step_in))) then
   min_time_step = 1000 ! i think this is in seconds - this would imply ice is moving at ~1 meter per second
  else
   min_time_step=min_time_step_in
  endif
  is = G%isc ; iec = G%iec ; js = G%jsc ; jec = G%jec

  ! NOTE: this relies on NE grid indexing
  ! dumtimeprint=time_type_to_real(Time)/spy
  if (is_root_pe()) print *, "TIME in ice shelf call, yrs: ", time_type_to_real(Time)/spy

  do while (time_step_remain .gt. 0.0)

  min_ratio = 1.0e16
  n=n+1
  do j=js,jec
    do i=is,iec

       local_u_max = 0 ; local_v_max = 0

       if (hmask (i,j) .eq. 1.0) then
         ! all 4 corners of the cell should have valid velocity values; otherwise something is wrong
        ! this is done by checking that umask and vmask are nonzero at all 4 corners
        do ki=1,2 ; do kj = 1,2

          local_u_max = max (local_u_max, abs(u_shelf(i-1+ki,j-1+kj)))
          local_v_max = max (local_v_max, abs(v_shelf(i-1+ki,j-1+kj)))

        enddo ; enddo

        ratio = min (G%areaT(i,j) / (local_u_max+1.0e-12), G%areaT(i,j) / (local_v_max+1.0e-12))
        min_ratio = min (min_ratio, ratio)

       endif
     enddo ! j loop
   enddo ! i loop

   ! solved velocities are in m/yr; we want m/s

   call mpp_min (min_ratio)

   time_step_int = min(CS%CFL_factor * min_ratio * (365*86400), time_step)

   if (time_step_int .lt. min_time_step) then
     call MOM_error (FATAL, "MOM_ice_shelf:solo_time_step: abnormally small timestep")
   else
     if (is_root_pe()) then
   write(*,*) "Ice model timestep: ", time_step_int, " seconds"
     endif
   endif

   if (time_step_int .ge. time_step_remain) then
     time_step_int = time_step_remain
     time_step_remain = 0.0
   else
     time_step_remain = time_step_remain - time_step_int
   endif

   write (stepnum,'(I4)') CS%velocity_update_sub_counter

   call ice_shelf_advect (CS, time_step_int, CS%lprec, Time)


   if (mpp_pe() .eq. 7) then
      call savearray2 ("hmask",CS%hmask,CS%write_output_to_file)
!!! OVS!!!
!      call savearray2 ("tshelf",CS%t_shelf,CS%write_output_to_file)
   endif


   do j=G%jsd,G%jed
     do i=G%isd,G%ied
       if ((CS%hmask(i,j) .eq. 1) .or. (CS%hmask(i,j) .eq. 2)) then
         CS%mass_shelf(i,j) = CS%h_shelf(i,j)*CS%density_ice
       endif
     enddo
   enddo

   ! if the last mini-timestep is a day or less, we cannot expect velocities to change by much.
   ! do not update them
   if (time_step_int .gt. 1000) then
     call update_velocity_masks (CS)

!     call savearray2 ("Umask"//"p"//trim(procnum)//"_"//trim(stepnum),CS%umask,CS%write_output_to_file)
!     call savearray2 ("Vmask"//"p"//trim(procnum)//"_"//trim(stepnum),CS%vmask,CS%write_output_to_file)

     call update_OD_ffrac_uncoupled (CS)
     call ice_shelf_solve_outer (CS, CS%u_shelf, CS%v_shelf, 1, iters, dummy)
   endif

!!! OVS!!!
   call ice_shelf_temp (CS, time_step_int, CS%lprec, Time)

  call enable_averaging(time_step,Time,CS%diag)
   if (CS%id_area_shelf_h > 0) call post_data(CS%id_area_shelf_h, CS%area_shelf_h, CS%diag)
   if (CS%id_col_thick > 0) call post_data(CS%id_col_thick, CS%OD_av, CS%diag)
   if (CS%id_h_shelf > 0) call post_data(CS%id_h_shelf,CS%h_shelf,CS%diag)
   if (CS%id_h_mask > 0) call post_data(CS%id_h_mask,CS%hmask,CS%diag)
   if (CS%id_u_mask > 0) call post_data(CS%id_u_mask,CS%umask,CS%diag)
   if (CS%id_v_mask > 0) call post_data(CS%id_v_mask,CS%vmask,CS%diag)
   if (CS%id_u_shelf > 0) call post_data(CS%id_u_shelf,CS%u_shelf,CS%diag)
   if (CS%id_v_shelf > 0) call post_data(CS%id_v_shelf,CS%v_shelf,CS%diag)
   if (CS%id_float_frac > 0) call post_data(CS%id_float_frac,CS%float_frac,CS%diag)
   if (CS%id_OD_av >0) call post_data(CS%id_OD_av,CS%OD_av,CS%diag)
   if (CS%id_float_frac_rt>0) call post_data(CS%id_float_frac_rt,CS%float_frac_rt,CS%diag)
!!! OVS!!!
!   if (CS%id_t_mask > 0)
   call post_data(CS%id_t_mask,CS%tmask,CS%diag)
!   if (CS%id_t_shelf > 0)
   call post_data(CS%id_t_shelf,CS%t_shelf,CS%diag)

  call disable_averaging(CS%diag)

 enddo

end subroutine solo_time_step

!!! OVS !!!
subroutine ice_shelf_temp(CS, time_step, melt_rate, Time)
  type(ice_shelf_CS),         pointer    :: CS
  real,                       intent(in) :: time_step
  real,pointer,dimension(:,:),intent(in) :: melt_rate
  type(time_type)             :: Time

! time_step: time step in sec
! melt_rate: basal melt rate in kg/m^2/s

! 5/23/12 OVS
! Arguments:
! CS - A structure containing the ice shelf state - including current velocities
! t0 - an array containing temperature at the beginning of the call
! t_after_uflux - an array containing the temperature after advection in u-direction
! t_after_vflux - similar
!
!    This subroutine takes the velocity (on the Bgrid) and timesteps (HT)_t = - div (uHT) + (adot Tsurd -bdot Tbot) once and then calculates T=HT/H
!
!    The flux overflows are included here. That is because they will be used to advect 3D scalars
!    into partial cells

  !
  ! flux_enter: this is to capture flow into partially covered cells; it gives the mass flux into a given
  ! cell across its boundaries.
  ! ###Perhaps flux_enter should be changed into u-face and v-face
  ! ###fluxes, which can then be used in halo updates, etc.
  !
  !   from left neighbor:   flux_enter (:,:,1)
  !   from right neighbor:  flux_enter (:,:,2)
  !   from bottom neighbor: flux_enter (:,:,3)
  !   from top neighbor:    flux_enter (:,:,4)
  !
  !  THESE ARE NOT CONSISTENT ==> FIND OUT WHAT YOU IMPLEMENTED

  ! flux_enter(isd:ied,jsd:jed,1:4): if cell is not ice-covered, gives flux of ice into cell from kth boundary
  !
  !   o--- (4) ---o
  !   |           |
  !  (1)         (2)
  !   |           |
  !   o--- (3) ---o
  !

  type(ocean_grid_type), pointer :: G
  real, dimension(size(CS%h_shelf,1),size(CS%h_shelf,2))   :: th_after_uflux, th_after_vflux, TH
  real, dimension(size(CS%h_shelf,1),size(CS%h_shelf,2),4) :: flux_enter
  integer                           :: isd, ied, jsd, jed, i, j, isc, iec, jsc, jec
  real                              :: rho, spy, t_bd, Tsurf, adot
  real, dimension(:,:), pointer     :: hmask, Tbot
  character(len=2)                  :: procnum

  hmask => CS%hmask
  G => CS%grid
  rho = CS%density_ice
  spy = 365 * 86400 ! seconds per year; is there a global constant for this?  No - it is dependent upon a calendar.

  adot = 0.1/spy ! for now adot and Tsurf are defined here adot=surf acc 0.1m/yr, Tsurf=-20oC, vary them later
  Tbot =>CS%Tfreeze
  Tsurf = -20.0

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  flux_enter (:,:,:) = 0.0

  th_after_uflux (:,:) = 0.0
  th_after_vflux (:,:) = 0.0
!   if (is_root_pe()) write(*,*) "ice_shelf_advect called"

  do j=jsd,jed
    do i=isd,ied
      t_bd = CS%t_boundary_values(i,j)
!      if (CS%hmask(i,j) .gt. 1) then
      if ((CS%hmask(i,j) .eq. 3) .or. (CS%hmask(i,j) .eq. -2)) then
          CS%t_shelf(i,j) = CS%t_boundary_values(i,j)
      endif
    enddo
  enddo

  do j=jsd,jed
    do i=isd,ied
        TH (i,j) = CS%t_shelf(i,j)*CS%h_shelf (i,j)
    enddo
  enddo


!  call enable_averaging(time_step,Time,CS%diag)
 ! call pass_var (h_after_uflux, G%domain)
!  if (CS%id_h_after_uflux > 0) call post_data(CS%id_h_after_uflux, h_after_uflux, CS%diag)
!  call disable_averaging(CS%diag)


!  call enable_averaging(time_step,Time,CS%diag)
!  call pass_var (h_after_vflux, G%domain)
!  if (CS%id_h_after_vflux > 0) call post_data(CS%id_h_after_vflux, h_after_vflux, CS%diag)
!  call disable_averaging(CS%diag)



  call ice_shelf_advect_temp_x (CS, time_step/spy, TH, th_after_uflux, flux_enter)
  call ice_shelf_advect_temp_y (CS, time_step/spy, th_after_uflux, th_after_vflux, flux_enter)

  do j=jsd,jed
    do i=isd,ied
!      if (CS%hmask(i,j) .eq. 1) then
      if (CS%h_shelf(i,j) .gt. 0.0) then
        CS%t_shelf (i,j) = th_after_vflux(i,j)/CS%h_shelf (i,j)
      else
          CS%t_shelf(i,j) = -10.0
      endif
    enddo
  enddo

  do j=jsd,jed
    do i=isd,ied
      t_bd = CS%t_boundary_values(i,j)
!      if (CS%hmask(i,j) .gt. 1) then
      if ((CS%hmask(i,j) .eq. 3) .or. (CS%hmask(i,j) .eq. -2)) then
          CS%t_shelf(i,j) = t_bd
!          CS%t_shelf(i,j) = -15.0
      endif
    enddo
  enddo

  do j=jsc,jec
    do i=isc,iec
      if ((CS%hmask(i,j) .eq. 1) .or. (CS%hmask(i,j) .eq. 2)) then
        if (CS%h_shelf(i,j) .gt. 0.0) then
!          CS%t_shelf (i,j) = CS%t_shelf (i,j) + time_step*(adot*Tsurf -melt_rate (i,j)*Tbot(i,j))/CS%h_shelf (i,j)
          CS%t_shelf (i,j) = CS%t_shelf (i,j) + time_step*(adot*Tsurf -3/spy*Tbot(i,j))/CS%h_shelf (i,j)
        else
          ! the ice is about to melt away
          ! in this case set thickness, area, and mask to zero
          ! NOTE: not mass conservative
          ! should maybe scale salt & heat flux for this cell

          CS%t_shelf(i,j) = -10.0
          CS%tmask(i,j) = 0.0
        endif
      endif
    enddo
  enddo

  call pass_var(CS%t_shelf, G%domain)
  call pass_var(CS%tmask, G%domain)

  if (CS%DEBUG) then
    call hchksum (CS%t_shelf, "temp after front", G%HI, haloshift=3)
  endif

end subroutine ice_shelf_temp


subroutine ice_shelf_advect_temp_x (CS, time_step, h0, h_after_uflux, flux_enter)
  type(ice_shelf_CS),         pointer    :: CS
  real,                       intent(in) :: time_step
  real, dimension(:,:), intent(in) :: h0
  real, dimension(:,:), intent(inout) :: h_after_uflux
  real, dimension(:,:,:), intent(inout) :: flux_enter

  ! use will be made of CS%hmask here - its value at the boundary will be zero, just like uncovered cells

  ! if there is an input bdry condition, the thickness there will be set in initialization

  ! flux_enter(isd:ied,jsd:jed,1:4): if cell is not ice-covered, gives flux of ice into cell from kth boundary
  !
  !   from left neighbor:   flux_enter (:,:,1)
  !   from right neighbor:  flux_enter (:,:,2)
  !   from bottom neighbor: flux_enter (:,:,3)
  !   from top neighbor:    flux_enter (:,:,4)
  !
  !        o--- (4) ---o
  !        |           |
  !       (1)         (2)
  !        |           |
  !        o--- (3) ---o
  !

  integer :: isym, i, j, is, ie, js, je, isd, ied, jsd, jed, gjed, gied
  integer :: i_off, j_off
  logical :: at_east_bdry, at_west_bdry, one_off_west_bdry, one_off_east_bdry
  type(ocean_grid_type), pointer :: G
  real, dimension(-2:2) :: stencil
  real, dimension(:,:), pointer  :: hmask, u_face_mask, u_flux_boundary_values,u_boundary_values,t_boundary
  real :: u_face, &  ! positive if out
      flux_diff_cell, phi, dxh, dyh, dxdyh

  character (len=1)        :: debug_str, procnum

!   if (CS%grid%symmetric) then
!     isym = 1
!   else
!     isym = 0
!   endif

  isym = 0

  G => CS%grid
  hmask => CS%hmask
  u_face_mask => CS%u_face_mask
  u_flux_boundary_values => CS%u_flux_boundary_values
  u_boundary_values => CS%u_shelf
!  h_boundaries => CS%h_shelf
  t_boundary => CS%t_boundary_values
  is = G%isc-2 ; ie = G%iec+2 ; js = G%jsc ; je = G%jec ; isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset

  do j=jsd+1,jed-1
    if (((j+j_off) .le. G%domain%njglobal+G%domain%njhalo) .AND. &
        ((j+j_off) .ge. G%domain%njhalo+1)) then ! based on mehmet's code - only if btw north & south boundaries

      stencil(:) = -1
!     if (i+i_off .eq. G%domain%nihalo+G%domain%nihalo)
      do i=is,ie

        if (((i+i_off) .le. G%domain%niglobal+G%domain%nihalo) .AND. &
             ((i+i_off) .ge. G%domain%nihalo+1)) then

          if (i+i_off .eq. G%domain%nihalo+1) then
            at_west_bdry=.true.
          else
            at_west_bdry=.false.
          endif

          if (i+i_off .eq. G%domain%niglobal+G%domain%nihalo) then
            at_east_bdry=.true.
          else
            at_east_bdry=.false.
          endif

          if (hmask(i,j) .eq. 1) then

            dxh = G%dxT(i,j) ; dyh = G%dyT(i,j) ; dxdyh = G%areaT(i,j)

            h_after_uflux(i,j) = h0(i,j)

            stencil(:) = h0(i-2:i+2,j)  ! fine as long has nx_halo >= 2

            flux_diff_cell = 0

            ! 1ST DO LEFT FACE

            if (u_face_mask (i-1,j) .eq. 4.) then

              flux_diff_cell = flux_diff_cell + dyh * time_step * u_flux_boundary_values (i-1,j) * &
                               t_boundary(i-1,j) / dxdyh
! assume no flux bc for temp
!               flux_diff_cell = flux_diff_cell + dyh * time_step * CS%u_shelf(i,j)*t_boundary(i-1,j) / dxdyh

            else

              ! get u-velocity at center of left face
              u_face = 0.5 * (CS%u_shelf(i-1,j-1) + CS%u_shelf(i-1,j))

  !            if (at_west_bdry .and. (i .eq. G%isc)) then
  !                print *, j, u_face, stencil(-1)
  !            endif

              if (u_face .gt. 0) then !flux is into cell - we need info from h(i-2), h(i-1) if available

              ! i may not cover all the cases.. but i cover the realistic ones

                if (at_west_bdry .AND. (hmask(i-1,j).eq.3)) then ! at western bdry but there is a thickness bdry condition,
                              ! and the stencil contains it
                  stencil (-1) = CS%t_boundary_values(i-1,j)*CS%h_shelf(i-1,j)
                  flux_diff_cell = flux_diff_cell + ABS(u_face) * dyh * time_step * stencil(-1) / dxdyh

                elseif (hmask(i-1,j) * hmask(i-2,j) .eq. 1) then  ! h(i-2) and h(i-1) are valid
                  phi = slope_limiter (stencil(-1)-stencil(-2), stencil(0)-stencil(-1))
                  flux_diff_cell = flux_diff_cell + ABS(u_face) * dyh* time_step / dxdyh * &
                           (stencil(-1) - phi * (stencil(-1)-stencil(0))/2)

                else                            ! h(i-1) is valid
                                    ! (o.w. flux would most likely be out of cell)
                                    !  but h(i-2) is not

                  flux_diff_cell = flux_diff_cell + ABS(u_face) * dyh * time_step / dxdyh * stencil(-1)

                endif

              elseif (u_face .lt. 0) then !flux is out of cell - we need info from h(i-1), h(i+1) if available
                if (hmask(i-1,j) * hmask(i+1,j) .eq. 1) then         ! h(i-1) and h(i+1) are both valid
                  phi = slope_limiter (stencil(0)-stencil(1), stencil(-1)-stencil(0))
                  flux_diff_cell = flux_diff_cell - ABS(u_face) * dyh * time_step / dxdyh * &
                             (stencil(0) - phi * (stencil(0)-stencil(-1))/2)

                else
                  flux_diff_cell = flux_diff_cell - ABS(u_face) * dyh * time_step / dxdyh * stencil(0)

                  if ((hmask(i-1,j) .eq. 0) .OR. (hmask(i-1,j) .eq. 2)) then
                    flux_enter(i-1,j,2) = ABS(u_face) * dyh * time_step * stencil(0)
                  endif
                endif
              endif
            endif

            ! NEXT DO RIGHT FACE

            ! get u-velocity at center of right face

            if (u_face_mask (i+1,j) .eq. 4.) then

              flux_diff_cell = flux_diff_cell + dyh * time_step * u_flux_boundary_values (i+1,j) *&
                               t_boundary(i+1,j)/ dxdyh
! assume no flux bc for temp
!               flux_diff_cell = flux_diff_cell + dyh * time_step *  CS%u_shelf(i,j)*t_boundary (i+1,j)/ dxdyh

            else

              u_face = 0.5 * (CS%u_shelf(i,j-1) + CS%u_shelf(i,j))

              if (u_face .lt. 0) then !flux is into cell - we need info from h(i+2), h(i+1) if available

                if (at_east_bdry .AND. (hmask(i+1,j).eq.3)) then ! at eastern bdry but there is a thickness bdry condition,
                                            ! and the stencil contains it

                  flux_diff_cell = flux_diff_cell + ABS(u_face) * dyh * time_step * stencil(1) / dxdyh

                elseif (hmask(i+1,j) * hmask(i+2,j) .eq. 1) then  ! h(i+2) and h(i+1) are valid

                  phi = slope_limiter (stencil(1)-stencil(2), stencil(0)-stencil(1))
                  flux_diff_cell = flux_diff_cell + ABS(u_face) * dyh * time_step / dxdyh * &
                      (stencil(1) - phi * (stencil(1)-stencil(0))/2)

                else                            ! h(i+1) is valid
                                            ! (o.w. flux would most likely be out of cell)
                                            !  but h(i+2) is not

                  flux_diff_cell = flux_diff_cell + ABS(u_face) * dyh * time_step / dxdyh * stencil(1)

                endif

              elseif (u_face .gt. 0) then !flux is out of cell - we need info from h(i-1), h(i+1) if available

                if (hmask(i-1,j) * hmask(i+1,j) .eq. 1) then         ! h(i-1) and h(i+1) are both valid

                  phi = slope_limiter (stencil(0)-stencil(-1), stencil(1)-stencil(0))
                  flux_diff_cell = flux_diff_cell - ABS(u_face) * dyh * time_step / dxdyh * &
                      (stencil(0) - phi * (stencil(0)-stencil(1))/2)

                else                            ! h(i+1) is valid
                                            ! (o.w. flux would most likely be out of cell)
                                            !  but h(i+2) is not

                  flux_diff_cell = flux_diff_cell - ABS(u_face) * dyh * time_step / dxdyh * stencil(0)

                  if ((hmask(i+1,j) .eq. 0) .OR. (hmask(i+1,j) .eq. 2)) then

                    flux_enter(i+1,j,1) = ABS(u_face) * dyh * time_step  * stencil(0)
                  endif

                endif

              endif

              h_after_uflux(i,j) = h_after_uflux(i,j) + flux_diff_cell

            endif

          elseif ((hmask(i,j) .eq. 0) .OR. (hmask(i,j) .eq. 2)) then

            if (at_west_bdry .AND. (hmask(i-1,j) .EQ. 3)) then
              u_face = 0.5 * (CS%u_shelf(i-1,j-1) + CS%u_shelf(i-1,j))
              flux_enter (i,j,1) = ABS(u_face) * G%dyT(i,j) * time_step * t_boundary(i-1,j)*CS%thickness_boundary_values(i+1,j)
            elseif (u_face_mask (i-1,j) .eq. 4.) then
              flux_enter (i,j,1) = G%dyT(i,j) * time_step * u_flux_boundary_values (i-1,j)*t_boundary(i-1,j)
!              flux_enter (i,j,1) = G%dyh(i,j) * time_step *  CS%u_shelf(i,j)*t_boundary (i-1,j)
! assume no flux bc for temp
            endif

            if (at_east_bdry .AND. (hmask(i+1,j) .EQ. 3)) then
              u_face = 0.5 * (CS%u_shelf(i,j-1) + CS%u_shelf(i,j))
              flux_enter(i,j,2) = ABS(u_face) * G%dyT(i,j) * time_step * t_boundary(i+1,j)*CS%thickness_boundary_values(i+1,j)
            elseif (u_face_mask (i+1,j) .eq. 4.) then
              flux_enter (i,j,2) = G%dyT(i,j) * time_step * u_flux_boundary_values (i+1,j) * t_boundary(i+1,j)
! assume no flux bc for temp
!              flux_enter (i,j,2) = G%dyh(i,j) * time_step *  CS%u_shelf(i,j)*t_boundary (i+1,j)
            endif

!            if ((i .eq. is) .AND. (hmask(i,j) .eq. 0) .AND. (hmask(i-1,j) .eq. 1)) then
              ! this is solely for the purposes of keeping the mask consistent while advancing the front without having
              ! to call pass_var - if cell is empty and cell to left is ice-covered then this cell will become partly covered

!              hmask(i,j) = 2
!            elseif ((i .eq. ie) .AND. (hmask(i,j) .eq. 0) .AND. (hmask(i+1,j) .eq. 1)) then
              ! this is solely for the purposes of keeping the mask consistent while advancing the front without having
              ! to call pass_var - if cell is empty and cell to left is ice-covered then this cell will become partly covered

!              hmask(i,j) = 2

!            endif

          endif

        endif

      enddo ! i loop

    endif

  enddo ! j loop

!  write (procnum,'(I1)') mpp_pe()

end subroutine ice_shelf_advect_temp_x

subroutine ice_shelf_advect_temp_y (CS, time_step, h_after_uflux, h_after_vflux, flux_enter)
  type(ice_shelf_CS),         pointer    :: CS
  real,                       intent(in) :: time_step
  real, dimension(:,:), intent(in) :: h_after_uflux
  real, dimension(:,:), intent(inout) :: h_after_vflux
  real, dimension(:,:,:), intent(inout) :: flux_enter

  ! use will be made of CS%hmask here - its value at the boundary will be zero, just like uncovered cells

  ! if there is an input bdry condition, the thickness there will be set in initialization

  ! flux_enter(isd:ied,jsd:jed,1:4): if cell is not ice-covered, gives flux of ice into cell from kth boundary
  !
  !   from left neighbor:   flux_enter (:,:,1)
  !   from right neighbor:  flux_enter (:,:,2)
  !   from bottom neighbor: flux_enter (:,:,3)
  !   from top neighbor:    flux_enter (:,:,4)
  !
  !        o--- (4) ---o
  !        |           |
  !       (1)         (2)
  !        |           |
  !        o--- (3) ---o
  !

  integer :: isym, i, j, is, ie, js, je, isd, ied, jsd, jed, gjed, gied
  integer :: i_off, j_off
  logical :: at_north_bdry, at_south_bdry, one_off_west_bdry, one_off_east_bdry
  type(ocean_grid_type), pointer :: G
  real, dimension(-2:2) :: stencil
  real, dimension(:,:), pointer  :: hmask, v_face_mask, v_flux_boundary_values,t_boundary,v_boundary_values
  real :: v_face, &  ! positive if out
      flux_diff_cell, phi, dxh, dyh, dxdyh
  character(len=1)        :: debug_str, procnum

!   if (CS%grid%symmetric) then
!     isym = 1
!   else
!     isym = 0
!   endif

  isym = 0

  G => CS%grid
  hmask => CS%hmask
  v_face_mask => CS%v_face_mask
  v_flux_boundary_values => CS%v_flux_boundary_values
  t_boundary => CS%t_boundary_values
  v_boundary_values => CS%v_shelf
  is = G%isc ; ie = G%iec ; js = G%jsc-1 ; je = G%jec+1 ; isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset

  do i=isd+2,ied-2
    if (((i+i_off) .le. G%domain%niglobal+G%domain%nihalo) .AND. &
       ((i+i_off) .ge. G%domain%nihalo+1)) then  ! based on mehmet's code - only if btw east & west boundaries

      stencil(:) = -1

      do j=js,je

        if (((j+j_off) .le. G%domain%njglobal+G%domain%njhalo) .AND. &
             ((j+j_off) .ge. G%domain%njhalo+1)) then

          if (j+j_off .eq. G%domain%njhalo+1) then
            at_south_bdry=.true.
          else
            at_south_bdry=.false.
          endif
          if (j+j_off .eq. G%domain%njglobal+G%domain%njhalo) then
            at_north_bdry=.true.
          else
            at_north_bdry=.false.
          endif

          if (hmask(i,j) .eq. 1) then
            dxh = G%dxT(i,j) ; dyh = G%dyT(i,j) ; dxdyh = G%areaT(i,j)
            h_after_vflux (i,j) = h_after_uflux (i,j)

            stencil (:) = h_after_uflux (i,j-2:j+2)  ! fine as long has ny_halo >= 2
            flux_diff_cell = 0

            ! 1ST DO south FACE

            if (v_face_mask (i,j-1) .eq. 4.) then

              flux_diff_cell = flux_diff_cell + dxh * time_step * v_flux_boundary_values (i,j-1) * t_boundary(i,j-1)/ dxdyh
! assume no flux bc for temp
!              flux_diff_cell = flux_diff_cell + dxh * time_step *  CS%v_shelf(i,j)*t_boundary (i,j-1) / dxdyh

            else

              ! get u-velocity at center of left face
              v_face = 0.5 * (CS%v_shelf(i-1,j-1) + CS%v_shelf(i,j-1))

              if (v_face .gt. 0) then !flux is into cell - we need info from h(j-2), h(j-1) if available

                ! i may not cover all the cases.. but i cover the realistic ones

                if (at_south_bdry .AND. (hmask(i,j-1).eq.3)) then ! at western bdry but there is a thickness bdry condition,
                                            ! and the stencil contains it
                  flux_diff_cell = flux_diff_cell + ABS(v_face) * dxh * time_step * stencil(-1) / dxdyh

                elseif (hmask(i,j-1) * hmask(i,j-2) .eq. 1) then  ! h(j-2) and h(j-1) are valid

                  phi = slope_limiter (stencil(-1)-stencil(-2), stencil(0)-stencil(-1))
                  flux_diff_cell = flux_diff_cell + ABS(v_face) * dxh * time_step / dxdyh * &
                      (stencil(-1) - phi * (stencil(-1)-stencil(0))/2)

                else     ! h(j-1) is valid
                         ! (o.w. flux would most likely be out of cell)
                         !  but h(j-2) is not
                  flux_diff_cell = flux_diff_cell + ABS(v_face) * dxh * time_step / dxdyh * stencil(-1)
                endif

              elseif (v_face .lt. 0) then !flux is out of cell - we need info from h(j-1), h(j+1) if available

                if (hmask(i,j-1) * hmask(i,j+1) .eq. 1) then  ! h(j-1) and h(j+1) are both valid
                  phi = slope_limiter (stencil(0)-stencil(1), stencil(-1)-stencil(0))
                  flux_diff_cell = flux_diff_cell - ABS(v_face) * dxh * time_step / dxdyh * &
                      (stencil(0) - phi * (stencil(0)-stencil(-1))/2)
                else
                  flux_diff_cell = flux_diff_cell - ABS(v_face) * dxh * time_step / dxdyh * stencil(0)

                  if ((hmask(i,j-1) .eq. 0) .OR. (hmask(i,j-1) .eq. 2)) then
                    flux_enter(i,j-1,4) = ABS(v_face) * dyh * time_step * stencil(0)
                  endif

                endif

              endif

            endif

            ! NEXT DO north FACE

            if (v_face_mask(i,j+1) .eq. 4.) then

              flux_diff_cell = flux_diff_cell + dxh * time_step * v_flux_boundary_values (i,j+1) *&
                               t_boundary(i,j+1)/ dxdyh
! assume no flux bc for temp
!              flux_diff_cell = flux_diff_cell + dxh * time_step *  CS%v_shelf(i,j)*t_boundary (i,j+1) / dxdyh

            else

            ! get u-velocity at center of right face
              v_face = 0.5 * (CS%v_shelf(i-1,j) + CS%v_shelf(i,j))

              if (v_face .lt. 0) then !flux is into cell - we need info from h(j+2), h(j+1) if available

                if (at_north_bdry .AND. (hmask(i,j+1).eq.3)) then ! at eastern bdry but there is a thickness bdry condition,
                                            ! and the stencil contains it
                  flux_diff_cell = flux_diff_cell + ABS(v_face) * dxh * time_step * stencil(1) / dxdyh
                elseif (hmask(i,j+1) * hmask(i,j+2) .eq. 1) then  ! h(j+2) and h(j+1) are valid
                  phi = slope_limiter (stencil(1)-stencil(2), stencil(0)-stencil(1))
                  flux_diff_cell = flux_diff_cell + ABS(v_face) * dxh * time_step / dxdyh * &
                      (stencil(1) - phi * (stencil(1)-stencil(0))/2)
                else     ! h(j+1) is valid
                         ! (o.w. flux would most likely be out of cell)
                         !  but h(j+2) is not
                  flux_diff_cell = flux_diff_cell + ABS(v_face) * dxh * time_step / dxdyh * stencil(1)
                endif

              elseif (v_face .gt. 0) then !flux is out of cell - we need info from h(j-1), h(j+1) if available

                if (hmask(i,j-1) * hmask(i,j+1) .eq. 1) then         ! h(j-1) and h(j+1) are both valid
                  phi = slope_limiter (stencil(0)-stencil(-1), stencil(1)-stencil(0))
                  flux_diff_cell = flux_diff_cell - ABS(v_face) * dxh * time_step / dxdyh * &
                      (stencil(0) - phi * (stencil(0)-stencil(1))/2)
                else   ! h(j+1) is valid
                       ! (o.w. flux would most likely be out of cell)
                       !  but h(j+2) is not
                  flux_diff_cell = flux_diff_cell - ABS(v_face) * dxh * time_step / dxdyh * stencil(0)
                  if ((hmask(i,j+1) .eq. 0) .OR. (hmask(i,j+1) .eq. 2)) then
                    flux_enter(i,j+1,3) = ABS(v_face) * dxh * time_step * stencil(0)
                  endif
                endif

              endif

            endif

            h_after_vflux (i,j) = h_after_vflux (i,j) + flux_diff_cell

          elseif ((hmask(i,j) .eq. 0) .OR. (hmask(i,j) .eq. 2)) then

            if (at_south_bdry .AND. (hmask(i,j-1) .EQ. 3)) then
              v_face = 0.5 * (CS%v_shelf(i-1,j-1) + CS%v_shelf(i,j-1))
              flux_enter (i,j,3) = ABS(v_face) * G%dxT(i,j) * time_step * t_boundary(i,j-1)*CS%thickness_boundary_values(i,j-1)
            elseif (v_face_mask(i,j-1) .eq. 4.) then
              flux_enter (i,j,3) = G%dxT(i,j) * time_step * v_flux_boundary_values (i,j-1)*t_boundary(i,j-1)
! assume no flux bc for temp
!              flux_enter (i,j,3) = G%dxh(i,j) * time_step *  CS%v_shelf(i,j)*t_boundary (i,j-1)

            endif

            if (at_north_bdry .AND. (hmask(i,j+1) .EQ. 3)) then
              v_face = 0.5 * (CS%v_shelf(i-1,j) + CS%v_shelf(i,j))
              flux_enter (i,j,4) = ABS(v_face) * G%dxT(i,j) * time_step * t_boundary(i,j+1)*CS%thickness_boundary_values(i,j+1)
            elseif (v_face_mask(i,j+1) .eq. 4.) then
              flux_enter (i,j,4) = G%dxT(i,j) * time_step * v_flux_boundary_values (i,j+1)*t_boundary(i,j+1)
! assume no flux bc for temp
!              flux_enter (i,j,4) = G%dxh(i,j) * time_step * CS%v_shelf(i,j)*t_boundary (i,j+1)
            endif

!            if ((j .eq. js) .AND. (hmask(i,j) .eq. 0) .AND. (hmask(i,j-1) .eq. 1)) then
 ! this is solely for the purposes of keeping the mask consistent while advancing the front without having
                ! to call pass_var - if cell is empty and cell to left is ice-covered then this cell will become partly covered
 !             hmask (i,j) = 2
 !           elseif ((j .eq. je) .AND. (hmask(i,j) .eq. 0) .AND. (hmask(i,j+1) .eq. 1)) then
                ! this is solely for the purposes of keeping the mask consistent while advancing the front without having
                ! to call pass_var - if cell is empty and cell to left is ice-covered then this cell will become partly covered
!              hmask (i,j) = 2
!            endif

          endif
        endif
      enddo ! j loop
    endif
  enddo ! i loop

  !write (procnum,'(I1)') mpp_pe()

end subroutine ice_shelf_advect_temp_y

end module MOM_ice_shelf
