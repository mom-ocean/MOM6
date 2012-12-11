module MOM
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
!*            The Modular Ocean Model, Version 6.0                     *
!*                               MOM                                   *
!*                                                                     *
!*  By Alistair Adcroft, Stephen Griffies, and Robert Hallberg         *
!*                                                                     *
!*  With software contributions from:                                  *
!*    Whit Anderson, Brian Arbic, Will Cooke, Anand Gnanadesikan,      *
!*    Matthew Harrison, Mehmet Ilicak, Laura Jackson, Jasmine John,    *
!*    Bonnie Samuels, Harper Simmons, Laurent White and Niki Zadeh     *
!*                                                                     *
!*  MOM ice-shelf code by Daniel Goldberg, Robert Hallberg             *
!*    Chris Little, and Olga Sergienko                                 *
!*                                                                     *
!*  This file was first relased in 2012.                               *
!*                                                                     *
!*    This program (MOM) simulates the ocean by numerically solving    *
!*  the hydrostatic primitive equations in generalized Lagrangian      *
!*  vertical coordinates, typically trackng stretched pressure (p*)    *
!*  surfaces or following isopycnals in the ocean's interior, and      *
!*  general orthogonal horizontal coordinates.  Unlike earlier verions *
!*  of MOM, in MOM6 these equations are horizontally discretized on an *
!*  Arakawa C-grid.  (It remains to be seen whether a B-grid dynamic   *
!*  core will be revived in MOM6 at a later date; for now applications *
!*  requiring a B-grid discretization should use MOM5.1.)  MOM6 offers *
!*  a range of options for the physical parameterizations, from those  *
!*  most appropriate to highly idealized models for geophysical fluid  *
!*  dynamics studies to a rich suite of processes appropriate for      *
!*  realistic ocean simulations.  The thermodynamic options typically  *
!*  use conservative temperature and preformed salinity as conservative*
!*  state variables and a full nonlinear equation of state, but there  *
!*  are also idealized adiabatic configurations of the model that use  *
!*  fixed density layers.  Version 6.0 of MOM continues in the long    *
!*  tradition of a commitment to climate-quality ocean simulations     *
!*  embodied in previous versions of MOM, even as it draws extensively *
!*  on the lessons learned in the development of the Generalized Ocean *
!*  Layered Dynamics (GOLD) ocean model, which was also primarily      *
!*  developed at NOAA/GFDL.  MOM has also benefitted tremendously from *
!*  the FMS infrastucture, which it utilizes and shares with other     *
!*  component models developed at GFDL.                                *
!*                                                                     *
!*    When run is isopycnal-coordinate mode, the uppermost few layers  *
!*  are often used to describe a bulk mixed layer, including the       *
!*  effects of penetrating shortwave radiation.  Either a split-       *
!*  explicit time stepping scheme or a non-split scheme may be used    *
!*  for the dynamics, while the time stepping may be split (and use    *
!*  different numbers of steps to cover the same interval) for the     *
!*  forcing, the thermodynamics, and for the dynamics.  Most of the    *
!*  numerics are second order accurate in space.  MOM can run with an  *
!*  absurdly thin minimum layer thickness. A variety of non-isopycnal  *
!*  vertical coordinate options are under development, but all exploit *
!*  the advantages of a Lagrangian vertical coordinate, as discussed   *
!*  in detail by Adcroft and Hallberg (Ocean Modelling, 2006).         *
!*                                                                     *
!*    Details of the numerics and physical parameterizations are       *
!*  provided in the appropriate source files.  All of the available    *
!*  options are selected at run-time by parsing the input files,       *
!*  usually MOM_input and MOM_override, and the options choices are    *
!*  then documented for each run in MOM_param_docs.                    *
!*                                                                     *
!*    MOM6 integrates the equations forward in time in three distinct  *
!*  phases.  In one phase, the dynamic equations for the velocities    *
!*  and layerthicknesses are advanced, capturing the propagation of    *
!*  external and internal inertia-gravity waves, Rossby waves, and     *
!*  other strictly adiabatic processes, including lateral stresses,    *
!*  vertical viscosity and momentum forcing, and interface height      *
!*  diffusion (commonly called Gent-McWilliams diffusion in depth-     *
!*  coordinate models).  In the second phase, all tracers are advected *
!*  and diffused along the layers.  The third phase applies diabatic   *
!*  processes, vertical mixing of water properties, and perhaps        *
!*  vertical remapping to cause the layers to track the desired        *
!*  vertical coordinate.                                               *
!*                                                                     *
!*    The present file (MOM.F90) orchestrates the main time stepping   *
!*  loops. One time integration option for the dynamics uses a split   *
!*  explicit time stepping scheme to rapidly step the barotropic       *
!*  pressure and velocity fields. The barotropic velocities are        *
!*  averaged over the baroclinic time step before they are used to     *
!*  advect thickness and determine the baroclinic accelerations.  As   *
!*  described in Hallberg and Adcroft (2009), a barotropic correction  *
!*  is applied to the time-mean layer velocities to ensure that the    *
!*  sum of the layer transports agrees with the time-mean barotropic   *
!*  transport, thereby ensuring that the estimates of the free surface *
!*  from the sum of the layer thinknesses agrees with the final free   *
!*  surface height as calculated by the barotropic solver.  The        *
!*  barotropic and baroclinic velocities are kept consistent by        *
!*  recalculating the barotropic velocities from the baroclinic        *
!*  transports each time step. This scheme is described in Hallberg,   *
!*  1997, J. Comp. Phys. 135, 54-65 and in Hallberg and Adcroft, 2009, *
!*  Ocean Modelling, 29, 15-26.                                        *
!*                                                                     *
!*    The other time integration options use non-split time stepping   *
!*  schemes based on the 3-step third order Runge-Kutta scheme         *
!*  described in Matsuno, 1966, J. Met. Soc. Japan, 44, 85-88, or on   *
!*  a two-step quasi-2nd order Runge-Kutta scheme.  These are much     *
!*  slower than the split time-stepping scheme, but they are useful    *
!*  for providing a more robust solution for debugging cases where the *
!*  more complicated split time-stepping scheme may be giving suspect  *
!*  solutions.                                                         *
!*                                                                     *
!*    There are a range of closure options available.  Horizontal      *
!*  velocities are subject to a combination of horizontal biharmonic   *
!*  and Laplacian friction (based on a stress tensor formalism) and a  *
!*  vertical Fickian viscosity (perhaps using the kinematic viscosity  *
!*  of water).  The horizontal viscosities may be constant, spatially  *
!*  varying or may be dynamically calculated using Smagorinsky's       *
!*  approach.  A diapycnal diffusion of density and thermodynamic      *
!*  quantities is also allowed, but not required, as is horizontal     *
!*  diffusion of interface heights (akin to the Gent-McWilliams        *
!*  closure of geopotential coordinate models).  The diapycnal mixing  *
!*  may use a fixed diffusivity or it may use the shear Richardson     *
!*  number dependent closure, like that described in Jackson et al.    *
!*  (JPO, 2008).  When there is diapycnal diffusion, it applies to     *
!*  momentum as well. As this is in addition to the vertical viscosity,*
!*  the vertical Prandtl always exceeds 1.  A refined bulk-mixed layer *
!*  is often used to describe the planetary boundary layer in realistic*
!*  ocean simulations.                                                 *
!*                                                                     *
!*    MOM has a number of noteworthy debugging capabilities.           *
!*  Excessively large velocities are truncated and MOM will stop       *
!*  itself after a number of such instances to keep the model from     *
!*  crashing altogether.  This is useful in diagnosing failures,       *
!*  or (by accepting some truncations) it may be useful for getting    *
!*  the model past the adjustment from an ill-balanced initial         *
!*  condition.  In addition, all of the accelerations in the columns   *
!*  with excessively large velocities may be directed to a text file.  *
!*  Parallelization errors may be diagnosed using the DEBUG option,    *
!*  which causes extensive checksums to be written out along with      *
!*  comments indicating where in the algorithm the sums originate and  *
!*  what variable is being summed.  The point where these checksums    *
!*  differ between runs is usually a good indication of where in the   *
!*  code the problem lies.  All of the test cases provided with MOM    *
!*  are routinely tested to ensure that they give bitwise identical    *
!*  results regardless of the domain decomposition, or whether they    *
!*  use static or dynamic memory allocation.                           *
!*                                                                     *
!*    About 115 other files of source code and 4 header files comprise *
!*  the MOM code, although there are several hundred more files that   *
!*  make up the FMS infrastructure upon which MOM is built.  Each of   *
!*  the MOM files contains comments documenting what it does, and      *
!*  most of the file names are fairly self-evident. In addition, all   *
!*  subroutines and data types are referenced via a module use, only   *
!*  statement, and the module names are consistent with the file names,*
!*  so it is not too hard to find the source file for a subroutine.    *
!*                                                                     *
!*    The typical MOM directory tree is as follows:                    *
!*        ../MOM                                                       *
!*        |-- config_src                                               *
!*        |   |-- coupled_driver                                       *
!*        |   |-- dynamic                                              *
!*        |   `-- solo_driver                                          *
!*        |-- examples                                                 *
!*        |   |-- CM2G                                                 *
!*        |   |-- ...                                                  *
!*        |   `-- torus_advection_test                                 *
!*        `-- src                                                      *
!*            |-- core                                                 *
!*            |-- diagnostics                                          *
!*            |-- equation_of_state                                    *
!*            |-- framework                                            *
!*            |-- ice_shelf                                            *
!*            |-- initialization                                       *
!*            |-- parameterizations                                    *
!*            |   |-- lateral                                          *
!*            |   `-- vertical                                         *
!*            |-- tracer                                               *
!*            `-- user                                                 *
!*  Rather than describing each file here, each directory's contents   *
!*  will be described to give a broad overview of the MOM code         *
!*  structure.                                                         *
!*                                                                     *
!*    The directories under config_src contain files that are used for *
!*  configuring the code, for instance for coupled or ocean-only runs. *
!*  Only one or two of these directories are used in compiling any,    *
!*  particular run.                                                    *
!*                                                                     *
!*  config_src/coupled_driver:                                         *
!*    The files here are used to couple MOM as a component in a larger *
!*    run driven by the FMS coupler.  This includes code that converts *
!*    various forcing fields into the code structures and flux and unit*
!*    conventions used by MOM, and converts the MOM surface fields     *
!*    back to the forms used by other FMS components.                  *
!*  config_src/dynamic:                                                *
!*    The only file here is the version of MOM_memory.h that is used   *
!*    for dynamic memory configurations of MOM.                        *
!*  config_src/solo_driver:                                            *
!*    The files here are include the _main driver that is used when    *
!*    MOM is configured as an ocean-only model, as well as the files   *
!*    that specify the surface forcing in this configuration.          *
!*                                                                     *
!*    The directories under examples provide a large number of working *
!*  configurations of MOM, along with reference solutions for several  *
!*  different compilers on GFDL's latest large computer.  The versions *
!*  of MOM_memory.h in these directories need not be used if dynamic   *
!*  memory allocation is desired, and the answers should be unchanged. *
!*                                                                     *
!*    The directories under src contain most of the MOM files.  These  *
!*  files are used in every configuration using MOM.                   *
!*                                                                     *
!*  src/core:                                                          *
!*    The files here constitute the MOM dynamic core.  This directory  *
!*    also includes files with the types that describe the model's     *
!*    lateral grid and have defined types that are shared across       *
!*    various MOM modules to allow for more succinct and flexible      *
!*    subroutine argument lists.                                       *
!*  src/diagnostics:                                                   *
!*    The files here calculate various diagnostics that are anciliary  *
!*    to the model itself.  While most of these diagnostics do not     *
!*    directly affect the model's solution, there are some, like the   *
!*    calculation of the deformation radius, that are used in some     *
!*    of the process parameterizations.                                *
!*  src/equation_of_state:                                             *
!*    These files describe the physical properties of sea-water,       *
!*    including both the equation of state and when it freezes.        *
!*  src/framework:                                                     *
!*    These files provide infrastructure utilities for MOM.  Many are  *
!*    simply wrappers for capabilities provided by FMS, although others*
!*    provide capabilities (like the file_parser) that are unique to   *
!*    MOM.  When MOM is adapted to use a different modeling            *
!*    infrastructure, most of the required changes are in this         *
!*    directory.                                                       *
!*  src/initialization:                                                *
!*    These are the files that are used to initialize the MOM grid     *
!*    or provide the initial physical state for MOM.  These files are  *
!*    not intended to be modified, but provide a means for calling     *
!*    user-specific initialization code like the examples in src/user. *
!*  src/parameterizations/lateral:                                     *
!*    These files implement a number of quasi-lateral (along-layer)    *
!*    process parameterizations, including lateral viscosities,        *
!*    parameterizations of eddy effects, and the calculation of tidal  *
!*    forcing.                                                         *
!*  src/parameterizations/vertical:                                    *
!*    These files implement a number of vertical mixing or diabatic    *
!*    processes, including the effects of vertical viscosity and       *
!*    code to parameterize the planetary boundary layer.  There is a   *
!*    separate driver that orchestrates this portion of the algorithm, *
!*    and there is a diversity of parameterizations to be found here.  *
!*  src/tracer:                                                        *
!*    These files handle the lateral transport and diffusion of        *
!*    tracers, or are the code to implement various passive tracer     *
!*    packages.  Additional tracer packages are readily accomodated.   *
!*  src/user:                                                          *
!*    These are either stub routines that a user could use to change   *
!*    the model's initial conditions or forcing, or are examples that  *
!*    implement specific test cases.  These files can easily  be hand  *
!*    edited to create new analytically specified configurations.      *
!*                                                                     *
!*                                                                     *
!*    Most simulations can be set up by modifying only the files       *
!*  MOM_input, and possibly one or two of the files in src/user.       *
!*  In addition, the diag_table (MOM_diag_table) will commonly be      *
!*  modified to tailor the output to the needs of the question at      *
!*  hand.  The FMS utility mkmf works with a file called path_names    *
!*  to build an appropriate makefile, and path_names should be edited  *
!*  to reflect the actual location of the desired source code.         *
!*                                                                     *
!*                                                                     *
!*    There are 3 publicly visible subroutines in this file (MOM.F90). *
!*  step_MOM steps MOM over a specified interval of time.              *
!*  MOM_initialize calls initialize and does other initialization      *
!*    that does not warrant user modification.                         *
!*  calculate_surface_state determines the surface (mixed layer)       *
!*    properties of the current model state and packages pointers      *
!*    to these fields into an exported structure.                      *
!*                                                                     *
!*    The remaining subroutines in this file (src/core/MOM.F90) are:   *
!*  find_total_transport determines the barotropic mass transport.     *
!*  register_diags registers many diagnostic fields for the dynamic    *
!*    solver, or of the main model variables.                          *
!*  MOM_timing_init initializes various CPU time clocks.               *
!*  write_static_fields writes out various time-invariant fields.      *
!*  set_restart_fields is used to specify those fields that are        *
!*    written to and read from the restart file.                       *
!*  smooth_SSH applies filtering of the sea-surface height field that  *
!*    is reported back to the calling routine; it is rarely used.      *
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
  surface, &      ! A structure containing pointers to various fields which
                  ! may be used describe the surface state of MOM, and
                  ! which will be returned to the calling program
  thermo_var_ptrs, & ! A structure containing pointers to an assortment of
                  ! thermodynamic fields that may be available, including
                  ! potential temperature, salinity and mixed layer density.
  ocean_internal_state  ! A structure containing pointers to most of the above.

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
use MOM_diabatic_driver, only : adiabatic, adiabatic_driver_init, diabatic_driver_end
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
use MOM_sponge, only : init_sponge_diags, sponge_CS
use MOM_CS_type, only : MOM_control_struct, MOM_thermo_chksum, MOM_state_chksum
use MOM_dynamics_unsplit, only : step_MOM_dyn_unsplit, register_restarts_dyn_unsplit
use MOM_dynamics_unsplit, only : initialize_dyn_unsplit
use MOM_dynamics_split_RK2, only : step_MOM_dyn_split_RK2, register_restarts_dyn_split_RK2
use MOM_dynamics_split_RK2, only : initialize_dyn_split_RK2
use MOM_dynamics_unsplit_RK2, only : step_MOM_dyn_unsplit_RK2, register_restarts_dyn_unsplit_RK2
use MOM_dynamics_unsplit_RK2, only : initialize_dyn_unsplit_RK2
use MOM_regridding, only : initialize_regridding, end_regridding, regridding_main

implicit none ; private

#include <MOM_memory.h>

public initialize_MOM, step_MOM, MOM_end, calculate_surface_state
public MOM_control_struct ! This is exported for MOM_driver.F90 to use. We could
                           ! instead have MOM_driver import MOM_cs_type directly?

integer :: id_clock_ocean, id_clock_dynamics, id_clock_thermo
integer :: id_clock_tracer, id_clock_diabatic
integer :: id_clock_continuity ! Also in dynamics s/r
integer :: id_clock_thick_diff ! Also in dynamics s/r
integer :: id_clock_diagnostics, id_clock_Z_diag
integer :: id_clock_init, id_clock_MOM_init
integer :: id_clock_pass, id_clock_pass_init ! Also in dynamics d/r

contains

function step_MOM(fluxes, state, Time_start, time_interval, CS)
!   This subroutine orchestrates the time stepping of MOM.  The adiabatic
! dynamics are stepped by one of the calls to one of the step_MOM_dyn_...
! routines.  The action of lateral processes on tracers occur in the calls to
! advect_tracer and tracer_hordiff.  The vertical mixing and possibly remapping
! occur inside of diabatic.
  integer                              :: step_MOM
  type(forcing), intent(inout)         :: fluxes
  type(surface), intent(inout)         :: state
  type(time_type), intent(in)          :: Time_start
  real,           intent(in)           :: time_interval
  type(MOM_control_struct), pointer    :: CS
! Arguments: fluxes - An intent in structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (out)     state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (in)      Time_start - The starting time of a run segment, as a time type.
!  (in)      time_interval - The interval of time over which to integrate in s.
!  (in)      CS - The control structure returned by a previous call to
!                 initialize_MOM.
  type(ocean_grid_type), pointer :: grid ! A pointer to a structure containing
                                  ! metrics and related information.
  integer, save :: nt = 1 ! The running number of iterations.
  integer :: ntstep ! The number of time steps between tracer updates
                    ! or diabatic forcing.
  integer :: n_max  ! The number of steps to take in this call.
  integer :: m = 1  ! The current time level (1, 2, or 3).
  integer :: mp     ! The previous value of m.
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, n
  integer :: isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq
  real :: dt        ! The baroclinic time step in s.
  real :: dtth      ! The time step used for thickness diffusion, in s.
  real :: dtnt      ! The elapsed time since updating the tracers and applying
                    ! diabatic processes, in s.
  real :: dtbt_reset_time ! The value of CS%rel_time when DTBT was last
                    ! calculated, in s.
  real :: wt_end, wt_beg
  real, dimension(SZI_(CS%grid),SZJ_(CS%grid)) :: &
    eta_av, &       ! The average sea surface height or column mass over
                    ! a time step, in m or kg m-2.
    ssh             ! The sea surface height based on eta_av, in m.
  real, allocatable, dimension(:,:) :: &
    sfc_speed, &    ! The sea surface speed at h-points, in m s-1.
    frazil_ave, &   ! The average heat frazil heat flux
                    ! required to keep the temperature above freezing, in W m-2.
    salt_deficit_ave, &  ! The average salt flux required to keep the
                    ! salinity above 0.01 PSU, in gSalt m-2 s-1.
    Heat_PmE_ave, & !   The average effective heat flux into the ocean due to
                    ! the exchange of water with other components, times the
                    ! heat capacity of water, in W m-2.   
    intern_heat_ave !   The average heat flux into the ocean from geothermal or
                    ! other internal heat sources, in W m-2.   
  real, pointer, dimension(:,:,:,:) :: &
    u, &                     ! u : Zonal velocity, in m s-1.
    v, &                     ! v : Meridional velocity, in m s-1.
    h                        ! h : Layer thickness, in m.
  real, dimension(SZI_(CS%grid),SZJ_(CS%grid),SZK_(CS%grid)+1) :: eta_predia
  real :: tot_wt_ssh, Itot_wt_ssh, I_time_int
  type(time_type) :: Time_local
  integer :: pid_tau, pid_ustar, pid_psurf, pid_u, pid_h
  integer :: pid_T, pid_S

  grid => CS%grid
  is = grid%isc ; ie = grid%iec ; js = grid%jsc ; je = grid%jec ; nz = grid%ke
  Isq = grid%Iscq ; Ieq = grid%Iecq ; Jsq = grid%Jscq ; Jeq = grid%Jecq
  isd = grid%isd ; ied = grid%ied ; jsd = grid%jsd ; jed = grid%jed
  Isdq = grid%Isdq ; Iedq = grid%Iedq ; Jsdq = grid%Jsdq ; Jedq = grid%Jedq
  u => CS%u ; v => CS%v ; h => CS%h

  call cpu_clock_begin(id_clock_ocean)
 !   First determine the time step that is consistent with this call.
 ! It is anticipated that the time step will almost always coincide
 ! with dt.  In addition, ntstep is determined, subject to the constraint
 ! that ntstep cannot exceed n_max.
  if (time_interval <= CS%dt) then
    n_max = 1
  else
    n_max = ceiling(time_interval/CS%dt - 0.001)
  endif

  dt = real(time_interval) / real(n_max)
  dtnt = 0.0
  ntstep = MAX(1,MIN(n_max,floor(CS%dt_therm/dt + 0.001)))

  CS%calc_bbl = .true.
  if (.not.ASSOCIATED(fluxes%p_surf)) CS%interp_p_surf = .false.

  call cpu_clock_begin(id_clock_pass)
  if (grid%nonblocking_updates) then
    pid_tau = pass_vector_start(fluxes%taux, fluxes%tauy, grid%Domain)
    if (ASSOCIATED(fluxes%ustar)) &
      pid_ustar = pass_var_start(fluxes%ustar(:,:), grid%Domain)
    if (ASSOCIATED(fluxes%p_surf)) &
      pid_psurf = pass_var_start(fluxes%p_surf(:,:), grid%Domain)
  else
    call pass_vector(fluxes%taux, fluxes%tauy, grid%Domain)
    if (ASSOCIATED(fluxes%ustar)) call pass_var(fluxes%ustar(:,:), grid%Domain)
    if (ASSOCIATED(fluxes%p_surf)) call pass_var(fluxes%p_surf(:,:), grid%Domain)
  endif
  call cpu_clock_end(id_clock_pass)
  if (ASSOCIATED(CS%tv%frazil)) CS%tv%frazil(:,:) = 0.0
  if (ASSOCIATED(CS%tv%salt_deficit)) CS%tv%salt_deficit(:,:) = 0.0   
  if (ASSOCIATED(CS%tv%TempxPmE)) CS%tv%TempxPmE(:,:) = 0.0
  if (ASSOCIATED(CS%tv%internal_heat)) CS%tv%internal_heat(:,:) = 0.0

  CS%rel_time = 0.0

  tot_wt_ssh = 0.0
  do j=js,je ; do i=is,ie ; CS%ave_ssh(i,j) = 0.0 ; enddo ; enddo

  if (CS%interp_p_surf) then
    if (.not.ASSOCIATED(CS%p_surf_end)) allocate(CS%p_surf_end(isd:ied,jsd:jed))
    if (.not.ASSOCIATED(CS%p_surf_begin)) allocate(CS%p_surf_begin(isd:ied,jsd:jed))
    
    if (.not.CS%p_surf_prev_set) then
      do j=jsd,jed ; do i=isd,ied
        CS%p_surf_prev(i,j) = fluxes%p_surf(i,j)
      enddo ; enddo
      CS%p_surf_prev_set = .true.
    endif
  else
    CS%p_surf_end  => fluxes%p_surf
  endif

  mp = 1 ; if (CS%split) mp = MOD(nt+1,2) + 1
  if (CS%debug) then
    call MOM_state_chksum("Before steps ", u(:,:,:,mp), v(:,:,:,mp), &
                          h(:,:,:,mp), CS%uh, CS%vh, grid)
    call check_redundant("Before steps mp ", u(:,:,:,mp), v(:,:,:,mp), grid)
  endif

  if (associated(CS%VarMix)) then
    call enable_averaging(time_interval, Time_start+set_time(int(time_interval)), &
                          CS%diag)
    call calc_resoln_function(h(:,:,:,mp), CS%tv, grid, CS%VarMix)
    call disable_averaging(CS%diag)
  endif

  if (grid%nonblocking_updates) then
    call cpu_clock_begin(id_clock_pass)
    call pass_vector_complete(pid_tau, fluxes%taux, fluxes%tauy, grid%Domain)
    if (ASSOCIATED(fluxes%ustar)) &
      call pass_var_complete(pid_ustar, fluxes%ustar(:,:), grid%Domain)
    if (ASSOCIATED(fluxes%p_surf)) &
      call pass_var_complete(pid_psurf, fluxes%p_surf(:,:), grid%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  do n=1,n_max
    nt = nt + 1
    ! Set the universally visible time to the middle of the time step
    CS%Time = Time_start + set_time(int(floor(CS%rel_time+0.5*dt+0.5)))
    CS%rel_time = CS%rel_time + dt
    ! Set the local time to the end of the time step.
    Time_local = Time_start + set_time(int(floor(CS%rel_time+0.5)))

    call cpu_clock_begin(id_clock_dynamics)
    call disable_averaging(CS%diag)

    if (CS%thickness_diffuse .and. CS%thickness_diffuse_first) then
      if (MOD(n-1,ntstep) == 0) then
        dtth = dt*min(ntstep,n_max-n+1)
        call enable_averaging(dtth,Time_local+set_time(int(floor(dtth-dt+0.5))), CS%diag)
        mp = 1 ; if (CS%split) mp = MOD(nt,2) + 1
        call cpu_clock_begin(id_clock_thick_diff)
        if (associated(CS%VarMix)) &
          call calc_slope_function(h(:,:,:,mp), CS%tv, grid, CS%VarMix)
        call thickness_diffuse(h(:,:,:,mp), CS%uhtr, CS%vhtr, CS%tv, dtth, grid, &
                               CS%MEKE, CS%VarMix, CS%thickness_diffuse_CSp)
        call cpu_clock_end(id_clock_thick_diff)
        call cpu_clock_begin(id_clock_pass)
        call pass_var(h(:,:,:,mp), grid%Domain)
        call cpu_clock_end(id_clock_pass)
        call disable_averaging(CS%diag)
      endif
    endif

    if (CS%interp_p_surf) then
      wt_end = real(n) / real(n_max)
      wt_beg = real(n-1) / real(n_max)
      do j=jsd,jed ; do i=isd,ied
        CS%p_surf_end(i,j) = wt_end * fluxes%p_surf(i,j) + &
                        (1.0-wt_end) * CS%p_surf_prev(i,j)
        CS%p_surf_begin(i,j) = wt_beg * fluxes%p_surf(i,j) + &
                        (1.0-wt_beg) * CS%p_surf_prev(i,j)
      enddo ; enddo
    endif

    if (CS%calc_bbl) &
      CS%bbl_calc_time_interval = dt*real(1+MIN(ntstep-MOD(n,ntstep),n_max-n))

    if (associated(CS%u_prev) .and. associated(CS%v_prev)) then
      mp = 1 ; if (CS%split) mp = MOD(nt,2) + 1
      do k=1,nz ; do j=jsd,jed ; do I=Isdq,Iedq
        CS%u_prev(I,j,k) = u(I,j,k,mp)
      enddo ; enddo ; enddo
      do k=1,nz ; do J=Jsdq,Jedq ; do i=isd,ied
        CS%v_prev(I,j,k) = u(I,j,k,mp)
      enddo ; enddo ; enddo
    endif

    if (CS%split) then !--------------------------- start SPLIT
!   This section uses a predictor corrector scheme, that is somewhere
! (determined by be) between the forward-backward (be=0.5) scheme and
! the backward Euler scheme (be=1.0) to time step the dynamic equations.
      mp = MOD(nt,2) + 1
      m  = 3 - mp

      CS%calc_dtbt = .false.
      if ((CS%dtbt_reset_period >= 0.0) .and. &
          ((n==1) .or. (CS%dtbt_reset_period == 0.0) .or. &
           (CS%rel_time >= dtbt_reset_time + 0.999*CS%dtbt_reset_period))) then
        CS%calc_dtbt = .true.
        dtbt_reset_time = CS%rel_time
      endif

      call step_MOM_dyn_split_RK2(u(:,:,:,mp), v(:,:,:,mp), h(:,:,:,mp), &
                    CS%eta, CS%uhbt_in, CS%vhbt_in, Time_local, dt, &
                    fluxes, CS%p_surf_begin, CS%p_surf_end, dtnt, dt*ntstep, &
                    CS%uh, CS%vh, CS%u_av, CS%v_av, CS%h_av, eta_av, &
                    u(:,:,:,m), v(:,:,:,m), h(:,:,:,m), CS%eta, grid, CS)

    else ! --------------------------------------------------- not SPLIT

      if (CS%use_RK2) then
        call step_MOM_dyn_unsplit_RK2(u(:,:,:,1), v(:,:,:,1), h(:,:,:,1), &
                    Time_local, dt, fluxes, CS%p_surf_begin, CS%p_surf_end, &
                    CS%uh, CS%vh, eta_av, grid, CS)
      else
        call step_MOM_dyn_unsplit(u(:,:,:,1), v(:,:,:,1), h(:,:,:,1), &
                    Time_local, dt, fluxes, CS%p_surf_begin, CS%p_surf_end, &
                    CS%uh, CS%vh, eta_av, grid, CS)
      endif

    endif ! -------------------------------------------------- end SPLIT

    call disable_averaging(CS%diag)
    call cpu_clock_end(id_clock_dynamics)

    dtnt = dtnt + dt
    if ((MOD(n,ntstep) == 0) .or. (n==n_max)) then
      if (CS%debug) then
        call uchksum(u(:,:,:,m),"Pre-advection u",grid,haloshift=2)
        call vchksum(v(:,:,:,m),"Pre-advection v",grid,haloshift=2)
        call hchksum(h(:,:,:,m),"Pre-advection h",grid,haloshift=1)
        call uchksum(CS%uhtr,"Pre-advection uh",grid,haloshift=0)
        call vchksum(CS%vhtr,"Pre-advection vh",grid,haloshift=0)
      ! call MOM_state_chksum("Pre-advection ", u(:,:,:,m), v(:,:,:,m), &
      !                       h(:,:,:,m), CS%uhtr, CS%vhtr, grid, haloshift=1)
          if (associated(CS%tv%T)) call hchksum(CS%tv%T, "Pre-advection T",grid,haloshift=1)
          if (associated(CS%tv%S)) call hchksum(CS%tv%S, "Pre-advection S",grid,haloshift=1)
          if (associated(CS%tv%frazil)) call hchksum(CS%tv%frazil, "Pre-advection frazil",grid,haloshift=0)
          if (associated(CS%tv%salt_deficit)) call hchksum(CS%tv%salt_deficit, "Pre-advection salt deficit",grid,haloshift=0)
      ! call MOM_thermo_chksum("Pre-advection ", CS%tv, grid)
        call check_redundant("Pre-advection ", u(:,:,:,m), v(:,:,:,m), grid)
      endif

      call cpu_clock_begin(id_clock_thermo)
      call enable_averaging(dtnt,Time_local, CS%diag)

      call cpu_clock_begin(id_clock_tracer)
      call advect_tracer(h(:,:,:,m), CS%uhtr, CS%vhtr, CS%OBC, dtnt, grid, &
                         CS%tracer_CSp)
      call tracer_hordiff(h(:,:,:,m), dtnt, CS%MEKE, CS%VarMix, grid, CS%tracer_CSp, &
                          CS%tv)
      call cpu_clock_end(id_clock_tracer)

      call cpu_clock_begin(id_clock_Z_diag)
      call calculate_Z_transport(CS%uhtr, CS%vhtr, h(:,:,:,m), dtnt, grid, &
                                 CS%diag_to_Z_CSp)
      call cpu_clock_end(id_clock_Z_diag)

      if (CS%id_u_predia > 0) call post_data(CS%id_u_predia, u(:,:,:,m), CS%diag)
      if (CS%id_v_predia > 0) call post_data(CS%id_v_predia, v(:,:,:,m), CS%diag)
      if (CS%id_h_predia > 0) call post_data(CS%id_h_predia, h(:,:,:,m), CS%diag)
      if (CS%id_T_predia > 0) call post_data(CS%id_T_predia, CS%tv%T, CS%diag)
      if (CS%id_S_predia > 0) call post_data(CS%id_S_predia, CS%tv%S, CS%diag)
      if (CS%id_e_predia > 0) then
        call find_eta(h(:,:,:,m), CS%tv, grid%g_Earth, grid, eta_predia)
        call post_data(CS%id_e_predia, eta_predia, CS%diag)
      endif

      if (.not.CS%adiabatic) then
        if (CS%debug) then
          call uchksum(u(:,:,:,m),"Pre-diabatic u",grid,haloshift=2)
          call vchksum(v(:,:,:,m),"Pre-diabatic v",grid,haloshift=2)
          call hchksum(h(:,:,:,m),"Pre-diabatic h",grid,haloshift=1)
          call uchksum(CS%uhtr,"Pre-diabatic uh",grid,haloshift=0)
          call vchksum(CS%vhtr,"Pre-diabatic vh",grid,haloshift=0)
        ! call MOM_state_chksum("Pre-diabatic ", u(:,:,:,m), v(:,:,:,m), &
        !                       h(:,:,:,m), CS%uhtr, CS%vhtr, grid)
          call MOM_thermo_chksum("Pre-diabatic ", CS%tv, grid,haloshift=0)
          call check_redundant("Pre-diabatic ", u(:,:,:,m), v(:,:,:,m), grid)
        endif

        if (CS%readjust_BT_trans .and. .not.CS%BT_include_udhdt) then
          call find_total_transport(u(:,:,:,m), v(:,:,:,m), h(:,:,:,m), &
                                    CS%uhbt_in, CS%vhbt_in, dt, grid, CS)
          CS%readjust_velocity = .true.
        endif

        call cpu_clock_begin(id_clock_diabatic)
        call diabatic(u(:,:,:,m),v(:,:,:,m),h(:,:,:,m),CS%tv,fluxes,CS%visc,dtnt, &
                      grid, CS%diabatic_CSp)
        call cpu_clock_end(id_clock_diabatic)

        ! Regridding is done here, at the end of the thermodynamical time step
        ! (that may comprise several dynamical time steps)
        ! The routine 'regridding_main' can be found in 'regridding.F90'.
        if ( CS%regridding_opts%use_regridding ) then 
          call regridding_main(grid, h(:,:,:,m), CS%h_aux(:,:,:), &
                               u(:,:,:,m), v(:,:,:,m), CS%tv, CS%regridding_opts )
!         call pass_vector(u(:,:,:,m), v(:,:,:,m), grid%Domain)
!         call pass_var(CS%tv%T, grid%Domain, complete=.false.)
!         call pass_var(CS%tv%S, grid%Domain, complete=.false.)
!         call pass_var(h(:,:,:,m), grid%Domain)
!         call pass_var(h_aux(:,:,:), grid%Domain)
        end if   

        call cpu_clock_begin(id_clock_pass)
        if (grid%nonblocking_updates) then        
          pid_u = pass_vector_start(u(:,:,:,m), v(:,:,:,m), grid%Domain)
          if (CS%use_temperature) then
            pid_T = pass_var_start(CS%tv%T, grid%Domain)
            pid_S = pass_var_start(CS%tv%S, grid%Domain)
          endif
          pid_h = pass_var_start(h(:,:,:,m), grid%Domain)

          call pass_vector_complete(pid_u, u(:,:,:,m), v(:,:,:,m), grid%Domain)
          if (CS%use_temperature) then
            call pass_var_complete(pid_T, CS%tv%T, grid%Domain)
            call pass_var_complete(pid_S, CS%tv%S, grid%Domain)
          endif
          call pass_var_complete(pid_h, h(:,:,:,m), grid%Domain)
        else 
          call pass_vector(u(:,:,:,m), v(:,:,:,m), grid%Domain)
          if (CS%use_temperature) then
            call pass_var(CS%tv%T, grid%Domain, complete=.false.)
            call pass_var(CS%tv%S, grid%Domain, complete=.false.)
          endif
          call pass_var(h(:,:,:,m), grid%Domain)
        endif
        call cpu_clock_end(id_clock_pass)

        if (CS%debug) then
          call uchksum(u(:,:,:,m),"Post-diabatic u",grid,haloshift=2)
          call vchksum(v(:,:,:,m),"Post-diabatic v",grid,haloshift=2)
          call hchksum(h(:,:,:,m),"Post-diabatic h",grid,haloshift=1)
          call uchksum(CS%uhtr,"Post-diabatic uh",grid,haloshift=0)
          call vchksum(CS%vhtr,"Post-diabatic vh",grid,haloshift=0)
        ! call MOM_state_chksum("Post-diabatic ", u(:,:,:,m), v(:,:,:,m), &
        !                       h(:,:,:,m), CS%uhtr, CS%vhtr, grid, haloshift=1)
          if (associated(CS%tv%T)) call hchksum(CS%tv%T, "Post-diabatic T",grid,haloshift=1)
          if (associated(CS%tv%S)) call hchksum(CS%tv%S, "Post-diabatic S",grid,haloshift=1)
          if (associated(CS%tv%frazil)) call hchksum(CS%tv%frazil, "Post-diabatic frazil",grid,haloshift=0)
          if (associated(CS%tv%salt_deficit)) call hchksum(CS%tv%salt_deficit, "Post-diabatic salt deficit",grid,haloshift=0)
        ! call MOM_thermo_chksum("Post-diabatic ", CS%tv, grid)
          call check_redundant("Post-diabatic ", u(:,:,:,m), v(:,:,:,m), grid)
        endif
      else

        call cpu_clock_begin(id_clock_diabatic)
        call adiabatic(h(:,:,:,m), CS%tv, fluxes, dtnt, grid, CS%diabatic_CSp)
        call cpu_clock_end(id_clock_diabatic)

        if (CS%use_temperature) then
          call pass_var(CS%tv%T, grid%Domain, complete=.false.)
          call pass_var(CS%tv%S, grid%Domain, complete=.true.)
          if (CS%debug) then
            if (associated(CS%tv%T)) call hchksum(CS%tv%T, "Post-diabatic T",grid,haloshift=1)
            if (associated(CS%tv%S)) call hchksum(CS%tv%S, "Post-diabatic S",grid,haloshift=1)
          endif
        endif
      endif                                                  ! ADIABATIC
      CS%uhtr(:,:,:) = 0.0
      CS%vhtr(:,:,:) = 0.0
      call cpu_clock_end(id_clock_thermo)

      call cpu_clock_begin(id_clock_diagnostics)
      if (CS%split) then
        call calculate_diagnostic_fields(u(:,:,:,m),v(:,:,:,m),h(:,:,:,m), &
                 CS%uh, CS%vh, m, CS%tv, dtnt, grid, CS%diagnostics_CSp, &
                 CS%eta)
      else
        call calculate_diagnostic_fields(u(:,:,:,m),v(:,:,:,m),h(:,:,:,m), &
                 CS%uh, CS%vh, m, CS%tv, dtnt, grid, CS%diagnostics_CSp)
      endif
      call cpu_clock_end(id_clock_diagnostics)

      if (CS%id_T > 0) call post_data(CS%id_T, CS%tv%T, CS%diag)
      if (CS%id_S > 0) call post_data(CS%id_S, CS%tv%S, CS%diag)

      if (CS%id_Tadx > 0) call post_data(CS%id_Tadx, CS%T_adx, CS%diag)
      if (CS%id_Tady > 0) call post_data(CS%id_Tady, CS%T_ady, CS%diag)
      if (CS%id_Tdiffx > 0) call post_data(CS%id_Tdiffx, CS%T_diffx, CS%diag)
      if (CS%id_Tdiffy > 0) call post_data(CS%id_Tdiffy, CS%T_diffy, CS%diag)

      if (CS%id_Sadx > 0) call post_data(CS%id_Sadx, CS%S_adx, CS%diag)
      if (CS%id_Sady > 0) call post_data(CS%id_Sady, CS%S_ady, CS%diag)
      if (CS%id_Sdiffx > 0) call post_data(CS%id_Sdiffx, CS%S_diffx, CS%diag)
      if (CS%id_Sdiffy > 0) call post_data(CS%id_Sdiffy, CS%S_diffy, CS%diag)

      if (CS%id_Tadx_2d > 0) call post_data(CS%id_Tadx_2d, CS%T_adx_2d, CS%diag)
      if (CS%id_Tady_2d > 0) call post_data(CS%id_Tady_2d, CS%T_ady_2d, CS%diag)
      if (CS%id_Tdiffx_2d > 0) call post_data(CS%id_Tdiffx_2d, CS%T_diffx_2d, CS%diag)
      if (CS%id_Tdiffy_2d > 0) call post_data(CS%id_Tdiffy_2d, CS%T_diffy_2d, CS%diag)

      if (CS%id_Sadx_2d > 0) call post_data(CS%id_Sadx_2d, CS%S_adx_2d, CS%diag)
      if (CS%id_Sady_2d > 0) call post_data(CS%id_Sady_2d, CS%S_ady_2d, CS%diag)
      if (CS%id_Sdiffx_2d > 0) call post_data(CS%id_Sdiffx_2d, CS%S_diffx_2d, CS%diag)
      if (CS%id_Sdiffy_2d > 0) call post_data(CS%id_Sdiffy_2d, CS%S_diffy_2d, CS%diag)

      call disable_averaging(CS%diag)

      call cpu_clock_begin(id_clock_Z_diag)

      if (Time_local + set_time(int(0.5*CS%dt_therm)) > CS%Z_diag_time) then
        call enable_averaging(real(time_type_to_real(CS%Z_diag_interval)), &
                              CS%Z_diag_time, CS%diag)
        call calculate_Z_diag_fields(u(:,:,:,m),v(:,:,:,m),h(:,:,:,m), dtnt, &
                                     grid, CS%diag_to_Z_CSp)
        CS%Z_diag_time = CS%Z_diag_time + CS%Z_diag_interval
        call disable_averaging(CS%diag)
      endif
      call cpu_clock_end(id_clock_Z_diag)

      dtnt = 0.0
      CS%calc_bbl = .true.
    else  ! It is not time to do thermodynamics.
      if (.not.CS%BT_include_udhdt) CS%readjust_velocity = .false.
    endif

    call enable_averaging(dt,Time_local, CS%diag)
    if (CS%id_u > 0) call post_data(CS%id_u, u(:,:,:,m), CS%diag)
    if (CS%id_v > 0) call post_data(CS%id_v, v(:,:,:,m), CS%diag)
    if (CS%id_h > 0) call post_data(CS%id_h, h(:,:,:,m), CS%diag)

    tot_wt_ssh = tot_wt_ssh + dt
    call find_eta(h(:,:,:,1), CS%tv, grid%g_Earth, grid, ssh, eta_av)
    do j=js,je ; do i=is,ie
      CS%ave_ssh(i,j) = CS%ave_ssh(i,j) + dt*ssh(i,j)
    enddo ; enddo
    if (CS%id_ssh_inst > 0) call post_data(CS%id_ssh_inst, ssh, CS%diag)
    call disable_averaging(CS%diag)

  enddo ! End of n loop

  Itot_wt_ssh = 1.0/tot_wt_ssh
  do j=js,je ; do i=is,ie
    CS%ave_ssh(i,j) = CS%ave_ssh(i,j)*Itot_wt_ssh
  enddo ; enddo

  call enable_averaging(dt*n_max,Time_local, CS%diag)
    I_time_int = 1.0/(dt*n_max)
    if (CS%id_ssh > 0) &
      call post_data(CS%id_ssh, CS%ave_ssh, CS%diag, mask=grid%hmask)
    if (ASSOCIATED(CS%tv%frazil) .and. (CS%id_fraz > 0)) then
      allocate(frazil_ave(grid%isd:grid%ied,grid%jsd:grid%jed))
      do j=js,je ; do i=is,ie
        frazil_ave(i,j) = CS%tv%frazil(i,j) * I_time_int
      enddo ; enddo
      call post_data(CS%id_fraz, frazil_ave, CS%diag, mask=grid%hmask)
      deallocate(frazil_ave)
    endif
    if (ASSOCIATED(CS%tv%salt_deficit) .and. (CS%id_salt_deficit > 0)) then
      allocate(salt_deficit_ave(grid%isd:grid%ied,grid%jsd:grid%jed))
      do j=js,je ; do i=is,ie
        salt_deficit_ave(i,j) = CS%tv%salt_deficit(i,j) * I_time_int
      enddo ; enddo
      call post_data(CS%id_salt_deficit, salt_deficit_ave, CS%diag, mask=grid%hmask)
      deallocate(salt_deficit_ave)
    endif
    if (ASSOCIATED(CS%tv%TempxPmE) .and. (CS%id_Heat_PmE > 0)) then
      allocate(Heat_PmE_ave(grid%isd:grid%ied,grid%jsd:grid%jed))
      do j=js,je ; do i=is,ie
        Heat_PmE_ave(i,j) = CS%tv%TempxPmE(i,j) * (CS%tv%C_p * I_time_int)
      enddo ; enddo
      call post_data(CS%id_Heat_PmE, Heat_PmE_ave, CS%diag, mask=grid%hmask)      
      deallocate(Heat_PmE_ave)
    endif
    if (ASSOCIATED(CS%tv%internal_heat) .and. (CS%id_intern_heat > 0)) then
      allocate(intern_heat_ave(grid%isd:grid%ied,grid%jsd:grid%jed))
      do j=js,je ; do i=is,ie
        intern_heat_ave(i,j) = CS%tv%internal_heat(i,j) * (CS%tv%C_p * I_time_int)
      enddo ; enddo
      call post_data(CS%id_intern_heat, intern_heat_ave, CS%diag, mask=grid%hmask)
      deallocate(intern_heat_ave)
    endif
  call disable_averaging(CS%diag)

  call calculate_surface_state(state, u(:,:,:,m), v(:,:,:,m), h(:,:,:,m), &
                               CS%ave_ssh, grid, CS, fluxes%p_surf_full)

  call enable_averaging(dt*n_max,Time_local, CS%diag)
    if (CS%id_sst > 0) &
      call post_data(CS%id_sst, state%SST, CS%diag, mask=grid%hmask)
    if (CS%id_sst_sq > 0) then
      do j=js,je ; do i=is,ie
        CS%SST_sq(i,j) = state%SST(i,j)*state%SST(i,j)
      enddo ; enddo
      call post_data(CS%id_sst_sq, CS%SST_sq, CS%diag, mask=grid%hmask)
    endif
    if (CS%id_sss > 0) &
      call post_data(CS%id_sss, state%SSS, CS%diag, mask=grid%hmask)
    if (CS%id_ssu > 0) &
      call post_data(CS%id_ssu, state%u, CS%diag, mask=grid%umask)
    if (CS%id_ssv > 0) &
      call post_data(CS%id_ssv, state%v, CS%diag, mask=grid%umask)
    if (CS%id_speed > 0) then
      allocate(sfc_speed(grid%isd:grid%ied,grid%jsd:grid%jed))
      do j=js,je ; do i=is,ie
        sfc_speed(i,j) = sqrt(0.5*(state%u(I-1,j)**2 + state%u(I,j)**2) + &
                              0.5*(state%v(i,J-1)**2 + state%v(i,J)**2))
      enddo ; enddo
      call post_data(CS%id_speed, sfc_speed, CS%diag, mask=grid%hmask)
      deallocate(sfc_speed)
    endif
  call disable_averaging(CS%diag)

  if (CS%interp_p_surf) then ; do j=jsd,jed ; do i=isd,ied
    CS%p_surf_prev(i,j) = fluxes%p_surf(i,j)
  enddo ; enddo ; endif

  call cpu_clock_end(id_clock_ocean)
  step_MOM = m

end function step_MOM

! ============================================================================

subroutine find_total_transport(u_in, v_in, h_in, uh_tot, vh_tot, dt, G, CS)
  real, dimension(NXMEMQ_,NYMEM_,NKMEM_), intent(in)    :: u_in
  real, dimension(NXMEM_,NYMEMQ_,NKMEM_), intent(in)    :: v_in
  real,                                   intent(in)    :: dt
  real, dimension(NXMEMQ_,NYMEM_),        intent(out)   :: uh_tot
  real, dimension(NXMEM_,NYMEMQ_),        intent(out)   :: vh_tot
  real, dimension(NXMEM_,NYMEM_,NKMEM_),  intent(in)    :: h_in
  type(ocean_grid_type),                  intent(inout) :: G
  type(MOM_control_struct),               pointer       :: CS
!   This subroutine determines the vertically summed transport based on input
! velocities and thicknesses.  The individual layers' transports are not retained.

! Arguments: u_in - The input zonal velocity, in m s-1. (Intent in.)
!  (in)      v_in - The input meridional velocity, in m s-1.
!  (in)      h_in - The input layer thicknesses, in m or kg m-2, depending on
!                   whether the Boussinesq approximation is made.
!  (out)     uh_tot - The vertically summed zonal and meridional volume or mass
!  (out)     vh_tot -  transports, in m3 s-1 or kg s-1.
!  (in)      dt - The time step in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure set up by initialize_MOM.

  ! Temporary arrays to contain layer thickness fluxes in m3 s-1 or kg s-1.
  real, dimension(SZIQ_(G),SZJ_(G),SZK_(G)) :: uh_temp 
  real, dimension(SZI_(G),SZJQ_(G),SZK_(G)) :: vh_temp
  ! A temporary array to contain layer projected thicknesses in m or kg m-2.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G))  :: h_temp
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call cpu_clock_begin(id_clock_continuity)
  call continuity(u_in, v_in, h_in, h_temp, uh_temp, vh_temp, dt, G, &
                  CS%continuity_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_continuity)

  do j=js,je ; do I=is-1,ie ; uh_tot(I,j) = uh_temp(I,j,1) ; enddo ; enddo
  do k=2,nz ; do j=js,je ; do I=is-1,ie
    uh_tot(I,j) = uh_tot(I,j) + uh_temp(I,j,k)
  enddo ; enddo ; enddo
  do J=js-1,je ; do i=is,ie ; vh_tot(i,J) = vh_temp(i,J,1) ; enddo ; enddo
  do k=2,nz ; do J=js-1,je ; do i=is,ie
    vh_tot(i,J) = vh_tot(i,J) + vh_temp(i,J,k)
  enddo ; enddo ; enddo

end subroutine find_total_transport

! ============================================================================

subroutine initialize_MOM(Time, param_file, dirs, CS, Time_in)
  type(time_type), target, intent(inout) :: Time
  type(param_file_type), intent(out)     :: param_file
  type(directories), intent(out)         :: dirs
  type(MOM_control_struct), pointer      :: CS
  type(time_type), optional, intent(in)  :: Time_in
! Arguments: Time - The model time, set in this routine.
!  (out)     param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      dirs - A structure containing several relevant directory paths.
!  (out)     CS - A pointer set in this routine to the MOM control structure.
!  (in)      Time_in - An optional time passed to MOM_initialize to use when
!                      the model is not being started from a restart file.
  type(ocean_grid_type), pointer :: grid ! A pointer to a structure containing
                                  ! metrics and related information.
  type(diag_ptrs), pointer :: diag
  character(len=4), parameter :: vers_num = 'v2.0'
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: Isdq, Iedq, Jsdq, Jedq
  real    :: dtbt
  real    :: Z_diag_int      ! The minimum interval between calculations of
                             ! Depth-space diagnostic quantities, in s.
  real, pointer, dimension(:,:,:,:) :: &
    u, &                     ! u : Zonal velocity, in m s-1.
    v, &                     ! v : Meridional velocity, in m s-1.
    h                        ! h : Layer thickness, in m or kg m-2.
  real, allocatable, dimension(:,:,:) :: e ! The interface heights in m.
  type(MOM_restart_CS),  pointer :: restart_CSp_tmp => NULL()
  integer :: nkml, nkbl, verbosity
  real    :: default_val     ! A default value for a parameter.
  logical :: new_sim
  logical :: use_geothermal  ! If true, apply geothermal heating.
  logical :: use_EOS         ! If true, density is calculated from T & S using
                             ! an equation of state.
  logical :: use_tides       ! If true, tidal momentum forcing is used.
  logical :: save_IC         ! If true, save the initial conditions.
  character(len=80) :: IC_file ! A file into which the initial conditions are
                             ! written in a new run if save_IC is true.
  type(vardesc) :: vd
  type(time_type) :: Start_time
  type(MOM_initialization_struct) :: init_CS
  type(ocean_internal_state) :: MOM_internal_state

  if (associated(CS)) then
    call MOM_error(WARNING, "initialize_MOM called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)
  grid => CS%grid ; CS%Time => Time
  diag => CS%diag

  id_clock_init = cpu_clock_id('Ocean Initialization', grain=CLOCK_SUBCOMPONENT)
  call cpu_clock_begin(id_clock_init)

  Start_time = Time ; if (present(Time_in)) Start_time = Time_in

  call Get_MOM_Input(param_file, dirs)

  verbosity = 2 ; call read_param(param_file, "VERBOSITY", verbosity)
  call MOM_set_verbosity(verbosity)

  call find_obsolete_params(param_file)
  call MOM_domains_init(grid%domain, param_file)
  call MOM_checksums_init(param_file)

  call MOM_io_init(param_file)
  call diag_mediator_init(param_file)
  call MOM_grid_init(grid, param_file)
  is = grid%isc ; ie = grid%iec ; js = grid%jsc ; je = grid%jec ; nz = grid%ke
  isd = grid%isd ; ied = grid%ied ; jsd = grid%jsd ; jed = grid%jed
  Isdq = grid%Isdq ; Iedq = grid%Iedq ; Jsdq = grid%Jsdq ; Jedq = grid%Jedq
  call set_diag_mediator_grid(grid, diag)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, "MOM", version, tagname, "")
  call get_param(param_file, "MOM", "VERBOSITY", verbosity,  &
                 "Integer controlling level of messaging\n" // &
                 "\t0 = Only FATAL messages\n" // &
                 "\t2 = Only FATAL, WARNING, NOTE [default]\n" // &
                 "\t9 = All)", default=2)
  call get_param(param_file, "MOM", "SPLIT", CS%split, &
                 "Use the split time stepping if true.", default=.true.)
  if (.not. CS%split) &
    call get_param(param_file, "MOM", "USE_RK2", CS%use_RK2, &
                 "If true, use RK2 instead of RK3 in the unsplit time stepping.", &
                 default=.false.)
  call get_param(param_file, "MOM", "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state \n"//&
                 "variables.", default=.true.)
  call get_param(param_file, "MOM", "USE_EOS", use_EOS, &
                 "If true,  density is calculated from temperature and \n"//&
                 "salinity with an equation of state.  If USE_EOS is \n"//&
                 "true, ENABLE_THERMODYNAMICS must be true as well.", &
                 default=CS%use_temperature)
  call get_param(param_file, "MOM", "ADIABATIC", CS%adiabatic, &
                 "There are no diapycnal mass fluxes if ADIABATIC is \n"//&
                 "true. This assumes that KD = KDML = 0.0 and that \n"//&
                 "there is no buoyancy forcing, but makes the model \n"//&
                 "faster by eliminating subroutine calls.", default=.false.)
  call get_param(param_file, "MOM", "BULKMIXEDLAYER", CS%bulkmixedlayer, &
                 "If true, use a Kraus-Turner-like bulk mixed layer \n"//&
                 "with transitional buffer layers.  Layers 1 through  \n"//&
                 "NKML+NKBL have variable densities. There must be at \n"//&
                 "least NKML+NKBL+1 layers if BULKMIXEDLAYER is true. \n"//&
                 "The default is the same setting as ENABLE_THERMODYNAMICS.", &
                 default=CS%use_temperature)
  call get_param(param_file, "MOM", "THICKNESSDIFFUSE", CS%thickness_diffuse, &
                 "If true, interfaces or isopycnal surfaces are diffused, \n"//&
                 "depending on the value of FULL_THICKNESSDIFFUSE.", &
                 default=.false.)
  call get_param(param_file, "MOM", "THICKNESSDIFFUSE_FIRST", &
                                     CS%thickness_diffuse_first, &
                 "If true, do thickness diffusion before dynamics.\n"//&
                 "This is only used if THICKNESSDIFFUSE is true.", &
                 default=.false.)
  call get_param(param_file, "MOM", "MIXEDLAYER_RESTRAT",CS%mixedlayer_restrat, &
                 "If true, a density-gradient dependent re-stratifying \n"//&
                 "flow is imposed in the mixed layer. \n"//&
                 "This is only used if BULKMIXEDLAYER is true.", default=.false.)

  call get_param(param_file, "MOM", "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", default=.false.)
  call get_param(param_file, "MOM", "DEBUG_TRUNCATIONS", CS%debug_truncations, &
                 "If true, calculate all diagnostics that are useful for \n"//&
                 "debugging truncations.", default=.false.)

  call get_param(param_file, "MOM", "DT", CS%dt, &
                 "The (baroclinic) dynamics time step.  The time-step that \n"//&
                 "is actually used will be an integer fraction of the \n"//&
                 "forcing time-step (DT_FORCING in ocean-only mode or the \n"//&
                 "coupling timestep in coupled mode.)", units="s", &
                 fail_if_missing=.true.)
  call get_param(param_file, "MOM", "DT_THERM", CS%dt_therm, &
                 "The thermodynamic and tracer advection time step. \n"//&
                 "Ideally DT_THERM should be an integer multiple of DT \n"//&
                 "and less than the forcing or coupling time-step. \n"//&
                 "By default DT_THERM is set to DT.", units="s", default=CS%dt)

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
  call get_param(param_file, "MOM", "HMIX", CS%Hmix, &
                 "If BULKMIXEDLAYER is false, HMIX is the depth over \n"//&
                 "which to average to find surface properties like SST \n"//&
                 "and SSS, and over which the vertical viscosity and \n"//&
                 "diapycnal diffusivity are elevated.  HMIX is only used \n"//&
                 "directly if BULKMIXEDLAYER is false, but provides a \n"//&
                 "default value for other variables if BULKMIXEDLAYER is \n"//&
                 "true.", units="m", default=1.0)
  call get_param(param_file, "MOM", "MIN_Z_DIAG_INTERVAL", Z_diag_int, &
                 "The minimum amount of time in seconds between \n"//&
                 "calculations of depth-space diagnostics. Making this \n"//&
                 "larger than DT_THERM reduces the  performance penalty \n"//&
                 "of regridding to depth online.", units="s", default=0.0)
  call get_param(param_file, "MOM", "FLUX_BT_COUPLING", CS%flux_BT_coupling, &
                 "If true, use mass fluxes to ensure consistency between \n"//&
                 "the baroclinic and barotropic modes. This is only used \n"//&
                 "if SPLIT is true.", default=.true.)
  call get_param(param_file, "MOM", "INTERPOLATE_P_SURF", CS%interp_p_surf, &
                 "If true, linearly interpolate the surface pressure \n"//&
                 "over the coupling time step, using the specified value \n"//&
                 "at the end of the step.", default=.false.)
  call get_param(param_file, "MOM", "SSH_SMOOTHING_PASSES", CS%smooth_ssh_passes, &
                 "The number of Laplacian smoothing passes to apply to the \n"//&
                 "the sea surface height that is reported to the sea-ice.", &
                 units="nondim", default=0.0)

  if (CS%split) then
    call get_param(param_file, "MOM", "DTBT", dtbt, default=-0.98)
    default_val = CS%dt_therm ; if (dtbt > 0.0) default_val = -1.0
    CS%dtbt_reset_period = -1.0
    call get_param(param_file, "MOM", "DTBT_RESET_PERIOD", CS%dtbt_reset_period, &
                 "The period between recalculations of DTBT (if DTBT <= 0). \n"//&
                 "If DTBT_RESET_PERIOD is negative, DTBT is set based \n"//&
                 "only on information available at initialization.  If \n"//&
                 "dynamic, DTBT will be set at least every forcing time \n"//&
                 "step, and if 0, every dynamics time step.  The default is \n"//&
                 "set by DT_THERM.  This is only used if SPLIT is true.", &
                 units="s", default=default_val, do_not_read=(dtbt > 0.0))

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
                 "in the barotropic continuity equation.", default=.false.)
  endif

  if (CS%split .and. CS%flux_BT_coupling) then
    call get_param(param_file, "MOM", "BT_INCLUDE_UDHDT", CS%BT_include_udhdt, &
                 "If true, include the barotropic transport tendancies \n"//&
                 "from sum(u dhdt) and sum(v dhdt) in the barotropic \n"//&
                 "solver.", default=.false.)
  endif
  if (.not.(CS%split .and. CS%flux_BT_coupling) .or. CS%adiabatic) &
    CS%readjust_BT_trans = .false.

  ! This is here in case these values are used inappropriately.
  CS%use_frazil = .false. ; CS%bound_salinity = .false. ; CS%tv%P_Ref = 2.0e7
  if (CS%use_temperature) then
    call get_param(param_file, "MOM", "FRAZIL", CS%use_frazil, &
                 "If true, water freezes if it gets too cold, and the \n"//&
                 "the accumulated heat deficit is returned in the \n"//&
                 "surface state.  FRAZIL is only used if \n"//&
                 "ENABLE_THERMODYNAMICS is true.", default=.false.)
    call get_param(param_file, "MOM", "DO_GEOTHERMAL", use_geothermal, &
                 "If true, apply geothermal heating.", default=.false.)
    call get_param(param_file, "MOM", "BOUND_SALINITY", CS%bound_salinity, &
                 "If true, limit salinity to being positive. (The sea-ice \n"//&
                 "model may ask for more salt than is available and \n"//&
                 "drive the salinity negative otherwise.)", default=.false.)
    call get_param(param_file, "MOM", "C_P", CS%tv%C_p, &
                 "The heat capacity of sea water, approximated as a \n"//&
                 "constant. This is only used if ENABLE_THERMODYNAMICS is \n"//&
                 "true. The default value is from the TEOS-10 definition \n"//&
                 "of conservative temperature.", units="J kg-1 K-1", &
                 default=3991.86795711963)
  endif
  if (use_EOS) call get_param(param_file, "MOM", "P_REF", CS%tv%P_Ref, &
                 "The pressure that is used for calculating the coordinate \n"//&
                 "density.  (1 Pa = 1e4 dbar, so 2e7 is commonly used.) \n"//&
                 "This is only used if USE_EOS and ENABLE_THERMODYNAMICS \n"//&
                 "are true.", units="Pa", default=2.0e7)

  call get_param(param_file, "MOM", "TIDES", use_tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  if (CS%bulkmixedlayer) then
    call get_param(param_file, "MOM", "NKML", nkml, &
                 "The number of sublayers within the mixed layer if \n"//&
                 "BULKMIXEDLAYER is true.", units="nondim", default=2)
    call get_param(param_file, "MOM", "NKBL", nkbl, &
                 "The number of layers that are used as variable density \n"//&
                 "buffer layers if BULKMIXEDLAYER is true.", units="nondim", &
                 default=2)
  endif

  call get_param(param_file, "MOM", "CHECK_BAD_SURFACE_VALS", &
                                     CS%check_bad_surface_vals, &
                 "If true, check the surface state for ridiculous values.", &
                 default=.false.)
  if (CS%check_bad_surface_vals) then
    call get_param(param_file, "MOM", "BAD_VAL_SSH_MAX", CS%bad_val_ssh_max, &
                 "The value of SSH above which a bad value message is \n"//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", units="m", &
                 default=20.0)
    call get_param(param_file, "MOM", "BAD_VAL_SSS_MAX", CS%bad_val_sss_max, &
                 "The value of SSS above which a bad value message is \n"//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", units="PSU", &
                 default=45.0)
    call get_param(param_file, "MOM", "BAD_VAL_SST_MAX", CS%bad_val_sst_max, &
                 "The value of SST above which a bad value message is \n"//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", &
                 units="deg C", default=45.0)
    call get_param(param_file, "MOM", "BAD_VAL_SST_MIN", CS%bad_val_sst_min, &
                 "The value of SST below which a bad value message is \n"//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", &
                 units="deg C", default=-2.1)
  endif

  call get_param(param_file, "MOM", "SAVE_INITIAL_CONDS", save_IC, &
                 "If true, write the initial conditions to a file given \n"//&
                 "by IC_OUTPUT_FILE.", default=.false.)
  call get_param(param_file, "MOM", "IC_OUTPUT_FILE", IC_file, &
                 "The file into which to write the initial conditions.", &
                 default="MOM_IC")

  ! Check for inconsistent settings.
  if (CS%adiabatic .and. CS%use_temperature) call MOM_error(WARNING, &
    "MOM: ADIABATIC and ENABLE_THERMODYNAMICS both defined is usually unwise.")
  if (use_EOS .and. .not.CS%use_temperature) call MOM_error(FATAL, &
    "MOM: ENABLE_THERMODYNAMICS must be defined to use USE_EOS.")
  if (CS%adiabatic .and. CS%bulkmixedlayer) call MOM_error(FATAL, &
    "MOM: ADIABATIC and BULKMIXEDLAYER can not both be defined.")
  if (CS%mixedlayer_restrat .and. .not.CS%bulkmixedlayer) call MOM_error(FATAL, &
    "MOM: MIXEDLAYER_RESTRAT true requires BULKMIXEDLAYER to be true to work.")

  ! Allocate the auxiliary non-symmetric domain for debugging or I/O purposes.
  if (CS%debug .or. grid%symmetric) &
    call MOM_domains_init(grid%Domain_aux, param_file, symmetric=.false.)

  call MOM_timing_init(CS)

  call advect_tracer_init(param_file, CS%tracer_CSp)

! Allocate and initialize space for the primary MOM variables.
  ALLOC(CS%u(Isdq:Iedq,jsd:jed,nz,2)) ; CS%u(:,:,:,:) = 0.0
  ALLOC(CS%v(isd:ied,Jsdq:Jedq,nz,2)) ; CS%v(:,:,:,:) = 0.0
  ALLOC(CS%h(isd:ied,jsd:jed,nz,2))   ; CS%h(:,:,:,:) = grid%Angstrom
  u => CS%u ; v => CS%v ; h => CS%h
  ALLOC(CS%uh(Isdq:Iedq,jsd:jed,nz))  ; CS%uh(:,:,:) = 0.0
  ALLOC(CS%vh(isd:ied,Jsdq:Jedq,nz))  ; CS%vh(:,:,:) = 0.0
  ALLOC(CS%diffu(Isdq:Iedq,jsd:jed,nz)) ; CS%diffu(:,:,:) = 0.0
  ALLOC(CS%diffv(isd:ied,Jsdq:Jedq,nz)) ; CS%diffv(:,:,:) = 0.0
  ALLOC(CS%CAu(Isdq:Iedq,jsd:jed,nz)) ; CS%CAu(:,:,:) = 0.0
  ALLOC(CS%CAv(isd:ied,Jsdq:Jedq,nz)) ; CS%CAv(:,:,:) = 0.0
  ALLOC(CS%PFu(Isdq:Iedq,jsd:jed,nz)) ; CS%PFu(:,:,:) = 0.0
  ALLOC(CS%PFv(isd:ied,Jsdq:Jedq,nz)) ; CS%PFv(:,:,:) = 0.0
  if (CS%use_temperature) then
    ALLOC(CS%T(isd:ied,jsd:jed,nz))   ; CS%T(:,:,:) = 0.0
    ALLOC(CS%S(isd:ied,jsd:jed,nz))   ; CS%S(:,:,:) = 0.0
    CS%tv%T => CS%T ; CS%tv%S => CS%S
    call register_tracer(CS%tv%T, "T", param_file, CS%tracer_CSp)
    call register_tracer(CS%tv%S, "S", param_file, CS%tracer_CSp)
  endif
  if (CS%use_frazil) then
    allocate(CS%tv%frazil(isd:ied,jsd:jed)) ; CS%tv%frazil(:,:) = 0.0
  endif
  if (CS%bound_salinity) then
    allocate(CS%tv%salt_deficit(isd:ied,jsd:jed)) ; CS%tv%salt_deficit(:,:)=0.0
  endif
  
  if (CS%bulkmixedlayer) then
    if (.not.use_EOS) call MOM_error(FATAL, &
      "initialize_MOM: A bulk mixed layer can only be used with T & S as "//&
      "state variables. Add #define USE_EOS.")
    grid%nkml = nkml
    grid%nk_rho_varies = nkml + nkbl
    allocate(CS%tv%Hml(isd:ied,jsd:jed)) ; CS%tv%Hml(:,:) = 0.0
  else
    grid%nkml = 0 ; grid%nk_rho_varies = 0
  endif

  ALLOC(CS%uhtr(Isdq:Iedq,jsd:jed,nz)) ; CS%uhtr(:,:,:) = 0.0
  ALLOC(CS%vhtr(isd:ied,Jsdq:Jedq,nz)) ; CS%vhtr(:,:,:) = 0.0

  if (CS%debug_truncations) then
    allocate(CS%u_prev(Isdq:Iedq,jsd:jed,nz)) ; CS%u_prev(:,:,:) = 0.0
    allocate(CS%v_prev(isd:ied,Jsdq:Jedq,nz)) ; CS%u_prev(:,:,:) = 0.0
  endif

  MOM_internal_state%u => u ; MOM_internal_state%v => v ; MOM_internal_state%h =>h
  MOM_internal_state%uh => CS%uh ; MOM_internal_state%vh => CS%vh
  MOM_internal_state%diffu => CS%diffu ; MOM_internal_state%diffv => CS%diffv
  MOM_internal_state%PFu => CS%PFu ; MOM_internal_state%PFv => CS%PFv
  MOM_internal_state%CAu => CS%CAu ; MOM_internal_state%CAv => CS%CAv
  if (CS%use_temperature) then
    MOM_internal_state%T => CS%T ; MOM_internal_state%S => CS%S
  endif
  if (CS%split) then
    ALLOC(CS%eta(isd:ied,jsd:jed))      ; CS%eta(:,:) = 0.0
    ALLOC(CS%uhbt(Isdq:Iedq,jsd:jed))   ; CS%uhbt(:,:) = 0.0
    ALLOC(CS%vhbt(isd:ied,Jsdq:Jedq))   ; CS%vhbt(:,:) = 0.0
    ALLOC(CS%uhbt_in(Isdq:Iedq,jsd:jed)) ; CS%uhbt_in(:,:) = 0.0
    ALLOC(CS%vhbt_in(isd:ied,Jsdq:Jedq)) ; CS%vhbt_in(:,:) = 0.0
    ALLOC(CS%u_av(Isdq:Iedq,jsd:jed,nz)); CS%u_av(:,:,:) = 0.0
    ALLOC(CS%v_av(isd:ied,Jsdq:Jedq,nz)); CS%v_av(:,:,:) = 0.0
    ALLOC(CS%h_av(isd:ied,jsd:jed,nz))  ; CS%h_av(:,:,:) = grid%Angstrom
    ALLOC(CS%eta_PF(isd:ied,jsd:jed))   ; CS%eta_PF(:,:) = 0.0
    ALLOC(CS%pbce(isd:ied,jsd:jed,nz))  ; CS%pbce(:,:,:) = 0.0
    ALLOC(CS%u_accel_bt(Isdq:Iedq,jsd:jed,nz)) ; CS%u_accel_bt(:,:,:) = 0.0
    ALLOC(CS%v_accel_bt(isd:ied,Jsdq:Jedq,nz)) ; CS%v_accel_bt(:,:,:) = 0.0
    ALLOC(CS%visc_rem_u(Isdq:Iedq,jsd:jed,nz)) ; CS%visc_rem_u(:,:,:) = 0.0
    ALLOC(CS%visc_rem_v(isd:ied,Jsdq:Jedq,nz)) ; CS%visc_rem_v(:,:,:) = 0.0

    MOM_internal_state%u_accel_bt => CS%u_accel_bt
    MOM_internal_state%v_accel_bt => CS%v_accel_bt
    MOM_internal_state%pbce => CS%pbce
    MOM_internal_state%u_av => CS%u_av
    MOM_internal_state%v_av => CS%v_av
  endif
  if (CS%interp_p_surf) then
    allocate(CS%p_surf_prev(isd:ied,jsd:jed)) ; CS%p_surf_prev(:,:) = 0.0
  endif
  allocate(CS%taux_bot(Isdq:Iedq,jsd:jed)) ; CS%taux_bot(:,:) = 0.0
  allocate(CS%tauy_bot(isd:ied,Jsdq:Jedq)) ; CS%tauy_bot(:,:) = 0.0

  ALLOC(CS%ave_ssh(isd:ied,jsd:jed)) ; CS%ave_ssh(:,:) = 0.0

! Use the Wright equation of state by default, unless otherwise specified
! Note: this line and the following block ought to be in a separate
! initialization routine for tv.
  if (use_EOS) call select_eqn_of_state(param_file,CS%tv%eqn_of_state)
  if (CS%use_temperature) then
    allocate(CS%tv%TempxPmE(isd:ied,jsd:jed))
    CS%tv%TempxPmE(:,:) = 0.0
    if (use_geothermal) then
      allocate(CS%tv%internal_heat(isd:ied,jsd:jed))
      CS%tv%internal_heat(:,:) = 0.0
    endif
  endif

!   Set the fields that are needed for bitwise identical restarting
! the time stepping scheme.
  call restart_init(param_file, CS%restart_CSp)
  call set_restart_fields(grid, param_file, CS)
  if (CS%split) then
    call register_restarts_dyn_split_RK2(grid, param_file, CS, CS%restart_CSp)
  else
    if (CS%use_RK2) then
      call register_restarts_dyn_unsplit_RK2(grid, param_file, CS, CS%restart_CSp)
    else
      call register_restarts_dyn_unsplit(grid, param_file, CS, CS%restart_CSp)
    endif
  endif
!   This subroutine calls user-specified tracer registration routines.
! Additional calls can be added to MOM_tracer_flow_control.F90.
  call call_tracer_register(grid, param_file, CS%tracer_flow_CSp, &
                            diag, CS%tracer_CSp, CS%restart_CSp)
  call MEKE_alloc_register_restart(grid, param_file, CS%MEKE, CS%restart_CSp)

!   Initialize all of the relevant fields.
  if (associated(CS%tracer_CSp)) &
    init_CS%advect_tracer_CSp => CS%tracer_CSp

  call cpu_clock_begin(id_clock_MOM_init)
  call MOM_initialize(u(:,:,:,1), v(:,:,:,1), h(:,:,:,1), CS%tv, Time, &
                       grid, param_file, dirs, CS%restart_CSp, init_CS, Time_in)
  call cpu_clock_end(id_clock_MOM_init)

  if (associated(init_CS%advect_tracer_CSp)) &
    CS%tracer_CSp => init_CS%advect_tracer_CSp
  if (associated(init_CS%OBC)) then
    CS%OBC => init_CS%OBC
    call open_boundary_init(Time, grid, param_file, diag, CS%open_boundary_CSp)
  endif
! if (associated(init_CS%sponge_CSp)) CS%sponge_CSp => init_CS%sponge_CSp

  if (use_tides) call tidal_forcing_init(Time, grid, param_file, CS%tides_CSp)

  call continuity_init(Time, grid, param_file, diag, CS%continuity_CSp)
  call CoriolisAdv_init(Time, grid, param_file, diag, CS%CoriolisAdv_CSp)
  call PressureForce_init(Time, grid, param_file, diag, CS%PressureForce_CSp, &
                          CS%tides_CSp)

  call hor_visc_init(Time, grid, param_file, diag, CS%hor_visc_CSp)
  call vertvisc_init(MOM_internal_state, Time, grid, param_file, diag, dirs, &
                     CS%ntrunc, CS%vertvisc_CSp)
  call set_visc_init(Time, grid, param_file, diag, CS%visc, CS%set_visc_CSp)
  call thickness_diffuse_init(Time, grid, param_file, diag, CS%thickness_diffuse_CSp)
  if (CS%mixedlayer_restrat) &
    call mixedlayer_restrat_init(Time, grid, param_file, diag, CS%mixedlayer_restrat_CSp)
  call MEKE_init(Time, grid, param_file, diag, CS%MEKE_CSp, CS%MEKE)
  call VarMix_init(Time, grid, param_file, diag, CS%VarMix)
  ! Need an ALE CS !!!! -AJA
  ALLOC(CS%h_aux(isd:ied,jsd:jed,nz)); CS%h_aux(:,:,:) = 0.
  call initialize_regridding(param_file, CS%regridding_opts, grid, &
                             h(:,:,:,:), CS%h_aux(:,:,:), u(:,:,:,1), v(:,:,:,1), CS%tv)
  call MOM_diagnostics_init(MOM_internal_state, Time, grid, param_file, &
                             diag, CS%diagnostics_CSp)

  CS%Z_diag_interval = set_time(int((CS%dt_therm) * &
       max(1,floor(0.01 + Z_diag_int/(CS%dt_therm)))))
  call MOM_diag_to_Z_init(Time, grid, param_file, diag, CS%diag_to_Z_CSp)
  CS%Z_diag_time = Start_time + CS%Z_diag_interval * (1 + &
    ((Time + set_time(int(CS%dt_therm))) - Start_time) / CS%Z_diag_interval)

  if (associated(init_CS%sponge_CSp)) &
    call init_sponge_diags(Time, grid, diag, init_CS%sponge_CSp)
  if (CS%adiabatic) then
    call adiabatic_driver_init(Time, grid, param_file, diag, CS%diabatic_CSp, &
                               CS%tracer_flow_CSp, CS%diag_to_Z_CSp)
  else
    call diabatic_driver_init(Time, grid, param_file, diag, CS%diabatic_CSp, &
                     CS%tracer_flow_CSp, init_CS%sponge_CSp, CS%diag_to_Z_CSp)
  endif

  call register_diags(Time, grid, CS)

  call advect_tracer_diag_init(Time, grid, diag, CS%tracer_CSp)
  if (CS%use_temperature) then
    ! If needed T_adx, etc., would have been allocated in register_diags.
    call add_tracer_diagnostics("T", CS%tracer_CSp, CS%T_adx, CS%T_ady, &
                                CS%T_diffx, CS%T_diffy)
    call add_tracer_diagnostics("S", CS%tracer_CSp, CS%S_adx, CS%S_ady, &
                                CS%S_diffx, CS%S_diffy)
    call add_tracer_2d_diagnostics("T", CS%tracer_CSp, CS%T_adx_2d, CS%T_ady_2d, &
                                   CS%T_diffx_2d, CS%T_diffy_2d)
    call add_tracer_2d_diagnostics("S", CS%tracer_CSp, CS%S_adx_2d, CS%S_ady_2d, &
                                    CS%S_diffx_2d, CS%S_diffy_2d)
    call register_Z_tracer(CS%tv%T, "temp_z", "Potential Temperature", "degC", Time, &
                           grid, CS%diag_to_Z_CSp)
    call register_Z_tracer(CS%tv%S, "salt_z", "Salinity", "PSU", Time, &
                           grid, CS%diag_to_Z_CSp)
  endif

  ! This subroutine initializes any tracer packages.
  new_sim = ((dirs%input_filename(1:1) == 'n') .and. &
             (LEN_TRIM(dirs%input_filename) == 1))
  call tracer_flow_control_init(.not.new_sim, Time, grid, h(:,:,:,1), CS%OBC, &
           CS%tracer_flow_CSp, init_CS%sponge_CSp, CS%diag_to_Z_CSp)

  call cpu_clock_begin(id_clock_pass_init)
  call pass_vector(u(:,:,:,1),v(:,:,:,1), grid%Domain)
  call pass_var(h(:,:,:,1), grid%Domain)

  if (CS%use_temperature) then
    call pass_var(CS%tv%T, grid%Domain)
    call pass_var(CS%tv%S, grid%Domain)
  endif
  call cpu_clock_end(id_clock_pass_init)

  if (CS%split) then
    call initialize_dyn_split_RK2(u(:,:,:,1), v(:,:,:,1), h(:,:,:,1), Time, &
                                  grid, param_file, diag, CS, CS%restart_CSp)
  else
    if (CS%use_RK2) then
      call initialize_dyn_unsplit_RK2(u(:,:,:,1), v(:,:,:,1), h(:,:,:,1), Time, &
                                  grid, param_file, diag, CS, CS%restart_CSp)
    else
      call initialize_dyn_unsplit(u(:,:,:,1), v(:,:,:,1), h(:,:,:,1), Time, &
                                  grid, param_file, diag, CS, CS%restart_CSp)
    endif
  endif

  call write_static_fields(grid, CS%diag)
  call enable_averaging(0.0,Time, CS%diag)

!  call calculate_diagnostic_fields(u(:,:,:,1),v(:,:,:,1),h(:,:,:,1), &
!            uh, vh, 1, CS%tv, 0.0, grid, CS%diagnostics_CSp)

!  if (CS%id_u > 0) call post_data(CS%id_u, CS%u(:,:,:,1), CS%diag)
!  if (CS%id_v > 0) call post_data(CS%id_v, CS%v(:,:,:,1), CS%diag)
!  if (CS%id_h > 0) call post_data(CS%id_h, CS%h(:,:,:,1), CS%diag)
!  if (CS%id_T > 0) call post_data(CS%id_T, CS%tv%T, CS%diag)
!  if (CS%id_S > 0) call post_data(CS%id_S, CS%tv%S, CS%diag)

  call disable_averaging(CS%diag)

  if (CS%use_frazil) then
    if (.not.query_initialized(CS%tv%frazil,"frazil",CS%restart_CSp)) &
      CS%tv%frazil(:,:) = 0.0
  endif

  if (CS%interp_p_surf) then 
    CS%p_surf_prev_set = &
      query_initialized(CS%p_surf_prev,"p_surf_prev",CS%restart_CSp)

    if (CS%p_surf_prev_set) call pass_var(CS%p_surf_prev,grid%domain)
  endif
  
  if (.not.query_initialized(CS%ave_ssh,"ave_ssh",CS%restart_CSp)) then
    if (CS%split) then
      call find_eta(h(:,:,:,1), CS%tv, grid%g_Earth, grid, CS%ave_ssh, CS%eta(:,:))
    else
      call find_eta(h(:,:,:,1), CS%tv, grid%g_Earth, grid, CS%ave_ssh)
    endif
  endif

  if (save_IC .and. .not.((dirs%input_filename(1:1) == 'r') .and. &
                          (LEN_TRIM(dirs%input_filename) == 1))) then
    allocate(restart_CSp_tmp)
    restart_CSp_tmp = CS%restart_CSp
    allocate(e(SZI_(grid),SZJ_(grid),SZK_(grid)+1))
    call find_eta(h(:,:,:,1), CS%tv, grid%g_Earth, grid, e)
    vd = vardesc("eta","Interface heights",'h','i','1',"meter", 'd')
    call register_restart_field(e, e, vd, .true., restart_CSp_tmp)
    
    call save_restart(dirs%output_directory, Time, 1, grid, &
                      restart_CSp_tmp, filename=IC_file)
    deallocate(e)
    deallocate(restart_CSp_tmp)
  endif

!  call calculate_surface_state(state, u(:,:,:,1), v(:,:,:,1), h(:,:,:,1), &
!                               CS%ave_ssh, grid, CS)

  call cpu_clock_end(id_clock_init)

end subroutine initialize_MOM

! ============================================================================

subroutine register_diags(Time, G, CS)
  type(time_type),           intent(in)    :: Time
  type(ocean_grid_type),     intent(inout) :: G
  type(MOM_control_struct), intent(inout) :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure set up by initialize_MOM.
  character(len=48) :: thickness_units, flux_units, T_flux_units, S_flux_units
  integer :: isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq, nz
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq

  thickness_units = get_thickness_units(G)
  flux_units = get_flux_units(G)
  T_flux_units = get_tr_flux_units(G, "Celsius")
  S_flux_units = get_tr_flux_units(G, "PSU")

  CS%id_u = register_diag_field('ocean_model', 'u', G%axesuL, Time, &
      'Zonal velocity', 'meter second-1')
  CS%id_v = register_diag_field('ocean_model', 'v', G%axesvL, Time, &
      'Meridional velocity', 'meter second-1')
  CS%id_h = register_diag_field('ocean_model', 'h', G%axeshL, Time, &
      'Layer Thickness', thickness_units)
  CS%id_ssh = register_diag_field('ocean_model', 'SSH', G%axesh1, Time, &
      'Sea Surface Height', 'meter', CS%missing)
  CS%id_ssh_inst = register_diag_field('ocean_model', 'SSH_inst', G%axesh1, Time, &
      'Instantaneous Sea Surface Height', 'meter', CS%missing)
  CS%id_ssu = register_diag_field('ocean_model', 'SSU', G%axesu1, Time, &
      'Sea Surface Zonal Velocity', 'meter second-1', CS%missing)
  CS%id_ssv = register_diag_field('ocean_model', 'SSV', G%axesv1, Time, &
      'Sea Surface Meridional Velocity', 'meter second-1', CS%missing)
  CS%id_speed = register_diag_field('ocean_model', 'speed', G%axesh1, Time, &
      'Sea Surface Speed', 'meter second-1', CS%missing)
  if (CS%use_temperature) then
    CS%id_T = register_diag_field('ocean_model', 'temp', G%axeshL, Time, &
        'Potential Temperature', 'Celsius')
    CS%id_S = register_diag_field('ocean_model', 'salt', G%axeshL, Time, &
        'Salinity', 'PSU')
    CS%id_sst = register_diag_field('ocean_model', 'SST', G%axesh1, Time, &
        'Sea Surface Temperature', 'Celsius', CS%missing)
    CS%id_sst_sq = register_diag_field('ocean_model', 'SST_sq', G%axesh1, Time, &
        'Sea Surface Temperature Squared', 'Celsius**2', CS%missing)    
    CS%id_sss = register_diag_field('ocean_model', 'SSS', G%axesh1, Time, &
        'Sea Surface Salinity', 'PSU', CS%missing)
    if (CS%id_sst_sq > 0) call safe_alloc_ptr(CS%SST_sq,isd,ied,jsd,jed)    
  endif
  if (CS%use_temperature .and. CS%use_frazil) then
    CS%id_fraz = register_diag_field('ocean_model', 'frazil', G%axesh1, Time, &
          'Heat sink from frazil formation', 'Watt meter-2')
  endif

  CS%id_salt_deficit = register_diag_field('ocean_model', 'salt_deficit', G%axesh1, Time, &
         'Salt sink in ocean due to ice flux', 'g Salt meter-2 s-1')
  CS%id_Heat_PmE = register_diag_field('ocean_model', 'Heat_PmE', G%axesh1, Time, &
         'Heat flux into ocean from mass flux into ocean', 'Watt meter-2')
  CS%id_intern_heat = register_diag_field('ocean_model', 'internal_heat', G%axesh1, Time, &
         'Heat flux into ocean from geothermal or other internal sources', &
         'Watt meter-2')

  CS%id_CAu = register_diag_field('ocean_model', 'CAu', G%axesuL, Time, &
      'Zonal Coriolis and Advective Acceleration', 'meter second-2')
  CS%id_CAv = register_diag_field('ocean_model', 'CAv', G%axesvL, Time, &
      'Meridional Coriolis and Advective Acceleration', 'meter second-2')
  CS%id_PFu = register_diag_field('ocean_model', 'PFu', G%axesuL, Time, &
      'Zonal Pressure Force Acceleration', 'meter second-2')
  CS%id_PFv = register_diag_field('ocean_model', 'PFv', G%axesvL, Time, &
      'Meridional Pressure Force Acceleration', 'meter second-2')
  if (CS%split) then
    if (CS%id_PFu > 0) call safe_alloc_ptr(CS%diag%PFu_tot,Isdq,Iedq,jsd,jed,nz)
    if (CS%id_PFv > 0) call safe_alloc_ptr(CS%diag%PFv_tot,isd,ied,Jsdq,Jedq,nz)
    if (CS%id_CAu > 0) call safe_alloc_ptr(CS%diag%CAu_tot,Isdq,Iedq,jsd,jed,nz)
    if (CS%id_CAv > 0) call safe_alloc_ptr(CS%diag%CAv_tot,isd,ied,Jsdq,Jedq,nz)
  endif
  CS%id_u_BT_accel = register_diag_field('ocean_model', 'u_BT_accel', G%axesul, Time, &
    'Barotropic Anomaly Zonal Acceleration', 'meter second-1')
  CS%id_v_BT_accel = register_diag_field('ocean_model', 'v_BT_accel', G%axesvl, Time, &
    'Barotropic Anomaly Meridional Acceleration', 'meter second-1')

  CS%id_Tadx = register_diag_field('ocean_model', 'T_adx', G%axesul, Time, &
      'Advective Zonal Flux of Potential Temperature', T_flux_units)
  CS%id_Tady = register_diag_field('ocean_model', 'T_ady', G%axesvl, Time, &
      'Advective Meridional Flux of Potential Temperature', T_flux_units)
  CS%id_Tdiffx = register_diag_field('ocean_model', 'T_diffx', G%axesul, Time, &
      'Diffusive Zonal Flux of Potential Temperature', T_flux_units)
  CS%id_Tdiffy = register_diag_field('ocean_model', 'T_diffy', G%axesvl, Time, &
      'Diffusive Meridional Flux of Potential Temperature', T_flux_units)
  if (CS%id_Tadx > 0)   call safe_alloc_ptr(CS%T_adx,Isdq,Iedq,jsd,jed,nz)
  if (CS%id_Tady > 0)   call safe_alloc_ptr(CS%T_ady,isd,ied,Jsdq,Jedq,nz)
  if (CS%id_Tdiffx > 0) call safe_alloc_ptr(CS%T_diffx,Isdq,Iedq,jsd,jed,nz)
  if (CS%id_Tdiffy > 0) call safe_alloc_ptr(CS%T_diffy,isd,ied,Jsdq,Jedq,nz)

  CS%id_Sadx = register_diag_field('ocean_model', 'S_adx', G%axesul, Time, &
      'Advective Zonal Flux of Salinity', S_flux_units)
  CS%id_Sady = register_diag_field('ocean_model', 'S_ady', G%axesvl, Time, &
      'Advective Meridional Flux of Salinity', S_flux_units)
  CS%id_Sdiffx = register_diag_field('ocean_model', 'S_diffx', G%axesul, Time, &
      'Diffusive Zonal Flux of Salinity', S_flux_units)
  CS%id_Sdiffy = register_diag_field('ocean_model', 'S_diffy', G%axesvl, Time, &
      'Diffusive Meridional Flux of Salinity', S_flux_units)
  if (CS%id_Sadx > 0)   call safe_alloc_ptr(CS%S_adx,Isdq,Iedq,jsd,jed,nz)
  if (CS%id_Sady > 0)   call safe_alloc_ptr(CS%S_ady,isd,ied,Jsdq,Jedq,nz)
  if (CS%id_Sdiffx > 0) call safe_alloc_ptr(CS%S_diffx,Isdq,Iedq,jsd,jed,nz)
  if (CS%id_Sdiffy > 0) call safe_alloc_ptr(CS%S_diffy,isd,ied,Jsdq,Jedq,nz)

  CS%id_Tadx_2d = register_diag_field('ocean_model', 'T_adx_2d', G%axesu1, Time, &
      'Vertically Integrated Advective Zonal Flux of Potential Temperature', T_flux_units)
  CS%id_Tady_2d = register_diag_field('ocean_model', 'T_ady_2d', G%axesv1, Time, &
      'Vertically Integrated Advective Meridional Flux of Potential Temperature', T_flux_units)
  CS%id_Tdiffx_2d = register_diag_field('ocean_model', 'T_diffx_2d', G%axesu1, Time, &
      'Vertically Integrated Diffusive Zonal Flux of Potential Temperature', T_flux_units)
  CS%id_Tdiffy_2d = register_diag_field('ocean_model', 'T_diffy_2d', G%axesv1, Time, &
      'Vertically Integrated Diffusive Meridional Flux of Potential Temperature', T_flux_units)
  if (CS%id_Tadx_2d > 0)   call safe_alloc_ptr(CS%T_adx_2d,Isdq,Iedq,jsd,jed)
  if (CS%id_Tady_2d > 0)   call safe_alloc_ptr(CS%T_ady_2d,isd,ied,Jsdq,Jedq)
  if (CS%id_Tdiffx_2d > 0) call safe_alloc_ptr(CS%T_diffx_2d,Isdq,Iedq,jsd,jed)
  if (CS%id_Tdiffy_2d > 0) call safe_alloc_ptr(CS%T_diffy_2d,isd,ied,Jsdq,Jedq)

  CS%id_Sadx_2d = register_diag_field('ocean_model', 'S_adx_2d', G%axesu1, Time, &
      'Vertically Integrated Advective Zonal Flux of Salinity', S_flux_units)
  CS%id_Sady_2d = register_diag_field('ocean_model', 'S_ady_2d', G%axesv1, Time, &
      'Vertically Integrated Advective Meridional Flux of Salinity', S_flux_units)
  CS%id_Sdiffx_2d = register_diag_field('ocean_model', 'S_diffx_2d', G%axesu1, Time, &
      'Vertically Integrated Diffusive Zonal Flux of Salinity', S_flux_units)
  CS%id_Sdiffy_2d = register_diag_field('ocean_model', 'S_diffy_2d', G%axesv1, Time, &
      'Vertically Integrated Diffusive Meridional Flux of Salinity', S_flux_units)
  if (CS%id_Sadx_2d > 0)   call safe_alloc_ptr(CS%S_adx_2d,Isdq,Iedq,jsd,jed)
  if (CS%id_Sady_2d > 0)   call safe_alloc_ptr(CS%S_ady_2d,isd,ied,Jsdq,Jedq)
  if (CS%id_Sdiffx_2d > 0) call safe_alloc_ptr(CS%S_diffx_2d,Isdq,Iedq,jsd,jed)
  if (CS%id_Sdiffy_2d > 0) call safe_alloc_ptr(CS%S_diffy_2d,isd,ied,Jsdq,Jedq)

  CS%id_uh = register_diag_field('ocean_model', 'uh', G%axesul, Time, &
      'Zonal Thickness Flux', flux_units)
  CS%id_vh = register_diag_field('ocean_model', 'vh', G%axesvl, Time, &
      'Meridional Thickness Flux', flux_units)
  CS%id_uav = register_diag_field('ocean_model', 'uav', G%axesul, Time, &
      'Barotropic-step Averaged Zonal Velocity', 'meter second-1')
  CS%id_vav = register_diag_field('ocean_model', 'vav', G%axesvl, Time, &
      'Barotropic-step Averaged Meridional Velocity', 'meter second-1')

  if (CS%debug_truncations) then
    if (CS%split .and. CS%flux_BT_coupling) then
      call safe_alloc_ptr(CS%diag%du_adj,Isdq,Iedq,jsd,jed,nz)
      call safe_alloc_ptr(CS%diag%dv_adj,isd,ied,Jsdq,Jedq,nz)
    endif
    if (CS%split .and. CS%flux_BT_coupling .and. CS%readjust_BT_trans) then
      call safe_alloc_ptr(CS%diag%du_adj2,Isdq,Iedq,jsd,jed,nz)
      call safe_alloc_ptr(CS%diag%dv_adj2,isd,ied,Jsdq,Jedq,nz)
    endif
    call safe_alloc_ptr(CS%diag%du_dt_visc,Isdq,Iedq,jsd,jed,nz)
    call safe_alloc_ptr(CS%diag%dv_dt_visc,isd,ied,Jsdq,Jedq,nz)
    if (.not.CS%adiabatic) then
      call safe_alloc_ptr(CS%diag%du_dt_dia,Isdq,Iedq,jsd,jed,nz)
      call safe_alloc_ptr(CS%diag%dv_dt_dia,isd,ied,Jsdq,Jedq,nz)
    endif
  endif

  if (CS%split .and. CS%flux_BT_coupling) then
    CS%id_du_adj = register_diag_field('ocean_model', 'du_adj', G%axesuL, Time, &
        'Zonal velocity Adjustment 1', 'meter second-1')
    CS%id_dv_adj = register_diag_field('ocean_model', 'dv_adj', G%axesvL, Time, &
        'Meridional velocity Adjustment 1', 'meter second-1')
    if (CS%id_du_adj > 0) call safe_alloc_ptr(CS%diag%du_adj,Isdq,Iedq,jsd,jed,nz)
    if (CS%id_dv_adj > 0) call safe_alloc_ptr(CS%diag%dv_adj,isd,ied,Jsdq,Jedq,nz)
    if (CS%readjust_BT_trans) then
      CS%id_du_adj2 = register_diag_field('ocean_model', 'du_adj2', G%axesuL, Time, &
          'Zonal velocity Adjustment 2', 'meter second-1')
      CS%id_dv_adj2 = register_diag_field('ocean_model', 'dv_adj2', G%axesvL, Time, &
          'Meridional velocity Adjustment 2', 'meter second-1')
      if (CS%id_du_adj2 > 0) call safe_alloc_ptr(CS%diag%du_adj2,Isdq,Iedq,jsd,jed,nz)
      if (CS%id_dv_adj2 > 0) call safe_alloc_ptr(CS%diag%dv_adj2,isd,ied,Jsdq,Jedq,nz)
    endif
    if (CS%BT_include_udhdt) then
      CS%id_h_dudt = register_diag_field('ocean_model', 'BT_u_dhdt', G%axesu1, Time, &
          'Barotropic zonal transport tendancy from sum(h du_dt)', 'meter3 second-2')
      CS%id_h_dvdt = register_diag_field('ocean_model', 'BT_v_dhdt', G%axesv1, Time, &
          'Barotropic meridional transport tendancy from sum(h du_dt)', 'meter3 second-2')
    endif
  endif


  CS%id_u_predia = register_diag_field('ocean_model', 'u_predia', G%axesuL, Time, &
      'Zonal velocity', 'meter second-1')
  CS%id_v_predia = register_diag_field('ocean_model', 'v_predia', G%axesvL, Time, &
      'Meridional velocity', 'meter second-1')
  CS%id_h_predia = register_diag_field('ocean_model', 'h_predia', G%axeshL, Time, &
      'Layer Thickness', thickness_units)
  CS%id_e_predia = register_diag_field('ocean_model', 'e_predia', G%axeshi, Time, &
      'Interface Heights', 'meter')
  if (CS%use_temperature) then
    CS%id_T_predia = register_diag_field('ocean_model', 'temp_predia', G%axeshL, Time, &
        'Potential Temperature', 'Celsius')
    CS%id_S_predia = register_diag_field('ocean_model', 'salt_predia', G%axeshL, Time, &
        'Salinity', 'PSU')
  endif


end subroutine register_diags

subroutine MOM_timing_init(CS)
  type(MOM_control_struct), intent(in)    :: CS
! Arguments: CS - The control structure set up by initialize_MOM.
  ! This subroutine sets up clock IDs for timing various subroutines.

 id_clock_ocean = cpu_clock_id('Ocean', grain=CLOCK_COMPONENT)
 id_clock_dynamics = cpu_clock_id('Ocean dynamics', grain=CLOCK_SUBCOMPONENT)
 id_clock_thermo = cpu_clock_id('Ocean thermodynamics and tracers', grain=CLOCK_SUBCOMPONENT)
 id_clock_tracer = cpu_clock_id('(Ocean tracer advection)', grain=CLOCK_MODULE_DRIVER)
 if (.not.CS%adiabatic) &
   id_clock_diabatic = cpu_clock_id('(Ocean diabatic driver)', grain=CLOCK_MODULE_DRIVER)

 id_clock_continuity = cpu_clock_id('(Ocean continuity equation *)', grain=CLOCK_MODULE)
 id_clock_pass = cpu_clock_id('(Ocean message passing *)', grain=CLOCK_MODULE)
 id_clock_MOM_init = cpu_clock_id('(Ocean MOM_initialize)', grain=CLOCK_MODULE)
 id_clock_pass_init = cpu_clock_id('(Ocean init message passing *)', grain=CLOCK_ROUTINE)
 if (CS%thickness_diffuse) &
   id_clock_thick_diff = cpu_clock_id('(Ocean thickness diffusion *)', grain=CLOCK_MODULE)
 id_clock_diagnostics = cpu_clock_id('(Ocean collective diagnostics)', grain=CLOCK_MODULE)
 id_clock_Z_diag = cpu_clock_id('(Ocean Z-space diagnostics)', grain=CLOCK_MODULE)

end subroutine MOM_timing_init

! ============================================================================

subroutine write_static_fields(G, diag)
  type(ocean_grid_type),   intent(in) :: G
  type(diag_ptrs), target, intent(in) :: diag
!   This subroutine offers the static fields in the ocean grid type
! for output via the diag_manager.
! Arguments: G - The ocean's grid structure.  Effectively intent in.
!  (in)      diag - A structure containing pointers to common diagnostic fields.

  ! The out_X arrays are needed because some of the elements of the grid
  ! type may be reduced rank macros.
  real :: out_h(SZI_(G),SZJ_(G))
  integer :: id, i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  out_h(:,:) = 0.0

  id = register_static_field('ocean_model', 'geolat', G%axesh1, &
        'Latitude of tracer (h) points', 'degrees_N')
  if (id > 0) call post_data(id, G%geolath, diag, .true.)

  id = register_static_field('ocean_model', 'geolon', G%axesh1, &
        'Longitude of tracer (h) points', 'degrees_E')
  if (id > 0) call post_data(id, G%geolonh, diag, .true.)

  id = register_static_field('ocean_model', 'geolat_c', G%axesq1, &
        'Latitude of corner (q) points', 'degrees_N')
  if (id > 0) call post_data(id, G%geolatq, diag, .true.)

  id = register_static_field('ocean_model', 'geolon_c', G%axesq1, &
        'Longitude of corner (q) points', 'degrees_E')
  if (id > 0) call post_data(id, G%geolonq, diag, .true.)

  id = register_static_field('ocean_model', 'geolat_v', G%axesv1, &
        'Latitude of meridional velocity (v) points', 'degrees_N')
  if (id > 0) call post_data(id, G%geolatv, diag, .true.)

  id = register_static_field('ocean_model', 'geolon_v', G%axesv1, &
        'Longitude of meridional velocity (v) points', 'degrees_E')
  if (id > 0) call post_data(id, G%geolonv, diag, .true.)

  id = register_static_field('ocean_model', 'geolat_u', G%axesu1, &
        'Latitude of zonal velocity (u) points', 'degrees_N')
  if (id > 0) call post_data(id, G%geolatu, diag, .true.)

  id = register_static_field('ocean_model', 'geolon_u', G%axesu1, &
        'Longitude of zonal velocity (u) points', 'degrees_E')
  if (id > 0) call post_data(id, G%geolonu, diag, .true.)

  id = register_static_field('ocean_model', 'area_t', G%axesh1, &
        'Surface area of tracer (h) cells', 'degrees_E')
  if (id > 0) then
    do j=js,je ; do i=is,ie ; out_h(i,j) = G%DXDYh(i,j) ; enddo ; enddo
    call post_data(id, out_h, diag, .true.)
  endif

  id = register_static_field('ocean_model', 'depth_ocean', G%axesh1, &
        'Depth of the ocean at tracer points', 'm', &
        standard_name='sea_floor_depth_below_geoid')
  if (id > 0) call post_data(id, G%D, diag, .true.)

  id = register_static_field('ocean_model', 'wet', G%axesh1, &
        '0 if land, 1 if ocean at tracer points', 'none')
  if (id > 0) call post_data(id, G%hmask, diag, .true.)

  id = register_static_field('ocean_model', 'wet_c', G%axesq1, &
        '0 if land, 1 if ocean at corner (q) points', 'none')
  if (id > 0) call post_data(id, G%qmask, diag, .true.)

  id = register_static_field('ocean_model', 'wet_u', G%axesu1, &
        '0 if land, 1 if ocean at zonal velocity (u) points', 'none')
  if (id > 0) call post_data(id, G%umask, diag, .true.)

  id = register_static_field('ocean_model', 'wet_v', G%axesv1, &
        '0 if land, 1 if ocean at meridional velocity (v) points', 'none')
  if (id > 0) call post_data(id, G%vmask, diag, .true.)

  id = register_static_field('ocean_model', 'Coriolis', G%axesq1, &
        'Coriolis parameter at corner (q) points', 's-1')
  if (id > 0) call post_data(id, G%f, diag, .true.)

end subroutine write_static_fields

! ============================================================================

subroutine set_restart_fields(grid, param_file, CS)
  type(ocean_grid_type),    intent(in) :: grid
  type(param_file_type),    intent(in) :: param_file
  type(MOM_control_struct), intent(in) :: CS
!   Set the fields that are needed for bitwise identical restarting
! the time stepping scheme.  In addition to those specified here
! directly, there may be fields related to the forcing or to the
! barotropic solver that are needed; these are specified in sub-
! routines that are called from this one.
!   This routine should be altered if there are any changes to the
! time stepping scheme.  The CHECK_RESTART facility may be used to
! confirm that all needed restart fields have been included.
!
! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      CS - The control structure set up by initialize_MOM.
  type(vardesc) :: vd
  character(len=48) :: thickness_units, flux_units

  thickness_units = get_thickness_units(grid)
  flux_units = get_flux_units(grid)

  vd = vardesc("u","Zonal velocity",'u','L','s',"meter second-1", 'd')
  call register_restart_field(CS%u(:,:,:,1), CS%u(:,:,:,2), vd, .true., CS%restart_CSp)

  vd = vardesc("v","Meridional velocity",'v','L','s',"meter second-1", 'd')
  call register_restart_field(CS%v(:,:,:,1), CS%v(:,:,:,2), vd, .true., CS%restart_CSp)

  vd = vardesc("h","Layer Thickness",'h','L','s',thickness_units, 'd')
  call register_restart_field(CS%h(:,:,:,1), CS%h(:,:,:,2), vd, .true., CS%restart_CSp)

  if (CS%use_temperature) then
    vd = vardesc("Temp","Potential Temperature",'h','L','s',"degC", 'd')
    call register_restart_field(CS%tv%T, CS%tv%T, vd, .true., CS%restart_CSp)

    vd = vardesc("Salt","Salinity",'h','L','s',"PSU", 'd')
    call register_restart_field(CS%tv%S, CS%tv%S, vd, .true., CS%restart_CSp)
  endif

  if (CS%use_frazil) then
    vd = vardesc("frazil","Frazil heat flux into ocean",'h','1','s',"J m-2", 'd')
    call register_restart_field(CS%tv%frazil, CS%tv%frazil, vd, .false., CS%restart_CSp)
  endif

  if (CS%interp_p_surf) then
    vd = vardesc("p_surf_prev","Previous ocean surface pressure",'h','1','s',"Pa", 'd')
    call register_restart_field(CS%p_surf_prev, CS%p_surf_prev, vd, .false., CS%restart_CSp)
  endif

  vd = vardesc("ave_ssh","Time average sea surface height",'h','1','s',"meter", 'd')
  call register_restart_field(CS%ave_ssh, CS%ave_ssh, vd, .false., CS%restart_CSp)

end subroutine set_restart_fields

! ============================================================================

subroutine calculate_surface_state(state, u, v, h, ssh, G, CS, p_atm)
  type(surface),                                  intent(inout) :: state
  real, target, dimension(NXMEMQ_,NYMEM_,NKMEM_), intent(in)    :: u
  real, target, dimension(NXMEM_,NYMEMQ_,NKMEM_), intent(in)    :: v
  real, target, dimension(NXMEM_,NYMEM_,NKMEM_),  intent(in)    :: h
  real, target, dimension(NXMEM_,NYMEM_),         intent(inout) :: ssh
  type(ocean_grid_type),                          intent(inout) :: G
  type(MOM_control_struct),                       intent(in)    :: CS
  real, optional, pointer, dimension(:,:)                       :: p_atm
!   This subroutine sets the surface (return) properties of the ocean
! model by setting the appropriate pointers in state.  Unused fields
! are set to NULL.
!
! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m.
!  (in/out)  ssh - Time mean sea surface hieght, in m.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure set up by initialize_MOM.
!  (in)      p_atm - The atmospheric pressure, in Pa.
!  (out)     state - A structure containing fields that describe the
!                    surface state of the ocean.
  real :: depth(SZI_(G))    ! The distance from the surface, in m.
  real :: depth_ml          ! The depth over which to average to
                            ! determine mixed layer properties, m.
  real :: dh                ! The thickness of a layer that is within the
                            ! mixed layer, in m.
  real :: mass              ! The mass per unit area of a layer, in kg m-2.
  real :: IgR0
  integer :: i, j, k, is, ie, js, je, nz, num_pnts, num_errs
  integer :: isd, ied, jsd, jed
  character(128) :: msg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  state%sea_lev => ssh

  if (present(p_atm)) then ; if (ASSOCIATED(p_atm)) then
    IgR0 = 1.0 / (G%Rho0 * G%g_Earth)
    do j=js,je ; do i=is,ie
      ssh(i,j) = ssh(i,j) + p_atm(i,j) * IgR0
    enddo ; enddo
  endif ; endif

  if (CS%smooth_ssh_passes > 0.0) then
    call smooth_SSH(ssh, G, CS%smooth_ssh_passes)
  endif

  if (CS%bulkmixedlayer) then
    state%SST => CS%tv%T(:,:,1)
    state%SSS => CS%tv%S(:,:,1)
    nullify(state%sfc_density)
    if (associated(CS%tv%Hml)) state%Hml => CS%tv%Hml
  else
    if (CS%use_temperature) then
      if (.not.associated(state%SST)) then
        allocate(state%SST(isd:ied,jsd:jed)) ; state%SST(:,:) = 0.0
      endif
      if (.not.associated(state%SSS)) then
        allocate(state%SSS(isd:ied,jsd:jed)) ; state%SSS(:,:) = 0.0
      endif
      nullify(state%sfc_density)
    else
      if (.not.associated(state%sfc_density)) then
        allocate(state%sfc_density(isd:ied,jsd:jed)) ; state%sfc_density(:,:) = 0.0
      endif
      nullify(state%SST) ; nullify(state%SSS)
    endif
    if (.not.associated(state%Hml)) allocate(state%Hml(isd:ied,jsd:jed))

    depth_ml = CS%Hmix
  !   Determine the mean properties of the uppermost depth_ml fluid.
    do j=js,je
      do i=is,ie
        depth(i) = 0.0
        if (CS%use_temperature) then
          state%SST(i,j) = 0.0 ; state%SSS(i,j) = 0.0
        else
          state%sfc_density(i,j) = 0.0
        endif
      enddo

      do k=1,nz ; do i=is,ie
        if (depth(i) + h(i,j,k) < depth_ml) then
          dh = h(i,j,k)
        elseif (depth(i) < depth_ml) then
          dh = depth_ml - depth(i)
        else
          dh = 0.0
        endif
        if (CS%use_temperature) then
          state%SST(i,j) = state%SST(i,j) + dh * CS%tv%T(i,j,k)
          state%SSS(i,j) = state%SSS(i,j) + dh * CS%tv%S(i,j,k)
        else
          state%sfc_density(i,j) = state%sfc_density(i,j) + dh * G%Rlay(k)
        endif
        depth(i) = depth(i) + dh
      enddo ; enddo
  ! Calculate the average properties of the mixed layer depth.
      do i=is,ie
        if (depth(i) < G%H_subroundoff) depth(i) = G%H_subroundoff
        if (CS%use_temperature) then
          state%SST(i,j) = state%SST(i,j) / depth(i)
          state%SSS(i,j) = state%SSS(i,j) / depth(i)
        else
          state%sfc_density(i,j) = state%sfc_density(i,j) / depth(i)
        endif
        state%Hml(:,:) = depth(i)
      enddo
    enddo ! end of j loop
  endif                                             ! end BULKMIXEDLAYER

  state%u => u(:,:,1)
  state%v => v(:,:,1)
  state%frazil => CS%tv%frazil
  state%TempxPmE => CS%tv%TempxPmE
  state%internal_heat => CS%tv%internal_heat

  if (associated(state%salt_deficit) .and. associated(CS%tv%salt_deficit)) then
    do j=js,je ; do i=is,ie
      ! Convert from gSalt to kgSalt
      state%salt_deficit(i,j) = 1000.0 * CS%tv%salt_deficit(i,j)
    enddo ; enddo
  endif

  ! Allocate structures for ocean_mass, ocean_heat, and ocean_salt.  This could
  ! be wrapped in a run-time flag to disable it for economy, since the 3-d
  ! sums are not negligible.
  if (.not.associated(state%ocean_mass)) then
    allocate(state%ocean_mass(isd:ied,jsd:jed)) ; state%ocean_mass(:,:) = 0.0
  endif
  if (CS%use_temperature) then
    if (.not.associated(state%ocean_heat)) then
      allocate(state%ocean_heat(isd:ied,jsd:jed)) ; state%ocean_heat(:,:) = 0.0
    endif
    if (.not.associated(state%ocean_salt)) then
      allocate(state%ocean_salt(isd:ied,jsd:jed)) ; state%ocean_salt(:,:) = 0.0
    endif
  endif

  if (associated(state%ocean_mass) .and. associated(state%ocean_heat) .and. &
      associated(state%ocean_salt)) then
    do j=js,je ; do i=is,ie
      state%ocean_mass(i,j) = 0.0
      state%ocean_heat(i,j) = 0.0 ; state%ocean_salt(i,j) = 0.0
    enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      mass = G%H_to_kg_m2*h(i,j,k)
      state%ocean_mass(i,j) = state%ocean_mass(i,j) + mass
      state%ocean_heat(i,j) = state%ocean_heat(i,j) + mass*CS%tv%T(i,j,k)
      state%ocean_salt(i,j) = state%ocean_salt(i,j) + &
                              mass * (1.0e-3*CS%tv%S(i,j,k))
    enddo ; enddo ; enddo
  else
    if (associated(state%ocean_mass)) then
      do j=js,je ; do i=is,ie ; state%ocean_mass(i,j) = 0.0 ; enddo ; enddo
      do k=1,nz ; do j=js,je ; do i=is,ie
        state%ocean_mass(i,j) = state%ocean_mass(i,j) + G%H_to_kg_m2*h(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (associated(state%ocean_heat)) then
      do j=js,je ; do i=is,ie ; state%ocean_heat(i,j) = 0.0 ; enddo ; enddo
      do k=1,nz ; do j=js,je ; do i=is,ie
        mass = G%H_to_kg_m2*h(i,j,k)
        state%ocean_heat(i,j) = state%ocean_heat(i,j) + mass*CS%tv%T(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (associated(state%ocean_salt)) then
      do j=js,je ; do i=is,ie ; state%ocean_salt(i,j) = 0.0 ; enddo ; enddo
      do k=1,nz ; do j=js,je ; do i=is,ie
        mass = G%H_to_kg_m2*h(i,j,k)
        state%ocean_salt(i,j) = state%ocean_salt(i,j) + &
                                mass * (1.0e-3*CS%tv%S(i,j,k))
      enddo ; enddo ; enddo
    endif
  endif

  if (associated(CS%visc%taux_shelf)) state%taux_shelf => CS%visc%taux_shelf
  if (associated(CS%visc%tauy_shelf)) state%tauy_shelf => CS%visc%tauy_shelf

  if (associated(CS%tracer_flow_CSp)) then
    if (.not.associated(state%tr_fields)) allocate(state%tr_fields)
    call call_tracer_surface_state(state, h, G, CS%tracer_flow_CSp)
  endif

  if (CS%check_bad_surface_vals) then
    num_errs=0 ! count number of errors
    num_pnts=0 ! count number of errors
    do j=js,je; do i=is,ie
      if (num_errs>99) exit ! If things get really bad, stop filling up the tty
      k=0 ! Num errors at this point
      if (G%hmask(i,j)>0.) then
        if (state%sea_lev(i,j)<=-G%D(i,j)) then
          k=k+1
          write(msg(1:128),'(2(a,i4,x),2(a,f8.3,x),a,es12.3)') &
             'Sea level < bathymetry at i=',i,'j=',j,'x=',G%geolonh(i,j),'y=',G%geolath(i,j),'SSH=',state%sea_lev(i,j)
          call MOM_error(WARNING, trim(msg), all_print=.true.)
        endif
        if (state%sea_lev(i,j)>=CS%bad_val_ssh_max) then
          k=k+1
          write(msg(1:128),'(2(a,i4,x),2(a,f8.3,x),a,es12.3)') &
             'Very high sea level at i=',i,'j=',j,'x=',G%geolonh(i,j),'y=',G%geolath(i,j),'SSH=',state%sea_lev(i,j)
          call MOM_error(WARNING, trim(msg), all_print=.true.)
        endif
        if (CS%use_temperature) then
          if (state%SSS(i,j)<0.) then
            k=k+1
            write(msg(1:128),'(2(a,i4,x),2(a,f8.3,x),a,es12.3)') &
               'Negative salinity at i=',i,'j=',j,'x=',G%geolonh(i,j),'y=',G%geolath(i,j),'SSS=',state%SSS(i,j)
            call MOM_error(WARNING, trim(msg), all_print=.true.)
          elseif (state%SSS(i,j)>=CS%bad_val_sss_max) then
            k=k+1
            write(msg(1:128),'(2(a,i4,x),2(a,f8.3,x),a,es12.3)') &
               'Very high salinity at i=',i,'j=',j,'x=',G%geolonh(i,j),'y=',G%geolath(i,j),'SSS=',state%SSS(i,j)
            call MOM_error(WARNING, trim(msg), all_print=.true.)
          endif
          if (state%SST(i,j)<CS%bad_val_sst_min) then
            k=k+1
            write(msg(1:128),'(2(a,i4,x),2(a,f8.3,x),a,es12.3)') &
               'Very cold SST at i=',i,'j=',j,'x=',G%geolonh(i,j),'y=',G%geolath(i,j),'SST=',state%SST(i,j)
            call MOM_error(WARNING, trim(msg), all_print=.true.)
          elseif (state%SST(i,j)>=CS%bad_val_sst_max) then
            k=k+1
            write(msg(1:128),'(2(a,i4,x),2(a,f8.3,x),a,es12.3)') &
               'Very hot SST at i=',i,'j=',j,'x=',G%geolonh(i,j),'y=',G%geolath(i,j),'SST=',state%SST(i,j)
            call MOM_error(WARNING, trim(msg), all_print=.true.)
          endif
        endif ! use_temperature
        num_pnts=num_pnts+min(1,k)
        num_errs=num_errs+k
      endif ! hmask
    enddo; enddo
    if (num_errs>0) then
      write(msg(1:128),'(3(a,i4,x))') 'There were',num_errs,'errors involving',num_pnts,'points'
      call MOM_error(FATAL, trim(msg))
    endif
  endif

end subroutine calculate_surface_state

! ============================================================================

subroutine smooth_SSH(ssh, G, smooth_passes)
  real, dimension(NXMEM_,NYMEM_),      intent(inout) :: ssh
  type(ocean_grid_type),               intent(inout) :: G
  real,                                intent(in)    :: smooth_passes
!   This subroutine applies a number of 2-D smoothing passes, each of which
! applies a nominal filter with the following weights:
!         1/8
!    1/8  1/2  1/8
!         1/8  
! Arguments: ssh - Time mean sea surface hieght, in m.  (Intent inout.)
!  (in)      G - The ocean's grid structure.
!  (in)      smooth_passes - the number of smoothing passes to apply.
  real, dimension(SZIQ_(G), SZJ_(G)) :: flux_x, area_x
  real, dimension(SZI_(G), SZJQ_(G)) :: flux_y, area_y
  real :: wt
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, isl, iel, jsl, jel, halo
  integer :: pass, tot_pass

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (smooth_passes <= 0.0) return

  tot_pass = ceiling(smooth_passes)
  wt = 0.125 * smooth_passes / real(tot_pass)
  halo = -1

  do j=jsd+1,jed-1 ; do I=isd,ied-1
    area_x(I,j) = min(G%dy_u(I,j)*G%dxu(I,j), G%dxdyh(i,j), G%dxdyh(i+1,j))
  enddo ; enddo

  do J=jsd,jed-1 ; do i=isd+1,ied-1
    area_y(i,J) = min(G%dx_v(i,J)*G%dyv(i,J), G%dxdyh(i,j), G%dxdyh(i,j+1))
  enddo ; enddo

  do pass=1,tot_pass
    if (halo < 0) then
      call pass_var(ssh, G%domain)
      halo = min(is-isd-1, ied-ie-1, js-jsd-1, jed-je-1, tot_pass-pass)
    endif
    isl = is-halo ; iel = ie+halo ; jsl = js-halo ; jel = je+halo

    do j=jsl,jel ; do I=isl-1,iel
      flux_x(I,j) =  (wt * area_x(I,j)) * (ssh(i,j) - ssh(i+1,j))
    enddo ; enddo

    do J=jsl-1,jel ; do i=isl,iel
      flux_y(i,J) =  (wt * area_y(i,J)) * (ssh(i,j) - ssh(i,j+1))
    enddo ; enddo
  
    do j=jsl,jel ; do i=isl,iel
      ssh(i,j) = ssh(i,j) + ((flux_x(I-1,j) - flux_x(I,j)) + &
                             (flux_y(i,J-1) - flux_y(i,J))) * G%Idxdyh(i,j)
    enddo ; enddo

    halo = halo - 1
  enddo

end subroutine smooth_SSH

! ============================================================================

subroutine MOM_end(CS)
  type(MOM_control_struct), pointer      :: CS

  call end_regridding(CS%regridding_opts)

  DEALLOC(CS%u) ; DEALLOC(CS%v) ; DEALLOC(CS%h)
  DEALLOC(CS%uh) ; DEALLOC(CS%vh)
  DEALLOC(CS%diffu) ; DEALLOC(CS%diffv)
  DEALLOC(CS%CAu) ; DEALLOC(CS%CAv)
  DEALLOC(CS%PFu) ; DEALLOC(CS%PFv)
  if (CS%use_temperature) then
    DEALLOC(CS%tv%T) ; DEALLOC(CS%tv%S)
  endif
  if (associated(CS%tv%frazil)) deallocate(CS%tv%frazil)
  if (associated(CS%tv%salt_deficit)) deallocate(CS%tv%salt_deficit)  
  if (associated(CS%tv%Hml)) deallocate(CS%tv%Hml)

  DEALLOC(CS%uhtr) ; DEALLOC(CS%vhtr)
  if (CS%split) then
    DEALLOC(CS%eta) ; DEALLOC(CS%uhbt) ; DEALLOC(CS%vhbt)
    DEALLOC(CS%uhbt_in) ; DEALLOC(CS%vhbt_in)
    DEALLOC(CS%h_av) ; DEALLOC(CS%u_av) ; DEALLOC(CS%v_av)
    DEALLOC(CS%eta_PF) ; DEALLOC(CS%pbce)
    DEALLOC(CS%u_accel_bt) ; DEALLOC(CS%v_accel_bt)
    DEALLOC(CS%visc_rem_u) ; DEALLOC(CS%visc_rem_v)
    call dealloc_BT_cont_type(CS%BT_cont)
  endif
  DEALLOC(CS%ave_ssh)

  deallocate(CS)

end subroutine MOM_end

end module MOM
