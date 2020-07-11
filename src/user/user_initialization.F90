!> A template of a user to code up customized initial conditions.
module user_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N
use MOM_open_boundary, only : OBC_DIRECTION_S
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
implicit none ; private

#include <MOM_memory.h>

public USER_set_coord, USER_initialize_topography, USER_initialize_thickness
public USER_initialize_velocity, USER_init_temperature_salinity
public USER_initialize_sponges, USER_set_OBC_data, USER_set_rotation

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> A module variable that should not be used.
!! \todo Move this module variable into a control structure.
logical :: first_call = .true.

contains

!> Set vertical coordinates.
subroutine USER_set_coord(Rlay, g_prime, GV, US, param_file, eqn_of_state)
  type(verticalGrid_type), intent(in)  :: GV         !< The ocean's vertical grid
                                                     !! structure.
  real, dimension(:),      intent(out) :: Rlay       !< Layer potential density [R ~> kg m-3].
  real, dimension(:),      intent(out) :: g_prime    !< The reduced gravity at
                                                     !! each interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),   intent(in)  :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file !< A structure indicating the
                                                     !! open file to parse for model
                                                     !! parameter values.
  type(EOS_type),          pointer     :: eqn_of_state !< Integer that selects the
                                                     !! equation of state.

  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_set_coord: " // &
    "Unmodified user routine called - you must edit the routine to use it")
  Rlay(:) = 0.0
  g_prime(:) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_set_coord

!> Initialize topography.
subroutine USER_initialize_topography(D, G, param_file, max_depth, US)
  type(dyn_horgrid_type),          intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                   intent(out) :: D !< Ocean bottom depth in m or Z if US is present
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth in the units of D
  type(unit_scale_type), optional, intent(in)  :: US !< A dimensional unit scaling type

  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_initialize_topography: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  D(:,:) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_topography

!> initialize thicknesses.
subroutine USER_initialize_thickness(h, G, GV, param_file, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h  !< The thicknesses being initialized [H ~> m or kg m-2].
  type(param_file_type),   intent(in)  :: param_file !< A structure indicating the open
                                             !! file to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                             !! only read parameters without changing h.

  logical :: just_read    ! If true, just read parameters but set nothing.

  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_initialize_thickness: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (just_read) return ! All run-time parameters have been read, so return.

  h(:,:,1) = 0.0 ! h should be set [H ~> m or kg m-2].

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_thickness

!> initialize velocities.
subroutine USER_initialize_velocity(u, v, G, US, param_file, just_read_params)
  type(ocean_grid_type),                       intent(in)  :: G !< Ocean grid structure.
  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)), intent(out) :: u !< i-component of velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G), SZJB_(G), SZK_(G)), intent(out) :: v !< j-component of velocity [L T-1 ~> m s-1]
  type(unit_scale_type),                       intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                       intent(in)  :: param_file !< A structure indicating the
                                                            !! open file to parse for model
                                                            !! parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.

  logical :: just_read    ! If true, just read parameters but set nothing.

  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_initialize_velocity: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (just_read) return ! All run-time parameters have been read, so return.

  u(:,:,1) = 0.0
  v(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_velocity

!> This function puts the initial layer temperatures and salinities
!! into T(:,:,:) and S(:,:,:).
subroutine USER_init_temperature_salinity(T, S, G, param_file, eqn_of_state, just_read_params)
  type(ocean_grid_type),                     intent(in)  :: G !< Ocean grid structure.
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T !< Potential temperature [degC].
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: S !< Salinity [ppt].
  type(param_file_type),                     intent(in)  :: param_file !< A structure indicating the
                                                            !! open file to parse for model
                                                            !! parameter values.
  type(EOS_type),                            pointer     :: eqn_of_state !< Integer that selects the
                                                            !! equation of state.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will only
                                                           !! read parameters without changing T & S.

  logical :: just_read    ! If true, just read parameters but set nothing.

  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_init_temperature_salinity: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (just_read) return ! All run-time parameters have been read, so return.

  T(:,:,1) = 0.0
  S(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_init_temperature_salinity

!> Set up the sponges.
subroutine USER_initialize_sponges(G, GV, use_temp, tv, param_file, CSp, h)
  type(ocean_grid_type),   intent(in) :: G             !< Ocean grid structure.
  type(verticalGrid_type), intent(in) :: GV            !< The ocean's vertical grid structure.
  logical,                 intent(in) :: use_temp      !< If true, temperature and salinity are state variables.
  type(thermo_var_ptrs),   intent(in) :: tv            !< A structure containing pointers
                                                       !! to any available thermodynamic
                                                       !! fields, potential temperature and
                                                       !! salinity or mixed layer density.
                                                       !! Absent fields have NULL ptrs.
  type(param_file_type),   intent(in) :: param_file    !< A structure indicating the
                                                       !! open file to parse for model
                                                       !! parameter values.
  type(sponge_CS),         pointer    :: CSp           !< A pointer to the sponge control structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h             !< Layer thicknesses [H ~> m or kg m-2].
  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_initialize_sponges: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_sponges

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine USER_set_OBC_data(OBC, tv, G, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC   !< This open boundary condition type specifies
                                                  !! whether, where, and what open boundary
                                                  !! conditions are used.
  type(thermo_var_ptrs),      intent(in) :: tv    !< A structure containing pointers to any
                                       !! available thermodynamic fields, including potential
                                       !! temperature and salinity or mixed layer density. Absent
                                       !! fields have NULL ptrs.
  type(ocean_grid_type),      intent(in) :: G     !< The ocean's grid structure.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.
  type(tracer_registry_type), pointer    :: tr_Reg !< Tracer registry.
!  call MOM_error(FATAL, &
!   "USER_initialization.F90, USER_set_OBC_data: " // &
!   "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_set_OBC_data

subroutine USER_set_rotation(G, param_file)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(param_file_type), intent(in)    :: param_file !< A structure to parse for run-time parameters
  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_set_rotation: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_set_rotation

!> Write output about the parameter values being used.
subroutine write_user_log(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "user_initialization" ! This module's name.

  call log_version(param_file, mdl, version)
  first_call = .false.

end subroutine write_user_log

!> \namespace user_initialization
!!
!!  This subroutine initializes the fields for the simulations.
!!  The one argument passed to initialize, Time, is set to the
!!  current time of the simulation.  The fields which might be initialized
!!  here are:
!!  - u - Zonal velocity [Z T-1 ~> m s-1].
!!  - v - Meridional velocity [Z T-1 ~> m s-1].
!!  - h - Layer thickness [H ~> m or kg m-2].  (Must be positive.)
!!  - G%bathyT - Basin depth [Z ~> m].  (Must be positive.)
!!  - G%CoriolisBu - The Coriolis parameter [T-1 ~> s-1].
!!  - GV%g_prime - The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2].
!!  - GV%Rlay - Layer potential density (coordinate variable) [R ~> kg m-3].
!!  If ENABLE_THERMODYNAMICS is defined:
!!  - T - Temperature [degC].
!!  - S - Salinity [psu].
!!  If BULKMIXEDLAYER is defined:
!!  - Rml - Mixed layer and buffer layer potential densities [R ~> kg m-3].
!!  If SPONGE is defined:
!!  - A series of subroutine calls are made to set up the damping
!!    rates and reference profiles for all variables that are damped
!!    in the sponge.
!!
!!  Any user provided tracer code is also first linked through this
!!  subroutine.
!!
!!  These variables are all set in the set of subroutines (in this
!!  file) USER_initialize_bottom_depth, USER_initialize_thickness,
!!  USER_initialize_velocity,  USER_initialize_temperature_salinity,
!!  USER_initialize_mixed_layer_density, USER_initialize_sponges,
!!  USER_set_coord, and USER_set_ref_profile.
!!
!!  The names of these subroutines should be self-explanatory. They
!!  start with "USER_" to indicate that they will likely have to be
!!  modified for each simulation to set the initial conditions and
!!  boundary conditions.  Most of these take two arguments: an integer
!!  argument specifying whether the fields are to be calculated
!!  internally or read from a NetCDF file; and a string giving the
!!  path to that file.  If the field is initialized internally, the
!!  path is ignored.

end module user_initialization
