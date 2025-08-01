!> Drives the generic version of tracers TOPAZ and CFC and other GFDL BGC components
module MOM_generic_tracer

! This file is part of MOM6. See LICENSE.md for the license.

#include <MOM_memory.h>

! The following macro is usually defined in <fms_platform.h> but since MOM6 should not directly
! include files from FMS we replicate the macro lines here:
#ifdef NO_F2000
#define _ALLOCATED associated
#else
#define _ALLOCATED allocated
#endif

! ### These imports should not reach into FMS directly ###

use MOM_ALE_sponge, only : ALE_sponge_CS
use MOM_coms, only : EFP_type
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : param_file_type
use MOM_forcing_type, only : forcing, optics_type
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : MOM_restart_CS
use MOM_sponge, only : sponge_CS
use MOM_time_manager, only : time_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface, thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

!> A state hidden in module data that is very much not allowed in MOM6
! ### This needs to be fixed
logical :: g_registered = .false.

public register_MOM_generic_tracer, initialize_MOM_generic_tracer
public MOM_generic_tracer_column_physics, MOM_generic_tracer_surface_state
public end_MOM_generic_tracer, MOM_generic_tracer_get
public MOM_generic_tracer_stock
public MOM_generic_flux_init
public MOM_generic_tracer_min_max
public MOM_generic_tracer_fluxes_accumulate
public register_MOM_generic_tracer_segments

!> Control structure for generic tracers
type, public :: MOM_generic_tracer_CS ; private
  character(len = 200) :: IC_file !< The file in which the generic tracer initial values can
                                  !! be found, or an empty string for internal initialization.
  logical :: Z_IC_file !< If true, the generic_tracer IC_file is in Z-space.  The default is false.
  real :: tracer_IC_val = 0.0    !< The initial value assigned to tracers, in
                                 !! concentration units [conc]
  real :: tracer_land_val = -1.0 !< The values of tracers used where land is masked out, in
                                 !! concentration units [conc]
  logical :: tracers_may_reinit  !< If true, tracers may go through the
                                 !! initialization code if they are not found in the restart files.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                             !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< Restart control structure
  type(ocean_OBC_type), pointer :: OBC => NULL() !<open boundary condition type
  !> Pointer to the first element of the linked list of generic tracers.
  !type(g_tracer_type), pointer :: g_tracer_list => NULL()

end type MOM_generic_tracer_CS

contains

!> Initializes the generic tracer packages and adds their tracers to the list
!! Adds the tracers in the list of generic tracers to the set of MOM tracers (i.e., MOM-register them)
!! Register these tracers for restart
function register_MOM_generic_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
!subroutine register_MOM_generic_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in)   :: HI         !< Horizontal index ranges
  type(verticalGrid_type),    intent(in)   :: GV         !< The ocean's vertical grid structure
  type(param_file_type),      intent(in)   :: param_file !< A structure to parse for run-time parameters
  type(MOM_generic_tracer_CS), pointer     :: CS         !< Pointer to the control structure for this module
  type(tracer_registry_type), pointer      :: tr_Reg     !< Pointer to the control structure for the tracer
                                                         !! advection and diffusion module.
  type(MOM_restart_CS), target, intent(inout)  :: restart_CS !< MOM restart control struct

  logical :: register_MOM_generic_tracer

  register_MOM_generic_tracer = .false.

  call MOM_error(FATAL, "register_MOM_generic_tracer should not be called with the stub code "// &
         "in MOM6/config_src/external, as it does nothing.  Recompile using the full MOM_generic_tracer package.")

end function register_MOM_generic_tracer

!> Register OBC segments for generic tracers
subroutine register_MOM_generic_tracer_segments(CS, GV, OBC, tr_Reg, param_file)
  type(MOM_generic_tracer_CS), pointer    :: CS         !< Pointer to the control structure for this module.
  type(verticalGrid_type),     intent(in) :: GV         !< The ocean's vertical grid structure
  type(ocean_OBC_type),        pointer    :: OBC        !< This open boundary condition type specifies whether,
                                                        !! where, and what open boundary conditions are used.
  type(tracer_registry_type),  pointer    :: tr_Reg     !< Pointer to the control structure for the tracer
                                                        !! advection and diffusion module.
  type(param_file_type),       intent(in) :: param_file !< A structure to parse for run-time parameters

end subroutine register_MOM_generic_tracer_segments

!>  Initialize phase II:  Initialize required variables for generic tracers
!!  There are some steps of initialization that cannot be done in register_MOM_generic_tracer
!!  This is the place and time to do them:
!!      Set the grid mask and initial time for all generic tracers.
!!      Diag_register them.
!!      Z_diag_register them.
!!
!!   This subroutine initializes the NTR tracer fields in tr(:,:,:,:)
!! and it sets up the tracer output.
subroutine initialize_MOM_generic_tracer(restart, day, G, GV, US, h, tv, param_file, diag, OBC, &
                                         CS, sponge_CSp, ALE_sponge_CSp)
  logical,                               intent(in) :: restart !< .true. if the fields have already been
                                                               !! read from a restart file.
  type(time_type), target,               intent(in) :: day     !< Time of the start of the run.
  type(ocean_grid_type),                 intent(inout) :: G    !< The ocean's grid structure
  type(verticalGrid_type),               intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                 intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                 intent(in) :: tv      !< A structure pointing to various thermodynamic
                                                               !! variables
  type(param_file_type),                 intent(in) :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl),               target, intent(in) :: diag    !< Regulates diagnostic output.
  type(ocean_OBC_type),                  pointer    :: OBC     !< This open boundary condition type specifies whether,
                                                               !! where, and what open boundary conditions are used.
  type(MOM_generic_tracer_CS),           pointer    :: CS      !< Pointer to the control structure for this module.
  type(sponge_CS),                       pointer    :: sponge_CSp !< Pointer to the control structure for the sponges.
  type(ALE_sponge_CS),                   pointer    :: ALE_sponge_CSp !< Pointer  to the control structure for the
                                                               !! ALE sponges.

end subroutine initialize_MOM_generic_tracer

!>  Column physics for generic tracers.
!!      Get the coupler values for generic tracers that exchange with atmosphere
!!      Update generic tracer concentration fields from sources and sinks.
!!      Vertically diffuse generic tracer concentration fields.
!!      Update generic tracers from bottom and their bottom reservoir.
!!
!!   This subroutine applies diapycnal diffusion and any other column
!! tracer physics or chemistry to the tracers from this file.
!! CFCs are relatively simple, as they are passive tracers. with only a surface
!! flux as a source.
subroutine MOM_generic_tracer_column_physics(h_old, h_new, ea, eb, fluxes, Hml, dt, G, GV, US, CS, tv, optics, &
      evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G     !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV    !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: ea    !< The amount of fluid entrained from the layer
                                               !! above during this call [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: eb    !< The amount of fluid entrained from the layer
                                               !! below during this call [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                               !! and tracer forcing fields.
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: Hml  !< Mixed layer depth [Z ~> m]
  real,                    intent(in) :: dt    !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US    !< A dimensional unit scaling type
  type(MOM_generic_tracer_CS), pointer :: CS   !< Pointer to the control structure for this module.
  type(thermo_var_ptrs),   intent(in) :: tv    !< A structure pointing to various thermodynamic variables
  type(optics_type),       intent(in) :: optics !< The structure containing optical properties.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                               !! be fluxed out of the top layer in a timestep [nondim]
                                               !   Stored previously in diabatic CS.
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which fluxes
                                               !!  can be applied [H ~> m or kg m-2]
                                               !   Stored previously in diabatic CS.

end subroutine MOM_generic_tracer_column_physics

!> This subroutine calculates mass-weighted integral on the PE either
!! of all available tracer concentrations, or of a tracer that is
!! being requested specifically, returning the number of stocks it has
!! calculated. If the stock_index is present, only the stock corresponding
!! to that coded index is returned.
function MOM_generic_tracer_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< Layer thicknesses [H ~> m or kg m-2]
  type(EFP_type), dimension(:),       intent(out)   :: stocks !< The mass-weighted integrated amount of each
                                                              !! tracer, in kg times concentration units [kg conc]
  type(MOM_generic_tracer_CS),        pointer       :: CS     !< Pointer to the control structure for this module.
  character(len=*), dimension(:),     intent(out)   :: names  !< The names of the stocks calculated.
  character(len=*), dimension(:),     intent(out)   :: units  !< The units of the stocks calculated.
  integer, optional,                  intent(in)    :: stock_index !< The coded index of a specific stock
                                                                   !! being sought.
  integer                                           :: MOM_generic_tracer_stock !< Return value, the
                                                                   !! number of stocks calculated here.

  MOM_generic_tracer_stock = 0

end function MOM_generic_tracer_stock

!> This subroutine finds the global min and max of either of all available
!! tracer concentrations, or of a tracer that is being requested specifically,
!! returning the number of tracers it has evaluated.
!! It also optionally returns the locations of the extrema.
function MOM_generic_tracer_min_max(ind_start, got_minmax, gmin, gmax, G, CS, names, units, &
                                    xgmin, ygmin, zgmin, xgmax, ygmax, zgmax)
  integer,                        intent(in)    :: ind_start !< The index of the tracer to start with
  logical, dimension(:),          intent(out)   :: got_minmax !< Indicates whether the global min and
                                                          !! max are found for each tracer
  real, dimension(:),             intent(out)   :: gmin   !< Global minimum of each tracer [conc]
  real, dimension(:),             intent(out)   :: gmax   !< Global maximum of each tracer [conc]
  type(ocean_grid_type),          intent(in)    :: G      !< The ocean's grid structure
  type(MOM_generic_tracer_CS),    pointer       :: CS     !< Pointer to the control structure for this module.
  character(len=*), dimension(:), intent(out)   :: names  !< The names of the stocks calculated.
  character(len=*), dimension(:), intent(out)   :: units  !< The units of the stocks calculated.
  real, dimension(:),   optional, intent(out)   :: xgmin  !< The x-position of the global minimum in the
                                                          !! units of G%geoLonT, often [degrees_E] or [km] or [m]
  real, dimension(:),   optional, intent(out)   :: ygmin  !< The y-position of the global minimum in the
                                                          !! units of G%geoLatT, often [degrees_N] or [km] or [m]
  real, dimension(:),   optional, intent(out)   :: zgmin  !< The z-position of the global minimum [layer]
  real, dimension(:),   optional, intent(out)   :: xgmax  !< The x-position of the global maximum in the
                                                          !! units of G%geoLonT, often [degrees_E] or [km] or [m]
  real, dimension(:),   optional, intent(out)   :: ygmax  !< The y-position of the global maximum in the
                                                          !! units of G%geoLatT, often [degrees_N] or [km] or [m]
  real, dimension(:),   optional, intent(out)   :: zgmax  !< The z-position of the global maximum [layer]
  integer                                       :: MOM_generic_tracer_min_max !< Return value, the
                                                          !! number of tracers done here.

  MOM_generic_tracer_min_max = 0

end function MOM_generic_tracer_min_max

!> This subroutine calculates the surface state and sets coupler values for
!! those generic tracers that have flux exchange with atmosphere.
!!
!! This subroutine sets up the fields that the coupler needs to calculate the
!! CFC fluxes between the ocean and atmosphere.
subroutine MOM_generic_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),                 intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),               intent(in)    :: GV   !< The ocean's vertical grid structure
  type(surface),                         intent(inout) :: sfc_state !< A structure containing fields that
                                                               !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(MOM_generic_tracer_CS),           pointer       :: CS   !< Pointer to the control structure for this module.

end subroutine MOM_generic_tracer_surface_state

!ALL PE subroutine on Ocean!  Due to otpm design the fluxes should be initialized like this on ALL PE's!
subroutine MOM_generic_flux_init(verbosity)
  integer, optional, intent(in) :: verbosity  !< A 0-9 integer indicating a level of verbosity.

end subroutine MOM_generic_flux_init

subroutine MOM_generic_tracer_fluxes_accumulate(flux_tmp, weight)
  type(forcing), intent(in)    :: flux_tmp  !< A structure containing pointers to
                                            !! thermodynamic and tracer forcing fields.
  real,          intent(in)    :: weight    !< A weight for accumulating this flux [nondim]

end subroutine MOM_generic_tracer_fluxes_accumulate

!> Copy the requested tracer into an array.
subroutine MOM_generic_tracer_get(name,member,array, CS)
  character(len=*),         intent(in)  :: name   !< Name of requested tracer.
  character(len=*),         intent(in)  :: member !< The tracer element to return.
  real, dimension(:,:,:),   intent(out) :: array  !< Array filled by this routine, in arbitrary units [A]
  type(MOM_generic_tracer_CS), pointer  :: CS     !< Pointer to the control structure for this module.

  ! Local variables
  real, dimension(:,:,:),   pointer :: array_ptr  ! The tracer in the generic tracer structures, in
                                                  ! arbitrary units [A]
  character(len=128), parameter :: sub_name = 'MOM_generic_tracer_get'

end subroutine MOM_generic_tracer_get

!> This subroutine deallocates the memory owned by this module.
subroutine end_MOM_generic_tracer(CS)
  type(MOM_generic_tracer_CS), pointer :: CS   !< Pointer to the control structure for this module.

end subroutine end_MOM_generic_tracer

!----------------------------------------------------------------
! <CONTACT EMAIL="Niki.Zadeh@noaa.gov"> Niki Zadeh
! </CONTACT>
!
! <REVIEWER EMAIL="William.Cooke@noaa.gov"> William Cooke
! </REVIEWER>
!
! <OVERVIEW>
!  This module drives the generic version of tracers TOPAZ and CFC
! </OVERVIEW>
!----------------------------------------------------------------

end module MOM_generic_tracer
