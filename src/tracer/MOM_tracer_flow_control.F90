!> Orchestrates the registration and calling of tracer packages
module MOM_tracer_flow_control

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type, close_param_file
use MOM_forcing_type, only : forcing, optics_type
use MOM_get_input, only : Get_MOM_input
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : MOM_restart_CS
use MOM_sponge, only : sponge_CS
use MOM_ALE_sponge, only : ALE_sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface, thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
#include <MOM_memory.h>

! Add references to other user-provide tracer modules here.
use USER_tracer_example, only : tracer_column_physics, USER_initialize_tracer, USER_tracer_stock
use USER_tracer_example, only : USER_register_tracer_example, USER_tracer_surface_state
use USER_tracer_example, only : USER_tracer_example_end, USER_tracer_example_CS
use DOME_tracer, only : register_DOME_tracer, initialize_DOME_tracer
use DOME_tracer, only : DOME_tracer_column_physics, DOME_tracer_surface_state
use DOME_tracer, only : DOME_tracer_end, DOME_tracer_CS
use ISOMIP_tracer, only : register_ISOMIP_tracer, initialize_ISOMIP_tracer
use ISOMIP_tracer, only : ISOMIP_tracer_column_physics, ISOMIP_tracer_surface_state
use ISOMIP_tracer, only : ISOMIP_tracer_end, ISOMIP_tracer_CS
use RGC_tracer, only : register_RGC_tracer, initialize_RGC_tracer
use RGC_tracer, only : RGC_tracer_column_physics
use RGC_tracer, only : RGC_tracer_end, RGC_tracer_CS
use ideal_age_example, only : register_ideal_age_tracer, initialize_ideal_age_tracer
use ideal_age_example, only : ideal_age_tracer_column_physics, ideal_age_tracer_surface_state
use ideal_age_example, only : ideal_age_stock, ideal_age_example_end, ideal_age_tracer_CS
use regional_dyes, only : register_dye_tracer, initialize_dye_tracer
use regional_dyes, only : dye_tracer_column_physics, dye_tracer_surface_state
use regional_dyes, only : dye_stock, regional_dyes_end, dye_tracer_CS
use MOM_OCMIP2_CFC, only : register_OCMIP2_CFC, initialize_OCMIP2_CFC, flux_init_OCMIP2_CFC
use MOM_OCMIP2_CFC, only : OCMIP2_CFC_column_physics, OCMIP2_CFC_surface_state
use MOM_OCMIP2_CFC, only : OCMIP2_CFC_stock, OCMIP2_CFC_end, OCMIP2_CFC_CS
use oil_tracer, only : register_oil_tracer, initialize_oil_tracer
use oil_tracer, only : oil_tracer_column_physics, oil_tracer_surface_state
use oil_tracer, only : oil_stock, oil_tracer_end, oil_tracer_CS
use advection_test_tracer, only : register_advection_test_tracer, initialize_advection_test_tracer
use advection_test_tracer, only : advection_test_tracer_column_physics, advection_test_tracer_surface_state
use advection_test_tracer, only : advection_test_stock, advection_test_tracer_end, advection_test_tracer_CS
use dyed_obc_tracer, only : register_dyed_obc_tracer, initialize_dyed_obc_tracer
use dyed_obc_tracer, only : dyed_obc_tracer_column_physics
use dyed_obc_tracer, only : dyed_obc_tracer_end, dyed_obc_tracer_CS
#ifdef _USE_GENERIC_TRACER
use MOM_generic_tracer, only : register_MOM_generic_tracer, initialize_MOM_generic_tracer
use MOM_generic_tracer, only : MOM_generic_tracer_column_physics, MOM_generic_tracer_surface_state
use MOM_generic_tracer, only : end_MOM_generic_tracer, MOM_generic_tracer_get, MOM_generic_flux_init
use MOM_generic_tracer, only : MOM_generic_tracer_stock, MOM_generic_tracer_min_max, MOM_generic_tracer_CS
#endif
use pseudo_salt_tracer, only : register_pseudo_salt_tracer, initialize_pseudo_salt_tracer
use pseudo_salt_tracer, only : pseudo_salt_tracer_column_physics, pseudo_salt_tracer_surface_state
use pseudo_salt_tracer, only : pseudo_salt_stock, pseudo_salt_tracer_end, pseudo_salt_tracer_CS
use boundary_impulse_tracer, only : register_boundary_impulse_tracer, initialize_boundary_impulse_tracer
use boundary_impulse_tracer, only : boundary_impulse_tracer_column_physics, boundary_impulse_tracer_surface_state
use boundary_impulse_tracer, only : boundary_impulse_stock, boundary_impulse_tracer_end
use boundary_impulse_tracer, only : boundary_impulse_tracer_CS

implicit none ; private

public call_tracer_register, tracer_flow_control_init, call_tracer_set_forcing
public call_tracer_column_fns, call_tracer_surface_state, call_tracer_stocks
public call_tracer_flux_init, get_chl_from_model, tracer_flow_control_end

!> The control structure for orchestrating the calling of tracer packages
type, public :: tracer_flow_control_CS ; private
  logical :: use_USER_tracer_example = .false.     !< If true, use the USER_tracer_example package
  logical :: use_DOME_tracer = .false.             !< If true, use the DOME_tracer package
  logical :: use_ISOMIP_tracer = .false.           !< If true, use the ISOMPE_tracer package
  logical :: use_RGC_tracer =.false.               !< If true, use the RGC_tracer package
  logical :: use_ideal_age = .false.               !< If true, use the ideal age tracer package
  logical :: use_regional_dyes = .false.           !< If true, use the regional dyes tracer package
  logical :: use_oil = .false.                     !< If true, use the oil tracer package
  logical :: use_advection_test_tracer = .false.   !< If true, use the advection_test_tracer package
  logical :: use_OCMIP2_CFC = .false.              !< If true, use the OCMIP2_CFC tracer package
  logical :: use_MOM_generic_tracer = .false.      !< If true, use the MOM_generic_tracer packages
  logical :: use_pseudo_salt_tracer = .false.      !< If true, use the psuedo_salt tracer  package
  logical :: use_boundary_impulse_tracer = .false. !< If true, use the boundary impulse tracer package
  logical :: use_dyed_obc_tracer = .false.         !< If true, use the dyed OBC tracer package
  !>@{ Pointers to the control strucures for the tracer packages
  type(USER_tracer_example_CS), pointer :: USER_tracer_example_CSp => NULL()
  type(DOME_tracer_CS), pointer :: DOME_tracer_CSp => NULL()
  type(ISOMIP_tracer_CS), pointer :: ISOMIP_tracer_CSp => NULL()
  type(RGC_tracer_CS), pointer :: RGC_tracer_CSp => NULL()
  type(ideal_age_tracer_CS), pointer :: ideal_age_tracer_CSp => NULL()
  type(dye_tracer_CS), pointer :: dye_tracer_CSp => NULL()
  type(oil_tracer_CS), pointer :: oil_tracer_CSp => NULL()
  type(advection_test_tracer_CS), pointer :: advection_test_tracer_CSp => NULL()
  type(OCMIP2_CFC_CS), pointer :: OCMIP2_CFC_CSp => NULL()
#ifdef _USE_GENERIC_TRACER
  type(MOM_generic_tracer_CS), pointer :: MOM_generic_tracer_CSp => NULL()
#endif
  type(pseudo_salt_tracer_CS), pointer :: pseudo_salt_tracer_CSp => NULL()
  type(boundary_impulse_tracer_CS), pointer :: boundary_impulse_tracer_CSp => NULL()
  type(dyed_obc_tracer_CS), pointer :: dyed_obc_tracer_CSp => NULL()
  !!@}
end type tracer_flow_control_CS

contains

!> This subroutine carries out a series of calls to initialize the air-sea
!! tracer fluxes, but it does not record the generated indicies, and it may
!! be called _before_ the ocean model has been initialized and may be called
!! on non-ocean PEs.  It is not necessary to call this routine for ocean-only
!! runs, because the same calls are made again inside of the routines called by
!! call_tracer_register
subroutine call_tracer_flux_init(verbosity)
  integer, optional, intent(in) :: verbosity !< A 0-9 integer indicating a level of verbosity.

  type(param_file_type) :: param_file ! A structure to parse for run-time parameters
  character(len=40)  :: mdl = "call_tracer_flux_init"  ! This module's name.
  logical :: use_OCMIP_CFCs, use_MOM_generic_tracer

  ! Determine which tracer routines with tracer fluxes are to be called.  Note
  ! that not every tracer package is required to have a flux_init call.
  call get_MOM_Input(param_file, check_params=.false.)

  call get_param(param_file, mdl, "USE_OCMIP2_CFC", use_OCMIP_CFCs, &
                 default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "USE_generic_tracer", use_MOM_generic_tracer,&
                 default=.false., do_not_log=.true.)
  call close_param_file(param_file, quiet_close=.true.)

  if (use_OCMIP_CFCs) call flux_init_OCMIP2_CFC(verbosity=verbosity)
  if (use_MOM_generic_tracer) then
#ifdef _USE_GENERIC_TRACER
    call MOM_generic_flux_init(verbosity=verbosity)
#else
    call MOM_error(FATAL, &
       "call_tracer_flux_init: use_MOM_generic_tracer=.true. but MOM6 was "//&
       "not compiled with _USE_GENERIC_TRACER")
#endif
  endif

end subroutine call_tracer_flux_init

! The following 5 subroutines and associated definitions provide the machinery to register and call
! the subroutines that initialize tracers and apply vertical column processes to tracers.

!> This subroutine determines which tracer packages are to be used and does the calls to
!! register their tracers to be advected, diffused, and read from restarts.
subroutine call_tracer_register(HI, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),         intent(in) :: HI         !< A horizontal index type structure.
  type(verticalGrid_type),      intent(in) :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),        intent(in) :: US         !< A dimensional unit scaling type
  type(param_file_type),        intent(in) :: param_file !< A structure to parse for run-time
                                                         !! parameters.
  type(tracer_flow_control_CS), pointer    :: CS         !< A pointer that is set to point to the
                                                         !! control structure for this module.
  type(tracer_registry_type),   pointer    :: tr_Reg     !< A pointer that is set to point to the
                                                         !! control structure for the tracer
                                                         !! advection and diffusion module.
  type(MOM_restart_CS),         pointer    :: restart_CS !< A pointer to the restart control
                                                         !! structure.


  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_tracer_flow_control" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "call_tracer_register called with an associated "// &
                            "control structure.")
    return
  else ; allocate(CS) ; endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "USE_USER_TRACER_EXAMPLE", &
                                CS%use_USER_tracer_example, &
                 "If true, use the USER_tracer_example tracer package.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_DOME_TRACER", CS%use_DOME_tracer, &
                 "If true, use the DOME_tracer tracer package.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_ISOMIP_TRACER", CS%use_ISOMIP_tracer, &
                 "If true, use the ISOMIP_tracer tracer package.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_RGC_TRACER", CS%use_RGC_tracer, &
                 "If true, use the RGC_tracer tracer package.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_IDEAL_AGE_TRACER", CS%use_ideal_age, &
                 "If true, use the ideal_age_example tracer package.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_REGIONAL_DYES", CS%use_regional_dyes, &
                 "If true, use the regional_dyes tracer package.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_OIL_TRACER", CS%use_oil, &
                 "If true, use the oil_tracer tracer package.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_ADVECTION_TEST_TRACER", CS%use_advection_test_tracer, &
                 "If true, use the advection_test_tracer tracer package.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_OCMIP2_CFC", CS%use_OCMIP2_CFC, &
                 "If true, use the MOM_OCMIP2_CFC tracer package.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_generic_tracer", CS%use_MOM_generic_tracer, &
                 "If true and _USE_GENERIC_TRACER is defined as a "//&
                 "preprocessor macro, use the MOM_generic_tracer packages.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_PSEUDO_SALT_TRACER", CS%use_pseudo_salt_tracer, &
                 "If true, use the pseudo salt tracer, typically run as a diagnostic.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_BOUNDARY_IMPULSE_TRACER", CS%use_boundary_impulse_tracer, &
                 "If true, use the boundary impulse tracer.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_DYED_OBC_TRACER", CS%use_dyed_obc_tracer, &
                 "If true, use the dyed_obc_tracer tracer package.", &
                 default=.false.)

#ifndef _USE_GENERIC_TRACER
  if (CS%use_MOM_generic_tracer) call MOM_error(FATAL, &
       "call_tracer_register: use_MOM_generic_tracer=.true. but MOM6 was "//&
       "not compiled with _USE_GENERIC_TRACER")
#endif

!    Add other user-provided calls to register tracers for restarting here. Each
!  tracer package registration call returns a logical false if it cannot be run
!  for some reason.  This then overrides the run-time selection from above.
  if (CS%use_USER_tracer_example) CS%use_USER_tracer_example = &
    USER_register_tracer_example(HI, GV, param_file, CS%USER_tracer_example_CSp, &
                                 tr_Reg, restart_CS)
  if (CS%use_DOME_tracer) CS%use_DOME_tracer = &
    register_DOME_tracer(HI, GV, param_file, CS%DOME_tracer_CSp, &
                         tr_Reg, restart_CS)
  if (CS%use_ISOMIP_tracer) CS%use_ISOMIP_tracer = &
    register_ISOMIP_tracer(HI, GV, param_file, CS%ISOMIP_tracer_CSp, &
                           tr_Reg, restart_CS)
  if (CS%use_RGC_tracer) CS%use_RGC_tracer = &
    register_RGC_tracer(HI, GV, param_file, CS%RGC_tracer_CSp, &
                           tr_Reg, restart_CS)
  if (CS%use_ideal_age) CS%use_ideal_age = &
    register_ideal_age_tracer(HI, GV, param_file,  CS%ideal_age_tracer_CSp, &
                              tr_Reg, restart_CS)
  if (CS%use_regional_dyes) CS%use_regional_dyes = &
    register_dye_tracer(HI, GV, US, param_file,  CS%dye_tracer_CSp, &
                        tr_Reg, restart_CS)
  if (CS%use_oil) CS%use_oil = &
    register_oil_tracer(HI, GV, US, param_file,  CS%oil_tracer_CSp, &
                        tr_Reg, restart_CS)
  if (CS%use_advection_test_tracer) CS%use_advection_test_tracer = &
    register_advection_test_tracer(HI, GV, param_file, CS%advection_test_tracer_CSp, &
                                   tr_Reg, restart_CS)
  if (CS%use_OCMIP2_CFC) CS%use_OCMIP2_CFC = &
    register_OCMIP2_CFC(HI, GV, param_file,  CS%OCMIP2_CFC_CSp, &
                        tr_Reg, restart_CS)
#ifdef _USE_GENERIC_TRACER
  if (CS%use_MOM_generic_tracer) CS%use_MOM_generic_tracer = &
    register_MOM_generic_tracer(HI, GV, param_file,  CS%MOM_generic_tracer_CSp, &
                                tr_Reg, restart_CS)
#endif
  if (CS%use_pseudo_salt_tracer) CS%use_pseudo_salt_tracer = &
    register_pseudo_salt_tracer(HI, GV, param_file,  CS%pseudo_salt_tracer_CSp, &
                                tr_Reg, restart_CS)
  if (CS%use_boundary_impulse_tracer) CS%use_boundary_impulse_tracer = &
    register_boundary_impulse_tracer(HI, GV, param_file,  CS%boundary_impulse_tracer_CSp, &
                                     tr_Reg, restart_CS)
  if (CS%use_dyed_obc_tracer) CS%use_dyed_obc_tracer = &
    register_dyed_obc_tracer(HI, GV, param_file, CS%dyed_obc_tracer_CSp, &
                             tr_Reg, restart_CS)


end subroutine call_tracer_register

!> This subroutine calls all registered tracer initialization
!! subroutines.
subroutine tracer_flow_control_init(restart, day, G, GV, US, h, param_file, diag, OBC, &
                                    CS, sponge_CSp, ALE_sponge_CSp, tv)
  logical,                               intent(in)    :: restart !< 1 if the fields have already
                                                                  !! been read from a restart file.
  type(time_type), target,               intent(in)    :: day     !< Time of the start of the run.
  type(ocean_grid_type),                 intent(inout) :: G       !< The ocean's grid structure.
  type(verticalGrid_type),               intent(in)    :: GV      !< The ocean's vertical grid
                                                                  !! structure.
  type(unit_scale_type),                 intent(in)    :: US      !< A dimensional unit scaling type
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h       !< Layer thicknesses [H ~> m or kg m-2]
  type(param_file_type),                 intent(in)    :: param_file !< A structure to parse for
                                                                  !! run-time parameters
  type(diag_ctrl), target,               intent(in)    :: diag    !< A structure that is used to
                                                                  !! regulate diagnostic output.
  type(ocean_OBC_type),                  pointer       :: OBC     !< This open boundary condition
                                                                  !! type specifies whether, where,
                                                                  !! and what open boundary
                                                                  !! conditions are used.
  type(tracer_flow_control_CS),          pointer       :: CS      !< The control structure returned
                                                                  !! by a previous call to
                                                                  !! call_tracer_register.
  type(sponge_CS),                       pointer       :: sponge_CSp     !< A pointer to the control
                                               !! structure for the sponges, if they are in use.
                                               !! Otherwise this may be unassociated.
  type(ALE_sponge_CS),                   pointer       :: ALE_sponge_CSp !< A pointer to the control
                                               !! structure for the ALE sponges, if they are in use.
                                               !! Otherwise this may be unassociated.
  type(thermo_var_ptrs),                 intent(in)    :: tv      !< A structure pointing to various
                                                                  !! thermodynamic variables

  if (.not. associated(CS)) call MOM_error(FATAL, "tracer_flow_control_init: "// &
         "Module must be initialized via call_tracer_register before it is used.")

!  Add other user-provided calls here.
  if (CS%use_USER_tracer_example) &
    call USER_initialize_tracer(restart, day, G, GV, h, diag, OBC, CS%USER_tracer_example_CSp, &
                                sponge_CSp)
  if (CS%use_DOME_tracer) &
    call initialize_DOME_tracer(restart, day, G, GV, US, h, diag, OBC, CS%DOME_tracer_CSp, &
                                sponge_CSp, param_file)
  if (CS%use_ISOMIP_tracer) &
    call initialize_ISOMIP_tracer(restart, day, G, GV, h, diag, OBC, CS%ISOMIP_tracer_CSp, &
                                ALE_sponge_CSp)
  if (CS%use_RGC_tracer) &
    call initialize_RGC_tracer(restart, day, G, GV, h, diag, OBC, &
                  CS%RGC_tracer_CSp, sponge_CSp, ALE_sponge_CSp)
  if (CS%use_ideal_age) &
    call initialize_ideal_age_tracer(restart, day, G, GV, US, h, diag, OBC, CS%ideal_age_tracer_CSp, &
                                     sponge_CSp)
  if (CS%use_regional_dyes) &
    call initialize_dye_tracer(restart, day, G, GV, h, diag, OBC, CS%dye_tracer_CSp, &
                                     sponge_CSp)
  if (CS%use_oil) &
    call initialize_oil_tracer(restart, day, G, GV, US, h, diag, OBC, CS%oil_tracer_CSp, &
                                     sponge_CSp)
  if (CS%use_advection_test_tracer) &
    call initialize_advection_test_tracer(restart, day, G, GV, h, diag, OBC, CS%advection_test_tracer_CSp, &
                                sponge_CSp)
  if (CS%use_OCMIP2_CFC) &
    call initialize_OCMIP2_CFC(restart, day, G, GV, US, h, diag, OBC, CS%OCMIP2_CFC_CSp, &
                                sponge_CSp)
#ifdef _USE_GENERIC_TRACER
  if (CS%use_MOM_generic_tracer) &
    call initialize_MOM_generic_tracer(restart, day, G, GV, US, h, param_file, diag, OBC, &
        CS%MOM_generic_tracer_CSp, sponge_CSp, ALE_sponge_CSp)
#endif
  if (CS%use_pseudo_salt_tracer) &
    call initialize_pseudo_salt_tracer(restart, day, G, GV, h, diag, OBC, CS%pseudo_salt_tracer_CSp, &
                                sponge_CSp, tv)
  if (CS%use_boundary_impulse_tracer) &
    call initialize_boundary_impulse_tracer(restart, day, G, GV, h, diag, OBC, CS%boundary_impulse_tracer_CSp, &
                                sponge_CSp, tv)
  if (CS%use_dyed_obc_tracer) &
    call initialize_dyed_obc_tracer(restart, day, G, GV, h, diag, OBC, CS%dyed_obc_tracer_CSp)

end subroutine tracer_flow_control_init

!> This subroutine extracts the chlorophyll concentrations from the model state, if possible
subroutine get_chl_from_model(Chl_array, G, CS)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), &
                                intent(out) :: Chl_array !< The array in which to store the model's
                                                         !! Chlorophyll-A concentrations in mg m-3.
  type(ocean_grid_type),        intent(in)  :: G         !< The ocean's grid structure.
  type(tracer_flow_control_CS), pointer     :: CS        !< The control structure returned by a
                                                         !! previous call to call_tracer_register.

#ifdef _USE_GENERIC_TRACER
  if (CS%use_MOM_generic_tracer) then
    call MOM_generic_tracer_get('chl','field',Chl_array, CS%MOM_generic_tracer_CSp)
  else
    call MOM_error(FATAL, "get_chl_from_model was called in a configuration "// &
             "that is unable to provide a sensible model-based value.\n"// &
             "CS%use_MOM_generic_tracer is false and no other viable options are on.")
  endif
#else
  call MOM_error(FATAL, "get_chl_from_model was called in a configuration "// &
           "that is unable to provide a sensible model-based value.\n"// &
           "_USE_GENERIC_TRACER is undefined and no other options "//&
           "are currently viable.")
#endif

end subroutine get_chl_from_model

!> This subroutine calls the individual tracer modules' subroutines to
!! specify or read quantities related to their surface forcing.
subroutine call_tracer_set_forcing(state, fluxes, day_start, day_interval, G, CS)

  type(surface),                intent(inout) :: state     !< A structure containing fields that
                                                           !! describe the surface state of the
                                                           !! ocean.
  type(forcing),                intent(inout) :: fluxes    !< A structure containing pointers to any
                                                           !! possible forcing fields. Unused fields
                                                           !! have NULL ptrs.
  type(time_type),              intent(in)    :: day_start !< Start time of the fluxes.
  type(time_type),              intent(in)    :: day_interval !< Length of time over which these
                                                           !! fluxes will be applied.
  type(ocean_grid_type),        intent(in)    :: G         !< The ocean's grid structure.
  type(tracer_flow_control_CS), pointer       :: CS        !< The control structure returned by a
                                                           !! previous call to call_tracer_register.

  if (.not. associated(CS)) call MOM_error(FATAL, "call_tracer_set_forcing"// &
         "Module must be initialized via call_tracer_register before it is used.")
!  if (CS%use_ideal_age) &
!    call ideal_age_tracer_set_forcing(state, fluxes, day_start, day_interval, &
!                                      G, CS%ideal_age_tracer_CSp)

end subroutine call_tracer_set_forcing

!> This subroutine calls all registered tracer column physics subroutines.
subroutine call_tracer_column_fns(h_old, h_new, ea, eb, fluxes, Hml, dt, G, GV, US, tv, optics, CS, &
                                  debug, evap_CFL_limit, minimum_forcing_depth)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: h_old  !< Layer thickness before entrainment
                                                              !! [H ~> m or kg m-2].
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: h_new  !< Layer thickness after entrainment
                                                              !! [H ~> m or kg m-2].
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: ea     !< an array to which the amount of
                                          !! fluid entrained from the layer above during this call
                                          !! will be added [H ~> m or kg m-2].
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: eb     !< an array to which the amount of
                                          !! fluid entrained from the layer below during this call
                                          !! will be added [H ~> m or kg m-2].
  type(forcing),                         intent(in) :: fluxes !< A structure containing pointers to
                                                              !! any possible forcing fields.
                                                              !! Unused fields have NULL ptrs.
  real, dimension(NIMEM_,NJMEM_),        intent(in) :: Hml    !< Mixed layer depth [H ~> m or kg m-2]
  real,                                  intent(in) :: dt     !< The amount of time covered by this
                                                              !! call [T ~> s]
  type(ocean_grid_type),                 intent(in) :: G      !< The ocean's grid structure.
  type(verticalGrid_type),               intent(in) :: GV     !< The ocean's vertical grid
                                                              !! structure.
  type(unit_scale_type),                 intent(in) :: US     !< A dimensional unit scaling type
  type(thermo_var_ptrs),                 intent(in) :: tv     !< A structure pointing to various
                                                              !! thermodynamic variables.
  type(optics_type),                     pointer    :: optics !< The structure containing optical
                                                              !! properties.
  type(tracer_flow_control_CS),          pointer    :: CS     !< The control structure returned by
                                                              !! a previous call to
                                                              !! call_tracer_register.
  logical,                               intent(in) :: debug  !< If true calculate checksums
  real,                        optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of
                                                              !! the water that can be fluxed out
                                                              !! of the top layer in a timestep [nondim]
  real,                        optional, intent(in) :: minimum_forcing_depth !< The smallest depth over
                                                              !! which fluxes can be applied [H ~> m or kg m-2]

  if (.not. associated(CS)) call MOM_error(FATAL, "call_tracer_column_fns: "// &
         "Module must be initialized via call_tracer_register before it is used.")

  ! Use the applyTracerBoundaryFluxesInOut to handle surface fluxes
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    ! Add calls to tracer column functions here.
    if (CS%use_USER_tracer_example) &
      call tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                 G, GV, US, CS%USER_tracer_example_CSp)
    if (CS%use_DOME_tracer) &
      call DOME_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                      G, GV, US, CS%DOME_tracer_CSp, &
                                      evap_CFL_limit=evap_CFL_limit, &
                                      minimum_forcing_depth=minimum_forcing_depth)
    if (CS%use_ISOMIP_tracer) &
      call ISOMIP_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                        G, GV, US, CS%ISOMIP_tracer_CSp, &
                                        evap_CFL_limit=evap_CFL_limit, &
                                        minimum_forcing_depth=minimum_forcing_depth)
    if (CS%use_RGC_tracer) &
      call RGC_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                        G, GV, US, CS%RGC_tracer_CSp, &
                                        evap_CFL_limit=evap_CFL_limit, &
                                        minimum_forcing_depth=minimum_forcing_depth)
    if (CS%use_ideal_age) &
      call ideal_age_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                           G, GV, US, CS%ideal_age_tracer_CSp, &
                                           evap_CFL_limit=evap_CFL_limit, &
                                           minimum_forcing_depth=minimum_forcing_depth)
    if (CS%use_regional_dyes) &
      call dye_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                     G, GV, US, CS%dye_tracer_CSp, &
                                     evap_CFL_limit=evap_CFL_limit, &
                                     minimum_forcing_depth=minimum_forcing_depth)
    if (CS%use_oil) &
      call oil_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                     G, GV, US, CS%oil_tracer_CSp, tv, &
                                     evap_CFL_limit=evap_CFL_limit, &
                                     minimum_forcing_depth=minimum_forcing_depth)

    if (CS%use_advection_test_tracer) &
      call advection_test_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                                G, GV, US, CS%advection_test_tracer_CSp, &
                                                evap_CFL_limit=evap_CFL_limit, &
                                                minimum_forcing_depth=minimum_forcing_depth)
    if (CS%use_OCMIP2_CFC) &
      call OCMIP2_CFC_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                     G, GV, US, CS%OCMIP2_CFC_CSp, &
                                     evap_CFL_limit=evap_CFL_limit, &
                                     minimum_forcing_depth=minimum_forcing_depth)
#ifdef _USE_GENERIC_TRACER
    if (CS%use_MOM_generic_tracer) &
      call MOM_generic_tracer_column_physics(h_old, h_new, ea, eb, fluxes, Hml, US%T_to_s*dt, &
                                             G, GV, CS%MOM_generic_tracer_CSp, tv, optics, &
                                             evap_CFL_limit=evap_CFL_limit, &
                                             minimum_forcing_depth=minimum_forcing_depth)
#endif
    if (CS%use_pseudo_salt_tracer) &
      call pseudo_salt_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                     G, GV, US, CS%pseudo_salt_tracer_CSp, tv, debug, &
                                     evap_CFL_limit=evap_CFL_limit, &
                                     minimum_forcing_depth=minimum_forcing_depth)
    if (CS%use_boundary_impulse_tracer) &
      call boundary_impulse_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                     G, GV, US, CS%boundary_impulse_tracer_CSp, tv, debug, &
                                     evap_CFL_limit=evap_CFL_limit, &
                                     minimum_forcing_depth=minimum_forcing_depth)
    if (CS%use_dyed_obc_tracer) &
      call dyed_obc_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                      G, GV, US, CS%dyed_obc_tracer_CSp, &
                                      evap_CFL_limit=evap_CFL_limit, &
                                      minimum_forcing_depth=minimum_forcing_depth)


  else ! Apply tracer surface fluxes using ea on the first layer
    if (CS%use_USER_tracer_example) &
      call tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                 G, GV, US, CS%USER_tracer_example_CSp)
    if (CS%use_DOME_tracer) &
      call DOME_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                      G, GV, US, CS%DOME_tracer_CSp)
    if (CS%use_ISOMIP_tracer) &
      call ISOMIP_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                      G, GV, US, CS%ISOMIP_tracer_CSp)
    if (CS%use_RGC_tracer) &
      call RGC_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                      G, GV, US, CS%RGC_tracer_CSp)
    if (CS%use_ideal_age) &
      call ideal_age_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                           G, GV, US, CS%ideal_age_tracer_CSp)
    if (CS%use_regional_dyes) &
      call dye_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                           G, GV, US, CS%dye_tracer_CSp)
    if (CS%use_oil) &
      call oil_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                     G, GV, US, CS%oil_tracer_CSp, tv)
    if (CS%use_advection_test_tracer) &
      call advection_test_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                      G, GV, US, CS%advection_test_tracer_CSp)
    if (CS%use_OCMIP2_CFC) &
      call OCMIP2_CFC_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                     G, GV, US, CS%OCMIP2_CFC_CSp)
#ifdef _USE_GENERIC_TRACER
    if (CS%use_MOM_generic_tracer) &
      call MOM_generic_tracer_column_physics(h_old, h_new, ea, eb, fluxes, Hml, US%T_to_s*dt, &
                                     G, GV, CS%MOM_generic_tracer_CSp, tv, optics)
#endif
    if (CS%use_pseudo_salt_tracer) &
      call pseudo_salt_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                     G, GV, US, CS%pseudo_salt_tracer_CSp, tv, debug)
    if (CS%use_boundary_impulse_tracer) &
      call boundary_impulse_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                     G, GV, US, CS%boundary_impulse_tracer_CSp, tv, debug)
    if (CS%use_dyed_obc_tracer) &
      call dyed_obc_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, &
                                      G, GV, US, CS%dyed_obc_tracer_CSp)

  endif


end subroutine call_tracer_column_fns

!> This subroutine calls all registered tracer packages to enable them to
!! add to the surface state returned to the coupler. These routines are optional.
subroutine call_tracer_stocks(h, stock_values, G, GV, CS, stock_names, stock_units, &
                              num_stocks, stock_index, got_min_max, global_min, global_max, &
                              xgmin, ygmin, zgmin, xgmax, ygmax, zgmax)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),    &
                                  intent(in)  :: h           !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(:),             intent(out) :: stock_values !< The integrated amounts of a tracer
                             !! on the current PE, usually in kg x concentration [kg conc].
  type(ocean_grid_type),          intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type),        intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(tracer_flow_control_CS),   pointer     :: CS          !< The control structure returned by a
                                                             !! previous call to
                                                             !! call_tracer_register.
  character(len=*), dimension(:), &
                        optional, intent(out) :: stock_names !< Diagnostic names to use for each stock.
  character(len=*), dimension(:), &
                        optional, intent(out) :: stock_units !< Units to use in the metadata for each stock.
  integer,              optional, intent(out) :: num_stocks  !< The number of tracer stocks being returned.
  integer,              optional, intent(in)  :: stock_index !< The integer stock index from
                             !! stocks_constants_mod of the stock to be returned.  If this is
                             !! present and greater than 0, only a single stock can be returned.
  logical, dimension(:), &
                      optional, intent(inout) :: got_min_max !< Indicates whether the global min and
                                                             !! max are found for each tracer
  real, dimension(:), optional, intent(out)   :: global_min  !< The global minimum of each tracer
  real, dimension(:), optional, intent(out)   :: global_max  !< The global maximum of each tracer
  real, dimension(:), optional, intent(out)   :: xgmin       !< The x-position of the global minimum
  real, dimension(:), optional, intent(out)   :: ygmin       !< The y-position of the global minimum
  real, dimension(:), optional, intent(out)   :: zgmin       !< The z-position of the global minimum
  real, dimension(:), optional, intent(out)   :: xgmax       !< The x-position of the global maximum
  real, dimension(:), optional, intent(out)   :: ygmax       !< The y-position of the global maximum
  real, dimension(:), optional, intent(out)   :: zgmax       !< The z-position of the global maximum

  ! Local variables
  character(len=200), dimension(MAX_FIELDS_) :: names, units
  character(len=200) :: set_pkg_name
  real, dimension(MAX_FIELDS_) :: values
  integer :: max_ns, ns_tot, ns, index, pkg, max_pkgs, nn

  if (.not. associated(CS)) call MOM_error(FATAL, "call_tracer_stocks: "// &
       "Module must be initialized via call_tracer_register before it is used.")

  index = -1 ; if (present(stock_index)) index = stock_index
  ns_tot = 0
  max_ns = size(stock_values)
  if (present(stock_names)) max_ns = min(max_ns,size(stock_names))
  if (present(stock_units)) max_ns = min(max_ns,size(stock_units))

!  Add other user-provided calls here.
  if (CS%use_USER_tracer_example) then
    ns = USER_tracer_stock(h, values, G, GV, CS%USER_tracer_example_CSp, &
                           names, units, stock_index)
    call store_stocks("tracer_example", ns, names, units, values, index, stock_values, &
                       set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  endif
! if (CS%use_DOME_tracer) then
!   ns = DOME_tracer_stock(h, values, G, GV, CS%DOME_tracer_CSp, &
!                          names, units, stock_index)
!   call store_stocks("DOME_tracer", ns, names, units, values, index, stock_values, &
!                      set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
! endif
  if (CS%use_ideal_age) then
    ns = ideal_age_stock(h, values, G, GV, CS%ideal_age_tracer_CSp, &
                         names, units, stock_index)
    call store_stocks("ideal_age_example", ns, names, units, values, index, &
           stock_values, set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  endif
  if (CS%use_regional_dyes) then
    ns = dye_stock(h, values, G, GV, CS%dye_tracer_CSp, &
                         names, units, stock_index)
    call store_stocks("regional_dyes", ns, names, units, values, index, &
           stock_values, set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  endif
  if (CS%use_oil) then
    ns = oil_stock(h, values, G, GV, CS%oil_tracer_CSp, &
                         names, units, stock_index)
    call store_stocks("oil_tracer", ns, names, units, values, index, &
           stock_values, set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  endif
  if (CS%use_OCMIP2_CFC) then
    ns = OCMIP2_CFC_stock(h, values, G, GV, CS%OCMIP2_CFC_CSp, names, units, stock_index)
    call store_stocks("MOM_OCMIP2_CFC", ns, names, units, values, index, stock_values, &
                       set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  endif

  if (CS%use_advection_test_tracer) then
    ns = advection_test_stock( h, values, G, GV, CS%advection_test_tracer_CSp, &
                         names, units, stock_index )
    call store_stocks("advection_test_tracer", ns, names, units, values, index, &
           stock_values, set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  endif

#ifdef _USE_GENERIC_TRACER
  if (CS%use_MOM_generic_tracer) then
    ns = MOM_generic_tracer_stock(h, values, G, GV, CS%MOM_generic_tracer_CSp, &
                                   names, units, stock_index)
    call store_stocks("MOM_generic_tracer", ns, names, units, values, index, stock_values, &
                       set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
    nn=ns_tot-ns+1
    nn=MOM_generic_tracer_min_max(nn, got_min_max, global_min,  global_max, &
                                  xgmin, ygmin, zgmin, xgmax, ygmax, zgmax ,&
                                  G, CS%MOM_generic_tracer_CSp,names, units)

  endif
#endif
  if (CS%use_pseudo_salt_tracer) then
    ns = pseudo_salt_stock(h, values, G, GV, CS%pseudo_salt_tracer_CSp, &
                         names, units, stock_index)
    call store_stocks("pseudo_salt_tracer", ns, names, units, values, index, &
           stock_values, set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  endif

  if (CS%use_boundary_impulse_tracer) then
    ns = boundary_impulse_stock(h, values, G, GV, CS%boundary_impulse_tracer_CSp, &
                         names, units, stock_index)
    call store_stocks("boundary_impulse_tracer", ns, names, units, values, index, &
           stock_values, set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  endif

  if (ns_tot == 0) stock_values(1) = 0.0

  if (present(num_stocks)) num_stocks = ns_tot

end subroutine call_tracer_stocks

!> This routine stores the stocks and does error handling for call_tracer_stocks.
subroutine store_stocks(pkg_name, ns, names, units, values, index, stock_values, &
                        set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
  character(len=*),   intent(in)    :: pkg_name !< The tracer package name
  integer,            intent(in)    :: ns      !< The number of stocks associated with this tracer package
  character(len=*), dimension(:), &
                      intent(in)    :: names   !< Diagnostic names to use for each stock.
  character(len=*), dimension(:), &
                      intent(in)    :: units   !< Units to use in the metadata for each stock.
  real, dimension(:), intent(in)    :: values  !< The values of the tracer stocks
  integer,            intent(in)    :: index   !< The integer stock index from
                             !! stocks_constants_mod of the stock to be returned.  If this is
                             !! present and greater than 0, only a single stock can be returned.
  real, dimension(:), intent(inout) :: stock_values !< The master list of stock values
  character(len=*),   intent(inout) :: set_pkg_name !< The name of the last tracer package whose
                                               !! stocks were stored for a specific index.  This is
                                               !! used to trigger an error if there are redundant stocks.
  integer,            intent(in)    :: max_ns  !< The maximum size of the master stock list
  integer,            intent(inout) :: ns_tot  !< The total number of stocks in the master list
  character(len=*), dimension(:), &
            optional, intent(inout) :: stock_names !< Diagnostic names to use for each stock in the master list
  character(len=*), dimension(:), &
            optional, intent(inout) :: stock_units !< Units to use in the metadata for each stock in the master list

! This routine stores the stocks and does error handling for call_tracer_stocks.
  character(len=16) :: ind_text, ns_text, max_text
  integer :: n

  if ((index > 0) .and. (ns > 0)) then
    write(ind_text,'(i8)') index
    if (ns > 1) then
      call MOM_error(FATAL,"Tracer package "//trim(pkg_name)//&
          " is not permitted to return more than one value when queried"//&
          " for specific stock index "//trim(adjustl(ind_text))//".")
    elseif (ns+ns_tot > 1) then
      call MOM_error(FATAL,"Tracer packages "//trim(pkg_name)//" and "//&
          trim(set_pkg_name)//" both attempted to set values for"//&
          " specific stock index "//trim(adjustl(ind_text))//".")
    else
      set_pkg_name = pkg_name
    endif
  endif

  if (ns_tot+ns > max_ns) then
    write(ns_text,'(i8)') ns_tot+ns ; write(max_text,'(i8)') max_ns
    call MOM_error(FATAL,"Attempted to return more tracer stock values (at least "//&
      trim(adjustl(ns_text))//") than the size "//trim(adjustl(max_text))//&
      "of the smallest value, name, or units array.")
  endif

  do n=1,ns
    stock_values(ns_tot+n) = values(n)
    if (present(stock_names)) stock_names(ns_tot+n) = names(n)
    if (present(stock_units)) stock_units(ns_tot+n) = units(n)
  enddo
  ns_tot = ns_tot + ns

end subroutine store_stocks

!> This subroutine calls all registered tracer packages to enable them to
!! add to the surface state returned to the coupler. These routines are optional.
subroutine call_tracer_surface_state(state, h, G, CS)
  type(surface),                intent(inout) :: state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  real, dimension(NIMEM_,NJMEM_,NKMEM_), &
                                intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  type(ocean_grid_type),        intent(in)    :: G     !< The ocean's grid structure.
  type(tracer_flow_control_CS), pointer       :: CS    !< The control structure returned by a
                                                       !! previous call to call_tracer_register.

  if (.not. associated(CS)) call MOM_error(FATAL, "call_tracer_surface_state: "// &
         "Module must be initialized via call_tracer_register before it is used.")

!  Add other user-provided calls here.
  if (CS%use_USER_tracer_example) &
    call USER_tracer_surface_state(state, h, G, CS%USER_tracer_example_CSp)
  if (CS%use_DOME_tracer) &
    call DOME_tracer_surface_state(state, h, G, CS%DOME_tracer_CSp)
  if (CS%use_ISOMIP_tracer) &
    call ISOMIP_tracer_surface_state(state, h, G, CS%ISOMIP_tracer_CSp)
  if (CS%use_ideal_age) &
    call ideal_age_tracer_surface_state(state, h, G, CS%ideal_age_tracer_CSp)
  if (CS%use_regional_dyes) &
    call dye_tracer_surface_state(state, h, G, CS%dye_tracer_CSp)
  if (CS%use_oil) &
    call oil_tracer_surface_state(state, h, G, CS%oil_tracer_CSp)
  if (CS%use_advection_test_tracer) &
    call advection_test_tracer_surface_state(state, h, G, CS%advection_test_tracer_CSp)
  if (CS%use_OCMIP2_CFC) &
    call OCMIP2_CFC_surface_state(state, h, G, CS%OCMIP2_CFC_CSp)
#ifdef _USE_GENERIC_TRACER
  if (CS%use_MOM_generic_tracer) &
    call MOM_generic_tracer_surface_state(state, h, G, CS%MOM_generic_tracer_CSp)
#endif

end subroutine call_tracer_surface_state

subroutine tracer_flow_control_end(CS)
  type(tracer_flow_control_CS), pointer :: CS    !< The control structure returned by a
                                                 !! previous call to call_tracer_register.

  if (CS%use_USER_tracer_example) &
    call USER_tracer_example_end(CS%USER_tracer_example_CSp)
  if (CS%use_DOME_tracer) call DOME_tracer_end(CS%DOME_tracer_CSp)
  if (CS%use_ISOMIP_tracer) call ISOMIP_tracer_end(CS%ISOMIP_tracer_CSp)
  if (CS%use_RGC_tracer) call RGC_tracer_end(CS%RGC_tracer_CSp)
  if (CS%use_ideal_age) call ideal_age_example_end(CS%ideal_age_tracer_CSp)
  if (CS%use_regional_dyes) call regional_dyes_end(CS%dye_tracer_CSp)
  if (CS%use_oil) call oil_tracer_end(CS%oil_tracer_CSp)
  if (CS%use_advection_test_tracer) call advection_test_tracer_end(CS%advection_test_tracer_CSp)
  if (CS%use_OCMIP2_CFC) call OCMIP2_CFC_end(CS%OCMIP2_CFC_CSp)
#ifdef _USE_GENERIC_TRACER
  if (CS%use_MOM_generic_tracer) call end_MOM_generic_tracer(CS%MOM_generic_tracer_CSp)
#endif
  if (CS%use_pseudo_salt_tracer) call pseudo_salt_tracer_end(CS%pseudo_salt_tracer_CSp)
  if (CS%use_boundary_impulse_tracer) call boundary_impulse_tracer_end(CS%boundary_impulse_tracer_CSp)
  if (CS%use_dyed_obc_tracer) call dyed_obc_tracer_end(CS%dyed_obc_tracer_CSp)

  if (associated(CS)) deallocate(CS)
end subroutine tracer_flow_control_end

!> \namespace MOM_tracer_flow_control
!!
!!  By Will Cooke, April 2003
!!  Edited by Elizabeth Yankovsky, May 2019
!!
!!    This module contains two subroutines into which calls to other
!!  tracer initialization (call_tracer_init_fns) and column physics
!!  routines (call_tracer_column_fns) can be inserted.
!!
end module MOM_tracer_flow_control
