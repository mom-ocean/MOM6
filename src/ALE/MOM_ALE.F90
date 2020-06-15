!> This module contains the main regridding routines.
!!
!! Regridding comprises two steps:
!! 1. Interpolation and creation of a new grid based on target interface
!!    densities (or any other criterion).
!! 2. Remapping of quantities between old grid and new grid.
!!
!! Original module written by Laurent White, 2008.06.09
module MOM_ALE

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,        only : check_column_integrals
use MOM_diag_mediator,    only : register_diag_field, post_data, diag_ctrl
use MOM_diag_mediator,    only : time_type, diag_update_remap_grids
use MOM_diag_vkernels,    only : interpolate_column, reintegrate_column
use MOM_domains,          only : create_group_pass, do_group_pass, group_pass_type
use MOM_EOS,              only : calculate_density
use MOM_domains,          only : create_group_pass, do_group_pass, group_pass_type
use MOM_error_handler,    only : MOM_error, FATAL, WARNING
use MOM_error_handler,    only : callTree_showQuery
use MOM_error_handler,    only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,      only : get_param, param_file_type, log_param
use MOM_io,               only : vardesc, var_desc, fieldtype, SINGLE_FILE
use MOM_io,               only : create_file, write_field, close_file
use MOM_interface_heights,only : find_eta
use MOM_open_boundary,    only : ocean_OBC_type, OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary,    only : OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_regridding,       only : initialize_regridding, regridding_main, end_regridding
use MOM_regridding,       only : uniformResolution
use MOM_regridding,       only : inflate_vanished_layers_old
use MOM_regridding,       only : set_target_densities_from_GV, set_target_densities
use MOM_regridding,       only : regriddingCoordinateModeDoc, DEFAULT_COORDINATE_MODE
use MOM_regridding,       only : regriddingInterpSchemeDoc, regriddingDefaultInterpScheme
use MOM_regridding,       only : regriddingDefaultBoundaryExtrapolation
use MOM_regridding,       only : regriddingDefaultMinThickness
use MOM_regridding,       only : regridding_CS, set_regrid_params
use MOM_regridding,       only : getCoordinateInterfaces, getCoordinateResolution
use MOM_regridding,       only : getCoordinateUnits, getCoordinateShortName
use MOM_regridding,       only : getStaticThickness
use MOM_remapping,        only : initialize_remapping, end_remapping
use MOM_remapping,        only : remapping_core_h, remapping_core_w
use MOM_remapping,        only : remappingSchemesDoc, remappingDefaultScheme
use MOM_remapping,        only : remapping_CS, dzFromH1H2
use MOM_string_functions, only : uppercase, extractWord, extract_integer
use MOM_tracer_registry,  only : tracer_registry_type, tracer_type, MOM_tracer_chkinv
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,     only : get_thickness_units, verticalGrid_type

!use regrid_consts,       only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts,        only : coordinateUnits, coordinateMode, state_dependent
use regrid_edge_values,   only : edge_values_implicit_h4
use PLM_functions,        only : PLM_reconstruction, PLM_boundary_extrapolation
use PLM_functions,        only : PLM_extrapolate_slope, PLM_monotonized_slope, PLM_slope_wa
use PPM_functions,        only : PPM_reconstruction, PPM_boundary_extrapolation

implicit none ; private
#include <MOM_memory.h>


!> ALE control structure
type, public :: ALE_CS ; private
  logical :: remap_uv_using_old_alg !< If true, uses the old "remapping via a delta z"
                                    !! method. If False, uses the new method that
                                    !! remaps between grids described by h.

  real :: regrid_time_scale !< The time-scale used in blending between the current (old) grid
                            !! and the target (new) grid [T ~> s]

  type(regridding_CS) :: regridCS !< Regridding parameters and work arrays
  type(remapping_CS)  :: remapCS  !< Remapping parameters and work arrays

  integer :: nk             !< Used only for queries, not directly by this module

  logical :: remap_after_initialization !< Indicates whether to regrid/remap after initializing the state.

  logical :: answers_2018   !< If true, use the order of arithmetic and expressions for remapping
                            !! that recover the answers from the end of 2018.  Otherwise, use more
                            !! robust and accurate forms of mathematically equivalent expressions.

  logical :: show_call_tree !< For debugging

  ! for diagnostics
  type(diag_ctrl), pointer           :: diag                          !< structure to regulate output
  integer, dimension(:), allocatable :: id_tracer_remap_tendency      !< diagnostic id
  integer, dimension(:), allocatable :: id_Htracer_remap_tendency     !< diagnostic id
  integer, dimension(:), allocatable :: id_Htracer_remap_tendency_2d  !< diagnostic id
  logical, dimension(:), allocatable :: do_tendency_diag              !< flag for doing diagnostics
  integer                            :: id_dzRegrid = -1              !< diagnostic id

  ! diagnostic for fields prior to applying ALE remapping
  integer :: id_u_preale = -1 !< diagnostic id for zonal velocity before ALE.
  integer :: id_v_preale = -1 !< diagnostic id for meridional velocity before ALE.
  integer :: id_h_preale = -1 !< diagnostic id for layer thicknesses before ALE.
  integer :: id_T_preale = -1 !< diagnostic id for temperatures before ALE.
  integer :: id_S_preale = -1 !< diagnostic id for salinities before ALE.
  integer :: id_e_preale = -1 !< diagnostic id for interface heights before ALE.
  integer :: id_vert_remap_h = -1      !< diagnostic id for layer thicknesses used for remapping
  integer :: id_vert_remap_h_tendency = -1 !< diagnostic id for layer thickness tendency due to ALE

end type

! Publicly available functions
public ALE_init
public ALE_end
public ALE_main
public ALE_main_offline
public ALE_offline_inputs
public ALE_offline_tracer_final
public ALE_build_grid
public ALE_regrid_accelerated
public ALE_remap_scalar
public ALE_PLM_edge_values
public TS_PLM_edge_values
public TS_PPM_edge_values
public adjustGridForIntegrity
public ALE_initRegridding
public ALE_getCoordinate
public ALE_getCoordinateUnits
public ALE_writeCoordinateFile
public ALE_updateVerticalGridType
public ALE_initThicknessToCoord
public ALE_update_regrid_weights
public ALE_remap_init_conds
public ALE_register_diags

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> This routine is typically called (from initialize_MOM in file MOM.F90)
!! before the main time integration loop to initialize the regridding stuff.
!! We read the MOM_input file to register the values of different
!! regridding/remapping parameters.
subroutine ALE_init( param_file, GV, US, max_depth, CS)
  type(param_file_type),   intent(in) :: param_file !< Parameter file
  type(verticalGrid_type), intent(in) :: GV         !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US         !< A dimensional unit scaling type
  real,                    intent(in) :: max_depth  !< The maximum depth of the ocean [Z ~> m].
  type(ALE_CS),            pointer    :: CS         !< Module control structure

  ! Local variables
  real, dimension(:), allocatable :: dz
  character(len=40)               :: mdl = "MOM_ALE" ! This module's name.
  character(len=80)               :: string ! Temporary strings
  real                            :: filter_shallow_depth, filter_deep_depth
  logical       :: default_2018_answers
  logical                         :: check_reconstruction
  logical                         :: check_remapping
  logical                         :: force_bounds_in_subcell
  logical                         :: local_logical
  logical                         :: remap_boundary_extrap

  if (associated(CS)) then
    call MOM_error(WARNING, "ALE_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%show_call_tree = callTree_showQuery()
  if (CS%show_call_tree) call callTree_enter("ALE_init(), MOM_ALE.F90")

  call get_param(param_file, mdl, "REMAP_UV_USING_OLD_ALG", CS%remap_uv_using_old_alg, &
                 "If true, uses the old remapping-via-a-delta-z method for "//&
                 "remapping u and v. If false, uses the new method that remaps "//&
                 "between grids described by an old and new thickness.", &
                 default=.false.)

  ! Initialize and configure regridding
  call ALE_initRegridding(GV, US, max_depth, param_file, mdl, CS%regridCS)

  ! Initialize and configure remapping
  call get_param(param_file, mdl, "REMAPPING_SCHEME", string, &
                 "This sets the reconstruction scheme used "//&
                 "for vertical remapping for all variables. "//&
                 "It can be one of the following schemes: "//&
                 trim(remappingSchemesDoc), default=remappingDefaultScheme)
  call get_param(param_file, mdl, "FATAL_CHECK_RECONSTRUCTIONS", check_reconstruction, &
                 "If true, cell-by-cell reconstructions are checked for "//&
                 "consistency and if non-monotonicity or an inconsistency is "//&
                 "detected then a FATAL error is issued.", default=.false.)
  call get_param(param_file, mdl, "FATAL_CHECK_REMAPPING", check_remapping, &
                 "If true, the results of remapping are checked for "//&
                 "conservation and new extrema and if an inconsistency is "//&
                 "detected then a FATAL error is issued.", default=.false.)
  call get_param(param_file, mdl, "REMAP_BOUND_INTERMEDIATE_VALUES", force_bounds_in_subcell, &
                 "If true, the values on the intermediate grid used for remapping "//&
                 "are forced to be bounded, which might not be the case due to "//&
                 "round off.", default=.false.)
  call get_param(param_file, mdl, "REMAP_BOUNDARY_EXTRAP", remap_boundary_extrap, &
                 "If true, values at the interfaces of boundary cells are "//&
                 "extrapolated instead of piecewise constant", default=.false.)
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  call get_param(param_file, mdl, "REMAPPING_2018_ANSWERS", CS%answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
  call initialize_remapping( CS%remapCS, string, &
                             boundary_extrapolation=remap_boundary_extrap, &
                             check_reconstruction=check_reconstruction, &
                             check_remapping=check_remapping, &
                             force_bounds_in_subcell=force_bounds_in_subcell, &
                             answers_2018=CS%answers_2018)

  call get_param(param_file, mdl, "REMAP_AFTER_INITIALIZATION", CS%remap_after_initialization, &
                 "If true, applies regridding and remapping immediately after "//&
                 "initialization so that the state is ALE consistent. This is a "//&
                 "legacy step and should not be needed if the initialization is "//&
                 "consistent with the coordinate mode.", default=.true.)

  call get_param(param_file, mdl, "REGRID_TIME_SCALE", CS%regrid_time_scale, &
                 "The time-scale used in blending between the current (old) grid "//&
                 "and the target (new) grid. A short time-scale favors the target "//&
                 "grid (0. or anything less than DT_THERM) has no memory of the old "//&
                 "grid. A very long time-scale makes the model more Lagrangian.", &
                 units="s", default=0., scale=US%s_to_T)
  call get_param(param_file, mdl, "REGRID_FILTER_SHALLOW_DEPTH", filter_shallow_depth, &
                 "The depth above which no time-filtering is applied. Above this depth "//&
                 "final grid exactly matches the target (new) grid.", &
                 units="m", default=0., scale=GV%m_to_H)
  call get_param(param_file, mdl, "REGRID_FILTER_DEEP_DEPTH", filter_deep_depth, &
                 "The depth below which full time-filtering is applied with time-scale "//&
                 "REGRID_TIME_SCALE. Between depths REGRID_FILTER_SHALLOW_DEPTH and "//&
                 "REGRID_FILTER_SHALLOW_DEPTH the filter weights adopt a cubic profile.", &
                 units="m", default=0., scale=GV%m_to_H)
  call set_regrid_params(CS%regridCS, depth_of_time_filter_shallow=filter_shallow_depth, &
                         depth_of_time_filter_deep=filter_deep_depth)
  call get_param(param_file, mdl, "REGRID_USE_OLD_DIRECTION", local_logical, &
                 "If true, the regridding ntegrates upwards from the bottom for "//&
                 "interface positions, much as the main model does. If false "//&
                 "regridding integrates downward, consistant with the remapping "//&
                 "code.", default=.true., do_not_log=.true.)
  call set_regrid_params(CS%regridCS, integrate_downward_for_e=.not.local_logical)

  ! Keep a record of values for subsequent queries
  CS%nk = GV%ke

  if (CS%show_call_tree) call callTree_leave("ALE_init()")
end subroutine ALE_init

!> Initialize diagnostics for the ALE module.
subroutine ALE_register_diags(Time, G, GV, US, diag, CS)
  type(time_type),target,     intent(in)  :: Time  !< Time structure
  type(ocean_grid_type),      intent(in)  :: G     !< Grid structure
  type(unit_scale_type),      intent(in)  :: US    !< A dimensional unit scaling type
  type(verticalGrid_type),    intent(in)  :: GV    !< Ocean vertical grid structure
  type(diag_ctrl), target,    intent(in)  :: diag  !< Diagnostics control structure
  type(ALE_CS), pointer                   :: CS    !< Module control structure

  CS%diag => diag

  ! These diagnostics of the state variables before ALE are useful for
  ! debugging the ALE code.
  CS%id_u_preale = register_diag_field('ocean_model', 'u_preale', diag%axesCuL, Time, &
      'Zonal velocity before remapping', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_v_preale = register_diag_field('ocean_model', 'v_preale', diag%axesCvL, Time, &
      'Meridional velocity before remapping', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_h_preale = register_diag_field('ocean_model', 'h_preale', diag%axesTL, Time, &
      'Layer Thickness before remapping', get_thickness_units(GV), &
      conversion=GV%H_to_MKS, v_extensive=.true.)
  CS%id_T_preale = register_diag_field('ocean_model', 'T_preale', diag%axesTL, Time, &
      'Temperature before remapping', 'degC')
  CS%id_S_preale = register_diag_field('ocean_model', 'S_preale', diag%axesTL, Time, &
      'Salinity before remapping', 'PSU')
  CS%id_e_preale = register_diag_field('ocean_model', 'e_preale', diag%axesTi, Time, &
      'Interface Heights before remapping', 'm', conversion=US%Z_to_m)

  CS%id_dzRegrid = register_diag_field('ocean_model','dzRegrid',diag%axesTi,Time, &
      'Change in interface height due to ALE regridding', 'm', &
      conversion=GV%H_to_m)
  cs%id_vert_remap_h = register_diag_field('ocean_model', 'vert_remap_h', &
      diag%axestl, time, 'layer thicknesses after ALE regridding and remapping', 'm', &
      conversion=GV%H_to_m, v_extensive=.true.)
  cs%id_vert_remap_h_tendency = register_diag_field('ocean_model','vert_remap_h_tendency',diag%axestl,time, &
      'Layer thicknesses tendency due to ALE regridding and remapping', 'm', &
      conversion=GV%H_to_m*US%s_to_T, v_extensive = .true.)

end subroutine ALE_register_diags

!> Crudely adjust (initial) grid for integrity.
!! This routine is typically called (from initialize_MOM in file MOM.F90)
!! before the main time integration loop to initialize the regridding stuff.
!! We read the MOM_input file to register the values of different
!! regridding/remapping parameters.
subroutine adjustGridForIntegrity( CS, G, GV, h )
  type(ALE_CS),                              pointer       :: CS  !< Regridding parameters and options
  type(ocean_grid_type),                     intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                   intent(in)    :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h   !< Current 3D grid thickness that
                                                                  !! are to be adjusted [H ~> m or kg-2]
  call inflate_vanished_layers_old( CS%regridCS, G, GV, h(:,:,:) )

end subroutine adjustGridForIntegrity


!> End of regridding (memory deallocation).
!! This routine is typically called (from MOM_end in file MOM.F90)
!! after the main time integration loop to deallocate the regridding stuff.
subroutine ALE_end(CS)
  type(ALE_CS), pointer :: CS  !< module control structure

  ! Deallocate memory used for the regridding
  call end_remapping( CS%remapCS )
  call end_regridding( CS%regridCS )

  deallocate(CS)

end subroutine ALE_end

!> Takes care of (1) building a new grid and (2) remapping all variables between
!! the old grid and the new grid. The creation of the new grid can be based
!! on z coordinates, target interface densities, sigma coordinates or any
!! arbitrary coordinate system.
subroutine ALE_main( G, GV, US, h, u, v, tv, Reg, CS, OBC, dt, frac_shelf_h)
  type(ocean_grid_type),                      intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                    intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h   !< Current 3D grid obtained after the
                                                                   !! last time step [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u   !< Zonal velocity field [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v   !< Meridional velocity field [L T-1 ~> m s-1]
  type(thermo_var_ptrs),                      intent(inout) :: tv  !< Thermodynamic variable structure
  type(tracer_registry_type),                 pointer       :: Reg !< Tracer registry structure
  type(ALE_CS),                               pointer       :: CS  !< Regridding parameters and options
  type(ocean_OBC_type),                       pointer       :: OBC !< Open boundary structure
  real,                             optional, intent(in)    :: dt  !< Time step between calls to ALE_main [T ~> s]
  real, dimension(:,:),             optional, pointer       :: frac_shelf_h !< Fractional ice shelf coverage
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: dzRegrid ! The change in grid interface positions
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: eta_preale
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_new ! New 3D grid obtained after last time step [H ~> m or kg-2]
  integer :: nk, i, j, k, isc, iec, jsc, jec
  logical :: ice_shelf

  nk = GV%ke; isc = G%isc; iec = G%iec; jsc = G%jsc; jec = G%jec

  ice_shelf = .false.
  if (present(frac_shelf_h)) then
    if (associated(frac_shelf_h)) ice_shelf = .true.
  endif

  if (CS%show_call_tree) call callTree_enter("ALE_main(), MOM_ALE.F90")

  ! These diagnostics of the state before ALE is applied are mostly used for debugging.
  if (CS%id_u_preale > 0) call post_data(CS%id_u_preale, u,    CS%diag)
  if (CS%id_v_preale > 0) call post_data(CS%id_v_preale, v,    CS%diag)
  if (CS%id_h_preale > 0) call post_data(CS%id_h_preale, h,    CS%diag)
  if (CS%id_T_preale > 0) call post_data(CS%id_T_preale, tv%T, CS%diag)
  if (CS%id_S_preale > 0) call post_data(CS%id_S_preale, tv%S, CS%diag)
  if (CS%id_e_preale > 0) then
    call find_eta(h, tv, G, GV, US, eta_preale)
    call post_data(CS%id_e_preale, eta_preale, CS%diag)
  endif

  if (present(dt)) then
    call ALE_update_regrid_weights( dt, CS )
  endif
  dzRegrid(:,:,:) = 0.0

  ! Build new grid. The new grid is stored in h_new. The old grid is h.
  ! Both are needed for the subsequent remapping of variables.
  if (ice_shelf) then
     call regridding_main( CS%remapCS, CS%regridCS, G, GV, h, tv, h_new, dzRegrid, frac_shelf_h)
  else
     call regridding_main( CS%remapCS, CS%regridCS, G, GV, h, tv, h_new, dzRegrid)
  endif

  call check_grid( G, GV, h, 0. )

  if (CS%show_call_tree) call callTree_waypoint("new grid generated (ALE_main)")

  ! The presence of dt is used for expediency to distinguish whether ALE_main is being called during init
  ! or in the main loop. Tendency diagnostics in remap_all_state_vars also rely on this logic.
  if (present(dt)) then
    call diag_update_remap_grids(CS%diag)
  endif
  ! Remap all variables from old grid h onto new grid h_new
  call remap_all_state_vars( CS%remapCS, CS, G, GV, h, h_new, Reg, OBC, -dzRegrid, &
                             u, v, CS%show_call_tree, dt )

  if (CS%show_call_tree) call callTree_waypoint("state remapped (ALE_main)")

  ! Override old grid with new one. The new grid 'h_new' is built in
  ! one of the 'build_...' routines above.
  !$OMP parallel do default(shared)
  do k = 1,nk ; do j = jsc-1,jec+1 ; do i = isc-1,iec+1
    h(i,j,k) = h_new(i,j,k)
  enddo ; enddo ; enddo

  if (CS%show_call_tree) call callTree_leave("ALE_main()")

  if (CS%id_dzRegrid>0 .and. present(dt)) call post_data(CS%id_dzRegrid, dzRegrid, CS%diag)


end subroutine ALE_main

!> Takes care of (1) building a new grid and (2) remapping all variables between
!! the old grid and the new grid. The creation of the new grid can be based
!! on z coordinates, target interface densities, sigma coordinates or any
!! arbitrary coordinate system.
subroutine ALE_main_offline( G, GV, h, tv, Reg, CS, OBC, dt)
  type(ocean_grid_type),                      intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                    intent(in)    :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h   !< Current 3D grid obtained after the
                                                                   !! last time step [H ~> m or kg-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv  !< Thermodynamic variable structure
  type(tracer_registry_type),                 pointer       :: Reg !< Tracer registry structure
  type(ALE_CS),                               pointer       :: CS  !< Regridding parameters and options
  type(ocean_OBC_type),                       pointer       :: OBC !< Open boundary structure
  real,                             optional, intent(in)    :: dt  !< Time step between calls to ALE_main [T ~> s]
  ! Local variables
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)+1) :: dzRegrid ! The change in grid interface positions
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_new ! New 3D grid obtained after last time step [H ~> m or kg-2]
  integer :: nk, i, j, k, isc, iec, jsc, jec

  nk = GV%ke; isc = G%isc; iec = G%iec; jsc = G%jsc; jec = G%jec

  if (CS%show_call_tree) call callTree_enter("ALE_main_offline(), MOM_ALE.F90")

  if (present(dt)) then
    call ALE_update_regrid_weights( dt, CS )
  endif
  dzRegrid(:,:,:) = 0.0

  ! Build new grid. The new grid is stored in h_new. The old grid is h.
  ! Both are needed for the subsequent remapping of variables.
  call regridding_main( CS%remapCS, CS%regridCS, G, GV, h, tv, h_new, dzRegrid )

  call check_grid( G, GV, h, 0. )

  if (CS%show_call_tree) call callTree_waypoint("new grid generated (ALE_main)")

  ! Remap all variables from old grid h onto new grid h_new

  call remap_all_state_vars(CS%remapCS, CS, G, GV, h, h_new, Reg, OBC, &
                            debug=CS%show_call_tree, dt=dt )

  if (CS%show_call_tree) call callTree_waypoint("state remapped (ALE_main)")

  ! Override old grid with new one. The new grid 'h_new' is built in
  ! one of the 'build_...' routines above.
  !$OMP parallel do default(shared)
  do k = 1,nk ; do j = jsc-1,jec+1 ; do i = isc-1,iec+1
    h(i,j,k) = h_new(i,j,k)
  enddo ; enddo ; enddo

  if (CS%show_call_tree) call callTree_leave("ALE_main()")
  if (CS%id_dzRegrid>0 .and. present(dt)) call post_data(CS%id_dzRegrid, dzRegrid, CS%diag)

end subroutine ALE_main_offline

!> Regrid/remap stored fields used for offline tracer integrations. These input fields are assumed to have
!! the same layer thicknesses at the end of the last offline interval (which should be a Zstar grid). This
!! routine builds a grid on the runtime specified vertical coordinate
subroutine ALE_offline_inputs(CS, G, GV, h, tv, Reg, uhtr, vhtr, Kd, debug, OBC)
  type(ALE_CS),                                 pointer       :: CS    !< Regridding parameters and options
  type(ocean_grid_type),                        intent(in   ) :: G     !< Ocean grid informations
  type(verticalGrid_type),                      intent(in   ) :: GV    !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(inout) :: h     !< Layer thicknesses
  type(thermo_var_ptrs),                        intent(inout) :: tv    !< Thermodynamic variable structure
  type(tracer_registry_type),                   pointer       :: Reg   !< Tracer registry structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),   intent(inout) :: uhtr  !< Zonal mass fluxes
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),   intent(inout) :: vhtr  !< Meridional mass fluxes
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(inout) :: Kd    !< Input diffusivites
  logical,                                      intent(in   ) :: debug !< If true, then turn checksums
  type(ocean_OBC_type),                         pointer       :: OBC   !< Open boundary structure
  ! Local variables
  integer :: nk, i, j, k, isc, iec, jsc, jec
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV))   :: h_new    ! Layer thicknesses after regridding
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)+1) :: dzRegrid ! The change in grid interface positions
  real, dimension(SZK_(GV)) :: h_src
  real, dimension(SZK_(GV)) :: h_dest, uh_dest
  real, dimension(SZK_(GV)) :: temp_vec

  nk = GV%ke; isc = G%isc; iec = G%iec; jsc = G%jsc; jec = G%jec
  dzRegrid(:,:,:) = 0.0
  h_new(:,:,:) = 0.0

  if (debug) call MOM_tracer_chkinv("Before ALE_offline_inputs", G, h, Reg%Tr, Reg%ntr)

  ! Build new grid from the Zstar state onto the requested vertical coordinate. The new grid is stored
  ! in h_new. The old grid is h. Both are needed for the subsequent remapping of variables. Convective
  ! adjustment right now is not used because it is unclear what to do with vanished layers
  call regridding_main( CS%remapCS, CS%regridCS, G, GV, h, tv, h_new, dzRegrid, conv_adjust = .false. )
  call check_grid( G, GV, h_new, 0. )
  if (CS%show_call_tree) call callTree_waypoint("new grid generated (ALE_offline_inputs)")

  ! Remap all variables from old grid h onto new grid h_new
  call remap_all_state_vars( CS%remapCS, CS, G, GV, h, h_new, Reg, OBC, debug=CS%show_call_tree )
  if (CS%show_call_tree) call callTree_waypoint("state remapped (ALE_inputs)")

  ! Reintegrate mass transports from Zstar to the offline vertical coordinate
  do j=jsc,jec ; do i=G%iscB,G%iecB
    if (G%mask2dCu(i,j)>0.) then
      h_src(:) = 0.5 * (h(i,j,:) + h(i+1,j,:))
      h_dest(:) = 0.5 * (h_new(i,j,:) + h_new(i+1,j,:))
      call reintegrate_column(nk, h_src, uhtr(I,j,:), nk, h_dest, 0., temp_vec)
      uhtr(I,j,:) = temp_vec
    endif
  enddo ; enddo
  do j=G%jscB,G%jecB ; do i=isc,iec
    if (G%mask2dCv(i,j)>0.) then
      h_src(:) = 0.5 * (h(i,j,:) + h(i,j+1,:))
      h_dest(:) = 0.5 * (h_new(i,j,:) + h_new(i,j+1,:))
      call reintegrate_column(nk, h_src, vhtr(I,j,:), nk, h_dest, 0., temp_vec)
      vhtr(I,j,:) = temp_vec
    endif
  enddo ; enddo

  do j = jsc,jec ; do i=isc,iec
    if (G%mask2dT(i,j)>0.) then
      if (check_column_integrals(nk, h_src, nk, h_dest)) then
        call MOM_error(FATAL, "ALE_offline_inputs: Kd interpolation columns do not match")
      endif
      call interpolate_column(nk, h(i,j,:), Kd(i,j,:), nk, h_new(i,j,:), 0., Kd(i,j,:))
    endif
  enddo ; enddo

  call ALE_remap_scalar(CS%remapCS, G, GV, nk, h, tv%T, h_new, tv%T, answers_2018=CS%answers_2018)
  call ALE_remap_scalar(CS%remapCS, G, GV, nk, h, tv%S, h_new, tv%S, answers_2018=CS%answers_2018)

  if (debug) call MOM_tracer_chkinv("After ALE_offline_inputs", G, h_new, Reg%Tr, Reg%ntr)

  ! Copy over the new layer thicknesses
  do k = 1,nk  ; do j = jsc-1,jec+1 ; do i = isc-1,iec+1
      h(i,j,k) = h_new(i,j,k)
  enddo ; enddo ; enddo

  if (CS%show_call_tree) call callTree_leave("ALE_offline_inputs()")
end subroutine ALE_offline_inputs


!> Remaps all tracers from h onto h_target. This is intended to be called when tracers
!! are done offline. In the case where transports don't quite conserve, we still want to
!! make sure that layer thicknesses offline do not drift too far away from the online model
subroutine ALE_offline_tracer_final( G, GV, h, tv, h_target, Reg, CS, OBC)
  type(ocean_grid_type),                      intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                    intent(in)    :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h   !< Current 3D grid obtained after the
                                                                   !! last time step [H ~> m or kg-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv  !< Thermodynamic variable structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h_target !< Current 3D grid obtained after
                                                                        !! last time step  [H ~> m or kg-2]
  type(tracer_registry_type),                 pointer       :: Reg !< Tracer registry structure
  type(ALE_CS),                               pointer       :: CS  !< Regridding parameters and options
  type(ocean_OBC_type),                       pointer       :: OBC !< Open boundary structure
  ! Local variables

  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)+1) :: dzRegrid !< The change in grid interface positions
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV))   :: h_new    !< Regridded target thicknesses
  integer :: nk, i, j, k, isc, iec, jsc, jec

  nk = GV%ke; isc = G%isc; iec = G%iec; jsc = G%jsc; jec = G%jec

  if (CS%show_call_tree) call callTree_enter("ALE_offline_tracer_final(), MOM_ALE.F90")
  ! Need to make sure that h_target is consistent with the current offline ALE confiuration
  call regridding_main( CS%remapCS, CS%regridCS, G, GV, h_target, tv, h_new, dzRegrid )
  call check_grid( G, GV, h_target, 0. )


  if (CS%show_call_tree) call callTree_waypoint("Source and target grids checked (ALE_offline_tracer_final)")

  ! Remap all variables from old grid h onto new grid h_new

  call remap_all_state_vars( CS%remapCS, CS, G, GV, h, h_new, Reg, OBC, debug=CS%show_call_tree )

  if (CS%show_call_tree) call callTree_waypoint("state remapped (ALE_offline_tracer_final)")

  ! Override old grid with new one. The new grid 'h_new' is built in
  ! one of the 'build_...' routines above.
  !$OMP parallel do default(shared)
  do k = 1,nk
    do j = jsc-1,jec+1 ; do i = isc-1,iec+1
      h(i,j,k) = h_new(i,j,k)
    enddo ; enddo
  enddo
  if (CS%show_call_tree) call callTree_leave("ALE_offline_tracer_final()")
end subroutine ALE_offline_tracer_final

!> Check grid for negative thicknesses
subroutine check_grid( G, GV, h, threshold )
  type(ocean_grid_type),                     intent(in) :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in) :: GV !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h  !< Current 3D grid obtained after the
                                                              !! last time step [H ~> m or kg m-2]
  real,                                      intent(in) :: threshold !< Value below which to flag issues,
                                                              !! [H ~> m or kg m-2]
  ! Local variables
  integer :: i, j

  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    if (G%mask2dT(i,j)>0.) then
      if (minval(h(i,j,:)) < threshold) then
        write(0,*) 'check_grid: i,j=',i,j,'h(i,j,:)=',h(i,j,:)
        if (threshold <= 0.) then
          call MOM_error(FATAL,"MOM_ALE, check_grid: negative thickness encountered.")
        else
          call MOM_error(FATAL,"MOM_ALE, check_grid: too tiny thickness encountered.")
        endif
      endif
    endif
  enddo ; enddo


end subroutine check_grid

!> Generates new grid
subroutine ALE_build_grid( G, GV, regridCS, remapCS, h, tv, debug, frac_shelf_h )
  type(ocean_grid_type),                   intent(in)    :: G        !< Ocean grid structure
  type(verticalGrid_type),                 intent(in)    :: GV       !< Ocean vertical grid structure
  type(regridding_CS),                     intent(in)    :: regridCS !< Regridding parameters and options
  type(remapping_CS),                      intent(in)    :: remapCS  !< Remapping parameters and options
  type(thermo_var_ptrs),                   intent(inout) :: tv       !< Thermodynamical variable structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)), intent(inout) :: h     !< Current 3D grid obtained after the
                                                                     !! last time step [H ~> m or kg-2]
  logical,                       optional, intent(in)    :: debug    !< If true, show the call tree
  real, dimension(:,:),          optional, pointer       :: frac_shelf_h !< Fractional ice shelf coverage
  ! Local variables
  integer :: nk, i, j, k
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)+1) :: dzRegrid ! The change in grid interface positions
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)) :: h_new ! The new grid thicknesses
  logical :: show_call_tree, use_ice_shelf

  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug
  if (show_call_tree) call callTree_enter("ALE_build_grid(), MOM_ALE.F90")
  use_ice_shelf = .false.
  if (present(frac_shelf_h)) then
    if (associated(frac_shelf_h)) use_ice_shelf = .true.
  endif

  ! Build new grid. The new grid is stored in h_new. The old grid is h.
  ! Both are needed for the subsequent remapping of variables.
  if (use_ice_shelf) then
     call regridding_main( remapCS, regridCS, G, GV, h, tv, h_new, dzRegrid, frac_shelf_h )
  else
     call regridding_main( remapCS, regridCS, G, GV, h, tv, h_new, dzRegrid )
  endif

  ! Override old grid with new one. The new grid 'h_new' is built in
  ! one of the 'build_...' routines above.
!$OMP parallel do default(none) shared(G,h,h_new)
  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    if (G%mask2dT(i,j)>0.) h(i,j,:) = h_new(i,j,:)
  enddo ; enddo

  if (show_call_tree) call callTree_leave("ALE_build_grid()")
end subroutine ALE_build_grid

!> For a state-based coordinate, accelerate the process of regridding by
!! repeatedly applying the grid calculation algorithm
subroutine ALE_regrid_accelerated(CS, G, GV, h, tv, n, u, v, OBC, Reg, dt, dzRegrid, initial)
  type(ALE_CS),            pointer       :: CS     !< ALE control structure
  type(ocean_grid_type),   intent(inout) :: G      !< Ocean grid
  type(verticalGrid_type), intent(in)    :: GV     !< Vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h      !< Original thicknesses [H ~> m or kg-2]
  type(thermo_var_ptrs),   intent(inout) :: tv     !< Thermo vars (T/S/EOS)
  integer,                 intent(in)    :: n      !< Number of times to regrid
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v      !< Meridional velocity [L T-1 ~> m s-1]
  type(ocean_OBC_type),    pointer       :: OBC    !< Open boundary structure
  type(tracer_registry_type), &
                 optional, pointer       :: Reg    !< Tracer registry to remap onto new grid
  real,          optional, intent(in)    :: dt     !< Model timestep to provide a timescale for regridding [T ~> s]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                 optional, intent(inout) :: dzRegrid !< Final change in interface positions
  logical,       optional, intent(in)    :: initial !< Whether we're being called from an initialization
                                                    !! routine (and expect diagnostics to work)

  ! Local variables
  integer :: i, j, k, nz
  type(thermo_var_ptrs) :: tv_local ! local/intermediate temp/salt
  type(group_pass_type) :: pass_T_S_h ! group pass if the coordinate has a stencil
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))         :: h_loc, h_orig ! A working copy of layer thicknesses
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target :: T, S ! local temporary state
  ! we have to keep track of the total dzInterface if for some reason
  ! we're using the old remapping algorithm for u/v
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: dzInterface, dzIntTotal

  nz = GV%ke

  ! initial total interface displacement due to successive regridding
  dzIntTotal(:,:,:) = 0.

  call create_group_pass(pass_T_S_h, T, G%domain)
  call create_group_pass(pass_T_S_h, S, G%domain)
  call create_group_pass(pass_T_S_h, h_loc, G%domain)

  ! copy original temp/salt and set our local tv_pointers to them
  tv_local = tv
  T(:,:,:) = tv%T(:,:,:)
  S(:,:,:) = tv%S(:,:,:)
  tv_local%T => T
  tv_local%S => S

  ! get local copy of thickness and save original state for remapping
  h_loc(:,:,:) = h(:,:,:)
  h_orig(:,:,:) = h(:,:,:)

  ! Apply timescale to regridding (for e.g. filtered_grid_motion)
  if (present(dt)) &
    call ALE_update_regrid_weights(dt, CS)

  do k = 1, n
    call do_group_pass(pass_T_S_h, G%domain)

    ! generate new grid
    call regridding_main(CS%remapCS, CS%regridCS, G, GV, h_loc, tv_local, h, dzInterface)
    dzIntTotal(:,:,:) = dzIntTotal(:,:,:) + dzInterface(:,:,:)

    ! remap from original grid onto new grid
    do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1
      call remapping_core_h(CS%remapCS, nz, h_orig(i,j,:), tv%S(i,j,:), nz, h(i,j,:), tv_local%S(i,j,:))
      call remapping_core_h(CS%remapCS, nz, h_orig(i,j,:), tv%T(i,j,:), nz, h(i,j,:), tv_local%T(i,j,:))
    enddo ; enddo

    ! starting grid for next iteration
    h_loc(:,:,:) = h(:,:,:)
  enddo

  ! remap all state variables (including those that weren't needed for regridding)
  call remap_all_state_vars(CS%remapCS, CS, G, GV, h_orig, h, Reg, OBC, dzIntTotal, u, v)

  ! save total dzregrid for diags if needed?
  if (present(dzRegrid)) dzRegrid(:,:,:) = dzIntTotal(:,:,:)
end subroutine ALE_regrid_accelerated

!> This routine takes care of remapping all variable between the old and the
!! new grids. When velocity components need to be remapped, thicknesses at
!! velocity points are taken to be arithmetic averages of tracer thicknesses.
!! This routine is called during initialization of the model at time=0, to
!! remap initiali conditions to the model grid.  It is also called during a
!! time step to update the state.
subroutine remap_all_state_vars(CS_remapping, CS_ALE, G, GV, h_old, h_new, Reg, OBC, &
                                dxInterface, u, v, debug, dt)
  type(remapping_CS),                        intent(in)    :: CS_remapping !< Remapping control structure
  type(ALE_CS),                              intent(in)    :: CS_ALE       !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G            !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV           !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_old        !< Thickness of source grid
                                                                           !! [H ~> m or kg-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_new        !< Thickness of destination grid
                                                                           !! [H ~> m or kg-2]
  type(tracer_registry_type),                pointer       :: Reg          !< Tracer registry structure
  type(ocean_OBC_type),                      pointer       :: OBC          !< Open boundary structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                                   optional, intent(in)    :: dxInterface  !< Change in interface position
                                                                           !! [H ~> m or kg-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                   optional, intent(inout) :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                   optional, intent(inout) :: v      !< Meridional velocity [L T-1 ~> m s-1]
  logical,                         optional, intent(in)    :: debug  !< If true, show the call tree
  real,                            optional, intent(in)    :: dt     !< time step for diagnostics [T ~> s]
  ! Local variables
  integer                                     :: i, j, k, m
  integer                                     :: nz, ntr
  real, dimension(GV%ke+1)                    :: dx
  real, dimension(GV%ke)                      :: h1, u_column
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)) :: work_conc
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)) :: work_cont
  real, dimension(SZI_(G), SZJ_(G))           :: work_2d
  real                                        :: Idt ! The inverse of the timestep [T-1 ~> s-1]
  real                                        :: ppt2mks
  real, dimension(GV%ke)                      :: h2
  real :: h_neglect, h_neglect_edge
  logical                                     :: show_call_tree
  type(tracer_type), pointer                  :: Tr => NULL()

  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug
  if (show_call_tree) call callTree_enter("remap_all_state_vars(), MOM_ALE.F90")

  ! If remap_uv_using_old_alg is .true. and u or v is requested, then we must have dxInterface. Otherwise,
  ! u and v can be remapped without dxInterface
  if ( .not. present(dxInterface) .and. (CS_ALE%remap_uv_using_old_alg .and. (present(u) .or. present(v))) ) then
    call MOM_error(FATAL, "remap_all_state_vars: dxInterface must be present if using old algorithm "// &
                          "and u/v are to be remapped")
  endif

  if (.not.CS_ALE%answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif

  nz      = GV%ke
  ppt2mks = 0.001

  ntr = 0 ; if (associated(Reg)) ntr = Reg%ntr

  if (present(dt)) then
    Idt = 1.0/dt
    work_conc(:,:,:) = 0.0
    work_cont(:,:,:) = 0.0
  endif

  ! Remap tracer
  if (ntr>0) then
    if (show_call_tree) call callTree_waypoint("remapping tracers (remap_all_state_vars)")
    !$OMP parallel do default(shared) private(h1,h2,u_column,Tr)
    do m=1,ntr ! For each tracer
      Tr => Reg%Tr(m)
      do j = G%jsc,G%jec ; do i = G%isc,G%iec ; if (G%mask2dT(i,j)>0.) then
        ! Build the start and final grids
        h1(:) = h_old(i,j,:)
        h2(:) = h_new(i,j,:)
        call remapping_core_h(CS_remapping, nz, h1, Tr%t(i,j,:), nz, h2, &
                              u_column, h_neglect, h_neglect_edge)

        ! Intermediate steps for tendency of tracer concentration and tracer content.
        if (present(dt)) then
          if (Tr%id_remap_conc > 0) then
            do k=1,GV%ke
              work_conc(i,j,k) = (u_column(k) - Tr%t(i,j,k)) * Idt
            enddo
          endif
          if (Tr%id_remap_cont > 0 .or. Tr%id_remap_cont_2d > 0) then
            do k=1,GV%ke
              work_cont(i,j,k) = (u_column(k)*h2(k) - Tr%t(i,j,k)*h1(k)) * Idt
            enddo
          endif
        endif
        ! update tracer concentration
        Tr%t(i,j,:) = u_column(:)
      endif ; enddo ; enddo

      ! tendency diagnostics.
      if (present(dt)) then
        if (Tr%id_remap_conc > 0) then
          call post_data(Tr%id_remap_conc, work_conc, CS_ALE%diag)
        endif
        if (Tr%id_remap_cont > 0) then
          call post_data(Tr%id_remap_cont, work_cont, CS_ALE%diag)
        endif
        if (Tr%id_remap_cont_2d > 0) then
          do j = G%jsc,G%jec ; do i = G%isc,G%iec
            work_2d(i,j) = 0.0
            do k = 1,GV%ke
              work_2d(i,j) = work_2d(i,j) + work_cont(i,j,k)
            enddo
          enddo ; enddo
          call post_data(Tr%id_remap_cont_2d, work_2d, CS_ALE%diag)
        endif
      endif
    enddo ! m=1,ntr

  endif   ! endif for ntr > 0

  if (show_call_tree) call callTree_waypoint("tracers remapped (remap_all_state_vars)")

  ! Remap u velocity component
  if ( present(u) ) then
    !$OMP parallel do default(shared) private(h1,h2,dx,u_column)
    do j = G%jsc,G%jec ; do I = G%iscB,G%iecB ; if (G%mask2dCu(I,j)>0.) then
      ! Build the start and final grids
      h1(:) = 0.5 * ( h_old(i,j,:) + h_old(i+1,j,:) )
      if (CS_ALE%remap_uv_using_old_alg) then
        dx(:) = 0.5 * ( dxInterface(i,j,:) + dxInterface(i+1,j,:) )
        do k = 1, nz
          h2(k) = max( 0., h1(k) + ( dx(k+1) - dx(k) ) )
        enddo
      else
        h2(:) = 0.5 * ( h_new(i,j,:) + h_new(i+1,j,:) )
      endif
      if (associated(OBC)) then
        if (OBC%segnum_u(I,j) .ne. 0) then
          if (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_E) then
            h1(:) = h_old(i,j,:)
            h2(:) = h_new(i,j,:)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_W)
            h1(:) = h_old(i+1,j,:)
            h2(:) = h_new(i+1,j,:)
          endif
        endif
      endif
      call remapping_core_h(CS_remapping, nz, h1, u(I,j,:), nz, h2, &
                            u_column, h_neglect, h_neglect_edge)
      u(I,j,:) = u_column(:)
    endif ; enddo ; enddo
  endif

  if (show_call_tree) call callTree_waypoint("u remapped (remap_all_state_vars)")

  ! Remap v velocity component
  if ( present(v) ) then
    !$OMP parallel do default(shared) private(h1,h2,dx,u_column)
    do J = G%jscB,G%jecB ; do i = G%isc,G%iec ; if (G%mask2dCv(i,j)>0.) then
      ! Build the start and final grids
      h1(:) = 0.5 * ( h_old(i,j,:) + h_old(i,j+1,:) )
      if (CS_ALE%remap_uv_using_old_alg) then
        dx(:) = 0.5 * ( dxInterface(i,j,:) + dxInterface(i,j+1,:) )
        do k = 1, nz
          h2(k) = max( 0., h1(k) + ( dx(k+1) - dx(k) ) )
        enddo
      else
        h2(:) = 0.5 * ( h_new(i,j,:) + h_new(i,j+1,:) )
      endif
      if (associated(OBC)) then
        if (OBC%segnum_v(i,J) .ne. 0) then
          if (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_N) then
            h1(:) = h_old(i,j,:)
            h2(:) = h_new(i,j,:)
          else ! (OBC%segment(n)%direction == OBC_DIRECTION_S)
            h1(:) = h_old(i,j+1,:)
            h2(:) = h_new(i,j+1,:)
          endif
        endif
      endif
      call remapping_core_h(CS_remapping, nz, h1, v(i,J,:), nz, h2, &
                            u_column, h_neglect, h_neglect_edge)
      v(i,J,:) = u_column(:)
    endif ; enddo ; enddo
  endif

  if (CS_ALE%id_vert_remap_h > 0) call post_data(CS_ALE%id_vert_remap_h, h_old, CS_ALE%diag)
  if ((CS_ALE%id_vert_remap_h_tendency > 0) .and. present(dt)) then
    do k = 1, nz ; do j = G%jsc,G%jec ; do i = G%isc,G%iec
      work_cont(i,j,k) = (h_new(i,j,k) - h_old(i,j,k))*Idt
    enddo ; enddo ; enddo
    call post_data(CS_ALE%id_vert_remap_h_tendency, work_cont, CS_ALE%diag)
  endif
  if (show_call_tree) call callTree_waypoint("v remapped (remap_all_state_vars)")
  if (show_call_tree) call callTree_leave("remap_all_state_vars()")

end subroutine remap_all_state_vars


!> Remaps a single scalar between grids described by thicknesses h_src and h_dst.
!! h_dst must be dimensioned as a model array with GV%ke layers while h_src can
!! have an arbitrary number of layers specified by nk_src.
subroutine ALE_remap_scalar(CS, G, GV, nk_src, h_src, s_src, h_dst, s_dst, all_cells, old_remap, answers_2018 )
  type(remapping_CS),                      intent(in)    :: CS        !< Remapping control structure
  type(ocean_grid_type),                   intent(in)    :: G         !< Ocean grid structure
  type(verticalGrid_type),                 intent(in)    :: GV        !< Ocean vertical grid structure
  integer,                                 intent(in)    :: nk_src    !< Number of levels on source grid
  real, dimension(SZI_(G),SZJ_(G),nk_src), intent(in)    :: h_src     !< Level thickness of source grid
                                                                      !! [H ~> m or kg-2]
  real, dimension(SZI_(G),SZJ_(G),nk_src), intent(in)    :: s_src     !< Scalar on source grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),intent(in)   :: h_dst     !< Level thickness of destination grid
                                                                      !! [H ~> m or kg-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),intent(inout) :: s_dst    !< Scalar on destination grid
  logical, optional,                       intent(in)    :: all_cells !< If false, only reconstruct for
                                                                      !! non-vanished cells. Use all vanished
                                                                      !! layers otherwise (default).
  logical, optional,                       intent(in)    :: old_remap !< If true, use the old "remapping_core_w"
                                                                      !! method, otherwise use "remapping_core_h".
  logical,                       optional, intent(in)    :: answers_2018 !< If true, use the order of arithmetic
                                                                      !! and expressions that recover the answers for
                                                                      !! remapping from the end of 2018. Otherwise,
                                                                      !! use more robust forms of the same expressions.
  ! Local variables
  integer :: i, j, k, n_points
  real :: dx(GV%ke+1)
  real :: h_neglect, h_neglect_edge
  logical :: ignore_vanished_layers, use_remapping_core_w, use_2018_remap

  ignore_vanished_layers = .false.
  if (present(all_cells)) ignore_vanished_layers = .not. all_cells
  use_remapping_core_w = .false.
  if (present(old_remap)) use_remapping_core_w = old_remap
  n_points = nk_src
  use_2018_remap = .true. ; if (present(answers_2018)) use_2018_remap = answers_2018

  if (.not.use_2018_remap) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif

  !$OMP parallel do default(shared) firstprivate(n_points,dx)
  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    if (G%mask2dT(i,j) > 0.) then
      if (ignore_vanished_layers) then
        n_points = 0
        do k = 1, nk_src
          if (h_src(i,j,k)>0.) n_points = n_points + 1
        enddo
        s_dst(i,j,:) = 0.
      endif
      if (use_remapping_core_w) then
        call dzFromH1H2( n_points, h_src(i,j,1:n_points), GV%ke, h_dst(i,j,:), dx )
        call remapping_core_w(CS, n_points, h_src(i,j,1:n_points), s_src(i,j,1:n_points), &
                              GV%ke, dx, s_dst(i,j,:), h_neglect, h_neglect_edge)
      else
        call remapping_core_h(CS, n_points, h_src(i,j,1:n_points), s_src(i,j,1:n_points), &
                              GV%ke, h_dst(i,j,:), s_dst(i,j,:), h_neglect, h_neglect_edge)
      endif
    else
      s_dst(i,j,:) = 0.
    endif
  enddo ; enddo

end subroutine ALE_remap_scalar


!> Calculate edge values (top and bottom of layer) for T and S consistent with a PLM reconstruction
!! in the vertical direction. Boundary reconstructions are PCM unless bdry_extrap is true.
subroutine TS_PLM_edge_values( CS, S_t, S_b, T_t, T_b, G, GV, tv, h, bdry_extrap )
  type(ocean_grid_type),   intent(in)    :: G    !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean vertical grid structure
  type(ALE_CS),            intent(inout) :: CS   !< module control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: S_t  !< Salinity at the top edge of each layer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: S_b  !< Salinity at the bottom edge of each layer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_t  !< Temperature at the top edge of each layer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_b  !< Temperature at the bottom edge of each layer
  type(thermo_var_ptrs),   intent(in)    :: tv   !< thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< layer thickness [H ~> m or kg m-2]
  logical,                 intent(in)    :: bdry_extrap !< If true, use high-order boundary
                                                 !! extrapolation within boundary cells

  call ALE_PLM_edge_values( CS, G, GV, h, tv%S, bdry_extrap, S_t, S_b )
  call ALE_PLM_edge_values( CS, G, GV, h, tv%T, bdry_extrap, T_t, T_b )

end subroutine TS_PLM_edge_values

!> Calculate edge values (top and bottom of layer) 3d scalar array.
!! Boundary reconstructions are PCM unless bdry_extrap is true.
subroutine ALE_PLM_edge_values( CS, G, GV, h, Q, bdry_extrap, Q_t, Q_b )
  type(ALE_CS),            intent(in)    :: CS   !< module control structure
  type(ocean_grid_type),   intent(in)    :: G    !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: Q    !< 3d scalar array
  logical,                 intent(in)    :: bdry_extrap !< If true, use high-order boundary
                                                 !! extrapolation within boundary cells
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: Q_t  !< Scalar at the top edge of each layer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: Q_b  !< Scalar at the bottom edge of each layer
  ! Local variables
  integer :: i, j, k
  real :: slp(GV%ke)
  real :: mslp
  real :: h_neglect

  if (.not.CS%answers_2018) then
    h_neglect = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30
  endif

  !$OMP parallel do default(shared) private(slp,mslp)
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1
    slp(1) = 0.
    do k = 2, GV%ke-1
      slp(k) = PLM_slope_wa(h(i,j,k-1), h(i,j,k), h(i,j,k+1), h_neglect, Q(i,j,k-1), Q(i,j,k), Q(i,j,k+1))
    enddo
    slp(GV%ke) = 0.

    do k = 2, GV%ke-1
      mslp = PLM_monotonized_slope(Q(i,j,k-1), Q(i,j,k), Q(i,j,k+1), slp(k-1), slp(k), slp(k+1))
      Q_t(i,j,k) = Q(i,j,k) - 0.5 * mslp
      Q_b(i,j,k) = Q(i,j,k) + 0.5 * mslp
    enddo
    if (bdry_extrap) then
      mslp = - PLM_extrapolate_slope(h(i,j,2), h(i,j,1), h_neglect, Q(i,j,2), Q(i,j,1))
      Q_t(i,j,1) = Q(i,j,1) - 0.5 * mslp
      Q_b(i,j,1) = Q(i,j,1) + 0.5 * mslp
      mslp = PLM_extrapolate_slope(h(i,j,GV%ke-1), h(i,j,GV%ke), h_neglect, Q(i,j,GV%ke-1), Q(i,j,GV%ke))
      Q_t(i,j,GV%ke) = Q(i,j,GV%ke) - 0.5 * mslp
      Q_b(i,j,GV%ke) = Q(i,j,GV%ke) + 0.5 * mslp
    else
      Q_t(i,j,1) = Q(i,j,1)
      Q_b(i,j,1) = Q(i,j,1)
      Q_t(i,j,GV%ke) = Q(i,j,GV%ke)
      Q_b(i,j,GV%ke) = Q(i,j,GV%ke)
    endif

  enddo ; enddo

end subroutine ALE_PLM_edge_values

!> Calculate edge values (top and bottom of layer) for T and S consistent with a PPM reconstruction
!! in the vertical direction. Boundary reconstructions are PCM unless bdry_extrap is true.
subroutine TS_PPM_edge_values( CS, S_t, S_b, T_t, T_b, G, GV, tv, h, bdry_extrap )
  type(ocean_grid_type),   intent(in)    :: G    !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean vertical grid structure
  type(ALE_CS),            intent(inout) :: CS   !< module control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: S_t  !< Salinity at the top edge of each layer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: S_b  !< Salinity at the bottom edge of each layer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_t  !< Temperature at the top edge of each layer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_b  !< Temperature at the bottom edge of each layer
  type(thermo_var_ptrs),   intent(in)    :: tv   !< thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< layer thicknesses [H ~> m or kg m-2]
  logical,                 intent(in)    :: bdry_extrap !< If true, use high-order boundary
                                                 !! extrapolation within boundary cells

  ! Local variables
  integer :: i, j, k
  real    :: hTmp(GV%ke) ! A 1-d copy of h [H ~> m or kg m-2]
  real    :: tmp(GV%ke)  ! A 1-d copy of a column of temperature [degC] or salinity [ppt]
  real, dimension(CS%nk,2) :: &
      ppol_E            ! Edge value of polynomial in [degC] or [ppt]
  real, dimension(CS%nk,3) :: &
      ppol_coefs        ! Coefficients of polynomial, all in [degC] or [ppt]
  real :: h_neglect, h_neglect_edge ! Tiny thicknesses [H ~> m or kg m-2]

  if (.not.CS%answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif

  ! Determine reconstruction within each column
  !$OMP parallel do default(shared) private(hTmp,tmp,ppol_E,ppol_coefs)
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1

    ! Build current grid
    hTmp(:) = h(i,j,:)
    tmp(:) = tv%S(i,j,:)

    ! Reconstruct salinity profile
    ppol_E(:,:) = 0.0
    ppol_coefs(:,:) = 0.0
    call edge_values_implicit_h4( GV%ke, hTmp, tmp, ppol_E, h_neglect=h_neglect_edge, &
                                  answers_2018=CS%answers_2018 )
    call PPM_reconstruction( GV%ke, hTmp, tmp, ppol_E, ppol_coefs, h_neglect, &
                                  answers_2018=CS%answers_2018 )
    if (bdry_extrap) &
      call PPM_boundary_extrapolation( GV%ke, hTmp, tmp, ppol_E, ppol_coefs, h_neglect )

    do k = 1,GV%ke
      S_t(i,j,k) = ppol_E(k,1)
      S_b(i,j,k) = ppol_E(k,2)
    enddo

    ! Reconstruct temperature profile
    ppol_E(:,:) = 0.0
    ppol_coefs(:,:) = 0.0
    tmp(:) = tv%T(i,j,:)
    if (CS%answers_2018) then
      call edge_values_implicit_h4( GV%ke, hTmp, tmp, ppol_E, h_neglect=1.0e-10*GV%m_to_H, &
                                  answers_2018=CS%answers_2018 )
    else
      call edge_values_implicit_h4( GV%ke, hTmp, tmp, ppol_E, h_neglect=GV%H_subroundoff, &
                                  answers_2018=CS%answers_2018 )
    endif
    call PPM_reconstruction( GV%ke, hTmp, tmp, ppol_E, ppol_coefs, h_neglect, &
                                  answers_2018=CS%answers_2018 )
    if (bdry_extrap) &
      call PPM_boundary_extrapolation(GV%ke, hTmp, tmp, ppol_E, ppol_coefs, h_neglect )

    do k = 1,GV%ke
      T_t(i,j,k) = ppol_E(k,1)
      T_b(i,j,k) = ppol_E(k,2)
    enddo

  enddo ; enddo

end subroutine TS_PPM_edge_values


!> Initializes regridding for the main ALE algorithm
subroutine ALE_initRegridding(GV, US, max_depth, param_file, mdl, regridCS)
  type(verticalGrid_type), intent(in)  :: GV         !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)  :: US         !< A dimensional unit scaling type
  real,                    intent(in)  :: max_depth  !< The maximum depth of the ocean [Z ~> m].
  type(param_file_type),   intent(in)  :: param_file !< parameter file
  character(len=*),        intent(in)  :: mdl        !< Name of calling module
  type(regridding_CS),     intent(out) :: regridCS   !< Regridding parameters and work arrays
  ! Local variables
  character(len=30) :: coord_mode

  call get_param(param_file, mdl, "REGRIDDING_COORDINATE_MODE", coord_mode, &
                 "Coordinate mode for vertical regridding. "//&
                 "Choose among the following possibilities: "//&
                 trim(regriddingCoordinateModeDoc), &
                 default=DEFAULT_COORDINATE_MODE, fail_if_missing=.true.)

  call initialize_regridding(regridCS, GV, US, max_depth, param_file, mdl, coord_mode, '', '')

end subroutine ALE_initRegridding

!> Query the target coordinate interfaces positions
function ALE_getCoordinate( CS )
  type(ALE_CS), pointer    :: CS                  !< module control structure

  real, dimension(CS%nk+1) :: ALE_getCoordinate
  ALE_getCoordinate(:) = getCoordinateInterfaces( CS%regridCS, undo_scaling=.true. )

end function ALE_getCoordinate


!> Query the target coordinate units
function ALE_getCoordinateUnits( CS )
  type(ALE_CS), pointer :: CS   !< module control structure

  character(len=20)     :: ALE_getCoordinateUnits

  ALE_getCoordinateUnits = getCoordinateUnits( CS%regridCS )

end function ALE_getCoordinateUnits


!> Returns true if initial conditions should be regridded and remapped
logical function ALE_remap_init_conds( CS )
  type(ALE_CS), pointer :: CS   !< module control structure

  ALE_remap_init_conds = .false.
  if (associated(CS)) ALE_remap_init_conds = CS%remap_after_initialization
end function ALE_remap_init_conds

!> Updates the weights for time filtering the new grid generated in regridding
subroutine ALE_update_regrid_weights( dt, CS )
  real,         intent(in) :: dt !< Time-step used between ALE calls [T ~> s]
  type(ALE_CS), pointer    :: CS !< ALE control structure
  ! Local variables
  real :: w  ! An implicit weighting estimate.

  if (associated(CS)) then
    w = 0.0
    if (CS%regrid_time_scale > 0.0) then
      w = CS%regrid_time_scale / (CS%regrid_time_scale + dt)
    endif
    call set_regrid_params(CS%regridCS, old_grid_weight=w)
  endif

end subroutine ALE_update_regrid_weights

!> Update the vertical grid type with ALE information.
!! This subroutine sets information in the verticalGrid_type to be
!! consistent with the use of ALE mode.
subroutine ALE_updateVerticalGridType(CS, GV)
  type(ALE_CS),            pointer :: CS  !< ALE control structure
  type(verticalGrid_type), pointer :: GV  !< vertical grid information

  integer :: nk

  nk = GV%ke
  GV%sInterface(1:nk+1) = getCoordinateInterfaces( CS%regridCS, undo_scaling=.true. )
  GV%sLayer(1:nk) = 0.5*( GV%sInterface(1:nk) + GV%sInterface(2:nk+1) )
  GV%zAxisUnits = getCoordinateUnits( CS%regridCS )
  GV%zAxisLongName = getCoordinateShortName( CS%regridCS )
  GV%direction = -1 ! Because of ferret in z* mode. Need method to set
                    ! as function of coordinae mode.

end subroutine ALE_updateVerticalGridType


!> Write the vertical coordinate information into a file.
!! This subroutine writes out a file containing any available data related
!! to the vertical grid used by the MOM ocean model when in ALE mode.
subroutine ALE_writeCoordinateFile( CS, GV, directory )
  type(ALE_CS),            pointer     :: CS         !< module control structure
  type(verticalGrid_type), intent(in)  :: GV         !< ocean vertical grid structure
  character(len=*),        intent(in)  :: directory  !< directory for writing grid info

  character(len=240) :: filepath
  type(vardesc)      :: vars(2)
  type(fieldtype)    :: fields(2)
  integer            :: unit
  real               :: ds(GV%ke), dsi(GV%ke+1)

  filepath    = trim(directory) // trim("Vertical_coordinate")
  ds(:)       = getCoordinateResolution( CS%regridCS, undo_scaling=.true. )
  dsi(1)      = 0.5*ds(1)
  dsi(2:GV%ke) = 0.5*( ds(1:GV%ke-1) + ds(2:GV%ke) )
  dsi(GV%ke+1) = 0.5*ds(GV%ke)

  vars(1) = var_desc('ds', getCoordinateUnits( CS%regridCS ), &
                    'Layer Coordinate Thickness','1','L','1')
  vars(2) = var_desc('ds_interface', getCoordinateUnits( CS%regridCS ), &
                    'Layer Center Coordinate Separation','1','i','1')

  call create_file(unit, trim(filepath), vars, 2, fields, SINGLE_FILE, GV=GV)
  call write_field(unit, fields(1), ds)
  call write_field(unit, fields(2), dsi)
  call close_file(unit)

end subroutine ALE_writeCoordinateFile


!> Set h to coordinate values for fixed coordinate systems
subroutine ALE_initThicknessToCoord( CS, G, GV, h )
  type(ALE_CS), intent(inout)                            :: CS  !< module control structure
  type(ocean_grid_type), intent(in)                      :: G   !< module grid structure
  type(verticalGrid_type), intent(in)                    :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: h   !< layer thickness [H ~> m or kg m-2]

  ! Local variables
  integer :: i, j, k

  do j = G%jsd,G%jed ; do i = G%isd,G%ied
    h(i,j,:) = GV%Z_to_H * getStaticThickness( CS%regridCS, 0., G%bathyT(i,j) )
  enddo ; enddo

end subroutine ALE_initThicknessToCoord

end module MOM_ALE
