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
use MOM_diag_mediator,    only : time_type, diag_update_remap_grids, query_averaging_enabled
use MOM_domains,          only : create_group_pass, do_group_pass, group_pass_type
use MOM_error_handler,    only : MOM_error, FATAL, WARNING
use MOM_error_handler,    only : callTree_showQuery
use MOM_error_handler,    only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_hybgen_unmix,     only : hybgen_unmix, init_hybgen_unmix, end_hybgen_unmix, hybgen_unmix_CS
use MOM_hybgen_regrid,    only : hybgen_regrid_CS
use MOM_file_parser,      only : get_param, param_file_type, log_param
use MOM_interface_heights,only : find_eta, calc_derived_thermo
use MOM_open_boundary,    only : ocean_OBC_type, OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary,    only : OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_regridding,       only : initialize_regridding, regridding_main, end_regridding
use MOM_regridding,       only : uniformResolution
use MOM_regridding,       only : inflate_vanished_layers_old
use MOM_regridding,       only : regridding_preadjust_reqs, convective_adjustment
use MOM_regridding,       only : set_target_densities_from_GV, set_target_densities
use MOM_regridding,       only : regriddingCoordinateModeDoc, DEFAULT_COORDINATE_MODE
use MOM_regridding,       only : regriddingInterpSchemeDoc, regriddingDefaultInterpScheme
use MOM_regridding,       only : regriddingDefaultBoundaryExtrapolation
use MOM_regridding,       only : regriddingDefaultMinThickness
use MOM_regridding,       only : regridding_CS, set_regrid_params, write_regrid_file
use MOM_regridding,       only : getCoordinateInterfaces
use MOM_regridding,       only : getCoordinateUnits, getCoordinateShortName
use MOM_regridding,       only : getStaticThickness
use MOM_remapping,        only : initialize_remapping, end_remapping
use MOM_remapping,        only : remapping_core_h, remapping_core_w
use MOM_remapping,        only : remappingSchemesDoc, remappingDefaultScheme
use MOM_remapping,        only : interpolate_column, reintegrate_column
use MOM_remapping,        only : remapping_CS, dzFromH1H2, remapping_set_param
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
  logical :: partial_cell_vel_remap !< If true, use partial cell thicknesses at velocity points
                                    !! that are masked out where they extend below the shallower
                                    !! of the neighboring bathymetry for remapping velocity.

  real :: regrid_time_scale !< The time-scale used in blending between the current (old) grid
                            !! and the target (new) grid [T ~> s]

  type(regridding_CS) :: regridCS !< Regridding parameters and work arrays
  type(remapping_CS)  :: remapCS  !< Remapping parameters and work arrays
  type(remapping_CS)  :: vel_remapCS  !< Remapping parameters for velocities and work arrays

  type(hybgen_unmix_CS), pointer :: hybgen_unmixCS => NULL() !< Parameters for hybgen remapping

  logical :: use_hybgen_unmix   !< If true, use the hybgen unmixing code before regridding
  logical :: do_conv_adj        !< If true, do convective adjustment before regridding

  integer :: nk             !< Used only for queries, not directly by this module
  real :: BBL_h_vel_mask    !< The thickness of a bottom boundary layer within which velocities in
                            !! thin layers are zeroed out after remapping, following practice with
                            !! Hybgen remapping, or a negative value to avoid such filtering
                            !! altogether, in [H ~> m or kg m-2].
  real :: h_vel_mask        !< A thickness at velocity points below which near-bottom layers are
                            !! zeroed out after remapping, following the practice with Hybgen
                            !! remapping, or a negative value to avoid such filtering altogether,
                            !! in [H ~> m or kg m-2].

  logical :: remap_after_initialization !< Indicates whether to regrid/remap after initializing the state.

  integer :: answer_date    !< The vintage of the expressions and order of arithmetic to use for
                            !! remapping. Values below 20190101 result in the use of older, less
                            !! accurate expressions that were in use at the end of 2018.  Higher
                            !! values result in the use of more robust and accurate forms of
                            !! mathematically equivalent expressions.

  logical :: conserve_ke    !< Apply a correction to the baroclinic velocity after remapping to
                            !! conserve KE.

  logical :: debug   !< If true, write verbose checksums for debugging purposes.
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
  integer :: id_remap_delta_integ_u2 = -1  !< Change in depth-integrated rho0*u**2/2
  integer :: id_remap_delta_integ_v2 = -1  !< Change in depth-integrated rho0*v**2/2

end type

! Publicly available functions
public ALE_init
public ALE_end
public ALE_regrid
public ALE_offline_inputs
public ALE_regrid_accelerated
public ALE_remap_scalar
public ALE_remap_tracers
public ALE_remap_velocities
public ALE_remap_set_h_vel, ALE_remap_set_h_vel_via_dz
public ALE_remap_interface_vals
public ALE_remap_vertex_vals
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
public pre_ALE_diagnostics
public pre_ALE_adjustments
public ALE_remap_init_conds
public ALE_register_diags
public ALE_set_extrap_boundaries

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
  character(len=40) :: mdl = "MOM_ALE" ! This module's name.
  character(len=80) :: string, vel_string ! Temporary strings
  real              :: filter_shallow_depth, filter_deep_depth ! Depth ranges of filtering [H ~> m or kg m-2]
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  logical           :: check_reconstruction
  logical           :: check_remapping
  logical           :: force_bounds_in_subcell
  logical           :: local_logical
  logical           :: remap_boundary_extrap
  logical           :: init_boundary_extrap
  logical           :: om4_remap_via_sub_cells
  type(hybgen_regrid_CS), pointer :: hybgen_regridCS => NULL() ! Control structure for hybgen regridding
                                                         ! for sharing parameters.
  real :: h_neglect, h_neglect_edge ! small thicknesses [H ~> m or kg m-2]

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
  call regridding_preadjust_reqs(CS%regridCS, CS%do_conv_adj, CS%use_hybgen_unmix, hybgen_CS=hybgen_regridCS)

  ! Initialize and configure remapping that is orchestrated by ALE.
  call get_param(param_file, mdl, "REMAPPING_SCHEME", string, &
                 "This sets the reconstruction scheme used "//&
                 "for vertical remapping for all variables. "//&
                 "It can be one of the following schemes: \n"//&
                 trim(remappingSchemesDoc), default=remappingDefaultScheme)
  call get_param(param_file, mdl, "VELOCITY_REMAPPING_SCHEME", vel_string, &
                 "This sets the reconstruction scheme used for vertical remapping "//&
                 "of velocities. By default it is the same as REMAPPING_SCHEME. "//&
                 "It can be one of the following schemes: \n"//&
                 trim(remappingSchemesDoc), default=trim(string))
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
  call get_param(param_file, mdl, "INIT_BOUNDARY_EXTRAP", init_boundary_extrap, &
                 "If true, values at the interfaces of boundary cells are "//&
                 "extrapolated instead of piecewise constant during initialization."//&
                 "Defaults to REMAP_BOUNDARY_EXTRAP.", default=remap_boundary_extrap)
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 "This selects the remapping algorithm used in OM4 that does not use "//&
                 "the full reconstruction for the top- and lower-most sub-layers, but instead "//&
                 "assumes they are always vanished (untrue) and so just uses their edge values. "//&
                 "We recommend setting this option to false.", default=.true.)
  call get_param(param_file, mdl, "REMAPPING_ANSWER_DATE", CS%answer_date, &
                 "The vintage of the expressions and order of arithmetic to use for remapping.  "//&
                 "Values below 20190101 result in the use of older, less accurate expressions "//&
                 "that were in use at the end of 2018.  Higher values result in the use of more "//&
                 "robust and accurate forms of mathematically equivalent expressions.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) CS%answer_date = max(CS%answer_date, 20230701)

  if (CS%answer_date >= 20190101) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H * 1.0e-30 ; h_neglect_edge = GV%m_to_H * 1.0e-10
  else
    h_neglect = GV%kg_m2_to_H * 1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H * 1.0e-10
  endif

  call initialize_remapping( CS%remapCS, string, nk=GV%ke, &
                             boundary_extrapolation=init_boundary_extrap, &
                             check_reconstruction=check_reconstruction, &
                             check_remapping=check_remapping, &
                             force_bounds_in_subcell=force_bounds_in_subcell, &
                             om4_remap_via_sub_cells=om4_remap_via_sub_cells, &
                             answer_date=CS%answer_date, &
                             h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping( CS%vel_remapCS, vel_string, nk=GV%ke, &
                             boundary_extrapolation=init_boundary_extrap, &
                             check_reconstruction=check_reconstruction, &
                             check_remapping=check_remapping, &
                             force_bounds_in_subcell=force_bounds_in_subcell, &
                             om4_remap_via_sub_cells=om4_remap_via_sub_cells, &
                             answer_date=CS%answer_date, &
                             h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)

  call get_param(param_file, mdl, "PARTIAL_CELL_VELOCITY_REMAP", CS%partial_cell_vel_remap, &
                 "If true, use partial cell thicknesses at velocity points that are masked out "//&
                 "where they extend below the shallower of the neighboring bathymetry for "//&
                 "remapping velocity.", default=.false.)

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
                 "If true, the regridding integrates upwards from the bottom for "//&
                 "interface positions, much as the main model does. If false "//&
                 "regridding integrates downward, consistent with the remapping code.", &
                 default=.true., do_not_log=.true.)
  call set_regrid_params(CS%regridCS, integrate_downward_for_e=.not.local_logical)

  call get_param(param_file, mdl, "REMAP_VEL_MASK_BBL_THICK", CS%BBL_h_vel_mask, &
                 "A thickness of a bottom boundary layer below which velocities in thin layers "//&
                 "are zeroed out after remapping, following practice with Hybgen remapping, "//&
                 "or a negative value to avoid such filtering altogether.", &
                 default=-0.001, units="m", scale=GV%m_to_H)
  call get_param(param_file, mdl, "REMAP_VEL_MASK_H_THIN", CS%h_vel_mask, &
                 "A thickness at velocity points below which near-bottom layers are zeroed out "//&
                 "after remapping, following practice with Hybgen remapping, or a negative value "//&
                 "to avoid such filtering altogether.", &
                 default=1.0e-6, units="m", scale=GV%m_to_H, do_not_log=(CS%BBL_h_vel_mask<=0.0))

  if (CS%use_hybgen_unmix) &
    call init_hybgen_unmix(CS%hybgen_unmixCS, GV, US, param_file, hybgen_regridCS)

  call get_param(param_file, mdl, "REMAP_VEL_CONSERVE_KE", CS%conserve_ke, &
                 "If true, a correction is applied to the baroclinic component of velocity "//&
                 "after remapping so that total KE is conserved. KE may not be conserved "//&
                 "when (CS%BBL_h_vel_mask > 0.0) .and. (CS%h_vel_mask > 0.0)", &
                 default=.false.)
  call get_param(param_file, "MOM", "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)

  ! Keep a record of values for subsequent queries
  CS%nk = GV%ke

  if (CS%show_call_tree) call callTree_leave("ALE_init()")
end subroutine ALE_init

!> Sets the boundary extrapolation set for the remapping type.
subroutine ALE_set_extrap_boundaries( param_file, CS)
  type(param_file_type),   intent(in) :: param_file !< Parameter file
  type(ALE_CS),            pointer    :: CS         !< Module control structure

  logical :: remap_boundary_extrap
  call get_param(param_file, "MOM_ALE", "REMAP_BOUNDARY_EXTRAP", remap_boundary_extrap, &
                 "If true, values at the interfaces of boundary cells are "//&
                 "extrapolated instead of piecewise constant", default=.false.)
  call remapping_set_param(CS%remapCS, boundary_extrapolation=remap_boundary_extrap)
end subroutine ALE_set_extrap_boundaries

!> Sets the remapping algorithm to that of OM4
!!
!! The remapping aglorithm used in OM4 made poor assumptions about the reconstructions
!! in the top/bottom layers, namely that they were always vanished and could be
!! represented solely by their upper/lower edge value respectively.
!! Passing .false. here uses the full reconstruction of those top and bottom layers
!! and properly sample those layers.
subroutine ALE_set_OM4_remap_algorithm( CS, om4_remap_via_sub_cells )
  type(ALE_CS), pointer :: CS !< Module control structure
  logical, intent(in)   :: om4_remap_via_sub_cells !< If true, use OM4 remapping algorithm

  call remapping_set_param(CS%remapCS, om4_remap_via_sub_cells=om4_remap_via_sub_cells )

end subroutine ALE_set_OM4_remap_algorithm

!> Initialize diagnostics for the ALE module.
subroutine ALE_register_diags(Time, G, GV, US, diag, CS)
  type(time_type),target,     intent(in)  :: Time  !< Time structure
  type(ocean_grid_type),      intent(in)  :: G     !< Grid structure
  type(unit_scale_type),      intent(in)  :: US    !< A dimensional unit scaling type
  type(verticalGrid_type),    intent(in)  :: GV    !< Ocean vertical grid structure
  type(diag_ctrl), target,    intent(in)  :: diag  !< Diagnostics control structure
  type(ALE_CS), pointer                   :: CS    !< Module control structure

  ! Local variables
  character(len=48)  :: thickness_units

  CS%diag => diag
  thickness_units = get_thickness_units(GV)

  ! These diagnostics of the state variables before ALE are useful for
  ! debugging the ALE code.
  CS%id_u_preale = register_diag_field('ocean_model', 'u_preale', diag%axesCuL, Time, &
      'Zonal velocity before remapping', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_v_preale = register_diag_field('ocean_model', 'v_preale', diag%axesCvL, Time, &
      'Meridional velocity before remapping', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_h_preale = register_diag_field('ocean_model', 'h_preale', diag%axesTL, Time, &
      'Layer Thickness before remapping', thickness_units, conversion=GV%H_to_MKS, &
      v_extensive=.true.)
  CS%id_T_preale = register_diag_field('ocean_model', 'T_preale', diag%axesTL, Time, &
      'Temperature before remapping', 'degC', conversion=US%C_to_degC)
  CS%id_S_preale = register_diag_field('ocean_model', 'S_preale', diag%axesTL, Time, &
      'Salinity before remapping', 'PSU', conversion=US%S_to_ppt)
  CS%id_e_preale = register_diag_field('ocean_model', 'e_preale', diag%axesTi, Time, &
      'Interface Heights before remapping', 'm', conversion=US%Z_to_m)

  CS%id_dzRegrid = register_diag_field('ocean_model', 'dzRegrid', diag%axesTi, Time, &
      'Change in interface height due to ALE regridding', 'm', conversion=GV%H_to_m)
  CS%id_vert_remap_h = register_diag_field('ocean_model', 'vert_remap_h', diag%axestl, Time, &
      'layer thicknesses after ALE regridding and remapping', &
      thickness_units, conversion=GV%H_to_MKS, v_extensive=.true.)
  CS%id_vert_remap_h_tendency = register_diag_field('ocean_model', &
      'vert_remap_h_tendency', diag%axestl, Time, &
      'Layer thicknesses tendency due to ALE regridding and remapping', &
      trim(thickness_units)//" s-1", conversion=GV%H_to_MKS*US%s_to_T, v_extensive=.true.)
  CS%id_remap_delta_integ_u2 = register_diag_field('ocean_model', 'ale_u2', diag%axesCu1, Time, &
      'Rate of change in half rho0 times depth integral of squared zonal'//&
      ' velocity by remapping. If REMAP_VEL_CONSERVE_KE is .true. then '//&
      ' this measures the change before the KE-conserving correction is applied.', &
      'W m-2', conversion=GV%H_to_kg_m2 * US%L_T_to_m_s**2 * US%s_to_T)
  CS%id_remap_delta_integ_v2 = register_diag_field('ocean_model', 'ale_v2', diag%axesCv1, Time, &
      'Rate of change in half rho0 times depth integral of squared meridional'//&
      ' velocity by remapping. If REMAP_VEL_CONSERVE_KE is .true. then '//&
      ' this measures the change before the KE-conserving correction is applied.', &
      'W m-2', conversion=GV%H_to_kg_m2 * US%L_T_to_m_s**2 * US%s_to_T)

end subroutine ALE_register_diags

!> Crudely adjust (initial) grid for integrity.
!! This routine is typically called (from initialize_MOM in file MOM.F90)
!! before the main time integration loop to initialize the regridding stuff.
!! We read the MOM_input file to register the values of different
!! regridding/remapping parameters.
subroutine adjustGridForIntegrity( CS, G, GV, h )
  type(ALE_CS),                              intent(in)    :: CS  !< Regridding parameters and options
  type(ocean_grid_type),                     intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                   intent(in)    :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h   !< Current 3D grid thickness that
                                                                  !! are to be adjusted [H ~> m or kg m-2]
  call inflate_vanished_layers_old( CS%regridCS, G, GV, h(:,:,:) )

end subroutine adjustGridForIntegrity


!> End of regridding (memory deallocation).
!! This routine is typically called (from MOM_end in file MOM.F90)
!! after the main time integration loop to deallocate the regridding stuff.
subroutine ALE_end(CS)
  type(ALE_CS), pointer :: CS  !< module control structure

  ! Deallocate memory used for the regridding
  call end_remapping( CS%remapCS )

  if (CS%use_hybgen_unmix) call end_hybgen_unmix( CS%hybgen_unmixCS )
  call end_regridding( CS%regridCS )

  deallocate(CS)

end subroutine ALE_end

!> Save any diagnostics of the state before ALE remapping.  These diagnostics are
!! mostly used for debugging.
subroutine pre_ALE_diagnostics(G, GV, US, h, u, v, tv, CS)
  type(ocean_grid_type),                      intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                    intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h   !< Current 3D grid obtained after the
                                                                   !! last time step [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u   !< Zonal velocity field [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v   !< Meridional velocity field [L T-1 ~> m s-1]
  type(thermo_var_ptrs),                      intent(inout) :: tv  !< Thermodynamic variable structure
  type(ALE_CS),                               pointer       :: CS  !< Regridding parameters and options

  ! Local variables
  real :: eta_preale(SZI_(G),SZJ_(G),SZK_(GV)+1)  ! Interface heights before remapping [Z ~> m]

  if (CS%id_u_preale > 0) call post_data(CS%id_u_preale, u,    CS%diag)
  if (CS%id_v_preale > 0) call post_data(CS%id_v_preale, v,    CS%diag)
  if (CS%id_h_preale > 0) call post_data(CS%id_h_preale, h,    CS%diag)
  if (CS%id_T_preale > 0) call post_data(CS%id_T_preale, tv%T, CS%diag)
  if (CS%id_S_preale > 0) call post_data(CS%id_S_preale, tv%S, CS%diag)
  if (CS%id_e_preale > 0) then
    call find_eta(h, tv, G, GV, US, eta_preale, dZref=G%Z_ref)
    call post_data(CS%id_e_preale, eta_preale, CS%diag)
  endif

end subroutine pre_ALE_diagnostics


!> Potentially do some preparatory work, such as convective adjustment, to clean up the model
!! state before regridding.
subroutine pre_ALE_adjustments(G, GV, US, h, tv, Reg, CS, u, v)
  type(ocean_grid_type),                      intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                    intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h   !< Current 3D grid obtained after the
                                                                   !! last time step [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv  !< Thermodynamic variable structure
  type(tracer_registry_type),                 pointer       :: Reg !< Tracer registry structure
  type(ALE_CS),                               pointer       :: CS  !< Regridding parameters and options
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                    optional, intent(inout) :: u   !< Zonal velocity field [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    optional, intent(inout) :: v   !< Meridional velocity field [L T-1 ~> m s-1]

  integer :: ntr

  ! Do column-wise convective adjustment.
  ! Tracers and velocities should probably also undergo consistent adjustments.
  if (CS%do_conv_adj) call convective_adjustment(G, GV, h, tv)

  if (CS%use_hybgen_unmix) then
    ntr = 0 ; if (associated(Reg)) ntr = Reg%ntr
    call hybgen_unmix(G, GV, US, CS%hybgen_unmixCS, tv, Reg, ntr, h)
  endif

end subroutine pre_ALE_adjustments

!> Takes care of building a new grid. The creation of the new grid can be based on z coordinates,
!! target interface densities, sigma coordinates or any arbitrary coordinate system.
subroutine ALE_regrid( G, GV, US, h, h_new, dzRegrid, tv, CS, frac_shelf_h, PCM_cell)
  type(ocean_grid_type),                      intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                    intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h   !< Layer thicknesses in 3D grid before
                                                                   !! regridding [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(out)   :: h_new !< Layer thicknesses in 3D grid after
                                                                   !! regridding [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(out)  :: dzRegrid !< The change in grid interface positions
                                                                   !! due to regridding, in the same units as
                                                                   !! thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv  !< Thermodynamic variable structure
  type(ALE_CS),                               pointer       :: CS  !< Regridding parameters and options
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)    :: frac_shelf_h !< Fractional ice shelf coverage [nondim]
  logical, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    optional, intent(out)   :: PCM_cell !< If true, use PCM remapping in a cell.

  ! Local variables
  logical :: showCallTree

  showCallTree = callTree_showQuery()

  if (showCallTree) call callTree_enter("ALE_regrid(), MOM_ALE.F90")

  ! Build the new grid and store it in h_new. The old grid is retained as h.
  ! Both are needed for the subsequent remapping of variables.
  dzRegrid(:,:,:) = 0.0
  call regridding_main( CS%remapCS, CS%regridCS, G, GV, US, h, tv, h_new, dzRegrid, &
                        frac_shelf_h=frac_shelf_h, PCM_cell=PCM_cell)

  if (CS%id_dzRegrid>0) then ; if (query_averaging_enabled(CS%diag)) then
    call post_data(CS%id_dzRegrid, dzRegrid, CS%diag, alt_h=h_new)
  endif ; endif

  if (showCallTree) call callTree_leave("ALE_regrid()")

end subroutine ALE_regrid

!> Regrid/remap stored fields used for offline tracer integrations. These input fields are assumed to have
!! the same layer thicknesses at the end of the last offline interval (which should be a Zstar grid). This
!! routine builds a grid on the runtime specified vertical coordinate
subroutine ALE_offline_inputs(CS, G, GV, US, h, tv, Reg, uhtr, vhtr, Kd, debug, OBC)
  type(ALE_CS),                                 pointer       :: CS    !< Regridding parameters and options
  type(ocean_grid_type),                        intent(in   ) :: G     !< Ocean grid informations
  type(verticalGrid_type),                      intent(in   ) :: GV    !< Ocean vertical grid structure
  type(unit_scale_type),                        intent(in   ) :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(inout) :: h     !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                        intent(inout) :: tv    !< Thermodynamic variable structure
  type(tracer_registry_type),                   pointer       :: Reg   !< Tracer registry structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),   intent(inout) :: uhtr  !< Zonal mass fluxes [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),   intent(inout) :: vhtr  !< Meridional mass fluxes [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(inout) :: Kd    !< Input diffusivities
                                                                       !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  logical,                                      intent(in   ) :: debug !< If true, then turn checksums
  type(ocean_OBC_type),                         pointer       :: OBC   !< Open boundary structure
  ! Local variables
  integer :: nk, i, j, k, isc, iec, jsc, jec
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV))   :: h_new    ! Layer thicknesses after regridding [H ~> m or kg m-2]
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)+1) :: dzRegrid ! The change in grid interface positions [H ~> m or kg m-2]
  real, dimension(SZK_(GV)) :: h_src   ! Source grid thicknesses at velocity points [H ~> m or kg m-2]
  real, dimension(SZK_(GV)) :: h_dest  ! Destination grid thicknesses at velocity points [H ~> m or kg m-2]
  real, dimension(SZK_(GV)) :: temp_vec ! Transports on the destination grid [H L2 ~> m3 or kg]

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nk = GV%ke
  dzRegrid(:,:,:) = 0.0
  h_new(:,:,:) = 0.0

  if (debug) call MOM_tracer_chkinv("Before ALE_offline_inputs", G, GV, h, Reg%Tr, Reg%ntr)

  ! Build new grid from the Zstar state onto the requested vertical coordinate. The new grid is stored
  ! in h_new. The old grid is h. Both are needed for the subsequent remapping of variables. Convective
  ! adjustment right now is not used because it is unclear what to do with vanished layers
  call regridding_main( CS%remapCS, CS%regridCS, G, GV, US, h, tv, h_new, dzRegrid)
  if (CS%show_call_tree) call callTree_waypoint("new grid generated (ALE_offline_inputs)")

  ! Remap all variables from old grid h onto new grid h_new
  call ALE_remap_tracers(CS, G, GV, h, h_new, Reg, debug=CS%show_call_tree)
  if (allocated(tv%SpV_avg)) tv%valid_SpV_halo = -1   ! Record that SpV_avg is no longer valid.
  if (CS%show_call_tree) call callTree_waypoint("state remapped (ALE_inputs)")

  ! Reintegrate mass transports from Zstar to the offline vertical coordinate
  do j=jsc,jec ; do i=G%iscB,G%iecB
    if (G%mask2dCu(i,j)>0.) then
      h_src(:) = 0.5 * (h(i,j,:) + h(i+1,j,:))
      h_dest(:) = 0.5 * (h_new(i,j,:) + h_new(i+1,j,:))
      call reintegrate_column(nk, h_src, uhtr(I,j,:), nk, h_dest, temp_vec)
      uhtr(I,j,:) = temp_vec
    endif
  enddo ; enddo
  do j=G%jscB,G%jecB ; do i=isc,iec
    if (G%mask2dCv(i,j)>0.) then
      h_src(:) = 0.5 * (h(i,j,:) + h(i,j+1,:))
      h_dest(:) = 0.5 * (h_new(i,j,:) + h_new(i,j+1,:))
      call reintegrate_column(nk, h_src, vhtr(I,j,:), nk, h_dest, temp_vec)
      vhtr(I,j,:) = temp_vec
    endif
  enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    if (G%mask2dT(i,j)>0.) then
      if (check_column_integrals(nk, h_src, nk, h_dest)) then
        call MOM_error(FATAL, "ALE_offline_inputs: Kd interpolation columns do not match")
      endif
      call interpolate_column(nk, h(i,j,:), Kd(i,j,:), nk, h_new(i,j,:), Kd(i,j,:), .true.)
    endif
  enddo ; enddo

  call ALE_remap_scalar(CS%remapCS, G, GV, nk, h, tv%T, h_new, tv%T)
  call ALE_remap_scalar(CS%remapCS, G, GV, nk, h, tv%S, h_new, tv%S)

  if (debug) call MOM_tracer_chkinv("After ALE_offline_inputs", G, GV, h_new, Reg%Tr, Reg%ntr)

  ! Copy over the new layer thicknesses
  do k = 1,nk  ; do j = jsc-1,jec+1 ; do i = isc-1,iec+1
    h(i,j,k) = h_new(i,j,k)
  enddo ; enddo ; enddo

  if (allocated(tv%SpV_avg)) tv%valid_SpV_halo = -1   ! Record that SpV_avg is no longer valid.

  if (CS%show_call_tree) call callTree_leave("ALE_offline_inputs()")
end subroutine ALE_offline_inputs


!> For a state-based coordinate, accelerate the process of regridding by
!! repeatedly applying the grid calculation algorithm
subroutine ALE_regrid_accelerated(CS, G, GV, US, h, tv, n_itt, u, v, OBC, Reg, dt, dzRegrid, initial)
  type(ALE_CS),            pointer       :: CS     !< ALE control structure
  type(ocean_grid_type),   intent(inout) :: G      !< Ocean grid
  type(verticalGrid_type), intent(in)    :: GV     !< Vertical grid
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h      !< Original thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv     !< Thermo vars (T/S/EOS)
  integer,                 intent(in)    :: n_itt  !< Number of times to regrid
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v      !< Meridional velocity [L T-1 ~> m s-1]
  type(ocean_OBC_type),    pointer       :: OBC    !< Open boundary structure
  type(tracer_registry_type), &
                 optional, pointer       :: Reg    !< Tracer registry to remap onto new grid
  real,          optional, intent(in)    :: dt     !< Model timestep to provide a timescale for regridding [T ~> s]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                 optional, intent(inout) :: dzRegrid !< Final change in interface positions [H ~> m or kg m-2]
  logical,       optional, intent(in)    :: initial !< Whether we're being called from an initialization
                                                    !! routine (and expect diagnostics to work)

  ! Local variables
  integer :: i, j, itt, nz
  type(thermo_var_ptrs) :: tv_local ! local/intermediate temp/salt
  type(group_pass_type) :: pass_T_S_h ! group pass if the coordinate has a stencil
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))         :: h_loc  ! A working copy of layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))         :: h_orig ! The original layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target :: T      ! local temporary temperatures [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target :: S      ! local temporary salinities [S ~> ppt]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV))        :: h_old_u ! Source grid thickness at zonal
                                                               ! velocity points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV))        :: h_old_v ! Source grid thickness at meridional
                                                               ! velocity points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV))        :: h_new_u ! Destination grid thickness at zonal
                                                               ! velocity points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV))        :: h_new_v ! Destination grid thickness at meridional
                                                               ! velocity points [H ~> m or kg m-2]

  ! we have to keep track of the total dzInterface if for some reason
  ! we're using the old remapping algorithm for u/v
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: dzInterface ! Interface height changes within
                                                             ! an iteration [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: dzIntTotal  ! Cumulative interface position changes [H ~> m or kg m-2]

  nz = GV%ke

  ! initial total interface displacement due to successive regridding
  if (CS%remap_uv_using_old_alg) &
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

  do itt = 1, n_itt

    call do_group_pass(pass_T_S_h, G%domain)

    ! generate new grid
    if (CS%do_conv_adj) call convective_adjustment(G, GV, h_loc, tv_local)

    ! Update the layer specific volumes if necessary
    if (allocated(tv_local%SpV_avg)) call calc_derived_thermo(tv_local, h, G, GV, US, halo=1)

    call regridding_main(CS%remapCS, CS%regridCS, G, GV, US, h_loc, tv_local, h, dzInterface)
    if (CS%remap_uv_using_old_alg) &
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
  call ALE_remap_tracers(CS, G, GV, h_orig, h, Reg)

  call ALE_remap_set_h_vel(CS, G, GV, h_orig, h_old_u, h_old_v, OBC)
  if (CS%remap_uv_using_old_alg) then
    call ALE_remap_set_h_vel_via_dz(CS, G, GV, h, h_new_u, h_new_v, OBC, h_orig, dzIntTotal)
  else
    call ALE_remap_set_h_vel(CS, G, GV, h, h_new_u, h_new_v, OBC)
  endif

  call ALE_remap_velocities(CS, G, GV, h_old_u, h_old_v, h_new_u, h_new_v, u, v)

  ! save total dzregrid for diags if needed?
  if (present(dzRegrid)) dzRegrid(:,:,:) = dzIntTotal(:,:,:)

  if (allocated(tv%SpV_avg)) tv%valid_SpV_halo = -1   ! Record that SpV_avg is no longer valid.

end subroutine ALE_regrid_accelerated

!> This routine takes care of remapping all tracer variables between the old and the
!! new grids. This routine is called during initialization of the model at time=0, to
!! remap initial conditions to the model grid.  It is also called during a
!! time step to update the state.
subroutine ALE_remap_tracers(CS, G, GV, h_old, h_new, Reg, debug, dt, PCM_cell)
  type(ALE_CS),                              intent(in)    :: CS           !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G            !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV           !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_old        !< Thickness of source grid
                                                                           !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_new        !< Thickness of destination grid
                                                                           !! [H ~> m or kg m-2]
  type(tracer_registry_type),                pointer       :: Reg          !< Tracer registry structure
  logical,                         optional, intent(in)    :: debug  !< If true, show the call tree
  real,                            optional, intent(in)    :: dt     !< time step for diagnostics [T ~> s]
  logical, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                   optional, intent(in)    :: PCM_cell !< Use PCM remapping in cells where true

  ! Local variables
  real :: tr_column(GV%ke)  ! A column of updated tracer concentrations [CU ~> Conc]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: work_conc ! The rate of change of concentrations [Conc T-1 ~> Conc s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: work_cont ! The rate of change of cell-integrated tracer
                                                       ! content [Conc H T-1 ~> Conc m s-1 or Conc kg m-2 s-1] or
                                                       ! cell thickness [H T-1 ~> m s-1 or kg m-2 s-1]
  real, dimension(SZI_(G),SZJ_(G))          :: work_2d ! The rate of change of column-integrated tracer
                                                       ! content [Conc H T-1 ~> Conc m s-1 or Conc kg m-2 s-1]
  logical :: PCM(GV%ke) ! If true, do PCM remapping from a cell.
  real :: Idt           ! The inverse of the timestep [T-1 ~> s-1]
  real :: h1(GV%ke)     ! A column of source grid layer thicknesses [H ~> m or kg m-2]
  real :: h2(GV%ke)     ! A column of target grid layer thicknesses [H ~> m or kg m-2]
  logical :: show_call_tree
  type(tracer_type), pointer :: Tr => NULL()
  integer :: i, j, k, m, nz, ntr

  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug

  if (show_call_tree) call callTree_enter("ALE_remap_tracers(), MOM_ALE.F90")

  nz = GV%ke

  ntr = 0 ; if (associated(Reg)) ntr = Reg%ntr

  if (present(dt)) then
    Idt = 1.0/dt
    work_conc(:,:,:) = 0.0
    work_cont(:,:,:) = 0.0
  endif

  ! Remap all registered tracers, including temperature and salinity.
  if (ntr>0) then
    if (show_call_tree) call callTree_waypoint("remapping tracers (ALE_remap_tracers)")
    !$OMP parallel do default(shared) private(h1,h2,tr_column,Tr,PCM,work_conc,work_cont,work_2d)
    do m=1,ntr ! For each tracer
      Tr => Reg%Tr(m)
      do j = G%jsc,G%jec ; do i = G%isc,G%iec ; if (G%mask2dT(i,j)>0.) then
        ! Build the start and final grids
        h1(:) = h_old(i,j,:)
        h2(:) = h_new(i,j,:)
        if (present(PCM_cell)) then
          PCM(:) = PCM_cell(i,j,:)
          call remapping_core_h(CS%remapCS, nz, h1, Tr%t(i,j,:), nz, h2, tr_column, PCM_cell=PCM)
        else
          call remapping_core_h(CS%remapCS, nz, h1, Tr%t(i,j,:), nz, h2, tr_column)
        endif

        ! Possibly underflow any very tiny tracer concentrations to 0.  Note that this is not conservative!
        if (Tr%conc_underflow > 0.0) then ; do k=1,GV%ke
          if (abs(tr_column(k)) < Tr%conc_underflow) tr_column(k) = 0.0
        enddo ; endif

        ! Intermediate steps for tendency of tracer concentration and tracer content.
        if (present(dt)) then
          if (Tr%id_remap_conc > 0) then
            do k=1,GV%ke
              work_conc(i,j,k) = (tr_column(k) - Tr%t(i,j,k)) * Idt
            enddo
          endif
          if (Tr%id_remap_cont > 0 .or. Tr%id_remap_cont_2d > 0) then
            do k=1,GV%ke
              work_cont(i,j,k) = (tr_column(k)*h2(k) - Tr%t(i,j,k)*h1(k)) * Idt
            enddo
          endif
        endif

        ! update tracer concentration
        Tr%t(i,j,:) = tr_column(:)
      endif ; enddo ; enddo

      ! tendency diagnostics.
      if (present(dt)) then
        if (Tr%id_remap_conc > 0) then
          call post_data(Tr%id_remap_conc, work_conc, CS%diag)
        endif
        if (Tr%id_remap_cont > 0) then
          call post_data(Tr%id_remap_cont, work_cont, CS%diag)
        endif

        if (Tr%id_remap_cont_2d > 0) then
          do j = G%jsc,G%jec ; do i = G%isc,G%iec
            work_2d(i,j) = 0.0
            do k = 1,GV%ke
              work_2d(i,j) = work_2d(i,j) + work_cont(i,j,k)
            enddo
          enddo ; enddo
          call post_data(Tr%id_remap_cont_2d, work_2d, CS%diag)
        endif
      endif
    enddo ! m=1,ntr

  endif  ! endif for ntr > 0


  if (CS%id_vert_remap_h > 0) call post_data(CS%id_vert_remap_h, h_old, CS%diag)
  if ((CS%id_vert_remap_h_tendency > 0) .and. present(dt)) then
    do k = 1, nz ; do j = G%jsc,G%jec ; do i = G%isc,G%iec
      work_cont(i,j,k) = (h_new(i,j,k) - h_old(i,j,k))*Idt
    enddo ; enddo ; enddo
    call post_data(CS%id_vert_remap_h_tendency, work_cont, CS%diag)
  endif

  if (show_call_tree) call callTree_leave("ALE_remap_tracers(), MOM_ALE.F90")

end subroutine ALE_remap_tracers

!> This routine sets the thicknesses at velocity points used for vertical remapping.
subroutine ALE_remap_set_h_vel(CS, G, GV, h_new, h_u, h_v, OBC, debug)
  type(ALE_CS),                              intent(in)    :: CS      !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G       !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV      !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_new   !< Thickness at tracer points of the
                                                                      !! grid being interpolated to velocity
                                                                      !! points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(inout) :: h_u     !< Grid thickness at zonal velocity
                                                                      !! points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(inout) :: h_v     !< Grid thickness at meridional velocity
                                                                      !! points [H ~> m or kg m-2]
  type(ocean_OBC_type),                      pointer       :: OBC     !< Open boundary structure
  logical,                         optional, intent(in)    :: debug   !< If true, show the call tree

  ! Local variables
  logical :: show_call_tree
  integer :: i, j, k

  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug
  if (show_call_tree) call callTree_enter("ALE_remap_set_h_vel()")

  ! Build the u- and v-velocity grid thicknesses for remapping.

  !$OMP parallel do default(shared)
  do k=1,GV%ke ; do j=G%jsc,G%jec ; do I=G%IscB,G%IecB ; if (G%mask2dCu(I,j)>0.) then
    h_u(I,j,k) = 0.5*(h_new(i,j,k) + h_new(i+1,j,k))
  endif ; enddo ; enddo ; enddo
  !$OMP parallel do default(shared)
  do k=1,GV%ke ; do J=G%JscB,G%JecB ; do i=G%isc,G%iec ; if (G%mask2dCv(i,J)>0.) then
    h_v(i,J,k) = 0.5*(h_new(i,j,k) + h_new(i,j+1,k))
  endif ; enddo ; enddo ; enddo

  ! Mask out blocked portions of velocity cells.
  if (CS%partial_cell_vel_remap) call ALE_remap_set_h_vel_partial(CS, G, GV, h_new, h_u, h_v)

  ! Take open boundary conditions into account.
  if (associated(OBC)) call ALE_remap_set_h_vel_OBC(G, GV, h_new, h_u, h_v, OBC)

  if (show_call_tree) call callTree_leave("ALE_remap_set_h_vel()")

end subroutine ALE_remap_set_h_vel

!> This routine sets the thicknesses at velocity points used for vertical remapping using a
!! combination of the old grid and interface movements.
subroutine ALE_remap_set_h_vel_via_dz(CS, G, GV, h_new, h_u, h_v, OBC, h_old, dzInterface, debug)
  type(ALE_CS),                              intent(in)    :: CS           !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G            !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV           !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_new        !< Thickness at tracer points of the
                                                                           !! grid being interpolated to velocity
                                                                           !! points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(inout) :: h_u          !< Grid thickness at zonal velocity
                                                                           !! points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(inout) :: h_v          !< Grid thickness at meridional velocity
                                                                           !! points [H ~> m or kg m-2]
  type(ocean_OBC_type),                      pointer       :: OBC          !< Open boundary structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                             intent(in)    :: h_old        !< Thickness of source grid when generating
                                                                           !! the destination grid via the old
                                                                           !! algorithm [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                                             intent(in)    :: dzInterface  !< Change in interface position
                                                                           !! [H ~> m or kg m-2]
  logical,                         optional, intent(in)    :: debug        !< If true, show the call tree

  ! Local variables
  logical :: show_call_tree
  integer :: i, j, k

  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug
  if (show_call_tree) call callTree_enter("ALE_remap_set_h_vel()")

  ! Build the u- and v-velocity grid thicknesses for remapping using the old grid and interface movement.

  !$OMP parallel do default(shared)
  do k=1,GV%ke ; do j=G%jsc,G%jec ; do I=G%IscB,G%IecB ; if (G%mask2dCu(I,j)>0.) then
    h_u(I,j,k) = max( 0., 0.5*(h_old(i,j,k) + h_old(i+1,j,k)) + &
            0.5 * (( dzInterface(i,j,k) + dzInterface(i+1,j,k) ) - &
                   ( dzInterface(i,j,k+1) + dzInterface(i+1,j,k+1) )) )
  endif ; enddo ; enddo ; enddo

  !$OMP parallel do default(shared)
  do k=1,GV%ke ; do J=G%JscB,G%JecB ; do i=G%isc,G%iec ; if (G%mask2dCv(i,J)>0.) then
    h_v(i,J,k) = max( 0., 0.5*(h_old(i,j,k) + h_old(i,j+1,k)) + &
            0.5 * (( dzInterface(i,j,k) + dzInterface(i,j+1,k) ) - &
                   ( dzInterface(i,j,k+1) + dzInterface(i,j+1,k+1) )) )
  endif ; enddo ; enddo ; enddo

  ! Mask out blocked portions of velocity cells.
  if (CS%partial_cell_vel_remap) call ALE_remap_set_h_vel_partial(CS, G, GV, h_old, h_u, h_v)

  ! Take open boundary conditions into account.
  if (associated(OBC)) call ALE_remap_set_h_vel_OBC(G, GV, h_new, h_u, h_v, OBC)

  if (show_call_tree) call callTree_leave("ALE_remap_set_h_vel()")

end subroutine ALE_remap_set_h_vel_via_dz

!> Mask out the thicknesses at velocity points where they are below the minimum depth
!! at adjacent tracer points
subroutine ALE_remap_set_h_vel_partial(CS, G, GV, h_mask, h_u, h_v)
  type(ALE_CS),                              intent(in)    :: CS           !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G            !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV           !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_mask       !< Thickness at tracer points
                                                                           !! used to apply the partial
                                                                           !! cell masking [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(inout) :: h_u          !< Grid thickness at zonal velocity
                                                                           !! points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(inout) :: h_v          !< Grid thickness at meridional velocity
                                                                           !! points [H ~> m or kg m-2]
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: h_tot  ! The vertically summed thicknesses [H ~> m or kg m-2]
  real :: h_mask_vel ! A depth below which the thicknesses at a velocity point are masked out [H ~> m or kg m-2]
  integer :: i, j, k

  h_tot(:,:) = 0.0
  do k=1,GV%ke ; do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
    h_tot(i,j) = h_tot(i,j) + h_mask(i,j,k)
  enddo ; enddo ; enddo

  !$OMP parallel do default(shared) private(h_mask_vel)
  do j=G%jsc,G%jec ; do I=G%IscB,G%IecB ; if (G%mask2dCu(I,j)>0.) then
    h_mask_vel = min(h_tot(i,j), h_tot(i+1,j))
    call apply_partial_cell_mask(h_u(I,j,:), h_mask_vel)
  endif ; enddo ; enddo

  !$OMP parallel do default(shared) private(h_mask_vel)
  do J=G%JscB,G%JecB ; do i=G%isc,G%iec ; if (G%mask2dCv(i,J)>0.) then
    h_mask_vel = min(h_tot(i,j), h_tot(i,j+1))
    call apply_partial_cell_mask(h_v(i,J,:), h_mask_vel)
  endif ; enddo ; enddo

end subroutine ALE_remap_set_h_vel_partial

! Reset thicknesses at velocity points on open boundary condition segments
subroutine ALE_remap_set_h_vel_OBC(G, GV, h_new, h_u, h_v, OBC)
  type(ocean_grid_type),                     intent(in)    :: G            !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV           !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_new        !< Thickness at tracer points of the
                                                                           !! grid being interpolated to velocity
                                                                           !! points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(inout) :: h_u          !< Grid thickness at zonal velocity
                                                                           !! points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(inout) :: h_v          !< Grid thickness at meridional velocity
                                                                           !! points [H ~> m or kg m-2]
  type(ocean_OBC_type),                      pointer       :: OBC          !< Open boundary structure

  ! Local variables
  integer :: i, j, k, nz, is_OBC, ie_OBC, js_OBC, je_OBC

  if (.not.associated(OBC)) return

  nz = GV%ke

  ! Take open boundary conditions into account.
  if (OBC%u_E_OBCs_on_PE) then
    js_OBC = max(G%jsc,  OBC%js_u_E_obc) ; je_OBC = min(G%jec,  OBC%je_u_E_obc)
    Is_OBC = max(G%IscB, OBC%Is_u_E_obc) ; Ie_OBC = min(G%IecB, OBC%Ie_u_E_obc)
    !$OMP parallel do default(shared)
    do j=js_OBC,je_OBC ; do I=Is_OBC,Ie_OBC ; if (OBC%segnum_u(I,j) > 0) then !  OBC_DIRECTION_E
      do k=1,nz ; h_u(I,j,k) = h_new(i,j,k) ; enddo
    endif ; enddo ; enddo
  endif
  if (OBC%u_W_OBCs_on_PE) then
    js_OBC = max(G%jsc,  OBC%js_u_W_obc) ; je_OBC = min(G%jec,  OBC%je_u_W_obc)
    Is_OBC = max(G%IscB, OBC%Is_u_W_obc) ; Ie_OBC = min(G%IecB, OBC%Ie_u_W_obc)
    !$OMP parallel do default(shared)
    do j=js_OBC,je_OBC ; do I=Is_OBC,Ie_OBC ; if (OBC%segnum_u(I,j) < 0) then !  OBC_DIRECTION_W
      do k=1,nz ; h_u(I,j,k) = h_new(i+1,j,k) ; enddo
    endif ; enddo ; enddo
  endif

  if (OBC%v_N_OBCs_on_PE) then
    Js_OBC = max(G%JscB, OBC%Js_v_N_obc) ; Je_OBC = min(G%JecB, OBC%Je_v_N_obc)
    is_OBC = max(G%isc,  OBC%is_v_N_obc) ; ie_OBC = min(G%iec,  OBC%ie_v_N_obc)
    !$OMP parallel do default(shared)
    do J=Js_OBC,Je_OBC ; do i=is_OBC,ie_OBC ; if (OBC%segnum_v(i,J) > 0) then !  OBC_DIRECTION_N
      do k=1,nz ; h_v(i,J,k) = h_new(i,j,k) ; enddo
    endif ; enddo ; enddo
  endif
  if (OBC%v_S_OBCs_on_PE) then
    Js_OBC = max(G%JscB, OBC%Js_v_S_obc) ; Je_OBC = min(G%JecB, OBC%Je_v_S_obc)
    is_OBC = max(G%isc,  OBC%is_v_S_obc) ; ie_OBC = min(G%iec,  OBC%ie_v_S_obc)
    !$OMP parallel do default(shared)
    do J=Js_OBC,Je_OBC ; do i=is_OBC,ie_OBC ; if (OBC%segnum_v(i,J) < 0) then !  OBC_DIRECTION_S
      do k=1,nz ; h_v(i,J,k) = h_new(i,j+1,k) ; enddo
    endif ; enddo ; enddo
  endif

end subroutine ALE_remap_set_h_vel_OBC

!> This routine remaps velocity components between the old and the new grids,
!! with thicknesses at velocity points taken to be arithmetic averages of tracer thicknesses.
!! This routine may be called during initialization of the model at time=0, to
!! remap initial conditions to the model grid.  It is also called during a
!! time step to update the state.
subroutine ALE_remap_velocities(CS, G, GV, h_old_u, h_old_v, h_new_u, h_new_v, u, v, debug, &
                                dt, allow_preserve_variance)
  type(ALE_CS),                              intent(in)    :: CS        !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G         !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV        !< Ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(in)    :: h_old_u   !< Source grid thickness at zonal
                                                                        !! velocity points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(in)    :: h_old_v   !< Source grid thickness at meridional
                                                                        !! velocity points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(in)    :: h_new_u   !< Destination grid thickness at zonal
                                                                        !! velocity points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(in)    :: h_new_v   !< Destination grid thickness at meridional
                                                                        !! velocity points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                             intent(inout) :: u         !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                             intent(inout) :: v         !< Meridional velocity [L T-1 ~> m s-1]
  logical,                         optional, intent(in)    :: debug     !< If true, show the call tree
  real,                            optional, intent(in)    :: dt        !< time step for diagnostics [T ~> s]
  logical,                         optional, intent(in)    :: allow_preserve_variance !< If true, enables ke-conserving
                                                                                      !! correction

  ! Local variables
  real :: h_mask_vel ! A depth below which the thicknesses at a velocity point are masked out [H ~> m or kg m-2]
  real :: u_src(GV%ke)  ! A column of u-velocities on the source grid [L T-1 ~> m s-1]
  real :: u_tgt(GV%ke)  ! A column of u-velocities on the target grid [L T-1 ~> m s-1]
  real :: v_src(GV%ke)  ! A column of v-velocities on the source grid [L T-1 ~> m s-1]
  real :: v_tgt(GV%ke)  ! A column of v-velocities on the target grid [L T-1 ~> m s-1]
  real :: h1(GV%ke)     ! A column of source grid layer thicknesses [H ~> m or kg m-2]
  real :: h2(GV%ke)     ! A column of target grid layer thicknesses [H ~> m or kg m-2]
  real :: rescale_coef  ! Factor that scales the baroclinic velocity to conserve ke [nondim]
  real :: u_bt, v_bt    ! Depth-averaged velocity components [L T-1 ~> m s-1]
  real :: ke_c_src, ke_c_tgt ! \int [u_c or v_c]^2 dz on src and tgt grids [H L2 T-2 ~> m3 s-2]
  real, dimension(SZIB_(G),SZJ_(G)) :: du2h_tot  ! The rate of change of vertically integrated
                                                 ! 0.5 * rho0 *  u**2 [R Z L2 T-3 ~> W m-2]
  real, dimension(SZI_(G),SZJB_(G)) :: dv2h_tot  ! The rate of change of vertically integrated
                                                 ! 0.5 * rho0 *  v**2 [R Z L2 T-3 ~> W m-2]
  real :: u2h_tot, v2h_tot   ! The vertically integrated u**2 and v**2 [H L2 T-2 ~> m3 s-2 or kg s-2]
  real :: I_dt               ! 1 / dt [T-1 ~> s-1]
  logical :: variance_option ! Contains the value of allow_preserve_variance when present, else false
  logical :: show_call_tree
  integer :: i, j, k, nz

  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug
  if (show_call_tree) call callTree_enter("ALE_remap_velocities()")

  ! Setup related to KE conservation
  variance_option = .false.
  if (present(allow_preserve_variance)) variance_option=allow_preserve_variance
  if (present(dt)) I_dt = 1.0 / dt

  if (CS%id_remap_delta_integ_u2>0) du2h_tot(:,:) = 0.
  if (CS%id_remap_delta_integ_v2>0) dv2h_tot(:,:) = 0.

  if (((CS%id_remap_delta_integ_u2>0) .or. (CS%id_remap_delta_integ_v2>0)) .and. .not.present(dt))&
    call MOM_error(FATAL, "ALE KE diagnostics requires passing dt into ALE_remap_velocities")

  nz = GV%ke

  ! --- Remap u profiles from the source vertical grid onto the new target grid.

  !$OMP parallel do default(shared) private(h1,h2,u_src,h_mask_vel,u_tgt, &
  !$OMP                                     u_bt,ke_c_src,ke_c_tgt,rescale_coef, &
  !$OMP                                     u2h_tot,v2h_tot)
  do j=G%jsc,G%jec ; do I=G%IscB,G%IecB ; if (G%mask2dCu(I,j)>0.) then
    ! Make a 1-d copy of the start and final grids and the source velocity
    do k=1,nz
      h1(k) = h_old_u(I,j,k)
      h2(k) = h_new_u(I,j,k)
      u_src(k) = u(I,j,k)
    enddo

    if (CS%id_remap_delta_integ_u2>0) then
      u2h_tot = 0.
      do k=1,nz
        u2h_tot = u2h_tot - h1(k) * (u_src(k)**2)
      enddo
    endif

    call remapping_core_h(CS%vel_remapCS, nz, h1, u_src, nz, h2, u_tgt)

    if (variance_option .and. CS%conserve_ke) then
    ! Conserve ke_u by correcting baroclinic component.
    ! Assumes total depth doesn't change during remap, and
    ! that \int u(z) dz doesn't change during remap.
      ! First get barotropic component
      u_bt = 0.0
      do k=1,nz
        u_bt = u_bt + h2(k) * u_tgt(k) ! Dimensions [H L T-1]
      enddo
      u_bt = u_bt / (sum(h2(1:nz)) + GV%H_subroundoff) ! Dimensions return to [L T-1]
      ! Next get baroclinic ke = \int (u-u_bt)^2 from source and target
      ke_c_src = 0.0
      ke_c_tgt = 0.0
      do k=1,nz
        ke_c_src = ke_c_src + h1(k) * (u_src(k) - u_bt)**2
        ke_c_tgt = ke_c_tgt + h2(k) * (u_tgt(k) - u_bt)**2
      enddo
      ! Next rescale baroclinic component on target grid to conserve ke
      ! The values 1.5625 = 1.25**2 and 1.25 below mean that the KE-conserving
      ! correction cannot amplify the baroclinic part of velocity by more
      ! than 25%. This threshold is somewhat arbitrary. It was added to
      ! prevent unstable behavior when the amplification factor is large.
      if (ke_c_src < 1.5625 * ke_c_tgt) then
        rescale_coef = sqrt(ke_c_src / ke_c_tgt)
      else
        rescale_coef = 1.25
      endif
      do k=1,nz
        u_tgt(k) = u_bt + rescale_coef * (u_tgt(k) - u_bt)
      enddo
    endif

    if (CS%id_remap_delta_integ_u2>0) then
      do k=1,nz
        u2h_tot = u2h_tot + h2(k) * (u_tgt(k)**2)
      enddo
      du2h_tot(I,j) = u2h_tot * I_dt
    endif

    if ((CS%BBL_h_vel_mask > 0.0) .and. (CS%h_vel_mask > 0.0)) &
      call mask_near_bottom_vel(u_tgt, h2, CS%BBL_h_vel_mask, CS%h_vel_mask, nz)

    ! Copy the column of new velocities back to the 3-d array
    do k=1,nz
      u(I,j,k) = u_tgt(k)
    enddo !k
  endif ; enddo ; enddo

  if (CS%id_remap_delta_integ_u2>0) call post_data(CS%id_remap_delta_integ_u2, du2h_tot, CS%diag)

  if (show_call_tree) call callTree_waypoint("u remapped (ALE_remap_velocities)")


  ! --- Remap v profiles from the source vertical grid onto the new target grid.

  !$OMP parallel do default(shared) private(h1,h2,v_src,h_mask_vel,v_tgt, &
  !$OMP                                     v_bt,ke_c_src,ke_c_tgt,rescale_coef, &
  !$OMP                                     u2h_tot,v2h_tot)
  do J=G%JscB,G%JecB ; do i=G%isc,G%iec ; if (G%mask2dCv(i,J)>0.) then

    do k=1,nz
      h1(k) = h_old_v(i,J,k)
      h2(k) = h_new_v(i,J,k)
      v_src(k) = v(i,J,k)
    enddo

    if (CS%id_remap_delta_integ_v2>0) then
      v2h_tot = 0.
      do k=1,nz
        v2h_tot = v2h_tot - h1(k) * (v_src(k)**2)
      enddo
    endif

    call remapping_core_h(CS%vel_remapCS, nz, h1, v_src, nz, h2, v_tgt)

    if (variance_option .and. CS%conserve_ke) then
    ! Conserve ke_v by correcting baroclinic component.
    ! Assumes total depth doesn't change during remap, and
    ! that \int v(z) dz doesn't change during remap.
      ! First get barotropic component
      v_bt = 0.0
      do k=1,nz
        v_bt = v_bt + h2(k) * v_tgt(k) ! Dimensions [H L T-1]
      enddo
      v_bt = v_bt / (sum(h2(1:nz)) + GV%H_subroundoff) ! Dimensions return to [L T-1]
      ! Next get baroclinic ke = \int (u-u_bt)^2 from source and target
      ke_c_src = 0.0
      ke_c_tgt = 0.0
      do k=1,nz
        ke_c_src = ke_c_src + h1(k) * (v_src(k) - v_bt)**2
        ke_c_tgt = ke_c_tgt + h2(k) * (v_tgt(k) - v_bt)**2
      enddo
      ! Next rescale baroclinic component on target grid to conserve ke
      if (ke_c_src < 1.5625 * ke_c_tgt) then
        rescale_coef = sqrt(ke_c_src / ke_c_tgt)
      else
        rescale_coef = 1.25
      endif
      do k=1,nz
        v_tgt(k) = v_bt + rescale_coef * (v_tgt(k) - v_bt)
      enddo
    endif

    if (CS%id_remap_delta_integ_v2>0) then
      do k=1,nz
        v2h_tot = v2h_tot + h2(k) * (v_tgt(k)**2)
      enddo
      dv2h_tot(I,j) = v2h_tot * I_dt
    endif

    if ((CS%BBL_h_vel_mask > 0.0) .and. (CS%h_vel_mask > 0.0)) then
      call mask_near_bottom_vel(v_tgt, h2, CS%BBL_h_vel_mask, CS%h_vel_mask, nz)
    endif

    ! Copy the column of new velocities back to the 3-d array
    do k=1,nz
      v(i,J,k) = v_tgt(k)
    enddo !k
  endif ; enddo ; enddo

  if (CS%id_remap_delta_integ_v2>0) call post_data(CS%id_remap_delta_integ_v2, dv2h_tot, CS%diag)

  if (show_call_tree) call callTree_waypoint("v remapped (ALE_remap_velocities)")
  if (show_call_tree) call callTree_leave("ALE_remap_velocities()")

end subroutine ALE_remap_velocities

!> Interpolate to find an updated array of values at interfaces after remapping.
subroutine ALE_remap_interface_vals(CS, G, GV, h_old, h_new, int_val)
  type(ALE_CS),                              intent(in)    :: CS       !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G        !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV       !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_old    !< Thickness of source grid
                                                                       !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_new    !< Thickness of destination grid
                                                                       !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                                             intent(inout) :: int_val  !< The interface values to interpolate [A]

  real :: val_src(GV%ke+1)  ! A column of interface values on the source grid [A]
  real :: val_tgt(GV%ke+1)  ! A column of interface values on the target grid [A]
  real :: h_src(GV%ke)      ! A column of source grid layer thicknesses [H ~> m or kg m-2]
  real :: h_tgt(GV%ke)      ! A column of target grid layer thicknesses [H ~> m or kg m-2]
  integer :: i, j, k, nz

  nz = GV%ke

  do j=G%jsc,G%jec ; do i=G%isc,G%iec ; if (G%mask2dT(i,j)>0.) then
    do k=1,nz
      h_src(k) = h_old(i,j,k)
      h_tgt(k) = h_new(i,j,k)
    enddo

    do K=1,nz+1
      val_src(K) = int_val(i,j,K)
    enddo

    call interpolate_column(nz, h_src, val_src, nz, h_tgt, val_tgt, .false.)

    do K=1,nz+1
      int_val(i,j,K) = val_tgt(K)
    enddo
  endif ; enddo ; enddo

end subroutine ALE_remap_interface_vals

!> Interpolate to find an updated array of values at vertices of tracer cells after remapping.
subroutine ALE_remap_vertex_vals(CS, G, GV, h_old, h_new, vert_val)
  type(ALE_CS),                              intent(in)    :: CS       !< ALE control structure
  type(ocean_grid_type),                     intent(in)    :: G        !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV       !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_old    !< Thickness of source grid
                                                                       !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_new    !< Thickness of destination grid
                                                                       !! [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1), &
                                             intent(inout) :: vert_val  !< The interface values to interpolate [A]

  real :: val_src(GV%ke+1)  ! A column of interface values on the source grid [A]
  real :: val_tgt(GV%ke+1)  ! A column of interface values on the target grid [A]
  real :: h_src(GV%ke)      ! A column of source grid layer thicknesses [H ~> m or kg m-2]
  real :: h_tgt(GV%ke)      ! A column of target grid layer thicknesses [H ~> m or kg m-2]
  real :: I_mask_sum        ! The inverse of the tracer point masks surrounding a corner [nondim]
  integer :: i, j, k, nz

  nz = GV%ke

  do J=G%JscB,G%JecB ; do I=G%IscB,G%IecB
    if ((G%mask2dT(i,j) + G%mask2dT(i+1,j+1)) + (G%mask2dT(i+1,j) + G%mask2dT(i,j+1)) > 0.0 ) then
      I_mask_sum = 1.0 / ((G%mask2dT(i,j) + G%mask2dT(i+1,j+1)) + (G%mask2dT(i+1,j) + G%mask2dT(i,j+1)))

    do k=1,nz
      h_src(k) = ((G%mask2dT(i,j) * h_old(i,j,k) + G%mask2dT(i+1,j+1) * h_old(i+1,j+1,k)) + &
                  (G%mask2dT(i+1,j) * h_old(i+1,j,k) + G%mask2dT(i,j+1) * h_old(i,j+1,k)) ) * I_mask_sum
      h_tgt(k) = ((G%mask2dT(i,j) * h_new(i,j,k) + G%mask2dT(i+1,j+1) * h_new(i+1,j+1,k)) + &
                  (G%mask2dT(i+1,j) * h_new(i+1,j,k) + G%mask2dT(i,j+1) * h_new(i,j+1,k)) ) * I_mask_sum
    enddo

    do K=1,nz+1
      val_src(K) = vert_val(I,J,K)
    enddo

    call interpolate_column(nz, h_src, val_src, nz, h_tgt, val_tgt, .false.)

    do K=1,nz+1
      vert_val(I,J,K) = val_tgt(K)
    enddo
  endif ; enddo ; enddo

end subroutine ALE_remap_vertex_vals

!> Mask out thicknesses to 0 when their running sum exceeds a specified value.
subroutine apply_partial_cell_mask(h1, h_mask)
  real, dimension(:), intent(inout) :: h1 !< A column of thicknesses to be masked out after their
                                          !! running vertical sum exceeds h_mask [H ~> m or kg m-2]
  real,               intent(in)    :: h_mask !< The depth after which the thicknesses in h1 are
                                          !! masked out [H ~> m or kg m-2]
  ! Local variables
  real :: h1_rsum  ! The running sum of h1 [H ~> m or kg m-2]
  integer :: k

  h1_rsum = 0.0
  do k=1,size(h1)
    if (h1(k) > h_mask - h1_rsum) then
      ! This thickness is reduced because it extends below the shallower neighboring bathymetry.
      h1(k) = max(h_mask - h1_rsum, 0.0)
      h1_rsum = h_mask
    else
      h1_rsum = h1_rsum + h1(k)
    endif
  enddo
end subroutine apply_partial_cell_mask


!> Zero out velocities in a column in very thin layers near the seafloor
subroutine mask_near_bottom_vel(vel, h, h_BBL, h_thin, nk)
  integer, intent(in)    :: nk      !< The number of layers in this column
  real,    intent(inout) :: vel(nk) !< The velocity component being zeroed out [L T-1 ~> m s-1]
  real,    intent(in)    :: h(nk)   !< The layer thicknesses at velocity points  [H ~> m or kg m-2]
  real,    intent(in)    :: h_BBL   !< The thickness of the near-bottom region over which to apply
                                    !! the filtering [H ~> m or kg m-2]
  real,    intent(in)    :: h_thin  !< A layer thickness below which the filtering is applied [H ~> m or kg m-2]

  ! Local variables
  real :: h_from_bot  ! The distance between the top of a layer and the seafloor [H ~> m or kg m-2]
  integer :: k

  if ((h_BBL < 0.0) .or. (h_thin < 0.0)) return

  h_from_bot = 0.0
  do k=nk,1,-1
    h_from_bot = h_from_bot + h(k)
    if (h_from_bot > h_BBL) return
    ! Set the velocity to zero in thin, near-bottom layers.
    if (h(k) <= h_thin) vel(k) = 0.0
  enddo !k

end subroutine mask_near_bottom_vel


!> Remaps a single scalar between grids described by thicknesses h_src and h_dst.
!! h_dst must be dimensioned as a model array with GV%ke layers while h_src can
!! have an arbitrary number of layers specified by nk_src.
subroutine ALE_remap_scalar(CS, G, GV, nk_src, h_src, s_src, h_dst, s_dst, all_cells, old_remap)
  type(remapping_CS),                      intent(in)    :: CS        !< Remapping control structure
  type(ocean_grid_type),                   intent(in)    :: G         !< Ocean grid structure
  type(verticalGrid_type),                 intent(in)    :: GV        !< Ocean vertical grid structure
  integer,                                 intent(in)    :: nk_src    !< Number of levels on source grid
  real, dimension(SZI_(G),SZJ_(G),nk_src), intent(in)    :: h_src     !< Level thickness of source grid
                                                                      !! [H ~> m or kg m-2] or other units
                                                                      !! if H_neglect is provided
  real, dimension(SZI_(G),SZJ_(G),nk_src), intent(in)    :: s_src     !< Scalar on source grid, in arbitrary units [A]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),intent(in)   :: h_dst     !< Level thickness of destination grid in the
                                                                      !! same units as h_src, often [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),intent(inout) :: s_dst    !< Scalar on destination grid, in the same
                                                                      !! arbitrary units as s_src [A]
  logical, optional,                       intent(in)    :: all_cells !< If false, only reconstruct for
                                                                      !! non-vanished cells. Use all vanished
                                                                      !! layers otherwise (default).
  logical, optional,                       intent(in)    :: old_remap !< If true, use the old "remapping_core_w"
                                                                      !! method, otherwise use "remapping_core_h".
   ! Local variables
  integer :: i, j, k, n_points
  real :: dx(GV%ke+1) ! Change in interface position [H ~> m or kg m-2]
  logical :: ignore_vanished_layers, use_remapping_core_w

  ignore_vanished_layers = .false.
  if (present(all_cells)) ignore_vanished_layers = .not. all_cells
  use_remapping_core_w = .false.
  if (present(old_remap)) use_remapping_core_w = old_remap
  n_points = nk_src

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
                              GV%ke, dx, s_dst(i,j,:))
      else
        call remapping_core_h(CS, n_points, h_src(i,j,1:n_points), s_src(i,j,1:n_points), &
                              GV%ke, h_dst(i,j,:), s_dst(i,j,:))
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
                           intent(inout) :: S_t  !< Salinity at the top edge of each layer [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: S_b  !< Salinity at the bottom edge of each layer [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_t  !< Temperature at the top edge of each layer [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_b  !< Temperature at the bottom edge of each layer [C ~> degC]
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
                           intent(in)    :: Q    !< 3d scalar array, in arbitrary units [A]
  logical,                 intent(in)    :: bdry_extrap !< If true, use high-order boundary
                                                 !! extrapolation within boundary cells
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: Q_t  !< Scalar at the top edge of each layer [A]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: Q_b  !< Scalar at the bottom edge of each layer [A]
  ! Local variables
  integer :: i, j, k
  real :: slp(GV%ke) ! Tracer slope times the cell width [A]
  real :: mslp       ! Monotonized tracer slope times the cell width [A]
  real :: h_neglect  ! Tiny thicknesses used in remapping [H ~> m or kg m-2]

  if (CS%answer_date >= 20190101) then
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
                           intent(inout) :: S_t  !< Salinity at the top edge of each layer [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: S_b  !< Salinity at the bottom edge of each layer [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_t  !< Temperature at the top edge of each layer [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: T_b  !< Temperature at the bottom edge of each layer [C ~> degC]
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

  if (CS%answer_date >= 20190101) then
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
                                  answer_date=CS%answer_date )
    call PPM_reconstruction( GV%ke, hTmp, tmp, ppol_E, ppol_coefs, h_neglect, &
                                  answer_date=CS%answer_date )
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
    if (CS%answer_date < 20190101) then
      call edge_values_implicit_h4( GV%ke, hTmp, tmp, ppol_E, h_neglect=1.0e-10*GV%m_to_H, &
                                  answer_date=CS%answer_date )
    else
      call edge_values_implicit_h4( GV%ke, hTmp, tmp, ppol_E, h_neglect=GV%H_subroundoff, &
                                  answer_date=CS%answer_date )
    endif
    call PPM_reconstruction( GV%ke, hTmp, tmp, ppol_E, ppol_coefs, h_neglect, &
                                  answer_date=CS%answer_date )
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

  real, dimension(CS%nk+1) :: ALE_getCoordinate !< The coordinate positions, in the appropriate units
                                                !! of the target coordinate, e.g. [Z ~> m] for z*,
                                                !! non-dimensional for sigma, etc.
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
  real :: w  ! An implicit weighting estimate [nondim]

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
                    ! as function of coordinate mode.

end subroutine ALE_updateVerticalGridType


!> Write the vertical coordinate information into a file.
!! This subroutine writes out a file containing any available data related
!! to the vertical grid used by the MOM ocean model when in ALE mode.
subroutine ALE_writeCoordinateFile( CS, GV, directory )
  type(ALE_CS),            pointer     :: CS         !< module control structure
  type(verticalGrid_type), intent(in)  :: GV         !< ocean vertical grid structure
  character(len=*),        intent(in)  :: directory  !< directory for writing grid info

  character(len=240) :: filepath

  filepath = trim(directory) // trim("Vertical_coordinate.nc")

  call write_regrid_file(CS%regridCS, GV, filepath)

end subroutine ALE_writeCoordinateFile

!> Set h to coordinate values for fixed coordinate systems
subroutine ALE_initThicknessToCoord( CS, G, GV, h, height_units )
  type(ALE_CS), intent(inout)                            :: CS  !< module control structure
  type(ocean_grid_type), intent(in)                      :: G   !< module grid structure
  type(verticalGrid_type), intent(in)                    :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: h   !< layer thickness in thickness units
                                                                !! [H ~> m or kg m-2] or height units [Z ~> m]
  logical,                          optional, intent(in) :: height_units !< If present and true, the
                                                                !! thicknesses are in height units

  ! Local variables
  real :: scale ! A scaling value for the thicknesses [nondim] or [H Z-1 ~> nondim or kg m-3]
  integer :: i, j

  scale = GV%Z_to_H
  if (present(height_units)) then ; if (height_units) scale = 1.0 ; endif
  do j = G%jsd,G%jed ; do i = G%isd,G%ied
    h(i,j,:) = scale * getStaticThickness( CS%regridCS, 0., G%bathyT(i,j)+G%Z_ref )
  enddo ; enddo

end subroutine ALE_initThicknessToCoord

end module MOM_ALE
