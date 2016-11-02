!> This module contains the main regridding routines. 
!! Regridding comprises two steps:
!! (1) Interpolation and creation of a new grid based on target interface 
!!     densities (or any other criterion).
!! (2) Remapping of quantities between old grid and new grid.
!! Original module written by Laurent White, 2008.06.09
module MOM_ALE

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator,    only : register_diag_field, post_data, diag_ctrl, time_type
use MOM_EOS,              only : calculate_density
use MOM_error_handler,    only : MOM_error, FATAL, WARNING
use MOM_error_handler,    only : callTree_showQuery
use MOM_error_handler,    only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,      only : get_param, param_file_type, log_param
use MOM_io,               only : file_exists, field_exists, MOM_read_data
use MOM_io,               only : vardesc, var_desc, fieldtype, SINGLE_FILE
use MOM_io,               only : create_file, write_field, close_file, slasher
use MOM_regridding,       only : initialize_regridding, regridding_main, end_regridding
use MOM_regridding,       only : uniformResolution
use MOM_regridding,       only : inflate_vanished_layers_old, setCoordinateResolution
use MOM_regridding,       only : set_target_densities_from_GV, set_target_densities
use MOM_regridding,       only : regriddingCoordinateModeDoc, DEFAULT_COORDINATE_MODE
use MOM_regridding,       only : regriddingInterpSchemeDoc, regriddingDefaultInterpScheme
use MOM_regridding,       only : regriddingDefaultBoundaryExtrapolation
use MOM_regridding,       only : regriddingDefaultMinThickness
use MOM_regridding,       only : set_regrid_max_depths
use MOM_regridding,       only : set_regrid_max_thickness
use MOM_regridding,       only : regridding_CS, set_regrid_params
use MOM_regridding,       only : getCoordinateInterfaces, getCoordinateResolution
use MOM_regridding,       only : getCoordinateUnits, getCoordinateShortName
use MOM_regridding,       only : getStaticThickness
use MOM_remapping,        only : initialize_remapping, end_remapping
use MOM_remapping,        only : remapping_core_h, remapping_core_w
use MOM_remapping,        only : remappingSchemesDoc, remappingDefaultScheme
use MOM_remapping,        only : remapping_CS, dzFromH1H2
use MOM_string_functions, only : uppercase, extractWord
use MOM_tracer_registry,  only : tracer_registry_type
use MOM_variables,        only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,     only : verticalGrid_type

use regrid_defs,          only : PRESSURE_RECONSTRUCTION_PLM
!use regrid_consts,       only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts,        only : coordinateUnits, coordinateMode
use regrid_consts,        only : REGRIDDING_ZSTAR, REGRIDDING_RHO
use regrid_consts,        only : REGRIDDING_HYCOM1, REGRIDDING_SLIGHT
use regrid_edge_values,   only : edge_values_implicit_h4
use PLM_functions,        only : PLM_reconstruction, PLM_boundary_extrapolation
use PPM_functions,        only : PPM_reconstruction, PPM_boundary_extrapolation
use P1M_functions,        only : P1M_interpolation,  P1M_boundary_extrapolation
use P3M_functions,        only : P3M_interpolation,  P3M_boundary_extrapolation


implicit none ; private
#include <MOM_memory.h>


!> ALE control structure 
type, public :: ALE_CS
  private

  logical :: boundary_extrapolation_for_pressure !< Indicate whether high-order boundary 
                                                 !! extrapolation should be used within boundary cells

  logical :: reconstructForPressure = .false.    !< Indicates whether integrals for FV 
                                                 !! pressure gradient calculation will
                                                 !! use reconstruction of T/S.
                                                 !! By default, it is true if regridding 
                                                 !! has been initialized, otherwise false.

  integer :: pressureReconstructionScheme        !< Form of the reconstruction of T/S 
                                                 !! for FV pressure gradient calculation.
                                                 !! By default, it is =1 (PLM)

  logical :: remap_uv_using_old_alg              !< If true, uses the old "remapping via a delta z"
                                                 !! method. If False, uses the new method that
                                                 !! remaps between grids described by h.

  real :: regrid_time_scale !< The time-scale used in blending between the current (old) grid
                            !! and the target (new) grid. (s)

  type(regridding_CS) :: regridCS !< Regridding parameters and work arrays
  type(remapping_CS)  :: remapCS  !< Remapping parameters and work arrays

  integer :: nk              !< Used only for queries, not directly by this module
  integer :: degree_linear=1 !< Degree of linear piecewise polynomial
  integer :: degree_parab=2  !< Degree of parabolic piecewise polynomial

  logical :: remap_after_initialization !<   Indicates whether to regrid/remap after initializing the state.

  logical :: show_call_tree !< For debugging
  real    :: C_p            !< seawater heat capacity (J/(kg deg C))

  ! for diagnostics 
  type(diag_ctrl), pointer           :: diag                          !< structure to regulate output
  integer, dimension(:), allocatable :: id_tracer_remap_tendency      !< diagnostic id 
  integer, dimension(:), allocatable :: id_Htracer_remap_tendency     !< diagnostic id 
  integer, dimension(:), allocatable :: id_Htracer_remap_tendency_2d  !< diagnostic id 
  logical, dimension(:), allocatable :: do_tendency_diag              !< flag for doing diagnostics
  integer                            :: id_dzRegrid

end type

! Publicly available functions
public ALE_init
public ALE_end
public ALE_main
public ALE_main_offline
public ALE_offline_tracer_final
public ALE_build_grid
public ALE_remap_scalar
public pressure_gradient_plm
public pressure_gradient_ppm
public usePressureReconstruction
public pressureReconstructionScheme
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

contains

!> This routine is typically called (from initialize_MOM in file MOM.F90)
!! before the main time integration loop to initialize the regridding stuff.
!! We read the MOM_input file to register the values of different
!! regridding/remapping parameters.
subroutine ALE_init( param_file, GV, max_depth, CS)
  type(param_file_type),   intent(in) :: param_file !< Parameter file
  type(verticalGrid_type), intent(in) :: GV         !< Ocean vertical grid structure
  real,                    intent(in) :: max_depth  !< The maximum depth of the ocean, in m.
  type(ALE_CS),            pointer    :: CS         !< Module control structure

  ! Local variables
  real, dimension(:), allocatable :: dz
  character(len=40)               :: mod = "MOM_ALE" ! This module's name.
  character(len=80)               :: string ! Temporary strings
  real                            :: filter_shallow_depth, filter_deep_depth
  logical                         :: check_reconstruction
  logical                         :: check_remapping
  logical                         :: force_bounds_in_subcell
  logical                         :: local_logical

  if (associated(CS)) then
    call MOM_error(WARNING, "ALE_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%show_call_tree = callTree_showQuery()
  if (CS%show_call_tree) call callTree_enter("ALE_init(), MOM_ALE.F90")

  ! --- BOUNDARY EXTRAPOLATION --
  ! This sets whether high-order (rather than PCM) reconstruction schemes
  ! should be used within boundary cells
  call get_param(param_file, mod, "BOUNDARY_EXTRAPOLATION_PRESSURE", &
                 CS%boundary_extrapolation_for_pressure, &
                 "When defined, the reconstruction is extrapolated\n"//&
                 "within boundary cells rather than assume PCM for the.\n"//&
                 "calculation of pressure. e.g. if PPM is used, a\n"//&
                 "PPM reconstruction will also be used within\n"//&
                 "boundary cells.", default=.true.)

  ! --- PRESSURE GRADIENT CALCULATION ---
  call get_param(param_file, mod, "RECONSTRUCT_FOR_PRESSURE", &
                 CS%reconstructForPressure , &
                 "If True, use vertical reconstruction of T/S within\n"//&
                 "the integrals of teh FV pressure gradient calculation.\n"//&
                 "If False, use the constant-by-layer algorithm.\n"//&
                 "By default, this is True when using ALE and False otherwise.", &
                 default=.true. )

  call get_param(param_file, mod, "PRESSURE_RECONSTRUCTION_SCHEME", &
                 CS%pressureReconstructionScheme, &
                 "Type of vertical reconstruction of T/S to use in integrals\n"//&
                 "within the FV pressure gradient calculation."//&
                 " 1: PLM reconstruction.\n"//&
                 " 2: PPM reconstruction.", default=PRESSURE_RECONSTRUCTION_PLM)

  call get_param(param_file, mod, "REMAP_UV_USING_OLD_ALG", &
                 CS%remap_uv_using_old_alg, &
                 "If true, uses the old remapping-via-a-delta-z method for\n"//&
                 "remapping u and v. If false, uses the new method that remaps\n"//&
                 "between grids described by an old and new thickness.", &
                 default=.true.)

  ! Initialize and configure regridding
  allocate( dz(GV%ke) )
  call ALE_initRegridding( GV, max_depth, param_file, mod, CS%regridCS, dz )
  deallocate( dz )

  ! Initialize and configure remapping
  call get_param(param_file, mod, "REMAPPING_SCHEME", string, &
                 "This sets the reconstruction scheme used\n"//&
                 "for vertical remapping for all variables.\n"//&
                 "It can be one of the following schemes:\n"//&
                 trim(remappingSchemesDoc), default=remappingDefaultScheme)
  call get_param(param_file, mod, "FATAL_CHECK_RECONSTRUCTIONS", check_reconstruction, &
                 "If true, cell-by-cell reconstructions are checked for\n"//&
                 "consistency and if non-monotonicty or an inconsistency is\n"//&
                 "detected then a FATAL error is issued.", default=.false.)
  call get_param(param_file, mod, "FATAL_CHECK_REMAPPING", check_remapping, &
                 "If true, the results of remapping are checked for\n"//&
                 "conservation and new extrema and if an inconsistency is\n"//&
                 "detected then a FATAL error is issued.", default=.false.)
  call get_param(param_file, mod, "REMAP_BOUND_INTERMEDIATE_VALUES", force_bounds_in_subcell, &
                 "If true, the values on the intermediate grid used for remapping\n"//&
                 "are forced to be bounded, which might not be the case due to\n"//&
                 "round off.", default=.false.)
  call initialize_remapping( CS%remapCS, string, &
                             boundary_extrapolation=.false., &
                             check_reconstruction=check_reconstruction, &
                             check_remapping=check_remapping, &
                             force_bounds_in_subcell=force_bounds_in_subcell)

  call get_param(param_file, mod, "REMAP_AFTER_INITIALIZATION", CS%remap_after_initialization, &
                 "If true, applies regridding and remapping immediately after\n"//&
                 "initialization so that the state is ALE consistent. This is a\n"//&
                 "legacy step and should not be needed if the initialization is\n"//&
                 "consistent with the coordinate mode.", default=.true.)

  call get_param(param_file, mod, "REGRID_TIME_SCALE", CS%regrid_time_scale, &
                 "The time-scale used in blending between the current (old) grid\n"//&
                 "and the target (new) grid. A short time-scale favors the target\n"//&
                 "grid (0. or anything less than DT_THERM) has no memory of the old\n"//&
                 "grid. A very long time-scale makes the model more Lagrangian.", &
                 units="s", default=0.)
  call get_param(param_file, mod, "REGRID_FILTER_SHALLOW_DEPTH", filter_shallow_depth, &
                 "The depth above which no time-filtering is applied. Above this depth\n"//&
                 "final grid exactly matches the target (new) grid.", units="m", default=0.)
  call get_param(param_file, mod, "REGRID_FILTER_DEEP_DEPTH", filter_deep_depth, &
                 "The depth below which full time-filtering is applied with time-scale\n"//&
                 "REGRID_TIME_SCALE. Between depths REGRID_FILTER_SHALLOW_DEPTH and\n"//&
                 "REGRID_FILTER_SHALLOW_DEPTH the filter wieghts adopt a cubic profile.", &
                 units="m", default=0.)
  call set_regrid_params(CS%regridCS, depth_of_time_filter_shallow=filter_shallow_depth*GV%m_to_H, &
                                      depth_of_time_filter_deep=filter_deep_depth*GV%m_to_H)
  call get_param(param_file, mod, "REGRID_USE_OLD_DIRECTION", local_logical, &
                 "If true, the regridding ntegrates upwards from the bottom for\n"//&
                 "interface positions, much as the main model does. If false\n"//&
                 "regridding integrates downward, consistant with the remapping\n"//&
                 "code.", default=.true., do_not_log=.true.)
  call set_regrid_params(CS%regridCS, integrate_downward_for_e=.not.local_logical)

  ! Keep a record of values for subsequent queries
  CS%nk = GV%ke

  if (CS%show_call_tree) call callTree_leave("ALE_init()")
end subroutine ALE_init

!> Initialize diagnostics for the ALE module. 
subroutine ALE_register_diags(Time, G, diag, C_p, Reg, CS)
  type(time_type),target,     intent(in)  :: Time  !< Time structure
  type(ocean_grid_type),      intent(in)  :: G     !< Grid structure
  type(diag_ctrl), target,    intent(in)  :: diag  !< Diagnostics control structure
  real,                       intent(in)  :: C_p   !< seawater heat capacity (J/(kg deg C))
  type(tracer_registry_type), pointer     :: Reg   !< Tracer registry
  type(ALE_CS), pointer                   :: CS    !< Module control structure 

  integer :: m, ntr, nsize 

  if (associated(Reg)) then
    ntr = Reg%ntr
  else
    ntr = 0
  endif
  nsize = max(1,ntr)

  CS%diag => diag
  CS%C_p  = C_p

  allocate(CS%id_tracer_remap_tendency(nsize))
  allocate(CS%id_Htracer_remap_tendency(nsize))
  allocate(CS%id_Htracer_remap_tendency_2d(nsize))
  allocate(CS%do_tendency_diag(nsize)) 
  CS%do_tendency_diag(:)             = .false.
  CS%id_tracer_remap_tendency(:)     = -1
  CS%id_Htracer_remap_tendency(:)    = -1
  CS%id_Htracer_remap_tendency_2d(:) = -1

  CS%id_dzRegrid = register_diag_field('ocean_model','dzRegrid',diag%axesTi,Time, &
      'Change in interface height due to ALE regridding', 'meter')

  if(ntr > 0) then 

    do m=1,ntr
      if(trim(Reg%Tr(m)%name) == 'T') then 

        CS%id_tracer_remap_tendency(m) = register_diag_field('ocean_model',                &
        trim(Reg%Tr(m)%name)//'_tendency_vert_remap', diag%axesTL, Time,                   &
        'Tendency from vertical remapping for tracer concentration '//trim(Reg%Tr(m)%name),&
        'degC/s')

        CS%id_Htracer_remap_tendency(m) = register_diag_field('ocean_model',&
        trim(Reg%Tr(m)%name)//'h_tendency_vert_remap', diag%axesTL, Time,   &
        'Tendency from vertical remapping for heat',                        &
         'W/m2')

        CS%id_Htracer_remap_tendency_2d(m) = register_diag_field('ocean_model',&
        trim(Reg%Tr(m)%name)//'h_tendency_vert_remap_2d', diag%axesT1, Time,   &
        'Vertical sum of tendency from vertical remapping for heat',           &
         'W/m2')

      else 

        CS%id_tracer_remap_tendency(m) = register_diag_field('ocean_model',                &
        trim(Reg%Tr(m)%name)//'_tendency_vert_remap', diag%axesTL, Time,                   & 
        'Tendency from vertical remapping for tracer concentration '//trim(Reg%Tr(m)%name),&
        'tracer conc / sec')

        CS%id_Htracer_remap_tendency(m) = register_diag_field('ocean_model',         &
        trim(Reg%Tr(m)%name)//'h_tendency_vert_remap', diag%axesTL, Time,            &
        'Tendency from vertical remapping for tracer content '//trim(Reg%Tr(m)%name),&
        'kg m-2 s-1')

        CS%id_Htracer_remap_tendency_2d(m) = register_diag_field('ocean_model',                      &
        trim(Reg%Tr(m)%name)//'h_tendency_vert_remap_2d', diag%axesT1, Time,                         &
        'Vertical sum of tendency from vertical remapping for tracer content '//trim(Reg%Tr(m)%name),&
        'kg m-2 s-1')

      endif 

      if(CS%id_tracer_remap_tendency(m)     > 0) CS%do_tendency_diag(m) = .true.
      if(CS%id_Htracer_remap_tendency(m)    > 0) CS%do_tendency_diag(m) = .true.
      if(CS%id_Htracer_remap_tendency_2d(m) > 0) CS%do_tendency_diag(m) = .true.

    enddo   ! m loop over tracers 

  endif ! ntr > 0  

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
                                                                   !! are to be adjusted (m or Pa)
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
subroutine ALE_main( G, GV, h, u, v, tv, Reg, CS, dt, frac_shelf_h)
  type(ocean_grid_type),                      intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                    intent(in)    :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h   !< Current 3D grid obtained after last time step (m or Pa)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u   !< Zonal velocity field (m/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v   !< Meridional velocity field (m/s)
  type(thermo_var_ptrs),                      intent(inout) :: tv  !< Thermodynamic variable structure
  type(tracer_registry_type),                 pointer       :: Reg !< Tracer registry structure
  type(ALE_CS),                               pointer       :: CS  !< Regridding parameters and options
  real,                             optional, intent(in)    :: dt  !< Time step between calls to ALE_main()
  real, dimension(:,:),             optional, pointer       :: frac_shelf_h !< Fractional ice shelf coverage
  ! Local variables
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)+1) :: dzRegrid ! The change in grid interface positions
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_new ! New 3D grid obtained after last time step (m or Pa)
  integer :: nk, i, j, k, isc, iec, jsc, jec
  logical :: ice_shelf

  nk = GV%ke; isc = G%isc; iec = G%iec; jsc = G%jsc; jec = G%jec

  ice_shelf = .false.
  if (present(frac_shelf_h)) then
    if (associated(frac_shelf_h)) ice_shelf = .true.
  endif

  if (CS%show_call_tree) call callTree_enter("ALE_main(), MOM_ALE.F90")

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

  ! Remap all variables from old grid h onto new grid h_new

  call remap_all_state_vars( CS%remapCS, CS, G, GV, h, h_new, -dzRegrid, Reg, &
                             u, v, CS%show_call_tree, dt )

  if (CS%show_call_tree) call callTree_waypoint("state remapped (ALE_main)")

  ! Override old grid with new one. The new grid 'h_new' is built in
  ! one of the 'build_...' routines above.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,nk,h,h_new,CS)
  do k = 1,nk
    do j = jsc-1,jec+1 ; do i = isc-1,iec+1
      h(i,j,k) = h_new(i,j,k)
    enddo ; enddo
  enddo

  if (CS%show_call_tree) call callTree_leave("ALE_main()")

  if (CS%id_dzRegrid>0 .and. present(dt)) call post_data(CS%id_dzRegrid, dzRegrid, CS%diag)


end subroutine ALE_main

!> Takes care of (1) building a new grid and (2) remapping all variables between
!! the old grid and the new grid. The creation of the new grid can be based
!! on z coordinates, target interface densities, sigma coordinates or any
!! arbitrary coordinate system.
subroutine ALE_main_offline( G, GV, h, tv, Reg, CS, dt)
  type(ocean_grid_type),                      intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                    intent(in)    :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h   !< Current 3D grid obtained after last time step (m or Pa)
  type(thermo_var_ptrs),                      intent(inout) :: tv  !< Thermodynamic variable structure
  type(tracer_registry_type),                 pointer       :: Reg !< Tracer registry structure
  type(ALE_CS),                               pointer       :: CS  !< Regridding parameters and options
  real,                             optional, intent(in)    :: dt  !< Time step between calls to ALE_main()
  ! Local variables
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)+1) :: dzRegrid ! The change in grid interface positions
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_new ! New 3D grid obtained after last time step (m or Pa)
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

  call remap_all_state_vars( CS%remapCS, CS, G, GV, h, h_new, -dzRegrid, Reg, &
                             debug=CS%show_call_tree, dt=dt )

  if (CS%show_call_tree) call callTree_waypoint("state remapped (ALE_main)")

  ! Override old grid with new one. The new grid 'h_new' is built in
  ! one of the 'build_...' routines above.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,nk,h,h_new,CS)
  do k = 1,nk
    do j = jsc-1,jec+1 ; do i = isc-1,iec+1
      h(i,j,k) = h_new(i,j,k)
    enddo ; enddo
  enddo

  if (CS%show_call_tree) call callTree_leave("ALE_main()")
  if (CS%id_dzRegrid>0 .and. present(dt)) call post_data(CS%id_dzRegrid, dzRegrid, CS%diag)

end subroutine ALE_main_offline

!> Remaps all tracers from h onto h_target. This is intended to be called when tracers
!! are done offline. In the case where transports don't quite conserve, we still want to
!! make sure that layer thicknesses offline do not drift too far away from the online model
subroutine ALE_offline_tracer_final( G, GV, h, h_target, Reg, CS)
  type(ocean_grid_type),                      intent(in)    :: G   !< Ocean grid informations
  type(verticalGrid_type),                    intent(in)    :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h   !< Current 3D grid obtained after last time step (m or Pa)
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h_target   !< Current 3D grid obtained after last time step (m or Pa)
  type(tracer_registry_type),                 pointer       :: Reg !< Tracer registry structure
  type(ALE_CS),                               pointer       :: CS  !< Regridding parameters and options
  ! Local variables

  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)+1) :: dzRegrid ! The change in grid interface positions
  integer :: nk, i, j, k, isc, iec, jsc, jec

  nk = GV%ke; isc = G%isc; iec = G%iec; jsc = G%jsc; jec = G%jec

  if (CS%show_call_tree) call callTree_enter("ALE_offline_tracer_final(), MOM_ALE.F90")
  
  ! It does not seem that remap_all_state_vars uses dzRegrid for tracers, only for u, v
  dzRegrid(:,:,:) = 0.0

  call check_grid( G, GV, h, 0. )
  call check_grid( G, GV, h_target, 0. )

  if (CS%show_call_tree) call callTree_waypoint("Source and target grids checked (ALE_offline_tracer)")

  ! Remap all variables from old grid h onto new grid h_new

  call remap_all_state_vars( CS%remapCS, CS, G, GV, h, h_target, -dzRegrid, Reg, &
                             debug=CS%show_call_tree )

  if (CS%show_call_tree) call callTree_waypoint("state remapped (ALE_offline_tracer)")

  ! Override old grid with new one. The new grid 'h_new' is built in
  ! one of the 'build_...' routines above.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,nk,h,h_target,CS)
  do k = 1,nk
    do j = jsc-1,jec+1 ; do i = isc-1,iec+1
      h(i,j,k) = h_target(i,j,k)
    enddo ; enddo
  enddo

  if (CS%show_call_tree) call callTree_leave("ALE_offline_tracer()")

end subroutine ALE_offline_tracer_final

!> Check grid for negative thicknesses
subroutine check_grid( G, GV, h, threshold )
  type(ocean_grid_type),                     intent(in) :: G !< Ocean grid structure
  type(verticalGrid_type),                   intent(in) :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< Current 3D grid obtained after the last time step (H units)
  real,                                      intent(in) :: threshold !< Value below which to flag issues (H units)
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
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)), intent(inout) :: h      !< Current 3D grid obtained after the last time step (m or Pa)
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

!> This routine takes care of remapping all variable between the old and the
!! new grids. When velocity components need to be remapped, thicknesses at
!! velocity points are taken to be arithmetic averages of tracer thicknesses.
!! This routine is called during initialization of the model at time=0, to 
!! remap initiali conditions to the model grid.  It is also called during a
!! time step to update the state.  
subroutine remap_all_state_vars(CS_remapping, CS_ALE, G, GV, h_old, h_new, dxInterface, Reg, u, v, debug, dt)
  type(remapping_CS),                               intent(in)    :: CS_remapping  !< Remapping control structure
  type(ALE_CS),                                     intent(in)    :: CS_ALE        !< ALE control structure 
  type(ocean_grid_type),                            intent(in)    :: G             !< Ocean grid structure
  type(verticalGrid_type),                          intent(in)    :: GV            !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),        intent(in)    :: h_old         !< Thickness of source grid (m or Pa)
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),        intent(in)    :: h_new         !< Thickness of destination grid (m or Pa)
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),      intent(in)    :: dxInterface   !< Change in interface position (Hm or Pa)
  type(tracer_registry_type),                       pointer       :: Reg           !< Tracer registry structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), optional, intent(inout) :: u          !< Zonal velocity component (m/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), optional, intent(inout) :: v          !< Meridional velocity component (m/s)
  logical,                                    optional, intent(in)    :: debug      !< If true, show the call tree
  real,                                       optional, intent(in)    :: dt         !< time step for diagnostics 
  ! Local variables
  integer                                     :: i, j, k, m
  integer                                     :: nz, ntr
  real, dimension(GV%ke+1)                    :: dx
  real, dimension(GV%ke)                      :: h1, u_column
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)) :: work_conc 
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)) :: work_cont
  real, dimension(SZI_(G), SZJ_(G))           :: work_2d 
  real                                        :: Idt, ppt2mks
  real, dimension(GV%ke)                      :: h2
  logical                                     :: show_call_tree

  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug
  if (show_call_tree) call callTree_enter("remap_all_state_vars(), MOM_ALE.F90")

  nz      = GV%ke
  ppt2mks = 0.001 

  if (associated(Reg)) then
    ntr = Reg%ntr
  else
    ntr = 0
  endif

  if(present(dt)) then   
    work_conc(:,:,:) = 0.0
    work_cont(:,:,:) = 0.0
    work_2d(:,:)     = 0.0 
    Idt              = 1.0/dt 
  endif 

  ! Remap tracer
!$OMP parallel default(none) shared(G,GV,h_old,h_new,dxInterface,CS_remapping,nz,Reg,u,v,ntr,show_call_tree, &
!$OMP                               dt,h2,CS_ALE,work_conc,work_cont,work_2d,Idt,ppt2mks) &
!$OMP                       private(h1,dx,u_column)
  if (ntr>0) then
    if (show_call_tree) call callTree_waypoint("remapping tracers (remap_all_state_vars)")
!$OMP do
    do m=1,ntr ! For each tracer 

      do j = G%jsc,G%jec
        do i = G%isc,G%iec

          if (G%mask2dT(i,j)>0.) then

            ! Build the start and final grids
            h1(:) = h_old(i,j,:)
            h2(:) = h_new(i,j,:)
            call remapping_core_h( nz, h1, Reg%Tr(m)%t(i,j,:), nz, h2, u_column, CS_remapping )

            ! Intermediate steps for tendency of tracer concentration and tracer content.
            ! Note: do not merge the two if-tests, since do_tendency_diag(:) is not 
            ! allocated during the time=0 initialization call to this routine. 
            if(present(dt)) then 
              if(CS_ALE%do_tendency_diag(m)) then 
                do k=1,GV%ke
                  work_conc(i,j,k) = (u_column(k)    - Reg%Tr(m)%t(i,j,k)      ) * Idt
                  work_cont(i,j,k) = (u_column(k)*h2(k) - Reg%Tr(m)%t(i,j,k)*h1(k)) * Idt * GV%H_to_kg_m2
                enddo 
              endif 
            endif 

            ! update tracer concentration 
            Reg%Tr(m)%t(i,j,:) = u_column(:)

          endif

        enddo ! i
      enddo ! j


      ! tendency diagnostics.
      ! Note: do not merge the two if-tests if(present(dt)) and 
      ! if(CS_ALE%do_tendency_diag(m)).  The reason is that 
      ! do_tendency_diag(:) is not allocated when this routine is called
      ! during initialization (time=0). So need to keep the if-tests split.  
      if(present(dt)) then 
        if(CS_ALE%do_tendency_diag(m)) then 

          if(CS_ALE%id_tracer_remap_tendency(m) > 0) then 
            call post_data(CS_ALE%id_tracer_remap_tendency(m), work_conc, CS_ALE%diag)
          endif 

          if (CS_ALE%id_Htracer_remap_tendency(m) > 0 .or. CS_ALE%id_Htracer_remap_tendency_2d(m) > 0) then 
            if(trim(Reg%Tr(m)%name) == 'T') then 
              do k=1,GV%ke
                do j = G%jsc,G%jec
                  do i = G%isc,G%iec
                    work_cont(i,j,k) = work_cont(i,j,k) * CS_ALE%C_p
                  enddo
                enddo
              enddo
            elseif(trim(Reg%Tr(m)%name) == 'S') then 
              do k=1,GV%ke
                do j = G%jsc,G%jec
                  do i = G%isc,G%iec
                    work_cont(i,j,k) = work_cont(i,j,k) * ppt2mks 
                  enddo
                enddo
              enddo
            endif 
          endif 

          if (CS_ALE%id_Htracer_remap_tendency(m) > 0) then 
            call post_data(CS_ALE%id_Htracer_remap_tendency(m), work_cont, CS_ALE%diag)
          endif 
          if (CS_ALE%id_Htracer_remap_tendency_2d(m) > 0) then 
            do j = G%jsc,G%jec
              do i = G%isc,G%iec
                work_2d(i,j) = 0.0
                do k = 1,GV%ke
                  work_2d(i,j) = work_2d(i,j) + work_cont(i,j,k)
                enddo 
              enddo 
            enddo 
            call post_data(CS_ALE%id_Htracer_remap_tendency_2d(m), work_2d, CS_ALE%diag)
          endif

        endif 
      endif 

    enddo ! m=1,ntr

  endif   ! endif for ntr > 0

  if (show_call_tree) call callTree_waypoint("tracers remapped (remap_all_state_vars)")

  ! Remap u velocity component
  if ( present(u) ) then
!$OMP do
    do j = G%jsc,G%jec
      do I = G%iscB,G%iecB
        if (G%mask2dCu(i,j)>0.) then
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
          call remapping_core_h( nz, h1, u(I,j,:), nz, h2, u_column, CS_remapping )
          u(I,j,:) = u_column(:)
        endif
      enddo
    enddo
  endif

  if (show_call_tree) call callTree_waypoint("u remapped (remap_all_state_vars)")

  ! Remap v velocity component
  if ( present(v) ) then
!$OMP do
    do J = G%jscB,G%jecB
      do i = G%isc,G%iec
        if (G%mask2dCv(i,j)>0.) then
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
          call remapping_core_h( nz, h1, v(i,J,:), nz, h2, u_column, CS_remapping )
          v(i,J,:) = u_column(:)
        endif
      enddo
    enddo
  endif
!$OMP end parallel

  if (show_call_tree) call callTree_waypoint("v remapped (remap_all_state_vars)")
  if (show_call_tree) call callTree_leave("remap_all_state_vars()")

end subroutine remap_all_state_vars


!> Remaps a single scalar between grids described by thicknesses h_src and h_dst.
!! h_dst must be dimensioned as a model array with GV%ke layers while h_src can
!! have an arbitrary number of layers specified by nk_src.
subroutine ALE_remap_scalar(CS, G, GV, nk_src, h_src, s_src, h_dst, s_dst, all_cells, old_remap )
  type(remapping_CS),                      intent(in)    :: CS        !< Remapping control structure
  type(ocean_grid_type),                   intent(in)    :: G         !< Ocean grid structure
  type(verticalGrid_type),                 intent(in)    :: GV        !< Ocean vertical grid structure
  integer,                                 intent(in)    :: nk_src    !< Number of levels on source grid
  real, dimension(SZI_(G),SZJ_(G),nk_src), intent(in)    :: h_src     !< Level thickness of source grid (m or Pa)
  real, dimension(SZI_(G),SZJ_(G),nk_src), intent(in)    :: s_src     !< Scalar on source grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),intent(in)    :: h_dst    !< Level thickness of destination grid (m or Pa)
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),intent(inout) :: s_dst    !< Scalar on destination grid
  logical, optional,                       intent(in)    :: all_cells !< If false, only reconstruct for
                                                                      !! non-vanished cells. Use all vanished
                                                                      !! layers otherwise (default).
  logical, optional,                       intent(in)    :: old_remap !< If true, use the old "remapping_core_w"
                                                                      !! method, otherwise use "remapping_core_h".
  ! Local variables
  integer :: i, j, k, n_points
  real :: dx(GV%ke+1)
  logical :: ignore_vanished_layers, use_remapping_core_w

  ignore_vanished_layers = .false.
  if (present(all_cells)) ignore_vanished_layers = .not. all_cells
  use_remapping_core_w = .false.
  if (present(old_remap)) use_remapping_core_w = old_remap
  n_points = nk_src

!$OMP parallel default(none) shared(CS,G,GV,h_src,s_src,h_dst,s_dst &
!$OMP                               ,ignore_vanished_layers, use_remapping_core_w, nk_src,dx ) &
!$OMP                        firstprivate(n_points)
!$OMP do
  do j = G%jsc,G%jec
    do i = G%isc,G%iec
      if (G%mask2dT(i,j)>0.) then
        if (ignore_vanished_layers) then
          n_points = 0
          do k = 1, nk_src
            if (h_src(i,j,k)>0.) n_points = n_points + 1
          enddo
          s_dst(i,j,:) = 0.
        endif
        if (use_remapping_core_w) then
          call dzFromH1H2( n_points, h_src(i,j,1:n_points), GV%ke, h_dst(i,j,:), dx )
          call remapping_core_w(CS, n_points, h_src(i,j,1:n_points), s_src(i,j,1:n_points), GV%ke, dx, s_dst(i,j,:))
        else
          call remapping_core_h(n_points, h_src(i,j,1:n_points), s_src(i,j,1:n_points), GV%ke, h_dst(i,j,:), s_dst(i,j,:), CS)
        endif
      else
        s_dst(i,j,:) = 0.
      endif
    enddo
  enddo
!$OMP end parallel

end subroutine ALE_remap_scalar


!> Use plm reconstruction for pressure gradient (determine edge values)
!! By using a PLM (limited piecewise linear method) reconstruction, this 
!! routine determines the edge values for the salinity and temperature 
!! within each layer. These edge values are returned and are used to compute 
!! the pressure gradient (by computing the densities).
subroutine pressure_gradient_plm( CS, S_t, S_b, T_t, T_b, G, GV, tv, h )
  type(ocean_grid_type),                     intent(in)    :: G    !< ocean grid structure 
  type(verticalGrid_type),                   intent(in)    :: GV   !< Ocean vertical grid structure
  type(ALE_CS),                              intent(inout) :: CS   !< module control structure 
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: S_t  !< Salinity at the top edge of each layer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: S_b  !< Salinity at the bottom edge of each layer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: T_t  !< Temperature at the top edge of each layer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: T_b  !< Temperature at the bottom edge of each layer
  type(thermo_var_ptrs),                     intent(in)    :: tv   !< thermodynamics structure 
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h    !< layer thickness 

  ! Local variables
  integer :: i, j, k
  real    :: hTmp(GV%ke)
  real    :: tmp(GV%ke)
  real, dimension(CS%nk,2)                  :: ppoly_linear_E            !Edge value of polynomial
  real, dimension(CS%nk,CS%degree_linear+1) :: ppoly_linear_coefficients !Coefficients of polynomial

  ! NOTE: the variables 'CS%grid_generic' and 'CS%ppoly_linear' are declared at
  ! the module level.

  ! Determine reconstruction within each column
!$OMP parallel do default(none) shared(G,GV,h,tv,CS,S_t,S_b,T_t,T_b)                  &
!$OMP                          private(hTmp,ppoly_linear_E,ppoly_linear_coefficients,tmp)
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1
      ! Build current grid
      hTmp(:) = h(i,j,:)*GV%H_to_m
      tmp(:) = tv%S(i,j,:)
      ! Reconstruct salinity profile    
      ppoly_linear_E = 0.0
      ppoly_linear_coefficients = 0.0
      call PLM_reconstruction( GV%ke, hTmp, tmp, ppoly_linear_E, ppoly_linear_coefficients )
      if (CS%boundary_extrapolation_for_pressure) call &
        PLM_boundary_extrapolation( GV%ke, hTmp, tmp, ppoly_linear_E, ppoly_linear_coefficients )
      
      do k = 1,GV%ke
        S_t(i,j,k) = ppoly_linear_E(k,1)
        S_b(i,j,k) = ppoly_linear_E(k,2)
      end do
      
      ! Reconstruct temperature profile 
      ppoly_linear_E = 0.0
      ppoly_linear_coefficients = 0.0
      tmp(:) = tv%T(i,j,:)
      call PLM_reconstruction( GV%ke, hTmp, tmp, ppoly_linear_E, ppoly_linear_coefficients )
      if (CS%boundary_extrapolation_for_pressure) call &
        PLM_boundary_extrapolation( GV%ke, hTmp, tmp, ppoly_linear_E, ppoly_linear_coefficients )
      
      do k = 1,GV%ke
        T_t(i,j,k) = ppoly_linear_E(k,1)
        T_b(i,j,k) = ppoly_linear_E(k,2)
      end do
      
    end do
  end do

end subroutine pressure_gradient_plm


!> Use ppm reconstruction for pressure gradient (determine edge values)
!> By using a PPM (limited piecewise linear method) reconstruction, this 
!> routine determines the edge values for the salinity and temperature 
!> within each layer. These edge values are returned and are used to compute 
!> the pressure gradient (by computing the densities).
subroutine pressure_gradient_ppm( CS, S_t, S_b, T_t, T_b, G, GV, tv, h )
  type(ocean_grid_type),                     intent(in)    :: G    !< ocean grid structure 
  type(verticalGrid_type),                   intent(in)    :: GV   !< Ocean vertical grid structure
  type(ALE_CS),                              intent(inout) :: CS   !< module control structure 
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: S_t  !< Salinity at top edge of each layer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: S_b  !< Salinity at bottom edge of each layer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: T_t  !< Temperature at the top edge of each layer
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: T_b  !< Temperature at the bottom edge of each layer
  type(thermo_var_ptrs),                     intent(in)    :: tv   !< ocean thermodynamics structure 
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h    !< layer thickness 

  ! Local variables
  integer :: i, j, k
  real    :: hTmp(GV%ke)
  real    :: tmp(GV%ke)
  real, dimension(CS%nk,2) :: &
      ppoly_parab_E            !Edge value of polynomial
  real, dimension(CS%nk,CS%degree_parab+1) :: &
      ppoly_parab_coefficients !Coefficients of polynomial


  ! NOTE: the variables 'CS%grid_generic' and 'CS%ppoly_parab' are declared at
  ! the module level.

  ! Determine reconstruction within each column
!$OMP parallel do default(none) shared(G,GV,h,tv,CS,S_t,S_b,T_t,T_b) &
!$OMP                          private(hTmp,tmp,ppoly_parab_E,ppoly_parab_coefficients)
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1
     
      ! Build current grid
      hTmp(:) = h(i,j,:) * GV%H_to_m
      tmp(:) = tv%S(i,j,:)
      
      ! Reconstruct salinity profile    
      ppoly_parab_E = 0.0
      ppoly_parab_coefficients = 0.0
      call edge_values_implicit_h4( GV%ke, hTmp, tmp, ppoly_parab_E )
      call PPM_reconstruction( GV%ke, hTmp, tmp, ppoly_parab_E, ppoly_parab_coefficients )
      if (CS%boundary_extrapolation_for_pressure) call &
        PPM_boundary_extrapolation( GV%ke, hTmp, tmp, ppoly_parab_E, ppoly_parab_coefficients )
      
      do k = 1,GV%ke
        S_t(i,j,k) = ppoly_parab_E(k,1)
        S_b(i,j,k) = ppoly_parab_E(k,2)
      end do
      
      ! Reconstruct temperature profile 
      ppoly_parab_E = 0.0
      ppoly_parab_coefficients = 0.0
      tmp(:) = tv%T(i,j,:)
      call edge_values_implicit_h4( GV%ke, hTmp, tmp, ppoly_parab_E )
      call PPM_reconstruction( GV%ke, hTmp, tmp, ppoly_parab_E, ppoly_parab_coefficients )
      if (CS%boundary_extrapolation_for_pressure) call &
        PPM_boundary_extrapolation( GV%ke, hTmp, tmp, ppoly_parab_E, ppoly_parab_coefficients )
      
      do k = 1,GV%ke
        T_t(i,j,k) = ppoly_parab_E(k,1)
        T_b(i,j,k) = ppoly_parab_E(k,2)
      end do
      
    end do
  end do

end subroutine pressure_gradient_ppm


!> pressure reconstruction logical 
logical function usePressureReconstruction(CS)
  type(ALE_CS), pointer :: CS  !< control structure 

  if (associated(CS)) then
    usePressureReconstruction=CS%reconstructForPressure
  else
    usePressureReconstruction=.false.
  endif

end function usePressureReconstruction


!> pressure reconstruction integer 
integer function pressureReconstructionScheme(CS)
  type(ALE_CS), pointer :: CS !< control structure 

  if (associated(CS)) then
    pressureReconstructionScheme=CS%pressureReconstructionScheme
  else
    pressureReconstructionScheme=-1
  endif

end function pressureReconstructionScheme


!> Initialize regridding module 
subroutine ALE_initRegridding(GV, max_depth, param_file, mod, regridCS, dz )
  type(verticalGrid_type), intent(in)  :: GV         !< Ocean vertical grid structure
  real,                    intent(in)  :: max_depth  !< The maximum depth of the ocean, in m.
  type(param_file_type),   intent(in)  :: param_file !< parameter file 
  character(len=*),        intent(in)  :: mod        !< Name of calling module
  type(regridding_CS),     intent(out) :: regridCS   !< Regridding parameters and work arrays
  real, dimension(:),      intent(out) :: dz         !< Resolution (thickness) in units of coordinate

  ! Local variables
  character(len=80)  :: string, varName ! Temporary strings
  character(len=40)  :: coordMode, interpScheme, coordUnits ! Temporary strings
  character(len=200) :: inputdir, fileName
  character(len=320) :: message ! Temporary strings
  integer :: K, ke
  logical :: tmpLogical, fix_haloclines, set_max, do_sum
  real :: filt_len, strat_tol, index_scale
  real :: tmpReal, compress_fraction
  real :: dz_fixed_sfc, Rho_avg_depth, nlay_sfc_int
  integer :: nz_fixed_sfc
  real :: rho_target(GV%ke+1) ! Target density used in HYBRID mode
  real, dimension(size(dz))   :: h_max  ! Maximum layer thicknesses, in m.
  real, dimension(size(dz))   :: dz_max ! Thicknesses used to find maximum interface depths, in m.
  real, dimension(size(dz)+1) :: z_max  ! Maximum tinterface depths, in m.

  ke = size(dz) ! Number of levels in resolution vector

  call get_param(param_file, mod, "REGRIDDING_COORDINATE_MODE", coordMode, &
                 "Coordinate mode for vertical regridding.\n"//&
                 "Choose among the following possibilities:\n"//&
                 trim(regriddingCoordinateModeDoc),&
                 default=DEFAULT_COORDINATE_MODE, fail_if_missing=.true.)
  call get_param(param_file, mod, "REGRIDDING_COORDINATE_UNITS", coordUnits, &
                 "Units of the regridding coordinuate.",&
                 default=coordinateUnits(coordMode))

  call get_param(param_file, mod, "INTERPOLATION_SCHEME", interpScheme, &
                 "This sets the interpolation scheme to use to\n"//&
                 "determine the new grid. These parameters are\n"//&
                 "only relevant when REGRIDDING_COORDINATE_MODE is\n"//&
                 "set to a function of state. Otherwise, it is not\n"//&
                 "used. It can be one of the following schemes:\n"//&
                 trim(regriddingInterpSchemeDoc),&
                 default=regriddingDefaultInterpScheme)

  call get_param(param_file, mod, "REGRID_COMPRESSIBILITY_FRACTION", compress_fraction, &
                 "When interpolating potential density profiles we can add\n"//&
                 "some artificial compressibility solely to make homogenous\n"//&
                 "regions appear stratified.", default=0.)

  call initialize_regridding( GV%ke, coordMode, interpScheme, regridCS, &
                              compressibility_fraction=compress_fraction )

  call get_param(param_file, mod, "ALE_COORDINATE_CONFIG", string, &
                 "Determines how to specify the coordinate\n"//&
                 "resolution. Valid options are:\n"//&
                 " PARAM       - use the vector-parameter ALE_RESOLUTION\n"//&
                 " UNIFORM     - uniformly distributed\n"//&
                 " FILE:string - read from a file. The string specifies\n"//&
                 "               the filename and variable name, separated\n"//&
                 "               by a comma or space, e.g. FILE:lev.nc,Z\n"//&
                 " FNC1:string - FNC1:dz_min,H_total,power,precision\n"//&
                 " HYBRID:string - read from a file. The string specifies\n"//&
                 "               the filename and two variable names, separated\n"//&
                 "               by a comma or space, for sigma-2 and dz. e.g.\n"//&
                 "               HYBRID:vgrid.nc,sigma2,dz",&
                 default='UNIFORM')
  message = "The distribution of vertical resolution for the target\n"//&
            "grid used for Eulerian-like coordinates. For example,\n"//&
            "in z-coordinate mode, the parameter is a list of level\n"//&
            "thicknesses (in m). In sigma-coordinate mode, the list\n"//&
            "is of non-dimensional fractions of the water column."
  select case ( trim(string) )
    case ("UNIFORM")
      dz(:) = uniformResolution(GV%ke, coordMode, max_depth, &
                 GV%Rlay(1)+0.5*(GV%Rlay(1)-GV%Rlay(2)), &
                 GV%Rlay(GV%ke)+0.5*(GV%Rlay(GV%ke)-GV%Rlay(GV%ke-1)) )
      call log_param(param_file, mod, "!ALE_RESOLUTION", dz, &
                   trim(message), units=trim(coordUnits))
    case ("PARAM")
      call get_param(param_file, mod, "ALE_RESOLUTION", dz, &
                   trim(message), units=trim(coordUnits), fail_if_missing=.true.)
    case default 
      if (index(trim(string),'FILE:')==1) then
        call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
        inputdir = slasher(inputdir)

        if (string(6:6)=='.' .or. string(6:6)=='/') then
          ! If we specified "FILE:./xyz" or "FILE:/xyz" then we have a relative or absolute path
          fileName = trim( extractWord(trim(string(6:80)), 1) )
        else
          ! Otherwise assume we should look for the file in INPUTDIR
          fileName = trim(inputdir) // trim( extractWord(trim(string(6:80)), 1) )
        endif
        if (.not. file_exists(fileName)) call MOM_error(FATAL,"ALE_initRegridding: "// &
          "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")

        varName = trim( extractWord(trim(string(6:)), 2) )
        if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,"ALE_initRegridding: "// &
          "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
        if (len_trim(varName)==0) then
          if (field_exists(fileName,'dz')) then; varName = 'dz'
          elseif (field_exists(fileName,'dsigma')) then; varName = 'dsigma'
          elseif (field_exists(fileName,'ztest')) then; varName = 'ztest'
          else ;  call MOM_error(FATAL,"ALE_initRegridding: "// &
            "Coordinate variable not specified and none could be guessed.")
          endif
        endif
        call MOM_read_data(trim(fileName), trim(varName), dz)
        call log_param(param_file, mod, "!ALE_RESOLUTION", dz, &
                   trim(message), units=coordinateUnits(coordMode))
      elseif (index(trim(string),'FNC1:')==1) then
        call dz_function1( trim(string(6:)), dz )
        call log_param(param_file, mod, "!ALE_RESOLUTION", dz, &
                   trim(message), units=coordinateUnits(coordMode))
      elseif (index(trim(string),'HYBRID:')==1) then
        call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
        inputdir = slasher(inputdir)

        ! The following assumes the FILE: syntax of above but without "FILE:" in the string
        fileName = trim( extractWord(trim(string(8:)), 1) )
        if (fileName(1:1)/='.' .and. filename(1:1)/='/') fileName = trim(inputdir) // trim( fileName )
        if (.not. file_exists(fileName)) call MOM_error(FATAL,"ALE_initRegridding: HYBRID "// &
          "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")
        varName = trim( extractWord(trim(string(8:)), 2) )
        if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,"ALE_initRegridding: HYBRID "// &
          "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
        call MOM_read_data(trim(fileName), trim(varName), rho_target)
        call set_target_densities( regridCS, rho_target )
        varName = trim( extractWord(trim(string(8:)), 3) )
        if (varName(1:5) == 'FNC1:') then ! Use FNC1 to calculate dz
          call dz_function1( trim(string((index(trim(string),'FNC1:')+5):)), dz )
        else ! Read dz from file
          if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,"ALE_initRegridding: HYBRID "// &
            "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
          call MOM_read_data(trim(fileName), trim(varName), dz)
        endif
        call log_param(param_file, mod, "!ALE_RESOLUTION", dz, &
                   trim(message), units=coordinateUnits(coordMode))
        call log_param(param_file, mod, "!TARGET_DENSITIES", rho_target, &
                   'HYBRID target densities for itnerfaces', units=coordinateUnits(coordMode))

      else
        call MOM_error(FATAL,"ALE_initRegridding: "// &
          "Unrecognized coordinate configuraiton"//trim(string))
      endif
  end select
  if (coordinateMode(coordMode) == REGRIDDING_ZSTAR .or. &
      coordinateMode(coordMode) == REGRIDDING_HYCOM1 .or. &
      coordinateMode(coordMode) == REGRIDDING_SLIGHT) then
    ! Adjust target grid to be consistent with max_depth
    ! This is a work around to the from_Z initialization...  ???
    tmpReal = sum( dz(:) )
    if (tmpReal < max_depth) then
      dz(ke) = dz(ke) + ( max_depth - tmpReal )
    elseif (tmpReal > max_depth) then
      if ( dz(ke) + ( max_depth - tmpReal ) > 0. ) then
        dz(ke) = dz(ke) + ( max_depth - tmpReal )
      else
        call MOM_error(FATAL,"ALE_initRegridding: "// &
          "MAXIMUM_DEPTH was too shallow to adjust bottom layer of DZ!"//trim(string))
      endif
    endif
  endif
  call setCoordinateResolution( dz, regridCS )
  if (coordinateMode(coordMode) == REGRIDDING_RHO) call set_target_densities_from_GV( GV, regridCS )

  call get_param(param_file, mod, "MIN_THICKNESS", tmpReal, &
                 "When regridding, this is the minimum layer\n"//&
                 "thickness allowed.", units="m",&
                 default=regriddingDefaultMinThickness )

  call get_param(param_file, mod, "BOUNDARY_EXTRAPOLATION", tmpLogical, &
                 "When defined, a proper high-order reconstruction\n"//&
                 "scheme is used within boundary cells rather\n"//&
                 "than PCM. E.g., if PPM is used for remapping, a\n"//&
                 "PPM reconstruction will also be used within\n"//&
                 "boundary cells.", default=regriddingDefaultBoundaryExtrapolation)
  call set_regrid_params( regridCS, min_thickness=tmpReal, boundary_extrapolation=tmpLogical )

  if (coordinateMode(coordMode) == REGRIDDING_SLIGHT) then
    ! Set SLight-specific regridding parameters.
    call get_param(param_file, mod, "SLIGHT_DZ_SURFACE", dz_fixed_sfc, &
                 "The nominal thickness of fixed thickness near-surface\n"//&
                 "layers with the SLight coordinate.", units="m", default=1.0)
    call get_param(param_file, mod, "SLIGHT_NZ_SURFACE_FIXED", nz_fixed_sfc, &
                 "The number of fixed-depth surface layers with the SLight\n"//&
                 "coordinate.", units="nondimensional", default=2)
    call get_param(param_file, mod, "SLIGHT_SURFACE_AVG_DEPTH", Rho_avg_depth, &
                 "The thickness of the surface region over which to average\n"//&
                 "when calculating the density to use to define the interior\n"//&
                 "with the SLight coordinate.", units="m", default=1.0)
    call get_param(param_file, mod, "SLIGHT_NLAY_TO_INTERIOR", nlay_sfc_int, &
                 "The number of layers to offset the surface density when\n"//&
                 "defining where the interior ocean starts with SLight.", &
                 units="nondimensional", default=2.0)
    call get_param(param_file, mod, "SLIGHT_FIX_HALOCLINES", fix_haloclines, &
                 "If true, identify regions above the reference pressure\n"//&
                 "where the reference pressure systematically underestimates\n"//&
                 "the stratification and use this in the definition of the\n"//&
                 "interior with the SLight coordinate.", default=.false.)
                 
    call set_regrid_params( regridCS, dz_min_surface=dz_fixed_sfc, &
                nz_fixed_surface=nz_fixed_sfc, Rho_ML_avg_depth=Rho_avg_depth, &
                nlay_ML_to_interior=nlay_sfc_int, fix_haloclines=fix_haloclines)
    if (fix_haloclines) then
      ! Set additional parameters related to SLIGHT_FIX_HALOCLINES.
      call get_param(param_file, mod, "HALOCLINE_FILTER_LENGTH", filt_len, &
                 "A length scale over which to smooth the temperature and\n"//&
                 "salinity before identifying erroneously unstable haloclines.", &
                 units="m", default=2.0)
      call get_param(param_file, mod, "HALOCLINE_STRAT_TOL", strat_tol, &
                 "A tolerance for the ratio of the stratification of the\n"//&
                 "apparent coordinate stratification to the actual value\n"//&
                 "that is used to identify erroneously unstable haloclines.\n"//&
                 "This ratio is 1 when they are equal, and sensible values \n"//&
                 "are between 0 and 0.5.", units="nondimensional", default=0.2)
      call set_regrid_params(regridCS, halocline_filt_len=filt_len, &
                             halocline_strat_tol=strat_tol)
    endif

  endif

  call get_param(param_file, mod, "MAXIMUM_INT_DEPTH_CONFIG", string, &
                 "Determines how to specify the maximum interface depths.\n"//&
                 "Valid options are:\n"//&
                 " NONE        - there are no maximum interface depths\n"//&
                 " PARAM       - use the vector-parameter MAXIMUM_INTERFACE_DEPTHS\n"//&
                 " FILE:string - read from a file. The string specifies\n"//&
                 "               the filename and variable name, separated\n"//&
                 "               by a comma or space, e.g. FILE:lev.nc,Z\n"//&
                 " FNC1:string - FNC1:dz_min,H_total,power,precision",&
                 default='NONE')
  message = "The list of maximum depths for each interface."
  if ( trim(string) == "NONE") then
    ! Do nothing.
  elseif ( trim(string) ==  "PARAM") then
    call get_param(param_file, mod, "MAXIMUM_INTERFACE_DEPTHS", z_max, &
                 trim(message), units="m", fail_if_missing=.true.)
    call set_regrid_max_depths( regridCS, z_max, GV%m_to_H )
  elseif (index(trim(string),'FILE:')==1) then
    call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)

    if (string(6:6)=='.' .or. string(6:6)=='/') then
      ! If we specified "FILE:./xyz" or "FILE:/xyz" then we have a relative or absolute path
      fileName = trim( extractWord(trim(string(6:80)), 1) )
    else
      ! Otherwise assume we should look for the file in INPUTDIR
      fileName = trim(inputdir) // trim( extractWord(trim(string(6:80)), 1) )
    endif
    if (.not. file_exists(fileName)) call MOM_error(FATAL,"ALE_initRegridding: "// &
      "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")

    do_sum = .false.
    varName = trim( extractWord(trim(string(6:)), 2) )
    if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,"ALE_initRegridding: "// &
      "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
    if (len_trim(varName)==0) then
      if (field_exists(fileName,'z_max')) then; varName = 'z_max'
      elseif (field_exists(fileName,'dz')) then; varName = 'dz' ; do_sum = .true.
      elseif (field_exists(fileName,'dz_max')) then; varName = 'dz_max' ; do_sum = .true.
      else ; call MOM_error(FATAL,"ALE_initRegridding: "// &
        "MAXIMUM_INT_DEPTHS variable not specified and none could be guessed.")
      endif
    endif
    if (do_sum) then
      call MOM_read_data(trim(fileName), trim(varName), dz_max)
      z_max(1) = 0.0 ; do K=1,ke ; z_max(K+1) = z_max(K) + dz_max(k) ; enddo
    else
      call MOM_read_data(trim(fileName), trim(varName), z_max)
    endif
    call log_param(param_file, mod, "!MAXIMUM_INT_DEPTHS", z_max, &
               trim(message), units=coordinateUnits(coordMode))
    call set_regrid_max_depths( regridCS, z_max, GV%m_to_H )
  elseif (index(trim(string),'FNC1:')==1) then
    call dz_function1( trim(string(6:)), dz_max )
    if ((coordinateMode(coordMode) == REGRIDDING_SLIGHT) .and. &
        (dz_fixed_sfc > 0.0)) then
      do k=1,nz_fixed_sfc ; dz_max(k) = dz_fixed_sfc ; enddo
    endif
    z_max(1) = 0.0 ; do K=1,ke ; z_max(K+1) = z_max(K) + dz_max(K) ; enddo
    call log_param(param_file, mod, "!MAXIMUM_INT_DEPTHS", z_max, &
               trim(message), units=coordinateUnits(coordMode))
    call set_regrid_max_depths( regridCS, z_max, GV%m_to_H )
  else
    call MOM_error(FATAL,"ALE_initRegridding: "// &
      "Unrecognized MAXIMUM_INT_DEPTH_CONFIG "//trim(string))
  endif

  ! Optionally specify maximum thicknesses for each layer, enforced by moving
  ! the interface below a layer downward.
  call get_param(param_file, mod, "MAX_LAYER_THICKNESS_CONFIG", string, &
                 "Determines how to specify the maximum layer thicknesses.\n"//&
                 "Valid options are:\n"//&
                 " NONE        - there are no maximum layer thicknesses\n"//&
                 " PARAM       - use the vector-parameter MAX_LAYER_THICKNESS\n"//&
                 " FILE:string - read from a file. The string specifies\n"//&
                 "               the filename and variable name, separated\n"//&
                 "               by a comma or space, e.g. FILE:lev.nc,Z\n"//&
                 " FNC1:string - FNC1:dz_min,H_total,power,precision",&
                 default='NONE')
  message = "The list of maximum thickness for each layer."
  if ( trim(string) == "NONE") then
    ! Do nothing.
  elseif ( trim(string) ==  "PARAM") then
    call get_param(param_file, mod, "MAX_LAYER_THICKNESS", h_max, &
                 trim(message), units="m", fail_if_missing=.true.)
    call set_regrid_max_thickness( regridCS, h_max, GV%m_to_H )
  elseif (index(trim(string),'FILE:')==1) then
    call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)

    if (string(6:6)=='.' .or. string(6:6)=='/') then
      ! If we specified "FILE:./xyz" or "FILE:/xyz" then we have a relative or absolute path
      fileName = trim( extractWord(trim(string(6:80)), 1) )
    else
      ! Otherwise assume we should look for the file in INPUTDIR
      fileName = trim(inputdir) // trim( extractWord(trim(string(6:80)), 1) )
    endif
    if (.not. file_exists(fileName)) call MOM_error(FATAL,"ALE_initRegridding: "// &
      "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")

    varName = trim( extractWord(trim(string(6:)), 2) )
    if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,"ALE_initRegridding: "// &
      "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
    if (len_trim(varName)==0) then
      if (field_exists(fileName,'h_max')) then; varName = 'h_max'
      elseif (field_exists(fileName,'dz_max')) then; varName = 'dz_max'
      else ; call MOM_error(FATAL,"ALE_initRegridding: "// &
        "MAXIMUM_INT_DEPTHS variable not specified and none could be guessed.")
      endif
    endif
    call MOM_read_data(trim(fileName), trim(varName), h_max)
    call log_param(param_file, mod, "!MAX_LAYER_THICKNESS", h_max, &
               trim(message), units=coordinateUnits(coordMode))
    call set_regrid_max_thickness( regridCS, h_max, GV%m_to_H )
  elseif (index(trim(string),'FNC1:')==1) then
    call dz_function1( trim(string(6:)), h_max )
    call log_param(param_file, mod, "!MAX_LAYER_THICKNESS", h_max, &
               trim(message), units=coordinateUnits(coordMode))
    call set_regrid_max_thickness( regridCS, h_max, GV%m_to_H )
  else
    call MOM_error(FATAL,"ALE_initRegridding: "// &
      "Unrecognized MAX_LAYER_THICKNESS_CONFIG "//trim(string))
  endif

end subroutine ALE_initRegridding


!> Parses a string and generates a dz(:) profile
subroutine dz_function1( string, dz )
  character(len=*),   intent(in)    :: string !< String with list of parameters
  real, dimension(:), intent(inout) :: dz     !< Profile of nominal thicknesses

  ! Local variables
  integer :: nk, k
  real    :: dz_min, power, prec, H_total

  nk = size(dz) ! Number of cells
  prec = -1024.
  read( string, *) dz_min, H_total, power, prec
  if (prec == -1024.) call MOM_error(FATAL,"dz_function1: "// &
          "Problem reading FNC1: string  ="//trim(string))
  ! Create profile of ( dz - dz_min )
  do k = 1, nk
    dz(k) = (real(k-1)/real(nk-1))**power
  enddo
  dz(:) = ( H_total - real(nk) * dz_min ) * ( dz(:) / sum(dz) ) ! Rescale to so total is H_total
  dz(:) = anint( dz(:) / prec ) * prec ! Rounds to precision prec
  dz(:) = ( H_total - real(nk) * dz_min ) * ( dz(:) / sum(dz) ) ! Rescale to so total is H_total
  dz(:) = anint( dz(:) / prec ) * prec ! Rounds to precision prec
  dz(nk) = dz(nk) + ( H_total - sum( dz(:) + dz_min ) ) ! Adjust bottommost layer
  dz(:) = anint( dz(:) / prec ) * prec ! Rounds to precision prec
  dz(:) = dz(:) + dz_min ! Finally add in the constant dz_min

end subroutine dz_function1


!> Query the target coordinate interfaces positions
function ALE_getCoordinate( CS )
  type(ALE_CS), pointer    :: CS                  !< module control structure 

  real, dimension(CS%nk+1) :: ALE_getCoordinate  
  ALE_getCoordinate(:) = getCoordinateInterfaces( CS%regridCS )

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
  real,         intent(in) :: dt !< Time-step used between ALE calls
  type(ALE_CS), pointer    :: CS !< ALE control structure
  ! Local variables
  real :: w  ! An implicit weighting estimate.

  if (associated(CS)) then
    w = 0.0
    if (CS%regrid_time_scale > 0.0) then
      w = CS%regrid_time_scale / (CS%regrid_time_scale + dt)
    endif
    call set_regrid_params( CS%regridCS, old_grid_weight=w )
  endif

end subroutine ALE_update_regrid_weights

!> Update the vertical grid type with ALE information.
!! This subroutine sets information in the verticalGrid_type to be
!! consistent with the use of ALE mode.
subroutine ALE_updateVerticalGridType( CS, GV )
  type(ALE_CS),            pointer :: CS  ! module control structure 
  type(verticalGrid_type), pointer :: GV  ! vertical grid information 

  integer :: nk

  nk = GV%ke
  GV%sInterface(1:nk+1) = getCoordinateInterfaces( CS%regridCS )
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
  ds(:)       = getCoordinateResolution( CS%regridCS )
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
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: h   !< layer thickness 

  ! Local variables
  integer :: i, j, k

  do j = G%jsd,G%jed ; do i = G%isd,G%ied
    h(i,j,:) = getStaticThickness( CS%regridCS, 0., G%bathyT(i,j) )
  enddo; enddo

end subroutine ALE_initThicknessToCoord

end module MOM_ALE
