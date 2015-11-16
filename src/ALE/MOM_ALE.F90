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
use MOM_regridding,       only : initialize_regridding, regridding_main , end_regridding
use MOM_regridding,       only : uniformResolution, set_old_grid_weight
use MOM_regridding,       only : inflate_vanished_layers_old, setCoordinateResolution
use MOM_regridding,       only : set_target_densities_from_G, set_target_densities
use MOM_regridding,       only : regriddingCoordinateModeDoc, DEFAULT_COORDINATE_MODE
use MOM_regridding,       only : regriddingInterpSchemeDoc, regriddingDefaultInterpScheme
use MOM_regridding,       only : setRegriddingBoundaryExtrapolation
use MOM_regridding,       only : regriddingDefaultBoundaryExtrapolation
use MOM_regridding,       only : set_regrid_min_thickness, regriddingDefaultMinThickness
use MOM_regridding,       only : check_remapping_grid
use MOM_regridding,       only : regridding_CS
use MOM_regridding,       only : getCoordinateInterfaces, getCoordinateResolution
use MOM_regridding,       only : getCoordinateUnits, getCoordinateShortName
use MOM_regridding,       only : getStaticThickness
use MOM_remapping,        only : initialize_remapping, remapping_core, end_remapping
use MOM_remapping,        only : remappingSchemesDoc, remappingDefaultScheme
use MOM_remapping,        only : remapDisableBoundaryExtrapolation, remapEnableBoundaryExtrapolation
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

  logical :: boundary_extrapolation_for_pressure !<  Indicate whether high-order boundary 
                                                 !!  extrapolation should be used within boundary cells

  logical :: reconstructForPressure = .false.    !< Indicates whether integrals for FV 
                                                 !! pressure gradient calculation will
                                                 !! use reconstruction of T/S.
                                                 !! By default, it is true if regridding 
                                                 !! has been initialized, otherwise false.

  integer :: pressureReconstructionScheme        !<  Form of the reconstruction of T/S 
                                                 !! for FV pressure gradient calculation.
                                                 !! By default, it is =1 (PLM)

  real :: regrid_time_scale !< The time-scale used in blending between the current (old) grid
                            !! and the target (new) grid. (s)

  type(regridding_CS) :: regridCS !< Regridding parameters and work arrays
  type(remapping_CS)  :: remapCS  !< Remapping parameters and work arrays

  integer :: nk              !< Used only for queries, not directly by this module
  integer :: degree_linear=1 !< Degree of linear piecewise polynomial
  integer :: degree_parab=2  !< Degree of parabolic piecewise polynomial

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NK_INTERFACE_) :: dzRegrid  !< Work space for communicating
                                                                       !! between regridding and remapping

  logical :: remap_after_initialization !<   Indicates whether to regrid/remap after initializing the state.

  logical :: show_call_tree !< For debugging
  real    :: C_p            !< seawater heat capacity (J/(kg deg C))

  ! for diagnostics 
  type(diag_ctrl), pointer           :: diag                          !< structure to regulate output
  integer, dimension(:), allocatable :: id_tracer_remap_tendency      !< diagnostic id 
  integer, dimension(:), allocatable :: id_Htracer_remap_tendency     !< diagnostic id 
  integer, dimension(:), allocatable :: id_Htracer_remap_tendency_2d  !< diagnostic id 
  logical, dimension(:), allocatable :: do_tendency_diag              !< flag for doing diagnostics 

end type


public initialize_ALE
public end_ALE
public ALE_main 
public regrid_only
public regrid_remap_T_S
public remap_scalar_h_to_h
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
public check_remapping_grid
public remap_init_conds
public register_diags_ALE 

! List of interpolation schemes
integer, parameter :: INTERPOLATION_P1M_H2     = 0 ! O(h^2)
integer, parameter :: INTERPOLATION_P1M_H4     = 1 ! O(h^2)
integer, parameter :: INTERPOLATION_P1M_IH4    = 2 ! O(h^2)
integer, parameter :: INTERPOLATION_PLM        = 3 ! O(h^2)
integer, parameter :: INTERPOLATION_PPM_H4     = 4 ! O(h^3)
integer, parameter :: INTERPOLATION_PPM_IH4    = 5 ! O(h^3)
integer, parameter :: INTERPOLATION_P3M_IH4IH3 = 6 ! O(h^4)
integer, parameter :: INTERPOLATION_P3M_IH6IH5 = 7 ! O(h^4)
integer, parameter :: INTERPOLATION_PQM_IH4IH3 = 8 ! O(h^4)
integer, parameter :: INTERPOLATION_PQM_IH6IH5 = 9 ! O(h^5)

contains


!> This routine is typically called (from initialize_MOM in file MOM.F90)
!! before the main time integration loop to initialize the regridding stuff.
!! We read the MOM_input file to register the values of different
!! regridding/remapping parameters.
subroutine initialize_ALE( param_file, G, CS)
  type(param_file_type),   intent(in) :: param_file !< Parameter file
  type(ocean_grid_type),   intent(in) :: G          !< Ocean grid structure
  type(ALE_CS),            pointer    :: CS         !< Module control structure

  ! Local variables
  real, dimension(:), allocatable :: dz
  character(len=40)               :: mod = "MOM_ALE" ! This module's name.
  character(len=80)               :: string ! Temporary strings

  if (associated(CS)) then
    call MOM_error(WARNING, "initialize_ALE called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%show_call_tree = callTree_showQuery()
  if (CS%show_call_tree) call callTree_enter("initialize_ALE(), MOM_ALE.F90")

  ! Memory allocation for regridding
  call ALE_memory_allocation( G, CS )

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

  ! Initialize and configure regridding
  allocate( dz(G%ke) )
  call ALE_initRegridding( G, param_file, mod, CS%regridCS, dz )
  deallocate( dz )

  ! Initialize and configure remapping
  call get_param(param_file, mod, "REMAPPING_SCHEME", string, &
                 "This sets the reconstruction scheme used\n"//&
                 "for vertical remapping for all variables.\n"//&
                 "It can be one of the following schemes:\n"//&
                 trim(remappingSchemesDoc), default=remappingDefaultScheme)
  call initialize_remapping( G%ke, string, CS%remapCS )
  call remapDisableBoundaryExtrapolation( CS%remapCS )

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

  ! Keep a record of values for subsequent queries
  CS%nk = G%ke

  if (CS%show_call_tree) call callTree_leave("initialize_ALE()")
end subroutine initialize_ALE


!> Initialize diagnostics for the ALE module. 
subroutine register_diags_ALE(Time, G, diag, C_p, Reg, CS)
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


end subroutine register_diags_ALE


!> Crudely adjust (initial) grid for integrity.
!! This routine is typically called (from initialize_MOM in file MOM.F90)
!! before the main time integration loop to initialize the regridding stuff.
!! We read the MOM_input file to register the values of different
!! regridding/remapping parameters.
subroutine adjustGridForIntegrity( CS, G, h )
  type(ALE_CS), pointer                                 :: CS
  type(ocean_grid_type), intent(in)                     :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout) :: h

  call inflate_vanished_layers_old( CS%regridCS, G, h(:,:,:) )

end subroutine adjustGridForIntegrity


!> End of regridding (memory deallocation).
!! This routine is typically called (from MOM_end in file MOM.F90)
!! after the main time integration loop to deallocate the regridding stuff.
subroutine end_ALE(CS)
  type(ALE_CS), pointer :: CS  !< module control structure 
  
  ! Deallocate memory used for the regridding
  call end_remapping( CS%remapCS )
  call end_regridding( CS%regridCS )
  call ALE_memory_deallocation( CS )

  deallocate(CS)

end subroutine end_ALE

!> Takes care of (1) building a new grid and (2) remapping all variables between
!! the old grid and the new grid. The creation of the new grid can be based
!! on z coordinates, target interface densities, sigma coordinates or any
!! arbitrary coordinate system.
subroutine ALE_main( G, h, u, v, tv, Reg, CS, dt)
  type(ocean_grid_type),                   intent(in)    :: G   !< Ocean grid informations
  real, dimension(NIMEM_,NJMEM_, NKMEM_),  intent(inout) :: h   !< Current 3D grid obtained after last time step (m or Pa)
  real, dimension(NIMEMB_,NJMEM_, NKMEM_), intent(inout) :: u   !< Zonal velocity field (m/s)
  real, dimension(NIMEM_,NJMEMB_, NKMEM_), intent(inout) :: v   !< Meridional velocity field (m/s)
  type(thermo_var_ptrs),                   intent(inout) :: tv  !< Thermodynamic variable structure
  type(tracer_registry_type),              pointer       :: Reg !< Tracer registry structure
  type(ALE_CS),                            pointer       :: CS  !< Regridding parameters and options
  real,                          optional, intent(in)    :: dt  !< Time step between calls to ALE_main()

  ! Local variables
  integer :: nk, i, j, k, isd, ied, jsd, jed

  if (CS%show_call_tree) call callTree_enter("ALE_main(), MOM_ALE.F90")

  if (present(dt)) then
    call ALE_update_regrid_weights( dt, CS )
  endif

  ! Build new grid. The new grid is stored in h_new. The old grid is h.
  ! Both are needed for the subsequent remapping of variables.
  call regridding_main( CS%remapCS, CS%regridCS, G, h, tv, CS%dzRegrid )

  call check_remapping_grid( G, h, CS%dzRegrid, 'in ALE_main()' )

  if (CS%show_call_tree) call callTree_waypoint("new grid generated (ALE_main)")

  ! Remap all variables from old grid h onto new grid h_new
  call remapping_main( CS%remapCS, CS, G, h, -CS%dzRegrid, Reg, u, v, CS%show_call_tree, dt)

  if (CS%show_call_tree) call callTree_waypoint("state remapped (ALE_main)")

  ! Override old grid with new one. The new grid 'h_new' is built in
  ! one of the 'build_...' routines above.
  nk = G%ke; isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,nk,h,CS)
  do k = 1,nk
    do j = jsd,jed ; do i = isd,ied
      h(i,j,k) = h(i,j,k) + ( CS%dzRegrid(i,j,k) - CS%dzRegrid(i,j,k+1) )
     enddo ; enddo
  enddo

  if (CS%show_call_tree) call callTree_leave("ALE_main()")
end subroutine ALE_main


!> Generates new grid
subroutine regrid_only( G, regridCS, remapCS, h, tv, debug )
  type(ocean_grid_type),                   intent(in)    :: G        !< Ocean grid structure 
  type(regridding_CS),                     intent(in)    :: regridCS !< Regridding parameters and options
  type(remapping_CS),                      intent(in)    :: remapCS  !< Remapping parameters and options
  type(thermo_var_ptrs),                   intent(inout) :: tv       !< Thermodynamical variable structure
  real, dimension(NIMEM_,NJMEM_, NKMEM_),  intent(inout) :: h        !< Current 3D grid obtained after the last time step (m or Pa)
  logical,                       optional, intent(in)    :: debug    !< If true, show the call tree

  ! Local variables
  integer :: nk, i, j, k
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)+1) :: dzRegrid ! The changein grid interface positions
  logical :: show_call_tree

  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug
  if (show_call_tree) call callTree_enter("regrid_only(), MOM_ALE.F90")

  ! Build new grid. The new grid is stored in h_new. The old grid is h.
  ! Both are needed for the subsequent remapping of variables.
  call regridding_main( remapCS, regridCS, G, h, tv, dzRegrid )

  call check_remapping_grid( G, h, dzRegrid, 'in regrid_only()' )

  ! Override old grid with new one. The new grid 'h_new' is built in
  ! one of the 'build_...' routines above.
!$OMP parallel do default(none) shared(G,h,dzRegrid)
  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    if (G%mask2dT(i,j)>0.) then
      do k = 1,G%ke
        h(i,j,k) = h(i,j,k) + ( dzRegrid(i,j,k) - dzRegrid(i,j,k+1) )
      enddo
    endif
  enddo ; enddo

  if (show_call_tree) call callTree_leave("regrid_only()")
end subroutine regrid_only


!> Generates new grid and remaps T and S
subroutine regrid_remap_T_S( G, regridCS, remapCS, h, tv, debug )
  type(ocean_grid_type),                   intent(in)    :: G        !< Ocean grid informations
  type(regridding_CS),                     intent(in)    :: regridCS !< Regridding parameters and options
  type(remapping_CS),                      intent(in)    :: remapCS  !< Remapping parameters and options
  real, dimension(NIMEM_,NJMEM_, NKMEM_),  intent(inout) :: h        !< Current 3D grid obtained after the last time step (m or Pa)
  type(thermo_var_ptrs),                   intent(inout) :: tv       !< Thermodynamical variable structure
  logical,                       optional, intent(in)    :: debug    !< If true, show the call tree

  ! Local variables
  integer :: nk, i, j, k, isd, ied, jsd, jed
  real, dimension(SZI_(G), SZJ_(G), SZK_(G))   :: h_old    ! source grid
  real, dimension(SZI_(G), SZJ_(G), SZK_(G))   :: scalar   ! source data
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)+1) :: dzRegrid ! change in grid interface positions
  logical :: show_call_tree

  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug
  if (show_call_tree) call callTree_enter("regrid_map_T_S(), MOM_ALE.F90")

  ! Build new grid. The new grid is stored in h_new. The old grid is h.
  ! Both are needed for the subsequent remapping of variables.
  call regridding_main( remapCS, regridCS, G, h, tv, dzRegrid )

  call check_remapping_grid( G, h, dzRegrid, 'in regrid_map_T_S()' )

  ! Override old grid with new one. The new grid 'h_new' is built in
  ! one of the 'build_...' routines above.
  nk = G%ke; isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,nk,h,h_old,dzRegrid)
  do k = 1,nk
    do j = jsd,jed ; do i = isd,ied
      h_old(i,j,k) = h(i,j,k)
      h(i,j,k) = h(i,j,k) + ( dzRegrid(i,j,k) - dzRegrid(i,j,k+1) )
    enddo ; enddo
  enddo

  if (show_call_tree) call callTree_waypoint("new grid generated (regrid_map_T_S)")

  ! Remap T and S from old grid h_old onto new grid h
  scalar(:,:,:) = tv%T(:,:,:)
  call remap_scalar_h_to_h( remapCS, G, nk, h_old, scalar, h, tv%T )

  scalar(:,:,:) = tv%S(:,:,:)
  call remap_scalar_h_to_h( remapCS, G, nk, h_old, scalar, h, tv%S )

  if (show_call_tree) call callTree_leave("regrid_map_T_S()")
end subroutine regrid_remap_T_S


!> This routine takes care of remapping all variable between the old and the
!! new grids. When velocity components need to be remapped, thicknesses at
!! velocity points are taken to be arithmetic averages of tracer thicknesses.
!! This routine is called during initialization of the model at time=0, to 
!! remap initiali conditions to the model grid.  It is also called during a
!! time step to update the state.  
subroutine remapping_main(CS_remapping, CS_ALE, G, h, dxInterface, Reg, u, v, debug, dt)
  type(remapping_CS),                               intent(in)    :: CS_remapping  !< Remapping control structure
  type(ALE_CS),                                     intent(in)    :: CS_ALE        !< ALE control structure 
  type(ocean_grid_type),                            intent(in)    :: G             !< Ocean grid structure
  real, dimension(NIMEM_,NJMEM_,NKMEM_),            intent(in)    :: h             !< Level thickness (m or Pa)
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_),     intent(in)    :: dxInterface   !< Change in interface position (Hm or Pa)
  type(tracer_registry_type),                       pointer       :: Reg           !< Tracer registry structure
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), optional, intent(inout) :: u             !< Zonal velocity component (m/s)
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), optional, intent(inout) :: v             !< Meridional velocity component (m/s)
  logical,                                optional, intent(in)    :: debug         !< If true, show the call tree
  real,                                   optional, intent(in)    :: dt            !< time step for diagnostics 

  ! Local variables
  integer                                     :: i, j, k, m
  integer                                     :: nz, ntr
  real, dimension(G%ke+1)                     :: dx
  real, dimension(G%ke)                       :: h1, u_column
  real, dimension(SZI_(G), SZJ_(G), SZK_(G))  :: work_conc 
  real, dimension(SZI_(G), SZJ_(G), SZK_(G))  :: work_cont
  real, dimension(SZI_(G), SZJ_(G))           :: work_2d 
  real                                        :: Idt, h2, ppt2mks   
  logical                                     :: show_call_tree

  show_call_tree = .false.
  if (present(debug)) show_call_tree = debug
  if (show_call_tree) call callTree_enter("remapping_main(), MOM_ALE.F90")

  nz      = G%ke
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
!$OMP parallel default(none) shared(G,h,dxInterface,CS_remapping,nz,Reg,u,v,ntr,show_call_tree) &
!$OMP                       private(h1,dx,u_column)
  if (ntr>0) then
    if (show_call_tree) call callTree_waypoint("remapping tracers (remapping_main)")
!$OMP do
    do m=1,ntr ! For each tracer 

      do j = G%jsc,G%jec
        do i = G%isc,G%iec

          if (G%mask2dT(i,j)>0.) then

            ! Build the start and final grids
            h1(:) = h(i,j,:)
            dx(:) = dxInterface(i,j,:)
            call remapping_core(CS_remapping, nz, h1, Reg%Tr(m)%t(i,j,:), nz, dx, u_column)

            ! Intermediate steps for tendency of tracer concentration and tracer content.
            ! Note: do not merge the two if-tests, since do_tendency_diag(:) is not 
            ! allocated during the time=0 initialization call to this routine. 
            if(present(dt)) then 
              if(CS_ALE%do_tendency_diag(m)) then 
                do k=1,G%ke
                  h2               = h1(k) - (dx(k)-dx(k+1))
                  work_conc(i,j,k) = (u_column(k)    - Reg%Tr(m)%t(i,j,k)      ) * Idt 
                  work_cont(i,j,k) = (u_column(k)*h2 - Reg%Tr(m)%t(i,j,k)*h1(k)) * Idt * G%H_to_kg_m2 
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
              do k=1,G%ke
                do j = G%jsc,G%jec
                  do i = G%isc,G%iec
                    work_cont(i,j,k) = work_cont(i,j,k) * CS_ALE%C_p
                  enddo
                enddo
              enddo
            elseif(trim(Reg%Tr(m)%name) == 'S') then 
              do k=1,G%ke
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
                do k = 1,G%ke
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

  if (show_call_tree) call callTree_waypoint("tracers remapped (remapping_main)")

  ! Remap u velocity component
  if ( present(u) ) then
!$OMP do
    do j = G%jsc,G%jec
      do i = G%iscB,G%iecB
        if (G%mask2dCu(i,j)>0.) then
          ! Build the start and final grids
          h1(:) = 0.5 * ( h(i,j,:) + h(i+1,j,:) )
          dx(:) = 0.5 * ( dxInterface(i,j,:) + dxInterface(i+1,j,:) )
          call remapping_core(CS_remapping, nz, h1, u(i,j,:), nz, dx, u_column)
          u(i,j,:) = u_column(:)
        endif
      enddo
    enddo
  endif

  if (show_call_tree) call callTree_waypoint("u remapped (remapping_main)")

  ! Remap v velocity component
  if ( present(v) ) then
!$OMP do
    do j = G%jscB,G%jecB
      do i = G%isc,G%iec
        if (G%mask2dCv(i,j)>0.) then
          ! Build the start and final grids
          h1(:) = 0.5 * ( h(i,j,:) + h(i,j+1,:) )
          dx(:) = 0.5 * ( dxInterface(i,j,:) + dxInterface(i,j+1,:) )
          call remapping_core(CS_remapping, nz, h1, v(i,j,:), nz, dx, u_column)
          v(i,j,:) = u_column(:)
        endif
      enddo
    enddo
  endif
!$OMP end parallel

  if (show_call_tree) call callTree_waypoint("u remapped (remapping_main)")
  if (show_call_tree) call callTree_leave("remapping_main()")

end subroutine remapping_main


!> Remaps a single scalar between grids described by thicknesses h_src and h_dst.
!! h_dst must be dimensioned as a model array with G%ke layers while h_src can
!! have an arbitrary number of layers specified by nk_src.
subroutine remap_scalar_h_to_h(CS, G, nk_src, h_src, s_src, h_dst, s_dst, all_cells )
  type(remapping_CS),                      intent(in)    :: CS        !< Remapping control structure
  type(ocean_grid_type),                   intent(in)    :: G         !< Ocean grid structure
  integer,                                 intent(in)    :: nk_src    !< Number of levels on source grid
  real, dimension(SZI_(G),SZJ_(G),nk_src), intent(in)    :: h_src     !< Level thickness of source grid (m or Pa)
  real, dimension(SZI_(G),SZJ_(G),nk_src), intent(in)    :: s_src     !< Scalar on source grid
  real, dimension(NIMEM_,NJMEM_,NKMEM_),   intent(in)    :: h_dst     !< Level thickness of destination grid (m or Pa)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),   intent(inout) :: s_dst     !< Scalar on destination grid
  logical, optional,                       intent(in)    :: all_cells !< If false, only reconstruct for
                                                                      !! non-vanished cells. Use all vanished
                                                                      !! layers otherwise (default).
  ! Local variables
  integer :: i, j, k, n_points
  real :: dx(G%ke+1)
  logical :: ignore_vanished_layers

  ignore_vanished_layers = .false.
  if (present(all_cells)) ignore_vanished_layers = .not. all_cells
  n_points = nk_src

!$OMP parallel default(none) shared(CS,G,h_src,s_src,h_dst,s_dst &
!$OMP                               ignore_vanished_layers, nk_src ) &
!$OMP                        private(d,n_pointx)
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
        call dzFromH1H2( n_points, h_src(i,j,1:n_points), G%ke, h_dst(i,j,:), dx )
        call remapping_core(CS, n_points, h_src(i,j,1:n_points), s_src(i,j,1:n_points), G%ke, dx, s_dst(i,j,:))
      else
        s_dst(i,j,:) = 0.
      endif
    enddo
  enddo
!$OMP end parallel

end subroutine remap_scalar_h_to_h


!> Use plm reconstruction for pressure gradient (determine edge values)
!! By using a PLM (limited piecewise linear method) reconstruction, this 
!! routine determines the edge values for the salinity and temperature 
!! within each layer. These edge values are returned and are used to compute 
!! the pressure gradient (by computing the densities).
subroutine pressure_gradient_plm( CS, S_t, S_b, T_t, T_b, G, tv, h )
  type(ALE_CS), intent(inout)                          :: CS   !< module control structure 
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: S_t  !< Salinity at the top edge of each layer
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: S_b  !< Salinity at the bottom edge of each layer
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: T_t  !< Temperature at the top edge of each layer
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: T_b  !< Temperature at the bottom edge of each layer
  type(ocean_grid_type), intent(in)                    :: G    !< ocean grid structure 
  type(thermo_var_ptrs), intent(in)                    :: tv   !< thermodynamics structure 
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h    !< layer thickness 

  ! Local variables
  integer :: i, j, k
  real    :: hTmp(G%ke)
  real    :: tmp(G%ke)
  real, dimension(CS%nk,2)                  :: ppoly_linear_E            !Edge value of polynomial
  real, dimension(CS%nk,CS%degree_linear+1) :: ppoly_linear_coefficients !Coefficients of polynomial

  ! NOTE: the variables 'CS%grid_generic' and 'CS%ppoly_linear' are declared at
  ! the module level. Memory is allocated once at the beginning of the run
  ! in 'ALE_memory_allocation'.

  ! Determine reconstruction within each column
!$OMP parallel do default(none) shared(G,h,tv,CS,S_t,S_b,T_t,T_b)                     &
!$OMP                          private(hTmp,ppoly_linear_E,ppoly_linear_coefficients,tmp)
  do j = G%jsc,G%jec+1
    do i = G%isc,G%iec+1
      ! Build current grid
      hTmp(:) = h(i,j,:)*G%H_to_m
      tmp(:) = tv%S(i,j,:)
      ! Reconstruct salinity profile    
      ppoly_linear_E = 0.0
      ppoly_linear_coefficients = 0.0
      call PLM_reconstruction( G%ke, hTmp, tmp, ppoly_linear_E, ppoly_linear_coefficients )
      if (CS%boundary_extrapolation_for_pressure) call &
        PLM_boundary_extrapolation( G%ke, hTmp, tmp, ppoly_linear_E, ppoly_linear_coefficients )
      
      do k = 1,G%ke
        S_t(i,j,k) = ppoly_linear_E(k,1)
        S_b(i,j,k) = ppoly_linear_E(k,2)
      end do
      
      ! Reconstruct temperature profile 
      ppoly_linear_E = 0.0
      ppoly_linear_coefficients = 0.0
      tmp(:) = tv%T(i,j,:)
      call PLM_reconstruction( G%ke, hTmp, tmp, ppoly_linear_E, ppoly_linear_coefficients )
      if (CS%boundary_extrapolation_for_pressure) call &
        PLM_boundary_extrapolation( G%ke, hTmp, tmp, ppoly_linear_E, ppoly_linear_coefficients )
      
      do k = 1,G%ke
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
subroutine pressure_gradient_ppm( CS, S_t, S_b, T_t, T_b, G, tv, h )
  type(ALE_CS), intent(inout)                          :: CS   !< module control structure 
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: S_t  !< Salinity at top edge of each layer
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: S_b  !< Salinity at bottom edge of each layer
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: T_t  !< Temperature at the top edge of each layer
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: T_b  !< Temperature at the bottom edge of each layer
  type(ocean_grid_type), intent(in)                    :: G    !< ocean grid structure 
  type(thermo_var_ptrs), intent(in)                    :: tv   !< ocean thermodynamics structure 
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h    !< layer thickness 

  ! Local variables
  integer :: i, j, k
  real    :: hTmp(G%ke)
  real    :: tmp(G%ke)
  real, dimension(CS%nk,2) :: &
      ppoly_parab_E            !Edge value of polynomial
  real, dimension(CS%nk,CS%degree_parab+1) :: &
      ppoly_parab_coefficients !Coefficients of polynomial


  ! NOTE: the variables 'CS%grid_generic' and 'CS%ppoly_parab' are declared at
  ! the module level. Memory is allocated once at the beginning of the run
  ! in 'ALE_memory_allocation'.

  ! Determine reconstruction within each column
!$OMP parallel do default(none) shared(G,h,tv,CS,S_t,S_b,T_t,T_b) &
!$OMP                          private(hTmp,tmp,ppoly_parab_E,ppoly_parab_coefficients)
  do j = G%jsc,G%jec+1
    do i = G%isc,G%iec+1
     
      ! Build current grid
      hTmp(:) = h(i,j,:) * G%H_to_m
      tmp(:) = tv%S(i,j,:)
      
      ! Reconstruct salinity profile    
      ppoly_parab_E = 0.0
      ppoly_parab_coefficients = 0.0
      call edge_values_implicit_h4( G%ke, hTmp, tmp, ppoly_parab_E )
      call PPM_reconstruction( G%ke, hTmp, tmp, ppoly_parab_E, ppoly_parab_coefficients )
      if (CS%boundary_extrapolation_for_pressure) call &
        PPM_boundary_extrapolation( G%ke, hTmp, tmp, ppoly_parab_E, ppoly_parab_coefficients )
      
      do k = 1,G%ke
        S_t(i,j,k) = ppoly_parab_E(k,1)
        S_b(i,j,k) = ppoly_parab_E(k,2)
      end do
      
      ! Reconstruct temperature profile 
      ppoly_parab_E = 0.0
      ppoly_parab_coefficients = 0.0
      tmp(:) = tv%T(i,j,:)
      call edge_values_implicit_h4( G%ke, hTmp, tmp, ppoly_parab_E )
      call PPM_reconstruction( G%ke, hTmp, tmp, ppoly_parab_E, ppoly_parab_coefficients )
      if (CS%boundary_extrapolation_for_pressure) call &
        PPM_boundary_extrapolation( G%ke, hTmp, tmp, ppoly_parab_E, ppoly_parab_coefficients )
      
      do k = 1,G%ke
        T_t(i,j,k) = ppoly_parab_E(k,1)
        T_b(i,j,k) = ppoly_parab_E(k,2)
      end do
      
    end do
  end do

end subroutine pressure_gradient_ppm


!> Allocate memory for regridding.
!! In this routine, we allocate the memory needed to carry out regridding
!! steps in the course of the simulation. 
!! For example, to compute implicit edge-value estimates, a tridiagonal system
!! must be solved. We allocate the needed memory at the beginning of the
!! simulation because the number of layers never changes.
subroutine ALE_memory_allocation( G, CS )
  type(ocean_grid_type), intent(in) :: G   !< ocean grid structure 
  type(ALE_CS), intent(inout)       :: CS  !< module control structure 

  ! Local variables
  integer :: nz
  nz = G%ke

  ! Work space
  ALLOC_(CS%dzRegrid(G%isd:G%ied,G%jsd:G%jed,nz+1)); CS%dzRegrid(:,:,:) = 0.
  
end subroutine ALE_memory_allocation


!> Deallocate memory for regridding.
!! In this routine, we reclaim the memory that was allocated for the regridding. 
subroutine ALE_memory_deallocation( CS )
  type(ALE_CS), intent(inout) :: CS  !< module control structure 
  
  ! Work space
  DEALLOC_(CS%dzRegrid)

end subroutine ALE_memory_deallocation

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
subroutine ALE_initRegridding( G, param_file, mod, regridCS, dz )
  type(ocean_grid_type), intent(in)  :: G            !< ocean grid structure 
  type(param_file_type), intent(in)  :: param_file   !< parameter file 
  character(len=*),      intent(in)  :: mod          !< Name of calling module
  type(regridding_CS),   intent(out) :: regridCS     !< Regridding parameters and work arrays
  real, dimension(:),    intent(out) :: dz           !< Resolution (thickness) in units of coordinate

  ! Local variables
  character(len=80)  :: string, varName ! Temporary strings
  character(len=40)  :: coordMode, interpScheme, coordUnits ! Temporary strings
  character(len=200) :: inputdir, fileName
  character(len=320) :: message ! Temporary strings
  integer :: ke
  logical :: tmpLogical
  real :: tmpReal, compress_fraction
  real :: rho_target(G%ke+1) ! Target density used in HYBRID mode

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

  call initialize_regridding( G%ke, coordMode, interpScheme, regridCS, &
                              compressibility_fraction=compress_fraction )

  call get_param(param_file, mod, "ALE_COORDINATE_CONFIG", string, &
                 "Determines how to specify the coordinate\n"//&
                 "resolution. Valid options are:\n"//&
                 " PARAM       - use the vector-parameter ALE_RESOLUTION\n"//&
                 " UNIFORM     - uniformly distributed\n"//&
                 " FILE:string - read from a file. The string specifies\n"//&
                 "               the filename and variable name, separated\n"//&
                 "               by a comma or space, e.g. FILE:lev.nc,Z\n"//&
                 " FNC1:string - FNC1:dz_min,H_total,power,precision",&
                 default='UNIFORM')
!                " FNC1:string - FNC1:dz_min,H_total,power,precision\n"//&
!                " HYBRID:strg  - read from a file. The string specifies\n"//&
!                "               the filename and two variable names, separated\n"//&
!                "               by a comma or space, for sigma-2 and dz. e.g.\n"//&
!                "               HYBRID:vgrid.nc,sigma2,dz",&
  message = "The distribution of vertical resolution for the target\n"//&
            "grid used for Eulerian-like coordinates. For example,\n"//&
            "in z-coordinate mode, the parameter is a list of level\n"//&
            "thicknesses (in m). In sigma-coordinate mode, the list\n"//&
            "is of non-dimensional fractions of the water column."
  select case ( trim(string) )
    case ("UNIFORM")
      dz(:) = uniformResolution(G%ke, coordMode, G%max_depth, &
                 G%Rlay(1)+0.5*(G%Rlay(1)-G%Rlay(2)), &
                 G%Rlay(G%ke)+0.5*(G%Rlay(G%ke)-G%Rlay(G%ke-1)) )
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
          endif
        endif
        if (len_trim(varName)==0) call MOM_error(FATAL,"ALE_initRegridding: "// &
          "Coordinate variable not specified and none could be guessed.")
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
        if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,"ALE_initRegridding: HYBRID "// &
          "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
        call MOM_read_data(trim(fileName), trim(varName), dz)
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
    ! Adjust target grid to be consistent with G%max_depth
    ! This is a work around to the from_Z initialization...  ???
    tmpReal = sum( dz(:) )
    if (tmpReal < G%max_depth) then
      dz(ke) = dz(ke) + ( G%max_depth - tmpReal )
    elseif (tmpReal > G%max_depth) then
      if ( dz(ke) + ( G%max_depth - tmpReal ) > 0. ) then
        dz(ke) = dz(ke) + ( G%max_depth - tmpReal )
      else
        call MOM_error(FATAL,"ALE_initRegridding: "// &
          "MAXIMUM_DEPTH was too shallow to adjust bottom layer of DZ!"//trim(string))
      endif
    endif
  endif
  call setCoordinateResolution( dz, regridCS )
  if (coordinateMode(coordMode) == REGRIDDING_RHO) call set_target_densities_from_G( G, regridCS )

  call get_param(param_file, mod, "MIN_THICKNESS", tmpReal, &
                 "When regridding, this is the minimum layer\n"//&
                 "thickness allowed.", units="m",&
                 default=regriddingDefaultMinThickness )
  call set_regrid_min_thickness( tmpReal, regridCS )

  call get_param(param_file, mod, "BOUNDARY_EXTRAPOLATION", tmpLogical, &
                 "When defined, a proper high-order reconstruction\n"//&
                 "scheme is used within boundary cells rather\n"//&
                 "than PCM. E.g., if PPM is used for remapping, a\n"//&
                 "PPM reconstruction will also be used within\n"//&
                 "boundary cells.", default=regriddingDefaultBoundaryExtrapolation)
  call setRegriddingBoundaryExtrapolation( tmpLogical, regridCS )

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
  dz(nk) = dz(nk) + ( H_total - sum( dz(:) + dz_min ) ) ! Adjust bottom most layer
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
logical function remap_init_conds( CS )
  type(ALE_CS), pointer :: CS   !< module control structure 

  remap_init_conds = .false.
  if (associated(CS)) remap_init_conds = CS%remap_after_initialization
end function remap_init_conds

!> Updates the weights for time filtering the new grid generated in regridding
subroutine ALE_update_regrid_weights( dt, CS )
  real,         intent(in) :: dt !< Time-step used between ALE calls
  type(ALE_CS), pointer    :: CS !< ALE control structure
  ! Local variables
  real :: w

  if (associated(CS)) then
    if (CS%regrid_time_scale <= dt) then
      w = 0.
    else
      w = ( CS%regrid_time_scale - dt ) / CS%regrid_time_scale
    endif
    call set_old_grid_weight( w, CS%regridCS )
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
subroutine ALE_writeCoordinateFile( CS, G, directory )
  type(ALE_CS),          pointer       :: CS         !< module control structure 
  type(ocean_grid_type), intent(inout) :: G          !< model grid structure 
  character(len=*),      intent(in)    :: directory  !< directory for writing grid info 

  character(len=120) :: filepath
  type(vardesc)      :: vars(2)
  type(fieldtype)    :: fields(2)
  integer            :: unit
  real               :: ds(G%ke), dsi(G%ke+1)

  filepath    = trim(directory) // trim("Vertical_coordinate")
  ds(:)       = getCoordinateResolution( CS%regridCS )
  dsi(1)      = 0.5*ds(1)
  dsi(2:G%ke) = 0.5*( ds(1:G%ke-1) + ds(2:G%ke) )
  dsi(G%ke+1) = 0.5*ds(G%ke)

  vars(1) = var_desc('ds', getCoordinateUnits( CS%regridCS ), &
                    'Layer Coordinate Thickness','1','L','1')
  vars(2) = var_desc('ds_interface', getCoordinateUnits( CS%regridCS ), &
                    'Layer Center Coordinate Separation','1','i','1')

  call create_file(unit, trim(filepath), vars, 2, G, fields, SINGLE_FILE)
  call write_field(unit, fields(1), ds)
  call write_field(unit, fields(2), dsi)
  call close_file(unit)

end subroutine ALE_writeCoordinateFile


!> Set h to coordinate values for fixed coordinate systems
subroutine ALE_initThicknessToCoord( CS, G, h )
  type(ALE_CS), intent(inout)                        :: CS  !< module control structure 
  type(ocean_grid_type), intent(in)                  :: G   !< module grid structure 
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(out) :: h   !< layer thickness 

  ! Local variables
  integer :: i, j, k

  do j = G%jsd,G%jed ; do i = G%isd,G%ied
    h(i,j,:) = getStaticThickness( CS%regridCS, 0., G%bathyT(i,j) )
  enddo; enddo

end subroutine ALE_initThicknessToCoord

end module MOM_ALE
