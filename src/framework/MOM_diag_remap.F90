!> This module is used for runtime remapping of diagnostics to z star, sigma and
!! rho vertical coordinates. It defines the diag_remap_ctrl type which
!! represents a remapping of diagnostics to a particular vertical coordinate.
!! The module is used by the diag mediator module in the following way:
!! 1) _init() is called to initialise a diag_remap_ctrl instance.
!! 2) _configure_axes() is called to read the configuration file and set up the
!!    vertical coordinate / axes definitions.
!! 3) _get_axes_info() returns information needed for the diag mediator to
!!    define new axes for the remapped diagnostics.
!! 4) _update() is called periodically (whenever h, T or S change) to either
!!    create or update the target remapping grids.
!! 5) _do_remap() is called from within a diag post() to do the remapping before
!!    the diagnostic is written out.

module MOM_diag_remap

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms,             only : sum_across_PEs
use MOM_error_handler,    only : MOM_error, FATAL, assert
use MOM_diag_vkernels,    only : interpolate_column, reintegrate_column
use MOM_file_parser,      only : get_param, log_param, param_file_type
use MOM_io,               only : slasher, mom_read_data
use MOM_io,               only : file_exists, field_size
use MOM_string_functions, only : lowercase, extractWord
use MOM_grid,             only : ocean_grid_type
use MOM_verticalGrid,     only : verticalGrid_type
use MOM_EOS,              only : EOS_type
use MOM_remapping,        only : remapping_CS, initialize_remapping
use MOM_remapping,        only : remapping_core_h
use MOM_regridding,       only : regridding_CS, initialize_regridding
use MOM_regridding,       only : build_zstar_column, build_rho_column, build_sigma_column
use MOM_regridding,       only : set_regrid_params, get_regrid_size, uniformResolution
use MOM_regridding,       only : getCoordinateInterfaces
use regrid_consts,        only : coordinateMode

use diag_axis_mod,     only : get_diag_axis_name
use diag_manager_mod,  only : diag_axis_init

implicit none ; private

public diag_remap_ctrl
public diag_remap_init, diag_remap_end, diag_remap_update, diag_remap_do_remap
public diag_remap_configure_axes, diag_remap_axes_configured
public diag_remap_get_axes_info, diag_remap_set_active
public diag_remap_diag_registration_closed
public vertically_reintegrate_diag_field
public vertically_interpolate_diag_field
public horizontally_average_diag_field

!> This type represents remapping of diagnostics to a particular vertical
!! coordinate.
!! There is one of these types for each vertical coordinate. The vertical axes
!! of a diagnostic will reference an instance of this type indicating how (or
!! if) the diagnostic should be vertically remapped when being posted.
type :: diag_remap_ctrl
  logical :: configured = .false. !< Whether vertical coordinate has been configured
  logical :: initialized = .false.  !< Whether remappping initialized
  logical :: used = .false.  !< Whether this coordinate actually gets used.
  integer :: vertical_coord = 0 !< The vertical coordinate that we remap to
  character(len=10) :: vertical_coord_name ='' !< The coordinate name as understood by ALE
  character(len=16) :: diag_coord_name = '' !< A name for the purpose of run-time parameters
  character(len=8) :: diag_module_suffix = '' !< The suffix for the module to appear in diag_table
  type(remapping_CS) :: remap_cs !< Remapping control structure use for this axes
  type(regridding_CS) :: regrid_cs !< Regridding control structure that defines the coordiantes for this axes
  integer :: nz = 0 !< Number of vertical levels used for remapping
  real, dimension(:,:,:), allocatable :: h !< Remap grid thicknesses
  real, dimension(:), allocatable :: dz !< Nominal layer thicknesses
  integer :: interface_axes_id = 0 !< Vertical axes id for remapping at interfaces
  integer :: layer_axes_id = 0 !< Vertical axes id for remapping on layers
end type diag_remap_ctrl

contains

!> Initialize a diagnostic remapping type with the given vertical coordinate.
subroutine diag_remap_init(remap_cs, coord_tuple)
  type(diag_remap_ctrl), intent(inout) :: remap_cs !< Diag remapping control structure
  character(len=*),      intent(in)    :: coord_tuple !< A string in form of
                                                      !! MODULE_SUFFIX PARAMETER_SUFFIX COORDINATE_NAME

  remap_cs%diag_module_suffix = trim(extractWord(coord_tuple, 1))
  remap_cs%diag_coord_name = trim(extractWord(coord_tuple, 2))
  remap_cs%vertical_coord_name = trim(extractWord(coord_tuple, 3))
  remap_cs%vertical_coord = coordinateMode(remap_cs%vertical_coord_name)
  remap_cs%configured = .false.
  remap_cs%initialized = .false.
  remap_cs%used = .false.
  remap_cs%nz = 0

end subroutine diag_remap_init

!> De-init a diagnostic remapping type.
!! Free allocated memory.
subroutine diag_remap_end(remap_cs)
  type(diag_remap_ctrl), intent(inout) :: remap_cs !< Diag remapping control structure

  if (allocated(remap_cs%h)) deallocate(remap_cs%h)
  if (allocated(remap_cs%dz)) deallocate(remap_cs%dz)
  remap_cs%configured = .false.
  remap_cs%initialized = .false.
  remap_cs%used = .false.
  remap_cs%nz = 0

end subroutine diag_remap_end

!> Inform that all diagnostics have been registered.
!! If _set_active() has not been called on the remapping control structure
!! will be disabled. This saves time in the case that a vertical coordinate was
!! configured but no diagnostics which use the coordinate appeared in the
!! diag_table.
subroutine diag_remap_diag_registration_closed(remap_cs)
  type(diag_remap_ctrl), intent(inout) :: remap_cs !< Diag remapping control structure

  if (.not. remap_cs%used) then
    call diag_remap_end(remap_cs)
  endif

end subroutine diag_remap_diag_registration_closed

!> Indicate that this remapping type is actually used by the diag manager.
!! If this is never called then the type will be disabled to save time.
!! See further explanation with diag_remap_registration_closed.
subroutine diag_remap_set_active(remap_cs)
  type(diag_remap_ctrl), intent(inout) :: remap_cs !< Diag remapping control structure

  remap_cs%used = .true.

end subroutine diag_remap_set_active

!> Configure the vertical axes for a diagnostic remapping control structure.
!! Reads a configuration parameters to determine coordinate generation.
subroutine diag_remap_configure_axes(remap_cs, GV, param_file)
  type(diag_remap_ctrl),   intent(inout) :: remap_cs !< Diag remap control structure
  type(verticalGrid_type),    intent(in) :: GV !< ocean vertical grid structure
  type(param_file_type),      intent(in) :: param_file !< Parameter file structure
  ! Local variables
  integer :: nzi(4), nzl(4), k
  character(len=200) :: inputdir, string, filename, int_varname, layer_varname
  character(len=40)  :: mod  = "MOM_diag_remap" ! This module's name.
  character(len=8)   :: units, expected_units
  character(len=34)  :: longname, string2

  character(len=256) :: err_msg
  logical :: ierr

  real, allocatable, dimension(:) :: interfaces, layers

  call initialize_regridding(remap_cs%regrid_cs, GV, GV%max_depth, param_file, mod, &
           trim(remap_cs%vertical_coord_name), "DIAG_COORD", trim(remap_cs%diag_coord_name))
  call set_regrid_params(remap_cs%regrid_cs, min_thickness=0., integrate_downward_for_e=.false.)

  remap_cs%nz = get_regrid_size(remap_cs%regrid_cs)

  if (remap_cs%vertical_coord == coordinateMode('SIGMA')) then
    units = 'nondim'
    longname = 'Fraction'
  elseif (remap_cs%vertical_coord == coordinateMode('RHO')) then
    units = 'kg m-3'
    longname = 'Target Potential Density'
  else
    units = 'meters'
    longname = 'Depth'
  endif

  ! Make axes objects
  allocate(interfaces(remap_cs%nz+1))
  allocate(layers(remap_cs%nz))

  interfaces(:) = getCoordinateInterfaces(remap_cs%regrid_cs)
  layers(:) = 0.5 * ( interfaces(1:remap_cs%nz) + interfaces(2:remap_cs%nz+1) )

  remap_cs%interface_axes_id = diag_axis_init(lowercase(trim(remap_cs%diag_coord_name))//'_i', &
                                              interfaces, trim(units), 'z', &
                                              trim(longname)//' at interface', direction=-1)
  remap_cs%layer_axes_id = diag_axis_init(lowercase(trim(remap_cs%diag_coord_name))//'_l', &
                                          layers, trim(units), 'z', &
                                          trim(longname)//' at cell center', direction=-1, &
                                          edges=remap_cs%interface_axes_id)

  ! Axes have now been configured.
  remap_cs%configured = .true.

  deallocate(interfaces)
  deallocate(layers)

end subroutine diag_remap_configure_axes

!> Get layer and interface axes ids for this coordinate
!! Needed when defining axes groups.
subroutine diag_remap_get_axes_info(remap_cs, nz, id_layer, id_interface)
  type(diag_remap_ctrl), intent(in) :: remap_cs !< Diagnostic coordinate control structure
  integer, intent(out) :: nz !< Number of vertical levels for the coordinate
  integer, intent(out) :: id_layer !< 1D-axes id for layer points
  integer, intent(out) :: id_interface !< 1D-axes id for interface points

  nz = remap_cs%nz
  id_layer = remap_cs%layer_axes_id
  id_interface = remap_cs%interface_axes_id

end subroutine diag_remap_get_axes_info


!> Whether or not the axes for this vertical coordinated has been configured.
!! Configuration is complete when diag_remap_configure_axes() has been
!! successfully called.
function diag_remap_axes_configured(remap_cs)
  type(diag_remap_ctrl), intent(in) :: remap_cs
  logical :: diag_remap_axes_configured

  diag_remap_axes_configured = remap_cs%configured

end function

!> Build/update target vertical grids for diagnostic remapping.
!! \note The target grids need to be updated whenever sea surface
!! height or layer thicknesses changes. In the case of density-based
!! coordinates then technically we should also regenerate the
!! target grid whenever T/S change.
subroutine diag_remap_update(remap_cs, G, h, T, S, eqn_of_state)
  type(diag_remap_ctrl), intent(inout) :: remap_cs !< Diagnostic coordinate control structure
  type(ocean_grid_type), pointer :: G !< The ocean's grid type
  real, dimension(:, :, :), intent(in) :: h, T, S !< New thickness, T and S
  type(EOS_type), pointer, intent(in) :: eqn_of_state !< A pointer to the equation of state

  ! Local variables
  integer :: i, j, k, nz
  real, dimension(remap_cs%nz + 1) :: zInterfaces
  real, dimension(remap_cs%nz) :: resolution

  ! Note that coordinateMode('LAYER') is never 'configured' so will
  ! always return here.
  if (.not. remap_cs%configured) then
    return
  endif

  nz = remap_cs%nz

  if (.not. remap_cs%initialized) then
    ! Initialize remapping and regridding on the first call
    call initialize_remapping(remap_cs%remap_cs, 'PPM_IH4', boundary_extrapolation=.false.)
    allocate(remap_cs%h(G%isd:G%ied,G%jsd:G%jed, nz))
    remap_cs%initialized = .true.
  endif

  ! Calculate remapping thicknesses for different target grids based on
  ! nominal/target interface locations. This happens for every call on the
  ! assumption that h, T, S has changed.
  do j=G%jsd, G%jed
    do i=G%isd, G%ied
      if (G%mask2dT(i,j)==0.) then
        remap_cs%h(i,j,:) = 0.
        cycle
      endif

      if (remap_cs%vertical_coord == coordinateMode('ZSTAR')) then
        call build_zstar_column(remap_cs%regrid_cs, nz, &
                              G%bathyT(i,j), sum(h(i,j,:)), zInterfaces)
      elseif (remap_cs%vertical_coord == coordinateMode('SIGMA')) then
        call build_sigma_column(remap_cs%regrid_cs, nz, &
                                G%bathyT(i,j), sum(h(i,j,:)), zInterfaces)
      elseif (remap_cs%vertical_coord == coordinateMode('RHO')) then
        call build_rho_column(remap_cs%regrid_cs, remap_cs%remap_cs, nz, &
                              G%bathyT(i,j), h(i,j,:), T(i, j, :), S(i, j, :), &
                              eqn_of_state, zInterfaces)
      elseif (remap_cs%vertical_coord == coordinateMode('SLIGHT')) then
!       call build_slight_column(remap_cs%regrid_cs,remap_cs%remap_cs, nz, &
!                             G%bathyT(i,j), sum(h(i,j,:)), zInterfaces)
        call MOM_error(FATAL,"diag_remap_update: SLIGHT coordinate not coded for diagnostics yet!")
      elseif (remap_cs%vertical_coord == coordinateMode('HYCOM1')) then
!       call build_hycom1_column(remap_cs%regrid_cs, nz, &
!                             G%bathyT(i,j), sum(h(i,j,:)), zInterfaces)
        call MOM_error(FATAL,"diag_remap_update: HYCOM1 coordinate not coded for diagnostics yet!")
      endif
      remap_cs%h(i,j,:) = zInterfaces(1:nz) - zInterfaces(2:nz+1)
    enddo
  enddo

end subroutine diag_remap_update

!> Remap diagnostic field to alternative vertical grid.
subroutine diag_remap_do_remap(remap_cs, G, h, staggered_in_x, staggered_in_y, &
                               mask, missing_value, field, remapped_field)
  type(diag_remap_ctrl),  intent(in) :: remap_cs !< Diagnostic coodinate control structure
  type(ocean_grid_type),  intent(in) :: G !< Ocean grid structure
  real, dimension(:,:,:), intent(in) :: h   !< The current thicknesses
  logical,                intent(in) :: staggered_in_x !< True is the x-axis location is at u or q points
  logical,                intent(in) :: staggered_in_y !< True is the y-axis location is at v or q points
  real, dimension(:,:,:), pointer    :: mask !< A mask for the field
  real,                   intent(in) :: missing_value !< A missing_value to assign land/vanished points
  real, dimension(:,:,:), intent(in) :: field(:,:,:) !< The diagnostic field to be remapped
  real, dimension(:,:,:), intent(inout) :: remapped_field !< Field remapped to new coordinate
  ! Local variables
  real, dimension(remap_cs%nz) :: h_dest
  real, dimension(size(h,3)) :: h_src
  logical :: mask_vanished_layers
  integer :: nz_src, nz_dest
  integer :: i, j, k

  call assert(remap_cs%initialized, 'diag_remap_do_remap: remap_cs not initialized.')
  call assert(size(field, 3) == size(h, 3), &
              'diag_remap_do_remap: Remap field and thickness z-axes do not match.')

  nz_src = size(field,3)
  nz_dest = remap_cs%nz
  mask_vanished_layers = (remap_cs%vertical_coord == coordinateMode('ZSTAR'))
  remapped_field(:,:,:) = missing_value

  if (staggered_in_x .and. .not. staggered_in_y) then
    ! U-points
    do j=G%jsc, G%jec
      do I=G%iscB, G%iecB
        if (associated(mask)) then
          if (mask(i,j,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i,j,:) + h(i+1,j,:))
        h_dest(:) = 0.5 * (remap_cs%h(i,j,:) + remap_cs%h(i+1,j,:))
        call remapping_core_h(nz_src, h_src(:), field(I,j,:), &
                              nz_dest, h_dest(:), remapped_field(I,j,:), &
                              remap_cs%remap_cs)
        if (mask_vanished_layers) then ! This only works for z-like output
          do k=1, nz_dest
            if (h_dest(k) == 0.) remapped_field(i, j, k:nz_dest) = missing_value
          enddo
        endif
      enddo
    enddo
  elseif (staggered_in_y .and. .not. staggered_in_x) then
    ! V-points
    do J=G%jscB, G%jecB
      do i=G%isc, G%iec
        if (associated(mask)) then
          if (mask(i,j,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i,j,:) + h(i,j+1,:))
        h_dest(:) = 0.5 * (remap_cs%h(i,j,:) + remap_cs%h(i,j+1,:) )
        call remapping_core_h(nz_src, h_src(:), field(i,J,:), &
                              nz_dest, h_dest(:), remapped_field(i,J,:), &
                              remap_cs%remap_cs)
        if (mask_vanished_layers) then ! This only works for z-like output
          do k=1, nz_dest
            if (h_dest(k) == 0.) remapped_field(i,j,k) = missing_value
          enddo
        endif
      enddo
    enddo
  elseif ((.not. staggered_in_x) .and. (.not. staggered_in_y)) then
    ! H-points
    do j=G%jsc, G%jec
      do i=G%isc, G%iec
        if (associated(mask)) then
          if (mask(i,j, 1) == 0.) cycle
        endif
        h_dest(:) = remap_cs%h(i,j,:)
        call remapping_core_h(nz_src, h(i,j,:), field(i,j,:), &
                              nz_dest, h_dest(:), remapped_field(i,j,:), &
                              remap_cs%remap_cs)
        if (mask_vanished_layers) then ! This only works for z-like output
          do k=1, nz_dest
            if (h_dest(k)==0.) remapped_field(i,j,k) = missing_value
          enddo
        endif
      enddo
    enddo
  else
    call assert(.false., 'diag_remap_do_remap: Unsupported axis combination')
  endif

end subroutine diag_remap_do_remap

!> Vertically re-grid an already vertically-integrated diagnostic field to alternative vertical grid.
subroutine vertically_reintegrate_diag_field(remap_cs, G, h, staggered_in_x, staggered_in_y, &
                                             mask, missing_value, field, reintegrated_field)
  type(diag_remap_ctrl),  intent(in) :: remap_cs !< Diagnostic coodinate control structure
  type(ocean_grid_type),  intent(in) :: G !< Ocean grid structure
  real, dimension(:,:,:), intent(in) :: h   !< The current thicknesses
  logical,                intent(in) :: staggered_in_x !< True is the x-axis location is at u or q points
  logical,                intent(in) :: staggered_in_y !< True is the y-axis location is at v or q points
  real, dimension(:,:,:), pointer    :: mask !< A mask for the field
  real,                   intent(in) :: missing_value !< A missing_value to assign land/vanished points
  real, dimension(:,:,:), intent(in) :: field !<  The diagnostic field to be remapped
  real, dimension(:,:,:), intent(inout) :: reintegrated_field !< Field argument remapped to alternative coordinate
  ! Local variables
  real, dimension(remap_cs%nz) :: h_dest
  real, dimension(size(h,3)) :: h_src
  integer :: nz_src, nz_dest
  integer :: i, j, k

  call assert(remap_cs%initialized, 'vertically_reintegrate_diag_field: remap_cs not initialized.')
  call assert(size(field, 3) == size(h, 3), &
              'vertically_reintegrate_diag_field: Remap field and thickness z-axes do not match.')

  nz_src = size(field,3)
  nz_dest = remap_cs%nz
  reintegrated_field(:,:,:) = missing_value

  if (staggered_in_x .and. .not. staggered_in_y) then
    ! U-points
    do j=G%jsc, G%jec
      do I=G%iscB, G%iecB
        if (associated(mask)) then
          if (mask(i,j,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i,j,:) + h(i+1,j,:))
        h_dest(:) = 0.5 * ( remap_cs%h(i,j,:) + remap_cs%h(i+1,j,:) )
        call reintegrate_column(nz_src, h_src, field(I,j,:), &
                                nz_dest, h_dest, missing_value, reintegrated_field(I,j,:))
      enddo
    enddo
  elseif (staggered_in_y .and. .not. staggered_in_x) then
    ! V-points
    do J=G%jscB, G%jecB
      do i=G%isc, G%iec
        if (associated(mask)) then
          if (mask(i,j,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i,j,:) + h(i,j+1,:))
        h_dest(:) = 0.5 * ( remap_cs%h(i,j,:) + remap_cs%h(i,j+1,:) )
        call reintegrate_column(nz_src, h_src, field(i,J,:), &
                                nz_dest, h_dest, missing_value, reintegrated_field(i,J,:))
      enddo
    enddo
  elseif ((.not. staggered_in_x) .and. (.not. staggered_in_y)) then
    ! H-points
    do j=G%jsc, G%jec
      do i=G%isc, G%iec
        if (associated(mask)) then
          if (mask(i,j, 1) == 0.) cycle
        endif
        h_src(:) = h(i,j,:)
        h_dest(:) = remap_cs%h(i,j,:)
        call reintegrate_column(nz_src, h_src, field(i,j,:), &
                                nz_dest, h_dest, missing_value, reintegrated_field(i,j,:))
      enddo
    enddo
  else
    call assert(.false., 'vertically_reintegrate_diag_field: Q point remapping is not coded yet.')
  endif

end subroutine vertically_reintegrate_diag_field

!> Vertically interpolate diagnostic field to alternative vertical grid.
subroutine vertically_interpolate_diag_field(remap_cs, G, h, staggered_in_x, staggered_in_y, &
                                             mask, missing_value, field, interpolated_field)
  type(diag_remap_ctrl),  intent(in) :: remap_cs !< Diagnostic coodinate control structure
  type(ocean_grid_type),  intent(in) :: G !< Ocean grid structure
  real, dimension(:,:,:), intent(in) :: h   !< The current thicknesses
  logical,                intent(in) :: staggered_in_x !< True is the x-axis location is at u or q points
  logical,                intent(in) :: staggered_in_y !< True is the y-axis location is at v or q points
  real, dimension(:,:,:), pointer    :: mask !< A mask for the field
  real,                   intent(in) :: missing_value !< A missing_value to assign land/vanished points
  real, dimension(:,:,:), intent(in) :: field !<  The diagnostic field to be remapped
  real, dimension(:,:,:), intent(inout) :: interpolated_field !< Field argument remapped to alternative coordinate
  ! Local variables
  real, dimension(remap_cs%nz) :: h_dest
  real, dimension(size(h,3)) :: h_src
  integer :: nz_src, nz_dest
  integer :: i, j, k

  call assert(remap_cs%initialized, 'vertically_interpolate_diag_field: remap_cs not initialized.')
  call assert(size(field, 3) == size(h, 3)+1, &
              'vertically_interpolate_diag_field: Remap field and thickness z-axes do not match.')

  interpolated_field(:,:,:) = missing_value

  nz_src = size(h,3)
  nz_dest = remap_cs%nz

  if (staggered_in_x .and. .not. staggered_in_y) then
    ! U-points
    do j=G%jsc, G%jec
      do I=G%iscB, G%iecB
        if (associated(mask)) then
          if (mask(i,j,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i,j,:) + h(i+1,j,:))
        h_dest(:) = 0.5 * ( remap_cs%h(i,j,:) + remap_cs%h(i+1,j,:) )
        call interpolate_column(nz_src, h_src, field(I,j,:), &
                                nz_dest, h_dest, missing_value, interpolated_field(I,j,:))
      enddo
    enddo
  elseif (staggered_in_y .and. .not. staggered_in_x) then
    ! V-points
    do J=G%jscB, G%jecB
      do i=G%isc, G%iec
        if (associated(mask)) then
          if (mask(i,j,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i,j,:) + h(i,j+1,:))
        h_dest(:) = 0.5 * ( remap_cs%h(i,j,:) + remap_cs%h(i,j+1,:) )
        call interpolate_column(nz_src, h_src, field(i,J,:), &
                                nz_dest, h_dest, missing_value, interpolated_field(i,J,:))
      enddo
    enddo
  elseif ((.not. staggered_in_x) .and. (.not. staggered_in_y)) then
    ! H-points
    do j=G%jsc, G%jec
      do i=G%isc, G%iec
        if (associated(mask)) then
          if (mask(i,j, 1) == 0.) cycle
        endif
        h_src(:) = h(i,j,:)
        h_dest(:) = remap_cs%h(i,j,:)
        call interpolate_column(nz_src, h_src, field(i,j,:), &
                                nz_dest, h_dest, missing_value, interpolated_field(i,j,:))
      enddo
    enddo
  else
    call assert(.false., 'vertically_interpolate_diag_field: Q point remapping is not coded yet.')
  endif

end subroutine vertically_interpolate_diag_field

!> Horizontally average field
subroutine horizontally_average_diag_field(G, h, staggered_in_x, staggered_in_y, &
                                           is_layer, is_extensive, &
                                           missing_value, field, averaged_field)
  type(ocean_grid_type),  intent(in) :: G !< Ocean grid structure
  real, dimension(:,:,:), intent(in) :: h !< The current thicknesses
  logical,                intent(in) :: staggered_in_x !< True if the x-axis location is at u or q points
  logical,                intent(in) :: staggered_in_y !< True if the y-axis location is at v or q points
  logical,                intent(in) :: is_layer !< True if the z-axis location is at h points
  logical,                intent(in) :: is_extensive !< True if the z-direction is spatially integrated (over layers)
  real,                   intent(in) :: missing_value !< A missing_value to assign land/vanished points
  real, dimension(:,:,:), intent(in) :: field !<  The diagnostic field to be remapped
  real, dimension(:),  intent(inout) :: averaged_field !< Field argument horizontally averaged
  ! Local variables
  real, dimension(size(field, 3)) :: vol_sum, stuff_sum ! nz+1 is needed for interface averages
  real :: v1, v2, total_volume, total_stuff, val
  integer :: i, j, k, nz

  nz = size(field, 3)

  if (staggered_in_x .and. .not. staggered_in_y) then
    if (is_layer) then
      ! U-points
      do k=1,nz
        vol_sum(k) = 0.
        stuff_sum(k) = 0.
        if (is_extensive) then
          do j=G%jsc, G%jec ; do i=G%isc, G%iec
            v1 = G%areaCu(I,j)
            v2 = G%areaCu(I-1,j)
            vol_sum(k) = vol_sum(k) + 0.5 * ( v1 + v2 ) * G%mask2dT(i,j)
            stuff_sum(k) = stuff_sum(k) + 0.5 * ( v1 * field(I,j,k) + v2 * field(I-1,j,k) ) * G%mask2dT(i,j)
          enddo ; enddo
        else ! Intensive
          do j=G%jsc, G%jec ; do i=G%isc, G%iec
            v1 = G%areaCu(I,j) * 0.5 * ( h(i,j,k) + h(i+1,j,k) )
            v2 = G%areaCu(I-1,j) * 0.5 * ( h(i,j,k) + h(i-1,j,k) )
            vol_sum(k) = vol_sum(k) + 0.5 * ( v1 + v2 ) * G%mask2dT(i,j)
            stuff_sum(k) = stuff_sum(k) + 0.5 * ( v1 * field(I,j,k) + v2 * field(I-1,j,k) ) * G%mask2dT(i,j)
          enddo ; enddo
        endif
      enddo
    else ! Interface
      do k=1,nz
        vol_sum(k) = 0.
        stuff_sum(k) = 0.
        do j=G%jsc, G%jec ; do i=G%isc, G%iec
          v1 = G%areaCu(I,j)
          v2 = G%areaCu(I-1,j)
          vol_sum(k) = vol_sum(k) + 0.5 * ( v1 + v2 ) * G%mask2dT(i,j)
          stuff_sum(k) = stuff_sum(k) + 0.5 * ( v1 * field(I,j,k) + v2 * field(I-1,j,k) ) * G%mask2dT(i,j)
        enddo ; enddo
      enddo
    endif
  elseif (staggered_in_y .and. .not. staggered_in_x) then
    if (is_layer) then
      ! V-points
      do k=1,nz
        vol_sum(k) = 0.
        stuff_sum(k) = 0.
        if (is_extensive) then
          do j=G%jsc, G%jec ; do i=G%isc, G%iec
            v1 = G%areaCv(i,J)
            v2 = G%areaCv(i,J-1)
            vol_sum(k) = vol_sum(k) + 0.5 * ( v1 + v2 ) * G%mask2dT(i,j)
            stuff_sum(k) = stuff_sum(k) + 0.5 * ( v1 * field(i,J,k) + v2 * field(i,J-1,k) ) * G%mask2dT(i,j)
          enddo ; enddo
        else ! Intensive
          do j=G%jsc, G%jec ; do i=G%isc, G%iec
            v1 = G%areaCv(i,J) * 0.5 * ( h(i,j,k) + h(i,j+1,k) )
            v2 = G%areaCv(i,J-1) * 0.5 * ( h(i,j,k) + h(i,j-1,k) )
            vol_sum(k) = vol_sum(k) + 0.5 * ( v1 + v2 ) * G%mask2dT(i,j)
            stuff_sum(k) = stuff_sum(k) + 0.5 * ( v1 * field(i,J,k) + v2 * field(i,J-1,k) ) * G%mask2dT(i,j)
          enddo ; enddo
        endif
      enddo
    else ! Interface
      do k=1,nz
        vol_sum(k) = 0.
        stuff_sum(k) = 0.
        do j=G%jsc, G%jec ; do i=G%isc, G%iec
          v1 = G%areaCv(i,J)
          v2 = G%areaCv(i,J-1)
          vol_sum(k) = vol_sum(k) + 0.5 * ( v1 + v2 ) * G%mask2dT(i,j)
          stuff_sum(k) = stuff_sum(k) + 0.5 * ( v1 * field(i,J,k) + v2 * field(i,J-1,k) ) * G%mask2dT(i,j)
        enddo ; enddo
      enddo
    endif
  elseif ((.not. staggered_in_x) .and. (.not. staggered_in_y)) then
    if (is_layer) then
      ! H-points
      do k=1,nz
        vol_sum(k) = 0.
        stuff_sum(k) = 0.
        if (is_extensive) then
          do j=G%jsc, G%jec ; do i=G%isc, G%iec
            if (G%mask2dT(i,j)>0. .and. h(i,j,k)>0.) then
              v1 = G%areaT(i,j)
              vol_sum(k) = vol_sum(k) + v1
              stuff_sum(k) = stuff_sum(k) + v1 * field(i,j,k)
            endif
          enddo ; enddo
        else ! Intensive
          do j=G%jsc, G%jec ; do i=G%isc, G%iec
            if (G%mask2dT(i,j)>0. .and. h(i,j,k)>0.) then
              v1 = G%areaT(i,j) * h(i,j,k)
              vol_sum(k) = vol_sum(k) + v1
              stuff_sum(k) = stuff_sum(k) + v1 * field(i,j,k)
            endif
          enddo ; enddo
        endif
      enddo
    else ! Interface
      do k=1,nz
        vol_sum(k) = 0.
        stuff_sum(k) = 0.
        do j=G%jsc, G%jec ; do i=G%isc, G%iec
          val = field(i,j,k)
          if (G%mask2dT(i,j)>0. .and. val/=missing_value) then
            v1 = G%areaT(i,j)
            vol_sum(k) = vol_sum(k) + v1
            stuff_sum(k) = stuff_sum(k) + v1 * field(i,j,k)
          endif
        enddo ; enddo
      enddo
    endif
  else
    call assert(.false., 'horizontally_average_diag_field: Q point averaging is not coded yet.')
  endif

  call sum_across_PEs(vol_sum, nz)
  call sum_across_PEs(stuff_sum, nz)

  do k=1,nz
    if (vol_sum(k)>0.) then
      averaged_field(k) = stuff_sum(k) / vol_sum(k)
    else
      averaged_field(k) = missing_value
    endif
  enddo

end subroutine horizontally_average_diag_field

end module MOM_diag_remap
