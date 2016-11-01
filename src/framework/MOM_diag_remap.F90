module MOM_diag_remap

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
!*    The subroutines here are used for runtime remapping of           *
!*    diagnostics.                                                     *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

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
use MOM_regridding,       only : regridding_CS, initialize_regridding, setCoordinateResolution
use MOM_regridding,       only : build_zstar_column, build_rho_column, build_sigma_column
use MOM_regridding,       only : set_regrid_params, uniformResolution
use regrid_consts,        only : coordinateMode, DEFAULT_COORDINATE_MODE, vertical_coord_strings
use regrid_consts,        only : REGRIDDING_ZSTAR

use diag_axis_mod,     only : get_diag_axis_name
use diag_manager_mod,  only : diag_axis_init

use netcdf

implicit none ; private

public diag_remap_ctrl
public diag_remap_init, diag_remap_end, diag_remap_update, diag_remap_do_remap
public diag_remap_set_vertical_axes, diag_remap_axes_setup_done, diag_remap_get_nz
public diag_remap_get_vertical_ids
public vertically_reintegrate_diag_field
public vertically_interpolate_diag_field

!> This type represents a remapping of a diagnostic to a particular vertical
!! coordinate.
!! There is one of these types for each vertical coordinate. The vertical axes
!! of a diagnostic will reference an instance of this type indicating how (or
!! if) the diagnostic should be vertically remapped when being posted.
type :: diag_remap_ctrl
  logical :: defined = .false. !< Whether a coordinate has been defined
  logical :: initialized = .false.  !< Whether remappping initialized
  integer :: vertical_coord = 0 !< The vertical coordinate that we remap to
  type(remapping_CS), pointer :: remap_cs => null() !< type for remapping using ALE module
  type(regridding_CS), pointer :: regrid_cs => null() !< type for regridding using ALE module
  integer :: nz = 0 !< Number of vertical levels used for remapping
  real, dimension(:,:,:), allocatable :: h !< Remap grid thicknesses
  real, dimension(:), allocatable :: dz !< Nominal layer thicknesses
  integer :: interface_axes_id = 0 !< Vertical axes id for remapping at interfaces
  integer :: layer_axes_id = 0 !< Vertical axes id for remapping on layers
end type diag_remap_ctrl

contains

!> Initialize a diagnostic remapping type with the given vertical coordinate.
!! The possible coordinates are defined in module regrid_consts
subroutine diag_remap_init(remap_cs, vertical_coord)
  type(diag_remap_ctrl), intent(inout) :: remap_cs !< Diag remapping control structure
  integer, intent(in) :: vertical_coord !< The vertical coordinate it represents

  remap_cs%vertical_coord = vertical_coord
  remap_cs%initialized = .false.
  remap_cs%defined = .false.
  remap_cs%nz = 0

end subroutine diag_remap_init

!> De-init a diagnostic remapping type.
!! Free allocated memory.
subroutine diag_remap_end(remap_cs)
  type(diag_remap_ctrl), intent(inout) :: remap_cs !< Diag remapping control structure

  if (allocated(remap_cs%h)) deallocate(remap_cs%h)
  if (allocated(remap_cs%dz)) deallocate(remap_cs%dz)
  remap_cs%nz = 0

end subroutine diag_remap_end

!> Configure the vertical axes for a diagnostic remapping control structure.
!! Reads a configuration file to determine nominal location of vertical
!! layers/interfaces.
subroutine diag_remap_set_vertical_axes(remap_cs, G, GV, param_file)
  type(diag_remap_ctrl), intent(inout) :: remap_cs !< Diag remap control structure
  type(ocean_grid_type), intent(inout) :: G !< Ocean grid structure
  type(verticalGrid_type), intent(in)  :: GV !< ocean vertical grid structure
  type(param_file_type), intent(in)    :: param_file !< Parameter file structure

  if (remap_cs%vertical_coord /= coordinateMode('LAYER')) then
    call setup_axes(remap_cs, G, GV, param_file)
  endif

end subroutine diag_remap_set_vertical_axes

!> Read grid definition spec to configure axes.
subroutine setup_axes(remap_cs, G, GV, param_file)
  type(diag_remap_ctrl),  intent(inout) :: remap_cs !< Diag remap control structure
  type(ocean_grid_type), intent(inout) :: G !< Ocean grid structure
  type(verticalGrid_type), intent(in)  :: GV !< ocean vertical grid structure
  type(param_file_type), intent(in)    :: param_file !< Parameter file structure

  integer :: nzi(4), nzl(4), k
  character(len=200) :: inputdir, string, filename, int_varname, layer_varname
  character(len=40)  :: mod  = "MOM_diag_mediator" ! This module's name.
  character(len=8)   :: units, expected_units
  character(len=34)  :: longname, coord_string

  character(len=256) :: err_msg
  logical :: ierr

  real, allocatable, dimension(:) :: interfaces, layers

  ! Read info needed for z-space remapping
  coord_string = trim(vertical_coord_strings(remap_cs%vertical_coord))
  if (coordinateMode(coord_string) == REGRIDDING_ZSTAR) then
    coord_string = 'Z' ! Special case to carry forward original use of _z and _Z for z*-coordinates
  endif
  call get_param(param_file, mod, &
                 'DIAG_REMAP_'//trim(coord_string)//'_GRID_DEF', &
                 string, &
                 "This sets the file and variable names that define the\n"//&
                 "vertical grid used for diagnostic output remapping. \n"//&
                 " It should look like:\n"//&
                 " FILE:<file>,<varI>,<varL> - where <file> is a file within\n"//&
                 "                             the INPUTDIR, <varI> is\n"//&
                 "                             the name of the variable that\n"//&
                 "                             contains interface positions,\n"//&
                 "                             <varL> is the name of the variable\n"//&
                 "                             that contains layer positions,\n"//&
                 "UNIFORM                    - vertical grid is uniform\n"//&
                 "                             between surface and max depth.\n",&
                 default="")
  if (len_trim(string) > 0) then

    if (trim(string) == 'UNIFORM') then
      nzi(1) = GV%ke + 1
      nzl(1) = GV%ke
      allocate(remap_cs%dz(nzl(1)))
      allocate(interfaces(nzi(1)))
      allocate(layers(nzl(1)))

      remap_cs%dz(:) = uniformResolution(GV%ke, vertical_coord_strings(remap_cs%vertical_coord), &
                                      G%max_depth, GV%Rlay(1)+0.5*(GV%Rlay(1)-GV%Rlay(2)), &
                                      GV%Rlay(GV%ke)+0.5*(GV%Rlay(GV%ke)-GV%Rlay(GV%ke-1)))

      if (remap_cs%vertical_coord == coordinateMode('RHO')) then
        interfaces(1) = GV%Rlay(1)
      else
        interfaces(1) = 0
      endif

      ! Calculate interface and layer positions
      do k=2,nzi(1)
        interfaces(k) = interfaces(k - 1) + remap_cs%dz(k - 1)
      enddo
      layers(:) = interfaces(1:nzl(1)) + remap_cs%dz(:) / 2

    elseif (index(trim(string), 'FILE:') == 1) then
      ! read coordinate information from a file
      if (string(6:6)=='.' .or. string(6:6)=='/') then
        inputdir = "."
        filename = trim(extractWord(trim(string(6:200)), 1))
      else
        call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
        inputdir = slasher(inputdir)
        filename = trim(inputdir) // trim(extractWord(trim(string(6:200)), 1))
      endif
      int_varname = trim(extractWord(trim(string(1:200)), 2))
      layer_varname = trim(extractWord(trim(string(1:200)), 3))

      if (.not. file_exists(trim(filename))) then
        call MOM_error(FATAL,"set_axes_info: Specified file not found: "//&
                             "Looking for '"//trim(filename)//"'")
      endif

      if (remap_cs%vertical_coord == coordinateMode('SIGMA')) then
        expected_units = 'nondim'
      elseif (remap_cs%vertical_coord == coordinateMode('RHO')) then
        expected_units = 'kg m-3'
      else
        expected_units = 'meters'
      endif

      ! Check that the vars have expected format, units etc.
      call check_grid_def(filename, int_varname, expected_units, err_msg, ierr)
      if (ierr) then
           call MOM_error(FATAL,"set_axes_info: Unsupported format in grid "//&
                                 "definition '"//trim(filename)//"'"//&
                                 ". Error message "//err_msg)
      endif

      call check_grid_def(filename, layer_varname, expected_units, err_msg, ierr)
      if (ierr) then
           call MOM_error(FATAL,"set_axes_info: Unsupported format in grid "//&
                                 "definition '"//trim(filename)//"'"//&
                                 ". Error message "//err_msg)
      endif

      ! Log the expanded result as a comment since it cannot be read back in
      call log_param(param_file, mod, "! Remapping diagnostics", &
                     trim(filename)//","//trim(int_varname)//","//trim(layer_varname))

      ! Get interfaces
      call field_size(trim(filename), trim(int_varname), nzi)
      call assert(nzi(1) /= 0, 'set_axes_info: bad interface dimension size')
      allocate(interfaces(nzi(1)))
      call MOM_read_data(filename, int_varname, interfaces)
      ! Always convert heights into depths
      interfaces(:) = abs(interfaces(:))

      ! Get layer dimensions
      allocate(layers(nzi(1)-1))
      call field_size(trim(filename), trim(layer_varname), nzl)
      if (trim(layer_varname) /= trim(int_varname)) then
        call assert(nzl(1) /= 0 .and. nzl(1) == nzi(1) - 1, 'set_axes_info: bad layer dimension size')
        call MOM_read_data(filename, trim(layer_varname), layers)
        ! Always convert heights into depths
        layers(:) = abs(layers(:))
      else
        nzl(1) = nzi(1)-1
        layers(:) = 0.5*( interfaces(1:nzi(1)-1) + interfaces(2:nzi(1)) )
      endif

      allocate(remap_cs%dz(nzl(1)))
      remap_cs%dz(:) = interfaces(2:) - interfaces(1:nzi(1)-1)

    else
      ! unsupported method
      call MOM_error(FATAL,"set_axes_info: "//&
         "Incorrect remapping grid specification. Only 'FILE:file,var' and"//&
         "'UNIFORM' are currently supported."//&
         "Found '"//trim(string)//"'")
    endif

    remap_cs%nz = nzl(1)

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
    remap_cs%defined = .true.
    remap_cs%layer_axes_id = diag_axis_init(lowercase(trim(vertical_coord_strings(remap_cs%vertical_coord)))//'_l', &
                                            layers, trim(units), 'z', &
                                            trim(longname)//' at cell center', direction=-1)
    remap_cs%interface_axes_id = diag_axis_init(lowercase(trim(vertical_coord_strings(remap_cs%vertical_coord)))//'_i', &
                                                interfaces, trim(units), 'z', &
                                                trim(longname)//' Depth at interface', direction=-1)
    deallocate(interfaces)
    deallocate(layers)
  else
    ! This coordinate is not in use
    remap_cs%layer_axes_id = -1
    remap_cs%interface_axes_id = -1
  endif

end subroutine setup_axes

subroutine check_grid_def(filename, varname, expected_units, msg, ierr)
  ! Do some basic checks on the vertical grid definition file, variable
  character(len=*),   intent(in)  :: filename
  character(len=*),   intent(in)  :: varname
  character(len=*),   intent(in)  :: expected_units

  logical, intent(out) :: ierr
  character(len=*), intent(inout) :: msg

  character (len=200) :: units, long_name
  integer :: ncid, status, intid, vid

  ierr = .false.
  status = NF90_OPEN(trim(filename), NF90_NOWRITE, ncid);
  if (status /= NF90_NOERR) then
    ierr = .true.
    msg = 'File not found: '//trim(filename)
    return
  endif

  status = NF90_INQ_VARID(ncid, trim(varname), vid)
  if (status /= NF90_NOERR) then
    ierr = .true.
    msg = 'Var not found: '//trim(varname)
    return
  endif

  status = NF90_GET_ATT(ncid, vid, "units", units)
  if (status /= NF90_NOERR) then
    ierr = .true.
    msg = 'Attribute not found: units'
    return
  endif

  if (trim(units) /= trim(expected_units)) then
    if (trim(expected_units) == "meters") then
      if (trim(units) /= "m") then
        ierr = .true.
      endif
    else
      ierr = .true.
    endif
  endif

  if (ierr) then
    msg = 'Units incorrect: '//trim(units)//' /= '//trim(expected_units)
  endif

end subroutine check_grid_def

!> Get layer and interface axes ids for this coordinate
subroutine diag_remap_get_vertical_ids(remap_cs, id_layer, id_interface)
  type(diag_remap_ctrl), intent(in) :: remap_cs !< Diagnostic coordinate control structure
  integer, intent(out) :: id_layer !< 1D-axes id for layer points
  integer, intent(out) :: id_interface !< 1D-axes id for interface points

  id_layer = remap_cs%layer_axes_id
  id_interface = remap_cs%interface_axes_id

end subroutine diag_remap_get_vertical_ids

!> Get the number of vertical levels for this coordinate.
!! This is needed when defining axes groups.
function diag_remap_get_nz(remap_cs)
  type(diag_remap_ctrl), intent(in) :: remap_cs
  integer :: diag_remap_get_nz

  diag_remap_get_nz = remap_cs%nz

end function

function diag_remap_axes_setup_done(remap_cs)
  type(diag_remap_ctrl), intent(in) :: remap_cs
  logical :: diag_remap_axes_setup_done

  if (allocated(remap_cs%dz)) then
    diag_remap_axes_setup_done = .true.
  else
    diag_remap_axes_setup_done = .false.
  endif

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
  logical :: checked
  real, dimension(remap_cs%nz + 1) :: zInterfaces
  real, dimension(remap_cs%nz) :: resolution

  if (remap_cs%vertical_coord == coordinateMode('LAYER') .or. &
      .not. diag_remap_axes_setup_done(remap_cs)) then
    return
  endif

  call assert(remap_cs%defined, 'diag_remap_update: Attempting to update an undefined coordinate!')

  checked = .false.
  nz = remap_cs%nz

  if (.not. remap_cs%initialized) then
    ! Initialize remapping and regridding on the first call
    allocate(remap_cs%remap_cs)
    call initialize_remapping(remap_cs%remap_cs, 'PPM_IH4', boundary_extrapolation=.false.)

    allocate(remap_cs%regrid_cs)
    call initialize_regridding(remap_cs%nz, &
          vertical_coord_strings(remap_cs%vertical_coord), 'PPM_IH4', remap_cs%regrid_cs)
    call set_regrid_params(remap_cs%regrid_cs, min_thickness=0., integrate_downward_for_e=.false.)
    call setCoordinateResolution(remap_cs%dz, remap_cs%regrid_cs)

    allocate(remap_cs%h(G%isd:G%ied,G%jsd:G%jed,remap_cs%nz))
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

  remap_cs%initialized = .true.

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


  call assert(remap_cs%defined, 'diag_remap_do_remap: remap_cs is for an undefined coordinate!')
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
          if (mask(i,j,1)+mask(i+1,j,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i,j,:) + h(i+1,j,:))
        h_dest(:) = 0.5 * (remap_cs%h(i,j,:) + remap_cs%h(i+1,j,:))
        call remapping_core_h(nz_src, h_src(:), field(I,j,:), &
                              nz_dest, h_dest(:), remapped_field(I,j,:), &
                              remap_cs%remap_cs)
        if (mask_vanished_layers) then
          do k=1, nz_dest
            if (h_dest(k) == 0.) then ! This only works for z-like output
              remapped_field(i,j,k:nz_dest) = missing_value
              exit
            endif
          enddo
        endif
      enddo
    enddo
  elseif (staggered_in_y .and. .not. staggered_in_x) then
    ! V-points
    do J=G%jscB, G%jecB
      do i=G%isc, G%iec
        if (associated(mask)) then
          if (mask(i,j,1)+mask(i,j+1,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i,j,:) + h(i,j+1,:))
        h_dest(:) = 0.5 * (remap_cs%h(i,j,:) + remap_cs%h(i,j+1,:) )
        call remapping_core_h(nz_src, h_src(:), field(i,J,:), &
                              nz_dest, h_dest(:), remapped_field(i,J,:), &
                              remap_cs%remap_cs)
        if (mask_vanished_layers) then
          do k=1, nz_dest
            if (h_dest(k) == 0.) then ! This only works for z-like output
              remapped_field(i,j,k:nz_dest) = missing_value
              exit
            endif
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
        if (mask_vanished_layers) then
          do k=1, nz_dest
            if (h_dest(k) == 0.) then ! This only works for z-like output
              remapped_field(i,j,k:nz_dest) = missing_value
              exit
            endif
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
          if (mask(i,j,1)+mask(i+1,j,1) == 0.) cycle
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
          if (mask(i,j,1)+mask(i,j+1,1) == 0.) cycle
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
          if (mask(i,j,1)+mask(i+1,j,1) == 0.) cycle
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
          if (mask(i,j,1)+mask(i,j+1,1) == 0.) cycle
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

end module MOM_diag_remap
