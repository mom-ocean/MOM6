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

use diag_axis_mod,     only : get_diag_axis_name
use diag_manager_mod,  only : diag_axis_init

use netcdf

implicit none ; private

public diag_remap_ctrl, diag_remap_set_diag_axes
public diag_remap_init, diag_remap_end, diag_remap_update, diag_remap_do_remap
public diag_remap_set_vertical_axes, diag_remap_axes_setup_done, diag_remap_get_nz

type :: diag_remap_ctrl
  ! Whether remappping initialized
  logical :: initialized
  ! The vertical coordinate that we remap to
  integer :: vertical_coord
  ! Remap and regrid types needed for remaping using ALE
  type(remapping_CS) :: remap_cs
  type(regridding_CS) :: regrid_cs
  ! Remap grid thicknesses
  real, dimension(:,:,:), allocatable :: h
  ! Number of vertical levels used for remapping
  integer :: nz
  ! Nominal layer thicknesses
  real, dimension(:), allocatable :: dz
  ! Axes id's for the above
  integer :: interface_axes_id
  integer :: layer_axes_id
end type diag_remap_ctrl

contains

subroutine diag_remap_init(remap_cs, vertical_coord)
  type(diag_remap_ctrl), intent(inout) :: remap_cs
  integer, intent(in) :: vertical_coord

  remap_cs%vertical_coord = vertical_coord
  remap_cs%initialized = .false.

end subroutine diag_remap_init

subroutine diag_remap_end(remap_cs)
  type(diag_remap_ctrl), intent(inout) :: remap_cs

  if (allocated(remap_cs%h)) deallocate(remap_cs%h)
  if (allocated(remap_cs%dz)) deallocate(remap_cs%dz)

end subroutine diag_remap_end

subroutine diag_remap_set_vertical_axes(remap_cs, G, GV, param_file)
  type(diag_remap_ctrl), intent(inout) :: remap_cs
  type(ocean_grid_type), intent(inout) :: G !< Ocean grid structure
  type(verticalGrid_type), intent(in)  :: GV !< ocean vertical grid structure
  type(param_file_type), intent(in)    :: param_file !< Parameter file structure

  if (remap_cs%vertical_coord /= coordinateMode('LAYER')) then
    call setup_axes(remap_cs, G, GV, param_file)
  endif

end subroutine diag_remap_set_vertical_axes

subroutine setup_axes(remap_cs, G, GV, param_file)
  type(diag_remap_ctrl),  intent(inout) :: remap_cs !< Diag remap control structure
  type(ocean_grid_type), intent(inout) :: G !< Ocean grid structure
  type(verticalGrid_type), intent(in)  :: GV !< ocean vertical grid structure
  type(param_file_type), intent(in)    :: param_file !< Parameter file structure

  integer :: nzi(4), nzl(4), k
  character(len=200) :: inputdir, string, filename, int_varname, layer_varname
  character(len=40)  :: mod  = "MOM_diag_mediator" ! This module's name.
  character(len=8)   :: units, expected_units
  character(len=34)  :: longname

  real, allocatable, dimension(:) :: interfaces, layers

  ! Read info needed for z-space remapping
  call get_param(param_file, mod, &
                 'DIAG_REMAP_'//trim(vertical_coord_strings(remap_cs%vertical_coord))//'_GRID_DEF', &
                 string, &
                 "This sets the file and variable names that define the\n"//&
                 "vertical grid used for diagnostic output remapping. \n"//&
                 " It should look like:\n"//&
                 " FILE:<file>,<variable> - where <file> is a file within\n"//&
                 "                          the INPUTDIR, <variable> is\n"//&
                 "                          the name of the variable that\n"//&
                 "                          contains interface positions.\n"//&
                 " UNIFORM                - vertical grid is uniform\n"//&
                 "                          between surface and max depth.\n",&
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
      if (.not. check_grid_def(trim(filename), trim(int_varname), &
                                trim(expected_units))) then
        call MOM_error(FATAL,"set_axes_info: Bad grid definition in "//&
                             "'"//trim(filename)//"'")
      endif
      if (.not. check_grid_def(trim(filename), trim(layer_varname), &
                                trim(expected_units))) then
        call MOM_error(FATAL,"set_axes_info: Bad grid definition in "//&
                             "'"//trim(filename)//"'")
      endif

      ! Log the expanded result as a comment since it cannot be read back in
      call log_param(param_file, mod, "! Remapping diagnostics", &
                     trim(inputdir)//"/"//trim(filename)//","//trim(int_varname)//","//trim(layer_varname))

      ! Get interfaces
      call field_size(filename, int_varname, nzi)
      call assert(nzi(1) /= 0, 'set_axes_info: bad interface dimension size')
      allocate(interfaces(nzi(1)))
      call MOM_read_data(filename, int_varname, interfaces)
      ! Always convert heights into depths
      interfaces(:) = abs(interfaces(:))

      ! Get layer dimensions
      call field_size(filename, layer_varname, nzl)
      call assert(nzl(1) /= 0 .and. nzl(1) == nzi(1) - 1, 'set_axes_info: bad layer dimension size')
      allocate(layers(nzl(1)))
      call MOM_read_data(filename, layer_varname, layers)
      ! Always convert heights into depths
      layers(:) = abs(layers(:))

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
    remap_cs%layer_axes_id = diag_axis_init(lowercase(trim(vertical_coord_strings(remap_cs%vertical_coord)))//'_l', &
                                            layers, trim(units), 'z', &
                                            trim(longname)//' at cell center', direction=-1)
    remap_cs%interface_axes_id = diag_axis_init(lowercase(trim(vertical_coord_strings(remap_cs%vertical_coord)))//'_i', &
                                                interfaces, trim(units), 'z', &
                                                trim(longname)//' Depth at interface', direction=-1)
    deallocate(interfaces)
    deallocate(layers)
  else
    ! In this case the axes associated with these will never be used, however
    ! they need to be positive otherwise FMS complains.
    remap_cs%layer_axes_id = 1
    remap_cs%interface_axes_id = 1
  endif

end subroutine setup_axes

function check_grid_def(filename, varname, expected_units)
  ! Do some basic checks on the vertical grid definition file, variable
  character(len=*),   intent(in)  :: filename
  character(len=*),   intent(in)  :: varname
  character(len=*),   intent(in)  :: expected_units
  logical :: check_grid_def, check_units

  character (len=200) :: units, long_name
  integer :: ncid, status, intid, vid

  check_units = .false. ! FIXME: remove this.
  check_grid_def = .true.
  status = NF90_OPEN(filename, NF90_NOWRITE, ncid);
  if (status /= NF90_NOERR) then
    check_grid_def = .false.
  endif

  status = NF90_INQ_VARID(ncid, varname, vid)
  if (status /= NF90_NOERR) then
    check_grid_def = .false.
  endif

  if (check_units) then
    status = NF90_GET_ATT(ncid, vid, "units", units)
    if (status /= NF90_NOERR) then
      check_grid_def = .false.
    endif
    if (trim(units) /= trim(expected_units)) then
      if (trim(expected_units) == "meters") then
        if (trim(units) /= "m") then
          check_grid_def = .false.
        endif
      else
        check_grid_def = .false.
      endif
    endif
  endif

end function check_grid_def

subroutine diag_remap_set_diag_axes(remap_cs, axes)
  type(diag_remap_ctrl), intent(inout) :: remap_cs
  integer, dimension(:), intent(inout) :: axes

  character(len=2) :: axis_name
  integer :: vertical_axis

  call assert(size(axes) == 3, &
              'diag_remap_update_diag_axes: Unexpected number of axes')

  call get_diag_axis_name(axes(3), axis_name)

  if (axis_name == 'zl') then
    axes(3) = remap_cs%layer_axes_id
  elseif (axis_name == 'zi') then
    axes(3) = remap_cs%interface_axes_id
  else
    call assert(.false., &
                'diag_remap_update_diag_axes: Unexpected vertical axis name')
  endif

end subroutine

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
!! height changes.
subroutine diag_remap_update(remap_cs, G, h, T, S, eqn_of_state)
  type(diag_remap_ctrl), intent(inout) :: remap_cs
  type(ocean_grid_type), pointer :: G
  real, dimension(:, :, :), intent(in) :: h, T, S
  type(EOS_type), pointer, intent(in) :: eqn_of_state

  ! Local variables
  integer :: i, j, k, nz
  logical :: checked
  real, dimension(remap_cs%nz + 1) :: zInterfaces
  real, dimension(remap_cs%nz) :: resolution

  if (remap_cs%vertical_coord == coordinateMode('LAYER') .or. &
      .not. diag_remap_axes_setup_done(remap_cs)) then
    return
  endif

  checked = .false.
  nz = remap_cs%nz

  if (.not. remap_cs%initialized) then
    ! Initialize remapping and regridding on the first call
    call initialize_remapping(remap_cs%remap_cs, 'PPM_IH4', boundary_extrapolation=.false.)

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
      endif
      remap_cs%h(i,j,:) = zInterfaces(1:nz) - zInterfaces(2:nz+1)
    enddo
  enddo

  remap_cs%initialized = .true.

end subroutine diag_remap_update

subroutine diag_remap_do_remap(remap_cs, G, h, axes, mask, missing_value, field, remapped_field)
  type(diag_remap_ctrl),  intent(in) :: remap_cs
  type(ocean_grid_type),  intent(in) :: G !< Ocean grid structure
  real, dimension(:,:,:), intent(in) :: h   !<  The current thicknesses
  integer, dimension(3),  intent(in) :: axes !<  Axes for the field
  real, dimension(:,:,:), intent(in) :: mask !<  A mask for the field
  real,                   intent(in) :: missing_value
  real, dimension(:,:,:), intent(in) :: field(:,:,:)
  real, dimension(:,:,:), intent(inout) :: remapped_field !< Field remapped to new coordinate

  ! Local variables
  real, dimension(G%isd:G%ied, G%jsd:G%jed, remap_cs%nz) :: h_dest
  real, dimension(size(h, 3)) :: h_src
  character(len=8) :: x_axis, y_axis, z_axis

  integer :: nz_src, nz_dest
  integer :: i, j, k

  nz_src = size(field, 3)
  nz_dest = remap_cs%nz

  call assert(remap_cs%initialized, 'diag_remap_do_remap: remap_cs not initialized.')
  call assert(size(field, 3) == size(h, 3), &
              'diag_remap_do_remap: Remap field and thickness z-axes do not match.')

  call get_diag_axis_name(axes(1), x_axis)
  call get_diag_axis_name(axes(2), y_axis)
  call get_diag_axis_name(axes(3), z_axis)

  call assert(.not. is_layer_axis(remap_cs, z_axis), &
              'diag_remap_do_remap: Remapping on interfaces not supported')

  remapped_field(:,:,:) = missing_value
  h_dest(:,:,:) = 0.

  if (is_u_axis(x_axis, y_axis)) then
    do j=G%jsc, G%jec
      do I=G%iscB, G%iecB
        if (mask(i,j,1)+mask(i+1,j,1) == 0.) cycle
        h_src(:) = 0.5 * (h(i,j,:) + h(i+1,j,:))
        h_dest(i, j, :) = 0.5 * (remap_cs%h(i,j,:) + remap_cs%h(i+1,j,:))
        call remapping_core_h(nz_src, h_src(:), field(I,j,:), &
                              nz_dest, h_dest(i, j, :), remapped_field(I,j,:), &
                              remap_cs%remap_cs)
      enddo
    enddo
  elseif (is_v_axis(x_axis, y_axis)) then
    do J=G%jscB, G%jecB
      do i=G%isc, G%iec
        if (mask(i,j,1)+mask(i,j+1,1) == 0.) cycle
        h_src(:) = 0.5 * (h(i,j,:) + h(i,j+1,:))
        h_dest(i, j, :) = 0.5 * (remap_cs%h(i,j,:) + remap_cs%h(i,j+1,:) )
        call remapping_core_h(nz_src, h_src(:), field(i,J,:), &
                              nz_dest, h_dest(i, j, :), remapped_field(i,J,:), &
                              remap_cs%remap_cs)
      enddo
    enddo
  elseif (is_B_axis(x_axis, y_axis)) then
    do j=G%jscB, G%jecB
      do I=G%iscB, G%iecB
        if (mask(i,j+1,1)+mask(i+1,j,1) == 0.) cycle
        h_src(:) = 0.5 * (h(i,j+1,:) + h(i+1,j,:))
        h_dest(i,j,:) = 0.5 * (remap_cs%h(i,j+1,:) + remap_cs%h(i+1,j,:))
        call remapping_core_h(nz_src, h_src(:), field(I,j,:), &
                              nz_dest, h_dest(i,j,:), remapped_field(I,j,:), &
                              remap_cs%remap_cs)
      enddo
    enddo
  elseif (is_T_axis(x_axis, y_axis)) then
    do j=G%jsc, G%jec
      do i=G%isc, G%iec
        if (mask(i,j, 1) == 0.) cycle
        h_dest(i,j,:) = remap_cs%h(i,j,:)
        call remapping_core_h(nz_src, h(i,j,:), field(i,j,:), &
                              nz_dest, h_dest(i,j,:), remapped_field(i,j,:), &
                              remap_cs%remap_cs)
      enddo
    enddo
  else
    call assert(.false., 'diag_remap_do_remap: Unsupported axis combination')
  endif

  ! Now mask out zero thickness layers, only applicable to zstar output.
  ! Notice we keep h_dest just for this purpose.
  if (remap_cs%vertical_coord == coordinateMode('ZSTAR')) then
    do j=G%jsd, G%jed
      do i=G%isd, G%ied
        do k=1, nz_dest
          ! Find first zero-thickness layer and mask from there to the bottom.
          if (h_dest(i, j, k) == 0.) then
            remapped_field(i,j,k:nz_dest) = missing_value
            exit
          endif
        enddo
      enddo
    enddo
  endif

end subroutine diag_remap_do_remap

function is_u_axis(x_axis, y_axis)
  character(len=*), intent(in) :: x_axis, y_axis
  logical :: is_u_axis

  if (trim(x_axis) == 'xq' .and. trim(y_axis) == 'yh') then
    is_u_axis = .true.
  else
    is_u_axis = .false.
  endif

end function is_u_axis

function is_v_axis(x_axis, y_axis)
  character(len=*), intent(in) :: x_axis, y_axis
  logical :: is_v_axis

  if (trim(x_axis) == 'xh' .and. trim(y_axis) == 'yq') then
    is_v_axis = .true.
  else
    is_v_axis = .false.
  endif

end function is_v_axis

function is_B_axis(x_axis, y_axis)
  character(len=*), intent(in) :: x_axis, y_axis
  logical :: is_B_axis

  if (trim(x_axis) == 'xq' .and. trim(y_axis) == 'yq') then
    is_B_axis = .true.
  else
    is_B_axis = .false.
  endif

end function is_B_axis

function is_T_axis(x_axis, y_axis)
  character(len=*), intent(in) :: x_axis, y_axis
  logical :: is_T_axis

  if (trim(x_axis) == 'xh' .and. trim(y_axis) == 'yh') then
    is_T_axis = .true.
  else
    is_T_axis = .false.
  endif

end function is_T_axis

function is_layer_axis(remap_cs, z_axis)
  type(diag_remap_ctrl), intent(in) :: remap_cs
  character(len=*), intent(in) :: z_axis
  logical :: is_layer_axis

  if (z_axis == trim(vertical_coord_strings(remap_cs%vertical_coord))//'_zl') then
    is_layer_axis = .true.
  else
    is_layer_axis = .false.
  endif

end function is_layer_axis

end module MOM_diag_remap
