!> provides runtime remapping of diagnostics to z star, sigma and
!! rho vertical coordinates.
!!
!! The diag_remap_ctrl type represents a remapping of diagnostics to a particular
!! vertical coordinate. The module is used by the diag mediator module in the
!! following way:
!! 1. diag_remap_init() is called to initialize a diag_remap_ctrl instance.
!! 2. diag_remap_configure_axes() is called to read the configuration file and set up the
!!    vertical coordinate / axes definitions.
!! 3. diag_remap_get_axes_info() returns information needed for the diag mediator to
!!    define new axes for the remapped diagnostics.
!! 4. diag_remap_update() is called periodically (whenever h, T or S change) to either
!!    create or update the target remapping grids.
!! 5. diag_remap_do_remap() is called from within a diag post() to do the remapping before
!!    the diagnostic is written out.


! NOTE: In the following functions, the fields are passed using 1-based
! indexing, which requires special handling within the grid index loops.
!
!   * diag_remap_do_remap
!   * vertically_reintegrate_diag_field
!   * vertically_interpolate_diag_field
!   * horizontally_average_diag_field
!
! Symmetric grids add an additional row of western and southern points to u-
! and v-grids.  Non-symmetric grids are 1-based and symmetric grids are
! zero-based, allowing the same expressions to be used when accessing the
! fields.  But if u- or v-points become 1-indexed, as in these functions, then
! the stencils must be re-assessed.
!
! For interpolation between h and u grids, we use the following relations:
!
!   h->u: f_u(ig) = 0.5 * (f_h( ig ) + f_h(ig+1))
!         f_u(i1) = 0.5 * (f_h(i1-1) + f_h( i1 ))
!
!   u->h: f_h(ig) = 0.5 * (f_u(ig-1) + f_u( ig ))
!         f_h(i1) = 0.5 * (f_u( i1 ) + f_u(i1+1))
!
! where ig is the grid index and i1 is the 1-based index.  That is, a 1-based
! u-point is ahead of its matching h-point in non-symmetric mode, but behind
! its matching h-point in non-symmetric mode.
!
! We can combine these expressions by applying to ig a -1 shift on u-grids and
! a +1 shift on h-grids in symmetric mode.
!
! We do not adjust the h-point indices, since they are assumed to be 1-based.
! This is only correct when global indexing is disabled.  If global indexing is
! enabled, then all indices will need to be defined relative to the data
! domain.
!
! Finally, note that the mask input fields are pointers to arrays which are
! zero-indexed, and do not need any corrections over grid index loops.


module MOM_diag_remap

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms,             only : reproducing_sum_EFP, EFP_to_real
use MOM_coms,             only : EFP_type, assignment(=), EFP_sum_across_PEs
use MOM_error_handler,    only : MOM_error, FATAL, assert, WARNING
use MOM_debugging,        only : check_column_integrals
use MOM_diag_manager_infra,only : MOM_diag_axis_init
use MOM_diag_vkernels,    only : interpolate_column, reintegrate_column
use MOM_file_parser,      only : get_param, log_param, param_file_type
use MOM_string_functions, only : lowercase, extractWord
use MOM_grid,             only : ocean_grid_type
use MOM_unit_scaling,     only : unit_scale_type
use MOM_verticalGrid,     only : verticalGrid_type
use MOM_EOS,              only : EOS_type
use MOM_remapping,        only : remapping_CS, initialize_remapping
use MOM_remapping,        only : remapping_core_h
use MOM_regridding,       only : regridding_CS, initialize_regridding
use MOM_regridding,       only : set_regrid_params, get_regrid_size
use MOM_regridding,       only : getCoordinateInterfaces
use MOM_regridding,       only : get_zlike_CS, get_sigma_CS, get_rho_CS
use regrid_consts,        only : coordinateMode
use coord_zlike,          only : build_zstar_column
use coord_sigma,          only : build_sigma_column
use coord_rho,            only : build_rho_column


implicit none ; private

public diag_remap_ctrl
public diag_remap_init, diag_remap_end, diag_remap_update, diag_remap_do_remap
public diag_remap_configure_axes, diag_remap_axes_configured
public diag_remap_calc_hmask
public diag_remap_get_axes_info, diag_remap_set_active
public diag_remap_diag_registration_closed
public vertically_reintegrate_diag_field
public vertically_interpolate_diag_field
public horizontally_average_diag_field

!> Represents remapping of diagnostics to a particular vertical coordinate.
!!
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
  type(regridding_CS) :: regrid_cs !< Regridding control structure that defines the coordinates for this axes
  integer :: nz = 0 !< Number of vertical levels used for remapping
  real, dimension(:,:,:), allocatable :: h !< Remap grid thicknesses [H ~> m or kg m-2]
  real, dimension(:,:,:), allocatable :: h_extensive !< Remap grid thicknesses for extensive
                                           !! variables [H ~> m or kg m-2]
  integer :: interface_axes_id = 0 !< Vertical axes id for remapping at interfaces
  integer :: layer_axes_id = 0 !< Vertical axes id for remapping on layers
  logical :: answers_2018      !< If true, use the order of arithmetic and expressions for remapping
                               !! that recover the answers from the end of 2018. Otherwise, use
                               !! updated more robust forms of the same expressions.
end type diag_remap_ctrl

contains

!> Initialize a diagnostic remapping type with the given vertical coordinate.
subroutine diag_remap_init(remap_cs, coord_tuple, answers_2018)
  type(diag_remap_ctrl), intent(inout) :: remap_cs !< Diag remapping control structure
  character(len=*),      intent(in)    :: coord_tuple !< A string in form of
                                                      !! MODULE_SUFFIX PARAMETER_SUFFIX COORDINATE_NAME
  logical,               intent(in)    :: answers_2018 !< If true, use the order of arithmetic and expressions
                                                      !! for remapping that recover the answers from the end of 2018.
                                                      !! Otherwise, use more robust forms of the same expressions.

  remap_cs%diag_module_suffix = trim(extractWord(coord_tuple, 1))
  remap_cs%diag_coord_name = trim(extractWord(coord_tuple, 2))
  remap_cs%vertical_coord_name = trim(extractWord(coord_tuple, 3))
  remap_cs%vertical_coord = coordinateMode(remap_cs%vertical_coord_name)
  remap_cs%configured = .false.
  remap_cs%initialized = .false.
  remap_cs%used = .false.
  remap_cs%answers_2018 = answers_2018
  remap_cs%nz = 0

end subroutine diag_remap_init

!> De-init a diagnostic remapping type.
!! Free allocated memory.
subroutine diag_remap_end(remap_cs)
  type(diag_remap_ctrl), intent(inout) :: remap_cs !< Diag remapping control structure

  if (allocated(remap_cs%h)) deallocate(remap_cs%h)
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
subroutine diag_remap_configure_axes(remap_cs, GV, US, param_file)
  type(diag_remap_ctrl),   intent(inout) :: remap_cs !< Diag remap control structure
  type(verticalGrid_type), intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file structure
  ! Local variables
  integer :: nzi(4), nzl(4), k
  character(len=200) :: inputdir, string, filename, int_varname, layer_varname
  character(len=40)  :: mod  = "MOM_diag_remap" ! This module's name.
  character(len=8)   :: units, expected_units
  character(len=34)  :: longname, string2

  character(len=256) :: err_msg
  logical :: ierr

  real, allocatable, dimension(:) :: interfaces, layers

  call initialize_regridding(remap_cs%regrid_cs, GV, US, GV%max_depth, param_file, mod, &
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

  remap_cs%interface_axes_id = MOM_diag_axis_init(lowercase(trim(remap_cs%diag_coord_name))//'_i', &
                                              interfaces, trim(units), 'z', &
                                              trim(longname)//' at interface', direction=-1)
  remap_cs%layer_axes_id = MOM_diag_axis_init(lowercase(trim(remap_cs%diag_coord_name))//'_l', &
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
  type(diag_remap_ctrl), intent(in) :: remap_cs !< Diagnostic coordinate control structure
  logical :: diag_remap_axes_configured

  diag_remap_axes_configured = remap_cs%configured

end function

!> Build/update target vertical grids for diagnostic remapping.
!! \note The target grids need to be updated whenever sea surface
!! height or layer thicknesses changes. In the case of density-based
!! coordinates then technically we should also regenerate the
!! target grid whenever T/S change.
subroutine diag_remap_update(remap_cs, G, GV, US, h, T, S, eqn_of_state, h_target)
  type(diag_remap_ctrl),   intent(inout) :: remap_cs !< Diagnostic coordinate control structure
  type(ocean_grid_type),   pointer    :: G  !< The ocean's grid type
  type(verticalGrid_type), intent(in) :: GV !< ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  real, dimension(:,:,:),  intent(in) :: h  !< New thickness [H ~> m or kg m-2]
  real, dimension(:,:,:),  intent(in) :: T  !< New temperatures [degC]
  real, dimension(:,:,:),  intent(in) :: S  !< New salinities [ppt]
  type(EOS_type),          pointer    :: eqn_of_state !< A pointer to the equation of state
  real, dimension(:,:,:),  intent(inout) :: h_target  !< The new diagnostic thicknesses [H ~> m or kg m-2]

  ! Local variables
  real, dimension(remap_cs%nz + 1) :: zInterfaces ! Interface positions [H ~> m or kg m-2]
  real :: h_neglect, h_neglect_edge ! Negligible thicknesses [H ~> m or kg m-2]
  integer :: i, j, k, nz

  ! Note that coordinateMode('LAYER') is never 'configured' so will
  ! always return here.
  if (.not. remap_cs%configured) then
    return
  endif

  if (.not.remap_cs%answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif
  nz = remap_cs%nz

  if (.not. remap_cs%initialized) then
    ! Initialize remapping and regridding on the first call
    call initialize_remapping(remap_cs%remap_cs, 'PPM_IH4', boundary_extrapolation=.false., &
                              answers_2018=remap_cs%answers_2018)
    remap_cs%initialized = .true.
  endif

  ! Calculate remapping thicknesses for different target grids based on
  ! nominal/target interface locations. This happens for every call on the
  ! assumption that h, T, S has changed.
  do j=G%jsc-1, G%jec+1 ; do i=G%isc-1, G%iec+1
    if (G%mask2dT(i,j)==0.) then
      h_target(i,j,:) = 0.
      cycle
    endif

    if (remap_cs%vertical_coord == coordinateMode('ZSTAR')) then
      call build_zstar_column(get_zlike_CS(remap_cs%regrid_cs), &
                              GV%Z_to_H*G%bathyT(i,j), sum(h(i,j,:)), &
                              zInterfaces, zScale=GV%Z_to_H)
    elseif (remap_cs%vertical_coord == coordinateMode('SIGMA')) then
      call build_sigma_column(get_sigma_CS(remap_cs%regrid_cs), &
                              GV%Z_to_H*G%bathyT(i,j), sum(h(i,j,:)), zInterfaces)
    elseif (remap_cs%vertical_coord == coordinateMode('RHO')) then
      call build_rho_column(get_rho_CS(remap_cs%regrid_cs), GV%ke, &
                            GV%Z_to_H*G%bathyT(i,j), h(i,j,:), T(i,j,:), S(i,j,:), &
                            eqn_of_state, zInterfaces, h_neglect, h_neglect_edge)
    elseif (remap_cs%vertical_coord == coordinateMode('SLIGHT')) then
!     call build_slight_column(remap_cs%regrid_cs,remap_cs%remap_cs, nz, &
!                           GV%Z_to_H*G%bathyT(i,j), sum(h(i,j,:)), zInterfaces)
      call MOM_error(FATAL,"diag_remap_update: SLIGHT coordinate not coded for diagnostics yet!")
    elseif (remap_cs%vertical_coord == coordinateMode('HYCOM1')) then
!     call build_hycom1_column(remap_cs%regrid_cs, nz, &
!                           GV%Z_to_H*G%bathyT(i,j), sum(h(i,j,:)), zInterfaces)
      call MOM_error(FATAL,"diag_remap_update: HYCOM1 coordinate not coded for diagnostics yet!")
    endif
    do k = 1,nz
      h_target(i,j,k) = zInterfaces(k) - zInterfaces(k+1)
    enddo
  enddo ; enddo

end subroutine diag_remap_update

!> Remap diagnostic field to alternative vertical grid.
subroutine diag_remap_do_remap(remap_cs, G, GV, h, staggered_in_x, staggered_in_y, &
                               mask, field, remapped_field)
  type(diag_remap_ctrl),   intent(in) :: remap_cs !< Diagnostic coodinate control structure
  type(ocean_grid_type),   intent(in) :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV !< ocean vertical grid structure
  real, dimension(:,:,:),  intent(in) :: h  !< The current thicknesses [H ~> m or kg m-2]
  logical,                 intent(in) :: staggered_in_x !< True is the x-axis location is at u or q points
  logical,                 intent(in) :: staggered_in_y !< True is the y-axis location is at v or q points
  real, dimension(:,:,:),  pointer    :: mask !< A mask for the field [nondim]
  real, dimension(:,:,:),  intent(in) :: field(:,:,:) !< The diagnostic field to be remapped [A]
  real, dimension(:,:,:),  intent(inout) :: remapped_field !< Field remapped to new coordinate [A]
  ! Local variables
  real, dimension(remap_cs%nz) :: h_dest ! Destination thicknesses [H ~> m or kg m-2]
  real, dimension(size(h,3)) :: h_src    ! A column of source thicknesses [H ~> m or kg m-2]
  real :: h_neglect, h_neglect_edge ! Negligible thicknesses [H ~> m or kg m-2]
  integer :: nz_src, nz_dest
  integer :: i, j, k                !< Grid index
  integer :: i1, j1                 !< 1-based index
  integer :: i_lo, i_hi, j_lo, j_hi !< (uv->h) interpolation indices
  integer :: shift                  !< Symmetric offset for 1-based indexing

  call assert(remap_cs%initialized, 'diag_remap_do_remap: remap_cs not initialized.')
  call assert(size(field, 3) == size(h, 3), &
              'diag_remap_do_remap: Remap field and thickness z-axes do not match.')

  if (.not.remap_cs%answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif

  nz_src = size(field,3)
  nz_dest = remap_cs%nz
  remapped_field(:,:,:) = 0.

  ! Symmetric grid offset under 1-based indexing; see header for details.
  shift = 0 ; if (G%symmetric) shift = 1

  if (staggered_in_x .and. .not. staggered_in_y) then
    ! U-points
    do j=G%jsc, G%jec
      do I=G%iscB, G%iecB
        I1 = I - G%isdB + 1
        i_lo = I1 - shift; i_hi = i_lo + 1
        if (associated(mask)) then
          if (mask(I,j,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i_lo,j,:) + h(i_hi,j,:))
        h_dest(:) = 0.5 * (remap_cs%h(i_lo,j,:) + remap_cs%h(i_hi,j,:))
        call remapping_core_h(remap_cs%remap_cs, &
                              nz_src, h_src(:), field(I1,j,:), &
                              nz_dest, h_dest(:), remapped_field(I1,j,:), &
                              h_neglect, h_neglect_edge)
      enddo
    enddo
  elseif (staggered_in_y .and. .not. staggered_in_x) then
    ! V-points
    do J=G%jscB, G%jecB
      J1 = J - G%jsdB + 1
      j_lo = J1 - shift; j_hi = j_lo + 1
      do i=G%isc, G%iec
        if (associated(mask)) then
          if (mask(i,J,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i,j_lo,:) + h(i,j_hi,:))
        h_dest(:) = 0.5 * (remap_cs%h(i,j_lo,:) + remap_cs%h(i,j_hi,:))
        call remapping_core_h(remap_cs%remap_cs, &
                              nz_src, h_src(:), field(i,J1,:), &
                              nz_dest, h_dest(:), remapped_field(i,J1,:), &
                              h_neglect, h_neglect_edge)
      enddo
    enddo
  elseif ((.not. staggered_in_x) .and. (.not. staggered_in_y)) then
    ! H-points
    do j=G%jsc, G%jec
      do i=G%isc, G%iec
        if (associated(mask)) then
          if (mask(i,j,1) == 0.) cycle
        endif
        h_src(:) = h(i,j,:)
        h_dest(:) = remap_cs%h(i,j,:)
        call remapping_core_h(remap_cs%remap_cs, &
                              nz_src, h_src(:), field(i,j,:), &
                              nz_dest, h_dest(:), remapped_field(i,j,:), &
                              h_neglect, h_neglect_edge)
      enddo
    enddo
  else
    call assert(.false., 'diag_remap_do_remap: Unsupported axis combination')
  endif

end subroutine diag_remap_do_remap

!> Calculate masks for target grid
subroutine diag_remap_calc_hmask(remap_cs, G, mask)
  type(diag_remap_ctrl),  intent(in)  :: remap_cs !< Diagnostic coodinate control structure
  type(ocean_grid_type),  intent(in)  :: G    !< Ocean grid structure
  real, dimension(:,:,:), intent(out) :: mask !< h-point mask for target grid [nondim]
  ! Local variables
  real, dimension(remap_cs%nz) :: h_dest ! Destination thicknesses [H ~> m or kg m-2]
  integer :: i, j, k
  logical :: mask_vanished_layers
  real :: h_tot      ! Sum of all thicknesses [H ~> m or kg m-2]
  real :: h_err      ! An estimate of a negligible thickness [H ~> m or kg m-2]

  call assert(remap_cs%initialized, 'diag_remap_calc_hmask: remap_cs not initialized.')

  ! Only z*-like diagnostic coordinates should have a 3d mask
  mask_vanished_layers = (remap_cs%vertical_coord == coordinateMode('ZSTAR'))
  mask(:,:,:) = 0.

  do j=G%jsc-1, G%jec+1 ; do i=G%isc-1, G%iec+1
    if (G%mask2dT(i,j)>0.) then
      if (mask_vanished_layers) then
        h_dest(:) = remap_cs%h(i,j,:)
        h_tot = 0.
        h_err = 0.
        do k=1, remap_cs%nz
          h_tot = h_tot + h_dest(k)
          ! This is an overestimate of how thick a vanished layer might be, that
          ! appears due to round-off.
          h_err = h_err + epsilon(h_tot) * h_tot
          ! Mask out vanished layers
          if (h_dest(k)<=8.*h_err) then
            mask(i,j,k) = 0.
          else
            mask(i,j,k) = 1.
          endif
        enddo
      else ! all layers might contain data
        mask(i,j,:) = 1.
      endif
    endif
  enddo ; enddo

end subroutine diag_remap_calc_hmask

!> Vertically re-grid an already vertically-integrated diagnostic field to alternative vertical grid.
subroutine vertically_reintegrate_diag_field(remap_cs, G, h, h_target, staggered_in_x, staggered_in_y, &
                                             mask, field, reintegrated_field)
  type(diag_remap_ctrl),  intent(in) :: remap_cs !< Diagnostic coodinate control structure
  type(ocean_grid_type),  intent(in) :: G        !< Ocean grid structure
  real, dimension(:,:,:), intent(in) :: h        !< The thicknesses of the source grid [H ~> m or kg m-2]
  real, dimension(:,:,:), intent(in) :: h_target !< The thicknesses of the target grid [H ~> m or kg m-2]
  logical,                intent(in) :: staggered_in_x !< True is the x-axis location is at u or q points
  logical,                intent(in) :: staggered_in_y !< True is the y-axis location is at v or q points
  real, dimension(:,:,:), pointer    :: mask     !< A mask for the field [nondim]
  real, dimension(:,:,:), intent(in) :: field    !<  The diagnostic field to be remapped [A]
  real, dimension(:,:,:), intent(inout) :: reintegrated_field !< Field argument remapped to alternative coordinate [A]
  ! Local variables
  real, dimension(remap_cs%nz) :: h_dest ! Destination thicknesses [H ~> m or kg m-2]
  real, dimension(size(h,3)) :: h_src    ! A column of source thicknesses [H ~> m or kg m-2]
  integer :: nz_src, nz_dest
  integer :: i, j, k                !< Grid index
  integer :: i1, j1                 !< 1-based index
  integer :: i_lo, i_hi, j_lo, j_hi !< (uv->h) interpolation indices
  integer :: shift                  !< Symmetric offset for 1-based indexing

  call assert(remap_cs%initialized, 'vertically_reintegrate_diag_field: remap_cs not initialized.')
  call assert(size(field, 3) == size(h, 3), &
              'vertically_reintegrate_diag_field: Remap field and thickness z-axes do not match.')

  nz_src = size(field,3)
  nz_dest = remap_cs%nz
  reintegrated_field(:,:,:) = 0.

  ! Symmetric grid offset under 1-based indexing; see header for details.
  shift = 0 ; if (G%symmetric) shift = 1

  if (staggered_in_x .and. .not. staggered_in_y) then
    ! U-points
    do j=G%jsc, G%jec
      do I=G%iscB, G%iecB
        I1 = I - G%isdB + 1
        i_lo = I1 - shift; i_hi = i_lo + 1
        if (associated(mask)) then
          if (mask(I,j,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i_lo,j,:) + h(i_hi,j,:))
        h_dest(:) = 0.5 * (h_target(i_lo,j,:) + h_target(i_hi,j,:))
        call reintegrate_column(nz_src, h_src, field(I1,j,:), &
                                nz_dest, h_dest, 0., reintegrated_field(I1,j,:))
      enddo
    enddo
  elseif (staggered_in_y .and. .not. staggered_in_x) then
    ! V-points
    do J=G%jscB, G%jecB
      J1 = J - G%jsdB + 1
      j_lo = J1 - shift; j_hi = j_lo + 1
      do i=G%isc, G%iec
        if (associated(mask)) then
          if (mask(i,J,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i,j_lo,:) + h(i,j_hi,:))
        h_dest(:) = 0.5 * (h_target(i,j_lo,:) + h_target(i,j_hi,:))
        call reintegrate_column(nz_src, h_src, field(i,J1,:), &
                                nz_dest, h_dest, 0., reintegrated_field(i,J1,:))
      enddo
    enddo
  elseif ((.not. staggered_in_x) .and. (.not. staggered_in_y)) then
    ! H-points
    do j=G%jsc, G%jec
      do i=G%isc, G%iec
        if (associated(mask)) then
          if (mask(i,j,1) == 0.) cycle
        endif
        h_src(:) = h(i,j,:)
        h_dest(:) = h_target(i,j,:)
        call reintegrate_column(nz_src, h_src, field(i,j,:), &
                                nz_dest, h_dest, 0., reintegrated_field(i,j,:))
      enddo
    enddo
  else
    call assert(.false., 'vertically_reintegrate_diag_field: Q point remapping is not coded yet.')
  endif

end subroutine vertically_reintegrate_diag_field

!> Vertically interpolate diagnostic field to alternative vertical grid.
subroutine vertically_interpolate_diag_field(remap_cs, G, h, staggered_in_x, staggered_in_y, &
                                             mask, field, interpolated_field)
  type(diag_remap_ctrl),  intent(in) :: remap_cs !< Diagnostic coodinate control structure
  type(ocean_grid_type),  intent(in) :: G   !< Ocean grid structure
  real, dimension(:,:,:), intent(in) :: h   !< The current thicknesses [H ~> m or kg m-2]
  logical,                intent(in) :: staggered_in_x !< True is the x-axis location is at u or q points
  logical,                intent(in) :: staggered_in_y !< True is the y-axis location is at v or q points
  real, dimension(:,:,:), pointer    :: mask !< A mask for the field [nondim]
  real, dimension(:,:,:), intent(in) :: field !<  The diagnostic field to be remapped [A]
  real, dimension(:,:,:), intent(inout) :: interpolated_field !< Field argument remapped to alternative coordinate [A]
  ! Local variables
  real, dimension(remap_cs%nz) :: h_dest ! Destination thicknesses [H ~> m or kg m-2]
  real, dimension(size(h,3)) :: h_src    ! A column of source thicknesses [H ~> m or kg m-2]
  integer :: nz_src, nz_dest
  integer :: i, j, k                !< Grid index
  integer :: i1, j1                 !< 1-based index
  integer :: i_lo, i_hi, j_lo, j_hi !< (uv->h) interpolation indices
  integer :: shift                  !< Symmetric offset for 1-based indexing

  call assert(remap_cs%initialized, 'vertically_interpolate_diag_field: remap_cs not initialized.')
  call assert(size(field, 3) == size(h, 3)+1, &
              'vertically_interpolate_diag_field: Remap field and thickness z-axes do not match.')

  interpolated_field(:,:,:) = 0.

  nz_src = size(h,3)
  nz_dest = remap_cs%nz

  ! Symmetric grid offset under 1-based indexing; see header for details.
  shift = 0 ; if (G%symmetric) shift = 1

  if (staggered_in_x .and. .not. staggered_in_y) then
    ! U-points
    do j=G%jsc, G%jec
      do I=G%iscB, G%iecB
        I1 = I - G%isdB + 1
        i_lo = I1 - shift; i_hi = i_lo + 1
        if (associated(mask)) then
          if (mask(I,j,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i_lo,j,:) + h(i_hi,j,:))
        h_dest(:) = 0.5 * (remap_cs%h(i_lo,j,:) + remap_cs%h(i_hi,j,:))
        call interpolate_column(nz_src, h_src, field(I1,j,:), &
                                nz_dest, h_dest, 0., interpolated_field(I1,j,:))
      enddo
    enddo
  elseif (staggered_in_y .and. .not. staggered_in_x) then
    ! V-points
    do J=G%jscB, G%jecB
      J1 = J - G%jsdB + 1
      j_lo = J1 - shift; j_hi = j_lo + 1
      do i=G%isc, G%iec
        if (associated(mask)) then
          if (mask(i,J,1) == 0.) cycle
        endif
        h_src(:) = 0.5 * (h(i,j_lo,:) + h(i,j_hi,:))
        h_dest(:) = 0.5 * (remap_cs%h(i,j_lo,:) + remap_cs%h(i,j_hi,:))
        call interpolate_column(nz_src, h_src, field(i,J1,:), &
                                nz_dest, h_dest, 0., interpolated_field(i,J1,:))
      enddo
    enddo
  elseif ((.not. staggered_in_x) .and. (.not. staggered_in_y)) then
    ! H-points
    do j=G%jsc, G%jec
      do i=G%isc, G%iec
        if (associated(mask)) then
          if (mask(i,j,1) == 0.) cycle
        endif
        h_src(:) = h(i,j,:)
        h_dest(:) = remap_cs%h(i,j,:)
        call interpolate_column(nz_src, h_src, field(i,j,:), &
                                nz_dest, h_dest, 0., interpolated_field(i,j,:))
      enddo
    enddo
  else
    call assert(.false., 'vertically_interpolate_diag_field: Q point remapping is not coded yet.')
  endif

end subroutine vertically_interpolate_diag_field

!> Horizontally average field
subroutine horizontally_average_diag_field(G, GV, h, staggered_in_x, staggered_in_y, &
                                           is_layer, is_extensive, &
                                           field, averaged_field, &
                                           averaged_mask)
  type(ocean_grid_type),  intent(in) :: G !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV !< The ocean vertical grid structure
  real, dimension(:,:,:), intent(in) :: h !< The current thicknesses [H ~> m or kg m-2]
  logical,                intent(in) :: staggered_in_x !< True if the x-axis location is at u or q points
  logical,                intent(in) :: staggered_in_y !< True if the y-axis location is at v or q points
  logical,                intent(in) :: is_layer !< True if the z-axis location is at h points
  logical,                intent(in) :: is_extensive !< True if the z-direction is spatially integrated (over layers)
  real, dimension(:,:,:), intent(in) :: field !<  The diagnostic field to be remapped [A]
  real, dimension(:),  intent(inout) :: averaged_field !< Field argument horizontally averaged [A]
  logical, dimension(:), intent(inout) :: averaged_mask  !< Mask for horizontally averaged field [nondim]

  ! Local variables
  real, dimension(G%isc:G%iec, G%jsc:G%jec, size(field,3)) :: volume, stuff
  real, dimension(size(field, 3)) :: vol_sum, stuff_sum ! nz+1 is needed for interface averages
  type(EFP_type), dimension(2*size(field,3)) :: sums_EFP ! Sums of volume or stuff by layer
  real :: height  ! An average thickness attributed to an velocity point [H ~> m or kg m-2]
  integer :: i, j, k, nz
  integer :: i1, j1                 !< 1-based index

  nz = size(field, 3)

  ! TODO: These averages could potentially be modified to use the function in
  !       the MOM_spatial_means module.
  ! NOTE: Reproducible sums must be computed in the original MKS units

  if (staggered_in_x .and. .not. staggered_in_y) then
    if (is_layer) then
      ! U-points
      do k=1,nz
        vol_sum(k) = 0.
        stuff_sum(k) = 0.
        if (is_extensive) then
          do j=G%jsc, G%jec ; do I=G%isc, G%iec
            I1 = I - G%isdB + 1
            volume(I,j,k) = (G%US%L_to_m**2 * G%areaCu(I,j)) * G%mask2dCu(I,j)
            stuff(I,j,k) = volume(I,j,k) * field(I1,j,k)
          enddo ; enddo
        else ! Intensive
          do j=G%jsc, G%jec ; do I=G%isc, G%iec
            I1 = i - G%isdB + 1
            height = 0.5 * (h(i,j,k) + h(i+1,j,k))
            volume(I,j,k) = (G%US%L_to_m**2 * G%areaCu(I,j)) &
                * (GV%H_to_m * height) * G%mask2dCu(I,j)
            stuff(I,j,k) = volume(I,j,k) * field(I1,j,k)
          enddo ; enddo
        endif
      enddo
    else ! Interface
      do k=1,nz
        do j=G%jsc, G%jec ; do I=G%isc, G%iec
          I1 = I - G%isdB + 1
          volume(I,j,k) = (G%US%L_to_m**2 * G%areaCu(I,j)) * G%mask2dCu(I,j)
          stuff(I,j,k) = volume(I,j,k) * field(I1,j,k)
        enddo ; enddo
      enddo
    endif
  elseif (staggered_in_y .and. .not. staggered_in_x) then
    if (is_layer) then
      ! V-points
      do k=1,nz
        if (is_extensive) then
          do J=G%jsc, G%jec ; do i=G%isc, G%iec
            J1 = J - G%jsdB + 1
            volume(i,J,k) = (G%US%L_to_m**2 * G%areaCv(i,J)) * G%mask2dCv(i,J)
            stuff(i,J,k) = volume(i,J,k) * field(i,J1,k)
          enddo ; enddo
        else ! Intensive
          do J=G%jsc, G%jec ; do i=G%isc, G%iec
            J1 = J - G%jsdB + 1
            height = 0.5 * (h(i,j,k) + h(i,j+1,k))
            volume(i,J,k) = (G%US%L_to_m**2 * G%areaCv(i,J)) &
                * (GV%H_to_m * height) * G%mask2dCv(i,J)
            stuff(i,J,k) = volume(i,J,k) * field(i,J1,k)
          enddo ; enddo
        endif
      enddo
    else ! Interface
      do k=1,nz
        do J=G%jsc, G%jec ; do i=G%isc, G%iec
          J1 = J - G%jsdB + 1
          volume(i,J,k) = (G%US%L_to_m**2 * G%areaCv(i,J)) * G%mask2dCv(i,J)
          stuff(i,J,k) = volume(i,J,k) * field(i,J1,k)
        enddo ; enddo
      enddo
    endif
  elseif ((.not. staggered_in_x) .and. (.not. staggered_in_y)) then
    if (is_layer) then
      ! H-points
      do k=1,nz
        if (is_extensive) then
          do j=G%jsc, G%jec ; do i=G%isc, G%iec
            if (h(i,j,k) > 0.) then
              volume(i,j,k) = (G%US%L_to_m**2 * G%areaT(i,j)) * G%mask2dT(i,j)
              stuff(i,j,k) = volume(i,j,k) * field(i,j,k)
            else
              volume(i,j,k) = 0.
              stuff(i,j,k) = 0.
            endif
          enddo ; enddo
        else ! Intensive
          do j=G%jsc, G%jec ; do i=G%isc, G%iec
            volume(i,j,k) = (G%US%L_to_m**2 * G%areaT(i,j)) &
                * (GV%H_to_m * h(i,j,k)) * G%mask2dT(i,j)
            stuff(i,j,k) = volume(i,j,k) * field(i,j,k)
          enddo ; enddo
        endif
      enddo
    else ! Interface
      do k=1,nz
        do j=G%jsc, G%jec ; do i=G%isc, G%iec
          volume(i,j,k) = (G%US%L_to_m**2 * G%areaT(i,j)) * G%mask2dT(i,j)
          stuff(i,j,k) = volume(i,j,k) * field(i,j,k)
        enddo ; enddo
      enddo
    endif
  else
    call assert(.false., 'horizontally_average_diag_field: Q point averaging is not coded yet.')
  endif

  ! Packing the sums into a single array with a single call to sum across PEs saves reduces
  ! the costs of communication.
  do k=1,nz
    sums_EFP(2*k-1) = reproducing_sum_EFP(volume(:,:,k), only_on_PE=.true.)
    sums_EFP(2*k)   = reproducing_sum_EFP(stuff(:,:,k), only_on_PE=.true.)
  enddo
  call EFP_sum_across_PEs(sums_EFP, 2*nz)
  do k=1,nz
    vol_sum(k) = EFP_to_real(sums_EFP(2*k-1))
    stuff_sum(k) = EFP_to_real(sums_EFP(2*k))
  enddo

  averaged_mask(:) = .true.
  do k=1,nz
    if (vol_sum(k) > 0.) then
      averaged_field(k) = stuff_sum(k) / vol_sum(k)
    else
      averaged_field(k) = 0.
      averaged_mask(k) = .false.
    endif
  enddo

end subroutine horizontally_average_diag_field

end module MOM_diag_remap
