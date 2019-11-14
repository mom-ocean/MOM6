!> Regrid columns for the HyCOM coordinate
module coord_hycom

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_EOS,           only : EOS_type, calculate_density
use regrid_interp,     only : interp_CS_type, build_and_interpolate_grid

implicit none ; private

!> Control structure containing required parameters for the HyCOM coordinate
type, public :: hycom_CS ; private

  !> Number of layers/levels in generated grid
  integer :: nk

  !> Nominal near-surface resolution
  real, allocatable, dimension(:) :: coordinateResolution

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:) :: target_density

  !> Density scaling factor [R m3 kg-1 ~> 1]
  real :: kg_m3_to_R

  !> Maximum depths of interfaces
  real, allocatable, dimension(:) :: max_interface_depths

  !> Maximum thicknesses of layers
  real, allocatable, dimension(:) :: max_layer_thickness

  !> Interpolation control structure
  type(interp_CS_type) :: interp_CS
end type hycom_CS

public init_coord_hycom, set_hycom_params, build_hycom1_column, end_coord_hycom

contains

!> Initialise a hycom_CS with pointers to parameters
subroutine init_coord_hycom(CS, nk, coordinateResolution, target_density, interp_CS, rho_scale)
  type(hycom_CS),       pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,              intent(in) :: nk !< Number of layers in generated grid
  real, dimension(nk),  intent(in) :: coordinateResolution !< Nominal near-surface resolution [m]
  real, dimension(nk+1),intent(in) :: target_density !< Interface target densities [R ~> kg m-3]
  type(interp_CS_type), intent(in) :: interp_CS !< Controls for interpolation
  real,       optional, intent(in) :: rho_scale !< A dimensional scaling factor for target_density

  if (associated(CS)) call MOM_error(FATAL, "init_coord_hycom: CS already associated!")
  allocate(CS)
  allocate(CS%coordinateResolution(nk))
  allocate(CS%target_density(nk+1))

  CS%nk                      = nk
  CS%coordinateResolution(:) = coordinateResolution(:)
  CS%target_density(:)       = target_density(:)
  CS%interp_CS               = interp_CS
  CS%kg_m3_to_R = 1.0 ; if (present(rho_scale)) CS%kg_m3_to_R = rho_scale

end subroutine init_coord_hycom

!> This subroutine deallocates memory in the control structure for the coord_hycom module
subroutine end_coord_hycom(CS)
  type(hycom_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
  if (.not. associated(CS)) return
  deallocate(CS%coordinateResolution)
  deallocate(CS%target_density)
  if (allocated(CS%max_interface_depths)) deallocate(CS%max_interface_depths)
  if (allocated(CS%max_layer_thickness)) deallocate(CS%max_layer_thickness)
  deallocate(CS)
end subroutine end_coord_hycom

!> This subroutine can be used to set the parameters for the coord_hycom module
subroutine set_hycom_params(CS, max_interface_depths, max_layer_thickness, interp_CS)
  type(hycom_CS),                 pointer    :: CS !< Coordinate control structure
  real, dimension(:),   optional, intent(in) :: max_interface_depths !< Maximum depths of interfaces in m
  real, dimension(:),   optional, intent(in) :: max_layer_thickness  !< Maximum thicknesses of layers in m
  type(interp_CS_type), optional, intent(in) :: interp_CS !< Controls for interpolation

  if (.not. associated(CS)) call MOM_error(FATAL, "set_hycom_params: CS not associated")

  if (present(max_interface_depths)) then
    if (size(max_interface_depths) /= CS%nk+1) &
      call MOM_error(FATAL, "set_hycom_params: max_interface_depths inconsistent size")
    allocate(CS%max_interface_depths(CS%nk+1))
    CS%max_interface_depths(:) = max_interface_depths(:)
  endif

  if (present(max_layer_thickness)) then
    if (size(max_layer_thickness) /= CS%nk) &
      call MOM_error(FATAL, "set_hycom_params: max_layer_thickness inconsistent size")
    allocate(CS%max_layer_thickness(CS%nk))
    CS%max_layer_thickness(:) = max_layer_thickness(:)
  endif

  if (present(interp_CS)) CS%interp_CS = interp_CS
end subroutine set_hycom_params

!> Build a HyCOM coordinate column
subroutine build_hycom1_column(CS, eqn_of_state, nz, depth, h, T, S, p_col, &
                               z_col, z_col_new, zScale, h_neglect, h_neglect_edge)
  type(hycom_CS),        intent(in)    :: CS !< Coordinate control structure
  type(EOS_type),        pointer       :: eqn_of_state !< Equation of state structure
  integer,               intent(in)    :: nz !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive [H ~> m or kg m-2])
  real, dimension(nz),   intent(in)    :: T !< Temperature of column [degC]
  real, dimension(nz),   intent(in)    :: S !< Salinity of column [ppt]
  real, dimension(nz),   intent(in)    :: h  !< Layer thicknesses, in [m] or [H ~> m or kg m-2]
  real, dimension(nz),   intent(in)    :: p_col !< Layer pressure [Pa]
  real, dimension(nz+1), intent(in)    :: z_col !< Interface positions relative to the surface [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: z_col_new !< Absolute positions of interfaces
  real, optional,        intent(in)    :: zScale !< Scaling factor from the input thicknesses in [m]
                                                 !! to desired units for zInterface, perhaps m_to_H.
  real,        optional, intent(in)    :: h_neglect !< A negligibly small width for the
                                             !! purpose of cell reconstructions
                                             !! in the same units as h.
  real,        optional, intent(in)    :: h_neglect_edge !< A negligibly small width
                                             !! for the purpose of edge value calculations
                                             !! in the same units as h0.

  ! Local variables
  integer   :: k
  real, dimension(nz) :: rho_col ! Layer densities in a column [R ~> kg m-3]
  real, dimension(CS%nk) :: h_col_new ! New layer thicknesses
  real :: z_scale
  real :: stretching ! z* stretching, converts z* to z.
  real :: nominal_z ! Nominal depth of interface when using z* [Z ~> m]
  real :: hNew
  logical :: maximum_depths_set ! If true, the maximum depths of interface have been set.
  logical :: maximum_h_set      ! If true, the maximum layer thicknesses have been set.

  maximum_depths_set = allocated(CS%max_interface_depths)
  maximum_h_set = allocated(CS%max_layer_thickness)

  z_scale = 1.0 ; if (present(zScale)) z_scale = zScale

  ! Work bottom recording potential density
  call calculate_density(T, S, p_col, rho_col, 1, nz, eqn_of_state, scale=CS%kg_m3_to_R)
  ! This ensures the potential density profile is monotonic
  ! although not necessarily single valued.
  do k = nz-1, 1, -1
    rho_col(k) = min( rho_col(k), rho_col(k+1) )
  enddo

  ! Interpolates for the target interface position with the rho_col profile
  ! Based on global density profile, interpolate to generate a new grid
  call build_and_interpolate_grid(CS%interp_CS, rho_col, nz, h(:), z_col, &
           CS%target_density, CS%nk, h_col_new, z_col_new, h_neglect, h_neglect_edge)

  ! Sweep down the interfaces and make sure that the interface is at least
  ! as deep as a nominal target z* grid
  nominal_z = 0.
  stretching = z_col(nz+1) / depth ! Stretches z* to z
  do k = 2, CS%nk+1
    nominal_z = nominal_z + (z_scale * CS%coordinateResolution(k-1)) * stretching
    z_col_new(k) = max( z_col_new(k), nominal_z )
    z_col_new(k) = min( z_col_new(k), z_col(nz+1) )
  enddo

  if (maximum_depths_set .and. maximum_h_set) then ; do k=2,CS%nk
    ! The loop bounds are 2 & nz so the top and bottom interfaces do not move.
    ! Recall that z_col_new is positive downward.
    z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K), &
                       z_col_new(K-1) + CS%max_layer_thickness(k-1))
  enddo ; elseif (maximum_depths_set) then ; do K=2,CS%nk
    z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K))
  enddo ; elseif (maximum_h_set) then ; do k=2,CS%nk
    z_col_new(K) = min(z_col_new(K), z_col_new(K-1) + CS%max_layer_thickness(k-1))
  enddo ; endif
end subroutine build_hycom1_column

end module coord_hycom
