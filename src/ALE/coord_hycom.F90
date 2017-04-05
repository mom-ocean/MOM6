module coord_hycom

use MOM_error_handler, only : MOM_error, FATAL
use MOM_EOS,           only : EOS_type, calculate_density
use regrid_interp,     only : interp_CS_type, build_and_interpolate_grid

implicit none ; private

type, public :: hycom_CS
  private

  type(interp_CS_type), pointer :: interp_CS
  real, dimension(:),   pointer :: target_density
  real, dimension(:),   pointer :: coordinateResolution
  real, dimension(:),   pointer :: max_interface_depths
  real, dimension(:),   pointer :: max_layer_thickness
end type hycom_CS

public init_coord_hycom, build_hycom1_column

contains

subroutine init_coord_hycom(CS, interp_CS, target_density, coordinateResolution, &
     max_interface_depths, max_layer_thickness)
  type(hycom_CS),       pointer :: CS
  type(interp_CS_type), target  :: interp_CS
  real, dimension(:),   target  :: target_density, coordinateResolution, &
       max_interface_depths, max_layer_thickness

  if (associated(CS)) call MOM_error(FATAL, "init_coord_hycom: CS already associated!")
  allocate(CS)

  CS%interp_CS            => interp_CS
  CS%target_density       => target_density
  CS%coordinateResolution => coordinateResolution
  CS%max_interface_depths => max_interface_depths
  CS%max_layer_thickness  => max_layer_thickness
end subroutine init_coord_hycom

subroutine build_hycom1_column(CS, eqn_of_state, nz, depth, h, T, S, p_col, z_col, z_col_new)
  type(hycom_CS),        intent(in)    :: CS !< Regridding control structure
  type(EOS_type),        pointer       :: eqn_of_state !< Equation of state structure
  integer,               intent(in)    :: nz !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive in H)
  real, dimension(nz),   intent(in)    :: T, S !< T and S for column
  real, dimension(nz),   intent(in)    :: h  !< Layer thicknesses, in m
  real, dimension(nz),   intent(in)    :: p_col !< Layer pressure in Pa
  real, dimension(nz+1), intent(in)    :: z_col ! Interface positions relative to the surface in H units (m or kg m-2)
  real, dimension(nz+1), intent(inout) :: z_col_new !< Absolute positions of interfaces

  ! Local variables
  integer   :: k
  real, dimension(nz) :: rho_col, h_col_new ! Layer quantities
  real :: stretching ! z* stretching, converts z* to z.
  real :: nominal_z ! Nominal depth of interface is using z* (m or Pa)
  real :: hNew
  logical :: maximum_depths_set ! If true, the maximum depths of interface have been set.
  logical :: maximum_h_set      ! If true, the maximum layer thicknesses have been set.

  maximum_depths_set = allocated(CS%max_interface_depths)
  maximum_h_set = allocated(CS%max_layer_thickness)

  ! Work bottom recording potential density
  call calculate_density(T, S, p_col, rho_col, 1, nz, eqn_of_state)
  ! This ensures the potential density profile is monotonic
  ! although not necessarily single valued.
  do k = nz-1, 1, -1
    rho_col(k) = min( rho_col(k), rho_col(k+1) )
  enddo

  ! Interpolates for the target interface position with the rho_col profile
  ! Based on global density profile, interpolate to generate a new grid
  call build_and_interpolate_grid(CS%interp_CS, rho_col, nz, h(:), z_col, &
       CS%target_density, nz, h_col_new, z_col_new)

  ! Sweep down the interfaces and make sure that the interface is at least
  ! as deep as a nominal target z* grid
  nominal_z = 0.
  stretching = z_col(nz+1) / depth ! Stretches z* to z
  do k = 2, nz+1
    nominal_z = nominal_z + CS%coordinateResolution(k-1) * stretching
    z_col_new(k) = max( z_col_new(k), nominal_z )
    z_col_new(k) = min( z_col_new(k), z_col(nz+1) )
  enddo

  if (maximum_depths_set .and. maximum_h_set) then ; do k=2,nz
    ! The loop bounds are 2 & nz so the top and bottom interfaces do not move.
    ! Recall that z_col_new is positive downward.
    z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K), &
                       z_col_new(K-1) + CS%max_layer_thickness(k-1))
  enddo ; elseif (maximum_depths_set) then ; do K=2,nz
    z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K))
  enddo ; elseif (maximum_h_set) then ; do k=2,nz
    z_col_new(K) = min(z_col_new(K), z_col_new(K-1) + CS%max_layer_thickness(k-1))
  enddo ; endif
end subroutine build_hycom1_column

end module coord_hycom
