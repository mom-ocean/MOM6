module coord_sigma

use MOM_error_handler, only : MOM_error, FATAL

implicit none ; private

type, public :: sigma_CS
  private

  real,               pointer :: min_thickness
  real, dimension(:), pointer :: coordinateResolution
end type sigma_CS

public init_coord_sigma, build_sigma_column

contains

subroutine init_coord_sigma(CS, min_thickness, coordinateResolution)
  type(sigma_CS),     pointer :: CS
  real,               target  :: min_thickness
  real, dimension(:), target  :: coordinateResolution

  if (associated(CS)) call MOM_error(FATAL, "init_coord_sigma: CS already associated!")
  allocate(CS)

  CS%min_thickness => min_thickness
  CS%coordinateResolution => coordinateResolution
end subroutine init_coord_sigma

subroutine build_sigma_column(CS, nz, depth, totalThickness, zInterface)
  type(sigma_CS),        intent(in)    :: CS !< Coordinate control structure
  integer,               intent(in)    :: nz !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive in m)
  real,                  intent(in)    :: totalThickness !< Column thickness (positive in m)
  real, dimension(nz+1), intent(inout) :: zInterface !< Absolute positions of interfaces

  ! Local variables
  integer :: k

  zInterface(nz+1) = -depth
  do k = nz,1,-1
    zInterface(k) = zInterface(k+1) + (totalThickness * CS%coordinateResolution(k))
    ! Adjust interface position to accomodate inflating layers
    ! without disturbing the interface above
    if (zInterface(k) < (zInterface(k+1) + CS%min_thickness)) then
      zInterface(k) = zInterface(k+1) + CS%min_thickness
    endif
  enddo
end subroutine build_sigma_column

end module coord_sigma
