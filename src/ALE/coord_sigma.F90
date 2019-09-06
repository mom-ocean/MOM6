!> Regrid columns for the sigma coordinate
module coord_sigma

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL

implicit none ; private

!> Control structure containing required parameters for the sigma coordinate
type, public :: sigma_CS ; private

  !> Number of levels
  integer :: nk

  !> Minimum thickness allowed for layers
  real :: min_thickness

  !> Target coordinate resolution, nondimensional
  real, allocatable, dimension(:) :: coordinateResolution
end type sigma_CS

public init_coord_sigma, set_sigma_params, build_sigma_column, end_coord_sigma

contains

!> Initialise a sigma_CS with pointers to parameters
subroutine init_coord_sigma(CS, nk, coordinateResolution)
  type(sigma_CS),     pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,            intent(in) :: nk !< Number of layers in the grid
  real, dimension(:), intent(in) :: coordinateResolution !< Nominal coordinate resolution [nondim]

  if (associated(CS)) call MOM_error(FATAL, "init_coord_sigma: CS already associated!")
  allocate(CS)
  allocate(CS%coordinateResolution(nk))

  CS%nk                   = nk
  CS%coordinateResolution = coordinateResolution
end subroutine init_coord_sigma

!> This subroutine deallocates memory in the control structure for the coord_sigma module
subroutine end_coord_sigma(CS)
  type(sigma_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
  if (.not. associated(CS)) return
  deallocate(CS%coordinateResolution)
  deallocate(CS)
end subroutine end_coord_sigma

!> This subroutine can be used to set the parameters for the coord_sigma module
subroutine set_sigma_params(CS, min_thickness)
  type(sigma_CS), pointer    :: CS !< Coordinate control structure
  real, optional, intent(in) :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]

  if (.not. associated(CS)) call MOM_error(FATAL, "set_sigma_params: CS not associated")

  if (present(min_thickness)) CS%min_thickness = min_thickness
end subroutine set_sigma_params


!> Build a sigma coordinate column
subroutine build_sigma_column(CS, depth, totalThickness, zInterface)
  type(sigma_CS),           intent(in)    :: CS !< Coordinate control structure
  real,                     intent(in)    :: depth !< Depth of ocean bottom (positive [H ~> m or kg m-2])
  real,                     intent(in)    :: totalThickness !< Column thickness (positive [H ~> m or kg m-2])
  real, dimension(CS%nk+1), intent(inout) :: zInterface !< Absolute positions of interfaces [H ~> m or kg m-2]

  ! Local variables
  integer :: k

  zInterface(CS%nk+1) = -depth
  do k = CS%nk,1,-1
    zInterface(k) = zInterface(k+1) + (totalThickness * CS%coordinateResolution(k))
    ! Adjust interface position to accommodate inflating layers
    ! without disturbing the interface above
    if (zInterface(k) < (zInterface(k+1) + CS%min_thickness)) then
      zInterface(k) = zInterface(k+1) + CS%min_thickness
    endif
  enddo
end subroutine build_sigma_column

end module coord_sigma
