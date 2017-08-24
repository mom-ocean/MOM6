!> Regrid columns for the sigma coordinate
module coord_sigma

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL

implicit none ; private

!> Control structure containing required parameters for the sigma coordinate
type, public :: sigma_CS
  private

  !> Number of levels
  integer :: nk

  !> Minimum thickness allowed for layers
  real :: min_thickness

  !> Target coordinate resolution
  real, allocatable, dimension(:) :: coordinateResolution
end type sigma_CS

public init_coord_sigma, set_sigma_params, build_sigma_column, end_coord_sigma

contains

!> Initialise a sigma_CS with pointers to parameters
subroutine init_coord_sigma(CS, nk, coordinateResolution)
  type(sigma_CS),     pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,            intent(in) :: nk
  real, dimension(:), intent(in) :: coordinateResolution

  if (associated(CS)) call MOM_error(FATAL, "init_coord_sigma: CS already associated!")
  allocate(CS)
  allocate(CS%coordinateResolution(nk))

  CS%nk                   = nk
  CS%coordinateResolution = coordinateResolution
end subroutine init_coord_sigma

subroutine end_coord_sigma(CS)
  type(sigma_CS), pointer :: CS

  ! nothing to do
  if (.not. associated(CS)) return
  deallocate(CS%coordinateResolution)
  deallocate(CS)
end subroutine end_coord_sigma

subroutine set_sigma_params(CS, min_thickness)
  type(sigma_CS), pointer    :: CS
  real, optional, intent(in) :: min_thickness

  if (.not. associated(CS)) call MOM_error(FATAL, "set_sigma_params: CS not associated")

  if (present(min_thickness)) CS%min_thickness = min_thickness
end subroutine set_sigma_params


!> Build a sigma coordinate column
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
