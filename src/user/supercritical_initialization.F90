module supercritical_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_time_manager,   only : time_type, set_time, time_type_to_real

implicit none ; private

#include <MOM_memory.h>

public supercritical_initialize_topography
public supercritical_initialize_velocity
public supercritical_set_OBC_data

! This include declares and sets the variable "version".
#include "version_variable.h"

contains

! -----------------------------------------------------------------------------
!> This subroutine sets up the supercritical topography
subroutine supercritical_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),           intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),            intent(in)  :: param_file !< Parameter file structure
  real,                             intent(in)  :: max_depth  !< Maximum depth of model in m
  ! Local variables
  character(len=40)  :: mod = "supercritical_initialize_topography" ! This subroutine's name.
  real :: min_depth ! The minimum and maximum depths in m.
  real :: PI ! 3.1415...
  real :: coast_offset, coast_angle
  integer :: i, j

  call MOM_mesg("  supercritical_initialization.F90, supercritical_initialize_topography: setting topography", 5)

  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)
  call get_param(param_file, mod, "SUPERCRITICAL_COAST_OFFSET", coast_offset, &
                 "The distance alond the southern boundary at which the coasts angles in.", &
                 units="km", default=10.0)
  call get_param(param_file, mod, "SUPERCRITICAL_COAST_ANGLE", coast_angle, &
                 "The angle of the southern bondary beyond X=SUPERCRITICAL_COAST_OFFSET.", &
                 units="degrees", default=8.95)

  coast_angle = coast_angle * (atan(1.0)/45.) ! Convert to radians

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    D(i,j)=max_depth
    if ((G%geoLonT(i,j) > coast_offset).AND. &
        (atan2(G%geoLatT(i,j),G%geoLonT(i,j)-10.0) < coast_angle)) then
      D(i,j)=0.5*min_depth
    endif

    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
  enddo ; enddo

end subroutine supercritical_initialize_topography
! -----------------------------------------------------------------------------
!> Initialization of u and v in the supercritical test
subroutine supercritical_initialize_velocity(u, v, h, G, param_file)
  type(ocean_grid_type),                     intent(in)  :: G  !< Grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out) :: u  !< i-component of velocity [m/s]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out) :: v  !< j-component of velocity [m/s]
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h  !< Thickness [H]
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  ! Local variables
  character(len=40)  :: mod = "supercritical_initialize_topography" ! This subroutine's name.
  real :: zonal_flow
  integer :: i, j, k

  call get_param(param_file, mod, "SUPERCRITICAL_ZONAL_FLOW", zonal_flow, &
                 "The intial and imposed zonal flow.", units="m s-1", default=8.57)

  do k = 1, G%ke
    do j = G%jsc,G%jec ; do I = G%IscB,G%IecB
      u(I,j,k) = zonal_flow
    enddo ; enddo
    do J = G%JscB,G%JecB ; do I = G%isc,G%iec
      v(i,J,k) = 0.
    enddo ; enddo
  enddo

end subroutine supercritical_initialize_velocity
! -----------------------------------------------------------------------------
!> This subroutine sets the properties of flow at open boundary conditions.
subroutine supercritical_set_OBC_data(OBC, G, param_file)
  type(ocean_OBC_type),   pointer    :: OBC  !< This open boundary condition type specifies
                                             !! whether, where, and what open boundary
                                             !! conditions are used.
  type(ocean_grid_type),  intent(in) :: G    !< The ocean's grid structure.
  type(param_file_type),  intent(in) :: param_file !< Parameter file structure
  ! Local variables
  character(len=40)  :: mod = "supercritical_set_OBC_data" ! This subroutine's name.
  real :: zonal_flow
  integer :: i, j, k

  if (.not.associated(OBC)) call MOM_error(FATAL, 'supercritical_initialization.F90: '// &
        'supercritical_set_OBC_data() was called but OBC type was not initialized!')

  call get_param(param_file, mod, "SUPERCRITICAL_ZONAL_FLOW", zonal_flow, &
                 default=8.57, do_not_log=.true.) ! Do not log avoids maintaining identical descriptions

  allocate(OBC%u(G%IsdB:G%IedB,G%jsd:G%jed,G%ke)) ; OBC%u(:,:,:) = 0.0
  allocate(OBC%uh(G%IsdB:G%IedB,G%jsd:G%jed,G%ke)) ; OBC%uh(:,:,:) = 0.0

  do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
    if (OBC%OBC_segment_u(I,j)>0) then
      do k=1,G%ke
        OBC%u(I,j,k) = 8.57
        OBC%uh(I,j,k) = 8.57 * G%dyCu(I,j)
      enddo
      OBC%ubt_outer(I,j) = 8.57
    endif
  enddo ; enddo

end subroutine supercritical_set_OBC_data

!> \class supercritical_initialization
!!
!! The module configures the model for the "supercritical" experiment.
!! https://marine.rutgers.edu/po/index.php?model=test-problems&title=supercritical
end module supercritical_initialization
