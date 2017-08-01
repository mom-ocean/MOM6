module supercritical_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE, OBC_segment_type
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_time_manager,   only : time_type, set_time, time_type_to_real

implicit none ; private

#include <MOM_memory.h>

public supercritical_set_OBC_data
public supercritical_initialize_topography

! This include declares and sets the variable "version".
#include "version_variable.h"

contains

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine supercritical_set_OBC_data(OBC, G, param_file)
  type(ocean_OBC_type),   pointer    :: OBC  !< This open boundary condition type specifies
                                             !! whether, where, and what open boundary
                                             !! conditions are used.
  type(ocean_grid_type),  intent(in) :: G    !< The ocean's grid structure.
  type(param_file_type),  intent(in) :: param_file !< Parameter file structure
  ! Local variables
  character(len=40)  :: mdl = "supercritical_set_OBC_data" ! This subroutine's name.
  real :: zonal_flow
  integer :: i, j, k, l
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment ! pointer to segment type list

  if (.not.associated(OBC)) call MOM_error(FATAL, 'supercritical_initialization.F90: '// &
        'supercritical_set_OBC_data() was called but OBC type was not initialized!')

  call get_param(param_file, mdl, "SUPERCRITICAL_ZONAL_FLOW", zonal_flow, &
                 "Constant zonal flow imposed at upstream open boundary.", &
                 units="m/s", default=8.57)

  do l=1, OBC%number_of_segments
    segment => OBC%segment(l)
    if (.not. segment%on_pe) cycle
    if (segment%gradient) cycle
    if (segment%oblique .and. .not. segment%nudged .and. .not. segment%Flather) cycle

    if (segment%is_E_or_W) then
      jsd = segment%HI%jsd ; jed = segment%HI%jed
      IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
      do k=1,G%ke
        do j=jsd,jed ; do I=IsdB,IedB
          if (segment%specified .or. segment%nudged) then
            segment%normal_vel(I,j,k) = zonal_flow
          endif
          if (segment%specified) then
            segment%normal_trans(I,j,k) = zonal_flow * G%dyCu(I,j)
          endif
        enddo ; enddo
      enddo
      do j=jsd,jed ; do I=IsdB,IedB
        segment%normal_vel_bt(I,j) = zonal_flow
      enddo ; enddo
    else
      isd = segment%HI%isd ; ied = segment%HI%ied
      JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
!     do k=1,G%ke
!       do J=JsdB,JedB ; do i=isd,ied
!         segment%normal_vel(i,J,k) = 0.0
!       enddo ; enddo
!     enddo
      do J=JsdB,JedB ; do i=isd,ied
        segment%normal_vel_bt(i,J) = 0.0
      enddo ; enddo
    endif
  enddo

end subroutine supercritical_set_OBC_data

!> This subroutine sets up the supercritical topography and land mask.
!! We were not able to get the shock wave to behave this way and are
!! now using an external file.
subroutine supercritical_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),           intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),            intent(in)  :: param_file !< Parameter file structure
  real,                             intent(in)  :: max_depth  !< Maximum depth of model in m
  ! Local variables
  character(len=40)  :: mdl = "supercritical_initialize_topography" ! This subroutine's name.
  real :: min_depth ! The minimum and maximum depths in m.
  real :: PI ! 3.1415...
  real :: coast_offset, coast_angle
  integer :: i, j

  call MOM_mesg("  supercritical_initialization.F90, supercritical_initialize_topography: setting topography", 5)

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)
  call get_param(param_file, mdl, "SUPERCRITICAL_COAST_OFFSET", coast_offset, &
                 "The distance along the southern boundary at which the coasts angles in.", &
                 units="km", default=10.0)
  call get_param(param_file, mdl, "SUPERCRITICAL_COAST_ANGLE", coast_angle, &
                 "The angle of the southern bondary beyond X=SUPERCRITICAL_COAST_OFFSET.", &
                 units="degrees", default=8.95)

  coast_angle = coast_angle * (atan(1.0)/45.) ! Convert to radians

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    D(i,j)=max_depth
    if ((G%geoLonT(i,j) > coast_offset).AND. &
        (atan2(G%geoLatT(i,j),G%geoLonT(i,j) - coast_offset) < coast_angle)) then
      D(i,j)=0.5*min_depth
    endif

    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
  enddo ; enddo

end subroutine supercritical_initialize_topography

!> \namespace supercritical_initialization
!!
!! The module configures the model for the "supercritical" experiment.
!! https://marine.rutgers.edu/po/index.php?model=test-problems&title=supercritical
end module supercritical_initialization
