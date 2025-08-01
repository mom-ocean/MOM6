!> The "super critical" configuration
module supercritical_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type, OBC_segment_type, rotate_OBC_segment_direction
use MOM_open_boundary,  only : OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_time_manager,   only : time_type, time_type_to_real
use MOM_unit_scaling,   only : unit_scale_type
use MOM_verticalGrid,   only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public supercritical_set_OBC_data

contains

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine supercritical_set_OBC_data(OBC, G, GV, US, param_file)
  type(ocean_OBC_type),    pointer    :: OBC  !< This open boundary condition type specifies
                                              !! whether, where, and what open boundary
                                              !! conditions are used.
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file structure
  ! Local variables
  character(len=40)  :: mdl = "supercritical_set_OBC_data" ! This subroutine's name.
  real :: zonal_flow ! Inflow speed [L T-1 ~> m s-1]
  integer :: unrot_dir ! The unrotated direction of the segment
  integer :: turns    ! Number of index quarter turns
  integer :: i, j, k, l
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list

  if (.not.associated(OBC)) call MOM_error(FATAL, 'supercritical_initialization.F90: '// &
        'supercritical_set_OBC_data() was called but OBC type was not initialized!')

  call get_param(param_file, mdl, "SUPERCRITICAL_ZONAL_FLOW", zonal_flow, &
                 "Constant zonal flow imposed at upstream open boundary.", &
                 units="m/s", default=8.57, scale=US%m_s_to_L_T)

  turns = modulo(G%HI%turns, 4)

  do l=1, OBC%number_of_segments
    segment => OBC%segment(l)
    if (.not. segment%on_pe) cycle
    if (segment%gradient) cycle
    if (segment%oblique .and. .not. segment%nudged .and. .not. segment%Flather) cycle

   unrot_dir = segment%direction
   if (turns /= 0) unrot_dir = rotate_OBC_segment_direction(segment%direction, -turns)

    if ((unrot_dir == OBC_DIRECTION_E) .or. (unrot_dir == OBC_DIRECTION_W)) then
      jsd = segment%HI%jsd ; jed = segment%HI%jed
      IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
      do k=1,GV%ke
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
      do J=JsdB,JedB ; do i=isd,ied
        segment%normal_vel_bt(i,J) = 0.0
      enddo ; enddo
    endif
  enddo

end subroutine supercritical_set_OBC_data

!> \namespace supercritical_initialization
!!
!! The module configures the model for the "supercritical" experiment.
!! https://marine.rutgers.edu/po/index.php?model=test-problems&title=supercritical
end module supercritical_initialization
