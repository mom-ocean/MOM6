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

public supercritical_set_OBC_data

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
  character(len=40)  :: mod = "supercritical_set_OBC_data" ! This subroutine's name.
  real :: zonal_flow
  integer :: i, j, k

  if (.not.associated(OBC)) call MOM_error(FATAL, 'supercritical_initialization.F90: '// &
        'supercritical_set_OBC_data() was called but OBC type was not initialized!')

  call get_param(param_file, mod, "SUPERCRITICAL_ZONAL_FLOW", zonal_flow, &
                 default=8.57, do_not_log=.true.) ! Do not log avoids maintaining identical descriptions

  allocate(OBC%u(G%IsdB:G%IedB,G%jsd:G%jed,G%ke)) ; OBC%u(:,:,:) = 0.0
  allocate(OBC%uh(G%IsdB:G%IedB,G%jsd:G%jed,G%ke)) ; OBC%uh(:,:,:) = 0.0
  if (.not.associated(OBC%ubt_outer)) then
    allocate(OBC%ubt_outer(G%IsdB:G%IedB,G%jsd:G%jed)) ; OBC%ubt_outer(:,:) = 0.0
  endif

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
