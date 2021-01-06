!> Configures the model for the "tidal_bay" experiment.
!! tidal_bay = Tidally resonant bay from Zygmunt Kowalik's class on tides.
module tidal_bay_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms,           only : reproducing_sum
use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary,  only : OBC_segment_type, register_OBC
use MOM_open_boundary,  only : OBC_registry_type
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_time_manager,   only : time_type, time_type_to_real

implicit none ; private

#include <MOM_memory.h>

public tidal_bay_set_OBC_data, tidal_bay_OBC_end
public register_tidal_bay_OBC

!> Control structure for tidal bay open boundaries.
type, public :: tidal_bay_OBC_CS ; private
  real :: tide_flow = 3.0e6         !< Maximum tidal flux.
end type tidal_bay_OBC_CS

contains

!> Add tidal bay to OBC registry.
function register_tidal_bay_OBC(param_file, CS, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(tidal_bay_OBC_CS),   pointer    :: CS         !< tidal bay control structure.
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< OBC registry.
  logical                              :: register_tidal_bay_OBC
  character(len=32)  :: casename = "tidal bay"       !< This case's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "register_tidal_bay_OBC called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Register the open boundaries.
  call register_OBC(casename, param_file, OBC_Reg)
  register_tidal_bay_OBC = .true.

end function register_tidal_bay_OBC

!> Clean up the tidal bay OBC from registry.
subroutine tidal_bay_OBC_end(CS)
  type(tidal_bay_OBC_CS), pointer    :: CS         !< tidal bay control structure.

  if (associated(CS)) then
    deallocate(CS)
  endif
end subroutine tidal_bay_OBC_end

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine tidal_bay_set_OBC_data(OBC, CS, G, GV, h, Time)
  type(ocean_OBC_type),    pointer    :: OBC  !< This open boundary condition type specifies
                                              !! whether, where, and what open boundary
                                              !! conditions are used.
  type(tidal_bay_OBC_CS),  pointer    :: CS   !< tidal bay control structure.
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< layer thickness.
  type(time_type),         intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the tidal_bay example.
  real :: time_sec, cff
  real :: my_flux, total_area
  real :: PI
  real, allocatable :: my_area(:,:)
  character(len=40)  :: mdl = "tidal_bay_set_OBC_data" ! This subroutine's name.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz, n
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment => NULL()

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  PI = 4.0*atan(1.0)

  if (.not.associated(OBC)) return

  allocate(my_area(1:1,js:je))

  time_sec = time_type_to_real(Time)
  cff = 0.1*sin(2.0*PI*time_sec/(12.0*3600.0))
  my_area=0.0
  my_flux=0.0
  segment => OBC%segment(1)

  do j=segment%HI%jsc,segment%HI%jec ; do I=segment%HI%IscB,segment%HI%IecB
    if (OBC%segnum_u(I,j) /= OBC_NONE) then
      do k=1,nz
        my_area(1,j) = my_area(1,j) + h(I,j,k)*G%US%L_to_m*G%dyCu(I,j)
      enddo
    endif
  enddo ; enddo
  total_area = reproducing_sum(my_area)
  my_flux = - CS%tide_flow*SIN(2.0*PI*time_sec/(12.0*3600.0))

  do n = 1, OBC%number_of_segments
    segment => OBC%segment(n)

    if (.not. segment%on_pe) cycle

    segment%normal_vel_bt(:,:) = G%US%m_s_to_L_T*my_flux/total_area
    segment%eta(:,:) = cff

  enddo ! end segment loop

end subroutine tidal_bay_set_OBC_data

end module tidal_bay_initialization
