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
use MOM_unit_scaling,   only : unit_scale_type
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_time_manager,   only : time_type, time_type_to_real

implicit none ; private

#include <MOM_memory.h>

public tidal_bay_set_OBC_data, tidal_bay_OBC_end
public register_tidal_bay_OBC

!> Control structure for tidal bay open boundaries.
type, public :: tidal_bay_OBC_CS ; private
  real :: tide_flow = 3.0e6         !< Maximum tidal flux [L2 Z T-1 ~> m3 s-1]
end type tidal_bay_OBC_CS

contains

!> Add tidal bay to OBC registry.
function register_tidal_bay_OBC(param_file, CS, US, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(tidal_bay_OBC_CS),   pointer    :: CS         !< tidal bay control structure.
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< OBC registry.
  logical                              :: register_tidal_bay_OBC
  character(len=32)  :: casename = "tidal bay"       !< This case's name.
  character(len=40)  :: mdl = "tidal_bay_initialization" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "register_tidal_bay_OBC called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  call get_param(param_file, mdl, "TIDAL_BAY_FLOW", CS%tide_flow, &
                 "Maximum total tidal volume flux.", &
                 units="m3 s-1", default=3.0d6, scale=US%m_s_to_L_T*US%m_to_L*US%m_to_Z)

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
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< layer thickness [H ~> m or kg m-2]
  type(time_type),         intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the tidal_bay example.
  real :: time_sec
  real :: cff_eta     ! The total column thickness anomalies associated with the inflow [H ~> m or kg m-2]
  real :: my_flux     ! The vlume flux through the face [L2 Z T-1 ~> m3 s-1]
  real :: total_area  ! The total face area of the OBCs [L Z ~> m2]
  real :: PI
  real :: flux_scale  ! A scaling factor for the areas [m2 H-1 L-1 ~> nondim or m3 kg-1]
  real, allocatable :: my_area(:,:) ! The total OBC inflow area [m2]
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

  flux_scale = GV%H_to_m*G%US%L_to_m

  time_sec = time_type_to_real(Time)
  cff_eta = 0.1*GV%m_to_H * sin(2.0*PI*time_sec/(12.0*3600.0))
  my_area=0.0
  my_flux=0.0
  segment => OBC%segment(1)

  do j=segment%HI%jsc,segment%HI%jec ; do I=segment%HI%IscB,segment%HI%IecB
    if (OBC%segnum_u(I,j) /= OBC_NONE) then
      do k=1,nz
        ! This area has to be in MKS units to work with reproducing_sum.
        my_area(1,j) = my_area(1,j) + h(I,j,k)*flux_scale*G%dyCu(I,j)
      enddo
    endif
  enddo ; enddo
  total_area = reproducing_sum(my_area)
  my_flux = - CS%tide_flow*SIN(2.0*PI*time_sec/(12.0*3600.0))

  do n = 1, OBC%number_of_segments
    segment => OBC%segment(n)

    if (.not. segment%on_pe) cycle

    segment%normal_vel_bt(:,:) = my_flux / (G%US%m_to_Z*G%US%m_to_L*total_area)
    segment%eta(:,:) = cff_eta

  enddo ! end segment loop

end subroutine tidal_bay_set_OBC_data

end module tidal_bay_initialization
