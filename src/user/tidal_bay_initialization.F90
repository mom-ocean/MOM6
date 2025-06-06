!> Configures the model for the "tidal_bay" experiment.
!! tidal_bay = Tidally resonant bay from Zygmunt Kowalik's class on tides.
module tidal_bay_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms,           only : reproducing_sum
use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type
use MOM_open_boundary,  only : OBC_segment_type, register_OBC
use MOM_open_boundary,  only : OBC_registry_type
use MOM_unit_scaling,   only : unit_scale_type
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_time_manager,   only : time_type, time_type_to_real

implicit none ; private

#include <MOM_memory.h>

public tidal_bay_set_OBC_data
public register_tidal_bay_OBC

!> Control structure for tidal bay open boundaries.
type, public :: tidal_bay_OBC_CS ; private
  real :: tide_flow = 3.0e6  !< Maximum tidal flux with the tidal bay configuration [L2 Z T-1 ~> m3 s-1]
  real :: tide_period        !< The period associated with the tidal bay configuration [T ~> s]
  real :: tide_ssh_amp       !< The magnitude of the sea surface height anomalies at the inflow
                             !! with the tidal bay configuration [Z ~> m]
end type tidal_bay_OBC_CS

contains

!> Add tidal bay to OBC registry.
function register_tidal_bay_OBC(param_file, CS, US, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(tidal_bay_OBC_CS),   intent(inout) :: CS      !< tidal bay control structure.
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< OBC registry.
  logical                              :: register_tidal_bay_OBC
  character(len=32)  :: casename = "tidal bay"       !< This case's name.
  character(len=40)  :: mdl = "tidal_bay_initialization" ! This module's name.

  call get_param(param_file, mdl, "TIDAL_BAY_FLOW", CS%tide_flow, &
                 "Maximum total tidal volume flux.", &
                 units="m3 s-1", default=3.0e6, scale=US%m_s_to_L_T*US%m_to_L*US%m_to_Z)
  call get_param(param_file, mdl, "TIDAL_BAY_PERIOD", CS%tide_period, &
                 "Period of the inflow in the tidal bay configuration.", &
                 units="s", default=12.0*3600.0, scale=US%s_to_T)
  call get_param(param_file, mdl, "TIDAL_BAY_SSH_ANOM", CS%tide_ssh_amp, &
                 "Magnitude of the sea surface height anomalies at the inflow with the "//&
                 "tidal bay configuration.", &
                 units="m", default=0.1, scale=US%m_to_Z)

  ! Register the open boundaries.
  call register_OBC(casename, param_file, OBC_Reg)
  register_tidal_bay_OBC = .true.

end function register_tidal_bay_OBC

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine tidal_bay_set_OBC_data(OBC, CS, G, GV, US, h, Time)
  type(ocean_OBC_type),    pointer    :: OBC  !< This open boundary condition type specifies
                                              !! whether, where, and what open boundary
                                              !! conditions are used.
  type(tidal_bay_OBC_CS),  intent(in) :: CS   !< tidal bay control structure.
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< layer thickness [H ~> m or kg m-2]
  type(time_type),         intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the tidal_bay example.
  real :: time_sec    ! Elapsed model time [T ~> s]
  real :: cff_eta     ! The sea surface height anomalies associated with the inflow [Z ~> m]
  real :: my_flux     ! The volume flux through the face [L2 Z T-1 ~> m3 s-1]
  real :: total_area  ! The total face area of the OBCs [L Z ~> m2]
  real :: normal_vel  ! The normal velocity through the inflow face [L T-1 ~> m s-1]
  real :: PI          ! The ratio of the circumference of a circle to its diameter [nondim]
  real, allocatable :: my_area(:,:) ! The total OBC inflow area [L Z ~> m2]
  integer :: turns    ! Number of index quarter turns
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, n
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment => NULL()

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  PI = 4.0*atan(1.0)

  turns = modulo(G%HI%turns, 4)

  if (.not.associated(OBC)) return

  time_sec = US%s_to_T*time_type_to_real(Time)
  cff_eta = CS%tide_ssh_amp * sin(2.0*PI*time_sec / CS%tide_period)

  segment => OBC%segment(1)

  if (turns == 0) then
    allocate(my_area(1:1,js:je), source=0.0)
    do j=segment%HI%jsc,segment%HI%jec ; do I=segment%HI%IscB,segment%HI%IecB
      if (OBC%segnum_u(I,j) > 0) then ! (segment%direction == OBC_DIRECTION_E)
        do k=1,nz
          my_area(1,j) = my_area(1,j) + h(i,j,k)*(GV%H_to_m*US%m_to_Z)*G%dyCu(I,j)
        enddo
      endif
    enddo ; enddo
  elseif (turns == 1) then
    allocate(my_area(is:ie,1:1), source=0.0)
    do J=segment%HI%JscB,segment%HI%JecB ; do i=segment%HI%isc,segment%HI%iec
      if (OBC%segnum_v(i,J) > 0) then ! (segment%direction == OBC_DIRECTION_N)
        do k=1,nz
          my_area(i,1) = my_area(i,1) + h(i,j,k)*(GV%H_to_m*US%m_to_Z)*G%dxCv(i,J)
        enddo
      endif
    enddo ; enddo
  elseif (turns == 2) then
    allocate(my_area(1:1,js:je), source=0.0)
    do j=segment%HI%jsc,segment%HI%jec ; do I=segment%HI%IscB,segment%HI%IecB
      if (OBC%segnum_u(I,j) < 0) then ! (segment%direction == OBC_DIRECTION_W)
        do k=1,nz
          my_area(1,j) = my_area(1,j) + h(i+1,j,k)*(GV%H_to_m*US%m_to_Z)*G%dyCu(I,j)
        enddo
      endif
    enddo ; enddo
  elseif (turns == 3) then
    allocate(my_area(is:ie,1:1), source=0.0)
    do J=segment%HI%JscB,segment%HI%JecB ; do i=segment%HI%isc,segment%HI%iec
      if (OBC%segnum_v(i,J) < 0) then ! (segment%direction == OBC_DIRECTION_S)
        do k=1,nz
          my_area(i,1) = my_area(i,1) + h(i,j+1,k)*(GV%H_to_m*US%m_to_Z)*G%dxCv(i,J)
        enddo
      endif
    enddo ; enddo
  endif

  total_area = reproducing_sum(my_area, unscale=US%Z_to_m*US%L_to_m)
  my_flux = - CS%tide_flow * SIN(2.0*PI*time_sec / CS%tide_period)
  normal_vel = my_flux / total_area
  if ((turns==2) .or. (turns==3)) normal_vel = -1.0 * normal_vel

  do n = 1, OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) cycle

    segment%normal_vel_bt(:,:) = normal_vel
    segment%SSH(:,:) = cff_eta

  enddo ! end segment loop

end subroutine tidal_bay_set_OBC_data

end module tidal_bay_initialization
