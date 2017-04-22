module shelfwave_initialization
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

use MOM_domains,        only : sum_across_PEs
use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary,  only : OBC_segment_type, register_OBC
use MOM_open_boundary,  only : OBC_registry_type
use MOM_time_manager,   only : time_type, set_time, time_type_to_real

implicit none ; private

#include <MOM_memory.h>

character(len=40) :: mod = "shelfwave_initialization" ! This module's name.

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public shelfwave_initialize_topography
public shelfwave_set_OBC_data
public register_shelfwave_OBC, shelfwave_OBC_end

!> Control structure for shelfwave open boundaries.
type, public :: shelfwave_OBC_CS ; private
  real :: Lx = 20e3         !< Long-shore length scale.
  real :: Ly = 10e3         !< Cross-shore length scale.
  real :: f0 = 1.e-4        !< Coriolis parameter.
end type shelfwave_OBC_CS

contains

!> Add shelfwave to OBC registry.
function register_shelfwave_OBC(param_file, CS, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(shelfwave_OBC_CS),   pointer    :: CS         !< shelfwave control structure.
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< OBC registry.
  logical                              :: register_shelfwave_OBC
  character(len=32)  :: casename = "shelfwave"       !< This case's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "register_shelfwave_OBC called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Register the tracer for horizontal advection & diffusion.
  call register_OBC(casename, param_file, OBC_Reg)
  call get_param(param_file,mod,"F_0",CS%f0, &
                 do_not_log=.true.)
  call get_param(param_file,mod,"SHELFWAVE_X_LENGTH_SCALE",CS%Lx, &
                 "Length scale of shelfwave in x-direction.\n"//&
                 "Set to zero make wave uniform in the x-direction.", &
                 units="Same as x,y", default=20.)
  call get_param(param_file,mod,"SHELFWAVE_Y_LENGTH_SCALE",CS%Ly, &
                 "Length scale of shelfwave bathymetry in y-direction.\n"//&
                 "Set to zero make topography uniform in the y-direction.", &
                 units="Same as x,y", default=10.)
  register_shelfwave_OBC = .true.

end function register_shelfwave_OBC

!> Clean up the shelfwave OBC from registry.
subroutine shelfwave_OBC_end(CS)
  type(shelfwave_OBC_CS), pointer    :: CS         !< shelfwave control structure.

  if (associated(CS)) then
    deallocate(CS)
  endif
end subroutine shelfwave_OBC_end

!> Initialization of topography.
subroutine shelfwave_initialize_topography ( D, G, param_file, max_depth )
  ! Arguments
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m

  ! Local variables
  integer   :: i, j
  real      :: x, y, delta, rLy, Ly, H0

  call get_param(param_file,mod,"SHELFWAVE_Y_LENGTH_SCALE",Ly, &
                 default=10., do_not_log=.true.)
  call get_param(param_file,mod,"MINIMUM_DEPTH",H0, &
                 default=10., do_not_log=.true.)

  Ly = Ly / G%len_lat
  rLy = 0. ; if (Ly>0.) rLy = 1. / Ly
  do i=G%isc,G%iec
    do j=G%jsc,G%jec
      ! Compute normalized zonal coordinates (x,y=0 at center of domain)
      x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon - 0.5
      y = ( G%geoLatT(i,j) - G%south_lat ) / G%len_lat - 0.5
      D(i,j) = H0 * exp(2 * rLy * y)
    enddo
  enddo

end subroutine shelfwave_initialize_topography

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine shelfwave_set_OBC_data(OBC, CS, G, h, Time)
  type(ocean_OBC_type),   pointer    :: OBC  !< This open boundary condition type specifies
                                             !! whether, where, and what open boundary
                                             !! conditions are used.
  type(shelfwave_OBC_CS), pointer    :: CS   !< tidal bay control structure.
  type(ocean_grid_type),  intent(in) :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h !< layer thickness.
  type(time_type),        intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the tidal_bay example.
  real :: my_amp, time_sec
  real :: PI, cos_wt, cos_ly, sin_wt, sin_ly, omega
  real :: x, y, delta, Lx, rLx, Ly, rLy, f0
  character(len=40)  :: mod = "shelfwave_set_OBC_data" ! This subroutine's name.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, jj, n
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) return

  PI = 4.0*atan(1.0) ;

  Lx = CS%Lx / G%len_lon
  Ly = CS%Ly / G%len_lat
  rLx = 0. ; if (Lx>0.) rLx = 1. / Lx
  rLy = 0. ; if (Ly>0.) rLy = 1. / Ly

  time_sec = time_type_to_real(Time)
  omega = 2.0 * rLy * Lx * f0
  my_amp = 1.0
  jj = 1
  do n = 1, OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) cycle
    if (segment%Is_obc /= is - 1) cycle

    IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
    jsd = segment%HI%jsd ; jed = segment%HI%jed
    do j=jsd,jed ; do I=IsdB,IedB
      x = ( 0.5 * (G%geoLonT(i,j) + G%geoLonT(i+1,j)) - G%west_lon ) / G%len_lon
      y = ( 0.5 * (G%geoLatT(i,j) + G%geoLatT(i+1,j)) - G%south_lat ) / G%len_lat
      sin_wt = sin(rLx*x - omega*time_sec)
      cos_wt = cos(rLx*x - omega*time_sec)
      sin_ly = sin(jj * PI * y)
      cos_ly = cos(jj * PI * y)
      segment%normal_vel_bt(I,j) = my_amp * exp(- rLy * y) * cos_wt * &
           (rLy * sin_ly + jj * PI * rLy * cos_ly)
!     segment%tangential_vel_bt(I,j) = my_amp * rLx * exp(- rLy * y) * sin_wt * sin_ly
!     segment%vorticity_bt(I,j) = my_amp * exp(- rLy * y) * cos_wt * &
!           (rLy**2 + rLx**2 + (jj*PI)**2)
    enddo ; enddo
  enddo

end subroutine shelfwave_set_OBC_data

!> \namespace shelfwave_initialization
!!
!! The module configures the model for the idealized shelfwave
!! test case.
end module shelfwave_initialization
