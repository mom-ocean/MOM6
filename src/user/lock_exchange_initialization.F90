module lock_exchange_initialization
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

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public lock_exchange_initialize_thickness

contains

!> This subroutine initializes layer thicknesses for the lock_exchange experiment.
! -----------------------------------------------------------------------------
subroutine lock_exchange_initialize_thickness(h, G, GV, param_file)
  type(ocean_grid_type),   intent(in) :: G                    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV                   !< The ocean's vertical grid structure.
  real, intent(out), dimension(SZI_(G),SZJ_(G), SZK_(G)) :: h !< The thickness that is being initialized.
  type(param_file_type),   intent(in) :: param_file           !< A structure indicating the open file to
                                                              !! parse for model parameter values.

  real :: e0(SZK_(G))     ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: e_pert(SZK_(G)) ! Interface height perturbations, positive     !
                          ! upward, in m.                                !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  real :: front_displacement ! Vertical displacement acrodd front
  real :: thermocline_thickness ! Thickness of stratified region
  character(len=40)  :: mod = "lock_exchange_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("  lock_exchange_initialization.F90, lock_exchange_initialize_thickness: setting thickness", 5)

  call get_param(param_file, mod, "FRONT_DISPLACEMENT", front_displacement, &
                 "The vertical displacement of interfaces across the front. \n"//&
                 "A value larger in magnitude that MAX_DEPTH is truncated,", &
                 units="m", fail_if_missing=.true.)
  call get_param(param_file, mod, "THERMOCLINE_THICKNESS", thermocline_thickness, &
                 "The thickness of the thermocline in the lock exchange \n"//&
                 "experiment.  A value of zero creates a two layer system \n"//&
                 "with vanished layers in between the two inflated layers.", &
                 default=0., units="m")

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    do k=2,nz
      eta1D(K) = -0.5 * G%max_depth & ! Middle of column
              - thermocline_thickness * ( (real(k-1))/real(nz) -0.5 ) ! Stratification
      if (G%geoLonT(i,j)-G%west_lon < 0.5 * G%len_lon) then
        eta1D(K)=eta1D(K) + 0.5 * front_displacement
      elseif (G%geoLonT(i,j)-G%west_lon > 0.5 * G%len_lon) then
        eta1D(K)=eta1D(K) - 0.5 * front_displacement
      endif
    enddo
    eta1D(nz+1) = -G%max_depth ! Force bottom interface to bottom
    do k=nz,2,-1 ! Make sure interfaces increase upwards
      eta1D(K) = max( eta1D(K), eta1D(K+1) + GV%Angstrom )
    enddo
    eta1D(1) = 0. ! Force bottom interface to bottom
    do k=2,nz ! Make sure interfaces decrease downwards
      eta1D(K) = min( eta1D(K), eta1D(K-1) - GV%Angstrom )
    enddo
    do k=nz,1,-1
      h(i,j,k) = eta1D(K) - eta1D(K+1)
    enddo
  enddo ; enddo

end subroutine lock_exchange_initialize_thickness
! -----------------------------------------------------------------------------

!> \class lock_exchange_initialization
!!
!! The module configures the model for the "lock_exchange" experiment.
!! lock_exchange = A 2-d density driven hydraulic exchange flow.
end module lock_exchange_initialization
