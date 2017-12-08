module lock_exchange_initialization

! This file is part of MOM6. See LICENSE.md for the license.

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
subroutine lock_exchange_initialize_thickness(h, G, GV, param_file, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized, in m.
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.

  real :: e0(SZK_(GV))     ! The resting interface heights, in m, usually !
                           ! negative because it is positive upward.      !
  real :: e_pert(SZK_(GV)) ! Interface height perturbations, positive     !
                           ! upward, in m.                                !
  real :: eta1D(SZK_(GV)+1)! Interface height relative to the sea surface !
                           ! positive upward, in m.                       !
  real :: front_displacement ! Vertical displacement acrodd front
  real :: thermocline_thickness ! Thickness of stratified region
  logical :: just_read    ! If true, just read parameters but set nothing.
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "lock_exchange_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) &
    call MOM_mesg("  lock_exchange_initialization.F90, lock_exchange_initialize_thickness: setting thickness", 5)

  if (.not.just_read) call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "FRONT_DISPLACEMENT", front_displacement, &
                 "The vertical displacement of interfaces across the front. \n"//&
                 "A value larger in magnitude that MAX_DEPTH is truncated,", &
                 units="m", fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "THERMOCLINE_THICKNESS", thermocline_thickness, &
                 "The thickness of the thermocline in the lock exchange \n"//&
                 "experiment.  A value of zero creates a two layer system \n"//&
                 "with vanished layers in between the two inflated layers.", &
                 default=0., units="m", do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

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
      eta1D(K) = max( eta1D(K), eta1D(K+1) + GV%Angstrom_Z )
    enddo
    eta1D(1) = 0. ! Force bottom interface to bottom
    do k=2,nz ! Make sure interfaces decrease downwards
      eta1D(K) = min( eta1D(K), eta1D(K-1) - GV%Angstrom_Z )
    enddo
    do k=nz,1,-1
      h(i,j,k) = eta1D(K) - eta1D(K+1)
    enddo
  enddo ; enddo

end subroutine lock_exchange_initialize_thickness
! -----------------------------------------------------------------------------

!> \namespace lock_exchange_initialization
!!
!! The module configures the model for the "lock_exchange" experiment.
!! lock_exchange = A 2-d density driven hydraulic exchange flow.
end module lock_exchange_initialization
