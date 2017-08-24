module external_gwave_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_variables, only : thermo_var_ptrs
implicit none ; private

#include <MOM_memory.h>

public external_gwave_initialize_thickness

contains

! -----------------------------------------------------------------------------
!> This subroutine initializes layer thicknesses for the external_gwave experiment.
subroutine external_gwave_initialize_thickness(h, G, param_file, just_read_params)
  type(ocean_grid_type), intent(in) :: G                      !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(out) :: h           !< The thickness that is being initialized, in m.
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.

  real :: e0(SZK_(G))     ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: e_pert(SZK_(G)) ! Interface height perturbations, positive     !
                          ! upward, in m.                                !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  real :: ssh_anomaly_height ! Vertical height of ssh anomaly
  real :: ssh_anomaly_width ! Lateral width of anomaly
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40)  :: mdl = "external_gwave_initialize_thickness" ! This subroutine's name.
! This include declares and sets the variable "version".
#include "version_variable.h"
  integer :: i, j, k, is, ie, js, je, nz
  real :: PI, Xnondim

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) &
    call MOM_mesg("  external_gwave_initialization.F90, external_gwave_initialize_thickness: setting thickness", 5)

  if (.not.just_read) call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "SSH_ANOMALY_HEIGHT", ssh_anomaly_height, &
                 "The vertical displacement of the SSH anomaly. ", units="m", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "SSH_ANOMALY_WIDTH", ssh_anomaly_width, &
                 "The lateral width of the SSH anomaly. ", units="coordinate", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  PI = 4.0*atan(1.0)
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    Xnondim = (G%geoLonT(i,j)-G%west_lon-0.5*G%len_lon) / ssh_anomaly_width
    Xnondim = min(1., abs(Xnondim))
    eta1D(1) = ssh_anomaly_height * 0.5 * ( 1. + cos(PI*Xnondim) ) ! Cosine bell
    do k=2,nz
      eta1D(K) = -G%max_depth & ! Stretch interior interfaces with SSH
              + (eta1D(1)+G%max_depth) * ( real(nz+1-k)/real(nz) ) ! Stratification
    enddo
    eta1D(nz+1) = -G%max_depth ! Force bottom interface to bottom
    do k=1,nz
      h(i,j,k) = eta1D(K) - eta1D(K+1)
    enddo
  enddo ; enddo

end subroutine external_gwave_initialize_thickness
! -----------------------------------------------------------------------------

!> \namespace external_gwave_initialization
!!
!! The module configures the model for the "external_gwave" experiment.
!! external_gwave = External Gravity Wave
end module external_gwave_initialization
