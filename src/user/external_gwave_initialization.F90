!> Initialization for the "external gravity wave wave" configuration
module external_gwave_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
implicit none ; private

#include <MOM_memory.h>

public external_gwave_initialize_thickness

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> This subroutine initializes layer thicknesses for the external_gwave experiment.
subroutine external_gwave_initialize_thickness(h, G, GV, US, param_file, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [H ~> m or kg m-2].
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
  real :: eta1D(SZK_(GV)+1)  ! Interface height relative to the sea surface
                             ! positive upward [Z ~> m].
  real :: ssh_anomaly_height ! Vertical height of ssh anomaly [Z ~> m]
  real :: ssh_anomaly_width ! Lateral width of anomaly [degrees]
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40)  :: mdl = "external_gwave_initialize_thickness" ! This subroutine's name.
! This include declares and sets the variable "version".
#include "version_variable.h"
  integer :: i, j, k, is, ie, js, je, nz
  real :: PI, Xnondim

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) &
    call MOM_mesg("  external_gwave_initialization.F90, external_gwave_initialize_thickness: setting thickness", 5)

  if (.not.just_read) call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "SSH_ANOMALY_HEIGHT", ssh_anomaly_height, &
                 "The vertical displacement of the SSH anomaly. ", units="m", &
                 fail_if_missing=.not.just_read, do_not_log=just_read, scale=US%m_to_Z)
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
      h(i,j,k) = GV%Z_to_H * (eta1D(K) - eta1D(K+1))
    enddo
  enddo ; enddo

end subroutine external_gwave_initialize_thickness

end module external_gwave_initialization
