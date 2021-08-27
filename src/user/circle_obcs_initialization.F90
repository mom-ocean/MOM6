!> Configures the model for the "circle_obcs" experiment which tests
!! Open Boundary Conditions radiating an SSH anomaly.
module circle_obcs_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

implicit none ; private

#include <MOM_memory.h>

public circle_obcs_initialize_thickness

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> This subroutine initializes layer thicknesses for the circle_obcs experiment.
subroutine circle_obcs_initialize_thickness(h, depth_tot, G, GV, param_file, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G   !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV  !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.

  real :: e0(SZK_(GV)+1)   ! The resting interface heights, in depth units [Z ~> m], usually
                           ! negative because it is positive upward.
  real :: eta1D(SZK_(GV)+1)! Interface height relative to the sea surface
                           ! positive upward, in depth units [Z ~> m].
  real :: IC_amp           ! The amplitude of the initial height displacement [H ~> m or kg m-2].
  real :: diskrad, rad, xCenter, xRadius, lonC, latC, xOffset
  logical :: just_read
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "circle_obcs_initialization"   ! This module's name.
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) &
    call MOM_mesg("  circle_obcs_initialization.F90, circle_obcs_initialize_thickness: setting thickness", 5)

  if (.not.just_read) call log_version(param_file, mdl, version, "")
  ! Parameters read by cartesian grid initialization
  call get_param(param_file, mdl, "DISK_RADIUS", diskrad, &
                 "The radius of the initially elevated disk in the "//&
                 "circle_obcs test case.", units=G%x_axis_units, &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "DISK_X_OFFSET", xOffset, &
                 "The x-offset of the initially elevated disk in the "//&
                 "circle_obcs test case.", units=G%x_axis_units, &
                 default = 0.0, do_not_log=just_read)
  call get_param(param_file, mdl, "DISK_IC_AMPLITUDE", IC_amp, &
                 "Initial amplitude of interface height displacements "//&
                 "in the circle_obcs test case.", &
                 units='m', default=5.0, scale=GV%m_to_H, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  do k=1,nz
    e0(K) = -G%max_depth * real(k-1) / real(nz)
  enddo

  ! Uniform thicknesses for base state
  do j=js,je ; do i=is,ie
    eta1D(nz+1) = -depth_tot(i,j)
    do k=nz,1,-1
      eta1D(K) = e0(K)
      if (eta1D(K) < (eta1D(K+1) + GV%Angstrom_Z)) then
        eta1D(K) = eta1D(K+1) + GV%Angstrom_Z
        h(i,j,k) = GV%Angstrom_H
      else
        h(i,j,k) = GV%Z_to_H * (eta1D(K) - eta1D(K+1))
      endif
    enddo
  enddo ; enddo

  ! Perturb base state by circular anomaly in center
  k=nz
  latC = G%south_lat + 0.5*G%len_lat
  lonC = G%west_lon + 0.5*G%len_lon + xOffset
  do j=js,je ; do i=is,ie
    rad = sqrt((G%geoLonT(i,j)-lonC)**2+(G%geoLatT(i,j)-latC)**2)/(diskrad)
    ! if (rad <= 6.*diskrad) h(i,j,k) = h(i,j,k)+10.0*exp( -0.5*( rad**2 ) )
    rad = min( rad, 1. ) ! Flatten outside radius of diskrad
    rad = rad*(2.*asin(1.)) ! Map 0-1 to 0-pi
    if (nz==1) then
      ! The model is barotropic
      h(i,j,k) = h(i,j,k) + IC_amp * 0.5*(1.+cos(rad)) ! cosine bell
    else
      ! The model is baroclinic
      do k = 1, nz
        h(i,j,k) = h(i,j,k) - 0.5*(1.+cos(rad)) * IC_amp * real( 2*k-nz )
      enddo
    endif
  enddo ; enddo

end subroutine circle_obcs_initialize_thickness

end module circle_obcs_initialization
