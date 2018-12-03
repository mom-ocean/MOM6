!> Initialization for the "Neverland" configuration
module Neverland_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
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

public Neverland_initialize_topography
public Neverland_initialize_thickness

contains

!> This subroutine sets up the Neverland test case topography.
subroutine Neverland_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m
  ! Local variables
  real :: PI                   ! 3.1415926... calculated as 4*atan(1)
  real :: D0                   ! A constant to make the maximum     !
                               ! basin depth MAXIMUM_DEPTH.         !
  real :: x, y
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "Neverland_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  real :: nl_roughness_amp
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call MOM_mesg("  Neverland_initialization.F90, Neverland_initialize_topography: setting topography", 5)

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "NL_ROUGHNESS_AMP", nl_roughness_amp, &
                 "Amplitude of wavy signal in bathymetry.", default=0.05)

  PI = 4.0*atan(1.0)

!  Calculate the depth of the bottom.
  do i=is,ie
  do j=js,je
    x=(G%geoLonT(i,j)-G%west_lon)/G%len_lon
    y=(G%geoLatT(i,j)-G%south_lat)/G%len_lat
!  This sets topography that has a reentrant channel to the south.

    D(i,j) = 1.0 - (1.2 * spike(x,0.2) + 1.2 * spike(x-1.0,0.2)) * spike(MIN(0.0,y-0.3),0.2) & !< South America
              -  1.2 * spike(x-0.5,0.2) * spike(MIN(0.0,y-0.55),0.2)       & !< Africa
              -  1.1 * spike(y-1,0.12) - 1.1 * spike(y,0.12)               & !< The great northern wall and Antarctica
              -  1.2 * (spike(x,0.12)  + spike(x-1,0.12)) * spike(MAX(0.0,y-0.06),0.12)    & !< Antarctic Peninsula
              -  0.1 * (cosbell(x,0.1) + cosbell(x-1,0.1))                 & !< Drake Passage ridge
              -  0.5 * cosbell(x-0.16,0.05) * (cosbell(y-0.18,0.13)**0.4)  & !< Scotia Arc East
              -  0.4 * (cosbell(x-0.09,0.08)**0.4) * cosbell(y-0.26,0.05)  & !< Scotia Arc North
              -  0.4 * (cosbell(x-0.08,0.08)**0.4) * cosbell(y-0.1,0.05)   & !< Scotia Arc South
              -  nl_roughness_amp * cos(14*PI*x) * sin(14*PI*y)            & !< roughness
              -  nl_roughness_amp * cos(20*PI*x) * cos(20*PI*y)              !< roughness

    if (D(i,j) < 0.0) D(i,j) = 0.0
    D(i,j) = D(i,j) * max_depth
  enddo
  enddo

end subroutine Neverland_initialize_topography
! -----------------------------------------------------------------------------

!> Returns the value of a cosine-bell function evaluated at x/L
real function cosbell(x, L)
  real , intent(in) :: x       !< non-dimensional position
  real , intent(in) :: L       !< non-dimensional width
  real              :: PI      !< 3.1415926... calculated as 4*atan(1)

  PI      = 4.0*atan(1.0)
  cosbell = 0.5 * (1 + cos(PI*MIN(ABS(x/L),1.0)))
end function cosbell

!> Returns the value of a sin-spike function evaluated at x/L
real function spike(x, L)

  real , intent(in) :: x       !< non-dimensional position
  real , intent(in) :: L       !< non-dimensional width
  real              :: PI      !< 3.1415926... calculated as 4*atan(1)

  PI    = 4.0*atan(1.0)
  spike = (1 - sin(PI*MIN(ABS(x/L),0.5)))
end function spike

!> This subroutine initializes layer thicknesses for the Neverland test case,
!! by finding the depths of interfaces in a specified latitude-dependent
!! temperature profile with an exponentially decaying thermocline on top of a
!! linear stratification.
subroutine Neverland_initialize_thickness(h, G, GV, param_file, eqn_of_state, P_ref)
  type(ocean_grid_type),   intent(in) :: G                    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV                   !< The ocean's vertical grid structure.
  real, intent(out), dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h !< The thickness that is being
                                                              !! initialized, in H.
  type(param_file_type),   intent(in) :: param_file           !< A structure indicating the open
                                                              !! file to parse for model
                                                              !! parameter values.
  type(EOS_type),          pointer    :: eqn_of_state         !< integer that selects the
                                                              !! equation of state.
  real,                    intent(in) :: P_Ref                !< The coordinate-density
                                                              !! reference pressure in Pa.
  ! Local variables
  real :: e0(SZK_(G)+1)     ! The resting interface heights, in depth units (Z),
                            ! usually negative because it is positive upward.
  real, dimension(SZK_(G)) :: h_profile ! Vector of initial thickness profile (Z)
  real :: e_interface ! Current interface position (m)
  character(len=40)  :: mdl = "Neverland_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, k1, is, ie, js, je, nz, itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("  Neverland_initialization.F90, Neverland_initialize_thickness: setting thickness", 5)
  call get_param(param_file, mdl, "INIT_THICKNESS_PROFILE", h_profile, &
                 "Profile of initial layer thicknesses.", units="m", scale=GV%m_to_Z, &
                 fail_if_missing=.true.)

  ! e0 is the notional position of interfaces
  e0(1) = 0. ! The surface
  do k=1,nz
    e0(k+1) = e0(k) - h_profile(k)
  enddo

  do j=js,je ; do i=is,ie
    e_interface = -G%bathyT(i,j)
    do k=nz,1,-1
      h(i,j,k) = max( GV%Angstrom_H, GV%Z_to_H * (e0(k) - e_interface) )
      e_interface = max( e0(k), e_interface - GV%H_to_Z * h(i,j,k) )
    enddo
  enddo ; enddo

end subroutine Neverland_initialize_thickness

end module Neverland_initialization
