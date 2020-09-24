!> Initialization of the boundary-forced-basing configuration
module BFB_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use MOM_verticalGrid, only : verticalGrid_type
implicit none ; private

#include <MOM_memory.h>

public BFB_set_coord
public BFB_initialize_sponges_southonly

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Unsafe model variable
!! \todo Remove this module variable
logical :: first_call = .true.

contains

!> This subroutine specifies the vertical coordinate in terms of temperature at the surface and at the bottom.
!! This case is set up in such a way that the temperature of the topmost layer is equal to the SST at the
!! southern edge of the domain. The temperatures are then converted to densities of the top and bottom layers
!! and linearly interpolated for the intermediate layers.
subroutine BFB_set_coord(Rlay, g_prime, GV, US, param_file, eqn_of_state)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters
  type(EOS_type),           pointer     :: eqn_of_state !< Equation of state structure
  ! Local variables
  real                                 :: drho_dt, SST_s, T_bot, rho_top, rho_bot
  integer                              :: k, nz
  character(len=40)  :: mdl = "BFB_set_coord" ! This subroutine's name.

  call get_param(param_file, mdl, "DRHO_DT", drho_dt, &
          "Rate of change of density with temperature.", &
           units="kg m-3 K-1", default=-0.2, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "SST_S", SST_s, &
          "SST at the suothern edge of the domain.", units="C", default=20.0)
  call get_param(param_file, mdl, "T_BOT", T_bot, &
                 "Bottom Temp", units="C", default=5.0)
  rho_top = GV%Rho0 + drho_dt*SST_s
  rho_bot = GV%Rho0 + drho_dt*T_bot
  nz = GV%ke

  do k = 1,nz
    Rlay(k) = (rho_bot - rho_top)/(nz-1)*real(k-1) + rho_top
    if (k >1) then
      g_prime(k) = (Rlay(k) - Rlay(k-1)) * GV%g_Earth / (GV%Rho0)
    else
      g_prime(k) = GV%g_Earth
    endif
    !Rlay(:) = 0.0
    !g_prime(:) = 0.0
  enddo

  if (first_call) call write_BFB_log(param_file)

end subroutine BFB_set_coord

!> This subroutine sets up the sponges for the southern bouundary of the domain. Maximum damping occurs
!! within 2 degrees lat of the boundary. The damping linearly decreases northward over the next 2 degrees.
subroutine BFB_initialize_sponges_southonly(G, GV, US, use_temperature, tv, param_file, CSp, h)
  type(ocean_grid_type),   intent(in) :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  logical,                 intent(in) :: use_temperature !< If true, temperature and salinity are used as
                                            !! state variables.
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(sponge_CS),         pointer    :: CSp  !< A pointer to the sponge control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]

  ! Local variables
  real :: eta(SZI_(G),SZJ_(G),SZK_(GV)+1) ! A temporary array for eta, in depth units [Z ~> m].
  real :: Idamp(SZI_(G),SZJ_(G))    ! The inverse damping rate [T-1 ~> s-1].
  real :: H0(SZK_(GV))              ! Resting layer thicknesses in depth units [Z ~> m].
  real :: min_depth                 ! The minimum ocean depth in depth units [Z ~> m].
  real :: slat, wlon, lenlat, lenlon, nlat
  real :: max_damping               ! The maximum damping rate [T-1 ~> s-1]
  character(len=40)  :: mdl = "BFB_initialize_sponges_southonly" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  eta(:,:,:) = 0.0 ; Idamp(:,:) = 0.0

!  Here the inverse damping time [T-1 ~> s-1], is set. Set Idamp to 0
!  wherever there is no sponge, and the subroutines that are called
!  will automatically set up the sponges only where Idamp is positive
!  and mask2dT is 1.

!   Set up sponges for DOME configuration
  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0, scale=US%m_to_Z)

  call get_param(param_file, mdl, "SOUTHLAT", slat, &
                 "The southern latitude of the domain.", units="degrees")
  call get_param(param_file, mdl, "LENLAT", lenlat, &
                 "The latitudinal length of the domain.", units="degrees")
  call get_param(param_file, mdl, "WESTLON", wlon, &
                 "The western longitude of the domain.", units="degrees", default=0.0)
  call get_param(param_file, mdl, "LENLON", lenlon, &
                 "The longitudinal length of the domain.", units="degrees")
  nlat = slat + lenlat
  do k=1,nz ; H0(k) = -G%max_depth * real(k-1) / real(nz) ; enddo

  ! Use for meridional thickness profile initialization
!  do k=1,nz ; H0(k) = -G%max_depth * real(k-1) / real(nz-1) ; enddo

  max_damping = 1.0  / (86400.0*US%s_to_T)

  do i=is,ie; do j=js,je
    if (G%bathyT(i,j) <= min_depth) then ; Idamp(i,j) = 0.0
    elseif (G%geoLatT(i,j) < slat+2.0) then ; Idamp(i,j) = max_damping
    elseif (G%geoLatT(i,j) < slat+4.0) then
      Idamp(i,j) = max_damping * (slat+4.0-G%geoLatT(i,j))/2.0
    else ; Idamp(i,j) = 0.0
    endif

    ! These will be streched inside of apply_sponge, so they can be in
    ! depth space for Boussinesq or non-Boussinesq models.

    ! This section is used for uniform thickness initialization
    do k = 1,nz; eta(i,j,k) = H0(k); enddo

    ! The below section is used for meridional temperature profile thickness initiation
    ! do k = 1,nz; eta(i,j,k) = H0(k); enddo
    ! if (G%geoLatT(i,j) > 40.0) then
    !   do k = 1,nz
    !     eta(i,j,k) = -G%Angstrom_Z*(k-1)
    !   enddo
    ! elseif (G%geoLatT(i,j) > 20.0) then
    !   do k = 1,nz
    !     eta(i,j,k) = min(H0(k) + (G%geoLatT(i,j) - 20.0)*(G%max_depth - nz*G%Angstrom_Z)/20.0, &
    !                      -(k-1)*G%Angstrom_Z)
    !   enddo
    ! endif
    eta(i,j,nz+1) = -G%max_depth

  enddo ; enddo

!  This call sets up the damping rates and interface heights.
!  This sets the inverse damping timescale fields in the sponges.    !
  call initialize_sponge(Idamp, eta, G, param_file, CSp, GV)

!   Now register all of the fields which are damped in the sponge.   !
! By default, momentum is advected vertically within the sponge, but !
! momentum is typically not damped within the sponge.                !

  if (first_call) call write_BFB_log(param_file)

end subroutine BFB_initialize_sponges_southonly

!> Write output about the parameter values being used.
subroutine write_BFB_log(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "BFB_initialization" ! This module's name.

  call log_version(param_file, mdl, version)
  first_call = .false.

end subroutine write_BFB_log

end module BFB_initialization
