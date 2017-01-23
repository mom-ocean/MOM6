module BFB_surface_forcing
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
!
!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  Rewritten by Robert Hallberg, June 2009                            *
!*                                                                     *
!*  This file contains subroutines for specifying surface buoyancy     *
!*  forcing for the buoyancy-forced basin (BFB) case.                  *
!*  BFB_buoyancy_forcing is used to restore the surface buoayncy to    *
!*  a linear meridional ramp of temperature. The extent of the ramp    *
!*  can be specified by LFR_SLAT (linear forcing ramp southern         *
!*  latitude) and LFR_NLAT. The temperatures at these edges of the     *
!*  ramp can be specified by SST_S and SST_N.                          *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**
use MOM_diag_mediator, only : post_data, query_averaging_enabled
use MOM_diag_mediator, only : register_diag_field, diag_ctrl
use MOM_domains, only : pass_var, pass_vector, AGRID
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, param_file_type, log_version
use MOM_forcing_type, only : forcing, allocate_forcing_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, read_data
use MOM_time_manager, only : time_type, operator(+), operator(/), get_time
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_variables, only : surface

implicit none ; private

public BFB_buoyancy_forcing, BFB_surface_forcing_init

type, public :: BFB_surface_forcing_CS ; private
  !   This control structure should be used to store any run-time variables
  ! associated with the user-specified forcing.  It can be readily modified
  ! for a specific case, and because it is private there will be no changes
  ! needed in other code (although they will have to be recompiled).
  !   The variables in the cannonical example are used for some common
  ! cases, but do not need to be used.

  logical :: use_temperature ! If true, temperature and salinity are used as
                             ! state variables.
  logical :: restorebuoy     ! If true, use restoring surface buoyancy forcing.
  real :: Rho0               !   The density used in the Boussinesq
                             ! approximation, in kg m-3.
  real :: G_Earth            !   The gravitational acceleration in m s-2.
  real :: Flux_const         !   The restoring rate at the surface, in m s-1.
  real :: gust_const         !   A constant unresolved background gustiness
                             ! that contributes to ustar, in Pa.
  real :: SST_s              ! SST at the southern edge of the linear
                             ! forcing ramp
  real :: SST_n              ! SST at the northern edge of the linear
                             ! forcing ramp
  real :: lfrslat              ! Southern latitude where the linear forcing ramp
                             ! begins
  real :: lfrnlat              ! Northern latitude where the linear forcing ramp
                             ! ends
  real :: drho_dt            ! Rate of change of density with temperature.
                             ! Note that temperature is being used as a dummy
                             ! variable here. All temperatures are converted
                             ! into density.

  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
end type BFB_surface_forcing_CS

contains

subroutine BFB_buoyancy_forcing(state, fluxes, day, dt, G, CS)
  type(surface),                 intent(inout) :: state
  type(forcing),                 intent(inout) :: fluxes
  type(time_type),               intent(in)    :: day
  real,                          intent(in)    :: dt
  type(ocean_grid_type),         intent(in)    :: G
  type(BFB_surface_forcing_CS), pointer       :: CS

!    This subroutine specifies the current surface fluxes of buoyancy or
!  temperature and fresh water.  It may also be modified to add
!  surface fluxes of user provided tracers.

!    When temperature is used, there are long list of fluxes that need to be
!  set - essentially the same as for a full coupled model, but most of these
!  can be simply set to zero.  The net fresh water flux should probably be
!  set in fluxes%evap and fluxes%lprec, with any salinity restoring
!  appearing in fluxes%vprec, and the other water flux components
!  (fprec, lrunoff and frunoff) left as arrays full of zeros.
!  Evap is usually negative and precip is usually positive.  All heat fluxes
!  are in W m-2 and positive for heat going into the ocean.  All fresh water
!  fluxes are in kg m-2 s-1 and positive for water moving into the ocean.

! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (out)     fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      day_start - Start time of the fluxes.
!  (in)      day_interval - Length of time over which these fluxes
!                           will be applied.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure returned by a previous
!                 call to user_surface_forcing_init

  real :: Temp_restore   ! The temperature that is being restored toward, in C.
  real :: Salin_restore  ! The salinity that is being restored toward, in PSU.
  real :: density_restore  ! The potential density that is being restored
                         ! toward, in kg m-3.
  real :: rhoXcp ! The mean density times the heat capacity, in J m-3 K-1.
  real :: buoy_rest_const  ! A constant relating density anomalies to the
                           ! restoring buoyancy flux, in m5 s-3 kg-1.
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  !   When modifying the code, comment out this error message.  It is here
  ! so that the original (unmodified) version is not accidentally used.
  ! call MOM_error(FATAL, "User_buoyancy_surface_forcing: " // &
  !   "User forcing routine called without modification." )

  ! Allocate and zero out the forcing arrays, as necessary.  This portion is
  ! usually not changed.
  if (CS%use_temperature) then
    call alloc_if_needed(fluxes%evap, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%lprec, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%fprec, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%lrunoff, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%frunoff, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%vprec, isd, ied, jsd, jed)

    call alloc_if_needed(fluxes%sw, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%lw, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%latent, isd, ied, jsd, jed)
    call alloc_if_needed(fluxes%sens, isd, ied, jsd, jed)
  else ! This is the buoyancy only mode.
    call alloc_if_needed(fluxes%buoy, isd, ied, jsd, jed)
  endif


  ! MODIFY THE CODE IN THE FOLLOWING LOOPS TO SET THE BUOYANCY FORCING TERMS.

  if ( CS%use_temperature ) then
    ! Set whichever fluxes are to be used here.  Any fluxes that
    ! are always zero do not need to be changed here.
    do j=js,je ; do i=is,ie
      ! Fluxes of fresh water through the surface are in units of kg m-2 s-1
      ! and are positive downward - i.e. evaporation should be negative.
      fluxes%evap(i,j) = -0.0 * G%mask2dT(i,j)
      fluxes%lprec(i,j) = 0.0 * G%mask2dT(i,j)

      ! vprec will be set later, if it is needed for salinity restoring.
      fluxes%vprec(i,j) = 0.0

      !   Heat fluxes are in units of W m-2 and are positive into the ocean.
      fluxes%lw(i,j) = 0.0 * G%mask2dT(i,j)
      fluxes%latent(i,j) = 0.0 * G%mask2dT(i,j)
      fluxes%sens(i,j) = 0.0 * G%mask2dT(i,j)
      fluxes%sw(i,j) = 0.0 * G%mask2dT(i,j)
    enddo ; enddo
  else ! This is the buoyancy only mode.
    do j=js,je ; do i=is,ie
      !   fluxes%buoy is the buoyancy flux into the ocean in m2 s-3.  A positive
      ! buoyancy flux is of the same sign as heating the ocean.
      fluxes%buoy(i,j) = 0.0 * G%mask2dT(i,j)
    enddo ; enddo
  endif

  if (CS%restorebuoy) then
    if (CS%use_temperature) then
      call alloc_if_needed(fluxes%heat_added, isd, ied, jsd, jed)
      !   When modifying the code, comment out this error message.  It is here
      ! so that the original (unmodified) version is not accidentally used.
      call MOM_error(FATAL, "User_buoyancy_surface_forcing: " // &
        "Temperature and salinity restoring used without modification." )

      rhoXcp = CS%Rho0 * fluxes%C_p
      do j=js,je ; do i=is,ie
        !   Set Temp_restore and Salin_restore to the temperature (in C) and
        ! salinity (in PSU) that are being restored toward.
        Temp_restore = 0.0
        Salin_restore = 0.0

        fluxes%heat_added(i,j) = (G%mask2dT(i,j) * (rhoXcp * CS%Flux_const)) * &
            (Temp_restore - state%SST(i,j))
        fluxes%vprec(i,j) = - (G%mask2dT(i,j) * (CS%Rho0*CS%Flux_const)) * &
            ((Salin_restore - state%SSS(i,j)) / &
             (0.5 * (Salin_restore + state%SSS(i,j))))
      enddo ; enddo
    else
      !   When modifying the code, comment out this error message.  It is here
      ! so that the original (unmodified) version is not accidentally used.
      ! call MOM_error(FATAL, "User_buoyancy_surface_forcing: " // &
      !   "Buoyancy restoring used without modification." )

      ! The -1 is because density has the opposite sign to buoyancy.
      buoy_rest_const = -1.0 * (CS%G_Earth * CS%Flux_const) / CS%Rho0
      Temp_restore = 0.0
      do j=js,je ; do i=is,ie
       !   Set density_restore to an expression for the surface potential
       ! density in kg m-3 that is being restored toward.
        if (G%geoLatT(i,j) < CS%lfrslat) then
            Temp_restore = CS%SST_s
        else if (G%geoLatT(i,j) > CS%lfrnlat) then
            Temp_restore = CS%SST_n
        else
            Temp_restore = (CS%SST_s - CS%SST_n)/(CS%lfrslat - CS%lfrnlat) * &
                    (G%geoLatT(i,j) - CS%lfrslat) + CS%SST_s
        end if

        density_restore = Temp_restore*CS%drho_dt + CS%Rho0

        fluxes%buoy(i,j) = G%mask2dT(i,j) * buoy_rest_const * &
                          (density_restore - state%sfc_density(i,j))
      enddo ; enddo
    endif
  endif                                             ! end RESTOREBUOY

end subroutine BFB_buoyancy_forcing

subroutine alloc_if_needed(ptr, isd, ied, jsd, jed)
  ! If ptr is not associated, this routine allocates it with the given size
  ! and zeros out its contents.  This is equivalent to safe_alloc_ptr in
  ! MOM_diag_mediator, but is here so as to be completely transparent.
  real, pointer :: ptr(:,:)
  integer :: isd, ied, jsd, jed
  if (.not.ASSOCIATED(ptr)) then
    allocate(ptr(isd:ied,jsd:jed))
    ptr(:,:) = 0.0
  endif
end subroutine alloc_if_needed

subroutine BFB_surface_forcing_init(Time, G, param_file, diag, CS)
  type(time_type),                            intent(in) :: Time
  type(ocean_grid_type),                      intent(in) :: G
  type(param_file_type),                      intent(in) :: param_file
  type(diag_ctrl), target,                    intent(in) :: diag
  type(BFB_surface_forcing_CS), pointer    :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "BFB_surface_forcing" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "BFB_surface_forcing_init called with an associated "// &
                             "control structure.")
    return
  endif
  allocate(CS)
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state \n"//&
                 "variables.", default=.true.)

  call get_param(param_file, mod, "G_EARTH", CS%G_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)
  call get_param(param_file, mod, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  call get_param(param_file, mod, "LFR_SLAT", CS%lfrslat, &
                 "Southern latitude where the linear forcing ramp begins.", &
                 units="degrees", default = 20.0)
  call get_param(param_file, mod, "LFR_NLAT", CS%lfrnlat, &
                 "Northern latitude where the linear forcing ramp ends.", &
                 units="degrees", default = 40.0)
  call get_param(param_file, mod, "SST_S", CS%SST_s, &
                 "SST at the southern edge of the linear forcing ramp.", &
                 units="C", default = 20.0)
  call get_param(param_file, mod, "SST_N", CS%SST_n, &
                 "SST at the northern edge of the linear forcing ramp.", &
                 units="C", default = 10.0)
  call get_param(param_file, mod, "DRHO_DT", CS%drho_dt, &
                 "The rate of change of density with temperature.", &
                 units="kg m-3 K-1", default = -0.2)
  call get_param(param_file, mod, "GUST_CONST", CS%gust_const, &
                 "The background gustiness in the winds.", units="Pa", &
                 default=0.02)

  call get_param(param_file, mod, "RESTOREBUOY", CS%restorebuoy, &
                 "If true, the buoyancy fluxes drive the model back \n"//&
                 "toward some specified surface state with a rate \n"//&
                 "given by FLUXCONST.", default= .false.)
  if (CS%restorebuoy) then
    call get_param(param_file, mod, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes \n"//&
                 "to the relative surface anomalies (akin to a piston \n"//&
                 "velocity).  Note the non-MKS units.", units="m day-1", &
                 fail_if_missing=.true.)
    ! Convert CS%Flux_const from m day-1 to m s-1.
    CS%Flux_const = CS%Flux_const / 86400.0
  endif

end subroutine BFB_surface_forcing_init

end module BFB_surface_forcing
