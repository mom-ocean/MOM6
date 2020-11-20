!> Template for user to code up surface forcing.
module user_surface_forcing

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, query_averaging_enabled
use MOM_diag_mediator, only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_domains, only : pass_var, pass_vector, AGRID
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, param_file_type, log_version
use MOM_forcing_type, only : forcing, mech_forcing
use MOM_forcing_type, only : allocate_forcing_type, allocate_mech_forcing
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, read_data
use MOM_time_manager, only : time_type, operator(+), operator(/)
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface

implicit none ; private

public USER_wind_forcing, USER_buoyancy_forcing, USER_surface_forcing_init

!>   This control structure should be used to store any run-time variables
!! associated with the user-specified forcing.
!!
!! It can be readily modified for a specific case, and because it is private there
!! will be no changes needed in other code (although they will have to be recompiled).
type, public :: user_surface_forcing_CS ; private
  !   The variables in the cannonical example are used for some common
  ! cases, but do not need to be used.

  logical :: use_temperature !< If true, temperature and salinity are used as state variables.
  logical :: restorebuoy     !< If true, use restoring surface buoyancy forcing.
  real :: Rho0               !< The density used in the Boussinesq approximation [R ~> kg m-3].
  real :: G_Earth            !< The gravitational acceleration [L2 Z-1 s-2 ~> m s-2].
  real :: Flux_const         !< The restoring rate at the surface [Z T-1 ~> m s-1].
  real :: gust_const         !< A constant unresolved background gustiness
                             !! that contributes to ustar [R L Z T-1 ~> Pa].

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                             !! timing of diagnostic output.
end type user_surface_forcing_CS

contains

!> This subroutine sets the surface wind stresses, forces%taux and forces%tauy, in [R Z L T-2 ~> Pa].
!! These are the stresses in the direction of the model grid (i.e. the same
!! direction as the u- and v- velocities).
subroutine USER_wind_forcing(sfc_state, forces, day, G, US, CS)
  type(surface),                 intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),            intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),               intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),         intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type),         intent(in)    :: US   !< A dimensional unit scaling type
  type(user_surface_forcing_CS), pointer       :: CS   !< A pointer to the control structure returned
                                                       !! by a previous call to user_surface_forcing_init

  ! Local variables
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq

  !   When modifying the code, comment out this error message.  It is here
  ! so that the original (unmodified) version is not accidentally used.
  call MOM_error(FATAL, "User_wind_surface_forcing: " // &
     "User forcing routine called without modification." )

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  ! Allocate the forcing arrays, if necessary.
  call allocate_mech_forcing(G, forces, stress=.true., ustar=.true.)

  !  Set the surface wind stresses, in units of [R L Z T-1 ~> Pa].  A positive taux
  !  accelerates the ocean to the (pseudo-)east.

  !  The i-loop extends to is-1 so that taux can be used later in the
  ! calculation of ustar - otherwise the lower bound would be Isq.
  do j=js,je ; do I=is-1,Ieq
    ! Change this to the desired expression.
    forces%taux(I,j) = G%mask2dCu(I,j) * 0.0*US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z
  enddo ; enddo
  do J=js-1,Jeq ; do i=is,ie
    forces%tauy(i,J) = G%mask2dCv(i,J) * 0.0  ! Change this to the desired expression.
  enddo ; enddo

  !    Set the surface friction velocity, in units of [Z T-1 ~> m s-1].  ustar
  !  is always positive.
  if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
    !  This expression can be changed if desired, but need not be.
   forces%ustar(i,j) = G%mask2dT(i,j) * sqrt((CS%gust_const + &
            sqrt(0.5*(forces%taux(I-1,j)**2 + forces%taux(I,j)**2) + &
                 0.5*(forces%tauy(i,J-1)**2 + forces%tauy(i,J)**2))) * (US%L_to_Z/CS%Rho0))
  enddo ; enddo ; endif

end subroutine USER_wind_forcing

!>    This subroutine specifies the current surface fluxes of buoyancy or
!!  temperature and fresh water.  It may also be modified to add
!!  surface fluxes of user provided tracers.
subroutine USER_buoyancy_forcing(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),                 intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(forcing),                 intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),               intent(in)    :: day  !< The time of the fluxes
  real,                          intent(in)    :: dt   !< The amount of time over which
                                                       !! the fluxes apply [s]
  type(ocean_grid_type),         intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),         intent(in)    :: US   !< A dimensional unit scaling type
  type(user_surface_forcing_CS), pointer       :: CS   !< A pointer to the control structure returned
                                                       !! by a previous call to user_surface_forcing_init

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
!  fluxes are in [R Z T-1 ~> kg m-2 s-1] and positive for water moving into the ocean.

  ! Local variables
  real :: Temp_restore   ! The temperature that is being restored toward [degC].
  real :: Salin_restore  ! The salinity that is being restored toward [ppt]
  real :: density_restore  ! The potential density that is being restored
                         ! toward [R ~> kg m-3].
  real :: rhoXcp         ! The mean density times the heat capacity [Q R degC-1 ~> J m-3 degC-1].
  real :: buoy_rest_const  ! A constant relating density anomalies to the
                           ! restoring buoyancy flux [L2 m3 T-3 kg-1 ~> m5 s-3 kg-1].

  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  !   When modifying the code, comment out this error message.  It is here
  ! so that the original (unmodified) version is not accidentally used.
  call MOM_error(FATAL, "User_buoyancy_surface_forcing: " // &
    "User forcing routine called without modification." )

  ! Allocate and zero out the forcing arrays, as necessary.  This portion is
  ! usually not changed.
  if (CS%use_temperature) then
    call safe_alloc_ptr(fluxes%evap, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%lprec, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%fprec, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%lrunoff, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%frunoff, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%vprec, isd, ied, jsd, jed)

    call safe_alloc_ptr(fluxes%sw, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%lw, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%latent, isd, ied, jsd, jed)
    call safe_alloc_ptr(fluxes%sens, isd, ied, jsd, jed)
  else ! This is the buoyancy only mode.
    call safe_alloc_ptr(fluxes%buoy, isd, ied, jsd, jed)
  endif


  ! MODIFY THE CODE IN THE FOLLOWING LOOPS TO SET THE BUOYANCY FORCING TERMS.

  if ( CS%use_temperature ) then
    ! Set whichever fluxes are to be used here.  Any fluxes that
    ! are always zero do not need to be changed here.
    do j=js,je ; do i=is,ie
      ! Fluxes of fresh water through the surface are in units of [R Z T-1 ~> kg m-2 s-1]
      ! and are positive downward - i.e. evaporation should be negative.
      fluxes%evap(i,j) = -0.0 * G%mask2dT(i,j)
      fluxes%lprec(i,j) = 0.0 * G%mask2dT(i,j)

      ! vprec will be set later, if it is needed for salinity restoring.
      fluxes%vprec(i,j) = 0.0

      !   Heat fluxes are in units of [Q R Z T-1 ~> W m-2] and are positive into the ocean.
      fluxes%lw(i,j) = 0.0 * G%mask2dT(i,j)
      fluxes%latent(i,j) = 0.0 * G%mask2dT(i,j)
      fluxes%sens(i,j) = 0.0 * G%mask2dT(i,j)
      fluxes%sw(i,j) = 0.0 * G%mask2dT(i,j)
    enddo ; enddo
  else ! This is the buoyancy only mode.
    do j=js,je ; do i=is,ie
      !   fluxes%buoy is the buoyancy flux into the ocean [L2 T-3 ~> m2 s-3].  A positive
      ! buoyancy flux is of the same sign as heating the ocean.
      fluxes%buoy(i,j) = 0.0 * G%mask2dT(i,j)
    enddo ; enddo
  endif

  if (CS%restorebuoy) then
    if (CS%use_temperature) then
      call safe_alloc_ptr(fluxes%heat_added, isd, ied, jsd, jed)
      !   When modifying the code, comment out this error message.  It is here
      ! so that the original (unmodified) version is not accidentally used.
      call MOM_error(FATAL, "User_buoyancy_surface_forcing: " // &
        "Temperature and salinity restoring used without modification." )

      rhoXcp = CS%Rho0 * fluxes%C_p
      do j=js,je ; do i=is,ie
        !   Set Temp_restore and Salin_restore to the temperature (in degC) and
        ! salinity (in PSU or ppt) that are being restored toward.
        Temp_restore = 0.0
        Salin_restore = 0.0

        fluxes%heat_added(i,j) = (G%mask2dT(i,j) * (rhoXcp * CS%Flux_const)) * &
            (Temp_restore - sfc_state%SST(i,j))
        fluxes%vprec(i,j) = - (G%mask2dT(i,j) * (CS%Rho0*CS%Flux_const)) * &
            ((Salin_restore - sfc_state%SSS(i,j)) / (0.5 * (Salin_restore + sfc_state%SSS(i,j))))
      enddo ; enddo
    else
      !   When modifying the code, comment out this error message.  It is here
      ! so that the original (unmodified) version is not accidentally used.
      call MOM_error(FATAL, "User_buoyancy_surface_forcing: " // &
        "Buoyancy restoring used without modification." )

      ! The -1 is because density has the opposite sign to buoyancy.
      buoy_rest_const = -1.0 * (CS%G_Earth * CS%Flux_const) / CS%Rho0
      do j=js,je ; do i=is,ie
       !   Set density_restore to an expression for the surface potential
       ! density [kg m-3] that is being restored toward.
        density_restore = 1030.0*US%kg_m3_to_R

        fluxes%buoy(i,j) = G%mask2dT(i,j) * buoy_rest_const * &
                          (density_restore - sfc_state%sfc_density(i,j))
      enddo ; enddo
    endif
  endif                                             ! end RESTOREBUOY

end subroutine USER_buoyancy_forcing

!> This subroutine initializes the USER_surface_forcing module
subroutine USER_surface_forcing_init(Time, G, US, param_file, diag, CS)
  type(time_type),               intent(in) :: Time !< The current model time
  type(ocean_grid_type),         intent(in) :: G    !< The ocean's grid structure
  type(unit_scale_type),         intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),         intent(in) :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,       intent(in) :: diag !< A structure that is used to regulate diagnostic output.
  type(user_surface_forcing_CS), pointer    :: CS   !< A pointer that is set to point to
                                                    !! the control structure for this module

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "user_surface_forcing" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "USER_surface_forcing_init called with an associated "// &
                             "control structure.")
    return
  endif
  allocate(CS)
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state "//&
                 "variables.", default=.true.)

  call get_param(param_file, mdl, "G_EARTH", CS%G_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80, scale=US%m_to_L**2*US%Z_to_m*US%T_to_s**2)
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINSEQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "GUST_CONST", CS%gust_const, &
                 "The background gustiness in the winds.", &
                 units="Pa", default=0.0, scale=US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z)

  call get_param(param_file, mdl, "RESTOREBUOY", CS%restorebuoy, &
                 "If true, the buoyancy fluxes drive the model back "//&
                 "toward some specified surface state with a rate "//&
                 "given by FLUXCONST.", default= .false.)
  if (CS%restorebuoy) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 default=0.0, units="m day-1", scale=US%m_to_Z/(86400.0*US%s_to_T))
  endif

end subroutine USER_surface_forcing_init

!! \namespace user_surface_forcing
!!
!! This file contains the subroutines that a user should modify to
!! to set the surface wind stresses and fluxes of buoyancy or
!! temperature and fresh water.  They are called when the run-time
!! parameters WIND_CONFIG or BUOY_CONFIG are set to "USER".  The
!! standard version has simple examples, along with run-time error
!! messages that will cause the model to abort if this code has no
!! been modified.  This code is intended for use with relatively
!! simple specifications of the forcing.  For more complicated forms,
!! it is probably a good idea to read the forcing from input files
!! using "file" for WIND_CONFIG and BUOY_CONFIG.
!!
!! USER_wind_forcing() should set the surface wind stresses (taux and
!! tauy) perhaps along with the surface friction velocity (ustar).
!!
!! USER_buoyancy() forcing is used to set the surface buoyancy
!! forcing, which may include a number of fresh water flux fields
!! (evap, lprec, fprec, lrunoff, frunoff, and
!! vprec) and the surface heat fluxes (sw, lw, latent and sens)
!! if temperature and salinity are state variables, or it may simply
!! be the buoyancy flux if it is not.  This routine also has coded a
!! restoring to surface values of temperature and salinity.

end module user_surface_forcing
