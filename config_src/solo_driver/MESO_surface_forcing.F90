!> Sets forcing for the MESO configuration
module MESO_surface_forcing

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, query_averaging_enabled
use MOM_diag_mediator, only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_domains, only : pass_var, pass_vector, AGRID
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, mech_forcing
use MOM_forcing_type, only : allocate_forcing_type, allocate_mech_forcing
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, MOM_read_data, slasher
use MOM_time_manager, only : time_type, operator(+), operator(/)
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface

implicit none ; private

public MESO_buoyancy_forcing, MESO_surface_forcing_init

!> This control structure is used to store parameters associated with the MESO forcing.
type, public :: MESO_surface_forcing_CS ; private

  logical :: use_temperature !< If true, temperature and salinity are used as state variables.
  logical :: restorebuoy     !< If true, use restoring surface buoyancy forcing.
  real :: Rho0               !< The density used in the Boussinesq approximation [R ~> kg m-3].
  real :: G_Earth            !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real :: Flux_const         !< The restoring rate at the surface [Z T-1 ~> m s-1].
  real :: gust_const         !< A constant unresolved background gustiness
                             !! that contributes to ustar [Pa].
  real, dimension(:,:), pointer :: &
    T_Restore(:,:) => NULL(), & !< The temperature to restore the SST toward [degC].
    S_Restore(:,:) => NULL(), & !< The salinity to restore the sea surface salnity toward [ppt]
    PmE(:,:) => NULL(), &       !< The prescribed precip minus evap [Z T-1 ~> m s-1].
    Solar(:,:) => NULL()        !< The shortwave forcing into the ocean [Q R Z T-1 ~> W m-2].
  real, dimension(:,:), pointer :: Heat(:,:) => NULL() !< The prescribed longwave, latent and sensible
                                !! heat flux into the ocean [Q R Z T-1 ~> W m-2].
  character(len=200) :: inputdir !< The directory where NetCDF input files are.
  character(len=200) :: salinityrestore_file !< The file with the target sea surface salinity
  character(len=200) :: SSTrestore_file !< The file with the target sea surface temperature
  character(len=200) :: Solar_file !< The file with the shortwave forcing
  character(len=200) :: heating_file !< The file with the longwave, latent, and sensible heating
  character(len=200) :: PmE_file !< The file with precipitation minus evaporation
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                             !! timing of diagnostic output.
end type MESO_surface_forcing_CS

logical :: first_call = .true. !< True until after the first call to the MESO forcing routines

contains

!> This subroutine sets up the MESO buoyancy forcing, which uses control-theory style
!! specification restorative buoyancy fluxes at large scales.
subroutine MESO_buoyancy_forcing(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),                 intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),                 intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),               intent(in)    :: day  !< The time of the fluxes
  real,                          intent(in)    :: dt   !< The amount of time over which
                                                       !! the fluxes apply [s]
  type(ocean_grid_type),         intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),         intent(in)    :: US   !< A dimensional unit scaling type
  type(MESO_surface_forcing_CS), pointer       :: CS   !< A pointer to the control structure returned by
                                                       !! a previous call to MESO_surface_forcing_init

!    When temperature is used, there are long list of fluxes that need to be
!  set - essentially the same as for a full coupled model, but most of these
!  can be simply set to zero.  The net fresh water flux should probably be
!  set in fluxes%evap and fluxes%lprec, with any salinity restoring
!  appearing in fluxes%vprec, and the other water flux components
!  (fprec, lrunoff and frunoff) left as arrays full of zeros.
!  Evap is usually negative and precip is usually positive.  All heat fluxes
!  are in W m-2 and positive for heat going into the ocean.  All fresh water
!  fluxes are in kg m-2 s-1 and positive for water moving into the ocean.

  real :: Temp_restore   ! The temperature that is being restored toward [degC].
  real :: Salin_restore  ! The salinity that is being restored toward [ppt]
  real :: density_restore  ! The potential density that is being restored toward [R ~> kg m-3].
  real :: rhoXcp ! The mean density times the heat capacity [Q R degC-1 ~> J m-3 degC-1].
  real :: buoy_rest_const  ! A constant relating density anomalies to the
                           ! restoring buoyancy flux [L2 T-3 R-1 ~> m5 s-3 kg-1].

  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  !   When modifying the code, comment out this error message.  It is here
  ! so that the original (unmodified) version is not accidentally used.

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
    call safe_alloc_ptr(fluxes%heat_content_lprec, isd, ied, jsd, jed)
  else ! This is the buoyancy only mode.
    call safe_alloc_ptr(fluxes%buoy, isd, ied, jsd, jed)
  endif


  ! MODIFY THE CODE IN THE FOLLOWING LOOPS TO SET THE BUOYANCY FORCING TERMS.
  if (CS%restorebuoy .and. first_call) then !#CTRL# .or. associated(CS%ctrl_forcing_CSp)) then
    call safe_alloc_ptr(CS%T_Restore, isd, ied, jsd, jed)
    call safe_alloc_ptr(CS%S_Restore, isd, ied, jsd, jed)
    call safe_alloc_ptr(CS%Heat, isd, ied, jsd, jed)
    call safe_alloc_ptr(CS%PmE, isd, ied, jsd, jed)
    call safe_alloc_ptr(CS%Solar, isd, ied, jsd, jed)

    call MOM_read_data(trim(CS%inputdir)//trim(CS%SSTrestore_file), "SST", &
             CS%T_Restore(:,:), G%Domain)
    call MOM_read_data(trim(CS%inputdir)//trim(CS%salinityrestore_file), "SAL", &
             CS%S_Restore(:,:), G%Domain)
    call MOM_read_data(trim(CS%inputdir)//trim(CS%heating_file), "Heat", &
             CS%Heat(:,:), G%Domain, scale=US%W_m2_to_QRZ_T)
    call MOM_read_data(trim(CS%inputdir)//trim(CS%PmE_file), "PmE", &
             CS%PmE(:,:), G%Domain, scale=US%m_to_Z*US%T_to_s)
    call MOM_read_data(trim(CS%inputdir)//trim(CS%Solar_file), "NET_SOL", &
             CS%Solar(:,:), G%Domain, scale=US%W_m2_to_QRZ_T)
    first_call = .false.
  endif

  if ( CS%use_temperature ) then
    ! Set whichever fluxes are to be used here.  Any fluxes that
    ! are always zero do not need to be changed here.
    do j=js,je ; do i=is,ie
      ! Fluxes of fresh water through the surface are in units of [kg m-2 s-1]
      ! and are positive downward - i.e. evaporation should be negative.
      fluxes%evap(i,j) = -0.0 * G%mask2dT(i,j)
      fluxes%lprec(i,j) =  CS%PmE(i,j) * CS%Rho0 * G%mask2dT(i,j)

      ! vprec will be set later, if it is needed for salinity restoring.
      fluxes%vprec(i,j) = 0.0

      !   Heat fluxes are in units of [Q R Z T-1 ~> W m-2] and are positive into the ocean.
      fluxes%lw(i,j)     = 0.0 * G%mask2dT(i,j)
      fluxes%latent(i,j) = 0.0 * G%mask2dT(i,j)
      fluxes%sens(i,j)   = CS%Heat(i,j) * G%mask2dT(i,j)
      fluxes%sw(i,j)     = CS%Solar(i,j) * G%mask2dT(i,j)
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
!      call MOM_error(FATAL, "MESO_buoyancy_surface_forcing: " // &
!        "Temperature and salinity restoring used without modification." )

      rhoXcp = CS%Rho0 * fluxes%C_p
      do j=js,je ; do i=is,ie
        !   Set Temp_restore and Salin_restore to the temperature (in degC) and
        ! salinity (in ppt or PSU) that are being restored toward.
        if (G%mask2dT(i,j) > 0) then
          fluxes%heat_added(i,j) = G%mask2dT(i,j) * &
              ((CS%T_Restore(i,j) - sfc_state%SST(i,j)) * rhoXcp * CS%Flux_const)
          fluxes%vprec(i,j) = - (CS%Rho0 * CS%Flux_const) * &
              (CS%S_Restore(i,j) - sfc_state%SSS(i,j)) / &
              (0.5*(sfc_state%SSS(i,j) + CS%S_Restore(i,j)))
        else
          fluxes%heat_added(i,j) = 0.0
          fluxes%vprec(i,j) = 0.0
        endif
      enddo ; enddo
    else
      !   When modifying the code, comment out this error message.  It is here
      ! so that the original (unmodified) version is not accidentally used.
      call MOM_error(FATAL, "MESO_buoyancy_surface_forcing: " // &
        "Buoyancy restoring used without modification." )

      ! The -1 is because density has the opposite sign to buoyancy.
      buoy_rest_const = -1.0 * (CS%G_Earth * CS%Flux_const) / CS%Rho0
      do j=js,je ; do i=is,ie
       !   Set density_restore to an expression for the surface potential
       ! density [R ~> kg m-3] that is being restored toward.
        density_restore = 1030.0 * US%kg_m3_to_R

        fluxes%buoy(i,j) = G%mask2dT(i,j) * buoy_rest_const * &
                           (density_restore - sfc_state%sfc_density(i,j))
      enddo ; enddo
    endif
  endif                                             ! end RESTOREBUOY

end subroutine MESO_buoyancy_forcing

!> Initialize the MESO surface forcing module
subroutine MESO_surface_forcing_init(Time, G, US, param_file, diag, CS)

  type(time_type),               intent(in)    :: Time !< The current model time
  type(ocean_grid_type),         intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),         intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),         intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,       intent(inout) :: diag !< structure used to regulate diagnostic output
  type(MESO_surface_forcing_CS), pointer       :: CS   !< A pointer that is set to point to the
                                                       !! control structure for this module

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MESO_surface_forcing" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "MESO_surface_forcing_init called with an associated "// &
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
                 "The background gustiness in the winds.", units="Pa", default=0.0)

  call get_param(param_file, mdl, "RESTOREBUOY", CS%restorebuoy, &
                 "If true, the buoyancy fluxes drive the model back "//&
                 "toward some specified surface state with a rate "//&
                 "given by FLUXCONST.", default= .false.)

  if (CS%restorebuoy) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 default=0.0, units="m day-1", scale=US%m_to_Z/(86400.0*US%s_to_T))

    call get_param(param_file, mdl, "SSTRESTORE_FILE", CS%SSTrestore_file, &
                 "The file with the SST toward which to restore in "//&
                 "variable TEMP.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "SALINITYRESTORE_FILE", CS%salinityrestore_file, &
                 "The file with the surface salinity toward which to "//&
                 "restore in variable SALT.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "SENSIBLEHEAT_FILE", CS%heating_file, &
                 "The file with the non-shortwave heat flux in "//&
                 "variable Heat.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "PRECIP_FILE", CS%PmE_file, &
                 "The file with the net precipiation minus evaporation "//&
                 "in variable PmE.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "SHORTWAVE_FILE", CS%Solar_file, &
                 "The file with the shortwave heat flux in "//&
                 "variable NET_SOL.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, default=".")
    CS%inputdir = slasher(CS%inputdir)

   endif

end subroutine MESO_surface_forcing_init

!> \namespace meso_surface_forcing
!!
!! Rewritten by Robert Hallberg, June 2009
!!
!!   This file contains the subroutines that a user should modify to
!! to set the surface wind stresses and fluxes of buoyancy or
!! temperature and fresh water.  They are called when the run-time
!! parameters WIND_CONFIG or BUOY_CONFIG are set to "USER".  The
!! standard version has simple examples, along with run-time error
!! messages that will cause the model to abort if this code has not
!! been modified.  This code is intended for use with relatively
!! simple specifications of the forcing.  For more complicated forms,
!! it is probably a good idea to read the forcing from input files
!! using "file" for WIND_CONFIG and BUOY_CONFIG.
!!
!!   MESO_buoyancy forcing is used to set the surface buoyancy
!! forcing, which may include a number of fresh water flux fields
!! (evap, liq_precip, froz_precip, liq_runoff, froz_runoff, and
!! vprec) and the surface heat fluxes (sw, lw, latent and sens)
!! if temperature and salinity are state variables, or it may simply
!! be the buoyancy flux if it is not.  This routine also has coded a
!! restoring to surface values of temperature and salinity.

end module MESO_surface_forcing
