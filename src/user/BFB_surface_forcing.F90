!> Surface forcing for the boundary-forced-basin (BFB) configuration
module BFB_surface_forcing

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, query_averaging_enabled
use MOM_diag_mediator, only : register_diag_field, diag_ctrl
use MOM_domains, only : pass_var, pass_vector, AGRID
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, param_file_type, log_version
use MOM_forcing_type, only : forcing, allocate_forcing_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, read_data
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/)
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface

implicit none ; private

public BFB_buoyancy_forcing, BFB_surface_forcing_init

!> Control structure for BFB_surface_forcing
type, public :: BFB_surface_forcing_CS ; private

  logical :: use_temperature !< If true, temperature and salinity are used as state variables.
  logical :: restorebuoy     !< If true, use restoring surface buoyancy forcing.
  real :: Rho0               !< The density used in the Boussinesq approximation [R ~> kg m-3].
  real :: G_Earth            !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real :: Flux_const         !< The restoring rate at the surface [Z T-1 ~> m s-1].
  real :: gust_const         !< A constant unresolved background gustiness
                             !! that contributes to ustar [Pa].
  real :: SST_s              !< SST at the southern edge of the linear forcing ramp [degC]
  real :: SST_n              !< SST at the northern edge of the linear forcing ramp [degC]
  real :: lfrslat            !< Southern latitude where the linear forcing ramp begins [degLat]
  real :: lfrnlat            !< Northern latitude where the linear forcing ramp ends [degLat]
  real :: drho_dt            !< Rate of change of density with temperature [R degC-1 ~> kg m-3 degC-1].
                             !!   Note that temperature is being used as a dummy variable here.
                             !! All temperatures are converted into density.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                             !! regulate the timing of diagnostic output.
end type BFB_surface_forcing_CS

contains

!> Bouyancy forcing for the boundary-forced-basin (BFB) configuration
subroutine BFB_buoyancy_forcing(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),                intent(inout) :: sfc_state  !< A structure containing fields that
                                                      !! describe the surface state of the ocean.
  type(forcing),                intent(inout) :: fluxes !< A structure containing pointers to any
                                                      !! possible forcing fields. Unused fields
                                                      !! have NULL ptrs.
  type(time_type),              intent(in)    :: day  !< Time of the fluxes.
  real,                         intent(in)    :: dt   !< The amount of time over which
                                                      !! the fluxes apply [s]
  type(ocean_grid_type),        intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),        intent(in)    :: US   !< A dimensional unit scaling type
  type(BFB_surface_forcing_CS), pointer       :: CS   !< A pointer to the control structure
                                                      !! returned by a previous call to
                                                      !! BFB_surface_forcing_init.
  ! Local variables
  real :: Temp_restore   ! The temperature that is being restored toward [degC].
  real :: Salin_restore  ! The salinity that is being restored toward [ppt].
  real :: density_restore  ! The potential density that is being restored
                         ! toward [R ~> kg m-3].
  real :: rhoXcp           ! Reference density times heat capacity times unit scaling
                           ! factors [Q R degC-1 ~> J m-3 degC-1]
  real :: buoy_rest_const  ! A constant relating density anomalies to the
                           ! restoring buoyancy flux [L2 T-3 R-1 ~> m5 s-3 kg-1].
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

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

      ! Heat fluxes are in units of [Q R Z T-1 ~> W m-2] and are positive into the ocean.
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
        ! salinity (in ppt) that are being restored toward.
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
      ! call MOM_error(FATAL, "User_buoyancy_surface_forcing: " // &
      !   "Buoyancy restoring used without modification." )

      ! The -1 is because density has the opposite sign to buoyancy.
      buoy_rest_const = -1.0 * (CS%G_Earth * CS%Flux_const) / CS%Rho0
      Temp_restore = 0.0
      do j=js,je ; do i=is,ie
       !   Set density_restore to an expression for the surface potential
       ! density [R ~> kg m-3] that is being restored toward.
        if (G%geoLatT(i,j) < CS%lfrslat) then
            Temp_restore = CS%SST_s
        elseif (G%geoLatT(i,j) > CS%lfrnlat) then
            Temp_restore = CS%SST_n
        else
            Temp_restore = (CS%SST_s - CS%SST_n)/(CS%lfrslat - CS%lfrnlat) * &
                    (G%geoLatT(i,j) - CS%lfrslat) + CS%SST_s
        endif

        density_restore = Temp_restore*CS%drho_dt + CS%Rho0

        fluxes%buoy(i,j) = G%mask2dT(i,j) * buoy_rest_const * &
                          (density_restore - sfc_state%sfc_density(i,j))
      enddo ; enddo
    endif
  endif                                             ! end RESTOREBUOY

end subroutine BFB_buoyancy_forcing

!> Initialization for forcing the boundary-forced-basin (BFB) configuration
subroutine BFB_surface_forcing_init(Time, G, US, param_file, diag, CS)
  type(time_type),              intent(in) :: Time !< The current model time.
  type(ocean_grid_type),        intent(in) :: G    !< The ocean's grid structure
  type(unit_scale_type),        intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),        intent(in) :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,      intent(in) :: diag !< A structure that is used to
                                                   !! regulate diagnostic output.
  type(BFB_surface_forcing_CS), pointer    :: CS   !< A pointer to the control structure for this module
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "BFB_surface_forcing" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "BFB_surface_forcing_init called with an associated "// &
                             "control structure.")
    return
  endif
  allocate(CS)
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state variables.", default=.true.)

  call get_param(param_file, mdl, "G_EARTH", CS%G_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80, scale=US%m_to_L**2*US%Z_to_m*US%T_to_s**2)
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINSEQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "LFR_SLAT", CS%lfrslat, &
                 "Southern latitude where the linear forcing ramp begins.", &
                 units="degrees", default=20.0)
  call get_param(param_file, mdl, "LFR_NLAT", CS%lfrnlat, &
                 "Northern latitude where the linear forcing ramp ends.", &
                 units="degrees", default=40.0)
  call get_param(param_file, mdl, "SST_S", CS%SST_s, &
                 "SST at the southern edge of the linear forcing ramp.", &
                 units="C", default=20.0)
  call get_param(param_file, mdl, "SST_N", CS%SST_n, &
                 "SST at the northern edge of the linear forcing ramp.", &
                 units="C", default=10.0)
  call get_param(param_file, mdl, "DRHO_DT", CS%drho_dt, &
                 "The rate of change of density with temperature.", &
                 units="kg m-3 K-1", default=-0.2, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "GUST_CONST", CS%gust_const, &
                 "The background gustiness in the winds.", units="Pa", &
                 default=0.02)

  call get_param(param_file, mdl, "RESTOREBUOY", CS%restorebuoy, &
                 "If true, the buoyancy fluxes drive the model back "//&
                 "toward some specified surface state with a rate "//&
                 "given by FLUXCONST.", default= .false.)
  if (CS%restorebuoy) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes "//&
                 "to the relative surface anomalies (akin to a piston "//&
                 "velocity).  Note the non-MKS units.", &
                 units="m day-1", scale=US%m_to_Z*US%T_to_s, fail_if_missing=.true.)
    ! Convert CS%Flux_const from m day-1 to m s-1.
    CS%Flux_const = CS%Flux_const / 86400.0
  endif

end subroutine BFB_surface_forcing_init

end module BFB_surface_forcing
