!> Surface forcing for the dumbbell test case
module dumbbell_surface_forcing

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, query_averaging_enabled
use MOM_diag_mediator, only : register_diag_field, diag_ctrl
use MOM_domains, only : pass_var, pass_vector, AGRID
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, param_file_type, log_version
use MOM_forcing_type, only : forcing, allocate_forcing_type
use MOM_grid, only : ocean_grid_type
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/), get_time
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface

implicit none ; private

public dumbbell_dynamic_forcing, dumbbell_buoyancy_forcing, dumbbell_surface_forcing_init

!> Control structure for the dumbbell test case forcing
type, public :: dumbbell_surface_forcing_CS ; private
  logical :: use_temperature !< If true, temperature and salinity are used as state variables.
  logical :: restorebuoy     !< If true, use restoring surface buoyancy forcing.
  real :: Rho0               !< The density used in the Boussinesq approximation [R ~> kg m-3].
  real :: G_Earth            !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real :: Flux_const         !< The restoring rate at the surface [Z T-1 ~> m s-1].
! real :: gust_const         !< A constant unresolved background gustiness
!                            !! that contributes to ustar [R L Z T-2 ~> Pa].
  real :: slp_amplitude      !< The amplitude of pressure loading [R L2 T-2 ~> Pa] applied
                             !! to the reservoirs
  real :: slp_period         !< Period of sinusoidal pressure wave [days]
  real, dimension(:,:), allocatable :: &
    forcing_mask             !< A mask regulating where forcing occurs
  real, dimension(:,:), allocatable :: &
    S_restore                !< The surface salinity field toward which to restore [ppt].
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                             !! timing of diagnostic output.
end type dumbbell_surface_forcing_CS

contains

!> Surface buoyancy (heat and fresh water) fluxes for the dumbbell test case
subroutine dumbbell_buoyancy_forcing(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),                 intent(inout) :: sfc_state  !< A structure containing fields that
                                                         !! describe the surface state of the ocean.
  type(forcing),                 intent(inout) :: fluxes !< A structure containing pointers to any
                                                         !! possible forcing fields. Unused fields
                                                         !! have NULL ptrs.
  type(time_type),               intent(in)    :: day    !< Time of the fluxes.
  real,                          intent(in)    :: dt     !< The amount of time over which
                                                         !! the fluxes apply [s]
  type(ocean_grid_type),         intent(in)    :: G      !< The ocean's grid structure
  type(unit_scale_type),         intent(in)    :: US     !< A dimensional unit scaling type
  type(dumbbell_surface_forcing_CS),  pointer  :: CS     !< A control structure returned by a previous
                                                         !! call to dumbbell_surface_forcing_init
  ! Local variables
  real :: Temp_restore   ! The temperature that is being restored toward [degC].
  real :: Salin_restore  ! The salinity that is being restored toward [ppt].
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed


  ! Allocate and zero out the forcing arrays, as necessary.
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

  if (CS%use_temperature .and. CS%restorebuoy) then
    do j=js,je ; do i=is,ie
      if (CS%forcing_mask(i,j)>0.) then
        fluxes%vprec(i,j) = - (G%mask2dT(i,j) * (CS%Rho0*CS%Flux_const)) * &
                ((CS%S_restore(i,j) - sfc_state%SSS(i,j)) /  (0.5 * (CS%S_restore(i,j) + sfc_state%SSS(i,j))))

      endif
    enddo ; enddo
  endif

end subroutine dumbbell_buoyancy_forcing

!> Dynamic forcing for the dumbbell test case
subroutine dumbbell_dynamic_forcing(sfc_state, fluxes, day, dt, G, CS)
  type(surface),                 intent(inout) :: sfc_state  !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(forcing),                 intent(inout) :: fluxes !< A structure containing pointers to any
                                                       !! possible forcing fields. Unused fields
                                                       !! have NULL ptrs.
  type(time_type),               intent(in)    :: day  !< Time of the fluxes.
  real,                          intent(in)    :: dt   !< The amount of time over which
                                                       !! the fluxes apply [s]
  type(ocean_grid_type),         intent(in)    :: G    !< The ocean's grid structure
  type(dumbbell_surface_forcing_CS),  pointer  :: CS   !< A control structure returned by a previous
                                                       !! call to dumbbell_surface_forcing_init
  ! Local variables
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed
  integer :: idays, isecs
  real :: deg_rad, rdays


  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  deg_rad = atan(1.0)*4.0/180.

  call get_time(day,isecs,idays)
  rdays = real(idays) + real(isecs)/8.64e4
  ! This could be:  rdays = time_type_to_real(day)/8.64e4

  ! Allocate and zero out the forcing arrays, as necessary.
  call safe_alloc_ptr(fluxes%p_surf, isd, ied, jsd, jed)
  call safe_alloc_ptr(fluxes%p_surf_full, isd, ied, jsd, jed)

  do j=js,je ; do i=is,ie
    fluxes%p_surf(i,j) = CS%forcing_mask(i,j)* CS%slp_amplitude * &
                         G%mask2dT(i,j) * sin(deg_rad*(rdays/CS%slp_period))
    fluxes%p_surf_full(i,j) = CS%forcing_mask(i,j) * CS%slp_amplitude * &
                         G%mask2dT(i,j) * sin(deg_rad*(rdays/CS%slp_period))
  enddo ; enddo

end subroutine dumbbell_dynamic_forcing

!> Reads and sets up the forcing for the dumbbell test case
subroutine dumbbell_surface_forcing_init(Time, G, US, param_file, diag, CS)
  type(time_type),              intent(in) :: Time !< The current model time.
  type(ocean_grid_type),        intent(in) :: G    !< The ocean's grid structure
  type(unit_scale_type),        intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),        intent(in) :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl),      target, intent(in) :: diag !< A structure that is used to
                                                   !! regulate diagnostic output.
  type(dumbbell_surface_forcing_CS), &
                                pointer    :: CS   !< A pointer to the control structure for this module
  ! Local variables
  real :: S_surf, S_range
  real :: x, y
  integer :: i, j
  logical :: dbrotate    ! If true, rotate the domain.
#include "version_variable.h"
  character(len=40)  :: mdl = "dumbbell_surface_forcing" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "dumbbell_surface_forcing_init called with an associated "// &
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
  call get_param(param_file, mdl, "DUMBBELL_SLP_AMP", CS%slp_amplitude, &
                 "Amplitude of SLP forcing in reservoirs.", &
                 units="Pa", default = 10000.0, scale=US%kg_m3_to_R*US%m_s_to_L_T**2)
  call get_param(param_file, mdl, "DUMBBELL_SLP_PERIOD", CS%slp_period, &
                 "Periodicity of SLP forcing in reservoirs.", &
                 units="days", default = 1.0)
  call get_param(param_file, mdl, "DUMBBELL_ROTATION", dbrotate, &
                'Logical for rotation of dumbbell domain.',&
                 units='nondim', default=.false., do_not_log=.true.)
  call get_param(param_file, mdl,"INITIAL_SSS", S_surf, &
                 "Initial surface salinity", units="1e-3", default=34.0, do_not_log=.true.)
  call get_param(param_file, mdl,"INITIAL_S_RANGE", S_range, &
                 "Initial salinity range (bottom - surface)", units="1e-3", &
                 default=2., do_not_log=.true.)

  call get_param(param_file, mdl, "RESTOREBUOY", CS%restorebuoy, &
                 "If true, the buoyancy fluxes drive the model back "//&
                 "toward some specified surface state with a rate "//&
                 "given by FLUXCONST.", default= .false.)
  if (CS%restorebuoy) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 default=0.0, units="m day-1", scale=US%m_to_Z*US%T_to_s)
    ! Convert CS%Flux_const from m day-1 to m s-1.
    CS%Flux_const = CS%Flux_const / 86400.0


    allocate(CS%forcing_mask(G%isd:G%ied, G%jsd:G%jed)); CS%forcing_mask(:,:)=0.0
    allocate(CS%S_restore(G%isd:G%ied, G%jsd:G%jed))

    do j=G%jsc,G%jec
      do i=G%isc,G%iec
        ! Compute normalized zonal coordinates (x,y=0 at center of domain)
!       x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon - 0.5
!       y = ( G%geoLatT(i,j) - G%south_lat ) / G%len_lat - 0.5
        if (dbrotate) then
          ! This is really y in the rotated case
          x = ( G%geoLatT(i,j) - G%south_lat ) / G%len_lat - 0.5
        else
          x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon - 0.5
        endif
        CS%forcing_mask(i,j)=0
        CS%S_restore(i,j) = S_surf
        if ((x>0.25)) then
          CS%forcing_mask(i,j) = 1
          CS%S_restore(i,j) = S_surf + S_range
        elseif ((x<-0.25)) then
          CS%forcing_mask(i,j) = 1
          CS%S_restore(i,j) = S_surf - S_range
        endif
      enddo
    enddo
  endif

end subroutine dumbbell_surface_forcing_init

end module dumbbell_surface_forcing
