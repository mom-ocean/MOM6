!> Wind and buoyancy forcing for the Neverland configurations
module Neverland_surface_forcing

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, query_averaging_enabled
use MOM_diag_mediator, only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_domains, only : pass_var, pass_vector, AGRID
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, mech_forcing
use MOM_forcing_type, only : allocate_forcing_type, allocate_mech_forcing
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, read_data, slasher
use MOM_time_manager, only : time_type, operator(+), operator(/)
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface

implicit none ; private

public Neverland_wind_forcing
public Neverland_buoyancy_forcing
public Neverland_surface_forcing_init

!> This control structure should be used to store any run-time variables
!! associated with the Neverland forcing.
!!
!! It can be readily modified for a specific case, and because it is private there
!! will be no changes needed in other code (although they will have to be recompiled).
type, public :: Neverland_surface_forcing_CS ; private

  logical :: use_temperature !< If true, use temperature and salinity.
  logical :: restorebuoy     !< If true, use restoring surface buoyancy forcing.
  real :: Rho0               !< The density used in the Boussinesq approximation [R ~> kg m-3].
  real :: G_Earth            !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real :: flux_const         !<  The restoring rate at the surface [Z T-1 ~> m s-1].
  real, dimension(:,:), pointer :: &
    buoy_restore(:,:) => NULL() !< The pattern to restore buoyancy to.
  character(len=200) :: inputdir !< The directory where NetCDF input files are.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.
  logical :: first_call = .true. !< True until Neverland_buoyancy_forcing has been called
end type Neverland_surface_forcing_CS

contains

!> Sets the surface wind stresses, forces%taux and forces%tauy for the
!! Neverland forcing configuration.
subroutine Neverland_wind_forcing(sfc_state, forces, day, G, US, CS)
  type(surface),                 intent(inout) :: sfc_state !< A structure containing fields that
                                                         !! describe the surface state of the ocean.
  type(mech_forcing),            intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),               intent(in)    :: day    !< Time used for determining the fluxes.
  type(ocean_grid_type),         intent(inout) :: G      !< Grid structure.
  type(unit_scale_type),         intent(in)    :: US     !< A dimensional unit scaling type
  type(Neverland_surface_forcing_CS), pointer  :: CS     !< Control structure for this module.

  ! Local variables
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real :: x, y
  real :: PI
  real :: tau_max  ! The magnitude of the wind stress [R Z L T-2 ~> Pa]
  real :: off

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! Allocate the forcing arrays, if necessary.
  call allocate_mech_forcing(G, forces, stress=.true.)

  !  Set the surface wind stresses, in units of Pa.  A positive taux
  !  accelerates the ocean to the (pseudo-)east.

  !  The i-loop extends to is-1 so that taux can be used later in the
  ! calculation of ustar - otherwise the lower bound would be Isq.
  PI = 4.0*atan(1.0)
  forces%taux(:,:) = 0.0
  tau_max = 0.2 * US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z
  off = 0.02
  do j=js,je ; do I=is-1,Ieq
!    x = (G%geoLonT(i,j)-G%west_lon)/G%len_lon
    y = (G%geoLatT(i,j)-G%south_lat)/G%len_lat
!    forces%taux(I,j) = G%mask2dCu(I,j) * 0.0

    if (y <= 0.29) then
      forces%taux(I,j) = forces%taux(I,j) + tau_max * ( (1/0.29)*y - ( 1/(2*PI) )*sin( (2*PI*y) / 0.29 ) )
    endif
    if ((y > 0.29) .and. (y <= (0.8-off))) then
      forces%taux(I,j) = forces%taux(I,j) + tau_max *(0.35+0.65*cos(PI*(y-0.29)/(0.51-off))  )
    endif
    if ((y > (0.8-off)) .and. (y <= (1-off))) then
      forces%taux(I,j) = forces%taux(I,j) + tau_max *( 1.5*( (y-1+off) - (0.1/PI)*sin(10.0*PI*(y-0.8+off)) ) )
    endif
  enddo ; enddo

  do J=js-1,Jeq ; do i=is,ie
    forces%tauy(i,J) = G%mask2dCv(i,J) * 0.0  ! Change this to the desired expression.
  enddo ; enddo

  !    Set the surface friction velocity, in units of m s-1.  ustar
  !  is always positive.
! if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
!   !  This expression can be changed if desired, but need not be.
!   forces%ustar(i,j) = G%mask2dT(i,j) * sqrt((CS%gust_const + &
!           sqrt(0.5*(forces%taux(I-1,j)**2 + forces%taux(I,j)**2) + &
!                0.5*(forces%tauy(i,J-1)**2 + forces%tauy(i,J)**2))) * &
!            (US%L_to_Z / CS%Rho0) )
! enddo ; enddo ; endif

end subroutine Neverland_wind_forcing

!> Returns the value of a cosine-bell function evaluated at x/L
real function cosbell(x,L)

  real , intent(in) :: x       !< non-dimensional position
  real , intent(in) :: L       !< non-dimensional width
  real              :: PI      !< 3.1415926... calculated as 4*atan(1)

  PI      = 4.0*atan(1.0)
  cosbell = 0.5 * (1 + cos(PI*MIN(ABS(x/L),1.0)))
end function cosbell

!> Returns the value of a sin-spike function evaluated at x/L
real function spike(x,L)

  real , intent(in) :: x       !< non-dimensional position
  real , intent(in) :: L       !< non-dimensional width
  real              :: PI      !< 3.1415926... calculated as 4*atan(1)

  PI    = 4.0*atan(1.0)
  spike = (1 - sin(PI*MIN(ABS(x/L),0.5)))
end function spike


!> Surface fluxes of buoyancy for the Neverland configurations.
subroutine Neverland_buoyancy_forcing(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),                 intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),                 intent(inout) :: fluxes !< Forcing fields.
  type(time_type),               intent(in)    :: day !< Time used for determining the fluxes.
  real,                          intent(in)    :: dt !< Forcing time step (s).
  type(ocean_grid_type),         intent(inout) :: G  !< Grid structure.
  type(unit_scale_type),         intent(in)    :: US !< A dimensional unit scaling type
  type(Neverland_surface_forcing_CS), pointer  :: CS !< Control structure for this module.
  ! Local variables
  real :: buoy_rest_const  ! A constant relating density anomalies to the
                           ! restoring buoyancy flux [L2 T-3 R-1 ~> m5 s-3 kg-1].
  real :: density_restore  ! De
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  !   When modifying the code, comment out this error message.  It is here
  ! so that the original (unmodified) version is not accidentally used.

  ! Allocate and zero out the forcing arrays, as necessary.  This portion is
  ! usually not changed.
  if (CS%use_temperature) then
    call MOM_error(FATAL, "Neverland_buoyancy_forcing: " // &
        "Temperature and salinity mode not coded!" )
  else
    ! This is the buoyancy only mode.
    call safe_alloc_ptr(fluxes%buoy, isd, ied, jsd, jed)
  endif


  ! MODIFY THE CODE IN THE FOLLOWING LOOPS TO SET THE BUOYANCY FORCING TERMS.
  if (CS%restorebuoy .and. CS%first_call) then
    call safe_alloc_ptr(CS%buoy_restore, isd, ied, jsd, jed)
    CS%first_call = .false.
    ! Set CS%buoy_restore(i,j) here
  endif

  if ( CS%use_temperature ) then
    call MOM_error(FATAL, "Neverland_buoyancy_surface_forcing: " // &
      "Temperature/salinity restoring not coded!" )
  else ! This is the buoyancy only mode.
    do j=js,je ; do i=is,ie
      !   fluxes%buoy is the buoyancy flux into the ocean [L2 T-3 ~> m2 s-3].  A positive
      ! buoyancy flux is of the same sign as heating the ocean.
      fluxes%buoy(i,j) = 0.0 * G%mask2dT(i,j)
    enddo ; enddo
  endif

  if (CS%restorebuoy) then
    if (CS%use_temperature) then
      call MOM_error(FATAL, "Neverland_buoyancy_surface_forcing: " // &
        "Temperature/salinity restoring not coded!" )
    else
      !   When modifying the code, comment out this error message.  It is here
      ! so that the original (unmodified) version is not accidentally used.

      ! The -1 is because density has the opposite sign to buoyancy.
      buoy_rest_const = -1.0 * (CS%G_Earth * CS%Flux_const) / CS%Rho0
      do j=js,je ; do i=is,ie
       !   Set density_restore to an expression for the surface potential
       ! density [kg m-3] that is being restored toward.
        density_restore = 1030.0

        fluxes%buoy(i,j) = G%mask2dT(i,j) * buoy_rest_const * &
                          US%kg_m3_to_R*(density_restore - sfc_state%sfc_density(i,j))
      enddo ; enddo
    endif
  endif                                             ! end RESTOREBUOY

end subroutine Neverland_buoyancy_forcing

!> Initializes the Neverland control structure.
subroutine Neverland_surface_forcing_init(Time, G, US, param_file, diag, CS)
  type(time_type),         intent(in) :: Time       !< The current model time.
  type(ocean_grid_type),   intent(in) :: G          !< The ocean's grid structure.
  type(unit_scale_type),   intent(in) :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< A structure indicating the open file to parse for
                                                    !! model parameter values.
  type(diag_ctrl), target, intent(in) :: diag       !< A structure that is used to regulate diagnostic output.
  type(Neverland_surface_forcing_CS), pointer :: CS !< A pointer that is set to point to the control structure
                                                    !! for this module
  ! This include declares and sets the variable "version".
#include "version_variable.h"
  ! Local variables
  character(len=40) :: mdl = "Neverland_surface_forcing" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "Neverland_surface_forcing_init called with an associated "// &
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
! call get_param(param_file, mdl, "GUST_CONST", CS%gust_const, &
!                "The background gustiness in the winds.", units="Pa", &
!                default=0.02)

  call get_param(param_file, mdl, "RESTOREBUOY", CS%restorebuoy, &
                 "If true, the buoyancy fluxes drive the model back "//&
                 "toward some specified surface state with a rate "//&
                 "given by FLUXCONST.", default= .false.)

  if (CS%restorebuoy) then
    call get_param(param_file, mdl, "FLUXCONST", CS%flux_const, &
                 "The constant that relates the restoring surface fluxes "//&
                 "to the relative surface anomalies (akin to a piston "//&
                 "velocity).  Note the non-MKS units.", &
                 units="m day-1", scale=US%m_to_Z*US%T_to_s, fail_if_missing=.true.)
    ! Convert CS%flux_const from m day-1 to m s-1.
    CS%flux_const = CS%flux_const / 86400.0
  endif

end subroutine Neverland_surface_forcing_init

end module Neverland_surface_forcing
