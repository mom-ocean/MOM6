module MOM_surface_forcing
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
!*  By Robert Hallberg, November 1998 - May 2002                       *
!*  Edited by Stephen Griffies June 2014                               *
!*                                                                     *
!*    This program contains the subroutines that calculate the         *
!*  surface wind stresses and fluxes of buoyancy or temperature and    *
!*  fresh water.  These subroutines will be called every time step,    *
!*  even if the wind stresses or buoyancy fluxes are constant in time  *
!*  - in that case these routines return quickly without doing         *
!*  anything.  In addition, any I/O of forcing fields is controlled    *
!*  by surface_forcing_init, located in this file.                     *
!*                                                                     *
!*    set_forcing is a small entry subroutine for the subroutines in   *
!*  this file.  It provides the external access to these subroutines.  *
!*                                                                     *
!*    wind_forcing determines the wind stresses and places them into   *
!*  taux[][] and tauy[][].  Often wind_forcing must be tailored for    *
!*  a particular application - either by specifying file and variable  *
!*  names or by providing appropriate internal expressions for the     *
!*  stresses.                                                          *
!*                                                                     *
!*    buoyancy_forcing determines the surface fluxes of buoyancy,      *
!*  temperature, and fresh water, as is appropriate.  A restoring      *
!*  boundary condition is implemented, but the code for any other      *
!*  boundary condition will usually be modified - either to specify    *
!*  file and variable names and which time level to read, or to set    *
!*  an internal expression for the variables.                          *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, tauy                                  *
!*    j    x ^ x ^ x   At >:  u, taux                                  *
!*    j    > o > o >   At o:  h, fluxes.                               *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**
!### use MOM_controlled_forcing, only : apply_ctrl_forcing, register_ctrl_forcing_restarts
!### use MOM_controlled_forcing, only : controlled_forcing_init, controlled_forcing_end
!### use MOM_controlled_forcing, only : ctrl_forcing_CS
use MOM_cpu_clock,           only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,           only : CLOCK_MODULE
use MOM_diag_mediator,       only : post_data, query_averaging_enabled
use MOM_diag_mediator,       only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_domains,             only : pass_var, pass_vector, AGRID, To_South, To_West, To_All
use MOM_error_handler,       only : callTree_enter, callTree_leave
use MOM_error_handler,       only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,         only : get_param, log_version, param_file_type
use MOM_string_functions,    only : uppercase
use MOM_forcing_type,        only : forcing, allocate_forcing_type, deallocate_forcing_type
use MOM_get_input,           only : Get_MOM_Input, directories
use MOM_grid,                only : ocean_grid_type
use MOM_io,                  only : file_exists, read_data, slasher
use MOM_restart,             only : register_restart_field, restart_init, MOM_restart_CS
use MOM_restart,             only : restart_init_end, save_restart, restore_state
use MOM_time_manager,        only : time_type, operator(+), operator(/), get_time, set_time
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_variables,           only : surface
! use MESO_surface_forcing,  only : MESO_wind_forcing, MESO_buoyancy_forcing
! use MESO_surface_forcing,  only : MESO_surface_forcing_init, MESO_surface_forcing_CS
use user_surface_forcing,    only : USER_wind_forcing, USER_buoyancy_forcing
use user_surface_forcing,    only : USER_surface_forcing_init, user_surface_forcing_CS
use user_revise_forcing,     only : user_alter_forcing, user_revise_forcing_init
use user_revise_forcing,     only : user_revise_forcing_CS

implicit none ; private

#include <MOM_memory.h>

public set_forcing
public surface_forcing_init
public forcing_diagnostics
public forcing_save_restart

! surface_forcing_CS is a structure containing pointers to the forcing fields
! which may be used to drive MOM.  All fluxes are positive into the ocean.
type, public :: surface_forcing_CS ; private

  logical :: use_temperature    ! if true, temp & salinity used as state variables
  logical :: restorebuoy        ! if true, use restoring surface buoyancy forcing
  logical :: adiabatic          ! if true, no diapycnal mass fluxes or surface buoyancy forcing
  logical :: variable_winds     ! if true, wind stresses vary with time
  logical :: variable_buoyforce ! if true, buoyancy forcing varies with time.
  real    :: south_lat          ! southern latitude of the domain
  real    :: len_lat            ! domain length in latitude

  real :: Rho0                  ! Boussinesq reference density (kg/m^3)
  real :: G_Earth               ! gravitational acceleration (m/s^2)
  real :: Flux_const            ! piston velocity for surface restoring (m/s)

  real    :: gust_const                 ! constant unresolved background gustiness for ustar (Pa)
  logical :: read_gust_2d               ! if true, use 2-dimensional gustiness supplied from a file
  real, pointer :: gust(:,:) => NULL()  ! spatially varying unresolved background gustiness (Pa)
                                        ! gust is used when read_gust_2d is true.

  real, pointer :: T_Restore(:,:)    => NULL()  ! temperature to damp (restore) the SST to (deg C)
  real, pointer :: S_Restore(:,:)    => NULL()  ! salinity to damp (restore) the SSS (g/kg)
  real, pointer :: Dens_Restore(:,:) => NULL()  ! density to damp (restore) surface density (kg/m^3)

  integer :: wind_last_lev_read = -1 ! The last time level read from the wind input files
  integer :: buoy_last_lev_read = -1 ! The last time level read from buoyancy input files

  real :: gyres_taux_const, gyres_taux_sin_amp, gyres_taux_cos_amp, gyres_taux_n_pis
                             ! if WIND_CONFIG=='gyres' then use
                             ! = A, B, C and n respectively for
                             ! taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L)

  real :: T_north, T_south   ! target temperatures at north and south used in
                             ! buoyancy_forcing_linear
  real :: S_north, S_south   ! target salinity at north and south used in
                             ! buoyancy_forcing_linear

  logical          :: first_call_set_forcing = .true.
  real             :: wind_scale    ! value by which wind-stresses are scaled (nondimensional)
  character(len=8) :: wind_stagger

  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp  => NULL()
!###  type(ctrl_forcing_CS), pointer    :: ctrl_forcing_CSp => NULL()
  type(MOM_restart_CS), pointer         :: restart_CSp      => NULL()

  type(diag_ctrl), pointer :: diag ! structure used to regulate timing of diagnostic output

  character(len=200) :: inputdir     ! The directory where NetCDF input files are.
  character(len=200) :: wind_config  ! Indicator for wind forcing type (2gyre, USER, FILE..)
  character(len=200) :: wind_file    ! If wind_config is "file", file to use
  character(len=200) :: buoy_config  ! Indicator for buoyancy forcing type
  character(len=200) :: longwavedown_file
  character(len=200) :: longwaveup_file
  character(len=200) :: evaporation_file
  character(len=200) :: sensibleheat_file
  character(len=200) :: shortwaveup_file
  character(len=200) :: shortwavedown_file
  character(len=200) :: snow_file
  character(len=200) :: precip_file
  character(len=200) :: freshdischarge_file
  character(len=200) :: SSTrestore_file
  character(len=200) :: salinityrestore_file
  character(len=80)  :: stress_x_var, stress_y_var

  type(user_revise_forcing_CS),   pointer :: urf_CS           => NULL()
  type(user_surface_forcing_CS),  pointer :: user_forcing_CSp => NULL()
!  type(MESO_surface_forcing_CS), pointer :: MESO_forcing_CSp => NULL()
end type surface_forcing_CS

integer :: id_clock_forcing

contains

subroutine set_forcing(state, fluxes, day_start, day_interval, G, CS)
  type(surface),         intent(inout) :: state
  type(forcing),         intent(inout) :: fluxes
  type(time_type),       intent(in)    :: day_start
  type(time_type),       intent(in)    :: day_interval
  type(ocean_grid_type), intent(inout) :: G
  type(surface_forcing_CS), pointer    :: CS

! This subroutine calls other subroutines in this file to get surface forcing fields.
! It also allocates and initializes the fields in the flux type.

! Arguments:
!  (inout)   state        = structure describing ocean surface state
!  (inout)   fluxes       = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day_start    = Start time of the fluxes
!  (in)      day_interval = Length of time over which these fluxes applied
!  (in)      G            = ocean grid structure
!  (in)      CS           = pointer to control struct returned by previous surface_forcing_init call

  real :: dt                     ! length of time in seconds over which fluxes applied
  type(time_type) :: day_center  ! central time of the fluxes.
  integer :: intdt

  call cpu_clock_begin(id_clock_forcing)

  day_center = day_start + day_interval/2
  call get_time(day_interval, intdt)
  dt = real(intdt)

  if (CS%variable_winds .or. CS%first_call_set_forcing) then
    if (trim(CS%wind_config) == "file") then
      call wind_forcing_from_file(state, fluxes, day_center, G, CS)
    elseif (trim(CS%wind_config) == "2gyre") then
      call wind_forcing_2gyre(state, fluxes, day_center, G, CS)
    elseif (trim(CS%wind_config) == "1gyre") then
      call wind_forcing_1gyre(state, fluxes, day_center, G, CS)
    elseif (trim(CS%wind_config) == "gyres") then
      call wind_forcing_gyres(state, fluxes, day_center, G, CS)
    elseif (trim(CS%wind_config) == "zero") then
      call wind_forcing_zero(state, fluxes, day_center, G, CS)
    elseif (trim(CS%wind_config) == "MESO") then
      call MOM_error(FATAL, "MESO forcing is not available with the ice-shelf"//&
               "version of MOM_surface_forcing.")
!      call MESO_wind_forcing(state, fluxes, day_center, G, CS%MESO_forcing_CSp)
    elseif (trim(CS%wind_config) == "USER") then
      call USER_wind_forcing(state, fluxes, day_center, G, CS%user_forcing_CSp)
    elseif (CS%variable_winds .and. .not.CS%first_call_set_forcing) then
      call MOM_error(FATAL, &
       "MOM_surface_forcing: Variable winds defined with no wind config")
    else
       call MOM_error(FATAL, &
       "MOM_surface_forcing:Unrecognized wind config "//trim(CS%wind_config))
    endif
  endif
  if ((CS%variable_buoyforce .or. CS%first_call_set_forcing) .and. &
      (.not.CS%adiabatic)) then
    if (trim(CS%buoy_config) == "file") then
      call buoyancy_forcing_from_files(state, fluxes, day_center, dt, G, CS)
    elseif (trim(CS%buoy_config) == "zero") then
      call buoyancy_forcing_zero(state, fluxes, day_center, dt, G, CS)
    elseif (trim(CS%buoy_config) == "linear") then
      call buoyancy_forcing_linear(state, fluxes, day_center, dt, G, CS)
    elseif (trim(CS%buoy_config) == "MESO") then
      call MOM_error(FATAL, "MESO forcing is not available with the ice-shelf"//&
               "version of MOM_surface_forcing.")
!      call MESO_buoyancy_forcing(state, fluxes, day_center, dt, G, CS%MESO_forcing_CSp)
    elseif (trim(CS%buoy_config) == "USER") then
      call USER_buoyancy_forcing(state, fluxes, day_center, dt, G, CS%user_forcing_CSp)
    elseif (trim(CS%buoy_config) == "NONE") then
      call MOM_mesg("MOM_surface_forcing: buoyancy forcing has been set to omitted.")
    elseif (CS%variable_buoyforce .and. .not.CS%first_call_set_forcing) then
      call MOM_error(FATAL, &
       "MOM_surface_forcing: Variable buoy defined with no buoy config.")
    else
       call MOM_error(FATAL, &
       "MOM_surface_forcing: Unrecognized buoy config "//trim(CS%buoy_config))
    endif
  endif
  if (associated(CS%tracer_flow_CSp)) then
    if (.not.associated(fluxes%tr_fluxes)) allocate(fluxes%tr_fluxes)
    call call_tracer_set_forcing(state, fluxes, day_start, day_interval, G, CS%tracer_flow_CSp)
  endif

  ! Allow for user-written code to alter the fluxes after all the above
  call user_alter_forcing(state, fluxes, day_center, G, CS%urf_CS)

  CS%first_call_set_forcing = .false.

  call cpu_clock_end(id_clock_forcing)
end subroutine set_forcing

subroutine buoyancy_forcing_allocate(fluxes, G, CS)
  type(forcing),         intent(inout) :: fluxes
  type(ocean_grid_type), intent(in)    :: G
  type(surface_forcing_CS), pointer    :: CS

!  This subroutine allocates arrays for buoyancy forcing.

! Arguments:
!  (inout)   fluxes  = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      G       = ocean grid structure
!  (in)      CS      = pointer to control struct returned by previous surface_forcing_init call

  integer :: isd, ied, jsd, jed
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed


  ! this array is zero for all options
  if (associated(fluxes%p_surf)) then
    fluxes%p_surf(:,:) = 0.0
  endif

  if ( CS%use_temperature ) then

    ! specify surface freshwater forcing by setting the following (kg/(m^2 * s))
    ! with convention that positive values for water entering ocean.
    if (.not.associated(fluxes%evap)) then
      allocate(fluxes%evap(isd:ied,jsd:jed))
      fluxes%evap(:,:) = 0.0
    endif
    if (.not.associated(fluxes%lprec)) then
      allocate(fluxes%lprec(isd:ied,jsd:jed))
      fluxes%lprec(:,:) = 0.0
    endif
    if (.not.associated(fluxes%fprec)) then
      allocate(fluxes%fprec(isd:ied,jsd:jed))
      fluxes%fprec(:,:) = 0.0
    endif
    if (.not.associated(fluxes%vprec)) then
      allocate(fluxes%vprec(isd:ied,jsd:jed))
      fluxes%vprec(:,:) = 0.0
    endif
    if (.not.associated(fluxes%lrunoff)) then
      allocate(fluxes%lrunoff(isd:ied,jsd:jed))
      fluxes%lrunoff(:,:) = 0.0
    endif
    if (.not.associated(fluxes%frunoff)) then
      allocate(fluxes%frunoff(isd:ied,jsd:jed))
      fluxes%frunoff(:,:) = 0.0
    endif
    if (.not.associated(fluxes%seaice_melt)) then
      allocate(fluxes%seaice_melt(isd:ied,jsd:jed))
      fluxes%seaice_melt(:,:) = 0.0
    endif

    ! specify surface heat fluxes by setting the following (Watts/m^2)
    ! with convention that positive values for heat fluxes into the ocean.
    if (.not.associated(fluxes%sw)) then
      allocate(fluxes%sw(isd:ied,jsd:jed))
      fluxes%sw(:,:) = 0.0
    endif
    if (.not.associated(fluxes%lw)) then
      allocate(fluxes%lw(isd:ied,jsd:jed))
      fluxes%lw(:,:) = 0.0
    endif
    if (.not.associated(fluxes%latent)) then
      allocate(fluxes%latent(isd:ied,jsd:jed))
      fluxes%latent(:,:) = 0.0
    endif
    if (.not.associated(fluxes%latent_evap_diag)) then
      allocate(fluxes%latent_evap_diag(isd:ied,jsd:jed))
      fluxes%latent_evap_diag(:,:) = 0.0
    endif
    if (.not.associated(fluxes%latent_fprec_diag)) then
      allocate(fluxes%latent_fprec_diag(isd:ied,jsd:jed))
      fluxes%latent_fprec_diag(:,:) = 0.0
    endif
    if (.not.associated(fluxes%latent_frunoff_diag)) then
      allocate(fluxes%latent_frunoff_diag(isd:ied,jsd:jed))
      fluxes%latent_frunoff_diag(:,:) = 0.0
    endif
    if (.not.associated(fluxes%sens)) then
      allocate(fluxes%sens(isd:ied,jsd:jed))
      fluxes%sens(:,:) = 0.0
    endif
    if (.not.associated(fluxes%heat_content_cond)) then
      allocate(fluxes%heat_content_cond(isd:ied,jsd:jed))
      fluxes%heat_content_cond(:,:) = 0.0
    endif
    if (.not.associated(fluxes%heat_content_lprec)) then
      allocate(fluxes%heat_content_lprec(isd:ied,jsd:jed))
      fluxes%heat_content_lprec(:,:) = 0.0
    endif
    if (.not.associated(fluxes%heat_content_fprec)) then
      allocate(fluxes%heat_content_fprec(isd:ied,jsd:jed))
      fluxes%heat_content_fprec(:,:) = 0.0
    endif
    if (.not.associated(fluxes%heat_content_vprec)) then
      allocate(fluxes%heat_content_vprec(isd:ied,jsd:jed))
      fluxes%heat_content_vprec(:,:) = 0.0
    endif
    if (.not.associated(fluxes%heat_content_lrunoff)) then
      allocate(fluxes%heat_content_lrunoff(isd:ied,jsd:jed))
      fluxes%heat_content_lrunoff(:,:) = 0.0
    endif
    if (.not.associated(fluxes%heat_content_frunoff)) then
      allocate(fluxes%heat_content_frunoff(isd:ied,jsd:jed))
      fluxes%heat_content_frunoff(:,:) = 0.0
    endif
    if (.not.associated(fluxes%heat_content_massout)) then
      allocate(fluxes%heat_content_massout(isd:ied,jsd:jed))
      fluxes%heat_content_massout(:,:) = 0.0
    endif

    ! surface restoring fields
    if (CS%restorebuoy) then
      if (.not.associated(CS%T_Restore)) then
        allocate(CS%T_Restore(isd:ied,jsd:jed))
        CS%T_Restore(:,:) = 0.0
      endif
      if (.not.associated(fluxes%heat_restore)) then
        allocate(fluxes%heat_restore(isd:ied,jsd:jed))
        fluxes%heat_restore(:,:) = 0.0
      endif
      if (.not.associated(CS%S_Restore)) then
        allocate(CS%S_Restore(isd:ied,jsd:jed))
        CS%S_Restore(:,:) = 0.0
      endif
    endif

  else ! CS%use_temperature false.

    if (.not.associated(fluxes%buoy)) then
      allocate(fluxes%buoy(isd:ied,jsd:jed))
      fluxes%buoy(:,:) = 0.0
    endif
    if (CS%restorebuoy .and. .not.associated(CS%Dens_Restore)) then
      allocate(CS%Dens_Restore(isd:ied,jsd:jed))
      CS%Dens_Restore(:,:) = 0.0
    endif

  endif  ! endif for  CS%use_temperature

end subroutine buoyancy_forcing_allocate


subroutine wind_forcing_zero(state, fluxes, day, G, CS)
  type(surface),            intent(inout) :: state
  type(forcing),            intent(inout) :: fluxes
  type(time_type),          intent(in)    :: day
  type(ocean_grid_type),    intent(in)    :: G
  type(surface_forcing_CS), pointer       :: CS

! subroutine sets the surface wind stresses to zero

! Arguments:
!            state        = structure describing ocean surface state
!  (out)     fluxes       = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day          = time of the fluxes
!  (in)      G            = ocean grid structure
!  (in)      CS           = pointer to control struct returned by previous surface_forcing_init call

  real :: PI
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  call callTree_enter("wind_forcing_zero, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call allocate_forcing_type(G, fluxes, stress=.true., ustar=.true.)

  !set steady surface wind stresses, in units of Pa.
  PI = 4.0*atan(1.0)

  do j=js,je ; do I=Isq,Ieq
    fluxes%taux(I,j) = 0.0
  enddo ; enddo

  do J=Jsq,Jeq ; do i=is,ie
    fluxes%tauy(i,J) = 0.0
  enddo ; enddo

  if (CS%read_gust_2d) then
    if (associated(fluxes%ustar)) then ; do j=js,je ; do i=is,ie
      fluxes%ustar(i,j) = sqrt(CS%gust(i,j)/CS%Rho0)
    enddo ; enddo ; endif
  else
    if (associated(fluxes%ustar)) then ; do j=js,je ; do i=is,ie
      fluxes%ustar(i,j) = sqrt(CS%gust_const/CS%Rho0)
    enddo ; enddo ; endif
  endif

  call callTree_leave("wind_forcing_zero")
end subroutine wind_forcing_zero


subroutine wind_forcing_2gyre(state, fluxes, day, G, CS)
  type(surface),            intent(inout) :: state
  type(forcing),            intent(inout) :: fluxes
  type(time_type),          intent(in)    :: day
  type(ocean_grid_type),    intent(in)    :: G
  type(surface_forcing_CS), pointer       :: CS

! This subroutine sets the surface wind stresses according to double gyre.

! Arguments:
!            state   = structure describing ocean surface state
!  (out)     fluxes  = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day     = time of the fluxes
!  (in)      G       = ocean grid structure
!  (in)      CS      = pointer to control struct returned by previous surface_forcing_init call

  real :: PI
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  call callTree_enter("wind_forcing_2gyre, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call allocate_forcing_type(G, fluxes, stress=.true., ustar=.true.)

  !set the steady surface wind stresses, in units of Pa.
  PI = 4.0*atan(1.0)

  do j=js,je ; do I=Isq,Ieq
    fluxes%taux(I,j) = 0.1*(1.0 - cos(2.0*PI*(G%geoLatCu(I,j)-CS%South_lat) / &
                                      CS%len_lat))
  enddo ; enddo

  do J=Jsq,Jeq ; do i=is,ie
    fluxes%tauy(i,J) = 0.0
  enddo ; enddo

  call callTree_leave("wind_forcing_2gyre")
end subroutine wind_forcing_2gyre


subroutine wind_forcing_1gyre(state, fluxes, day, G, CS)
  type(surface),            intent(inout) :: state
  type(forcing),            intent(inout) :: fluxes
  type(time_type),          intent(in)    :: day
  type(ocean_grid_type),    intent(in)    :: G
  type(surface_forcing_CS), pointer       :: CS

! This subroutine sets the surface wind stresses according to single gyre.

! Arguments:
!            state   = structure describing ocean surface state
!  (out)     fluxes  = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day     = time of the fluxes
!  (in)      G       = ocean grid structure
!  (in)      CS      = pointer to control struct returned by previous surface_forcing_init call

  real :: PI
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  call callTree_enter("wind_forcing_1gyre, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call allocate_forcing_type(G, fluxes, stress=.true., ustar=.true.)

  ! set the steady surface wind stresses, in units of Pa.
  PI = 4.0*atan(1.0)

  do j=js,je ; do I=Isq,Ieq
    fluxes%taux(I,j) =-0.2*cos(PI*(G%geoLatCu(I,j)-CS%South_lat)/CS%len_lat)
  enddo ; enddo

  do J=Jsq,Jeq ; do i=is,ie
    fluxes%tauy(i,J) = 0.0
  enddo ; enddo

  call callTree_leave("wind_forcing_1gyre")
end subroutine wind_forcing_1gyre


subroutine wind_forcing_gyres(state, fluxes, day, G, CS)
  type(surface),            intent(inout) :: state
  type(forcing),            intent(inout) :: fluxes
  type(time_type),          intent(in)    :: day
  type(ocean_grid_type),    intent(in)    :: G
  type(surface_forcing_CS), pointer       :: CS

! This subroutine sets the surface wind stresses according to gyres.

! Arguments:
!            state  = structure describing ocean surface state
!  (out)     fluxes = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day    = time of the fluxes
!  (in)      G      = ocean grid structure
!  (in)      CS     = pointer to control struct returned by previous surface_forcing_init call

  real :: PI, y
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  call callTree_enter("wind_forcing_gyres, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call allocate_forcing_type(G, fluxes, stress=.true., ustar=.true.)

  ! steady surface wind stresses (Pa)
  PI = 4.0*atan(1.0)

  do j=jsd,jed ; do I=IsdB,IedB
    y = (G%geoLatCu(I,j)-CS%South_lat)/CS%len_lat
    fluxes%taux(I,j) = CS%gyres_taux_const +                            &
             (   CS%gyres_taux_sin_amp*sin(CS%gyres_taux_n_pis*PI*y)    &
               + CS%gyres_taux_cos_amp*cos(CS%gyres_taux_n_pis*PI*y) )
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    fluxes%tauy(i,J) = 0.0
  enddo ; enddo

  ! set the friction velocity
  do j=js,je ; do i=is,ie
    fluxes%ustar(i,j) = sqrt(sqrt(0.5*(fluxes%tauy(i,j-1)*fluxes%tauy(i,j-1) + &
      fluxes%tauy(i,j)*fluxes%tauy(i,j) + fluxes%taux(i-1,j)*fluxes%taux(i-1,j) + &
      fluxes%taux(i,j)*fluxes%taux(i,j)))/CS%Rho0 + (CS%gust_const/CS%Rho0))
  enddo ; enddo

  call callTree_leave("wind_forcing_gyres")
end subroutine wind_forcing_gyres


subroutine wind_forcing_from_file(state, fluxes, day, G, CS)
  type(surface),            intent(inout) :: state
  type(forcing),            intent(inout) :: fluxes
  type(time_type),          intent(in)    :: day
  type(ocean_grid_type),    intent(inout) :: G
  type(surface_forcing_CS), pointer       :: CS

! This subroutine sets the surface wind stresses.

! Arguments:
!            state  = structure describing ocean surface state
!  (out)     fluxes = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day    = time of the fluxes
!  (in)      G      = ocean grid structure
!  (in)      CS     = pointer to control struct returned by previous surface_forcing_init call

  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: time_lev             ! With fields from a file, this must
                                  ! be reset, depending on the time.
  character(len=200) :: filename  ! The name of the input file.
  real :: temp_x(SZI_(G),SZJ_(G)) ! Pseudo-zonal and psuedo-meridional
  real :: temp_y(SZI_(G),SZJ_(G)) ! wind stresses at h-points, in Pa.
  integer :: days, seconds

  call callTree_enter("wind_forcing_from_file, MOM_surface_forcing.F90")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call allocate_forcing_type(G, fluxes, stress=.true., ustar=.true.)
  call get_time(day,seconds,days)
  time_lev = days - 365*floor(real(days) / 365.0) +1

  if (time_lev /= CS%wind_last_lev_read) then
    filename = trim(CS%inputdir) // trim(CS%wind_file)
!    if (is_root_pe()) &
!      write(*,'("Wind_forcing Reading time level ",I," last was ",I,".")')&
!           time_lev-1,CS%wind_last_lev_read-1
    select case ( uppercase(CS%wind_stagger(1:1)) )
    case ("A")
      temp_x(:,:) = 0.0 ; temp_y(:,:) = 0.0
      call read_data(filename,CS%stress_x_var,temp_x(:,:), &
                     domain=G%Domain%mpp_domain,timelevel=time_lev)
      call read_data(filename,CS%stress_y_var,temp_y(:,:), &
                     domain=G%Domain%mpp_domain,timelevel=time_lev)

      call pass_vector(temp_x, temp_y, G%Domain, To_All, AGRID)
      do j=js,je ; do I=Isq,Ieq
        fluxes%taux(I,j) = 0.5 * CS%wind_scale * (temp_x(i,j) + temp_x(i+1,j))
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        fluxes%tauy(i,J) = 0.5 * CS%wind_scale * (temp_y(i,j) + temp_y(i,j+1))
      enddo ; enddo

      if (CS%read_gust_2d) then
        do j=js,je ; do i=is,ie
          fluxes%ustar(i,j) = sqrt((sqrt(temp_x(i,j)*temp_x(i,j) + &
              temp_y(i,j)*temp_y(i,j)) + CS%gust(i,j)) / CS%Rho0)
        enddo ; enddo
      else
        do j=js,je ; do i=is,ie
          fluxes%ustar(i,j) = sqrt(sqrt(temp_x(i,j)*temp_x(i,j) + &
              temp_y(i,j)*temp_y(i,j))/CS%Rho0 + (CS%gust_const/CS%Rho0))
        enddo ; enddo
      endif
    case ("C")
      call read_data(filename,CS%stress_x_var,fluxes%taux(:,:), &
                     domain=G%Domain%mpp_domain,timelevel=time_lev)
      call read_data(filename,CS%stress_y_var,fluxes%tauy(:,:), &
                     domain=G%Domain%mpp_domain,timelevel=time_lev)
      if (CS%wind_scale /= 1.0) then
        do j=js,je ; do I=Isq,Ieq
          fluxes%taux(I,j) = CS%wind_scale * fluxes%taux(I,j)
        enddo ; enddo
        do J=Jsq,Jeq ; do i=is,ie
          fluxes%tauy(i,J) = CS%wind_scale * fluxes%tauy(i,J)
        enddo ; enddo
      endif

      call pass_vector(fluxes%taux, fluxes%tauy, G%Domain, To_All)
      if (CS%read_gust_2d) then
        do j=js, je ; do i=is, ie
          fluxes%ustar(i,j) = sqrt((sqrt(0.5*((fluxes%tauy(i,j-1)**2 + &
            fluxes%tauy(i,j)**2) + (fluxes%taux(i-1,j)**2 + &
            fluxes%taux(i,j)**2))) + CS%gust(i,j)) / CS%Rho0 )
        enddo ; enddo
      else
        do j=js, je ; do i=is, ie
          fluxes%ustar(i,j) = sqrt(sqrt(0.5*((fluxes%tauy(i,j-1)**2 + &
            fluxes%tauy(i,j)**2) + (fluxes%taux(i-1,j)**2 + &
            fluxes%taux(i,j)**2)))/CS%Rho0 + (CS%gust_const/CS%Rho0))
        enddo ; enddo
      endif
    case default
      call MOM_error(FATAL, "wind_forcing_from_file: Unrecognized stagger "//&
                      trim(CS%wind_stagger)//" is not 'A' or 'C'.")
    end select
    CS%wind_last_lev_read = time_lev
  endif ! time_lev /= CS%wind_last_lev_read

  call callTree_leave("wind_forcing_from_file")
end subroutine wind_forcing_from_file


subroutine buoyancy_forcing_from_files(state, fluxes, day, dt, G, CS)
  type(surface),         intent(inout) :: state
  type(forcing),         intent(inout) :: fluxes
  type(time_type),       intent(in)    :: day
  real,                  intent(in)    :: dt
  type(ocean_grid_type), intent(inout) :: G
  type(surface_forcing_CS), pointer    :: CS

!  This subroutine specifies the current surface fluxes of buoyancy
!  temperature and fresh water.  It may also be modified to add
!  surface fluxes of user provided tracers.
!
! Arguments:
!            state   = structure describing ocean surface state
!  (out)     fluxes  = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day     = time of the fluxes
!  (in)      dt      = amount of time over which the fluxes apply
!  (in)      G       = ocean grid structure
!  (in)      CS      = pointer to control struct returned by previous surface_forcing_init call

  real :: rhoXcp ! mean density times the heat capacity, in J m-3 K-1.
  real :: Irho0  ! inverse Boussinesq reference density, in m3 kg-1.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed

  integer :: time_lev           ! With fields from a file, this must
                                ! be reset, depending on the time.
  integer :: time_lev_monthly   ! With fields from a file, this must
                                ! be reset, depending on the time.
  integer :: days, seconds
  real, dimension(SZI_(G),SZJ_(G)) :: &
    temp, &       ! A 2-d temporary work array with various units.
    SST_anom, &   ! Instantaneous sea surface temperature anomalies from a
                  ! target (observed) value, in deg C.
    SSS_anom, &   ! Instantaneous sea surface salinity anomalies from a target
                  ! (observed) value, in g kg-1.
    SSS_mean      ! A (mean?) salinity about which to normalize local salinity
                  ! anomalies when calculating restorative precipitation
                  ! anomalies, in g kg-1.

  call callTree_enter("buoyancy_forcing_from_files, MOM_surface_forcing.F90")

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ! allocate and initialize arrays
  call buoyancy_forcing_allocate(fluxes, G, CS)

  if (CS%use_temperature) rhoXcp = CS%Rho0 * fluxes%C_p
  Irho0 = 1.0/CS%Rho0

  ! Read the file containing the buoyancy forcing.
  call get_time(day,seconds,days)

   time_lev = days - 365*floor(real(days) / 365.0)

  if (time_lev < 31) then ; time_lev_monthly = 0
  else if (time_lev < 59)  then ; time_lev_monthly = 1
  else if (time_lev < 90)  then ; time_lev_monthly = 2
  else if (time_lev < 120) then ; time_lev_monthly = 3
  else if (time_lev < 151) then ; time_lev_monthly = 4
  else if (time_lev < 181) then ; time_lev_monthly = 5
  else if (time_lev < 212) then ; time_lev_monthly = 6
  else if (time_lev < 243) then ; time_lev_monthly = 7
  else if (time_lev < 273) then ; time_lev_monthly = 8
  else if (time_lev < 304) then ; time_lev_monthly = 9
  else if (time_lev < 334) then ; time_lev_monthly = 10
  else ; time_lev_monthly = 11
  endif

  time_lev = time_lev+1
  time_lev_monthly = time_lev_monthly+1

  if (time_lev /= CS%buoy_last_lev_read) then

!   if (is_root_pe()) &
!     write(*,'("buoyancy_forcing : Reading time level ",I3,", last was ",I3,".")')&
!          time_lev,CS%buoy_last_lev_read


    call read_data(trim(CS%inputdir)//trim(CS%longwavedown_file), "lwdn_sfc", &
             fluxes%LW(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev)
    call read_data(trim(CS%inputdir)//trim(CS%longwaveup_file), "lwup_sfc", &
             temp(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev)
    do j=js,je ; do i=is,ie ; fluxes%LW(i,j) = fluxes%LW(i,j) - temp(i,j) ; enddo ; enddo

    call read_data(trim(CS%inputdir)//trim(CS%evaporation_file), "evap", &
             temp(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev)
    do j=js,je ; do i=is,ie
      fluxes%latent(i,j)           = -hlv*temp(i,j)
      fluxes%evap(i,j)             = -temp(i,j)
      fluxes%latent_evap_diag(i,j) = fluxes%latent(i,j)

    enddo ; enddo

    call read_data(trim(CS%inputdir)//trim(CS%sensibleheat_file), "shflx", &
             temp(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev)
    do j=js,je ; do i=is,ie ; fluxes%sens(i,j) = -temp(i,j) ; enddo ; enddo

    call read_data(trim(CS%inputdir)//trim(CS%shortwavedown_file), "swdn_sfc", &
             fluxes%sw(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev)
    call read_data(trim(CS%inputdir)//trim(CS%shortwaveup_file), "swup_sfc", &
             temp(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev)
    do j=js,je ; do i=is,ie
      fluxes%sw(i,j) = fluxes%sw(i,j) - temp(i,j)
    enddo ; enddo

    call read_data(trim(CS%inputdir)//trim(CS%snow_file), "snow", &
             fluxes%fprec(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev)
    call read_data(trim(CS%inputdir)//trim(CS%precip_file), "precip", &
             fluxes%lprec(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev)
    do j=js,je ; do i=is,ie
      fluxes%lprec(i,j) = fluxes%lprec(i,j) - fluxes%fprec(i,j)
    enddo ; enddo

    call read_data(trim(CS%inputdir)//trim(CS%freshdischarge_file), "disch_w", &
             temp(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev_monthly)
    do j=js,je ; do i=is,ie
      fluxes%lrunoff(i,j) = temp(i,j)*G%IareaT(i,j)
    enddo ; enddo
    call read_data(trim(CS%inputdir)//trim(CS%freshdischarge_file), "disch_s", &
              temp(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev_monthly)
    do j=js,je ; do i=is,ie
      fluxes%frunoff(i,j) = temp(i,j)*G%IareaT(i,j)
    enddo ; enddo

!     Read the SST and SSS fields for damping.
    if (CS%restorebuoy) then !### .or. associated(CS%ctrl_forcing_CSp)) then
      call read_data(trim(CS%inputdir)//trim(CS%SSTrestore_file), "TEMP", &
               CS%T_Restore(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev_monthly)
      call read_data(trim(CS%inputdir)//trim(CS%salinityrestore_file), "SALT", &
               CS%S_Restore(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev_monthly)
    endif
    CS%buoy_last_lev_read = time_lev

    ! mask out land points and compute heat content of water fluxes
    ! assume liquid precip enters ocean at SST
    ! assume frozen precip enters ocean at 0degC
    ! assume liquid runoff enters ocean at SST
    ! assume solid runoff (calving) enters ocean at 0degC
    do j=js,je ; do i=is,ie
      fluxes%evap(i,j)                 = fluxes%evap(i,j)             * G%mask2dT(i,j)
      fluxes%lprec(i,j)                = fluxes%lprec(i,j)            * G%mask2dT(i,j)
      fluxes%fprec(i,j)                = fluxes%fprec(i,j)            * G%mask2dT(i,j)
      fluxes%lrunoff(i,j)              = fluxes%lrunoff(i,j)          * G%mask2dT(i,j)
      fluxes%frunoff(i,j)              = fluxes%frunoff(i,j)          * G%mask2dT(i,j)
      fluxes%LW(i,j)                   = fluxes%LW(i,j)               * G%mask2dT(i,j)
      fluxes%sens(i,j)                 = fluxes%sens(i,j)             * G%mask2dT(i,j)
      fluxes%sw(i,j)                   = fluxes%sw(i,j)               * G%mask2dT(i,j)
      fluxes%latent(i,j)               = fluxes%latent(i,j)           * G%mask2dT(i,j)

      fluxes%heat_content_lrunoff(i,j) = fluxes%C_p*fluxes%lrunoff(i,j)*state%SST(i,j)
      fluxes%latent_evap_diag(i,j)     = fluxes%latent_evap_diag(i,j) * G%mask2dT(i,j)
      fluxes%latent_fprec_diag(i,j)    = -fluxes%fprec(i,j)*hlf
      fluxes%latent_frunoff_diag(i,j)  = -fluxes%frunoff(i,j)*hlf
    enddo ; enddo

  endif ! time_lev /= CS%buoy_last_lev_read

  if (CS%restorebuoy) then
    if (CS%use_temperature) then
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0) then
          fluxes%heat_restore(i,j) = G%mask2dT(i,j) * &
              ((CS%T_Restore(i,j) - state%SST(i,j)) * rhoXcp * CS%Flux_const)
          fluxes%vprec(i,j) = - (CS%Rho0*CS%Flux_const) * &
              (CS%S_Restore(i,j) - state%SSS(i,j)) / &
              (0.5*(state%SSS(i,j) + CS%S_Restore(i,j)))
        else
          fluxes%heat_restore(i,j) = 0.0
          fluxes%vprec(i,j) = 0.0
        endif
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0) then
          fluxes%buoy(i,j) = (CS%Dens_Restore(i,j) - state%sfc_density(i,j)) * &
                             (CS%G_Earth*CS%Flux_const/CS%Rho0)
        else
          fluxes%buoy(i,j) = 0.0
        endif
      enddo ; enddo
    endif
  else                                              ! not RESTOREBUOY
    if (.not.CS%use_temperature) then
      call MOM_error(FATAL, "buoyancy_forcing in MOM_surface_forcing: "// &
                     "The fluxes need to be defined without RESTOREBUOY.")
    endif
  endif                                             ! end RESTOREBUOY

!### if (associated(CS%ctrl_forcing_CSp)) then
!###   do j=js,je ; do i=is,ie
!###     SST_anom(i,j) = state%SST(i,j) - CS%T_Restore(i,j)
!###     SSS_anom(i,j) = state%SSS(i,j) - CS%S_Restore(i,j)
!###     SSS_mean(i,j) = 0.5*(state%SSS(i,j) + CS%S_Restore(i,j))
!###   enddo ; enddo
!###   call apply_ctrl_forcing(SST_anom, SSS_anom, SSS_mean, fluxes%heat_restore, &
!###                           fluxes%vprec, day, dt, G, CS%ctrl_forcing_CSp)
!### endif

  call callTree_leave("buoyancy_forcing_from_files")
end subroutine buoyancy_forcing_from_files


subroutine buoyancy_forcing_zero(state, fluxes, day, dt, G, CS)
  type(surface),         intent(inout) :: state
  type(forcing),         intent(inout) :: fluxes
  type(time_type),       intent(in)    :: day
  real,                  intent(in)    :: dt
  type(ocean_grid_type), intent(in)    :: G
  type(surface_forcing_CS), pointer    :: CS

!  This subroutine specifies the current surface fluxes of buoyancy
!  temperature and fresh water.  It may also be modified to add
!  surface fluxes of user provided tracers.
!  This case has zero surface buoyancy forcing.

! Arguments:
!  (inout)   state   = structure describing ocean surface state
!  (inout)   fluxes  = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day     = time of the fluxes
!  (in)      dt      = amount of time over which the fluxes apply
!  (in)      G       = ocean grid structure
!  (in)      CS      = pointer to control struct returned by previous surface_forcing_init call

  integer :: i, j, is, ie, js, je

  call callTree_enter("buoyancy_forcing_zero, MOM_surface_forcing.F90")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec


  ! allocate and initialize arrays
  call buoyancy_forcing_allocate(fluxes, G, CS)

  if (CS%use_temperature) then
    do j=js,je ; do i=is,ie
      fluxes%evap(i,j)                 = 0.0
      fluxes%lprec(i,j)                = 0.0
      fluxes%fprec(i,j)                = 0.0
      fluxes%lrunoff(i,j)              = 0.0
      fluxes%frunoff(i,j)              = 0.0
      fluxes%lw(i,j)                   = 0.0
      fluxes%latent(i,j)               = 0.0
      fluxes%sens(i,j)                 = 0.0
      fluxes%sw(i,j)                   = 0.0
      fluxes%heat_content_lrunoff(i,j) = 0.0
      fluxes%latent_evap_diag(i,j)     = 0.0
      fluxes%latent_fprec_diag(i,j)    = 0.0
      fluxes%latent_frunoff_diag(i,j)  = 0.0
    enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      fluxes%buoy(i,j) = 0.0
    enddo ; enddo
  endif

  call callTree_leave("buoyancy_forcing_zero")
end subroutine buoyancy_forcing_zero


subroutine buoyancy_forcing_linear(state, fluxes, day, dt, G, CS)
  type(surface),         intent(inout) :: state
  type(forcing),         intent(inout) :: fluxes
  type(time_type),       intent(in)    :: day
  real,                  intent(in)    :: dt
  type(ocean_grid_type), intent(in)    :: G
  type(surface_forcing_CS), pointer    :: CS

!    This subroutine specifies the current surface fluxes of buoyancy
!  temperature and fresh water.  It may also be modified to add
!  surface fluxes of user provided tracers.
!
! Arguments:
!  (inout)   state   = structure describing ocean surface state
!  (inout)   fluxes  = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day     = time of the fluxes
!  (in)      dt      = amount of time over which the fluxes apply
!  (in)      G       = ocean grid structure
!  (in)      CS      = pointer to control struct returned by previous surface_forcing_init call

  real :: y, T_restore, S_restore
  integer :: i, j, is, ie, js, je

  call callTree_enter("buoyancy_forcing_linear, MOM_surface_forcing.F90")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  ! allocate and initialize arrays
  call buoyancy_forcing_allocate(fluxes, G, CS)

  ! This case has no surface buoyancy forcing.
  if (CS%use_temperature) then
    do j=js,je ; do i=is,ie
      fluxes%evap(i,j)                 = 0.0
      fluxes%lprec(i,j)                = 0.0
      fluxes%fprec(i,j)                = 0.0
      fluxes%lrunoff(i,j)              = 0.0
      fluxes%frunoff(i,j)              = 0.0
      fluxes%lw(i,j)                   = 0.0
      fluxes%latent(i,j)               = 0.0
      fluxes%sens(i,j)                 = 0.0
      fluxes%sw(i,j)                   = 0.0
      fluxes%heat_content_lrunoff(i,j) = 0.0
      fluxes%latent_evap_diag(i,j)     = 0.0
      fluxes%latent_fprec_diag(i,j)    = 0.0
      fluxes%latent_frunoff_diag(i,j)  = 0.0
    enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      fluxes%buoy(i,j) = 0.0
    enddo ; enddo
  endif

  if (CS%restorebuoy) then
    if (CS%use_temperature) then
      do j=js,je ; do i=is,ie
        y = (G%geoLatCu(I,j)-CS%South_lat)/CS%len_lat
        T_restore = CS%T_south + (CS%T_north-CS%T_south)*y
        S_restore = CS%S_south + (CS%S_north-CS%S_south)*y
        if (G%mask2dT(i,j) > 0) then
          fluxes%heat_restore(i,j) = G%mask2dT(i,j) * &
              ((T_Restore - state%SST(i,j)) * ((CS%Rho0 * fluxes%C_p) * CS%Flux_const))
          fluxes%vprec(i,j) = - (CS%Rho0*CS%Flux_const) * &
              (S_Restore - state%SSS(i,j)) / &
              (0.5*(state%SSS(i,j) + S_Restore))
        else
          fluxes%heat_restore(i,j) = 0.0
          fluxes%vprec(i,j) = 0.0
        endif
      enddo ; enddo
    else
      call MOM_error(FATAL, "buoyancy_forcing_linear in MOM_surface_forcing: "// &
                     "RESTOREBUOY to linear not written yet.")
     !do j=js,je ; do i=is,ie
     !  if (G%mask2dT(i,j) > 0) then
     !    fluxes%buoy(i,j) = (CS%Dens_Restore(i,j) - state%sfc_density(i,j)) * &
     !                       (CS%G_Earth*CS%Flux_const/CS%Rho0)
     !  else
     !    fluxes%buoy(i,j) = 0.0
     !  endif
     !enddo ; enddo
    endif
  else                                              ! not RESTOREBUOY
    if (.not.CS%use_temperature) then
      call MOM_error(FATAL, "buoyancy_forcing_linear in MOM_surface_forcing: "// &
                     "The fluxes need to be defined without RESTOREBUOY.")
    endif
  endif                                             ! end RESTOREBUOY

  call callTree_leave("buoyancy_forcing_linear")
end subroutine buoyancy_forcing_linear


subroutine forcing_save_restart(CS, G, Time, directory, time_stamped, &
                                filename_suffix)
  type(surface_forcing_CS),   pointer       :: CS
  type(ocean_grid_type),      intent(inout) :: G
  type(time_type),            intent(in)    :: Time
  character(len=*),           intent(in)    :: directory
  logical,          optional, intent(in)    :: time_stamped
  character(len=*), optional, intent(in)    :: filename_suffix
! Arguments: CS - A pointer to the control structure returned by a previous
!                 call to surface_forcing_init.
!  (in)      G - The ocean's grid structure.
!  (in)      Time - The model time at this call.  This is needed for mpp_write calls.
!  (in, opt) directory - An optional directory into which to write these restart files.
!  (in, opt) time_stamped - If true, the restart file names include
!                           a unique time stamp.  The default is false.
!  (in, opt) filename_suffix - An optional suffix (e.g., a time-stamp) to append
!                              to the restart file names.

  if (.not.associated(CS)) return
  if (.not.associated(CS%restart_CSp)) return

  call save_restart(directory, Time, 1, G, CS%restart_CSp, time_stamped)

end subroutine forcing_save_restart

subroutine surface_forcing_init(Time, G, param_file, diag, CS, tracer_flow_CSp)
  type(time_type),           intent(in) :: Time
  type(ocean_grid_type),     intent(in) :: G
  type(param_file_type),     intent(in) :: param_file
  type(diag_ctrl), target,   intent(in) :: diag
  type(surface_forcing_CS),  pointer    :: CS
  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      tracer_flow_CSp - A pointer to the control structure of the tracer
!                              flow control module.
  type(directories)  :: dirs
  logical            :: new_sim
  type(time_type)    :: Time_frc
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_surface_forcing" ! This module's name.
  character(len=60)  :: axis_units
  character(len=200) :: filename, gust_file ! The name of the gustiness input file.

  if (associated(CS)) then
    call MOM_error(WARNING, "surface_forcing_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  id_clock_forcing=cpu_clock_id('(Ocean surface forcing)', grain=CLOCK_MODULE)
  call cpu_clock_begin(id_clock_forcing)

  CS%diag => diag
  if (associated(tracer_flow_CSp)) CS%tracer_flow_CSp => tracer_flow_CSp

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state \n"//&
                 "variables.", default=.true.)
  call get_param(param_file, mod, "INPUTDIR", CS%inputdir, &
                 "The directory in which all input files are found.", &
                 default=".")
  CS%inputdir = slasher(CS%inputdir)

  call get_param(param_file, mod, "ADIABATIC", CS%adiabatic, &
                 "There are no diapycnal mass fluxes if ADIABATIC is \n"//&
                 "true. This assumes that KD = KDML = 0.0 and that \n"//&
                 "there is no buoyancy forcing, but makes the model \n"//&
                 "faster by eliminating subroutine calls.", default=.false.)
  call get_param(param_file, mod, "VARIABLE_WINDS", CS%variable_winds, &
                 "If true, the winds vary in time after the initialization.", &
                 default=.true.)
  call get_param(param_file, mod, "VARIABLE_BUOYFORCE", CS%variable_buoyforce, &
                 "If true, the buoyancy forcing varies in time after the \n"//&
                 "initialization of the model.", default=.true.)

  call get_param(param_file, mod, "BUOY_CONFIG", CS%buoy_config, &
                 "The character string that indicates how buoyancy forcing \n"//&
                 "is specified. Valid options include (file), (zero), \n"//&
                 "(linear), (USER), and (NONE).", fail_if_missing=.true.)
  if (trim(CS%buoy_config) == "file") then
    call get_param(param_file, mod, "LONGWAVEDOWN_FILE", CS%longwavedown_file, &
                 "The file with the downward longwave heat flux, in \n"//&
                 "variable lwdn_sfc.", fail_if_missing=.true.)
    call get_param(param_file, mod, "LONGWAVEUP_FILE", CS%longwaveup_file, &
                 "The file with the upward longwave heat flux, in \n"//&
                 "variable lwup_sfc.", fail_if_missing=.true.)
    call get_param(param_file, mod, "EVAPORATION_FILE", CS%evaporation_file, &
                 "The file with the evaporative moisture flux, in \n"//&
                 "variable evap.", fail_if_missing=.true.)
    call get_param(param_file, mod, "SENSIBLEHEAT_FILE", CS%sensibleheat_file, &
                 "The file with the sensible heat flux, in \n"//&
                 "variable shflx.", fail_if_missing=.true.)
    call get_param(param_file, mod, "SHORTWAVEUP_FILE", CS%shortwaveup_file, &
                 "The file with the upward shortwave heat flux.", &
                 fail_if_missing=.true.)
    call get_param(param_file, mod, "SHORTWAVEDOWN_FILE", CS%shortwavedown_file, &
                 "The file with the downward shortwave heat flux.", &
                 fail_if_missing=.true.)
    call get_param(param_file, mod, "SNOW_FILE", CS%snow_file, &
                 "The file with the downward frozen precip flux, in \n"//&
                 "variable snow.", fail_if_missing=.true.)
    call get_param(param_file, mod, "PRECIP_FILE", CS%precip_file, &
                 "The file with the downward total precip flux, in \n"//&
                 "variable precip.", fail_if_missing=.true.)
    call get_param(param_file, mod, "FRESHDISCHARGE_FILE", CS%freshdischarge_file, &
                 "The file with the fresh and frozen runoff/calving fluxes, \n"//&
                 "invariables disch_w and disch_s.", fail_if_missing=.true.)
    call get_param(param_file, mod, "SSTRESTORE_FILE", CS%SSTrestore_file, &
                 "The file with the SST toward which to restore in \n"//&
                 "variable TEMP.", fail_if_missing=.true.)
    call get_param(param_file, mod, "SALINITYRESTORE_FILE", CS%salinityrestore_file, &
                 "The file with the surface salinity toward which to \n"//&
                 "restore in variable SALT.", fail_if_missing=.true.)
  endif
  call get_param(param_file, mod, "WIND_CONFIG", CS%wind_config, &
                 "The character string that indicates how wind forcing \n"//&
                 "is specified. Valid options include (file), (2gyre), \n"//&
                 "(1gyre), (gyres), (zero), and (USER).", fail_if_missing=.true.)
  if (trim(CS%wind_config) == "file") then
    call get_param(param_file, mod, "WIND_FILE", CS%wind_file, &
                 "The file in which the wind stresses are found in \n"//&
                 "variables STRESS_X and STRESS_Y.", fail_if_missing=.true.)
    call get_param(param_file, mod, "WINDSTRESS_X_VAR",CS%stress_x_var, &
                 "The name of the x-wind stress variable in WIND_FILE.", &
                 default="STRESS_X")
    call get_param(param_file, mod, "WINDSTRESS_Y_VAR", CS%stress_y_var, &
                 "The name of the y-wind stress variable in WIND_FILE.", &
                 default="STRESS_Y")
    call get_param(param_file, mod, "WINDSTRESS_STAGGER",CS%wind_stagger, &
                 "A character indicating how the wind stress components \n"//&
                 "are staggered in WIND_FILE.  This may be A or C for now.", &
                 default="A")
    call get_param(param_file, mod, "WINDSTRESS_SCALE", CS%wind_scale, &
                 "A value by which the wind stresses in WIND_FILE are rescaled.", &
                 default=1.0, units="nondim")
  endif
  if (trim(CS%wind_config) == "gyres") then
    call get_param(param_file, mod, "TAUX_CONST", CS%gyres_taux_const, &
                 "With the gyres wind_config, the constant offset in the \n"//&
                 "zonal wind stress profile: \n"//&
                 "  A in taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L).", &
                 units="Pa", default=0.0)
    call get_param(param_file, mod, "TAUX_SIN_AMP",CS%gyres_taux_sin_amp, &
                 "With the gyres wind_config, the sine amplitude in the \n"//&
                 "zonal wind stress profile: \n"//&
                 "  B in taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L).", &
                 units="Pa", default=0.0)
    call get_param(param_file, mod, "TAUX_COS_AMP",CS%gyres_taux_cos_amp, &
                 "With the gyres wind_config, the cosine amplitude in \n"//&
                 "the zonal wind stress profile: \n"//&
                 "  C in taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L).", &
                 units="Pa", default=0.0)
    call get_param(param_file, mod, "TAUX_N_PIS",CS%gyres_taux_n_pis, &
                 "With the gyres wind_config, the number of gyres in \n"//&
                 "the zonal wind stress profile: \n"//&
                 "  n in taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L).", &
                 units="nondim", default=0.0)
  endif
  call get_param(param_file, mod, "SOUTHLAT", CS%south_lat, &
                 "The southern latitude of the domain or the equivalent \n"//&
                 "starting value for the y-axis.", units=axis_units, default=0.)
  call get_param(param_file, mod, "LENLAT", CS%len_lat, &
                 "The latitudinal or y-direction length of the domain.", &
                 units=axis_units, fail_if_missing=.true.)
   call get_param(param_file, mod, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
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
    if (trim(CS%buoy_config) == "linear") then
      call get_param(param_file, mod, "SST_NORTH", CS%T_north, &
                 "With buoy_config linear, the sea surface temperature \n"//&
                 "at the northern end of the domain toward which to \n"//&
                 "to restore.", units="deg C", default=0.0)
      call get_param(param_file, mod, "SST_SOUTH", CS%T_south, &
                 "With buoy_config linear, the sea surface temperature \n"//&
                 "at the southern end of the domain toward which to \n"//&
                 "to restore.", units="deg C", default=0.0)
      call get_param(param_file, mod, "SSS_NORTH", CS%S_north, &
                 "With buoy_config linear, the sea surface salinity \n"//&
                 "at the northern end of the domain toward which to \n"//&
                 "to restore.", units="PSU", default=35.0)
      call get_param(param_file, mod, "SSS_SOUTH", CS%S_south, &
                 "With buoy_config linear, the sea surface salinity \n"//&
                 "at the southern end of the domain toward which to \n"//&
                 "to restore.", units="PSU", default=35.0)
    endif
  endif
  call get_param(param_file, mod, "G_EARTH", CS%G_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)

  call get_param(param_file, mod, "GUST_CONST", CS%gust_const, &
                 "The background gustiness in the winds.", units="Pa", &
                 default=0.02)
  call get_param(param_file, mod, "READ_GUST_2D", CS%read_gust_2d, &
                 "If true, use a 2-dimensional gustiness supplied from \n"//&
                 "an input file", default=.false.)
  if (CS%read_gust_2d) then
    call get_param(param_file, mod, "GUST_2D_FILE", gust_file, &
                 "The file in which the wind gustiness is found in \n"//&
                 "variable gustiness.", fail_if_missing=.true.)
    call safe_alloc_ptr(CS%gust,G%isd,G%ied,G%jsd,G%jed) ; CS%gust(:,:) = 0.0
    filename = trim(CS%inputdir) // trim(gust_file)
    call read_data(filename,'gustiness',CS%gust,domain=G%domain%mpp_domain, &
                   timelevel=1) ! units should be Pa
  endif
  call get_param(param_file, mod, "AXIS_UNITS", axis_units, default="degrees")

!  All parameter settings are now known.

  if (trim(CS%wind_config) == "USER" .or. trim(CS%buoy_config) == "USER" ) then
    call USER_surface_forcing_init(Time, G, param_file, diag, CS%user_forcing_CSp)
  elseif (trim(CS%wind_config) == "MESO" .or. trim(CS%buoy_config) == "MESO" ) then
    call MOM_error(FATAL, "MESO forcing is not available with the ice-shelf"//&
               "version of MOM_surface_forcing.")
!    call MESO_surface_forcing_init(Time, G, param_file, diag, CS%MESO_forcing_CSp)
  endif

  ! Set up any restart fields associated with the forcing.
  call restart_init(G, param_file, CS%restart_CSp, "MOM_forcing.res")
!###  call register_ctrl_forcing_restarts(G, param_file, CS%ctrl_forcing_CSp, &
!###                                      CS%restart_CSp)
  call restart_init_end(CS%restart_CSp)

  if (associated(CS%restart_CSp)) then
    call Get_MOM_Input(dirs=dirs)

    new_sim = .false.
    if ((dirs%input_filename(1:1) == 'n') .and. &
        (LEN_TRIM(dirs%input_filename) == 1)) new_sim = .true.
    if (.not.new_sim) then
      call restore_state(dirs%input_filename, dirs%restart_input_dir, Time_frc, &
                         G, CS%restart_CSp)
    endif
  endif

!###  call controlled_forcing_init(Time, G, param_file, diag, CS%ctrl_forcing_CSp)

  call user_revise_forcing_init(param_file, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)
end subroutine surface_forcing_init


subroutine surface_forcing_end(CS, fluxes)
  type(surface_forcing_CS), pointer       :: CS
  type(forcing), optional,  intent(inout) :: fluxes
! Arguments:  CS - A pointer to the control structure returned by a previous
!                  call to surface_forcing_init, it will be deallocated here.
!  (inout)    fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.

  if (present(fluxes)) call deallocate_forcing_type(fluxes)

!###  call controlled_forcing_end(CS%ctrl_forcing_CSp)

  if (associated(CS)) deallocate(CS)
  CS => NULL()

end subroutine surface_forcing_end

end module MOM_surface_forcing
