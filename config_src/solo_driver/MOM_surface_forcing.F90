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
!*  Edited by Stephen Griffies, June 2014                              *
!*                                                                     *
!*    This program contains the subroutines that calculate the         *
!*  surface wind stresses and fluxes of buoyancy or temp/saln and      *
!*  fresh water.  These subroutines are called every time step,        *
!*  even if the wind stresses or buoyancy fluxes are constant in time  *
!*  - in that case these routines return quickly without doing         *
!*  anything.  In addition, any I/O of forcing fields is controlled    *
!*  by surface_forcing_init, located in this file.                     *
!*                                                                     *
!*    set_forcing is a small entry subroutine for the subroutines in   *
!*  this file.  It provides the external access to these subroutines.  *
!*                                                                     *
!*    wind_forcing determines the wind stresses and places them into   *
!*  fluxes%taux and fluxes%tauy.  Often wind_forcing must be tailored  *
!*  for a particular application - either by specifying file and input *
!*  variable names or by providing appropriate internal expressions    *
!*  for the stresses within a modified version of USER_wind_forcing.   *
!*                                                                     *
!*    buoyancy_forcing determines the surface fluxes of heat, fresh    *
!*  water and salt, as appropriate.  A restoring boundary              *
!*  condition plus a specified flux from a file is implemented here,   *
!*  but a user-provided internal expression can be set by modifying    *
!*  and calling USER_buoyancy_forcing.                                 *
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
use MOM_constants,           only : hlv, hlf
use MOM_cpu_clock,           only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,           only : CLOCK_MODULE
use MOM_diag_mediator,       only : post_data, query_averaging_enabled
use MOM_diag_mediator,       only : diag_ctrl, safe_alloc_ptr
use MOM_domains,             only : pass_var, pass_vector, AGRID, To_South, To_West, To_All
use MOM_domains,             only : fill_symmetric_edges, CGRID_NE
use MOM_error_handler,       only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_error_handler,       only : callTree_enter, callTree_leave
use MOM_file_parser,         only : get_param, log_version, param_file_type
use MOM_string_functions,    only : uppercase
use MOM_forcing_type,        only : forcing, forcing_diags, register_forcing_type_diags, deallocate_forcing_type
use MOM_forcing_type,        only : allocate_forcing_type, deallocate_forcing_type
use MOM_grid,                only : ocean_grid_type
use MOM_get_input,           only : Get_MOM_Input, directories
use MOM_io,                  only : file_exists, read_data, slasher, num_timelevels
use MOM_io,                  only : EAST_FACE, NORTH_FACE
use MOM_restart,             only : register_restart_field, restart_init, MOM_restart_CS
use MOM_restart,             only : restart_init_end, save_restart, restore_state
use MOM_time_manager,        only : time_type, operator(+), operator(/), get_time, set_time
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_variables,           only : surface
use MESO_surface_forcing,    only : MESO_wind_forcing, MESO_buoyancy_forcing
use MESO_surface_forcing,    only : MESO_surface_forcing_init, MESO_surface_forcing_CS
use user_surface_forcing,    only : USER_wind_forcing, USER_buoyancy_forcing
use user_surface_forcing,    only : USER_surface_forcing_init, user_surface_forcing_CS
use user_revise_forcing,     only : user_alter_forcing, user_revise_forcing_init
use user_revise_forcing,     only : user_revise_forcing_CS
use SCM_idealized_hurricane, only : SCM_idealized_hurricane_wind_init
use SCM_idealized_hurricane, only : SCM_idealized_hurricane_wind_forcing
use SCM_idealized_hurricane, only : SCM_idealized_hurricane_CS
use SCM_CVmix_tests,         only : SCM_CVmix_tests_surface_forcing_init
use SCM_CVmix_tests,         only : SCM_CVmix_tests_wind_forcing
use SCM_CVmix_tests,         only : SCM_CVmix_tests_buoyancy_forcing
use SCM_CVmix_tests,         only : SCM_CVmix_tests_CS
use BFB_surface_forcing,    only : BFB_buoyancy_forcing
use BFB_surface_forcing,    only : BFB_surface_forcing_init, BFB_surface_forcing_CS

use data_override_mod, only : data_override_init, data_override

implicit none ; private

#include <MOM_memory.h>

public set_forcing
public surface_forcing_init
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
  real :: latent_heat_fusion    ! latent heat of fusion (J/kg)
  real :: latent_heat_vapor     ! latent heat of vaporization (J/kg)
  real :: tau_x0, tau_y0        ! Constant wind stresses used in the WIND_CONFIG="const" forcing

  real    :: gust_const                 ! constant unresolved background gustiness for ustar (Pa)
  logical :: read_gust_2d               ! if true, use 2-dimensional gustiness supplied from a file
  real, pointer :: gust(:,:) => NULL()  ! spatially varying unresolved background gustiness (Pa)
                                        ! gust is used when read_gust_2d is true.

  real, pointer :: T_Restore(:,:)    => NULL()  ! temperature to damp (restore) the SST to (deg C)
  real, pointer :: S_Restore(:,:)    => NULL()  ! salinity to damp (restore) the SSS (g/kg)
  real, pointer :: Dens_Restore(:,:) => NULL()  ! density to damp (restore) surface density (kg/m^3)

  integer :: buoy_last_lev_read = -1 ! The last time level read from buoyancy input files

  real :: gyres_taux_const, gyres_taux_sin_amp, gyres_taux_cos_amp, gyres_taux_n_pis
                             ! if WIND_CONFIG=='gyres' then use
                             ! = A, B, C and n respectively for
                             ! taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L)

  real :: T_north, T_south   ! target temperatures at north and south used in
                             ! buoyancy_forcing_linear
  real :: S_north, S_south   ! target salinity at north and south used in
                             ! buoyancy_forcing_linear

  logical :: first_call_set_forcing = .true.
  logical :: archaic_OMIP_file = .true.
  logical :: dataOverrideIsInitialized = .false.

  real :: wind_scale          ! value by which wind-stresses are scaled, ND.
  real :: constantHeatForcing ! value used for sensible heat flux when buoy_config="const"

  character(len=8)   :: wind_stagger
  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp => NULL()
!###  type(ctrl_forcing_CS), pointer :: ctrl_forcing_CSp => NULL()
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()

  type(diag_ctrl), pointer :: diag ! structure used to regulate timing of diagnostic output

  character(len=200) :: inputdir    ! directory where NetCDF input files are.
  character(len=200) :: wind_config ! indicator for wind forcing type (2gyre, USER, FILE..)
  character(len=200) :: wind_file   ! if wind_config is "file", file to use
  character(len=200) :: buoy_config ! indicator for buoyancy forcing type

  character(len=200) :: longwave_file     = ''
  character(len=200) :: shortwave_file    = ''
  character(len=200) :: evaporation_file  = ''
  character(len=200) :: sensibleheat_file = ''
  character(len=200) :: latentheat_file   = ''

  character(len=200) :: rain_file   = ''
  character(len=200) :: snow_file   = ''
  character(len=200) :: runoff_file = ''

  character(len=200) :: longwaveup_file  = ''
  character(len=200) :: shortwaveup_file = ''

  character(len=200) :: SSTrestore_file      = ''
  character(len=200) :: salinityrestore_file = ''

  character(len=80)  :: & ! Variable names in the input files
    stress_x_var    = '', stress_y_var    = '', ustar_var   = '', &
    LW_var          = '', SW_var          = '', latent_var  = '', sens_var    = '', evap_var = '', &
    rain_var        = '', snow_var        = '', lrunoff_var = '', frunoff_var = '', &
    SST_restore_var = '', SSS_restore_var = ''

  ! These variables give the number of time levels in the various forcing files.
  integer :: SW_nlev   = -1, LW_nlev   = -1, latent_nlev = -1, sens_nlev   = -1
  integer :: wind_nlev = -1, evap_nlev = -1, precip_nlev = -1, runoff_nlev = -1
  integer :: SST_nlev  = -1, SSS_nlev  = -1

  ! These variables give the last time level read for the various forcing files.
  integer :: wind_last_lev = -1
  integer :: SW_last_lev   = -1, LW_last_lev = -1, latent_last_lev = -1
  integer :: sens_last_lev = -1
  integer :: evap_last_lev = -1, precip_last_lev = -1, runoff_last_lev = -1
  integer :: SST_last_lev  = -1, SSS_last_lev = -1

  ! Diagnostics handles
  type(forcing_diags), public :: handles

  type(user_revise_forcing_CS),  pointer :: urf_CS => NULL()
  type(user_surface_forcing_CS), pointer :: user_forcing_CSp => NULL()
  type(BFB_surface_forcing_CS), pointer :: BFB_forcing_CSp => NULL()
  type(MESO_surface_forcing_CS), pointer :: MESO_forcing_CSp => NULL()
  type(SCM_idealized_hurricane_CS), pointer :: SCM_idealized_hurricane_CSp => NULL()
  type(SCM_CVmix_tests_CS),      pointer :: SCM_CVmix_tests_CSp => NULL()

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
  call callTree_enter("set_forcing, MOM_surface_forcing.F90")

  day_center = day_start + day_interval/2
  call get_time(day_interval, intdt)
  dt = real(intdt)

  ! calls to various wind options
  if (CS%variable_winds .or. CS%first_call_set_forcing) then
    if (CS%first_call_set_forcing) call allocate_forcing_type(G, fluxes, stress=.true., ustar=.true.)
    if (trim(CS%wind_config) == "file") then
      call wind_forcing_from_file(state, fluxes, day_center, G, CS)
    elseif (trim(CS%wind_config) == "data_override") then
      call wind_forcing_by_data_override(state, fluxes, day_center, G, CS)
    elseif (trim(CS%wind_config) == "2gyre") then
      call wind_forcing_2gyre(state, fluxes, day_center, G, CS)
    elseif (trim(CS%wind_config) == "1gyre") then
      call wind_forcing_1gyre(state, fluxes, day_center, G, CS)
    elseif (trim(CS%wind_config) == "gyres") then
      call wind_forcing_gyres(state, fluxes, day_center, G, CS)
    elseif (trim(CS%wind_config) == "zero") then
      call wind_forcing_const(state, fluxes, 0., 0., day_center, G, CS)
    elseif (trim(CS%wind_config) == "const") then
      call wind_forcing_const(state, fluxes, CS%tau_x0, CS%tau_y0, day_center, G, CS)
    elseif (trim(CS%wind_config) == "MESO") then
      call MESO_wind_forcing(state, fluxes, day_center, G, CS%MESO_forcing_CSp)
    elseif (trim(CS%wind_config) == "SCM_ideal_hurr") then
      call SCM_idealized_hurricane_wind_forcing(state, fluxes, day_center, G, CS%SCM_idealized_hurricane_CSp)
    elseif (trim(CS%wind_config) == "SCM_CVmix_tests") then
      call SCM_CVmix_tests_wind_forcing(state, fluxes, day_center, G, CS%SCM_CVmix_tests_CSp)
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

! ! calls to various buoyancy forcing options
! if ((CS%variable_buoyforce .or. CS%first_call_set_forcing) .and. &
!     (.not.CS%adiabatic)) then
!   if (trim(CS%buoy_config) == "file") then
!     if (CS%first_call_set_forcing) call buoyancy_forcing_allocate(fluxes, G, CS)
!   elseif (trim(CS%wind_config) == "MESO") then
!     call MESO_wind_forcing(state, fluxes, day_center, G, CS%MESO_forcing_CSp)
!   elseif (trim(CS%wind_config) == "SCM_ideal_hurr") then
!     call SCM_idealized_hurricane_wind_forcing(state, fluxes, day_center, G, CS%SCM_idealized_hurricane_CSp)
!   elseif (trim(CS%wind_config) == "USER") then
!     call USER_wind_forcing(state, fluxes, day_center, G, CS%user_forcing_CSp)
!   elseif (CS%variable_winds .and. .not.CS%first_call_set_forcing) then
!     call MOM_error(FATAL, &
!      "MOM_surface_forcing: Variable winds defined with no wind config")
!   else
!      call MOM_error(FATAL, &
!      "MOM_surface_forcing:Unrecognized wind config "//trim(CS%wind_config))
!   endif
! endif

  ! calls to various buoyancy forcing options
  if ((CS%variable_buoyforce .or. CS%first_call_set_forcing) .and. &
      (.not.CS%adiabatic)) then
    if (trim(CS%buoy_config) == "file") then
      if (CS%first_call_set_forcing) call buoyancy_forcing_allocate(fluxes, G, CS)
      call buoyancy_forcing_from_files(state, fluxes, day_center, dt, G, CS)
    elseif (trim(CS%buoy_config) == "data_override") then
      if (CS%first_call_set_forcing) call buoyancy_forcing_allocate(fluxes, G, CS)
      call buoyancy_forcing_from_data_override(state, fluxes, day_center, dt, G, CS)
    elseif (trim(CS%buoy_config) == "zero") then
      if (CS%first_call_set_forcing) call buoyancy_forcing_allocate(fluxes, G, CS)
      call buoyancy_forcing_zero(state, fluxes, day_center, dt, G, CS)
    elseif (trim(CS%buoy_config) == "const") then
      if (CS%first_call_set_forcing) call buoyancy_forcing_allocate(fluxes, G, CS)
      call buoyancy_forcing_const(state, fluxes, day_center, dt, G, CS)
    elseif (trim(CS%buoy_config) == "linear") then
      if (CS%first_call_set_forcing) call buoyancy_forcing_allocate(fluxes, G, CS)
      call buoyancy_forcing_linear(state, fluxes, day_center, dt, G, CS)
    elseif (trim(CS%buoy_config) == "MESO") then
      call MESO_buoyancy_forcing(state, fluxes, day_center, dt, G, CS%MESO_forcing_CSp)
    elseif (trim(CS%buoy_config) == "SCM_CVmix_tests") then
      call SCM_CVmix_tests_buoyancy_forcing(state, fluxes, day_center, G, CS%SCM_CVmix_tests_CSp)
    elseif (trim(CS%buoy_config) == "USER") then
      call USER_buoyancy_forcing(state, fluxes, day_center, dt, G, CS%user_forcing_CSp)
    elseif (trim(CS%buoy_config) == "BFB") then
      call BFB_buoyancy_forcing(state, fluxes, day_center, dt, G, CS%BFB_forcing_CSp)
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
  call callTree_leave("set_forcing")

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
    ! specify surface heat fluxes by setting the following (Watts/m^2)
    ! with convention that positive values for heat fluxes into the ocean.
    call allocate_forcing_type(G, fluxes, water=.true., heat=.true.)

    ! surface restoring fields
    if (CS%restorebuoy) then
      call safe_alloc_ptr(CS%T_Restore,isd,ied,jsd,jed)
      call safe_alloc_ptr(fluxes%heat_added,isd,ied,jsd,jed)
      call safe_alloc_ptr(CS%S_Restore,isd,ied,jsd,jed)
    endif

  else ! CS%use_temperature false.

    call safe_alloc_ptr(fluxes%buoy,isd,ied,jsd,jed)
    if (CS%restorebuoy) call safe_alloc_ptr(CS%Dens_Restore,isd,ied,jsd,jed)

  endif  ! endif for  CS%use_temperature

end subroutine buoyancy_forcing_allocate


subroutine wind_forcing_const(state, fluxes, tau_x0, tau_y0, day, G, CS)
  type(surface),            intent(inout) :: state
  type(forcing),            intent(inout) :: fluxes
  real,                     intent(in)    :: tau_x0
  real,                     intent(in)    :: tau_y0
  type(time_type),          intent(in)    :: day
  type(ocean_grid_type),    intent(in)    :: G
  type(surface_forcing_CS), pointer       :: CS

! subroutine sets the surface wind stresses to zero

! Arguments:
!            state        = structure describing ocean surface state
!  (out)     fluxes       = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day          = time of the fluxes
!  (in)      G            = ocean grid structure
!  (in)      CS           = pointer to control returned by previous surface_forcing_init call

  real :: mag_tau
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  call callTree_enter("wind_forcing_const, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  !set steady surface wind stresses, in units of Pa.
  mag_tau = sqrt( tau_x0**2 + tau_y0**2)

  do j=js,je ; do I=Isq,Ieq
    fluxes%taux(I,j) = tau_x0
  enddo ; enddo

  do J=Jsq,Jeq ; do i=is,ie
    fluxes%tauy(i,J) = tau_y0
  enddo ; enddo

  if (CS%read_gust_2d) then
    if (associated(fluxes%ustar)) then ; do j=js,je ; do i=is,ie
      fluxes%ustar(i,j) = sqrt( ( mag_tau + CS%gust(i,j) ) / CS%Rho0 )
    enddo ; enddo ; endif
  else
    if (associated(fluxes%ustar)) then ; do j=js,je ; do i=is,ie
      fluxes%ustar(i,j) = sqrt( ( mag_tau + CS%gust_const ) / CS%Rho0 )
    enddo ; enddo ; endif
  endif

  call callTree_leave("wind_forcing_const")
end subroutine wind_forcing_const


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

  character(len=200) :: filename  ! The name of the input file.
  real    :: temp_x(SZI_(G),SZJ_(G)) ! Pseudo-zonal and psuedo-meridional
  real    :: temp_y(SZI_(G),SZJ_(G)) ! wind stresses at h-points, in Pa.
  integer :: time_lev_daily          ! The time levels to read for fields with
  integer :: time_lev_monthly        ! daily and montly cycles.
  integer :: time_lev                ! The time level that is used for a field.
  integer :: days, seconds
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  logical :: read_Ustar

  call callTree_enter("wind_forcing_from_file, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call get_time(day,seconds,days)
  time_lev_daily = days - 365*floor(real(days) / 365.0)

  if (time_lev_daily < 31) then ; time_lev_monthly = 0
  else if (time_lev_daily < 59)  then ; time_lev_monthly = 1
  else if (time_lev_daily < 90)  then ; time_lev_monthly = 2
  else if (time_lev_daily < 120) then ; time_lev_monthly = 3
  else if (time_lev_daily < 151) then ; time_lev_monthly = 4
  else if (time_lev_daily < 181) then ; time_lev_monthly = 5
  else if (time_lev_daily < 212) then ; time_lev_monthly = 6
  else if (time_lev_daily < 243) then ; time_lev_monthly = 7
  else if (time_lev_daily < 273) then ; time_lev_monthly = 8
  else if (time_lev_daily < 304) then ; time_lev_monthly = 9
  else if (time_lev_daily < 334) then ; time_lev_monthly = 10
  else ; time_lev_monthly = 11
  endif

  time_lev_daily = time_lev_daily+1
  time_lev_monthly = time_lev_monthly+1

  select case (CS%wind_nlev)
    case (12) ; time_lev = time_lev_monthly
    case (365) ; time_lev = time_lev_daily
    case default ; time_lev = 1
  end select

  if (time_lev /= CS%wind_last_lev) then
    filename = trim(CS%wind_file)
    read_Ustar = (len_trim(CS%ustar_var) > 0)
!    if (is_root_pe()) &
!      write(*,'("Wind_forcing Reading time level ",I," last was ",I,".")')&
!           time_lev-1,CS%wind_last_lev-1
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

      if (.not.read_Ustar) then
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
      endif
    case ("C")
      if (G%symmetric) then
        if (.not.associated(G%Domain_aux)) call MOM_error(FATAL, &
          " wind_forcing_from_file with C-grid input and symmetric memory "//&
          " called with a non-associated auxiliary domain in the grid type.")
        !   Read the data as though symmetric memory were not being used, and
        ! then translate it appropriately.
        temp_x(:,:) = 0.0 ; temp_y(:,:) = 0.0
        call read_data(filename, CS%stress_x_var, temp_x(:,:), position=EAST_FACE, &
                       domain=G%Domain_aux%mpp_domain, timelevel=time_lev)
        call read_data(filename, CS%stress_y_var, temp_y(:,:), position=NORTH_FACE, &
                       domain=G%Domain_aux%mpp_domain, timelevel=time_lev)

        do j=js,je ; do i=is,ie
          fluxes%taux(I,j) = CS%wind_scale * temp_x(I,j)
          fluxes%tauy(i,J) = CS%wind_scale * temp_y(i,J)
        enddo ; enddo
        call fill_symmetric_edges(fluxes%taux, fluxes%tauy, G%Domain, stagger=CGRID_NE)
      else
        call read_data(filename, CS%stress_x_var, fluxes%taux(:,:), &
                       domain=G%Domain%mpp_domain, position=EAST_FACE, &
                       timelevel=time_lev)
        call read_data(filename, CS%stress_y_var, fluxes%tauy(:,:), &
                       domain=G%Domain%mpp_domain, position=NORTH_FACE, &
                       timelevel=time_lev)

        if (CS%wind_scale /= 1.0) then
          do j=js,je ; do I=Isq,Ieq
            fluxes%taux(I,j) = CS%wind_scale * fluxes%taux(I,j)
          enddo ; enddo
          do J=Jsq,Jeq ; do i=is,ie
            fluxes%tauy(i,J) = CS%wind_scale * fluxes%tauy(i,J)
          enddo ; enddo
        endif
      endif

      call pass_vector(fluxes%taux, fluxes%tauy, G%Domain, To_All)
      if (.not.read_Ustar) then
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
      endif
    case default
      call MOM_error(FATAL, "wind_forcing_from_file: Unrecognized stagger "//&
                      trim(CS%wind_stagger)//" is not 'A' or 'C'.")
    end select

    if (read_Ustar) then
      call read_data(filename, CS%Ustar_var, fluxes%ustar(:,:), &
                     domain=G%Domain%mpp_domain, timelevel=time_lev)
    endif

    CS%wind_last_lev = time_lev

  endif ! time_lev /= CS%wind_last_lev

  call callTree_leave("wind_forcing_from_file")
end subroutine wind_forcing_from_file


subroutine wind_forcing_by_data_override(state, fluxes, day, G, CS)
  type(surface),            intent(inout) :: state
  type(forcing),            intent(inout) :: fluxes
  type(time_type),          intent(in)    :: day
  type(ocean_grid_type),    intent(inout) :: G
  type(surface_forcing_CS), pointer       :: CS
! This subroutine sets the surface wind stresses

! Arguments:
!            state  = structure describing ocean surface state
!  (out)     fluxes = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day    = time of the fluxes
!  (in)      G      = ocean grid structure
!  (in)      CS     = pointer to control struct returned by previous surface_forcing_init call

  real :: temp_x(SZI_(G),SZJ_(G)) ! Pseudo-zonal and psuedo-meridional
  real :: temp_y(SZI_(G),SZJ_(G)) ! wind stresses at h-points, in Pa.
  integer :: i, j, is_in, ie_in, js_in, je_in
  logical :: read_uStar

  call callTree_enter("wind_forcing_by_data_override, MOM_surface_forcing.F90")

  if (.not.CS%dataOverrideIsInitialized) then
    call allocate_forcing_type(G, fluxes, stress=.true., ustar=.true.)
    call data_override_init(Ocean_domain_in=G%Domain%mpp_domain)
    CS%dataOverrideIsInitialized = .True.
  endif

  is_in = G%isc - G%isd + 1
  ie_in = G%iec - G%isd + 1
  js_in = G%jsc - G%jsd + 1
  je_in = G%jec - G%jsd + 1

  temp_x(:,:) = 0.0 ; temp_y(:,:) = 0.0
  call data_override('OCN', 'taux', temp_x, day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
  call data_override('OCN', 'tauy', temp_y, day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
  call pass_vector(temp_x, temp_y, G%Domain, To_All, AGRID)
  ! Ignore CS%wind_scale when using data_override ?????
  do j=G%jsc,G%jec ; do I=G%IscB,G%IecB
    fluxes%taux(I,j) = 0.5 * (temp_x(i,j) + temp_x(i+1,j))
  enddo ; enddo
  do J=G%JscB,G%JecB ; do i=G%isc,G%iec
    fluxes%tauy(i,J) = 0.5 * (temp_y(i,j) + temp_y(i,j+1))
  enddo ; enddo

  read_Ustar = (len_trim(CS%ustar_var) > 0) ! Need better control higher up ????
  if (read_Ustar) then
    call data_override('OCN', 'ustar', fluxes%ustar, day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
  else
    if (CS%read_gust_2d) then
      call data_override('OCN', 'gust', CS%gust, day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        fluxes%ustar(i,j) = sqrt((sqrt(temp_x(i,j)*temp_x(i,j) + &
            temp_y(i,j)*temp_y(i,j)) + CS%gust(i,j)) / CS%Rho0)
      enddo ; enddo
    else
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        fluxes%ustar(i,j) = sqrt(sqrt(temp_x(i,j)*temp_x(i,j) + &
            temp_y(i,j)*temp_y(i,j))/CS%Rho0 + (CS%gust_const/CS%Rho0))
      enddo ; enddo
    endif
  endif

  call pass_vector(fluxes%taux, fluxes%tauy, G%Domain, To_All)
! call pass_var(fluxes%ustar, G%Domain, To_All)     Not needed  ?????

  call callTree_leave("wind_forcing_by_data_override")
end subroutine wind_forcing_by_data_override


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
!  This case has surface buoyancy forcing from input files.

! Arguments:
!            state   = structure describing ocean surface state
!  (out)     fluxes  = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day     = time of the fluxes
!  (in)      dt      = amount of time over which the fluxes apply
!  (in)      G       = ocean grid structure
!  (in)      CS      = pointer to control struct returned by previous surface_forcing_init call

  real, dimension(SZI_(G),SZJ_(G)) :: &
    temp, &       ! A 2-d temporary work array with various units.
    SST_anom, &   ! Instantaneous sea surface temperature anomalies from a
                  ! target (observed) value, in deg C.
    SSS_anom, &   ! Instantaneous sea surface salinity anomalies from a target
                  ! (observed) value, in g kg-1.
    SSS_mean      ! A (mean?) salinity about which to normalize local salinity
                  ! anomalies when calculating restorative precipitation
                  ! anomalies, in g kg-1.

  real :: rhoXcp ! reference density times heat capacity (J/(m^3 * K))
  real :: Irho0  ! inverse of the Boussinesq reference density (m^3/kg)

  integer :: time_lev_daily     ! time levels to read for fields with daily cycle
  integer :: time_lev_monthly   ! time levels to read for fields with monthly cycle
  integer :: time_lev           ! time level that for a field

  integer :: days, seconds
  integer :: i, j, is, ie, js, je

  call callTree_enter("buoyancy_forcing_from_files, MOM_surface_forcing.F90")

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je = G%jec

  if (CS%use_temperature) rhoXcp = CS%Rho0 * fluxes%C_p
  Irho0 = 1.0/CS%Rho0

  ! Read the buoyancy forcing file
  call get_time(day,seconds,days)

  time_lev_daily = days - 365*floor(real(days) / 365.0)

  if (time_lev_daily < 31) then ; time_lev_monthly = 0
  else if (time_lev_daily < 59)  then ; time_lev_monthly = 1
  else if (time_lev_daily < 90)  then ; time_lev_monthly = 2
  else if (time_lev_daily < 120) then ; time_lev_monthly = 3
  else if (time_lev_daily < 151) then ; time_lev_monthly = 4
  else if (time_lev_daily < 181) then ; time_lev_monthly = 5
  else if (time_lev_daily < 212) then ; time_lev_monthly = 6
  else if (time_lev_daily < 243) then ; time_lev_monthly = 7
  else if (time_lev_daily < 273) then ; time_lev_monthly = 8
  else if (time_lev_daily < 304) then ; time_lev_monthly = 9
  else if (time_lev_daily < 334) then ; time_lev_monthly = 10
  else ; time_lev_monthly = 11
  endif

  time_lev_daily   = time_lev_daily  +1
  time_lev_monthly = time_lev_monthly+1

  if (time_lev_daily /= CS%buoy_last_lev_read) then

    ! longwave
    select case (CS%LW_nlev)
      case (12)    ; time_lev = time_lev_monthly
      case (365)   ; time_lev = time_lev_daily
      case default ; time_lev = 1
    end select
    call read_data(CS%longwave_file, CS%LW_var, fluxes%LW(:,:), &
                   domain=G%Domain%mpp_domain, timelevel=time_lev)
    if (CS%archaic_OMIP_file) then
      call read_data(CS%longwaveup_file, &
       "lwup_sfc",temp(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev)
      do j=js,je ; do i=is,ie ; fluxes%LW(i,j) = fluxes%LW(i,j) - temp(i,j) ; enddo ; enddo
    endif
    CS%LW_last_lev = time_lev

    ! evaporation
    select case (CS%evap_nlev)
      case (12)    ; time_lev = time_lev_monthly
      case (365)   ; time_lev = time_lev_daily
      case default ; time_lev = 1
    end select
    if (CS%archaic_OMIP_file) then
      call read_data(CS%evaporation_file, CS%evap_var, temp(:,:), &
                     domain=G%Domain%mpp_domain, timelevel=time_lev)
      do j=js,je ; do i=is,ie
        fluxes%latent(i,j)           = -CS%latent_heat_vapor*temp(i,j)
        fluxes%evap(i,j)             = -temp(i,j)
        fluxes%latent_evap_diag(i,j) = fluxes%latent(i,j)
      enddo ; enddo
    else
      call read_data(CS%evaporation_file, CS%evap_var, fluxes%evap(:,:), &
                     domain=G%Domain%mpp_domain, timelevel=time_lev)
    endif
    CS%evap_last_lev = time_lev

    select case (CS%latent_nlev)
      case (12)    ; time_lev = time_lev_monthly
      case (365)   ; time_lev = time_lev_daily
      case default ; time_lev = 1
    end select
    if (.not.CS%archaic_OMIP_file) then
      call read_data(CS%latentheat_file, CS%latent_var, fluxes%latent(:,:), &
                     domain=G%Domain%mpp_domain, timelevel=time_lev)
      do j=js,je ; do i=is,ie
        fluxes%latent_evap_diag(i,j) = fluxes%latent(i,j)
      enddo ; enddo
    endif
    CS%latent_last_lev = time_lev

    select case (CS%sens_nlev)
      case (12)    ; time_lev = time_lev_monthly
      case (365)   ; time_lev = time_lev_daily
      case default ; time_lev = 1
    end select
    if (CS%archaic_OMIP_file) then
      call read_data(CS%sensibleheat_file, CS%sens_var, temp(:,:), &
                     domain=G%Domain%mpp_domain, timelevel=time_lev)
      do j=js,je ; do i=is,ie ; fluxes%sens(i,j) = -temp(i,j) ; enddo ; enddo
    else
      call read_data(CS%sensibleheat_file, CS%sens_var, fluxes%sens(:,:), &
                     domain=G%Domain%mpp_domain, timelevel=time_lev)
    endif
    CS%sens_last_lev = time_lev

    select case (CS%SW_nlev)
      case (12)    ; time_lev = time_lev_monthly
      case (365)   ; time_lev = time_lev_daily
      case default ; time_lev = 1
    end select
    call read_data(CS%shortwave_file, CS%SW_var, fluxes%sw(:,:), &
             domain=G%Domain%mpp_domain, timelevel=time_lev)
    if (CS%archaic_OMIP_file) then
      call read_data(CS%shortwaveup_file, "swup_sfc", temp(:,:), &
               domain=G%Domain%mpp_domain, timelevel=time_lev)
      do j=js,je ; do i=is,ie
        fluxes%sw(i,j) = fluxes%sw(i,j) - temp(i,j)
      enddo ; enddo
    endif
    CS%SW_last_lev = time_lev

    select case (CS%precip_nlev)
      case (12)    ; time_lev = time_lev_monthly
      case (365)   ; time_lev = time_lev_daily
      case default ; time_lev = 1
    end select
    call read_data(CS%snow_file, CS%snow_var, &
             fluxes%fprec(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev)
    call read_data(CS%rain_file, CS%rain_var, &
             fluxes%lprec(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev)
    if (CS%archaic_OMIP_file) then
      do j=js,je ; do i=is,ie
        fluxes%lprec(i,j) = fluxes%lprec(i,j) - fluxes%fprec(i,j)
      enddo ; enddo
    endif
    CS%precip_last_lev = time_lev

    select case (CS%runoff_nlev)
      case (12)    ; time_lev = time_lev_monthly
      case (365)   ; time_lev = time_lev_daily
      case default ; time_lev = 1
    end select
    if (CS%archaic_OMIP_file) then
      call read_data(CS%runoff_file, CS%lrunoff_var, temp(:,:), &
                     domain=G%Domain%mpp_domain, timelevel=time_lev)
      do j=js,je ; do i=is,ie
        fluxes%lrunoff(i,j) = temp(i,j)*G%IareaT(i,j)
      enddo ; enddo
      call read_data(CS%runoff_file, CS%frunoff_var, temp(:,:), &
                     domain=G%Domain%mpp_domain, timelevel=time_lev)
      do j=js,je ; do i=is,ie
        fluxes%frunoff(i,j) = temp(i,j)*G%IareaT(i,j)
      enddo ; enddo
    else
      call read_data(CS%runoff_file, CS%lrunoff_var, fluxes%lrunoff(:,:), &
                     domain=G%Domain%mpp_domain, timelevel=time_lev)
      call read_data(CS%runoff_file, CS%frunoff_var, fluxes%frunoff(:,:), &
                     domain=G%Domain%mpp_domain, timelevel=time_lev)
    endif
    CS%runoff_last_lev = time_lev

!     Read the SST and SSS fields for damping.
    if (CS%restorebuoy) then !### .or. associated(CS%ctrl_forcing_CSp)) then
      select case (CS%SST_nlev)
        case (12)    ; time_lev = time_lev_monthly
        case (365)   ; time_lev = time_lev_daily
        case default ; time_lev = 1
      end select
      call read_data(CS%SSTrestore_file, CS%SST_restore_var, &
               CS%T_Restore(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev)
      CS%SST_last_lev = time_lev

      select case (CS%SSS_nlev)
        case (12)    ; time_lev = time_lev_monthly
        case (365)   ; time_lev = time_lev_daily
        case default ; time_lev = 1
      end select
      call read_data(CS%salinityrestore_file, CS%SSS_restore_var, &
               CS%S_Restore(:,:), domain=G%Domain%mpp_domain, timelevel=time_lev)
      CS%SSS_last_lev = time_lev
    endif
    CS%buoy_last_lev_read = time_lev_daily

    ! mask out land points and compute heat content of water fluxes
    ! assume liquid precip enters ocean at SST
    ! assume frozen precip enters ocean at 0degC
    ! assume liquid runoff enters ocean at SST
    ! assume solid runoff (calving) enters ocean at 0degC
    ! mass leaving the ocean has heat_content determined in MOM_diabatic_driver.F90
    do j=js,je ; do i=is,ie
      fluxes%evap(i,j)    = fluxes%evap(i,j)    * G%mask2dT(i,j)
      fluxes%lprec(i,j)   = fluxes%lprec(i,j)   * G%mask2dT(i,j)
      fluxes%fprec(i,j)   = fluxes%fprec(i,j)   * G%mask2dT(i,j)
      fluxes%lrunoff(i,j) = fluxes%lrunoff(i,j) * G%mask2dT(i,j)
      fluxes%frunoff(i,j) = fluxes%frunoff(i,j) * G%mask2dT(i,j)
      fluxes%LW(i,j)      = fluxes%LW(i,j)      * G%mask2dT(i,j)
      fluxes%sens(i,j)    = fluxes%sens(i,j)    * G%mask2dT(i,j)
      fluxes%sw(i,j)      = fluxes%sw(i,j)      * G%mask2dT(i,j)
      fluxes%latent(i,j)  = fluxes%latent(i,j)  * G%mask2dT(i,j)

      fluxes%latent_evap_diag(i,j)     = fluxes%latent_evap_diag(i,j) * G%mask2dT(i,j)
      fluxes%latent_fprec_diag(i,j)    = -fluxes%fprec(i,j)*CS%latent_heat_fusion
      fluxes%latent_frunoff_diag(i,j)  = -fluxes%frunoff(i,j)*CS%latent_heat_fusion
    enddo ; enddo

  endif ! time_lev /= CS%buoy_last_lev_read


  ! restoring surface boundary fluxes
  if (CS%restorebuoy) then

    if (CS%use_temperature) then
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0) then
          fluxes%heat_added(i,j) = G%mask2dT(i,j) * &
              ((CS%T_Restore(i,j) - state%SST(i,j)) * rhoXcp * CS%Flux_const)
          fluxes%vprec(i,j) = - (CS%Rho0*CS%Flux_const) * &
              (CS%S_Restore(i,j) - state%SSS(i,j)) / &
              (0.5*(state%SSS(i,j) + CS%S_Restore(i,j)))
        else
          fluxes%heat_added(i,j) = 0.0
          fluxes%vprec(i,j)        = 0.0
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
!###   call apply_ctrl_forcing(SST_anom, SSS_anom, SSS_mean, fluxes%heat_added, &
!###                           fluxes%vprec, day, dt, G, CS%ctrl_forcing_CSp)
!### endif

  call callTree_leave("buoyancy_forcing_from_files")
end subroutine buoyancy_forcing_from_files


subroutine buoyancy_forcing_from_data_override(state, fluxes, day, dt, G, CS)
  type(surface),         intent(inout) :: state
  type(forcing),         intent(inout) :: fluxes
  type(time_type),       intent(in)    :: day
  real,                  intent(in)    :: dt
  type(ocean_grid_type), intent(inout) :: G
  type(surface_forcing_CS), pointer    :: CS

!  This subroutine specifies the current surface fluxes of buoyancy
!  temperature and fresh water.  It may also be modified to add
!  surface fluxes of user provided tracers.
!  This case has surface buoyancy forcing from data over-ride.

! Arguments:
!            state   = structure describing ocean surface state
!  (out)     fluxes  = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day     = time of the fluxes
!  (in)      dt      = amount of time over which the fluxes apply
!  (in)      G       = ocean grid structure
!  (in)      CS      = pointer to control struct returned by previous surface_forcing_init call

  real, dimension(SZI_(G),SZJ_(G)) :: &
    temp, &       ! A 2-d temporary work array with various units.
    SST_anom, &   ! Instantaneous sea surface temperature anomalies from a
                  ! target (observed) value, in deg C.
    SSS_anom, &   ! Instantaneous sea surface salinity anomalies from a target
                  ! (observed) value, in g kg-1.
    SSS_mean      ! A (mean?) salinity about which to normalize local salinity
                  ! anomalies when calculating restorative precipitation
                  ! anomalies, in g kg-1.
  real :: rhoXcp ! The mean density times the heat capacity, in J m-3 K-1.
  real :: Irho0  ! The inverse of the Boussinesq density, in m3 kg-1.

  integer :: time_lev_daily     ! The time levels to read for fields with
  integer :: time_lev_monthly   ! daily and montly cycles.
  integer :: itime_lev           ! The time level that is used for a field.

  integer :: days, seconds
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  integer :: is_in, ie_in, js_in, je_in

  call callTree_enter("buoyancy_forcing_from_data_override, MOM_surface_forcing.F90")

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (CS%use_temperature) rhoXcp = CS%Rho0 * fluxes%C_p
  Irho0 = 1.0/CS%Rho0

  if (.not.CS%dataOverrideIsInitialized) then
    call data_override_init(Ocean_domain_in=G%Domain%mpp_domain)
    CS%dataOverrideIsInitialized = .True.
  endif

  is_in = G%isc - G%isd + 1
  ie_in = G%iec - G%isd + 1
  js_in = G%jsc - G%jsd + 1
  je_in = G%jec - G%jsd + 1

  call data_override('OCN', 'lw', fluxes%LW(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
  call data_override('OCN', 'evap', fluxes%evap(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)

  ! note the sign convention
  do j=js,je ; do i=is,ie
     fluxes%evap(i,j) = -fluxes%evap(i,j)  ! Normal convention is positive into the ocean
                                           ! but evap is normally a positive quantity in the files
     fluxes%latent(i,j)           = CS%latent_heat_vapor*fluxes%evap(i,j)
     fluxes%latent_evap_diag(i,j) = fluxes%latent(i,j)
  enddo; enddo

  call data_override('OCN', 'sens', fluxes%sens(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)

  ! note the sign convention
  do j=js,je ; do i=is,ie
     fluxes%sens(i,j) = -fluxes%sens(i,j)  ! Normal convention is positive into the ocean
                                           ! but sensible is normally a positive quantity in the files
  enddo; enddo

  call data_override('OCN', 'sw', fluxes%sw(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)

  call data_override('OCN', 'snow', fluxes%fprec(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)

  call data_override('OCN', 'rain', fluxes%lprec(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)

  call data_override('OCN', 'runoff', fluxes%lrunoff(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)

  call data_override('OCN', 'calving', fluxes%frunoff(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)

!     Read the SST and SSS fields for damping.
  if (CS%restorebuoy) then !### .or. associated(CS%ctrl_forcing_CSp)) then
     call data_override('OCN', 'SST_restore', CS%T_restore(:,:), day, &
          is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)

     call data_override('OCN', 'SSS_restore', CS%S_restore(:,:), day, &
          is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)

  endif

  ! restoring boundary fluxes
  if (CS%restorebuoy) then
    if (CS%use_temperature) then
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0) then
          fluxes%heat_added(i,j) = G%mask2dT(i,j) * &
              ((CS%T_Restore(i,j) - state%SST(i,j)) * rhoXcp * CS%Flux_const)
          fluxes%vprec(i,j) = - (CS%Rho0*CS%Flux_const) * &
              (CS%S_Restore(i,j) - state%SSS(i,j)) / &
              (0.5*(state%SSS(i,j) + CS%S_Restore(i,j)))
        else
          fluxes%heat_added(i,j) = 0.0
          fluxes%vprec(i,j)        = 0.0
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


  ! mask out land points and compute heat content of water fluxes
  ! assume liquid precip enters ocean at SST
  ! assume frozen precip enters ocean at 0degC
  ! assume liquid runoff enters ocean at SST
  ! assume solid runoff (calving) enters ocean at 0degC
  ! mass leaving ocean has heat_content determined in MOM_diabatic_driver.F90
  do j=js,je ; do i=is,ie
    fluxes%evap(i,j)    = fluxes%evap(i,j)    * G%mask2dT(i,j)
    fluxes%lprec(i,j)   = fluxes%lprec(i,j)   * G%mask2dT(i,j)
    fluxes%fprec(i,j)   = fluxes%fprec(i,j)   * G%mask2dT(i,j)
    fluxes%lrunoff(i,j) = fluxes%lrunoff(i,j) * G%mask2dT(i,j)
    fluxes%frunoff(i,j) = fluxes%frunoff(i,j) * G%mask2dT(i,j)
    fluxes%LW(i,j)      = fluxes%LW(i,j)      * G%mask2dT(i,j)
    fluxes%latent(i,j)  = fluxes%latent(i,j)  * G%mask2dT(i,j)
    fluxes%sens(i,j)    = fluxes%sens(i,j)    * G%mask2dT(i,j)
    fluxes%sw(i,j)      = fluxes%sw(i,j)      * G%mask2dT(i,j)

    fluxes%latent_evap_diag(i,j)     = fluxes%latent_evap_diag(i,j) * G%mask2dT(i,j)
    fluxes%latent_fprec_diag(i,j)    = -fluxes%fprec(i,j)*CS%latent_heat_fusion
    fluxes%latent_frunoff_diag(i,j)  = -fluxes%frunoff(i,j)*CS%latent_heat_fusion
  enddo ; enddo


!### if (associated(CS%ctrl_forcing_CSp)) then
!###   do j=js,je ; do i=is,ie
!###     SST_anom(i,j) = state%SST(i,j) - CS%T_Restore(i,j)
!###     SSS_anom(i,j) = state%SSS(i,j) - CS%S_Restore(i,j)
!###     SSS_mean(i,j) = 0.5*(state%SSS(i,j) + CS%S_Restore(i,j))
!###   enddo ; enddo
!###   call apply_ctrl_forcing(SST_anom, SSS_anom, SSS_mean, fluxes%heat_added, &
!###                           fluxes%vprec, day, dt, G, CS%ctrl_forcing_CSp)
!### endif

  call callTree_leave("buoyancy_forcing_from_data_override")
end subroutine buoyancy_forcing_from_data_override


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

  if (CS%use_temperature) then
    do j=js,je ; do i=is,ie
      fluxes%evap(i,j)                 = 0.0
      fluxes%lprec(i,j)                = 0.0
      fluxes%fprec(i,j)                = 0.0
      fluxes%vprec(i,j)                = 0.0
      fluxes%lrunoff(i,j)              = 0.0
      fluxes%frunoff(i,j)              = 0.0
      fluxes%lw(i,j)                   = 0.0
      fluxes%latent(i,j)               = 0.0
      fluxes%sens(i,j)                 = 0.0
      fluxes%sw(i,j)                   = 0.0
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


subroutine buoyancy_forcing_const(state, fluxes, day, dt, G, CS)
  type(surface),         intent(inout) :: state
  type(forcing),         intent(inout) :: fluxes
  type(time_type),       intent(in)    :: day
  real,                  intent(in)    :: dt
  type(ocean_grid_type), intent(in)    :: G
  type(surface_forcing_CS), pointer    :: CS

!  This subroutine specifies the current surface fluxes of buoyancy
!  temperature and fresh water.  It may also be modified to add
!  surface fluxes of user provided tracers.
!  We here define a constant surface heat flux.

! Arguments:
!  (inout)   state   = structure describing ocean surface state
!  (inout)   fluxes  = structure with pointers to forcing fields; unused have NULL ptrs
!  (in)      day     = time of the fluxes
!  (in)      dt      = amount of time over which the fluxes apply
!  (in)      G       = ocean grid structure
!  (in)      CS      = pointer to control struct returned by previous surface_forcing_init call

  integer :: i, j, is, ie, js, je
  call callTree_enter("buoyancy_forcing_const, MOM_surface_forcing.F90")
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec

  if (CS%use_temperature) then
    do j=js,je ; do i=is,ie
      fluxes%evap(i,j)                 = 0.0
      fluxes%lprec(i,j)                = 0.0
      fluxes%fprec(i,j)                = 0.0
      fluxes%vprec(i,j)                = 0.0
      fluxes%lrunoff(i,j)              = 0.0
      fluxes%frunoff(i,j)              = 0.0
      fluxes%lw(i,j)                   = 0.0
      fluxes%latent(i,j)               = 0.0
      fluxes%sens(i,j)                 = CS%constantHeatForcing * G%mask2dT(i,j)
      fluxes%sw(i,j)                   = 0.0
      fluxes%latent_evap_diag(i,j)     = 0.0
      fluxes%latent_fprec_diag(i,j)    = 0.0
      fluxes%latent_frunoff_diag(i,j)  = 0.0
    enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      fluxes%buoy(i,j) = 0.0
    enddo ; enddo
  endif

  call callTree_leave("buoyancy_forcing_const")
end subroutine buoyancy_forcing_const


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

  ! This case has no surface buoyancy forcing.
  if (CS%use_temperature) then
    do j=js,je ; do i=is,ie
      fluxes%evap(i,j)                 = 0.0
      fluxes%lprec(i,j)                = 0.0
      fluxes%fprec(i,j)                = 0.0
      fluxes%vprec(i,j)                = 0.0
      fluxes%lrunoff(i,j)              = 0.0
      fluxes%frunoff(i,j)              = 0.0
      fluxes%lw(i,j)                   = 0.0
      fluxes%latent(i,j)               = 0.0
      fluxes%sens(i,j)                 = 0.0
      fluxes%sw(i,j)                   = 0.0
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
          fluxes%heat_added(i,j) = G%mask2dT(i,j) * &
              ((T_Restore - state%SST(i,j)) * ((CS%Rho0 * fluxes%C_p) * CS%Flux_const))
          fluxes%vprec(i,j) = - (CS%Rho0*CS%Flux_const) * &
              (S_Restore - state%SSS(i,j)) / &
              (0.5*(state%SSS(i,j) + S_Restore))
        else
          fluxes%heat_added(i,j) = 0.0
          fluxes%vprec(i,j)        = 0.0
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

! Arguments:
!            CS              = pointer to control structure from previous surface_forcing_init call
!  (in)      G               = ocean grid structure
!  (in)      Time            = model time at this call; needed for mpp_write calls
!  (in, opt) directory       = optional directory into which to write these restart files
!  (in, opt) time_stamped    = if true, the restart file names include a unique time stamp
!                              default is false.
!  (in, opt) filename_suffix = optional suffix (e.g., a time-stamp) to append to the restart fname

  if (.not.associated(CS)) return
  if (.not.associated(CS%restart_CSp)) return

  call save_restart(directory, Time, G, CS%restart_CSp, time_stamped)

end subroutine forcing_save_restart


subroutine surface_forcing_init(Time, G, param_file, diag, CS, tracer_flow_CSp)
  type(time_type),              intent(in)    :: Time
  type(ocean_grid_type),        intent(in)    :: G
  type(param_file_type),        intent(in)    :: param_file
  type(diag_ctrl), target,      intent(inout) :: diag
  type(surface_forcing_CS),     pointer       :: CS
  type(tracer_flow_control_CS), pointer       :: tracer_flow_CSp

! Arguments:
!            Time            = current model time
!  (in)      G               = ocean grid structure
!  (in)      param_file      = structure indicating the open file to parse for model parameter values
!  (in)      diag            = structure used to regulate diagnostic output
!  (in/out)  CS              = pointer set to point to the control structure for this module
!  (in)      tracer_flow_CSp = pointer to the control structure of the tracer flow control module

  type(directories)  :: dirs
  logical            :: new_sim
  type(time_type)    :: Time_frc
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_surface_forcing" ! This module's name.
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
  call log_version(param_file, mod, version, '')
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
                 "(linear), (USER), (BFB) and (NONE).", fail_if_missing=.true.)
  if (trim(CS%buoy_config) == "file") then
    call get_param(param_file, mod, "ARCHAIC_OMIP_FORCING_FILE", CS%archaic_OMIP_file, &
                 "If true, use the forcing variable decomposition from \n"//&
                 "the old German OMIP prescription that predated CORE. If \n"//&
                 "false, use the variable groupings available from MOM \n"//&
                 "output diagnostics of forcing variables.", default=.true.)
    if (CS%archaic_OMIP_file) then
      call get_param(param_file, mod, "LONGWAVEDOWN_FILE", CS%longwave_file, &
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
      call get_param(param_file, mod, "SHORTWAVEDOWN_FILE", CS%shortwave_file, &
                 "The file with the downward shortwave heat flux.", &
                 fail_if_missing=.true.)
      call get_param(param_file, mod, "SNOW_FILE", CS%snow_file, &
                 "The file with the downward frozen precip flux, in \n"//&
                 "variable snow.", fail_if_missing=.true.)
      call get_param(param_file, mod, "PRECIP_FILE", CS%rain_file, &
                 "The file with the downward total precip flux, in \n"//&
                 "variable precip.", fail_if_missing=.true.)
      call get_param(param_file, mod, "FRESHDISCHARGE_FILE", CS%runoff_file, &
                 "The file with the fresh and frozen runoff/calving fluxes, \n"//&
                 "invariables disch_w and disch_s.", fail_if_missing=.true.)

      ! These variable names are hard-coded, per the archaic OMIP conventions.
      CS%latentheat_file = CS%evaporation_file ; CS%latent_var = "evap"
      CS%LW_var = "lwdn_sfc"; CS%SW_var = "swdn_sfc"; CS%sens_var = "shflx"
      CS%evap_var = "evap"; CS%rain_var = "precip"; CS%snow_var = "snow"
      CS%lrunoff_var = "disch_w"; CS%frunoff_var = "disch_s"

    else
      call get_param(param_file, mod, "LONGWAVE_FILE", CS%longwave_file, &
                 "The file with the longwave heat flux, in the variable \n"//&
                 "given by LONGWAVE_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mod, "LONGWAVE_FORCING_VAR", CS%LW_var, &
                 "The variable with the longwave forcing field.", default="LW")

      call get_param(param_file, mod, "SHORTWAVE_FILE", CS%shortwave_file, &
                 "The file with the shortwave heat flux, in the variable \n"//&
                 "given by SHORTWAVE_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mod, "SHORTWAVE_FORCING_VAR", CS%SW_var, &
                 "The variable with the shortwave forcing field.", default="SW")

      call get_param(param_file, mod, "EVAPORATION_FILE", CS%evaporation_file, &
                 "The file with the evaporative moisture flux, in the \n"//&
                 "variable given by EVAP_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mod, "EVAP_FORCING_VAR", CS%evap_var, &
                 "The variable with the evaporative moisture flux.", &
                 default="evap")

      call get_param(param_file, mod, "LATENTHEAT_FILE", CS%latentheat_file, &
                 "The file with the latent heat flux, in the variable \n"//&
                 "given by LATENT_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mod, "LATENT_FORCING_VAR", CS%latent_var, &
                 "The variable with the latent heat flux.", default="latent")

      call get_param(param_file, mod, "SENSIBLEHEAT_FILE", CS%sensibleheat_file, &
                 "The file with the sensible heat flux, in the variable \n"//&
                 "given by SENSIBLE_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mod, "SENSIBLE_FORCING_VAR", CS%sens_var, &
                 "The variable with the sensible heat flux.", default="sensible")

      call get_param(param_file, mod, "RAIN_FILE", CS%rain_file, &
                 "The file with the liquid precipitation flux, in the \n"//&
                 "variable given by RAIN_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mod, "RAIN_FORCING_VAR", CS%rain_var, &
                 "The variable with the liquid precipitation flux.", &
                 default="liq_precip")
      call get_param(param_file, mod, "SNOW_FILE", CS%snow_file, &
                 "The file with the frozen precipitation flux, in the \n"//&
                 "variable given by SNOW_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mod, "SNOW_FORCING_VAR", CS%snow_var, &
                 "The variable with the frozen precipitation flux.", &
                 default="froz_precip")

      call get_param(param_file, mod, "RUNOFF_FILE", CS%runoff_file, &
                 "The file with the fresh and frozen runoff/calving \n"//&
                 "fluxes, in variables given by LIQ_RUNOFF_FORCING_VAR \n"//&
                 "and FROZ_RUNOFF_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mod, "LIQ_RUNOFF_FORCING_VAR", CS%lrunoff_var, &
                 "The variable with the liquid runoff flux.", &
                 default="liq_runoff")
      call get_param(param_file, mod, "FROZ_RUNOFF_FORCING_VAR", CS%frunoff_var, &
                 "The variable with the frozen runoff flux.", &
                 default="froz_runoff")
    endif

    call get_param(param_file, mod, "SSTRESTORE_FILE", CS%SSTrestore_file, &
                 "The file with the SST toward which to restore in the \n"//&
                 "variable given by SST_RESTORE_VAR.", fail_if_missing=.true.)
    call get_param(param_file, mod, "SALINITYRESTORE_FILE", CS%salinityrestore_file, &
                 "The file with the surface salinity toward which to \n"//&
                 "restore in the variable given by SSS_RESTORE_VAR.", &
                 fail_if_missing=.true.)

    if (CS%archaic_OMIP_file) then
      CS%SST_restore_var = "TEMP" ; CS%SSS_restore_var = "SALT"
    else
      call get_param(param_file, mod, "SST_RESTORE_VAR", CS%SST_restore_var, &
                 "The variable with the SST toward which to restore.", &
                 default="SST")
      call get_param(param_file, mod, "SSS_RESTORE_VAR", CS%SSS_restore_var, &
                 "The variable with the SSS toward which to restore.", &
                 default="SSS")
    endif

    ! Add inputdir to the file names.
    CS%shortwave_file = trim(CS%inputdir)//trim(CS%shortwave_file)
    CS%longwave_file = trim(CS%inputdir)//trim(CS%longwave_file)
    CS%sensibleheat_file = trim(CS%inputdir)//trim(CS%sensibleheat_file)
    CS%latentheat_file = trim(CS%inputdir)//trim(CS%latentheat_file)
    CS%evaporation_file = trim(CS%inputdir)//trim(CS%evaporation_file)
    CS%snow_file = trim(CS%inputdir)//trim(CS%snow_file)
    CS%rain_file = trim(CS%inputdir)//trim(CS%rain_file)
    CS%runoff_file = trim(CS%inputdir)//trim(CS%runoff_file)

    CS%shortwaveup_file = trim(CS%inputdir)//trim(CS%shortwaveup_file)
    CS%longwaveup_file = trim(CS%inputdir)//trim(CS%longwaveup_file)

    CS%SSTrestore_file = trim(CS%inputdir)//trim(CS%SSTrestore_file)
    CS%salinityrestore_file = trim(CS%inputdir)//trim(CS%salinityrestore_file)
  elseif (trim(CS%buoy_config) == "const") then
    call get_param(param_file, mod, "SENSIBLE_HEAT_FLUX", CS%constantHeatForcing, &
                 "A constant heat forcing (positive into ocean) applied \n"//&
                 "through the sensible heat flux field. ", &
                 units='W/m2', fail_if_missing=.true.)
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
    call get_param(param_file, mod, "USTAR_FORCING_VAR", CS%ustar_var, &
                 "The name of the friction velocity variable in WIND_FILE \n"//&
                 "or blank to get ustar from the wind stresses plus the \n"//&
                 "gustiness.", default=" ", units="nondim")
    CS%wind_file = trim(CS%inputdir) // trim(CS%wind_file)
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
  if ((trim(CS%wind_config) == "2gyre") .or. &
      (trim(CS%wind_config) == "1gyre") .or. &
      (trim(CS%wind_config) == "gyres") .or. &
      (trim(CS%buoy_config) == "linear")) then
    CS%south_lat = G%south_lat
    CS%len_lat = G%len_lat
  endif
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
  call get_param(param_file, mod, "LATENT_HEAT_FUSION", CS%latent_heat_fusion, &
                 "The latent heat of fusion.", units="J/kg", default=hlf)
  call get_param(param_file, mod, "LATENT_HEAT_VAPORIZATION", CS%latent_heat_vapor, &
                 "The latent heat of fusion.", units="J/kg", default=hlv)
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
    call safe_alloc_ptr(CS%gust,G%isd,G%ied,G%jsd,G%jed)
    filename = trim(CS%inputdir) // trim(gust_file)
    call read_data(filename,'gustiness',CS%gust,domain=G%domain%mpp_domain, &
                   timelevel=1) ! units should be Pa
  endif

!  All parameter settings are now known.

  if (trim(CS%wind_config) == "USER" .or. trim(CS%buoy_config) == "USER" ) then
    call USER_surface_forcing_init(Time, G, param_file, diag, CS%user_forcing_CSp)
  elseif (trim(CS%buoy_config) == "BFB" ) then
    call BFB_surface_forcing_init(Time, G, param_file, diag, CS%BFB_forcing_CSp)
  elseif (trim(CS%wind_config) == "MESO" .or. trim(CS%buoy_config) == "MESO" ) then
    call MESO_surface_forcing_init(Time, G, param_file, diag, CS%MESO_forcing_CSp)
  elseif (trim(CS%wind_config) == "SCM_ideal_hurr") then
    call SCM_idealized_hurricane_wind_init(Time, G, param_file, CS%SCM_idealized_hurricane_CSp)
  elseif (trim(CS%wind_config) == "const") then
    call get_param(param_file, mod, "CONST_WIND_TAUX", CS%tau_x0, &
                 "With wind_config const, this is the constant zonal\n"//&
                 "wind-stress", units="Pa", fail_if_missing=.true.)
    call get_param(param_file, mod, "CONST_WIND_TAUY", CS%tau_y0, &
                 "With wind_config const, this is the constant zonal\n"//&
                 "wind-stress", units="Pa", fail_if_missing=.true.)
  elseif (trim(CS%wind_config) == "SCM_CVmix_tests" .or. &
          trim(CS%buoy_config) == "SCM_CVmix_tests") then
    call SCM_CVmix_tests_surface_forcing_init(Time, G, param_file, CS%SCM_CVmix_tests_CSp)
    CS%SCM_CVmix_tests_CSp%Rho0 = CS%Rho0 !copy reference density for pass
  endif

  call register_forcing_type_diags(Time, diag, CS%use_temperature, CS%handles)

  ! Set up any restart fields associated with the forcing.
  call restart_init(param_file, CS%restart_CSp, "MOM_forcing.res")
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

  ! Determine how many time levels are in each forcing variable.
  if (trim(CS%buoy_config) == "file") then
    CS%SW_nlev = num_timelevels(CS%shortwave_file, CS%SW_var, min_dims=3)
    CS%LW_nlev = num_timelevels(CS%longwave_file, CS%LW_var, min_dims=3)
    CS%latent_nlev = num_timelevels(CS%latentheat_file, CS%latent_var, 3)
    CS%sens_nlev = num_timelevels(CS%sensibleheat_file, CS%sens_var, min_dims=3)

    CS%evap_nlev = num_timelevels(CS%evaporation_file, CS%evap_var, min_dims=3)
    CS%precip_nlev = num_timelevels(CS%rain_file, CS%rain_var, min_dims=3)
    CS%runoff_nlev = num_timelevels(CS%runoff_file, CS%lrunoff_var, 3)

    CS%SST_nlev = num_timelevels(CS%SSTrestore_file, CS%SST_restore_var, 3)
    CS%SSS_nlev = num_timelevels(CS%salinityrestore_file, CS%SSS_restore_var, 3)
  endif

  if (trim(CS%wind_config) == "file") &
    CS%wind_nlev = num_timelevels(CS%wind_file, CS%stress_x_var, min_dims=3)

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
