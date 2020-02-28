module MOM_surface_forcing

! This file is part of MOM6. See LICENSE.md for the license.

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

use MOM_constants,           only : hlv, hlf
use MOM_cpu_clock,           only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,           only : CLOCK_MODULE
use MOM_diag_mediator,       only : post_data, query_averaging_enabled
use MOM_diag_mediator,       only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_domains,             only : pass_var, pass_vector, AGRID, To_South, To_West, To_All
use MOM_error_handler,       only : callTree_enter, callTree_leave
use MOM_error_handler,       only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,         only : get_param, log_version, param_file_type
use MOM_string_functions,    only : uppercase
use MOM_forcing_type,        only : forcing, mech_forcing
use MOM_forcing_type,        only : forcing_diags, mech_forcing_diags, register_forcing_type_diags
use MOM_forcing_type,        only : set_net_mass_forcing, copy_common_forcing_fields
use MOM_forcing_type,        only : set_derived_forcing_fields
use MOM_forcing_type,        only : allocate_forcing_type, deallocate_forcing_type
use MOM_forcing_type,        only : allocate_mech_forcing, deallocate_mech_forcing
use MOM_get_input,           only : Get_MOM_Input, directories
use MOM_grid,                only : ocean_grid_type
use MOM_io,                  only : file_exists, MOM_read_data, MOM_read_vector, slasher
use MOM_restart,             only : register_restart_field, restart_init, MOM_restart_CS
use MOM_restart,             only : restart_init_end, save_restart, restore_state
use MOM_time_manager,        only : time_type, operator(+), operator(/), get_time, set_time
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_unit_scaling,        only : unit_scale_type
use MOM_variables,           only : surface
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

  logical :: use_temperature    !< if true, temp & salinity used as state variables
  logical :: restorebuoy        !< if true, use restoring surface buoyancy forcing
  logical :: adiabatic          !< if true, no diapycnal mass fluxes or surface buoyancy forcing
  logical :: variable_winds     !< if true, wind stresses vary with time
  logical :: variable_buoyforce !< if true, buoyancy forcing varies with time.
  real    :: south_lat          !< southern latitude of the domain
  real    :: len_lat            !< domain length in latitude

  real :: Rho0                  !< Boussinesq reference density [R ~> kg m-3]
  real :: G_Earth               !< gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real :: Flux_const            !< piston velocity for surface restoring [Z T-1 ~> m s-1]
  real :: latent_heat_fusion    !< latent heat of fusion times [Q ~> J kg-1]
  real :: latent_heat_vapor     !< latent heat of vaporization [Q ~> J kg-1]

  real    :: gust_const                 !< constant unresolved background gustiness for ustar [R L Z T-1 ~> Pa]
  logical :: read_gust_2d               !< if true, use 2-dimensional gustiness supplied from a file
  real, pointer :: gust(:,:) => NULL()  !< spatially varying unresolved background gustiness [R L Z T-1 ~> Pa]
                                        !< gust is used when read_gust_2d is true.

  real, pointer :: T_Restore(:,:)    => NULL()  !< temperature to damp (restore) the SST to [degC]
  real, pointer :: S_Restore(:,:)    => NULL()  !< salinity to damp (restore) the SSS [ppt]
  real, pointer :: Dens_Restore(:,:) => NULL()  !< density to damp (restore) surface density [kg m-3]

  integer :: wind_last_lev_read = -1 !< The last time level read from the wind input files
  integer :: buoy_last_lev_read = -1 !< The last time level read from buoyancy input files

  ! if WIND_CONFIG=='gyres' then use the following as  = A, B, C and n respectively for
  ! taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L)
  real :: gyres_taux_const   !< A constant wind stress [Pa].
  real :: gyres_taux_sin_amp !< The amplitude of cosine wind stress gyres [Pa], if WIND_CONFIG=='gyres'.
  real :: gyres_taux_cos_amp !< The amplitude of cosine wind stress gyres [Pa], if WIND_CONFIG=='gyres'.
  real :: gyres_taux_n_pis   !< The number of sine lobes in the basin if  if WIND_CONFIG=='gyres'

  real :: T_north   !< target temperatures at north used in buoyancy_forcing_linear
  real :: T_south   !< target temperatures at south used in buoyancy_forcing_linear
  real :: S_north   !< target salinity at north used in buoyancy_forcing_linear
  real :: S_south   !< target salinity at south used in buoyancy_forcing_linear

  logical :: first_call_set_forcing = .true. !< True until after the first call to set_forcing

  real :: wind_scale          !< value by which wind-stresses are scaled, ND.
  character(len=8)   :: wind_stagger !< A character indicating how the wind stress components
                              !! are staggered in WIND_FILE.  Valid values are A or C for now.

  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp => NULL() !< A pointer to the structure
                              !! that is used to orchestrate the calling of tracer packages
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(diag_ctrl), pointer :: diag !< structure used to regulate timing of diagnostic output

  character(len=200) :: inputdir    !< directory where NetCDF input files are.
  character(len=200) :: wind_config !< indicator for wind forcing type (2gyre, USER, FILE..)
  character(len=200) :: wind_file   !< if wind_config is "file", file to use
  character(len=200) :: buoy_config !< indicator for buoyancy forcing type

  character(len=200) :: longwavedown_file  = '' !< The file from which the downward longwave heat flux is read
  character(len=200) :: shortwavedown_file = '' !< The file from which the downward shortwave heat flux is read
  character(len=200) :: evaporation_file  = '' !< The file from which the evaporation is read
  character(len=200) :: sensibleheat_file = '' !< The file from which the sensible heat flux is read
  character(len=200) :: latentheat_file   = '' !< The file from which the latent heat flux is read

  character(len=200) :: precip_file   = '' !< The file from which the rainfall is read
  character(len=200) :: snow_file   = '' !< The file from which the snowfall is read
  character(len=200) :: freshdischarge_file = '' !< The file from which the runoff and calving are read

  character(len=200) :: longwaveup_file  = '' !< The file from which the upward longwave heat flux is read
  character(len=200) :: shortwaveup_file = '' !< The file from which the upward shorwave heat flux is read

  character(len=200) :: SSTrestore_file      = '' !< The file from which to read the sea surface
                                                  !! temperature to restore toward
  character(len=200) :: salinityrestore_file = '' !< The file from which to read the sea surface
                                                  !! salinity to restore toward

  character(len=80)  :: stress_x_var  = '' !< X-windstress variable name in the input file
  character(len=80)  :: stress_y_var  = '' !< Y-windstress variable name in the input file

  type(forcing_diags), public :: handles !< A structure with diagnostics handles

  !>@{ Control structures for named forcing packages
  type(user_revise_forcing_CS),  pointer :: urf_CS => NULL()
  type(user_surface_forcing_CS), pointer :: user_forcing_CSp => NULL()
  ! type(MESO_surface_forcing_CS), pointer :: MESO_forcing_CSp => NULL()
  !!@}
end type surface_forcing_CS

integer :: id_clock_forcing

contains

!> This subroutine calls other subroutines in this file to get surface forcing fields.
!! It also allocates and initializes the fields in the flux type.
subroutine set_forcing(sfc_state, forcing, fluxes, day_start, day_interval, G, US, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(mech_forcing),    intent(inout) :: forces !< A structure with the driving mechanical forces
  type(forcing),         intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),       intent(in)    :: day_start !< The start time of the fluxes
  type(time_type),       intent(in)    :: day_interval !< Length of time over which these fluxes applied
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer    :: CS   !< A pointer to the control structure returned by a
                                               !! previous surface_forcing_init call

  ! Local variables
  real :: dt                     ! length of time over which fluxes applied [s]
  type(time_type) :: day_center  ! central time of the fluxes.
  integer :: intdt
  integer :: isd, ied, jsd, jed
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call cpu_clock_begin(id_clock_forcing)

  day_center = day_start + day_interval/2
  call get_time(day_interval, intdt)
  dt = real(intdt)

  if (CS%first_call_set_forcing) then
    ! Allocate memory for the mechanical and thermodyanmic forcing fields.
    call allocate_mech_forcing(G, forces, stress=.true., ustar=.true., press=.true.)

    call allocate_forcing_type(G, fluxes, ustar=.true.)
    if (trim(CS%buoy_config) /= "NONE") then
      if ( CS%use_temperature ) then
        call allocate_forcing_type(G, fluxes, water=.true., heat=.true., press=.true.)
        if (CS%restorebuoy) then
          call safe_alloc_ptr(CS%T_Restore,isd, ied, jsd, jed)
          call safe_alloc_ptr(fluxes%heat_added, isd, ied, jsd, jed)
          call safe_alloc_ptr(CS%S_Restore, isd, ied, jsd, jed)
        endif
      else ! CS%use_temperature false.
        call safe_alloc_ptr(fluxes%buoy, isd, ied, jsd, jed)

        if (CS%restorebuoy) call safe_alloc_ptr(CS%Dens_Restore, isd, ied, jsd, jed)
      endif  ! endif for  CS%use_temperature
    endif
  endif

  ! calls to various wind options
  if (CS%variable_winds .or. CS%first_call_set_forcing) then
    if (trim(CS%wind_config) == "file") then
      call wind_forcing_from_file(sfc_state, forces, day_center, G, CS)
    elseif (trim(CS%wind_config) == "2gyre") then
      call wind_forcing_2gyre(sfc_state, forces, day_center, G, CS)
    elseif (trim(CS%wind_config) == "1gyre") then
      call wind_forcing_1gyre(sfc_state, forces, day_center, G, CS)
    elseif (trim(CS%wind_config) == "gyres") then
      call wind_forcing_gyres(sfc_state, forces, day_center, G, CS)
    elseif (trim(CS%wind_config) == "zero") then
      call wind_forcing_zero(sfc_state, forces, day_center, G, CS)
    elseif (trim(CS%wind_config) == "MESO") then
      call MOM_error(FATAL, "MESO forcing is not available with the ice-shelf"//&
               "version of MOM_surface_forcing.")
!      call MESO_wind_forcing(sfc_state, forces, day_center, G, CS%MESO_forcing_CSp)
    elseif (trim(CS%wind_config) == "USER") then
      call USER_wind_forcing(sfc_state, forces, day_center, G, CS%user_forcing_CSp)
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
      call buoyancy_forcing_from_files(sfc_state, fluxes, day_center, dt, G, US, CS)
    elseif (trim(CS%buoy_config) == "zero") then
      call buoyancy_forcing_zero(sfc_state, fluxes, day_center, dt, G, CS)
    elseif (trim(CS%buoy_config) == "linear") then
      call buoyancy_forcing_linear(sfc_state, fluxes, day_center, dt, G, US, CS)
    elseif (trim(CS%buoy_config) == "MESO") then
      call MOM_error(FATAL, "MESO forcing is not available with the ice-shelf"//&
               "version of MOM_surface_forcing.")
!      call MESO_buoyancy_forcing(sfc_state, fluxes, day_center, dt, G, CS%MESO_forcing_CSp)
    elseif (trim(CS%buoy_config) == "USER") then
      call USER_buoyancy_forcing(sfc_state, fluxes, day_center, dt, G, CS%user_forcing_CSp)
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
    call call_tracer_set_forcing(sfc_state, fluxes, day_start, day_interval, G, CS%tracer_flow_CSp)
  endif

  ! Allow for user-written code to alter the fluxes after all the above
  call user_alter_forcing(sfc_state, fluxes, day_center, G, CS%urf_CS)

  ! Fields that exist in both the forcing and mech_forcing types must be copied.
  if (CS%variable_winds .or. CS%first_call_set_forcing) then
    call copy_common_forcing_fields(forces, fluxes, G)
    call set_derived_forcing_fields(forces, fluxes, G, US, CS%Rho0)
  endif

  if ((CS%variable_buoyforce .or. CS%first_call_set_forcing) .and. &
      (.not.CS%adiabatic)) then
    call set_net_mass_forcing(fluxes, forces, G, US)
  endif

  CS%first_call_set_forcing = .false.

  call cpu_clock_end(id_clock_forcing)
end subroutine set_forcing

!> This subroutine allocates arrays for buoyancy forcing.
subroutine buoyancy_forcing_allocate(fluxes, G, CS)
  type(forcing),         intent(inout) :: fluxes !< A structure with pointers to thermodynamic
                                               !! forcing fields that will be allocated here
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(surface_forcing_CS), pointer    :: CS   !< A pointer to the control structure returned by a
                                               !! previous surface_forcing_init call

  integer :: isd, ied, jsd, jed
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if ( CS%use_temperature ) then
    call allocate_forcing_type(G, fluxes, water=.true., heat=.true., &
                               ustar=.true., press=.true.)

    ! surface restoring fields
    if (CS%restorebuoy) then
      call safe_alloc_ptr(CS%T_Restore,isd,ied,jsd,jed)
      call safe_alloc_ptr(fluxes%heat_added,isd,ied,jsd,jed)
      call safe_alloc_ptr(CS%S_Restore,isd,ied,jsd,jed)
    endif

  else ! CS%use_temperature false.
    call safe_alloc_ptr(fluxes%buoy,isd,ied,jsd,jed)

    call allocate_forcing_type(G, fluxes, water=.true., heat=.true., &
                               ustar=.true., press=.true.)

    if (CS%restorebuoy) call safe_alloc_ptr(CS%Dens_Restore,isd,ied,jsd,jed)

  endif  ! endif for  CS%use_temperature

end subroutine buoyancy_forcing_allocate


! This subroutine sets the surface wind stresses to zero
subroutine wind_forcing_zero(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day    !< Time used for determining the fluxes.
  type(ocean_grid_type),    intent(in)    :: G      !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS     !< A pointer to the control structure returned by a
                                                    !! previous surface_forcing_init call

  real :: PI
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  call callTree_enter("wind_forcing_zero, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  !set steady surface wind stresses, in units of Pa.
  PI = 4.0*atan(1.0)

  do j=js,je ; do I=Isq,Ieq
    forces%taux(I,j) = 0.0
  enddo ; enddo

  do J=Jsq,Jeq ; do i=is,ie
    forces%tauy(i,J) = 0.0
  enddo ; enddo

  if (CS%read_gust_2d) then
    if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
      forces%ustar(i,j) = sqrt(US%L_to_Z*CS%gust(i,j)/CS%Rho0)
    enddo ; enddo ; endif
  else
    if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
      forces%ustar(i,j) = sqrt(US%L_to_Z*CS%gust_const/CS%Rho0)
    enddo ; enddo ; endif
  endif

  call callTree_leave("wind_forcing_zero")
end subroutine wind_forcing_zero


!> This subroutine sets the surface wind stresses according to double gyre.
subroutine wind_forcing_2gyre(sfc_state, forces, day, G, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day    !< Time used for determining the fluxes.
  type(ocean_grid_type),    intent(in)    :: G      !< The ocean's grid structure
  type(surface_forcing_CS), pointer       :: CS     !< A pointer to the control structure returned by a
                                                    !! previous surface_forcing_init call

  ! Local variables
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
    forces%taux(I,j) = 0.1*US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z * &
                       (1.0 - cos(2.0*PI*(G%geoLatCu(I,j)-CS%South_lat) / CS%len_lat))
  enddo ; enddo

  do J=Jsq,Jeq ; do i=is,ie
    forces%tauy(i,J) = 0.0
  enddo ; enddo

  call callTree_leave("wind_forcing_2gyre")
end subroutine wind_forcing_2gyre


!> This subroutine sets the surface wind stresses according to single gyre.
subroutine wind_forcing_1gyre(sfc_state, forces, day, G, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day    !< Time used for determining the fluxes.
  type(ocean_grid_type),    intent(in)    :: G      !< The ocean's grid structure
  type(surface_forcing_CS), pointer       :: CS     !< A pointer to the control structure returned by a
                                                    !! previous surface_forcing_init call

  ! Local variables
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
    forces%taux(I,j) = -0.2*US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z * &
                       cos(PI*(G%geoLatCu(I,j)-CS%South_lat)/CS%len_lat)
  enddo ; enddo

  do J=Jsq,Jeq ; do i=is,ie
    forces%tauy(i,J) = 0.0
  enddo ; enddo

  call callTree_leave("wind_forcing_1gyre")
end subroutine wind_forcing_1gyre


!> This subroutine sets the surface wind stresses according to gyres.
subroutine wind_forcing_gyres(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day    !< Time used for determining the fluxes.
  type(ocean_grid_type),    intent(in)    :: G      !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS     !< A pointer to the control structure returned by a
                                                    !! previous surface_forcing_init call

  ! Local variables
  real :: PI, y
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  call callTree_enter("wind_forcing_gyres, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! steady surface wind stresses [Pa]
  PI = 4.0*atan(1.0)

  do j=jsd,jed ; do I=IsdB,IedB
    y = (G%geoLatCu(I,j)-CS%South_lat)/CS%len_lat
    forces%taux(I,j) = US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z * (CS%gyres_taux_const + &
             (   CS%gyres_taux_sin_amp*sin(CS%gyres_taux_n_pis*PI*y)    &
               + CS%gyres_taux_cos_amp*cos(CS%gyres_taux_n_pis*PI*y) ))
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    forces%tauy(i,J) = 0.0
  enddo ; enddo

  ! set the friction velocity
  do j=js,je ; do i=is,ie
    forces%ustar(i,j) = sqrt(US%L_to_S * (CS%gust_const/CS%Rho0 + &
            sqrt(0.5*(forces%tauy(i,j-1)*forces%tauy(i,j-1) + forces%tauy(i,j)*forces%tauy(i,j) + &
                      forces%taux(i-1,j)*forces%taux(i-1,j) + forces%taux(i,j)*forces%taux(i,j)))/CS%Rho0) )
  enddo ; enddo

  call callTree_leave("wind_forcing_gyres")
end subroutine wind_forcing_gyres

!> This subroutine sets the surface wind stresses by reading a file.
subroutine wind_forcing_from_file(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day    !< Time used for determining the fluxes.
  type(ocean_grid_type),    intent(inout) :: G      !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS     !< A pointer to the control structure returned by a
                                                    !! previous surface_forcing_init call

  ! Local variables
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: time_lev             ! With fields from a file, this must
                                  ! be reset, depending on the time.
  character(len=200) :: filename  ! The name of the input file.
  real :: temp_x(SZI_(G),SZJ_(G)) ! Pseudo-zonal and psuedo-meridional
  real :: temp_y(SZI_(G),SZJ_(G)) ! wind stresses at h-points [Pa].
  real    :: Pa_conversion           ! A unit conversion factor from Pa to the internal wind stress
                                     ! units [R Z L T-2 Pa-1 ~> 1]
  integer :: days, seconds

  call callTree_enter("wind_forcing_from_file, MOM_surface_forcing.F90")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  Pa_conversion = US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z

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
      call MOM_read_vector(filename, CS%stress_x_var, CS%stress_y_var, &
                           temp_x(:,:), temp_y(:,:), G%Domain, stagger=AGRID, &
                           timelevel=time_lev, scale=Pa_conversion)

      call pass_vector(temp_x, temp_y, G%Domain, To_All, AGRID)
      do j=js,je ; do I=Isq,Ieq
        forces%taux(I,j) = 0.5 * CS%wind_scale * (temp_x(i,j) + temp_x(i+1,j))
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        forces%tauy(i,J) = 0.5 * CS%wind_scale * (temp_y(i,j) + temp_y(i,j+1))
      enddo ; enddo

      if (CS%read_gust_2d) then
        do j=js,je ; do i=is,ie
          forces%ustar(i,j) = sqrt(US%L_to_Z * (CS%gust(i,j) + &
                  sqrt(temp_x(i,j)*temp_x(i,j) + temp_y(i,j)*temp_y(i,j)) ) / CS%Rho0)
        enddo ; enddo
      else
        do j=js,je ; do i=is,ie
          forces%ustar(i,j) = sqrt(US%L_to_Z * (CS%gust_const/CS%Rho0 + &
                  sqrt(temp_x(i,j)*temp_x(i,j) + temp_y(i,j)*temp_y(i,j)) / CS%Rho0) )
        enddo ; enddo
      endif
    case ("C")
      call MOM_read_vector(filename,CS%stress_x_var, CS%stress_y_var, &
                     forces%taux(:,:), forces%tauy(:,:), &
                     G%Domain, timelevel=time_lev, &
                     scale=Pa_conversion)
      if (CS%wind_scale /= 1.0) then
        do j=js,je ; do I=Isq,Ieq
          forces%taux(I,j) = CS%wind_scale * forces%taux(I,j)
        enddo ; enddo
        do J=Jsq,Jeq ; do i=is,ie
          forces%tauy(i,J) = CS%wind_scale * forces%tauy(i,J)
        enddo ; enddo
      endif

      call pass_vector(forces%taux, forces%tauy, G%Domain, To_All)
      if (CS%read_gust_2d) then
        do j=js, je ; do i=is, ie
          forces%ustar(i,j) = sqrt( (CS%gust(i,j) + &
                  sqrt(0.5*((forces%tauy(i,j-1)**2 + forces%tauy(i,j)**2) + &
                            (forces%taux(i-1,j)**2 +  forces%taux(i,j)**2)))) * US%L_to_Z / CS%Rho0 )
        enddo ; enddo
      else
        do j=js, je ; do i=is, ie
          forces%ustar(i,j) = sqrt(US%L_to_Z * (CS%gust_const/CS%Rho0 + &
                 sqrt(0.5*((forces%tauy(i,j-1)**2 + forces%tauy(i,j)**2) + &
                           (forces%taux(i-1,j)**2 + forces%taux(i,j)**2))) / CS%Rho0) )
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


!> This subroutine specifies the current surface fluxes of buoyancy, temperature and fresh water
!! by reading a file. It may also be modified to add surface fluxes of user provided tracers.
subroutine buoyancy_forcing_from_files(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),         intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),       intent(in)    :: day    !< Time used for determining the fluxes.
  real,                  intent(in)    :: dt   !< The amount of time over which
                                               !! the fluxes apply [s]
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer    :: CS   !< A pointer to the control structure returned by a
                                               !! previous surface_forcing_init call

  real :: rhoXcp ! mean density times the heat capacity [Q R degC-1 ~> J m-3 degC-1].
  real :: Irho0  ! inverse Boussinesq reference density [m3 kg-1].
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed

  integer :: time_lev           ! With fields from a file, this must
                                ! be reset, depending on the time.
  integer :: time_lev_monthly   ! With fields from a file, this must
                                ! be reset, depending on the time.
  integer :: days, seconds
  real, dimension(SZI_(G),SZJ_(G)) :: &
    temp, &       ! A 2-d temporary work array with various units.
    SST_anom, &   ! Instantaneous sea surface temperature anomalies from a
                  ! target (observed) value [degC].
    SSS_anom, &   ! Instantaneous sea surface salinity anomalies from a target
                  ! (observed) value [ppt].
    SSS_mean      ! A (mean?) salinity about which to normalize local salinity
                  ! anomalies when calculating restorative precipitation
                  ! anomalies [ppt].

  call callTree_enter("buoyancy_forcing_from_files, MOM_surface_forcing.F90")

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ! allocate and initialize arrays
  call buoyancy_forcing_allocate(fluxes, G, CS)

  if (CS%use_temperature) rhoXcp = CS%Rho0 * fluxes%C_p
  Irho0 = 1.0/(US%R_to_kg_m3*CS%Rho0)

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


    call MOM_read_data(trim(CS%inputdir)//trim(CS%longwavedown_file), "lwdn_sfc", &
             fluxes%LW(:,:), G%Domain, timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
    call MOM_read_data(trim(CS%inputdir)//trim(CS%longwaveup_file), "lwup_sfc", &
             temp(:,:), G%Domain, timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
    do j=js,je ; do i=is,ie ; fluxes%LW(i,j) = fluxes%LW(i,j) - temp(i,j) ; enddo ; enddo

    call MOM_read_data(trim(CS%inputdir)//trim(CS%evaporation_file), "evap", &
             fluxes%evap(:,:), G%Domain, timelevel=time_lev, scale=-US%kg_m2s_to_RZ_T)
    do j=js,je ; do i=is,ie
      fluxes%latent(i,j)           = CS%latent_heat_vapor*fluxes%evap(i,j)
      fluxes%latent_evap_diag(i,j) = fluxes%latent(i,j)
    enddo ; enddo

    call MOM_read_data(trim(CS%inputdir)//trim(CS%sensibleheat_file), "shflx", &
             fluxes%sens(:,:), G%Domain, timelevel=time_lev, scale=-US%W_m2_to_QRZ_T)

    call MOM_read_data(trim(CS%inputdir)//trim(CS%shortwavedown_file), "swdn_sfc", &
             fluxes%sw(:,:), G%Domain, timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
    call MOM_read_data(trim(CS%inputdir)//trim(CS%shortwaveup_file), "swup_sfc", &
             temp(:,:), G%Domain, timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
    do j=js,je ; do i=is,ie
      fluxes%sw(i,j) = fluxes%sw(i,j) - temp(i,j)
    enddo ; enddo

    call MOM_read_data(trim(CS%inputdir)//trim(CS%snow_file), "snow", &
             fluxes%fprec(:,:), G%Domain, timelevel=time_lev, scale=US%kg_m2s_to_RZ_T)
    call MOM_read_data(trim(CS%inputdir)//trim(CS%precip_file), "precip", &
             fluxes%lprec(:,:), G%Domain, timelevel=time_lev, scale=US%kg_m2s_to_RZ_T)
    do j=js,je ; do i=is,ie
      fluxes%lprec(i,j) = fluxes%lprec(i,j) - fluxes%fprec(i,j)
    enddo ; enddo

    call MOM_read_data(trim(CS%inputdir)//trim(CS%freshdischarge_file), "disch_w", &
             temp(:,:), G%Domain, timelevel=time_lev_monthly, scale=US%kg_m2s_to_RZ_T)
    do j=js,je ; do i=is,ie
      fluxes%lrunoff(i,j) = temp(i,j)*US%m_to_L**2*G%IareaT(i,j)
    enddo ; enddo
    call MOM_read_data(trim(CS%inputdir)//trim(CS%freshdischarge_file), "disch_s", &
              temp(:,:), G%Domain, timelevel=time_lev_monthly, scale=US%kg_m2s_to_RZ_T)
    do j=js,je ; do i=is,ie
      fluxes%frunoff(i,j) = temp(i,j)*US%m_to_L**2*G%IareaT(i,j)
    enddo ; enddo

!     Read the SST and SSS fields for damping.
    if (CS%restorebuoy) then
      call MOM_read_data(trim(CS%inputdir)//trim(CS%SSTrestore_file), "TEMP", &
               CS%T_Restore(:,:), G%Domain, timelevel=time_lev_monthly)
      call MOM_read_data(trim(CS%inputdir)//trim(CS%salinityrestore_file), "SALT", &
               CS%S_Restore(:,:), G%Domain, timelevel=time_lev_monthly)
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

      fluxes%heat_content_lrunoff(i,j) = fluxes%C_p * fluxes%lrunoff(i,j)*sfc_state%SST(i,j)
      fluxes%latent_evap_diag(i,j)     = fluxes%latent_evap_diag(i,j) * G%mask2dT(i,j)
      fluxes%latent_fprec_diag(i,j)    = -fluxes%fprec(i,j)*CS%latent_heat_fusion
      fluxes%latent_frunoff_diag(i,j)  = -fluxes%frunoff(i,j)*CS%latent_heat_fusion
    enddo ; enddo

  endif ! time_lev /= CS%buoy_last_lev_read

  if (CS%restorebuoy) then
    if (CS%use_temperature) then
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0) then
          fluxes%heat_added(i,j) = G%mask2dT(i,j) * &
              ((CS%T_Restore(i,j) - sfc_state%SST(i,j)) * rhoXcp * CS%Flux_const)
          fluxes%vprec(i,j) = - (CS%Rho0*CS%Flux_const) * &
              (CS%S_Restore(i,j) - sfc_state%SSS(i,j)) / &
              (0.5*(sfc_state%SSS(i,j) + CS%S_Restore(i,j)))
        else
          fluxes%heat_added(i,j) = 0.0
          fluxes%vprec(i,j) = 0.0
        endif
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0) then
          fluxes%buoy(i,j) = (CS%Dens_Restore(i,j) - sfc_state%sfc_density(i,j)) * &
                             (CS%G_Earth * CS%Flux_const/(US%R_to_kg_m3*CS%Rho0))
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

  call callTree_leave("buoyancy_forcing_from_files")
end subroutine buoyancy_forcing_from_files


!> This subroutine specifies the current surface fluxes of buoyancy, temperature and fresh water.
!! It may also be modified to add surface fluxes of user provided tracers.
!! This case has zero surface buoyancy forcing.
subroutine buoyancy_forcing_zero(sfc_state, fluxes, day, dt, G, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),         intent(inout) :: fluxes !< A structure with pointers to thermodynamic forcing fields
  type(time_type),       intent(in)    :: day    !< Time used for determining the fluxes.
  real,                  intent(in)    :: dt   !< The amount of time over which
                                               !! the fluxes apply [s]
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(surface_forcing_CS), pointer    :: CS   !< A pointer to the control structure returned by a
                                               !! previous surface_forcing_init call

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

!> This subroutine specifies the current surface fluxes of buoyancy, temperature and fresh water.
!! It may also be modified to add surface fluxes of user provided tracers.
subroutine buoyancy_forcing_linear(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),         intent(inout) :: fluxes !< A structure with pointers to thermodynamic forcing fields
  type(time_type),       intent(in)    :: day    !< Time used for determining the fluxes.
  real,                  intent(in)    :: dt   !< The amount of time over which
                                               !! the fluxes apply, in s
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer    :: CS   !< A pointer to the control structure returned by a
                                               !! previous surface_forcing_init call

  ! Local variables
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
          fluxes%heat_added(i,j) = G%mask2dT(i,j) * &
              ((T_Restore - sfc_state%SST(i,j)) * ((CS%Rho0 * fluxes%C_p) * CS%Flux_const))
          fluxes%vprec(i,j) = - (CS%Rho0*CS%Flux_const) * &
              (S_Restore - sfc_state%SSS(i,j)) / &
              (0.5*(sfc_state%SSS(i,j) + S_Restore))
        else
          fluxes%heat_added(i,j) = 0.0
          fluxes%vprec(i,j) = 0.0
        endif
      enddo ; enddo
    else
      call MOM_error(FATAL, "buoyancy_forcing_linear in MOM_surface_forcing: "// &
                     "RESTOREBUOY to linear not written yet.")
     !do j=js,je ; do i=is,ie
     !  if (G%mask2dT(i,j) > 0) then
     !   fluxes%buoy(i,j) = US%kg_m3_to_R*(CS%Dens_Restore(i,j) - sfc_state%sfc_density(i,j)) * &
     !                       (CS%G_Earth * CS%Flux_const / CS%Rho0)
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

!> Save any restart files associated with the surface forcing.
subroutine forcing_save_restart(CS, G, Time, directory, time_stamped, &
                                filename_suffix)
  type(surface_forcing_CS),   pointer       :: CS   !< A pointer to the control structure returned
                                                    !! by a previous call to surface_forcing_init
  type(ocean_grid_type),      intent(inout) :: G    !< The ocean's grid structure
  type(time_type),            intent(in)    :: Time !< The current model time
  character(len=*),           intent(in)    :: directory !< The directory into which to write the
                                                    !! restart files
  logical,          optional, intent(in)    :: time_stamped !< If true, the restart file names include
                                                    !! a unique time stamp.  The default is false.
  character(len=*), optional, intent(in)    :: filename_suffix !< An optional suffix (e.g., a time-
                                                    !! stamp) to append to the restart file names.

  if (.not.associated(CS)) return
  if (.not.associated(CS%restart_CSp)) return

  call save_restart(directory, Time, 1, G, CS%restart_CSp, time_stamped)

end subroutine forcing_save_restart

!> Initialize the surface forcing, including setting parameters and allocating permanent memory.
subroutine surface_forcing_init(Time, G, US, param_file, diag, CS, tracer_flow_CSp)
  type(time_type),           intent(in) :: Time !< The current model time
  type(ocean_grid_type),     intent(in) :: G    !< The ocean's grid structure
  type(unit_scale_type),     intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),     intent(in) :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,   intent(in) :: diag !< A structure that is used to regulate diagnostic output.
  type(surface_forcing_CS),  pointer    :: CS   !< A pointer that is set to point to the control structure
                                                !! for this module
  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp !< A pointer to the control structure of
                                                !! the tracer flow control module.

  ! Local variables
  type(directories)  :: dirs
  logical            :: new_sim
  type(time_type)    :: Time_frc
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_surface_forcing" ! This module's name.
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
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state "//&
                 "variables.", default=.true.)
  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, &
                 "The directory in which all input files are found.", &
                 default=".")
  CS%inputdir = slasher(CS%inputdir)

  call get_param(param_file, mdl, "ADIABATIC", CS%adiabatic, &
                 "There are no diapycnal mass fluxes if ADIABATIC is "//&
                 "true. This assumes that KD = KDML = 0.0 and that "//&
                 "there is no buoyancy forcing, but makes the model "//&
                 "faster by eliminating subroutine calls.", default=.false.)
  call get_param(param_file, mdl, "VARIABLE_WINDS", CS%variable_winds, &
                 "If true, the winds vary in time after the initialization.", &
                 default=.true.)
  call get_param(param_file, mdl, "VARIABLE_BUOYFORCE", CS%variable_buoyforce, &
                 "If true, the buoyancy forcing varies in time after the "//&
                 "initialization of the model.", default=.true.)

  call get_param(param_file, mdl, "BUOY_CONFIG", CS%buoy_config, &
                 "The character string that indicates how buoyancy forcing "//&
                 "is specified. Valid options include (file), (zero), "//&
                 "(linear), (USER), and (NONE).", fail_if_missing=.true.)
  if (trim(CS%buoy_config) == "file") then
    call get_param(param_file, mdl, "LONGWAVEDOWN_FILE", CS%longwavedown_file, &
                 "The file with the downward longwave heat flux, in "//&
                 "variable lwdn_sfc.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "LONGWAVEUP_FILE", CS%longwaveup_file, &
                 "The file with the upward longwave heat flux, in "//&
                 "variable lwup_sfc.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "EVAPORATION_FILE", CS%evaporation_file, &
                 "The file with the evaporative moisture flux, in "//&
                 "variable evap.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "SENSIBLEHEAT_FILE", CS%sensibleheat_file, &
                 "The file with the sensible heat flux, in "//&
                 "variable shflx.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "SHORTWAVEUP_FILE", CS%shortwaveup_file, &
                 "The file with the upward shortwave heat flux.", &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, "SHORTWAVEDOWN_FILE", CS%shortwavedown_file, &
                 "The file with the downward shortwave heat flux.", &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, "SNOW_FILE", CS%snow_file, &
                 "The file with the downward frozen precip flux, in "//&
                 "variable snow.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "PRECIP_FILE", CS%precip_file, &
                 "The file with the downward total precip flux, in "//&
                 "variable precip.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "FRESHDISCHARGE_FILE", CS%freshdischarge_file, &
                 "The file with the fresh and frozen runoff/calving fluxes, "//&
                 "invariables disch_w and disch_s.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "SSTRESTORE_FILE", CS%SSTrestore_file, &
                 "The file with the SST toward which to restore in "//&
                 "variable TEMP.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "SALINITYRESTORE_FILE", CS%salinityrestore_file, &
                 "The file with the surface salinity toward which to "//&
                 "restore in variable SALT.", fail_if_missing=.true.)
  endif
  call get_param(param_file, mdl, "WIND_CONFIG", CS%wind_config, &
                 "The character string that indicates how wind forcing "//&
                 "is specified. Valid options include (file), (2gyre), "//&
                 "(1gyre), (gyres), (zero), and (USER).", fail_if_missing=.true.)
  if (trim(CS%wind_config) == "file") then
    call get_param(param_file, mdl, "WIND_FILE", CS%wind_file, &
                 "The file in which the wind stresses are found in "//&
                 "variables STRESS_X and STRESS_Y.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "WINDSTRESS_X_VAR",CS%stress_x_var, &
                 "The name of the x-wind stress variable in WIND_FILE.", &
                 default="STRESS_X")
    call get_param(param_file, mdl, "WINDSTRESS_Y_VAR", CS%stress_y_var, &
                 "The name of the y-wind stress variable in WIND_FILE.", &
                 default="STRESS_Y")
    call get_param(param_file, mdl, "WINDSTRESS_STAGGER",CS%wind_stagger, &
                 "A character indicating how the wind stress components "//&
                 "are staggered in WIND_FILE.  This may be A or C for now.", &
                 default="A")
    call get_param(param_file, mdl, "WINDSTRESS_SCALE", CS%wind_scale, &
                 "A value by which the wind stresses in WIND_FILE are rescaled.", &
                 default=1.0, units="nondim")
  endif
  if (trim(CS%wind_config) == "gyres") then
    call get_param(param_file, mdl, "TAUX_CONST", CS%gyres_taux_const, &
                 "With the gyres wind_config, the constant offset in the "//&
                 "zonal wind stress profile: "//&
                 "  A in taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L).", &
                 units="Pa", default=0.0)
    call get_param(param_file, mdl, "TAUX_SIN_AMP",CS%gyres_taux_sin_amp, &
                 "With the gyres wind_config, the sine amplitude in the "//&
                 "zonal wind stress profile: "//&
                 "  B in taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L).", &
                 units="Pa", default=0.0)
    call get_param(param_file, mdl, "TAUX_COS_AMP",CS%gyres_taux_cos_amp, &
                 "With the gyres wind_config, the cosine amplitude in "//&
                 "the zonal wind stress profile: "//&
                 "  C in taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L).", &
                 units="Pa", default=0.0)
    call get_param(param_file, mdl, "TAUX_N_PIS",CS%gyres_taux_n_pis, &
                 "With the gyres wind_config, the number of gyres in "//&
                 "the zonal wind stress profile: "//&
                 "  n in taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L).", &
                 units="nondim", default=0.0)
  endif
  call get_param(param_file, mdl, "SOUTHLAT", CS%south_lat, &
                 "The southern latitude of the domain or the equivalent "//&
                 "starting value for the y-axis.", units=axis_units, default=0.)
  call get_param(param_file, mdl, "LENLAT", CS%len_lat, &
                 "The latitudinal or y-direction length of the domain.", &
                 units=axis_units, fail_if_missing=.true.)
   call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINSEQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "RESTOREBUOY", CS%restorebuoy, &
                 "If true, the buoyancy fluxes drive the model back "//&
                 "toward some specified surface state with a rate "//&
                 "given by FLUXCONST.", default= .false.)
  call get_param(param_file, mdl, "LATENT_HEAT_FUSION", CS%latent_heat_fusion, &
                 "The latent heat of fusion.", default=hlf, &
                 units="J/kg", scale=US%J_kg_to_Q)
  call get_param(param_file, mdl, "LATENT_HEAT_VAPORIZATION", CS%latent_heat_vapor, &
                 "The latent heat of fusion.", default=hlv, units="J/kg", scale=US%J_kg_to_Q)
  if (CS%restorebuoy) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes "//&
                 "to the relative surface anomalies (akin to a piston "//&
                 "velocity).  Note the non-MKS units.", &
                 units="m day-1", scale=US%m_to_Z*US%T_to_s/86400.0, fail_if_missing=.true.)
    if (trim(CS%buoy_config) == "linear") then
      call get_param(param_file, mdl, "SST_NORTH", CS%T_north, &
                 "With buoy_config linear, the sea surface temperature "//&
                 "at the northern end of the domain toward which to "//&
                 "to restore.", units="deg C", default=0.0)
      call get_param(param_file, mdl, "SST_SOUTH", CS%T_south, &
                 "With buoy_config linear, the sea surface temperature "//&
                 "at the southern end of the domain toward which to "//&
                 "to restore.", units="deg C", default=0.0)
      call get_param(param_file, mdl, "SSS_NORTH", CS%S_north, &
                 "With buoy_config linear, the sea surface salinity "//&
                 "at the northern end of the domain toward which to "//&
                 "to restore.", units="PSU", default=35.0)
      call get_param(param_file, mdl, "SSS_SOUTH", CS%S_south, &
                 "With buoy_config linear, the sea surface salinity "//&
                 "at the southern end of the domain toward which to "//&
                 "to restore.", units="PSU", default=35.0)
    endif
  endif
  call get_param(param_file, mdl, "G_EARTH", CS%G_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80, scale=US%m_to_L**2*US%Z_to_m*US%T_to_s**2)

  call get_param(param_file, mdl, "GUST_CONST", CS%gust_const, &
                 "The background gustiness in the winds.", &
                 units="Pa", default=0.02, scale=US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z)
  call get_param(param_file, mdl, "READ_GUST_2D", CS%read_gust_2d, &
                 "If true, use a 2-dimensional gustiness supplied from "//&
                 "an input file", default=.false.)
  if (CS%read_gust_2d) then
    call get_param(param_file, mdl, "GUST_2D_FILE", gust_file, &
                 "The file in which the wind gustiness is found in "//&
                 "variable gustiness.", fail_if_missing=.true.)
    call safe_alloc_ptr(CS%gust,G%isd,G%ied,G%jsd,G%jed) ; CS%gust(:,:) = 0.0
    filename = trim(CS%inputdir) // trim(gust_file)
    call MOM_read_data(filename,'gustiness',CS%gust,G%domain, timelevel=1, &
                 scale=US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z) ! units in file should be Pa
  endif
  call get_param(param_file, mdl, "AXIS_UNITS", axis_units, default="degrees")

!  All parameter settings are now known.

  if (trim(CS%wind_config) == "USER" .or. trim(CS%buoy_config) == "USER" ) then
    call USER_surface_forcing_init(Time, G, param_file, diag, CS%user_forcing_CSp)
  elseif (trim(CS%wind_config) == "MESO" .or. trim(CS%buoy_config) == "MESO" ) then
    call MOM_error(FATAL, "MESO forcing is not available with the ice-shelf"//&
               "version of MOM_surface_forcing.")
  endif

  call register_forcing_type_diags(Time, diag, US, CS%use_temperature, CS%handles)

  ! Set up any restart fields associated with the forcing.
  call restart_init(G, param_file, CS%restart_CSp, "MOM_forcing.res")
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

  call user_revise_forcing_init(param_file, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)
end subroutine surface_forcing_init

!> Clean up and deallocate any memory associated with this module and its children.
subroutine surface_forcing_end(CS, fluxes)
  type(surface_forcing_CS), pointer       :: CS     !< A pointer to the control structure returned
                                                    !! by a previous surface_forcing_init call
                                                    !! that will be deallocated here.
  type(forcing), optional,  intent(inout) :: fluxes !< A structure containing pointers to any possible
                                                    !! forcing fields that will be deallocated here.

  if (present(fluxes)) call deallocate_forcing_type(fluxes)

  if (associated(CS)) deallocate(CS)
  CS => NULL()

end subroutine surface_forcing_end

end module MOM_surface_forcing
