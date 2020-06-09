!> Functions that calculate the surface wind stresses and fluxes of buoyancy
!! or temperature/salinity andfresh water, in ocean-only (solo) mode.
!!
!! These functions are called every time step, even if the wind stresses
!! or buoyancy fluxes are constant in time - in that case these routines
!! return quickly without doing anything.  In addition, any I/O of forcing
!! fields is controlled by surface_forcing_init, located in this file.
module MOM_surface_forcing

! This file is part of MOM6. See LICENSE.md for the license.

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
use MOM_forcing_type,        only : forcing, mech_forcing
use MOM_forcing_type,        only : set_net_mass_forcing, copy_common_forcing_fields
use MOM_forcing_type,        only : set_derived_forcing_fields
use MOM_forcing_type,        only : forcing_diags, mech_forcing_diags, register_forcing_type_diags
use MOM_forcing_type,        only : allocate_forcing_type, deallocate_forcing_type
use MOM_forcing_type,        only : allocate_mech_forcing, deallocate_mech_forcing
use MOM_grid,                only : ocean_grid_type
use MOM_get_input,           only : Get_MOM_Input, directories
use MOM_io,                  only : file_exists, MOM_read_data, MOM_read_vector, slasher
use MOM_io,                  only : EAST_FACE, NORTH_FACE, num_timelevels
use MOM_restart,             only : register_restart_field, restart_init, MOM_restart_CS
use MOM_restart,             only : restart_init_end, save_restart, restore_state
use MOM_time_manager,        only : time_type, operator(+), operator(/), get_time, time_type_to_real
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_unit_scaling,        only : unit_scale_type
use MOM_variables,           only : surface
use MESO_surface_forcing,    only : MESO_buoyancy_forcing
use MESO_surface_forcing,    only : MESO_surface_forcing_init, MESO_surface_forcing_CS
use user_surface_forcing,    only : USER_wind_forcing, USER_buoyancy_forcing
use user_surface_forcing,    only : USER_surface_forcing_init, user_surface_forcing_CS
use user_revise_forcing,     only : user_alter_forcing, user_revise_forcing_init
use user_revise_forcing,     only : user_revise_forcing_CS
use idealized_hurricane, only : idealized_hurricane_wind_init
use idealized_hurricane, only : idealized_hurricane_wind_forcing, SCM_idealized_hurricane_wind_forcing
use idealized_hurricane, only : idealized_hurricane_CS
use SCM_CVmix_tests,         only : SCM_CVmix_tests_surface_forcing_init
use SCM_CVmix_tests,         only : SCM_CVmix_tests_wind_forcing
use SCM_CVmix_tests,         only : SCM_CVmix_tests_buoyancy_forcing
use SCM_CVmix_tests,         only : SCM_CVmix_tests_CS
use BFB_surface_forcing,    only : BFB_buoyancy_forcing
use BFB_surface_forcing,    only : BFB_surface_forcing_init, BFB_surface_forcing_CS
use dumbbell_surface_forcing,    only : dumbbell_surface_forcing_init, dumbbell_surface_forcing_CS
use dumbbell_surface_forcing, only    : dumbbell_buoyancy_forcing
use data_override_mod, only : data_override_init, data_override

implicit none ; private

#include <MOM_memory.h>

public set_forcing
public surface_forcing_init
public forcing_save_restart

!> Structure containing pointers to the forcing fields that may be used to drive MOM.
!!  All fluxes are positive into the ocean.
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
  real :: Flux_const_T          !< piston velocity for surface temperature restoring [m s-1]
  real :: Flux_const_S          !< piston velocity for surface salinity restoring [Z T-1 ~> m s-1]
  real :: latent_heat_fusion    !< latent heat of fusion times [Q ~> J kg-1]
  real :: latent_heat_vapor     !< latent heat of vaporization [Q ~> J kg-1]
  real :: tau_x0                !< Constant zonal wind stress used in the WIND_CONFIG="const" forcing
  real :: tau_y0                !< Constant meridional wind stress used in the WIND_CONFIG="const" forcing

  real    :: gust_const                 !< constant unresolved background gustiness for ustar [R L Z T-1 ~> Pa]
  logical :: read_gust_2d               !< if true, use 2-dimensional gustiness supplied from a file
  real, pointer :: gust(:,:) => NULL()  !< spatially varying unresolved background gustiness [R L Z T-1 ~> Pa]
                                        !! gust is used when read_gust_2d is true.

  real, pointer :: T_Restore(:,:)    => NULL()  !< temperature to damp (restore) the SST to [degC]
  real, pointer :: S_Restore(:,:)    => NULL()  !< salinity to damp (restore) the SSS [ppt]
  real, pointer :: Dens_Restore(:,:) => NULL()  !< density to damp (restore) surface density [R ~> kg m-3]

  integer :: buoy_last_lev_read = -1 !< The last time level read from buoyancy input files

  ! if WIND_CONFIG=='gyres' then use the following as  = A, B, C and n respectively for
  ! taux = A + B*sin(n*pi*y/L) + C*cos(n*pi*y/L)
  real :: gyres_taux_const   !< A constant wind stress [Pa].
  real :: gyres_taux_sin_amp !< The amplitude of cosine wind stress gyres [Pa], if WIND_CONFIG=='gyres'.
  real :: gyres_taux_cos_amp !< The amplitude of cosine wind stress gyres [Pa], if WIND_CONFIG=='gyres'.
  real :: gyres_taux_n_pis   !< The number of sine lobes in the basin if  if WIND_CONFIG=='gyres'
  logical :: answers_2018    !< If true, use the order of arithmetic and expressions that recover
                             !! the answers from the end of 2018.  Otherwise, use a form of the gyre
                             !! wind stresses that are rotationally invariant and more likely to be
                             !! the same between compilers.
  logical :: fix_ustar_gustless_bug         !< If true correct a bug in the time-averaging of the
                                            !! gustless wind friction velocity.
  ! if WIND_CONFIG=='scurves' then use the following to define a piecwise scurve profile
  real :: scurves_ydata(20) = 90. !< Latitudes of scurve nodes [degreesN]
  real :: scurves_taux(20) = 0.   !< Zonal wind stress values at scurve nodes [Pa]

  real :: T_north   !< target temperatures at north used in buoyancy_forcing_linear
  real :: T_south   !< target temperatures at south used in buoyancy_forcing_linear
  real :: S_north   !< target salinity at north used in buoyancy_forcing_linear
  real :: S_south   !< target salinity at south used in buoyancy_forcing_linear

  logical :: first_call_set_forcing = .true. !< True until after the first call to set_forcing
  logical :: archaic_OMIP_file = .true. !< If true use the variable names and data fields from
                                        !! a very old version of the OMIP forcing
  logical :: dataOverrideIsInitialized = .false. !< If true, data override has been initialized

  real :: wind_scale          !< value by which wind-stresses are scaled, ND.
  real :: constantHeatForcing !< value used for sensible heat flux when buoy_config="const" [Q R Z T-1 ~> W m-2]

  character(len=8)   :: wind_stagger !< A character indicating how the wind stress components
                              !! are staggered in WIND_FILE.  Valid values are A or C for now.
  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp => NULL() !< A pointer to the structure
                              !! that is used to orchestrate the calling of tracer packages
!#CTRL#  type(ctrl_forcing_CS), pointer :: ctrl_forcing_CSp => NULL()
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(diag_ctrl), pointer :: diag !< structure used to regulate timing of diagnostic output

  character(len=200) :: inputdir    !< directory where NetCDF input files are.
  character(len=200) :: wind_config !< indicator for wind forcing type (2gyre, USER, FILE..)
  character(len=200) :: wind_file   !< if wind_config is "file", file to use
  character(len=200) :: buoy_config !< indicator for buoyancy forcing type

  character(len=200) :: longwave_file     = '' !< The file from which the longwave heat flux is read
  character(len=200) :: shortwave_file    = '' !< The file from which the shortwave heat flux is read
  character(len=200) :: evaporation_file  = '' !< The file from which the evaporation is read
  character(len=200) :: sensibleheat_file = '' !< The file from which the sensible heat flux is read
  character(len=200) :: latentheat_file   = '' !< The file from which the latent heat flux is read

  character(len=200) :: rain_file   = '' !< The file from which the rainfall is read
  character(len=200) :: snow_file   = '' !< The file from which the snowfall is read
  character(len=200) :: runoff_file = '' !< The file from which the runoff is read

  character(len=200) :: longwaveup_file  = '' !< The file from which the upward longwave heat flux is read
  character(len=200) :: shortwaveup_file = '' !< The file from which the upward shorwave heat flux is read

  character(len=200) :: SSTrestore_file      = '' !< The file from which to read the sea surface
                                                  !! temperature to restore toward
  character(len=200) :: salinityrestore_file = '' !< The file from which to read the sea surface
                                                  !! salinity to restore toward

  character(len=80)  :: stress_x_var  = '' !< X-windstress variable name in the input file
  character(len=80)  :: stress_y_var  = '' !< Y-windstress variable name in the input file
  character(len=80)  :: ustar_var     = '' !< ustar variable name in the input file
  character(len=80)  :: LW_var        = '' !< lonngwave heat flux variable name in the input file
  character(len=80)  :: SW_var        = '' !< shortwave heat flux variable name in the input file
  character(len=80)  :: latent_var    = '' !< latent heat flux variable name in the input file
  character(len=80)  :: sens_var      = '' !< sensible heat flux variable name in the input file
  character(len=80)  :: evap_var      = '' !< evaporation variable name in the input file
  character(len=80)  :: rain_var      = '' !< rainfall variable name in the input file
  character(len=80)  :: snow_var      = '' !< snowfall variable name in the input file
  character(len=80)  :: lrunoff_var   = '' !< liquid runoff variable name in the input file
  character(len=80)  :: frunoff_var   = '' !< frozen runoff variable name in the input file
  character(len=80)  :: SST_restore_var = '' !< target sea surface temeperature variable name in the input file
  character(len=80)  :: SSS_restore_var = '' !< target sea surface salinity variable name in the input file

  ! These variables give the number of time levels in the various forcing files.
  integer :: wind_nlev = -1   !< The number of time levels in the file of wind stress
  integer :: SW_nlev   = -1   !< The number of time levels in the file of shortwave heat flux
  integer :: LW_nlev = -1     !< The number of time levels in the file of longwave heat flux
  integer :: latent_nlev = -1 !< The number of time levels in the file of latent heat flux
  integer :: sens_nlev = -1   !< The number of time levels in the file of sensible heat flux
  integer :: evap_nlev = -1   !< The number of time levels in the file of evaporation
  integer :: precip_nlev = -1 !< The number of time levels in the file of precipitation
  integer :: runoff_nlev = -1 !< The number of time levels in the file of runoff
  integer :: SST_nlev  = -1   !< The number of time levels in the file of target SST
  integer :: SSS_nlev = -1    !< The number of time levels in the file of target SSS

  ! These variables give the last time level read for the various forcing files.
  integer :: wind_last_lev = -1   !< The last time level read of wind stress
  integer :: SW_last_lev   = -1   !< The last time level read of shortwave heat flux
  integer :: LW_last_lev = -1     !< The last time level read of longwave heat flux
  integer :: latent_last_lev = -1 !< The last time level read of latent heat flux
  integer :: sens_last_lev = -1   !< The last time level read of sensible heat flux
  integer :: evap_last_lev = -1   !< The last time level read of evaporation
  integer :: precip_last_lev = -1 !< The last time level read of precipitation
  integer :: runoff_last_lev = -1 !< The last time level read of runoff
  integer :: SST_last_lev  = -1   !< The last time level read of target SST
  integer :: SSS_last_lev = -1    !< The last time level read of target SSS

  type(forcing_diags), public :: handles !< A structure with diagnostics handles

  !>@{ Control structures for named forcing packages
  type(user_revise_forcing_CS),  pointer :: urf_CS => NULL()
  type(user_surface_forcing_CS), pointer :: user_forcing_CSp => NULL()
  type(BFB_surface_forcing_CS), pointer :: BFB_forcing_CSp => NULL()
  type(dumbbell_surface_forcing_CS), pointer :: dumbbell_forcing_CSp => NULL()
  type(MESO_surface_forcing_CS), pointer :: MESO_forcing_CSp => NULL()
  type(idealized_hurricane_CS), pointer :: idealized_hurricane_CSp => NULL()
  type(SCM_CVmix_tests_CS),      pointer :: SCM_CVmix_tests_CSp => NULL()
  !>@}

end type surface_forcing_CS

integer :: id_clock_forcing !< A CPU time clock

contains

!> Calls subroutines in this file to get surface forcing fields.
!!
!! It also allocates and initializes the fields in the forcing and mech_forcing types
!! the first time it is called.
subroutine set_forcing(sfc_state, forces, fluxes, day_start, day_interval, G, US, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(mech_forcing),    intent(inout) :: forces !< A structure with the driving mechanical forces
  type(forcing),         intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),       intent(in)    :: day_start !< The start time of the fluxes
  type(time_type),       intent(in)    :: day_interval !< Length of time over which these fluxes applied
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer    :: CS   !< pointer to control struct returned by
                                               !! a previous surface_forcing_init call
  ! Local variables
  real :: dt                     ! length of time over which fluxes applied [s]
  type(time_type) :: day_center  ! central time of the fluxes.
  integer :: isd, ied, jsd, jed
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call cpu_clock_begin(id_clock_forcing)
  call callTree_enter("set_forcing, MOM_surface_forcing.F90")

  day_center = day_start + day_interval/2
  dt = time_type_to_real(day_interval)

  if (CS%first_call_set_forcing) then
    ! Allocate memory for the mechanical and thermodyanmic forcing fields.
    call allocate_mech_forcing(G, forces, stress=.true., ustar=.true., press=.true.)

    call allocate_forcing_type(G, fluxes, ustar=.true., fix_accum_bug=CS%fix_ustar_gustless_bug)
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
      call wind_forcing_from_file(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "data_override") then
      call wind_forcing_by_data_override(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "2gyre") then
      call wind_forcing_2gyre(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "1gyre") then
      call wind_forcing_1gyre(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "gyres") then
      call wind_forcing_gyres(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "zero") then
      call wind_forcing_const(sfc_state, forces, 0., 0., day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "const") then
      call wind_forcing_const(sfc_state, forces, CS%tau_x0, CS%tau_y0, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "Neverworld" .or. trim(CS%wind_config) == "Neverland") then
      call Neverworld_wind_forcing(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "scurves") then
      call scurve_wind_forcing(sfc_state, forces, day_center, G, US, CS)
    elseif (trim(CS%wind_config) == "ideal_hurr") then
      call idealized_hurricane_wind_forcing(sfc_state, forces, day_center, G, US, CS%idealized_hurricane_CSp)
    elseif (trim(CS%wind_config) == "SCM_ideal_hurr") then
      call SCM_idealized_hurricane_wind_forcing(sfc_state, forces, day_center, G, US, CS%idealized_hurricane_CSp)
    elseif (trim(CS%wind_config) == "SCM_CVmix_tests") then
      call SCM_CVmix_tests_wind_forcing(sfc_state, forces, day_center, G, US, CS%SCM_CVmix_tests_CSp)
    elseif (trim(CS%wind_config) == "USER") then
      call USER_wind_forcing(sfc_state, forces, day_center, G, US, CS%user_forcing_CSp)
    elseif (CS%variable_winds .and. .not.CS%first_call_set_forcing) then
      call MOM_error(FATAL, &
       "MOM_surface_forcing: Variable winds defined with no wind config")
    else
       call MOM_error(FATAL, &
       "MOM_surface_forcing:Unrecognized wind config "//trim(CS%wind_config))
    endif
  endif

  ! calls to various buoyancy forcing options
  if ((CS%variable_buoyforce .or. CS%first_call_set_forcing) .and. &
      (.not.CS%adiabatic)) then
    if (trim(CS%buoy_config) == "file") then
      call buoyancy_forcing_from_files(sfc_state, fluxes, day_center, dt, G, US, CS)
    elseif (trim(CS%buoy_config) == "data_override") then
      call buoyancy_forcing_from_data_override(sfc_state, fluxes, day_center, dt, G, US, CS)
    elseif (trim(CS%buoy_config) == "zero") then
      call buoyancy_forcing_zero(sfc_state, fluxes, day_center, dt, G, CS)
    elseif (trim(CS%buoy_config) == "const") then
      call buoyancy_forcing_const(sfc_state, fluxes, day_center, dt, G, US, CS)
    elseif (trim(CS%buoy_config) == "linear") then
      call buoyancy_forcing_linear(sfc_state, fluxes, day_center, dt, G, US, CS)
    elseif (trim(CS%buoy_config) == "MESO") then
      call MESO_buoyancy_forcing(sfc_state, fluxes, day_center, dt, G, US, CS%MESO_forcing_CSp)
    elseif (trim(CS%buoy_config) == "SCM_CVmix_tests") then
      call SCM_CVmix_tests_buoyancy_forcing(sfc_state, fluxes, day_center, G, US, CS%SCM_CVmix_tests_CSp)
    elseif (trim(CS%buoy_config) == "USER") then
      call USER_buoyancy_forcing(sfc_state, fluxes, day_center, dt, G, US, CS%user_forcing_CSp)
    elseif (trim(CS%buoy_config) == "BFB") then
      call BFB_buoyancy_forcing(sfc_state, fluxes, day_center, dt, G, US, CS%BFB_forcing_CSp)
    elseif (trim(CS%buoy_config) == "dumbbell") then
      call dumbbell_buoyancy_forcing(sfc_state, fluxes, day_center, dt, G, US, CS%dumbbell_forcing_CSp)
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
  call callTree_leave("set_forcing")

end subroutine set_forcing

!> Sets the surface wind stresses to constant values
subroutine wind_forcing_const(sfc_state, forces, tau_x0, tau_y0, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  real,                     intent(in)    :: tau_x0 !< The zonal wind stress [Pa]
  real,                     intent(in)    :: tau_y0 !< The meridional wind stress [Pa]
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control struct returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables
  real :: Pa_conversion ! A unit conversion factor from Pa to the internal units [R Z L T-2 Pa-1 ~> 1]
  real :: mag_tau
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq

  call callTree_enter("wind_forcing_const, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  Pa_conversion = US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z

  !set steady surface wind stresses, in units of Pa.
  mag_tau = Pa_conversion * sqrt( tau_x0**2 + tau_y0**2)

  do j=js,je ; do I=is-1,Ieq
    forces%taux(I,j) = tau_x0 * Pa_conversion
  enddo ; enddo

  do J=js-1,Jeq ; do i=is,ie
    forces%tauy(i,J) = tau_y0 * Pa_conversion
  enddo ; enddo

  if (CS%read_gust_2d) then
    if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
      forces%ustar(i,j) = sqrt( US%L_to_Z * ( mag_tau + CS%gust(i,j) ) / CS%Rho0 )
    enddo ; enddo ; endif
  else
    if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
      forces%ustar(i,j) = sqrt( US%L_to_Z * ( mag_tau + CS%gust_const ) / CS%Rho0 )
    enddo ; enddo ; endif
  endif

  call callTree_leave("wind_forcing_const")
end subroutine wind_forcing_const


!> Sets the surface wind stresses to set up two idealized gyres.
subroutine wind_forcing_2gyre(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control struct returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables
  real :: PI
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq

  call callTree_enter("wind_forcing_2gyre, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB

  !set the steady surface wind stresses, in units of Pa.
  PI = 4.0*atan(1.0)

  do j=js,je ; do I=is-1,Ieq
    forces%taux(I,j) = 0.1*US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z * &
                      (1.0 - cos(2.0*PI*(G%geoLatCu(I,j)-CS%South_lat) / CS%len_lat))
  enddo ; enddo

  do J=js-1,Jeq ; do i=is,ie
    forces%tauy(i,J) = 0.0
  enddo ; enddo

  call callTree_leave("wind_forcing_2gyre")
end subroutine wind_forcing_2gyre


!> Sets the surface wind stresses to set up a single idealized gyre.
subroutine wind_forcing_1gyre(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control struct returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables
  real :: PI
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq

  call callTree_enter("wind_forcing_1gyre, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB

  ! set the steady surface wind stresses, in units of Pa.
  PI = 4.0*atan(1.0)

  do j=js,je ; do I=is-1,Ieq
    forces%taux(I,j) = -0.2*US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z * &
                       cos(PI*(G%geoLatCu(I,j)-CS%South_lat)/CS%len_lat)
  enddo ; enddo

  do J=js-1,Jeq ; do i=is,ie
    forces%tauy(i,J) = 0.0
  enddo ; enddo

  call callTree_leave("wind_forcing_1gyre")
end subroutine wind_forcing_1gyre

!> Sets the surface wind stresses to set up idealized gyres.
subroutine wind_forcing_gyres(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control struct returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables
  real :: PI, y, I_rho
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq

  call callTree_enter("wind_forcing_gyres, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB

  ! steady surface wind stresses [Pa]
  PI = 4.0*atan(1.0)

  do j=js-1,je+1 ; do I=is-1,Ieq
    y = (G%geoLatCu(I,j)-CS%South_lat) / CS%len_lat
    forces%taux(I,j) = US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z * &
            (CS%gyres_taux_const +                            &
             (   CS%gyres_taux_sin_amp*sin(CS%gyres_taux_n_pis*PI*y)    &
               + CS%gyres_taux_cos_amp*cos(CS%gyres_taux_n_pis*PI*y) ))
  enddo ; enddo

  do J=js-1,Jeq ; do i=is-1,ie+1
    forces%tauy(i,J) = 0.0
  enddo ; enddo

  ! set the friction velocity
  if (CS%answers_2018) then
    do j=js,je ; do i=is,ie
      forces%ustar(i,j) = sqrt(US%L_to_Z * ((CS%gust_const/CS%Rho0) + &
              sqrt(0.5*(forces%tauy(i,j-1)*forces%tauy(i,j-1) + forces%tauy(i,j)*forces%tauy(i,j) + &
                        forces%taux(i-1,j)*forces%taux(i-1,j) + forces%taux(i,j)*forces%taux(i,j)))/CS%Rho0) )
    enddo ; enddo
  else
    I_rho = US%L_to_Z / CS%Rho0
    do j=js,je ; do i=is,ie
      forces%ustar(i,j) = sqrt( (CS%gust_const + &
            sqrt(0.5*((forces%tauy(i,J-1)**2 + forces%tauy(i,J)**2) + &
                      (forces%taux(I-1,j)**2 + forces%taux(I,j)**2))) ) * I_rho )
    enddo ; enddo
  endif

  call callTree_leave("wind_forcing_gyres")
end subroutine wind_forcing_gyres

!> Sets the surface wind stresses, forces%taux and forces%tauy for the
!! Neverworld forcing configuration.
subroutine Neverworld_wind_forcing(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day    !< Time used for determining the fluxes.
  type(ocean_grid_type),    intent(inout) :: G      !< Grid structure.
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS     !< pointer to control struct returned by
                                                    !! a previous surface_forcing_init call
  ! Local variables
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real :: PI, I_rho, y
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
    y = (G%geoLatT(i,j)-G%south_lat)/G%len_lat

    if (y <= 0.29) then
      forces%taux(I,j) = forces%taux(I,j) + tau_max * ( (1/0.29)*y - ( 1/(2*PI) )*sin( (2*PI*y) / 0.29 ) )
    endif
    if ((y > 0.29) .and. (y <= (0.8-off))) then
      forces%taux(I,j) = forces%taux(I,j) + tau_max *(0.35+0.65*cos(PI*(y-0.29)/(0.51-off))  )
    endif
    if ((y > (0.8-off)) .and. (y <= (1-off))) then
      forces%taux(I,j) = forces%taux(I,j) + tau_max *( 1.5*( (y-1+off) - (0.1/PI)*sin(10.0*PI*(y-0.8+off)) ) )
    endif
    forces%taux(I,j) = G%mask2dCu(I,j) * forces%taux(I,j)
  enddo ; enddo

  do J=js-1,Jeq ; do i=is,ie
    forces%tauy(i,J) = G%mask2dCv(i,J) * 0.0
  enddo ; enddo

  ! Set the surface friction velocity, in units of m s-1.  ustar is always positive.
  if (associated(forces%ustar)) then
    I_rho = US%L_to_Z / CS%Rho0
    do j=js,je ; do i=is,ie
      forces%ustar(i,j) = sqrt( (CS%gust_const + &
            sqrt(0.5*((forces%tauy(i,J-1)**2 + forces%tauy(i,J)**2) + &
                      (forces%taux(I-1,j)**2 + forces%taux(I,j)**2))) ) * I_rho )
    enddo ; enddo
  endif

end subroutine Neverworld_wind_forcing

!> Sets the zonal wind stresses to a piecewise series of s-curves.
subroutine scurve_wind_forcing(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day    !< Time used for determining the fluxes.
  type(ocean_grid_type),    intent(inout) :: G      !< Grid structure.
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS     !< pointer to control struct returned by
                                                    !! a previous surface_forcing_init call
  ! Local variables
  integer :: i, j, kseg
  real :: lon, lat, I_rho, y, L
! real :: ydata(7) = (/ -70., -45., -15., 0., 15., 45., 70. /)
! real :: taudt(7) = (/ 0., 0.2, -0.1, -0.02, -0.1, 0.1, 0. /)

  ! Allocate the forcing arrays, if necessary.
  call allocate_mech_forcing(G, forces, stress=.true.)

  kseg = 1
  do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
    lon = G%geoLonCu(I,j)
    lat = G%geoLatCu(I,j)

    ! Find segment k s.t. ydata(k)<= lat < ydata(k+1)
    do while (lat>=CS%scurves_ydata(kseg+1) .and. kseg<6)
      kseg = kseg+1
    enddo
    do while (lat<CS%scurves_ydata(kseg) .and. kseg>1)
      kseg = kseg-1
    enddo

    y = lat - CS%scurves_ydata(kseg)
    L = CS%scurves_ydata(kseg+1) - CS%scurves_ydata(kseg)
    forces%taux(I,j) = CS%scurves_taux(kseg) +  &
                       ( CS%scurves_taux(kseg+1) - CS%scurves_taux(kseg) ) * scurve(y, L)
    forces%taux(I,j) = G%mask2dCu(I,j) * forces%taux(I,j)
  enddo ; enddo

  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
    forces%tauy(i,J) = G%mask2dCv(i,J) * 0.0
  enddo ; enddo

  ! Set the surface friction velocity, in units of m s-1.  ustar is always positive.
  if (associated(forces%ustar)) then
    I_rho = US%L_to_Z / CS%Rho0
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      forces%ustar(i,j) = sqrt( (CS%gust_const + &
            sqrt(0.5*((forces%tauy(i,J-1)**2 + forces%tauy(i,J)**2) + &
                      (forces%taux(I-1,j)**2 + forces%taux(I,j)**2))) ) * I_rho )
    enddo ; enddo
  endif

end subroutine scurve_wind_forcing

!> Returns the value of a cosine-bell function evaluated at x/L
real function scurve(x,L)
  real , intent(in) :: x       !< non-dimensional position
  real , intent(in) :: L       !< non-dimensional width
  real :: s

  s = x/L
  scurve = (3. - 2.*s) * (s*s)
end function scurve

! Sets the surface wind stresses from input files.
subroutine wind_forcing_from_file(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),    intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control struct returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables
  character(len=200) :: filename  ! The name of the input file.
  real    :: temp_x(SZI_(G),SZJ_(G)) ! Pseudo-zonal and psuedo-meridional
  real    :: temp_y(SZI_(G),SZJ_(G)) ! wind stresses at h-points [R L Z T-1 ~> Pa].
  real    :: Pa_conversion           ! A unit conversion factor from Pa to the internal wind stress
                                     ! units [R Z L T-2 Pa-1 ~> 1]
  integer :: time_lev_daily          ! The time levels to read for fields with
  integer :: time_lev_monthly        ! daily and montly cycles.
  integer :: time_lev                ! The time level that is used for a field.
  integer :: days, seconds
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  logical :: read_Ustar

  call callTree_enter("wind_forcing_from_file, MOM_surface_forcing.F90")
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  Pa_conversion = US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z

  call get_time(day, seconds, days)
  time_lev_daily = days - 365*floor(real(days) / 365.0)

  if (time_lev_daily < 31) then ; time_lev_monthly = 0
  elseif (time_lev_daily < 59)  then ; time_lev_monthly = 1
  elseif (time_lev_daily < 90)  then ; time_lev_monthly = 2
  elseif (time_lev_daily < 120) then ; time_lev_monthly = 3
  elseif (time_lev_daily < 151) then ; time_lev_monthly = 4
  elseif (time_lev_daily < 181) then ; time_lev_monthly = 5
  elseif (time_lev_daily < 212) then ; time_lev_monthly = 6
  elseif (time_lev_daily < 243) then ; time_lev_monthly = 7
  elseif (time_lev_daily < 273) then ; time_lev_monthly = 8
  elseif (time_lev_daily < 304) then ; time_lev_monthly = 9
  elseif (time_lev_daily < 334) then ; time_lev_monthly = 10
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
      call MOM_read_vector(filename, CS%stress_x_var, CS%stress_y_var, &
                           temp_x(:,:), temp_y(:,:), G%Domain, stagger=AGRID, &
                           timelevel=time_lev, scale=Pa_conversion)

      call pass_vector(temp_x, temp_y, G%Domain, To_All, AGRID)
      do j=js,je ; do I=is-1,Ieq
        forces%taux(I,j) = 0.5 * CS%wind_scale * (temp_x(i,j) + temp_x(i+1,j))
      enddo ; enddo
      do J=js-1,Jeq ; do i=is,ie
        forces%tauy(i,J) = 0.5 * CS%wind_scale * (temp_y(i,j) + temp_y(i,j+1))
      enddo ; enddo

      if (.not.read_Ustar) then
        if (CS%read_gust_2d) then
          do j=js,je ; do i=is,ie
            forces%ustar(i,j) = sqrt((CS%gust(i,j) + &
                    sqrt(temp_x(i,j)*temp_x(i,j) + temp_y(i,j)*temp_y(i,j))) * US%L_to_Z / CS%Rho0)
          enddo ; enddo
        else
          do j=js,je ; do i=is,ie
            forces%ustar(i,j) = sqrt(US%L_to_Z * (CS%gust_const/CS%Rho0 + &
                    sqrt(temp_x(i,j)*temp_x(i,j) + temp_y(i,j)*temp_y(i,j)) / CS%Rho0) )
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
        call MOM_read_vector(filename, CS%stress_x_var, CS%stress_y_var, &
                             temp_x(:,:), temp_y(:,:), &
                             G%Domain_aux, stagger=CGRID_NE, timelevel=time_lev, &
                             scale=Pa_conversion)
        do j=js,je ; do i=is,ie
          forces%taux(I,j) = CS%wind_scale * temp_x(I,j)
          forces%tauy(i,J) = CS%wind_scale * temp_y(i,J)
        enddo ; enddo
        call fill_symmetric_edges(forces%taux, forces%tauy, G%Domain, stagger=CGRID_NE)
      else
        call MOM_read_vector(filename, CS%stress_x_var, CS%stress_y_var, &
                             forces%taux(:,:), forces%tauy(:,:), &
                             G%Domain, stagger=CGRID_NE, timelevel=time_lev, &
                             scale=Pa_conversion)

        if (CS%wind_scale /= 1.0) then
          do j=js,je ; do I=Isq,Ieq
            forces%taux(I,j) = CS%wind_scale * forces%taux(I,j)
          enddo ; enddo
          do J=Jsq,Jeq ; do i=is,ie
            forces%tauy(i,J) = CS%wind_scale * forces%tauy(i,J)
          enddo ; enddo
        endif
      endif

      call pass_vector(forces%taux, forces%tauy, G%Domain, To_All)
      if (.not.read_Ustar) then
        if (CS%read_gust_2d) then
          do j=js, je ; do i=is, ie
            forces%ustar(i,j) = sqrt((CS%gust(i,j) + &
                    sqrt(0.5*((forces%tauy(i,j-1)**2 + forces%tauy(i,j)**2) + &
                              (forces%taux(i-1,j)**2 + forces%taux(i,j)**2))) ) * US%L_to_Z / CS%Rho0 )
          enddo ; enddo
        else
          do j=js, je ; do i=is, ie
            forces%ustar(i,j) = sqrt(US%L_to_Z * ( (CS%gust_const/CS%Rho0) + &
                    sqrt(0.5*((forces%tauy(i,j-1)**2 + forces%tauy(i,j)**2) + &
                              (forces%taux(i-1,j)**2 + forces%taux(i,j)**2)))/CS%Rho0))
          enddo ; enddo
        endif
      endif
    case default
      call MOM_error(FATAL, "wind_forcing_from_file: Unrecognized stagger "//&
                      trim(CS%wind_stagger)//" is not 'A' or 'C'.")
    end select

    if (read_Ustar) then
      call MOM_read_data(filename, CS%Ustar_var, forces%ustar(:,:), &
                         G%Domain, timelevel=time_lev, scale=US%m_to_Z*US%T_to_s)
    endif

    CS%wind_last_lev = time_lev

  endif ! time_lev /= CS%wind_last_lev

  call callTree_leave("wind_forcing_from_file")
end subroutine wind_forcing_from_file


! Sets the surface wind stresses via the data override facility.
subroutine wind_forcing_by_data_override(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                       !! describe the surface state of the ocean.
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  type(ocean_grid_type),    intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control struct returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables
  real :: temp_x(SZI_(G),SZJ_(G)) ! Pseudo-zonal and psuedo-meridional
  real :: temp_y(SZI_(G),SZJ_(G)) ! wind stresses at h-points [Pa].
  real :: temp_ustar(SZI_(G),SZJ_(G)) ! ustar [m s-1] (not rescaled).
  real :: Pa_conversion ! A unit conversion factor from Pa to the internal units [R Z L T-2 Pa-1 ~> 1]
  integer :: i, j, is_in, ie_in, js_in, je_in
  logical :: read_uStar

  call callTree_enter("wind_forcing_by_data_override, MOM_surface_forcing.F90")

  if (.not.CS%dataOverrideIsInitialized) then
    call allocate_mech_forcing(G, forces, stress=.true., ustar=.true., press=.true.)
    call data_override_init(Ocean_domain_in=G%Domain%mpp_domain)
    CS%dataOverrideIsInitialized = .True.
  endif

  is_in = G%isc - G%isd + 1 ; ie_in = G%iec - G%isd + 1
  js_in = G%jsc - G%jsd + 1 ; je_in = G%jec - G%jsd + 1
  Pa_conversion = US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z

  temp_x(:,:) = 0.0 ; temp_y(:,:) = 0.0
  call data_override('OCN', 'taux', temp_x, day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
  call data_override('OCN', 'tauy', temp_y, day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
  call pass_vector(temp_x, temp_y, G%Domain, To_All, AGRID)
  ! Ignore CS%wind_scale when using data_override ?????
  do j=G%jsc,G%jec ; do I=G%isc-1,G%IecB
    forces%taux(I,j) = Pa_conversion * 0.5 * (temp_x(i,j) + temp_x(i+1,j))
  enddo ; enddo
  do J=G%jsc-1,G%JecB ; do i=G%isc,G%iec
    forces%tauy(i,J) = Pa_conversion * 0.5 * (temp_y(i,j) + temp_y(i,j+1))
  enddo ; enddo

  read_Ustar = (len_trim(CS%ustar_var) > 0) ! Need better control higher up ????
  if (read_Ustar) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec ; temp_ustar(i,j) = US%Z_to_m*US%s_to_T*forces%ustar(i,j) ; enddo ; enddo
    call data_override('OCN', 'ustar', temp_ustar, day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
    do j=G%jsc,G%jec ; do i=G%isc,G%iec ; forces%ustar(i,j) = US%m_to_Z*US%T_to_s*temp_ustar(i,j) ; enddo ; enddo
  else
    if (CS%read_gust_2d) then
      call data_override('OCN', 'gust', CS%gust, day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        forces%ustar(i,j) = sqrt((Pa_conversion * sqrt(temp_x(i,j)*temp_x(i,j) + &
            temp_y(i,j)*temp_y(i,j)) + CS%gust(i,j)) * US%L_to_Z / CS%Rho0)
      enddo ; enddo
    else
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        forces%ustar(i,j) = sqrt(US%L_to_Z * (Pa_conversion*sqrt(temp_x(i,j)*temp_x(i,j) + &
            temp_y(i,j)*temp_y(i,j))/CS%Rho0 + CS%gust_const/CS%Rho0 ))
      enddo ; enddo
    endif
  endif

  call pass_vector(forces%taux, forces%tauy, G%Domain, To_All)
! call pass_var(forces%ustar, G%Domain, To_All)     Not needed  ?????

  call callTree_leave("wind_forcing_by_data_override")
end subroutine wind_forcing_by_data_override


!> Specifies zero surface bouyancy fluxes from input files.
subroutine buoyancy_forcing_from_files(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),         intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),       intent(in)    :: day  !< The time of the fluxes
  real,                  intent(in)    :: dt   !< The amount of time over which
                                               !! the fluxes apply [s]
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer    :: CS   !< pointer to control struct returned by
                                               !! a previous surface_forcing_init call
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    temp, &       ! A 2-d temporary work array with various units.
    SST_anom, &   ! Instantaneous sea surface temperature anomalies from a
                  ! target (observed) value [degC].
    SSS_anom, &   ! Instantaneous sea surface salinity anomalies from a target
                  ! (observed) value [ppt].
    SSS_mean      ! A (mean?) salinity about which to normalize local salinity
                  ! anomalies when calculating restorative precipitation
                  ! anomalies [ppt].

  real :: kg_m2_s_conversion  ! A combination of unit conversion factors for rescaling
                              ! mass fluxes [R Z s m2 kg-1 T-1 ~> 1].
  real :: rhoXcp ! reference density times heat capacity [Q R degC-1 ~> J m-3 degC-1]

  integer :: time_lev_daily     ! time levels to read for fields with daily cycle
  integer :: time_lev_monthly   ! time levels to read for fields with monthly cycle
  integer :: time_lev           ! time level that for a field

  integer :: days, seconds
  integer :: i, j, is, ie, js, je

  call callTree_enter("buoyancy_forcing_from_files, MOM_surface_forcing.F90")

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je = G%jec
  kg_m2_s_conversion = US%kg_m2s_to_RZ_T

  if (CS%use_temperature) rhoXcp = CS%Rho0 * fluxes%C_p

  ! Read the buoyancy forcing file
  call get_time(day, seconds, days)

  time_lev_daily = days - 365*floor(real(days) / 365.0)

  if (time_lev_daily < 31) then ; time_lev_monthly = 0
  elseif (time_lev_daily < 59)  then ; time_lev_monthly = 1
  elseif (time_lev_daily < 90)  then ; time_lev_monthly = 2
  elseif (time_lev_daily < 120) then ; time_lev_monthly = 3
  elseif (time_lev_daily < 151) then ; time_lev_monthly = 4
  elseif (time_lev_daily < 181) then ; time_lev_monthly = 5
  elseif (time_lev_daily < 212) then ; time_lev_monthly = 6
  elseif (time_lev_daily < 243) then ; time_lev_monthly = 7
  elseif (time_lev_daily < 273) then ; time_lev_monthly = 8
  elseif (time_lev_daily < 304) then ; time_lev_monthly = 9
  elseif (time_lev_daily < 334) then ; time_lev_monthly = 10
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
    call MOM_read_data(CS%longwave_file, CS%LW_var, fluxes%lw(:,:), &
                       G%Domain, timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
    if (CS%archaic_OMIP_file) then
      call MOM_read_data(CS%longwaveup_file, "lwup_sfc", temp(:,:), G%Domain, &
                         timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
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
      call MOM_read_data(CS%evaporation_file, CS%evap_var, fluxes%evap(:,:), &
                     G%Domain, timelevel=time_lev, scale=-kg_m2_s_conversion)
      do j=js,je ; do i=is,ie
        fluxes%latent(i,j)           = CS%latent_heat_vapor*fluxes%evap(i,j)
        fluxes%latent_evap_diag(i,j) = fluxes%latent(i,j)
      enddo ; enddo
    else
      call MOM_read_data(CS%evaporation_file, CS%evap_var, fluxes%evap(:,:), &
                     G%Domain, timelevel=time_lev, scale=kg_m2_s_conversion)
    endif
    CS%evap_last_lev = time_lev

    select case (CS%latent_nlev)
      case (12)    ; time_lev = time_lev_monthly
      case (365)   ; time_lev = time_lev_daily
      case default ; time_lev = 1
    end select
    if (.not.CS%archaic_OMIP_file) then
      call MOM_read_data(CS%latentheat_file, CS%latent_var, fluxes%latent(:,:), &
                     G%Domain, timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
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
      call MOM_read_data(CS%sensibleheat_file, CS%sens_var, fluxes%sens(:,:), &
                     G%Domain, timelevel=time_lev, scale=-US%W_m2_to_QRZ_T)
    else
      call MOM_read_data(CS%sensibleheat_file, CS%sens_var, fluxes%sens(:,:), &
                     G%Domain, timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
    endif
    CS%sens_last_lev = time_lev

    select case (CS%SW_nlev)
      case (12)    ; time_lev = time_lev_monthly
      case (365)   ; time_lev = time_lev_daily
      case default ; time_lev = 1
    end select
    call MOM_read_data(CS%shortwave_file, CS%SW_var, fluxes%sw(:,:), G%Domain, &
                       timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
    if (CS%archaic_OMIP_file) then
      call MOM_read_data(CS%shortwaveup_file, "swup_sfc", temp(:,:), G%Domain, &
                         timelevel=time_lev, scale=US%W_m2_to_QRZ_T)
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
    call MOM_read_data(CS%snow_file, CS%snow_var, &
             fluxes%fprec(:,:), G%Domain, timelevel=time_lev, scale=kg_m2_s_conversion)
    call MOM_read_data(CS%rain_file, CS%rain_var, &
             fluxes%lprec(:,:), G%Domain, timelevel=time_lev, scale=kg_m2_s_conversion)
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
      call MOM_read_data(CS%runoff_file, CS%lrunoff_var, temp(:,:), &
                     G%Domain, timelevel=time_lev, scale=kg_m2_s_conversion)
      do j=js,je ; do i=is,ie
        fluxes%lrunoff(i,j) = temp(i,j)*US%m_to_L**2*G%IareaT(i,j)
      enddo ; enddo
      call MOM_read_data(CS%runoff_file, CS%frunoff_var, temp(:,:), &
                     G%Domain, timelevel=time_lev, scale=kg_m2_s_conversion)
      do j=js,je ; do i=is,ie
        fluxes%frunoff(i,j) = temp(i,j)*US%m_to_L**2*G%IareaT(i,j)
      enddo ; enddo
    else
      call MOM_read_data(CS%runoff_file, CS%lrunoff_var, fluxes%lrunoff(:,:), &
                     G%Domain, timelevel=time_lev, scale=kg_m2_s_conversion)
      call MOM_read_data(CS%runoff_file, CS%frunoff_var, fluxes%frunoff(:,:), &
                     G%Domain, timelevel=time_lev, scale=kg_m2_s_conversion)
    endif
    CS%runoff_last_lev = time_lev

!     Read the SST and SSS fields for damping.
    if (CS%restorebuoy) then !#CTRL# .or. associated(CS%ctrl_forcing_CSp)) then
      select case (CS%SST_nlev)
        case (12)    ; time_lev = time_lev_monthly
        case (365)   ; time_lev = time_lev_daily
        case default ; time_lev = 1
      end select
      call MOM_read_data(CS%SSTrestore_file, CS%SST_restore_var, &
               CS%T_Restore(:,:), G%Domain, timelevel=time_lev)
      CS%SST_last_lev = time_lev

      select case (CS%SSS_nlev)
        case (12)    ; time_lev = time_lev_monthly
        case (365)   ; time_lev = time_lev_daily
        case default ; time_lev = 1
      end select
      call MOM_read_data(CS%salinityrestore_file, CS%SSS_restore_var, &
               CS%S_Restore(:,:), G%Domain, timelevel=time_lev)
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
      fluxes%lw(i,j)      = fluxes%lw(i,j)      * G%mask2dT(i,j)
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
              ((CS%T_Restore(i,j) - sfc_state%SST(i,j)) * rhoXcp * CS%Flux_const_T)
          fluxes%vprec(i,j) = - (CS%Rho0*CS%Flux_const_S) * &
              (CS%S_Restore(i,j) - sfc_state%SSS(i,j)) / &
              (0.5*(sfc_state%SSS(i,j) + CS%S_Restore(i,j)))
        else
          fluxes%heat_added(i,j) = 0.0
          fluxes%vprec(i,j)      = 0.0
        endif
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0) then
          fluxes%buoy(i,j) = (CS%Dens_Restore(i,j) - sfc_state%sfc_density(i,j)) * &
                             (CS%G_Earth * CS%Flux_const / CS%Rho0)
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

!#CTRL# if (associated(CS%ctrl_forcing_CSp)) then
!#CTRL#   do j=js,je ; do i=is,ie
!#CTRL#     SST_anom(i,j) = sfc_state%SST(i,j) - CS%T_Restore(i,j)
!#CTRL#     SSS_anom(i,j) = sfc_state%SSS(i,j) - CS%S_Restore(i,j)
!#CTRL#     SSS_mean(i,j) = 0.5*(sfc_state%SSS(i,j) + CS%S_Restore(i,j))
!#CTRL#   enddo ; enddo
!#CTRL#   call apply_ctrl_forcing(SST_anom, SSS_anom, SSS_mean, fluxes%heat_added, &
!#CTRL#                           fluxes%vprec, day, dt, G, CS%ctrl_forcing_CSp)
!#CTRL# endif

  call callTree_leave("buoyancy_forcing_from_files")
end subroutine buoyancy_forcing_from_files

!> Specifies zero surface bouyancy fluxes from data over-ride.
subroutine buoyancy_forcing_from_data_override(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),            intent(inout) :: sfc_state !< A structure containing fields that
                                                  !! describe the surface state of the ocean.
  type(forcing),            intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),          intent(in)    :: day  !< The time of the fluxes
  real,                     intent(in)    :: dt   !< The amount of time over which
                                                  !! the fluxes apply [s]
  type(ocean_grid_type),    intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS   !< pointer to control struct returned by
                                                  !! a previous surface_forcing_init call
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    temp, &       ! A 2-d temporary work array with various units.
    SST_anom, &   ! Instantaneous sea surface temperature anomalies from a
                  ! target (observed) value [degC].
    SSS_anom, &   ! Instantaneous sea surface salinity anomalies from a target
                  ! (observed) value [ppt].
    SSS_mean      ! A (mean?) salinity about which to normalize local salinity
                  ! anomalies when calculating restorative precipitation
                  ! anomalies [ppt].
  real :: kg_m2_s_conversion  ! A combination of unit conversion factors for rescaling
                              ! mass fluxes [R Z s m2 kg-1 T-1 ~> 1].
  real :: rhoXcp ! The mean density times the heat capacity [Q R degC-1 ~> J m-3 degC-1].

  integer :: time_lev_daily     ! The time levels to read for fields with
  integer :: time_lev_monthly   ! daily and montly cycles.
  integer :: itime_lev           ! The time level that is used for a field.

  integer :: days, seconds
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  integer :: is_in, ie_in, js_in, je_in

  call callTree_enter("buoyancy_forcing_from_data_override, MOM_surface_forcing.F90")

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  kg_m2_s_conversion = US%kg_m2s_to_RZ_T

  if (CS%use_temperature) rhoXcp = CS%Rho0 * fluxes%C_p

  if (.not.CS%dataOverrideIsInitialized) then
    call data_override_init(Ocean_domain_in=G%Domain%mpp_domain)
    CS%dataOverrideIsInitialized = .True.
  endif

  is_in = G%isc - G%isd + 1
  ie_in = G%iec - G%isd + 1
  js_in = G%jsc - G%jsd + 1
  je_in = G%jec - G%jsd + 1

  call data_override('OCN', 'lw', fluxes%lw(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in) ! scale=US%W_m2_to_QRZ_T
  if (US%QRZ_T_to_W_m2 /= 1.0) then ; do j=js,je ; do i=is,ie
    fluxes%lw(i,j) = fluxes%lw(i,j) * US%W_m2_to_QRZ_T
  enddo ; enddo ; endif
  call data_override('OCN', 'evap', fluxes%evap(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)

  ! note the sign convention
  do j=js,je ; do i=is,ie
    ! The normal convention is that fluxes%evap positive into the ocean
    ! but evap is normally a positive quantity in the files
    ! This conversion is dangerous because it is not clear whether the data files have been read!
    fluxes%evap(i,j) = -kg_m2_s_conversion*fluxes%evap(i,j)
    fluxes%latent(i,j)           = CS%latent_heat_vapor*fluxes%evap(i,j)
    fluxes%latent_evap_diag(i,j) = fluxes%latent(i,j)
  enddo ; enddo

  call data_override('OCN', 'sens', fluxes%sens(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)

  ! note the sign convention
  do j=js,je ; do i=is,ie
     fluxes%sens(i,j) = -US%W_m2_to_QRZ_T * fluxes%sens(i,j)  ! Normal convention is positive into the ocean
                                           ! but sensible is normally a positive quantity in the files
  enddo ; enddo

  call data_override('OCN', 'sw', fluxes%sw(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in) ! scale=US%W_m2_to_QRZ_T
  if (US%QRZ_T_to_W_m2 /= 1.0) then ; do j=js,je ; do i=is,ie
    fluxes%sw(i,j) = fluxes%sw(i,j) * US%W_m2_to_QRZ_T
  enddo ; enddo ; endif

  call data_override('OCN', 'snow', fluxes%fprec(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in) ! scale=kg_m2_s_conversion

  call data_override('OCN', 'rain', fluxes%lprec(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in) ! scale=kg_m2_s_conversion

  call data_override('OCN', 'runoff', fluxes%lrunoff(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in) ! scale=kg_m2_s_conversion

  call data_override('OCN', 'calving', fluxes%frunoff(:,:), day, &
       is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in) ! scale=kg_m2_s_conversion

  if (kg_m2_s_conversion /= 1.0) then ; do j=js,je ; do i=is,ie
    fluxes%lprec(i,j) = fluxes%lprec(i,j) * kg_m2_s_conversion
    fluxes%fprec(i,j) = fluxes%fprec(i,j) * kg_m2_s_conversion
    fluxes%lrunoff(i,j) = fluxes%lrunoff(i,j) * kg_m2_s_conversion
    fluxes%frunoff(i,j) = fluxes%frunoff(i,j) * kg_m2_s_conversion
  enddo ; enddo ; endif

!     Read the SST and SSS fields for damping.
  if (CS%restorebuoy) then !#CTRL# .or. associated(CS%ctrl_forcing_CSp)) then
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
              ((CS%T_Restore(i,j) - sfc_state%SST(i,j)) * rhoXcp * CS%Flux_const_T)
          fluxes%vprec(i,j) = - (CS%Rho0*CS%Flux_const_S) * &
              (CS%S_Restore(i,j) - sfc_state%SSS(i,j)) / &
              (0.5*(sfc_state%SSS(i,j) + CS%S_Restore(i,j)))
        else
          fluxes%heat_added(i,j) = 0.0
          fluxes%vprec(i,j)      = 0.0
        endif
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0) then
          fluxes%buoy(i,j) = (CS%Dens_Restore(i,j) - sfc_state%sfc_density(i,j)) * &
                             (CS%G_Earth * CS%Flux_const / CS%Rho0)
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
    fluxes%lw(i,j)      = fluxes%lw(i,j)      * G%mask2dT(i,j)
    fluxes%latent(i,j)  = fluxes%latent(i,j)  * G%mask2dT(i,j)
    fluxes%sens(i,j)    = fluxes%sens(i,j)    * G%mask2dT(i,j)
    fluxes%sw(i,j)      = fluxes%sw(i,j)      * G%mask2dT(i,j)

    fluxes%latent_evap_diag(i,j)     = fluxes%latent_evap_diag(i,j) * G%mask2dT(i,j)
    fluxes%latent_fprec_diag(i,j)    = -fluxes%fprec(i,j)*CS%latent_heat_fusion
    fluxes%latent_frunoff_diag(i,j)  = -fluxes%frunoff(i,j)*CS%latent_heat_fusion
  enddo ; enddo


!#CTRL# if (associated(CS%ctrl_forcing_CSp)) then
!#CTRL#   do j=js,je ; do i=is,ie
!#CTRL#     SST_anom(i,j) = sfc_state%SST(i,j) - CS%T_Restore(i,j)
!#CTRL#     SSS_anom(i,j) = sfc_state%SSS(i,j) - CS%S_Restore(i,j)
!#CTRL#     SSS_mean(i,j) = 0.5*(sfc_state%SSS(i,j) + CS%S_Restore(i,j))
!#CTRL#   enddo ; enddo
!#CTRL#   call apply_ctrl_forcing(SST_anom, SSS_anom, SSS_mean, fluxes%heat_added, &
!#CTRL#                           fluxes%vprec, day, dt, G, CS%ctrl_forcing_CSp)
!#CTRL# endif

  call callTree_leave("buoyancy_forcing_from_data_override")
end subroutine buoyancy_forcing_from_data_override

!> This subroutine specifies zero surface bouyancy fluxes
subroutine buoyancy_forcing_zero(sfc_state, fluxes, day, dt, G, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),         intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),       intent(in)    :: day  !< The time of the fluxes
  real,                  intent(in)    :: dt   !< The amount of time over which
                                               !! the fluxes apply [s]
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(surface_forcing_CS), pointer    :: CS   !< pointer to control struct returned by
                                               !! a previous surface_forcing_init call
  ! Local variables
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


!> Sets up spatially and temporally constant surface heat fluxes.
subroutine buoyancy_forcing_const(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),         intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),       intent(in)    :: day  !< The time of the fluxes
  real,                  intent(in)    :: dt   !< The amount of time over which
                                               !! the fluxes apply [s]
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer    :: CS   !< pointer to control struct returned by
                                               !! a previous surface_forcing_init call
  ! Local variables
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

!> Sets surface fluxes of heat and salinity by restoring to temperature and
!! salinity profiles that vary linearly with latitude.
subroutine buoyancy_forcing_linear(sfc_state, fluxes, day, dt, G, US, CS)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),         intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
  type(time_type),       intent(in)    :: day  !< The time of the fluxes
  real,                  intent(in)    :: dt   !< The amount of time over which
                                               !! the fluxes apply [s]
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer    :: CS   !< pointer to control struct returned by
                                               !! a previous surface_forcing_init call
  ! Local variables
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
              ((T_Restore - sfc_state%SST(i,j)) * ((CS%Rho0 * fluxes%C_p) * CS%Flux_const))
          fluxes%vprec(i,j) = - (CS%Rho0*CS%Flux_const) * &
              (S_Restore - sfc_state%SSS(i,j)) / &
              (0.5*(sfc_state%SSS(i,j) + S_Restore))
        else
          fluxes%heat_added(i,j) = 0.0
          fluxes%vprec(i,j)      = 0.0
        endif
      enddo ; enddo
    else
      call MOM_error(FATAL, "buoyancy_forcing_linear in MOM_surface_forcing: "// &
                     "RESTOREBUOY to linear not written yet.")
     !do j=js,je ; do i=is,ie
     !  if (G%mask2dT(i,j) > 0) then
     !    fluxes%buoy(i,j) = (CS%Dens_Restore(i,j) - sfc_state%sfc_density(i,j)) * &
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

!> Save a restart file for the forcing fields
subroutine forcing_save_restart(CS, G, Time, directory, time_stamped, &
                                filename_suffix)
  type(surface_forcing_CS),   pointer       :: CS   !< pointer to control struct returned by
                                                    !! a previous surface_forcing_init call
  type(ocean_grid_type),      intent(inout) :: G    !< The ocean's grid structure
  type(time_type),            intent(in)    :: Time !< model time at this call; needed for mpp_write calls
  character(len=*),           intent(in)    :: directory !< directory into which to write these restart files
  logical,          optional, intent(in)    :: time_stamped !< If true, the restart file names
                                                    !! include a unique time stamp; the  default is false.
  character(len=*), optional, intent(in)    :: filename_suffix !< optional suffix (e.g., a time-stamp)
                                                    !! to append to the restart fname

  if (.not.associated(CS)) return
  if (.not.associated(CS%restart_CSp)) return

  call save_restart(directory, Time, G, CS%restart_CSp, time_stamped)

end subroutine forcing_save_restart

!> Initialize the surface forcing module
subroutine surface_forcing_init(Time, G, US, param_file, diag, CS, tracer_flow_CSp)
  type(time_type),              intent(in)    :: Time !< The current model time
  type(ocean_grid_type),        intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),        intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),        intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,      intent(inout) :: diag !< structure used to regulate diagnostic output
  type(surface_forcing_CS),     pointer       :: CS   !< pointer to control struct returned by
                                                      !! a previous surface_forcing_init call
  type(tracer_flow_control_CS), pointer       :: tracer_flow_CSp !< Forcing for tracers?

  ! Local variables
  type(directories)  :: dirs
  logical            :: new_sim
  type(time_type)    :: Time_frc
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  real :: flux_const_default ! The unscaled value of FLUXCONST [m day-1]
  logical :: default_2018_answers
  character(len=40)  :: mdl = "MOM_surface_forcing" ! This module's name.
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
  call log_version(param_file, mdl, version, '')
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
                 "(linear), (USER), (BFB) and (NONE).", default="zero")
  if (trim(CS%buoy_config) == "file") then
    call get_param(param_file, mdl, "ARCHAIC_OMIP_FORCING_FILE", CS%archaic_OMIP_file, &
                 "If true, use the forcing variable decomposition from "//&
                 "the old German OMIP prescription that predated CORE. If "//&
                 "false, use the variable groupings available from MOM "//&
                 "output diagnostics of forcing variables.", default=.true.)
    if (CS%archaic_OMIP_file) then
      call get_param(param_file, mdl, "LONGWAVEDOWN_FILE", CS%longwave_file, &
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
      call get_param(param_file, mdl, "SHORTWAVEDOWN_FILE", CS%shortwave_file, &
                 "The file with the downward shortwave heat flux.", &
                 fail_if_missing=.true.)
      call get_param(param_file, mdl, "SNOW_FILE", CS%snow_file, &
                 "The file with the downward frozen precip flux, in "//&
                 "variable snow.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "PRECIP_FILE", CS%rain_file, &
                 "The file with the downward total precip flux, in "//&
                 "variable precip.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "FRESHDISCHARGE_FILE", CS%runoff_file, &
                 "The file with the fresh and frozen runoff/calving fluxes, "//&
                 "invariables disch_w and disch_s.", fail_if_missing=.true.)

      ! These variable names are hard-coded, per the archaic OMIP conventions.
      CS%latentheat_file = CS%evaporation_file ; CS%latent_var = "evap"
      CS%LW_var = "lwdn_sfc"; CS%SW_var = "swdn_sfc"; CS%sens_var = "shflx"
      CS%evap_var = "evap"; CS%rain_var = "precip"; CS%snow_var = "snow"
      CS%lrunoff_var = "disch_w"; CS%frunoff_var = "disch_s"

    else
      call get_param(param_file, mdl, "LONGWAVE_FILE", CS%longwave_file, &
                 "The file with the longwave heat flux, in the variable "//&
                 "given by LONGWAVE_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "LONGWAVE_FORCING_VAR", CS%LW_var, &
                 "The variable with the longwave forcing field.", default="LW")

      call get_param(param_file, mdl, "SHORTWAVE_FILE", CS%shortwave_file, &
                 "The file with the shortwave heat flux, in the variable "//&
                 "given by SHORTWAVE_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "SHORTWAVE_FORCING_VAR", CS%SW_var, &
                 "The variable with the shortwave forcing field.", default="SW")

      call get_param(param_file, mdl, "EVAPORATION_FILE", CS%evaporation_file, &
                 "The file with the evaporative moisture flux, in the "//&
                 "variable given by EVAP_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "EVAP_FORCING_VAR", CS%evap_var, &
                 "The variable with the evaporative moisture flux.", &
                 default="evap")

      call get_param(param_file, mdl, "LATENTHEAT_FILE", CS%latentheat_file, &
                 "The file with the latent heat flux, in the variable "//&
                 "given by LATENT_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "LATENT_FORCING_VAR", CS%latent_var, &
                 "The variable with the latent heat flux.", default="latent")

      call get_param(param_file, mdl, "SENSIBLEHEAT_FILE", CS%sensibleheat_file, &
                 "The file with the sensible heat flux, in the variable "//&
                 "given by SENSIBLE_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "SENSIBLE_FORCING_VAR", CS%sens_var, &
                 "The variable with the sensible heat flux.", default="sensible")

      call get_param(param_file, mdl, "RAIN_FILE", CS%rain_file, &
                 "The file with the liquid precipitation flux, in the "//&
                 "variable given by RAIN_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "RAIN_FORCING_VAR", CS%rain_var, &
                 "The variable with the liquid precipitation flux.", &
                 default="liq_precip")
      call get_param(param_file, mdl, "SNOW_FILE", CS%snow_file, &
                 "The file with the frozen precipitation flux, in the "//&
                 "variable given by SNOW_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "SNOW_FORCING_VAR", CS%snow_var, &
                 "The variable with the frozen precipitation flux.", &
                 default="froz_precip")

      call get_param(param_file, mdl, "RUNOFF_FILE", CS%runoff_file, &
                 "The file with the fresh and frozen runoff/calving "//&
                 "fluxes, in variables given by LIQ_RUNOFF_FORCING_VAR "//&
                 "and FROZ_RUNOFF_FORCING_VAR.", fail_if_missing=.true.)
      call get_param(param_file, mdl, "LIQ_RUNOFF_FORCING_VAR", CS%lrunoff_var, &
                 "The variable with the liquid runoff flux.", &
                 default="liq_runoff")
      call get_param(param_file, mdl, "FROZ_RUNOFF_FORCING_VAR", CS%frunoff_var, &
                 "The variable with the frozen runoff flux.", &
                 default="froz_runoff")
    endif

    call get_param(param_file, mdl, "SSTRESTORE_FILE", CS%SSTrestore_file, &
                 "The file with the SST toward which to restore in the "//&
                 "variable given by SST_RESTORE_VAR.", fail_if_missing=.true.)
    call get_param(param_file, mdl, "SALINITYRESTORE_FILE", CS%salinityrestore_file, &
                 "The file with the surface salinity toward which to "//&
                 "restore in the variable given by SSS_RESTORE_VAR.", &
                 fail_if_missing=.true.)

    if (CS%archaic_OMIP_file) then
      CS%SST_restore_var = "TEMP" ; CS%SSS_restore_var = "SALT"
    else
      call get_param(param_file, mdl, "SST_RESTORE_VAR", CS%SST_restore_var, &
                 "The variable with the SST toward which to restore.", &
                 default="SST")
      call get_param(param_file, mdl, "SSS_RESTORE_VAR", CS%SSS_restore_var, &
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
    call get_param(param_file, mdl, "SENSIBLE_HEAT_FLUX", CS%constantHeatForcing, &
                 "A constant heat forcing (positive into ocean) applied "//&
                 "through the sensible heat flux field. ", &
                 units='W/m2', scale=US%W_m2_to_QRZ_T, fail_if_missing=.true.)
  endif
  call get_param(param_file, mdl, "WIND_CONFIG", CS%wind_config, &
                 "The character string that indicates how wind forcing "//&
                 "is specified. Valid options include (file), (2gyre), "//&
                 "(1gyre), (gyres), (zero), and (USER).", default="zero")
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
    call get_param(param_file, mdl, "WIND_STAGGER",CS%wind_stagger, &
                 "A character indicating how the wind stress components "//&
                 "are staggered in WIND_FILE.  This may be A or C for now.", &
                 default="C")
    call get_param(param_file, mdl, "WINDSTRESS_SCALE", CS%wind_scale, &
                 "A value by which the wind stresses in WIND_FILE are rescaled.", &
                 default=1.0, units="nondim")
    call get_param(param_file, mdl, "USTAR_FORCING_VAR", CS%ustar_var, &
                 "The name of the friction velocity variable in WIND_FILE "//&
                 "or blank to get ustar from the wind stresses plus the "//&
                 "gustiness.", default=" ", units="nondim")
    CS%wind_file = trim(CS%inputdir) // trim(CS%wind_file)
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
    call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
    call get_param(param_file, mdl, "WIND_GYRES_2018_ANSWERS", CS%answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the answers "//&
                 "from the end of 2018.  Otherwise, use expressions for the gyre friction velocities "//&
                 "that are rotationally invariant and more likely to be the same between compilers.", &
                 default=default_2018_answers)
  else
    CS%answers_2018 = .false.
  endif
  if (trim(CS%wind_config) == "scurves") then
    call get_param(param_file, mdl, "WIND_SCURVES_LATS", CS%scurves_ydata, &
                 "A list of latitudes defining a piecewise scurve profile "//&
                 "for zonal wind stress.", &
                 units="degrees N", fail_if_missing=.true.)
    call get_param(param_file, mdl, "WIND_SCURVES_TAUX", CS%scurves_taux, &
                 "A list of zonal wind stress values at latitudes "//&
                 "WIND_SCURVES_LATS defining a piecewise scurve profile.", &
                 units="Pa", fail_if_missing=.true.)
  endif
  if ((trim(CS%wind_config) == "2gyre") .or. &
      (trim(CS%wind_config) == "1gyre") .or. &
      (trim(CS%wind_config) == "gyres") .or. &
      (trim(CS%buoy_config) == "linear")) then
    CS%south_lat = G%south_lat
    CS%len_lat = G%len_lat
  endif

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
    ! These three variables use non-standard time units, but are rescaled as they are read.
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes to the relative "//&
                 "surface anomalies (akin to a piston velocity).  Note the non-MKS units.", &
                 default=0.0, units="m day-1", scale=US%m_to_Z*US%T_to_s/86400.0, &
                 unscaled=flux_const_default)

    if (CS%use_temperature) then
      call get_param(param_file, mdl, "FLUXCONST_T", CS%Flux_const_T, &
           "The constant that relates the restoring surface temperature "//&
           "flux to the relative surface anomaly (akin to a piston "//&
           "velocity).  Note the non-MKS units.", &
           units="m day-1", scale=US%m_to_Z*US%T_to_s/86400.0, &
           default=flux_const_default)
      call get_param(param_file, mdl, "FLUXCONST_S", CS%Flux_const_S, &
           "The constant that relates the restoring surface salinity "//&
           "flux to the relative surface anomaly (akin to a piston "//&
           "velocity).  Note the non-MKS units.", &
           units="m day-1", scale=US%m_to_Z*US%T_to_s/86400.0, &
           default=flux_const_default)
    endif

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
                 units="Pa", default=0.0, scale=US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z)
  call get_param(param_file, mdl, "FIX_USTAR_GUSTLESS_BUG", CS%fix_ustar_gustless_bug, &
                 "If true correct a bug in the time-averaging of the gustless wind friction velocity", &
                 default=.true.)
  call get_param(param_file, mdl, "READ_GUST_2D", CS%read_gust_2d, &
                 "If true, use a 2-dimensional gustiness supplied from "//&
                 "an input file", default=.false.)
  if (CS%read_gust_2d) then
    call get_param(param_file, mdl, "GUST_2D_FILE", gust_file, &
                 "The file in which the wind gustiness is found in "//&
                 "variable gustiness.", fail_if_missing=.true.)
    call safe_alloc_ptr(CS%gust,G%isd,G%ied,G%jsd,G%jed)
    filename = trim(CS%inputdir) // trim(gust_file)
    call MOM_read_data(filename,'gustiness',CS%gust,G%domain, timelevel=1, &
                   scale=US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z) ! units in file should be Pa
  endif

!  All parameter settings are now known.

  if (trim(CS%wind_config) == "USER" .or. trim(CS%buoy_config) == "USER" ) then
    call USER_surface_forcing_init(Time, G, US, param_file, diag, CS%user_forcing_CSp)
  elseif (trim(CS%buoy_config) == "BFB" ) then
    call BFB_surface_forcing_init(Time, G, US, param_file, diag, CS%BFB_forcing_CSp)
  elseif (trim(CS%buoy_config) == "dumbbell" ) then
    call dumbbell_surface_forcing_init(Time, G, US, param_file, diag, CS%dumbbell_forcing_CSp)
  elseif (trim(CS%wind_config) == "MESO" .or. trim(CS%buoy_config) == "MESO" ) then
    call MESO_surface_forcing_init(Time, G, US, param_file, diag, CS%MESO_forcing_CSp)
  elseif (trim(CS%wind_config) == "ideal_hurr" .or.&
          trim(CS%wind_config) == "SCM_ideal_hurr") then
    call idealized_hurricane_wind_init(Time, G, US, param_file, CS%idealized_hurricane_CSp)
  elseif (trim(CS%wind_config) == "const") then
    call get_param(param_file, mdl, "CONST_WIND_TAUX", CS%tau_x0, &
                 "With wind_config const, this is the constant zonal "//&
                 "wind-stress", units="Pa", fail_if_missing=.true.)
    call get_param(param_file, mdl, "CONST_WIND_TAUY", CS%tau_y0, &
                 "With wind_config const, this is the constant meridional "//&
                 "wind-stress", units="Pa", fail_if_missing=.true.)
  elseif (trim(CS%wind_config) == "SCM_CVmix_tests" .or. &
          trim(CS%buoy_config) == "SCM_CVmix_tests") then
    call SCM_CVmix_tests_surface_forcing_init(Time, G, param_file, CS%SCM_CVmix_tests_CSp)
  endif

  call register_forcing_type_diags(Time, diag, US, CS%use_temperature, CS%handles)

  ! Set up any restart fields associated with the forcing.
  call restart_init(param_file, CS%restart_CSp, "MOM_forcing.res")
!#CTRL#  call register_ctrl_forcing_restarts(G, param_file, CS%ctrl_forcing_CSp, &
!#CTRL#                                      CS%restart_CSp)
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

!#CTRL#  call controlled_forcing_init(Time, G, param_file, diag, CS%ctrl_forcing_CSp)

  call user_revise_forcing_init(param_file, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)
end subroutine surface_forcing_init


!> Deallocate memory associated with the surface forcing module
subroutine surface_forcing_end(CS, fluxes)
  type(surface_forcing_CS), pointer    :: CS   !< pointer to control struct returned by
                                               !! a previous surface_forcing_init call
  type(forcing), optional,  intent(inout) :: fluxes !< A structure containing thermodynamic forcing fields
! Arguments:  CS - A pointer to the control structure returned by a previous
!                  call to surface_forcing_init, it will be deallocated here.
!  (inout)    fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.

  if (present(fluxes)) call deallocate_forcing_type(fluxes)

!#CTRL#  call controlled_forcing_end(CS%ctrl_forcing_CSp)

  if (associated(CS)) deallocate(CS)
  CS => NULL()

end subroutine surface_forcing_end

end module MOM_surface_forcing
