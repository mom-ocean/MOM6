!> Initialize state variables, u, v, h, T and S.
module MOM_state_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging, only : hchksum, qchksum, uchksum, vchksum
use MOM_coms, only : max_across_PEs, min_across_PEs
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only :  CLOCK_ROUTINE, CLOCK_LOOP
use MOM_domains, only : pass_var, pass_vector, sum_across_PEs, broadcast
use MOM_domains, only : root_PE, To_All, SCALAR_PAIR, CGRID_NE, AGRID
use MOM_EOS, only : find_depth_of_pressure_in_cell
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, param_file_type
use MOM_file_parser, only : log_version
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type, isPointInCell
use MOM_interface_heights, only : find_eta
use MOM_io, only : close_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE, MULTIPLE
use MOM_io, only : slasher, vardesc, write_field
use MOM_io, only : EAST_FACE, NORTH_FACE
use MOM_open_boundary, only : ocean_OBC_type, open_boundary_init
use MOM_open_boundary, only : OBC_NONE, OBC_SIMPLE
use MOM_open_boundary, only : open_boundary_query, set_Flather_data
!use MOM_open_boundary, only : set_3D_OBC_data
use MOM_grid_initialize, only : initialize_masks, set_grid_metrics
use MOM_restart, only : restore_state, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, set_up_sponge_ML_density
use MOM_sponge, only : initialize_sponge, sponge_CS
use MOM_ALE_sponge, only : set_up_ALE_sponge_field, initialize_ALE_sponge
use MOM_ALE_sponge, only : ALE_sponge_CS
use MOM_string_functions, only : uppercase, lowercase
use MOM_time_manager, only : time_type, set_time
use MOM_tracer_registry, only : add_tracer_OBC_values, tracer_registry_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : setVerticalGridAxes, verticalGrid_type
use MOM_ALE, only : pressure_gradient_plm
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use MOM_EOS, only : int_specific_vol_dp, convert_temp_salt_for_TEOS10
use user_initialization, only : user_initialize_thickness, user_initialize_velocity
use user_initialization, only : user_init_temperature_salinity
use user_initialization, only : user_set_OBC_data
use user_initialization, only : user_initialize_sponges
use DOME_initialization, only : DOME_initialize_thickness
use DOME_initialization, only : DOME_set_OBC_data
use DOME_initialization, only : DOME_initialize_sponges
use ISOMIP_initialization, only : ISOMIP_initialize_thickness
use ISOMIP_initialization, only : ISOMIP_initialize_sponges
use ISOMIP_initialization, only : ISOMIP_initialize_temperature_salinity
use baroclinic_zone_initialization, only : baroclinic_zone_init_temperature_salinity
use benchmark_initialization, only : benchmark_initialize_thickness
use benchmark_initialization, only : benchmark_init_temperature_salinity
use circle_obcs_initialization, only : circle_obcs_initialize_thickness
use lock_exchange_initialization, only : lock_exchange_initialize_thickness
use external_gwave_initialization, only : external_gwave_initialize_thickness
use DOME2d_initialization, only : DOME2d_initialize_thickness
use DOME2d_initialization, only : DOME2d_initialize_temperature_salinity
use DOME2d_initialization, only : DOME2d_initialize_sponges
use adjustment_initialization, only : adjustment_initialize_thickness
use adjustment_initialization, only : adjustment_initialize_temperature_salinity
use sloshing_initialization, only : sloshing_initialize_thickness
use sloshing_initialization, only : sloshing_initialize_temperature_salinity
use seamount_initialization, only : seamount_initialize_thickness
use seamount_initialization, only : seamount_initialize_temperature_salinity
use Phillips_initialization, only : Phillips_initialize_thickness
use Phillips_initialization, only : Phillips_initialize_velocity
use Phillips_initialization, only : Phillips_initialize_sponges
use Rossby_front_2d_initialization, only : Rossby_front_initialize_thickness
use Rossby_front_2d_initialization, only : Rossby_front_initialize_temperature_salinity
use Rossby_front_2d_initialization, only : Rossby_front_initialize_velocity
use SCM_idealized_hurricane, only : SCM_idealized_hurricane_TS_init
use SCM_CVmix_tests, only: SCM_CVmix_tests_TS_init
use supercritical_initialization, only : supercritical_set_OBC_data
use soliton_initialization, only : soliton_initialize_velocity
use soliton_initialization, only : soliton_initialize_thickness
use BFB_initialization, only : BFB_initialize_sponges_southonly

use midas_vertmap, only : find_interfaces, tracer_Z_init
use midas_vertmap, only : determine_temperature

use MOM_ALE, only : ALE_initRegridding, ALE_CS, ALE_initThicknessToCoord
use MOM_ALE, only : ALE_remap_scalar, ALE_build_grid
use MOM_regridding, only : regridding_CS, set_regrid_params, getCoordinateResolution
use MOM_remapping, only : remapping_CS, initialize_remapping
use MOM_remapping, only : remapping_core_h
use MOM_tracer_initialization_from_Z, only : horiz_interp_and_extrap_tracer

implicit none ; private

#include <MOM_memory.h>

public MOM_initialize_state

character(len=40)  :: mod = "MOM_state_initialization" ! This module's name.

contains

! -----------------------------------------------------------------------------
subroutine MOM_initialize_state(u, v, h, tv, Time, G, GV, PF, dirs, &
                                restart_CS, ALE_CSp, tracer_Reg, sponge_CSp, &
                                ALE_sponge_CSp, OBC, Time_in)
  type(ocean_grid_type),                     intent(inout) :: G
  type(verticalGrid_type),                   intent(in)    :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out)   :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out)   :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(out)   :: h
  type(thermo_var_ptrs),                     intent(inout) :: tv
  type(time_type),                           intent(inout) :: Time
  type(param_file_type),                     intent(in)    :: PF
  type(directories),                         intent(in)    :: dirs
  type(MOM_restart_CS),                      pointer       :: restart_CS
  type(ALE_CS),                              pointer       :: ALE_CSp
  type(tracer_registry_type),                pointer       :: tracer_Reg
  type(sponge_CS),                           pointer       :: sponge_CSp
  type(ALE_sponge_CS),                       pointer       :: ALE_sponge_CSp
  type(ocean_OBC_type),                      pointer       :: OBC
  type(time_type), optional,                 intent(in)    :: Time_in
! Arguments: u  - Zonal velocity, in m s-1.
!  (out)     v  - Meridional velocity, in m s-1.
!  (out)     h  - Layer thickness, in m.
!  (out)     tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (out)     Time    - Time at the start of the run segment.
!  (inout)   G       - The ocean's grid structure.
!  (in)      GV      - The ocean's vertical grid structure.
!  (in)      PF      - A structure indicating the open file to parse for
!                      model parameter values.
!  (in)      dirs    - A structure containing several relevant directory paths.
!  (inout)   restart_CS - A pointer to the restart control structure.
!  (inout)   CS      - A structure of pointers to be exchanged with MOM.F90.
!  (in)      Time_in - Time at the start of the run segment. Time_in overrides
!                      any value set for Time.

  character(len=200) :: filename   ! The name of an input file.
  character(len=200) :: filename2  ! The name of an input files.
  character(len=200) :: inputdir   ! The directory where NetCDF input files are.
  character(len=200) :: config
  logical :: from_Z_file, useALE
  logical :: new_sim
  integer :: write_geom
  logical :: use_temperature, use_sponge
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  logical :: depress_sfc ! If true, remove the mass that would be displaced
                         ! by a large surface pressure by squeezing the column.
  logical :: trim_ic_for_p_surf ! If true, remove the mass that would be displaced
                         ! by a large surface pressure, such as with an ice sheet.
  logical :: Analytic_FV_PGF, obsol_test
  logical :: convert
  type(EOS_type), pointer :: eos => NULL()
  logical :: debug    ! indicates whether to write debugging output
! This include declares and sets the variable "version".
#include "version_variable.h"
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call callTree_enter("MOM_initialize_state(), MOM_state_initialization.F90")
  call log_version(PF, mod, version, "")
  call get_param(PF, mod, "DEBUG", debug, default=.false.)

  new_sim = .false.
  if ((dirs%input_filename(1:1) == 'n') .and. &
      (LEN_TRIM(dirs%input_filename) == 1)) new_sim = .true.

  call get_param(PF, mod, "INPUTDIR", inputdir, &
         "The directory in which input files are found.", default=".")
  inputdir = slasher(inputdir)

  use_temperature = ASSOCIATED(tv%T)
  useALE = associated(ALE_CSp)
  use_EOS = associated(tv%eqn_of_state)
  if (use_EOS) eos => tv%eqn_of_state

!====================================================================
!    Initialize temporally evolving fields, either as initial
!  conditions or by reading them from a restart (or saves) file.
!====================================================================

  if (new_sim) then
!  This block initializes all of the fields internally.              !
    call MOM_mesg("Run initialized internally.", 3)

    if (present(Time_in)) Time = Time_in
    ! Otherwise leave Time at its input value.

    call get_param(PF, mod, "INIT_LAYERS_FROM_Z_FILE", from_Z_file, &
               "If true, intialize the layer thicknesses, temperatures, \n"//&
               "and salnities from a Z-space file on a latitude- \n"//&
               "longitude grid.", default=.false.)
    ! h will be converted from m to H below
    h(:,:,:) = GV%Angstrom_z

    if (from_Z_file) then
!     Initialize thickness and T/S from z-coordinate data in a file.
      if (.NOT.use_temperature) call MOM_error(FATAL,"MOM_initialize_state : "//&
         "use_temperature must be true if INIT_LAYERS_FROM_Z_FILE is true")

      call MOM_temp_salt_initialize_from_Z(h, tv, G, GV, PF, dirs)
      call pass_var(h, G%Domain)

    else
!     Initialize thickness, h.
      call get_param(PF, mod, "THICKNESS_CONFIG", config, &
               "A string that determines how the initial layer \n"//&
               "thicknesses are specified for a new run: \n"//&
               " \t file - read interface heights from the file specified \n"//&
               " \t thickness_file - read thicknesses from the file specified \n"//&
               " \t\t by (THICKNESS_FILE).\n"//&
               " \t coord - determined by ALE coordinate.\n"//&
               " \t uniform - uniform thickness layers evenly distributed \n"//&
               " \t\t between the surface and MAXIMUM_DEPTH. \n"//&
               " \t DOME - use a slope and channel configuration for the \n"//&
               " \t\t DOME sill-overflow test case. \n"//&
               " \t ISOMIP - use a configuration for the \n"//&
               " \t\t ISOMIP test case. \n"//&
               " \t benchmark - use the benchmark test case thicknesses. \n"//&
               " \t search - search a density profile for the interface \n"//&
               " \t\t densities. This is not yet implemented. \n"//&
               " \t circle_obcs - the circle_obcs test case is used. \n"//&
               " \t DOME2D - 2D version of DOME initialization. \n"//&
               " \t adjustment2d - TBD AJA. \n"//&
               " \t sloshing - TBD AJA. \n"//&
               " \t seamount - TBD AJA. \n"//&
               " \t soliton - Equatorial Rossby soliton. \n"//&
               " \t rossby_front - a mixed layer front in thermal wind balance.\n"//&
               " \t USER - call a user modified routine.", &
               fail_if_missing=.true.)
      select case (trim(config))
         case ("file"); call initialize_thickness_from_file(h, G, GV, PF, .false.)
         case ("thickness_file"); call initialize_thickness_from_file(h, G, GV, PF, .true.)
         case ("coord")
           if (useALE) then
             call ALE_initThicknessToCoord( ALE_CSp, G, GV, h )
           else
             call MOM_error(FATAL, "MOM_initialize_state: USE_REGRIDDING must be True "//&
                                   "for THICKNESS_CONFIG of 'coord'")
           endif
         case ("uniform"); call initialize_thickness_uniform(h, G, GV, PF)
         case ("DOME"); call DOME_initialize_thickness(h, G, GV, PF)
         case ("ISOMIP"); call ISOMIP_initialize_thickness(h, G, GV, PF, tv)
         case ("benchmark"); call benchmark_initialize_thickness(h, G, GV, PF, &
                                 tv%eqn_of_state, tv%P_Ref)
         case ("search"); call initialize_thickness_search
         case ("circle_obcs"); call circle_obcs_initialize_thickness(h, G, GV, PF)
         case ("lock_exchange"); call lock_exchange_initialize_thickness(h, G, GV, PF)
         case ("external_gwave"); call external_gwave_initialize_thickness(h, G, PF)
         case ("DOME2D"); call DOME2d_initialize_thickness(h, G, GV, PF)
         case ("adjustment2d"); call adjustment_initialize_thickness(h, G, GV, PF)
         case ("sloshing"); call sloshing_initialize_thickness(h, G, GV, PF)
         case ("seamount"); call seamount_initialize_thickness(h, G, GV, PF)
         case ("soliton"); call soliton_initialize_thickness(h, G)
         case ("phillips"); call Phillips_initialize_thickness(h, G, GV, PF)
         case ("rossby_front"); call Rossby_front_initialize_thickness(h, G, GV, PF)
         case ("USER"); call user_initialize_thickness(h, G, PF, tv%T)
         case default ; call MOM_error(FATAL,  "MOM_initialize_state: "//&
              "Unrecognized layer thickness configuration "//trim(config))
      end select
      call pass_var(h, G%Domain)

!     Initialize temperature and salinity (T and S).
      if ( use_temperature ) then
        call get_param(PF, mod, "TS_CONFIG", config, &
               "A string that determines how the initial tempertures \n"//&
               "and salinities are specified for a new run: \n"//&
               " \t file - read velocities from the file specified \n"//&
               " \t\t by (TS_FILE). \n"//&
               " \t fit - find the temperatures that are consistent with \n"//&
               " \t\t the layer densities and salinity S_REF. \n"//&
               " \t TS_profile - use temperature and salinity profiles \n"//&
               " \t\t (read from TS_FILE) to set layer densities. \n"//&
               " \t benchmark - use the benchmark test case T & S. \n"//&
               " \t linear - linear in logical layer space. \n"//&
               " \t DOME2D - 2D DOME initialization. \n"//&
               " \t ISOMIP - ISOMIP initialization. \n"//&
               " \t adjustment2d - TBD AJA. \n"//&
               " \t sloshing - TBD AJA. \n"//&
               " \t seamount - TBD AJA. \n"//&
               " \t rossby_front - a mixed layer front in thermal wind balance.\n"//&
               " \t SCM_ideal_hurr - used in the SCM idealized hurricane test.\n"//&
               " \t SCM_CVmix_tests - used in the SCM CVmix tests.\n"//&
               " \t USER - call a user modified routine.", &
               fail_if_missing=.true.)
!              " \t baroclinic_zone - an analytic baroclinic zone. \n"//&
        select case (trim(config))
          case ("fit"); call initialize_temp_salt_fit(tv%T, tv%S, G, GV, PF, eos, tv%P_Ref)
          case ("file"); call initialize_temp_salt_from_file(tv%T, tv%S, G, PF)
          case ("benchmark"); call benchmark_init_temperature_salinity(tv%T, tv%S, &
                                   G, GV, PF, eos, tv%P_Ref)
          case ("TS_profile") ; call initialize_temp_salt_from_profile(tv%T, tv%S, G, PF)
          case ("linear"); call initialize_temp_salt_linear(tv%T, tv%S, G, PF)
          case ("DOME2D"); call DOME2d_initialize_temperature_salinity ( tv%T, &
                                tv%S, h, G, PF, eos)
          case ("ISOMIP"); call ISOMIP_initialize_temperature_salinity ( tv%T, &
                                tv%S, h, G, GV, PF, eos)
          case ("adjustment2d"); call adjustment_initialize_temperature_salinity ( tv%T, &
                                      tv%S, h, G, PF, eos)
          case ("baroclinic_zone"); call baroclinic_zone_init_temperature_salinity( tv%T, &
                                tv%S, h, G, PF)
          case ("sloshing"); call sloshing_initialize_temperature_salinity(tv%T, &
                                  tv%S, h, G, PF, eos)
          case ("seamount"); call seamount_initialize_temperature_salinity(tv%T, &
                                  tv%S, h, G, GV, PF, eos)
          case ("rossby_front"); call Rossby_front_initialize_temperature_salinity ( tv%T, &
                                tv%S, h, G, PF, eos)
          case ("SCM_ideal_hurr"); call SCM_idealized_hurricane_TS_init ( tv%T, &
                                tv%S, h, G, GV, PF)
          case ("SCM_CVmix_tests"); call SCM_CVmix_tests_TS_init (tv%T, &
                                tv%S, h, G, GV, PF)
          case ("USER"); call user_init_temperature_salinity(tv%T, tv%S, G, PF, eos)
          case default ; call MOM_error(FATAL,  "MOM_initialize_state: "//&
                 "Unrecognized Temp & salt configuration "//trim(config))
        end select
      endif
    endif  ! not from_Z_file.

!   Initialize velocity components, u and v
    call get_param(PF, mod, "VELOCITY_CONFIG", config, &
         "A string that determines how the initial velocities \n"//&
         "are specified for a new run: \n"//&
         " \t file - read velocities from the file specified \n"//&
         " \t\t by (VELOCITY_FILE). \n"//&
         " \t zero - the fluid is initially at rest. \n"//&
         " \t uniform - the flow is uniform (determined by\n"//&
         " \t\t parameters INITIAL_U_CONST and INITIAL_V_CONST).\n"//&
         " \t rossby_front - a mixed layer front in thermal wind balance.\n"//&
         " \t soliton - Equatorial Rossby soliton.\n"//&
         " \t USER - call a user modified routine.", default="zero")
    select case (trim(config))
       case ("file"); call initialize_velocity_from_file(u, v, G, PF)
       case ("zero"); call initialize_velocity_zero(u, v, G, PF)
       case ("uniform"); call initialize_velocity_uniform(u, v, G, PF)
       case ("circular"); call initialize_velocity_circular(u, v, G, PF)
       case ("phillips"); call Phillips_initialize_velocity(u, v, G, GV, PF)
       case ("rossby_front"); call Rossby_front_initialize_velocity(u, v, h, G, GV, PF)
       case ("soliton"); call soliton_initialize_velocity(u, v, h, G)
       case ("USER"); call user_initialize_velocity(u, v, G, PF)
       case default ; call MOM_error(FATAL,  "MOM_initialize_state: "//&
            "Unrecognized velocity configuration "//trim(config))
    end select

    call pass_vector(u, v, G%Domain)
    if (debug) call uchksum(u, "MOM_initialize_state: u ", G%HI, haloshift=1)
    if (debug) call vchksum(v, "MOM_initialize_state: v ", G%HI, haloshift=1)

!   Optionally convert the thicknesses from m to kg m-2.  This is particularly
! useful in a non-Boussinesq model.
    call get_param(PF, mod, "CONVERT_THICKNESS_UNITS", convert, &
                 "If true,  convert the thickness initial conditions from \n"//&
                 "units of m to kg m-2 or vice versa, depending on whether \n"//&
                 "BOUSSINESQ is defined. This does not apply if a restart \n"//&
                 "file is read.", default=.false.)
    if (convert .and. .not. GV%Boussinesq) then
      ! Convert h from m to kg m-2 then to thickness units (H)
      call convert_thickness(h, G, GV, PF, tv)
    elseif (GV%Boussinesq) then
      ! Convert h from m to thickness units (H)
      h(:,:,:) = h(:,:,:)*GV%m_to_H
    else
      h(:,:,:) = h(:,:,:)*GV%kg_m2_to_H
    endif

!  Remove the mass that would be displaced by an ice shelf or inverse barometer.
    call get_param(PF, mod, "DEPRESS_INITIAL_SURFACE", depress_sfc, &
                 "If true,  depress the initial surface to avoid huge \n"//&
                 "tsunamis when a large surface pressure is applied.", &
                 default=.false.)
    call get_param(PF, mod, "TRIM_IC_FOR_P_SURF", trim_ic_for_p_surf, &
                 "If true, cuts way the top of the column for initial conditions\n"//&
                 "at the depth where the hydrostatic presure matches the imposed\n"//&
                 "surface pressure which is read from file.", default=.false.)
    if (depress_sfc .and. trim_ic_for_p_surf) call MOM_error(FATAL, "MOM_initialize_state: "//&
             "DEPRESS_INITIAL_SURFACE and TRIM_IC_FOR_P_SURF are exclusive and cannot both be True")
    if (depress_sfc) call depress_surface(h, G, GV, PF, tv)
    if (trim_ic_for_p_surf) call trim_for_ice(PF, G, GV, ALE_CSp, tv, h)

  else ! Previous block for new_sim=.T., this block restores state
!    This line calls a subroutine that reads the initial conditions  !
!  from a previously generated file.                                 !
    call restore_state(dirs%input_filename, dirs%restart_input_dir, Time, &
                       G, restart_CS)
    if (present(Time_in)) Time = Time_in
  endif

  if ( use_temperature ) then
    call pass_var(tv%T, G%Domain, complete=.false.)
    call pass_var(tv%S, G%Domain, complete=.false.)
  endif
  call pass_var(h, G%Domain)

  if (debug) then
    call hchksum(h*GV%H_to_m, "MOM_initialize_state: h ", G%HI, haloshift=1)
    if ( use_temperature ) call hchksum(tv%T, "MOM_initialize_state: T ", G%HI, haloshift=1)
    if ( use_temperature ) call hchksum(tv%S, "MOM_initialize_state: S ", G%HI, haloshift=1)
  endif

  call get_param(PF, mod, "SPONGE", use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. \n"//&
                 "The exact location and properties of those sponges are \n"//&
                 "specified via SPONGE_CONFIG.", default=.false.)
  if ( use_sponge ) then
! The 3 arguments here are (1) a flag indicating whether the sponge  !
! values are to be read from files, (2) the name of a file containing!
! the state toward which the model is damped, and (3) the file in    !
! which the 2-D damping rate field can be found.                     !
    call get_param(PF, mod, "SPONGE_CONFIG", config, &
                 "A string that sets how the sponges are configured: \n"//&
                 " \t file - read sponge properties from the file \n"//&
                 " \t\t specified by (SPONGE_FILE).\n"//&
                 " \t ISOMIP - apply ale sponge in the ISOMIP case \n"//&
                 " \t DOME - use a slope and channel configuration for the \n"//&
                 " \t\t DOME sill-overflow test case. \n"//&
                 " \t BFB - Sponge at the southern boundary of the domain\n"//&
                 " \t\t for buoyancy-forced basin case.\n"//&
                 " \t USER - call a user modified routine.", default="file")
    select case (trim(config))
      case ("DOME"); call DOME_initialize_sponges(G, GV, tv, PF, sponge_CSp)
      case ("DOME2D"); call DOME2d_initialize_sponges(G, GV, tv, PF, useALE, &
                                                      sponge_CSp, ALE_sponge_CSp)
      case ("ISOMIP"); call ISOMIP_initialize_sponges(G, GV, tv, PF, useALE, &
                                                     sponge_CSp, ALE_sponge_CSp)
      case ("USER"); call user_initialize_sponges(G, use_temperature, tv, &
                                               PF, sponge_CSp, h)
      case ("BFB"); call BFB_initialize_sponges_southonly(G, use_temperature, tv, &
                                               PF, sponge_CSp, h)
      case ("phillips"); call Phillips_initialize_sponges(G, use_temperature, tv, &
                                               PF, sponge_CSp, h)
      case ("file"); call initialize_sponges_file(G, GV, use_temperature, tv, &
                                               PF, sponge_CSp)
      case default ; call MOM_error(FATAL,  "MOM_initialize_state: "//&
             "Unrecognized sponge configuration "//trim(config))
    end select
  endif

  ! Reads OBC parameters not pertaining to the location of the boundaries
  call open_boundary_init(G, PF, OBC)

  ! This controls user code for setting open boundary data
  if (associated(OBC)) then
    call get_param(PF, mod, "OBC_USER_CONFIG", config, &
                 "A string that sets how the user code is invoked to set open\n"//&
                 " boundary data: \n"//&
                 "   DOME - specified inflow on northern boundary\n"//&
                 "   tidal_bay - Flather with tidal forcing on eastern boundary\n"//&
                 "   supercritical - now only needed here for the allocations\n"//&
                 "   USER - user specified", default="none")
    if (trim(config) /= "none") OBC%OBC_user_config = trim(config)
    if (open_boundary_query(OBC, apply_specified_OBC=.true.)) then
      if (trim(config) == "DOME") then
        call DOME_set_OBC_data(OBC, tv, G, GV, PF, tracer_Reg)
      elseif (lowercase(trim(config)) == "supercritical") then
        call supercritical_set_OBC_data(OBC, G, PF)
      elseif (trim(config) == "tidal_bay") then
        OBC%update_OBC = .true.
      elseif (trim(config) == "USER") then
        call user_set_OBC_data(OBC, tv, G, PF, tracer_Reg)
      elseif (.not. trim(config) == "none") then
        call MOM_error(FATAL, "The open boundary conditions specified by "//&
                "OBC_USER_CONFIG = "//trim(config)//" have not been fully implemented.")
      endif
    endif
    if (open_boundary_query(OBC, apply_open_OBC=.true.)) then
      call set_Flather_data(OBC, tv, h, G, PF, tracer_Reg)
    endif
  endif
! if (open_boundary_query(OBC, apply_nudged_OBC=.true.)) then
!   call set_3D_OBC_data(OBC, tv, h, G, PF, tracer_Reg)
! endif
  ! Still need a way to specify the boundary values
  if (debug.and.associated(OBC)) then
    call hchksum(G%mask2dT, 'MOM_initialize_state: mask2dT ', G%HI)
    call uchksum(G%mask2dCu, 'MOM_initialize_state: mask2dCu ', G%HI)
    call vchksum(G%mask2dCv, 'MOM_initialize_state: mask2dCv ', G%HI)
    call qchksum(G%mask2dBu, 'MOM_initialize_state: mask2dBu ', G%HI)
  endif

  call callTree_leave('MOM_initialize_state()')

end subroutine MOM_initialize_state
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_thickness_from_file(h, G, GV, param_file, file_has_thickness)
  type(ocean_grid_type),                  intent(in)  :: G
  type(verticalGrid_type),                intent(in)  :: GV
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: h
  type(param_file_type),                  intent(in)  :: param_file
  logical,                                intent(in)  :: file_has_thickness
! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      file_has_thickness - If true, this file contains thicknesses;
!                                 otherwise it contains interface heights.

!  This subroutine reads the layer thicknesses from file.
  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1)
  integer :: inconsistent = 0
  real :: dilate     ! The amount by which each layer is dilated to agree
                     ! with the bottom depth and free surface height, nondim.
  logical :: correct_thickness
  character(len=40)  :: mod = "initialize_thickness_from_file" ! This subroutine's name.
  character(len=200) :: filename, thickness_file, inputdir, mesg ! Strings for file/path
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mod, "THICKNESS_FILE", thickness_file, &
                 "The name of the thickness file.", fail_if_missing=.true.)

  filename = trim(inputdir)//trim(thickness_file)
  call log_param(param_file, mod, "INPUTDIR/THICKNESS_FILE", filename)

  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
         " initialize_thickness_from_file: Unable to open "//trim(filename))

  if (file_has_thickness) then
    call read_data(filename,"h",h(:,:,:),domain=G%Domain%mpp_domain)
  else
    call get_param(param_file, mod, "ADJUST_THICKNESS", correct_thickness, &
                 "If true, all mass below the bottom removed if the \n"//&
                 "topography is shallower than the thickness input file \n"//&
                 "would indicate.", default=.false.)

    call read_data(filename,"eta",eta(:,:,:),domain=G%Domain%mpp_domain)

    if (correct_thickness) then
      call adjustEtaToFitBathymetry(G, GV, eta, h)
    else
      do k=nz,1,-1 ; do j=js,je ; do i=is,ie
        if (eta(i,j,K) < (eta(i,j,K+1) + GV%Angstrom_z)) then
          eta(i,j,K) = eta(i,j,K+1) + GV%Angstrom_z
          h(i,j,k) = GV%Angstrom_z
        else
          h(i,j,k) = eta(i,j,K) - eta(i,j,K+1)
        endif
      enddo ; enddo ; enddo

      do j=js,je ; do i=is,ie
        if (abs(eta(i,j,nz+1) + G%bathyT(i,j)) > 1.0) &
          inconsistent = inconsistent + 1
      enddo ; enddo
      call sum_across_PEs(inconsistent)

      if ((inconsistent > 0) .and. (is_root_pe())) then
        write(mesg,'("Thickness initial conditions are inconsistent ",'// &
                 '"with topography in ",I8," places.")') inconsistent
        call MOM_error(WARNING, mesg)
      endif
    endif

  endif
  call callTree_leave(trim(mod)//'()')
end subroutine initialize_thickness_from_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Adjust interface heights to fit the bathymetry and diagnose layer thickness.
!! If the bottom most interface is below the topography then the bottom-most
!! layers are contracted to GV%Angstrom_z.
!! If the bottom most interface is above the topography then the entire column
!! is dilated (expanded) to fill the void.
!!   @remark{There is a (hard-wired) "tolerance" parameter such that the
!! criteria for adjustment must equal or exceed 10cm.}
!!   @param[in]     G   Grid type
!!   @param[in,out] eta Interface heights
!!   @param[out]    h   Layer thicknesses
subroutine adjustEtaToFitBathymetry(G, GV, eta, h)
  type(ocean_grid_type),                          intent(in)    :: G
  type(verticalGrid_type),                        intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)+1),    intent(inout) :: eta
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)),      intent(inout) :: h
  ! Local variables
  integer :: i, j, k, is, ie, js, je, nz, contractions, dilations
  real, parameter :: hTolerance = 0.1 !<  Tolerance to exceed adjustment criteria (m)
  real :: hTmp, eTmp, dilate
  character(len=100) :: mesg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  contractions = 0
  do j=js,je ; do i=is,ie
    if (-eta(i,j,nz+1) > G%bathyT(i,j) + hTolerance) then
      eta(i,j,nz+1) = -G%bathyT(i,j)
      contractions = contractions + 1
    endif
  enddo ; enddo
  call sum_across_PEs(contractions)
  if ((contractions > 0) .and. (is_root_pe())) then
    write(mesg,'("Thickness initial conditions were contracted ",'// &
               '"to fit topography in ",I8," places.")') contractions
    call MOM_error(WARNING, 'adjustEtaToFitBathymetry: '//mesg)
  endif

  do k=nz,1,-1 ; do j=js,je ; do i=is,ie
    ! Collapse layers to thinnest possible if the thickness less than
    ! the thinnest possible (or negative).
    if (eta(i,j,K) < (eta(i,j,K+1) + GV%Angstrom_z)) then
      eta(i,j,K) = eta(i,j,K+1) + GV%Angstrom_z
      h(i,j,k) = GV%Angstrom_z
    else
      h(i,j,k) = eta(i,j,K) - eta(i,j,K+1)
    endif
  enddo ; enddo ; enddo

  dilations = 0
  do j=js,je ; do i=is,ie
    !   The whole column is dilated to accommodate deeper topography than
    ! the bathymetry would indicate.
    ! This should be...  if ((G%mask2dt(i,j)*(eta(i,j,1)-eta(i,j,nz+1)) > 0.0) .and. &
    if (-eta(i,j,nz+1) < G%bathyT(i,j) - hTolerance) then
      dilations = dilations + 1
      if (eta(i,j,1) <= eta(i,j,nz+1)) then
        do k=1,nz ; h(i,j,k) = (eta(i,j,1)+G%bathyT(i,j)) / real(nz) ; enddo
      else
        dilate = (eta(i,j,1)+G%bathyT(i,j)) / (eta(i,j,1)-eta(i,j,nz+1))
        do k=1,nz ; h(i,j,k) = h(i,j,k) * dilate ; enddo
      endif
      do k=nz, 2, -1; eta(i,j,K) = eta(i,j,K+1) + h(i,j,k); enddo
    endif
  enddo ; enddo
  call sum_across_PEs(dilations)
  if ((dilations > 0) .and. (is_root_pe())) then
    write(mesg,'("Thickness initial conditions were dilated ",'// &
               '"to fit topography in ",I8," places.")') dilations
    call MOM_error(WARNING, 'adjustEtaToFitBathymetry: '//mesg)
  endif

end subroutine adjustEtaToFitBathymetry
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_thickness_uniform(h, G, GV, param_file)
  type(ocean_grid_type),                  intent(in)  :: G
  type(verticalGrid_type),                intent(in)  :: GV
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: h
  type(param_file_type),                  intent(in)  :: param_file

! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!  This subroutine initializes the layer thicknesses to be uniform.
  character(len=40)  :: mod = "initialize_thickness_uniform" ! This subroutine's name.
  real :: e0(SZK_(G)+1)   ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_state_initialization.F90")

  if (G%max_depth<=0.) call MOM_error(FATAL,"initialize_thickness_uniform: "// &
      "MAXIMUM_DEPTH has a non-sensical value! Was it set?")

  do k=1,nz
    e0(K) = -G%max_depth * real(k-1) / real(nz)
  enddo

  do j=js,je ; do i=is,ie
!    This sets the initial thickness (in m) of the layers.  The      !
!  thicknesses are set to insure that: 1.  each layer is at least an !
!  Angstrom thick, and 2.  the interfaces are where they should be   !
!  based on the resting depths and interface height perturbations,   !
!  as long at this doesn't interfere with 1.                         !
    eta1D(nz+1) = -1.0*G%bathyT(i,j)
    do k=nz,1,-1
      eta1D(K) = e0(K)
      if (eta1D(K) < (eta1D(K+1) + GV%Angstrom_z)) then
        eta1D(K) = eta1D(K+1) + GV%Angstrom_z
        h(i,j,k) = GV%Angstrom_z
      else
        h(i,j,k) = eta1D(K) - eta1D(K+1)
      endif
    enddo
  enddo ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine initialize_thickness_uniform
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_thickness_search
! search density space for location of layers
  call MOM_error(FATAL,"  MOM_state_initialization.F90, initialize_thickness_search: NOT IMPLEMENTED")
end subroutine initialize_thickness_search
! -----------------------------------------------------------------------------

subroutine convert_thickness(h, G, GV, param_file, tv)
  type(ocean_grid_type),                  intent(in)    :: G
  type(verticalGrid_type),                intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(inout) :: h
  type(param_file_type),                  intent(in)    :: param_file
  type(thermo_var_ptrs),                  intent(in)    :: tv
! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    p_top, p_bot
  real :: dz_geo(SZI_(G),SZJ_(G))      ! The change in geopotential height
                                       ! across a layer, in m2 s-2.
  real :: rho(SZI_(G))
  real :: I_gEarth
  logical :: Boussinesq
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: itt, max_itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  max_itt = 10
  Boussinesq = GV%Boussinesq
  I_gEarth = 1.0 / GV%g_Earth

  if (Boussinesq) then
    call MOM_error(FATAL,"Not yet converting thickness with Boussinesq approx.")
  else
    if (associated(tv%eqn_of_state)) then
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        p_bot(i,j) = 0.0 ; p_top(i,j) = 0.0
      enddo ; enddo
      do k=1,nz
        do j=js,je
          do i=is,ie ; p_top(i,j) = p_bot(i,j) ; enddo
          call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_top(:,j), rho, &
                                 is, ie-is+1, tv%eqn_of_state)
          do i=is,ie
            p_bot(i,j) = p_top(i,j) + GV%g_Earth * h(i,j,k) * rho(i)
          enddo
        enddo

        do itt=1,max_itt
          call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p_top, p_bot, &
                                   0.0, G%HI, tv%eqn_of_state, dz_geo)
          if (itt < max_itt) then ; do j=js,je
            call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_bot(:,j), rho, &
                                   is, ie-is+1, tv%eqn_of_state)
            ! Use Newton's method to correct the bottom value.
            !   The hydrostatic equation is linear to such a
            ! high degree that no bounds-checking is needed.
            do i=is,ie
              p_bot(i,j) = p_bot(i,j) + rho(i) * (GV%g_Earth*h(i,j,k) - dz_geo(i,j))
            enddo
          enddo ; endif
        enddo

        do j=js,je ; do i=is,ie
          h(i,j,k) = (p_bot(i,j) - p_top(i,j)) * GV%kg_m2_to_H * I_gEarth
        enddo ; enddo
      enddo
    else
      do k=1,nz ; do j=js,je ; do i=is,ie
        h(i,j,k) = h(i,j,k) * GV%Rlay(k) * GV%kg_m2_to_H
      enddo ; enddo ; enddo
    endif
  endif

end subroutine convert_thickness

subroutine depress_surface(h, G, GV, param_file, tv)
  type(ocean_grid_type),                  intent(in)    :: G
  type(verticalGrid_type),                intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(inout) :: h
  type(param_file_type),                  intent(in)    :: param_file
  type(thermo_var_ptrs),                  intent(in)    :: tv
! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

  real, dimension(SZI_(G),SZJ_(G)) :: &
    eta_sfc  ! The free surface height that the model should use, in m.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: &
    eta  ! The free surface height that the model should use, in m.
  real :: dilate  ! A ratio by which layers are dilated, nondim.
  real :: scale_factor ! A scaling factor for the eta_sfc values that are read
                       ! in, which can be used to change units, for example.
  character(len=40)  :: mod = "depress_surface" ! This subroutine's name.
  character(len=200) :: inputdir, eta_srf_file ! Strings for file/path
  character(len=200) :: filename, eta_srf_var  ! Strings for file/path
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! Read the surface height (or pressure) from a file.

  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mod, "SURFACE_HEIGHT_IC_FILE", eta_srf_file,&
                 "The initial condition file for the surface height.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "SURFACE_HEIGHT_IC_VAR", eta_srf_var, &
                 "The initial condition variable for the surface height.",&
                 default="SSH")
  filename = trim(inputdir)//trim(eta_srf_file)
  call log_param(param_file,  mod, "INPUTDIR/SURFACE_HEIGHT_IC_FILE", filename)

  call read_data(filename,eta_srf_var,eta_sfc,domain=G%Domain%mpp_domain)
  call get_param(param_file, mod, "SURFACE_HEIGHT_IC_SCALE", scale_factor, &
                 "A scaling factor to convert SURFACE_HEIGHT_IC_VAR into \n"//&
                 "units of m", units="variable", default=1.0)

  if (scale_factor /= 1.0) then ; do j=js,je ; do i=is,ie
    eta_sfc(i,j) = eta_sfc(i,j) * scale_factor
  enddo ; enddo ; endif

  ! Convert thicknesses to interface heights.
  call find_eta(h, tv, GV%g_Earth, G, GV, eta)

  do j=js,je ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
!    if (eta_sfc(i,j) < eta(i,j,nz+1)) then
      ! Issue a warning?
!    endif
    if (eta_sfc(i,j) > eta(i,j,1)) then
      ! Dilate the water column to agree, but only up to 10-fold.
      if (eta_sfc(i,j) - eta(i,j,nz+1) > 10.0*(eta(i,j,1) - eta(i,j,nz+1))) then
        dilate = 10.0
        call MOM_error(WARNING, "Free surface height dilation attempted "//&
               "to exceed 10-fold.", all_print=.true.)
      else
        dilate = (eta_sfc(i,j) - eta(i,j,nz+1)) / (eta(i,j,1) - eta(i,j,nz+1))
      endif
      do k=1,nz ; h(i,j,k) = h(i,j,k) * dilate ; enddo
    elseif (eta(i,j,1) > eta_sfc(i,j)) then
      ! Remove any mass that is above the target free surface.
      do k=1,nz
        if (eta(i,j,K) <= eta_sfc(i,j)) exit
        if (eta(i,j,K+1) >= eta_sfc(i,j)) then
          h(i,j,k) = GV%Angstrom
        else
          h(i,j,k) = max(GV%Angstrom, h(i,j,k) * &
              (eta_sfc(i,j) - eta(i,j,K+1)) / (eta(i,j,K) - eta(i,j,K+1)) )
        endif
      enddo
    endif
  endif ; enddo ; enddo

end subroutine depress_surface

!> Adjust the layer thicknesses by cutting away the top of each model column at the depth
!! where the hydrostatic pressure matches an imposed surface pressure read from file.
subroutine trim_for_ice(PF, G, GV, ALE_CSp, tv, h)
  type(param_file_type),                    intent(in)    :: PF !< Parameter file structure
  type(ocean_grid_type),                    intent(in)    :: G !< Ocean grid structure
  type(verticalGrid_type),                  intent(in)    :: GV !< Vertical grid structure
  type(ALE_CS),                             pointer       :: ALE_CSp !< ALE control structure
  type(thermo_var_ptrs),                    intent(inout) :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: h !< Layer thickness (H units, m or Pa)
  ! Local variables
  character(len=200) :: mod = "trim_for_ice"
  real, dimension(SZI_(G),SZJ_(G)) :: p_surf ! Imposed pressure on ocean at surface (Pa)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: S_t, S_b, T_t, T_b ! Top and bottom edge values for reconstructions
                                                                 ! of salinity and temperature within each layer.
  character(len=200) :: inputdir, filename, p_surf_file, p_surf_var ! Strings for file/path
  real :: scale_factor, min_thickness
  integer :: i, j, k
  logical :: use_remapping
  type(remapping_CS), pointer :: remap_CS => NULL()

  call get_param(PF, mod, "SURFACE_PRESSURE_FILE", p_surf_file, &
                 "The initial condition file for the surface height.", &
                 fail_if_missing=.true.)
  call get_param(PF, mod, "SURFACE_PRESSURE_VAR", p_surf_var, &
                 "The initial condition variable for the surface height.", &
                 units="kg m-2", default="")
  call get_param(PF, mod, "INPUTDIR", inputdir, default=".", do_not_log=.true.)
  filename = trim(slasher(inputdir))//trim(p_surf_file)
  call log_param(PF,  mod, "!INPUTDIR/SURFACE_HEIGHT_IC_FILE", filename)

  call read_data(filename, p_surf_var, p_surf, domain=G%Domain%mpp_domain)
  call get_param(PF, mod, "SURFACE_PRESSURE_SCALE", scale_factor, &
                 "A scaling factor to convert SURFACE_PRESSURE_VAR from\n"//&
                 "file SURFACE_PRESSURE_FILE into a surface pressure.", &
                 units="file dependent", default=1.)
  if (scale_factor /= 1.) p_surf(:,:) = scale_factor * p_surf(:,:)
  call get_param(PF, mod, "MIN_THICKNESS", min_thickness, 'Minimum layer thickness', &
                 units='m', default=1.e-3)
  call get_param(PF, mod, "TRIMMING_USES_REMAPPING", use_remapping, &
                 'When trimming the column, also remap T and S.', &
                 default=.false.)
  if (use_remapping) then
    allocate(remap_CS)
    call initialize_remapping(remap_CS, 'PLM', boundary_extrapolation=.true.)
  endif

  ! Find edge values of T and S used in reconstructions
  if ( associated(ALE_CSp) ) then ! This should only be associated if we are in ALE mode
!   if ( PRScheme == PRESSURE_RECONSTRUCTION_PLM ) then
      call pressure_gradient_plm(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h)
!   elseif ( PRScheme == PRESSURE_RECONSTRUCTION_PPM ) then
!     call pressure_gradient_ppm(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h)
!   endif
  else
!    call MOM_error(FATAL, "trim_for_ice: Does not work without ALE mode")
    do k=1,G%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      T_t(i,j,k) = tv%T(i,j,k) ; T_b(i,j,k) = tv%T(i,j,k)
      S_t(i,j,k) = tv%S(i,j,k) ; S_b(i,j,k) = tv%S(i,j,k)
    enddo ; enddo ; enddo
  endif

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    call cut_off_column_top(GV%ke, tv, GV%Rho0, GV%g_Earth, G%bathyT(i,j), min_thickness, &
               tv%T(i,j,:), T_t(i,j,:), T_b(i,j,:), tv%S(i,j,:), S_t(i,j,:), S_b(i,j,:), &
               p_surf(i,j), h(i,j,:), remap_CS)
  enddo ; enddo

end subroutine trim_for_ice

!> Adjust the layer thicknesses by cutting away the top at the depth where the hydrostatic
!! pressure matches p_surf
subroutine cut_off_column_top(nk, tv, Rho0, G_earth, depth, min_thickness, &
                              T, T_t, T_b, S, S_t, S_b, p_surf, h, remap_CS)
  integer,               intent(in)    :: nk !< Number of layers
  type(thermo_var_ptrs), intent(in)    :: tv !< Thermodynamics structure
  real,                  intent(in)    :: Rho0 !< Reference density (kg/m3)
  real,                  intent(in)    :: G_earth !< Gravitational acceleration (m/s2)
  real,                  intent(in)    :: depth !< Depth of ocean column (m)
  real,                  intent(in)    :: min_thickness !< Smallest thickness allowed (m)
  real, dimension(nk),   intent(inout) :: T !< Layer mean temperature
  real, dimension(nk),   intent(in)    :: T_t !< Temperature at top of layer
  real, dimension(nk),   intent(in)    :: T_b !< Temperature at bottom of layer
  real, dimension(nk),   intent(inout) :: S !< Layer mean salinity
  real, dimension(nk),   intent(in)    :: S_t !< Salinity at top of layer
  real, dimension(nk),   intent(in)    :: S_b !< Salinity at bottom of layer
  real,                  intent(in)    :: p_surf !< Imposed pressure on ocean at surface (Pa)
  real, dimension(nk),   intent(inout) :: h !< Layer thickness (H units, m or Pa)
  type(remapping_CS),    pointer       :: remap_CS ! Remapping structure for remapping T and S, if associated
  ! Local variables
  real, dimension(nk+1) :: e ! Top and bottom edge values for reconstructions
  real, dimension(nk) :: h0, S0, T0, h1, S1, T1
  real :: P_t, P_b, z_out, e_top
  integer :: k

  ! Calculate original interface positions
  e(nk+1) = -depth
  do k=nk,1,-1
    e(K) = e(K+1) + h(k)
    h0(k) = h(nk+1-k) ! Keep a copy to use in remapping
  enddo

  P_t = 0.
  e_top = e(1)
  do k=1,nk
    call find_depth_of_pressure_in_cell(T_t(k), T_b(k), S_t(k), S_b(k), e(K), e(K+1), &
                                        P_t, p_surf, Rho0, G_earth, tv%eqn_of_state, P_b, z_out)
    if (z_out>=e(K)) then
      ! Imposed pressure was less that pressure at top of cell
      exit
    elseif (z_out<=e(K+1)) then
      ! Imposed pressure was greater than pressure at bottom of cell
      e_top = e(K+1)
    else
      ! Imposed pressure was fell between pressures at top and bottom of cell
      e_top = z_out
      exit
    endif
    P_t = P_b
  enddo
  if (e_top<e(1)) then
    ! Clip layers from the top down, if at all
    do K=1,nk
      if (e(K)>e_top) then
        ! Original e(K) is too high
        e(K) = e_top
        e_top = e_top - min_thickness ! Next interface must be at least this deep
      endif
      ! This layer needs trimming
      h(k) = max( min_thickness, e(K) - e(K+1) )
      if (e(K)<e_top) exit ! No need to go further
    enddo
  endif

  ! Now we need to remap but remapping assumes the surface is at the
  ! same place in the two columns so we turn the column upside down.
  if (associated(remap_CS)) then
    do k=1,nk
      S0(k) = S(nk+1-k)
      T0(k) = T(nk+1-k)
      h1(k) = h(nk+1-k)
    enddo
    call remapping_core_h(nk, h0, T0, nk, h1, T1, remap_CS )
    call remapping_core_h(nk, h0, S0, nk, h1, S1, remap_CS )
    do k=1,nk
      S(k) = S1(nk+1-k)
      T(k) = T1(nk+1-k)
    enddo
  endif

end subroutine cut_off_column_top

! -----------------------------------------------------------------------------
subroutine initialize_velocity_from_file(u, v, G, param_file)
  type(ocean_grid_type),                   intent(in)  :: G
  real, dimension(SZIB_(G),SZJ_(G), SZK_(G)), intent(out) :: u
  real, dimension(SZI_(G),SZJB_(G), SZK_(G)), intent(out) :: v
  type(param_file_type),                   intent(in)  :: param_file
! Arguments: u - The zonal velocity that is being initialized.
!  (out)     v - The meridional velocity that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file -  parameter file type

!   This subroutine reads the initial velocity components from file
  character(len=40)  :: mod = "initialize_velocity_from_file" ! This subroutine's name.
  character(len=200) :: filename,velocity_file,inputdir ! Strings for file/path

  call callTree_enter(trim(mod)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mod, "VELOCITY_FILE", velocity_file, &
                 "The name of the velocity initial condition file.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  filename = trim(inputdir)//trim(velocity_file)
  call log_param(param_file, mod, "INPUTDIR/VELOCITY_FILE", filename)

  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
         " initialize_velocity_from_file: Unable to open "//trim(filename))

  !  Read the velocities from a netcdf file.
  call read_data(filename,"u",u(:,:,:),domain=G%Domain%mpp_domain,position=EAST_FACE)
  call read_data(filename,"v",v(:,:,:),domain=G%Domain%mpp_domain,position=NORTH_FACE)

  call callTree_leave(trim(mod)//'()')
end subroutine initialize_velocity_from_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_velocity_zero(u, v, G, param_file)
  type(ocean_grid_type),                   intent(in)  :: G
  real, dimension(SZIB_(G),SZJ_(G), SZK_(G)), intent(out) :: u
  real, dimension(SZI_(G),SZJB_(G), SZK_(G)), intent(out) :: v
  type(param_file_type),                   intent(in)  :: param_file
! Arguments: u - The zonal velocity that is being initialized.
!  (out)     v - The meridional velocity that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file -  parameter file type

!   This subroutine sets the initial velocity components to zero
  character(len=200) :: mod = "initialize_velocity_zero" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  call callTree_enter(trim(mod)//"(), MOM_state_initialization.F90")

  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u(I,j,k) = 0.0
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v(i,J,k) = 0.0
  enddo ; enddo ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine initialize_velocity_zero
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_velocity_uniform(u, v, G, param_file)
  type(ocean_grid_type),                   intent(in)  :: G
  real, dimension(SZIB_(G),SZJ_(G), SZK_(G)), intent(out) :: u
  real, dimension(SZI_(G),SZJB_(G), SZK_(G)), intent(out) :: v
  type(param_file_type),                   intent(in)  :: param_file
! Arguments: u - The zonal velocity that is being initialized.
!  (out)     v - The meridional velocity that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file -  parameter file type

!   This subroutine sets the initial velocity components to uniform
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  real    :: initial_u_const, initial_v_const
  character(len=200) :: mod = "initialize_velocity_uniform" ! This subroutine's name.
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  call get_param(param_file, mod, "INITIAL_U_CONST", initial_u_const, &
                 "A initial uniform value for the zonal flow.", &
                 units="m s-1", fail_if_missing=.true.)
  call get_param(param_file, mod, "INITIAL_V_CONST", initial_v_const, &
                 "A initial uniform value for the meridional flow.", &
                 units="m s-1", fail_if_missing=.true.)

  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u(I,j,k) = initial_u_const
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v(i,J,k) = initial_v_const
  enddo ; enddo ; enddo

end subroutine initialize_velocity_uniform
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_velocity_circular(u, v, G, param_file)
  type(ocean_grid_type),                   intent(in)  :: G
  real, dimension(SZIB_(G),SZJ_(G), SZK_(G)), intent(out) :: u
  real, dimension(SZI_(G),SZJB_(G), SZK_(G)), intent(out) :: v
  type(param_file_type),                   intent(in)  :: param_file
! Arguments: u - The zonal velocity that is being initialized.
!  (out)     v - The meridional velocity that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file -  parameter file type

!   This subroutine sets the initial velocity components to be circular with
! no flow at edges of domain and center.
  character(len=200) :: mod = "initialize_velocity_circular"
  real :: circular_max_u
  real :: dpi, psi1, psi2
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  call get_param(param_file, mod, "CIRCULAR_MAX_U", circular_max_u, &
                 "The amplitude of zonal flow from which to scale the\n"// &
                 "circular stream function (m/s).", &
                 units="m s-1", default=0.)
  dpi=acos(0.0)*2.0 ! pi

  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    psi1 = my_psi(I,j)
    psi2 = my_psi(I,j-1)
    u(I,j,k) = (psi1-psi2)/G%dy_Cu(I,j)! *(circular_max_u*G%len_lon/(2.0*dpi))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    psi1 = my_psi(i,J)
    psi2 = my_psi(i-1,J)
    v(i,J,k) = (psi2-psi1)/G%dx_Cv(i,J)! *(circular_max_u*G%len_lon/(2.0*dpi))
  enddo ; enddo ; enddo

  contains

  real function my_psi(ig,jg) ! in-line function
    integer :: ig, jg
    real :: x, y, r
    x = 2.0*(G%geoLonBu(ig,jg)-G%west_lon)/G%len_lon-1.0  ! -1<x<1
    y = 2.0*(G%geoLatBu(ig,jg)-G%south_lat)/G%len_lat-1.0 ! -1<y<1
    r = sqrt( x**2 + y**2 ) ! Circulat stream fn nis fn of radius only
    r = min(1.0,r) ! Flatten stream function in corners of box
    my_psi = 0.5*(1.0 - cos(dpi*r))
    my_psi = my_psi * (circular_max_u*G%len_lon*1e3/dpi) ! len_lon is in km
  end function my_psi

end subroutine initialize_velocity_circular
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_temp_salt_from_file(T, S, G, param_file)
  type(ocean_grid_type),                  intent(in)  :: G
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T, S
  type(param_file_type),                  intent(in)  :: param_file
!  This function puts the initial layer temperatures and salinities  !
! into T(:,:,:) and S(:,:,:).                                        !

! Arguments: T - The potential temperature that is being initialized.
!  (out)     S - The salinity that is being initialized.
!  (in)      from_file - .true. if the variables that are set here are to
!                        be read from a file; .false. to be set internally.
!  (in)      filename - The name of the file to read.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
  character(len=200) :: filename, ts_file, salt_file, inputdir ! Strings for file/path
  character(len=40)  :: mod = "initialize_temp_salt_from_file"
  character(len=64)  :: temp_var, salt_var

  call callTree_enter(trim(mod)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mod, "TS_FILE", ts_file, &
                 "The initial condition file for temperature.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  filename = trim(inputdir)//trim(ts_file)
  call log_param(param_file, mod, "INPUTDIR/TS_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
     " initialize_temp_salt_from_file: Unable to open "//trim(filename))

  call get_param(param_file, mod, "TEMP_IC_VAR", temp_var, &
                 "The initial condition variable for potential temperature.", &
                 default="PTEMP")
  call get_param(param_file, mod, "SALT_IC_VAR", salt_var, &
                 "The initial condition variable for salinity.", default="SALT")

! Read the temperatures and salinities from a netcdf file.           !
  call read_data(filename, temp_var, T(:,:,:), domain=G%Domain%mpp_domain)

  call get_param(param_file, mod, "SALT_FILE", salt_file, &
                 "The initial condition file for salinity.", default=trim(ts_file))
  filename = trim(inputdir)//trim(ts_file)
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
     " initialize_temp_salt_from_file: Unable to open "//trim(filename))

  call read_data(filename, salt_var, S(:,:,:), domain=G%Domain%mpp_domain)

  call callTree_leave(trim(mod)//'()')
end subroutine initialize_temp_salt_from_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_temp_salt_from_profile(T, S, G, param_file)
  type(ocean_grid_type),                  intent(in)  :: G
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T, S
  type(param_file_type),                  intent(in)  :: param_file
!  This function puts the initial layer temperatures and salinities  !
! into T(:,:,:) and S(:,:,:).                                        !

! Arguments: T - The potential temperature that is being initialized.
!  (out)     S - The salinity that is being initialized.
!  (in)      from_file - .true. if the variables that are set here are to
!                        be read from a file; .false. to be set internally.
!  (in)      filename - The name of the file to read.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
  real, dimension(SZK_(G)) :: T0, S0
  integer :: i, j, k
  character(len=200) :: filename, ts_file, inputdir ! Strings for file/path
  character(len=40)  :: mod = "initialize_temp_salt_from_profile"

  call callTree_enter(trim(mod)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mod, "TS_FILE", ts_file, &
                 "The file with the reference profiles for temperature \n"//&
                 "and salinity.", fail_if_missing=.true.)
  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  filename = trim(inputdir)//trim(ts_file)
  call log_param(param_file, mod, "INPUTDIR/TS_FILE", filename)
  if (.not.file_exists(filename)) call MOM_error(FATAL, &
     " initialize_temp_salt_from_profile: Unable to open "//trim(filename))

! Read the temperatures and salinities from a netcdf file.           !
  call read_data(filename,"PTEMP",T0(:),domain=G%Domain%mpp_domain)
  call read_data(filename,"SALT", S0(:),domain=G%Domain%mpp_domain)

  do k=1,G%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    T(i,j,k) = T0(k) ; S(i,j,k) = S0(k)
  enddo ; enddo ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine initialize_temp_salt_from_profile
! -----------------------------------------------------------------------------


! -----------------------------------------------------------------------------
subroutine initialize_temp_salt_fit(T, S, G, GV, param_file, eqn_of_state, P_Ref)
  type(ocean_grid_type),                  intent(in)  :: G
  type(verticalGrid_type),                intent(in)  :: GV
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T, S
  type(param_file_type),                  intent(in)  :: param_file
  type(EOS_type),                         pointer     :: eqn_of_state
  real,                                   intent(in)  :: P_Ref
!  This function puts the initial layer temperatures and salinities  !
! into T(:,:,:) and S(:,:,:).                                        !

! Arguments: T - The potential temperature that is being initialized.
!  (out)     S - The salinity that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      eqn_of_state - integer that selects the equatio of state
!  (in)      P_Ref - The coordinate-density reference pressure in Pa.
  real :: T0(SZK_(G)), S0(SZK_(G))
  real :: T_Ref         ! Reference Temperature
  real :: S_Ref         ! Reference Salinity
  real :: pres(SZK_(G))      ! An array of the reference pressure in Pa.
  real :: drho_dT(SZK_(G))   ! Derivative of density with temperature in kg m-3 K-1.                              !
  real :: drho_dS(SZK_(G))   ! Derivative of density with salinity in kg m-3 PSU-1.                             !
  real :: rho_guess(SZK_(G)) ! Potential density at T0 & S0 in kg m-3.
  logical :: fit_salin       ! If true, accept the prescribed temperature and fit the salinity.
  character(len=40)  :: mod = "initialize_temp_salt_fit" ! This subroutine's name.
  integer :: i, j, k, itt, nz
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mod, "T_REF", T_Ref, &
                 "A reference temperature used in initialization.", &
                 units="degC", fail_if_missing=.true.)
  call get_param(param_file, mod, "S_REF", S_Ref, &
                 "A reference salinity used in initialization.", units="PSU", &
                 default=35.0)
  call get_param(param_file, mod, "FIT_SALINITY", fit_salin, &
                 "If true, accept the prescribed temperature and fit the \n"//&
                 "salinity; otherwise take salinity and fit temperature.", &
                 default=.false.)
  do k=1,nz
    pres(k) = P_Ref ; S0(k) = S_Ref
    T0(k) = T_Ref
  enddo

  call calculate_density(T0(1),S0(1),pres(1),rho_guess(1),eqn_of_state)
  call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,1,eqn_of_state)

  if (fit_salin) then
! A first guess of the layers' temperatures.
    do k=nz,1,-1
      S0(k) = max(0.0, S0(1) + (GV%Rlay(k) - rho_guess(1)) / drho_dS(1))
    enddo
! Refine the guesses for each layer.
    do itt=1,6
      call calculate_density(T0,S0,pres,rho_guess,1,nz,eqn_of_state)
      call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,nz,eqn_of_state)
      do k=1,nz
        S0(k) = max(0.0, S0(k) + (GV%Rlay(k) - rho_guess(k)) / drho_dS(k))
      enddo
    enddo
  else
! A first guess of the layers' temperatures.
    do k=nz,1,-1
      T0(k) = T0(1) + (GV%Rlay(k) - rho_guess(1)) / drho_dT(1)
    enddo
    do itt=1,6
      call calculate_density(T0,S0,pres,rho_guess,1,nz,eqn_of_state)
      call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,nz,eqn_of_state)
      do k=1,nz
        T0(k) = T0(k) + (GV%Rlay(k) - rho_guess(k)) / drho_dT(k)
      enddo
    enddo
  endif

  do k=1,nz ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
    T(i,j,k) = T0(k) ; S(i,j,k) = S0(k)
  enddo ; enddo ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine initialize_temp_salt_fit
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_temp_salt_linear(T, S, G, param_file)
  type(ocean_grid_type),                  intent(in)  :: G
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T, S
  type(param_file_type),                  intent(in)  :: param_file
  ! This subroutine initializes linear profiles for T and S according to
  ! reference surface layer salinity and temperature and a specified range.
  ! Note that the linear distribution is set up with respect to the layer
  ! number, not the physical position).
  integer :: k;
  real  :: delta_S, delta_T
  real  :: S_top, T_top ! Reference salinity and temerature within surface layer
  real  :: S_range, T_range ! Range of salinities and temperatures over the vertical
  real  :: delta
  character(len=40)  :: mod = "initialize_temp_salt_linear" ! This subroutine's name.

  call callTree_enter(trim(mod)//"(), MOM_state_initialization.F90")
  call get_param(param_file, mod, "T_TOP", T_top, &
                 "Initial temperature of the top surface.", &
                 units="degC", fail_if_missing=.true.)
  call get_param(param_file, mod, "T_RANGE", T_range, &
                 "Initial temperature difference (top-bottom).", &
                 units="degC", fail_if_missing=.true.)
  call get_param(param_file, mod, "S_TOP", S_top, &
                 "Initial salinity of the top surface.", &
                 units="PSU", fail_if_missing=.true.)
  call get_param(param_file, mod, "S_RANGE", S_range, &
                 "Initial salinity difference (top-bottom).", &
                 units="PSU", fail_if_missing=.true.)

! ! Prescribe salinity
! delta_S = S_range / ( G%ke - 1.0 );
! S(:,:,1) = S_top;
! do k = 2,G%ke
!   S(:,:,k) = S(:,:,k-1) + delta_S;
! end do
  do k = 1,G%ke
    S(:,:,k) = S_top - S_range*((real(k)-0.5)/real(G%ke))
    T(:,:,k) = T_top - T_range*((real(k)-0.5)/real(G%ke))
  end do

! ! Prescribe temperature
! delta_T = T_range / ( G%ke - 1.0 );
! T(:,:,1) = T_top;
! do k = 2,G%ke
!   T(:,:,k) = T(:,:,k-1) + delta_T;
! end do
! delta = 1;
! T(:,:,G%ke/2 - (delta-1):G%ke/2 + delta) = 1.0;

  call callTree_leave(trim(mod)//'()')
end subroutine initialize_temp_salt_linear
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_sponges_file(G, GV, use_temperature, tv, param_file, CSp)
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  logical,               intent(in) :: use_temperature
  type(thermo_var_ptrs), intent(in) :: tv
  type(param_file_type), intent(in) :: param_file
  type(sponge_CS),       pointer    :: CSp
!   This subroutine sets the inverse restoration time (Idamp), and   !
! the values towards which the interface heights and an arbitrary    !
! number of tracers should be restored within each sponge. The       !
! interface height is always subject to damping, and must always be  !
! the first registered field.                                        !

! Arguments: from_file - .true. if the variables that are used here are to
!                        be read from a file; .false. to be set internally.
!  (in)      filename - The name of the file to read for all fields
!                       except the inverse damping rate.
!  (in)      damp_file - The name of the file from which to read the
!                        inverse damping rate.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      use_temperature - If true, T & S are state variables.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CSp - A pointer that is set to point to the control structure
!                  for this module

  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! The target interface heights, in m.
  real, dimension (SZI_(G),SZJ_(G),SZK_(G)) :: &
    tmp, tmp2 ! A temporary array for tracers.
  real, dimension (SZI_(G),SZJ_(G)) :: &
    tmp_2d ! A temporary array for tracers.

  real :: Idamp(SZI_(G),SZJ_(G))    ! The inverse damping rate, in s-1.
  real :: pres(SZI_(G))     ! An array of the reference pressure, in Pa.

  integer :: i, j, k, is, ie, js, je, nz
  character(len=40) :: potemp_var, salin_var, Idamp_var, eta_var
  character(len=40) :: mod = "initialize_sponges_file"
  character(len=200) :: damping_file, state_file  ! Strings for filenames
  character(len=200) :: filename, inputdir ! Strings for file/path and path.
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  pres(:) = 0.0 ; eta(:,:,:) = 0.0 ; tmp(:,:,:) = 0.0 ; Idamp(:,:) = 0.0

  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mod, "SPONGE_DAMPING_FILE", damping_file, &
                 "The name of the file with the sponge damping rates.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "SPONGE_STATE_FILE", state_file, &
                 "The name of the file with the state to damp toward.", &
                 default=damping_file)

  call get_param(param_file, mod, "SPONGE_PTEMP_VAR", potemp_var, &
                 "The name of the potential temperature variable in \n"//&
                 "SPONGE_STATE_FILE.", default="PTEMP")
  call get_param(param_file, mod, "SPONGE_SALT_VAR", salin_var, &
                 "The name of the salinity variable in \n"//&
                 "SPONGE_STATE_FILE.", default="SALT")
  call get_param(param_file, mod, "SPONGE_ETA_VAR", eta_var, &
                 "The name of the interface height variable in \n"//&
                 "SPONGE_STATE_FILE.", default="ETA")
  call get_param(param_file, mod, "SPONGE_IDAMP_VAR", Idamp_var, &
                 "The name of the inverse damping rate variable in \n"//&
                 "SPONGE_DAMPING_FILE.", default="IDAMP")

  filename = trim(inputdir)//trim(damping_file)
  call log_param(param_file, mod, "INPUTDIR/SPONGE_DAMPING_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) &
    call MOM_error(FATAL, " initialize_sponges: Unable to open "//trim(filename))


  call read_data(filename,"Idamp",Idamp(:,:), domain=G%Domain%mpp_domain)

! Now register all of the fields which are damped in the sponge.     !
! By default, momentum is advected vertically within the sponge, but !
! momentum is typically not damped within the sponge.                !

  filename = trim(inputdir)//trim(state_file)
  call log_param(param_file, mod, "INPUTDIR/SPONGE_STATE_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) &
    call MOM_error(FATAL, " initialize_sponges: Unable to open "//trim(filename))


!  The first call to set_up_sponge_field is for the interface height.!
  call read_data(filename, eta_var, eta(:,:,:), domain=G%Domain%mpp_domain)

  do j=js,je ; do i=is,ie
    eta(i,j,nz+1) = -G%bathyT(i,j)
  enddo ; enddo
  do k=nz,1,-1 ; do j=js,je ; do i=is,ie
    if (eta(i,j,K) < (eta(i,j,K+1) + GV%Angstrom_z)) &
      eta(i,j,K) = eta(i,j,K+1) + GV%Angstrom_z
  enddo ; enddo ; enddo
! Set the inverse damping rates so that the model will know where to !
! apply the sponges, along with the interface heights.               !
  call initialize_sponge(Idamp, eta, G, param_file, CSp)

!   Now register all of the tracer fields which are damped in the    !
! sponge. By default, momentum is advected vertically within the     !
! sponge, but momentum is typically not damped within the sponge.    !

  if ( GV%nkml>0 ) then
!   This call to set_up_sponge_ML_density registers the target values of the
! mixed layer density, which is used in determining which layers can be
! inflated without causing static instabilities.
    do i=is-1,ie ; pres(i) = tv%P_Ref ; enddo

    call read_data(filename, potemp_var, tmp(:,:,:), domain=G%Domain%mpp_domain)
    call read_data(filename, salin_var, tmp2(:,:,:), domain=G%Domain%mpp_domain)

    do j=js,je
      call calculate_density(tmp(:,j,1), tmp2(:,j,1), pres, tmp_2d(:,j), &
                             is, ie-is+1, tv%eqn_of_state)
    enddo

    call set_up_sponge_ML_density(tmp_2d, G, CSp)
  endif

!  The remaining calls to set_up_sponge_field can be in any order.   !
  if ( use_temperature ) then
    call read_data(filename, potemp_var, tmp(:,:,:), domain=G%Domain%mpp_domain)
    call set_up_sponge_field(tmp, tv%T, G, nz, CSp)
    call read_data(filename, salin_var, tmp(:,:,:), domain=G%Domain%mpp_domain)
    call set_up_sponge_field(tmp, tv%S, G, nz, CSp)
  endif


end subroutine initialize_sponges_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_velocity_depth_max(G)
  type(ocean_grid_type), intent(inout) :: G
  ! This subroutine sets the 4 bottom depths at velocity points to be the
  ! maximum of the adjacent depths.
  integer :: i, j

  do I=G%isd,G%ied-1 ; do j=G%jsd,G%jed
    G%Dblock_u(I,j) = G%mask2dCu(I,j) * max(G%bathyT(i,j), G%bathyT(i+1,j))
    G%Dopen_u(I,j) = G%Dblock_u(I,j)
  enddo ; enddo
  do i=G%isd,G%ied ; do J=G%jsd,G%jed-1
    G%Dblock_v(I,J) = G%mask2dCv(i,J) * max(G%bathyT(i,j), G%bathyT(i,j+1))
    G%Dopen_v(I,J) = G%Dblock_v(I,J)
  enddo ; enddo
end subroutine set_velocity_depth_max
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine compute_global_grid_integrals(G)
  type(ocean_grid_type), intent(inout) :: G
  ! Subroutine to pre-compute global integrals of grid quantities for
  ! later use in reporting diagnostics
  integer :: i,j

  G%areaT_global = 0.0 ; G%IareaT_global = 0.0
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    G%areaT_global = G%areaT_global + ( G%areaT(i,j) * G%mask2dT(i,j) )
  enddo ; enddo
  call sum_across_PEs( G%areaT_global )
  G%IareaT_global = 1. / G%areaT_global
end subroutine compute_global_grid_integrals

! -----------------------------------------------------------------------------
subroutine set_velocity_depth_min(G)
  type(ocean_grid_type), intent(inout) :: G
  ! This subroutine sets the 4 bottom depths at velocity points to be the
  ! minimum of the adjacent depths.
  integer :: i, j

  do I=G%isd,G%ied-1 ; do j=G%jsd,G%jed
    G%Dblock_u(I,j) = G%mask2dCu(I,j) * min(G%bathyT(i,j), G%bathyT(i+1,j))
    G%Dopen_u(I,j) = G%Dblock_u(I,j)
  enddo ; enddo
  do i=G%isd,G%ied ; do J=G%jsd,G%jed-1
    G%Dblock_v(I,J) = G%mask2dCv(i,J) * min(G%bathyT(i,j), G%bathyT(i,j+1))
    G%Dopen_v(I,J) = G%Dblock_v(I,J)
  enddo ; enddo
end subroutine set_velocity_depth_min
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine MOM_temp_salt_initialize_from_Z(h, tv, G, GV, PF, dirs)

! Determines the isopycnal interfaces and layer potential
! temperatures and salinities directly from a z-space file on a latitude-
! longitude grid.
!
! This subroutine was written by M. Harrison, with input from R. Hallberg.
! and A. Adcroft.
!
! Arguments:
!  (out)     h  - Layer thickness, in m.
!  (out)     tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (inout)   G       - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      PF      - A structure indicating the open file to parse for
!                      model parameter values.
!  (in)      dirs    - A structure containing several relevant directory paths.

  type(ocean_grid_type),                 intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out)   :: h
  type(thermo_var_ptrs),                 intent(inout) :: tv
  type(verticalGrid_type),               intent(in)    :: GV
  type(param_file_type),                 intent(in)    :: PF
  type(directories),                     intent(in)    :: dirs

  character(len=200) :: filename   ! The name of an input file containing temperature
                                   ! and salinity in z-space; also used for  ice shelf area.
  character(len=200) :: inputdir ! The directory where NetCDF input files are.
  character(len=200) :: mesg, area_varname, ice_shelf_file

  type(EOS_type), pointer :: eos => NULL()

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_initialize_layers_from_Z" ! This module's name.

  integer :: is, ie, js, je, nz ! compute domain indices
  integer :: isc,iec,jsc,jec    ! global compute domain indices
  integer :: isg, ieg, jsg, jeg ! global extent
  integer :: isd, ied, jsd, jed ! data domain indices

  integer :: i, j, k, ks, np, ni, nj
  integer :: idbg, jdbg
  integer :: nkml, nkbl         ! number of mixed and buffer layers

  integer :: kd, inconsistent
  real    :: PI_180             ! for conversion from degrees to radians

  real, dimension(:,:), pointer :: shelf_area
  real    :: min_depth
  real    :: dilate
  real    :: missing_value_temp, missing_value_salt
  logical :: new_sim
  logical :: correct_thickness
  character(len=40) :: potemp_var, salin_var
  character(len=8)  :: laynum

  integer, parameter :: niter=10   ! number of iterations for t/s adjustment to layer density
  logical            :: adjust_temperature = .true.  ! fit t/s to target densities
  real, parameter    :: missing_value = -1.e20
  real, parameter    :: temp_land_fill = 0.0, salt_land_fill = 35.0
  logical :: reentrant_x, tripolar_n,dbg
  logical :: debug = .false.  ! manually set this to true for verbose output


  !data arrays
  real, dimension(:), allocatable :: z_edges_in, z_in, Rb
  real, dimension(:,:,:), allocatable, target :: temp_z, salt_z, mask_z
  real, dimension(:,:,:), allocatable :: rho_z
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: zi
  real, dimension(SZI_(G),SZJ_(G))  :: nlevs
  real, dimension(SZI_(G))   :: press


  ! Local variables for ALE remapping
  real, dimension(:), allocatable :: hTarget
  real, dimension(:,:), allocatable :: area_shelf_h
  real, dimension(:,:), allocatable, target  :: frac_shelf_h
  real, dimension(:,:,:), allocatable :: tmpT1dIn, tmpS1dIn, h1, tmp_mask_in
  real :: zTopOfCell, zBottomOfCell
  type(regridding_CS) :: regridCS ! Regridding parameters and work arrays
  type(remapping_CS) :: remapCS ! Remapping parameters and work arrays

  logical :: homogenize, useALEremapping, remap_full_column, remap_general, remap_old_alg
  logical :: use_ice_shelf
  character(len=10) :: remappingScheme
  real :: tempAvg, saltAvg
  integer :: nPoints, ans
  integer :: id_clock_routine, id_clock_read, id_clock_interp, id_clock_fill, id_clock_ALE

  id_clock_routine = cpu_clock_id('(Initialize from Z)', grain=CLOCK_ROUTINE)
  id_clock_ALE = cpu_clock_id('(Initialize from Z) ALE', grain=CLOCK_LOOP)

  call cpu_clock_begin(id_clock_routine)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  PI_180=atan(1.0)/45.

  call callTree_enter(trim(mod)//"(), MOM_state_initialization.F90")
  call log_version(PF, mod, version, "")

  new_sim = .false.
  if ((dirs%input_filename(1:1) == 'n') .and. &
       (LEN_TRIM(dirs%input_filename) == 1)) new_sim = .true.

  inputdir = "." ;  call get_param(PF, mod, "INPUTDIR", inputdir)
  inputdir = slasher(inputdir)

  eos => tv%eqn_of_state

! call mpp_get_compute_domain(G%domain%mpp_domain,isc,iec,jsc,jec)

  reentrant_x = .false. ; call get_param(PF, mod, "REENTRANT_X", reentrant_x,default=.true.)
  tripolar_n = .false. ;  call get_param(PF, mod, "TRIPOLAR_N", tripolar_n, default=.false.)
  call get_param(PF, mod, "MINIMUM_DEPTH", min_depth, default=0.0)

  call get_param(PF, mod, "NKML",nkml,default=0)
  call get_param(PF, mod, "NKBL",nkbl,default=0)

  call get_param(PF, mod, "TEMP_SALT_Z_INIT_FILE",filename, &
                 "The name of the z-space input file used to initialize \n"//&
                 "the layer thicknesses, temperatures and salinities.", &
                 default="temp_salt_z.nc")
  filename = trim(inputdir)//trim(filename)
  call get_param(PF, mod, "Z_INIT_FILE_PTEMP_VAR", potemp_var, &
                 "The name of the potential temperature variable in \n"//&
                 "TEMP_SALT_Z_INIT_FILE.", default="ptemp")
  call get_param(PF, mod, "Z_INIT_FILE_SALT_VAR", salin_var, &
                 "The name of the salinity variable in \n"//&
                 "TEMP_SALT_Z_INIT_FILE.", default="salt")
  call get_param(PF, mod, "Z_INIT_HOMOGENIZE", homogenize, &
                 "If True, then horizontally homogenize the interpolated \n"//&
                 "initial conditions.", default=.false.)
  call get_param(PF, mod, "Z_INIT_ALE_REMAPPING", useALEremapping, &
                 "If True, then remap straight to model coordinate from file.",&
                 default=.false.)
  call get_param(PF, mod, "Z_INIT_REMAPPING_SCHEME", remappingScheme, &
                 "The remapping scheme to use if using Z_INIT_ALE_REMAPPING\n"//&
                 "is True.", default="PPM_IH4")
  call get_param(PF, mod, "Z_INIT_REMAP_GENERAL", remap_general, &
                 "If false, only initializes to z* coordinates.\n"//&
                 "If true, allows initialization directly to general coordinates.",&
                 default=.false.)
  call get_param(PF, mod, "Z_INIT_REMAP_FULL_COLUMN", remap_full_column, &
                 "If false, only reconstructs profiles for valid data points.\n"//&
                 "If true, inserts vanished layers below the valid data.",&
                 default=remap_general)
  call get_param(PF, mod, "Z_INIT_REMAP_OLD_ALG", remap_old_alg, &
                 "If false, uses the preferred remapping algorithm for initialization.\n"//&
                 "If true, use an older, less robust algorithm for remapping.",&
                 default=.true.)

!   Read input grid coordinates for temperature and salinity field
!   in z-coordinate dataset. The file is REQUIRED to contain the
!   following:
!
!   dimension variables:
!            lon (degrees_E), lat (degrees_N), depth(meters)
!   variables:
!            ptemp(lon,lat,depth) : degC, potential temperature
!            salt (lon,lat,depth) : PSU, salinity
!
!   The first record will be read if there are multiple time levels.
!   The observation grid MUST tile the model grid. If the model grid extends
!   to the North/South Pole past the limits of the input data, they are extrapolated using the average
!   value at the northernmost/southernmost latitude.


  call horiz_interp_and_extrap_tracer(filename, potemp_var,1.0,1, &
       G, temp_z, mask_z, z_in, z_edges_in, missing_value_temp, reentrant_x, tripolar_n, homogenize)

  call horiz_interp_and_extrap_tracer(filename, salin_var,1.0,1, &
       G, salt_z, mask_z, z_in, z_edges_in, missing_value_salt, reentrant_x, tripolar_n, homogenize)

  kd = size(z_in,1)

  allocate(rho_z(isd:ied,jsd:jed,kd))
  allocate(area_shelf_h(isd:ied,jsd:jed))
  allocate(frac_shelf_h(isd:ied,jsd:jed))

  press(:)=tv%p_ref

  !Convert T&S to Absolute Salinity and Conservative Temperature if using TEOS10 or NEMO
  call convert_temp_salt_for_TEOS10(temp_z,salt_z, press, G, kd, mask_z, eos)

  do k=1, kd
    do j=js,je
      call calculate_density(temp_z(:,j,k),salt_z(:,j,k), press, rho_z(:,j,k), is, ie, eos)
    enddo
  enddo ! kd

  call pass_var(temp_z,G%Domain)
  call pass_var(salt_z,G%Domain)
  call pass_var(mask_z,G%Domain)
  call pass_var(rho_z,G%Domain)

  ! This is needed for building an ALE grid under ice shelves
  call get_param(PF, mod, "ICE_SHELF", use_ice_shelf, default=.false.)
  if (use_ice_shelf) then
     call get_param(PF, mod, "ICE_THICKNESS_FILE", ice_shelf_file, &
                    "The file from which the ice bathymetry and area are read.", &
                    fail_if_missing=.true.)
     filename = trim(inputdir)//trim(ice_shelf_file)
     call log_param(PF, mod, "INPUTDIR/THICKNESS_FILE", filename)
     call get_param(PF, mod, "ICE_AREA_VARNAME", area_varname, &
                    "The name of the area variable in ICE_THICKNESS_FILE.", &
                    fail_if_missing=.true.)
     if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       "MOM_temp_salt_initialize_from_Z: Unable to open "//trim(filename))

     call read_data(filename,trim(area_varname),area_shelf_h,domain=G%Domain%mpp_domain)

     ! initialize frac_shelf_h with zeros (open water everywhere)
     frac_shelf_h(:,:) = 0.0
     ! compute fractional ice shelf coverage of h
     do j=jsd,jed ; do i=isd,ied
         if (G%areaT(i,j) > 0.0) &
           frac_shelf_h(i,j) = area_shelf_h(i,j) / G%areaT(i,j)
     enddo ; enddo
     ! pass to the pointer
     shelf_area => frac_shelf_h

  endif

! Done with horizontal interpolation.
! Now remap to model coordinates
  if (useALEremapping) then
    call cpu_clock_begin(id_clock_ALE)
    ! The regridding tools (grid generation) are coded to work on model arrays of the same
    ! vertical shape. We need to re-write the regridding if the model has fewer layers
    ! than the data. -AJA
    if (kd>nz) call MOM_error(FATAL,"MOM_initialize_state, MOM_temp_salt_initialize_from_Z(): "//&
         "Data has more levels than the model - this has not been coded yet!")
    ! Build the source grid and copy data onto model-shaped arrays with vanished layers
    allocate( tmp_mask_in(isd:ied,jsd:jed,nz) ) ; tmp_mask_in(:,:,:) = 0.
    allocate( h1(isd:ied,jsd:jed,nz) ) ; h1(:,:,:) = 0.
    allocate( tmpT1dIn(isd:ied,jsd:jed,nz) ) ; tmpT1dIn(:,:,:) = 0.
    allocate( tmpS1dIn(isd:ied,jsd:jed,nz) ) ; tmpS1dIn(:,:,:) = 0.
    do j = js, je ; do i = is, ie
      if (G%mask2dT(i,j)>0.) then
        zTopOfCell = 0. ; zBottomOfCell = 0. ; nPoints = 0
        tmp_mask_in(i,j,1:kd) = mask_z(i,j,:)
        do k = 1, nz
          if (tmp_mask_in(i,j,k)>0. .and. k<=kd) then
            zBottomOfCell = -min( z_edges_in(k+1), G%bathyT(i,j) )
            tmpT1dIn(i,j,k) = temp_z(i,j,k)
            tmpS1dIn(i,j,k) = salt_z(i,j,k)
          elseif (k>1) then
            zBottomOfCell = -G%bathyT(i,j)
            tmpT1dIn(i,j,k) = tmpT1dIn(i,j,k-1)
            tmpS1dIn(i,j,k) = tmpS1dIn(i,j,k-1)
          else ! This next block should only ever be reached over land
            tmpT1dIn(i,j,k) = -99.9
            tmpS1dIn(i,j,k) = -99.9
          endif
          h1(i,j,k) = zTopOfCell - zBottomOfCell
          if (h1(i,j,k)>0.) nPoints = nPoints + 1
          zTopOfCell = zBottomOfCell ! Bottom becomes top for next value of k
        enddo
        h1(i,j,kd) = h1(i,j,kd) + ( zTopOfCell + G%bathyT(i,j) ) ! In case data is deeper than model
      endif ! mask2dT
    enddo ; enddo
    deallocate( tmp_mask_in )
    call pass_var(h1, G%Domain)
    call pass_var(tmpT1dIn, G%Domain)
    call pass_var(tmpS1dIn, G%Domain)

    ! Build the target grid (and set the model thickness to it)
    ! This call can be more general but is hard-coded for z* coordinates...  ????
    call ALE_initRegridding( GV, G%max_depth, PF, mod, regridCS ) ! sets regridCS

    if (.not. remap_general) then
      ! This is the old way of initializing to z* coordinates only
      allocate( hTarget(nz) )
      hTarget = getCoordinateResolution( regridCS )
      do j = js, je ; do i = is, ie
        h(i,j,:) = 0.
        if (G%mask2dT(i,j)>0.) then
          ! Build the target grid combining hTarget and topography
          zTopOfCell = 0. ; zBottomOfCell = 0.
          do k = 1, nz
            zBottomOfCell = max( zTopOfCell - hTarget(k), -G%bathyT(i,j) )
            h(i,j,k) = zTopOfCell - zBottomOfCell
            zTopOfCell = zBottomOfCell ! Bottom becomes top for next value of k
          enddo
        else
          h(i,j,:) = 0.
        endif ! mask2dT
      enddo ; enddo
      call pass_var(h, G%Domain)
      deallocate( hTarget )
    endif

    ! Now remap from source grid to target grid
    call initialize_remapping( remapCS, remappingScheme, boundary_extrapolation=.false. ) ! Reconstruction parameters
    if (remap_general) then
      call set_regrid_params( regridCS, min_thickness=0. )
      h(:,:,:) = h1(:,:,:) ; tv%T(:,:,:) = tmpT1dIn(:,:,:) ; tv%S(:,:,:) = tmpS1dIn(:,:,:)
      do j = js, je ; do i = is, ie
        if (G%mask2dT(i,j)==0.) then ! Ensure there are no nonsense values on land
          h(i,j,:) = 0. ; tv%T(i,j,:) = 0. ; tv%S(i,j,:) = 0.
        endif
      enddo ; enddo
      call pass_var(h, G%Domain)    ! Regridding might eventually use spatial information and
      call pass_var(tv%T, G%Domain) ! thus needs to be up to date in the halo regions even though
      call pass_var(tv%S, G%Domain) ! ALE_build_grid() only updates h on the computational domain.

      if (use_ice_shelf) then
         call ALE_build_grid( G, GV, regridCS, remapCS, h, tv, .true., shelf_area)
      else
         call ALE_build_grid( G, GV, regridCS, remapCS, h, tv, .true. )
      endif
    endif
    call ALE_remap_scalar( remapCS, G, GV, nz, h1, tmpT1dIn, h, tv%T, all_cells=remap_full_column, old_remap=remap_old_alg )
    call ALE_remap_scalar( remapCS, G, GV, nz, h1, tmpS1dIn, h, tv%S, all_cells=remap_full_column, old_remap=remap_old_alg )
    deallocate( h1 )
    deallocate( tmpT1dIn )
    deallocate( tmpS1dIn )

    call cpu_clock_end(id_clock_ALE)

  else ! remap to isopycnal layer space

! next find interface positions using local arrays
! nlevs contains the number of valid data points in each column
    nlevs = sum(mask_z,dim=3)

! Rb contains the layer interface densities
    allocate(Rb(nz+1))
    do k=2,nz ; Rb(k)=0.5*(GV%Rlay(k-1)+GV%Rlay(k)) ; enddo
    Rb(1) = 0.0 ;  Rb(nz+1) = 2.0*GV%Rlay(nz) - GV%Rlay(nz-1)

    zi(is:ie,js:je,:) = find_interfaces(rho_z(is:ie,js:je,:), z_in, Rb, G%bathyT(is:ie,js:je), &
                         nlevs(is:ie,js:je), nkml, nkbl, min_depth)

    call get_param(PF, mod, "ADJUST_THICKNESS", correct_thickness, &
                 "If true, all mass below the bottom removed if the \n"//&
                 "topography is shallower than the thickness input file \n"//&
                 "would indicate.", default=.false.)

    call get_param(PF, mod, "FIT_TO_TARGET_DENSITY_IC", adjust_temperature, &
                 "If true, all the interior layers are adjusted to \n"//&
                 "their target densities using mostly temperature \n"//&
                 "This approach can be problematic, particularly in the \n"//&
                 "high latitudes.", default=.true.)

    if (correct_thickness) then
      call adjustEtaToFitBathymetry(G, GV, zi, h)
    else
      do k=nz,1,-1 ; do j=js,je ; do i=is,ie
        if (zi(i,j,K) < (zi(i,j,K+1) + GV%Angstrom_z)) then
          zi(i,j,K) = zi(i,j,K+1) + GV%Angstrom_z
          h(i,j,k) = GV%Angstrom_z
        else
          h(i,j,k) = zi(i,j,K) - zi(i,j,K+1)
        endif
      enddo ; enddo ; enddo
      inconsistent=0
      do j=js,je ; do i=is,ie
        if (abs(zi(i,j,nz+1) + G%bathyT(i,j)) > 1.0) &
          inconsistent = inconsistent + 1
      enddo ; enddo
      call sum_across_PEs(inconsistent)

      if ((inconsistent > 0) .and. (is_root_pe())) then
        write(mesg,'("Thickness initial conditions are inconsistent ",'// &
                 '"with topography in ",I5," places.")') inconsistent
        call MOM_error(WARNING, mesg)
      endif
    endif

    tv%T(is:ie,js:je,:) = tracer_z_init(temp_z(is:ie,js:je,:),-1.0*z_edges_in,zi(is:ie,js:je,:), &
                                        nkml,nkbl,missing_value,G%mask2dT(is:ie,js:je),nz, &
                                        nlevs(is:ie,js:je),dbg,idbg,jdbg)
    tv%S(is:ie,js:je,:) = tracer_z_init(salt_z(is:ie,js:je,:),-1.0*z_edges_in,zi(is:ie,js:je,:), &
                                        nkml,nkbl,missing_value,G%mask2dT(is:ie,js:je),nz, &
                                        nlevs(is:ie,js:je))

    do k=1,nz

       nPoints = 0 ; tempAvg = 0. ; saltAvg = 0.
       do j=js,je
          do i=is,ie
             if (G%mask2dT(i,j) .ge. 1.0) then
                nPoints = nPoints + 1
                tempAvg = tempAvg + tv%T(i,j,k)
                saltAvg =saltAvg + tv%S(i,j,k)
             endif
          enddo
       enddo

    ! Horizontally homogenize data to produce perfectly "flat" initial conditions
       if (homogenize) then
          call sum_across_PEs(nPoints)
          call sum_across_PEs(tempAvg)
          call sum_across_PEs(saltAvg)
          if (nPoints>0) then
             tempAvg = tempAvg/real(nPoints)
             saltAvg = saltAvg/real(nPoints)
          endif
          tv%T(:,:,k) = tempAvg
          tv%S(:,:,k) = saltAvg
       endif

    enddo

  endif ! useALEremapping

! Fill land values
  do k=1,nz ; do j=js,je ; do i=is,ie
    if (tv%T(i,j,k) == missing_value) then
      tv%T(i,j,k)=temp_land_fill
      tv%S(i,j,k)=salt_land_fill
    endif
  enddo ; enddo ; enddo

! Finally adjust to target density
  ks=max(0,nkml)+max(0,nkbl)+1

  if (adjust_temperature .and. .not. useALEremapping) then
    call determine_temperature(tv%T(is:ie,js:je,:), tv%S(is:ie,js:je,:), &
            GV%Rlay(1:nz), tv%p_ref, niter, missing_value, h(is:ie,js:je,:), ks, eos)

  endif

  deallocate(z_in,z_edges_in,temp_z,salt_z,mask_z)

  call pass_var(h, G%Domain)
  call pass_var(tv%T, G%Domain)
  call pass_var(tv%S, G%Domain)

  call callTree_leave(trim(mod)//'()')
  call cpu_clock_end(id_clock_routine)

end subroutine MOM_temp_salt_initialize_from_Z

!> Run simple unit tests
subroutine MOM_state_init_tests(G, GV, tv)
  type(ocean_grid_type),     intent(inout) :: G
  type(verticalGrid_type),   intent(in)    :: GV
  type(thermo_var_ptrs),     intent(in)    :: tv !< Thermodynamics structure
  ! Local variables
  integer, parameter :: nk=5
  real, dimension(nk) :: T, T_t, T_b, S, S_t, S_b, rho, h, z
  real, dimension(nk+1) :: e
  integer :: k
  real :: P_tot, P_t, P_b, z_out
  type(remapping_CS), pointer :: remap_CS => NULL()

  do k = 1, nk
    h(k) = 100.
  enddo
  e(1) = 0.
  do K = 1, nk
    e(K+1) = e(K) - h(k)
  enddo
  P_tot = 0.
  do k = 1, nk
    z(k) = 0.5 * ( e(K) + e(K+1) )
    T_t(k) = 20.+(0./500.)*e(k)
    T(k)   = 20.+(0./500.)*z(k)
    T_b(k) = 20.+(0./500.)*e(k+1)
    S_t(k) = 35.-(0./500.)*e(k)
    S(k)   = 35.+(0./500.)*z(k)
    S_b(k) = 35.-(0./500.)*e(k+1)
    call calculate_density(0.5*(T_t(k)+T_b(k)), 0.5*(S_t(k)+S_b(k)), -GV%Rho0*GV%g_Earth*z(k), rho(k), tv%eqn_of_state)
    P_tot = P_tot + GV%g_Earth * rho(k) * h(k)
  enddo

  P_t = 0.
  do k = 1, nk
    call find_depth_of_pressure_in_cell(T_t(k), T_b(k), S_t(k), S_b(k), e(K), e(K+1), &
                                        P_t, 0.5*P_tot, GV%Rho0, GV%g_Earth, tv%eqn_of_state, P_b, z_out)
    write(0,*) k,P_t,P_b,0.5*P_tot,e(K),e(K+1),z_out
    P_t = P_b
  enddo
  write(0,*) P_b,P_tot

  write(0,*) ''
  write(0,*) ' ==================================================================== '
  write(0,*) ''
  write(0,*) h
  call cut_off_column_top(nk, tv, GV%Rho0, GV%g_Earth, -e(nk+1), GV%Angstrom, &
               T, T_t, T_b, S, S_t, S_b, 0.5*P_tot, h, remap_CS)
  write(0,*) h

end subroutine MOM_state_init_tests

end module MOM_state_initialization
