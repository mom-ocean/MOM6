!> Initialization functions for state variables, u, v, h, T and S.
module MOM_state_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging, only : hchksum, qchksum, uvchksum
use MOM_density_integrals, only : int_specific_vol_dp
use MOM_density_integrals, only : find_depth_of_pressure_in_cell
use MOM_coms, only : max_across_PEs, min_across_PEs, reproducing_sum
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only :  CLOCK_ROUTINE, CLOCK_LOOP
use MOM_domains, only : pass_var, pass_vector, sum_across_PEs, broadcast
use MOM_domains, only : root_PE, To_All, SCALAR_PAIR, CGRID_NE, AGRID
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, param_file_type
use MOM_file_parser, only : log_version
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type, isPointInCell
use MOM_interface_heights, only : find_eta
use MOM_io, only : file_exists, field_size, MOM_read_data, MOM_read_vector, slasher
use MOM_open_boundary, only : ocean_OBC_type, open_boundary_init, set_tracer_data
use MOM_open_boundary, only : OBC_NONE, OBC_SIMPLE
use MOM_open_boundary, only : open_boundary_query
use MOM_open_boundary, only : set_tracer_data, initialize_segment_data
use MOM_open_boundary, only : open_boundary_test_extern_h
use MOM_open_boundary, only : fill_temp_salt_segments
use MOM_open_boundary, only : update_OBC_segment_data
!use MOM_open_boundary, only : set_3D_OBC_data
use MOM_grid_initialize, only : initialize_masks, set_grid_metrics
use MOM_restart, only : restore_state, determine_is_new_run, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, set_up_sponge_ML_density
use MOM_sponge, only : initialize_sponge, sponge_CS
use MOM_ALE_sponge, only : set_up_ALE_sponge_field, initialize_ALE_sponge, ALE_sponge_CS
use MOM_string_functions, only : uppercase, lowercase
use MOM_time_manager, only : time_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : setVerticalGridAxes, verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type, EOS_domain
use MOM_EOS, only : convert_temp_salt_for_TEOS10
use user_initialization, only : user_initialize_thickness, user_initialize_velocity
use user_initialization, only : user_init_temperature_salinity, user_set_OBC_data
use user_initialization, only : user_initialize_sponges
use DOME_initialization, only : DOME_initialize_thickness
use DOME_initialization, only : DOME_set_OBC_data
use DOME_initialization, only : DOME_initialize_sponges
use ISOMIP_initialization, only : ISOMIP_initialize_thickness
use ISOMIP_initialization, only : ISOMIP_initialize_sponges
use ISOMIP_initialization, only : ISOMIP_initialize_temperature_salinity
use RGC_initialization, only : RGC_initialize_sponges
use baroclinic_zone_initialization, only : baroclinic_zone_init_temperature_salinity
use benchmark_initialization, only : benchmark_initialize_thickness
use benchmark_initialization, only : benchmark_init_temperature_salinity
use Neverworld_initialization, only : Neverworld_initialize_thickness
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
use dumbbell_initialization, only : dumbbell_initialize_thickness
use dumbbell_initialization, only : dumbbell_initialize_temperature_salinity
use Phillips_initialization, only : Phillips_initialize_thickness
use Phillips_initialization, only : Phillips_initialize_velocity
use Phillips_initialization, only : Phillips_initialize_sponges
use Rossby_front_2d_initialization, only : Rossby_front_initialize_thickness
use Rossby_front_2d_initialization, only : Rossby_front_initialize_temperature_salinity
use Rossby_front_2d_initialization, only : Rossby_front_initialize_velocity
use SCM_CVMix_tests, only: SCM_CVMix_tests_TS_init
use dyed_channel_initialization, only : dyed_channel_set_OBC_tracer_data
use dyed_obcs_initialization, only : dyed_obcs_set_OBC_data
use supercritical_initialization, only : supercritical_set_OBC_data
use soliton_initialization, only : soliton_initialize_velocity
use soliton_initialization, only : soliton_initialize_thickness
use BFB_initialization, only : BFB_initialize_sponges_southonly
use dense_water_initialization, only : dense_water_initialize_TS
use dense_water_initialization, only : dense_water_initialize_sponges
use dumbbell_initialization, only : dumbbell_initialize_sponges
use MOM_tracer_Z_init, only : find_interfaces, tracer_Z_init_array, determine_temperature
use MOM_ALE, only : ALE_initRegridding, ALE_CS, ALE_initThicknessToCoord
use MOM_ALE, only : ALE_remap_scalar, ALE_build_grid, ALE_regrid_accelerated
use MOM_ALE, only : TS_PLM_edge_values
use MOM_regridding, only : regridding_CS, set_regrid_params, getCoordinateResolution
use MOM_regridding, only : regridding_main
use MOM_remapping, only : remapping_CS, initialize_remapping
use MOM_remapping, only : remapping_core_h
use MOM_horizontal_regridding, only : horiz_interp_and_extrap_tracer

implicit none ; private

#include <MOM_memory.h>

public MOM_initialize_state

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

character(len=40)  :: mdl = "MOM_state_initialization" !< This module's name.

contains

!> Initialize temporally evolving fields, either as initial
!! conditions or by reading them from a restart (or saves) file.
subroutine MOM_initialize_state(u, v, h, tv, Time, G, GV, US, PF, dirs, &
                                restart_CS, ALE_CSp, tracer_Reg, sponge_CSp, &
                                ALE_sponge_CSp, OBC, Time_in)
  type(ocean_grid_type),      intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                              intent(out)   :: u    !< The zonal velocity that is being
                                                    !! initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                              intent(out)   :: v    !< The meridional velocity that is being
                                                    !! initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                              intent(out)   :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),      intent(inout) :: tv   !< A structure pointing to various thermodynamic
                                                    !! variables
  type(time_type),            intent(inout) :: Time !< Time at the start of the run segment.
  type(param_file_type),      intent(in)    :: PF   !< A structure indicating the open file to parse
                                                    !! for model parameter values.
  type(directories),          intent(in)    :: dirs !< A structure containing several relevant
                                                    !! directory paths.
  type(MOM_restart_CS),       pointer       :: restart_CS !< A pointer to the restart control
                                                    !! structure.
  type(ALE_CS),               pointer       :: ALE_CSp !< The ALE control structure for remapping
  type(tracer_registry_type), pointer       :: tracer_Reg !< A pointer to the tracer registry
  type(sponge_CS),            pointer       :: sponge_CSp !< The layerwise sponge control structure.
  type(ALE_sponge_CS),        pointer       :: ALE_sponge_CSp !< The ALE sponge control structure.
  type(ocean_OBC_type),       pointer       :: OBC   !< The open boundary condition control structure.
  type(time_type), optional,  intent(in)    :: Time_in !< Time at the start of the run segment.
                                                     !! Time_in overrides any value set for Time.
  ! Local variables
  character(len=200) :: filename   ! The name of an input file.
  character(len=200) :: filename2  ! The name of an input files.
  character(len=200) :: inputdir   ! The directory where NetCDF input files are.
  character(len=200) :: config
  real :: H_rescale   ! A rescaling factor for thicknesses from the representation in
                      ! a restart file to the internal representation in this run.
  real :: vel_rescale ! A rescaling factor for velocities from the representation in
                      ! a restart file to the internal representation in this run.
  real :: dt          ! The baroclinic dynamics timestep for this run [T ~> s].
  logical :: from_Z_file, useALE
  logical :: new_sim
  integer :: write_geom
  logical :: use_temperature, use_sponge, use_OBC
  logical :: use_EOS     ! If true, density is calculated from T & S using an equation of state.
  logical :: depress_sfc ! If true, remove the mass that would be displaced
                         ! by a large surface pressure by squeezing the column.
  logical :: trim_ic_for_p_surf ! If true, remove the mass that would be displaced
                         ! by a large surface pressure, such as with an ice sheet.
  logical :: regrid_accelerate
  integer :: regrid_iterations
!  logical :: Analytic_FV_PGF, obsol_test
  logical :: convert
  logical :: just_read  ! If true, only read the parameters because this
                        ! is a run from a restart file; this option
                        ! allows the use of Fatal unused parameters.
  type(EOS_type), pointer :: eos => NULL()
  logical :: debug      ! If true, write debugging output.
  logical :: debug_obc  ! If true, do debugging calls related to OBCs.
  logical :: debug_layers = .false.
  character(len=80) :: mesg
! This include declares and sets the variable "version".
#include "version_variable.h"
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call callTree_enter("MOM_initialize_state(), MOM_state_initialization.F90")
  call log_version(PF, mdl, version, "")
  call get_param(PF, mdl, "DEBUG", debug, default=.false.)
  call get_param(PF, mdl, "DEBUG_OBC", debug_obc, default=.false.)

  new_sim = determine_is_new_run(dirs%input_filename, dirs%restart_input_dir, &
                                 G, restart_CS)
  just_read = .not.new_sim

  call get_param(PF, mdl, "INPUTDIR", inputdir, &
         "The directory in which input files are found.", default=".")
  inputdir = slasher(inputdir)

  use_temperature = associated(tv%T)
  useALE = associated(ALE_CSp)
  use_EOS = associated(tv%eqn_of_state)
  use_OBC = associated(OBC)
  if (use_EOS) eos => tv%eqn_of_state

  !====================================================================
  !    Initialize temporally evolving fields, either as initial
  !  conditions or by reading them from a restart (or saves) file.
  !====================================================================

  if (new_sim) then
    call MOM_mesg("Run initialized internally.", 3)

    if (present(Time_in)) Time = Time_in
    ! Otherwise leave Time at its input value.

    ! This initialization should not be needed. Certainly restricting it
    ! to the computational domain helps detect possible uninitialized
    ! data in halos which should be covered by the pass_var(h) later.
    !do k = 1, nz; do j = js, je; do i = is, ie
    !  h(i,j,k) = 0.
    !enddo
  endif

  ! The remaining initialization calls are done, regardless of whether the
  ! fields are actually initialized here (if just_read=.false.) or whether it
  ! is just to make sure that all valid parameters are read to enable the
  ! detection of unused parameters.
  call get_param(PF, mdl, "INIT_LAYERS_FROM_Z_FILE", from_Z_file, &
             "If true, initialize the layer thicknesses, temperatures, "//&
             "and salinities from a Z-space file on a latitude-longitude "//&
             "grid.", default=.false., do_not_log=just_read)

  if (from_Z_file) then
    ! Initialize thickness and T/S from z-coordinate data in a file.
    if (.NOT.use_temperature) call MOM_error(FATAL,"MOM_initialize_state : "//&
       "use_temperature must be true if INIT_LAYERS_FROM_Z_FILE is true")

    call MOM_temp_salt_initialize_from_Z(h, tv, G, GV, US, PF, just_read_params=just_read)

  else
    ! Initialize thickness, h.
    call get_param(PF, mdl, "THICKNESS_CONFIG", config, &
             "A string that determines how the initial layer "//&
             "thicknesses are specified for a new run: \n"//&
             " \t file - read interface heights from the file specified \n"//&
             " \t thickness_file - read thicknesses from the file specified \n"//&
             " \t\t by (THICKNESS_FILE).\n"//&
             " \t coord - determined by ALE coordinate.\n"//&
             " \t uniform - uniform thickness layers evenly distributed \n"//&
             " \t\t between the surface and MAXIMUM_DEPTH. \n"//&
             " \t list - read a list of positive interface depths. \n"//&
             " \t DOME - use a slope and channel configuration for the \n"//&
             " \t\t DOME sill-overflow test case. \n"//&
             " \t ISOMIP - use a configuration for the \n"//&
             " \t\t ISOMIP test case. \n"//&
             " \t benchmark - use the benchmark test case thicknesses. \n"//&
             " \t Neverworld - use the Neverworld test case thicknesses. \n"//&
             " \t search - search a density profile for the interface \n"//&
             " \t\t densities. This is not yet implemented. \n"//&
             " \t circle_obcs - the circle_obcs test case is used. \n"//&
             " \t DOME2D - 2D version of DOME initialization. \n"//&
             " \t adjustment2d - 2D lock exchange thickness ICs. \n"//&
             " \t sloshing - sloshing gravity thickness ICs. \n"//&
             " \t seamount - no motion test with seamount ICs. \n"//&
             " \t dumbbell - sloshing channel ICs. \n"//&
             " \t soliton - Equatorial Rossby soliton. \n"//&
             " \t rossby_front - a mixed layer front in thermal wind balance.\n"//&
             " \t USER - call a user modified routine.", &
             default="uniform", do_not_log=just_read)
    select case (trim(config))
       case ("file")
         call initialize_thickness_from_file(h, G, GV, US, PF, .false., just_read_params=just_read)
       case ("thickness_file")
         call initialize_thickness_from_file(h, G, GV, US, PF, .true., just_read_params=just_read)
       case ("coord")
         if (new_sim .and. useALE) then
           call ALE_initThicknessToCoord( ALE_CSp, G, GV, h )
         elseif (new_sim) then
           call MOM_error(FATAL, "MOM_initialize_state: USE_REGRIDDING must be True "//&
                                 "for THICKNESS_CONFIG of 'coord'")
         endif
       case ("uniform"); call initialize_thickness_uniform(h, G, GV, PF, &
                                  just_read_params=just_read)
       case ("list"); call initialize_thickness_list(h, G, GV, US, PF, &
                                  just_read_params=just_read)
       case ("DOME"); call DOME_initialize_thickness(h, G, GV, PF, &
                               just_read_params=just_read)
       case ("ISOMIP"); call ISOMIP_initialize_thickness(h, G, GV, US, PF, tv, &
                                 just_read_params=just_read)
       case ("benchmark"); call benchmark_initialize_thickness(h, G, GV, US, PF, &
                                    tv%eqn_of_state, tv%P_Ref, just_read_params=just_read)
       case ("Neverwoorld","Neverland"); call Neverworld_initialize_thickness(h, G, GV, US, PF, &
                                 tv%eqn_of_state, tv%P_Ref)
       case ("search"); call initialize_thickness_search
       case ("circle_obcs"); call circle_obcs_initialize_thickness(h, G, GV, PF, &
                                      just_read_params=just_read)
       case ("lock_exchange"); call lock_exchange_initialize_thickness(h, G, GV, US, &
                                        PF, just_read_params=just_read)
       case ("external_gwave"); call external_gwave_initialize_thickness(h, G, GV, US, &
                                         PF, just_read_params=just_read)
       case ("DOME2D"); call DOME2d_initialize_thickness(h, G, GV, US, PF, &
                                 just_read_params=just_read)
       case ("adjustment2d"); call adjustment_initialize_thickness(h, G, GV, US, &
                                       PF, just_read_params=just_read)
       case ("sloshing"); call sloshing_initialize_thickness(h, G, GV, US, PF, &
                                   just_read_params=just_read)
       case ("seamount"); call seamount_initialize_thickness(h, G, GV, US, PF, &
                                   just_read_params=just_read)
       case ("dumbbell"); call dumbbell_initialize_thickness(h, G, GV, US, PF, &
                                   just_read_params=just_read)
       case ("soliton"); call soliton_initialize_thickness(h, G, GV, US)
       case ("phillips"); call Phillips_initialize_thickness(h, G, GV, US, PF, &
                                   just_read_params=just_read)
       case ("rossby_front"); call Rossby_front_initialize_thickness(h, G, GV, US, &
                                       PF, just_read_params=just_read)
       case ("USER"); call user_initialize_thickness(h, G, GV, PF, &
                               just_read_params=just_read)
       case default ; call MOM_error(FATAL,  "MOM_initialize_state: "//&
            "Unrecognized layer thickness configuration "//trim(config))
    end select

    ! Initialize temperature and salinity (T and S).
    if ( use_temperature ) then
      call get_param(PF, mdl, "TS_CONFIG", config, &
             "A string that determines how the initial tempertures "//&
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
             " \t adjustment2d - 2d lock exchange T/S ICs. \n"//&
             " \t sloshing - sloshing mode T/S ICs. \n"//&
             " \t seamount - no motion test with seamount ICs. \n"//&
             " \t dumbbell - sloshing channel ICs. \n"//&
             " \t rossby_front - a mixed layer front in thermal wind balance.\n"//&
             " \t SCM_CVMix_tests - used in the SCM CVMix tests.\n"//&
             " \t USER - call a user modified routine.", &
             fail_if_missing=new_sim, do_not_log=just_read)
!            " \t baroclinic_zone - an analytic baroclinic zone. \n"//&
      select case (trim(config))
        case ("fit"); call initialize_temp_salt_fit(tv%T, tv%S, G, GV, US, PF, &
                               eos, tv%P_Ref, just_read_params=just_read)
        case ("file"); call initialize_temp_salt_from_file(tv%T, tv%S, G, &
                                PF, just_read_params=just_read)
        case ("benchmark"); call benchmark_init_temperature_salinity(tv%T, tv%S, &
                                     G, GV, US, PF, eos, tv%P_Ref, just_read_params=just_read)
        case ("TS_profile") ; call initialize_temp_salt_from_profile(tv%T, tv%S, &
                                       G, PF, just_read_params=just_read)
        case ("linear"); call initialize_temp_salt_linear(tv%T, tv%S, G, PF, &
                                  just_read_params=just_read)
        case ("DOME2D"); call DOME2d_initialize_temperature_salinity ( tv%T, &
                                  tv%S, h, G, GV, PF, eos, just_read_params=just_read)
        case ("ISOMIP"); call ISOMIP_initialize_temperature_salinity ( tv%T, &
                                  tv%S, h, G, GV, US, PF, eos, just_read_params=just_read)
        case ("adjustment2d"); call adjustment_initialize_temperature_salinity ( tv%T, &
                                        tv%S, h, G, GV, PF, eos, just_read_params=just_read)
        case ("baroclinic_zone"); call baroclinic_zone_init_temperature_salinity( tv%T, &
                                           tv%S, h, G, GV, US, PF, just_read_params=just_read)
        case ("sloshing"); call sloshing_initialize_temperature_salinity(tv%T, &
                                    tv%S, h, G, GV, PF, eos, just_read_params=just_read)
        case ("seamount"); call seamount_initialize_temperature_salinity(tv%T, &
                                    tv%S, h, G, GV, PF, eos, just_read_params=just_read)
        case ("dumbbell"); call dumbbell_initialize_temperature_salinity(tv%T, &
                                    tv%S, h, G, GV, PF, eos, just_read_params=just_read)
        case ("rossby_front"); call Rossby_front_initialize_temperature_salinity ( tv%T, &
                                        tv%S, h, G, GV, PF, eos, just_read_params=just_read)
        case ("SCM_CVMix_tests"); call SCM_CVMix_tests_TS_init(tv%T, tv%S, h, &
                                           G, GV, US, PF, just_read_params=just_read)
        case ("dense"); call dense_water_initialize_TS(G, GV, PF, eos, tv%T, tv%S, &
                                 h, just_read_params=just_read)
        case ("USER"); call user_init_temperature_salinity(tv%T, tv%S, G, PF, eos, &
                                just_read_params=just_read)
        case default ; call MOM_error(FATAL,  "MOM_initialize_state: "//&
               "Unrecognized Temp & salt configuration "//trim(config))
      end select
    endif
  endif  ! not from_Z_file.
  if (use_temperature .and. use_OBC) &
    call fill_temp_salt_segments(G, OBC, tv)

  ! The thicknesses in halo points might be needed to initialize the velocities.
  if (new_sim) call pass_var(h, G%Domain)

  ! Initialize velocity components, u and v
  call get_param(PF, mdl, "VELOCITY_CONFIG", config, &
       "A string that determines how the initial velocities "//&
       "are specified for a new run: \n"//&
       " \t file - read velocities from the file specified \n"//&
       " \t\t by (VELOCITY_FILE). \n"//&
       " \t zero - the fluid is initially at rest. \n"//&
       " \t uniform - the flow is uniform (determined by\n"//&
       " \t\t parameters INITIAL_U_CONST and INITIAL_V_CONST).\n"//&
       " \t rossby_front - a mixed layer front in thermal wind balance.\n"//&
       " \t soliton - Equatorial Rossby soliton.\n"//&
       " \t USER - call a user modified routine.", default="zero", &
       do_not_log=just_read)
  select case (trim(config))
     case ("file"); call initialize_velocity_from_file(u, v, G, US, PF, &
                             just_read_params=just_read)
     case ("zero"); call initialize_velocity_zero(u, v, G, PF, &
                             just_read_params=just_read)
     case ("uniform"); call initialize_velocity_uniform(u, v, G, US, PF, &
                                just_read_params=just_read)
     case ("circular"); call initialize_velocity_circular(u, v, G, US, PF, &
                                 just_read_params=just_read)
     case ("phillips"); call Phillips_initialize_velocity(u, v, G, GV, US, PF, &
                                 just_read_params=just_read)
     case ("rossby_front"); call Rossby_front_initialize_velocity(u, v, h, &
                                     G, GV, US, PF, just_read_params=just_read)
     case ("soliton"); call soliton_initialize_velocity(u, v, h, G, US)
     case ("USER"); call user_initialize_velocity(u, v, G, US, PF, &
                             just_read_params=just_read)
     case default ; call MOM_error(FATAL,  "MOM_initialize_state: "//&
          "Unrecognized velocity configuration "//trim(config))
  end select

  if (new_sim) call pass_vector(u, v, G%Domain)
  if (debug .and. new_sim) then
    call uvchksum("MOM_initialize_state [uv]", u, v, G%HI, haloshift=1, scale=US%m_s_to_L_T)
  endif

  ! Optionally convert the thicknesses from m to kg m-2.  This is particularly
  ! useful in a non-Boussinesq model.
  call get_param(PF, mdl, "CONVERT_THICKNESS_UNITS", convert, &
               "If true,  convert the thickness initial conditions from "//&
               "units of m to kg m-2 or vice versa, depending on whether "//&
               "BOUSSINESQ is defined. This does not apply if a restart "//&
               "file is read.", default=.not.GV%Boussinesq, do_not_log=just_read)

  if (new_sim .and. convert .and. .not.GV%Boussinesq) &
    ! Convert thicknesses from geomtric distances to mass-per-unit-area.
    call convert_thickness(h, G, GV, US, tv)

  ! Remove the mass that would be displaced by an ice shelf or inverse barometer.
  call get_param(PF, mdl, "DEPRESS_INITIAL_SURFACE", depress_sfc, &
               "If true,  depress the initial surface to avoid huge "//&
               "tsunamis when a large surface pressure is applied.", &
               default=.false., do_not_log=just_read)
  call get_param(PF, mdl, "TRIM_IC_FOR_P_SURF", trim_ic_for_p_surf, &
               "If true, cuts way the top of the column for initial conditions "//&
               "at the depth where the hydrostatic pressure matches the imposed "//&
               "surface pressure which is read from file.", default=.false., &
               do_not_log=just_read)
  if (depress_sfc .and. trim_ic_for_p_surf) call MOM_error(FATAL, "MOM_initialize_state: "//&
           "DEPRESS_INITIAL_SURFACE and TRIM_IC_FOR_P_SURF are exclusive and cannot both be True")
  if (new_sim .and. debug .and. (depress_sfc .or. trim_ic_for_p_surf)) &
    call hchksum(h, "Pre-depress: h ", G%HI, haloshift=1, scale=GV%H_to_m)
  if (depress_sfc) call depress_surface(h, G, GV, US, PF, tv, just_read_params=just_read)
  if (trim_ic_for_p_surf) call trim_for_ice(PF, G, GV, US, ALE_CSp, tv, h, just_read_params=just_read)

  ! Perhaps we want to run the regridding coordinate generator for multiple
  ! iterations here so the initial grid is consistent with the coordinate
  if (useALE) then
    call get_param(PF, mdl, "REGRID_ACCELERATE_INIT", regrid_accelerate, &
         "If true, runs REGRID_ACCELERATE_ITERATIONS iterations of the regridding "//&
         "algorithm to push the initial grid to be consistent with the initial "//&
         "condition. Useful only for state-based and iterative coordinates.", &
         default=.false., do_not_log=just_read)
    if (regrid_accelerate) then
      call get_param(PF, mdl, "REGRID_ACCELERATE_ITERATIONS", regrid_iterations, &
           "The number of regridding iterations to perform to generate "//&
           "an initial grid that is consistent with the initial conditions.", &
           default=1, do_not_log=just_read)

      call get_param(PF, mdl, "DT", dt, "Timestep", fail_if_missing=.true., scale=US%s_to_T)

      if (new_sim .and. debug) &
        call hchksum(h, "Pre-ALE_regrid: h ", G%HI, haloshift=1, scale=GV%H_to_m)
      call ALE_regrid_accelerated(ALE_CSp, G, GV, h, tv, regrid_iterations, u, v, OBC, tracer_Reg, &
                                  dt=dt, initial=.true.)
    endif
  endif
  ! This is the end of the block of code that might have initialized fields
  ! internally at the start of a new run.

  if (.not.new_sim) then ! This block restores the state from a restart file.
    !    This line calls a subroutine that reads the initial conditions
    !  from a previously generated file.
    call restore_state(dirs%input_filename, dirs%restart_input_dir, Time, &
                       G, restart_CS)
    if (present(Time_in)) Time = Time_in
    if ((GV%m_to_H_restart /= 0.0) .and. (GV%m_to_H_restart /= GV%m_to_H)) then
      H_rescale = GV%m_to_H / GV%m_to_H_restart
      do k=1,nz ; do j=js,je ; do i=is,ie ; h(i,j,k) = H_rescale * h(i,j,k) ; enddo ; enddo ; enddo
    endif
    if ( (US%s_to_T_restart * US%m_to_L_restart /= 0.0) .and. &
         ((US%m_to_L * US%s_to_T_restart) /= (US%m_to_L_restart * US%s_to_T)) ) then
      vel_rescale = (US%m_to_L * US%s_to_T_restart) /  (US%m_to_L_restart * US%s_to_T)
      do k=1,nz ; do j=jsd,jed ; do I=IsdB,IeDB ; u(I,j,k) = vel_rescale * u(I,j,k) ; enddo ; enddo ; enddo
      do k=1,nz ; do J=JsdB,JedB ; do i=isd,ied ; v(i,J,k) = vel_rescale * v(i,J,k) ; enddo ; enddo ; enddo
    endif
  endif

  if ( use_temperature ) then
    call pass_var(tv%T, G%Domain, complete=.false.)
    call pass_var(tv%S, G%Domain, complete=.false.)
  endif
  call pass_var(h, G%Domain)

  if (debug) then
    call hchksum(h, "MOM_initialize_state: h ", G%HI, haloshift=1, scale=GV%H_to_m)
    if ( use_temperature ) call hchksum(tv%T, "MOM_initialize_state: T ", G%HI, haloshift=1)
    if ( use_temperature ) call hchksum(tv%S, "MOM_initialize_state: S ", G%HI, haloshift=1)
    if ( use_temperature .and. debug_layers) then ; do k=1,nz
      write(mesg,'("MOM_IS: T[",I2,"]")') k
      call hchksum(tv%T(:,:,k), mesg, G%HI, haloshift=1)
      write(mesg,'("MOM_IS: S[",I2,"]")') k
      call hchksum(tv%S(:,:,k), mesg, G%HI, haloshift=1)
    enddo ; endif
  endif

  call get_param(PF, mdl, "SPONGE", use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified via SPONGE_CONFIG.", default=.false.)
  if ( use_sponge ) then
    call get_param(PF, mdl, "SPONGE_CONFIG", config, &
                 "A string that sets how the sponges are configured: \n"//&
                 " \t file - read sponge properties from the file \n"//&
                 " \t\t specified by (SPONGE_FILE).\n"//&
                 " \t ISOMIP - apply ale sponge in the ISOMIP case \n"//&
                 " \t RGC - apply sponge in the rotating_gravity_current case \n"//&
                 " \t DOME - use a slope and channel configuration for the \n"//&
                 " \t\t DOME sill-overflow test case. \n"//&
                 " \t BFB - Sponge at the southern boundary of the domain\n"//&
                 " \t\t for buoyancy-forced basin case.\n"//&
                 " \t USER - call a user modified routine.", default="file")
    select case (trim(config))
      case ("DOME"); call DOME_initialize_sponges(G, GV, US, tv, PF, sponge_CSp)
      case ("DOME2D"); call DOME2d_initialize_sponges(G, GV, US, tv, PF, useALE, &
                                                      sponge_CSp, ALE_sponge_CSp)
      case ("ISOMIP"); call ISOMIP_initialize_sponges(G, GV, US, tv, PF, useALE, &
                                                      sponge_CSp, ALE_sponge_CSp)
      case("RGC"); call RGC_initialize_sponges(G, GV, US, tv, u, v, PF, useALE, &
                                                     sponge_CSp, ALE_sponge_CSp)
      case ("USER"); call user_initialize_sponges(G, GV, use_temperature, tv, PF, sponge_CSp, h)
      case ("BFB"); call BFB_initialize_sponges_southonly(G, GV, US, use_temperature, tv, PF, &
                                                          sponge_CSp, h)
      case ("DUMBBELL"); call dumbbell_initialize_sponges(G, GV, US, tv, PF, useALE, &
                                                          sponge_CSp, ALE_sponge_CSp)
      case ("phillips"); call Phillips_initialize_sponges(G, GV, US, tv, PF, sponge_CSp, h)
      case ("dense"); call dense_water_initialize_sponges(G, GV, US, tv, PF, useALE, &
                                                          sponge_CSp, ALE_sponge_CSp)
      case ("file"); call initialize_sponges_file(G, GV, US, use_temperature, tv, PF, &
                                                  sponge_CSp, ALE_sponge_CSp, Time)
      case default ; call MOM_error(FATAL,  "MOM_initialize_state: "//&
             "Unrecognized sponge configuration "//trim(config))
    end select
  endif

  ! Reads OBC parameters not pertaining to the location of the boundaries
  call open_boundary_init(G, GV, US, PF, OBC, restart_CS)

  ! This controls user code for setting open boundary data
  if (associated(OBC)) then
     call initialize_segment_data(G, OBC, PF) !   call initialize_segment_data(G, OBC, param_file)
!     call open_boundary_config(G, US, PF, OBC)
    ! Call this once to fill boundary arrays from fixed values
    if (.not. OBC%needs_IO_for_data)  &
      call update_OBC_segment_data(G, GV, US, OBC, tv, h, Time)

    call get_param(PF, mdl, "OBC_USER_CONFIG", config, &
                 "A string that sets how the user code is invoked to set open boundary data: \n"//&
                 "   DOME - specified inflow on northern boundary\n"//&
                 "   dyed_channel - supercritical with dye on the inflow boundary\n"//&
                 "   dyed_obcs - circle_obcs with dyes on the open boundaries\n"//&
                 "   Kelvin - barotropic Kelvin wave forcing on the western boundary\n"//&
                 "   shelfwave - Flather with shelf wave forcing on western boundary\n"//&
                 "   supercritical - now only needed here for the allocations\n"//&
                 "   tidal_bay - Flather with tidal forcing on eastern boundary\n"//&
                 "   USER - user specified", default="none")
    if (trim(config) == "DOME") then
      call DOME_set_OBC_data(OBC, tv, G, GV, US, PF, tracer_Reg)
    elseif (trim(config) == "dyed_channel") then
      call dyed_channel_set_OBC_tracer_data(OBC, G, GV, PF, tracer_Reg)
      OBC%update_OBC = .true.
    elseif (trim(config) == "dyed_obcs") then
      call dyed_obcs_set_OBC_data(OBC, G, GV, PF, tracer_Reg)
    elseif (trim(config) == "Kelvin") then
      OBC%update_OBC = .true.
    elseif (trim(config) == "shelfwave") then
      OBC%update_OBC = .true.
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
    if (open_boundary_query(OBC, apply_open_OBC=.true.)) then
      call set_tracer_data(OBC, tv, h, G, PF, tracer_Reg)
    endif
  endif
! if (open_boundary_query(OBC, apply_nudged_OBC=.true.)) then
!   call set_3D_OBC_data(OBC, tv, h, G, PF, tracer_Reg)
! endif
  ! Still need a way to specify the boundary values
  if (debug.and.associated(OBC)) then
    call hchksum(G%mask2dT, 'MOM_initialize_state: mask2dT ', G%HI)
    call uvchksum('MOM_initialize_state: mask2dC[uv]', G%mask2dCu,  &
                  G%mask2dCv, G%HI)
    call qchksum(G%mask2dBu, 'MOM_initialize_state: mask2dBu ', G%HI)
  endif

  if (debug_OBC) call open_boundary_test_extern_h(G, GV, OBC, h)
  call callTree_leave('MOM_initialize_state()')

end subroutine MOM_initialize_state

!> Reads the layer thicknesses or interface heights from a file.
subroutine initialize_thickness_from_file(h, G, GV, US, param_file, file_has_thickness, &
                                          just_read_params)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h    !< The thickness that is being initialized [H ~> m or kg m-2].
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                               !! to parse for model parameter values.
  logical,                 intent(in)  :: file_has_thickness !< If true, this file contains layer
                                               !! thicknesses; otherwise it contains
                                               !! interface heights.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                               !! only read parameters without changing h.

  ! Local variables
  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1)  ! Interface heights, in depth units.
  integer :: inconsistent = 0
  logical :: correct_thickness
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40)  :: mdl = "initialize_thickness_from_file" ! This subroutine's name.
  character(len=200) :: filename, thickness_file, inputdir, mesg ! Strings for file/path
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) &
    call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".", do_not_log=just_read)
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "THICKNESS_FILE", thickness_file, &
                 "The name of the thickness file.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)

  filename = trim(inputdir)//trim(thickness_file)
  if (.not.just_read) call log_param(param_file, mdl, "INPUTDIR/THICKNESS_FILE", filename)

  if ((.not.just_read) .and. (.not.file_exists(filename, G%Domain))) call MOM_error(FATAL, &
         " initialize_thickness_from_file: Unable to open "//trim(filename))

  if (file_has_thickness) then
    !### Consider adding a parameter to use to rescale h.
    if (just_read) return ! All run-time parameters have been read, so return.
    call MOM_read_data(filename, "h", h(:,:,:), G%Domain, scale=GV%m_to_H)
  else
    call get_param(param_file, mdl, "ADJUST_THICKNESS", correct_thickness, &
                 "If true, all mass below the bottom removed if the "//&
                 "topography is shallower than the thickness input file "//&
                 "would indicate.", default=.false., do_not_log=just_read)
    if (just_read) return ! All run-time parameters have been read, so return.

    call MOM_read_data(filename, "eta", eta(:,:,:), G%Domain, scale=US%m_to_Z)

    if (correct_thickness) then
      call adjustEtaToFitBathymetry(G, GV, US, eta, h)
    else
      do k=nz,1,-1 ; do j=js,je ; do i=is,ie
        if (eta(i,j,K) < (eta(i,j,K+1) + GV%Angstrom_Z)) then
          eta(i,j,K) = eta(i,j,K+1) + GV%Angstrom_Z
          h(i,j,k) = GV%Angstrom_H
        else
          h(i,j,k) = GV%Z_to_H * (eta(i,j,K) - eta(i,j,K+1))
        endif
      enddo ; enddo ; enddo

      do j=js,je ; do i=is,ie
        if (abs(eta(i,j,nz+1) + G%bathyT(i,j)) > 1.0*US%m_to_Z) &
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
  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_thickness_from_file

!> Adjust interface heights to fit the bathymetry and diagnose layer thickness.
!!
!! If the bottom most interface is below the topography then the bottom-most
!! layers are contracted to GV%Angstrom_m.
!! If the bottom most interface is above the topography then the entire column
!! is dilated (expanded) to fill the void.
!!   @remark{There is a (hard-wired) "tolerance" parameter such that the
!! criteria for adjustment must equal or exceed 10cm.}
subroutine adjustEtaToFitBathymetry(G, GV, US, eta, h)
  type(ocean_grid_type),                      intent(in)    :: G   !< The ocean's grid structure
  type(verticalGrid_type),                    intent(in)    :: GV  !< The ocean's vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: eta !< Interface heights [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(inout) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  ! Local variables
  integer :: i, j, k, is, ie, js, je, nz, contractions, dilations
  real :: hTolerance = 0.1 !<  Tolerance to exceed adjustment criteria [Z ~> m]
  real :: hTmp, eTmp, dilate
  character(len=100) :: mesg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  hTolerance = 0.1*US%m_to_Z

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

  !   To preserve previous answers in non-Boussinesq cases, delay converting
  ! thicknesses to units of H until the end of this routine.
  do k=nz,1,-1 ; do j=js,je ; do i=is,ie
    ! Collapse layers to thinnest possible if the thickness less than
    ! the thinnest possible (or negative).
    if (eta(i,j,K) < (eta(i,j,K+1) + GV%Angstrom_Z)) then
      eta(i,j,K) = eta(i,j,K+1) + GV%Angstrom_Z
      h(i,j,k) = GV%Angstrom_Z
    else
      h(i,j,k) = (eta(i,j,K) - eta(i,j,K+1))
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
        do k=1,nz ; h(i,j,k) = (eta(i,j,1) + G%bathyT(i,j)) / real(nz) ; enddo
      else
        dilate = (eta(i,j,1) + G%bathyT(i,j)) / (eta(i,j,1) - eta(i,j,nz+1))
        do k=1,nz ; h(i,j,k) = h(i,j,k) * dilate ; enddo
      endif
      do k=nz,2,-1 ; eta(i,j,K) = eta(i,j,K+1) + h(i,j,k) ; enddo
    endif
  enddo ; enddo

  ! Now convert thicknesses to units of H.
  do k=1,nz ; do j=js,je ; do i=is,ie
    h(i,j,k) = h(i,j,k)*GV%Z_to_H
  enddo ; enddo ; enddo

  call sum_across_PEs(dilations)
  if ((dilations > 0) .and. (is_root_pe())) then
    write(mesg,'("Thickness initial conditions were dilated ",'// &
               '"to fit topography in ",I8," places.")') dilations
    call MOM_error(WARNING, 'adjustEtaToFitBathymetry: '//mesg)
  endif

end subroutine adjustEtaToFitBathymetry

!> Initializes thickness to be uniform
subroutine initialize_thickness_uniform(h, G, GV, param_file, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [H ~> m or kg m-2].
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
  character(len=40)  :: mdl = "initialize_thickness_uniform" ! This subroutine's name.
  real :: e0(SZK_(G)+1)   ! The resting interface heights, in depth units, usually
                          ! negative because it is positive upward.
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface
                          ! positive upward, in depth units.
  logical :: just_read    ! If true, just read parameters but set nothing.
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (just_read) return ! This subroutine has no run-time parameters.

  call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  if (G%max_depth<=0.) call MOM_error(FATAL,"initialize_thickness_uniform: "// &
      "MAXIMUM_DEPTH has a non-sensical value! Was it set?")

  do k=1,nz
    e0(K) = -G%max_depth * real(k-1) / real(nz)
  enddo

  do j=js,je ; do i=is,ie
    ! This sets the initial thickness (in m) of the layers.  The
    ! thicknesses are set to insure that: 1.  each layer is at least an
    ! Angstrom thick, and 2.  the interfaces are where they should be
    ! based on the resting depths and interface height perturbations,
    ! as long at this doesn't interfere with 1.
    eta1D(nz+1) = -G%bathyT(i,j)
    do k=nz,1,-1
      eta1D(K) = e0(K)
      if (eta1D(K) < (eta1D(K+1) + GV%Angstrom_Z)) then
        eta1D(K) = eta1D(K+1) + GV%Angstrom_Z
        h(i,j,k) = GV%Angstrom_H
      else
        h(i,j,k) = GV%Z_to_H * (eta1D(K) - eta1D(K+1))
      endif
    enddo
  enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_thickness_uniform

!> Initialize thickness from a 1D list
subroutine initialize_thickness_list(h, G, GV, US, param_file, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [H ~> m or kg m-2].
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
  character(len=40)  :: mdl = "initialize_thickness_list" ! This subroutine's name.
  real :: e0(SZK_(G)+1)   ! The resting interface heights, in depth units [Z ~> m],
                          ! usually negative because it is positive upward.
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface
                          ! positive upward, in depth units [Z ~> m].
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=200) :: filename, eta_file, inputdir ! Strings for file/path
  character(len=72)  :: eta_var
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  call get_param(param_file, mdl, "INTERFACE_IC_FILE", eta_file, &
                 "The file from which horizontal mean initial conditions "//&
                 "for interface depths can be read.", fail_if_missing=.true.)
  call get_param(param_file, mdl, "INTERFACE_IC_VAR", eta_var, &
                 "The variable name for horizontal mean initial conditions "//&
                 "for interface depths relative to mean sea level.", &
                 default="eta")

  if (just_read) return

  call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  call get_param(param_file,  mdl, "INPUTDIR", inputdir, default=".")
  filename = trim(slasher(inputdir))//trim(eta_file)
  call log_param(param_file, mdl, "INPUTDIR/INTERFACE_IC_FILE", filename)

  e0(:) = 0.0
  call MOM_read_data(filename, eta_var, e0(:), scale=US%m_to_Z)

  if ((abs(e0(1)) - 0.0) > 0.001) then
    ! This list probably starts with the interior interface, so shift it up.
    do k=nz+1,2,-1 ; e0(K) = e0(K-1) ; enddo
    e0(1) = 0.0
  endif

  if (e0(2) > e0(1)) then ! Switch to the convention for interface heights increasing upward.
    do k=1,nz ; e0(K) = -e0(K) ; enddo
  endif

  do j=js,je ; do i=is,ie
    ! This sets the initial thickness (in m) of the layers.  The
    ! thicknesses are set to insure that: 1.  each layer is at least an
    ! Angstrom thick, and 2.  the interfaces are where they should be
    ! based on the resting depths and interface height perturbations,
    ! as long at this doesn't interfere with 1.
    eta1D(nz+1) = -G%bathyT(i,j)
    do k=nz,1,-1
      eta1D(K) = e0(K)
      if (eta1D(K) < (eta1D(K+1) + GV%Angstrom_Z)) then
        eta1D(K) = eta1D(K+1) + GV%Angstrom_Z
        h(i,j,k) = GV%Angstrom_H
      else
        h(i,j,k) = GV%Z_to_H * (eta1D(K) - eta1D(K+1))
      endif
    enddo
  enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_thickness_list

!> Search density space for location of layers (not implemented!)
subroutine initialize_thickness_search
  call MOM_error(FATAL,"  MOM_state_initialization.F90, initialize_thickness_search: NOT IMPLEMENTED")
end subroutine initialize_thickness_search

!> Converts thickness from geometric to pressure units
subroutine convert_thickness(h, G, GV, US, tv)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: h  !< Input geometric layer thicknesses being converted
                                               !! to layer pressure [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv !< A structure pointing to various
                                               !! thermodynamic variables
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    p_top, p_bot                  ! Pressure at the interfaces above and below a layer [R L2 T-2 ~> Pa]
  real :: dz_geo(SZI_(G),SZJ_(G)) ! The change in geopotential height across a layer [L2 T-2 ~> m2 s-2]
  real :: rho(SZI_(G))            ! The in situ density [R ~> kg m-3]
  real :: I_gEarth      ! Unit conversion factors divided by the gravitational acceleration
                        ! [H T2 R-1 L-2 ~> s2 m2 kg-1 or s2 m-1]
  real :: HR_to_pres    ! A conversion factor from the input geometric thicknesses times the layer
                        ! densities into pressure units [L2 T-2 H-1 ~> m s-2 or m4 kg-1 s-2].
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: itt, max_itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  max_itt = 10

  if (GV%Boussinesq) then
    call MOM_error(FATAL,"Not yet converting thickness with Boussinesq approx.")
  else
    I_gEarth = GV%RZ_to_H / GV%g_Earth
    HR_to_pres = GV%g_Earth * GV%H_to_Z

    if (associated(tv%eqn_of_state)) then
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        p_bot(i,j) = 0.0 ; p_top(i,j) = 0.0
      enddo ; enddo
      EOSdom(:) = EOS_domain(G%HI)
      do k=1,nz
        do j=js,je
          do i=is,ie ; p_top(i,j) = p_bot(i,j) ; enddo
          call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_top(:,j), rho, &
                                 tv%eqn_of_state, EOSdom)
          do i=is,ie
            p_bot(i,j) = p_top(i,j) + HR_to_pres * (h(i,j,k) * rho(i))
          enddo
        enddo

        do itt=1,max_itt
          call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p_top, p_bot, 0.0, G%HI, &
                                   tv%eqn_of_state, US, dz_geo)
          if (itt < max_itt) then ; do j=js,je
            call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_bot(:,j), rho, &
                                   tv%eqn_of_state, EOSdom)
            ! Use Newton's method to correct the bottom value.
            ! The hydrostatic equation is sufficiently linear that no bounds-checking is needed.
            do i=is,ie
              p_bot(i,j) = p_bot(i,j) + rho(i) * (HR_to_pres*h(i,j,k) - dz_geo(i,j))
            enddo
          enddo ; endif
        enddo

        do j=js,je ; do i=is,ie
          h(i,j,k) = (p_bot(i,j) - p_top(i,j)) * I_gEarth
        enddo ; enddo
      enddo
    else
      do k=1,nz ; do j=js,je ; do i=is,ie
        h(i,j,k) = h(i,j,k) * (GV%Rlay(k) / GV%Rho0)
      enddo ; enddo ; enddo
    endif
  endif

end subroutine convert_thickness

!> Depress the sea-surface based on an initial condition file
subroutine depress_surface(h, G, GV, US, param_file, tv, just_read_params)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure pointing to various thermodynamic variables
  logical,       optional, intent(in)    :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    eta_sfc  ! The free surface height that the model should use [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: &
    eta  ! The free surface height that the model should use [Z ~> m].
  real :: dilate  ! A ratio by which layers are dilated [nondim].
  real :: scale_factor ! A scaling factor for the eta_sfc values that are read
                       ! in, which can be used to change units, for example.
  character(len=40)  :: mdl = "depress_surface" ! This subroutine's name.
  character(len=200) :: inputdir, eta_srf_file ! Strings for file/path
  character(len=200) :: filename, eta_srf_var  ! Strings for file/path
  logical :: just_read    ! If true, just read parameters but set nothing.
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  ! Read the surface height (or pressure) from a file.

  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "SURFACE_HEIGHT_IC_FILE", eta_srf_file,&
                 "The initial condition file for the surface height.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "SURFACE_HEIGHT_IC_VAR", eta_srf_var, &
                 "The initial condition variable for the surface height.",&
                 default="SSH", do_not_log=just_read)
  filename = trim(inputdir)//trim(eta_srf_file)
  if (.not.just_read) &
    call log_param(param_file,  mdl, "INPUTDIR/SURFACE_HEIGHT_IC_FILE", filename)

  call get_param(param_file, mdl, "SURFACE_HEIGHT_IC_SCALE", scale_factor, &
                 "A scaling factor to convert SURFACE_HEIGHT_IC_VAR into units of m", &
                 units="variable", default=1.0, scale=US%m_to_Z, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  call MOM_read_data(filename, eta_srf_var, eta_sfc, G%Domain, scale=scale_factor)

  ! Convert thicknesses to interface heights.
  call find_eta(h, tv, G, GV, US, eta)

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
          h(i,j,k) = GV%Angstrom_H
        else
          h(i,j,k) = max(GV%Angstrom_H, h(i,j,k) * &
              (eta_sfc(i,j) - eta(i,j,K+1)) / (eta(i,j,K) - eta(i,j,K+1)) )
        endif
      enddo
    endif
  endif ; enddo ; enddo

end subroutine depress_surface

!> Adjust the layer thicknesses by cutting away the top of each model column at the depth
!! where the hydrostatic pressure matches an imposed surface pressure read from file.
subroutine trim_for_ice(PF, G, GV, US, ALE_CSp, tv, h, just_read_params)
  type(param_file_type),   intent(in)    :: PF !< Parameter file structure
  type(ocean_grid_type),   intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(ALE_CS),            pointer       :: ALE_CSp !< ALE control structure
  type(thermo_var_ptrs),   intent(inout) :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  logical,       optional, intent(in)    :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
  character(len=200) :: mdl = "trim_for_ice"
  real, dimension(SZI_(G),SZJ_(G)) :: p_surf ! Imposed pressure on ocean at surface [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: S_t, S_b ! Top and bottom edge values for reconstructions
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: T_t, T_b ! of salinity [ppt] and temperature [degC] within each layer.
  character(len=200) :: inputdir, filename, p_surf_file, p_surf_var ! Strings for file/path
  real :: scale_factor   ! A file-dependent scaling factor for the input pressure.
  real :: min_thickness  ! The minimum layer thickness, recast into Z units [Z ~> m].
  integer :: i, j, k
  logical :: default_2018_answers, remap_answers_2018
  logical :: just_read    ! If true, just read parameters but set nothing.
  logical :: use_remapping ! If true, remap the initial conditions.
  type(remapping_CS), pointer :: remap_CS => NULL()

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  call get_param(PF, mdl, "SURFACE_PRESSURE_FILE", p_surf_file, &
                 "The initial condition file for the surface pressure exerted by ice.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(PF, mdl, "SURFACE_PRESSURE_VAR", p_surf_var, &
                 "The initial condition variable for the surface pressure exerted by ice.", &
                 units="Pa", default="", do_not_log=just_read)
  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".", do_not_log=.true.)
  filename = trim(slasher(inputdir))//trim(p_surf_file)
  if (.not.just_read) call log_param(PF,  mdl, "!INPUTDIR/SURFACE_HEIGHT_IC_FILE", filename)

  call get_param(PF, mdl, "SURFACE_PRESSURE_SCALE", scale_factor, &
                 "A scaling factor to convert SURFACE_PRESSURE_VAR from "//&
                 "file SURFACE_PRESSURE_FILE into a surface pressure.", &
                 units="file dependent", default=1., do_not_log=just_read)
  call get_param(PF, mdl, "MIN_THICKNESS", min_thickness, 'Minimum layer thickness', &
                 units='m', default=1.e-3, do_not_log=just_read, scale=US%m_to_Z)
  call get_param(PF, mdl, "TRIMMING_USES_REMAPPING", use_remapping, &
                 'When trimming the column, also remap T and S.', &
                 default=.false., do_not_log=just_read)
  remap_answers_2018 = .true.
  if (use_remapping) then
    call get_param(PF, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
    call get_param(PF, mdl, "REMAPPING_2018_ANSWERS", remap_answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
  endif

  if (just_read) return ! All run-time parameters have been read, so return.

  call MOM_read_data(filename, p_surf_var, p_surf, G%Domain, &
                     scale=scale_factor*US%kg_m3_to_R*US%m_s_to_L_T**2)

  if (use_remapping) then
    allocate(remap_CS)
    call initialize_remapping(remap_CS, 'PLM', boundary_extrapolation=.true.)
  endif

  ! Find edge values of T and S used in reconstructions
  if ( associated(ALE_CSp) ) then ! This should only be associated if we are in ALE mode
    call TS_PLM_edge_values(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, .true.)
  else
!    call MOM_error(FATAL, "trim_for_ice: Does not work without ALE mode")
    do k=1,G%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      T_t(i,j,k) = tv%T(i,j,k) ; T_b(i,j,k) = tv%T(i,j,k)
      S_t(i,j,k) = tv%S(i,j,k) ; S_b(i,j,k) = tv%S(i,j,k)
    enddo ; enddo ; enddo
  endif

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    call cut_off_column_top(GV%ke, tv, GV, US, GV%g_Earth, G%bathyT(i,j), &
               min_thickness, tv%T(i,j,:), T_t(i,j,:), T_b(i,j,:), &
               tv%S(i,j,:), S_t(i,j,:), S_b(i,j,:), p_surf(i,j), h(i,j,:), remap_CS, &
               z_tol=1.0e-5*US%m_to_Z, remap_answers_2018=remap_answers_2018)
  enddo ; enddo

end subroutine trim_for_ice


!> Adjust the layer thicknesses by removing the top of the water column above the
!! depth where the hydrostatic pressure matches p_surf
subroutine cut_off_column_top(nk, tv, GV, US, G_earth, depth, min_thickness, T, T_t, T_b, &
                              S, S_t, S_b, p_surf, h, remap_CS, z_tol, remap_answers_2018)
  integer,               intent(in)    :: nk  !< Number of layers
  type(thermo_var_ptrs), intent(in)    :: tv  !< Thermodynamics structure
  type(verticalGrid_type), intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US  !< A dimensional unit scaling type
  real,                  intent(in)    :: G_earth !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real,                  intent(in)    :: depth !< Depth of ocean column [Z ~> m].
  real,                  intent(in)    :: min_thickness !< Smallest thickness allowed [Z ~> m].
  real, dimension(nk),   intent(inout) :: T   !< Layer mean temperature [degC]
  real, dimension(nk),   intent(in)    :: T_t !< Temperature at top of layer [degC]
  real, dimension(nk),   intent(in)    :: T_b !< Temperature at bottom of layer [degC]
  real, dimension(nk),   intent(inout) :: S   !< Layer mean salinity [ppt]
  real, dimension(nk),   intent(in)    :: S_t !< Salinity at top of layer [ppt]
  real, dimension(nk),   intent(in)    :: S_b !< Salinity at bottom of layer [ppt]
  real,                  intent(in)    :: p_surf !< Imposed pressure on ocean at surface [R L2 T-2 ~> Pa]
  real, dimension(nk),   intent(inout) :: h   !< Layer thickness [H ~> m or kg m-2]
  type(remapping_CS),    pointer       :: remap_CS !< Remapping structure for remapping T and S,
                                                   !! if associated
  real,        optional, intent(in)    :: z_tol !< The tolerance with which to find the depth
                                                !! matching the specified pressure [Z ~> m].
  logical,     optional, intent(in)    :: remap_answers_2018 !< If true, use the order of arithmetic
                                                !! and expressions that recover the answers for remapping
                                                !! from the end of 2018. Otherwise, use more robust
                                                !! forms of the same expressions.

  ! Local variables
  real, dimension(nk+1) :: e ! Top and bottom edge values for reconstructions [Z ~> m]
  real, dimension(nk) :: h0, S0, T0, h1, S1, T1
  real :: P_t, P_b  ! Top and bottom pressures [R L2 T-2 ~> Pa]
  real :: z_out, e_top
  logical :: answers_2018
  integer :: k

  answers_2018 = .true. ; if (present(remap_answers_2018)) answers_2018 = remap_answers_2018

  ! Calculate original interface positions
  e(nk+1) = -depth
  do k=nk,1,-1
    e(K) = e(K+1) + GV%H_to_Z*h(k)
    h0(k) = h(nk+1-k) ! Keep a copy to use in remapping
  enddo

  P_t = 0.
  e_top = e(1)
  do k=1,nk
    call find_depth_of_pressure_in_cell(T_t(k), T_b(k), S_t(k), S_b(k), e(K), e(K+1), &
                                        P_t, p_surf, GV%Rho0, G_earth, tv%eqn_of_state, &
                                        US, P_b, z_out, z_tol=z_tol)
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
      if (e(K) > e_top) then
        ! Original e(K) is too high
        e(K) = e_top
        e_top = e_top - min_thickness ! Next interface must be at least this deep
      endif
      ! This layer needs trimming
      h(k) = GV%Z_to_H * max( min_thickness, e(K) - e(K+1) )
      if (e(K) < e_top) exit ! No need to go further
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
    if (answers_2018) then
      call remapping_core_h(remap_CS, nk, h0, T0, nk, h1, T1, 1.0e-30*GV%m_to_H, 1.0e-10*GV%m_to_H)
      call remapping_core_h(remap_CS, nk, h0, S0, nk, h1, S1, 1.0e-30*GV%m_to_H, 1.0e-10*GV%m_to_H)
    else
      call remapping_core_h(remap_CS, nk, h0, T0, nk, h1, T1, GV%H_subroundoff, GV%H_subroundoff)
      call remapping_core_h(remap_CS, nk, h0, S0, nk, h1, S1, GV%H_subroundoff, GV%H_subroundoff)
    endif
    do k=1,nk
      S(k) = S1(nk+1-k)
      T(k) = T1(nk+1-k)
    enddo
  endif

end subroutine cut_off_column_top

!> Initialize horizontal velocity components from file
subroutine initialize_velocity_from_file(u, v, G, US, param_file, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(out) :: u  !< The zonal velocity that is being initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(out) :: v  !< The meridional velocity that is being initialized [L T-1 ~> m s-1]
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file to
                                                      !! parse for modelparameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
  character(len=40)  :: mdl = "initialize_velocity_from_file" ! This subroutine's name.
  character(len=200) :: filename,velocity_file,inputdir ! Strings for file/path
  logical :: just_read    ! If true, just read parameters but set nothing.

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mdl, "VELOCITY_FILE", velocity_file, &
                 "The name of the velocity initial condition file.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  if (just_read) return ! All run-time parameters have been read, so return.

  filename = trim(inputdir)//trim(velocity_file)
  call log_param(param_file, mdl, "INPUTDIR/VELOCITY_FILE", filename)

  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
         " initialize_velocity_from_file: Unable to open "//trim(filename))

  !  Read the velocities from a netcdf file.
  call MOM_read_vector(filename, "u", "v", u(:,:,:), v(:,:,:), G%Domain, scale=US%m_s_to_L_T)

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_velocity_from_file

!> Initialize horizontal velocity components to zero.
subroutine initialize_velocity_zero(u, v, G, param_file, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(out) :: u  !< The zonal velocity that is being initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(out) :: v  !< The meridional velocity that is being initialized [L T-1 ~> m s-1]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file to
                                                      !! parse for modelparameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
  character(len=200) :: mdl = "initialize_velocity_zero" ! This subroutine's name.
  logical :: just_read    ! If true, just read parameters but set nothing.
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  if (just_read) return ! All run-time parameters have been read, so return.

  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u(I,j,k) = 0.0
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v(i,J,k) = 0.0
  enddo ; enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_velocity_zero

!> Sets the initial velocity components to uniform
subroutine initialize_velocity_uniform(u, v, G, US, param_file, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(out) :: u  !< The zonal velocity that is being initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(out) :: v  !< The meridional velocity that is being initialized [L T-1 ~> m s-1]
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file to
                                                      !! parse for modelparameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  real    :: initial_u_const, initial_v_const
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=200) :: mdl = "initialize_velocity_uniform" ! This subroutine's name.
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  call get_param(param_file, mdl, "INITIAL_U_CONST", initial_u_const, &
                 "A initial uniform value for the zonal flow.", &
                 default=0.0, units="m s-1", scale=US%m_s_to_L_T, do_not_log=just_read)
  call get_param(param_file, mdl, "INITIAL_V_CONST", initial_v_const, &
                 "A initial uniform value for the meridional flow.", &
                 default=0.0, units="m s-1", scale=US%m_s_to_L_T, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    u(I,j,k) = initial_u_const
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    v(i,J,k) = initial_v_const
  enddo ; enddo ; enddo

end subroutine initialize_velocity_uniform

!> Sets the initial velocity components to be circular with
!! no flow at edges of domain and center.
subroutine initialize_velocity_circular(u, v, G, US, param_file, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(out) :: u  !< The zonal velocity that is being initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(out) :: v  !< The meridional velocity that is being initialized [L T-1 ~> m s-1]
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file to
                                                      !! parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
  character(len=200) :: mdl = "initialize_velocity_circular"
  real :: circular_max_u ! The amplitude of the zonal flow [L T-1 ~> m s-1]
  real :: dpi        ! A local variable storing pi = 3.14159265358979...
  real :: psi1, psi2 ! Values of the streamfunction at two points [L2 T-1 ~> m2 s-1]
  logical :: just_read    ! If true, just read parameters but set nothing.
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  call get_param(param_file, mdl, "CIRCULAR_MAX_U", circular_max_u, &
                 "The amplitude of zonal flow from which to scale the "// &
                 "circular stream function [m s-1].", &
                 units="m s-1", default=0., scale=US%m_s_to_L_T, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  dpi=acos(0.0)*2.0 ! pi

  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    psi1 = my_psi(I,j)
    psi2 = my_psi(I,j-1)
    u(I,j,k) = (psi1 - psi2) / G%dy_Cu(I,j) ! *(circular_max_u*G%len_lon/(2.0*dpi))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    psi1 = my_psi(i,J)
    psi2 = my_psi(i-1,J)
    v(i,J,k) = (psi2 - psi1) / G%dx_Cv(i,J) ! *(circular_max_u*G%len_lon/(2.0*dpi))
  enddo ; enddo ; enddo

  contains

  !> Returns the value of a circular stream function at (ig,jg) in [L2 T-1 ~> m2 s-1]
  real function my_psi(ig,jg)
    integer :: ig !< Global i-index
    integer :: jg !< Global j-index
    ! Local variables
    real :: x, y, r ! [nondim]

    x = 2.0*(G%geoLonBu(ig,jg)-G%west_lon) / G%len_lon - 1.0  ! -1<x<1
    y = 2.0*(G%geoLatBu(ig,jg)-G%south_lat) / G%len_lat - 1.0 ! -1<y<1
    r = sqrt( x**2 + y**2 ) ! Circular stream function is a function of radius only
    r = min(1.0, r) ! Flatten stream function in corners of box
    my_psi = 0.5*(1.0 - cos(dpi*r))
    my_psi = my_psi * (circular_max_u * G%US%m_to_L*G%len_lon*1e3 / dpi) ! len_lon is in km
  end function my_psi

end subroutine initialize_velocity_circular

!> Initializes temperature and salinity from file
subroutine initialize_temp_salt_from_file(T, S, G, param_file, just_read_params)
  type(ocean_grid_type),                  intent(in)  :: G    !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out) :: T !< The potential temperature that is
                                                             !! being initialized [degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out) :: S !< The salinity that is
                                                             !! being initialized [ppt]
  type(param_file_type),                  intent(in)  :: param_file !< A structure to parse for run-time parameters
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                           !! only read parameters without changing h.
  ! Local variables
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=200) :: filename, salt_filename ! Full paths to input files
  character(len=200) :: ts_file, salt_file, inputdir ! Strings for file/path
  character(len=40)  :: mdl = "initialize_temp_salt_from_file"
  character(len=64)  :: temp_var, salt_var ! Temperature and salinity names in files

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mdl, "TS_FILE", ts_file, &
                 "The initial condition file for temperature.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  filename = trim(inputdir)//trim(ts_file)
  if (.not.just_read) call log_param(param_file, mdl, "INPUTDIR/TS_FILE", filename)
  call get_param(param_file, mdl, "TEMP_IC_VAR", temp_var, &
                 "The initial condition variable for potential temperature.", &
                 default="PTEMP", do_not_log=just_read)
  call get_param(param_file, mdl, "SALT_IC_VAR", salt_var, &
                 "The initial condition variable for salinity.", &
                 default="SALT", do_not_log=just_read)
  call get_param(param_file, mdl, "SALT_FILE", salt_file, &
                 "The initial condition file for salinity.", &
                 default=trim(ts_file), do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
     " initialize_temp_salt_from_file: Unable to open "//trim(filename))

  ! Read the temperatures and salinities from netcdf files.
  call MOM_read_data(filename, temp_var, T(:,:,:), G%Domain)

  salt_filename = trim(inputdir)//trim(salt_file)
  if (.not.file_exists(salt_filename, G%Domain)) call MOM_error(FATAL, &
     " initialize_temp_salt_from_file: Unable to open "//trim(salt_filename))

  call MOM_read_data(salt_filename, salt_var, S(:,:,:), G%Domain)

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_temp_salt_from_file

!> Initializes temperature and salinity from a 1D profile
subroutine initialize_temp_salt_from_profile(T, S, G, param_file, just_read_params)
  type(ocean_grid_type),                    intent(in)  :: G !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out) :: T !< The potential temperature that is
                                                             !! being initialized [degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out) :: S !< The salinity that is
                                                             !! being initialized [ppt]
  type(param_file_type),                    intent(in)  :: param_file !< A structure to parse for run-time parameters
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                           !! only read parameters without changing h.
  ! Local variables
  real, dimension(SZK_(G)) :: T0, S0
  integer :: i, j, k
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=200) :: filename, ts_file, inputdir ! Strings for file/path
  character(len=40)  :: mdl = "initialize_temp_salt_from_profile"

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mdl, "TS_FILE", ts_file, &
                 "The file with the reference profiles for temperature "//&
                 "and salinity.", fail_if_missing=.not.just_read, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  filename = trim(inputdir)//trim(ts_file)
  call log_param(param_file, mdl, "INPUTDIR/TS_FILE", filename)
  if (.not.file_exists(filename)) call MOM_error(FATAL, &
     " initialize_temp_salt_from_profile: Unable to open "//trim(filename))

  ! Read the temperatures and salinities from a netcdf file.
  call MOM_read_data(filename, "PTEMP", T0(:))
  call MOM_read_data(filename, "SALT",  S0(:))

  do k=1,G%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    T(i,j,k) = T0(k) ; S(i,j,k) = S0(k)
  enddo ; enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_temp_salt_from_profile

!> Initializes temperature and salinity by fitting to density
subroutine initialize_temp_salt_fit(T, S, G, GV, US, param_file, eqn_of_state, P_Ref, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G            !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV           !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out) :: T !< The potential temperature that is
                                                       !! being initialized [degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out) :: S !< The salinity that is being
                                                       !! initialized [ppt].
  type(unit_scale_type),   intent(in)  :: US           !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file   !< A structure to parse for run-time
                                                       !! parameters.
  type(EOS_type),          pointer     :: eqn_of_state !< Integer that selects the equatio of state.
  real,                    intent(in)  :: P_Ref        !< The coordinate-density reference pressure
                                                       !! [R L2 T-2 ~> Pa].
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                       !! only read parameters without changing h.
  ! Local variables
  real :: T0(SZK_(G))   ! Layer potential temperatures [degC]
  real :: S0(SZK_(G))   ! Layer salinities [degC]
  real :: T_Ref         ! Reference Temperature [degC]
  real :: S_Ref         ! Reference Salinity [ppt]
  real :: pres(SZK_(G))      ! An array of the reference pressure [R L2 T-2 ~> Pa].
  real :: drho_dT(SZK_(G))   ! Derivative of density with temperature [R degC-1 ~> kg m-3 degC-1].
  real :: drho_dS(SZK_(G))   ! Derivative of density with salinity [R ppt-1 ~> kg m-3 ppt-1].
  real :: rho_guess(SZK_(G)) ! Potential density at T0 & S0 [R ~> kg m-3].
  logical :: fit_salin       ! If true, accept the prescribed temperature and fit the salinity.
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40)  :: mdl = "initialize_temp_salt_fit" ! This subroutine's name.
  integer :: i, j, k, itt, nz
  nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mdl, "T_REF", T_Ref, &
                 "A reference temperature used in initialization.", &
                 units="degC", fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "S_REF", S_Ref, &
                 "A reference salinity used in initialization.", units="PSU", &
                 default=35.0, do_not_log=just_read)
  call get_param(param_file, mdl, "FIT_SALINITY", fit_salin, &
                 "If true, accept the prescribed temperature and fit the "//&
                 "salinity; otherwise take salinity and fit temperature.", &
                 default=.false., do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  do k=1,nz
    pres(k) = P_Ref ; S0(k) = S_Ref
    T0(k) = T_Ref
  enddo

  call calculate_density(T0(1), S0(1), pres(1), rho_guess(1), eqn_of_state)
  call calculate_density_derivs(T0, S0, pres, drho_dT, drho_dS, eqn_of_state, (/1,1/) )

  if (fit_salin) then
    ! A first guess of the layers' temperatures.
    do k=nz,1,-1
      S0(k) = max(0.0, S0(1) + (GV%Rlay(k) - rho_guess(1)) / drho_dS(1))
    enddo
    ! Refine the guesses for each layer.
    do itt=1,6
      call calculate_density(T0, S0, pres, rho_guess, eqn_of_state)
      call calculate_density_derivs(T0, S0, pres, drho_dT, drho_dS, eqn_of_state)
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
      call calculate_density(T0, S0, pres, rho_guess, eqn_of_state)
      call calculate_density_derivs(T0, S0, pres, drho_dT, drho_dS, eqn_of_state)
      do k=1,nz
        T0(k) = T0(k) + (GV%Rlay(k) - rho_guess(k)) / drho_dT(k)
      enddo
    enddo
  endif

  do k=1,nz ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
    T(i,j,k) = T0(k) ; S(i,j,k) = S0(k)
  enddo ; enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_temp_salt_fit

!> Initializes T and S with linear profiles according to reference surface
!! layer salinity and temperature and a specified range.
!!
!! \remark Note that the linear distribution is set up with respect to the layer
!! number, not the physical position).
subroutine initialize_temp_salt_linear(T, S, G, param_file, just_read_params)
  type(ocean_grid_type),                    intent(in)  :: G          !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out) :: T !< The potential temperature that is
                                                             !! being initialized [degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out) :: S !< The salinity that is
                                                             !! being initialized [ppt]
  type(param_file_type),                    intent(in)  :: param_file !< A structure to parse for
                                                                     !! run-time parameters
  logical,                        optional, intent(in)  :: just_read_params !< If present and true,
                                                                      !! this call will only read
                                                                      !! parameters without
                                                                      !! changing h.

  integer :: k
  real  :: delta_S, delta_T
  real  :: S_top, T_top ! Reference salinity and temerature within surface layer
  real  :: S_range, T_range ! Range of salinities and temperatures over the vertical
  real  :: delta
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40)  :: mdl = "initialize_temp_salt_linear" ! This subroutine's name.

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")
  call get_param(param_file, mdl, "T_TOP", T_top, &
                 "Initial temperature of the top surface.", &
                 units="degC", fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "T_RANGE", T_range, &
                 "Initial temperature difference (top-bottom).", &
                 units="degC", fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "S_TOP", S_top, &
                 "Initial salinity of the top surface.", &
                 units="PSU", fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "S_RANGE", S_range, &
                 "Initial salinity difference (top-bottom).", &
                 units="PSU", fail_if_missing=.not.just_read, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  ! Prescribe salinity
! delta_S = S_range / ( G%ke - 1.0 )
! S(:,:,1) = S_top
! do k = 2,G%ke
!   S(:,:,k) = S(:,:,k-1) + delta_S
! enddo
  do k = 1,G%ke
    S(:,:,k) = S_top - S_range*((real(k)-0.5)/real(G%ke))
    T(:,:,k) = T_top - T_range*((real(k)-0.5)/real(G%ke))
  enddo

  ! Prescribe temperature
! delta_T = T_range / ( G%ke - 1.0 )
! T(:,:,1) = T_top
! do k = 2,G%ke
!   T(:,:,k) = T(:,:,k-1) + delta_T
! enddo
! delta = 1
! T(:,:,G%ke/2 - (delta-1):G%ke/2 + delta) = 1.0

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_temp_salt_linear

!> This subroutine sets the inverse restoration time (Idamp), and
!! the values towards which the interface heights and an arbitrary
!! number of tracers should be restored within each sponge. The
!! interface height is always subject to damping, and must always be
!! the first registered field.
subroutine initialize_sponges_file(G, GV, US, use_temperature, tv, param_file, CSp, ALE_CSp, Time)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  logical,                 intent(in) :: use_temperature !< If true, T & S are state variables.
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic
                                              !! variables.
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters.
  type(sponge_CS),         pointer    :: CSp  !< A pointer that is set to point to the control
                                              !! structure for this module (in layered mode).
  type(ALE_sponge_CS),     pointer    :: ALE_CSp  !< A pointer that is set to point to the control
                                                  !! structure for this module (in ALE mode).
  type(time_type),         intent(in) :: Time !< Time at the start of the run segment. Time_in
                                              !! overrides any value set for Time.
  ! Local variables
  real, allocatable, dimension(:,:,:) :: eta ! The target interface heights [Z ~> m].
  real, allocatable, dimension(:,:,:) :: h   ! The target interface thicknesses [H ~> m or kg m-2].

  real, dimension (SZI_(G),SZJ_(G),SZK_(G)) :: &
    tmp, tmp2 ! A temporary array for tracers.
  real, dimension (SZI_(G),SZJ_(G)) :: &
    tmp_2d ! A temporary array for tracers.

  real :: Idamp(SZI_(G),SZJ_(G))    ! The inverse damping rate [T-1 ~> s-1].
  real :: pres(SZI_(G))     ! An array of the reference pressure [R L2 T-2 ~> Pa].

  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, nz
  integer :: isd, ied, jsd, jed
  integer, dimension(4) :: siz
  integer :: nz_data  ! The size of the sponge source grid
  character(len=40) :: potemp_var, salin_var, Idamp_var, eta_var
  character(len=40) :: mdl = "initialize_sponges_file"
  character(len=200) :: damping_file, state_file  ! Strings for filenames
  character(len=200) :: filename, inputdir ! Strings for file/path and path.

  logical :: use_ALE ! True if ALE is being used, False if in layered mode
  logical :: new_sponges ! True if using the newer sponges which do not
                         ! need to reside on the model horizontal grid.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  pres(:) = 0.0 ; tmp(:,:,:) = 0.0 ; Idamp(:,:) = 0.0

  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "SPONGE_DAMPING_FILE", damping_file, &
                 "The name of the file with the sponge damping rates.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "SPONGE_STATE_FILE", state_file, &
                 "The name of the file with the state to damp toward.", &
                 default=damping_file)
  call get_param(param_file, mdl, "SPONGE_PTEMP_VAR", potemp_var, &
                 "The name of the potential temperature variable in "//&
                 "SPONGE_STATE_FILE.", default="PTEMP")
  call get_param(param_file, mdl, "SPONGE_SALT_VAR", salin_var, &
                 "The name of the salinity variable in "//&
                 "SPONGE_STATE_FILE.", default="SALT")
  call get_param(param_file, mdl, "SPONGE_ETA_VAR", eta_var, &
                 "The name of the interface height variable in "//&
                 "SPONGE_STATE_FILE.", default="ETA")
  call get_param(param_file, mdl, "SPONGE_IDAMP_VAR", Idamp_var, &
                 "The name of the inverse damping rate variable in "//&
                 "SPONGE_DAMPING_FILE.", default="IDAMP")
  call get_param(param_file, mdl, "USE_REGRIDDING", use_ALE, do_not_log = .true.)

  call get_param(param_file, mdl, "NEW_SPONGES", new_sponges, &
                 "Set True if using the newer sponging code which "//&
                 "performs on-the-fly regridding in lat-lon-time.",&
                 "of sponge restoring data.", default=.false.)

!  if (use_ALE) then
!    call get_param(param_file, mdl, "SPONGE_RESTORE_ETA", restore_eta, &
!                 "If true, then restore the interface positions towards "//&
!                 "target values (in ALE mode)", default = .false.)
!  endif

  filename = trim(inputdir)//trim(damping_file)
  call log_param(param_file, mdl, "INPUTDIR/SPONGE_DAMPING_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) &
    call MOM_error(FATAL, " initialize_sponges: Unable to open "//trim(filename))

  if (new_sponges .and. .not. use_ALE) &
    call MOM_error(FATAL, " initialize_sponges: Newer sponges are currently unavailable in layered mode ")

  call MOM_read_data(filename, "Idamp", Idamp(:,:), G%Domain, scale=US%T_to_s)

  ! Now register all of the fields which are damped in the sponge.
  ! By default, momentum is advected vertically within the sponge, but
  ! momentum is typically not damped within the sponge.

  filename = trim(inputdir)//trim(state_file)
  call log_param(param_file, mdl, "INPUTDIR/SPONGE_STATE_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) &
    call MOM_error(FATAL, " initialize_sponges: Unable to open "//trim(filename))

  ! The first call to set_up_sponge_field is for the interface heights if in layered mode.!

  if (.not. use_ALE) then
    allocate(eta(isd:ied,jsd:jed,nz+1)); eta(:,:,:) = 0.0
    call MOM_read_data(filename, eta_var, eta(:,:,:), G%Domain, scale=US%m_to_Z)

    do j=js,je ; do i=is,ie
      eta(i,j,nz+1) = -G%bathyT(i,j)
    enddo ; enddo
    do k=nz,1,-1 ; do j=js,je ; do i=is,ie
      if (eta(i,j,K) < (eta(i,j,K+1) + GV%Angstrom_Z)) &
        eta(i,j,K) = eta(i,j,K+1) + GV%Angstrom_Z
    enddo ; enddo ; enddo
    ! Set the inverse damping rates so that the model will know where to
    ! apply the sponges, along with the interface heights.
    call initialize_sponge(Idamp, eta, G, param_file, CSp, GV)
    deallocate(eta)
  elseif (.not. new_sponges) then ! ALE mode

    call field_size(filename,eta_var,siz,no_domain=.true.)
    if (siz(1) /= G%ieg-G%isg+1 .or. siz(2) /= G%jeg-G%jsg+1) &
      call MOM_error(FATAL,"initialize_sponge_file: Array size mismatch for sponge data.")

!   ALE_CSp%time_dependent_target = .false.
!   if (siz(4) > 1) ALE_CSp%time_dependent_target = .true.
    nz_data = siz(3)-1
    allocate(eta(isd:ied,jsd:jed,nz_data+1))
    allocate(h(isd:ied,jsd:jed,nz_data))

    call MOM_read_data(filename, eta_var, eta(:,:,:), G%Domain, scale=US%m_to_Z)

    do j=js,je ; do i=is,ie
      eta(i,j,nz+1) = -G%bathyT(i,j)
    enddo ; enddo

    do k=nz,1,-1 ; do j=js,je ; do i=is,ie
      if (eta(i,j,K) < (eta(i,j,K+1) + GV%Angstrom_Z)) &
        eta(i,j,K) = eta(i,j,K+1) + GV%Angstrom_Z
    enddo ; enddo ; enddo
    do k=1,nz; do j=js,je ; do i=is,ie
      h(i,j,k) = GV%Z_to_H*(eta(i,j,k)-eta(i,j,k+1))
    enddo ; enddo ; enddo
    call initialize_ALE_sponge(Idamp, G, param_file, ALE_CSp, h, nz_data)
    deallocate(eta)
    deallocate(h)
  else
    ! Initialize sponges without supplying sponge grid
    call initialize_ALE_sponge(Idamp, G, param_file, ALE_CSp)
  endif

  ! Now register all of the tracer fields which are damped in the
  ! sponge. By default, momentum is advected vertically within the
  ! sponge, but momentum is typically not damped within the sponge.

  if ( GV%nkml>0 .and. .not. new_sponges) then
    ! This call to set_up_sponge_ML_density registers the target values of the
    ! mixed layer density, which is used in determining which layers can be
    ! inflated without causing static instabilities.
    do i=is-1,ie ; pres(i) = tv%P_Ref ; enddo
    EOSdom(:) = EOS_domain(G%HI)

    call MOM_read_data(filename, potemp_var, tmp(:,:,:), G%Domain)
    call MOM_read_data(filename, salin_var, tmp2(:,:,:), G%Domain)

    do j=js,je
      call calculate_density(tmp(:,j,1), tmp2(:,j,1), pres, tmp_2d(:,j), tv%eqn_of_state, EOSdom)
    enddo

    call set_up_sponge_ML_density(tmp_2d, G, CSp)
  endif

  ! The remaining calls to set_up_sponge_field can be in any order.
  if ( use_temperature .and. .not. new_sponges) then
    call MOM_read_data(filename, potemp_var, tmp(:,:,:), G%Domain)
    call set_up_sponge_field(tmp, tv%T, G, nz, CSp)
    call MOM_read_data(filename, salin_var, tmp(:,:,:), G%Domain)
    call set_up_sponge_field(tmp, tv%S, G, nz, CSp)
  elseif (use_temperature) then
    call set_up_ALE_sponge_field(filename, potemp_var, Time, G, GV, US, tv%T, ALE_CSp)
    call set_up_ALE_sponge_field(filename, salin_var, Time, G, GV, US, tv%S, ALE_CSp)
  endif

end subroutine initialize_sponges_file

!> This subroutine sets the 4 bottom depths at velocity points to be the
!! maximum of the adjacent depths.
subroutine set_velocity_depth_max(G)
  type(ocean_grid_type), intent(inout) :: G !< The ocean's grid structure
  ! Local variables
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

!> Subroutine to pre-compute global integrals of grid quantities for
!! later use in reporting diagnostics
subroutine compute_global_grid_integrals(G, US)
  type(ocean_grid_type), intent(inout) :: G !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US !< A dimensional unit scaling type
  ! Local variables
  real, dimension(G%isc:G%iec, G%jsc:G%jec) :: tmpForSumming
  real :: area_scale
  integer :: i,j

  area_scale = US%L_to_m**2
  tmpForSumming(:,:) = 0.
  G%areaT_global = 0.0 ; G%IareaT_global = 0.0
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    tmpForSumming(i,j) = area_scale*G%areaT(i,j) * G%mask2dT(i,j)
  enddo ; enddo
  G%areaT_global = reproducing_sum(tmpForSumming)
  G%IareaT_global = 1. / (G%areaT_global)
end subroutine compute_global_grid_integrals

!> This subroutine sets the 4 bottom depths at velocity points to be the
!! minimum of the adjacent depths.
subroutine set_velocity_depth_min(G)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  ! Local variables
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

!> This subroutine determines the isopycnal or other coordinate interfaces and
!! layer potential temperatures and salinities directly from a z-space file on
!! a latitude-longitude grid.
subroutine MOM_temp_salt_initialize_from_Z(h, tv, G, GV, US, PF, just_read_params)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(out)   :: h    !< Layer thicknesses being initialized [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv   !< A structure pointing to various thermodynamic
                                                 !! variables including temperature and salinity
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: PF   !< A structure indicating the open file
                                                 !! to parse for model parameter values.
  logical,       optional, intent(in)    :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.

  ! Local variables
  character(len=200) :: filename   !< The name of an input file containing temperature
                                   !! and salinity in z-space; also used for  ice shelf area.
  character(len=200) :: tfilename  !< The name of an input file containing only temperature
                                   !! in z-space.
  character(len=200) :: sfilename  !< The name of an input file containing only salinity
                                   !! in z-space.
  character(len=200) :: shelf_file !< The name of an input file used for  ice shelf area.
  character(len=200) :: inputdir   !! The directory where NetCDF input filesare.
  character(len=200) :: mesg, area_varname, ice_shelf_file

  type(EOS_type), pointer :: eos => NULL()
  type(thermo_var_ptrs) :: tv_loc   ! A temporary thermo_var container
  type(verticalGrid_type) :: GV_loc ! A temporary vertical grid structure
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_initialize_layers_from_Z" ! This module's name.

  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: is, ie, js, je, nz ! compute domain indices
  integer :: isc,iec,jsc,jec    ! global compute domain indices
  integer :: isg, ieg, jsg, jeg ! global extent
  integer :: isd, ied, jsd, jed ! data domain indices

  integer :: i, j, k, ks, np, ni, nj
  integer :: nkml     ! The number of layers in the mixed layer.

  integer :: kd, inconsistent
  integer :: nkd      ! number of levels to use for regridding input arrays
  real    :: eps_Z    ! A negligibly thin layer thickness [Z ~> m].
  real    :: eps_rho  ! A negligibly small density difference [R ~> kg m-3].
  real    :: PI_180   ! for conversion from degrees to radians

  real, dimension(:,:), pointer :: shelf_area => NULL()
  real    :: Hmix_default ! The default initial mixed layer depth [m].
  real    :: Hmix_depth   ! The mixed layer depth in the initial condition [Z ~> m].
  real    :: dilate       ! A dilation factor to match topography [nondim]
  real    :: missing_value_temp, missing_value_salt
  logical :: correct_thickness
  character(len=40) :: potemp_var, salin_var
  character(len=8)  :: laynum

  integer, parameter :: niter=10   ! number of iterations for t/s adjustment to layer density
  logical :: just_read    ! If true, just read parameters but set nothing.
  logical            :: adjust_temperature = .true.  ! fit t/s to target densities
  real, parameter    :: missing_value = -1.e20
  real, parameter    :: temp_land_fill = 0.0, salt_land_fill = 35.0
  logical :: reentrant_x, tripolar_n,dbg
  logical :: debug = .false.  ! manually set this to true for verbose output

  ! data arrays
  real, dimension(:), allocatable :: z_edges_in, z_in ! Interface heights [Z ~> m]
  real, dimension(:), allocatable :: Rb  ! Interface densities [R ~> kg m-3]
  real, dimension(:,:,:), allocatable, target :: temp_z, salt_z, mask_z
  real, dimension(:,:,:), allocatable :: rho_z ! Densities in Z-space [R ~> kg m-3]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: zi   ! Interface heights [Z ~> m].
  integer, dimension(SZI_(G),SZJ_(G))  :: nlevs
  real, dimension(SZI_(G))   :: press  ! Pressures [R L2 T-2 ~> Pa].

  ! Local variables for ALE remapping
  real, dimension(:), allocatable :: hTarget ! Target thicknesses [Z ~> m].
  real, dimension(:,:), allocatable :: area_shelf_h ! Shelf-covered area per grid cell [L2 ~> m2]
  real, dimension(:,:), allocatable, target :: frac_shelf_h ! Fractional shelf area per grid cell [nondim]
  real, dimension(:,:,:), allocatable, target :: tmpT1dIn, tmpS1dIn
  real, dimension(:,:,:), allocatable :: tmp_mask_in
  real, dimension(:,:,:), allocatable :: h1 ! Thicknesses [H ~> m or kg m-2].
  real, dimension(:,:,:), allocatable :: dz_interface ! Change in position of interface due to regridding
  real :: zTopOfCell, zBottomOfCell ! Heights in Z units [Z ~> m].
  type(regridding_CS) :: regridCS ! Regridding parameters and work arrays
  type(remapping_CS) :: remapCS ! Remapping parameters and work arrays

  logical :: homogenize, useALEremapping, remap_full_column, remap_general, remap_old_alg
  logical :: answers_2018, default_2018_answers, hor_regrid_answers_2018
  logical :: use_ice_shelf
  logical :: pre_gridded
  logical :: separate_mixed_layer  ! If true, handle the mixed layers differently.
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

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")
  if (.not.just_read) call log_version(PF, mdl, version, "")

  inputdir = "." ;  call get_param(PF, mdl, "INPUTDIR", inputdir)
  inputdir = slasher(inputdir)

  eos => tv%eqn_of_state

  reentrant_x = .false. ; call get_param(PF, mdl, "REENTRANT_X", reentrant_x, default=.true.)
  tripolar_n = .false. ;  call get_param(PF, mdl, "TRIPOLAR_N", tripolar_n, default=.false.)

  call get_param(PF, mdl, "TEMP_SALT_Z_INIT_FILE",filename, &
                 "The name of the z-space input file used to initialize "//&
                 "temperatures (T) and salinities (S). If T and S are not "//&
                 "in the same file, TEMP_Z_INIT_FILE and SALT_Z_INIT_FILE "//&
                 "must be set.",default="temp_salt_z.nc",do_not_log=just_read)
  call get_param(PF, mdl, "TEMP_Z_INIT_FILE",tfilename, &
                 "The name of the z-space input file used to initialize "//&
                 "temperatures, only.", default=trim(filename),do_not_log=just_read)
  call get_param(PF, mdl, "SALT_Z_INIT_FILE",sfilename, &
                 "The name of the z-space input file used to initialize "//&
                 "temperatures, only.", default=trim(filename),do_not_log=just_read)
  filename = trim(inputdir)//trim(filename)
  tfilename = trim(inputdir)//trim(tfilename)
  sfilename = trim(inputdir)//trim(sfilename)
  call get_param(PF, mdl, "Z_INIT_FILE_PTEMP_VAR", potemp_var, &
                 "The name of the potential temperature variable in "//&
                 "TEMP_Z_INIT_FILE.", default="ptemp",do_not_log=just_read)
  call get_param(PF, mdl, "Z_INIT_FILE_SALT_VAR", salin_var, &
                 "The name of the salinity variable in "//&
                 "SALT_Z_INIT_FILE.", default="salt",do_not_log=just_read)
  call get_param(PF, mdl, "Z_INIT_HOMOGENIZE", homogenize, &
                 "If True, then horizontally homogenize the interpolated "//&
                 "initial conditions.", default=.false., do_not_log=just_read)
  call get_param(PF, mdl, "Z_INIT_ALE_REMAPPING", useALEremapping, &
                 "If True, then remap straight to model coordinate from file.", &
                 default=.false., do_not_log=just_read)
  call get_param(PF, mdl, "Z_INIT_REMAPPING_SCHEME", remappingScheme, &
                 "The remapping scheme to use if using Z_INIT_ALE_REMAPPING "//&
                 "is True.", default="PPM_IH4", do_not_log=just_read)
  call get_param(PF, mdl, "Z_INIT_REMAP_GENERAL", remap_general, &
                 "If false, only initializes to z* coordinates. "//&
                 "If true, allows initialization directly to general coordinates.",&
                 default=.false., do_not_log=just_read)
  call get_param(PF, mdl, "Z_INIT_REMAP_FULL_COLUMN", remap_full_column, &
                 "If false, only reconstructs profiles for valid data points. "//&
                 "If true, inserts vanished layers below the valid data.", &
                 default=remap_general, do_not_log=just_read)
  call get_param(PF, mdl, "Z_INIT_REMAP_OLD_ALG", remap_old_alg, &
                 "If false, uses the preferred remapping algorithm for initialization. "//&
                 "If true, use an older, less robust algorithm for remapping.", &
                 default=.false., do_not_log=just_read)
  call get_param(PF, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  call get_param(PF, mdl, "TEMP_SALT_INIT_VERTICAL_REMAP_ONLY", pre_gridded, &
                 "If true, initial conditions are on the model horizontal grid. " //&
                 "Extrapolation over missing ocean values is done using an ICE-9 "//&
                 "procedure with vertical ALE remapping .", &
                 default=.false.)
  if (useALEremapping) then
    call get_param(PF, mdl, "REMAPPING_2018_ANSWERS", answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
  endif
  call get_param(PF, mdl, "HOR_REGRID_2018_ANSWERS", hor_regrid_answers_2018, &
                 "If true, use the order of arithmetic for horizonal regridding that recovers "//&
                 "the answers from the end of 2018.  Otherwise, use rotationally symmetric "//&
                 "forms of the same expressions.", default=default_2018_answers)
  call get_param(PF, mdl, "ICE_SHELF", use_ice_shelf, default=.false.)
  if (use_ice_shelf) then
    call get_param(PF, mdl, "ICE_THICKNESS_FILE", ice_shelf_file, &
                 "The file from which the ice bathymetry and area are read.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
    shelf_file = trim(inputdir)//trim(ice_shelf_file)
    if (.not.just_read) call log_param(PF, mdl, "INPUTDIR/THICKNESS_FILE", shelf_file)
    call get_param(PF, mdl, "ICE_AREA_VARNAME", area_varname, &
                 "The name of the area variable in ICE_THICKNESS_FILE.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  endif
  if (.not.useALEremapping) then
    call get_param(PF, mdl, "ADJUST_THICKNESS", correct_thickness, &
                 "If true, all mass below the bottom removed if the "//&
                 "topography is shallower than the thickness input file "//&
                 "would indicate.", default=.false., do_not_log=just_read)

    call get_param(PF, mdl, "FIT_TO_TARGET_DENSITY_IC", adjust_temperature, &
                 "If true, all the interior layers are adjusted to "//&
                 "their target densities using mostly temperature "//&
                 "This approach can be problematic, particularly in the "//&
                 "high latitudes.", default=.true., do_not_log=just_read)
    call get_param(PF, mdl, "Z_INIT_SEPARATE_MIXED_LAYER", separate_mixed_layer, &
                 "If true, distribute the topmost Z_INIT_HMIX_DEPTH of water over NKML layers, "//&
                 "and do not correct the density of the topmost NKML+NKBL layers.  Otherwise "//&
                 "all layers are initialized based on the depths of their target densities.", &
                 default=.false., do_not_log=just_read.or.(GV%nkml==0))
    if (GV%nkml == 0) separate_mixed_layer = .false.
    call get_param(PF, mdl, "MINIMUM_DEPTH", Hmix_default, default=0.0)
    call get_param(PF, mdl, "Z_INIT_HMIX_DEPTH", Hmix_depth, &
                 "The mixed layer depth in the initial conditions when Z_INIT_SEPARATE_MIXED_LAYER "//&
                 "is set to true.", default=Hmix_default, units="m", scale=US%m_to_Z, &
                 do_not_log=(just_read .or. .not.separate_mixed_layer))
    ! Reusing MINIMUM_DEPTH for the default mixed layer depth may be a strange choice, but
    ! it reproduces previous answers.
  endif
  if (just_read) then
    call cpu_clock_end(id_clock_routine)
    return ! All run-time parameters have been read, so return.
  endif

  eps_z = GV%Angstrom_Z
  eps_rho = 1.0e-10*US%kg_m3_to_R

  ! Read input grid coordinates for temperature and salinity field
  ! in z-coordinate dataset. The file is REQUIRED to contain the
  ! following:
  !
  ! dimension variables:
  !          lon (degrees_E), lat (degrees_N), depth(meters)
  ! variables:
  !          ptemp(lon,lat,depth) : degC, potential temperature
  !          salt (lon,lat,depth) : ppt, salinity
  !
  ! The first record will be read if there are multiple time levels.
  ! The observation grid MUST tile the model grid. If the model grid extends
  ! to the North/South Pole past the limits of the input data, they are extrapolated using the average
  ! value at the northernmost/southernmost latitude.

  call horiz_interp_and_extrap_tracer(tfilename, potemp_var, 1.0, 1, &
       G, temp_z, mask_z, z_in, z_edges_in, missing_value_temp, reentrant_x, &
       tripolar_n, homogenize, m_to_Z=US%m_to_Z, answers_2018=hor_regrid_answers_2018, ongrid=pre_gridded)

  call horiz_interp_and_extrap_tracer(sfilename, salin_var, 1.0, 1, &
       G, salt_z, mask_z, z_in, z_edges_in, missing_value_salt, reentrant_x, &
       tripolar_n, homogenize, m_to_Z=US%m_to_Z, answers_2018=hor_regrid_answers_2018, ongrid=pre_gridded)

  kd = size(z_in,1)

  ! Convert the sign convention of Z_edges_in.
  do k=1,size(Z_edges_in,1) ; Z_edges_in(k) = -Z_edges_in(k) ; enddo

  allocate(rho_z(isd:ied,jsd:jed,kd))
  allocate(area_shelf_h(isd:ied,jsd:jed))
  allocate(frac_shelf_h(isd:ied,jsd:jed))

  ! Convert T&S to Absolute Salinity and Conservative Temperature if using TEOS10 or NEMO
  call convert_temp_salt_for_TEOS10(temp_z, salt_z, G%HI, kd, mask_z, eos)

  press(:) = tv%P_Ref
  EOSdom(:) = EOS_domain(G%HI)
  do k=1,kd ; do j=js,je
    call calculate_density(temp_z(:,j,k), salt_z(:,j,k), press, rho_z(:,j,k), eos, EOSdom)
  enddo ; enddo

  call pass_var(temp_z,G%Domain)
  call pass_var(salt_z,G%Domain)
  call pass_var(mask_z,G%Domain)
  call pass_var(rho_z,G%Domain)

  ! This is needed for building an ALE grid under ice shelves
  if (use_ice_shelf) then
    if (.not.file_exists(shelf_file, G%Domain)) call MOM_error(FATAL, &
      "MOM_temp_salt_initialize_from_Z: Unable to open shelf file "//trim(shelf_file))

    call MOM_read_data(shelf_file, trim(area_varname), area_shelf_h, G%Domain, scale=US%m_to_L**2)

    ! Initialize frac_shelf_h with zeros (open water everywhere)
    frac_shelf_h(:,:) = 0.0
    ! Compute fractional ice shelf coverage of h
    do j=jsd,jed ; do i=isd,ied
      if (G%areaT(i,j) > 0.0) &
        frac_shelf_h(i,j) = area_shelf_h(i,j) / G%areaT(i,j)
    enddo ; enddo
    ! Pass to the pointer for use as an argument to regridding_main
    shelf_area => frac_shelf_h

  endif

  ! Done with horizontal interpolation.
  ! Now remap to model coordinates
  if (useALEremapping) then
    call cpu_clock_begin(id_clock_ALE)
    nkd = max(GV%ke, kd)

    ! Build the source grid and copy data onto model-shaped arrays with vanished layers
    allocate( tmp_mask_in(isd:ied,jsd:jed,nkd) ) ; tmp_mask_in(:,:,:) = 0.
    allocate( h1(isd:ied,jsd:jed,nkd) ) ; h1(:,:,:) = 0.
    allocate( tmpT1dIn(isd:ied,jsd:jed,nkd) ) ; tmpT1dIn(:,:,:) = 0.
    allocate( tmpS1dIn(isd:ied,jsd:jed,nkd) ) ; tmpS1dIn(:,:,:) = 0.
    do j = js, je ; do i = is, ie
      if (G%mask2dT(i,j)>0.) then
        zTopOfCell = 0. ; zBottomOfCell = 0.
        tmp_mask_in(i,j,1:kd) = mask_z(i,j,:)
        do k = 1, nkd
          if (tmp_mask_in(i,j,k)>0. .and. k<=kd) then
            zBottomOfCell = max( z_edges_in(k+1), -G%bathyT(i,j) )
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
          h1(i,j,k) = GV%Z_to_H * (zTopOfCell - zBottomOfCell)
          zTopOfCell = zBottomOfCell ! Bottom becomes top for next value of k
        enddo
        h1(i,j,kd) = h1(i,j,kd) + GV%Z_to_H * max(0., zTopOfCell + G%bathyT(i,j) )
        ! The max here is in case the data data is shallower than model
      endif ! mask2dT
    enddo ; enddo
    deallocate( tmp_mask_in )
    call pass_var(h1, G%Domain)
    call pass_var(tmpT1dIn, G%Domain)
    call pass_var(tmpS1dIn, G%Domain)

    ! Build the target grid (and set the model thickness to it)
    ! This call can be more general but is hard-coded for z* coordinates...  ????
    call ALE_initRegridding( GV, US, G%max_depth, PF, mdl, regridCS ) ! sets regridCS

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
            h(i,j,k) = GV%Z_to_H * (zTopOfCell - zBottomOfCell)
            zTopOfCell = zBottomOfCell ! Bottom becomes top for next value of k
          enddo
        else
          h(i,j,:) = 0.
        endif ! mask2dT
      enddo ; enddo
      call pass_var(h, G%Domain)
      deallocate( hTarget )
    endif

    ! Now remap from source grid to target grid, first setting reconstruction parameters
    call initialize_remapping( remapCS, remappingScheme, boundary_extrapolation=.false., answers_2018=answers_2018 )
    if (remap_general) then
      call set_regrid_params( regridCS, min_thickness=0. )
      tv_loc = tv
      tv_loc%T => tmpT1dIn
      tv_loc%S => tmpS1dIn
      GV_loc = GV
      GV_loc%ke = nkd
      allocate( dz_interface(isd:ied,jsd:jed,nkd+1) ) ! Need for argument to regridding_main() but is not used
      if (use_ice_shelf) then
        call regridding_main( remapCS, regridCS, G, GV_loc, h1, tv_loc, h, dz_interface, shelf_area )
      else
        call regridding_main( remapCS, regridCS, G, GV_loc, h1, tv_loc, h, dz_interface )
      endif
      deallocate( dz_interface )
    endif
    call ALE_remap_scalar(remapCS, G, GV, nkd, h1, tmpT1dIn, h, tv%T, all_cells=remap_full_column, &
                          old_remap=remap_old_alg, answers_2018=answers_2018 )
    call ALE_remap_scalar(remapCS, G, GV, nkd, h1, tmpS1dIn, h, tv%S, all_cells=remap_full_column, &
                          old_remap=remap_old_alg, answers_2018=answers_2018 )
    deallocate( h1 )
    deallocate( tmpT1dIn )
    deallocate( tmpS1dIn )

    call cpu_clock_end(id_clock_ALE)

  else ! remap to isopycnal layer space

    ! Next find interface positions using local arrays
    ! nlevs contains the number of valid data points in each column
    nlevs = int(sum(mask_z,dim=3))

    ! Rb contains the layer interface densities
    allocate(Rb(nz+1))
    do k=2,nz ; Rb(k) = 0.5*(GV%Rlay(k-1)+GV%Rlay(k)) ; enddo
    Rb(1) = 0.0 ;  Rb(nz+1) = 2.0*GV%Rlay(nz) - GV%Rlay(nz-1)

    nkml = 0 ; if (separate_mixed_layer) nkml = GV%nkml

    call find_interfaces(rho_z, z_in, kd, Rb, G%bathyT, zi, G, US, &
                         nlevs, nkml, hml=Hmix_depth, eps_z=eps_z, eps_rho=eps_rho)

    if (correct_thickness) then
      call adjustEtaToFitBathymetry(G, GV, US, zi, h)
    else
      do k=nz,1,-1 ; do j=js,je ; do i=is,ie
        if (zi(i,j,K) < (zi(i,j,K+1) + GV%Angstrom_Z)) then
          zi(i,j,K) = zi(i,j,K+1) + GV%Angstrom_Z
          h(i,j,k) = GV%Angstrom_H
        else
          h(i,j,k) = GV%Z_to_H * (zi(i,j,K) - zi(i,j,K+1))
        endif
      enddo ; enddo ; enddo
      inconsistent=0
      do j=js,je ; do i=is,ie
        if (abs(zi(i,j,nz+1) + G%bathyT(i,j)) > 1.0*US%m_to_Z) &
          inconsistent = inconsistent + 1
      enddo ; enddo
      call sum_across_PEs(inconsistent)

      if ((inconsistent > 0) .and. (is_root_pe())) then
        write(mesg, '("Thickness initial conditions are inconsistent ",'// &
                    '"with topography in ",I5," places.")') inconsistent
        call MOM_error(WARNING, mesg)
      endif
    endif

    call tracer_z_init_array(temp_z, z_edges_in, kd, zi, missing_value, G, nz, nlevs, eps_z, tv%T)
    call tracer_z_init_array(salt_z, z_edges_in, kd, zi, missing_value, G, nz, nlevs, eps_z, tv%S)

    do k=1,nz
      nPoints = 0 ; tempAvg = 0. ; saltAvg = 0.
      do j=js,je ; do i=is,ie ; if (G%mask2dT(i,j) >= 1.0) then
        nPoints = nPoints + 1
        tempAvg = tempAvg + tv%T(i,j,k)
        saltAvg = saltAvg + tv%S(i,j,k)
      endif ; enddo ; enddo

      ! Horizontally homogenize data to produce perfectly "flat" initial conditions
      if (homogenize) then
        call sum_across_PEs(nPoints)
        call sum_across_PEs(tempAvg)
        call sum_across_PEs(saltAvg)
        if (nPoints>0) then
          tempAvg = tempAvg / real(nPoints)
          saltAvg = saltAvg / real(nPoints)
        endif
        tv%T(:,:,k) = tempAvg
        tv%S(:,:,k) = saltAvg
      endif
    enddo

  endif ! useALEremapping

  ! Fill land values
  do k=1,nz ; do j=js,je ; do i=is,ie
    if (tv%T(i,j,k) == missing_value) then
      tv%T(i,j,k) = temp_land_fill
      tv%S(i,j,k) = salt_land_fill
    endif
  enddo ; enddo ; enddo


  if (adjust_temperature .and. .not. useALEremapping) then
    ! Finally adjust to target density
    ks = 1 ; if (separate_mixed_layer) ks = GV%nk_rho_varies + 1
    call determine_temperature(tv%T, tv%S, GV%Rlay(1:nz), tv%P_Ref, niter, &
                               missing_value, h, ks, G, US, eos)
  endif

  deallocate(z_in, z_edges_in, temp_z, salt_z, mask_z)
  deallocate(rho_z, area_shelf_h, frac_shelf_h)

  call pass_var(h, G%Domain)
  call pass_var(tv%T, G%Domain)
  call pass_var(tv%S, G%Domain)

  call callTree_leave(trim(mdl)//'()')
  call cpu_clock_end(id_clock_routine)

end subroutine MOM_temp_salt_initialize_from_Z

!> Run simple unit tests
subroutine MOM_state_init_tests(G, GV, US, tv)
  type(ocean_grid_type),     intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),     intent(in)    :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),     intent(in)    :: tv   !< Thermodynamics structure.

  ! Local variables
  integer, parameter :: nk=5
  real, dimension(nk) :: T, T_t, T_b ! Temperatures [degC]
  real, dimension(nk) :: S, S_t, S_b ! Salinities [ppt]
  real, dimension(nk) :: rho ! Layer density [R ~> kg m-3]
  real, dimension(nk) :: h   ! Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nk) :: z   ! Height of layer center [Z ~> m]
  real, dimension(nk+1) :: e ! Interface heights [Z ~> m]
  integer :: k
  real :: P_tot, P_t, P_b    ! Pressures [R L2 T-2 ~> Pa]
  real :: z_out              ! Output height [Z ~> m]
  real :: I_z_scale          ! The inverse of the height scale for prescribed gradients [Z-1 ~> m-1]
  type(remapping_CS), pointer :: remap_CS => NULL()

  I_z_scale = 1.0 / (500.0*US%m_to_Z)
  do k = 1, nk
    h(k) = 100.0*GV%m_to_H
  enddo
  e(1) = 0.
  do K = 1, nk
    e(K+1) = e(K) - GV%H_to_Z * h(k)
  enddo
  P_tot = 0.
  do k = 1, nk
    z(k) = 0.5 * ( e(K) + e(K+1) )
    T_t(k) = 20. + (0. * I_z_scale) * e(k)
    T(k)   = 20. + (0. * I_z_scale)*z(k)
    T_b(k) = 20. + (0. * I_z_scale)*e(k+1)
    S_t(k) = 35. - (0. * I_z_scale)*e(k)
    S(k)   = 35. + (0. * I_z_scale)*z(k)
    S_b(k) = 35. - (0. * I_z_scale)*e(k+1)
    call calculate_density(0.5*(T_t(k)+T_b(k)), 0.5*(S_t(k)+S_b(k)), -GV%Rho0*GV%g_Earth*US%m_to_Z*z(k), &
                           rho(k), tv%eqn_of_state)
    P_tot = P_tot + GV%g_Earth * rho(k) * GV%H_to_Z*h(k)
  enddo

  P_t = 0.
  do k = 1, nk
    call find_depth_of_pressure_in_cell(T_t(k), T_b(k), S_t(k), S_b(k), e(K), e(K+1), P_t, 0.5*P_tot, &
                                        GV%Rho0, GV%g_Earth, tv%eqn_of_state, US, P_b, z_out)
    write(0,*) k, US%RL2_T2_to_Pa*P_t, US%RL2_T2_to_Pa*P_b, 0.5*US%RL2_T2_to_Pa*P_tot, &
               US%Z_to_m*e(K), US%Z_to_m*e(K+1), US%Z_to_m*z_out
    P_t = P_b
  enddo
  write(0,*) US%RL2_T2_to_Pa*P_b, US%RL2_T2_to_Pa*P_tot

  write(0,*) ''
  write(0,*) ' ==================================================================== '
  write(0,*) ''
  write(0,*) GV%H_to_m*h
  call cut_off_column_top(nk, tv, GV, US, GV%g_Earth, -e(nk+1), GV%Angstrom_Z, &
                          T, T_t, T_b, S, S_t, S_b, 0.5*P_tot, h, remap_CS)
  write(0,*) GV%H_to_m*h

end subroutine MOM_state_init_tests

end module MOM_state_initialization
