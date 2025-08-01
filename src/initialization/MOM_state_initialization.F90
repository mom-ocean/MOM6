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
use MOM_interface_heights, only : find_eta, dz_to_thickness, dz_to_thickness_simple
use MOM_interface_heights, only : calc_derived_thermo
use MOM_io, only : file_exists, field_size, MOM_read_data, MOM_read_vector, slasher
use MOM_open_boundary, only : ocean_OBC_type, open_boundary_test_extern_h
use MOM_open_boundary, only : fill_temp_salt_segments, setup_OBC_tracer_reservoirs
use MOM_open_boundary, only : set_initialized_OBC_tracer_reservoirs
use MOM_grid_initialize, only : initialize_masks, set_grid_metrics
use MOM_restart, only : restore_state, is_new_run, copy_restart_var, copy_restart_vector
use MOM_restart, only : restart_registry_lock, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, set_up_sponge_ML_density
use MOM_sponge, only : initialize_sponge, sponge_CS
use MOM_ALE_sponge, only : set_up_ALE_sponge_field, set_up_ALE_sponge_vel_field
use MOM_ALE_sponge, only : ALE_sponge_CS, initialize_ALE_sponge
use MOM_string_functions, only : uppercase, lowercase
use MOM_time_manager, only : time_type, operator(/=)
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
use MOM_tracer_Z_init, only : tracer_Z_init_array, determine_temperature
use MOM_ALE, only : ALE_initRegridding, ALE_CS, ALE_initThicknessToCoord
use MOM_ALE, only : ALE_remap_scalar, ALE_regrid_accelerated, TS_PLM_edge_values
use MOM_regridding, only : regridding_CS, set_regrid_params, getCoordinateResolution
use MOM_regridding, only : regridding_main, regridding_preadjust_reqs, convective_adjustment
use MOM_regridding, only : set_dz_neglect, set_h_neglect
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h
use MOM_horizontal_regridding, only : horiz_interp_and_extrap_tracer, homogenize_field
use MOM_oda_incupd, only: oda_incupd_CS, initialize_oda_incupd_fixed, initialize_oda_incupd
use MOM_oda_incupd, only: set_up_oda_incupd_field, set_up_oda_incupd_vel_field
use MOM_oda_incupd, only: calc_oda_increments, output_oda_incupd_inc

implicit none ; private

#include <MOM_memory.h>

public MOM_initialize_state, MOM_initialize_OBCs

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
                                ALE_sponge_CSp, oda_incupd_CSp, OBC, Time_in, frac_shelf_h, mass_shelf)
  type(ocean_grid_type),      intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                              intent(out)   :: u    !< The zonal velocity that is being
                                                    !! initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                              intent(out)   :: v    !< The meridional velocity that is being
                                                    !! initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(out)   :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),      intent(inout) :: tv   !< A structure pointing to various thermodynamic
                                                    !! variables
  type(time_type),            intent(inout) :: Time !< Time at the start of the run segment.
  type(param_file_type),      intent(in)    :: PF   !< A structure indicating the open file to parse
                                                    !! for model parameter values.
  type(directories),          intent(in)    :: dirs !< A structure containing several relevant
                                                    !! directory paths.
  type(MOM_restart_CS),       intent(inout) :: restart_CS !< MOM restart control structure
  type(ALE_CS),               pointer       :: ALE_CSp !< The ALE control structure for remapping
  type(tracer_registry_type), pointer       :: tracer_Reg !< A pointer to the tracer registry
  type(sponge_CS),            pointer       :: sponge_CSp !< The layerwise sponge control structure.
  type(ALE_sponge_CS),        pointer       :: ALE_sponge_CSp !< The ALE sponge control structure.
  type(ocean_OBC_type),       pointer       :: OBC   !< The open boundary condition control structure.
                          ! OBC is only used in MOM_initialize_state if OBC_RESERVOIR_INIT_BUG is true.
  type(oda_incupd_CS),        pointer       :: oda_incupd_CSp !< The oda_incupd control structure.
  type(time_type), optional,  intent(in)    :: Time_in !< Time at the start of the run segment.
  real, dimension(SZI_(G),SZJ_(G)), &
                     optional, intent(in)   :: frac_shelf_h    !< The fraction of the grid cell covered
                                                               !! by a floating ice shelf [nondim].
  real, dimension(SZI_(G),SZJ_(G)), &
                     optional, intent(in)   :: mass_shelf      !< The mass per unit area of the overlying
                                                               !! ice shelf [ R Z ~> kg m-2 ]
  ! Local variables
  real :: depth_tot(SZI_(G),SZJ_(G))   ! The nominal total depth of the ocean [Z ~> m]
  real :: dz(SZI_(G),SZJ_(G),SZK_(GV)) ! The layer thicknesses in geopotential (z) units [Z ~> m]
  character(len=200) :: inputdir   ! The directory where NetCDF input files are.
  character(len=200) :: config, h_config
  real :: H_rescale   ! A rescaling factor for thicknesses from the representation in
                      ! a restart file to the internal representation in this run [various units ~> 1]
  real :: dt          ! The baroclinic dynamics timestep for this run [T ~> s].

  logical :: from_Z_file, useALE
  logical :: new_sim, rotate_index
  logical :: use_temperature, use_sponge, use_oda_incupd
  logical :: verify_restart_time
  logical :: OBC_reservoir_init_bug  ! If true, set the OBC tracer reservoirs at the startup of a new
                         ! run from the interior tracer concentrations regardless of properties that
                         ! may be explicitly specified for the reservoir concentrations.
  logical :: use_EOS     ! If true, density is calculated from T & S using an equation of state.
  logical :: depress_sfc ! If true, remove the mass that would be displaced
                         ! by a large surface pressure by squeezing the column.
  logical :: trim_ic_for_p_surf ! If true, remove the mass that would be displaced
                         ! by a large surface pressure, such as with an ice sheet.
  logical :: regrid_accelerate
  integer :: regrid_iterations
  logical :: convert
  logical :: just_read  ! If true, only read the parameters because this
                        ! is a run from a restart file; this option
                        ! allows the use of Fatal unused parameters.
  type(EOS_type), pointer :: eos => NULL()
  logical :: enable_bugs  ! If true, the defaults for recently added bug-fix flags are set to
                          ! recreate the bugs, or if false bugs are only used if actively selected.
  logical :: debug      ! If true, write debugging output.
  logical :: debug_layers = .false.
  logical :: use_ice_shelf
  character(len=80) :: mesg
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call callTree_enter("MOM_initialize_state(), MOM_state_initialization.F90")
  call log_version(PF, mdl, version, "")
  call get_param(PF, mdl, "DEBUG", debug, default=.false.)

  new_sim = is_new_run(restart_CS)
  just_read = .not.new_sim

  call get_param(PF, mdl, "INPUTDIR", inputdir, &
         "The directory in which input files are found.", default=".")
  inputdir = slasher(inputdir)

  use_temperature = associated(tv%T)
  useALE = associated(ALE_CSp)
  use_EOS = associated(tv%eqn_of_state)
  if (use_EOS) eos => tv%eqn_of_state
  use_ice_shelf = PRESENT(frac_shelf_h)

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
    !do k=1,nz ; do j=js,je ; do i=is,ie
    !  h(i,j,k) = 0.
    !enddo

    ! Initialize the layer thicknesses.
    dz(:,:,:) = 0.0
  endif

  ! Set the nominal depth of the ocean, which might be different from the bathymetric
  ! geopotential height, for use by the various initialization routines.  G%bathyT has
  ! already been initialized in previous calls.
  do j=jsd,jed ; do i=isd,ied
    depth_tot(i,j) = G%bathyT(i,j) + G%Z_ref
  enddo ; enddo

  call get_param(PF, mdl, "FATAL_INCONSISTENT_RESTART_TIME", verify_restart_time, &
                 "If true and a time_in value is provided to MOM_initialize_state, verify that "//&
                 "the time read from a restart file is the same as time_in, and issue a fatal "//&
                 "error if it is not.  Otherwise, simply set the time to time_in if present.", &
                 default=.false.)

  ! The remaining initialization calls are done, regardless of whether the
  ! fields are actually initialized here (if just_read=.false.) or whether it
  ! is just to make sure that all valid parameters are read to enable the
  ! detection of unused parameters.
  call get_param(PF, mdl, "INIT_LAYERS_FROM_Z_FILE", from_Z_file, &
             "If true, initialize the layer thicknesses, temperatures, and "//&
             "salinities from a Z-space file on a latitude-longitude grid.", &
             default=.false., do_not_log=just_read)

  convert = new_sim  ! Thicknesses are initialized in height units in most cases.
  if (from_Z_file) then
    ! Initialize thickness and T/S from z-coordinate data in a file.
    if (.NOT.use_temperature) call MOM_error(FATAL,"MOM_initialize_state : "//&
       "use_temperature must be true if INIT_LAYERS_FROM_Z_FILE is true")

    call MOM_temp_salt_initialize_from_Z(h, tv, depth_tot, G, GV, US, PF, &
                                         just_read=just_read, frac_shelf_h=frac_shelf_h)
    convert = .false.
  else
    ! Initialize thickness, h.
    call get_param(PF, mdl, "THICKNESS_CONFIG", h_config, &
             "A string that determines how the initial layer "//&
             "thicknesses are specified for a new run: \n"//&
             " \t file - read interface heights from the file specified \n"//&
             " \t\t by (THICKNESS_FILE).\n"//&
             " \t thickness_file - read thicknesses from the file specified \n"//&
             " \t\t by (THICKNESS_FILE).\n"//&
             " \t mass_file - read thicknesses in units of mass per unit area from the file \n"//&
             " \t\t specified by (THICKNESS_FILE).\n"//&
             " \t coord - determined by ALE coordinate.\n"//&
             " \t uniform - uniform thickness layers evenly distributed \n"//&
             " \t\t between the surface and MAXIMUM_DEPTH. \n"//&
             " \t list - read a list of positive interface depths. \n"//&
             " \t param - use thicknesses from parameter THICKNESS_INIT_VALUES. \n"//&
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
    select case (trim(h_config))
      case ("file")
        call initialize_thickness_from_file(dz, depth_tot, G, GV, US, PF, file_has_thickness=.false., &
                                            mass_file=.false., just_read=just_read)
      case ("thickness_file")
        call initialize_thickness_from_file(dz, depth_tot, G, GV, US, PF, file_has_thickness=.true., &
                                            mass_file=.false., just_read=just_read)
      case ("mass_file")
        call initialize_thickness_from_file(h, depth_tot, G, GV, US, PF, file_has_thickness=.true., &
                                            mass_file=.true., just_read=just_read)
        convert = .false.
      case ("coord")
        if (new_sim .and. useALE) then
          call ALE_initThicknessToCoord( ALE_CSp, G, GV, dz, height_units=.true. )
        elseif (new_sim) then
          call MOM_error(FATAL, "MOM_initialize_state: USE_REGRIDDING must be True "//&
                                "for THICKNESS_CONFIG of 'coord'")
        endif
      case ("uniform"); call initialize_thickness_uniform(dz, depth_tot, G, GV, PF, &
                                 just_read=just_read)
      case ("list"); call initialize_thickness_list(dz, depth_tot, G, GV, US, PF, &
                                 just_read=just_read)
      case ("param"); call initialize_thickness_param(dz, depth_tot, G, GV, US, PF, &
                                 just_read=just_read)
      case ("DOME"); call DOME_initialize_thickness(dz, depth_tot, G, GV, PF, &
                              just_read=just_read)
      case ("ISOMIP"); call ISOMIP_initialize_thickness(dz, depth_tot, G, GV, US, PF, tv, &
                                just_read=just_read)
      case ("benchmark"); call benchmark_initialize_thickness(dz, depth_tot, G, GV, US, PF, &
                                   tv%eqn_of_state, tv%P_Ref, just_read=just_read)
      case ("Neverworld","Neverland"); call Neverworld_initialize_thickness(dz, depth_tot, &
                                   G, GV, US, PF, tv%P_Ref)
      case ("search"); call initialize_thickness_search()
      case ("circle_obcs"); call circle_obcs_initialize_thickness(dz, depth_tot, G, GV, US, PF, &
                                     just_read=just_read)
      case ("lock_exchange"); call lock_exchange_initialize_thickness(dz, G, GV, US, &
                                       PF, just_read=just_read)
      case ("external_gwave"); call external_gwave_initialize_thickness(dz, G, GV, US, &
                                        PF, just_read=just_read)
      case ("DOME2D"); call DOME2d_initialize_thickness(dz, depth_tot, G, GV, US, PF, &
                                just_read=just_read)
      case ("adjustment2d"); call adjustment_initialize_thickness(dz, G, GV, US, &
                                      PF, just_read=just_read)
      case ("sloshing"); call sloshing_initialize_thickness(dz, depth_tot, G, GV, US, PF, &
                                  just_read=just_read)
      case ("seamount"); call seamount_initialize_thickness(dz, depth_tot, G, GV, US, PF, &
                                  just_read=just_read)
      case ("dumbbell"); call dumbbell_initialize_thickness(dz, depth_tot, G, GV, US, PF, &
                                  just_read=just_read)
      case ("soliton"); call soliton_initialize_thickness(dz, depth_tot, G, GV, US, PF, &
                                  just_read=just_read)
      case ("phillips"); call Phillips_initialize_thickness(dz, depth_tot, G, GV, US, PF, &
                                  just_read=just_read)
      case ("rossby_front")
        call Rossby_front_initialize_thickness(h, G, GV, US, PF, just_read=just_read)
        convert = .false.  ! Rossby_front initialization works directly in thickness units.
      case ("USER"); call user_initialize_thickness(dz, G, GV, PF, &
                              just_read=just_read)
      case default ; call MOM_error(FATAL,  "MOM_initialize_state: "//&
           "Unrecognized layer thickness configuration "//trim(h_config))
    end select

    ! Initialize temperature and salinity (T and S).
    if ( use_temperature ) then
      call get_param(PF, mdl, "TS_CONFIG", config, &
             "A string that determines how the initial temperatures "//&
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

      ! Check for incompatible THICKNESS_CONFIG and TS_CONFIG settings
      if (new_sim .and. (.not.convert)) then ; select case (trim(config))
        case ("DOME2D", "ISOMIP", "adjustment2d", "baroclinic_zone", "sloshing", &
              "seamount", "dumbbell", "SCM_CVMix_tests", "dense")
          call MOM_error(FATAL, "TS_CONFIG = "//trim(config)//" does not work with thicknesses "//&
              "that have already been converted to thickness units, as is the case with "//&
              "THICKNESS_CONFIG = "//trim(h_config)//".")
      end select ; endif

      select case (trim(config))
        case ("fit"); call initialize_temp_salt_fit(tv%T, tv%S, G, GV, US, PF, &
                               eos, tv%P_Ref, just_read=just_read)
        case ("file"); call initialize_temp_salt_from_file(tv%T, tv%S, G, GV, US, &
                                PF, just_read=just_read)
        case ("benchmark"); call benchmark_init_temperature_salinity(tv%T, tv%S, &
                                     G, GV, US, PF, eos, tv%P_Ref, just_read=just_read)
        case ("TS_profile") ; call initialize_temp_salt_from_profile(tv%T, tv%S, &
                                       G, GV, US, PF, just_read=just_read)
        case ("linear"); call initialize_temp_salt_linear(tv%T, tv%S, G, GV, US, PF, &
                                  just_read=just_read)
        case ("DOME2D"); call DOME2d_initialize_temperature_salinity (tv%T, tv%S, dz, &
                                  G, GV, US, PF, just_read=just_read)
        case ("ISOMIP"); call ISOMIP_initialize_temperature_salinity (tv%T, tv%S, dz, &
                                  depth_tot, G, GV, US, PF, eos, just_read=just_read)
        case ("adjustment2d"); call adjustment_initialize_temperature_salinity ( tv%T, &
                                        tv%S, dz, depth_tot, G, GV, US, PF, just_read=just_read)
        case ("baroclinic_zone"); call baroclinic_zone_init_temperature_salinity( tv%T, &
                                           tv%S, dz, depth_tot, G, GV, US, PF, just_read=just_read)
        case ("sloshing"); call sloshing_initialize_temperature_salinity(tv%T, &
                                    tv%S, dz, G, GV, US, PF, just_read=just_read)
        case ("seamount"); call seamount_initialize_temperature_salinity(tv%T, &
                                    tv%S, dz, G, GV, US, PF, just_read=just_read)
        case ("dumbbell"); call dumbbell_initialize_temperature_salinity(tv%T, &
                                    tv%S, dz, G, GV, US, PF, just_read=just_read)
        case ("rossby_front")
          if (convert .and. .not.just_read) call dz_to_thickness(dz, tv, h, G, GV, US)
          call Rossby_front_initialize_temperature_salinity ( tv%T, tv%S, h, &
                                        G, GV, US, PF, just_read=just_read)
        case ("SCM_CVMix_tests"); call SCM_CVMix_tests_TS_init(tv%T, tv%S, dz, &
                                           G, GV, US, PF, just_read=just_read)
        case ("dense"); call dense_water_initialize_TS(G, GV, US, PF, tv%T, tv%S, &
                                 dz, just_read=just_read)
        case ("USER"); call user_init_temperature_salinity(tv%T, tv%S, G, GV, PF, &
                                just_read=just_read)
        case default ; call MOM_error(FATAL,  "MOM_initialize_state: "//&
               "Unrecognized Temp & salt configuration "//trim(config))
      end select
    endif
  endif  ! not from_Z_file.

  if (use_temperature .and. associated(OBC)) then
    call get_param(PF, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
    ! Log this parameter later with the other OBC parameters.
    call get_param(PF, mdl, "OBC_RESERVOIR_INIT_BUG", OBC_reservoir_init_bug, &
                 "If true, set the OBC tracer reservoirs at the startup of a new run from the "//&
                 "interior tracer concentrations regardless of properties that may be explicitly "//&
                 "specified for the reservoir concentrations.", default=enable_bugs, do_not_log=.true.)
    if (OBC_reservoir_init_bug) then
      ! These calls should be moved down to join the OBC code, but doing so changes answers because
      ! the temperatures and salinities can change due to the remapping and reading from the restarts.
      call pass_var(tv%T, G%Domain, complete=.false.)
      call pass_var(tv%S, G%Domain, complete=.true.)
      call fill_temp_salt_segments(G, GV, US, OBC, tv)
    endif
  endif

  ! Convert thicknesses from geometric distances in depth units to thickness units or mass-per-unit-area.
  if (new_sim .and. convert) call dz_to_thickness(dz, tv, h, G, GV, US)

  ! Handle the initial surface displacement under ice shelf
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
    call hchksum(h, "Pre-depress: h ", G%HI, haloshift=1, unscale=GV%H_to_MKS)

  ! Remove the mass that would be displaced by an ice shelf or inverse barometer.
  if (depress_sfc) then
    call depress_surface(h, G, GV, US, PF, tv, just_read=just_read)
  elseif (trim_ic_for_p_surf) then
    call trim_for_ice(PF, G, GV, US, ALE_CSp, tv, h, just_read=just_read)
  elseif (new_sim .and. use_ice_shelf .and. present(mass_shelf)) then
    call calc_sfc_displacement(PF, G, GV, US, mass_shelf, tv, h)
  endif

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

      call get_param(PF, mdl, "DT", dt, "Timestep", &
                     units="s", scale=US%s_to_T, fail_if_missing=.true.)

      if (new_sim .and. debug) &
        call hchksum(h, "Pre-ALE_regrid: h ", G%HI, haloshift=1, unscale=GV%H_to_MKS)
      ! In this call, OBC is only used for the directions of OBCs when setting thicknesses at
      ! velocity points.
      call ALE_regrid_accelerated(ALE_CSp, G, GV, US, h, tv, regrid_iterations, u, v, OBC, tracer_Reg, &
                                  dt=dt, initial=.true.)
    endif
  endif

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
    case ("file"); call initialize_velocity_from_file(u, v, G, GV, US, PF, just_read)
    case ("zero"); call initialize_velocity_zero(u, v, G, GV, PF, just_read)
    case ("uniform"); call initialize_velocity_uniform(u, v, G, GV, US, PF, just_read)
    case ("circular"); call initialize_velocity_circular(u, v, G, GV, US, PF, just_read)
    case ("phillips"); call Phillips_initialize_velocity(u, v, G, GV, US, PF, just_read)
    case ("rossby_front"); call Rossby_front_initialize_velocity(u, v, h, &
                                     G, GV, US, PF, just_read)
    case ("soliton"); call soliton_initialize_velocity(u, v, G, GV, US, PF, just_read)
    case ("USER"); call user_initialize_velocity(u, v, G, GV, US, PF, just_read)
    case default ; call MOM_error(FATAL,  "MOM_initialize_state: "//&
          "Unrecognized velocity configuration "//trim(config))
  end select

  if (new_sim) call pass_vector(u, v, G%Domain)
  if (debug .and. new_sim) then
    call uvchksum("MOM_initialize_state [uv]", u, v, G%HI, haloshift=1, unscale=US%L_T_to_m_s)
  endif

  ! This is the end of the block of code that might have initialized fields
  ! internally at the start of a new run.

  ! Initialized assimilative incremental update (oda_incupd) structure and
  ! register restart.
  call get_param(PF, mdl, "ODA_INCUPD", use_oda_incupd, &
                 "If true, oda incremental updates will be applied "//&
                 "everywhere in the domain.", default=.false.)
  if (use_oda_incupd) then
    call restart_registry_lock(restart_CS, unlocked=.true.)
    call initialize_oda_incupd_fixed(G, GV, US, oda_incupd_CSp, restart_CS)
    call restart_registry_lock(restart_CS)
  endif

  if (.not.new_sim) then ! This block restores the state from a restart file.
    !    This line calls a subroutine that reads the initial conditions
    !  from a previously generated file.
    call restore_state(dirs%input_filename, dirs%restart_input_dir, Time, G, restart_CS)
    if (present(Time_in)) then
      if (verify_restart_time .and. (Time /= Time_in)) call MOM_error(FATAL, &
        "MOM6 attempted to restart from a file from a different time than given by Time_in.")
      Time = Time_in
    endif
    call get_param(PF, mdl, "ROTATE_INDEX", rotate_index, &
                 "Enable rotation of the horizontal indices.", &
                 default=.false., debuggingParam=.true., do_not_log=.true.)
    if (rotate_index) then
      ! This model is using a rotated grid, so the unrotated variables used here have not been set yet.
      call copy_restart_var(h, "h", restart_CS, .true.)
      call copy_restart_vector(u, v, "u", "v", restart_CS, .true.)
      if ( use_temperature ) then
        call copy_restart_var(tv%T, "Temp", restart_CS, .true.)
        call copy_restart_var(tv%S, "Salt", restart_CS, .true.)
      endif
    endif
  endif

  if ( use_temperature ) then
    call pass_var(tv%T, G%Domain, complete=.false.)
    call pass_var(tv%S, G%Domain, complete=.false.)
  endif
  call pass_var(h, G%Domain)

  if (debug) then
    call hchksum(h, "MOM_initialize_state: h ", G%HI, haloshift=1, unscale=GV%H_to_MKS)
    if ( use_temperature ) call hchksum(tv%T, "MOM_initialize_state: T ", G%HI, haloshift=1, unscale=US%C_to_degC)
    if ( use_temperature ) call hchksum(tv%S, "MOM_initialize_state: S ", G%HI, haloshift=1, unscale=US%S_to_ppt)
    if ( use_temperature .and. debug_layers) then ; do k=1,nz
      write(mesg,'("MOM_IS: T[",I2,"]")') k
      call hchksum(tv%T(:,:,k), mesg, G%HI, haloshift=1, unscale=US%C_to_degC)
      write(mesg,'("MOM_IS: S[",I2,"]")') k
      call hchksum(tv%S(:,:,k), mesg, G%HI, haloshift=1, unscale=US%S_to_ppt)
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
      case ("DOME"); call DOME_initialize_sponges(G, GV, US, tv, depth_tot, PF, sponge_CSp)
      case ("DOME2D"); call DOME2d_initialize_sponges(G, GV, US, tv, depth_tot, PF, useALE, &
                                                      sponge_CSp, ALE_sponge_CSp)
      case ("ISOMIP"); call ISOMIP_initialize_sponges(G, GV, US, tv, depth_tot, PF, useALE, &
                                                      sponge_CSp, ALE_sponge_CSp)
      case("RGC"); call RGC_initialize_sponges(G, GV, US, tv, u, v, depth_tot, PF, useALE, &
                                                     sponge_CSp, ALE_sponge_CSp)
      case ("USER"); call user_initialize_sponges(G, GV, use_temperature, tv, PF, sponge_CSp, h)
      case ("BFB"); call BFB_initialize_sponges_southonly(G, GV, US, use_temperature, tv, depth_tot, PF, &
                                                          sponge_CSp, h)
      case ("DUMBBELL"); call dumbbell_initialize_sponges(G, GV, US, tv, h, depth_tot, PF, useALE, &
                                                          sponge_CSp, ALE_sponge_CSp)
      case ("phillips"); call Phillips_initialize_sponges(G, GV, US, tv, PF, sponge_CSp, h)
      case ("dense"); call dense_water_initialize_sponges(G, GV, US, tv, depth_tot, PF, useALE, &
                                                          sponge_CSp, ALE_sponge_CSp)
      case ("file"); call initialize_sponges_file(G, GV, US, use_temperature, tv, u, v, depth_tot, PF, &
                                                  sponge_CSp, ALE_sponge_CSp, Time)
      case default ; call MOM_error(FATAL,  "MOM_initialize_state: "//&
             "Unrecognized sponge configuration "//trim(config))
    end select
  endif

  ! Set-up of data Assimilation with incremental update
  if (use_oda_incupd) then
    call initialize_oda_incupd_file(G, GV, US, use_temperature, tv, h, u, v, &
                                    PF, oda_incupd_CSp, restart_CS, Time)
  endif

  call callTree_leave('MOM_initialize_state()')

end subroutine MOM_initialize_state

subroutine MOM_initialize_OBCs(h, tv, OBC, Time, G, GV, US, PF, restart_CS, tracer_Reg)
  type(ocean_grid_type),      intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(out)   :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),      intent(inout) :: tv   !< A structure pointing to various thermodynamic
                                                    !! variables
  type(ocean_OBC_type),       pointer       :: OBC   !< The open boundary condition control structure.
  type(time_type),            intent(in)    :: Time !< Time at the start of the run segment.
  type(param_file_type),      intent(in)    :: PF   !< A structure indicating the open file to parse
                                                    !! for model parameter values.
  type(MOM_restart_CS),       intent(inout) :: restart_CS !< MOM restart control structure
  type(tracer_registry_type), pointer       :: tracer_Reg !< A pointer to the tracer registry

  ! Local variables
  character(len=200) :: config
  logical :: enable_bugs  ! If true, the defaults for recently added bug-fix flags are set to
                          ! recreate the bugs, or if false bugs are only used if actively selected.
  logical :: debug      ! If true, write debugging output.
  logical :: debug_obc  ! If true, do additional calls resetting values to help debug the correctness
                        ! of the open boundary condition code.
  logical :: OBC_reservoir_init_bug  ! If true, set the OBC tracer reservoirs at the startup of a new
                        ! run from the interior tracer concentrations regardless of properties that
                        ! may be explicitly specified for the reservoir concentrations.

  call callTree_enter('MOM_initialize_OBCs()')
  if (associated(OBC)) then
    call get_param(PF, mdl, "DEBUG", debug, default=.false.)
    call get_param(PF, mdl, "OBC_DEBUGGING_TESTS", debug_obc, &
                 "If true, do additional calls resetting values to help verify the correctness "//&
                 "of the open boundary condition code.", default=.false.,  &
                 do_not_log=.true., old_name="DEBUG_OBC", debuggingParam=.true.)
    call get_param(PF, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
    call get_param(PF, mdl, "OBC_RESERVOIR_INIT_BUG", OBC_reservoir_init_bug, &
                 "If true, set the OBC tracer reservoirs at the startup of a new run from the "//&
                 "interior tracer concentrations regardless of properties that may be explicitly "//&
                 "specified for the reservoir concentrations.", default=enable_bugs)
    if (associated(tv%T)) then
      if (OBC_reservoir_init_bug) then
        if (is_new_run(restart_CS)) then
          ! Set up OBC%trex_x and OBC%tres_y as they have not been read from a restart file.
          call setup_OBC_tracer_reservoirs(G, GV, OBC)
          ! Ensure that the values of the tracer reservoirs that have just been set will not be revised.
          call set_initialized_OBC_tracer_reservoirs(G, OBC, restart_CS)
        endif
      else
        ! Store the updated temperatures and salinities at the open boundaries, noting that they may
        ! still be updated by the calls in the next 50 lines, so the code setting the tracer
        ! reservoir values will come later in the calling routine.
        call fill_temp_salt_segments(G, GV, US, OBC, tv)
      endif
    endif

    ! This controls user code for setting open boundary data
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
      call supercritical_set_OBC_data(OBC, G, GV, US, PF)
    elseif (trim(config) == "tidal_bay") then
      OBC%update_OBC = .true.
    elseif (trim(config) == "USER") then
      call user_set_OBC_data(OBC, tv, G, GV, PF, tracer_Reg)
    elseif (.not. trim(config) == "none") then
      call MOM_error(FATAL, "The open boundary conditions specified by "//&
              "OBC_USER_CONFIG = "//trim(config)//" have not been fully implemented.")
    endif

    if (debug) then
      call hchksum(G%mask2dT, 'MOM_initialize_OBCs: mask2dT ', G%HI)
      call uvchksum('MOM_initialize_OBCs: mask2dC[uv]', G%mask2dCu, G%mask2dCv, G%HI)
      call qchksum(G%mask2dBu, 'MOM_initialize_OBCs: mask2dBu ', G%HI)
    endif
    if (debug_OBC) call open_boundary_test_extern_h(G, GV, OBC, h)
  endif

  call callTree_leave('MOM_initialize_OBCs()')

end subroutine MOM_initialize_OBCs

!> Reads the layer thicknesses or interface heights from a file.
subroutine initialize_thickness_from_file(h, depth_tot, G, GV, US, param_file, file_has_thickness, &
                                          just_read, mass_file)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h    !< The thickness that is being initialized, in height
                                               !! or thickness units, depending on the value of
                                               !! mass_file [Z ~> m] or [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file !< A structure indicating the open file
                                               !! to parse for model parameter values.
  logical,                 intent(in)  :: file_has_thickness !< If true, this file contains layer
                                               !! thicknesses; otherwise it contains
                                               !! interface heights.
  logical,                 intent(in)  :: just_read !< If true, this call will only read
                                               !! parameters without changing h.
  logical,                 intent(in)  :: mass_file !< If true, this file contains layer thicknesses in
                                               !! units of mass per unit area.

  ! Local variables
  real :: eta(SZI_(G),SZJ_(G),SZK_(GV)+1) ! Interface heights, in depth units [Z ~> m].
  real :: h_rescale   ! A factor by which to rescale the initial thickness variable in the input
                      ! file to convert it to units of m [various]
  real :: eta_rescale ! A factor by which to rescale the initial interface heights to convert
                      ! them to units of m or correct sign conventions to positive upward [various]
  real :: h_tolerance ! A parameter that controls the tolerance when adjusting the
                      ! thickness to fit the bathymetry [Z ~> m].
  real :: tol_dz_bot  ! A tolerance for detecting inconsistent bottom depths when
                      ! correct_thickness is false [Z ~> m]
  integer :: inconsistent ! The total number of cells with in consistent topography and layer thicknesses.
  logical :: correct_thickness
  character(len=40)  :: mdl = "initialize_thickness_from_file" ! This subroutine's name.
  character(len=200) :: filename, thickness_file, inputdir, mesg ! Strings for file/path
  character(len=80)  :: eta_var ! The interface height variable name in the input file
  character(len=80)  :: h_var   ! The thickness variable name in the input file
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.just_read) &
    call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".", do_not_log=just_read)
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "THICKNESS_FILE", thickness_file, &
                 "The name of the thickness file.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)

  filename = trim(thickness_file)
  if (scan(thickness_file, "/") == 0) then ! prepend inputdir if only a filename is given
    filename = trim(inputdir)//trim(thickness_file)
  endif
  if (.not.just_read) call log_param(param_file, mdl, "INPUTDIR/THICKNESS_FILE", filename)

  if ((.not.just_read) .and. (.not.file_exists(filename, G%Domain))) call MOM_error(FATAL, &
         " initialize_thickness_from_file: Unable to open "//trim(filename))

  if (file_has_thickness) then
    call get_param(param_file, mdl, "THICKNESS_IC_VAR", h_var, &
                 "The variable name for layer thickness initial conditions.", &
                 default="h", do_not_log=just_read)
    call get_param(param_file, mdl, "THICKNESS_IC_RESCALE", h_rescale, &
                 'A factor by which to rescale the initial thicknesses in the input file to '//&
                 'convert them to units of kg/m2 (if THICKNESS_CONFIG="mass_file") or m.', &
                 default=1.0, units="various", do_not_log=just_read)
    if (just_read) return ! All run-time parameters have been read, so return.

    if (mass_file) then
      h_rescale = h_rescale*GV%kg_m2_to_H
    else
      h_rescale = h_rescale*US%m_to_Z
    endif
    call MOM_read_data(filename, h_var, h(:,:,:), G%Domain, scale=h_rescale)
  else
    call get_param(param_file, mdl, "ADJUST_THICKNESS", correct_thickness, &
                 "If true, all mass below the bottom removed if the "//&
                 "topography is shallower than the thickness input file "//&
                 "would indicate.", default=.false., do_not_log=just_read)
    if (correct_thickness) then
      call get_param(param_file, mdl, "THICKNESS_TOLERANCE", h_tolerance, &
                 "A parameter that controls the tolerance when adjusting the "//&
                 "thickness to fit the bathymetry. Used when ADJUST_THICKNESS=True.", &
                 units="m", default=0.1, scale=US%m_to_Z, do_not_log=just_read)
    endif
    call get_param(param_file, mdl, "DZ_BOTTOM_TOLERANCE", tol_dz_bot, &
                 "A tolerance for detecting inconsistent topography and input layer "//&
                 "thicknesses when ADJUST_THICKNESS is false.", &
                 units="m", default=1.0, scale=US%m_to_Z, &
                 do_not_log=(just_read.or.correct_thickness))
    call get_param(param_file, mdl, "INTERFACE_IC_VAR", eta_var, &
                 "The variable name for initial conditions for interface heights "//&
                 "relative to mean sea level, positive upward unless otherwise rescaled.", &
                 default="eta", do_not_log=just_read)
    call get_param(param_file, mdl, "INTERFACE_IC_RESCALE", eta_rescale, &
                 "A factor by which to rescale the initial interface heights to convert "//&
                 "them to units of m or correct sign conventions to positive upward.", &
                 default=1.0, units="various", do_not_log=just_read)
    if (just_read) return ! All run-time parameters have been read, so return.

    call MOM_read_data(filename, eta_var, eta(:,:,:), G%Domain, scale=US%m_to_Z*eta_rescale)

    if (correct_thickness) then
      call adjustEtaToFitBathymetry(G, GV, US, eta, h, h_tolerance, dZ_ref_eta=G%Z_ref)
    else
      do k=nz,1,-1 ; do j=js,je ; do i=is,ie
        if (eta(i,j,K) < (eta(i,j,K+1) + GV%Angstrom_Z)) then
          eta(i,j,K) = eta(i,j,K+1) + GV%Angstrom_Z
          h(i,j,k) = GV%Angstrom_Z
        else
          h(i,j,k) = eta(i,j,K) - eta(i,j,K+1)
        endif
      enddo ; enddo ; enddo

      inconsistent = 0
      do j=js,je ; do i=is,ie
        if (abs(eta(i,j,nz+1) + depth_tot(i,j)) > tol_dz_bot) &
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
!! layers are contracted to ANGSTROM thickness (which may be 0).
!! If the bottom most interface is above the topography then the entire column
!! is dilated (expanded) to fill the void.
subroutine adjustEtaToFitBathymetry(G, GV, US, eta, h, ht, dZ_ref_eta)
  type(ocean_grid_type),                       intent(in)    :: G   !< The ocean's grid structure
  type(verticalGrid_type),                     intent(in)    :: GV  !< The ocean's vertical grid structure
  type(unit_scale_type),                       intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: eta !< Interface heights [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(inout) :: h   !< Layer thicknesses [Z ~> m]
  real,                                        intent(in)    :: ht  !< Tolerance to exceed adjustment
                                                                    !! criteria [Z ~> m]
  real,                              optional, intent(in)    :: dZ_ref_eta !< The difference between the
                                                                    !! reference heights for bathyT and
                                                                    !! eta [Z ~> m], 0 by default.
  ! Local variables
  integer :: i, j, k, is, ie, js, je, nz, contractions, dilations
  real :: dilate ! A factor by which the column is dilated [nondim]
  real :: dZ_ref ! The difference in the reference heights for G%bathyT and eta [Z ~> m]
  character(len=100) :: mesg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  dZ_ref = 0.0 ; if (present(dZ_ref_eta)) dZ_ref = dZ_ref_eta

  contractions = 0
  do j=js,je ; do i=is,ie
    if (-eta(i,j,nz+1) > (G%bathyT(i,j) + dZ_ref) + ht) then
      eta(i,j,nz+1) = -(G%bathyT(i,j) + dZ_ref)
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
    if (-eta(i,j,nz+1) < (G%bathyT(i,j) + dZ_ref) - ht) then
      dilations = dilations + 1
      if (eta(i,j,1) <= eta(i,j,nz+1)) then
        do k=1,nz ; h(i,j,k) = (eta(i,j,1) + (G%bathyT(i,j) + dZ_ref)) / real(nz) ; enddo
      else
        dilate = (eta(i,j,1) + (G%bathyT(i,j) + dZ_ref)) / (eta(i,j,1) - eta(i,j,nz+1))
        do k=1,nz ; h(i,j,k) = h(i,j,k) * dilate ; enddo
      endif
      do k=nz,2,-1 ; eta(i,j,K) = eta(i,j,K+1) + h(i,j,k) ; enddo
    endif
  enddo ; enddo


  call sum_across_PEs(dilations)
  if ((dilations > 0) .and. (is_root_pe())) then
    write(mesg,'("Thickness initial conditions were dilated ",'// &
               '"to fit topography in ",I8," places.")') dilations
    call MOM_error(WARNING, 'adjustEtaToFitBathymetry: '//mesg)
  endif

end subroutine adjustEtaToFitBathymetry

!> Initializes thickness to be uniform
subroutine initialize_thickness_uniform(h, depth_tot, G, GV, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.
  ! Local variables
  character(len=40)  :: mdl = "initialize_thickness_uniform" ! This subroutine's name.
  real :: e0(SZK_(GV)+1)  ! The resting interface heights [Z ~> m], usually
                          ! negative because it is positive upward.
  real :: eta1D(SZK_(GV)+1)! Interface height relative to the sea surface,
                          ! positive upward [Z ~> m].
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (just_read) return ! This subroutine has no run-time parameters.

  call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  if (G%max_depth<=0.) call MOM_error(FATAL,"initialize_thickness_uniform: "// &
      "MAXIMUM_DEPTH has a nonsensical value! Was it set?")

  do k=1,nz
    e0(K) = -G%max_depth * real(k-1) / real(nz)
  enddo

  do j=js,je ; do i=is,ie
    ! This sets the initial thickness (in m) of the layers.  The
    ! thicknesses are set to insure that: 1.  each layer is at least an
    ! Angstrom thick, and 2.  the interfaces are where they should be
    ! based on the resting depths and interface height perturbations,
    ! as long at this doesn't interfere with 1.
    eta1D(nz+1) = -depth_tot(i,j)
    do k=nz,1,-1
      eta1D(K) = e0(K)
      if (eta1D(K) < (eta1D(K+1) + GV%Angstrom_Z)) then
        eta1D(K) = eta1D(K+1) + GV%Angstrom_Z
        h(i,j,k) = GV%Angstrom_Z
      else
        h(i,j,k) = eta1D(K) - eta1D(K+1)
      endif
    enddo
  enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_thickness_uniform

!> Initialize thickness from a 1D list
subroutine initialize_thickness_list(h, depth_tot, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.
  ! Local variables
  character(len=40)  :: mdl = "initialize_thickness_list" ! This subroutine's name.
  real :: e0(SZK_(GV)+1)  ! The resting interface heights, in depth units [Z ~> m],
                          ! usually negative because it is positive upward.
  real :: eta1D(SZK_(GV)+1)! Interface height relative to the sea surface
                          ! positive upward, in depth units [Z ~> m].
  character(len=200) :: filename, eta_file, inputdir ! Strings for file/path
  character(len=72)  :: eta_var ! The interface height variable name in the input file
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

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
    eta1D(nz+1) = -depth_tot(i,j)
    do k=nz,1,-1
      eta1D(K) = e0(K)
      if (eta1D(K) < (eta1D(K+1) + GV%Angstrom_Z)) then
        eta1D(K) = eta1D(K+1) + GV%Angstrom_Z
        h(i,j,k) = GV%Angstrom_Z
      else
        h(i,j,k) = eta1D(K) - eta1D(K+1)
      endif
    enddo
  enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_thickness_list

!> Initializes thickness based on a run-time parameter with nominal thickness
!! for each layer
subroutine initialize_thickness_param(h, depth_tot, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.
  ! Local variables
  character(len=40)  :: mdl = "initialize_thickness_param" ! This subroutine's name.
  real :: e0(SZK_(GV)+1)  ! The resting interface heights [Z ~> m], usually
                          ! negative because it is positive upward.
  real :: eta1D(SZK_(GV)+1)! Interface height relative to the sea surface,
                          ! positive upward [Z ~> m].
  real :: dz(SZK_(GV))    ! The nominal initial layer thickness [Z ~> m], usually
  real :: h0_def(SZK_(GV)) ! Uniform default values for dz [Z ~> m], usually
  integer :: i, j, k, is, ie, js, je, nz

  call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")
  if (G%max_depth<=0.) call MOM_error(FATAL, "initialize_thickness_param: "// &
      "MAXIMUM_DEPTH has a nonsensical value! Was it set?")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  h0_def(:) = ( G%max_depth / real(nz) ) * US%Z_to_m
  call get_param(param_file, mdl, "THICKNESS_INIT_VALUES", dz, &
                 "A list of nominal thickness for each layer to initialize with", &
                 units="m", scale=US%m_to_Z, defaults=h0_def, do_not_log=just_read)
  if (just_read) return ! This subroutine has no run-time parameters.

  e0(nz+1) = -G%max_depth
  do k=nz, 1, -1
    e0(K) = e0(K+1) + dz(k)
  enddo

  do j=js,je ; do i=is,ie
    ! This sets the initial thickness (in m) of the layers.  The
    ! thicknesses are set to insure that: 1.  each layer is at least an
    ! Angstrom thick, and 2.  the interfaces are where they should be
    ! based on the resting depths and interface height perturbations,
    ! as long at this doesn't interfere with 1.
    eta1D(nz+1) = -depth_tot(i,j)
    do k=nz,1,-1
      eta1D(K) = e0(K)
      if (eta1D(K) < (eta1D(K+1) + GV%Angstrom_Z)) then
        eta1D(K) = eta1D(K+1) + GV%Angstrom_Z
        h(i,j,k) = GV%Angstrom_Z
      else
        h(i,j,k) = eta1D(K) - eta1D(K+1)
      endif
    enddo
  enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_thickness_param

!> Search density space for location of layers (not implemented!)
subroutine initialize_thickness_search
  call MOM_error(FATAL,"  MOM_state_initialization.F90, initialize_thickness_search: NOT IMPLEMENTED")
end subroutine initialize_thickness_search

!> Depress the sea-surface based on an initial condition file
subroutine depress_surface(h, G, GV, US, param_file, tv, just_read, z_top_shelf)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure pointing to various thermodynamic variables
  logical,                 intent(in)    :: just_read !< If true, this call will only read
                                                      !! parameters without changing h.
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: z_top_shelf    !< Top interface position under ice shelf [Z ~> m]
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    eta_sfc  ! The free surface height that the model should use [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    eta  ! The free surface height that the model should use [Z ~> m].
  real :: dilate  ! A ratio by which layers are dilated [nondim].
  real :: scale_factor ! A scaling factor for the eta_sfc values that are read in,
                       ! which can be used to change units, for example, often [Z m-1 ~> 1].
  character(len=40)  :: mdl = "depress_surface" ! This subroutine's name.
  character(len=200) :: inputdir, eta_srf_file ! Strings for file/path
  character(len=200) :: filename, eta_srf_var  ! Strings for file/path
  integer :: i, j, k, is, ie, js, je, nz
  logical :: use_z_shelf

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  use_z_shelf = present(z_top_shelf)


  if (.not. use_z_shelf) then
  ! Read the surface height (or pressure) from a file.

    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    call get_param(param_file, mdl, "SURFACE_HEIGHT_IC_FILE", eta_srf_file, &
                   "The initial condition file for the surface height.", &
                   fail_if_missing=.not.just_read, do_not_log=just_read)
    call get_param(param_file, mdl, "SURFACE_HEIGHT_IC_VAR", eta_srf_var, &
                   "The initial condition variable for the surface height.", &
                   default="SSH", do_not_log=just_read)
    filename = trim(inputdir)//trim(eta_srf_file)
    if (.not.just_read) &
      call log_param(param_file,  mdl, "INPUTDIR/SURFACE_HEIGHT_IC_FILE", filename)

    call get_param(param_file, mdl, "SURFACE_HEIGHT_IC_SCALE", scale_factor, &
                   "A scaling factor to convert SURFACE_HEIGHT_IC_VAR into units of m", &
                   units="variable", default=1.0, scale=US%m_to_Z, do_not_log=just_read)

    if (just_read) return ! All run-time parameters have been read, so return.

    call MOM_read_data(filename, eta_srf_var, eta_sfc, G%Domain, scale=scale_factor)
  else
    do j=js,je ; do i=is,ie
      eta_sfc(i,j) = z_top_shelf(i,j)
    enddo; enddo
  endif

  ! Convert thicknesses to interface heights.
  call find_eta(h, tv, G, GV, US, eta, dZref=G%Z_ref)

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
subroutine trim_for_ice(PF, G, GV, US, ALE_CSp, tv, h, just_read)
  type(param_file_type),   intent(in)    :: PF !< Parameter file structure
  type(ocean_grid_type),   intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(ALE_CS),            pointer       :: ALE_CSp !< ALE control structure
  type(thermo_var_ptrs),   intent(inout) :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  logical,                 intent(in)    :: just_read !< If true, this call will only read
                                                      !! parameters without changing h.
  ! Local variables
  character(len=200) :: mdl = "trim_for_ice"
  real, dimension(SZI_(G),SZJ_(G)) :: p_surf ! Imposed pressure on ocean at surface [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: S_t, S_b ! Top and bottom edge values for reconstructions
                                                        ! of salinity within each layer [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: T_t, T_b ! Top and bottom edge values for reconstructions
                                                        ! of temperature within each layer [C ~> degC]
  character(len=200) :: inputdir, filename, p_surf_file, p_surf_var ! Strings for file/path
  real :: scale_factor   ! A file-dependent scaling factor for the input pressure [various].
  real :: min_thickness  ! The minimum layer thickness [H ~> m or kg m-2].
  real :: z_tolerance    ! The tolerance with which to find the depth matching a specified pressure [Z ~> m].
  integer :: i, j, k
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  integer :: remap_answer_date    ! The vintage of the order of arithmetic and expressions to use
                                  ! for remapping.  Values below 20190101 recover the remapping
                                  ! answers from 2018, while higher values use more robust
                                  ! forms of the same remapping expressions.
  logical :: use_remapping ! If true, remap the initial conditions.
  logical :: use_frac_dp_bugfix   ! If true, use bugfix. Otherwise, pressure input to EOS is negative.
  type(remapping_CS), pointer :: remap_CS => NULL()

  call get_param(PF, mdl, "SURFACE_PRESSURE_FILE", p_surf_file, &
                 "The initial condition file for the surface pressure exerted by ice.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(PF, mdl, "SURFACE_PRESSURE_VAR", p_surf_var, &
                 "The initial condition variable for the surface pressure exerted by ice.", &
                 default="", do_not_log=just_read)
  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".", do_not_log=.true.)
  filename = trim(slasher(inputdir))//trim(p_surf_file)
  if (.not.just_read) call log_param(PF,  mdl, "!INPUTDIR/SURFACE_HEIGHT_IC_FILE", filename)

  call get_param(PF, mdl, "SURFACE_PRESSURE_SCALE", scale_factor, &
                 "A scaling factor to convert SURFACE_PRESSURE_VAR from "//&
                 "file SURFACE_PRESSURE_FILE into a surface pressure.", &
                 units="file dependent", default=1., do_not_log=just_read)
  call get_param(PF, mdl, "MIN_THICKNESS", min_thickness, 'Minimum layer thickness', &
                 units='m', default=1.e-3, scale=GV%m_to_H, do_not_log=just_read)
  call get_param(PF, mdl, "TRIM_IC_Z_TOLERANCE", z_tolerance, &
                 "The tolerance with which to find the depth matching the specified "//&
                 "surface pressure with TRIM_IC_FOR_P_SURF.", &
                 units="m", default=1.0e-5, scale=US%m_to_Z, do_not_log=just_read)
  call get_param(PF, mdl, "FRAC_DP_AT_POS_NEGATIVE_P_BUGFIX", use_frac_dp_bugfix, &
                 "If true, use bugfix in ice shelf TRIM_IC initialization. "//&
                 "Otherwise, pressure input to density EOS is negative.", &
                 default=.false., do_not_log=just_read)
  call get_param(PF, mdl, "TRIMMING_USES_REMAPPING", use_remapping, &
                 'When trimming the column, also remap T and S.', &
                 default=.false., do_not_log=just_read)
  if (use_remapping) then
    call get_param(PF, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231, do_not_log=just_read)
    call get_param(PF, mdl, "REMAPPING_ANSWER_DATE", remap_answer_date, &
                 "The vintage of the expressions and order of arithmetic to use for remapping.  "//&
                 "Values below 20190101 result in the use of older, less accurate expressions "//&
                 "that were in use at the end of 2018.  Higher values result in the use of more "//&
                 "robust and accurate forms of mathematically equivalent expressions.", &
                 default=default_answer_date, do_not_log=just_read.or.(.not.GV%Boussinesq))
    if (.not.GV%Boussinesq) remap_answer_date = max(remap_answer_date, 20230701)
  else
    remap_answer_date = 20181231
    if (.not.GV%Boussinesq) remap_answer_date = 20230701
  endif

  if (just_read) return ! All run-time parameters have been read, so return.

  call MOM_read_data(filename, p_surf_var, p_surf, G%Domain, &
                     scale=scale_factor*US%Pa_to_RL2_T2)

  if (use_remapping) then
    allocate(remap_CS)
    if (remap_answer_date < 20190101) then
      call initialize_remapping(remap_CS, 'PLM', boundary_extrapolation=.true., &
                                h_neglect=1.0e-30*GV%m_to_H, h_neglect_edge=1.0e-10*GV%m_to_H)
    else
      call initialize_remapping(remap_CS, 'PLM', boundary_extrapolation=.true., &
                                h_neglect=GV%H_subroundoff,  h_neglect_edge=GV%H_subroundoff)
    endif
  endif

  ! Find edge values of T and S used in reconstructions
  if ( associated(ALE_CSp) ) then ! This should only be associated if we are in ALE mode
    call TS_PLM_edge_values(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, .true.)
  else
!    call MOM_error(FATAL, "trim_for_ice: Does not work without ALE mode")
    do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      T_t(i,j,k) = tv%T(i,j,k) ; T_b(i,j,k) = tv%T(i,j,k)
      S_t(i,j,k) = tv%S(i,j,k) ; S_b(i,j,k) = tv%S(i,j,k)
    enddo ; enddo ; enddo
  endif

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    call cut_off_column_top(GV%ke, tv, GV, US, GV%g_Earth, G%bathyT(i,j)+G%Z_ref, min_thickness, &
               tv%T(i,j,:), T_t(i,j,:), T_b(i,j,:), tv%S(i,j,:), S_t(i,j,:), S_b(i,j,:), &
               p_surf(i,j), h(i,j,:), remap_CS, z_tol=z_tolerance, &
               frac_dp_bugfix=use_frac_dp_bugfix)
  enddo ; enddo

end subroutine trim_for_ice

!> Calculate the hydrostatic equilibrium position of the surface under an ice shelf
subroutine calc_sfc_displacement(PF, G, GV, US, mass_shelf, tv, h)
  type(param_file_type),   intent(in)    :: PF !< Parameter file structure
  type(ocean_grid_type),   intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: mass_shelf  !< Ice shelf mass [R Z ~> kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]

  real :: z_top_shelf(SZI_(G),SZJ_(G))  ! The depth of the top interface under ice shelves [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
                                   eta  ! The free surface height that the model should use [Z ~> m].
  ! temporary arrays
  real, dimension(SZK_(GV)) :: rho_col   ! potential density in the column for use in ice [R ~> kg m-3]
  real, dimension(SZK_(GV)) :: rho_h     ! potential density multiplied by thickness [R Z ~> kg m-2]
  real, dimension(SZK_(GV)) :: h_tmp     ! temporary storage for thicknesses [H ~> m]
  real, dimension(SZK_(GV)) :: p_ref     ! pressure for density [R Z ~> kg m-2]
  real, dimension(SZK_(GV)+1) :: ei_tmp, ei_orig ! temporary storage for interface positions [Z ~> m]
  real :: z_top     ! An estimate of the height of the ice-ocean interface [Z ~> m]
  real :: mass_disp ! The net mass of sea water that has been displaced by the shelf [R Z ~> kg m-2]
  real :: residual  ! The difference between the displaced ocean mass and the ice shelf
                    ! mass [R Z ~> kg m-2]
  real :: tol       ! The initialization tolerance for ice shelf initialization [Z ~> m]
  integer :: is, ie, js, je, k, nz, i, j, max_iter, iter

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call get_param(PF, mdl, "ICE_SHELF_INITIALIZATION_Z_TOLERANCE", tol, &
                "A initialization tolerance for the calculation of the static "// &
                "ice shelf displacement (m) using initial temperature and salinity profile.", &
                 default=0.001, units="m", scale=US%m_to_Z)
  max_iter = 1e3
  call MOM_mesg("Started calculating initial interface position under ice shelf ")
  ! Convert thicknesses to interface heights.
  call find_eta(h, tv, G, GV, US, eta, dZref=G%Z_ref)
  do j=js,je ; do i=is,ie
    iter = 1
    z_top_shelf(i,j) = 0.0
    p_ref(:) = tv%p_ref
    if ((G%mask2dT(i,j) > 0.) .and. (mass_shelf(i,j) > 0.)) then
      call calculate_density(tv%T(i,j,:), tv%S(i,j,:), P_Ref, rho_col, tv%eqn_of_state)
      z_top = min(max(-1.0*mass_shelf(i,j)/rho_col(1), -G%bathyT(i,j)), 0.)
      h_tmp(:) = 0.0
      ei_tmp(1:nz+1) = eta(i,j,1:nz+1)
      ei_orig(1:nz+1) = eta(i,j,1:nz+1)
      do k=1,nz+1
        if (ei_tmp(k) < z_top) ei_tmp(k) = z_top
      enddo
      mass_disp = 0.0
      do k=1,nz
        h_tmp(k) = max(ei_tmp(k)-ei_tmp(k+1), GV%Angstrom_H)
        rho_h(k) = h_tmp(k) * rho_col(k)
        mass_disp = mass_disp + rho_h(k)
      enddo
      residual = mass_shelf(i,j) - mass_disp
      do while ((abs(residual) > tol) .and. (z_top > -G%bathyT(i,j)) .and. (iter < max_iter))
        z_top = min(max(z_top-(residual*0.5e-3), -G%bathyT(i,j)), 0.0)
        h_tmp(:) = 0.0
        ei_tmp(1:nz+1) = ei_orig(1:nz+1)
        do k=1,nz+1
          if (ei_tmp(k) < z_top) ei_tmp(k) = z_top
        enddo
        mass_disp = 0.0
        do k=1,nz
          h_tmp(k) = max(ei_tmp(k)-ei_tmp(k+1), GV%Angstrom_H)
          rho_h(k) = h_tmp(k) * rho_col(k)
          mass_disp = mass_disp + rho_h(k)
        enddo
        residual = mass_shelf(i,j) - mass_disp
        iter = iter+1
      end do
      if (iter >= max_iter) call MOM_mesg("Warning: calc_sfc_displacement too many iterations.")
      z_top_shelf(i,j) = z_top
    endif
  enddo ; enddo
  call MOM_mesg("Calling depress_surface ")
  call depress_surface(h, G, GV, US, PF, tv, just_read=.false.,z_top_shelf=z_top_shelf)
  call MOM_mesg("Finishing calling depress_surface ")
end subroutine calc_sfc_displacement

!> Adjust the layer thicknesses by removing the top of the water column above the
!! depth where the hydrostatic pressure matches p_surf
subroutine cut_off_column_top(nk, tv, GV, US, G_earth, depth, min_thickness, T, T_t, T_b, &
                              S, S_t, S_b, p_surf, h, remap_CS, z_tol, frac_dp_bugfix)
  integer,               intent(in)    :: nk  !< Number of layers
  type(thermo_var_ptrs), intent(in)    :: tv  !< Thermodynamics structure
  type(verticalGrid_type), intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US  !< A dimensional unit scaling type
  real,                  intent(in)    :: G_earth !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real,                  intent(in)    :: depth !< Depth of ocean column [Z ~> m].
  real,                  intent(in)    :: min_thickness !< Smallest thickness allowed [H ~> m or kg m-2].
  real, dimension(nk),   intent(inout) :: T   !< Layer mean temperature [C ~> degC]
  real, dimension(nk),   intent(in)    :: T_t !< Temperature at top of layer [C ~> degC]
  real, dimension(nk),   intent(in)    :: T_b !< Temperature at bottom of layer [C ~> degC]
  real, dimension(nk),   intent(inout) :: S   !< Layer mean salinity [S ~> ppt]
  real, dimension(nk),   intent(in)    :: S_t !< Salinity at top of layer [S ~> ppt]
  real, dimension(nk),   intent(in)    :: S_b !< Salinity at bottom of layer [S ~> ppt]
  real,                  intent(in)    :: p_surf !< Imposed pressure on ocean at surface [R L2 T-2 ~> Pa]
  real, dimension(nk),   intent(inout) :: h   !< Layer thickness [H ~> m or kg m-2]
  type(remapping_CS),    pointer       :: remap_CS !< Remapping structure for remapping T and S,
                                                   !! if associated
  real,                  intent(in)    :: z_tol !< The tolerance with which to find the depth
                                                !! matching the specified pressure [Z ~> m].
  logical,               intent(in)    :: frac_dp_bugfix !< If true, use bugfix in frac_dp_at_pos

  ! Local variables
  real, dimension(nk+1) :: e ! Top and bottom edge positions for reconstructions [Z ~> m]
  real, dimension(nk) :: h0, h1 ! Initial and remapped layer thicknesses [H ~> m or kg m-2]
  real, dimension(nk) :: S0, S1 ! Initial and remapped layer salinities [S ~> ppt]
  real, dimension(nk) :: T0, T1 ! Initial and remapped layer temperatures [C ~> degC]
  real :: P_t, P_b     ! Top and bottom pressures [R L2 T-2 ~> Pa]
  real :: z_out, e_top ! Interface height positions [Z ~> m]
  real :: min_dz       ! The minimum thickness in depth units [Z ~> m]
  real :: dh_surf_rem  ! The remaining thickness to remove in non-Bousinesq mode [H ~> kg m-2]
  integer :: k

  ! Keep a copy of the initial thicknesses in reverse order to use in remapping
  do k=1,nk ; h0(k) = h(nk+1-k) ; enddo

  if (GV%Boussinesq) then
    min_dz = GV%H_to_Z * min_thickness
    ! Calculate original interface positions
    e(nk+1) = -depth
    do k=nk,1,-1
      e(K) = e(K+1) + GV%H_to_Z*h(k)
    enddo

    P_t = 0.
    e_top = e(1)
    do k=1,nk
      call find_depth_of_pressure_in_cell(T_t(k), T_b(k), S_t(k), S_b(k), e(K), e(K+1), &
                                          P_t, p_surf, GV%Rho0, G_earth, tv%eqn_of_state, &
                                          US, P_b, z_out, z_tol=z_tol, &
                                          frac_dp_bugfix=frac_dp_bugfix)
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
          e_top = e_top - min_dz ! Next interface must be at least this deep
        endif
        ! This layer needs trimming
        h(k) = max( min_thickness, GV%Z_to_H * (e(K) - e(K+1)) )
        if (e(K) < e_top) exit ! No need to go further
      enddo
    endif
  else
    ! In non-Bousinesq mode, we are already in mass units so the calculation is much easier.
    if (p_surf > 0.0) then
      dh_surf_rem = p_surf * GV%RZ_to_H / G_earth
      do k=1,nk
        if (h(k) <= min_thickness) then  ! This layer has no mass to remove.
          cycle
        elseif ((h(k) - min_thickness) < dh_surf_rem) then  ! This layer should be removed entirely.
          dh_surf_rem = dh_surf_rem - (h(k) - min_thickness)
          h(k) = min_thickness
        else  ! This is the last layer that should be removed.
          h(k) = h(k) - dh_surf_rem
          dh_surf_rem = 0.0
          exit
        endif
      enddo
    endif
  endif

  ! Now we need to remap but remapping assumes the surface is at the
  ! same place in the two columns so we turn the column upside down.
  if (associated(remap_CS)) then
    do k=1,nk
      S0(k) = S(nk+1-k)
      T0(k) = T(nk+1-k)
      h1(k) = h(nk+1-k)
    enddo
    call remapping_core_h(remap_CS, nk, h0, T0, nk, h1, T1)
    call remapping_core_h(remap_CS, nk, h0, S0, nk, h1, S1)
    do k=1,nk
      S(k) = S1(nk+1-k)
      T(k) = T1(nk+1-k)
    enddo
  endif

end subroutine cut_off_column_top

!> Initialize horizontal velocity components from file
subroutine initialize_velocity_from_file(u, v, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: u  !< The zonal velocity that is being initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out) :: v  !< The meridional velocity that is being initialized [L T-1 ~> m s-1]
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file to
                                                      !! parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing u or v.
  ! Local variables
  character(len=40)  :: mdl = "initialize_velocity_from_file" ! This subroutine's name.
  character(len=200) :: filename, velocity_file, inputdir ! Strings for file/path
  character(len=64)  :: u_IC_var, v_IC_var ! Velocity component names in files

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mdl, "VELOCITY_FILE", velocity_file, &
                 "The name of the velocity initial condition file.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  filename = trim(velocity_file)
  if (scan(velocity_file, '/')== 0) then ! prepend inputdir if only a filename is given
    filename = trim(inputdir)//trim(velocity_file)
  endif
  if (.not.just_read) call log_param(param_file, mdl, "INPUTDIR/VELOCITY_FILE", filename)

  call get_param(param_file, mdl, "U_IC_VAR", u_IC_var, &
                 "The initial condition variable for zonal velocity in VELOCITY_FILE.", &
                 default="u")
  call get_param(param_file, mdl, "V_IC_VAR", v_IC_var, &
                 "The initial condition variable for meridional velocity in VELOCITY_FILE.", &
                 default="v")

  if (just_read) return ! All run-time parameters have been read, so return.

  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
         " initialize_velocity_from_file: Unable to open "//trim(filename))

  !  Read the velocities from a netcdf file.
  call MOM_read_vector(filename, u_IC_var, v_IC_var, u(:,:,:), v(:,:,:), G%Domain, scale=US%m_s_to_L_T)

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_velocity_from_file

!> Initialize horizontal velocity components to zero.
subroutine initialize_velocity_zero(u, v, G, GV, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: u  !< The zonal velocity that is being initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out) :: v  !< The meridional velocity that is being initialized [L T-1 ~> m s-1]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file to
                                                      !! parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.
  ! Local variables
  character(len=200) :: mdl = "initialize_velocity_zero" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

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
subroutine initialize_velocity_uniform(u, v, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: u  !< The zonal velocity that is being initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out) :: v  !< The meridional velocity that is being initialized [L T-1 ~> m s-1]
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file to
                                                      !! parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing u or v.
  ! Local variables
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  real    :: initial_u_const, initial_v_const ! Constant initial velocities [L T-1 ~> m s-1]
  character(len=200) :: mdl = "initialize_velocity_uniform" ! This subroutine's name.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

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
subroutine initialize_velocity_circular(u, v, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: u  !< The zonal velocity that is being initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out) :: v  !< The meridional velocity that is being initialized [L T-1 ~> m s-1]
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file to
                                                      !! parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing u or v.
  ! Local variables
  character(len=200) :: mdl = "initialize_velocity_circular"
  real :: circular_max_u ! The amplitude of the zonal flow [L T-1 ~> m s-1]
  real :: dpi        ! A local variable storing pi = 3.14159265358979... [nondim]
  real :: psi1, psi2 ! Values of the streamfunction at two points [L2 T-1 ~> m2 s-1]
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  call get_param(param_file, mdl, "CIRCULAR_MAX_U", circular_max_u, &
                 "The amplitude of zonal flow from which to scale the "// &
                 "circular stream function [m s-1].", &
                 units="m s-1", default=0., scale=US%m_s_to_L_T, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  if (G%grid_unit_to_L <= 0.) call MOM_error(FATAL, "MOM_state_initialization.F90: "//&
          "initialize_velocity_circular() is only set to work with Cartesian axis units.")

  dpi = acos(0.0)*2.0 ! pi

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
    r = sqrt( (x**2) + (y**2) ) ! Circular stream function is a function of radius only
    r = min(1.0, r) ! Flatten stream function in corners of box
    my_psi = 0.5*(1.0 - cos(dpi*r))
    my_psi = my_psi * (circular_max_u * G%len_lon * G%grid_unit_to_L / dpi) ! len_lon is in km
  end function my_psi

end subroutine initialize_velocity_circular

!> Initializes temperature and salinity from file
subroutine initialize_temp_salt_from_file(T, S, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< The potential temperature that is
                                                               !! being initialized [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< The salinity that is
                                                               !! being initialized [S ~> ppt]
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file !< A structure to parse for run-time parameters
  logical,                                   intent(in)  :: just_read !< If true, this call will only
                                                           !! read parameters without changing T or S.
  ! Local variables
  character(len=200) :: filename, salt_filename ! Full paths to input files
  character(len=200) :: ts_file, salt_file, inputdir ! Strings for file/path
  character(len=40)  :: mdl = "initialize_temp_salt_from_file"
  character(len=64)  :: temp_var, salt_var ! Temperature and salinity names in files

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mdl, "TS_FILE", ts_file, &
                 "The initial condition file for temperature.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  filename = trim(ts_file)
  if (scan(ts_file, '/')== 0) then ! prepend inputdir if only a filename is given
    filename = trim(inputdir)//trim(ts_file)
  endif
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
  call MOM_read_data(filename, temp_var, T(:,:,:), G%Domain, scale=US%degC_to_C)

  salt_filename = trim(salt_file)
  if (scan(salt_file, '/')== 0) then ! prepend inputdir if only a filename is given
    salt_filename = trim(inputdir)//trim(salt_file)
  endif
  if (.not.file_exists(salt_filename, G%Domain)) call MOM_error(FATAL, &
     " initialize_temp_salt_from_file: Unable to open "//trim(salt_filename))

  call MOM_read_data(salt_filename, salt_var, S(:,:,:), G%Domain, scale=US%ppt_to_S)

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_temp_salt_from_file

!> Initializes temperature and salinity from a 1D profile
subroutine initialize_temp_salt_from_profile(T, S, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< The potential temperature that is
                                                               !! being initialized [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< The salinity that is
                                                               !! being initialized [S ~> ppt]
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file !< A structure to parse for run-time parameters
  logical,                                   intent(in)  :: just_read !< If true, this call will only read
                                                               !! parameters without changing T or S.
  ! Local variables
  real, dimension(SZK_(GV)) :: T0 ! The profile of temperatures [C ~> degC]
  real, dimension(SZK_(GV)) :: S0 ! The profile of salinities [S ~> ppt]
  integer :: i, j, k
  character(len=200) :: filename, ts_file, inputdir ! Strings for file/path
  character(len=64)  :: temp_var, salt_var ! Temperature and salinity names in files
  character(len=40)  :: mdl = "initialize_temp_salt_from_profile"

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mdl, "TS_FILE", ts_file, &
                 "The file with the reference profiles for temperature and salinity.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "TEMP_IC_VAR", temp_var, &
                 "The initial condition variable for potential temperature.", &
                 default="PTEMP", do_not_log=just_read)
  call get_param(param_file, mdl, "SALT_IC_VAR", salt_var, &
                 "The initial condition variable for salinity.", &
                 default="SALT", do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  filename = trim(inputdir)//trim(ts_file)
  call log_param(param_file, mdl, "INPUTDIR/TS_FILE", filename)
   if (.not.file_exists(filename)) call MOM_error(FATAL, &
     " initialize_temp_salt_from_profile: Unable to open "//trim(filename))

  ! Read the temperatures and salinities from a netcdf file.
  call MOM_read_data(filename, temp_var, T0(:), scale=US%degC_to_C)
  call MOM_read_data(filename, salt_var, S0(:), scale=US%ppt_to_S)

  do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    T(i,j,k) = T0(k) ; S(i,j,k) = S0(k)
  enddo ; enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_temp_salt_from_profile

!> Initializes temperature and salinity by fitting to density
subroutine initialize_temp_salt_fit(T, S, G, GV, US, param_file, eqn_of_state, P_Ref, just_read)
  type(ocean_grid_type),   intent(in)  :: G            !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV           !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T !< The potential temperature that is
                                                       !! being initialized [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S !< The salinity that is being
                                                       !! initialized [S ~> ppt].
  type(unit_scale_type),   intent(in)  :: US           !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file   !< A structure to parse for run-time
                                                       !! parameters.
  type(EOS_type),          intent(in)  :: eqn_of_state !< Equation of state structure
  real,                    intent(in)  :: P_Ref        !< The coordinate-density reference pressure
                                                       !! [R L2 T-2 ~> Pa].
  logical,                 intent(in)  :: just_read    !< If true, this call will only read
                                                       !! parameters without changing T or S.
  ! Local variables
  real :: T0(SZK_(GV))  ! Layer potential temperatures [C ~> degC]
  real :: S0(SZK_(GV))  ! Layer salinities [S ~> ppt]
  real :: T_Ref         ! Reference Temperature [C ~> degC]
  real :: S_Ref         ! Reference Salinity [S ~> ppt]
  real :: pres(SZK_(GV))     ! An array of the reference pressure [R L2 T-2 ~> Pa].
  real :: drho_dT(SZK_(GV))  ! Derivative of density with temperature [R C-1 ~> kg m-3 degC-1].
  real :: drho_dS(SZK_(GV))  ! Derivative of density with salinity [R S-1 ~> kg m-3 ppt-1].
  real :: rho_guess(SZK_(GV)) ! Potential density at T0 & S0 [R ~> kg m-3].
  logical :: fit_salin       ! If true, accept the prescribed temperature and fit the salinity.
  character(len=40)  :: mdl = "initialize_temp_salt_fit" ! This subroutine's name.
  integer :: i, j, k, itt, nz
  nz = GV%ke

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  call get_param(param_file, mdl, "T_REF", T_Ref, &
                 "A reference temperature used in initialization.", &
                 units="degC", scale=US%degC_to_C, fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "S_REF", S_Ref, &
                 "A reference salinity used in initialization.", &
                 units="ppt", default=35.0, scale=US%ppt_to_S, do_not_log=just_read)
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
subroutine initialize_temp_salt_linear(T, S, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G  !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< The potential temperature that is
                                                               !! being initialized [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< The salinity that is
                                                               !! being initialized [S ~> ppt]
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file !< A structure to parse for
                                                               !! run-time parameters
  logical,                                   intent(in)  :: just_read !< If present and true,
                                                               !! this call will only read parameters
                                                               !! without changing T or S.

  ! Local variables
  real :: S_top, S_range ! Reference salinity in the surface layer and its vertical range [S ~> ppt]
  real :: T_top, T_range ! Reference temperature in the surface layer and its vertical range [C ~> degC]
  character(len=40)  :: mdl = "initialize_temp_salt_linear" ! This subroutine's name.
  integer :: k

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")
  call get_param(param_file, mdl, "T_TOP", T_top, &
                 "Initial temperature of the top surface.", &
                 units="degC", scale=US%degC_to_C, fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "T_RANGE", T_range, &
                 "Initial temperature difference (top-bottom).", &
                 units="degC", scale=US%degC_to_C, fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "S_TOP", S_top, &
                 "Initial salinity of the top surface.", &
                 units="ppt", scale=US%ppt_to_S, fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "S_RANGE", S_range, &
                 "Initial salinity difference (top-bottom).", &
                 units="ppt", scale=US%ppt_to_S, fail_if_missing=.not.just_read, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  ! Prescribe salinity and temperature, with the extrapolated top interface value prescribed.
  do k=1,GV%ke
    S(:,:,k) = S_top - S_range*((real(k)-0.5)/real(GV%ke))
    T(:,:,k) = T_top - T_range*((real(k)-0.5)/real(GV%ke))
  enddo

  ! Prescribe salinity and temperature, but with the top layer value matching the surface value.
  ! S(:,:,1) = S_top ; T(:,:,1) = T_top
  ! do k=2,GV%ke
  !   S(:,:,k) = S_top - S_range * (real(k-1) / real(GV%ke-1))
  !   T(:,:,k) = T_top - T_range * (real(k-1) / real(GV%ke-1))
  ! enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_temp_salt_linear

!> This subroutine sets the inverse restoration time (Idamp), and
!! the values towards which the interface heights and an arbitrary
!! number of tracers should be restored within each sponge. The
!! interface height is always subject to damping, and must always be
!! the first registered field.
subroutine initialize_sponges_file(G, GV, US, use_temperature, tv, u, v, depth_tot, param_file, &
                                   Layer_CSp, ALE_CSp, Time)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  logical,                 intent(in) :: use_temperature !< If true, T & S are state variables.
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic
                                              !! variables.
  real, target, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: u    !< The zonal velocity that is being
                                              !! initialized [L T-1 ~> m s-1]
  real, target, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in) :: v    !< The meridional velocity that is being
                                              !! initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters.
  type(sponge_CS),         pointer    :: Layer_CSp  !< A pointer that is set to point to the control
                                              !! structure for this module (in layered mode).
  type(ALE_sponge_CS),     pointer    :: ALE_CSp  !< A pointer that is set to point to the control
                                                  !! structure for this module (in ALE mode).
  type(time_type),         intent(in) :: Time !< Time at the start of the run segment. Time_in
                                              !! overrides any value set for Time.
  ! Local variables
  real, allocatable, dimension(:,:,:) :: eta ! The target interface heights [Z ~> m].
  real, allocatable, dimension(:,:,:) :: dz  ! The target interface thicknesses in height units [Z ~> m]
  real, allocatable, dimension(:,:,:) :: h   ! The target interface thicknesses [H ~> m or kg m-2].

  real, dimension (SZI_(G),SZJ_(G),SZK_(GV)) :: &
    tmp, &    ! A temporary array for temperatures [C ~> degC] or other tracers.
    tmp2      ! A temporary array for salinities [S ~> ppt]
  real, dimension (SZI_(G),SZJ_(G)) :: &
    tmp_2d    ! A temporary array for mixed layer densities [R ~> kg m-3]
  real, allocatable, dimension(:,:,:) :: tmp_T ! A temporary array for reading sponge target temperatures
                                    ! on the vertical grid of the input file  [C ~> degC]
  real, allocatable, dimension(:,:,:) :: tmp_S ! A temporary array for reading sponge target salinities
                                    ! on the vertical grid of the input file [S ~> ppt]
  real, allocatable, dimension(:,:,:) :: tmp_u ! Temporary array for reading sponge target zonal
                                    ! velocities on the vertical grid of the input file [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:) :: tmp_v ! Temporary array for reading sponge target meridional
                                    ! velocities on the vertical grid of the input file [L T-1 ~> m s-1]

  real :: Idamp(SZI_(G),SZJ_(G))    ! The sponge damping rate [T-1 ~> s-1]
  real :: Idamp_u(SZIB_(G),SZJ_(G)) ! The sponge damping rate for velocity fields [T-1 ~> s-1]
  real :: Idamp_v(SZI_(G),SZJB_(G)) ! The sponge damping rate for velocity fields [T-1 ~> s-1]
  real :: pres(SZI_(G))             ! An array of the reference pressure [R L2 T-2 ~> Pa]

  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, nz
  integer :: isd, ied, jsd, jed
  integer, dimension(4) :: siz
  integer :: nz_data  ! The size of the sponge source grid
  logical :: sponge_uv ! Apply sponges in u and v, in addition to tracers.
  character(len=40) :: potemp_var, salin_var, u_var, v_var, Idamp_var, Idamp_u_var, Idamp_v_var, eta_var
  character(len=40) :: mdl = "initialize_sponges_file"
  character(len=200) :: damping_file, uv_damping_file, state_file, state_uv_file  ! Strings for filenames
  character(len=200) :: filename, inputdir ! Strings for file/path and path.

  logical :: use_ALE ! True if ALE is being used, False if in layered mode
  logical :: time_space_interp_sponge ! If true use sponge data that need to be interpolated in both
                              ! the horizontal dimension and in time prior to vertical remapping.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  pres(:) = 0.0 ; tmp(:,:,:) = 0.0 ; Idamp(:,:) = 0.0 ; Idamp_u(:,:) = 0.0 ; Idamp_v(:,:) = 0.0

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
  call get_param(param_file, mdl, "SPONGE_UV", sponge_uv, &
                 "Apply sponges in u and v, in addition to tracers.", &
                 default=.false.)
  if (sponge_uv) then
    call get_param(param_file, mdl, "SPONGE_UV_STATE_FILE", state_uv_file, &
                 "The name of the file with the state to damp UV toward.", &
                 default=damping_file)
    call get_param(param_file, mdl, "SPONGE_U_VAR", u_var, &
                 "The name of the zonal velocity variable in "//&
                 "SPONGE_UV_STATE_FILE.", default="UVEL")
    call get_param(param_file, mdl, "SPONGE_V_VAR", v_var, &
                 "The name of the vertical velocity variable in "//&
                 "SPONGE_UV_STATE_FILE.", default="VVEL")
  endif
  call get_param(param_file, mdl, "SPONGE_ETA_VAR", eta_var, &
                 "The name of the interface height variable in "//&
                 "SPONGE_STATE_FILE.", default="ETA")
  call get_param(param_file, mdl, "SPONGE_IDAMP_VAR", Idamp_var, &
                 "The name of the inverse damping rate variable in "//&
                 "SPONGE_DAMPING_FILE.", default="Idamp")
  if (sponge_uv) then
    call get_param(param_file, mdl, "SPONGE_UV_DAMPING_FILE", uv_damping_file, &
                   "The name of the file with sponge damping rates for the velocity variables.", &
                   default=damping_file)
    call get_param(param_file, mdl, "SPONGE_IDAMP_U_var", Idamp_u_var, &
                   "The name of the inverse damping rate variable in "//&
                   "SPONGE_UV_DAMPING_FILE for the velocities.", default=Idamp_var)
    call get_param(param_file, mdl, "SPONGE_IDAMP_V_var", Idamp_v_var, &
                   "The name of the inverse damping rate variable in "//&
                   "SPONGE_UV_DAMPING_FILE for the velocities.", default=Idamp_var)
  endif
  call get_param(param_file, mdl, "USE_REGRIDDING", use_ALE, default=.false., do_not_log=.true.)

  call get_param(param_file, mdl, "INTERPOLATE_SPONGE_TIME_SPACE", time_space_interp_sponge, &
                 "If True, perform on-the-fly regridding in lat-lon-time of sponge restoring data.", &
                 default=.false.)

  ! Read in sponge damping rate for tracers
  filename = trim(damping_file)
  if (scan(damping_file, '/')== 0) then ! prepend inputdir if only a filename is given
    filename = trim(inputdir)//trim(damping_file)
  endif
  call log_param(param_file, mdl, "INPUTDIR/SPONGE_DAMPING_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) &
    call MOM_error(FATAL, " initialize_sponges: Unable to open "//trim(filename))

  if (time_space_interp_sponge .and. .not. use_ALE) &
    call MOM_error(FATAL, " initialize_sponges: Time-varying sponges are currently unavailable in layered mode ")

  call MOM_read_data(filename, Idamp_var, Idamp(:,:), G%Domain, scale=US%T_to_s)

  ! Read in sponge damping rate for velocities
  if (sponge_uv) then
    if (separate_idamp_for_uv()) then
      filename = trim(inputdir)//trim(uv_damping_file)
      call log_param(param_file, mdl, "INPUTDIR/SPONGE_UV_DAMPING_FILE", filename)

      if (.not.file_exists(filename, G%Domain)) &
        call MOM_error(FATAL, " initialize_sponges: Unable to open "//trim(filename))

      call MOM_read_vector(filename, Idamp_u_var,Idamp_v_var,Idamp_u(:,:),Idamp_v(:,:), G%Domain, scale=US%T_to_s)
    else
      !      call MOM_error(FATAL, "Must provide SPONGE_IDAMP_U_var and SPONGE_IDAMP_V_var")
      call pass_var(Idamp,G%Domain)
      do j=G%jsc,G%jec
        do i=G%iscB,G%iecB
          Idamp_u(I,j) = 0.5*(Idamp(i,j)+Idamp(i+1,j))
        enddo
      enddo
      do j=G%jscB,G%jecB
        do i=G%isc,G%iec
          Idamp_v(i,J) = 0.5*(Idamp(i,j)+Idamp(i,j+1))
        enddo
      enddo
    endif
  endif

  ! Now register all of the fields which are damped in the sponge.
  ! By default, momentum is advected vertically within the sponge, but
  ! momentum is typically not damped within the sponge.

  filename = trim(inputdir)//trim(state_file)
  call log_param(param_file, mdl, "INPUTDIR/SPONGE_STATE_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) &
    call MOM_error(FATAL, " initialize_sponges: Unable to open "//trim(filename))


  if (.not. use_ALE) then
    ! The first call to set_up_sponge_field is for the interface heights if in layered mode.
    allocate(eta(isd:ied,jsd:jed,nz+1), source=0.0)
    call MOM_read_data(filename, eta_var, eta(:,:,:), G%Domain, scale=US%m_to_Z)

    do j=js,je ; do i=is,ie
      eta(i,j,nz+1) = -depth_tot(i,j)
    enddo ; enddo
    do k=nz,1,-1 ; do j=js,je ; do i=is,ie
      if (eta(i,j,K) < (eta(i,j,K+1) + GV%Angstrom_Z)) &
        eta(i,j,K) = eta(i,j,K+1) + GV%Angstrom_Z
    enddo ; enddo ; enddo
    ! Set the sponge damping rates so that the model will know where to
    ! apply the sponges, along with the interface heights.
    call initialize_sponge(Idamp, eta, G, param_file, Layer_CSp, GV)
    deallocate(eta)

    if ( GV%nkml>0) then
      ! This call to set_up_sponge_ML_density registers the target values of the
      ! mixed layer density, which is used in determining which layers can be
      ! inflated without causing static instabilities.
      do i=is,ie ; pres(i) = tv%P_Ref ; enddo
      EOSdom(:) = EOS_domain(G%HI)

      call MOM_read_data(filename, potemp_var, tmp(:,:,:), G%Domain, scale=US%degC_to_C)
      call MOM_read_data(filename, salin_var, tmp2(:,:,:), G%Domain, scale=US%ppt_to_S)

      do j=js,je
        call calculate_density(tmp(:,j,1), tmp2(:,j,1), pres, tmp_2d(:,j), tv%eqn_of_state, EOSdom)
      enddo

      call set_up_sponge_ML_density(tmp_2d, G, Layer_CSp)
    endif

   ! Now register all of the tracer fields which are damped in the
   ! sponge. By default, momentum is advected vertically within the
   ! sponge, but momentum is typically not damped within the sponge.


    ! The remaining calls to set_up_sponge_field can be in any order.
    if ( use_temperature) then
      call MOM_read_data(filename, potemp_var, tmp(:,:,:), G%Domain, scale=US%degC_to_C)
      call set_up_sponge_field(tmp, tv%T, G, GV, nz, Layer_CSp)
      call MOM_read_data(filename, salin_var, tmp2(:,:,:), G%Domain, scale=US%ppt_to_S)
      call set_up_sponge_field(tmp2, tv%S, G, GV, nz, Layer_CSp)
    endif

!  else
    ! Initialize sponges without supplying sponge grid
!    if (sponge_uv) then
!      call initialize_ALE_sponge(Idamp, G, GV, US, param_file, ALE_CSp, Idamp_u, Idamp_v)
!    else
!      call initialize_ALE_sponge(Idamp, G, GV, US, param_file, ALE_CSp)
!    endif
  endif


  if  (use_ALE) then ! ALE mode
    if (.not. time_space_interp_sponge) then
      call field_size(filename,eta_var,siz,no_domain=.true.)
      if (siz(1) /= G%ieg-G%isg+1 .or. siz(2) /= G%jeg-G%jsg+1) &
        call MOM_error(FATAL,"initialize_sponge_file: Array size mismatch for sponge data.")
      nz_data = siz(3)-1
      allocate(eta(isd:ied,jsd:jed,nz_data+1))
      allocate(dz(isd:ied,jsd:jed,nz_data))
      call MOM_read_data(filename, eta_var, eta(:,:,:), G%Domain, scale=US%m_to_Z)
      do j=js,je ; do i=is,ie
        eta(i,j,nz_data+1) = -depth_tot(i,j)
      enddo ; enddo
      do k=nz_data,1,-1 ; do j=js,je ; do i=is,ie
        if (eta(i,j,K) < (eta(i,j,K+1) + GV%Angstrom_Z)) &
          eta(i,j,K) = eta(i,j,K+1) + GV%Angstrom_Z
      enddo ; enddo ; enddo
      do k=1,nz_data ; do j=js,je ; do i=is,ie
        dz(i,j,k) = eta(i,j,k)-eta(i,j,k+1)
      enddo; enddo ; enddo
      deallocate(eta)

      if (use_temperature) then
        allocate(tmp_T(isd:ied,jsd:jed,nz_data))
        allocate(tmp_S(isd:ied,jsd:jed,nz_data))
        call MOM_read_data(filename, potemp_var, tmp_T(:,:,:), G%Domain, scale=US%degC_to_C)
        call MOM_read_data(filename, salin_var, tmp_S(:,:,:), G%Domain, scale=US%ppt_to_S)
      endif

      if (sponge_uv) then
        call initialize_ALE_sponge(Idamp, G, GV, param_file, ALE_CSp, dz, nz_data, Idamp_u, Idamp_v, &
                                   data_h_is_Z=.true.)
      else
        call initialize_ALE_sponge(Idamp, G, GV, param_file, ALE_CSp, dz, nz_data, &
                                   data_h_is_Z=.true.)
      endif
      if (use_temperature) then
        call set_up_ALE_sponge_field(tmp_T, G, GV, tv%T, ALE_CSp, 'temp', &
               sp_long_name='temperature', sp_unit='degC s-1')
        call set_up_ALE_sponge_field(tmp_S, G, GV, tv%S, ALE_CSp, 'salt', &
               sp_long_name='salinity', sp_unit='g kg-1 s-1')
        deallocate(tmp_S)
        deallocate(tmp_T)
      endif
      deallocate(dz)

      if (sponge_uv) then
        filename = trim(inputdir)//trim(state_uv_file)
        call log_param(param_file, mdl, "INPUTDIR/SPONGE_STATE_UV_FILE", filename)
        if (.not.file_exists(filename, G%Domain)) &
             call MOM_error(FATAL, " initialize_sponges: Unable to open "//trim(filename))
        allocate(tmp_u(G%IsdB:G%IedB,jsd:jed,nz_data))
        allocate(tmp_v(isd:ied,G%JsdB:G%JedB,nz_data))
        call MOM_read_vector(filename, u_var, v_var, tmp_u(:,:,:), tmp_v(:,:,:), G%Domain, scale=US%m_s_to_L_T)
        call set_up_ALE_sponge_vel_field(tmp_u, tmp_v, G, GV, u, v, ALE_CSp)
        deallocate(tmp_u,tmp_v)
      endif
    else
      ! Initialize sponges without supplying sponge grid
      if (sponge_uv) then
        call initialize_ALE_sponge(Idamp, G, GV, US, param_file, ALE_CSp, Idamp_u, Idamp_v)
      else
        call initialize_ALE_sponge(Idamp, G, GV, US, param_file, ALE_CSp)
      endif
      ! The remaining calls to set_up_sponge_field can be in any order.
      if ( use_temperature) then
        call set_up_ALE_sponge_field(filename, potemp_var, Time, G, GV, US, tv%T, ALE_CSp, &
               'temp', sp_long_name='temperature', sp_unit='degC s-1', scale=US%degC_to_C)
        call set_up_ALE_sponge_field(filename, salin_var, Time, G, GV, US, tv%S, ALE_CSp, &
               'salt', sp_long_name='salinity', sp_unit='g kg-1 s-1', scale=US%ppt_to_S)
      endif
      if (sponge_uv) then
        filename = trim(inputdir)//trim(state_uv_file)
        call log_param(param_file, mdl, "INPUTDIR/SPONGE_STATE_UV_FILE", filename)
        if (.not.file_exists(filename, G%Domain)) &
             call MOM_error(FATAL, " initialize_sponges: Unable to open "//trim(filename))
        call set_up_ALE_sponge_vel_field(filename, u_var, filename, v_var, Time, G, GV, US, &
                                         ALE_CSp, u, v, scale=US%m_s_to_L_T)
      endif
    endif
  endif

  if (sponge_uv .and. .not. use_ALE) call MOM_error(FATAL,'initialize_sponges_file: '// &
                       'UV damping to target values only available in ALE mode')


  contains

  ! returns true if a separate idamp is provided for u and/or v
  logical function separate_idamp_for_uv()
    separate_idamp_for_uv = (lowercase(damping_file)/=lowercase(uv_damping_file) .or. &
         lowercase(Idamp_var)/=lowercase(Idamp_u_var) .or. lowercase(Idamp_var)/=lowercase(Idamp_v_var))
  end function separate_idamp_for_uv

end subroutine initialize_sponges_file

subroutine initialize_oda_incupd_file(G, GV, US, use_temperature, tv, h, u, v, param_file, &
                                      oda_incupd_CSp, restart_CS, Time)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  logical,                 intent(in)    :: use_temperature !< If true, T & S are state variables.
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure pointing to various thermodynamic
                                                 !! variables.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2] (in)

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                             intent(in) :: u     !< The zonal velocity that is being
                                                 !! initialized [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                             intent(in) :: v     !< The meridional velocity that is being
                                                 !! initialized [L T-1 ~> m s-1]
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters.
  type(oda_incupd_CS),     pointer    :: oda_incupd_CSp  !< A pointer that is set to point to the control
                                                 !! structure for this module.
  type(MOM_restart_CS),    intent(in) :: restart_CS !< MOM restart control structure
  type(time_type),         intent(in) :: Time    !< Time at the start of the run segment. Time_in
                                                 !! overrides any value set for Time.
  ! Local variables
  real, allocatable, dimension(:,:,:) :: hoda ! The layer thickness increment and oda layer thickness [H ~> m or kg m-2]
  real, allocatable, dimension(:,:,:) :: tmp_tr ! A temporary array for reading oda tracer increments
                                    ! on the vertical grid of the input file, used for both
                                    ! temperatures [C ~> degC] and salinities [S ~> ppt]
  real, allocatable, dimension(:,:,:) :: tmp_u ! Temporary array for reading oda zonal velocity
                                    ! increments on the vertical grid of the input file [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:) :: tmp_v ! Temporary array for reading oda meridional velocity
                                    ! increments on the vertical grid of the input file [L T-1 ~> m s-1]

  integer :: is, ie, js, je, nz
  integer :: isd, ied, jsd, jed

  integer, dimension(4) :: siz
  integer :: nz_data  ! The size of the sponge source grid
  logical :: oda_inc  ! input files are increments (true) or full fields (false)
  logical :: save_inc ! save increments if using full fields
  logical :: uv_inc   ! use u and v increments
  logical :: reset_ncount ! reset ncount to zero if true

  character(len=40)  :: tempinc_var, salinc_var, uinc_var, vinc_var, h_var
  character(len=40)  :: mdl = "initialize_oda_incupd_file"
  character(len=200) :: inc_file, uv_inc_file  ! Strings for filenames
  character(len=200) :: filename, inputdir ! Strings for file/path and path.

!  logical :: use_ALE ! True if ALE is being used, False if in layered mode

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  call get_param(param_file, mdl, "ODA_INCUPD_FILE", inc_file, &
                 "The name of the file with the T,S,h increments.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "ODA_INCUPD_INC", oda_inc, &
                 "INCUPD files are increments and not full fields.", &
                 default=.true.)
  if (.not.oda_inc) then
    call get_param(param_file, mdl, "ODA_INCUPD_SAVE", save_inc, &
                   "If true, save the increments when using full fields.", &
                   default=.false.)
  endif
  call get_param(param_file, mdl, "ODA_INCUPD_RESET_NCOUNT", reset_ncount, &
                 "If True, reinitialize number of updates already done, ncount.", &
                 default=.true.)
  if (.not.oda_inc .and. .not.reset_ncount) &
    call MOM_error(FATAL, " initialize_oda_incupd: restarting during update "// &
                 "necessitates increments input file")

  call get_param(param_file, mdl, "ODA_TEMPINC_VAR", tempinc_var, &
                 "The name of the potential temperature inc. variable in "//&
                 "ODA_INCUPD_FILE.", default="ptemp_inc")
  call get_param(param_file, mdl, "ODA_SALTINC_VAR", salinc_var, &
                 "The name of the salinity inc. variable in "//&
                 "ODA_INCUPD_FILE.", default="sal_inc")
  call get_param(param_file, mdl, "ODA_THK_VAR", h_var, &
                 "The name of the layer thickness variable in "//&
                 "ODA_INCUPD_FILE.", default="h")
  call get_param(param_file, mdl, "ODA_INCUPD_UV", uv_inc, &
                 "use U,V increments.", &
                 default=.true.)
  call get_param(param_file, mdl, "ODA_INCUPD_UV_FILE", uv_inc_file, &
                 "The name of the file with the U,V increments.", &
                 default=inc_file)
  call get_param(param_file, mdl, "ODA_UINC_VAR", uinc_var, &
                 "The name of the zonal vel. inc. variable in "//&
                 "ODA_INCUPD_FILE.", default="u_inc")
  call get_param(param_file, mdl, "ODA_VINC_VAR", vinc_var, &
                 "The name of the meridional vel. inc. variable in "//&
                 "ODA_INCUPD_FILE.", default="v_inc")

!  call get_param(param_file, mdl, "USE_REGRIDDING", use_ALE, default=.false., do_not_log=.true.)

  ! Read in incremental update for tracers
  filename = trim(inc_file)
  if (scan(inc_file, '/')== 0) then ! prepend inputdir if only a filename is given
    filename = trim(inputdir)//trim(inc_file)
  endif
  call log_param(param_file, mdl, "INPUTDIR/ODA_INCUPD_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) &
    call MOM_error(FATAL, " initialize_oda_incupd: Unable to open "//trim(filename))

  call field_size(filename,h_var,siz,no_domain=.true.)
  if (siz(1) /= G%ieg-G%isg+1 .or. siz(2) /= G%jeg-G%jsg+1) &
         call MOM_error(FATAL,"initialize_oda_incupd_file: Array size mismatch for oda data.")
  nz_data = siz(3)
  ! get h increments
  allocate(hoda(isd:ied,jsd:jed,nz_data))
  call MOM_read_data(filename, h_var   , hoda(:,:,:), G%Domain, scale=US%m_to_Z)
  call initialize_oda_incupd( G, GV, US, param_file, oda_incupd_CSp, hoda, nz_data, restart_CS)
  deallocate(hoda)

  ! set-up T and S increments arrays
  if (use_temperature) then
    allocate(tmp_tr(isd:ied,jsd:jed,nz_data))
    ! temperature inc. in array Inc(1)
    call MOM_read_data(filename, tempinc_var, tmp_tr(:,:,:), G%Domain, scale=US%degC_to_C)
    call set_up_oda_incupd_field(tmp_tr, G, GV, oda_incupd_CSp)
    ! salinity inc. in array Inc(2)
    call MOM_read_data(filename, salinc_var, tmp_tr(:,:,:), G%Domain, scale=US%ppt_to_S)
    call set_up_oda_incupd_field(tmp_tr, G, GV, oda_incupd_CSp)
    deallocate(tmp_tr)
  endif

  ! set-up U and V increments arrays
  if (uv_inc) then
    filename = trim(inputdir)//trim(uv_inc_file)
    call log_param(param_file, mdl, "INPUTDIR/ODA_INCUPD_UV_FILE", filename)
    if (.not.file_exists(filename, G%Domain)) &
            call MOM_error(FATAL, " initialize_oda_incupd_uv: Unable to open "//trim(filename))
    allocate(tmp_u(G%IsdB:G%IedB,jsd:jed,nz_data), source=0.0)
    allocate(tmp_v(isd:ied,G%JsdB:G%JedB,nz_data), source=0.0)
    call MOM_read_vector(filename, uinc_var, vinc_var, tmp_u, tmp_v, G%Domain, scale=US%m_s_to_L_T)
    call set_up_oda_incupd_vel_field(tmp_u, tmp_v, G, GV, oda_incupd_CSp)
    deallocate(tmp_u, tmp_v)
  endif

  ! calculate increments if input are full fields
  if (oda_inc) then ! input are increments
    if (is_root_pe()) call MOM_mesg("incupd using increments fields ")
  else ! inputs are full fields
    if (is_root_pe()) call MOM_mesg("incupd using full fields ")
    call calc_oda_increments(h, tv, u, v, G, GV, US, oda_incupd_CSp)
    if (save_inc) then
      call output_oda_incupd_inc(Time, G, GV, param_file, oda_incupd_CSp, US)
    endif
  endif  ! not oda_inc

end subroutine initialize_oda_incupd_file


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
subroutine MOM_temp_salt_initialize_from_Z(h, tv, depth_tot, G, GV, US, PF, just_read, frac_shelf_h)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: h    !< Layer thicknesses being initialized [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv   !< A structure pointing to various thermodynamic
                                                 !! variables including temperature and salinity
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: PF   !< A structure indicating the open file
                                                 !! to parse for model parameter values.
  logical,                 intent(in)    :: just_read !< If true, this call will only read
                                                 !! parameters without changing T or S.
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: frac_shelf_h  !< The fraction of the grid cell covered
                                                 !! by a floating ice shelf [nondim].

  ! Local variables
  character(len=200) :: filename   !< The name of an input file containing temperature
                                   !! and salinity in z-space; by default it is also used for ice shelf area.
  character(len=200) :: tfilename  !< The name of an input file containing temperature in z-space.
  character(len=200) :: sfilename  !< The name of an input file containing salinity in z-space.
  character(len=200) :: inputdir   !! The directory where NetCDF input files are.
  character(len=200) :: mesg

  type(EOS_type), pointer :: eos => NULL()
  type(thermo_var_ptrs) :: tv_loc   ! A temporary thermo_var container
  type(verticalGrid_type) :: GV_loc ! A temporary vertical grid structure
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_initialize_layers_from_Z" ! This module's name.

  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: is, ie, js, je, nz ! compute domain indices
  integer :: isg, ieg, jsg, jeg ! global extent
  integer :: isd, ied, jsd, jed ! data domain indices

  integer :: i, j, k, ks
  integer :: nkml     ! The number of layers in the mixed layer.

  integer :: inconsistent ! The total number of cells with in consistent topography and layer thicknesses.
  integer :: kd       ! The number of levels in the input data
  integer :: nkd      ! number of levels to use for regridding input arrays
  real    :: eps_Z    ! A negligibly thin layer thickness [Z ~> m].
  real    :: eps_rho  ! A negligibly small density difference [R ~> kg m-3].
  real    :: PI_180   ! for conversion from degrees to radians [radian degree-1]
  real    :: Hmix_default ! The default initial mixed layer depth [Z ~> m].
  real    :: Hmix_depth   ! The mixed layer depth in the initial condition [Z ~> m].
  real    :: missing_value_temp  ! The missing value in the input temperature field [C ~> degC]
  real    :: missing_value_salt  ! The missing value in the input salinity field [S ~> ppt]
  real    :: tol_temp ! The tolerance for changes in temperature during the horizontal
                      ! interpolation from an input dataset [C ~> degC]
  real    :: tol_sal  ! The tolerance for changes in salinity during the horizontal
                      ! interpolation from an input dataset [S ~> ppt]
  logical :: correct_thickness  ! If true, correct the column thicknesses to match the topography
  real    :: h_tolerance ! A parameter that controls the tolerance when adjusting the
                         ! thickness to fit the bathymetry [Z ~> m].
  real    :: tol_dz_bot  ! A tolerance for detecting inconsistent bottom depths when
                         ! correct_thickness is false [Z ~> m]
  character(len=40) :: potemp_var, salin_var

  integer, parameter :: niter=10   ! number of iterations for t/s adjustment to layer density
  logical            :: adjust_temperature = .true.  ! fit t/s to target densities
  real    :: temp_land_fill  ! A temperature value to use for land points [C ~> degC]
  real    :: salt_land_fill  ! A salinity value to use for land points [C ~> degC]

  ! data arrays
  real, dimension(:), allocatable :: z_edges_in ! Input data interface heights or depths [Z ~> m]
  real, dimension(:), allocatable :: z_in       ! Input data cell heights or depths [Z ~> m]
  real, dimension(:), allocatable :: Rb         ! Interface densities [R ~> kg m-3]
  real, dimension(:,:,:), allocatable, target :: temp_z ! Input temperatures [C ~> degC]
  real, dimension(:,:,:), allocatable, target :: salt_z ! Input salinities [S ~> ppt]
  real, dimension(:,:,:), allocatable, target :: mask_z ! 1 for valid data points [nondim]
  real, dimension(:,:,:), allocatable :: rho_z  ! Densities in Z-space [R ~> kg m-3]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: zi   ! Interface heights [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: dz     ! Layer thicknesses in height units [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)) :: Z_bottom  ! The (usually negative) height of the seafloor
                                                ! relative to the surface [Z ~> m].
  integer, dimension(SZI_(G),SZJ_(G))  :: nlevs ! The number of levels in each column with valid data
  real, dimension(SZI_(G))   :: press  ! Pressures [R L2 T-2 ~> Pa].

  ! Local variables for ALE remapping
  real, dimension(:), allocatable :: hTarget ! Target thicknesses [Z ~> m].
  real, dimension(:,:,:), allocatable, target :: tmpT1dIn ! Input temperatures on a model-sized grid [C ~> degC]
  real, dimension(:,:,:), allocatable, target :: tmpS1dIn ! Input salinities on a model-sized grid [S ~> ppt]
  real, dimension(:,:,:), allocatable :: tmp_mask_in      ! The valid data mask on a model-sized grid [nondim]
  real, dimension(:,:,:), allocatable :: dz1 ! Input grid thicknesses in depth units [Z ~> m]
  real, dimension(:,:,:), allocatable :: h1  ! Thicknesses on the input grid [H ~> m or kg m-2].
  real, dimension(:,:,:), allocatable :: dz_interface ! Change in position of interface due to
                                    ! regridding [H ~> m or kg m-2]
  real :: dz_neglect                ! A negligibly small vertical layer extent used in
                                    ! remapping cell reconstructions [Z ~> m]
  real :: dz_neglect_edge           ! A negligibly small vertical layer extent used in
                                    ! remapping edge value calculations [Z ~> m]
  real :: zTopOfCell, zBottomOfCell ! Heights in Z units [Z ~> m].
  type(regridding_CS) :: regridCS ! Regridding parameters and work arrays
  type(remapping_CS) :: remapCS ! Remapping parameters and work arrays

  logical :: homogenize, useALEremapping, remap_full_column, remap_general, remap_old_alg
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  integer :: remap_answer_date    ! The vintage of the order of arithmetic and expressions to use
                                  ! for remapping.  Values below 20190101 recover the remapping
                                  ! answers from 2018, while higher values use more robust
                                  ! forms of the same remapping expressions.
  integer :: hor_regrid_answer_date  ! The vintage of the order of arithmetic and expressions to use
                                  ! for horizontal regridding.  Values below 20190101 recover the
                                  ! answers from 2018, while higher values use expressions that have
                                  ! been rearranged for rotational invariance.
  logical :: pre_gridded
  logical :: separate_mixed_layer  ! If true, handle the mixed layers differently.
  logical :: density_extrap_bug    ! If true use an expression with a vertical indexing bug for
                                   ! extrapolating the densities at the bottom of unstable profiles
                                   ! from data when finding the initial interface locations in
                                   ! layered mode from a dataset of T and S.
  character(len=64) :: remappingScheme
  real :: tempAvg  ! Spatially averaged temperatures on a layer [C ~> degC]
  real :: saltAvg  ! Spatially averaged salinities on a layer [S ~> ppt]
  logical :: om4_remap_via_sub_cells ! If true, use the OM4 remapping algorithm (only used if useALEremapping)
  logical :: do_conv_adj, ignore
  integer :: nPoints
  integer :: id_clock_routine, id_clock_ALE

  id_clock_routine = cpu_clock_id('(Initialize from Z)', grain=CLOCK_ROUTINE)
  id_clock_ALE = cpu_clock_id('(Initialize from Z) ALE', grain=CLOCK_LOOP)

  call cpu_clock_begin(id_clock_routine)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  PI_180=atan(1.0)/45.

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")
  if (.not.just_read) call log_version(PF, mdl, version, "")

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  eos => tv%eqn_of_state

  call get_param(PF, mdl, "TEMP_SALT_Z_INIT_FILE", filename, &
                 "The name of the z-space input file used to initialize "//&
                 "temperatures (T) and salinities (S). If T and S are not "//&
                 "in the same file, TEMP_Z_INIT_FILE and SALT_Z_INIT_FILE "//&
                 "must be set.", default="temp_salt_z.nc", do_not_log=just_read)
  call get_param(PF, mdl, "TEMP_Z_INIT_FILE", tfilename, &
                 "The name of the z-space input file used to initialize "//&
                 "temperatures, only.", default=trim(filename), do_not_log=just_read)
  call get_param(PF, mdl, "SALT_Z_INIT_FILE", sfilename, &
                 "The name of the z-space input file used to initialize "//&
                 "temperatures, only.", default=trim(filename), do_not_log=just_read)
  filename = trim(inputdir)//trim(filename)
  tfilename = trim(inputdir)//trim(tfilename)
  sfilename = trim(inputdir)//trim(sfilename)
  call get_param(PF, mdl, "Z_INIT_FILE_PTEMP_VAR", potemp_var, &
                 "The name of the potential temperature variable in "//&
                 "TEMP_Z_INIT_FILE.", default="ptemp", do_not_log=just_read)
  call get_param(PF, mdl, "Z_INIT_FILE_SALT_VAR", salin_var, &
                 "The name of the salinity variable in "//&
                 "SALT_Z_INIT_FILE.", default="salt", do_not_log=just_read)
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
                 "If true, allows initialization directly to general coordinates.", &
                 default=.not.(GV%Boussinesq.or.GV%semi_Boussinesq) , do_not_log=just_read)
  call get_param(PF, mdl, "Z_INIT_REMAP_FULL_COLUMN", remap_full_column, &
                 "If false, only reconstructs profiles for valid data points. "//&
                 "If true, inserts vanished layers below the valid data.", &
                 default=remap_general, do_not_log=just_read)
  call get_param(PF, mdl, "Z_INIT_REMAP_OLD_ALG", remap_old_alg, &
                 "If false, uses the preferred remapping algorithm for initialization. "//&
                 "If true, use an older, less robust algorithm for remapping.", &
                 default=.false., do_not_log=just_read)
  call get_param(PF, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231, do_not_log=just_read)
  call get_param(PF, mdl, "TEMP_SALT_INIT_VERTICAL_REMAP_ONLY", pre_gridded, &
                 "If true, initial conditions are on the model horizontal grid. " //&
                 "Extrapolation over missing ocean values is done using an ICE-9 "//&
                 "procedure with vertical ALE remapping .", &
                 default=.false., do_not_log=just_read)
  if (useALEremapping) then
    call get_param(PF, mdl, "REMAPPING_ANSWER_DATE", remap_answer_date, &
                 "The vintage of the expressions and order of arithmetic to use for remapping.  "//&
                 "Values below 20190101 result in the use of older, less accurate expressions "//&
                 "that were in use at the end of 2018.  Higher values result in the use of more "//&
                 "robust and accurate forms of mathematically equivalent expressions.", &
                 default=default_answer_date, do_not_log=just_read.or.(.not.GV%Boussinesq))
    call get_param(PF, mdl, "REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                   do_not_log=.true., default=.true.)
    call get_param(PF, mdl, "Z_INIT_REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 "If true, use the OM4 remapping-via-subcells algorithm for initialization. "//&
                 "See REMAPPING_USE_OM4_SUBCELLS for more details. "//&
                 "We recommend setting this option to false.", default=om4_remap_via_sub_cells)
    if (.not.GV%Boussinesq) remap_answer_date = max(remap_answer_date, 20230701)
  endif
  call get_param(PF, mdl, "HOR_REGRID_ANSWER_DATE", hor_regrid_answer_date, &
                 "The vintage of the order of arithmetic for horizontal regridding.  "//&
                 "Dates before 20190101 give the same answers as the code did in late 2018, "//&
                 "while later versions add parentheses for rotational symmetry.  "//&
                 "Dates after 20230101 use reproducing sums for global averages.", &
                 default=default_answer_date, do_not_log=just_read.or.(.not.GV%Boussinesq))
  if (.not.GV%Boussinesq) hor_regrid_answer_date = max(hor_regrid_answer_date, 20230701)

  if (.not.useALEremapping) then
    call get_param(PF, mdl, "ADJUST_THICKNESS", correct_thickness, &
                 "If true, all mass below the bottom removed if the "//&
                 "topography is shallower than the thickness input file "//&
                 "would indicate.", default=.false., do_not_log=just_read)
    call get_param(PF, mdl, "THICKNESS_TOLERANCE", h_tolerance, &
                 "A parameter that controls the tolerance when adjusting the "//&
                 "thickness to fit the bathymetry. Used when ADJUST_THICKNESS=True.", &
                 units="m", default=0.1, scale=US%m_to_Z, &
                 do_not_log=(just_read.or..not.correct_thickness))
    call get_param(PF, mdl, "DZ_BOTTOM_TOLERANCE", tol_dz_bot, &
                 "A tolerance for detecting inconsistent topography and input layer "//&
                 "thicknesses when ADJUST_THICKNESS is false.", &
                 units="m", default=1.0, scale=US%m_to_Z, &
                 do_not_log=(just_read.or.correct_thickness))

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
    call get_param(PF, mdl, "MINIMUM_DEPTH", Hmix_default, &
                 units="m", default=0.0, scale=US%m_to_Z)
    call get_param(PF, mdl, "Z_INIT_HMIX_DEPTH", Hmix_depth, &
                 "The mixed layer depth in the initial conditions when Z_INIT_SEPARATE_MIXED_LAYER "//&
                 "is set to true.", units="m", default=US%Z_to_m*Hmix_default, scale=US%m_to_Z, &
                 do_not_log=(just_read .or. .not.separate_mixed_layer))
    ! Reusing MINIMUM_DEPTH for the default mixed layer depth may be a strange choice, but
    ! it reproduces previous answers.
    call get_param(PF, mdl, "DENSITY_INTERP_TOLERANCE", eps_rho, &
                 "A small density tolerance used when finding depths in a density profile.", &
                 units="kg m-3", default=1.0e-10, scale=US%kg_m3_to_R, &
                 do_not_log=useALEremapping.or.just_read)
    call get_param(PF, mdl, "LAYER_Z_INIT_IC_EXTRAP_BUG", density_extrap_bug, &
                 "If true use an expression with a vertical indexing bug for extrapolating the "//&
                 "densities at the bottom of unstable profiles from data when finding the "//&
                 "initial interface locations in layered mode from a dataset of T and S.", &
                 default=.false., do_not_log=just_read)
  endif
  call get_param(PF, mdl, "LAND_FILL_TEMP", temp_land_fill, &
                 "A value to use to fill in ocean temperatures on land points.", &
                 units="degC", default=0.0, scale=US%degC_to_C, do_not_log=just_read)
  call get_param(PF, mdl, "LAND_FILL_SALIN", salt_land_fill, &
                 "A value to use to fill in ocean salinities on land points.", &
                 units="ppt", default=35.0, scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(PF, mdl, "HORIZ_INTERP_TOL_TEMP", tol_temp, &
                 "The tolerance in temperature changes between iterations when interpolating "//&
                 "from an input dataset using horiz_interp_and_extrap_tracer.  This routine "//&
                 "converges slowly, so an overly small tolerance can get expensive.", &
                 units="degC", default=1.0e-3, scale=US%degC_to_C, do_not_log=just_read)
  call get_param(PF, mdl, "HORIZ_INTERP_TOL_SALIN", tol_sal, &
                 "The tolerance in salinity changes between iterations when interpolating "//&
                 "from an input dataset using horiz_interp_and_extrap_tracer.  This routine "//&
                 "converges slowly, so an overly small tolerance can get expensive.", &
                 units="ppt", default=1.0e-3, scale=US%ppt_to_S, do_not_log=just_read)

  if (just_read) then
    if ((.not.useALEremapping) .and. adjust_temperature) &
      ! This call is just here to read and log the determine_temperature parameters
      call determine_temperature(tv%T, tv%S, GV%Rlay(1:nz), eos, tv%P_Ref, 0, &
                                 0, G, GV, US, PF, just_read=.true.)
    call cpu_clock_end(id_clock_routine)
    return ! All run-time parameters have been read, so return.
  endif

  eps_z = GV%Angstrom_Z

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

  call horiz_interp_and_extrap_tracer(tfilename, potemp_var, 1, &
            G, temp_z, mask_z, z_in, z_edges_in, missing_value_temp, &
            scale=US%degC_to_C, homogenize=homogenize, m_to_Z=US%m_to_Z, &
            answer_date=hor_regrid_answer_date, ongrid=pre_gridded, tr_iter_tol=tol_temp)

  call horiz_interp_and_extrap_tracer(sfilename, salin_var, 1, &
            G, salt_z, mask_z, z_in, z_edges_in, missing_value_salt, &
            scale=US%ppt_to_S, homogenize=homogenize, m_to_Z=US%m_to_Z, &
            answer_date=hor_regrid_answer_date, ongrid=pre_gridded, tr_iter_tol=tol_sal)

  kd = size(z_in,1)

  ! Convert the sign convention of Z_edges_in.
  do k=1,size(Z_edges_in,1) ; Z_edges_in(k) = -Z_edges_in(k) ; enddo

  ! Convert T&S to Absolute Salinity and Conservative Temperature if using TEOS10 or NEMO
  call convert_temp_salt_for_TEOS10(temp_z, salt_z, G%HI, kd, mask_z, eos)

  do j=js,je ; do i=is,ie
    Z_bottom(i,j) = -depth_tot(i,j)
  enddo ; enddo

  ! Done with horizontal interpolation.
  ! Now remap to model coordinates
  if (useALEremapping) then
    call cpu_clock_begin(id_clock_ALE)
    nkd = max(GV%ke, kd)

    ! Build the source grid and copy data onto model-shaped arrays with vanished layers
    allocate( tmp_mask_in(isd:ied,jsd:jed,nkd), source=0.0 )
    allocate( dz1(isd:ied,jsd:jed,nkd), source=0.0 )
    allocate( h1(isd:ied,jsd:jed,nkd), source=0.0 )
    allocate( tmpT1dIn(isd:ied,jsd:jed,nkd), source=0.0 )
    allocate( tmpS1dIn(isd:ied,jsd:jed,nkd), source=0.0 )
    do j = js, je ; do i = is, ie
      if (G%mask2dT(i,j) > 0.) then
        zTopOfCell = 0. ; zBottomOfCell = 0.
        tmp_mask_in(i,j,1:kd) = mask_z(i,j,:)
        do k = 1, nkd
          if ((tmp_mask_in(i,j,k) > 0.) .and. (k <= kd)) then
            zBottomOfCell = max( z_edges_in(k+1), Z_bottom(i,j))
            tmpT1dIn(i,j,k) = temp_z(i,j,k)
            tmpS1dIn(i,j,k) = salt_z(i,j,k)
          elseif (k>1) then
            zBottomOfCell = Z_bottom(i,j)
            tmpT1dIn(i,j,k) = tmpT1dIn(i,j,k-1)
            tmpS1dIn(i,j,k) = tmpS1dIn(i,j,k-1)
          else ! This next block should only ever be reached over land
            tmpT1dIn(i,j,k) = temp_land_fill
            tmpS1dIn(i,j,k) = salt_land_fill
          endif
          dz1(i,j,k) = (zTopOfCell - zBottomOfCell)
          zTopOfCell = zBottomOfCell ! Bottom becomes top for next value of k
        enddo
        dz1(i,j,kd) = dz1(i,j,kd) + max(0., zTopOfCell - Z_bottom(i,j) )
        ! The max here is in case the data data is shallower than model
      endif ! mask2dT
    enddo ; enddo
    deallocate( tmp_mask_in )

    ! Convert input thicknesses to units of H.  In non-Boussinesq mode this is done by inverting
    ! integrals of specific volume in pressure, so it can be expensive.
    tv_loc = tv
    tv_loc%T => tmpT1dIn
    tv_loc%S => tmpS1dIn
    GV_loc = GV
    GV_loc%ke = nkd
    call dz_to_thickness(dz1, tv_loc, h1, G, GV_loc, US)

    ! Build the target grid (and set the model thickness to it)

    call ALE_initRegridding( GV, US, G%max_depth, PF, mdl, regridCS ) ! sets regridCS
    if (remap_general) then
      dz_neglect = set_h_neglect(GV, remap_answer_date, dz_neglect_edge)
    else
      dz_neglect = set_dz_neglect(GV, US, remap_answer_date, dz_neglect_edge)
    endif
    call initialize_remapping( remapCS, remappingScheme, boundary_extrapolation=.false., &
                               om4_remap_via_sub_cells=om4_remap_via_sub_cells, answer_date=remap_answer_date, &
                               h_neglect=dz_neglect, h_neglect_edge=dz_neglect_edge)

    ! Now remap from source grid to target grid, first setting reconstruction parameters
    if (remap_general) then
      call set_regrid_params( regridCS, min_thickness=0. )
      allocate( dz_interface(isd:ied,jsd:jed,nkd+1) ) ! Need for argument to regridding_main() but is not used

      call regridding_preadjust_reqs(regridCS, do_conv_adj, ignore)
      if (do_conv_adj) call convective_adjustment(G, GV_loc, h1, tv_loc)
      call regridding_main( remapCS, regridCS, G, GV_loc, US, h1, tv_loc, h, dz_interface, &
                            frac_shelf_h=frac_shelf_h )

      deallocate( dz_interface )

      call ALE_remap_scalar(remapCS, G, GV, nkd, h1, tmpT1dIn, h, tv%T, all_cells=remap_full_column, &
                            old_remap=remap_old_alg )
      call ALE_remap_scalar(remapCS, G, GV, nkd, h1, tmpS1dIn, h, tv%S, all_cells=remap_full_column, &
                            old_remap=remap_old_alg )
    else
      ! This is the old way of initializing to z* coordinates only
      allocate( hTarget(nz) )
      hTarget = getCoordinateResolution( regridCS )
      do j = js, je ; do i = is, ie
        dz(i,j,:) = 0.
        if (G%mask2dT(i,j) > 0.) then
          ! Build the target grid combining hTarget and topography
          zTopOfCell = 0. ; zBottomOfCell = 0.
          do k = 1, nz
            zBottomOfCell = max( zTopOfCell - hTarget(k), Z_bottom(i,j))
            dz(i,j,k) = zTopOfCell - zBottomOfCell
            zTopOfCell = zBottomOfCell ! Bottom becomes top for next value of k
          enddo
        else
          dz(i,j,:) = 0.
        endif ! mask2dT
      enddo ; enddo
      deallocate( hTarget )

      dz_neglect = set_dz_neglect(GV, US, remap_answer_date, dz_neglect_edge)
      call ALE_remap_scalar(remapCS, G, GV, nkd, dz1, tmpT1dIn, dz, tv%T, all_cells=remap_full_column, &
                            old_remap=remap_old_alg)
      call ALE_remap_scalar(remapCS, G, GV, nkd, dz1, tmpS1dIn, dz, tv%S, all_cells=remap_full_column, &
                            old_remap=remap_old_alg)

      if (GV%Boussinesq .or. GV%semi_Boussinesq) then
        ! This is a simple conversion of the target grid to thickness units that is not
        ! appropriate in non-Boussinesq mode.
        call dz_to_thickness_simple(dz, h, G, GV, US)
      else
        ! Convert dz into thicknesses in units of H using the equation of state as appropriate.
        call dz_to_thickness(dz, tv, h, G, GV, US)
      endif
    endif

    deallocate( dz1 )
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
    Rb(1) = 0.0
    if (nz>1) then
      Rb(nz+1) = 2.0*GV%Rlay(nz) - GV%Rlay(nz-1)
    else
      Rb(nz+1) = 2.0 * GV%Rlay(1)
    endif

    nkml = 0 ; if (separate_mixed_layer) nkml = GV%nkml

    press(:) = tv%P_Ref
    EOSdom(:) = EOS_domain(G%HI)
    allocate(rho_z(isd:ied,jsd:jed,kd))
    do k=1,kd ; do j=js,je
      call calculate_density(temp_z(:,j,k), salt_z(:,j,k), press, rho_z(:,j,k), eos, EOSdom)
    enddo ; enddo

    call find_interfaces(rho_z, z_in, kd, Rb, Z_bottom, zi, G, GV, US, nlevs, nkml, &
                         Hmix_depth, eps_z, eps_rho, density_extrap_bug)

    deallocate(rho_z)

    dz(:,:,:) = 0.0
    if (correct_thickness) then
      call adjustEtaToFitBathymetry(G, GV, US, zi, dz, h_tolerance, dZ_ref_eta=G%Z_ref)
    else
      do k=nz,1,-1 ; do j=js,je ; do i=is,ie
        if (zi(i,j,K) < (zi(i,j,K+1) + GV%Angstrom_Z)) then
          zi(i,j,K) = zi(i,j,K+1) + GV%Angstrom_Z
          dz(i,j,k) = GV%Angstrom_Z
        else
          dz(i,j,k) = zi(i,j,K) - zi(i,j,K+1)
        endif
      enddo ; enddo ; enddo
      inconsistent = 0
      do j=js,je ; do i=is,ie
        if (abs(zi(i,j,nz+1) - Z_bottom(i,j)) > tol_dz_bot) &
          inconsistent = inconsistent + 1
      enddo ; enddo
      call sum_across_PEs(inconsistent)

      if ((inconsistent > 0) .and. (is_root_pe())) then
        write(mesg, '("Thickness initial conditions are inconsistent ",'// &
                    '"with topography in ",I5," places.")') inconsistent
        call MOM_error(WARNING, mesg)
      endif
    endif

    call tracer_z_init_array(temp_z, z_edges_in, kd, zi, temp_land_fill, G, nz, nlevs, eps_z, tv%T)
    call tracer_z_init_array(salt_z, z_edges_in, kd, zi, salt_land_fill, G, nz, nlevs, eps_z, tv%S)

    if (homogenize) then
      ! Horizontally homogenize data to produce perfectly "flat" initial conditions
      do k=1,nz
        call homogenize_field(tv%T(:,:,k), G, tmp_scale=US%C_to_degC, answer_date=hor_regrid_answer_date)
        call homogenize_field(tv%S(:,:,k), G, tmp_scale=US%S_to_ppt, answer_date=hor_regrid_answer_date)
      enddo
    endif

    if (adjust_temperature) then
      ! Finally adjust to target density
      ks = 1 ; if (separate_mixed_layer) ks = GV%nk_rho_varies + 1
      call determine_temperature(tv%T, tv%S, GV%Rlay(1:nz), eos, tv%P_Ref, niter, &
                                 ks, G, GV, US, PF, just_read)
    endif

    ! Now convert dz into thicknesses in units of H.
    call dz_to_thickness(dz, tv, h, G, GV, US)

  endif ! useALEremapping

  deallocate(z_in, z_edges_in, temp_z, salt_z, mask_z)

  call pass_var(h, G%Domain)
  call pass_var(tv%T, G%Domain)
  call pass_var(tv%S, G%Domain)

  call callTree_leave(trim(mdl)//'()')
  call cpu_clock_end(id_clock_routine)

end subroutine MOM_temp_salt_initialize_from_Z


!> Find interface positions corresponding to interpolated depths in a density profile
subroutine find_interfaces(rho, zin, nk_data, Rb, Z_bot, zi, G, GV, US, nlevs, nkml, hml, &
                           eps_z, eps_rho, density_extrap_bug)
  type(ocean_grid_type),      intent(in)  :: G     !< The ocean's grid structure
  type(verticalGrid_type),    intent(in)  :: GV    !< The ocean's vertical grid structure
  integer,                    intent(in)  :: nk_data !< The number of levels in the input data
  real, dimension(SZI_(G),SZJ_(G),nk_data), &
                              intent(in)  :: rho   !< Potential density in z-space [R ~> kg m-3]
  real, dimension(nk_data),   intent(in)  :: zin   !< Input data levels [Z ~> m].
  real, dimension(SZK_(GV)+1), intent(in) :: Rb    !< target interface densities [R ~> kg m-3]
  real, dimension(SZI_(G),SZJ_(G)), &
                              intent(in)  :: Z_bot !< The (usually negative) height of the seafloor
                                                   !! relative to the surface [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                              intent(out) :: zi    !< The returned interface heights [Z ~> m]
  type(unit_scale_type),      intent(in)  :: US    !< A dimensional unit scaling type
  integer, dimension(SZI_(G),SZJ_(G)), &
                              intent(in)  :: nlevs !< number of valid points in each column
  integer,                    intent(in)  :: nkml  !< number of mixed layer pieces to distribute over
                                                   !! a depth of hml.
  real,                       intent(in)  :: hml   !< mixed layer depth [Z ~> m].
  real,                       intent(in)  :: eps_z !< A negligibly small layer thickness [Z ~> m].
  real,                       intent(in)  :: eps_rho !< A negligibly small density difference [R ~> kg m-3].
  logical,                    intent(in)  :: density_extrap_bug !< If true use an expression with an
                                                   !! indexing bug for projecting the densities at
                                                   !! the bottom of unstable profiles from data when
                                                   !! finding the initial interface locations in
                                                   !! layered mode from a dataset of T and S.

  ! Local variables
  real, dimension(nk_data) :: rho_ ! A column of densities [R ~> kg m-3]
  real, dimension(SZK_(GV)+1) :: zi_ ! A column interface heights (negative downward) [Z ~> m].
  real    :: slope      ! The rate of change of height with density [Z R-1 ~> m4 kg-1]
  real    :: drhodz     ! A local vertical density gradient [R Z-1 ~> kg m-4]
  real, parameter :: zoff = 0.999 ! A small fractional adjustment to the density differences [nondim]
  logical :: unstable   ! True if the column is statically unstable anywhere.
  integer :: nlevs_data ! The number of data values in a column.
  logical :: work_down  ! This indicates whether this pass goes up or down the water column.
  integer :: k_int, lo_int, hi_int, mid
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  zi(:,:,:) = 0.0

  do j=js,je ; do i=is,ie
    nlevs_data = nlevs(i,j)
    do k=1,nlevs_data ; rho_(k) = rho(i,j,k) ; enddo

    unstable = .true.
    work_down = .true.
    do while (unstable)
      ! Modify the input profile until it no longer has densities that decrease with depth.
      unstable = .false.
      if (work_down) then
        do k=2,nlevs_data-1 ; if (rho_(k) - rho_(k-1) < 0.0) then
          if (k == 2) then
            rho_(k-1) = rho_(k) - eps_rho
          else
            drhodz = (rho_(k+1)-rho_(k-1)) / (zin(k+1)-zin(k-1))
            if (drhodz < 0.0) unstable = .true.
            rho_(k) = rho_(k-1) + drhodz*zoff*(zin(k)-zin(k-1))
          endif
        endif ; enddo
        work_down = .false.
      else
        do k=nlevs_data-1,2,-1 ;  if (rho_(k+1) - rho_(k) < 0.0) then
          if (k == nlevs_data-1) then
            if (density_extrap_bug) then
              rho_(k+1) = rho_(k-1) + eps_rho
            else
              rho_(k+1) = rho_(k) + eps_rho
            endif
          else
            drhodz = (rho_(k+1)-rho_(k-1)) / (zin(k+1)-zin(k-1))
            if (drhodz < 0.0) unstable = .true.
            rho_(k) = rho_(k+1) - drhodz*(zin(k+1)-zin(k))
          endif
        endif ; enddo
        work_down = .true.
      endif
    enddo

    ! Find and store the interface depths.
    zi_(1) = 0.0
    if (nlevs_data < 1) then
      ! There is no data to use, so set the interfaces at the bottom.
      do K=2,nz ; zi_(K) = Z_bot(i,j) ; enddo
    elseif (nlevs_data == 1) then
      ! There is data for only one input layer, so set the interfaces at the bottom or top,
      ! depending on how their target densities compare with the one data point.
      do K=2,nz
        if (Rb(K) < rho_(1)) then ; zi_(K) = 0.0
        else ; zi_(K) = Z_bot(i,j) ; endif
      enddo
    else
      do K=2,nz
        ! Find the value of k_int in the list of rho_ where rho_(k_int) <= Rb(K) < rho_(k_int+1).
        ! This might be made a little faster by exploiting the fact that Rb is
        ! monotonically increasing and not resetting lo_int back to 1 inside the K loop.
        lo_int = 1 ; hi_int = nlevs_data
        do while (lo_int < hi_int)
          mid = (lo_int+hi_int) / 2
          if (Rb(K) < rho_(mid)) then ; hi_int = mid
          else ; lo_int = mid+1 ; endif
        enddo
        k_int = max(1, lo_int-1)

        ! Linearly interpolate to find the depth, zi_, where Rb would be found.
        slope = (zin(k_int+1) - zin(k_int)) / max(rho_(k_int+1) - rho_(k_int), eps_rho)
        zi_(K) = -1.0*(zin(k_int) + slope*(Rb(K)-rho_(k_int)))
        zi_(K) = min(max(zi_(K), Z_bot(i,j)), -1.0*hml)
      enddo
    endif
    zi_(nz+1) = Z_bot(i,j)
    if (nkml > 0) then ; do K=2,nkml+1
      zi_(K) = max(hml*((1.0-real(K))/real(nkml)), Z_bot(i,j))
    enddo ; endif
    do K=nz,max(nkml+2,2),-1
      if (zi_(K) < zi_(K+1) + eps_Z) zi_(K) = zi_(K+1) + eps_Z
      if (zi_(K) > -1.0*hml)  zi_(K) = max(-1.0*hml, Z_bot(i,j))
    enddo

    do K=1,nz+1
      zi(i,j,K) = zi_(K)
    enddo
  enddo ; enddo ! i- and j- loops

end subroutine find_interfaces

!> Run simple unit tests
subroutine MOM_state_init_tests(G, GV, US, tv)
  type(ocean_grid_type),     intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),     intent(in)    :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),     intent(in)    :: tv   !< Thermodynamics structure.

  ! Local variables
  integer, parameter :: nk=5
  real, dimension(nk) :: T, T_t, T_b ! Temperatures [C ~> degC]
  real, dimension(nk) :: S, S_t, S_b ! Salinities [S ~> ppt]
  real, dimension(nk) :: rho ! Layer density [R ~> kg m-3]
  real, dimension(nk) :: h   ! Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nk) :: z   ! Height of layer center [Z ~> m]
  real, dimension(nk+1) :: e ! Interface heights [Z ~> m]
  real :: T_ref              ! A reference temperature [C ~> degC]
  real :: S_ref              ! A reference salinity [S ~> ppt]
  real :: P_tot, P_t, P_b    ! Pressures [R L2 T-2 ~> Pa]
  real :: z_out              ! Output height [Z ~> m]
  real :: I_z_scale          ! The inverse of the height scale for prescribed gradients [Z-1 ~> m-1]
  real :: z_tol              ! The tolerance with which to find the depth matching a specified pressure [Z ~> m].
  integer :: k
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
  T_ref = 20.0*US%degC_to_C
  S_ref = 35.0*US%ppt_to_S
  z_tol = 1.0e-5*US%m_to_Z
  do k = 1, nk
    z(k) = 0.5 * ( e(K) + e(K+1) )
    T_t(k) = T_ref + (0. * I_z_scale) * e(k)
    T(k)   = T_ref + (0. * I_z_scale)*z(k)
    T_b(k) = T_ref + (0. * I_z_scale)*e(k+1)
    S_t(k) = S_ref - (0. * I_z_scale)*e(k)
    S(k)   = S_ref + (0. * I_z_scale)*z(k)
    S_b(k) = S_ref - (0. * I_z_scale)*e(k+1)
    call calculate_density(0.5*(T_t(k)+T_b(k)), 0.5*(S_t(k)+S_b(k)), -GV%Rho0*GV%g_Earth*z(k), &
                           rho(k), tv%eqn_of_state)
    P_tot = P_tot + GV%g_Earth * rho(k) * GV%H_to_Z*h(k)
  enddo

  P_t = 0.
  do k = 1, nk
    call find_depth_of_pressure_in_cell(T_t(k), T_b(k), S_t(k), S_b(k), e(K), e(K+1), P_t, 0.5*P_tot, &
                                        GV%Rho0, GV%g_Earth, tv%eqn_of_state, US, P_b, z_out, z_tol=z_tol, &
                                        frac_dp_bugfix=.false.)
    write(0,*) k, US%RL2_T2_to_Pa*P_t, US%RL2_T2_to_Pa*P_b, 0.5*US%RL2_T2_to_Pa*P_tot, &
               US%Z_to_m*e(K), US%Z_to_m*e(K+1), US%Z_to_m*z_out
    P_t = P_b
  enddo
  write(0,*) US%RL2_T2_to_Pa*P_b, US%RL2_T2_to_Pa*P_tot

  write(0,*) ''
  write(0,*) ' ==================================================================== '
  write(0,*) ''
  write(0,*) GV%H_to_m*h(:)

  ! For consistency with the usual call, add the following:
  ! if (use_remapping) then
  !   allocate(remap_CS)
  !   call initialize_remapping(remap_CS, 'PLM', boundary_extrapolation=.true., &
  !                             h_neglect=GV%H_subroundoff, h_neglect_edge=GV%H_subroundoff)
  ! endif
  call cut_off_column_top(nk, tv, GV, US, GV%g_Earth, -e(nk+1), GV%Angstrom_H, &
                          T, T_t, T_b, S, S_t, S_b, 0.5*P_tot, h, remap_CS, z_tol=z_tol, &
                          frac_dp_bugfix=.false.)
  write(0,*) GV%H_to_m*h(:)
  if (associated(remap_CS)) deallocate(remap_CS)

end subroutine MOM_state_init_tests

end module MOM_state_initialization
