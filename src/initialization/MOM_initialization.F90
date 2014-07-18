module MOM_initialization
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

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, April 1994 - June 2002                         *
!*                                                                     *
!*    This subroutine initializes the fields for the simulations.      *
!*  The one argument passed to initialize, Time, is set to the         *
!*  current time of the simulation.  The fields which are initialized  *
!*  here are:                                                          *
!*    u - Zonal velocity in m s-1.                                     *
!*    v - Meridional velocity in m s-1.                                *
!*    h - Layer thickness in m.  (Must be positive.)                   *
!*    G%bathyT - Basin depth in m.  (Must be positive.)                *
!*    G%CoriolisBu - The Coriolis parameter, in s-1.                   *
!*    G%g_prime - The reduced gravity at each interface, in m s-2.     *
!*    G%Rlay - Layer potential density (coordinate variable) in kg m-3.*
!*  If ENABLE_THERMODYNAMICS is defined:                               *
!*    tv%T - Temperature in C.                                         *
!*    tv%S - Salinity in psu.                                          *
!*  If SPONGE is defined:                                              *
!*    A series of subroutine calls are made to set up the damping      *
!*    rates and reference profiles for all variables that are damped   *
!*    in the sponge.                                                   *
!*  Any user provided tracer code is also first linked through this    *
!*  subroutine.                                                        *
!*                                                                     *
!*    Forcing-related fields (taux, tauy, buoy, ustar, etc.) are set   *
!*  in MOM_surface_forcing.F90.                                        *
!*                                                                     *
!*    These variables are all set in the set of subroutines (in this   *
!*  file): initialize_topography, initialize_thickness,                *
!*  initialize_velocity, initialize_temp_sal, initialize_sponges, and  *
!*  set_coordinate.                                                    *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, CoriolisBu                            *
!*    j+1  > o > o >   At ^:  v, tauy                                  *
!*    j    x ^ x ^ x   At >:  u, taux                                  *
!*    j    > o > o >   At o:  h, bathyT, buoy, tr, T, S, Rml, ustar    *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**


use MOM_checksums, only : hchksum, qchksum, uchksum, vchksum, chksum
use MOM_coms, only : max_across_PEs, min_across_PEs
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only :  CLOCK_ROUTINE, CLOCK_LOOP
use MOM_domains, only : pass_var, pass_vector, sum_across_PEs, broadcast
use MOM_domains, only : root_PE, To_All, SCALAR_PAIR, CGRID_NE
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, param_file_type
use MOM_file_parser, only : log_version
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type, isPointInCell
use MOM_interface_heights, only : find_eta
use MOM_io, only : close_file, create_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE, MULTIPLE
use MOM_io, only : slasher, vardesc, write_field
use MOM_io, only : EAST_FACE, NORTH_FACE
use MOM_grid_initialize, only : initialize_masks, set_grid_metrics
use MOM_restart, only : restore_state, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, set_up_sponge_ML_density
use MOM_sponge, only : initialize_sponge, sponge_CS
use MOM_string_functions, only : uppercase
use MOM_time_manager, only : time_type, set_time
use MOM_tracer_registry, only : add_tracer_OBC_values, tracer_registry_type
use MOM_variables, only : thermo_var_ptrs, ocean_OBC_type
use MOM_variables, only : OBC_NONE, OBC_SIMPLE, OBC_FLATHER_E, OBC_FLATHER_W
use MOM_variables, only : OBC_FLATHER_N, OBC_FLATHER_S
use MOM_verticalGrid, only : setVerticalGridAxes
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use MOM_EOS, only : int_specific_vol_dp
use user_initialization, only : user_set_coord, user_initialize_topography
use user_initialization, only : user_initialize_thickness, user_initialize_velocity
use user_initialization, only : user_init_temperature_salinity
use user_initialization, only : user_set_Open_Bdry_Conds, user_initialize_sponges
use DOME_initialization, only : DOME_initialize_thickness
use DOME_initialization, only : DOME_initialize_topography
use DOME_initialization, only : DOME_set_Open_Bdry_Conds
use DOME_initialization, only : DOME_initialize_sponges
use benchmark_initialization, only : benchmark_initialize_thickness
use benchmark_initialization, only : benchmark_initialize_topography
use benchmark_initialization, only : benchmark_init_temperature_salinity
use circle_obcs_initialization, only : circle_obcs_initialize_thickness
use lock_exchange_initialization, only : lock_exchange_initialize_thickness
use external_gwave_initialization, only : external_gwave_initialize_thickness
use DOME2d_initialization, only : DOME2d_initialize_topography
use DOME2d_initialization, only : DOME2d_initialize_thickness
use DOME2d_initialization, only : DOME2d_initialize_temperature_salinity
use adjustment_initialization, only : adjustment_initialize_thickness
use adjustment_initialization, only : adjustment_initialize_temperature_salinity
use sloshing_initialization, only : sloshing_initialize_topography
use sloshing_initialization, only : sloshing_initialize_thickness
use sloshing_initialization, only : sloshing_initialize_temperature_salinity
use seamount_initialization, only : seamount_initialize_topography
use seamount_initialization, only : seamount_initialize_thickness
use seamount_initialization, only : seamount_initialize_temperature_salinity
use Phillips_initialization, only : Phillips_initialize_thickness
use Phillips_initialization, only : Phillips_initialize_velocity
use Phillips_initialization, only : Phillips_initialize_sponges
use Rossby_front_2d_initialization, only : Rossby_front_initialize_thickness
use Rossby_front_2d_initialization, only : Rossby_front_initialize_temperature_salinity
use Rossby_front_2d_initialization, only : Rossby_front_initialize_velocity

use midas_vertmap, only : find_interfaces, tracer_Z_init, meshgrid
use midas_vertmap, only : determine_temperature

use MOM_ALE, only : ALE_initRegridding
use MOM_regridding, only : regridding_CS
use MOM_remapping, only : remapping_CS, remapping_core, initialize_remapping
use MOM_remapping, only : dzFromH1H2, remapDisableBoundaryExtrapolation

use mpp_domains_mod, only  : mpp_global_field, mpp_get_compute_domain
use mpp_mod, only          : mpp_broadcast,mpp_root_pe,mpp_sync,mpp_sync_self


use horiz_interp_mod, only : horiz_interp_new, horiz_interp,horiz_interp_type
use horiz_interp_mod, only : horiz_interp_init, horiz_interp_del

use netcdf

implicit none ; private

#include <MOM_memory.h>

public MOM_initialize, MOM_initialize_rotation, MOM_initialize_topography

! This structure is to simplify communication with the calling code.
type, public :: MOM_initialization_struct
  type(tracer_registry_type), pointer :: tracer_Reg => NULL()
  type(sponge_CS), pointer :: sponge_CSp => NULL()
  type(ocean_OBC_type), pointer :: OBC => NULL()
end type MOM_initialization_struct

contains

! -----------------------------------------------------------------------------
subroutine MOM_initialize(u, v, h, tv, Time, G, PF, dirs, &
                              restart_CS, CS, Time_in)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(out)   :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(out)   :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(out)   :: h
  type(thermo_var_ptrs),                  intent(inout) :: tv
  type(time_type),                        intent(inout) :: Time
  type(ocean_grid_type),                  intent(inout) :: G
  type(param_file_type),                  intent(in)    :: PF
  type(directories),                      intent(in)    :: dirs
  type(MOM_restart_CS),                   pointer       :: restart_CS
  type(MOM_initialization_struct),        intent(inout) :: CS
  type(time_type), optional,              intent(in)    :: Time_in
! Arguments: u  - Zonal velocity, in m s-1.
!  (out)     v  - Meridional velocity, in m s-1.
!  (out)     h  - Layer thickness, in m.
!  (out)     tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (out)     Time    - Time at the start of the run segment.
!  (inout)   G       - The ocean's grid structure.
!  (in)      PF      - A structure indicating the open file to parse for
!                      model parameter values.
!  (in)      dirs    - A structure containing several relevant directory paths.
!  (inout)   restart_CS - A pointer to the restart control structure.
!  (inout)   CS      - A structure of pointers to be exchanged with MOM.F90.
!  (in)      Time_in - Time at the start of the run segment. Time_in overrides
!                      any value set for Time.

  character(len=200) :: filename   ! The name of an input file.
  character(len=200) :: filename2  ! The name of an input files.
  character(len = 200) :: inputdir ! The directory where NetCDF input files are.
  character(len=200) :: config
  logical :: from_Z_file
  logical :: new_sim
  integer :: write_geom
  logical :: use_temperature, use_sponge
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  logical :: depress_sfc ! If true, remove the mass that would be displaced
                         ! by a large surface pressure, such as with an ice
                         ! sheet.
  logical :: Analytic_FV_PGF, obsol_test
  logical :: apply_OBC_u, apply_OBC_v
  logical :: apply_OBC_u_flather_east, apply_OBC_u_flather_west
  logical :: apply_OBC_v_flather_north, apply_OBC_v_flather_south
  logical :: convert
  type(EOS_type), pointer :: eos => NULL()
  logical :: debug    ! indicates whether to write debugging output
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_initialization" ! This module's name.
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call callTree_enter("MOM_initialize(), MOM_initialization.F90")
  call log_version(PF, mod, version)

  new_sim = .false.
  if ((dirs%input_filename(1:1) == 'n') .and. &
      (LEN_TRIM(dirs%input_filename) == 1)) new_sim = .true.

  call get_param(PF, mod, "INPUTDIR", inputdir, &
         "The directory in which input files are found.", default=".")
  inputdir = slasher(inputdir)

  use_temperature = ASSOCIATED(tv%T)
  use_EOS = associated(tv%eqn_of_state)
  if (use_EOS) eos => tv%eqn_of_state

  call get_param(PF, mod, "DEBUG", debug, default=.false.)

! ====================================================================
!    Initialize fields that are time invariant - metrics, topography,
!  masks, vertical coordinate, Coriolis parameter.
! ====================================================================

! Set-up the layer densities, G%Rlay, and reduced gravities, G%g_prime.
  call get_param(PF, mod, "COORD_CONFIG", config, &
                 "This specifies how layers are to be defined: \n"//&
                 " \t file - read coordinate information from the file \n"//&
                 " \t\t specified by (COORD_FILE).\n"//&
                 " \t linear - linear based on interfaces not layesrs. \n"//&
                 " \t ts_ref - use reference temperature and salinity \n"//&
                 " \t ts_range - use range of temperature and salinity \n"//&
                 " \t\t (T_REF and S_REF) to determine surface density \n"//&
                 " \t\t and GINT calculate internal densities. \n"//&
                 " \t gprime - use reference density (RHO_0) for surface \n"//&
                 " \t\t density and GINT calculate internal densities. \n"//&
                 " \t ts_profile - use temperature and salinity profiles \n"//&
                 " \t\t (read from COORD_FILE) to set layer densities. \n"//&
                 " \t USER - call a user modified routine.", &
                 fail_if_missing=.true.)
  select case ( trim(config) )
    case ("gprime")
      call set_coord_from_gprime(G%Rlay, G%g_prime, G, PF)
    case ("layer_ref")
      call set_coord_from_layer_density(G%Rlay, G%g_prime, G, PF)
    case ("linear")
      call set_coord_linear(G%Rlay, G%g_prime, G, PF)
    case ("ts_ref")
      call set_coord_from_ts_ref(G%Rlay, G%g_prime, G, PF, eos, tv%P_Ref)
    case ("ts_profile")
      call set_coord_from_TS_profile(G%Rlay, G%g_prime, G, PF, eos, tv%P_Ref)
    case ("ts_range")
      call set_coord_from_TS_range(G%Rlay, G%g_prime, G, PF, eos, tv%P_Ref)
    case ("file")
      call set_coord_from_file(G%Rlay, G%g_prime, G, PF)
    case ("USER")
      call user_set_coord(G%Rlay, G%g_prime, G, PF, eos)
    case ("none")
    case default ; call MOM_error(FATAL,"MOM_initialize: "// &
      "Unrecognized coordinate setup"//trim(config))
  end select
  if (debug) call chksum(G%Rlay, "MOM_initialize: Rlay ", 1, nz)
  if (debug) call chksum(G%g_prime, "MOM_initialize: g_prime ", 1, nz)
  call setVerticalGridAxes( G%Rlay, G%GV )

! Set up the parameters of the physical domain (i.e. the grid), G
  call set_grid_metrics(G, PF)

! Set up the bottom depth, G%bathyT either analytically or from file
  call MOM_initialize_topography(G%bathyT, G%max_depth, G, PF)

!    This call sets seamasks that prohibit flow over any point with  !
!  a bottom that is shallower than min_depth from PF.                !
  call initialize_masks(G, PF)
  if (debug) then
    call hchksum(G%bathyT, 'MOM_initialize: depth ', G, haloshift=1)
    call hchksum(G%mask2dT, 'MOM_initialize: mask2dT ', G)
    call uchksum(G%mask2dCu, 'MOM_initialize: mask2dCu ', G)
    call vchksum(G%mask2dCv, 'MOM_initialize: mask2dCv ', G)
    call qchksum(G%mask2dBu, 'MOM_initialize: mask2dBu ', G)
  endif

! Modulate geometric scales according to geography.
  call get_param(PF, mod, "CHANNEL_CONFIG", config, &
                 "A parameter that determines which set of channels are \n"//&
                 "restricted to specific  widths.  Options are:\n"//&
                 " \t none - All channels have the grid width.\n"//&
                 " \t global_1deg - Sets 16 specific channels appropriate \n"//&
                 " \t\t for a 1-degree model, as used in CM2G.\n"//&
                 " \t list - Read the channel locations and widths from a \n"//&
                 " \t\t text file, like MOM_channel_list in the MOM_SIS \n"//&
                 " \t\t test case.\n"//&
                 " \t file - Read open face widths everywhere from a \n"//&
                 " \t\t NetCDF file on the model grid.", &
                 default="none")
  select case ( trim(config) )
    case ("none")
    case ("list") ; call reset_face_lengths_list(G, PF)
    case ("file") ; call reset_face_lengths_file(G, PF)
    case ("global_1deg") ; call reset_face_lengths_named(G, PF, trim(config))
    case default ; call MOM_error(FATAL, "MOM_initialize: "// &
      "Unrecognized channel configuration "//trim(config))
  end select

!   This call sets the topography at velocity points.
  if (G%bathymetry_at_vel) then
    call get_param(PF, mod, "VELOCITY_DEPTH_CONFIG", config, &
                   "A string that determines how the topography is set at \n"//&
                   "velocity points. This may be 'min' or 'max'.", &
                   default="max")
    select case ( trim(config) )
      case ("max") ; call set_velocity_depth_max(G)
      case ("min") ; call set_velocity_depth_min(G)
      case default ; call MOM_error(FATAL, "MOM_initialize: "// &
        "Unrecognized velocity depth configuration "//trim(config))
    end select
  endif

!    Calculate the value of the Coriolis parameter at the latitude   !
!  of the q grid points, in s-1.
  call MOM_initialize_rotation(G%CoriolisBu, G, PF)
  if (debug) then
    call qchksum(G%CoriolisBu, "MOM_initialize: f ", G)
  endif

! Compute global integrals of grid values for later use in scalar diagnostics !
  call compute_global_grid_integrals(G)

! Write out all of the grid data used by this run.
  call get_param(PF, mod, "WRITE_GEOM", write_geom, &
                 "If =0, never write the geometry and vertical grid files.\n"//&
                 "If =1, write the geometry and vertical grid files only for\n"//&
                 "a new simulation. If =2, always write the geometry and\n"//&
                 "vertical grid files. Other values are invalid.", default=1)
  if (write_geom<0 .or. write_geom>2) call MOM_error(FATAL,"MOM_initialize: "//&
         "WRITE_GEOM must be equal to 0, 1 or 2.")
  if ((write_geom==1 .and. new_sim) .or. write_geom==2) then
    call write_ocean_geometry_file(G, PF, dirs%output_directory)
    call write_vertgrid_file(G, PF, dirs%output_directory)
  endif

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

    if (from_Z_file) then
!     Initialize thickness and T/S from z-coordinate data in a file.
      if (.NOT.use_temperature) call MOM_error(FATAL,"MOM_initialize : "//&
         "use_temperature must be true if INIT_LAYERS_FROM_Z_FILE is true")

      call MOM_temp_salt_initialize_from_Z(h, tv, G, PF, dirs)
      call pass_var(h, G%Domain)
    else
!     Initialize thickness, h.
      call get_param(PF, mod, "THICKNESS_CONFIG", config, &
               "A string that determines how the initial layer \n"//&
               "thicknesses are specified for a new run: \n"//&
               " \t file - read interface heights from the file specified \n"//&
               " \t thickness_file - read thicknesses from the file specified \n"//&
               " \t\t by (THICKNESS_FILE).\n"//&
               " \t uniform - uniform thickness layers evenly distributed \n"//&
               " \t\t between the surface and MAXIMUM_DEPTH. \n"//&
               " \t DOME - use a slope and channel configuration for the \n"//&
               " \t\t DOME sill-overflow test case. \n"//&
               " \t benchmark - use the benchmark test case thicknesses. \n"//&
               " \t search - search a density profile for the interface \n"//&
               " \t\t densities. This is not yet implemented. \n"//&
               " \t circle_obcs - the circle_obcs test case is used. \n"//&
               " \t DOME2D - 2D version of DOME initialization. \n"//&
               " \t adjustment2d - TBD AJA. \n"//&
               " \t sloshing - TBD AJA. \n"//&
               " \t seamount - TBD AJA. \n"//&
               " \t rossby_front - a mixed layer front in thermal wind balance.\n"//&
               " \t USER - call a user modified routine.", &
               fail_if_missing=.true.)
      select case (trim(config))
         case ("file"); call initialize_thickness_from_file(h, G, PF, .false.)
         case ("thickness_file"); call initialize_thickness_from_file(h, G, PF, .true.)
         case ("uniform"); call initialize_thickness_uniform(h, G, PF)
         case ("DOME"); call DOME_initialize_thickness(h, G, PF)
         case ("benchmark"); call benchmark_initialize_thickness(h, G, PF, &
                                 tv%eqn_of_state, tv%P_Ref)
         case ("search"); call initialize_thickness_search
         case ("circle_obcs"); call circle_obcs_initialize_thickness(h, G, PF)
         case ("lock_exchange"); call lock_exchange_initialize_thickness(h, G, PF)
         case ("external_gwave"); call external_gwave_initialize_thickness(h, G, PF)
         case ("DOME2D"); call DOME2d_initialize_thickness(h, G, PF)
         case ("adjustment2d"); call adjustment_initialize_thickness(h, G, PF)
         case ("sloshing"); call sloshing_initialize_thickness(h, G, PF)
         case ("seamount"); call seamount_initialize_thickness(h, G, PF)
         case ("phillips"); call Phillips_initialize_thickness(h, G, PF)
         case ("rossby_front"); call Rossby_front_initialize_thickness(h, G, PF)
         case ("USER"); call user_initialize_thickness(h, G, PF, tv%T)
         case default ; call MOM_error(FATAL,  "MOM_initialize: "//&
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
               " \t adjustment2d - TBD AJA. \n"//&
               " \t sloshing - TBD AJA. \n"//&
               " \t seamount - TBD AJA. \n"//&
               " \t rossby_front - a mixed layer front in thermal wind balance.\n"//&
               " \t USER - call a user modified routine.", &
               fail_if_missing=.true.)
        select case (trim(config))
          case ("fit"); call initialize_temp_salt_fit(tv%T, tv%S, G, PF, eos, tv%P_Ref)
          case ("file"); call initialize_temp_salt_from_file(tv%T, tv%S, G, PF)
          case ("benchmark"); call benchmark_init_temperature_salinity(tv%T, tv%S, &
                                   G, PF, eos, tv%P_Ref)
          case ("TS_profile") ; call initialize_temp_salt_from_profile(tv%T, tv%S, G, PF)
          case ("linear"); call initialize_temp_salt_linear(tv%T, tv%S, G, PF)
          case ("DOME2D"); call DOME2d_initialize_temperature_salinity ( tv%T, &
                                tv%S, h, G, PF, eos )
          case ("adjustment2d"); call adjustment_initialize_temperature_salinity ( tv%T, &
                                      tv%S, h, G, PF, eos )
          case ("sloshing"); call sloshing_initialize_temperature_salinity(tv%T, &
                                  tv%S, h, G, PF, eos )
          case ("seamount"); call seamount_initialize_temperature_salinity(tv%T, &
                                  tv%S, h, G, PF, eos )
          case ("rossby_front"); call Rossby_front_initialize_temperature_salinity ( tv%T, &
                                tv%S, h, G, PF, eos )
          case ("USER"); call user_init_temperature_salinity(tv%T, tv%S, G, PF, eos)
          case default ; call MOM_error(FATAL,  "MOM_initialize: "//&
                 "Unrecognized Temp & salt configuration "//trim(config))
        end select
      endif
    endif  ! not from_Z_file.

    if (debug) then
      call hchksum(h, "MOM_initialize: h ", G, haloshift=1)
      if ( use_temperature ) call hchksum(tv%T, "MOM_initialize: T ", G, haloshift=1)
      if ( use_temperature ) call hchksum(tv%S, "MOM_initialize: S ", G, haloshift=1)
    endif

!   Initialize velocity components, u and v
    call get_param(PF, mod, "VELOCITY_CONFIG", config, &
         "A string that determines how the initial velocities \n"//&
         "are specified for a new run: \n"//&
         " \t file - read velocities from the file specified \n"//&
         " \t\t by (VELOCITY_FILE). \n"//&
         " \t zero - the fluid is initially at rest. \n"//&
         " \t uniform - the flow is uniform (determined by\n"//&
         " \t\t parameters TORUS_U and TORUS_V).\n"//&
         " \t rossby_front - a mixed layer front in thermal wind balance.\n"//&
         " \t USER - call a user modified routine.", default="zero")
    select case (trim(config))
       case ("file"); call initialize_velocity_from_file(u, v, G, PF)
       case ("zero"); call initialize_velocity_zero(u, v, G, PF)
       case ("uniform"); call initialize_velocity_uniform(u, v, G, PF)
       case ("circular"); call initialize_velocity_circular(u, v, G, PF)
       case ("phillips"); call Phillips_initialize_velocity(u, v, G, PF)
       case ("rossby_front"); call Rossby_front_initialize_velocity(u, v, h, G, PF)
       case ("USER"); call user_initialize_velocity(u, v, G, PF)
       case default ; call MOM_error(FATAL,  "MOM_initialize: "//&
            "Unrecognized velocity configuration "//trim(config))
    end select

    call pass_vector(u, v, G%Domain)
    if (debug) call uchksum(u, "MOM_initialize: u ", G, haloshift=1)
    if (debug) call vchksum(v, "MOM_initialize: v ", G, haloshift=1)


!   Optionally convert the thicknesses from m to kg m-2.  This is particularly
! useful in a non-Boussinesq model.
    call get_param(PF, mod, "CONVERT_THICKNESS_UNITS", convert, &
                 "If true,  convert the thickness initial conditions from \n"//&
                 "units of m to kg m-2 or vice versa, depending on whether \n"//&
                 "BOUSSINESQ is defined. This does not apply if a restart \n"//&
                 "file is read.", default=.false.)
    if (convert) call convert_thickness(h, G, PF, tv)

!  Remove the mass that would be displaced by an ice shelf or inverse barometer.
    call get_param(PF, mod, "DEPRESS_INITIAL_SURFACE", depress_sfc, &
                 "If true,  depress the initial surface to avoid huge \n"//&
                 "tsunamis when a large surface pressure is applied.", &
                 default=.false.)
    if (depress_sfc) call depress_surface(h, G, PF, tv)

  else ! Previous block for new_sim=.T., this block restores state
!    This line calls a subroutine that reads the initial conditions  !
!  from a previously generated file.                                 !
    call restore_state(dirs%input_filename, dirs%restart_input_dir, Time, &
                       G, restart_CS)
    if (present(Time_in)) Time = Time_in
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
                 " \t DOME - use a slope and channel configuration for the \n"//&
                 " \t\t DOME sill-overflow test case. \n"//&
                 " \t USER - call a user modified routine.", default="file")

    select case (trim(config))
      case ("DOME"); call DOME_initialize_sponges(G, tv, PF, CS%sponge_CSp)
      case ("USER"); call user_initialize_sponges(G, use_temperature, tv, &
                                                  PF, CS%sponge_CSp, h)
      case ("phillips"); call Phillips_initialize_sponges(G, use_temperature, tv, &
                                                  PF, CS%sponge_CSp, h)
      case ("file"); call initialize_sponges_file(G, use_temperature, tv, &
                                                  PF, CS%sponge_CSp)
      case default ; call MOM_error(FATAL,  "MOM_initialize: "//&
             "Unrecognized sponge configuration "//trim(config))
    end select

  endif

! This subroutine call sets optional open boundary conditions.
  call get_param(PF, mod, "APPLY_OBC_U", apply_OBC_u, &
                 "If true, open boundary conditions may be set at some \n"//&
                 "u-points, with the configuration controlled by OBC_CONFIG", &
                 default=.false.)
  call get_param(PF, mod, "APPLY_OBC_V", apply_OBC_v, &
                 "If true, open boundary conditions may be set at some \n"//&
                 "v-points, with the configuration controlled by OBC_CONFIG", &
                 default=.false.)
  if (apply_OBC_u .or. apply_OBC_v) then 
    call get_param(PF, mod, "OBC_CONFIG", config, &
                 "A string that sets how the open boundary conditions are \n"//&
                 " configured: \n"//&
                 " \t DOME - use a slope and channel configuration for the \n"//&
                 " \t\t DOME sill-overflow test case. \n"//&
                 " \t USER - call a user modified routine.", default="file", &
                 fail_if_missing=.true.)
    if (trim(config) == "DOME") then
      call DOME_set_Open_Bdry_Conds(CS%OBC, tv, G, PF, CS%tracer_Reg)
    elseif (trim(config) == "USER") then
      call user_set_Open_Bdry_Conds(CS%OBC, tv, G, PF, CS%tracer_Reg)
    else
      call MOM_error(FATAL, "The open boundary conditions specified by "//&
              "OBC_CONFIG = "//trim(config)//" have not been fully implemented.")
      call set_Open_Bdry_Conds(CS%OBC, tv, G, PF, CS%tracer_Reg)
    endif
  endif

  call get_param(PF, mod, "APPLY_OBC_U_FLATHER_EAST", apply_OBC_u_flather_east,&
                 "Apply a Flather open boundary condition on the eastern \n"//&
                 "side of the global domain", default=.false.)
  call get_param(PF, mod, "APPLY_OBC_U_FLATHER_WEST", apply_OBC_u_flather_west,&
                 "Apply a Flather open boundary condition on the western \n"//&
                 "side of the global domain", default=.false.)
  call get_param(PF, mod, "APPLY_OBC_V_FLATHER_NORTH", apply_OBC_v_flather_north,&
                 "Apply a Flather open boundary condition on the northern \n"//&
                 "side of the global domain", default=.false.)
  call get_param(PF, mod, "APPLY_OBC_V_FLATHER_SOUTH", apply_OBC_v_flather_south,&
                 "Apply a Flather open boundary condition on the southern \n"//&
                 "side of the global domain", default=.false.)
  if (apply_OBC_u_flather_east .or. apply_OBC_u_flather_west .or. apply_OBC_v_flather_north .or. apply_OBC_v_flather_south) then
    call set_Flather_Bdry_Conds(CS%OBC, tv, h, G, PF, CS%tracer_Reg)
  endif

  call callTree_leave('MOM_initialize()')

end subroutine MOM_initialize
! -----------------------------------------------------------------------------


! -----------------------------------------------------------------------------
subroutine set_coord_from_gprime(Rlay, g_prime, G, param_file)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(param_file_type), intent(in)  :: param_file
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real :: g_int   ! Reduced gravities across the internal interfaces, in m s-2.
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  character(len=40)  :: mod = "set_coord_from_gprime" ! This subroutine's name.
  integer :: k, nz
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, mod, "GFS" , g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call get_param(param_file, mod, "GINT", g_int, &
                 "The reduced gravity across internal interfaces.", &
                 units="m s-2", fail_if_missing=.true.)

  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = g_int ; enddo
  Rlay(1) = G%Rho0
  do k=2,nz ; Rlay(k) = Rlay(k-1) + g_prime(k)*(G%Rho0/G%g_Earth) ; enddo

  call callTree_leave(trim(mod)//'()')

end subroutine set_coord_from_gprime
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_layer_density(Rlay, g_prime, G, param_file)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(param_file_type), intent(in)  :: param_file
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  real :: Rlay_Ref! The surface layer's target density, in kg m-3.
  real :: RLay_range ! The range of densities, in kg m-3.
  character(len=40)  :: mod = "set_coord_from_layer_density" ! This subroutine's name.
  integer :: k, nz
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call get_param(param_file, mod, "LIGHTEST_DENSITY", Rlay_Ref, &
                 "The reference potential density used for layer 1.", &
                 units="kg m-3", default=G%Rho0)
  call get_param(param_file, mod, "DENSITY_RANGE", Rlay_range, &
                 "The range of reference potential densities in the layers.", &
                 units="kg m-3", default=2.0)

  g_prime(1) = g_fs
  Rlay(1) = Rlay_Ref
  do k=2,nz
     Rlay(k) = Rlay(k-1) + RLay_range/(real(nz-1))
  enddo
!    These statements set the interface reduced gravities.           !
  do k=2,nz
     g_prime(k) = (G%g_Earth/G%Rho0) * (Rlay(k) - Rlay(k-1))
  enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_coord_from_layer_density
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_TS_ref(Rlay, g_prime, G, param_file, eqn_of_state, &
                                 P_Ref)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(param_file_type), intent(in)  :: param_file
  type(EOS_type),        pointer     :: eqn_of_state
  real,                  intent(in)  :: P_Ref
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      eqn_of_state - integer selecting the equation of state
!  (in)      P_Ref - The coordinate-density reference pressure in Pa.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real :: T_ref   ! Reference temperature
  real :: S_ref   ! Reference salinity
  real :: g_int   ! Reduced gravities across the internal interfaces, in m s-2.
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  character(len=40)  :: mod = "set_coord_from_TS_ref" ! This subroutine's name.
  integer :: k, nz
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, mod, "T_REF", T_Ref, &
                 "The initial temperature of the lightest layer.", units="degC", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "S_REF", S_Ref, &
                 "The initial salinities.", units="PSU", default=35.0)
  call get_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call get_param(param_file, mod, "GINT", g_int, &
                 "The reduced gravity across internal interfaces.", &
                 units="m s-2", fail_if_missing=.true.)
                                      !
!    These statements set the interface reduced gravities.           !
  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = g_int ; enddo

!    The uppermost layer's density is set here.  Subsequent layers'  !
!  densities are determined from this value and the g values.        !
!        T0 = 28.228 ; S0 = 34.5848 ; Pref = P_Ref
  call calculate_density(T_ref, S_ref, P_ref, Rlay(1), eqn_of_state)

!    These statements set the layer densities.                       !
  do k=2,nz ; Rlay(k) = Rlay(k-1) + g_prime(k)*(G%Rho0/G%g_Earth) ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_coord_from_TS_ref
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_TS_profile(Rlay, g_prime, G, param_file, &
                                     eqn_of_state, P_Ref)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(param_file_type), intent(in)  :: param_file
  type(EOS_type),        pointer     :: eqn_of_state
  real,                  intent(in)  :: P_Ref
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      eqn_of_state - integer that selects equation of state
!  (in)      P_Ref - The coordinate-density reference pressure in Pa.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real, dimension(SZK_(G)) :: T0, S0,  Pref
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  integer :: k, nz
  character(len=40)  :: mod = "set_coord_from_TS_profile" ! This subroutine's name.
  character(len=200) :: filename, coord_file, inputdir ! Strings for file/path
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call get_param(param_file, mod, "COORD_FILE", coord_file, &
                 "The file from which the coordinate temperatures and \n"//&
                 "salnities are read.", fail_if_missing=.true.)

  call get_param(param_file,  mod, "INPUTDIR", inputdir, default=".")
  filename = trim(slasher(inputdir))//trim(coord_file)
  call log_param(param_file, mod, "INPUTDIR/COORD_FILE", filename)

  call read_data(filename,"PTEMP",T0(:),domain=G%Domain%mpp_domain)
  call read_data(filename,"SALT",S0(:),domain=G%Domain%mpp_domain)

  if (.not.file_exists(filename)) call MOM_error(FATAL, &
      " set_coord_from_TS_profile: Unable to open " //trim(filename))
!    These statements set the interface reduced gravities.           !
  g_prime(1) = g_fs
  do k=1,nz ; Pref(k) = P_ref ; enddo
  call calculate_density(T0, S0, Pref, Rlay, 1,nz,eqn_of_state)
  do k=2,nz; g_prime(k) = (G%g_Earth/G%Rho0) * (Rlay(k) - Rlay(k-1)); enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_coord_from_TS_profile
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_TS_range(Rlay, g_prime, G, param_file, &
                                     eqn_of_state, P_Ref)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(param_file_type), intent(in)  :: param_file
  type(EOS_type),        pointer     :: eqn_of_state
  real,                  intent(in)  :: P_Ref
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      eqn_of_state - integer that selects equation of state
!  (in)      P_Ref - The coordinate-density reference pressure in Pa.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real, dimension(SZK_(G)) :: T0, S0,  Pref
  real :: S_Ref, S_Light, S_Dense ! Salnity range parameters in PSU.
  real :: T_Ref, T_Light, T_Dense ! Temperature range parameters in dec C.
  real :: res_rat ! The ratio of density space resolution in the denser part
                  ! of the range to that in the lighter part of the range.
                  ! Setting this greater than 1 increases the resolution for
                  ! the denser water.
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  real :: a1, frac_dense, k_frac
  integer :: k, nz, k_light
  character(len=40)  :: mod = "set_coord_from_TS_range" ! This subroutine's name.
  character(len=200) :: filename, coord_file, inputdir ! Strings for file/path
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, mod, "T_REF", T_Ref, &
                 "The default initial temperatures.", units="degC", default=10.0)
  call get_param(param_file, mod, "TS_RANGE_T_LIGHT", T_Light, &
                 "The initial temperature of the lightest layer when \n"//&
                 "COORD_CONFIG is set to ts_range.", units="degC", default=T_Ref)
  call get_param(param_file, mod, "TS_RANGE_T_DENSE", T_Dense, &
                 "The initial temperature of the densest layer when \n"//&
                 "COORD_CONFIG is set to ts_range.", units="degC", default=T_Ref)

  call get_param(param_file, mod, "S_REF", S_Ref, &
                 "The default initial salinities.", units="PSU", default=35.0)
  call get_param(param_file, mod, "TS_RANGE_S_LIGHT", S_Light, &
                 "The initial lightest salinities when COORD_CONFIG \n"//&
                 "is set to ts_range.", default = S_Ref, units="PSU")
  call get_param(param_file, mod, "TS_RANGE_S_DENSE", S_Dense, &
                 "The initial densest salinities when COORD_CONFIG \n"//&
                 "is set to ts_range.", default = S_Ref, units="PSU")
 
  call get_param(param_file, mod, "TS_RANGE_RESOLN_RATIO", res_rat, &
                 "The ratio of density space resolution in the densest \n"//&
                 "part of the range to that in the lightest part of the \n"//&
                 "range when COORD_CONFIG is set to ts_range. Values \n"//&
                 "greater than 1 increase the resolution of the denser water.",&
                 default=1.0, units="nondim")

  call get_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)

  k_light = G%nk_rho_varies + 1

  ! Set T0(k) to range from T_LIGHT to T_DENSE, and simliarly for S0(k).
  T0(k_light) = T_light ; S0(k_light) = S_light
  a1 = 2.0 * res_rat / (1.0 + res_rat)
  do k=k_light+1,nz
    k_frac = real(k-k_light)/real(nz-k_light)
    frac_dense = a1 * k_frac + (1.0 - a1) * k_frac**2
    T0(k) = frac_dense * (T_Dense - T_Light) + T_Light
    S0(k) = frac_dense * (S_Dense - S_Light) + S_Light
  enddo

  g_prime(1) = g_fs
  do k=1,nz ; Pref(k) = P_ref ; enddo
  call calculate_density(T0, S0, Pref, Rlay, k_light,nz-k_light+1,eqn_of_state)
  ! Extrapolate target densities for the variable density mixed and buffer layers.
  do k=k_light-1,1,-1
    Rlay(k) = 2.0*Rlay(k+1) -  Rlay(k+2)
  enddo
  do k=2,nz; g_prime(k) = (G%g_Earth/G%Rho0) * (Rlay(k) - Rlay(k-1)); enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_coord_from_TS_range
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_file(Rlay, g_prime, G, param_file)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(param_file_type), intent(in)  :: param_file
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  integer :: k, nz
  character(len=40)  :: mod = "set_coord_from_file" ! This subroutine's name.
  character(len=40)  :: coord_var
  character(len=200) :: filename,coord_file,inputdir ! Strings for file/path
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mod, "COORD_FILE", coord_file, &
                 "The file from which the coordinate densities are read.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "COORD_VAR", coord_var, &
                 "The variable in COORD_FILE that is to be used for the \n"//&
                 "coordinate densities.", default="Layer")
  filename = trim(inputdir)//trim(coord_file)
  call log_param(param_file, mod, "INPUTDIR/COORD_FILE", filename)
  if (.not.file_exists(filename)) call MOM_error(FATAL, &
      " set_coord_from_file: Unable to open "//trim(filename))

  call read_axis_data(filename, coord_var, Rlay)
  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = (G%g_Earth/G%Rho0) * (Rlay(k) - Rlay(k-1)) ; enddo
  do k=1,nz ; if (g_prime(k) <= 0.0) then
    call MOM_error(FATAL, "MOM_initialization set_coord_from_file: "//&
       "Zero or negative g_primes read from variable "//"Layer"//" in file "//&
       trim(filename))
  endif ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_coord_from_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_linear(Rlay, g_prime, G, param_file)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(param_file_type), intent(in)  :: param_file
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets the layer densities (Rlay) and the interface  
! reduced gravities (g) according to a linear profile starting at a
! reference surface layer density and spanning a range of densities 
! defined by the parameter RLAY_RANGE (defaulted to 2 if not defined) 
  character(len=40)  :: mod = "set_coord_linear" ! This subroutine
  real :: Rlay_ref, Rlay_range, g_fs
  integer :: k, nz
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, mod, "LIGHTEST_DENSITY", Rlay_Ref, &
                 "The reference potential density used for layer 1.", &
                 units="kg m-3", default=G%Rho0)
  call get_param(param_file, mod, "DENSITY_RANGE", Rlay_range, &
                 "The range of reference potential densities in the layers.", &
                 units="kg m-3", default=2.0)
  call get_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
! Rlay(1) = Rlay_Ref
! do k=2,nz
!    Rlay(k) = Rlay(k-1) + Rlay_range/(real(nz-1))
! enddo
  ! This following sets the target layer densities such that a the
  ! surface interface has density Rlay_ref and the bottom
  ! is Rlay_range larger
  do k=1,nz
     Rlay(k) = Rlay_Ref + RLay_range*((real(k)-0.5)/real(nz))
  enddo
!    These statements set the interface reduced gravities.           !
  g_prime(1) = g_fs
  do k=2,nz 
     g_prime(k) = (G%g_Earth/G%Rho0) * (Rlay(k) - Rlay(k-1))
  enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_coord_linear
! -----------------------------------------------------------------------------

subroutine MOM_initialize_rotation(f, G, PF)
  type(ocean_grid_type),                        intent(in)  :: G
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), intent(out) :: f
  type(param_file_type),                        intent(in)  :: PF
! Arguments: f  - the Coriolis parameter in s-1. Intent out.
!  (in)      G  - The ocean's grid structure.
!  (in)      PF - A structure indicating the open file to parse for
!                         model parameter values.

!   This subroutine makes the appropriate call to set up the Coriolis parameter.
! This is a separate subroutine so that it can be made public and shared with
! the ice-sheet code or other components.
! Set up the Coriolis parameter, f, either analytically or from file.
  character(len=40)  :: mod = "MOM_initialize_rotation" ! This subroutine's name.
  character(len=200) :: config

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")
  call get_param(PF, mod, "ROTATION", config, &
                 "This specifies how the Coriolis parameter is specified: \n"//&
                 " \t 2omegasinlat - Use twice the planetary rotation rate \n"//&
                 " \t\t times the sine of latitude.\n"//&
                 " \t betaplane - Use a beta-plane or f-plane. \n"//&
                 " \t USER - call a user modified routine.", &
                 default="2omegasinlat")
  select case (trim(config))
    case ("2omegasinlat"); call set_rotation_planetary(f, G, PF)
    case ("beta"); call set_rotation_beta_plane(f, G, PF)
    case ("betaplane"); call set_rotation_beta_plane(f, G, PF)
   !case ("nonrotating") ! Note from AJA: Missing case?
    case default ; call MOM_error(FATAL,"MOM_initialize: "// &
      "Unrecognized rotation setup "//trim(config))
  end select
  call callTree_leave(trim(mod)//'()')
end subroutine MOM_initialize_rotation

subroutine MOM_initialize_topography(D, max_depth, G, PF)
  real, dimension(NIMEM_,NJMEM_), intent(out) :: D
  real,                           intent(out) :: max_depth
  type(ocean_grid_type),          intent(in)  :: G
  type(param_file_type),          intent(in)  :: PF
! Arguments: D  - the bottom depth in m. Intent out.
!  (in)      G  - The ocean's grid structure.
!  (in)      PF - A structure indicating the open file to parse for
!                         model parameter values.

!  This subroutine makes the appropriate call to set up the bottom depth.
!  This is a separate subroutine so that it can be made public and shared with
!  the ice-sheet code or other components.
! Set up the bottom depth, G%bathyT either analytically or from file
  character(len=40)  :: mod = "MOM_initialize_topography" ! This subroutine's name.
  character(len=200) :: config

  call get_param(PF, mod, "TOPO_CONFIG", config, &
                 "This specifies how bathymetry is specified: \n"//&
                 " \t file - read bathymetric information from the file \n"//&
                 " \t\t specified by (TOPO_FILE).\n"//&
                 " \t flat - flat bottom set to MAXIMUM_DEPTH. \n"//&
                 " \t bowl - an analytically specified bowl-shaped basin \n"//&
                 " \t\t ranging between MAXIMUM_DEPTH and MINIMUM_DEPTH. \n"//&
                 " \t spoon - a similar shape to 'bowl', but with an vertical \n"//&
                 " \t\t wall at the southern face. \n"//&
                 " \t halfpipe - a zonally uniform channel with a half-sine \n"//&
                 " \t\t profile in the meridional direction. \n"//&
                 " \t benchmark - use the benchmark test case topography. \n"//&
                 " \t DOME - use a slope and channel configuration for the \n"//&
                 " \t\t DOME sill-overflow test case. \n"//&
                 " \t DOME2D - use a shelf and slope configuration for the \n"//&
                 " \t\t DOME2D gravity current/overflow test case. \n"//&
                 " \t seamount - Gaussian bump for spontaneous motion test case.\n"//&
                 " \t USER - call a user modified routine.", &
                 fail_if_missing=.true.)
  max_depth = -1.e9; call read_param(PF, "MAXIMUM_DEPTH", max_depth)
  select case ( trim(config) )
    case ("file");      call initialize_topography_from_file(D, G, PF)
    case ("flat");      call initialize_topography_named(D, G, PF, config, max_depth)
    case ("spoon");     call initialize_topography_named(D, G, PF, config, max_depth)
    case ("bowl");      call initialize_topography_named(D, G, PF, config, max_depth)
    case ("halfpipe");  call initialize_topography_named(D, G, PF, config, max_depth)
    case ("DOME");      call DOME_initialize_topography(D, G, PF, max_depth)
    case ("benchmark"); call benchmark_initialize_topography(D, G, PF, max_depth)
    case ("DOME2D");    call DOME2d_initialize_topography(D, G, PF, max_depth)
    case ("sloshing");  call sloshing_initialize_topography(D, G, PF, max_depth)
    case ("seamount");  call seamount_initialize_topography(D, G, PF, max_depth)
    case ("USER");      call user_initialize_topography(D, G, PF)
    case default ;      call MOM_error(FATAL,"MOM_initialize_topography: "// &
      "Unrecognized topography setup '"//trim(config)//"'")
  end select
  if (max_depth>0.) then
    call log_param(PF, mod, "MAXIMUM_DEPTH", max_depth, &
                   "The maximum depth of the ocean.", units="m")
  else
    max_depth = diagnoseMaximumDepth(D,G)
    call log_param(PF, mod, "!MAXIMUM_DEPTH", max_depth, &
                   "The (diagnosed) maximum depth of the ocean.", units="m")
  endif
  if (trim(config) .ne. "DOME") then
    call limit_topography(D, G, PF, max_depth)
  endif
  
end subroutine MOM_initialize_topography

! -----------------------------------------------------------------------------
function diagnoseMaximumDepth(D,G)
  real, dimension(NIMEM_,NJMEM_), intent(in) :: D
  type(ocean_grid_type),          intent(in) :: G
  real :: diagnoseMaximumDepth
  ! Local variables
  integer :: i,j
  diagnoseMaximumDepth=D(G%isc,G%jsc)
  do j=G%jsc, G%jec
    do i=G%isc, G%iec
      diagnoseMaximumDepth=max(diagnoseMaximumDepth,D(i,j))
    enddo
  enddo
  call max_across_PEs(diagnoseMaximumDepth)
end function diagnoseMaximumDepth
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_topography_from_file(D, G, param_file )
  real, dimension(NIMEM_,NJMEM_), intent(out) :: D
  type(ocean_grid_type),          intent(in)  :: G
  type(param_file_type),          intent(in)  :: param_file
! Arguments: D          - the bottom depth in m. Intent out.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!  This subroutine reads depths from a file and puts it into D(:,:) in m.
  character(len=200) :: filename, topo_file, inputdir ! Strings for file/path
  character(len=200) :: topo_varname                  ! Variable name in file
  character(len=40)  :: mod = "initialize_topography_from_file" ! This subroutine's name.

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mod, "TOPO_FILE", topo_file, &
                 "The file from which the bathymetry is read.", &
                 default="topog.nc")
  call get_param(param_file, mod, "TOPO_VARNAME", topo_varname, &
                 "The name of the bathymetry variable in TOPO_FILE.", &
                 default="depth")

  filename = trim(inputdir)//trim(topo_file)
  call log_param(param_file, mod, "INPUTDIR/TOPO_FILE", filename)

  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topography_from_file: Unable to open "//trim(filename))

  call read_data(filename,trim(topo_varname),D,domain=G%Domain%mpp_domain)

  call callTree_leave(trim(mod)//'()')
end subroutine initialize_topography_from_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_topography_named(D, G, param_file, topog_config, max_depth)
  real, dimension(NIMEM_,NJMEM_), intent(out) :: D
  type(ocean_grid_type),          intent(in)  :: G
  type(param_file_type),          intent(in)  :: param_file
  character(len=*),               intent(in)  :: topog_config
  real,                           intent(in)  :: max_depth
! Arguments: D          - the bottom depth in m. Intent out.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      topog_config - The name of an idealized topographic configuration.
!  (in)      max_depth  - The maximum depth in m.

! This subroutine places the bottom depth in m into D(:,:), shaped in a spoon
  real :: min_depth            ! The minimum depth in m.
  real :: PI                   ! 3.1415926... calculated as 4*atan(1)
  real :: D0                   ! A constant to make the maximum     !
                               ! basin depth MAXIMUM_DEPTH.         !
  real :: expdecay             ! A decay scale of associated with   !
                               ! the sloping boundaries, in m.      !
  real :: Dedge                ! The depth in m at the basin edge.  !
! real :: south_lat, west_lon, len_lon, len_lat, Rad_earth
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  character(len=40)  :: mod = "initialize_topography_named" ! This subroutine's name.
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")
  call MOM_mesg("  MOM_initialization.F90, initialize_topography_named: "//&
                 "TOPO_CONFIG = "//trim(topog_config), 5)

  call get_param(param_file, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)
  if (max_depth<=0.) call MOM_error(FATAL,"initialize_topography_named: "// &
      "MAXIMUM_DEPTH has a non-sensical value! Was it set?")

  if (trim(topog_config) /= "flat") then
    call get_param(param_file, mod, "EDGE_DEPTH", Dedge, &
                   "The depth at the edge of one of the named topographies.", &
                   units="m", default=100.0)
!   call get_param(param_file, mod, "SOUTHLAT", south_lat, &
!                  "The southern latitude of the domain.", units="degrees", &
!                  fail_if_missing=.true.)
!   call get_param(param_file, mod, "LENLAT", len_lat, &
!                  "The latitudinal length of the domain.", units="degrees", &
!                  fail_if_missing=.true.)
!   call get_param(param_file, mod, "WESTLON", west_lon, &
!                  "The western longitude of the domain.", units="degrees", &
!                  default=0.0)
!   call get_param(param_file, mod, "LENLON", len_lon, &
!                  "The longitudinal length of the domain.", units="degrees", &
!                  fail_if_missing=.true.)
!   call get_param(param_file, mod, "RAD_EARTH", Rad_Earth, &
!                  "The radius of the Earth.", units="m", default=6.378e6)
    call get_param(param_file, mod, "TOPOG_SLOPE_SCALE", expdecay, &
                   "The exponential decay scale used in defining some of \n"//&
                   "the named topographies.", units="m", default=400000.0)
  endif


  PI = 4.0*atan(1.0)

  if (trim(topog_config) == "flat") then
    do i=is,ie ; do j=js,je ; D(i,j) = max_depth ; enddo ; enddo
  elseif (trim(topog_config) == "spoon") then
    D0 = (max_depth - Dedge) / &
             ((1.0 - exp(-0.5*G%len_lat*G%Rad_earth*PI/(180.0 *expdecay))) * &
              (1.0 - exp(-0.5*G%len_lat*G%Rad_earth*PI/(180.0 *expdecay))))
    do i=is,ie ; do j=js,je
  !  This sets a bowl shaped (sort of) bottom topography, with a       !
  !  maximum depth of max_depth.                                   !
      D(i,j) =  Dedge + D0 * &
             (sin(PI * (G%geoLonT(i,j) - (G%west_lon)) / G%len_lon) * &
           (1.0 - exp((G%geoLatT(i,j) - (G%south_lat+G%len_lat))*G%Rad_earth*PI / &
                      (180.0*expdecay)) ))
    enddo ; enddo
  elseif (trim(topog_config) == "bowl") then
    D0 = (max_depth - Dedge) / &
             ((1.0 - exp(-0.5*G%len_lat*G%Rad_earth*PI/(180.0 *expdecay))) * &
              (1.0 - exp(-0.5*G%len_lat*G%Rad_earth*PI/(180.0 *expdecay))))

  !  This sets a bowl shaped (sort of) bottom topography, with a
  !  maximum depth of max_depth.
    do i=is,ie ; do j=js,je
      D(i,j) =  Dedge + D0 * &
             (sin(PI * (G%geoLonT(i,j) - G%west_lon) / G%len_lon) * &
             ((1.0 - exp(-(G%geoLatT(i,j) - G%south_lat)*G%Rad_Earth*PI/ &
                          (180.0*expdecay))) * &
             (1.0 - exp((G%geoLatT(i,j) - (G%south_lat+G%len_lat))* &
                         G%Rad_Earth*PI/(180.0*expdecay)))))
    enddo ; enddo
  elseif (trim(topog_config) == "halfpipe") then
    D0 = max_depth - Dedge
    do i=is,ie ; do j=js,je
      D(i,j) =  Dedge + D0 * ABS(sin(PI*(G%geoLatT(i,j) - G%south_lat)/G%len_lat))
    enddo ; enddo
  else
    call MOM_error(FATAL,"initialize_topography_named: "// &
      "Unrecognized topography name "//trim(topog_config))
  endif

  ! This is here just for safety.  Hopefully it doesn't do anything.
  do i=is,ie ; do j=js,je
    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
  enddo ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine initialize_topography_named
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine limit_topography(D, G, param_file, max_depth)
  real, dimension(NIMEM_,NJMEM_), intent(inout) :: D
  type(ocean_grid_type),          intent(in)    :: G
  type(param_file_type),          intent(in)    :: param_file
  real,                           intent(in)    :: max_depth
! Arguments: D          - the bottom depth in m. Intent in/out.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      max_depth  - The maximum depth in m.

! This subroutine ensures that    min_depth < D(x,y) < max_depth
  integer :: i, j
  character(len=40)  :: mod = "limit_topography" ! This subroutine's name.
  real :: min_depth, mask_depth

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, "MOM_grid_init initialize_masks", "MINIMUM_DEPTH", min_depth, &
                 "If MASKING_DEPTH is unspecified, then anything shallower than\n"//&
                 "MINIMUM_DEPTH is assumed to be land and all fluxes are masked out.\n"//&
                 "If MASKING_DEPTH is specified, then all depths shallower than\n"//&
                 "MINIMUM_DEPTH but depper than MASKING_DEPTH are rounded to MINIMUM_DEPTH.", &
                 units="m", default=0.0)
  call get_param(param_file, mod, "MASKING_DEPTH", mask_depth, &
                 "The depth below which to mask the ocean as land.", units="m", &
                 default=-9999.0, do_not_log=.true.)

! Make sure that min_depth < D(x,y) < max_depth
  if (mask_depth<-9990.) then
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      D(i,j) = min( max( D(i,j), 0.5*min_depth ), max_depth )
    enddo ; enddo
  else
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      if (D(i,j)>0.) then
        D(i,j) = min( max( D(i,j), min_depth ), max_depth )
      else
        D(i,j) = 0.
      endif
    enddo ; enddo
  endif

  call callTree_leave(trim(mod)//'()')
end subroutine limit_topography
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_rotation_planetary(f, G, param_file)
  type(ocean_grid_type),                        intent(in)  :: G
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), intent(out) :: f
  type(param_file_type),                        intent(in)  :: param_file
! Arguments: f          - Coriolis parameter (vertical component) in s^-1
!     (in)   G          - grid type
!     (in)   param_file - parameter file type

! This subroutine sets up the Coriolis parameter for a sphere
  character(len=30) :: mod = "set_rotation_planetary" ! This subroutine's name.
  integer :: I, J
  real    :: PI, omega

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, "set_rotation_planetary", "OMEGA", omega, &
                 "The rotation rate of the earth.", units="s-1", &
                 default=7.2921e-5)
  PI = 4.0*atan(1.0)

  do I=G%IsdB,G%IedB ; do J=G%JsdB,G%JedB
    f(I,J) = ( 2.0 * omega ) * sin( ( PI * G%geoLatBu(I,J) ) / 180.)
  enddo ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_rotation_planetary
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_rotation_beta_plane(f, G, param_file)
  type(ocean_grid_type),                        intent(in)  :: G
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), intent(out) :: f
  type(param_file_type),                        intent(in)  :: param_file
! Arguments: f          - Coriolis parameter (vertical component) in s^-1
!     (in)   G          - grid type
!     (in)   param_file - parameter file type

! This subroutine sets up the Coriolis parameter for a beta-plane
  integer :: I, J
  real    :: f_0, beta, y_scl, Rad_Earth, PI
  character(len=40)  :: mod = "set_rotation_beta_plane" ! This subroutine's name.
  character(len=200) :: axis_units

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, mod, "F_0", f_0, &
                 "The reference value of the Coriolis parameter with the \n"//&
                 "betaplane option.", units="s-1", default=0.0)
  call get_param(param_file, mod, "BETA", beta, &
                 "The northward gradient of the Coriolis parameter with \n"//&
                 "the betaplane option.", units="m-1 s-1", default=0.0)
  call get_param(param_file, mod, "AXIS_UNITS", axis_units, default="degrees")

  PI = 4.0*atan(1.0)
  select case (axis_units(1:1))
    case ("d")
      call get_param(param_file, mod, "RAD_EARTH", Rad_Earth, &
                   "The radius of the Earth.", units="m", default=6.378e6)
      y_scl = Rad_Earth/PI
    case ("k"); y_scl = 1.E3
    case ("m"); y_scl = 1.
    case ("c"); y_scl = 1.E-2
    case default ; call MOM_error(FATAL, &
      " set_rotation_beta_plane: unknown AXIS_UNITS = "//trim(axis_units))
  end select

  do I=G%IsdB,G%IedB ; do J=G%JsdB,G%JedB
    f(I,J) = f_0 + beta * ( G%geoLatBu(I,J) * y_scl )
  enddo ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_rotation_beta_plane
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_thickness_from_file(h, G, param_file, file_has_thickness)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: h
  type(ocean_grid_type),                  intent(in)  :: G
  type(param_file_type),                  intent(in)  :: param_file
  logical,                                intent(in)  :: file_has_thickness
! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
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

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

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
      call adjustEtaToFitBathymetry(G, eta, h)
    else
      do k=nz,1,-1 ; do j=js,je ; do i=is,ie
        if (eta(i,j,K) < (eta(i,j,K+1) + G%Angstrom_z)) then
          eta(i,j,K) = eta(i,j,K+1) + G%Angstrom_z
          h(i,j,k) = G%Angstrom_z
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
!! layers are contracted to G%Angstrom_z.
!! If the bottom most interface is above the topography then the entire column
!! is dilated (expanded) to fill the void.
!!   @remark{There is a (hard-wired) "tolerance" parameter such that the
!! criteria for adjustment must equal or exceed 10cm.}
!!   @param[in]     G   Grid type
!!   @param[in,out] eta Interface heights
!!   @param[out]    h   Layer thicknesses
subroutine adjustEtaToFitBathymetry(G, eta, h)
  type(ocean_grid_type),                          intent(in)    :: G
  real, dimension(NIMEM_,NJMEM_, NK_INTERFACE_ ), intent(inout) :: eta
  real, dimension(NIMEM_,NJMEM_, NKMEM_),         intent(inout) :: h
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
    if (eta(i,j,K) < (eta(i,j,K+1) + G%Angstrom_z)) then
      eta(i,j,K) = eta(i,j,K+1) + G%Angstrom_z
      h(i,j,k) = G%Angstrom_z
    else
      h(i,j,k) = eta(i,j,K) - eta(i,j,K+1)
    endif
  enddo ; enddo ; enddo

  dilations = 0
  do j=js,je ; do i=is,ie
    !   The whole column is dilated to accomodate deeper topography than
    ! the bathymetry would indicate.
    if (-eta(i,j,nz+1) < G%bathyT(i,j) - hTolerance) then
      dilations = dilations + 1
      dilate = (eta(i,j,1)+G%bathyT(i,j)) / (eta(i,j,1)-eta(i,j,nz+1))
      do k=1,nz ; h(i,j,k) = h(i,j,k) * dilate ; enddo
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
subroutine initialize_thickness_uniform(h, G, param_file)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: h
  type(ocean_grid_type),                  intent(in)  :: G
  type(param_file_type),                  intent(in)  :: param_file

! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
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

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

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
      if (eta1D(K) < (eta1D(K+1) + G%Angstrom_z)) then
        eta1D(K) = eta1D(K+1) + G%Angstrom_z
        h(i,j,k) = G%Angstrom_z
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
  call MOM_error(FATAL,"  MOM_initialization.F90, initialize_thickness_search: NOT IMPLEMENTED")
end subroutine initialize_thickness_search
! -----------------------------------------------------------------------------

subroutine convert_thickness(h, G, param_file, tv)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout) :: h
  type(ocean_grid_type),                  intent(in)    :: G
  type(param_file_type),                  intent(in)    :: param_file
  type(thermo_var_ptrs),                  intent(in)    :: tv
! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
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
  Boussinesq = G%Boussinesq
  I_gEarth = 1.0 / G%g_Earth

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
            p_bot(i,j) = p_top(i,j) + G%g_Earth * h(i,j,k) * rho(i)
          enddo
        enddo

        do itt=1,max_itt
          call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p_top, p_bot, &
                                   0.0, G, tv%eqn_of_state, dz_geo)
          if (itt < max_itt) then ; do j=js,je
            call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_bot(:,j), rho, &
                                   is, ie-is+1, tv%eqn_of_state)
            ! Use Newton's method to correct the bottom value.
            !   The hydrostatic equation is linear to such a
            ! high degree that no bounds-checking is needed.
            do i=is,ie
              p_bot(i,j) = p_bot(i,j) + rho(i) * (G%g_Earth*h(i,j,k) - dz_geo(i,j))
            enddo
          enddo ; endif
        enddo

        do j=js,je ; do i=is,ie
          h(i,j,k) = (p_bot(i,j) - p_top(i,j)) * G%kg_m2_to_H * I_gEarth
        enddo ; enddo
      enddo
    else
      do k=1,nz ; do j=js,je ; do i=is,ie
        h(i,j,k) = h(i,j,k) * G%Rlay(k) * G%kg_m2_to_H
      enddo ; enddo ; enddo
    endif
  endif

end subroutine convert_thickness

subroutine depress_surface(h, G, param_file, tv)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout) :: h
  type(ocean_grid_type),                  intent(in)    :: G
  type(param_file_type),                  intent(in)    :: param_file
  type(thermo_var_ptrs),                  intent(in)    :: tv
! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
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
  call find_eta(h, tv, G%g_Earth, G, eta)
  
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
          h(i,j,k) = G%Angstrom
        else
          h(i,j,k) = max(G%Angstrom, h(i,j,k) * &
              (eta_sfc(i,j) - eta(i,j,K+1)) / (eta(i,j,K) - eta(i,j,K+1)) )
        endif
      enddo
    endif
  endif ; enddo ; enddo

end subroutine depress_surface

! -----------------------------------------------------------------------------
subroutine initialize_velocity_from_file(u, v, G, param_file)
  real, dimension(NIMEMB_,NJMEM_, NKMEM_), intent(out) :: u
  real, dimension(NIMEM_,NJMEMB_, NKMEM_), intent(out) :: v
  type(ocean_grid_type),                   intent(in)  :: G
  type(param_file_type),                   intent(in)  :: param_file
! Arguments: u - The zonal velocity that is being initialized.
!  (out)     v - The meridional velocity that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file -  parameter file type

!   This subroutine reads the initial velocity components from file
  character(len=40)  :: mod = "initialize_velocity_from_file" ! This subroutine's name.
  character(len=200) :: filename,velocity_file,inputdir ! Strings for file/path

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

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
  real, dimension(NIMEMB_,NJMEM_, NKMEM_), intent(out) :: u
  real, dimension(NIMEM_,NJMEMB_, NKMEM_), intent(out) :: v
  type(ocean_grid_type),                   intent(in)  :: G
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

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

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
  real, dimension(NIMEMB_,NJMEM_, NKMEM_), intent(out) :: u
  real, dimension(NIMEM_,NJMEMB_, NKMEM_), intent(out) :: v
  type(ocean_grid_type),                   intent(in)  :: G
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
  real, dimension(NIMEMB_,NJMEM_, NKMEM_), intent(out) :: u
  real, dimension(NIMEM_,NJMEMB_, NKMEM_), intent(out) :: v
  type(ocean_grid_type),                   intent(in)  :: G
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
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: T, S
  type(ocean_grid_type),                  intent(in)  :: G
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

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

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
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: T, S
  type(ocean_grid_type),                  intent(in)  :: G
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

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

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
subroutine initialize_temp_salt_fit(T, S, G, param_file, eqn_of_state, P_Ref)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: T, S
  type(ocean_grid_type),                  intent(in)  :: G
  type(param_file_type),                  intent(in)  :: param_file
  type(EOS_type),                         pointer     :: eqn_of_state
  real,                                   intent(in)  :: P_Ref
!  This function puts the initial layer temperatures and salinities  !
! into T(:,:,:) and S(:,:,:).                                        !

! Arguments: T - The potential temperature that is being initialized.
!  (out)     S - The salinity that is being initialized.
!  (in)      G - The ocean's grid structure.
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
  character(len=40)  :: mod = "initialize_temp_salt_fit" ! This subroutine's name.
  integer :: i, j, k, itt, nz
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, mod, "T_REF", T_Ref, &
                 "A reference temperature used in initialization.", &
                 units="degC", fail_if_missing=.true.)
  call get_param(param_file, mod, "S_REF", S_Ref, &
                 "A reference salinity used in initialization.", units="PSU", &
                 default=35.0)
  do k=1,nz
    pres(k) = P_Ref ; S0(k) = S_Ref
  enddo
  T0(1) = T_Ref

  call calculate_density(T0(1),S0(1),pres(1),rho_guess(1),eqn_of_state)
  call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,1,eqn_of_state)

! A first guess of the layers' temperatures.                         !
  do k=nz,1,-1
    T0(k) = T0(1) + (G%Rlay(k) - rho_guess(1)) / drho_dT(1)
  enddo

! Refine the guesses for each layer.                                 !
  do itt=1,6
    call calculate_density(T0,S0,pres,rho_guess,1,nz,eqn_of_state)
    call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,nz,eqn_of_state)
    do k=1,nz
      T0(k) = T0(k) + (G%Rlay(k) - rho_guess(k)) / drho_dT(k)
    enddo
  enddo

  do k=1,nz ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
    T(i,j,k) = T0(k) ; S(i,j,k) = S0(k)
  enddo ; enddo ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine initialize_temp_salt_fit
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine initialize_temp_salt_linear(T, S, G, param_file)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: T, S
  type(ocean_grid_type),                  intent(in)  :: G
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
  
  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")
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
subroutine initialize_sponges_file(G, use_temperature, tv, param_file, CSp)
  type(ocean_grid_type), intent(in) :: G
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
    if (eta(i,j,K) < (eta(i,j,K+1) + G%Angstrom_z)) &
      eta(i,j,K) = eta(i,j,K+1) + G%Angstrom_z
  enddo ; enddo ; enddo
! Set the inverse damping rates so that the model will know where to !
! apply the sponges, along with the interface heights.               !
  call initialize_sponge(Idamp, eta, G, param_file, CSp)

!   Now register all of the tracer fields which are damped in the    !
! sponge. By default, momentum is advected vertically within the     !
! sponge, but momentum is typically not damped within the sponge.    !

  if ( G%nkml>0 ) then
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

    call set_up_sponge_ML_density(tmp_2d, CSp)
  endif

!  The remaining calls to set_up_sponge_field can be in any order.   !
  if ( use_temperature ) then
    call read_data(filename, potemp_var, tmp(:,:,:), domain=G%Domain%mpp_domain)
    call set_up_sponge_field(tmp, tv%T, nz, CSp)
    call read_data(filename, salin_var, tmp(:,:,:), domain=G%Domain%mpp_domain)
    call set_up_sponge_field(tmp, tv%S, nz, CSp)
  endif


end subroutine initialize_sponges_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_Open_Bdry_Conds(OBC, tv, G, param_file, tracer_Reg)
  type(ocean_OBC_type),       pointer    :: OBC
  type(thermo_var_ptrs),      intent(in) :: tv
  type(ocean_grid_type),      intent(in) :: G
  type(param_file_type),      intent(in) :: param_file
  type(tracer_registry_type), pointer    :: tracer_Reg
!   This subroutine sets the properties of flow at open boundary conditions.
! This particular example is for the DOME inflow describe in Legg et al. 2006.

! Arguments: OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (out)     tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

  logical :: any_OBC        ! Set to true if any points in this subdomain use
                            ! open boundary conditions.
  logical, pointer, dimension(:,:) :: &
    OBC_mask_u => NULL(), & ! These arrays are true at zonal or meridional
    OBC_mask_v => NULL()    ! velocity points that have prescribed open boundary
                            ! conditions.
  real, pointer, dimension(:,:,:) :: &
    OBC_T_u => NULL(), &    ! These arrays should be allocated and set to
    OBC_T_v => NULL(), &    ! specify the values of T and S that should come
    OBC_S_u => NULL(), &    ! in through u- and v- points through the open
    OBC_S_v => NULL()       ! boundary conditions, in C and psu.
  logical :: apply_OBC_u = .false., apply_OBC_v = .false.
  ! The following variables are used to set the target temperature and salinity.
  real :: T0(SZK_(G)), S0(SZK_(G))
  real :: pres(SZK_(G))      ! An array of the reference pressure in Pa.
  real :: drho_dT(SZK_(G))   ! Derivative of density with temperature in kg m-3 K-1.                              !
  real :: drho_dS(SZK_(G))   ! Derivative of density with salinity in kg m-3 PSU-1.                             !
  real :: rho_guess(SZK_(G)) ! Potential density at T0 & S0 in kg m-3.
  character(len=40) :: mod = "set_Open_Bdry_Conds"
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call get_param(param_file, mod, "APPLY_OBC_U", apply_OBC_u, default=.false.)
  call get_param(param_file, mod, "APPLY_OBC_V", apply_OBC_v, default=.false.)

  if (apply_OBC_u) then
    ! Determine where u points are applied.
    allocate(OBC_mask_u(IsdB:IedB,jsd:jed)) ; OBC_mask_u(:,:) = .false.
    any_OBC = .false.
    do j=jsd,jed ; do I=IsdB,IedB
    ! if (SOME_TEST_FOR_U_OPEN_BCS) then
    !   OBC_mask_u(I,j) = .true. ; any_OBC = .true.
    ! endif
    enddo ; enddo
    if (.not.any_OBC) then
      ! This processor has no u points at which open boundary conditions are
      ! to be applied.
      apply_OBC_u = .false.
      deallocate(OBC_mask_u)
    endif
  endif
  if (apply_OBC_v) then
    ! Determine where v points are applied.
    allocate(OBC_mask_v(isd:ied,JsdB:JedB)) ; OBC_mask_v(:,:) = .false.
    any_OBC = .false.
    do J=JsdB,JedB ; do i=isd,ied
    ! if (SOME_TEST_FOR_V_OPEN_BCS) then
    !   OBC_mask_v(i,J) = .true. ; any_OBC = .true.
    ! endif
    enddo ; enddo
    if (.not.any_OBC) then
      ! This processor has no v points at which open boundary conditions are
      ! to be applied.
      apply_OBC_v = .false.
      deallocate(OBC_mask_v)
    endif
  endif

  if (.not.(apply_OBC_u .or. apply_OBC_v)) return

  if (.not.associated(OBC)) allocate(OBC)

  if (apply_OBC_u) then
    OBC%apply_OBC_u = .true.
    OBC%OBC_mask_u => OBC_mask_u
    allocate(OBC%u(IsdB:IedB,jsd:jed,nz)) ; OBC%u(:,:,:) = 0.0
    allocate(OBC%uh(IsdB:IedB,jsd:jed,nz)) ; OBC%uh(:,:,:) = 0.0
    allocate(OBC%OBC_kind_u(IsdB:IedB,jsd:jed)) ; OBC%OBC_kind_u(:,:) = OBC_NONE
    do j=jsd,jed ; do I=IsdB,IedB
      if (OBC%OBC_mask_u(I,j)) OBC%OBC_kind_u(I,j) = OBC_SIMPLE
    enddo ; enddo
  endif
  if (apply_OBC_v) then
    OBC%apply_OBC_v = .true.
    OBC%OBC_mask_v => OBC_mask_v
    allocate(OBC%v(isd:ied,JsdB:JedB,nz)) ; OBC%v(:,:,:) = 0.0
    allocate(OBC%vh(isd:ied,JsdB:JedB,nz)) ; OBC%vh(:,:,:) = 0.0
    allocate(OBC%OBC_kind_v(isd:ied,JsdB:JedB)) ; OBC%OBC_kind_v(:,:) = OBC_NONE
    do J=JsdB,JedB ; do i=isd,ied
      if (OBC%OBC_mask_v(i,J)) OBC%OBC_kind_v(i,J) = OBC_SIMPLE
    enddo ; enddo
  endif

  if (apply_OBC_v) then
    do k=1,nz ; do J=Jsd,Jed ; do i=isd,ied
      if (OBC_mask_v(i,J)) then
        ! An appropriate expression for the meridional inflow velocities and
        ! transports should go here.
        OBC%vh(i,J,k) = 0.0 * G%m_to_H ; OBC%v(i,J,k) = 0.0
      else
        OBC%vh(i,J,k) = 0.0 ; OBC%v(i,J,k) = 0.0
      endif
    enddo ; enddo ; enddo
  endif

  if (apply_OBC_u) then
    do k=1,nz ; do j=jsd,jed ; do I=IsdB,IedB
      if (OBC_mask_u(I,j)) then
        ! An appropriate expression for the zonal inflow velocities and
        ! transports should go here.
        OBC%uh(I,j,k) = 0.0 * G%m_to_H ; OBC%u(I,j,k) = 0.0
      else
        OBC%uh(I,j,k) = 0.0 ; OBC%u(I,j,k) = 0.0
      endif
    enddo ; enddo ; enddo
  endif

  !   The inflow values of temperature and salinity also need to be set here if
  ! these variables are used.  The following code is just a naive example.
  if (apply_OBC_u .or. apply_OBC_v) then
    if (associated(tv%S)) then
      ! In this example, all S inflows have values of 35 psu.
      call add_tracer_OBC_values("S", tracer_Reg, OBC_inflow=35.0)
    endif
    if (associated(tv%T)) then
      ! In this example, the T values are set to be consistent with the layer
      ! target density and a salinity of 35 psu.  This code is taken from
      !  initialize_temp_sal.
      pres(:) = tv%P_Ref ; S0(:) = 35.0 ; T0(1) = 25.0
      call calculate_density(T0(1),S0(1),pres(1),rho_guess(1),tv%eqn_of_state)
      call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,1,tv%eqn_of_state)

      do k=1,nz ; T0(k) = T0(1) + (G%Rlay(k)-rho_guess(1)) / drho_dT(1) ; enddo
      do itt=1,6
        call calculate_density(T0,S0,pres,rho_guess,1,nz,tv%eqn_of_state)
        call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,nz,tv%eqn_of_state)
        do k=1,nz ; T0(k) = T0(k) + (G%Rlay(k)-rho_guess(k)) / drho_dT(k) ; enddo
      enddo

      if (apply_OBC_u) then
        allocate(OBC_T_u(IsdB:IedB,jsd:jed,nz))
        do k=1,nz ; do j=jsd,jed ; do I=IsdB,IedB
          OBC_T_u(I,j,k) = T0(k)
        enddo ; enddo ; enddo
      endif
      if (apply_OBC_v) then
        allocate(OBC_T_v(isd:ied,JsdB:JedB,nz))
        do k=1,nz ; do J=JsdB,JedB ; do i=isd,ied
          OBC_T_v(i,J,k) = T0(k)
        enddo ; enddo ; enddo
      endif
      call add_tracer_OBC_values("T", tracer_Reg, OBC_in_u=OBC_T_u, &
                                                  OBC_in_v=OBC_T_v)
    endif
  endif

end subroutine set_Open_Bdry_Conds
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_Flather_Bdry_Conds(OBC, tv, h, G, PF, tracer_Reg)
  type(ocean_OBC_type),                   pointer    :: OBC
  type(thermo_var_ptrs),                  intent(inout) :: tv
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout) :: h
  type(ocean_grid_type),                  intent(inout) :: G
  type(param_file_type),                  intent(in) :: PF
  type(tracer_registry_type),             pointer    :: tracer_Reg
!   This subroutine sets the initial definitions of the characteristic open boundary
!   conditions. Written by Mehmet Ilicak

! Arguments: OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (out)     tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in)      G - The ocean's grid structure.
!  (in)      PF - A structure indicating the open file to parse for
!                         model parameter values.

  logical :: any_OBC        ! Set to true if any points in this subdomain use
                            ! open boundary conditions.

  logical :: apply_OBC_u_flather_east = .false., apply_OBC_u_flather_west = .false.
  logical :: apply_OBC_v_flather_north = .false., apply_OBC_v_flather_south = .false.  
  logical :: read_OBC_eta = .false.
  logical :: read_OBC_uv = .false.
  logical :: read_OBC_TS = .false.

  integer :: isd_global, jsd_global
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: isd_off, jsd_off
  integer :: IsdB, IedB, JsdB, JedB
  integer :: east_boundary, west_boundary, north_boundary, south_boundary
  character(len=40)  :: mod = "set_Flather_Bdry_Conds" ! This subroutine's name.
  character(len=200) :: filename, OBC_file, inputdir ! Strings for file/path

  real :: temp_u(G%domain%niglobal+1,G%domain%njglobal)
  real :: temp_v(G%domain%niglobal,G%domain%njglobal+1)

  real, pointer, dimension(:,:,:) :: &
    OBC_T_u => NULL(), &    ! These arrays should be allocated and set to
    OBC_T_v => NULL(), &    ! specify the values of T and S that should come
    OBC_S_u => NULL(), & 
    OBC_S_v => NULL()     

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  
  isd_global = G%isd_global
  jsd_global = G%jsd_global
  call get_param(PF, mod, "APPLY_OBC_U_FLATHER_EAST", apply_OBC_u_flather_east,&
                 "Apply a Flather open boundary condition on the eastern \n"//&
                 "side of the global domain", default=.false.)
  call get_param(PF, mod, "APPLY_OBC_U_FLATHER_WEST", apply_OBC_u_flather_west,&
                 "Apply a Flather open boundary condition on the western \n"//&
                 "side of the global domain", default=.false.)
  call get_param(PF, mod, "APPLY_OBC_V_FLATHER_NORTH", apply_OBC_v_flather_north,&
                 "Apply a Flather open boundary condition on the northern \n"//&
                 "side of the global domain", default=.false.)
  call get_param(PF, mod, "APPLY_OBC_V_FLATHER_SOUTH", apply_OBC_v_flather_south,&
                 "Apply a Flather open boundary condition on the southern \n"//&
                 "side of the global domain", default=.false.)

  if (.not.(apply_OBC_u_flather_east  .or. apply_OBC_u_flather_west .or. &            
            apply_OBC_v_flather_north .or. apply_OBC_v_flather_south)) return
  
  if (.not.associated(OBC)) allocate(OBC)
       
  OBC%apply_OBC_u_flather_east = apply_OBC_u_flather_east
  OBC%apply_OBC_u_flather_west = apply_OBC_u_flather_west 
  OBC%apply_OBC_v_flather_north = apply_OBC_v_flather_north 
  OBC%apply_OBC_v_flather_south = apply_OBC_v_flather_south

  call get_param(PF, mod, "READ_OBC_UV", read_OBC_uv, &
                 "If true, read the values for the velocity open boundary \n"//&
                 "conditions from the file specified by OBC_FILE.", &
                 default=.false.)
  call get_param(PF, mod, "READ_OBC_ETA", read_OBC_eta, &
                 "If true, read the values for the sea surface height \n"//&
                 "open boundary conditions from the file specified by \n"//&
                 "OBC_FILE.", default=.false.)
  call get_param(PF, mod, "READ_OBC_TS", read_OBC_TS, &
                 "If true, read the values for the temperature and \n"//&
                 "salinity open boundary conditions from the file \n"//&
                 "specified by OBC_FILE.", default=.false.)
  if (read_OBC_uv .or. read_OBC_eta .or. read_OBC_TS) then
    call get_param(PF, mod, "OBC_FILE", OBC_file, &
                 "The file from which the appropriate open boundary \n"//&
                 "condition values are read.", default="MOM_OBC_FILE.nc")
    call get_param(PF, mod, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    filename = trim(inputdir)//trim(OBC_file)
    call log_param(PF, mod, "INPUTDIR/OBC_FILE", filename)
  endif

  if (G%symmetric) then
    east_boundary = G%ieg
    west_boundary = G%isg-1
    north_boundary = G%jeg
    south_boundary = G%jsg-1
  else
    ! I am not entirely sure that this works properly. -RWH
    east_boundary = G%ieg-1
    west_boundary = G%isg
    north_boundary = G%jeg-1
    south_boundary = G%jsg
  endif

  if (.not.associated(OBC%OBC_mask_u)) then
    allocate(OBC%OBC_mask_u(IsdB:IedB,jsd:jed)) ; OBC%OBC_mask_u(:,:) = .false.
  endif
  if (.not.associated(OBC%OBC_kind_u)) then
    allocate(OBC%OBC_kind_u(IsdB:IedB,jsd:jed)) ; OBC%OBC_kind_u(:,:) = OBC_NONE
  endif
  if (.not.associated(OBC%OBC_mask_v)) then
    allocate(OBC%OBC_mask_v(isd:ied,JsdB:JedB)) ; OBC%OBC_mask_v(:,:) = .false.
  endif
  if (.not.associated(OBC%OBC_kind_v)) then
    allocate(OBC%OBC_kind_v(isd:ied,JsdB:JedB)) ; OBC%OBC_kind_v(:,:) = OBC_NONE
  endif

  if (.not.associated(OBC%vbt_outer)) then
    allocate(OBC%vbt_outer(isd:ied,JsdB:JedB)) ; OBC%vbt_outer(:,:) = 0.0
  endif

  if (.not.associated(OBC%ubt_outer)) then
    allocate(OBC%ubt_outer(IsdB:IedB,jsd:jed)) ; OBC%ubt_outer(:,:) = 0.0
  endif

  if (.not.associated(OBC%eta_outer_u)) then
    allocate(OBC%eta_outer_u(IsdB:IedB,jsd:jed)) ; OBC%eta_outer_u(:,:) = 0.0
  endif

  if (.not.associated(OBC%eta_outer_v)) then
    allocate(OBC%eta_outer_v(isd:ied,JsdB:JedB)) ; OBC%eta_outer_v(:,:) = 0.0
  endif
  
  if (read_OBC_uv) then
    call read_data(filename, 'ubt', OBC%ubt_outer, &
                   domain=G%Domain%mpp_domain, position=EAST_FACE)
    call read_data(filename, 'vbt', OBC%vbt_outer, &
                   domain=G%Domain%mpp_domain, position=NORTH_FACE)
  endif

  if (read_OBC_eta) then
    call read_data(filename, 'eta_outer_u', OBC%eta_outer_u, &
                   domain=G%Domain%mpp_domain, position=EAST_FACE)
    call read_data(filename, 'eta_outer_v', OBC%eta_outer_v, &
                   domain=G%Domain%mpp_domain, position=NORTH_FACE)
  endif

  call pass_vector(OBC%eta_outer_u,OBC%eta_outer_v,G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
  call pass_vector(OBC%ubt_outer,OBC%vbt_outer,G%Domain)

  ! This code should be modified to allow OBCs to be applied anywhere.

  if (apply_OBC_u_flather_east) then
    ! Determine where u points are applied at east side 
    do j=jsd,jed ; do I=IsdB,IedB
      if ((I+isd_global-isd) == east_boundary) then !eastern side
        OBC%OBC_mask_u(I,j) = .true.
        OBC%OBC_kind_u(I,j) = OBC_FLATHER_E
        if ((i+1>isd) .and. (i+1<ied) .and. (J>JsdB) .and. (J<JedB)) then
          OBC%OBC_mask_v(i+1,J) = .true.
          if (OBC%OBC_kind_v(i+1,J) == OBC_NONE) OBC%OBC_kind_v(i+1,J) = OBC_FLATHER_E
        endif
        if ((i+1>isd) .and. (i+1<ied) .and. (J-1>JsdB) .and. (J-1<JedB)) then
          OBC%OBC_mask_v(i+1,J-1) = .true.
          if (OBC%OBC_kind_v(i+1,J-1) == OBC_NONE) OBC%OBC_kind_v(i+1,J-1) = OBC_FLATHER_E
        endif
      endif
    enddo ; enddo
  endif

  if (apply_OBC_u_flather_west) then
    ! Determine where u points are applied at west side 
    do j=jsd,jed ; do I=IsdB,IedB
      if ((I+isd_global-isd) == west_boundary) then !western side
        OBC%OBC_mask_u(I,j) = .true.
        OBC%OBC_kind_u(I,j) = OBC_FLATHER_W
        if ((i>isd) .and. (i<ied) .and. (J>JsdB) .and. (J<JedB)) then
          OBC%OBC_mask_v(i,J) = .true.
          if (OBC%OBC_kind_v(i,J) == OBC_NONE) OBC%OBC_kind_v(i,J) = OBC_FLATHER_W
        endif
        if ((i>isd) .and. (i<ied) .and. (J-1>JsdB) .and. (J-1<JedB)) then
          OBC%OBC_mask_v(i,J-1) = .true.
          if (OBC%OBC_kind_v(i,J-1) == OBC_NONE) OBC%OBC_kind_v(i,J-1) = OBC_FLATHER_W
        endif
      endif
    enddo ; enddo
  endif


  if (apply_OBC_v_flather_north) then
    ! Determine where v points are applied at north side 
    do J=JsdB,JedB ; do i=isd,ied
      if ((J+jsd_global-jsd) == north_boundary) then         !northern side
        OBC%OBC_mask_v(i,J) = .true.
        OBC%OBC_kind_v(i,J) = OBC_FLATHER_N
        if ((I>IsdB) .and. (I<IedB) .and. (j+1>jsd) .and. (j+1<jed)) then
          OBC%OBC_mask_u(I,j+1) = .true.
          if (OBC%OBC_kind_u(I,j+1) == OBC_NONE) OBC%OBC_kind_u(I,j+1) = OBC_FLATHER_N
        endif
        if ((I-1>IsdB) .and. (I-1<IedB) .and. (j+1>jsd) .and. (j+1<jed)) then
          OBC%OBC_mask_u(I-1,j+1) = .true.
          if (OBC%OBC_kind_u(I-1,j+1) == OBC_NONE) OBC%OBC_kind_u(I-1,j+1) = OBC_FLATHER_N
        endif
     endif
    enddo ; enddo
  endif
  
  if (apply_OBC_v_flather_south) then
    ! Determine where v points are applied at south side 
    do J=JsdB,JedB ; do i=isd,ied
      if ((J+jsd_global-jsd) == south_boundary) then         !southern side
        OBC%OBC_mask_v(i,J) = .true.
        OBC%OBC_kind_v(i,J) = OBC_FLATHER_S
        if ((I>IsdB) .and. (I<IedB) .and. (j>jsd) .and. (j<jed)) then
          OBC%OBC_mask_u(I,j) = .true.
          if (OBC%OBC_kind_u(I,j) == OBC_NONE) OBC%OBC_kind_u(I,j) = OBC_FLATHER_S
        endif
        if ((I-1>IsdB) .and. (I-1<IedB) .and. (j>jsd) .and. (j<jed)) then
          OBC%OBC_mask_u(I-1,j) = .true.
          if (OBC%OBC_kind_u(I-1,j) == OBC_NONE) OBC%OBC_kind_u(I-1,j) = OBC_FLATHER_S
        endif
      endif
    enddo ; enddo
  endif

  !   If there are no OBC points on this PE, there is no reason to keep the OBC
  ! type, and it could be deallocated.


  ! Define radiation coefficients r[xy]_old_[uvh] as needed.  For now, there are
  ! no radiation conditions applied to the thicknesses, since the thicknesses
  ! might not be physically motivated.  Instead, sponges should be used to
  ! enforce the near-boundary layer structure.
  if (apply_OBC_u_flather_west .or. apply_OBC_u_flather_east) then
    allocate(OBC%rx_old_u(IsdB:IedB,jsd:jed,nz)) ; OBC%rx_old_u(:,:,:) = 0.0
 !   allocate(OBC%rx_old_h(Isd:Ied,jsd:jed,nz))   ; OBC%rx_old_h(:,:,:) = 0.0
  endif
  if (apply_OBC_v_flather_south .or. apply_OBC_v_flather_north) then
    allocate(OBC%ry_old_v(isd:ied,JsdB:JedB,nz)) ; OBC%ry_old_v(:,:,:) = 0.0
 !   allocate(OBC%ry_old_h(isd:ied,Jsd:Jed,nz))   ; OBC%ry_old_h(:,:,:) = 0.0
  endif


  if (associated(tv%T)) then
    allocate(OBC_T_u(IsdB:IedB,jsd:jed,nz)) ; OBC_T_u(:,:,:) = 0.0
    allocate(OBC_S_u(IsdB:IedB,jsd:jed,nz)) ; OBC_S_u(:,:,:) = 0.0
    allocate(OBC_T_v(isd:ied,JsdB:JedB,nz)) ; OBC_T_v(:,:,:) = 0.0
    allocate(OBC_S_v(isd:ied,JsdB:JedB,nz)) ; OBC_S_v(:,:,:) = 0.0

    if (read_OBC_TS) then
      call read_data(filename, 'OBC_T_u', OBC_T_u, &
                     domain=G%Domain%mpp_domain, position=EAST_FACE)
      call read_data(filename, 'OBC_S_u', OBC_S_u, &
                     domain=G%Domain%mpp_domain, position=EAST_FACE)

      call read_data(filename, 'OBC_T_v', OBC_T_v, &
                     domain=G%Domain%mpp_domain, position=NORTH_FACE)
      call read_data(filename, 'OBC_S_v', OBC_S_v, &
                     domain=G%Domain%mpp_domain, position=NORTH_FACE)
    else
      call pass_var(tv%T, G%Domain)
      call pass_var(tv%S, G%Domain)
      do k=1,nz ; do j=js,je ; do I=is-1,ie
        if (OBC%OBC_mask_u(I,j)) then
          if (OBC%OBC_kind_u(I,j) == OBC_FLATHER_E) then
            OBC_T_u(I,j,k) = tv%T(i,j,k)
            OBC_S_u(I,j,k) = tv%S(i,j,k)
          elseif (OBC%OBC_kind_u(I,j) == OBC_FLATHER_W) then
            OBC_T_u(I,j,k) = tv%T(i+1,j,k)
            OBC_S_u(I,j,k) = tv%S(i+1,j,k)
          elseif (G%mask2dT(i,j) + G%mask2dT(i+1,j) > 0) then
            OBC_T_u(I,j,k) = (G%mask2dT(i,j)*tv%T(i,j,k) + G%mask2dT(i+1,j)*tv%T(i+1,j,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i+1,j))
            OBC_S_u(I,j,k) = (G%mask2dT(i,j)*tv%S(i,j,k) + G%mask2dT(i+1,j)*tv%S(i+1,j,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i+1,j))
          else ! This probably shouldn't happen or maybe it doesn't matter?
            OBC_T_u(I,j,k) = 0.5*(tv%T(i,j,k)+tv%T(i+1,j,k))
            OBC_S_u(I,j,k) = 0.5*(tv%S(i,j,k)+tv%S(i+1,j,k))
          endif
        else
          OBC_T_u(I,j,k) = 0.5*(tv%T(i,j,k)+tv%T(i+1,j,k))
          OBC_S_u(I,j,k) = 0.5*(tv%S(i,j,k)+tv%S(i+1,j,k))
        endif
      enddo; enddo ; enddo

      do k=1,nz ; do J=js-1,je ; do i=is,ie
        if (OBC%OBC_mask_v(i,J)) then
          if (OBC%OBC_kind_v(i,J) == OBC_FLATHER_N) then
            OBC_T_v(i,J,k) = tv%T(i,j,k)
            OBC_S_v(i,J,k) = tv%S(i,j,k)
          elseif (OBC%OBC_kind_v(i,J) == OBC_FLATHER_S) then
            OBC_T_v(i,J,k) = tv%T(i,j+1,k)
            OBC_S_v(i,J,k) = tv%S(i,j+1,k)
          elseif (G%mask2dT(i,j) + G%mask2dT(i,j+1) > 0) then
            OBC_T_v(i,J,k) = (G%mask2dT(i,j)*tv%T(i,j,k) + G%mask2dT(i,j+1)*tv%T(i,j+1,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i,j+1))
            OBC_S_v(i,J,k) = (G%mask2dT(i,j)*tv%S(i,j,k) + G%mask2dT(i,j+1)*tv%S(i,j+1,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i,j+1))
          else ! This probably shouldn't happen or maybe it doesn't matter?
            OBC_T_v(i,J,k) = 0.5*(tv%T(i,j,k)+tv%T(i,j+1,k))
            OBC_S_v(i,J,k) = 0.5*(tv%S(i,j,k)+tv%S(i,j+1,k))
          endif
        else
          OBC_T_v(i,J,k) = 0.5*(tv%T(i,j,k)+tv%T(i,j+1,k))
          OBC_S_v(i,J,k) = 0.5*(tv%S(i,j,k)+tv%S(i,j+1,k))
        endif
      enddo; enddo ; enddo
    endif

    call pass_vector(OBC_T_u, OBC_T_v, G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
    call pass_vector(OBC_S_u, OBC_S_v, G%Domain, To_All+SCALAR_PAIR, CGRID_NE)

    call add_tracer_OBC_values("T", tracer_Reg, OBC_in_u=OBC_T_u, &
                                                OBC_in_v=OBC_T_v)
    call add_tracer_OBC_values("S", tracer_Reg, OBC_in_u=OBC_S_u, &
                                                OBC_in_v=OBC_S_v)
    do k=1,nz ; do j=js,je ; do I=is-1,ie
      if (OBC%OBC_kind_u(I,j) == OBC_FLATHER_E) then
        tv%T(i+1,j,k) = tv%T(i,j,k) ; tv%S(i+1,j,k) = tv%S(i,j,k)
      elseif (OBC%OBC_kind_u(I,j) == OBC_FLATHER_W) then
        tv%T(i,j,k) = tv%T(i+1,j,k) ; tv%S(i,j,k) = tv%S(i+1,j,k)
      endif
    enddo ; enddo ; enddo
    do k=1,nz ; do J=js-1,je ; do i=is,ie
      if (OBC%OBC_kind_v(i,J) == OBC_FLATHER_N) then
        tv%T(i,j+1,k) = tv%T(i,j,k) ; tv%S(i,j+1,k) = tv%S(i,j,k)
      elseif (OBC%OBC_kind_v(i,J) == OBC_FLATHER_S) then
        tv%T(i,j,k) = tv%T(i,j+1,k) ; tv%S(i,j,k) = tv%S(i,j+1,k)
      endif
    enddo ; enddo ; enddo
  endif

  do k=1,nz ; do j=js,je ; do I=is-1,ie
    if (OBC%OBC_kind_u(I,j) == OBC_FLATHER_E) h(i+1,j,k) = h(i,j,k)
    if (OBC%OBC_kind_u(I,j) == OBC_FLATHER_W) h(i,j,k) = h(i+1,j,k)
  enddo ; enddo ; enddo
  do k=1,nz ; do J=js-1,je ; do i=is,ie
    if (OBC%OBC_kind_v(i,J) == OBC_FLATHER_N) h(i,j+1,k) = h(i,j,k)
    if (OBC%OBC_kind_v(i,J) == OBC_FLATHER_S) h(i,j,k) = h(i,j+1,k)
  enddo ; enddo ; enddo

end subroutine set_Flather_Bdry_Conds   
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine reset_face_lengths_named(G, param_file, name)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
  character(len=*),      intent(in)    :: name
!   This subroutine sets the open face lengths at selected points to restrict
! passages to their observed widths.

! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      name - The name for the set of face lengths.
  character(len=256) :: mesg    ! Message for error messages.
  real    :: dx_2 = -1.0, dy_2 = -1.0
  real    :: pi_180
  integer :: option = -1
  integer :: i, j, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  pi_180 = (4.0*atan(1.0))/180.0

  select case ( trim(name) )
    case ("global_1deg")    ; option = 1 ; dx_2 = 0.5*1.0
    case default ; call MOM_error(FATAL, "reset_face_lengths_named: "//&
      "Unrecognized channel configuration name "//trim(name))
  end select

  if (option==1) then ! 1-degree settings.
    do j=jsd,jed ; do I=IsdB,IedB  ! Change any u-face lengths within this loop.
      dy_2 = dx_2 * G%dyCu(I,j)*G%IdxCu(I,j) * cos(pi_180 * G%geoLatCu(I,j))

      if ((abs(G%geoLatCu(I,j)-35.5) < dy_2) .and. (G%geoLonCu(I,j) < -4.5) .and. &
          (G%geoLonCu(I,j) > -6.5)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*12000.0   ! Gibraltar

      if ((abs(G%geoLatCu(I,j)-12.5) < dy_2) .and. (abs(G%geoLonCu(I,j)-43.0) < dx_2)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*10000.0   ! Red Sea

      if ((abs(G%geoLatCu(i,j)-40.5) < dy_2) .and. (abs(G%geoLonCu(i,j)-26.0) < dx_2)) &
        G%dy_Cu(i,j) = G%mask2dCu(i,j)*5000.0   ! Dardanelles

      if ((abs(G%geoLatCu(I,j)-41.5) < dy_2) .and. (abs(G%geoLonCu(I,j)+220.0) < dx_2)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*35000.0   ! Tsugaru strait at 140.0e

      if ((abs(G%geoLatCu(I,j)-45.5) < dy_2) .and. (abs(G%geoLonCu(I,j)+217.5) < 0.9)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*15000.0   ! Betw Hokkaido and Sakhalin at 217&218 = 142e


      ! Greater care needs to be taken in the tripolar region.
      if ((abs(G%geoLatCu(I,j)-80.84) < 0.2) .and. (abs(G%geoLonCu(I,j)+64.9) < 0.8)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*38000.0   ! Smith Sound in Canadian Arch - tripolar region

    enddo ; enddo

    do J=JsdB,JedB ; do i=isd,ied  ! Change any v-face lengths within this loop.
      dy_2 = dx_2 * G%dyCv(i,J)*G%IdxCv(i,J) * cos(pi_180 * G%geoLatCv(i,J))
      if ((abs(G%geoLatCv(i,J)-41.0) < dy_2) .and. (abs(G%geoLonCv(i,J)-28.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*2500.0   ! Bosporus - should be 1000.0 m wide.

      if ((abs(G%geoLatCv(i,J)-13.0) < dy_2) .and. (abs(G%geoLonCv(i,J)-42.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*10000.0   ! Red Sea

      if ((abs(G%geoLatCv(i,J)+2.8) < 0.8) .and. (abs(G%geoLonCv(i,J)+241.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*40000.0   ! Makassar Straits at 241.5 W = 118.5 E

      if ((abs(G%geoLatCv(i,J)-0.56) < 0.5) .and. (abs(G%geoLonCv(i,J)+240.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*80000.0   ! entry to Makassar Straits at 240.5 W = 119.5 E

      if ((abs(G%geoLatCv(i,J)-0.19) < 0.5) .and. (abs(G%geoLonCv(i,J)+230.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*25000.0   ! Channel betw N Guinea and Halmahara 230.5 W = 129.5 E

      if ((abs(G%geoLatCv(i,J)-0.19) < 0.5) .and. (abs(G%geoLonCv(i,J)+229.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*25000.0   ! Channel betw N Guinea and Halmahara 229.5 W = 130.5 E

      if ((abs(G%geoLatCv(i,J)-0.0) < 0.25) .and. (abs(G%geoLonCv(i,J)+228.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*25000.0   ! Channel betw N Guinea and Halmahara 228.5 W = 131.5 E

      if ((abs(G%geoLatCv(i,J)+8.5) < 0.5) .and. (abs(G%geoLonCv(i,J)+244.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*20000.0   ! Lombok Straits at 244.5 W = 115.5 E

      if ((abs(G%geoLatCv(i,J)+8.5) < 0.5) .and. (abs(G%geoLonCv(i,J)+235.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*20000.0   ! Timor Straits at 235.5 W = 124.5 E

      if ((abs(G%geoLatCv(i,J)-52.5) < dy_2) .and. (abs(G%geoLonCv(i,J)+218.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*2500.0    ! Russia and Sakhalin Straits at 218.5 W = 141.5 E

      ! Greater care needs to be taken in the tripolar region.
      if ((abs(G%geoLatCv(i,J)-76.8) < 0.06) .and. (abs(G%geoLonCv(i,J)+88.7) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*8400.0    ! Jones Sound in Canadian Arch - tripolar region

    enddo ; enddo
  endif

  ! These checks apply regardless of the chosen option.

  do j=jsd,jed ; do I=IsdB,IedB
    if (G%dy_Cu(I,j) > G%dyCu(I,j)) then
      write(mesg,'("dy_Cu of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4," at lon/lat of ", ES11.4, ES11.4)') &
                   G%dy_Cu(I,j), G%dyCu(I,j), G%dy_Cu(I,j)-G%dyCu(I,j), &
                   G%geoLonCu(I,j), G%geoLatCu(I,j)
      call MOM_error(FATAL,"reset_face_lengths_named "//mesg)
    endif
    G%areaCu(I,j) = G%dxCu(I,j)*G%dy_Cu(I,j)
    G%IareaCu(I,j) = 0.0
    if (G%areaCu(I,j) > 0.0) G%IareaCu(I,j) = G%mask2dCu(I,j) / G%areaCu(I,j)
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    if (G%dx_Cv(i,J) > G%dxCv(i,J)) then
      write(mesg,'("dx_Cv of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4, " at lon/lat of ", ES11.4, ES11.4)') &
                   G%dx_Cv(i,J), G%dxCv(i,J), G%dx_Cv(i,J)-G%dxCv(i,J), &
                   G%geoLonCv(i,J), G%geoLatCv(i,J)

      call MOM_error(FATAL,"reset_face_lengths_named "//mesg)
    endif
    G%areaCv(i,J) = G%dyCv(i,J)*G%dx_Cv(i,J)
    G%IareaCv(i,J) = 0.0
    if (G%areaCv(i,J) > 0.0) G%IareaCv(i,J) = G%mask2dCv(i,J) / G%areaCv(i,J)
  enddo ; enddo

end subroutine reset_face_lengths_named
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------

subroutine reset_face_lengths_file(G, param_file)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
!   This subroutine sets the open face lengths at selected points to restrict
! passages to their observed widths.

! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
  character(len=40)  :: mod = "reset_face_lengths_file" ! This subroutine's name.
  character(len=256) :: mesg    ! Message for error messages.
  character(len=200) :: filename, chan_file, inputdir ! Strings for file/path
  integer :: i, j, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  ! These checks apply regardless of the chosen option.

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, mod, "CHANNEL_WIDTH_FILE", chan_file, &
                 "The file from which the list of narrowed channels is read.", &
                 default="ocean_geometry.nc")
  call get_param(param_file,  mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  filename = trim(inputdir)//trim(chan_file)
  call log_param(param_file, mod, "INPUTDIR/CHANNEL_WIDTH_FILE", filename)

  if (is_root_pe()) then ; if (.not.file_exists(filename)) &
    call MOM_error(FATAL," reset_face_lengths_file: Unable to open "//&
                           trim(filename))
  endif

  call read_data(filename,"dyCuo",G%dy_Cu,domain=G%Domain%mpp_domain)
  call read_data(filename,"dxCvo",G%dx_Cv,domain=G%Domain%mpp_domain)
  call pass_vector(G%dy_Cu, G%dx_Cv, G%Domain, To_All+SCALAR_PAIR, CGRID_NE)

  do j=jsd,jed ; do I=IsdB,IedB
    if (G%dy_Cu(I,j) > G%dyCu(I,j)) then
      write(mesg,'("dy_Cu of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4," at lon/lat of ", ES11.4, ES11.4)') &
                   G%dy_Cu(I,j), G%dyCu(I,j), G%dy_Cu(I,j)-G%dyCu(I,j), &
                   G%geoLonCu(I,j), G%geoLatCu(I,j)
      call MOM_error(FATAL,"reset_face_lengths_file "//mesg)
    endif
    G%areaCu(I,j) = G%dxCu(I,j)*G%dy_Cu(I,j)
    G%IareaCu(I,j) = 0.0
    if (G%areaCu(I,j) > 0.0) G%IareaCu(I,j) = G%mask2dCu(I,j) / G%areaCu(I,j)
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    if (G%dx_Cv(i,J) > G%dxCv(i,J)) then
      write(mesg,'("dx_Cv of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4, " at lon/lat of ", ES11.4, ES11.4)') &
                   G%dx_Cv(i,J), G%dxCv(i,J), G%dx_Cv(i,J)-G%dxCv(i,J), &
                   G%geoLonCv(i,J), G%geoLatCv(i,J)

      call MOM_error(FATAL,"reset_face_lengths_file "//mesg)
    endif
    G%areaCv(i,J) = G%dyCv(i,J)*G%dx_Cv(i,J)
    G%IareaCv(i,J) = 0.0
    if (G%areaCv(i,J) > 0.0) G%IareaCv(i,J) = G%mask2dCv(i,J) / G%areaCv(i,J)
  enddo ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine reset_face_lengths_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine reset_face_lengths_list(G, param_file)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
!   This subroutine sets the open face lengths at selected points to restrict
! passages to their observed widths.

! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
  character(len=120), pointer, dimension(:) :: lines => NULL()
  character(len=120) :: line
  character(len=200) :: filename, chan_file, inputdir ! Strings for file/path
  character(len=40)  :: mod = "reset_face_lengths_list" ! This subroutine's name.
  real, pointer, dimension(:,:) :: &
    u_lat => NULL(), u_lon => NULL(), v_lat => NULL(), v_lon => NULL()
  real, pointer, dimension(:) :: &
    u_width => NULL(), v_width => NULL()
  real    :: lat, lon     ! The latitude and longitude of a point.
  real    :: lon_p, lon_m ! The longitude of a point shifted by 360 degrees.
  logical :: check_360    ! If true, check for longitudes that are shifted by
                          ! +/- 360 degrees from the specified range of values.
  logical :: found_u, found_v
  logical :: unit_in_use
  integer :: ios, iounit, isu, isv
  integer :: last, num_lines, nl_read, ln, npt, u_pt, v_pt
  integer :: i, j, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")

  call get_param(param_file, mod, "CHANNEL_LIST_FILE", chan_file, &
                 "The file from which the list of narrowed channels is read.", &
                 default="MOM_channel_list")
  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  filename = trim(inputdir)//trim(chan_file)
  call log_param(param_file, mod, "INPUTDIR/CHANNEL_LIST_FILE", filename)
  call get_param(param_file, mod, "CHANNEL_LIST_360_LON_CHECK", check_360, &
                 "If true, the channel configuration list works for any \n"//&
                 "longitudes in the range of -360 to 360.", default=.true.)

  if (is_root_pe()) then
    ! Open the input file.
    if (.not.file_exists(filename)) call MOM_error(FATAL, &
        " reset_face_lengths_list: Unable to open "//trim(filename))

    ! Find an unused unit number.
    do iounit=10,512
      INQUIRE(iounit,OPENED=unit_in_use) ; if (.not.unit_in_use) exit
    enddo
    if (iounit >= 512) call MOM_error(FATAL, &
        "reset_face_lengths_list: No unused file unit could be found.")

    ! Open the parameter file.
    open(iounit, file=trim(filename), access='SEQUENTIAL', &
         form='FORMATTED', action='READ', position='REWIND', iostat=ios)
    if (ios /= 0) call MOM_error(FATAL, &
            "reset_face_lengths_list: Error opening "//trim(filename))

    ! Count the number of u_width and v_width entries.
    call read_face_length_list(iounit, filename, num_lines, lines)
  endif

  ! Broadcast the number of lines and allocate the required space.
  call broadcast(num_lines, root_PE())
  u_pt = 0 ; v_pt = 0
  if (num_lines > 0) then
    allocate (lines(num_lines))
    if (num_lines > 0) then
      allocate(u_lat(2,num_lines)) ; u_lat(:,:) = -1e34
      allocate(u_lon(2,num_lines)) ; u_lon(:,:) = -1e34
      allocate(u_width(num_lines)) ; u_width(:) = -1e34

      allocate(v_lat(2,num_lines)) ; v_lat(:,:) = -1e34
      allocate(v_lon(2,num_lines)) ; v_lon(:,:) = -1e34
      allocate(v_width(num_lines)) ; v_width(:) = -1e34
    endif

    ! Actually read the lines.
    if (is_root_pe()) then
      call read_face_length_list(iounit, filename, nl_read, lines)
      if (nl_read /= num_lines) &
        call MOM_error(FATAL, 'reset_face_lengths_list : Found different '// &
                  'number of valid lines on second reading of '//trim(filename))
      close(iounit) ; iounit = -1
    endif

    ! Broadcast the lines.
    call broadcast(lines, 120, root_PE())

    ! Populate the u_width, etc., data.
    do ln=1,num_lines
      line = lines(ln)
      ! Detect keywords
      found_u = .false.; found_v = .false.
      isu = index(uppercase(line), "U_WIDTH" ); if (isu > 0) found_u = .true.
      isv = index(uppercase(line), "V_WIDTH" ); if (isv > 0) found_v = .true.
      
      ! Store and check the relevant values.
      if (found_u) then
        u_pt = u_pt + 1
        read(line(isu+8:),*) u_lon(1:2,u_pt), u_lat(1:2,u_pt), u_width(u_pt) 
        if (is_root_PE()) then
          if (check_360) then
            if ((abs(u_lon(1,u_pt)) > 360.0) .or. (abs(u_lon(2,u_pt)) > 360.0)) &
              call MOM_error(WARNING, "reset_face_lengths_list : Out-of-bounds "//&
                 "u-longitude found when reading line "//trim(line)//" from file "//&
                 trim(filename))
            if ((abs(u_lat(1,u_pt)) > 180.0) .or. (abs(u_lat(2,u_pt)) > 180.0)) &
              call MOM_error(WARNING, "reset_face_lengths_list : Out-of-bounds "//&
                 "u-latitude found when reading line "//trim(line)//" from file "//&
                 trim(filename))
          endif
          if (u_lat(1,u_pt) > u_lat(2,u_pt)) &
            call MOM_error(WARNING, "reset_face_lengths_list : Out-of-order "//&
               "u-face latitudes found when reading line "//trim(line)//" from file "//&
               trim(filename))
          if (u_lon(1,u_pt) > u_lon(2,u_pt)) &
            call MOM_error(WARNING, "reset_face_lengths_list : Out-of-order "//&
               "u-face longitudes found when reading line "//trim(line)//" from file "//&
               trim(filename))
          if (u_width(u_pt) < 0.0) &
            call MOM_error(WARNING, "reset_face_lengths_list : Negative "//&
               "u-width found when reading line "//trim(line)//" from file "//&
               trim(filename))
        endif
      elseif (found_v) then
        v_pt = v_pt + 1
        read(line(isv+8:),*) v_lon(1:2,v_pt), v_lat(1:2,v_pt), v_width(v_pt)
        if (is_root_PE()) then
          if (check_360) then
            if ((abs(v_lon(1,v_pt)) > 360.0) .or. (abs(v_lon(2,v_pt)) > 360.0)) &
              call MOM_error(WARNING, "reset_face_lengths_list : Out-of-bounds "//&
                 "v-longitude found when reading line "//trim(line)//" from file "//&
                 trim(filename))
            if ((abs(v_lat(1,v_pt)) > 180.0) .or. (abs(v_lat(2,v_pt)) > 180.0)) &
              call MOM_error(WARNING, "reset_face_lengths_list : Out-of-bounds "//&
                 "v-latitude found when reading line "//trim(line)//" from file "//&
                 trim(filename))
          endif
          if (v_lat(1,v_pt) > v_lat(2,v_pt)) &
            call MOM_error(WARNING, "reset_face_lengths_list : Out-of-order "//&
               "v-face latitudes found when reading line "//trim(line)//" from file "//&
               trim(filename))
          if (v_lon(1,v_pt) > v_lon(2,v_pt)) &
            call MOM_error(WARNING, "reset_face_lengths_list : Out-of-order "//&
               "v-face longitudes found when reading line "//trim(line)//" from file "//&
               trim(filename))
          if (v_width(v_pt) < 0.0) &
            call MOM_error(WARNING, "reset_face_lengths_list : Negative "//&
               "v-width found when reading line "//trim(line)//" from file "//&
               trim(filename))
        endif
      endif
    enddo
    
    deallocate(lines)
  endif

  do j=jsd,jed ; do I=IsdB,IedB
    lat = G%geoLatCu(I,j) ; lon = G%geoLonCu(I,j)
    if (check_360) then ; lon_p = lon+360.0 ; lon_m = lon-360.0
    else ; lon_p = lon ; lon_m = lon ; endif

    do npt=1,u_pt
      if (((lat >= u_lat(1,npt)) .and. (lat <= u_lat(2,npt))) .and. &
          (((lon >= u_lon(1,npt)) .and. (lon <= u_lon(2,npt))) .or. &
           ((lon_p >= u_lon(1,npt)) .and. (lon_p <= u_lon(2,npt))) .or. &
           ((lon_m >= u_lon(1,npt)) .and. (lon_m <= u_lon(2,npt)))) ) &
  
      G%dy_Cu(I,j) = G%mask2dCu(I,j) * min(G%dyCu(I,j), max(u_width(npt), 0.0))
    enddo

    G%areaCu(I,j) = G%dxCu(I,j)*G%dy_Cu(I,j)
    G%IareaCu(I,j) = 0.0
    if (G%areaCu(I,j) > 0.0) G%IareaCu(I,j) = G%mask2dCu(I,j) / G%areaCu(I,j)
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    lat = G%geoLatCv(i,J) ; lon = G%geoLonCv(i,J)
    if (check_360) then ; lon_p = lon+360.0 ; lon_m = lon-360.0
    else ; lon_p = lon ; lon_m = lon ; endif

    do npt=1,v_pt
      if (((lat >= v_lat(1,npt)) .and. (lat <= v_lat(2,npt))) .and. &
          (((lon >= v_lon(1,npt)) .and. (lon <= v_lon(2,npt))) .or. &
           ((lon_p >= v_lon(1,npt)) .and. (lon_p <= v_lon(2,npt))) .or. &
           ((lon_m >= v_lon(1,npt)) .and. (lon_m <= v_lon(2,npt)))) ) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J) * min(G%dxCv(i,J), max(v_width(npt), 0.0))
    enddo

    G%areaCv(i,J) = G%dyCv(i,J)*G%dx_Cv(i,J)
    G%IareaCv(i,J) = 0.0
    if (G%areaCv(i,J) > 0.0) G%IareaCv(i,J) = G%mask2dCv(i,J) / G%areaCv(i,J)
  enddo ; enddo

  if (num_lines > 0) then
    deallocate(u_lat) ; deallocate(u_lon) ; deallocate(u_width)
    deallocate(v_lat) ; deallocate(v_lon) ; deallocate(v_width)
  endif

  call callTree_leave(trim(mod)//'()')
end subroutine reset_face_lengths_list

subroutine read_face_length_list(iounit, filename, num_lines, lines)
  integer,                          intent(in)  :: iounit
  character(len=*),                 intent(in)  :: filename
  integer,                          intent(out) :: num_lines
  character(len=120), dimension(:), pointer     :: lines

  !   This subroutine reads and counts the non-blank lines in the face length
  ! list file, after removing comments.
  character(len=120) :: line, line_up
  logical :: found_u, found_v
  integer :: isu, isv, icom, verbose
  integer :: last
  
  num_lines = 0

  if (iounit <= 0) return
  rewind(iounit)
  do while(.true.)
    read(iounit, '(a)', end=8, err=9) line
    last = len_trim(line)
    ! Eliminate either F90 or C comments from the line.
    icom = index(line(:last), "!") ; if (icom > 0) last = icom-1
    icom = index(line(:last), "/*") ; if (icom > 0) last = icom-1
    if (last < 1) cycle

    ! Detect keywords
    line_up = uppercase(line)
    found_u = .false.; found_v = .false.
    isu = index(line_up(:last), "U_WIDTH" ); if (isu > 0) found_u = .true.
    isv = index(line_up(:last), "V_WIDTH" ); if (isv > 0) found_v = .true.

    if (found_u .and. found_v) call MOM_error(FATAL, &
      "read_face_length_list : both U_WIDTH and V_WIDTH found when "//&
      "reading the line "//trim(line(:last))//" in file "//trim(filename))
    if (found_u .or. found_v) then
      num_lines = num_lines + 1
      if (associated(lines)) then
        lines(num_lines) = line(1:last)
      endif
    endif
  enddo ! while (.true.)

8 continue
  return

9 call MOM_error(FATAL, "read_face_length_list : "//&
                  "Error while reading file "//trim(filename))
 
end subroutine read_face_length_list

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
subroutine write_ocean_geometry_file(G, param_file, directory)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
  character(len=*),      intent(in)    :: directory
!   This subroutine writes out a file containing all of the ocean geometry
! and grid data uses by the MOM ocean model.
! Arguments: G - The ocean's grid structure.  Effectively intent in.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      directory - The directory into which to place the file.
  character(len=120) :: filepath
  character(len=40)  :: mod = "write_ocean_geometry_file"
  integer, parameter :: nFlds=23
  type(vardesc) :: vars(nFlds)
  type(fieldtype) :: fields(nFlds)
  integer :: unit
  integer :: file_threading
  integer :: nFlds_used
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  logical :: multiple_files
  real :: out_h(SZI_(G),SZJ_(G))
  real :: out_u(SZIB_(G),SZJ_(G))
  real :: out_v(SZI_(G),SZJB_(G))
  real :: out_q(SZIB_(G),SZJB_(G))
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

!   vardesc is a structure defined in MOM_io.F90.  The elements of
! this structure, in order, are:
! (1) the variable name for the NetCDF file
! (2) the variable's long name
! (3) a character indicating the  horizontal grid, which may be '1' (column),
!     'h', 'q', 'u', or 'v', for the corresponding C-grid variable
! (4) a character indicating the vertical grid, which may be 'L' (layer),
!     'i' (interface), or '1' (no vertical location)
! (5) a character indicating the time levels of the field, which may be
!    's' (snap-shot), 'p' (periodic), or '1' (no time variation)
! (6) the variable's units
  vars(1) = vardesc("geolatb","latitude at corner (Bu) points",'q','1','1',"degree")
  vars(2) = vardesc("geolonb","longitude at corner (Bu) points",'q','1','1',"degree")
  vars(3) = vardesc("geolat", "latitude at tracer (T) points", 'h','1','1',"degree")
  vars(4) = vardesc("geolon","longitude at tracer (T) points",'h','1','1',"degree")
  vars(5) = vardesc("D","Basin Depth",'h','1','1',"meter")
  vars(6) = vardesc("f","Coriolis Parameter",'q','1','1',"second-1")
  vars(7) = vardesc("dxCv","Zonal grid spacing at v points",'v','1','1',"m")
  vars(8) = vardesc("dyCu","Meridional grid spacing at u points",'u','1','1',"m")
  vars(9) = vardesc("dxCu","Zonal grid spacing at u points",'u','1','1',"m")
  vars(10)= vardesc("dyCv","Meridional grid spacing at v points",'v','1','1',"m")
  vars(11)= vardesc("dxT","Zonal grid spacing at h points",'h','1','1',"m")
  vars(12)= vardesc("dyT","Meridional grid spacing at h points",'h','1','1',"m")
  vars(13)= vardesc("dxBu","Zonal grid spacing at q points",'q','1','1',"m")
  vars(14)= vardesc("dyBu","Meridional grid spacing at q points",'q','1','1',"m")
  vars(15)= vardesc("Ah","Area of h cells",'h','1','1',"m2")
  vars(16)= vardesc("Aq","Area of q cells",'q','1','1',"m2")

  vars(17)= vardesc("dxCvo","Open zonal grid spacing at v points",'v','1','1',"m")
  vars(18)= vardesc("dyCuo","Open meridional grid spacing at u points",'u','1','1',"m")
  vars(19)= vardesc("wet", "land or ocean?", 'h','1','1',"none")

  vars(20) = vardesc("Dblock_u","Blocked depth at u points",'u','1','1',"meter")
  vars(21) = vardesc("Dopen_u","Open depth at u points",'u','1','1',"meter")
  vars(22) = vardesc("Dblock_v","Blocked depth at v points",'v','1','1',"meter")
  vars(23) = vardesc("Dopen_v","Open depth at v points",'v','1','1',"meter")

  nFlds_used = 19 ; if (G%bathymetry_at_vel) nFlds_used = 23

  filepath = trim(directory) // "ocean_geometry"

  out_h(:,:) = 0.0
  out_u(:,:) = 0.0
  out_v(:,:) = 0.0
  out_q(:,:) = 0.0

  call get_param(param_file, mod, "PARALLEL_RESTARTFILES", multiple_files, &
                 "If true, each processor writes its own restart file, \n"//&
                 "otherwise a single restart file is generated", &
                 default=.false.)
  file_threading = SINGLE_FILE
  if (multiple_files) file_threading = MULTIPLE

  call create_file(unit, trim(filepath), vars, nFlds_used, G, fields, file_threading)

  do J=Jsq,Jeq; do I=Isq,Ieq; out_q(I,J) = G%geoLatBu(I,J); enddo; enddo
  call write_field(unit, fields(1), G%Domain%mpp_domain, out_q)
  do J=Jsq,Jeq; do I=Isq,Ieq; out_q(I,J) = G%geoLonBu(I,J); enddo; enddo
  call write_field(unit, fields(2), G%Domain%mpp_domain, out_q)
  call write_field(unit, fields(3), G%Domain%mpp_domain, G%geoLatT)
  call write_field(unit, fields(4), G%Domain%mpp_domain, G%geoLonT)

  call write_field(unit, fields(5), G%Domain%mpp_domain, G%bathyT)
  call write_field(unit, fields(6), G%Domain%mpp_domain, G%CoriolisBu)

  do J=Jsq,Jeq; do i=is,ie; out_v(i,J) = G%dxCv(i,J); enddo; enddo
  call write_field(unit, fields(7), G%Domain%mpp_domain, out_v)
  do j=js,je; do I=Isq,Ieq; out_u(I,j) = G%dyCu(I,j); enddo; enddo
  call write_field(unit, fields(8), G%Domain%mpp_domain, out_u)

  do J=Jsq,Jeq; do i=is,ie; out_u(i,J) = G%dxCu(i,J); enddo; enddo
  call write_field(unit, fields(9), G%Domain%mpp_domain, out_u)
  do j=js,je; do I=Isq,Ieq; out_v(I,j) = G%dyCv(I,j); enddo; enddo
  call write_field(unit, fields(10), G%Domain%mpp_domain, out_v)

  do J=Jsq,Jeq; do i=is,ie; out_h(i,J) = G%dxT(i,J); enddo; enddo
  call write_field(unit, fields(11), G%Domain%mpp_domain, out_h)
  do j=js,je; do I=Isq,Ieq; out_h(I,j) = G%dyT(I,j); enddo; enddo
  call write_field(unit, fields(12), G%Domain%mpp_domain, out_h)

  do J=Jsq,Jeq; do i=is,ie; out_q(i,J) = G%dxBu(i,J); enddo; enddo
  call write_field(unit, fields(13), G%Domain%mpp_domain, out_q)
  do j=js,je; do I=Isq,Ieq; out_q(I,j) = G%dyBu(I,j); enddo; enddo
  call write_field(unit, fields(14), G%Domain%mpp_domain, out_q)

  do j=js,je; do i=is,ie; out_h(i,j) = G%areaT(i,j); enddo; enddo
  call write_field(unit, fields(15), G%Domain%mpp_domain, out_h)
  do j=js,je; do i=is,ie; out_q(i,j) = G%areaBu(i,j); enddo; enddo
  call write_field(unit, fields(16), G%Domain%mpp_domain, out_q)

!  do J=Jsq,Jeq; do i=is,ie; out_v(i,J) = G%dx_Cv(i,J); enddo; enddo
  call write_field(unit, fields(17), G%Domain%mpp_domain, G%dx_Cv)
!  do j=js,je; do I=Isq,Ieq; out_u(I,j) = G%dy_Cu(I,j); enddo; enddo
  call write_field(unit, fields(18), G%Domain%mpp_domain, G%dy_Cu)
  call write_field(unit, fields(19), G%Domain%mpp_domain, G%mask2dT)

  if (G%bathymetry_at_vel) then
    call write_field(unit, fields(20), G%Domain%mpp_domain, G%Dblock_u)
    call write_field(unit, fields(21), G%Domain%mpp_domain, G%Dopen_u)
    call write_field(unit, fields(22), G%Domain%mpp_domain, G%Dblock_v)
    call write_field(unit, fields(23), G%Domain%mpp_domain, G%Dopen_v)
  endif

  call close_file(unit)

end subroutine write_ocean_geometry_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine write_vertgrid_file(G, param_file, directory)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
  character(len=*),      intent(in)    :: directory
!   This subroutine writes out a file containing any available data related
! to the vertical grid used by the MOM ocean model.
! Arguments: G - The ocean's grid structure.  Effectively intent in.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      directory - The directory into which to place the file.
  character(len=120) :: filepath
  type(vardesc) :: vars(2)
  type(fieldtype) :: fields(2)
  integer :: unit

  filepath = trim(directory) // trim("Vertical_coordinate")

  vars(1) = vardesc("R","Target Potential Density",'1','L','1',"kilogram meter-3")
  vars(2) = vardesc("g","Reduced gravity",'1','L','1',"meter second-2")

  call create_file(unit, trim(filepath), vars, 2, G, fields, SINGLE_FILE)

  call write_field(unit, fields(1), G%Rlay)
  call write_field(unit, fields(2), G%g_prime)

  call close_file(unit)

end subroutine write_vertgrid_file
! -----------------------------------------------------------------------------

subroutine MOM_temp_salt_initialize_from_Z(h, tv, G, PF, dirs)

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
!  (in)      PF      - A structure indicating the open file to parse for
!                      model parameter values.
!  (in)      dirs    - A structure containing several relevant directory paths.


  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(out)   :: h    
  type(thermo_var_ptrs),                 intent(inout) :: tv
  type(ocean_grid_type),                 intent(inout)    :: G
  type(param_file_type),                 intent(in)    :: PF
  type(directories),                     intent(in)    :: dirs

  character(len=200) :: filename   ! The name of an input file containing temperature
                                   ! and salinity in z-space.
  character(len=200) :: inputdir ! The directory where NetCDF input files are.
  character(len=200) :: mesg

  type(EOS_type), pointer :: eos => NULL()

! This include declares and sets the variable "version".
#include "version_variable.h"


  character(len=40)  :: mod = "MOM_initialize_layers_from_Z" ! This module's name.



  integer :: is, ie, js, je, nz ! compute domain indices
  integer :: isc,iec,jsc,jec    ! global compute domain indices
  integer :: isg, ieg, jsg, jeg ! global extent
  integer :: isd, ied, jsd, jed ! data domain indices

  integer :: rcode, no_fill
  integer :: ndims                 
  integer :: i, j, k, ks, np, ni, nj
  integer :: idbg, jdbg
  integer :: nkml, nkbl         ! number of mixed and buffer layers
  integer :: i_offset, j_offset ! Offsets between the global grid and the local
                                ! 1-indexed version of the global grid.
  integer :: ncid, varid_t, varid_s
  integer :: id, jd, kd, jdp, inconsistent
  integer :: im,jm
  real    :: PI_180             ! for conversion from degrees to radians
  real    :: npole, pole
  real    :: max_lat, min_depth, max_depth
  real    :: dilate
  real    :: missing_value_temp, missing_value_salt    
  logical :: new_sim
  logical :: correct_thickness
  type(horiz_interp_type) :: Interp
  character(len=40) :: potemp_var, salin_var
  character(len=8)  :: laynum

  integer, parameter :: niter=10   ! number of iterations for t/s adjustment to layer density
  logical            :: adjust_temperature = .true.  ! fit t/s to target densities
  real, parameter    :: missing_value = -1.e20
  real, parameter    :: temp_land_fill = 0.0, salt_land_fill = 35.0

  !data arrays
  real, dimension(:,:), allocatable :: x_in, y_in
  real, dimension(:), allocatable :: lon_in, lon_in_p, lat_in, lat_in_p, last_row
  real, dimension(:), allocatable :: z_in, z_edges_in, Rb
  real, dimension(:,:), allocatable :: temp_in, salt_in, mask_in


  integer, dimension(4) :: start, count, dims    
  
  real, dimension(:,:), allocatable :: tmp_in ! A 2-d array for holding input data.  



  real, dimension(:,:,:), allocatable :: temp_z, salt_z, mask_z, rho_z


  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: zi

  real, dimension(:,:), allocatable :: Depth

  real, dimension(SZI_(G),SZJ_(G))  :: temp2,salt2,good2,fill2,temp_prev2, salt_prev2
  real, dimension(SZI_(G),SZJ_(G))  :: temp_out, salt_out, rho_out, mask_out
  real, dimension(SZI_(G),SZJ_(G)) :: lon_out, lat_out, good, fill
  real, dimension(SZI_(G),SZJ_(G))  :: nlevs
  real, dimension(SZI_(G))   :: press


  logical :: reentrant_x, tripolar_n, add_np,dbg
  logical :: debug = .false.  ! manually set this to true for verbose output

  ! Local variables for ALE remapping
  real, dimension(:), allocatable :: h1, h2, hTarget, deltaE, tmpT1d, tmpS1d
  real, dimension(:), allocatable :: tmpT1dIn, tmpS1dIn
  real :: zTopOfCell, zBottomOfCell
  type(regridding_CS) :: regridCS ! Regridding parameters and work arrays
  type(remapping_CS) :: remapCS ! Remapping parameters and work arrays

  real, dimension(:,:,:), allocatable :: tmp1
  logical :: homogenize, useALEremapping
  character(len=10) :: remappingScheme
  real :: tempAvg, saltAvg
  integer :: nPoints, ans
  integer :: id_clock_routine, id_clock_read, id_clock_interp, id_clock_fill, id_clock_ALE

  id_clock_routine = cpu_clock_id('(Initialize from Z)', grain=CLOCK_ROUTINE)
  id_clock_read = cpu_clock_id('(Initialize from Z) read', grain=CLOCK_LOOP)
  id_clock_interp = cpu_clock_id('(Initialize from Z) interp', grain=CLOCK_LOOP)
  id_clock_fill = cpu_clock_id('(Initialize from Z) fill', grain=CLOCK_LOOP)
  id_clock_ALE = cpu_clock_id('(Initialize from Z) ALE', grain=CLOCK_LOOP)

  call cpu_clock_begin(id_clock_routine)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  PI_180=atan(1.0)/45.

  call callTree_enter(trim(mod)//"(), MOM_initialization.F90")
  call log_version(PF, mod, version)

  new_sim = .false.
  if ((dirs%input_filename(1:1) == 'n') .and. &
       (LEN_TRIM(dirs%input_filename) == 1)) new_sim = .true.

  inputdir = "." ;  call get_param(PF, mod, "INPUTDIR", inputdir)
  inputdir = slasher(inputdir)    

  eos => tv%eqn_of_state

  call mpp_get_compute_domain(G%domain%mpp_domain,isc,iec,jsc,jec)

  reentrant_x = .false. ;  call get_param(PF, mod, "REENTRANT_X", reentrant_x,default=.true.)
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

  call cpu_clock_begin(id_clock_read)

  rcode = NF90_OPEN(filename, NF90_NOWRITE, ncid)
  if (rcode .ne. 0) call MOM_error(FATAL,"error opening file "//trim(filename)//&
                           " in MOM_initialize_layers_from_Z")
  rcode = NF90_INQ_VARID(ncid, potemp_var, varid_t)
  if (rcode .ne. 0) call MOM_error(FATAL,"error finding variable "//trim(potemp_var)//&
                                 " in file "//trim(filename)//" in MOM_initialize_layers_from_Z")        
  rcode = NF90_INQ_VARID(ncid, salin_var, varid_s)
  if (rcode .ne. 0) call MOM_error(FATAL,"error finding variable "//trim(salin_var)//&
                                 " in file "//trim(filename)//" in MOM_initialize_layers_from_Z")        
  rcode = NF90_INQUIRE_VARIABLE(ncid, varid_t, ndims=ndims, dimids=dims)
  if (rcode .ne. 0) call MOM_error(FATAL,'error inquiring dimensions MOM_initialize_layers_from_Z')            
  if (ndims < 3) call MOM_error(FATAL,"Variable "//trim(potemp_var)//" in file "// &
              trim(filename)//" has too few dimensions.")            
  rcode = NF90_INQUIRE_DIMENSION(ncid, dims(1), len=id)
  if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 1 data for "// &
                trim(potemp_var)//" in file "// trim(filename)//" in MOM_initialize_layers_from_Z")               
  rcode = NF90_INQUIRE_DIMENSION(ncid, dims(2), len=jd)
  if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 2 data for "// &
                trim(potemp_var)//" in file "// trim(filename)//" in MOM_initialize_layers_from_Z")               
  rcode = NF90_INQUIRE_DIMENSION(ncid, dims(3), len=kd)
  if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 3 data for "// &
                trim(potemp_var)//" in file "// trim(filename)//" in MOM_initialize_layers_from_Z")               

  missing_value_temp=0.0 ; missing_value_salt=0.0
  rcode = NF90_GET_ATT(ncid, varid_t, "_FillValue", missing_value_temp)
  if (rcode .ne. 0) call MOM_error(FATAL,"error finding missing value for "//&
       trim(potemp_var)//" in file "// trim(filename)//" in MOM_initialize_layers_from_Z")    
  rcode = NF90_GET_ATT(ncid, varid_s, "_FillValue", missing_value_salt)
  if (rcode .ne. 0) call MOM_error(FATAL,"error finding missing value for "//&
       trim(salin_var)//" in file "// trim(filename)//" in MOM_initialize_layers_from_Z")    

  allocate(lon_in(id),lat_in(jd),z_in(kd),z_edges_in(kd+1))
  allocate(temp_z(isd:ied,jsd:jed,kd), salt_z(isd:ied,jsd:jed,kd), rho_z(isd:ied,jsd:jed,kd),&
       & mask_z(isd:ied,jsd:jed,kd))

  start = 1; count = 1; count(1) = id
  rcode = NF90_GET_VAR(ncid, dims(1), lon_in, start, count)
  if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 1 values for "// &
                trim(potemp_var)//" in file "// trim(filename)//" in MOM_initialize_layers_from_Z")               
  start = 1; count = 1; count(1) = jd
  rcode = NF90_GET_VAR(ncid, dims(2), lat_in, start, count)
  if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 2 values for "// &
                trim(potemp_var)//" in file "// trim(filename)//" in MOM_initialize_layers_from_Z")               
  start = 1; count = 1; count(1) = kd
  rcode = NF90_GET_VAR(ncid, dims(3), z_in, start, count)
  if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 3 values for "// &
                trim(potemp_var)//" in file "// trim(filename)//" in MOM_initialize_layers_from_Z")               

! extrapolate the input data to the north pole using the northerm-most latitude

  max_lat = maxval(lat_in)
  add_np=.false.
  if (max_lat < 90.0) then
    add_np=.true.
    jdp=jd+1
    allocate(lat_in_p(jdp))
    lat_in_p(1:jd)=lat_in(:)
    lat_in_p(jd+1)=90.0
    deallocate(lat_in)
    allocate(lat_in(1:jdp))
    lat_in(:)=lat_in_p(:)
  else
    jdp=jd
  endif

! construct level cell boundaries as the mid-point between adjacent centers

  z_edges_in(1) = 0.0
  do k=2,kd
   z_edges_in(k)=0.5*(z_in(k-1)+z_in(k))
  enddo
  z_edges_in(kd+1)=2.0*z_in(kd) - z_in(kd-1)

  call horiz_interp_init()

  lon_in = lon_in*PI_180
  lat_in = lat_in*PI_180
  allocate(x_in(id,jdp),y_in(id,jdp))        
  call meshgrid(lon_in,lat_in, x_in, y_in)

  temp_prev2 = 0.0; salt_prev2 = 0.0

  lon_out(:,:) = G%geoLonT(:,:)*PI_180
  lat_out(:,:) = G%geoLatT(:,:)*PI_180


  allocate(tmp_in(id,jd)) ; tmp_in(:,:)=0.0
  allocate(temp_in(id,jdp)) ; temp_in(:,:)=0.0
  allocate(salt_in(id,jdp)) ; salt_in(:,:)=0.0
  allocate(mask_in(id,jdp)) ; mask_in(:,:)=0.0
  allocate(last_row(id))    ; last_row(:)=0.0

  ni=ieg-isg+1 ; nj = jeg-jsg+1
  allocate(Depth(ni,nj))

  press(:)=tv%p_ref

! get the global depth array

  call mpp_global_field(G%domain%mpp_domain, G%bathyT, Depth)    

  max_depth = maxval(Depth)
  if (z_edges_in(kd+1)<max_depth) z_edges_in(kd+1)=max_depth


! loop through each data level and interpolate to model grid.
! after interpolating, fill in points which will be needed
! to define the layers

  call cpu_clock_end(id_clock_read)
  do k=1,kd
    call cpu_clock_begin(id_clock_read)
    write(laynum,'(I8)') k ; laynum = adjustl(laynum)

    if (is_root_pe()) then
      start = 1; start(3) = k; count = 1; count(1) = id; count(2) = jd
      rcode = NF90_GET_VAR(ncid,varid_t, tmp_in, start, count)
      if (rcode .ne. 0) call MOM_error(FATAL,"MOM_initialize_layers_from_Z: "//&
           "error reading level "//trim(laynum)//" of variable "//&
           trim(potemp_var)//" in file "// trim(filename))

      if (add_np) then
         last_row(:)=tmp_in(:,jd); pole=0.0;npole=0.0
         do i=1,id
            if (abs(tmp_in(i,jd)-missing_value_temp) .gt. abs(G%Angstrom_Z*missing_value_temp)) then
               pole = pole+last_row(i)
               npole = npole+1.0
            endif
         enddo
         if (npole > 0) then
            pole=pole/npole
         else
            pole=missing_value
         endif
         temp_in(:,1:jd) = tmp_in(:,:)
         temp_in(:,jdp) = pole
      else
         temp_in(:,:) = tmp_in(:,:)
      endif
      
      rcode = NF90_GET_VAR(ncid, varid_s, tmp_in, start, count)
      if (rcode .ne. 0) call MOM_error(FATAL,"MOM_initialize_layers_from_Z: "//&
           "error reading level "//trim(laynum)//" of variable "//&
           trim(salin_var)//" in file "// trim(filename))
      
      if (add_np) then
         last_row(:)=tmp_in(:,jd); pole=0.0;npole=0.0
         do i=1,id
            if (abs(tmp_in(i,jd)-missing_value_salt) .gt. abs(G%Angstrom_Z*missing_value_salt)) then              
               pole = pole+last_row(i)
               npole = npole+1.0
            endif
         enddo
         if (npole > 0) then
            pole = pole/npole
         else
            pole = missing_value
         endif
         salt_in(:,1:jd) = tmp_in(:,:)
         salt_in(:,jdp) = pole
      else
         salt_in(:,:) = tmp_in(:,:)
      endif
    endif
  
    call mpp_sync()
    call mpp_broadcast(temp_in,id*jdp,root_PE())
    call mpp_broadcast(salt_in,id*jdp,root_PE())
    call mpp_sync_self ()

    mask_in=0.0

    do j=1,jdp
      do i=1,id
         if (abs(temp_in(i,j)-missing_value_temp) .gt. abs(G%Angstrom_Z*missing_value_temp)) then                           
           mask_in(i,j)=1.0
         else
           temp_in(i,j)=missing_value
           salt_in(i,j)=missing_value
         endif
      enddo 
    enddo


    call cpu_clock_end(id_clock_read)
    call cpu_clock_begin(id_clock_interp)
    
! call fms routine horiz_interp to interpolate input level data to model horizontal grid

    if (k == 1) then
      call horiz_interp_new(Interp,x_in,y_in,lon_out(is:ie,js:je),lat_out(is:ie,js:je), &
               interp_method='bilinear',src_modulo=reentrant_x)
    endif

    if (debug) then
       call myStats(temp_in,missing_value, k,'Temp from file')
       call myStats(salt_in,missing_value, k,'Salt from file')
    endif

    temp_out(:,:) = 0.0
    salt_out(:,:) = 0.0

    call horiz_interp(Interp,temp_in,temp_out(is:ie,js:je), missing_value=missing_value, new_missing_handle=.true.)

    mask_out=1.0
    do j=js,je
      do i=is,ie
        if (abs(temp_out(i,j)-missing_value) .lt. abs(G%Angstrom_Z*missing_value)) mask_out(i,j)=0.
      enddo
    enddo

    call horiz_interp(Interp,salt_in,salt_out(is:ie,js:je), missing_value=missing_value, new_missing_handle=.true.)

    if (debug) then
      call hchksum(temp_out,'temperature after hinterp ',G)
      call hchksum(salt_out,'salinity after hinterp ',G)
    endif


    call cpu_clock_end(id_clock_interp)
    call cpu_clock_begin(id_clock_fill)
    fill = 0.0; good = 0.0

    nPoints = 0 ; tempAvg = 0. ; saltAvg = 0.
    do j=js,je
      do i=is,ie
        if (mask_out(i,j) .lt. 1.0) then
          temp_out(i,j)=missing_value
          salt_out(i,j)=missing_value
        else
          good(i,j)=1.0
          nPoints = nPoints + 1
          tempAvg = tempAvg + temp_out(i,j)
          saltAvg = saltAvg + salt_out(i,j)
        endif
        if (G%mask2dT(i,j) == 1.0 .and. z_edges_in(k) <= G%bathyT(i,j) .and. mask_out(i,j) .lt. 1.0) fill(i,j)=1.0
      enddo
    enddo
    call pass_var(fill,G%Domain)
    call pass_var(good,G%Domain)

    if (debug) then
      call myStats(temp_out,missing_value, k,'Temp from horiz_interp()')
      call myStats(salt_out,missing_value, k,'Salt from horiz_interp()')
    endif

    ! Horizontally homogenize data to produce perfectly "flat" initial conditions
    if (homogenize) then
      call sum_across_PEs(nPoints)
      call sum_across_PEs(tempAvg)
      call sum_across_PEs(saltAvg)
      if (nPoints>0) then
        tempAvg = tempAvg/real(nPoints)
        saltAvg = saltAvg/real(nPoints)
      endif
      temp_out(:,:) = tempAvg
      salt_out(:,:) = saltAvg
    endif

! temp_out,salt_out contain input z-space data on the model grid with missing values
! now fill in missing values using "ICE-nine" algorithm. 

    temp2(:,:)=temp_out(:,:);good2(:,:)=good(:,:); fill2(:,:)=fill(:,:)
    salt2(:,:)=salt_out(:,:) 

    call fill_miss_2d(temp2,good2,fill2,temp_prev2,G,smooth=.true.)
    call myStats(temp2,missing_value,k,'Temp from fill_miss_2d()')
    call fill_miss_2d(salt2,good2,fill2,salt_prev2,G,smooth=.true.)
    call myStats(salt2,missing_value,k,'Salt from fill_miss_2d()')

    temp_z(:,:,k) = temp2(:,:)*G%mask2dT(:,:)
    salt_z(:,:,k) = salt2(:,:)*G%mask2dT(:,:)
    mask_z(:,:,k) = good2(:,:)+fill2(:,:)

    temp_prev2(:,:)=temp_z(:,:,k)
    salt_prev2(:,:)=salt_z(:,:,k)
    call cpu_clock_end(id_clock_fill)

    if (debug) then
      call hchksum(temp2,'temperature after fill ',G)
      call hchksum(salt2,'salinity after fill ',G)
    endif

! next use the equation of state to create a potential density field using filled z-data

    do j=js,je
      call calculate_density(temp_z(:,j,k),salt_z(:,j,k), press, rho_z(:,j,k), is, ie, eos)
    enddo

  enddo ! kd

  call pass_var(temp_z,G%Domain)
  call pass_var(salt_z,G%Domain)
  call pass_var(mask_z,G%Domain)
  call pass_var(rho_z,G%Domain)

  call horiz_interp_del(Interp)    

! Done with horizontal interpolation.    
! Now remap to model coordinates
  if (useALEremapping) then
    call cpu_clock_begin(id_clock_ALE)
    ! First we reserve a work space for reconstructions of the source data
    allocate( h1(kd) )
    allocate( tmpT1dIn(kd) )
    allocate( tmpS1dIn(kd) )
    call initialize_remapping( kd, remappingScheme, remapCS ) ! Data for reconstructions
    call remapDisableBoundaryExtrapolation( remapCS )
    ! Next we initialize the regridding package so that it knows about the target grid
    allocate( hTarget(nz) )
    allocate( h2(nz) )
    allocate( tmpT1d(nz) )
    allocate( tmpS1d(nz) )
    allocate( deltaE(nz+1) )
    ! This call can be more general but is hard-coded for z* coordinates...  ????
    call ALE_initRegridding( G, PF, mod, regridCS, hTarget ) ! sets regridCS and hTarget(1:nz)
    ! For each column ...
    do j = js, je ; do i = is, ie
      if (G%mask2dT(i,j)>0.) then
        ! Build the source grid
        zTopOfCell = 0. ; zBottomOfCell = 0. ; nPoints = 0
        do k = 1, kd
          if (mask_z(i,j,k) > 0.) then
            zBottomOfCell = -min( z_edges_in(k+1), G%bathyT(i,j) )
            tmpT1dIn(k) = temp_z(i,j,k)
            tmpS1dIn(k) = salt_z(i,j,k)
          elseif (k>1) then
            zBottomOfCell = -G%bathyT(i,j) 
            tmpT1dIn(k) = tmpT1dIn(k-1)
            tmpS1dIn(k) = tmpS1dIn(k-1)
          else ! This next block should only ever be reached over land
            tmpT1dIn(k) = -99.9
            tmpS1dIn(k) = -99.9
          endif
          h1(k) = zTopOfCell - zBottomOfCell
          if (h1(k)>0.) nPoints = nPoints + 1
          zTopOfCell = zBottomOfCell ! Bottom becomes top for next value of k
        enddo
        h1(kd) = h1(kd) + ( zTopOfCell + G%bathyT(i,j) ) ! In case data is deeper than model
        ! Build the target grid combining hTarget and topography
        zTopOfCell = 0. ; zBottomOfCell = 0.
        do k = 1, nz
          zBottomOfCell = max( zTopOfCell - hTarget(k), -G%bathyT(i,j) )
          h2(k) = zTopOfCell - zBottomOfCell
          zTopOfCell = zBottomOfCell ! Bottom becomes top for next value of k
        enddo
        ! Calcaulate an effectiveadisplacement, deltaE
        call dzFromH1H2( nPoints, h1, nz, h2, deltaE ) ! sets deltaE
        ! Now remap from h1 to h2=h1+div.deltaE
        call remapping_core( remapCS, nPoints, h1, tmpT1dIn, nz, deltaE, tmpT1d ) ! sets tmpT1d
        call remapping_core( remapCS, nPoints, h1, tmpS1dIn, nz, deltaE, tmpS1d ) ! sets tmpS1d
        h(i,j,:) = h2(:)
        tv%T(i,j,:) = tmpT1d(:)
        tv%S(i,j,:) = tmpS1d(:)
      else
        tv%T(i,j,:) = 0.
        tv%S(i,j,:) = 0.
        h(i,j,:) = 0.
      endif ! mask2dT
    enddo ; enddo
    deallocate( h1 )
    deallocate( h2 )
    deallocate( hTarget )
    deallocate( tmpT1d )
    deallocate( tmpS1d )
    deallocate( tmpT1dIn )
    deallocate( tmpS1dIn )
    deallocate( deltaE )

    do k=1,nz
      call myStats(tv%T(is:ie,js:je,k),missing_value,k,'Temp from ALE()')
    enddo
    call cpu_clock_end(id_clock_ALE)
  else ! remap to isopycnal layer space
! next find interface positions using local arrays
! nlevs contains the number of valid data points in each column

    nlevs = sum(mask_z,dim=3)


! Rb contains the layer interface densities
    allocate(Rb(nz+1))
    do k=2,nz
       Rb(k)=0.5*(G%Rlay(k-1)+G%Rlay(k))
    enddo
    Rb(1)=0.0
    Rb(nz+1)=2.0*G%Rlay(nz) - G%Rlay(nz-1)

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
      call adjustEtaToFitBathymetry(G, zi, h)
    else
      do k=nz,1,-1 ; do j=js,je ; do i=is,ie
        if (zi(i,j,K) < (zi(i,j,K+1) + G%Angstrom_z)) then
          zi(i,j,K) = zi(i,j,K+1) + G%Angstrom_z
          h(i,j,k) = G%Angstrom_z
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



    tv%T(is:ie,js:je,:) = tracer_z_init(temp_z(is:ie,js:je,:),-1.0*z_edges_in,zi(is:ie,js:je,:),nkml,nkbl,missing_value,G%mask2dT(is:ie,js:je),nz,nlevs(is:ie,js:je),dbg,idbg,jdbg)
    tv%S(is:ie,js:je,:) = tracer_z_init(salt_z(is:ie,js:je,:),-1.0*z_edges_in,zi(is:ie,js:je,:),nkml,nkbl,missing_value,G%mask2dT(is:ie,js:je),nz,nlevs(is:ie,js:je))

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
            G%Rlay(1:nz), tv%p_ref, niter, missing_value, h(is:ie,js:je,:), ks, eos)

  endif

  deallocate(lon_in,lat_in,x_in,y_in)
  deallocate(z_in,z_edges_in)
  deallocate(tmp_in,temp_in,salt_in,mask_in,last_row)
  deallocate(Depth)

  call callTree_leave(trim(mod)//'()')
  call cpu_clock_end(id_clock_routine)

  contains
  subroutine myStats(array, missing, k, mesg)
  real, dimension(:,:), intent(in) :: array
  real, intent(in) :: missing
  integer :: k
  character(len=*) :: mesg
  ! Local variables
  real :: minA, maxA
  integer :: i,j,is,ie,js,je
  logical :: found
  character(len=120) :: lMesg
  minA = 9.E24 ; maxA = -9.E24 ; found = .false.

  is = G%isc;ie=G%iec;js=G%jsc;je=G%jec

  do j = js, je
    do i = is, ie
      if (array(i,j) /= array(i,j)) stop 'Nan!'
      if (abs(array(i,j)-missing)>1.e-6*abs(missing)) then
        if (found) then
          minA = min(minA, array(i,j))
          maxA = max(maxA, array(i,j))
        else
          found = .true.
          minA = array(i,j)
          maxA = array(i,j)
        endif
      endif
    enddo
  enddo
  call min_across_PEs(minA)
  call max_across_PEs(maxA)
  if (is_root_pe()) then
    write(lMesg(1:120),'(2(a,es12.4),a,i3,x,a)') &
       'init_from_Z: min=',minA,' max=',maxA,' Level=',k,trim(mesg)
    call MOM_mesg(lMesg,8)
  endif
  end subroutine myStats
  
  subroutine fill_miss_2d(aout,good,fill,prev,G,smooth,num_pass,relc,crit,keep_bug,debug)
!
!# Use ICE-9 algorithm to populate points (fill=1) with 
!# valid data (good=1). If no information is available,
!# Then use a previous guess (prev). Optionally (smooth) 
!# blend the filled points to achieve a more desirable result.
!
!  (in)        a   : input 2-d array with missing values 
!  (in)     good   : valid data mask for incoming array (1==good data; 0==missing data)
!  (in)     fill   : same shape array of points which need filling (1==please fill;0==leave it alone)   
!  (in)     prev   : first guess where isolated holes exist,
!

    use MOM_coms, only : sum_across_PEs

    real, dimension(NIMEM_,NJMEM_), intent(inout) :: aout
    real, dimension(NIMEM_,NJMEM_), intent(in) :: good,fill
    real, dimension(NIMEM_,NJMEM_), optional, intent(in) :: prev
    type(ocean_grid_type), intent(inout)  :: G
    logical, intent(in), optional :: smooth
    integer, intent(in), optional :: num_pass
    real, intent(in), optional    :: relc,crit
    logical, intent(in), optional :: keep_bug, debug
    
    
    real, dimension(SZI_(G),SZJ_(G)) :: b,r
    real, dimension(SZI_(G),SZJ_(G)) :: fill_pts,good_,good_new   
    
    integer :: i,j,k
    real    :: east,west,north,south,sor
    real    :: ge,gw,gn,gs,ngood
    logical :: do_smooth,siena_bug
    real    :: nfill, nfill_prev
    integer, parameter :: num_pass_default = 10000
    real, parameter :: relc_default = 0.25, crit_default = 1.e-3
    
    integer :: npass
    integer :: is, ie, js, je, nz
    real    :: relax_coeff, acrit, ares
    logical :: debug_it

    debug_it=.false.
    if (PRESENT(debug)) debug_it=debug

    is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
   
    npass = num_pass_default
    if (PRESENT(num_pass)) npass = num_pass

    relax_coeff = relc_default
    if (PRESENT(relc)) relax_coeff = relc

    acrit = crit_default
    if (PRESENT(crit)) acrit = crit

    siena_bug=.false.
    if (PRESENT(keep_bug)) siena_bug = keep_bug

    do_smooth=.false.
    if (PRESENT(smooth)) do_smooth=smooth
   
    fill_pts(:,:)=fill(:,:)

    nfill = sum(fill(is:ie,js:je))
    call sum_across_PEs(nfill)

    nfill_prev = nfill
    good_(:,:)=good(:,:)
    r(:,:)=0.0

    do while (nfill > 0.0)

       call pass_var(good_,G%Domain)
       call pass_var(aout,G%Domain)

       b(:,:)=aout(:,:)
       good_new(:,:)=good_(:,:)

       do j=js,je
          i_loop: do i=is,ie

             if (good_(i,j) .eq. 1.0 .or. fill(i,j) .eq. 0.) cycle i_loop

             ge=good_(i+1,j);gw=good_(i-1,j)
             gn=good_(i,j+1);gs=good_(i,j-1)
             east=0.0;west=0.0;north=0.0;south=0.0
             if (ge.eq.1.0) east=aout(i+1,j)*ge
             if (gw.eq.1.0) west=aout(i-1,j)*gw
             if (gn.eq.1.0) north=aout(i,j+1)*gn
             if (gs.eq.1.0) south=aout(i,j-1)*gs     

             ngood = ge+gw+gn+gs
             if (ngood > 0.) then
                b(i,j)=(east+west+north+south)/ngood
                fill_pts(i,j)=0.0
                good_new(i,j)=1.0
             endif
          enddo i_loop
       enddo

       aout(is:ie,js:je)=b(is:ie,js:je)
       good_(is:ie,js:je)=good_new(is:ie,js:je)
       nfill_prev = nfill
       nfill = sum(fill_pts(is:ie,js:je))
       call sum_across_PEs(nfill)

       if (nfill == nfill_prev .and. PRESENT(prev)) then
          do j=js,je
             do i=is,ie
                if (fill_pts(i,j).eq.1.0) then
                   aout(i,j)=prev(i,j)
                   fill_pts(i,j)=0.0
                endif
             enddo
          enddo
       else if (nfill .eq. nfill_prev) then
          print *,&
               'Unable to fill missing points using either data at the same vertical level from a connected basin'//&
               'or using a point from a previous vertical level.  Make sure that the original data has some valid'//&
               'data in all basins.'
          print *,'nfill=',nfill
       endif

       nfill = sum(fill_pts(is:ie,js:je))
       call sum_across_PEs(nfill)
       
    end do

    if (do_smooth) then
       do k=1,npass
          call pass_var(aout,G%Domain)
          do j=js,je
             do i=is,ie
                if (fill(i,j) .eq. 1) then
                east=max(good(i+1,j),fill(i+1,j));west=max(good(i-1,j),fill(i-1,j))
                north=max(good(i,j+1),fill(i,j+1));south=max(good(i,j-1),fill(i,j-1))
                r(i,j) = relax_coeff*(south*aout(i,j-1)+north*aout(i,j+1)+west*aout(i-1,j)+east*aout(i+1,j) - (south+north+west+east)*aout(i,j))
               else
                r(i,j) = 0.
               endif
             enddo
          enddo
          aout(is:ie,js:je)=r(is:ie,js:je)+aout(is:ie,js:je)
          ares = maxval(abs(r))
          call max_across_PEs(ares)
          if (ares <= acrit) exit
       enddo
    endif

    do j=js,je
       do i=is,ie
          if (good_(i,j).eq.0.0 .and. fill_pts(i,j) .eq. 1.0) then
             print *,'in fill_miss, fill, good,i,j= ',fill_pts(i,j),good_(i,j),i,j
             call MOM_error(FATAL,"MOM_initialize: "// &
                  "fill is true and good is false after fill_miss, how did this happen? ")
          endif
       enddo
    enddo
    
    return
    
  end subroutine fill_miss_2d

end subroutine MOM_temp_salt_initialize_from_Z

end module MOM_initialization
