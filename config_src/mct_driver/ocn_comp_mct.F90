!> This is the main driver for MOM6 in CIME
module ocn_comp_mct

! This file is part of MOM6. See LICENSE.md for the license.

! mct modules
use ESMF,                only: ESMF_clock, ESMF_time, ESMF_timeInterval, ESMF_TimeInc
use ESMF,                only: ESMF_ClockGet, ESMF_TimeGet, ESMF_TimeIntervalGet
use seq_cdata_mod,       only: seq_cdata, seq_cdata_setptrs
use seq_flds_mod,        only: ice_ncat, seq_flds_i2o_per_cat
use mct_mod,             only: mct_gsMap, mct_gsmap_init, mct_gsMap_lsize, &
                               mct_gsmap_orderedpoints
use mct_mod,             only: mct_aVect, mct_aVect_init, mct_aVect_zero, &
                               mct_aVect_nRattr
use mct_mod,             only: mct_gGrid, mct_gGrid_init, mct_gGrid_importRAttr, &
                               mct_gGrid_importIAttr
use mct_mod,             only: mct_avect_indexra, mct_aVect_clean
use seq_flds_mod,        only: seq_flds_x2o_fields, seq_flds_o2x_fields, seq_flds_dom_coord, &
                               seq_flds_dom_other
use seq_infodata_mod,    only: seq_infodata_type, seq_infodata_GetData, &
                               seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                               seq_infodata_start_type_brnch, seq_infodata_PutData
use seq_comm_mct,        only: seq_comm_name, seq_comm_inst, seq_comm_suffix
use seq_timemgr_mod,     only: seq_timemgr_EClockGetData, seq_timemgr_RestartAlarmIsOn
use perf_mod,            only: t_startf, t_stopf
use shr_kind_mod,        only: shr_kind_r8
use shr_file_mod,        only: shr_file_getUnit, shr_file_freeUnit, shr_file_setIO, &
                               shr_file_getLogUnit, shr_file_getLogLevel, &
                               shr_file_setLogUnit, shr_file_setLogLevel

! MOM6 modules
use MOM_coms,             only : reproducing_sum
use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,        only : CLOCK_SUBCOMPONENT
use MOM,                  only: initialize_MOM, step_MOM, MOM_control_struct, MOM_end
use MOM,                  only: calculate_surface_state, allocate_surface_state
use MOM,                  only: finish_MOM_initialization, step_offline
use MOM_forcing_type,     only: forcing, forcing_diags, register_forcing_type_diags
use MOM_forcing_type,     only: allocate_forcing_type, deallocate_forcing_type
use MOM_forcing_type,     only: mech_forcing_diags, forcing_accumulate, forcing_diagnostics
use MOM_restart,          only: save_restart
use MOM_domains,          only: MOM_infra_init, num_pes, root_pe, pe_here
use MOM_domains,          only: pass_vector, BGRID_NE, CGRID_NE, To_All
use MOM_domains,          only: pass_var, AGRID, fill_symmetric_edges
use MOM_grid,             only: ocean_grid_type, get_global_grid_size
use MOM_verticalGrid,     only: verticalGrid_type
use MOM_variables,        only: surface
use MOM_error_handler,    only: MOM_error, FATAL, is_root_pe, WARNING
use MOM_error_handler,    only: callTree_enter, callTree_leave
use MOM_time_manager,     only: time_type, set_date, set_time, set_calendar_type, NOLEAP, get_date
use MOM_time_manager,     only: operator(+), operator(-), operator(*), operator(/)
use MOM_time_manager,     only: operator(/=), operator(>), get_time
use MOM_file_parser,      only: get_param, log_version, param_file_type, close_param_file
use MOM_get_input,        only: Get_MOM_Input, directories
use MOM_diag_mediator,    only: diag_ctrl, enable_averaging, disable_averaging
use MOM_diag_mediator,    only: diag_mediator_close_registration, diag_mediator_end
use MOM_diag_mediator,    only: safe_alloc_ptr
use MOM_ice_shelf,        only: initialize_ice_shelf, shelf_calc_flux, ice_shelf_CS
use MOM_ice_shelf,        only: ice_shelf_end, ice_shelf_save_restart
use MOM_sum_output,       only: MOM_sum_output_init, sum_output_CS
use MOM_sum_output,       only: write_energy, accumulate_net_input
use MOM_string_functions, only: uppercase
use MOM_constants,        only: CELSIUS_KELVIN_OFFSET, hlf, hlv
use MOM_EOS,              only: gsw_sp_from_sr, gsw_pt_from_ct
use user_revise_forcing,  only : user_alter_forcing, user_revise_forcing_init
use user_revise_forcing,  only : user_revise_forcing_CS
use MOM_restart,          only : restart_init, MOM_restart_CS
use MOM_restart,          only : restart_init_end, save_restart, restore_state
use data_override_mod,    only : data_override_init, data_override
use MOM_io,               only : slasher, write_version_number
use MOM_spatial_means,    only : adjust_area_mean_to_zero

! FMS modules
use mpp_domains_mod,          only : domain2d, mpp_get_layout, mpp_get_global_domain
use mpp_domains_mod,          only : mpp_define_domains, mpp_get_compute_domain
use time_interp_external_mod, only : init_external_field, time_interp_external
use time_interp_external_mod, only : time_interp_external_init
use fms_mod,                  only : read_data

! GFDL coupler modules
use coupler_types_mod,   only : coupler_1d_bc_type, coupler_2d_bc_type
use coupler_types_mod,   only : coupler_type_spawn
use coupler_types_mod,   only : coupler_type_initialized, coupler_type_copy_data

! By default make data private
implicit none; private

#include <MOM_memory.h>

! Public member functions
public :: ocn_init_mct
public :: ocn_run_mct
public :: ocn_final_mct
! Flag for debugging
logical, parameter :: debug=.true.

!> Structure with MCT attribute vectors and indices
type cpl_indices

  ! ocean to coupler
  integer :: o2x_So_t         !< Surface potential temperature (deg C)
  integer :: o2x_So_u         !< Surface zonal velocity (m/s)
  integer :: o2x_So_v         !< Surface meridional velocity (m/s)
  integer :: o2x_So_s         !< Surface salinity (PSU)
  integer :: o2x_So_dhdx      !< Zonal slope in the sea surface height
  integer :: o2x_So_dhdy      !< Meridional lope in the sea surface height
  integer :: o2x_So_bldepth   !< Boundary layer depth (m)
  integer :: o2x_Fioo_q       !< Heat flux?
  integer :: o2x_Faoo_fco2_ocn!< CO2 flux
  integer :: o2x_Faoo_fdms_ocn!< DMS flux

  ! coupler to ocean
  integer :: x2o_Si_ifrac        !< Fractional ice wrt ocean
  integer :: x2o_So_duu10n       !< 10m wind speed squared (m^2/s^2)
  integer :: x2o_Sa_pslv         !< Sea-level pressure (Pa)
  integer :: x2o_Sa_co2prog      !< Bottom atm level prognostic CO2
  integer :: x2o_Sa_co2diag      !< Bottom atm level diagnostic CO2
  integer :: x2o_Sw_lamult       !< Wave model langmuir multiplier
  integer :: x2o_Sw_ustokes      !< Surface Stokes drift, x-component
  integer :: x2o_Sw_vstokes      !< Surface Stokes drift, y-component
  integer :: x2o_Foxx_taux       !< Zonal wind stress (W/m2)
  integer :: x2o_Foxx_tauy       !< Meridonal wind stress (W/m2)
  integer :: x2o_Foxx_swnet      !< Net short-wave heat flux (W/m2)
  integer :: x2o_Foxx_sen        !< Sensible heat flux (W/m2)
  integer :: x2o_Foxx_lat        !< Latent heat flux  (W/m2)
  integer :: x2o_Foxx_lwup       !< Longwave radiation, up (W/m2)
  integer :: x2o_Faxa_lwdn       !< Longwave radiation, down (W/m2)
  integer :: x2o_Fioi_melth      !< Heat flux from snow & ice melt (W/m2)
  integer :: x2o_Fioi_meltw      !< Snow melt flux (kg/m2/s)
  integer :: x2o_Fioi_bcpho      !< Black Carbon hydrophobic release
                                 !! from sea ice component
  integer :: x2o_Fioi_bcphi      !< Black Carbon hydrophilic release from
                                 !! sea ice component
  integer :: x2o_Fioi_flxdst     !< Dust release from sea ice component
  integer :: x2o_Fioi_salt       !< Salt flux    (kg(salt)/m2/s)
  integer :: x2o_Foxx_evap       !< Evaporation flux  (kg/m2/s)
  integer :: x2o_Faxa_prec       !< Total precipitation flux (kg/m2/s)
  integer :: x2o_Faxa_snow       !< Water flux due to snow (kg/m2/s)
  integer :: x2o_Faxa_rain       !< Water flux due to rain (kg/m2/s)
  integer :: x2o_Faxa_bcphidry   !< Black   Carbon hydrophilic dry deposition
  integer :: x2o_Faxa_bcphodry   !< Black   Carbon hydrophobic dry deposition
  integer :: x2o_Faxa_bcphiwet   !< Black   Carbon hydrophilic wet deposition
  integer :: x2o_Faxa_ocphidry   !< Organic Carbon hydrophilic dry deposition
  integer :: x2o_Faxa_ocphodry   !< Organic Carbon hydrophobic dry deposition
  integer :: x2o_Faxa_ocphiwet   !< Organic Carbon hydrophilic dry deposition
  integer :: x2o_Faxa_dstwet1    !< Size 1 dust -- wet deposition
  integer :: x2o_Faxa_dstwet2    !< Size 2 dust -- wet deposition
  integer :: x2o_Faxa_dstwet3    !< Size 3 dust -- wet deposition
  integer :: x2o_Faxa_dstwet4    !< Size 4 dust -- wet deposition
  integer :: x2o_Faxa_dstdry1    !< Size 1 dust -- dry deposition
  integer :: x2o_Faxa_dstdry2    !< Size 2 dust -- dry deposition
  integer :: x2o_Faxa_dstdry3    !< Size 3 dust -- dry deposition
  integer :: x2o_Faxa_dstdry4    !< Size 4 dust -- dry deposition
  integer :: x2o_Foxx_rofl       !< River runoff flux (kg/m2/s)
  integer :: x2o_Foxx_rofi       !< Ice runoff flux (kg/m2/s)

  ! optional per thickness category fields
  integer, dimension(:), allocatable :: x2o_frac_col !< Fraction of ocean cell,
                                                     !! per column
  integer, dimension(:), allocatable :: x2o_fracr_col!< Fraction of ocean cell used
                                                     !! in radiation computations,
                                                     !! per column
  integer, dimension(:), allocatable :: x2o_qsw_fracr_col !< qsw * fracr, per column
  end type cpl_indices

!> This type is used for communication with other components via the FMS coupler.
! The element names and types can be changed only with great deliberation, hence
! the persistnce of things like the cutsy element name "avg_kount".
type, public ::  ocean_public_type
  type(domain2d) :: Domain    !< The domain for the surface fields.
  logical :: is_ocean_pe      !! .true. on processors that run the ocean model.
  character(len=32) :: instance_name = '' !< A name that can be used to identify
                              !! this instance of an ocean model, for example
                              !! in ensembles when writing messages.
  integer, pointer, dimension(:) :: pelist => NULL()   !< The list of ocean PEs.
  logical, pointer, dimension(:,:) :: maskmap =>NULL() !< A pointer to an array
                              !! indicating which logical processors are actually
                              !! used for the ocean code. The other logical
                              !! processors would be all land points and are not
                              !! assigned to actual processors. This need not be
                              !! assigned if all logical processors are used.

  integer :: stagger = -999 !< The staggering relative to the tracer points
                            !! of the two velocity components. Valid entries
                    !! include AGRID, BGRID_NE, CGRID_NE, BGRID_SW, and CGRID_SW,
                    !! corresponding to the community-standard Arakawa notation.
                    !! (These are named integers taken from mpp_parameter_mod.)
                    !! Following MOM, this is BGRID_NE by default when the ocean
                    !! is initialized, but here it is set to -999 so that a
                    !! global max across ocean and non-ocean processors  can be
                    !! used to determine its value.
  real, pointer, dimension(:,:)  :: &
    t_surf => NULL(), &  !< SST on t-cell (degrees Kelvin)
    s_surf => NULL(), &  !< SSS on t-cell (psu)
    u_surf => NULL(), &  !< i-velocity at the locations indicated by stagger, m/s.
    v_surf => NULL(), &  !< j-velocity at the locations indicated by stagger, m/s.
    sea_lev => NULL(), & !< Sea level in m after correction for surface pressure,
                         !! i.e. dzt(1) + eta_t + patm/rho0/grav (m)
    frazil =>NULL(), &   !< Accumulated heating (in Joules/m^2) from frazil
                         !! formation in the ocean.
    area => NULL()       !< cell area of the ocean surface, in m2.
  type(coupler_2d_bc_type) :: fields    !< A structure that may contain an
                                        !! array of named tracer-related fields.
  integer                  :: avg_kount !< Used for accumulating averages of this type.
  integer, dimension(2)    :: axes = 0  !< Axis numbers that are available
                                        ! for I/O using this surface data.
end type ocean_public_type

!> Contains pointers to the forcing fields which may be used to drive MOM.
!! All fluxes are positive downward.
type, public :: surface_forcing_CS ; private
  integer :: wind_stagger       !< AGRID, BGRID_NE, or CGRID_NE (integer values
                                !! from MOM_domains) to indicate the staggering of
                                !! the winds that are being provided in calls to
                                !! update_ocean_model. CIME uses AGRID, so this option
                                !! is being hard coded for now.
  logical :: use_temperature    !< If true, temp and saln used as state variables
  real :: wind_stress_multiplier!< A multiplier applied to incoming wind stress (nondim).
  ! smg: remove when have A=B code reconciled
  logical :: bulkmixedlayer     !< If true, model based on bulk mixed layer code
  real :: Rho0                  !< Boussinesq reference density (kg/m^3)
  real :: area_surf = -1.0      !< total ocean surface area (m^2)
  real :: latent_heat_fusion    !< latent heat of fusion (J/kg)
  real :: latent_heat_vapor     !< latent heat of vaporization (J/kg)
  real :: max_p_surf            !< maximum surface pressure that can be
                                !! exerted by the atmosphere and floating sea-ice,
                                !! in Pa.  This is needed because the FMS coupling
                                !! structure does not limit the water that can be
                                !! frozen out of the ocean and the ice-ocean heat
                                !! fluxes are treated explicitly.
  logical :: use_limited_P_SSH  !< If true, return the sea surface height with
                                !! the correction for the atmospheric (and sea-ice)
                                !! pressure limited by max_p_surf instead of the
                                !! full atmospheric pressure.  The default is true.
  real :: gust_const            !< constant unresolved background gustiness for ustar (Pa)
  logical :: read_gust_2d       !< If true, use a 2-dimensional gustiness supplied
                                !! from an input file.
  real, pointer, dimension(:,:) :: &
    TKE_tidal => NULL(), &      !< turbulent kinetic energy introduced to the
                                !! bottom boundary layer by drag on the tidal flows,
                                !! in W m-2.
    gust => NULL(), &           !< spatially varying unresolved background
                                !! gustiness that contributes to ustar (Pa).
                                !! gust is used when read_gust_2d is true.
    ustar_tidal => NULL()     !< tidal contribution to the bottom friction velocity (m/s)
  real :: cd_tides            !< drag coefficient that applies to the tides (nondimensional)
  real :: utide               !< constant tidal velocity to use if read_tideamp
                              !! is false, in m s-1.
  logical :: read_tideamp     !< If true, spatially varying tidal amplitude read from a file.
  logical :: rigid_sea_ice    !< If true, sea-ice exerts a rigidity that acts
                              !! to damp surface deflections (especially surface
                              !! gravity waves).  The default is false.
  real    :: Kv_sea_ice       !< viscosity in sea-ice that resists sheared vertical motions (m^2/s)
  real    :: density_sea_ice  !< typical density of sea-ice (kg/m^3). The value is
                              !! only used to convert the ice pressure into
                              !! appropriate units for use with Kv_sea_ice.
  real    :: rigid_sea_ice_mass !< A mass per unit area of sea-ice beyond which
                                !! sea-ice viscosity becomes effective, in kg m-2,
                                !! typically of order 1000 kg m-2.
  logical :: allow_flux_adjustments !< If true, use data_override to obtain flux adjustments
  real    :: Flux_const             !< piston velocity for surface restoring (m/s)
  logical :: salt_restore_as_sflux  !< If true, SSS restore as salt flux instead of water flux
  logical :: adjust_net_srestore_to_zero !< adjust srestore to zero (for both salt_flux or vprec)
  logical :: adjust_net_srestore_by_scaling !< adjust srestore w/o moving zero contour
  logical :: adjust_net_fresh_water_to_zero !< adjust net surface fresh-water (w/ restoring) to zero
  logical :: adjust_net_fresh_water_by_scaling !< adjust net surface fresh-water  w/o moving zero contour
  logical :: mask_srestore_under_ice        !< If true, use an ice mask defined by frazil
                                            !! criteria for salinity restoring.
  real    :: ice_salt_concentration         !< salt concentration for sea ice (kg/kg)
  logical :: mask_srestore_marginal_seas    !< if true, then mask SSS restoring in marginal seas
  real    :: max_delta_srestore             !< maximum delta salinity used for restoring
  real    :: max_delta_trestore             !< maximum delta sst used for restoring
  real, pointer, dimension(:,:) :: basin_mask => NULL() !< mask for SSS restoring
  type(diag_ctrl), pointer :: diag          !< structure to regulate diagnostic output timing
  character(len=200)       :: inputdir      !< directory where NetCDF input files are
  character(len=200)       :: salt_restore_file !< filename for salt restoring data
  character(len=30)        :: salt_restore_var_name !< name of surface salinity in salt_restore_file
  character(len=200)       :: temp_restore_file     !< filename for sst restoring data
  character(len=30)        :: temp_restore_var_name !< name of surface temperature in temp_restore_file
  integer :: id_srestore = -1     !< id number for time_interp_external.
  integer :: id_trestore = -1     !< id number for time_interp_external.
  type(forcing_diags), public :: handles !< diagnostics handles
  !###  type(ctrl_forcing_CS), pointer :: ctrl_forcing_CSp => NULL()
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()   !< restart pointer
  type(user_revise_forcing_CS), pointer :: urf_CS => NULL()!< user revise pointer
end type surface_forcing_CS

!> Contains information about the ocean state, although it is not necessary that
!! this is implemented with all models. This type is private, and can therefore vary
!! between different ocean models.
type, public :: ocean_state_type ; private
  logical :: is_ocean_PE = .false. !< True if this is an ocean PE.
  type(time_type) :: Time          !< The ocean model's time and master clock.
  integer :: Restart_control !< An integer that is bit-tested to determine whether
                             !! incremental restart files are saved and whether they
                             !! have a time stamped name.  +1 (bit 0) for generic
                             !! files and +2 (bit 1) for time-stamped files.  A
                             !! restart file is saved at the end of a run segment
                             !! unless Restart_control is negative.
  type(time_type) :: energysavedays  !< The interval between writing the energies
                                     !! and other integral quantities of the run.
  type(time_type) :: write_energy_time !< The next time to write to the energy file.
  integer :: nstep = 0        !< The number of calls to update_ocean.
  logical :: use_ice_shelf    !< If true, the ice shelf model is enabled.
  logical :: icebergs_apply_rigid_boundary !< If true, the icebergs can change ocean bd condition.
  real :: kv_iceberg          !< The viscosity of the icebergs in m2/s (for ice rigidity)
  real :: berg_area_threshold !< Fraction of grid cell which iceberg must occupy
                              !! so that fluxes below are set to zero. (0.5 is a
                              !! good value to use. Not applied for negative values.
  real :: latent_heat_fusion  !< Latent heat of fusion
  real :: density_iceberg     !< A typical density of icebergs in kg/m3 (for ice rigidity)
  type(ice_shelf_CS), pointer :: Ice_shelf_CSp => NULL() !< ice shelf structure.
  logical :: restore_salinity !< If true, the coupled MOM driver adds a term to
                              !! restore salinity to a specified value.
  logical :: restore_temp     !< If true, the coupled MOM driver adds a term to
                              !! restore sst to a specified value.
  real :: press_to_z          !< A conversion factor between pressure and ocean
                              !! depth in m, usually 1/(rho_0*g), in m Pa-1.
  real :: C_p                 !< The heat capacity of seawater, in J K-1 kg-1.
  type(directories) :: dirs   !< A structure containing several relevant directory paths.
  type(forcing)   :: fluxes   !< A structure containing pointers to
                              !! the ocean forcing fields.
  type(forcing)   :: flux_tmp !< A secondary structure containing pointers to the
                              !! ocean forcing fields for when multiple coupled
                              !! timesteps are taken per thermodynamic step.
  type(surface)   :: state    !< A structure containing pointers to
                              !! the ocean surface state fields.
  type(ocean_grid_type), pointer :: grid => NULL() !< A pointer to a grid structure
                              !! containing metrics and related information.
  type(verticalGrid_type), pointer :: GV => NULL() !< A pointer to a vertical grid
                              !! structure containing metrics and related information.
  type(MOM_control_struct), pointer :: MOM_CSp => NULL()
  type(surface_forcing_CS), pointer :: forcing_CSp => NULL()
  type(sum_output_CS),      pointer :: sum_output_CSp => NULL()
end type ocean_state_type

!> Control structure for this module
type MCT_MOM_Data

  type(ocean_state_type), pointer  :: ocn_state => NULL()   !< The private state of ocean
  type(ocean_public_type), pointer :: ocn_public => NULL()  !< The public state of ocean
  type(ocean_grid_type), pointer   :: grid => NULL()        !< The grid structure
  type(surface), pointer           :: ocn_surface => NULL() !< The ocean surface state
  type(forcing)                    :: fluxes                !< Structure that contains pointers to the
                                                            !! boundary forcing used to drive the liquid
                                                            !! ocean simulated by MOM.
  type(seq_infodata_type), pointer :: infodata              !< The input info type
  type(cpl_indices), public        :: ind                   !< Variable IDs
  ! runtime params
  logical :: sw_decomp   !< Controls whether shortwave is decomposed into four components
  real    :: c1, c2, c3, c4 !< Coeffs. used in the shortwave decomposition
  ! i/o
  character(len=384) :: pointer_filename !< Name of the ascii file that contains the path
                                         !! and filename of the latest restart file.
  integer :: stdout !< standard output unit. (by default, it should point to ocn.log.* file)
end type MCT_MOM_Data

type(MCT_MOM_Data) :: glb               !< global structure
integer :: id_clock_forcing

contains

!> This subroutine initializes MOM6.
subroutine ocn_init_mct( EClock, cdata_o, x2o_o, o2x_o, NLFilename )
  type(ESMF_Clock),             intent(inout) :: EClock      !< Time and time step ? \todo Why must this
                                                             !! be intent(inout)?
  type(seq_cdata)             , intent(inout) :: cdata_o     !< Input parameters
  type(mct_aVect)             , intent(inout) :: x2o_o       !< Fluxes from coupler to ocean, computed by ocean
  type(mct_aVect)             , intent(inout) :: o2x_o       !< Fluxes from ocean to coupler, computed by ocean
  character(len=*), optional  , intent(in)    :: NLFilename  !< Namelist filename

  !  local variables
  type(time_type)     :: time_init         !< Start time of coupled model's calendar
  type(time_type)     :: time_in           !< Time at the beginning of the first ocn coupling interval
  type(ESMF_time)     :: current_time      !< Current time
  type(ESMF_time)     :: time_in_ESMF      !< Initial time for ocean
  type(ESMF_timeInterval) :: ocn_cpl_interval !< Ocean coupling interval
  integer             :: year, month, day, hour, minute, seconds, seconds_n, seconds_d, rc
  character(len=240)  :: runid             !< Run ID
  character(len=240)  :: runtype           !< Run type
  character(len=240)  :: restartfile       !< Path/Name of restart file
  integer             :: nu                !< i/o unit to read pointer file
  character(len=240)  :: restart_pointer_file !< File name for restart pointer file
  character(len=240)  :: restartpath       !< Path of the restart file
  character(len=32)   :: starttype         !< infodata start type
  integer             :: mpicom_ocn        !< MPI ocn communicator
  integer             :: npes, pe0         !< # of processors and current processor
  integer             :: i, errorCode
  integer             :: lsize, nsend, nrecv
  logical             :: ldiag_cpl = .false.
  integer             :: isc, iec, jsc, jec, ni, nj     !< Indices for the start and end of the domain
                                                        !! in the x and y dir., respectively.
  ! runi-time params
  type(param_file_type) :: param_file           !< A structure to parse for run-time parameters
  type(directories)     :: dirs_tmp             !< A structure containing several relevant directory paths
  character(len=40)     :: mdl = "ocn_comp_mct" !< This module's name.

  ! mct variables (these are local for now)
  integer                   :: MOM_MCT_ID
  type(mct_gsMap), pointer  :: MOM_MCT_gsMap => NULL() !< 2d, points to cdata
  type(mct_gGrid), pointer  :: MOM_MCT_dom => NULL()   !< 2d, points to cdata
  type(mct_gsMap)           :: MOM_MCT_gsMap3d         !< for 3d streams, local
  type(mct_gGrid)           :: MOM_MCT_dom3d           !< for 3d streams, local

  ! time management
  integer                   :: ocn_cpl_dt   !< one ocn coupling interval in seconds. (to be received from cesm)
  real (kind=8)             :: mom_cpl_dt   !< one ocn coupling interval in seconds. (internal)
  real (kind=8), parameter  ::        &
      seconds_in_minute =    60.0d0, &
      seconds_in_hour   =  3600.0d0, &
      seconds_in_day    = 86400.0d0, &
      minutes_in_hour   =    60.0d0

  character(len=99) :: ocn_modelio_name !< ocn model input namelist filename
  integer           :: shrlogunit       !< original log file unit
  integer           :: shrloglev        !< original log level

  integer(kind=4)     :: inst_index     !< instance control vars (these are local for now)
  character(len=16)   :: inst_name
  character(len=16)   :: inst_suffix

  !!!DANGER!!!: change the following vars with the corresponding MOM6 vars
  integer :: km=1                   !< Number of vertical levels
  integer :: nx_block=0, ny_block=0 !< Size of block domain in x,y dir including ghost cells
  integer :: max_blocks_clinic=0    !< Max. number of blocks per processor in each distribution
  integer :: ncouple_per_day
  logical :: lsend_precip_fact      !< If T,send precip_fact to cpl for use in fw balance
                                    !! (partially-coupled option)
  character(len=128) :: err_msg     !< Error message

  ! set (actually, get from mct) the cdata pointers:
  call seq_cdata_setptrs(cdata_o, id=MOM_MCT_ID, mpicom=mpicom_ocn, &
                         gsMap=MOM_MCT_gsMap, dom=MOM_MCT_dom, infodata=glb%infodata)

  ! Determine attribute vector indices
  call coupler_indices_init(glb%ind)

  call seq_infodata_GetData( glb%infodata, case_name=runid )

  call seq_infodata_GetData( glb%infodata, start_type=starttype)

  if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
     runtype = "initial"
  else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
     runtype = "continue"
  else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
     runtype = "branch"
  else
     write(glb%stdout,*) 'ocn_comp_mct ERROR: unknown starttype'
     call exit(0)
  end if

  ! instance control
  inst_name   = seq_comm_name(MOM_MCT_ID)
  inst_index  = seq_comm_inst(MOM_MCT_ID)
  inst_suffix = seq_comm_suffix(MOM_MCT_ID)

  call t_startf('MOM_init')

  ! Initialize MOM6 comm
  call MOM_infra_init(mpicom_ocn)

  ! initialize ocn log file
  if (is_root_pe()) then

    ! get original log file properties
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)

    glb%stdout = shr_file_getUnit() ! get an unused unit number

    ! open the ocn_modelio.nml file and then open a log file associated with stdout
    ocn_modelio_name = 'ocn_modelio.nml' // trim(inst_suffix)
    call shr_file_setIO(ocn_modelio_name,glb%stdout)

    !  set the shr log io unit number
    call shr_file_setLogUnit(glb%stdout)
  end if

  call set_calendar_type(NOLEAP)  !TODO: confirm this

  ! Get the ESMF clock instance (assigned by CESM for MOM6)
  call ESMF_ClockGet(EClock, currTime=current_time, rc=rc)

  ! Get the initial CESM time
  call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
  time_init = set_date(year, month, day, hour, minute, seconds, err_msg=err_msg)

  ! Compute time_in: time at the beginning of the first ocn coupling interval
  call ESMF_ClockGet(EClock, TimeStep=ocn_cpl_interval, rc=rc)
  if (runtype /= "continue") then
    ! In startup runs, take the one ocn coupling interval lag into account to
    ! compute the initial ocn time.  (time_in = time_init + ocn_cpl_interval)
    time_in_ESMF = ESMF_TimeInc(current_time, ocn_cpl_interval)
  else
    time_in_ESMF = current_time
  endif
  call ESMF_TimeGet(time_in_ESMF, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
  time_in = set_date(year, month, day, hour, minute, seconds, err_msg=err_msg)

  ! Debugging clocks
  if (debug .and. is_root_pe()) then
    write(glb%stdout,*) 'ocn_init_mct, current time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, StartTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_init_mct, start time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, StopTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_init_mct, stop time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, PrevTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_init_mct, previous time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, TimeStep=ocn_cpl_interval, rc=rc)
    call ESMF_TimeIntervalGet(ocn_cpl_interval, yy=year, mm=month, d=day, s=seconds, sn=seconds_n, sd=seconds_d, rc=rc)
    write(glb%stdout,*) 'ocn_init_mct, time step: y,m,d-',year,month,day,'s,sn,sd=',seconds,seconds_n,seconds_d
  endif

  npes = num_pes()
  pe0 = root_pe()

  allocate(glb%ocn_public)
  glb%ocn_public%is_ocean_PE = .true.
  allocate(glb%ocn_public%pelist(npes))
  glb%ocn_public%pelist(:) = (/(i,i=pe0,pe0+npes)/)
  ! \todo Set other bits of glb$ocn_public

  ! This include declares and sets the variable "version".
  ! read useful runtime params
  call get_MOM_Input(param_file, dirs_tmp, check_params=.false.)
  !call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "POINTER_FILENAME", glb%pointer_filename, &
                 "Name of the ascii file that contains the path and filename of" // &
                 " the latest restart file.", default='rpointer.ocn')
  call get_param(param_file, mdl, "SW_DECOMP", glb%sw_decomp, &
                 "If True, read coeffs c1, c2, c3 and c4 and decompose" // &
                 "the net shortwave radiation (SW) into four components:\n" // &
                 "visible, direct shortwave  = c1 * SW \n" // &
                 "visible, diffuse shortwave = c2 * SW \n" // &
                 "near-IR, direct shortwave  = c3 * SW \n" // &
                 "near-IR, diffuse shortwave = c4 * SW", default=.true.)
  if (glb%sw_decomp) then
    call get_param(param_file, mdl, "SW_c1", glb%c1, &
                  "Coeff. used to convert net shortwave rad. into \n"//&
                  "visible, direct shortwave.", units="nondim", default=0.285)
    call get_param(param_file, mdl, "SW_c2", glb%c2, &
                  "Coeff. used to convert net shortwave rad. into \n"//&
                  "visible, diffuse shortwave.", units="nondim", default=0.285)
    call get_param(param_file, mdl, "SW_c3", glb%c3, &
                  "Coeff. used to convert net shortwave rad. into \n"//&
                  "near-IR, direct shortwave.", units="nondim", default=0.215)
    call get_param(param_file, mdl, "SW_c4", glb%c4, &
                  "Coeff. used to convert net shortwave rad. into \n"//&
                  "near-IR, diffuse shortwave.", units="nondim", default=0.215)
  else
    glb%c1 = 0.0; glb%c2 = 0.0; glb%c3 = 0.0; glb%c4 = 0.0
  endif

  ! Initialize the MOM6 model
  if (runtype == "initial") then ! startup (new run) - 'n' is needed below since we don't
                                 ! specify input_filename in input.nml
    call ocean_model_init(glb%ocn_public, glb%ocn_state, time_init, time_in, input_restart_file = 'n')
  else                           ! hybrid or branch or continuos runs
    ! output path root
    call seq_infodata_GetData( glb%infodata, outPathRoot=restartpath )
    ! read name of restart file in the pointer file
    nu = shr_file_getUnit()
    restart_pointer_file = trim(glb%pointer_filename)
    if (is_root_pe()) write(glb%stdout,*) 'Reading ocn pointer file: ',restart_pointer_file
    open(nu, file=restart_pointer_file, form='formatted', status='unknown')
    read(nu,'(a)') restartfile
    close(nu)
    !restartfile = trim(restartpath) // trim(restartfile)
    if (is_root_pe()) write(glb%stdout,*) 'Reading restart file: ',trim(restartfile)
    !endif
    call shr_file_freeUnit(nu)
    call ocean_model_init(glb%ocn_public, glb%ocn_state, time_init, time_in, input_restart_file=trim(restartfile))
  endif

  ! Initialize ocn_state%state out of sight
  call ocean_model_init_sfc(glb%ocn_state, glb%ocn_public)

  ! store pointers to components inside MOM
  call get_state_pointers(glb%ocn_state, grid=glb%grid)

  call t_stopf('MOM_init')

  ! Initialize MCT attribute vectors and indices
  call t_startf('MOM_mct_init')

  if (debug .and. root_pe().eq.pe_here()) print *, "calling ocn_SetGSMap_mct"

  ! Set mct global seg maps:

  call ocn_SetGSMap_mct(mpicom_ocn, MOM_MCT_ID, MOM_MCT_GSMap, MOM_MCT_GSMap3d)
  lsize = mct_gsMap_lsize(MOM_MCT_gsmap, mpicom_ocn)

  ! Initialize mct ocn domain (needs ocn initialization info)

  if (debug .and. root_pe().eq.pe_here()) print *, "calling ocn_domain_mct"
  call ocn_domain_mct(lsize, MOM_MCT_gsmap, MOM_MCT_dom)
  call ocn_domain_mct(lsize*km, MOM_MCT_gsmap3d, MOM_MCT_dom3d) !TODO: this is not used

  ! Inialize mct attribute vectors

  if (debug .and. root_pe().eq.pe_here()) print *, "calling mct_avect_init a"

  ! Initialize the mct attribute vector x2o_o, given Attribute list and length:
  call mct_aVect_init(x2o_o, rList=seq_flds_x2o_fields, lsize=lsize)
  ! set the mct attribute vector x2o_o to zero:
  call mct_aVect_zero(x2o_o)

  if (debug .and. root_pe().eq.pe_here()) print *, "calling mct_avect_init b"

  ! Initialize the mct attribute vector o2x_o, given Attribute list and length:
  call mct_aVect_init(o2x_o, rList=seq_flds_o2x_fields, lsize=lsize)
  ! set the mct attribute vector o2x_o to zero:
  call mct_aVect_zero(o2x_o)

  ! allocate send buffer
  nsend = mct_avect_nRattr(o2x_o)
  nrecv = mct_avect_nRattr(x2o_o)

  ! initialize necessary coupling info

  if (debug .and. root_pe().eq.pe_here()) print *, "calling seq_timemgr_eclockgetdata"

  call seq_timemgr_EClockGetData(EClock, dtime=ocn_cpl_dt)

  ! \todo Need interface to get dt from MOM6
  ncouple_per_day = seconds_in_day / ocn_cpl_dt
  mom_cpl_dt = seconds_in_day / ncouple_per_day
  if (mom_cpl_dt /= ocn_cpl_dt) then
     write(glb%stdout,*) 'ERROR mom_cpl_dt and ocn_cpl_dt must be identical'
     call exit(0)
  end if

  ! send initial state to driver

  !TODO:
  ! if ( lsend_precip_fact )  then
  !    call seq_infodata_PutData( infodata, precip_fact=precip_fact)
  ! end if


  if (debug .and. root_pe().eq.pe_here()) print *, "calling ocn_export"
  call ocn_export(glb%ind, glb%ocn_public, glb%grid, o2x_o%rattr)

  call t_stopf('MOM_mct_init')

  if (debug .and. root_pe().eq.pe_here()) print *, "calling get_state_pointers"

  ! Size of global domain
  call get_global_grid_size(glb%grid, ni, nj)

  if (debug .and. root_pe().eq.pe_here()) print *, "calling seq_infodata_putdata"

   call seq_infodata_PutData( glb%infodata, &
        ocn_nx = ni , ocn_ny = nj)
   call seq_infodata_PutData( glb%infodata, &
        ocn_prognostic=.true., ocnrof_prognostic=.true.)

  if (debug .and. root_pe().eq.pe_here()) print *, "leaving ocean_init_mct"

  ! Reset shr logging to original values
  if (is_root_pe()) then
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
  end if

end subroutine ocn_init_mct

!> Determines attribute vector indices
subroutine coupler_indices_init(ind)

  type(cpl_indices), intent(inout) :: ind !< Structure with coupler indices
                                          !! and vectors

  ! Local Variables
  type(mct_aVect) :: o2x      !< Array with ocean to coupler data
  type(mct_aVect) :: x2o      !< Array with coupler to ocean data

  integer          :: ncat    !< Thickness category index
  character(len=2) :: cncat   !< Character version of ncat
  integer          :: ncol    !< Column index
  integer          :: mcog_ncols !< Number of ice thickness categories?
  integer          :: lmcog_flds_sent !< Used to convert per thickness
                                      !! category fields?

  ! create temporary attribute vectors
  call mct_aVect_init(x2o, rList=seq_flds_x2o_fields, lsize=1)
  call mct_aVect_init(o2x, rList=seq_flds_o2x_fields, lsize=1)

  ! ocean to coupler
  ind%o2x_So_t          = mct_avect_indexra(o2x,'So_t')
  ind%o2x_So_u          = mct_avect_indexra(o2x,'So_u')
  ind%o2x_So_v          = mct_avect_indexra(o2x,'So_v')
  ind%o2x_So_s          = mct_avect_indexra(o2x,'So_s')
  ind%o2x_So_dhdx       = mct_avect_indexra(o2x,'So_dhdx')
  ind%o2x_So_dhdy       = mct_avect_indexra(o2x,'So_dhdy')
  ! QL, 150526, to wav, boundary layer depth
  ind%o2x_So_bldepth    = mct_avect_indexra(o2x,'So_bldepth')
  ind%o2x_Fioo_q        = mct_avect_indexra(o2x,'Fioo_q')
  ind%o2x_Faoo_fco2_ocn = mct_avect_indexra(o2x,'Faoo_fco2_ocn',perrWith='quiet')
  ind%o2x_Faoo_fdms_ocn = mct_avect_indexra(o2x,'Faoo_fdms_ocn',perrWith='quiet')

  ! coupler to ocean
  ind%x2o_Si_ifrac      = mct_avect_indexra(x2o,'Si_ifrac')
  ind%x2o_Sa_pslv       = mct_avect_indexra(x2o,'Sa_pslv')
  ind%x2o_So_duu10n     = mct_avect_indexra(x2o,'So_duu10n')
  ! QL, 150526, from wav
  ind%x2o_Sw_lamult     = mct_avect_indexra(x2o,'Sw_lamult')
  ind%x2o_Sw_ustokes    = mct_avect_indexra(x2o,'Sw_ustokes')
  ind%x2o_Sw_vstokes    = mct_avect_indexra(x2o,'Sw_vstokes')
  ind%x2o_Foxx_tauy     = mct_avect_indexra(x2o,'Foxx_tauy')
  ind%x2o_Foxx_taux     = mct_avect_indexra(x2o,'Foxx_taux')
  ind%x2o_Foxx_swnet    = mct_avect_indexra(x2o,'Foxx_swnet')
  ind%x2o_Foxx_lat      = mct_avect_indexra(x2o,'Foxx_lat')
  ind%x2o_Foxx_sen      = mct_avect_indexra(x2o,'Foxx_sen')
  ind%x2o_Foxx_lwup     = mct_avect_indexra(x2o,'Foxx_lwup')
  ind%x2o_Faxa_lwdn     = mct_avect_indexra(x2o,'Faxa_lwdn')
  ind%x2o_Fioi_melth    = mct_avect_indexra(x2o,'Fioi_melth')
  ind%x2o_Fioi_meltw    = mct_avect_indexra(x2o,'Fioi_meltw')
  ind%x2o_Fioi_salt     = mct_avect_indexra(x2o,'Fioi_salt')
  ind%x2o_Fioi_bcpho    = mct_avect_indexra(x2o,'Fioi_bcpho')
  ind%x2o_Fioi_bcphi    = mct_avect_indexra(x2o,'Fioi_bcphi')
  ind%x2o_Fioi_flxdst   = mct_avect_indexra(x2o,'Fioi_flxdst')
  ind%x2o_Faxa_prec     = mct_avect_indexra(x2o,'Faxa_prec')
  ind%x2o_Faxa_snow     = mct_avect_indexra(x2o,'Faxa_snow')
  ind%x2o_Faxa_rain     = mct_avect_indexra(x2o,'Faxa_rain')
  ind%x2o_Foxx_evap     = mct_avect_indexra(x2o,'Foxx_evap')
  ind%x2o_Foxx_rofl     = mct_avect_indexra(x2o,'Foxx_rofl')
  ind%x2o_Foxx_rofi     = mct_avect_indexra(x2o,'Foxx_rofi')
  ind%x2o_Faxa_bcphidry = mct_avect_indexra(x2o,'Faxa_bcphidry')
  ind%x2o_Faxa_bcphodry = mct_avect_indexra(x2o,'Faxa_bcphodry')
  ind%x2o_Faxa_bcphiwet = mct_avect_indexra(x2o,'Faxa_bcphiwet')
  ind%x2o_Faxa_ocphidry = mct_avect_indexra(x2o,'Faxa_ocphidry')
  ind%x2o_Faxa_ocphodry = mct_avect_indexra(x2o,'Faxa_ocphodry')
  ind%x2o_Faxa_ocphiwet = mct_avect_indexra(x2o,'Faxa_ocphiwet')
  ind%x2o_Faxa_dstdry1  = mct_avect_indexra(x2o,'Faxa_dstdry1')
  ind%x2o_Faxa_dstdry2  = mct_avect_indexra(x2o,'Faxa_dstdry2')
  ind%x2o_Faxa_dstdry3  = mct_avect_indexra(x2o,'Faxa_dstdry3')
  ind%x2o_Faxa_dstdry4  = mct_avect_indexra(x2o,'Faxa_dstdry4')
  ind%x2o_Faxa_dstwet1  = mct_avect_indexra(x2o,'Faxa_dstwet1')
  ind%x2o_Faxa_dstwet2  = mct_avect_indexra(x2o,'Faxa_dstwet2')
  ind%x2o_Faxa_dstwet3  = mct_avect_indexra(x2o,'Faxa_dstwet3')
  ind%x2o_Faxa_dstwet4  = mct_avect_indexra(x2o,'Faxa_dstwet4')
  ind%x2o_Sa_co2prog    = mct_avect_indexra(x2o,'Sa_co2prog',perrWith='quiet')
  ind%x2o_Sa_co2diag    = mct_avect_indexra(x2o,'Sa_co2diag',perrWith='quiet')
  ! optional per thickness category fields

  ! convert cpl indices to mcog column indices
  ! this implementation only handles columns due to ice thickness categories
  lmcog_flds_sent = seq_flds_i2o_per_cat

  if (seq_flds_i2o_per_cat) then
    mcog_ncols = ice_ncat+1
    allocate(ind%x2o_frac_col(mcog_ncols))
    allocate(ind%x2o_fracr_col(mcog_ncols))
    allocate(ind%x2o_qsw_fracr_col(mcog_ncols))
    ncol = 1
    ind%x2o_frac_col(ncol)        = mct_avect_indexra(x2o,'Sf_afrac')
    ind%x2o_fracr_col(ncol)       = mct_avect_indexra(x2o,'Sf_afracr')
    ind%x2o_qsw_fracr_col(ncol)   = mct_avect_indexra(x2o,'Foxx_swnet_afracr')

    do ncat = 1, ice_ncat
      write(cncat,'(i2.2)') ncat
      ncol = ncat+1
      ind%x2o_frac_col(ncol)      = mct_avect_indexra(x2o,'Si_ifrac_'//cncat)
      ind%x2o_fracr_col(ncol)     = ind%x2o_frac_col(ncol)
      ind%x2o_qsw_fracr_col(ncol) = mct_avect_indexra(x2o,'PFioi_swpen_ifrac_'//cncat)
    enddo
  else
    mcog_ncols = 1
  endif

  call mct_aVect_clean(x2o)
  call mct_aVect_clean(o2x)

end subroutine coupler_indices_init

!> Initializes the ocean model, including registering fields
!! for restarts and reading restart files if appropriate.
subroutine ocean_model_init(Ocean_sfc, OS, Time_init, Time_in, gas_fields_ocn, input_restart_file)
  type(ocean_public_type), target, &
                       intent(inout) :: Ocean_sfc !< A structure containing various
                                !! publicly visible ocean surface properties after initialization,
                                !! the data in this type is intent(out).
  type(ocean_state_type), pointer    :: OS        !< A structure whose internal
                                !! contents are private to ocean_model_mod that may be used to
                                !! contain all information about the ocean's interior state.
  type(time_type),     intent(in)    :: Time_init !< The start time for the coupled model's calendar
  type(time_type),     intent(in)    :: Time_in   !< The time at which to initialize the ocean model.
  type(coupler_1d_bc_type), &
             optional, intent(in)    :: gas_fields_ocn !< If present, this type describes the
                                              !! ocean and surface-ice fields that will participate
                                              !! in the calculation of additional gas or other
                                              !! tracer fluxes, and can be used to spawn related
                                              !! internal variables in the ice model.
  character(len=*), optional, intent(in) :: input_restart_file !< If present, name of restart file to read

! This subroutine initializes both the ocean state and the ocean surface type.
! Because of the way that indicies and domains are handled, Ocean_sfc must have
! been used in a previous call to initialize_ocean_type.

  real :: Time_unit   !< The time unit in seconds for ENERGYSAVEDAYS.
  real :: Rho0        !< The Boussinesq ocean density, in kg m-3.
  real :: G_Earth     !< The gravitational acceleration in m s-2.
                      !! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "ocean_model_init"  !< This module's name.
  character(len=48)  :: stagger
  integer :: secs, days
  type(param_file_type) :: param_file !< A structure to parse for run-time parameters
  logical :: offline_tracer_mode

  call callTree_enter("ocean_model_init(), ocn_comp_mct.F90")
  if (associated(OS)) then
    call MOM_error(WARNING, "ocean_model_init called with an associated "// &
                    "ocean_state_type structure. Model is already initialized.")
    return
  endif
  allocate(OS)

  OS%is_ocean_pe = Ocean_sfc%is_ocean_pe
  if (.not.OS%is_ocean_pe) return

  OS%Time = Time_in
  call initialize_MOM(OS%Time, param_file, OS%dirs, OS%MOM_CSp, Time_in, &
      offline_tracer_mode=offline_tracer_mode, input_restart_file=input_restart_file)
  OS%grid => OS%MOM_CSp%G ; OS%GV => OS%MOM_CSp%GV
  OS%C_p = OS%MOM_CSp%tv%C_p
  OS%fluxes%C_p = OS%MOM_CSp%tv%C_p

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "RESTART_CONTROL", OS%Restart_control, &
                 "An integer whose bits encode which restart files are \n"//&
                 "written. Add 2 (bit 1) for a time-stamped file, and odd \n"//&
                 "(bit 0) for a non-time-stamped file.  A restart file \n"//&
                 "will be saved at the end of the run segment for any \n"//&
                 "non-negative value.", default=1)
  call get_param(param_file, mdl, "TIMEUNIT", Time_unit, &
                 "The time unit for ENERGYSAVEDAYS.", &
                 units="s", default=86400.0)
  call get_param(param_file, mdl, "ENERGYSAVEDAYS",OS%energysavedays, &
                 "The interval in units of TIMEUNIT between saves of the \n"//&
                 "energies of the run and other globally summed diagnostics.", &
                 default=set_time(0,days=1), timeunit=Time_unit)

  call get_param(param_file, mdl, "OCEAN_SURFACE_STAGGER", stagger, &
                 "A case-insensitive character string to indicate the \n"//&
                 "staggering of the surface velocity field that is \n"//&
                 "returned to the coupler.  Valid values include \n"//&
                 "'A', 'B', or 'C'.", default="C")
  if (uppercase(stagger(1:1)) == 'A') then ; Ocean_sfc%stagger = AGRID
  elseif (uppercase(stagger(1:1)) == 'B') then ; Ocean_sfc%stagger = BGRID_NE
  elseif (uppercase(stagger(1:1)) == 'C') then ; Ocean_sfc%stagger = CGRID_NE
  else ; call MOM_error(FATAL,"ocean_model_init: OCEAN_SURFACE_STAGGER = "// &
                        trim(stagger)//" is invalid.") ; endif

  call get_param(param_file, mdl, "RESTORE_SALINITY",OS%restore_salinity, &
                 "If true, the coupled driver will add a globally-balanced \n"//&
                 "fresh-water flux that drives sea-surface salinity \n"//&
                 "toward specified values.", default=.false.)
  call get_param(param_file, mdl, "RESTORE_TEMPERATURE",OS%restore_temp, &
                 "If true, the coupled driver will add a  \n"//&
                 "heat flux that drives sea-surface temperauture \n"//&
                 "toward specified values.", default=.false.)
  call get_param(param_file, mdl, "RHO_0", Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  call get_param(param_file, mdl, "G_EARTH", G_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)

  call get_param(param_file, mdl, "ICE_SHELF",  OS%use_ice_shelf, &
                 "If true, enables the ice shelf model.", default=.false.)

  call get_param(param_file, mdl, "ICEBERGS_APPLY_RIGID_BOUNDARY",  OS%icebergs_apply_rigid_boundary, &
                 "If true, allows icebergs to change boundary condition felt by ocean", default=.false.)

  if (OS%icebergs_apply_rigid_boundary) then
    call get_param(param_file, mdl, "KV_ICEBERG",  OS%kv_iceberg, &
                 "The viscosity of the icebergs",  units="m2 s-1",default=1.0e10)
    call get_param(param_file, mdl, "DENSITY_ICEBERGS",  OS%density_iceberg, &
                 "A typical density of icebergs.", units="kg m-3", default=917.0)
    call get_param(param_file, mdl, "LATENT_HEAT_FUSION", OS%latent_heat_fusion, &
                 "The latent heat of fusion.", units="J/kg", default=hlf)
    call get_param(param_file, mdl, "BERG_AREA_THRESHOLD", OS%berg_area_threshold, &
                 "Fraction of grid cell which iceberg must occupy, so that fluxes \n"//&
                 "below berg are set to zero. Not applied for negative \n"//&
                 " values.", units="non-dim", default=-1.0)
  endif

  OS%press_to_z = 1.0/(Rho0*G_Earth)

  !   Consider using a run-time flag to determine whether to do the diagnostic
  ! vertical integrals, since the related 3-d sums are not negligible in cost.
  call allocate_surface_state(OS%state, OS%grid, OS%MOM_CSp%use_temperature, &
                              do_integrals=.true., gas_fields_ocn=gas_fields_ocn)

  call surface_forcing_init(Time_in, OS%grid, param_file, OS%MOM_CSp%diag, &
                            OS%forcing_CSp, OS%restore_salinity, OS%restore_temp)

  if (OS%use_ice_shelf)  then
    call initialize_ice_shelf(param_file, OS%grid, OS%Time, OS%ice_shelf_CSp, &
                              OS%MOM_CSp%diag, OS%fluxes)
  endif
  if (OS%icebergs_apply_rigid_boundary)  then
    !call allocate_forcing_type(OS%grid, OS%fluxes, iceberg=.true.)
    !This assumes that the iceshelf and ocean are on the same grid. I hope this is true
    if (.not. OS%use_ice_shelf) call allocate_forcing_type(OS%grid, OS%fluxes, ustar=.true., shelf=.true.)
  endif

  call MOM_sum_output_init(OS%grid, param_file, OS%dirs%output_directory, &
                            OS%MOM_CSp%ntrunc, Time_init, OS%sum_output_CSp)

  ! This call has been moved into the first call to update_ocean_model.
!  call write_energy(OS%MOM_CSp%u, OS%MOM_CSp%v, OS%MOM_CSp%h, OS%MOM_CSp%tv, &
!             OS%Time, 0, OS%grid, OS%GV, OS%sum_output_CSp, OS%MOM_CSp%tracer_flow_CSp)

  ! write_energy_time is the next integral multiple of energysavedays.
  OS%write_energy_time = Time_init + OS%energysavedays * &
                         (1 + (OS%Time - Time_init) / OS%energysavedays)

  if (ASSOCIATED(OS%grid%Domain%maskmap)) then
    call initialize_ocean_public_type(OS%grid%Domain%mpp_domain, Ocean_sfc, &
                                      OS%MOM_CSp%diag, maskmap=OS%grid%Domain%maskmap, &
                                      gas_fields_ocn=gas_fields_ocn)
  else
    call initialize_ocean_public_type(OS%grid%Domain%mpp_domain, Ocean_sfc, &
                                      OS%MOM_CSp%diag, gas_fields_ocn=gas_fields_ocn)
  endif

  ! This call can only occur here if the coupler_bc_type variables have been
  ! initialized already using the information from gas_fields_ocn.
  if (present(gas_fields_ocn)) then
    call calculate_surface_state(OS%state, OS%MOM_CSp%u, &
             OS%MOM_CSp%v, OS%MOM_CSp%h, OS%MOM_CSp%ave_ssh,&
             OS%grid, OS%GV, OS%MOM_CSp)

    call convert_state_to_ocean_type(OS%state, Ocean_sfc, OS%grid, &
                                     OS%MOM_CSp%use_conT_absS)
  endif

  call close_param_file(param_file)
  call diag_mediator_close_registration(OS%MOM_CSp%diag)

  if (is_root_pe()) &
    write(glb%stdout,'(/12x,a/)') '======== COMPLETED MOM INITIALIZATION ========'

  call callTree_leave("ocean_model_init(")
end subroutine ocean_model_init

!> Extracts the surface properties from the ocean's internal
!! state and stores them in the ocean type returned to the calling ice model.
!! It has to be separate from the ocean_initialization call because the coupler
!! module allocates the space for some of these variables.
subroutine ocean_model_init_sfc(OS, Ocean_sfc)
  type(ocean_state_type),  pointer       :: OS
  type(ocean_public_type), intent(inout) :: Ocean_sfc

  integer :: is, ie, js, je

  is = OS%grid%isc ; ie = OS%grid%iec ; js = OS%grid%jsc ; je = OS%grid%jec
  call coupler_type_spawn(Ocean_sfc%fields, OS%state%tr_fields, &
                          (/is,is,ie,ie/), (/js,js,je,je/), as_needed=.true.)

  call calculate_surface_state(OS%state, OS%MOM_CSp%u, &
           OS%MOM_CSp%v, OS%MOM_CSp%h, OS%MOM_CSp%ave_ssh,&
           OS%grid, OS%GV, OS%MOM_CSp)

  call convert_state_to_ocean_type(OS%state, Ocean_sfc, OS%grid, &
                                   OS%MOM_CSp%use_conT_absS)

end subroutine ocean_model_init_sfc

!> Initializes surface forcing: get relevant parameters and allocate arrays.
subroutine surface_forcing_init(Time, G, param_file, diag, CS, restore_salt, restore_temp)
  type(time_type),          intent(in)    :: Time !< The current model time
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(param_file_type),    intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,  intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(surface_forcing_CS), pointer       :: CS !< A pointer that is set to point to the
                                                !! control structure for this module
  logical, optional,        intent(in)    :: restore_salt, restore_temp !< If present and true,
                                             !! temp/salt restoring will be applied

  ! local variables
  real :: utide  !< The RMS tidal velocity, in m s-1.
  type(directories)  :: dirs
  logical            :: new_sim, iceberg_flux_diags
  type(time_type)    :: Time_frc
  character(len=200) :: TideAmp_file, gust_file, salt_file, temp_file ! Input file names.
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "ocn_comp_mct"  ! This module's name.
  character(len=48)  :: stagger
  character(len=240) :: basin_file
  integer :: i, j, isd, ied, jsd, jed

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (associated(CS)) then
    call MOM_error(WARNING, "surface_forcing_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  id_clock_forcing=cpu_clock_id('Ocean surface forcing', grain=CLOCK_SUBCOMPONENT)
  call cpu_clock_begin(id_clock_forcing)

  CS%diag => diag

  call write_version_number (version)
  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, &
                 "The directory in which all input files are found.", &
                 default=".")
  CS%inputdir = slasher(CS%inputdir)
  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state \n"//&
                 "variables.", default=.true.)
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  call get_param(param_file, mdl, "LATENT_HEAT_FUSION", CS%latent_heat_fusion, &
                 "The latent heat of fusion.", units="J/kg", default=hlf)
  call get_param(param_file, mdl, "LATENT_HEAT_VAPORIZATION", CS%latent_heat_vapor, &
                 "The latent heat of fusion.", units="J/kg", default=hlv)
  call get_param(param_file, mdl, "MAX_P_SURF", CS%max_p_surf, &
                 "The maximum surface pressure that can be exerted by the \n"//&
                 "atmosphere and floating sea-ice or ice shelves. This is \n"//&
                 "needed because the FMS coupling structure does not \n"//&
                 "limit the water that can be frozen out of the ocean and \n"//&
                 "the ice-ocean heat fluxes are treated explicitly.  No \n"//&
                 "limit is applied if a negative value is used.", units="Pa", &
                 default=-1.0)
  call get_param(param_file, mdl, "ADJUST_NET_SRESTORE_TO_ZERO", &
                 CS%adjust_net_srestore_to_zero, &
                 "If true, adjusts the salinity restoring seen to zero\n"//&
                 "whether restoring is via a salt flux or virtual precip.",&
                 default=restore_salt)
  call get_param(param_file, mdl, "ADJUST_NET_SRESTORE_BY_SCALING", &
                 CS%adjust_net_srestore_by_scaling, &
                 "If true, adjustments to salt restoring to achieve zero net are\n"//&
                 "made by scaling values without moving the zero contour.",&
                 default=.false.)
  call get_param(param_file, mdl, "ADJUST_NET_FRESH_WATER_TO_ZERO", &
                 CS%adjust_net_fresh_water_to_zero, &
                 "If true, adjusts the net fresh-water forcing seen \n"//&
                 "by the ocean (including restoring) to zero.", default=.false.)
  call get_param(param_file, mdl, "ADJUST_NET_FRESH_WATER_BY_SCALING", &
                 CS%adjust_net_fresh_water_by_scaling, &
                 "If true, adjustments to net fresh water to achieve zero net are\n"//&
                 "made by scaling values without moving the zero contour.",&
                 default=.false.)
  call get_param(param_file, mdl, "ICE_SALT_CONCENTRATION", &
                 CS%ice_salt_concentration, &
                 "The assumed sea-ice salinity needed to reverse engineer the \n"//&
                 "melt flux (or ice-ocean fresh-water flux).", &
                 units="kg/kg", default=0.005)
  call get_param(param_file, mdl, "USE_LIMITED_PATM_SSH", CS%use_limited_P_SSH, &
                 "If true, return the sea surface height with the \n"//&
                 "correction for the atmospheric (and sea-ice) pressure \n"//&
                 "limited by max_p_surf instead of the full atmospheric \n"//&
                 "pressure.", default=.true.)

! smg: should get_param call should be removed when have A=B code reconciled.
! this param is used to distinguish how to diagnose surface heat content from water.
  call get_param(param_file, mdl, "BULKMIXEDLAYER", CS%bulkmixedlayer, &
                 default=CS%use_temperature,do_not_log=.true.)

  call get_param(param_file, mdl, "WIND_STAGGER", stagger, &
                 "A case-insensitive character string to indicate the \n"//&
                 "staggering of the input wind stress field.  Valid \n"//&
                 "values are 'A', 'B', or 'C'.", default="C")
  if (uppercase(stagger(1:1)) == 'A') then ; CS%wind_stagger = AGRID
  elseif (uppercase(stagger(1:1)) == 'B') then ; CS%wind_stagger = BGRID_NE
  elseif (uppercase(stagger(1:1)) == 'C') then ; CS%wind_stagger = CGRID_NE
  else ; call MOM_error(FATAL,"surface_forcing_init: WIND_STAGGER = "// &
                        trim(stagger)//" is invalid.") ; endif
  call get_param(param_file, mdl, "WIND_STRESS_MULTIPLIER", CS%wind_stress_multiplier, &
                 "A factor multiplying the wind-stress given to the ocean by the\n"//&
                 "coupler. This is used for testing and should be =1.0 for any\n"//&
                 "production runs.", default=1.0)

  if (restore_salt) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes \n"//&
                 "to the relative surface anomalies (akin to a piston \n"//&
                 "velocity).  Note the non-MKS units.", units="m day-1", &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, "SALT_RESTORE_FILE", CS%salt_restore_file, &
                 "A file in which to find the surface salinity to use for restoring.", &
                 default="salt_restore.nc")
    call get_param(param_file, mdl, "SALT_RESTORE_VARIABLE", CS%salt_restore_var_name, &
                 "The name of the surface salinity variable to read from "//&
                 "SALT_RESTORE_FILE for restoring salinity.", &
                 default="salt")
! Convert CS%Flux_const from m day-1 to m s-1.
    CS%Flux_const = CS%Flux_const / 86400.0

    call get_param(param_file, mdl, "SRESTORE_AS_SFLUX", CS%salt_restore_as_sflux, &
                 "If true, the restoring of salinity is applied as a salt \n"//&
                 "flux instead of as a freshwater flux.", default=.false.)
    call get_param(param_file, mdl, "MAX_DELTA_SRESTORE", CS%max_delta_srestore, &
                 "The maximum salinity difference used in restoring terms.", &
                 units="PSU or g kg-1", default=999.0)
    call get_param(param_file, mdl, "MASK_SRESTORE_UNDER_ICE", &
                 CS%mask_srestore_under_ice, &
                 "If true, disables SSS restoring under sea-ice based on a frazil\n"//&
                 "criteria (SST<=Tf). Only used when RESTORE_SALINITY is True.",      &
                 default=.false.)
    call get_param(param_file, mdl, "MASK_SRESTORE_MARGINAL_SEAS", &
                 CS%mask_srestore_marginal_seas, &
                 "If true, disable SSS restoring in marginal seas. Only used when\n"//&
                 "RESTORE_SALINITY is True.", default=.false.)
    call get_param(param_file, mdl, "BASIN_FILE", basin_file, &
                 "A file in which to find the basin masks, in variable 'basin'.", &
                 default="basin.nc")
    basin_file = trim(CS%inputdir) // trim(basin_file)
    call safe_alloc_ptr(CS%basin_mask,isd,ied,jsd,jed) ; CS%basin_mask(:,:) = 1.0
    if (CS%mask_srestore_marginal_seas) then
      call read_data(basin_file,'basin',CS%basin_mask,domain=G%domain%mpp_domain,timelevel=1)
      do j=jsd,jed ; do i=isd,ied
        if (CS%basin_mask(i,j) >= 6.0) then ; CS%basin_mask(i,j) = 0.0
        else ; CS%basin_mask(i,j) = 1.0 ; endif
      enddo ; enddo
    endif
  endif

  if (restore_temp) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes \n"//&
                 "to the relative surface anomalies (akin to a piston \n"//&
                 "velocity).  Note the non-MKS units.", units="m day-1", &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, "SST_RESTORE_FILE", CS%temp_restore_file, &
                 "A file in which to find the surface temperature to use for restoring.", &
                 default="temp_restore.nc")
    call get_param(param_file, mdl, "SST_RESTORE_VARIABLE", CS%temp_restore_var_name, &
                 "The name of the surface temperature variable to read from "//&
                 "SST_RESTORE_FILE for restoring sst.", &
                 default="temp")
! Convert CS%Flux_const from m day-1 to m s-1.
    CS%Flux_const = CS%Flux_const / 86400.0

    call get_param(param_file, mdl, "MAX_DELTA_TRESTORE", CS%max_delta_trestore, &
                 "The maximum sst difference used in restoring terms.", &
                 units="degC ", default=999.0)

  endif

! Optionally read tidal amplitude from input file (m s-1) on model grid.
! Otherwise use default tidal amplitude for bottom frictionally-generated
! dissipation. Default cd_tides is chosen to yield approx 1 TWatt of
! work done against tides globally using OSU tidal amplitude.
  call get_param(param_file, mdl, "CD_TIDES", CS%cd_tides, &
                 "The drag coefficient that applies to the tides.", &
                 units="nondim", default=1.0e-4)
  call get_param(param_file, mdl, "READ_TIDEAMP", CS%read_TIDEAMP, &
                 "If true, read a file (given by TIDEAMP_FILE) containing \n"//&
                 "the tidal amplitude with INT_TIDE_DISSIPATION.", default=.false.)
  if (CS%read_TIDEAMP) then
    call get_param(param_file, mdl, "TIDEAMP_FILE", TideAmp_file, &
                 "The path to the file containing the spatially varying \n"//&
                 "tidal amplitudes with INT_TIDE_DISSIPATION.", &
                 default="tideamp.nc")
    CS%utide=0.0
  else
    call get_param(param_file, mdl, "UTIDE", CS%utide, &
                 "The constant tidal amplitude used with INT_TIDE_DISSIPATION.", &
                 units="m s-1", default=0.0)
  endif

  call safe_alloc_ptr(CS%TKE_tidal,isd,ied,jsd,jed)
  call safe_alloc_ptr(CS%ustar_tidal,isd,ied,jsd,jed)

  if (CS%read_TIDEAMP) then
    TideAmp_file = trim(CS%inputdir) // trim(TideAmp_file)
    call read_data(TideAmp_file,'tideamp',CS%TKE_tidal,domain=G%domain%mpp_domain,timelevel=1)
    do j=jsd, jed; do i=isd, ied
      utide = CS%TKE_tidal(i,j)
      CS%TKE_tidal(i,j) = G%mask2dT(i,j)*CS%Rho0*CS%cd_tides*(utide*utide*utide)
      CS%ustar_tidal(i,j)=sqrt(CS%cd_tides)*utide
    enddo ; enddo
  else
    do j=jsd,jed; do i=isd,ied
      utide=CS%utide
      CS%TKE_tidal(i,j) = CS%Rho0*CS%cd_tides*(utide*utide*utide)
      CS%ustar_tidal(i,j)=sqrt(CS%cd_tides)*utide
    enddo ; enddo
  endif

  call time_interp_external_init

! Optionally read a x-y gustiness field in place of a global
! constant.

  call get_param(param_file, mdl, "READ_GUST_2D", CS%read_gust_2d, &
                 "If true, use a 2-dimensional gustiness supplied from \n"//&
                 "an input file", default=.false.)
  call get_param(param_file, mdl, "GUST_CONST", CS%gust_const, &
               "The background gustiness in the winds.", units="Pa", &
               default=0.02)
  if (CS%read_gust_2d) then
    call get_param(param_file, mdl, "GUST_2D_FILE", gust_file, &
                 "The file in which the wind gustiness is found in \n"//&
                 "variable gustiness.")

    call safe_alloc_ptr(CS%gust,isd,ied,jsd,jed)
    gust_file = trim(CS%inputdir) // trim(gust_file)
    call read_data(gust_file,'gustiness',CS%gust,domain=G%domain%mpp_domain, &
                   timelevel=1) ! units should be Pa
  endif

! See whether sufficiently thick sea ice should be treated as rigid.
  call get_param(param_file, mdl, "USE_RIGID_SEA_ICE", CS%rigid_sea_ice, &
                 "If true, sea-ice is rigid enough to exert a \n"//&
                 "nonhydrostatic pressure that resist vertical motion.", &
                 default=.false.)
  if (CS%rigid_sea_ice) then
    call get_param(param_file, mdl, "SEA_ICE_MEAN_DENSITY", CS%density_sea_ice, &
                 "A typical density of sea ice, used with the kinematic \n"//&
                 "viscosity, when USE_RIGID_SEA_ICE is true.", units="kg m-3", &
                 default=900.0)
    call get_param(param_file, mdl, "SEA_ICE_VISCOSITY", CS%Kv_sea_ice, &
                 "The kinematic viscosity of sufficiently thick sea ice \n"//&
                 "for use in calculating the rigidity of sea ice.", &
                 units="m2 s-1", default=1.0e9)
    call get_param(param_file, mdl, "SEA_ICE_RIGID_MASS", CS%rigid_sea_ice_mass, &
                 "The mass of sea-ice per unit area at which the sea-ice \n"//&
                 "starts to exhibit rigidity", units="kg m-2", default=1000.0)
  endif

  call get_param(param_file, mdl, "ALLOW_ICEBERG_FLUX_DIAGNOSTICS", iceberg_flux_diags, &
                 "If true, makes available diagnostics of fluxes from icebergs\n"//&
                 "as seen by MOM6.", default=.false.)
  call register_forcing_type_diags(Time, diag, CS%use_temperature, CS%handles, &
                                   use_berg_fluxes=iceberg_flux_diags)

  call get_param(param_file, mdl, "ALLOW_FLUX_ADJUSTMENTS", CS%allow_flux_adjustments, &
                 "If true, allows flux adjustments to specified via the \n"//&
                 "data_table using the component name 'OCN'.", default=.false.)
  if (CS%allow_flux_adjustments) then
    call data_override_init(Ocean_domain_in=G%Domain%mpp_domain)
  endif

  if (present(restore_salt)) then ; if (restore_salt) then
    salt_file = trim(CS%inputdir) // trim(CS%salt_restore_file)
    CS%id_srestore = init_external_field(salt_file, CS%salt_restore_var_name, domain=G%Domain%mpp_domain)
  endif ; endif

  if (present(restore_temp)) then ; if (restore_temp) then
    temp_file = trim(CS%inputdir) // trim(CS%temp_restore_file)
    CS%id_trestore = init_external_field(temp_file, CS%temp_restore_var_name, domain=G%Domain%mpp_domain)
  endif ; endif

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

!###  call controlled_forcing_init(Time, G, param_file, diag, CS%ctrl_forcing_CSp)

  call user_revise_forcing_init(param_file, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)
end subroutine surface_forcing_init

!> Initializes domain and state variables contained in the ocean public type.
subroutine initialize_ocean_public_type(input_domain, Ocean_sfc, diag, maskmap, &
                                        gas_fields_ocn)
  type(domain2D), intent(in)             :: input_domain !< The FMS domain for the input structure
  type(ocean_public_type), intent(inout) :: Ocean_sfc    !< Ocean surface state
  type(diag_ctrl), intent(in)            :: diag         !< A structure used to control diagnostics.
  logical, intent(in), optional          :: maskmap(:,:) !< A pointer to an array indicating which
                                         !! logical processors are actually used for the ocean code.
  type(coupler_1d_bc_type), &
                 optional, intent(in)    :: gas_fields_ocn !< If present, this type describes the
                                         !! ocean and surface-ice fields that will participate
                                         !! in the calculation of additional gas or other
                                         !! tracer fluxes.
  ! local variables
  integer :: xsz, ysz, layout(2)
  integer :: isc, iec, jsc, jec

  call mpp_get_layout(input_domain,layout)
  call mpp_get_global_domain(input_domain, xsize=xsz, ysize=ysz)
  if(PRESENT(maskmap)) then
     call mpp_define_domains((/1,xsz,1,ysz/),layout,Ocean_sfc%Domain, maskmap=maskmap)
  else
     call mpp_define_domains((/1,xsz,1,ysz/),layout,Ocean_sfc%Domain)
  endif
  call mpp_get_compute_domain(Ocean_sfc%Domain, isc, iec, jsc, jec)

  allocate ( Ocean_sfc%t_surf (isc:iec,jsc:jec), &
             Ocean_sfc%s_surf (isc:iec,jsc:jec), &
             Ocean_sfc%u_surf (isc:iec,jsc:jec), &
             Ocean_sfc%v_surf (isc:iec,jsc:jec), &
             Ocean_sfc%sea_lev(isc:iec,jsc:jec), &
             Ocean_sfc%area   (isc:iec,jsc:jec), &
             Ocean_sfc%frazil (isc:iec,jsc:jec))

  Ocean_sfc%t_surf  = 0.0  ! time averaged sst (Kelvin) passed to atmosphere/ice model
  Ocean_sfc%s_surf  = 0.0  ! time averaged sss (psu) passed to atmosphere/ice models
  Ocean_sfc%u_surf  = 0.0  ! time averaged u-current (m/sec) passed to atmosphere/ice models
  Ocean_sfc%v_surf  = 0.0  ! time averaged v-current (m/sec)  passed to atmosphere/ice models
  Ocean_sfc%sea_lev = 0.0  ! time averaged thickness of top model grid cell (m) plus patm/rho0/grav
  Ocean_sfc%frazil  = 0.0  ! time accumulated frazil (J/m^2) passed to ice model
  Ocean_sfc%area    = 0.0
  Ocean_sfc%axes    = diag%axesT1%handles !diag axes to be used by coupler tracer flux diagnostics

  if (present(gas_fields_ocn)) then
    call coupler_type_spawn(gas_fields_ocn, Ocean_sfc%fields, (/isc,isc,iec,iec/), &
                              (/jsc,jsc,jec,jec/), suffix = '_ocn', as_needed=.true.)
  endif

end subroutine initialize_ocean_public_type

!> Translates the coupler's ocean_data_type into MOM6's surface state variable.
!! This may eventually be folded into the MOM6's code that calculates the
!! surface state in the first place.
subroutine convert_state_to_ocean_type(state, Ocean_sfc, G, use_conT_absS, &
                                       patm, press_to_z)
  type(surface),           intent(inout) :: state
  type(ocean_public_type), target, intent(inout) :: Ocean_sfc !< Ocean surface state
  type(ocean_grid_type),   intent(inout) :: G                 !< The ocean's grid structure
  logical,                 intent(in)    :: use_conT_absS     !< If true, , the prognostics
                                         !! T&S are the conservative temperature
  real,          optional, intent(in)    :: patm(:,:)         !< Atmospheric pressure.
  real,          optional, intent(in)    :: press_to_z        !< Factor to tranform atmospheric
                                         !! pressure to z?

  ! local variables
  real :: IgR0
  character(len=48)  :: val_str
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd
  integer :: i, j, i0, j0, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  call pass_vector(state%u,state%v,G%Domain)

  call mpp_get_compute_domain(Ocean_sfc%Domain, isc_bnd, iec_bnd, &
                              jsc_bnd, jec_bnd)
  if (present(patm)) then
    ! Check that the inidicies in patm are (isc_bnd:iec_bnd,jsc_bnd:jec_bnd).
    if (.not.present(press_to_z)) call MOM_error(FATAL, &
        'convert_state_to_ocean_type: press_to_z must be present if patm is.')
  endif

  i0 = is - isc_bnd ; j0 = js - jsc_bnd
  !If directed convert the surface T&S
  !from conservative T to potential T and
  !from absolute (reference) salinity to practical salinity
  !
  if(use_conT_absS) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%s_surf(i,j) = gsw_sp_from_sr(state%SSS(i+i0,j+j0))
      Ocean_sfc%t_surf(i,j) = gsw_pt_from_ct(state%SSS(i+i0,j+j0),state%SST(i+i0,j+j0)) + CELSIUS_KELVIN_OFFSET
    enddo ; enddo
  else
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%t_surf(i,j) = state%SST(i+i0,j+j0) + CELSIUS_KELVIN_OFFSET
      Ocean_sfc%s_surf(i,j) = state%SSS(i+i0,j+j0)
    enddo ; enddo
  endif

  do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
    Ocean_sfc%sea_lev(i,j) = state%sea_lev(i+i0,j+j0)
    if (present(patm)) &
      Ocean_sfc%sea_lev(i,j) = Ocean_sfc%sea_lev(i,j) + patm(i,j) * press_to_z
      if (associated(state%frazil)) &
      Ocean_sfc%frazil(i,j) = state%frazil(i+i0,j+j0)
    Ocean_sfc%area(i,j)   =  G%areaT(i+i0,j+j0)
  enddo ; enddo

  if (Ocean_sfc%stagger == AGRID) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%u_surf(i,j) = G%mask2dT(i+i0,j+j0)*0.5*(state%u(I+i0,j+j0)+state%u(I-1+i0,j+j0))
      Ocean_sfc%v_surf(i,j) = G%mask2dT(i+i0,j+j0)*0.5*(state%v(i+i0,J+j0)+state%v(i+i0,J-1+j0))
    enddo ; enddo
  elseif (Ocean_sfc%stagger == BGRID_NE) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%u_surf(i,j) = G%mask2dBu(I+i0,J+j0)*0.5*(state%u(I+i0,j+j0)+state%u(I+i0,j+j0+1))
      Ocean_sfc%v_surf(i,j) = G%mask2dBu(I+i0,J+j0)*0.5*(state%v(i+i0,J+j0)+state%v(i+i0+1,J+j0))
    enddo ; enddo
  elseif (Ocean_sfc%stagger == CGRID_NE) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%u_surf(i,j) = G%mask2dCu(I+i0,j+j0)*state%u(I+i0,j+j0)
      Ocean_sfc%v_surf(i,j) = G%mask2dCv(i+i0,J+j0)*state%v(i+i0,J+j0)
    enddo ; enddo
  else
    write(val_str, '(I8)') Ocean_sfc%stagger
    call MOM_error(FATAL, "convert_state_to_ocean_type: "//&
      "Ocean_sfc%stagger has the unrecognized value of "//trim(val_str))
  endif

  if (coupler_type_initialized(state%tr_fields)) then
    if (.not.coupler_type_initialized(Ocean_sfc%fields)) then
      call MOM_error(FATAL, "convert_state_to_ocean_type: "//&
               "Ocean_sfc%fields has not been initialized.")
    endif
    call coupler_type_copy_data(state%tr_fields, Ocean_sfc%fields)
  endif

end subroutine convert_state_to_ocean_type

!> Returns pointers to objects within ocean_state_type
subroutine get_state_pointers(OS, grid, surf)
  type(ocean_state_type),          pointer :: OS !< Ocean state type
  type(ocean_grid_type), optional, pointer :: grid !< Ocean grid
  type(surface), optional, pointer         :: surf !< Ocean surface state

  if (present(grid)) grid => OS%grid
  if (present(surf)) surf=> OS%state

end subroutine get_state_pointers

!> Maps outgoing ocean data to MCT buffer.
!! See \ref section_ocn_export for a summary of the data
!! that is transferred from MOM6 to MCT.
subroutine ocn_export(ind, ocn_public, grid, o2x)
  type(cpl_indices),       intent(inout) :: ind        !< Structure with coupler
                                                       !! indices and vectors
  type(ocean_public_type), intent(in)    :: ocn_public !< Ocean surface state
  type(ocean_grid_type),   intent(in)    :: grid       !< Ocean model grid
  real(kind=8),            intent(inout) :: o2x(:,:)   !< MCT outgoing bugger
  ! Local variables
  real, dimension(grid%isd:grid%ied,grid%jsd:grid%jed) :: ssh !< Local copy of sea_lev with updated halo
  integer :: i, j, n, ig, jg  !< Grid indices
  real :: slp_L, slp_R, slp_C, slope, u_min, u_max

  ! Copy from ocn_public to o2x. ocn_public uses global indexing with no halos.
  ! The mask comes from "grid" that uses the usual MOM domain that has halos
  ! and does not use global indexing.
  n = 0
  do j=grid%jsc, grid%jec
    jg = j + grid%jdg_offset
    do i=grid%isc,grid%iec
      n = n+1
      ig = i + grid%idg_offset
      ! surface temperature in Kelvin
      o2x(ind%o2x_So_t, n) = ocn_public%t_surf(ig,jg)  * grid%mask2dT(i,j)
      o2x(ind%o2x_So_s, n) = ocn_public%s_surf(ig,jg) * grid%mask2dT(i,j)
      o2x(ind%o2x_So_u, n) = ocn_public%u_surf(ig,jg) * grid%mask2dT(i,j)
      o2x(ind%o2x_So_v, n) = ocn_public%v_surf(ig,jg) * grid%mask2dT(i,j)
      ! Make a copy of ssh in order to do a halo update. We use the usual MOM domain
      ! in order to update halos. i.e. does not use global indexing.
      ssh(i,j) = ocn_public%sea_lev(ig,jg)
    end do
  end do

  ! Update halo of ssh so we can calculate gradients
  call pass_var(ssh, grid%domain)

  ! d/dx ssh
  n = 0
  do j=grid%jsc, grid%jec ; do i=grid%isc,grid%iec
    n = n+1
    ! This is a simple second-order difference
    !o2x(ind%o2x_So_dhdx, n) = 0.5 * (ssh(i+1,j) - ssh(i-1,j)) * grid%IdxT(i,j) * grid%mask2dT(i,j)
    ! This is a PLM slope which might be less prone to the A-grid null mode
    slp_L = (ssh(I,j) - ssh(I-1,j)) * grid%mask2dCu(I-1,j)
    if (grid%mask2dCu(I-1,j)==0.) slp_L = 0.
    slp_R = (ssh(I+1,j) - ssh(I,j)) * grid%mask2dCu(I,j)
    if (grid%mask2dCu(I+1,j)==0.) slp_R = 0.
    slp_C = 0.5 * (slp_L + slp_R)
    if ( (slp_L * slp_R) > 0.0 ) then
      ! This limits the slope so that the edge values are bounded by the
      ! two cell averages spanning the edge.
      u_min = min( ssh(i-1,j), ssh(i,j), ssh(i+1,j) )
      u_max = max( ssh(i-1,j), ssh(i,j), ssh(i+1,j) )
      slope = sign( min( abs(slp_C), 2.*min( ssh(i,j) - u_min, u_max - ssh(i,j) ) ), slp_C )
    else
      ! Extrema in the mean values require a PCM reconstruction avoid generating
      ! larger extreme values.
      slope = 0.0
    end if
    o2x(ind%o2x_So_dhdx, n) = slope * grid%IdxT(i,j) * grid%mask2dT(i,j)
    if (grid%mask2dT(i,j)==0.) o2x(ind%o2x_So_dhdx, n) = 0.0
  end do; end do

  ! d/dy ssh
  n = 0
  do j=grid%jsc, grid%jec ; do i=grid%isc,grid%iec
    n = n+1
    ! This is a simple second-order difference
    !o2x(ind%o2x_So_dhdy, n) = 0.5 * (ssh(i,j+1) - ssh(i,j-1)) * grid%IdyT(i,j) * grid%mask2dT(i,j)
    ! This is a PLM slope which might be less prone to the A-grid null mode
    slp_L = ssh(i,J) - ssh(i,J-1) * grid%mask2dCv(i,J-1)
    if (grid%mask2dCv(i,J-1)==0.) slp_L = 0.
    slp_R = ssh(i,J+1) - ssh(i,J) * grid%mask2dCv(i,J)
    if (grid%mask2dCv(i,J+1)==0.) slp_R = 0.
    slp_C = 0.5 * (slp_L + slp_R)
    !write(6,*)'slp_L, slp_R,i,j,slp_L*slp_R', slp_L, slp_R,i,j,slp_L*slp_R
    if ((slp_L * slp_R) > 0.0) then
      ! This limits the slope so that the edge values are bounded by the
      ! two cell averages spanning the edge.
      u_min = min( ssh(i,j-1), ssh(i,j), ssh(i,j+1) )
      u_max = max( ssh(i,j-1), ssh(i,j), ssh(i,j+1) )
      slope = sign( min( abs(slp_C), 2.*min( ssh(i,j) - u_min, u_max - ssh(i,j) ) ), slp_C )
    else
      ! Extrema in the mean values require a PCM reconstruction avoid generating
      ! larger extreme values.
      slope = 0.0
    end if
    o2x(ind%o2x_So_dhdy, n) = slope * grid%IdyT(i,j) * grid%mask2dT(i,j)
    if (grid%mask2dT(i,j)==0.) o2x(ind%o2x_So_dhdy, n) = 0.0
  end do; end do

end subroutine ocn_export

!> Step forward ocean model for coupling interval
subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o)
  type(ESMF_Clock), intent(inout) :: EClock  !< Time and time step ? \todo Why must this be intent(inout)?
  type(seq_cdata),  intent(inout) :: cdata_o !< Input parameters
  type(mct_aVect),  intent(inout) :: x2o_o   !< Fluxes from coupler to ocean, computed by ocean
  type(mct_aVect),  intent(inout) :: o2x_o   !< Fluxes from ocean to coupler, computed by ocean
  ! Local variables
  type(ESMF_time) :: current_time             !< Time to be reached at the end of ocean cpl interval
  type(ESMF_time) :: time_start_ESMF          !< Time at the start of the coupling interval
  type(ESMF_timeInterval) :: ocn_cpl_interval !< The length of one ocean coupling interval
  integer :: year, month, day, hour, minute, seconds, seconds_n, seconds_d, rc
  logical :: write_restart_at_eod      !< Controls if restart files must be written
  logical :: debug=.false.
  type(time_type) :: time_start        !< Start of coupled time interval to pass to MOM6
  type(time_type) :: coupling_timestep !< Coupled time interval to pass to MOM6
  character(len=128) :: err_msg        !< Error message
  character(len=32)  :: timestamp      !< Name of intermediate restart file
  character(len=384) :: restartname    !< The restart file name (no dir)
  character(len=384) :: restart_pointer_file !< File name for restart pointer file
  character(len=384) :: runid                !< Run ID
  integer            :: nu                   !< i/o unit to write pointer file
  integer            :: shrlogunit ! original log file unit
  integer            :: shrloglev  ! original log level

  ! reset shr logging to ocn log file:
  if (is_root_pe()) then
    call shr_file_getLogUnit(shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit(glb%stdout)
  endif

  ! Compute the time at the start of this coupling interval
  call ESMF_ClockGet(EClock, PrevTime=time_start_ESMF, rc=rc)
  call ESMF_TimeGet(time_start_ESMF, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
  time_start = set_date(year, month, day, hour, minute, seconds, err_msg=err_msg)

  ! Debugging clocks
  if (debug .and. is_root_pe()) then
    call ESMF_ClockGet(EClock, CurrTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_run_mct, current time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, StartTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_run_mct, start time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, StopTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_run_mct, stop time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, PrevTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_run_mct, previous time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, TimeStep=ocn_cpl_interval, rc=rc)
    call ESMF_TimeIntervalGet(ocn_cpl_interval, yy=year, mm=month, d=day, s=seconds, sn=seconds_n, sd=seconds_d, rc=rc)
    write(glb%stdout,*) 'ocn_init_mct, time step: y,m,d-',year,month,day,'s,sn,sd=',seconds,seconds_n,seconds_d
  endif

  ! Translate the coupling time interval
  call ESMF_ClockGet(EClock, TimeStep=ocn_cpl_interval, rc=rc)
  call ESMF_TimeIntervalGet(ocn_cpl_interval, yy=year, mm=month, d=day, s=seconds, sn=seconds_n, sd=seconds_d, rc=rc)
  coupling_timestep = set_time(seconds, days=day, err_msg=err_msg)

  ! set (actually, get from mct) the cdata pointers:
  ! \todo this was done in _init_, is it needed again. Does this infodata need to be in glb%?
  ! GMM, check if  this is needed!
  call seq_cdata_setptrs(cdata_o, infodata=glb%infodata)

  call update_ocean_model(glb%ocn_state, glb%ocn_public, time_start, coupling_timestep, &
                         x2o_o%rattr, glb%ind, glb%sw_decomp, glb%c1, glb%c2, glb%c3, glb%c4)

  ! return export state to driver
  call ocn_export(glb%ind, glb%ocn_public, glb%grid, o2x_o%rattr)

  !--- write out intermediate restart file when needed.
  ! Check alarms for flag to write restart at end of day
  write_restart_at_eod = seq_timemgr_RestartAlarmIsOn(EClock)
  if (debug .and. is_root_pe()) write(glb%stdout,*) 'ocn_run_mct, write_restart_at_eod=', write_restart_at_eod

  if (write_restart_at_eod) then
    ! case name
    call seq_infodata_GetData( glb%infodata, case_name=runid )
    ! add time stamp to the restart filename
    call ESMF_ClockGet(EClock, CurrTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    seconds = seconds + hour*3600 + minute*60
    write(restartname,'(A,".mom6.r.",I4.4,"-",I2.2,"-",I2.2,"-",I5.5)') trim(runid), year, month, day, seconds

    call save_restart(glb%ocn_state%dirs%restart_output_dir, glb%ocn_state%Time, glb%grid, &
                     glb%ocn_state%MOM_CSp%restart_CSp, .false., filename=restartname,GV=glb%ocn_state%GV)

    ! write name of restart file in the rpointer file
    nu = shr_file_getUnit()
    if (is_root_pe()) then
      restart_pointer_file = trim(glb%pointer_filename)
      open(nu, file=restart_pointer_file, form='formatted', status='unknown')
      write(nu,'(a)') trim(restartname) //'.nc'
      close(nu)
      write(glb%stdout,*) 'ocn restart pointer file written: ',trim(restartname)
    endif
    call shr_file_freeUnit(nu)

    ! Is this needed?
    call forcing_save_restart(glb%ocn_state%forcing_CSp, glb%grid, glb%ocn_state%Time, &
                              glb%ocn_state%dirs%restart_output_dir, .true.)
    ! Once we start using the ice shelf module, the following will be needed
    if (glb%ocn_state%use_ice_shelf) then
      call ice_shelf_save_restart(glb%ocn_state%Ice_shelf_CSp, glb%ocn_state%Time, glb%ocn_state%dirs%restart_output_dir, .true.)
    endif

  endif

  ! reset shr logging to original values
  if (is_root_pe()) then
    call shr_file_setLogUnit(shrlogunit)
    call shr_file_setLogLevel(shrloglev)
  endif

end subroutine ocn_run_mct

!> Saves restart fields associated with the forcing
subroutine forcing_save_restart(CS, G, Time, directory, time_stamped, &
                                filename_suffix)
  type(surface_forcing_CS),   pointer       :: CS   !< pointer to the control structure
                                                    !! returned by a previous call to
                                                    !! surface_forcing_init
  type(ocean_grid_type),      intent(inout) :: G    !< The ocean's grid structure
  type(time_type),            intent(in)    :: Time !<  model time at this call
  character(len=*),           intent(in)    :: directory !< optional directory into which
                                                         !! to write these restart files
  logical,          optional, intent(in)    :: time_stamped !< If true, the restart file
                                                            !! names include a unique time
                                                            !! stamp
  character(len=*), optional, intent(in)    :: filename_suffix !< optional suffix
                                               !! (e.g., a time-stamp) to append to the
                                               !! restart file names
  if (.not.associated(CS)) return
  if (.not.associated(CS%restart_CSp)) return
  call save_restart(directory, Time, G, CS%restart_CSp, time_stamped)

end subroutine forcing_save_restart

!> Updates the ocean model fields.  This code wraps the call to step_MOM with MOM6's call.
!! It uses the forcing to advance the ocean model's state from the
!! input value of Ocean_state (which must be for time time_start_update) for a time interval
!! of Ocean_coupling_time_step, returning the publicly visible ocean surface properties in
!! Ocean_sfc and storing the new ocean properties in Ocean_state.
subroutine update_ocean_model(OS, Ocean_sfc, time_start_update, &
                              Ocean_coupling_time_step, x2o_o, ind, sw_decomp,      &
                              c1, c2, c3, c4)
  type(ocean_state_type),        pointer       :: OS !< Structure containing the internal ocean state
  type(ocean_public_type),       intent(inout) :: Ocean_sfc !< Structure containing all the publicly
                                               !! visible ocean surface fields after a coupling time step
  type(time_type), intent(in)                  :: time_start_update !< Time at the beginning of the update step
  type(time_type), intent(in)                  :: Ocean_coupling_time_step !< Amount of time over which to
                                               !! advance the ocean
  real(kind=8),           intent(in)    :: x2o_o(:,:) !< Fluxes from coupler to ocean, computed by ocean
  type(cpl_indices),      intent(inout) :: ind !< Structure with MCT attribute vectors and indices
  logical,                intent(in)    :: sw_decomp !< controls if shortwave is
                                                     !!decomposed into four components
  real(kind=8),           intent(in), optional :: c1, c2, c3, c4 !< Coeffs. used in the shortwave decomposition

  ! local variables
  type(time_type) :: Master_time !< This allows step_MOM to temporarily change
                                 !! the time that is seen by internal modules.
  type(time_type) :: Time1       !< The value of the ocean model's time at the
                                 !! start of a call to step_MOM.
  real :: weight                 !< Flux accumulation weight
  real :: time_step              !< The time step of a call to step_MOM in seconds.
  integer :: secs, days
  integer :: is, ie, js, je

  call callTree_enter("update_ocean_model(), ocn_comp_mct.F90")
  call get_time(Ocean_coupling_time_step, secs, days)
  time_step = 86400.0*real(days) + real(secs)

  if (time_start_update /= OS%Time) then
    call MOM_error(WARNING, "update_ocean_model: internal clock does not "//&
                            "agree with time_start_update argument.")
  endif

  if (.not.associated(OS)) then
    call MOM_error(FATAL, "update_ocean_model called with an unassociated "// &
                    "ocean_state_type structure. ocean_model_init must be "//  &
                    "called first to allocate this structure.")
    return
  endif

  ! This is benign but not necessary if ocean_model_init_sfc was called or if
  ! OS%state%tr_fields was spawnded in ocean_model_init.  Consider removing it.
  is = OS%grid%isc ; ie = OS%grid%iec ; js = OS%grid%jsc ; je = OS%grid%jec
  call coupler_type_spawn(Ocean_sfc%fields, OS%state%tr_fields, &
                          (/is,is,ie,ie/), (/js,js,je,je/), as_needed=.true.)

  weight = 1.0

  if (OS%fluxes%fluxes_used) then
    ! GMM, is enable_averaging needed now?
    call enable_averaging(time_step, OS%Time + Ocean_coupling_time_step, OS%MOM_CSp%diag)
    call ocn_import(OS%fluxes, OS%Time, OS%grid, OS%forcing_CSp, OS%state, x2o_o, ind, sw_decomp, &
                    c1, c2, c3, c4, OS%restore_salinity,OS%restore_temp)
#ifdef _USE_GENERIC_TRACER
    call MOM_generic_tracer_fluxes_accumulate(OS%fluxes, weight) !here weight=1, just saving the current fluxes
#endif

    ! Add ice shelf fluxes
    if (OS%use_ice_shelf) then
      call shelf_calc_flux(OS%State, OS%fluxes, OS%Time, time_step, OS%Ice_shelf_CSp)
    endif

    ! GMM, check ocean_model_MOM.F90 to enable the following option
    !if (OS%icebergs_apply_rigid_boundary)  then
      !This assumes that the iceshelf and ocean are on the same grid. I hope this is true
    !  call add_berg_flux_to_shelf(OS%grid, OS%fluxes,OS%use_ice_shelf,OS%density_iceberg,OS%kv_iceberg, OS%latent_heat_fusion, OS%State, time_step, OS%berg_area_threshold)
    !endif

    ! Indicate that there are new unused fluxes.
    OS%fluxes%fluxes_used = .false.
    OS%fluxes%dt_buoy_accum = time_step
  else
    OS%flux_tmp%C_p = OS%fluxes%C_p
    ! Import fluxes from coupler to ocean. Also, perform do SST and SSS restoring, if needed.
    call ocn_import(OS%flux_tmp, OS%Time, OS%grid, OS%forcing_CSp, &
                       OS%state, x2o_o, ind, sw_decomp, c1, c2, c3, c4, &
                       OS%restore_salinity,OS%restore_temp)

    if (OS%use_ice_shelf) then
      call shelf_calc_flux(OS%State, OS%flux_tmp, OS%Time, time_step, OS%Ice_shelf_CSp)
    endif

    ! GMM, check ocean_model_MOM.F90 to enable the following option
    !if (OS%icebergs_apply_rigid_boundary)  then
     !This assumes that the iceshelf and ocean are on the same grid. I hope this is true
    ! call add_berg_flux_to_shelf(OS%grid, OS%flux_tmp, OS%use_ice_shelf,OS%density_iceberg,OS%kv_iceberg, OS%latent_heat_fusion, OS%State, time_step, OS%berg_area_threshold)
    !endif

    ! Accumulate the forcing over time steps
    call forcing_accumulate(OS%flux_tmp, OS%fluxes, time_step, OS%grid, weight)
#ifdef _USE_GENERIC_TRACER
    call MOM_generic_tracer_fluxes_accumulate(OS%flux_tmp, weight) !weight of the current flux in the running average
#endif
  endif

  if (OS%nstep==0) then
    call finish_MOM_initialization(OS%Time, OS%dirs, OS%MOM_CSp, OS%fluxes)

    call write_energy(OS%MOM_CSp%u, OS%MOM_CSp%v, OS%MOM_CSp%h, OS%MOM_CSp%tv, &
                      OS%Time, 0, OS%grid, OS%GV, OS%sum_output_CSp, &
                      OS%MOM_CSp%tracer_flow_CSp)
  endif

  call disable_averaging(OS%MOM_CSp%diag)
  Master_time = OS%Time ; Time1 = OS%Time

  if(OS%MOM_Csp%offline_tracer_mode) then
    call step_offline(OS%fluxes, OS%state, Time1, time_step, OS%MOM_CSp)
  else
    call step_MOM(OS%fluxes, OS%state, Time1, time_step, OS%MOM_CSp)
  endif

  OS%Time = Master_time + Ocean_coupling_time_step
  OS%nstep = OS%nstep + 1

  call enable_averaging(time_step, OS%Time, OS%MOM_CSp%diag)
  call mech_forcing_diags(OS%fluxes, time_step, OS%grid, &
                          OS%MOM_CSp%diag, OS%forcing_CSp%handles)
  call disable_averaging(OS%MOM_CSp%diag)

  if (OS%fluxes%fluxes_used) then
    call enable_averaging(OS%fluxes%dt_buoy_accum, OS%Time, OS%MOM_CSp%diag)
    call forcing_diagnostics(OS%fluxes, OS%state, OS%fluxes%dt_buoy_accum, &
                             OS%grid, OS%MOM_CSp%diag, OS%forcing_CSp%handles)
    call accumulate_net_input(OS%fluxes, OS%state, OS%fluxes%dt_buoy_accum, &
                              OS%grid, OS%sum_output_CSp)
    call disable_averaging(OS%MOM_CSp%diag)
  endif

!  See if it is time to write out the energy.
  if ((OS%Time + ((Ocean_coupling_time_step)/2) > OS%write_energy_time) .and. &
      (OS%MOM_CSp%t_dyn_rel_adv==0.0)) then
    call write_energy(OS%MOM_CSp%u, OS%MOM_CSp%v, OS%MOM_CSp%h, OS%MOM_CSp%tv, &
                      OS%Time, OS%nstep, OS%grid, OS%GV, OS%sum_output_CSp, &
                      OS%MOM_CSp%tracer_flow_CSp)
    OS%write_energy_time = OS%write_energy_time + OS%energysavedays
  endif

! Translate state into Ocean.
!  call convert_state_to_ocean_type(OS%state, Ocean_sfc, OS%grid, &
!                                   Ice_ocean_boundary%p, OS%press_to_z)
  call convert_state_to_ocean_type(OS%state, Ocean_sfc, OS%grid, &
                                   OS%MOM_CSp%use_conT_absS)

  call callTree_leave("update_ocean_model()")
end subroutine update_ocean_model

!> This function has a few purposes: 1) it allocates and initializes the data
!! in the fluxes structure; 2) it imports surface fluxes using data from
!! the coupler; and 3) it can apply restoring in SST and SSS.
!! See \ref section_ocn_import for a summary of the surface fluxes that are
!! passed from MCT to MOM6, including fluxes that need to be included in
!! the future.
subroutine ocn_import(fluxes, Time, G, CS, state, x2o_o, ind, sw_decomp, &
                             c1, c2, c3, c4, restore_salt, restore_temp)
  type(forcing),              intent(inout)        :: fluxes !< Surface fluxes structure
  type(time_type),            intent(in)           :: Time !< Model time structure
  type(ocean_grid_type),      intent(inout)        :: G  !< The ocean's grid structure
  type(surface_forcing_CS),   pointer              :: CS !< control structure returned by
                                                   !! a previous call to surface_forcing_init
  type(surface),              intent(in)           :: state !< control structure to ocean
                                                   !! surface state fields.
  real(kind=8),               intent(in)           :: x2o_o(:,:)!< Fluxes from coupler to ocean, computed by ocean
  type(cpl_indices),          intent(inout)        :: ind !< Structure with MCT attribute vectors and indices
  logical,                    intent(in)           :: sw_decomp !< controls if shortwave is
                                                   !!decomposed into four components
  real(kind=8),               intent(in), optional :: c1, c2, c3, c4 !< Coeffs. used in the shortwave decomposition
  logical, optional,          intent(in)            :: restore_salt, restore_temp !< Controls if salt and temp are
                                                   !! restored

  ! local variables
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    taux_at_q, &     ! Zonal wind stresses at q points (Pa)
    tauy_at_q        ! Meridional wind stresses at q points (Pa)

  real, dimension(SZI_(G),SZJ_(G)) :: &
    taux_at_h,     & ! Zonal wind stresses at h points (Pa)
    tauy_at_h,     & ! Meridional wind stresses at h points (Pa)
    data_restore,  & ! The surface value toward which to restore (g/kg or degC)
    SST_anom,      & ! Instantaneous sea surface temperature anomalies from a target value (deg C)
    SSS_anom,      & ! Instantaneous sea surface salinity anomalies from a target value (g/kg)
    SSS_mean,      & ! A (mean?) salinity about which to normalize local salinity
                     ! anomalies when calculating restorative precipitation anomalies (g/kg)
    PmE_adj,       & ! The adjustment to PminusE that will cause the salinity
                     ! to be restored toward its target value (kg/(m^2 * s))
    net_FW,        & ! The area integrated net freshwater flux into the ocean (kg/s)
    net_FW2,       & ! The area integrated net freshwater flux into the ocean (kg/s)
    work_sum,      & ! A 2-d array that is used as the work space for a global
                     ! sum, used with units of m2 or (kg/s)
    open_ocn_mask    ! a binary field indicating where ice is present based on frazil criteria

  real :: gustiness     ! unresolved gustiness that contributes to ustar (Pa)
  real :: Irho0         ! inverse of the mean density in (m^3/kg)
  real :: taux2, tauy2  ! squared wind stresses (Pa^2)
  real :: tau_mag       ! magnitude of the wind stress (Pa)
  real :: I_GEarth      ! 1.0 / G%G_Earth  (s^2/m)
  real :: Kv_rho_ice    ! (CS%kv_sea_ice / CS%density_sea_ice) ( m^5/(s*kg) )
  real :: mass_ice      ! mass of sea ice at a face (kg/m^2)
  real :: mass_eff      ! effective mass of sea ice for rigidity (kg/m^2)

  integer :: wind_stagger  ! AGRID, BGRID_NE, or CGRID_NE (integers from MOM_domains)
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, i0, j0
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isr, ier, jsr, jer
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd

  logical :: restore_salinity ! local copy of the argument restore_salt, if it
                              ! is present, or false (no restoring) otherwise.
  logical :: restore_sst      ! local copy of the argument restore_temp, if it
                              ! is present, or false (no restoring) otherwise.
  real :: delta_sss           ! temporary storage for sss diff from restoring value
  real :: delta_sst           ! temporary storage for sst diff from restoring value

  real :: C_p                 ! heat capacity of seawater ( J/(K kg) )

  call cpu_clock_begin(id_clock_forcing)

  is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
  Isq  = G%IscB  ; Ieq  = G%IecB   ; Jsq  = G%JscB  ; Jeq  = G%JecB
  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB
  isr = is-isd+1 ; ier  = ie-isd+1 ; jsr = js-jsd+1 ; jer = je-jsd+1

  C_p                    = fluxes%C_p
  Irho0                  = 1.0/CS%Rho0
  open_ocn_mask(:,:)     = 1.0
  pme_adj(:,:)           = 0.0
  fluxes%vPrecGlobalAdj  = 0.0
  fluxes%vPrecGlobalScl  = 0.0
  fluxes%saltFluxGlobalAdj = 0.0
  fluxes%saltFluxGlobalScl = 0.0
  fluxes%netFWGlobalAdj = 0.0
  fluxes%netFWGlobalScl = 0.0

  restore_salinity = .false.
  if (present(restore_salt)) restore_salinity = restore_salt
  restore_sst = .false.
  if (present(restore_temp)) restore_sst = restore_temp

  ! if true, allocation and initialization
  if (fluxes%dt_buoy_accum < 0) then
    call allocate_forcing_type(G, fluxes, stress=.true., ustar=.true., &
                               water=.true., heat=.true.)

    call safe_alloc_ptr(fluxes%sw_vis_dir,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_vis_dif,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_nir_dir,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_nir_dif,isd,ied,jsd,jed)

    call safe_alloc_ptr(fluxes%p_surf,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%p_surf_full,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%p_surf_SSH,isd,ied,jsd,jed)

    call safe_alloc_ptr(fluxes%salt_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%salt_flux_in,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%salt_flux_added,isd,ied,jsd,jed)

    call safe_alloc_ptr(fluxes%TKE_tidal,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%ustar_tidal,isd,ied,jsd,jed)

    if (CS%allow_flux_adjustments) then
      call safe_alloc_ptr(fluxes%heat_added,isd,ied,jsd,jed)
      call safe_alloc_ptr(fluxes%salt_flux_added,isd,ied,jsd,jed)
    endif

    do j=js-2,je+2 ; do i=is-2,ie+2
      fluxes%TKE_tidal(i,j)   = CS%TKE_tidal(i,j)
      fluxes%ustar_tidal(i,j) = CS%ustar_tidal(i,j)
    enddo; enddo

    if (CS%rigid_sea_ice) then
      call safe_alloc_ptr(fluxes%rigidity_ice_u,IsdB,IedB,jsd,jed)
      call safe_alloc_ptr(fluxes%rigidity_ice_v,isd,ied,JsdB,JedB)
    endif

    if (restore_temp) call safe_alloc_ptr(fluxes%heat_added,isd,ied,jsd,jed)

    fluxes%dt_buoy_accum = 0.0
  endif   ! endif for allocation and initialization

  if (CS%allow_flux_adjustments) then
   fluxes%heat_added(:,:)=0.0
   fluxes%salt_flux_added(:,:)=0.0
  endif

  if (CS%area_surf < 0.0) then
    do j=js,je ; do i=is,ie
      work_sum(i,j) = G%areaT(i,j) * G%mask2dT(i,j)
    enddo ; enddo
    CS%area_surf = reproducing_sum(work_sum, isr, ier, jsr, jer)
  endif    ! endif for allocation and initialization

  do j=js,je ; do i=is,ie
    fluxes%salt_flux(i,j) = 0.0
    fluxes%vprec(i,j) = 0.0
  enddo ; enddo

  ! Salinity restoring logic
  if (restore_salinity) then
    call time_interp_external(CS%id_srestore,Time,data_restore)
    ! open_ocn_mask indicates where to restore salinity (1 means restore, 0 does not)
    open_ocn_mask(:,:) = 1.0
    if (CS%mask_srestore_under_ice) then ! Do not restore under sea-ice
      do j=js,je ; do i=is,ie
        if (state%SST(i,j) .le. -0.0539*state%SSS(i,j)) open_ocn_mask(i,j)=0.0
      enddo; enddo
    endif
    if (CS%salt_restore_as_sflux) then
      do j=js,je ; do i=is,ie
        delta_sss = data_restore(i,j)- state%SSS(i,j)
        delta_sss = sign(1.0,delta_sss)*min(abs(delta_sss),CS%max_delta_srestore)
        fluxes%salt_flux(i,j) = 1.e-3*G%mask2dT(i,j) * (CS%Rho0*CS%Flux_const)* &
                  (CS%basin_mask(i,j)*open_ocn_mask(i,j)) *delta_sss  ! kg Salt m-2 s-1
      enddo; enddo
      if (CS%adjust_net_srestore_to_zero) then
        if (CS%adjust_net_srestore_by_scaling) then
          call adjust_area_mean_to_zero(fluxes%salt_flux, G, fluxes%saltFluxGlobalScl)
          fluxes%saltFluxGlobalAdj = 0.
        else
          work_sum(is:ie,js:je) = G%areaT(is:ie,js:je)*fluxes%salt_flux(is:ie,js:je)
          fluxes%saltFluxGlobalAdj = reproducing_sum(work_sum(:,:), isr,ier, jsr,jer)/CS%area_surf
          fluxes%salt_flux(is:ie,js:je) = fluxes%salt_flux(is:ie,js:je) - fluxes%saltFluxGlobalAdj
        endif
      endif
      fluxes%salt_flux_added(is:ie,js:je) = fluxes%salt_flux(is:ie,js:je) ! Diagnostic
    else
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0.5) then
          delta_sss = state%SSS(i,j) - data_restore(i,j)
          delta_sss = sign(1.0,delta_sss)*min(abs(delta_sss),CS%max_delta_srestore)
          fluxes%vprec(i,j) = (CS%basin_mask(i,j)*open_ocn_mask(i,j))* &
                      (CS%Rho0*CS%Flux_const) * &
                      delta_sss / (0.5*(state%SSS(i,j) + data_restore(i,j)))
        endif
      enddo; enddo
      if (CS%adjust_net_srestore_to_zero) then
        if (CS%adjust_net_srestore_by_scaling) then
          call adjust_area_mean_to_zero(fluxes%vprec, G, fluxes%vPrecGlobalScl)
          fluxes%vPrecGlobalAdj = 0.
        else
          work_sum(is:ie,js:je) = G%areaT(is:ie,js:je)*fluxes%vprec(is:ie,js:je)
          fluxes%vPrecGlobalAdj = reproducing_sum(work_sum(:,:), isr, ier, jsr, jer) / CS%area_surf
          do j=js,je ; do i=is,ie
            fluxes%vprec(i,j) = ( fluxes%vprec(i,j) - fluxes%vPrecGlobalAdj ) * G%mask2dT(i,j)
          enddo; enddo
        endif
      endif
    endif
  endif

  ! SST restoring logic
  if (restore_sst) then
    call time_interp_external(CS%id_trestore,Time,data_restore)
    do j=js,je ; do i=is,ie
       delta_sst = data_restore(i,j)- state%SST(i,j)
       delta_sst = sign(1.0,delta_sst)*min(abs(delta_sst),CS%max_delta_trestore)
       fluxes%heat_added(i,j) = G%mask2dT(i,j) * (CS%Rho0*fluxes%C_p) * delta_sst * CS%Flux_const   ! W m-2
    enddo; enddo
  endif

  ! GMM, CIME uses AGRID. All the BGRID_NE code can be cleaned later
  wind_stagger = AGRID

  if (wind_stagger == BGRID_NE) then
    ! This is necessary to fill in the halo points.
    taux_at_q(:,:) = 0.0 ; tauy_at_q(:,:) = 0.0
  endif
  if (wind_stagger == AGRID) then
    ! This is necessary to fill in the halo points.
    taux_at_h(:,:) = 0.0 ; tauy_at_h(:,:) = 0.0
  endif

  k = 0
  do j=js,je ; do i=is,ie
    k = k + 1 ! Increment position within gindex

    if (wind_stagger == BGRID_NE) then
      taux_at_q(I,J) = x2o_o(ind%x2o_Foxx_taux,k) * CS%wind_stress_multiplier
      tauy_at_q(I,J) = x2o_o(ind%x2o_Foxx_tauy,k) * CS%wind_stress_multiplier
    ! GMM, cime uses AGRID
    elseif (wind_stagger == AGRID) then
      taux_at_h(i,j) = x2o_o(ind%x2o_Foxx_taux,k) * CS%wind_stress_multiplier
      tauy_at_h(i,j) = x2o_o(ind%x2o_Foxx_tauy,k) * CS%wind_stress_multiplier
    else ! C-grid wind stresses.
      fluxes%taux(I,j) = x2o_o(ind%x2o_Foxx_taux,k) * CS%wind_stress_multiplier
      fluxes%tauy(i,J) = x2o_o(ind%x2o_Foxx_tauy,k) * CS%wind_stress_multiplier
    endif

    ! liquid precipitation (rain)
    if (ASSOCIATED(fluxes%lprec)) &
      fluxes%lprec(i,j) = x2o_o(ind%x2o_Faxa_rain,k) * G%mask2dT(i,j)

    ! frozen precipitation (snow)
    if (ASSOCIATED(fluxes%fprec)) &
      fluxes%fprec(i,j) = x2o_o(ind%x2o_Faxa_snow,k) * G%mask2dT(i,j)

    ! evaporation
    if (ASSOCIATED(fluxes%evap)) &
      fluxes%evap(i,j) = x2o_o(ind%x2o_Foxx_evap,k) * G%mask2dT(i,j)

    ! river runoff flux
    if (ASSOCIATED(fluxes%lrunoff)) &
      fluxes%lrunoff(i,j) = x2o_o(ind%x2o_Foxx_rofl,k) * G%mask2dT(i,j)

    ! ice runoff flux
    if (ASSOCIATED(fluxes%frunoff)) &
      fluxes%frunoff(i,j) = x2o_o(ind%x2o_Foxx_rofi,k) * G%mask2dT(i,j)

    ! GMM, we don't have an icebergs yet so the following is not needed
    !if (((ASSOCIATED(IOB%ustar_berg) .and. (.not. ASSOCIATED(fluxes%ustar_berg)))   &
    !  .or. (ASSOCIATED(IOB%area_berg) .and. (.not. ASSOCIATED(fluxes%area_berg)))) &
    !  .or. (ASSOCIATED(IOB%mass_berg) .and. (.not. ASSOCIATED(fluxes%mass_berg)))) &
    !  call allocate_forcing_type(G, fluxes, iceberg=.true.)
    !if (ASSOCIATED(IOB%ustar_berg)) &
    !  fluxes%ustar_berg(i,j) = IOB%ustar_berg(i-i0,j-j0) * G%mask2dT(i,j)
    !if (ASSOCIATED(IOB%area_berg)) &
    !  fluxes%area_berg(i,j) = IOB%area_berg(i-i0,j-j0) * G%mask2dT(i,j)
    !if (ASSOCIATED(IOB%mass_berg)) &
    !  fluxes%mass_berg(i,j) = IOB%mass_berg(i-i0,j-j0) * G%mask2dT(i,j)

    ! GMM, cime does not not have an equivalent for heat_content_lrunoff and
    ! heat_content_frunoff. I am seeting these to zero for now.
    if (ASSOCIATED(fluxes%heat_content_lrunoff)) &
      fluxes%heat_content_lrunoff(i,j) = 0.0 * G%mask2dT(i,j)

    if (ASSOCIATED(fluxes%heat_content_frunoff)) &
      fluxes%heat_content_frunoff(i,j) = 0.0 * G%mask2dT(i,j)

    ! longwave radiation, sum up and down (W/m2)
    if (ASSOCIATED(fluxes%LW)) &
      fluxes%LW(i,j) = (x2o_o(ind%x2o_Faxa_lwdn,k) + x2o_o(ind%x2o_Foxx_lwup,k)) * G%mask2dT(i,j)

    ! sensible heat flux (W/m2)
    if (ASSOCIATED(fluxes%sens)) &
      fluxes%sens(i,j) =  x2o_o(ind%x2o_Foxx_sen,k) * G%mask2dT(i,j)

    ! latent heat flux (W/m^2)
    if (ASSOCIATED(fluxes%latent)) &
      fluxes%latent(i,j) = x2o_o(ind%x2o_Foxx_lat,k) * G%mask2dT(i,j)

    if (sw_decomp) then
      ! Use runtime coefficients to decompose net short-wave heat flux into 4 components
      ! 1) visible, direct shortwave (W/m2)
      if (ASSOCIATED(fluxes%sw_vis_dir)) &
        fluxes%sw_vis_dir(i,j) = G%mask2dT(i,j) * x2o_o(ind%x2o_Foxx_swnet,k)*c1
      ! 2) visible, diffuse shortwave (W/m2)
      if (ASSOCIATED(fluxes%sw_vis_dif)) &
        fluxes%sw_vis_dif(i,j) = G%mask2dT(i,j) * x2o_o(ind%x2o_Foxx_swnet,k)*c2
      ! 3) near-IR, direct shortwave (W/m2)
      if (ASSOCIATED(fluxes%sw_nir_dir)) &
        fluxes%sw_nir_dir(i,j) = G%mask2dT(i,j) * x2o_o(ind%x2o_Foxx_swnet,k)*c3
      ! 4) near-IR, diffuse shortwave (W/m2)
      if (ASSOCIATED(fluxes%sw_nir_dif)) &
        fluxes%sw_nir_dif(i,j) = G%mask2dT(i,j) * x2o_o(ind%x2o_Foxx_swnet,k)*c4

      fluxes%sw(i,j) = fluxes%sw_vis_dir(i,j) + fluxes%sw_vis_dif(i,j) + &
                     fluxes%sw_nir_dir(i,j) + fluxes%sw_nir_dif(i,j)
    else
        call MOM_error(FATAL,"fill_ice_ocean_bnd: this option has not been implemented yet."// &
                       "Shortwave must be decomposed using coeffs. c1, c2, c3, c4.");
    endif

    ! applied surface pressure from atmosphere and cryosphere
    ! sea-level pressure (Pa)
    if (ASSOCIATED(fluxes%p_surf_full) .and. ASSOCIATED(fluxes%p_surf)) then
      fluxes%p_surf_full(i,j) = G%mask2dT(i,j) * x2o_o(ind%x2o_Sa_pslv,k)
      if (CS%max_p_surf >= 0.0) then
          fluxes%p_surf(i,j) = MIN(fluxes%p_surf_full(i,j),CS%max_p_surf)
      else
          fluxes%p_surf(i,j) = fluxes%p_surf_full(i,j)
      endif

      if (CS%use_limited_P_SSH) then
        fluxes%p_surf_SSH(i,j) = fluxes%p_surf(i,j)
      else
        fluxes%p_surf_SSH(i,j) = fluxes%p_surf_full(i,j)
      endif

    endif

    ! salt flux
    ! more salt restoring logic
    if (ASSOCIATED(fluxes%salt_flux)) &
      fluxes%salt_flux(i,j) = G%mask2dT(i,j)*(x2o_o(ind%x2o_Fioi_salt,k) + fluxes%salt_flux(i,j))

    if (ASSOCIATED(fluxes%salt_flux_in)) &
      fluxes%salt_flux_in(i,j) = G%mask2dT(i,j)*x2o_o(ind%x2o_Fioi_salt,k)

  enddo ; enddo
  ! ############################ END OF MCT to MOM ##############################

  ! adjust the NET fresh-water flux to zero, if flagged
  if (CS%adjust_net_fresh_water_to_zero) then
    do j=js,je ; do i=is,ie
      net_FW(i,j) = (((fluxes%lprec(i,j)   + fluxes%fprec(i,j)) + &
                      (fluxes%lrunoff(i,j) + fluxes%frunoff(i,j))) + &
                      (fluxes%evap(i,j)    + fluxes%vprec(i,j)) ) * G%areaT(i,j)
      !   The following contribution appears to be calculating the volume flux of sea-ice
      ! melt. This calculation is clearly WRONG if either sea-ice has variable
      ! salinity or the sea-ice is completely fresh.
      !   Bob thinks this is trying ensure the net fresh-water of the ocean + sea-ice system
      ! is constant.
      !   To do this correctly we will need a sea-ice melt field added to IOB. -AJA
      if (ASSOCIATED(fluxes%salt_flux) .and. (CS%ice_salt_concentration>0.0)) &
        net_FW(i,j) = net_FW(i,j) - G%areaT(i,j) * &
                     (fluxes%salt_flux(i,j) / CS%ice_salt_concentration)
      net_FW2(i,j) = net_FW(i,j)
    enddo ; enddo

    if (CS%adjust_net_fresh_water_by_scaling) then
      call adjust_area_mean_to_zero(net_FW2, G, fluxes%netFWGlobalScl)
      do j=js,je ; do i=is,ie
        fluxes%vprec(i,j) = fluxes%vprec(i,j) + (net_FW2(i,j) - net_FW(i,j)) * G%mask2dT(i,j)
      enddo; enddo
    else
      fluxes%netFWGlobalAdj = reproducing_sum(net_FW(:,:), isr, ier, jsr, jer) / CS%area_surf
      do j=js,je ; do i=is,ie
        fluxes%vprec(i,j) = ( fluxes%vprec(i,j) - fluxes%netFWGlobalAdj ) * G%mask2dT(i,j)
      enddo; enddo
    endif

  endif

  ! surface momentum stress related fields as function of staggering
  if (wind_stagger == BGRID_NE) then
    if (G%symmetric) &
      call fill_symmetric_edges(taux_at_q, tauy_at_q, G%Domain, stagger=BGRID_NE)
    call pass_vector(taux_at_q, tauy_at_q, G%Domain, stagger=BGRID_NE)

    do j=js,je ; do I=Isq,Ieq
      fluxes%taux(I,j) = 0.0
      if ((G%mask2dBu(I,J) + G%mask2dBu(I,J-1)) > 0) &
        fluxes%taux(I,j) = (G%mask2dBu(I,J)*taux_at_q(I,J) + &
                            G%mask2dBu(I,J-1)*taux_at_q(I,J-1)) / &
                           (G%mask2dBu(I,J) + G%mask2dBu(I,J-1))
    enddo ; enddo

    do J=Jsq,Jeq ; do i=is,ie
      fluxes%tauy(i,J) = 0.0
      if ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J)) > 0) &
        fluxes%tauy(i,J) = (G%mask2dBu(I,J)*tauy_at_q(I,J) + &
                            G%mask2dBu(I-1,J)*tauy_at_q(I-1,J)) / &
                           (G%mask2dBu(I,J) + G%mask2dBu(I-1,J))
    enddo ; enddo

    ! ustar is required for the bulk mixed layer formulation. The background value
    ! of 0.02 Pa is a relatively small value intended to give reasonable behavior
    ! in regions of very weak winds.

    do j=js,je ; do i=is,ie
      tau_mag = 0.0 ; gustiness = CS%gust_const
      if (((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + &
           (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) > 0) then
        tau_mag = sqrt(((G%mask2dBu(I,J)*(taux_at_q(I,J)**2 + tauy_at_q(I,J)**2) + &
            G%mask2dBu(I-1,J-1)*(taux_at_q(I-1,J-1)**2 + tauy_at_q(I-1,J-1)**2)) + &
           (G%mask2dBu(I,J-1)*(taux_at_q(I,J-1)**2 + tauy_at_q(I,J-1)**2) + &
            G%mask2dBu(I-1,J)*(taux_at_q(I-1,J)**2 + tauy_at_q(I-1,J)**2)) ) / &
          ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) )
        if (CS%read_gust_2d) gustiness = CS%gust(i,j)
      endif
      fluxes%ustar(i,j) = sqrt(gustiness*Irho0 + Irho0*tau_mag)
    enddo ; enddo

  elseif (wind_stagger == AGRID) then
    call pass_vector(taux_at_h, tauy_at_h, G%Domain,stagger=AGRID)

    do j=js,je ; do I=Isq,Ieq
      fluxes%taux(I,j) = 0.0
      if ((G%mask2dT(i,j) + G%mask2dT(i+1,j)) > 0) &
        fluxes%taux(I,j) = (G%mask2dT(i,j)*taux_at_h(i,j) + &
                            G%mask2dT(i+1,j)*taux_at_h(i+1,j)) / &
                           (G%mask2dT(i,j) + G%mask2dT(i+1,j))
    enddo ; enddo

    do J=Jsq,Jeq ; do i=is,ie
      fluxes%tauy(i,J) = 0.0
      if ((G%mask2dT(i,j) + G%mask2dT(i,j+1)) > 0) &
        fluxes%tauy(i,J) = (G%mask2dT(i,j)*tauy_at_h(i,j) + &
                            G%mask2dT(i,J+1)*tauy_at_h(i,j+1)) / &
                           (G%mask2dT(i,j) + G%mask2dT(i,j+1))
    enddo ; enddo

    do j=js,je ; do i=is,ie
      gustiness = CS%gust_const
      if (CS%read_gust_2d .and. (G%mask2dT(i,j) > 0)) gustiness = CS%gust(i,j)
      fluxes%ustar(i,j) = sqrt(gustiness*Irho0 + Irho0 * G%mask2dT(i,j) * &
                               sqrt(taux_at_h(i,j)**2 + tauy_at_h(i,j)**2))
    enddo ; enddo

  else ! C-grid wind stresses.
    if (G%symmetric) &
      call fill_symmetric_edges(fluxes%taux, fluxes%tauy, G%Domain)
    call pass_vector(fluxes%taux, fluxes%tauy, G%Domain)

    do j=js,je ; do i=is,ie
      taux2 = 0.0
      if ((G%mask2dCu(I-1,j) + G%mask2dCu(I,j)) > 0) &
        taux2 = (G%mask2dCu(I-1,j)*fluxes%taux(I-1,j)**2 + &
                 G%mask2dCu(I,j)*fluxes%taux(I,j)**2) / (G%mask2dCu(I-1,j) + G%mask2dCu(I,j))

      tauy2 = 0.0
      if ((G%mask2dCv(i,J-1) + G%mask2dCv(i,J)) > 0) &
        tauy2 = (G%mask2dCv(i,J-1)*fluxes%tauy(i,J-1)**2 + &
                 G%mask2dCv(i,J)*fluxes%tauy(i,J)**2) / (G%mask2dCv(i,J-1) + G%mask2dCv(i,J))

      if (CS%read_gust_2d) then
        fluxes%ustar(i,j) = sqrt(CS%gust(i,j)*Irho0 + Irho0*sqrt(taux2 + tauy2))
      else
        fluxes%ustar(i,j) = sqrt(CS%gust_const*Irho0 + Irho0*sqrt(taux2 + tauy2))
      endif
    enddo ; enddo

  endif   ! endif for wind related fields


  ! sea ice related fields
  if (CS%rigid_sea_ice) then
    ! The commented out code here and in the following lines is the correct
    ! version, but the incorrect version is being retained temporarily to avoid
    ! changing answers.
    call pass_var(fluxes%p_surf_full, G%Domain)
    I_GEarth = 1.0 / G%G_Earth
    Kv_rho_ice = (CS%kv_sea_ice / CS%density_sea_ice)
    do I=isd,ied-1 ; do j=jsd,jed
      mass_ice = min(fluxes%p_surf_full(i,j), fluxes%p_surf_full(i+1,j)) * I_GEarth
      mass_eff = 0.0
      if (mass_ice > CS%rigid_sea_ice_mass) then
        mass_eff = (mass_ice - CS%rigid_sea_ice_mass) **2 / &
                   (mass_ice + CS%rigid_sea_ice_mass)
      endif
      ! CAUTION: with both rigid_sea_ice and ice shelves, we will need to make this
      ! a maximum for the second call.
      fluxes%rigidity_ice_u(I,j) = Kv_rho_ice * mass_eff
    enddo ; enddo
    do i=isd,ied ; do J=jsd,jed-1
      mass_ice = min(fluxes%p_surf_full(i,j), fluxes%p_surf_full(i,j+1)) * I_GEarth
      mass_eff = 0.0
      if (mass_ice > CS%rigid_sea_ice_mass) then
        mass_eff = (mass_ice - CS%rigid_sea_ice_mass) **2 / &
                   (mass_ice + CS%rigid_sea_ice_mass)
      endif
      fluxes%rigidity_ice_v(i,J) = Kv_rho_ice * mass_eff
    enddo ; enddo
  endif

  if (CS%allow_flux_adjustments) then
    ! Apply adjustments to fluxes
    call apply_flux_adjustments(G, CS, Time, fluxes)
  endif

  ! Allow for user-written code to alter fluxes after all the above
  call user_alter_forcing(state, fluxes, Time, G, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)

end subroutine  ocn_import

!> Adds flux adjustments obtained via data_override
!! Component name is 'OCN'
!! Available adjustments are:
!! - taux_adj (Zonal wind stress delta, positive to the east, in Pa)
!! - tauy_adj (Meridional wind stress delta, positive to the north, in Pa)
subroutine apply_flux_adjustments(G, CS, Time, fluxes)
  type(ocean_grid_type),    intent(inout) :: G !< Ocean grid structure
  type(surface_forcing_CS), pointer       :: CS !< Surface forcing control structure
  type(time_type),          intent(in)    :: Time !< Model time structure
  type(forcing), optional,  intent(inout) :: fluxes !< Surface fluxes structure

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: tempx_at_h ! Delta to zonal wind stress at h points (Pa)
  real, dimension(SZI_(G),SZJ_(G)) :: tempy_at_h ! Delta to meridional wind stress at h points (Pa)
  real, dimension(SZI_(G),SZJ_(G)) :: temp_at_h ! Fluxes at h points (W m-2 or kg m-2 s-1)

  integer :: isc, iec, jsc, jec, i, j
  real :: dLonDx, dLonDy, rDlon, cosA, sinA, zonal_tau, merid_tau
  logical :: overrode_x, overrode_y, overrode_h

  isc = G%isc; iec = G%iec
  jsc = G%jsc; jec = G%jec

  overrode_h = .false.
  call data_override('OCN', 'hflx_adj', temp_at_h(isc:iec,jsc:jec), Time, override=overrode_h)

  if (overrode_h) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      fluxes%heat_added(i,j) = fluxes%heat_added(i,j) + temp_at_h(i,j)* G%mask2dT(i,j)
    enddo; enddo
  endif

  call pass_var(fluxes%heat_added, G%Domain)

  overrode_h = .false.
  call data_override('OCN', 'sflx_adj', temp_at_h(isc:iec,jsc:jec), Time, override=overrode_h)

  if (overrode_h) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      fluxes%salt_flux_added(i,j) = fluxes%salt_flux_added(i,j) + temp_at_h(i,j)* G%mask2dT(i,j)
    enddo; enddo
  endif

  call pass_var(fluxes%salt_flux_added, G%Domain)
  overrode_h = .false.

  call data_override('OCN', 'prcme_adj', temp_at_h(isc:iec,jsc:jec), Time, override=overrode_h)

  if (overrode_h) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      fluxes%vprec(i,j) = fluxes%vprec(i,j) + temp_at_h(i,j)* G%mask2dT(i,j)
    enddo; enddo
  endif

  call pass_var(fluxes%vprec, G%Domain)


  tempx_at_h(:,:) = 0.0 ; tempy_at_h(:,:) = 0.0
  ! Either reads data or leaves contents unchanged
  overrode_x = .false. ; overrode_y = .false.
  call data_override('OCN', 'taux_adj', tempx_at_h(isc:iec,jsc:jec), Time, override=overrode_x)
  call data_override('OCN', 'tauy_adj', tempy_at_h(isc:iec,jsc:jec), Time, override=overrode_y)

  if (overrode_x .or. overrode_y) then
    if (.not. (overrode_x .and. overrode_y)) call MOM_error(FATAL,"apply_flux_adjustments: "//&
            "Both taux_adj and tauy_adj must be specified, or neither, in data_table")

    ! Rotate winds
    call pass_vector(tempx_at_h, tempy_at_h, G%Domain, To_All, AGRID)
    do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
      dLonDx = G%geoLonCu(I,j) - G%geoLonCu(I-1,j)
      dLonDy = G%geoLonCv(i,J) - G%geoLonCv(i,J-1)
      rDlon = sqrt( dLonDx * dLonDx + dLonDy * dLonDy )
      if (rDlon > 0.) rDlon = 1. / rDlon
      cosA = dLonDx * rDlon
      sinA = dLonDy * rDlon
      zonal_tau = tempx_at_h(i,j)
      merid_tau = tempy_at_h(i,j)
      tempx_at_h(i,j) = cosA * zonal_tau - sinA * merid_tau
      tempy_at_h(i,j) = sinA * zonal_tau + cosA * merid_tau
    enddo ; enddo

    ! Average to C-grid locations
    do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
      fluxes%taux(I,j) = fluxes%taux(I,j) + 0.5 * ( tempx_at_h(i,j) + tempx_at_h(i+1,j) )
    enddo ; enddo

    do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
      fluxes%tauy(i,J) = fluxes%tauy(i,J) + 0.5 * ( tempy_at_h(i,j) + tempy_at_h(i,j+1) )
    enddo ; enddo
  endif ! overrode_x .or. overrode_y

end subroutine apply_flux_adjustments

!> Finalizes MOM6
!!
!! \todo This needs to be done here.
subroutine ocn_final_mct( EClock, cdata_o, x2o_o, o2x_o)
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o  !< Fluxes from coupler to ocean, computed by ocean
    type(mct_aVect)             , intent(inout) :: o2x_o  !< Fluxes from ocean to coupler, computed by ocean

    call ocean_model_end(glb%ocn_public, glb%ocn_state, glb%ocn_state%Time)

end subroutine ocn_final_mct

!> Terminates the model run, saving the ocean state in a
!! restart file and deallocating any data associated with the ocean.
subroutine ocean_model_end(Ocean_sfc, Ocean_state, Time)
  type(ocean_public_type),     intent(inout) :: Ocean_sfc !< An ocean_public_type structure that is to be
                                                          !! deallocated upon termination.
  type(ocean_state_type),      pointer       :: Ocean_state!<  pointer to the structure containing the internal
                                 !                        !! ocean state to be deallocated upon termination.
  type(time_type),             intent(in)    :: Time      !< The model time, used for writing restarts.

  !if (debug .and. is_root_pe()) write(glb%stdout,*)'Here 1'
  !GMM call save_restart(Ocean_state, Time)
  call diag_mediator_end(Time, Ocean_state%MOM_CSp%diag)
  call MOM_end(Ocean_state%MOM_CSp)
  if (Ocean_state%use_ice_shelf) call ice_shelf_end(Ocean_state%Ice_shelf_CSp)
  !if (debug .and. is_root_pe()) write(glb%stdout,*)'Here 2'

end subroutine ocean_model_end

!> Sets mct global segment maps for the MOM decomposition.
!!
!! \todo Find out if we should only provide indirect indexing for ocean points and not land.
subroutine ocn_SetGSMap_mct(mpicom_ocn, MOM_MCT_ID, gsMap_ocn, gsMap3d_ocn)
  integer,         intent(in)    :: mpicom_ocn  !< MPI communicator
  integer,         intent(in)    :: MOM_MCT_ID  !< MCT component ID
  type(mct_gsMap), intent(inout) :: gsMap_ocn   !< MCT global segment map for 2d data
  type(mct_gsMap), intent(inout) :: gsMap3d_ocn !< MCT global segment map for 3d data
  ! Local variables
  integer :: lsize   !< Local size of indirect indexing array
  integer :: i, j, k !< Local indices
  integer :: ni, nj  !< Declared sizes of h-point arrays
  integer :: ig, jg  !< Global indices
  type(ocean_grid_type), pointer :: grid => NULL() !< A pointer to a grid structure
  integer, allocatable           :: gindex(:)      !< Indirect indices

  grid => glb%grid ! for convenience
  if (.not. associated(grid)) call MOM_error(FATAL, 'ocn_comp_mct.F90, ocn_SetGSMap_mct():' // &
      'grid returned from get_state_pointers() was not associated!')

  ! Size of computational domain
  lsize = ( grid%iec - grid%isc + 1 ) * ( grid%jec - grid%jsc + 1 )

  ! Size of global domain
  call get_global_grid_size(grid, ni, nj)

  ! Create indirect indices for the computational domain
  allocate(gindex(lsize))

  ! Set indirect indices in gindex
  k = 0
  do j = grid%jsc, grid%jec
    jg = j + grid%jdg_offset ! TODO: check this calculation
    do i = grid%isc, grid%iec
      ig = i + grid%idg_offset ! TODO: check this calculation
      k = k + 1 ! Increment position within gindex
      gindex(k) = ni * (jg - 1) + ig
    enddo
  enddo

  ! Tell MCT how to indirectly index into the 2d buffer
  call mct_gsMap_init(gsMap_ocn, gindex, mpicom_ocn, MOM_MCT_ID, lsize, ni * nj)

  deallocate(gindex)

end subroutine ocn_SetGSMap_mct

!> Sets MCT global segment maps for the MOM6 decomposition
subroutine ocn_domain_mct( lsize, gsMap_ocn, dom_ocn)
  integer        , intent(in)    :: lsize      !< Size of attr. vector
  type(mct_gsMap), intent(in)    :: gsMap_ocn  !< MCT global segment map for 2d data
  type(mct_ggrid), intent(inout) :: dom_ocn    !< WHAT IS THIS?

  ! Local Variables
  integer, parameter              :: SHR_REAL_R8 = selected_real_kind(12)
  integer, pointer                :: idata(:)
  integer                         :: i,j,k
  real(kind=SHR_REAL_R8), pointer :: data(:)
  real(kind=SHR_REAL_R8)          :: m2_to_rad2
  type(ocean_grid_type), pointer :: grid => NULL() ! A pointer to a grid structure

  grid => glb%grid ! for convenience

  ! set coords to lat and lon, and areas to rad^2
  call mct_gGrid_init(GGrid=dom_ocn, CoordChars=trim(seq_flds_dom_coord), &
                      OtherChars=trim(seq_flds_dom_other), lsize=lsize )

  call mct_avect_zero(dom_ocn%data)
  allocate(data(lsize))

  ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
  k = pe_here()
  call mct_gsMap_orderedPoints(gsMap_ocn, k, idata)
  call mct_gGrid_importIAttr(dom_ocn,'GlobGridNum',idata,lsize)

  !initialization
  data(:) = -9999.0
  call mct_gGrid_importRAttr(dom_ocn,"lat"  ,data,lsize)
  call mct_gGrid_importRAttr(dom_ocn,"lon"  ,data,lsize)
  call mct_gGrid_importRAttr(dom_ocn,"area" ,data,lsize)
  call mct_gGrid_importRAttr(dom_ocn,"aream",data,lsize)
  data(:) = 0.0
  call mct_gGrid_importRAttr(dom_ocn,"mask",data,lsize)
  call mct_gGrid_importRAttr(dom_ocn,"frac",data,lsize)

  k = 0
  do j = grid%jsc, grid%jec
    do i = grid%isc, grid%iec
      k = k + 1 ! Increment position within gindex
      data(k) = grid%geoLonT(i,j)
    enddo
  enddo
  call mct_gGrid_importRattr(dom_ocn,"lon",data,lsize)

  k = 0
  do j = grid%jsc, grid%jec
    do i = grid%isc, grid%iec
      k = k + 1 ! Increment position within gindex
      data(k) = grid%geoLatT(i,j)
    enddo
  enddo
  call mct_gGrid_importRattr(dom_ocn,"lat",data,lsize)

  k = 0
  m2_to_rad2 = 1./grid%Rad_Earth**2
  do j = grid%jsc, grid%jec
    do i = grid%isc, grid%iec
      k = k + 1 ! Increment position within gindex
      data(k) = grid%AreaT(i,j) * m2_to_rad2
    enddo
  enddo
  call mct_gGrid_importRattr(dom_ocn,"area",data,lsize)

  k = 0
  do j = grid%jsc, grid%jec
    do i = grid%isc, grid%iec
      k = k + 1 ! Increment position within gindex
      data(k) = grid%mask2dT(i,j)
    enddo
  enddo
  call mct_gGrid_importRattr(dom_ocn,"mask",data,lsize)
  call mct_gGrid_importRattr(dom_ocn,"frac",data,lsize)

  deallocate(data)
  deallocate(idata)

end subroutine ocn_domain_mct

!> \namespace ocn_comp_mct
!!
!! \section section_ocn_import Fluxes imported from the coupler (MCT) to MOM6
!! The following summarizes the mismatches between MCT and MOM6 in terms
!! of ice ocean fluxes.
!!
!! Redundancies:
!! x2o_Faxa_prec = x2o_Faxa_rain + x2o_Faxa_snow
!!
!! Variables whose units and sign  **could not** be verified so far:
!! x2o_Foxx_rofl
!! x2o_Foxx_rof
!!
!! Variables in MOM6 fluxes that are **NOT** filled by the coupler:
!! ustar_berg, frictional velocity beneath icebergs (m/s)
!! area_berg, area covered by icebergs(m2/m2)
!! mass_berg, mass of icebergs(kg/m2)
!! runoff_hflx, heat content of liquid runoff (W/m2)
!! calving_hflx, heat content of frozen runoff (W/m2)
!! mi, mass of ice (kg/m2)
!!
!! Variables in the coupler that are **NOT** used in MOM6 (i.e., no corresponding field in fluxes):
!! x2o_Fioi_melth, heat flux from snow & ice melt (W/m2)
!! x2o_Fioi_meltw, snow melt flux (kg/m2/s)
!! x2o_Si_ifrac, fractional ice wrt ocean
!! x2o_So_duu10n, 10m wind speed squared (m^2/s^2)
!! x2o_Sa_co2prog, bottom atm level prognostic CO2
!! x2o_Sa_co2diag, bottom atm level diagnostic CO2
!!
!! \TODO Langmuir related fields:
!! surface Stokes drift, x-comp. (x2o_Sw_ustokes)
!! surface Stokes drift, y-comp. (x2o_Sw_vstokes)
!! wave model langmuir multiplier (x2o_Sw_lamult)
!!
!! \TODO Biogeochemistry:
!! x2o_Fioi_bcpho, Black Carbon hydrophobic release from sea ice component
!! x2o_Fioi_bcphi, Black Carbon hydrophilic release from sea ice component
!! x2o_Fioi_flxdst, Dust release from sea ice component
!! x2o_Faxa_bcphidry, Black Carbon hydrophilic dry deposition
!! x2o_Faxa_bcphodry, Black Carbon hydrophobic dry deposition
!! x2o_Faxa_bcphiwet, Black Carbon hydrophobic wet deposition
!! x2o_Faxa_ocphidry, Organic Carbon hydrophilic dry deposition
!! x2o_Faxa_ocphodry, Organic Carbon hydrophobic dry deposition
!! x2o_Faxa_ocphiwet, Organic Carbon hydrophilic dry deposition
!! x2o_Faxa_dstwet, Sizes 1 to 4 dust - wet deposition
!! x2o_Faxa_dstdry, Sizes 1 to 4 dust - dry deposition
!!
!! \section section_ocn_export Fluxes exported from MOM6 to the coupler (MCT)
!!
!! Variables that are currently being exported:
!!
!! Surface temperature (Kelvin)
!! Surface salinity (psu)
!! Surface eastward velocity (m/s)
!! Surface northward velocity (m/s)
!! Zonal slope in the sea surface height
!! Meridional slope in the sea surface height
!!
!! \TODO Variables that **are not** currently being exported:
!!
!! Boundary layer depth
!! CO2
!! DMS
!! o2x_Fioo_q       !< Heat flux?

end module ocn_comp_mct
