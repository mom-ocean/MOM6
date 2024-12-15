!> Interface for surface waves
module MOM_wave_interface

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_data_override, only : data_override_init, data_override
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : diag_ctrl
use MOM_domains,       only : pass_var, pass_vector, AGRID
use MOM_domains,       only : To_South, To_West, To_All
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_forcing_type,  only : mech_forcing
use MOM_grid,          only : ocean_grid_type
use MOM_hor_index,     only : hor_index_type
use MOM_io,            only : file_exists, get_var_sizes, read_variable
use MOM_io,            only : vardesc, var_desc
use MOM_safe_alloc,    only : safe_alloc_ptr
use MOM_time_manager,  only : time_type, operator(+), operator(/)
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs, surface
use MOM_verticalgrid,  only : verticalGrid_type
use MOM_restart,       only : register_restart_pair, MOM_restart_CS

implicit none ; private

#include <MOM_memory.h>

public MOM_wave_interface_init ! Public interface to fully initialize the wave routines.
public query_wave_properties ! Public interface to obtain information from the waves control structure.
public Update_Surface_Waves ! Public interface to update wave information at the
                            ! coupler/driver level.
public Update_Stokes_Drift ! Public interface to update the Stokes drift profiles
                           ! called in step_mom.
public get_Langmuir_Number ! Public interface to compute Langmuir number called from
                           ! ePBL or KPP routines.
public Stokes_PGF ! Public interface to compute Stokes-shear induced pressure gradient force anomaly
public StokesMixing ! NOT READY - Public interface to add down-Stokes gradient
                    ! momentum mixing (e.g. the approach of Harcourt 2013/2015)
public CoriolisStokes ! NOT READY - Public interface to add Coriolis-Stokes acceleration
                      ! of the mean currents, needed for comparison with LES.  It is
                      ! presently advised against implementing in non-1d settings without
                      ! serious consideration of the full 3d wave-averaged Navier-Stokes
                      ! CL2 effects.
public Waves_end ! public interface to deallocate and free wave related memory.
public get_wave_method ! public interface to obtain the wave method string
public waves_register_restarts ! public interface to register wave restart fields

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Container for all surface wave related parameters
type, public :: wave_parameters_CS ; private

  ! Main surface wave options and publicly visible variables
  logical, public :: UseWaves = .false.     !< Flag to enable surface gravity wave feature
  logical, public :: Stokes_VF = .false.    !< True if Stokes vortex force is used
  logical, public :: Passive_Stokes_VF = .false. !< Computes Stokes VF, but doesn't affect dynamics
  logical, public :: Stokes_PGF = .false.   !< True if Stokes shear pressure Gradient force is used
  logical, public :: robust_Stokes_PGF = .false.  !< If true, use expressions to calculate the
                                            !! Stokes-induced pressure gradient anomalies that are
                                            !! more accurate in the limit of thin layers.
  logical, public :: Passive_Stokes_PGF = .false. !< Keeps Stokes_PGF on, but doesn't affect dynamics
  logical, public :: Stokes_DDT = .false.   !< Developmental:
                                            !! True if Stokes d/dt is used
  logical, public :: Passive_Stokes_DDT = .false.   !< Keeps Stokes_DDT on, but doesn't affect dynamics

  real, allocatable, dimension(:,:,:), public :: &
    Us_x               !< 3d zonal Stokes drift profile [L T-1 ~> m s-1]
                       !! Horizontal -> U points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    Us_y               !< 3d meridional Stokes drift profile [L T-1 ~> m s-1]
                       !! Horizontal -> V points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    ddt_Us_x           !< 3d time tendency of zonal Stokes drift profile [L T-2 ~> m s-2]
                       !! Horizontal -> U points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    ddt_Us_y           !< 3d time tendency of meridional Stokes drift profile [L T-2 ~> m s-2]
                       !! Horizontal -> V points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    Us_x_from_ddt      !< Check of 3d zonal Stokes drift profile [L T-1 ~> m s-1]
                       !! Horizontal -> U points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    Us_y_from_ddt      !< Check of 3d meridional Stokes drift profile [L T-1 ~> m s-1]
                       !! Horizontal -> V points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    Us_x_prev          !< 3d zonal Stokes drift profile, previous dynamics call [L T-1 ~> m s-1]
                       !! Horizontal -> U points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    Us_y_prev          !< 3d meridional Stokes drift profile, previous dynamics call [L T-1 ~> m s-1]
                       !! Horizontal -> V points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    KvS                !< Viscosity for Stokes Drift shear [H Z T-1 ~> m2 s-1 or Pa s]
  real, allocatable, dimension(:), public :: &
    WaveNum_Cen        !< Wavenumber bands for read/coupled [Z-1 ~> m-1]
  real, allocatable, dimension(:,:,:), public :: &
    UStk_Hb            !< Surface Stokes Drift spectrum (zonal) [L T-1 ~> m s-1]
                       !! Horizontal -> H-points
                       !! 3rd dimension -> Freq/Wavenumber
  real, allocatable, dimension(:,:,:), public :: &
    VStk_Hb            !< Surface Stokes Drift spectrum (meridional) [L T-1 ~> m s-1]
                       !! Horizontal -> H-points
                       !! 3rd dimension -> Freq/Wavenumber
  real, allocatable, dimension(:,:), public :: &
    Omega_w2x          !< wind direction ccw from model x- axis   [nondim radians]
  integer, public :: NumBands = 0   !< Number of wavenumber/frequency partitions
                                    !! Must match the number of bands provided
                                    !! via either coupling or file.

  ! The remainder of this control structure is private
  integer :: WaveMethod = -99 !< Options for including wave information
                              !! Valid (tested) choices are:
                              !!   0 - Test Profile
                              !!   1 - Surface Stokes Drift Bands
                              !!   2 - DHH85
                              !!   3 - LF17
                              !! -99 - No waves computed, but empirical Langmuir number used.
  logical :: LagrangianMixing !< This feature is in development and not ready
                              !! True if Stokes drift is present and mixing
                              !! should be applied to Lagrangian current
                              !! (mean current + Stokes drift).
                              !! See Reichl et al., 2016 KPP-LT approach
  logical :: StokesMixing     !< This feature is in development and not ready.
                              !! True if vertical mixing of momentum
                              !! should be applied directly to Stokes current
                              !! (with separate mixing parameter for Eulerian
                              !! mixing contribution).
                              !! See Harcourt 2013, 2015 Second-Moment approach
  logical :: CoriolisStokes   !< This feature is in development and not ready.
                              ! True if Coriolis-Stokes acceleration should be applied.
  real :: Stokes_min_thick_avg !< A layer thickness below which the cell-center Stokes drift is
                              !! used instead of the cell average [Z ~> m].  This is only used if
                              !! WAVE_INTERFACE_ANSWER_DATE < 20230101.
  integer :: answer_date      !< The vintage of the order of arithmetic and expressions in the
                              !! surface wave calculations.  Values below 20230101 recover the
                              !! answers from the end of 2022, while higher values use updated
                              !! and more robust forms of the same expressions.

  ! Options if WaveMethod is Surface Stokes Drift Bands (1)
  integer :: PartitionMode  !< Method for partition mode (meant to check input)
                            !! 0 - wavenumbers
                            !! 1 - frequencies
  integer :: DataSource !< Integer that specifies where the model Looks for data
                        !! Valid choices are:
                        !! 1 - FMS DataOverride Routine
                        !! 2 - Reserved For Coupler
                        !! 3 - User input (fixed values, useful for 1d testing)

  ! Options if using FMS DataOverride Routine
  character(len=40)  :: SurfBandFileName !< Filename if using DataOverride
  real :: land_speed    !< A large Stokes velocity that can be used to indicate land values in
                        !! a data override file [L T-1 ~> m s-1].  Stokes drift components larger
                        !! than this are set to zero in data override calls for the Stokes drift.
  logical :: DataOver_initialized !< Flag for DataOverride Initialization

  ! Options for computing Langmuir number
  real :: LA_FracHBL         !< Fraction of OSBL for averaging Langmuir number [nondim]
  real :: LA_HBL_min         !< Minimum boundary layer depth for averaging Langmuir number [Z ~> m]
  logical :: LA_Misalignment = .false. !< Flag to use misalignment in Langmuir number
  logical :: LA_misalign_bug = .false. !< Flag to use code with a sign error when calculating the
                       !! misalignment between the shear and waves in the Langmuir number calculation.
  real :: g_Earth      !< The gravitational acceleration, equivalent to GV%g_Earth but with
                       !! different dimensional rescaling appropriate for deep-water gravity
                       !! waves [Z T-2 ~> m s-2]
  real :: I_g_Earth    !< The inverse of the gravitational acceleration, with dimensional rescaling
                       !! appropriate for deep-water gravity waves [T2 Z-1 ~> s2 m-1]
  ! Surface Wave Dependent 1d/2d/3d vars
  real, allocatable, dimension(:) :: &
    Freq_Cen           !< Central frequency for wave bands, including a factor of 2*pi [T-1 ~> s-1]
  real, allocatable, dimension(:) :: &
    PrescribedSurfStkX !< Surface Stokes drift if prescribed [L T-1 ~> m s-1]
  real, allocatable, dimension(:) :: &
    PrescribedSurfStkY !< Surface Stokes drift if prescribed [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:) :: &
    La_Turb            !< Aligned Turbulent Langmuir number [nondim]
                       !! Horizontal -> H points
  real, allocatable, dimension(:,:) :: &
    US0_x              !< Surface Stokes Drift (zonal) [L T-1 ~> m s-1]
                       !! Horizontal -> U points
  real, allocatable, dimension(:,:) :: &
    US0_y              !< Surface Stokes Drift (meridional) [L T-1 ~> m s-1]
                       !! Horizontal -> V points
  real, allocatable, dimension(:,:,:) :: &
    STKx0              !< Stokes Drift spectrum (zonal) [L T-1 ~> m s-1]
                       !! Horizontal -> U points
                       !! 3rd dimension -> Freq/Wavenumber
  real, allocatable, dimension(:,:,:) :: &
    STKy0              !< Stokes Drift spectrum (meridional) [L T-1 ~> m s-1]
                       !! Horizontal -> V points
                       !! 3rd dimension -> Freq/Wavenumber

  real :: La_min       !< An arbitrary lower-bound on the Langmuir number [nondim].
                       !! Langmuir number is sqrt(u_star/u_stokes).  When both are small
                       !! but u_star is orders of magnitude smaller, the Langmuir number could
                       !! have unintended consequences.  Since both are small it can be safely
                       !! capped to avoid such consequences.
  real :: La_Stk_backgnd !< A small background Stokes velocity used in the denominator of
                       !! some expressions for the Langmuir number [L T-1 ~> m s-1]

  ! Parameters used in estimating the wind speed or wave properties from the friction velocity
  real :: VonKar = -1.0 !< The von Karman coefficient as used in the MOM_wave_interface module [nondim]
  real :: rho_air  !< A typical density of air at sea level, as used in wave calculations [R ~> kg m-3]
  real :: nu_air   !< The viscosity of air, as used in wave calculations [Z2 T-1 ~> m2 s-1]
  real :: rho_ocn  !< A typical surface density of seawater, as used in wave calculations in
                   !! comparison with the density of air [R ~> kg m-3].  The default is RHO_0.
  real :: SWH_from_u10sq !< A factor for converting the square of the 10 m wind speed to the
                   !! significant wave height [Z T2 L-2 ~> s2 m-1]
  real :: Charnock_min !< The minimum value of the Charnock coefficient, which relates the square of
                   !! the air friction velocity divided by the gravitational acceleration to the
                   !! wave roughness length [nondim]
  real :: Charnock_slope_U10 !< The partial derivative of the Charnock coefficient with the 10 m wind
                   !! speed [T L-1 ~> s m-1].   Note that in eq. 13 of the Edson et al. 2013 describing
                   !! the COARE 3.5 bulk flux algorithm, this slope is given as 0.017.  However, 0.0017
                   !! reproduces the curve in their figure 6, so that is the default value used in MOM6.
  real :: Charnock_intercept !< The intercept of the fit for the Charnock coefficient in the limit of
                   !! no wind [nondim].  Note that this can be negative because CHARNOCK_MIN will keep
                   !! the final value for the Charnock coefficient from being from being negative.

  ! Options used with the test profile
  real    :: TP_STKX0     !< Test profile x-stokes drift amplitude [L T-1 ~> m s-1]
  real    :: TP_STKY0     !< Test profile y-stokes drift amplitude [L T-1 ~> m s-1]
  real    :: TP_WVL       !< Test profile wavelength [Z ~> m]

  ! Options for use with the Donelan et al., 1985 (DHH85) spectrum
  logical :: WaveAgePeakFreq !< Flag to use wave age to determine the peak frequency with DHH85
  logical :: StaticWaves  !< Flag to disable updating DHH85 Stokes drift
  logical :: DHH85_is_set !< The if the wave properties have been set when WaveMethod = DHH85.
  real    :: WaveAge      !< The fixed wave age used with the DHH85 spectrum [nondim]
  real    :: WaveWind     !< Wind speed for the DHH85 spectrum [L T-1 ~> m s-1]
  real    :: omega_min    !< Minimum wave frequency with the DHH85 spectrum [T-1 ~> s-1]
  real    :: omega_max    !< Maximum wave frequency with the DHH85 spectrum [T-1 ~> s-1]

  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  !>@{ Diagnostic handles
  integer, public :: id_PFu_Stokes = -1 , id_PFv_Stokes = -1
  integer, public :: id_3dstokes_x_from_ddt = -1 , id_3dstokes_y_from_ddt = -1
  integer :: id_P_deltaStokes_L = -1, id_P_deltaStokes_i = -1
  integer :: id_surfacestokes_x = -1 , id_surfacestokes_y = -1
  integer :: id_3dstokes_x = -1 , id_3dstokes_y = -1
  integer :: id_ddt_3dstokes_x = -1 , id_ddt_3dstokes_y = -1
  integer :: id_La_turb = -1
  !>@}

end type wave_parameters_CS

! Switches needed in import_stokes_drift
!>@{ Enumeration values for the wave method
integer, parameter :: TESTPROF = 0, SURFBANDS = 1, DHH85 = 2, LF17 = 3, EFACTOR = 4, NULL_WaveMethod = -99
!>@}
!>@{ Enumeration values for the wave data source
integer, parameter :: DATAOVR = 1, COUPLER = 2, INPUT = 3
!>@}

! Strings for the wave method
character*(5), parameter  :: NULL_STRING      = "EMPTY"         !< null wave method string
character*(12), parameter :: TESTPROF_STRING  = "TEST_PROFILE"  !< test profile string
character*(13), parameter :: SURFBANDS_STRING = "SURFACE_BANDS" !< surface bands string
character*(5), parameter  :: DHH85_STRING     = "DHH85"         !< DHH85 wave method string
character*(4), parameter  :: LF17_STRING      = "LF17"          !< LF17 wave method string
character*(7), parameter  :: EFACTOR_STRING   = "EFACTOR"       !< EFACTOR (based on vr12-ma) wave method string

contains

!> Initializes parameters related to MOM_wave_interface
subroutine MOM_wave_interface_init(time, G, GV, US, param_file, CS, diag)
  type(time_type), target, intent(in)    :: Time       !< Model time
  type(ocean_grid_type),   intent(inout) :: G          !< Grid structure
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Input parameter structure
  type(wave_parameters_CS), pointer      :: CS         !< Wave parameter control structure
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostic Pointer

  ! Local variables
  character(len=40)  :: mdl = "MOM_wave_interface" !< This module's name.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character*(13) :: TMPSTRING1, TMPSTRING2
  character*(12), parameter :: DATAOVR_STRING   = "DATAOVERRIDE"
  character*(7), parameter  :: COUPLER_STRING   = "COUPLER"
  character*(5), parameter  :: INPUT_STRING     = "INPUT"
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags
  logical :: use_waves
  logical :: StatisticalWaves

  ! Dummy Check
  if (.not. associated(CS)) then
    call MOM_error(FATAL, "wave_interface_init called without an associated control structure.")
    return
  endif

  call get_param(param_file, mdl, "USE_WAVES", use_waves, &
       "If true, enables surface wave modules.", default=.false.)

  ! Check if using LA_LI2016
  call get_param(param_file, mdl, "USE_LA_LI2016", StatisticalWaves, &
                 do_not_log=.true.,default=.false.)

  if (.not.(use_waves .or. StatisticalWaves)) return

  CS%UseWaves = use_waves
  CS%diag => diag
  CS%Time => Time

  CS%g_Earth = GV%g_Earth_Z_T2
  CS%I_g_Earth = 1.0 / CS%g_Earth

  ! Add any initializations needed here
  CS%DataOver_initialized = .false.

  call log_version(param_file, mdl, version)

  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)

  call get_param(param_file, mdl, "WAVE_INTERFACE_ANSWER_DATE", CS%answer_date, &
                 "The vintage of the order of arithmetic and expressions in the surface wave "//&
                 "calculations.  Values below 20230101 recover the answers from the end of 2022, "//&
                 "while higher values use updated and more robust forms of the same expressions:\n"//&
                 "\t <  20230101 - Original answers for wave interface routines\n"//&
                 "\t >= 20230101 - More robust expressions for Update_Stokes_Drift\n"//&
                 "\t >= 20230102 - More robust expressions for get_StokesSL_LiFoxKemper\n"//&
                 "\t >= 20230103 - More robust expressions for ust_2_u10_coare3p5", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) CS%answer_date = max(CS%answer_date, 20230701)

  ! Langmuir number Options
  call get_param(param_file, mdl, "LA_DEPTH_RATIO", CS%LA_FracHBL, &
                 "The depth (normalized by BLD) to average Stokes drift over in "//&
                 "Langmuir number calculation, where La = sqrt(ust/Stokes).", &
                 units="nondim", default=0.04)
  call get_param(param_file, mdl, "LA_DEPTH_MIN", CS%LA_HBL_min, &
                 "The minimum depth over which to average the Stokes drift in the Langmuir "//&
                 "number calculation.", units="m", default=0.1, scale=US%m_to_Z)

  if (StatisticalWaves) then
    CS%WaveMethod = LF17
    call set_LF17_wave_params(param_file, mdl, GV, US, CS)
    if (.not.use_waves) return
  else
    CS%WaveMethod = NULL_WaveMethod
  end if

  ! Wave modified physics
  !  Presently these are all in research mode
  call get_param(param_file, mdl, "LAGRANGIAN_MIXING", CS%LagrangianMixing, &
                 "Flag to use Lagrangian Mixing of momentum", default=.false., &
                 do_not_log=.not.use_waves)
  if (CS%LagrangianMixing) then
    ! Force Code Intervention
    call MOM_error(FATAL,"Should you be enabling Lagrangian Mixing? Code not ready.")
  endif
  call get_param(param_file, mdl, "STOKES_MIXING", CS%StokesMixing, &
                 "Flag to use Stokes Mixing of momentum", default=.false., &
                 do_not_log=.not.use_waves)
  if (CS%StokesMixing) then
    ! Force Code Intervention
    call MOM_error(FATAL, "Should you be enabling Stokes Mixing? Code not ready.")
  endif
  call get_param(param_file, mdl, "CORIOLIS_STOKES", CS%CoriolisStokes, &
                 "Flag to use Coriolis Stokes acceleration", default=.false., &
                 do_not_log=.not.use_waves)
  if (CS%CoriolisStokes) then
    ! Force Code Intervention
    call MOM_error(FATAL, "Should you be enabling Coriolis-Stokes? Code not ready.")
  endif

  call get_param(param_file, mdl, "STOKES_VF", CS%Stokes_VF, &
       "Flag to use Stokes vortex force", &
       default=.false.)
  call get_param(param_file, mdl, "PASSIVE_STOKES_VF", CS%Passive_Stokes_VF, &
       "Flag to make Stokes vortex force diagnostic only.", &
       default=.false.)
  call get_param(param_file, mdl, "STOKES_PGF", CS%Stokes_PGF, &
       "Flag to use Stokes-induced pressure gradient anomaly", &
       default=.false.)
  call get_param(param_file, mdl, "ROBUST_STOKES_PGF", CS%robust_Stokes_PGF, &
       "If true, use expressions to calculate the Stokes-induced pressure gradient "//&
       "anomalies that are more accurate in the limit of thin layers.", &
       default=.false., do_not_log=.not.CS%Stokes_PGF)
       !### Change the default for ROBUST_STOKES_PGF to True.
  call get_param(param_file, mdl, "PASSIVE_STOKES_PGF", CS%Passive_Stokes_PGF, &
       "Flag to make Stokes-induced pressure gradient anomaly diagnostic only.", &
       default=.false.)
  call get_param(param_file, mdl, "STOKES_DDT", CS%Stokes_DDT, &
       "Flag to use Stokes d/dt", &
       default=.false.)
  call get_param(param_file, mdl, "PASSIVE_STOKES_DDT", CS%Passive_Stokes_DDT, &
       "Flag to make Stokes d/dt diagnostic only", &
       default=.false.)

  ! Get Wave Method and write to integer WaveMethod
  call get_param(param_file,mdl,"WAVE_METHOD",TMPSTRING1,             &
       "Choice of wave method, valid options include: \n"//           &
       " TEST_PROFILE  - Prescribed from surface Stokes drift \n"//   &
       "                 and a decay wavelength.\n"//                 &
       " SURFACE_BANDS - Computed from multiple surface values \n"//  &
       "                 and decay wavelengths.\n"//                  &
       " DHH85         - Uses Donelan et al. 1985 empirical \n"//     &
       "                 wave spectrum with prescribed values. \n"//  &
       " LF17          - Infers Stokes drift profile from wind \n"//  &
       "                 speed following Li and Fox-Kemper 2017.\n"// &
       " EFACTOR       - Applies an enhancement factor to the KPP\n"//&
       "                 turbulent velocity scale received \n"//      &
       "                 directly from WW3 and is based on the \n"//  &
       "                 surface layer and projected Langmuir \n"//   &
       "                 number (Li 2016)\n", &
       default=NULL_STRING)
  select case (TRIM(TMPSTRING1))
  case (NULL_STRING)! No Waves
    call MOM_error(FATAL, "wave_interface_init called with no specified "//&
                           "WAVE_METHOD.")
  case (TESTPROF_STRING)! Test Profile
    CS%WaveMethod = TESTPROF
    call get_param(param_file, mdl, "TP_STKX_SURF", CS%TP_STKX0, &
         'Surface Stokes (x) for test profile', &
         units='m/s', default=0.1, scale=US%m_s_to_L_T)
    call get_param(param_file, mdl, "TP_STKY_SURF", CS%TP_STKY0, &
         'Surface Stokes (y) for test profile', &
         units='m/s', default=0.0, scale=US%m_s_to_L_T)
    call get_param(param_file,mdl, "TP_WVL", CS%TP_WVL, &
         'Wavelength for test profile', &
         units='m', default=50.0, scale=US%m_to_Z)
  case (SURFBANDS_STRING)! Surface Stokes Drift Bands
    CS%WaveMethod = SURFBANDS
    call get_param(param_file, mdl, "SURFBAND_MIN_THICK_AVG", CS%Stokes_min_thick_avg, &
                 "A layer thickness below which the cell-center Stokes drift is used instead of "//&
                 "the cell average.  This is only used if WAVE_INTERFACE_ANSWER_DATE < 20230101.", &
                 units="m", default=0.1, scale=US%m_to_Z, do_not_log=(CS%answer_date>=20230101))
    call get_param(param_file, mdl, "SURFBAND_SOURCE", TMPSTRING2, &
                 "Choice of SURFACE_BANDS data mode, valid options include: \n"//&
                 " DATAOVERRIDE  - Read from NetCDF using FMS DataOverride. \n"//&
                 " COUPLER       - Look for variables from coupler pass \n"//&
                 " INPUT         - Testing with fixed values.", default=NULL_STRING)
    select case (TRIM(TMPSTRING2))
    case (NULL_STRING)! Default
      call MOM_error(FATAL, "wave_interface_init called with SURFACE_BANDS"//&
                           " but no SURFBAND_SOURCE.")
    case (DATAOVR_STRING)! Using Data Override
      CS%DataSource = DATAOVR
      call get_param(param_file, mdl, "SURFBAND_FILENAME", CS%SurfBandFileName, &
                 "Filename of surface Stokes drift input band data.", default="StkSpec.nc")
      call get_param(param_file, mdl, "SURFBAND_OVERRIDE_LAND_SPEED", CS%land_speed, &
                 "A large Stokes velocity that can be used to indicate land values in "//&
                 "a data override file.  Stokes drift components larger than this are "//&
                 "set to zero in data override calls for the Stokes drift.", &
                 units="m s-1", default=10.0, scale=US%m_s_to_L_T)
    case (COUPLER_STRING)! Reserved for coupling
      CS%DataSource = COUPLER
      ! This is just to make something work, but it needs to be read from the wavemodel.
      call get_param(param_file, mdl, "STK_BAND_COUPLER",CS%NumBands, &
                 "STK_BAND_COUPLER is the number of Stokes drift bands in the coupler. "//&
                 "This has to be consistent with the number of Stokes drift bands in WW3, "//&
                 "or the model will fail.", default=1)
      allocate( CS%WaveNum_Cen(CS%NumBands), source=0.0 )
      allocate( CS%STKx0(G%isdB:G%iedB,G%jsd:G%jed,CS%NumBands), source=0.0 )
      allocate( CS%STKy0(G%isd:G%ied,G%jsdB:G%jedB,CS%NumBands), source=0.0 )
      allocate( CS%UStk_Hb(G%isc:G%iec,G%jsc:G%jec,CS%NumBands), source=0.0 )
      allocate( CS%VStk_Hb(G%isc:G%iec,G%jsc:G%jec,CS%NumBands), source=0.0 )
      allocate( CS%Omega_w2x(G%isc:G%iec,G%jsc:G%jec)          , source=0.0 )
      CS%PartitionMode = 0
      call get_param(param_file, mdl, "SURFBAND_WAVENUMBERS", CS%WaveNum_Cen, &
           "Central wavenumbers for surface Stokes drift bands.", &
           units='rad/m', default=0.12566, scale=US%Z_to_m)
    case (INPUT_STRING)! A method to input the Stokes band (globally uniform)
      CS%DataSource = INPUT
      call get_param(param_file, mdl, "SURFBAND_NB", CS%NumBands, &
                 "Prescribe number of wavenumber bands for Stokes drift. "//&
                 "Make sure this is consistnet w/ WAVENUMBERS, STOKES_X, and "//&
                 "STOKES_Y, there are no safety checks in the code.", default=1)
      allocate( CS%WaveNum_Cen(1:CS%NumBands), source=0.0 )
      allocate( CS%PrescribedSurfStkX(1:CS%NumBands), source=0.0 )
      allocate( CS%PrescribedSurfStkY(1:CS%NumBands), source=0.0 )
      allocate( CS%STKx0(G%isdB:G%iedB,G%jsd:G%jed,1:CS%NumBands), source=0.0 )
      allocate( CS%STKy0(G%isd:G%ied,G%jsdB:G%jedB,1:CS%NumBands), source=0.0 )

      CS%PartitionMode = 0
      call get_param(param_file, mdl, "SURFBAND_WAVENUMBERS", CS%WaveNum_Cen, &
           "Central wavenumbers for surface Stokes drift bands.", &
           units='rad/m', default=0.12566, scale=US%Z_to_m)
      call get_param(param_file, mdl, "SURFBAND_STOKES_X", CS%PrescribedSurfStkX, &
           "X-direction surface Stokes drift for bands.", &
           units='m/s', default=0.15, scale=US%m_s_to_L_T)
      call get_param(param_file, mdl, "SURFBAND_STOKES_Y", CS%PrescribedSurfStkY, &
           "Y-direction surface Stokes drift for bands.", &
           units='m/s', default=0.0, scale=US%m_s_to_L_T)
    case default! No method provided
      call MOM_error(FATAL,'Check WAVE_METHOD.')
    end select

  case (DHH85_STRING) !Donelan et al., 1985 spectrum
    CS%WaveMethod = DHH85
    call MOM_error(WARNING,"DHH85 only ever set-up for uniform cases w/"//&
                           " Stokes drift in x-direction.")
    call get_param(param_file, mdl, "DHH85_AGE_FP", CS%WaveAgePeakFreq, &
         "Choose true to use waveage in peak frequency.", default=.false.)
    call get_param(param_file, mdl, "DHH85_AGE", CS%WaveAge, &
         "Wave Age for DHH85 spectrum.", &
         units='nondim', default=1.2)
    call get_param(param_file, mdl, "DHH85_WIND", CS%WaveWind, &
         "Wind speed for DHH85 spectrum.", &
         units='m s-1', default=10.0, scale=US%m_s_to_L_T)
    call get_param(param_file, mdl, "DHH85_MIN_WAVE_FREQ", CS%omega_min, &
                 "Minimum wave frequency for the DHH85 spectrum.", &
                 units='s-1', default=0.1, scale=US%T_to_s)
    call get_param(param_file, mdl, "DHH85_MAX_WAVE_FREQ", CS%omega_max, &
                 "Maximum wave frequency for the DHH85 spectrum.", &
                 units='s-1', default=10.0, scale=US%T_to_s) ! The default is about a 30 cm cutoff wavelength.
    call get_param(param_file, mdl, "STATIC_DHH85", CS%StaticWaves, &
         "Flag to disable updating DHH85 Stokes drift.", default=.false.)
  case (LF17_STRING) !Li and Fox-Kemper 17 wind-sea Langmuir number
    CS%WaveMethod = LF17
    call set_LF17_wave_params(param_file, mdl, GV, US, CS)
  case (EFACTOR_STRING) !Li and Fox-Kemper 16
    CS%WaveMethod = EFACTOR
  case default
    call MOM_error(FATAL,'Check WAVE_METHOD.')
  end select

  ! Langmuir number Options  (Note that CS%LA_FracHBL is set above.)
  call get_param(param_file, mdl, "LA_MISALIGNMENT", CS%LA_Misalignment, &
         "Flag (logical) if using misalignment between shear and waves in LA", &
         default=.false.)
  call get_param(param_file, mdl, "LA_MISALIGNMENT_BUG", CS%LA_misalign_bug, &
         "If true, use a code with a sign error when calculating the misalignment between "//&
         "the shear and waves when LA_MISALIGNMENT is true.", &
         default=.false., do_not_log=.not.CS%LA_Misalignment)
  call get_param(param_file, mdl, "MIN_LANGMUIR", CS%La_min,    &
         "A minimum value for all Langmuir numbers that is not physical, "//&
         "but is likely only encountered when the wind is very small and "//&
         "therefore its effects should be mostly benign.", &
         units="nondim", default=0.05)
  call get_param(param_file, mdl, "LANGMUIR_STOKES_BACKGROUND", CS%La_Stk_backgnd, &
         "A small background Stokes velocity used in the denominator of some "//&
         "expressions for the Langmuir number.", &
         units="m s-1", default=1.0e-10, scale=US%m_s_to_L_T, do_not_log=(CS%WaveMethod==LF17))

  ! Allocate and initialize
  ! a. Stokes driftProfiles
  allocate(CS%Us_x(G%isdB:G%IedB,G%jsd:G%jed,G%ke), source=0.0)
  allocate(CS%Us_y(G%isd:G%Ied,G%jsdB:G%jedB,G%ke), source=0.0)
  if (CS%Stokes_DDT) then
    !allocate(CS%Us_x_prev(G%isdB:G%IedB,G%jsd:G%jed,G%ke), source=0.0)
    !allocate(CS%Us_y_prev(G%isd:G%Ied,G%jsdB:G%jedB,G%ke), source=0.0)
    allocate(CS%ddt_Us_x(G%isdB:G%IedB,G%jsd:G%jed,G%ke), source=0.0)
    allocate(CS%ddt_Us_y(G%isd:G%Ied,G%jsdB:G%jedB,G%ke), source=0.0)
    allocate(CS%Us_x_from_ddt(G%isdB:G%IedB,G%jsd:G%jed,G%ke), source=0.0)
    allocate(CS%Us_y_from_ddt(G%isd:G%Ied,G%jsdB:G%jedB,G%ke), source=0.0)
  endif
  ! b. Surface Values
  allocate(CS%US0_x(G%isdB:G%iedB,G%jsd:G%jed), source=0.0)
  allocate(CS%US0_y(G%isd:G%ied,G%jsdB:G%jedB), source=0.0)
  ! c. Langmuir number
  allocate(CS%La_turb(G%isc:G%iec,G%jsc:G%jec), source=0.0)
  ! d. Viscosity for Stokes drift
  if (CS%StokesMixing) then
    allocate(CS%KvS(G%isd:G%Ied,G%jsd:G%jed,GV%ke), source=0.0)
  endif

  ! Initialize Wave related outputs
  CS%id_surfacestokes_y = register_diag_field('ocean_model','surface_stokes_y', &
       CS%diag%axesCv1,Time,'Surface Stokes drift (y)', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_surfacestokes_x = register_diag_field('ocean_model','surface_stokes_x', &
       CS%diag%axesCu1,Time,'Surface Stokes drift (x)', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_3dstokes_y = register_diag_field('ocean_model','3d_stokes_y', &
       CS%diag%axesCvL,Time,'3d Stokes drift (y)', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_3dstokes_x = register_diag_field('ocean_model','3d_stokes_x', &
       CS%diag%axesCuL,Time,'3d Stokes drift (x)', 'm s-1', conversion=US%L_T_to_m_s)
  if (CS%Stokes_DDT) then
    CS%id_ddt_3dstokes_y = register_diag_field('ocean_model','dvdt_Stokes', &
         CS%diag%axesCvL,Time,'d/dt Stokes drift (meridional)', 'm s-2', conversion=US%L_T2_to_m_s2)
    CS%id_ddt_3dstokes_x = register_diag_field('ocean_model','dudt_Stokes', &
         CS%diag%axesCuL,Time,'d/dt Stokes drift (zonal)', 'm s-2', conversion=US%L_T2_to_m_s2)
    CS%id_3dstokes_y_from_ddt = register_diag_field('ocean_model','3d_stokes_y_from_ddt', &
         CS%diag%axesCvL,Time,'3d Stokes drift from ddt (y)', 'm s-1', conversion=US%L_T_to_m_s)
    CS%id_3dstokes_x_from_ddt = register_diag_field('ocean_model','3d_stokes_x_from_ddt', &
         CS%diag%axesCuL,Time,'3d Stokes drift from ddt (x)', 'm s-1', conversion=US%L_T_to_m_s)
  endif
  CS%id_PFv_Stokes = register_diag_field('ocean_model','PFv_Stokes', &
       CS%diag%axesCvL,Time,'PF from Stokes drift (meridional)','m s-2',conversion=US%L_T2_to_m_s2)
  CS%id_PFu_Stokes = register_diag_field('ocean_model','PFu_Stokes', &
       CS%diag%axesCuL,Time,'PF from Stokes drift (zonal)','m s-2',conversion=US%L_T2_to_m_s2)
  CS%id_P_deltaStokes_i = register_diag_field('ocean_model','P_deltaStokes_i', &
       CS%diag%axesTi,Time,'Interfacial pressure anomaly from Stokes drift used in PFu_Stokes',&
       'm2 s-2',conversion=US%L_T_to_m_s**2)
  CS%id_P_deltaStokes_L = register_diag_field('ocean_model','P_deltaStokes_L', &
       CS%diag%axesTL,Time,'Layer averaged pressure anomaly from Stokes drift used in PFu_Stokes',&
       'm2 s-2',conversion=US%L_T_to_m_s**2)
  CS%id_La_turb = register_diag_field('ocean_model','La_turbulent', &
       CS%diag%axesT1,Time,'Surface (turbulent) Langmuir number','nondim')

end subroutine MOM_wave_interface_init

!> Set the parameters that are used to determine the averaged Stokes drift and Langmuir numbers
subroutine set_LF17_wave_params(param_file, mdl, GV, US, CS)
  type(param_file_type),   intent(in)    :: param_file !< Input parameter structure
  character(len=*),        intent(in)    :: mdl        !< A module name to use in the get_param calls
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(wave_parameters_CS), pointer      :: CS         !< Wave parameter control structure

  ! A separate routine is used to set these parameters because there are multiple ways that the
  ! underlying parameterizations are enabled.

  call get_param(param_file, mdl, "VISCOSITY_AIR", CS%nu_air, &
                 "A typical viscosity of air at sea level, as used in wave calculations", &
                 units="m2 s-1", default=1.0e-6, scale=US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "VON_KARMAN_WAVES", CS%vonKar, &
                 "The value the von Karman constant as used for surface wave calculations.", &
                 units="nondim", default=0.40)  ! The default elsewhere in MOM6 is usually 0.41.
  call get_param(param_file, mdl, "RHO_AIR", CS%rho_air, &
                 "A typical density of air at sea level, as used in wave calculations", &
                 units="kg m-3", default=1.225, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "RHO_SFC_WAVES", CS%Rho_ocn, &
                 "A typical surface density of seawater, as used in wave calculations in "//&
                 "comparison with the density of air.  The default is RHO_0.", &
                 units="kg m-3", default=GV%Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "WAVE_HEIGHT_SCALE_FACTOR", CS%SWH_from_u10sq, &
                 "A factor relating the square of the 10 m wind speed to the significant "//&
                 "wave height, with a default value based on the Pierson-Moskowitz spectrum.", &
                 units="s2 m-1", default=0.0246, scale=US%m_to_Z*US%L_T_to_m_s**2)
  call get_param(param_file, mdl, "CHARNOCK_MIN", CS%Charnock_min, &
                 "The minimum value of the Charnock coefficient, which relates the square of "//&
                 "the air friction velocity divided by the gravitational acceleration to the "//&
                 "wave roughness length.", units="nondim", default=0.028)
  call get_param(param_file, mdl, "CHARNOCK_SLOPE_U10", CS%Charnock_slope_U10, &
                 "The partial derivative of the Charnock coefficient with the 10 m wind speed.  "//&
                 "Note that in eq. 13 of the Edson et al. 2013 describing the COARE 3.5 bulk "//&
                 "flux algorithm, this slope is given as 0.017.  However, 0.0017 reproduces "//&
                 "the curve in their figure 6, so that is the default value used in MOM6.", &
                 units="s m-1", default=0.0017, scale=US%L_T_to_m_s)
  call get_param(param_file, mdl, "CHARNOCK_0_WIND_INTERCEPT", CS%Charnock_intercept, &
                 "The intercept of the fit for the Charnock coefficient in the limit of no wind.  "//&
                 "Note that this can be negative because CHARNOCK_MIN will keep the final "//&
                 "value for the Charnock coefficient from being from being negative.", &
                 units="nondim", default=-0.005)

end subroutine set_LF17_wave_params

!> This interface provides the caller with information from the waves control structure.
subroutine query_wave_properties(CS, NumBands, WaveNumbers, US)
  type(wave_parameters_CS),        pointer     :: CS   !< Wave parameter Control structure
  integer,               optional, intent(out) :: NumBands    !< If present, this returns the number of
                                                       !!< wavenumber partitions in the wave discretization
  real, dimension(:),    optional, intent(out) :: Wavenumbers !< If present this returns the characteristic
                                                       !! wavenumbers of the wave discretization [m-1] or [Z-1 ~> m-1]
  type(unit_scale_type), optional, intent(in)  :: US   !< A dimensional unit scaling type that is used to undo
                                                       !! the dimensional scaling of the output variables, if present
  integer :: n

  if (present(NumBands)) NumBands = CS%NumBands
  if (present(Wavenumbers)) then
    if (size(Wavenumbers) < CS%NumBands) call MOM_error(FATAL, "query_wave_properties called "//&
                                "with a Wavenumbers array that is smaller than the number of bands.")
    if (present(US)) then
      do n=1,CS%NumBands ; Wavenumbers(n) = US%m_to_Z * CS%WaveNum_Cen(n) ; enddo
    else
      do n=1,CS%NumBands ; Wavenumbers(n) = CS%WaveNum_Cen(n) ; enddo
    endif
  endif

end subroutine query_wave_properties

!> Subroutine that handles updating of surface wave/Stokes drift related properties
subroutine Update_Surface_Waves(G, GV, US, Time_present, dt, CS, forces)
  type(wave_parameters_CS), pointer    :: CS  !< Wave parameter Control structure
  type(ocean_grid_type), intent(inout) :: G   !< Grid structure
  type(verticalGrid_type), intent(in)  :: GV  !< Vertical grid structure
  type(unit_scale_type),   intent(in)  :: US  !< A dimensional unit scaling type
  type(time_type),         intent(in)  :: Time_present !< Model Time
  type(time_type),         intent(in)  :: dt  !< Time increment as a time-type
  type(mech_forcing),      intent(in), optional  :: forces !< MOM_forcing_type
  ! Local variables
  type(time_type) :: Stokes_Time
  integer :: i, j, b

  if (CS%WaveMethod == TESTPROF) then
    ! Do nothing
  elseif (CS%WaveMethod == SURFBANDS) then
    if (CS%DataSource == DATAOVR) then
      ! Updating Stokes drift time to center of time increment.
      !  This choice makes sense for the thermodynamics, but for the
      !  dynamics it may be more useful to update to the end of the
      !  time increment.
      Stokes_Time = Time_present + dt/2
      call Surface_Bands_by_data_override(Stokes_Time, G, GV, US, CS)
    elseif (CS%DataSource == COUPLER) then
      if (.not.present(FORCES)) then
        call MOM_error(FATAL,"The option SURFBAND = COUPLER can not be used with "//&
             "this driver. If you are using a coupled driver with a wave model then "//&
             "check the arguments in the subroutine call to Update_Surface_Waves, "//&
             "otherwise select another option for SURFBAND_SOURCE.")
      endif
      if (size(CS%WaveNum_Cen) /= size(forces%stk_wavenumbers)) then
        call MOM_error(FATAL, "Number of wavenumber bands in WW3 does not match that in MOM6. "//&
             "Make sure that STK_BAND_COUPLER in MOM6 input is equal to the number of bands in "//&
             "ww3_grid.inp, and that your mod_def.ww3 is up to date.")
      endif

      do b=1,CS%NumBands
        CS%WaveNum_Cen(b) = forces%stk_wavenumbers(b)
        !Interpolate from a grid to c grid
        do j=G%jsc,G%jec
          do I=G%iscB,G%iecB
            CS%STKx0(I,j,b) = 0.5*(forces%UStkb(i,j,b)+forces%UStkb(i+1,j,b))
          enddo
        enddo
        do J=G%jscB,G%jecB
          do i=G%isc,G%iec
            CS%STKY0(i,J,b) = 0.5*(forces%VStkb(i,j,b)+forces%VStkb(i,j+1,b))
          enddo
        enddo
        call pass_vector(CS%STKx0(:,:,b),CS%STKy0(:,:,b), G%Domain)
      enddo
      do j=G%jsc,G%jec
        do i=G%isc,G%iec
          CS%Omega_w2x(i,j)   = forces%omega_w2x(i,j)
          do b=1,CS%NumBands
            CS%UStk_Hb(i,j,b) = forces%UStkb(i,j,b)
            CS%VStk_Hb(i,j,b) = forces%VStkb(i,j,b)
          enddo
        enddo
      enddo
    elseif (CS%DataSource == INPUT) then
      do b=1,CS%NumBands
        do j=G%jsd,G%jed
          do I=G%isdB,G%iedB
            CS%STKx0(I,j,b) = CS%PrescribedSurfStkX(b)
          enddo
        enddo
        do J=G%jsdB, G%jedB
          do i=G%isd,G%ied
            CS%STKY0(i,J,b) = CS%PrescribedSurfStkY(b)
          enddo
        enddo
      enddo
    endif
  endif

end subroutine Update_Surface_Waves

!> Constructs the Stokes Drift profile on the model grid based on
!! desired coupling options
subroutine Update_Stokes_Drift(G, GV, US, CS, dz, ustar, dt, dynamics_step)
  type(wave_parameters_CS), pointer       :: CS    !< Wave parameter Control structure
  type(ocean_grid_type),    intent(inout) :: G     !< Grid structure
  type(verticalGrid_type),  intent(in)    :: GV    !< Vertical grid structure
  type(unit_scale_type),    intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: dz    !< Thickness in height units [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                            intent(in)    :: ustar !< Wind friction velocity [Z T-1 ~> m s-1].
  real, intent(in)                        :: dt    !< Time-step for computing Stokes-tendency [T ~> s]
  logical, intent(in)                     :: dynamics_step !< True if this call is on a dynamics step

  ! Local Variables
  real    :: Top, MidPoint, Bottom ! Positions within the layer [Z ~> m]
  real    :: level_thick ! The thickness of each layer [Z ~> m]
  real    :: DecayScale ! A vertical decay scale in the test profile [Z-1 ~> m-1]
  real    :: CMN_FAC  ! A nondimensional factor [nondim]
  real    :: WN       ! Model wavenumber [Z-1 ~> m-1]
  real    :: UStokes  ! A Stokes drift velocity [L T-1 ~> m s-1]
  real    :: PI       ! 3.1415926535... [nondim]
  real    :: La       ! The local Langmuir number [nondim]
  integer :: i, j, k, b
  real    :: I_dt     ! The inverse of the time step [T-1 ~> s-1]

  if (CS%WaveMethod==EFACTOR) return

  if (allocated(CS%US_x) .and. allocated(CS%US_y)) then
    call pass_vector(CS%US_x(:,:,:),CS%US_y(:,:,:), G%Domain)
  endif

  ! 1. If Test Profile Option is chosen
  !    Computing mid-point value from surface value and decay wavelength
  if (CS%WaveMethod==TESTPROF) then
    PI = 4.0*atan(1.0)
    DecayScale = 4.*PI / CS%TP_WVL !4pi
    do j=G%jsc,G%jec
      do I=G%iscB,G%iecB
        Bottom = 0.0
        MidPoint = 0.0
        do k = 1,GV%ke
          Top = Bottom
          MidPoint = Bottom - 0.25*(dz(i,j,k)+dz(i+1,j,k))
          Bottom = Bottom - 0.5*(dz(i,j,k)+dz(i+1,j,k))
          CS%Us_x(I,j,k) = CS%TP_STKX0*exp(MidPoint*DecayScale)
        enddo
      enddo
    enddo
    do J=G%jscB,G%jecB
      do i=G%isc,G%iec
        Bottom = 0.0
        MidPoint = 0.0
        do k = 1,GV%ke
          Top = Bottom
          MidPoint = Bottom - 0.25*(dz(i,j,k)+dz(i,j+1,k))
          Bottom = Bottom - 0.5*(dz(i,j,k)+dz(i,j+1,k))
          CS%Us_y(i,J,k) = CS%TP_STKY0*exp(MidPoint*DecayScale)
        enddo
      enddo
    enddo
    call pass_vector(CS%Us_x(:,:,:),CS%Us_y(:,:,:), G%Domain, To_All)
  ! 2. If Surface Bands is chosen
  !    In wavenumber mode compute integral for layer averaged Stokes drift.
  !    In frequency mode compuate value at midpoint.
  elseif (CS%WaveMethod==SURFBANDS) then
    CS%Us_x(:,:,:) = 0.0
    CS%Us_y(:,:,:) = 0.0
    CS%Us0_x(:,:) = 0.0
    CS%Us0_y(:,:) = 0.0
    ! Computing X direction Stokes drift
    do j=G%jsc,G%jec
      do I=G%iscB,G%iecB
        ! 1. First compute the surface Stokes drift
        !    by summing over the partitions.
        do b = 1,CS%NumBands
          CS%US0_x(I,j) = CS%US0_x(I,j) + CS%STKx0(I,j,b)
        enddo
        ! 2. Second compute the level averaged Stokes drift
        bottom = 0.0
        do k = 1,GV%ke
          Top = Bottom
          level_thick = 0.5*(dz(i,j,k)+dz(i+1,j,k))
          MidPoint = Top - 0.5*level_thick
          Bottom = Top - level_thick

          if (CS%answer_date >= 20230101) then
            ! Use more accurate and numerically stable expressions that work even for vanished layers.
            do b = 1,CS%NumBands
              if (CS%PartitionMode == 0) then
                ! Average over a layer using the bin's central wavenumber.
                CMN_FAC = exp(2.*CS%WaveNum_Cen(b)*Top) * one_minus_exp_x(2.*CS%WaveNum_Cen(b)*level_thick)
              else
                ! Use an analytic expression for the average of an exponential over a layer
                WN = CS%Freq_Cen(b)**2 * CS%I_g_Earth
                CMN_FAC = exp(2.*WN*Top) * one_minus_exp_x(2.*WN*level_thick)
              endif
              CS%US_x(I,j,k) = CS%US_x(I,j,k) + CS%STKx0(I,j,b)*CMN_FAC
            enddo

          elseif (level_thick > CS%Stokes_min_thick_avg) then
            ! -> Stokes drift in thin layers not averaged.
            do b = 1,CS%NumBands
              if (CS%PartitionMode == 0) then
                ! In wavenumber we are averaging over level
                CMN_FAC = (exp(Top*2.*CS%WaveNum_Cen(b))-exp(Bottom*2.*CS%WaveNum_Cen(b))) &
                          / ((Top-Bottom)*(2.*CS%WaveNum_Cen(b)))
              else
                ! Use a numerical integration and then divide by layer thickness
                WN = CS%Freq_Cen(b)**2 / CS%g_Earth !bgr bug-fix missing g
                CMN_FAC = (exp(2.*WN*Top)-exp(2.*WN*Bottom)) / (2.*WN*(Top-Bottom))
              endif
              CS%US_x(I,j,k) = CS%US_x(I,j,k) + CS%STKx0(I,j,b)*CMN_FAC
            enddo
          else ! Take the value at the midpoint
            do b = 1,CS%NumBands
              if (CS%PartitionMode == 0) then
                CMN_FAC = exp(MidPoint * 2. * CS%WaveNum_Cen(b))
              else
                CMN_FAC = exp(MidPoint * 2. * CS%Freq_Cen(b)**2 / CS%g_Earth)
              endif
              CS%US_x(I,j,k) = CS%US_x(I,j,k) + CS%STKx0(I,j,b)*CMN_FAC
            enddo
          endif
        enddo
      enddo
    enddo

    ! Computing Y direction Stokes drift
    do J=G%jscB,G%jecB
      do i=G%isc,G%iec
        ! Set the surface value to that at z=0
        do b = 1,CS%NumBands
          CS%US0_y(i,J) = CS%US0_y(i,J) + CS%STKy0(i,J,b)
        enddo
        ! Compute the level averages.
        bottom = 0.0
        do k = 1,GV%ke
          Top = Bottom
          level_thick = 0.5*(dz(i,j,k)+dz(i,j+1,k))
          MidPoint = Top - 0.5*level_thick
          Bottom = Top - level_thick

          if (CS%answer_date >= 20230101) then
            ! Use more accurate and numerically stable expressions that work even for vanished layers.
            do b = 1,CS%NumBands
              if (CS%PartitionMode == 0) then
                ! Average over a layer using the bin's central wavenumber.
                CMN_FAC = exp(2.*CS%WaveNum_Cen(b)*Top) * one_minus_exp_x(2.*CS%WaveNum_Cen(b)*level_thick)
              else
                ! Use an analytic expression for the average of an exponential over a layer
                WN = CS%Freq_Cen(b)**2 * CS%I_g_Earth
                CMN_FAC = exp(2.*WN*Top) * one_minus_exp_x(2.*WN*level_thick)
              endif
              CS%US_y(i,J,k) = CS%US_y(i,J,k) + CS%STKy0(i,J,b)*CMN_FAC
            enddo
          elseif (level_thick > CS%Stokes_min_thick_avg) then
            ! -> Stokes drift in thin layers not averaged.
            do b = 1,CS%NumBands
              if (CS%PartitionMode == 0) then
                ! In wavenumber we are averaging over level
                CMN_FAC = (exp(Top*2.*CS%WaveNum_Cen(b))-exp(Bottom*2.*CS%WaveNum_Cen(b))) &
                          / ((Top-Bottom)*(2.*CS%WaveNum_Cen(b)))
              else
                ! Use a numerical integration and then divide by layer thickness
                WN = CS%Freq_Cen(b)**2 / CS%g_Earth !bgr bug-fix missing g
                CMN_FAC = (exp(2.*WN*Top)-exp(2.*WN*Bottom)) / (2.*WN*(Top-Bottom))
              endif
              CS%US_y(i,J,k) = CS%US_y(i,J,k) + CS%STKy0(i,J,b)*CMN_FAC
            enddo
          else ! Take the value at the midpoint
            do b = 1,CS%NumBands
              if (CS%PartitionMode == 0) then
                CMN_FAC = exp(MidPoint*2.*CS%WaveNum_Cen(b))
              else
                CMN_FAC = exp(MidPoint * 2. * CS%Freq_Cen(b)**2 / CS%g_Earth)
              endif
              CS%US_y(i,J,k) = CS%US_y(i,J,k) + CS%STKy0(i,J,b)*CMN_FAC
            enddo
          endif
        enddo
      enddo
    enddo
    call pass_vector(CS%Us_x(:,:,:),CS%Us_y(:,:,:), G%Domain, To_All)
    call pass_vector(CS%Us0_x(:,:),CS%Us0_y(:,:), G%Domain)
  elseif (CS%WaveMethod == DHH85) then
    if (.not.(CS%StaticWaves .and. CS%DHH85_is_set)) then
      do j=G%jsc,G%jec
        do I=G%iscB,G%iecB
          bottom = 0.0
          do k = 1,GV%ke
            Top = Bottom
            MidPoint = Top - 0.25*(dz(i,j,k)+dz(i+1,j,k))
            Bottom = Top - 0.5*(dz(i,j,k)+dz(i+1,j,k))
            !bgr note that this is using a u-point I on h-point ustar
            !    this code has only been previous used for uniform
            !    grid cases.  This needs fixed if DHH85 is used for non
            !    uniform cases.
            call DHH85_mid(GV, US, CS, MidPoint, UStokes)
            ! Putting into x-direction (no option for direction
            CS%US_x(I,j,k) = UStokes
          enddo
        enddo
      enddo
      do J=G%jscB,G%jecB
        do i=G%isc,G%iec
          Bottom = 0.0
          do k = 1,GV%ke
            Top = Bottom
            MidPoint = Bottom - 0.25*(dz(i,j,k)+dz(i,j+1,k))
            Bottom = Bottom - 0.5*(dz(i,j,k)+dz(i,j+1,k))
            !bgr note that this is using a v-point J on h-point ustar
            !    this code has only been previous used for uniform
            !    grid cases.  This needs fixed if DHH85 is used for non
            !    uniform cases.
            ! call DHH85_mid(GV, US, CS, Midpoint, UStokes)
            ! Putting into x-direction, so setting y direction to 0
            CS%US_y(i,J,k) = 0.0
            ! For rotational symmetry there should be the option for this to become = UStokes
            !    bgr - see note above, but this is true
            !          if this is used for anything
            !          other than simple LES comparison
          enddo
        enddo
      enddo
      CS%DHH85_is_set = .true.
    endif
    call pass_vector(CS%Us_x(:,:,:),CS%Us_y(:,:,:), G%Domain)
  else! Keep this else, fallback to 0 Stokes drift
    CS%Us_x(:,:,:) = 0.
    CS%Us_y(:,:,:) = 0.
  endif

  ! Turbulent Langmuir number is computed here and available to use anywhere.
  ! SL Langmuir number requires mixing layer depth, and therefore is computed
  ! in the routine it is needed by (e.g. KPP or ePBL).
  do j=G%jsc, G%jec
    do i=G%isc,G%iec
      call get_Langmuir_Number( La, G, GV, US, dz(i,j,1), ustar(i,j), i, j, &
             dz(i,j,:), CS, Override_MA=.false.)
      CS%La_turb(i,j) = La
    enddo
  enddo

  ! Finding tendency of Stokes drift over the time step to apply
  !  as an acceleration to the models current.
  if ( dynamics_step .and. CS%Stokes_DDT ) then
    I_dt = 1.0 / dt
    CS%ddt_us_x(:,:,:) = (CS%US_x(:,:,:) - CS%US_x_prev(:,:,:)) * I_dt
    CS%ddt_us_y(:,:,:) = (CS%US_y(:,:,:) - CS%US_y_prev(:,:,:)) * I_dt
    CS%US_x_prev(:,:,:) = CS%US_x(:,:,:)
    CS%US_y_prev(:,:,:) = CS%US_y(:,:,:)
  endif

  ! Output any desired quantities
  if (CS%id_surfacestokes_y>0) &
    call post_data(CS%id_surfacestokes_y, CS%us0_y, CS%diag)
  if (CS%id_surfacestokes_x>0) &
    call post_data(CS%id_surfacestokes_x, CS%us0_x, CS%diag)
  if (CS%id_3dstokes_y>0) &
    call post_data(CS%id_3dstokes_y, CS%us_y, CS%diag)
  if (CS%id_3dstokes_x>0) &
    call post_data(CS%id_3dstokes_x, CS%us_x, CS%diag)
  if (CS%Stokes_DDT) then
    if (CS%id_ddt_3dstokes_x>0) &
      call post_data(CS%id_ddt_3dstokes_x, CS%ddt_us_x, CS%diag)
    if (CS%id_ddt_3dstokes_y>0) &
      call post_data(CS%id_ddt_3dstokes_y, CS%ddt_us_y, CS%diag)
    if (CS%id_3dstokes_x_from_ddt>0) &
      call post_data(CS%id_3dstokes_x_from_ddt, CS%us_x_from_ddt, CS%diag)
    if (CS%id_3dstokes_y_from_ddt>0) &
      call post_data(CS%id_3dstokes_y_from_ddt, CS%us_y_from_ddt, CS%diag)
  endif
  if (CS%id_La_turb>0) &
    call post_data(CS%id_La_turb, CS%La_turb, CS%diag)

end subroutine Update_Stokes_Drift

!> Return the value of (1 - exp(-x))/x [nondim], using an accurate expression for small values of x.
real function one_minus_exp_x(x)
  real, intent(in) :: x !< The argument of the function ((1 - exp(-x))/x) [nondim]
  real, parameter :: C1_6 = 1.0/6.0  ! A rational fraction [nondim]
  if (abs(x) <= 2.0e-5) then
    ! The Taylor series expression for exp(-x) gives a more accurate expression for 64-bit reals.
    one_minus_exp_x = 1.0 - x * (0.5 - C1_6*x)
  else
    one_minus_exp_x = (1.0 - exp(-x)) / x
  endif
end function one_minus_exp_x

!> Return the value of (1 - exp(-x)) [nondim], using an accurate expression for small values of x.
real function one_minus_exp(x)
  real, intent(in) :: x !< The argument of the function ((1 - exp(-x))/x) [nondim]
  real, parameter :: C1_6 = 1.0/6.0  ! A rational fraction [nondim]
  if (abs(x) <= 2.0e-5) then
    ! The Taylor series expression for exp(-x) gives a more accurate expression for 64-bit reals.
    one_minus_exp = x * (1.0 - x * (0.5 - C1_6*x))
  else
    one_minus_exp = 1.0 - exp(-x)
  endif
end function one_minus_exp

!> A subroutine to fill the Stokes drift from a NetCDF file
!! using the data_override procedures.
subroutine Surface_Bands_by_data_override(Time, G, GV, US, CS)
  type(time_type),          intent(in) :: Time       !< Time to get Stokes drift bands
  type(wave_parameters_CS), pointer    :: CS         !< Wave structure
  type(ocean_grid_type), intent(inout) :: G          !< Grid structure
  type(verticalGrid_type),  intent(in) :: GV         !< Vertical grid structure
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type

  ! Local variables
  real    :: temp_x(SZI_(G),SZJ_(G)) ! Pseudo-zonal Stokes drift of band at h-points [L T-1 ~> m s-1]
  real    :: temp_y(SZI_(G),SZJ_(G)) ! Pseudo-meridional Stokes drift of band at h-points [L T-1 ~> m s-1]
  integer, dimension(4) :: sizes    ! The sizes of the various dimensions of the variable.
  character(len=48) :: dim_name(4)  ! The names of the dimensions of the variable.
  character(len=20) :: varname      ! The name of an input variable for data override.
  real :: PI       ! 3.1415926535... [nondim]
  logical :: wavenumber_exists
  integer :: ndims, b, i, j

  if (.not.CS%DataOver_initialized) then
    call data_override_init(G%Domain)
    CS%DataOver_initialized = .true.

    if (.not.file_exists(CS%SurfBandFileName)) &
      call MOM_error(FATAL, "MOM_wave_interface is unable to find file "//trim(CS%SurfBandFileName))

    ! Check if input has wavenumber or frequency variables.

    ! Read the number of wavenumber bands in the file, if the variable 'wavenumber' exists.
    call get_var_sizes(CS%SurfBandFileName, 'wavenumber', ndims, sizes, dim_names=dim_name)
    wavenumber_exists = (ndims > -1)

    if (.not.wavenumber_exists) then
      ! Read the number of frequency bands in the file, if the variable 'frequency' exists.
      call get_var_sizes(CS%SurfBandFileName, 'frequency', ndims, sizes, dim_names=dim_name)
      if (ndims < 0) &
        call MOM_error(FATAL, "error finding variable 'wavenumber' or 'frequency' in file "//&
                              trim(CS%SurfBandFileName)//" in MOM_wave_interface.")
    endif

    CS%NUMBANDS = sizes(1)
    ! Allocate the wavenumber bins
    allocate( CS%WaveNum_Cen(CS%NUMBANDS), source=0.0 )

    if (wavenumber_exists) then
      ! Wavenumbers found, so this file uses the old method:
      CS%PartitionMode = 0

      ! Reading wavenumber bins
      call read_variable(CS%SurfBandFileName, dim_name(1), CS%WaveNum_Cen, scale=US%Z_to_m)

    else
      ! Frequencies found, so this file uses the newer method:
      CS%PartitionMode = 1

      ! Allocate the frequency bins
      allocate( CS%Freq_Cen(CS%NUMBANDS), source=0.0 )

      ! Reading frequencies
      PI = 4.0*atan(1.0)
      call read_variable(CS%SurfBandFileName, dim_name(1), CS%Freq_Cen, scale=2.*PI*US%T_to_s)

      do b = 1,CS%NumBands
        CS%WaveNum_Cen(b) = CS%Freq_Cen(b)**2 / CS%g_Earth
      enddo
    endif

    if (.not.allocated(CS%STKx0)) then
      allocate( CS%STKx0(G%isdB:G%iedB,G%jsd:G%jed,CS%NUMBANDS), source=0.0 )
    endif
    if (.not.allocated(CS%STKy0)) then
      allocate( CS%STKy0(G%isd:G%ied,G%jsdB:G%jedB,CS%NUMBANDS), source=0.0 )
    endif
  endif

  do b = 1,CS%NumBands
    temp_x(:,:) = 0.0
    temp_y(:,:) = 0.0
    varname = '                    '
    write(varname, "(A3,I0)") 'Usx', b
    call data_override(G%Domain, trim(varname), temp_x, Time, scale=US%m_s_to_L_T)
    varname = '                    '
    write(varname, "(A3,I0)") 'Usy', b
    call data_override(G%Domain, trim(varname), temp_y, Time, scale=US%m_s_to_L_T)
    ! Update halo on h-grid
    call pass_vector(temp_x, temp_y, G%Domain, To_All, AGRID)
    ! Filter land values
    do j = G%jsd,G%jed
      do i = G%Isd,G%Ied
        if ((abs(temp_x(i,j)) > CS%land_speed) .or. (abs(temp_y(i,j)) > CS%land_speed)) then
          ! Assume land-mask and zero out
          temp_x(i,j) = 0.0
          temp_y(i,j) = 0.0
        endif
      enddo
    enddo

    ! Interpolate to u/v grids
    do j = G%jsc,G%jec
      do I = G%IscB,G%IecB
        CS%STKx0(I,j,b) = 0.5 * (temp_x(i,j) + temp_x(i+1,j))
      enddo
    enddo
    do J = G%JscB,G%JecB
      do i = G%isc,G%iec
        CS%STKy0(i,J,b) = 0.5 * (temp_y(i,j) + temp_y(i,j+1))
      enddo
    enddo
  enddo !Closes b-loop

  ! Update halo on u/v grids
  call pass_vector(CS%STKx0(:,:,:), CS%STKy0(:,:,:), G%Domain, To_ALL)

end subroutine Surface_Bands_by_data_override

!> Interface to get Langmuir number based on options stored in wave structure
!!
!! Note this can be called with an unallocated Waves pointer, which is okay if we
!!  want the wind-speed only dependent Langmuir number.  Therefore, we need to be
!!  careful about what we try to access here.
subroutine get_Langmuir_Number( LA, G, GV, US, HBL, ustar, i, j, dz, Waves, &
                                U_H, V_H, Override_MA )
  type(ocean_grid_type),     intent(in)  :: G     !< Ocean grid structure
  type(verticalGrid_type),   intent(in)  :: GV    !< Ocean vertical grid structure
  real,                      intent(out) :: LA    !< Langmuir number [nondim]
  type(unit_scale_type),     intent(in)  :: US    !< A dimensional unit scaling type
  real,                      intent(in)  :: HBL   !< (Positive) thickness of boundary layer [Z ~> m]
  real,                      intent(in)  :: ustar !< Friction velocity [Z T-1 ~> m s-1]
  integer,                   intent(in)  :: i     !< Meridional index of h-point
  integer,                   intent(in)  :: j     !< Zonal index of h-point
  real, dimension(SZK_(GV)), intent(in)  :: dz    !< Grid layer thickness [Z ~> m]
  type(Wave_parameters_CS),  pointer     :: Waves !< Surface wave control structure.
  real, dimension(SZK_(GV)), &
                   optional, intent(in)  :: U_H   !< Zonal velocity at H point [L T-1 ~> m s-1] or [m s-1]
  real, dimension(SZK_(GV)), &
                   optional, intent(in)  :: V_H   !< Meridional velocity at H point [L T-1 ~> m s-1] or [m s-1]
  logical,         optional, intent(in)  :: Override_MA !< Override to use misalignment in LA
                                                  !! calculation. This can be used if diagnostic
                                                  !! LA outputs are desired that are different than
                                                  !! those used by the dynamical model.


!Local Variables
  real :: Top, Bottom, MidPoint  ! Positions within each layer [Z ~> m]
  real :: Dpt_LASL         ! Averaging depth for Stokes drift [Z ~> m]
  real :: ShearDirection   ! Shear angular direction from atan2 [radians]
  real :: WaveDirection    ! Wave angular direction from atan2 [radians]
  real :: LA_STKx, LA_STKy, LA_STK ! Stokes velocities in [L T-1 ~> m s-1]
  logical :: ContinueLoop, USE_MA
  real, dimension(SZK_(GV)) :: US_H, VS_H ! Profiles of Stokes velocities [L T-1 ~> m s-1]
  real, allocatable :: StkBand_X(:), StkBand_Y(:) ! Stokes drifts by band [L T-1 ~> m s-1]
  integer :: k, BB

  ! Compute averaging depth for Stokes drift (negative)
  Dpt_LASL = -1.0*max(Waves%LA_FracHBL*HBL, Waves%LA_HBL_min)

  USE_MA = Waves%LA_Misalignment
  if (present(Override_MA)) USE_MA = Override_MA

  ! If requesting to use misalignment in the Langmuir number compute the Shear Direction
  if (USE_MA) then
    if (.not.(present(U_H).and.present(V_H))) call MOM_error(FATAL, &
        "Get_LA_waves requested to consider misalignment, but velocities were not provided.")
    ContinueLoop = .true.
    bottom = 0.0
    do k = 1,GV%ke
      Top = Bottom
      MidPoint = Bottom + 0.5*dz(k)
      Bottom = Bottom + dz(k)

      if (Waves%LA_Misalign_bug) then
        ! Given the sign convention that Dpt_LASL is negative, the next line has a bug.
        if (MidPoint > Dpt_LASL .and. k > 1 .and. ContinueLoop) then
          ShearDirection = atan2(V_H(1)-V_H(k), U_H(1)-U_H(k))
          ContinueLoop = .false.
        endif
      else ! This version avoids the bug in the version above.
        if (MidPoint > abs(Dpt_LASL) .and. (k > 1) .and. ContinueLoop) then
          ShearDirection = atan2(V_H(1)-V_H(k), U_H(1)-U_H(k))
          ContinueLoop = .false.
        endif
      endif
    enddo
  endif

  if (Waves%WaveMethod==TESTPROF) then
    do k = 1,GV%ke
      US_H(k) = 0.5*(Waves%US_X(I,j,k)+Waves%US_X(I-1,j,k))
      VS_H(k) = 0.5*(Waves%US_Y(i,J,k)+Waves%US_Y(i,J-1,k))
    enddo
    call Get_SL_Average_Prof( GV, Dpt_LASL, dz, US_H, LA_STKx)
    call Get_SL_Average_Prof( GV, Dpt_LASL, dz, VS_H, LA_STKy)
    LA_STK = sqrt((LA_STKX*LA_STKX) + (LA_STKY*LA_STKY))
  elseif (Waves%WaveMethod==SURFBANDS) then
    allocate(StkBand_X(Waves%NumBands), StkBand_Y(Waves%NumBands))
    do bb = 1,Waves%NumBands
      StkBand_X(bb) = 0.5*(Waves%STKx0(I,j,bb)+Waves%STKx0(I-1,j,bb))
      StkBand_Y(bb) = 0.5*(Waves%STKy0(i,J,bb)+Waves%STKy0(i,J-1,bb))
    enddo
    call Get_SL_Average_Band(GV, Dpt_LASL, Waves%NumBands, Waves%WaveNum_Cen, StkBand_X, LA_STKx )
    call Get_SL_Average_Band(GV, Dpt_LASL, Waves%NumBands, Waves%WaveNum_Cen, StkBand_Y, LA_STKy )
    LA_STK = sqrt((LA_STKX**2) + (LA_STKY**2))
    deallocate(StkBand_X, StkBand_Y)
  elseif (Waves%WaveMethod==DHH85) then
    ! Temporarily integrating profile rather than spectrum for simplicity
    do k = 1,GV%ke
      US_H(k) = 0.5*(Waves%US_X(I,j,k)+Waves%US_X(I-1,j,k))
      VS_H(k) = 0.5*(Waves%US_Y(i,J,k)+Waves%US_Y(i,J-1,k))
    enddo
    call Get_SL_Average_Prof( GV, Dpt_LASL, dz, US_H, LA_STKx)
    call Get_SL_Average_Prof( GV, Dpt_LASL, dz, VS_H, LA_STKy)
    LA_STK = sqrt((LA_STKX**2) + (LA_STKY**2))
  elseif (Waves%WaveMethod==LF17) then
    call get_StokesSL_LiFoxKemper(ustar, HBL*Waves%LA_FracHBL, GV, US, Waves, LA_STK, LA)
  elseif (Waves%WaveMethod==Null_WaveMethod) then
    call MOM_error(FATAL, "Get_Langmuir_number called without defining a WaveMethod. "//&
                          "Suggest to make sure USE_LT is set/overridden to False or choose "//&
                          "a wave method (or set USE_LA_LI2016 to use statistical waves).")
  endif

  if (.not.(Waves%WaveMethod==LF17)) then
    ! This expression uses an arbitrary lower bound on Langmuir number.
    ! We shouldn't expect values lower than this, but there is also no good reason to cap it here
    ! other than to prevent large enhancements in unconstrained parts of the curve fit parameterizations.
    LA = max(Waves%La_min, sqrt(US%Z_to_L*ustar / (LA_STK + Waves%La_Stk_backgnd)))
  endif

  if (Use_MA) then
    WaveDirection = atan2(LA_STKy, LA_STKx)
    LA = LA / sqrt(max(1.e-8, cos( WaveDirection - ShearDirection)))
  endif

end subroutine get_Langmuir_Number

!> function to return the wave method string set in the param file
function get_wave_method(CS)
  character(:), allocatable :: get_wave_method
  type(wave_parameters_CS), pointer :: CS !< Control structure

  if (associated(CS)) then
    select case(CS%WaveMethod)
      case (Null_WaveMethod)
        get_wave_method = NULL_STRING
      case (TESTPROF)
        get_wave_method = TESTPROF_STRING
      case (SURFBANDS)
        get_wave_method = SURFBANDS_STRING
      case (DHH85)
        get_wave_method = DHH85_STRING
      case (LF17)
        get_wave_method = LF17_STRING
      case (EFACTOR)
        get_wave_method = EFACTOR_STRING
    end select
  else
    get_wave_method = NULL_STRING
  endif
end function get_wave_method

!> Get SL averaged Stokes drift from Li/FK 17 method
!!
!! Original description:
!! - This function returns the enhancement factor, given the 10-meter
!!   wind [m s-1], friction velocity [m s-1] and the boundary layer depth [m].
!!
!! Update (Jan/25):
!! - Converted from function to subroutine, now returns Langmuir number.
!! - Compute 10m wind internally, so only ustar and hbl need passed to
!!   subroutine.
!!
!! Qing Li, 160606
!! - BGR port from CVMix to MOM6 Jan/25/2017
!! - BGR change output to LA from Efactor
!! - BGR remove u10 input
!! - BGR note: fixed parameter values should be changed to "get_params"
subroutine get_StokesSL_LiFoxKemper(ustar, hbl, GV, US, CS, UStokes_SL, LA)
  real, intent(in)  :: ustar !< water-side surface friction velocity [Z T-1 ~> m s-1].
  real, intent(in)  :: hbl   !< boundary layer depth [Z ~> m].
  type(verticalGrid_type), intent(in) :: GV !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  type(wave_parameters_CS), pointer   :: CS  !< Wave parameter Control structure
  real, intent(out) :: UStokes_SL !< Surface layer averaged Stokes drift [L T-1 ~> m s-1]
  real, intent(out) :: LA    !< Langmuir number [nondim]
  ! Local variables
  ! parameters
  real, parameter :: u19p5_to_u10 = 1.075 ! ratio of U19.5 to U10 (Holthuijsen, 2007) [nondim]
  real, parameter :: fm_into_fp = 1.296 ! ratio of mean frequency to peak frequency for
                   ! Pierson-Moskowitz spectrum (Webb, 2011) [nondim]
  real, parameter :: us_to_u10 = 0.0162 ! ratio of surface Stokes drift to U10 [nondim]
  real, parameter :: r_loss = 0.667 ! loss ratio of Stokes transport [nondim]
  real :: UStokes  ! The surface Stokes drift [L T-1 ~> m s-1]
  real :: hm0      ! The significant wave height [Z ~> m]
  real :: fm       ! The mean wave frequency [T-1 ~> s-1]
  real :: fp       ! The peak wave frequency [T-1 ~> s-1]
  real :: kphil    ! A peak wavenumber in the Phillips spectrum [Z-1 ~> m-1]
  real :: kstar    ! A rescaled wavenumber? [Z-1 ~> m-1]
  real :: vstokes  ! The total Stokes transport [Z L T-1 ~> m2 s-1]
  real :: z0       ! The boundary layer depth [Z ~> m]
  real :: z0i      ! The inverse of the boundary layer depth [Z-1 ~> m-1]
  real :: r1, r2, r3, r4  ! Nondimensional ratios [nondim]
  real :: r5       ! A single expression that combines r2 and r4 [nondim]
  real :: root_2kz ! The square root of twice the peak wavenumber times the
                   ! boundary layer depth [nondim]
  real :: u10      ! The 10 m wind speed [L T-1 ~> m s-1]
  real :: PI       ! 3.1415926535... [nondim]

  PI = 4.0*atan(1.0)
  UStokes_sl = 0.0
  LA = 1.e8
  if (ustar > 0.0) then
    ! This code should be revised to minimize the number of divisions and cancel out common factors.

    ! Computing u10 based on u_star and COARE 3.5 relationships
    call ust_2_u10_coare3p5(ustar*sqrt(CS%rho_ocn/CS%rho_air), u10, GV, US, CS)
    ! surface Stokes drift
    UStokes = us_to_u10*u10
    !
    ! significant wave height from Pierson-Moskowitz spectrum (Bouws, 1998)
    hm0 = CS%SWH_from_u10sq * u10**2
    !
    ! peak frequency (PM, Bouws, 1998)
    fp = 0.877 * (US%L_to_Z*GV%g_Earth) / (2.0 * PI * u19p5_to_u10 * u10)
    !
    ! mean frequency
    fm = fm_into_fp * fp
    !
    ! total Stokes transport (a factor r_loss is applied to account
    !  for the effect of directional spreading, multidirectional waves
    !  and the use of PM peak frequency and PM significant wave height
    !  on estimating the Stokes transport)
    vstokes = 0.125 * PI * r_loss * US%Z_to_L * fm * hm0**2
    !
    ! the general peak wavenumber for Phillips' spectrum
    ! (Breivik et al., 2016) with correction of directional spreading
    kphil = 0.176 * UStokes / vstokes

    ! Combining all of the expressions above gives kPhil as the following
    ! where the first two lines are just a constant:
    ! kphil = ((0.176 * us_to_u10 * u19p5_to_u10) / &
    !          (0.5*0.125 * r_loss * fm_into_fp * 0.877 * CS%SWH_from_u10sq**2)) / &
    !         (GV%g_Earth * u10**2)

    ! surface layer
    z0 = abs(hbl)

    if (CS%answer_date < 20230102) then
      z0i = 1.0 / z0

      ! Surface layer averaged Stokes drift with Stokes drift profile
      ! estimated from Phillips' spectrum (Breivik et al., 2016)
      ! The directional spreading effect from Webb and Fox-Kemper, 2015 is also included.
      kstar = kphil * 2.56

      ! Terms 1 to 4, as written in the appendix of Li et al. (2017)
      r1 = ( 0.151 / kphil * z0i - 0.84 ) * &
           ( 1.0 - exp(-2.0 * kphil * z0) )
      r2 = -( 0.84 + 0.0591 / kphil * z0i ) * &
           sqrt( 2.0 * PI * kphil * z0 ) * &
           erfc( sqrt( 2.0 * kphil * z0 ) )
      r3 = ( 0.0632 / kstar * z0i + 0.125 ) * &
           (1.0 - exp(-2.0 * kstar * z0) )
      r4 = ( 0.125 + 0.0946 / kstar * z0i ) * &
           sqrt( 2.0 * PI * kstar * z0) * &
           erfc( sqrt( 2.0 * kstar * z0 ) )
      UStokes_sl = UStokes * (0.715 + r1 + r2 + r3 + r4)
    else
      ! The following is equivalent to the code above, but avoids singularities
      r1 = ( 0.302 - 1.68*(kphil*z0) ) * one_minus_exp_x(2.0 * (kphil * z0))
      r3 = ( 0.1264 + 0.64*(kphil*z0) ) * one_minus_exp_x(5.12 * (kphil * z0))

      root_2kz = sqrt(2.0 * kphil * z0)
      ! r2 = -( 0.84 + 0.0591*2.0 / (root_2kz**2) ) * sqrt(PI) * root_2kz * erfc( root_2kz )
      ! r4 = ( 0.2 + 0.059125*2.0 / (root_2kz**2) ) * sqrt(PI) * root_2kz * erfc( 1.6 * root_2kz )

      ! r5 = r2 + r4 (with a small correction to one coefficient to avoid a singularity when z0 = 0):
      ! The correction leads to <1% relative differences in (r2+r4) for root_2kz > 0.05, but without
      ! it the values of r2 + r4 are qualitatively wrong (>50% errors) for root_2kz < 0.0015 .
      !   It has been verified that these two expressions for r5 are the same to 6 decimal places for
      ! root_2kz  between 1e-10 and 1e-3, but that the first one degrades for smaller values.
      if (root_2kz > 1e-3) then
        r5 = sqrt(PI) * (root_2kz * (-0.84 * erfc(root_2kz) + 0.2 * erfc(1.6*root_2kz)) + &
                         0.1182 * (erfc(1.6*root_2kz) - erfc(root_2kz)) / root_2kz)
      else
        ! It is more accurate to replace erf with the first two terms of its Taylor series
        !  erf(z) = (2/sqrt(pi)) * z * (1. - (1/3)*z**2 + (1/10)*z**4 - (1/42)*z**6 + ...)
        ! and then cancel or combine common terms and drop negligibly small terms.
        r5 = -0.64*sqrt(PI)*root_2kz + (-0.14184 + 1.0839648 * root_2kz**2)
      endif
      UStokes_sl = UStokes * (0.715 + ((r1 + r3) + r5))
    endif

    if (UStokes_sl /= 0.0) LA = sqrt(US%Z_to_L*ustar / UStokes_sl)
  endif

end subroutine Get_StokesSL_LiFoxKemper

!> Get SL Averaged Stokes drift from a Stokes drift Profile
subroutine Get_SL_Average_Prof( GV, AvgDepth, dz, Profile, Average )
  type(verticalGrid_type),  &
       intent(in)   :: GV       !< Ocean vertical grid structure
  real, intent(in)  :: AvgDepth !< Depth to average over (negative) [Z ~> m]
  real, dimension(SZK_(GV)), &
       intent(in)   :: dz       !< Grid thickness [Z ~> m]
  real, dimension(SZK_(GV)), &
       intent(in)   :: Profile  !< Profile of quantity to be averaged in arbitrary units [A]
                                !! (used here for Stokes drift)
  real, intent(out) :: Average  !< Output quantity averaged over depth AvgDepth [A]
                                !! (used here for Stokes drift)
  !Local variables
  real :: Top, Bottom ! Depths, negative downward [Z ~> m]
  real :: Sum  ! The depth weighted vertical sum of a quantity [A Z ~> A m]
  integer :: k

  ! Initializing sum
  Sum = 0.0

  ! Integrate
  bottom = 0.0
  do k = 1, GV%ke
    Top = Bottom
    Bottom = Bottom - dz(k)
    if (AvgDepth < Bottom) then ! The whole cell is within H_LA
      Sum = Sum + Profile(k) * dz(k)
    elseif (AvgDepth < Top) then ! A partial cell is within H_LA
      Sum = Sum + Profile(k) * (Top-AvgDepth)
      exit
    else
      exit
    endif
  enddo

  ! Divide by AvgDepth or the depth in the column, whichever is smaller.
  if (abs(AvgDepth) <= abs(Bottom)) then
    Average = Sum / abs(AvgDepth)
  elseif (abs(Bottom) > 0.0) then
    Average = Sum / abs(Bottom)
  else
    Average = 0.0
  endif

end subroutine Get_SL_Average_Prof

!> Get SL averaged Stokes drift from the banded Spectrum method
subroutine Get_SL_Average_Band( GV, AvgDepth, NB, WaveNumbers, SurfStokes, Average )
  type(verticalGrid_type),  &
       intent(in)     :: GV          !< Ocean vertical grid
  real, intent(in)    :: AvgDepth    !< Depth to average over [Z ~> m].
  integer, intent(in) :: NB          !< Number of bands used
  real, dimension(NB), &
       intent(in)     :: WaveNumbers !< Wavenumber corresponding to each band [Z-1 ~> m-1]
  real, dimension(NB), &
       intent(in)     :: SurfStokes  !< Surface Stokes drift for each band [L T-1 ~> m s-1]
  real, intent(out)   :: Average     !< Output average Stokes drift over depth AvgDepth [L T-1 ~> m s-1]

  ! Local variables
  integer :: bb

  ! Loop over bands
  Average = 0.0
  do bb = 1, NB
    ! Factor includes analytical integration of e(2kz)
    !  - divided by (-H_LA) to get average from integral.
    Average = Average + SurfStokes(BB) * &
              (1.-EXP(-abs(AvgDepth * 2.0 * WaveNumbers(BB)))) / &
              abs(AvgDepth * 2.0 * WaveNumbers(BB))

    ! For accuracy when AvgDepth is small change the above to:
    ! Average = Average + SurfStokes(BB) * one_minus_exp_x(abs(AvgDepth * 2.0 * WaveNumbers(BB)))
  enddo

end subroutine Get_SL_Average_Band

!> Compute the Stokes drift at a given depth
!!
!! Taken from Qing Li (Brown)
!! use for comparing MOM6 simulation to his LES
!! computed at z mid point (I think) and not depth averaged.
!! Should be fine to integrate in frequency from 0.1 to sqrt(-0.2*grav*2pi/dz
subroutine DHH85_mid(GV, US, CS, zpt, UStokes)
  type(verticalGrid_type), intent(in)  :: GV  !< Ocean vertical grid
  type(unit_scale_type),   intent(in)  :: US  !< A dimensional unit scaling type
  type(wave_parameters_CS), pointer    :: CS  !< Wave parameter Control structure
  real, intent(in)  :: zpt   !< Depth to get Stokes drift [Z ~> m].
  real, intent(out) :: UStokes !< Stokes drift [L T-1 ~> m s-1]
  !
  real :: ann, Bnn, Snn, Cnn, Dnn ! Nondimensional factors [nondim]
  real :: omega_peak ! The peak wave frequency [T-1 ~> s-1]
  real :: omega      ! The average frequency in the band [T-1 ~> s-1]
  real :: domega     ! The width in frequency of the band [T-1 ~> s-1]
  real :: u10        ! The wind speed for this spectrum [Z T-1 ~> m s-1]
  real :: wavespec   ! The wave spectrum [L Z T ~> m2 s]
  real :: Stokes     ! The Stokes displacement per cycle [L ~> m]
  real :: PI         ! 3.1415926535... [nondim]
  integer :: Nomega  ! The number of wavenumber bands
  integer :: OI

  u10 = CS%WaveWind*US%L_to_Z

  !/
  NOmega = 1000
  domega = (CS%omega_max - CS%omega_min) / real(NOmega)

  !
  if (CS%WaveAgePeakFreq) then
    omega_peak = CS%g_Earth / (CS%WaveAge * u10)
  else
    PI = 4.0*atan(1.0)
    omega_peak = 2. * PI * 0.13 * CS%g_Earth / u10
  endif
  !/
  Ann = 0.006 * CS%WaveAge**(-0.55)
  Bnn = 1.0
  Snn = 0.08 * (1.0 + 4.0 * CS%WaveAge**3)
  Cnn = 1.7
  if (CS%WaveAge < 1.) then
    Cnn = Cnn - 6.0*log10(CS%WaveAge)
  endif
  !/
  UStokes = 0.0
  omega = CS%omega_min + 0.5*domega
  do oi = 1,nomega-1
    Dnn = exp ( -0.5 * (omega-omega_peak)**2 / (Snn**2 * omega_peak**2) )
    ! wavespec units [L Z T ~> m2 s]
    wavespec = US%Z_to_L * (Ann * CS%g_Earth**2 / (omega_peak*omega**4 ) ) * &
               exp(-bnn*(omega_peak/omega)**4)*Cnn**Dnn
    ! Stokes units [L ~> m] (multiply by frequency range for units of [L T-1 ~> m s-1])
    Stokes = 2.0 * wavespec * omega**3 * &
         exp( 2.0 * omega**2 * zpt / CS%g_Earth) / CS%g_Earth
    UStokes = UStokes + Stokes*domega
    omega = omega + domega
  enddo

end subroutine DHH85_mid

!> Explicit solver for Stokes mixing.
!! Still in development do not use.
subroutine StokesMixing(G, GV, dt, h, dz, u, v, Waves )
  type(ocean_grid_type), &
       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type), &
       intent(in)    :: GV    !< Ocean vertical grid
  real, intent(in)   :: dt    !< Time step of MOM6 [T ~> s] for explicit solver
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
       intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
       intent(in)    :: dz    !< Vertical distance between interfaces around a layer [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
       intent(inout) :: u     !< Velocity i-component [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
       intent(inout) :: v     !< Velocity j-component [L T-1 ~> m s-1]
  type(Wave_parameters_CS), &
       pointer       :: Waves !< Surface wave related control structure.
  ! Local variables
  real :: dTauUp, dTauDn ! Vertical momentum fluxes [H L T-2 ~> m2 s-2 or Pa]
  real :: h_lay  ! The layer thickness at a velocity point [H ~> m or kg m-2]
  real :: dz_lay ! The distance between interfaces at a velocity point [Z ~> m]
  integer :: i, j, k

! This is a template to think about down-Stokes mixing.
! This is not ready for use...

  do k = 1, GV%ke
    do j = G%jsc, G%jec
      do I = G%iscB, G%iecB
        h_lay = 0.5*(h(i,j,k)+h(i+1,j,k))
        dz_lay = 0.5*(dz(i,j,k)+dz(i+1,j,k))
        dTauUp = 0.0
        if (k > 1) &
          dTauUp = (0.5*(waves%Kvs(i,j,k)+waves%Kvs(i+1,j,k))) * &
               (waves%us_x(i,j,k-1)-waves%us_x(i,j,k)) / &
               (0.5*(dz_lay + 0.5*(dz(i,j,k-1)+dz(i+1,j,k-1)) ))
        dTauDn = 0.0
        if (k < GV%ke-1) &
          dTauDn = (0.5*(waves%Kvs(i,j,k+1)+waves%Kvs(i+1,j,k+1))) * &
               (waves%us_x(i,j,k)-waves%us_x(i,j,k+1)) / &
               (0.5*(dz_lay + 0.5*(dz(i,j,k+1)+dz(i+1,j,k+1)) ))
        u(i,j,k) = u(i,j,k) + dt * (dTauUp-dTauDn) / h_lay
      enddo
    enddo
  enddo

  do k = 1, GV%ke
    do J = G%jscB, G%jecB
      do i = G%isc, G%iec
        h_lay = 0.5*(h(i,j,k)+h(i,j+1,k))
        dz_lay = 0.5*(dz(i,j,k)+dz(i,j+1,k))
        dTauUp = 0.
        if (k > 1) &
          dTauUp = (0.5*(waves%Kvs(i,j,k)+waves%Kvs(i,j+1,k))) * &
               (waves%us_y(i,j,k-1)-waves%us_y(i,j,k)) / &
               (0.5*(dz_lay + 0.5*(dz(i,j,k-1)+dz(i,j+1,k-1)) ))
        dTauDn = 0.0
        if (k < GV%ke-1) &
          dTauDn = (0.5*(waves%Kvs(i,j,k+1)+waves%Kvs(i,j+1,k+1))) * &
               (waves%us_y(i,j,k)-waves%us_y(i,j,k+1)) / &
               (0.5*(dz_lay + 0.5*(dz(i,j,k+1)+dz(i,j+1,k+1)) ))
        v(i,J,k) = v(i,J,k) + dt * (dTauUp-dTauDn) / h_lay
      enddo
    enddo
  enddo

end subroutine StokesMixing

!> Solver to add Coriolis-Stokes to model
!! Still in development and not meant for general use.
!! Can be activated (with code intervention) for LES comparison
!! CHECK THAT RIGHT TIMESTEP IS PASSED IF YOU USE THIS**
!!
!! Not accessed in the standard code.
subroutine CoriolisStokes(G, GV, dt, h, u, v, Waves)
  type(ocean_grid_type), &
       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type), &
       intent(in)   :: GV     !< Ocean vertical grid
  real, intent(in)  :: dt     !< Time step of MOM6 [T ~> s]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
       intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
       intent(inout) :: u     !< Velocity i-component [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
       intent(inout) :: v     !< Velocity j-component [L T-1 ~> m s-1]
  type(Wave_parameters_CS), &
       pointer       :: Waves !< Surface wave related control structure.

  ! Local variables
  real :: DVel ! A rescaled velocity change [L T-2 ~> m s-2]
  integer :: i, j, k

  do k = 1, GV%ke
    do j = G%jsc, G%jec
      do I = G%iscB, G%iecB
        DVel = 0.25*((Waves%us_y(i,J-1,k)+Waves%us_y(i+1,J-1,k)) * G%CoriolisBu(I,J-1)) + &
               0.25*((Waves%us_y(i,J,k)+Waves%us_y(i+1,J,k)) * G%CoriolisBu(I,J))
        u(I,j,k) = u(I,j,k) + DVEL*dt
      enddo
    enddo
  enddo

  do k = 1, GV%ke
    do J = G%jscB, G%jecB
      do i = G%isc, G%iec
        DVel = 0.25*((Waves%us_x(I-1,j,k)+Waves%us_x(I-1,j+1,k)) * G%CoriolisBu(I-1,j)) + &
               0.25*((Waves%us_x(I,j,k)+Waves%us_x(I,j+1,k)) * G%CoriolisBu(I,J))
        v(i,J,k) = v(i,j,k) - DVEL*dt
      enddo
    enddo
  enddo
end subroutine CoriolisStokes

!> Computes tendency due to Stokes pressure gradient force anomaly
!! including analytical integration of Stokes shear using multiple-exponential decay
!! Stokes drift profile and vertical integration of the resulting pressure
!! anomaly to the total pressure gradient force
subroutine Stokes_PGF(G, GV, US, dz, u, v, PFu_Stokes, PFv_Stokes, CS )
  type(ocean_grid_type), &
       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type), &
       intent(in)    :: GV    !< Ocean vertical grid
  type(unit_scale_type), &
       intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),&
       intent(in)    :: dz      !< Layer thicknesses in height units [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
       intent(in) :: u          !< Lagrangian Velocity i-component [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
       intent(in) :: v          !< Lagrangian Velocity j-component [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
       intent(out) :: PFu_Stokes !< PGF Stokes-shear i-component [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
       intent(out) :: PFv_Stokes !< PGF Stokes-shear j-component [L T-2 ~> m s-2]
  type(Wave_parameters_CS), &
       pointer       :: CS !< Surface wave related control structure.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: P_deltaStokes_L ! The Stokes induced pressure anomaly,
                                                              ! layer averaged [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: P_deltaStokes_i ! The Stokes induced pressure anomaly
                                                              ! at interfaces [L2 T-2 ~> m2 s-2]
  real :: P_Stokes_l, P_Stokes_r ! Stokes-induced pressure anomaly over layer (left/right of point) [L2 T-2 ~> m2 s-2]
  real :: P_Stokes_l0, P_Stokes_r0 ! Stokes-induced pressure anomaly at interface
                                   ! (left/right of point) [L2 T-2 ~> m2 s-2]
  real :: dP_Stokes_l_dz, dP_Stokes_r_dz ! Contribution of layer to integrated Stokes pressure anomaly for summation
                                         ! (left/right of point) [Z L2 T-2 ~> m3 s-2]
  real :: dP_lay_Stokes_l, dP_lay_Stokes_r ! Contribution of layer to integrated Stokes pressure anomaly for summation
                                         ! (left/right of point) [L2 T-2 ~> m2 s-2]
  real :: dP_Stokes_l, dP_Stokes_r ! Net increment of Stokes pressure anomaly across layer for summation
                                   ! (left/right of point) [L2 T-2 ~> m2 s-2]
  real :: uE_l, uE_r, vE_l, vE_r ! Eulerian velocity components (left/right of point) [L T-1 ~> m s-1]
  real :: uS0_l, uS0_r, vS0_l, vS0_r ! Surface Stokes velocity components (left/right of point) [L T-1 ~> m s-1]
  real :: zi_l(SZK_(G)+1), zi_r(SZK_(G)+1)   ! The height of the edges of the cells (left/right of point) [Z ~> m].
  real :: idz_l(SZK_(G)), idz_r(SZK_(G)) ! The inverse thickness of the cells (left/right of point) [Z-1 ~> m-1]
  real :: h_l, h_r   ! The thickness of the cell (left/right of point) [Z ~> m].
  real :: exp_top    ! The decay of the surface stokes drift to the interface atop a layer [nondim]
  real :: dexp2kzL, dexp4kzL, dexp2kzR, dexp4kzR ! Analytical evaluation of multi-exponential decay
                                              ! contribution to Stokes pressure anomalies [nondim].
  real :: TwoK, FourK   ! Wavenumbers multiplied by a factor [Z-1 ~> m-1]
  real :: iTwoK, iFourK ! Inverses of wavenumbers [Z ~> m]

  integer :: i, j, k, l

  !---------------------------------------------------------------
  ! Compute the Stokes contribution to the pressure gradient force
  !---------------------------------------------------------------
  ! Notes on the algorithm/code:
  ! This code requires computing velocities at bounding h points
  ! of the u/v points to get the pressure-gradient. In this
  ! implementation there are several redundant calculations as the
  ! left/right points are computed at each cell while integrating
  ! in the vertical, requiring about twice the calculations.  The
  ! velocities at the tracer points could be precomputed and
  ! stored, but this would require more memory and cycling through
  ! large 3d arrays while computing the pressures. This could be
  ! explored as a way to speed up this code.
  !---------------------------------------------------------------

  PFu_Stokes(:,:,:) = 0.0
  PFv_Stokes(:,:,:) = 0.0
  if (CS%id_P_deltaStokes_i > 0) P_deltaStokes_i(:,:,:) = 0.0
  if (CS%id_P_deltaStokes_L > 0) P_deltaStokes_L(:,:,:) = 0.0

  ! First compute PGFu.  The Stokes-induced pressure anomaly diagnostic is stored from this calculation.
  ! > Seeking PGFx at (I,j), meaning we need to compute pressure at h-points (i,j) and (i+1,j).
  !   UL(i,j)   -> found as average of I-1 & I on j
  !   UR(i+1,j) -> found as average of I & I+1 on j
  !   VL(i,j)   -> found on i as average of J-1 & J
  !   VR(i+1,j) -> found on i+1 as average of J-1 & J
  !
  do j = G%jsc, G%jec ; do I = G%iscB, G%iecB
    if (G%mask2dCu(I,j)>0.5) then
      P_Stokes_l0 = 0.0
      P_Stokes_r0 = 0.0
      ! We don't need to precompute the grid in physical space arrays and could have done this during
      ! the next loop, but this gives flexibility if the loop directions (integrations) are performed
      ! upwards instead of downwards (it seems downwards is the better approach).
      zi_l(1) = 0.0
      zi_r(1) = 0.0
      do k = 1, G%ke
        h_l = dz(i,j,k)
        h_r = dz(i+1,j,k)
        zi_l(k+1) = zi_l(k) - h_l
        zi_r(k+1) = zi_r(k) - h_r
        if (.not.CS%robust_Stokes_PGF) then
          ! When the code is properly refactored, the following hard-coded constants are unnecessary.
          Idz_l(k) = 1./max(0.1*US%m_to_Z, h_l)
          Idz_r(k) = 1./max(0.1*US%m_to_Z, h_r)
        endif
      enddo
      do k = 1,G%ke
        ! Computing (left/right) Eulerian velocities assuming the velocity passed to this routine is the
        ! Lagrangian velocity.  This requires the wave acceleration terms to be activated together.
        uE_l = 0.5*((u(I-1,j,k)-CS%Us_x(I-1,j,k))*G%mask2dCu(I-1,j) + &
                    (u(I,j,k)-CS%Us_x(I-1,j,k))*G%mask2dCu(I,j))
        uE_r = 0.5*((u(I,j,k)-CS%Us_x(I,j,k))*G%mask2dCu(I,j) + &
                    (u(I+1,j,k)-CS%Us_x(I+1,j,k))*G%mask2dCu(I+1,j))
        vE_l = 0.5*((v(i,J-1,k)-CS%Us_y(i,J-1,k))*G%mask2dCv(i,J-1) + &
                    (v(i,J,k)-CS%Us_y(i,J,k))*G%mask2dCv(i,J))
        vE_r = 0.5*((v(i+1,J-1,k)-CS%Us_y(i+1,J-1,k))*G%mask2dCv(i+1,J-1) + &
                    (v(i+1,J,k)-CS%Us_y(i+1,J,k))*G%mask2dCv(i+1,J))

        dP_Stokes_l_dz = 0.0
        dP_Stokes_r_dz = 0.0
        dP_Stokes_l = 0.0
        dP_Stokes_r = 0.0

        do l = 1, CS%numbands

          ! Computing (left/right) surface Stokes drift velocities at wavenumber band
          uS0_l = 0.5*(CS%Stkx0(I-1,j,l)*G%mask2dCu(I-1,j) + &
                       CS%Stkx0(I,j,l)*G%mask2dCu(I,j))
          uS0_r = 0.5*(CS%Stkx0(I,j,l)*G%mask2dCu(I,j) + &
                       CS%Stkx0(I+1,j,l)*G%mask2dCu(I+1,j))
          vS0_l = 0.5*(CS%Stky0(i,J-1,l)*G%mask2dCv(i,J-1) + &
                       CS%Stky0(i,J,l)*G%mask2dCv(i,J))
          vS0_r = 0.5*(CS%Stky0(i+1,J-1,l)*G%mask2dCv(i+1,J-1) + &
                       CS%Stky0(i+1,J,l)*G%mask2dCv(i+1,J))

          ! Wavenumber terms that are useful to simplify the pressure calculations
          TwoK = 2.*CS%WaveNum_Cen(l)
          FourK = 2.*TwoK
          if (.not.CS%robust_Stokes_PGF) then
            iTwoK = 1. / TwoK
            iFourK = 1. / FourK
          endif

          ! Compute Pressure at interface and integrated over layer on left/right bounding points.
          ! These are summed over wavenumber bands.
          if (G%mask2dT(i,j)>0.5) then
            if (.not.CS%robust_Stokes_PGF) then
              dexp2kzL = exp(TwoK*zi_l(k))-exp(TwoK*zi_l(k+1))
              dexp4kzL = exp(FourK*zi_l(k))-exp(FourK*zi_l(k+1))
              dP_Stokes_l_dz = dP_Stokes_l_dz + &
                               ((uE_l*uS0_l+vE_l*vS0_l)*iTwoK*dexp2kzL + 0.5*(uS0_l*uS0_l+vS0_l*vS0_l)*iFourK*dexp4kzL)
              dP_Stokes_l = dP_Stokes_l + (uE_l*uS0_l+vE_l*vS0_l)*dexp2kzL + 0.5*(uS0_l*uS0_l+vS0_l*vS0_l)*dexp4kzL
            else  ! These expressions are equivalent to those above for thick layers, but more accurate for thin layers.
              exp_top = exp(TwoK*zi_l(k))
              dP_lay_Stokes_l = dP_lay_Stokes_l + &
                  ((((uE_l*uS0_l)+(vE_l*vS0_l)) * exp_top) * one_minus_exp_x(TwoK*dz(i,j,k)) + &
                   (0.5*((uS0_l**2)+(vS0_l**2)) * exp_top**2) * one_minus_exp_x(FourK*dz(i,j,k)) )
              dP_Stokes_l = dP_Stokes_l + &
                  ((((uE_l*uS0_l)+(vE_l*vS0_l)) * exp_top) * one_minus_exp(TwoK*dz(i,j,k)) + &
                   (0.5*((uS0_l**2)+(vS0_l**2)) * exp_top**2) * one_minus_exp(FourK*dz(i,j,k)) )
            endif
          endif
          if (G%mask2dT(i+1,j)>0.5) then
            if (.not.CS%robust_Stokes_PGF) then
              dexp2kzR = exp(TwoK*zi_r(k))-exp(TwoK*zi_r(k+1))
              dexp4kzR = exp(FourK*zi_r(k))-exp(FourK*zi_r(k+1))
              dP_Stokes_r_dz = dP_Stokes_r_dz + &
                               ((uE_r*uS0_r+vE_r*vS0_r)*iTwoK*dexp2kzR + 0.5*(uS0_l*uS0_l+vS0_l*vS0_l)*iFourK*dexp4kzR)
              dP_Stokes_r = dP_Stokes_r + (uE_r*uS0_r+vE_r*vS0_r)*dexp2kzR + 0.5*(uS0_l*uS0_l+vS0_l*vS0_l)*dexp4kzR
            else  ! These expressions are equivalent to those above for thick layers, but more accurate for thin layers.
              exp_top = exp(TwoK*zi_r(k))
              dP_lay_Stokes_r = dP_lay_Stokes_r + &
                  ((((uE_r*uS0_r)+(vE_r*vS0_r)) * exp_top) * one_minus_exp_x(TwoK*dz(i+1,j,k)) + &
                   (0.5*((uS0_r**2)+(vS0_r**2)) * exp_top**2) * one_minus_exp_x(FourK*dz(i+1,j,k)) )
              dP_Stokes_r = dP_Stokes_r + &
                  ((((uE_r*uS0_r)+(vE_r*vS0_r)) * exp_top) * one_minus_exp(TwoK*dz(i+1,j,k)) + &
                   (0.5*((uS0_r**2)+(vS0_r**2)) * exp_top**2) * one_minus_exp(FourK*dz(i+1,j,k)) )
            endif
          endif
        enddo

        ! Summing PF over bands
        ! > Increment the Layer averaged pressure
        if (.not.CS%robust_Stokes_PGF) then
          P_Stokes_l = P_Stokes_l0 + dP_Stokes_l_dz*Idz_l(k)
          P_Stokes_r = P_Stokes_r0 + dP_Stokes_r_dz*Idz_r(k)
        else
          P_Stokes_l = P_Stokes_l0 + dP_lay_Stokes_l
          P_Stokes_r = P_Stokes_r0 + dP_lay_Stokes_r
        endif

        ! > Increment the Interface pressure
        P_Stokes_l0 = P_Stokes_l0 + dP_Stokes_l
        P_Stokes_r0 = P_Stokes_r0 + dP_Stokes_r

        ! Pressure force anomaly is finite difference across the cell.
        PFu_Stokes(I,j,k) = (P_Stokes_r - P_Stokes_l)*G%IdxCu(I,j)

        ! Choose to output the pressure delta on the h-points from the PFu calculation.
        if (G%mask2dT(i,j)>0.5 .and. CS%id_P_deltaStokes_L > 0) P_deltaStokes_L(i,j,k) = P_Stokes_l
        if (G%mask2dT(i,j)>0.5 .and. CS%id_P_deltaStokes_i > 0) P_deltaStokes_i(i,j,k+1) = P_Stokes_l0

      enddo
    endif
  enddo ;  enddo

  ! Next compute PGFv.  The Stokes-induced pressure anomaly diagnostic is stored from this calculation.
  ! > Seeking PGFy at (i,J), meaning we need to compute pressure at h-points (i,j) and (i,j+1).
  !   UL(i,j)   -> found as average of I-1 & I on j
  !   UR(i,j+1) -> found as average of I-1 & I on j+1
  !   VL(i,j)   -> found on i as average of J-1 & J
  !   VR(i,j+1) -> found on i as average of J & J+1
  !
  do J = G%jscB, G%jecB ; do i = G%isc, G%iec
    if (G%mask2dCv(i,J)>0.5) then
      P_Stokes_l0 = 0.0
      P_Stokes_r0 = 0.0
      zi_l(1) = 0.0
      zi_r(1) = 0.0
      do k = 1, G%ke
        h_l = dz(i,j,k)
        h_r = dz(i,j+1,k)
        zi_l(k+1) = zi_l(k) - h_l
        zi_r(k+1) = zi_r(k) - h_r
        if (.not.CS%robust_Stokes_PGF) then
          ! When the code is properly refactored, the following hard-coded constants are unnecessary.
          Idz_l(k) = 1. / max(0.1*US%m_to_Z, h_l)
          Idz_r(k) = 1. / max(0.1*US%m_to_Z, h_r)
        endif
      enddo
      do k = 1,G%ke
        ! Computing (left/right) Eulerian velocities assuming the velocity passed to this routine is the
        ! Lagrangian velocity.  This requires the wave acceleration terms to be activated together.
        uE_l = 0.5*((u(I-1,j,k)-CS%Us_x(I-1,j,k))*G%mask2dCu(I-1,j) + &
                    (u(I,j,k)-CS%Us_x(I,j,k))*G%mask2dCu(I,j))
        uE_r = 0.5*((u(I-1,j+1,k)-CS%Us_x(I-1,j+1,k))*G%mask2dCu(I-1,j+1) + &
                    (u(I,j+1,k)-CS%Us_x(I,j+1,k))*G%mask2dCu(I,j+1))
        vE_l = 0.5*((v(i,J-1,k)-CS%Us_y(i,J-1,k))*G%mask2dCv(i,J-1) + &
                    (v(i,J,k)-CS%Us_y(i,J,k))*G%mask2dCv(i,J))
        vE_r = 0.5*((v(i,J,k)-CS%Us_y(i,J,k))*G%mask2dCv(i,J) + &
                    (v(i,J+1,k)-CS%Us_y(i,J+1,k))*G%mask2dCv(i,J+1))

        dP_Stokes_l_dz = 0.0
        dP_Stokes_r_dz = 0.0
        dP_Stokes_l = 0.0
        dP_Stokes_r = 0.0

        do l = 1, CS%numbands

          ! Computing (left/right) surface Stokes drift velocities at wavenumber band
          uS0_l = 0.5*(CS%Stkx0(I-1,j,l)*G%mask2dCu(I-1,j) + &
                       CS%Stkx0(I,j,l)*G%mask2dCu(I,j))
          uS0_r = 0.5*(CS%Stkx0(I-1,j+1,l)*G%mask2dCu(I-1,j+1) + &
                       CS%Stkx0(I,j+1,l)*G%mask2dCu(I,j+1))
          vS0_l = 0.5*(CS%Stky0(i,J-1,l)*G%mask2dCv(i,J-1) + &
                       CS%Stky0(i,J,l)*G%mask2dCv(i,J))
          vS0_r = 0.5*(CS%Stky0(i,J,l)*G%mask2dCv(i,J) + &
                       CS%Stky0(i,J+1,l)*G%mask2dCv(i,J+1))

          ! Wavenumber terms that are useful to simplify the pressure calculations
          TwoK = 2.*CS%WaveNum_Cen(l)
          FourK = 2.*TwoK
          if (.not.CS%robust_Stokes_PGF) then
            iTwoK = 1. / TwoK
            iFourK = 1. / FourK
          endif

          ! Compute Pressure at interface and integrated over layer on left/right bounding points.
          ! These are summed over wavenumber bands.
          if (G%mask2dT(i,j)>0.5) then
            if (.not.CS%robust_Stokes_PGF) then
              dexp2kzL = exp(TwoK*zi_l(k))-exp(TwoK*zi_l(k+1))
              dexp4kzL = exp(FourK*zi_l(k))-exp(FourK*zi_l(k+1))
              dP_Stokes_l_dz = dP_Stokes_l_dz + &
                               ((uE_l*uS0_l+vE_l*vS0_l)*iTwoK*dexp2kzL + 0.5*(uS0_l*uS0_l+vS0_l*vS0_l)*iFourK*dexp4kzL)
              dP_Stokes_l = dP_Stokes_l + (uE_l*uS0_l+vE_l*vS0_l)*dexp2kzL + 0.5*(uS0_l*uS0_l+vS0_l*vS0_l)*dexp4kzL
            else  ! These expressions are equivalent to those above for thick layers, but more accurate for thin layers.
              exp_top = exp(TwoK*zi_l(k))
              dP_lay_Stokes_l = dP_lay_Stokes_l + &
                  ((((uE_l*uS0_l)+(vE_l*vS0_l)) * exp_top) * one_minus_exp_x(TwoK*dz(i,j,k)) + &
                   (0.5*((uS0_l**2)+(vS0_l**2)) * exp_top**2) * one_minus_exp_x(FourK*dz(i,j,k)) )
              dP_Stokes_l = dP_Stokes_l + &
                  ((((uE_l*uS0_l)+(vE_l*vS0_l)) * exp_top) * one_minus_exp(TwoK*dz(i,j,k)) + &
                   (0.5*((uS0_l**2)+(vS0_l**2)) * exp_top**2) * one_minus_exp(FourK*dz(i,j,k)) )
            endif
          endif
          if (G%mask2dT(i,j+1)>0.5) then
            if (.not.CS%robust_Stokes_PGF) then
              dexp2kzR = exp(TwoK*zi_r(k))-exp(TwoK*zi_r(k+1))
              dexp4kzR = exp(FourK*zi_r(k))-exp(FourK*zi_r(k+1))
              dP_Stokes_r_dz = dP_Stokes_r_dz + &
                               ((uE_r*uS0_r+vE_r*vS0_r)*iTwoK*dexp2kzR + 0.5*(uS0_l*uS0_l+vS0_l*vS0_l)*iFourK*dexp4kzR)
              dP_Stokes_r = dP_Stokes_r + (uE_r*uS0_r+vE_r*vS0_r)*dexp2kzR + 0.5*(uS0_l*uS0_l+vS0_l*vS0_l)*dexp4kzR
            else  ! These expressions are equivalent to those above for thick layers, but more accurate for thin layers.
              exp_top = exp(TwoK*zi_r(k))
              dP_lay_Stokes_r = dP_lay_Stokes_r + &
                  ((((uE_r*uS0_r)+(vE_r*vS0_r)) * exp_top) * one_minus_exp_x(TwoK*dz(i,j+1,k)) + &
                   (0.5*((uS0_r**2)+(vS0_r**2)) * exp_top**2) * one_minus_exp_x(FourK*dz(i,j+1,k)) )
              dP_Stokes_r = dP_Stokes_r + &
                  ((((uE_r*uS0_r)+(vE_r*vS0_r)) * exp_top) * one_minus_exp(TwoK*dz(i,j+1,k)) + &
                   (0.5*((uS0_r**2)+(vS0_r**2)) * exp_top**2) * one_minus_exp(FourK*dz(i,j+1,k)) )
            endif
          endif
        enddo

        ! Summing PF over bands
        ! > Increment the Layer averaged pressure
        if (.not.CS%robust_Stokes_PGF) then
          P_Stokes_l = P_Stokes_l0 + dP_Stokes_l_dz*Idz_l(k)
          P_Stokes_r = P_Stokes_r0 + dP_Stokes_r_dz*Idz_r(k)
        else
          P_Stokes_l = P_Stokes_l0 + dP_lay_Stokes_l
          P_Stokes_r = P_Stokes_r0 + dP_lay_Stokes_r
        endif

        ! > Increment the Interface pressure
        P_Stokes_l0 = P_Stokes_l0 + dP_Stokes_l
        P_Stokes_r0 = P_Stokes_r0 + dP_Stokes_r

        ! Pressure force anomaly is finite difference across the cell.
        PFv_Stokes(I,j,k) = (P_Stokes_r - P_Stokes_l)*G%IdyCv(i,J)

      enddo
    endif
  enddo ; enddo

  if (CS%id_PFv_Stokes>0) &
    call post_data(CS%id_PFv_Stokes, PFv_Stokes, CS%diag)
  if (CS%id_PFu_Stokes>0) &
    call post_data(CS%id_PFu_Stokes, PFu_Stokes, CS%diag)
  if (CS%id_P_deltaStokes_L>0) &
    call post_data(CS%id_P_deltaStokes_L, P_deltaStokes_L, CS%diag)
  if (CS%id_P_deltaStokes_i>0) &
    call post_data(CS%id_P_deltaStokes_i, P_deltaStokes_i, CS%diag)

end subroutine Stokes_PGF

!> Computes wind speed from ustar_air based on COARE 3.5 Cd relationship
!! Probably doesn't belong in this module, but it is used here to estimate
!! wind speed for wind-wave relationships.  Should be a fine way to estimate
!! the neutral wind-speed as written here.
subroutine ust_2_u10_coare3p5(USTair, U10, GV, US, CS)
  real, intent(in)                    :: USTair !< Wind friction velocity [Z T-1 ~> m s-1]
  real, intent(out)                   :: U10    !< 10-m neutral wind speed [L T-1 ~> m s-1]
  type(verticalGrid_type), intent(in) :: GV     !< vertical grid type
  type(unit_scale_type),   intent(in) :: US     !< A dimensional unit scaling type
  type(wave_parameters_CS), pointer   :: CS     !< Wave parameter Control structure

  ! Local variables
  real :: z0sm, z0, z0rough  ! Roughness lengths [Z ~> m]
  real :: ten_m_scale ! The 10 m reference height, in rescaled units [Z ~> m]
  real :: I_ten_m_scale ! The inverse of the 10 m reference height, in rescaled units [Z-1 ~> m-1]
  real :: u10a  ! The previous guess for u10 [L T-1 ~> m s-1]
  real :: alpha ! The Charnock coeffient relating the wind friction velocity squared to the
                ! roughness length [nondim]
  real :: Cd    ! The drag coefficient [nondim]
  real :: I_sqrtCd  ! The inverse of the square root of the drag coefficient [nondim]
  real :: I_vonKar  ! The inverse of the von Karman coefficient [nondim]
  integer :: CT

  ! Uses empirical formula for z0 to convert ustar_air to u10 based on the
  !  COARE 3.5 paper (Edson et al., 2013)
  ! alpha=m*U10+b
  ! Note in Edson et al. 2013, eq. 13 m is given as 0.017.  However,
  ! m=0.0017 reproduces the curve in their figure 6.

  if (CS%vonKar < 0.0) call MOM_error(FATAL, &
    "ust_2_u10_coare3p5 called with a negative value of Waves%vonKar")

  z0sm = 0.11 * CS%nu_air / USTair ! Compute z0smooth from ustar guess
  u10a = 1000.0*US%m_s_to_L_T ! An insanely large upper bound for u10.

  if (CS%answer_date < 20230103) then
    u10 = US%Z_to_L*USTair / sqrt(0.001)  ! Guess for u10
    ten_m_scale = 10.0*US%m_to_Z
    CT=0
    do while (abs(u10a/u10 - 1.) > 0.001)
      CT=CT+1
      u10a = u10
      alpha = min(CS%Charnock_min, CS%Charnock_slope_U10 * u10 + CS%Charnock_intercept)
      z0rough = alpha * (US%Z_to_L*USTair)**2 / GV%g_Earth ! Compute z0rough from ustar guess
      z0 = z0sm + z0rough
      Cd = ( CS%vonKar / log(ten_m_scale / z0) )**2 ! Compute Cd from derived roughness
      u10 = US%Z_to_L*USTair/sqrt(Cd)  ! Compute new u10 from derived Cd, while loop
                             ! ends and checks for convergence...CT counter
                             ! makes sure loop doesn't run away if function
                             ! doesn't converge.  This code was produced offline
                             ! and converged rapidly (e.g. 2 cycles)
                             ! for ustar=0.0001:0.0001:10.
      if (CT>20) then
        u10 = US%Z_to_L*USTair/sqrt(0.0015) ! I don't expect to get here, but just
                                !  in case it will output a reasonable value.
        exit
      endif
    enddo

  else ! Use more efficient expressions that are mathematically equivalent to those above.
    u10 = US%Z_to_L*USTair * sqrt(1000.0) ! First guess for u10.
    ! In the line above 1000 is the inverse of a plausible first guess of the drag coefficient.
    I_vonKar = 1.0 / CS%vonKar
    I_ten_m_scale = 0.1*US%Z_to_m

    do CT=1,20
      if (abs(u10a - u10) <= 0.001*u10) exit   ! Check for convergence.
      u10a = u10
      alpha = min(CS%Charnock_min, CS%Charnock_slope_U10 * u10 + CS%Charnock_intercept)
      z0rough = alpha * (CS%I_g_Earth * USTair**2) ! Compute z0rough from ustar guess
      z0 = z0sm + z0rough
      I_sqrtCd = abs(log(z0 * I_ten_m_scale)) * I_vonKar ! Compute Cd from derived roughness
      u10 = US%Z_to_L*USTair * I_sqrtCd  ! Compute new u10 from the derived Cd.
    enddo

    ! Output a reasonable estimate of u10 if the iteration has not converged. The hard-coded
    ! number 25.82 is 1/sqrt(0.0015) to 4 decimal places, but the exact value should not matter.
    if (abs(u10a - u10) > 0.001*u10) u10 = US%Z_to_L*USTair * 25.82
  endif

end subroutine ust_2_u10_coare3p5

!> Clear pointers, deallocate memory
subroutine Waves_end(CS)
  type(wave_parameters_CS), pointer :: CS !< Control structure

  if (allocated(CS%WaveNum_Cen)) deallocate( CS%WaveNum_Cen )
  if (allocated(CS%Freq_Cen))    deallocate( CS%Freq_Cen )
  if (allocated(CS%Us_x))        deallocate( CS%Us_x )
  if (allocated(CS%Us_y))        deallocate( CS%Us_y )
  if (allocated(CS%La_turb))     deallocate( CS%La_turb )
  if (allocated(CS%STKx0))       deallocate( CS%STKx0 )
  if (allocated(CS%STKy0))       deallocate( CS%STKy0 )
  if (allocated(CS%UStk_Hb))     deallocate( CS%UStk_Hb )
  if (allocated(CS%VStk_Hb))     deallocate( CS%VStk_Hb )
  if (allocated(CS%Omega_w2x))   deallocate( CS%Omega_w2x )
  if (allocated(CS%KvS))         deallocate( CS%KvS )
  if (allocated(CS%Us0_y))       deallocate( CS%Us0_y )
  if (allocated(CS%Us0_x))       deallocate( CS%Us0_x )

  deallocate( CS )

end subroutine Waves_end

!> Register wave restart fields. To be called before MOM_wave_interface_init
subroutine waves_register_restarts(CS, HI, GV, US, param_file, restart_CSp)
  type(wave_parameters_CS), pointer       :: CS           !< Wave parameter Control structure
  type(hor_index_type),     intent(inout) :: HI           !< Grid structure
  type(verticalGrid_type),  intent(in)    :: GV           !< Vertical grid structure
  type(unit_scale_type),    intent(in)    :: US           !< A dimensional unit scaling type
  type(param_file_type),    intent(in)    :: param_file   !< Input parameter structure
  type(MOM_restart_CS),     pointer       :: restart_CSp  !< Restart structure, data intent(inout)
  ! Local variables
  type(vardesc) :: vd(2)
  logical :: use_waves
  logical :: StatisticalWaves
  logical :: time_tendency_term
  character(len=40)  :: mdl = "MOM_wave_interface" !< This module's name.

  if (associated(CS)) then
    call MOM_error(FATAL, "waves_register_restarts: Called with initialized waves control structure")
  endif
  allocate(CS)

  call get_param(param_file, mdl, "USE_WAVES", use_waves, &
       "If true, enables surface wave modules.", do_not_log=.true., default=.false.)

  ! Check if using LA_LI2016
  call get_param(param_file,mdl,"USE_LA_LI2016",StatisticalWaves,     &
                 do_not_log=.true.,default=.false.)

  if (.not.(use_waves .or. StatisticalWaves)) return

  call get_param(param_file, mdl, "STOKES_DDT", time_tendency_term, do_not_log=.true., default=.false.)

  if (time_tendency_term) then
    ! Allocate wave fields needed for restart file
    allocate(CS%Us_x_prev(HI%isdB:HI%IedB,HI%jsd:HI%jed,GV%ke), source=0.0)
    allocate(CS%Us_y_prev(HI%isd:HI%Ied,HI%jsdB:HI%jedB,GV%ke), source=0.0)

    ! Register to restart files.  If these are not found in a restart file, they stay 0.
    vd(1) = var_desc("Us_x_prev", "m s-1", "3d zonal Stokes drift profile",&
                     hor_grid='u', z_grid='L')
    vd(2) = var_desc("Us_y_prev", "m s-1", "3d meridional Stokes drift profile",&
                      hor_grid='v', z_grid='L')
    call register_restart_pair(CS%US_x_prev, CS%US_y_prev, vd(1), vd(2), .false., &
                               restart_CSp, conversion=US%L_T_to_m_s)
  endif

end subroutine waves_register_restarts

!> \namespace  mom_wave_interface
!!
!! \author Brandon Reichl, 2018.
!!
!! This module should be moved as wave coupling progresses and
!! likely will should mirror the iceberg or sea-ice model set-up.
!!
!! This module is meant to contain the routines to read in and
!! interpret surface wave data for MOM6. In its original form, the
!! capabilities include setting the Stokes drift in the model (from a
!! variety of sources including prescribed, empirical, and input
!! files).  In short order, the plan is to also amend the subroutine
!! to accept Stokes drift information from an external coupler.
!! Eventually, it will be necessary to break this file apart so that
!! general wave information may be stored in the control structure
!! and the Stokes drift effect can be isolated from processes such as
!! sea-state dependent momentum fluxes, gas fluxes, and other wave
!! related air-sea interaction and boundary layer phenomenon.
!!
!! The Stokes drift are stored on the C-grid with the conventional
!! protocol to interpolate to the h-grid to compute Langmuir number,
!! the primary quantity needed for Langmuir turbulence
!! parameterizations in both the ePBL and KPP approach.  This module
!! also computes full 3d Stokes drift profiles, which will be useful
!! if second-order type boundary layer parameterizations are
!! implemented (perhaps via GOTM, work in progress).

end module MOM_wave_interface
