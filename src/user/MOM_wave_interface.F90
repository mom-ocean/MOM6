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
use MOM_safe_alloc,    only : safe_alloc_ptr
use MOM_time_manager,  only : time_type, operator(+), operator(/)
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs, surface
use MOM_verticalgrid,  only : verticalGrid_type

use netcdf, only : NF90_open, NF90_inq_varid, NF90_inquire_variable, NF90_get_var
use netcdf, only : NF90_inquire_dimension, NF90_close, NF90_NOWRITE, NF90_NOERR

implicit none ; private

#include <MOM_memory.h>

public MOM_wave_interface_init ! Public interface to fully initialize the wave routines.
public MOM_wave_interface_init_lite ! Public interface to quick initialize this module.
public Update_Surface_Waves ! Public interface to update wave information at the
                            ! coupler/driver level.
public Update_Stokes_Drift ! Public interface to update the Stokes drift profiles
                           ! called in step_mom.
public get_Langmuir_Number ! Public interface to compute Langmuir number called from
                           ! ePBL or KPP routines.
public StokesMixing ! NOT READY - Public interface to add down-Stokes gradient
                    ! momentum mixing (e.g. the approach of Harcourt 2013/2015)
public CoriolisStokes ! NOT READY - Public interface to add Coriolis-Stokes acceleration
                      ! of the mean currents, needed for comparison with LES.  It is
                      ! presently advised against implementing in non-1d settings without
                      ! serious consideration of the full 3d wave-averaged Navier-Stokes
                      ! CL2 effects.
public Waves_end ! public interface to deallocate and free wave related memory.

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Container for all surface wave related parameters
type, public :: wave_parameters_CS ; private

  !Main surface wave options
  logical, public :: UseWaves         !< Flag to enable surface gravity wave feature
  logical, public :: LagrangianMixing !< This feature is in development and not ready
                                      !! True if Stokes drift is present and mixing
                                      !! should be applied to Lagrangian current
                                      !! (mean current + Stokes drift).
                                      !! See Reichl et al., 2016 KPP-LT approach
  logical, public :: StokesMixing     !< This feature is in development and not ready.
                                      !! True if vertical mixing of momentum
                                      !! should be applied directly to Stokes current
                                      !! (with separate mixing parameter for Eulerian
                                      !! mixing contribution).
                                      !! See Harcourt 2013, 2015 Second-Moment approach
  logical, public :: CoriolisStokes   !< This feature is in development and not ready.
                                      ! True if Coriolis-Stokes acceleration should be applied.
  integer, public :: StkLevelMode=1   !< Sets if Stokes drift is defined at mid-points
                                      !! or layer averaged.  Set to 0 if mid-point and set to
                                      !! 1 if average value of Stokes drift over level.
                                      !! If advecting with Stokes transport, 1 is the correct
                                      !! approach.

  ! Surface Wave Dependent 1d/2d/3d vars
  integer, public :: NumBands =0 !< Number of wavenumber/frequency partitions to receive
                               !! This needs to match the number of bands provided
                               !! via either coupling or file.
  real, allocatable, dimension(:), public :: &
    WaveNum_Cen        !< Wavenumber bands for read/coupled [m-1]
  real, allocatable, dimension(:), public :: &
    Freq_Cen           !< Frequency bands for read/coupled [s-1]
  real, allocatable, dimension(:), public :: &
    PrescribedSurfStkX !< Surface Stokes drift if prescribed [m s-1]
  real, allocatable, dimension(:), public :: &
    PrescribedSurfStkY !< Surface Stokes drift if prescribed [m s-1]
  real, allocatable, dimension(:,:,:), public :: &
    Us_x               !< 3d zonal Stokes drift profile [m s-1]
                       !! Horizontal -> U points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
    Us_y               !< 3d meridional Stokes drift profile [m s-1]
                       !! Horizontal -> V points
                       !! Vertical -> Mid-points
  real, allocatable, dimension(:,:), public :: &
    La_SL,&            !< SL Langmuir number (directionality factored later)
                       !! Horizontal -> H points
    La_Turb            !< Aligned Turbulent Langmuir number
                       !! Horizontal -> H points
  real, allocatable, dimension(:,:), public :: &
    US0_x              !< Surface Stokes Drift (zonal, m/s)
                       !! Horizontal -> U points
  real, allocatable, dimension(:,:), public :: &
    US0_y              !< Surface Stokes Drift (meridional, m/s)
                       !! Horizontal -> V points
  real, allocatable, dimension(:,:,:), public :: &
    STKx0              !< Stokes Drift spectrum (zonal, m/s)
                       !! Horizontal -> U points
                       !! 3rd dimension -> Freq/Wavenumber
  real, allocatable, dimension(:,:,:), public :: &
    STKy0              !< Stokes Drift spectrum (meridional, m/s)
                       !! Horizontal -> V points
                       !! 3rd dimension -> Freq/Wavenumber
  real, allocatable, dimension(:,:,:), public :: &
    KvS                !< Viscosity for Stokes Drift shear [Z2 T-1 ~> m2 s-1]

  ! Pointers to auxiliary fields
  type(time_type), pointer, public :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer, public :: diag !< A structure that is used to regulate the
                                           !! timing of diagnostic output.

  !> An arbitrary lower-bound on the Langmuir number.  Run-time parameter.
  !! Langmuir number is sqrt(u_star/u_stokes). When both are small
  !! but u_star is orders of magnitude smaller the Langmuir number could
  !! have unintended consequences.  Since both are small it can be safely capped
  !! to avoid such consequences.
  real :: La_min = 0.05

  !>@{ Diagnostic handles
  integer, public :: id_surfacestokes_x = -1 , id_surfacestokes_y = -1
  integer, public :: id_3dstokes_x = -1 , id_3dstokes_y = -1
  integer, public :: id_La_turb = -1
  !>@}

end type wave_parameters_CS

! Options not needed outside of this module

integer :: WaveMethod=-99 !< Options for including wave information
                          !! Valid (tested) choices are:
                          !!   0 - Test Profile
                          !!   1 - Surface Stokes Drift Bands
                          !!   2 - DHH85
                          !!   3 - LF17
                          !! -99 - No waves computed, but empirical Langmuir number used.
                          !! \todo Module variable! Move into a control structure.

! Options if WaveMethod is Surface Stokes Drift Bands (1)
integer, public :: PartitionMode !< Method for partition mode (meant to check input)
                                 !! 0 - wavenumbers
                                 !! 1 - frequencies
                                 !! \todo Module variable! Move into a control structure.
integer :: DataSource !< Integer that specifies where the Model Looks for Data
                      !! Valid choices are:
                      !! 1 - FMS DataOverride Routine
                      !! 2 - Reserved For Coupler
                      !! 3 - User input (fixed values, useful for 1d testing)
                      !! \todo Module variable! Move into a control structure.

! Options if using FMS DataOverride Routine
character(len=40)  :: SurfBandFileName !< Filename if using DataOverride
                                       !! \todo Module variable! Move into a control structure.
logical :: dataoverrideisinitialized !< Flag for DataOverride Initialization
                                     !! \todo Module variable! Move into a control structure.

! Options for computing Langmuir number
real :: LA_FracHBL         !< Fraction of OSBL for averaging Langmuir number
                           !! \todo Module variable! Move into a control structure.
logical :: LA_Misalignment = .false. !< Flag to use misalignment in Langmuir number
                                     !! \todo Module variable! Move into a control structure.

! This include declares and sets the variable "version".
#include "version_variable.h"

character(len=40)  :: mdl = "MOM_wave_interface" !< This module's name.

!>@{ Undocumented parameters.
!! \todo These module variables need to be documented as static/private variables or moved
!! into a control structure.
! Switches needed in import_stokes_drift
integer, parameter :: TESTPROF = 0, SURFBANDS = 1, &
                      DHH85 = 2, LF17 = 3, NULL_WaveMethod=-99, &
                      DATAOVR = 1, COUPLER = 2, INPUT = 3

! Options For Test Prof
Real    :: TP_STKX0, TP_STKY0, TP_WVL
logical :: WaveAgePeakFreq ! Flag to use W
logical :: StaticWaves, DHH85_Is_Set
real    :: WaveAge, WaveWind
real    :: PI
!>@}

contains

!> Initializes parameters related to MOM_wave_interface
subroutine MOM_wave_interface_init(time, G, GV, US, param_file, CS, diag )
  type(time_type), target, intent(in)    :: Time       !< Model time
  type(ocean_grid_type),   intent(inout) :: G          !< Grid structure
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Input parameter structure
  type(wave_parameters_CS), pointer      :: CS         !< Wave parameter control structure
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostic Pointer
  ! Local variables
  ! I/O
  character*(13) :: TMPSTRING1,TMPSTRING2
  character*(5), parameter  :: NULL_STRING      = "EMPTY"
  character*(12), parameter :: TESTPROF_STRING  = "TEST_PROFILE"
  character*(13), parameter :: SURFBANDS_STRING = "SURFACE_BANDS"
  character*(5), parameter  :: DHH85_STRING     = "DHH85"
  character*(4), parameter  :: LF17_STRING      = "LF17"
  character*(12), parameter :: DATAOVR_STRING   = "DATAOVERRIDE"
  character*(7), parameter  :: COUPLER_STRING   = "COUPLER"
  character*(5), parameter  :: INPUT_STRING     = "INPUT"

  ! Dummy Check
  if (associated(CS)) then
    call MOM_error(FATAL, "wave_interface_init called with an associated control structure.")
    return
  endif

  PI=4.0*atan(1.0)

  ! Allocate CS and set pointers
  allocate(CS)

  CS%diag => diag
  CS%Time => Time

  ! Add any initializations needed here
  dataOverrideIsInitialized = .false.

  ! The only way to get here is with UseWaves enabled.
  CS%UseWaves=.true.

  call log_version(param_file, mdl, version)

  ! Wave modified physics
  !  Presently these are all in research mode
  call get_param(param_file, mdl, "LAGRANGIAN_MIXING", CS%LagrangianMixing, &
       "Flag to use Lagrangian Mixing of momentum", units="", &
       Default=.false.)
  if (CS%LagrangianMixing) then
    ! Force Code Intervention
    call MOM_error(FATAL,"Should you be enabling Lagrangian Mixing? Code not ready.")
  endif
  call get_param(param_file, mdl, "STOKES_MIXING", CS%StokesMixing, &
       "Flag to use Stokes Mixing of momentum", units="", &
       Default=.false.)
  if (CS%StokesMixing) then
    ! Force Code Intervention
    call MOM_error(FATAL,"Should you be enabling Stokes Mixing? Code not ready.")
  endif
  call get_param(param_file, mdl, "CORIOLIS_STOKES", CS%CoriolisStokes, &
       "Flag to use Coriolis Stokes acceleration", units="", &
       Default=.false.)
  if (CS%CoriolisStokes) then
    ! Force Code Intervention
    call MOM_error(FATAL,"Should you be enabling Coriolis-Stokes? Code not ready.")
  endif

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
       "                 speed following Li and Fox-Kemper 2017.\n",  &
       units='', default=NULL_STRING)
  select case (TRIM(TMPSTRING1))
  case (NULL_STRING)! No Waves
    call MOM_error(FATAL, "wave_interface_init called with no specified "//&
                           "WAVE_METHOD.")
  case (TESTPROF_STRING)! Test Profile
    WaveMethod = TESTPROF
    call get_param(param_file,mdl,"TP_STKX_SURF",TP_STKX0,&
         'Surface Stokes (x) for test profile',&
         units='m/s',default=0.1)
    call get_param(param_file,mdl,"TP_STKY_SURF",TP_STKY0,&
         'Surface Stokes (y) for test profile',&
         units='m/s',default=0.0)
    call get_param(param_file,mdl,"TP_WVL",TP_WVL,&
         units='m', default=50.0, scale=US%m_to_Z)
  case (SURFBANDS_STRING)! Surface Stokes Drift Bands
    WaveMethod = SURFBANDS
    call get_param(param_file, mdl, "SURFBAND_SOURCE",TMPSTRING2,       &
       "Choice of SURFACE_BANDS data mode, valid options include: \n"// &
       " DATAOVERRIDE  - Read from NetCDF using FMS DataOverride. \n"// &
       " COUPLER       - Look for variables from coupler pass \n"//     &
       " INPUT         - Testing with fixed values.",                   &
       units='', default=NULL_STRING)
    select case (TRIM(TMPSTRING2))
    case (NULL_STRING)! Default
      call MOM_error(FATAL, "wave_interface_init called with SURFACE_BANDS"//&
                           " but no SURFBAND_SOURCE.")
    case (DATAOVR_STRING)! Using Data Override
      DataSource = DATAOVR
      call get_param(param_file, mdl, "SURFBAND_FILENAME", SurfBandFileName, &
           "Filename of surface Stokes drift input band data.", default="StkSpec.nc")
    case (COUPLER_STRING)! Reserved for coupling
      DataSource = Coupler
      ! This is just to make something work, but it needs to be read from the wavemodel.
      call get_param(param_file,mdl,"STK_BAND_COUPLER",CS%NumBands,                &
         "STK_BAND_COUPLER is the number of Stokes drift bands in the coupler. "// &
         "This has to be consistent with the number of Stokes drift bands in WW3, "//&
         "or the model will fail.",units='', default=1)
      allocate( CS%WaveNum_Cen(CS%NumBands) )
      allocate( CS%STKx0(G%isdB:G%iedB,G%jsd:G%jed,CS%NumBands))
      allocate( CS%STKy0(G%isdB:G%iedB,G%jsd:G%jed,CS%NumBands))
      CS%WaveNum_Cen(:) = 0.0
      CS%STKx0(:,:,:) = 0.0
      CS%STKy0(:,:,:) = 0.0
      partitionmode = 0
    case (INPUT_STRING)! A method to input the Stokes band (globally uniform)
      DataSource = Input
      call get_param(param_file,mdl,"SURFBAND_NB",CS%NumBands,              &
         "Prescribe number of wavenumber bands for Stokes drift. "//      &
         "Make sure this is consistnet w/ WAVENUMBERS, STOKES_X, and "// &
         "STOKES_Y, there are no safety checks in the code.",              &
         units='', default=1)
      allocate( CS%WaveNum_Cen(1:CS%NumBands) )
      CS%WaveNum_Cen(:) = 0.0
      allocate( CS%PrescribedSurfStkX(1:CS%NumBands))
      CS%PrescribedSurfStkX(:) = 0.0
      allocate( CS%PrescribedSurfStkY(1:CS%NumBands))
      CS%PrescribedSurfStkY(:) = 0.0
      allocate( CS%STKx0(G%isdB:G%iedB,G%jsd:G%jed,1:CS%NumBands))
      CS%STKx0(:,:,:) = 0.0
      allocate( CS%STKy0(G%isd:G%ied,G%jsdB:G%jedB,1:CS%NumBands))
      CS%STKy0(:,:,:) = 0.0
      partitionmode=0
      call get_param(param_file,mdl,"SURFBAND_WAVENUMBERS",CS%WaveNum_Cen,      &
           "Central wavenumbers for surface Stokes drift bands.",units='rad/m', &
           default=0.12566)
      call get_param(param_file,mdl,"SURFBAND_STOKES_X",CS%PrescribedSurfStkX,      &
           "X-direction surface Stokes drift for bands.",units='m/s', &
           default=0.15)
      call get_param(param_file,mdl,"SURFBAND_STOKES_Y",CS%PrescribedSurfStkY,      &
           "Y-direction surface Stokes drift for bands.",units='m/s', &
           default=0.0)
    case default! No method provided
      call MOM_error(FATAL,'Check WAVE_METHOD.')
    end select

  case (DHH85_STRING)!Donelan et al., 1985 spectrum
    WaveMethod = DHH85
    call MOM_error(WARNING,"DHH85 only ever set-up for uniform cases w/"//&
                           " Stokes drift in x-direction.")
    call get_param(param_file,mdl,"DHH85_AGE_FP",WaveAgePeakFreq,   &
         "Choose true to use waveage in peak frequency.", &
         units='', default=.false.)
    call get_param(param_file,mdl,"DHH85_AGE",WaveAge,   &
         "Wave Age for DHH85 spectrum.", &
         units='', default=1.2)
    call get_param(param_file,mdl,"DHH85_WIND",WaveWind,   &
         "Wind speed for DHH85 spectrum.", &
         units='', default=10.0)
    call get_param(param_file,mdl,"STATIC_DHH85",StaticWaves,   &
         "Flag to disable updating DHH85 Stokes drift.", &
          default=.false.)
  case (LF17_STRING)!Li and Fox-Kemper 17 wind-sea Langmuir number
    WaveMethod = LF17
  case default
    call MOM_error(FATAL,'Check WAVE_METHOD.')
  end select

  ! Langmuir number Options
  call get_param(param_file, mdl, "LA_DEPTH_RATIO", LA_FracHBL,              &
         "The depth (normalized by BLD) to average Stokes drift over in "//&
         "Langmuir number calculation, where La = sqrt(ust/Stokes).",       &
         units="nondim",default=0.04)
  call get_param(param_file, mdl, "LA_MISALIGNMENT", LA_Misalignment,    &
         "Flag (logical) if using misalignment bt shear and waves in LA",&
         default=.false.)
  call get_param(param_file, mdl, "MIN_LANGMUIR", CS%La_min,    &
         "A minimum value for all Langmuir numbers that is not physical, "//&
         "but is likely only encountered when the wind is very small and "//&
         "therefore its effects should be mostly benign.",units="nondim",&
         default=0.05)

  ! Allocate and initialize
  ! a. Stokes driftProfiles
  allocate(CS%Us_x(G%isdB:G%IedB,G%jsd:G%jed,GV%ke))
  CS%Us_x(:,:,:) = 0.0
  allocate(CS%Us_y(G%isd:G%Ied,G%jsdB:G%jedB,GV%ke))
  CS%Us_y(:,:,:) = 0.0
  ! b. Surface Values
  allocate(CS%US0_x(G%isdB:G%iedB,G%jsd:G%jed))
  CS%US0_x(:,:) = 0.0
  allocate(CS%US0_y(G%isd:G%ied,G%jsdB:G%jedB))
  CS%US0_y(:,:) = 0.0
  ! c. Langmuir number
  allocate(CS%La_SL(G%isc:G%iec,G%jsc:G%jec))
  allocate(CS%La_turb(G%isc:G%iec,G%jsc:G%jec))
  CS%La_SL(:,:) = 0.0
  CS%La_turb (:,:) = 0.0
  ! d. Viscosity for Stokes drift
  if (CS%StokesMixing) then
    allocate(CS%KvS(G%isd:G%Ied,G%jsd:G%jed,GV%ke))
    CS%KvS(:,:,:) = 0.0
  endif

  ! Initialize Wave related outputs
  CS%id_surfacestokes_y = register_diag_field('ocean_model','surface_stokes_y', &
       CS%diag%axesCu1,Time,'Surface Stokes drift (y)','m s-1')
  CS%id_surfacestokes_x = register_diag_field('ocean_model','surface_stokes_x', &
       CS%diag%axesCv1,Time,'Surface Stokes drift (x)','m s-1')
  CS%id_3dstokes_y = register_diag_field('ocean_model','3d_stokes_y', &
       CS%diag%axesCvL,Time,'3d Stokes drift (y)','m s-1')
  CS%id_3dstokes_x = register_diag_field('ocean_model','3d_stokes_x', &
       CS%diag%axesCuL,Time,'3d Stokes drift (y)','m s-1')
  CS%id_La_turb = register_diag_field('ocean_model','La_turbulent',&
       CS%diag%axesT1,Time,'Surface (turbulent) Langmuir number','nondim')

  return
end subroutine MOM_wave_interface_init

!> A 'lite' init subroutine to initialize a few inputs needed if using wave information
!! with the wind-speed dependent Stokes drift formulation of LF17
subroutine MOM_wave_interface_init_lite(param_file)
  type(param_file_type), intent(in) :: param_file !< Input parameter structure
  character*(5), parameter  :: NULL_STRING      = "EMPTY"
  character*(4), parameter  :: LF17_STRING      = "LF17"
  character*(13) :: TMPSTRING1
  logical :: StatisticalWaves

  ! Langmuir number Options
  call get_param(param_file, mdl, "LA_DEPTH_RATIO", LA_FracHBL,              &
       "The depth (normalized by BLD) to average Stokes drift over in "//&
       "Langmuir number calculation, where La = sqrt(ust/Stokes).",       &
       units="nondim",default=0.04)

  ! Check if using LA_LI2016
  call get_param(param_file,mdl,"USE_LA_LI2016",StatisticalWaves,     &
                 do_not_log=.true.,default=.false.)
  if (StatisticalWaves) then
    WaveMethod = LF17
    PI=4.0*atan(1.0)
  else
    WaveMethod = NULL_WaveMethod
  end if

  return
end subroutine MOM_wave_interface_init_lite

!> Subroutine that handles updating of surface wave/Stokes drift related properties
subroutine Update_Surface_Waves(G, GV, US, Day, dt, CS, forces)
  type(wave_parameters_CS), pointer    :: CS  !< Wave parameter Control structure
  type(ocean_grid_type), intent(inout) :: G   !< Grid structure
  type(verticalGrid_type), intent(in)  :: GV  !< Vertical grid structure
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  type(time_type),         intent(in)  :: Day !< Current model time
  type(time_type),         intent(in)  :: dt  !< Timestep as a time-type
  type(mech_forcing),      intent(in), optional  :: forces !< MOM_forcing_type
  ! Local variables
  integer :: ii, jj, kk, b
  type(time_type) :: Day_Center

  ! Computing central time of time step
  Day_Center = Day + DT/2

  if (WaveMethod == TESTPROF) then
    ! Do nothing
  elseif (WaveMethod==SURFBANDS) then
    if (DataSource==DATAOVR) then
      call Surface_Bands_by_data_override(day_center, G, GV, US, CS)
    elseif (DataSource==Coupler) then
      if (.not.present(FORCES)) then
        call MOM_error(FATAL,"The option SURFBAND = COUPLER can not be used with "//&
             "this driver. If you are using a coupled driver with a wave model then "//&
             "check the arguments in the subroutine call to Update_Surface_Waves, "//&
             "otherwise select another option for SURFBAND_SOURCE.")
      endif
      if (size(CS%WaveNum_Cen).ne.size(forces%stk_wavenumbers)) then
        call MOM_error(FATAL, "Number of wavenumber bands in WW3 does not match that in MOM6. "//&
             "Make sure that STK_BAND_COUPLER in MOM6 input is equal to the number of bands in "//&
             "ww3_grid.inp, and that your mod_def.ww3 is up to date.")
      endif

      do b=1,CS%NumBands
        CS%WaveNum_Cen(b) = forces%stk_wavenumbers(b)
        !Interpolate from a grid to c grid
        do jj=G%jsc,G%jec
          do II=G%iscB,G%iecB
            CS%STKx0(II,jj,b) = 0.5*(forces%UStkb(ii,jj,b)+forces%UStkb(ii+1,jj,b))
          enddo
        enddo
        do JJ=G%jscB, G%jecB
          do ii=G%isc,G%iec
            CS%STKY0(ii,JJ,b) = 0.5*(forces%VStkb(ii,jj,b)+forces%VStkb(ii,jj+1,b))
          enddo
        enddo
        call pass_vector(CS%STKx0(:,:,b),CS%STKy0(:,:,b), G%Domain)
      enddo
    elseif (DataSource==Input) then
      do b=1,CS%NumBands
        do jj=G%jsd,G%jed
          do II=G%isdB,G%iedB
            CS%STKx0(II,jj,b) = CS%PrescribedSurfStkX(b)
          enddo
        enddo
        do JJ=G%jsdB, G%jedB
          do ii=G%isd,G%ied
            CS%STKY0(ii,JJ,b) = CS%PrescribedSurfStkY(b)
          enddo
        enddo
      enddo
    endif
  endif

  return
end subroutine Update_Surface_Waves

!> Constructs the Stokes Drift profile on the model grid based on
!! desired coupling options
subroutine Update_Stokes_Drift(G, GV, US, CS, h, ustar)
  type(wave_parameters_CS), pointer       :: CS    !< Wave parameter Control structure
  type(ocean_grid_type),    intent(inout) :: G     !< Grid structure
  type(verticalGrid_type),  intent(in)    :: GV    !< Vertical grid structure
  type(unit_scale_type),    intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h     !< Thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), &
                            intent(in)    :: ustar !< Wind friction velocity [Z T-1 ~> m s-1].
  ! Local Variables
  real    :: Top, MidPoint, Bottom, one_cm, level_thick, min_level_thick_avg
  real    :: DecayScale
  real    :: CMN_FAC, WN, UStokes
  real    :: La
  integer :: ii, jj, kk, b, iim1, jjm1

  one_cm = 0.01*US%m_to_Z
  min_level_thick_avg = 1.e-3*US%m_to_Z

  ! 1. If Test Profile Option is chosen
  !    Computing mid-point value from surface value and decay wavelength
  if (WaveMethod==TESTPROF) then
    DecayScale = 4.*PI / TP_WVL !4pi
    do jj = G%jsd,G%jed
      do II = G%isdB,G%iedB
        IIm1 = max(1,II-1)
        Bottom = 0.0
        MidPoint = 0.0
        do kk = 1,GV%ke
          Top = Bottom
          MidPoint = Bottom - GV%H_to_Z*0.25*(h(II,jj,kk)+h(IIm1,jj,kk))
          Bottom = Bottom - GV%H_to_Z*0.5*(h(II,jj,kk)+h(IIm1,jj,kk))
          CS%Us_x(II,jj,kk) = TP_STKX0*exp(MidPoint*DecayScale)
        enddo
      enddo
    enddo
    do JJ = G%jsdB,G%jedB
      do ii = G%isd,G%ied
        JJm1 = max(1,JJ-1)
        Bottom = 0.0
        MidPoint = 0.0
        do kk = 1,GV%ke
          Top = Bottom
          MidPoint = Bottom - GV%H_to_Z*0.25*(h(ii,JJ,kk)+h(ii,JJm1,kk))
          Bottom = Bottom - GV%H_to_Z*0.5*(h(ii,JJ,kk)+h(ii,JJm1,kk))
          CS%Us_y(ii,JJ,kk) = TP_STKY0*exp(MidPoint*DecayScale)
        enddo
      enddo
    enddo
  ! 2. If Surface Bands is chosen
  !    In wavenumber mode compute integral for layer averaged Stokes drift.
  !    In frequency mode compuate value at midpoint.
  elseif (WaveMethod==SURFBANDS) then
    CS%Us_x(:,:,:) = 0.0
    CS%Us_y(:,:,:) = 0.0
    CS%Us0_x(:,:) = 0.0
    CS%Us0_y(:,:) = 0.0
    ! Computing X direction Stokes drift
    do jj = G%jsd,G%jed
      do II = G%isdB,G%iedB
        ! 1. First compute the surface Stokes drift
        !    by integrating over the partitionas.
        do b = 1,CS%NumBands
          if (PartitionMode==0) then
            ! In wavenumber we are averaging over (small) level
            CMN_FAC = (1.0-exp(-one_cm*2.*CS%WaveNum_Cen(b))) / &
                      (one_cm*2.*CS%WaveNum_Cen(b))
          elseif (PartitionMode==1) then
             ! In frequency we are not averaging over level and taking top
            CMN_FAC = 1.0
          endif
          CS%US0_x(II,jj) = CS%US0_x(II,jj) + CS%STKx0(II,jj,b)*CMN_FAC
        enddo
        ! 2. Second compute the level averaged Stokes drift
        bottom = 0.0
        do kk = 1,GV%ke
          Top = Bottom
          IIm1 = max(II-1,1)
          level_thick = 0.5*GV%H_to_Z*(h(II,jj,kk)+h(IIm1,jj,kk))
          MidPoint = Bottom - 0.5*level_thick
          Bottom = Bottom - level_thick
          ! -> Stokes drift in thin layers not averaged.
          if (level_thick>min_level_thick_avg) then
            do b = 1,CS%NumBands
              if (PartitionMode==0) then
              ! In wavenumber we are averaging over level
                CMN_FAC = (exp(Top*2.*CS%WaveNum_Cen(b))-exp(Bottom*2.*CS%WaveNum_Cen(b)))&
                          / ((Top-Bottom)*(2.*CS%WaveNum_Cen(b)))
              elseif (PartitionMode==1) then
                if (CS%StkLevelMode==0) then
                  ! Take the value at the midpoint
                  CMN_FAC = exp(MidPoint*2.*(2.*PI*CS%Freq_Cen(b)*US%T_to_s)**2/(US%L_to_Z**2*GV%g_Earth))
                elseif (CS%StkLevelMode==1) then
                  ! Use a numerical integration and then
                  ! divide by layer thickness
                  WN = (2.*PI*CS%Freq_Cen(b)*US%T_to_s)**2 / (US%L_to_Z**2*GV%g_Earth) !bgr bug-fix missing g
                  CMN_FAC = (exp(2.*WN*Top)-exp(2.*WN*Bottom)) / (2.*WN*(Top-Bottom))
                endif
              endif
              CS%US_x(II,jj,kk) = CS%US_x(II,jj,kk) + CS%STKx0(II,jj,b)*CMN_FAC
            enddo
          else
            ! Take the value at the midpoint
            do b = 1,CS%NumBands
              if (PartitionMode==0) then
                CMN_FAC = exp(MidPoint*2.*CS%WaveNum_Cen(b))
              elseif (PartitionMode==1) then
                 CMN_FAC = exp(MidPoint*2.*(2.*PI*CS%Freq_Cen(b)*US%T_to_s)**2/(US%L_to_Z**2*GV%g_Earth))
              endif
              CS%US_x(II,jj,kk) = CS%US_x(II,jj,kk) + CS%STKx0(II,jj,b)*CMN_FAC
            enddo
          endif
        enddo
      enddo
    enddo
    ! Computing Y direction Stokes drift
    do JJ = G%jsdB,G%jedB
      do ii = G%isd,G%ied
        ! Compute the surface values.
        do b = 1,CS%NumBands
          if (PartitionMode==0) then
            ! In wavenumber we are averaging over (small) level
            CMN_FAC = (1.0-exp(-one_cm*2.*CS%WaveNum_Cen(b))) / &
                      (one_cm*2.*CS%WaveNum_Cen(b))
          elseif (PartitionMode==1) then
            ! In frequency we are not averaging over level and taking top
            CMN_FAC = 1.0
          endif
          CS%US0_y(ii,JJ) = CS%US0_y(ii,JJ) + CS%STKy0(ii,JJ,b)*CMN_FAC
        enddo
        ! Compute the level averages.
        bottom = 0.0
        do kk = 1,GV%ke
          Top = Bottom
          JJm1 = max(JJ-1,1)
          level_thick = 0.5*GV%H_to_Z*(h(ii,JJ,kk)+h(ii,JJm1,kk))
          MidPoint = Bottom - 0.5*level_thick
          Bottom = Bottom - level_thick
          ! -> Stokes drift in thin layers not averaged.
          if (level_thick>min_level_thick_avg) then
            do b = 1,CS%NumBands
              if (PartitionMode==0) then
              ! In wavenumber we are averaging over level
                CMN_FAC = (exp(Top*2.*CS%WaveNum_Cen(b))-exp(Bottom*2.*CS%WaveNum_Cen(b)))&
                          / ((Top-Bottom)*(2.*CS%WaveNum_Cen(b)))
              elseif (PartitionMode==1) then
                if (CS%StkLevelMode==0) then
                  ! Take the value at the midpoint
                  CMN_FAC = exp(MidPoint*2.*(2.*PI*CS%Freq_Cen(b)*US%T_to_s)**2/(US%L_to_Z**2*GV%g_Earth))
                elseif (CS%StkLevelMode==1) then
                  ! Use a numerical integration and then
                  ! divide by layer thickness
                  WN = (2.*PI*CS%Freq_Cen(b)*US%T_to_s)**2 / (US%L_to_Z**2*GV%g_Earth) !bgr bug-fix missing g
                  CMN_FAC = (exp(2.*WN*Top)-exp(2.*WN*Bottom)) / (2.*WN*(Top-Bottom))
                endif
              endif
              CS%US_y(ii,JJ,kk) = CS%US_y(ii,JJ,kk) + CS%STKy0(ii,JJ,b)*CMN_FAC
            enddo
          else
            ! Take the value at the midpoint
            do b = 1,CS%NumBands
              if (PartitionMode==0) then
                CMN_FAC = exp(MidPoint*2.*CS%WaveNum_Cen(b))
              elseif (PartitionMode==1) then
                 CMN_FAC = exp(MidPoint*2.*(2.*PI*CS%Freq_Cen(b)*US%T_to_s)**2/(US%L_to_Z**2*GV%g_Earth))
              endif
              CS%US_y(ii,JJ,kk) = CS%US_y(ii,JJ,kk) + CS%STKy0(ii,JJ,b)*CMN_FAC
            enddo
          endif
        enddo
      enddo
    enddo
  elseif (WaveMethod==DHH85) then
    if (.not.(StaticWaves .and. DHH85_is_set)) then
      do jj = G%jsd,G%jed
        do II = G%isdB,G%iedB
          bottom = 0.0
          do kk = 1,GV%ke
            Top = Bottom
            IIm1 = max(II-1,1)
            MidPoint = Bottom - GV%H_to_Z*0.25*(h(II,jj,kk)+h(IIm1,jj,kk))
            Bottom = Bottom - GV%H_to_Z*0.5*(h(II,jj,kk)+h(IIm1,jj,kk))
            !bgr note that this is using a u-point ii on h-point ustar
            !    this code has only been previous used for uniform
            !    grid cases.  This needs fixed if DHH85 is used for non
            !    uniform cases.
            call DHH85_mid(GV, US, MidPoint, UStokes)
            ! Putting into x-direction (no option for direction
            CS%US_x(II,jj,kk) = UStokes
          enddo
        enddo
      enddo
      do JJ = G%jsdB,G%jedB
        do ii = G%isd,G%ied
          Bottom = 0.0
          do kk=1, GV%ke
            Top = Bottom
            JJm1 = max(JJ-1,1)
            MidPoint = Bottom - GV%H_to_Z*0.25*(h(ii,JJ,kk)+h(ii,JJm1,kk))
            Bottom = Bottom - GV%H_to_Z*0.5*(h(ii,JJ,kk)+h(ii,JJm1,kk))
            !bgr note that this is using a v-point jj on h-point ustar
            !    this code has only been previous used for uniform
            !    grid cases.  This needs fixed if DHH85 is used for non
            !    uniform cases.
            ! call DHH85_mid(GV, US, Midpoint, UStokes)
            ! Putting into x-direction, so setting y direction to 0
            CS%US_y(ii,JJ,kk) = 0.0
            ! For rotational symmetry there should be the option for this to become = UStokes
            !    bgr - see note above, but this is true
            !          if this is used for anything
            !          other than simple LES comparison
          enddo
        enddo
      enddo
      DHH85_is_set = .true.
    endif
  else! Keep this else, fallback to 0 Stokes drift
    do kk= 1,GV%ke
      do jj = G%jsd,G%jed
        do II = G%isdB,G%iedB
          CS%Us_x(II,jj,kk) = 0.
        enddo
      enddo
      do JJ = G%jsdB,G%jedB
        do ii = G%isd,G%ied
          CS%Us_y(ii,JJ,kk) = 0.
        enddo
      enddo
    enddo
  endif

  ! Turbulent Langmuir number is computed here and available to use anywhere.
  ! SL Langmuir number requires mixing layer depth, and therefore is computed
  ! in the routine it is needed by (e.g. KPP or ePBL).
  do jj = G%jsc, G%jec
    do ii = G%isc,G%iec
      Top = h(ii,jj,1)*GV%H_to_Z
      call get_Langmuir_Number( La, G, GV, US, Top, ustar(ii,jj), ii, jj, &
             H(ii,jj,:),Override_MA=.false.,WAVES=CS)
      CS%La_turb(ii,jj) = La
    enddo
  enddo

  ! Output any desired quantities
  if (CS%id_surfacestokes_y>0) &
    call post_data(CS%id_surfacestokes_y, CS%us0_y, CS%diag)
  if (CS%id_surfacestokes_x>0) &
    call post_data(CS%id_surfacestokes_x, CS%us0_x, CS%diag)
  if (CS%id_3dstokes_y>0) &
    call post_data(CS%id_3dstokes_y, CS%us_y, CS%diag)
  if (CS%id_3dstokes_x>0) &
    call post_data(CS%id_3dstokes_x, CS%us_x, CS%diag)
  if (CS%id_La_turb>0) &
    call post_data(CS%id_La_turb, CS%La_turb, CS%diag)

end subroutine Update_Stokes_Drift

!> A subroutine to fill the Stokes drift from a NetCDF file
!! using the data_override procedures.
subroutine Surface_Bands_by_data_override(day_center, G, GV, US, CS)
  type(time_type),          intent(in) :: day_center !< Center of timestep
  type(wave_parameters_CS), pointer    :: CS         !< Wave structure
  type(ocean_grid_type), intent(inout) :: G          !< Grid structure
  type(verticalGrid_type),  intent(in) :: GV         !< Vertical grid structure
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type

  ! Local variables
  real    :: temp_x(SZI_(G),SZJ_(G)) ! Pseudo-zonal Stokes drift of band at h-points [m s-1]
  real    :: temp_y(SZI_(G),SZJ_(G)) ! Psuedo-meridional Stokes drift of band at h-points [m s-1]
  real    :: Top, MidPoint
  integer :: b
  integer :: i, j
  integer, dimension(4) :: start, counter, dims, dim_id
  character(len=12)  :: dim_name(4)
  character(20) :: varname, varread1, varread2
  integer :: rcode_fr, rcode_wn, ncid, varid_fr, varid_wn, id, ndims

  if (.not.dataOverrideIsInitialized) then
    call data_override_init(G%Domain)
    dataOverrideIsInitialized = .true.

    ! Read in number of wavenumber bands in file to set number to be read in
    ! Hardcoded filename/variables
    varread1 = 'wavenumber' !Old method gives wavenumber
    varread2 = 'frequency'  !New method gives frequency
    rcode_wn = NF90_OPEN(trim(SurfBandFileName), NF90_NOWRITE, ncid)
    if (rcode_wn /= 0) then
      call MOM_error(FATAL,"error opening file "//trim(SurfBandFileName)//&
            " in MOM_wave_interface.")
    endif

    ! Check if rcode_wn or rcode_fr is 0 (checks if input has wavenumber or frequency)
    rcode_wn = NF90_INQ_VARID(ncid, varread1, varid_wn)
    rcode_fr = NF90_INQ_VARID(ncid, varread2, varid_fr)

    if (rcode_wn /= 0 .and. rcode_fr /= 0) then
      call MOM_error(FATAL,"error finding variable "//trim(varread1)//&
         " or "//trim(varread2)//" in file "//trim(SurfBandFileName)//" in MOM_wave_interface.")

    elseif (rcode_wn == 0) then
      ! wavenumbers found:
      PartitionMode = 0
      rcode_wn = NF90_INQUIRE_VARIABLE(ncid, varid_wn, ndims=ndims, &
           dimids=dims)
      if (rcode_wn /= 0) then
        call MOM_error(FATAL, &
             'error inquiring dimensions MOM_wave_interface.')
      endif
      rcode_wn = NF90_INQUIRE_DIMENSION(ncid, dims(1), dim_name(1), len=id)
      if (rcode_wn /= 0) then
        call MOM_error(FATAL,"error reading dimension 1 data for "// &
             trim(varread1)//" in file "// trim(SurfBandFileName)//          &
             " in MOM_wave_interface.")
      endif
      rcode_wn = NF90_INQ_VARID(ncid, dim_name(1), dim_id(1))
      if (rcode_wn /= 0) then
        call MOM_error(FATAL,"error finding variable "//trim(dim_name(1))//&
          " in file "//trim(SurfBandFileName)//" in MOM_wave_interace.")
      endif
      ! Allocating size of wavenumber bins
      allocate( CS%WaveNum_Cen(1:id) )
      CS%WaveNum_Cen(:) = 0.0
    elseif (rcode_fr == 0) then
      ! frequencies found:
      PartitionMode = 1
      rcode_fr = NF90_INQUIRE_VARIABLE(ncid, varid_fr, ndims=ndims, &
           dimids=dims)
      if (rcode_fr /= 0) then
        call MOM_error(FATAL,&
             'error inquiring dimensions MOM_wave_interface.')
      endif
      rcode_fr = NF90_INQUIRE_DIMENSION(ncid, dims(1), dim_name(1), len=id)
      if (rcode_fr /= 0) then
        call MOM_error(FATAL,"error reading dimension 1 data for "// &
             trim(varread2)//" in file "// trim(SurfBandFileName)// &
             " in MOM_wave_interface.")
      endif
      rcode_fr = NF90_INQ_VARID(ncid, dim_name(1), dim_id(1))
      if (rcode_fr /= 0) then
        call MOM_error(FATAL,"error finding variable "//trim(dim_name(1))//&
             " in file "//trim(SurfBandFileName)//" in MOM_wave_interace.")
      endif
      ! Allocating size of frequency bins
      allocate( CS%Freq_Cen(1:id) )
      CS%Freq_Cen(:) = 0.0
      ! Allocating size of wavenumber bins
      allocate( CS%WaveNum_Cen(1:id) )
      CS%WaveNum_Cen(:) = 0.0
      allocate( CS%STKx0(G%isdB:G%iedB,G%jsd:G%jed,1:id))
      CS%STKx0(:,:,:) = 0.0
      allocate( CS%STKy0(G%isd:G%ied,G%jsdB:G%jedB,1:id))
      CS%STKy0(:,:,:) = 0.0
    endif

    ! Reading wavenumber bins/Frequencies
    start(:) = 1 ! Set all start to 1
    counter(:) = 1 ! Set all counter to 1
    counter(1) = id ! Set counter(1) to id (number of frequency bins)
    if (PartitionMode==0) then
      rcode_wn = NF90_GET_VAR(ncid, dim_id(1), CS%WaveNum_Cen, start, counter)
      if (rcode_wn /= 0) then
        call MOM_error(FATAL,&
             "error reading dimension 1 values for var_name "// &
             trim(varread1)//",dim_name "//trim(dim_name(1))//  &
             " in file "// trim(SurfBandFileName)//" in MOM_wave_interface")
      endif
      CS%NUMBANDS = ID
      do B = 1,CS%NumBands ; CS%WaveNum_Cen(b) = US%Z_to_m*CS%WaveNum_Cen(b) ; enddo
    elseif (PartitionMode==1) then
      rcode_fr = NF90_GET_VAR(ncid, dim_id(1), CS%Freq_Cen, start, counter)
      if (rcode_fr /= 0) then
        call MOM_error(FATAL,&
             "error reading dimension 1 values for var_name "// &
             trim(varread2)//",dim_name "//trim(dim_name(1))//  &
             " in file "// trim(SurfBandFileName)//" in MOM_wave_interface")
      endif
      CS%NUMBANDS = ID
      do B = 1,CS%NumBands
        CS%WaveNum_Cen(b) = (2.*PI*CS%Freq_Cen(b)*US%T_to_s)**2 / (US%L_to_Z**2*GV%g_Earth)
      enddo
    endif

    rcode_wn = NF90_close(ncid)
    if (rcode_wn /= 0) call MOM_error(WARNING, &
            "Error closing file "//trim(SurfBandFileName)//" in MOM_wave_interface.")

  endif

  do b = 1,CS%NumBands
    temp_x(:,:) = 0.0
    temp_y(:,:) = 0.0
    varname = '                    '
    write(varname,"(A3,I0)")'Usx',b
    call data_override('OCN',trim(varname), temp_x, day_center)
    varname = '                    '
    write(varname,'(A3,I0)')'Usy',b
    call data_override('OCN',trim(varname), temp_y, day_center)
    ! Disperse into halo on h-grid
    call pass_vector(temp_x, temp_y, G%Domain, To_All, AGRID)
    !Filter land values
    do j = G%jsd,G%jed
      do i = G%Isd,G%Ied
        if (abs(temp_x(i,j)) > 10. .or. abs(temp_y(i,j)) > 10.) then
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
    ! Disperse into halo on u/v grids
    call pass_vector(CS%STKx0(:,:,b),CS%STKy0(:,:,b), G%Domain, To_ALL)
  enddo !Closes b-loop

end subroutine Surface_Bands_by_data_override

!> Interface to get Langmuir number based on options stored in wave structure
!!
!! Note this can be called with an unallocated Waves pointer, which is okay if we
!!  want the wind-speed only dependent Langmuir number.  Therefore, we need to be
!!  careful about what we try to access here.
subroutine get_Langmuir_Number( LA, G, GV, US, HBL, ustar, i, j, &
                                H, U_H, V_H, Override_MA, Waves )
  type(ocean_grid_type),   intent(in) :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  integer, intent(in) :: i      !< Meridional index of h-point
  integer, intent(in) :: j      !< Zonal index of h-point
  real, intent(in)    :: ustar  !< Friction velocity [Z T-1 ~> m s-1].
  real, intent(in)    :: HBL    !< (Positive) thickness of boundary layer [Z ~> m].
  logical, optional,       intent(in) :: Override_MA !< Override to use misalignment in LA
                                !! calculation. This can be used if diagnostic
                                !! LA outputs are desired that are different than
                                !! those used by the dynamical model.
  real, dimension(SZK_(GV)), optional, &
       intent(in)      :: H     !< Grid layer thickness [H ~> m or kg m-2]
  real, dimension(SZK_(GV)), optional, &
       intent(in)      :: U_H   !< Zonal velocity at H point [m s-1]
  real, dimension(SZK_(GV)), optional, &
       intent(in)      :: V_H   !< Meridional velocity at H point [m s-1]
  type(Wave_parameters_CS), &
       pointer         :: Waves !< Surface wave control structure.

  real, intent(out)    :: LA    !< Langmuir number

!Local Variables
  real :: Top, bottom, midpoint
  real :: Dpt_LASL, ShearDirection, WaveDirection
  real :: LA_STKx, LA_STKy, LA_STK ! Stokes velocities in [m s-1]
  logical :: ContinueLoop, USE_MA
  real, dimension(SZK_(GV)) :: US_H, VS_H
  real, allocatable :: StkBand_X(:), StkBand_Y(:)
  integer :: KK, BB


 ! Compute averaging depth for Stokes drift (negative)
  Dpt_LASL = min(-0.1*US%m_to_Z, -LA_FracHBL*HBL)

  USE_MA = LA_Misalignment
  if (present(Override_MA)) USE_MA = Override_MA

  ! If requesting to use misalignment in the Langmuir number compute the Shear Direction
  if (USE_MA) then
    if (.not.(present(H).and.present(U_H).and.present(V_H))) then
      call MOM_error(Fatal,'Get_LA_waves requested to consider misalignment.')
    endif
    ContinueLoop = .true.
    bottom = 0.0
    do kk = 1,GV%ke
      Top = Bottom
      MidPoint = Bottom + GV%H_to_Z*0.5*h(kk)
      Bottom = Bottom + GV%H_to_Z*h(kk)
      if (MidPoint > Dpt_LASL .and. kk > 1 .and. ContinueLoop) then
        ShearDirection = atan2(V_H(1)-V_H(kk),U_H(1)-U_H(kk))
        ContinueLoop = .false.
      endif
    enddo
  endif

  if (WaveMethod==TESTPROF) then
    do kk = 1,GV%ke
      US_H(kk) = 0.5*(WAVES%US_X(I,j,kk)+WAVES%US_X(I-1,j,kk))
      VS_H(kk) = 0.5*(WAVES%US_Y(i,J,kk)+WAVES%US_Y(i,J-1,kk))
    enddo
    call Get_SL_Average_Prof( GV, Dpt_LASL, H, US_H, LA_STKx)
    call Get_SL_Average_Prof( GV, Dpt_LASL, H, VS_H, LA_STKy)
    LA_STK = sqrt(LA_STKX*LA_STKX+LA_STKY*LA_STKY)
  elseif (WaveMethod==SURFBANDS) then
    allocate(StkBand_X(WAVES%NumBands), StkBand_Y(WAVES%NumBands))
    do bb = 1,WAVES%NumBands
      StkBand_X(bb) = 0.5*(WAVES%STKx0(I,j,bb)+WAVES%STKx0(I-1,j,bb))
      StkBand_Y(bb) = 0.5*(WAVES%STKy0(i,J,bb)+WAVES%STKy0(i,J-1,bb))
    enddo
    call Get_SL_Average_Band(GV, Dpt_LASL, WAVES%NumBands, WAVES%WaveNum_Cen, StkBand_X, LA_STKx )
    call Get_SL_Average_Band(GV, Dpt_LASL, WAVES%NumBands, WAVES%WaveNum_Cen, StkBand_Y, LA_STKy )
    LA_STK = sqrt(LA_STKX**2 + LA_STKY**2)
    deallocate(StkBand_X, StkBand_Y)
  elseif (WaveMethod==DHH85) then
    ! Temporarily integrating profile rather than spectrum for simplicity
    do kk = 1,GV%ke
      US_H(kk) = 0.5*(WAVES%US_X(I,j,kk)+WAVES%US_X(I-1,j,kk))
      VS_H(kk) = 0.5*(WAVES%US_Y(i,J,kk)+WAVES%US_Y(i,J-1,kk))
    enddo
    call Get_SL_Average_Prof( GV, Dpt_LASL, H, US_H, LA_STKx)
    call Get_SL_Average_Prof( GV, Dpt_LASL, H, VS_H, LA_STKy)
    LA_STK = sqrt(LA_STKX**2 + LA_STKY**2)
  elseif (WaveMethod==LF17) then
    call get_StokesSL_LiFoxKemper(ustar, hbl*LA_FracHBL, GV, US, LA_STK, LA)
  elseif (WaveMethod==Null_WaveMethod) then
    call MOM_error(FATAL, "Get_Langmuir_number called without defining a WaveMethod. "//&
                          "Suggest to make sure USE_LT is set/overridden to False or "//&
                          "choose a wave method (or set USE_LA_LI2016 to use statistical "//&
                          "waves.")
  endif

  if (.not.(WaveMethod==LF17)) then
    ! This is an arbitrary lower bound on Langmuir number.
    ! We shouldn't expect values lower than this, but
    ! there is also no good reason to cap it here other then
    ! to prevent large enhancements in unconstrained parts of
    ! the curve fit parameterizations.
    ! Note the dimensional constant background Stokes velocity of 10^-10 m s-1.
    LA = max(WAVES%La_min, sqrt(US%Z_to_m*US%s_to_T*ustar / (LA_STK+1.e-10)))
  endif

  if (Use_MA) then
    WaveDirection = atan2(LA_STKy, LA_STKx)
    LA = LA / sqrt(max(1.e-8, cos( WaveDirection - ShearDirection)))
  endif

  return
end subroutine get_Langmuir_Number

!> Get SL averaged Stokes drift from Li/FK 17 method
!!
!! Original description:
!! - This function returns the enhancement factor, given the 10-meter
!!   wind [m s-1], friction velocity [m s-1] and the boundary layer depth [m].
!!
!! Update (Jan/25):
!! - Converted from function to subroutine, now returns Langmuir number.
!! - Computs 10m wind internally, so only ustar and hbl need passed to
!!   subroutine.
!!
!! Qing Li, 160606
!! - BGR port from CVMix to MOM6 Jan/25/2017
!! - BGR change output to LA from Efactor
!! - BGR remove u10 input
!! - BGR note: fixed parameter values should be changed to "get_params"
subroutine get_StokesSL_LiFoxKemper(ustar, hbl, GV, US, UStokes_SL, LA)
  real, intent(in)  :: ustar !< water-side surface friction velocity [Z T-1 ~> m s-1].
  real, intent(in)  :: hbl   !< boundary layer depth [Z ~> m].
  type(verticalGrid_type), intent(in) :: GV !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  real, intent(out) :: UStokes_SL !< Surface layer averaged Stokes drift [m s-1]
  real, intent(out) :: LA    !< Langmuir number
  ! Local variables
  ! parameters
  real, parameter :: &
       ! ratio of U19.5 to U10 (Holthuijsen, 2007)
       u19p5_to_u10 = 1.075, &
       ! ratio of mean frequency to peak frequency for
       ! Pierson-Moskowitz spectrum (Webb, 2011)
       fm_into_fp = 1.296, &
       ! ratio of surface Stokes drift to U10
       us_to_u10 = 0.0162, &
       ! loss ratio of Stokes transport
       r_loss = 0.667
  real :: UStokes, hm0, fm, fp, vstokes, kphil, kstar
  real :: z0, z0i, r1, r2, r3, r4, tmp, lasl_sqr_i
  real :: u10

  UStokes_sl = 0.0
  LA=1.e8
  if (ustar > 0.0) then
    ! Computing u10 based on u_star and COARE 3.5 relationships
    call ust_2_u10_coare3p5(US%Z_to_m*US%s_to_T*ustar*sqrt(US%R_to_kg_m3*GV%Rho0/1.225), u10, GV, US)
    ! surface Stokes drift
    UStokes = us_to_u10*u10
    !
    ! significant wave height from Pierson-Moskowitz
    ! spectrum (Bouws, 1998)
    hm0 = 0.0246 *u10**2
    !
    ! peak frequency (PM, Bouws, 1998)
    tmp = 2.0 * PI * u19p5_to_u10 * u10
    fp = 0.877 * US%L_T_to_m_s**2*US%m_to_Z * GV%g_Earth / tmp
    !
    ! mean frequency
    fm = fm_into_fp * fp
    !
    ! total Stokes transport (a factor r_loss is applied to account
    !  for the effect of directional spreading, multidirectional waves
    !  and the use of PM peak frequency and PM significant wave height
    !  on estimating the Stokes transport)
    vstokes = 0.125 * PI * r_loss * fm * hm0**2
    !
    ! the general peak wavenumber for Phillips' spectrum
    ! (Breivik et al., 2016) with correction of directional spreading
    kphil = 0.176 * UStokes / vstokes
    !
    ! surface layer averaged Stokes dirft with Stokes drift profile
    ! estimated from Phillips' spectrum (Breivik et al., 2016)
    ! the directional spreading effect from Webb and Fox-Kemper, 2015
    ! is also included
    kstar = kphil * 2.56
    ! surface layer
    z0 = abs(US%Z_to_m*hbl)
    z0i = 1.0 / z0
    ! term 1 to 4
    r1 = ( 0.151 / kphil * z0i -0.84 ) * &
         ( 1.0 - exp(-2.0 * kphil * z0) )
    r2 = -( 0.84 + 0.0591 / kphil * z0i ) * &
         sqrt( 2.0 * PI * kphil * z0 ) * &
         erfc( sqrt( 2.0 * kphil * z0 ) )
    r3 = ( 0.0632 / kstar * z0i + 0.125 ) * &
         (1.0 - exp(-2.0 * kstar * z0) )
    r4 = ( 0.125 + 0.0946 / kstar * z0i ) * &
         sqrt( 2.0 * PI *kstar * z0) * &
         erfc( sqrt( 2.0 * kstar * z0 ) )
    UStokes_sl = UStokes * (0.715 + r1 + r2 + r3 + r4)
    if(UStokes_sl .ne. 0.0)LA = sqrt(US%Z_to_m*US%s_to_T*ustar / UStokes_sl)
  endif

end subroutine Get_StokesSL_LiFoxKemper

!> Get SL Averaged Stokes drift from a Stokes drift Profile
subroutine Get_SL_Average_Prof( GV, AvgDepth, H, Profile, Average )
  type(verticalGrid_type),  &
       intent(in)   :: GV       !< Ocean vertical grid structure
  real, intent(in)  :: AvgDepth !< Depth to average over (negative) [Z ~> m].
  real, dimension(SZK_(GV)), &
       intent(in)   :: H        !< Grid thickness [H ~> m or kg m-2]
  real, dimension(SZK_(GV)), &
       intent(in)   :: Profile  !< Profile of quantity to be averaged [arbitrary]
                                !! (used here for Stokes drift)
  real, intent(out) :: Average  !< Output quantity averaged over depth AvgDepth [arbitrary]
                                !! (used here for Stokes drift)
  !Local variables
  real :: top, midpoint, bottom ! Depths, negative downward [Z ~> m].
  real :: Sum
  integer :: kk

  ! Initializing sum
  Sum = 0.0

  ! Integrate
  bottom = 0.0
  do kk = 1, GV%ke
    Top = Bottom
    MidPoint = Bottom - GV%H_to_Z * 0.5*h(kk)
    Bottom = Bottom - GV%H_to_Z * h(kk)
    if (AvgDepth < Bottom) then ! The whole cell is within H_LA
      Sum = Sum + Profile(kk) * (GV%H_to_Z * H(kk))
    elseif (AvgDepth < Top) then ! A partial cell is within H_LA
      Sum = Sum + Profile(kk) * (Top-AvgDepth)
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
       intent(in)     :: SurfStokes  !< Surface Stokes drift for each band [m s-1]
  real, intent(out)   :: Average     !< Output average Stokes drift over depth AvgDepth [m s-1]

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
  enddo

  return
end subroutine Get_SL_Average_Band

!> Compute the Stokes drift at a given depth
!!
!! Taken from Qing Li (Brown)
!! use for comparing MOM6 simulation to his LES
!! computed at z mid point (I think) and not depth averaged.
!! Should be fine to integrate in frequency from 0.1 to sqrt(-0.2*grav*2pi/dz
subroutine DHH85_mid(GV, US, zpt, UStokes)
  type(verticalGrid_type), intent(in)  :: GV  !< Ocean vertical grid
  type(unit_scale_type),   intent(in)  :: US  !< A dimensional unit scaling type
  real, intent(in)  :: zpt   !< Depth to get Stokes drift [Z ~> m].
  real, intent(out) :: UStokes !< Stokes drift [m s-1]
  !
  real :: ann, Bnn, Snn, Cnn, Dnn
  real :: omega_peak, omega, u10, WA, domega
  real :: omega_min, omega_max, wavespec, Stokes
  real :: g_Earth ! Gravitational acceleration [m s-2]
  integer :: Nomega, OI

  WA = WaveAge
  u10 = WaveWind
  g_Earth = US%L_T_to_m_s**2*US%m_to_Z * GV%g_Earth

  !/
  omega_min = 0.1 ! Hz
  ! Cut off at 30cm for now...
  omega_max = 10. ! ~sqrt(0.2*g_Earth*2*pi/0.3)
  NOmega = 1000
  domega = (omega_max-omega_min)/real(NOmega)

  !
  if (WaveAgePeakFreq) then
    omega_peak = g_Earth / (WA * u10)
  else
    omega_peak = 2. * pi * 0.13 * g_Earth / U10
  endif
  !/
  Ann = 0.006 * WaveAge**(-0.55)
  Bnn = 1.0
  Snn = 0.08 * (1.0 + 4.0 * WaveAge**3)
  Cnn = 1.7
  if (WA < 1.) then
    Cnn = Cnn - 6.0*log10(WA)
  endif
  !/
  UStokes = 0.0
  omega = omega_min + 0.5*domega
  do oi = 1,nomega-1
    Dnn = exp ( -0.5 * (omega-omega_peak)**2 / (Snn**2 * omega_peak**2) )
    ! wavespec units = m2s
    wavespec = (Ann * g_Earth**2 / (omega_peak*omega**4 ) ) * &
               exp(-bnn*(omega_peak/omega)**4)*Cnn**Dnn
    ! Stokes units m  (multiply by frequency range for units of m/s)
    Stokes = 2.0 * wavespec * omega**3 * &
         exp( 2.0 * omega**2 * US%Z_to_m*zpt / g_Earth) / g_Earth
    UStokes = UStokes + Stokes*domega
    omega = omega + domega
  enddo

  return
end subroutine DHH85_mid

!> Explicit solver for Stokes mixing.
!! Still in development do not use.
subroutine StokesMixing(G, GV, dt, h, u, v, Waves )
  type(ocean_grid_type), &
       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type), &
       intent(in)    :: GV    !< Ocean vertical grid
  real, intent(in)   :: dt    !< Time step of MOM6 [T ~> s] for explicit solver
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
       intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
       intent(inout) :: u     !< Velocity i-component [m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
       intent(inout) :: v     !< Velocity j-component [m s-1]
  type(Wave_parameters_CS), &
       pointer       :: Waves !< Surface wave related control structure.
  ! Local variables
  real :: dTauUp, dTauDn ! Vertical momentum fluxes [Z T-1 m s-1]
  real :: h_Lay  ! The layer thickness at a velocity point [Z ~> m].
  integer :: i,j,k

! This is a template to think about down-Stokes mixing.
! This is not ready for use...

  do k = 1, GV%ke
    do j = G%jsc, G%jec
      do I = G%iscB, G%iecB
        h_lay = GV%H_to_Z*0.5*(h(i,j,k)+h(i+1,j,k))
        dTauUp = 0.0
        if (k > 1) &
          dTauUp = 0.5*(waves%Kvs(i,j,k)+waves%Kvs(i+1,j,k)) * &
               (waves%us_x(i,j,k-1)-waves%us_x(i,j,k)) / &
               (0.5*(h_lay + GV%H_to_Z*0.5*(h(i,j,k-1)+h(i+1,j,k-1)) ))
        dTauDn = 0.0
        if (k < GV%ke-1) &
          dTauDn = 0.5*(waves%Kvs(i,j,k+1)+waves%Kvs(i+1,j,k+1)) * &
               (waves%us_x(i,j,k)-waves%us_x(i,j,k+1)) / &
               (0.5*(h_lay + GV%H_to_Z*0.5*(h(i,j,k+1)+h(i+1,j,k+1)) ))
        u(i,j,k) = u(i,j,k) + dt * (dTauUp-dTauDn) / h_Lay
      enddo
    enddo
  enddo

  do k = 1, GV%ke
    do J = G%jscB, G%jecB
      do i = G%isc, G%iec
        h_Lay = GV%H_to_Z*0.5*(h(i,j,k)+h(i,j+1,k))
        dTauUp = 0.
        if (k > 1) &
          dTauUp = 0.5*(waves%Kvs(i,j,k)+waves%Kvs(i,j+1,k)) * &
               (waves%us_y(i,j,k-1)-waves%us_y(i,j,k)) / &
               (0.5*(h_lay + GV%H_to_Z*0.5*(h(i,j,k-1)+h(i,j+1,k-1)) ))
        dTauDn = 0.0
        if (k < GV%ke-1) &
          dTauDn =0.5*(waves%Kvs(i,j,k+1)+waves%Kvs(i,j+1,k+1)) * &
               (waves%us_y(i,j,k)-waves%us_y(i,j,k+1)) / &
               (0.5*(h_lay + GV%H_to_Z*0.5*(h(i,j,k+1)+h(i,j+1,k+1)) ))
        v(i,J,k) = v(i,J,k) + dt * (dTauUp-dTauDn) / h_Lay
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
subroutine CoriolisStokes(G, GV, DT, h, u, v, WAVES, US)
  type(ocean_grid_type), &
       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type), &
       intent(in)   :: GV     !< Ocean vertical grid
  real, intent(in)  :: Dt     !< Time step of MOM6 [s] CHECK IF PASSING RIGHT TIMESTEP
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
       intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
       intent(inout) :: u     !< Velocity i-component [m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
       intent(inout) :: v     !< Velocity j-component [m s-1]
  type(Wave_parameters_CS), &
       pointer       :: Waves !< Surface wave related control structure.
  type(unit_scale_type),   intent(in) :: US     !< A dimensional unit scaling type
  ! Local variables
  real :: DVel ! A rescaled velocity change [m s-1 T-1 ~> m s-2]
  integer :: i,j,k

  do k = 1, GV%ke
    do j = G%jsc, G%jec
      do I = G%iscB, G%iecB
        DVel = 0.25*(WAVES%us_y(i,j+1,k)+WAVES%us_y(i-1,j+1,k))*G%CoriolisBu(i,j+1) + &
               0.25*(WAVES%us_y(i,j,k)+WAVES%us_y(i-1,j,k))*G%CoriolisBu(i,j)
        u(I,j,k) = u(I,j,k) + DVEL*US%s_to_T*DT
      enddo
    enddo
  enddo

  do k = 1, GV%ke
    do J = G%jscB, G%jecB
      do i = G%isc, G%iec
        DVel = 0.25*(WAVES%us_x(i+1,j,k)+WAVES%us_x(i+1,j-1,k))*G%CoriolisBu(i+1,j) + &
               0.25*(WAVES%us_x(i,j,k)+WAVES%us_x(i,j-1,k))*G%CoriolisBu(i,j)
        v(i,J,k) = v(i,j,k) - DVEL*US%s_to_T*DT
      enddo
    enddo
  enddo
end subroutine CoriolisStokes

!> Computes wind speed from ustar_air based on COARE 3.5 Cd relationship
!! Probably doesn't belong in this module, but it is used here to estimate
!! wind speed for wind-wave relationships.  Should be a fine way to estimate
!! the neutral wind-speed as written here.
subroutine ust_2_u10_coare3p5(USTair, U10, GV, US)
  real, intent(in)                    :: USTair !< Wind friction velocity [m s-1]
  real, intent(out)                   :: U10    !< 10-m neutral wind speed [m s-1]
  type(verticalGrid_type), intent(in) :: GV     !< vertical grid type
  type(unit_scale_type),   intent(in) :: US     !< A dimensional unit scaling type

  ! Local variables
  real, parameter :: vonkar = 0.4 ! Should access a get_param von karman
  real, parameter :: nu=1e-6 ! Should access a get_param air-viscosity
  real :: z0sm, z0, z0rough, u10a, alpha, CD
  integer :: CT

  ! Uses empirical formula for z0 to convert ustar_air to u10 based on the
  !  COARE 3.5 paper (Edson et al., 2013)
  ! alpha=m*U10+b
  ! Note in Edson et al. 2013, eq. 13 m is given as 0.017.  However,
  ! m=0.0017 reproduces the curve in their figure 6.

  z0sm = 0.11 * nu * US%m_to_Z / USTair !Compute z0smooth from ustar guess
  u10 = USTair/sqrt(0.001)  !Guess for u10
  u10a = 1000

  CT=0
  do while (abs(u10a/u10-1.) > 0.001)
    CT=CT+1
    u10a = u10
    alpha = min(0.028, 0.0017 * u10 - 0.005)
    z0rough = alpha * (US%m_s_to_L_T*USTair)**2 / GV%g_Earth ! Compute z0rough from ustar guess
    z0 = z0sm + z0rough
    CD = ( vonkar / log(10.*US%m_to_Z / z0) )**2 ! Compute CD from derived roughness
    u10 = USTair/sqrt(CD)  ! Compute new u10 from derived CD, while loop
                           ! ends and checks for convergence...CT counter
                           ! makes sure loop doesn't run away if function
                           ! doesn't converge.  This code was produced offline
                           ! and converged rapidly (e.g. 2 cycles)
                           ! for ustar=0.0001:0.0001:10.
    if (CT>20) then
      u10 = USTair/sqrt(0.0015) ! I don't expect to get here, but just
                              !  in case it will output a reasonable value.
      exit
    endif
  enddo
  return
end subroutine ust_2_u10_coare3p5

!> Clear pointers, deallocate memory
subroutine Waves_end(CS)
  type(wave_parameters_CS), pointer :: CS !< Control structure

  if (allocated(CS%WaveNum_Cen)) deallocate( CS%WaveNum_Cen )
  if (allocated(CS%Freq_Cen))    deallocate( CS%Freq_Cen )
  if (allocated(CS%Us_x))        deallocate( CS%Us_x )
  if (allocated(CS%Us_y))        deallocate( CS%Us_y )
  if (allocated(CS%La_SL))       deallocate( CS%La_SL )
  if (allocated(CS%La_turb))     deallocate( CS%La_turb )
  if (allocated(CS%STKx0))       deallocate( CS%STKx0 )
  if (allocated(CS%STKy0))       deallocate( CS%STKy0 )
  if (allocated(CS%KvS))         deallocate( CS%KvS )
  if (allocated(CS%Us0_y))       deallocate( CS%Us0_y )
  if (allocated(CS%Us0_x))       deallocate( CS%Us0_x )

  deallocate( CS )

  return
end subroutine Waves_end

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
!! files).  In short order, the plan is to also ammend the subroutine
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
