!> Provides the K-Profile Parameterization (KPP) of Large et al., 1994, via CVMix.
module MOM_CVMix_KPP

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms,           only : max_across_PEs
use MOM_debugging,      only : hchksum, is_NaN
use MOM_diag_mediator,  only : time_type, diag_ctrl, safe_alloc_ptr, post_data
use MOM_diag_mediator,  only : query_averaging_enabled, register_diag_field
use MOM_error_handler,  only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_PE
use MOM_EOS,            only : EOS_type, calculate_density
use MOM_file_parser,    only : get_param, log_param, log_version, param_file_type
use MOM_file_parser,    only : openParameterBlock, closeParameterBlock
use MOM_grid,           only : ocean_grid_type, isPointInCell
use MOM_interface_heights, only : thickness_to_dz
use MOM_restart,        only : MOM_restart_CS, register_restart_field
use MOM_unit_scaling,   only : unit_scale_type
use MOM_variables,      only : thermo_var_ptrs
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS, Get_Langmuir_Number, get_wave_method
use MOM_domains,        only : pass_var
use MOM_cpu_clock,      only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,      only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_tracer_types,   only : tracer_type

use CVMix_kpp, only : CVMix_init_kpp, CVMix_put_kpp, CVMix_get_kpp_real
use CVMix_kpp, only : CVMix_coeffs_kpp
use CVMix_kpp, only : CVMix_kpp_compute_OBL_depth
use CVMix_kpp, only : CVMix_kpp_compute_turbulent_scales
use CVMix_kpp, only : CVMix_kpp_compute_bulk_Richardson
use CVMix_kpp, only : CVMix_kpp_compute_unresolved_shear
use CVMix_kpp, only : CVMix_kpp_params_type
use CVMix_kpp, only : CVMix_kpp_compute_kOBL_depth
use CVMix_kpp, only : CVMix_kpp_compute_StokesXi

implicit none ; private

#include "MOM_memory.h"

public :: register_KPP_restarts
public :: KPP_init
public :: KPP_compute_BLD
public :: KPP_calculate
public :: KPP_end
public :: KPP_NonLocalTransport_temp
public :: KPP_NonLocalTransport_saln
public :: KPP_NonLocalTransport
public :: KPP_get_BLD

! Enumerated constants
integer, private, parameter :: NLT_SHAPE_CVMix     = 0 !< Use the CVMix profile
integer, private, parameter :: NLT_SHAPE_LINEAR    = 1 !< Linear, \f$ G(\sigma) = 1-\sigma \f$
integer, private, parameter :: NLT_SHAPE_PARABOLIC = 2 !< Parabolic, \f$ G(\sigma) = (1-\sigma)^2 \f$
integer, private, parameter :: NLT_SHAPE_CUBIC     = 3 !< Cubic, \f$ G(\sigma) = 1 + (2\sigma-3) \sigma^2\f$
integer, private, parameter :: NLT_SHAPE_CUBIC_LMD = 4 !< Original shape,
                                                       !!    \f$ G(\sigma) = \frac{27}{4} \sigma (1-\sigma)^2 \f$

integer, private, parameter :: SW_METHOD_ALL_SW = 0 !< Use all shortwave radiation
integer, private, parameter :: SW_METHOD_MXL_SW = 1 !< Use shortwave radiation absorbed in mixing layer
integer, private, parameter :: SW_METHOD_LV1_SW = 2 !< Use shortwave radiation absorbed in layer 1
integer, private, parameter :: LT_K_CONSTANT = 1,        & !< Constant enhance K through column
                               LT_K_SCALED = 2,          & !< Enhance K scales with G(sigma)
                               LT_K_MODE_CONSTANT = 1,   & !< Prescribed enhancement for K
                               LT_K_MODE_VR12 = 2,       & !< Enhancement for K based on
                                                           !! Van Roekel et al., 2012
                               LT_K_MODE_RW16 = 3,       & !< Enhancement for K based on
                                                           !! Reichl et al., 2016
                               LT_VT2_MODE_CONSTANT = 1, & !< Prescribed enhancement for Vt2
                               LT_VT2_MODE_VR12 = 2,     & !< Enhancement for Vt2 based on
                                                           !! Van Roekel et al., 2012
                               LT_VT2_MODE_RW16 = 3,     & !< Enhancement for Vt2 based on
                                                           !! Reichl et al., 2016
                               LT_VT2_MODE_LF17 = 4        !< Enhancement for Vt2 based on
                                                           !! Li and Fox-Kemper, 2017

!> Control structure for containing KPP parameters/data
type, public :: KPP_CS ; private

  ! Parameters
  real    :: Ri_crit                   !< Critical bulk Richardson number (defines OBL depth) [nondim]
  real    :: vonKarman                 !< von Karman constant (dimensionless) [nondim]
  real    :: cs                        !< Parameter for computing velocity scale function (dimensionless) [nondim]
  real    :: cs2                       !< Parameter for multiplying by non-local term [nondim]
                                       !   This is active for NLT_SHAPE_CUBIC_LMD only
  logical :: enhance_diffusion         !< If True, add enhanced diffusivity at base of boundary layer.
  character(len=32) :: interpType      !< Type of interpolation to compute bulk Richardson number
  character(len=32) :: interpType2     !< Type of interpolation to compute diff and visc at OBL_depth
  logical :: StokesMOST                !< If True, use Stokes similarity package
  logical :: computeEkman              !< If True, compute Ekman depth limit for OBLdepth
  logical :: computeMoninObukhov       !< If True, compute Monin-Obukhov limit for OBLdepth
  logical :: passiveMode               !< If True, makes KPP passive meaning it does NOT alter the diffusivity
  real    :: deepOBLoffset             !< If non-zero, is a distance from the bottom that the OBL can not
                                       !! penetrate through [Z ~> m]
  real    :: minOBLdepth               !< If non-zero, is a minimum depth for the OBL [Z ~> m]
  real    :: surf_layer_ext            !< Fraction of OBL depth considered in the surface layer [nondim]
  real    :: minVtsqr                  !< Min for the squared unresolved velocity used in Rib CVMix
                                       !! calculation [L2 T-2 ~> m2 s-2]
  logical :: fixedOBLdepth             !< If True, will fix the OBL depth at fixedOBLdepth_value
  real    :: fixedOBLdepth_value       !< value for the fixed OBL depth when fixedOBLdepth==True [Z ~> m]
  logical :: debug                     !< If True, calculate checksums and write debugging information
  character(len=30) :: MatchTechnique  !< Method used in CVMix for setting diffusivity and NLT profile functions
  integer :: NLT_shape                 !< MOM6 over-ride of CVMix NLT shape function
  logical :: applyNonLocalTrans        !< If True, apply non-local transport to all tracers
  integer :: n_smooth                  !< Number of times smoothing operator is applied on OBLdepth.
  logical :: deepen_only               !< If true, apply OBLdepth smoothing at a cell only if the OBLdepth gets deeper.
  logical :: KPPzeroDiffusivity        !< If True, will set diffusivity and viscosity from KPP to zero
                                       !! for testing purposes.
  logical :: KPPisAdditive             !< If True, will add KPP diffusivity to initial diffusivity.
                                       !! If False, will replace initial diffusivity wherever KPP diffusivity
                                       !! is non-zero.
  real    :: min_thickness             !< A minimum thickness used to avoid division by small numbers
                                       !! in the vicinity of vanished layers [Z ~> m]
  integer :: SW_METHOD                 !< Sets method for using shortwave radiation in surface buoyancy flux
  logical :: LT_K_Enhancement          !< Flags if enhancing mixing coefficients due to LT
  integer :: LT_K_Shape                !< Integer for constant or shape function enhancement
  integer :: LT_K_Method               !< Integer for mixing coefficients LT method
  real    :: KPP_CVt2                  !< Parameter for Stokes MOST convection entrainment [nondim]
  real    :: KPP_K_ENH_FAC             !< Factor to multiply by K if Method is CONSTANT [nondim]
  logical :: LT_Vt2_Enhancement        !< Flags if enhancing Vt2 due to LT
  integer :: LT_VT2_METHOD             !< Integer for Vt2 LT method
  real    :: KPP_VT2_ENH_FAC           !< Factor to multiply by VT2 if Method is CONSTANT [nondim]
  real    :: MLD_guess_min             !< The minimum estimate of the mixed layer depth used to
                                       !! calculate the Langmuir number for Langmuir turbulence
                                       !! enhancement with KPP [Z ~> m]
  logical :: STOKES_MIXING             !< Flag if model is mixing down Stokes gradient
                                       !! This is relevant for which current to use in RiB
  integer :: answer_date               !< The vintage of the order of arithmetic in the CVMix KPP
                                       !! calculations.  Values below 20240501 recover the answers
                                       !! from early in 2024, while higher values use expressions
                                       !! that have been refactored for rotational symmetry.

  !> CVMix parameters
  type(CVMix_kpp_params_type), pointer :: KPP_params => NULL()

  type(diag_ctrl), pointer :: diag => NULL() !< Pointer to diagnostics control structure
  !>@{ Diagnostic handles
  integer :: id_OBLdepth = -1, id_BulkRi   = -1
  integer :: id_N        = -1, id_N2       = -1
  integer :: id_Ws       = -1, id_Vt2      = -1
  integer :: id_BulkUz2  = -1, id_BulkDrho = -1
  integer :: id_uStar    = -1, id_buoyFlux = -1
  integer :: id_sigma    = -1, id_Kv_KPP   = -1
  integer :: id_Kt_KPP   = -1, id_Ks_KPP   = -1
  integer :: id_Tsurf    = -1, id_Ssurf    = -1
  integer :: id_Usurf    = -1, id_Vsurf    = -1
  integer :: id_Kd_in    = -1
  integer :: id_NLTt     = -1
  integer :: id_NLTs     = -1
  integer :: id_EnhK     = -1, id_EnhVt2   = -1
  integer :: id_EnhW     = -1
  integer :: id_La_SL    = -1
  integer :: id_OBLdepth_original = -1
  integer :: id_StokesXI = -1
  integer :: id_Lam2     = -1
  !>@}

  ! Diagnostics arrays
  real, pointer,     dimension(:,:)   :: OBLdepth  !< Depth (positive) of ocean boundary layer (OBL) [Z ~> m]
  real, allocatable, dimension(:,:)   :: OBLdepth_original  !< Depth (positive) of OBL without smoothing [Z ~> m]
  real, allocatable, dimension(:,:)   :: StokesParXI !< Stokes similarity parameter [nondim]
  real, allocatable, dimension(:,:)   :: Lam2      !< La^(-2) = Ustk0/u* [nondim]
  real, allocatable, dimension(:,:)   :: kOBL      !< Level (+fraction) of OBL extent [nondim]
  real, allocatable, dimension(:,:)   :: OBLdepthprev !< previous Depth (positive) of OBL [Z ~> m]
  real, allocatable, dimension(:,:)   :: La_SL     !< Langmuir number used in KPP [nondim]
  real, allocatable, dimension(:,:,:) :: dRho      !< Bulk difference in density [R ~> kg m-3]
  real, allocatable, dimension(:,:,:) :: Uz2       !< Square of bulk difference in resolved velocity [L2 T-2 ~> m2 s-2]
  real, allocatable, dimension(:,:,:) :: BulkRi    !< Bulk Richardson number for each layer [nondim]
  real, allocatable, dimension(:,:,:) :: sigma     !< Sigma coordinate (dimensionless) [nondim]
  real, allocatable, dimension(:,:,:) :: Ws        !< Turbulent velocity scale for scalars [Z T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:) :: N         !< Brunt-Vaisala frequency [T-1 ~> s-1]
  real, allocatable, dimension(:,:,:) :: N2        !< Squared Brunt-Vaisala frequency [T-2 ~> s-2]
  real, allocatable, dimension(:,:,:) :: Vt2       !< Unresolved squared turbulence velocity for
                                                   !! bulk Ri [Z2 T-2 ~> m2 s-2]
  real, allocatable, dimension(:,:,:) :: Kt_KPP    !< Temp diffusivity from KPP [Z2 T-1 ~> m2 s-1]
  real, allocatable, dimension(:,:,:) :: Ks_KPP    !< Scalar diffusivity from KPP [Z2 T-1 ~> m2 s-1]
  real, allocatable, dimension(:,:,:) :: Kv_KPP    !< Viscosity due to KPP [Z2 T-1 ~> m2 s-1]
  real, allocatable, dimension(:,:)   :: Tsurf     !< Temperature of surface layer [C ~> degC]
  real, allocatable, dimension(:,:)   :: Ssurf     !< Salinity of surface layer [S ~> ppt]
  real, allocatable, dimension(:,:)   :: Usurf     !< i-velocity of surface layer [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:)   :: Vsurf     !< j-velocity of surface layer [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:) :: EnhK      !< Enhancement for mixing coefficient [nondim]
  real, allocatable, dimension(:,:,:) :: EnhVt2    !< Enhancement for Vt2 [nondim]

end type KPP_CS

!>@{ CPU time clocks
integer :: id_clock_KPP_calc, id_clock_KPP_compute_BLD, id_clock_KPP_smoothing
!>@}

#define __DO_SAFETY_CHECKS__

contains

!> Routine to register restarts, pass-through to children modules
subroutine register_KPP_restarts(G, param_file, restart_CSp, CS)
  type(ocean_grid_type), intent(in)    :: G           !< The ocean's grid structure
  type(param_file_type), intent(in)    :: param_file  !< A structure to parse for run-time parameters
  type(MOM_restart_CS),  pointer       :: restart_CSp  !< MOM restart control structure
  type(KPP_CS),         pointer        :: CS           !< module control structure

  character(len=40) :: mdl = 'MOM_CVMix_KPP' !< name of this module
  logical :: use_kpp, fpmix

  if (associated(CS)) call MOM_error(FATAL, 'MOM_CVMix_KPP, register_KPP_restarts: '// &
           'Control structure has already been initialized')
  call get_param(param_file, mdl, "USE_KPP", use_kpp, default=.false., do_not_log=.true.)
  ! Forego remainder of initialization if not using this scheme
  if (.not. use_kpp) return
  allocate(CS)

  allocate(CS%OBLdepth(SZI_(G),SZJ_(G)), source=0.0)

  ! FPMIX is needed to decide if boundary layer depth should be added to restart file
  call get_param(param_file, '', "FPMIX", fpmix, &
                 "If true, add non-local momentum flux increments and diffuse down the Eulerian gradient.", &
                 default=.false., do_not_log=.true.)
    if (fpmix) call register_restart_field(CS%OBLdepth, 'KPP_OBLdepth', .false., restart_CSp)

end subroutine register_KPP_restarts

!> Initialize the CVMix KPP module and set up diagnostics
!! Returns True if KPP is to be used, False otherwise.
logical function KPP_init(paramFile, G, GV, US, diag, Time, CS, passive)

  ! Arguments
  type(param_file_type),   intent(in)    :: paramFile !< File parser
  type(ocean_grid_type),   intent(in)    :: G         !< Ocean grid
  type(verticalGrid_type), intent(in)    :: GV        !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US        !< A dimensional unit scaling type
  type(diag_ctrl), target, intent(in)    :: diag      !< Diagnostics
  type(time_type),         intent(in)    :: Time      !< Model time
  type(KPP_CS),            pointer       :: CS        !< Control structure
  logical,       optional, intent(out)   :: passive   !< Copy of %passiveMode

  ! Local variables
# include "version_variable.h"
  character(len=40) :: mdl = 'MOM_CVMix_KPP' !< name of this module
  character(len=20) :: string          !< local temporary string
  character(len=20) :: langmuir_mixing_opt = 'NONE' !< Langmuir mixing option to be passed to CVMix, e.g., LWF16
  character(len=20) :: langmuir_entrainment_opt = 'NONE' !< Langmuir entrainment option to be
                                       !! passed to CVMix, e.g., LWF16
  integer :: default_answer_date       ! The default setting for the various ANSWER_DATE flags.
  logical :: CS_IS_ONE=.false.         !< Logical for setting Cs based on Non-local
  logical :: lnoDGat1=.false.          !< True => G'(1) = 0 (shape function)
                                       !! False => compute G'(1) as in LMD94
  ! Read parameters
  call get_param(paramFile, mdl, "USE_KPP", KPP_init, default=.false., do_not_log=.true.)
  call log_version(paramFile, mdl, version, 'This is the MOM wrapper to CVMix:KPP\n' // &
            'See http://cvmix.github.io/', all_default=.not.KPP_init)
  call get_param(paramFile, mdl, "USE_KPP", KPP_init, &
                 "If true, turns on the [CVMix] KPP scheme of Large et al., 1994, "// &
                 "to calculate diffusivities and non-local transport in the OBL.",     &
                 default=.false.)
  ! Forego remainder of initialization if not using this scheme
  if (.not. KPP_init) return

  call get_param(paramFile, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231, do_not_log=.true.)

  call openParameterBlock(paramFile,'KPP')
  call get_param(paramFile, mdl, 'PASSIVE', CS%passiveMode,           &
                 'If True, puts KPP into a passive-diagnostic mode.', &
                  default=.False.)
  !BGR: Note using PASSIVE for KPP creates warning for PASSIVE from Convection
  !     should we create a separate flag?
  if (present(passive)) passive=CS%passiveMode ! This is passed back to the caller so
                                               ! the caller knows to not use KPP output
  call get_param(paramFile, mdl, 'APPLY_NONLOCAL_TRANSPORT', CS%applyNonLocalTrans,  &
                 'If True, applies the non-local transport to all tracers. '//  &
                 'If False, calculates the non-local transport and tendencies but '//&
                 'purely for diagnostic purposes.',                                   &
                 default=.not. CS%passiveMode)
  call get_param(paramFile, mdl, 'N_SMOOTH', CS%n_smooth,  &
                 'The number of times the 1-1-4-1-1 Laplacian filter is applied on OBL depth.', &
                 default=0)
  if (CS%n_smooth > G%domain%nihalo) then
    call MOM_error(FATAL,'KPP smoothing number (N_SMOOTH) cannot be greater than NIHALO.')
  elseif (CS%n_smooth > G%domain%njhalo) then
    call MOM_error(FATAL,'KPP smoothing number (N_SMOOTH) cannot be greater than NJHALO.')
  endif
  if (CS%n_smooth > 0) then
    call get_param(paramFile, mdl, 'DEEPEN_ONLY_VIA_SMOOTHING', CS%deepen_only,  &
                   'If true, apply OBLdepth smoothing at a cell only if the OBLdepth '// &
                   'gets deeper via smoothing.',   &
                   default=.false.)
    id_clock_KPP_smoothing = cpu_clock_id('(Ocean KPP BLD smoothing)', grain=CLOCK_ROUTINE)
  endif
  call get_param(paramFile, mdl, 'RI_CRIT', CS%Ri_crit,                            &
                 'Critical bulk Richardson number used to define depth of the '// &
                 'surface Ocean Boundary Layer (OBL).',                            &
                 units='nondim', default=0.3)
  call get_param(paramFile, mdl, 'VON_KARMAN', CS%vonKarman, &
                 'von Karman constant.',                     &
                 units='nondim', default=0.40)
  call get_param(paramFile, mdl, 'ENHANCE_DIFFUSION', CS%enhance_diffusion,              &
                 'If True, adds enhanced diffusion at the based of the boundary layer.', &
                 default=.true.)
  call get_param(paramFile, mdl, 'INTERP_TYPE', CS%interpType,           &
                 'Type of interpolation to determine the OBL depth.\n'// &
                 'Allowed types are: linear, quadratic, cubic.',         &
                 default='quadratic')
  call get_param(paramFile, mdl, 'INTERP_TYPE2', CS%interpType2,           &
                 'Type of interpolation to compute diff and visc at OBL_depth.\n'// &
                 'Allowed types are: linear, quadratic, cubic or LMD94.',         &
                 default='LMD94')
  call get_param(paramFile, mdl, 'STOKES_MOST', CS%StokesMOST,             &
                 'If True, use Stokes Similarity package.', &
                 default=.False.)
  call get_param(paramFile, mdl, 'COMPUTE_EKMAN', CS%computeEkman,             &
                 'If True, limit OBL depth to be no deeper than Ekman depth.', &
                 default=.False.)
  call get_param(paramFile, mdl, 'COMPUTE_MONIN_OBUKHOV', CS%computeMoninObukhov, &
                 'If True, limit the OBL depth to be no deeper than '//          &
                 'Monin-Obukhov depth.',                                          &
                 default=.False.)
  call get_param(paramFile, mdl, 'CS', CS%cs,                        &
                 'Parameter for computing velocity scale function.', &
                 units='nondim', default=98.96)
  call get_param(paramFile, mdl, 'CS2', CS%cs2,                        &
                 'Parameter for computing non-local term.', &
                 units='nondim', default=6.32739901508)
  call get_param(paramFile, mdl, 'DEEP_OBL_OFFSET', CS%deepOBLoffset,                             &
                 'If non-zero, the distance above the bottom to which the OBL is clipped '//     &
                 'if it would otherwise reach the bottom. The smaller of this and 0.1D is used.', &
                 units='m', default=0., scale=US%m_to_Z)
  call get_param(paramFile, mdl, 'FIXED_OBLDEPTH', CS%fixedOBLdepth,       &
                 'If True, fix the OBL depth to FIXED_OBLDEPTH_VALUE '//  &
                 'rather than using the OBL depth from CVMix. '//         &
                 'This option is just for testing purposes.',              &
                 default=.False.)
  call get_param(paramFile, mdl, 'FIXED_OBLDEPTH_VALUE', CS%fixedOBLdepth_value,  &
                 'Value for the fixed OBL depth when fixedOBLdepth==True. '//   &
                 'This parameter is for just for testing purposes. '//          &
                 'It will over-ride the OBLdepth computed from CVMix.',           &
                 units='m', default=30.0, scale=US%m_to_Z)
  call get_param(paramFile, mdl, 'SURF_LAYER_EXTENT', CS%surf_layer_ext,   &
                 'Fraction of OBL depth considered in the surface layer.', &
                 units='nondim', default=0.10)
  call get_param(paramFile, mdl, 'MINIMUM_OBL_DEPTH', CS%minOBLdepth,                            &
                 'If non-zero, a minimum depth to use for KPP OBL depth. Independent of '//     &
                 'this parameter, the OBL depth is always at least as deep as the first layer.', &
                 units='m', default=0., scale=US%m_to_Z)
  call get_param(paramFile, mdl, 'MINIMUM_VT2', CS%minVtsqr,                                   &
                 'Min of the unresolved velocity Vt2 used in Rib CVMix calculation.\n'//  &
                 'Scaling: MINIMUM_VT2 = const1*d*N*ws, with d=1m, N=1e-5/s, ws=1e-6 m/s.',    &
                 units='m2/s2', default=1e-10, scale=US%m_s_to_L_T**2)

  call get_param(paramFile, mdl, 'NLT_SHAPE', string, &
                 'MOM6 method to set nonlocal transport profile. '//                          &
                 'Over-rides the result from CVMix.  Allowed values are: \n'//                 &
                 '\t CVMix     - Uses the profiles from CVMix specified by MATCH_TECHNIQUE\n'//&
                 '\t LINEAR    - A linear profile, 1-sigma\n'//                                &
                 '\t PARABOLIC - A parablic profile, (1-sigma)^2\n'//                          &
                 '\t CUBIC     - A cubic profile, (1-sigma)^2(1+2*sigma)\n'//                  &
                 '\t CUBIC_LMD - The original KPP profile',                                    &
                 default='CVMix')
  select case ( trim(string) )
    case ("CVMix")     ; CS%NLT_shape = NLT_SHAPE_CVMix
    case ("LINEAR")    ; CS%NLT_shape = NLT_SHAPE_LINEAR
    case ("PARABOLIC") ; CS%NLT_shape = NLT_SHAPE_PARABOLIC
    case ("CUBIC")     ; CS%NLT_shape = NLT_SHAPE_CUBIC
    case ("CUBIC_LMD") ; CS%NLT_shape = NLT_SHAPE_CUBIC_LMD
    case default ; call MOM_error(FATAL,"KPP_init: "// &
                   "Unrecognized NLT_SHAPE option"//trim(string))
  end select
  call get_param(paramFile, mdl, 'MATCH_TECHNIQUE', CS%MatchTechnique,                                    &
                 'CVMix method to set profile function for diffusivity and NLT, '//                      &
                 'as well as matching across OBL base. Allowed values are: \n'//                          &
                 '\t SimpleShapes      = sigma*(1-sigma)^2 for both diffusivity and NLT\n'//              &
                 '\t MatchGradient     = sigma*(1-sigma)^2 for NLT; diffusivity profile from matching\n'//&
                 '\t MatchBoth         = match gradient for both diffusivity and NLT\n'//                 &
                 '\t ParabolicNonLocal = sigma*(1-sigma)^2 for diffusivity; (1-sigma)^2 for NLT',         &
                 default='SimpleShapes')
  if (CS%MatchTechnique == 'ParabolicNonLocal') then
    ! This forces Cs2 (Cs in non-local computation) to equal 1 for parabolic non-local option.
    !  May be used during CVMix initialization.
    Cs_is_one=.true.
  endif
  if (CS%MatchTechnique == 'ParabolicNonLocal' .or. CS%MatchTechnique == 'SimpleShapes') then
    ! if gradient won't be matched, lnoDGat1=.true.
    lnoDGat1=.true.
  endif

  ! safety check to avoid negative diff/visc
  if (CS%MatchTechnique == 'MatchBoth' .and. (CS%interpType2 == 'cubic' .or. &
      CS%interpType2 == 'quadratic')) then
    call MOM_error(FATAL,"If MATCH_TECHNIQUE=MatchBoth, INTERP_TYPE2 must be set to \n"//&
               "linear or LMD94 (recommended) to avoid negative viscosity and diffusivity.\n"//&
               "Please select one of these valid options." )
  endif

  call get_param(paramFile, mdl, 'KPP_ZERO_DIFFUSIVITY', CS%KPPzeroDiffusivity,            &
                 'If True, zeroes the KPP diffusivity and viscosity; for testing purpose.',&
                 default=.False.)
  call get_param(paramFile, mdl, 'KPP_IS_ADDITIVE', CS%KPPisAdditive,                &
                 'If true, adds KPP diffusivity to diffusivity from other schemes.\n'//&
                 'If false, KPP is the only diffusivity wherever KPP is non-zero.',  &
                 default=.True.)
  call get_param(paramFile, mdl, 'KPP_SHORTWAVE_METHOD',string,                      &
                 'Determines contribution of shortwave radiation to KPP surface '// &
                 'buoyancy flux.  Options include:\n'//                             &
                 '  ALL_SW: use total shortwave radiation\n'//                      &
                 '  MXL_SW: use shortwave radiation absorbed by mixing layer\n'//  &
                 '  LV1_SW: use shortwave radiation absorbed by top model layer',  &
                 default='MXL_SW')
  select case ( trim(string) )
    case ("ALL_SW") ; CS%SW_METHOD = SW_METHOD_ALL_SW
    case ("MXL_SW") ; CS%SW_METHOD = SW_METHOD_MXL_SW
    case ("LV1_SW") ; CS%SW_METHOD = SW_METHOD_LV1_SW
    case default ; call MOM_error(FATAL,"KPP_init: "// &
                   "Unrecognized KPP_SHORTWAVE_METHOD option"//trim(string))
  end select
  call get_param(paramFile, mdl, 'CVMix_ZERO_H_WORK_AROUND', CS%min_thickness,                           &
                 'A minimum thickness used to avoid division by small numbers in the vicinity '//       &
                 'of vanished layers. This is independent of MIN_THICKNESS used in other parts of MOM.', &
                 units='m', default=0., scale=US%m_to_Z)

!/BGR: New options for including Langmuir effects
!/ 1. Options related to enhancing the mixing coefficient
  call get_param(paramFile, mdl, "USE_KPP_LT_K", CS%LT_K_Enhancement, &
       'Flag for Langmuir turbulence enhancement of turbulent'//&
       'mixing coefficient.', Default=.false.)
  call get_param(paramFile, mdl, "STOKES_MIXING", CS%Stokes_Mixing, &
       'Flag for Langmuir turbulence enhancement of turbulent'//&
       'mixing coefficient.', Default=.false.)
  if (CS%LT_K_Enhancement) then
    call get_param(paramFile, mdl, 'KPP_LT_K_SHAPE', string,                 &
                 'Vertical dependence of LT enhancement of mixing. '//     &
                 'Valid options are: \n'//                                   &
                 '\t CONSTANT = Constant value for full OBL\n'//             &
                 '\t SCALED   = Varies based on normalized shape function.', &
                 default='CONSTANT')
    select case ( trim(string))
      case ("CONSTANT") ; CS%LT_K_SHAPE = LT_K_CONSTANT
      case ("SCALED")   ; CS%LT_K_SHAPE = LT_K_SCALED
      case default ; call MOM_error(FATAL,"KPP_init: "//&
                    "Unrecognized KPP_LT_K_SHAPE option: "//trim(string))
    end select
    call get_param(paramFile, mdl, "KPP_LT_K_METHOD", string ,                   &
                   'Method to enhance mixing coefficient in KPP. '//             &
                   'Valid options are: \n'//                                     &
                   '\t CONSTANT = Constant value (KPP_K_ENH_FAC) \n'//           &
                   '\t VR12     = Function of Langmuir number based on VR12\n'// &
                   '\t            (Van Roekel et al. 2012)\n'//                  &
                   '\t            (Li et al. 2016, OM) \n'//                     &
                   '\t RW16     = Function of Langmuir number based on RW16\n'// &
                   '\t            (Reichl et al., 2016, JPO)',    &
                   default='CONSTANT')
    select case ( trim(string))
      case ("CONSTANT")
        CS%LT_K_METHOD = LT_K_MODE_CONSTANT
        langmuir_mixing_opt = 'LWF16'
      case ("VR12")
        CS%LT_K_METHOD = LT_K_MODE_VR12
        langmuir_mixing_opt = 'LWF16'
      case ("RW16")
        CS%LT_K_METHOD = LT_K_MODE_RW16
        langmuir_mixing_opt = 'RWHGK16'
      case default
        call MOM_error(FATAL,"KPP_init: "//&
                    "Unrecognized KPP_LT_K_METHOD option: "//trim(string))
    end select
    if (CS%LT_K_METHOD==LT_K_MODE_CONSTANT) then
      call get_param(paramFile, mdl, "KPP_K_ENH_FAC", CS%KPP_K_ENH_FAC,    &
                   'Constant value to enhance mixing coefficient in KPP.', &
                   units="nondim", default=1.0)
    endif
  endif
!/ 2. Options related to enhancing the unresolved Vt2/entrainment in Rib
  call get_param(paramFile, mdl, "USE_KPP_LT_VT2", CS%LT_Vt2_Enhancement, &
       'Flag for Langmuir turbulence enhancement of Vt2'//&
       'in Bulk Richardson Number.', Default=.false.)
  if (CS%LT_Vt2_Enhancement) then
    call get_param(paramFile, mdl, "KPP_LT_VT2_METHOD",string ,                  &
                   'Method to enhance Vt2 in KPP. '//                            &
                   'Valid options are: \n'//                                     &
                   '\t CONSTANT = Constant value (KPP_VT2_ENH_FAC) \n'//         &
                   '\t VR12     = Function of Langmuir number based on VR12\n'// &
                   '\t            (Van Roekel et al., 2012) \n'//                &
                   '\t            (Li et al. 2016, OM) \n'//                     &
                   '\t RW16     = Function of Langmuir number based on RW16\n'// &
                   '\t            (Reichl et al., 2016, JPO) \n'//               &
                   '\t LF17     = Function of Langmuir number based on LF17\n'// &
                   '\t            (Li and Fox-Kemper, 2017, JPO)',    &
                   default='CONSTANT')
    select case ( trim(string))
      case ("CONSTANT")
        CS%LT_VT2_METHOD = LT_VT2_MODE_CONSTANT
        langmuir_entrainment_opt = 'LWF16'
      case ("VR12")
        CS%LT_VT2_METHOD = LT_VT2_MODE_VR12
        langmuir_entrainment_opt = 'LWF16'
      case ("RW16")
        CS%LT_VT2_METHOD = LT_VT2_MODE_RW16
        langmuir_entrainment_opt = 'RWHGK16'
      case ("LF17")
        CS%LT_VT2_METHOD = LT_VT2_MODE_LF17
        langmuir_entrainment_opt = 'LF17'
      case default
        call MOM_error(FATAL,"KPP_init: "//&
          "Unrecognized KPP_LT_VT2_METHOD option: "//trim(string))
    end select
    if (CS%LT_VT2_METHOD==LT_VT2_MODE_CONSTANT) then
      call get_param(paramFile, mdl, "KPP_VT2_ENH_FAC", CS%KPP_VT2_ENH_FAC,     &
                   'Constant value to enhance VT2 in KPP.',  &
                   units="nondim", default=1.0)
    endif
  endif

  if (CS%LT_K_ENHANCEMENT .or. CS%LT_VT2_ENHANCEMENT) then
    call get_param(paramFile, mdl, "KPP_LT_MLD_GUESS_MIN", CS%MLD_guess_min,     &
                   "The minimum estimate of the mixed layer depth used to calculate "//&
                   "the Langmuir number for Langmuir turbulence enhancement with KPP.", &
                   units="m", default=1.0, scale=US%m_to_Z)
  endif

  call get_param(paramFile, mdl, "KPP_CVt2", CS%KPP_CVt2, &
                 'Parameter for Stokes MOST convection entrainment', &
                 units="nondim", default=1.6)

  call get_param(paramFile, mdl, "ANSWER_DATE", CS%answer_date, &
                 "The vintage of the order of arithmetic in the CVMix KPP calculations.  Values "//&
                 "below 20240501 recover the answers from early in 2024, while higher values "//&
                 "use expressions that have been refactored for rotational symmetry.", &
                 default=default_answer_date)

  call closeParameterBlock(paramFile)

  call get_param(paramFile, mdl, 'DEBUG', CS%debug, default=.False., do_not_log=.True.)

  call CVMix_init_kpp( Ri_crit=CS%Ri_crit,                 &
                       minOBLdepth=US%Z_to_m*CS%minOBLdepth, &
                       minVtsqr=US%L_T_to_m_s**2*CS%minVtsqr, &
                       vonKarman=CS%vonKarman,             &
                       surf_layer_ext=CS%surf_layer_ext,   &
                       CVt2=CS%KPP_CVt2,                   &
                       interp_type=CS%interpType,          &
                       interp_type2=CS%interpType2,        &
                       lEkman=CS%computeEkman,             &
                       lStokesMOST=CS%StokesMOST,          &
                       lMonOb=CS%computeMoninObukhov,      &
                       MatchTechnique=CS%MatchTechnique,   &
                       lenhanced_diff=CS%enhance_diffusion,&
                       lnonzero_surf_nonlocal=Cs_is_one   ,&
                       lnoDGat1=lnoDGat1                  ,&
                       langmuir_mixing_str=langmuir_mixing_opt,&
                       langmuir_entrainment_str=langmuir_entrainment_opt,&
                       CVMix_kpp_params_user=CS%KPP_params )

  ! Register diagnostics
  CS%diag => diag
  CS%id_OBLdepth = register_diag_field('ocean_model', 'KPP_OBLdepth', diag%axesT1, Time, &
      'Thickness of the surface Ocean Boundary Layer calculated by [CVMix] KPP', &
      'meter', conversion=US%Z_to_m, &
      cmor_field_name='oml', cmor_long_name='ocean_mixed_layer_thickness_defined_by_mixing_scheme', &
      cmor_units='m', cmor_standard_name='Ocean Mixed Layer Thickness Defined by Mixing Scheme')
      ! CMOR names are placeholders; must be modified by time period
      ! for CMOR compliance. Diag manager will be used for omlmax and
      ! omldamax.
  if (CS%n_smooth > 0) then
    CS%id_OBLdepth_original = register_diag_field('ocean_model', 'KPP_OBLdepth_original', diag%axesT1, Time, &
        'Thickness of the surface Ocean Boundary Layer without smoothing calculated by [CVMix] KPP', &
        'meter', conversion=US%Z_to_m, &
        cmor_field_name='oml', cmor_long_name='ocean_mixed_layer_thickness_defined_by_mixing_scheme', &
        cmor_units='m', cmor_standard_name='Ocean Mixed Layer Thickness Defined by Mixing Scheme')
  endif
  if( CS%StokesMOST ) then
  CS%id_StokesXI = register_diag_field('ocean_model', 'StokesXI', diag%axesT1, Time, &
      'Stokes Similarity Parameter', 'nondim')
  CS%id_Lam2     = register_diag_field('ocean_model', 'Lam2',  diag%axesT1, Time, &
      'Ustk0_ustar', 'nondim')
  endif
  CS%id_BulkDrho = register_diag_field('ocean_model', 'KPP_BulkDrho', diag%axesTL, Time, &
      'Bulk difference in density used in Bulk Richardson number, as used by [CVMix] KPP', &
      'kg/m3', conversion=US%R_to_kg_m3)
  CS%id_BulkUz2 = register_diag_field('ocean_model', 'KPP_BulkUz2', diag%axesTL, Time, &
      'Square of bulk difference in resolved velocity used in Bulk Richardson number via [CVMix] KPP', &
      'm2/s2', conversion=US%L_T_to_m_s**2)
  CS%id_BulkRi = register_diag_field('ocean_model', 'KPP_BulkRi', diag%axesTL, Time, &
      'Bulk Richardson number used to find the OBL depth used by [CVMix] KPP', 'nondim')
  CS%id_Sigma = register_diag_field('ocean_model', 'KPP_sigma', diag%axesTi, Time, &
      'Sigma coordinate used by [CVMix] KPP', 'nondim')
  CS%id_Ws = register_diag_field('ocean_model', 'KPP_Ws', diag%axesTL, Time, &
      'Turbulent vertical velocity scale for scalars used by [CVMix] KPP', &
      'm/s', conversion=US%Z_to_m*US%s_to_T)
  CS%id_N = register_diag_field('ocean_model', 'KPP_N', diag%axesTi, Time, &
      '(Adjusted) Brunt-Vaisala frequency used by [CVMix] KPP', '1/s', conversion=US%s_to_T)
  CS%id_N2 = register_diag_field('ocean_model', 'KPP_N2', diag%axesTi, Time, &
      'Square of Brunt-Vaisala frequency used by [CVMix] KPP', '1/s2', conversion=US%s_to_T**2)
  CS%id_Vt2 = register_diag_field('ocean_model', 'KPP_Vt2', diag%axesTL, Time, &
      'Unresolved shear turbulence used by [CVMix] KPP', 'm2/s2', conversion=US%Z_to_m**2*US%s_to_T**2)
  CS%id_uStar = register_diag_field('ocean_model', 'KPP_uStar', diag%axesT1, Time, &
      'Friction velocity, u*, as used by [CVMix] KPP', 'm/s', conversion=US%Z_to_m*US%s_to_T)
  CS%id_buoyFlux = register_diag_field('ocean_model', 'KPP_buoyFlux', diag%axesTi, Time, &
      'Surface (and penetrating) buoyancy flux, as used by [CVMix] KPP', &
      'm2/s3', conversion=US%L_to_m**2*US%s_to_T**3)
  CS%id_Kt_KPP = register_diag_field('ocean_model', 'KPP_Kheat', diag%axesTi, Time, &
      'Heat diffusivity due to KPP, as calculated by [CVMix] KPP', &
      'm2/s', conversion=US%Z2_T_to_m2_s)
  CS%id_Kd_in = register_diag_field('ocean_model', 'KPP_Kd_in', diag%axesTi, Time, &
      'Diffusivity passed to KPP', 'm2/s', conversion=GV%HZ_T_to_m2_s)
  CS%id_Ks_KPP = register_diag_field('ocean_model', 'KPP_Ksalt', diag%axesTi, Time, &
      'Salt diffusivity due to KPP, as calculated by [CVMix] KPP', &
      'm2/s', conversion=US%Z2_T_to_m2_s)
  CS%id_Kv_KPP = register_diag_field('ocean_model', 'KPP_Kv', diag%axesTi, Time, &
      'Vertical viscosity due to KPP, as calculated by [CVMix] KPP', &
      'm2/s', conversion=US%Z2_T_to_m2_s)
  CS%id_NLTt = register_diag_field('ocean_model', 'KPP_NLtransport_heat', diag%axesTi, Time, &
      'Non-local transport (Cs*G(sigma)) for heat, as calculated by [CVMix] KPP', 'nondim')
  CS%id_NLTs = register_diag_field('ocean_model', 'KPP_NLtransport_salt', diag%axesTi, Time, &
      'Non-local tranpsort (Cs*G(sigma)) for scalars, as calculated by [CVMix] KPP', 'nondim')
  CS%id_Tsurf = register_diag_field('ocean_model', 'KPP_Tsurf', diag%axesT1, Time, &
      'Temperature of surface layer (10% of OBL depth) as passed to [CVMix] KPP', &
      'C', conversion=US%C_to_degC)
  CS%id_Ssurf = register_diag_field('ocean_model', 'KPP_Ssurf', diag%axesT1, Time, &
      'Salinity of surface layer (10% of OBL depth) as passed to [CVMix] KPP', &
      'ppt', conversion=US%S_to_ppt)
  CS%id_Usurf = register_diag_field('ocean_model', 'KPP_Usurf', diag%axesCu1, Time, &
      'i-component flow of surface layer (10% of OBL depth) as passed to [CVMix] KPP', &
      'm/s', conversion=US%L_T_to_m_s)
  CS%id_Vsurf = register_diag_field('ocean_model', 'KPP_Vsurf', diag%axesCv1, Time, &
      'j-component flow of surface layer (10% of OBL depth) as passed to [CVMix] KPP', &
      'm/s', conversion=US%L_T_to_m_s)
  CS%id_EnhK = register_diag_field('ocean_model', 'EnhK', diag%axesTI, Time, &
      'Langmuir number enhancement to K as used by [CVMix] KPP','nondim')
  CS%id_EnhVt2 = register_diag_field('ocean_model', 'EnhVt2', diag%axesTL, Time, &
      'Langmuir number enhancement to Vt2 as used by [CVMix] KPP','nondim')
  CS%id_La_SL = register_diag_field('ocean_model', 'KPP_La_SL', diag%axesT1, Time, &
      'Surface-layer Langmuir number computed in [CVMix] KPP','nondim')

  allocate( CS%N( SZI_(G), SZJ_(G), SZK_(GV)+1 ), source=0. )
  allocate( CS%StokesParXI( SZI_(G), SZJ_(G) ), source=0. )
  allocate( CS%Lam2    ( SZI_(G), SZJ_(G) ), source=0. )
  allocate( CS%kOBL( SZI_(G), SZJ_(G) ), source=0. )
  allocate( CS%La_SL( SZI_(G), SZJ_(G) ), source=0. )
  allocate( CS%Vt2( SZI_(G), SZJ_(G), SZK_(GV) ), source=0. )
  if (CS%id_OBLdepth_original > 0) allocate( CS%OBLdepth_original( SZI_(G), SZJ_(G) ) )

  allocate( CS%OBLdepthprev( SZI_(G), SZJ_(G) ), source=0.0 )
  if (CS%id_BulkDrho > 0) allocate( CS%dRho( SZI_(G), SZJ_(G), SZK_(GV) ), source=0. )
  if (CS%id_BulkUz2 > 0)  allocate( CS%Uz2( SZI_(G), SZJ_(G), SZK_(GV) ), source=0. )
  if (CS%id_BulkRi > 0)   allocate( CS%BulkRi( SZI_(G), SZJ_(G), SZK_(GV) ), source=0. )
  if (CS%id_Sigma > 0)    allocate( CS%sigma( SZI_(G), SZJ_(G), SZK_(GV)+1 ), source=0. )
  if (CS%id_Ws > 0)       allocate( CS%Ws( SZI_(G), SZJ_(G), SZK_(GV) ), source=0. )
  if (CS%id_N2 > 0)       allocate( CS%N2( SZI_(G), SZJ_(G), SZK_(GV)+1 ), source=0. )
  if (CS%id_Kt_KPP > 0)   allocate( CS%Kt_KPP( SZI_(G), SZJ_(G), SZK_(GV)+1 ), source=0. )
  if (CS%id_Ks_KPP > 0)   allocate( CS%Ks_KPP( SZI_(G), SZJ_(G), SZK_(GV)+1 ), source=0. )
  if (CS%id_Kv_KPP > 0)   allocate( CS%Kv_KPP( SZI_(G), SZJ_(G), SZK_(GV)+1 ), source=0. )
  if (CS%id_Tsurf > 0)    allocate( CS%Tsurf( SZI_(G), SZJ_(G) ), source=0. )
  if (CS%id_Ssurf > 0)    allocate( CS%Ssurf( SZI_(G), SZJ_(G) ), source=0. )
  if (CS%id_Usurf > 0)    allocate( CS%Usurf( SZIB_(G), SZJ_(G) ), source=0. )
  if (CS%id_Vsurf > 0)    allocate( CS%Vsurf( SZI_(G), SZJB_(G) ), source=0. )
  if (CS%id_EnhVt2 > 0)   allocate( CS%EnhVt2( SZI_(G), SZJ_(G), SZK_(GV) ), source=0. )
  if (CS%id_EnhK > 0)     allocate( CS%EnhK( SZI_(G), SZJ_(G), SZK_(GV)+1 ), source=0. )

  id_clock_KPP_calc = cpu_clock_id('Ocean KPP calculate)', grain=CLOCK_MODULE)
  id_clock_KPP_compute_BLD = cpu_clock_id('(Ocean KPP comp BLD)', grain=CLOCK_ROUTINE)

end function KPP_init

!> KPP vertical diffusivity/viscosity and non-local tracer transport
subroutine KPP_calculate(CS, G, GV, US, h, tv, uStar, buoyFlux, Kt, Ks, Kv, &
                         nonLocalTransHeat, nonLocalTransScalar, Waves, lamult)

  ! Arguments
  type(KPP_CS),                                pointer       :: CS    !< Control structure
  type(ocean_grid_type),                       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type),                     intent(in)    :: GV    !< Ocean vertical grid
  type(unit_scale_type),                       intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                       intent(in)    :: tv    !< Thermodynamics structure.
  real, dimension(SZI_(G),SZJ_(G)),            intent(in)    :: uStar !< Surface friction velocity [Z T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: buoyFlux !< Surface buoyancy flux [L2 T-3 ~> m2 s-3]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kt  !< (in)  Vertical diffusivity of heat w/o KPP
                                                                    !! (out) Vertical diffusivity including KPP
                                                                    !!       [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Ks  !< (in)  Vertical diffusivity of salt w/o KPP
                                                                    !! (out) Vertical diffusivity including KPP
                                                                    !!       [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kv  !< (in)  Vertical viscosity w/o KPP
                                                                    !! (out) Vertical viscosity including KPP
                                                                    !!       [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: nonLocalTransHeat   !< Temp non-local transport [nondim]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: nonLocalTransScalar !< scalar non-local trans. [nondim]
  type(wave_parameters_CS),                    pointer       :: Waves   !< Wave CS for Langmuir turbulence
  real, dimension(SZI_(G),SZJ_(G)),  optional, intent(in)    :: lamult  !< Langmuir enhancement multiplier [nondim]

  ! Local variables
  integer :: i, j, k                            ! Loop indices
  real, dimension(SZI_(G),SZK_(GV)) :: dz       ! Height change across layers [Z ~> m]
  real, dimension( GV%ke )     :: cellHeight    ! Cell center heights referenced to surface [Z ~> m] (negative in ocean)
  real, dimension( GV%ke+1 )   :: iFaceHeight   ! Interface heights referenced to surface [Z ~> m] (negative in ocean)
  real, dimension( GV%ke )     :: z_cell        ! Cell center heights referenced to surface [m] (negative in ocean)
  real, dimension( GV%ke+1 )   :: z_inter       ! Cell interface heights referenced to surface [m] (negative in ocean)
  real, dimension( GV%ke+1, 2) :: Kdiffusivity  ! Vertical diffusivity at interfaces in MKS units [m2 s-1]
  real, dimension( GV%ke+1 )   :: Kviscosity    ! Vertical viscosity at interfaces in MKS units [m2 s-1]
  real, dimension( GV%ke+1, 2) :: nonLocalTrans ! Non-local transport for heat/salt at interfaces [nondim]

  real :: surfFricVel   ! Surface friction velocity in MKS units [m s-1]
  real :: surfBuoyFlux  ! Surface buoyancy flux in MKS units [m2 s-3]
  real :: sigma      ! Fractional vertical position within the boundary layer [nondim]
  real :: sigmaRatio ! A cubic function of sigma [nondim]
  real :: buoy_scale ! A unit conversion factor for buoyancy fluxes [m2 T3 L-2 s-3 ~> 1]
  real :: dh    ! The local thickness used for calculating interface positions [Z ~> m]
  real :: hcorr ! A cumulative correction arising from inflation of vanished layers [Z ~> m]

  ! For Langmuir Calculations
  real :: LangEnhK     ! Langmuir enhancement for mixing coefficient [nondim]

  if (CS%Stokes_Mixing .and. .not.associated(Waves)) call MOM_error(FATAL, &
      "KPP_calculate: The Waves control structure must be associated if STOKES_MIXING is True.")

  if (CS%debug) then
    call hchksum(h, "KPP in: h", G%HI, haloshift=0, unscale=GV%H_to_m)
    call hchksum(uStar, "KPP in: uStar", G%HI, haloshift=0, unscale=US%Z_to_m*US%s_to_T)
    call hchksum(buoyFlux, "KPP in: buoyFlux", G%HI, haloshift=0, unscale=US%L_to_m**2*US%s_to_T**3)
    call hchksum(Kt, "KPP in: Kt", G%HI, haloshift=0, unscale=GV%HZ_T_to_m2_s)
    call hchksum(Ks, "KPP in: Ks", G%HI, haloshift=0, unscale=GV%HZ_T_to_m2_s)
  endif

  nonLocalTrans(:,:) = 0.0

  if (CS%id_Kd_in > 0) call post_data(CS%id_Kd_in, Kt, CS%diag)

  call cpu_clock_begin(id_clock_KPP_calc)
  buoy_scale = US%L_to_m**2*US%s_to_T**3

  !$OMP parallel do default(none) firstprivate(nonLocalTrans)                               &
  !$OMP                           private(surfFricVel, iFaceHeight, hcorr, dh, dz, cellHeight,  &
  !$OMP                           surfBuoyFlux, Kdiffusivity, Kviscosity, LangEnhK, sigma,  &
  !$OMP                           sigmaRatio, z_inter, z_cell)                              &
  !$OMP                           shared(G, GV, CS, US, tv, uStar, h, buoy_scale, buoyFlux, Kt, &
  !$OMP                           Ks, Kv, nonLocalTransHeat, nonLocalTransScalar, Waves, lamult)
  ! loop over horizontal points on processor
  do j = G%jsc, G%jec

    ! Find the vertical distances across layers.
    call thickness_to_dz(h, tv, dz, j, G, GV)

    do i = G%isc, G%iec ; if (G%mask2dT(i,j) > 0.0) then

      ! things independent of position within the column
      surfFricVel = US%Z_to_m*US%s_to_T * uStar(i,j)

      iFaceHeight(1) = 0.0 ! BBL is all relative to the surface
      hcorr = 0.
      do k=1,GV%ke

        ! cell center and cell bottom in meters (negative values in the ocean)
        dh = dz(i,k)    ! Nominal thickness to use for increment
        dh = dh + hcorr ! Take away the accumulated error (could temporarily make dh<0)
        hcorr = min( dh - CS%min_thickness, 0. ) ! If inflating then hcorr<0
        dh = max( dh, CS%min_thickness ) ! Limit increment dh>=min_thickness
        cellHeight(k)    = iFaceHeight(k) - 0.5 * dh
        iFaceHeight(k+1) = iFaceHeight(k) - dh

      enddo ! k-loop finishes

      surfBuoyFlux = buoy_scale*buoyFlux(i,j,1) ! This is only used in kpp_compute_OBL_depth to limit
                                     ! h to Monin-Obukhov (default is false, ie. not used)

      ! Call CVMix/KPP to obtain OBL diffusivities, viscosities and non-local transports

      ! Unlike LMD94, we do not match to interior diffusivities. If using the original
      ! LMD94 shape function, not matching is equivalent to matching to a zero diffusivity.

      !BGR/ Add option for use of surface buoyancy flux with total sw flux.
      if (CS%SW_METHOD == SW_METHOD_ALL_SW) then
         surfBuoyFlux = buoy_scale * buoyFlux(i,j,1)
      elseif (CS%SW_METHOD == SW_METHOD_MXL_SW) then
         ! We know the actual buoyancy flux into the OBL
         surfBuoyFlux  = buoy_scale * (buoyFlux(i,j,1) - buoyFlux(i,j,int(CS%kOBL(i,j))+1))
      elseif (CS%SW_METHOD == SW_METHOD_LV1_SW) then
         surfBuoyFlux  = buoy_scale * (buoyFlux(i,j,1) - buoyFlux(i,j,2))
      endif

      ! If option "MatchBoth" is selected in CVMix, MOM should be capable of matching.
      if (.not. (CS%MatchTechnique == 'MatchBoth')) then
         Kdiffusivity(:,:) = 0. ! Diffusivities for heat and salt [m2 s-1]
         Kviscosity(:)     = 0. ! Viscosity [m2 s-1]
      else
         Kdiffusivity(:,1) = GV%HZ_T_to_m2_s * Kt(i,j,:)
         Kdiffusivity(:,2) = GV%HZ_T_to_m2_s * Ks(i,j,:)
         Kviscosity(:) = GV%HZ_T_to_m2_s * Kv(i,j,:)
      endif

      IF (CS%LT_K_ENHANCEMENT) then
        if (CS%LT_K_METHOD==LT_K_MODE_CONSTANT) then
           LangEnhK = CS%KPP_K_ENH_FAC
        elseif (CS%LT_K_METHOD==LT_K_MODE_VR12) then
          if (present(lamult)) then
            LangEnhK = lamult(i,j)
          else
            LangEnhK = sqrt(1.+(1.5*CS%La_SL(i,j))**(-2) + &
                (5.4*CS%La_SL(i,j))**(-4))
          endif
        elseif (CS%LT_K_METHOD==LT_K_MODE_RW16) then
          !This maximum value is proposed in Reichl et al., 2016 JPO formula
          LangEnhK = min(2.25, 1. + 1./CS%La_SL(i,j))
        else
           !This shouldn't be reached.
           !call MOM_error(WARNING,"Unexpected behavior in MOM_CVMix_KPP, see error in LT_K_ENHANCEMENT")
           LangEnhK = 1.0
        endif

        ! diffusivities don't need to be enhanced below anymore since LangEnhK is applied within CVMix.
        ! todo: need to double check if the distinction between the two different options of LT_K_SHAPE may need to be
        ! treated specially.
        do k=1,GV%ke
          if (CS%LT_K_SHAPE== LT_K_CONSTANT) then
            if (CS%id_EnhK > 0) CS%EnhK(i,j,:) = LangEnhK
            !Kdiffusivity(k,1) = Kdiffusivity(k,1) * LangEnhK
            !Kdiffusivity(k,2) = Kdiffusivity(k,2) * LangEnhK
            !Kviscosity(k)     = Kviscosity(k)   * LangEnhK
          elseif (CS%LT_K_SHAPE == LT_K_SCALED) then
            sigma = min(1.0,-iFaceHeight(k)/CS%OBLdepth(i,j))
            SigmaRatio = sigma * (1. - sigma)**2 / 0.148148037
            if (CS%id_EnhK > 0) CS%EnhK(i,j,k) = (1.0 + (LangEnhK - 1.)*sigmaRatio)
            !Kdiffusivity(k,1) = Kdiffusivity(k,1) * ( 1. + &
            !                    ( LangEnhK - 1.)*sigmaRatio)
            !Kdiffusivity(k,2) = Kdiffusivity(k,2) * ( 1. + &
            !                    ( LangEnhK - 1.)*sigmaRatio)
            !Kviscosity(k) = Kviscosity(k) * ( 1. + &
            !                    ( LangEnhK - 1.)*sigmaRatio)
          endif
        enddo
      endif

      ! Convert columns to MKS units for passing to CVMix
      do k = 1, GV%ke
        z_cell(k) = US%Z_to_m*cellHeight(k)
      enddo
      do K = 1, GV%ke+1
        z_inter(K) = US%Z_to_m*iFaceHeight(K)
      enddo

      call CVMix_coeffs_kpp(Kviscosity(:),     & ! (inout) Total viscosity [m2 s-1]
                            Kdiffusivity(:,1), & ! (inout) Total heat diffusivity [m2 s-1]
                            Kdiffusivity(:,2), & ! (inout) Total salt diffusivity [m2 s-1]
                            z_inter(:),        & ! (in) Height of interfaces [m]
                            z_cell(:),         & ! (in) Height of level centers [m]
                            Kviscosity(:),     & ! (in) Original viscosity [m2 s-1]
                            Kdiffusivity(:,1), & ! (in) Original heat diffusivity [m2 s-1]
                            Kdiffusivity(:,2), & ! (in) Original salt diffusivity [m2 s-1]
                            US%Z_to_m*CS%OBLdepth(i,j),  & ! (in) OBL depth [m]
                            CS%kOBL(i,j),      & ! (in) level (+fraction) of OBL extent
                            nonLocalTrans(:,1),& ! (out) Non-local heat transport [nondim]
                            nonLocalTrans(:,2),& ! (out) Non-local salt transport [nondim]
                            surfFricVel,       & ! (in) Turbulent friction velocity at surface [m s-1]
                            surfBuoyFlux,      & ! (in) Buoyancy flux at surface [m2 s-3]
                            GV%ke,             & ! (in) Number of levels to compute coeffs for
                            GV%ke,             & ! (in) Number of levels in array shape
                            Langmuir_EFactor=LangEnhK,& ! Langmuir enhancement multiplier
                            StokesXi = CS%StokesParXI(i,j), & ! Stokes forcing parameter
                            CVMix_kpp_params_user=CS%KPP_params )

      ! safety check, Kviscosity and Kdiffusivity must be >= 0
      do k=1, GV%ke+1
        if (Kviscosity(k) < 0. .or. Kdiffusivity(k,1) < 0.) then
          write(*,'(a,3i3)')    'interface, i, j, k = ',j, j, k
          write(*,'(a,2f12.5)') 'lon,lat=', G%geoLonT(i,j), G%geoLatT(i,j)
          write(*,'(a,es12.4)') 'depth, z_inter(k) =',z_inter(k)
          write(*,'(a,es12.4)') 'Kviscosity(k) =',Kviscosity(k)
          write(*,'(a,es12.4)') 'Kdiffusivity(k,1) =',Kdiffusivity(k,1)
          write(*,'(a,es12.4)') 'Kdiffusivity(k,2) =',Kdiffusivity(k,2)
          write(*,'(a,es12.4)') 'OBLdepth =',US%Z_to_m*CS%OBLdepth(i,j)
          write(*,'(a,f8.4)')   'kOBL =',CS%kOBL(i,j)
          write(*,'(a,es12.4)') 'u* =',surfFricVel
          write(*,'(a,es12.4)') 'bottom, z_inter(GV%ke+1) =',z_inter(GV%ke+1)
          write(*,'(a,es12.4)') 'CS%La_SL(i,j) =',CS%La_SL(i,j)
          write(*,'(a,es12.4)') 'LangEnhK =',LangEnhK
          if (present(lamult)) write(*,'(a,es12.4)') 'lamult(i,j) =',lamult(i,j)
          write(*,*) 'Kviscosity(:) =',Kviscosity(:)
          write(*,*) 'Kdiffusivity(:,1) =',Kdiffusivity(:,1)

          call MOM_error(FATAL,"KPP_calculate, after CVMix_coeffs_kpp: "// &
                   "Negative vertical viscosity or diffusivity has been detected. " // &
                   "This is likely related to the choice of MATCH_TECHNIQUE and INTERP_TYPE2." //&
                   "You might consider using the default options for these parameters." )
        endif
      enddo

      ! Over-write CVMix NLT shape function with one of the following choices.
      ! The CVMix code has yet to update for thse options, so we compute in MOM6.
      ! Note that nonLocalTrans = Cs * G(sigma) (LMD94 notation), with
      ! Cs = 6.32739901508.
      ! Start do-loop at k=2, since k=1 is ocean surface (sigma=0)
      ! and we do not wish to double-count the surface forcing.
      ! Only compute nonlocal transport for 0 <= sigma <= 1.
      ! MOM6 recommended shape is the parabolic; it gives deeper boundary layer
      ! and no spurious extrema.
      if (surfBuoyFlux < 0.0) then
        if (CS%NLT_shape == NLT_SHAPE_CUBIC) then
          do k = 2, GV%ke
            sigma = min(1.0,-iFaceHeight(k)/CS%OBLdepth(i,j))
            nonLocalTrans(k,1) = (1.0 - sigma)**2 * (1.0 + 2.0*sigma) !*
            nonLocalTrans(k,2) = nonLocalTrans(k,1)
          enddo
        elseif (CS%NLT_shape == NLT_SHAPE_PARABOLIC) then
          do k = 2, GV%ke
            sigma = min(1.0,-iFaceHeight(k)/CS%OBLdepth(i,j))
            nonLocalTrans(k,1) = (1.0 - sigma)**2 !*CS%CS2
            nonLocalTrans(k,2) = nonLocalTrans(k,1)
          enddo
        elseif (CS%NLT_shape == NLT_SHAPE_LINEAR) then
          do k = 2, GV%ke
            sigma = min(1.0,-iFaceHeight(k)/CS%OBLdepth(i,j))
            nonLocalTrans(k,1) = (1.0 - sigma)!*CS%CS2
            nonLocalTrans(k,2) = nonLocalTrans(k,1)
          enddo
        elseif (CS%NLT_shape == NLT_SHAPE_CUBIC_LMD) then
          ! Sanity check (should agree with CVMix result using simple matching)
          do k = 2, GV%ke
            sigma = min(1.0,-iFaceHeight(k)/CS%OBLdepth(i,j))
            nonLocalTrans(k,1) = CS%CS2 * sigma*(1.0 -sigma)**2
            nonLocalTrans(k,2) = nonLocalTrans(k,1)
          enddo
        endif
      endif

      ! we apply nonLocalTrans in subroutines
      ! KPP_NonLocalTransport_temp and KPP_NonLocalTransport_saln
      nonLocalTransHeat(i,j,:)   = nonLocalTrans(:,1) ! temperature
      nonLocalTransScalar(i,j,:) = nonLocalTrans(:,2) ! salinity

      ! set the KPP diffusivity and viscosity to zero for testing purposes
      if (CS%KPPzeroDiffusivity) then
         Kdiffusivity(:,1) = 0.0
         Kdiffusivity(:,2) = 0.0
         Kviscosity(:)     = 0.0
      endif

      ! Copy 1d data into 3d diagnostic arrays
      !/ grabbing obldepth_0d for next time step.
      CS%OBLdepthprev(i,j) = CS%OBLdepth(i,j)
      if (CS%id_sigma > 0) then
        CS%sigma(i,j,:)  = 0.
        if (CS%OBLdepth(i,j)>0.) CS%sigma(i,j,:)  = -iFaceHeight(:)/CS%OBLdepth(i,j)
      endif
      if (CS%id_Kt_KPP > 0) CS%Kt_KPP(i,j,:) = US%m2_s_to_Z2_T * Kdiffusivity(:,1)
      if (CS%id_Ks_KPP > 0) CS%Ks_KPP(i,j,:) = US%m2_s_to_Z2_T * Kdiffusivity(:,2)
      if (CS%id_Kv_KPP > 0) CS%Kv_KPP(i,j,:) = US%m2_s_to_Z2_T * Kviscosity(:)

      ! Update output of routine
      if (.not. CS%passiveMode) then
        if (CS%KPPisAdditive) then
          do k=1, GV%ke+1
            Kt(i,j,k) = Kt(i,j,k) + GV%m2_s_to_HZ_T * Kdiffusivity(k,1)
            Ks(i,j,k) = Ks(i,j,k) + GV%m2_s_to_HZ_T * Kdiffusivity(k,2)
            Kv(i,j,k) = Kv(i,j,k) + GV%m2_s_to_HZ_T * Kviscosity(k)
            if (CS%Stokes_Mixing) Waves%KvS(i,j,k) = Kv(i,j,k)
          enddo
        else ! KPP replaces prior diffusivity when former is non-zero
          do k=1, GV%ke+1
            if (Kdiffusivity(k,1) /= 0.) Kt(i,j,k) = GV%m2_s_to_HZ_T * Kdiffusivity(k,1)
            if (Kdiffusivity(k,2) /= 0.) Ks(i,j,k) = GV%m2_s_to_HZ_T * Kdiffusivity(k,2)
            if (Kviscosity(k) /= 0.) Kv(i,j,k) = GV%m2_s_to_HZ_T * Kviscosity(k)
            if (CS%Stokes_Mixing) Waves%KvS(i,j,k) = Kv(i,j,k)
          enddo
        endif
      endif


    ! end of the horizontal do-loops over the vertical columns
    endif ; enddo ! i
  enddo ! j

  call cpu_clock_end(id_clock_KPP_calc)

  if (CS%debug) then
    call hchksum(Kt, "KPP out: Kt", G%HI, haloshift=0, unscale=GV%HZ_T_to_m2_s)
    call hchksum(Ks, "KPP out: Ks", G%HI, haloshift=0, unscale=GV%HZ_T_to_m2_s)
  endif

  ! send diagnostics to post_data
  if (CS%id_OBLdepth > 0) call post_data(CS%id_OBLdepth, CS%OBLdepth,        CS%diag)
  if (CS%id_OBLdepth_original > 0) call post_data(CS%id_OBLdepth_original,CS%OBLdepth_original,CS%diag)
  if (CS%id_sigma    > 0) call post_data(CS%id_sigma,    CS%sigma,           CS%diag)
  if (CS%id_Ws       > 0) call post_data(CS%id_Ws,       CS%Ws,              CS%diag)
  if (CS%id_uStar    > 0) call post_data(CS%id_uStar,    uStar,              CS%diag)
  if (CS%id_buoyFlux > 0) call post_data(CS%id_buoyFlux, buoyFlux,           CS%diag)
  if (CS%id_Kt_KPP   > 0) call post_data(CS%id_Kt_KPP,   CS%Kt_KPP,          CS%diag)
  if (CS%id_Ks_KPP   > 0) call post_data(CS%id_Ks_KPP,   CS%Ks_KPP,          CS%diag)
  if (CS%id_Kv_KPP   > 0) call post_data(CS%id_Kv_KPP,   CS%Kv_KPP,          CS%diag)
  if (CS%id_NLTt     > 0) call post_data(CS%id_NLTt,     nonLocalTransHeat,  CS%diag)
  if (CS%id_NLTs     > 0) call post_data(CS%id_NLTs,     nonLocalTransScalar,CS%diag)


end subroutine KPP_calculate


!> Compute OBL depth
subroutine KPP_compute_BLD(CS, G, GV, US, h, Temp, Salt, u, v, tv, uStar, buoyFlux, Waves, lamult)

  ! Arguments
  type(KPP_CS),                               pointer       :: CS    !< Control structure
  type(ocean_grid_type),                      intent(inout) :: G     !< Ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV    !< Ocean vertical grid
  type(unit_scale_type),                      intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: Temp  !< potential/cons temp [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: Salt  !< Salinity [S ~> ppt]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u     !< Velocity i-component [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v     !< Velocity j-component [L T-1 ~> m s-1]
  type(thermo_var_ptrs),                      intent(in)    :: tv    !< Thermodynamics structure.
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: uStar !< Surface friction velocity [Z T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)   :: buoyFlux !< Surface buoyancy flux [L2 T-3 ~> m2 s-3]
  type(wave_parameters_CS),                   pointer       :: Waves !< Wave CS for Langmuir turbulence
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)    :: lamult !< Langmuir enhancement factor [nondim]

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) ::  dz  ! Height change across layers [Z ~> m]
  ! Variables for passing to CVMix routines, often in MKS units
  real, dimension( GV%ke )   :: Ws_1d          ! Profile of vertical velocity scale for scalars in MKS units [m s-1]
  real, dimension( GV%ke )   :: deltaRho       ! delta Rho in numerator of Bulk Ri number [R ~> kg m-3]
  real, dimension( GV%ke )   :: deltaBuoy      ! Change in Buoyancy based on deltaRho [m s-2]
  real, dimension( GV%ke )   :: deltaU2        ! square of delta U (shear) in denominator of Bulk Ri [m2 s-2]
  real, dimension( GV%ke )   :: surfBuoyFlux2  ! Surface buoyancy flux in MKS units [m2 s-3]
  real, dimension( GV%ke )   :: BulkRi_1d      ! Bulk Richardson number for each layer [nondim]
  real, dimension( GV%ke )   :: Vt2_1d         ! Unresolved squared turbulence velocity for bulk Ri [m2 s-2]
  real, dimension( GV%ke )   :: z_cell         ! Cell center heights referenced to surface [m] (negative in ocean)
  real, dimension( GV%ke )   :: OBL_depth      ! Cell center depths referenced to surface [m] (positive in ocean)
  real, dimension( GV%ke+1 ) :: z_inter        ! Cell interface heights referenced to surface [m] (negative in ocean)
  real, dimension( GV%ke+1 ) :: N_col          ! A column of buoyancy frequencies at interfaces in MKS units [s-1]
  real :: surfFricVel           ! Surface friction velocity in MKS units [m s-1]
  real :: surfBuoyFlux          ! Surface buoyancy flux in MKS units [m2 s-3]
  real :: Coriolis              ! Coriolis parameter at tracer points in MKS units [s-1]
  real :: KPP_OBL_depth         ! Boundary layer depth calculated by CVMix_kpp_compute_OBL_depth in MKS units [m]

  ! Variables for EOS calculations
  real, dimension( 3*GV%ke )   :: rho_1D   ! A column of densities [R ~> kg m-3]
  real, dimension( 3*GV%ke )   :: pres_1D  ! A column of pressures [R L2 T-2 ~> Pa]
  real, dimension( 3*GV%ke )   :: Temp_1D  ! A column of temperatures [C ~> degC]
  real, dimension( 3*GV%ke )   :: Salt_1D  ! A column of salinities [S ~> ppt]

  real, dimension( GV%ke )     :: cellHeight   ! Cell center heights referenced to surface [Z ~> m] (negative in ocean)
  real, dimension( GV%ke+1 )   :: iFaceHeight  ! Interface heights referenced to surface [Z ~> m] (negative in ocean)
  real, dimension( GV%ke+1 )   :: N2_1d        ! Brunt-Vaisala frequency squared, at interfaces [T-2 ~> s-2]
  real :: zBottomMinusOffset    ! Height of bottom plus a little bit [Z ~> m]
  real :: GoRho         ! Gravitational acceleration in MKS units divided by density [m s-2 R-1 ~> m4 kg-1 s-2]
  real :: GoRho_Z_L2    ! Gravitational acceleration, perhaps divided by density, times aspect ratio
                        ! rescaling [H T-2 R-1 ~> m4 kg-1 s-2 or m s-2]
  real :: pRef          ! The interface pressure [R L2 T-2 ~> Pa]
  real :: Uk, Vk        ! Layer velocities relative to their averages in the surface layer [L T-1 ~> m s-1]
  real :: SLdepth_0d    ! Surface layer depth = surf_layer_ext*OBLdepth [Z ~> m]
  real :: hTot          ! Running sum of thickness used in the surface layer average [Z ~> m]
  real :: I_hTot        ! The inverse of hTot [Z-1 ~> m-1]
  real :: buoy_scale    ! A unit conversion factor for buoyancy fluxes [m2 T3 L-2 s-3 ~> 1]
  real :: delH          ! Thickness of a layer [Z ~> m]
  real :: surfTemp      ! Average of temperature over the surface layer [C ~> degC]
  real :: surfHtemp     ! Integral of temperature over the surface layer [Z C ~> m degC]
  real :: surfSalt      ! Average of salinity over the surface layer [S ~> ppt]
  real :: surfHsalt     ! Integral of salinity over the surface layer [Z S ~> m ppt]
  real :: surfHu, surfHv  ! Integral of u and v over the surface layer [Z L T-1 ~> m2 s-1]
  real :: surfU, surfV  ! Average of u and v over the surface layer [Z T-1 ~> m s-1]
  real :: dh            ! The local thickness used for calculating interface positions [Z ~> m]
  real :: hcorr         ! A cumulative correction arising from inflation of vanished layers [Z ~> m]

  ! For Langmuir Calculations
  real :: Vt_layer     ! non-dimensional extent contribution to unresolved shear
  real :: LangEnhVt2   ! Langmuir enhancement for unresolved shear [nondim]
  real, dimension(GV%ke) :: U_H, V_H ! Velocities at tracer points [L T-1 ~> m s-1]
  real :: MLD_guess    ! A guess at the mixed layer depth for calculating the Langmuir number [Z ~> m]
  real :: LA           ! The local Langmuir number [nondim]
  real :: surfHuS, surfHvS ! Stokes drift velocities integrated over the boundary layer [Z L T-1 ~> m2 s-1]
  real :: surfUs, surfVs   ! Stokes drift velocities averaged over the boundary layer [Z T-1 ~> m s-1]

  integer :: i, j, k, km1, kk, ksfc, ktmp    ! Loop indices

  real, dimension(GV%ke)   :: uE_H, vE_H        ! Eulerian velocities h-points, centers [L T-1 ~> m s-1]
  real, dimension(GV%ke)   :: uS_H, vS_H        ! Stokes drift components h-points, centers [L T-1 ~> m s-1]
  real, dimension(GV%ke)   :: uSbar_H, vSbar_H  ! Cell Average Stokes drift h-points [L T-1 ~> m s-1]
  real, dimension(GV%ke+1) :: uS_Hi, vS_Hi      ! Stokes Drift components at interfaces [L T-1 ~> m s-1]
  real :: uS_SLD , vS_SLD, uS_SLC , vS_SLC, uSbar_SLD, vSbar_SLD ! Stokes at/to to Surface Layer Extent
                                                                 ! [L T-1 ~> m s-1]
  real :: StokesXI     ! Stokes similarity parameter [nondim]
  real, dimension( GV%ke )     :: StokesXI_1d ,   StokesVt_1d    !  Parameters of TKE production ratio [nondim]
  real :: Llimit ! Stable boundary Layer Limit =  vonk Lstar [Z ~> m]
  integer :: kbl ! index of cell containing boundary layer depth

  if (CS%Stokes_Mixing .and. .not.associated(Waves)) call MOM_error(FATAL, &
      "KPP_compute_BLD: The Waves control structure must be associated if STOKES_MIXING is True.")

  if (CS%debug) then
    call hchksum(Salt, "KPP in: S", G%HI, haloshift=0, unscale=US%S_to_ppt)
    call hchksum(Temp, "KPP in: T", G%HI, haloshift=0, unscale=US%C_to_degC)
    call hchksum(u, "KPP in: u", G%HI, haloshift=0, unscale=US%L_T_to_m_s)
    call hchksum(v, "KPP in: v", G%HI, haloshift=0, unscale=US%L_T_to_m_s)
  endif

  call cpu_clock_begin(id_clock_KPP_compute_BLD)

  ! some constants
  GoRho = US%Z_to_m*US%s_to_T**2 * (GV%g_Earth_Z_T2 / GV%Rho0)
  if (GV%Boussinesq) then
    GoRho_Z_L2 = GV%Z_to_H * GV%g_Earth_Z_T2 / GV%Rho0
  else
    GoRho_Z_L2 = GV%g_Earth_Z_T2 * GV%RZ_to_H
  endif
  buoy_scale = US%L_to_m**2*US%s_to_T**3

  ! Find the vertical distances across layers.
  call thickness_to_dz(h, tv, dz, G, GV, US)

  ! loop over horizontal points on processor
  !$OMP parallel do default(none) private(surfFricVel, iFaceHeight, hcorr, dh, cellHeight,  &
  !$OMP                           surfBuoyFlux, U_H, V_H, Coriolis, pRef, SLdepth_0d, vt2_1d, &
  !$OMP                           ksfc, surfHtemp, surfHsalt, surfHu, surfHv, surfHuS,      &
  !$OMP                           surfHvS, hTot, I_hTot, delH, surftemp, surfsalt, surfu, surfv,    &
  !$OMP                           surfUs, surfVs, Uk, Vk, deltaU2, km1, kk, pres_1D, N_col, &
  !$OMP                           Temp_1D, salt_1D, surfBuoyFlux2, MLD_guess, LA, rho_1D,   &
  !$OMP                           deltarho, deltaBuoy, N2_1d, ws_1d, LangEnhVT2,KPP_OBL_depth, z_cell, &
  !$OMP                           z_inter, OBL_depth, BulkRi_1d, zBottomMinusOffset, uE_H, vE_H, &
  !$OMP                           uS_H, vS_H, uSbar_H, vSbar_H , uS_Hi, vS_Hi, &
  !$OMP                           uS_SLD, vS_SLD, uS_SLC, vS_SLC, uSbar_SLD, vSbar_SLD, &
  !$OMP                           StokesXI, StokesXI_1d, StokesVt_1d, kbl) &
  !$OMP                           shared(G, GV, CS, US, uStar, h, dz, buoy_scale, buoyFlux, &
  !$OMP                           Temp, Salt, waves, tv, GoRho, GoRho_Z_L2, u, v, lamult, Vt_layer)
  do j = G%jsc, G%jec
    do i = G%isc, G%iec ; if (G%mask2dT(i,j) > 0.0) then

      do k=1,GV%ke
        U_H(k) = 0.5 * (u(I,j,k)+u(I-1,j,k))
        V_H(k) = 0.5 * (v(i,J,k)+v(i,J-1,k))
      enddo
      if (CS%StokesMOST) then
        do k=1,GV%ke
          uE_H(k) = 0.5 * (u(I,j,k)+u(I-1,j,k)-Waves%US_x(I,j,k)-Waves%US_x(I-1,j,k))
          vE_H(k) = 0.5 * (v(i,J,k)+v(i,J-1,k)-Waves%US_y(i,J,k)-Waves%US_y(i,J-1,k))
        enddo
      endif
      ! things independent of position within the column
      Coriolis = 0.25*US%s_to_T*( (G%CoriolisBu(i,j)   + G%CoriolisBu(i-1,j-1)) + &
                                  (G%CoriolisBu(i-1,j) + G%CoriolisBu(i,j-1)) )
      surfFricVel = US%Z_to_m*US%s_to_T * uStar(i,j)

      ! Bulk Richardson number computed for each cell in a column,
      ! assuming OBLdepth = grid cell depth. After Rib(k) is
      ! known for the column, then CVMix interpolates to find
      ! the actual OBLdepth. This approach avoids need to iterate
      ! on the OBLdepth calculation. It follows that used in MOM5
      ! and POP.
      iFaceHeight(1) = 0.0 ! BBL is all relative to the surface
      pRef = 0. ; if (associated(tv%p_surf)) pRef = tv%p_surf(i,j)
      hcorr = 0.

      if (CS%StokesMOST)  call Compute_StokesDrift( i, j, h(i,j,1) , iFaceHeight(1),  &
             uS_Hi(1), vS_Hi(1), uS_H(1), vS_H(1), uSbar_H(1), vSbar_H(1), Waves)

      do k=1,GV%ke
        ! cell center and cell bottom in meters (negative values in the ocean)
        dh = dz(i,j,k) ! Nominal thickness to use for increment
        dh = dh + hcorr ! Take away the accumulated error (could temporarily make dh<0)
        hcorr = min( dh - CS%min_thickness, 0. ) ! If inflating then hcorr<0
        dh = max( dh, CS%min_thickness ) ! Limit increment dh>=min_thickness
        cellHeight(k)    = iFaceHeight(k) - 0.5 * dh
        iFaceHeight(k+1) = iFaceHeight(k) - dh

        ! find ksfc for cell where "surface layer" sits
        SLdepth_0d = CS%surf_layer_ext*max( max(-cellHeight(k),-iFaceHeight(2) ), CS%minOBLdepth )
        ksfc = k
        do ktmp = 1,k
          if (-1.0*iFaceHeight(ktmp+1) >= SLdepth_0d) then
            ksfc = ktmp
            exit
          endif
        enddo

        if (CS%StokesMOST) then
          ! if k=1, want buoyFlux(i,j,1) - buoyFlux(i,j,2), otherwise
          ! subtract average of buoyFlux(i,j,k) and buoyFlux(i,j,k+1)
          surfBuoyFlux   = buoy_scale * &
                          (buoyFlux(i,j,1) - 0.5*(buoyFlux(i,j,max(2,k))+buoyFlux(i,j,k+1)) )
          surfBuoyFlux2(k) = surfBuoyFlux
          call Compute_StokesDrift(i,j, iFaceHeight(k),iFaceHeight(k+1), &
              uS_Hi(k+1), vS_Hi(k+1), uS_H(k), vS_H(k), uSbar_H(k), vSbar_H(k), Waves)
          call Compute_StokesDrift(i,j, iFaceHeight(ksfc) , -SLdepth_0d,  &
              uS_SLD  , vS_SLD, uS_SLC , vS_SLC,  uSbar_SLD, vSbar_SLD, Waves)
          call cvmix_kpp_compute_StokesXi( iFaceHeight,CellHeight,ksfc ,SLdepth_0d,surfBuoyFlux, &
               surfFricVel,waves%omega_w2x(i,j), uE_H, vE_H, uS_Hi, vS_Hi, uSbar_H, vSbar_H, uS_SLD,&
               vS_SLD, uSbar_SLD, vSbar_SLD, StokesXI, CVMix_kpp_params_user=CS%KPP_params )
          StokesXI_1d(k) = StokesXI
          StokesVt_1d(k) = 0.0  ! StokesXI

          ! average temperature, salinity, u and v over surface layer starting at ksfc
          delH      = SLdepth_0d + iFaceHeight(ksfc)
          surfHtemp = Temp(i,j,ksfc) * delH
          surfHsalt = Salt(i,j,ksfc) * delH
          surfHu    = (uE_H(ksfc) + uSbar_SLD)  * delH
          surfHv    = (vE_H(ksfc) + vSbar_SLD)  * delH
          hTot      = delH
          do ktmp = 1,ksfc-1                            ! if ksfc >=2
            delH = h(i,j,ktmp)*GV%H_to_Z
            hTot = hTot + delH
            surfHtemp = surfHtemp + Temp(i,j,ktmp) * delH
            surfHsalt = surfHsalt + Salt(i,j,ktmp) * delH
            surfHu    = surfHu + (uE_H(ktmp) + uSbar_H(ktmp)) * delH
            surfHv    = surfHv + (vE_H(ktmp) + vSbar_H(ktmp)) * delH
          enddo
          I_hTot = 1./hTot
          surfTemp = surfHtemp * I_hTot
          surfSalt = surfHsalt * I_hTot
          surfU    = surfHu    * I_hTot
          surfV    = surfHv    * I_hTot

          Uk       = uE_H(k) + uS_H(k) - surfU
          Vk       = vE_H(k) + vS_H(k) - surfV

        else   !not StokesMOST
          StokesXI_1d(k) = 0.0
          ! average temperature, salinity, u and v over surface layer
          ! use C-grid average to get u and v on T-points.
          surfHtemp = 0.0
          surfHsalt = 0.0
          surfHu    = 0.0
          surfHv    = 0.0
          surfHuS   = 0.0
          surfHvS   = 0.0
          hTot      = 0.0
          do ktmp = 1,ksfc

            ! SLdepth_0d can be between cell interfaces
            delH = min( max(0.0, SLdepth_0d - hTot), dz(i,j,ktmp) )

            ! surface layer thickness
            hTot = hTot + delH

            ! surface averaged fields
            surfHtemp = surfHtemp + Temp(i,j,ktmp) * delH
            surfHsalt = surfHsalt + Salt(i,j,ktmp) * delH
            surfHu    = surfHu + 0.5*(u(i,j,ktmp)+u(i-1,j,ktmp)) * delH
            surfHv    = surfHv + 0.5*(v(i,j,ktmp)+v(i,j-1,ktmp)) * delH
            if (CS%Stokes_Mixing) then
              surfHus = surfHus + 0.5*(Waves%US_x(i,j,ktmp)+Waves%US_x(i-1,j,ktmp)) * delH
              surfHvs = surfHvs + 0.5*(Waves%US_y(i,j,ktmp)+Waves%US_y(i,j-1,ktmp)) * delH
            endif

          enddo
          !I_hTot = 1./hTot
          !surfTemp = surfHtemp * I_hTot
          !surfSalt = surfHsalt * I_hTot
          !surfU    = surfHu    * I_hTot
          !surfV    = surfHv    * I_hTot
          !surfUs   = surfHus   * I_hTot
          !surfVs   = surfHvs   * I_hTot

          surfTemp = surfHtemp / hTot
          surfSalt = surfHsalt / hTot
          surfU    = surfHu    / hTot
          surfV    = surfHv    / hTot
          surfUs   = surfHus   / hTot
          surfVs   = surfHvs   / hTot
          ! vertical shear between present layer and surface layer averaged surfU and surfV.
          ! C-grid average to get Uk and Vk on T-points.
          Uk         = 0.5*(u(i,j,k)+u(i-1,j,k)) - surfU
          Vk         = 0.5*(v(i,j,k)+v(i,j-1,k)) - surfV

          if (CS%Stokes_Mixing) then
            ! If momentum is mixed down the Stokes drift gradient, then
            !  the Stokes drift must be included in the bulk Richardson number
            !  calculation.
            Uk =  Uk + (0.5*(Waves%Us_x(i,j,k)+Waves%US_x(i-1,j,k)) - surfUs )
            Vk =  Vk + (0.5*(Waves%Us_y(i,j,k)+Waves%Us_y(i,j-1,k)) - surfVs )
          endif

          ! this difference accounts for penetrating SW
          surfBuoyFlux     = buoy_scale *  (buoyFlux(i,j,1) - buoyFlux(i,j,k+1))
          surfBuoyFlux2(k) = surfBuoyFlux

        endif     ! StokesMOST

        deltaU2(k) = US%L_T_to_m_s**2 * ((Uk**2) + (Vk**2))

        ! pressure, temperature, and salinity for calling the equation of state
        ! kk+1 = surface fields
        ! kk+2 = k fields
        ! kk+3 = km1 fields
        km1  = max(1, k-1)
        kk   = 3*(k-1)
        pres_1D(kk+1) = pRef
        pres_1D(kk+2) = pRef
        pres_1D(kk+3) = pRef
        Temp_1D(kk+1) = surfTemp
        Temp_1D(kk+2) = Temp(i,j,k)
        Temp_1D(kk+3) = Temp(i,j,km1)
        Salt_1D(kk+1) = surfSalt
        Salt_1D(kk+2) = Salt(i,j,k)
        Salt_1D(kk+3) = Salt(i,j,km1)

        ! pRef is pressure at interface between k and km1 [R L2 T-2 ~> Pa].
        ! iterate pRef for next pass through k-loop.
        pRef = pRef + (GV%g_Earth * GV%H_to_RZ) * h(i,j,k)

      enddo ! k-loop finishes

      if ( (CS%LT_K_ENHANCEMENT .or. CS%LT_VT2_ENHANCEMENT) .and. .not. present(lamult)) then
        MLD_guess = max( CS%MLD_guess_min, abs(CS%OBLdepthprev(i,j) ) )
        call get_Langmuir_Number(LA, G, GV, US, MLD_guess, uStar(i,j), i, j, &
                                 dz=dz(i,j,:), U_H=U_H, V_H=V_H, WAVES=WAVES)
        CS%La_SL(i,j) = LA
      endif


      ! compute in-situ density
      call calculate_density(Temp_1D, Salt_1D, pres_1D, rho_1D, tv%eqn_of_state)

      ! N2 (can be negative) and N (non-negative) on interfaces.
      ! deltaRho is non-local rho difference used for bulk Richardson number.
      ! CS%N is local N (with floor) used for unresolved shear calculation.
      do k = 1, GV%ke
        km1 = max(1, k-1)
        kk = 3*(k-1)
        deltaRho(k) = rho_1D(kk+2) - rho_1D(kk+1)
        if (GV%Boussinesq .or. GV%semi_Boussinesq) then
          deltaBuoy(k) = GoRho*(rho_1D(kk+2) - rho_1D(kk+1))
        else
          deltaBuoy(k) = (US%Z_to_m*US%s_to_T**2) * GV%g_Earth_Z_T2 * &
              ( (rho_1D(kk+2) - rho_1D(kk+1)) / (0.5 * (rho_1D(kk+2) + rho_1D(kk+1))) )
        endif
        N2_1d(k)    = (GoRho_Z_L2 * (rho_1D(kk+2) - rho_1D(kk+3)) ) / &
                      ((0.5*(h(i,j,km1) + h(i,j,k))+GV%H_subroundoff))
        CS%N(i,j,k)     = sqrt( max( N2_1d(k), 0.) )
      enddo
      N2_1d(GV%ke+1 ) = 0.0
      CS%N(i,j,GV%ke+1 )  = 0.0

      ! Convert columns to MKS units for passing to CVMix
      do k = 1, GV%ke
        OBL_depth(k) = -US%Z_to_m * cellHeight(k)
        z_cell(k) = US%Z_to_m*cellHeight(k)
      enddo
      do K = 1, GV%ke+1
        N_col(K) = US%s_to_T*CS%N(i,j,K)
        z_inter(K) = US%Z_to_m*iFaceHeight(K)
      enddo

      ! CVMix_kpp_compute_turbulent_scales_1d_OBL computes w_s velocity scale at cell centers for
      ! CVmix_kpp_compute_bulk_Richardson call to CVmix_kpp_compute_unresolved_shear
      ! at sigma=Vt_layer (CS%surf_layer_ext or 1.0) for this calculation.
      ! StokesVt_1d controls Stokes enhancement  (= 0 for none)
      Vt_layer = 1.0      ! CS%surf_layer_ext
      call CVMix_kpp_compute_turbulent_scales( &                             ! 1d_OBL
              Vt_layer,          & ! (in)  Boundary layer extent contributing to unresolved shear
              OBL_depth,         & ! (in)  OBL depth [m]
              surfBuoyFlux2,     & ! (in)  Buoyancy flux at surface [m2 s-3]
              surfFricVel,       & ! (in)  Turbulent friction velocity at surface [m s-1]
              xi=StokesVt_1d,    & ! (in)  Stokes similarity parameter-->1/CHI(xi) enhance of Vt
              w_s=Ws_1d,         & ! (out) Turbulent velocity scale profile [m s-1]
              CVMix_kpp_params_user=CS%KPP_params )

      ! Determine the enhancement factor for unresolved shear
      IF (CS%LT_VT2_ENHANCEMENT) then
        IF (CS%LT_VT2_METHOD==LT_VT2_MODE_CONSTANT) then
          LangEnhVT2 = CS%KPP_VT2_ENH_FAC
        elseif (CS%LT_VT2_METHOD==LT_VT2_MODE_VR12) then
          !Introduced minimum value for La_SL, so maximum value for enhvt2 is removed.
          if (present(lamult)) then
            LangEnhVT2 = lamult(i,j)
          else
            LangEnhVT2 = sqrt(1.+(1.5*CS%La_SL(i,j))**(-2) + &
                      (5.4*CS%La_SL(i,j))**(-4))
          endif
        else
          ! for other methods (e.g., LT_VT2_MODE_RW16, LT_VT2_MODE_LF17), the enhancement factor is
          ! computed internally within CVMix using LaSL, bfsfc, and ustar to be passed to CVMix.
          LangEnhVT2 = 1.0
        endif
      else
        LangEnhVT2 = 1.0
      endif

      surfBuoyFlux = buoy_scale * buoyFlux(i,j,1)

      ! Calculate Bulk Richardson number from eq (21) of LMD94
      BulkRi_1d = CVmix_kpp_compute_bulk_Richardson( &
                  zt_cntr=z_cell,                    & ! Depth of cell center [m]
                  delta_buoy_cntr=deltaBuoy,         & ! Bulk buoyancy difference, Br-B(z) [m s-2]
                  delta_Vsqr_cntr=deltaU2,           & ! Square of resolved velocity difference [m2 s-2]
                  ws_cntr=Ws_1d,                     & ! Turbulent velocity scale profile [m s-1]
                  N_iface=N_col,                     & ! Buoyancy frequency [s-1]
                  EFactor=LangEnhVT2,                & ! Langmuir enhancement factor [nondim]
                  LaSL=CS%La_SL(i,j),                & ! surface layer averaged Langmuir number [nondim]
                  bfsfc=surfBuoyFlux2,               & ! surface buoyancy flux [m2 s-3]
                  uStar=surfFricVel,                 & ! surface friction velocity [m s-1]
                  CVMix_kpp_params_user=CS%KPP_params ) ! KPP parameters

!      ! A hack to avoid KPP reaching the bottom. It was needed during development
!      ! because KPP was unable to handle vanishingly small layers near the bottom.
!      if (CS%deepOBLoffset>0.) then
!        zBottomMinusOffset = iFaceHeight(GV%ke+1) + min(CS%deepOBLoffset, -0.1*iFaceHeight(GV%ke+1))
!        CS%OBLdepth(i,j) = min( CS%OBLdepth(i,j), -zBottomMinusOffset )
!      endif
      zBottomMinusOffset = iFaceHeight(GV%ke+1) + min(CS%deepOBLoffset,-0.1*iFaceHeight(GV%ke+1))

      call CVMix_kpp_compute_OBL_depth( &
            BulkRi_1d,              & ! (in) Bulk Richardson number
            z_inter,                & ! (in) Height of interfaces [m]
            KPP_OBL_depth,          & ! (out) OBL depth [m]
            CS%kOBL(i,j),           & ! (out) level (+fraction) of OBL extent
            zt_cntr=z_cell,         & ! (in) Height of cell centers [m]
            surf_fric=surfFricVel,  & ! (in) Turbulent friction velocity at surface [m s-1]
            surf_buoy=surfBuoyFlux2, & ! (in) Buoyancy flux at surface [m2 s-3]
            Coriolis=Coriolis,      & ! (in) Coriolis parameter [s-1]
            Xi = StokesXI_1d,       & ! (in) Stokes similarity parameter Lmob limit (1-Xi)
            zBottom = zBottomMinusOffset,     & ! (in) Numerical limit on OBLdepth
            CVMix_kpp_params_user=CS%KPP_params ) ! KPP parameters
      CS%OBLdepth(i,j) = US%m_to_Z * KPP_OBL_depth

    if (CS%StokesMOST) then
      kbl = int(CS%kOBL(i,j))
      SLdepth_0d = CS%surf_layer_ext*CS%OBLdepth(i,j)
      surfBuoyFlux = surfBuoyFlux2(kbl)
        ! find ksfc for cell where "surface layer" sits
      ksfc = kbl
      do ktmp = 1, kbl
        if (-1.0*iFaceHeight(ktmp+1) >= SLdepth_0d) then
          ksfc = ktmp
          exit
        endif
      enddo

      call Compute_StokesDrift(i,j, iFaceHeight(ksfc) , -SLdepth_0d,  &
              uS_SLD  , vS_SLD, uS_SLC , vS_SLC,  uSbar_SLD, vSbar_SLD, Waves)
      call cvmix_kpp_compute_StokesXi( iFaceHeight,CellHeight,ksfc ,SLdepth_0d,  &
               surfBuoyFlux, surfFricVel,waves%omega_w2x(i,j), uE_H, vE_H, uS_Hi, &
               vS_Hi, uSbar_H, vSbar_H, uS_SLD, vS_SLD, uSbar_SLD, vSbar_SLD,     &
               StokesXI, CVMix_kpp_params_user=CS%KPP_params )
      CS%StokesParXI(i,j) = StokesXI
      CS%Lam2(i,j)        = sqrt(US_Hi(1)**2+VS_Hi(1)**2) / MAX(surfFricVel,0.0002)

    else                         !.not Stokes_MOST
      CS%StokesParXI(i,j) = 10.0
      CS%Lam2(i,j)        = sqrt(US_Hi(1)**2+VS_Hi(1)**2) / MAX(surfFricVel,0.0002)

      ! A hack to avoid KPP reaching the bottom. It was needed during development
      ! because KPP was unable to handle vanishingly small layers near the bottom.
      if (CS%deepOBLoffset>0.) then
        zBottomMinusOffset = iFaceHeight(GV%ke+1) + min(CS%deepOBLoffset, -0.1*iFaceHeight(GV%ke+1))
        CS%OBLdepth(i,j) = min( CS%OBLdepth(i,j), -zBottomMinusOffset )
      endif

      ! apply some constraints on OBLdepth
      if (CS%fixedOBLdepth) CS%OBLdepth(i,j) = CS%fixedOBLdepth_value
      CS%OBLdepth(i,j) = max( CS%OBLdepth(i,j), -iFaceHeight(2) )       ! no shallower than top layer
      CS%OBLdepth(i,j) = min( CS%OBLdepth(i,j), -iFaceHeight(GV%ke+1) ) ! no deeper than bottom
      CS%kOBL(i,j)     = CVMix_kpp_compute_kOBL_depth( iFaceHeight, cellHeight, CS%OBLdepth(i,j) )

    endif                              !Stokes_MOST

      ! compute unresolved squared velocity for diagnostics
      if (CS%id_Vt2 > 0) then
        Vt2_1d(:) = CVmix_kpp_compute_unresolved_shear( &
                    z_cell,             & ! Depth of cell center [m]
                    ws_cntr=Ws_1d,      & ! Turbulent velocity scale profile, at centers [m s-1]
                    N_iface=N_col,      & ! Buoyancy frequency at interface [s-1]
                    EFactor=LangEnhVT2, & ! Langmuir enhancement factor [nondim]
                    LaSL=CS%La_SL(i,j), & ! surface layer averaged Langmuir number [nondim]
                    bfsfc=surfBuoyFlux2, & ! surface buoyancy flux [m2 s-3]
                    uStar=surfFricVel,  & ! surface friction velocity [m s-1]
                    CVmix_kpp_params_user=CS%KPP_params ) ! KPP parameters
        CS%Vt2(i,j,:) = US%m_to_Z**2*US%T_to_s**2 * Vt2_1d(:)
      endif

      ! recompute wscale for diagnostics, now that we in fact know boundary layer depth
      !BGR consider if LTEnhancement is wanted for diagnostics
      if (CS%id_Ws > 0) then
        call CVMix_kpp_compute_turbulent_scales( &
            -cellHeight(:)/CS%OBLdepth(i,j),       & ! (in)  Normalized boundary layer coordinate [nondim]
            US%Z_to_m*CS%OBLdepth(i,j),            & ! (in)  OBL depth [m]
            surfBuoyFlux,                          & ! (in)  Buoyancy flux at surface [m2 s-3]
            surfFricVel,                           & ! (in)  Turbulent friction velocity at surface [m s-1]
            xi=StokesXI,                           & ! (in) Stokes similarity parameter-->1/CHI(xi) enhance
            w_s=Ws_1d,                             & ! (out) Turbulent velocity scale profile [m s-1]
            CVMix_kpp_params_user=CS%KPP_params)     !       KPP parameters
        CS%Ws(i,j,:) = US%m_to_Z*US%T_to_s*Ws_1d(:)
      endif

      ! Diagnostics
      if (CS%id_N2     > 0)   CS%N2(i,j,:)     = N2_1d(:)
      if (CS%id_BulkDrho > 0) CS%dRho(i,j,:)   = deltaRho(:)
      if (CS%id_BulkRi > 0)   CS%BulkRi(i,j,:) = BulkRi_1d(:)
      if (CS%id_BulkUz2 > 0)  CS%Uz2(i,j,:)    = US%m_s_to_L_T**2 * deltaU2(:)
      if (CS%id_Tsurf  > 0)   CS%Tsurf(i,j)    = surfTemp
      if (CS%id_Ssurf  > 0)   CS%Ssurf(i,j)    = surfSalt
      if (CS%id_Usurf  > 0)   CS%Usurf(i,j)    = surfU
      if (CS%id_Vsurf  > 0)   CS%Vsurf(i,j)    = surfV

    endif ; enddo
  enddo

  call cpu_clock_end(id_clock_KPP_compute_BLD)

  ! send diagnostics to post_data
  if (CS%id_BulkRi   > 0) call post_data(CS%id_BulkRi,   CS%BulkRi,          CS%diag)
  if (CS%id_N        > 0) call post_data(CS%id_N,        CS%N,               CS%diag)
  if (CS%id_N2       > 0) call post_data(CS%id_N2,       CS%N2,              CS%diag)
  if (CS%id_Tsurf    > 0) call post_data(CS%id_Tsurf,    CS%Tsurf,           CS%diag)
  if (CS%id_Ssurf    > 0) call post_data(CS%id_Ssurf,    CS%Ssurf,           CS%diag)
  if (CS%id_Usurf    > 0) call post_data(CS%id_Usurf,    CS%Usurf,           CS%diag)
  if (CS%id_Vsurf    > 0) call post_data(CS%id_Vsurf,    CS%Vsurf,           CS%diag)
  if (CS%id_BulkDrho > 0) call post_data(CS%id_BulkDrho, CS%dRho,            CS%diag)
  if (CS%id_BulkUz2  > 0) call post_data(CS%id_BulkUz2,  CS%Uz2,             CS%diag)
  if (CS%id_EnhK     > 0) call post_data(CS%id_EnhK,     CS%EnhK,            CS%diag)
  if (CS%id_EnhVt2   > 0) call post_data(CS%id_EnhVt2,   CS%EnhVt2,          CS%diag)
  if (CS%id_La_SL    > 0) call post_data(CS%id_La_SL,    CS%La_SL,           CS%diag)
  if (CS%id_Vt2      > 0) call post_data(CS%id_Vt2,      CS%Vt2,             CS%diag)

  if (CS%StokesMOST) then
    if (CS%id_StokesXI > 0) call post_data(CS%id_StokesXI, CS%StokesParXI,     CS%diag)
    if (CS%id_Lam2     > 0) call post_data(CS%id_Lam2    , CS%Lam2    ,        CS%diag)
  endif

  ! BLD smoothing:
  if (CS%n_smooth > 0) call KPP_smooth_BLD(CS, G, GV, US, dz)

end subroutine KPP_compute_BLD


!> Apply a 1-1-4-1-1 Laplacian filter one time on BLD to reduce any horizontal two-grid-point noise
subroutine KPP_smooth_BLD(CS, G, GV, US, dz)
  ! Arguments
  type(KPP_CS),                           pointer       :: CS   !< Control structure
  type(ocean_grid_type),                  intent(inout) :: G    !< Ocean grid
  type(verticalGrid_type),                intent(in)    :: GV   !< Ocean vertical grid
  type(unit_scale_type),                  intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: dz   !< Layer thicknesses [Z ~> m]

  ! local variables
  real, dimension(SZI_(G),SZJ_(G)) :: OBLdepth_prev     ! OBLdepth before s.th smoothing iteration [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)) :: total_depth       ! The total depth of the water column, adjusted
                                                        ! for the minimum layer thickness [Z ~> m]
  real, dimension( GV%ke )         :: cellHeight        ! Cell center heights referenced to surface [Z ~> m]
                                                        ! (negative in the ocean)
  real, dimension( GV%ke+1 )       :: iFaceHeight       ! Interface heights referenced to surface [Z ~> m]
                                                        ! (negative in the ocean)
  real :: wc, ww, we, wn, ws ! averaging weights for smoothing [nondim]
  real :: dh                 ! The local thickness used for calculating interface positions [Z ~> m]
  real :: h_cor(SZI_(G))     ! A cumulative correction arising from inflation of vanished layers [Z ~> m]
  real :: hcorr              ! A cumulative correction arising from inflation of vanished layers [Z ~> m]
  integer :: i, j, k, s, halo

  call cpu_clock_begin(id_clock_KPP_smoothing)

  ! Find the total water column thickness first, as it is reused for each smoothing pass.
  total_depth(:,:) = 0.0

  !$OMP parallel do default(shared) private(dh, h_cor)
  do j = G%jsc, G%jec
    h_cor(:) = 0.
    do k=1,GV%ke
      do i=G%isc,G%iec ; if (G%mask2dT(i,j) > 0.0) then
        ! This code replicates the interface height calculations below.  It could be simpler, as shown below.
        dh = dz(i,j,k)   ! Nominal thickness to use for increment
        dh = dh + h_cor(i) ! Take away the accumulated error (could temporarily make dh<0)
        h_cor(i) = min( dh - CS%min_thickness, 0. ) ! If inflating then hcorr<0
        dh = max( dh, CS%min_thickness ) ! Limit increment dh>=min_thickness
        total_depth(i,j) = total_depth(i,j) + dh
      endif ; enddo
    enddo
  enddo
  ! A much simpler (but answer changing) version of the total_depth calculation would be
  ! do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
  !   total_depth(i,j) = total_depth(i,j) + dz(i,j,k)
  ! enddo ; enddo ; enddo

  ! Update halos once, then march inward for each iteration
  if (CS%n_smooth > 1) call pass_var(total_depth, G%Domain, halo=CS%n_smooth, complete=.false.)
  call pass_var(CS%OBLdepth, G%Domain, halo=CS%n_smooth)

  if (CS%id_OBLdepth_original > 0) CS%OBLdepth_original(:,:) = CS%OBLdepth(:,:)

  do s=1,CS%n_smooth

    OBLdepth_prev(:,:) = CS%OBLdepth(:,:)
    halo = CS%n_smooth - s

    ! apply smoothing on OBL depth
    !$OMP parallel do default(none) shared(G, GV, CS, OBLdepth_prev, total_depth, halo) &
    !$OMP                           private(wc, ww, we, wn, ws)
    do j = G%jsc-halo, G%jec+halo
      do i = G%isc-halo, G%iec+halo ; if (G%mask2dT(i,j) > 0.0) then
        ! compute weights
        ww = 0.125 * G%mask2dT(i-1,j)
        we = 0.125 * G%mask2dT(i+1,j)
        ws = 0.125 * G%mask2dT(i,j-1)
        wn = 0.125 * G%mask2dT(i,j+1)
        wc = 1.0 - (ww+we+wn+ws)

        if (CS%answer_date < 20240501) then
          CS%OBLdepth(i,j) =  wc * OBLdepth_prev(i,j)   &
                            + ww * OBLdepth_prev(i-1,j) &
                            + we * OBLdepth_prev(i+1,j) &
                            + ws * OBLdepth_prev(i,j-1) &
                            + wn * OBLdepth_prev(i,j+1)
        else
          CS%OBLdepth(i,j) =  wc * OBLdepth_prev(i,j) &
                            + ((ww * OBLdepth_prev(i-1,j) + we * OBLdepth_prev(i+1,j)) &
                             + (ws * OBLdepth_prev(i,j-1) + wn * OBLdepth_prev(i,j+1)))
        endif

        ! Apply OBLdepth smoothing at a cell only if the OBLdepth gets deeper via smoothing.
        if (CS%deepen_only) CS%OBLdepth(i,j) = max(CS%OBLdepth(i,j), OBLdepth_prev(i,j))

        ! prevent OBL depths deeper than the bathymetric depth
        CS%OBLdepth(i,j) = min( CS%OBLdepth(i,j), total_depth(i,j) ) ! no deeper than bottom
      endif ; enddo
    enddo

  enddo ! s-loop

  ! Determine the fractional index of the bottom of the boundary layer.
  !$OMP parallel do default(none) shared(G, GV, CS, dz) &
  !$OMP                           private(dh, hcorr, cellHeight, iFaceHeight)
  do j=G%jsc,G%jec ; do i=G%isc,G%iec ; if (G%mask2dT(i,j) > 0.0) then

    iFaceHeight(1) = 0.0 ! BBL is all relative to the surface
    hcorr = 0.
    do k=1,GV%ke
      ! cell center and cell bottom in meters (negative values in the ocean)
      dh = dz(i,j,k)   ! Nominal thickness to use for increment
      dh = dh + hcorr  ! Take away the accumulated error (could temporarily make dh<0)
      hcorr = min( dh - CS%min_thickness, 0. ) ! If inflating then hcorr<0
      dh = max( dh, CS%min_thickness ) ! Limit increment dh>=min_thickness
      cellHeight(k)    = iFaceHeight(k) - 0.5 * dh
      iFaceHeight(k+1) = iFaceHeight(k) - dh
    enddo

    CS%kOBL(i,j) = CVMix_kpp_compute_kOBL_depth( iFaceHeight, cellHeight, CS%OBLdepth(i,j) )
  endif ; enddo ; enddo

  call cpu_clock_end(id_clock_KPP_smoothing)

end subroutine KPP_smooth_BLD


!> Copies KPP surface boundary layer depth into BLD, in units of [Z ~> m] unless other units are specified.
subroutine KPP_get_BLD(CS, BLD, G, US, m_to_BLD_units)
  type(KPP_CS),                     pointer     :: CS  !< Control structure for
                                                       !! this module
  type(ocean_grid_type),            intent(in)  :: G   !< Grid structure
  type(unit_scale_type),            intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: BLD !< Boundary layer depth [Z ~> m] or other units
  real,                   optional, intent(in)  :: m_to_BLD_units !< A conversion factor from meters
                                                       !! to the desired units for BLD [various]
  ! Local variables
  real :: scale  ! A dimensional rescaling factor in [nondim] or other units.
  integer :: i,j

  scale = 1.0 ; if (present(m_to_BLD_units)) scale = US%Z_to_m*m_to_BLD_units

  !$OMP parallel do default(none) shared(BLD, CS, G, scale)
  do j = G%jsc, G%jec ; do i = G%isc, G%iec
    BLD(i,j) = scale * CS%OBLdepth(i,j)
  enddo ; enddo

end subroutine KPP_get_BLD

!> Apply KPP non-local transport of surface fluxes for a given tracer
subroutine KPP_NonLocalTransport(CS, G, GV, h, nonLocalTrans, surfFlux, &
                                 dt, diag, tr_ptr, scalar, flux_scale)
  type(KPP_CS),                               intent(in)    :: CS            !< Control structure
  type(ocean_grid_type),                      intent(in)    :: G             !< Ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV            !< Ocean vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h             !< Layer/level thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)   :: nonLocalTrans !< Non-local transport [nondim]
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: surfFlux      !< Surface flux of scalar
                                                                        !! [conc H T-1 ~> conc m s-1 or conc kg m-2 s-1]
  real,                                       intent(in)    :: dt            !< Time-step [T ~> s]
  type(diag_ctrl), target,                    intent(in)    :: diag          !< Diagnostics
  type(tracer_type), pointer,                 intent(in)    :: tr_ptr        !< tracer_type has diagnostic ids on it
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: scalar        !< Scalar (scalar units [conc])
  real, optional,                             intent(in)    :: flux_scale    !< Scale factor to get surfFlux
                                                                             !! into proper units [various]

  integer :: i, j, k
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: dtracer ! Rate of tracer change [conc T-1 ~> conc s-1]
  real, dimension(SZI_(G),SZJ_(G)) :: surfFlux_loc ! An optionally rescaled surface flux of the scalar
                                                   ! in [conc H T-1 ~> conc m s-1 or conc kg m-2 s-1] or other units

  ! term used to scale
  if (present(flux_scale)) then
    do j = G%jsc, G%jec ; do i = G%isc, G%iec
      surfFlux_loc(i,j) = surfFlux(i,j) * flux_scale
    enddo ; enddo
  else
    surfFlux_loc(:,:) = surfFlux(:,:)
  endif

  ! Post surface flux diagnostic
  if (tr_ptr%id_net_surfflux > 0) call post_data(tr_ptr%id_net_surfflux, surfFlux_loc(:,:), diag)

  ! Only continue if we are applying the nonlocal tendency
  ! or the nonlocal tendency diagnostic has been requested
  if ((tr_ptr%id_NLT_tendency > 0) .or. (CS%applyNonLocalTrans)) then

    !$OMP parallel do default(none) shared(dtracer, nonLocalTrans, h, G, GV, surfFlux_loc)
    do k = 1, GV%ke ; do j = G%jsc, G%jec ; do i = G%isc, G%iec
      dtracer(i,j,k) = ( nonLocalTrans(i,j,k) - nonLocalTrans(i,j,k+1) ) / &
                       ( h(i,j,k) + GV%H_subroundoff ) * surfFlux_loc(i,j)
    enddo ; enddo ; enddo

    !  Update tracer due to non-local redistribution of surface flux
    if (CS%applyNonLocalTrans) then
      !$OMP parallel do default(none) shared(G, GV, dt, scalar, dtracer)
      do k = 1, GV%ke ; do j = G%jsc, G%jec ; do i = G%isc, G%iec
        scalar(i,j,k) = scalar(i,j,k) + dt * dtracer(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (tr_ptr%id_NLT_tendency > 0) call post_data(tr_ptr%id_NLT_tendency, dtracer,  diag)

  endif


  if (tr_ptr%id_NLT_budget > 0) then
    !$OMP parallel do default(none) shared(G, GV, dtracer, nonLocalTrans, surfFlux_loc)
    do k = 1, GV%ke ; do j = G%jsc, G%jec ; do i = G%isc, G%iec
      ! Here dtracer has units of [Q R Z T-1 ~> W m-2].
      dtracer(i,j,k) = (nonLocalTrans(i,j,k) - nonLocalTrans(i,j,k+1)) * surfFlux_loc(i,j)
    enddo ; enddo ; enddo
    call post_data(tr_ptr%id_NLT_budget, dtracer(:,:,:), diag)
  endif

end subroutine KPP_NonLocalTransport


!> Apply KPP non-local transport of surface fluxes for temperature.
subroutine KPP_NonLocalTransport_temp(CS, G, GV, h, nonLocalTrans, surfFlux, dt, tr_ptr, scalar, C_p)
  type(KPP_CS),                               intent(in)    :: CS     !< Control structure
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h      !< Layer/level thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)   :: nonLocalTrans !< Non-local transport [nondim]
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: surfFlux  !< Surface flux of temperature
                                                                      !! [C H T-1 ~> degC m s-1 or degC kg m-2 s-1]
  real,                                       intent(in)    :: dt     !< Time-step [T ~> s]
  type(tracer_type), pointer,                 intent(in)    :: tr_ptr !< tracer_type has diagnostic ids on it
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: scalar !< temperature [C ~> degC]
  real,                                       intent(in)    :: C_p    !< Seawater specific heat capacity
                                                                      !! [Q C-1 ~> J kg-1 degC-1]

  call KPP_NonLocalTransport(CS, G, GV, h, nonLocalTrans, surfFlux, dt, CS%diag, &
                             tr_ptr, scalar)

end subroutine KPP_NonLocalTransport_temp


!> Apply KPP non-local transport of surface fluxes for salinity.
subroutine KPP_NonLocalTransport_saln(CS, G, GV, h, nonLocalTrans, surfFlux, dt, tr_ptr, scalar)
  type(KPP_CS),                               intent(in)    :: CS            !< Control structure
  type(ocean_grid_type),                      intent(in)    :: G             !< Ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV            !< Ocean vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h             !< Layer/level thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)   :: nonLocalTrans !< Non-local transport [nondim]
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: surfFlux      !< Surface flux of salt
                                                                             !! [S H T-1 ~> ppt m s-1 or ppt kg m-2 s-1]
  real,                                       intent(in)    :: dt            !< Time-step [T ~> s]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: scalar        !< Salinity [S ~> ppt]
  type(tracer_type), pointer,                 intent(in)    :: tr_ptr        !< tracer_type has diagnostic ids on it

  call KPP_NonLocalTransport(CS, G, GV, h, nonLocalTrans, surfFlux, dt, CS%diag, &
                             tr_ptr, scalar)

end subroutine KPP_NonLocalTransport_saln

!> Compute Stokes Drift components at zbot < ztop <= 0  and at k=0.5*(ztop+zbot) and
!! average components from ztop to zbot <= 0
subroutine Compute_StokesDrift(i ,j, ztop, zbot, uS_i, vS_i, uS_k, vS_k, uSbar, vSbar, waves)

  type(wave_parameters_CS), pointer  :: waves  !< Wave CS for Langmuir turbulence
  real,                intent(in)    :: ztop   !< cell top
  real,                intent(in)    :: zbot   !< cell bottom
  real,                intent(inout) :: uS_i   !< Stokes u velocity at zbot interface
  real,                intent(inout) :: vS_i   !< Stokes v velocity at zbot interface
  real,                intent(inout) :: uS_k   !< Stokes u velocity at zk center
  real,                intent(inout) :: vS_k   !< Stokes v at zk =0.5(ztop+zbot)
  real,                intent(inout) :: uSbar  !< mean Stokes u (ztop to zbot)
  real,                intent(inout) :: vSbar  !< mean Stokes v (ztop to zbot)
  integer,             intent(in)    :: i      !< Meridional index of H-point
  integer,             intent(in)    :: j      !< Zonal index of H-point

  ! local variables
  integer                            ::   b     !< wavenumber band index
  real                               :: fexp    !< an exponential function
  real                               :: WaveNum !< Wavenumber

  uS_i  = 0.0
  vS_i  = 0.0
  uS_k  = 0.0
  vS_k  = 0.0
  uSbar = 0.0
  vSbar = 0.0
  do b  = 1, waves%NumBands
    WaveNum =  waves%WaveNum_Cen(b)
    fexp  =  exp(2. * WaveNum * zbot)
    uS_i  =  uS_i + waves%Ustk_Hb(i,j,b) * fexp
    vS_i  =  vS_i + waves%Vstk_Hb(i,j,b) * fexp
    fexp  =  exp( WaveNum * (ztop + zbot) )
    uS_k  =  uS_k+ waves%Ustk_Hb(i,j,b) * fexp
    vS_k  =  vS_k+ waves%Vstk_Hb(i,j,b) * fexp
    fexp  =  exp(2. * WaveNum * ztop) - exp(2. * WaveNum * zbot)
    uSbar =  uSbar + 0.5 * waves%Ustk_Hb(i,j,b) * fexp / WaveNum
    vSbar =  vSbar + 0.5 * waves%Vstk_Hb(i,j,b) * fexp / WaveNum
  enddo
  uSbar = uSbar / (ztop-zbot)
  vSbar = vSbar / (ztop-zbot)

end subroutine Compute_StokesDrift

!> Clear pointers, deallocate memory
subroutine KPP_end(CS)
  type(KPP_CS), pointer :: CS !< Control structure

  if (.not.associated(CS)) return

  deallocate(CS)

end subroutine KPP_end

end module MOM_CVMix_KPP
