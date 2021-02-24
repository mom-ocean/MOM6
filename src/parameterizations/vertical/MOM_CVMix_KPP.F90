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
use MOM_unit_scaling,   only : unit_scale_type
use MOM_variables,      only : thermo_var_ptrs
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS, Get_Langmuir_Number
use MOM_domains,        only : pass_var
use MOM_cpu_clock,      only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,      only : CLOCK_MODULE, CLOCK_ROUTINE

use CVMix_kpp, only : CVMix_init_kpp, CVMix_put_kpp, CVMix_get_kpp_real
use CVMix_kpp, only : CVMix_coeffs_kpp
use CVMix_kpp, only : CVMix_kpp_compute_OBL_depth
use CVMix_kpp, only : CVMix_kpp_compute_turbulent_scales
use CVMix_kpp, only : CVMix_kpp_compute_bulk_Richardson
use CVMix_kpp, only : CVMix_kpp_compute_unresolved_shear
use CVMix_kpp, only : CVMix_kpp_params_type
use CVMix_kpp, only : CVMix_kpp_compute_kOBL_depth

implicit none ; private

#include "MOM_memory.h"

public :: KPP_init
public :: KPP_compute_BLD
public :: KPP_calculate
public :: KPP_end
public :: KPP_NonLocalTransport_temp
public :: KPP_NonLocalTransport_saln
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
  real    :: Ri_crit                   !< Critical bulk Richardson number (defines OBL depth)
  real    :: vonKarman                 !< von Karman constant (dimensionless)
  real    :: cs                        !< Parameter for computing velocity scale function (dimensionless)
  real    :: cs2                       !< Parameter for multiplying by non-local term
                                       !   This is active for NLT_SHAPE_CUBIC_LMD only
  logical :: enhance_diffusion         !< If True, add enhanced diffusivity at base of boundary layer.
  character(len=10) :: interpType      !< Type of interpolation to compute bulk Richardson number
  character(len=10) :: interpType2     !< Type of interpolation to compute diff and visc at OBL_depth
  logical :: computeEkman              !< If True, compute Ekman depth limit for OBLdepth
  logical :: computeMoninObukhov       !< If True, compute Monin-Obukhov limit for OBLdepth
  logical :: passiveMode               !< If True, makes KPP passive meaning it does NOT alter the diffusivity
  real    :: deepOBLoffset             !< If non-zero, is a distance from the bottom that the OBL can not
                                       !! penetrate through [m]
  real    :: minOBLdepth               !< If non-zero, is a minimum depth for the OBL [m]
  real    :: surf_layer_ext            !< Fraction of OBL depth considered in the surface layer [nondim]
  real    :: minVtsqr                  !< Min for the squared unresolved velocity used in Rib CVMix calculation [m2 s-2]
  logical :: fixedOBLdepth             !< If True, will fix the OBL depth at fixedOBLdepth_value
  real    :: fixedOBLdepth_value       !< value for the fixed OBL depth when fixedOBLdepth==True.
  logical :: debug                     !< If True, calculate checksums and write debugging information
  character(len=30) :: MatchTechnique  !< Method used in CVMix for setting diffusivity and NLT profile functions
  integer :: NLT_shape                 !< MOM6 over-ride of CVMix NLT shape function
  logical :: applyNonLocalTrans        !< If True, apply non-local transport to heat and scalars
  integer :: n_smooth                  !< Number of times smoothing operator is applied on OBLdepth.
  logical :: deepen_only               !< If true, apply OBLdepth smoothing at a cell only if the OBLdepth gets deeper.
  logical :: KPPzeroDiffusivity        !< If True, will set diffusivity and viscosity from KPP to zero
                                       !! for testing purposes.
  logical :: KPPisAdditive             !< If True, will add KPP diffusivity to initial diffusivity.
                                       !! If False, will replace initial diffusivity wherever KPP diffusivity
                                       !! is non-zero.
  real    :: min_thickness             !< A minimum thickness used to avoid division by small numbers
                                       !! in the vicinity of vanished layers.
  ! smg: obsolete below
  logical :: correctSurfLayerAvg       !< If true, applies a correction to the averaging of surface layer properties
  real    :: surfLayerDepth            !< A guess at the depth of the surface layer (which should 0.1 of OBLdepth) [m]
  ! smg: obsolete above
  integer :: SW_METHOD                 !< Sets method for using shortwave radiation in surface buoyancy flux
  logical :: LT_K_Enhancement          !< Flags if enhancing mixing coefficients due to LT
  integer :: LT_K_Shape                !< Integer for constant or shape function enhancement
  integer :: LT_K_Method               !< Integer for mixing coefficients LT method
  real    :: KPP_K_ENH_FAC             !< Factor to multiply by K if Method is CONSTANT
  logical :: LT_Vt2_Enhancement        !< Flags if enhancing Vt2 due to LT
  integer :: LT_VT2_METHOD             !< Integer for Vt2 LT method
  real    :: KPP_VT2_ENH_FAC           !< Factor to multiply by VT2 if Method is CONSTANT
  logical :: STOKES_MIXING             !< Flag if model is mixing down Stokes gradient
                                       !! This is relavent for which current to use in RiB

  !> CVMix parameters
  type(CVMix_kpp_params_type), pointer :: KPP_params => NULL()

  type(diag_ctrl), pointer :: diag => NULL() !< Pointer to diagnostics control structure
  !>@{ Diagnostic handles
  integer :: id_OBLdepth = -1, id_BulkRi   = -1
  integer :: id_N        = -1, id_N2       = -1
  integer :: id_Ws       = -1, id_Vt2      = -1
  integer :: id_BulkUz2  = -1, id_BulkDrho = -1
  integer :: id_uStar    = -1, id_buoyFlux = -1
  integer :: id_QminusSW = -1, id_netS     = -1
  integer :: id_sigma    = -1, id_Kv_KPP   = -1
  integer :: id_Kt_KPP   = -1, id_Ks_KPP   = -1
  integer :: id_Tsurf    = -1, id_Ssurf    = -1
  integer :: id_Usurf    = -1, id_Vsurf    = -1
  integer :: id_Kd_in    = -1
  integer :: id_NLTt     = -1
  integer :: id_NLTs     = -1
  integer :: id_NLT_dSdt = -1
  integer :: id_NLT_dTdt = -1
  integer :: id_NLT_temp_budget = -1
  integer :: id_NLT_saln_budget = -1
  integer :: id_EnhK     = -1, id_EnhVt2   = -1
  integer :: id_EnhW     = -1
  integer :: id_La_SL    = -1
  integer :: id_OBLdepth_original = -1
  !>@}

  ! Diagnostics arrays
  real, allocatable, dimension(:,:)   :: OBLdepth  !< Depth (positive) of OBL [m]
  real, allocatable, dimension(:,:)   :: OBLdepth_original  !< Depth (positive) of OBL [m] without smoothing
  real, allocatable, dimension(:,:)   :: kOBL      !< Level (+fraction) of OBL extent
  real, allocatable, dimension(:,:)   :: OBLdepthprev !< previous Depth (positive) of OBL [m]
  real, allocatable, dimension(:,:)   :: La_SL     !< Langmuir number used in KPP
  real, allocatable, dimension(:,:,:) :: dRho      !< Bulk difference in density [R ~> kg m-3]
  real, allocatable, dimension(:,:,:) :: Uz2       !< Square of bulk difference in resolved velocity [m2 s-2]
  real, allocatable, dimension(:,:,:) :: BulkRi    !< Bulk Richardson number for each layer (dimensionless)
  real, allocatable, dimension(:,:,:) :: sigma     !< Sigma coordinate (dimensionless)
  real, allocatable, dimension(:,:,:) :: Ws        !< Turbulent velocity scale for scalars [m s-1]
  real, allocatable, dimension(:,:,:) :: N         !< Brunt-Vaisala frequency [s-1]
  real, allocatable, dimension(:,:,:) :: N2        !< Squared Brunt-Vaisala frequency [s-2]
  real, allocatable, dimension(:,:,:) :: Vt2       !< Unresolved squared turbulence velocity for bulk Ri [m2 s-2]
  real, allocatable, dimension(:,:,:) :: Kt_KPP    !< Temp diffusivity from KPP [m2 s-1]
  real, allocatable, dimension(:,:,:) :: Ks_KPP    !< Scalar diffusivity from KPP [m2 s-1]
  real, allocatable, dimension(:,:,:) :: Kv_KPP    !< Viscosity due to KPP [m2 s-1]
  real, allocatable, dimension(:,:)   :: Tsurf     !< Temperature of surface layer [degC]
  real, allocatable, dimension(:,:)   :: Ssurf     !< Salinity of surface layer [ppt]
  real, allocatable, dimension(:,:)   :: Usurf     !< i-velocity of surface layer [m s-1]
  real, allocatable, dimension(:,:)   :: Vsurf     !< j-velocity of surface layer [m s-1]
  real, allocatable, dimension(:,:,:) :: EnhK      !< Enhancement for mixing coefficient
  real, allocatable, dimension(:,:,:) :: EnhVt2    !< Enhancement for Vt2

end type KPP_CS

!>@{ CPU time clocks
integer :: id_clock_KPP_calc, id_clock_KPP_compute_BLD, id_clock_KPP_smoothing
!>@}

#define __DO_SAFETY_CHECKS__

contains

!> Initialize the CVMix KPP module and set up diagnostics
!! Returns True if KPP is to be used, False otherwise.
logical function KPP_init(paramFile, G, GV, US, diag, Time, CS, passive, Waves)

  ! Arguments
  type(param_file_type),   intent(in)    :: paramFile !< File parser
  type(ocean_grid_type),   intent(in)    :: G         !< Ocean grid
  type(verticalGrid_type), intent(in)    :: GV        !< Vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US        !< A dimensional unit scaling type
  type(diag_ctrl), target, intent(in)    :: diag      !< Diagnostics
  type(time_type),         intent(in)    :: Time      !< Model time
  type(KPP_CS),            pointer       :: CS        !< Control structure
  logical,       optional, intent(out)   :: passive   !< Copy of %passiveMode
  type(wave_parameters_CS), optional, pointer :: Waves !< Wave CS

  ! Local variables
# include "version_variable.h"
  character(len=40) :: mdl = 'MOM_CVMix_KPP' !< name of this module
  character(len=20) :: string          !< local temporary string
  logical :: CS_IS_ONE=.false.         !< Logical for setting Cs based on Non-local
  logical :: lnoDGat1=.false.          !< True => G'(1) = 0 (shape function)
                                       !! False => compute G'(1) as in LMD94
  if (associated(CS)) call MOM_error(FATAL, 'MOM_CVMix_KPP, KPP_init: '// &
           'Control structure has already been initialized')

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
  allocate(CS)

  call openParameterBlock(paramFile,'KPP')
  call get_param(paramFile, mdl, 'PASSIVE', CS%passiveMode,           &
                 'If True, puts KPP into a passive-diagnostic mode.', &
                  default=.False.)
  !BGR: Note using PASSIVE for KPP creates warning for PASSIVE from Convection
  !     should we create a separate flag?
  if (present(passive)) passive=CS%passiveMode ! This is passed back to the caller so
                                               ! the caller knows to not use KPP output
  call get_param(paramFile, mdl, 'APPLY_NONLOCAL_TRANSPORT', CS%applyNonLocalTrans,  &
                 'If True, applies the non-local transport to heat and scalars. '//  &
                 'If False, calculates the non-local transport and tendencies but '//&
                 'purely for diagnostic purposes.',                                   &
                 default=.not. CS%passiveMode)
  call get_param(paramFile, mdl, 'N_SMOOTH', CS%n_smooth,  &
                 'The number of times the 1-1-4-1-1 Laplacian filter is applied on '//  &
                 'OBL depth.',   &
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
                 units='m',default=0.)
  call get_param(paramFile, mdl, 'FIXED_OBLDEPTH', CS%fixedOBLdepth,       &
                 'If True, fix the OBL depth to FIXED_OBLDEPTH_VALUE '//  &
                 'rather than using the OBL depth from CVMix. '//         &
                 'This option is just for testing purposes.',              &
                 default=.False.)
  call get_param(paramFile, mdl, 'FIXED_OBLDEPTH_VALUE', CS%fixedOBLdepth_value,  &
                 'Value for the fixed OBL depth when fixedOBLdepth==True. '//   &
                 'This parameter is for just for testing purposes. '//          &
                 'It will over-ride the OBLdepth computed from CVMix.',           &
                 units='m',default=30.0)
  call get_param(paramFile, mdl, 'SURF_LAYER_EXTENT', CS%surf_layer_ext,   &
                 'Fraction of OBL depth considered in the surface layer.', &
                 units='nondim',default=0.10)
  call get_param(paramFile, mdl, 'MINIMUM_OBL_DEPTH', CS%minOBLdepth,                            &
                 'If non-zero, a minimum depth to use for KPP OBL depth. Independent of '//     &
                 'this parameter, the OBL depth is always at least as deep as the first layer.', &
                 units='m',default=0.)
  call get_param(paramFile, mdl, 'MINIMUM_VT2', CS%minVtsqr,                                   &
                 'Min of the unresolved velocity Vt2 used in Rib CVMix calculation.\n'//  &
                 'Scaling: MINIMUM_VT2 = const1*d*N*ws, with d=1m, N=1e-5/s, ws=1e-6 m/s.',    &
                 units='m2/s2',default=1e-10)

! smg: for removal below
  call get_param(paramFile, mdl, 'CORRECT_SURFACE_LAYER_AVERAGE', CS%correctSurfLayerAvg,   &
                 'If true, applies a correction step to the averaging of surface layer '// &
                 'properties. This option is obsolete.', default=.False.)
  if (CS%correctSurfLayerAvg) &
    call MOM_error(FATAL,'Correct surface layer average disabled in code.  To recover \n'// &
                       ' feature will require code intervention.')
  call get_param(paramFile, mdl, 'FIRST_GUESS_SURFACE_LAYER_DEPTH', CS%surfLayerDepth,              &
                 'The first guess at the depth of the surface layer used for averaging '//         &
                 'the surface layer properties. If =0, the top model level properties '//          &
                 'will be used for the surface layer. If CORRECT_SURFACE_LAYER_AVERAGE=True, a '// &
                 'subsequent correction is applied. This parameter is obsolete', units='m', default=0.)
! smg: for removal above

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
                 units='m', default=0.)

!/BGR: New options for including Langmuir effects
!/ 1. Options related to enhancing the mixing coefficient
  call get_param(paramFile, mdl, "USE_KPP_LT_K", CS%LT_K_Enhancement, &
       'Flag for Langmuir turbulence enhancement of turbulent'//&
       'mixing coefficient.', units="", Default=.false.)
  call get_param(paramFile, mdl, "STOKES_MIXING", CS%STOKES_MIXING, &
       'Flag for Langmuir turbulence enhancement of turbulent'//&
       'mixing coefficient.', units="", Default=.false.)
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
                   'Method to enhance mixing coefficient in KPP. '//           &
                   'Valid options are: \n'//                                     &
                   '\t CONSTANT = Constant value (KPP_K_ENH_FAC) \n'//           &
                   '\t VR12     = Function of Langmuir number based on VR12\n'// &
                   '\t RW16     = Function of Langmuir number based on RW16',    &
                   default='CONSTANT')
    select case ( trim(string))
      case ("CONSTANT") ; CS%LT_K_METHOD = LT_K_MODE_CONSTANT
      case ("VR12")     ; CS%LT_K_METHOD = LT_K_MODE_VR12
      case ("RW16")     ; CS%LT_K_METHOD = LT_K_MODE_RW16
      case default      ; call MOM_error(FATAL,"KPP_init: "//&
                    "Unrecognized KPP_LT_K_METHOD option: "//trim(string))
    end select
    if (CS%LT_K_METHOD==LT_K_MODE_CONSTANT) then
      call get_param(paramFile, mdl, "KPP_K_ENH_FAC",CS%KPP_K_ENH_FAC ,     &
                   'Constant value to enhance mixing coefficient in KPP.',  &
                   default=1.0)
    endif
  endif
!/ 2. Options related to enhancing the unresolved Vt2/entrainment in Rib
  call get_param(paramFile, mdl, "USE_KPP_LT_VT2", CS%LT_Vt2_Enhancement, &
       'Flag for Langmuir turbulence enhancement of Vt2'//&
       'in Bulk Richardson Number.', units="", Default=.false.)
  if (CS%LT_Vt2_Enhancement) then
    call get_param(paramFile, mdl, "KPP_LT_VT2_METHOD",string ,                  &
                   'Method to enhance Vt2 in KPP. '//                          &
                   'Valid options are: \n'//                                     &
                   '\t CONSTANT = Constant value (KPP_VT2_ENH_FAC) \n'//         &
                   '\t VR12     = Function of Langmuir number based on VR12\n'// &
                   '\t RW16     = Function of Langmuir number based on RW16\n'// &
                   '\t LF17     = Function of Langmuir number based on LF17',    &
                   default='CONSTANT')
    select case ( trim(string))
      case ("CONSTANT") ; CS%LT_VT2_METHOD = LT_VT2_MODE_CONSTANT
      case ("VR12")     ; CS%LT_VT2_METHOD = LT_VT2_MODE_VR12
      case ("RW16")     ; CS%LT_VT2_METHOD = LT_VT2_MODE_RW16
      case ("LF17")     ; CS%LT_VT2_METHOD = LT_VT2_MODE_LF17
      case default      ; call MOM_error(FATAL,"KPP_init: "//&
                    "Unrecognized KPP_LT_VT2_METHOD option: "//trim(string))
    end select
    if (CS%LT_VT2_METHOD==LT_VT2_MODE_CONSTANT) then
      call get_param(paramFile, mdl, "KPP_VT2_ENH_FAC",CS%KPP_VT2_ENH_FAC ,     &
                   'Constant value to enhance VT2 in KPP.',  &
                   default=1.0)
    endif
  endif

  call closeParameterBlock(paramFile)
  call get_param(paramFile, mdl, 'DEBUG', CS%debug, default=.False., do_not_log=.True.)

  call CVMix_init_kpp( Ri_crit=CS%Ri_crit,                 &
                       minOBLdepth=CS%minOBLdepth,         &
                       minVtsqr=CS%minVtsqr,               &
                       vonKarman=CS%vonKarman,             &
                       surf_layer_ext=CS%surf_layer_ext,   &
                       interp_type=CS%interpType,          &
                       interp_type2=CS%interpType2,        &
                       lEkman=CS%computeEkman,             &
                       lMonOb=CS%computeMoninObukhov,      &
                       MatchTechnique=CS%MatchTechnique,   &
                       lenhanced_diff=CS%enhance_diffusion,&
                       lnonzero_surf_nonlocal=Cs_is_one   ,&
                       lnoDGat1=lnoDGat1                  ,&
                       CVMix_kpp_params_user=CS%KPP_params )

  ! Register diagnostics
  CS%diag => diag
  CS%id_OBLdepth = register_diag_field('ocean_model', 'KPP_OBLdepth', diag%axesT1, Time, &
      'Thickness of the surface Ocean Boundary Layer calculated by [CVMix] KPP', 'meter', &
      cmor_field_name='oml', cmor_long_name='ocean_mixed_layer_thickness_defined_by_mixing_scheme', &
      cmor_units='m', cmor_standard_name='Ocean Mixed Layer Thickness Defined by Mixing Scheme')
      ! CMOR names are placeholders; must be modified by time period
      ! for CMOR compliance. Diag manager will be used for omlmax and
      ! omldamax.
  if (CS%n_smooth > 0) then
    CS%id_OBLdepth_original = register_diag_field('ocean_model', 'KPP_OBLdepth_original', diag%axesT1, Time, &
        'Thickness of the surface Ocean Boundary Layer without smoothing calculated by [CVMix] KPP', 'meter', &
        cmor_field_name='oml', cmor_long_name='ocean_mixed_layer_thickness_defined_by_mixing_scheme', &
        cmor_units='m', cmor_standard_name='Ocean Mixed Layer Thickness Defined by Mixing Scheme')
  endif
  CS%id_BulkDrho = register_diag_field('ocean_model', 'KPP_BulkDrho', diag%axesTL, Time, &
      'Bulk difference in density used in Bulk Richardson number, as used by [CVMix] KPP', &
      'kg/m3', conversion=US%R_to_kg_m3)
  CS%id_BulkUz2 = register_diag_field('ocean_model', 'KPP_BulkUz2', diag%axesTL, Time, &
      'Square of bulk difference in resolved velocity used in Bulk Richardson number via [CVMix] KPP', 'm2/s2')
  CS%id_BulkRi = register_diag_field('ocean_model', 'KPP_BulkRi', diag%axesTL, Time, &
      'Bulk Richardson number used to find the OBL depth used by [CVMix] KPP', 'nondim')
  CS%id_Sigma = register_diag_field('ocean_model', 'KPP_sigma', diag%axesTi, Time, &
      'Sigma coordinate used by [CVMix] KPP', 'nondim')
  CS%id_Ws = register_diag_field('ocean_model', 'KPP_Ws', diag%axesTL, Time, &
      'Turbulent vertical velocity scale for scalars used by [CVMix] KPP', 'm/s')
  CS%id_N = register_diag_field('ocean_model', 'KPP_N', diag%axesTi, Time, &
      '(Adjusted) Brunt-Vaisala frequency used by [CVMix] KPP', '1/s')
  CS%id_N2 = register_diag_field('ocean_model', 'KPP_N2', diag%axesTi, Time, &
      'Square of Brunt-Vaisala frequency used by [CVMix] KPP', '1/s2')
  CS%id_Vt2 = register_diag_field('ocean_model', 'KPP_Vt2', diag%axesTL, Time, &
      'Unresolved shear turbulence used by [CVMix] KPP', 'm2/s2')
  CS%id_uStar = register_diag_field('ocean_model', 'KPP_uStar', diag%axesT1, Time, &
      'Friction velocity, u*, as used by [CVMix] KPP', 'm/s', conversion=US%Z_to_m*US%s_to_T)
  CS%id_buoyFlux = register_diag_field('ocean_model', 'KPP_buoyFlux', diag%axesTi, Time, &
      'Surface (and penetrating) buoyancy flux, as used by [CVMix] KPP', 'm2/s3', conversion=US%L_to_m**2*US%s_to_T**3)
  CS%id_QminusSW = register_diag_field('ocean_model', 'KPP_QminusSW', diag%axesT1, Time, &
      'Net temperature flux ignoring short-wave, as used by [CVMix] KPP', 'K m/s')
  CS%id_netS = register_diag_field('ocean_model', 'KPP_netSalt', diag%axesT1, Time, &
      'Effective net surface salt flux, as used by [CVMix] KPP', 'ppt m/s')
  CS%id_Kt_KPP = register_diag_field('ocean_model', 'KPP_Kheat', diag%axesTi, Time, &
      'Heat diffusivity due to KPP, as calculated by [CVMix] KPP', 'm2/s')
  CS%id_Kd_in = register_diag_field('ocean_model', 'KPP_Kd_in', diag%axesTi, Time, &
      'Diffusivity passed to KPP', 'm2/s', conversion=US%Z2_T_to_m2_s)
  CS%id_Ks_KPP = register_diag_field('ocean_model', 'KPP_Ksalt', diag%axesTi, Time, &
      'Salt diffusivity due to KPP, as calculated by [CVMix] KPP', 'm2/s')
  CS%id_Kv_KPP = register_diag_field('ocean_model', 'KPP_Kv', diag%axesTi, Time, &
      'Vertical viscosity due to KPP, as calculated by [CVMix] KPP', 'm2/s')
  CS%id_NLTt = register_diag_field('ocean_model', 'KPP_NLtransport_heat', diag%axesTi, Time, &
      'Non-local transport (Cs*G(sigma)) for heat, as calculated by [CVMix] KPP', 'nondim')
  CS%id_NLTs = register_diag_field('ocean_model', 'KPP_NLtransport_salt', diag%axesTi, Time, &
      'Non-local tranpsort (Cs*G(sigma)) for scalars, as calculated by [CVMix] KPP', 'nondim')
  CS%id_NLT_dTdt = register_diag_field('ocean_model', 'KPP_NLT_dTdt', diag%axesTL, Time, &
      'Temperature tendency due to non-local transport of heat, as calculated by [CVMix] KPP', 'K/s')
  CS%id_NLT_dSdt = register_diag_field('ocean_model', 'KPP_NLT_dSdt', diag%axesTL, Time, &
      'Salinity tendency due to non-local transport of salt, as calculated by [CVMix] KPP', 'ppt/s')
  CS%id_NLT_temp_budget = register_diag_field('ocean_model', 'KPP_NLT_temp_budget', diag%axesTL, Time, &
      'Heat content change due to non-local transport, as calculated by [CVMix] KPP', 'W/m^2')
  CS%id_NLT_saln_budget = register_diag_field('ocean_model', 'KPP_NLT_saln_budget', diag%axesTL, Time, &
      'Salt content change due to non-local transport, as calculated by [CVMix] KPP', 'kg/(sec*m^2)')
  CS%id_Tsurf = register_diag_field('ocean_model', 'KPP_Tsurf', diag%axesT1, Time, &
      'Temperature of surface layer (10% of OBL depth) as passed to [CVMix] KPP', 'C')
  CS%id_Ssurf = register_diag_field('ocean_model', 'KPP_Ssurf', diag%axesT1, Time, &
      'Salinity of surface layer (10% of OBL depth) as passed to [CVMix] KPP', 'ppt')
  CS%id_Usurf = register_diag_field('ocean_model', 'KPP_Usurf', diag%axesCu1, Time, &
      'i-component flow of surface layer (10% of OBL depth) as passed to [CVMix] KPP', 'm/s')
  CS%id_Vsurf = register_diag_field('ocean_model', 'KPP_Vsurf', diag%axesCv1, Time, &
      'j-component flow of surface layer (10% of OBL depth) as passed to [CVMix] KPP', 'm/s')
  CS%id_EnhK = register_diag_field('ocean_model', 'EnhK', diag%axesTI, Time, &
      'Langmuir number enhancement to K as used by [CVMix] KPP','nondim')
  CS%id_EnhVt2 = register_diag_field('ocean_model', 'EnhVt2', diag%axesTL, Time, &
      'Langmuir number enhancement to Vt2 as used by [CVMix] KPP','nondim')
  CS%id_La_SL = register_diag_field('ocean_model', 'KPP_La_SL', diag%axesT1, Time, &
      'Surface-layer Langmuir number computed in [CVMix] KPP','nondim')

  allocate( CS%N( SZI_(G), SZJ_(G),SZK_(GV)+1 ) )
  CS%N(:,:,:) = 0.
  allocate( CS%OBLdepth( SZI_(G), SZJ_(G) ) )
  CS%OBLdepth(:,:) = 0.
  allocate( CS%kOBL( SZI_(G), SZJ_(G) ) )
  CS%kOBL(:,:) = 0.
  allocate( CS%La_SL( SZI_(G), SZJ_(G) ) )
  CS%La_SL(:,:) = 0.
  allocate( CS%Vt2( SZI_(G), SZJ_(G),SZK_(GV) ) )
  CS%Vt2(:,:,:) = 0.
  if (CS%id_OBLdepth_original > 0) allocate( CS%OBLdepth_original( SZI_(G), SZJ_(G) ) )

  allocate( CS%OBLdepthprev( SZI_(G), SZJ_(G) ) ) ; CS%OBLdepthprev(:,:) = 0.0
  if (CS%id_BulkDrho > 0) allocate( CS%dRho( SZI_(G), SZJ_(G),SZK_(GV) ) )
  if (CS%id_BulkDrho > 0) CS%dRho(:,:,:) = 0.
  if (CS%id_BulkUz2 > 0)  allocate( CS%Uz2( SZI_(G), SZJ_(G),SZK_(GV) ) )
  if (CS%id_BulkUz2 > 0)  CS%Uz2(:,:,:) = 0.
  if (CS%id_BulkRi > 0)   allocate( CS%BulkRi( SZI_(G), SZJ_(G),SZK_(GV) ) )
  if (CS%id_BulkRi > 0)   CS%BulkRi(:,:,:) = 0.
  if (CS%id_Sigma > 0)    allocate( CS%sigma( SZI_(G), SZJ_(G),SZK_(GV)+1 ) )
  if (CS%id_Sigma > 0)    CS%sigma(:,:,:) = 0.
  if (CS%id_Ws > 0)       allocate( CS%Ws( SZI_(G), SZJ_(G),SZK_(GV) ) )
  if (CS%id_Ws > 0)       CS%Ws(:,:,:) = 0.
  if (CS%id_N2 > 0)       allocate( CS%N2( SZI_(G), SZJ_(G),SZK_(GV)+1 ) )
  if (CS%id_N2 > 0)       CS%N2(:,:,:) = 0.
  if (CS%id_Kt_KPP > 0)   allocate( CS%Kt_KPP( SZI_(G), SZJ_(G),SZK_(GV)+1 ) )
  if (CS%id_Kt_KPP > 0)   CS%Kt_KPP(:,:,:) = 0.
  if (CS%id_Ks_KPP > 0)   allocate( CS%Ks_KPP( SZI_(G), SZJ_(G),SZK_(GV)+1 ) )
  if (CS%id_Ks_KPP > 0)   CS%Ks_KPP(:,:,:) = 0.
  if (CS%id_Kv_KPP > 0)   allocate( CS%Kv_KPP( SZI_(G), SZJ_(G),SZK_(GV)+1 ) )
  if (CS%id_Kv_KPP > 0)   CS%Kv_KPP(:,:,:) = 0.
  if (CS%id_Tsurf > 0)    allocate( CS%Tsurf( SZI_(G), SZJ_(G)) )
  if (CS%id_Tsurf > 0)    CS%Tsurf(:,:) = 0.
  if (CS%id_Ssurf > 0)    allocate( CS%Ssurf( SZI_(G), SZJ_(G)) )
  if (CS%id_Ssurf > 0)    CS%Ssurf(:,:) = 0.
  if (CS%id_Usurf > 0)    allocate( CS%Usurf( SZIB_(G), SZJ_(G)) )
  if (CS%id_Usurf > 0)    CS%Usurf(:,:) = 0.
  if (CS%id_Vsurf > 0)    allocate( CS%Vsurf( SZI_(G), SZJB_(G)) )
  if (CS%id_Vsurf > 0)    CS%Vsurf(:,:) = 0.
  if (CS%id_EnhVt2 > 0)    allocate( CS%EnhVt2( SZI_(G), SZJ_(G),SZK_(GV)) )
  if (CS%id_EnhVt2 > 0)    CS%EnhVt2(:,:,:) = 0.
  if (CS%id_EnhK > 0)    allocate( CS%EnhK( SZI_(G), SZJ_(G),SZK_(GV)+1 ) )
  if (CS%id_EnhK > 0)    CS%EnhK(:,:,:) = 0.

  id_clock_KPP_calc = cpu_clock_id('Ocean KPP calculate)', grain=CLOCK_MODULE)
  id_clock_KPP_compute_BLD = cpu_clock_id('(Ocean KPP comp BLD)', grain=CLOCK_ROUTINE)

end function KPP_init

!> KPP vertical diffusivity/viscosity and non-local tracer transport
subroutine KPP_calculate(CS, G, GV, US, h, uStar, &
                         buoyFlux, Kt, Ks, Kv, nonLocalTransHeat,&
                         nonLocalTransScalar, waves)

  ! Arguments
  type(KPP_CS),                                pointer       :: CS    !< Control structure
  type(ocean_grid_type),                       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type),                     intent(in)    :: GV    !< Ocean vertical grid
  type(unit_scale_type),                       intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h     !< Layer/level thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),            intent(in)    :: uStar !< Surface friction velocity [Z T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: buoyFlux !< Surface buoyancy flux [L2 T-3 ~> m2 s-3]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kt  !< (in)  Vertical diffusivity of heat w/o KPP
                                                                    !! (out) Vertical diffusivity including KPP
                                                                    !!       [Z2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Ks  !< (in)  Vertical diffusivity of salt w/o KPP
                                                                    !! (out) Vertical diffusivity including KPP
                                                                    !!       [Z2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kv  !< (in)  Vertical viscosity w/o KPP
                                                                    !! (out) Vertical viscosity including KPP
                                                                    !!       [Z2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: nonLocalTransHeat   !< Temp non-local transport [m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: nonLocalTransScalar !< scalar non-local trans. [m s-1]
  type(wave_parameters_CS),          optional, pointer       :: Waves !< Wave CS

! Local variables
  integer :: i, j, k                             ! Loop indices
  real, dimension( GV%ke )     :: cellHeight     ! Cell center heights referenced to surface [m] (negative in ocean)
  real, dimension( GV%ke+1 )   :: iFaceHeight    ! Interface heights referenced to surface [m] (negative in ocean)
  real, dimension( GV%ke+1, 2) :: Kdiffusivity   ! Vertical diffusivity at interfaces [m2 s-1]
  real, dimension( GV%ke+1 )   :: Kviscosity     ! Vertical viscosity at interfaces [m2 s-1]
  real, dimension( GV%ke+1, 2) :: nonLocalTrans  ! Non-local transport for heat/salt at interfaces [nondim]

  real :: surfFricVel, surfBuoyFlux
  real :: sigma, sigmaRatio
  real :: buoy_scale ! A unit conversion factor for buoyancy fluxes [m2 T3 L-2 s-3 ~> nondim]
  real :: dh    ! The local thickness used for calculating interface positions [m]
  real :: hcorr ! A cumulative correction arising from inflation of vanished layers [m]

  ! For Langmuir Calculations
  real :: LangEnhK     ! Langmuir enhancement for mixing coefficient


  if (CS%debug) then
    call hchksum(h, "KPP in: h",G%HI,haloshift=0, scale=GV%H_to_m)
    call hchksum(uStar, "KPP in: uStar",G%HI,haloshift=0, scale=US%Z_to_m*US%s_to_T)
    call hchksum(buoyFlux, "KPP in: buoyFlux",G%HI,haloshift=0)
    call hchksum(Kt, "KPP in: Kt",G%HI,haloshift=0, scale=US%Z2_T_to_m2_s)
    call hchksum(Ks, "KPP in: Ks",G%HI,haloshift=0, scale=US%Z2_T_to_m2_s)
  endif

  nonLocalTrans(:,:) = 0.0

  if (CS%id_Kd_in > 0) call post_data(CS%id_Kd_in, Kt, CS%diag)

  call cpu_clock_begin(id_clock_KPP_calc)
  buoy_scale = US%L_to_m**2*US%s_to_T**3

  !$OMP parallel do default(none) firstprivate(nonLocalTrans)                               &
  !$OMP                           private(surfFricVel, iFaceHeight, hcorr, dh, cellHeight,  &
  !$OMP                           surfBuoyFlux, Kdiffusivity, Kviscosity, LangEnhK, sigma,  &
  !$OMP                           sigmaRatio)                                               &
  !$OMP                           shared(G, GV, CS, US, uStar, h, buoy_scale, buoyFlux, Kt, &
  !$OMP                           Ks, Kv, nonLocalTransHeat, nonLocalTransScalar, waves)
  ! loop over horizontal points on processor
  do j = G%jsc, G%jec
    do i = G%isc, G%iec

      ! skip calling KPP for land points
      if (G%mask2dT(i,j)==0.) cycle

      ! things independent of position within the column
      surfFricVel = US%Z_to_m*US%s_to_T * uStar(i,j)

      iFaceHeight(1) = 0.0 ! BBL is all relative to the surface
      hcorr = 0.
      do k=1,GV%ke

        ! cell center and cell bottom in meters (negative values in the ocean)
        dh = h(i,j,k) * GV%H_to_m ! Nominal thickness to use for increment
        dh = dh + hcorr ! Take away the accumulated error (could temporarily make dh<0)
        hcorr = min( dh - CS%min_thickness, 0. ) ! If inflating then hcorr<0
        dh = max( dh, CS%min_thickness ) ! Limit increment dh>=min_thickness
        cellHeight(k)    = iFaceHeight(k) - 0.5 * dh
        iFaceHeight(k+1) = iFaceHeight(k) - dh

      enddo ! k-loop finishes

      surfBuoyFlux = buoy_scale*buoyFlux(i,j,1) ! This is only used in kpp_compute_OBL_depth to limit
                                     ! h to Monin-Obukov (default is false, ie. not used)

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
         Kdiffusivity(:,1) = US%Z2_T_to_m2_s * Kt(i,j,:)
         Kdiffusivity(:,2) = US%Z2_T_to_m2_s * Ks(i,j,:)
         Kviscosity(:) = US%Z2_T_to_m2_s * Kv(i,j,:)
      endif

      call CVMix_coeffs_kpp(Kviscosity(:),     & ! (inout) Total viscosity [m2 s-1]
                            Kdiffusivity(:,1), & ! (inout) Total heat diffusivity [m2 s-1]
                            Kdiffusivity(:,2), & ! (inout) Total salt diffusivity [m2 s-1]
                            iFaceHeight,       & ! (in) Height of interfaces [m]
                            cellHeight,        & ! (in) Height of level centers [m]
                            Kviscosity(:),     & ! (in) Original viscosity [m2 s-1]
                            Kdiffusivity(:,1), & ! (in) Original heat diffusivity [m2 s-1]
                            Kdiffusivity(:,2), & ! (in) Original salt diffusivity [m2 s-1]
                            CS%OBLdepth(i,j),  & ! (in) OBL depth [m]
                            CS%kOBL(i,j),      & ! (in) level (+fraction) of OBL extent
                            nonLocalTrans(:,1),& ! (out) Non-local heat transport [nondim]
                            nonLocalTrans(:,2),& ! (out) Non-local salt transport [nondim]
                            surfFricVel,       & ! (in) Turbulent friction velocity at surface [m s-1]
                            surfBuoyFlux,      & ! (in) Buoyancy flux at surface [m2 s-3]
                            GV%ke,             & ! (in) Number of levels to compute coeffs for
                            GV%ke,             & ! (in) Number of levels in array shape
                            CVMix_kpp_params_user=CS%KPP_params )

      ! safety check, Kviscosity and Kdiffusivity must be >= 0
      do k=1, GV%ke+1
        if (Kviscosity(k) < 0. .or. Kdiffusivity(k,1) < 0.) then
          call MOM_error(FATAL,"KPP_calculate, after CVMix_coeffs_kpp: "// &
                   "Negative vertical viscosity or diffusivity has been detected. " // &
                   "This is likely related to the choice of MATCH_TECHNIQUE and INTERP_TYPE2." //&
                   "You might consider using the default options for these parameters." )
        endif
      enddo

      IF (CS%LT_K_ENHANCEMENT) then
        if (CS%LT_K_METHOD==LT_K_MODE_CONSTANT) then
           LangEnhK = CS%KPP_K_ENH_FAC
        elseif (CS%LT_K_METHOD==LT_K_MODE_VR12) then
           ! Added minimum value for La_SL, so removed maximum value for LangEnhK.
           LangEnhK = sqrt(1.+(1.5*CS%La_SL(i,j))**(-2) + &
                (5.4*CS%La_SL(i,j))**(-4))
        elseif (CS%LT_K_METHOD==LT_K_MODE_RW16) then
          !This maximum value is proposed in Reichl et al., 2016 JPO formula
          LangEnhK = min(2.25, 1. + 1./CS%La_SL(i,j))
        else
           !This shouldn't be reached.
           !call MOM_error(WARNING,"Unexpected behavior in MOM_CVMix_KPP, see error in LT_K_ENHANCEMENT")
           LangEnhK = 1.0
        endif
        do k=1,GV%ke
          if (CS%LT_K_SHAPE== LT_K_CONSTANT) then
            if (CS%id_EnhK > 0) CS%EnhK(i,j,:) = LangEnhK
            Kdiffusivity(k,1) = Kdiffusivity(k,1) * LangEnhK
            Kdiffusivity(k,2) = Kdiffusivity(k,2) * LangEnhK
            Kviscosity(k)     = Kviscosity(k)   * LangEnhK
          elseif (CS%LT_K_SHAPE == LT_K_SCALED) then
            sigma = min(1.0,-iFaceHeight(k)/CS%OBLdepth(i,j))
            SigmaRatio = sigma * (1. - sigma)**2 / 0.148148037
            if (CS%id_EnhK > 0) CS%EnhK(i,j,k) = (1.0 + (LangEnhK - 1.)*sigmaRatio)
            Kdiffusivity(k,1) = Kdiffusivity(k,1) * ( 1. + &
                                ( LangEnhK - 1.)*sigmaRatio)
            Kdiffusivity(k,2) = Kdiffusivity(k,2) * ( 1. + &
                                ( LangEnhK - 1.)*sigmaRatio)
            Kviscosity(k) = Kviscosity(k) * ( 1. + &
                                ( LangEnhK - 1.)*sigmaRatio)
          endif
        enddo
      endif

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
      nonLocalTransHeat(i,j,:)   = nonLocalTrans(:,1) ! temp
      nonLocalTransScalar(i,j,:) = nonLocalTrans(:,2) ! saln

      ! set the KPP diffusivity and viscosity to zero for testing purposes
      if (CS%KPPzeroDiffusivity) then
         Kdiffusivity(:,1) = 0.0
         Kdiffusivity(:,2) = 0.0
         Kviscosity(:)     = 0.0
      endif


      ! compute unresolved squared velocity for diagnostics
      if (CS%id_Vt2 > 0) then
!BGR Now computing VT2 above so can modify for LT
!    therefore, don't repeat this operation here
!        CS%Vt2(i,j,:) = CVmix_kpp_compute_unresolved_shear( &
!                    cellHeight(1:GV%ke),                & ! Depth of cell center [m]
!                    ws_cntr=Ws_1d,                      & ! Turbulent velocity scale profile, at centers [m s-1]
!                    N_iface=CS%N(i,j,:),                & ! Buoyancy frequency at interface [s-1]
!                    CVmix_kpp_params_user=CS%KPP_params ) ! KPP parameters
      endif

      ! Copy 1d data into 3d diagnostic arrays
      !/ grabbing obldepth_0d for next time step.
      CS%OBLdepthprev(i,j)=CS%OBLdepth(i,j)
      if (CS%id_sigma > 0) then
        CS%sigma(i,j,:)  = 0.
        if (CS%OBLdepth(i,j)>0.)   CS%sigma(i,j,:)  = -iFaceHeight/CS%OBLdepth(i,j)
      endif
      if (CS%id_Kt_KPP > 0)   CS%Kt_KPP(i,j,:) = Kdiffusivity(:,1)
      if (CS%id_Ks_KPP > 0)   CS%Ks_KPP(i,j,:) = Kdiffusivity(:,2)
      if (CS%id_Kv_KPP > 0)   CS%Kv_KPP(i,j,:) = Kviscosity(:)

      ! Update output of routine
      if (.not. CS%passiveMode) then
        if (CS%KPPisAdditive) then
          do k=1, GV%ke+1
            Kt(i,j,k) = Kt(i,j,k) + US%m2_s_to_Z2_T * Kdiffusivity(k,1)
            Ks(i,j,k) = Ks(i,j,k) + US%m2_s_to_Z2_T * Kdiffusivity(k,2)
            Kv(i,j,k) = Kv(i,j,k) + US%m2_s_to_Z2_T * Kviscosity(k)
            if (CS%Stokes_Mixing) Waves%KvS(i,j,k) = Kv(i,j,k)
          enddo
        else ! KPP replaces prior diffusivity when former is non-zero
          do k=1, GV%ke+1
            if (Kdiffusivity(k,1) /= 0.) Kt(i,j,k) = US%m2_s_to_Z2_T * Kdiffusivity(k,1)
            if (Kdiffusivity(k,2) /= 0.) Ks(i,j,k) = US%m2_s_to_Z2_T * Kdiffusivity(k,2)
            if (Kviscosity(k) /= 0.) Kv(i,j,k) = US%m2_s_to_Z2_T * Kviscosity(k)
            if (CS%Stokes_Mixing) Waves%KvS(i,j,k) = Kv(i,j,k)
          enddo
        endif
      endif


    ! end of the horizontal do-loops over the vertical columns
    enddo ! i
  enddo ! j

  call cpu_clock_end(id_clock_KPP_calc)

  if (CS%debug) then
    call hchksum(Kt, "KPP out: Kt", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
    call hchksum(Ks, "KPP out: Ks", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
  endif

  ! send diagnostics to post_data
  if (CS%id_OBLdepth > 0) call post_data(CS%id_OBLdepth, CS%OBLdepth,        CS%diag)
  if (CS%id_OBLdepth_original > 0) call post_data(CS%id_OBLdepth_original,CS%OBLdepth_original,CS%diag)
  if (CS%id_sigma    > 0) call post_data(CS%id_sigma,    CS%sigma,           CS%diag)
  if (CS%id_Ws       > 0) call post_data(CS%id_Ws,       CS%Ws,              CS%diag)
  if (CS%id_Vt2      > 0) call post_data(CS%id_Vt2,      CS%Vt2,             CS%diag)
  if (CS%id_uStar    > 0) call post_data(CS%id_uStar,    uStar,              CS%diag)
  if (CS%id_buoyFlux > 0) call post_data(CS%id_buoyFlux, buoyFlux,           CS%diag)
  if (CS%id_Kt_KPP   > 0) call post_data(CS%id_Kt_KPP,   CS%Kt_KPP,          CS%diag)
  if (CS%id_Ks_KPP   > 0) call post_data(CS%id_Ks_KPP,   CS%Ks_KPP,          CS%diag)
  if (CS%id_Kv_KPP   > 0) call post_data(CS%id_Kv_KPP,   CS%Kv_KPP,          CS%diag)
  if (CS%id_NLTt     > 0) call post_data(CS%id_NLTt,     nonLocalTransHeat,  CS%diag)
  if (CS%id_NLTs     > 0) call post_data(CS%id_NLTs,     nonLocalTransScalar,CS%diag)


end subroutine KPP_calculate


!> Compute OBL depth
subroutine KPP_compute_BLD(CS, G, GV, US, h, Temp, Salt, u, v, tv, uStar, buoyFlux, Waves)

  ! Arguments
  type(KPP_CS),                               pointer       :: CS    !< Control structure
  type(ocean_grid_type),                      intent(inout) :: G     !< Ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV    !< Ocean vertical grid
  type(unit_scale_type),                      intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h     !< Layer/level thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: Temp  !< potential/cons temp [degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: Salt  !< Salinity [ppt]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u     !< Velocity i-component [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v     !< Velocity j-component [L T-1 ~> m s-1]
  type(thermo_var_ptrs),                      intent(in)    :: tv    !< Thermodynamics structure.
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: uStar !< Surface friction velocity [Z T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: buoyFlux !< Surface buoyancy flux [L2 T-3 ~> m2 s-3]
  type(wave_parameters_CS),         optional, pointer       :: Waves !< Wave CS

  ! Local variables
  integer :: i, j, k, km1                        ! Loop indices
  real, dimension( GV%ke )     :: cellHeight     ! Cell center heights referenced to surface [m] (negative in ocean)
  real, dimension( GV%ke+1 )   :: iFaceHeight    ! Interface heights referenced to surface [m] (negative in ocean)
  real, dimension( GV%ke+1 )   :: N2_1d          ! Brunt-Vaisala frequency squared, at interfaces [s-2]
  real, dimension( GV%ke )     :: Ws_1d          ! Profile of vertical velocity scale for scalars [m s-1]
  real, dimension( GV%ke )     :: deltaRho       ! delta Rho in numerator of Bulk Ri number [R ~> kg m-3]
  real, dimension( GV%ke )     :: deltaU2        ! square of delta U (shear) in denominator of Bulk Ri [m2 s-2]
  real, dimension( GV%ke )     :: surfBuoyFlux2
  real, dimension( GV%ke )     :: BulkRi_1d      ! Bulk Richardson number for each layer [nondim]

  ! for EOS calculation
  real, dimension( 3*GV%ke )   :: rho_1D   ! A column of densities [R ~> kg m-3]
  real, dimension( 3*GV%ke )   :: pres_1D  ! A column of pressures [R L2 T-2 ~> Pa]
  real, dimension( 3*GV%ke )   :: Temp_1D  ! A column of temperatures [degC]
  real, dimension( 3*GV%ke )   :: Salt_1D  ! A column of salinities [ppt]

  real :: surfFricVel, surfBuoyFlux, Coriolis
  real :: GoRho  ! Gravitational acceleration divided by density in MKS units [m R-1 s-2 ~> m4 kg-1 s-2]
  real :: pRef   ! The interface pressure [R L2 T-2 ~> Pa]
  real :: rho1, rhoK, Uk, Vk, sigma, sigmaRatio

  real :: zBottomMinusOffset   ! Height of bottom plus a little bit [m]
  real :: SLdepth_0d           ! Surface layer depth = surf_layer_ext*OBLdepth.
  real :: hTot                 ! Running sum of thickness used in the surface layer average [m]
  real :: buoy_scale           ! A unit conversion factor for buoyancy fluxes [m2 T3 L-2 s-3 ~> nondim]
  real :: delH                 ! Thickness of a layer [m]
  real :: surfHtemp, surfTemp  ! Integral and average of temp over the surface layer
  real :: surfHsalt, surfSalt  ! Integral and average of saln over the surface layer
  real :: surfHu, surfU        ! Integral and average of u over the surface layer
  real :: surfHv, surfV        ! Integral and average of v over the surface layer
  real :: dh    ! The local thickness used for calculating interface positions [m]
  real :: hcorr ! A cumulative correction arising from inflation of vanished layers [m]
  integer :: kk, ksfc, ktmp

  ! For Langmuir Calculations
  real :: LangEnhW     ! Langmuir enhancement for turbulent velocity scale
  real, dimension(GV%ke) :: LangEnhVt2   ! Langmuir enhancement for unresolved shear
  real, dimension(GV%ke) :: U_H, V_H
  real :: MLD_GUESS, LA
  real :: surfHuS, surfHvS, surfUs, surfVs, wavedir, currentdir
  real :: VarUp, VarDn, M, VarLo, VarAvg
  real :: H10pct, H20pct,CMNFACT, USx20pct, USy20pct, enhvt2
  integer :: B
  real :: WST


  if (CS%debug) then
    call hchksum(Salt, "KPP in: S",G%HI,haloshift=0)
    call hchksum(Temp, "KPP in: T",G%HI,haloshift=0)
    call hchksum(u, "KPP in: u",G%HI,haloshift=0)
    call hchksum(v, "KPP in: v",G%HI,haloshift=0)
  endif

  call cpu_clock_begin(id_clock_KPP_compute_BLD)

  ! some constants
  GoRho = US%L_T_to_m_s**2*US%m_to_Z * GV%g_Earth / GV%Rho0
  buoy_scale = US%L_to_m**2*US%s_to_T**3

  ! loop over horizontal points on processor
  !$OMP parallel do default(none) private(surfFricVel, iFaceHeight, hcorr, dh, cellHeight,  &
  !$OMP                           surfBuoyFlux, U_H, V_H, Coriolis, pRef, SLdepth_0d,       &
  !$OMP                           ksfc, surfHtemp, surfHsalt, surfHu, surfHv, surfHuS,      &
  !$OMP                           surfHvS, hTot, delH, surftemp, surfsalt, surfu, surfv,    &
  !$OMP                           surfUs, surfVs, Uk, Vk, deltaU2, km1, kk, pres_1D,        &
  !$OMP                           Temp_1D, salt_1D, surfBuoyFlux2, MLD_GUESS, LA, rho_1D,   &
  !$OMP                           deltarho, N2_1d, ws_1d, LangEnhVT2, enhvt2, wst,          &
  !$OMP                           BulkRi_1d, zBottomMinusOffset) &
  !$OMP                           shared(G, GV, CS, US, uStar, h, buoy_scale, buoyFlux,     &
  !$OMP                           Temp, Salt, waves, tv, GoRho, u, v)
  do j = G%jsc, G%jec
    do i = G%isc, G%iec

      ! skip calling KPP for land points
      if (G%mask2dT(i,j)==0.) cycle

      do k=1,GV%ke
        U_H(k) = 0.5 * US%L_T_to_m_s*(u(i,j,k)+u(i-1,j,k))
        V_H(k) = 0.5 * US%L_T_to_m_s*(v(i,j,k)+v(i,j-1,k))
      enddo

      ! things independent of position within the column
      Coriolis = 0.25*US%s_to_T*( (G%CoriolisBu(i,j)   + G%CoriolisBu(i-1,j-1)) + &
                                  (G%CoriolisBu(i-1,j) + G%CoriolisBu(i,j-1)) )
      surfFricVel = US%Z_to_m*US%s_to_T * uStar(i,j)

      ! Bullk Richardson number computed for each cell in a column,
      ! assuming OBLdepth = grid cell depth. After Rib(k) is
      ! known for the column, then CVMix interpolates to find
      ! the actual OBLdepth. This approach avoids need to iterate
      ! on the OBLdepth calculation. It follows that used in MOM5
      ! and POP.
      iFaceHeight(1) = 0.0 ! BBL is all relative to the surface
      pRef = 0. ; if (associated(tv%p_surf)) pRef = tv%p_surf(i,j)
      hcorr = 0.
      do k=1,GV%ke

        ! cell center and cell bottom in meters (negative values in the ocean)
        dh = h(i,j,k) * GV%H_to_m ! Nominal thickness to use for increment
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

        ! average temp, saln, u, v over surface layer
        ! use C-grid average to get u,v on T-points.
        surfHtemp=0.0
        surfHsalt=0.0
        surfHu   =0.0
        surfHv   =0.0
        surfHuS  =0.0
        surfHvS  =0.0
        hTot     =0.0
        do ktmp = 1,ksfc

          ! SLdepth_0d can be between cell interfaces
          delH = min( max(0.0, SLdepth_0d - hTot), h(i,j,ktmp)*GV%H_to_m )

          ! surface layer thickness
          hTot = hTot + delH

          ! surface averaged fields
          surfHtemp = surfHtemp + Temp(i,j,ktmp) * delH
          surfHsalt = surfHsalt + Salt(i,j,ktmp) * delH
          surfHu    = surfHu + 0.5*US%L_T_to_m_s*(u(i,j,ktmp)+u(i-1,j,ktmp)) * delH
          surfHv    = surfHv + 0.5*US%L_T_to_m_s*(v(i,j,ktmp)+v(i,j-1,ktmp)) * delH
          if (CS%Stokes_Mixing) then
            surfHus = surfHus + 0.5*(WAVES%US_x(i,j,ktmp)+WAVES%US_x(i-1,j,ktmp)) * delH
            surfHvs = surfHvs + 0.5*(WAVES%US_y(i,j,ktmp)+WAVES%US_y(i,j-1,ktmp)) * delH
          endif

        enddo
        surfTemp = surfHtemp / hTot
        surfSalt = surfHsalt / hTot
        surfU    = surfHu    / hTot
        surfV    = surfHv    / hTot
        surfUs   = surfHus   / hTot
        surfVs   = surfHvs   / hTot

        ! vertical shear between present layer and
        ! surface layer averaged surfU,surfV.
        ! C-grid average to get Uk and Vk on T-points.
        Uk         = 0.5*US%L_T_to_m_s*(u(i,j,k)+u(i-1,j,k)) - surfU
        Vk         = 0.5*US%L_T_to_m_s*(v(i,j,k)+v(i,j-1,k)) - surfV

        if (CS%Stokes_Mixing) then
          ! If momentum is mixed down the Stokes drift gradient, then
          !  the Stokes drift must be included in the bulk Richardson number
          !  calculation.
          Uk =  Uk + (0.5*(Waves%Us_x(i,j,k)+Waves%US_x(i-1,j,k)) -surfUs )
          Vk =  Vk + (0.5*(Waves%Us_y(i,j,k)+Waves%Us_y(i,j-1,k)) -surfVs )
        endif

        deltaU2(k) = Uk**2 + Vk**2

        ! pressure, temp, and saln for EOS
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

        ! this difference accounts for penetrating SW
        surfBuoyFlux2(k) = buoy_scale * (buoyFlux(i,j,1) - buoyFlux(i,j,k+1))

      enddo ! k-loop finishes

      if (CS%LT_K_ENHANCEMENT .or. CS%LT_VT2_ENHANCEMENT) then
        MLD_GUESS = max( 1.*US%m_to_Z, abs(US%m_to_Z*CS%OBLdepthprev(i,j) ) )
        call get_Langmuir_Number(LA, G, GV, US, MLD_guess, uStar(i,j), i, j, &
                                 H=H(i,j,:), U_H=U_H, V_H=V_H, WAVES=WAVES)
        CS%La_SL(i,j)=LA
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
        N2_1d(k)    = (GoRho * (rho_1D(kk+2) - rho_1D(kk+3)) ) / &
                      ((0.5*(h(i,j,km1) + h(i,j,k))+GV%H_subroundoff)*GV%H_to_m)
        CS%N(i,j,k)     = sqrt( max( N2_1d(k), 0.) )
      enddo
      N2_1d(GV%ke+1 ) = 0.0
      CS%N(i,j,GV%ke+1 )  = 0.0

      ! turbulent velocity scales w_s and w_m computed at the cell centers.
      ! Note that if sigma > CS%surf_layer_ext, then CVMix_kpp_compute_turbulent_scales
      ! computes w_s and w_m velocity scale at sigma=CS%surf_layer_ext. So we only pass
      ! sigma=CS%surf_layer_ext for this calculation.
      call CVMix_kpp_compute_turbulent_scales( &
        CS%surf_layer_ext, & ! (in)  Normalized surface layer depth; sigma = CS%surf_layer_ext
        -cellHeight,       & ! (in)  Assume here that OBL depth [m] = -cellHeight(k)
        surfBuoyFlux2,     & ! (in)  Buoyancy flux at surface [m2 s-3]
        surfFricVel,       & ! (in)  Turbulent friction velocity at surface [m s-1]
        w_s=Ws_1d,         & ! (out) Turbulent velocity scale profile [m s-1]
        CVMix_kpp_params_user=CS%KPP_params )

      !Compute CVMix VT2
      CS%Vt2(i,j,:) = CVmix_kpp_compute_unresolved_shear( &
                      zt_cntr=cellHeight(1:GV%ke),        & ! Depth of cell center [m]
                      ws_cntr=Ws_1d,                      & ! Turbulent velocity scale profile, at centers [m s-1]
                      N_iface=CS%N(i,j,:),                & ! Buoyancy frequency at interface [s-1]
                    CVmix_kpp_params_user=CS%KPP_params ) ! KPP parameters

      !Modify CVMix VT2
      IF (CS%LT_VT2_ENHANCEMENT) then
        IF (CS%LT_VT2_METHOD==LT_VT2_MODE_CONSTANT) then
          do k=1,GV%ke
            LangEnhVT2(k) = CS%KPP_VT2_ENH_FAC
          enddo
        elseif (CS%LT_VT2_METHOD==LT_VT2_MODE_VR12) then
          !Introduced minimum value for La_SL, so maximum value for enhvt2 is removed.
          enhvt2 = sqrt(1.+(1.5*CS%La_SL(i,j))**(-2) + &
                   (5.4*CS%La_SL(i,j))**(-4))
          do k=1,GV%ke
             LangEnhVT2(k) = enhvt2
          enddo
        elseif (CS%LT_VT2_METHOD==LT_VT2_MODE_RW16) then
          !Introduced minimum value for La_SL, so maximum value for enhvt2 is removed.
          enhvt2 = 1. + 2.3*CS%La_SL(i,j)**(-0.5)
          do k=1,GV%ke
            LangEnhVT2(k) = enhvt2
          enddo
        elseif (CS%LT_VT2_METHOD==LT_VT2_MODE_LF17) then
          CS%CS=cvmix_get_kpp_real('c_s',CS%KPP_params)
          do k=1,GV%ke
            WST = (max(0.,-buoy_scale*buoyflux(i,j,1))*(-cellHeight(k)))**(1./3.)
            LangEnhVT2(k) = sqrt((0.15*WST**3. + 0.17*surfFricVel**3.* &
                 (1.+0.49*CS%La_SL(i,j)**(-2.)))  / &
                 (0.2*ws_1d(k)**3/(CS%cs*CS%surf_layer_ext*CS%vonKarman**4.)))
          enddo
        else
           !This shouldn't be reached.
           !call MOM_error(WARNING,"Unexpected behavior in MOM_CVMix_KPP, see error in Vt2")
           LangEnhVT2(:) = 1.0
        endif
      else
        LangEnhVT2(:) = 1.0
      endif

      do k=1,GV%ke
        CS%Vt2(i,j,k)=CS%Vt2(i,j,k)*LangEnhVT2(k)
        if (CS%id_EnhVt2 > 0) CS%EnhVt2(i,j,k)=LangEnhVT2(k)
      enddo

      ! Calculate Bulk Richardson number from eq (21) of LMD94
      BulkRi_1d = CVmix_kpp_compute_bulk_Richardson( &
                  zt_cntr = cellHeight(1:GV%ke),     & ! Depth of cell center [m]
                  delta_buoy_cntr=GoRho*deltaRho,    & ! Bulk buoyancy difference, Br-B(z) [s-1]
                  delta_Vsqr_cntr=deltaU2,           & ! Square of resolved velocity difference [m2 s-2]
                  Vt_sqr_cntr=CS%Vt2(i,j,:),         &
                  ws_cntr=Ws_1d,                     & ! Turbulent velocity scale profile [m s-1]
                  N_iface=CS%N(i,j,:))               ! Buoyancy frequency [s-1]


      surfBuoyFlux = buoy_scale * buoyFlux(i,j,1) ! This is only used in kpp_compute_OBL_depth to limit
                                     ! h to Monin-Obukov (default is false, ie. not used)

      call CVMix_kpp_compute_OBL_depth( &
        BulkRi_1d,              & ! (in) Bulk Richardson number
        iFaceHeight,            & ! (in) Height of interfaces [m]
        CS%OBLdepth(i,j),       & ! (out) OBL depth [m]
        CS%kOBL(i,j),           & ! (out) level (+fraction) of OBL extent
        zt_cntr=cellHeight,     & ! (in) Height of cell centers [m]
        surf_fric=surfFricVel,  & ! (in) Turbulent friction velocity at surface [m s-1]
        surf_buoy=surfBuoyFlux, & ! (in) Buoyancy flux at surface [m2 s-3]
        Coriolis=Coriolis,      & ! (in) Coriolis parameter [s-1]
        CVMix_kpp_params_user=CS%KPP_params ) ! KPP parameters

      ! A hack to avoid KPP reaching the bottom. It was needed during development
      ! because KPP was unable to handle vanishingly small layers near the bottom.
      if (CS%deepOBLoffset>0.) then
        zBottomMinusOffset = iFaceHeight(GV%ke+1) + min(CS%deepOBLoffset,-0.1*iFaceHeight(GV%ke+1))
        CS%OBLdepth(i,j) = min( CS%OBLdepth(i,j), -zBottomMinusOffset )
      endif

      ! apply some constraints on OBLdepth
      if(CS%fixedOBLdepth)  CS%OBLdepth(i,j) = CS%fixedOBLdepth_value
      CS%OBLdepth(i,j) = max( CS%OBLdepth(i,j), -iFaceHeight(2) )       ! no shallower than top layer
      CS%OBLdepth(i,j) = min( CS%OBLdepth(i,j), -iFaceHeight(GV%ke+1) ) ! no deeper than bottom
      CS%kOBL(i,j)     = CVMix_kpp_compute_kOBL_depth( iFaceHeight, cellHeight, CS%OBLdepth(i,j) )


      ! recompute wscale for diagnostics, now that we in fact know boundary layer depth
      !BGR consider if LTEnhancement is wanted for diagnostics
      if (CS%id_Ws > 0) then
          call CVMix_kpp_compute_turbulent_scales( &
            -CellHeight/CS%OBLdepth(i,j),          & ! (in)  Normalized boundary layer coordinate
            CS%OBLdepth(i,j),                      & ! (in)  OBL depth [m]
            surfBuoyFlux,                          & ! (in)  Buoyancy flux at surface [m2 s-3]
            surfFricVel,                           & ! (in)  Turbulent friction velocity at surface [m s-1]
            w_s=Ws_1d,                             & ! (out) Turbulent velocity scale profile [m s-1]
            CVMix_kpp_params_user=CS%KPP_params)     !       KPP parameters
          CS%Ws(i,j,:) = Ws_1d(:)
      endif

      ! Diagnostics
      if (CS%id_N2     > 0)   CS%N2(i,j,:)     = N2_1d(:)
      if (CS%id_BulkDrho > 0) CS%dRho(i,j,:)   = deltaRho(:)
      if (CS%id_BulkRi > 0)   CS%BulkRi(i,j,:) = BulkRi_1d(:)
      if (CS%id_BulkUz2 > 0)  CS%Uz2(i,j,:)    = deltaU2(:)
      if (CS%id_Tsurf  > 0)   CS%Tsurf(i,j)    = surfTemp
      if (CS%id_Ssurf  > 0)   CS%Ssurf(i,j)    = surfSalt
      if (CS%id_Usurf  > 0)   CS%Usurf(i,j)    = surfU
      if (CS%id_Vsurf  > 0)   CS%Vsurf(i,j)    = surfv

    enddo
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

  ! BLD smoothing:
  if (CS%n_smooth > 0) call KPP_smooth_BLD(CS,G,GV,h)

end subroutine KPP_compute_BLD


!> Apply a 1-1-4-1-1 Laplacian filter one time on BLD to reduce any horizontal two-grid-point noise
subroutine KPP_smooth_BLD(CS,G,GV,h)
  ! Arguments
  type(KPP_CS),                           pointer       :: CS   !< Control structure
  type(ocean_grid_type),                  intent(inout) :: G    !< Ocean grid
  type(verticalGrid_type),                intent(in)    :: GV   !< Ocean vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h    !< Layer/level thicknesses [H ~> m or kg m-2]

  ! local
  real, dimension(SZI_(G),SZJ_(G)) :: OBLdepth_prev     ! OBLdepth before s.th smoothing iteration
  real, dimension( GV%ke )         :: cellHeight        ! Cell center heights referenced to surface [m]
                                                        ! (negative in the ocean)
  real, dimension( GV%ke+1 )       :: iFaceHeight       ! Interface heights referenced to surface [m]
                                                        ! (negative in the ocean)
  real :: wc, ww, we, wn, ws ! averaging weights for smoothing
  real :: dh                 ! The local thickness used for calculating interface positions [m]
  real :: hcorr              ! A cumulative correction arising from inflation of vanished layers [m]
  integer :: i, j, k, s

  call cpu_clock_begin(id_clock_KPP_smoothing)

  ! Update halos
  call pass_var(CS%OBLdepth, G%Domain, halo=CS%n_smooth)

  if (CS%id_OBLdepth_original > 0) CS%OBLdepth_original = CS%OBLdepth

  do s=1,CS%n_smooth

    OBLdepth_prev = CS%OBLdepth

    ! apply smoothing on OBL depth
    !$OMP parallel do default(none) shared(G, GV, CS, h, OBLdepth_prev) &
    !$OMP                           private(wc, ww, we, wn, ws, dh, hcorr, cellHeight, iFaceHeight)
    do j = G%jsc, G%jec
      do i = G%isc, G%iec

         ! skip land points
        if (G%mask2dT(i,j)==0.) cycle

        iFaceHeight(1) = 0.0 ! BBL is all relative to the surface
        hcorr = 0.
        do k=1,GV%ke

          ! cell center and cell bottom in meters (negative values in the ocean)
          dh = h(i,j,k) * GV%H_to_m ! Nominal thickness to use for increment
          dh = dh + hcorr ! Take away the accumulated error (could temporarily make dh<0)
          hcorr = min( dh - CS%min_thickness, 0. ) ! If inflating then hcorr<0
          dh = max( dh, CS%min_thickness ) ! Limit increment dh>=min_thickness
          cellHeight(k)    = iFaceHeight(k) - 0.5 * dh
          iFaceHeight(k+1) = iFaceHeight(k) - dh
        enddo

        ! compute weights
        ww = 0.125 * G%mask2dT(i-1,j)
        we = 0.125 * G%mask2dT(i+1,j)
        ws = 0.125 * G%mask2dT(i,j-1)
        wn = 0.125 * G%mask2dT(i,j+1)
        wc = 1.0 - (ww+we+wn+ws)

        CS%OBLdepth(i,j) =  wc * OBLdepth_prev(i,j)   &
                          + ww * OBLdepth_prev(i-1,j) &
                          + we * OBLdepth_prev(i+1,j) &
                          + ws * OBLdepth_prev(i,j-1) &
                          + wn * OBLdepth_prev(i,j+1)

        ! Apply OBLdepth smoothing at a cell only if the OBLdepth gets deeper via smoothing.
        if (CS%deepen_only) CS%OBLdepth(i,j) = max(CS%OBLdepth(i,j), OBLdepth_prev(i,j))

        ! prevent OBL depths deeper than the bathymetric depth
        CS%OBLdepth(i,j) = min( CS%OBLdepth(i,j), -iFaceHeight(GV%ke+1) ) ! no deeper than bottom
        CS%kOBL(i,j)     = CVMix_kpp_compute_kOBL_depth( iFaceHeight, cellHeight, CS%OBLdepth(i,j) )
      enddo
    enddo

  enddo ! s-loop

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
                                                       !! to the desired units for BLD
  ! Local variables
  real :: scale  ! A dimensional rescaling factor
  integer :: i,j

  scale = US%m_to_Z ; if (present(m_to_BLD_units)) scale = m_to_BLD_units

  !$OMP parallel do default(none) shared(BLD, CS, G, scale)
  do j = G%jsc, G%jec ; do i = G%isc, G%iec
    BLD(i,j) = scale * CS%OBLdepth(i,j)
  enddo ; enddo

end subroutine KPP_get_BLD

!> Apply KPP non-local transport of surface fluxes for temperature.
subroutine KPP_NonLocalTransport_temp(CS, G, GV, h, nonLocalTrans, surfFlux, &
                                      dt, scalar, C_p)

  type(KPP_CS),                               intent(in)    :: CS     !< Control structure
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h      !< Layer/level thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)   :: nonLocalTrans !< Non-local transport [nondim]
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: surfFlux  !< Surface flux of scalar
                                                                      !! [conc H s-1 ~> conc m s-1 or conc kg m-2 s-1]
  real,                                       intent(in)    :: dt     !< Time-step [s]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: scalar !< temperature
  real,                                       intent(in)    :: C_p    !< Seawater specific heat capacity [J kg-1 degC-1]

  integer :: i, j, k
  real, dimension( SZI_(G), SZJ_(G),SZK_(GV) ) :: dtracer


  dtracer(:,:,:) = 0.0
  !$OMP parallel do default(none) shared(dtracer, nonLocalTrans, h, G, GV, surfFlux)
  do k = 1, GV%ke
    do j = G%jsc, G%jec
      do i = G%isc, G%iec
        dtracer(i,j,k) = ( nonLocalTrans(i,j,k) - nonLocalTrans(i,j,k+1) ) / &
                         ( h(i,j,k) + GV%H_subroundoff ) * surfFlux(i,j)
      enddo
    enddo
  enddo

  !  Update tracer due to non-local redistribution of surface flux
  if (CS%applyNonLocalTrans) then
    !$OMP parallel do default(none) shared(dt, scalar, dtracer, G, GV)
    do k = 1, GV%ke
      do j = G%jsc, G%jec
        do i = G%isc, G%iec
          scalar(i,j,k) = scalar(i,j,k) + dt * dtracer(i,j,k)
        enddo
      enddo
    enddo
  endif

  ! Diagnostics
  if (CS%id_QminusSW        > 0) call post_data(CS%id_QminusSW, surfFlux, CS%diag)
  if (CS%id_NLT_dTdt        > 0) call post_data(CS%id_NLT_dTdt, dtracer,  CS%diag)
  if (CS%id_NLT_temp_budget > 0) then
    dtracer(:,:,:) = 0.0
    !$OMP parallel do default(none) shared(dtracer, nonLocalTrans, surfFlux, C_p, G, GV)
    do k = 1, GV%ke
      do j = G%jsc, G%jec
        do i = G%isc, G%iec
          dtracer(i,j,k) = (nonLocalTrans(i,j,k) - nonLocalTrans(i,j,k+1)) * &
                           surfFlux(i,j) * C_p * GV%H_to_kg_m2
        enddo
      enddo
    enddo
    call post_data(CS%id_NLT_temp_budget, dtracer, CS%diag)
  endif

end subroutine KPP_NonLocalTransport_temp


!> Apply KPP non-local transport of surface fluxes for salinity.
!> This routine is a useful prototype for other material tracers.
subroutine KPP_NonLocalTransport_saln(CS, G, GV, h, nonLocalTrans, surfFlux, dt, scalar)

  type(KPP_CS),                               intent(in)    :: CS            !< Control structure
  type(ocean_grid_type),                      intent(in)    :: G             !< Ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV            !< Ocean vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h             !< Layer/level thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)   :: nonLocalTrans !< Non-local transport [nondim]
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: surfFlux      !< Surface flux of scalar
                                                                        !! [conc H s-1 ~> conc m s-1 or conc kg m-2 s-1]
  real,                                       intent(in)    :: dt            !< Time-step [s]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: scalar        !< Scalar (scalar units [conc])

  integer :: i, j, k
  real, dimension( SZI_(G), SZJ_(G),SZK_(GV) ) :: dtracer


  dtracer(:,:,:) = 0.0
  !$OMP parallel do default(none) shared(dtracer, nonLocalTrans, h, G, GV, surfFlux)
  do k = 1, GV%ke
    do j = G%jsc, G%jec
      do i = G%isc, G%iec
        dtracer(i,j,k) = ( nonLocalTrans(i,j,k) - nonLocalTrans(i,j,k+1) ) / &
                         ( h(i,j,k) + GV%H_subroundoff ) * surfFlux(i,j)
      enddo
    enddo
  enddo

  !  Update tracer due to non-local redistribution of surface flux
  if (CS%applyNonLocalTrans) then
    !$OMP parallel do default(none) shared(G, GV, dt, scalar, dtracer)
    do k = 1, GV%ke
      do j = G%jsc, G%jec
        do i = G%isc, G%iec
          scalar(i,j,k) = scalar(i,j,k) + dt * dtracer(i,j,k)
        enddo
      enddo
    enddo
  endif

  ! Diagnostics
  if (CS%id_netS            > 0) call post_data(CS%id_netS,     surfFlux, CS%diag)
  if (CS%id_NLT_dSdt        > 0) call post_data(CS%id_NLT_dSdt, dtracer,  CS%diag)
  if (CS%id_NLT_saln_budget > 0) then
    dtracer(:,:,:) = 0.0
    !$OMP parallel do default(none) shared(G, GV, dtracer, nonLocalTrans, surfFlux)
    do k = 1, GV%ke
      do j = G%jsc, G%jec
        do i = G%isc, G%iec
          dtracer(i,j,k) = (nonLocalTrans(i,j,k) - nonLocalTrans(i,j,k+1)) * &
                           surfFlux(i,j) * GV%H_to_kg_m2
        enddo
      enddo
    enddo
    call post_data(CS%id_NLT_saln_budget, dtracer, CS%diag)
  endif

end subroutine KPP_NonLocalTransport_saln


!> Clear pointers, deallocate memory
subroutine KPP_end(CS)
  type(KPP_CS), pointer :: CS !< Control structure

  if (.not.associated(CS)) return

  deallocate(CS)

end subroutine KPP_end

end module MOM_CVMix_KPP
