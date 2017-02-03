!> Provides the K-Profile Parameterization (KPP) of Large et al., 1994, via CVMix.
module MOM_KPP

! License goes here?

use MOM_coms,          only : max_across_PEs
use MOM_debugging,     only : hchksum, is_NaN
use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_ptr, post_data
use MOM_diag_mediator, only : query_averaging_enabled, register_diag_field
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_PE
use MOM_EOS,           only : EOS_type, calculate_density
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_file_parser,   only : openParameterBlock, closeParameterBlock
use MOM_grid,          only : ocean_grid_type, isPointInCell
use MOM_verticalGrid,  only : verticalGrid_type

use CVmix_kpp, only : CVmix_init_kpp, CVmix_put_kpp, CVmix_get_kpp_real
use CVmix_kpp, only : CVmix_coeffs_kpp
use CVmix_kpp, only : CVmix_kpp_compute_OBL_depth
use CVmix_kpp, only : CVmix_kpp_compute_turbulent_scales
use CVmix_kpp, only : CVmix_kpp_compute_bulk_Richardson
use CVmix_kpp, only : CVmix_kpp_compute_unresolved_shear
use CVmix_kpp, only : CVmix_kpp_params_type
use CVmix_kpp, only : CVmix_kpp_compute_kOBL_depth

implicit none ; private

#include "MOM_memory.h"

public :: KPP_init
public :: KPP_calculate
public :: KPP_end
public :: KPP_NonLocalTransport_temp
public :: KPP_NonLocalTransport_saln

! Enumerated constants
integer, private, parameter :: NLT_SHAPE_CVMIX     = 0 !< Use the CVmix profile
integer, private, parameter :: NLT_SHAPE_LINEAR    = 1 !< Linear, \f$ G(\sigma) = 1-\sigma \f$
integer, private, parameter :: NLT_SHAPE_PARABOLIC = 2 !< Parabolic, \f$ G(\sigma) = (1-\sigma)^2 \f$
integer, private, parameter :: NLT_SHAPE_CUBIC     = 3 !< Cubic, \f$ G(\sigma) = 1 + (2\sigma-3) \sigma^2\f$
integer, private, parameter :: NLT_SHAPE_CUBIC_LMD = 4 !< Original shape, \f$ G(\sigma) = \frac{27}{4} \sigma (1-\sigma)^2 \f$

integer, private, parameter :: SW_METHOD_ALL_SW = 0 !< Use all shortwave radiation
integer, private, parameter :: SW_METHOD_MXL_SW = 1 !< Use shortwave radiation absorbed in mixing layer
integer, private, parameter :: SW_METHOD_LV1_SW = 2 !< Use shortwave radiation absorbed in layer 1

!> Control structure for containing KPP parameters/data
type, public :: KPP_CS ; private

  ! Parameters
  real    :: Ri_crit                   !< Critical bulk Richardson number (defines OBL depth)
  real    :: vonKarman                 !< von Karman constant (dimensionless)
  real    :: cs                        !< Parameter for computing velocity scale function (dimensionless)
  real    :: cs2                       !< Parameter for multiplying by non-local term
                                       !   This is active for NLT_SHAPE_CUBIC_LMD only
  logical :: enhance_diffusion         !< If True, add enhanced diffusivity at base of boundary layer.
  character(len=10) :: interpType      !< Type of interpolation in determining OBL depth
  logical :: computeEkman              !< If True, compute Ekman depth limit for OBLdepth
  logical :: computeMoninObukhov       !< If True, compute Monin-Obukhov limit for OBLdepth
  logical :: passiveMode               !< If True, makes KPP passive meaning it does NOT alter the diffusivity
  real    :: deepOBLoffset             !< If non-zero, is a distance from the bottom that the OBL can not penetrate through (m)
  real    :: minOBLdepth               !< If non-zero, is a minimum depth for the OBL (m)
  real    :: surf_layer_ext            !< Fraction of OBL depth considered in the surface layer (nondim)
  real    :: minVtsqr                  !< Min for the squared unresolved velocity used in Rib CVMix calculation (m2/s2)
  logical :: fixedOBLdepth             !< If True, will fix the OBL depth at fixedOBLdepth_value
  real    :: fixedOBLdepth_value       !< value for the fixed OBL depth when fixedOBLdepth==True.
  logical :: debug                     !< If True, calculate checksums and write debugging information
  character(len=30) :: MatchTechnique  !< Method used in CVMix for setting diffusivity and NLT profile functions
  integer :: NLT_shape                 !< MOM6 over-ride of CVMix NLT shape function
  logical :: applyNonLocalTrans        !< If True, apply non-local transport to heat and scalars
  logical :: KPPzeroDiffusivity        !< If True, will set diffusivity and viscosity from KPP to zero; for testing purposes.
  logical :: KPPisAdditive             !< If True, will add KPP diffusivity to initial diffusivity.
                                       !! If False, will replace initial diffusivity wherever KPP diffusivity is non-zero.
  real    :: min_thickness             !< A minimum thickness used to avoid division by small numbers in the vicinity of vanished layers.
  ! smg: obsolete below
  logical :: correctSurfLayerAvg       !< If true, applies a correction to the averaging of surface layer properties
  real    :: surfLayerDepth            !< A guess at the depth of the surface layer (which should 0.1 of OBLdepth) (m)
  ! smg: obsolete above
  integer :: SW_METHOD                 !<Sets method for using shortwave radiation in surface buoyancy flux

  !> CVmix parameters
  type(CVmix_kpp_params_type), pointer :: KPP_params => NULL()

  ! Diagnostic handles and pointers
  type(diag_ctrl), pointer :: diag => NULL()
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

  ! Diagnostics arrays
  real, allocatable, dimension(:,:)   :: OBLdepth  !< Depth (positive) of OBL (m)
  real, allocatable, dimension(:,:,:) :: dRho      !< Bulk difference in density (kg/m3)
  real, allocatable, dimension(:,:,:) :: Uz2       !< Square of bulk difference in resolved velocity (m2/s2)
  real, allocatable, dimension(:,:,:) :: BulkRi    !< Bulk Richardson number for each layer (dimensionless)
  real, allocatable, dimension(:,:,:) :: sigma     !< Sigma coordinate (dimensionless)
  real, allocatable, dimension(:,:,:) :: Ws        !< Turbulent velocity scale for scalars (m/s)
  real, allocatable, dimension(:,:,:) :: N         !< Brunt-Vaisala frequency (1/s)
  real, allocatable, dimension(:,:,:) :: N2        !< Squared Brunt-Vaisala frequency (1/s2)
  real, allocatable, dimension(:,:,:) :: Vt2       !< Unresolved squared turbulence velocity for bulk Ri (m2/s2)
  real, allocatable, dimension(:,:,:) :: Kt_KPP    !< Temp diffusivity from KPP (m2/s)
  real, allocatable, dimension(:,:,:) :: Ks_KPP    !< Scalar diffusivity from KPP (m2/s)
  real, allocatable, dimension(:,:,:) :: Kv_KPP    !< Viscosity due to KPP (m2/s)
  real, allocatable, dimension(:,:)   :: Tsurf     !< Temperature of surface layer (C)
  real, allocatable, dimension(:,:)   :: Ssurf     !< Salinity of surface layer (ppt)
  real, allocatable, dimension(:,:)   :: Usurf     !< i-velocity of surface layer (m/s)
  real, allocatable, dimension(:,:)   :: Vsurf     !< j-velocity of surface layer (m/s)

end type KPP_CS

! Module data used for debugging only
logical, parameter :: verbose = .False.
#define __DO_SAFETY_CHECKS__

contains

!> Initialize the CVmix KPP module and set up diagnostics
!! Returns True if KPP is to be used, False otherwise.
logical function KPP_init(paramFile, G, diag, Time, CS, passive)

  ! Arguments
  type(param_file_type),   intent(in)    :: paramFile !< File parser
  type(ocean_grid_type),   intent(in)    :: G         !< Ocean grid
  type(diag_ctrl), target, intent(in)    :: diag      !< Diagnostics
  type(time_type),         intent(in)    :: Time      !< Time
  type(KPP_CS),            pointer       :: CS        !< Control structure
  logical, optional,       intent(out)   :: passive   !< Copy of %passiveMode

  ! Local variables
#include "version_variable.h"
  character(len=40) :: mod = 'MOM_KPP' ! name of this module
  character(len=20) :: string          ! local temporary string
  logical :: CS_IS_ONE=.false.

  if (associated(CS)) call MOM_error(FATAL, 'MOM_KPP, KPP_init: '// &
           'Control structure has already been initialized')
  allocate(CS)

  ! Read parameters
  call log_version(paramFile, mod, version, 'This is the MOM wrapper to CVmix:KPP\n' // &
            'See http://code.google.com/p/cvmix/')
  call get_param(paramFile, mod, "USE_KPP", KPP_init, &
                 "If true, turns on the [CVmix] KPP scheme of Large et al., 1994,\n"// &
                 "to calculate diffusivities and non-local transport in the OBL.",     &
                 default=.false.)
  ! Forego remainder of initialization if not using this scheme
  if (.not. KPP_init) return

  call openParameterBlock(paramFile,'KPP')
  call get_param(paramFile, mod, 'PASSIVE', CS%passiveMode,           &
                 'If True, puts KPP into a passive-diagnostic mode.', &
                  default=.False.)
  if (present(passive)) passive=CS%passiveMode ! This is passed back to the caller so
                                               ! the caller knows to not use KPP output
  call get_param(paramFile, mod, 'APPLY_NONLOCAL_TRANSPORT', CS%applyNonLocalTrans,  &
                 'If True, applies the non-local transport to heat and scalars.\n'//  &
                 'If False, calculates the non-local transport and tendencies but\n'//&
                 'purely for diagnostic purposes.',                                   &
                 default=.not. CS%passiveMode)
  call get_param(paramFile, mod, 'RI_CRIT', CS%Ri_crit,                            &
                 'Critical bulk Richardson number used to define depth of the\n'// &
                 'surface Ocean Boundary Layer (OBL).',                            &
                 units='nondim', default=0.3)
  call get_param(paramFile, mod, 'VON_KARMAN', CS%vonKarman, &
                 'von Karman constant.',                     &
                 units='nondim', default=0.40)
  call get_param(paramFile, mod, 'ENHANCE_DIFFUSION', CS%enhance_diffusion,              &
                 'If True, adds enhanced diffusion at the based of the boundary layer.', &
                 default=.true.)
  call get_param(paramFile, mod, 'INTERP_TYPE', CS%interpType,           &
                 'Type of interpolation to determine the OBL depth.\n'// &
                 'Allowed types are: linear, quadratic, cubic.',         &
                 default='cubic')
  call get_param(paramFile, mod, 'COMPUTE_EKMAN', CS%computeEkman,             &
                 'If True, limit OBL depth to be no deeper than Ekman depth.', &
                 default=.False.)
  call get_param(paramFile, mod, 'COMPUTE_MONIN_OBUKHOV', CS%computeMoninObukhov, &
                 'If True, limit the OBL depth to be no deeper than\n'//          &
                 'Monin-Obukhov depth.',                                          &
                 default=.False.)
  call get_param(paramFile, mod, 'CS', CS%cs,                        &
                 'Parameter for computing velocity scale function.', &
                 units='nondim', default=98.96)
  call get_param(paramFile, mod, 'CS2', CS%cs2,                        &
                 'Parameter for computing non-local term.', &
                 units='nondim', default=6.32739901508)
  call get_param(paramFile, mod, 'DEEP_OBL_OFFSET', CS%deepOBLoffset,                             &
                 'If non-zero, the distance above the bottom to which the OBL is clipped\n'//     &
                 'if it would otherwise reach the bottom. The smaller of this and 0.1D is used.', &
                 units='m',default=0.)
  call get_param(paramFile, mod, 'FIXED_OBLDEPTH', CS%fixedOBLdepth,       &
                 'If True, fix the OBL depth to FIXED_OBLDEPTH_VALUE\n'//  &
                 'rather than using the OBL depth from CVMix.\n'//         &
                 'This option is just for testing purposes.',              &
                 default=.False.)
  call get_param(paramFile, mod, 'FIXED_OBLDEPTH_VALUE', CS%fixedOBLdepth_value,  &
                 'Value for the fixed OBL depth when fixedOBLdepth==True. \n'//   &
                 'This parameter is for just for testing purposes. \n'//          &
                 'It will over-ride the OBLdepth computed from CVMix.',           &
                 units='m',default=30.0)
  call get_param(paramFile, mod, 'SURF_LAYER_EXTENT', CS%surf_layer_ext,   &
                 'Fraction of OBL depth considered in the surface layer.', &
                 units='nondim',default=0.10)
  call get_param(paramFile, mod, 'MINIMUM_OBL_DEPTH', CS%minOBLdepth,                            &
                 'If non-zero, a minimum depth to use for KPP OBL depth. Independent of\n'//     &
                 'this parameter, the OBL depth is always at least as deep as the first layer.', &
                 units='m',default=0.)
  call get_param(paramFile, mod, 'MINIMUM_VT2', CS%minVtsqr,                                   &
                 'Min of the unresolved velocity Vt2 used in Rib CVMix calculation.     \n'//  &
                 'Scaling: MINIMUM_VT2 = const1*d*N*ws, with d=1m, N=1e-5/s, ws=1e-6 m/s.',    &
                 units='m2/s2',default=1e-10)

! smg: for removal below
  call get_param(paramFile, mod, 'CORRECT_SURFACE_LAYER_AVERAGE', CS%correctSurfLayerAvg,   &
                 'If true, applies a correction step to the averaging of surface layer\n'// &
                 'properties. This option is obsolete.', default=.False.)
  call get_param(paramFile, mod, 'FIRST_GUESS_SURFACE_LAYER_DEPTH', CS%surfLayerDepth,              &
                 'The first guess at the depth of the surface layer used for averaging\n'//         &
                 'the surface layer properties. If =0, the top model level properties\n'//          &
                 'will be used for the surface layer. If CORRECT_SURFACE_LAYER_AVERAGE=True, a\n'// &
                 'subsequent correction is applied. This parameter is obsolete', units='m', default=0.)
! smg: for removal above

  call get_param(paramFile, mod, 'NLT_SHAPE', string, &
                 'MOM6 method to set nonlocal transport profile.\n'//                          &
                 'Over-rides the result from CVMix.  Allowed values are: \n'//                 &
                 '\t CVMIX     - Uses the profiles from CVmix specified by MATCH_TECHNIQUE\n'//&
                 '\t LINEAR    - A linear profile, 1-sigma\n'//                                &
                 '\t PARABOLIC - A parablic profile, (1-sigma)^2\n'//                          &
                 '\t CUBIC     - A cubic profile, (1-sigma)^2(1+2*sigma)\n'//                  &
                 '\t CUBIC_LMD - The original KPP profile',                                    &
                 default='CVMIX')
  select case ( trim(string) )
    case ("CVMIX")     ; CS%NLT_shape = NLT_SHAPE_CVMIX
    case ("LINEAR")    ; CS%NLT_shape = NLT_SHAPE_LINEAR
    case ("PARABOLIC") ; CS%NLT_shape = NLT_SHAPE_PARABOLIC
    case ("CUBIC")     ; CS%NLT_shape = NLT_SHAPE_CUBIC
    case ("CUBIC_LMD") ; CS%NLT_shape = NLT_SHAPE_CUBIC_LMD
    case default ; call MOM_error(FATAL,"KPP_init: "// &
                   "Unrecognized NLT_SHAPE option"//trim(string))
  end select
  call get_param(paramFile, mod, 'MATCH_TECHNIQUE', CS%MatchTechnique,                                    &
                 'CVMix method to set profile function for diffusivity and NLT,\n'//                      &
                 'as well as matching across OBL base. Allowed values are: \n'//                          &
                 '\t SimpleShapes      = sigma*(1-sigma)^2 for both diffusivity and NLT\n'//              &
                 '\t MatchGradient     = sigma*(1-sigma)^2 for NLT; diffusivity profile from matching\n'//&
                 '\t MatchBoth         = match gradient for both diffusivity and NLT\n'//                 &
                 '\t ParabolicNonLocal = sigma*(1-sigma)^2 for diffusivity; (1-sigma)^2 for NLT',         &
                 default='SimpleShapes')
  if (CS%MatchTechnique.eq.'ParabolicNonLocal') then
     ! This forces Cs2 (Cs in non-local computation) to equal 1 for parabolic non-local option.
     !  May be used during CVmix initialization.
     Cs_is_one=.true.
  endif
  call get_param(paramFile, mod, 'KPP_ZERO_DIFFUSIVITY', CS%KPPzeroDiffusivity,            &
                 'If True, zeroes the KPP diffusivity and viscosity; for testing purpose.',&
                 default=.False.)
  call get_param(paramFile, mod, 'KPP_IS_ADDITIVE', CS%KPPisAdditive,                &
                 'If true, adds KPP diffusivity to diffusivity from other schemes.'//&
                 'If false, KPP is the only diffusivity wherever KPP is non-zero.',  &
                 default=.True.)
  call get_param(paramFile, mod, 'KPP_SHORTWAVE_METHOD',string,                      &
                 'Determines contribution of shortwave radiation to KPP surface '// &
                 'buoyancy flux.  Options include:\n'//                             &
                 '  ALL_SW: use total shortwave radiation\n'//                      &
                 '  MXL_SW:  use shortwave radiation absorbed by mixing layer\n'//  &
                 '  LV1_SW:  use shortwave radiation absorbed by top model layer',  &
                 default='MXL_SW')
  select case ( trim(string) )
    case ("ALL_SW") ; CS%SW_METHOD = SW_METHOD_ALL_SW
    case ("MXL_SW") ; CS%SW_METHOD = SW_METHOD_MXL_SW
    case ("LV1_SW") ; CS%SW_METHOD = SW_METHOD_LV1_SW
    case default ; call MOM_error(FATAL,"KPP_init: "// &
                   "Unrecognized KPP_SHORTWAVE_METHOD option"//trim(string))
  end select
  call get_param(paramFile, mod, 'CVMIX_ZERO_H_WORK_AROUND', CS%min_thickness,                           &
                 'A minimum thickness used to avoid division by small numbers in the vicinity\n'//       &
                 'of vanished layers. This is independent of MIN_THICKNESS used in other parts of MOM.', &
                 units='m', default=0.)

  call closeParameterBlock(paramFile)
  call get_param(paramFile, mod, 'DEBUG', CS%debug, default=.False., do_not_log=.True.)

  call CVmix_init_kpp( Ri_crit=CS%Ri_crit,                 &
                       minOBLdepth=CS%minOBLdepth,         &
                       minVtsqr=CS%minVtsqr,               &
                       vonKarman=CS%vonKarman,             &
                       surf_layer_ext=CS%surf_layer_ext,   &
                       interp_type=CS%interpType,          &
                       lEkman=CS%computeEkman,             &
                       lMonOb=CS%computeMoninObukhov,      &
                       MatchTechnique=CS%MatchTechnique,   &
                       lenhanced_diff=CS%enhance_diffusion,&
                       lnonzero_surf_nonlocal=Cs_is_one   ,&
                       CVmix_kpp_params_user=CS%KPP_params )

  ! Register diagnostics
  CS%diag => diag
  CS%id_OBLdepth = register_diag_field('ocean_model', 'KPP_OBLdepth', diag%axesT1, Time, &
      'Thickness of the surface Ocean Boundary Layer calculated by [CVmix] KPP', 'meter', &
      cmor_field_name='oml', cmor_long_name='ocean_mixed_layer_thickness_defined_by_mixing_scheme', &
      cmor_units='m', cmor_standard_name='Ocean Mixed Layer Thickness Defined by Mixing Scheme')
      ! CMOR names are placeholders; must be modified by time period
      ! for CMOR compliance. Diag manager will be used for omlmax and
      ! omldamax.
  CS%id_BulkDrho = register_diag_field('ocean_model', 'KPP_BulkDrho', diag%axesTL, Time, &
      'Bulk difference in density used in Bulk Richardson number, as used by [CVmix] KPP', 'kg/m3')
  CS%id_BulkUz2 = register_diag_field('ocean_model', 'KPP_BulkUz2', diag%axesTL, Time, &
      'Square of bulk difference in resolved velocity used in Bulk Richardson number via [CVmix] KPP', 'm2/s2')
  CS%id_BulkRi = register_diag_field('ocean_model', 'KPP_BulkRi', diag%axesTL, Time, &
      'Bulk Richardson number used to find the OBL depth used by [CVmix] KPP', 'nondim')
  CS%id_Sigma = register_diag_field('ocean_model', 'KPP_sigma', diag%axesTi, Time, &
      'Sigma coordinate used by [CVmix] KPP', 'nondim')
  CS%id_Ws = register_diag_field('ocean_model', 'KPP_Ws', diag%axesTL, Time, &
      'Turbulent vertical velocity scale for scalars used by [CVmix] KPP', 'm/s')
  CS%id_N = register_diag_field('ocean_model', 'KPP_N', diag%axesTi, Time, &
      '(Adjusted) Brunt-Vaisala frequency used by [CVmix] KPP', '1/s')
  CS%id_N2 = register_diag_field('ocean_model', 'KPP_N2', diag%axesTi, Time, &
      'Square of Brunt-Vaisala frequency used by [CVmix] KPP', '1/s2')
  CS%id_Vt2 = register_diag_field('ocean_model', 'KPP_Vt2', diag%axesTL, Time, &
      'Unresolved shear turbulence used by [CVmix] KPP', 'm2/s2')
  CS%id_uStar = register_diag_field('ocean_model', 'KPP_uStar', diag%axesT1, Time, &
      'Friction velocity, u*, as used by [CVmix] KPP', 'm/s')
  CS%id_buoyFlux = register_diag_field('ocean_model', 'KPP_buoyFlux', diag%axesTi, Time, &
      'Surface (and penetrating) buoyancy flux, as used by [CVmix] KPP', 'm2/s3')
  CS%id_QminusSW = register_diag_field('ocean_model', 'KPP_QminusSW', diag%axesT1, Time, &
      'Net temperature flux ignoring short-wave, as used by [CVmix] KPP', 'K m/s')
  CS%id_netS = register_diag_field('ocean_model', 'KPP_netSalt', diag%axesT1, Time, &
      'Effective net surface salt flux, as used by [CVmix] KPP', 'ppt m/s')
  CS%id_Kt_KPP = register_diag_field('ocean_model', 'KPP_Kheat', diag%axesTi, Time, &
      'Heat diffusivity due to KPP, as calculated by [CVmix] KPP', 'm2/s')
  CS%id_Kd_in = register_diag_field('ocean_model', 'KPP_Kd_in', diag%axesTi, Time, &
      'Diffusivity passed to KPP', 'm2/s')
  CS%id_Ks_KPP = register_diag_field('ocean_model', 'KPP_Ksalt', diag%axesTi, Time, &
      'Salt diffusivity due to KPP, as calculated by [CVmix] KPP', 'm2/s')
  CS%id_Kv_KPP = register_diag_field('ocean_model', 'KPP_Kv', diag%axesTi, Time, &
      'Vertical viscosity due to KPP, as calculated by [CVmix] KPP', 'm2/s')
  CS%id_NLTt = register_diag_field('ocean_model', 'KPP_NLtransport_heat', diag%axesTi, Time, &
      'Non-local transport (Cs*G(sigma)) for heat, as calculated by [CVmix] KPP', 'nondim')
  CS%id_NLTs = register_diag_field('ocean_model', 'KPP_NLtransport_salt', diag%axesTi, Time, &
      'Non-local tranpsort (Cs*G(sigma)) for scalars, as calculated by [CVmix] KPP', 'nondim')
  CS%id_NLT_dTdt = register_diag_field('ocean_model', 'KPP_NLT_dTdt', diag%axesTL, Time, &
      'Temperature tendency due to non-local transport of heat, as calculated by [CVmix] KPP', 'K/s')
  CS%id_NLT_dSdt = register_diag_field('ocean_model', 'KPP_NLT_dSdt', diag%axesTL, Time, &
      'Salinity tendency due to non-local transport of salt, as calculated by [CVmix] KPP', 'ppt/s')
  CS%id_NLT_temp_budget = register_diag_field('ocean_model', 'KPP_NLT_temp_budget', diag%axesTL, Time, &
      'Heat content change due to non-local transport, as calculated by [CVmix] KPP', 'W/m^2')
  CS%id_NLT_saln_budget = register_diag_field('ocean_model', 'KPP_NLT_saln_budget', diag%axesTL, Time, &
      'Salt content change due to non-local transport, as calculated by [CVmix] KPP', 'kg/(sec*m^2)')
  CS%id_Tsurf = register_diag_field('ocean_model', 'KPP_Tsurf', diag%axesT1, Time, &
      'Temperature of surface layer (10% of OBL depth) as passed to [CVmix] KPP', 'C')
  CS%id_Ssurf = register_diag_field('ocean_model', 'KPP_Ssurf', diag%axesT1, Time, &
      'Salinity of surface layer (10% of OBL depth) as passed to [CVmix] KPP', 'ppt')
  CS%id_Usurf = register_diag_field('ocean_model', 'KPP_Usurf', diag%axesCu1, Time, &
      'i-component flow of surface layer (10% of OBL depth) as passed to [CVmix] KPP', 'm/s')
  CS%id_Vsurf = register_diag_field('ocean_model', 'KPP_Vsurf', diag%axesCv1, Time, &
      'j-component flow of surface layer (10% of OBL depth) as passed to [CVmix] KPP', 'm/s')

  if (CS%id_OBLdepth > 0) allocate( CS%OBLdepth( SZI_(G), SZJ_(G) ) )
  if (CS%id_OBLdepth > 0) CS%OBLdepth(:,:) = 0.
  if (CS%id_BulkDrho > 0) allocate( CS%dRho( SZI_(G), SZJ_(G), SZK_(G) ) )
  if (CS%id_BulkDrho > 0) CS%dRho(:,:,:) = 0.
  if (CS%id_BulkUz2 > 0)  allocate( CS%Uz2( SZI_(G), SZJ_(G), SZK_(G) ) )
  if (CS%id_BulkUz2 > 0)  CS%Uz2(:,:,:) = 0.
  if (CS%id_BulkRi > 0)   allocate( CS%BulkRi( SZI_(G), SZJ_(G), SZK_(G) ) )
  if (CS%id_BulkRi > 0)   CS%BulkRi(:,:,:) = 0.
  if (CS%id_Sigma > 0)    allocate( CS%sigma( SZI_(G), SZJ_(G), SZK_(G)+1 ) )
  if (CS%id_Sigma > 0)    CS%sigma(:,:,:) = 0.
  if (CS%id_Ws > 0)       allocate( CS%Ws( SZI_(G), SZJ_(G), SZK_(G) ) )
  if (CS%id_Ws > 0)       CS%Ws(:,:,:) = 0.
  if (CS%id_N > 0)        allocate( CS%N( SZI_(G), SZJ_(G), SZK_(G)+1 ) )
  if (CS%id_N > 0)        CS%N(:,:,:) = 0.
  if (CS%id_N2 > 0)       allocate( CS%N2( SZI_(G), SZJ_(G), SZK_(G)+1 ) )
  if (CS%id_N2 > 0)       CS%N2(:,:,:) = 0.
  if (CS%id_Vt2 > 0)      allocate( CS%Vt2( SZI_(G), SZJ_(G), SZK_(G) ) )
  if (CS%id_Vt2 > 0)      CS%Vt2(:,:,:) = 0.
  if (CS%id_Kt_KPP > 0)   allocate( CS%Kt_KPP( SZI_(G), SZJ_(G), SZK_(G)+1 ) )
  if (CS%id_Kt_KPP > 0)   CS%Kt_KPP(:,:,:) = 0.
  if (CS%id_Ks_KPP > 0)   allocate( CS%Ks_KPP( SZI_(G), SZJ_(G), SZK_(G)+1 ) )
  if (CS%id_Ks_KPP > 0)   CS%Ks_KPP(:,:,:) = 0.
  if (CS%id_Kv_KPP > 0)   allocate( CS%Kv_KPP( SZI_(G), SZJ_(G), SZK_(G)+1 ) )
  if (CS%id_Kv_KPP > 0)   CS%Kv_KPP(:,:,:) = 0.
  if (CS%id_Tsurf > 0)    allocate( CS%Tsurf( SZI_(G), SZJ_(G)) )
  if (CS%id_Tsurf > 0)    CS%Tsurf(:,:) = 0.
  if (CS%id_Ssurf > 0)    allocate( CS%Ssurf( SZI_(G), SZJ_(G)) )
  if (CS%id_Ssurf > 0)    CS%Ssurf(:,:) = 0.
  if (CS%id_Usurf > 0)    allocate( CS%Usurf( SZIB_(G), SZJ_(G)) )
  if (CS%id_Usurf > 0)    CS%Usurf(:,:) = 0.
  if (CS%id_Vsurf > 0)    allocate( CS%Vsurf( SZI_(G), SZJB_(G)) )
  if (CS%id_Vsurf > 0)    CS%Vsurf(:,:) = 0.

end function KPP_init



!> KPP vertical diffusivity/viscosity and non-local tracer transport
subroutine KPP_calculate(CS, G, GV, h, Temp, Salt, u, v, EOS, uStar, &
                         buoyFlux, Kt, Ks, Kv, nonLocalTransHeat,&
                         nonLocalTransScalar)

  ! Arguments
  type(KPP_CS),                           pointer       :: CS             !< Control structure
  type(ocean_grid_type),                  intent(in)    :: G              !< Ocean grid
  type(verticalGrid_type),                intent(in)    :: GV             !< Ocean vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h              !< Layer/level thicknesses (units of H)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: Temp           !< potential/cons temp (deg C)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: Salt           !< Salinity (ppt)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u              !< Velocity i-component (m/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v              !< Velocity j-component (m/s)
  type(EOS_type),                         pointer       :: EOS            !< Equation of state
  real, dimension(SZI_(G),SZJ_(G)),         intent(in)    :: uStar          !< Surface friction velocity (m/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(in)    :: buoyFlux !< Surface buoyancy flux (m2/s3)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Kt       !< (in)  Vertical diffusivity of heat w/o KPP (m2/s)
                                                                          !< (out) Vertical diffusivity including KPP (m2/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Ks       !< (in)  Vertical diffusivity of salt w/o KPP (m2/s)
                                                                          !< (out) Vertical diffusivity including KPP (m2/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Kv       !< (in)  Vertical viscosity w/o KPP (m2/s)
                                                                          !< (out) Vertical viscosity including KPP (m2/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: nonLocalTransHeat   !< Temp non-local transport (m/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: nonLocalTransScalar !< scalar non-local transport (m/s)

  ! Local variables
  integer :: i, j, k, km1,kp1                    ! Loop indices
  real, dimension( G%ke )     :: cellHeight      ! Cell center heights referenced to surface (m) (negative in ocean)
  real, dimension( G%ke+1 )   :: iFaceHeight     ! Interface heights referenced to surface (m) (negative in ocean)
  real, dimension( G%ke+1 )   :: N2_1d           ! Brunt-Vaisala frequency squared, at interfaces (1/s2)
  real, dimension( G%ke+1 )   :: N_1d            ! Brunt-Vaisala frequency at interfaces (1/s) (floored at 0)
  real, dimension( G%ke )     :: Ws_1d           ! Profile of vertical velocity scale for scalars (m/s)
  real, dimension( G%ke )     :: Wm_1d           ! Profile of vertical velocity scale for momentum (m/s)
  real, dimension( G%ke )     :: Vt2_1d          ! Unresolved velocity for bulk Ri calculation/diagnostic (m2/s2)
  real, dimension( G%ke )     :: BulkRi_1d       ! Bulk Richardson number for each layer
  real, dimension( G%ke )     :: deltaRho        ! delta Rho in numerator of Bulk Ri number
  real, dimension( G%ke )     :: deltaU2         ! square of delta U (shear) in denominator of Bulk Ri (m2/s2)
  real, dimension( G%ke+1, 2) :: Kdiffusivity    ! Vertical diffusivity at interfaces (m2/s)
  real, dimension( G%ke+1 )   :: Kviscosity      ! Vertical viscosity at interfaces (m2/s)
  real, dimension( G%ke+1, 2) :: nonLocalTrans   ! Non-local transport for heat/salt at interfaces (non-dimensional)
  real, dimension( G%ke )     :: surfBuoyFlux2

  ! for EOS calculation
  real, dimension( 3*G%ke )   :: rho_1D
  real, dimension( 3*G%ke )   :: pres_1D
  real, dimension( 3*G%ke )   :: Temp_1D
  real, dimension( 3*G%ke )   :: Salt_1D

  real :: kOBL, OBLdepth_0d, surfFricVel, surfBuoyFlux, Coriolis
  real :: GoRho, pRef, rho1, rhoK, rhoKm1, Uk, Vk, sigma

  real :: zBottomMinusOffset   ! Height of bottom plus a little bit (m)
  real :: SLdepth_0d           ! Surface layer depth = surf_layer_ext*OBLdepth.
  real :: hTot                 ! Running sum of thickness used in the surface layer average (m)
  real :: delH                 ! Thickness of a layer (m)
  real :: surfHtemp, surfTemp  ! Integral and average of temp over the surface layer
  real :: surfHsalt, surfSalt  ! Integral and average of saln over the surface layer
  real :: surfHu, surfU        ! Integral and average of u over the surface layer
  real :: surfHv, surfV        ! Integral and average of v over the surface layer
  real :: dh    ! The local thickness used for calculating interface positions (m)
  real :: hcorr ! A cumulative correction arising from inflation of vanished layers (m)
  integer :: kk, ksfc, ktmp

#ifdef __DO_SAFETY_CHECKS__
  if (CS%debug) then
    call hchksum(h*GV%H_to_m, "KPP in: h",G%HI,haloshift=0)
    call hchksum(Temp, "KPP in: T",G%HI,haloshift=0)
    call hchksum(Salt, "KPP in: S",G%HI,haloshift=0)
    call hchksum(u, "KPP in: u",G%HI,haloshift=0)
    call hchksum(v, "KPP in: v",G%HI,haloshift=0)
    call hchksum(uStar, "KPP in: uStar",G%HI,haloshift=0)
    call hchksum(buoyFlux, "KPP in: buoyFlux",G%HI,haloshift=0)
    call hchksum(Kt, "KPP in: Kt",G%HI,haloshift=0)
    call hchksum(Ks, "KPP in: Ks",G%HI,haloshift=0)
  endif
#endif

  ! some constants
  GoRho = GV%g_Earth / GV%Rho0
  nonLocalTrans(:,:) = 0.0

  if (CS%id_Kd_in > 0) call post_data(CS%id_Kd_in, Kt, CS%diag)

!$OMP parallel do default(none) shared(G,GV,CS,EOS,uStar,Temp,Salt,u,v,h,GoRho,       &
!$OMP                                  buoyFlux, nonLocalTransHeat,                   &
!$OMP                                  nonLocalTransScalar,Kt,Ks,Kv)                  &
!$OMP                     firstprivate(nonLocalTrans)                                 &
!$OMP                          private(Coriolis,surfFricVel,SLdepth_0d,hTot,surfTemp, &
!$OMP                                  surfHtemp,surfSalt,surfHsalt,surfU,            &
!$OMP                                  surfHu,surfV,surfHv,iFaceHeight,               &
!$OMP                                  pRef,km1,cellHeight,Uk,Vk,deltaU2,             &
!$OMP                                  rho1,rhoK,rhoKm1,deltaRho,N2_1d,N_1d,delH,     &
!$OMP                                  surfBuoyFlux,Ws_1d,Vt2_1d,BulkRi_1d,           &
!$OMP                                  OBLdepth_0d,zBottomMinusOffset,Kdiffusivity,   &
!$OMP                                  Kviscosity,sigma,kOBL,kk,pres_1D,Temp_1D,      &
!$OMP                                  Salt_1D,rho_1D,surfBuoyFlux2,ksfc,dh,hcorr)

  ! loop over horizontal points on processor
  do j = G%jsc, G%jec
    do i = G%isc, G%iec

      ! skip calling KPP for land points
      if (G%mask2dT(i,j)==0.) cycle

      ! things independent of position within the column
      Coriolis = 0.25*( (G%CoriolisBu(i,j)   + G%CoriolisBu(i-1,j-1)) &
                       +(G%CoriolisBu(i-1,j) + G%CoriolisBu(i,j-1)) )
      surfFricVel = uStar(i,j)

      ! Bullk Richardson number computed for each cell in a column,
      ! assuming OBLdepth = grid cell depth. After Rib(k) is
      ! known for the column, then CVMix interpolates to find
      ! the actual OBLdepth. This approach avoids need to iterate
      ! on the OBLdepth calculation. It follows that used in MOM5
      ! and POP.
      iFaceHeight(1) = 0.0 ! BBL is all relative to the surface
      pRef = 0.
      hcorr = 0.
      do k=1,G%ke

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
        hTot     =0.0
        do ktmp = 1,ksfc

          ! SLdepth_0d can be between cell interfaces
          delH = min( max(0.0, SLdepth_0d - hTot), h(i,j,ktmp)*GV%H_to_m )

          ! surface layer thickness
          hTot = hTot + delH

          ! surface averaged fields
          surfHtemp = surfHtemp + Temp(i,j,ktmp) * delH
          surfHsalt = surfHsalt + Salt(i,j,ktmp) * delH
          surfHu    = surfHu + 0.5*(u(i,j,ktmp)+u(i-1,j,ktmp)) * delH
          surfHv    = surfHv + 0.5*(v(i,j,ktmp)+v(i,j-1,ktmp)) * delH

        enddo
        surfTemp = surfHtemp / hTot
        surfSalt = surfHsalt / hTot
        surfU    = surfHu    / hTot
        surfV    = surfHv    / hTot

        ! vertical shear between present layer and
        ! surface layer averaged surfU,surfV.
        ! C-grid average to get Uk and Vk on T-points.
        Uk         = 0.5*(u(i,j,k)+u(i-1,j,k)) - surfU
        Vk         = 0.5*(v(i,j,k)+v(i,j-1,k)) - surfV
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

        ! pRef is pressure at interface between k and km1.
        ! iterate pRef for next pass through k-loop.
        pRef = pRef + GV%H_to_Pa * h(i,j,k)

        ! this difference accounts for penetrating SW
        surfBuoyFlux2(k) = buoyFlux(i,j,1) - buoyFlux(i,j,k+1)

      enddo ! k-loop finishes

      ! compute in-situ density
      call calculate_density(Temp_1D, Salt_1D, pres_1D, rho_1D, 1, 3*G%ke, EOS)

      ! N2 (can be negative) and N (non-negative) on interfaces.
      ! deltaRho is non-local rho difference used for bulk Richardson number.
      ! N_1d is local N (with floor) used for unresolved shear calculation.
      do k = 1, G%ke
        km1 = max(1, k-1)
        kk = 3*(k-1)
        deltaRho(k) = rho_1D(kk+2) - rho_1D(kk+1)
        N2_1d(k)    = (GoRho * (rho_1D(kk+2) - rho_1D(kk+3)) ) / &
                      ((0.5*(h(i,j,km1) + h(i,j,k))+GV%H_subroundoff)*GV%H_to_m)
        N_1d(k)     = sqrt( max( N2_1d(k), 0.) )
      enddo
      N2_1d(G%ke+1 ) = 0.0
      N_1d(G%ke+1 )  = 0.0

      ! turbulent velocity scales w_s and w_m computed at the cell centers.
      ! Note that if sigma > CS%surf_layer_ext, then CVmix_kpp_compute_turbulent_scales
      ! computes w_s and w_m velocity scale at sigma=CS%surf_layer_ext. So we only pass
      ! sigma=CS%surf_layer_ext for this calculation.
      call CVmix_kpp_compute_turbulent_scales( &
        CS%surf_layer_ext, & ! (in)  Normalized surface layer depth; sigma = CS%surf_layer_ext
        -cellHeight,       & ! (in)  Assume here that OBL depth (m) = -cellHeight(k)
        surfBuoyFlux2,     & ! (in)  Buoyancy flux at surface (m2/s3)
        surfFricVel,       & ! (in)  Turbulent friction velocity at surface (m/s)
        w_s=Ws_1d,         & ! (out) Turbulent velocity scale profile (m/s)
        CVmix_kpp_params_user=CS%KPP_params )

      ! Calculate Bulk Richardson number from eq (21) of LMD94
      BulkRi_1d = CVmix_kpp_compute_bulk_Richardson( &
                  cellHeight(1:G%ke),                & ! Depth of cell center (m)
                  GoRho*deltaRho,                    & ! Bulk buoyancy difference, Br-B(z) (1/s)
                  deltaU2,                           & ! Square of resolved velocity difference (m2/s2)
                  ws_cntr=Ws_1d,                     & ! Turbulent velocity scale profile (m/s)
                  N_iface=N_1d)                        ! Buoyancy frequency (1/s)


      surfBuoyFlux = buoyFlux(i,j,1) ! This is only used in kpp_compute_OBL_depth to limit
                                     ! h to Monin-Obukov (default is false, ie. not used)

      call CVmix_kpp_compute_OBL_depth( &
        BulkRi_1d,              & ! (in) Bulk Richardson number
        iFaceHeight,            & ! (in) Height of interfaces (m)
        OBLdepth_0d,            & ! (out) OBL depth (m)
        kOBL,                   & ! (out) level (+fraction) of OBL extent
        zt_cntr=cellHeight,     & ! (in) Height of cell centers (m)
        surf_fric=surfFricVel,  & ! (in) Turbulent friction velocity at surface (m/s)
        surf_buoy=surfBuoyFlux, & ! (in) Buoyancy flux at surface (m2/s3)
        Coriolis=Coriolis,      & ! (in) Coriolis parameter (1/s)
        CVmix_kpp_params_user=CS%KPP_params ) ! KPP parameters

      ! A hack to avoid KPP reaching the bottom. It was needed during development
      ! because KPP was unable to handle vanishingly small layers near the bottom.
      if (CS%deepOBLoffset>0.) then
        zBottomMinusOffset = iFaceHeight(G%ke+1) + min(CS%deepOBLoffset,-0.1*iFaceHeight(G%ke+1))
        OBLdepth_0d = min( OBLdepth_0d, -zBottomMinusOffset )
      endif

      ! apply some constraints on OBLdepth
      if(CS%fixedOBLdepth)  OBLdepth_0d = CS%fixedOBLdepth_value
      OBLdepth_0d = max( OBLdepth_0d, -iFaceHeight(2) )      ! no shallower than top layer
      OBLdepth_0d = min( OBLdepth_0d, -iFaceHeight(G%ke+1) ) ! no deeper than bottom
      kOBL        = CVmix_kpp_compute_kOBL_depth( iFaceHeight, cellHeight, OBLdepth_0d )

!*************************************************************************
! smg: remove code below

! Following "correction" step has been found to be unnecessary.
! Code should be removed after further testing.
      if (CS%correctSurfLayerAvg) then

        SLdepth_0d = CS%surf_layer_ext * OBLdepth_0d
        hTot      = h(i,j,1)
        surfTemp  = Temp(i,j,1) ; surfHtemp = surfTemp * hTot
        surfSalt  = Salt(i,j,1) ; surfHsalt = surfSalt * hTot
        surfU     = 0.5*(u(i,j,1)+u(i-1,j,1)) ; surfHu = surfU * hTot
        surfV     = 0.5*(v(i,j,1)+v(i,j-1,1)) ; surfHv = surfV * hTot
        pRef      = 0.0

        do k = 2, G%ke

          ! Recalculate differences with surface layer
          Uk = 0.5*(u(i,j,k)+u(i-1,j,k)) - surfU
          Vk = 0.5*(v(i,j,k)+v(i,j-1,k)) - surfV
          deltaU2(k) = Uk**2 + Vk**2
          pRef = pRef + GV%H_to_Pa * h(i,j,k)
          call calculate_density(surfTemp, surfSalt, pRef, rho1, EOS)
          call calculate_density(Temp(i,j,k), Salt(i,j,k), pRef, rhoK, EOS)
          deltaRho(k) = rhoK - rho1

          ! Surface layer averaging (needed for next k+1 iteration of this loop)
          if (hTot < SLdepth_0d) then
            delH = min( max(0., SLdepth_0d - hTot), h(i,j,k)*GV%H_to_m )
            hTot = hTot + delH
            surfHtemp = surfHtemp + Temp(i,j,k) * delH ; surfTemp = surfHtemp / hTot
            surfHsalt = surfHsalt + Salt(i,j,k) * delH ; surfSalt = surfHsalt / hTot
            surfHu = surfHu + 0.5*(u(i,j,k)+u(i-1,j,k)) * delH ; surfU = surfHu / hTot
            surfHv = surfHv + 0.5*(v(i,j,k)+v(i,j-1,k)) * delH ; surfV = surfHv / hTot
          endif

        enddo

        BulkRi_1d = CVmix_kpp_compute_bulk_Richardson( &
                    cellHeight(1:G%ke),                & ! Depth of cell center (m)
                    GoRho*deltaRho,                    & ! Bulk buoyancy difference, Br-B(z) (1/s)
                    deltaU2,                           & ! Square of resolved velocity difference (m2/s2)
                    ws_cntr=Ws_1d,                     & ! Turbulent velocity scale profile (m/s)
                    N_iface=N_1d )                       ! Buoyancy frequency (1/s)

        surfBuoyFlux = buoyFlux(i,j,1) ! This is only used in kpp_compute_OBL_depth to limit
                                       ! h to Monin-Obukov (default is false, ie. not used)

        call CVmix_kpp_compute_OBL_depth( &
          BulkRi_1d,              & ! (in) Bulk Richardson number
          iFaceHeight,            & ! (in) Height of interfaces (m)
          OBLdepth_0d,            & ! (out) OBL depth (m)
          kOBL,                   & ! (out) level (+fraction) of OBL extent
          zt_cntr=cellHeight,     & ! (in) Height of cell centers (m)
          surf_fric=surfFricVel,  & ! (in) Turbulent friction velocity at surface (m/s)
          surf_buoy=surfBuoyFlux, & ! (in) Buoyancy flux at surface (m2/s3)
          Coriolis=Coriolis,      & ! (in) Coriolis parameter (1/s)
          CVmix_kpp_params_user=CS%KPP_params ) ! KPP parameters

        if (CS%deepOBLoffset>0.) then
          zBottomMinusOffset = iFaceHeight(G%ke+1) + min(CS%deepOBLoffset,-0.1*iFaceHeight(G%ke+1))
          OBLdepth_0d = min( OBLdepth_0d, -zBottomMinusOffset )
          kOBL = CVmix_kpp_compute_kOBL_depth( iFaceHeight, cellHeight, OBLdepth_0d )
        endif

        ! apply some constraints on OBLdepth
        if(CS%fixedOBLdepth)  OBLdepth_0d = CS%fixedOBLdepth_value
        OBLdepth_0d = max( OBLdepth_0d, -iFaceHeight(2) )      ! no shallower than top layer
        OBLdepth_0d = min( OBLdepth_0d, -iFaceHeight(G%ke+1) ) ! no deep than bottom
        kOBL        = CVmix_kpp_compute_kOBL_depth( iFaceHeight, cellHeight, OBLdepth_0d )

      endif   ! endif for "correction" step

! smg: remove code above
! **********************************************************************


      ! Call CVMix/KPP to obtain OBL diffusivities, viscosities and non-local transports

      ! Unlike LMD94, we do not match to interior diffusivities. If using the original
      ! LMD94 shape function, not matching is equivalent to matching to a zero diffusivity.

      !BGR/ Add option for use of surface buoyancy flux with total sw flux.
      if (CS%SW_METHOD .eq. SW_METHOD_ALL_SW) then
         surfBuoyFlux = buoyFlux(i,j,1)
      elseif (CS%SW_METHOD .eq. SW_METHOD_MXL_SW) then
         surfBuoyFlux  = buoyFlux(i,j,1) - buoyFlux(i,j,int(kOBL)+1) ! We know the actual buoyancy flux into the OBL
      elseif (CS%SW_METHOD .eq. SW_METHOD_LV1_SW) then
         surfBuoyFlux  = buoyFlux(i,j,1) - buoyFlux(i,j,2)
      endif

      ! If option "MatchBoth" is selected in CVMix, MOM should be capable of matching.
      if (.not. (CS%MatchTechnique.eq.'MatchBoth')) then
         Kdiffusivity(:,:) = 0. ! Diffusivities for heat and salt (m2/s)
         Kviscosity(:)     = 0. ! Viscosity (m2/s)
      else
         Kdiffusivity(:,1) = Kt(i,j,:)
         Kdiffusivity(:,2) = Ks(i,j,:)
         Kviscosity(:)=Kv(i,j,:)
      endif

      call cvmix_coeffs_kpp(Kviscosity,        & ! (inout) Total viscosity (m2/s)
                            Kdiffusivity(:,1), & ! (inout) Total heat diffusivity (m2/s)
                            Kdiffusivity(:,2), & ! (inout) Total salt diffusivity (m2/s)
                            iFaceHeight,       & ! (in) Height of interfaces (m)
                            cellHeight,        & ! (in) Height of level centers (m)
                            Kviscosity,        & ! (in) Original viscosity (m2/s)
                            Kdiffusivity(:,1), & ! (in) Original heat diffusivity (m2/s)
                            Kdiffusivity(:,2), & ! (in) Original salt diffusivity (m2/s)
                            OBLdepth_0d,       & ! (in) OBL depth (m)
                            kOBL,              & ! (in) level (+fraction) of OBL extent
                            nonLocalTrans(:,1),& ! (out) Non-local heat transport (non-dimensional)
                            nonLocalTrans(:,2),& ! (out) Non-local salt transport (non-dimensional)
                            surfFricVel,       & ! (in) Turbulent friction velocity at surface (m/s)
                            surfBuoyFlux,      & ! (in) Buoyancy flux at surface (m2/s3)
                            G%ke,              & ! (in) Number of levels to compute coeffs for
                            G%ke,              & ! (in) Number of levels in array shape
                            CVmix_kpp_params_user=CS%KPP_params )


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
          do k = 2, G%ke
            sigma = min(1.0,-iFaceHeight(k)/OBLdepth_0d)
            nonLocalTrans(k,1) = (1.0 - sigma)**2 * (1.0 + 2.0*sigma) !*
            nonLocalTrans(k,2) = nonLocalTrans(k,1)
          enddo
        elseif (CS%NLT_shape == NLT_SHAPE_PARABOLIC) then
          do k = 2, G%ke
            sigma = min(1.0,-iFaceHeight(k)/OBLdepth_0d)
            nonLocalTrans(k,1) = (1.0 - sigma)**2 !*CS%CS2
            nonLocalTrans(k,2) = nonLocalTrans(k,1)
          enddo
        elseif (CS%NLT_shape == NLT_SHAPE_LINEAR) then
          do k = 2, G%ke
            sigma = min(1.0,-iFaceHeight(k)/OBLdepth_0d)
            nonLocalTrans(k,1) = (1.0 - sigma)!*CS%CS2
            nonLocalTrans(k,2) = nonLocalTrans(k,1)
          enddo
        elseif (CS%NLT_shape == NLT_SHAPE_CUBIC_LMD) then
          ! Sanity check (should agree with CVMix result using simple matching)
          do k = 2, G%ke
            sigma = min(1.0,-iFaceHeight(k)/OBLdepth_0d)
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
      if(CS%KPPzeroDiffusivity) then
         Kdiffusivity(:,1) = 0.0
         Kdiffusivity(:,2) = 0.0
         Kviscosity(:)     = 0.0
      endif

      ! recompute wscale for diagnostics, now that we in fact know boundary layer depth
      if (CS%id_Ws > 0) then
          call CVmix_kpp_compute_turbulent_scales( &
            -CellHeight/OBLdepth_0d,               & ! (in)  Normalized boundary layer coordinate
            OBLdepth_0d,                           & ! (in)  OBL depth (m)
            surfBuoyFlux,                          & ! (in)  Buoyancy flux at surface (m2/s3)
            surfFricVel,                           & ! (in)  Turbulent friction velocity at surface (m/s)
            w_s=Ws_1d,                             & ! (out) Turbulent velocity scale profile (m/s)
            CVmix_kpp_params_user=CS%KPP_params    & !       KPP parameters
            )
          CS%Ws(i,j,:) = Ws_1d(:)
      endif

      ! compute unresolved squared velocity for diagnostics
      if (CS%id_Vt2 > 0) then
        Vt2_1d(:) = CVmix_kpp_compute_unresolved_shear( &
                    cellHeight(1:G%ke),                 & ! Depth of cell center (m)
                    ws_cntr=Ws_1d,                      & ! Turbulent velocity scale profile, at centers (m/s)
                    N_iface=N_1d,                       & ! Buoyancy frequency at interface (1/s)
                    CVmix_kpp_params_user=CS%KPP_params ) ! KPP parameters
        CS%Vt2(i,j,:) = Vt2_1d(:)
      endif

      ! Copy 1d data into 3d diagnostic arrays
      if (CS%id_OBLdepth > 0) CS%OBLdepth(i,j) = OBLdepth_0d
      if (CS%id_BulkDrho > 0) CS%dRho(i,j,:)   = deltaRho(:)
      if (CS%id_BulkUz2 > 0)  CS%Uz2(i,j,:)    = deltaU2(:)
      if (CS%id_BulkRi > 0)   CS%BulkRi(i,j,:) = BulkRi_1d(:)
      if (CS%id_sigma > 0) then
        CS%sigma(i,j,:)  = 0.
        if (OBLdepth_0d>0.)   CS%sigma(i,j,:)  = -iFaceHeight/OBLdepth_0d
      endif
      if (CS%id_N      > 0)   CS%N(i,j,:)      = N_1d(:)
      if (CS%id_N2     > 0)   CS%N2(i,j,:)     = N2_1d(:)
      if (CS%id_Kt_KPP > 0)   CS%Kt_KPP(i,j,:) = Kdiffusivity(:,1)
      if (CS%id_Ks_KPP > 0)   CS%Ks_KPP(i,j,:) = Kdiffusivity(:,2)
      if (CS%id_Kv_KPP > 0)   CS%Kv_KPP(i,j,:) = Kviscosity(:)
      if (CS%id_Tsurf  > 0)   CS%Tsurf(i,j)    = surfTemp
      if (CS%id_Ssurf  > 0)   CS%Ssurf(i,j)    = surfSalt
      if (CS%id_Usurf  > 0)   CS%Usurf(i,j)    = surfU
      if (CS%id_Vsurf  > 0)   CS%Vsurf(i,j)    = surfv


       ! Update output of routine
      if (.not. CS%passiveMode) then
        if (CS%KPPisAdditive) then
          do k=1, G%ke+1
            Kt(i,j,k) = Kt(i,j,k) + Kdiffusivity(k,1)
            Ks(i,j,k) = Ks(i,j,k) + Kdiffusivity(k,2)
            Kv(i,j,k) = Kv(i,j,k) + Kviscosity(k)
          enddo
        else ! KPP replaces prior diffusivity when former is non-zero
          do k=1, G%ke+1
            if (Kdiffusivity(k,1) /= 0.) Kt(i,j,k) = Kdiffusivity(k,1)
            if (Kdiffusivity(k,2) /= 0.) Ks(i,j,k) = Kdiffusivity(k,2)
            if (Kviscosity(k) /= 0.) Kv(i,j,k) = Kviscosity(k)
          enddo
        endif
      endif


    ! end of the horizontal do-loops over the vertical columns
    enddo ! i
  enddo ! j


#ifdef __DO_SAFETY_CHECKS__
  if (CS%debug) then
    call hchksum(Kt, "KPP out: Kt",G%HI,haloshift=0)
    call hchksum(Ks, "KPP out: Ks",G%HI,haloshift=0)
  endif
#endif

  ! send diagnostics to post_data
  if (CS%id_OBLdepth > 0) call post_data(CS%id_OBLdepth, CS%OBLdepth,        CS%diag)
  if (CS%id_BulkDrho > 0) call post_data(CS%id_BulkDrho, CS%dRho,            CS%diag)
  if (CS%id_BulkUz2  > 0) call post_data(CS%id_BulkUz2,  CS%Uz2,             CS%diag)
  if (CS%id_BulkRi   > 0) call post_data(CS%id_BulkRi,   CS%BulkRi,          CS%diag)
  if (CS%id_sigma    > 0) call post_data(CS%id_sigma,    CS%sigma,           CS%diag)
  if (CS%id_Ws       > 0) call post_data(CS%id_Ws,       CS%Ws,              CS%diag)
  if (CS%id_N        > 0) call post_data(CS%id_N,        CS%N,               CS%diag)
  if (CS%id_N2       > 0) call post_data(CS%id_N2,       CS%N2,              CS%diag)
  if (CS%id_Vt2      > 0) call post_data(CS%id_Vt2,      CS%Vt2,             CS%diag)
  if (CS%id_uStar    > 0) call post_data(CS%id_uStar,    uStar,              CS%diag)
  if (CS%id_buoyFlux > 0) call post_data(CS%id_buoyFlux, buoyFlux,           CS%diag)
  if (CS%id_Kt_KPP   > 0) call post_data(CS%id_Kt_KPP,   CS%Kt_KPP,          CS%diag)
  if (CS%id_Ks_KPP   > 0) call post_data(CS%id_Ks_KPP,   CS%Ks_KPP,          CS%diag)
  if (CS%id_Kv_KPP   > 0) call post_data(CS%id_Kv_KPP,   CS%Kv_KPP,          CS%diag)
  if (CS%id_NLTt     > 0) call post_data(CS%id_NLTt,     nonLocalTransHeat,  CS%diag)
  if (CS%id_NLTs     > 0) call post_data(CS%id_NLTs,     nonLocalTransScalar,CS%diag)
  if (CS%id_Tsurf    > 0) call post_data(CS%id_Tsurf,    CS%Tsurf,           CS%diag)
  if (CS%id_Ssurf    > 0) call post_data(CS%id_Ssurf,    CS%Ssurf,           CS%diag)
  if (CS%id_Usurf    > 0) call post_data(CS%id_Usurf,    CS%Usurf,           CS%diag)
  if (CS%id_Vsurf    > 0) call post_data(CS%id_Vsurf,    CS%Vsurf,           CS%diag)

end subroutine KPP_calculate



!> Apply KPP non-local transport of surface fluxes for temperature.
subroutine KPP_NonLocalTransport_temp(CS, G, GV, h, nonLocalTrans, surfFlux, &
                                      dt, scalar, C_p)

  type(KPP_CS),                               intent(in)    :: CS            !< Control structure
  type(ocean_grid_type),                      intent(in)    :: G             !< Ocean grid
  type(verticalGrid_type),                    intent(in)    :: GV            !< Ocean vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: h             !< Layer/level thickness (units of H)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(in)    :: nonLocalTrans !< Non-local transport (non-dimensional)
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: surfFlux      !< Surface flux of scalar (H/s * scalar)
  real,                                       intent(in)    :: dt            !< Time-step (s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(inout) :: scalar        !< temperature
  real,                                       intent(in)    :: C_p           !< Seawater specific heat capacity (J/(kg*K))

  integer :: i, j, k
  real, dimension( SZI_(G), SZJ_(G), SZK_(G) ) :: dtracer


  dtracer(:,:,:) = 0.0
!$OMP parallel do default(none) shared(G,GV,dtracer,nonLocalTrans,h,surfFlux,CS,scalar,dt)
  do k = 1, G%ke
    do j = G%jsc, G%jec
      do i = G%isc, G%iec
        dtracer(i,j,k) = ( nonLocalTrans(i,j,k) - nonLocalTrans(i,j,k+1) ) / &
                         ( h(i,j,k) + GV%H_subroundoff ) * surfFlux(i,j)
      enddo
    enddo
  enddo

  !  Update tracer due to non-local redistribution of surface flux
  if (CS%applyNonLocalTrans) then
    do k = 1, G%ke
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
    do k = 1, G%ke
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
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)    :: h             !< Layer/level thickness (units of H)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(in)    :: nonLocalTrans !< Non-local transport (non-dimensional)
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: surfFlux      !< Surface flux of scalar (H/s * scalar)
  real,                                       intent(in)    :: dt            !< Time-step (s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(inout) :: scalar        !< Scalar (scalar units)

  integer :: i, j, k
  real, dimension( SZI_(G), SZJ_(G), SZK_(G) ) :: dtracer


  dtracer(:,:,:) = 0.0
!$OMP parallel do default(none) shared(G,GV,dtracer,nonLocalTrans,h,surfFlux,CS,scalar,dt)
  do k = 1, G%ke
    do j = G%jsc, G%jec
      do i = G%isc, G%iec
        dtracer(i,j,k) = ( nonLocalTrans(i,j,k) - nonLocalTrans(i,j,k+1) ) / &
                         ( h(i,j,k) + GV%H_subroundoff ) * surfFlux(i,j)
      enddo
    enddo
  enddo

  !  Update tracer due to non-local redistribution of surface flux
  if (CS%applyNonLocalTrans) then
    do k = 1, G%ke
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
    do k = 1, G%ke
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

  deallocate(CS)
end subroutine KPP_end

!> \namespace mom_kpp
!!
!! \section section_KPP The K-Profile Parameterization
!!
!! The K-Profile Parameterization (KPP) of Large et al., 1994, (http://dx.doi.org/10.1029/94RG01872) is
!! implemented via the Community Vertical Mixing package, [CVmix](https://code.google.com/p/cvmix),
!! which is called directly by this module.
!!
!! The formulation and implementation of KPP is described in great detail in the
!! [CVMix manual](https://cvmix.googlecode.com/svn/trunk/manual/cvmix.pdf) (written by our own Stephen Griffies).
!!
!! \subsection section_KPP_nutshell KPP in a nutshell
!!
!! Large et al., 1994, decompose the parameterized boundary layer turbulent flux of a scalar, \f$ s \f$, as
!! \f[ \overline{w^\prime s^\prime} = -K \partial_z s + K \gamma_s(\sigma), \f]
!! where \f$ \sigma = -z/h \f$ is a non-dimensional coordinate within the boundary layer of depth \f$ h \f$.
!! \f$ K \f$ is the eddy diffusivity and is a function of position within the boundary layer as well as a
!! function of the surface forcing:
!! \f[ K = h w_s(\sigma) G(\sigma) . \f]
!! Here, \f$ w_s \f$ is the vertical velocity scale of the boundary layer turbulence and \f$ G(\sigma) \f$ is
!! a "shape function" which is described later.
!! The last term is the "non-local transport" which involves a function \f$ \gamma_s(\sigma) \f$ that is matched
!! to the forcing but is not actually needed in the final implementation.
!! Instead, the entire non-local transport term can be equivalently written
!! \f[ K \gamma_s(\sigma) = C_s G(\sigma) Q_s \f]
!! where \f$ Q_s \f$ is the surface flux of \f$ s \f$ and \f$ C_s \f$ is a constant.
!! The vertical structure of the redistribution (non-local) term is solely due  to the shape function, \f$ G(\sigma) \f$.
!! In our implementation of KPP, we allow the shape functions used for \f$ K \f$ and for the non-local transport
!! to be chosen independently.
!!
!! [google_thread_NLT]: https://groups.google.com/forum/#!msg/cvmix-dev/i6rF-eHOtKI/Ti8BeyksrhAJ "Extreme values of non-local transport"
!!
!! The particular shape function most widely used in the atmospheric community is
!! \f[ G(\sigma) = \sigma (1-\sigma)^2 \f]
!! which satisfies the boundary conditions
!!  \f$ G(0) = 0 \f$,
!!  \f$ G(1) = 0 \f$,
!!  \f$ G^\prime(0) = 1 \f$, and
!!  \f$ G^\prime(1) = 0 \f$.
!! Large et al, 1994, alter the function so as to match interior diffusivities but we have found that this leads
!! to inconsistencies within the formulation (see google groups thread [Extreme values of non-local transport][google_thread_NLT]).
!! Instead, we use either the above form, or even simpler forms that use alternative upper boundary conditions.
!!
!! The KPP boundary layer depth is a function of the bulk Richardson number, Rib.
!! But to compute Rib, we need the boundary layer depth.  To address this circular
!! logic, we compute Rib for each vertical cell in a column, assuming the BL depth
!! equals to the depth of the given grid cell.  Once we have a vertical array of Rib(k),
!! we then call the OBLdepth routine from CVMix to compute the actual
!! OBLdepth. We optionally then "correct" the OBLdepth by cycling through once more,
!! this time knowing the OBLdepth from the first pass. This "correction" step is not
!! used by NCAR. It has been found in idealized MOM6 tests to not be necessary.
!!
!! \sa
!! kpp_calculate(), kpp_applynonlocaltransport()
end module MOM_KPP
