!> Provides the K-Profile Parameterization (KPP) of Large et al., 1994, via CVMix.
module MOM_KPP

! License goes here?

use MOM_coms, only : max_across_PEs
use MOM_checksums, only : hchksum, is_NaN
use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_ptr, post_data
use MOM_diag_mediator, only : query_averaging_enabled, register_diag_field
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_PE
use MOM_EOS, only : EOS_type, calculate_density
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_file_parser, only : openParameterBlock, closeParameterBlock
use MOM_grid, only : ocean_grid_type, isPointInCell

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

public :: KPP_init, KPP_calculate, KPP_end
public :: KPP_applyNonLocalTransport

! Enumerated constants
integer, private, parameter :: NLT_SHAPE_CVMIX     = 0 !< Use the CVmix profile
integer, private, parameter :: NLT_SHAPE_LINEAR    = 1 !< Linear, \f$ G(\sigma) = 1-\sigma \f$
integer, private, parameter :: NLT_SHAPE_PARABOLIC = 2 !< Parabolic, \f$ G(\sigma) = (1-\sigma)^2 \f$
integer, private, parameter :: NLT_SHAPE_CUBIC     = 3 !< Cubic, \f$ G(\sigma) = 1 + (2\sigma-3) \sigma^2\f$ 
integer, private, parameter :: NLT_SHAPE_CUBIC_LMD = 4 !< Original shape, \f$ G(\sigma) = \frac{27}{4} \sigma (1-\sigma)^2 \f$

!> Control structure for containing KPP parameters/data
type, public :: KPP_CS ; private

  ! Parameters
  real    :: Ri_crit              !< Critical bulk Richardson number (defines OBL depth)
  real    :: vonKarman            !< von Karman constant (dimensionless)
  real    :: cs                   !< Parameter for computing velocity scale function (dimensionless)
  character(len=10) :: interpType !< Type of interpolation to use in determining OBL
  logical :: computeEkman         !< If True, compute Ekman depth limit
  logical :: computeMoninObukhov  !< If True, compute Monin-Obukhov limit
  logical :: passiveMode          !< If True, makes KPP passive meaning it does NOT alter the diffusivity
  logical :: applyNonLocalTrans   !< If True, apply non-local transport to heat and scalars
  real    :: deepOBLoffset        !< If non-zero, is a distance from the bottom that the OBL can not penetrate through (m)
  real    :: minOBLdepth          !< If non-zero, is a minimum depth for the OBL (m)
  logical :: debug                !< If True, calculate checksums and write debugging information
  logical :: correctSurfLayerAvg  !< If true, applies a correction to the averaging of surface layer properties
  real    :: surfLayerDepth       !< A guess at the depth of the surface layer (which should 0.1 of OBLdepth) (m)
  integer :: NLT_shape            !< Determines the shape function for nonlocal transport.
  logical :: KPPzeroDiffusivity   !< If True, will set diffusivity and viscosity from KPP to zero; for testing purposes. 
  logical :: KPPisAdditive        !< If True, will add KPP diffusivity to initial diffusivity.
                                  !! If False, will replace initial diffusivity wherever KPP diffusivity is non-zero.

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
  integer :: id_NLTt     = -1, id_NLTs     = -1
  integer :: id_Tsurf    = -1, id_Ssurf    = -1
  integer :: id_Usurf    = -1, id_Vsurf    = -1
  integer :: id_dSdt     = -1, id_dTdt     = -1
  integer :: id_Kd_in    = -1

  ! Diagnostics arrays
  real, allocatable, dimension(:,:)   :: &
 OBLdepth  !< Depth (positive) of OBL (m)
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
  real, allocatable, dimension(:,:)   :: Usurf     !< i-component of velocity for surface layer (m/s)
  real, allocatable, dimension(:,:)   :: Vsurf     !< j-component of velocity for surface layer (m/s)

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
  character(len=40) :: mod = 'MOM_KPP' ! This module's name.
  character(len=20) :: string ! A local temporary string

  if (associated(CS)) call MOM_error(FATAL, 'MOM_KPP, KPP_init: '// &
           'Control structure has already been initialized')
  allocate(CS)

! Read parameters
  call log_version(paramFile, mod, version, 'This is the MOM wrapper to CVmix:KPP\n' // &
            'See http://code.google.com/p/cvmix/')
  call get_param(paramFile, mod, "USE_KPP", KPP_init, &
                 "If true, turns on the [CVmix] KPP scheme of Large et al., 1984,\n"// &
                 "to calculate diffusivities and non-local transport in the OBL.", &
                 default=.false.)

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
  call get_param(paramFile, mod, 'RI_CRIT', CS%Ri_crit,                       &
                 'Critical bulk Richardson number used to define depth of the\n'// &
                 'Oceab Boundary Layer (OBL).',                               &
                 units='nondim', default=0.3)
  call get_param(paramFile, mod, 'VON_KARMAN', CS%vonKarman, &
                 'von Karman constant.',                     &
                 units='nondim', default=0.40)
  call get_param(paramFile, mod, 'INTERP_TYPE', CS%interpType,                  &
                 'Type of interpolation to use to determine the OBL depth.\n'// &
                 'Allowed types are: linear, quadratic, cubic.',                &
                 default='quadratic')
  call get_param(paramFile, mod, 'COMPUTE_EKMAN', CS%computeEkman,                     &
                 'If True, limit the OBL depth to be shallower than the Ekman depth.', &
                 default=.False.)
  call get_param(paramFile, mod, 'COMPUTE_MONIN_OBUKHOV', CS%computeMoninObukhov,  &
                 'If True, limit the OBL depth to be shallower than the\n'//       &
                 'Monin-Obukhov depth.',                                           &
                 default=.False.)
  call get_param(paramFile, mod, 'CS', CS%cs, &
                 'Parameter for computing velocity scale function.', &
                 units='nondim', default=98.96)
  call get_param(paramFile, mod, 'DEEP_OBL_OFFSET', CS%deepOBLoffset, &
                 'If non-zero, the distance above the bottom to which the OBL is clipped\n'// &
                 'if it would otherwise reach the bottom. The smaller of this and 0.1D is used.', &
                 units='m',default=0.)
  call get_param(paramFile, mod, 'MINIMUM_OBL_DEPTH', CS%minOBLdepth, &
                 'If non-zero, a minimum depth to use for KPP OBL depth. Independent of\n'// &
                 'this parameter, the OBL depth is always at least as deep as the first layer.', &
                 units='m',default=0.)
  call get_param(paramFile, mod, 'CORRECT_SURFACE_LAYER_AVERAGE', CS%correctSurfLayerAvg, &
                 'If true, applies a correction step to the averaging of surface layer\n'// &
                 'properties.', default=.False.)
  call get_param(paramFile, mod, 'FIRST_GUESS_SURFACE_LAYER_DEPTH', CS%surfLayerDepth, &
                 'The first guess at the depth of the surface layer used for averaging\n'// &
                 'the surface layer properties. If =0, the top model level properties\n'//&
                 'will be used for the surface layer. If CORRECT_SURFACE_LAYER_AVERAGE=True, a\n'// &
                 'subsequent correction is applied.', units='m', default=0.)
  call get_param(paramFile, mod, 'NLT_SHAPE', string, &
                 'The shape of the nonlocal transport (or redistribution of surface\n'//&
                 'forcina. Allowed values are: \n'//&
                 '\t CVMIX     - Uses the profile from CVmix\n'//&
                 '\t LINEAR    - A linear profile, 1-sigma\n'//&
                 '\t PARABOLIC - A paroblic profile, (1-sigma)^2\n'//&
                 '\t CUBIC     - A cubic profile, (1-sigma)^2(1+2*sigma)\n'//&
                 '\t CUBIC_LMD - The original KPP profile', &
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
  call get_param(paramFile, mod, 'KPP_ZERO_DIFFUSIVITY', CS%KPPzeroDiffusivity, &
                 'If true, sets both the diffusivity and viscosity from KPP to zero; for testing.',&
                 default=.False.)
  call get_param(paramFile, mod, 'KPP_IS_ADDITIVE', CS%KPPisAdditive, &
                 'If true, adds KPP diffusivity to the existing diffusivity. If false, replaces '//&
                 'exisiting diffusivity with KPP diffusivity wherever the latter is non-zero.',&
                 default=.False.)

  call closeParameterBlock(paramFile)
  call get_param(paramFile, mod, 'DEBUG', CS%debug, default=.False., do_not_log=.True.)

! Forego remainder of initialization if not using this scheme
  if (.not. KPP_init) return

  call CVmix_init_kpp( Ri_crit=CS%Ri_crit,                 &
                       vonKarman=CS%vonKarman,             &
                       interp_type=CS%interpType,          &
                       lEkman=CS%computeEkman,             &
                       lMonOb=CS%computeMoninObukhov,      &
                       MatchTechnique='SimpleShapes',      &
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
      'Square of bulk difference in resolved velocity used in Bulk Richardson number, as used by [CVmix] KPP', 'm2/s2')
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
      'Frictional velocity, u*, as used by [CVmix] KPP', 'm/s')
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
  CS%id_dTdt = register_diag_field('ocean_model', 'KPP_dTdt', diag%axesTL, Time, &
      'Temperature tendency due to non-local transport of heat, as calculated by [CVmix] KPP', 'K/s')
  CS%id_dSdt = register_diag_field('ocean_model', 'KPP_dSdt', diag%axesTL, Time, &
      'Salinity tendency due to non-local transport of heat, as calculated by [CVmix] KPP', 'ppt/s')
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
  if (CS%id_Usurf > 0)    CS%Tsurf(:,:) = 0.
  if (CS%id_Vsurf > 0)    allocate( CS%Vsurf( SZI_(G), SZJB_(G)) )
  if (CS%id_Vsurf > 0)    CS%Ssurf(:,:) = 0.

end function KPP_init


!> Calculates diffusivity and non-local transport according to KPP parameterization 
subroutine KPP_calculate(CS, G, h, Temp, Salt, u, v, EOS, uStar, buoyFlux, Kt, Ks, Kv, nonLocalTransHeat, nonLocalTransScalar)

! Arguments
  type(KPP_CS),                           pointer       :: CS    !< Control structure
  type(ocean_grid_type),                  intent(in)    :: G     !< Ocean grid
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: h     !< Layer/level thicknesses (units of H)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: Temp  !< Pot. temperature (degrees C)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: Salt  !< Salinity (ppt)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in)    :: u     !< Velocity components (m/s)
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in)    :: v     !< Velocity components (m/s)
  type(EOS_type),                         pointer       :: EOS   !< Equation of state
  real, dimension(NIMEM_,NJMEM_),         intent(in)    :: uStar !< Piston velocity (m/s)
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_), intent(in)    :: buoyFlux !< Forcing buoyancy flux (m2/s3)
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_), intent(inout) :: Kt !< (in) Vertical diffusivity of heat in interior (m2/s)
                                                                    !< (out) Vertical diffusivity including KPP (m2/s)
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_), intent(inout) :: Ks !< (in) Vertical diffusivity of salt in interior (m2/s)
                                                                    !< (out) Vertical diffusivity including KPP (m2/s)
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_), intent(inout) :: Kv !< (in) Vertical viscosity in interior (m2/s)
                                                                    !< (out) Vertical viscosity including KPP (m2/s)
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_), intent(inout) :: nonLocalTransHeat   !< Temp non-local transport (m/s)
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_), intent(inout) :: nonLocalTransScalar !< scalar non-local transport (m/s)

  ! Local variables
  integer :: i, j, k, km1                               ! Loop indices
  real, dimension( G%ke )     :: cellHeight             ! Cell center heights referenced to surface (m)
  real, dimension( G%ke+1 )   :: iFaceHeight            ! Interface heights referenced to surface (m)
  real, dimension( G%ke+1 )   :: N2_1d                  ! Brunt-Vaisala frequency squared, at interfaces (1/s2)
  real, dimension( G%ke+1 )   :: N_1d                   ! (Adjusted) Brunt-Vaisala frequency, at interfaces (1/s)
  real, dimension( G%ke )     :: Ws_1d, Wm_1d           ! Profiles of vertical velocity scale for scalars/momentum (m/s)
  real, dimension( G%ke )     :: Vt2_1d                 ! Unresolved shear turbulence, at interfaces (m2/s2)
  real, dimension( G%ke )     :: BulkRi_1d              ! Bulk Richardson number for each layer
  real, dimension( G%ke )     :: deltaRho               ! delta Rho in numerator of Bulk Ri number
  real, dimension( G%ke )     :: deltaU2                ! square of delta U (shear) in denominator of Bulk Ri (m2/s2)
  real, dimension( G%ke+1, 2) :: Kdiffusivity           ! Vertical diffusivity at interfaces (m2/s)
  real, dimension( G%ke+1 )   :: Kviscosity             ! Vertical viscosity at interfaces (m2/s)
  real, dimension( G%ke+1, 2) :: nonLocalTrans          ! Non-local transport for heat/salt at interfaces (non-dimensional)
  real, dimension( G%ke )     :: surfBuoyFlux2
  real, dimension( 3*G%ke )   :: rho_1D, pres_1D, Temp_1D, Salt_1D 
  real :: kOBL, OBLdepth_0d, surfFricVel, surfBuoyFlux, Coriolis
  real :: GoRho, pRef, rho1, rhoK, rhoKm1, Uk, Vk, const1, Cv, sigma
  real :: zBottomMinusOffset                  ! Height of bottom plus a little bit (m)
  real, parameter :: eps = 0.1                ! Non-dimensional extent of Monin-Obukov surface layer. Used for const1 below.
  real, parameter :: BetaT = -0.2             ! Ratio of entrainment flux to surface buoyancy flux. Used for const1 below.
  real, parameter :: minimumVt2 = 1.e-11      ! A small number added to unresolved velocity Vt2 to avoid divide by zero.
                                              ! This value should be larger than roundoff for sensible behavior 
                                              ! with zero vertical stratification and zero resolved velocity.  
                                              ! We compute 1e-11 by the following rough scaling: 
                                              ! minimumVt2 = const1*depth*N*ws, with depth=1m, N = 1e-5 s^{-1}, ws = 1e-6 m/s
  real :: SLdepth_0d                 ! Surface layer depth. This is meant to be 0.1*OBLdepth but has to be guessed at the beginning (m)
  real :: hTot, delH                 ! The running total of thickness used in the surface layer average (m), the thickness from this layer (m)
  real :: surfHtemp, surfTemp        ! Integral and average of temperature over the surface layer
  real :: surfHsalt, surfSalt        ! Integral and average of salinity over the surface layer
  real :: hTotU, surfHu, surfU       ! Thickness, integral and average of U over the surface layer
  real :: hTotV, surfHv, surfV       ! Thickness, integral and average of V over the surface layer
  real :: hTotUm1, surfHum1, surfUm1 ! Same fo i-1,j or U
  real :: hTotVm1, surfHvm1, surfVm1 ! Same fo i,j-1 of V
  integer :: kk
#ifdef __DO_SAFETY_CHECKS__
  if (CS%debug) then
    call hchksum(h, "KPP in: h",G,haloshift=0)
    call hchksum(Temp, "KPP in: T",G,haloshift=0)
    call hchksum(Salt, "KPP in: S",G,haloshift=0)
    call hchksum(u, "KPP in: u",G,haloshift=0)
    call hchksum(v, "KPP in: v",G,haloshift=0)
    call hchksum(uStar, "KPP in: uStar",G,haloshift=0)
    call hchksum(buoyFlux, "KPP in: buoyFlux",G,haloshift=0)
    call hchksum(Kt, "KPP in: Kt",G,haloshift=0)
    call hchksum(Ks, "KPP in: Ks",G,haloshift=0)
  endif
#endif

  GoRho = G%g_Earth / G%Rho0
  ! const1 is a constant factor in unresolved squared velocity, Vt2 (eq. 23 in LMD94)
  const1 = sqrt( abs(BetaT) / (CS%cs * eps) )/( CS%Ri_crit * (CS%vonKarman**2) )
  nonLocalTrans(:,:) = 0.0

  if (CS%id_Kd_in > 0) call post_data(CS%id_Kd_in, Kt, CS%diag)

!$OMP parallel do default(none) shared(G,CS,EOS,uStar,Temp,Salt,u,v,h,GoRho,buoyFlux, &
!$OMP                                  const1,nonLocalTransHeat,                      &
!$OMP                                  nonLocalTransScalar,Kt,Ks,Kv)                  &
!$OMP                     firstprivate(nonLocalTrans)                                 &
!$OMP                          private(Coriolis,surfFricVel,SLdepth_0d,hTot,surfTemp, &
!$OMP                                  surfHtemp,surfSalt,surfHsalt,hTotU,surfU,      &
!$OMP                                  surfHu,hTotV,surfV,surfHv,hTotUm1,surfUm1,     &
!$OMP                                  surfHum1,hTotVm1,surfVm1,surfHvm1,iFaceHeight, &
!$OMP                                  pRef,km1,cellHeight,Uk,Vk,deltaU2,             &
!$OMP                                  rho1,rhoK,rhoKm1,deltaRho,N2_1d,N_1d,delH,     &
!$OMP                                  surfBuoyFlux,Ws_1d,Cv,Vt2_1d,BulkRi_1d,        &
!$OMP                                  OBLdepth_0d,zBottomMinusOffset,Kdiffusivity,   &
!$OMP                                  Kviscosity,sigma,kOBL,kk,pres_1D,Temp_1D,      &
!$OMP                                  Salt_1D,rho_1D,surfBuoyFlux2)
  do j = G%jsc, G%jec
    do i = G%isc, G%iec
      if (G%mask2dT(i,j)==0.) cycle ! Skip calling KPP for land points

      ! Things that are independent of position within the column
      Coriolis = 0.25*( (G%CoriolisBu(i,j)   + G%CoriolisBu(i-1,j-1)) &
                       +(G%CoriolisBu(i-1,j) + G%CoriolisBu(i,j-1)) )
      surfFricVel = uStar(i,j)

      ! Initialize the surface properties to layer k=1 values and initialize the running integrals
      SLdepth_0d = CS%surfLayerDepth ! This is a first guess at the surface layer depth (which we do not know yet)
      hTot = 1.e-16 ! We initialize to non-zero to avoid divide by zero within the k-loop
      surfTemp = Temp(i,j,1) ; surfHtemp = surfTemp * hTot
      surfSalt = Salt(i,j,1) ; surfHsalt = surfSalt * hTot
      hTotU = hTot ; surfU = u(i,j,1) ; surfHu = surfU * hTotU
      hTotV = hTot ; surfV = v(i,j,1) ; surfHv = surfV * hTotV
      hTotUm1 = hTot ; surfUm1 = u(i-1,j,1) ; surfHum1 = surfUm1 * hTotUm1
      hTotVm1 = hTot ; surfVm1 = v(i,j-1,1) ; surfHvm1 = surfVm1 * hTotVm1

      ! This k-loop calculates quantities that will be passed to KPP
      iFaceHeight(1) = 0.
      pRef = 0.
      do k = 1, G%ke
        km1 = max(1, k-1)
        kk = 3*(k-1)
        pres_1D(kk+1) = pRef
        pres_1D(kk+2) = pRef
        pres_1D(kk+3) = pRef
        Temp_1D(kk+1) = surfTemp
        Temp_1D(kk+2) = Temp(i,j,k)
        Temp_1D(kk+3) = Temp(i,j,km1)
        Salt_1D(kk+1) = surfsalt
        Salt_1D(kk+2) = Salt(i,j,k)
        Salt_1D(kk+3) = Salt(i,j,km1)
        ! Compute heights, referenced to the surface (z=0)
        cellHeight(k) = iFaceHeight(k) - 0.5 * h(i,j,k) * G%H_to_m ! cell center in meters 
        iFaceHeight(k+1) = iFaceHeight(k) - h(i,j,k) * G%H_to_m ! cell bottom in meters

        ! Compute Bulk Richardson number
        ! rho1 is meant to be the average over the surface layer (0.1*OBLdepth) which
        ! we do not know at this point. So we compute a running average down to a
        ! prescribed depth, as a first guess.
        ! In z-mode, this will typically just be the top level. but a proper integral
        ! will be needed for fine vertical resolution or arbitrary coordinates.   ???????

        ! Compute shear between surface layer and this layer for use in the Bulk Richardson number
        Uk = 0.5 * ( abs( u(i,j,k) - surfU ) + abs( u(i-1,j,k) - surfUm1 ) ) ! delta_k U  w/ C-grid average
        Vk = 0.5 * ( abs( v(i,j,k) - surfV ) + abs( v(i,j-1,k) - surfVm1 ) ) ! delta_k V  w/ C-grid average
        deltaU2(k) = Uk**2 + Vk**2

        ! Pressure at bottom of level k will become pressure at top of level on next iteration
        pRef = pRef + G%g_Earth * G%Rho0 * h(i,j,k) * G%H_to_m ! Boussinesq approximation!!!! ?????

        ! Surface layer averaging (needed for next k+1 iteration of this loop)
        if (hTot < SLdepth_0d) then
          delH = min( max(0., SLdepth_0d - hTot), h(i,j,k)*G%H_to_m )
          hTot = hTot + delH
          surfHtemp = surfHtemp + Temp(i,j,k) * delH ; surfTemp = surfHtemp / hTot
          surfHsalt = surfHsalt + Salt(i,j,k) * delH ; surfSalt = surfHsalt / hTot
        endif
        if (hTotU < SLdepth_0d) then
          delH = min( max(0., SLdepth_0d - hTotU), 0.5*(h(i,j,k)+h(i+1,j,k))*G%H_to_m )
          hTotU = hTotU + delH
          surfHu = surfHu + u(i,j,k) * delH ; surfU = surfHu / hTotU
        endif
        if (hTotV < SLdepth_0d) then
          delH = min( max(0., SLdepth_0d - hTotV), 0.5*(h(i,j,k)+h(i,j+1,k))*G%H_to_m )
          hTotV = hTotV + delH
          surfHv = surfHv + v(i,j,k) * delH ; surfV = surfHv / hTotV
        endif
        if (hTotUm1 < SLdepth_0d) then
          delH = min( max(0., SLdepth_0d - hTotUm1), 0.5*(h(i-1,j,k)+h(i,j,k))*G%H_to_m )
          hTotUm1 = hTotUm1 + delH
          surfHum1 = surfHum1 + u(i-1,j,k) * delH ; surfUm1 = surfHum1 / hTotUm1
        endif
        if (hTotVm1 < SLdepth_0d) then
          delH = min( max(0., SLdepth_0d - hTotVm1), 0.5*(h(i,j-1,k)+h(i,j,k))*G%H_to_m )
          hTotVm1 = hTotVm1 + delH
          surfHvm1 = surfHvm1 + v(i,j-1,k) * delH ; surfVm1 = surfHvm1 / hTotVm1
        endif
      enddo ! k
      call calculate_density(Temp_1D, Salt_1D, pres_1D, rho_1D, 1, 3*G%ke, EOS)
      do k = 1, G%ke
        km1 = max(1, k-1)
        kk = 3*(k-1)
        deltaRho(k) = rho_1D(kk+2) - rho_1D(kk+1)
        ! N2 is on interfaces.
        N2_1d(k) = ( GoRho * (rho_1D(kk+2) - rho_1D(kk+3)) ) / (0.5*(h(i,j,km1) + h(i,j,k))+G%H_subroundoff) ! Can be negative
        ! N = sqrt(N^2) but because N^2 can be negative, we clip N^2 before taking the square root   ??????
        N_1d(k) = sqrt( max( N2_1d(k), 0.) )
   !    N_1d(k) = sqrt( abs( N2_1d(k)) ) ! From MOM4p1 ?????
   !    N_1d(k) = sign( sqrt( abs( N2_1d(k) ) ), N2_1d(k) ) ! Suggestion ????
      enddo
      N2_1d( G%ke+1 ) = 0.
      N_1d( G%ke+1 ) = 0.

      do k = 1, G%ke
        surfBuoyFlux2(k) = buoyFlux(i,j,1) - buoyFlux(i,j,k+1) ! This difference accounts for penetrating of SW
      enddo

      call CVmix_kpp_compute_turbulent_scales( &
        eps,            & ! (in)  Normalized boundary layer depth; sigma = eps
        -cellHeight, & ! (in)  Guess that OBL depth (m) = -cellHeight(k)
        surfBuoyFlux2,   & ! (in)  Buoyancy flux at surface (m2/s3)
        surfFricVel,    & ! (in)  Turbulent friction velocity at surface (m/s)
        w_s=Ws_1d,   & ! (out) Turbulent velocity scale profile (m/s)
        CVmix_kpp_params_user=CS%KPP_params ) ! KPP parameters

      ! Estimate Ws in order to estimate Vt^2  (eq. 23 in LMD94)
      do k = 1, G%ke

        ! Note that if sigma > eps, then CVmix_kpp_compute_turbulent_scales 
        ! computes w_s and w_m velocity scale at sigma=eps. So we only pass
        ! sigma=eps for this calculation.    

        ! Unresolved squared velocity, Vt^2, eq 23 from LMD94
        Cv = max( 1.7, 2.1 - 200. * N_1d(k) ) ! Cv from eq A3 of Danabasoglu et al. 2006
    !   Cv = 1.8 ! MOM5 
        ! The calculation is for Vt^2 at level center but uses N from the interface below
        ! (and depth of lower interface) to bias towards higher estimates.  One would
        ! otherwise use d=-cellHeight(k) and a vertical average of N.  ?????
        Vt2_1d(k) = minimumVt2 + const1 * Cv * ( -iFaceHeight(k+1) ) * N_1d(k+1) * Ws_1d(k)

      enddo ! k

    ! The following call gives a similar answer to the above but is much less efficient
    ! Vt2_1d(:) = CVmix_kpp_compute_unresolved_shear( &
    !               iFaceHeight(2:G%ke+1), & ! Height of level centers (m) NOTE DISCREPANCY ????
    !               N_1d(2:G%ke+1),        & ! Buoyancy frequency at centers (1/s) NOTE DISCREPANCY ????
    !               Ws_1d,                 & ! Turbulent velocity scale profile, at centers (m/s)
    !               CVmix_kpp_params_user=CS%KPP_params ) ! KPP parameters
    ! Vt2_1d(:) = Vt2_1d(:) + minimumVt2

      ! Calculate Bulk Richardson number, eq 21 of LMD94
      BulkRi_1d = CVmix_kpp_compute_bulk_Richardson( &
                    iFaceHeight(1:G%ke), & ! Height of level centers (m) NOTE DISCREPANCY ????
                    GoRho*deltaRho,      & ! Bulk buoyancy difference, Br -B(z) (1/s)
                    deltaU2,             & ! Square of bulk shear (m/s)
                    Vt2_1d )               ! Square of unresolved turbulence (m2/s2)
    ! do k = 1, G%ke
    ! ! Notes:
    ! ! o BulRi(k=1)=0 because rho1=rhoK
    !   BulkRi_1d(k) = ( ( GoRho * deltaRho(k) ) * ( cellHeight(1)-cellHeight(k) ) ) / ( deltaU2(k) + Vt2_1d(k) )
    !   ! The distance here between the surface and the interface at which a stability
    !   ! calculation (and the pressure used) would take place, ie. the upper interface     ?????
    !   BulkRi_1d(k) = ( ( GoRho * deltaRho(k) ) * ( -iFaceHeight(k) ) ) / ( deltaU2(k) + Vt2_1d(k) )
    ! enddo ! k
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
      OBLdepth_0d = max( OBLdepth_0d, -iFaceHeight(2) ) ! Keep at least as deep as top layer
      OBLdepth_0d = max( OBLdepth_0d, CS%minOBLdepth ) ! Forc the OBL depth to be at least this deep
      if (CS%deepOBLoffset>0.) then
        ! This is a hack to avoid KPP reaching the bottom. It was needed during development
        ! because KPP was unable to handle vanishingly small layers near the bottom.
        zBottomMinusOffset = iFaceHeight(G%ke+1) + min(CS%deepOBLoffset,-0.1*iFaceHeight(G%ke+1))
        OBLdepth_0d = min( OBLdepth_0d, -zBottomMinusOffset )
        kOBL = CVmix_kpp_compute_kOBL_depth( iFaceHeight, cellHeight, OBLdepth_0d )
      endif
      OBLdepth_0d = min( OBLdepth_0d, -iFaceHeight(G%ke+1) ) ! Do not let the OBL pass the bottom

      ! The surface temp/salt was estimated as equal to the top layer in the first guess.
      ! Now we have an estimate of BL_depth, we can do a proper average for surfTemp and surfSalt
      if (CS%correctSurfLayerAvg) then
        SLdepth_0d = eps * OBLdepth_0d ! This is the corrected depth of the surface layer depth (m)
        hTot = h(i,j,1) ! We initialize to first layer
        surfTemp = Temp(i,j,1) ; surfHtemp = surfTemp * hTot
        surfSalt = Salt(i,j,1) ; surfHsalt = surfSalt * hTot
        hTotU = hTot ; surfU = u(i,j,1) ; surfHu = surfU * hTotU
        hTotV = hTot ; surfV = v(i,j,1) ; surfHv = surfV * hTotV
        hTotUm1 = hTot ; surfUm1 = u(i-1,j,1) ; surfHum1 = surfUm1 * hTotUm1
        hTotVm1 = hTot ; surfVm1 = v(i,j-1,1) ; surfHvm1 = surfVm1 * hTotVm1
        pRef = 0.
        do k = 2, G%ke
          ! Recalculate differences with surface layer
          Uk = 0.5 * ( abs( u(i,j,k) - surfU ) + abs( u(i-1,j,k) - surfUm1 ) ) ! delta_k U
          Vk = 0.5 * ( abs( v(i,j,k) - surfV ) + abs( v(i,j-1,k) - surfVm1 ) ) ! delta_k V
          deltaU2(k) = Uk**2 + Vk**2
          pRef = pRef + G%g_Earth * G%Rho0 * h(i,j,k) * G%H_to_m ! Boussinesq approximation!!!! ?????
          call calculate_density(surfTemp, surfSalt, pRef, rho1, EOS)
          call calculate_density(Temp(i,j,k), Salt(i,j,k), pRef, rhoK, EOS)
          deltaRho(k) = rhoK - rho1

          ! Surface layer averaging (needed for next k+1 iteration of this loop)
          if (hTot < SLdepth_0d) then
            delH = min( max(0., SLdepth_0d - hTot), h(i,j,k)*G%H_to_m )
            hTot = hTot + delH
            surfHtemp = surfHtemp + Temp(i,j,k) * delH ; surfTemp = surfHtemp / hTot
            surfHsalt = surfHsalt + Salt(i,j,k) * delH ; surfSalt = surfHsalt / hTot
          endif
          if (hTotU < SLdepth_0d) then
            delH = min( max(0., SLdepth_0d - hTotU), 0.5*(h(i,j,k)+h(i+1,j,k))*G%H_to_m )
            hTotU = hTotU + delH
            surfHu = surfHu + u(i,j,k) * delH ; surfU = surfHu / hTotU
          endif
          if (hTotV < SLdepth_0d) then
            delH = min( max(0., SLdepth_0d - hTotV), 0.5*(h(i,j,k)+h(i,j+1,k))*G%H_to_m )
            hTotV = hTotV + delH
            surfHv = surfHv + v(i,j,k) * delH ; surfV = surfHv / hTotV
          endif
          if (hTotUm1 < SLdepth_0d) then
            delH = min( max(0., SLdepth_0d - hTotUm1), 0.5*(h(i-1,j,k)+h(i,j,k))*G%H_to_m )
            hTotUm1 = hTotUm1 + delH
            surfHum1 = surfHum1 + u(i-1,j,k) * delH ; surfUm1 = surfHum1 / hTotUm1
          endif
          if (hTotVm1 < SLdepth_0d) then
            delH = min( max(0., SLdepth_0d - hTotVm1), 0.5*(h(i,j-1,k)+h(i,j,k))*G%H_to_m )
            hTotVm1 = hTotVm1 + delH
            surfHvm1 = surfHvm1 + v(i,j-1,k) * delH ; surfVm1 = surfHvm1 / hTotVm1
          else
            exit ! Avoid unnecessary calculations
          endif
        enddo

        ! Calculate Bulk Richardson number, eq 21 of LMD94
        BulkRi_1d = CVmix_kpp_compute_bulk_Richardson( &
                      iFaceHeight(1:G%ke), & ! Height of level centers (m) NOTE DISCREPANCY ????
                      GoRho*deltaRho,      & ! Bulk buoyancy difference, Br -B(z) (1/s)
                      deltaU2,             & ! Square of bulk shear (m/s)
                      Vt2_1d )               ! Square of unresolved turbulence (m2/s2)
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
        OBLdepth_0d = max( OBLdepth_0d, -iFaceHeight(2) ) ! Keep at least as deep as top layer
        if (CS%deepOBLoffset>0.) then
          zBottomMinusOffset = iFaceHeight(G%ke+1) + min(CS%deepOBLoffset,-0.1*iFaceHeight(G%ke+1))
          OBLdepth_0d = min( OBLdepth_0d, -zBottomMinusOffset )
          kOBL = CVmix_kpp_compute_kOBL_depth( iFaceHeight, cellHeight, OBLdepth_0d )
        endif
      endif

      ! Now call KPP proper to obtain BL diffusivities, viscosities and non-local transports

      ! Unlike LMD94, we do not match to interior diffusivities. If using the original
      ! LMD94 shape function, not matching is equivalent to matching to a zero diffusivity.
      Kdiffusivity(:,:) = 0. ! Diffusivities for heat and salt (m2/s)
      Kviscosity(:)     = 0. ! Viscosity (m2/s)
      surfBuoyFlux  = buoyFlux(i,j,1) - buoyFlux(i,j,int(kOBL)+1) ! We know the actual buoyancy flux into the OBL
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
                            CVmix_kpp_params_user=CS%KPP_params )


      ! Over-write CVMix shape function with one of the following choices. 
      ! Note that nonLocalTrans = Cs * G(sigma) (LMD94 notation), with 
      ! Cs = 6.32739901508.
      ! Start do-loop at k=2, since k=1 is ocean surface (sigma=0) 
      ! and we do not wish to double-count the surface forcing.  
      ! Only compute nonlocal transport for 0 <= sigma <= 1. 
      ! Recommended shape is the parabolic; it gives deeper boundary layer. 
      if (surfBuoyFlux < 0.0) then
        if (CS%NLT_shape == NLT_SHAPE_CUBIC) then
          do k = 2, G%ke
            sigma = min(1.0,-iFaceHeight(k)/OBLdepth_0d)
           !nonLocalTrans(k,1) = 1.0 + (2.0*sigma-3)*sigma**2
            nonLocalTrans(k,1) = (1.0 - sigma)**2 * (1.0 + 2.0*sigma)
            nonLocalTrans(k,2) = nonLocalTrans(k,1)
          enddo 
        elseif (CS%NLT_shape == NLT_SHAPE_PARABOLIC) then
          do k = 2, G%ke
            sigma = min(1.0,-iFaceHeight(k)/OBLdepth_0d)
            nonLocalTrans(k,1) = (1.0 - sigma)**2
            nonLocalTrans(k,2) = nonLocalTrans(k,1)
          enddo 
        elseif (CS%NLT_shape == NLT_SHAPE_LINEAR) then
          do k = 2, G%ke
            sigma = min(1.0,-iFaceHeight(k)/OBLdepth_0d)
            nonLocalTrans(k,1) = (1.0 - sigma)
            nonLocalTrans(k,2) = nonLocalTrans(k,1)
          enddo 
        elseif (CS%NLT_shape == NLT_SHAPE_CUBIC_LMD) then
          ! Sanity check (should agree with CVMix result using simple matching)
          do k = 2, G%ke
            sigma = min(1.0,-iFaceHeight(k)/OBLdepth_0d)
            nonLocalTrans(k,1) = 6.32739901508 * sigma*(1.0 -sigma)**2
            nonLocalTrans(k,2) = nonLocalTrans(k,1)
          enddo 
        endif 
      endif 

      nonLocalTransHeat(i,j,:)   = nonLocalTrans(:,1) ! correct index ???
      nonLocalTransScalar(i,j,:) = nonLocalTrans(:,2) ! correct index ???

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

      ! Copy 1d data into 3d diagnostic arrays
      if (CS%id_OBLdepth > 0) CS%OBLdepth(i,j) = OBLdepth_0d
      if (CS%id_BulkDrho > 0) CS%dRho(i,j,:)   = deltaRho(:)
      if (CS%id_BulkUz2 > 0)  CS%Uz2(i,j,:)    = deltaU2(:)
      if (CS%id_BulkRi > 0)   CS%BulkRi(i,j,:) = BulkRi_1d(:)
      if (CS%id_sigma > 0)    CS%sigma(i,j,:)  = -iFaceHeight/OBLdepth_0d
      if (CS%id_N > 0)        CS%N(i,j,:)      = N_1d(:)
      if (CS%id_N2 > 0)       CS%N2(i,j,:)     = N2_1d(:)
      if (CS%id_Vt2 > 0)      CS%Vt2(i,j,:)    = Vt2_1d(:)
      if (CS%id_Kt_KPP > 0)   CS%Kt_KPP(i,j,:) = Kdiffusivity(:,1)
      if (CS%id_Ks_KPP > 0)   CS%Ks_KPP(i,j,:) = Kdiffusivity(:,2)
      if (CS%id_Kv_KPP > 0)   CS%Kv_KPP(i,j,:) = Kviscosity(:)
      if (CS%id_Tsurf > 0)    CS%Tsurf(i,j)    = surfTemp
      if (CS%id_Ssurf > 0)    CS%Ssurf(i,j)    = surfSalt
      if (CS%id_Usurf > 0)    CS%Usurf(i,j)    = surfU
      if (CS%id_Vsurf > 0)    CS%Vsurf(i,j)    = surfv

!if (abs(G%geoLonT(i,j)+80.99621)+abs(G%geoLatT(i,j)-82.64066)<0.5) then
! write(0,*) G%geoLonT(i,j), G%geoLatT(i,j), isPointInCell(G,i,j,-80.99621,82.64066)
!endif
!     if (isPointInCell(G,i,j,-80.99621,82.64066)) then
      if (maxval(abs(nonLocalTrans(:,1)))>8.) then
        write(*,'(3(a,es10.3,x),es10.3)') 'Lon=',G%geoLonT(i,j),'Lat=',G%geoLatT(i,j),'Depth=',G%bathyT(i,j),-iFaceHeight(G%ke+1)
        write(*,'(7(a,1x,1es10.3,1x))') 'OBL_depth=',OBLdepth_0d,'u*=',surfFricVel,'buoySurf=',surfBuoyFlux,'kOBL=',kOBL
        write(*,'(a4,5x,12a13)') 'k','zw','N^2','Vt^2','Kt (in)','K (out)','NLT'
        write(*,'(a4,12a13)') 'k','zT','T','S','dB*d','dU2','Ws','Rib','div.NLT'
        do k=1,min(G%ke,int(kOBL)+3)
        write(*,'(f4.1,5x,12(1x,1es12.3))') float(k)-0.5,iFaceHeight(k), &
           N2_1d(k),Vt2_1d(k),Kt(i,j,k),Kdiffusivity(k,1),nonLocalTrans(k,1)
        write(*,'(f4.1,12(1x,1es12.3))') float(k),cellHeight(k), &
           Temp(i,j,k),Salt(i,j,k),GoRho*deltaRho(k)*(cellHeight(1)-cellHeight(k)),deltaU2(k),Ws_1d(k), &
           BulkRi_1d(k),(nonLocalTrans(k,1)-nonLocalTrans(k+1,1))/h(i,j,k)
        enddo
        write(*,'(f4.1,5x,12(1x,1es12.3))') float(k)-0.5,iFaceHeight(k), &
           N2_1d(k),Vt2_1d(k),Kt(i,j,k),Kdiffusivity(k,1),nonLocalTrans(k,1)
      endif

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

    enddo ! i
  enddo ! j

#ifdef __DO_SAFETY_CHECKS__
  if (CS%debug) then
    call hchksum(Kt, "KPP out: Kt",G,haloshift=0)
    call hchksum(Ks, "KPP out: Ks",G,haloshift=0)
  endif
#endif

  if (CS%id_OBLdepth > 0) call post_data(CS%id_OBLdepth, CS%OBLdepth, CS%diag)
  if (CS%id_BulkDrho > 0) call post_data(CS%id_BulkDrho, CS%dRho, CS%diag)
  if (CS%id_BulkUz2 > 0)  call post_data(CS%id_BulkUz2, CS%Uz2, CS%diag)
  if (CS%id_BulkRi > 0)   call post_data(CS%id_BulkRi, CS%BulkRi, CS%diag)
  if (CS%id_sigma > 0)    call post_data(CS%id_sigma, CS%sigma, CS%diag)
  if (CS%id_Ws > 0)       call post_data(CS%id_Ws, CS%Ws, CS%diag)
  if (CS%id_N > 0)        call post_data(CS%id_N, CS%N, CS%diag)
  if (CS%id_N2 > 0)       call post_data(CS%id_N2, CS%N2, CS%diag)
  if (CS%id_Vt2 > 0)      call post_data(CS%id_Vt2, CS%Vt2, CS%diag)
  if (CS%id_uStar > 0)    call post_data(CS%id_uStar, uStar, CS%diag)
  if (CS%id_buoyFlux > 0) call post_data(CS%id_buoyFlux, buoyFlux, CS%diag)
  if (CS%id_Kt_KPP > 0)   call post_data(CS%id_Kt_KPP, CS%Kt_KPP, CS%diag)
  if (CS%id_Ks_KPP > 0)   call post_data(CS%id_Ks_KPP, CS%Ks_KPP, CS%diag)
  if (CS%id_Kv_KPP > 0)   call post_data(CS%id_Kv_KPP, CS%Kv_KPP, CS%diag)
  if (CS%id_NLTt > 0)     call post_data(CS%id_NLTt, nonLocalTransHeat, CS%diag)
  if (CS%id_NLTs > 0)     call post_data(CS%id_NLTs, nonLocalTransScalar, CS%diag)
  if (CS%id_Tsurf > 0)    call post_data(CS%id_Tsurf, CS%Tsurf, CS%diag)
  if (CS%id_Ssurf > 0)    call post_data(CS%id_Ssurf, CS%Ssurf, CS%diag)
  if (CS%id_Usurf > 0)    call post_data(CS%id_Usurf, CS%Usurf, CS%diag)
  if (CS%id_Vsurf > 0)    call post_data(CS%id_Vsurf, CS%Vsurf, CS%diag)

end subroutine KPP_calculate


!subroutine fixNLTamplitude( h, NLT )
!! Arguments
!  real, dimension(:), intent(in)    :: h
!  real, dimension(:), intent(inout) :: NLT
!! Local variables
!  real :: maxDelta
!  integer :: n
!  n=size(NLT,1)
!  maxDelta = maxval( abs( (NLT(1:n-1)-NLT(2:n)) )/h(:) )
!  if (maxDelta>1.) NLT = NLT / maxDelta
!end subroutine fixNLTamplitude


!> Applies the KPP non-local transport of surface fluxes; only available for tracers
subroutine KPP_applyNonLocalTransport(CS, G, h, nonLocalTrans, surfFlux, dt, scalar, isHeat, isSalt)

  type(KPP_CS),                                 intent(in)    :: CS            !< Control structure
  type(ocean_grid_type),                        intent(in)    :: G             !< Ocean grid
  real, dimension(NIMEM_,NJMEM_,NKMEM_),        intent(in)    :: h             !< Layer/level thicknesses (units of H)
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_), intent(in)    :: nonLocalTrans !< Non-local transport (non-dimensional)
  real, dimension(NIMEM_,NJMEM_),               intent(in)    :: surfFlux      !< Surface source of scalar (m/s * scalar)
  real,                                         intent(in)    :: dt            !< Time-step (s)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),        intent(inout) :: scalar        !< Scalar or temperature (scalar units)
  logical, optional,                            intent(in)    :: isHeat        !< Indicates scalar is heat for diagnostics
  logical, optional,                            intent(in)    :: isSalt        !< Indicates scalar is salt for diagnostics

  integer :: i, j, k
  logical :: diagHeat, diagSalt
  logical :: debugColumn
  real, dimension( SZI_(G), SZJ_(G), SZK_(G) ) :: dSdt ! Tendency in scalar due to non-local transport (scalar/s)

  diagHeat = .False.
  if (present(isHeat)) then
    if (isHeat) then
      if (CS%id_dTdt > 0) diagHeat = .True.
      if (CS%id_QminusSW > 0) call post_data(CS%id_QminusSW, surfFlux, CS%diag)
    endif
  endif
  diagSalt = .False.
  if (present(isSalt)) then
    if (isSalt) then
      if (CS%id_dSdt > 0) diagSalt = .True.
      if (CS%id_netS > 0) call post_data(CS%id_netS, surfFlux, CS%diag)
    endif
  endif

  if (diagHeat .or. diagSalt) dSdt(:,:,:) = 0. ! Zero halos for diagnostics ???
!$OMP parallel do default(none) shared(G,dSdt,nonLocalTrans,h,surfFlux,CS,scalar,dt)
  do k = 1, G%ke
    do j = G%jsc, G%jec
      do i = G%isc, G%iec
        ! Tendency due to non-local transport of scalar
        dSdt(i,j,k) = ( nonLocalTrans(i,j,k) - nonLocalTrans(i,j,k+1) ) / ( h(i,j,k) + G%H_subroundoff ) * surfFlux(i,j)
        ! Update the scalar
        if (CS%applyNonLocalTrans) scalar(i,j,k) = scalar(i,j,k) + dt * dSdt(i,j,k)
      enddo ! i
    enddo ! j
  enddo ! k

  if (diagHeat) call post_data(CS%id_dTdt, dSdt, CS%diag)
  if (diagSalt) call post_data(CS%id_dSdt, dSdt, CS%diag)

end subroutine KPP_applyNonLocalTransport


!> Clear pointers, deallocate memory
subroutine KPP_end(CS)
  type(KPP_CS), pointer :: CS !< Control structure

  deallocate(CS)
end subroutine KPP_end

!> \class mom_kpp
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
!! \sa
!! kpp_calculate(), kpp_applynonlocaltransport()
end module MOM_KPP
