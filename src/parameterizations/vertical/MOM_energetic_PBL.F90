!> Energetically consistent planetary boundary layer parameterization
module MOM_energetic_PBL

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,      only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_coms,           only : EFP_type, real_to_EFP, EFP_to_real, operator(+), assignment(=), EFP_sum_across_PEs
use MOM_debugging,      only : hchksum
use MOM_diag_mediator,  only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator,  only : time_type, diag_ctrl
use MOM_domains,        only : create_group_pass, do_group_pass, group_pass_type
use MOM_error_handler,  only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_file_parser,    only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,   only : forcing
use MOM_grid,           only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_intrinsic_functions, only : cuberoot
use MOM_string_functions, only : uppercase
use MOM_unit_scaling,   only : unit_scale_type
use MOM_variables,      only : thermo_var_ptrs, vertvisc_type
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS, Get_Langmuir_Number
use MOM_stochastics,    only : stochastic_CS

implicit none ; private

#include <MOM_memory.h>

public energetic_PBL, energetic_PBL_init, energetic_PBL_end
public energetic_PBL_get_MLD

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> This control structure holds parameters for the MOM_energetic_PBL module
type, public :: energetic_PBL_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.

  !/ Constants
  real    :: VonKar          !< The von Karman coefficient as used in the ePBL module [nondim]
  real    :: omega           !< The Earth's rotation rate [T-1 ~> s-1].
  real    :: omega_frac      !< When setting the decay scale for turbulence, use this fraction of
                             !! the absolute rotation rate blended with the local value of f, as
                             !! sqrt((1-omega_frac)*f^2 + omega_frac*4*omega^2) [nondim].

  !/ Convection related terms
  real    :: nstar           !< The fraction of the TKE input to the mixed layer available to drive
                             !! entrainment [nondim]. This quantity is the vertically integrated
                             !! buoyancy production minus the vertically integrated dissipation of
                             !! TKE produced by buoyancy.

  !/ Mixing Length terms
  logical :: Use_MLD_iteration !< If true, use the proximity to the bottom of the actively turbulent
                             !! surface boundary layer to constrain the mixing lengths.
  logical :: MLD_iteration_guess !< False to default to guessing half the
                             !! ocean depth for the first iteration.
  logical :: MLD_bisection   !< If true, use bisection with the iterative determination of the
                             !! self-consistent mixed layer depth.  Otherwise use the false position
                             !! after a maximum and minimum bound have been evaluated and the
                             !! returned value from the previous guess or bisection before this.
  logical :: MLD_iter_bug    !< If true use buggy logic that gives the wrong bounds for the next
                             !! iteration when successive guesses increase by exactly EPBL_MLD_TOLERANCE.
  integer :: max_MLD_its     !< The maximum number of iterations that can be used to find a
                             !! self-consistent mixed layer depth with Use_MLD_iteration.
  real    :: MixLenExponent  !< Exponent in the mixing length shape-function [nondim].
                             !! 1 is law-of-the-wall at top and bottom,
                             !! 2 is more KPP like.
  real    :: MKE_to_TKE_effic !< The efficiency with which mean kinetic energy released by
                             !!  mechanically forced entrainment of the mixed layer is converted to
                             !!  TKE, times conversion factors between the natural units of mean
                             !!  kinetic energy and those used for TKE [Z2 L-2 ~> nondim].
  logical :: direct_calc     !< If true and there is no conversion from mean kinetic energy to ePBL
                             !! turbulent kinetic energy, use a direct calculation of the
                             !! diffusivity that is supported by a given energy input instead of the
                             !! more general but slower iterative solver.
  real    :: ustar_min       !< A minimum value of ustar to avoid numerical problems [Z T-1 ~> m s-1].
                             !! If the value is small enough, this should not affect the solution.
  real    :: Ekman_scale_coef !< A nondimensional scaling factor controlling the inhibition of the
                             !! diffusive length scale by rotation [nondim].  Making this larger decreases
                             !! the diffusivity in the planetary boundary layer.
  real    :: transLay_scale  !< A scale for the mixing length in the transition layer
                             !! at the edge of the boundary layer as a fraction of the
                             !! boundary layer thickness [nondim].  The default is 0, but a
                             !! value of 0.1 might be better justified by observations.
  real    :: MLD_tol         !< A tolerance for determining the boundary layer thickness when
                             !! Use_MLD_iteration is true [Z ~> m].
  real    :: min_mix_len     !< The minimum mixing length scale that will be used by ePBL [Z ~> m].
                             !! The default (0) does not set a minimum.

  !/ Velocity scale terms
  integer :: wT_scheme       !< An enumerated value indicating the method for finding the turbulent
                             !! velocity scale.  There are currently two options:
                             !! wT_mwT_from_cRoot_TKE is the original (TKE_remaining)^1/3
                             !! wT_from_RH18 is the version described by Reichl and Hallberg, 2018
  real    :: wstar_ustar_coef !< A ratio relating the efficiency with which convectively released
                             !! energy is converted to a turbulent velocity, relative to
                             !! mechanically forced turbulent kinetic energy [nondim].
                             !! Making this larger increases the diffusivity.
  real    :: vstar_surf_fac  !< If (wT_scheme == wT_from_RH18) this is the proportionality coefficient between
                             !! ustar and the surface mechanical contribution to vstar [nondim]
  real    :: vstar_scale_fac !< An overall nondimensional scaling factor for vstar [nondim].  Making
                             !! this larger increases the diffusivity.

  !mstar related options
  integer :: mstar_scheme    !< An encoded integer to determine which formula is used to set mstar
  integer :: BBL_mstar_scheme !< An encoded integer to determine which formula is used to set mstar
  real    :: mstar_cap       !< Since mstar is restoring undissipated energy to mixing,
                             !! there must be a cap on how large it can be [nondim].  This
                             !! is definitely a function of latitude (Ekman limit),
                             !! but will be taken as constant for now.

  !/ vertical decay related options
  real    :: TKE_decay       !< The ratio of the natural Ekman depth to the TKE decay scale [nondim].

  !/ mstar_scheme == 0
  real    :: fixed_mstar     !< mstar is the ratio of the friction velocity cubed to the TKE available to
                             !! drive entrainment [nondim]. This quantity is the vertically
                             !! integrated shear production minus the vertically integrated
                             !! dissipation of TKE produced by shear.  This value is used if the option
                             !! for using a fixed mstar is used.
  real    :: BBL_fixed_mstar !< Similar to fixed_mstar, but for the bottom boundary layer

  !/ mstar_scheme == 2
  real :: C_Ek = 0.17        !< mstar Coefficient in rotation limit for EPBL_MSTAR_SCHEME=OM4 [nondim]
  real :: mstar_coef = 0.3   !< mstar coefficient in rotation/stabilizing balance for EPBL_MSTAR_SCHEME=OM4 [nondim]

  !/ mstar_scheme == 3
  real    :: RH18_mstar_cN1  !< mstar_N coefficient 1 (outer-most coefficient for fit) [nondim].
                             !! Value of 0.275 in RH18.  Increasing this
                             !! coefficient increases mechanical mixing for all values of Hf/ust,
                             !! but is most effective at low values (weakly developed OSBLs).
  real    :: RH18_mstar_cN2  !< mstar_N coefficient 2 (coefficient outside of exponential decay) [nondim].
                             !! Value of 8.0 in RH18.  Increasing this coefficient increases mstar
                             !! for all values of HF/ust, with a consistent affect across
                             !! a wide range of Hf/ust.
  real    :: RH18_mstar_cN3  !< mstar_N coefficient 3 (exponential decay coefficient) [nondim]. Value of
                             !! -5.0 in RH18.  Increasing this increases how quickly the value
                             !! of mstar decreases as Hf/ust increases.
  real    :: RH18_mstar_cS1  !< mstar_S coefficient for RH18 in stabilizing limit [nondim].
                             !! Value of 0.2 in RH18.
  real    :: RH18_mstar_cS2  !< mstar_S exponent for RH18 in stabilizing limit [nondim].
                             !! Value of 0.4 in RH18.

  !/ Coefficient for shear/convective turbulence interaction
  real :: mstar_convect_coef !< Factor to reduce mstar when statically unstable [nondim].

  !/ Langmuir turbulence related parameters
  logical :: Use_LT = .false. !< Flag for using LT in Energy calculation
  integer :: LT_enhance_form !< Integer for Enhancement functional form (various options)
  real    :: LT_enhance_coef !< Coefficient in fit for Langmuir Enhancement [nondim]
  real    :: LT_enhance_exp  !< Exponent in fit for Langmuir Enhancement [nondim]
  real :: LaC_MLD_Ek         !< Coefficient for Langmuir number modification based on the ratio of
                             !! the mixed layer depth over the Ekman depth [nondim].
  real :: LaC_MLD_Ob_stab    !< Coefficient for Langmuir number modification based on the ratio of
                             !! the mixed layer depth over the Obukhov depth with stabilizing forcing [nondim].
  real :: LaC_Ek_Ob_stab     !< Coefficient for Langmuir number modification based on the ratio of
                             !! the Ekman depth over the Obukhov depth with stabilizing forcing [nondim].
  real :: LaC_MLD_Ob_un      !< Coefficient for Langmuir number modification based on the ratio of
                             !! the mixed layer depth over the Obukhov depth with destabilizing forcing [nondim].
  real :: LaC_Ek_Ob_un       !< Coefficient for Langmuir number modification based on the ratio of
                             !! the Ekman depth over the Obukhov depth with destabilizing forcing [nondim].
  real :: Max_Enhance_M = 5. !< The maximum allowed LT enhancement to the mixing [nondim].

  !/ Machine learned equation discovery model paramters
  logical :: eqdisc       !< Uses machine learned shape function
  logical :: eqdisc_v0    !< Uses machine learned velocity scale
  logical :: eqdisc_v0h   !< Uses machine learned velocity scale that uses boundary layer depth as input
  real :: v0_lower_cap    !< Lower cap to prevent v0 from attaining anomlously low values [Z T-1 ~> m s-1]
  real :: v0_upper_cap    !< Upper cap to prevent v0 from attaining anomlously high values [Z T-1 ~> m s-1]
  real :: f_lower !< Lower cap of |f| i.e. absolute of Coriolis parameter [T-1 ~> s-1]
                  !! Used only in get_eqdisc_v0 subroutine. Default is 0.1deg Lat
  real :: bflux_lower_cap !< Lower cap for capping blfux [Z2 T-3 ~> m2 s-3]
  real :: bflux_upper_cap !< Upper cap for capping blfux [Z2 T-3 ~> m2 s-3]
  real :: sigma_max_lower_cap    !< Lower cap to prevent sigma_max from attaining low values [nondim]
  real :: sigma_max_upper_cap    !< Upper cap to prevent sigma_max from attaining high values [nondim]
  real :: Eh_upper_cap !< Upper cap to prevent Eh = hf/(u__*) from attaining high values [nondim]
  real :: Lh_cap       !< Cap to prevent Lh = h/Monin_Obukhov_depth from attaining beyond extreme values [nondim]
  real, allocatable, dimension(:) :: shape_function !< shape function used in machine learned diffusivity [nondim]
  !/ Coefficients used for Machine learned diffusivity
  real :: ML_c(18) !< Array of non-dimensional constants used in machine learned (ML) diffusivity [nondim]
  real :: shape_function_epsilon !< An small value of shape_function below the boundary layer depth [nondim]

  !/ Bottom boundary layer mixing related options
  real :: ePBL_BBL_effic     !< The efficiency of bottom boundary layer mixing via ePBL driven by
                             !! the bottom drag dissipation of mean kinetic energy, times
                             !! conversion factors between the natural units of mean kinetic energy
                             !! and those used for TKE [Z2 L-2 ~> nondim].
  real :: ePBL_tidal_effic   !< The efficiency of bottom boundary layer mixing via ePBL driven by
                             !! the bottom drag dissipation of tides, times conversion factors
                             !! between the natural units of mean kinetic energy and those used for
                             !! TKE [Z2 L-2 ~> nondim].
  logical :: Use_BBLD_iteration !< If true, use the proximity to the top of the actively turbulent
                             !! bottom boundary layer to constrain the mixing lengths.
  real    :: TKE_decay_BBL   !< The ratio of the natural Ekman depth to the TKE decay scale for
                             !! bottom boundary layer mixing [nondim]
  real    :: min_BBL_mix_len !< The minimum mixing length scale that will be used by ePBL in the bottom
                             !! boundary layer mixing [Z ~> m].  The default (0) does not set a minimum.
  real    :: MixLenExponent_BBL !< Exponent in the bottom boundary layer mixing length shape-function [nondim].
                             !! 1 is law-of-the-wall at top and bottom,
                             !! 2 is more KPP like.
  real    :: BBLD_tol        !< The tolerance for the iteratively determined bottom boundary layer depth [Z ~> m].
                             !! This is only used with USE_MLD_ITERATION.
  integer :: max_BBLD_its    !< The maximum number of iterations that can be used to find a self-consistent
                             !! bottom boundary layer depth.
  integer :: wT_scheme_BBL   !< An enumerated value indicating the method for finding the bottom boundary
                             !! layer turbulent velocity scale.  There are currently two options:
                             !! wT_mwT_from_cRoot_TKE is the original (TKE_remaining)^1/3
                             !! wT_from_RH18 is the version described by Reichl and Hallberg, 2018
  real :: vstar_scale_fac_BBL !< An overall nondimensional scaling factor for wT in the bottom boundary layer [nondim].
                             !! Making this larger increases the bottom boundary layer diffusivity.", &
  real :: vstar_surf_fac_BBL !< If (wT_scheme_BBL == wT_from_RH18) this is the proportionality coefficient between
                             !! ustar and the bottom boundayer layer mechanical contribution to vstar [nondim]
  real :: Ekman_scale_coef_BBL !< A nondimensional scaling factor controlling the inhibition of the
                             !! diffusive length scale by rotation in the bottom boundary layer [nondim].
                             !! Making this larger decreases the bottom boundary layer diffusivity.
  logical :: decay_adjusted_BBL_TKE !< If true, include an adjustment factor in the bottom boundary layer
                             !! energetics that accounts for an exponential decay of TKE from a
                             !! near-bottom source and an assumed piecewise linear linear profile
                             !! of the buoyancy flux response to a change in a diffusivity.
  logical :: BBL_effic_bug   !< If true, overestimate the efficiency of the non-tidal ePBL bottom boundary
                             !! layer diffusivity by a factor of 1/sqrt(CDRAG), which is often a factor of
                             !! about 18.3.
  logical :: ePBL_BBL_use_mstar !< If true, use an mstar*ustar^3 paramaterization to get the TKE available
                             !! to drive mixing in the bottom boundary layer version of ePBL.  Otherwise,
                             !! use the meanflow energy loss to bottom drag scaled by a constant efficiency.

  !/ Options for documenting differences from parameter choices
  integer :: options_diff    !< If positive, this is a coded integer indicating a pair of
                             !! settings whose differences are diagnosed in a passive diagnostic mode
                             !! via extra calls to ePBL_column.  If this is 0 or negative no extra
                             !! calls occur.

  !/ Others
  type(time_type), pointer :: Time=>NULL() !< A pointer to the ocean model's clock.

  logical :: TKE_diagnostics = .false. !< If true, diagnostics of the TKE budget are being calculated.
  integer :: answer_date     !< The vintage of the order of arithmetic and expressions in the ePBL
                             !! calculations.  Values below 20190101 recover the answers from the
                             !! end of 2018, while higher values use updated and more robust forms
                             !! of the same expressions.  Values below 20240101 use A**(1./3.) to
                             !! estimate the cube root of A in several expressions, while higher
                             !! values use the integer root function cuberoot(A) and therefore
                             !! can work with scaled variables.
  logical :: orig_PE_calc    !< If true, the ePBL code uses the original form of the
                             !! potential energy change code.  Otherwise, it uses a newer version
                             !! that can work with successive increments to the diffusivity in
                             !! upward or downward passes.
  logical :: debug           !< If true, write verbose checksums for debugging purposes.
  type(diag_ctrl), pointer :: diag=>NULL() !< A structure that is used to regulate the
                             !! timing of diagnostic output.

  real, allocatable, dimension(:,:) :: &
    ML_depth                 !< The mixed layer depth determined by active mixing in ePBL, which may
                             !! be used for the first guess in the next time step [H ~> m or kg m-2]
  real, allocatable, dimension(:,:) :: &
    BBL_depth                !< The bottom boundary layer depth determined by active mixing in ePBL [H ~> m or kg m-2]

  type(EFP_type), dimension(2) :: sum_its !< The total number of iterations and columns worked on
  type(EFP_type), dimension(2) :: sum_its_BBL !< The total number of iterations and columns worked on

  !>@{ Diagnostic IDs
  integer :: id_ML_depth = -1, id_hML_depth = -1, id_TKE_wind = -1, id_TKE_mixing = -1
  integer :: id_ustar_ePBL = -1, id_bflx_ePBL = -1
  integer :: id_TKE_MKE = -1, id_TKE_conv = -1, id_TKE_forcing = -1
  integer :: id_TKE_mech_decay = -1, id_TKE_conv_decay = -1
  integer :: id_Mixing_Length = -1, id_Velocity_Scale = -1
  integer :: id_Kd_BBL = -1, id_BBL_Mix_Length = -1, id_BBL_Vel_Scale = -1
  integer :: id_TKE_BBL = -1, id_TKE_BBL_mixing = -1, id_TKE_BBL_decay = -1
  integer :: id_ustar_BBL = -1, id_bflx_BBL = -1, id_BBL_decay_scale = -1, id_BBL_depth = -1
  integer :: id_mstar_sfc = -1, id_mstar_BBL = -1, id_LA_mod = -1, id_LA = -1, id_mstar_LT = -1
  ! The next options are used when passively diagnosing sensitivities from parameter choices
  integer :: id_opt_diff_Kd_ePBL = -1, id_opt_maxdiff_Kd_ePBL = -1, id_opt_diff_hML_depth = -1
  !>@}
end type energetic_PBL_CS

!>@{ Enumeration values for mstar_scheme
integer, parameter :: Use_Fixed_mstar = 0  !< The value of mstar_scheme to use a constant mstar
integer, parameter :: mstar_from_Ekman = 2 !< The value of mstar_scheme to base mstar on the ratio
                                           !! of the Ekman layer depth to the Obukhov depth
integer, parameter :: mstar_from_RH18 = 3  !< The value of mstar_scheme to base mstar of of RH18
integer, parameter :: No_Langmuir = 0      !< The value of LT_enhance_form not use Langmuir turbulence.
integer, parameter :: Langmuir_rescale = 2 !< The value of LT_enhance_form to use a multiplicative
                                           !! rescaling of mstar to account for Langmuir turbulence.
integer, parameter :: Langmuir_add = 3     !< The value of LT_enhance_form to add a contribution to
                                           !! mstar from Langmuir turbulence to other contributions.
integer, parameter :: wT_from_cRoot_TKE = 0 !< Use a constant times the cube root of remaining TKE
                                           !! to calculate the turbulent velocity.
integer, parameter :: wT_from_RH18 = 1     !< Use a scheme based on a combination of w* and v* as
                                           !! documented in Reichl & Hallberg (2018) to calculate
                                           !! the turbulent velocity.
character*(20), parameter :: CONSTANT_STRING = "CONSTANT"
character*(20), parameter :: OM4_STRING = "OM4"
character*(20), parameter :: RH18_STRING = "REICHL_H18"
character*(20), parameter :: ROOT_TKE_STRING = "CUBE_ROOT_TKE"
character*(20), parameter :: NONE_STRING = "NONE"
character*(20), parameter :: RESCALED_STRING = "RESCALE"
character*(20), parameter :: ADDITIVE_STRING = "ADDITIVE"
!>@}

logical :: report_avg_its = .false.  !< Report the average number of ePBL iterations for debugging.

!> A type for conveniently passing around ePBL diagnostics for a column.
type, public :: ePBL_column_diags ; private
  !>@{ Local column copies of energy change diagnostics, all in [R Z3 T-3 ~> W m-2].
  real :: dTKE_conv, dTKE_forcing, dTKE_wind, dTKE_mixing ! Local column diagnostics [R Z3 T-3 ~> W m-2]
  real :: dTKE_MKE, dTKE_mech_decay, dTKE_conv_decay      ! Local column diagnostics [R Z3 T-3 ~> W m-2]
  real :: dTKE_BBL, dTKE_BBL_decay, dTKE_BBL_mixing       ! Local column diagnostics [R Z3 T-3 ~> W m-2]
  !>@}
  real :: LA        !< The value of the Langmuir number [nondim]
  real :: LAmod     !< The modified Langmuir number by convection [nondim]
  real :: mstar     !< The value of mstar used in ePBL [nondim]
  real :: mstar_BBL !< The value of mstar used in ePBL BBL [nondim]
  real :: mstar_LT  !< The portion of mstar due to Langmuir turbulence [nondim]
  integer :: OBL_its !< The number of iterations used to find a self-consistent surface boundary layer depth
  integer :: BBL_its !< The number of iterations used to find a self-consistent bottom boundary layer depth
end type ePBL_column_diags

contains

!>    This subroutine determines the diffusivities from the integrated energetics
!!  mixed layer model.  It assumes that heating, cooling and freshwater fluxes
!!  have already been applied.  All calculations are done implicitly, and there
!!  is no stability limit on the time step.
subroutine energetic_PBL(h_3d, u_3d, v_3d, tv, fluxes, visc, dt, Kd_int, G, GV, US, CS, &
                         stoch_CS, dSV_dT, dSV_dS, TKE_forced, buoy_flux, BBL_buoy_flux, Waves )
  type(ocean_grid_type),   intent(inout) :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h_3d   !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u_3d   !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: v_3d   !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dSV_dT !< The partial derivative of in-situ specific
                                                   !! volume with potential temperature
                                                   !! [R-1 C-1 ~> m3 kg-1 degC-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dSV_dS !< The partial derivative of in-situ specific
                                                   !! volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: TKE_forced !< The forcing requirements to homogenize the
                                                   !! forcing that has been applied to each layer
                                                   !! [R Z3 T-2 ~> J m-2].
  type(thermo_var_ptrs),   intent(inout) :: tv     !< A structure containing pointers to any
                                                   !! available thermodynamic fields. Absent fields
                                                   !! have NULL ptrs.
  type(forcing),           intent(inout) :: fluxes !< A structure containing pointers to any
                                                   !! possible forcing fields. Unused fields have
                                                   !! NULL ptrs.
  type(vertvisc_type),     intent(in)    :: visc   !< Structure with vertical viscosities,
                                                   !! BBL properties and related fields
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(out)   :: Kd_int !< The diagnosed diffusivities at interfaces
                                                   !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  type(energetic_PBL_CS),  intent(inout) :: CS     !< Energetic PBL control structure
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: buoy_flux !< The surface buoyancy flux [Z2 T-3 ~> m2 s-3].
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: BBL_buoy_flux !< The bottom buoyancy flux [Z2 T-3 ~> m2 s-3].
  type(wave_parameters_CS), pointer      :: Waves  !< Waves control structure for Langmuir turbulence
  type(stochastic_CS),     pointer       :: stoch_CS  !< The control structure returned by a previous

!    This subroutine determines the diffusivities from the integrated energetics
!  mixed layer model.  It assumes that heating, cooling and freshwater fluxes
!  have already been applied.  All calculations are done implicitly, and there
!  is no stability limit on the time step.
!
!    For each interior interface, first discard the TKE to account for mixing
! of shortwave radiation through the next denser cell.  Next drive mixing based
! on the local? values of ustar + wstar, subject to available energy.  This
! step sets the value of Kd(K).  Any remaining energy is then subject to decay
! before being handed off to the next interface.  mech_TKE and conv_PErel are treated
! separately for the purposes of decay, but are used proportionately to drive
! mixing.
!
!   The key parameters for the mixed layer are found in the control structure.
!   To use the classic constant mstar mixed layers choose EPBL_MSTAR_SCHEME=CONSTANT.
! The key parameters then include mstar, nstar, TKE_decay, and conv_decay.
! For the Oberhuber (1993) mixed layer,the values of these are:
!      mstar = 1.25,  nstar = 1, TKE_decay = 2.5, conv_decay = 0.5
! TKE_decay is 1/kappa in eq. 28 of Oberhuber (1993), while conv_decay is 1/mu.
! For a traditional Kraus-Turner mixed layer, the values are:
!      mstar = 1.25, nstar = 0.4, TKE_decay = 0.0, conv_decay = 0.0

  ! Local variables
  real, dimension(SZI_(G),SZK_(GV)) :: &
    h_2d, &         ! A 2-d slice of the layer thickness [H ~> m or kg m-2].
    dz_2d, &        ! A 2-d slice of the vertical distance across layers [Z ~> m].
    T_2d, &         ! A 2-d slice of the layer temperatures [C ~> degC].
    S_2d, &         ! A 2-d slice of the layer salinities [S ~> ppt].
    TKE_forced_2d, & ! A 2-d slice of TKE_forced [R Z3 T-2 ~> J m-2].
    dSV_dT_2d, &    ! A 2-d slice of dSV_dT [R-1 C-1 ~> m3 kg-1 degC-1].
    dSV_dS_2d, &    ! A 2-d slice of dSV_dS [R-1 S-1 ~> m3 kg-1 ppt-1].
    u_2d, &         ! A 2-d slice of the zonal velocity [L T-1 ~> m s-1].
    v_2d            ! A 2-d slice of the meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK_(GV)+1) :: &
    Kd_2d           ! A 2-d version of the diapycnal diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZK_(GV)) :: &
    h, &            ! The layer thickness [H ~> m or kg m-2].
    dz, &           ! The vertical distance across layers [Z ~> m].
    T0, &           ! The initial layer temperatures [C ~> degC].
    S0, &           ! The initial layer salinities [S ~> ppt].
    dSV_dT_1d, &    ! The partial derivatives of specific volume with temperature [R-1 C-1 ~> m3 kg-1 degC-1].
    dSV_dS_1d, &    ! The partial derivatives of specific volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1].
    TKE_forcing, &  ! Forcing of the TKE in the layer coming from TKE_forced [R Z3 T-2 ~> J m-2].
    u, &            ! The zonal velocity [L T-1 ~> m s-1].
    v               ! The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZK_(GV)+1) :: &
    Kd, &           ! The diapycnal diffusivity due to ePBL [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
    mixvel, &       ! A turbulent mixing velocity [Z T-1 ~> m s-1].
    mixlen, &       ! A turbulent mixing length [Z ~> m].
    mixvel_BBL, &   ! A bottom boundary layer turbulent mixing velocity [Z T-1 ~> m s-1].
    mixlen_BBL, &   ! A bottom boundary layer turbulent mixing length [Z ~> m].
    Kd_BBL, &       ! The bottom boundary layer diapycnal diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
    SpV_dt, &       ! Specific volume interpolated to interfaces divided by dt or 1.0 / (dt * Rho0),
                    ! in [R-1 T-1 ~> m3 kg-1 s-1], used to convert local TKE into a turbulence velocity cubed.
    SpV_dt_cf       ! Specific volume interpolated to interfaces divided by dt or 1.0 / (dt * Rho0)
                    ! times conversion factors for answer dates before 20240101 in
                    ! [m3 Z-3 R-1 T2 s-3 ~> m3 kg-1 s-1] or without the conversion factors for
                    ! answer dates of 20240101 and later in [R-1 T-1 ~> m3 kg-1 s-1], used to
                    ! convert local TKE into a turbulence velocity cubed.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].

  real :: absf      ! The absolute value of f [T-1 ~> s-1].
  real :: U_star    ! The surface friction velocity [Z T-1 ~> m s-1].
  real :: U_Star_Mean ! The surface friction without gustiness [Z T-1 ~> m s-1].
  real :: mech_TKE  ! The mechanically generated turbulent kinetic energy available for mixing over a
                    ! timestep before the application of the efficiency in mstar [R Z3 T-2 ~> J m-2]
  real :: u_star_BBL ! The bottom boundary layer friction velocity [H T-1 ~> m s-1 or kg m-2 s-1].
  real :: u_star_BBL_z_t ! The bottom boundary layer friction velocity converted to Z T-1 [Z T-1 ~> m s-1].
  real :: BBL_TKE   ! The mechanically generated turbulent kinetic energy available for bottom
                    ! boundary layer mixing within a timestep [R Z3 T-2 ~> J m-2]
  real :: I_rho     ! The inverse of the Boussinesq reference density [R-1 ~> m3 kg-1]
  real :: I_dt      ! The Adcroft reciprocal of the timestep [T-1 ~> s-1]
  real :: I_rho0dt  ! The inverse of the Boussinesq reference density times the time
                    ! step [R-1 T-1 ~> m3 kg-1 s-1]
  real :: B_Flux    ! The surface buoyancy flux [Z2 T-3 ~> m2 s-3]
  real :: MLD_io    ! The mixed layer depth found by ePBL_column [Z ~> m]
  real :: BBLD_io   ! The bottom boundary layer thickness found by ePBL_BBL_column [Z ~> m]
  real :: MLD_in    ! The first guess at the mixed layer depth [Z ~> m]
  real :: BBLD_in   ! The first guess at the bottom boundary layer thickness [Z ~> m]

  type(ePBL_column_diags) :: eCD ! A container for passing around diagnostics.

  ! The following variables are used for diagnostics
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    diag_Velocity_Scale, & ! The velocity scale used in getting Kd [Z T-1 ~> m s-1]
    diag_Mixing_Length, &  ! The length scale used in getting Kd [Z ~> m]
    Kd_BBL_3d, &           ! The bottom boundary layer diffusivities [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
    BBL_Vel_Scale, &       ! The velocity scale used in getting the BBL part of Kd [Z T-1 ~> m s-1]
    BBL_Mix_Length         ! The length scale used in getting the BBL part of Kd [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    ! The next 7 diagnostics are terms in the mixed layer TKE budget, all in [R Z3 T-3 ~> W m-2 = kg s-3].
    diag_TKE_wind, &   ! The wind source of TKE [R Z3 T-3 ~> W m-2]
    diag_TKE_MKE, &    ! The resolved KE source of TKE [R Z3 T-3 ~> W m-2]
    diag_TKE_conv, &   ! The convective source of TKE [R Z3 T-3 ~> W m-2]
    diag_TKE_forcing, & ! The TKE sink required to mix surface penetrating shortwave heating [R Z3 T-3 ~> W m-2]
    diag_TKE_mech_decay, & ! The decay of mechanical TKE [R Z3 T-3 ~> W m-2]
    diag_TKE_conv_decay, & ! The decay of convective TKE [R Z3 T-3 ~> W m-2]
    diag_TKE_mixing, & ! The work done by TKE to deepen the mixed layer [R Z3 T-3 ~> W m-2]
    diag_TKE_BBL, &    ! The source of TKE to the bottom boundary layer [R Z3 T-3 ~> W m-2].
    diag_TKE_BBL_mixing, & ! The work done by TKE to thicken the bottom boundary layer [R Z3 T-3 ~> W m-2].
    diag_TKE_BBL_decay, & ! The work lost to decy of mechanical TKE in the bottom boundary
                       ! layer [R Z3 T-3 ~> W m-2].
    diag_ustar_BBL, &  ! The bottom boundary layer friction velocity [H T-1 ~> m s-1 or kg m-2 s-1]
    diag_BBL_decay_scale, & ! The bottom boundary layer TKE decay length scale [H ~> m]
    diag_mstar_sfc, &  ! mstar used in EPBL [nondim]
    diag_mstar_BBL, &  ! mstar used in EPBL BBL [nondim]
    diag_mstar_LT, &   ! mstar due to Langmuir turbulence [nondim]
    diag_LA, &         ! Langmuir number [nondim]
    diag_LA_mod, &     ! Modified Langmuir number [nondim]
    diag_ustar, &      ! The surface boundary layer friction velocity [Z T-1 ~> m s-1]
    diag_bflx          ! The surface boundary layer buoyancy flux  [Z2 T-3 ~> m2 s-3]

  ! The following variables are only used for diagnosing sensitivities to ePBL settings
  real, dimension(SZK_(GV)+1) :: &
    Kd_1, Kd_2      ! Diapycnal diffusivities found with different ePBL options [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real :: diff_Kd(SZI_(G),SZJ_(G),SZK_(GV)+1) ! The change in diapycnal diffusivities found with different
                        ! ePBL options [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real :: max_abs_diff_Kd(SZI_(G),SZJ_(G))  ! The column maximum magnitude of the change in diapycnal
                        ! diffusivities found with different ePBL options [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real :: diff_hML_depth(SZI_(G),SZJ_(G))   ! The change in diagnosed active mixing layer depth with
                        ! different ePBL options [Z ~> m]
  real :: BLD_1, BLD_2  ! Surface or bottom boundary layer depths found with different ePBL_column options [Z ~> m]
  real :: SpV_scale1    ! A factor that accounts for the varying scaling of SpV_dt with answer date
                        ! [nondim] or [T3 m3 Z-3 s-3 ~> 1]
  real :: SpV_scale2    ! A factor that accounts for the varying scaling of SpV_dt with answer date
                        ! [nondim] or [Z3 s3 T-3 m-3 ~> 1]
  real :: SpV_dt_tmp(SZK_(GV)+1)  ! Specific volume interpolated to interfaces divided by dt or 1.0 / (dt * Rho0)
                        ! times conversion factors for answer dates before 20240101 in
                        ! [m3 Z-3 R-1 T2 s-3 ~> m3 kg-1 s-1] or without the conversion factors for
                        ! answer dates of 20240101 and later in [R-1 T-1 ~> m3 kg-1 s-1], used to
                        ! convert local TKE into a turbulence velocity cubed.
  type(ePBL_column_diags) :: eCD_tmp   ! A container for not passing around diagnostics.
  type(energetic_PBL_CS)  :: CS_tmp1, CS_tmp2 ! Copies of the energetic PBL control structure that
                                       ! can be modified to test for sensitivities
  logical :: BBL_mixing ! If true, there is bottom boundary layer mixing.
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not. CS%initialized) call MOM_error(FATAL, "energetic_PBL: "//&
         "Module must be initialized before it is used.")
  if (.not. associated(tv%eqn_of_state)) call MOM_error(FATAL, &
      "energetic_PBL: Temperature, salinity and an equation of state "//&
      "must now be used.")
  if (.not.(associated(fluxes%ustar) .or. associated(fluxes%tau_mag))) call MOM_error(FATAL, &
      "energetic_PBL: No surface friction velocity (ustar or tau_mag) defined in fluxes type.")
  if ((.not.GV%Boussinesq) .and. (.not.associated(fluxes%tau_mag))) call MOM_error(FATAL, &
      "energetic_PBL: No surface wind stress magnitude defined in fluxes type in non-Boussinesq mode.")
  if (CS%use_LT .and. .not.associated(Waves)) call MOM_error(FATAL, &
      "energetic_PBL: The Waves control structure must be associated if CS%use_LT "//&
      "(i.e., USE_LA_LI2016 or EPBL_LT) is True.")


  h_neglect = GV%H_subroundoff
  I_rho = GV%H_to_Z * GV%RZ_to_H ! == 1.0 / GV%Rho0 ! This is not used when fully non-Boussinesq.
  I_dt = 0.0 ; if (dt > 0.0) I_dt = 1.0 / dt
  I_rho0dt = 1.0 / (GV%Rho0 * dt)  ! This is not used when fully non-Boussinesq.
  BBL_mixing = ((CS%ePBL_BBL_effic > 0.0) .or. (CS%ePBL_tidal_effic > 0.0) .or. CS%ePBL_BBL_use_mstar)

  ! Zero out diagnostics before accumulation.
  if (CS%TKE_diagnostics) then
    !!OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      diag_TKE_wind(i,j) = 0.0 ; diag_TKE_MKE(i,j) = 0.0
      diag_TKE_conv(i,j) = 0.0 ; diag_TKE_forcing(i,j) = 0.0
      diag_TKE_mixing(i,j) = 0.0 ; diag_TKE_mech_decay(i,j) = 0.0
      diag_TKE_conv_decay(i,j) = 0.0 !; diag_TKE_unbalanced(i,j) = 0.0
    enddo ; enddo
    if (BBL_mixing) then
      !!OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        diag_TKE_BBL(i,j) = 0.0 ; diag_TKE_BBL_mixing(i,j) = 0.0
        diag_TKE_BBL_decay(i,j) = 0.0
      enddo ; enddo
    endif
  endif
  if (CS%debug .or. (CS%id_Mixing_Length>0)) diag_Mixing_Length(:,:,:) = 0.0
  if (CS%debug .or. (CS%id_Velocity_Scale>0)) diag_Velocity_Scale(:,:,:) = 0.0
  if (BBL_mixing) then
    if (CS%debug .or. (CS%id_BBL_Mix_Length>0)) BBL_Mix_Length(:,:,:) = 0.0
    if (CS%debug .or. (CS%id_BBL_Vel_Scale>0)) BBL_Vel_Scale(:,:,:) = 0.0
    if (CS%id_Kd_BBL > 0) Kd_BBL_3d(:,:,:) = 0.0
    if (CS%id_ustar_BBL > 0) diag_ustar_BBL(:,:) = 0.0
    if (CS%id_BBL_decay_scale > 0) diag_BBL_decay_scale(:,:) = 0.0
  endif

  ! CS_tmp is used to test sensitivity to parameter setting changes.
  if (CS%options_diff > 0) then
    CS_tmp1 = CS ; CS_tmp2 = CS
    SpV_scale1 = 1.0 ; SpV_scale2 = 1.0

    if (CS%options_diff == 1) then
      CS_tmp1%orig_PE_calc = .true. ; CS_tmp2%orig_PE_calc = .false.
    elseif (CS%options_diff == 2) then
      CS_tmp1%answer_date = 20181231 ; CS_tmp2%answer_date = 20240101
    elseif (CS%options_diff == 3) then
      CS_tmp1%direct_calc = .true.   ; CS_tmp2%direct_calc = .false.
      CS_tmp1%MKE_to_TKE_effic = 0.0 ; CS_tmp2%MKE_to_TKE_effic = 0.0
      CS_tmp1%orig_PE_calc = .false. ; CS_tmp2%orig_PE_calc = .false.
    elseif (CS%options_diff == 4) then
      CS_tmp1%direct_calc = .true.   ; CS_tmp2%direct_calc = .false.
      CS_tmp1%MKE_to_TKE_effic = 0.0 ; CS_tmp2%MKE_to_TKE_effic = 0.0
      CS_tmp1%ePBL_BBL_effic = 0.2   ; CS_tmp2%ePBL_BBL_effic = 0.2
    elseif (CS%options_diff == 5) then
      CS_tmp1%decay_adjusted_BBL_TKE = .true. ; CS_tmp2%decay_adjusted_BBL_TKE = .false.
      CS_tmp1%MKE_to_TKE_effic = 0.0 ; CS_tmp2%MKE_to_TKE_effic = 0.0
      CS_tmp1%ePBL_BBL_effic = 0.2   ; CS_tmp2%ePBL_BBL_effic = 0.2
    endif
    ! This logic is needed because the scaling of SpV_dt changes with answer date.
    if (CS_tmp1%answer_date < 20240101) SpV_scale1 = US%m_to_Z**3 * US%T_to_s**3
    if (CS_tmp2%answer_date < 20240101) SpV_scale2 = US%m_to_Z**3 * US%T_to_s**3
    if (CS%id_opt_diff_Kd_ePBL > 0)    diff_Kd(:,:,:) = 0.0
    if (CS%id_opt_maxdiff_Kd_ePBL > 0) max_abs_diff_Kd(:,:) = 0.0
    if (CS%id_opt_diff_hML_depth > 0)  diff_hML_depth(:,:) = 0.0
  endif

  !!OMP parallel do default(private) shared(js,je,nz,is,ie,h_3d,u_3d,v_3d,tv,dt,I_dt,BBL_mixing, &
  !!OMP                                  CS,G,GV,US,fluxes,TKE_forced,dSV_dT,dSV_dS,Kd_int)
  do j=js,je
    ! Copy the thicknesses and other fields to 2-d arrays.
    do k=1,nz ; do i=is,ie
      h_2d(i,k) = h_3d(i,j,k) ; u_2d(i,k) = u_3d(i,j,k) ; v_2d(i,k) = v_3d(i,j,k)
      T_2d(i,k) = tv%T(i,j,k) ; S_2d(i,k) = tv%S(i,j,k)
      TKE_forced_2d(i,k) = TKE_forced(i,j,k)
      dSV_dT_2d(i,k) = dSV_dT(i,j,k) ; dSV_dS_2d(i,k) = dSV_dS(i,j,k)
    enddo ; enddo
    call thickness_to_dz(h_3d, tv, dz_2d, j, G, GV)

    ! Set the inverse density used to translating local TKE into a turbulence velocity
    SpV_dt(:) = 0.0
    if ((dt > 0.0) .and. GV%Boussinesq .or. .not.allocated(tv%SpV_avg)) then
      if (CS%answer_date < 20240101) then
        do K=1,nz+1
          SpV_dt(K) = 1.0 / (dt*GV%Rho0)
        enddo
      else
        do K=1,nz+1
          SpV_dt(K) = I_rho0dt
        enddo
      endif
    endif

    !   Determine the initial mech_TKE and conv_PErel, including the energy required
    ! to mix surface heating through the topmost cell, the energy released by mixing
    ! surface cooling & brine rejection down through the topmost cell, and
    ! homogenizing the shortwave heating within that cell.  This sets the energy
    ! and ustar and wstar available to drive mixing at the first interior
    ! interface.
    do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then

      ! Copy the thicknesses and other fields to 1-d arrays.
      do k=1,nz
        h(k) = h_2d(i,k) + GV%H_subroundoff ; dz(k) = dz_2d(i,k) + GV%dZ_subroundoff
        u(k) = u_2d(i,k) ; v(k) = v_2d(i,k)
        T0(k) = T_2d(i,k) ; S0(k) = S_2d(i,k) ; TKE_forcing(k) =  TKE_forced_2d(i,k)
        dSV_dT_1d(k) = dSV_dT_2d(i,k) ; dSV_dS_1d(k) = dSV_dS_2d(i,k)
      enddo
      do K=1,nz+1 ; Kd(K) = 0.0 ; enddo

      ! Make local copies of surface forcing and process them.
      if (associated(fluxes%ustar) .and. (GV%Boussinesq .or. .not.associated(fluxes%tau_mag))) then
        u_star = fluxes%ustar(i,j)
        u_star_Mean = fluxes%ustar_gustless(i,j)
        mech_TKE = dt * GV%Rho0 * u_star**3
      elseif (allocated(tv%SpV_avg)) then
        u_star = sqrt(fluxes%tau_mag(i,j) * tv%SpV_avg(i,j,1))
        u_star_Mean = sqrt(fluxes%tau_mag_gustless(i,j) * tv%SpV_avg(i,j,1))
        mech_TKE = dt * u_star * fluxes%tau_mag(i,j)
      else
        u_star = sqrt(fluxes%tau_mag(i,j) * I_rho)
        u_star_Mean = sqrt(fluxes%tau_mag_gustless(i,j) * I_rho)
        mech_TKE = dt * GV%Rho0 * u_star**3
        ! The line above is equivalent to: mech_TKE = dt * u_star * fluxes%tau_mag(i,j)
      endif
      diag_ustar(i,j) = u_star

      if (allocated(tv%SpV_avg) .and. .not.GV%Boussinesq) then
        SpV_dt(1) = tv%SpV_avg(i,j,1) * I_dt
        do K=2,nz
          SpV_dt(K) = 0.5*(tv%SpV_avg(i,j,k-1) + tv%SpV_avg(i,j,k)) * I_dt
        enddo
        SpV_dt(nz+1) = tv%SpV_avg(i,j,nz) * I_dt
      endif

      B_flux = buoy_flux(i,j)
      if (associated(fluxes%ustar_shelf) .and. associated(fluxes%frac_shelf_h)) then
        if (fluxes%frac_shelf_h(i,j) > 0.0) &
          u_star = (1.0 - fluxes%frac_shelf_h(i,j)) * u_star + &
                   fluxes%frac_shelf_h(i,j) * fluxes%ustar_shelf(i,j)
      endif
      if (u_star < CS%ustar_min) u_star = CS%ustar_min
      if (CS%omega_frac >= 1.0) then
        absf = 2.0*CS%omega
      else
        absf = 0.25*((abs(G%CoriolisBu(I,J)) + abs(G%CoriolisBu(I-1,J-1))) + &
                     (abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I-1,J))))
        if (CS%omega_frac > 0.0) &
          absf = sqrt(CS%omega_frac*4.0*CS%omega**2 + (1.0-CS%omega_frac)*absf**2)
      endif

      ! Perhaps provide a first guess for MLD based on a stored previous value.
      MLD_io = -1.0
      if (CS%MLD_iteration_guess .and. (CS%ML_depth(i,j) > 0.0))  MLD_io = CS%ML_depth(i,j)
      BBLD_io = 0.0

      ! Store the initial guesses at the boundary layer depths for testing sensitivities.
      MLD_in = MLD_io

      if (CS%answer_date < 20240101) then
        do K=1,nz+1 ; SpV_dt_cf(K) = (US%Z_to_m**3*US%s_to_T**3) * SpV_dt(K) ; enddo
      else
        do K=1,nz+1 ; SpV_dt_cf(K) = SpV_dt(K) ; enddo
      endif
      if (stoch_CS%pert_epbl) then ! stochastics are active
        call ePBL_column(h, dz, u, v, T0, S0, dSV_dT_1d, dSV_dS_1d, SpV_dt_cf, TKE_forcing, B_flux, absf, &
                         u_star, u_star_mean, mech_TKE, dt, MLD_io, Kd, mixvel, mixlen, GV, &
                         US, CS, eCD, Waves, G, i, j, &
                         TKE_gen_stoch=stoch_CS%epbl1_wts(i,j), TKE_diss_stoch=stoch_CS%epbl2_wts(i,j))
      else
        call ePBL_column(h, dz, u, v, T0, S0, dSV_dT_1d, dSV_dS_1d, SpV_dt_cf, TKE_forcing, B_flux, absf, &
                         u_star, u_star_mean, mech_TKE, dt, MLD_io, Kd, mixvel, mixlen, GV, &
                         US, CS, eCD, Waves, G, i, j)
      endif

      ! Add the diffusivity due to bottom boundary layer mixing, if there is energy to drive this mixing.
      if (BBL_mixing) then
        if (CS%MLD_iteration_guess .and. (CS%BBL_depth(i,j) > 0.0)) BBLD_io = CS%BBL_depth(i,j)
        BBLD_in = BBLD_io
        u_star_BBL = max(visc%ustar_BBL(i,j), CS%ustar_min*GV%Z_to_H)  ! units are H T-1
        if (GV%Boussinesq) then
          u_star_BBL_z_t = u_star_BBL*GV%H_to_Z
        else
          u_star_BBL_z_t = u_star_BBL*GV%H_to_RZ*tv%SpV_avg(i,j,1)
        endif

        if (CS%ePBL_BBL_use_mstar) then
          BBL_TKE = dt * ((u_star_BBL*GV%H_to_RZ) * u_star_BBL_z_t**2)
        else
          if (CS%BBL_effic_bug) then
            BBL_TKE = CS%ePBL_BBL_effic * GV%H_to_RZ * dt * visc%BBL_meanKE_loss_sqrtCd(i,j)
          else
            BBL_TKE = CS%ePBL_BBL_effic * GV%H_to_RZ * dt * visc%BBL_meanKE_loss(i,j)
          endif
          ! Add in tidal dissipation energy at the bottom, noting that fluxes%BBL_tidal_dis is
          ! in [R Z L2 T-3 ~> W m-2], unlike visc%BBL_meanKE_loss.
          if ((CS%ePBL_tidal_effic > 0.0) .and. associated(fluxes%BBL_tidal_dis)) &
            BBL_TKE = BBL_TKE + CS%ePBL_tidal_effic * dt * fluxes%BBL_tidal_dis(i,j)
        endif

        call ePBL_BBL_column(h, dz, u, v, T0, S0, dSV_dT_1d, dSV_dS_1d, SpV_dt, absf, dt, Kd, BBL_TKE, &
                             u_star_BBL, u_star_BBL_z_t, BBL_buoy_flux(i,j), Kd_BBL, BBLD_io, mixvel_BBL, mixlen_BBL, &
                             GV, US, CS, eCD)

        do K=1,nz+1 ; Kd(K) = Kd(K) + Kd_BBL(K) ; enddo
        if (CS%id_Kd_BBL > 0) then ; do K=1,nz+1
          Kd_BBL_3d(i,j,K) = Kd_BBL(K)
        enddo ; endif
        if (CS%id_ustar_BBL > 0) diag_ustar_BBL(i,j) = u_star_BBL
        if ((CS%id_BBL_decay_scale > 0) .and. (CS%TKE_decay * absf > 0)) &
          diag_BBL_decay_scale(i,j) = u_star_BBL / (CS%TKE_decay * absf)
      endif

      ! Copy the diffusivities to a 2-d array.
      do K=1,nz+1
        Kd_2d(i,K) = Kd(K)
      enddo
      CS%ML_depth(i,j) = MLD_io
      CS%BBL_depth(i,j) = BBLD_io

      if (CS%TKE_diagnostics) then
        diag_TKE_MKE(i,j) = diag_TKE_MKE(i,j) + eCD%dTKE_MKE
        diag_TKE_conv(i,j) = diag_TKE_conv(i,j) + eCD%dTKE_conv
        diag_TKE_forcing(i,j) = diag_TKE_forcing(i,j) + eCD%dTKE_forcing
        diag_TKE_wind(i,j) = diag_TKE_wind(i,j) + eCD%dTKE_wind
        diag_TKE_mixing(i,j) = diag_TKE_mixing(i,j) + eCD%dTKE_mixing
        diag_TKE_mech_decay(i,j) = diag_TKE_mech_decay(i,j) + eCD%dTKE_mech_decay
        diag_TKE_conv_decay(i,j) = diag_TKE_conv_decay(i,j) + eCD%dTKE_conv_decay
       ! diag_TKE_unbalanced(i,j) = diag_TKE_unbalanced(i,j) + eCD%dTKE_unbalanced
      endif
      ! Write mixing length and velocity scale to 3-D arrays for diagnostic output
      if (CS%debug .or. (CS%id_Mixing_Length > 0)) then ; do K=1,nz+1
        diag_Mixing_Length(i,j,K) = mixlen(K)
      enddo ; endif
      if (CS%debug .or. (CS%id_Velocity_Scale > 0)) then ; do K=1,nz+1
        diag_Velocity_Scale(i,j,K) = mixvel(K)
      enddo ; endif
      if (BBL_mixing) then
        if (CS%debug .or. (CS%id_BBL_Mix_Length>0)) then ; do k=1,nz
          BBL_Mix_Length(i,j,k) = mixlen_BBL(k)
        enddo ; endif
        if (CS%debug .or. (CS%id_BBL_Vel_Scale>0)) then ; do k=1,nz
          BBL_Vel_Scale(i,j,k) = mixvel_BBL(k)
        enddo ; endif
        if (CS%id_TKE_BBL>0) &
          diag_TKE_BBL(i,j) = diag_TKE_BBL(i,j) + BBL_TKE
      endif
      if (CS%id_mstar_sfc > 0) diag_mstar_sfc(i,j) = eCD%mstar
      if (CS%id_mstar_bbl > 0) diag_mstar_BBL(i,j) = eCD%mstar_BBL
      if (CS%id_mstar_LT > 0) diag_mstar_lt(i,j) = eCD%mstar_LT
      if (CS%id_LA > 0) diag_LA(i,j) = eCD%LA
      if (CS%id_LA_mod > 0) diag_LA_mod(i,j) = eCD%LAmod
      if (report_avg_its) then
        CS%sum_its(1) = CS%sum_its(1) + real_to_EFP(real(eCD%OBL_its))
        CS%sum_its(2) = CS%sum_its(2) + real_to_EFP(1.0)
        if (BBL_mixing) then
          CS%sum_its_BBL(1) = CS%sum_its_BBL(1) + real_to_EFP(real(eCD%BBL_its))
          CS%sum_its_BBL(2) = CS%sum_its_BBL(2) + real_to_EFP(1.0)
        endif
      endif

      if (CS%options_diff > 0) then
        ! Call ePBL_column of ePBL_BBL_column with different parameter settings to diagnose sensitivities.
        ! These do not change the model state, and are only used for diagnostic purposes.
        if (CS%options_diff < 4) then
          BLD_1 = MLD_in ; BLD_2 = MLD_in
          do K=1,nz+1 ; SpV_dt_tmp(K) = SpV_scale1 * SpV_dt(K) ; enddo
          call ePBL_column(h, dz, u, v, T0, S0, dSV_dT_1d, dSV_dS_1d, SpV_dt_tmp, TKE_forcing, &
                           B_flux, absf, u_star, u_star_mean, mech_TKE, dt, BLD_1, Kd_1, &
                           mixvel, mixlen, GV, US, CS_tmp1, eCD_tmp, Waves, G, i, j)
          do K=1,nz+1 ; SpV_dt_tmp(K) = SpV_scale2 * SpV_dt(K) ; enddo
          call ePBL_column(h, dz, u, v, T0, S0, dSV_dT_1d, dSV_dS_1d, SpV_dt_tmp, TKE_forcing, &
                           B_flux, absf, u_star, u_star_mean, mech_TKE, dt, BLD_2, Kd_2, &
                           mixvel, mixlen, GV, US, CS_tmp2, eCD_tmp, Waves, G, i, j)
        else
          BLD_1 = BBLD_in ; BLD_2 = BBLD_in
          BBL_TKE = CS%ePBL_BBL_effic * GV%H_to_RZ * dt * visc%BBL_meanKE_loss(i,j)
          if ((CS%ePBL_tidal_effic > 0.0) .and. associated(fluxes%BBL_tidal_dis)) &
            BBL_TKE = BBL_TKE + CS%ePBL_tidal_effic * dt * fluxes%BBL_tidal_dis(i,j)
          u_star_BBL = max(visc%ustar_BBL(i,j), CS%ustar_min*GV%Z_to_H)
          u_star_BBL_z_t = u_star_bbl*GV%H_to_Z
          call ePBL_BBL_column(h, dz, u, v, T0, S0, dSV_dT_1d, dSV_dS_1d, SpV_dt, absf, dt, Kd, BBL_TKE, &
                               u_star_BBL, u_star_BBL_z_t, BBL_buoy_flux(i,j), Kd_1, BLD_1, mixvel_BBL, mixlen_BBL, &
                               GV, US, CS_tmp1, eCD_tmp)
          call ePBL_BBL_column(h, dz, u, v, T0, S0, dSV_dT_1d, dSV_dS_1d, SpV_dt, absf, dt, Kd, BBL_TKE, &
                               u_star_BBL, u_star_BBL_z_t, BBL_buoy_flux(i,j), Kd_2, BLD_2, mixvel_BBL, mixlen_BBL, &
                               GV, US, CS_tmp2, eCD_tmp)
        endif

        if (CS%id_opt_diff_Kd_ePBL > 0) then
          do K=1,nz+1 ; diff_Kd(i,j,K) = Kd_1(K) - Kd_2(K) ; enddo
        endif
        if (CS%id_opt_maxdiff_Kd_ePBL > 0) then
          max_abs_diff_Kd(i,j) = 0.0
          do K=1,nz+1 ; max_abs_diff_Kd(i,j) = max(max_abs_diff_Kd(i,j), abs(Kd_1(K) - Kd_2(K))) ; enddo
        endif
        if (CS%id_opt_diff_hML_depth > 0) diff_hML_depth(i,j) = BLD_1 - BLD_2
      endif

    else ! End of the ocean-point part of the i-loop
      ! For masked points, Kd_int must still be set (to 0) because it has intent out.
      do K=1,nz+1 ; Kd_2d(i,K) = 0. ; enddo
      CS%ML_depth(i,j) = 0.0
      CS%BBL_depth(i,j) = 0.0
    endif ; enddo ! Close of i-loop - Note the unusual loop order, with k-loops inside i-loops.

    do K=1,nz+1 ; do i=is,ie ; Kd_int(i,j,K) = Kd_2d(i,K) ; enddo ; enddo

  enddo ! j-loop

  if (CS%debug .and. BBL_mixing) then
    call hchksum(visc%BBL_meanKE_loss, "ePBL visc%BBL_meanKE_loss", G%HI, &
                 unscale=GV%H_to_MKS*US%L_T_to_m_s**2*US%s_to_T)
    call hchksum(visc%ustar_BBL, "ePBL visc%ustar_BBL", G%HI, unscale=GV%H_to_MKS*US%s_to_T)
    call hchksum(Kd_int, "End of ePBL Kd_int", G%HI, unscale=GV%H_to_MKS*US%Z_to_m*US%s_to_T)
    call hchksum(diag_Velocity_Scale, "ePBL Velocity_Scale", G%HI, unscale=US%Z_to_m*US%s_to_T)
    call hchksum(diag_Mixing_Length, "ePBL Mixing_Length", G%HI, unscale=US%Z_to_m)
    call hchksum(BBL_Vel_Scale, "ePBL BBL_Vel_Scale", G%HI, unscale=US%Z_to_m*US%s_to_T)
    call hchksum(BBL_Mix_Length, "ePBL BBL_Mix_Length", G%HI, unscale=US%Z_to_m)
  endif

  if (CS%id_ML_depth > 0) call post_data(CS%id_ML_depth, CS%ML_depth, CS%diag)
  if (CS%id_ustar_ePBL > 0) call post_data(CS%id_ustar_ePBL, diag_ustar, CS%diag)
  if (CS%id_bflx_ePBL > 0) call post_data(CS%id_bflx_ePBL, buoy_flux, CS%diag)
  if (CS%id_hML_depth > 0) call post_data(CS%id_hML_depth, CS%ML_depth, CS%diag)
  if (CS%id_TKE_wind > 0) call post_data(CS%id_TKE_wind, diag_TKE_wind, CS%diag)
  if (CS%id_TKE_MKE > 0)  call post_data(CS%id_TKE_MKE, diag_TKE_MKE, CS%diag)
  if (CS%id_TKE_conv > 0) call post_data(CS%id_TKE_conv, diag_TKE_conv, CS%diag)
  if (CS%id_TKE_forcing > 0) call post_data(CS%id_TKE_forcing, diag_TKE_forcing, CS%diag)
  if (CS%id_TKE_mixing > 0) call post_data(CS%id_TKE_mixing, diag_TKE_mixing, CS%diag)
  if (CS%id_TKE_mech_decay > 0) &
    call post_data(CS%id_TKE_mech_decay, diag_TKE_mech_decay, CS%diag)
  if (CS%id_TKE_conv_decay > 0) &
    call post_data(CS%id_TKE_conv_decay, diag_TKE_conv_decay, CS%diag)
  if (CS%id_Mixing_Length > 0) call post_data(CS%id_Mixing_Length, diag_Mixing_Length, CS%diag)
  if (CS%id_Velocity_Scale >0) call post_data(CS%id_Velocity_Scale, diag_Velocity_Scale, CS%diag)
  if (CS%id_mstar_sfc > 0)     call post_data(CS%id_mstar_sfc, diag_mstar_sfc, CS%diag)
  if (BBL_mixing) then
    if (CS%id_Kd_BBL > 0) call post_data(CS%id_Kd_BBL, Kd_BBL_3d, CS%diag)
    if (CS%id_BBL_Mix_Length > 0) call post_data(CS%id_BBL_Mix_Length, BBL_Mix_Length, CS%diag)
    if (CS%id_BBL_Vel_Scale > 0) call post_data(CS%id_BBL_Vel_Scale, BBL_Vel_Scale, CS%diag)
    if (CS%id_ustar_BBL > 0) call post_data(CS%id_ustar_BBL, diag_ustar_BBL, CS%diag)
    if (CS%id_BBL_decay_scale > 0) call post_data(CS%id_BBL_decay_scale, diag_BBL_decay_scale, CS%diag)
    if (CS%id_TKE_BBL > 0) call post_data(CS%id_TKE_BBL, diag_TKE_BBL, CS%diag)
    if (CS%id_TKE_BBL_mixing > 0) call post_data(CS%id_TKE_BBL_mixing, diag_TKE_BBL_mixing, CS%diag)
    if (CS%id_TKE_BBL_decay > 0) call post_data(CS%id_TKE_BBL_decay, diag_TKE_BBL_decay, CS%diag)
    if (CS%id_BBL_depth > 0) call post_data(CS%id_BBL_depth, CS%BBL_depth, CS%diag)
    if (CS%id_mstar_BBL > 0)     call post_data(CS%id_mstar_BBL, diag_mstar_BBL, CS%diag)
  endif
  if (CS%id_LA > 0)       call post_data(CS%id_LA, diag_LA, CS%diag)
  if (CS%id_LA_mod > 0)   call post_data(CS%id_LA_mod, diag_LA_mod, CS%diag)
  if (CS%id_mstar_LT > 0) call post_data(CS%id_mstar_LT, diag_mstar_LT, CS%diag)
  if (stoch_CS%pert_epbl) then
    if (stoch_CS%id_epbl1_wts > 0) call post_data(stoch_CS%id_epbl1_wts, stoch_CS%epbl1_wts, CS%diag)
    if (stoch_CS%id_epbl2_wts > 0) call post_data(stoch_CS%id_epbl2_wts, stoch_CS%epbl2_wts, CS%diag)
  endif

  if (CS%options_diff > 0) then
    ! These diagnostics are only for determining sensitivities to different ePBL settings.
    if (CS%id_opt_diff_Kd_ePBL > 0)    call post_data(CS%id_opt_diff_Kd_ePBL, diff_Kd, CS%diag)
    if (CS%id_opt_maxdiff_Kd_ePBL > 0) call post_data(CS%id_opt_maxdiff_Kd_ePBL, max_abs_diff_Kd, CS%diag)
    if (CS%id_opt_diff_hML_depth > 0)  call post_data(CS%id_opt_diff_hML_depth, diff_hML_depth, CS%diag)
  endif

end subroutine energetic_PBL



!> This subroutine determines the diffusivities from the integrated energetics
!!  mixed layer model for a single column of water.
subroutine ePBL_column(h, dz, u, v, T0, S0, dSV_dT, dSV_dS, SpV_dt, TKE_forcing, B_flux, absf, &
                       u_star, u_star_mean, mech_TKE_in, dt, MLD_io, Kd, mixvel, mixlen, GV, US, CS, eCD, &
                       Waves, G, i, j, TKE_gen_stoch, TKE_diss_stoch)
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZK_(GV)), intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZK_(GV)), intent(in)  :: dz     !< The vertical distance across layers [Z ~> m].
  real, dimension(SZK_(GV)), intent(in)  :: u      !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZK_(GV)), intent(in)  :: v      !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZK_(GV)), intent(in)  :: T0     !< The initial layer temperatures [C ~> degC].
  real, dimension(SZK_(GV)), intent(in)  :: S0     !< The initial layer salinities [S ~> ppt].

  real, dimension(SZK_(GV)), intent(in)  :: dSV_dT !< The partial derivative of in-situ specific
                                                   !! volume with potential temperature
                                                   !! [R-1 C-1 ~> m3 kg-1 degC-1].
  real, dimension(SZK_(GV)), intent(in)  :: dSV_dS !< The partial derivative of in-situ specific
                                                   !! volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZK_(GV)+1), intent(in) :: SpV_dt !< Specific volume interpolated to interfaces
                                                   !! divided by dt or 1.0 / (dt * Rho0), times conversion
                                                   !! factors for answer dates before 20240101 in
                                                   !! [m3 Z-3 R-1 T2 s-3 ~> m3 kg-1 s-1] or without
                                                   !! the conversion factors for answer dates of
                                                   !! 20240101 and later in [R-1 T-1 ~> m3 kg-1 s-1],
                                                   !! used to convert local TKE into a turbulence
                                                   !! velocity cubed.
  real, dimension(SZK_(GV)), intent(in)  :: TKE_forcing !< The forcing requirements to homogenize the
                                                   !! forcing that has been applied to each layer
                                                   !! [R Z3 T-2 ~> J m-2].
  real,                    intent(in)    :: B_flux !< The surface buoyancy flux [Z2 T-3 ~> m2 s-3]
  real,                    intent(in)    :: absf   !< The absolute value of the Coriolis parameter [T-1 ~> s-1].
  real,                    intent(in)    :: u_star !< The surface friction velocity [Z T-1 ~> m s-1].
  real,                    intent(in)    :: u_star_mean !< The surface friction velocity without any
                                                   !! contribution from unresolved gustiness  [Z T-1 ~> m s-1].
  real,                    intent(in)    :: mech_TKE_in !< The mechanically generated turbulent
                                                   !! kinetic energy available for mixing over a time
                                                   !! step before the application of the efficiency
                                                   !! in mstar. [R Z3 T-2 ~> J m-2].
  real,                    intent(inout) :: MLD_io !< A first guess at the mixed layer depth on input, and
                                                   !! the calculated mixed layer depth on output [Z ~> m]
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: Kd     !< The diagnosed diffusivities at interfaces
                                                   !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: mixvel !< The mixing velocity scale used in Kd
                                                   !! [Z T-1 ~> m s-1].
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: mixlen !< The mixing length scale used in Kd [Z ~> m].
  type(energetic_PBL_CS),  intent(in)    :: CS     !< Energetic PBL control structure
  type(ePBL_column_diags), intent(inout) :: eCD    !< A container for passing around diagnostics.
  type(wave_parameters_CS), pointer      :: Waves  !< Waves control structure for Langmuir turbulence
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  integer,                 intent(in)    :: i      !< The i-index to work on (used for Waves)
  integer,                 intent(in)    :: j      !< The j-index to work on (used for Waves)
  real,          optional, intent(in)    :: TKE_gen_stoch  !< random factor used to perturb TKE generation [nondim]
  real,          optional, intent(in)    :: TKE_diss_stoch !< random factor used to perturb TKE dissipation [nondim]

!    This subroutine determines the diffusivities in a single column from the integrated energetics
!  planetary boundary layer (ePBL) model.  It assumes that heating, cooling and freshwater fluxes
!  have already been applied.  All calculations are done implicitly, and there
!  is no stability limit on the time step.
!
!    For each interior interface, first discard the TKE to account for mixing
! of shortwave radiation through the next denser cell.  Next drive mixing based
! on the local? values of ustar + wstar, subject to available energy.  This
! step sets the value of Kd(K).  Any remaining energy is then subject to decay
! before being handed off to the next interface.  mech_TKE and conv_PErel are treated
! separately for the purposes of decay, but are used proportionately to drive
! mixing.

  ! Local variables
  real, dimension(SZK_(GV)+1) :: &
    pres_Z, &       ! Interface pressures with a rescaling factor to convert interface height
                    ! movements into changes in column potential energy [R Z2 T-2 ~> kg m-1 s-2].
    hb_hs           ! The distance from the bottom over the thickness of the
                    ! water column [nondim].
  real :: mech_TKE  !   The mechanically generated turbulent kinetic energy
                    ! available for mixing over a time step [R Z3 T-2 ~> J m-2].
  real :: conv_PErel ! The potential energy that has been convectively released
                    ! during this timestep [R Z3 T-2 ~> J m-2]. A portion nstar_FC
                    ! of conv_PErel is available to drive mixing.
  real :: htot      ! The total thickness of the layers above an interface [H ~> m or kg m-2].
  real :: dztot     ! The total depth of the layers above an interface [Z ~> m].
  real :: uhtot     ! The depth integrated zonal velocities in the layers above [H L T-1 ~> m2 s-1 or kg m-1 s-1]
  real :: vhtot     ! The depth integrated meridional velocities in the layers above [H L T-1 ~> m2 s-1 or kg m-1 s-1]
  real :: Idecay_len_TKE  ! The inverse of a turbulence decay length scale [H-1 ~> m-1 or m2 kg-1].
  real :: dz_sum    ! The total thickness of the water column [Z ~> m].

  real, dimension(SZK_(GV)) :: &
    dT_to_dColHt, & ! Partial derivative of the total column height with the temperature changes
                    ! within a layer [Z C-1 ~> m degC-1].
    dS_to_dColHt, & ! Partial derivative of the total column height with the salinity changes
                    ! within a layer  [Z S-1 ~> m ppt-1].
    dT_to_dPE, &    ! Partial derivatives of column potential energy with the temperature
                    ! changes within a layer, in [R Z3 T-2 C-1 ~> J m-2 degC-1].
    dS_to_dPE, &    ! Partial derivatives of column potential energy with the salinity changes
                    ! within a layer, in [R Z3 T-2 S-1 ~> J m-2 ppt-1].
    dT_to_dColHt_a, & ! Partial derivative of the total column height with the temperature changes
                    ! within a layer, including the implicit effects  of mixing with layers higher
                    ! in the water column [Z C-1 ~> m degC-1].
    dS_to_dColHt_a, & ! Partial derivative of the total column height with the salinity changes
                    ! within a layer, including the implicit effects  of mixing with layers higher
                    ! in the water column [Z S-1 ~> m ppt-1].
    dT_to_dPE_a, &  ! Partial derivatives of column potential energy with the temperature changes
                    ! within a layer, including the implicit effects of mixing with layers higher
                    ! in the water column [R Z3 T-2 C-1 ~> J m-2 degC-1].
    dS_to_dPE_a, &  ! Partial derivative of column potential energy with the salinity changes
                    ! within a layer, including the implicit effects of mixing with layers higher
                    ! in the water column [R Z3 T-2 S-1 ~> J m-2 ppt-1].
    c1, &           ! c1 is used by the tridiagonal solver [nondim].
    Te, &           ! Estimated final values of T in the column [C ~> degC].
    Se, &           ! Estimated final values of S in the column [S ~> ppt].
    dTe, &          ! Running (1-way) estimates of temperature change [C ~> degC].
    dSe, &          ! Running (1-way) estimates of salinity change [S ~> ppt].
    hp_a, &         ! An effective pivot thickness of the layer including the effects
                    ! of coupling with layers above [H ~> m or kg m-2].  This is the first term
                    ! in the denominator of b1 in a downward-oriented tridiagonal solver.
    Th_a, &         ! An effective temperature times a thickness in the layer above, including implicit
                    ! mixing effects with other yet higher layers [C H ~> degC m or degC kg m-2].
    Sh_a, &         ! An effective salinity times a thickness in the layer above, including implicit
                    ! mixing effects with other yet higher layers [S H ~> ppt m or ppt kg m-2].
    Th_b, &         ! An effective temperature times a thickness in the layer below, including implicit
                    ! mixing effects with other yet lower layers [C H ~> degC m or degC kg m-2].
    Sh_b            ! An effective salinity times a thickness in the layer below, including implicit
                    ! mixing effects with other yet lower layers [S H ~> ppt m or ppt kg m-2].
  real, dimension(SZK_(GV)+1) :: &
    MixLen_shape, & ! A nondimensional shape factor for the mixing length that
                    ! gives it an appropriate asymptotic value at the bottom of
                    ! the boundary layer [nondim].
    h_dz_int, &     ! The ratio of the layer thicknesses over the vertical distances
                    ! across the layers surrounding an interface [H Z-1 ~> nondim or kg m-3]
    Kddt_h          ! The total diapycnal diffusivity at an interface times a timestep divided by the
                    ! average thicknesses around an interface [H ~> m or kg m-2].
  real :: b1        ! b1 is inverse of the pivot used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: dz_neglect ! A vertical distance that is so small it is usually lost
                    ! in roundoff and can be neglected [Z ~> m].
  real :: dMass     ! The mass per unit area within a layer [Z R ~> kg m-2].
  real :: dPres     ! The hydrostatic pressure change across a layer [R Z2 T-2 ~> Pa = J m-3].
  real :: dMKE_max  ! The maximum amount of mean kinetic energy that could be
                    ! converted to turbulent kinetic energy if the velocity in
                    ! the layer below an interface were homogenized with all of
                    ! the water above the interface [R Z3 T-2 ~> J m-2].
  real :: MKE2_Hharm ! Twice the inverse of the harmonic mean of the thickness
                    ! of a layer and the thickness of the water above, used in
                    ! the MKE conversion equation [H-1 ~> m-1 or m2 kg-1].

  real :: dt_h      ! The timestep divided by the averages of the vertical distances around
                    ! a layer [T Z-1 ~> s m-1].
  real :: dz_bot    ! The distance from the bottom [Z ~> m].
  real :: dz_rsum   ! The running sum of dz from the top [Z ~> m].
  real :: I_dzsum   ! The inverse of dz_sum [Z-1 ~> m-1].
  real :: I_MLD     ! The inverse of the current value of MLD [Z-1 ~> m-1].
  real :: dz_tt     ! The distance from the surface or up to the next interface
                    ! that did not exhibit turbulent mixing from this scheme plus
                    ! a surface mixing roughness length given by dz_tt_min [Z ~> m].
  real :: dz_tt_min  ! A surface roughness length [Z ~> m].

  real :: C1_3      ! = 1/3  [nondim]
  real :: vstar     ! An in-situ turbulent velocity [Z T-1 ~> m s-1].
  real :: mstar_total ! The value of mstar used in ePBL [nondim]
  real :: mstar_LT  ! An addition to mstar due to Langmuir turbulence [nondim] (output for diagnostic)
  real :: MLD_output ! The mixed layer depth output from this routine [Z ~> m]
  real :: LA        ! The value of the Langmuir number [nondim]
  real :: LAmod     ! The modified Langmuir number by convection [nondim]
  real :: hbs_here  ! The local minimum of hb_hs and MixLen_shape [nondim]
  real :: nstar_FC  ! The fraction of conv_PErel that can be converted to mixing [nondim].
  real :: TKE_reduc ! The fraction by which TKE and other energy fields are
                    ! reduced to support mixing [nondim]. between 0 and 1.
  real :: tot_TKE   ! The total TKE available to support mixing at interface K [R Z3 T-2 ~> J m-2].
  real :: TKE_here  ! The total TKE at this point in the algorithm [R Z3 T-2 ~> J m-2].
  real :: dT_km1_t2 ! A diffusivity-independent term related to the temperature
                    ! change in the layer above the interface [C ~> degC].
  real :: dS_km1_t2 ! A diffusivity-independent term related to the salinity
                    ! change in the layer above the interface [S ~> ppt].
  real :: dTe_term  ! A diffusivity-independent term related to the temperature
                    ! change in the layer below the interface [C H ~> degC m or degC kg m-2].
  real :: dSe_term  ! A diffusivity-independent term related to the salinity
                    ! change in the layer above the interface [S H ~> ppt m or ppt kg m-2].
  real :: dTe_t2    ! A part of dTe_term [C H ~> degC m or degC kg m-2].
  real :: dSe_t2    ! A part of dSe_term [S H ~> ppt m or ppt kg m-2].
  real :: dPE_conv  ! The convective change in column potential energy [R Z3 T-2 ~> J m-2].
  real :: MKE_src   ! The mean kinetic energy source of TKE due to Kddt_h(K) [R Z3 T-2 ~> J m-2].
  real :: dMKE_src_dK  ! The partial derivative of MKE_src with Kddt_h(K) [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real :: Kd_guess0    ! A first guess of the diapycnal diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real :: PE_chg_g0    ! The potential energy change when Kd is Kd_guess0 [R Z3 T-2 ~> J m-2]
  real :: Kddt_h_g0    ! The first guess of the change in diapycnal diffusivity times a timestep
                       ! divided by the average thicknesses around an interface [H ~> m or kg m-2].
  real :: PE_chg_max   ! The maximum PE change for very large values of Kddt_h(K) [R Z3 T-2 ~> J m-2].
  real :: dPEc_dKd_Kd0 ! The partial derivative of PE change with Kddt_h(K)
                       ! for very small values of Kddt_h(K) [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real :: PE_chg    ! The change in potential energy due to mixing at an
                    ! interface [R Z3 T-2 ~> J m-2], positive for the column increasing
                    ! in potential energy (i.e., consuming TKE).
  real :: TKE_left  ! The amount of turbulent kinetic energy left for the most
                    ! recent guess at Kddt_h(K) [R Z3 T-2 ~> J m-2].
  real :: dPEc_dKd  ! The partial derivative of PE_chg with Kddt_h(K) [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real :: TKE_left_min, TKE_left_max ! Maximum and minimum values of TKE_left [R Z3 T-2 ~> J m-2].
  real :: Kddt_h_max, Kddt_h_min ! Maximum and minimum values of Kddt_h(K) [H ~> m or kg m-2].
  real :: Kddt_h_guess ! A guess at the value of Kddt_h(K) [H ~> m or kg m-2].
  real :: Kddt_h_next  ! The next guess at the value of Kddt_h(K) [H ~> m or kg m-2].
  real :: dKddt_h      ! The change between guesses at Kddt_h(K) [H ~> m or kg m-2].
  real :: dKddt_h_Newt ! The change between guesses at Kddt_h(K) with Newton's method [H ~> m or kg m-2].
  real :: Kddt_h_newt  ! The Newton's method next guess for Kddt_h(K) [H ~> m or kg m-2].
  real :: exp_kh    ! The nondimensional decay of TKE across a layer [nondim].
  real :: vstar_unit_scale ! A unit conversion factor for turbulent velocities [Z T-1 s m-1 ~> 1]
  logical :: use_Newt  ! Use Newton's method for the next guess at Kddt_h(K).
  logical :: convectively_stable ! If true the water column is convectively stable at this interface.
  logical :: sfc_connected   ! If true the ocean is actively turbulent from the present
                    ! interface all the way up to the surface.
  logical :: sfc_disconnect ! If true, any turbulence has become disconnected
                    ! from the surface.

  ! The following is only used for diagnostics.
  real :: I_dtdiag  !  = 1.0 / dt [T-1 ~> s-1].

  !----------------------------------------------------------------------
  !/BGR added Aug24,2016 for adding iteration to get boundary layer depth
  !    - needed to compute new mixing length.
  real :: MLD_guess, MLD_found ! Mixing Layer depth guessed/found for iteration [Z ~> m]
  real :: min_MLD, max_MLD ! Iteration bounds on MLD [Z ~> m], which are adjusted at each step
                    !  - These are initialized based on surface/bottom
                    !  1. The iteration guesses a value (possibly from prev step or neighbor).
                    !  2. The iteration checks if value is converged, too shallow, or too deep.
                    !  3. Based on result adjusts the Max/Min and searches through the water column.
                    !  - If using an accurate guess the iteration is very quick (e.g. if MLD doesn't
                    !    change over timestep).  Otherwise it takes 5-10 passes, but has a high
                    !    convergence rate.  Other iteration may be tried, but this method seems to
                    !    fail very rarely and the added cost is likely not significant.
                    !    Additionally, when it fails to converge it does so in a reasonable
                    !    manner giving a usable guess. When it does fail, it is due to convection
                    !    within the boundary layer.  Likely, a new method e.g. surface_disconnect,
                    !    can improve this.
  real :: dMLD_min  ! The change in diagnosed mixed layer depth when the guess is min_MLD [Z ~> m]
  real :: dMLD_max  ! The change in diagnosed mixed layer depth when the guess is max_MLD [Z ~> m]
  integer :: OBL_it        ! Iteration counter

  real :: TKE_used  ! The TKE used to support mixing at an interface [R Z3 T-2 ~> J m-2].
  ! real :: Kd_add    ! The additional diffusivity at an interface [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real :: frac_in_BL ! The fraction of the energy required to support dKd_max that is suppiled by
                     ! max_PE_chg, used here to determine a fractional layer contribution to the
                     ! boundary layer thickness [nondim]
  real :: Surface_Scale ! Surface decay scale for vstar [nondim]
  logical :: calc_Te    ! If true calculate the expected final temperature and salinity values.
  logical :: no_MKE_conversion ! If true, there is no conversion from MKE to TKE, so a simpler solver can be used.
  logical :: debug      ! This is used as a hard-coded value for debugging.
  logical :: convectively_unstable  ! If true, there is convective instability at an interface.

  !  The following arrays are used only for debugging purposes.
  real :: dPE_debug     ! An estimate of the potential energy change [R Z3 T-2 ~> J m-2]
  real :: mixing_debug  ! An estimate of the rate of change of potential energy due to mixing [R Z3 T-3 ~> W m-2]
  real, dimension(20) :: TKE_left_itt   ! The value of TKE_left after each iteration [R Z3 T-2 ~> J m-2]
  real, dimension(20) :: PE_chg_itt     ! The value of PE_chg after each iteration [R Z3 T-2 ~> J m-2]
  real, dimension(20) :: Kddt_h_itt     ! The value of Kddt_h_guess after each iteration [H ~> m or kg m-2]
  real, dimension(20) :: dPEa_dKd_itt   ! The value of dPEc_dKd after each iteration [R Z3 T-2 H-1 ~> J m-3 or J kg-1]
  real, dimension(20) :: MKE_src_itt    ! The value of MKE_src after each iteration [R Z3 T-2 ~> J m-2]
  real, dimension(SZK_(GV)) :: mech_TKE_k  ! The mechanically generated turbulent kinetic energy
                    ! available for mixing over a time step for each layer [R Z3 T-2 ~> J m-2].
  real, dimension(SZK_(GV)) :: conv_PErel_k ! The potential energy that has been convectively released
                    ! during this timestep for each layer [R Z3 T-2 ~> J m-2].
  real, dimension(SZK_(GV)) :: nstar_k   ! The fraction of conv_PErel that can be converted to mixing
                    ! for each layer [nondim].
  real, dimension(SZK_(GV)) :: dT_expect ! Expected temperature changes [C ~> degC]
  real, dimension(SZK_(GV)) :: dS_expect ! Expected salinity changes [S ~> ppt]
  integer, dimension(SZK_(GV)) :: num_itts

  integer :: k, nz, itt, max_itt

  ! variables for ML based diffusivity
  real :: v0_ML_turb_vel_scale ! turbulence vel scale from ML in diffusivity [Z T-1 ~> m s-1]

  nz = GV%ke

  debug = .false.  ! Change this hard-coded value for debugging.
  calc_Te = (debug .or. (.not.CS%orig_PE_calc))
  no_MKE_conversion = ((CS%direct_calc) .and. (CS%MKE_to_TKE_effic == 0.0))

  h_neglect = GV%H_subroundoff
  dz_neglect = GV%dZ_subroundoff

  C1_3 = 1.0 / 3.0
  I_dtdiag = 1.0 / dt
  max_itt = 20

  dz_tt_min = 0.0
  if (CS%answer_date < 20240101) vstar_unit_scale = US%m_to_Z * US%T_to_s

  MLD_guess = MLD_io

!   Determine the initial mech_TKE and conv_PErel, including the energy required
! to mix surface heating through the topmost cell, the energy released by mixing
! surface cooling & brine rejection down through the topmost cell, and
! homogenizing the shortwave heating within that cell.  This sets the energy
! and ustar and wstar available to drive mixing at the first interior
! interface.

  do K=1,nz+1 ; Kd(K) = 0.0 ; enddo

  pres_Z(1) = 0.0
  do k=1,nz
    dMass = GV%H_to_RZ * h(k)
    dPres = GV%g_Earth_Z_T2 * dMass
    dT_to_dPE(k) = (dMass * (pres_Z(K) + 0.5*dPres)) * dSV_dT(k)
    dS_to_dPE(k) = (dMass * (pres_Z(K) + 0.5*dPres)) * dSV_dS(k)
    dT_to_dColHt(k) = dMass * dSV_dT(k)
    dS_to_dColHt(k) = dMass * dSV_dS(k)

    pres_Z(K+1) = pres_Z(K) + dPres
  enddo

  ! Determine the total thickness (dz_sum) and the fractional distance from the bottom (hb_hs).
  dz_sum = dz_neglect ; do k=1,nz ; dz_sum = dz_sum + dz(k) ; enddo
  I_dzsum = 0.0 ; if (dz_sum > 0.0) I_dzsum = 1.0 / dz_sum
  dz_bot = 0.0
  hb_hs(nz+1) = 0.0
  do k=nz,1,-1
    dz_bot = dz_bot + dz(k)
    hb_hs(K) = dz_bot * I_dzsum
  enddo

  MLD_output = dz(1)

  !/The following lines are for the iteration over MLD
  ! max_MLD will initialized as ocean bottom depth
  max_MLD = 0.0 ; do k=1,nz ; max_MLD = max_MLD + dz(k) ; enddo
  ! min_MLD will be initialized to 0.
  min_MLD = 0.0
  ! Set values of the wrong signs to indicate that these changes are not based on valid estimates
  dMLD_min = -1.0*US%m_to_Z ; dMLD_max = 1.0*US%m_to_Z

  ! If no first guess is provided for MLD, try the middle of the water column
  if (MLD_guess <= min_MLD) MLD_guess = 0.5 * (min_MLD + max_MLD)

  if (GV%Boussinesq) then
    do K=1,nz+1 ; h_dz_int(K) = GV%Z_to_H ; enddo
  else
    h_dz_int(1) = (h(1) + h_neglect) / (dz(1) + dz_neglect)
    do K=2,nz
      h_dz_int(K) = (h(k-1) + h(k) + h_neglect) / (dz(k-1) + dz(k) + dz_neglect)
    enddo
    h_dz_int(nz+1) = (h(nz) + h_neglect) / (dz(nz) + dz_neglect)
  endif

  ! Iterate to determine a converged EPBL depth.
  do OBL_it=1,CS%Max_MLD_Its

    if (debug) then ; mech_TKE_k(:) = 0.0 ; conv_PErel_k(:) = 0.0 ; endif

    ! Reset ML_depth
    MLD_output = dz(1)
    sfc_connected = .true.

    !/ Here we get mstar, which is the ratio of convective TKE driven mixing to UStar**3
    if (CS%Use_LT) then
      call get_Langmuir_Number(LA, G, GV, US, abs(MLD_guess), u_star_mean, i, j, dz, Waves, &
                               U_H=u, V_H=v)
      call find_mstar(CS, US, B_flux, u_star, MLD_guess, absf, .false., &
                      mstar_total, Langmuir_Number=La, Convect_Langmuir_Number=LAmod,&
                      mstar_LT=mstar_LT)
    else
      call find_mstar(CS, US, B_flux, u_star, MLD_guess, absf, .false., mstar_total)
    endif

    !/ Apply mstar to get mech_TKE
    if ((CS%answer_date < 20190101) .and. (CS%mstar_scheme==Use_Fixed_mstar)) then
      mech_TKE = (dt*mstar_total*GV%Rho0) * u_star**3
    else
      mech_TKE = mstar_total * mech_TKE_in
      ! mech_TKE = mstar_total * (dt*GV%Rho0* u_star**3)
    endif
    ! stochastically perturb mech_TKE in the UFS
    if (present(TKE_gen_stoch)) mech_TKE = mech_TKE*TKE_gen_stoch

    if (CS%TKE_diagnostics) then
      eCD%dTKE_conv = 0.0 ; eCD%dTKE_mixing = 0.0
      eCD%dTKE_MKE = 0.0 ; eCD%dTKE_mech_decay = 0.0 ; eCD%dTKE_conv_decay = 0.0

      eCD%dTKE_wind = mech_TKE * I_dtdiag
      if (TKE_forcing(1) <= 0.0) then
        eCD%dTKE_forcing = max(-mech_TKE, TKE_forcing(1)) * I_dtdiag
        ! eCD%dTKE_unbalanced = min(0.0, TKE_forcing(1) + mech_TKE) * I_dtdiag
      else
        eCD%dTKE_forcing = CS%nstar*TKE_forcing(1) * I_dtdiag
        ! eCD%dTKE_unbalanced = 0.0
      endif
    endif

    if (TKE_forcing(1) <= 0.0) then
      mech_TKE = mech_TKE + TKE_forcing(1)
      if (mech_TKE < 0.0) mech_TKE = 0.0
      conv_PErel = 0.0
    else
      conv_PErel = TKE_forcing(1)
    endif


    ! Store in 1D arrays for output.
    do K=1,nz+1 ; mixvel(K) = 0.0 ; mixlen(K) = 0.0 ; enddo

    ! Determine the mixing shape function MixLen_shape.
    if ((.not.CS%Use_MLD_iteration) .or. &
        (CS%transLay_scale >= 1.0) .or. (CS%transLay_scale < 0.0) ) then
      do K=1,nz+1
        MixLen_shape(K) = 1.0
      enddo
    elseif (MLD_guess <= 0.0) then
      if (CS%transLay_scale > 0.0) then ; do K=1,nz+1
        MixLen_shape(K) = CS%transLay_scale
      enddo ; else ; do K=1,nz+1
        MixLen_shape(K) = 1.0
      enddo ; endif
    else
      ! Reduce the mixing length based on MLD, with a quadratic
      ! expression that follows KPP.
      I_MLD = 1.0 / MLD_guess
      dz_rsum = 0.0
      MixLen_shape(1) = 1.0
      if (CS%eqdisc) then ! update Kd as per Machine Learning equation discovery
        call kappa_eqdisc(MixLen_shape, CS, GV, h, absf, B_flux, u_star, MLD_guess)
      else
        do K=2,nz+1
          dz_rsum = dz_rsum + dz(k-1)
          if (CS%MixLenExponent==2.0) then
            MixLen_shape(K) = CS%transLay_scale + (1.0 - CS%transLay_scale) * &
               (max(0.0, (MLD_guess - dz_rsum)*I_MLD) )**2 ! CS%MixLenExponent
          else
            MixLen_shape(K) = CS%transLay_scale + (1.0 - CS%transLay_scale) * &
               (max(0.0, (MLD_guess - dz_rsum)*I_MLD) )**CS%MixLenExponent
          endif
        enddo
      endif
    endif

    v0_ML_turb_vel_scale = 0.0 ! a variable that gets passed on to get_eqdisc_v0 & get_eqdisc_v0h
    if (CS%eqdisc_v0) then
      call get_eqdisc_v0(CS,absf,B_flux,u_star,v0_ML_turb_vel_scale)
    elseif (CS%eqdisc_v0h) then
      call get_eqdisc_v0h(CS,B_flux,u_star,MLD_guess,v0_ML_turb_vel_scale)
    endif

    Kd(1) = 0.0 ; Kddt_h(1) = 0.0
    hp_a(1) = h(1)
    dT_to_dPE_a(1) = dT_to_dPE(1) ; dT_to_dColHt_a(1) = dT_to_dColHt(1)
    dS_to_dPE_a(1) = dS_to_dPE(1) ; dS_to_dColHt_a(1) = dS_to_dColHt(1)

    htot = h(1) ; dztot = dz(1) ; uhtot = u(1)*h(1) ; vhtot = v(1)*h(1)

    if (debug) then
      mech_TKE_k(1) = mech_TKE ; conv_PErel_k(1) = conv_PErel
      nstar_k(:) = 0.0 ; nstar_k(1) = CS%nstar ; num_itts(:) = -1
    endif

    do K=2,nz
      ! Apply dissipation to the TKE, here applied as an exponential decay
      ! due to 3-d turbulent energy being lost to inefficient rotational modes.

      !   There should be several different "flavors" of TKE that decay at
      ! different rates.  The following form is often used for mechanical
      ! stirring from the surface, perhaps due to breaking surface gravity
      ! waves and wind-driven turbulence.
      if (GV%Boussinesq) then
        Idecay_len_TKE = (CS%TKE_decay * absf / u_star) * GV%H_to_Z
      else
        Idecay_len_TKE = (CS%TKE_decay * absf) / (h_dz_int(K) * u_star)
      endif
      exp_kh = 1.0
      if (Idecay_len_TKE > 0.0) exp_kh = exp(-h(k-1)*Idecay_len_TKE)
      if (CS%TKE_diagnostics) &
        eCD%dTKE_mech_decay = eCD%dTKE_mech_decay + (exp_kh-1.0) * mech_TKE * I_dtdiag
      if (present(TKE_diss_stoch)) then ! perturb the TKE destruction
        mech_TKE = mech_TKE * (1.0 + (exp_kh-1.0) * TKE_diss_stoch)
      else
        mech_TKE = mech_TKE * exp_kh
      endif

      !   Accumulate any convectively released potential energy to contribute
      ! to wstar and to drive penetrating convection.
      if (TKE_forcing(k) > 0.0) then
        conv_PErel = conv_PErel + TKE_forcing(k)
        if (CS%TKE_diagnostics) &
          eCD%dTKE_forcing = eCD%dTKE_forcing + CS%nstar*TKE_forcing(k) * I_dtdiag
      endif

      if (debug) then
        mech_TKE_k(K) = mech_TKE ; conv_PErel_k(K) = conv_PErel
      endif

      !  Determine the total energy
      nstar_FC = CS%nstar
      if (CS%nstar * conv_PErel > 0.0) then
        ! Here nstar is a function of the natural Rossby number 0.2/(1+0.2/Ro), based
        ! on a curve fit from the data of Wang (GRL, 2003).
        ! Note:         Ro = 1.0 / sqrt(0.5 * dt * Rho0 * (absf*dztot)**3 / conv_PErel)
        if (GV%Boussinesq) then
          nstar_FC = CS%nstar * conv_PErel / (conv_PErel + 0.2 * &
                     sqrt(0.5 * dt * GV%Rho0 * (absf*dztot)**3 * conv_PErel))
        else
          nstar_FC = CS%nstar * conv_PErel / (conv_PErel + 0.2 * &
                     sqrt(0.5 * dt * GV%H_to_RZ * (absf**3 * (dztot**2 * htot)) * conv_PErel))
        endif
      endif

      if (debug) nstar_k(K) = nstar_FC

      tot_TKE = mech_TKE + nstar_FC * conv_PErel

      !   For each interior interface, first discard the TKE to account for
      ! mixing of shortwave radiation through the next denser cell.
      if (TKE_forcing(k) < 0.0) then
        if (TKE_forcing(k) + tot_TKE < 0.0) then
          ! The shortwave requirements deplete all the energy in this layer.
          if (CS%TKE_diagnostics) then
            eCD%dTKE_mixing = eCD%dTKE_mixing + tot_TKE * I_dtdiag
            eCD%dTKE_forcing = eCD%dTKE_forcing - tot_TKE * I_dtdiag
            ! eCD%dTKE_unbalanced = eCD%dTKE_unbalanced + (TKE_forcing(k) + tot_TKE) * I_dtdiag
            eCD%dTKE_conv_decay = eCD%dTKE_conv_decay + (CS%nstar-nstar_FC) * conv_PErel * I_dtdiag
          endif
          tot_TKE = 0.0 ; mech_TKE = 0.0 ; conv_PErel = 0.0
        else
          ! Reduce the mechanical and convective TKE proportionately.
          TKE_reduc = (tot_TKE + TKE_forcing(k)) / tot_TKE
          if (CS%TKE_diagnostics) then
            eCD%dTKE_mixing = eCD%dTKE_mixing - TKE_forcing(k) * I_dtdiag
            eCD%dTKE_forcing = eCD%dTKE_forcing + TKE_forcing(k) * I_dtdiag
            eCD%dTKE_conv_decay = eCD%dTKE_conv_decay + &
                (1.0-TKE_reduc)*(CS%nstar-nstar_FC) * conv_PErel * I_dtdiag
          endif
          tot_TKE = TKE_reduc*tot_TKE   ! = tot_TKE + TKE_forcing(k)
          mech_TKE = TKE_reduc*mech_TKE
          conv_PErel = TKE_reduc*conv_PErel
        endif
      endif

      ! Precalculate some temporary expressions that are independent of Kddt_h(K).
      if (CS%orig_PE_calc) then
        if (K==2) then
          dTe_t2 = 0.0 ; dSe_t2 = 0.0
        else
          dTe_t2 = Kddt_h(K-1) * ((T0(k-2) - T0(k-1)) + dTe(k-2))
          dSe_t2 = Kddt_h(K-1) * ((S0(k-2) - S0(k-1)) + dSe(k-2))
        endif
      endif
      dt_h = dt / max(0.5*(dz(k-1)+dz(k)), 1e-15*dz_sum)

      !   This tests whether the layers above and below this interface are in
      ! a convectively stable configuration, without considering any effects of
      ! mixing at higher interfaces.  It is an approximation to the more
      ! complete test dPEc_dKd_Kd0 >= 0.0, that would include the effects of
      ! mixing across interface K-1.  The dT_to_dColHt here are effectively
      ! mass-weighted estimates of dSV_dT.
      Convectively_stable = ( 0.0 <= &
           ( (dT_to_dColHt(k) + dT_to_dColHt(k-1) ) * (T0(k-1)-T0(k)) + &
             (dS_to_dColHt(k) + dS_to_dColHt(k-1) ) * (S0(k-1)-S0(k)) ) )

      if ((mech_TKE + conv_PErel) <= 0.0 .and. Convectively_stable) then
        ! Energy is already exhausted, so set Kd = 0 and cycle or exit?
        tot_TKE = 0.0 ; mech_TKE = 0.0 ; conv_PErel = 0.0
        Kd(K) = 0.0 ; Kddt_h(K) = 0.0
        sfc_disconnect = .true.
        ! if (.not.debug) exit

       !   The estimated properties for layer k-1 can be calculated, using
       ! greatly simplified expressions when Kddt_h = 0.  This enables the
       ! tridiagonal solver for the whole column to be completed for debugging
       ! purposes, and also allows for something akin to convective adjustment
       ! in unstable interior regions?
        b1 = 1.0 / hp_a(k-1)
        c1(K) = 0.0
        if (CS%orig_PE_calc) then
          dTe(k-1) = b1 * ( dTe_t2 )
          dSe(k-1) = b1 * ( dSe_t2 )
        endif

        hp_a(k) = h(k)
        dT_to_dPE_a(k) = dT_to_dPE(k)
        dS_to_dPE_a(k) = dS_to_dPE(k)
        dT_to_dColHt_a(k) = dT_to_dColHt(k)
        dS_to_dColHt_a(k) = dS_to_dColHt(k)

      else ! tot_TKE > 0.0 or this is a potentially convectively unstable profile.
        sfc_disconnect = .false.

        ! Precalculate some more temporary expressions that are independent of
        ! Kddt_h(K).
        if (CS%orig_PE_calc) then
          if (K==2) then
            dT_km1_t2 = (T0(k)-T0(k-1))
            dS_km1_t2 = (S0(k)-S0(k-1))
          else
            dT_km1_t2 = (T0(k)-T0(k-1)) - &
                  (Kddt_h(K-1) / hp_a(k-1)) * ((T0(k-2) - T0(k-1)) + dTe(k-2))
            dS_km1_t2 = (S0(k)-S0(k-1)) - &
                  (Kddt_h(K-1) / hp_a(k-1)) * ((S0(k-2) - S0(k-1)) + dSe(k-2))
          endif
          dTe_term = dTe_t2 + hp_a(k-1) * (T0(k-1)-T0(k))
          dSe_term = dSe_t2 + hp_a(k-1) * (S0(k-1)-S0(k))
        else
          if (K<=2) then
            Th_a(k-1) = h(k-1) * T0(k-1) ; Sh_a(k-1) = h(k-1) * S0(k-1)
          else
            Th_a(k-1) = h(k-1) * T0(k-1) + Kddt_h(K-1) * Te(k-2)
            Sh_a(k-1) = h(k-1) * S0(k-1) + Kddt_h(K-1) * Se(k-2)
          endif
          Th_b(k) = h(k) * T0(k) ; Sh_b(k) = h(k) * S0(k)
        endif

        !   Using Pr=1 and the diffusivity at the bottom interface (once it is
        ! known), determine how much resolved mean kinetic energy (MKE) will be
        ! extracted within a timestep and add a fraction CS%MKE_to_TKE_effic of
        ! this to the mTKE budget available for mixing in the next layer.

        if ((CS%MKE_to_TKE_effic > 0.0) .and. (htot*h(k) > 0.0)) then
          ! This is the energy that would be available from homogenizing the
          ! velocities between layer k and the layers above.
          dMKE_max = (GV%H_to_RZ * CS%MKE_to_TKE_effic) * 0.5 * &
              (h(k) / ((htot + h(k))*htot)) * &
              (((uhtot-u(k)*htot)**2) + ((vhtot-v(k)*htot)**2))
          ! A fraction (1-exp(Kddt_h*MKE2_Hharm)) of this energy would be
          ! extracted by mixing with a finite viscosity.
          MKE2_Hharm = (htot + h(k) + 2.0*h_neglect) / &
                       ((htot+h_neglect) * (h(k)+h_neglect))
        else
          dMKE_max = 0.0
          MKE2_Hharm = 0.0
        endif

        ! At this point, Kddt_h(K) will be unknown because its value may depend
        ! on how much energy is available.  mech_TKE might be negative due to
        ! contributions from TKE_forced.
        dz_tt = dztot + dz_tt_min
        TKE_here = mech_TKE + CS%wstar_ustar_coef*conv_PErel
        if (TKE_here > 0.0) then
          if (CS%answer_date < 20240101) then
            if (CS%wT_scheme==wT_from_cRoot_TKE) then
              vstar = CS%vstar_scale_fac * vstar_unit_scale * (SpV_dt(K)*TKE_here)**C1_3
            elseif (CS%wT_scheme==wT_from_RH18) then
              Surface_Scale = max(0.05, 1.0 - dztot / MLD_guess)
              vstar = CS%vstar_scale_fac * Surface_Scale * (CS%vstar_surf_fac*u_star + &
                        vstar_unit_scale * (CS%wstar_ustar_coef*conv_PErel*SpV_dt(K))**C1_3)
            endif
          else
            if (CS%wT_scheme==wT_from_cRoot_TKE) then
              vstar = CS%vstar_scale_fac * cuberoot(SpV_dt(K)*TKE_here)
            elseif (CS%wT_scheme==wT_from_RH18) then
              Surface_Scale = max(0.05, 1.0 - dztot / MLD_guess)
              vstar = (CS%vstar_scale_fac * Surface_Scale) * ( CS%vstar_surf_fac*u_star + &
                        cuberoot((CS%wstar_ustar_coef*conv_PErel) * SpV_dt(K)) )
            endif
          endif
          hbs_here = min(hb_hs(K), MixLen_shape(K))
          mixlen(K) = max(CS%min_mix_len, ((dz_tt*hbs_here)*vstar) / &
              ((CS%Ekman_scale_coef * absf) * (dz_tt*hbs_here) + vstar))
          !Note setting Kd_guess0 to vstar * CS%vonKar * mixlen(K) here will
          ! change the answers.  Therefore, skipping that.
          if (.not.CS%Use_MLD_iteration) then
            Kd_guess0 = (h_dz_int(K)*vstar) * CS%vonKar * ((dz_tt*hbs_here)*vstar) / &
              ((CS%Ekman_scale_coef * absf) * (dz_tt*hbs_here) + vstar)
          elseif (CS%eqdisc) then  ! ML-eqdisc line1/2
            Kd_guess0 = MixLen_shape(K) * v0_ML_turb_vel_scale * MLD_guess ! ML-eqdisc
          else
            Kd_guess0 = (h_dz_int(K)*vstar) * CS%vonKar * mixlen(K)
          endif
        else
          vstar = 0.0 ; Kd_guess0 = 0.0
        endif
        mixvel(K) = vstar ! Track vstar
        Kddt_h_g0 = Kd_guess0 * dt_h

        if (no_MKE_conversion) then
          ! Without conversion from MKE to TKE, the updated diffusivity can be determined directly.
          ! Replace h(k) with hp_b(k) = h(k), and dT_to_dPE with dT_to_dPE_b, etc., for a 2-direction solver.
          call find_Kd_from_PE_chg(0.0, Kd_guess0, dt_h, tot_TKE, hp_a(k-1), h(k), &
                      Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                      dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE(k), dS_to_dPE(k), &
                      pres_Z(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                      dT_to_dColHt(k), dS_to_dColHt(k), Kd_add=Kd(K), PE_chg=TKE_used, &
                      dPE_max=PE_chg_max, frac_dKd_max_PE=frac_in_BL)
          convectively_unstable = (PE_chg_max < 0.0)
          PE_chg_g0 = TKE_used  ! This is only used in the convectively unstable limit.
          MKE_src = 0.0
        elseif (CS%orig_PE_calc) then
          call find_PE_chg_orig(Kddt_h_g0, h(k), hp_a(k-1), dTe_term, dSe_term, &
                   dT_km1_t2, dS_km1_t2, dT_to_dPE(k), dS_to_dPE(k), &
                   dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), &
                   pres_Z(K), dT_to_dColHt(k), dS_to_dColHt(k), &
                   dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                   PE_chg=PE_chg_g0, dPE_max=PE_chg_max, dPEc_dKd_0=dPEc_dKd_Kd0 )
          convectively_unstable =  (PE_chg_g0 < 0.0) .or. ((vstar == 0.0) .and. (dPEc_dKd_Kd0 < 0.0))
          MKE_src = dMKE_max*(1.0 - exp(-Kddt_h_g0 * MKE2_Hharm))
        else
          call find_PE_chg(0.0, Kddt_h_g0, hp_a(k-1), h(k), &
                   Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                   dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE(k), dS_to_dPE(k), &
                   pres_Z(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                   dT_to_dColHt(k), dS_to_dColHt(k), &
                   PE_chg=PE_chg_g0, dPE_max=PE_chg_max, dPEc_dKd_0=dPEc_dKd_Kd0 )
          convectively_unstable =  (PE_chg_g0 < 0.0) .or. ((vstar == 0.0) .and. (dPEc_dKd_Kd0 < 0.0))
          MKE_src = dMKE_max*(1.0 - exp(-Kddt_h_g0 * MKE2_Hharm))
        endif

        ! This block checks out different cases to determine Kd at the present interface.
        if (convectively_unstable) then
          ! This column is convectively unstable.
          if (PE_chg_max <= 0.0) then
            ! Does MKE_src need to be included in the calculation of vstar here?
            TKE_here = mech_TKE + CS%wstar_ustar_coef*(conv_PErel-PE_chg_max)
            if (TKE_here > 0.0) then
              if (CS%answer_date < 20240101) then
                if (CS%wT_scheme==wT_from_cRoot_TKE) then
                  vstar = CS%vstar_scale_fac * vstar_unit_scale * (SpV_dt(K)*TKE_here)**C1_3
                elseif (CS%wT_scheme==wT_from_RH18) then
                  Surface_Scale = max(0.05, 1. - dztot / MLD_guess)
                  vstar = CS%vstar_scale_fac * Surface_Scale * (CS%vstar_surf_fac*u_star + &
                                  vstar_unit_scale * (CS%wstar_ustar_coef*conv_PErel*SpV_dt(K))**C1_3)
                endif
              else
                if (CS%wT_scheme==wT_from_cRoot_TKE) then
                  vstar = CS%vstar_scale_fac * cuberoot(SpV_dt(K)*TKE_here)
                elseif (CS%wT_scheme==wT_from_RH18) then
                  Surface_Scale = max(0.05, 1. - dztot / MLD_guess)
                  vstar = (CS%vstar_scale_fac * Surface_Scale) * ( CS%vstar_surf_fac*u_star + &
                                  cuberoot((CS%wstar_ustar_coef*conv_PErel) * SpV_dt(K)) )
                endif
              endif
              hbs_here = min(hb_hs(K), MixLen_shape(K))
              mixlen(K) = max(CS%min_mix_len, ((dz_tt*hbs_here)*vstar) / &
                  ((CS%Ekman_scale_coef * absf) * (dz_tt*hbs_here) + vstar))
              if (.not.CS%Use_MLD_iteration) then
              ! Note again (as prev) that using mixlen here
              !  instead of redoing the computation will change answers...
                Kd(K) = (h_dz_int(K)*vstar) * CS%vonKar *  ((dz_tt*hbs_here)*vstar) / &
                      ((CS%Ekman_scale_coef * absf) * (dz_tt*hbs_here) + vstar)
              elseif (CS%eqdisc)  then  ! ML-eqdisc line2/2
                Kd(K) = MixLen_shape(K) * v0_ML_turb_vel_scale * MLD_guess ! ML-eqdisc
              else
                Kd(K) = (h_dz_int(K)*vstar) * CS%vonKar * mixlen(K)
              endif
            else
              vstar = 0.0 ; Kd(K) = 0.0
            endif
            mixvel(K) = vstar

            if (CS%orig_PE_calc) then
              call find_PE_chg_orig(Kd(K)*dt_h, h(k), hp_a(k-1), dTe_term, dSe_term, &
                       dT_km1_t2, dS_km1_t2, dT_to_dPE(k), dS_to_dPE(k), &
                       dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), &
                       pres_Z(K), dT_to_dColHt(k), dS_to_dColHt(k), &
                       dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                       PE_chg=dPE_conv)
            else
              call find_PE_chg(0.0, Kd(K)*dt_h, hp_a(k-1), h(k), &
                       Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                       dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE(k), dS_to_dPE(k), &
                       pres_Z(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                       dT_to_dColHt(k), dS_to_dColHt(k), &
                       PE_chg=dPE_conv)
            endif
            ! Should this be iterated to convergence for Kd?
            if (dPE_conv > 0.0) then
              Kd(K) = Kd_guess0 ; dPE_conv = PE_chg_g0
            else
              MKE_src = dMKE_max*(1.0 - exp(-(Kd(K)*dt_h) * MKE2_Hharm))
            endif
          else
            ! The energy change does not vary monotonically with Kddt_h.  Find the maximum?
            Kd(K) = Kd_guess0 ; dPE_conv = PE_chg_g0
          endif

          conv_PErel = conv_PErel - dPE_conv
          mech_TKE = mech_TKE + MKE_src
          if (CS%TKE_diagnostics) then
            eCD%dTKE_conv = eCD%dTKE_conv - CS%nstar*dPE_conv * I_dtdiag
            eCD%dTKE_MKE = eCD%dTKE_MKE + MKE_src * I_dtdiag
          endif
          if (sfc_connected) then
            MLD_output = MLD_output + dz(k)
          endif

          Kddt_h(K) = Kd(K) * dt_h

        elseif (no_MKE_conversion) then  ! (PE_chg_max >= 0.0) and use the diffusivity from find_Kd_from_PE_chg.
          ! Kd(K) and TKE_used were already set by find_Kd_from_PE_chg.

          ! frac_in_BL = min((TKE_used / PE_chg_g0), 1.0)
          if (sfc_connected)  MLD_output = MLD_output + frac_in_BL*dz(k)
          if (frac_in_BL < 1.0) sfc_disconnect = .true.

          ! Reduce the mechanical and convective TKE proportionately.
          TKE_reduc = 0.0   ! tot_TKE could be 0 if Convectively_stable is false.
          if ((tot_TKE > 0.0) .and. (tot_TKE > TKE_used)) TKE_reduc = (tot_TKE - TKE_used) / tot_TKE

          ! All TKE should have been consumed.
          if (CS%TKE_diagnostics) then
            eCD%dTKE_mixing = eCD%dTKE_mixing - TKE_used * I_dtdiag
            eCD%dTKE_conv_decay = eCD%dTKE_conv_decay + &
                (1.0-TKE_reduc)*(CS%nstar-nstar_FC) * conv_PErel * I_dtdiag
          endif

          tot_TKE = tot_TKE - TKE_used
          mech_TKE = TKE_reduc*mech_TKE
          conv_PErel = TKE_reduc*conv_PErel

          Kddt_h(K) = Kd(K) * dt_h

        elseif (tot_TKE + (MKE_src - PE_chg_g0) >= 0.0) then
          ! This column is convectively stable and there is energy to support the suggested
          ! mixing.  Keep that estimate.
          Kd(K) = Kd_guess0
          Kddt_h(K) = Kddt_h_g0

          ! Reduce the mechanical and convective TKE proportionately.
          tot_TKE = tot_TKE + MKE_src
          TKE_reduc = 0.0   ! tot_TKE could be 0 if Convectively_stable is false.
          if (tot_TKE > 0.0) TKE_reduc = (tot_TKE - PE_chg_g0) / tot_TKE
          if (CS%TKE_diagnostics) then
            eCD%dTKE_mixing = eCD%dTKE_mixing - PE_chg_g0 * I_dtdiag
            eCD%dTKE_MKE = eCD%dTKE_MKE + MKE_src * I_dtdiag
            eCD%dTKE_conv_decay = eCD%dTKE_conv_decay + &
                (1.0-TKE_reduc)*(CS%nstar-nstar_FC) * conv_PErel * I_dtdiag
          endif
          tot_TKE = TKE_reduc*tot_TKE
          mech_TKE = TKE_reduc*(mech_TKE + MKE_src)
          conv_PErel = TKE_reduc*conv_PErel
          if (sfc_connected) then
            MLD_output = MLD_output + dz(k)
          endif

        elseif (tot_TKE == 0.0) then
          ! This can arise if nstar_FC = 0, but it is not common.
          Kd(K) = 0.0 ; Kddt_h(K) = 0.0
          tot_TKE = 0.0 ; conv_PErel = 0.0 ; mech_TKE = 0.0
          sfc_disconnect = .true.
        else
          ! There is not enough energy to support the mixing, so reduce the
          ! diffusivity to what can be supported.
          Kddt_h_max = Kddt_h_g0 ; Kddt_h_min = 0.0
          TKE_left_max = tot_TKE + (MKE_src - PE_chg_g0)
          TKE_left_min = tot_TKE

          ! As a starting guess, take the minimum of a false position estimate
          ! and a Newton's method estimate starting from Kddt_h = 0.0.
          Kddt_h_guess = tot_TKE * Kddt_h_max / max( PE_chg_g0  - MKE_src, &
                           Kddt_h_max * (dPEc_dKd_Kd0 - dMKE_max * MKE2_Hharm) )
          ! The above expression is mathematically the same as the following
          ! except it is not susceptible to division by zero when
          !   dPEc_dKd_Kd0 = dMKE_max = 0 .
          !  Kddt_h_guess = tot_TKE * min( Kddt_h_max / (PE_chg_g0 - MKE_src), &
          !                      1.0 / (dPEc_dKd_Kd0 - dMKE_max * MKE2_Hharm) )
          if (debug) then
            TKE_left_itt(:) = 0.0 ; dPEa_dKd_itt(:) = 0.0 ; PE_chg_itt(:) = 0.0
            MKE_src_itt(:) = 0.0 ; Kddt_h_itt(:) = 0.0
          endif
          do itt=1,max_itt
            if (CS%orig_PE_calc) then
              call find_PE_chg_orig(Kddt_h_guess, h(k), hp_a(k-1), dTe_term, dSe_term, &
                       dT_km1_t2, dS_km1_t2, dT_to_dPE(k), dS_to_dPE(k), &
                       dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), &
                       pres_Z(K), dT_to_dColHt(k), dS_to_dColHt(k), &
                       dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                       PE_chg=PE_chg, dPEc_dKd=dPEc_dKd )
            else
              call find_PE_chg(0.0, Kddt_h_guess, hp_a(k-1), h(k), &
                       Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                       dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE(k), dS_to_dPE(k), &
                       pres_Z(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                       dT_to_dColHt(k), dS_to_dColHt(k), &
                       PE_chg=PE_chg, dPEc_dKd=dPEc_dKd)
            endif
            MKE_src = dMKE_max * (1.0 - exp(-MKE2_Hharm * Kddt_h_guess))
            dMKE_src_dK = dMKE_max * MKE2_Hharm * exp(-MKE2_Hharm * Kddt_h_guess)

            TKE_left = tot_TKE + (MKE_src - PE_chg)
            if (debug .and. itt<=20) then
              Kddt_h_itt(itt) = Kddt_h_guess ; MKE_src_itt(itt) = MKE_src
              PE_chg_itt(itt) = PE_chg ; dPEa_dKd_itt(itt) = dPEc_dKd
              TKE_left_itt(itt) = TKE_left
            endif
            ! Store the new bounding values, bearing in mind that min and max
            ! here refer to Kddt_h and dTKE_left/dKddt_h < 0:
            if (TKE_left >= 0.0) then
              Kddt_h_min = Kddt_h_guess ; TKE_left_min = TKE_left
            else
              Kddt_h_max = Kddt_h_guess ; TKE_left_max = TKE_left
            endif

            ! Try to use Newton's method, but if it would go outside the bracketed
            ! values use the false-position method instead.
            use_Newt = .true.
            if (dPEc_dKd - dMKE_src_dK <= 0.0) then
              use_Newt = .false.
            else
              dKddt_h_Newt = TKE_left / (dPEc_dKd - dMKE_src_dK)
              Kddt_h_Newt = Kddt_h_guess + dKddt_h_Newt
              if ((Kddt_h_Newt > Kddt_h_max) .or. (Kddt_h_Newt < Kddt_h_min)) &
                use_Newt = .false.
            endif

            if (use_Newt) then
              Kddt_h_next = Kddt_h_guess + dKddt_h_Newt
              dKddt_h = dKddt_h_Newt
            else
              Kddt_h_next = (TKE_left_max * Kddt_h_min - Kddt_h_max * TKE_left_min) / &
                            (TKE_left_max - TKE_left_min)
              dKddt_h = Kddt_h_next - Kddt_h_guess
            endif

            if ((abs(dKddt_h) < 1e-9*Kddt_h_guess) .or. (itt==max_itt)) then
              ! Use the old value so that the energy calculation does not need to be repeated.
              if (debug) num_itts(K) = itt
              exit
            else
              Kddt_h_guess = Kddt_h_next
            endif
          enddo ! Inner iteration loop on itt.
          Kd(K) = Kddt_h_guess / dt_h ; Kddt_h(K) = Kd(K) * dt_h

          ! All TKE should have been consumed.
          if (CS%TKE_diagnostics) then
            eCD%dTKE_mixing = eCD%dTKE_mixing - (tot_TKE + MKE_src) * I_dtdiag
            eCD%dTKE_MKE = eCD%dTKE_MKE + MKE_src * I_dtdiag
            eCD%dTKE_conv_decay = eCD%dTKE_conv_decay + &
                (CS%nstar-nstar_FC) * conv_PErel * I_dtdiag
          endif

          if (sfc_connected) MLD_output = MLD_output + (PE_chg / (PE_chg_g0)) * dz(k)

          tot_TKE = 0.0 ; mech_TKE = 0.0 ; conv_PErel = 0.0
          sfc_disconnect = .true.
        endif ! End of convective or forced mixing cases to determine Kd.

        Kddt_h(K) = Kd(K) * dt_h
        !   At this point, the final value of Kddt_h(K) is known, so the
        ! estimated properties for layer k-1 can be calculated.
        b1 = 1.0 / (hp_a(k-1) + Kddt_h(K))
        c1(K) = Kddt_h(K) * b1
        if (CS%orig_PE_calc) then
          dTe(k-1) = b1 * ( Kddt_h(K)*(T0(k)-T0(k-1)) + dTe_t2 )
          dSe(k-1) = b1 * ( Kddt_h(K)*(S0(k)-S0(k-1)) + dSe_t2 )
        endif

        hp_a(k) = h(k) + (hp_a(k-1) * b1) * Kddt_h(K)
        dT_to_dPE_a(k) = dT_to_dPE(k) + c1(K)*dT_to_dPE_a(k-1)
        dS_to_dPE_a(k) = dS_to_dPE(k) + c1(K)*dS_to_dPE_a(k-1)
        dT_to_dColHt_a(k) = dT_to_dColHt(k) + c1(K)*dT_to_dColHt_a(k-1)
        dS_to_dColHt_a(k) = dS_to_dColHt(k) + c1(K)*dS_to_dColHt_a(k-1)

      endif  ! tot_TKT > 0.0 branch.  Kddt_h(K) has been set.

      ! Store integrated velocities and thicknesses for MKE conversion calculations.
      if (sfc_disconnect) then
        ! There is no turbulence at this interface, so zero out the running sums.
        uhtot = u(k)*h(k)
        vhtot = v(k)*h(k)
        htot  = h(k)
        dztot = dz(k)
        sfc_connected = .false.
      else
        uhtot = uhtot + u(k)*h(k)
        vhtot = vhtot + v(k)*h(k)
        htot  = htot + h(k)
        dztot = dztot + dz(k)
      endif

      if (calc_Te) then
        if (K==2) then
          Te(1) = b1*(h(1)*T0(1))
          Se(1) = b1*(h(1)*S0(1))
        else
          Te(k-1) = b1 * (h(k-1) * T0(k-1) + Kddt_h(K-1) * Te(k-2))
          Se(k-1) = b1 * (h(k-1) * S0(k-1) + Kddt_h(K-1) * Se(k-2))
        endif
      endif
    enddo
    Kd(nz+1) = 0.0

    if (debug) then
      ! Complete the tridiagonal solve for Te.
      b1 = 1.0 / hp_a(nz)
      Te(nz) = b1 * (h(nz) * T0(nz) + Kddt_h(nz) * Te(nz-1))
      Se(nz) = b1 * (h(nz) * S0(nz) + Kddt_h(nz) * Se(nz-1))
      dT_expect(nz) = Te(nz) - T0(nz) ; dS_expect(nz) = Se(nz) - S0(nz)
      do k=nz-1,1,-1
        Te(k) = Te(k) + c1(K+1)*Te(k+1)
        Se(k) = Se(k) + c1(K+1)*Se(k+1)
        dT_expect(k) = Te(k) - T0(k) ; dS_expect(k) = Se(k) - S0(k)
      enddo
    endif

    if (debug) then
      dPE_debug = 0.0
      do k=1,nz
        dPE_debug = dPE_debug + (dT_to_dPE(k) * (Te(k) - T0(k)) + &
                                 dS_to_dPE(k) * (Se(k) - S0(k)))
      enddo
      mixing_debug = dPE_debug * I_dtdiag
    endif

    if (OBL_it >= CS%Max_MLD_Its) exit

    ! The following lines are used for the iteration.  Note the iteration has been altered
    ! to use the value predicted by the TKE threshold (ML_depth).  This is because the mstar
    ! is now dependent on the ML, and therefore the ML needs to be estimated more precisely
    ! than the grid spacing.

    ! New method uses ML_depth as computed in ePBL routine
    MLD_found = MLD_output

    ! Find out whether to do another iteration and the new bounds on it.
    if (CS%MLD_iter_bug) then
      ! There is a bug in the logic here if (MLD_found - MLD_guess == CS%MLD_tol).
      if (MLD_found - MLD_guess > CS%MLD_tol) then
        min_MLD = MLD_guess ; dMLD_min = MLD_found - MLD_guess
      elseif (abs(MLD_found - MLD_guess) < CS%MLD_tol) then
        exit ! Break the MLD convergence loop
      else ! We know this guess was too deep
        max_MLD = MLD_guess ; dMLD_max = MLD_found - MLD_guess ! < -CS%MLD_tol
      endif
    else
      if (abs(MLD_found - MLD_guess) < CS%MLD_tol) then
        exit ! Break the MLD convergence loop
      elseif (MLD_found > MLD_guess) then  ! This guess was too shallow
        min_MLD = MLD_guess ; dMLD_min = MLD_found - MLD_guess
      else ! We know this guess was too deep
        max_MLD = MLD_guess ; dMLD_max = MLD_found - MLD_guess ! < -CS%MLD_tol
      endif
    endif

    if (OBL_it < CS%Max_MLD_Its) then
      if (CS%MLD_bisection) then
        ! For the next pass, guess the average of the minimum and maximum values.
        MLD_guess = 0.5*(min_MLD + max_MLD)
      else ! Try using the false position method or the returned value instead of simple bisection.
        ! Taking the occasional step with MLD_output empirically helps to converge faster.
        if ((dMLD_min > 0.0) .and. (dMLD_max < 0.0) .and. (OBL_it > 2) .and. (mod(OBL_it-1,4) > 0)) then
          ! Both bounds have valid change estimates and are probably in the range of possible outputs.
          MLD_guess = (dMLD_min*max_MLD - dMLD_max*min_MLD) / (dMLD_min - dMLD_max)
        elseif ((MLD_found > min_MLD) .and. (MLD_found < max_MLD)) then
          ! The output MLD_found is an interesting guess, as it is likely to bracket the true solution
          ! along with the previous value of MLD_guess and to be close to the solution.
          MLD_guess = MLD_found
        else ! Bisect if the other guesses would be out-of-bounds.  This does not happen much.
          MLD_guess = 0.5*(min_MLD + max_MLD)
        endif
      endif
    endif

  enddo ! Iteration loop for converged boundary layer thickness.

  eCD%OBL_its = min(OBL_it, CS%Max_MLD_Its)
  if (CS%Use_LT) then
    eCD%LA = LA ; eCD%LAmod = LAmod ; eCD%mstar = mstar_total ; eCD%mstar_LT = mstar_LT
  else
    eCD%LA = 0.0 ; eCD%LAmod = 0.0 ; eCD%mstar = mstar_total ; eCD%mstar_LT = 0.0
  endif

  MLD_io = MLD_output

end subroutine ePBL_column


!> This subroutine determines the diffusivities from a bottom boundary layer version of
!! the integrated energetics mixed layer model for a single column of water.
subroutine ePBL_BBL_column(h, dz, u, v, T0, S0, dSV_dT, dSV_dS, SpV_dt, absf, &
                           dt, Kd, BBL_TKE_in, u_star_BBL, u_star_BBL_z_t, b_flux_BBL, Kd_BBL, BBLD_io, mixvel_BBL, &
                           mixlen_BBL, GV, US, CS, eCD)
  type(verticalGrid_type),   intent(in)  :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZK_(GV)), intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZK_(GV)), intent(in)  :: dz     !< The vertical distance across layers [Z ~> m].
  real, dimension(SZK_(GV)), intent(in)  :: u      !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZK_(GV)), intent(in)  :: v      !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZK_(GV)), intent(in)  :: T0     !< The initial layer temperatures [C ~> degC].
  real, dimension(SZK_(GV)), intent(in)  :: S0     !< The initial layer salinities [S ~> ppt].

  real, dimension(SZK_(GV)), intent(in)  :: dSV_dT !< The partial derivative of in-situ specific
                                                   !! volume with potential temperature
                                                   !! [R-1 C-1 ~> m3 kg-1 degC-1].
  real, dimension(SZK_(GV)), intent(in)  :: dSV_dS !< The partial derivative of in-situ specific
                                                   !! volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZK_(GV)+1), intent(in) :: SpV_dt !< Specific volume interpolated to interfaces
                                                   !! divided by dt (if non-Boussinesq) or
                                                   !! 1.0 / (dt * Rho0), in [R-1 T-1 ~> m3 kg-1 s-1],
                                                   !! used to convert local TKE into a turbulence
                                                   !! velocity cubed.
  real,                    intent(in)    :: absf   !< The absolute value of the Coriolis parameter [T-1 ~> s-1].
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  real, dimension(SZK_(GV)+1), &
                           intent(in)    :: Kd     !< The diffusivities at interfaces due to previously
                                                   !! applied mixing processes [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real,                    intent(in)    :: BBL_TKE_in !< The mechanically generated turbulent
                                                   !! kinetic energy available for bottom boundary
                                                   !! layer mixing within a time step [R Z3 T-2 ~> J m-2].
  real,                    intent(in)    :: u_star_BBL !< The bottom boundary layer friction velocity
                                                       !! in thickness flux units [H T-1 ~> m s-1 or kg m-2 s-1]
  real,                    intent(in)    :: u_star_BBL_z_t !< The bottom boundary layer friction velocity
                                                       !! converted to length flux units [Z T-1 ~> m s-1]
  real,                    intent(in)    :: b_flux_BBL !< The bottom boundary layer buoyancy flux
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: Kd_BBL !< The bottom boundary layer contribution to diffusivities
                                                   !! at interfaces [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real,                    intent(inout) :: BBLD_io !< A first guess at the bottom boundary layer depth on input, and
                                                   !! the calculated bottom boundary layer depth on output [Z ~> m]
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: mixvel_BBL !< The profile of boundary layer turbulent mixing
                                                   !! velocities [Z T-1 ~> m s-1]
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: mixlen_BBL !< The profile of bottom boundary layer turbulent
                                                   !! mixing lengths [Z ~> m]
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(energetic_PBL_CS),  intent(in)    :: CS     !< Energetic PBL control structure
  type(ePBL_column_diags), intent(inout) :: eCD    !< A container for passing around diagnostics.

!    This subroutine determines the contributions from diffusivities in a single column from a
!  bottom-boundary layer adaptation of the integrated energetics planetary boundary layer (ePBL)
!  model.  It accounts for the possibility that the surface boundary diffusivities have already
!  been determined.  All calculations are done implicitly, and there is no stability limit on the
!  time step.  Only mechanical mixing in the bottom boundary layer is considered.  (Geothermal heat
!  fluxes are addressed elsewhere in the MOM6 code, and convection throughout the water column is
!  handled by the surface version of ePBL.)  There is no conversion of released mean kinetic energy
!  into bottom boundary layer turbulent kinetic energy (at least for now), apart from the explicit
!  energy that is supplied as an argument to this routine.

  ! Local variables
  real, dimension(SZK_(GV)+1) :: &
    pres_Z, &       ! Interface pressures with a rescaling factor to convert interface height
                    ! movements into changes in column potential energy [R Z2 T-2 ~> kg m-1 s-2].
    dztop_dztot     ! The distance from the surface divided by the thickness of the
                    ! water column [nondim].
  real :: mech_BBL_TKE ! The mechanically generated turbulent kinetic energy available for
                    ! bottom boundary layer mixing within a time step [R Z3 T-2 ~> J m-2].
  real :: TKE_eff_avail ! The turbulent kinetic energy that is effectively available to drive mixing
                    ! after any effects of exponentially decay have been taken into account
                    ! [R Z3 T-2 ~> J m-2]
  real :: TKE_eff_used ! The amount of TKE_eff_avail that has been used to drive mixing [R Z3 T-2 ~> J m-2]
  real :: htot      ! The total thickness of the layers above an interface [H ~> m or kg m-2].
  real :: dztot     ! The total depth of the layers above an interface [Z ~> m].
  real :: uhtot     ! The depth integrated zonal velocities in the layers above [H L T-1 ~> m2 s-1 or kg m-1 s-1]
  real :: vhtot     ! The depth integrated meridional velocities in the layers above [H L T-1 ~> m2 s-1 or kg m-1 s-1]
  real :: Idecay_len_TKE  ! The inverse of a turbulence decay length scale [H-1 ~> m-1 or m2 kg-1].
  real :: dz_sum    ! The total thickness of the water column [Z ~> m].

  real, dimension(SZK_(GV)) :: &
    dT_to_dColHt, & ! Partial derivative of the total column height with the temperature changes
                    ! within a layer [Z C-1 ~> m degC-1].
    dS_to_dColHt, & ! Partial derivative of the total column height with the salinity changes
                    ! within a layer  [Z S-1 ~> m ppt-1].
    dT_to_dPE, &    ! Partial derivatives of column potential energy with the temperature
                    ! changes within a layer, in [R Z3 T-2 C-1 ~> J m-2 degC-1].
    dS_to_dPE, &    ! Partial derivatives of column potential energy with the salinity changes
                    ! within a layer, in [R Z3 T-2 S-1 ~> J m-2 ppt-1].
    dT_to_dColHt_a, & ! Partial derivative of the total column height with the temperature changes
                    ! within a layer, including the implicit effects  of mixing with layers higher
                    ! in the water column [Z C-1 ~> m degC-1].
    dS_to_dColHt_a, & ! Partial derivative of the total column height with the salinity changes
                    ! within a layer, including the implicit effects  of mixing with layers higher
                    ! in the water column [Z S-1 ~> m ppt-1].
    dT_to_dPE_a, &  ! Partial derivatives of column potential energy with the temperature changes
                    ! within a layer, including the implicit effects of mixing with layers higher
                    ! in the water column [R Z3 T-2 C-1 ~> J m-2 degC-1].
    dS_to_dPE_a, &  ! Partial derivative of column potential energy with the salinity changes
                    ! within a layer, including the implicit effects of mixing with layers higher
                    ! in the water column [R Z3 T-2 S-1 ~> J m-2 ppt-1].
    dT_to_dColHt_b, & ! Partial derivative of the total column height with the temperature changes
                    ! within a layer, including the implicit effects of mixing with layers deeper
                    ! in the water column [Z C-1 ~> m degC-1].
    dS_to_dColHt_b, & ! Partial derivative of the total column height with the salinity changes
                    ! within a layer, including the implicit effects of mixing with layers deeper
                    ! in the water column [Z S-1 ~> m ppt-1].
    dT_to_dPE_b, &  ! Partial derivatives of column potential energy with the temperature changes
                    ! within a layer, including the implicit effects of mixing with layers deeper
                    ! in the water column [R Z3 T-2 C-1 ~> J m-2 degC-1].
    dS_to_dPE_b, &  ! Partial derivative of column potential energy with the salinity changes
                    ! within a layer, including the implicit effects of mixing with layers deeper
                    ! in the water column [R Z3 T-2 S-1 ~> J m-2 ppt-1].
    c1, &           ! c1 is used by the tridiagonal solver [nondim].
    Te, &           ! Estimated final values of T in the column [C ~> degC].
    Se, &           ! Estimated final values of S in the column [S ~> ppt].
    dTe, &          ! Running (1-way) estimates of temperature change [C ~> degC].
    dSe, &          ! Running (1-way) estimates of salinity change [S ~> ppt].
    hp_a, &         ! An effective pivot thickness of the layer including the effects
                    ! of coupling with layers above [H ~> m or kg m-2].  This is the first term
                    ! in the denominator of b1 in a downward-oriented tridiagonal solver.
    hp_b, &         ! An effective pivot thickness of the layer including the effects
                    ! of coupling with layers below [H ~> m or kg m-2].  This is the first term
                    ! in the denominator of b1 in an upward-oriented tridiagonal solver.
    Th_a, &         ! An effective temperature times a thickness in the layer above, including implicit
                    ! mixing effects with other yet higher layers [C H ~> degC m or degC kg m-2].
    Sh_a, &         ! An effective salinity times a thickness in the layer above, including implicit
                    ! mixing effects with other yet higher layers [S H ~> ppt m or ppt kg m-2].
    Th_b, &         ! An effective temperature times a thickness in the layer below, including implicit
                    ! mixing effects with other yet lower layers [C H ~> degC m or degC kg m-2].
    Sh_b            ! An effective salinity times a thickness in the layer below, including implicit
                    ! mixing effects with other yet lower layers [S H ~> ppt m or ppt kg m-2].
  real, dimension(SZK_(GV)+1) :: &
    MixLen_shape, & ! A nondimensional shape factor for the mixing length that
                    ! gives it an appropriate asymptotic value at the bottom of
                    ! the boundary layer [nondim].
    h_dz_int, &     ! The ratio of the layer thicknesses over the vertical distances
                    ! across the layers surrounding an interface [H Z-1 ~> nondim or kg m-3]
    Kddt_h          ! The total diapycnal diffusivity at an interface times a timestep divided by the
                    ! average thicknesses around an interface [H ~> m or kg m-2].
  real :: b1        ! b1 is inverse of the pivot used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: dz_neglect ! A vertical distance that is so small it is usually lost
                    ! in roundoff and can be neglected [Z ~> m].
  real :: dMass     ! The mass per unit area within a layer [Z R ~> kg m-2].
  real :: dPres     ! The hydrostatic pressure change across a layer [R Z2 T-2 ~> Pa = J m-3].

  real :: dt_h      ! The timestep divided by the averages of the vertical distances around
                    ! a layer [T Z-1 ~> s m-1].
  real :: dz_top    ! The distance from the surface [Z ~> m].
  real :: dz_rsum   ! The distance from the seafloor [Z ~> m].
  real :: I_dzsum   ! The inverse of dz_sum [Z-1 ~> m-1].
  real :: I_BBLD    ! The inverse of the current value of BBLD [Z-1 ~> m-1].
  real :: dz_tt     ! The distance from the surface or up to the next interface
                    ! that did not exhibit turbulent mixing from this scheme plus
                    ! a surface mixing roughness length given by dz_tt_min [Z ~> m].
  real :: dz_tt_min ! A surface roughness length [Z ~> m].
  real :: C1_3      ! = 1/3  [nondim]
  real :: vstar     ! An in-situ turbulent velocity [Z T-1 ~> m s-1].
  real :: BBLD_output ! The bottom boundary layer depth output from this routine [Z ~> m]
  real :: hbs_here  ! The local minimum of dztop_dztot and MixLen_shape [nondim]
  real :: TKE_used  ! The TKE used to support mixing at an interface [R Z3 T-2 ~> J m-2].
  real :: Kd_guess0    ! A first guess of the diapycnal diffusivity [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real :: PE_chg_g0    ! The potential energy change when Kd is Kd_guess0 [R Z3 T-2 ~> J m-2]
  real :: Kddt_h_g0    ! The first guess of the change in diapycnal diffusivity times a timestep
                       ! divided by the average thicknesses around an interface [H ~> m or kg m-2].
  real :: Kddt_h_prev  ! The diapycnal diffusivity before modification times a timestep divided
                       ! by the average thicknesses around an interface [H ~> m or kg m-2].
  real :: dPEc_dKd_Kd0 ! The partial derivative of PE change with Kddt_h(K)
                       ! for very small values of Kddt_h(K) [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real :: PE_chg    ! The change in potential energy due to mixing at an
                    ! interface [R Z3 T-2 ~> J m-2], positive for the column increasing
                    ! in potential energy (i.e., consuming TKE).
  real :: TKE_left  ! The amount of turbulent kinetic energy left for the most
                    ! recent guess at Kddt_h(K) [R Z3 T-2 ~> J m-2].
  real :: dPEc_dKd  ! The partial derivative of PE_chg with Kddt_h(K) [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real :: TKE_left_min, TKE_left_max ! Maximum and minimum values of TKE_left [R Z3 T-2 ~> J m-2].
  real :: Kddt_h_max, Kddt_h_min ! Maximum and minimum values of Kddt_h(K) [H ~> m or kg m-2].
  real :: Kddt_h_guess ! A guess at the value of Kddt_h(K) [H ~> m or kg m-2].
  real :: Kddt_h_next  ! The next guess at the value of Kddt_h(K) [H ~> m or kg m-2].
  real :: dKddt_h      ! The change between guesses at Kddt_h(K) [H ~> m or kg m-2].
  real :: dKddt_h_Newt ! The change between guesses at Kddt_h(K) with Newton's method [H ~> m or kg m-2].
  real :: Kddt_h_newt  ! The Newton's method next guess for Kddt_h(K) [H ~> m or kg m-2].
  real :: exp_kh       ! The nondimensional decay of TKE across a layer [nondim].
  real :: frac_in_BL   ! The fraction of the energy required to support dKd_max that is suppiled by
                       ! max_PE_chg, used here to determine a fractional layer contribution to the
                       ! boundary layer thickness [nondim]
  real :: TKE_rescale  ! The effective fractional increase in energy available to
                       ! mixing at this interface once its exponential decay is accounted for [nondim]
  logical :: use_Newt  ! Use Newton's method for the next guess at Kddt_h(K).
  logical :: convectively_stable ! If true the water column is convectively stable at this interface.
  logical :: bot_connected   ! If true the ocean is actively turbulent from the present
                    ! interface all the way down to the seafloor.
  logical :: bot_disconnect ! If true, any turbulence has become disconnected
                    ! from the bottom.

  ! The following is only used for diagnostics.
  real :: I_dtdiag  !  = 1.0 / dt [T-1 ~> s-1].

  real :: BBLD_guess, BBLD_found ! Bottom boundary layer depth guessed/found for iteration [Z ~> m]
  real :: min_BBLD, max_BBLD ! Iteration bounds on BBLD [Z ~> m], which are adjusted at each step
  real :: dBBLD_min  ! The change in diagnosed mixed layer depth when the guess is min_BLD [Z ~> m]
  real :: dBBLD_max  ! The change in diagnosed mixed layer depth when the guess is max_BLD [Z ~> m]
  logical :: BBL_converged ! Flag for convergence of BBLD
  integer :: BBL_it        ! Iteration counter

  real :: Surface_Scale ! Surface decay scale for vstar [nondim]
  logical :: debug      ! This is used as a hard-coded value for debugging.
  logical :: no_MKE_conversion  ! If true, there is conversion of MKE to TKE in this routine.
  real :: mstar_BBL !< the value of mstar for the bottom boundary layer [nondim]

  !  The following arrays are used only for debugging purposes.
  real :: dPE_debug     ! An estimate of the potential energy change [R Z3 T-2 ~> J m-2]
  real :: mixing_debug  ! An estimate of the rate of change of potential energy due to mixing [R Z3 T-3 ~> W m-2]
  real, dimension(20) :: TKE_left_itt   ! The value of TKE_left after each iteration [R Z3 T-2 ~> J m-2]
  real, dimension(20) :: PE_chg_itt     ! The value of PE_chg after each iteration [R Z3 T-2 ~> J m-2]
  real, dimension(20) :: Kddt_h_itt     ! The value of Kddt_h_guess after each iteration [H ~> m or kg m-2]
  real, dimension(20) :: dPEa_dKd_itt   ! The value of dPEc_dKd after each iteration [R Z3 T-2 H-1 ~> J m-3 or J kg-1]
!  real, dimension(20) :: MKE_src_itt    ! The value of MKE_src after each iteration [R Z3 T-2 ~> J m-2]
  real, dimension(SZK_(GV)) :: dT_expect !< Expected temperature changes [C ~> degC]
  real, dimension(SZK_(GV)) :: dS_expect !< Expected salinity changes [S ~> ppt]
  real, dimension(SZK_(GV)) :: mech_BBL_TKE_k  ! The mechanically generated turbulent kinetic energy
                    ! available for bottom boundary mixing over a time step for each layer [R Z3 T-2 ~> J m-2].
  integer, dimension(SZK_(GV)) :: num_itts

  integer :: k, nz, itt, max_itt

  nz = GV%ke

  debug = .false.  ! Change this hard-coded value for debugging.
  no_MKE_conversion = ((CS%direct_calc) ) ! .and. (CS%MKE_to_TKE_effic == 0.0))

  ! Add bottom boundary layer mixing if there is energy to support it.
  if (((CS%ePBL_BBL_effic <= 0.0) .and. (CS%ePBL_tidal_effic <= 0.0) .and. (.not.CS%ePBL_BBL_use_mstar)) &
      .or. (BBL_TKE_in <= 0.0)) then
    ! There is no added bottom boundary layer mixing.
    BBLD_io = 0.0
    Kd_BBL(:) = 0.0
    mixvel_BBL(:) = 0.0 ; mixlen_BBL(:) = 0.0
    eCD%BBL_its = 0
    if (CS%TKE_diagnostics) then
      eCD%dTKE_BBL_mixing = 0.0 ; eCD%dTKE_BBL_decay = 0.0 ; eCD%dTKE_BBL = 0.0
      ! eCD%dTKE_BBL_MKE = 0.0
    endif
    return
  else
    ! There will be added bottom boundary layer mixing.

    h_neglect = GV%H_subroundoff
    dz_neglect = GV%dZ_subroundoff

    C1_3 = 1.0 / 3.0
    I_dtdiag = 1.0 / dt
    max_itt = 20
    dz_tt_min = 0.0

    ! The next two blocks of code could be shared with ePBL_column.

    ! Set up fields relating a layer's temperature and salinity changes to potential energy changes.
    pres_Z(1) = 0.0
    do k=1,nz
      dMass = GV%H_to_RZ * h(k)
      dPres = GV%g_Earth_Z_T2 * dMass
      dT_to_dPE(k) = (dMass * (pres_Z(K) + 0.5*dPres)) * dSV_dT(k)
      dS_to_dPE(k) = (dMass * (pres_Z(K) + 0.5*dPres)) * dSV_dS(k)
      dT_to_dColHt(k) = dMass * dSV_dT(k)
      dS_to_dColHt(k) = dMass * dSV_dS(k)

      pres_Z(K+1) = pres_Z(K) + dPres
    enddo

    if (GV%Boussinesq) then
      do K=1,nz+1 ; h_dz_int(K) = GV%Z_to_H ; enddo
    else
      h_dz_int(1) = (h(1) + h_neglect) / (dz(1) + dz_neglect)
      do K=2,nz
        h_dz_int(K) = (h(k-1) + h(k) + h_neglect) / (dz(k-1) + dz(k) + dz_neglect)
      enddo
      h_dz_int(nz+1) = (h(nz) + h_neglect) / (dz(nz) + dz_neglect)
    endif
    ! The two previous blocks of code could be shared with ePBL_column.

    ! Determine the total thickness (dz_sum) and the fractional distance from the top (dztop_dztot).
    dz_sum = 0.0 ; do k=1,nz ; dz_sum = dz_sum + dz(k) ; enddo
    I_dzsum = 0.0 ; if (dz_sum > 0.0) I_dzsum = 1.0 / dz_sum
    dz_top = 0.0
    dztop_dztot(nz+1) = 0.0
    do k=1,nz
      dz_top = dz_top + dz(k)
      dztop_dztot(K) = dz_top * I_dzsum
    enddo

    ! Set terms from a tridiagonal solver based on the previously determined diffusivities.
    Kddt_h(1) = 0.0
    hp_a(1) = h(1)
    dT_to_dPE_a(1) = dT_to_dPE(1) ; dT_to_dColHt_a(1) = dT_to_dColHt(1)
    dS_to_dPE_a(1) = dS_to_dPE(1) ; dS_to_dColHt_a(1) = dS_to_dColHt(1)
    do K=2,nz
      dt_h = dt / max(0.5*(dz(k-1)+dz(k)), 1e-15*dz_sum)
      Kddt_h(K) = Kd(K) * dt_h
      b1 = 1.0 / (hp_a(k-1) + Kddt_h(K))
      c1(K) = Kddt_h(K) * b1
      hp_a(k) = h(k) + (hp_a(k-1) * b1) * Kddt_h(K)
      dT_to_dPE_a(k) = dT_to_dPE(k) + c1(K)*dT_to_dPE_a(k-1)
      dS_to_dPE_a(k) = dS_to_dPE(k) + c1(K)*dS_to_dPE_a(k-1)
      dT_to_dColHt_a(k) = dT_to_dColHt(k) + c1(K)*dT_to_dColHt_a(k-1)
      dS_to_dColHt_a(k) = dS_to_dColHt(k) + c1(K)*dS_to_dColHt_a(k-1)
      if (K<=2) then
        Te(k-1) = b1*(h(k-1)*T0(k-1)) ; Se(k-1) = b1*(h(k-1)*S0(k-1))
        Th_a(k-1) = h(k-1) * T0(k-1) ; Sh_a(k-1) = h(k-1) * S0(k-1)
      else
        Te(k-1) = b1 * (h(k-1) * T0(k-1) + Kddt_h(K-1) * Te(k-2))
        Se(k-1) = b1 * (h(k-1) * S0(k-1) + Kddt_h(K-1) * Se(k-2))
        Th_a(k-1) = h(k-1) * T0(k-1) + Kddt_h(K-1) * Te(k-2)
        Sh_a(k-1) = h(k-1) * S0(k-1) + Kddt_h(K-1) * Se(k-2)
      endif
    enddo
    Kddt_h(nz+1) = 0.0
    if (debug) then
      ! Complete the tridiagonal solve for Te and Se, which may be useful for debugging.
      b1 = 1.0 / hp_a(nz)
      Te(nz) = b1 * (h(nz) * T0(nz) + Kddt_h(nz) * Te(nz-1))
      Se(nz) = b1 * (h(nz) * S0(nz) + Kddt_h(nz) * Se(nz-1))
      do k=nz-1,1,-1
        Te(k) = Te(k) + c1(K+1)*Te(k+1)
        Se(k) = Se(k) + c1(K+1)*Se(k+1)
      enddo
    endif

    BBLD_guess = BBLD_io

    !/The following lines are for the iteration over BBLD
    ! max_BBLD will initialized as ocean bottom depth
    max_BBLD = 0.0 ; do k=1,nz ; max_BBLD = max_BBLD + dz(k) ; enddo
    ! min_BBLD will be initialized to 0.
    min_BBLD = 0.0
    ! Set values of the wrong signs to indicate that these changes are not based on valid estimates
    dBBLD_min = -1.0*US%m_to_Z ; dBBLD_max = 1.0*US%m_to_Z

    ! If no first guess is provided for BBLD, try the middle of the water column
    if (BBLD_guess <= min_BBLD) BBLD_guess = 0.5 * (min_BBLD + max_BBLD)

    ! Iterate to determine a converged EPBL bottom boundary layer depth.
    do BBL_it=1,CS%max_BBLD_its

      if (debug) then ; mech_BBL_TKE_k(:) = 0.0 ; endif

      ! Reset BBL_depth
      BBLD_output = dz(nz)
      bot_connected = .true.

      if (CS%ePBL_BBL_use_mstar) then
        call find_mstar(CS, US, B_flux_BBL, u_star_BBL_z_t, BBLD_guess, absf, .true., mstar_BBL)
        eCD%mstar_BBL = mstar_BBL
        mech_BBL_TKE = mstar_BBL * BBL_TKE_in
      else
        mech_BBL_TKE = BBL_TKE_in
        eCD%mstar_BBL = 0.0
      endif
      if (CS%TKE_diagnostics) then
        ! eCD%dTKE_BBL_MKE = 0.0
        eCD%dTKE_BBL_mixing = 0.0
        eCD%dTKE_BBL_decay = 0.0
        eCD%dTKE_BBL = mech_BBL_TKE * I_dtdiag
      endif

      ! Store in 1D arrays for output.
      do K=1,nz+1 ; mixvel_BBL(K) = 0.0 ; mixlen_BBL(K) = 0.0 ; enddo

      ! Determine the mixing shape function MixLen_shape.
      if ((.not.CS%Use_BBLD_iteration) .or. &
          (CS%transLay_scale >= 1.0) .or. (CS%transLay_scale < 0.0) ) then
        do K=1,nz+1
          MixLen_shape(K) = 1.0
        enddo
      elseif (BBLD_guess <= 0.0) then
        if (CS%transLay_scale > 0.0) then ; do K=1,nz+1
          MixLen_shape(K) = CS%transLay_scale
        enddo ; else ; do K=1,nz+1
          MixLen_shape(K) = 1.0
        enddo ; endif
      else
        ! Reduce the mixing length based on BBLD, with a quadratic
        ! expression that follows KPP.
        I_BBLD = 1.0 / BBLD_guess
        dz_rsum = 0.0
        MixLen_shape(nz+1) = 1.0
        if (CS%MixLenExponent_BBL==2.0) then
          do K=nz,1,-1
            dz_rsum = dz_rsum + dz(k)
            MixLen_shape(K) = CS%transLay_scale + (1.0 - CS%transLay_scale) * &
                 (max(0.0, (BBLD_guess - dz_rsum)*I_BBLD) )**2
          enddo
        elseif (CS%MixLenExponent_BBL==1.0) then
          do K=nz,1,-1
            dz_rsum = dz_rsum + dz(k)
            MixLen_shape(K) = CS%transLay_scale + (1.0 - CS%transLay_scale) * &
                 (max(0.0, (BBLD_guess - dz_rsum)*I_BBLD) )
          enddo
        else ! (CS%MixLenExponent_BBL /= 2.0 or 1.0) then
          do K=nz,1,-1
            dz_rsum = dz_rsum + dz(k)
            MixLen_shape(K) = CS%transLay_scale + (1.0 - CS%transLay_scale) * &
                 (max(0.0, (BBLD_guess - dz_rsum)*I_BBLD) )**CS%MixLenExponent_BBL
          enddo
        endif
      endif

      Kd_BBL(nz+1) = 0.0 ; Kddt_h(nz+1) = 0.0
      hp_b(nz) = h(nz)
      dT_to_dPE_b(nz) = dT_to_dPE(nz) ; dT_to_dColHt_b(nz) = dT_to_dColHt(nz)
      dS_to_dPE_b(nz) = dS_to_dPE(nz) ; dS_to_dColHt_b(nz) = dS_to_dColHt(nz)

      htot = h(nz) ; dztot = dz(nz) ; uhtot = u(nz)*h(nz) ; vhtot = v(nz)*h(nz)

      if (debug) then
        mech_BBL_TKE_k(nz) = mech_BBL_TKE
        num_itts(:) = -1
      endif

      Idecay_len_TKE = (CS%TKE_decay_BBL * absf) / u_star_BBL
      do K=nz,2,-1
        !   Apply dissipation to the TKE, here applied as an exponential decay
        ! due to 3-d turbulent energy being lost to inefficient rotational modes.
        ! The following form is often used for mechanical stirring from the surface.
        ! There could be several different "flavors" of TKE that decay at different rates.
        exp_kh = 1.0
        if (Idecay_len_TKE > 0.0) exp_kh = exp(-h(k)*Idecay_len_TKE)
        if (CS%TKE_diagnostics) &
          eCD%dTKE_BBL_decay = eCD%dTKE_BBL_decay + (exp_kh-1.0) * mech_BBL_TKE * I_dtdiag
        mech_BBL_TKE = mech_BBL_TKE * exp_kh

        if (debug) then
          mech_BBL_TKE_k(K) = mech_BBL_TKE
        endif

        ! Precalculate some temporary expressions that are independent of Kddt_h(K).
        dt_h = dt / max(0.5*(dz(k-1)+dz(k)), 1e-15*dz_sum)

        !   This tests whether the layers above and below this interface are in
        ! a convectively stable configuration, without considering any effects of
        ! mixing at higher interfaces.  It is an approximation to the more
        ! complete test dPEc_dKd_Kd0 >= 0.0, that would include the effects of
        ! mixing across interface K+1.  The dT_to_dColHt here are effectively
        ! mass-weighted estimates of dSV_dT.
        Convectively_stable = ( 0.0 <= &
             ( (dT_to_dColHt(k) + dT_to_dColHt(k-1) ) * (T0(k-1)-T0(k)) + &
               (dS_to_dColHt(k) + dS_to_dColHt(k-1) ) * (S0(k-1)-S0(k)) ) )

        if ((mech_BBL_TKE <= 0.0) .and. Convectively_stable) then
          ! Energy is already exhausted, so set Kd_BBL = 0 and cycle or exit?
          mech_BBL_TKE = 0.0
          Kd_BBL(K) = 0.0 ; Kddt_h(K) = Kd(K) * dt_h
          bot_disconnect = .true.
          ! if (.not.debug) exit

        else ! mech_BBL_TKE > 0.0 or this is a potentially convectively unstable profile.
          bot_disconnect = .false.

          ! Precalculate some more temporary expressions that are independent of Kddt_h(K).
          if (K>=nz) then
            Th_b(k) = h(k) * T0(k) ; Sh_b(k) = h(k) * S0(k)
          else
            Th_b(k) = h(k) * T0(k) + Kddt_h(K+1) * Te(k+1)
            Sh_b(k) = h(k) * S0(k) + Kddt_h(K+1) * Se(k+1)
          endif

          !   Using Pr=1 and the diffusivity at the upper interface (once it is
          ! known), determine how much resolved mean kinetic energy (MKE) will be
          ! extracted within a timestep and add a fraction CS%MKE_to_TKE_effic of
          ! this to the mTKE budget available for mixing in the next layer.
          ! This is not enabled yet for BBL mixing.
         !  if ((CS%MKE_to_TKE_effic > 0.0) .and. (htot*h(k-1) > 0.0)) then
         !    ! This is the energy that would be available from homogenizing the
         !    ! velocities between layer k-1 and the layers below.
         !    dMKE_max = (GV%H_to_RZ * CS%MKE_to_TKE_effic) * 0.5 * &
         !        (h(k-1) / ((htot + h(k-1))*htot)) * &
         !        ((uhtot-u(k-1)*htot)**2 + (vhtot-v(k-1)*htot)**2)
         !    ! A fraction (1-exp(Kddt_h*MKE2_Hharm)) of this energy would be
         !    ! extracted by mixing with a finite viscosity.
         !    MKE2_Hharm = (htot + h(k-1) + 2.0*h_neglect) / &
         !                 ((htot+h_neglect) * (h(k-1)+h_neglect))
         !  else
         !    dMKE_max = 0.0
         !    MKE2_Hharm = 0.0
         !  endif

          ! At this point, Kddt_h(K) will be unknown because its value may depend
          ! on how much energy is available.
          dz_tt = dztot + dz_tt_min
          if (mech_BBL_TKE > 0.0) then
            if (CS%wT_scheme_BBL==wT_from_cRoot_TKE) then
              vstar = CS%vstar_scale_fac_BBL * cuberoot(SpV_dt(K)*mech_BBL_TKE)
            elseif (CS%wT_scheme_BBL==wT_from_RH18) then
              Surface_Scale = max(0.05, 1.0 - dztot / BBLD_guess)
              vstar = (CS%vstar_scale_fac_BBL * Surface_Scale) * ( CS%vstar_surf_fac_BBL*u_star_BBL/h_dz_int(K) )
            endif
            hbs_here = min(dztop_dztot(K), MixLen_shape(K))
            mixlen_BBL(K) = max(CS%min_BBL_mix_len, ((dz_tt*hbs_here)*vstar) / &
                ((CS%Ekman_scale_coef_BBL * absf) * (dz_tt*hbs_here) + vstar))
            Kd_guess0 = (h_dz_int(K)*vstar) * CS%vonKar * mixlen_BBL(K)
          else
            vstar = 0.0 ; Kd_guess0 = 0.0
          endif
          mixvel_BBL(K) = vstar ! Track vstar

          TKE_rescale = 1.0
          if (CS%decay_adjusted_BBL_TKE) then
            ! Add a scaling factor that accounts for the exponential decay of TKE from a
            ! near-bottom source and the assumption that an increase in the diffusivity at an
            ! interface causes a linearly increasing buoyancy flux going from 0 at the bottom
            ! to a peak at the interface, and then going back to 0 atop the layer above.
            TKE_rescale = exp_decay_TKE_adjust(htot, h(k-1), Idecay_len_TKE)
          endif

          TKE_eff_avail = TKE_rescale*mech_BBL_TKE

          if (no_MKE_conversion) then
            ! Without conversion from MKE to TKE, the updated diffusivity can be determined directly.
            call find_Kd_from_PE_chg(Kd(K), Kd_guess0, dt_h, TKE_eff_avail, hp_a(k-1), hp_b(k), &
                        Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                        dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                        pres_Z(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                        dT_to_dColHt_b(k), dS_to_dColHt_b(k), Kd_add=Kd_BBL(K), PE_chg=TKE_eff_used, &
                        frac_dKd_max_PE=frac_in_BL)

            ! Do not add energy if the column is convectively unstable.  This was handled previously
            ! for mixing from the surface.
            if (TKE_eff_used < 0.0) TKE_eff_used = 0.0

            ! Convert back to the TKE that has actually been used.
            if (CS%decay_adjusted_BBL_TKE) then
              if (TKE_rescale == 0.0) then  ! This probably never occurs, even at roundoff.
                TKE_used = mech_BBL_TKE  ! All the energy was dissipated before it could mix.
              else
                TKE_used = TKE_eff_used / TKE_rescale
              endif
            else
              TKE_used = TKE_eff_used
            endif

            if (bot_connected)  BBLD_output = BBLD_output + frac_in_BL*dz(k-1)
            if (frac_in_BL < 1.0) bot_disconnect = .true.

            if (CS%TKE_diagnostics) then
              eCD%dTKE_BBL_mixing = eCD%dTKE_BBL_mixing - TKE_eff_used * I_dtdiag
              eCD%dTKE_BBL_decay = eCD%dTKE_BBL_decay - (TKE_used-TKE_eff_used) * I_dtdiag
            endif

            mech_BBL_TKE = mech_BBL_TKE - TKE_used

            Kddt_h(K) = (Kd(K) + Kd_BBL(K)) * dt_h

          else
            Kddt_h_prev = Kd(K) * dt_h
            Kddt_h_g0 = Kd_guess0 * dt_h
            ! Find the change in PE with the guess at the added bottom boundary layer mixing.
            call find_PE_chg(Kddt_h_prev, Kddt_h_g0, hp_a(k-1), hp_b(k), &
                       Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                       dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                       pres_Z(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                       dT_to_dColHt_b(k), dS_to_dColHt_b(k), &
                       PE_chg=PE_chg_g0, dPEc_dKd_0=dPEc_dKd_Kd0 )

            ! MKE_src = 0.0 ! Enable later?: = dMKE_max*(1.0 - exp(-Kddt_h_g0 * MKE2_Hharm))

            ! Do not add energy if the column is convectively unstable.  This was handled previously
            ! for mixing from the surface.
            if (PE_chg_g0 < 0.0) PE_chg_g0 = 0.0

            ! This block checks out different cases to determine Kd at the present interface.
            ! if (mech_BBL_TKE*TKE_rescale + (MKE_src - PE_chg_g0) >= 0.0) then
            if (TKE_eff_avail  - PE_chg_g0 >= 0.0) then
              ! This column is convectively stable and there is energy to support the suggested
              ! mixing, or it is convectively unstable.  Keep this first estimate of Kd.
              Kd_BBL(K) = Kd_guess0
              Kddt_h(K) = Kddt_h_prev + Kddt_h_g0

              TKE_used = PE_chg_g0 / TKE_rescale

              if (CS%TKE_diagnostics) then
                eCD%dTKE_BBL_mixing = eCD%dTKE_BBL_mixing - PE_chg_g0 * I_dtdiag
!                eCD%dTKE_MKE = eCD%dTKE_MKE + MKE_src * I_dtdiag
                eCD%dTKE_BBL_decay = eCD%dTKE_BBL_decay - (TKE_used - PE_chg_g0) * I_dtdiag
              endif

              ! mech_BBL_TKE = mech_BBL_TKE + MKE_src - TKE_used
              mech_BBL_TKE = mech_BBL_TKE - TKE_used
              if (bot_connected) then
                BBLD_output = BBLD_output + dz(k-1)
              endif

            elseif (TKE_eff_avail == 0.0) then
              ! This can arise if there is no energy input to drive mixing or if there
              ! is such strong decay that the mech_BBL_TKE becomes 0 via an underflow.
              Kd_BBL(K) = 0.0 ; Kddt_h(K) = Kddt_h_prev
              if (CS%TKE_diagnostics) then
                eCD%dTKE_BBL_decay = eCD%dTKE_BBL_decay - mech_BBL_TKE * I_dtdiag
              endif
              mech_BBL_TKE = 0.0
              bot_disconnect = .true.
            else
              ! There is not enough energy to support the mixing, so reduce the
              ! diffusivity to what can be supported.
              Kddt_h_max = Kddt_h_g0 ; Kddt_h_min = 0.0
              ! TKE_left_max = TKE_eff_avail + (MKE_src - PE_chg_g0)
              TKE_left_max = TKE_eff_avail - PE_chg_g0
              TKE_left_min = TKE_eff_avail

              ! As a starting guess, take the minimum of a false position estimate
              ! and a Newton's method estimate starting from dKddt_h = 0.0
              ! Enable conversion from MKE to TKE in the bottom boundary layer later?
              ! Kddt_h_guess = TKE_eff_avail * Kddt_h_max / max( PE_chg_g0  - MKE_src, &
              !                  Kddt_h_max * (dPEc_dKd_Kd0 - dMKE_max * MKE2_Hharm) )
              Kddt_h_guess = TKE_eff_avail * Kddt_h_max / max( PE_chg_g0, Kddt_h_max * dPEc_dKd_Kd0 )
              ! The above expression is mathematically the same as the following
              ! except it is not susceptible to division by zero when
              !   dPEc_dKd_Kd0 = dMKE_max = 0 .
              !  Kddt_h_guess = TKE_eff_avail * min( Kddt_h_max / (PE_chg_g0 - MKE_src), &
              !                      1.0 / (dPEc_dKd_Kd0 - dMKE_max * MKE2_Hharm) )
              if (debug) then
                TKE_left_itt(:) = 0.0 ; dPEa_dKd_itt(:) = 0.0 ; PE_chg_itt(:) = 0.0
                Kddt_h_itt(:) = 0.0 ! ; MKE_src_itt(:) = 0.0
              endif
              do itt=1,max_itt
                call find_PE_chg(Kddt_h_prev, Kddt_h_guess, hp_a(k-1), hp_b(k), &
                           Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                           dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE_b(k), dS_to_dPE_b(k), &
                           pres_Z(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                           dT_to_dColHt_b(k), dS_to_dColHt_b(k), &
                           PE_chg=PE_chg, dPEc_dKd=dPEc_dKd)
                ! Enable conversion from MKE to TKE in the bottom boundary layer later?
                ! MKE_src = dMKE_max * (1.0 - exp(-MKE2_Hharm * Kddt_h_guess))
                ! dMKE_src_dK = dMKE_max * MKE2_Hharm * exp(-MKE2_Hharm * Kddt_h_guess)

                ! TKE_left = TKE_eff_avail + (MKE_src - PE_chg)
                TKE_left = TKE_eff_avail - PE_chg
                if (debug .and. itt<=20) then
                  Kddt_h_itt(itt) = Kddt_h_guess ! ; MKE_src_itt(itt) = MKE_src
                  PE_chg_itt(itt) = PE_chg ; dPEa_dKd_itt(itt) = dPEc_dKd
                  TKE_left_itt(itt) = TKE_left
                endif
                ! Store the new bounding values, bearing in mind that min and max
                ! here refer to Kddt_h and dTKE_left/dKddt_h < 0:
                if (TKE_left >= 0.0) then
                  Kddt_h_min = Kddt_h_guess ; TKE_left_min = TKE_left
                else
                  Kddt_h_max = Kddt_h_guess ; TKE_left_max = TKE_left
                endif

                ! Try to use Newton's method, but if it would go outside the bracketed
                ! values use the false-position method instead.
                use_Newt = .true.
                ! if (dPEc_dKd - dMKE_src_dK <= 0.0) then
                if (dPEc_dKd <= 0.0) then
                  use_Newt = .false.
                else
                  ! dKddt_h_Newt = TKE_left / (dPEc_dKd - dMKE_src_dK)
                  dKddt_h_Newt = TKE_left / dPEc_dKd
                  Kddt_h_Newt = Kddt_h_guess + dKddt_h_Newt
                  if ((Kddt_h_Newt > Kddt_h_max) .or. (Kddt_h_Newt < Kddt_h_min)) &
                    use_Newt = .false.
                endif

                if (use_Newt) then
                  Kddt_h_next = Kddt_h_guess + dKddt_h_Newt
                  dKddt_h = dKddt_h_Newt
                else
                  Kddt_h_next = (TKE_left_max * Kddt_h_min - Kddt_h_max * TKE_left_min) / &
                                (TKE_left_max - TKE_left_min)
                  dKddt_h = Kddt_h_next - Kddt_h_guess
                endif

                if ((abs(dKddt_h) < 1e-9*Kddt_h_guess) .or. (itt==max_itt)) then
                  ! Use the old value so that the energy calculation does not need to be repeated.
                  if (debug) num_itts(K) = itt
                  exit
                else
                  Kddt_h_guess = Kddt_h_next
                endif
              enddo ! Inner iteration loop on itt.
              Kd_BBL(K) = Kddt_h_guess / dt_h
              Kddt_h(K) = (Kd(K) + Kd_BBL(K)) * dt_h

              ! All TKE should have been consumed.
              if (CS%TKE_diagnostics) then
                ! eCD%dTKE_BBL_mixing = eCD%dTKE_BBL_mixing - (TKE_eff_avail + MKE_src) * I_dtdiag
                ! eCD%dTKE_BBL_MKE = eCD%dTKE_BBL_MKE + MKE_src * I_dtdiag
                eCD%dTKE_BBL_mixing = eCD%dTKE_BBL_mixing - TKE_eff_avail * I_dtdiag
                eCD%dTKE_BBL_decay = eCD%dTKE_BBL_decay - (mech_BBL_TKE-TKE_eff_avail) * I_dtdiag
              endif

              if (bot_connected) BBLD_output = BBLD_output + (PE_chg / PE_chg_g0) * dz(k-1)

              mech_BBL_TKE = 0.0
              bot_disconnect = .true.
            endif ! End of convective or forced mixing cases to determine Kd.
          endif

          Kddt_h(K) = (Kd(K) + Kd_BBL(K)) * dt_h
        endif  ! tot_TKT > 0.0 branch.  Kddt_h(K) has been set.

        !   At this point, the final value of Kddt_h(K) is known, so the
        ! estimated properties for layer k can be calculated.
        b1 = 1.0 / (hp_b(k) + Kddt_h(K))
        c1(K) = Kddt_h(K) * b1

        hp_b(k-1) = h(k-1) + (hp_b(k) * b1) * Kddt_h(K)
        dT_to_dPE_b(k-1) = dT_to_dPE(k-1) + c1(K)*dT_to_dPE_b(k)
        dS_to_dPE_b(k-1) = dS_to_dPE(k-1) + c1(K)*dS_to_dPE_b(k)
        dT_to_dColHt_b(k-1) = dT_to_dColHt(k-1) + c1(K)*dT_to_dColHt_b(k)
        dS_to_dColHt_b(k-1) = dS_to_dColHt(k-1) + c1(K)*dS_to_dColHt_b(k)

        ! Store integrated velocities and thicknesses for MKE conversion calculations.
        if (bot_disconnect) then
          ! There is no turbulence at this interface, so restart the running sums.
          uhtot = u(k-1)*h(k-1)
          vhtot = v(k-1)*h(k-1)
          htot  = h(k-1)
          dztot = dz(k-1)
          bot_connected = .false.
        else
          uhtot = uhtot + u(k-1)*h(k-1)
          vhtot = vhtot + v(k-1)*h(k-1)
          htot  = htot + h(k-1)
          dztot = dztot + dz(k-1)
        endif

        if (K==nz) then
          Te(k) = b1*(h(k)*T0(k))
          Se(k) = b1*(h(k)*S0(k))
        else
          Te(k) = b1 * (h(k) * T0(k) + Kddt_h(K+1) * Te(k+1))
          Se(k) = b1 * (h(k) * S0(k) + Kddt_h(K+1) * Se(k+1))
        endif
      enddo
      Kd_BBL(1) = 0.0

      if (debug) then
        ! Complete the tridiagonal solve for Te with a downward pass.
        b1 = 1.0 / hp_b(1)
        Te(1) = b1 * (h(1) * T0(1) + Kddt_h(2) * Te(2))
        Se(1) = b1 * (h(1) * S0(1) + Kddt_h(2) * Se(2))
        dT_expect(1) = Te(1) - T0(1) ; dS_expect(1) = Se(1) - S0(1)
        do k=2,nz
          Te(k) = Te(k) + c1(K)*Te(k-1)
          Se(k) = Se(k) + c1(K)*Se(k-1)
          dT_expect(k) = Te(k) - T0(k) ; dS_expect(k) = Se(k) - S0(k)
        enddo

        dPE_debug = 0.0
        do k=1,nz
          dPE_debug = dPE_debug + (dT_to_dPE(k) * (Te(k) - T0(k)) + &
                                   dS_to_dPE(k) * (Se(k) - S0(k)))
        enddo
        mixing_debug = dPE_debug * I_dtdiag
      endif

      ! Skip the rest of the contents of the do loop if there are no more BBL depth iterations.
      if (BBL_it >= CS%max_BBLD_its) exit

      ! The following lines are used for the iteration to determine the boundary layer depth.
      ! Note that the iteration uses the value predicted by the TKE threshold (BBL_DEPTH),
      ! because the mixing length shape is dependent on the BBL depth, and therefore the BBL depth
      ! should be estimated more precisely than the grid spacing.

      ! New method uses BBL_DEPTH as computed in ePBL routine
      BBLD_found = BBLD_output
      if (abs(BBLD_found - BBLD_guess) < CS%BBLD_tol) then
        exit ! Break the BBL depth convergence loop
      elseif (BBLD_found > BBLD_guess) then
        min_BBLD = BBLD_guess ; dBBLD_min = BBLD_found - BBLD_guess
      else ! We know this guess was too deep
        max_BBLD = BBLD_guess ; dBBLD_max = BBLD_found - BBLD_guess ! <= -CS%BBLD_tol
      endif

      ! Try using the false position method or the returned value instead of simple bisection.
      ! Taking the occasional step with BBLD_output empirically helps to converge faster.
      if ((dBBLD_min > 0.0) .and. (dBBLD_max < 0.0) .and. (BBL_it > 2) .and. (mod(BBL_it-1,4) > 0)) then
        ! Both bounds have valid change estimates and are probably in the range of possible outputs.
        BBLD_guess = (dBBLD_min*max_BBLD - dBBLD_max*min_BBLD) / (dBBLD_min - dBBLD_max)
      elseif ((BBLD_found > min_BBLD) .and. (BBLD_found < max_BBLD)) then
        ! The output BBLD_found is an interesting guess, as it is likely to bracket the true solution
        ! along with the previous value of BBLD_guess and to be close to the solution.
        BBLD_guess = BBLD_found
      else ! Bisect if the other guesses would be out-of-bounds.  This does not happen much.
        BBLD_guess = 0.5*(min_BBLD + max_BBLD)
      endif

    enddo ! Iteration loop for converged boundary layer thickness.

    eCD%BBL_its = min(BBL_it, CS%max_BBLD_its)
    BBLD_io = BBLD_output
  endif

end subroutine ePBL_BBL_column

!> Gives shape function that sets the vertical structure of OSBL diffusivity
!! as described in Sane et al. 2025
subroutine kappa_eqdisc(shape_func, CS, GV, dz, absf, B_flux, u_star, MLD_guess)

  type(verticalGrid_type), intent(in) :: GV     !< The ocean's vertical grid structure.
  type(energetic_PBL_CS),  intent(in) :: CS     !< Energetic PBL control struct
  real, dimension(SZK_(GV)+1), intent(inout) :: shape_func  !< shape function, [nondim]
  real, intent(in) :: absf      !< The absolute value of f [T-1 ~> s-1]
  real, intent(in) :: u_star    !< The surface friction velocity [Z T-1 ~> m s-1]
  real, intent(in) :: B_Flux    !< The surface buoyancy flux [Z2 T-3 ~> m2 s-3]
  real, dimension(SZK_(GV)), intent(in)  :: dz     !< The vertical distance across layers [Z ~> m]
  real, intent(in) :: MLD_guess !< Mixing Layer depth guessed/found for iteration [Z ~> m].
  real, dimension(SZK_(GV)+1) :: hz !< depth variable, only used in this routine [H ~> m]

  ! local variables for this subroutine
  integer :: nz
  integer :: K, n ! integers for looping
  real :: Lh ! ((B_flux * h))/(u_star^3), boundary layer depth by M-O depth, [nondim]
  real :: Eh ! ((h f)/u_star ),  boundary layer depth by Ekman depth, [nondim]
  real :: sm ! sigma_max: location of maximum of shape function in sigma coordinate [nondim]
  real :: hbl ! Boundary layer depth, same as MLD_guess [Z ~> m]
  real :: F ! function, used in asymptotic model for sm, Equation 7 in Sane et al. 2024 [nondim]
  real :: F_Eh ! F multiplied by Eh [nondim]
  real :: u_star_I  ! inverse of u_star [Z-1 T ~> m-1 s]

  ! variables used for optimizing computations:
  real :: sm_h     ! sigma_max multiplied by boundary layer depth [Z ~> m]
  real :: sm_h_I   ! inverse of sm_h,[Z-1 ~> m-1]
  real :: sm_h_I2  ! An inverse variable given by 1.0/(h - sm_h), [Z-1 ~> m-1]
  real :: hz_n     ! z depth to avoid calling hz multiple times [Z ~> m]
  real :: z_minus_sm_h  ! depth z minus \sigma_m * MLD_Guess [Z ~> m]
  real :: z_minus_sm_h2 ! (depth z minus \sigma_m * MLD_Guess)^2 [Z2 ~> m2]
  real :: z_minus_sm_h3 ! (depth z minus \sigma_m * MLD_Guess)^3 [Z3 ~> m3]
  real :: h_minus_smh_I ! inverse of (MLD_Guess - \sigma_m * MLD_Guess)  [Z-1 ~> m-1]
  real :: h_minus_smh_I2 ! inverse of (MLD_Guess - \sigma_m * MLD_Guess) ^ 2 [Z-2 ~> m-2]
  real :: h_minus_smh_I3 ! inverse of (MLD_Guess - \sigma_m * MLD_Guess) ^ 3 [Z-3 ~> m-3]
  real :: z_sm_h_I      ! depth divided by (\sigma_m * MLD_Guess) [nondim]
  real :: coef_c2       ! = 2.98 * h_minus_smh_I2 !  [Z-2 ~> m-2]
  real :: coef_c3       ! = 2.98 * h_minus_smh_I2 !  [Z-3 ~> m-3]

  nz = SZK_(GV)+1
  hz(1) = 0.0
  do K=2,nz
    hz(K) = hz(K-1) + dz(K-1)
  end do
  hbl = MLD_Guess ! hbl is boundary layer depth.

  u_star_I = 1.0/u_star
  Lh = (-B_flux * hbl) * ((u_star_I * u_star_I) * u_star_I) ! Boundary layer depth divided by Monin-Obukhov depth
  Eh = (hbl * absf) * u_star_I   ! Boundary layer depth divided by Ekman depth

  ! B_flux given negative sign to follow convention used in Sane et al. 2023
  ! Lh < 0 --> surface stabilizing i.e. heating, and Lh > 0 --> surface destabilizing i.e. cooling
  ! This capping does not matter because these equations have asymptotes. Not sensitive beyond the caps.
  Eh = min(Eh, CS%Eh_upper_cap) ! capping p1 to less than 2.0. It is always >0.0.
  Lh = min(max(Lh, -CS%Lh_cap), CS%Lh_cap) ! capping Lh between -8 and 8

  ! Empirical model to predict sm:
  ! F is Equation (6) in Sane et al. 2025, and needs to be computed before sigma_m:
  ! \mathcal{F} = \frac{1}{c_3 + c_4 \cdot e^{-\left( \text{sgn}(B) \cdot {c_5} \cdot {{L_h}^3} \right)}} + c_6
  ! Equation (5) in Sane et al. 2025:
  ! \sigma_{m} = \frac{1}{c_1 + \frac{c_2}{\mathcal{F} \cdot E_h}}
  ! Note: Lh over here is ((Bh)/ustar^3), whereas in Sane et al. 2025, L_h = (((Bh)^{1/3})/(ustar))

  F = (1.0/ ( CS%ML_c(3) + CS%ML_c(4) * exp(-CS%ML_c(5) * Lh) ) ) + CS%ML_c(6)
  F_Eh = F * Eh
  sm = F_Eh / (CS%ML_c(1)*F_Eh +CS%ML_c(2))
  sm = min(max(sm, CS%sigma_max_lower_cap), CS%sigma_max_upper_cap) ! makes sure 0.1<sm<0.7
                                                                    ! true sm range is (approx) 0.2 to 0.60

  sm_h = sm * hbl
  sm_h_I = 1.0/sm_h                                 ! 1.0 /  (sm x hbl)
  h_minus_smh_I  = 1.0/(hbl-sm_h)                   ! 1.0 /  (hbl-sm_h)
  h_minus_smh_I2 = h_minus_smh_I * h_minus_smh_I    !  (1.0 / (hbl - sm*hbl))^2
  h_minus_smh_I3 = h_minus_smh_I2 * h_minus_smh_I   !  (1.0 / (hbl - sm*hbl))^3

  ! The coefficients coef_c3 and coef_c2 are dependent on CS%shape_function_epsilon.
  ! Above depth sm_h, shape_func is quadratic, and below sm_h, it is cubic.
  ! For iterative ePBL solver, shape_func should not be zero below hbl, so that it has been set to a small value
  ! set by CS%shape_function_epsilon. To make the cubic part of shapefunc behave smoothly, the below two coefficients
  ! are used that depend on CS%shape_function_epsilon. The numbers 1.0, 2.0, 3.0 below are constants,
  ! and should not be changed.

  coef_c3 = ( 2.0 * ( 1.0 - CS%shape_function_epsilon ) ) * h_minus_smh_I3
  coef_c2 = ( 3.0 * ( CS%shape_function_epsilon - 1.0 ) ) * h_minus_smh_I2

  ! gives the shape, quadratic above sm, cubic below sm in sigma coordinate
  ! see Equation 3 in Sane et al. 2024
  ! interpolates a quadratic function from z=0 to z=sm_h, and then a cubic from z=sm_h to z=hbl

  shape_func(1) = 0.0  ! initializing the first element of shape function array
  do n = 2,nz
    hz_n = hz(n) ! calls hz(n) once to avoid calling it multiple times below

    if  (hz_n <= sm_h) then
      ! Eq.3a in Sane et al. 2025: -(\frac{z}{\sigma_m \cdot h})^2+2(\frac{z}{\sigma_m h}) : Eq. (3) in Sane et al. 2025

      z_sm_h_I = hz_n * sm_h_I ! pre multiplying
      shape_func(n) = -z_sm_h_I*z_sm_h_I + 2.0 * z_sm_h_I

    elseif  (hz_n <= hbl) then
      ! Eq.3b in Sane et al. 2025: 2\left(\frac{\s - \sm}{1 - \sm} \right)^3 -
      ! 3\left(\frac{\s - \sm}{1 - \sm} \right)^2 + 1

      z_minus_sm_h  = (hz_n - sm_h)
      z_minus_sm_h2 = z_minus_sm_h * z_minus_sm_h
      z_minus_sm_h3 = z_minus_sm_h * z_minus_sm_h2

      shape_func(n) = (coef_c3 * z_minus_sm_h3 + coef_c2 * z_minus_sm_h2) + 1.0

    elseif (hz(n) > hbl) then
      shape_func(n) = CS%shape_function_epsilon ! set an arbitrary low constant value below hbl, default 0.01
    endif
  end do
end subroutine kappa_eqdisc

!> Gives velocity scale (v_0) using equations that approximate neural network of Sane et al. 2023
subroutine get_eqdisc_v0(CS, absf, B_flux, u_star, v0_dummy)
  type(energetic_PBL_CS),  intent(in) :: CS     !< Energetic PBL control struct
  real, intent(in) :: B_flux !< The surface buoyancy flux [Z2 T-3 ~> m2 s-3]
  real, intent(in) :: u_star !< The surface friction velocity [Z T-1 ~> m s-1]
  real, intent(in) :: absf  !< The absolute value of f [T-1 ~> s-1].
  real, intent(inout) :: v0_dummy   !< velocity scale v0, local variable [Z T-1 ~> m s-1]

  ! local variables for this subroutine
  real :: bflux_c  ! capped bflux [Z2 T-3 ~> m2 s-3]
  real :: absf_c   ! capped absf [T-1 ~> s-1]
  real :: root_b_f ! square root of (abs(B_flux) * Coriolis) [Z T-2 ~> m s-2]
  real :: f_u2     ! Coriolis X ustar^2 [Z2 T-3 ~> m2 s-3]
  real :: den      ! denominator, units iof buuyancy flux [Z2 T-3 ~> m2 s-3]
  real :: root_B_by_Omega ! sqrt( B / Omega )   [Z T-1 ~> m s-1]
  real :: f_prime  ! Coriolis divided by Earth's rotation [nondim]
  real :: omega_I  ! Inverse of the Earth's rotation rate, 1 divided by omega [T ~> s]

  if (B_flux <= CS%bflux_lower_cap) then
    bflux_c = CS%bflux_lower_cap
  elseif (B_flux >= CS%bflux_upper_cap) then
    bflux_c = CS%bflux_upper_cap
  else
    bflux_c = B_flux
  endif

  if (absf <= CS%f_lower) then   !
    absf_c = CS%f_lower    ! 0.1 deg Latitude, cap avoids zero coriolis, solution insensitive below 0.1 deg.
  else
    absf_c = absf
  endif

  f_u2 = absf_c * (u_star * u_star) ! pre-computing

  ! setting v0_dummy here:
  ! \lambda = (1/ustar) \sqrt(bflux_c/absf_c)

  if (bflux_c >= 0.0) then ! surface heating and neutral conditions
  ! Equation 7 in Sane et al. 2025:
  ! \frac{v_0}{u_*} = \frac{c_{7}}{\lambda + c_{8} + \frac{c_{9}^2}{\lambda + c_{9}} }

    root_b_f = sqrt( bflux_c  * absf_c)
    den = bflux_c + (CS%ML_c(8) + CS%ML_c(9)) * u_star * root_b_f  + &
          (CS%ML_c(8) * CS%ML_c(9) + CS%ML_c(9)**2) * f_u2
    v0_dummy = ( ( CS%ML_c(7)*( (u_star * root_b_f) + (CS%ML_c(9)*f_u2) ) ) * u_star) / den

  else ! surface cooling
  ! Equation 8 in Sane et al. 2025:
  ! \frac{v_0}{u_*}=\frac{c_{10} \cdot \lambda \cdot \sqrt{f'} }{1 +
  ! \frac{(c_{11} e^{(-c_{12} \cdot f')} + c_{13}) }{\lambda ^2} } + c_{14}

    omega_I = 1.0 / CS%omega
    f_prime = absf_c * omega_I  ! Coriolis divided by Earth's rotation
    root_B_by_Omega = sqrt( -bflux_c * omega_I  )
    den = ( -bflux_c + CS%ML_c(11) * f_u2 * exp(-f_prime * CS%ML_c(12) ) ) + CS%ML_c(13)*f_u2
    v0_dummy = ( CS%ML_c(10) * (-bflux_c * root_B_by_Omega) / den  ) + ( CS%ML_c(14) * u_star )

  endif

  v0_dummy = min( max(v0_dummy, CS%v0_lower_cap), CS%v0_upper_cap )
  ! upper cap kept for safety, but has never hit this cap.

  ! v0_lower_cap has been set to 0.0001 as data below that values does not exist in the training
  ! solution was tested for lower cap of 0.00001 and was found to be insensitive.
  ! sensitivity arises when lower cap is 0.0. That is when diffusivity attains extremely low values and
  ! they go near molecular diffusivity. Boundary layers might become "sub-grid" i.e. < 1 metre
  ! some cause issues such as anomlous surface warming.
  ! this needs further investigation, our choices are motivated by practicallity for now.
end subroutine get_eqdisc_v0

!> Gives velocity scale (v_0^h) using equations that with using boundary layer depth as one of its inputs
!! These equations are different than those set in get_eqdisc_v0 subroutine
subroutine get_eqdisc_v0h(CS, B_flux, u_star, MLD_guess, v0_dummy)
  type(energetic_PBL_CS),  intent(in) :: CS     !< Energetic PBL control struct
  real, intent(in) :: B_flux !< The surface buoyancy flux [Z2 T-3 ~> m2 s-3]
  real, intent(in) :: u_star !< The surface friction velocity [Z T-1 ~> m s-1]
  real, intent(in) :: MLD_guess !< boundary layer depth guessed/found for iteration [Z ~> m]

  real, intent(inout) :: v0_dummy   !< velocity scale v0, local variable [Z T-1 ~> m s-1]

  ! local variables for this subroutine
  real :: bflux_c  ! capped bflux [Z2 T-3 ~> m2 s-3]
  real :: B_h, den ! Surface buoyancy flux multiplied by boundary layer depth, den is a denominator [Z3 T-3 ~> m3 s-3]
  real :: B_h_power1by3 ! cuberoot of (Surface buoyancy flux multiplied by boundary layer depth) [Z T-1 ~> m s-1]
  real :: u_star_2      ! u_star squared, [Z2 T-2 ~> m2 s-2]
  real :: u_star_3      ! u_star cubed,   [Z3 T-3 ~> m3 s-3]

  u_star_2 = u_star * u_star ! pre-multiplying to get ustar ^ 2
  u_star_3 = u_star_2 * u_star ! ustar ^ 3.0

  if (B_flux <= CS%bflux_lower_cap) then
    bflux_c = CS%bflux_lower_cap
  elseif (B_flux >= CS%bflux_upper_cap) then
    bflux_c = CS%bflux_upper_cap
  else
    bflux_c = B_flux
  endif

  B_h = abs(bflux_c) * MLD_guess
  B_h_power1by3 = cuberoot(B_h)

  ! setting v0_dummy here:

  if (bflux_c >= 0.0) then ! surface heating and neutral conditions
    ! Equation 9 in Sane et al. 2025:
    ! \frac{v_0^h}{u_*} = \frac{C_{14}}{ c_{15} L_h^3 + c_{16} L_h^2  + 1 }

    den = ( CS%ML_c(15) * B_h + CS%ML_c(16)* u_star*(B_h_power1by3*B_h_power1by3)) &
           + (u_star*u_star_2)
    v0_dummy = ( CS%ML_c(14) * (u_star_2 * u_star_2)) / den

  else
    ! Equation 10 in Sane et al. 2025:
    ! \frac{v_0^h}{u_*} = \frac{L_h}{c_{17} + \frac{c_{18}}{L_h ^2}}  + c_{14}
    den = CS%ML_c(17) * (B_h_power1by3*B_h_power1by3) + CS%ML_c(18) * u_star_2
    v0_dummy = (B_h / den ) + CS%ML_c(14) * u_star
  endif

  v0_dummy = min( max(v0_dummy, CS%v0_lower_cap), CS%v0_upper_cap )
  ! upper cap kept for safety, but has never hit this cap.

  ! v0_lower_cap has been set to 0.0001 as data below that values does not exist in the training
  ! solution was tested for lower cap of 0.00001 and was found to be insensitive.
  ! sensitivity arises when lower cap is 0.0. That is when diffusivity attains extremely low values and
  ! they go near molecular diffusivity. Boundary layers might become "sub-grid" i.e. < 1 metre
  ! some cause issues such as anomlous surface warming.
  ! this needs further investigation, our choices are motivated by practicallity for now.
end subroutine get_eqdisc_v0h

!> Determine a scaling factor that accounts for the exponential decay of turbulent kinetic energy
!! from a boundary source and the assumption that an increase in the diffusivity at an interface
!! causes a linearly increasing buoyancy flux going from 0 at the bottom to a peak at the interface,
!! and then going back to 0 atop the layer above.  Where this factor increases the available mixing
!! TKE, it is only compensating for the fact  that the TKE has already been reduced by the same
!! exponential decay rate.  ha and hb must be non-negative, and this function generally increases
!! with hb and decreases with ha.
!!
!!   Exp_decay_TKE_adjust is coded to have a lower bound of 1e-30 on the return value. For large
!! values of ha*Idecay, the return value is about 0.5*ka*(ha+hb)*Idecay**2 * exp(-ha*Idecay), but
!! return values of less than 1e-30 are deliberately reset to 1e-30.  For relatively large values
!! of hb*Idecay, the return value increases linearly with hb.  When Idecay ~= 0, the return value
!! is close to 1.
function exp_decay_TKE_adjust(hb, ha, Idecay) result(TKE_to_PE_scale)
  real, intent(in) :: hb   !< The thickness over which the buoyancy flux varies on the
                           !! near-boundary side of an interface (e.g., a well-mixed bottom
                           !! boundary layer thickness) [H ~> m or kg m-2]
  real, intent(in) :: ha   !< The thickness of the layer on the opposite side of an interface from
                           !! the boundary [H ~> m or kg m-2]
  real, intent(in) :: Idecay !< The inverse of a turbulence decay length scale [H-1 ~> m-1 or m2 kg-1]
  real             :: TKE_to_PE_scale !< The effective fractional change in energy available to
                           !! drive mixing at this interface once the exponential decay of TKE
                           !! is accounted for [nondim].  TKE_to_PE_scale is always positive.

  real :: khb  ! The thickness on the boundary side times the TKE decay rate [nondim]
  real :: kha  ! The thickness away from from the boundary times the TKE decay rate [nondim]
  real, parameter :: C1_3 = 1.0/3.0  ! A rational constant [nondim]

  khb = abs(hb*Idecay)
  kha = abs(ha*Idecay)

  !  For large enough kha that exp(kha) > 1.0e17*kha:
  !    TKE_to_PE_scale = (0.5 * (khb + kha) * kha) * exp(-kha) > (0.5 * kha**2) * exp(-kha)
  !  To keep TKE_to_PE_scale > -1e30 and avoid overflow in the exp(), keep kha < kha_max_30, where:
  !    kha_max_30 = ln(0.5*1e30) + 2.0 * ln(kha_max_30) ~= 68.3844 + 2.0 * ln(68.3844+8.6895))
  !    If kha_max = 77.0739, (0.5 * kha_max**2) * exp(-kha_max) = 1.0e-30.

  if (kha > 77.0739) then
    TKE_to_PE_scale = 1.0e-30
  elseif ((kha > 2.2e-4) .and. (khb > 2.2e-4)) then
    ! This is the usual case, derived from integrals of z exp(z) over the layers above and below.
    ! TKE_to_PE_scale = (0.5 * (khb + kha)) / &
    !                   ((exp(-khb) - (1.0 - khb)) / khb + (exp(kha) - (1.0 + kha)) / kha)
    TKE_to_PE_scale = (0.5 * (khb + kha) * (kha * khb)) / &
                      (kha * (exp(-khb) - (1.0 - khb)) + khb * (exp(kha) - (1.0 + kha)))
  elseif (khb > 2.2e-4) then
    ! For small values of kha, approximate (exp(kha) - (1.0 + hha)) by the first two
    ! terms of its Taylor series: 0.5*kha**2 + C1_6*kha**3 + ... + kha**n/n! + ...
    ! which is more accurate when kha**4/24. < 1e-16 or kha < ~ 2.21e-4.
    TKE_to_PE_scale = (0.5 * (khb + kha) * khb) / &
                      ((exp(-khb) - (1.0 - khb)) + 0.5*(khb * kha) * (1.0 + C1_3*kha))
  elseif (kha > 2.2e-4) then
    ! Use a Taylor series expansion for small values of khb
    TKE_to_PE_scale = (0.5 * (khb + kha) * kha) / &
                      (0.5 * (kha * khb) * (1.0 - C1_3*Khb) + (exp(kha) - (1.0 + kha)))
  else ! (kha < 2.2e-4) .and. (khb < 2.2e-4) - use Taylor series approximations for both
    TKE_to_PE_scale = 1.0 / (1.0 + C1_3*(kha - khb))
  endif

  if (TKE_to_PE_scale < 1.0e-30) TKE_to_PE_scale = 1.0e-30

  !  For kha >> 1:
  !    TKE_to_PE_scale = (0.5 * (khb + kha) * kha) * exp(-kha)

  !  For khb >> 1:
  !    TKE_to_PE_scale = (0.5 * (khb + kha) * (kha * khb)) / &
  !                      (khb * exp(kha) - (kha + khb)))
  !  For khb >> 1 and khb >> kha:
  !    TKE_to_PE_scale = (0.5 * (kha * khb)) / (exp(kha) - 1.0))

end function exp_decay_TKE_adjust

!> This subroutine calculates the change in potential energy and or derivatives
!! for several changes in an interface's diapycnal diffusivity times a timestep.
subroutine find_PE_chg(Kddt_h0, dKddt_h, hp_a, hp_b, Th_a, Sh_a, Th_b, Sh_b, &
                       dT_to_dPE_a, dS_to_dPE_a, dT_to_dPE_b, dS_to_dPE_b, &
                       pres_Z, dT_to_dColHt_a, dS_to_dColHt_a, dT_to_dColHt_b, dS_to_dColHt_b, &
                       PE_chg, dPEc_dKd, dPE_max, dPEc_dKd_0, PE_ColHt_cor)
  real, intent(in)  :: Kddt_h0  !< The previously used diffusivity at an interface times
                                !! the time step and divided by the average of the
                                !! thicknesses around the interface [H ~> m or kg m-2].
  real, intent(in)  :: dKddt_h  !< The trial change in the diffusivity at an interface times
                                !! the time step and divided by the average of the
                                !! thicknesses around the interface [H ~> m or kg m-2].
  real, intent(in)  :: hp_a     !< The effective pivot thickness of the layer above the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above [H ~> m or kg m-2].
  real, intent(in)  :: hp_b     !< The effective pivot thickness of the layer below the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface below [H ~> m or kg m-2].
  real, intent(in)  :: Th_a     !< An effective temperature times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_a     !< An effective salinity times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: Th_b     !< An effective temperature times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_b     !< An effective salinity times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: dT_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers above [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers above [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: dT_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers below [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers below [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: pres_Z   !< The rescaled hydrostatic interface pressure, which relates
                                !! the changes in column thickness to the energy that is radiated
                                !! as gravity waves and unavailable to drive mixing [R Z2 T-2 ~> J m-3].
  real, intent(in)  :: dT_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers above [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers above [Z S-1 ~> m ppt-1].
  real, intent(in)  :: dT_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers below [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers below [Z S-1 ~> m ppt-1].

  real, intent(out) :: PE_chg   !< The change in column potential energy from applying
                                !! dKddt_h at the present interface [R Z3 T-2 ~> J m-2].
  real, optional, intent(out) :: dPEc_dKd !< The partial derivative of PE_chg with dKddt_h
                                          !! [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real, optional, intent(out) :: dPE_max  !< The maximum change in column potential energy that could
                                          !! be realized by applying a huge value of dKddt_h at the
                                          !! present interface [R Z3 T-2 ~> J m-2].
  real, optional, intent(out) :: dPEc_dKd_0 !< The partial derivative of PE_chg with dKddt_h in the
                                            !! limit where dKddt_h = 0 [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real, optional, intent(out) :: PE_ColHt_cor !< The correction to PE_chg that is made due to a net
                                            !! change in the column height [R Z3 T-2 ~> J m-2].

  ! Local variables
  real :: hps ! The sum of the two effective pivot thicknesses [H ~> m or kg m-2].
  real :: bdt1 ! A product of the two pivot thicknesses plus a diffusive term [H2 ~> m2 or kg2 m-4].
  real :: dT_c ! The core term in the expressions for the temperature changes [C H2 ~> degC m2 or degC kg2 m-4].
  real :: dS_c ! The core term in the expressions for the salinity changes [S H2 ~> ppt m2 or ppt kg2 m-4].
  real :: PEc_core ! The diffusivity-independent core term in the expressions
                   ! for the potential energy changes [R Z2 T-2 ~> J m-3].
  real :: ColHt_core ! The diffusivity-independent core term in the expressions
                     ! for the column height changes [H Z ~> m2 or kg m-1].
  real :: ColHt_chg  ! The change in the column height [H ~> m or kg m-2].
  real :: y1_3 ! A local temporary term in [H-3 ~> m-3 or m6 kg-3].
  real :: y1_4 ! A local temporary term in [H-4 ~> m-4 or m8 kg-4].

  !   The expression for the change in potential energy used here is derived
  ! from the expression for the final estimates of the changes in temperature
  ! and salinities, and then extensively manipulated to get it into its most
  ! succinct form. The derivation is not necessarily obvious, but it demonstrably
  ! works by comparison with separate calculations of the energy changes after
  ! the tridiagonal solver for the final changes in temperature and salinity are
  ! applied.

  hps = hp_a + hp_b
  bdt1 = hp_a * hp_b + Kddt_h0 * hps
  dT_c = hp_a * Th_b - hp_b * Th_a
  dS_c = hp_a * Sh_b - hp_b * Sh_a
  PEc_core = hp_b * (dT_to_dPE_a * dT_c + dS_to_dPE_a * dS_c) - &
             hp_a * (dT_to_dPE_b * dT_c + dS_to_dPE_b * dS_c)
  ColHt_core = hp_b * (dT_to_dColHt_a * dT_c + dS_to_dColHt_a * dS_c) - &
               hp_a * (dT_to_dColHt_b * dT_c + dS_to_dColHt_b * dS_c)

  ! Find the change in column potential energy due to the change in the
  ! diffusivity at this interface by dKddt_h.
  y1_3 = dKddt_h / (bdt1 * (bdt1 + dKddt_h * hps))
  PE_chg = PEc_core * y1_3
  ColHt_chg = ColHt_core * y1_3
  if (ColHt_chg < 0.0) PE_chg = PE_chg - pres_Z * ColHt_chg

  if (present(PE_ColHt_cor)) PE_ColHt_cor = -pres_Z * min(ColHt_chg, 0.0)

  if (present(dPEc_dKd)) then
    ! Find the derivative of the potential energy change with dKddt_h.
    y1_4 = 1.0 / (bdt1 + dKddt_h * hps)**2
    dPEc_dKd = PEc_core * y1_4
    ColHt_chg = ColHt_core * y1_4
    if (ColHt_chg < 0.0) dPEc_dKd = dPEc_dKd - pres_Z * ColHt_chg
  endif

  if (present(dPE_max)) then
    ! This expression is the limit of PE_chg for infinite dKddt_h.
    y1_3 = 1.0 / (bdt1 * hps)
    dPE_max = PEc_core * y1_3
    ColHt_chg = ColHt_core * y1_3
    if (ColHt_chg < 0.0) dPE_max = dPE_max - pres_Z * ColHt_chg
  endif

  if (present(dPEc_dKd_0)) then
    ! This expression is the limit of dPEc_dKd for dKddt_h = 0.
    y1_4 = 1.0 / bdt1**2
    dPEc_dKd_0 = PEc_core * y1_4
    ColHt_chg = ColHt_core * y1_4
    if (ColHt_chg < 0.0) dPEc_dKd_0 = dPEc_dKd_0 - pres_Z * ColHt_chg
  endif

end subroutine find_PE_chg


!> This subroutine directly calculates the an increment in the diapycnal diffusivity based on the
!! change in potential energy within a timestep, subject to bounds on the possible change in
!! diffusivity, returning both the added diffusivity and the realized potential energy change, and
!! optionally also the maximum change in potential energy that would be realized for an infinitely
!! large diffusivity.
subroutine find_Kd_from_PE_chg(Kd_prev, dKd_max, dt_h, max_PE_chg, hp_a, hp_b, Th_a, Sh_a, Th_b, Sh_b, &
                       dT_to_dPE_a, dS_to_dPE_a, dT_to_dPE_b, dS_to_dPE_b, pres_Z, &
                       dT_to_dColHt_a, dS_to_dColHt_a, dT_to_dColHt_b, dS_to_dColHt_b, &
                       Kd_add, PE_chg, dPE_max, frac_dKd_max_PE)
  real, intent(in)  :: Kd_prev  !< The previously used diffusivity at an interface
                                !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, intent(in)  :: dKd_max  !< The maximum change in the diffusivity at an interface
                                !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, intent(in)  :: dt_h     !< The time step and divided by the average of the
                                !! thicknesses around the interface [T Z-1 ~> s m-1].
  real, intent(in)  :: max_PE_chg !< The maximum change in the column potential energy due to
                                !! additional mixing at an interface [R Z3 T-2 ~> J m-2].

  real, intent(in)  :: hp_a     !< The effective pivot thickness of the layer above the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above [H ~> m or kg m-2].
  real, intent(in)  :: hp_b     !< The effective pivot thickness of the layer below the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface below [H ~> m or kg m-2].
  real, intent(in)  :: Th_a     !< An effective temperature times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_a     !< An effective salinity times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: Th_b     !< An effective temperature times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_b     !< An effective salinity times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: dT_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers above [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers above [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: dT_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers below [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers below [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: pres_Z   !< The rescaled hydrostatic interface pressure, which relates
                                !! the changes in column thickness to the energy that is radiated
                                !! as gravity waves and unavailable to drive mixing [R Z2 T-2 ~> J m-3].
  real, intent(in)  :: dT_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers above [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers above [Z S-1 ~> m ppt-1].
  real, intent(in)  :: dT_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers below [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers below [Z S-1 ~> m ppt-1].
  real, intent(out) :: Kd_add   !< The additional diffusivity at an interface
                                !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, intent(out) :: PE_chg   !< The realized change in the column potential energy due to
                                !! additional mixing at an interface [R Z3 T-2 ~> J m-2].
  real, optional, &
        intent(out) :: dPE_max  !< The maximum change in column potential energy that could
                                !! be realized by applying a huge value of dKddt_h at the
                                !! present interface [R Z3 T-2 ~> J m-2].
  real, optional, &
        intent(out) :: frac_dKd_max_PE !< The fraction of the energy required to support dKd_max
                                !! that is supplied by max_PE_chg [nondim]

  ! Local variables
  real :: Kddt_h0   ! The previously used diffusivity at an interface times the time step
                    ! and divided by the average of the thicknesses around the
                    ! interface [H ~> m or kg m-2].
  real :: dKddt_h   ! The upper bound on the change in the diffusivity at an interface times
                    ! the time step and divided by the average of the thicknesses around
                    ! the interface [H ~> m or kg m-2].
  real :: hps       ! The sum of the two effective pivot thicknesses [H ~> m or kg m-2].
  real :: bdt1      ! A product of the two pivot thicknesses plus a diffusive term [H2 ~> m2 or kg2 m-4].
  real :: dT_c      ! The core term in the expressions for the temperature changes [C H2 ~> degC m2 or degC kg2 m-4].
  real :: dS_c      ! The core term in the expressions for the salinity changes [S H2 ~> ppt m2 or ppt kg2 m-4].
  real :: PEc_core  ! The diffusivity-independent core term in the expressions
                    ! for the potential energy changes [R Z2 T-2 ~> J m-3].
  real :: ColHt_core ! The diffusivity-independent core term in the expressions
                     ! for the column height changes [H Z ~> m2 or kg m-1].

  ! The expression for the change in potential energy used here is derived from the expression
  ! for the final estimates of the changes in temperature and salinities, which is then
  ! extensively manipulated to get it into its most succinct form.  It is the same as the
  ! expression that appears in find_PE_chg.

  Kddt_h0 = Kd_prev * dt_h
  hps = hp_a + hp_b
  bdt1 = hp_a * hp_b + Kddt_h0 * hps
  dT_c = hp_a * Th_b - hp_b * Th_a
  dS_c = hp_a * Sh_b - hp_b * Sh_a
  PEc_core = hp_b * (dT_to_dPE_a * dT_c + dS_to_dPE_a * dS_c) - &
             hp_a * (dT_to_dPE_b * dT_c + dS_to_dPE_b * dS_c)
  ColHt_core = hp_b * (dT_to_dColHt_a * dT_c + dS_to_dColHt_a * dS_c) - &
               hp_a * (dT_to_dColHt_b * dT_c + dS_to_dColHt_b * dS_c)
  if (ColHt_core < 0.0) PEc_core = PEc_core - pres_Z * ColHt_core

  ! Find the change in column potential energy due to the change in the
  ! diffusivity at this interface by dKd_max, and use this to dermine which limit applies.
  dKddt_h = dKd_max * dt_h
  if ( (PEc_core * dKddt_h <= max_PE_chg * (bdt1 * (bdt1 + dKddt_h * hps))) .or. (PEc_core <= 0.0) ) then
    ! There is more than enough energy available to support the maximum permitted diffusivity.
    Kd_add = dKd_max
    PE_chg = PEc_core * dKddt_h / (bdt1 * (bdt1 + dKddt_h * hps))
    if (present(frac_dKd_max_PE)) frac_dKd_max_PE = 1.0
  else
    ! Mixing is constrained by the available energy, so solve the following for Kd_add:
    !   max_PE_chg = PEc_core * Kd_add * dt_h / (bdt1 * (bdt1 + Kd_add * dt_h * hps))
    ! It has been verified that the two branches are continuous.
    Kd_add = (bdt1**2 * max_PE_chg) / (dt_h * (PEc_core - bdt1 * hps * max_PE_chg))
    PE_chg = max_PE_chg
    if (present(frac_dKd_max_PE)) &
      frac_dKd_max_PE = (PE_chg * (bdt1 * (bdt1 + dKddt_h * hps))) / (PEc_core * dKddt_h)
  endif

  ! Note that the derivative of PE_chg with dKddt_h is monotonic:
  !   dPE_chg_dKd = PEc_core * ( (bdt1 * (bdt1 + dKddt_h * hps)) - bdtl * hps * dKddt_h ) / &
  !                              (bdt1 * (bdt1 + dKddt_h * hps))**2
  !   dPE_chg_dKd = PEc_core / (bdt1 + dKddt_h * hps)**2

  ! This expression is the limit of PE_chg for infinite dKddt_h.
  if (present(dPE_max)) dPE_max = PEc_core / (bdt1 * hps)

end subroutine find_Kd_from_PE_chg


!> This subroutine calculates the change in potential energy and or derivatives
!! for several changes in an interface's diapycnal diffusivity times a timestep
!! using the original form used in the first version of ePBL.
subroutine find_PE_chg_orig(Kddt_h, h_k, b_den_1, dTe_term, dSe_term, &
                       dT_km1_t2, dS_km1_t2, dT_to_dPE_k, dS_to_dPE_k, &
                       dT_to_dPEa, dS_to_dPEa, pres_Z, dT_to_dColHt_k, &
                       dS_to_dColHt_k, dT_to_dColHta, dS_to_dColHta, PE_chg, &
                       dPEc_dKd, dPE_max, dPEc_dKd_0)
  real, intent(in)  :: Kddt_h   !< The diffusivity at an interface times the time step and
                                !! divided by the average of the thicknesses around the
                                !! interface [H ~> m or kg m-2].
  real, intent(in)  :: h_k      !< The thickness of the layer below the interface [H ~> m or kg m-2].
  real, intent(in)  :: b_den_1  !< The first term in the denominator of the pivot
                                !! for the tridiagonal solver, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above [H ~> m or kg m-2].
  real, intent(in)  :: dTe_term !< A diffusivity-independent term related to the temperature change
                                !! in the layer below the interface [C H ~> degC m or degC kg m-2].
  real, intent(in)  :: dSe_term !< A diffusivity-independent term related to the salinity change
                                !! in the layer below the interface [S H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: dT_km1_t2 !< A diffusivity-independent term related to the
                                 !! temperature change in the layer above the interface [C ~> degC].
  real, intent(in)  :: dS_km1_t2 !< A diffusivity-independent term related to the
                                 !! salinity change in the layer above the interface [S ~> ppt].
  real, intent(in)  :: pres_Z    !< The rescaled hydrostatic interface pressure, which relates
                                 !! the changes in column thickness to the energy that is radiated
                                 !! as gravity waves and unavailable to drive mixing [R Z2 T-2 ~> J m-3].
  real, intent(in)  :: dT_to_dPE_k !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                 !! a layer's temperature change to the change in column potential
                                 !! energy, including all implicit diffusive changes in the
                                 !! temperatures of all the layers below [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_k !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                 !! a layer's salinity change to the change in column potential
                                 !! energy, including all implicit diffusive changes in the
                                 !! in the salinities of all the layers below [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: dT_to_dPEa !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                 !! a layer's temperature change to the change in column potential
                                 !! energy, including all implicit diffusive changes in the
                                 !! temperatures of all the layers above [R Z3 T-2 C-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPEa !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                 !! a layer's salinity change to the change in column potential
                                 !! energy, including all implicit diffusive changes in the
                                 !! salinities of all the layers above [R Z3 T-2 S-1 ~> J m-2 ppt-1].
  real, intent(in)  :: dT_to_dColHt_k !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! height, including all implicit diffusive changes in the
                                 !! temperatures of all the layers below [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_k !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the salinities of all the layers below [Z S-1 ~> m ppt-1].
  real, intent(in)  :: dT_to_dColHta !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the temperatures of all the layers above [Z C-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHta !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the salinities of all the layers above [Z S-1 ~> m ppt-1].

  real, intent(out) :: PE_chg    !< The change in column potential energy from applying
                                 !! Kddt_h at the present interface [R Z3 T-2 ~> J m-2].
  real, optional, intent(out) :: dPEc_dKd !< The partial derivative of PE_chg with Kddt_h
                                          !! [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real, optional, intent(out) :: dPE_max  !< The maximum change in column potential energy that could
                                          !! be realized by applying a huge value of Kddt_h at the
                                          !! present interface [R Z3 T-2 ~> J m-2].
  real, optional, intent(out) :: dPEc_dKd_0 !< The partial derivative of PE_chg with Kddt_h in the
                                          !! limit where Kddt_h = 0 [R Z3 T-2 H-1 ~> J m-3 or J kg-1].

!   This subroutine determines the total potential energy change due to mixing
! at an interface, including all of the implicit effects of the prescribed
! mixing at interfaces above.  Everything here is derived by careful manipulation
! of the robust tridiagonal solvers used for tracers by MOM6.  The results are
! positive for mixing in a stably stratified environment.
!   The comments describing these arguments are for a downward mixing pass, but
! this routine can also be used for an upward pass with the sense of direction
! reversed.

  ! Local variables
  real :: b1            ! b1 is used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
  real :: b1Kd          ! Temporary array [nondim]
  real :: ColHt_chg     ! The change in column thickness [Z ~> m].
  real :: dColHt_max    ! The change in column thickness for infinite diffusivity [Z ~> m].
  real :: dColHt_dKd    ! The partial derivative of column thickness with Kddt_h [Z H-1 ~> nondim or m3 kg-1]
  real :: dT_k, dT_km1  ! Temperature changes in layers k and k-1 [C ~> degC]
  real :: dS_k, dS_km1  ! Salinity changes in layers k and k-1 [S ~> ppt]
  real :: I_Kr_denom    ! Temporary array [H-2 ~> m-2 or m4 kg-2]
  real :: dKr_dKd       ! Temporary array [H-2 ~> m-2 or m4 kg-2]
  real :: ddT_k_dKd, ddT_km1_dKd ! Temporary arrays indicating the temperature changes
                        ! per unit change in Kddt_h [C H-1 ~> degC m-1 or degC m2 kg-1]
  real :: ddS_k_dKd, ddS_km1_dKd ! Temporary arrays indicating the salinity changes
                        ! per unit change in Kddt_h [S H-1 ~> ppt m-1 or ppt m2 kg-1]

  b1 = 1.0 / (b_den_1 + Kddt_h)
  b1Kd = Kddt_h*b1

  ! Start with the temperature change in layer k-1 due to the diffusivity at
  ! interface K without considering the effects of changes in layer k.

  ! Calculate the change in PE due to the diffusion at interface K
  ! if Kddt_h(K+1) = 0.
  I_Kr_denom = 1.0 / (h_k*b_den_1 + (b_den_1 + h_k)*Kddt_h)

  dT_k = (Kddt_h*I_Kr_denom) * dTe_term
  dS_k = (Kddt_h*I_Kr_denom) * dSe_term

  ! Find the change in energy due to diffusion with strength Kddt_h at this interface.
  ! Increment the temperature changes in layer k-1 due the changes in layer k.
  dT_km1 = b1Kd * ( dT_k + dT_km1_t2 )
  dS_km1 = b1Kd * ( dS_k + dS_km1_t2 )
  PE_chg = (dT_to_dPE_k * dT_k + dT_to_dPEa * dT_km1) + &
           (dS_to_dPE_k * dS_k + dS_to_dPEa * dS_km1)
  ColHt_chg = (dT_to_dColHt_k * dT_k + dT_to_dColHta * dT_km1) + &
              (dS_to_dColHt_k * dS_k + dS_to_dColHta * dS_km1)
  if (ColHt_chg < 0.0) PE_chg = PE_chg - pres_Z * ColHt_chg

  if (present(dPEc_dKd)) then
    ! Find the derivatives of the temperature and salinity changes with Kddt_h.
    dKr_dKd = (h_k*b_den_1) * I_Kr_denom**2

    ddT_k_dKd = dKr_dKd * dTe_term
    ddS_k_dKd = dKr_dKd * dSe_term
    ddT_km1_dKd = (b1**2 * b_den_1) * ( dT_k + dT_km1_t2 ) + b1Kd * ddT_k_dKd
    ddS_km1_dKd = (b1**2 * b_den_1) * ( dS_k + dS_km1_t2 ) + b1Kd * ddS_k_dKd

    ! Calculate the partial derivative of Pe_chg with Kddt_h.
    dPEc_dKd = (dT_to_dPE_k * ddT_k_dKd + dT_to_dPEa * ddT_km1_dKd) + &
               (dS_to_dPE_k * ddS_k_dKd + dS_to_dPEa * ddS_km1_dKd)
    dColHt_dKd = (dT_to_dColHt_k * ddT_k_dKd + dT_to_dColHta * ddT_km1_dKd) + &
                 (dS_to_dColHt_k * ddS_k_dKd + dS_to_dColHta * ddS_km1_dKd)
    if (dColHt_dKd < 0.0) dPEc_dKd = dPEc_dKd - pres_Z * dColHt_dKd
  endif

  if (present(dPE_max)) then
    ! This expression is the limit of PE_chg for infinite Kddt_h.
    dPE_max = (dT_to_dPEa * dT_km1_t2 + dS_to_dPEa * dS_km1_t2) + &
              ((dT_to_dPE_k + dT_to_dPEa) * dTe_term + &
               (dS_to_dPE_k + dS_to_dPEa) * dSe_term) / (b_den_1 + h_k)
    dColHt_max = (dT_to_dColHta * dT_km1_t2 + dS_to_dColHta * dS_km1_t2) + &
              ((dT_to_dColHt_k + dT_to_dColHta) * dTe_term + &
               (dS_to_dColHt_k + dS_to_dColHta) * dSe_term) / (b_den_1 + h_k)
    if (dColHt_max < 0.0) dPE_max = dPE_max - pres_Z*dColHt_max
  endif

  if (present(dPEc_dKd_0)) then
    ! This expression is the limit of dPEc_dKd for Kddt_h = 0.
    dPEc_dKd_0 = (dT_to_dPEa * dT_km1_t2 + dS_to_dPEa * dS_km1_t2) / (b_den_1) + &
                 (dT_to_dPE_k * dTe_term + dS_to_dPE_k * dSe_term) / (h_k*b_den_1)
    dColHt_dKd = (dT_to_dColHta * dT_km1_t2 + dS_to_dColHta * dS_km1_t2) / (b_den_1) + &
                 (dT_to_dColHt_k * dTe_term + dS_to_dColHt_k * dSe_term) / (h_k*b_den_1)
    if (dColHt_dKd < 0.0) dPEc_dKd_0 = dPEc_dKd_0 - pres_Z*dColHt_dKd
  endif

end subroutine find_PE_chg_orig

!> This subroutine finds the mstar value for ePBL
subroutine find_mstar(CS, US, Buoyancy_Flux, UStar, &
                      BLD, Abs_Coriolis, Is_BBL, mstar, &
                      Langmuir_Number, mstar_LT, Convect_Langmuir_Number)
  type(energetic_PBL_CS), intent(in) :: CS    !< Energetic PBL control structure
  type(unit_scale_type), intent(in)  :: US    !< A dimensional unit scaling type
  real,                  intent(in)  :: UStar !< ustar including gustiness [Z T-1 ~> m s-1]
  real,                  intent(in)  :: Abs_Coriolis !< absolute value of the Coriolis parameter [T-1 ~> s-1]
  real,                  intent(in)  :: Buoyancy_Flux !< Buoyancy flux [Z2 T-3 ~> m2 s-3]
  real,                  intent(in)  :: BLD   !< boundary layer depth [Z ~> m]
  logical,               intent(in)  :: Is_BBL !< Logcal flag to indicate if bottom boundary layer mode
  real,                  intent(out) :: mstar !< Output mstar (Mixing/ustar**3) [nondim]
  real,        optional, intent(in)  :: Langmuir_Number !< Langmuir number [nondim]
  real,        optional, intent(out) :: mstar_LT !< mstar increase due to Langmuir turbulence [nondim]
  real,        optional, intent(out) :: Convect_Langmuir_number !< Langmuir number including buoyancy flux [nondim]

  !/ Variables used in computing mstar
  real :: MSN_term       ! Temporary terms [nondim]
  real :: MSCR_term1, MSCR_term2 ! Temporary terms [Z3 T-3 ~> m3 s-3]
  real :: mstar_Conv_Red ! Adjustment made to mstar due to convection reducing mechanical mixing [nondim]
  real :: mstar_S, mstar_N ! mstar in (S)tabilizing/(N)ot-stabilizing buoyancy flux [nondim]
  integer :: mstar_scheme ! Toggles between surface and bottom boundary layer mstar scheme from control structure

  !/  Integer options for how to find mstar

  !/

  if (Is_BBL) then
    mstar_scheme = CS%BBL_mstar_scheme
  else
    mstar_scheme = CS%mstar_scheme
  endif

  if (mstar_scheme == Use_Fixed_mstar) then
    if (Is_BBL) then
      mstar = CS%BBL_Fixed_mstar
    else
      mstar = CS%Fixed_mstar
    endif
  !/ 1. Get mstar
  elseif (mstar_scheme == mstar_from_Ekman) then

    if (CS%answer_date < 20190101) then
      ! The limit for the balance of rotation and stabilizing is f(L_Ekman,L_Obukhov)
      mstar_S = CS%mstar_coef*sqrt(max(0.0,Buoyancy_Flux) / UStar**2 / &
                    (Abs_Coriolis + 1.e-10*US%T_to_s) )
      ! The limit for rotation (Ekman length) limited mixing
      mstar_N =  CS%C_Ek * log( max( 1., UStar / (Abs_Coriolis + 1.e-10*US%T_to_s) / BLD ) )
    else
      ! The limit for the balance of rotation and stabilizing is f(L_Ekman,L_Obukhov)
      mstar_S = CS%mstar_coef*sqrt(max(0.0, Buoyancy_Flux) / (UStar**2 * max(Abs_Coriolis, 1.e-20*US%T_to_s)))
      ! The limit for rotation (Ekman length) limited mixing
      mstar_N = 0.0
      if (UStar > Abs_Coriolis * BLD) mstar_N = CS%C_Ek * log(UStar / (Abs_Coriolis * BLD))
    endif

    ! Here 1.25 is about .5/von Karman, which gives the Obukhov limit.
    mstar = max(mstar_S, min(1.25, mstar_N))
    if (CS%mstar_Cap > 0.0) mstar = min( CS%mstar_Cap,mstar )
  elseif ( mstar_scheme == mstar_from_RH18 ) then
    if (CS%answer_date < 20190101) then
      mstar_N = CS%RH18_mstar_cn1 * ( 1.0 - 1.0 / ( 1. + CS%RH18_mstar_cn2 * &
                exp( CS%RH18_mstar_CN3 * BLD * Abs_Coriolis / UStar) ) )
    else
      MSN_term = CS%RH18_mstar_cn2 * exp( CS%RH18_mstar_CN3 * BLD * Abs_Coriolis / UStar)
      mstar_N = (CS%RH18_mstar_cn1 *  MSN_term) / ( 1. + MSN_term)
    endif
    mstar_S = CS%RH18_mstar_CS1 * ( max(0.0, Buoyancy_Flux)**2 * BLD / &
             ( UStar**5 * max(Abs_Coriolis,1.e-20*US%T_to_s) ) )**CS%RH18_mstar_cs2
    mstar = mstar_N + mstar_S
  endif

  !/ 2. Adjust mstar to account for convective turbulence
  if (CS%answer_date < 20190101) then
    mstar_Conv_Red = 1. - CS%mstar_Convect_coef * (-min(0.0,Buoyancy_Flux) + 1.e-10*US%T_to_s**3*US%m_to_Z**2) / &
                         ( (-min(0.0,Buoyancy_Flux) + 1.e-10*US%T_to_s**3*US%m_to_Z**2) + &
                         2.0 *mstar * UStar**3 / BLD )
  else
    MSCR_term1 = -BLD * min(0.0, Buoyancy_Flux)
    MSCR_term2 = 2.0*mstar * UStar**3
    if ( abs(MSCR_term2) > 0.0) then
      mstar_Conv_Red = ((1.-CS%mstar_convect_coef) * MSCR_term1 + MSCR_term2) / (MSCR_term1 + MSCR_term2)
    else
      mstar_Conv_Red = 1.-CS%mstar_convect_coef
    endif
  endif

  !/3. Combine various mstar terms to get final value
  mstar = mstar * mstar_Conv_Red

  if ((.not.Is_BBL) .and. (present(Langmuir_Number))) then
    call mstar_Langmuir(CS, US, Abs_Coriolis, Buoyancy_Flux, UStar, BLD, Langmuir_Number, mstar, &
                        mstar_LT, Convect_Langmuir_Number)
  endif

end subroutine Find_mstar

!> This subroutine modifies the mstar value if the Langmuir number is present
subroutine mstar_Langmuir(CS, US, Abs_Coriolis, Buoyancy_Flux, UStar, BLD, Langmuir_Number, &
                          mstar, mstar_LT, Convect_Langmuir_Number)
  type(energetic_PBL_CS), intent(in) :: CS    !< Energetic PBL control structure
  type(unit_scale_type), intent(in)  :: US    !< A dimensional unit scaling type
  real,                  intent(in)  :: Abs_Coriolis !< Absolute value of the Coriolis parameter [T-1 ~> s-1]
  real,                  intent(in)  :: Buoyancy_Flux !< Buoyancy flux [Z2 T-3 ~> m2 s-3]
  real,                  intent(in)  :: UStar !< Surface friction velocity with? gustiness [Z T-1 ~> m s-1]
  real,                  intent(in)  :: BLD   !< boundary layer depth [Z ~> m]
  real,                  intent(inout) :: mstar !< Input/output mstar (Mixing/ustar**3) [nondim]
  real,                  intent(in)  :: Langmuir_Number !< Langmuir number [nondim]
  real,                  intent(out) :: mstar_LT !< mstar increase due to Langmuir turbulence [nondim]
  real,                  intent(out) :: Convect_Langmuir_number !< Langmuir number including buoyancy flux [nondim]

  !/
  real, parameter :: Max_ratio = 1.0e16  ! The maximum value of a nondimensional ratio [nondim].
  real :: enhance_mstar ! A multiplicative scaling of mstar due to Langmuir turbulence [nondim].
  real :: mstar_LT_add ! A value that is added to mstar due to Langmuir turbulence [nondim].
  real :: iL_Ekman    ! Inverse of Ekman length scale [Z-1 ~> m-1].
  real :: iL_Obukhov  ! Inverse of Obukhov length scale [Z-1 ~> m-1].
  real :: I_ustar     ! The Adcroft reciprocal of ustar [T Z-1 ~> s m-1]
  real :: I_f         ! The Adcroft reciprocal of the Coriolis parameter [T ~> s]
  real :: MLD_Ekman          ! The ratio of the mixed layer depth to the Ekman layer depth [nondim].
  real :: Ekman_Obukhov      ! The Ekman layer thickness divided by the Obukhov depth [nondim].
  real :: MLD_Obukhov        ! The mixed layer depth divided by the Obukhov depth [nondim].
  real :: MLD_Obukhov_stab   ! The mixed layer depth divided by the Obukhov depth under stable
                             ! conditions or 0 under unstable conditions [nondim].
  real :: Ekman_Obukhov_stab ! The Ekman layer thickness divided by the Obukhov depth under stable
                             ! conditions or 0 under unstable conditions [nondim].
  real :: MLD_Obukhov_un     ! The mixed layer depth divided by the Obukhov depth under unstable
                             ! conditions or 0 under stable conditions [nondim].
  real :: Ekman_Obukhov_un   ! The Ekman layer thickness divided by the Obukhov depth under unstable
                             ! conditions or 0 under stable conditions [nondim].

  ! Set default values for no Langmuir effects.
  enhance_mstar = 1.0 ; mstar_LT_add = 0.0

  if (CS%LT_enhance_form /= No_Langmuir) then
    ! a. Get parameters for modified LA
    if (CS%answer_date < 20190101) then
      iL_Ekman   = Abs_Coriolis / Ustar
      iL_Obukhov = Buoyancy_Flux*CS%vonkar / Ustar**3
      Ekman_Obukhov_stab = abs(max(0., iL_Obukhov / (iL_Ekman + 1.e-10*US%Z_to_m)))
      Ekman_Obukhov_un = abs(min(0., iL_Obukhov / (iL_Ekman + 1.e-10*US%Z_to_m)))
      MLD_Obukhov_stab = abs(max(0., BLD*iL_Obukhov))
      MLD_Obukhov_un = abs(min(0., BLD*iL_Obukhov))
      MLD_Ekman = abs( BLD*iL_Ekman )
    else
      Ekman_Obukhov = Max_ratio ; MLD_Obukhov = Max_ratio ; MLD_Ekman = Max_ratio
      I_f = 0.0 ; if (abs(abs_Coriolis) > 0.0) I_f = 1.0 / abs_Coriolis
      I_ustar = 0.0 ; if (abs(Ustar) > 0.0) I_ustar = 1.0 / Ustar
      if (abs(Buoyancy_Flux*CS%vonkar) < Max_ratio*(abs_Coriolis * Ustar**2)) &
        Ekman_Obukhov = abs(Buoyancy_Flux*CS%vonkar) * (I_f * I_Ustar**2)
      if (abs(BLD*Buoyancy_Flux*CS%vonkar) < Max_ratio*Ustar**3) &
        MLD_Obukhov = abs(BLD*Buoyancy_Flux*CS%vonkar) * I_Ustar**3
      if (BLD*Abs_Coriolis < Max_ratio*Ustar) &
        MLD_Ekman = BLD*Abs_Coriolis * I_Ustar

      if (Buoyancy_Flux > 0.0) then
        Ekman_Obukhov_stab = Ekman_Obukhov ; Ekman_Obukhov_un = 0.0
        MLD_Obukhov_stab = MLD_Obukhov ; MLD_Obukhov_un = 0.0
      else
        Ekman_Obukhov_un = Ekman_Obukhov ; Ekman_Obukhov_stab = 0.0
        MLD_Obukhov_un = MLD_Obukhov ; MLD_Obukhov_stab = 0.0
      endif
    endif

    ! b. Adjust LA based on various parameters.
    !    Assumes linear factors based on length scale ratios to adjust LA
    !    Note when these coefficients are set to 0 recovers simple LA.
    Convect_Langmuir_Number = Langmuir_Number * &
                    ( (1.0 + max(-0.5, CS%LaC_MLD_Ek * MLD_Ekman)) + &
                   ((CS%LaC_Ek_Ob_stab * Ekman_Obukhov_stab + CS%LaC_Ek_Ob_un * Ekman_Obukhov_un) + &
                    (CS%LaC_MLD_Ob_stab * MLD_Obukhov_stab  + CS%LaC_MLD_Ob_un * MLD_Obukhov_un)) )

    if (CS%LT_enhance_form == Langmuir_rescale) then
      ! Enhancement is multiplied (added mst_lt set to 0)
      Enhance_mstar = min(CS%Max_Enhance_M, &
                          (1. + CS%LT_enhance_coef * Convect_Langmuir_Number**CS%LT_enhance_exp) )
    elseif (CS%LT_enhance_form == Langmuir_add) then
      ! or Enhancement is additive (multiplied enhance_m set to 1)
      mstar_LT_add = CS%LT_enhance_coef * Convect_Langmuir_Number**CS%LT_enhance_exp
    endif
  endif

  mstar_LT = (enhance_mstar - 1.0)*mstar + mstar_LT_add  ! Diagnose the full increase in mstar.
  mstar = mstar*enhance_mstar + mstar_LT_add

end subroutine mstar_Langmuir


!> Copies the ePBL active mixed layer depth into MLD, in units of [Z ~> m] unless other units are specified.
subroutine energetic_PBL_get_MLD(CS, MLD, G, US, m_to_MLD_units)
  type(energetic_PBL_CS),           intent(in)  :: CS  !< Energetic PBL control structure
  type(ocean_grid_type),            intent(in)  :: G   !< Grid structure
  type(unit_scale_type),            intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: MLD !< Depth of ePBL active mixing layer [Z ~> m]
                                                       !! or other units
  real,                   optional, intent(in)  :: m_to_MLD_units !< A conversion factor from meters
                                                       !! to the desired units for MLD, sometimes [Z m-1 ~> 1]
  ! Local variables
  real :: scale  ! A dimensional rescaling factor, often [nondim] or [m Z-1 ~> 1]
  integer :: i, j

  scale = 1.0 ; if (present(m_to_MLD_units)) scale = US%Z_to_m * m_to_MLD_units

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    MLD(i,j) = scale*CS%ML_depth(i,j)
  enddo ; enddo

end subroutine energetic_PBL_get_MLD


!> This subroutine initializes the energetic_PBL module
subroutine energetic_PBL_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< The current model time
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(energetic_PBL_CS),  intent(inout) :: CS   !< Energetic PBL control structure

  ! Local variables
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_energetic_PBL"  ! This module's name.
  character(len=20)  :: tmpstr  ! A string that is parsed for parameter settings
  character(len=20)  :: mstar_scheme ! A string that is parsed for mstar parameter settings
  character(len=20)  :: vel_scale_str ! A string that is parsed for velocity scale parameter settings
  character(len=120) :: diff_text ! A clause describing parameter setting that differ.
  real :: omega_frac_dflt  ! The default for omega_frac [nondim]
  integer :: isd, ied, jsd, jed
  integer :: mstar_mode, LT_enhance, wT_mode
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  logical :: enable_bugs  ! If true, the defaults for recently added bug-fix flags are set to
                          ! recreate the bugs, or if false bugs are only used if actively selected.
  logical :: use_omega
  logical :: no_BBL  ! If true, EPBL_BBL_EFFIC < 0 and EPBL_BBL_TIDAL_EFFIC < 0, so
                     ! bottom boundary layer mixing is not enabled.
  logical :: use_la_windsea
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%initialized = .true.
  CS%diag => diag
  CS%Time => Time

! Set default, read and log parameters
  call log_version(param_file, mdl, version, "")


!/1. General ePBL settings
  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "OMEGA", CS%omega, &
                 "The rotation rate of the earth.", &
                 units="s-1", default=7.2921e-5, scale=US%T_to_S)
  call get_param(param_file, mdl, "ML_USE_OMEGA", use_omega, &
                 "If true, use the absolute rotation rate instead of the "//&
                 "vertical component of rotation when setting the decay "//&
                 "scale for turbulence.", default=.false., do_not_log=.true.)
  omega_frac_dflt = 0.0
  if (use_omega) then
    call MOM_error(WARNING, "ML_USE_OMEGA is deprecated; use ML_OMEGA_FRAC=1.0 instead.")
    omega_frac_dflt = 1.0
  endif
  call get_param(param_file, mdl, "ML_OMEGA_FRAC", CS%omega_frac, &
                 "When setting the decay scale for turbulence, use this "//&
                 "fraction of the absolute rotation rate blended with the "//&
                 "local value of f, as sqrt((1-of)*f^2 + of*4*omega^2).", &
                 units="nondim", default=omega_frac_dflt)
  call get_param(param_file, mdl, "EKMAN_SCALE_COEF", CS%Ekman_scale_coef, &
                 "A nondimensional scaling factor controlling the inhibition "//&
                 "of the diffusive length scale by rotation. Making this larger "//&
                 "decreases the PBL diffusivity.", units="nondim", default=1.0)
  call get_param(param_file, mdl, 'VON_KARMAN_CONST', CS%vonKar, &
                 'The value the von Karman constant as used for mixed layer viscosity.', &
                 units='nondim', default=0.41)
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "EPBL_ANSWER_DATE", CS%answer_date, &
                 "The vintage of the order of arithmetic and expressions in the energetic "//&
                 "PBL calculations.  Values below 20190101 recover the answers from the "//&
                 "end of 2018, while higher values use updated and more robust forms of the "//&
                 "same expressions.  Values below 20240101 use A**(1./3.) to estimate the cube "//&
                 "root of A in several expressions, while higher values use the integer root "//&
                 "function cuberoot(A) and therefore can work with scaled variables.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) CS%answer_date = max(CS%answer_date, 20230701)

  call get_param(param_file, mdl, "EPBL_ORIGINAL_PE_CALC", CS%orig_PE_calc, &
                 "If true, the ePBL code uses the original form of the potential energy change "//&
                 "code.  Otherwise, the newer version that can work with successive increments "//&
                 "to the diffusivity in upward or downward passes is used.", &
                 default=.true.) ! Change the default to .false.?

  call get_param(param_file, mdl, "MKE_TO_TKE_EFFIC", CS%MKE_to_TKE_effic, &
                 "The efficiency with which mean kinetic energy released "//&
                 "by mechanically forced entrainment of the mixed layer "//&
                 "is converted to turbulent kinetic energy.", &
                 units="nondim", default=0.0, scale=US%L_to_Z**2)
  call get_param(param_file, mdl, "TKE_DECAY", CS%TKE_decay, &
                 "TKE_DECAY relates the vertical rate of decay of the TKE available "//&
                 "for mechanical entrainment to the natural Ekman depth.", &
                 units="nondim", default=2.5)
  call get_param(param_file, mdl, "DIRECT_EPBL_MIXING_CALC", CS%direct_calc, &
                 "If true and there is no conversion from mean kinetic energy to ePBL turbulent "//&
                 "kinetic energy, use a direct calculation of the diffusivity that is supported "//&
                 "by a given energy input instead of the more general but slower iterative solver.", &
                 default=.false., do_not_log=(CS%MKE_to_TKE_effic>0.0))


!/2. Options related to setting mstar

  call get_param(param_file, mdl, "EPBL_MSTAR_SCHEME", mstar_scheme, &
                 "EPBL_MSTAR_SCHEME selects the method for setting mstar.  Valid values are: \n"//&
                 "\t CONSTANT   - Use a fixed mstar given by MSTAR \n"//&
                 "\t OM4        - Use L_Ekman/L_Obukhov in the stabilizing limit, as in OM4 \n"//&
                 "\t REICHL_H18 - Use the scheme documented in Reichl & Hallberg, 2018.", &
                 default=CONSTANT_STRING, do_not_log=.true.)
  call get_param(param_file, mdl, "MSTAR_MODE", mstar_mode, default=-1)
  if (mstar_mode == 0) then
    mstar_scheme = CONSTANT_STRING
    call MOM_error(WARNING, "Use EPBL_MSTAR_SCHEME = CONSTANT instead of the archaic MSTAR_MODE = 0.")
  elseif (mstar_mode == 1) then
    call MOM_error(FATAL, "You are using a legacy mstar mode in ePBL that has been phased out. "//&
                          "If you need to use this setting please report this error.  Also use "//&
                          "EPBL_MSTAR_SCHEME to specify the scheme for mstar.")
  elseif (mstar_mode == 2) then
    mstar_scheme = OM4_STRING
    call MOM_error(WARNING, "Use EPBL_MSTAR_SCHEME = OM4 instead of the archaic MSTAR_MODE = 2.")
  elseif (mstar_mode == 3) then
    mstar_scheme = RH18_STRING
    call MOM_error(WARNING, "Use EPBL_MSTAR_SCHEME = REICHL_H18 instead of the archaic MSTAR_MODE = 3.")
  elseif (mstar_mode > 3) then
    call MOM_error(FATAL, "An unrecognized value of the obsolete parameter MSTAR_MODE was specified.")
  endif
  call log_param(param_file, mdl, "EPBL_MSTAR_SCHEME", mstar_scheme, &
                 "EPBL_MSTAR_SCHEME selects the method for setting mstar.  Valid values are: \n"//&
                 "\t CONSTANT   - Use a fixed mstar given by MSTAR \n"//&
                 "\t OM4        - Use L_Ekman/L_Obukhov in the stabilizing limit, as in OM4 \n"//&
                 "\t REICHL_H18 - Use the scheme documented in Reichl & Hallberg, 2018.", &
                 default=CONSTANT_STRING)
  mstar_scheme = uppercase(mstar_scheme)
  select case (mstar_scheme)
    case (CONSTANT_STRING)
      CS%mstar_scheme = Use_Fixed_mstar
    case (OM4_STRING)
      CS%mstar_scheme = mstar_from_Ekman
    case (RH18_STRING)
      CS%mstar_scheme = mstar_from_RH18
    case default
      call MOM_mesg('energetic_PBL_init: EPBL_MSTAR_SCHEME ="'//trim(mstar_scheme)//'"', 0)
      call MOM_error(FATAL, "energetic_PBL_init: Unrecognized setting "// &
            "EPBL_MSTAR_SCHEME = "//trim(mstar_scheme)//" found in input file.")
  end select
  call get_param(param_file, mdl, "MSTAR", CS%fixed_mstar, &
                 "The ratio of the friction velocity cubed to the TKE input to the "//&
                 "surface boundary layer.  This option is used if EPBL_MSTAR_SCHEME = CONSTANT.", &
                 units="nondim", default=1.2, do_not_log=(CS%mstar_scheme/=Use_Fixed_mstar))

  call get_param(param_file, mdl, "MSTAR_CAP", CS%mstar_cap, &
                 "If this value is positive, it sets the maximum value of mstar "//&
                 "allowed in ePBL.  (This is not used if EPBL_mstar_scheme = CONSTANT).", &
                 units="nondim", default=-1.0, do_not_log=(CS%mstar_scheme==Use_Fixed_mstar))
  ! mstar_scheme==mstar_from_Ekman options
  call get_param(param_file, mdl, "MSTAR2_COEF1", CS%mstar_coef, &
                 "Coefficient in computing mstar when rotation and stabilizing "//&
                 "effects are both important (used if EPBL_mstar_scheme = OM4).", &
                 units="nondim", default=0.3, do_not_log=(CS%mstar_scheme/=mstar_from_Ekman))
  call get_param(param_file, mdl, "MSTAR2_COEF2", CS%C_Ek, &
                 "Coefficient in computing mstar when only rotation limits "// &
                 "the total mixing (used if EPBL_MSTAR_SCHEME = OM4)", &
                 units="nondim", default=0.085, do_not_log=(CS%mstar_scheme/=mstar_from_Ekman))
  ! mstar_scheme==mstar_from_RH18 options
  call get_param(param_file, mdl, "RH18_MSTAR_CN1", CS%RH18_mstar_cn1,&
                 "MSTAR_N coefficient 1 (outer-most coefficient for fit). "//&
                 "The value of 0.275 is given in RH18.  Increasing this "//&
                 "coefficient increases mstar for all values of Hf/ust, but more "//&
                 "effectively at low values (weakly developed OSBLs).", &
                 units="nondim", default=0.275, do_not_log=(CS%mstar_scheme/=mstar_from_RH18))
  call get_param(param_file, mdl, "RH18_MSTAR_CN2", CS%RH18_mstar_cn2,&
                 "MSTAR_N coefficient 2 (coefficient outside of exponential decay). "//&
                 "The value of 8.0 is given in RH18.  Increasing this coefficient "//&
                 "increases mstar for all values of HF/ust, with a much more even "//&
                 "effect across a wide range of Hf/ust than CN1.", &
                 units="nondim", default=8.0, do_not_log=(CS%mstar_scheme/=mstar_from_RH18))
  call get_param(param_file, mdl, "RH18_MSTAR_CN3", CS%RH18_mstar_CN3,&
                 "MSTAR_N coefficient 3 (exponential decay coefficient). "//&
                 "The value of -5.0 is given in RH18.  Increasing this increases how "//&
                 "quickly the value of mstar decreases as Hf/ust increases.", &
                  units="nondim", default=-5.0, do_not_log=(CS%mstar_scheme/=mstar_from_RH18))
  call get_param(param_file, mdl, "RH18_MSTAR_CS1", CS%RH18_mstar_cs1,&
                 "MSTAR_S coefficient for RH18 in stabilizing limit. "//&
                 "The value of 0.2 is given in RH18 and increasing it increases "//&
                 "mstar in the presence of a stabilizing surface buoyancy flux.", &
                 units="nondim", default=0.2, do_not_log=(CS%mstar_scheme/=mstar_from_RH18))
  call get_param(param_file, mdl, "RH18_MSTAR_CS2", CS%RH18_mstar_cs2,&
                 "MSTAR_S exponent for RH18 in stabilizing limit. "//&
                 "The value of 0.4 is given in RH18 and increasing it increases mstar "//&
                 "exponentially in the presence of a stabilizing surface buoyancy flux.", &
                 Units="nondim", default=0.4, do_not_log=(CS%mstar_scheme/=mstar_from_RH18))
!/ BBL mstar related options
  call get_param(param_file, mdl, "EPBL_BBL_USE_MSTAR", CS%ePBL_BBL_use_mstar, &
                 "A logical to use mstar in the calculation of TKE in the ePBL BBL scheme", &
                 units="nondim", default=.false.)
  if (CS%ePBL_BBL_use_mstar) then
    call get_param(param_file, mdl, "EPBL_BBL_MSTAR_SCHEME", tmpstr, &
                   "EPBL_BBL_MSTAR_SCHEME selects the method for setting mstar in the BBL.  Valid values are: \n"//&
                   "\t CONSTANT   - Use a fixed mstar given by MSTAR_BBL \n"//&
                   "\t OM4        - Use L_Ekman/L_Obukhov in the stabilizing limit, as in OM4 \n"//&
                   "\t REICHL_H18 - Use the scheme documented in Reichl & Hallberg, 2018.", &
                   default=mstar_scheme)
    tmpstr = uppercase(tmpstr)
    select case (tmpstr)
      case (CONSTANT_STRING)
        CS%BBL_mstar_scheme = Use_Fixed_mstar
      case (OM4_STRING)
        CS%BBL_mstar_scheme = mstar_from_Ekman
      case (RH18_STRING)
        CS%BBL_mstar_scheme = mstar_from_RH18
      case default
        call MOM_mesg('energetic_PBL_init: EPBL_BBL_MSTAR_SCHEME ="'//trim(tmpstr)//'"', 0)
        call MOM_error(FATAL, "energetic_PBL_init: Unrecognized setting "// &
              "EPBL_BBL_MSTAR_SCHEME = "//trim(tmpstr)//" found in input file.")
    end select
    call get_param(param_file, mdl, "MSTAR_BBL", CS%BBL_fixed_mstar, &
                   "The ratio of the friction velocity cubed to the TKE input to the "//&
                   "bottom boundary layer.  This option is used if EPBL_BBL_MSTAR_SCHEME = CONSTANT.", &
                   units="nondim", default=1.2, do_not_log=(CS%BBL_mstar_scheme/=Use_Fixed_mstar))
  endif

!/ Convective turbulence related options
  call get_param(param_file, mdl, "NSTAR", CS%nstar, &
                 "The portion of the buoyant potential energy imparted by "//&
                 "surface fluxes that is available to drive entrainment "//&
                 "at the base of mixed layer when that energy is positive.", &
                 units="nondim", default=0.2)
  call get_param(param_file, mdl, "MSTAR_CONV_ADJ", CS%mstar_convect_coef, &
                 "Coefficient used for reducing mstar during convection "//&
                 "due to reduction of stable density gradient.", &
                 units="nondim", default=0.0)

!/ Mixing Length Options
  call get_param(param_file, mdl, "USE_MLD_ITERATION", CS%Use_MLD_iteration, &
                 "A logical that specifies whether or not to use the "//&
                 "distance to the bottom of the actively turbulent boundary "//&
                 "layer to help set the EPBL length scale.", default=.true.)
  call get_param(param_file, mdl, "EPBL_TRANSITION_SCALE", CS%transLay_scale, &
                 "A scale for the mixing length in the transition layer "//&
                 "at the edge of the boundary layer as a fraction of the "//&
                 "boundary layer thickness.", units="nondim", default=0.1)
  if ( CS%Use_MLD_iteration .and. abs(CS%transLay_scale-0.5) >= 0.5) then
    call MOM_error(FATAL, "If flag USE_MLD_ITERATION is true, then "//&
                 "EPBL_TRANSITION should be greater than 0 and less than 1.")
  endif

  call get_param(param_file, mdl, "MLD_ITERATION_GUESS", CS%MLD_ITERATION_GUESS, &
                 "If true, use the previous timestep MLD as a first guess in the MLD iteration, "//&
                 "otherwise use half the ocean depth as the first guess of the boundary layer "//&
                 "depth.  The default is false to facilitate reproducibility.", &
                 default=.false., do_not_log=.not.CS%Use_MLD_iteration)
  call get_param(param_file, mdl, "EPBL_MLD_TOLERANCE", CS%MLD_tol, &
                 "The tolerance for the iteratively determined mixed "//&
                 "layer depth.  This is only used with USE_MLD_ITERATION.", &
                 units="meter", default=1.0, scale=US%m_to_Z, do_not_log=.not.CS%Use_MLD_iteration)
  call get_param(param_file, mdl, "EPBL_MLD_BISECTION", CS%MLD_bisection, &
                 "If true, use bisection with the iterative determination of the self-consistent "//&
                 "mixed layer depth.  Otherwise use the false position after a maximum and minimum "//&
                 "bound have been evaluated and the returned value or bisection before this.", &
                 default=.false., do_not_log=.not.CS%Use_MLD_iteration)
   call get_param(param_file, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
   call get_param(param_file, mdl, "EPBL_MLD_ITER_BUG", CS%MLD_iter_bug, &
                 "If true, use buggy logic that gives the wrong bounds for the next iteration "//&
                 "when successive guesses increase by exactly EPBL_MLD_TOLERANCE.", &
                 default=enable_bugs, do_not_log=.not.CS%Use_MLD_iteration)
  call get_param(param_file, mdl, "EPBL_MLD_MAX_ITS", CS%max_MLD_its, &
                 "The maximum number of iterations that can be used to find a self-consistent "//&
                 "mixed layer depth.  If EPBL_MLD_BISECTION is true, the maximum number "//&
                 "of iterations needed is set by Depth/2^MAX_ITS < EPBL_MLD_TOLERANCE.", &
                 default=20, do_not_log=.not.CS%Use_MLD_iteration)
  if (.not.CS%Use_MLD_iteration) CS%Max_MLD_Its = 1
  call get_param(param_file, mdl, "EPBL_MIN_MIX_LEN", CS%min_mix_len, &
                 "The minimum mixing length scale that will be used "//&
                 "by ePBL.  The default (0) does not set a minimum.", &
                 units="meter", default=0.0, scale=US%m_to_Z)

  call get_param(param_file, mdl, "MIX_LEN_EXPONENT", CS%MixLenExponent, &
                 "The exponent applied to the ratio of the distance to the MLD "//&
                 "and the MLD depth which determines the shape of the mixing length. "//&
                 "This is only used if USE_MLD_ITERATION is True.", &
                 units="nondim", default=2.0)

!/ Turbulent velocity scale in mixing coefficient
  call get_param(param_file, mdl, "EPBL_VEL_SCALE_SCHEME", vel_scale_str, &
                 "Selects the method for translating TKE into turbulent velocities. "//&
                 "Valid values are: \n"//&
                 "\t CUBE_ROOT_TKE  - A constant times the cube root of remaining TKE. \n"//&
                 "\t REICHL_H18 - Use the scheme based on a combination of w* and v* as \n"//&
                 "\t              documented in Reichl & Hallberg, 2018.", &
                 default=ROOT_TKE_STRING, do_not_log=.true.)
  call get_param(param_file, mdl, "EPBL_VEL_SCALE_MODE", wT_mode, default=-1)
  if (wT_mode == 0) then
    vel_scale_str = ROOT_TKE_STRING
    call MOM_error(WARNING, "Use EPBL_VEL_SCALE_SCHEME = CUBE_ROOT_TKE instead of the archaic EPBL_VEL_SCALE_MODE = 0.")
  elseif (wT_mode == 1) then
    vel_scale_str = RH18_STRING
    call MOM_error(WARNING, "Use EPBL_VEL_SCALE_SCHEME = REICHL_H18 instead of the archaic EPBL_VEL_SCALE_MODE = 1.")
  elseif (wT_mode >= 2) then
    call MOM_error(FATAL, "An unrecognized value of the obsolete parameter EPBL_VEL_SCALE_MODE was specified.")
  endif
  call log_param(param_file, mdl, "EPBL_VEL_SCALE_SCHEME", vel_scale_str, &
                 "Selects the method for translating TKE into turbulent velocities. "//&
                 "Valid values are: \n"//&
                 "\t CUBE_ROOT_TKE  - A constant times the cube root of remaining TKE. \n"//&
                 "\t REICHL_H18 - Use the scheme based on a combination of w* and v* as \n"//&
                 "\t              documented in Reichl & Hallberg, 2018.", &
                 default=ROOT_TKE_STRING)
  vel_scale_str = uppercase(vel_scale_str)
  select case (vel_scale_str)
    case (ROOT_TKE_STRING)
      CS%wT_scheme = wT_from_cRoot_TKE
    case (RH18_STRING)
      CS%wT_scheme = wT_from_RH18
    case default
      call MOM_mesg('energetic_PBL_init: EPBL_VEL_SCALE_SCHEME ="'//trim(vel_scale_str)//'"', 0)
      call MOM_error(FATAL, "energetic_PBL_init: Unrecognized setting "// &
            "EPBL_VEL_SCALE_SCHEME = "//trim(vel_scale_str)//" found in input file.")
  end select

  call get_param(param_file, mdl, "WSTAR_USTAR_COEF", CS%wstar_ustar_coef, &
                 "A ratio relating the efficiency with which convectively "//&
                 "released energy is converted to a turbulent velocity, "//&
                 "relative to mechanically forced TKE. Making this larger "//&
                 "increases the BL diffusivity", units="nondim", default=1.0)
  call get_param(param_file, mdl, "EPBL_VEL_SCALE_FACTOR", CS%vstar_scale_fac, &
                 "An overall nondimensional scaling factor for wT. "//&
                 "Making this larger increases the PBL diffusivity.", &
                 units="nondim", default=1.0)
  call get_param(param_file, mdl, "VSTAR_SURF_FAC", CS%vstar_surf_fac,&
                 "The proportionality times ustar to set vstar at the surface.", &
                 units="nondim", default=1.2)

  !/ Bottom boundary layer mixing related options
  call get_param(param_file, mdl, "EPBL_BBL_EFFIC", CS%ePBL_BBL_effic, &
                 "The efficiency of bottom boundary layer mixing via ePBL.  Setting this to a "//&
                 "value that is greater than 0 to enable bottom boundary layer mixing from EPBL.", &
                 units="nondim", default=0.0, scale=US%L_to_Z**2)
  call get_param(param_file, mdl, "EPBL_BBL_TIDAL_EFFIC", CS%ePBL_tidal_effic, &
                 "The efficiency of bottom boundary layer mixing via ePBL driven by the "//&
                 "bottom drag dissipation of tides, as provided in fluxes%BBL_tidal_dis.", &
                 units="nondim", default=0.0, scale=US%L_to_Z**2) !### Change the default to follow EPBL_BBL_EFFIC?
  no_BBL = ((CS%ePBL_BBL_effic <= 0.0) .and. (CS%ePBL_tidal_effic <= 0.0))

  call get_param(param_file, mdl, "USE_BBLD_ITERATION", CS%Use_BBLD_iteration, &
                 "A logical that specifies whether or not to use the distance to the top of the "//&
                 "actively turbulent bottom boundary layer to help set the EPBL length scale.", &
                 default=.true., do_not_log=no_BBL)
  call get_param(param_file, mdl, "TKE_DECAY_BBL", CS%TKE_decay_BBL, &
                 "TKE_DECAY_BBL relates the vertical rate of decay of the TKE available for "//&
                 "mechanical entrainment in the bottom boundary layer to the natural Ekman depth.", &
                 units="nondim", default=CS%TKE_decay, do_not_log=no_BBL)
  call get_param(param_file, mdl, "MIX_LEN_EXPONENT_BBL", CS%MixLenExponent_BBL, &
                 "The exponent applied to the ratio of the distance to the top of the BBL "//&
                 "and the total BBL depth which determines the shape of the mixing length. "//&
                 "This is only used if USE_MLD_ITERATION is True.", &
                 units="nondim", default=2.0, do_not_log=(no_BBL.or.(.not.CS%Use_BBLD_iteration)))
  call get_param(param_file, mdl, "EPBL_MIN_BBL_MIX_LEN", CS%min_BBL_mix_len, &
                 "The minimum mixing length scale that will be used by ePBL for bottom boundary "//&
                 "layer mixing.  Choosing (0) does not set a minimum.", &
                 units="meter", default=CS%min_mix_len, scale=US%m_to_Z, do_not_log=no_BBL)
  call get_param(param_file, mdl, "EPBL_BBLD_TOLERANCE", CS%BBLD_tol, &
                 "The tolerance for the iteratively determined bottom boundary layer depth.  "//&
                 "This is only used with USE_MLD_ITERATION.", &
                 units="meter", default=US%Z_to_m*CS%MLD_tol, scale=US%m_to_Z, &
                 do_not_log=(no_BBL.or.(.not.CS%Use_MLD_iteration)))
  call get_param(param_file, mdl, "EPBL_BBLD_MAX_ITS", CS%max_BBLD_its, &
                 "The maximum number of iterations that can be used to find a self-consistent "//&
                 "bottom boundary layer depth.", &
                 default=CS%max_MLD_its, do_not_log=(no_BBL.or.(.not.CS%Use_MLD_iteration)))
  if (.not.CS%Use_MLD_iteration) CS%max_BBLD_its = 1

  call get_param(param_file, mdl, "EPBL_BBL_VEL_SCALE_SCHEME", tmpstr, &
                 "Selects the method for translating bottom boundary layer TKE into turbulent velocities. "//&
                 "Valid values are: \n"//&
                 "\t CUBE_ROOT_TKE  - A constant times the cube root of remaining BBL TKE. \n"//&
                 "\t REICHL_H18 - Use the scheme based on a combination of w* and v* as \n"//&
                 "\t              documented in Reichl & Hallberg, 2018.", &
                 default=vel_scale_str, do_not_log=no_BBL)
  select case (tmpstr)
    case (ROOT_TKE_STRING)
      CS%wT_scheme_BBL = wT_from_cRoot_TKE
    case (RH18_STRING)
      CS%wT_scheme_BBL = wT_from_RH18
    case default
      call MOM_mesg('energetic_PBL_init: EPBL_BBL_VEL_SCALE_SCHEME ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "energetic_PBL_init: Unrecognized setting "// &
            "EPBL_BBL_VEL_SCALE_SCHEME = "//trim(tmpstr)//" found in input file.")
  end select
  call get_param(param_file, mdl, "EPBL_BBL_VEL_SCALE_FACTOR", CS%vstar_scale_fac_BBL, &
                 "An overall nondimensional scaling factor for wT in the bottom boundary layer. "//&
                 "Making this larger increases the bottom boundary layer diffusivity.", &
                 units="nondim", default=CS%vstar_scale_fac, do_not_log=no_BBL)
  call get_param(param_file, mdl, "VSTAR_BBL_SURF_FAC", CS%vstar_surf_fac_BBL,&
                 "The proportionality times ustar to set vstar in the bottom boundary layer.", &
                 units="nondim", default=CS%vstar_surf_fac, do_not_log=(no_BBL.or.(CS%wT_scheme_BBL/=wT_from_RH18)))
  call get_param(param_file, mdl, "EKMAN_SCALE_COEF_BBL", CS%Ekman_scale_coef_BBL, &
                 "A nondimensional scaling factor controlling the inhibition of the diffusive "//&
                 "length scale by rotation in the bottom boundary layer.  Making this larger "//&
                 "decreases the bottom boundary layer diffusivity.", &
                 units="nondim", default=CS%Ekman_scale_coef, do_not_log=no_BBL)
  call get_param(param_file, mdl, "EPBL_BBL_EFFIC_BUG", CS%BBL_effic_bug, &
                 "If true, overestimate the efficiency of the non-tidal ePBL bottom boundary "//&
                 "layer diffusivity by a factor of 1/sqrt(CDRAG), which is often a factor of "//&
                 "about 18.3.", default=.false., do_not_log=(CS%ePBL_BBL_effic<=0.0))

  call get_param(param_file, mdl, "DECAY_ADJUSTED_BBL_TKE", CS%decay_adjusted_BBL_TKE, &
                 "If true, include an adjustment factor in the bottom boundary layer energetics "//&
                 "that accounts for an exponential decay of TKE from a near-bottom source and "//&
                 "an assumed piecewise linear profile of the buoyancy flux response to a change "//&
                 "in a diffusivity.", &
                 default=.false., do_not_log=no_BBL)

  !/ Options related to Langmuir turbulence
  call get_param(param_file, mdl, "USE_LA_LI2016", use_LA_Windsea, &
       "A logical to use the Li et al. 2016 (submitted) formula to "//&
       "determine the Langmuir number.", default=.false.)
  ! Note this can be activated in other ways, but this preserves the old method.
  if (use_LA_windsea) then
    CS%use_LT = .true.
  else
    call get_param(param_file, mdl, "EPBL_LT", CS%use_LT, &
                 "A logical to use a LT parameterization.", default=.false.)
  endif
  if (CS%use_LT) then
    call get_param(param_file, mdl, "EPBL_LANGMUIR_SCHEME", tmpstr, &
                 "EPBL_LANGMUIR_SCHEME selects the method for including Langmuir turbulence. "//&
                 "Valid values are: \n"//&
                 "\t NONE     - Do not do any extra mixing due to Langmuir turbulence \n"//&
                 "\t RESCALE  - Use a multiplicative rescaling of mstar to account for Langmuir turbulence \n"//&
                 "\t ADDITIVE - Add a Langmuir turbulence contribution to mstar to other contributions", &
                 default=NONE_STRING, do_not_log=.true.)
    call get_param(param_file, mdl, "LT_ENHANCE", LT_enhance, default=-1)
    if (LT_enhance == 0) then
      tmpstr = NONE_STRING
      call MOM_error(WARNING, "Use EPBL_LANGMUIR_SCHEME = NONE instead of the archaic LT_ENHANCE = 0.")
    elseif (LT_enhance == 1) then
      call MOM_error(FATAL, "You are using a legacy LT_ENHANCE mode in ePBL that has been phased out. "//&
                            "If you need to use this setting please report this error.  Also use "//&
                            "EPBL_LANGMUIR_SCHEME to specify the scheme for mstar.")
    elseif (LT_enhance == 2) then
      tmpstr = RESCALED_STRING
      call MOM_error(WARNING, "Use EPBL_LANGMUIR_SCHEME = RESCALE instead of the archaic LT_ENHANCE = 2.")
    elseif (LT_enhance == 3) then
      tmpstr = ADDITIVE_STRING
      call MOM_error(WARNING, "Use EPBL_LANGMUIR_SCHEME = ADDITIVE instead of the archaic LT_ENHANCE = 3.")
    elseif (LT_enhance > 3) then
      call MOM_error(FATAL, "An unrecognized value of the obsolete parameter LT_ENHANCE was specified.")
    endif
    call log_param(param_file, mdl, "EPBL_LANGMUIR_SCHEME", tmpstr, &
                 "EPBL_LANGMUIR_SCHEME selects the method for including Langmuir turbulence. "//&
                 "Valid values are: \n"//&
                 "\t NONE     - Do not do any extra mixing due to Langmuir turbulence \n"//&
                 "\t RESCALE  - Use a multiplicative rescaling of mstar to account for Langmuir turbulence \n"//&
                 "\t ADDITIVE - Add a Langmuir turbulence contribution to mstar to other contributions", &
                 default=NONE_STRING)
    tmpstr = uppercase(tmpstr)
    select case (tmpstr)
      case (NONE_STRING)
        CS%LT_enhance_form = No_Langmuir
      case (RESCALED_STRING)
        CS%LT_enhance_form = Langmuir_rescale
      case (ADDITIVE_STRING)
        CS%LT_enhance_form = Langmuir_add
      case default
        call MOM_mesg('energetic_PBL_init: EPBL_LANGMUIR_SCHEME ="'//trim(tmpstr)//'"', 0)
        call MOM_error(FATAL, "energetic_PBL_init: Unrecognized setting "// &
              "EPBL_LANGMUIR_SCHEME = "//trim(tmpstr)//" found in input file.")
    end select

    call get_param(param_file, mdl, "LT_ENHANCE_COEF", CS%LT_enhance_coef, &
                 "Coefficient for Langmuir enhancement of mstar", &
                 units="nondim", default=0.447, do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_ENHANCE_EXP", CS%LT_enhance_exp, &
                 "Exponent for Langmuir enhancement of mstar", &
                 units="nondim", default=-1.33,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC1", CS%LaC_MLD_Ek, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "MLD approaching Ekman depth.", &
                 units="nondim", default=-0.87,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC2", CS%LaC_MLD_Ob_stab, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "MLD approaching stable Obukhov depth.", &
                 units="nondim", default=0.0,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC3", CS%LaC_MLD_Ob_un, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "MLD approaching unstable Obukhov depth.", &
                 units="nondim", default=0.0,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC4", CS%Lac_Ek_Ob_stab, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "ratio of Ekman to stable Obukhov depth.", &
                 units="nondim", default=0.95,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC5", CS%Lac_Ek_Ob_un, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "ratio of Ekman to unstable Obukhov depth.", &
                 units="nondim", default=0.95,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
  endif

  !/Options related to Machine Learning Equation Discovery
  ! Logial flags for using shape function from equation discovery - machine learning
  ! EPBL_EQD_DIFFUSIVITY : EPBL + Equation Discovery Diffusivity parameters

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_SHAPE", CS%eqdisc, &
                 "Logical flag for activating ML equation for shape function "// &
                 "that uses forcing to change its structure.", &
                 units="nondim", default=.false.)

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_VELOCITY", CS%eqdisc_v0, &
                   "Logical flag for activating ML equation discovery for velocity scale", &
                   units="nondim", default=.false.)

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_VELOCITY_H", CS%eqdisc_v0h, &
                   "Logical flag for activating ML equation discovery for velocity scale with h as input", &
                   units="nondim", default=.false.)


  ! sets a  lower cap for abs_f (Coriolis parameter) required in equation for v_0.
  ! Small value, solution not sensitive below 1 deg Latitute
  ! Default value of 2.5384E-07 corresponds to 0.1 deg.
  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_CORIOLIS_LOWER_CAP", CS%f_lower, &
                       "value of lower limit cap for v0, default is for 0.1 deg, insensitive below 1deg", &
                       units="s-1", default=2.5384E-07, scale=US%T_to_S, &
                       do_not_log=.not.CS%eqdisc_v0)

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_V0_LOWER_CAP", CS%v0_lower_cap, &
                       "value of lower limit cap for Coriolis in v0", &
                       units="m s-1", default=0.0001, scale=US%m_to_Z*US%T_to_s, &
                       do_not_log=.not.(CS%eqdisc_v0.or.CS%eqdisc_v0h))

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_V0_UPPER_CAP", CS%v0_upper_cap, &
                       "value of upper limit cap for Coriolis in v0", &
                       units="m s-1", default=0.1, scale=US%m_to_Z*US%T_to_s, &
                       do_not_log=.not.(CS%eqdisc_v0.or.CS%eqdisc_v0h))

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_BFLUX_LOWER_CAP", CS%bflux_lower_cap, &
                       "value of lower limit cap for Bflux used in setting in v0", &
                       units="m2 s-3", default=-7.0E-07, scale=(US%m_to_L**2)*(US%T_to_s**3), &
                       do_not_log=.not.(CS%eqdisc_v0.or.CS%eqdisc_v0h))

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_BFLUX_UPPER_CAP", CS%bflux_upper_cap, &
                       "value of upper limit cap for Bflux used in setting in v0", &
                       units="m2 s-3", default=7.0E-07, scale=(US%m_to_L**2)*(US%T_to_s**3), &
                       do_not_log=.not.(CS%eqdisc_v0.or.CS%eqdisc_v0h))

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_SIGMA_MAX_LOWER_CAP", CS%sigma_max_lower_cap, &
                       "value of lower limit cap for sigma coordinate of maximum for diffusivity", &
                       units="nondim", default=0.1, do_not_log=.not.CS%eqdisc)

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_SIGMA_MAX_UPPER_CAP", CS%sigma_max_upper_cap, &
                       "value of upper limit cap for sigma coordinate of maximum for diffusivity", &
                       units="nondim", default=0.7, do_not_log=.not.CS%eqdisc)

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_EH_UPPER_CAP", CS%Eh_upper_cap, &
                       "value of upper limit cap for boundary layer depth by Ekman depth hf/u", &
                       units="nondim", default=2.0, do_not_log=.not.CS%eqdisc)

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_LH_CAP", CS%Lh_cap, &
                       "value of upper limit cap for boundary layer depth by Monin-Obukhov depth hB/u^3", &
                       units="nondim", default=8.0, do_not_log=.not.CS%eqdisc)

  ! The coefficients used for machine learned diffusivity
  ! c1 to c6 used for sigma_m,
  !  7 to 9 v_0 surface heating, 10 to 14 v_0 surface cooling (ML velocity scale without h as input)
  ! 14, 15, & 16 for v_0h surface heating, 17, 18, & 14 for v_0h surface cooling (ML velocity scale with h as input)
  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_COEFFS", CS%ML_c, &
                 "Coefficient used for ML diffusivity 1 to 18 ", units="nondim", &
                  defaults=(/1.7908 , 0.6904, 0.0712, 0.4380, 2.6821, 1.5845, 0.1550,  1.1120,  0.8616, 0.0984, &
                             45.0,    2.8570, 3.290,  0.0785, 0.650,  0.0944, 6.0277, 15.7292 /), &
                  do_not_log=.not.(CS%eqdisc .or. CS%eqdisc_v0 .or. CS%eqdisc_v0h))

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_SHAPE_FUNCTION_EPSILON", CS%shape_function_epsilon, &
                 "Constant value of OSBL shape function below the boundary layer", &
                 units="nondim", default=0.01, do_not_log=.not.CS%eqdisc)

  !/ options end for Machine Learning Equation Discovery

  !/ Options for documenting differences from parameter choices
  call get_param(param_file, mdl, "EPBL_OPTIONS_DIFF", CS%options_diff, &
                 "If positive, this is a coded integer indicating a pair of settings whose "//&
                 "differences are diagnosed in a passive diagnostic mode  via extra calls to "//&
                 "ePBL_column.  If this is 0 or negative no extra calls occur.", &
                 default=0)
  if (CS%options_diff > 0) then
    if (CS%options_diff == 1) then
      diff_text = "EPBL_ORIGINAL_PE_CALC settings"
    elseif (CS%options_diff == 2) then
      diff_text = "EPBL_ANSWER_DATE settings"
    elseif (CS%options_diff == 3) then
      diff_text = "DIRECT_EPBL_MIXING_CALC settings"
    elseif (CS%options_diff == 4) then
      diff_text = "BBL DIRECT_EPBL_MIXING_CALC settings"
    elseif (CS%options_diff == 5) then
      diff_text = "BBL DECAY_ADJUSTED_BBL_TKE settings"
    else
      diff_text = "unchanged settings"
    endif
  endif

!/ Logging parameters
  ! This gives a minimum decay scale that is typically much less than Angstrom.
  CS%ustar_min = 2e-4*CS%omega*(GV%Angstrom_Z + GV%dZ_subroundoff)
  call log_param(param_file, mdl, "!EPBL_USTAR_MIN", CS%ustar_min, &
                 "The (tiny) minimum friction velocity used within the "//&
                 "ePBL code, derived from OMEGA and ANGSTROM.", &
                 units="m s-1", unscale=US%Z_to_m*US%s_to_T, &
                 like_default=.true.)


!/ Checking output flags
  CS%id_ML_depth = register_diag_field('ocean_model', 'ePBL_h_ML', diag%axesT1, &
      Time, 'Surface boundary layer depth', units='m', conversion=US%Z_to_m, &
      cmor_long_name='Ocean Mixed Layer Thickness Defined by Mixing Scheme')
  ! This is an alias for the same variable as ePBL_h_ML
  CS%id_hML_depth = register_diag_field('ocean_model', 'h_ML', diag%axesT1, &
      Time, 'Surface mixed layer depth based on active turbulence', units='m', conversion=US%Z_to_m)
  CS%id_ustar_ePBL = register_diag_field('ocean_model', 'ePBL_ustar', diag%axesT1, &
      Time, 'Surface friction in ePBL', units='m s-1', conversion=US%Z_to_m*US%s_to_T)
  CS%id_bflx_ePBL = register_diag_field('ocean_model', 'ePBL_bflx', diag%axesT1, &
      Time, 'Surface buoyancy flux in ePBL', units='m2 s-3', conversion=US%Z_to_m**2*US%s_to_T**3)
  CS%id_TKE_wind = register_diag_field('ocean_model', 'ePBL_TKE_wind', diag%axesT1, &
      Time, 'Wind-stirring source of mixed layer TKE', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_TKE_MKE = register_diag_field('ocean_model', 'ePBL_TKE_MKE', diag%axesT1, &
      Time, 'Mean kinetic energy source of mixed layer TKE', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_TKE_conv = register_diag_field('ocean_model', 'ePBL_TKE_conv', diag%axesT1, &
      Time, 'Convective source of mixed layer TKE', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_TKE_forcing = register_diag_field('ocean_model', 'ePBL_TKE_forcing', diag%axesT1, &
      Time, 'TKE consumed by mixing surface forcing or penetrative shortwave radation '//&
            'through model layers', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_TKE_mixing = register_diag_field('ocean_model', 'ePBL_TKE_mixing', diag%axesT1, &
      Time, 'TKE consumed by mixing that deepens the mixed layer', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_TKE_mech_decay = register_diag_field('ocean_model', 'ePBL_TKE_mech_decay', diag%axesT1, &
      Time, 'Mechanical energy decay sink of mixed layer TKE', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_TKE_conv_decay = register_diag_field('ocean_model', 'ePBL_TKE_conv_decay', diag%axesT1, &
      Time, 'Convective energy decay sink of mixed layer TKE', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_Mixing_Length = register_diag_field('ocean_model', 'Mixing_Length', diag%axesTi, &
      Time, 'Mixing Length that is used', units='m', conversion=US%Z_to_m)
  CS%id_Velocity_Scale = register_diag_field('ocean_model', 'Velocity_Scale', diag%axesTi, &
      Time, 'Velocity Scale that is used.', units='m s-1', conversion=US%Z_to_m*US%s_to_T)
  CS%id_mstar_sfc = register_diag_field('ocean_model', 'MSTAR', diag%axesT1, &
      Time, 'Total mstar that is used.', 'nondim')
  if ((CS%ePBL_BBL_effic > 0.0) .or. (CS%ePBL_tidal_effic > 0.0) .or. CS%ePBL_BBL_use_mstar) then
    CS%id_Kd_BBL = register_diag_field('ocean_model', 'Kd_ePBL_BBL', diag%axesTi, &
        Time, 'ePBL bottom boundary layer diffusivity', units='m2 s-1', conversion=GV%HZ_T_to_m2_s)
    CS%id_BBL_Mix_Length = register_diag_field('ocean_model', 'BBL_Mixing_Length', diag%axesTi, &
        Time, 'ePBL bottom boundary layer mixing length', units='m', conversion=US%Z_to_m)
    CS%id_BBL_Vel_Scale = register_diag_field('ocean_model', 'BBL_Velocity_Scale', diag%axesTi, &
        Time, 'ePBL bottom boundary layer velocity scale', units='m s-1', conversion=US%Z_to_m*US%s_to_T)
    CS%id_BBL_depth = register_diag_field('ocean_model', 'h_BBL', diag%axesT1, &
        Time, 'Bottom boundary layer depth based on active turbulence', units='m', conversion=US%Z_to_m)
    CS%id_ustar_BBL = register_diag_field('ocean_model', 'ePBL_ustar_BBL', diag%axesT1, &
        Time, 'The bottom boundary layer friction velocity', units='m s-1', conversion=GV%H_to_m*US%s_to_T)
    CS%id_BBL_decay_scale = register_diag_field('ocean_model', 'BBL_decay_scale', diag%axesT1, &
        Time, 'The bottom boundary layer TKE decay lengthscale', units='m', conversion=GV%H_to_m)
    CS%id_TKE_BBL = register_diag_field('ocean_model', 'ePBL_BBL_TKE', diag%axesT1, &
        Time, 'The source of TKE for the bottom boundary layer', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
    CS%id_TKE_BBL_mixing = register_diag_field('ocean_model', 'ePBL_BBL_TKE_mixing', diag%axesT1, &
        Time, 'TKE consumed by mixing that thickens the bottom boundary layer', &
        units='W m-2', conversion=US%RZ3_T3_to_W_m2)
    CS%id_TKE_BBL_decay = register_diag_field('ocean_model', 'ePBL_BBL_TKE_decay', diag%axesT1, &
        Time, 'Energy decay sink of mixed layer TKE in the bottom boundary layer', &
        units='W m-2', conversion=US%RZ3_T3_to_W_m2)
    CS%id_mstar_BBL = register_diag_field('ocean_model', 'MSTAR_BBL', diag%axesT1, &
        Time, 'Total BBL mstar that is used.', 'nondim')
  endif
  if (CS%use_LT) then
    CS%id_LA = register_diag_field('ocean_model', 'LA', diag%axesT1, &
        Time, 'Langmuir number.', 'nondim')
    CS%id_LA_mod = register_diag_field('ocean_model', 'LA_MOD', diag%axesT1, &
        Time, 'Modified Langmuir number.', 'nondim')
    CS%id_mstar_LT = register_diag_field('ocean_model', 'MSTAR_LT', diag%axesT1, &
        Time, 'Increase in mstar due to Langmuir Turbulence.', 'nondim')
  endif

  if (CS%options_diff > 0) then
    CS%id_opt_diff_Kd_ePBL = register_diag_field('ocean_model', 'ePBL_opt_diff_Kd_ePBL', diag%axesTi, &
        Time, 'Change in ePBL diapycnal diffusivity at interfaces due to '//trim(diff_text), &
        units='m2 s-1', conversion=GV%HZ_T_to_m2_s)
    CS%id_opt_maxdiff_Kd_ePBL = register_diag_field('ocean_model', 'ePBL_opt_maxdiff_Kd_ePBL', diag%axesT1, &
        Time, 'Column maximum change in ePBL diapycnal diffusivity at interfaces due to '//trim(diff_text), &
        units='m2 s-1', conversion=GV%HZ_T_to_m2_s)
    CS%id_opt_diff_hML_depth = register_diag_field('ocean_model', 'ePBL_opt_diff_h_ML', diag%axesT1, Time, &
        'Change in surface or bottom boundary layer depth based on active turbulence due to '//trim(diff_text), &
        units='m', conversion=US%Z_to_m)
  endif

  if (report_avg_its) then
    CS%sum_its(1) = real_to_EFP(0.0) ; CS%sum_its(2) = real_to_EFP(0.0)
    CS%sum_its_BBL(1) = real_to_EFP(0.0) ; CS%sum_its_BBL(2) = real_to_EFP(0.0)
  endif

  CS%TKE_diagnostics = (max(CS%id_TKE_wind, CS%id_TKE_MKE, CS%id_TKE_conv, &
                            CS%id_TKE_mixing, CS%id_TKE_mech_decay, CS%id_TKE_forcing, &
                            CS%id_TKE_conv_decay) > 0)
  if ((CS%ePBL_BBL_effic > 0.0) .or. (CS%ePBL_tidal_effic > 0.0) .or. CS%ePBL_BBL_use_mstar) then
    CS%TKE_diagnostics = CS%TKE_diagnostics .or. &
        (max(CS%id_TKE_BBL, CS%id_TKE_BBL_mixing, CS%id_TKE_BBL_decay) > 0)
  endif

  call safe_alloc_alloc(CS%ML_depth, isd, ied, jsd, jed)
  call safe_alloc_alloc(CS%BBL_depth, isd, ied, jsd, jed)

end subroutine energetic_PBL_init

!> Clean up and deallocate memory associated with the energetic_PBL module.
subroutine energetic_PBL_end(CS)
  type(energetic_PBL_CS), intent(inout) :: CS !< Energetic_PBL control structure

  character(len=256) :: mesg
  real :: avg_its ! The averaged number of iterations used by ePBL [nondim]

  if (allocated(CS%ML_depth))            deallocate(CS%ML_depth)
  if (allocated(CS%BBL_depth))           deallocate(CS%BBL_depth)

  if (report_avg_its) then
    call EFP_sum_across_PEs(CS%sum_its, 2)
    avg_its = EFP_to_real(CS%sum_its(1)) / EFP_to_real(CS%sum_its(2))
    write (mesg,*) "Average ePBL iterations = ", avg_its
    call MOM_mesg(mesg)

    if ((CS%ePBL_BBL_effic > 0.0) .or. (CS%ePBL_tidal_effic > 0.0) .or. CS%ePBL_BBL_use_mstar) then
      call EFP_sum_across_PEs(CS%sum_its_BBL, 2)
      avg_its = EFP_to_real(CS%sum_its_BBL(1)) / EFP_to_real(CS%sum_its_BBL(2))
      write (mesg,*) "Average ePBL BBL iterations = ", avg_its
      call MOM_mesg(mesg)
    endif
  endif
end subroutine energetic_PBL_end

!> \namespace MOM_energetic_PBL
!!
!! By Robert Hallberg, 2015.
!!
!!   This file contains the subroutine (energetic_PBL) that uses an
!! integrated boundary layer energy budget (like a bulk- or refined-
!! bulk mixed layer scheme), but instead of homogenizing this model
!! calculates a finite diffusivity and viscosity, which in this
!! regard is conceptually similar to what is done with KPP or various
!! two-equation closures.  However, the scheme that is implemented
!! here has the big advantage that is entirely implicit, but is
!! simple enough that it requires only a single vertical pass to
!! determine the diffusivity. The development of bulk mixed layer
!! models stems from the work of various people, as described in the
!! review paper by \cite niiler1977. The work here draws in
!! with particular on the form for TKE decay proposed by
!! \cite oberhuber1993, with an extension to a refined bulk mixed
!! layer as described in Hallberg (\cite muller2003).  The physical
!! processes portrayed in this subroutine include convectively driven
!! mixing and mechanically driven mixing.  Unlike boundary-layer
!! mixing, stratified shear mixing is not a one-directional turbulent
!! process, and it is dealt with elsewhere in the MOM6 code within
!! the module MOM_kappa_shear.F90.  It is assumed that the heat,
!! mass, and salt fluxes have been applied elsewhere, but that their
!! implications for the integrated TKE budget have been captured in
!! an array that is provided as an argument to this subroutine.  This
!! is a full 3-d array due to the effects of penetrating shortwave
!! radiation.

end module MOM_energetic_PBL
