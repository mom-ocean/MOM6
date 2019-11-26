!> Energetically consistent planetary boundary layer parameterization
module MOM_energetic_PBL

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_domains,       only : create_group_pass, do_group_pass, group_pass_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_string_functions, only : uppercase
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_wave_interface, only: wave_parameters_CS, Get_Langmuir_Number

! use MOM_EOS, only : calculate_density, calculate_density_derivs

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

  !/ Constants
  real    :: VonKar = 0.41   !<   The von Karman coefficient.  This should be runtime, but because
                             !!  it is runtime in KPP and set to 0.4 it might change answers.
  real    :: omega           !<   The Earth's rotation rate [T-1 ~> s-1].
  real    :: omega_frac      !<   When setting the decay scale for turbulence, use this fraction of
                             !!  the absolute rotation rate blended with the local value of f, as
                             !!  sqrt((1-of)*f^2 + of*4*omega^2) [nondim].

  !/ Convection related terms
  real    :: nstar           !< The fraction of the TKE input to the mixed layer available to drive
                             !! entrainment [nondim]. This quantity is the vertically integrated
                             !! buoyancy production minus the vertically integrated dissipation of
                             !! TKE produced by buoyancy.

  !/ Mixing Length terms
  logical :: Use_MLD_iteration=.false. !< False to use old ePBL method.
  logical :: MLD_iteration_guess=.false. !< False to default to guessing half the
                             !! ocean depth for the iteration.
  integer :: max_MLD_its     !< The maximum number of iterations that can be used to find a
                             !! self-consistent mixed layer depth with Use_MLD_iteration.
  real    :: MixLenExponent  !< Exponent in the mixing length shape-function.
                             !! 1 is law-of-the-wall at top and bottom,
                             !! 2 is more KPP like.
  real    :: MKE_to_TKE_effic !< The efficiency with which mean kinetic energy released by
                             !!  mechanically forced entrainment of the mixed layer is converted to
                             !!  TKE [nondim].
  real    :: ustar_min       !< A minimum value of ustar to avoid numerical problems [Z T-1 ~> m s-1].
                             !! If the value is small enough, this should not affect the solution.
  real    :: Ekman_scale_coef !< A nondimensional scaling factor controlling the inhibition of the
                             !! diffusive length scale by rotation.  Making this larger decreases
                             !! the diffusivity in the planetary boundary layer.
  real    :: transLay_scale  !< A scale for the mixing length in the transition layer
                             !! at the edge of the boundary layer as a fraction of the
                             !! boundary layer thickness.  The default is 0, but a
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
  real    :: vstar_scale_fac !< An overall nondimensional scaling factor for vstar times a unit
                             !! conversion factor [Z s T-1 m-1 ~> nondim].  Making this larger increases
                             !! the diffusivity.

  !mstar related options
  integer :: mstar_scheme    !< An encoded integer to determine which formula is used to set mstar
  logical :: MSTAR_FLATCAP=.true. !< Set false to use asymptotic mstar cap.
  real    :: mstar_cap       !< Since MSTAR is restoring undissipated energy to mixing,
                             !! there must be a cap on how large it can be.  This
                             !! is definitely a function of latitude (Ekman limit),
                             !! but will be taken as constant for now.

  !/ vertical decay related options
  real    :: TKE_decay       !< The ratio of the natural Ekman depth to the TKE decay scale [nondim].

  !/ mstar_scheme == 0
  real    :: fixed_mstar     !< Mstar is the ratio of the friction velocity cubed to the TKE available to
                             !! drive entrainment, nondimensional. This quantity is the vertically
                             !! integrated shear production minus the vertically integrated
                             !! dissipation of TKE produced by shear.  This value is used if the option
                             !! for using a fixed mstar is used.

  !/ mstar_scheme == 2
  real :: C_EK = 0.17        !< MSTAR Coefficient in rotation limit for mstar_scheme=OM4
  real :: MSTAR_COEF = 0.3   !< MSTAR coefficient in rotation/stabilizing balance for mstar_scheme=OM4

  !/ mstar_scheme == 3
  real    :: RH18_mstar_cN1  !< MSTAR_N coefficient 1 (outter-most coefficient for fit).
                             !! Value of 0.275 in RH18.  Increasing this
                             !! coefficient increases mechanical mixing for all values of Hf/ust,
                             !! but is most effective at low values (weakly developed OSBLs).
  real    :: RH18_mstar_cN2  !< MSTAR_N coefficient 2 (coefficient outside of exponential decay).
                             !! Value of 8.0 in RH18.  Increasing this coefficient increases MSTAR
                             !! for all values of HF/ust, with a consistent affect across
                             !! a wide range of Hf/ust.
  real    :: RH18_mstar_cN3  !< MSTAR_N coefficient 3 (exponential decay coefficient). Value of
                             !! -5.0 in RH18.  Increasing this increases how quickly the value
                             !! of MSTAR decreases as Hf/ust increases.
  real    :: RH18_mstar_cS1  !< MSTAR_S coefficient for RH18 in stabilizing limit.
                             !! Value of 0.2 in RH18.
  real    :: RH18_mstar_cS2  !< MSTAR_S exponent for RH18 in stabilizing limit.
                             !! Value of 0.4 in RH18.

  !/ Coefficient for shear/convective turbulence interaction
  real :: mstar_convect_coef !< Factor to reduce mstar when statically unstable.

  !/ Langmuir turbulence related parameters
  logical :: Use_LT = .false. !< Flag for using LT in Energy calculation
  integer :: LT_ENHANCE_FORM !< Integer for Enhancement functional form (various options)
  real    :: LT_ENHANCE_COEF !< Coefficient in fit for Langmuir Enhancment
  real    :: LT_ENHANCE_EXP  !< Exponent in fit for Langmuir Enhancement
  real :: LaC_MLDoEK         !< Coefficient for Langmuir number modification based on the ratio of
                             !! the mixed layer depth over the Ekman depth.
  real :: LaC_MLDoOB_stab    !< Coefficient for Langmuir number modification based on the ratio of
                             !! the mixed layer depth over the Obukov depth with stablizing forcing.
  real :: LaC_EKoOB_stab     !< Coefficient for Langmuir number modification based on the ratio of
                             !! the Ekman depth over the Obukov depth with stablizing forcing.
  real :: LaC_MLDoOB_un      !< Coefficient for Langmuir number modification based on the ratio of
                             !! the mixed layer depth over the Obukov depth with destablizing forcing.
  real :: LaC_EKoOB_un       !< Coefficient for Langmuir number modification based on the ratio of
                             !! the Ekman depth over the Obukov depth with destablizing forcing.
  real :: Max_Enhance_M = 5. !< The maximum allowed LT enhancement to the mixing.

  !/ Others
  type(time_type), pointer :: Time=>NULL() !< A pointer to the ocean model's clock.

  logical :: TKE_diagnostics = .false. !< If true, diagnostics of the TKE budget are being calculated.
  logical :: answers_2018    !< If true, use the order of arithmetic and expressions that recover the
                             !! answers from the end of 2018.  Otherwise, use updated and more robust
                             !! forms of the same expressions.
  logical :: orig_PE_calc    !< If true, the ePBL code uses the original form of the
                             !! potential energy change code.  Otherwise, it uses a newer version
                             !! that can work with successive increments to the diffusivity in
                             !! upward or downward passes.
  type(diag_ctrl), pointer :: diag=>NULL() !< A structure that is used to regulate the
                             !! timing of diagnostic output.

  real, allocatable, dimension(:,:) :: &
    ML_depth            !< The mixed layer depth determined by active mixing in ePBL [Z ~> m].

  ! These are terms in the mixed layer TKE budget, all in [R Z3 T-3 ~> W m-2 = kg s-3].
  real, allocatable, dimension(:,:) :: &
    diag_TKE_wind, &   !< The wind source of TKE [R Z3 T-3 ~> W m-2].
    diag_TKE_MKE, &    !< The resolved KE source of TKE [R Z3 T-3 ~> W m-2].
    diag_TKE_conv, &   !< The convective source of TKE [R Z3 T-3 ~> W m-2].
    diag_TKE_forcing, & !< The TKE sink required to mix surface penetrating shortwave heating
                       !! [R Z3 T-3 ~> W m-2].
    diag_TKE_mech_decay, & !< The decay of mechanical TKE [R Z3 T-3 ~> W m-2].
    diag_TKE_conv_decay, & !< The decay of convective TKE [R Z3 T-3 ~> W m-2].
    diag_TKE_mixing, & !< The work done by TKE to deepen the mixed layer [R Z3 T-3 ~> W m-2].
    ! These additional diagnostics are also 2d.
    MSTAR_MIX, &       !< Mstar used in EPBL [nondim]
    MSTAR_LT, &        !< Mstar due to Langmuir turbulence [nondim]
    LA, &              !< Langmuir number [nondim]
    LA_MOD             !< Modified Langmuir number [nondim]

  real, allocatable, dimension(:,:,:) :: &
    Velocity_Scale, & !< The velocity scale used in getting Kd [Z T-1 ~> m s-1]
    Mixing_Length     !< The length scale used in getting Kd [Z ~> m]
  !>@{ Diagnostic IDs
  integer :: id_ML_depth = -1, id_TKE_wind = -1, id_TKE_mixing = -1
  integer :: id_TKE_MKE = -1, id_TKE_conv = -1, id_TKE_forcing = -1
  integer :: id_TKE_mech_decay = -1, id_TKE_conv_decay = -1
  integer :: id_Mixing_Length = -1, id_Velocity_Scale = -1
  integer :: id_MSTAR_mix = -1, id_LA_mod = -1, id_LA = -1, id_MSTAR_LT = -1
  !!@}
end type energetic_PBL_CS

!>@{ Enumeration values for mstar_Scheme
integer, parameter :: Use_Fixed_MStar = 0  !< The value of mstar_scheme to use a constant mstar
integer, parameter :: MStar_from_Ekman = 2 !< The value of mstar_scheme to base mstar on the ratio
                                           !! of the Ekman layer depth to the Obukhov depth
integer, parameter :: MStar_from_RH18 = 3  !< The value of mstar_scheme to base mstar of of RH18
integer, parameter :: No_Langmuir = 0      !< The value of LT_ENHANCE_FORM not use Langmuir turbolence.
integer, parameter :: Langmuir_rescale = 2 !< The value of LT_ENHANCE_FORM to use a multiplicative
                                           !! rescaling of mstar to account for Langmuir turbulence.
integer, parameter :: Langmuir_add = 3     !< The value of LT_ENHANCE_FORM to add a contribution to
                                           !! mstar from Langmuir turblence to other contributions.
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
!!@}

!> A type for conveniently passing around ePBL diagnostics for a column.
type, public :: ePBL_column_diags ; private
  !>@{ Local column copies of energy change diagnostics, all in [R Z3 T-3 ~> W m-2].
  real :: dTKE_conv, dTKE_forcing, dTKE_wind, dTKE_mixing
  real :: dTKE_MKE, dTKE_mech_decay, dTKE_conv_decay
  !!@}
  real :: LA        !< The value of the Langmuir number [nondim]
  real :: LAmod     !< The modified Langmuir number by convection [nondim]
  real :: mstar     !< The value of mstar used in ePBL [nondim]
  real :: mstar_LT  !< The portion of mstar due to Langmuir turbulence [nondim]
  real, allocatable, dimension(:) :: dT_expect !< Expected temperature changes [degC]
  real, allocatable, dimension(:) :: dS_expect !< Expected salinity changes [ppt]
end type ePBL_column_diags

contains

!>    This subroutine determines the diffusivities from the integrated energetics
!!  mixed layer model.  It assumes that heating, cooling and freshwater fluxes
!!  have already been applied.  All calculations are done implicitly, and there
!!  is no stability limit on the time step.
subroutine energetic_PBL(h_3d, u_3d, v_3d, tv, fluxes, dt, Kd_int, G, GV, US, CS, &
                         dSV_dT, dSV_dS, TKE_forced, buoy_flux, dt_diag, last_call, &
                         dT_expected, dS_expected, Waves )
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
                                                   !! [R-1 degC-1 ~> m3 kg-1 degC-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dSV_dS !< The partial derivative of in-situ specific
                                                   !! volume with salinity [R-1 ppt-1 ~> m3 kg-1 ppt-1].
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
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(out)   :: Kd_int !< The diagnosed diffusivities at interfaces
                                                   !! [Z2 s-1 ~> m2 s-1].
  type(energetic_PBL_CS),  pointer       :: CS     !< The control structure returned by a previous
                                                   !! call to mixedlayer_init.
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: buoy_flux !< The surface buoyancy flux [Z2 T-3 ~> m2 s-3].
  real,          optional, intent(in)    :: dt_diag   !< The diagnostic time step, which may be less
                                                   !! than dt if there are two calls to mixedlayer [T ~> s].
  logical,       optional, intent(in)    :: last_call !< If true, this is the last call to
                                                   !! mixedlayer in the current time step, so
                                                   !! diagnostics will be written. The default
                                                   !! is .true.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: dT_expected !< The values of temperature change that
                                                   !! should be expected when the returned
                                                   !! diffusivities are applied [degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: dS_expected !< The values of salinity change that
                                                   !! should be expected when the returned
                                                   !! diffusivities are applied [ppt].
  type(wave_parameters_CS), &
                 optional, pointer       :: Waves  !< Wave CS

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
!   To use the classic constant mstar mixied layers choose MSTAR_SCHEME=CONSTANT.
! The key parameters then include mstar, nstar, TKE_decay, and conv_decay.
! For the Oberhuber (1993) mixed layer,the values of these are:
!      mstar = 1.25,  nstar = 1, TKE_decay = 2.5, conv_decay = 0.5
! TKE_decay is 1/kappa in eq. 28 of Oberhuber (1993), while conv_decay is 1/mu.
! For a traditional Kraus-Turner mixed layer, the values are:
!      mstar = 1.25, nstar = 0.4, TKE_decay = 0.0, conv_decay = 0.0

  ! Local variables
  real, dimension(SZI_(G),SZK_(GV)) :: &
    h_2d, &         ! A 2-d slice of the layer thickness [H ~> m or kg m-2].
    T_2d, &         ! A 2-d slice of the layer temperatures [degC].
    S_2d, &         ! A 2-d slice of the layer salinities [ppt].
    TKE_forced_2d, & ! A 2-d slice of TKE_forced [R Z3 T-2 ~> J m-2].
    dSV_dT_2d, &    ! A 2-d slice of dSV_dT [R-1 degC-1 ~> m3 kg-1 degC-1].
    dSV_dS_2d, &    ! A 2-d slice of dSV_dS [R-1 ppt-1 ~> m3 kg-1 ppt-1].
    u_2d, &         ! A 2-d slice of the zonal velocity [L T-1 ~> m s-1].
    v_2d            ! A 2-d slice of the meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK_(GV)+1) :: &
    Kd_2d           ! A 2-d version of the diapycnal diffusivity [Z2 T-1 ~> m2 s-1].
  real, dimension(SZK_(GV)) :: &
    h, &            ! The layer thickness [H ~> m or kg m-2].
    T0, &           ! The initial layer temperatures [degC].
    S0, &           ! The initial layer salinities [ppt].
    dSV_dT_1d, &    ! The partial derivatives of specific volume with temperature [R-1 degC-1 ~> m3 kg-1 degC-1].
    dSV_dS_1d, &    ! The partial derivatives of specific volume with salinity [R-1 ppt-1 ~> m3 kg-1 ppt-1].
    TKE_forcing, &  ! Forcing of the TKE in the layer coming from TKE_forced [R Z3 T-2 ~> J m-2].
    u, &            ! The zonal velocity [L T-1 ~> m s-1].
    v               ! The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZK_(GV)+1) :: &
    Kd, &           ! The diapycnal diffusivity [Z2 T-1 ~> m2 s-1].
    mixvel, &       ! A turbulent mixing veloxity [Z T-1 ~> m s-1].
    mixlen          ! A turbulent mixing length [Z ~> m].
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].

  real :: absf      ! The absolute value of f [T-1 ~> s-1].
  real :: U_star    ! The surface friction velocity [Z T-1 ~> m s-1].
  real :: U_Star_Mean ! The surface friction without gustiness [Z T-1 ~> m s-1].
  real :: B_Flux    ! The surface buoyancy flux [Z2 T-3 ~> m2 s-3]
  real :: MLD_io    ! The mixed layer depth found by ePBL_column [Z ~> m].

! The following are only used for diagnostics.
  real :: dt__diag  ! A copy of dt_diag (if present) or dt [T ~> s].
  logical :: write_diags  ! If true, write out diagnostics with this step.
  logical :: reset_diags  ! If true, zero out the accumulated diagnostics.

  logical :: debug=.false.  ! Change this hard-coded value for debugging.
  type(ePBL_column_diags) :: eCD ! A container for passing around diagnostics.

  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "energetic_PBL: "//&
         "Module must be initialized before it is used.")

  if (.not. associated(tv%eqn_of_state)) call MOM_error(FATAL, &
      "energetic_PBL: Temperature, salinity and an equation of state "//&
      "must now be used.")
  if (.NOT. associated(fluxes%ustar)) call MOM_error(FATAL, &
      "energetic_PBL: No surface TKE fluxes (ustar) defined in mixedlayer!")
  debug = .false. ; if (present(dT_expected) .or. present(dS_expected)) debug = .true.

  if (debug) allocate(eCD%dT_expect(nz), eCD%dS_expect(nz))

  h_neglect = GV%H_subroundoff

  dt__diag = dt ; if (present(dt_diag)) dt__diag = dt_diag
  write_diags = .true. ; if (present(last_call)) write_diags = last_call


  ! Determine whether to zero out diagnostics before accumulation.
  reset_diags = .true.
  if (present(dt_diag) .and. write_diags .and. (dt__diag > dt)) &
    reset_diags = .false.  ! This is the second call to mixedlayer.

  if (reset_diags) then
    if (CS%TKE_diagnostics) then
!!OMP parallel do default(none) shared(is,ie,js,je,CS)
      do j=js,je ; do i=is,ie
        CS%diag_TKE_wind(i,j) = 0.0 ; CS%diag_TKE_MKE(i,j) = 0.0
        CS%diag_TKE_conv(i,j) = 0.0 ; CS%diag_TKE_forcing(i,j) = 0.0
        CS%diag_TKE_mixing(i,j) = 0.0 ; CS%diag_TKE_mech_decay(i,j) = 0.0
        CS%diag_TKE_conv_decay(i,j) = 0.0 !; CS%diag_TKE_unbalanced(i,j) = 0.0
      enddo ; enddo
    endif
  endif
  ! if (CS%id_Mixing_Length>0) CS%Mixing_Length(:,:,:) = 0.0
  ! if (CS%id_Velocity_Scale>0) CS%Velocity_Scale(:,:,:) = 0.0

!!OMP parallel do default(private) shared(js,je,nz,is,ie,h_3d,u_3d,v_3d,tv,dt, &
!!OMP                                  CS,G,GV,US,fluxes,debug, &
!!OMP                                  TKE_forced,dSV_dT,dSV_dS,Kd_int)
  do j=js,je
    ! Copy the thicknesses and other fields to 2-d arrays.
    do k=1,nz ; do i=is,ie
      h_2d(i,k) = h_3d(i,j,k) ; u_2d(i,k) = u_3d(i,j,k) ; v_2d(i,k) = v_3d(i,j,k)
      T_2d(i,k) = tv%T(i,j,k) ; S_2d(i,k) = tv%S(i,j,k)
      TKE_forced_2d(i,k) = TKE_forced(i,j,k)
      dSV_dT_2d(i,k) = dSV_dT(i,j,k) ; dSV_dS_2d(i,k) = dSV_dS(i,j,k)
    enddo ; enddo

    !   Determine the initial mech_TKE and conv_PErel, including the energy required
    ! to mix surface heating through the topmost cell, the energy released by mixing
    ! surface cooling & brine rejection down through the topmost cell, and
    ! homogenizing the shortwave heating within that cell.  This sets the energy
    ! and ustar and wstar available to drive mixing at the first interior
    ! interface.
    do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then

      ! Copy the thicknesses and other fields to 1-d arrays.
      do k=1,nz
        h(k) = h_2d(i,k) + GV%H_subroundoff ; u(k) = u_2d(i,k) ; v(k) = v_2d(i,k)
        T0(k) = T_2d(i,k) ; S0(k) = S_2d(i,k) ; TKE_forcing(k) =  TKE_forced_2d(i,k)
        dSV_dT_1d(k) = dSV_dT_2d(i,k) ; dSV_dS_1d(k) = dSV_dS_2d(i,k)
      enddo
      do K=1,nz+1 ; Kd(K) = 0.0 ; enddo

      ! Make local copies of surface forcing and process them.
      u_star = fluxes%ustar(i,j)
      u_star_Mean = fluxes%ustar_gustless(i,j)
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
      if (CS%MLD_iteration_guess .and. (CS%ML_Depth(i,j) > 0.0))  MLD_io = CS%ML_Depth(i,j)

      call ePBL_column(h, u, v, T0, S0, dSV_dT_1d, dSV_dS_1d, TKE_forcing, B_flux, absf, &
                       u_star, u_star_mean, dt, MLD_io, Kd, mixvel, mixlen, GV, &
                       US, CS, eCD, dt_diag=dt_diag, Waves=Waves, G=G, i=i, j=j)


      ! Copy the diffusivities to a 2-d array.
      do K=1,nz+1
        Kd_2d(i,K) = Kd(K)
      enddo
      CS%ML_depth(i,j) = MLD_io

      if (present(dT_expected)) then
        do k=1,nz ; dT_expected(i,j,k) = eCD%dT_expect(k) ; enddo
      endif
      if (present(dS_expected)) then
        do k=1,nz ; dS_expected(i,j,k) = eCD%dS_expect(k) ; enddo
      endif

      if (CS%TKE_diagnostics) then
        CS%diag_TKE_MKE(i,j) = CS%diag_TKE_MKE(i,j) + eCD%dTKE_MKE
        CS%diag_TKE_conv(i,j) = CS%diag_TKE_conv(i,j) + eCD%dTKE_conv
        CS%diag_TKE_forcing(i,j) = CS%diag_TKE_forcing(i,j) + eCD%dTKE_forcing
        CS%diag_TKE_wind(i,j) = CS%diag_TKE_wind(i,j) + eCD%dTKE_wind
        CS%diag_TKE_mixing(i,j) = CS%diag_TKE_mixing(i,j) + eCD%dTKE_mixing
        CS%diag_TKE_mech_decay(i,j) = CS%diag_TKE_mech_decay(i,j) + eCD%dTKE_mech_decay
        CS%diag_TKE_conv_decay(i,j) = CS%diag_TKE_conv_decay(i,j) + eCD%dTKE_conv_decay
       ! CS%diag_TKE_unbalanced(i,j) = CS%diag_TKE_unbalanced(i,j) + eCD%dTKE_unbalanced
      endif
      ! Write to 3-D for outputing Mixing length and velocity scale.
      if (CS%id_Mixing_Length>0) then ; do k=1,nz
        CS%Mixing_Length(i,j,k) = mixlen(k)
      enddo ; endif
      if (CS%id_Velocity_Scale>0) then ; do k=1,nz
        CS%Velocity_Scale(i,j,k) = mixvel(k)
      enddo ; endif
      if (allocated(CS%mstar_mix)) CS%mstar_mix(i,j) = eCD%mstar
      if (allocated(CS%mstar_lt)) CS%mstar_lt(i,j) = eCD%mstar_LT
      if (allocated(CS%La)) CS%La(i,j) = eCD%LA
      if (allocated(CS%La_mod)) CS%La_mod(i,j) = eCD%LAmod
    else ! End of the ocean-point part of the i-loop
      ! For masked points, Kd_int must still be set (to 0) because it has intent out.
      do K=1,nz+1 ; Kd_2d(i,K) = 0. ; enddo
      CS%ML_depth(i,j) = 0.0

      if (present(dT_expected)) then
        do k=1,nz ; dT_expected(i,j,k) = 0.0 ; enddo
      endif
      if (present(dS_expected)) then
        do k=1,nz ; dS_expected(i,j,k) = 0.0 ; enddo
      endif
    endif ; enddo ! Close of i-loop - Note unusual loop order!

    do K=1,nz+1 ; do i=is,ie ; Kd_int(i,j,K) = Kd_2d(i,K) ; enddo ; enddo

  enddo ! j-loop

  if (write_diags) then
    if (CS%id_ML_depth > 0) call post_data(CS%id_ML_depth, CS%ML_depth, CS%diag)
    if (CS%id_TKE_wind > 0) call post_data(CS%id_TKE_wind, CS%diag_TKE_wind, CS%diag)
    if (CS%id_TKE_MKE > 0)  call post_data(CS%id_TKE_MKE, CS%diag_TKE_MKE, CS%diag)
    if (CS%id_TKE_conv > 0) call post_data(CS%id_TKE_conv, CS%diag_TKE_conv, CS%diag)
    if (CS%id_TKE_forcing > 0) call post_data(CS%id_TKE_forcing, CS%diag_TKE_forcing, CS%diag)
    if (CS%id_TKE_mixing > 0) call post_data(CS%id_TKE_mixing, CS%diag_TKE_mixing, CS%diag)
    if (CS%id_TKE_mech_decay > 0) &
      call post_data(CS%id_TKE_mech_decay, CS%diag_TKE_mech_decay, CS%diag)
    if (CS%id_TKE_conv_decay > 0) &
      call post_data(CS%id_TKE_conv_decay, CS%diag_TKE_conv_decay, CS%diag)
    if (CS%id_Mixing_Length > 0) call post_data(CS%id_Mixing_Length, CS%Mixing_Length, CS%diag)
    if (CS%id_Velocity_Scale >0) call post_data(CS%id_Velocity_Scale, CS%Velocity_Scale, CS%diag)
    if (CS%id_MSTAR_MIX > 0)     call post_data(CS%id_MSTAR_MIX, CS%MSTAR_MIX, CS%diag)
    if (CS%id_LA > 0)       call post_data(CS%id_LA, CS%LA, CS%diag)
    if (CS%id_LA_MOD > 0)   call post_data(CS%id_LA_MOD, CS%LA_MOD, CS%diag)
    if (CS%id_MSTAR_LT > 0) call post_data(CS%id_MSTAR_LT, CS%MSTAR_LT, CS%diag)
  endif

  if (debug) deallocate(eCD%dT_expect, eCD%dS_expect)

end subroutine energetic_PBL



!> This subroutine determines the diffusivities from the integrated energetics
!!  mixed layer model for a single column of water.
subroutine ePBL_column(h, u, v, T0, S0, dSV_dT, dSV_dS, TKE_forcing, B_flux, absf, &
                       u_star, u_star_mean, dt, MLD_io, Kd, mixvel, mixlen, GV, US, CS, eCD, &
                       dt_diag, Waves, G, i, j)
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZK_(GV)), intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZK_(GV)), intent(in)  :: u      !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZK_(GV)), intent(in)  :: v      !< Zonal velocities interpolated to h points
                                                   !! [L T-1 ~> m s-1].
  real, dimension(SZK_(GV)), intent(in)  :: T0     !< The initial layer temperatures [degC].
  real, dimension(SZK_(GV)), intent(in)  :: S0     !< The initial layer salinities [ppt].

  real, dimension(SZK_(GV)), intent(in)  :: dSV_dT !< The partial derivative of in-situ specific
                                                   !! volume with potential temperature
                                                   !! [R-1 degC-1 ~> m3 kg-1 degC-1].
  real, dimension(SZK_(GV)), intent(in)  :: dSV_dS !< The partial derivative of in-situ specific
                                                   !! volume with salinity [R-1 ppt-1 ~> m3 kg-1 ppt-1].
  real, dimension(SZK_(GV)), intent(in)  :: TKE_forcing !< The forcing requirements to homogenize the
                                                   !! forcing that has been applied to each layer
                                                   !! [R Z3 T-2 ~> J m-2].
  real,                    intent(in)    :: B_flux !< The surface buoyancy flux [Z2 T-3 ~> m2 s-3]
  real,                    intent(in)    :: absf   !< The absolute value of the Coriolis parameter [T-1 ~> s-1].
  real,                    intent(in)    :: u_star !< The surface friction velocity [Z T-1 ~> m s-1].
  real,                    intent(in)    :: u_star_mean !< The surface friction velocity without any
                                                   !! contribution from unresolved gustiness  [Z T-1 ~> m s-1].
  real,                    intent(inout) :: MLD_io !< A first guess at the mixed layer depth on input, and
                                                   !! the calculated mixed layer depth on output [Z ~> m].
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: Kd     !< The diagnosed diffusivities at interfaces
                                                   !! [Z2 T-1 ~> m2 s-1].
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: mixvel !< The mixing velocity scale used in Kd
                                                   !! [Z T-1 ~> m s-1].
  real, dimension(SZK_(GV)+1), &
                           intent(out)   :: mixlen !< The mixing length scale used in Kd [Z ~> m].
  type(energetic_PBL_CS),  pointer       :: CS     !< The control structure returned by a previous
                                                   !! call to mixedlayer_init.
  type(ePBL_column_diags), intent(inout) :: eCD    !< A container for passing around diagnostics.
  real,          optional, intent(in)    :: dt_diag   !< The diagnostic time step, which may be less
                                                   !! than dt if there are two calls to mixedlayer [T ~> s].
  type(wave_parameters_CS), &
                 optional, pointer       :: Waves  !< Wave CS for Langmuir turbulence
  type(ocean_grid_type), &
                 optional, intent(inout) :: G      !< The ocean's grid structure.
  integer,       optional, intent(in)    :: i      !< The i-index to work on (used for Waves)
  integer,       optional, intent(in)    :: j      !< The i-index to work on (used for Waves)

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
  real :: htot      !   The total depth of the layers above an interface [H ~> m or kg m-2].
  real :: uhtot     !   The depth integrated zonal and meridional velocities in the
  real :: vhtot     ! layers above [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real :: Idecay_len_TKE  ! The inverse of a turbulence decay length scale [H-1 ~> m-1 or m2 kg-1].
  real :: h_sum     ! The total thickness of the water column [H ~> m or kg m-2].

  real, dimension(SZK_(GV)) :: &
    dT_to_dColHt, & ! Partial derivative of the total column height with the temperature changes
                    ! within a layer [Z degC-1 ~> m degC-1].
    dS_to_dColHt, & ! Partial derivative of the total column height with the salinity changes
                    ! within a layer  [Z ppt-1 ~> m ppt-1].
    dT_to_dPE, &    ! Partial derivatives of column potential energy with the temperature
                    ! changes within a layer, in [R Z3 T-2 degC-1 ~> J m-2 degC-1].
    dS_to_dPE, &    ! Partial derivatives of column potential energy with the salinity changes
                    ! within a layer, in [R Z3 T-2 ppt-1 ~> J m-2 ppt-1].
    dT_to_dColHt_a, & ! Partial derivative of the total column height with the temperature changes
                    ! within a layer, including the implicit effects  of mixing with layers higher
                    ! in the water column [Z degC-1 ~> m degC-1].
    dS_to_dColHt_a, & ! Partial derivative of the total column height with the salinity changes
                    ! within a layer, including the implicit effects  of mixing with layers higher
                    ! in the water column [Z ppt-1 ~> m ppt-1].
    dT_to_dPE_a, &  ! Partial derivatives of column potential energy with the temperature changes
                    ! within a layer, including the implicit effects of mixing with layers higher
                    ! in the water column [R Z3 T-2 degC-1 ~> J m-2 degC-1].
    dS_to_dPE_a, &  ! Partial derivative of column potential energy with the salinity changes
                    ! within a layer, including the implicit effects of mixing with layers higher
                    ! in the water column [R Z3 T-2 ppt-1 ~> J m-2 ppt-1].
    c1, &           ! c1 is used by the tridiagonal solver [nondim].
    Te, &           ! Estimated final values of T in the column [degC].
    Se, &           ! Estimated final values of S in the column [ppt].
    dTe, &          ! Running (1-way) estimates of temperature change [degC].
    dSe, &          ! Running (1-way) estimates of salinity change [ppt].
    Th_a, &         ! An effective temperature times a thickness in the layer above, including implicit
                    ! mixing effects with other yet higher layers [degC H ~> degC m or degC kg m-2].
    Sh_a, &         ! An effective salinity times a thickness in the layer above, including implicit
                    ! mixing effects with other yet higher layers [ppt H ~> ppt m or ppt kg m-2].
    Th_b, &         ! An effective temperature times a thickness in the layer below, including implicit
                    ! mixing effects with other yet lower layers [degC H ~> degC m or degC kg m-2].
    Sh_b            ! An effective salinity times a thickness in the layer below, including implicit
                    ! mixing effects with other yet lower layers [ppt H ~> ppt m or ppt kg m-2].
  real, dimension(SZK_(GV)+1) :: &
    MixLen_shape, & ! A nondimensional shape factor for the mixing length that
                    ! gives it an appropriate assymptotic value at the bottom of
                    ! the boundary layer.
    Kddt_h          ! The diapycnal diffusivity times a timestep divided by the
                    ! average thicknesses around a layer [H ~> m or kg m-2].
  real :: b1        ! b1 is inverse of the pivot used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
  real :: hp_a      ! An effective pivot thickness of the layer including the effects
                    ! of coupling with layers above [H ~> m or kg m-2].  This is the first term
                    ! in the denominator of b1 in a downward-oriented tridiagonal solver.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: dMass     ! The mass per unit area within a layer [Z R ~> kg m-2].
  real :: dPres     ! The hydrostatic pressure change across a layer [R Z2 T-2 ~> kg m-1 s-2 = Pa = J m-3].
  real :: dMKE_max  ! The maximum amount of mean kinetic energy that could be
                    ! converted to turbulent kinetic energy if the velocity in
                    ! the layer below an interface were homogenized with all of
                    ! the water above the interface [R Z3 T-2 ~> J m-2].
  real :: MKE2_Hharm ! Twice the inverse of the harmonic mean of the thickness
                    ! of a layer and the thickness of the water above, used in
                    ! the MKE conversion equation [H-1 ~> m-1 or m2 kg-1].

  real :: dt_h      ! The timestep divided by the averages of the thicknesses around
                    ! a layer, times a thickness conversion factor [H T m-2 ~> s m-1 or kg s m-4].
  real :: h_bot     ! The distance from the bottom [H ~> m or kg m-2].
  real :: h_rsum    ! The running sum of h from the top [Z ~> m].
  real :: I_hs      ! The inverse of h_sum [H-1 ~> m-1 or m2 kg-1].
  real :: I_MLD     ! The inverse of the current value of MLD [Z-1 ~> m-1].
  real :: h_tt      ! The distance from the surface or up to the next interface
                    ! that did not exhibit turbulent mixing from this scheme plus
                    ! a surface mixing roughness length given by h_tt_min [H ~> m or kg m-2].
  real :: h_tt_min  ! A surface roughness length [H ~> m or kg m-2].

  real :: C1_3      ! = 1/3.
  real :: I_dtrho   ! 1.0 / (dt * Rho0) times conversion factors in [m3 Z-3 R-1 T2 s-3 ~> m3 kg-1 s-1].
                    ! This is used convert TKE back into ustar^3.
  real :: vstar     ! An in-situ turbulent velocity [Z T-1 ~> m s-1].
  real :: mstar_total ! The value of mstar used in ePBL [nondim]
  real :: mstar_LT  ! An addition to mstar due to Langmuir turbulence [nondim] (output for diagnostic)
  real :: MLD_output ! The mixed layer depth output from this routine [Z ~> m].
  real :: LA        ! The value of the Langmuir number [nondim]
  real :: LAmod     ! The modified Langmuir number by convection [nondim]
  real :: hbs_here  ! The local minimum of hb_hs and MixLen_shape, times a
                    ! conversion factor from H to Z [Z H-1 ~> 1 or m3 kg-1].
  real :: nstar_FC  ! The fraction of conv_PErel that can be converted to mixing [nondim].
  real :: TKE_reduc ! The fraction by which TKE and other energy fields are
                    ! reduced to support mixing [nondim]. between 0 and 1.
  real :: tot_TKE   ! The total TKE available to support mixing at interface K [R Z3 T-2 ~> J m-2].
  real :: TKE_here  ! The total TKE at this point in the algorithm [R Z3 T-2 ~> J m-2].
  real :: dT_km1_t2 ! A diffusivity-independent term related to the temperature
                    ! change in the layer above the interface [degC].
  real :: dS_km1_t2 ! A diffusivity-independent term related to the salinity
                    ! change in the layer above the interface [ppt].
  real :: dTe_term  ! A diffusivity-independent term related to the temperature
                    ! change in the layer below the interface [degC H ~> degC m or degC kg m-2].
  real :: dSe_term  ! A diffusivity-independent term related to the salinity
                    ! change in the layer above the interface [ppt H ~> ppt m or ppt kg m-2].
  real :: dTe_t2    ! A part of dTe_term [degC H ~> degC m or degC kg m-2].
  real :: dSe_t2    ! A part of dSe_term [ppt H ~> ppt m or ppt kg m-2].
  real :: dPE_conv  ! The convective change in column potential energy [R Z3 T-2 ~> J m-2].
  real :: MKE_src   ! The mean kinetic energy source of TKE due to Kddt_h(K) [R Z3 T-2 ~> J m-2].
  real :: dMKE_src_dK  ! The partial derivative of MKE_src with Kddt_h(K) [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real :: Kd_guess0    ! A first guess of the diapycnal diffusivity [Z2 T-1 ~> m2 s-1].
  real :: PE_chg_g0    ! The potential energy change when Kd is Kd_guess0 [R Z3 T-2 ~> J m-2]
  !### The following might be unused.
  real :: dPEa_dKd_g0  ! The derivative of the change in the potential energy of the column above an interface
                       ! with the diffusivity when the Kd is Kd_guess0 [R Z T-1 ~> J s m-4]
  real :: Kddt_h_g0    ! The first guess diapycnal diffusivity times a timestep divided
                       ! by the average thicknesses around a layer [H ~> m or kg m-2].
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
  real :: vstar_unit_scale ! A unit converion factor for turbulent velocities [Z T-1 s m-1 ~> 1]
  logical :: use_Newt  ! Use Newton's method for the next guess at Kddt_h(K).
  logical :: convectively_stable ! If true the water column is convectively stable at this interface.
  logical :: sfc_connected   ! If true the ocean is actively turbulent from the present
                    ! interface all the way up to the surface.
  logical :: sfc_disconnect ! If true, any turbulence has become disconnected
                    ! from the surface.

! The following are only used for diagnostics.
  real :: dt__diag  ! A copy of dt_diag (if present) or dt [T ~> s].
  real :: I_dtdiag  !  = 1.0 / dt__diag [T-1 ~> s-1].

  !----------------------------------------------------------------------
  !/BGR added Aug24,2016 for adding iteration to get boundary layer depth
  !    - needed to compute new mixing length.
  real :: MLD_guess, MLD_found ! Mixing Layer depth guessed/found for iteration [Z ~> m].
  real :: min_MLD   ! Iteration bounds [Z ~> m], which are adjusted at each step
  real :: max_MLD   !  - These are initialized based on surface/bottom
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
  logical :: FIRST_OBL     ! Flag for computing "found" Mixing layer depth
  logical :: OBL_converged ! Flag for convergence of MLD
  integer :: OBL_it        ! Iteration counter

  real :: Surface_Scale ! Surface decay scale for vstar

  logical :: debug=.false.  ! Change this hard-coded value for debugging.

  !  The following arrays are used only for debugging purposes.
  real :: dPE_debug, mixing_debug, taux2, tauy2
  real, dimension(20) :: TKE_left_itt, PE_chg_itt, Kddt_h_itt, dPEa_dKd_itt, MKE_src_itt
  real, dimension(SZK_(GV)) :: mech_TKE_k, conv_PErel_k, nstar_k
  integer, dimension(SZK_(GV)) :: num_itts

  integer :: k, nz, itt, max_itt

  nz = GV%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "energetic_PBL: "//&
         "Module must be initialized before it is used.")

  debug = .false. ; if (allocated(eCD%dT_expect) .or. allocated(eCD%dS_expect)) debug = .true.

  h_neglect = GV%H_subroundoff

  C1_3 = 1.0 / 3.0
  dt__diag = dt ; if (present(dt_diag)) dt__diag = dt_diag
  I_dtdiag = 1.0 / dt__diag
  max_itt = 20

  h_tt_min = 0.0
  I_dtrho = 0.0 ; if (dt*GV%Rho0 > 0.0) I_dtrho = (US%Z_to_m**3*US%s_to_T**3) / (dt*GV%Rho0)
  vstar_unit_scale = US%m_to_Z * US%T_to_s

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
    dPres = US%L_to_Z**2 * GV%g_Earth * dMass  ! Equivalent to GV%H_to_Pa * h(k) with rescaling
    dT_to_dPE(k) = (dMass * (pres_Z(K) + 0.5*dPres)) * dSV_dT(k)
    dS_to_dPE(k) = (dMass * (pres_Z(K) + 0.5*dPres)) * dSV_dS(k)
    dT_to_dColHt(k) = dMass * dSV_dT(k)
    dS_to_dColHt(k) = dMass * dSV_dS(k)

    pres_Z(K+1) = pres_Z(K) + dPres
  enddo

  ! Determine the total thickness (h_sum) and the fractional distance from the bottom (hb_hs).
  h_sum = H_neglect ; do k=1,nz ; h_sum = h_sum + h(k) ; enddo
  I_hs = 0.0 ; if (h_sum > 0.0) I_hs = 1.0 / h_sum
  h_bot = 0.0
  hb_hs(nz+1) = 0.0
  do k=nz,1,-1
    h_bot = h_bot + h(k)
    hb_hs(K) = h_bot * I_hs
  enddo

  MLD_output = h(1)*GV%H_to_Z

  !/The following lines are for the iteration over MLD
  ! max_MLD will initialized as ocean bottom depth
  max_MLD = 0.0 ; do k=1,nz ; max_MLD = max_MLD + h(k)*GV%H_to_Z ; enddo
  !min_MLD will initialize as 0.
  min_MLD = 0.0

  ! If no first guess is provided for MLD, try the middle of the water column
  if (MLD_guess <= min_MLD) MLD_guess = 0.5 * (min_MLD + max_MLD)

  ! Iterate to determine a converged EPBL depth.
  OBL_converged = .false.
  do OBL_it=1,CS%Max_MLD_Its

    if (.not. OBL_converged) then
      ! If not using MLD_Iteration flag loop to only execute once.
      if (.not.CS%Use_MLD_iteration) OBL_converged = .true.

      if (debug) then ; mech_TKE_k(:) = 0.0 ; conv_PErel_k(:) = 0.0 ; endif

      ! Reset ML_depth
      MLD_output = h(1)*GV%H_to_Z
      sfc_connected = .true.

      !/ Here we get MStar, which is the ratio of convective TKE driven mixing to UStar**3
      if (CS%Use_LT) then
        call get_Langmuir_Number(LA, G, GV, US, abs(MLD_guess), u_star_mean, i, j, &
                                 H=h, U_H=u, V_H=v, Waves=Waves)
        call find_mstar(CS, US, B_flux, u_star, u_star_Mean, MLD_Guess, absf, &
                        MStar_total, Langmuir_Number=La, Convect_Langmuir_Number=LAmod,&
                        mstar_LT=mstar_LT)
      else
        call find_mstar(CS, US, B_flux, u_star, u_star_mean, MLD_guess, absf, mstar_total)
      endif

      !/ Apply MStar to get mech_TKE
      if ((CS%answers_2018) .and. (CS%mstar_scheme==Use_Fixed_MStar)) then
        mech_TKE = (dt*MSTAR_total*GV%Rho0) * u_star**3
      else
        mech_TKE = MSTAR_total * (dt*GV%Rho0* u_star**3)
      endif

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
        h_rsum = 0.0
        MixLen_shape(1) = 1.0
        do K=2,nz+1
          h_rsum = h_rsum + h(k-1)*GV%H_to_Z
          if (CS%MixLenExponent==2.0) then
            MixLen_shape(K) = CS%transLay_scale + (1.0 - CS%transLay_scale) * &
                 (max(0.0, (MLD_guess - h_rsum)*I_MLD) )**2 ! CS%MixLenExponent
          else
            MixLen_shape(K) = CS%transLay_scale + (1.0 - CS%transLay_scale) * &
                 (max(0.0, (MLD_guess - h_rsum)*I_MLD) )**CS%MixLenExponent
          endif
        enddo
      endif

      Kd(1) = 0.0 ; Kddt_h(1) = 0.0
      hp_a = h(1)
      dT_to_dPE_a(1) = dT_to_dPE(1) ; dT_to_dColHt_a(1) = dT_to_dColHt(1)
      dS_to_dPE_a(1) = dS_to_dPE(1) ; dS_to_dColHt_a(1) = dS_to_dColHt(1)

      htot = h(1) ; uhtot = u(1)*h(1) ; vhtot = v(1)*h(1)

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
        Idecay_len_TKE = (CS%TKE_decay * absf / u_star) * GV%H_to_Z
        exp_kh = 1.0
        if (Idecay_len_TKE > 0.0) exp_kh = exp(-h(k-1)*Idecay_len_TKE)
        if (CS%TKE_diagnostics) &
          eCD%dTKE_mech_decay = eCD%dTKE_mech_decay + (exp_kh-1.0) * mech_TKE * I_dtdiag
        mech_TKE = mech_TKE * exp_kh

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
          ! Note:         Ro = 1.0 / sqrt(0.5 * dt * Rho0 * (absf*htot)**3 / conv_PErel)
          nstar_FC = CS%nstar * conv_PErel / (conv_PErel + 0.2 * &
                     sqrt(0.5 * dt * GV%Rho0 * (absf*(htot*GV%H_to_Z))**3 * conv_PErel))
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
        dt_h = (GV%Z_to_H**2*dt) / max(0.5*(h(k-1)+h(k)), 1e-15*h_sum)

        !   This tests whether the layers above and below this interface are in
        ! a convetively stable configuration, without considering any effects of
        ! mixing at higher interfaces.  It is an approximation to the more
        ! complete test dPEc_dKd_Kd0 >= 0.0, that would include the effects of
        ! mixing across interface K-1.  The dT_to_dColHt here are effectively
        ! mass-weigted estimates of dSV_dT.
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
          b1 = 1.0 / hp_a
          c1(K) = 0.0
          if (CS%orig_PE_calc) then
            dTe(k-1) = b1 * ( dTe_t2 )
            dSe(k-1) = b1 * ( dSe_t2 )
          endif

          hp_a = h(k)
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
                    (Kddt_h(K-1) / hp_a) * ((T0(k-2) - T0(k-1)) + dTe(k-2))
              dS_km1_t2 = (S0(k)-S0(k-1)) - &
                    (Kddt_h(K-1) / hp_a) * ((S0(k-2) - S0(k-1)) + dSe(k-2))
            endif
            dTe_term = dTe_t2 + hp_a * (T0(k-1)-T0(k))
            dSe_term = dSe_t2 + hp_a * (S0(k-1)-S0(k))
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
            dMKE_max = (US%L_to_Z**2*GV%H_to_RZ * CS%MKE_to_TKE_effic) * 0.5 * &
                (h(k) / ((htot + h(k))*htot)) * &
                ((uhtot-u(k)*htot)**2 + (vhtot-v(k)*htot)**2)
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
          h_tt = htot + h_tt_min
          TKE_here = mech_TKE + CS%wstar_ustar_coef*conv_PErel
          if (TKE_here > 0.0) then
            if (CS%wT_scheme==wT_from_cRoot_TKE) then
              vstar = CS%vstar_scale_fac * vstar_unit_scale * (I_dtrho*TKE_here)**C1_3
            elseif (CS%wT_scheme==wT_from_RH18) then
              Surface_Scale = max(0.05, 1.0 - htot/MLD_guess)
              vstar = CS%vstar_scale_fac * Surface_Scale * (CS%vstar_surf_fac*u_star + &
                        vstar_unit_scale * (CS%wstar_ustar_coef*conv_PErel*I_dtrho)**C1_3)
            endif
            hbs_here = GV%H_to_Z * min(hb_hs(K), MixLen_shape(K))
            mixlen(K) = MAX(CS%min_mix_len, ((h_tt*hbs_here)*vstar) / &
                ((CS%Ekman_scale_coef * absf) * (h_tt*hbs_here) + vstar))
            !Note setting Kd_guess0 to vstar * CS%vonKar * mixlen(K) here will
            ! change the answers.  Therefore, skipping that.
            if (.not.CS%Use_MLD_iteration) then
              Kd_guess0 = vstar * CS%vonKar * ((h_tt*hbs_here)*vstar) / &
                ((CS%Ekman_scale_coef * absf) * (h_tt*hbs_here) + vstar)
            else
              Kd_guess0 = vstar * CS%vonKar * mixlen(K)
            endif
          else
            vstar = 0.0 ; Kd_guess0 = 0.0
          endif
          mixvel(K) = vstar ! Track vstar
          Kddt_h_g0 = Kd_guess0 * dt_h

          if (CS%orig_PE_calc) then
            call find_PE_chg_orig(Kddt_h_g0, h(k), hp_a, dTe_term, dSe_term, &
                     dT_km1_t2, dS_km1_t2, dT_to_dPE(k), dS_to_dPE(k), &
                     dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), &
                     pres_Z(K), dT_to_dColHt(k), dS_to_dColHt(k), &
                     dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                     PE_chg=PE_chg_g0, dPEc_dKd=dPEa_dKd_g0, dPE_max=PE_chg_max, &
                     dPEc_dKd_0=dPEc_dKd_Kd0 )
          else
            call find_PE_chg(0.0, Kddt_h_g0, hp_a, h(k), &
                     Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                     dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE(k), dS_to_dPE(k), &
                     pres_Z(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                     dT_to_dColHt(k), dS_to_dColHt(k), &
                     PE_chg=PE_chg_g0, dPEc_dKd=dPEa_dKd_g0, dPE_max=PE_chg_max, &
                     dPEc_dKd_0=dPEc_dKd_Kd0 )
          endif

          MKE_src = dMKE_max*(1.0 - exp(-Kddt_h_g0 * MKE2_Hharm))

          ! This block checks out different cases to determine Kd at the present interface.
          if ((PE_chg_g0 < 0.0) .or. ((vstar == 0.0) .and. (dPEc_dKd_Kd0 < 0.0))) then
            ! This column is convectively unstable.
            if (PE_chg_max <= 0.0) then
              ! Does MKE_src need to be included in the calculation of vstar here?
              TKE_here = mech_TKE + CS%wstar_ustar_coef*(conv_PErel-PE_chg_max)
              if (TKE_here > 0.0) then
                if (CS%wT_scheme==wT_from_cRoot_TKE) then
                  vstar = CS%vstar_scale_fac * vstar_unit_scale * (I_dtrho*TKE_here)**C1_3
                elseif (CS%wT_scheme==wT_from_RH18) then
                  Surface_Scale = max(0.05, 1. - htot/MLD_guess)
                  vstar = CS%vstar_scale_fac * Surface_Scale * (CS%vstar_surf_fac*u_star + &
                                  vstar_unit_scale * (CS%wstar_ustar_coef*conv_PErel*I_dtrho)**C1_3)
                endif
                hbs_here = GV%H_to_Z * min(hb_hs(K), MixLen_shape(K))
                mixlen(K) = max(CS%min_mix_len, ((h_tt*hbs_here)*vstar) / &
                    ((CS%Ekman_scale_coef * absf) * (h_tt*hbs_here) + vstar))
                if (.not.CS%Use_MLD_iteration) then
                ! Note again (as prev) that using mixlen here
                !  instead of redoing the computation will change answers...
                  Kd(K) = vstar * CS%vonKar *  ((h_tt*hbs_here)*vstar) / &
                        ((CS%Ekman_scale_coef * absf) * (h_tt*hbs_here) + vstar)
                else
                  Kd(K) = vstar * CS%vonKar * mixlen(K)
                endif
              else
                vstar = 0.0 ; Kd(K) = 0.0
              endif
              mixvel(K) = vstar

              if (CS%orig_PE_calc) then
                call find_PE_chg_orig(Kd(K)*dt_h, h(k), hp_a, dTe_term, dSe_term, &
                         dT_km1_t2, dS_km1_t2, dT_to_dPE(k), dS_to_dPE(k), &
                         dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), &
                         pres_Z(K), dT_to_dColHt(k), dS_to_dColHt(k), &
                         dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                         PE_chg=dPE_conv)
              else
                call find_PE_chg(0.0, Kd(K)*dt_h, hp_a, h(k), &
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
              MLD_output = MLD_output + GV%H_to_Z * h(k)
            endif

            Kddt_h(K) = Kd(K) * dt_h
          elseif (tot_TKE + (MKE_src - PE_chg_g0) >= 0.0) then
            ! This column is convctively stable and there is energy to support the suggested
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
              MLD_output = MLD_output + GV%H_to_Z * h(k)
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
                call find_PE_chg_orig(Kddt_h_guess, h(k), hp_a, dTe_term, dSe_term, &
                         dT_km1_t2, dS_km1_t2, dT_to_dPE(k), dS_to_dPE(k), &
                         dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), &
                         pres_Z(K), dT_to_dColHt(k), dS_to_dColHt(k), &
                         dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                         PE_chg=PE_chg, dPEc_dKd=dPEc_dKd )
              else
                call find_PE_chg(0.0, Kddt_h_guess, hp_a, h(k), &
                         Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                         dT_to_dPE_a(k-1), dS_to_dPE_a(k-1), dT_to_dPE(k), dS_to_dPE(k), &
                         pres_Z(K), dT_to_dColHt_a(k-1), dS_to_dColHt_a(k-1), &
                         dT_to_dColHt(k), dS_to_dColHt(k), &
                         PE_chg=dPE_conv)
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

            if (sfc_connected) MLD_output = MLD_output + &
                 (PE_chg / (PE_chg_g0)) * GV%H_to_Z * h(k)

            tot_TKE = 0.0 ; mech_TKE = 0.0 ; conv_PErel = 0.0
            sfc_disconnect = .true.
          endif ! End of convective or forced mixing cases to determine Kd.

          Kddt_h(K) = Kd(K) * dt_h
          !   At this point, the final value of Kddt_h(K) is known, so the
          ! estimated properties for layer k-1 can be calculated.
          b1 = 1.0 / (hp_a + Kddt_h(K))
          c1(K) = Kddt_h(K) * b1
          if (CS%orig_PE_calc) then
            dTe(k-1) = b1 * ( Kddt_h(K)*(T0(k)-T0(k-1)) + dTe_t2 )
            dSe(k-1) = b1 * ( Kddt_h(K)*(S0(k)-S0(k-1)) + dSe_t2 )
          endif

          hp_a = h(k) + (hp_a * b1) * Kddt_h(K)
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
          sfc_connected = .false.
        else
          uhtot = uhtot + u(k)*h(k)
          vhtot = vhtot + v(k)*h(k)
          htot  = htot + h(k)
        endif

        if (debug) then
          if (k==2) then
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
        b1 = 1.0 / hp_a
        Te(nz) = b1 * (h(nz) * T0(nz) + Kddt_h(nz) * Te(nz-1))
        Se(nz) = b1 * (h(nz) * S0(nz) + Kddt_h(nz) * Se(nz-1))
        eCD%dT_expect(nz) = Te(nz) - T0(nz) ; eCD%dS_expect(nz) = Se(nz) - S0(nz)
        do k=nz-1,1,-1
          Te(k) = Te(k) + c1(K+1)*Te(k+1)
          Se(k) = Se(k) + c1(K+1)*Se(k+1)
          eCD%dT_expect(k) = Te(k) - T0(k) ; eCD%dS_expect(k) = Se(k) - S0(k)
        enddo

        dPE_debug = 0.0
        do k=1,nz
          dPE_debug = dPE_debug + (dT_to_dPE(k) * (Te(k) - T0(k)) + &
                                   dS_to_dPE(k) * (Se(k) - S0(k)))
        enddo
        mixing_debug = dPE_debug * I_dtdiag
      endif
      k = nz ! This is here to allow a breakpoint to be set.
      !/BGR
      ! The following lines are used for the iteration
      ! note the iteration has been altered to use the value predicted by
      ! the TKE threshold (ML_DEPTH).  This is because the MSTAR
      ! is now dependent on the ML, and therefore the ML needs to be estimated
      ! more precisely than the grid spacing.

      !New method uses ML_DEPTH as computed in ePBL routine
      MLD_found = MLD_output
      if (MLD_found - CS%MLD_tol > MLD_guess) then
        min_MLD = MLD_guess
      elseif (abs(MLD_guess - MLD_found) < CS%MLD_tol) then
        OBL_converged = .true. ! Break convergence loop
      else
        max_MLD = MLD_guess ! We know this guess was too deep
      endif

      ! For next pass, guess average of minimum and maximum values.
      !### We should try using the false position method instead of simple bisection.
      MLD_guess = 0.5*(min_MLD + max_MLD)
    endif
  enddo ! Iteration loop for converged boundary layer thickness.
  if (CS%Use_LT) then
    eCD%LA = LA ; eCD%LAmod = LAmod ; eCD%mstar = mstar_total ; eCD%mstar_LT = mstar_LT
  else
    eCD%LA = 0.0 ; eCD%LAmod = 0.0 ; eCD%mstar = mstar_total ; eCD%mstar_LT = 0.0
  endif

  MLD_io = MLD_output

end subroutine ePBL_column

!> This subroutine calculates the change in potential energy and or derivatives
!! for several changes in an interfaces's diapycnal diffusivity times a timestep.
subroutine find_PE_chg(Kddt_h0, dKddt_h, hp_a, hp_b, Th_a, Sh_a, Th_b, Sh_b, &
                       dT_to_dPE_a, dS_to_dPE_a, dT_to_dPE_b, dS_to_dPE_b, &
                       pres_Z, dT_to_dColHt_a, dS_to_dColHt_a, dT_to_dColHt_b, dS_to_dColHt_b, &
                       PE_chg, dPEc_dKd, dPE_max, dPEc_dKd_0, PE_ColHt_cor)
  real, intent(in)  :: Kddt_h0  !< The previously used diffusivity at an interface times
                                !! the time step and  divided by the average of the
                                !! thicknesses around the interface [H ~> m or kg m-2].
  real, intent(in)  :: dKddt_h  !< The trial change in the diffusivity at an interface times
                                !! the time step and  divided by the average of the
                                !! thicknesses around the interface [H ~> m or kg m-2].
  real, intent(in)  :: hp_a     !< The effective pivot thickness of the layer above the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above [H ~> m or kg m-2].
  real, intent(in)  :: hp_b     !< The effective pivot thickness of the layer below the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above [H ~> m or kg m-2].
  real, intent(in)  :: Th_a     !< An effective temperature times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [degC H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_a     !< An effective salinity times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers [degC H ~> degC m or degC kg m-2].
  real, intent(in)  :: Th_b     !< An effective temperature times a thickness in the layer
                                !! below, including implicit mixfing effects with other
                                !! yet lower layers [degC H ~> degC m or degC kg m-2].
  real, intent(in)  :: Sh_b     !< An effective salinity times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers [degC H ~> degC m or degC kg m-2].
  real, intent(in)  :: dT_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers above [R Z3 T-2 degC-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers above [R Z3 T-2 ppt-1 ~> J m-2 ppt-1].
  real, intent(in)  :: dT_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! temperatures of all the layers below [R Z3 T-2 degC-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column potential
                                !! energy, including all implicit diffusive changes in the
                                !! salinities of all the layers below [R Z3 T-2 ppt-1 ~> J m-2 ppt-1].
  real, intent(in)  :: pres_Z   !< The rescaled hydrostatic interface pressure, which relates
                                !! the changes in column thickness to the energy that is radiated
                                !! as gravity waves and unavailable to drive mixing [R Z2 T-2 ~> J m-3].
  real, intent(in)  :: dT_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers above [Z degC-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers above [Z ppt-1 ~> m ppt-1].
  real, intent(in)  :: dT_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers below [Z degC-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers below [Z ppt-1 ~> m ppt-1].

  real, optional, intent(out) :: PE_chg   !< The change in column potential energy from applying
                                          !! Kddt_h at the present interface [R Z3 T-2 ~> J m-2].
  real, optional, intent(out) :: dPEc_dKd !< The partial derivative of PE_chg with Kddt_h
                                          !! [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real, optional, intent(out) :: dPE_max  !< The maximum change in column potential energy that could
                                          !! be realizedd by applying a huge value of Kddt_h at the
                                          !! present interface [R Z3 T-2 ~> J m-2].
  real, optional, intent(out) :: dPEc_dKd_0 !< The partial derivative of PE_chg with Kddt_h in the
                                            !! limit where Kddt_h = 0 [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real, optional, intent(out) :: PE_ColHt_cor !< The correction to PE_chg that is made due to a net
                                            !! change in the column height [R Z3 T-2 ~> J m-2].

  real :: hps ! The sum of the two effective pivot thicknesses [H ~> m or kg m-2].
  real :: bdt1 ! A product of the two pivot thicknesses plus a diffusive term [H2 ~> m2 or kg2 m-4].
  real :: dT_c ! The core term in the expressions for the temperature changes [degC H2 ~> degC m2 or degC kg2 m-4].
  real :: dS_c ! The core term in the expressions for the salinity changes [ppt H2 ~> ppt m2 or ppt kg2 m-4].
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
  ! succint form. The derivation is not necessarily obvious, but it demonstrably
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

  if (present(PE_chg)) then
    ! Find the change in column potential energy due to the change in the
    ! diffusivity at this interface by dKddt_h.
    y1_3 = dKddt_h / (bdt1 * (bdt1 + dKddt_h * hps))
    PE_chg = PEc_core * y1_3
    ColHt_chg = ColHt_core * y1_3
    if (ColHt_chg < 0.0) PE_chg = PE_chg - pres_Z * ColHt_chg
    if (present(PE_ColHt_cor)) PE_ColHt_cor = -pres_Z * min(ColHt_chg, 0.0)
  elseif (present(PE_ColHt_cor)) then
    y1_3 = dKddt_h / (bdt1 * (bdt1 + dKddt_h * hps))
    PE_ColHt_cor = -pres_Z * min(ColHt_core * y1_3, 0.0)
  endif

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

!> This subroutine calculates the change in potential energy and or derivatives
!! for several changes in an interfaces's diapycnal diffusivity times a timestep
!! using the original form used in the first version of ePBL.
subroutine find_PE_chg_orig(Kddt_h, h_k, b_den_1, dTe_term, dSe_term, &
                       dT_km1_t2, dS_km1_t2, dT_to_dPE_k, dS_to_dPE_k, &
                       dT_to_dPEa, dS_to_dPEa, pres_Z, dT_to_dColHt_k, &
                       dS_to_dColHt_k, dT_to_dColHta, dS_to_dColHta, &
                       PE_chg, dPEc_dKd, dPE_max, dPEc_dKd_0)
  real, intent(in)  :: Kddt_h   !< The diffusivity at an interface times the time step and
                                !! divided by the average of the thicknesses around the
                                !! interface [H ~> m or kg m-2].
  real, intent(in)  :: h_k      !< The thickness of the layer below the interface [H ~> m or kg m-2].
  real, intent(in)  :: b_den_1  !< The first term in the denominator of the pivot
                                !! for the tridiagonal solver, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above [H ~> m or kg m-2].
  real, intent(in)  :: dTe_term !< A diffusivity-independent term related to the temperature change
                                !! in the layer below the interface [degC H ~> degC m or degC kg m-2].
  real, intent(in)  :: dSe_term !< A diffusivity-independent term related to the salinity change
                                !! in the layer below the interface [ppt H ~> ppt m or ppt kg m-2].
  real, intent(in)  :: dT_km1_t2 !< A diffusivity-independent term related to the
                                 !! temperature change in the layer above the interface [degC].
  real, intent(in)  :: dS_km1_t2 !< A diffusivity-independent term related to the
                                 !! salinity change in the layer above the interface [ppt].
  real, intent(in)  :: pres_Z    !< The rescaled hydrostatic interface pressure, which relates
                                 !! the changes in column thickness to the energy that is radiated
                                 !! as gravity waves and unavailable to drive mixing [R Z2 T-2 ~> J m-3].
  real, intent(in)  :: dT_to_dPE_k !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                 !! a layer's temperature change to the change in column potential
                                 !! energy, including all implicit diffusive changes in the
                                 !! temperatures of all the layers below [R Z3 T-2 degC-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPE_k !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                 !! a layer's salinity change to the change in column potential
                                 !! energy, including all implicit diffusive changes in the
                                 !! in the salinities of all the layers below [R Z3 T-2 ppt-1 ~> J m-2 ppt-1].
  real, intent(in)  :: dT_to_dPEa !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                 !! a layer's temperature change to the change in column potential
                                 !! energy, including all implicit diffusive changes in the
                                 !! temperatures of all the layers above [R Z3 T-2 degC-1 ~> J m-2 degC-1].
  real, intent(in)  :: dS_to_dPEa !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                 !! a layer's salinity change to the change in column potential
                                 !! energy, including all implicit diffusive changes in the
                                 !! salinities of all the layers above [R Z3 T-2 ppt-1 ~> J m-2 ppt-1].
  real, intent(in)  :: dT_to_dColHt_k !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! height, including all implicit diffusive changes in the
                                 !! temperatures of all the layers below [Z degC-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHt_k !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the salinities of all the layers below [Z ppt-1 ~> m ppt-1].
  real, intent(in)  :: dT_to_dColHta !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the temperatures of all the layers above [Z degC-1 ~> m degC-1].
  real, intent(in)  :: dS_to_dColHta !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the salinities of all the layers above [Z ppt-1 ~> m ppt-1].

  real, optional, intent(out) :: PE_chg   !< The change in column potential energy from applying
                                          !! Kddt_h at the present interface [R Z3 T-2 ~> J m-2].
  real, optional, intent(out) :: dPEc_dKd !< The partial derivative of PE_chg with Kddt_h
                                          !! [R Z3 T-2 H-1 ~> J m-3 or J kg-1].
  real, optional, intent(out) :: dPE_max  !< The maximum change in column potential energy that could
                                          !! be realizedd by applying a huge value of Kddt_h at the
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

  real :: b1            ! b1 is used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
  real :: b1Kd          ! Temporary array [nondim]
  real :: ColHt_chg     ! The change in column thickness [Z ~> m].
  real :: dColHt_max    ! The change in column thickness for infinite diffusivity [Z ~> m].
  real :: dColHt_dKd    ! The partial derivative of column thickness with diffusivity [s Z-1 ~> s m-1].
  real :: dT_k, dT_km1  ! Temporary arrays [degC].
  real :: dS_k, dS_km1  ! Temporary arrays [ppt].
  real :: I_Kr_denom    ! Temporary array [H-2 ~> m-2 or m4 kg-2]
  real :: dKr_dKd       ! Nondimensional temporary array.
  real :: ddT_k_dKd, ddT_km1_dKd ! Temporary arrays [degC H-1 ~> m-1 or m2 kg-1].
  real :: ddS_k_dKd, ddS_km1_dKd ! Temporary arrays [ppt H-1 ~> ppt m-1 or ppt m2 kg-1].

  b1 = 1.0 / (b_den_1 + Kddt_h)
  b1Kd = Kddt_h*b1

  ! Start with the temperature change in layer k-1 due to the diffusivity at
  ! interface K without considering the effects of changes in layer k.

  ! Calculate the change in PE due to the diffusion at interface K
  ! if Kddt_h(K+1) = 0.
  I_Kr_denom = 1.0 / (h_k*b_den_1 + (b_den_1 + h_k)*Kddt_h)

  dT_k = (Kddt_h*I_Kr_denom) * dTe_term
  dS_k = (Kddt_h*I_Kr_denom) * dSe_term

  if (present(PE_chg)) then
    ! Find the change in energy due to diffusion with strength Kddt_h at this interface.
    ! Increment the temperature changes in layer k-1 due the changes in layer k.
    dT_km1 = b1Kd * ( dT_k + dT_km1_t2 )
    dS_km1 = b1Kd * ( dS_k + dS_km1_t2 )
    PE_chg = (dT_to_dPE_k * dT_k + dT_to_dPEa * dT_km1) + &
             (dS_to_dPE_k * dS_k + dS_to_dPEa * dS_km1)
    ColHt_chg = (dT_to_dColHt_k * dT_k + dT_to_dColHta * dT_km1) + &
                (dS_to_dColHt_k * dS_k + dS_to_dColHta * dS_km1)
    if (ColHt_chg < 0.0) PE_chg = PE_chg - pres_Z * ColHt_chg
  endif

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

!> This subroutine finds the Mstar value for ePBL
subroutine find_mstar(CS, US, Buoyancy_Flux, UStar, UStar_Mean,&
                      BLD, Abs_Coriolis, MStar, Langmuir_Number,&
                      MStar_LT, Convect_Langmuir_Number)
  type(energetic_PBL_CS), pointer    :: CS    !< Energetic_PBL control structure.
  type(unit_scale_type), intent(in)  :: US    !< A dimensional unit scaling type
  real,                  intent(in)  :: UStar !< ustar w/ gustiness [Z T-1 ~> m s-1]
  real,                  intent(in)  :: UStar_Mean !< ustar w/o gustiness [Z T-1 ~> m s-1]
  real,                  intent(in)  :: Abs_Coriolis !< abolute value of the Coriolis parameter [T-1 ~> s-1]
  real,                  intent(in)  :: Buoyancy_Flux !< Buoyancy flux [Z2 T-3 ~> m2 s-3]
  real,                  intent(in)  :: BLD   !< boundary layer depth [Z ~> m]
  real,                  intent(out) :: Mstar !< Ouput mstar (Mixing/ustar**3) [nondim]
  real,        optional, intent(in)  :: Langmuir_Number !< Langmuir number [nondim]
  real,        optional, intent(out) :: MStar_LT !< Mstar increase due to Langmuir turbulence [nondim]
  real,        optional, intent(out) :: Convect_Langmuir_number !< Langmuir number including buoyancy flux [nondim]

  !/ Variables used in computing mstar
  real :: MSN_term       ! Temporary terms [nondim]
  real :: MSCR_term1, MSCR_term2 ! Temporary terms [Z3 T-3 ~> m3 s-3]
  real :: MStar_Conv_Red ! Adjustment made to mstar due to convection reducing mechanical mixing [nondim]
  real :: MStar_S, MStar_N ! Mstar in (S)tabilizing/(N)ot-stabilizing buoyancy flux [nondim]

  !/  Integer options for how to find mstar

  !/

  if (CS%mstar_scheme == Use_Fixed_MStar) then
    MStar = CS%Fixed_MStar
  !/ 1. Get mstar
  elseif (CS%mstar_scheme == MStar_from_Ekman) then

    if (CS%answers_2018) then
      ! The limit for the balance of rotation and stabilizing is f(L_Ekman,L_Obukhov)
      MStar_S = CS%MStar_coef*sqrt(max(0.0,Buoyancy_Flux) / UStar**2 / &
                    (Abs_Coriolis + 1.e-10*US%T_to_s) )
      ! The limit for rotation (Ekman length) limited mixing
      MStar_N =  CS%C_Ek * log( max( 1., UStar / (Abs_Coriolis + 1.e-10*US%T_to_s) / BLD ) )
    else
      ! The limit for the balance of rotation and stabilizing is f(L_Ekman,L_Obukhov)
      MStar_S = CS%MSTAR_COEF*sqrt(max(0.0, Buoyancy_Flux) / (UStar**2 * max(Abs_Coriolis, 1.e-20*US%T_to_s)))
      ! The limit for rotation (Ekman length) limited mixing
      MStar_N = 0.0
      if (UStar > Abs_Coriolis * BLD) Mstar_N = CS%C_EK * log(UStar / (Abs_Coriolis * BLD))
    endif

    ! Here 1.25 is about .5/von Karman, which gives the Obukhov limit.
    MStar = max(MStar_S, min(1.25, MStar_N))
    if (CS%MStar_Cap > 0.0) MStar = min( CS%MStar_Cap,MStar )
  elseif ( CS%mstar_scheme == MStar_from_RH18 ) then
    if (CS%answers_2018) then
      MStar_N = CS%RH18_MStar_cn1 * ( 1.0 - 1.0 / ( 1. + CS%RH18_MStar_cn2 * &
                exp( CS%RH18_mstar_CN3 * BLD * Abs_Coriolis / UStar) ) )
    else
      MSN_term = CS%RH18_MStar_cn2 * exp( CS%RH18_mstar_CN3 * BLD * Abs_Coriolis / UStar)
      MStar_N = (CS%RH18_MStar_cn1 *  MSN_term) / ( 1. + MSN_term)
    endif
    MStar_S = CS%RH18_MStar_CS1 * ( max(0.0, Buoyancy_Flux)**2 * BLD / &
             ( UStar**5 * max(Abs_Coriolis,1.e-20*US%T_to_s) ) )**CS%RH18_mstar_cs2
    MStar = MStar_N + MStar_S
  endif

  !/ 2. Adjust mstar to account for convective turbulence
  if (CS%answers_2018) then
    MStar_Conv_Red = 1. - CS%MStar_Convect_coef * (-min(0.0,Buoyancy_Flux) + 1.e-10*US%T_to_s**3*US%m_to_Z**2) / &
                         ( (-min(0.0,Buoyancy_Flux) + 1.e-10*US%T_to_s**3*US%m_to_Z**2) + &
                         2.0 *MStar * UStar**3 / BLD )
  else
    MSCR_term1 = -BLD * min(0.0, Buoyancy_Flux)
    MSCR_term2 = 2.0*MStar * UStar**3
    if ( abs(MSCR_term2) > 0.0) then
      MStar_Conv_Red = ((1.-CS%mstar_convect_coef) * MSCR_term1 + MSCR_term2) / (MSCR_term1 + MSCR_term2)
    else
      MStar_Conv_Red = 1.-CS%mstar_convect_coef
    endif
  endif

  !/3. Combine various mstar terms to get final value
  MStar = MStar * MStar_Conv_Red

  if (present(Langmuir_Number)) then
    !### In this call, ustar was previously ustar_mean.  Is this change deliberate?
    call mstar_Langmuir(CS, US, Abs_Coriolis, Buoyancy_Flux, UStar, BLD, Langmuir_Number, MStar, &
                        MStar_LT, Convect_Langmuir_Number)
  endif

end subroutine Find_Mstar

!> This subroutine modifies the Mstar value if the Langmuir number is present
subroutine Mstar_Langmuir(CS, US, Abs_Coriolis, Buoyancy_Flux, UStar, BLD, Langmuir_Number, &
                          Mstar, MStar_LT, Convect_Langmuir_Number)
  type(energetic_PBL_CS), pointer    :: CS    !< Energetic_PBL control structure.
  type(unit_scale_type), intent(in)  :: US    !< A dimensional unit scaling type
  real,                  intent(in)  :: Abs_Coriolis !< Absolute value of the Coriolis parameter [T-1 ~> s-1]
  real,                  intent(in)  :: Buoyancy_Flux !< Buoyancy flux [Z2 T-3 ~> m2 s-3]
  real,                  intent(in)  :: UStar !< Surface friction velocity with? gustiness [Z T-1 ~> m s-1]
  real,                  intent(in)  :: BLD   !< boundary layer depth [Z ~> m]
  real,                  intent(inout) :: Mstar !< Input/output mstar (Mixing/ustar**3) [nondim]
  real,                  intent(in)  :: Langmuir_Number !< Langmuir number [nondim]
  real,                  intent(out) :: MStar_LT !< Mstar increase due to Langmuir turbulence [nondim]
  real,                  intent(out) :: Convect_Langmuir_number !< Langmuir number including buoyancy flux [nondim]

  !/
  real, parameter :: Max_ratio = 1.0e16  ! The maximum value of a nondimensional ratio.
  real :: enhance_mstar ! A multiplicative scaling of mstar due to Langmuir turbulence.
  real :: mstar_LT_add ! A value that is added to mstar due to Langmuir turbulence.
  real :: iL_Ekman    ! Inverse of Ekman length scale [Z-1 ~> m-1].
  real :: iL_Obukhov  ! Inverse of Obukhov length scale [Z-1 ~> m-1].
  real :: I_ustar     ! The Adcroft reciprocal of ustar [T Z-1 ~> s m-1]
  real :: I_f         ! The Adcroft reciprocal of the Coriolis parameter [T ~> s]
  real :: MLD_Ekman          ! The ratio of the mixed layer depth to the Ekman layer depth [nondim].
  real :: Ekman_Obukhov      ! The Ekman layer thickness divided by the Obukhov depth [nondim].
  real :: MLD_Obukhov        ! The mixed layer depth divided by the Obukhov depth [nondim].
  real :: MLD_Obukhov_stab   ! Ratios of length scales where MLD is boundary layer depth [nondim].
  real :: Ekman_Obukhov_stab ! >
  real :: MLD_Obukhov_un     ! Ratios of length scales where MLD is boundary layer depth
  real :: Ekman_Obukhov_un   ! >

  ! Set default values for no Langmuir effects.
  enhance_mstar = 1.0 ; mstar_LT_add = 0.0

  if (CS%LT_Enhance_Form /= No_Langmuir) then
    ! a. Get parameters for modified LA
    if (CS%answers_2018) then
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
                    ( (1.0 + max(-0.5, CS%LaC_MLDoEK * MLD_Ekman)) + &
                   ((CS%LaC_EKoOB_stab * Ekman_Obukhov_stab + CS%LaC_EKoOB_un * Ekman_Obukhov_un) + &
                    (CS%LaC_MLDoOB_stab * MLD_Obukhov_stab  + CS%LaC_MLDoOB_un * MLD_Obukhov_un)) )

    if (CS%LT_Enhance_Form == Langmuir_rescale) then
      ! Enhancement is multiplied (added mst_lt set to 0)
      Enhance_mstar = min(CS%Max_Enhance_M, &
                          (1. + CS%LT_ENHANCE_COEF * Convect_Langmuir_Number**CS%LT_ENHANCE_EXP) )
    elseif (CS%LT_ENHANCE_Form == Langmuir_add) then
      ! or Enhancement is additive (multiplied enhance_m set to 1)
      mstar_LT_add = CS%LT_ENHANCE_COEF * Convect_Langmuir_Number**CS%LT_ENHANCE_EXP
    endif
  endif

  mstar_LT = (enhance_mstar - 1.0)*mstar + mstar_LT_add  ! Diagnose the full increase in mstar.
  mstar = mstar*enhance_mstar + mstar_LT_add

end subroutine Mstar_Langmuir


!> Copies the ePBL active mixed layer depth into MLD
subroutine energetic_PBL_get_MLD(CS, MLD, G, US, m_to_MLD_units)
  type(energetic_PBL_CS),           pointer     :: CS  !< Control structure for ePBL
  type(ocean_grid_type),            intent(in)  :: G   !< Grid structure
  type(unit_scale_type),            intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: MLD !< Depth of ePBL active mixing layer [m or other units]
  real,                   optional, intent(in)  :: m_to_MLD_units !< A conversion factor to the
                                                       !! desired units for MLD
  ! Local variables
  real :: scale  ! A dimensional rescaling factor
  integer :: i,j

  scale = US%Z_to_m ; if (present(m_to_MLD_units)) scale = scale * m_to_MLD_units

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    MLD(i,j) = scale*CS%ML_Depth(i,j)
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
  type(energetic_PBL_CS),  pointer       :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module
  ! Local variables
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_energetic_PBL"  ! This module's name.
  character(len=20)  :: tmpstr
  real :: omega_frac_dflt
  real :: R_Z3_T3_to_kg_s3 ! A conversion factor for work diagnostics [kg T3 R-1 Z-3 s-3 ~> nondim]
  integer :: isd, ied, jsd, jed
  integer :: mstar_mode, LT_enhance, wT_mode
  logical :: default_2018_answers
  logical :: use_temperature, use_omega
  logical :: use_la_windsea
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (associated(CS)) then
    call MOM_error(WARNING, "mixedlayer_init called with an associated"//&
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag
  CS%Time => Time

! Set default, read and log parameters
  call log_version(param_file, mdl, version, "")


!/1. General ePBL settings
  call get_param(param_file, mdl, "OMEGA", CS%omega, &
                 "The rotation rate of the earth.", units="s-1", &
                 default=7.2921e-5, scale=US%T_to_S)
  call get_param(param_file, mdl, "ML_USE_OMEGA", use_omega, &
                 "If true, use the absolute rotation rate instead of the "//&
                 "vertical component of rotation when setting the decay "//&
                 "scale for turbulence.", default=.false., do_not_log=.true.)
  omega_frac_dflt = 0.0
  if (use_omega) then
    call MOM_error(WARNING, "ML_USE_OMEGA is depricated; use ML_OMEGA_FRAC=1.0 instead.")
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
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.true.)
  call get_param(param_file, mdl, "EPBL_2018_ANSWERS", CS%answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)


  call get_param(param_file, mdl, "EPBL_ORIGINAL_PE_CALC", CS%orig_PE_calc, &
                 "If true, the ePBL code uses the original form of the "//&
                 "potential energy change code.  Otherwise, the newer "//&
                 "version that can work with successive increments to the "//&
                 "diffusivity in upward or downward passes is used.", default=.true.)

  call get_param(param_file, mdl, "MKE_TO_TKE_EFFIC", CS%MKE_to_TKE_effic, &
                 "The efficiency with which mean kinetic energy released "//&
                 "by mechanically forced entrainment of the mixed layer "//&
                 "is converted to turbulent kinetic energy.", units="nondim", &
                 default=0.0)
  call get_param(param_file, mdl, "TKE_DECAY", CS%TKE_decay, &
                 "TKE_DECAY relates the vertical rate of decay of the "//&
                 "TKE available for mechanical entrainment to the natural "//&
                 "Ekman depth.", units="nondim", default=2.5)


!/2. Options related to setting MSTAR
  call get_param(param_file, mdl, "EPBL_MSTAR_SCHEME", tmpstr, &
                 "EPBL_MSTAR_SCHEME selects the method for setting mstar.  Valid values are: \n"//&
                 "\t CONSTANT   - Use a fixed mstar given by MSTAR \n"//&
                 "\t OM4        - Use L_Ekman/L_Obukhov in the sabilizing limit, as in OM4 \n"//&
                 "\t REICHL_H18 - Use the scheme documented in Reichl & Hallberg, 2018.", &
                 default=CONSTANT_STRING, do_not_log=.true.)
  call get_param(param_file, mdl, "MSTAR_MODE", mstar_mode, default=-1)
  if (mstar_mode == 0) then
    tmpstr = CONSTANT_STRING
    call MOM_error(WARNING, "Use EPBL_MSTAR_SCHEME = CONSTANT instead of the archaic MSTAR_MODE = 0.")
  elseif (mstar_mode == 1) then
    call MOM_error(FATAL, "You are using a legacy mstar mode in ePBL that has been phased out. "//&
                          "If you need to use this setting please report this error.  Also use "//&
                          "EPBL_MSTAR_SCHEME to specify the scheme for mstar.")
  elseif (mstar_mode == 2) then
    tmpstr = OM4_STRING
    call MOM_error(WARNING, "Use EPBL_MSTAR_SCHEME = OM4 instead of the archaic MSTAR_MODE = 2.")
  elseif (mstar_mode == 3) then
    tmpstr = RH18_STRING
    call MOM_error(WARNING, "Use EPBL_MSTAR_SCHEME = REICHL_H18 instead of the archaic MSTAR_MODE = 3.")
  elseif (mstar_mode > 3) then
    call MOM_error(FATAL, "An unrecognized value of the obsolete parameter MSTAR_MODE was specified.")
  endif
  call log_param(param_file, mdl, "EPBL_MSTAR_SCHEME", tmpstr, &
                 "EPBL_MSTAR_SCHEME selects the method for setting mstar.  Valid values are: \n"//&
                 "\t CONSTANT   - Use a fixed mstar given by MSTAR \n"//&
                 "\t OM4        - Use L_Ekman/L_Obukhov in the sabilizing limit, as in OM4 \n"//&
                 "\t REICHL_H18 - Use the scheme documented in Reichl & Hallberg, 2018.", &
                 default=CONSTANT_STRING)
  tmpstr = uppercase(tmpstr)
  select case (tmpstr)
    case (CONSTANT_STRING)
      CS%mstar_Scheme = Use_Fixed_MStar
    case (OM4_STRING)
      CS%mstar_Scheme = MStar_from_Ekman
    case (RH18_STRING)
      CS%mstar_Scheme = MStar_from_RH18
    case default
      call MOM_mesg('energetic_PBL_init: EPBL_MSTAR_SCHEME ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "energetic_PBL_init: Unrecognized setting "// &
            "EPBL_MSTAR_SCHEME = "//trim(tmpstr)//" found in input file.")
  end select

  call get_param(param_file, mdl, "MSTAR", CS%fixed_mstar, &
                 "The ratio of the friction velocity cubed to the TKE input to the "//&
                 "mixed layer.  This option is used if EPBL_MSTAR_SCHEME = CONSTANT.", &
                 units="nondim", default=1.2, do_not_log=(CS%mstar_scheme/=Use_Fixed_MStar))
  call get_param(param_file, mdl, "MSTAR_CAP", CS%mstar_cap, &
                 "If this value is positive, it sets the maximum value of mstar "//&
                 "allowed in ePBL.  (This is not used if EPBL_MSTAR_SCHEME = CONSTANT).", &
                 units="nondim", default=-1.0, do_not_log=(CS%mstar_scheme==Use_Fixed_MStar))
  ! mstar_scheme==MStar_from_Ekman options
  call get_param(param_file, mdl, "MSTAR2_COEF1", CS%MSTAR_COEF, &
                 "Coefficient in computing mstar when rotation and stabilizing "//&
                 "effects are both important (used if EPBL_MSTAR_SCHEME = OM4).", &
                 units="nondim", default=0.3, do_not_log=(CS%mstar_scheme/=MStar_from_Ekman))
  call get_param(param_file, mdl, "MSTAR2_COEF2", CS%C_EK, &
                 "Coefficient in computing mstar when only rotation limits "// &
                 "the total mixing (used if EPBL_MSTAR_SCHEME = OM4)", &
                 units="nondim", default=0.085, do_not_log=(CS%mstar_scheme/=MStar_from_Ekman))
  ! mstar_scheme==MStar_from_RH18 options
  call get_param(param_file, mdl, "RH18_MSTAR_CN1", CS%RH18_mstar_cn1,&
                 "MSTAR_N coefficient 1 (outter-most coefficient for fit). "//&
                 "The value of 0.275 is given in RH18.  Increasing this "//&
                 "coefficient increases MSTAR for all values of Hf/ust, but more "//&
                 "effectively at low values (weakly developed OSBLs).", &
                 units="nondim", default=0.275, do_not_log=(CS%mstar_scheme/=MStar_from_RH18))
  call get_param(param_file, mdl, "RH18_MSTAR_CN2", CS%RH18_mstar_cn2,&
                 "MSTAR_N coefficient 2 (coefficient outside of exponential decay). "//&
                 "The value of 8.0 is given in RH18.  Increasing this coefficient "//&
                 "increases MSTAR for all values of HF/ust, with a much more even "//&
                 "effect across a wide range of Hf/ust than CN1.", &
                 units="nondim", default=8.0, do_not_log=(CS%mstar_scheme/=MStar_from_RH18))
  call get_param(param_file, mdl, "RH18_MSTAR_CN3", CS%RH18_mstar_CN3,&
                 "MSTAR_N coefficient 3 (exponential decay coefficient). "//&
                 "The value of -5.0 is given in RH18.  Increasing this increases how "//&
                 "quickly the value of MSTAR decreases as Hf/ust increases.", &
                  units="nondim", default=-5.0, do_not_log=(CS%mstar_scheme/=MStar_from_RH18))
  call get_param(param_file, mdl, "RH18_MSTAR_CS1", CS%RH18_mstar_cs1,&
                 "MSTAR_S coefficient for RH18 in stabilizing limit. "//&
                 "The value of 0.2 is given in RH18 and increasing it increases "//&
                 "MSTAR in the presence of a stabilizing surface buoyancy flux.", &
                 units="nondim", default=0.2, do_not_log=(CS%mstar_scheme/=MStar_from_RH18))
  call get_param(param_file, mdl, "RH18_MSTAR_CS2", CS%RH18_mstar_cs2,&
                 "MSTAR_S exponent for RH18 in stabilizing limit. "//&
                 "The value of 0.4 is given in RH18 and increasing it increases MSTAR "//&
                 "exponentially in the presence of a stabilizing surface buoyancy flux.", &
                 Units="nondim", default=0.4, do_not_log=(CS%mstar_scheme/=MStar_from_RH18))


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
  !### THIS DEFAULT SHOULD BECOME TRUE.
  call get_param(param_file, mdl, "USE_MLD_ITERATION", CS%Use_MLD_iteration, &
                 "A logical that specifies whether or not to use the "//&
                 "distance to the bottom of the actively turbulent boundary "//&
                 "layer to help set the EPBL length scale.", default=.false.)
  call get_param(param_file, mdl, "EPBL_TRANSITION_SCALE", CS%transLay_scale, &
                 "A scale for the mixing length in the transition layer "//&
                 "at the edge of the boundary layer as a fraction of the "//&
                 "boundary layer thickness.", units="nondim", default=0.1)
  if ( CS%Use_MLD_iteration .and. abs(CS%transLay_scale-0.5) >= 0.5) then
    call MOM_error(FATAL, "If flag USE_MLD_ITERATION is true, then "//&
                 "EPBL_TRANSITION should be greater than 0 and less than 1.")
  endif

  call get_param(param_file, mdl, "MLD_ITERATION_GUESS", CS%MLD_ITERATION_GUESS, &
                 "A logical that specifies whether or not to use the "//&
                 "previous timestep MLD as a first guess in the MLD iteration. "//&
                 "The default is false to facilitate reproducibility.", default=.false.)
  call get_param(param_file, mdl, "EPBL_MLD_TOLERANCE", CS%MLD_tol, &
                 "The tolerance for the iteratively determined mixed "//&
                 "layer depth.  This is only used with USE_MLD_ITERATION.", &
                 units="meter", default=1.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "EPBL_MLD_MAX_ITS", CS%max_MLD_its, &
                 "The maximum number of iterations that can be used to find a self-consistent "//&
                 "mixed layer depth.  For now, due to the use of bisection, the maximum number "//&
                 "iteractions needed is set by Depth/2^MAX_ITS < EPBL_MLD_TOLERANCE.", &
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
  call get_param(param_file, mdl, "EPBL_VEL_SCALE_SCHEME", tmpstr, &
                 "Selects the method for translating TKE into turbulent velocities. "//&
                 "Valid values are: \n"//&
                 "\t CUBE_ROOT_TKE  - A constant times the cube root of remaining TKE. \n"//&
                 "\t REICHL_H18 - Use the scheme based on a combination of w* and v* as \n"//&
                 "\t              documented in Reichl & Hallberg, 2018.", &
                 default=ROOT_TKE_STRING, do_not_log=.true.)
  call get_param(param_file, mdl, "EPBL_VEL_SCALE_MODE", wT_mode, default=-1)
  if (wT_mode == 0) then
    tmpstr = ROOT_TKE_STRING
    call MOM_error(WARNING, "Use EPBL_VEL_SCALE_SCHEME = CUBE_ROOT_TKE instead of the archaic EPBL_VEL_SCALE_MODE = 0.")
  elseif (wT_mode == 1) then
    tmpstr = RH18_STRING
    call MOM_error(WARNING, "Use EPBL_VEL_SCALE_SCHEME = REICHL_H18 instead of the archaic EPBL_VEL_SCALE_MODE = 1.")
  elseif (wT_mode >= 2) then
    call MOM_error(FATAL, "An unrecognized value of the obsolete parameter EPBL_VEL_SCALE_MODE was specified.")
  endif
  call log_param(param_file, mdl, "EPBL_VEL_SCALE_SCHEME", tmpstr, &
                 "Selects the method for translating TKE into turbulent velocities. "//&
                 "Valid values are: \n"//&
                 "\t CUBE_ROOT_TKE  - A constant times the cube root of remaining TKE. \n"//&
                 "\t REICHL_H18 - Use the scheme based on a combination of w* and v* as \n"//&
                 "\t              documented in Reichl & Hallberg, 2018.", &
                 default=ROOT_TKE_STRING)
  tmpstr = uppercase(tmpstr)
  select case (tmpstr)
    case (ROOT_TKE_STRING)
      CS%wT_scheme = wT_from_cRoot_TKE
    case (RH18_STRING)
      CS%wT_scheme = wT_from_RH18
    case default
      call MOM_mesg('energetic_PBL_init: EPBL_VEL_SCALE_SCHEME ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "energetic_PBL_init: Unrecognized setting "// &
            "EPBL_VEL_SCALE_SCHEME = "//trim(tmpstr)//" found in input file.")
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

  !/ Options related to Langmuir turbulence
  call get_param(param_file, mdl, "USE_LA_LI2016", use_LA_Windsea, &
       "A logical to use the Li et al. 2016 (submitted) formula to "//&
       "determine the Langmuir number.", units="nondim", default=.false.)
  ! Note this can be activated in other ways, but this preserves the old method.
  if (use_LA_windsea) then
    CS%USE_LT = .true.
  else
    call get_param(param_file, mdl, "EPBL_LT", CS%USE_LT, &
                 "A logical to use a LT parameterization.", &
                 units="nondim", default=.false.)
  endif
  if (CS%USE_LT) then
    call get_param(param_file, mdl, "EPBL_LANGMUIR_SCHEME", tmpstr, &
                 "EPBL_LANGMUIR_SCHEME selects the method for including Langmuir turbulence. "//&
                 "Valid values are: \n"//&
                 "\t NONE     - Do not do any extra mixing due to Langmuir turbulence \n"//&
                 "\t RESCALE  - Use a multiplicative rescaling of mstar to account for Langmuir turbulence \n"//&
                 "\t ADDITIVE - Add a Langmuir turblence contribution to mstar to other contributions", &
                 default=NONE_STRING, do_not_log=.true.)
    call get_param(param_file, mdl, "LT_ENHANCE", LT_enhance, default=-1)
    if (LT_ENHANCE == 0) then
      tmpstr = NONE_STRING
      call MOM_error(WARNING, "Use EPBL_LANGMUIR_SCHEME = NONE instead of the archaic LT_ENHANCE = 0.")
    elseif (LT_ENHANCE == 1) then
      call MOM_error(FATAL, "You are using a legacy LT_ENHANCE mode in ePBL that has been phased out. "//&
                            "If you need to use this setting please report this error.  Also use "//&
                            "EPBL_LANGMUIR_SCHEME to specify the scheme for mstar.")
    elseif (LT_ENHANCE == 2) then
      tmpstr = RESCALED_STRING
      call MOM_error(WARNING, "Use EPBL_LANGMUIR_SCHEME = RESCALE instead of the archaic LT_ENHANCE = 2.")
    elseif (LT_ENHANCE == 3) then
      tmpstr = ADDITIVE_STRING
      call MOM_error(WARNING, "Use EPBL_LANGMUIR_SCHEME = ADDITIVE instead of the archaic LT_ENHANCE = 3.")
    elseif (LT_ENHANCE > 3) then
      call MOM_error(FATAL, "An unrecognized value of the obsolete parameter LT_ENHANCE was specified.")
    endif
    call log_param(param_file, mdl, "EPBL_LANGMUIR_SCHEME", tmpstr, &
                 "EPBL_LANGMUIR_SCHEME selects the method for including Langmuir turbulence. "//&
                 "Valid values are: \n"//&
                 "\t NONE     - Do not do any extra mixing due to Langmuir turbulence \n"//&
                 "\t RESCALE  - Use a multiplicative rescaling of mstar to account for Langmuir turbulence \n"//&
                 "\t ADDITIVE - Add a Langmuir turblence contribution to mstar to other contributions", &
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

    call get_param(param_file, mdl, "LT_ENHANCE_COEF", CS%LT_ENHANCE_COEF, &
                 "Coefficient for Langmuir enhancement of mstar", &
                 units="nondim", default=0.447, do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_ENHANCE_EXP", CS%LT_ENHANCE_EXP, &
                 "Exponent for Langmuir enhancementt of mstar", &
                 units="nondim", default=-1.33,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC1", CS%LaC_MLDoEK, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "MLD approaching Ekman depth.", &
                 units="nondim", default=-0.87,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC2", CS%LaC_MLDoOB_stab, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "MLD approaching stable Obukhov depth.", &
                 units="nondim", default=0.0,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC3", CS%LaC_MLDoOB_un, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "MLD approaching unstable Obukhov depth.", &
                 units="nondim", default=0.0,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC4", CS%Lac_EKoOB_stab, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "ratio of Ekman to stable Obukhov depth.", &
                 units="nondim", default=0.95,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC5", CS%Lac_EKoOB_un, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "ratio of Ekman to unstable Obukhov depth.", &
                 units="nondim", default=0.95,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
  endif


!/ Logging parameters
  ! This gives a minimum decay scale that is typically much less than Angstrom.
  CS%ustar_min = 2e-4*CS%omega*(GV%Angstrom_Z + GV%H_to_Z*GV%H_subroundoff)
  call log_param(param_file, mdl, "EPBL_USTAR_MIN", CS%ustar_min*US%Z_to_m*US%s_to_T, &
                 "The (tiny) minimum friction velocity used within the "//&
                 "ePBL code, derived from OMEGA and ANGSTROM.", units="m s-1")


!/ Checking output flags
  R_Z3_T3_to_kg_s3 = US%R_to_kg_m3 * US%Z_to_m**3  * US%s_to_T**3
  CS%id_ML_depth = register_diag_field('ocean_model', 'ePBL_h_ML', diag%axesT1, &
      Time, 'Surface boundary layer depth', 'm', conversion=US%Z_to_m, &
      cmor_long_name='Ocean Mixed Layer Thickness Defined by Mixing Scheme')
  CS%id_TKE_wind = register_diag_field('ocean_model', 'ePBL_TKE_wind', diag%axesT1, &
      Time, 'Wind-stirring source of mixed layer TKE', 'm3 s-3', conversion=R_Z3_T3_to_kg_s3)
  CS%id_TKE_MKE = register_diag_field('ocean_model', 'ePBL_TKE_MKE', diag%axesT1, &
      Time, 'Mean kinetic energy source of mixed layer TKE', 'm3 s-3', conversion=R_Z3_T3_to_kg_s3)
  CS%id_TKE_conv = register_diag_field('ocean_model', 'ePBL_TKE_conv', diag%axesT1, &
      Time, 'Convective source of mixed layer TKE', 'm3 s-3', conversion=R_Z3_T3_to_kg_s3)
  CS%id_TKE_forcing = register_diag_field('ocean_model', 'ePBL_TKE_forcing', diag%axesT1, &
      Time, 'TKE consumed by mixing surface forcing or penetrative shortwave radation '//&
            'through model layers', 'm3 s-3', conversion=R_Z3_T3_to_kg_s3)
  CS%id_TKE_mixing = register_diag_field('ocean_model', 'ePBL_TKE_mixing', diag%axesT1, &
      Time, 'TKE consumed by mixing that deepens the mixed layer', 'm3 s-3', conversion=R_Z3_T3_to_kg_s3)
  CS%id_TKE_mech_decay = register_diag_field('ocean_model', 'ePBL_TKE_mech_decay', diag%axesT1, &
      Time, 'Mechanical energy decay sink of mixed layer TKE', 'm3 s-3', conversion=R_Z3_T3_to_kg_s3)
  CS%id_TKE_conv_decay = register_diag_field('ocean_model', 'ePBL_TKE_conv_decay', diag%axesT1, &
      Time, 'Convective energy decay sink of mixed layer TKE', 'm3 s-3', conversion=R_Z3_T3_to_kg_s3)
  CS%id_Mixing_Length = register_diag_field('ocean_model', 'Mixing_Length', diag%axesTi, &
      Time, 'Mixing Length that is used', 'm', conversion=US%Z_to_m)
  CS%id_Velocity_Scale = register_diag_field('ocean_model', 'Velocity_Scale', diag%axesTi, &
      Time, 'Velocity Scale that is used.', 'm s-1', conversion=US%Z_to_m*US%s_to_T)
  CS%id_MSTAR_mix = register_diag_field('ocean_model', 'MSTAR', diag%axesT1, &
      Time, 'Total mstar that is used.', 'nondim')

  if (CS%use_LT) then
    CS%id_LA = register_diag_field('ocean_model', 'LA', diag%axesT1, &
        Time, 'Langmuir number.', 'nondim')
    CS%id_LA_mod = register_diag_field('ocean_model', 'LA_MOD', diag%axesT1, &
        Time, 'Modified Langmuir number.', 'nondim')
    CS%id_MSTAR_LT = register_diag_field('ocean_model', 'MSTAR_LT', diag%axesT1, &
        Time, 'Increase in mstar due to Langmuir Turbulence.', 'nondim')
  endif

  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", use_temperature, &
                 "If true, temperature and salinity are used as state "//&
                 "variables.", default=.true.)

  if (max(CS%id_TKE_wind, CS%id_TKE_MKE, CS%id_TKE_conv, &
          CS%id_TKE_mixing, CS%id_TKE_mech_decay, CS%id_TKE_forcing, &
          CS%id_TKE_conv_decay) > 0) then
    call safe_alloc_alloc(CS%diag_TKE_wind, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%diag_TKE_MKE, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%diag_TKE_conv, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%diag_TKE_forcing, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%diag_TKE_mixing, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%diag_TKE_mech_decay, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%diag_TKE_conv_decay, isd, ied, jsd, jed)

    CS%TKE_diagnostics = .true.
  endif
  if (CS%id_Velocity_Scale>0) call safe_alloc_alloc(CS%Velocity_Scale, isd, ied, jsd, jed, GV%ke+1)
  if (CS%id_Mixing_Length>0) call safe_alloc_alloc(CS%Mixing_Length, isd, ied, jsd, jed, GV%ke+1)

  call safe_alloc_alloc(CS%ML_depth, isd, ied, jsd, jed)
  if (max(CS%id_mstar_mix, CS%id_LA, CS%id_LA_mod, CS%id_MSTAR_LT ) >0) then
    call safe_alloc_alloc(CS%Mstar_mix, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%LA, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%LA_MOD, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%MSTAR_LT, isd, ied, jsd, jed)
  endif

end subroutine energetic_PBL_init

!> Clean up and deallocate memory associated with the energetic_PBL module.
subroutine energetic_PBL_end(CS)
  type(energetic_PBL_CS), pointer :: CS !< Energetic_PBL control structure that
                                        !! will be deallocated in this subroutine.

  if (.not.associated(CS)) return

  if (allocated(CS%ML_depth))            deallocate(CS%ML_depth)
  if (allocated(CS%LA))                  deallocate(CS%LA)
  if (allocated(CS%LA_MOD))              deallocate(CS%LA_MOD)
  if (allocated(CS%MSTAR_MIX))           deallocate(CS%MSTAR_MIX)
  if (allocated(CS%MSTAR_LT))            deallocate(CS%MSTAR_LT)
  if (allocated(CS%diag_TKE_wind))       deallocate(CS%diag_TKE_wind)
  if (allocated(CS%diag_TKE_MKE))        deallocate(CS%diag_TKE_MKE)
  if (allocated(CS%diag_TKE_conv))       deallocate(CS%diag_TKE_conv)
  if (allocated(CS%diag_TKE_forcing))    deallocate(CS%diag_TKE_forcing)
  if (allocated(CS%diag_TKE_mixing))     deallocate(CS%diag_TKE_mixing)
  if (allocated(CS%diag_TKE_mech_decay)) deallocate(CS%diag_TKE_mech_decay)
  if (allocated(CS%diag_TKE_conv_decay)) deallocate(CS%diag_TKE_conv_decay)
  if (allocated(CS%Mixing_Length))       deallocate(CS%Mixing_Length)
  if (allocated(CS%Velocity_Scale))      deallocate(CS%Velocity_Scale)

  deallocate(CS)

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
!! review paper by Niiler and Kraus (1979). The work here draws in
!! with particular on the form for TKE decay proposed by Oberhuber
!! (JPO, 1993, 808-829), with an extension to a refined bulk mixed
!! layer as described in Hallberg (Aha Huliko'a, 2003).  The physical
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
