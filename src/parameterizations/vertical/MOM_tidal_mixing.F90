!> Interface to vertical tidal mixing schemes including CVMix tidal mixing.
module MOM_tidal_mixing

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator,      only : diag_ctrl, time_type, register_diag_field
use MOM_diag_mediator,      only : safe_alloc_ptr, post_data
use MOM_debugging,          only : hchksum
use MOM_EOS,                only : calculate_density
use MOM_error_handler,      only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,        only : openParameterBlock, closeParameterBlock
use MOM_file_parser,        only : get_param, log_param, log_version, param_file_type
use MOM_grid,               only : ocean_grid_type
use MOM_io,                 only : slasher, MOM_read_data, field_size
use MOM_remapping,          only : remapping_CS, initialize_remapping, remapping_core_h
use MOM_string_functions,   only : uppercase, lowercase
use MOM_unit_scaling,       only : unit_scale_type
use MOM_variables,          only : thermo_var_ptrs, p3d
use MOM_verticalGrid,       only : verticalGrid_type
use CVMix_tidal,            only : CVMix_init_tidal, CVMix_compute_Simmons_invariant
use CVMix_tidal,            only : CVMix_coeffs_tidal, CVMix_tidal_params_type
use CVMix_tidal,            only : CVMix_compute_Schmittner_invariant, CVMix_compute_SchmittnerCoeff
use CVMix_tidal,            only : CVMix_coeffs_tidal_schmittner
use CVMix_kinds_and_types,  only : CVMix_global_params_type
use CVMix_put_get,          only : CVMix_put

implicit none ; private

#include <MOM_memory.h>

public tidal_mixing_init
public setup_tidal_diagnostics
public calculate_tidal_mixing
public post_tidal_diagnostics
public tidal_mixing_h_amp
public tidal_mixing_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Containers for tidal mixing diagnostics
type, public :: tidal_mixing_diags ; private
  real, pointer, dimension(:,:,:) :: &
    Kd_itidal             => NULL(),& !< internal tide diffusivity at interfaces [Z2 T-1 ~> m2 s-1].
    Fl_itidal             => NULL(),& !< vertical flux of tidal turbulent dissipation [Z3 T-3 ~> m3 s-3]
    Kd_Niku               => NULL(),& !< lee-wave diffusivity at interfaces [Z2 T-1 ~> m2 s-1].
    Kd_Niku_work          => NULL(),& !< layer integrated work by lee-wave driven mixing [R Z3 T-3 ~> W m-2]
    Kd_Itidal_Work        => NULL(),& !< layer integrated work by int tide driven mixing [R Z3 T-3 ~> W m-2]
    Kd_Lowmode_Work       => NULL(),& !< layer integrated work by low mode driven mixing [R Z3 T-3 ~> W m-2]
    N2_int                => NULL(),& !< Bouyancy frequency squared at interfaces [T-2 ~> s-2]
    vert_dep_3d           => NULL(),& !< The 3-d mixing energy deposition [W m-3]
    Schmittner_coeff_3d   => NULL()   !< The coefficient in the Schmittner et al mixing scheme, in UNITS?
  real, pointer, dimension(:,:,:) :: tidal_qe_md => NULL() !< Input tidal energy dissipated locally,
                                           !! interpolated to model vertical coordinate [W m-3?]
  real, pointer, dimension(:,:,:) :: Kd_lowmode => NULL() !< internal tide diffusivity at interfaces
                                           !! due to propagating low modes [Z2 T-1 ~> m2 s-1].
  real, pointer, dimension(:,:,:) :: Fl_lowmode => NULL() !< vertical flux of tidal turbulent
                                           !! dissipation due to propagating low modes [Z3 T-3 ~> m3 s-3]
  real, pointer, dimension(:,:) :: &
    TKE_itidal_used           => NULL(),& !< internal tide TKE input at ocean bottom [R Z3 T-3 ~> W m-2]
    N2_bot                    => NULL(),& !< bottom squared buoyancy frequency [T-2 ~> s-2]
    N2_meanz                  => NULL(),& !< vertically averaged buoyancy frequency [T-2 ~> s-2]
    Polzin_decay_scale_scaled => NULL(),& !< vertical scale of decay for tidal dissipation [Z ~> m]
    Polzin_decay_scale        => NULL(),& !< vertical decay scale for tidal diss with Polzin [Z ~> m]
    Simmons_coeff_2d          => NULL()   !< The Simmons et al mixing coefficient

end type

!> Control structure with parameters for the tidal mixing module.
type, public :: tidal_mixing_cs ; private
  logical :: debug = .true.   !< If true, do more extensive debugging checks.  This is hard-coded.

  ! Parameters
  logical :: int_tide_dissipation = .false. !< Internal tide conversion (from barotropic)
                              !! with the schemes of St Laurent et al (2002) & Simmons et al (2004)

  integer :: Int_tide_profile !< A coded integer indicating the vertical profile
                              !! for dissipation of the internal waves.  Schemes that are
                              !! currently encoded are St Laurent et al (2002) and Polzin (2009).
  logical :: Lee_wave_dissipation = .false. !< Enable lee-wave driven mixing, following
                              !! Nikurashin (2010), with a vertical energy
                              !! deposition profile specified by Lee_wave_profile to be
                              !! St Laurent et al (2002) or Simmons et al (2004) scheme

  integer :: Lee_wave_profile !< A coded integer indicating the vertical profile
                              !! for dissipation of the lee waves.  Schemes that are
                              !! currently encoded are St Laurent et al (2002) and
                              !! Polzin (2009).
  real :: Int_tide_decay_scale !< decay scale for internal wave TKE [Z ~> m].

  real :: Mu_itides           !< efficiency for conversion of dissipation
                              !! to potential energy [nondim]

  real :: Gamma_itides        !< fraction of local dissipation [nondim]

  real :: Gamma_lee           !< fraction of local dissipation for lee waves
                              !! (Nikurashin's energy input) [nondim]
  real :: Decay_scale_factor_lee !< Scaling factor for the decay scale of lee
                              !! wave energy dissipation [nondim]

  real :: min_zbot_itides     !< minimum depth for internal tide conversion [Z ~> m].
  logical :: Lowmode_itidal_dissipation = .false.  !< If true, consider mixing due to breaking low
                              !! modes that have been remotely generated using an internal tidal
                              !! dissipation scheme to specify the vertical profile of the energy
                              !! input to drive diapycnal mixing, along the lines of St. Laurent
                              !! et al. (2002) and Simmons et al. (2004).

  real :: Nu_Polzin           !< The non-dimensional constant used in Polzin form of
                              !! the vertical scale of decay of tidal dissipation [nondim]

  real :: Nbotref_Polzin      !< Reference value for the buoyancy frequency at the
                              !! ocean bottom used in Polzin formulation of the
                              !! vertical scale of decay of tidal dissipation [T-1 ~> s-1]
  real :: Polzin_decay_scale_factor !< Scaling factor for the decay length scale
                              !! of the tidal dissipation profile in Polzin [nondim]
  real :: Polzin_decay_scale_max_factor !< The decay length scale of tidal dissipation
                              !! profile in Polzin formulation should not exceed
                              !! Polzin_decay_scale_max_factor * depth of the ocean [nondim].
  real :: Polzin_min_decay_scale !< minimum decay scale of the tidal dissipation
                              !! profile in Polzin formulation [Z ~> m].

  real :: TKE_itide_max       !< maximum internal tide conversion [R Z3 T-3 ~> W m-2]
                              !! available to mix above the BBL

  real :: utide               !< constant tidal amplitude [Z T-1 ~> m s-1] if READ_TIDEAMP is false.
  real :: kappa_itides        !< topographic wavenumber and non-dimensional scaling [Z-1 ~> m-1].
  real :: kappa_h2_factor     !< factor for the product of wavenumber * rms sgs height
  character(len=200) :: inputdir !< The directory in which to find input files

  logical :: use_CVMix_tidal = .false. !< true if CVMix is to be used for determining
                              !! diffusivity due to tidal mixing

  real :: min_thickness       !< Minimum thickness allowed [m]

  ! CVMix-specific parameters
  integer                         :: CVMix_tidal_scheme = -1  !< 1 for Simmons, 2 for Schmittner
  type(CVMix_tidal_params_type)   :: CVMix_tidal_params !< A CVMix-specific type with parameters for tidal mixing
  type(CVMix_global_params_type)  :: CVMix_glb_params   !< CVMix-specific for Prandtl number only
  real                            :: tidal_max_coef     !< CVMix-specific maximum allowable tidal diffusivity. [m^2/s]
  real                            :: tidal_diss_lim_tc  !< CVMix-specific dissipation limit depth for
                                                        !! tidal-energy-constituent data [Z ~> m].
  type(remapping_CS)              :: remap_CS           !< The control structure for remapping
  logical :: remap_answers_2018 = .true.  !< If true, use the order of arithmetic and expressions that
                                       !! recover the remapping answers from 2018.  If false, use more
                                       !! robust forms of the same remapping expressions.

  ! Data containers
  real, pointer, dimension(:,:) :: TKE_Niku    => NULL() !< Lee wave driven Turbulent Kinetic Energy input
                                                         !! [R Z3 T-3 ~> W m-2]
  real, pointer, dimension(:,:) :: TKE_itidal  => NULL() !< The internal Turbulent Kinetic Energy input divided
                                                         !! by the bottom stratfication [R Z3 T-2 ~> J m-2].
  real, pointer, dimension(:,:) :: Nb          => NULL() !< The near bottom buoyancy frequency [T-1 ~> s-1].
  real, pointer, dimension(:,:) :: mask_itidal => NULL() !< A mask of where internal tide energy is input
  real, pointer, dimension(:,:) :: h2          => NULL() !< Squared bottom depth variance [Z2 ~> m2].
  real, pointer, dimension(:,:) :: tideamp     => NULL() !< RMS tidal amplitude [Z T-1 ~> m s-1]
  real, allocatable, dimension(:)     :: h_src           !< tidal constituent input layer thickness [m]
  real, allocatable, dimension(:,:)   :: tidal_qe_2d     !< Tidal energy input times the local dissipation
                                                         !! fraction, q*E(x,y), with the CVMix implementation
                                                         !! of Jayne et al tidal mixing [W m-2].
                                                         !! TODO: make this E(x,y) only
  real, allocatable, dimension(:,:,:) :: tidal_qe_3d_in  !< q*E(x,y,z) with the Schmittner parameterization [W m-3?]

  logical :: answers_2018   !< If true, use the order of arithmetic and expressions that recover the
                            !! answers from the end of 2018.  Otherwise, use updated and more robust
                            !! forms of the same expressions.

  ! Diagnostics
  type(diag_ctrl),          pointer :: diag => NULL() !< structure to regulate diagnostic output timing
  type(tidal_mixing_diags), pointer :: dd => NULL() !< A pointer to a structure of diagnostic arrays

  !>@{ Diagnostic identifiers
  integer :: id_TKE_itidal                = -1
  integer :: id_TKE_leewave               = -1
  integer :: id_Kd_itidal                 = -1
  integer :: id_Kd_Niku                   = -1
  integer :: id_Kd_lowmode                = -1
  integer :: id_Kd_Itidal_Work            = -1
  integer :: id_Kd_Niku_Work              = -1
  integer :: id_Kd_Lowmode_Work           = -1
  integer :: id_Nb                        = -1
  integer :: id_N2_bot                    = -1
  integer :: id_N2_meanz                  = -1
  integer :: id_Fl_itidal                 = -1
  integer :: id_Fl_lowmode                = -1
  integer :: id_Polzin_decay_scale        = -1
  integer :: id_Polzin_decay_scale_scaled = -1
  integer :: id_N2_int                    = -1
  integer :: id_Simmons_coeff             = -1
  integer :: id_Schmittner_coeff          = -1
  integer :: id_tidal_qe_md               = -1
  integer :: id_vert_dep                  = -1
  !>@}

end type tidal_mixing_cs

!>@{ Coded parmameters for specifying mixing schemes
character*(20), parameter :: STLAURENT_PROFILE_STRING   = "STLAURENT_02"
character*(20), parameter :: POLZIN_PROFILE_STRING      = "POLZIN_09"
integer,        parameter :: STLAURENT_02 = 1
integer,        parameter :: POLZIN_09    = 2
character*(20), parameter :: SIMMONS_SCHEME_STRING      = "SIMMONS"
character*(20), parameter :: SCHMITTNER_SCHEME_STRING   = "SCHMITTNER"
integer,        parameter :: SIMMONS   = 1
integer,        parameter :: SCHMITTNER   = 2
!>@}

contains

!> Initializes internal tidal dissipation scheme for diapycnal mixing
logical function tidal_mixing_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type),          intent(in)    :: Time       !< The current time.
  type(ocean_grid_type),    intent(in)    :: G          !< Grid structure.
  type(verticalGrid_type),  intent(in)    :: GV         !< Vertical grid structure.
  type(unit_scale_type),    intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),    intent(in)    :: param_file !< Run-time parameter file handle
  type(diag_ctrl), target,  intent(inout) :: diag       !< Diagnostics control structure.
  type(tidal_mixing_cs),    pointer       :: CS         !< This module's control structure.

  ! Local variables
  logical :: read_tideamp
  logical :: default_2018_answers
  character(len=20)  :: tmpstr, int_tide_profile_str
  character(len=20)  :: CVMix_tidal_scheme_str, tidal_energy_type
  character(len=200) :: filename, h2_file, Niku_TKE_input_file
  character(len=200) :: tidal_energy_file, tideamp_file
  real :: utide, hamp, prandtl_tidal, max_frac_rough
  real :: Niku_scale ! local variable for scaling the Nikurashin TKE flux data
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_tidal_mixing"     !< This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "tidal_mixing_init called when control structure "// &
                            "is already associated.")
    return
  endif
  allocate(CS)
  allocate(CS%dd)

  CS%debug = CS%debug.and.is_root_pe()

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%diag => diag

  ! Read parameters
  call get_param(param_file, mdl, "USE_CVMix_TIDAL", CS%use_CVMix_tidal, &
                 default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "INT_TIDE_DISSIPATION", CS%int_tide_dissipation, &
                 default=CS%use_CVMix_tidal, do_not_log=.true.)
  call log_version(param_file, mdl, version, &
                 "Vertical Tidal Mixing Parameterization", &
                 all_default=.not.(CS%use_CVMix_tidal .or. CS%int_tide_dissipation))
  call get_param(param_file, mdl, "USE_CVMix_TIDAL", CS%use_CVMix_tidal, &
                 "If true, turns on tidal mixing via CVMix", &
                 default=.false.)

  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, default=".",do_not_log=.true.)
  CS%inputdir = slasher(CS%inputdir)
  call get_param(param_file, mdl, "INT_TIDE_DISSIPATION", CS%int_tide_dissipation, &
                 "If true, use an internal tidal dissipation scheme to "//&
                 "drive diapycnal mixing, along the lines of St. Laurent "//&
                 "et al. (2002) and Simmons et al. (2004).", default=CS%use_CVMix_tidal)

  ! return if tidal mixing is inactive
  tidal_mixing_init = CS%int_tide_dissipation
  if (.not. tidal_mixing_init) return

  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDAL_MIXING_2018_ANSWERS", CS%answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
  call get_param(param_file, mdl, "REMAPPING_2018_ANSWERS", CS%remap_answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)

  if (CS%int_tide_dissipation) then

    ! Read in CVMix tidal scheme if CVMix tidal mixing is on
    if (CS%use_CVMix_tidal) then
      call get_param(param_file, mdl, "CVMIX_TIDAL_SCHEME", CVMix_tidal_scheme_str, &
                 "CVMIX_TIDAL_SCHEME selects the CVMix tidal mixing "//&
                 "scheme with INT_TIDE_DISSIPATION. Valid values are:\n"//&
                 "\t SIMMONS - Use the Simmons et al (2004) tidal \n"//&
                 "\t                mixing scheme.\n"//&
                 "\t SCHMITTNER - Use the Schmittner et al (2014) tidal \n"//&
                 "\t                mixing scheme.", &
                 default=SIMMONS_SCHEME_STRING)
      CVMix_tidal_scheme_str = uppercase(CVMix_tidal_scheme_str)

      select case (CVMix_tidal_scheme_str)
        case (SIMMONS_SCHEME_STRING)    ; CS%CVMix_tidal_scheme = SIMMONS
        case (SCHMITTNER_SCHEME_STRING) ; CS%CVMix_tidal_scheme = SCHMITTNER
        case default
        call MOM_error(FATAL, "tidal_mixing_init: Unrecognized setting "// &
            "#define CVMIX_TIDAL_SCHEME "//trim(CVMix_tidal_scheme_str)//" found in input file.")
      end select
    endif ! CS%use_CVMix_tidal

    ! Read in vertical profile of tidal energy dissipation
    if ( CS%CVMix_tidal_scheme.eq.SCHMITTNER .or. .not. CS%use_CVMix_tidal) then
      call get_param(param_file, mdl, "INT_TIDE_PROFILE", int_tide_profile_str, &
                   "INT_TIDE_PROFILE selects the vertical profile of energy "//&
                   "dissipation with INT_TIDE_DISSIPATION. Valid values are:\n"//&
                   "\t STLAURENT_02 - Use the St. Laurent et al exponential \n"//&
                   "\t                decay profile.\n"//&
                   "\t POLZIN_09 - Use the Polzin WKB-stretched algebraic \n"//&
                   "\t                decay profile.", &
                   default=STLAURENT_PROFILE_STRING)
      int_tide_profile_str = uppercase(int_tide_profile_str)

      select case (int_tide_profile_str)
        case (STLAURENT_PROFILE_STRING)   ; CS%int_tide_profile = STLAURENT_02
        case (POLZIN_PROFILE_STRING)      ; CS%int_tide_profile = POLZIN_09
        case default
          call MOM_error(FATAL, "tidal_mixing_init: Unrecognized setting "// &
              "#define INT_TIDE_PROFILE "//trim(int_tide_profile_str)//" found in input file.")
      end select
    endif

  elseif (CS%use_CVMix_tidal) then
        call MOM_error(FATAL, "tidal_mixing_init: Cannot set INT_TIDE_DISSIPATION to False "// &
            "when USE_CVMix_TIDAL is set to True.")
  endif

  call get_param(param_file, mdl, "LEE_WAVE_DISSIPATION", CS%Lee_wave_dissipation, &
                 "If true, use an lee wave driven dissipation scheme to "//&
                 "drive diapycnal mixing, along the lines of Nikurashin "//&
                 "(2010) and using the St. Laurent et al. (2002) "//&
                 "and Simmons et al. (2004) vertical profile", default=.false.)
  if (CS%lee_wave_dissipation) then
    if (CS%use_CVMix_tidal) then
        call MOM_error(FATAL, "tidal_mixing_init: Lee wave driven dissipation scheme cannot "// &
            "be used when CVMix tidal mixing scheme is active.")
    endif
    call get_param(param_file, mdl, "LEE_WAVE_PROFILE", tmpstr, &
                 "LEE_WAVE_PROFILE selects the vertical profile of energy "//&
                 "dissipation with LEE_WAVE_DISSIPATION. Valid values are:\n"//&
                 "\t STLAURENT_02 - Use the St. Laurent et al exponential \n"//&
                 "\t                decay profile.\n"//&
                 "\t POLZIN_09 - Use the Polzin WKB-stretched algebraic \n"//&
                 "\t                decay profile.", &
                 default=STLAURENT_PROFILE_STRING)
    tmpstr = uppercase(tmpstr)
    select case (tmpstr)
      case (STLAURENT_PROFILE_STRING) ; CS%lee_wave_profile = STLAURENT_02
      case (POLZIN_PROFILE_STRING) ; CS%lee_wave_profile = POLZIN_09
      case default
        call MOM_error(FATAL, "tidal_mixing_init: Unrecognized setting "// &
            "#define LEE_WAVE_PROFILE "//trim(tmpstr)//" found in input file.")
    end select
  endif

  call get_param(param_file, mdl, "INT_TIDE_LOWMODE_DISSIPATION", CS%Lowmode_itidal_dissipation, &
                 "If true, consider mixing due to breaking low modes that "//&
                 "have been remotely generated; as with itidal drag on the "//&
                 "barotropic tide, use an internal tidal dissipation scheme to "//&
                 "drive diapycnal mixing, along the lines of St. Laurent "//&
                 "et al. (2002) and Simmons et al. (2004).", default=.false.)

  if ((CS%Int_tide_dissipation .and. (CS%int_tide_profile == POLZIN_09)) .or. &
      (CS%lee_wave_dissipation .and. (CS%lee_wave_profile == POLZIN_09))) then
    if (CS%use_CVMix_tidal) then
        call MOM_error(FATAL, "tidal_mixing_init: Polzin scheme cannot "// &
            "be used when CVMix tidal mixing scheme is active.")
    endif
    call get_param(param_file, mdl, "NU_POLZIN", CS%Nu_Polzin, &
                 "When the Polzin decay profile is used, this is a "//&
                 "non-dimensional constant in the expression for the "//&
                 "vertical scale of decay for the tidal energy dissipation.", &
                 units="nondim", default=0.0697)
    call get_param(param_file, mdl, "NBOTREF_POLZIN", CS%Nbotref_Polzin, &
                 "When the Polzin decay profile is used, this is the "//&
                 "reference value of the buoyancy frequency at the ocean "//&
                 "bottom in the Polzin formulation for the vertical "//&
                 "scale of decay for the tidal energy dissipation.", &
                 units="s-1", default=9.61e-4, scale=US%T_to_s)
    call get_param(param_file, mdl, "POLZIN_DECAY_SCALE_FACTOR", &
                 CS%Polzin_decay_scale_factor, &
                 "When the Polzin decay profile is used, this is a "//&
                 "scale factor for the vertical scale of decay of the tidal "//&
                 "energy dissipation.", default=1.0, units="nondim")
    call get_param(param_file, mdl, "POLZIN_SCALE_MAX_FACTOR", &
                 CS%Polzin_decay_scale_max_factor, &
                 "When the Polzin decay profile is used, this is a factor "//&
                 "to limit the vertical scale of decay of the tidal "//&
                 "energy dissipation to POLZIN_DECAY_SCALE_MAX_FACTOR "//&
                 "times the depth of the ocean.", units="nondim", default=1.0)
    call get_param(param_file, mdl, "POLZIN_MIN_DECAY_SCALE", CS%Polzin_min_decay_scale, &
                 "When the Polzin decay profile is used, this is the "//&
                 "minimum vertical decay scale for the vertical profile\n"//&
                 "of internal tide dissipation with the Polzin (2009) formulation", &
                 units="m", default=0.0, scale=US%m_to_Z)
  endif

  if (CS%Int_tide_dissipation .or. CS%Lee_wave_dissipation) then
    call get_param(param_file, mdl, "INT_TIDE_DECAY_SCALE", CS%Int_tide_decay_scale, &
                 "The decay scale away from the bottom for tidal TKE with "//&
                 "the new coding when INT_TIDE_DISSIPATION is used.", &
                 !units="m", default=0.0)
                 units="m", default=500.0, scale=US%m_to_Z)  ! TODO: confirm this new default
    call get_param(param_file, mdl, "MU_ITIDES", CS%Mu_itides, &
                 "A dimensionless turbulent mixing efficiency used with "//&
                 "INT_TIDE_DISSIPATION, often 0.2.", units="nondim", default=0.2)
    call get_param(param_file, mdl, "GAMMA_ITIDES", CS%Gamma_itides, &
                 "The fraction of the internal tidal energy that is "//&
                 "dissipated locally with INT_TIDE_DISSIPATION. "//&
                 "THIS NAME COULD BE BETTER.", &
                 units="nondim", default=0.3333)
    call get_param(param_file, mdl, "MIN_ZBOT_ITIDES", CS%min_zbot_itides, &
                 "Turn off internal tidal dissipation when the total "//&
                 "ocean depth is less than this value.", units="m", default=0.0, scale=US%m_to_Z)
  endif

  if ( (CS%Int_tide_dissipation .or. CS%Lee_wave_dissipation) .and. &
        .not. CS%use_CVMix_tidal) then

    call safe_alloc_ptr(CS%Nb,isd,ied,jsd,jed)
    call safe_alloc_ptr(CS%h2,isd,ied,jsd,jed)
    call safe_alloc_ptr(CS%TKE_itidal,isd,ied,jsd,jed)
    call safe_alloc_ptr(CS%mask_itidal,isd,ied,jsd,jed) ; CS%mask_itidal(:,:) = 1.0

    call get_param(param_file, mdl, "KAPPA_ITIDES", CS%kappa_itides, &
                 "A topographic wavenumber used with INT_TIDE_DISSIPATION. "//&
                 "The default is 2pi/10 km, as in St.Laurent et al. 2002.", &
                 units="m-1", default=8.e-4*atan(1.0), scale=US%Z_to_m)

    call get_param(param_file, mdl, "UTIDE", CS%utide, &
                 "The constant tidal amplitude used with INT_TIDE_DISSIPATION.", &
                 units="m s-1", default=0.0, scale=US%m_to_Z*US%T_to_s)
    call safe_alloc_ptr(CS%tideamp,is,ie,js,je) ; CS%tideamp(:,:) = CS%utide

    call get_param(param_file, mdl, "KAPPA_H2_FACTOR", CS%kappa_h2_factor, &
                 "A scaling factor for the roughness amplitude with "//&
                 "INT_TIDE_DISSIPATION.",  units="nondim", default=1.0)
    call get_param(param_file, mdl, "TKE_ITIDE_MAX", CS%TKE_itide_max, &
                 "The maximum internal tide energy source available to mix "//&
                 "above the bottom boundary layer with INT_TIDE_DISSIPATION.", &
                 units="W m-2", default=1.0e3, scale=US%W_m2_to_RZ3_T3)

    call get_param(param_file, mdl, "READ_TIDEAMP", read_tideamp, &
                 "If true, read a file (given by TIDEAMP_FILE) containing "//&
                 "the tidal amplitude with INT_TIDE_DISSIPATION.", default=.false.)
    if (read_tideamp) then
      if (CS%use_CVMix_tidal) then
          call MOM_error(FATAL, "tidal_mixing_init: Tidal amplitude files are "// &
              "not compatible with CVMix tidal mixing. ")
      endif
      call get_param(param_file, mdl, "TIDEAMP_FILE", tideamp_file, &
                 "The path to the file containing the spatially varying "//&
                 "tidal amplitudes with INT_TIDE_DISSIPATION.", default="tideamp.nc")
      filename = trim(CS%inputdir) // trim(tideamp_file)
      call log_param(param_file, mdl, "INPUTDIR/TIDEAMP_FILE", filename)
      call MOM_read_data(filename, 'tideamp', CS%tideamp, G%domain, scale=US%m_to_Z*US%T_to_s)
    endif

    call get_param(param_file, mdl, "H2_FILE", h2_file, &
                 "The path to the file containing the sub-grid-scale "//&
                 "topographic roughness amplitude with INT_TIDE_DISSIPATION.", &
                 fail_if_missing=(.not.CS%use_CVMix_tidal))
    filename = trim(CS%inputdir) // trim(h2_file)
    call log_param(param_file, mdl, "INPUTDIR/H2_FILE", filename)
    call MOM_read_data(filename, 'h2', CS%h2, G%domain, scale=US%m_to_Z**2)

    call get_param(param_file, mdl, "FRACTIONAL_ROUGHNESS_MAX", max_frac_rough, &
                 "The maximum topographic roughness amplitude as a fraction of the mean depth, "//&
                 "or a negative value for no limitations on roughness.", &
                 units="nondim", default=0.1)

    do j=js,je ; do i=is,ie
      if (G%bathyT(i,j) < CS%min_zbot_itides) CS%mask_itidal(i,j) = 0.0
      CS%tideamp(i,j) = CS%tideamp(i,j) * CS%mask_itidal(i,j) * G%mask2dT(i,j)

      ! Restrict rms topo to a fraction (often 10 percent) of the column depth.
      if (CS%answers_2018 .and. (max_frac_rough >= 0.0)) then
        hamp = min(max_frac_rough*G%bathyT(i,j), sqrt(CS%h2(i,j)))
        CS%h2(i,j) = hamp*hamp
      else
        if (max_frac_rough >= 0.0) &
          CS%h2(i,j) = min((max_frac_rough*G%bathyT(i,j))**2, CS%h2(i,j))
      endif

      utide = CS%tideamp(i,j)
      ! Compute the fixed part of internal tidal forcing.
      ! The units here are [R Z3 T-2 ~> J m-2 = kg s-2] here.
      CS%TKE_itidal(i,j) = 0.5 * CS%kappa_h2_factor * GV%Rho0 * &
           CS%kappa_itides * CS%h2(i,j) * utide*utide
    enddo ; enddo

  endif

  if (CS%Lee_wave_dissipation) then

    call get_param(param_file, mdl, "NIKURASHIN_TKE_INPUT_FILE",Niku_TKE_input_file, &
                 "The path to the file containing the TKE input from lee "//&
                 "wave driven mixing. Used with LEE_WAVE_DISSIPATION.", &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, "NIKURASHIN_SCALE",Niku_scale, &
                 "A non-dimensional factor by which to scale the lee-wave "//&
                 "driven TKE input. Used with LEE_WAVE_DISSIPATION.", &
                 units="nondim", default=1.0)

    filename = trim(CS%inputdir) // trim(Niku_TKE_input_file)
    call log_param(param_file, mdl, "INPUTDIR/NIKURASHIN_TKE_INPUT_FILE", &
                   filename)
    call safe_alloc_ptr(CS%TKE_Niku,is,ie,js,je) ; CS%TKE_Niku(:,:) = 0.0
    call MOM_read_data(filename, 'TKE_input', CS%TKE_Niku, G%domain, timelevel=1, &  ! ??? timelevel -aja
                       scale=Niku_scale*US%W_m2_to_RZ3_T3)

    call get_param(param_file, mdl, "GAMMA_NIKURASHIN",CS%Gamma_lee, &
                 "The fraction of the lee wave energy that is dissipated "//&
                 "locally with LEE_WAVE_DISSIPATION.", units="nondim", &
                 default=0.3333)
    call get_param(param_file, mdl, "DECAY_SCALE_FACTOR_LEE",CS%Decay_scale_factor_lee, &
                 "Scaling for the vertical decay scaleof the local "//&
                 "dissipation of lee waves dissipation.", units="nondim", &
                 default=1.0)
  else
    CS%Decay_scale_factor_lee = -9.e99 ! This should never be used if CS%Lee_wave_dissipation = False
  endif

  ! Configure CVMix
  if (CS%use_CVMix_tidal) then

    ! Read in CVMix params
    !call openParameterBlock(param_file,'CVMix_TIDAL')
    call get_param(param_file, mdl, "TIDAL_MAX_COEF", CS%tidal_max_coef, &
                   "largest acceptable value for tidal diffusivity", &
                   units="m^2/s", default=50e-4) ! the default is 50e-4 in CVMix, 100e-4 in POP.
    call get_param(param_file, mdl, "TIDAL_DISS_LIM_TC", CS%tidal_diss_lim_tc, &
                   "Min allowable depth for dissipation for tidal-energy-constituent data. "//&
                   "No dissipation contribution is applied above TIDAL_DISS_LIM_TC.", &
                   units="m", default=0.0, scale=US%m_to_Z)
    call get_param(param_file, mdl, "TIDAL_ENERGY_FILE",tidal_energy_file, &
                 "The path to the file containing tidal energy "//&
                 "dissipation. Used with CVMix tidal mixing schemes.", &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, 'MIN_THICKNESS', CS%min_thickness, default=0.001, &
                   do_not_log=.True.)
    call get_param(param_file, mdl, "PRANDTL_TIDAL", prandtl_tidal, &
                   "Prandtl number used by CVMix tidal mixing schemes "//&
                   "to convert vertical diffusivities into viscosities.", &
                    units="nondim", default=1.0, &
                   do_not_log=.true.)
    call CVMix_put(CS%CVMix_glb_params,'Prandtl',prandtl_tidal)

    tidal_energy_file = trim(CS%inputdir) // trim(tidal_energy_file)
    call get_param(param_file, mdl, "TIDAL_ENERGY_TYPE",tidal_energy_type, &
                 "The type of input tidal energy flux dataset. Valid values are"//&
                   "\t Jayne\n"//&
                   "\t ER03 \n",&
                 fail_if_missing=.true.)
    ! Check whether tidal energy input format and CVMix tidal mixing scheme are consistent
    if ( .not. ( &
          (uppercase(tidal_energy_type(1:4)).eq.'JAYN' .and. CS%CVMix_tidal_scheme.eq.SIMMONS).or. &
          (uppercase(tidal_energy_type(1:4)).eq.'ER03' .and. CS%CVMix_tidal_scheme.eq.SCHMITTNER) ) )then
        call MOM_error(FATAL, "tidal_mixing_init: Tidal energy file type ("//&
                      trim(tidal_energy_type)//") is incompatible with CVMix tidal "//&
                      " mixing scheme: "//trim(CVMix_tidal_scheme_str) )
    endif
    CVMix_tidal_scheme_str = lowercase(CVMix_tidal_scheme_str)

    ! Set up CVMix
    call CVMix_init_tidal(CVmix_tidal_params_user = CS%CVMix_tidal_params,    &
                          mix_scheme              = CVMix_tidal_scheme_str,   &
                          efficiency              = CS%Mu_itides,             &
                          vertical_decay_scale    = CS%int_tide_decay_scale*US%Z_to_m,  &
                          max_coefficient         = CS%tidal_max_coef,        &
                          local_mixing_frac       = CS%Gamma_itides,          &
                          depth_cutoff            = CS%min_zbot_itides*US%Z_to_m)

    call read_tidal_energy(G, US, tidal_energy_type, tidal_energy_file, CS)

    !call closeParameterBlock(param_file)

  endif ! CVMix on

  ! Register Diagnostics fields

  if (CS%Int_tide_dissipation .or. CS%Lee_wave_dissipation .or. &
      CS%Lowmode_itidal_dissipation) then

    CS%id_Kd_itidal = register_diag_field('ocean_model','Kd_itides',diag%axesTi,Time, &
         'Internal Tide Driven Diffusivity', 'm2 s-1', conversion=US%Z2_T_to_m2_s)

    if (CS%use_CVMix_tidal) then
      CS%id_N2_int = register_diag_field('ocean_model','N2_int',diag%axesTi,Time, &
          'Bouyancy frequency squared, at interfaces', 's-2', conversion=US%s_to_T**2)
      !> TODO: add units
      CS%id_Simmons_coeff = register_diag_field('ocean_model','Simmons_coeff',diag%axesT1,Time, &
           'time-invariant portion of the tidal mixing coefficient using the Simmons', '')
      CS%id_Schmittner_coeff = register_diag_field('ocean_model','Schmittner_coeff',diag%axesTL,Time, &
           'time-invariant portion of the tidal mixing coefficient using the Schmittner', '')
      CS%id_tidal_qe_md = register_diag_field('ocean_model','tidal_qe_md',diag%axesTL,Time, &
           'input tidal energy dissipated locally interpolated to model vertical coordinates', '')
      CS%id_vert_dep = register_diag_field('ocean_model','vert_dep',diag%axesTi,Time, &
           'vertical deposition function needed for Simmons et al tidal  mixing', '')

    else
      CS%id_TKE_itidal = register_diag_field('ocean_model','TKE_itidal',diag%axesT1,Time, &
          'Internal Tide Driven Turbulent Kinetic Energy', &
          'W m-2', conversion=US%RZ3_T3_to_W_m2)
      CS%id_Nb = register_diag_field('ocean_model','Nb',diag%axesT1,Time, &
           'Bottom Buoyancy Frequency', 's-1', conversion=US%s_to_T)

      CS%id_Kd_lowmode = register_diag_field('ocean_model','Kd_lowmode',diag%axesTi,Time, &
           'Internal Tide Driven Diffusivity (from propagating low modes)', &
           'm2 s-1', conversion=US%Z2_T_to_m2_s)

      CS%id_Fl_itidal = register_diag_field('ocean_model','Fl_itides',diag%axesTi,Time, &
          'Vertical flux of tidal turbulent dissipation', &
          'm3 s-3', conversion=(US%Z_to_m**3*US%s_to_T**3))

      CS%id_Fl_lowmode = register_diag_field('ocean_model','Fl_lowmode',diag%axesTi,Time, &
           'Vertical flux of tidal turbulent dissipation (from propagating low modes)', &
           'm3 s-3', conversion=(US%Z_to_m**3*US%s_to_T**3))

      CS%id_Polzin_decay_scale = register_diag_field('ocean_model','Polzin_decay_scale',diag%axesT1,Time, &
           'Vertical decay scale for the tidal turbulent dissipation with Polzin scheme', &
           'm', conversion=US%Z_to_m)

      CS%id_Polzin_decay_scale_scaled = register_diag_field('ocean_model', &
           'Polzin_decay_scale_scaled', diag%axesT1, Time, &
           'Vertical decay scale for the tidal turbulent dissipation with Polzin scheme, '// &
           'scaled by N2_bot/N2_meanz', 'm', conversion=US%Z_to_m)

      CS%id_N2_bot = register_diag_field('ocean_model','N2_b',diag%axesT1,Time, &
           'Bottom Buoyancy frequency squared', 's-2', conversion=US%s_to_T**2)

      CS%id_N2_meanz = register_diag_field('ocean_model','N2_meanz', diag%axesT1, Time, &
           'Buoyancy frequency squared averaged over the water column', 's-2', conversion=US%s_to_T**2)

      CS%id_Kd_Itidal_Work = register_diag_field('ocean_model','Kd_Itidal_Work',diag%axesTL,Time, &
           'Work done by Internal Tide Diapycnal Mixing', &
           'W m-2', conversion=US%RZ3_T3_to_W_m2)

      CS%id_Kd_Niku_Work = register_diag_field('ocean_model','Kd_Nikurashin_Work',diag%axesTL,Time, &
           'Work done by Nikurashin Lee Wave Drag Scheme', &
           'W m-2', conversion=US%RZ3_T3_to_W_m2)

      CS%id_Kd_Lowmode_Work = register_diag_field('ocean_model','Kd_Lowmode_Work',diag%axesTL,Time, &
           'Work done by Internal Tide Diapycnal Mixing (low modes)', &
           'W m-2', conversion=US%RZ3_T3_to_W_m2)

      if (CS%Lee_wave_dissipation) then
        CS%id_TKE_leewave = register_diag_field('ocean_model','TKE_leewave',diag%axesT1,Time, &
            'Lee wave Driven Turbulent Kinetic Energy', &
            'W m-2', conversion=US%RZ3_T3_to_W_m2)
        CS%id_Kd_Niku = register_diag_field('ocean_model','Kd_Nikurashin',diag%axesTi,Time, &
            'Lee Wave Driven Diffusivity', 'm2 s-1', conversion=US%Z2_T_to_m2_s)
      endif
    endif ! S%use_CVMix_tidal
  endif

end function tidal_mixing_init


!> Depending on whether or not CVMix is active, calls the associated subroutine to compute internal
!! tidal dissipation and to add the effect of internal-tide-driven mixing to the layer or interface
!! diffusivities.
subroutine calculate_tidal_mixing(h, N2_bot, j, TKE_to_Kd, max_TKE, G, GV, US, CS, &
                                  N2_lay, N2_int, Kd_lay, Kd_int, Kd_max, Kv)
  type(ocean_grid_type),            intent(in)    :: G      !< The ocean's grid structure
  type(verticalGrid_type),          intent(in)    :: GV     !< The ocean's vertical grid structure
  type(unit_scale_type),            intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G)),         intent(in)    :: N2_bot !< The near-bottom squared buoyancy
                                                            !! frequency [T-2 ~> s-2].
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: N2_lay !< The squared buoyancy frequency of the
                                                            !! layers [T-2 ~> s-2].
  real, dimension(SZI_(G),SZK_(GV)+1), intent(in) :: N2_int !< The squared buoyancy frequency at the
                                                            !! interfaces [T-2 ~> s-2].
  integer,                          intent(in)    :: j      !< The j-index to work on
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: TKE_to_Kd !< The conversion rate between the TKE
                                                            !! dissipated within a layer and the
                                                            !! diapycnal diffusivity within that layer,
                                                            !! usually (~Rho_0 / (G_Earth * dRho_lay))
                                                            !! [Z2 T-1 / Z3 T-3 = T2 Z-1 ~> s2 m-1]
  real, dimension(SZI_(G),SZK_(GV)), intent(in)   :: max_TKE !< The energy required to for a layer to entrain
                                                            !! to its maximum realizable thickness [Z3 T-3 ~> m3 s-3]
  type(tidal_mixing_cs),            pointer       :: CS     !< The control structure for this module
  real, dimension(SZI_(G),SZK_(GV)), &
                          optional, intent(inout) :: Kd_lay !< The diapycnal diffusivity in layers [Z2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1), &
                          optional, intent(inout) :: Kd_int !< The diapycnal diffusivity at interfaces,
                                                            !! [Z2 T-1 ~> m2 s-1].
  real,                             intent(in)    :: Kd_max !< The maximum increment for diapycnal
                                                            !! diffusivity due to TKE-based processes,
                                                            !! [Z2 T-1 ~> m2 s-1].
                                                            !! Set this to a negative value to have no limit.
  real, dimension(:,:,:),           pointer       :: Kv     !< The "slow" vertical viscosity at each interface
                                                            !! (not layer!) [Z2 T-1 ~> m2 s-1].

  if (CS%Int_tide_dissipation .or. CS%Lee_wave_dissipation .or. CS%Lowmode_itidal_dissipation) then
    if (CS%use_CVMix_tidal) then
      call calculate_CVMix_tidal(h, j, G, GV, US, CS, N2_int, Kd_lay, Kd_int, Kv)
    else
      call add_int_tide_diffusivity(h, N2_bot, j, TKE_to_Kd, max_TKE, &
                                    G, GV, US, CS, N2_lay, Kd_lay, Kd_int, Kd_max)
    endif
  endif
end subroutine calculate_tidal_mixing


!> Calls the CVMix routines to compute tidal dissipation and to add the effect of internal-tide-driven
!! mixing to the interface diffusivities.
subroutine calculate_CVMix_tidal(h, j, G, GV, US, CS, N2_int, Kd_lay, Kd_int, Kv)
  integer,                 intent(in)    :: j     !< The j-index to work on
  type(ocean_grid_type),   intent(in)    :: G     !< Grid structure.
  type(verticalGrid_type), intent(in)    :: GV    !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US    !< A dimensional unit scaling type
  type(tidal_mixing_cs),   pointer       :: CS    !< This module's control structure.
  real, dimension(SZI_(G),SZK_(GV)+1), intent(in) :: N2_int !< The squared buoyancy
                                                  !! frequency at the interfaces [T-2 ~> s-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)), &
                 optional, intent(inout) :: Kd_lay!< The diapycnal diffusivity in the layers [Z2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1), &
                 optional, intent(inout) :: Kd_int!< The diapycnal diffusivity at interfaces [Z2 T-1 ~> m2 s-1].
  real, dimension(:,:,:),  pointer       :: Kv    !< The "slow" vertical viscosity at each interface
                                                  !! (not layer!) [Z2 T-1 ~> m2 s-1].
  ! Local variables
  real, dimension(SZK_(GV)+1) :: Kd_tidal    ! tidal diffusivity [m2 s-1]
  real, dimension(SZK_(GV)+1) :: Kv_tidal    ! tidal viscosity [m2 s-1]
  real, dimension(SZK_(GV)+1) :: vert_dep    ! vertical deposition
  real, dimension(SZK_(GV)+1) :: iFaceHeight ! Height of interfaces [m]
  real, dimension(SZK_(GV)+1) :: SchmittnerSocn
  real, dimension(SZK_(GV))   :: cellHeight  ! Height of cell centers [m]
  real, dimension(SZK_(GV))   :: tidal_qe_md ! Tidal dissipation energy interpolated from 3d input
                                             ! to model coordinates
  real, dimension(SZK_(GV)+1) :: N2_int_i    ! De-scaled interface buoyancy frequency [s-2]
  real, dimension(SZK_(GV))   :: Schmittner_coeff
  real, dimension(SZK_(GV))   :: h_m         ! Cell thickness [m]
  real, allocatable, dimension(:,:) :: exp_hab_zetar

  integer :: i, k, is, ie
  real :: dh, hcorr, Simmons_coeff
  real, parameter :: rho_fw = 1000.0 ! fresh water density [kg/m^3]
                                     ! TODO: when coupled, get this from CESM (SHR_CONST_RHOFW)
  type(tidal_mixing_diags), pointer :: dd => NULL()

  is  = G%isc ; ie  = G%iec
  dd => CS%dd

  select case (CS%CVMix_tidal_scheme)
  case (SIMMONS)
    do i=is,ie

      if (G%mask2dT(i,j)<1) cycle

      iFaceHeight = 0.0 ! BBL is all relative to the surface
      hcorr = 0.0
      do k=1,GV%ke
        ! cell center and cell bottom in meters (negative values in the ocean)
        dh = h(i,j,k) * GV%H_to_m ! Nominal thickness to use for increment, rescaled to m for use by CVMix.
        dh = dh + hcorr ! Take away the accumulated error (could temporarily make dh<0)
        hcorr = min( dh - CS%min_thickness, 0. ) ! If inflating then hcorr<0
        dh = max( dh, CS%min_thickness ) ! Limit increment dh>=min_thickness
        cellHeight(k)    = iFaceHeight(k) - 0.5 * dh
        iFaceHeight(k+1) = iFaceHeight(k) - dh
      enddo

      call CVMix_compute_Simmons_invariant( nlev                    = GV%ke,               &
                                            energy_flux             = CS%tidal_qe_2d(i,j), &
                                            rho                     = rho_fw,              &
                                            SimmonsCoeff            = Simmons_coeff,       &
                                            VertDep                 = vert_dep,            &
                                            zw                      = iFaceHeight,         &
                                            zt                      = cellHeight,          &
                                            CVMix_tidal_params_user = CS%CVMix_tidal_params)

      ! Since we pass tidal_qe_2d=(CS%Gamma_itides)*tidal_energy_flux_2d, and not tidal_energy_flux_2d in
      ! above subroutine call, we divide Simmons_coeff by CS%Gamma_itides as a corrective step:
      ! TODO: (CS%Gamma_itides)*tidal_energy_flux_2d is unnecessary, directly use tidal_energy_flux_2d
      Simmons_coeff = Simmons_coeff / CS%Gamma_itides


      ! XXX: Temporary de-scaling of N2_int(i,:) into a temporary variable
      do K=1,GV%ke+1
        N2_int_i(K) = US%s_to_T**2 * N2_int(i,K)
      enddo

      call CVMix_coeffs_tidal( Mdiff_out               = Kv_tidal,             &
                               Tdiff_out               = Kd_tidal,             &
                               Nsqr                    = N2_int_i,             &
                               OceanDepth              = -iFaceHeight(GV%ke+1),&
                               SimmonsCoeff            = Simmons_coeff,        &
                               vert_dep                = vert_dep,             &
                               nlev                    = GV%ke,                &
                               max_nlev                = GV%ke,                &
                               CVMix_params            = CS%CVMix_glb_params,  &
                               CVMix_tidal_params_user = CS%CVMix_tidal_params)

      ! Update diffusivity
      if (present(Kd_lay)) then
        do k=1,GV%ke
          Kd_lay(i,k) = Kd_lay(i,k) + 0.5 * US%m2_s_to_Z2_T * (Kd_tidal(k) + Kd_tidal(k+1))
        enddo
      endif
      if (present(Kd_int)) then
        do K=1,GV%ke+1
          Kd_int(i,K) = Kd_int(i,K) +  (US%m2_s_to_Z2_T * Kd_tidal(K))
        enddo
      endif
      ! Update viscosity with the proper unit conversion.
      if (associated(Kv)) then
        do K=1,GV%ke+1
          Kv(i,j,K) = Kv(i,j,K) + US%m2_s_to_Z2_T * Kv_tidal(K)  ! Rescale from m2 s-1 to Z2 T-1.
        enddo
      endif

      ! diagnostics
      if (associated(dd%Kd_itidal)) then
        dd%Kd_itidal(i,j,:) = US%m2_s_to_Z2_T*Kd_tidal(:)
      endif
      if (associated(dd%N2_int)) then
        dd%N2_int(i,j,:) = N2_int(i,:)
      endif
      if (associated(dd%Simmons_coeff_2d)) then
        dd%Simmons_coeff_2d(i,j) = Simmons_coeff
      endif
      if (associated(dd%vert_dep_3d)) then
        dd%vert_dep_3d(i,j,:) = vert_dep(:)
      endif

    enddo ! i=is,ie

  case (SCHMITTNER)

    ! TODO: correct exp_hab_zetar shapes in CVMix_compute_Schmittner_invariant
    ! and CVMix_compute_SchmittnerCoeff low subroutines

    allocate(exp_hab_zetar(GV%ke+1,GV%ke+1))

    do i=is,ie

      if (G%mask2dT(i,j)<1) cycle

      iFaceHeight = 0.0 ! BBL is all relative to the surface
      hcorr = 0.0
      do k=1,GV%ke
        h_m(k) = h(i,j,k)*GV%H_to_m  ! Rescale thicknesses to m for use by CVmix.
        ! cell center and cell bottom in meters (negative values in the ocean)
        dh = h_m(k) + hcorr ! Nominal thickness less the accumulated error (could temporarily make dh<0)
        hcorr = min( dh - CS%min_thickness, 0. ) ! If inflating then hcorr<0
        dh = max( dh, CS%min_thickness ) ! Limit increment dh>=min_thickness
        cellHeight(k)    = iFaceHeight(k) - 0.5 * dh
        iFaceHeight(k+1) = iFaceHeight(k) - dh
      enddo

      SchmittnerSocn = 0.0 ! TODO: compute this

      ! form the time-invariant part of Schmittner coefficient term
      call CVMix_compute_Schmittner_invariant(nlev                    = GV%ke,          &
                                              VertDep                 = vert_dep,       &
                                              efficiency              = CS%Mu_itides,   &
                                              rho                     = rho_fw,         &
                                              exp_hab_zetar           = exp_hab_zetar,  &
                                              zw                      = iFaceHeight,    &
                                              CVmix_tidal_params_user = CS%CVMix_tidal_params)
                  !TODO: in above call, there is no need to pass efficiency, since it gets
                  ! passed via CVMix_init_tidal and stored in CVMix_tidal_params. Change
                  ! CVMix API to prevent this redundancy.

      ! remap from input z coordinate to model coordinate:
      tidal_qe_md = 0.0
      call remapping_core_h(CS%remap_cs, size(CS%h_src), CS%h_src, CS%tidal_qe_3d_in(i,j,:), &
                            GV%ke, h_m, tidal_qe_md)

      ! form the Schmittner coefficient that is based on 3D q*E, which is formed from
      ! summing q_i*TidalConstituent_i over the number of constituents.
      call CVMix_compute_SchmittnerCoeff( nlev                    = GV%ke,              &
                                          energy_flux             = tidal_qe_md(:),     &
                                          SchmittnerCoeff         = Schmittner_coeff,   &
                                          exp_hab_zetar           = exp_hab_zetar,      &
                                          CVmix_tidal_params_user = CS%CVMix_tidal_params)

      ! XXX: Temporary de-scaling of N2_int(i,:) into a temporary variable
      do k=1,GV%ke+1
        N2_int_i(k) = US%s_to_T**2 * N2_int(i,k)
      enddo

      call CVMix_coeffs_tidal_schmittner( Mdiff_out               = Kv_tidal,             &
                                          Tdiff_out               = Kd_tidal,             &
                                          Nsqr                    = N2_int_i,             &
                                          OceanDepth              = -iFaceHeight(GV%ke+1), &
                                          nlev                    = GV%ke,                &
                                          max_nlev                = GV%ke,                &
                                          SchmittnerCoeff         = Schmittner_coeff,     &
                                          SchmittnerSouthernOcean = SchmittnerSocn,       &
                                          CVmix_params            = CS%CVMix_glb_params,  &
                                          CVmix_tidal_params_user = CS%CVMix_tidal_params)

      ! Update diffusivity
      if (present(Kd_lay)) then
        do k=1,GV%ke
          Kd_lay(i,k) = Kd_lay(i,k) + 0.5 * US%m2_s_to_Z2_T * (Kd_tidal(k) + Kd_tidal(k+1))
        enddo
      endif
      if (present(Kd_int)) then
        do K=1,GV%ke+1
          Kd_int(i,K) = Kd_int(i,K) +  (US%m2_s_to_Z2_T * Kd_tidal(K))
        enddo
      endif

      ! Update viscosity
      if (associated(Kv)) then
        do K=1,GV%ke+1
          Kv(i,j,K) = Kv(i,j,K) + US%m2_s_to_Z2_T * Kv_tidal(K)   ! Rescale from m2 s-1 to Z2 T-1.
        enddo
      endif

      ! diagnostics
      if (associated(dd%Kd_itidal)) then
        dd%Kd_itidal(i,j,:) = US%m2_s_to_Z2_T*Kd_tidal(:)
      endif
      if (associated(dd%N2_int)) then
        dd%N2_int(i,j,:) = N2_int(i,:)
      endif
      if (associated(dd%Schmittner_coeff_3d)) then
        dd%Schmittner_coeff_3d(i,j,:) = Schmittner_coeff(:)
      endif
      if (associated(dd%tidal_qe_md)) then
        dd%tidal_qe_md(i,j,:) = tidal_qe_md(:)
      endif
      if (associated(dd%vert_dep_3d)) then
        dd%vert_dep_3d(i,j,:) = vert_dep(:)
      endif
    enddo ! i=is,ie

    deallocate(exp_hab_zetar)

  case default
    call MOM_error(FATAL, "tidal_mixing_init: Unrecognized setting "// &
         "#define CVMIX_TIDAL_SCHEME found in input file.")
  end select

end subroutine calculate_CVMix_tidal


!> This subroutine adds the effect of internal-tide-driven mixing to the layer diffusivities.
!! The mechanisms considered are (1) local dissipation of internal waves generated by the
!! barotropic flow ("itidal"), (2) local dissipation of internal waves generated by the propagating
!! low modes (rays) of the internal tide ("lowmode"), and (3) local dissipation of internal lee waves.
!! Will eventually need to add diffusivity due to other wave-breaking processes (e.g. Bottom friction,
!! Froude-number-depending breaking, PSI, etc.).
subroutine add_int_tide_diffusivity(h, N2_bot, j, TKE_to_Kd, max_TKE, G, GV, US, CS, &
                                    N2_lay, Kd_lay, Kd_int, Kd_max)
  type(ocean_grid_type),             intent(in)    :: G      !< The ocean's grid structure
  type(verticalGrid_type),           intent(in)    :: GV     !< The ocean's vertical grid structure
  type(unit_scale_type),             intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                     intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G)),          intent(in)    :: N2_bot !< The near-bottom squared buoyancy frequency
                                                             !! frequency [T-2 ~> s-2].
  real, dimension(SZI_(G),SZK_(GV)), intent(in)    :: N2_lay !< The squared buoyancy frequency of the
                                                             !! layers [T-2 ~> s-2].
  integer,                           intent(in)    :: j      !< The j-index to work on
  real, dimension(SZI_(G),SZK_(GV)), intent(in)    :: TKE_to_Kd !< The conversion rate between the TKE
                                                             !! dissipated within a layer and the
                                                             !! diapycnal diffusivity within that layer,
                                                             !! usually (~Rho_0 / (G_Earth * dRho_lay))
                                                             !! [Z2 T-1 / Z3 T-3 = T2 Z-1 ~> s2 m-1]
  real, dimension(SZI_(G),SZK_(GV)), intent(in)    :: max_TKE !< The energy required to for a layer to entrain
                                                             !! to its maximum realizable thickness [Z3 T-3 ~> m3 s-3]
  type(tidal_mixing_cs),             pointer       :: CS     !< The control structure for this module
  real, dimension(SZI_(G),SZK_(GV)), &
                           optional, intent(inout) :: Kd_lay !< The diapycnal diffusivity in layers [Z2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZK_(GV)+1), &
                           optional, intent(inout) :: Kd_int !< The diapycnal diffusivity at interfaces
                                                             !! [Z2 T-1 ~> m2 s-1].
  real,                              intent(in)    :: Kd_max !< The maximum increment for diapycnal
                                                             !! diffusivity due to TKE-based processes
                                                             !! [Z2 T-1 ~> m2 s-1].
                                                             !! Set this to a negative value to have no limit.

  ! local

  real, dimension(SZI_(G)) :: &
    htot,             & ! total thickness above or below a layer, or the
                        ! integrated thickness in the BBL [Z ~> m].
    htot_WKB,         & ! WKB scaled distance from top to bottom [Z ~> m].
    TKE_itidal_bot,   & ! internal tide TKE at ocean bottom [Z3 T-3 ~> m3 s-3]
    TKE_Niku_bot,     & ! lee-wave TKE at ocean bottom [Z3 T-3 ~> m3 s-3]
    TKE_lowmode_bot,  & ! internal tide TKE at ocean bottom lost from all remote low modes [Z3 T-3 ~> m3 s-3] (BDM)
    Inv_int,          & ! inverse of TKE decay for int tide over the depth of the ocean [nondim]
    Inv_int_lee,      & ! inverse of TKE decay for lee waves over the depth of the ocean [nondim]
    Inv_int_low,      & ! inverse of TKE decay for low modes over the depth of the ocean [nondim] (BDM)
    z0_Polzin,        & ! TKE decay scale in Polzin formulation [Z ~> m].
    z0_Polzin_scaled, & ! TKE decay scale in Polzin formulation [Z ~> m].
                        ! multiplied by N2_bot/N2_meanz to be coherent with the WKB scaled z
                        ! z*=int(N2/N2_bot) * N2_bot/N2_meanz = int(N2/N2_meanz)
                        ! z0_Polzin_scaled = z0_Polzin * N2_bot/N2_meanz
    N2_meanz,         & ! vertically averaged squared buoyancy frequency [T-2] for WKB scaling
    TKE_itidal_rem,   & ! remaining internal tide TKE (from barotropic source) [Z3 T-3 ~> m3 s-3]
    TKE_Niku_rem,     & ! remaining lee-wave TKE [Z3 T-3 ~> m3 s-3]
    TKE_lowmode_rem,  & ! remaining internal tide TKE (from propagating low mode source) [Z3 T-3 ~> m3 s-3] (BDM)
    TKE_frac_top,     & ! fraction of bottom TKE that should appear at top of a layer [nondim]
    TKE_frac_top_lee, & ! fraction of bottom TKE that should appear at top of a layer [nondim]
    TKE_frac_top_lowmode, &
                        ! fraction of bottom TKE that should appear at top of a layer [nondim] (BDM)
    z_from_bot,       & ! distance from bottom [Z ~> m].
    z_from_bot_WKB      ! WKB scaled distance from bottom [Z ~> m].

  real :: I_rho0        ! Inverse of the Boussinesq reference density, i.e. 1 / RHO0 [R-1 ~> m3 kg-1]
  real :: Kd_add        ! diffusivity to add in a layer [Z2 T-1 ~> m2 s-1].
  real :: TKE_itide_lay ! internal tide TKE imparted to a layer (from barotropic) [Z3 T-3 ~> m3 s-3]
  real :: TKE_Niku_lay  ! lee-wave TKE imparted to a layer [Z3 T-3 ~> m3 s-3]
  real :: TKE_lowmode_lay ! internal tide TKE imparted to a layer (from low mode) [Z3 T-3 ~> m3 s-3] (BDM)
  real :: frac_used     ! fraction of TKE that can be used in a layer [nondim]
  real :: Izeta         ! inverse of TKE decay scale [Z-1 ~> m-1].
  real :: Izeta_lee     ! inverse of TKE decay scale for lee waves [Z-1 ~> m-1].
  real :: z0Ps_num      ! The numerator of the unlimited z0_Polzin_scaled [Z T-3 ~> m s-3].
  real :: z0Ps_denom    ! The denominator of the unlimited z0_Polzin_scaled [T-3 ~> s-3].
  real :: z0_psl        ! temporary variable [Z ~> m].
  real :: TKE_lowmode_tot ! TKE from all low modes [R Z3 T-3 ~> W m-2] (BDM)

  logical :: use_Polzin, use_Simmons
  character(len=160) :: mesg  ! The text of an error message
  integer :: i, k, is, ie, nz
  integer :: a, fr, m
  type(tidal_mixing_diags), pointer :: dd => NULL()

  is = G%isc ; ie = G%iec ; nz = GV%ke
  dd => CS%dd

  if (.not.(CS%Int_tide_dissipation .or. CS%Lee_wave_dissipation)) return

  do i=is,ie ; htot(i) = 0.0 ; Inv_int(i) = 0.0 ; Inv_int_lee(i) = 0.0 ; Inv_int_low(i) = 0.0 ;enddo
  do k=1,nz ; do i=is,ie
    htot(i) = htot(i) + GV%H_to_Z*h(i,j,k)
  enddo ; enddo

  I_Rho0 = 1.0 / (GV%Rho0)

  use_Polzin = ((CS%Int_tide_dissipation .and. (CS%int_tide_profile == POLZIN_09)) .or. &
                (CS%lee_wave_dissipation .and. (CS%lee_wave_profile == POLZIN_09)) .or. &
                (CS%Lowmode_itidal_dissipation .and. (CS%int_tide_profile == POLZIN_09)))
  use_Simmons = ((CS%Int_tide_dissipation .and. (CS%int_tide_profile == STLAURENT_02)) .or. &
                 (CS%lee_wave_dissipation .and. (CS%lee_wave_profile == STLAURENT_02)) .or. &
                 (CS%Lowmode_itidal_dissipation .and. (CS%int_tide_profile == STLAURENT_02)))

  ! Calculate parameters for vertical structure of dissipation
  ! Simmons:
  if ( use_Simmons ) then
    Izeta = 1.0 / max(CS%Int_tide_decay_scale, GV%H_subroundoff*GV%H_to_Z)
    Izeta_lee = 1.0 / max(CS%Int_tide_decay_scale*CS%Decay_scale_factor_lee, &
                          GV%H_subroundoff*GV%H_to_Z)
    do i=is,ie
      CS%Nb(i,j) = sqrt(N2_bot(i))
      if (associated(dd%N2_bot)) dd%N2_bot(i,j) = N2_bot(i)
      if ( CS%Int_tide_dissipation ) then
        if (Izeta*htot(i) > 1.0e-14) then ! L'Hospital's version of Adcroft's reciprocal rule.
          Inv_int(i) = 1.0 / (1.0 - exp(-Izeta*htot(i)))
        endif
      endif
      if ( CS%Lee_wave_dissipation ) then
        if (Izeta_lee*htot(i) > 1.0e-14) then  ! L'Hospital's version of Adcroft's reciprocal rule.
          Inv_int_lee(i) = 1.0 / (1.0 - exp(-Izeta_lee*htot(i)))
        endif
      endif
      if ( CS%Lowmode_itidal_dissipation) then
        if (Izeta*htot(i) > 1.0e-14) then ! L'Hospital's version of Adcroft's reciprocal rule.
          Inv_int_low(i) = 1.0 / (1.0 - exp(-Izeta*htot(i)))
        endif
      endif
      z_from_bot(i) = GV%H_to_Z*h(i,j,nz)
    enddo
  endif ! Simmons

  ! Polzin:
  if ( use_Polzin ) then
    ! WKB scaling of the vertical coordinate
    do i=is,ie ; N2_meanz(i) = 0.0 ; enddo
    do k=1,nz ; do i=is,ie
      N2_meanz(i) = N2_meanz(i) + N2_lay(i,k) * GV%H_to_Z * h(i,j,k)
    enddo ; enddo
    do i=is,ie
      N2_meanz(i) = N2_meanz(i) / (htot(i) + GV%H_subroundoff*GV%H_to_Z)
      if (associated(dd%N2_meanz))  dd%N2_meanz(i,j) = N2_meanz(i)
    enddo

    ! WKB scaled z*(z=H) z* at the surface using the modified Polzin WKB scaling
    do i=is,ie ; htot_WKB(i) = htot(i) ; enddo
!    do i=is,ie ; htot_WKB(i) = 0.0 ; enddo
!    do k=1,nz ; do i=is,ie
!      htot_WKB(i) = htot_WKB(i) + GV%H_to_Z*h(i,j,k) * N2_lay(i,k) / N2_meanz(i)
!    enddo ; enddo
    ! htot_WKB(i) = htot(i) ! Nearly equivalent and simpler

    do i=is,ie
      CS%Nb(i,j) = sqrt(N2_bot(i))
      if (CS%answers_2018) then
        if ((CS%tideamp(i,j) > 0.0) .and. &
            (CS%kappa_itides**2 * CS%h2(i,j) * CS%Nb(i,j)**3 > 1.0e-14*US%T_to_s**3) ) then
          z0_polzin(i) = CS%Polzin_decay_scale_factor * CS%Nu_Polzin * &
                         CS%Nbotref_Polzin**2 * CS%tideamp(i,j) / &
                       ( CS%kappa_itides**2 * CS%h2(i,j) * CS%Nb(i,j)**3 )
          if (z0_polzin(i) < CS%Polzin_min_decay_scale) &
            z0_polzin(i) = CS%Polzin_min_decay_scale
          if (N2_meanz(i) > 1.0e-14*US%T_to_s**2  ) then
            z0_polzin_scaled(i) = z0_polzin(i)*CS%Nb(i,j)**2 / N2_meanz(i)
          else
            z0_polzin_scaled(i) = CS%Polzin_decay_scale_max_factor * htot(i)
          endif
          if (z0_polzin_scaled(i) > (CS%Polzin_decay_scale_max_factor * htot(i)) ) &
            z0_polzin_scaled(i) = CS%Polzin_decay_scale_max_factor * htot(i)
        else
          z0_polzin(i) = CS%Polzin_decay_scale_max_factor * htot(i)
          z0_polzin_scaled(i) = CS%Polzin_decay_scale_max_factor * htot(i)
        endif
      else
        z0Ps_num = (CS%Polzin_decay_scale_factor * CS%Nu_Polzin * CS%Nbotref_Polzin**2) * CS%tideamp(i,j)
        z0Ps_denom = ( CS%kappa_itides**2 * CS%h2(i,j) * CS%Nb(i,j) * N2_meanz(i) )
        if ((CS%tideamp(i,j) > 0.0) .and. &
            (z0Ps_num < z0Ps_denom * CS%Polzin_decay_scale_max_factor * htot(i))) then
          z0_polzin_scaled(i) = z0Ps_num / z0Ps_denom

          if (abs(N2_meanz(i) * z0_polzin_scaled(i)) < &
              CS%Nb(i,j)**2 * (CS%Polzin_decay_scale_max_factor * htot(i))) then
            z0_polzin(i) = z0_polzin_scaled(i) * (N2_meanz(i) / CS%Nb(i,j)**2)
          else
            z0_polzin(i) = CS%Polzin_decay_scale_max_factor * htot(i)
          endif
        else
          z0_polzin(i) = CS%Polzin_decay_scale_max_factor * htot(i)
          z0_polzin_scaled(i) = CS%Polzin_decay_scale_max_factor * htot(i)
        endif
      endif

      if (associated(dd%Polzin_decay_scale)) &
        dd%Polzin_decay_scale(i,j) = z0_polzin(i)
      if (associated(dd%Polzin_decay_scale_scaled)) &
        dd%Polzin_decay_scale_scaled(i,j) = z0_polzin_scaled(i)
      if (associated(dd%N2_bot)) dd%N2_bot(i,j) = CS%Nb(i,j)*CS%Nb(i,j)

      if (CS%answers_2018) then
        ! These expressions use dimensional constants to avoid NaN values.
        if ( CS%Int_tide_dissipation .and. (CS%int_tide_profile == POLZIN_09) ) then
          if (htot_WKB(i) > 1.0e-14*US%m_to_Z) &
            Inv_int(i) = ( z0_polzin_scaled(i) / htot_WKB(i) ) + 1.0
        endif
        if ( CS%lee_wave_dissipation .and. (CS%lee_wave_profile == POLZIN_09) ) then
          if (htot_WKB(i) > 1.0e-14*US%m_to_Z) &
            Inv_int_lee(i) = ( z0_polzin_scaled(i)*CS%Decay_scale_factor_lee / htot_WKB(i) ) + 1.0
        endif
        if ( CS%Lowmode_itidal_dissipation .and. (CS%int_tide_profile == POLZIN_09) ) then
          if (htot_WKB(i) > 1.0e-14*US%m_to_Z) &
            Inv_int_low(i) = ( z0_polzin_scaled(i) / htot_WKB(i) ) + 1.0
        endif
      else
        ! These expressions give values of Inv_int < 10^14 using a variant of Adcroft's reciprocal rule.
        Inv_int(i) = 0.0 ; Inv_int_lee(i) = 0.0 ; Inv_int_low(i) = 0.0
        if ( CS%Int_tide_dissipation .and. (CS%int_tide_profile == POLZIN_09) ) then
          if (z0_polzin_scaled(i) < 1.0e14 * htot_WKB(i)) &
            Inv_int(i) = ( z0_polzin_scaled(i) / htot_WKB(i) ) + 1.0
        endif
        if ( CS%lee_wave_dissipation .and. (CS%lee_wave_profile == POLZIN_09) ) then
          if (z0_polzin_scaled(i) < 1.0e14 * htot_WKB(i)) &
            Inv_int_lee(i) = ( z0_polzin_scaled(i)*CS%Decay_scale_factor_lee / htot_WKB(i) ) + 1.0
        endif
        if ( CS%Lowmode_itidal_dissipation .and. (CS%int_tide_profile == POLZIN_09) ) then
          if (z0_polzin_scaled(i) < 1.0e14 * htot_WKB(i)) &
            Inv_int_low(i) = ( z0_polzin_scaled(i) / htot_WKB(i) ) + 1.0
        endif
      endif

      z_from_bot(i) = GV%H_to_Z*h(i,j,nz)
      ! Use the new formulation for WKB scaling.  N2 is referenced to its vertical mean.
      if (CS%answers_2018) then
        if (N2_meanz(i) > 1.0e-14*US%T_to_s**2 ) then
          z_from_bot_WKB(i) = GV%H_to_Z*h(i,j,nz) * N2_lay(i,nz) / N2_meanz(i)
        else ; z_from_bot_WKB(i) = 0 ; endif
      else
        if (GV%H_to_Z*h(i,j,nz) * N2_lay(i,nz) < N2_meanz(i) * (1.0e14 * htot_WKB(i))) then
          z_from_bot_WKB(i) = GV%H_to_Z*h(i,j,nz) * N2_lay(i,nz) / N2_meanz(i)
        else ; z_from_bot_WKB(i) = 0 ; endif
      endif
    enddo
  endif  ! Polzin

  ! Calculate/get dissipation values at bottom
  ! Both Polzin and Simmons:
  do i=is,ie
    ! Dissipation of locally trapped internal tide (non-propagating high modes)
    TKE_itidal_bot(i) = min(CS%TKE_itidal(i,j)*CS%Nb(i,j), CS%TKE_itide_max)
    if (associated(dd%TKE_itidal_used)) &
      dd%TKE_itidal_used(i,j) = TKE_itidal_bot(i)
    TKE_itidal_bot(i) = (I_rho0 * CS%Mu_itides * CS%Gamma_itides) * TKE_itidal_bot(i)
    ! Dissipation of locally trapped lee waves
    TKE_Niku_bot(i) = 0.0
    if (CS%Lee_wave_dissipation) then
      TKE_Niku_bot(i) = (I_rho0 * CS%Mu_itides * CS%Gamma_lee) * CS%TKE_Niku(i,j)
    endif
    ! Dissipation of propagating internal tide (baroclinic low modes; rays) (BDM)
    TKE_lowmode_tot    = 0.0
    TKE_lowmode_bot(i) = 0.0
    if (CS%Lowmode_itidal_dissipation) then
      ! get loss rate due to wave drag on low modes (already multiplied by q)

      ! TODO: uncomment the following call and fix it
      !call get_lowmode_loss(i,j,G,CS%int_tide_CSp,"WaveDrag",TKE_lowmode_tot)
      write (mesg,*) "========", __FILE__, __LINE__
      call MOM_error(FATAL,trim(mesg)//": this block not supported yet. (aa)")

      TKE_lowmode_bot(i) = CS%Mu_itides * I_rho0 * TKE_lowmode_tot
    endif
    ! Vertical energy flux at bottom
    TKE_itidal_rem(i)  = Inv_int(i)     * TKE_itidal_bot(i)
    TKE_Niku_rem(i)    = Inv_int_lee(i) * TKE_Niku_bot(i)
    TKE_lowmode_rem(i) = Inv_int_low(i) * TKE_lowmode_bot(i)

    if (associated(dd%Fl_itidal)) dd%Fl_itidal(i,j,nz) = TKE_itidal_rem(i) !why is this here? BDM
  enddo

  ! Estimate the work that would be done by mixing in each layer.
  ! Simmons:
  if ( use_Simmons ) then
    do k=nz-1,2,-1 ; do i=is,ie
      if (max_TKE(i,k) <= 0.0) cycle
      z_from_bot(i) = z_from_bot(i) + GV%H_to_Z*h(i,j,k)

      ! Fraction of bottom flux predicted to reach top of this layer
      TKE_frac_top(i)         = Inv_int(i)     * exp(-Izeta * z_from_bot(i))
      TKE_frac_top_lee(i)     = Inv_int_lee(i) * exp(-Izeta_lee * z_from_bot(i))
      TKE_frac_top_lowmode(i) = Inv_int_low(i) * exp(-Izeta * z_from_bot(i))

      ! Actual influx at bottom of layer minus predicted outflux at top of layer to give
      ! predicted power expended
      TKE_itide_lay   = TKE_itidal_rem(i)  - TKE_itidal_bot(i) * TKE_frac_top(i)
      TKE_Niku_lay    = TKE_Niku_rem(i)    - TKE_Niku_bot(i)   * TKE_frac_top_lee(i)
      TKE_lowmode_lay = TKE_lowmode_rem(i) - TKE_lowmode_bot(i)* TKE_frac_top_lowmode(i)

      ! Actual power expended may be less than predicted if stratification is weak; adjust
      if (TKE_itide_lay + TKE_Niku_lay + TKE_lowmode_lay > (max_TKE(i,k))) then
        frac_used = (max_TKE(i,k)) / (TKE_itide_lay + TKE_Niku_lay + TKE_lowmode_lay)
        TKE_itide_lay   = frac_used * TKE_itide_lay
        TKE_Niku_lay    = frac_used * TKE_Niku_lay
        TKE_lowmode_lay = frac_used * TKE_lowmode_lay
      endif

      ! Calculate vertical flux available to bottom of layer above
      TKE_itidal_rem(i)  = TKE_itidal_rem(i)  - TKE_itide_lay
      TKE_Niku_rem(i)    = TKE_Niku_rem(i)    - TKE_Niku_lay
      TKE_lowmode_rem(i) = TKE_lowmode_rem(i) - TKE_lowmode_lay

      ! Convert power to diffusivity
      Kd_add  = TKE_to_Kd(i,k) * (TKE_itide_lay + TKE_Niku_lay + TKE_lowmode_lay)

      if (Kd_max >= 0.0) Kd_add = min(Kd_add, Kd_max)
      if (present(Kd_lay)) then
        Kd_lay(i,k) = Kd_lay(i,k) + Kd_add
      endif

      if (present(Kd_int)) then
        Kd_int(i,K)   = Kd_int(i,K)   + 0.5 * Kd_add
        Kd_int(i,K+1) = Kd_int(i,K+1) + 0.5 * Kd_add
      endif

      ! diagnostics
      if (associated(dd%Kd_itidal)) then
        ! If at layers, dd%Kd_itidal is just TKE_to_Kd(i,k) * TKE_itide_lay
        ! The following sets the interface diagnostics.
        Kd_add = TKE_to_Kd(i,k) * TKE_itide_lay
        if (Kd_max >= 0.0) Kd_add = min(Kd_add, Kd_max)
        if (k>1)  dd%Kd_itidal(i,j,K)   = dd%Kd_itidal(i,j,K)   + 0.5*Kd_add
        if (k<nz) dd%Kd_itidal(i,j,K+1) = dd%Kd_itidal(i,j,K+1) + 0.5*Kd_add
      endif
      if (associated(dd%Kd_Itidal_work)) &
        dd%Kd_itidal_work(i,j,k) = GV%Rho0 * TKE_itide_lay
      if (associated(dd%Fl_itidal)) dd%Fl_itidal(i,j,k) = TKE_itidal_rem(i)

      if (associated(dd%Kd_Niku)) then
        ! If at layers, dd%Kd_Niku(i,j,K) is just TKE_to_Kd(i,k) * TKE_Niku_lay
        ! The following sets the interface diagnostics.
        Kd_add = TKE_to_Kd(i,k) * TKE_Niku_lay
        if (Kd_max >= 0.0) Kd_add = min(Kd_add, Kd_max)
        if (k>1) dd%Kd_Niku(i,j,K)    = dd%Kd_Niku(i,j,K)   + 0.5*Kd_add
        if (k<nz) dd%Kd_Niku(i,j,K+1) = dd%Kd_Niku(i,j,K+1) + 0.5*Kd_add
      endif
!     if (associated(dd%Kd_Niku)) dd%Kd_Niku(i,j,K) = TKE_to_Kd(i,k) * TKE_Niku_lay
      if (associated(dd%Kd_Niku_work)) &
        dd%Kd_Niku_work(i,j,k) = GV%Rho0 * TKE_Niku_lay

      if (associated(dd%Kd_lowmode)) then
        ! If at layers, dd%Kd_lowmode is just TKE_to_Kd(i,k) * TKE_lowmode_lay
        ! The following sets the interface diagnostics.
        Kd_add = TKE_to_Kd(i,k) * TKE_lowmode_lay
        if (Kd_max >= 0.0) Kd_add = min(Kd_add, Kd_max)
        if (k>1)  dd%Kd_lowmode(i,j,K)   = dd%Kd_lowmode(i,j,K)   + 0.5*Kd_add
        if (k<nz) dd%Kd_lowmode(i,j,K+1) = dd%Kd_lowmode(i,j,K+1) + 0.5*Kd_add
      endif
      if (associated(dd%Kd_lowmode_work)) &
        dd%Kd_lowmode_work(i,j,k) = GV%Rho0 * TKE_lowmode_lay
      if (associated(dd%Fl_lowmode)) dd%Fl_lowmode(i,j,k) = TKE_lowmode_rem(i)

    enddo ; enddo
  endif ! Simmons

  ! Polzin:
  if ( use_Polzin ) then
    do k=nz-1,2,-1 ; do i=is,ie
      if (max_TKE(i,k) <= 0.0) cycle
      z_from_bot(i) = z_from_bot(i) + GV%H_to_Z*h(i,j,k)
      if (CS%answers_2018) then
        if (N2_meanz(i) > 1.0e-14*US%T_to_s**2 ) then
          z_from_bot_WKB(i) = z_from_bot_WKB(i) &
              + GV%H_to_Z * h(i,j,k) * N2_lay(i,k) / N2_meanz(i)
        else ; z_from_bot_WKB(i) = 0 ; endif
      else
        if (GV%H_to_Z*h(i,j,k) * N2_lay(i,k) < (1.0e14 * htot_WKB(i)) * N2_meanz(i)) then
          z_from_bot_WKB(i) = z_from_bot_WKB(i) + &
             GV%H_to_Z*h(i,j,k) * N2_lay(i,k) / N2_meanz(i)
        endif
      endif

      ! Fraction of bottom flux predicted to reach top of this layer
      TKE_frac_top(i)     = ( Inv_int(i) * z0_polzin_scaled(i) ) / &
                            ( z0_polzin_scaled(i) + z_from_bot_WKB(i) )
      z0_psl = z0_polzin_scaled(i)*CS%Decay_scale_factor_lee
      TKE_frac_top_lee(i) = (Inv_int_lee(i) * z0_psl) / (z0_psl + z_from_bot_WKB(i))
      TKE_frac_top_lowmode(i) = ( Inv_int_low(i) * z0_polzin_scaled(i) ) / &
                            ( z0_polzin_scaled(i) + z_from_bot_WKB(i) )

      ! Actual influx at bottom of layer minus predicted outflux at top of layer to give
      ! predicted power expended
      TKE_itide_lay   = TKE_itidal_rem(i)  - TKE_itidal_bot(i) *TKE_frac_top(i)
      TKE_Niku_lay    = TKE_Niku_rem(i)    - TKE_Niku_bot(i)   * TKE_frac_top_lee(i)
      TKE_lowmode_lay = TKE_lowmode_rem(i) - TKE_lowmode_bot(i)*TKE_frac_top_lowmode(i)

      ! Actual power expended may be less than predicted if stratification is weak; adjust
      if (TKE_itide_lay + TKE_Niku_lay + TKE_lowmode_lay > (max_TKE(i,k))) then
        frac_used = (max_TKE(i,k)) / (TKE_itide_lay + TKE_Niku_lay + TKE_lowmode_lay)
        TKE_itide_lay   = frac_used * TKE_itide_lay
        TKE_Niku_lay    = frac_used * TKE_Niku_lay
        TKE_lowmode_lay = frac_used * TKE_lowmode_lay
      endif

      ! Calculate vertical flux available to bottom of layer above
      TKE_itidal_rem(i)  = TKE_itidal_rem(i)  - TKE_itide_lay
      TKE_Niku_rem(i)    = TKE_Niku_rem(i)    - TKE_Niku_lay
      TKE_lowmode_rem(i) = TKE_lowmode_rem(i) - TKE_lowmode_lay

      ! Convert power to diffusivity
      Kd_add  = TKE_to_Kd(i,k) * (TKE_itide_lay + TKE_Niku_lay + TKE_lowmode_lay)

      if (Kd_max >= 0.0) Kd_add = min(Kd_add, Kd_max)
      if (present(Kd_lay)) then
        Kd_lay(i,k) = Kd_lay(i,k) + Kd_add
      endif

      if (present(Kd_int)) then
        Kd_int(i,K)   = Kd_int(i,K)   + 0.5 * Kd_add
        Kd_int(i,K+1) = Kd_int(i,K+1) + 0.5 * Kd_add
      endif

      ! diagnostics
      if (associated(dd%Kd_itidal)) then
        ! If at layers, this is just dd%Kd_itidal(i,j,K) = TKE_to_Kd(i,k) * TKE_itide_lay
        ! The following sets the interface diagnostics.
        Kd_add = TKE_to_Kd(i,k) * TKE_itide_lay
        if (Kd_max >= 0.0) Kd_add = min(Kd_add, Kd_max)
        if (k>1)  dd%Kd_itidal(i,j,K)   = dd%Kd_itidal(i,j,K)   + 0.5*Kd_add
        if (k<nz) dd%Kd_itidal(i,j,K+1) = dd%Kd_itidal(i,j,K+1) + 0.5*Kd_add
      endif
      if (associated(dd%Kd_Itidal_work)) &
        dd%Kd_itidal_work(i,j,k) = GV%Rho0 * TKE_itide_lay
      if (associated(dd%Fl_itidal)) dd%Fl_itidal(i,j,k) = TKE_itidal_rem(i)

      if (associated(dd%Kd_Niku)) then
        ! If at layers, this is just dd%Kd_Niku(i,j,K) = TKE_to_Kd(i,k) * TKE_Niku_lay
        ! The following sets the interface diagnostics.
        Kd_add = TKE_to_Kd(i,k) * TKE_Niku_lay
        if (Kd_max >= 0.0) Kd_add = min(Kd_add, Kd_max)
        if (k>1) dd%Kd_Niku(i,j,K)    = dd%Kd_Niku(i,j,K)   + 0.5*Kd_add
        if (k<nz) dd%Kd_Niku(i,j,K+1) = dd%Kd_Niku(i,j,K+1) + 0.5*Kd_add
      endif
   !  if (associated(dd%Kd_Niku)) dd%Kd_Niku(i,j,K) = TKE_to_Kd(i,k) * TKE_Niku_lay
      if (associated(dd%Kd_Niku_work)) dd%Kd_Niku_work(i,j,k) = GV%Rho0 * TKE_Niku_lay

      if (associated(dd%Kd_lowmode)) then
        ! If at layers, dd%Kd_lowmode is just TKE_to_Kd(i,k) * TKE_lowmode_lay
        ! The following sets the interface diagnostics.
        Kd_add = TKE_to_Kd(i,k) * TKE_lowmode_lay
        if (Kd_max >= 0.0) Kd_add = min(Kd_add, Kd_max)
        if (k>1)  dd%Kd_lowmode(i,j,K)   = dd%Kd_lowmode(i,j,K)   + 0.5*Kd_add
        if (k<nz) dd%Kd_lowmode(i,j,K+1) = dd%Kd_lowmode(i,j,K+1) + 0.5*Kd_add
      endif
      if (associated(dd%Kd_lowmode_work)) &
        dd%Kd_lowmode_work(i,j,k) = GV%Rho0 * TKE_lowmode_lay
      if (associated(dd%Fl_lowmode)) dd%Fl_lowmode(i,j,k) = TKE_lowmode_rem(i)

    enddo ; enddo
  endif ! Polzin

end subroutine add_int_tide_diffusivity

!> Sets up diagnostics arrays for tidal mixing.
subroutine setup_tidal_diagnostics(G, GV, CS)
  type(ocean_grid_type),   intent(in) :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV !< The ocean's vertical grid structure
  type(tidal_mixing_cs),   pointer    :: CS !< The control structure for this module

  ! local
  integer :: isd, ied, jsd, jed, nz
  type(tidal_mixing_diags), pointer :: dd => NULL()

  isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed; nz = GV%ke
  dd => CS%dd

  if ((CS%id_Kd_itidal > 0) .or. (CS%id_Kd_Itidal_work > 0)) then
    allocate(dd%Kd_itidal(isd:ied,jsd:jed,nz+1)) ; dd%Kd_itidal(:,:,:) = 0.0
  endif
  if ((CS%id_Kd_lowmode > 0) .or. (CS%id_Kd_lowmode_work > 0)) then
    allocate(dd%Kd_lowmode(isd:ied,jsd:jed,nz+1)) ; dd%Kd_lowmode(:,:,:) = 0.0
  endif
  if ( (CS%id_Fl_itidal > 0) ) then
    allocate(dd%Fl_itidal(isd:ied,jsd:jed,nz+1)) ; dd%Fl_itidal(:,:,:) = 0.0
  endif
  if ( (CS%id_Fl_lowmode > 0) ) then
    allocate(dd%Fl_lowmode(isd:ied,jsd:jed,nz+1)) ; dd%Fl_lowmode(:,:,:) = 0.0
  endif
  if ( (CS%id_Polzin_decay_scale > 0) ) then
    allocate(dd%Polzin_decay_scale(isd:ied,jsd:jed))
    dd%Polzin_decay_scale(:,:) = 0.0
  endif
  if ( (CS%id_N2_bot > 0) ) then
    allocate(dd%N2_bot(isd:ied,jsd:jed)) ; dd%N2_bot(:,:) = 0.0
  endif
  if ( (CS%id_N2_meanz > 0) ) then
    allocate(dd%N2_meanz(isd:ied,jsd:jed)) ; dd%N2_meanz(:,:) = 0.0
  endif
  if ( (CS%id_Polzin_decay_scale_scaled > 0) ) then
    allocate(dd%Polzin_decay_scale_scaled(isd:ied,jsd:jed))
    dd%Polzin_decay_scale_scaled(:,:) = 0.0
  endif
  if ((CS%id_Kd_Niku > 0) .or. (CS%id_Kd_Niku_work > 0)) then
    allocate(dd%Kd_Niku(isd:ied,jsd:jed,nz+1)) ; dd%Kd_Niku(:,:,:) = 0.0
  endif
  if (CS%id_Kd_Niku_work > 0) then
    allocate(dd%Kd_Niku_work(isd:ied,jsd:jed,nz)) ; dd%Kd_Niku_work(:,:,:) = 0.0
  endif
  if (CS%id_Kd_Itidal_work > 0) then
    allocate(dd%Kd_Itidal_work(isd:ied,jsd:jed,nz))
    dd%Kd_Itidal_work(:,:,:) = 0.0
  endif
  if (CS%id_Kd_Lowmode_Work > 0) then
    allocate(dd%Kd_Lowmode_Work(isd:ied,jsd:jed,nz))
    dd%Kd_Lowmode_Work(:,:,:) = 0.0
  endif
  if (CS%id_TKE_itidal > 0) then
    allocate(dd%TKE_Itidal_used(isd:ied,jsd:jed)) ; dd%TKE_Itidal_used(:,:) = 0.
  endif
  ! additional diags for CVMix
  if (CS%id_N2_int > 0) then
    allocate(dd%N2_int(isd:ied,jsd:jed,nz+1)) ; dd%N2_int(:,:,:) = 0.0
  endif
  if (CS%id_Simmons_coeff > 0) then
    if (CS%CVMix_tidal_scheme .ne. SIMMONS) then
      call MOM_error(FATAL, "setup_tidal_diagnostics: Simmons_coeff diagnostics is available "//&
                            "only when CVMix_tidal_scheme is Simmons")
    endif
    allocate(dd%Simmons_coeff_2d(isd:ied,jsd:jed)) ; dd%Simmons_coeff_2d(:,:) = 0.0
  endif
  if (CS%id_vert_dep > 0) then
    allocate(dd%vert_dep_3d(isd:ied,jsd:jed,nz+1)) ; dd%vert_dep_3d(:,:,:) = 0.0
  endif
  if (CS%id_Schmittner_coeff > 0) then
    if (CS%CVMix_tidal_scheme .ne. SCHMITTNER) then
      call MOM_error(FATAL, "setup_tidal_diagnostics: Schmittner_coeff diagnostics is available "//&
                            "only when CVMix_tidal_scheme is Schmittner.")
    endif
    allocate(dd%Schmittner_coeff_3d(isd:ied,jsd:jed,nz)) ; dd%Schmittner_coeff_3d(:,:,:) = 0.0
  endif
  if (CS%id_tidal_qe_md > 0) then
    if (CS%CVMix_tidal_scheme .ne. SCHMITTNER) then
      call MOM_error(FATAL, "setup_tidal_diagnostics: tidal_qe_md diagnostics is available "//&
                            "only when CVMix_tidal_scheme is Schmittner.")
    endif
    allocate(dd%tidal_qe_md(isd:ied,jsd:jed,nz)) ; dd%tidal_qe_md(:,:,:) = 0.0
  endif
end subroutine setup_tidal_diagnostics

!> This subroutine offers up diagnostics of the tidal mixing.
subroutine post_tidal_diagnostics(G, GV, h ,CS)
  type(ocean_grid_type),    intent(in)   :: G   !< The ocean's grid structure
  type(verticalGrid_type),  intent(in)   :: GV  !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                            intent(in)   :: h   !< Layer thicknesses [H ~> m or kg m-2].
  type(tidal_mixing_cs),    pointer      :: CS  !< The control structure for this module

  ! local
  type(tidal_mixing_diags), pointer :: dd => NULL()

  dd => CS%dd

  if (CS%Int_tide_dissipation .or. CS%Lee_wave_dissipation .or. CS%Lowmode_itidal_dissipation) then
    if (CS%id_TKE_itidal  > 0) call post_data(CS%id_TKE_itidal,  dd%TKE_itidal_used, CS%diag)
    if (CS%id_TKE_leewave > 0) call post_data(CS%id_TKE_leewave, CS%TKE_Niku,        CS%diag)
    if (CS%id_Nb          > 0) call post_data(CS%id_Nb,      CS%Nb,      CS%diag)
    if (CS%id_N2_bot      > 0) call post_data(CS%id_N2_bot,  dd%N2_bot,  CS%diag)
    if (CS%id_N2_meanz    > 0) call post_data(CS%id_N2_meanz,dd%N2_meanz,CS%diag)

    if (CS%id_Fl_itidal > 0) call post_data(CS%id_Fl_itidal, dd%Fl_itidal, CS%diag)
    if (CS%id_Kd_itidal > 0) call post_data(CS%id_Kd_itidal, dd%Kd_itidal, CS%diag)
    if (CS%id_Kd_Niku   > 0) call post_data(CS%id_Kd_Niku,   dd%Kd_Niku,   CS%diag)
    if (CS%id_Kd_lowmode> 0) call post_data(CS%id_Kd_lowmode, dd%Kd_lowmode, CS%diag)
    if (CS%id_Fl_lowmode> 0) call post_data(CS%id_Fl_lowmode, dd%Fl_lowmode, CS%diag)

    if (CS%id_N2_int> 0)        call post_data(CS%id_N2_int, dd%N2_int, CS%diag)
    if (CS%id_vert_dep> 0)      call post_data(CS%id_vert_dep, dd%vert_dep_3d, CS%diag)
    if (CS%id_Simmons_coeff> 0) call post_data(CS%id_Simmons_coeff, dd%Simmons_coeff_2d, CS%diag)
    if (CS%id_Schmittner_coeff> 0) call post_data(CS%id_Schmittner_coeff, dd%Schmittner_coeff_3d, CS%diag)
    if (CS%id_tidal_qe_md> 0) call post_data(CS%id_tidal_qe_md, dd%tidal_qe_md, CS%diag)

    if (CS%id_Kd_Itidal_Work > 0) &
      call post_data(CS%id_Kd_Itidal_Work, dd%Kd_Itidal_Work, CS%diag)
    if (CS%id_Kd_Niku_Work > 0) call post_data(CS%id_Kd_Niku_Work, dd%Kd_Niku_Work, CS%diag)
    if (CS%id_Kd_Lowmode_Work > 0) &
      call post_data(CS%id_Kd_Lowmode_Work, dd%Kd_Lowmode_Work, CS%diag)

    if (CS%id_Polzin_decay_scale > 0 ) &
      call post_data(CS%id_Polzin_decay_scale, dd%Polzin_decay_scale, CS%diag)
    if (CS%id_Polzin_decay_scale_scaled > 0 ) &
      call post_data(CS%id_Polzin_decay_scale_scaled, dd%Polzin_decay_scale_scaled, CS%diag)
  endif

  if (associated(dd%Kd_itidal)) deallocate(dd%Kd_itidal)
  if (associated(dd%Kd_lowmode)) deallocate(dd%Kd_lowmode)
  if (associated(dd%Fl_itidal)) deallocate(dd%Fl_itidal)
  if (associated(dd%Fl_lowmode)) deallocate(dd%Fl_lowmode)
  if (associated(dd%Polzin_decay_scale)) deallocate(dd%Polzin_decay_scale)
  if (associated(dd%Polzin_decay_scale_scaled)) deallocate(dd%Polzin_decay_scale_scaled)
  if (associated(dd%N2_bot)) deallocate(dd%N2_bot)
  if (associated(dd%N2_meanz)) deallocate(dd%N2_meanz)
  if (associated(dd%Kd_Niku)) deallocate(dd%Kd_Niku)
  if (associated(dd%Kd_Niku_work)) deallocate(dd%Kd_Niku_work)
  if (associated(dd%Kd_Itidal_Work))  deallocate(dd%Kd_Itidal_Work)
  if (associated(dd%Kd_Lowmode_Work)) deallocate(dd%Kd_Lowmode_Work)
  if (associated(dd%TKE_itidal_used)) deallocate(dd%TKE_itidal_used)
  if (associated(dd%N2_int)) deallocate(dd%N2_int)
  if (associated(dd%vert_dep_3d)) deallocate(dd%vert_dep_3d)
  if (associated(dd%Simmons_coeff_2d)) deallocate(dd%Simmons_coeff_2d)
  if (associated(dd%Schmittner_coeff_3d)) deallocate(dd%Schmittner_coeff_3d)
  if (associated(dd%tidal_qe_md)) deallocate(dd%tidal_qe_md)

end subroutine post_tidal_diagnostics

!> This subroutine returns a zonal slice of the topographic roughness amplitudes
subroutine tidal_mixing_h_amp(h_amp, G, j, CS)
  type(ocean_grid_type),    intent(in)  :: G     !< The ocean's grid structure
  real, dimension(SZI_(G)), intent(out) :: h_amp !< The topographic roughness amplitude [Z ~> m]
  integer,                  intent(in)  :: j     !< j-index of the row to work on
  type(tidal_mixing_cs),    pointer     :: CS    !< The control structure for this module

  integer :: i

  h_amp(:) = 0.0
  if ( CS%Int_tide_dissipation .and. .not. CS%use_CVMix_tidal ) then
    do i=G%isc,G%iec
      h_amp(i) = sqrt(CS%h2(i,j))
    enddo
  endif

end subroutine tidal_mixing_h_amp

! TODO: move this subroutine to MOM_internal_tide_input module (?)
!> This subroutine read tidal energy inputs from a file.
subroutine read_tidal_energy(G, US, tidal_energy_type, tidal_energy_file, CS)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  character(len=20),       intent(in) :: tidal_energy_type !< The type of tidal energy inputs to read
  character(len=200),      intent(in) :: tidal_energy_file !< The file from which to read tidalinputs
  type(tidal_mixing_cs),   pointer    :: CS   !< The control structure for this module
  ! local
  integer :: i, j, isd, ied, jsd, jed
  real, allocatable, dimension(:,:) :: tidal_energy_flux_2d ! input tidal energy flux at T-grid points [W m-2]

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  select case (uppercase(tidal_energy_type(1:4)))
  case ('JAYN') ! Jayne 2009
    if (.not. allocated(CS%tidal_qe_2d)) allocate(CS%tidal_qe_2d(isd:ied,jsd:jed))
    allocate(tidal_energy_flux_2d(isd:ied,jsd:jed))
    call MOM_read_data(tidal_energy_file,'wave_dissipation',tidal_energy_flux_2d, G%domain)
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      CS%tidal_qe_2d(i,j) = CS%Gamma_itides * tidal_energy_flux_2d(i,j)
    enddo ; enddo
    deallocate(tidal_energy_flux_2d)
  case ('ER03') ! Egbert & Ray 2003
    call read_tidal_constituents(G, US, tidal_energy_file, CS)
  case default
    call MOM_error(FATAL, "read_tidal_energy: Unknown tidal energy file type.")
  end select

end subroutine read_tidal_energy

!> This subroutine reads tidal input energy from a file by constituent.
subroutine read_tidal_constituents(G, US, tidal_energy_file, CS)
  type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in) :: US   !< A dimensional unit scaling type
  character(len=200),    intent(in) :: tidal_energy_file !< The file from which to read tidal energy inputs
  type(tidal_mixing_cs), pointer    :: CS   !< The control structure for this module

  ! local variables
  real, parameter :: C1_3 = 1.0/3.0
  real, dimension(SZI_(G),SZJ_(G)) :: &
    tidal_qk1, &  ! qk1 coefficient used in Schmittner & Egbert
    tidal_qo1     ! qo1 coefficient used in Schmittner & Egbert
  real, allocatable, dimension(:) :: &
    z_t, &        ! depth from surface to midpoint of input layer [Z]
    z_w           ! depth from surface to top of input layer [Z]
  real, allocatable, dimension(:,:,:) :: &
    tc_m2, &      ! input lunar semidiurnal tidal energy flux [W/m^2]
    tc_s2, &      ! input solar semidiurnal tidal energy flux [W/m^2]
    tc_k1, &      ! input lunar diurnal tidal energy flux [W/m^2]
    tc_o1         ! input lunar diurnal tidal energy flux [W/m^2]
  integer, dimension(4) :: nz_in
  integer               :: k, is, ie, js, je, isd, ied, jsd, jed, i, j

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ! get number of input levels:
  call field_size(tidal_energy_file, 'z_t', nz_in)

  ! allocate local variables
  allocate(z_t(nz_in(1)), z_w(nz_in(1)) )
  allocate(tc_m2(isd:ied,jsd:jed,nz_in(1)), &
           tc_s2(isd:ied,jsd:jed,nz_in(1)), &
           tc_k1(isd:ied,jsd:jed,nz_in(1)), &
           tc_o1(isd:ied,jsd:jed,nz_in(1)) )

  ! allocate CS variables associated with 3d tidal energy dissipation
  if (.not. allocated(CS%tidal_qe_3d_in)) allocate(CS%tidal_qe_3d_in(isd:ied,jsd:jed,nz_in(1)))
  if (.not. allocated(CS%h_src))          allocate(CS%h_src(nz_in(1)))

  ! read in tidal constituents
  call MOM_read_data(tidal_energy_file, 'M2', tc_m2, G%domain)
  call MOM_read_data(tidal_energy_file, 'S2', tc_s2, G%domain)
  call MOM_read_data(tidal_energy_file, 'K1', tc_k1, G%domain)
  call MOM_read_data(tidal_energy_file, 'O1', tc_o1, G%domain)
  ! Note the hard-coded assumption that z_t and z_w in the file are in centimeters.
  call MOM_read_data(tidal_energy_file, 'z_t', z_t, scale=100.0*US%m_to_Z)
  call MOM_read_data(tidal_energy_file, 'z_w', z_w, scale=100.0*US%m_to_Z)

  do j=js,je ; do i=is,ie
    if (abs(G%geoLatT(i,j)) < 30.0) then
      tidal_qk1(i,j) = C1_3
      tidal_qo1(i,j) = C1_3
    else
      tidal_qk1(i,j) = 1.0
      tidal_qo1(i,j) = 1.0
    endif
  enddo ; enddo

  CS%tidal_qe_3d_in(:,:,:) = 0.0
  do k=1,nz_in(1)
    ! Store the input cell thickness in m for use with CVmix.
    CS%h_src(k) = US%Z_to_m*(z_t(k)-z_w(k))*2.0
    ! form tidal_qe_3d_in from weighted tidal constituents
    do j=js,je ; do i=is,ie
      if ((z_t(k) <= G%bathyT(i,j)) .and. (z_w(k) > CS%tidal_diss_lim_tc)) &
        CS%tidal_qe_3d_in(i,j,k) = C1_3*tc_m2(i,j,k) + C1_3*tc_s2(i,j,k) + &
                tidal_qk1(i,j)*tc_k1(i,j,k) + tidal_qo1(i,j)*tc_o1(i,j,k)
    enddo ; enddo
  enddo

  !open(unit=1905,file="out_1905.txt",access="APPEND")
  !do j=G%jsd,G%jed
  !  do i=isd,ied
  !    if ( i+G%idg_offset .eq. 90 .and. j+G%jdg_offset .eq. 126) then
  !      write(1905,*) "-------------------------------------------"
  !      do k=50,nz_in(1)
  !          write(1905,*) i,j,k
  !          write(1905,*) CS%tidal_qe_3d_in(i,j,k), tc_m2(i,j,k)
  !          write(1905,*) z_t(k), G%bathyT(i,j), z_w(k),CS%tidal_diss_lim_tc
  !      end do
  !    endif
  !  enddo
  !enddo
  !close(1905)

  ! test if qE is positive
  if (any(CS%tidal_qe_3d_in<0.0)) then
    call MOM_error(FATAL, "read_tidal_constituents: Negative tidal_qe_3d_in terms.")
  endif

  !! collapse 3D q*E to 2D q*E
  !CS%tidal_qe_2d(:,:) = 0.0
  !do k=1,nz_in(1) ; do j=js,je ; do i=is,ie
  !  if (z_t(k) <= G%bathyT(i,j)) &
  !    CS%tidal_qe_2d(i,j) = CS%tidal_qe_2d(i,j) + CS%tidal_qe_3d_in(i,j,k)
  !enddo ; enddo ; enddo

  ! initialize input remapping:
  call initialize_remapping(CS%remap_cs, remapping_scheme="PLM", &
                            boundary_extrapolation=.false., check_remapping=CS%debug, &
                            answers_2018=CS%remap_answers_2018)

  deallocate(tc_m2)
  deallocate(tc_s2)
  deallocate(tc_k1)
  deallocate(tc_o1)
  deallocate(z_t)
  deallocate(z_w)

end subroutine read_tidal_constituents

!> Clear pointers and deallocate memory
subroutine tidal_mixing_end(CS)
  type(tidal_mixing_cs), pointer :: CS !< This module's control structure, which
                                       !! will be deallocated in this routine.

  if (.not.associated(CS)) return

  !TODO deallocate all the dynamically allocated members here ...
  if (allocated(CS%tidal_qe_2d))    deallocate(CS%tidal_qe_2d)
  if (allocated(CS%tidal_qe_3d_in)) deallocate(CS%tidal_qe_3d_in)
  if (allocated(CS%h_src))          deallocate(CS%h_src)
  deallocate(CS%dd)
  deallocate(CS)

end subroutine tidal_mixing_end

end module MOM_tidal_mixing
