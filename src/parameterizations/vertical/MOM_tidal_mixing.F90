!> Interface to vertical tidal mixing schemes including CVMix tidal mixing.
module MOM_tidal_mixing

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator,       only : diag_ctrl, time_type, register_diag_field
use MOM_diag_mediator,       only : safe_alloc_ptr, post_data
use MOM_EOS,                 only : calculate_density
use MOM_variables,           only : thermo_var_ptrs
use MOM_error_handler,       only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_debugging,           only : hchksum
use MOM_grid,                only : ocean_grid_type
use MOM_verticalGrid,        only : verticalGrid_type
use MOM_file_parser,         only : openParameterBlock, closeParameterBlock
use MOM_file_parser,         only : get_param, log_param, log_version, param_file_type
use MOM_string_functions,    only : uppercase
use MOM_io,                  only : slasher, MOM_read_data
use cvmix_tidal,             only : cvmix_init_tidal

implicit none ; private

#include <MOM_memory.h>

public tidal_mixing_init
public calculate_cvmix_tidal
public add_int_tide_diffusivity
public tidal_mixing_end

!> Control structure including parameters for tidal mixing.
type, public :: tidal_mixing_cs
  logical :: debug = .true.

  ! Parameters
  logical :: int_tide_dissipation  ! Internal tide conversion (from barotropic) with
                                   ! the schemes of St Laurent et al (2002)/
                                   ! Simmons et al (2004)
  integer :: Int_tide_profile ! A coded integer indicating the vertical profile
                              ! for dissipation of the internal waves.  Schemes that
                              ! are currently encoded are St Laurent et al (2002) and
                              ! Polzin (2009).
  logical :: Lee_wave_dissipation ! Enable lee-wave driven mixing, following
                                  ! Nikurashin (2010), with a vertical energy
                                  ! deposition profile specified by Lee_wave_profile.
                                  ! St Laurent et al (2002) or
                                  ! Simmons et al (2004) scheme
  integer :: Lee_wave_profile ! A coded integer indicating the vertical profile
                              ! for dissipation of the lee waves.  Schemes that are
                              ! currently encoded are St Laurent et al (2002) and
                              ! Polzin (2009).
  real :: Int_tide_decay_scale ! decay scale for internal wave TKE (meter)
  real :: Mu_itides          ! efficiency for conversion of dissipation
                             ! to potential energy (nondimensional)
  real :: Gamma_itides       ! fraction of local dissipation (nondimensional)
  real :: Gamma_lee          ! fraction of local dissipation for lee waves
                             ! (Nikurashin's energy input) (nondimensional)
  real :: Decay_scale_factor_lee ! Scaling factor for the decay scale of lee
                             ! wave energy dissipation (nondimensional)
  real :: min_zbot_itides    ! minimum depth for internal tide conversion (meter)
  logical :: Lowmode_itidal_dissipation ! Internal tide conversion (from low modes) with
                                        ! the schemes of St Laurent et al (2002)/
                                        ! Simmons et al (2004) !BDM
  real :: Nu_Polzin      ! The non-dimensional constant used in Polzin form of
                         ! the vertical scale of decay of tidal dissipation
  real :: Nbotref_Polzin ! Reference value for the buoyancy frequency at the
                         ! ocean bottom used in Polzin formulation of the
                         ! vertical scale of decay of tidal dissipation (1/s)
  real :: Polzin_decay_scale_factor ! Scaling factor for the decay length scale
                                    ! of the tidal dissipation profile in Polzin
                                    ! (nondimensional)
  real :: Polzin_decay_scale_max_factor  ! The decay length scale of tidal
                         ! dissipation profile in Polzin formulation should not
                         ! exceed Polzin_decay_scale_max_factor * depth of the
                         ! ocean (nondimensional).
  real :: Polzin_min_decay_scale ! minimum decay scale of the tidal dissipation
                                 ! profile in Polzin formulation (meter)
  real :: TKE_itide_max       ! maximum internal tide conversion (W m-2)
                              ! available to mix above the BBL
  real :: utide               ! constant tidal amplitude (m s-1) used if
                              ! tidal amplitude file is not present
  real :: kappa_itides        ! topographic wavenumber and non-dimensional scaling
  real :: kappa_h2_factor     ! factor for the product of wavenumber * rms sgs height
  character(len=200) :: inputdir

  real :: tidal_max_coef    !< maximum allowable tidal diffusivity. [m^2/s]

  real, pointer, dimension(:,:) :: TKE_Niku    => NULL()
  real, pointer, dimension(:,:) :: TKE_itidal  => NULL()
  real, pointer, dimension(:,:) :: Nb          => NULL()
  real, pointer, dimension(:,:) :: mask_itidal => NULL()
  real, pointer, dimension(:,:) :: h2          => NULL()
  real, pointer, dimension(:,:) :: tideamp     => NULL() ! RMS tidal amplitude (m/s)

end type tidal_mixing_cs

type, public :: tidal_mixing_diags
  real, pointer, dimension(:,:,:) :: &
    Kd_itidal      => NULL(),& ! internal tide diffusivity at interfaces (m2/s)
    Fl_itidal      => NULL(),& ! vertical flux of tidal turbulent dissipation (m3/s3)
    Kd_lowmode     => NULL(),& ! internal tide diffusivity at interfaces
                               ! due to propagating low modes (m2/s) (BDM)
    Fl_lowmode     => NULL(),& ! vertical flux of tidal turbulent dissipation
                               ! due to propagating low modes (m3/s3) (BDM)
    Kd_Niku        => NULL(),& ! lee-wave diffusivity at interfaces (m2/s)
    Kd_Niku_work   => NULL(),& ! layer integrated work by lee-wave driven mixing (W/m2)
    Kd_Itidal_Work => NULL(),& ! layer integrated work by int tide driven mixing (W/m2)
    Kd_Lowmode_Work=> NULL()   ! layer integrated work by low mode driven mixing (W/m2) BDM

  real, pointer, dimension(:,:) :: &
    TKE_itidal_used           => NULL(),& ! internal tide TKE input at ocean bottom (W/m2)
    N2_bot                    => NULL(),& ! bottom squared buoyancy frequency (1/s2)
    N2_meanz                  => NULL(),& ! vertically averaged buoyancy frequency (1/s2)
    Polzin_decay_scale_scaled => NULL(),& ! vertical scale of decay for tidal dissipation
    Polzin_decay_scale        => NULL()   ! vertical decay scale for tidal diss with Polzin (meter)

end type

character(len=40)         :: mdl = "MOM_tidal_mixing"     !< This module's name.
character*(20), parameter :: STLAURENT_PROFILE_STRING = "STLAURENT_02"
character*(20), parameter :: POLZIN_PROFILE_STRING = "POLZIN_09"
integer,        parameter :: STLAURENT_02 = 1
integer,        parameter :: POLZIN_09    = 2

contains

logical function tidal_mixing_init(Time, G, GV, param_file, diag, CS)

  type(time_type),          intent(in)    :: Time       !< The current time.
  type(ocean_grid_type),    intent(in)    :: G          !< Grid structure.
  type(verticalGrid_type),  intent(in)    :: GV         !< Vertical grid structure.
  type(param_file_type),    intent(in)    :: param_file !< Run-time parameter file handle
  type(diag_ctrl), target,  intent(inout) :: diag       !< Diagnostics control structure.
  type(tidal_mixing_cs),     pointer       :: CS         !< This module's control structure.

  ! Local variables
  logical :: read_tideamp
  character(len=20)  :: tmpstr
  character(len=200) :: filename, tideamp_file, h2_file, Niku_TKE_input_file
  real :: utide, zbot, hamp
  real :: Niku_scale ! local variable for scaling the Nikurashin TKE flux data
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed

! This include declares and sets the variable "version".
#include "version_variable.h"

  if (associated(CS)) then
    call MOM_error(WARNING, "tidal_mixing_init called when control structure "// &
                            "is already associated.")
    return
  endif
  allocate(CS)

  CS%debug = CS%debug.and.is_root_pe()

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ! Read parameters
  call log_version(param_file, mdl, version, &
    "Vertical Tidal Mixing Parameterization")
  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, default=".",do_not_log=.true.)
  CS%inputdir = slasher(CS%inputdir)
  call get_param(param_file, mdl, "INT_TIDE_DISSIPATION", CS%int_tide_dissipation, &
                 "If true, use an internal tidal dissipation scheme to \n"//&
                 "drive diapycnal mixing, along the lines of St. Laurent \n"//&
                 "et al. (2002) and Simmons et al. (2004).", default=.false.)
  if (CS%int_tide_dissipation) then
    call get_param(param_file, mdl, "INT_TIDE_PROFILE", tmpstr, &
                 "INT_TIDE_PROFILE selects the vertical profile of energy \n"//&
                 "dissipation with INT_TIDE_DISSIPATION. Valid values are:\n"//&
                 "\t STLAURENT_02 - Use the St. Laurent et al exponential \n"//&
                 "\t                decay profile.\n"//&
                 "\t POLZIN_09 - Use the Polzin WKB-streched algebraic \n"//&
                 "\t                decay profile.", &
                 default=STLAURENT_PROFILE_STRING)
    tmpstr = uppercase(tmpstr)
    select case (tmpstr)
      case (STLAURENT_PROFILE_STRING) ; CS%int_tide_profile = STLAURENT_02
      case (POLZIN_PROFILE_STRING) ; CS%int_tide_profile = POLZIN_09
      case default
        call MOM_error(FATAL, "set_diffusivity_init: Unrecognized setting "// &
            "#define INT_TIDE_PROFILE "//trim(tmpstr)//" found in input file.")
    end select
  endif

  call get_param(param_file, mdl, "LEE_WAVE_DISSIPATION", CS%Lee_wave_dissipation, &
                 "If true, use an lee wave driven dissipation scheme to \n"//&
                 "drive diapycnal mixing, along the lines of Nikurashin \n"//&
                 "(2010) and using the St. Laurent et al. (2002) \n"//&
                 "and Simmons et al. (2004) vertical profile", default=.false.)
  if (CS%lee_wave_dissipation) then
    call get_param(param_file, mdl, "LEE_WAVE_PROFILE", tmpstr, &
                 "LEE_WAVE_PROFILE selects the vertical profile of energy \n"//&
                 "dissipation with LEE_WAVE_DISSIPATION. Valid values are:\n"//&
                 "\t STLAURENT_02 - Use the St. Laurent et al exponential \n"//&
                 "\t                decay profile.\n"//&
                 "\t POLZIN_09 - Use the Polzin WKB-streched algebraic \n"//&
                 "\t                decay profile.", &
                 default=STLAURENT_PROFILE_STRING)
    tmpstr = uppercase(tmpstr)
    select case (tmpstr)
      case (STLAURENT_PROFILE_STRING) ; CS%lee_wave_profile = STLAURENT_02
      case (POLZIN_PROFILE_STRING) ; CS%lee_wave_profile = POLZIN_09
      case default
        call MOM_error(FATAL, "set_diffusivity_init: Unrecognized setting "// &
            "#define LEE_WAVE_PROFILE "//trim(tmpstr)//" found in input file.")
    end select
  endif

  call get_param(param_file, mdl, "INT_TIDE_LOWMODE_DISSIPATION", CS%Lowmode_itidal_dissipation, &
                 "If true, consider mixing due to breaking low modes that \n"//&
                 "have been remotely generated; as with itidal drag on the \n"//&
                 "barotropic tide, use an internal tidal dissipation scheme to \n"//&
                 "drive diapycnal mixing, along the lines of St. Laurent \n"//&
                 "et al. (2002) and Simmons et al. (2004).", default=.false.)

  if ((CS%Int_tide_dissipation .and. (CS%int_tide_profile == POLZIN_09)) .or. &
      (CS%lee_wave_dissipation .and. (CS%lee_wave_profile == POLZIN_09))) then
    call get_param(param_file, mdl, "NU_POLZIN", CS%Nu_Polzin, &
                 "When the Polzin decay profile is used, this is a \n"//&
                 "non-dimensional constant in the expression for the \n"//&
                 "vertical scale of decay for the tidal energy dissipation.", &
                 units="nondim", default=0.0697)
    call get_param(param_file, mdl, "NBOTREF_POLZIN", CS%Nbotref_Polzin, &
                 "When the Polzin decay profile is used, this is the \n"//&
                 "Rreference value of the buoyancy frequency at the ocean \n"//&
                 "bottom in the Polzin formulation for the vertical \n"//&
                 "scale of decay for the tidal energy dissipation.", &
                 units="s-1", default=9.61e-4)
    call get_param(param_file, mdl, "POLZIN_DECAY_SCALE_FACTOR", &
                 CS%Polzin_decay_scale_factor, &
                 "When the Polzin decay profile is used, this is a \n"//&
                 "scale factor for the vertical scale of decay of the tidal \n"//&
                 "energy dissipation.", default=1.0, units="nondim")
    call get_param(param_file, mdl, "POLZIN_SCALE_MAX_FACTOR", &
                 CS%Polzin_decay_scale_max_factor, &
                 "When the Polzin decay profile is used, this is a factor \n"//&
                 "to limit the vertical scale of decay of the tidal \n"//&
                 "energy dissipation to POLZIN_DECAY_SCALE_MAX_FACTOR \n"//&
                 "times the depth of the ocean.", units="nondim", default=1.0)
    call get_param(param_file, mdl, "POLZIN_MIN_DECAY_SCALE", CS%Polzin_min_decay_scale, &
                 "When the Polzin decay profile is used, this is the \n"//&
                 "minimum vertical decay scale for the vertical profile\n"//&
                 "of internal tide dissipation with the Polzin (2009) formulation", &
                 units="m", default=0.0)
  endif

  if (CS%Int_tide_dissipation .or. CS%Lee_wave_dissipation) then
    call get_param(param_file, mdl, "INT_TIDE_DECAY_SCALE", CS%Int_tide_decay_scale, &
                 "The decay scale away from the bottom for tidal TKE with \n"//&
                 "the new coding when INT_TIDE_DISSIPATION is used.", &
                 !units="m", default=0.0)
                 units="m", default=500.0)  ! TODO: confirm this new default
    call get_param(param_file, mdl, "MU_ITIDES", CS%Mu_itides, &
                 "A dimensionless turbulent mixing efficiency used with \n"//&
                 "INT_TIDE_DISSIPATION, often 0.2.", units="nondim", default=0.2)
    call get_param(param_file, mdl, "GAMMA_ITIDES", CS%Gamma_itides, &
                 "The fraction of the internal tidal energy that is \n"//&
                 "dissipated locally with INT_TIDE_DISSIPATION.  \n"//&
                 "THIS NAME COULD BE BETTER.", &
                 units="nondim", default=0.3333)
    call get_param(param_file, mdl, "MIN_ZBOT_ITIDES", CS%min_zbot_itides, &
                 "Turn off internal tidal dissipation when the total \n"//&
                 "ocean depth is less than this value.", units="m", default=0.0)

    call safe_alloc_ptr(CS%Nb,isd,ied,jsd,jed)
    call safe_alloc_ptr(CS%h2,isd,ied,jsd,jed)
    call safe_alloc_ptr(CS%TKE_itidal,isd,ied,jsd,jed)
    call safe_alloc_ptr(CS%mask_itidal,isd,ied,jsd,jed) ; CS%mask_itidal(:,:) = 1.0

    call get_param(param_file, mdl, "KAPPA_ITIDES", CS%kappa_itides, &
                 "A topographic wavenumber used with INT_TIDE_DISSIPATION. \n"//&
                 "The default is 2pi/10 km, as in St.Laurent et al. 2002.", &
                 units="m-1", default=8.e-4*atan(1.0))

    call get_param(param_file, mdl, "UTIDE", CS%utide, &
                 "The constant tidal amplitude used with INT_TIDE_DISSIPATION.", &
                 units="m s-1", default=0.0)
    call safe_alloc_ptr(CS%tideamp,is,ie,js,je) ; CS%tideamp(:,:) = CS%utide

    call get_param(param_file, mdl, "KAPPA_H2_FACTOR", CS%kappa_h2_factor, &
                 "A scaling factor for the roughness amplitude with n"//&
                 "INT_TIDE_DISSIPATION.",  units="nondim", default=1.0)
    call get_param(param_file, mdl, "TKE_ITIDE_MAX", CS%TKE_itide_max, &
                 "The maximum internal tide energy source availble to mix \n"//&
                 "above the bottom boundary layer with INT_TIDE_DISSIPATION.", &
                 units="W m-2",  default=1.0e3)

    call get_param(param_file, mdl, "READ_TIDEAMP", read_tideamp, &
                 "If true, read a file (given by TIDEAMP_FILE) containing \n"//&
                 "the tidal amplitude with INT_TIDE_DISSIPATION.", default=.false.)
    if (read_tideamp) then
      call get_param(param_file, mdl, "TIDEAMP_FILE", tideamp_file, &
                 "The path to the file containing the spatially varying \n"//&
                 "tidal amplitudes with INT_TIDE_DISSIPATION.", default="tideamp.nc")
      filename = trim(CS%inputdir) // trim(tideamp_file)
      call log_param(param_file, mdl, "INPUTDIR/TIDEAMP_FILE", filename)
      call MOM_read_data(filename, 'tideamp', CS%tideamp, G%domain, timelevel=1)
    endif

    call get_param(param_file, mdl, "H2_FILE", h2_file, &
                 "The path to the file containing the sub-grid-scale \n"//&
                 "topographic roughness amplitude with INT_TIDE_DISSIPATION.", &
                 fail_if_missing=.true.)
    filename = trim(CS%inputdir) // trim(h2_file)
    call log_param(param_file, mdl, "INPUTDIR/H2_FILE", filename)
    call MOM_read_data(filename, 'h2', CS%h2, G%domain, timelevel=1)

    do j=js,je ; do i=is,ie
      if (G%bathyT(i,j) < CS%min_zbot_itides) CS%mask_itidal(i,j) = 0.0
      CS%tideamp(i,j) = CS%tideamp(i,j) * CS%mask_itidal(i,j) * G%mask2dT(i,j)

      ! Restrict rms topo to 10 percent of column depth.
      zbot = G%bathyT(i,j)
      hamp = sqrt(CS%h2(i,j))
      hamp = min(0.1*zbot,hamp)
      CS%h2(i,j) = hamp*hamp

      utide = CS%tideamp(i,j)
      ! Compute the fixed part of internal tidal forcing; units are [kg s-2] here.
      CS%TKE_itidal(i,j) = 0.5*CS%kappa_h2_factor*GV%Rho0*&
           CS%kappa_itides*CS%h2(i,j)*utide*utide
    enddo; enddo

  endif

  if (CS%Lee_wave_dissipation) then

    call get_param(param_file, mdl, "NIKURASHIN_TKE_INPUT_FILE",Niku_TKE_input_file, &
                 "The path to the file containing the TKE input from lee \n"//&
                 "wave driven mixing. Used with LEE_WAVE_DISSIPATION.", &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, "NIKURASHIN_SCALE",Niku_scale, &
                 "A non-dimensional factor by which to scale the lee-wave \n"//&
                 "driven TKE input. Used with LEE_WAVE_DISSIPATION.", &
                 units="nondim", default=1.0)

    filename = trim(CS%inputdir) // trim(Niku_TKE_input_file)
    call log_param(param_file, mdl, "INPUTDIR/NIKURASHIN_TKE_INPUT_FILE", &
                   filename)
    call safe_alloc_ptr(CS%TKE_Niku,is,ie,js,je); CS%TKE_Niku(:,:) = 0.0
    call MOM_read_data(filename, 'TKE_input', CS%TKE_Niku, G%domain, timelevel=1 ) ! ??? timelevel -aja
    CS%TKE_Niku(:,:) = Niku_scale * CS%TKE_Niku(:,:)

    call get_param(param_file, mdl, "GAMMA_NIKURASHIN",CS%Gamma_lee, &
                 "The fraction of the lee wave energy that is dissipated \n"//&
                 "locally with LEE_WAVE_DISSIPATION.", units="nondim", &
                 default=0.3333)
    call get_param(param_file, mdl, "DECAY_SCALE_FACTOR_LEE",CS%Decay_scale_factor_lee, &
                 "Scaling for the vertical decay scaleof the local \n"//&
                 "dissipation of lee waves dissipation.", units="nondim", &
                 default=1.0)
  else
    CS%Decay_scale_factor_lee = -9.e99 ! This should never be used if CS%Lee_wave_dissipation = False
  endif

  call get_param(param_file, mdl, "USE_CVMIX_TIDAL", tidal_mixing_init, &
                 "If true, turns on tidal mixing scheme via CVMix", &
                 default=.false.)

  if (tidal_mixing_init) then 

    ! Read in CVMix params
    call openParameterBlock(param_file,'CVMIX_TIDAL')
    call get_param(param_file, mdl, "TIDAL_MAX_COEF", CS%tidal_max_coef, &
                   "largest acceptable value for tidal diffusivity", &
                   units="m^2/s", default=100e-4) ! the default is 50e-4 in CVMIX, 100e-4 in POP.
    call closeParameterBlock(param_file)


    ! Set up CVMix
    call cvmix_init_tidal(mix_scheme            = 'Simmons',                &
                          efficiency            = CS%Mu_itides,             &
                          vertical_decay_scale  = cs%int_tide_decay_scale,  &
                          max_coefficient       = cs%tidal_max_coef,        &
                          local_mixing_frac     = cs%Gamma_itides,          &
                          depth_cutoff          = 0.0)
                          
  endif ! cvmix on  

end function tidal_mixing_init

  !> This subroutine adds the effect of internal-tide-driven mixing to the layer diffusivities.
  !! The mechanisms considered are (1) local dissipation of internal waves generated by the
  !! barotropic flow ("itidal"), (2) local dissipation of internal waves generated by the propagating
  !! low modes (rays) of the internal tide ("lowmode"), and (3) local dissipation of internal lee waves.
  !! Will eventually need to add diffusivity due to other wave-breaking processes (e.g. Bottom friction,
  !! Froude-number-depending breaking, PSI, etc.).
subroutine add_int_tide_diffusivity(h, N2_bot, j, TKE_to_Kd, max_TKE, G, GV, CS, &
                                    dd, N2_lay, Kd, Kd_int, Kd_max)
  type(ocean_grid_type),                    intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                  intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h    !< Layer thicknesses, in H (usually m or kg m-2)
  real, dimension(SZI_(G)),                 intent(in)    :: N2_bot
  real, dimension(SZI_(G),SZK_(G)),         intent(in)    :: N2_lay
  integer,                                  intent(in)    :: j
  real, dimension(SZI_(G),SZK_(G)),         intent(in)    :: TKE_to_Kd, max_TKE
  type(tidal_mixing_cs),                    pointer       :: CS
  type(tidal_mixing_diags),                 intent(inout) :: dd
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: Kd
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), optional, intent(inout) :: Kd_int
  real,                                     intent(inout) :: Kd_max


  ! This subroutine adds the effect of internal-tide-driven mixing to the layer diffusivities.
  ! The mechanisms considered are (1) local dissipation of internal waves generated by the
  ! barotropic flow ("itidal"), (2) local dissipation of internal waves generated by the propagating
  ! low modes (rays) of the internal tide ("lowmode"), and (3) local dissipation of internal lee waves.
  ! Will eventually need to add diffusivity due to other wave-breaking processes (e.g. Bottom friction,
  ! Froude-number-depending breaking, PSI, etc.).

  real, dimension(SZI_(G)) :: &
    htot,             & ! total thickness above or below a layer, or the
                        ! integrated thickness in the BBL (meter)
    htot_WKB,         & ! distance from top to bottom (meter) WKB scaled
    TKE_itidal_bot,   & ! internal tide TKE at ocean bottom (m3/s3)
    TKE_Niku_bot,     & ! lee-wave TKE at ocean bottom (m3/s3)
    TKE_lowmode_bot,  & ! internal tide TKE at ocean bottom lost from all remote low modes (m3/s3) (BDM)
    Inv_int,          & ! inverse of TKE decay for int tide over the depth of the ocean (nondim)
    Inv_int_lee,      & ! inverse of TKE decay for lee waves over the depth of the ocean (nondim)
    Inv_int_low,      & ! inverse of TKE decay for low modes over the depth of the ocean (nondim) (BDM)
    z0_Polzin,        & ! TKE decay scale in Polzin formulation (meter)
    z0_Polzin_scaled, & ! TKE decay scale in Polzin formulation (meter)
                        ! multiplied by N2_bot/N2_meanz to be coherent with the WKB scaled z
                        ! z*=int(N2/N2_bot) * N2_bot/N2_meanz = int(N2/N2_meanz)
                        ! z0_Polzin_scaled = z0_Polzin * N2_bot/N2_meanz
    N2_meanz,         & ! vertically averaged squared buoyancy frequency (1/s2) for WKB scaling
    TKE_itidal_rem,   & ! remaining internal tide TKE (from barotropic source)
    TKE_Niku_rem,     & ! remaining lee-wave TKE
    TKE_lowmode_rem,  & ! remaining internal tide TKE (from propagating low mode source) (BDM)
    TKE_frac_top,     & ! fraction of bottom TKE that should appear at top of a layer (nondim)
    TKE_frac_top_lee, & ! fraction of bottom TKE that should appear at top of a layer (nondim)
    TKE_frac_top_lowmode, &
                        ! fraction of bottom TKE that should appear at top of a layer (nondim) (BDM)
    z_from_bot,       & ! distance from bottom (meter)
    z_from_bot_WKB      ! distance from bottom (meter), WKB scaled

  real :: I_rho0        ! 1 / RHO0, (m3/kg)
  real :: Kd_add        ! diffusivity to add in a layer (m2/sec)
  real :: TKE_itide_lay ! internal tide TKE imparted to a layer (from barotropic) (m3/s3)
  real :: TKE_Niku_lay  ! lee-wave TKE imparted to a layer (m3/s3)
  real :: TKE_lowmode_lay ! internal tide TKE imparted to a layer (from low mode) (m3/s3) (BDM)
  real :: frac_used     ! fraction of TKE that can be used in a layer (nondim)
  real :: Izeta         ! inverse of TKE decay scale (1/meter)
  real :: Izeta_lee     ! inverse of TKE decay scale for lee waves (1/meter)
  real :: z0_psl        ! temporary variable with units of meter
  real :: TKE_lowmode_tot ! TKE from all low modes (W/m2) (BDM)

  logical :: use_Polzin, use_Simmons
  integer :: i, k, is, ie, nz
  integer :: a, fr, m
  is = G%isc ; ie = G%iec ; nz = G%ke

  if (.not.(CS%Int_tide_dissipation .or. CS%Lee_wave_dissipation)) return

  do i=is,ie ; htot(i) = 0.0 ; Inv_int(i) = 0.0 ; Inv_int_lee(i) = 0.0 ; Inv_int_low(i) = 0.0 ;enddo
  do k=1,nz ; do i=is,ie
    htot(i) = htot(i) + GV%H_to_m*h(i,j,k)
  enddo ; enddo

  I_Rho0 = 1.0/GV%Rho0

  use_Polzin = ((CS%Int_tide_dissipation .and. (CS%int_tide_profile == POLZIN_09)) .or. &
                (CS%lee_wave_dissipation .and. (CS%lee_wave_profile == POLZIN_09)) .or. &
                (CS%Lowmode_itidal_dissipation .and. (CS%int_tide_profile == POLZIN_09)))
  use_Simmons = ((CS%Int_tide_dissipation .and. (CS%int_tide_profile == STLAURENT_02)) .or. &
                 (CS%lee_wave_dissipation .and. (CS%lee_wave_profile == STLAURENT_02)) .or. &
                 (CS%Lowmode_itidal_dissipation .and. (CS%int_tide_profile == STLAURENT_02)))

  ! Calculate parameters for vertical structure of dissipation
  ! Simmons:
  if ( use_Simmons ) then
    Izeta = 1.0 / max(CS%Int_tide_decay_scale, GV%H_subroundoff*GV%H_to_m)
    Izeta_lee = 1.0 / max(CS%Int_tide_decay_scale*CS%Decay_scale_factor_lee, &
                          GV%H_subroundoff*GV%H_to_m)
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
      z_from_bot(i) = GV%H_to_m*h(i,j,nz)
    enddo
  endif ! Simmons

  ! Polzin:
  if ( use_Polzin ) then
    ! WKB scaling of the vertical coordinate
    do i=is,ie ; N2_meanz(i)=0.0 ; enddo
    do k=1,nz ; do i=is,ie
      N2_meanz(i) = N2_meanz(i) + N2_lay(i,k)*GV%H_to_m*h(i,j,k)
    enddo ; enddo
    do i=is,ie
      N2_meanz(i) = N2_meanz(i) / (htot(i) + GV%H_subroundoff*GV%H_to_m)
      if (associated(dd%N2_meanz))  dd%N2_meanz(i,j) = N2_meanz(i)
    enddo

    ! WKB scaled z*(z=H) z* at the surface using the modified Polzin WKB scaling
    do i=is,ie ; htot_WKB(i) = htot(i) ; enddo
!    do i=is,ie ; htot_WKB(i) = 0.0 ; enddo
!    do k=1,nz ; do i=is,ie
!      htot_WKB(i) = htot_WKB(i) + GV%H_to_m*h(i,j,k)*N2_lay(i,k) / N2_meanz(i)
!    enddo ; enddo
    ! htot_WKB(i) = htot(i) ! Nearly equivalent and simpler

    do i=is,ie
      CS%Nb(i,j) = sqrt(N2_bot(i))
      if ((CS%tideamp(i,j) > 0.0) .and. &
          (CS%kappa_itides**2 * CS%h2(i,j) * CS%Nb(i,j)**3 > 1.0e-14) ) then
        z0_polzin(i) = CS%Polzin_decay_scale_factor * CS%Nu_Polzin * &
                       CS%Nbotref_Polzin**2 * CS%tideamp(i,j) / &
                     ( CS%kappa_itides**2 * CS%h2(i,j) * CS%Nb(i,j)**3 )
        if (z0_polzin(i) < CS%Polzin_min_decay_scale) &
          z0_polzin(i) = CS%Polzin_min_decay_scale
        if (N2_meanz(i) > 1.0e-14  ) then
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

      if (associated(dd%Polzin_decay_scale)) &
        dd%Polzin_decay_scale(i,j) = z0_polzin(i)
      if (associated(dd%Polzin_decay_scale_scaled)) &
        dd%Polzin_decay_scale_scaled(i,j) = z0_polzin_scaled(i)
      if (associated(dd%N2_bot)) dd%N2_bot(i,j) = CS%Nb(i,j)*CS%Nb(i,j)

      if ( CS%Int_tide_dissipation .and. (CS%int_tide_profile == POLZIN_09) ) then
        ! For the Polzin formulation, this if loop prevents the vertical
        ! flux of energy dissipation from having NaN values
        if (htot_WKB(i) > 1.0e-14) then
          Inv_int(i) = ( z0_polzin_scaled(i) / htot_WKB(i) ) + 1
        endif
      endif
      if ( CS%lee_wave_dissipation .and. (CS%lee_wave_profile == POLZIN_09) ) then
        ! For the Polzin formulation, this if loop prevents the vertical
        ! flux of energy dissipation from having NaN values
        if (htot_WKB(i) > 1.0e-14) then
          Inv_int_lee(i) = ( z0_polzin_scaled(i)*CS%Decay_scale_factor_lee / htot_WKB(i) ) + 1
        endif
      endif
      if ( CS%Lowmode_itidal_dissipation .and. (CS%int_tide_profile == POLZIN_09) ) then
        ! For the Polzin formulation, this if loop prevents the vertical
        ! flux of energy dissipation from having NaN values
        if (htot_WKB(i) > 1.0e-14) then
          Inv_int_low(i) = ( z0_polzin_scaled(i) / htot_WKB(i) ) + 1
        endif
      endif

      z_from_bot(i) = GV%H_to_m*h(i,j,nz)
      ! Use the new formulation for WKB scaling.  N2 is referenced to its
      ! vertical mean.
      if (N2_meanz(i) > 1.0e-14 ) then
        z_from_bot_WKB(i) = GV%H_to_m*h(i,j,nz)*N2_lay(i,nz) / N2_meanz(i)
      else ; z_from_bot_WKB(i) = 0 ; endif
    enddo
  endif  ! Polzin

  ! Calculate/get dissipation values at bottom
  ! Both Polzin and Simmons:
  do i=is,ie
    ! Dissipation of locally trapped internal tide (non-propagating high modes)
    TKE_itidal_bot(i) = min(CS%TKE_itidal(i,j)*CS%Nb(i,j),CS%TKE_itide_max)
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
      print *, "========", __FILE__, __LINE__
      call MOM_error(FATAL,"this block not supported yet. (aa)")

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
      z_from_bot(i) = z_from_bot(i) + GV%H_to_m*h(i,j,k)

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
      if (TKE_itide_lay + TKE_Niku_lay + TKE_lowmode_lay > max_TKE(i,k)) then
         frac_used = max_TKE(i,k) / (TKE_itide_lay + TKE_Niku_lay + TKE_lowmode_lay)
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
      Kd(i,j,k) = Kd(i,j,k) + Kd_add

      if (present(Kd_int)) then
        Kd_int(i,j,K)   = Kd_int(i,j,K)   + 0.5*Kd_add
        Kd_int(i,j,K+1) = Kd_int(i,j,K+1) + 0.5*Kd_add
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

    enddo ; enddo ;
  endif ! Simmons

  ! Polzin:
  if ( use_Polzin ) then
    do k=nz-1,2,-1 ; do i=is,ie
      if (max_TKE(i,k) <= 0.0) cycle
      z_from_bot(i) = z_from_bot(i) + GV%H_to_m*h(i,j,k)
      if (N2_meanz(i) > 1.0e-14 ) then
        z_from_bot_WKB(i) = z_from_bot_WKB(i) + GV%H_to_m*h(i,j,k)*N2_lay(i,k)/N2_meanz(i)
      else ; z_from_bot_WKB(i) = 0 ; endif

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
      if (TKE_itide_lay + TKE_Niku_lay + TKE_lowmode_lay > max_TKE(i,k)) then
        frac_used = max_TKE(i,k) / (TKE_itide_lay + TKE_Niku_lay + TKE_lowmode_lay)
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
      Kd(i,j,k) = Kd(i,j,k) + Kd_add

      if (present(Kd_int)) then
        Kd_int(i,j,K)   = Kd_int(i,j,K)   + 0.5*Kd_add
        Kd_int(i,j,K+1) = Kd_int(i,j,K+1) + 0.5*Kd_add
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

    enddo ; enddo;
  endif ! Polzin

end subroutine add_int_tide_diffusivity


!TODO:
subroutine calculate_cvmix_tidal()
  continue
end subroutine calculate_cvmix_tidal


!> Clear pointers and deallocate memory
subroutine tidal_mixing_end(CS)
  type(tidal_mixing_cs), pointer :: CS ! This module's control structure

  !TODO deallocate all the dynamically allocated members here ...
  deallocate(CS)
end subroutine tidal_mixing_end


end module MOM_tidal_mixing
