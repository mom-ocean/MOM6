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


  !real :: local_mixing_frac !< fraction of wave energy dissipated locally.
  !real :: mixing_efficiency !< The efficiency that mechanical energy dissipation translates into mixing
  !                          !! that can be parameterized by a diffusivity acting on vertical stratification.
  !real :: vert_decay_scale  !< zeta in the Simmons paper (to compute the vertical deposition function). [m]
  !real :: tidal_max_coef    !< maximum allowable tidel diffusivity. [m^2/s]

  real, pointer, dimension(:,:) :: TKE_Niku    => NULL()
  real, pointer, dimension(:,:) :: TKE_itidal  => NULL()
  real, pointer, dimension(:,:) :: Nb          => NULL()
  real, pointer, dimension(:,:) :: mask_itidal => NULL()
  real, pointer, dimension(:,:) :: h2          => NULL()
  real, pointer, dimension(:,:) :: tideamp     => NULL() ! RMS tidal amplitude (m/s)

end type tidal_mixing_cs

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
                 units="m", default=0.0)
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
  !call get_param(param_file, mdl, "USE_CVMIX_TIDAL", tidal_mixing_init, &
  !               "If true, turns on tidal mixing scheme via CVMix", &
  !               default=.false.)
  !!call openParameterBlock(param_file,'CVMIX_TIDAL')
  !call get_param(param_file, mdl, "LOCAL_MIXING_FRAC", CS%local_mixing_frac, &
  !               "Fraction of wave energy dissipated locally.", &
  !               units="nondim", default=0.33)
  !call get_param(param_file, mdl, "MIXING_EFFICIENCY", CS%mixing_efficiency, &
  !               "Gamma in Simmons, 2004", &
  !               units="nondim", default=0.20)
  !!TODO: make sure GAMMA_ITIDES (same as LOCAL_MIXING_FRAC
  !call get_param(param_file, mdl, "VERTICAL_DECAY_SCALE", CS%vert_decay_scale, &
  !               "zeta in Simmons, 2004. Used to compute the vertical deposition function", &
  !               units="m", default=500.0)
  !!TODO: make sure int_tide_decay scale (same as VERTICAL_DECAY_SCALE is removed from code).
  !call get_param(param_file, mdl, "TIDAL_MAX_COEF", CS%tidal_max_coef, &
  !               "largest acceptable value for tidal diffusivity", &
  !               units="m^2/s", default=100e-4) ! the default is 50e-4 in CVMIX, 100e-4 in POP.
  !!call closeParameterBlock(param_file)


  !if (CS%debug) print *, __FILE__, __LINE__, tidal_mixing_init

  !! Set up CVMix
  !call cvmix_init_tidal(mix_scheme            = 'Simmons',            &
  !                      efficiency            = cs%mixing_efficiency, &
  !                      vertical_decay_scale  = cs%vert_decay_scale,  &
  !                      max_coefficient       = cs%tidal_max_coef,    &
  !                      local_mixing_frac     = cs%local_mixing_frac, &
  !                      depth_cutoff          = 0.0)
  !                      
  !! TODO: read in energy

end function tidal_mixing_init


!> ....
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
