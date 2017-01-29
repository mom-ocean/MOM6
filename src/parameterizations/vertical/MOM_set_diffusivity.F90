module MOM_set_diffusivity
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, September 1997 - June 2007                     *
!*                                                                     *
!*    This file contains the subroutines that sets the diapycnal       *
!*  diffusivity, perhaps adding up pieces that are calculated in other *
!*  files and passed in via the vertvisc type argument.                *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, buoy, ustar, T, S, Kd, ea, eb, etc.   *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_cpu_clock,           only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,           only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,       only : diag_ctrl, time_type
use MOM_diag_mediator,       only : safe_alloc_ptr, post_data, register_diag_field
use MOM_diag_to_Z,           only : diag_to_Z_CS, register_Zint_diag, calc_Zint_diags
use MOM_debugging,           only : hchksum, uchksum, vchksum
use MOM_EOS,                 only : calculate_density, calculate_density_derivs
use MOM_error_handler,       only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_error_handler,       only : callTree_showQuery
use MOM_error_handler,       only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,         only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,        only : forcing, optics_type
use MOM_grid,                only : ocean_grid_type
use MOM_internal_tides,      only : int_tide_CS, get_lowmode_loss
use MOM_intrinsic_functions, only : invcosh
use MOM_io,                  only : slasher, vardesc, var_desc
use MOM_kappa_shear,         only : calculate_kappa_shear, kappa_shear_init, Kappa_shear_CS
use MOM_cvmix_shear,         only : calculate_cvmix_shear, cvmix_shear_init, cvmix_shear_CS
use MOM_string_functions,    only : uppercase
use MOM_thickness_diffuse,   only : vert_fill_TS
use MOM_variables,           only : thermo_var_ptrs, vertvisc_type, p3d
use MOM_verticalGrid, only : verticalGrid_type

use user_change_diffusivity, only : user_change_diff, user_change_diff_init
use user_change_diffusivity, only : user_change_diff_end, user_change_diff_CS

use fms_mod,                 only : read_data

implicit none ; private

#include <MOM_memory.h>

public set_diffusivity
public set_BBL_TKE
public set_diffusivity_init
public set_diffusivity_end

type, public :: set_diffusivity_CS ; private
  logical :: debug           ! If true, write verbose checksums for debugging.

  logical :: bulkmixedlayer  ! If true, a refined bulk mixed layer is used with
                             ! GV%nk_rho_varies variable density mixed & buffer
                             ! layers.
  real    :: FluxRi_max      ! The flux Richardson number where the stratification is
                             ! large enough that N2 > omega2.  The full expression for
                             ! the Flux Richardson number is usually
                             ! FLUX_RI_MAX*N2/(N2+OMEGA2). The default is 0.2.
  logical :: Henyey_IGW_background  ! If true, use a simplified variant of the
                             ! Henyey et al, JGR (1986) latitudinal scaling for
                             ! the background diapycnal diffusivity, which gives
                             ! a marked decrease in the diffusivity near the
                             ! equator.  The simplification here is to assume
                             ! that the in-situ stratification is the same as
                             ! the reference stratificaiton.
  logical :: Henyey_IGW_background_new ! same as Henyey_IGW_background
                             ! but incorporate the effect of
                             ! stratification on TKE dissipation,
                             !
                             ! e = f/f_0 * acosh(N/f) / acosh(N_0/f_0) * e_0
                             !
                             ! where e is the TKE dissipation, and N_0 and f_0 are the
                             ! reference buoyancy frequency and inertial frequencies respectively.
                             ! e_0 is the reference dissipation at (N_0,f_0). In the
                             ! previous version, N=N_0.
                             !
                             ! Additionally, the squared inverse relationship between
                             ! diapycnal diffusivities and stratification is included
                             !
                             ! kd = e/N^2
                             !
                             ! where kd is the diapycnal diffusivity.
                             ! This approach assumes that work done
                             ! against gravity is uniformly distributed
                             ! throughout the column. Whereas, kd=kd_0*e,
                             ! as in the original version, concentrates buoyancy
                             ! work in regions of strong stratification.

  logical :: Kd_tanh_lat_fn  ! If true, use the tanh dependence of Kd_sfc on
                             ! latitude, like GFDL CM2.1/CM2M.  There is no physical
                             ! justification for this form, and it can not be
                             ! used with Henyey_IGW_background.
  real :: Kd_tanh_lat_scale  ! A nondimensional scaling for the range of
                             ! diffusivities with Kd_tanh_lat_fn. Valid values
                             ! are in the range of -2 to 2; 0.4 reproduces CM2M.

  logical :: bottomdraglaw   ! If true, the  bottom stress is calculated with a
                             ! drag law c_drag*|u|*u.
  logical :: BBL_mixing_as_max !  If true, take the maximum of the diffusivity
                             ! from the BBL mixing and the other diffusivities.
                             ! Otherwise, diffusivities from the BBL_mixing is
                             ! added.
  logical :: use_LOTW_BBL_diffusivity ! If true, use simpler/less precise, BBL diffusivity.
  logical :: LOTW_BBL_use_omega ! If true, use simpler/less precise, BBL diffusivity.
  real    :: BBL_effic       ! efficiency with which the energy extracted
                             ! by bottom drag drives BBL diffusion (nondim)
  real    :: cdrag           ! quadratic drag coefficient (nondim)
  real    :: IMax_decay      ! inverse of a maximum decay scale for
                             ! bottom-drag driven turbulence, (1/m)

  real    :: Kd              ! interior diapycnal diffusivity (m2/s)
  real    :: Kd_min          ! minimum diapycnal diffusivity (m2/s)
  real    :: Kd_max          ! maximum increment for diapycnal diffusivity (m2/s)
                             ! Set to a negative value to have no limit.
  real    :: Kd_add          ! uniform diffusivity added everywhere without
                             ! filtering or scaling (m2/s)
  real    :: Kv              ! interior vertical viscosity (m2/s)
  real    :: Kdml            ! mixed layer diapycnal diffusivity (m2/s)
                             ! when bulkmixedlayer==.false.
  real    :: Hmix            ! mixed layer thickness (meter) when
                             ! bulkmixedlayer==.false.

  logical :: Bryan_Lewis_diffusivity ! If true, background vertical diffusivity
                                     ! uses Bryan-Lewis (1979) like tanh profile.
  real    :: Kd_Bryan_Lewis_deep     ! abyssal value of Bryan-Lewis profile (m2/s)
  real    :: Kd_Bryan_Lewis_surface  ! surface value of Bryan-Lewis profile (m2/s)
  real    :: Bryan_Lewis_depth_cent  ! center of transition depth in Bryan-Lewis (meter)
  real    :: Bryan_Lewis_width_trans ! width of transition for Bryan-Lewis (meter)

  real    :: N0_2Omega       ! ratio of the typical Buoyancy frequency to
                             ! twice the Earth's rotation period, used with the
                             ! Henyey scaling from the mixing
  real    :: N2_FLOOR_IOMEGA2 ! floor applied to N2(k) scaled by Omega^2
                              ! If =0., N2(k) is positive definite
                              ! If =1., N2(k) > Omega^2 everywhere

  type(diag_ctrl), pointer :: diag ! structure to regulate diagn output timing

  real :: Int_tide_decay_scale ! decay scale for internal wave TKE (meter)
  real :: Mu_itides          ! efficiency for conversion of dissipation
                             ! to potential energy (nondimensional)
  real :: Gamma_itides       ! fraction of local dissipation (nondimensional)
  real :: Gamma_lee          ! fraction of local dissipation for lee waves
                             ! (Nikurashin's energy input) (nondimensional)
  real :: Decay_scale_factor_lee ! Scaling factor for the decay scale of lee
                             ! wave energy dissipation (nondimensional)
  real :: min_zbot_itides    ! minimum depth for internal tide conversion (meter)
  logical :: Int_tide_dissipation  ! Internal tide conversion (from barotropic) with
                                   ! the schemes of St Laurent et al (2002)/
                                   ! Simmons et al (2004)
  logical :: Lowmode_itidal_dissipation ! Internal tide conversion (from low modes) with
                                        ! the schemes of St Laurent et al (2002)/
                                        ! Simmons et al (2004) !BDM
  integer :: Int_tide_profile ! A coded integer indicating the vertical profile
                              ! for dissipation of the internal waves.  Schemes that
                              ! are currently encoded are St Laurent et al (2002) and
                              ! Polzin (2009).
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
  logical :: Lee_wave_dissipation ! Enable lee-wave driven mixing, following
                                  ! Nikurashin (2010), with a vertical energy
                                  ! deposition profile specified by Lee_wave_profile.
                                  ! St Laurent et al (2002) or
                                  ! Simmons et al (2004) scheme
  integer :: Lee_wave_profile ! A coded integer indicating the vertical profile
                              ! for dissipation of the lee waves.  Schemes that are
                              ! currently encoded are St Laurent et al (2002) and
                              ! Polzin (2009).
  logical :: limit_dissipation ! If enabled, dissipation is limited to be larger
                               ! than the following:
  real :: dissip_min    ! Minimum dissipation (W/m3)
  real :: dissip_N0     ! Coefficient a in minimum dissipation = a+b*N (W/m3)
  real :: dissip_N1     ! Coefficient b in minimum dissipation = a+b*N (J/m3)
  real :: dissip_N2     ! Coefficient c in minimum dissipation = c*N2 (W m-3 s2)
  real :: dissip_Kd_min ! Minimum Kd (m2/s) with dissipatio Rho0*Kd_min*N^2

  real :: TKE_itide_max       ! maximum internal tide conversion (W m-2)
                              ! available to mix above the BBL
  real :: omega               ! Earth's rotation frequency (s-1)
  real :: utide               ! constant tidal amplitude (m s-1) used if
                              ! tidal amplitude file is not present
  real :: kappa_itides        ! topographic wavenumber and non-dimensional scaling
  real :: kappa_h2_factor     ! factor for the product of wavenumber * rms sgs height
  logical :: ML_radiation     ! allow a fraction of TKE available from wind work
                              ! to penetrate below mixed layer base with a vertical
                              ! decay scale determined by the minimum of
                              ! (1) The depth of the mixed layer, or
                              ! (2) An Ekman length scale.
                              ! Energy availble to drive mixing below the mixed layer is
                              ! given by E = ML_RAD_COEFF*MSTAR*USTAR**3.  Optionally, if
                              ! ML_rad_TKE_decay is true, this is further reduced by a factor
                              ! of exp(-h_ML*Idecay_len_TkE), where Idecay_len_TKE is
                              ! calculated the same way as in the mixed layer code.
                              ! The diapycnal diffusivity is KD(k) = E/(N2(k)+OMEGA2),
                              ! where N2 is the squared buoyancy frequency (s-2) and OMEGA2
                              ! is the rotation rate of the earth squared.
  real :: ML_rad_kd_max       ! Maximum diapycnal diffusivity due to turbulence
                              ! radiated from the base of the mixed layer (m2/s)
  real :: ML_rad_efold_coeff  ! non-dim coefficient to scale penetration depth
  real :: ML_rad_coeff        ! coefficient, which scales MSTAR*USTAR^3 to
                              ! obtain energy available for mixing below
                              ! mixed layer base (nondimensional)
  logical :: ML_rad_TKE_decay ! If true, apply same exponential decay
                              ! to ML_rad as applied to the other surface
                              ! sources of TKE in the mixed layer code.
  real    :: ustar_min        ! A minimum value of ustar to avoid numerical
                              ! problems (m/s).  If the value is small enough,
                              ! this parameter should not affect the solution.
  real    :: TKE_decay        ! ratio of natural Ekman depth to TKE decay scale (nondim)
  real    :: mstar            ! ratio of friction velocity cubed to
                              ! TKE input to the mixed layer (nondim)
  logical :: ML_use_omega     ! If true, use absolute rotation rate instead
                              ! of the vertical component of rotation when
                              ! setting the decay scale for mixed layer turbulence.
  real    :: ML_omega_frac    !   When setting the decay scale for turbulence, use
                              ! this fraction of the absolute rotation rate blended
                              ! with the local value of f, as f^2 ~= (1-of)*f^2 + of*4*omega^2.
  logical :: user_change_diff ! If true, call user-defined code to change diffusivity.
  logical :: useKappaShear    ! If true, use the kappa_shear module to find the
                              ! shear-driven diapycnal diffusivity.

  logical :: useCVmix         ! If true, use one of the CVMix modules to find
                              ! shear-driven diapycnal diffusivity.

  logical :: double_diffusion           ! If true, enable double-diffusive mixing.
  logical :: simple_TKE_to_Kd ! If true, uses a simple estimate of Kd/TKE that
                              ! does not rely on a layer-formulation.
  real    :: Max_Rrho_salt_fingers      ! max density ratio for salt fingering
  real    :: Max_salt_diff_salt_fingers ! max salt diffusivity for salt fingers (m2/s)
  real    :: Kv_molecular               ! molecular visc for double diff convect (m2/s)

  real, pointer, dimension(:,:) :: TKE_Niku    => NULL()
  real, pointer, dimension(:,:) :: TKE_itidal  => NULL()
  real, pointer, dimension(:,:) :: Nb          => NULL()
  real, pointer, dimension(:,:) :: mask_itidal => NULL()
  real, pointer, dimension(:,:) :: h2          => NULL()
  real, pointer, dimension(:,:) :: tideamp     => NULL() ! RMS tidal amplitude (m/s)

  character(len=200)                 :: inputdir
  type(user_change_diff_CS), pointer :: user_change_diff_CSp => NULL()
  type(diag_to_Z_CS),        pointer :: diag_to_Z_CSp        => NULL()
  type(Kappa_shear_CS),      pointer :: kappaShear_CSp       => NULL()
  type(CVMix_shear_CS),      pointer :: CVMix_Shear_CSp      => NULL()
  type(int_tide_CS),         pointer :: int_tide_CSp         => NULL()

  integer :: id_TKE_itidal  = -1
  integer :: id_TKE_leewave = -1
  integer :: id_maxTKE      = -1
  integer :: id_TKE_to_Kd   = -1

  integer :: id_Kd_itidal      = -1
  integer :: id_Kd_Niku        = -1
  integer :: id_Kd_lowmode     = -1
  integer :: id_Kd_user        = -1
  integer :: id_Kd_layer       = -1
  integer :: id_Kd_BBL         = -1
  integer :: id_Kd_BBL_z       = -1
  integer :: id_Kd_itidal_z    = -1
  integer :: id_Kd_Niku_z      = -1
  integer :: id_Kd_lowmode_z   = -1
  integer :: id_Kd_user_z      = -1
  integer :: id_Kd_Work        = -1
  integer :: id_Kd_Itidal_Work = -1
  integer :: id_Kd_Niku_Work   = -1
  integer :: id_Kd_Lowmode_Work= -1

  integer :: id_Fl_itidal                 = -1
  integer :: id_Fl_lowmode                = -1
  integer :: id_Polzin_decay_scale        = -1
  integer :: id_Polzin_decay_scale_scaled = -1

  integer :: id_Nb       = -1
  integer :: id_N2       = -1
  integer :: id_N2_z     = -1
  integer :: id_N2_bot   = -1
  integer :: id_N2_meanz = -1

  integer :: id_KT_extra   = -1
  integer :: id_KS_extra   = -1
  integer :: id_KT_extra_z = -1
  integer :: id_KS_extra_z = -1

end type set_diffusivity_CS

type diffusivity_diags
  real, pointer, dimension(:,:,:) :: &
    N2_3d          => NULL(),& ! squared buoyancy frequency at interfaces (1/s2)
    Kd_itidal      => NULL(),& ! internal tide diffusivity at interfaces (m2/s)
    Fl_itidal      => NULL(),& ! vertical flux of tidal turbulent dissipation (m3/s3)
    Kd_lowmode     => NULL(),& ! internal tide diffusivity at interfaces
                               ! due to propagating low modes (m2/s) (BDM)
    Fl_lowmode     => NULL(),& ! vertical flux of tidal turbulent dissipation
                               ! due to propagating low modes (m3/s3) (BDM)
    Kd_Niku        => NULL(),& ! lee-wave diffusivity at interfaces (m2/s)
    Kd_user        => NULL(),& ! user-added diffusivity at interfaces (m2/s)
    Kd_BBL         => NULL(),& ! BBL diffusivity at interfaces (m2/s)
    Kd_work        => NULL(),& ! layer integrated work by diapycnal mixing (W/m2)
    Kd_Niku_work   => NULL(),& ! layer integrated work by lee-wave driven mixing (W/m2)
    Kd_Itidal_Work => NULL(),& ! layer integrated work by int tide driven mixing (W/m2)
    Kd_Lowmode_Work=> NULL(),& ! layer integrated work by low mode driven mixing (W/m2) BDM
    maxTKE         => NULL(),& ! energy required to entrain to h_max (m3/s3)
    TKE_to_Kd      => NULL(),& ! conversion rate (~1.0 / (G_Earth + dRho_lay))
                               ! between TKE dissipated within a layer and Kd
                               ! in that layer, in m2 s-1 / m3 s-3 = s2 m-1
    KT_extra       => NULL(),& ! double diffusion diffusivity for temp (m2/s)
    KS_extra       => NULL()   ! double diffusion diffusivity for saln (m2/s)

  real, pointer, dimension(:,:) :: &
    TKE_itidal_used           => NULL(),& ! internal tide TKE input at ocean bottom (W/m2)
    N2_bot                    => NULL(),& ! bottom squared buoyancy frequency (1/s2)
    N2_meanz                  => NULL(),& ! vertically averaged buoyancy frequency (1/s2)
    Polzin_decay_scale_scaled => NULL(),& ! vertical scale of decay for tidal dissipation
    Polzin_decay_scale        => NULL()   ! vertical decay scale for tidal diss with Polzin (meter)

end type diffusivity_diags

character*(20), parameter :: STLAURENT_PROFILE_STRING = "STLAURENT_02"
character*(20), parameter :: POLZIN_PROFILE_STRING = "POLZIN_09"
integer,        parameter :: STLAURENT_02 = 1
integer,        parameter :: POLZIN_09    = 2

! Clocks
integer :: id_clock_kappaShear

contains

subroutine set_diffusivity(u, v, h, u_h, v_h, tv, fluxes, optics, visc, dt, &
                           G, GV, CS, Kd, Kd_int)
  type(ocean_grid_type),                  intent(in)    :: G
  type(verticalGrid_type),                intent(in)    :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in) :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in) :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h, u_h, v_h
  type(thermo_var_ptrs),                  intent(inout) :: tv  ! out is for tv%TempxPmE
  type(forcing),                          intent(in)    :: fluxes
  type(optics_type),                      pointer       :: optics
  type(vertvisc_type),                    intent(inout) :: visc
  real,                                   intent(in)    :: dt
  type(set_diffusivity_CS),               pointer       :: CS
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out) :: Kd
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), optional, intent(out) :: Kd_int

! Arguments:
!  (in)      u      - zonal velocity (m/s)
!  (in)      v      - meridional velocity (m/s)
!  (in)      h      - Layer thickness (m or kg/m2)
!  (in)      tv     - structure with pointers to thermodynamic fields
!  (in)      fluxes - structure of surface fluxes that may be used
!  (in)      visc   - structure containing vertical viscosities, bottom boundary
!                     layer properies, and related fields
!  (in)      dt     - time increment (sec)
!  (in)      G      - ocean grid structure
!  (in)      GV     - The ocean's vertical grid structure.
!  (in)      CS     - module control structure
!  (in)      j      - meridional index upon which to work
!  (out)     Kd     - diapycnal diffusivity of each layer (m2/sec)
!  (out,opt) Kd_int - diapycnal diffusivity at each interface (m2/sec)

  real, dimension(SZI_(G)) :: &
    depth, &      ! distance from surface of an interface (meter)
    N2_bot        ! bottom squared buoyancy frequency (1/s2)
  real, dimension(SZI_(G), SZJ_(G)) :: &
    Kd_sfc        ! surface value of the diffusivity (m2/s)

  type(diffusivity_diags) :: dd ! structure w/ arrays of pointers to avail diags

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    T_f, S_f      ! temperature and salinity (deg C and ppt);
                  ! massless layers filled vertically by diffusion.

  real, dimension(SZI_(G),SZK_(G)) :: &
    N2_lay, &     ! squared buoyancy frequency associated with layers (1/s2)
    maxTKE, &     ! energy required to entrain to h_max (m3/s3)
    TKE_to_Kd     ! conversion rate (~1.0 / (G_Earth + dRho_lay)) between
                  ! TKE dissipated within a layer and Kd in that layer, in
                  ! m2 s-1 / m3 s-3 = s2 m-1.

  real, dimension(SZI_(G),SZK_(G)+1) :: &
    N2_int,   &   ! squared buoyancy frequency associated at interfaces (1/s2)
    dRho_int, &   ! locally ref potential density difference across interfaces (in s-2) smg: or kg/m3?
    KT_extra, &   ! double difusion diffusivity on temperature (m2/sec)
    KS_extra      ! double difusion diffusivity on salinity (m2/sec)

  real :: I_trans       ! inverse of the transitional for Bryan-Lewis (1/m)
  real :: depth_c       ! depth of the center of a layer (meter)
  real :: I_Hmix        ! inverse of fixed mixed layer thickness (1/m)
  real :: I_Rho0        ! inverse of Boussinesq density (m3/kg)
  real :: I_x30         ! 2/acos(2) = 1/(sin(30 deg) * acosh(1/sin(30 deg)))
  real :: abs_sin       ! absolute value of sine of latitude (nondim)
  real :: atan_fn_sfc   ! surface value of Bryan-Lewis profile (nondim)
  real :: atan_fn_lay   ! value of Bryan-Lewis profile in layer middle (nondim)
  real :: I_atan_fn     ! inverse of change in Bryan-Lewis profile from surface to infinite depth (nondim)
  real :: deg_to_rad    ! factor converting degrees to radians, pi/180.
  real :: dissip        ! local variable for dissipation calculations (W/m3)
  real :: Omega2        ! squared absolute rotation rate (1/s2)
  real :: I_2Omega      ! 1/(2 Omega) (sec)
  real :: N_2Omega
  real :: N02_N2
  real :: epsilon

  logical   :: use_EOS      ! If true, compute density from T/S using equation of state.
  type(p3d) :: z_ptrs(6)    ! pointers to diagns to be interpolated into depth space
  integer   :: kb(SZI_(G))  ! The index of the lightest layer denser than the
                            ! buffer layer.
  integer   :: num_z_diags  ! number of diagns to be interpolated to depth space
  integer   :: z_ids(6)     ! id numbers of diagns to be interpolated to depth space
  logical   :: showCallTree ! If true, show the call tree.

  integer :: i, j, k, is, ie, js, je, nz
  integer :: isd, ied, jsd, jed

  real      :: kappa_fill   ! diffusivity used to fill massless layers
  real      :: dt_fill      ! timestep used to fill massless layers

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("set_diffusivity(), MOM_set_diffusivity.F90")

  if (.not.associated(CS)) call MOM_error(FATAL,"set_diffusivity: "//&
         "Module must be initialized before it is used.")

  I_Rho0     = 1.0/GV%Rho0
  kappa_fill = 1.e-3 ! m2 s-1
  dt_fill    = 7200.
  deg_to_rad = atan(1.0)/45.0 ! = PI/180
  Omega2     = CS%Omega*CS%Omega
  I_2Omega   = 0.5/CS%Omega
  epsilon    = 1.e-10

  use_EOS = associated(tv%eqn_of_state)

  if ((CS%double_diffusion) .and. &
      .not.(associated(visc%Kd_extra_T) .and. associated(visc%Kd_extra_S)) ) &
    call MOM_error(FATAL, "set_diffusivity: visc%Kd_extra_T and "//&
         "visc%Kd_extra_S must be associated when DOUBLE_DIFFUSION is true.")

  ! Set up arrays for diagnostics.

  if ((CS%id_N2 > 0) .or. (CS%id_N2_z > 0)) then
    allocate(dd%N2_3d(isd:ied,jsd:jed,nz+1)) ; dd%N2_3d(:,:,:) = 0.0
  endif
  if ((CS%id_Kd_itidal > 0) .or. (CS%id_Kd_itidal_z > 0) .or. &
      (CS%id_Kd_Itidal_work > 0)) then
    allocate(dd%Kd_itidal(isd:ied,jsd:jed,nz+1)) ; dd%Kd_itidal(:,:,:) = 0.0
  endif
  if ((CS%id_Kd_lowmode > 0) .or. (CS%id_Kd_lowmode_z > 0) .or. &
      (CS%id_Kd_lowmode_work > 0)) then
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
  if ( (CS%id_Polzin_decay_scale_scaled > 0) ) then
    allocate(dd%Polzin_decay_scale_scaled(isd:ied,jsd:jed))
    dd%Polzin_decay_scale_scaled(:,:) = 0.0
  endif
  if ( (CS%id_N2_bot > 0) ) then
    allocate(dd%N2_bot(isd:ied,jsd:jed)) ; dd%N2_bot(:,:) = 0.0
  endif
  if ( (CS%id_N2_meanz > 0) ) then
    allocate(dd%N2_meanz(isd:ied,jsd:jed)) ; dd%N2_meanz(:,:) = 0.0
  endif
  if ((CS%id_Kd_Niku > 0) .or. (CS%id_Kd_Niku_z > 0) .or. &
      (CS%id_Kd_Niku_work > 0)) then
    allocate(dd%Kd_Niku(isd:ied,jsd:jed,nz+1)) ; dd%Kd_Niku(:,:,:) = 0.0
  endif
  if ((CS%id_Kd_user > 0) .or. (CS%id_Kd_user_z > 0)) then
    allocate(dd%Kd_user(isd:ied,jsd:jed,nz+1)) ; dd%Kd_user(:,:,:) = 0.0
  endif
  if (CS%id_Kd_work > 0) then
    allocate(dd%Kd_work(isd:ied,jsd:jed,nz)) ; dd%Kd_work(:,:,:) = 0.0
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
  if (CS%id_maxTKE > 0) then
    allocate(dd%maxTKE(isd:ied,jsd:jed,nz)) ; dd%maxTKE(:,:,:) = 0.0
  endif
  if (CS%id_TKE_to_Kd > 0) then
    allocate(dd%TKE_to_Kd(isd:ied,jsd:jed,nz)) ; dd%TKE_to_Kd(:,:,:) = 0.0
  endif
  if ((CS%id_KT_extra > 0) .or. (CS%id_KT_extra_z > 0)) then
    allocate(dd%KT_extra(isd:ied,jsd:jed,nz+1)) ; dd%KT_extra(:,:,:) = 0.0
  endif
  if ((CS%id_KS_extra > 0) .or. (CS%id_KS_extra_z > 0)) then
    allocate(dd%KS_extra(isd:ied,jsd:jed,nz+1)) ; dd%KS_extra(:,:,:) = 0.0
  endif
  if ((CS%id_Kd_BBL > 0) .or. (CS%id_Kd_BBL_z > 0)) then
    allocate(dd%Kd_BBL(isd:ied,jsd:jed,nz+1)) ; dd%Kd_BBL(:,:,:) = 0.0
  endif

  ! Smooth the properties through massless layers.
  if (use_EOS) then
    if (CS%debug) then
      call hchksum(tv%T, "before vert_fill_TS tv%T",G%HI)
      call hchksum(tv%S, "before vert_fill_TS tv%S",G%HI)
      call hchksum(h*GV%H_to_m, "before vert_fill_TS h",G%HI)
    endif
    call vert_fill_TS(h, tv%T, tv%S, kappa_fill, dt_fill, T_f, S_f, G, GV)
    if (CS%debug) then
      call hchksum(tv%T, "after vert_fill_TS tv%T",G%HI)
      call hchksum(tv%S, "after vert_fill_TS tv%S",G%HI)
      call hchksum(h*GV%H_to_m, "after vert_fill_TS h",G%HI)
    endif
  endif

  if (CS%useKappaShear) then
    if (CS%debug) then
      call hchksum(u_h, "before calc_KS u_h",G%HI)
      call hchksum(v_h, "before calc_KS v_h",G%HI)
    endif
    call cpu_clock_begin(id_clock_kappaShear)
    ! Changes: visc%Kd_turb, visc%TKE_turb (not clear that TKE_turb is used as input ????)
    ! Sets visc%Kv_turb
    call calculate_kappa_shear(u_h, v_h, h, tv, fluxes%p_surf, visc%Kd_turb, visc%TKE_turb, &
                               visc%Kv_turb, dt, G, GV, CS%kappaShear_CSp)
    call cpu_clock_end(id_clock_kappaShear)
    if (CS%debug) then
      call hchksum(visc%Kd_turb, "after calc_KS visc%Kd_turb",G%HI)
      call hchksum(visc%Kv_turb, "after calc_KS visc%Kv_turb",G%HI)
      call hchksum(visc%TKE_turb, "after calc_KS visc%TKE_turb",G%HI)
    endif
    if (showCallTree) call callTree_waypoint("done with calculate_kappa_shear (set_diffusivity)")
 elseif (CS%useCVMix) then
    !NOTE{BGR}: this needs cleaned up.  Works in 1D case, not tested outside.
    call calculate_cvmix_shear(u_h, v_h, h, tv, visc%Kd_turb, visc%Kv_turb,G,GV,CS%CVMix_shear_CSp)
  elseif (associated(visc%Kv_turb)) then
    visc%Kv_turb(:,:,:) = 0. ! needed if calculate_kappa_shear is not enabled
  endif

!   Calculate the diffusivity, Kd, for each layer.  This would be
! the appropriate place to add a depth-dependent parameterization or
! another explicit parameterization of Kd.

  if (CS%Bryan_Lewis_diffusivity) then
!$OMP parallel do default(none) shared(is,ie,js,je,CS,Kd_sfc)
    do j=js,je ; do i=is,ie
      Kd_sfc(i,j) = CS%Kd_Bryan_Lewis_surface
    enddo ; enddo
  else
!$OMP parallel do default(none) shared(is,ie,js,je,CS,Kd_sfc)
    do j=js,je ; do i=is,ie
      Kd_sfc(i,j) = CS%Kd
    enddo ; enddo
  endif
  if (CS%Henyey_IGW_background) then
    I_x30 = 2.0 / invcosh(CS%N0_2Omega*2.0) ! This is evaluated at 30 deg.
!$OMP parallel do default(none) shared(is,ie,js,je,Kd_sfc,CS,G,deg_to_rad,epsilon,I_x30) &
!$OMP                          private(abs_sin)
    do j=js,je ; do i=is,ie
      abs_sin = abs(sin(G%geoLatT(i,j)*deg_to_rad))
      Kd_sfc(i,j) = max(CS%Kd_min, Kd_sfc(i,j) * &
           ((abs_sin * invcosh(CS%N0_2Omega/max(epsilon,abs_sin))) * I_x30) )
    enddo ; enddo
  elseif (CS%Kd_tanh_lat_fn) then
!$OMP parallel do default(none) shared(is,ie,js,je,Kd_sfc,CS,G)
    do j=js,je ; do i=is,ie
      !   The transition latitude and latitude range are hard-scaled here, since
      ! this is not really intended for wide-spread use, but rather for
      ! comparison with CM2M / CM2.1 settings.
      Kd_sfc(i,j) = max(CS%Kd_min, Kd_sfc(i,j) * (1.0 + &
          CS%Kd_tanh_lat_scale * 0.5*tanh((abs(G%geoLatT(i,j)) - 35.0)/5.0) ))
    enddo ; enddo
  endif

  if (CS%debug) call hchksum(Kd_sfc,"Kd_sfc",G%HI,haloshift=0)
!$OMP parallel do default(none) shared(is,ie,js,je,nz,G,GV,CS,h,tv,T_f,S_f,fluxes,dd, &
!$OMP                                  Kd,Kd_sfc,epsilon,deg_to_rad,I_2Omega,visc,    &
!$OMP                                  Kd_int,dt,u,v,Omega2)   &
!$OMP                          private(dRho_int,I_trans,atan_fn_sfc,I_atan_fn,atan_fn_lay, &
!$OMP                                  I_Hmix,depth_c,depth,N2_lay, N2_int, N2_bot,        &
!$OMP                                  I_x30,abs_sin,N_2Omega,N02_N2,KT_extra, KS_extra,   &
!$OMP                                  TKE_to_Kd,maxTKE,dissip,kb)
  do j=js,je
    ! Set up variables related to the stratification.
    call find_N2(h, tv, T_f, S_f, fluxes, j, G, GV, CS, dRho_int, N2_lay, N2_int, N2_bot)
    if (associated(dd%N2_3d)) then
      do K=1,nz+1 ; do i=is,ie ; dd%N2_3d(i,j,K) = N2_int(i,K) ; enddo ; enddo
    endif

    ! Set up the background diffusivity.
    if (CS%Bryan_Lewis_diffusivity) then
      I_trans = 1.0 / CS%Bryan_Lewis_width_trans
      atan_fn_sfc = atan(CS%Bryan_Lewis_depth_cent*I_trans)
      I_atan_fn = 1.0 / (2.0*atan(1.0) + atan_fn_sfc)
      do i=is,ie ; depth(i) = 0.0 ; enddo
      do k=1,nz ; do i=is,ie
        atan_fn_lay = atan((CS%Bryan_Lewis_depth_cent - &
                            (depth(i)+0.5*GV%H_to_m*h(i,j,k)))*I_trans)
        Kd(i,j,k) = Kd_sfc(i,j) + (CS%Kd_Bryan_Lewis_deep - Kd_sfc(i,j)) * &
                                  (atan_fn_sfc - atan_fn_lay) * I_atan_fn
        depth(i) = depth(i) + GV%H_to_m*h(i,j,k)
      enddo ; enddo
    elseif ((.not.CS%bulkmixedlayer) .and. (CS%Kd /= CS%Kdml)) then
      I_Hmix = 1.0 / CS%Hmix
      do i=is,ie ; depth(i) = 0.0 ; enddo
      do k=1,nz ; do i=is,ie
        depth_c = depth(i) + 0.5*GV%H_to_m*h(i,j,k)

        if (depth_c <= CS%Hmix) then ; Kd(i,j,k) = CS%Kdml
        elseif (depth_c >= 2.0*CS%Hmix) then ; Kd(i,j,k) = Kd_sfc(i,j)
        else
          Kd(i,j,k) = ((Kd_sfc(i,j) - CS%Kdml) * I_Hmix) * depth_c + &
                      (2.0*CS%Kdml - Kd_sfc(i,j))
        endif

        depth(i) = depth(i) + GV%H_to_m*h(i,j,k)
      enddo ; enddo
    elseif (CS%Henyey_IGW_background_new) then
      I_x30 = 2.0 / invcosh(CS%N0_2Omega*2.0) ! This is evaluated at 30 deg.
      do k=1,nz ; do i=is,ie
        abs_sin = max(epsilon,abs(sin(G%geoLatT(i,j)*deg_to_rad)))
        N_2Omega = max(abs_sin,sqrt(N2_lay(i,k))*I_2Omega)
        N02_N2 = (CS%N0_2Omega/N_2Omega)**2
        Kd(i,j,k) = max(CS%Kd_min, Kd_sfc(i,j) * &
             ((abs_sin * invcosh(N_2Omega/abs_sin)) * I_x30)*N02_N2)
      enddo ; enddo
    else
      do k=1,nz ; do i=is,ie
        Kd(i,j,k) = Kd_sfc(i,j)
      enddo ; enddo
    endif

    if (CS%double_diffusion) then
      call double_diffusion(tv, h, T_f, S_f, j, G, GV, CS, KT_extra, KS_extra)
      do K=2,nz ; do i=is,ie
        if (KS_extra(i,K) > KT_extra(i,K)) then ! salt fingering
          Kd(i,j,k-1) = Kd(i,j,k-1) + 0.5*KT_extra(i,K)
          Kd(i,j,k)   = Kd(i,j,k)   + 0.5*KT_extra(i,K)
          visc%Kd_extra_S(i,j,k) = KS_extra(i,K)-KT_extra(i,K)
          visc%Kd_extra_T(i,j,k) = 0.0
        elseif (KT_extra(i,K) > 0.0) then ! double-diffusive convection
          Kd(i,j,k-1) = Kd(i,j,k-1) + 0.5*KS_extra(i,K)
          Kd(i,j,k)   = Kd(i,j,k)   + 0.5*KS_extra(i,K)
          visc%Kd_extra_T(i,j,k) = KT_extra(i,K)-KS_extra(i,K)
          visc%Kd_extra_S(i,j,k) = 0.0
        else ! There is no double diffusion at this interface.
          visc%Kd_extra_T(i,j,k) = 0.0
          visc%Kd_extra_S(i,j,k) = 0.0
        endif
      enddo; enddo
      if (associated(dd%KT_extra)) then ; do K=1,nz+1 ; do i=is,ie
        dd%KT_extra(i,j,K) = KT_extra(i,K)
      enddo ; enddo ; endif

      if (associated(dd%KS_extra)) then ; do K=1,nz+1 ; do i=is,ie
        dd%KS_extra(i,j,K) = KS_extra(i,K)
      enddo ; enddo ; endif
    endif

  ! Add the input turbulent diffusivity.
    if (CS%useKappaShear .or. CS%useCVMix) then
      if (present(Kd_int)) then
        do K=2,nz ; do i=is,ie
          Kd_int(i,j,K) = visc%Kd_turb(i,j,K) + 0.5*(Kd(i,j,k-1) + Kd(i,j,k))
        enddo ; enddo
        do i=is,ie
          Kd_int(i,j,1) = visc%Kd_turb(i,j,1) ! This isn't actually used. It could be 0.
          Kd_int(i,j,nz+1) = 0.0
        enddo
      endif
      do k=1,nz ; do i=is,ie
        Kd(i,j,k) = Kd(i,j,k) + 0.5*(visc%Kd_turb(i,j,K) + visc%Kd_turb(i,j,K+1))
      enddo ; enddo
    else
      if (present(Kd_int)) then
        do i=is,ie
          Kd_int(i,j,1) = Kd(i,j,1) ; Kd_int(i,j,nz+1) = 0.0
        enddo
        do K=2,nz ; do i=is,ie
          Kd_int(i,j,K) = 0.5*(Kd(i,j,k-1) + Kd(i,j,k))
        enddo ; enddo
      endif
    endif

    call find_TKE_to_Kd(h, tv, dRho_int, N2_lay, j, dt, G, GV, CS, TKE_to_Kd, &
                        maxTKE, kb)
    if (associated(dd%maxTKE)) then ; do k=1,nz ; do i=is,ie
      dd%maxTKE(i,j,k) = maxTKE(i,k)
    enddo ; enddo ; endif
    if (associated(dd%TKE_to_Kd)) then ; do k=1,nz ; do i=is,ie
      dd%TKE_to_Kd(i,j,k) = TKE_to_Kd(i,k)
    enddo ; enddo ; endif

    ! Add the ML_Rad diffusivity.
    if (CS%ML_radiation) &
      call add_MLrad_diffusivity(h, fluxes, j, G, GV, CS, Kd, TKE_to_Kd, Kd_int)

    ! Add the Nikurashin and / or tidal bottom-driven mixing
    if (CS%Int_tide_dissipation .or. CS%Lee_wave_dissipation .or. CS%Lowmode_itidal_dissipation) &
      call add_int_tide_diffusivity(h, N2_bot, j, TKE_to_Kd, maxTKE, G, GV, CS, &
                                    dd, N2_lay, Kd, Kd_int)

    ! This adds the diffusion sustained by the energy extracted from the flow
    ! by the bottom drag.
    if (CS%bottomdraglaw .and. (CS%BBL_effic>0.0)) then
      if (CS%use_LOTW_BBL_diffusivity) then
        call add_LOTW_BBL_diffusivity(h, u, v, tv, fluxes, visc, j, N2_int, G, GV, CS,  &
                                      Kd, Kd_int, dd%Kd_BBL)
      else
        call add_drag_diffusivity(h, u, v, tv, fluxes, visc, j, TKE_to_Kd, &
                                  maxTKE, kb, G, GV, CS, Kd, Kd_int, dd%Kd_BBL)
      endif
    endif

    if (CS%limit_dissipation) then
      do k=2,nz-1 ; do i=is,ie
      ! This calculates the dissipation ONLY from Kd calculated in this routine
      ! dissip has units of W/m3 (kg/m3 * m2/s * 1/s2 = J/s/m3)
      !   1) a global constant,
      !   2) a dissipation proportional to N (aka Gargett) and
      !   3) dissipation corresponding to a (nearly) constant diffusivity.
        dissip = max( CS%dissip_min, &   ! Const. floor on dissip.
                      CS%dissip_N0 + CS%dissip_N1 * sqrt(N2_lay(i,k)), & ! Floor aka Gargett
                      CS%dissip_N2 * N2_lay(i,k) ) ! Floor of Kd_min*rho0/F_Ri
        Kd(i,j,k) = max( Kd(i,j,k) , &  ! Apply floor to Kd
           dissip * (CS%FluxRi_max / (GV%Rho0 * (N2_lay(i,k) + Omega2))) )
      enddo ; enddo

      if (present(Kd_int)) then ; do K=2,nz ; do i=is,ie
      ! This calculates the dissipation ONLY from Kd calculated in this routine
      ! dissip has units of W/m3 (kg/m3 * m2/s * 1/s2 = J/s/m3)
      !   1) a global constant,
      !   2) a dissipation proportional to N (aka Gargett) and
      !   3) dissipation corresponding to a (nearly) constant diffusivity.
        dissip = max( CS%dissip_min, &   ! Const. floor on dissip.
                      CS%dissip_N0 + CS%dissip_N1 * sqrt(N2_int(i,K)), & ! Floor aka Gargett
                      CS%dissip_N2 * N2_int(i,K) ) ! Floor of Kd_min*rho0/F_Ri
        Kd_int(i,j,K) = max( Kd_int(i,j,K) , &  ! Apply floor to Kd
           dissip * (CS%FluxRi_max / (GV%Rho0 * (N2_int(i,K) + Omega2))) )
      enddo ; enddo ; endif
    endif

    if (associated(dd%Kd_work)) then
      do k=1,nz ; do i=is,ie
        dd%Kd_Work(i,j,k)  = GV%Rho0 * Kd(i,j,k) * N2_lay(i,k) * &
                             GV%H_to_m*h(i,j,k)  ! Watt m-2 s or kg s-3
      enddo ; enddo
    endif
  enddo ! j-loop

  if (CS%debug) then
    call hchksum(Kd,"BBL Kd",G%HI,haloshift=0)
    if (CS%useKappaShear) call hchksum(visc%Kd_turb,"Turbulent Kd",G%HI,haloshift=0)
    if (associated(visc%kv_bbl_u) .and. associated(visc%kv_bbl_v)) then
      call uchksum(visc%kv_bbl_u,"BBL Kv_bbl_u",G%HI,haloshift=1)
      call vchksum(visc%kv_bbl_v,"BBL Kv_bbl_v",G%HI,haloshift=1)
    endif
    if (associated(visc%bbl_thick_u) .and. associated(visc%bbl_thick_v)) then
      call uchksum(visc%bbl_thick_u,"BBL bbl_thick_u",G%HI,haloshift=1)
      call vchksum(visc%bbl_thick_v,"BBL bbl_thick_v",G%HI,haloshift=1)
    endif
    if (associated(visc%Ray_u) .and. associated(visc%Ray_v)) then
      call uchksum(visc%Ray_u,"Ray_u",G%HI)
      call vchksum(visc%Ray_v,"Ray_v",G%HI)
    endif
  endif

  if (CS%Kd_add > 0.0) then
    if (present(Kd_int)) then
!$OMP parallel do default(none) shared(is,ie,js,je,nz,Kd_int,CS,Kd)
      do k=1,nz ; do j=js,je ; do i=is,ie
        Kd_int(i,j,K) = Kd_int(i,j,K) + CS%Kd_add
        Kd(i,j,k) = Kd(i,j,k) + CS%Kd_add
      enddo ; enddo ; enddo
    else
!$OMP parallel do default(none) shared(is,ie,js,je,nz,CS,Kd)
      do k=1,nz ; do j=js,je ; do i=is,ie
        Kd(i,j,k) = Kd(i,j,k) + CS%Kd_add
      enddo ; enddo ; enddo
    endif
  endif

  if (CS%user_change_diff) then
    call user_change_diff(h, tv, G, CS%user_change_diff_CSp, Kd, Kd_int, &
                          T_f, S_f, dd%Kd_user)
  endif

  if (CS%id_Kd_layer > 0) call post_data(CS%id_Kd_layer, Kd, CS%diag)

  num_z_diags = 0
  if (CS%Int_tide_dissipation .or. CS%Lee_wave_dissipation .or. CS%Lowmode_itidal_dissipation) then
    if (CS%id_TKE_itidal  > 0) call post_data(CS%id_TKE_itidal,  dd%TKE_itidal_used, CS%diag)
    if (CS%id_TKE_leewave > 0) call post_data(CS%id_TKE_leewave, CS%TKE_Niku,        CS%diag)
    if (CS%id_Nb          > 0) call post_data(CS%id_Nb,      CS%Nb,      CS%diag)
    if (CS%id_N2          > 0) call post_data(CS%id_N2,      dd%N2_3d,   CS%diag)
    if (CS%id_N2_bot      > 0) call post_data(CS%id_N2_bot,  dd%N2_bot,  CS%diag)
    if (CS%id_N2_meanz    > 0) call post_data(CS%id_N2_meanz,dd%N2_meanz,CS%diag)

    if (CS%id_Fl_itidal > 0) call post_data(CS%id_Fl_itidal, dd%Fl_itidal, CS%diag)
    if (CS%id_Kd_itidal > 0) call post_data(CS%id_Kd_itidal, dd%Kd_itidal, CS%diag)
    if (CS%id_Kd_Niku   > 0) call post_data(CS%id_Kd_Niku,   dd%Kd_Niku,   CS%diag)
    if (CS%id_Kd_lowmode> 0) call post_data(CS%id_Kd_lowmode, dd%Kd_lowmode, CS%diag)
    if (CS%id_Fl_lowmode> 0) call post_data(CS%id_Fl_lowmode, dd%Fl_lowmode, CS%diag)
    if (CS%id_Kd_user   > 0) call post_data(CS%id_Kd_user,   dd%Kd_user,   CS%diag)
    if (CS%id_Kd_Work   > 0) call post_data(CS%id_Kd_Work,   dd%Kd_Work,   CS%diag)
    if (CS%id_Kd_Itidal_Work > 0) &
      call post_data(CS%id_Kd_Itidal_Work, dd%Kd_Itidal_Work, CS%diag)
    if (CS%id_Kd_Niku_Work > 0) call post_data(CS%id_Kd_Niku_Work, dd%Kd_Niku_Work, CS%diag)
    if (CS%id_Kd_Lowmode_Work > 0) &
      call post_data(CS%id_Kd_Lowmode_Work, dd%Kd_Lowmode_Work, CS%diag)
    if (CS%id_maxTKE > 0) call post_data(CS%id_maxTKE, dd%maxTKE, CS%diag)
    if (CS%id_TKE_to_Kd > 0) call post_data(CS%id_TKE_to_Kd, dd%TKE_to_Kd, CS%diag)

    if (CS%id_Polzin_decay_scale > 0 ) &
      call post_data(CS%id_Polzin_decay_scale, dd%Polzin_decay_scale, CS%diag)
    if (CS%id_Polzin_decay_scale_scaled > 0 ) &
      call post_data(CS%id_Polzin_decay_scale_scaled, dd%Polzin_decay_scale_scaled, CS%diag)

    if (CS%id_Kd_itidal_z > 0) then
      num_z_diags        = num_z_diags + 1
      z_ids(num_z_diags) = CS%id_Kd_itidal_z
      z_ptrs(num_z_diags)%p => dd%Kd_itidal
    endif

    if (CS%id_Kd_Niku_z > 0) then
      num_z_diags        = num_z_diags + 1
      z_ids(num_z_diags) = CS%id_Kd_Niku_z
      z_ptrs(num_z_diags)%p => dd%Kd_Niku
    endif

    if (CS%id_Kd_lowmode_z > 0) then
      num_z_diags        = num_z_diags + 1
      z_ids(num_z_diags) = CS%id_Kd_lowmode_z
      z_ptrs(num_z_diags)%p => dd%Kd_lowmode
    endif

    if (CS%id_N2_z > 0) then
      num_z_diags = num_z_diags + 1
      z_ids(num_z_diags) = CS%id_N2_z
      z_ptrs(num_z_diags)%p => dd%N2_3d
    endif

    if (CS%id_Kd_user_z > 0) then
      num_z_diags        = num_z_diags + 1
      z_ids(num_z_diags) = CS%id_Kd_user_z
      z_ptrs(num_z_diags)%p => dd%Kd_user
    endif

  endif

  if (CS%id_KT_extra > 0) call post_data(CS%id_KT_extra, dd%KT_extra, CS%diag)
  if (CS%id_KS_extra > 0) call post_data(CS%id_KS_extra, dd%KS_extra, CS%diag)
  if (CS%id_Kd_BBL > 0)   call post_data(CS%id_Kd_BBL, dd%Kd_BBL, CS%diag)

  if (CS%id_KT_extra_z > 0) then
      num_z_diags = num_z_diags + 1
      z_ids(num_z_diags) = CS%id_KT_extra_z
      z_ptrs(num_z_diags)%p => dd%KT_extra
  endif

  if (CS%id_KS_extra_z > 0) then
      num_z_diags = num_z_diags + 1
      z_ids(num_z_diags) = CS%id_KS_extra_z
      z_ptrs(num_z_diags)%p => dd%KS_extra
  endif

  if (CS%id_Kd_BBL_z > 0) then
      num_z_diags = num_z_diags + 1
      z_ids(num_z_diags) = CS%id_Kd_BBL_z
      z_ptrs(num_z_diags)%p => dd%KS_extra
  endif

  if (num_z_diags > 0) &
    call calc_Zint_diags(h, z_ptrs, z_ids, num_z_diags, G, GV, CS%diag_to_Z_CSp)

  if (associated(dd%N2_3d)) deallocate(dd%N2_3d)
  if (associated(dd%Kd_itidal)) deallocate(dd%Kd_itidal)
  if (associated(dd%Kd_lowmode)) deallocate(dd%Kd_lowmode)
  if (associated(dd%Fl_itidal)) deallocate(dd%Fl_itidal)
  if (associated(dd%Fl_lowmode)) deallocate(dd%Fl_lowmode)
  if (associated(dd%Polzin_decay_scale)) deallocate(dd%Polzin_decay_scale)
  if (associated(dd%Polzin_decay_scale_scaled)) deallocate(dd%Polzin_decay_scale_scaled)
  if (associated(dd%N2_bot)) deallocate(dd%N2_bot)
  if (associated(dd%N2_meanz)) deallocate(dd%N2_meanz)
  if (associated(dd%Kd_work)) deallocate(dd%Kd_work)
  if (associated(dd%Kd_user)) deallocate(dd%Kd_user)
  if (associated(dd%Kd_Niku)) deallocate(dd%Kd_Niku)
  if (associated(dd%Kd_Niku_work)) deallocate(dd%Kd_Niku_work)
  if (associated(dd%Kd_Itidal_Work))  deallocate(dd%Kd_Itidal_Work)
  if (associated(dd%Kd_Lowmode_Work)) deallocate(dd%Kd_Lowmode_Work)
  if (associated(dd%TKE_itidal_used)) deallocate(dd%TKE_itidal_used)
  if (associated(dd%maxTKE)) deallocate(dd%maxTKE)
  if (associated(dd%TKE_to_Kd)) deallocate(dd%TKE_to_Kd)
  if (associated(dd%KT_extra)) deallocate(dd%KT_extra)
  if (associated(dd%KS_extra)) deallocate(dd%KS_extra)
  if (associated(dd%Kd_BBL)) deallocate(dd%Kd_BBL)

  if (showCallTree) call callTree_leave("set_diffusivity()")

end subroutine set_diffusivity

subroutine find_TKE_to_Kd(h, tv, dRho_int, N2_lay, j, dt, G, GV, CS, &
                          TKE_to_Kd, maxTKE, kb)
  type(ocean_grid_type),                   intent(in)    :: G
  type(verticalGrid_type),                 intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)   :: h
  type(thermo_var_ptrs),                   intent(in)    :: tv
  real, dimension(SZI_(G),SZK_(G)+1),      intent(in)    :: dRho_int
  real, dimension(SZI_(G),SZK_(G)),        intent(in)    :: N2_lay
  integer,                                 intent(in)    :: j
  real,                                    intent(in)    :: dt
  type(set_diffusivity_CS),                pointer       :: CS
  real, dimension(SZI_(G),SZK_(G)),        intent(out)   :: TKE_to_Kd, maxTKE
  integer, dimension(SZI_(G)),             intent(out)   :: kb

  real, dimension(SZI_(G),SZK_(G)) :: &
    ds_dsp1, &    ! coordinate variable (sigma-2) difference across an
                  ! interface divided by the difference across the interface
                  ! below it (nondimensional)
    dsp1_ds, &    ! inverse coordinate variable (sigma-2) difference
                  ! across an interface times the difference across the
                  ! interface above it (nondimensional)
    rho_0,   &    ! Layer potential densities relative to surface pressure (kg/m3)
    maxEnt        ! maxEnt is the maximum value of entrainment from below (with
                  ! compensating entrainment from above to keep the layer
                  ! density from changing) that will not deplete all of the
                  ! layers above or below a layer within a timestep (meter)
  real, dimension(SZI_(G)) :: &
    htot,    &    ! total thickness above or below a layer, or the
                  ! integrated thickness in the BBL (meter)
    mFkb,    &    ! total thickness in the mixed and buffer layers
                  ! times ds_dsp1 (meter)
    p_ref,   &    ! array of tv%P_Ref pressures
    Rcv_kmb, &    ! coordinate density in the lowest buffer layer
    p_0           ! An array of 0 pressures

  real :: dh_max      ! maximum amount of entrainment a layer could
                      ! undergo before entraining all fluid in the layers
                      ! above or below (meter)
  real :: dRho_lay    ! density change across a layer (kg/m3)
  real :: Omega2      ! rotation rate squared (1/s2)
  real :: G_Rho0      ! gravitation accel divided by Bouss ref density (m4 s-2 kg-1)
  real :: I_Rho0      ! inverse of Boussinesq reference density (m3/kg)
  real :: I_dt        ! 1/dt (1/sec)
  real :: H_neglect   ! negligibly small thickness (units as h)
  real :: hN2pO2      ! h * (N^2 + Omega^2), in m s-2.
  logical :: do_i(SZI_(G))

  integer :: i, k, is, ie, nz, i_rem, kmb, kb_min
  is = G%isc ; ie = G%iec ; nz = G%ke

  I_dt      = 1.0/dt
  Omega2    = CS%Omega**2
  G_Rho0    = GV%g_Earth / GV%Rho0
  H_neglect = GV%H_subroundoff
  I_Rho0    = 1.0/GV%Rho0

  ! Simple but coordinate-independent estimate of Kd/TKE
  if (CS%simple_TKE_to_Kd) then
    do k=1,nz ; do i=is,ie
      hN2pO2 = ( GV%H_to_m * h(i,j,k) ) * ( N2_lay(i,k) + Omega2 ) ! Units of m s-2.
      if (hN2pO2>0.) then
        TKE_to_Kd(i,k) = 1./ hN2pO2 ! Units of s2 m-1.
      else; TKE_to_Kd(i,k) = 0.; endif
      ! The maximum TKE conversion we allow is really a statement
      ! about the upper diffusivity we allow. Kd_max must be set.
      maxTKE(i,k) = hN2pO2 * CS%Kd_max ! Units of m3 s-3.
    enddo ; enddo
    kb(is:ie) = -1 ! kb should not be used by any code in non-layered mode -AJA
    return
  endif

  ! Determine kb - the index of the shallowest active interior layer.
  if (CS%bulkmixedlayer) then
    kmb = GV%nk_rho_varies
    do i=is,ie ; p_0(i) = 0.0 ; p_ref(i) = tv%P_Ref ; enddo
    do k=1,nz
      call calculate_density(tv%T(:,j,k),tv%S(:,j,k),p_0,rho_0(:,k),&
                             is,ie-is+1,tv%eqn_of_state)
    enddo
    call calculate_density(tv%T(:,j,kmb),tv%S(:,j,kmb),p_ref,Rcv_kmb,&
                           is,ie-is+1,tv%eqn_of_state)

    kb_min = kmb+1
    do i=is,ie
      !   Determine the next denser layer than the buffer layer in the
      ! coordinate density (sigma-2).
      do k=kmb+1,nz-1 ; if (Rcv_kmb(i) <= GV%Rlay(k)) exit ; enddo
      kb(i) = k

    !   Backtrack, in case there are massive layers above that are stable
    ! in sigma-0.
      do k=kb(i)-1,kmb+1,-1
        if (rho_0(i,kmb) > rho_0(i,k)) exit
        if (h(i,j,k)>2.0*GV%Angstrom) kb(i) = k
      enddo
    enddo

    call set_density_ratios(h, tv, kb, G, GV, CS, j, ds_dsp1, rho_0)
  else ! not bulkmixedlayer
    kb_min = 2 ; kmb = 0
    do i=is,ie ; kb(i) = 1 ; enddo
    call set_density_ratios(h, tv, kb, G, GV, CS, j, ds_dsp1)
  endif

  ! Determine maxEnt - the maximum permitted entrainment from below by each
  ! interior layer.
  do k=2,nz-1 ; do i=is,ie
    dsp1_ds(i,k) = 1.0 / ds_dsp1(i,k)
  enddo ; enddo
  do i=is,ie ; dsp1_ds(i,nz) = 0.0 ; enddo

  if (CS%bulkmixedlayer) then
    kmb = GV%nk_rho_varies
    do i=is,ie
      htot(i) = GV%H_to_m*h(i,j,kmb)
      mFkb(i) = 0.0
      if (kb(i) < nz) &
        mFkb(i) = ds_dsp1(i,kb(i)) * (GV%H_to_m*(h(i,j,kmb) - GV%Angstrom))
    enddo
    do k=1,kmb-1 ; do i=is,ie
      htot(i) = htot(i) + GV%H_to_m*h(i,j,k)
      mFkb(i) = mFkb(i) + ds_dsp1(i,k+1)*(GV%H_to_m*(h(i,j,k) - GV%Angstrom))
    enddo ; enddo
  else
    do i=is,i
      maxEnt(i,1) = 0.0 ; htot(i) = GV%H_to_m*(h(i,j,1) - GV%Angstrom)
    enddo
  endif
  do k=kb_min,nz-1 ; do i=is,ie
    if (k == kb(i)) then
      maxEnt(i,kb(i))= mFkb(i)
    elseif (k > kb(i)) then
      maxEnt(i,k) = (1.0/dsp1_ds(i,k))*(maxEnt(i,k-1) + htot(i))
!        maxEnt(i,k) = ds_dsp1(i,k)*(maxEnt(i,k-1) + htot(i)) ! BITWISE CHG
      htot(i) = htot(i) + GV%H_to_m*(h(i,j,k) - GV%Angstrom)
    endif
  enddo ; enddo

  do i=is,ie
    htot(i) = GV%H_to_m*(h(i,j,nz) - GV%Angstrom) ; maxEnt(i,nz) = 0.0
    do_i(i) = (G%mask2dT(i,j) > 0.5)
  enddo
  do k=nz-1,kb_min,-1
    i_rem = 0
    do i=is,ie ; if (do_i(i)) then
      if (k<kb(i)) then ; do_i(i) = .false. ; cycle ; endif
      i_rem = i_rem + 1  ! Count the i-rows that are still being worked on.
      maxEnt(i,k) = MIN(maxEnt(i,k),dsp1_ds(i,k+1)*maxEnt(i,k+1) + htot(i))
      htot(i) = htot(i) + GV%H_to_m*(h(i,j,k) - GV%Angstrom)
    endif ; enddo
    if (i_rem == 0) exit
  enddo ! k-loop

  ! Now set maxTKE and TKE_to_Kd.
  do i=is,ie
    maxTKE(i,1) = 0.0 ; TKE_to_Kd(i,1) = 0.0
    maxTKE(i,nz) = 0.0 ; TKE_to_Kd(i,nz) = 0.0
  enddo
  do k=2,kmb ; do i=is,ie
    maxTKE(i,k) = 0.0
    TKE_to_Kd(i,k) = 1.0 / ((N2_lay(i,k) + Omega2) * &
                            (GV%H_to_m*(h(i,j,k) + H_neglect)))
  enddo ; enddo
  do k=kmb+1,kb_min-1 ; do i=is,ie
    !   These are the properties in the deeper mixed and buffer layers, and
    ! should perhaps be revisited.
    maxTKE(i,k) = 0.0 ; TKE_to_Kd(i,k) = 0.0
  enddo ; enddo
  do k=kb_min,nz-1 ; do i=is,ie
    if (k<kb(i)) then
      maxTKE(i,k) = 0.0
      TKE_to_Kd(i,k) = 0.0
    else
      ! maxTKE is found by determining the kappa that gives maxEnt.
      ! ### This should be 1 / G_Earth * (delta rho_InSitu)
      !  kappa_max = I_dt * dRho_int(i,K+1) * maxEnt(i,k) * &
      !             (GV%H_to_m*h(i,j,k) + dh_max) / dRho_lay
      !  maxTKE(i,k) = GV%g_Earth * dRho_lay * kappa_max
      ! dRho_int should already be non-negative, so the max is redundant?
      dh_max = maxEnt(i,k) * (1.0 + dsp1_ds(i,k))
      dRho_lay = 0.5 * max(dRho_int(i,K) + dRho_int(i,K+1), 0.0)
      maxTKE(i,k) = I_dt * ((GV%g_Earth * I_Rho0) * &
          (0.5*max(dRho_int(i,K+1) + dsp1_ds(i,k)*dRho_int(i,K),0.0))) * &
                   ((GV%H_to_m*h(i,j,k) + dh_max) * maxEnt(i,k))
      TKE_to_Kd(i,k) = 1.0 / (G_Rho0 * dRho_lay + &
                              CS%Omega**2 * GV%H_to_m*(h(i,j,k) + H_neglect))
    endif
  enddo ; enddo

end subroutine find_TKE_to_Kd

subroutine find_N2(h, tv, T_f, S_f, fluxes, j, G, GV, CS, dRho_int, &
                   N2_lay, N2_int, N2_bot)
  type(ocean_grid_type),                    intent(in)   :: G
  type(verticalGrid_type),                  intent(in)   :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)   :: h
  type(thermo_var_ptrs),                    intent(in)   :: tv
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)   :: T_f, S_f
  type(forcing),                            intent(in)   :: fluxes
  integer,                                  intent(in)   :: j
  type(set_diffusivity_CS),                 pointer      :: CS
  real, dimension(SZI_(G),SZK_(G)+1),       intent(out)  :: dRho_int, N2_int
  real, dimension(SZI_(G),SZK_(G)),         intent(out)  :: N2_lay
  real, dimension(SZI_(G)),                 intent(out)  :: N2_bot

  real, dimension(SZI_(G),SZK_(G)+1) :: &
    dRho_int_unfilt, & ! unfiltered density differences across interfaces
    dRho_dT,         & ! partial derivative of density wrt temp (kg m-3 degC-1)
    dRho_dS            ! partial derivative of density wrt saln (kg m-3 PPT-1)

  real, dimension(SZI_(G)) :: &
    pres,      &  ! pressure at each interface (Pa)
    Temp_int,  &  ! temperature at each interface (degC)
    Salin_int, &  ! salinity at each interface (PPT)
    drho_bot,  &
    h_amp,     &
    hb,        &
    z_from_bot

  real :: Rml_base  ! density of the deepest variable density layer
  real :: dz_int    ! thickness associated with an interface (meter)
  real :: G_Rho0    ! gravitation acceleration divided by Bouss reference density (m4 s-2 kg-1)
  real :: H_neglect ! negligibly small thickness, in the same units as h.

  logical :: do_i(SZI_(G)), do_any
  integer :: i, k, is, ie, nz

  is = G%isc ; ie = G%iec ; nz = G%ke
  G_Rho0    = GV%g_Earth / GV%Rho0
  H_neglect = GV%H_subroundoff

  ! Find the (limited) density jump across each interface.
  do i=is,ie
    dRho_int(i,1) = 0.0 ; dRho_int(i,nz+1) = 0.0
    dRho_int_unfilt(i,1) = 0.0 ; dRho_int_unfilt(i,nz+1) = 0.0
  enddo
  if (associated(tv%eqn_of_state)) then
    if (associated(fluxes%p_surf)) then
      do i=is,ie ; pres(i) = fluxes%p_surf(i,j) ; enddo
    else
      do i=is,ie ; pres(i) = 0.0 ; enddo
    endif
    do K=2,nz
      do i=is,ie
        pres(i) = pres(i) + GV%H_to_Pa*h(i,j,k-1)
        Temp_Int(i) = 0.5 * (T_f(i,j,k) + T_f(i,j,k-1))
        Salin_Int(i) = 0.5 * (S_f(i,j,k) + S_f(i,j,k-1))
      enddo
      call calculate_density_derivs(Temp_int, Salin_int, pres, &
               dRho_dT(:,K), dRho_dS(:,K), is, ie-is+1, tv%eqn_of_state)
      do i=is,ie
        dRho_int(i,K) = max(dRho_dT(i,K)*(T_f(i,j,k) - T_f(i,j,k-1)) + &
                            dRho_dS(i,K)*(S_f(i,j,k) - S_f(i,j,k-1)), 0.0)
        dRho_int_unfilt(i,K) = max(dRho_dT(i,K)*(tv%T(i,j,k) - tv%T(i,j,k-1)) + &
                            dRho_dS(i,K)*(tv%S(i,j,k) - tv%S(i,j,k-1)), 0.0)
      enddo
    enddo
  else
    do K=2,nz ; do i=is,ie
      dRho_int(i,K) = GV%Rlay(k) - GV%Rlay(k-1)
    enddo ; enddo
  endif

  ! Set the buoyancy frequencies.
  do k=1,nz ; do i=is,ie
    N2_lay(i,k) = G_Rho0 * 0.5*(dRho_int(i,K) + dRho_int(i,K+1)) / &
                  (GV%H_to_m*(h(i,j,k) + H_neglect))
  enddo ; enddo
  do i=is,ie ; N2_int(i,1) = 0.0 ; N2_int(i,nz+1) = 0.0 ; enddo
  do K=2,nz ; do i=is,ie
    N2_int(i,K) = G_Rho0 * dRho_int(i,K) / &
                  (0.5*GV%H_to_m*(h(i,j,k-1) + h(i,j,k) + H_neglect))
  enddo ; enddo

  ! Find the bottom boundary layer stratification, and use this in the deepest layers.
  do i=is,ie
    hb(i) = 0.0 ; dRho_bot(i) = 0.0
    z_from_bot(i) = 0.5*GV%H_to_m*h(i,j,nz)
    do_i(i) = (G%mask2dT(i,j) > 0.5)

    if (CS%Int_tide_dissipation .or. CS%Lee_wave_dissipation) then
      h_amp(i) = sqrt(CS%h2(i,j)) ! for computing Nb
    else
      h_amp(i) = 0.0
    endif
  enddo

  do k=nz,2,-1
    do_any = .false.
    do i=is,ie ; if (do_i(i)) then
      dz_int = 0.5*GV%H_to_m*(h(i,j,k) + h(i,j,k-1))
      z_from_bot(i) = z_from_bot(i) + dz_int ! middle of the layer above

      hb(i) = hb(i) + dz_int
      dRho_bot(i) = dRho_bot(i) + dRho_int(i,K)

      if (z_from_bot(i) > h_amp(i)) then
        if (k>2) then
          ! Always include at least one full layer.
          hb(i) = hb(i) + 0.5*GV%H_to_m*(h(i,j,k-1) + h(i,j,k-2))
          dRho_bot(i) = dRho_bot(i) + dRho_int(i,K-1)
        endif
        do_i(i) = .false.
      else
        do_any = .true.
      endif
    endif ; enddo
    if (.not.do_any) exit
  enddo

  do i=is,ie
    if (hb(i) > 0.0) then
      N2_bot(i) = (G_Rho0 * dRho_bot(i)) / hb(i)
    else ;  N2_bot(i) = 0.0 ; endif
    z_from_bot(i) = 0.5*GV%H_to_m*h(i,j,nz)
    do_i(i) = (G%mask2dT(i,j) > 0.5)
  enddo

  do k=nz,2,-1
    do_any = .false.
    do i=is,ie ; if (do_i(i)) then
      dz_int = 0.5*GV%H_to_m*(h(i,j,k) + h(i,j,k-1))
      z_from_bot(i) = z_from_bot(i) + dz_int ! middle of the layer above

      N2_int(i,K) = N2_bot(i)
      if (k>2) N2_lay(i,k-1) = N2_bot(i)

      if (z_from_bot(i) > h_amp(i)) then
        if (k>2) N2_int(i,K-1) = N2_bot(i)
        do_i(i) = .false.
      else
        do_any = .true.
      endif
    endif ; enddo
    if (.not.do_any) exit
  enddo

  if (associated(tv%eqn_of_state)) then
    do K=1,nz+1 ; do i=is,ie
      dRho_int(i,K) = dRho_int_unfilt(i,K)
    enddo ; enddo
  endif

end subroutine find_N2


subroutine double_diffusion(tv, h, T_f, S_f, j, G, GV, CS, Kd_T_dd, Kd_S_dd)
  type(ocean_grid_type),                    intent(in)  :: G
  type(verticalGrid_type),                  intent(in)  :: GV
  type(thermo_var_ptrs),                    intent(in)  :: tv
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)  :: h
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)  :: T_f, S_f
  integer,                                  intent(in)  :: j
  type(set_diffusivity_CS),                 pointer     :: CS
  real, dimension(SZI_(G),SZK_(G)+1),       intent(out) :: Kd_T_dd, Kd_S_dd

! Arguments:
!  (in)      tv      - structure containing pointers to any available
!                      thermodynamic fields; absent fields have NULL ptrs
!  (in)      h       - layer thickness (m or kg m-2)
!  (in)      T_f     - layer temp in C with the values in massless layers
!                      filled vertically by diffusion
!  (in)      S_f     - layer salinities in PPT with values in massless layers
!                      filled vertically by diffusion
!  (in)      G       - ocean grid structure
!  (in)      GV      - The ocean's vertical grid structure.
!  (in)      CS      - module control structure
!  (in)      j       - meridional index upon which to work
!  (out)     Kd_T_dd - interface double diffusion diapycnal diffusivity for temp (m2/sec)
!  (out)     Kd_S_dd - interface double diffusion diapycnal diffusivity for saln (m2/sec)

! This subroutine sets the additional diffusivities of temperature and
! salinity due to double diffusion, using the same functional form as is
! used in MOM4.1, and taken from an NCAR technical note (###REF?) that updates
! what was in Large et al. (1994).  All the coefficients here should probably
! be made run-time variables rather than hard-coded constants.

  real, dimension(SZI_(G)) :: &
    dRho_dT,  &    ! partial derivatives of density wrt temp (kg m-3 degC-1)
    dRho_dS,  &    ! partial derivatives of density wrt saln (kg m-3 PPT-1)
    pres,     &    ! pressure at each interface (Pa)
    Temp_int, &    ! temp and saln at interfaces
    Salin_int

  real ::  alpha_dT ! density difference between layers due to temp diffs (kg/m3)
  real ::  beta_dS  ! density difference between layers due to saln diffs (kg/m3)

  real  :: Rrho    ! vertical density ratio
  real  :: diff_dd ! factor for double-diffusion
  real  :: prandtl ! flux ratio for diffusive convection regime

  real, parameter :: Rrho0  = 1.9          ! limit for double-diffusive density ratio
  real, parameter :: dsfmax = 1.e-4        ! max diffusivity in case of salt fingering
  real, parameter :: Kv_molecular = 1.5e-6 ! molecular viscosity  (m2/sec)

  integer :: i, k, is, ie, nz
  is = G%isc ; ie = G%iec ; nz = G%ke

  if (associated(tv%eqn_of_state)) then
    do i=is,ie
      pres(i) = 0.0 ; Kd_T_dd(i,1) = 0.0 ; Kd_S_dd(i,1) = 0.0
      Kd_T_dd(i,nz+1) = 0.0 ; Kd_S_dd(i,nz+1) = 0.0
    enddo
    do K=2,nz
      do i=is,ie
        pres(i) = pres(i) + GV%H_to_Pa*h(i,j,k-1)
        Temp_Int(i) = 0.5 * (T_f(i,j,k-1) + T_f(i,j,k))
        Salin_Int(i) = 0.5 * (S_f(i,j,k-1) + S_f(i,j,k))
      enddo
      call calculate_density_derivs(Temp_int, Salin_int, pres, &
             dRho_dT(:), dRho_dS(:), is, ie-is+1, tv%eqn_of_state)

      do i=is,ie
        alpha_dT = -1.0*dRho_dT(i) * (T_f(i,j,k-1) - T_f(i,j,k))
        beta_dS  = dRho_dS(i) * (S_f(i,j,k-1) - S_f(i,j,k))

        if ((alpha_dT > beta_dS) .and. (beta_dS > 0.0)) then  ! salt finger case
          Rrho = min(alpha_dT/beta_dS,Rrho0)
          diff_dd = 1.0 - ((RRho-1.0)/(RRho0-1.0))
          diff_dd = dsfmax*diff_dd*diff_dd*diff_dd
          Kd_T_dd(i,K) = 0.7*diff_dd
          Kd_S_dd(i,K) = diff_dd
        elseif ((alpha_dT < 0.) .and. (beta_dS < 0.) .and. (alpha_dT > beta_dS)) then ! diffusive convection
          Rrho = alpha_dT/beta_dS
          diff_dd = Kv_molecular*0.909*exp(4.6*exp(-0.54*(1/Rrho-1)))
          prandtl = 0.15*Rrho
          if (Rrho > 0.5) prandtl = (1.85-0.85/Rrho)*Rrho
          Kd_T_dd(i,K) = diff_dd
          Kd_S_dd(i,K) = prandtl*diff_dd
        else
          Kd_T_dd(i,K) = 0.0 ; Kd_S_dd(i,K) = 0.0
        endif
      enddo
    enddo
  endif

end subroutine double_diffusion

subroutine add_drag_diffusivity(h, u, v, tv, fluxes, visc, j, TKE_to_Kd, &
                                maxTKE, kb, G, GV, CS, Kd, Kd_int, Kd_BBL)
  type(ocean_grid_type),                     intent(in)    :: G
  type(verticalGrid_type),                   intent(in)    :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h
  type(thermo_var_ptrs),                     intent(in)    :: tv
  type(forcing),                             intent(in)    :: fluxes
  type(vertvisc_type),                       intent(in)    :: visc
  integer,                                   intent(in)    :: j
  real, dimension(SZI_(G),SZK_(G)),          intent(in)    :: TKE_to_Kd, maxTKE
  integer, dimension(SZI_(G)),               intent(in)    :: kb
  type(set_diffusivity_CS),                  pointer       :: CS
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: Kd
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Kd_int
  real, dimension(:,:,:),                    pointer       :: Kd_BBL

! This routine adds diffusion sustained by flow energy extracted by bottom drag.

  real, dimension(SZK_(G)+1) :: &
    Rint          ! coordinate density of an interface (kg/m3)
  real, dimension(SZI_(G)) :: &
    htot, &       ! total thickness above or below a layer, or the
                  ! integrated thickness in the BBL (meter)
    rho_htot, &   ! running integral with depth of density (kg/m2)
    gh_sum_top, & ! BBL value of g'h that can be supported by
                  ! the local ustar, times R0_g (kg/m2)
    Rho_top, &    ! density at top of the BBL (kg/m3)
    TKE, &        ! turbulent kinetic energy available to drive
                  ! bottom-boundary layer mixing in a layer (m3/s3)
    I2decay       ! inverse of twice the TKE decay scale (1/m)

  real    :: TKE_to_layer   ! TKE used to drive mixing in a layer (m3/s3)
  real    :: TKE_Ray        ! TKE from layer Rayleigh drag used to drive mixing in layer (m3/s3)
  real    :: TKE_here       ! TKE that goes into mixing in this layer (m3/s3)
  real    :: dRl, dRbot     ! temporaries holding density differences (kg/m3)
  real    :: cdrag_sqrt     ! square root of the drag coefficient (nondimensional)
  real    :: ustar_h        ! value of ustar at a thickness point (m/s)
  real    :: absf           ! average absolute Coriolis parameter around a thickness point (1/s)
  real    :: R0_g           ! Rho0 / G_Earth (kg s2 m-2)
  real    :: I_rho0         ! 1 / RHO0
  real    :: delta_Kd       ! increment to Kd from the bottom boundary layer mixing (m2/s)
  logical :: Rayleigh_drag  ! Set to true if Rayleigh drag velocities
                            ! defined in visc, on the assumption that this
                            ! extracted energy also drives diapycnal mixing.

  logical :: domore, do_i(SZI_(G))
  logical :: do_diag_Kd_BBL

  integer :: i, k, is, ie, nz, i_rem, kb_min
  is = G%isc ; ie = G%iec ; nz = G%ke

  do_diag_Kd_BBL = associated(Kd_BBL)

  if (.not.(CS%bottomdraglaw .and. (CS%BBL_effic>0.0))) return

  cdrag_sqrt = sqrt(CS%cdrag)
  TKE_Ray = 0.0 ; Rayleigh_drag = .false.
  if (associated(visc%Ray_u) .and. associated(visc%Ray_v)) Rayleigh_drag = .true.

  I_Rho0 = 1.0/GV%Rho0
  R0_g = GV%Rho0/GV%g_Earth

  do K=2,nz ; Rint(K) = 0.5*(GV%Rlay(k-1)+GV%Rlay(k)) ; enddo

  kb_min = max(GV%nk_rho_varies+1,2)

  ! The turbulence decay scale is 0.5*ustar/f from K&E & MOM_vertvisc.F90
  ! Any turbulence that makes it into the mixed layers is assumed
  ! to be relatively small and is discarded.
  do i=is,ie
    ustar_h = visc%ustar_BBL(i,j)
    if (ASSOCIATED(fluxes%ustar_tidal)) &
      ustar_h = ustar_h + fluxes%ustar_tidal(i,j)
    absf = 0.25*((abs(G%CoriolisBu(I-1,J-1)) + abs(G%CoriolisBu(I,J))) + &
                 (abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J-1))))
    if ((ustar_h > 0.0) .and. (absf > 0.5*CS%IMax_decay*ustar_h))  then
      I2decay(i) = absf / ustar_h
    else
      ! The maximum decay scale should be something of order 200 m.
      ! If ustar_h = 0, this is land so this value doesn't matter.
      I2decay(i) = 0.5*CS%IMax_decay
    endif
    TKE(i) = ((CS%BBL_effic * cdrag_sqrt) * &
              exp(-I2decay(i)*(GV%H_to_m*h(i,j,nz))) ) * &
             visc%TKE_BBL(i,j)

    if (ASSOCIATED(fluxes%TKE_tidal)) &
      TKE(i) = TKE(i) + fluxes%TKE_tidal(i,j) * I_Rho0 * &
           (CS%BBL_effic * exp(-I2decay(i)*(GV%H_to_m*h(i,j,nz))))

    ! Distribute the work over a BBL of depth 20^2 ustar^2 / g' following
    ! Killworth & Edwards (1999) and Zilitikevich & Mironov (1996).
    ! Rho_top is determined by finding the density where
    ! integral(bottom, Z) (rho(z') - rho(Z)) dz' = rho_0 400 ustar^2 / g

    gh_sum_top(i) = R0_g * 400.0 * ustar_h**2

    do_i(i) = (G%mask2dT(i,j) > 0.5)
    htot(i) = GV%H_to_m*h(i,j,nz)
    rho_htot(i) = GV%Rlay(nz)*(GV%H_to_m*h(i,j,nz))
    Rho_top(i) = GV%Rlay(1)
    if (CS%bulkmixedlayer .and. do_i(i)) Rho_top(i) = GV%Rlay(kb(i)-1)
  enddo

  do k=nz-1,2,-1 ; domore = .false.
    do i=is,ie ; if (do_i(i)) then
      htot(i) = htot(i) + GV%H_to_m*h(i,j,k)
      rho_htot(i) = rho_htot(i) + GV%Rlay(k)*(GV%H_to_m*h(i,j,k))
      if (htot(i)*GV%Rlay(k-1) <= (rho_htot(i) - gh_sum_top(i))) then
        ! The top of the mixing is in the interface atop the current layer.
        Rho_top(i) = (rho_htot(i) - gh_sum_top(i)) / htot(i)
        do_i(i) = .false.
      elseif (k <= kb(i)) then ; do_i(i) = .false.
      else ; domore = .true. ; endif
    endif ; enddo
    if (.not.domore) exit
  enddo ! k-loop

  do i=is,ie ; do_i(i) = (G%mask2dT(i,j) > 0.5) ; enddo
  do k=nz-1,kb_min,-1
    i_rem = 0
    do i=is,ie ; if (do_i(i)) then
      if (k<kb(i)) then ; do_i(i) = .false. ; cycle ; endif
      i_rem = i_rem + 1  ! Count the i-rows that are still being worked on.
      !   Apply vertical decay of the turbulent energy.  This energy is
      ! simply lost.
      TKE(i) = TKE(i) * exp(-I2decay(i) * (GV%H_to_m*(h(i,j,k) + h(i,j,k+1))))

!      if (maxEnt(i,k) <= 0.0) cycle
      if (maxTKE(i,k) <= 0.0) cycle

  ! This is an analytic integral where diffusity is a quadratic function of
  ! rho that goes asymptotically to 0 at Rho_top (vaguely following KPP?).
      if (TKE(i) > 0.0) then
        if (Rint(K) <= Rho_top(i)) then
          TKE_to_layer = TKE(i)
        else
          dRl = Rint(K+1) - Rint(K) ; dRbot = Rint(K+1) - Rho_top(i)
          TKE_to_layer = TKE(i) * dRl * (3.0*dRbot*(Rint(K) - Rho_top(i)) +&
                                         dRl**2) / dRbot**3
        endif
      else ; TKE_to_layer = 0.0 ; endif

      ! TKE_Ray has been initialized to 0 above.
      if (Rayleigh_drag) TKE_Ray = 0.5*CS%BBL_effic * G%IareaT(i,j) * &
            ((G%areaCu(I-1,j) * visc%Ray_u(I-1,j,k) * u(I-1,j,k)**2 + &
              G%areaCu(I,j)   * visc%Ray_u(I,j,k)   * u(I,j,k)**2) + &
             (G%areaCv(i,J-1) * visc%Ray_v(i,J-1,k) * v(i,J-1,k)**2 + &
              G%areaCv(i,J)   * visc%Ray_v(i,J,k)   * v(i,J,k)**2))

      if (TKE_to_layer + TKE_Ray > 0.0) then
        if (CS%BBL_mixing_as_max) then
          if (TKE_to_layer + TKE_Ray > maxTKE(i,k)) &
              TKE_to_layer = maxTKE(i,k) - TKE_Ray

          TKE(i) = TKE(i) - TKE_to_layer

          if (Kd(i,j,k) < (TKE_to_layer+TKE_Ray)*TKE_to_Kd(i,k)) then
            delta_Kd = (TKE_to_layer+TKE_Ray)*TKE_to_Kd(i,k) - Kd(i,j,k)
            if ((CS%Kd_max >= 0.0) .and. (delta_Kd > CS%Kd_max)) then
              delta_Kd = CS%Kd_Max
              Kd(i,j,k) = Kd(i,j,k) + delta_Kd
            else
              Kd(i,j,k) = (TKE_to_layer+TKE_Ray)*TKE_to_Kd(i,k)
            endif
            Kd_int(i,j,K) = Kd_int(i,j,K) + 0.5*delta_Kd
            Kd_int(i,j,K+1) = Kd_int(i,j,K+1) + 0.5*delta_Kd
            if (do_diag_Kd_BBL) then
              Kd_BBL(i,j,K) = Kd_BBL(i,j,K) + 0.5*delta_Kd
              Kd_BBL(i,j,K+1) = Kd_BBL(i,j,K+1) + 0.5*delta_Kd
            endif
          endif
        else
          if (Kd(i,j,k) >= maxTKE(i,k)*TKE_to_Kd(i,k)) then
            TKE_here = 0.0
            TKE(i) = TKE(i) + TKE_Ray
          elseif (Kd(i,j,k) + (TKE_to_layer+TKE_Ray)*TKE_to_Kd(i,k) > &
                  maxTKE(i,k)*TKE_to_Kd(i,k)) then
            TKE_here = ((TKE_to_layer+TKE_Ray) + Kd(i,j,k)/TKE_to_Kd(i,k)) - &
                       maxTKE(i,k)
            TKE(i) = TKE(i) - TKE_here + TKE_Ray
          else
            TKE_here = TKE_to_layer + TKE_ray
            TKE(i) = TKE(i) - TKE_to_Layer
          endif
          if (TKE(i) < 0.0) TKE(i) = 0.0 ! This should be unnecessary?

          if (TKE_here > 0.0) then
            delta_Kd = TKE_here*TKE_to_Kd(i,k)
            if (CS%Kd_max >= 0.0) delta_Kd = min(delta_Kd, CS%Kd_max)
            Kd(i,j,k) = Kd(i,j,k) + delta_Kd
            Kd_int(i,j,K) = Kd_int(i,j,K) + 0.5*delta_Kd
            Kd_int(i,j,K+1) = Kd_int(i,j,K+1) + 0.5*delta_Kd
            if (do_diag_Kd_BBL) then
              Kd_BBL(i,j,K) = Kd_BBL(i,j,K) + 0.5*delta_Kd
              Kd_BBL(i,j,K+1) = Kd_BBL(i,j,K+1) + 0.5*delta_Kd
            endif
          endif
        endif
      endif

      ! This may be risky - in the case that there are exactly zero
      ! velocities at 4 neighboring points, but nonzero velocities
      ! above the iterations would stop too soon. I don't see how this
      ! could happen in practice. RWH
      if ((TKE(i)<= 0.0) .and. (TKE_Ray == 0.0)) then
        do_i(i) = .false. ; i_rem = i_rem - 1
      endif

    endif ; enddo
    if (i_rem == 0) exit
  enddo ! k-loop

end subroutine add_drag_diffusivity

!> Calculates a BBL diffusivity use a Prandtl number 1 diffusivitiy with a law of the
!! wall turbulent viscosity, up to a BBL height where the energy used for mixing has
!! consumed the mechanical TKE input.
subroutine add_LOTW_BBL_diffusivity(h, u, v, tv, fluxes, visc, j, N2_int, &
                                    G, GV, CS, Kd, Kd_int, Kd_BBL)
  type(ocean_grid_type),                        intent(in)    :: G !< Grid structure
  type(verticalGrid_type),                      intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)),    intent(in)    :: u !< u component of flow (m s-1)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)),    intent(in)    :: v !< v component of flow (m s-1)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),     intent(in)    :: h !< Layer thickness (m or kg m-2)
  type(thermo_var_ptrs),                        intent(in)    :: tv !< Thermodynamic variables structure
  type(forcing),                                intent(in)    :: fluxes !< Surface fluxes structure
  type(vertvisc_type),                          intent(in)    :: visc !< Vertical viscosity structure
  integer,                                      intent(in)    :: j !< j-index of row to work on
  real, dimension(SZI_(G),SZK_(G)+1),           intent(in)    :: N2_int !< Square of Brunt-Vaisala at interfaces (s-2)
  type(set_diffusivity_CS),                     pointer       :: CS !< Diffusivity control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),     intent(inout) :: Kd !< Layer net diffusivity (m2 s-1)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1),   intent(inout) :: Kd_int !< Interface net diffusivity (m2 s-1)
  real, dimension(:,:,:),                       pointer       :: Kd_BBL !< Interface BBL diffusivity (m2 s-1)

  ! Local variables
  real :: TKE_column       ! net TKE input into the column (m3 s-3)
  real :: TKE_to_layer     ! TKE used to drive mixing in a layer (m3 s-3)
  real :: TKE_Ray          ! TKE from a layer Rayleigh drag used to drive mixing in that layer (m3 s-3)
  real :: TKE_remaining    ! remaining TKE available for mixing in this layer and above (m3 s-3)
  real :: TKE_consumed     ! TKE used for mixing in this layer (m3 s-3)
  real :: TKE_Kd_wall      ! TKE associated with unlimited law of the wall mixing (m3 s-3)
  real :: cdrag_sqrt       ! square root of the drag coefficient (nondimensional)
  real :: ustar            ! value of ustar at a thickness point (m/s)
  real :: ustar2           ! square of ustar, for convenience (m2/s2)
  real :: absf             ! average absolute value of Coriolis parameter around a thickness point (1/sec)
  real :: dh, dhm1         ! thickness of layers k and k-1, respecitvely (meter)
  real :: z                ! distance to interface k from bottom (meter)
  real :: D_minus_z        ! distance to interface k from surface (meter)
  real :: total_thickness  ! total thickness of water column (meter)
  real :: Idecay           ! inverse of decay scale used for "Joule heating" loss of TKE with height (1/m)
  real :: Kd_wall          ! Law of the wall diffusivity (m2/s)
  real :: Kd_lower         ! diffusivity for lower interface (m2/sec)
  real :: ustar_D          ! u* x D  (m2/s)
  real :: I_Rho0           ! 1 / rho0
  real :: N2_min           ! Minimum value of N2 to use in calculation of TKE_Kd_wall (1/s2)
  logical :: Rayleigh_drag ! Set to true if there are Rayleigh drag velocities defined in visc, on
                           ! the assumption that this extracted energy also drives diapycnal mixing.
  integer :: i, k, km1
  real, parameter :: von_karm = 0.41 ! Von Karman constant (http://en.wikipedia.org/wiki/Von_Karman_constant)
  logical :: do_diag_Kd_BBL

  if (.not.(CS%bottomdraglaw .and. (CS%BBL_effic>0.0))) return
  do_diag_Kd_BBL = associated(Kd_BBL)

  N2_min = 0.
  if (CS%LOTW_BBL_use_omega) N2_min = (CS%omega**2)

  ! Determine whether to add Rayleigh drag contribution to TKE
  Rayleigh_drag = .false.
  if (associated(visc%Ray_u) .and. associated(visc%Ray_v)) Rayleigh_drag = .true.
  I_Rho0 = 1.0/GV%Rho0
  cdrag_sqrt = sqrt(CS%cdrag)

  TKE_Ray = 0. ! In case Rayleigh_drag is not used.
  do i=G%isc,G%iec ! Developed in single-column mode

    ! Column-wise parameters.
    absf = 0.25*((abs(G%CoriolisBu(I-1,J-1)) + abs(G%CoriolisBu(I,J))) + &
                 (abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J-1)))) ! Non-zero on equator!

    ! u* at the bottom, in m s-1.
    ustar = visc%ustar_BBL(i,j)
    ustar2 = ustar**2
    ! In add_drag_diffusivity(), fluxes%ustar_tidal is added in. This might be double counting
    ! since ustar_BBL should already include all contributions to u*? -AJA
    if (ASSOCIATED(fluxes%ustar_tidal)) ustar = ustar + fluxes%ustar_tidal(i,j)

    ! The maximum decay scale should be something of order 200 m. We use the smaller of u*/f and
    ! (IMax_decay)^-1 as the decay scale. If ustar = 0, this is land so this value doesn't matter.
    Idecay = CS%IMax_decay
    if ((ustar > 0.0) .and. (absf > CS%IMax_decay*ustar)) Idecay = absf / ustar

    ! Energy input at the bottom, in m3 s-3.
    ! (Note that visc%TKE_BBL is in m3 s-3, set in set_BBL_TKE().)
    ! I am still unsure about sqrt(cdrag) in this expressions - AJA
    TKE_column = cdrag_sqrt * visc%TKE_BBL(i,j)
    ! Add in tidal dissipation energy at the bottom, in m3 s-3.
    ! Note that TKE_tidal is in W m-2.
    if (ASSOCIATED(fluxes%TKE_tidal)) TKE_column = TKE_column + fluxes%TKE_tidal(i,j) * I_Rho0
    TKE_column = CS%BBL_effic * TKE_column ! Only use a fraction of the mechanical dissipation for mixing.

    TKE_remaining = TKE_column
    total_thickness = ( sum(h(i,j,:)) + GV%H_subroundoff )* GV%H_to_m ! Total column thickness, in m.
    ustar_D = ustar * total_thickness
    z = 0.
    Kd_lower = 0. ! Diffusivity on bottom boundary.

    ! Work upwards from the bottom, accumulating work used until it exceeds the available TKE input
    ! at the bottom.
    do k=G%ke,2,-1
      dh = GV%H_to_m * h(i,j,k) ! Thickness of this level in m.
      km1 = max(k-1, 1)
      dhm1 = GV%H_to_m * h(i,j,km1) ! Thickness of level above in m.

      ! Add in additional energy input from bottom-drag against slopes (sides)
      if (Rayleigh_drag) TKE_remaining = TKE_remaining + 0.5*CS%BBL_effic * G%IareaT(i,j) * &
            ((G%areaCu(I-1,j) * visc%Ray_u(I-1,j,k) * u(I-1,j,k)**2 + &
              G%areaCu(I,j)   * visc%Ray_u(I,j,k)   * u(I,j,k)**2) + &
             (G%areaCv(i,J-1) * visc%Ray_v(i,J-1,k) * v(i,J-1,k)**2 + &
              G%areaCv(i,J)   * visc%Ray_v(i,J,k)   * v(i,J,k)**2))

      ! Exponentially decay TKE across the thickness of the layer.
      ! This is energy loss in addition to work done as mixing, apparently to Joule heating.
      TKE_remaining = exp(-Idecay*dh) * TKE_remaining

      z = z + h(i,j,k)*GV%H_to_m ! Distance between upper interface of layer and the bottom, in m.
      D_minus_z = max(total_thickness - z, 0.) ! Thickness above layer, m.

      ! Diffusivity using law of the wall, limited by rotation, at height z, in m2/s.
      ! This calculation is at the upper interface of the layer
      if ( ustar_D + absf * ( z * D_minus_z ) == 0.) then
        Kd_wall = 0.
      else
        Kd_wall = ( ( von_karm * ustar2 ) * ( z * D_minus_z ) )/( ustar_D + absf * ( z * D_minus_z ) )
      endif

      ! TKE associated with Kd_wall, in m3 s-2.
      ! This calculation if for the volume spanning the interface.
      TKE_Kd_wall = Kd_wall * 0.5 * (dh + dhm1) * max(N2_int(i,k), N2_min)

      ! Now bound Kd such that the associated TKE is no greater than available TKE for mixing.
      if (TKE_Kd_wall>0.) then
        TKE_consumed = min(TKE_Kd_wall, TKE_remaining)
        Kd_wall = (TKE_consumed/TKE_Kd_wall) * Kd_wall ! Scale Kd so that only TKE_consumed is used.
      else
        ! Either N2=0 or dh = 0.
        if (TKE_remaining>0.) then
          Kd_wall = CS%Kd_max
        else
          Kd_wall = 0.
        endif
        TKE_consumed = 0.
      endif

      ! Now use up the appropriate about of TKE associated with the diffusivity chosen
      TKE_remaining = TKE_remaining - TKE_consumed ! Note this will be non-negative

      ! Add this BBL diffusivity to the model net diffusivity.
      Kd_int(i,j,k) = Kd_int(i,j,k) + Kd_wall
      Kd(i,j,k) = Kd(i,j,k) + 0.5*(Kd_wall + Kd_lower)
      Kd_lower = Kd_wall ! Store for next level up.
      if (do_diag_Kd_BBL) Kd_BBL(i,j,k) = Kd_wall
    enddo ! k
  enddo ! i

end subroutine add_LOTW_BBL_diffusivity

subroutine add_MLrad_diffusivity(h, fluxes, j, G, GV, CS, Kd, TKE_to_Kd, Kd_int)
  type(ocean_grid_type),                    intent(in)    :: G
  type(verticalGrid_type),                  intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h
  type(forcing),                            intent(in)    :: fluxes
  integer,                                  intent(in)    :: j
  type(set_diffusivity_CS),                 pointer       :: CS
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: Kd
  real, dimension(SZI_(G),SZK_(G)),         intent(in)    :: TKE_to_Kd
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), optional, intent(inout) :: Kd_int

! This routine adds effects of mixed layer radiation to the layer diffusivities.

  real, dimension(SZI_(G)) ::         &
                         h_ml,        &
                         TKE_ml_flux, &
                         I_decay,     &
                         Kd_mlr_ml

  real :: f_sq, h_ml_sq, ustar_sq, Kd_mlr, C1_6
  real :: Omega2            ! rotation rate squared (1/s2)
  real :: z1                ! layer thickness times I_decay (nondim)
  real :: dzL               ! thickness converted to meter
  real :: I_decay_len2_TKE  ! squared inverse decay lengthscale for
                            ! TKE, as used in the mixed layer code (1/m2)
  real :: h_neglect         ! negligibly small thickness (meter)

  logical :: do_any, do_i(SZI_(G))
  integer :: i, k, is, ie, nz, kml
  is = G%isc ; ie = G%iec ; nz = G%ke

  Omega2    = CS%Omega**2
  C1_6      = 1.0 / 6.0
  kml       = GV%nkml
  h_neglect = GV%H_subroundoff*GV%H_to_m

  if (.not.CS%ML_radiation) return

  do i=is,ie ; h_ml(i) = 0.0 ; do_i(i) = (G%mask2dT(i,j) > 0.5) ; enddo
  do k=1,kml ; do i=is,ie ; h_ml(i) = h_ml(i) + GV%H_to_m*h(i,j,k) ; enddo ; enddo

  do i=is,ie ; if (do_i(i)) then
    if (CS%ML_omega_frac >= 1.0) then
      f_sq = 4.0*Omega2
    else
      f_sq = 0.25*((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                   (G%CoriolisBu(I,J-1)**2 + G%CoriolisBu(I-1,J)**2))
      if (CS%ML_omega_frac > 0.0) &
        f_sq = CS%ML_omega_frac*4.0*Omega2 + (1.0-CS%ML_omega_frac)*f_sq
    endif

    ustar_sq = max(fluxes%ustar(i,j), CS%ustar_min)**2

    TKE_ml_flux(i) = (CS%mstar*CS%ML_rad_coeff)*(ustar_sq*fluxes%ustar(i,j))
    I_decay_len2_TKE = CS%TKE_decay**2 * (f_sq / ustar_sq)

    if (CS%ML_rad_TKE_decay) &
      TKE_ml_flux(i) = TKE_ml_flux(i) * exp(-h_ml(i) * sqrt(I_decay_len2_TKE))

    ! Calculate the inverse decay scale
    h_ml_sq = (CS%ML_rad_efold_coeff * (h_ml(i)+h_neglect))**2
    I_decay(i) = sqrt((I_decay_len2_TKE * h_ml_sq + 1.0) / h_ml_sq)

    ! Average the dissipation layer kml+1, using
    ! a more accurate Taylor series approximations for very thin layers.
    z1 = (GV%H_to_m*h(i,j,kml+1)) * I_decay(i)
    if (z1 > 1e-5) then
      Kd_mlr = (TKE_ml_flux(i) * TKE_to_Kd(i,kml+1)) * &
               (1.0 - exp(-z1))
    else
      Kd_mlr = (TKE_ml_flux(i) * TKE_to_Kd(i,kml+1)) * &
               (z1 * (1.0 - z1 * (0.5 - C1_6*z1)))
    endif
    Kd_mlr_ml(i) = min(Kd_mlr,CS%ML_rad_kd_max)
    TKE_ml_flux(i) = TKE_ml_flux(i) * exp(-z1)
  endif ; enddo

  do k=1,kml+1 ; do i=is,ie ; if (do_i(i)) then
    Kd(i,j,k) = Kd(i,j,k) + Kd_mlr_ml(i)
  endif ; enddo ; enddo
  if (present(Kd_int)) then
    do K=2,kml+1 ; do i=is,ie ; if (do_i(i)) then
      Kd_int(i,j,K) = Kd_int(i,j,K) + Kd_mlr_ml(i)
    endif ; enddo ; enddo
    if (kml<=nz-1) then ; do i=is,ie ; if (do_i(i)) then
      Kd_int(i,j,Kml+2) = Kd_int(i,j,Kml+2) + 0.5*Kd_mlr_ml(i)
    endif ; enddo ; endif
  endif

  do k=kml+2,nz-1
    do_any = .false.
    do i=is,ie ; if (do_i(i)) then
      dzL = GV%H_to_m*h(i,j,k) ;  z1 = dzL*I_decay(i)
      if (z1 > 1e-5) then
        Kd_mlr = (TKE_ml_flux(i) * TKE_to_Kd(i,k)) * &
                 ((1.0 - exp(-z1)) / dzL)
      else
        Kd_mlr = (TKE_ml_flux(i) * TKE_to_Kd(i,k)) * &
                 (I_decay(i) * (1.0 - z1 * (0.5 - C1_6*z1)))
      endif
      Kd_mlr = min(Kd_mlr,CS%ML_rad_kd_max)
      Kd(i,j,k) =  Kd(i,j,k) + Kd_mlr
      if (present(Kd_int)) then
        Kd_int(i,j,K) = Kd_int(i,j,K) + 0.5*Kd_mlr
        Kd_int(i,j,K+1) = Kd_int(i,j,K+1) + 0.5*Kd_mlr
      endif

      TKE_ml_flux(i) = TKE_ml_flux(i) * exp(-z1)
      if (TKE_ml_flux(i) * I_decay(i) < 0.1*CS%Kd_min*Omega2) then
        do_i(i) = .false.
      else ; do_any = .true. ; endif
    endif ; enddo
    if (.not.do_any) exit
  enddo

end subroutine add_MLrad_diffusivity

subroutine add_int_tide_diffusivity(h, N2_bot, j, TKE_to_Kd, max_TKE, G, GV, CS, &
                                    dd, N2_lay, Kd, Kd_int )
  type(ocean_grid_type),                    intent(in)    :: G
  type(verticalGrid_type),                  intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h
  real, dimension(SZI_(G)),                 intent(in)    :: N2_bot
  real, dimension(SZI_(G),SZK_(G)),         intent(in)    :: N2_lay
  integer,                                  intent(in)    :: j
  real, dimension(SZI_(G),SZK_(G)),         intent(in)    :: TKE_to_Kd, max_TKE
  type(set_diffusivity_CS),                 pointer       :: CS
  type(diffusivity_diags),                  intent(inout) :: dd
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: Kd
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), optional, intent(inout) :: Kd_int

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
      call get_lowmode_loss(i,j,G,CS%int_tide_CSp,"WaveDrag",TKE_lowmode_tot)
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

      if (CS%Kd_max >= 0.0) Kd_add = min(Kd_add, CS%Kd_max)
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
        if (CS%Kd_max >= 0.0) Kd_add = min(Kd_add, CS%Kd_max)
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
        if (CS%Kd_max >= 0.0) Kd_add = min(Kd_add, CS%Kd_max)
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
        if (CS%Kd_max >= 0.0) Kd_add = min(Kd_add, CS%Kd_max)
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

      if (CS%Kd_max >= 0.0) Kd_add = min(Kd_add, CS%Kd_max)
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
        if (CS%Kd_max >= 0.0) Kd_add = min(Kd_add, CS%Kd_max)
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
        if (CS%Kd_max >= 0.0) Kd_add = min(Kd_add, CS%Kd_max)
        if (k>1) dd%Kd_Niku(i,j,K)    = dd%Kd_Niku(i,j,K)   + 0.5*Kd_add
        if (k<nz) dd%Kd_Niku(i,j,K+1) = dd%Kd_Niku(i,j,K+1) + 0.5*Kd_add
      endif
   !  if (associated(dd%Kd_Niku)) dd%Kd_Niku(i,j,K) = TKE_to_Kd(i,k) * TKE_Niku_lay
      if (associated(dd%Kd_Niku_work)) dd%Kd_Niku_work(i,j,k) = GV%Rho0 * TKE_Niku_lay

      if (associated(dd%Kd_lowmode)) then
        ! If at layers, dd%Kd_lowmode is just TKE_to_Kd(i,k) * TKE_lowmode_lay
        ! The following sets the interface diagnostics.
        Kd_add = TKE_to_Kd(i,k) * TKE_lowmode_lay
        if (CS%Kd_max >= 0.0) Kd_add = min(Kd_add, CS%Kd_max)
        if (k>1)  dd%Kd_lowmode(i,j,K)   = dd%Kd_lowmode(i,j,K)   + 0.5*Kd_add
        if (k<nz) dd%Kd_lowmode(i,j,K+1) = dd%Kd_lowmode(i,j,K+1) + 0.5*Kd_add
      endif
      if (associated(dd%Kd_lowmode_work)) &
        dd%Kd_lowmode_work(i,j,k) = GV%Rho0 * TKE_lowmode_lay
      if (associated(dd%Fl_lowmode)) dd%Fl_lowmode(i,j,k) = TKE_lowmode_rem(i)

    enddo ; enddo;
  endif ! Polzin

end subroutine add_int_tide_diffusivity

subroutine set_BBL_TKE(u, v, h, fluxes, visc, G, GV, CS)
  type(ocean_grid_type),                     intent(in)    :: G
  type(verticalGrid_type),                   intent(in)    :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h
  type(forcing),                             intent(in)    :: fluxes
  type(vertvisc_type),                       intent(inout) :: visc
  type(set_diffusivity_CS),                  pointer       :: CS

  ! This subroutine calculates several properties related to bottom
  ! boundary layer turbulence.

  real, dimension(SZI_(G)) :: &
    htot          ! total thickness above or below a layer, or the
                  ! integrated thickness in the BBL (meter)

  real, dimension(SZIB_(G)) :: &
    uhtot, &      ! running integral of u in the BBL (m2/s)
    ustar, &      ! bottom boundary layer turbulence speed (m/s)
    u2_bbl        ! square of the mean zonal velocity in the BBL (m2/s2)

  real :: vhtot(SZI_(G)) ! running integral of v in the BBL (m2/sec)

  real, dimension(SZI_(G),SZJB_(G)) :: &
    vstar, & ! ustar at at v-points in 2 j-rows (m/s)
    v2_bbl   ! square of average meridional velocity in BBL (m2/s2)

  real :: cdrag_sqrt  ! square root of the drag coefficient (nondim)
  real :: hvel        ! thickness at velocity points (meter)

  logical :: domore, do_i(SZI_(G))
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(CS)) call MOM_error(FATAL,"set_BBL_TKE: "//&
         "Module must be initialized before it is used.")

  if (.not.CS%bottomdraglaw .or. (CS%BBL_effic<=0.0)) then
    if (associated(visc%ustar_BBL)) then
      do j=js,je ; do i=is,ie ; visc%ustar_BBL(i,j) = 0.0 ; enddo ; enddo
    endif
    if (associated(visc%TKE_BBL)) then
      do j=js,je ; do i=is,ie ; visc%TKE_BBL(i,j) = 0.0 ; enddo ; enddo
    endif
    return
  endif

  cdrag_sqrt = sqrt(CS%cdrag)

!$OMP parallel default(none) shared(cdrag_sqrt,is,ie,js,je,nz,visc,CS,G,GV,vstar,h,v, &
!$OMP                               v2_bbl,u) &
!$OMP                       private(do_i,vhtot,htot,domore,hvel,uhtot,ustar,u2_bbl)
!$OMP do
  do J=js-1,je
    ! Determine ustar and the square magnitude of the velocity in the
    ! bottom boundary layer. Together these give the TKE source and
    ! vertical decay scale.
    do i=is,ie ; if ((G%mask2dCv(i,J) > 0.5) .and. (cdrag_sqrt*visc%bbl_thick_v(i,J) > 0.0)) then
      do_i(i) = .true. ; vhtot(i) = 0.0 ; htot(i) = 0.0
      vstar(i,J) = visc%kv_bbl_v(i,J)/(cdrag_sqrt*visc%bbl_thick_v(i,J))
    else
      do_i(i) = .false. ; vstar(i,J) = 0.0 ; htot(i) = 0.0
    endif ; enddo
    do k=nz,1,-1
      domore = .false.
      do i=is,ie ; if (do_i(i)) then
        hvel = 0.5*GV%H_to_m*(h(i,j,k) + h(i,j+1,k))
        if ((htot(i) + hvel) >= visc%bbl_thick_v(i,J)) then
          vhtot(i) = vhtot(i) + (visc%bbl_thick_v(i,J) - htot(i))*v(i,J,k)
          htot(i) = visc%bbl_thick_v(i,J)
          do_i(i) = .false.
        else
          vhtot(i) = vhtot(i) + hvel*v(i,J,k)
          htot(i) = htot(i) + hvel
          domore = .true.
        endif
      endif ; enddo
      if (.not.domore) exit
    enddo
    do i=is,ie ; if ((G%mask2dCv(i,J) > 0.5) .and. (htot(i) > 0.0)) then
      v2_bbl(i,J) = (vhtot(i)*vhtot(i))/(htot(i)*htot(i))
    else
      v2_bbl(i,J) = 0.0
    endif ; enddo
  enddo
!$OMP do
  do j=js,je
    do I=is-1,ie ; if ((G%mask2dCu(I,j) > 0.5) .and. (cdrag_sqrt*visc%bbl_thick_u(I,j) > 0.0))  then
      do_i(I) = .true. ; uhtot(I) = 0.0 ; htot(I) = 0.0
      ustar(I) = visc%kv_bbl_u(I,j)/(cdrag_sqrt*visc%bbl_thick_u(I,j))
    else
      do_i(I) = .false. ; ustar(I) = 0.0 ; htot(I) = 0.0
    endif ; enddo
    do k=nz,1,-1 ; domore = .false.
      do I=is-1,ie ; if (do_i(I)) then
        hvel = 0.5*GV%H_to_m*(h(i,j,k) + h(i+1,j,k))
        if ((htot(I) + hvel) >= visc%bbl_thick_u(I,j)) then
          uhtot(I) = uhtot(I) + (visc%bbl_thick_u(I,j) - htot(I))*u(I,j,k)
          htot(I) = visc%bbl_thick_u(I,j)
          do_i(I) = .false.
        else
          uhtot(I) = uhtot(I) + hvel*u(I,j,k)
          htot(I) = htot(I) + hvel
          domore = .true.
        endif
      endif ; enddo
      if (.not.domore) exit
    enddo
    do I=is-1,ie ; if ((G%mask2dCu(I,j) > 0.5) .and. (htot(i) > 0.0)) then
      u2_bbl(I) = (uhtot(I)*uhtot(I))/(htot(I)*htot(I))
    else
      u2_bbl(I) = 0.0
    endif ; enddo

    do i=is,ie
      visc%ustar_BBL(i,j) = sqrt(0.5*G%IareaT(i,j) * &
                ((G%areaCu(I-1,j)*(ustar(I-1)*ustar(I-1)) + &
                  G%areaCu(I,j)*(ustar(I)*ustar(I))) + &
                 (G%areaCv(i,J-1)*(vstar(i,J-1)*vstar(i,J-1)) + &
                  G%areaCv(i,J)*(vstar(i,J)*vstar(i,J))) ) )
      visc%TKE_BBL(i,j) = (((G%areaCu(I-1,j)*(ustar(I-1)*u2_bbl(I-1)) + &
                    G%areaCu(I,j) * (ustar(I)*u2_bbl(I))) + &
                   (G%areaCv(i,J-1)*(vstar(i,J-1)*v2_bbl(i,J-1)) + &
                    G%areaCv(i,J) * (vstar(i,J)*v2_bbl(i,J))))*G%IareaT(i,j))
    enddo
  enddo
!$OMP end parallel

end subroutine set_BBL_TKE

subroutine set_density_ratios(h, tv, kb, G, GV, CS, j, ds_dsp1, rho_0)
  type(ocean_grid_type),                 intent(in)    :: G
  type(verticalGrid_type),               intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
  type(thermo_var_ptrs),                 intent(in)    :: tv
  integer, dimension(SZI_(G)),            intent(in)   :: kb
  type(set_diffusivity_CS),              pointer       :: CS
  integer,                               intent(in)    :: j
  real, dimension(SZI_(G),SZK_(G)),      intent(out)   :: ds_dsp1
  real, dimension(SZI_(G),SZK_(G)), optional, intent(in) :: rho_0

! Arguments:
!  (in)      h       - layer thickness (meter)
!  (in)      tv      - structure containing pointers to any available
!                      thermodynamic fields; absent fields have NULL ptrs
!  (in)      kb      - index of lightest layer denser than the buffer layer
!  (in)      G       - ocean grid structure
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS      - control structure returned by previous call to diabatic_entrain_init
!  (in)      j       - meridional index upon which to work
!  (in)      ds_dsp1 - coordinate variable (sigma-2) difference across an
!                      interface divided by the difference across the interface
!                      below it (nondimensional)
!  (in)      rho_0   - layer potential densities relative to surface press (kg/m3)

  real :: g_R0                     ! g_R0 is g/Rho (m4 kg-1 s-2)
  real :: eps, tmp                 ! nondimensional temproray variables
  real :: a(SZK_(G)), a_0(SZK_(G)) ! nondimensional temporary variables
  real :: p_ref(SZI_(G))           ! an array of tv%P_Ref pressures
  real :: Rcv(SZI_(G),SZK_(G))     ! coordinate density in the mixed and buffer layers (kg/m3)
  real :: I_Drho                   ! temporary variable (m3/kg)

  integer :: i, k, k3, is, ie, nz, kmb
  is = G%isc ; ie = G%iec ; nz = G%ke

  do k=2,nz-1
    if (GV%g_prime(k+1)/=0.) then
      do i=is,ie
        ds_dsp1(i,k) = GV%g_prime(k) / GV%g_prime(k+1)
      enddo
    else
      do i=is,ie
        ds_dsp1(i,k) = 1.
      enddo
    endif
  enddo

  if (CS%bulkmixedlayer) then
    g_R0 = GV%g_Earth/GV%Rho0
    kmb = GV%nk_rho_varies
    eps = 0.1
    do i=is,ie ; p_ref(i) = tv%P_Ref ; enddo
    do k=1,kmb
      call calculate_density(tv%T(:,j,k),tv%S(:,j,k),p_ref,Rcv(:,k),&
                             is,ie-is+1,tv%eqn_of_state)
    enddo
    do i=is,ie
      if (kb(i) <= nz-1) then
!   Set up appropriately limited ratios of the reduced gravities of the
! interfaces above and below the buffer layer and the next denser layer.
        k = kb(i)

        I_Drho = g_R0 / GV%g_prime(k+1)
        ! The indexing convention for a is appropriate for the interfaces.
        do k3=1,kmb
          a(k3+1) = (GV%Rlay(k) - Rcv(i,k3)) * I_Drho
        enddo
        if ((present(rho_0)) .and. (a(kmb+1) < 2.0*eps*ds_dsp1(i,k))) then
!   If the buffer layer nearly matches the density of the layer below in the
! coordinate variable (sigma-2), use the sigma-0-based density ratio if it is
! greater (and stable).
          if ((rho_0(i,k) > rho_0(i,kmb)) .and. &
              (rho_0(i,k+1) > rho_0(i,k))) then
            I_Drho = 1.0 / (rho_0(i,k+1)-rho_0(i,k))
            a_0(kmb+1) = min((rho_0(i,k)-rho_0(i,kmb)) * I_Drho, ds_dsp1(i,k))
            if (a_0(kmb+1) > a(kmb+1)) then
              do k3=2,kmb
                a_0(k3) = a_0(kmb+1) + (rho_0(i,kmb)-rho_0(i,k3-1)) * I_Drho
              enddo
              if (a(kmb+1) <= eps*ds_dsp1(i,k)) then
                do k3=2,kmb+1 ; a(k3) = a_0(k3) ; enddo
              else
! Alternative...  tmp = 0.5*(1.0 - cos(PI*(a(K2+1)/(eps*ds_dsp1(i,k)) - 1.0)) )
                tmp = a(kmb+1)/(eps*ds_dsp1(i,k)) - 1.0
                do k3=2,kmb+1 ; a(k3) = tmp*a(k3) + (1.0-tmp)*a_0(k3) ; enddo
              endif
            endif
          endif
        endif

        ds_dsp1(i,k) = MAX(a(kmb+1),1e-5)

        do k3=2,kmb
!           ds_dsp1(i,k3) = MAX(a(k3),1e-5)
          ! Deliberately treat convective instabilies of the upper mixed
          ! and buffer layers with respect to the deepest buffer layer as
          ! though they don't exist.  They will be eliminated by the upcoming
          ! call to the mixedlayer code anyway.
          ! The indexing convention is appropriate for the interfaces.
          ds_dsp1(i,k3) = MAX(a(k3),ds_dsp1(i,k))
        enddo
      endif ! (kb(i) <= nz-1)
    enddo ! I-loop.
  endif ! bulkmixedlayer

end subroutine set_density_ratios

subroutine set_diffusivity_init(Time, G, GV, param_file, diag, CS, diag_to_Z_CSp, int_tide_CSp)
  type(time_type),          intent(in)    :: Time
  type(ocean_grid_type),    intent(inout) :: G
  type(verticalGrid_type),  intent(in)    :: GV
  type(param_file_type),    intent(in)    :: param_file
  type(diag_ctrl), target,  intent(inout) :: diag
  type(set_diffusivity_CS), pointer       :: CS
  type(diag_to_Z_CS),       pointer       :: diag_to_Z_CSp
  type(int_tide_CS),        pointer       :: int_tide_CSp

! Arguments:
!  (in)      Time          - current model time
!  (in)      G             - ocean grid structure
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file    - structure indicating open file to parse for params
!  (in)      diag          - structure used to regulate diagnostic output
!  (in/out)  CS            - pointer set to point to the module control structure
!  (in)      diag_to_Z_CSp - pointer to the Z-diagnostics control structure
!  (in)      int_tide_CSp  - pointer to the internal tides control structure (BDM)

  real :: decay_length, utide, zbot, hamp
  type(vardesc) :: vd
  logical :: read_tideamp, ML_use_omega
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_set_diffusivity"  ! This module's name.
  character(len=20)  :: tmpstr
  character(len=200) :: filename, tideamp_file, h2_file, Niku_TKE_input_file
  real :: Niku_scale ! local variable for scaling the Nikurashin TKE flux data
  real :: omega_frac_dflt
  integer :: i, j, is, ie, js, je
  integer :: isd, ied, jsd, jed

  if (associated(CS)) then
    call MOM_error(WARNING, "diabatic_entrain_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%diag => diag
  if (associated(int_tide_CSp))  CS%int_tide_CSp  => int_tide_CSp
  if (associated(diag_to_Z_CSp)) CS%diag_to_Z_CSp => diag_to_Z_CSp

  ! These default values always need to be set.
  CS%BBL_mixing_as_max = .true.
  CS%Kdml = 0.0 ; CS%cdrag = 0.003 ; CS%BBL_effic = 0.0 ;
  CS%bulkmixedlayer = (GV%nkml > 0)


  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")

  call get_param(param_file, mod, "INPUTDIR", CS%inputdir, default=".")
  CS%inputdir = slasher(CS%inputdir)
  call get_param(param_file, mod, "FLUX_RI_MAX", CS%FluxRi_max, &
                 "The flux Richardson number where the stratification is \n"//&
                 "large enough that N2 > omega2.  The full expression for \n"//&
                 "the Flux Richardson number is usually \n"//&
                 "FLUX_RI_MAX*N2/(N2+OMEGA2).", default=0.2)
  call get_param(param_file, mod, "OMEGA", CS%omega, &
                 "The rotation rate of the earth.", units="s-1", &
                 default=7.2921e-5)

  call get_param(param_file, mod, "ML_RADIATION", CS%ML_radiation, &
                 "If true, allow a fraction of TKE available from wind \n"//&
                 "work to penetrate below the base of the mixed layer \n"//&
                 "with a vertical decay scale determined by the minimum \n"//&
                 "of: (1) The depth of the mixed layer, (2) an Ekman \n"//&
                 "length scale.", default=.false.)
  if (CS%ML_radiation) then
    ! This give a minimum decay scale that is typically much less than Angstrom.
    CS%ustar_min = 2e-4*CS%omega*(GV%Angstrom + GV%H_subroundoff)

    call get_param(param_file, mod, "ML_RAD_EFOLD_COEFF", CS%ML_rad_efold_coeff, &
                 "A coefficient that is used to scale the penetration \n"//&
                 "depth for turbulence below the base of the mixed layer. \n"//&
                 "This is only used if ML_RADIATION is true.", units="nondim", &
                 default=0.2)
    call get_param(param_file, mod, "ML_RAD_KD_MAX", CS%ML_rad_kd_max, &
                 "The maximum diapycnal diffusivity due to turbulence \n"//&
                 "radiated from the base of the mixed layer. \n"//&
                 "This is only used if ML_RADIATION is true.", units="m2 s-1", &
                 default=1.0e-3)
    call get_param(param_file, mod, "ML_RAD_COEFF", CS%ML_rad_coeff, &
                 "The coefficient which scales MSTAR*USTAR^3 to obtain \n"//&
                 "the energy available for mixing below the base of the \n"//&
                 "mixed layer. This is only used if ML_RADIATION is true.", &
                 units="nondim", default=0.2)
    call get_param(param_file, mod, "ML_RAD_APPLY_TKE_DECAY", CS%ML_rad_TKE_decay, &
                 "If true, apply the same exponential decay to ML_rad as \n"//&
                 "is applied to the other surface sources of TKE in the \n"//&
                 "mixed layer code. This is only used if ML_RADIATION is true.",&
                 default=.true.)
    call get_param(param_file, mod, "MSTAR", CS%mstar, &
                 "The ratio of the friction velocity cubed to the TKE \n"//&
                 "input to the mixed layer.", "units=nondim", default=1.2)
    call get_param(param_file, mod, "TKE_DECAY", CS%TKE_decay, &
                 "The ratio of the natural Ekman depth to the TKE decay scale.", &
                 units="nondim", default=2.5)
    call get_param(param_file, mod, "ML_USE_OMEGA", ML_use_omega, &
                 "If true, use the absolute rotation rate instead of the \n"//&
                 "vertical component of rotation when setting the decay \n"//&
                 "scale for turbulence.", default=.false., do_not_log=.true.)
    omega_frac_dflt = 0.0
    if (ML_use_omega) then
      call MOM_error(WARNING, "ML_USE_OMEGA is depricated; use ML_OMEGA_FRAC=1.0 instead.")
      omega_frac_dflt = 1.0
    endif
    call get_param(param_file, mod, "ML_OMEGA_FRAC", CS%ML_omega_frac, &
                   "When setting the decay scale for turbulence, use this \n"//&
                   "fraction of the absolute rotation rate blended with the \n"//&
                   "local value of f, as sqrt((1-of)*f^2 + of*4*omega^2).", &
                   units="nondim", default=omega_frac_dflt)
  endif

  call get_param(param_file, mod, "BOTTOMDRAGLAW", CS%bottomdraglaw, &
                 "If true, the bottom stress is calculated with a drag \n"//&
                 "law of the form c_drag*|u|*u. The velocity magnitude \n"//&
                 "may be an assumed value or it may be based on the \n"//&
                 "actual velocity in the bottommost HBBL, depending on \n"//&
                 "LINEAR_DRAG.", default=.true.)
  if  (CS%bottomdraglaw) then
    call get_param(param_file, mod, "CDRAG", CS%cdrag, &
                 "The drag coefficient relating the magnitude of the \n"//&
                 "velocity field to the bottom stress. CDRAG is only used \n"//&
                 "if BOTTOMDRAGLAW is true.", units="nondim", default=0.003)
    call get_param(param_file, mod, "BBL_EFFIC", CS%BBL_effic, &
                 "The efficiency with which the energy extracted by \n"//&
                 "bottom drag drives BBL diffusion.  This is only \n"//&
                 "used if BOTTOMDRAGLAW is true.", units="nondim", default=0.20)
    call get_param(param_file, mod, "BBL_MIXING_MAX_DECAY", decay_length, &
                 "The maximum decay scale for the BBL diffusion, or 0 \n"//&
                 "to allow the mixing to penetrate as far as \n"//&
                 "stratification and rotation permit.  The default is 0. \n"//&
                 "This is only used if BOTTOMDRAGLAW is true.", units="m", &
                 default=0.0)

    CS%IMax_decay = 1.0/200.0
    if (decay_length > 0.0) CS%IMax_decay = 1.0/decay_length
    call get_param(param_file, mod, "BBL_MIXING_AS_MAX", CS%BBL_mixing_as_max, &
                 "If true, take the maximum of the diffusivity from the \n"//&
                 "BBL mixing and the other diffusivities. Otherwise, \n"//&
                 "diffusiviy from the BBL_mixing is simply added.", &
                 default=.true.)
    call get_param(param_file, mod, "USE_LOTW_BBL_DIFFUSIVITY", CS%use_LOTW_BBL_diffusivity, &
                 "If true, uses a simple, imprecise but non-coordinate dependent, model\n"//&
                 "of BBL mixing diffusivity based on Law of the Wall. Otherwise, uses\n"//&
                 "the original BBL scheme.", default=.false.)
    if (CS%use_LOTW_BBL_diffusivity) then
      call get_param(param_file, mod, "LOTW_BBL_USE_OMEGA", CS%LOTW_BBL_use_omega, &
                 "If true, use the maximum of Omega and N for the TKE to diffusion\n"//&
                 "calculation. Otherwise, N is N.", default=.true.)
    endif
  else
    CS%use_LOTW_BBL_diffusivity = .false. ! This parameterization depends on a u* from viscous BBL
  endif
  CS%id_Kd_BBL = register_diag_field('ocean_model','Kd_BBL',diag%axesTi,Time, &
       'Bottom Boundary Layer Diffusivity', 'meter2 sec-1')
  call get_param(param_file, mod, "SIMPLE_TKE_TO_KD", CS%simple_TKE_to_Kd, &
                 "If true, uses a simple estimate of Kd/TKE that will\n"//&
                 "work for arbitrary vertical coordinates. If false,\n"//&
                 "calculates Kd/TKE and bounds based on exact energetics/n"//&
                 "for an isopycnal layer-formulation.", &
                 default=.false.)

  call get_param(param_file, mod, "BRYAN_LEWIS_DIFFUSIVITY", &
                                CS%Bryan_Lewis_diffusivity, &
                 "If true, use a Bryan & Lewis (JGR 1979) like tanh \n"//&
                 "profile of background diapycnal diffusivity with depth.", &
                 default=.false.)
  if (CS%Bryan_Lewis_diffusivity) then
    call get_param(param_file, mod, "KD_BRYAN_LEWIS_DEEP", &
                                  CS%Kd_Bryan_Lewis_deep, &
                 "The abyssal value of a Bryan-Lewis diffusivity profile. \n"//&
                 "KD_BRYAN_LEWIS_DEEP is only used if \n"//&
                 "BRYAN_LEWIS_DIFFUSIVITY is true.", units="m2 s-1", &
                 fail_if_missing=.true.)
    call get_param(param_file, mod, "KD_BRYAN_LEWIS_SURFACE", &
                                  CS%Kd_Bryan_Lewis_surface, &
                 "The surface value of a Bryan-Lewis diffusivity profile. \n"//&
                 "KD_BRYAN_LEWIS_SURFACE is only used if \n"//&
                 "BRYAN_LEWIS_DIFFUSIVITY is true.", units="m2 s-1", &
                 fail_if_missing=.true.)
    call get_param(param_file, mod, "BRYAN_LEWIS_DEPTH_CENT", &
                                  CS%Bryan_Lewis_depth_cent, &
                 "The depth about which the transition in the Bryan-Lewis \n"//&
                 "profile is centered. BRYAN_LEWIS_DEPTH_CENT is only \n"//&
                 "used if BRYAN_LEWIS_DIFFUSIVITY is true.", units="m", &
                 fail_if_missing=.true.)
    call get_param(param_file, mod, "BRYAN_LEWIS_WIDTH_TRANS", &
                                  CS%Bryan_Lewis_width_trans, &
                 "The width of the transition in the Bryan-Lewis \n"//&
                 "profile. BRYAN_LEWIS_WIDTH_TRANS is only \n"//&
                 "used if BRYAN_LEWIS_DIFFUSIVITY is true.", units="m", &
                 fail_if_missing=.true.)
  endif

  call get_param(param_file, mod, "HENYEY_IGW_BACKGROUND", &
                                CS%Henyey_IGW_background, &
                 "If true, use a latitude-dependent scaling for the near \n"//&
                 "surface background diffusivity, as described in \n"//&
                 "Harrison & Hallberg, JPO 2008.", default=.false.)
  call get_param(param_file, mod, "HENYEY_IGW_BACKGROUND_NEW", &
                                CS%Henyey_IGW_background_new, &
                 "If true, use a better latitude-dependent scaling for the\n"//&
                 "background diffusivity, as described in \n"//&
                 "Harrison & Hallberg, JPO 2008.", default=.false.)
  if (CS%Henyey_IGW_background .and. CS%Henyey_IGW_background_new) call MOM_error(FATAL, &
                 "set_diffusivity_init: HENYEY_IGW_BACKGROUND and HENYEY_IGW_BACKGROUND_NEW "// &
                 "are mutually exclusive. Set only one or none.")
  if (CS%Henyey_IGW_background) &
    call get_param(param_file, mod, "HENYEY_N0_2OMEGA", CS%N0_2Omega, &
                  "The ratio of the typical Buoyancy frequency to twice \n"//&
                  "the Earth's rotation period, used with the Henyey \n"//&
                  "scaling from the mixing.", units="nondim", default=20.0)
  call get_param(param_file, mod, "N2_FLOOR_IOMEGA2", CS%N2_FLOOR_IOMEGA2, &
                  "The floor applied to N2(k) scaled by Omega^2:\n"//&
                  "\tIf =0., N2(k) is simply positive definite.\n"//&
                  "\tIf =1., N2(k) > Omega^2 everywhere.", units="nondim", &
                  default=1.0)

  call get_param(param_file, mod, "KD_TANH_LAT_FN", &
                                  CS%Kd_tanh_lat_fn, &
                 "If true, use a tanh dependence of Kd_sfc on latitude, \n"//&
                 "like CM2.1/CM2M.  There is no physical justification \n"//&
                 "for this form, and it can not be used with \n"//&
                 "HENYEY_IGW_BACKGROUND.", default=.false.)
  if (CS%Kd_tanh_lat_fn) &
    call get_param(param_file, mod, "KD_TANH_LAT_SCALE", &
                                  CS%Kd_tanh_lat_scale, &
                 "A nondimensional scaling for the range ofdiffusivities \n"//&
                 "with KD_TANH_LAT_FN. Valid values are in the range of \n"//&
                 "-2 to 2; 0.4 reproduces CM2M.", units="nondim", default=0.0)

  call get_param(param_file, mod, "KV", CS%Kv, &
                 "The background kinematic viscosity in the interior. \n"//&
                 "The molecular value, ~1e-6 m2 s-1, may be used.", &
                 units="m2 s-1", fail_if_missing=.true.)

  call get_param(param_file, mod, "KD", CS%Kd, &
                 "The background diapycnal diffusivity of density in the \n"//&
                 "interior. Zero or the molecular value, ~1e-7 m2 s-1, \n"//&
                 "may be used.", units="m2 s-1", fail_if_missing=.true.)
  call get_param(param_file, mod, "KD_MIN", CS%Kd_min, &
                 "The minimum diapycnal diffusivity.", &
                 units="m2 s-1", default=0.01*CS%Kd)
  call get_param(param_file, mod, "KD_MAX", CS%Kd_max, &
                 "The maximum permitted increment for the diapycnal \n"//&
                 "diffusivity from TKE-based parameterizations, or a \n"//&
                 "negative value for no limit.", units="m2 s-1", default=-1.0)
  if (CS%simple_TKE_to_Kd .and. CS%Kd_max<=0.) call MOM_error(FATAL, &
         "set_diffusivity_init: To use SIMPLE_TKE_TO_KD, KD_MAX must be set to >0.")
  call get_param(param_file, mod, "KD_ADD", CS%Kd_add, &
                 "A uniform diapycnal diffusivity that is added \n"//&
                 "everywhere without any filtering or scaling.", &
                 units="m2 s-1", default=0.0)
  if (CS%use_LOTW_BBL_diffusivity .and. CS%Kd_max<=0.) call MOM_error(FATAL, &
                 "set_diffusivity_init: KD_MAX must be set (positive) when "// &
                 "USE_LOTW_BBL_DIFFUSIVITY=True.")

  if (CS%bulkmixedlayer) then
    ! Check that Kdml is not set when using bulk mixed layer
    call get_param(param_file, mod, "KDML", CS%Kdml, default=-1.)
    if (CS%Kdml>0.) call MOM_error(FATAL, &
                 "set_diffusivity_init: KDML cannot be set when using"// &
                 "bulk mixed layer.")
    CS%Kdml = CS%Kd ! This is not used with a bulk mixed layer, but also
                    ! cannot be a NaN.
  else
    call get_param(param_file, mod, "KDML", CS%Kdml, &
                 "If BULKMIXEDLAYER is false, KDML is the elevated \n"//&
                 "diapycnal diffusivity in the topmost HMIX of fluid. \n"//&
                 "KDML is only used if BULKMIXEDLAYER is false.", &
                 units="m2 s-1", default=CS%Kd)
    call get_param(param_file, mod, "HMIX_FIXED", CS%Hmix, &
                 "The prescribed depth over which the near-surface \n"//&
                 "viscosity and diffusivity are elevated when the bulk \n"//&
                 "mixed layer is not used.", units="m", fail_if_missing=.true.)
  endif
  call get_param(param_file, mod, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", default=.false.)

  call get_param(param_file, mod, "INT_TIDE_DISSIPATION", CS%Int_tide_dissipation, &
                 "If true, use an internal tidal dissipation scheme to \n"//&
                 "drive diapycnal mixing, along the lines of St. Laurent \n"//&
                 "et al. (2002) and Simmons et al. (2004).", default=.false.)
  if (CS%Int_tide_dissipation) then
    call get_param(param_file, mod, "INT_TIDE_PROFILE", tmpstr, &
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

  call get_param(param_file, mod, "LEE_WAVE_DISSIPATION", CS%Lee_wave_dissipation, &
                 "If true, use an lee wave driven dissipation scheme to \n"//&
                 "drive diapycnal mixing, along the lines of Nikurashin \n"//&
                 "(2010) and using the St. Laurent et al. (2002) \n"//&
                 "and Simmons et al. (2004) vertical profile", default=.false.)
  if (CS%lee_wave_dissipation) then
    call get_param(param_file, mod, "LEE_WAVE_PROFILE", tmpstr, &
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

  call get_param(param_file, mod, "INT_TIDE_LOWMODE_DISSIPATION", CS%Lowmode_itidal_dissipation, &
                 "If true, consider mixing due to breaking low modes that \n"//&
                 "have been remotely generated; as with itidal drag on the \n"//&
                 "barotropic tide, use an internal tidal dissipation scheme to \n"//&
                 "drive diapycnal mixing, along the lines of St. Laurent \n"//&
                 "et al. (2002) and Simmons et al. (2004).", default=.false.)

  if ((CS%Int_tide_dissipation .and. (CS%int_tide_profile == POLZIN_09)) .or. &
      (CS%lee_wave_dissipation .and. (CS%lee_wave_profile == POLZIN_09))) then
    call get_param(param_file, mod, "NU_POLZIN", CS%Nu_Polzin, &
                 "When the Polzin decay profile is used, this is a \n"//&
                 "non-dimensional constant in the expression for the \n"//&
                 "vertical scale of decay for the tidal energy dissipation.", &
                 units="nondim", default=0.0697)
    call get_param(param_file, mod, "NBOTREF_POLZIN", CS%Nbotref_Polzin, &
                 "When the Polzin decay profile is used, this is the \n"//&
                 "Rreference value of the buoyancy frequency at the ocean \n"//&
                 "bottom in the Polzin formulation for the vertical \n"//&
                 "scale of decay for the tidal energy dissipation.", &
                 units="s-1", default=9.61e-4)
    call get_param(param_file, mod, "POLZIN_DECAY_SCALE_FACTOR", &
                 CS%Polzin_decay_scale_factor, &
                 "When the Polzin decay profile is used, this is a \n"//&
                 "scale factor for the vertical scale of decay of the tidal \n"//&
                 "energy dissipation.", default=1.0, units="nondim")
    call get_param(param_file, mod, "POLZIN_SCALE_MAX_FACTOR", &
                 CS%Polzin_decay_scale_max_factor, &
                 "When the Polzin decay profile is used, this is a factor \n"//&
                 "to limit the vertical scale of decay of the tidal \n"//&
                 "energy dissipation to POLZIN_DECAY_SCALE_MAX_FACTOR \n"//&
                 "times the depth of the ocean.", units="nondim", default=1.0)
    call get_param(param_file, mod, "POLZIN_MIN_DECAY_SCALE", CS%Polzin_min_decay_scale, &
                 "When the Polzin decay profile is used, this is the \n"//&
                 "minimum vertical decay scale for the vertical profile\n"//&
                 "of internal tide dissipation with the Polzin (2009) formulation", &
                 units="m", default=0.0)
  endif
  call get_param(param_file, mod, "USER_CHANGE_DIFFUSIVITY", CS%user_change_diff, &
                 "If true, call user-defined code to change the diffusivity.", &
                 default=.false.)

  call get_param(param_file, mod, "DISSIPATION_MIN", CS%dissip_min, &
                 "The minimum dissipation by which to determine a lower \n"//&
                 "bound of Kd (a floor).", units="W m-3", default=0.0)
  call get_param(param_file, mod, "DISSIPATION_N0", CS%dissip_N0, &
                 "The intercept when N=0 of the N-dependent expression \n"//&
                 "used to set a minimum dissipation by which to determine \n"//&
                 "a lower bound of Kd (a floor): A in eps_min = A + B*N.", &
                 units="W m-3", default=0.0)
  call get_param(param_file, mod, "DISSIPATION_N1", CS%dissip_N1, &
                 "The coefficient multiplying N, following Gargett, used to \n"//&
                 "set a minimum dissipation by which to determine a lower \n"//&
                 "bound of Kd (a floor): B in eps_min = A + B*N", &
                 units="J m-3", default=0.0)
  call get_param(param_file, mod, "DISSIPATION_KD_MIN", CS%dissip_Kd_min, &
                 "The minimum vertical diffusivity applied as a floor.", &
                 units="m2 s-1", default=0.0)

  CS%limit_dissipation = (CS%dissip_min>0.) .or. (CS%dissip_N1>0.) .or. &
                         (CS%dissip_N0>0.) .or. (CS%dissip_Kd_min>0.)
  CS%dissip_N2 = 0.0
  if (CS%FluxRi_max > 0.0) &
    CS%dissip_N2 = CS%dissip_Kd_min * GV%Rho0 / CS%FluxRi_max

  if (CS%Int_tide_dissipation .or. CS%Lee_wave_dissipation) then
    call get_param(param_file, mod, "INT_TIDE_DECAY_SCALE", CS%Int_tide_decay_scale, &
                 "The decay scale away from the bottom for tidal TKE with \n"//&
                 "the new coding when INT_TIDE_DISSIPATION is used.", &
                 units="m", default=0.0)
    call get_param(param_file, mod, "MU_ITIDES", CS%Mu_itides, &
                 "A dimensionless turbulent mixing efficiency used with \n"//&
                 "INT_TIDE_DISSIPATION, often 0.2.", units="nondim", default=0.2)
    call get_param(param_file, mod, "GAMMA_ITIDES", CS%Gamma_itides, &
                 "The fraction of the internal tidal energy that is \n"//&
                 "dissipated locally with INT_TIDE_DISSIPATION.  \n"//&
                 "THIS NAME COULD BE BETTER.", &
                 units="nondim", default=0.3333)
    call get_param(param_file, mod, "MIN_ZBOT_ITIDES", CS%min_zbot_itides, &
                 "Turn off internal tidal dissipation when the total \n"//&
                 "ocean depth is less than this value.", units="m", default=0.0)

    call safe_alloc_ptr(CS%Nb,isd,ied,jsd,jed)
    call safe_alloc_ptr(CS%h2,isd,ied,jsd,jed)
    call safe_alloc_ptr(CS%TKE_itidal,isd,ied,jsd,jed)
    call safe_alloc_ptr(CS%mask_itidal,isd,ied,jsd,jed) ; CS%mask_itidal(:,:) = 1.0

    call get_param(param_file, mod, "KAPPA_ITIDES", CS%kappa_itides, &
                 "A topographic wavenumber used with INT_TIDE_DISSIPATION. \n"//&
                 "The default is 2pi/10 km, as in St.Laurent et al. 2002.", &
                 units="m-1", default=8.e-4*atan(1.0))

    call get_param(param_file, mod, "UTIDE", CS%utide, &
                 "The constant tidal amplitude used with INT_TIDE_DISSIPATION.", &
                 units="m s-1", default=0.0)
    call safe_alloc_ptr(CS%tideamp,is,ie,js,je) ; CS%tideamp(:,:) = CS%utide

    call get_param(param_file, mod, "KAPPA_H2_FACTOR", CS%kappa_h2_factor, &
                 "A scaling factor for the roughness amplitude with n"//&
                 "INT_TIDE_DISSIPATION.",  units="nondim", default=1.0)
    call get_param(param_file, mod, "TKE_ITIDE_MAX", CS%TKE_itide_max, &
                 "The maximum internal tide energy source availble to mix \n"//&
                 "above the bottom boundary layer with INT_TIDE_DISSIPATION.", &
                 units="W m-2",  default=1.0e3)

    call get_param(param_file, mod, "READ_TIDEAMP", read_tideamp, &
                 "If true, read a file (given by TIDEAMP_FILE) containing \n"//&
                 "the tidal amplitude with INT_TIDE_DISSIPATION.", default=.false.)
    if (read_tideamp) then
      call get_param(param_file, mod, "TIDEAMP_FILE", tideamp_file, &
                 "The path to the file containing the spatially varying \n"//&
                 "tidal amplitudes with INT_TIDE_DISSIPATION.", default="tideamp.nc")
      filename = trim(CS%inputdir) // trim(tideamp_file)
      call log_param(param_file, mod, "INPUTDIR/TIDEAMP_FILE", filename)
      call read_data(filename, 'tideamp', CS%tideamp, &
                     domain=G%domain%mpp_domain, timelevel=1)
    endif

    call get_param(param_file, mod, "H2_FILE", h2_file, &
                 "The path to the file containing the sub-grid-scale \n"//&
                 "topographic roughness amplitude with INT_TIDE_DISSIPATION.", &
                 fail_if_missing=.true.)
    filename = trim(CS%inputdir) // trim(h2_file)
    call log_param(param_file, mod, "INPUTDIR/H2_FILE", filename)
    call read_data(filename, 'h2', CS%h2, domain=G%domain%mpp_domain, &
                   timelevel=1)

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

  CS%id_Kd_layer = register_diag_field('ocean_model', 'Kd_layer', diag%axesTL, Time, &
      'Diapycnal diffusivity of layers (as set)', 'meter2 second-1')

  if (CS%Lee_wave_dissipation) then

    call get_param(param_file, mod, "NIKURASHIN_TKE_INPUT_FILE",Niku_TKE_input_file, &
                 "The path to the file containing the TKE input from lee \n"//&
                 "wave driven mixing. Used with LEE_WAVE_DISSIPATION.", &
                 fail_if_missing=.true.)
    call get_param(param_file, mod, "NIKURASHIN_SCALE",Niku_scale, &
                 "A non-dimensional factor by which to scale the lee-wave \n"//&
                 "driven TKE input. Used with LEE_WAVE_DISSIPATION.", &
                 units="nondim", default=1.0)

    filename = trim(CS%inputdir) // trim(Niku_TKE_input_file)
    call log_param(param_file, mod, "INPUTDIR/NIKURASHIN_TKE_INPUT_FILE", &
                   filename)
    call safe_alloc_ptr(CS%TKE_Niku,is,ie,js,je); CS%TKE_Niku(:,:) = 0.0
    call read_data(filename, 'TKE_input', CS%TKE_Niku, &
                   domain=G%domain%mpp_domain, timelevel=1 ) ! ??? timelevel -aja
    CS%TKE_Niku(:,:) = Niku_scale * CS%TKE_Niku(:,:)

    call get_param(param_file, mod, "GAMMA_NIKURASHIN",CS%Gamma_lee, &
                 "The fraction of the lee wave energy that is dissipated \n"//&
                 "locally with LEE_WAVE_DISSIPATION.", units="nondim", &
                 default=0.3333)
    call get_param(param_file, mod, "DECAY_SCALE_FACTOR_LEE",CS%Decay_scale_factor_lee, &
                 "Scaling for the vertical decay scaleof the local \n"//&
                 "dissipation of lee waves dissipation.", units="nondim", &
                 default=1.0)

    CS%id_TKE_leewave = register_diag_field('ocean_model','TKE_leewave',diag%axesT1,Time, &
        'Lee wave Driven Turbulent Kinetic Energy', 'Watt meter-2')
    CS%id_Kd_Niku = register_diag_field('ocean_model','Kd_Nikurashin',diag%axesTi,Time, &
         'Lee Wave Driven Diffusivity', 'meter2 sec-1')
  else
    CS%Decay_scale_factor_lee = -9.e99 ! This should never be used if CS%Lee_wave_dissipation = False
  endif

  if (CS%Int_tide_dissipation .or. CS%Lee_wave_dissipation .or. &
      CS%Lowmode_itidal_dissipation) then

    CS%id_TKE_itidal = register_diag_field('ocean_model','TKE_itidal',diag%axesT1,Time, &
        'Internal Tide Driven Turbulent Kinetic Energy', 'Watt meter-2')
    CS%id_maxTKE = register_diag_field('ocean_model','maxTKE',diag%axesTL,Time, &
           'Maximum layer TKE', 'meter3 second-3')
    CS%id_TKE_to_Kd = register_diag_field('ocean_model','TKE_to_Kd',diag%axesTL,Time, &
           'Convert TKE to Kd', 'second2 meter')

    CS%id_Nb = register_diag_field('ocean_model','Nb',diag%axesT1,Time, &
         'Bottom Buoyancy Frequency', 'sec-1')

    CS%id_Kd_itidal = register_diag_field('ocean_model','Kd_itides',diag%axesTi,Time, &
         'Internal Tide Driven Diffusivity', 'meter2 sec-1')

    CS%id_Kd_lowmode = register_diag_field('ocean_model','Kd_lowmode',diag%axesTi,Time, &
         'Internal Tide Driven Diffusivity (from propagating low modes)', 'meter2 sec-1')

    CS%id_Fl_itidal = register_diag_field('ocean_model','Fl_itides',diag%axesTi,Time, &
        'Vertical flux of tidal turbulent dissipation', 'meter3 sec-3')

    CS%id_Fl_lowmode = register_diag_field('ocean_model','Fl_lowmode',diag%axesTi,Time, &
         'Vertical flux of tidal turbulent dissipation (from propagating low modes)', 'meter3 sec-3')

    CS%id_Polzin_decay_scale = register_diag_field('ocean_model','Polzin_decay_scale',diag%axesT1,Time, &
         'Vertical decay scale for the tidal turbulent dissipation with Polzin scheme', 'meter')

    CS%id_Polzin_decay_scale_scaled = register_diag_field('ocean_model','Polzin_decay_scale_scaled',diag%axesT1,Time, &
         'Vertical decay scale for the tidal turbulent dissipation with Polzin scheme, scaled by N2_bot/N2_meanz', 'meter')

    CS%id_N2_bot = register_diag_field('ocean_model','N2_b',diag%axesT1,Time, &
         'Bottom Buoyancy frequency squared', 's-2')

    CS%id_N2_meanz = register_diag_field('ocean_model','N2_meanz',diag%axesT1,Time, &
         'Buoyancy frequency squared averaged over the water column', 's-2')

    CS%id_Kd_Work = register_diag_field('ocean_model','Kd_Work',diag%axesTL,Time, &
         'Work done by Diapycnal Mixing', 'Watts m-2')

    CS%id_Kd_Itidal_Work = register_diag_field('ocean_model','Kd_Itidal_Work',diag%axesTL,Time, &
         'Work done by Internal Tide Diapycnal Mixing', 'Watts m-2')

    CS%id_Kd_Niku_Work = register_diag_field('ocean_model','Kd_Nikurashin_Work',diag%axesTL,Time, &
         'Work done by Nikurashin Lee Wave Drag Scheme', 'Watts m-2')

    CS%id_Kd_Lowmode_Work = register_diag_field('ocean_model','Kd_Lowmode_Work',diag%axesTL,Time, &
         'Work done by Internal Tide Diapycnal Mixing (low modes)', 'Watts m-2')

    CS%id_N2 = register_diag_field('ocean_model','N2',diag%axesTi,Time,            &
         'Buoyancy frequency squared', 'sec-2', cmor_field_name='obvfsq',          &
          cmor_units='s-2', cmor_long_name='Square of seawater buoyancy frequency',&
          cmor_standard_name='square_of_brunt_vaisala_frequency_in_sea_water')

    if (CS%user_change_diff) &
      CS%id_Kd_user = register_diag_field('ocean_model','Kd_user',diag%axesTi,Time, &
           'User-specified Extra Diffusivity', 'meter2 sec-1')

    if (associated(diag_to_Z_CSp)) then
      vd = var_desc("N2", "second-2",&
                    "Buoyancy frequency, interpolated to z", z_grid='z')
      CS%id_N2_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
      vd = var_desc("Kd_itides","meter2 second-1", &
                    "Internal Tide Driven Diffusivity, interpolated to z", z_grid='z')
      CS%id_Kd_itidal_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
      if (CS%Lee_wave_dissipation) then
         vd = var_desc("Kd_Nikurashin", "meter2 second-1", &
                       "Lee Wave Driven Diffusivity, interpolated to z", z_grid='z')
         CS%id_Kd_Niku_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
      endif
      if (CS%Lowmode_itidal_dissipation) then
        vd = var_desc("Kd_lowmode","meter2 second-1", &
                  "Internal Tide Driven Diffusivity (from low modes), interpolated to z",&
                  z_grid='z')
        CS%id_Kd_lowmode_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
      endif
      if (CS%user_change_diff) &
        CS%id_Kd_user_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    endif
  endif

  call get_param(param_file, mod, "DOUBLE_DIFFUSION", CS%double_diffusion, &
                 "If true, increase diffusivitives for temperature or salt \n"//&
                 "based on double-diffusive paramaterization from MOM4/KPP.", &
                 default=.false.)
  if (CS%double_diffusion) then
    call get_param(param_file, mod, "MAX_RRHO_SALT_FINGERS", CS%Max_Rrho_salt_fingers, &
                 "Maximum density ratio for salt fingering regime.", &
                 default=2.55, units="nondim")
    call get_param(param_file, mod, "MAX_SALT_DIFF_SALT_FINGERS", CS%Max_salt_diff_salt_fingers, &
                 "Maximum salt diffusivity for salt fingering regime.", &
                 default=1.e-4, units="m2 s-1")
    call get_param(param_file, mod, "KV_MOLECULAR", CS%Kv_molecular, &
                 "Molecular viscosity for calculation of fluxes under \n"//&
                 "double-diffusive convection.", default=1.5e-6, units="m2 s-1")
    ! The default molecular viscosity follows the CCSM4.0 and MOM4p1 defaults.

    CS%id_KT_extra = register_diag_field('ocean_model','KT_extra',diag%axesTi,Time, &
         'Double-diffusive diffusivity for temperature', 'meter2 sec-1')

    CS%id_KS_extra = register_diag_field('ocean_model','KS_extra',diag%axesTi,Time, &
         'Double-diffusive diffusivity for salinity', 'meter2 sec-1')

    if (associated(diag_to_Z_CSp)) then
      vd = var_desc("KT_extra", "meter2 second-1", &
                    "Double-Diffusive Temperature Diffusivity, interpolated to z", &
                    z_grid='z')
      CS%id_KT_extra_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
      vd = var_desc("KS_extra", "meter2 second-1", &
                    "Double-Diffusive Salinity Diffusivity, interpolated to z",&
                    z_grid='z')
      CS%id_KS_extra_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
      vd = var_desc("Kd_BBL", "meter2 second-1", &
                    "Bottom Boundary Layer Diffusivity", z_grid='z')
      CS%id_Kd_BBL_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    endif
  endif

  if (CS%Int_tide_dissipation .and. CS%Bryan_Lewis_diffusivity) &
    call MOM_error(FATAL,"MOM_Set_Diffusivity: "// &
         "Bryan-Lewis and internal tidal dissipation are both enabled. Choose one.")
  if (CS%Henyey_IGW_background .and. CS%Kd_tanh_lat_fn) call MOM_error(FATAL, &
    "Set_diffusivity: KD_TANH_LAT_FN can not be used with HENYEY_IGW_BACKGROUND.")

  if (CS%user_change_diff) then
    call user_change_diff_init(Time, G, param_file, diag, CS%user_change_diff_CSp)
  endif

  CS%useKappaShear = kappa_shear_init(Time, G, GV, param_file, CS%diag, CS%kappaShear_CSp)
  if (CS%useKappaShear) &
    id_clock_kappaShear = cpu_clock_id('(Ocean kappa_shear)', grain=CLOCK_MODULE)
  CS%useCVMix = CVMix_shear_init(Time, G, GV, param_file, CS%diag, CS%CVMix_shear_CSp)



end subroutine set_diffusivity_init

subroutine set_diffusivity_end(CS)
  type(set_diffusivity_CS), pointer :: CS

  if (CS%user_change_diff) &
    call user_change_diff_end(CS%user_change_diff_CSp)

  if (associated(CS)) deallocate(CS)

end subroutine set_diffusivity_end

end module MOM_set_diffusivity
