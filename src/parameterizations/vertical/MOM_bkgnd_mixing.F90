!> Interface to background mixing schemes, including the Bryan and Lewis (1979)
!! which is applied via CVMix.

module MOM_bkgnd_mixing

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator,   only : diag_ctrl, time_type, register_diag_field
use MOM_diag_mediator,   only : post_data
use MOM_EOS,             only : calculate_density, calculate_density_derivs
use MOM_variables,       only : thermo_var_ptrs
use MOM_forcing_type,    only : forcing
use MOM_error_handler,   only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,     only : openParameterBlock, closeParameterBlock
use MOM_debugging,       only : hchksum
use MOM_grid,            only : ocean_grid_type
use MOM_verticalGrid,    only : verticalGrid_type
use MOM_file_parser,     only : get_param, log_version, param_file_type
use cvmix_background,    only : cvmix_init_bkgnd, cvmix_coeffs_bkgnd
use MOM_intrinsic_functions, only : invcosh

implicit none ; private

#include <MOM_memory.h>

public bkgnd_mixing_init, bkgnd_mixing_end, calculate_bkgnd_mixing

!> Control structure including parameters for this module.
type, public :: bkgnd_mixing_cs

  ! Parameters
  real    :: Kd_Bryan_Lewis_deep    !< The abyssal value of a Bryan-Lewis diffusivity profile
                                    !! (m2/s)
  real    :: Kd_Bryan_Lewis_surface !< "The surface value of a Bryan-Lewis diffusivity profile
                                    !! (m2/s)
  real    :: Bryan_Lewis_depth_cent !<  The depth about which the transition in the Bryan-Lewis
                                    !! is centered (m)
  real    :: Bryan_Lewis_width_trans!< The width of the transition in the Bryan-Lewis profile (m)
  real    :: Kd_min                 !< minimum diapycnal diffusivity (m2/s)
  real    :: Kd                     !< interior diapycnal diffusivity (m2/s)
  real    :: N0_2Omega              !< ratio of the typical Buoyancy frequency to
                                    !! twice the Earth's rotation period, used with the
                                    !! Henyey scaling from the mixing
  real    :: prandtl_turb           !< Turbulent Prandtl number
  real    :: Kd_tanh_lat_scale      !< A nondimensional scaling for the range of
                                    !! diffusivities with Kd_tanh_lat_fn. Valid values
                                    !! are in the range of -2 to 2; 0.4 reproduces CM2M.
  real    :: Kdml                   !< mixed layer diapycnal diffusivity (m2/s)
                                    !! when bulkmixedlayer==.false.
  real    :: Hmix                   !< mixed layer thickness (meter) when
                                    !! bulkmixedlayer==.false.
  logical :: Kd_tanh_lat_fn         !< If true, use the tanh dependence of Kd_sfc on
                                    !! latitude, like GFDL CM2.1/CM2M.  There is no
                                    !! physical justification for this form, and it can
                                    !! not be used with Henyey_IGW_background.
  logical :: Bryan_Lewis_diffusivity!< If true, background vertical diffusivity
                                    !! uses Bryan-Lewis (1979) like tanh profile.
  logical :: Henyey_IGW_background  !< If true, use a simplified variant of the
             !! Henyey et al, JGR (1986) latitudinal scaling for the background diapycnal diffusivity,
             !! which gives a marked decrease in the diffusivity near the equator.  The simplification
             !! here is to assume that the in-situ stratification is the same as the reference stratificaiton.
  logical :: Henyey_IGW_background_new !< same as Henyey_IGW_background
             !! but incorporate the effect of stratification on TKE dissipation,
             !! e = f/f_0 * acosh(N/f) / acosh(N_0/f_0) * e_0
             !! where e is the TKE dissipation, and N_0 and f_0
             !! are the reference buoyancy frequency and inertial frequencies respectively.
             !! e_0 is the reference dissipation at (N_0,f_0). In the previous version, N=N_0.
             !! Additionally, the squared inverse relationship between  diapycnal diffusivities
             !! and stratification is included:
             !!
             !! kd = e/N^2
             !!
             !! where kd is the diapycnal diffusivity. This approach assumes that work done
             !! against gravity is uniformly distributed throughout the column. Whereas, kd=kd_0*e,
             !! as in the original version, concentrates buoyancy work in regions of strong stratification.
  logical :: bulkmixedlayer !< If true, a refined bulk mixed layer scheme is used
  logical :: debug !< If true, turn on debugging in this module
  ! Daignostic handles and pointers
  type(diag_ctrl), pointer :: diag => NULL()
  integer :: id_kd_bkgnd = -1, id_kv_bkgnd = -1

  ! Diagnostics arrays
  real, allocatable, dimension(:,:,:) :: kd_bkgnd !< Background diffusivity (m2/s)
  real, allocatable, dimension(:,:,:) :: kv_bkgnd !< Background viscosity  (m2/s)

end type bkgnd_mixing_cs

character(len=40)  :: mdl = "MOM_bkgnd_mixing" !< This module's name.

contains

!> Initialize the background mixing routine.
subroutine bkgnd_mixing_init(Time, G, GV, param_file, diag, CS)

  type(time_type),         intent(in)    :: Time       !< The current time.
  type(ocean_grid_type),   intent(in)    :: G          !< Grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure.
  type(param_file_type),   intent(in)    :: param_file !< Run-time parameter file handle
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics control structure.
  type(bkgnd_mixing_cs),    pointer        :: CS        !< This module's control structure.

  ! Local variables

! This include declares and sets the variable "version".
#include "version_variable.h"

  if (associated(CS)) then
    call MOM_error(WARNING, "bkgnd_mixing_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  ! Read parameters
  call log_version(param_file, mdl, version, &
    "Adding static vertical background mixing coefficients")

  call get_param(param_file, mdl, "KD", CS%Kd, &
                 "The background diapycnal diffusivity of density in the \n"//&
                 "interior. Zero or the molecular value, ~1e-7 m2 s-1, \n"//&
                 "may be used.", units="m2 s-1", fail_if_missing=.true.)

  call get_param(param_file, mdl, "KD_MIN", CS%Kd_min, &
                 "The minimum diapycnal diffusivity.", &
                 units="m2 s-1", default=0.01*CS%Kd)

  ! The following is needed to set one of the choices of vertical background mixing
  call get_param(param_file, mdl, "BULKMIXEDLAYER", CS%bulkmixedlayer, &
                 do_not_log=.true.)
  if (CS%bulkmixedlayer) then
    ! Check that Kdml is not set when using bulk mixed layer
    call get_param(param_file, mdl, "KDML", CS%Kdml, default=-1.)
    if (CS%Kdml>0.) call MOM_error(FATAL, &
                 "bkgnd_mixing_init: KDML cannot be set when using"// &
                 "bulk mixed layer.")
    CS%Kdml = CS%Kd ! This is not used with a bulk mixed layer, but also
                    ! cannot be a NaN.
  else
    call get_param(param_file, mdl, "KDML", CS%Kdml, &
                 "If BULKMIXEDLAYER is false, KDML is the elevated \n"//&
                 "diapycnal diffusivity in the topmost HMIX of fluid. \n"//&
                 "KDML is only used if BULKMIXEDLAYER is false.", &
                 units="m2 s-1", default=CS%Kd)
    call get_param(param_file, mdl, "HMIX_FIXED", CS%Hmix, &
                 "The prescribed depth over which the near-surface \n"//&
                 "viscosity and diffusivity are elevated when the bulk \n"//&
                 "mixed layer is not used.", units="m", fail_if_missing=.true.)
  endif


  call get_param(param_file, mdl, 'DEBUG', CS%debug, default=.False., do_not_log=.True.)

  call get_param(param_file, mdl, "PRANDTL_TURB", CS%prandtl_turb, &
                 units="nondim", default=1.0, do_not_log=.true.)

  call openParameterBlock(param_file,'MOM_BACKGROUND_MIXING')

  call get_param(param_file, mdl, "BRYAN_LEWIS_DIFFUSIVITY", &
                                CS%Bryan_Lewis_diffusivity, &
                 "If true, use a Bryan & Lewis (JGR 1979) like tanh \n"//&
                 "profile of background diapycnal diffusivity with depth. \n"//&
                 "This is done via CVMix.", default=.false.)

  if (CS%Bryan_Lewis_diffusivity) then

    call get_param(param_file, mdl, "KD_BRYAN_LEWIS_DEEP", &
                   CS%Kd_Bryan_Lewis_deep, &
                   "The abyssal value of a Bryan-Lewis diffusivity profile.", &
                   units="m2 s-1", fail_if_missing=.true.)

    call get_param(param_file, mdl, "KD_BRYAN_LEWIS_SURFACE", &
                   CS%Kd_Bryan_Lewis_surface, &
                   "The surface value of a Bryan-Lewis diffusivity profile.", &
                   units="m2 s-1", fail_if_missing=.true.)

    call get_param(param_file, mdl, "BRYAN_LEWIS_DEPTH_CENT", &
                   CS%Bryan_Lewis_depth_cent, &
                   "The depth about which the transition in the Bryan-Lewis.", &
                   units="m", fail_if_missing=.true.)

    call get_param(param_file, mdl, "BRYAN_LEWIS_WIDTH_TRANS", &
                   CS%Bryan_Lewis_width_trans, &
                   "The width of the transition in the Bryan-Lewis.",&
                   units="m", fail_if_missing=.true.)

  endif ! CS%Bryan_Lewis_diffusivity

  call get_param(param_file, mdl, "HENYEY_IGW_BACKGROUND", &
                                CS%Henyey_IGW_background, &
                 "If true, use a latitude-dependent scaling for the near \n"//&
                 "surface background diffusivity, as described in \n"//&
                 "Harrison & Hallberg, JPO 2008.", default=.false.)

  call get_param(param_file, mdl, "HENYEY_IGW_BACKGROUND_NEW", &
                                CS%Henyey_IGW_background_new, &
                 "If true, use a better latitude-dependent scaling for the\n"//&
                 "background diffusivity, as described in \n"//&
                 "Harrison & Hallberg, JPO 2008.", default=.false.)

  if (CS%Henyey_IGW_background .and. CS%Henyey_IGW_background_new) &
     call MOM_error(FATAL, "set_diffusivity_init: HENYEY_IGW_BACKGROUND and \n"//&
          "HENYEY_IGW_BACKGROUND_NEW are mutually exclusive. Set only one or none.")

  if (CS%Henyey_IGW_background) &
    call get_param(param_file, mdl, "HENYEY_N0_2OMEGA", CS%N0_2Omega, &
                  "The ratio of the typical Buoyancy frequency to twice \n"//&
                  "the Earth's rotation period, used with the Henyey \n"//&
                  "scaling from the mixing.", units="nondim", default=20.0)

  call get_param(param_file, mdl, "KD_TANH_LAT_FN", &
                  CS%Kd_tanh_lat_fn, &
                 "If true, use a tanh dependence of Kd_sfc on latitude, \n"//&
                 "like CM2.1/CM2M.  There is no physical justification \n"//&
                 "for this form, and it can not be used with \n"//&
                 "HENYEY_IGW_BACKGROUND.", default=.false.)

  if (CS%Kd_tanh_lat_fn) &
  call get_param(param_file, mdl, "KD_TANH_LAT_SCALE", &
                 CS%Kd_tanh_lat_scale, &
                 "A nondimensional scaling for the range ofdiffusivities \n"//&
                 "with KD_TANH_LAT_FN. Valid values are in the range of \n"//&
                 "-2 to 2; 0.4 reproduces CM2M.", units="nondim", default=0.0)

  if (CS%Henyey_IGW_background .and. CS%Kd_tanh_lat_fn) call MOM_error(FATAL, &
    "MOM_bkgnd_mixing: KD_TANH_LAT_FN can not be used with HENYEY_IGW_BACKGROUND.")

  call closeParameterBlock(param_file)

  ! allocate arrays and set them to zero
  allocate(CS%kd_bkgnd(SZI_(G), SZJ_(G), SZK_(G)+1)); CS%kd_bkgnd(:,:,:) = 0.
  allocate(CS%kv_bkgnd(SZI_(G), SZJ_(G), SZK_(G)+1)); CS%kv_bkgnd(:,:,:) = 0.

  ! Register diagnostics
  CS%diag => diag
  CS%id_kd_bkgnd = register_diag_field('ocean_model', 'bkgnd_kd', diag%axesTi, Time, &
      'Background diffusivity added by MOM_bkgnd_mixing module', 'm2/s')
  CS%id_kv_bkgnd = register_diag_field('ocean_model', 'bkgnd_kv', diag%axesTi, Time, &
      'Background viscosity added by MOM_bkgnd_mixing module', 'm2/s')

end subroutine bkgnd_mixing_init

!> Subroutine for calculating vertical background diffusivities/viscosities
subroutine calculate_bkgnd_mixing(h, tv, T_f, S_f, fluxes, G, GV, CS)

  type(ocean_grid_type),                      intent(in)  :: G  !< Grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV !< Vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)  :: h  !< Layer thickness, in m or kg m-2.
  type(thermo_var_ptrs),                      intent(in)  :: tv !< Thermodynamics structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)  :: T_f!< temperature (deg C), after massless
                                                                !! layers filled vertically by diffusion.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)  :: S_f!< salinity, after massless
                                                                !! layers filled vertically by diffusion.
  type(forcing),                        intent(in)    :: fluxes !< Structure of surface fluxes that may be
                                                                !! used.
  type(bkgnd_mixing_cs),                          pointer :: CS !< The control structure returned by
                                                                !! a previous call to bkgnd_mixing_init.

  ! local variables
  real, dimension(SZI_(G), SZJ_(G)) ::  Kd_sfc !< surface value of the diffusivity (m2/s)
  real, dimension(SZI_(G)) :: &
    depth_i, &   !< distance from surface of an interface (meter)
    N2_bot       !< bottom squared buoyancy frequency (1/s2)
  real, dimension(SZI_(G),SZK_(G)+1) :: &
    N2_int,   &  !< squared buoyancy frequency associated at interfaces (1/s2)
    dRho_int     !< locally ref potential density difference across interfaces (in s-2) smg: or kg/m3?
  real, dimension(SZI_(G),SZK_(G)) :: &
    N2_lay       !< squared buoyancy frequency associated with layers (1/s2)

  real, dimension(SZK_(G)) :: depth_k !< distance from surface of an interface (meter)
  real :: I_x30  !< 2/acos(2) = 1/(sin(30 deg) * acosh(1/sin(30 deg)))
  real :: deg_to_rad !< factor converting degrees to radians, pi/180.
  real :: abs_sin    !< absolute value of sine of latitude (nondim)
  real :: depth_c    !< depth of the center of a layer (meter)
  real :: I_Hmix     !< inverse of fixed mixed layer thickness (1/m)
  real :: epsilon
  real :: I_2Omega   !< 1/(2 Omega) (sec)
  real :: N_2Omega
  real :: N02_N2
  integer :: i, j, k, is, ie, js, je, nz

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = G%ke

  ! set some parameters
  deg_to_rad = atan(1.0)/45.0 ! = PI/180
  epsilon = 1.e-10

  if (CS%Bryan_Lewis_diffusivity) then
!$OMP parallel do default(none) shared(is,ie,js,je,nz, CS)
!$OMP                          private(cvmix_init_bkgnd,cvmix_coeffs_bkgnd)

    ! Bryan & Lewis is computed via CVMix
    do j=js,je; do i=is,ie

      depth_k(:) = 0.0
      do k=1,nz
        depth_k(k) = depth_k(k) + GV%H_to_m*h(i,j,k)
      enddo

        call cvmix_init_bkgnd(max_nlev=nz, &
                              zw = depth_k(:), &  !< interface depth, must be positive.
                              bl1 = CS%Kd_Bryan_Lewis_deep, &
                              bl2 = CS%Kd_Bryan_Lewis_surface, &
                              bl3 = 1.0/CS%Bryan_Lewis_depth_cent , &
                              bl4 = CS%Bryan_Lewis_width_trans, &
                              prandtl = CS%prandtl_turb)

        call cvmix_coeffs_bkgnd(Mdiff_out=CS%kv_bkgnd(i,j,:), &
                                Tdiff_out=CS%kd_bkgnd(i,j,:), &
                                nlev=nz, &
                                max_nlev=nz)
    enddo; enddo
  else
!$OMP parallel do default(none) shared(is,ie,js,je,CS,Kd_sfc)
    do j=js,je ; do i=is,ie
      Kd_sfc(i,j) = CS%Kd
    enddo ; enddo
  endif

  if (CS%Henyey_IGW_background) then
    I_x30 = 2.0 / invcosh(CS%N0_2Omega*2.0) ! This is evaluated at 30 deg.
!$OMP parallel do default(none)
!shared(is,ie,js,je,Kd_sfc,CS,G,deg_to_rad,epsilon,I_x30) &
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

!$OMP parallel do default(none)
!shared(is,ie,js,je,nz,G,GV,CS,h,tv,T_f,S_f,fluxes,dd, &
!$OMP
!Kd,Kd_sfc,epsilon,deg_to_rad,I_2Omega,visc,    &
!$OMP                                  Kd_int,dt,u,v,Omega2)   &
!$OMP
!private(dRho_int,I_trans,atan_fn_sfc,I_atan_fn,atan_fn_lay, &
!$OMP                                  I_Hmix,depth_c,depth,N2_lay, N2_int,
!N2_bot,        &
!$OMP                                  I_x30,abs_sin,N_2Omega,N02_N2,KT_extra,
!KS_extra,   &
!$OMP                                  TKE_to_Kd,maxTKE,dissip,kb)
  do j=js,je
    ! Set up variables related to the stratification.
    call find_N2(h, tv, T_f, S_f, fluxes, j, G, GV, dRho_int, N2_lay, N2_int, N2_bot)
    !if (associated(dd%N2_3d)) then
    !  do K=1,nz+1 ; do i=is,ie ; dd%N2_3d(i,j,K) = N2_int(i,K) ; enddo ; enddo
    !endif

    ! Set up the background diffusivity.
    if ((.not. CS%Bryan_Lewis_diffusivity) .and. (.not.CS%bulkmixedlayer) .and. &
       (CS%Kd/= CS%Kdml)) then

      I_Hmix = 1.0 / CS%Hmix
      do i=is,ie ; depth_i(i) = 0.0 ; enddo
      do k=1,nz ; do i=is,ie
        depth_c = depth_i(i) + 0.5*GV%H_to_m*h(i,j,k)

        if (depth_c <= CS%Hmix) then ; CS%kd_bkgnd(i,j,k) = CS%Kdml
        elseif (depth_c >= 2.0*CS%Hmix) then ; CS%kd_bkgnd(i,j,k) = Kd_sfc(i,j)
        else
          CS%kd_bkgnd(i,j,k) = ((Kd_sfc(i,j) - CS%Kdml) * I_Hmix) * depth_c + &
                      (2.0*CS%Kdml - Kd_sfc(i,j))
        endif

        depth_i(i) = depth_i(i) + GV%H_to_m*h(i,j,k)
      enddo ; enddo
    elseif (CS%Henyey_IGW_background_new) then
      I_x30 = 2.0 / invcosh(CS%N0_2Omega*2.0) ! This is evaluated at 30 deg.
      do k=1,nz ; do i=is,ie
        abs_sin = max(epsilon,abs(sin(G%geoLatT(i,j)*deg_to_rad)))
        N_2Omega = max(abs_sin,sqrt(N2_lay(i,k))*I_2Omega)
        N02_N2 = (CS%N0_2Omega/N_2Omega)**2
        CS%kd_bkgnd(i,j,k) = max(CS%Kd_min, Kd_sfc(i,j) * &
             ((abs_sin * invcosh(N_2Omega/abs_sin)) * I_x30)*N02_N2)
      enddo ; enddo

    else
      do k=1,nz ; do i=is,ie
        CS%kd_bkgnd(i,j,k) = Kd_sfc(i,j)
      enddo ; enddo
    endif
  enddo ! j-loop

  if (CS%debug) then
    call hchksum(Kd_sfc,"Kd_sfc",G%HI,haloshift=0)
    call hchksum(CS%kd_bkgnd, "MOM_bkgnd_mixing: kd_bkgnd",G%HI,haloshift=0)
    call hchksum(CS%kv_bkgnd, "MOM_bkgnd_mixing: kv_bkgnd",G%HI,haloshift=0)
  endif

  ! send diagnostics to post_data
  if (CS%id_kd_bkgnd > 0) call post_data(CS%id_kd_bkgnd, CS%kd_bkgnd, CS%diag)
  if (CS%id_kv_bkgnd > 0) call post_data(CS%id_kv_bkgnd, CS%kv_bkgnd, CS%diag)

end subroutine calculate_bkgnd_mixing


!> Computes N2
subroutine find_N2(h, tv, T_f, S_f, fluxes, j, G, GV, dRho_int, &
                   N2_lay, N2_int, N2_bot)
  type(ocean_grid_type),                    intent(in)   :: G    !< The ocean's grid structure
  type(verticalGrid_type),                  intent(in)   :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)   :: h    !< Layer thicknesses, in H (usually m or kg m-2)
  type(thermo_var_ptrs),                    intent(in)   :: tv
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)   :: T_f, S_f
  type(forcing),                            intent(in)   :: fluxes
  integer,                                  intent(in)   :: j
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
    h_amp(i) = 0.0
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

!> Reads the parameter "USE_CVMIX_BACKGROUND" and returns state.
!! This function allows other modules to know whether this parameterization will
!! be used without needing to duplicate the log entry.
logical function cvmix_bkgnd_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  call get_param(param_file, mdl, "USE_CVMIX_BACKGROUND", cvmix_bkgnd_is_used, &
                 default=.false., do_not_log = .true.)

end function cvmix_bkgnd_is_used

!> Clear pointers and dealocate memory
subroutine bkgnd_mixing_end(CS)
  type(bkgnd_mixing_cs), pointer :: CS ! Control structure

  deallocate(CS%kd_bkgnd)
  deallocate(CS%kv_bkgnd)
  deallocate(CS)

end subroutine bkgnd_mixing_end


end module MOM_bkgnd_mixing
