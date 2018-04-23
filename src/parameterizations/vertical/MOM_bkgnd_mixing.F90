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
use MOM_error_handler,   only : is_root_pe
use MOM_file_parser,     only : openParameterBlock, closeParameterBlock
use MOM_debugging,       only : hchksum
use MOM_grid,            only : ocean_grid_type
use MOM_verticalGrid,    only : verticalGrid_type
use MOM_file_parser,     only : get_param, log_version, param_file_type
use CVMix_background,    only : CVMix_init_bkgnd, CVMix_coeffs_bkgnd
use MOM_variables,       only : vertvisc_type
use MOM_intrinsic_functions, only : invcosh

implicit none ; private

#include <MOM_memory.h>

public bkgnd_mixing_init
public bkgnd_mixing_end
public calculate_bkgnd_mixing
public sfc_bkgnd_mixing

!> Control structure including parameters for this module.
type, public :: bkgnd_mixing_cs

  ! Parameters
  real    :: Bryan_Lewis_c1         !< The vertical diffusivity values for  Bryan-Lewis profile
                                    !! at |z|=D (m2/s)
  real    :: Bryan_Lewis_c2         !< The amplitude of variation in diffusivity for the
                                    !! Bryan-Lewis diffusivity profile (m2/s)
  real    :: Bryan_Lewis_c3         !< The inverse length scale for transition region in the
                                    !! Bryan-Lewis diffusivity profile (1/m)
  real    :: Bryan_Lewis_c4         !< The depth where diffusivity is Bryan_Lewis_bl1 in the
                                    !! Bryan-Lewis profile (m)
  real    :: Kd_min                 !< minimum diapycnal diffusivity (m2/s)
  real    :: Kd                     !< interior diapycnal diffusivity (m2/s)
  real    :: N0_2Omega              !< ratio of the typical Buoyancy frequency to
                                    !! twice the Earth's rotation period, used with the
                                    !! Henyey scaling from the mixing
  real    :: prandtl_bkgnd          !< Turbulent Prandtl number used to convert
                                    !! vertical background diffusivity into viscosity
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

  real, allocatable, dimension(:,:)   ::  Kd_sfc !< surface value of the diffusivity (m2/s)
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

  ! BULKMIXEDLAYER is not always defined (e.g., CM2G63L), so the following by pass
  ! the need to include BULKMIXEDLAYER in MOM_input
  CS%bulkmixedlayer = (GV%nkml > 0)
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

  call get_param(param_file, mdl, "PRANDTL_BKGND", CS%prandtl_bkgnd, &
                 "Turbulent Prandtl number used to convert vertical \n"//&
                 "background diffusivities into viscosities.", &
                 units="nondim", default=1.0)

!  call openParameterBlock(param_file,'MOM_BACKGROUND_MIXING')

  call get_param(param_file, mdl, "BRYAN_LEWIS_DIFFUSIVITY", &
                                CS%Bryan_Lewis_diffusivity, &
                 "If true, use a Bryan & Lewis (JGR 1979) like tanh \n"//&
                 "profile of background diapycnal diffusivity with depth. \n"//&
                 "This is done via CVMix.", default=.false.)

  if (CS%Bryan_Lewis_diffusivity) then

    call get_param(param_file, mdl, "BRYAN_LEWIS_C1", &
                   CS%Bryan_Lewis_c1, &
                   "The vertical diffusivity values for Bryan-Lewis profile at |z|=D.", &
                   units="m2 s-1", fail_if_missing=.true.)

    call get_param(param_file, mdl, "BRYAN_LEWIS_C2", &
                   CS%Bryan_Lewis_c2, &
                   "The amplitude of variation in diffusivity for the Bryan-Lewis profile", &
                   units="m2 s-1", fail_if_missing=.true.)

    call get_param(param_file, mdl, "BRYAN_LEWIS_C3", &
                   CS%Bryan_Lewis_c3, &
                   "The inverse length scale for transition region in the Bryan-Lewis profile", &
                   units="m-1", fail_if_missing=.true.)

    call get_param(param_file, mdl, "BRYAN_LEWIS_C4", &
                   CS%Bryan_Lewis_c4, &
                   "The depth where diffusivity is BRYAN_LEWIS_C1 in the Bryan-Lewis profile",&
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

!  call closeParameterBlock(param_file)

  ! allocate arrays and set them to zero
  allocate(CS%kd_bkgnd(SZI_(G), SZJ_(G), SZK_(G)+1)); CS%kd_bkgnd(:,:,:) = 0.
  allocate(CS%kv_bkgnd(SZI_(G), SZJ_(G), SZK_(G)+1)); CS%kv_bkgnd(:,:,:) = 0.
  allocate(CS%Kd_sfc(SZI_(G), SZJ_(G))); CS%Kd_sfc(:,:) = 0.

  ! Register diagnostics
  CS%diag => diag
  CS%id_kd_bkgnd = register_diag_field('ocean_model', 'Kd_bkgnd', diag%axesTi, Time, &
      'Background diffusivity added by MOM_bkgnd_mixing module', 'm2/s')
  CS%id_kv_bkgnd = register_diag_field('ocean_model', 'Kv_bkgnd', diag%axesTi, Time, &
      'Background viscosity added by MOM_bkgnd_mixing module', 'm2/s')

end subroutine bkgnd_mixing_init

!> Get surface vertical background diffusivities/viscosities.
subroutine sfc_bkgnd_mixing(G, CS)

  type(ocean_grid_type),                      intent(in)  :: G  !< Grid structure.
  type(bkgnd_mixing_cs),           pointer, intent(inout) :: CS !< The control structure returned by
                                                                !! a previous call to bkgnd_mixing_init.
  ! local variables
  real :: I_x30  !< 2/acos(2) = 1/(sin(30 deg) * acosh(1/sin(30 deg)))
  real :: deg_to_rad !< factor converting degrees to radians, pi/180.
  real :: abs_sin    !< absolute value of sine of latitude (nondim)
  real :: epsilon
  integer :: i, j, k, is, ie, js, je

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec

  ! set some parameters
  deg_to_rad = atan(1.0)/45.0 ! = PI/180
  epsilon = 1.e-10


  if (.not. CS%Bryan_Lewis_diffusivity) then
!$OMP parallel do default(none) shared(is,ie,js,je,CS,Kd_sfc)
    do j=js,je ; do i=is,ie
      CS%Kd_sfc(i,j) = CS%Kd
    enddo ; enddo
  endif

  if (CS%Henyey_IGW_background) then
    I_x30 = 2.0 / invcosh(CS%N0_2Omega*2.0) ! This is evaluated at 30 deg.
!$OMP parallel do default(none)
!shared(is,ie,js,je,Kd_sfc,CS,G,deg_to_rad,epsilon,I_x30) &
!$OMP                          private(abs_sin)
    do j=js,je ; do i=is,ie
      abs_sin = abs(sin(G%geoLatT(i,j)*deg_to_rad))
      CS%Kd_sfc(i,j) = max(CS%Kd_min, CS%Kd_sfc(i,j) * &
           ((abs_sin * invcosh(CS%N0_2Omega/max(epsilon,abs_sin))) * I_x30) )
    enddo ; enddo
  elseif (CS%Kd_tanh_lat_fn) then
!$OMP parallel do default(none) shared(is,ie,js,je,Kd_sfc,CS,G)
    do j=js,je ; do i=is,ie
      !   The transition latitude and latitude range are hard-scaled here, since
      ! this is not really intended for wide-spread use, but rather for
      ! comparison with CM2M / CM2.1 settings.
      CS%Kd_sfc(i,j) = max(CS%Kd_min, CS%Kd_sfc(i,j) * (1.0 + &
          CS%Kd_tanh_lat_scale * 0.5*tanh((abs(G%geoLatT(i,j)) - 35.0)/5.0) ))
    enddo ; enddo
  endif

  if (CS%debug) call hchksum(CS%Kd_sfc,"After sfc_bkgnd_mixing: Kd_sfc",G%HI,haloshift=0)

end subroutine sfc_bkgnd_mixing


!> Calculates the vertical background diffusivities/viscosities
subroutine calculate_bkgnd_mixing(h, tv, N2_lay, kd_lay, kv, j, G, GV, CS)

  type(ocean_grid_type),                      intent(in)  :: G   !< Grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)  :: h   !< Layer thickness, in m or kg m-2.
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamics structure.
  real, dimension(SZI_(G),SZK_(G)),           intent(in)  :: N2_lay!< squared buoyancy frequency associated
                                                                 !! with layers (1/s2)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: kd_lay!< Diapycnal diffusivity of each layer m2 s-1.
  real, dimension(:,:,:),                     pointer     :: kv  !< The "slow" vertical viscosity at each interface
                                                                 !! (not layer!) in m2 s-1.
  integer,                                   intent(in)   :: j   !< Meridional grid indice.
  type(bkgnd_mixing_cs),                          pointer :: CS  !< The control structure returned by
                                                                 !! a previous call to bkgnd_mixing_init.

  ! local variables
  real, dimension(SZI_(G), SZK_(G)+1) :: depth_2d  !< distance from surface of an interface (m)
  real, dimension(SZI_(G)) :: &
        depth        !< distance from surface of an interface (meter)
  real :: depth_c    !< depth of the center of a layer (meter)
  real :: I_Hmix     !< inverse of fixed mixed layer thickness (1/m)
  real :: I_2Omega   !< 1/(2 Omega) (sec)
  real :: N_2Omega
  real :: N02_N2
  real :: I_x30  !< 2/acos(2) = 1/(sin(30 deg) * acosh(1/sin(30 deg)))
  real :: deg_to_rad !< factor converting degrees to radians, pi/180.
  real :: abs_sin    !< absolute value of sine of latitude (nondim)
  real :: epsilon
  integer :: i, k, is, ie, js, je, nz

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = G%ke

  ! set some parameters
  deg_to_rad = atan(1.0)/45.0 ! = PI/180
  epsilon = 1.e-10

  depth_2d(:,:) = 0.0
  ! Set up the background diffusivity.
  if (CS%Bryan_Lewis_diffusivity) then

    do i=is,ie
      do k=2,nz+1
        depth_2d(i,k) = depth_2d(i,k-1) + GV%H_to_m*h(i,j,k-1)
      enddo
      ! if (is_root_pe()) write(*,*)'depth_3d(i,j,:)',depth_3d(i,j,:)

      call CVMix_init_bkgnd(max_nlev=nz, &
                            zw = depth_2d(i,:), &  !< interface depth, must bepositive.
                            bl1 = CS%Bryan_Lewis_c1, &
                            bl2 = CS%Bryan_Lewis_c2, &
                            bl3 = CS%Bryan_Lewis_c3, &
                            bl4 = CS%Bryan_Lewis_c4, &
                            prandtl = CS%prandtl_bkgnd)

      call CVMix_coeffs_bkgnd(Mdiff_out=CS%kv_bkgnd(i,j,:), &
                              Tdiff_out=CS%kd_bkgnd(i,j,:), &
                              nlev=nz, &
                              max_nlev=nz)

      ! Update Kd
      do k=1,nz
        kd_lay(i,j,k) = kd_lay(i,j,k) + 0.5*(CS%kd_bkgnd(i,j,K) + CS%kd_bkgnd(i,j,K+1))
      enddo
    enddo ! i loop

  elseif ((.not. CS%Bryan_Lewis_diffusivity) .and. (.not.CS%bulkmixedlayer) .and. &
     (CS%Kd/= CS%Kdml)) then
    I_Hmix = 1.0 / CS%Hmix
    do i=is,ie ; depth(i) = 0.0 ; enddo
    do k=1,nz ; do i=is,ie
      depth_c = depth(i) + 0.5*GV%H_to_m*h(i,j,k)
      if (depth_c <= CS%Hmix) then ; CS%kd_bkgnd(i,j,k) = CS%Kdml
      elseif (depth_c >= 2.0*CS%Hmix) then ; CS%kd_bkgnd(i,j,k) = CS%Kd_sfc(i,j)
      else
        kd_lay(i,j,k) = ((CS%Kd_sfc(i,j) - CS%Kdml) * I_Hmix) * depth_c + &
                    (2.0*CS%Kdml - CS%Kd_sfc(i,j))
      endif

      depth(i) = depth(i) + GV%H_to_m*h(i,j,k)
    enddo ; enddo

  elseif (CS%Henyey_IGW_background_new) then
    I_x30 = 2.0 / invcosh(CS%N0_2Omega*2.0) ! This is evaluated at 30 deg.
    do k=1,nz ; do i=is,ie
      abs_sin = max(epsilon,abs(sin(G%geoLatT(i,j)*deg_to_rad)))
      N_2Omega = max(abs_sin,sqrt(N2_lay(i,k))*I_2Omega)
      N02_N2 = (CS%N0_2Omega/N_2Omega)**2
      kd_lay(i,j,k) = max(CS%Kd_min, CS%Kd_sfc(i,j) * &
           ((abs_sin * invcosh(N_2Omega/abs_sin)) * I_x30)*N02_N2)
    enddo ; enddo

  else
    do k=1,nz ; do i=is,ie
      kd_lay(i,j,k) = CS%Kd_sfc(i,j)
    enddo ; enddo
  endif

  ! Update CS%kd_bkgnd and CS%kv_bkgnd for diagnostic purposes
  if (.not. CS%Bryan_Lewis_diffusivity) then
    do i=is,ie
      CS%kd_bkgnd(i,j,1) = 0.0; CS%kv_bkgnd(i,j,1) = 0.0
      CS%kd_bkgnd(i,j,nz+1) = 0.0; CS%kv_bkgnd(i,j,nz+1) = 0.0
      do k=2,nz
        CS%kd_bkgnd(i,j,k) = CS%kd_bkgnd(i,j,k) + 0.5*(kd_lay(i,j,K-1) + kd_lay(i,j,K))
        CS%kv_bkgnd(i,j,k) = CS%kd_bkgnd(i,j,k) * CS%prandtl_bkgnd
      enddo
    enddo
  endif

  ! Update kv
  if (associated(kv)) then
    do i=is,ie
      do k=1,nz+1
        kv(i,j,k) = kv(i,j,k) + CS%kv_bkgnd(i,j,k)
      enddo
    enddo
  endif

end subroutine calculate_bkgnd_mixing

!> Reads the parameter "USE_CVMix_BACKGROUND" and returns state.
!! This function allows other modules to know whether this parameterization will
!! be used without needing to duplicate the log entry.
logical function CVMix_bkgnd_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  call get_param(param_file, mdl, "USE_CVMix_BACKGROUND", CVMix_bkgnd_is_used, &
                 default=.false., do_not_log = .true.)

end function CVMix_bkgnd_is_used

!> Clear pointers and dealocate memory
subroutine bkgnd_mixing_end(CS)
  type(bkgnd_mixing_cs), pointer :: CS ! Control structure

  deallocate(CS%kd_bkgnd)
  deallocate(CS%kv_bkgnd)
  deallocate(CS)

end subroutine bkgnd_mixing_end


end module MOM_bkgnd_mixing
