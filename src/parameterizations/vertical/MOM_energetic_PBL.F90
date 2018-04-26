module MOM_energetic_PBL

! This file is part of MOM6. See LICENSE.md for the license.

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, 2015.                                          *
!*                                                                     *
!*    This file contains the subroutine (energetic_PBL) that uses an   *
!*  integrated boundary layer energy budget (like a bulk- or refined-  *
!*  bulk mixed layer scheme), but instead of homogenizing this model   *
!*  calculates a finite diffusivity and viscosity, which in this       *
!*  regard is conceptually similar to what is done with KPP or various *
!*  two-equation closures.  However, the scheme that is implemented    *
!*  here has the big advantage that is entirely implicit, but is       *
!*  simple enough that it requires only a single vertical pass to      *
!*  determine the diffusivity. The development of bulk mixed layer     *
!*  models stems from the work of various people, as described in the  *
!*  review paper by Niiler and Kraus (1979). The work here draws in    *
!*  with particular on the form for TKE decay proposed by Oberhuber    *
!*  (JPO, 1993, 808-829), with an extension to a refined bulk mixed    *
!*  layer as described in Hallberg (Aha Huliko'a, 2003).  The physical *
!*  processes portrayed in this subroutine include convectively driven *
!*  mixing and mechanically driven mixing.  Unlike boundary-layer      *
!*  mixing, stratified shear mixing is not a one-directional turbulent *
!*  process, and it is dealt with elsewhere in the MOM6 code within    *
!*  the module MOM_kappa_shear.F90.  It is assumed that the heat,      *
!*  mass, and salt fluxes have been applied elsewhere, but that their  *
!*  implications for the integrated TKE budget have been captured in   *
!*  an array that is provided as an argument to this subroutine.  This *
!*  is a full 3-d array due to the effects of penetrating shortwave    *
!*  radiation.                                                         *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, T, S, Kd, etc.                        *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_domains,       only : create_group_pass, do_group_pass, group_pass_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
! use MOM_EOS, only : calculate_density, calculate_density_derivs
use cvmix_energetic_pbl, only : cvmix_energetic_PBL_CS, cvmix_epbl_init, cvmix_epbl_column, cvmix_epbl_end

implicit none ; private

#include <MOM_memory.h>

public energetic_PBL, energetic_PBL_init, energetic_PBL_end
public energetic_PBL_get_MLD
public energetic_PBL_CS

type, public :: energetic_PBL_CS ; private
  real    :: mstar           ! The ratio of the friction velocity cubed to the
                             ! TKE available to drive entrainment, nondimensional.
                             ! This quantity is the vertically integrated
                             ! shear production minus the vertically integrated
                             ! dissipation of TKE produced by shear.
  real    :: nstar           ! The fraction of the TKE input to the mixed layer
                             ! available to drive entrainment, nondim.
                             ! This quantity is the vertically integrated
                             ! buoyancy production minus the vertically integrated
                             ! dissipation of TKE produced by buoyancy.
  real    :: MixLenExponent  ! Exponent in the mixing length shape-function.
                             ! 1 is law-of-the-wall at top and bottom,
                             ! 2 is more KPP like.
  real    :: TKE_decay       ! The ratio of the natural Ekman depth to the TKE
                             ! decay scale, nondimensional.
  real    :: MKE_to_TKE_effic ! The efficiency with which mean kinetic energy
                             ! released by mechanically forced entrainment of
                             ! the mixed layer is converted to TKE, nondim.
!  real    :: Hmix_min        ! The minimum mixed layer thickness in m.
  real    :: ustar_min       ! A minimum value of ustar to avoid numerical
                             ! problems, in m s-1.  If the value is small enough,
                             ! this should not affect the solution.
  real    :: omega           !   The Earth's rotation rate, in s-1.
  real    :: omega_frac      !   When setting the decay scale for turbulence, use
                             ! this fraction of the absolute rotation rate blended
                             ! with the local value of f, as sqrt((1-of)*f^2 + of*4*omega^2).
  real    :: wstar_ustar_coef ! A ratio relating the efficiency with which
                             ! convectively released energy is converted to a
                             ! turbulent velocity, relative to mechanically
                             ! forced turbulent kinetic energy, nondim. Making
                             ! this larger increases the diffusivity.
  real    :: vstar_scale_fac ! An overall nondimensional scaling factor
                             ! for vstar.  Making this larger increases the
                             ! diffusivity.
  real    :: Ekman_scale_coef ! A nondimensional scaling factor controlling
                             ! the inhibition of the diffusive length scale by
                             ! rotation.  Making this larger decreases the
                             ! diffusivity in the planetary boundary layer.
  real    :: transLay_scale  ! A scale for the mixing length in the transition layer
                             ! at the edge of the boundary layer as a fraction of the
                             ! boundary layer thickness.  The default is 0, but a
                             ! value of 0.1 might be better justified by observations.
  real    :: MLD_tol         ! A tolerance for determining the boundary layer
                             ! thickness when Use_MLD_iteration is true, in m.
  real    :: min_mix_len     ! The minimum mixing length scale that will be
                             ! used by ePBL, in m.  The default (0) does not
                             ! set a minimum.
  real    :: N2_Dissipation_Scale_Neg
  real    :: N2_Dissipation_Scale_Pos
                             ! A nondimensional scaling factor controlling the
                             ! loss of TKE due to enhanced dissipation in the presence
                             ! of stratification.  This dissipation is applied to the
                             ! available TKE which includes both that generated at the
                             ! surface and that generated at depth.  It may be important
                             ! to distinguish which TKE flavor that this dissipation
                             ! applies to in subsequent revisions of this code.
                             ! "_Neg" and "_Pos" refer to which scale is applied as a
                             ! function of negative or positive local buoyancy.
  real    :: MSTAR_CAP       ! Since MSTAR is restoring undissipated energy to mixing,
                             ! there must be a cap on how large it can be.  This
                             ! is definitely a function of latitude (Ekman limit),
                             ! but will be taken as constant for now.
  real    :: MSTAR_SLOPE     ! Slope of the function which relates the shear production
                             ! to the mixing layer depth, Ekman depth, and Monin-Obukhov
                             ! depth.
  real    :: MSTAR_XINT      ! Value where MSTAR function transitions from linear
                             ! to decay toward MSTAR->0 at fully developed Ekman depth.
  real    :: MSTAR_XINT_UP   ! Similar but for transition to asymptotic cap.
  real    :: MSTAR_AT_XINT   ! Intercept value of MSTAR at value where function
                             ! changes to linear transition.
  integer :: LT_ENHANCE_FORM ! Integer for Enhancement functional form (various options)
  real    :: LT_ENHANCE_COEF ! Coefficient in fit for Langmuir Enhancment
  real    :: LT_ENHANCE_EXP  ! Exponent in fit for Langmuir Enhancement
  real :: MSTAR_N = -2.      ! Exponent in decay at negative and positive limits of MLD_over_STAB
  real :: MSTAR_A,MSTAR_A2   ! MSTAR_A and MSTAR_B are coefficients in asymptote toward limits.
  real :: MSTAR_B,MSTAR_B2   !  These are computed to match the function value and slope at both
                             !  ends of the linear fit within the well constrained region.
  real :: C_EK = 0.17        ! MSTAR Coefficient in rotation limit for mstar_mode=2
  real :: MSTAR_COEF = 0.3   ! MSTAR coefficient in rotation/stabilizing balance for mstar_mode=2
  real :: LaC_MLDoEK         ! Coefficients for Langmuir number modification based on
  real :: LaC_MLDoOB_stab    !  length scale ratios, MLD is boundary, EK is Ekman,
  real :: LaC_EKoOB_stab     !  and OB is Obukhov, the "o" in the name is for division.
  real :: LaC_MLDoOB_un      !  Stab/un are for stable (pos) and unstable (neg) Obukhov depths
  real :: LaC_EKoOB_un       !   ...
  real :: LaDepthRatio=0.04  ! The ratio of OBL depth to average Stokes drift over
  real :: Max_Enhance_M = 5. ! The maximum allowed LT enhancement to the mixing.
  real :: CNV_MST_FAC        ! Factor to reduce mstar when statically unstable.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.

  integer :: MSTAR_MODE = 0  ! An integer to determine which formula is used to
                             !  set mstar
  integer :: CONST_MSTAR=0,MLD_o_OBUKHOV=1,EKMAN_o_OBUKHOV=2
  logical :: MSTAR_FLATCAP=.true. !Set false to use asymptotic mstar cap.
  logical :: TKE_diagnostics = .false.
  logical :: Use_LA_windsea = .false.
  logical :: orig_PE_calc = .true.
  logical :: Use_MLD_iteration=.false. ! False to use old ePBL method.
  logical :: Orig_MLD_iteration=.false. ! False to use old MLD value
  logical :: MLD_iteration_guess=.false. ! False to default to guessing half the
                                         ! ocean depth for the iteration.
  logical :: Mixing_Diagnostics = .false. ! Will be true when outputing mixing
                                          !  length and velocity scale
  logical :: MSTAR_Diagnostics=.false.
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.

! These are terms in the mixed layer TKE budget, all in J m-2 = kg s-2.
  real, allocatable, dimension(:,:) :: &
    diag_TKE_wind, &   ! The wind source of TKE.
    diag_TKE_MKE, &    ! The resolved KE source of TKE.
    diag_TKE_conv, &   ! The convective source of TKE.
    diag_TKE_forcing, & ! The TKE sink required to mix surface
                       ! penetrating shortwave heating.
    diag_TKE_mech_decay, & ! The decay of mechanical TKE.
    diag_TKE_conv_decay, & ! The decay of convective TKE.
    diag_TKE_mixing,&  ! The work done by TKE to deepen
                       ! the mixed layer.
    ! Additional output parameters also 2d
    ML_depth, &        ! The mixed layer depth in m. (result after iteration step)
    ML_depth2, &       ! The mixed layer depth in m. (guess for iteration step)
    Enhance_M, &       ! The enhancement to the turbulent velocity scale (non-dim)
    MSTAR_MIX, &       ! Mstar used in EPBL
    MSTAR_LT, &        ! Mstar for Langmuir turbulence
    MLD_EKMAN, &       ! MLD over Ekman length
    MLD_OBUKHOV, &     ! MLD over Obukhov length
    EKMAN_OBUKHOV, &   ! Ekman over Obukhov length
    LA, &              ! Langmuir number
    LA_MOD             ! Modified Langmuir number

  real, allocatable, dimension(:,:,:) :: &
    Velocity_Scale, & ! The velocity scale used in getting Kd
    Mixing_Length     ! The length scale used in getting Kd
  integer :: id_ML_depth = -1, id_TKE_wind = -1, id_TKE_mixing = -1
  integer :: id_TKE_MKE = -1, id_TKE_conv = -1, id_TKE_forcing = -1
  integer :: id_TKE_mech_decay = -1, id_TKE_conv_decay = -1
  integer :: id_Hsfc_used = -1
  integer :: id_Mixing_Length = -1, id_Velocity_Scale = -1
  integer :: id_OSBL = -1, id_LT_Enhancement = -1, id_MSTAR_mix = -1
  integer :: id_mld_ekman = -1, id_mld_obukhov = -1, id_ekman_obukhov = -1
  integer :: id_LA_mod = -1, id_LA = -1, id_MSTAR_LT = -1

  type(cvmix_energetic_PBL_CS), pointer  :: cvmix_CS => null() ! A structure that is used to drive the cvmix kernel for kappa_shear
  
end type energetic_PBL_CS

type(diag_ctrl), pointer :: MY_diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
integer :: num_msg = 0, max_msg = 2

contains

!>    This subroutine determines the diffusivities from the integrated energetics
!!  mixed layer model.  It assumes that heating, cooling and freshwater fluxes
!!  have already been applied.  All calculations are done implicitly, and there
!!  is no stability limit on the time step.
subroutine energetic_PBL(h_3d, u_3d, v_3d, tv, fluxes, dt, Kd_int, G, GV, CS, &
                         dSV_dT, dSV_dS, TKE_forced, Buoy_Flux, dt_diag, last_call, &
                         dT_expected, dS_expected )
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h_3d   !< Layer thickness, in m or kg m-2.
                                                   !! (Intent in/out) The units of h are referred
                                                   !! to as H below.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u_3d   !< Zonal velocities interpolated to h points,
                                                   !! m s-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: v_3d   !< Zonal velocities interpolated to h points,
                                                   !! m s-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dSV_dT !< The partial derivative of in-situ specific
                                                   !! volume with potential temperature,
                                                   !! in m3 kg-1 K-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dSV_dS !< The partial derivative of in-situ specific
                                                   !! volume with salinity, in m3 kg-1 ppt-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: TKE_forced !< The forcing requirements to homogenize the
                                                   !! forcing that has been applied to each layer
                                                   !! through each layer, in J m-2.
  type(thermo_var_ptrs),   intent(in)    :: tv     !< A structure containing pointers to any
                                                   !! available thermodynamic fields. Absent fields
                                                   !! have NULL ptrs.
  type(forcing),           intent(in)    :: fluxes !< A structure containing pointers to any
                                                   !! possible forcing fields. Unused fields have
                                                   !! NULL ptrs.
  real,                    intent(in)    :: dt     !< Time increment, in s.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(out)   :: Kd_int !< The diagnosed diffusivities at interfaces,
                                                   !! in m2 s-1.
  type(energetic_PBL_CS),  pointer       :: CS     !< The control structure returned by a previous
                                                   !! call to mixedlayer_init.
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: Buoy_Flux !< The surface buoyancy flux. in m2/s3.
  real,          optional, intent(in)    :: dt_diag   !< The diagnostic time step, which may be less
                                                   !! than dt if there are two callse to
                                                   !! mixedlayer, in s.
  logical,       optional, intent(in)    :: last_call !< If true, this is the last call to
                                                   !! mixedlayer in the current time step, so
                                                   !! diagnostics will be written. The default
                                                   !! is .true.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: dT_expected, dS_expected

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
! These include mstar, nstar, TKE_decay, and conv_decay.  For the Oberhuber (1993) mixed layer,
! the values of these are:
!      mstar = 1.25,  nstar = 1, TKE_decay = 2.5, conv_decay = 0.5
!  TKE_decay is 1/kappa in eq. 28 of Oberhuber (1993), while
!  conv_decay is 1/mu.
!    For a traditional Kraus-Turner mixed layer, the values are:
!      mstar = 1.25, nstar = 0.4, TKE_decay = 0.0, conv_decay = 0.0

! Arguments: h_3d - Layer thickness, in m or kg m-2. (Intent in/out)
!                   The units of h are referred to as H below.
!  (in)      u_3d - Zonal velocities interpolated to h points, m s-1.
!  (in)      v_3d - Zonal velocities interpolated to h points, m s-1.
!  (in/out)  tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      dt - Time increment, in s.
!  (out)     Kd_int - The diagnosed diffusivities at interfaces, in m2 s-1.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 mixedlayer_init.
!  (in)      dSV_dT - The partial derivative of in-situ specific volume with
!                     potential temperature, in m3 kg-1 K-1.
!  (in)      dSV_dS - The partial derivative of in-situ specific volume with
!                     salinity, in m3 kg-1 ppt-1.
!  (in)      TKE_forced - The forcing requirements to homogenize the forcing
!                 that has been applied to each layer through each layer, in J m-2.
!  (in)      Buoy_Flux - The surface buoyancy flux. in m2/s3.
!  (in,opt)  dt_diag - The diagnostic time step, which may be less than dt
!                      if there are two callse to mixedlayer, in s.
!  (in,opt)  last_call - if true, this is the last call to mixedlayer in the
!                        current time step, so diagnostics will be written.
!                        The default is .true.
  real, dimension(SZI_(G),SZK_(GV)+1) :: Kd  ! The diapycnal diffusivity, in m2 s-1.
  real, dimension(SZI_(G)) :: &
    mech_TKE, &     !   The mechanically generated turbulent kinetic energy
                    ! available for mixing over a time step, in J m-2 = kg s-2.
    conv_PErel, &   ! The potential energy that has been convectively released
                    ! during this timestep, in J m-2 = kg s-2. A portion nstar_FC
                    ! of conv_PErel is available to drive mixing.
    absf            ! The absolute value of f, in s-1.
  real, dimension(SZI_(G),SZK_(GV)) :: &
    h, &            !   The layer thickness, in H (usually m or kg m-2).
    T, &            !   The layer temperatures, in deg C.
    S, &            !   The layer salinities, in psu.
    u, &            !   The zonal velocity, in m s-1.
    v               !   The meridional velocity, in m s-1.
  real, dimension(SZI_(G),SZJ_(G)) ::   Hsfc_used     !   The thickness of the surface region after buffer layer
  real    :: U_Star,U_Star_Mean
! The following is only used as a diagnostic.
  real :: dt__diag  ! A copy of dt_diag (if present) or dt, in s.
  real :: IdtdR0    !  = 1.0 / (dt__diag * Rho0), in m3 kg-1 s-1.
  logical :: write_diags  ! If true, write out diagnostics with this step.
  logical :: reset_diags  ! If true, zero out the accumulated diagnostics.
                ! detrainment, in units of m.

  integer :: i, j, k, is, ie, js, je, nz, itt, max_itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "energetic_PBL: "//&
         "Module must be initialized before it is used.")

  if (.not. ASSOCIATED(tv%eqn_of_state)) call MOM_error(FATAL, &
      "energetic_PBL: Temperature, salinity and an equation of state "//&
      "must now be used.")
  if (.NOT. ASSOCIATED(fluxes%ustar)) call MOM_error(FATAL, &
      "energetic_PBL: No surface TKE fluxes (ustar) defined in mixedlayer!")

  dt__diag = dt ; if (present(dt_diag)) dt__diag = dt_diag
  IdtdR0 = 1.0 / (dt__diag * GV%Rho0)

  write_diags = .true. ; if (present(last_call)) write_diags = last_call

  ! Determine whether to zero out diagnostics before accumulation.
  reset_diags = .true.
  if (present(dt_diag) .and. write_diags .and. (dt__diag > dt)) &
    reset_diags = .false.  ! This is the second call to mixedlayer.

  if (reset_diags) then
    if (CS%TKE_diagnostics) then
      do j=js,je ; do i=is,ie
        CS%diag_TKE_wind(i,j) = 0.0 ; CS%diag_TKE_MKE(i,j) = 0.0
        CS%diag_TKE_conv(i,j) = 0.0 ; CS%diag_TKE_forcing(i,j) = 0.0
        CS%diag_TKE_mixing(i,j) = 0.0 ; CS%diag_TKE_mech_decay(i,j) = 0.0
        CS%diag_TKE_conv_decay(i,j) = 0.0 !; CS%diag_TKE_unbalanced_forcing(i,j) = 0.0
      enddo ; enddo
    endif
    if (CS%Mixing_Diagnostics) then
      CS%Mixing_Length(:,:,:) = 0.0
      CS%Velocity_Scale(:,:,:) = 0.0
    endif
  endif

  if (present(dT_expected)) then
     dT_expected(:,:,:) = 0.0
  endif
  if (present(dS_expected)) then
     dS_expected(:,:,:) = 0.0
  endif

  do j=js,je
    ! Copy the thicknesses and other fields to 2-d arrays.
    do k=1,nz ; do i=is,ie
      h(i,k) = h_3d(i,j,k) + GV%H_subroundoff ; u(i,k) = u_3d(i,j,k) ; v(i,k) = v_3d(i,j,k)
      T(i,k) = tv%T(i,j,k) ; S(i,k) = tv%S(i,j,k)
      Kd(i,K) = 0.0
    enddo ; enddo
    do i=is,ie
      CS%ML_depth(i,j) = h(i,1)*GV%H_to_m
      !CS%ML_depth2(i,j) = h(i,1)*GV%H_to_m
!?      sfc_connected(i) = .true. !Niki Does this need to be initialized here?
    enddo

!?    if (debug) then
!?      mech_TKE_k(:,:) = 0.0 ; conv_PErel_k(:,:) = 0.0
!?    endif

    !   Determine the initial mech_TKE and conv_PErel, including the energy required
    ! to mix surface heating through the topmost cell, the energy released by mixing
    ! surface cooling & brine rejection down through the topmost cell, and
    ! homogenizing the shortwave heating within that cell.  This sets the energy
    ! and ustar and wstar available to drive mixing at the first interior
    ! interface.
    do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
      U_Star = fluxes%ustar(i,j)
      U_Star_Mean = fluxes%ustar_gustless(i,j)
      if (associated(fluxes%ustar_shelf) .and. associated(fluxes%frac_shelf_h)) then
        if (fluxes%frac_shelf_h(i,j) > 0.0) &
          U_Star = (1.0 - fluxes%frac_shelf_h(i,j)) * U_star + &
                    fluxes%frac_shelf_h(i,j) * fluxes%ustar_shelf(i,j)
      endif
      if (U_Star < CS%ustar_min) U_Star = CS%ustar_min
      if (CS%omega_frac >= 1.0) then ; absf(i) = 2.0*CS%omega
      else
        absf(i) = 0.25*((abs(G%CoriolisBu(I,J)) + abs(G%CoriolisBu(I-1,J-1))) + &
                     (abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I-1,J))))
        if (CS%omega_frac > 0.0) &
          absf(i) = sqrt(CS%omega_frac*4.0*CS%omega**2 + (1.0-CS%omega_frac)*absf(i)**2)
      endif
      if (CS%Mstar_Mode == CS%CONST_MSTAR) then
        mech_TKE(i) = (dt*CS%mstar*GV%Rho0)*((U_Star**3))
        conv_PErel(i) = 0.0

        if (CS%TKE_diagnostics) then
          CS%diag_TKE_wind(i,j) = CS%diag_TKE_wind(i,j) + mech_TKE(i) * IdtdR0
          if (TKE_forced(i,j,1) <= 0.0) then
            CS%diag_TKE_forcing(i,j) = CS%diag_TKE_forcing(i,j) + &
                 max(-mech_TKE(i), TKE_forced(i,j,1)) * IdtdR0
            ! CS%diag_TKE_unbalanced_forcing(i,j) = CS%diag_TKE_unbalanced_forcing(i,j) + &
            !     min(0.0, TKE_forced(i,j,1) + mech_TKE(i)) * IdtdR0
          else
            CS%diag_TKE_forcing(i,j) = CS%diag_TKE_forcing(i,j) + CS%nstar*TKE_forced(i,j,1) * IdtdR0
          endif
        endif

        if (TKE_forced(i,j,1) <= 0.0) then
          mech_TKE(i) = mech_TKE(i) + TKE_forced(i,j,1)
          if (mech_TKE(i) < 0.0) mech_TKE(i) = 0.0
        else
          conv_PErel(i) = conv_PErel(i) + TKE_forced(i,j,1)
        endif

      endif

      call cvmix_epbl_column(GV%ke, Kd(i,:), h(i,:),u(i,:),v(i,:),T(i,:),S(i,:), &
                       dSV_dT(i,j,:),dSV_dS(i,j,:),TKE_forced(i,j,:),dt, IdtdR0, &
                       GV%m_to_H, GV%H_to_m,GV%H_to_kg_m2,GV%g_Earth,GV%Rho0,GV%H_subroundoff,&
                       Buoy_Flux(i,j),U_Star,U_Star_Mean,absf(i),mech_TKE(i),conv_PErel(i),&
                       CS%ML_Depth(i,j),CS%ML_Depth2(i,j),&
                       i,j, & !?Niki: just to be able to set the CS% arrays. Can we get rid of them?
                       CS%cvmix_CS) !?Niki: container type for the subroutine
!                       dT_expected(i,j,:), dS_expected(i,j,:)) !?Niki: Cannot use these, cause crash

    else
      ! For masked points, Kd_int must still be set (to 0) because it has intent(out).
      do K=1,nz+1
        Kd(i,K) = 0.
      enddo
    endif ; enddo ; ! Close of i-loop - Note unusual loop order!

    if (CS%id_Hsfc_used > 0) then
      do i=is,ie ; Hsfc_used(i,j) = h(i,1)*GV%H_to_m ; enddo
      do k=2,nz ; do i=is,ie
        if (Kd(i,K) > 0.0) Hsfc_used(i,j) = Hsfc_used(i,j) + h(i,k)*GV%H_to_m
      enddo ; enddo
    endif

    do K=1,nz+1 ; do i=is,ie
      Kd_int(i,j,K) = Kd(i,K)
    enddo ; enddo

  enddo ! j-loop

  if (write_diags) then
    if (CS%id_ML_depth > 0) &
      call post_data(CS%id_ML_depth, CS%ML_depth, MY_diag)
    if (CS%id_TKE_wind > 0) &
      call post_data(CS%id_TKE_wind, CS%diag_TKE_wind, MY_diag)
    if (CS%id_TKE_MKE > 0) &
      call post_data(CS%id_TKE_MKE, CS%diag_TKE_MKE, MY_diag)
    if (CS%id_TKE_conv > 0) &
      call post_data(CS%id_TKE_conv, CS%diag_TKE_conv, MY_diag)
    if (CS%id_TKE_forcing > 0) &
      call post_data(CS%id_TKE_forcing, CS%diag_TKE_forcing, MY_diag)
    if (CS%id_TKE_mixing > 0) &
      call post_data(CS%id_TKE_mixing, CS%diag_TKE_mixing, MY_diag)
    if (CS%id_TKE_mech_decay > 0) &
      call post_data(CS%id_TKE_mech_decay, CS%diag_TKE_mech_decay, MY_diag)
    if (CS%id_TKE_conv_decay > 0) &
      call post_data(CS%id_TKE_conv_decay, CS%diag_TKE_conv_decay, MY_diag)
    if (CS%id_Hsfc_used > 0) &
      call post_data(CS%id_Hsfc_used, Hsfc_used, MY_diag)
    if (CS%id_Mixing_Length > 0) &
      call post_data(CS%id_Mixing_Length, CS%Mixing_Length, MY_diag)
    if (CS%id_Velocity_Scale >0) &
      call post_data(CS%id_Velocity_Scale, CS%Velocity_Scale, MY_diag)
    if (CS%id_OSBL >0) &
      call post_data(CS%id_OSBL, CS%ML_Depth2, MY_diag)
    if (CS%id_LT_Enhancement >0) &
      call post_data(CS%id_LT_Enhancement, CS%Enhance_M, MY_diag)
    if (CS%id_MSTAR_MIX >0) &
      call post_data(CS%id_MSTAR_MIX, CS%MSTAR_MIX, MY_diag)
    if (CS%id_MLD_OBUKHOV >0) &
      call post_data(CS%id_MLD_Obukhov, CS%MLD_OBUKHOV, MY_diag)
    if (CS%id_MLD_EKMAN >0) &
      call post_data(CS%id_MLD_Ekman, CS%MLD_EKMAN, MY_diag)
    if (CS%id_Ekman_Obukhov >0) &
      call post_data(CS%id_Ekman_Obukhov, CS%Ekman_Obukhov, MY_diag)
    if (CS%id_LA >0) &
      call post_data(CS%id_LA, CS%LA, MY_diag)
    if (CS%id_LA_MOD >0) &
      call post_data(CS%id_LA_MOD, CS%LA_MOD, MY_diag)

  endif
end subroutine energetic_PBL

!> Copies the ePBL active mixed layer depth into MLD
subroutine energetic_PBL_get_MLD(CS, MLD, G)
  type(energetic_PBL_CS),           pointer     :: CS  !< Control structure for ePBL
  type(ocean_grid_type),            intent(in)  :: G   !< Grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: MLD !< Depth of ePBL active mixing layer
  ! Local variables
  integer :: i,j
  do j = G%jsc, G%jec ; do i = G%isc, G%iec
    MLD(i,j) = CS%ML_depth(i,j)
  enddo ; enddo
end subroutine energetic_PBL_get_MLD


subroutine energetic_PBL_init(Time, G, GV, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target, intent(inout) :: diag
  type(energetic_PBL_CS), pointer        :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                  for this module
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_energetic_PBL"  ! This module's name.
  real :: omega_frac_dflt
  integer :: isd, ied, jsd, jed
  logical :: use_temperature, use_omega
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (associated(CS)) then
    call MOM_error(WARNING, "mixedlayer_init called with an associated"// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  MY_diag => diag

! Set default, read and log parameters
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "MSTAR_MODE", CS%mstar_mode, &
                 "An integer switch for how to compute MSTAR. \n"//&
                 "    0 for constant MSTAR\n"//&
                 "    1 for MSTAR w/ MLD in stabilizing limit\n"//&
                 "    2 for MSTAR w/ L_E/L_O in stabilizing limit.",&
                 "units=nondim",default=0)
  call get_param(param_file, mdl, "MSTAR", CS%mstar, &
                 "The ratio of the friction velocity cubed to the TKE \n"//&
                 "input to the mixed layer.", "units=nondim", default=1.2)
  call get_param(param_file, mdl, "MIX_LEN_EXPONENT", CS%MixLenExponent, &
                 "The exponent applied to the ratio of the distance to the MLD \n"//&
                 "and the MLD depth which determines the shape of the mixing length.",&
                 "units=nondim", default=2.0)
  call get_param(param_file, mdl, "MSTAR_CAP", CS%mstar_cap, &
                 "Maximum value of mstar allowed in model if non-negative\n"//&
                 "(used if MSTAR_MODE>0).",&
                 "units=nondim", default=-1.0)
  call get_param(param_file, mdl, "MSTAR_CONV_ADJ", CS%cnv_mst_fac, &
                 "Factor used for reducing mstar during convection \n"//&
                 " due to reduction of stable density gradient.",&
                 "units=nondim", default=0.0)
  call get_param(param_file, mdl, "MSTAR_SLOPE", CS%mstar_slope, &
                 "The slope of the linear relationship between mstar \n"//&
                 "and the length scale ratio (used if MSTAR_MODE=1).",&
                 "units=nondim", default=0.85)
  call get_param(param_file, mdl, "MSTAR_XINT", CS%mstar_xint, &
                 "The value of the length scale ratio where the mstar \n"//&
                 "is linear above (used if MSTAR_MODE=1).",&
                 "units=nondim", default=-0.3)
  call get_param(param_file, mdl, "MSTAR_AT_XINT", CS%mstar_at_xint, &
                 "The value of mstar at MSTAR_XINT \n"//&
                 "(used if MSTAR_MODE=1).",&
                 "units=nondim", default=0.095)
  call get_param(param_file, mdl, "MSTAR_FLATCAP", CS%MSTAR_FLATCAP, &
                 "Set false to use asymptotic cap, defaults to true.\n"//&
                 "(used only if MSTAR_MODE=1)"&
                 ,"units=nondim",default=.true.)
  call get_param(param_file, mdl, "MSTAR2_COEF1", CS%MSTAR_COEF, &
                 "Coefficient in computing mstar when rotation and \n"//&
                 " stabilizing effects are both important (used if MSTAR_MODE=2)"&
                  ,"units=nondim",default=0.3)
  call get_param(param_file, mdl, "MSTAR2_COEF2", CS%C_EK, &
                 "Coefficient in computing mstar when only rotation limits \n"//&
                 " the total mixing. (used only if MSTAR_MODE=2)"&
                  ,"units=nondim",default=0.085)
  call get_param(param_file, mdl, "NSTAR", CS%nstar, &
                 "The portion of the buoyant potential energy imparted by \n"//&
                 "surface fluxes that is available to drive entrainment \n"//&
                 "at the base of mixed layer when that energy is positive.", &
                 units="nondim", default=0.2)
  call get_param(param_file, mdl, "MKE_TO_TKE_EFFIC", CS%MKE_to_TKE_effic, &
                 "The efficiency with which mean kinetic energy released \n"//&
                 "by mechanically forced entrainment of the mixed layer \n"//&
                 "is converted to turbulent kinetic energy.", units="nondim", &
                 default=0.0)
  call get_param(param_file, mdl, "TKE_DECAY", CS%TKE_decay, &
                 "TKE_DECAY relates the vertical rate of decay of the \n"//&
                 "TKE available for mechanical entrainment to the natural \n"//&
                 "Ekman depth.", units="nondim", default=2.5)
!  call get_param(param_file, mdl, "HMIX_MIN", CS%Hmix_min, &
!                 "The minimum mixed layer depth if the mixed layer depth \n"//&
!                 "is determined dynamically.", units="m", default=0.0)

  call get_param(param_file, mdl, "OMEGA",CS%omega,              &
                 "The rotation rate of the earth.", units="s-1", &
                 default=7.2921e-5)
  call get_param(param_file, mdl, "ML_USE_OMEGA", use_omega,                  &
                 "If true, use the absolute rotation rate instead of the \n"//&
                 "vertical component of rotation when setting the decay \n"// &
                 "scale for turbulence.", default=.false., do_not_log=.true.)
  omega_frac_dflt = 0.0
  if (use_omega) then
    call MOM_error(WARNING, "ML_USE_OMEGA is depricated; use ML_OMEGA_FRAC=1.0 instead.")
    omega_frac_dflt = 1.0
  endif
  call get_param(param_file, mdl, "ML_OMEGA_FRAC", CS%omega_frac,              &
                 "When setting the decay scale for turbulence, use this \n"//  &
                 "fraction of the absolute rotation rate blended with the \n"//&
                 "local value of f, as sqrt((1-of)*f^2 + of*4*omega^2).",      &
                 units="nondim", default=omega_frac_dflt)
  call get_param(param_file, mdl, "WSTAR_USTAR_COEF", CS%wstar_ustar_coef,     &
                 "A ratio relating the efficiency with which convectively \n"//&
                 "released energy is converted to a turbulent velocity, \n"//  &
                 "relative to mechanically forced TKE. Making this larger \n"//&
                 "increases the BL diffusivity", units="nondim", default=1.0)
  call get_param(param_file, mdl, "VSTAR_SCALE_FACTOR", CS%vstar_scale_fac, &
                 "An overall nondimensional scaling factor for v*. \n"//    &
                 "Making this larger decreases the PBL diffusivity.",       &
                 units="nondim", default=1.0)
  call get_param(param_file, mdl, "EKMAN_SCALE_COEF", CS%Ekman_scale_coef,           &
                 "A nondimensional scaling factor controlling the inhibition \n"//   &
                 "of the diffusive length scale by rotation. Making this larger \n"//&
                 "decreases the PBL diffusivity.", units="nondim", default=1.0)
  call get_param(param_file, mdl, "USE_MLD_ITERATION", CS%USE_MLD_ITERATION,    &
                 "A logical that specifies whether or not to use the \n"//      &
                 "distance to the bottom of the actively turblent boundary \n"//&
                 "layer to help set the EPBL length scale.", default=.false.)
  call get_param(param_file, mdl, "ORIG_MLD_ITERATION", CS%ORIG_MLD_ITERATION,  &
                 "A logical that specifies whether or not to use the \n"//      &
                 "old method for determining MLD depth in iteration, which \n"//&
                 "is limited to resolution.", default=.true.)
  call get_param(param_file, mdl, "MLD_ITERATION_GUESS", CS%MLD_ITERATION_GUESS,       &
                 "A logical that specifies whether or not to use the \n"//             &
                 "previous timestep MLD as a first guess in the MLD iteration.\n"//    &
                 "The default is false to facilitate reproducibility.", default=.false.)
  call get_param(param_file, mdl, "EPBL_MLD_TOLERANCE", CS%MLD_tol,         &
                 "The tolerance for the iteratively determined mixed \n"//  &
                 "layer depth.  This is only used with USE_MLD_ITERATION.", &
                 units="meter", default=1.0)
  call get_param(param_file, mdl, "EPBL_MIN_MIX_LEN", CS%min_mix_len,    &
                 "The minimum mixing length scale that will be used \n"//&
                 "by ePBL.  The default (0) does not set a minimum.",    &
                 units="meter", default=0.0)
  call get_param(param_file, mdl, "EPBL_ORIGINAL_PE_CALC", CS%orig_PE_calc,         &
                 "If true, the ePBL code uses the original form of the \n"//        &
                 "potential energy change code.  Otherwise, the newer \n"//         &
                 "version that can work with successive increments to the \n"//     &
                 "diffusivity in upward or downward passes is used.", default=.true.)
  call get_param(param_file, mdl, "EPBL_TRANSITION_SCALE", CS%transLay_scale, &
                 "A scale for the mixing length in the transition layer \n"// &
                 "at the edge of the boundary layer as a fraction of the \n"//&
                 "boundary layer thickness.  The default is 0.1.", &
                 units="nondim", default=0.1)
  if ( CS%USE_MLD_ITERATION .and. abs(CS%transLay_scale-0.5).ge.0.5) then
    call MOM_error(FATAL, "If flag USE_MLD_ITERATION is true, then "//       &
                 "EPBL_TRANSITION should be greater than 0 and less than 1.")
  endif
  call get_param(param_file, mdl, "N2_DISSIPATION_POS", CS%N2_Dissipation_Scale_Pos, &
                 "A scale for the dissipation of TKE due to stratification \n"//     &
                 "in the boundary layer, applied when local stratification \n"//     &
                 "is positive.  The default is 0, but should probably be ~0.4.",     &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "N2_DISSIPATION_NEG", CS%N2_Dissipation_Scale_Neg,&
                 "A scale for the dissipation of TKE due to stratification \n"//    &
                 "in the boundary layer, applied when local stratification \n"//    &
                 "is negative.  The default is 0, but should probably be ~1.",      &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "USE_LA_LI2016", CS%USE_LA_Windsea,            &
                 "A logical to use the Li et al. 2016 (submitted) formula to \n"//&
                 " determine the Langmuir number.",                               &
                 units="nondim", default=.false.)
  call get_param(param_file, mdl, "LA_DEPTH_RATIO", CS%LaDepthRatio,                &
                 "The depth (normalized by BLD) to average Stokes drift over in \n"//&
                 " Lanmguir number calculation, where La = sqrt(ust/Stokes).",       &
                 units="nondim",default=0.04)
  call get_param(param_file, mdl, "LT_ENHANCE", CS%LT_ENHANCE_FORM,        &
                 "Integer for Langmuir number mode. \n"//                   &
                 " *Requires USE_LA_LI2016 to be set to True. \n"//         &
                 "Options: 0 - No Langmuir \n"//                            &
                 "         1 - Van Roekel et al. 2014/Li et al., 2016  \n"//&
                 "         2 - Multiplied w/ adjusted La. \n"//             &
                 "         3 - Added w/ adjusted La.",                      &
                 units="nondim", default=0)
  call get_param(param_file, mdl, "LT_ENHANCE_COEF", CS%LT_ENHANCE_COEF, &
                 "Coefficient for Langmuir enhancement if LT_ENHANCE > 1",&
                 units="nondim", default=0.447)
  call get_param(param_file, mdl, "LT_ENHANCE_EXP", CS%LT_ENHANCE_EXP, &
                 "Exponent for Langmuir enhancement if LT_ENHANCE > 1", &
                 units="nondim", default=-1.33)
  call get_param(param_file, mdl, "LT_MOD_LAC1", CS%LaC_MLDoEK,             &
                 "Coefficient for modification of Langmuir number due to\n"//&
                 " MLD approaching Ekman depth if LT_ENHANCE=2.",            &
                 units="nondim", default=-0.87)
  call get_param(param_file, mdl, "LT_MOD_LAC2", CS%LaC_MLDoOB_stab,        &
                 "Coefficient for modification of Langmuir number due to\n"//&
                 " MLD approaching stable Obukhov depth if LT_ENHANCE=2.",   &
                  units="nondim", default=0.0)
  call get_param(param_file, mdl, "LT_MOD_LAC3", CS%LaC_MLDoOB_un,          &
                 "Coefficient for modification of Langmuir number due to\n"//&
                 " MLD approaching unstable Obukhov depth if LT_ENHANCE=2.", &
                  units="nondim", default=0.0)
  call get_param(param_file, mdl, "LT_MOD_LAC4", CS%Lac_EKoOB_stab,         &
                 "Coefficient for modification of Langmuir number due to\n"//&
                 " ratio of Ekman to stable Obukhov depth if LT_ENHANCE=2.", &
                  units="nondim", default=0.95)
  call get_param(param_file, mdl, "LT_MOD_LAC5", CS%Lac_EKoOB_un,            &
                 "Coefficient for modification of Langmuir number due to\n"// &
                 " ratio of Ekman to unstable Obukhov depth if LT_ENHANCE=2.",&
                  units="nondim", default=0.95)
  if (CS%LT_ENHANCE_FORM>0 .and. (.not.CS%USE_LA_Windsea)) then
    call MOM_error(FATAL, "If flag USE_LA_LI2016 is false, LT_ENHANCE must be 0.")
  endif
  ! This gives a minimum decay scale that is typically much less than Angstrom.
  CS%ustar_min = 2e-4*CS%omega*(GV%Angstrom_z + GV%H_to_m*GV%H_subroundoff)
  call log_param(param_file, mdl, "EPBL_USTAR_MIN", CS%ustar_min, &
                 "The (tiny) minimum friction velocity used within the \n"//&
                 "ePBL code, derived from OMEGA and ANGSTROM.", units="m s-1")

  CS%id_ML_depth = register_diag_field('ocean_model', 'ePBL_h_ML', diag%axesT1, &
      Time, 'Surface boundary layer depth', 'm',                            &
      cmor_long_name='Ocean Mixed Layer Thickness Defined by Mixing Scheme')
  CS%id_TKE_wind = register_diag_field('ocean_model', 'ePBL_TKE_wind', diag%axesT1, &
      Time, 'Wind-stirring source of mixed layer TKE', 'm3 s-3')
  CS%id_TKE_MKE = register_diag_field('ocean_model', 'ePBL_TKE_MKE', diag%axesT1, &
      Time, 'Mean kinetic energy source of mixed layer TKE', 'm3 s-3')
  CS%id_TKE_conv = register_diag_field('ocean_model', 'ePBL_TKE_conv', diag%axesT1, &
      Time, 'Convective source of mixed layer TKE', 'm3 s-3')
  CS%id_TKE_forcing = register_diag_field('ocean_model', 'ePBL_TKE_forcing', diag%axesT1, &
      Time, 'TKE consumed by mixing surface forcing or penetrative shortwave radation'//&
            ' through model layers', 'm3 s-3')
  CS%id_TKE_mixing = register_diag_field('ocean_model', 'ePBL_TKE_mixing', diag%axesT1, &
      Time, 'TKE consumed by mixing that deepens the mixed layer', 'm3 s-3')
  CS%id_TKE_mech_decay = register_diag_field('ocean_model', 'ePBL_TKE_mech_decay', diag%axesT1, &
      Time, 'Mechanical energy decay sink of mixed layer TKE', 'm3 s-3')
  CS%id_TKE_conv_decay = register_diag_field('ocean_model', 'ePBL_TKE_conv_decay', diag%axesT1, &
      Time, 'Convective energy decay sink of mixed layer TKE', 'm3 s-3')
  CS%id_Hsfc_used = register_diag_field('ocean_model', 'ePBL_Hs_used', diag%axesT1, &
      Time, 'Surface region thickness that is used', 'm')
  CS%id_Mixing_Length = register_diag_field('ocean_model', 'Mixing_Length', diag%axesTi, &
      Time, 'Mixing Length that is used', 'm')
  CS%id_Velocity_Scale = register_diag_field('ocean_model', 'Velocity_Scale', diag%axesTi, &
      Time, 'Velocity Scale that is used.', 'm s-1')
  CS%id_LT_enhancement = register_diag_field('ocean_model', 'LT_Enhancement', diag%axesT1, &
      Time, 'LT enhancement that is used.', 'nondim')
  CS%id_MSTAR_mix = register_diag_field('ocean_model', 'MSTAR', diag%axesT1, &
      Time, 'MSTAR that is used.', 'nondim')
  CS%id_OSBL = register_diag_field('ocean_model', 'ePBL_OSBL', diag%axesT1, &
      Time, 'ePBL Surface Boundary layer depth.', 'm')
  ! BGR (9/21/2017) Note that ePBL_OSBL is the guess for iteration step while ePBL_h_ML is
  !                 result from iteration step.
  CS%id_mld_ekman = register_diag_field('ocean_model', 'MLD_EKMAN', diag%axesT1, &
      Time, 'Boundary layer depth over Ekman length.', 'm')
  CS%id_mld_obukhov = register_diag_field('ocean_model', 'MLD_OBUKHOV', diag%axesT1, &
      Time, 'Boundary layer depth over Obukhov length.', 'm')
  CS%id_ekman_obukhov = register_diag_field('ocean_model', 'EKMAN_OBUKHOV', diag%axesT1, &
      Time, 'Ekman length over Obukhov length.', 'm')
  CS%id_LA = register_diag_field('ocean_model', 'LA', diag%axesT1, &
      Time, 'Langmuir number.', 'nondim')
  CS%id_LA_mod = register_diag_field('ocean_model', 'LA_MOD', diag%axesT1, &
      Time, 'Modified Langmuir number.', 'nondim')


  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", use_temperature, &
                 "If true, temperature and salinity are used as state \n"//&
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
  if ((CS%id_Mixing_Length>0) .or. (CS%id_Velocity_Scale>0)) then
    call safe_alloc_alloc(CS%Velocity_Scale,isd,ied,jsd,jed,GV%ke+1)
    call safe_alloc_alloc(CS%Mixing_Length,isd,ied,jsd,jed,GV%ke+1)
    CS%Velocity_Scale(:,:,:) = 0.0
    CS%Mixing_Length(:,:,:) = 0.0
    CS%mixing_diagnostics = .true.
  endif
  call safe_alloc_alloc(CS%ML_depth, isd, ied, jsd, jed)
  call safe_alloc_alloc(CS%ML_depth2, isd, ied, jsd, jed)
  if (max(CS%id_LT_Enhancement, CS%id_mstar_mix,CS%id_mld_ekman, &
       CS%id_ekman_obukhov, CS%id_mld_obukhov, CS%id_LA, CS%id_LA_mod)>0) then
    call safe_alloc_alloc(CS%Mstar_mix, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%Enhance_M, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%MLD_EKMAN, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%MLD_OBUKHOV, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%EKMAN_OBUKHOV, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%LA, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%LA_MOD, isd, ied, jsd, jed)
  endif

  !Fitting coefficients to asymptote twoard 0 as MLD -> Ekman depth
  CS%MSTAR_A = CS%MSTAR_AT_XINT**(1./CS%MSTAR_N)
  CS%MSTAR_B = CS%MSTAR_SLOPE / (CS%MSTAR_N*CS%MSTAR_A**(CS%MSTAR_N-1.))
  !Fitting coefficients to asymptote toward MSTAR_CAP
  !*Fixed to begin asymptote at MSTAR_CAP-0.5 toward MSTAR_CAP
  CS%MSTAR_A2 = 0.5**(1./CS%MSTAR_N)
  CS%MSTAR_B2 = -CS%MSTAR_SLOPE / (CS%MSTAR_N*CS%MSTAR_A2**(CS%MSTAR_N-1))
  !Compute value of X (referenced to MSTAR_XINT) where transition
  ! to asymptotic regime based on value of X where MSTAR=MSTAR_CAP-0.5
  CS%MSTAR_XINT_UP = (CS%MSTAR_CAP-0.5-CS%MSTAR_AT_XINT)/CS%MSTAR_SLOPE

  call cvmix_epbl_init(CS%mstar,CS%nstar,CS%MixLenExponent,CS%TKE_decay,CS%MKE_to_TKE_effic,CS%omega,CS%omega_frac,&
CS%wstar_ustar_coef,CS%vstar_scale_fac,CS%transLay_scale,CS%MLD_tol,CS%min_mix_len,&
CS%N2_Dissipation_Scale_Neg,CS%N2_Dissipation_Scale_Pos,CS%MSTAR_CAP,CS%MSTAR_SLOPE,CS%MSTAR_XINT,CS%MSTAR_AT_XINT,&
CS%LT_ENHANCE_COEF,CS%LT_ENHANCE_EXP,CS%MSTAR_N,CS%C_EK,CS%MSTAR_COEF,CS%MSTAR_A,CS%MSTAR_B,CS%MSTAR_A2,CS%MSTAR_B2,&
CS%LaC_MLDoEK,CS%LaC_MLDoOB_stab,CS%LaC_EKoOB_stab,&
CS%LaC_MLDoOB_un,CS%LaC_EKoOB_un,CS%Max_Enhance_M,CS%CNV_MST_FAC,CS%LT_Enhance_Form,CS%MSTAR_MODE,&
CS%CONST_MSTAR,CS%MLD_o_OBUKHOV,CS%EKMAN_o_OBUKHOV,CS%MSTAR_FLATCAP,CS%TKE_diagnostics,CS%Use_LA_windsea,&
CS%orig_PE_calc,CS%Use_MLD_iteration,CS%Orig_MLD_iteration,CS%MLD_iteration_guess,CS%Mixing_Diagnostics,CS%cvmix_CS)

end subroutine energetic_PBL_init

subroutine energetic_PBL_end(CS)
  type(energetic_PBL_CS), pointer        :: CS

  call cvmix_epbl_end(CS%cvmix_CS)

  deallocate(CS)

end subroutine energetic_PBL_end

end module MOM_energetic_PBL
