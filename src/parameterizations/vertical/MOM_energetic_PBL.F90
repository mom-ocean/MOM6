module MOM_energetic_PBL
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

implicit none ; private

#include <MOM_memory.h>

public energetic_PBL, energetic_PBL_init, energetic_PBL_end
public energetic_PBL_get_MLD

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
  real    :: LT_ENHANCE_COEF ! Coefficient in fit for Langmuir Enhancment
  real    :: LT_ENHANCE_EXP  ! Exponent in fit for Langmuir Enhancement
  real :: MSTAR_N = -2.      ! Exponent in decay at negative and positive limits of MLD_over_STAB
  real :: MSTAR_A,MSTAR_A2   ! MSTAR_A and MSTAR_B are coefficients in asymptote toward limits.
  real :: MSTAR_B,MSTAR_B2   !  These are computed to match the function value and slope at both
                             !  ends of the linear fit within the well constrained region.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  integer :: LT_Enhance_Form = 0 ! Option for Langmuir enhancement function
  logical :: Use_Mstar_Fixed = .true. ! A logical to revert to a fixed m*
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
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.

! These are terms in the mixed layer TKE budget, all in J m-2 = kg s-2.
  real, allocatable, dimension(:,:) :: &
    ML_depth, &        ! The mixed layer depth in m.
    ML_depth2, &       ! The mixed layer depth in m.
    Enhance_V, &       ! The enhancement to the turbulent velocity scale (non-dim)
    MSTAR_MIX, &       ! Mstar used in EPBL
    diag_TKE_wind, &   ! The wind source of TKE.
    diag_TKE_MKE, &    ! The resolved KE source of TKE.
    diag_TKE_conv, &   ! The convective source of TKE.
    diag_TKE_forcing, & ! The TKE sink required to mix surface
                       ! penetrating shortwave heating.
    diag_TKE_mech_decay, & ! The decay of mechanical TKE.
    diag_TKE_conv_decay, & ! The decay of convective TKE.
    diag_TKE_mixing    ! The work done by TKE to deepen
                       ! the mixed layer.
  real, allocatable, dimension(:,:,:) :: &
    Velocity_Scale, & ! The velocity scale used in getting Kd
    Mixing_Length     ! The length scale used in getting Kd
  integer :: id_ML_depth = -1, id_TKE_wind = -1, id_TKE_mixing = -1
  integer :: id_TKE_MKE = -1, id_TKE_conv = -1, id_TKE_forcing = -1
  integer :: id_TKE_mech_decay = -1, id_TKE_conv_decay = -1
  integer :: id_Hsfc_used = -1
  integer :: id_Mixing_Length = -1, id_Velocity_Scale = -1
  integer :: id_OSBL = -1, id_LT_Enhancement = -1, id_MSTAR_mix = -1
end type energetic_PBL_CS

integer :: num_msg = 0, max_msg = 2

contains

subroutine energetic_PBL(h_3d, u_3d, v_3d, tv, fluxes, dt, Kd_int, G, GV, CS, &
                         dSV_dT, dSV_dS, TKE_forced, Buoy_Flux, dt_diag, last_call, &
                         dT_expected, dS_expected )
  type(ocean_grid_type),                     intent(inout) :: G
  type(verticalGrid_type),                   intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h_3d
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u_3d, v_3d, dSV_dT, dSV_dS
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: TKE_forced
  type(thermo_var_ptrs),                     intent(inout) :: tv
  type(forcing),                             intent(inout) :: fluxes
  real,                                      intent(in)    :: dt
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(out) :: Kd_int
  type(energetic_PBL_CS),                    pointer       :: CS
  real, dimension(SZI_(G),SZJ_(G)), intent(in)             :: Buoy_Flux
  real,                            optional, intent(in)    :: dt_diag
  logical,                         optional, intent(in)    :: last_call
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

  real, dimension(SZI_(G),SZK_(GV)) :: &
    h, &            !   The layer thickness, in H (usually m or kg m-2).
    T, &            !   The layer temperatures, in deg C.
    S, &            !   The layer salinities, in psu.
    u, &            !   The zonal velocity, in m s-1.
    v               !   The meridional velocity, in m s-1.
  real, dimension(SZI_(G),SZK_(GV)+1) :: &
    Kd, &           ! The diapycnal diffusivity, in m2 s-1.
    pres, &         ! Interface pressures in Pa.
    hb_hs           ! The distance from the bottom over the thickness of the
                    ! water column, nondim.
  real, dimension(SZI_(G)) :: &
    mech_TKE, &     !   The mechanically generated turbulent kinetic energy
                    ! available for mixing over a time step, in J m-2 = kg s-2.
    conv_PErel, &   ! The potential energy that has been convectively released
                    ! during this timestep, in J m-2 = kg s-2. A portion nstar_FC
                    ! of conv_PErel is available to drive mixing.
    htot, &         !   The total depth of the layers above an interface, in H.
    uhtot, &        !   The depth integrated zonal and meridional velocities in the
    vhtot, &        ! layers above, in H m s-1.
    mech_TKE_top, & !    The value of mech_TKE at the top of the column, in J m-2.
    conv_PErel_top, & !  The value of conv_PErel at the top of the column, in J m-2.

    Idecay_len_TKE, &  ! The inverse of a turbulence decay length scale, in H-1.
    h_sum, &        ! The total thickness of the water column, in H.
    absf            ! The absolute value of f, in s-1.


  real, dimension(SZI_(G),SZK_(GV)) :: &
    dT_to_dColHt, & ! Partial derivatives of the total column height with the temperature
    dS_to_dColHt, & ! and salinity changes within a layer, in m K-1 and m ppt-1.
    dT_to_dPE, &    ! Partial derivatives of column potential energy with the temperature
    dS_to_dPE, &    ! and salinity changes within a layer, in J m-2 K-1 and J m-2 ppt-1.
    dT_to_dColHt_a, & ! Partial derivatives of the total column height with the temperature
    dS_to_dColHt_a, & ! and salinity changes within a layer, including the implicit effects
                    ! of mixing with layers higher in the water colun, in m K-1 and m ppt-1.
    dT_to_dPE_a, &  ! Partial derivatives of column potential energy with the temperature
    dS_to_dPE_a     ! and salinity changes within a layer, including the implicit effects
                    ! of mixing with layers higher in the water column, in
                    ! units of J m-2 K-1 and J m-2 ppt-1.
  real, dimension(SZK_(GV)) :: &
    T0, S0, &       ! Initial values of T and S in the column, in K and ppt.
    Te, Se, &       ! Estimated final values of T and S in the column, in K and ppt.
    c1, &           ! c1 is used by the tridiagonal solver, ND.
    dTe, dSe        ! Running (1-way) estimates of temperature and salinity change.
  real, dimension(SZK_(GV)) :: &
    Th_a, &         ! An effective temperature times a thickness in the layer above,
                    ! including implicit mixing effects with other yet higher layers, in K H.
    Sh_a, &         ! An effective salinity times a thickness in the layer above,
                    ! including implicit mixing effects with other yet higher layers, in K H.
    Th_b, &         ! An effective temperature times a thickness in the layer below,
                    ! including implicit mixing effects with other yet lower layers, in K H.
    Sh_b            ! An effective salinity times a thickness in the layer below,
                    ! including implicit mixing effects with other yet lower layers, in K H.
  real, dimension(SZI_(G)) :: &
    hp_a            ! An effective pivot thickness of the layer including the effects
                    ! of coupling with layers above, in H.  This is the first term
                    ! in the denominator of b1 in a downward-oriented tridiagonal solver.
  real, dimension(SZK_(GV)+1) :: &
    MixLen_shape, & ! A nondimensional shape factor for the mixing length that
                    ! gives it an appropriate assymptotic value at the bottom of
                    ! the boundary layer.
    Kddt_h          ! The diapycnal diffusivity times a timestep divided by the
                    ! average thicknesses around a layer, in H (m or kg m-2).
  real :: b1        ! b1 is inverse of the pivot used by the tridiagonal solver, in H-1.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in H.
  real :: dMass     ! The mass per unit area within a layer, in kg m-2.
  real :: dPres     ! The hydrostatic pressure change across a layer, in Pa.
  real :: dMKE_max  ! The maximum amount of mean kinetic energy that could be
                    ! converted to turbulent kinetic energy if the velocity in
                    ! the layer below an interface were homogenized with all of
                    ! the water above the interface, in J m-2 = kg s-2.
  real :: MKE2_Hharm ! Twice the inverse of the harmonic mean of the thickness
                    ! of a layer and the thickness of the water above, used in
                    ! the MKE conversion equation, in H-1.

  real :: dt_h      ! The timestep divided by the averages of the thicknesses around
                    ! a layer, times a thickness conversion factor, in H s m-2.
  real :: h_bot     ! The distance from the bottom, in H.
  real :: h_rsum    ! The running sum of h from the top, in H.
  real :: I_hs      ! The inverse of h_sum, in H-1.
  real :: I_mld     ! The inverse of the current value of MLD, in H-1.
  real :: h_tt      ! The distance from the surface or up to the next interface
                    ! that did not exhibit turbulent mixing from this scheme plus
                    ! a surface mixing roughness length given by h_tt_min, in H.
  real :: h_tt_min  ! A surface roughness length, in H.

  real :: C1_3      ! = 1/3.
  real :: vonKar    ! The vonKarman constant.
  real :: I_dtrho   ! 1.0 / (dt * Rho0) in m3 kg-1 s-1.  This is
                    ! used convert TKE back into ustar^3.
  real :: U_star    ! The surface friction velocity, in m s-1.
  real :: U_Star_Mean ! The surface friction without gustiness in m s-1.
  real :: vstar     ! An in-situ turbulent velocity, in m s-1.
  real :: Enhance_V ! An enhancement factor for vstar, based here on Langmuir impact.
  real :: LA        ! The Langmuir number (non-dim)
  real :: hbs_here  ! The local minimum of hb_hs and MixLen_shape, times a
                    ! conversion factor from H to M, in m H-1.
  real :: nstar_FC  ! The fraction of conv_PErel that can be converted to mixing, nondim.
  real :: TKE_reduc ! The fraction by which TKE and other energy fields are
                    ! reduced to support mixing, nondim. between 0 and 1.
  real :: tot_TKE   ! The total TKE available to support mixing at interface K, in J m-2.
  real :: TKE_here  ! The total TKE at this point in the algorithm, in J m-2.
  real :: dT_km1_t2 ! A diffusivity-independent term related to the temperature
                    ! change in the layer above the interface, in K.
  real :: dS_km1_t2 ! A diffusivity-independent term related to the salinity
                    ! change in the layer above the interface, in ppt.
  real :: dTe_term  ! A diffusivity-independent term related to the temperature
                    ! change in the layer below the interface, in K H.
  real :: dSe_term  ! A diffusivity-independent term related to the salinity
                    ! change in the layer above the interface, in ppt H.
  real :: dTe_t2    ! A part of dTe_term, in K H.
  real :: dSe_t2    ! A part of dSe_term, in ppt H.
  real :: dPE_conv  ! The convective change in column potential energy, in J m-2.
  real :: MKE_src   ! The mean kinetic energy source of TKE due to Kddt_h(K), in J m-2.
  real :: dMKE_src_dK  ! The partial derivative of MKE_src with Kddt_h(K), in J m-2 H-1.
  real :: Kd_guess0, PE_chg_g0, dPEa_dKd_g0, Kddt_h_g0
  real :: PE_chg_max   ! The maximum PE change for very large values of Kddt_h(K).
  real :: dPEc_dKd_Kd0 ! The partial derivative of PE change with Kddt_h(K)
                       ! for very small values of Kddt_h(K), in J m-2 H-1.
  real :: PE_chg    ! The change in potential energy due to mixing at an
                    ! interface, in J m-2, positive for the column increasing
                    ! in potential energy (i.e., consuming TKE).
  real :: TKE_left  ! The amount of turbulent kinetic energy left for the most
                    ! recent guess at Kddt_h(K), in J m-2.
  real :: dPEc_dKd  ! The partial derivative of PE_chg with Kddt_h(K), in J m-2 H-1.
  real :: TKE_left_min, TKE_left_max, Kddt_h_max, Kddt_h_min
  real :: Kddt_h_guess ! A guess at the value of Kddt_h(K), in H.
  real :: Kddt_h_next  ! The next guess at the value of Kddt_h(K), in H.
  real :: dKddt_h      ! The change between guesses at Kddt_h(K), in H.
  real :: dKddt_h_Newt ! The change between guesses at Kddt_h(K) with Newton's method, in H.
  real :: Kddt_h_newt  ! The Newton's method next guess for Kddt_h(K), in H.
  real :: exp_kh    ! The nondimensional decay of TKE across a layer, ND.
  logical :: use_Newt  ! Use Newton's method for the next guess at Kddt_h(K).
  logical :: convectively_stable
  logical, dimension(SZI_(G)) :: &
    sfc_connected   ! If true the ocean is actively turbulent from the present
                    ! interface all the way up to the surface.
  logical :: sfc_disconnect ! If true, any turbulence has become disconnected
                    ! from the surface.

! The following is only used as a diagnostic.
  real :: dt__diag  ! A copy of dt_diag (if present) or dt, in s.
  real :: IdtdR0    !  = 1.0 / (dt__diag * Rho0), in m3 kg-1 s-1.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    Hsfc_used   ! The thickness of the surface region after buffer layer
  logical :: write_diags  ! If true, write out diagnostics with this step.
  logical :: reset_diags  ! If true, zero out the accumulated diagnostics.
                ! detrainment, in units of m.
  ! Local column copies of energy change diagnostics, all in J m-2.
  real :: dTKE_conv, dTKE_forcing, dTKE_mixing
  real :: dTKE_MKE, dTKE_mech_decay, dTKE_conv_decay
  !----------------------------------------------------------------------
  !/BGR added Aug24,2016 for adding iteration to get boundary layer depth
  !    - needed to compute new mixing length.
  real :: MLD_guess, MLD_found ! Mixing Layer depth guessed/found for iteration, in m.
  real :: max_MLD, min_MLD ! Iteration bounds, in m, which are adjusted at each step
                           !  - These are initialized based on surface/bottom
                           !  1. The iteration guesses a value (possibly from
                           !     prev step or neighbor).
                           !  2. The iteration checks if value is converged,
                           !     too shallow, or too deep.
                           !  3. Based on result adjusts the Max/Min
                           !     and searches through the water column.
                           !  - If using an accurate guess the iteration
                           !    is very quick (e.g. if MLD doesn't change
                           !    over timestep).  Otherwise it takes 5-10
                           !    passes, but has a high convergence rate.
                           !    Other iteration may be tried, but this
                           !    method seems to rarely fail and the added
                           !    cost is likely not significant.  Additionally,
                           !    when it fails it does so in a reasonable
                           !    manner giving a usable guess. When it
                           !    does fail, it is due to convection within
                           !    the boundary.  Likely, a new method e.g.
                           !    surface_disconnect, can improve this.
  logical :: FIRST_OBL     ! Flag for computing "found" Mixing layer depth
  logical :: OBL_CONVERGED ! Flag for convergence of MLD
  integer :: OBL_IT        ! Iteration counter
!### These need to be made into run-time parameters.
  integer :: MAX_OBL_IT=20 ! Set maximum number of iterations.  Probably
                           !  best as an input parameter, but then may want
                           !  to use allocatable arrays if storing
                           !  guess/found (as diagnostic); skipping for now.
                           !  In reality, the maximum number of guesses
                           !  needed is set by:
                           !    DEPTH/2^M < DZ
                           !    where M is the number of guesses
                           !    e.g. M=12 for DEPTH=4000m and DZ=1m
  real, dimension(SZK_(GV)+1) :: Vstar_Used, &      ! 1D arrays used to store
                               Mixing_Length_Used   ! Vstar and Mixing_Length
  !/BGR - remaining variables are related to tracking iteration statistics.
  logical :: OBL_IT_STATS=.false. ! Flag for computing OBL iteration statistics
  REAL :: ITguess(20), ITresult(20),ITmax(20),ITmin(20) ! Flag for storing guess/result
                                                        ! should have dim=MAX_OBL_IT
  integer, save :: MAXIT=0   ! Stores maximum number of iterations
  integer, save :: MINIT=1e8 ! Stores minimum number of iterations
  integer, save :: SUMIT=0   ! Stores total iterations (summed over all)
  integer, save :: NUMIT=0   ! Stores number of times iterated
                             !e.g. Average iterations = SUMIT/NUMIT
  integer, save :: CONVERGED!
  integer, save :: NOTCONVERGED!
  !-End BGR iteration parameters-----------------------------------------
  real :: N2_dissipation
  real :: B_STAR ! Buoyancy flux over ustar
  real :: STAB_SCALE ! Composite of Stabilizing length scales:
                     !  Ekman scale and Monin-Obukhov scale.
  real :: C_MO = 1. ! Constant in STAB_SCALE for Monin-Obukhov
  real :: C_EK = 2. ! Constant in STAB_SCALE for Ekman length
  real :: MLD_over_STAB ! Mixing layer depth divided by STAB_SCALE
  real :: MSTAR_MIX! The value of mstar (Proportionality of TKE to drive mixing to ustar
                    ! cubed) which is computed as a function of latitude, boundary layer depth,
                    ! and the Monin-Obhukov depth.
  logical :: debug=.false.  ! Change this hard-coded value for debugging.
!  The following arrays are used only for debugging purposes.
  real :: dPE_debug, mixing_debug, taux2, tauy2
  real, dimension(20) :: TKE_left_itt, PE_chg_itt, Kddt_h_itt, dPEa_dKd_itt, MKE_src_itt
  real, dimension(SZI_(G),SZK_(GV)) :: &
    mech_TKE_k, conv_PErel_k
  real, dimension(SZK_(GV)) :: nstar_k
  integer, dimension(SZK_(GV)) :: num_itts

  integer :: i, j, k, is, ie, js, je, nz, itt, max_itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "energetic_PBL: "//&
         "Module must be initialized before it is used.")

  if (.not. ASSOCIATED(tv%eqn_of_state)) call MOM_error(FATAL, &
      "energetic_PBL: Temperature, salinity and an equation of state "//&
      "must now be used.")
  if (.NOT. ASSOCIATED(fluxes%ustar)) call MOM_error(FATAL, &
      "energetic_PBL: No surface TKE fluxes (ustar) defined in mixedlayer!")
  if (present(dT_expected) .or. present(dS_expected)) debug = .true.

  h_neglect = GV%H_subroundoff

  if(.not.CS%Use_MLD_Iteration) MAX_OBL_IT=1
  C1_3 = 1.0 / 3.0
  dt__diag = dt ; if (present(dt_diag)) dt__diag = dt_diag
  IdtdR0 = 1.0 / (dt__diag * GV%Rho0)
  write_diags = .true. ; if (present(last_call)) write_diags = last_call
  max_itt = 20

  h_tt_min = 0.0
  vonKar = 0.41
  mstar_mix=CS%MSTAR!Initialize to mstar
  I_dtrho = 0.0 ; if (dt*GV%Rho0 > 0.0) I_dtrho = 1.0 / (dt*GV%Rho0)

  ! Determine whether to zero out diagnostics before accumulation.
  reset_diags = .true.
  if (present(dt_diag) .and. write_diags .and. (dt__diag > dt)) &
    reset_diags = .false.  ! This is the second call to mixedlayer.

  if (reset_diags) then
!!OMP parallel default(none) shared(is,ie,js,je,CS)
    if (CS%TKE_diagnostics) then
!!OMP do
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
!!OMP end parallel
  endif


!!OMP parallel do default(none) shared(js,je,nz,is,ie,h_3d,u_3d,v_3d,tv,dt,      &
!!OMP                                  CS,G,GV,fluxes,IdtdR0,                    &
!!OMP                                  TKE_forced,debug,H_neglect,dSV_dT,        &
!!OMP                                  dSV_dS,I_dtrho,C1_3,h_tt_min,vonKar,     &
!!OMP                                  max_itt,Kd_int)                           &
!!OMP                           private(i,j,k,h,u,v,T,S,Kd,mech_TKE_k,conv_PErel_k,    &
!!OMP                                   U_Star,absf,mech_TKE,conv_PErel,nstar_k, &
!!OMP                                   h_sum,I_hs,h_bot,hb_hs,T0,S0,num_itts,   &
!!OMP                                   pres,dMass,dPres,dT_to_dPE,dS_to_dPE,    &
!!OMP                                   dT_to_dColHt,dS_to_dColHt,Kddt_h,hp_a,   &
!!OMP                                   Th_a,Sh_a,Th_b,Sh_b,dT_to_dPE_a,htot,    &
!!OMP                                   dT_to_dColHt_a,dS_to_dColHt_a,uhtot,vhtot, &
!!OMP                                   Idecay_len_TKE,exp_kh,nstar_FC,tot_TKE,  &
!!OMP                                   TKE_reduc,dTe_t2,dSe_t2,dTe,dSe,dt_h,    &
!!OMP                                   Convectively_stable,sfc_disconnect,b1,   &
!!OMP                                   c1,dT_km1_t2,dS_km1_t2,dTe_term,         &
!!OMP                                   dSe_term,MKE2_Hharm,vstar,h_tt,h_rsum,   &
!!OMP                                   Kd_guess0,Kddt_h_g0,dPEc_dKd_Kd0,        &
!!OMP                                   PE_chg_max,dPEa_dKd_g0,PE_chg_g0,        &
!!OMP                                   MKE_src,dPE_conv,Kddt_h_max,Kddt_h_min,  &
!!OMP                                   dTKE_conv, dTKE_forcing, dTKE_mixing,    &
!!OMP                                   dTKE_MKE,dTKE_mech_decay,dTKE_conv_decay,&
!!OMP                                   TKE_left_max,TKE_left_min,Kddt_h_guess,  &
!!OMP                                   TKE_left_itt,dPEa_dKd_itt,PE_chg_itt,    &
!!OMP                                   MKE_src_itt,Kddt_h_itt,dPEc_dKd,PE_chg,  &
!!OMP                                   dMKE_src_dK,TKE_left,use_Newt,           &
!!OMP                                   dKddt_h_Newt,Kddt_h_Newt,Kddt_h_next,    &
!!OMP                                   dKddt_h,Te,Se,Hsfc_used,dS_to_dPE_a,     &
!!OMP                                   dMKE_max,sfc_connected,TKE_here)
  do j=js,je
    ! Copy the thicknesses and other fields to 2-d arrays.
    do k=1,nz ; do i=is,ie
      h(i,k) = h_3d(i,j,k) + h_neglect ; u(i,k) = u_3d(i,j,k) ; v(i,k) = v_3d(i,j,k)
      T(i,k) = tv%T(i,j,k) ; S(i,k) = tv%S(i,j,k)
      Kd(i,K) = 0.0
    enddo ; enddo
    do i=is,ie
       CS%ML_depth(i,j) = h(i,1)*GV%H_to_m
       !CS%ML_depth2(i,j) = h(i,1)*GV%H_to_m
       sfc_connected(i) = .true.
    enddo

    if (debug) then
      mech_TKE_k(:,:) = 0.0 ; conv_PErel_k(:,:) = 0.0
    endif

    !   Determine the initial mech_TKE and conv_PErel, including the energy required
    ! to mix surface heating through the topmost cell, the energy released by mixing
    ! surface cooling & brine rejection down through the topmost cell, and
    ! homogenizing the shortwave heating within that cell.  This sets the energy
    ! and ustar and wstar available to drive mixing at the first interior
    ! interface.
    do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then


      U_Star = fluxes%ustar(i,j)
      taux2 = 0.
      if ((G%mask2dCu(I-1,j) + G%mask2dCu(I,j)) > 0) &
        taux2 = (G%mask2dCu(I-1,j)*fluxes%taux(I-1,j)**2 + &
                 G%mask2dCu(I,j)*fluxes%taux(I,j)**2) / (G%mask2dCu(I-1,j) + G%mask2dCu(I,j))
      tauy2 = 0.0
      if ((G%mask2dCv(i,J-1) + G%mask2dCv(i,J)) > 0) &
        tauy2 = (G%mask2dCv(i,J-1)*fluxes%tauy(i,J-1)**2 + &
                 G%mask2dCv(i,J)*fluxes%tauy(i,J)**2) / (G%mask2dCv(i,J-1) + G%mask2dCv(i,J))
      U_Star_Mean = sqrt(sqrt(taux2 + tauy2)/GV%rho0)
      if (associated(fluxes%ustar_shelf) .and. associated(fluxes%frac_shelf_h)) then
        if (fluxes%frac_shelf_h(i,j) > 0.0) &
          U_Star = (1.0 - fluxes%frac_shelf_h(i,j)) * U_star + &
                    fluxes%frac_shelf_h(i,j) * fluxes%ustar_shelf(i,j)
      endif
      if (U_Star < CS%ustar_min) U_Star = CS%ustar_min
      ! Computing B_Star, or the Buoyancy flux over the friction velocity.
      B_Star = min(0.0,-buoy_Flux(i,j)/U_Star)
      if (CS%omega_frac >= 1.0) then ; absf(i) = 2.0*CS%omega
      else
        absf(i) = 0.25*((abs(G%CoriolisBu(I,J)) + abs(G%CoriolisBu(I-1,J-1))) + &
                     (abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I-1,J))))
        if (CS%omega_frac > 0.0) &
          absf(i) = sqrt(CS%omega_frac*4.0*CS%omega**2 + (1.0-CS%omega_frac)*absf(i)**2)
      endif
      ! Computing stability scale which correlates with TKE for mixing, where
      ! TKE for mixing = TKE production minus TKE dissipation
      Stab_Scale = -u_star**2 / ( VonKar * ( C_MO * B_STAR +  C_EK * u_star * absf(i)))

      if (CS%Use_Mstar_Fixed) then
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

!    endif ; enddo

!    do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then

      h_sum(i) = H_neglect ; do k=1,nz ; h_sum(i) = h_sum(i) + h(i,k) ; enddo
      I_hs = 0.0 ; if (h_sum(i) > 0.0) I_hs = 1.0 / h_sum(i)

      h_bot = 0.0 ; hb_hs(i,nz+1) = 0.0
      do k=nz,1,-1
        h_bot = h_bot + h(i,k)
        hb_hs(i,K) = h_bot * I_hs
      enddo

      pres(i,1) = 0.0
      do k=1,nz
        dMass = GV%H_to_kg_m2 * h(i,k)
        dPres = GV%g_Earth * dMass
        dT_to_dPE(i,k) = (dMass * (pres(i,K) + 0.5*dPres)) * dSV_dT(i,j,k)
        dS_to_dPE(i,k) = (dMass * (pres(i,K) + 0.5*dPres)) * dSV_dS(i,j,k)
        dT_to_dColHt(i,k) = dMass * dSV_dT(i,j,k)
        dS_to_dColHt(i,k) = dMass * dSV_dS(i,j,k)

        pres(i,K+1) = pres(i,K) + dPres
      enddo

!    endif ; enddo

    ! Note the outer i-loop and inner k-loop loop order!!!
!    do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
      do k=1,nz ; T0(k) = T(i,k) ; S0(k) = S(i,k) ; enddo

      ! Store the initial mechanical TKE and convectively released PE to
      ! enable multiple iterations.
      mech_TKE_top(i) = mech_TKE(i) ; conv_PErel_top(i) = conv_PErel(i)

      !/The following lines are for the iteration over MLD
      !{
      ! max_MLD will initialized as ocean bottom depth
      max_MLD = 0.0 ; do k=1,nz ; max_MLD = max_MLD + h(i,k)*GV%H_to_m ; enddo
      min_MLD = 0.0 !min_MLD will initialize as 0.
      !/BGR: May add user-input bounds for max/min MLD

      !/BGR: Add MLD_guess based on stored previous value.
      !      note that this is different from ML_Depth already
      !      computed by EPBL, need to figure out why.
      if (CS%MLD_iteration_guess .and. CS%ML_Depth2(i,j) > 1.) then
        !If prev value is present use for guess.
        MLD_guess=CS%ML_Depth2(i,j)
      else
        !Otherwise guess middle of water column (or stab_scale if smaller).
        MLD_guess = 0.5 * (min_MLD+max_MLD)
      endif

      ! Iterate up to MAX_OBL_IT times to determine a converged EPBL depth.
      OBL_CONVERGED = .false.
      ! Initialize ENHANCE_V to 1
      ENHANCE_V=1.e0
      do OBL_IT=1,MAX_OBL_IT ; if (.not. OBL_CONVERGED) then
        if (CS%Use_LA_windsea) then
          call get_LA_windsea( u_star_mean, MLD_guess, GV, LA)
          if (CS%LT_Enhance_Form==1) then
            !Original w'/ust scaling
            ENHANCE_V = (1+(1.4*LA)**(-2)+(5.4*LA)**(-4))**(1.5)
          elseif (CS%LT_Enhance_Form==2) then
            !New Mix_LT/Mix_ST scaling
            ENHANCE_V = (1.+CS%LT_ENHANCE_COEF*LA**CS%LT_ENHANCE_EXP)
          endif
        endif
        CS%ML_depth(i,j) = h(i,1)*GV%H_to_m
        !CS%ML_depth2(i,j) = h(i,1)*GV%H_to_m

        sfc_connected(i) = .true.

        if (.not.CS%Use_Mstar_Fixed) then
          ! Note the value of mech_TKE(i) now must be iterated over, so it is moved here
          ! First solve for the TKE to PE length scale
          MLD_over_Stab = MLD_guess / Stab_Scale - CS%MSTAR_XINT
          if ((MLD_over_Stab) .le. 0.0) then
            !Asymptote to 0 as MLD_over_Stab -> -infinity (always)
            MSTAR_mix = (CS%MSTAR_B*(MLD_over_Stab)+CS%MSTAR_A)**(CS%MSTAR_N)
          else
            if (CS%MSTAR_CAP>=0.) then
              if (CS%MSTAR_FLATCAP .OR. (MLD_over_Stab .le.CS%MSTAR_XINT_UP)) then
                !If using flat cap (or if using asymptotic cap but within linear regime we
                ! can make use of same code)
                MSTAR_mix = min(CS%MSTAR_CAP, &
                     CS%MSTAR_SLOPE*(MLD_over_Stab)+CS%MSTAR_AT_XINT)
              else
                !Asymptote to MSTAR_CAP as MLD_over_Stab -> infinity
                MSTAR_mix = CS%MSTAR_CAP - &
                     (CS%MSTAR_B2*(MLD_over_Stab-CS%MSTAR_XINT_UP)+CS%MSTAR_A2)**(CS%MSTAR_N)
              endif
            else
              !No cap if negative cap value given.
              MSTAR_mix = CS%MSTAR_SLOPE*(MLD_over_Stab)+CS%MSTAR_AT_XINT
            endif
          endif

          !Reset mech_tke and conv_perel values (based on new mstar)
          mech_TKE(i) = (dt*MSTAR_mix*ENHANCE_V*GV%Rho0)*((U_Star**3))
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
        else
          mech_TKE(i) = mech_TKE_top(i)*ENHANCE_V ; conv_PErel(i) = conv_PErel_top(i)
        endif

        if (CS%TKE_diagnostics) then
          dTKE_conv = 0.0 ; dTKE_forcing = 0.0 ; dTKE_mixing = 0.0
          dTKE_MKE = 0.0 ; dTKE_mech_decay = 0.0 ; dTKE_conv_decay = 0.0
        endif

        ! Store in 1D arrays cleared out each iteration.  Only write in
        !  3D arrays after convergence.
        do k=1,nz
          Vstar_Used(k) = 0.0 ; Mixing_Length_Used(k) = 0.0
        enddo
        if (.not.CS%Use_MLD_Iteration) OBL_CONVERGED = .true.

        if ((.not.CS%Use_MLD_Iteration) .or. &
            (CS%transLay_scale >= 1.0) .or. (CS%transLay_scale < 0.0) ) then
          do K=1,nz+1 ; MixLen_shape(K) = 1.0 ; enddo
        elseif (MLD_guess <= 0.0) then
          if (CS%transLay_scale > 0.0) then
            do K=1,nz+1 ; MixLen_shape(K) = CS%transLay_scale ; enddo
          else
            do K=1,nz+1 ; MixLen_shape(K) = 1.0 ; enddo
          endif
        else
          ! Reduce the mixing length based on MLD, with a quadratic
          ! expression that follows KPP.
          I_MLD = 1.0 / MLD_guess ; h_rsum = 0.0
          MixLen_shape(1) = 1.0
          do K=2,nz+1
            h_rsum = h_rsum + h(i,k-1)
            if (CS%MixLenExponent==2.0)then
              MixLen_shape(K) = CS%transLay_scale + (1.0 - CS%transLay_scale) * &
                   (max(0.0, (MLD_guess - h_rsum)*I_MLD) )**2!CS%MixLenExponent
            else
              MixLen_shape(K) = CS%transLay_scale + (1.0 - CS%transLay_scale) * &
                   (max(0.0, (MLD_guess - h_rsum)*I_MLD) )**CS%MixLenExponent
            endif
          enddo
        endif


        Kd(i,1) = 0.0 ; Kddt_h(1) = 0.0
        hp_a(i) = h(i,1)
        dT_to_dPE_a(i,1) = dT_to_dPE(i,1) ; dT_to_dColHt_a(i,1) = dT_to_dColHt(i,1)
        dS_to_dPE_a(i,1) = dS_to_dPE(i,1) ; dS_to_dColHt_a(i,1) = dS_to_dColHt(i,1)

        htot(i) = h(i,1) ; uhtot(i) = u(i,1)*h(i,1) ; vhtot(i) = v(i,1)*h(i,1)

        if (debug) then
          mech_TKE_k(i,1) = mech_TKE(i) ; conv_PErel_k(i,1) = conv_PErel(i)
          nstar_k(:) = 0.0 ; nstar_k(1) = CS%nstar ; num_itts(:) = -1
        endif

        do K=2,nz
          ! Apply dissipation to the TKE, here applied as an exponential decay
          ! due to 3-d turbulent energy being lost to inefficient rotational modes.

          !   There should be several different "flavors" of TKE that decay at
          ! different rates.  The following form is often used for mechanical
          ! stirring from the surface, perhaps due to breaking surface gravity
          ! waves and wind-driven turbulence.
          Idecay_len_TKE(i) = (CS%TKE_decay * absf(i) / U_Star) * GV%H_to_m
          exp_kh = 1.0
          if (Idecay_len_TKE(i) > 0.0) exp_kh = exp(-h(i,k-1)*Idecay_len_TKE(i))
          if (CS%TKE_diagnostics) &
            dTKE_mech_decay = dTKE_mech_decay + (exp_kh-1.0) * mech_TKE(i) * IdtdR0
          mech_TKE(i) = mech_TKE(i) * exp_kh

          !   Accumulate any convectively released potential energy to contribute
          ! to wstar and to drive penetrating convection.
          if (TKE_forced(i,j,k) > 0.0) then
            conv_PErel(i) = conv_PErel(i) + TKE_forced(i,j,k)
            if (CS%TKE_diagnostics) &
              dTKE_forcing = dTKE_forcing + CS%nstar*TKE_forced(i,j,k) * IdtdR0
          endif

          if (debug) then
            mech_TKE_k(i,K) = mech_TKE(i) ; conv_PErel_k(i,K) = conv_PErel(i)
          endif

          !  Determine the total energy
          nstar_FC = CS%nstar
          if (CS%nstar * conv_PErel(i) > 0.0) then
            ! Here nstar is a function of the natural Rossby number 0.2/(1+0.2/Ro), based
            ! on a curve fit from the data of Wang (GRL, 2003).
            ! Note:         Ro = 1.0 / sqrt(0.5 * dt * Rho0 * (absf*htot(i))**3 / conv_PErel(i))
            nstar_FC = CS%nstar * conv_PErel(i) / (conv_PErel(i) + 0.2 * &
                            sqrt(0.5 * dt * GV%Rho0 * (absf(i)*(htot(i)*GV%H_to_m))**3 * conv_PErel(i)))
          endif
          if (debug) nstar_k(K) = nstar_FC

          tot_TKE = mech_TKE(i) + nstar_FC * conv_PErel(i)

          !   For each interior interface, first discard the TKE to account for
          ! mixing of shortwave radiation through the next denser cell.
          if (TKE_forced(i,j,k) < 0.0) then
            if (TKE_forced(i,j,k) + tot_TKE < 0.0) then
              ! The shortwave requirements deplete all the energy in this layer.
              if (CS%TKE_diagnostics) then
                dTKE_mixing = dTKE_mixing + tot_TKE * IdtdR0
                dTKE_forcing = dTKE_forcing - tot_TKE * IdtdR0
                ! dTKE_unbalanced_forcing = dTKE_unbalanced_forcing + &
                !     (TKE_forced(i,j,k) + tot_TKE) * IdtdR0
                dTKE_conv_decay = dTKE_conv_decay + &
                        (CS%nstar-nstar_FC) * conv_PErel(i) * IdtdR0
              endif
              tot_TKE = 0.0 ; mech_TKE(i) = 0.0 ; conv_PErel(i) = 0.0
            else
              ! Reduce the mechanical and convective TKE proportionately.
              TKE_reduc = (tot_TKE + TKE_forced(i,j,k)) / tot_TKE
              if (CS%TKE_diagnostics) then
                dTKE_mixing = dTKE_mixing - TKE_forced(i,j,k) * IdtdR0
                dTKE_forcing = dTKE_forcing + TKE_forced(i,j,k) * IdtdR0
                dTKE_conv_decay = dTKE_conv_decay + &
                    (1.0-TKE_reduc)*(CS%nstar-nstar_FC) * conv_PErel(i) * IdtdR0
              endif
              tot_TKE = TKE_reduc*tot_TKE   ! = tot_TKE + TKE_forced(i,j,k)
              mech_TKE(i) = TKE_reduc*mech_TKE(i)
              conv_PErel(i) = TKE_reduc*conv_PErel(i)
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
          dt_h = (GV%m_to_H**2*dt) / max(0.5*(h(i,k-1)+h(i,k)), 1e-15*h_sum(i))

          !   This tests whether the layers above and below this interface are in
          ! a convetively stable configuration, without considering any effects of
          ! mixing at higher interfaces.  It is an approximation to the more
          ! complete test dPEc_dKd_Kd0 >= 0.0, that would include the effects of
          ! mixing across interface K-1.  The dT_to_dColHt here are effectively
          ! mass-weigted estimates of dSV_dT.
          Convectively_stable = ( 0.0 <= &
            ( (dT_to_dColHt(i,k) + dT_to_dColHt(i,k-1) ) * (T0(k-1)-T0(k)) + &
              (dS_to_dColHt(i,k) + dS_to_dColHt(i,k-1) ) * (S0(k-1)-S0(k)) ) )

          if ((mech_TKE(i) + conv_PErel(i)) <= 0.0 .and. Convectively_stable) then
            ! Energy is already exhausted, so set Kd = 0 and cycle or exit?
            tot_TKE = 0.0 ; mech_TKE(i) = 0.0 ; conv_PErel(i) = 0.0
            Kd(i,K) = 0.0 ; Kddt_h(K) = 0.0
            sfc_disconnect = .true.
            ! if (.not.debug) exit

           !   The estimated properties for layer k-1 can be calculated, using
           ! greatly simplified expressions when Kddt_h = 0.  This enables the
           ! tridiagonal solver for the whole column to be completed for debugging
           ! purposes, and also allows for something akin to convective adjustment
           ! in unstable interior regions?
            b1 = 1.0 / hp_a(i)
            c1(K) = 0.0
            if (CS%orig_PE_calc) then
              dTe(k-1) = b1 * ( dTe_t2 )
              dSe(k-1) = b1 * ( dSe_t2 )
            endif

            hp_a(i) = h(i,k)
            dT_to_dPE_a(i,k) = dT_to_dPE(i,k)
            dS_to_dPE_a(i,k) = dS_to_dPE(i,k)
            dT_to_dColHt_a(i,k) = dT_to_dColHt(i,k)
            dS_to_dColHt_a(i,k) = dS_to_dColHt(i,k)

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
                      (Kddt_h(K-1) / hp_a(i)) * ((T0(k-2) - T0(k-1)) + dTe(k-2))
                dS_km1_t2 = (S0(k)-S0(k-1)) - &
                      (Kddt_h(K-1) / hp_a(i)) * ((S0(k-2) - S0(k-1)) + dSe(k-2))
              endif
              dTe_term = dTe_t2 + hp_a(i) * (T0(k-1)-T0(k))
              dSe_term = dSe_t2 + hp_a(i) * (S0(k-1)-S0(k))
            else
              if (K<=2) then
                Th_a(k-1) = h(i,k-1) * T0(k-1) ; Sh_a(k-1) = h(i,k-1) * S0(k-1)
              else
                Th_a(k-1) = h(i,k-1) * T0(k-1) + Kddt_h(K-1) * Te(k-2)
                Sh_a(k-1) = h(i,k-1) * S0(k-1) + Kddt_h(K-1) * Se(k-2)
              endif
              Th_b(k) = h(i,k) * T0(k) ; Sh_b(k) = h(i,k) * S0(k)
            endif

            !   Using Pr=1 and the diffusivity at the bottom interface (once it is
            ! known), determine how much resolved mean kinetic energy (MKE) will be
            ! extracted within a timestep and add a fraction CS%MKE_to_TKE_effic of
            ! this to the mTKE budget available for mixing in the next layer.

            if ((CS%MKE_to_TKE_effic > 0.0) .and. (htot(i)*h(i,k) > 0.0)) then
              ! This is the energy that would be available from homogenizing the
              ! velocities between layer k and the layers above.
              dMKE_max = (GV%H_to_kg_m2 * CS%MKE_to_TKE_effic) * 0.5 * &
                  (h(i,k) / ((htot(i) + h(i,k))*htot(i))) * &
                  ((uhtot(i)-u(i,k)*htot(i))**2 + (vhtot(i)-v(i,k)*htot(i))**2)
              ! A fraction (1-exp(Kddt_h*MKE2_Hharm)) of this energy would be
              ! extracted by mixing with a finite viscosity.
              MKE2_Hharm = (htot(i) + h(i,k) + 2.0*h_neglect) / &
                           ((htot(i)+h_neglect) * (h(i,k)+h_neglect))
            else
              dMKE_max = 0.0 ; MKE2_Hharm = 0.0
            endif

            ! At this point, Kddt_h(K) will be unknown because its value may depend
            ! on how much energy is available.  mech_TKE might be negative due to
            ! contributions from TKE_forced.
            h_tt = htot(i) + h_tt_min
            TKE_here = mech_TKE(i) + CS%wstar_ustar_coef*conv_PErel(i)
            if (TKE_here > 0.0) then
              vstar = CS%vstar_scale_fac * (I_dtrho*TKE_here)**C1_3
              hbs_here = GV%H_to_m * min(hb_hs(i,K), MixLen_shape(K))
              Mixing_Length_Used(k) = MAX(CS%min_mix_len,((h_tt*hbs_here)*vstar) / &
                  ((CS%Ekman_scale_coef * absf(i)) * (h_tt*hbs_here) + vstar))
              !Note setting Kd_guess0 to Mixing_Length_Used(K) here will
              ! change the answers.  Therefore, skipping that.
              if (.not.CS%Use_MLD_Iteration) then
                 Kd_guess0 = vstar * vonKar *  ((h_tt*hbs_here)*vstar) / &
                  ((CS%Ekman_scale_coef * absf(i)) * (h_tt*hbs_here) + vstar)
              else
                 Kd_guess0 = vstar * vonKar * Mixing_Length_Used(k)
              endif
            else
              vstar = 0.0 ; Kd_guess0 = 0.0
            endif
            Vstar_Used(k) = vstar ! Track vstar
            Kddt_h_g0 = Kd_guess0*dt_h

            if (CS%orig_PE_calc) then
              call find_PE_chg_orig(Kddt_h_g0, h(i,k), hp_a(i), dTe_term, dSe_term, &
                       dT_km1_t2, dS_km1_t2, dT_to_dPE(i,k), dS_to_dPE(i,k), &
                       dT_to_dPE_a(i,k-1), dS_to_dPE_a(i,k-1), &
                       pres(i,K), dT_to_dColHt(i,k), dS_to_dColHt(i,k), &
                       dT_to_dColHt_a(i,k-1), dS_to_dColHt_a(i,k-1), &
                       PE_chg=PE_chg_g0, dPEc_dKd=dPEa_dKd_g0, dPE_max=PE_chg_max, &
                       dPEc_dKd_0=dPEc_dKd_Kd0 )
            else
              call find_PE_chg(0.0, Kddt_h_g0, hp_a(i), h(i,k), &
                         Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                         dT_to_dPE_a(i,k-1), dS_to_dPE_a(i,k-1), dT_to_dPE(i,k), dS_to_dPE(i,k), &
                         pres(i,K), dT_to_dColHt_a(i,k-1), dS_to_dColHt_a(i,k-1), &
                         dT_to_dColHt(i,k), dS_to_dColHt(i,k), &
                         PE_chg=PE_chg_g0, dPEc_dKd=dPEa_dKd_g0, dPE_max=PE_chg_max, &
                         dPEc_dKd_0=dPEc_dKd_Kd0 )
            endif

            MKE_src = dMKE_max*(1.0 - exp(-Kddt_h_g0 * MKE2_Hharm))

            if (pe_chg_g0 .gt. 0.0) then
              !Negative buoyancy (increases PE)
              N2_dissipation = 1.+CS%N2_DISSIPATION_SCALE_NEG
            else
              !Positive buoyancy (decreases PE)
              N2_dissipation = 1.+CS%N2_DISSIPATION_SCALE_POS
            endif

            if ((PE_chg_g0 < 0.0) .or. ((vstar == 0.0) .and. (dPEc_dKd_Kd0 < 0.0))) then
              ! This column is convectively unstable.
              if (PE_chg_max <= 0.0) then
                ! Does MKE_src need to be included in the calculation of vstar here?
                TKE_here = mech_TKE(i) + CS%wstar_ustar_coef*(conv_PErel(i)-PE_chg_max)
                if (TKE_here > 0.0) then
                  vstar = CS%vstar_scale_fac * (I_dtrho*TKE_here)**C1_3
                  hbs_here = GV%H_to_m * min(hb_hs(i,K), MixLen_shape(K))
                  Mixing_Length_Used(k) = max(CS%min_mix_len,((h_tt*hbs_here)*vstar) / &
                      ((CS%Ekman_scale_coef * absf(i)) * (h_tt*hbs_here) + vstar))
                  if (.not.CS%Use_MLD_Iteration) then
                  ! Note again (as prev) that using Mixing_Length_Used here
                  !  instead of redoing the computation will change answers...
                     Kd(i,k) = vstar * vonKar *  ((h_tt*hbs_here)*vstar) / &
                          ((CS%Ekman_scale_coef * absf(i)) * (h_tt*hbs_here) + vstar)
                  else
                     Kd(i,k) = vstar * vonKar * Mixing_Length_Used(k)
                  endif
                else
                  vstar = 0.0 ; Kd(i,k) = 0.0
                endif
                Vstar_Used(k) = vstar

                if (CS%orig_PE_calc) then
                  call find_PE_chg_orig(Kd(i,k)*dt_h, h(i,k), hp_a(i), dTe_term, dSe_term, &
                           dT_km1_t2, dS_km1_t2, dT_to_dPE(i,k), dS_to_dPE(i,k), &
                           dT_to_dPE_a(i,k-1), dS_to_dPE_a(i,k-1), &
                           pres(i,K), dT_to_dColHt(i,k), dS_to_dColHt(i,k), &
                           dT_to_dColHt_a(i,k-1), dS_to_dColHt_a(i,k-1), &
                           PE_chg=dPE_conv)
                else
                  call find_PE_chg(0.0, Kd(i,k)*dt_h, hp_a(i), h(i,k), &
                           Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                           dT_to_dPE_a(i,k-1), dS_to_dPE_a(i,k-1), dT_to_dPE(i,k), dS_to_dPE(i,k), &
                           pres(i,K), dT_to_dColHt_a(i,k-1), dS_to_dColHt_a(i,k-1), &
                           dT_to_dColHt(i,k), dS_to_dColHt(i,k), &
                           PE_chg=dPE_conv)
                endif
                ! Should this be iterated to convergence for Kd?
                if (dPE_conv > 0.0) then
                  Kd(i,k) = Kd_guess0 ; dPE_conv = PE_chg_g0
                else
                  MKE_src = dMKE_max*(1.0 - exp(-(Kd(i,k)*dt_h) * MKE2_Hharm))
                endif
              else
                ! The energy change does not vary monotonically with Kddt_h.  Find the maximum?
                Kd(i,k) = Kd_guess0 ; dPE_conv = PE_chg_g0
              endif
              conv_PErel(i) = conv_PErel(i) - dPE_conv
              mech_TKE(i) = mech_TKE(i) + MKE_src
              if (CS%TKE_diagnostics) then
                dTKE_conv = dTKE_conv - CS%nstar*dPE_conv * IdtdR0
                dTKE_MKE = dTKE_MKE + MKE_src * IdtdR0
              endif
              if (sfc_connected(i)) then
                CS%ML_depth(i,J) = CS%ML_depth(i,J) + GV%H_to_m * h(i,k)
                !CS%ML_depth2(i,j) = CS%ML_depth2(i,J) + GV%H_to_m * h(i,k)
              endif

              Kddt_h(K) = Kd(i,k)*dt_h
            elseif (tot_TKE + (MKE_src - N2_DISSIPATION*PE_chg_g0) >= 0.0) then
              ! There is energy to support the suggested mixing.  Keep that estimate.
              Kd(i,k) = Kd_guess0
              Kddt_h(K) = Kddt_h_g0

              ! Reduce the mechanical and convective TKE proportionately.
              tot_TKE = tot_TKE + MKE_src
              TKE_reduc = 0.0   ! tot_TKE could be 0 if Convectively_stable is false.
              if (tot_TKE > 0.0) TKE_reduc = (tot_TKE - N2_DISSIPATION*PE_chg_g0) &
                                             / tot_TKE
              if (CS%TKE_diagnostics) then
                dTKE_mixing = dTKE_mixing - PE_chg_g0 * IdtdR0
                dTKE_MKE = dTKE_MKE + MKE_src * IdtdR0
                dTKE_conv_decay = dTKE_conv_decay + &
                    (1.0-TKE_reduc)*(CS%nstar-nstar_FC) * conv_PErel(i) * IdtdR0
              endif
              tot_TKE = TKE_reduc*tot_TKE
              mech_TKE(i) = TKE_reduc*(mech_TKE(i) + MKE_src)
              conv_PErel(i) = TKE_reduc*conv_PErel(i)
              if (sfc_connected(i)) then
                CS%ML_depth(i,J) = CS%ML_depth(i,J) + GV%H_to_m * h(i,k)
                !CS%ML_depth2(i,J) = CS%ML_depth2(i,J) + GV%H_to_m * h(i,k)
              endif
            elseif (tot_TKE == 0.0) then
              ! This can arise if nstar_FC = 0.
              Kd(i,k) = 0.0 ; Kddt_h(K) = 0.0
              tot_TKE = 0.0 ; conv_PErel(i) = 0.0 ; mech_TKE(i) = 0.0
              sfc_disconnect = .true.
            else
              ! There is not enough energy to support the mixing, so reduce the
              ! diffusivity to what can be supported.
              Kddt_h_max = Kddt_h_g0 ; Kddt_h_min = 0.0
              TKE_left_max = tot_TKE + (MKE_src - N2_DISSIPATION*PE_chg_g0) ;
              TKE_left_min = tot_TKE

              ! As a starting guess, take the minimum of a false position estimate
              ! and a Newton's method estimate starting from Kddt_h = 0.0.
              Kddt_h_guess = tot_TKE * Kddt_h_max / max( N2_DISSIPATION*PE_chg_g0 &
                                 - MKE_src, Kddt_h_max * (dPEc_dKd_Kd0 - dMKE_max *        &
                                 MKE2_Hharm) )
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
                  call find_PE_chg_orig(Kddt_h_guess, h(i,k), hp_a(i), dTe_term, dSe_term, &
                           dT_km1_t2, dS_km1_t2, dT_to_dPE(i,k), dS_to_dPE(i,k), &
                           dT_to_dPE_a(i,k-1), dS_to_dPE_a(i,k-1), &
                           pres(i,K), dT_to_dColHt(i,k), dS_to_dColHt(i,k), &
                           dT_to_dColHt_a(i,k-1), dS_to_dColHt_a(i,k-1), &
                           PE_chg=PE_chg, dPEc_dKd=dPEc_dKd )
                else
                  call find_PE_chg(0.0, Kddt_h_guess, hp_a(i), h(i,k), &
                           Th_a(k-1), Sh_a(k-1), Th_b(k), Sh_b(k), &
                           dT_to_dPE_a(i,k-1), dS_to_dPE_a(i,k-1), dT_to_dPE(i,k), dS_to_dPE(i,k), &
                           pres(i,K), dT_to_dColHt_a(i,k-1), dS_to_dColHt_a(i,k-1), &
                           dT_to_dColHt(i,k), dS_to_dColHt(i,k), &
                           PE_chg=dPE_conv)
                endif
                MKE_src = dMKE_max * (1.0 - exp(-MKE2_Hharm * Kddt_h_guess))
                dMKE_src_dK = dMKE_max * MKE2_Hharm * exp(-MKE2_Hharm * Kddt_h_guess)

                TKE_left = tot_TKE + (MKE_src - N2_DISSIPATION*PE_chg)
                if (debug) then
                  Kddt_h_itt(itt) = Kddt_h_guess ; MKE_src_itt(itt) = MKE_src
                  PE_chg_itt(itt) = N2_DISSIPATION*PE_chg
                  TKE_left_itt(itt) = TKE_left
                  dPEa_dKd_itt(itt) = dPEc_dKd
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
                if (dPEc_dKd*N2_DISSIPATION - dMKE_src_dK <= 0.0) then
                  use_Newt = .false.
                else
                  dKddt_h_Newt = TKE_left / (dPEc_dKd*N2_DISSIPATION - dMKE_src_dK)
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
              enddo
              Kd(i,K) = Kddt_h_guess / dt_h ; Kddt_h(K) = Kd(i,K)*dt_h

              ! All TKE should have been consumed.
              if (CS%TKE_diagnostics) then
                dTKE_mixing = dTKE_mixing - (tot_TKE + MKE_src) * IdtdR0
                dTKE_MKE = dTKE_MKE + MKE_src * IdtdR0
                dTKE_conv_decay = dTKE_conv_decay + &
                    (CS%nstar-nstar_FC) * conv_PErel(i) * IdtdR0
              endif

              if (sfc_connected(i)) CS%ML_depth(i,J) = CS%ML_depth(i,J) + &
                   (PE_chg / PE_chg_g0) * GV%H_to_m * h(i,k)
              tot_TKE = 0.0 ; mech_TKE(i) = 0.0 ; conv_PErel(i) = 0.0
              sfc_disconnect = .true.
            endif

            Kddt_h(K) = Kd(i,K)*dt_h
           !   At this point, the final value of Kddt_h(K) is known, so the
           ! estimated properties for layer k-1 can be calculated.
            b1 = 1.0 / (hp_a(i) + Kddt_h(K))
            c1(K) = Kddt_h(K) * b1
            if (CS%orig_PE_calc) then
              dTe(k-1) = b1 * ( Kddt_h(K)*(T0(k)-T0(k-1)) + dTe_t2 )
              dSe(k-1) = b1 * ( Kddt_h(K)*(S0(k)-S0(k-1)) + dSe_t2 )
            endif

            hp_a(i) = h(i,k) + (hp_a(i) * b1) * Kddt_h(K)
            dT_to_dPE_a(i,k) = dT_to_dPE(i,k) + c1(K)*dT_to_dPE_a(i,k-1)
            dS_to_dPE_a(i,k) = dS_to_dPE(i,k) + c1(K)*dS_to_dPE_a(i,k-1)
            dT_to_dColHt_a(i,k) = dT_to_dColHt(i,k) + c1(K)*dT_to_dColHt_a(i,k-1)
            dS_to_dColHt_a(i,k) = dS_to_dColHt(i,k) + c1(K)*dS_to_dColHt_a(i,k-1)

          endif  ! tot_TKT > 0.0 branch.  Kddt_h(K) has been set.

          ! Store integrated velocities and thicknesses for MKE conversion calculations.
          if (sfc_disconnect) then
            ! There is no turbulence at this interface, so zero out the running sums.
            uhtot(i) = u(i,k)*h(i,k)
            vhtot(i) = v(i,k)*h(i,k)
            htot(i)  = h(i,k)
            sfc_connected(i) = .false.
          else
            uhtot(i) = uhtot(i) + u(i,k)*h(i,k)
            vhtot(i) = vhtot(i) + v(i,k)*h(i,k)
            htot(i)  = htot(i) + h(i,k)
          endif

          if (debug) then
            if (k==2) then
              Te(1) = b1*(h(i,1)*T0(1))
              Se(1) = b1*(h(i,1)*S0(1))
            else
              Te(k-1) = b1 * (h(i,k-1) * T0(k-1) + Kddt_h(K-1) * Te(k-2))
              Se(k-1) = b1 * (h(i,k-1) * S0(k-1) + Kddt_h(K-1) * Se(k-2))
            endif
          endif
        enddo
        Kd(i,nz+1) = 0.0

        if (debug) then
          ! Complete the tridiagonal solve for Te.
          b1 = 1.0 / hp_a(i)
          Te(nz) = b1 * (h(i,nz) * T0(nz) + Kddt_h(nz) * Te(nz-1))
          Se(nz) = b1 * (h(i,nz) * S0(nz) + Kddt_h(nz) * Se(nz-1))
          do k=nz-1,1,-1
            Te(k) = Te(k) + c1(K+1)*Te(k+1)
            Se(k) = Se(k) + c1(K+1)*Se(k+1)
          enddo
        endif
        if (present(dT_expected)) then
          do k=1,nz ; dT_expected(i,j,k) = Te(k) - T0(k) ; enddo
        endif
        if (present(dS_expected)) then
          do k=1,nz ; dS_expected(i,j,k) = Se(k) - S0(k) ; enddo
        endif
        if (debug) then
          dPE_debug = 0.0
          do k=1,nz
            dPE_debug = dPE_debug + (dT_to_dPE(i,k) * (Te(k) - T0(k)) + &
                                     dS_to_dPE(i,k) * (Se(k) - S0(k)))
          enddo
          mixing_debug = dPE_debug * IdtdR0
        endif
        k = nz ! This is here to allow a breakpoint to be set.
        !/BGR
        ! The following lines are used for the iteration
        ! note the iteration has been altered to use the value predicted by
        ! the TKE threshold (ML_DEPTH).  This is because the MSTAR
        ! is now dependent on the ML, and therefore the ML needs to be estimated
        ! more precisely than the grid spacing.
        !/
        ITmax(obl_it) = max_MLD       ! Track max    }
        ITmin(obl_it) = min_MLD       ! Track min    } For debug purpose
        ITguess(obl_it) = MLD_guess   ! Track guess  }
        !/
        MLD_FOUND=0.0 ; FIRST_OBL=.true.
        if (CS%Orig_MLD_iteration) then
          !This is how the iteration was original conducted
          do k=2,nz
            if (FIRST_OBL) then !Breaks when OBL found
              if (Vstar_Used(k) > 1.e-10 .and. k < nz) then
                MLD_FOUND = MLD_FOUND+h(i,k-1)*GV%H_to_m
              else
                FIRST_OBL = .false.
                if (MLD_FOUND-CS%MLD_tol > MLD_guess) then
                  min_MLD = MLD_guess
                elseif ((MLD_guess-MLD_FOUND) < max(CS%MLD_tol,h(i,k-1)*GV%H_to_m)) then
                  OBL_CONVERGED = .true.!Break convergence loop
                  if (OBL_IT_STATS) then !Compute iteration statistics
                    MAXIT = max(MAXIT,obl_it)
                    MINIT = min(MINIT,obl_it)
                    SUMIT = SUMIT+obl_it
                    NUMIT = NUMIT+1
                    print*,MAXIT,MINIT,SUMIT/NUMIT
                  endif
                  CS%ML_Depth2(i,j) = MLD_guess
                else
                  max_MLD = MLD_guess !We know this guess was too deep
                endif
              endif
            endif
          enddo
        else
          !New method uses ML_DEPTH as computed in ePBL routine
          MLD_FOUND=CS%ML_DEPTH(i,j)
          if (MLD_FOUND-CS%MLD_tol > MLD_guess) then
            min_MLD = MLD_guess
          elseif (abs(MLD_guess-MLD_FOUND) < (CS%MLD_tol)) then
            OBL_CONVERGED = .true.!Break convergence loop
            if (OBL_IT_STATS) then !Compute iteration statistics
              MAXIT = max(MAXIT,obl_it)
              MINIT = min(MINIT,obl_it)
              SUMIT = SUMIT+obl_it
              NUMIT = NUMIT+1
              print*,MAXIT,MINIT,SUMIT/NUMIT
            endif
            CS%ML_Depth2(i,j) = MLD_guess
          else
            max_MLD = MLD_guess !We know this guess was too deep
          endif
        endif
        ! For next pass, guess average of minimum and maximum values.
        MLD_guess = min_MLD*0.5 + max_MLD*0.5
        ITresult(obl_it) = MLD_FOUND
      endif ; enddo ! Iteration loop for converged boundary layer thickness.
      if (.not.OBL_CONVERGED) then
        !/Temp output, warn that EPBL didn't converge
        !/Print guess/found for every iteration step
        !print*,'EPBL MLD DID NOT CONVERGE'
        NOTCONVERGED=NOTCONVERGED+1
        !do obl_it=1,max_obl_it
        !   print*,ITmin(obl_it),ITmax(obl_it)
        !   print*,obl_it,ITguess(obl_it),ITresult(obl_it)
        !enddo
        !Activate to print out some output when not converged
        !{
        !print*,'Min/Max: ',ITmin(50),ITmax(50)
        !print*,'Guess/result: ',ITguess(50),ITresult(50)
        !print*,'Stats on CPU: ',CONVERGED,NOTCONVERGED,&
        !     real(NOTCONVERGED)/real(CONVERGED)
        !}
        !stop !Kill if not converged during testing.
      else
        CONVERGED=CONVERGED+1
      endif

      if (CS%TKE_diagnostics) then
        CS%diag_TKE_MKE(i,j) = CS%diag_TKE_MKE(i,j) + dTKE_MKE
        CS%diag_TKE_conv(i,j) = CS%diag_TKE_conv(i,j) + dTKE_conv
        CS%diag_TKE_forcing(i,j) = CS%diag_TKE_forcing(i,j) + dTKE_forcing
        CS%diag_TKE_mixing(i,j) = CS%diag_TKE_mixing(i,j) + dTKE_mixing
        CS%diag_TKE_mech_decay(i,j) = CS%diag_TKE_mech_decay(i,j) + dTKE_mech_decay
        CS%diag_TKE_conv_decay(i,j) = CS%diag_TKE_conv_decay(i,j) + dTKE_conv_decay
       ! CS%diag_TKE_unbalanced_forcing(i,j) = CS%diag_TKE_unbalanced_forcing(i,j) + dTKE_unbalanced
      endif
      if (CS%Mixing_Diagnostics) then
        !Write to 3-D for outputing Mixing length and
        !  velocity scale.
        do k=1,nz
          CS%Mixing_Length(i,j,k) = Mixing_Length_Used(k)
          CS%Velocity_Scale(i,j,k) = Vstar_Used(k)
        enddo
      endif
      CS%Enhance_V(i,j) = Enhance_V
      CS%MSTAR_MIX(i,j) = MSTAR_MIX
    else
      ! For masked points, Kd_int must still be set (to 0) because it has intent(out).
      do K=1,nz+1
        Kd(i,K) = 0.
      enddo
      if (present(dT_expected)) then
        do k=1,nz ; dT_expected(i,j,k) = 0.0 ; enddo
      endif
      if (present(dS_expected)) then
        do k=1,nz ; dS_expected(i,j,k) = 0.0 ; enddo
      endif
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
      call post_data(CS%id_ML_depth, CS%ML_depth, CS%diag)
    if (CS%id_TKE_wind > 0) &
      call post_data(CS%id_TKE_wind, CS%diag_TKE_wind, CS%diag)
    if (CS%id_TKE_MKE > 0) &
      call post_data(CS%id_TKE_MKE, CS%diag_TKE_MKE, CS%diag)
    if (CS%id_TKE_conv > 0) &
      call post_data(CS%id_TKE_conv, CS%diag_TKE_conv, CS%diag)
    if (CS%id_TKE_forcing > 0) &
      call post_data(CS%id_TKE_forcing, CS%diag_TKE_forcing, CS%diag)
    if (CS%id_TKE_mixing > 0) &
      call post_data(CS%id_TKE_mixing, CS%diag_TKE_mixing, CS%diag)
    if (CS%id_TKE_mech_decay > 0) &
      call post_data(CS%id_TKE_mech_decay, CS%diag_TKE_mech_decay, CS%diag)
    if (CS%id_TKE_conv_decay > 0) &
      call post_data(CS%id_TKE_conv_decay, CS%diag_TKE_conv_decay, CS%diag)
    if (CS%id_Hsfc_used > 0) &
      call post_data(CS%id_Hsfc_used, Hsfc_used, CS%diag)
    if (CS%id_Mixing_Length > 0) &
      call post_data(CS%id_Mixing_Length, CS%Mixing_Length, CS%diag)
    if (CS%id_Velocity_Scale >0) &
      call post_data(CS%id_Velocity_Scale, CS%Velocity_Scale, CS%diag)
    if (CS%id_OSBL >0) &
      call post_data(CS%id_OSBL, CS%ML_Depth2, CS%diag)
    if (CS%id_LT_Enhancement >0) &
      call post_data(CS%id_LT_Enhancement, CS%Enhance_V, CS%diag)
    if (CS%id_MSTAR_MIX >0) &
      call post_data(CS%id_MSTAR_MIX, CS%MSTAR_MIX, CS%diag)
  endif

end subroutine energetic_PBL

!> This subroutine calculates the change in potential energy and or derivatives
!! for several changes in an interfaces's diapycnal diffusivity times a timestep.
subroutine find_PE_chg(Kddt_h0, dKddt_h, hp_a, hp_b, Th_a, Sh_a, Th_b, Sh_b, &
                       dT_to_dPE_a, dS_to_dPE_a, dT_to_dPE_b, dS_to_dPE_b, &
                       pres, dT_to_dColHt_a, dS_to_dColHt_a, dT_to_dColHt_b, dS_to_dColHt_b, &
                       PE_chg, dPEc_dKd, dPE_max, dPEc_dKd_0, ColHt_cor)
  real, intent(in)  :: Kddt_h0  !< The previously used diffusivity at an interface times
                                !! the time step and  divided by the average of the
                                !! thicknesses around the interface, in units of H (m or kg-2).
  real, intent(in)  :: dKddt_h  !< The trial change in the diffusivity at an interface times
                                !! the time step and  divided by the average of the
                                !! thicknesses around the interface, in units of H (m or kg-2).
  real, intent(in)  :: hp_a     !< The effective pivot thickness of the layer above the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above, in H.
  real, intent(in)  :: hp_b     !< The effective pivot thickness of the layer below the
                                !! interface, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above, in H.
  real, intent(in)  :: Th_a     !< An effective temperature times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers, in K H.
  real, intent(in)  :: Sh_a     !< An effective salinity times a thickness in the layer
                                !! above, including implicit mixing effects with other
                                !! yet higher layers, in K H.
  real, intent(in)  :: Th_b     !< An effective temperature times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers, in K H.
  real, intent(in)  :: Sh_b     !< An effective salinity times a thickness in the layer
                                !! below, including implicit mixing effects with other
                                !! yet lower layers, in K H.
  real, intent(in)  :: dT_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! potential energy, including all implicit diffusive changes
                                !! in the temperatures of all the layers above, in J m-2 K-1.
  real, intent(in)  :: dS_to_dPE_a !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! potential energy, including all implicit diffusive changes
                                !! in the salinities of all the layers above, in J m-2 ppt-1.
  real, intent(in)  :: dT_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! potential energy, including all implicit diffusive changes
                                !! in the temperatures of all the layers below, in J m-2 K-1.
  real, intent(in)  :: dS_to_dPE_b !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! potential energy, including all implicit diffusive changes
                                !! in the salinities of all the layers below, in J m-2 ppt-1.
  real, intent(in)  :: pres     !< The hydrostatic interface pressure, which is used to relate
                                !! the changes in column thickness to the energy that is radiated
                                !! as gravity waves and unavailable to drive mixing, in Pa.
  real, intent(in)  :: dT_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers above, in m K-1.
  real, intent(in)  :: dS_to_dColHt_a !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers above, in m ppt-1.
  real, intent(in)  :: dT_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                !! a layer's temperature change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the temperatures of all the layers below, in m K-1.
  real, intent(in)  :: dS_to_dColHt_b !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                !! a layer's salinity change to the change in column
                                !! height, including all implicit diffusive changes
                                !! in the salinities of all the layers below, in m ppt-1.

  real, optional, intent(out) :: PE_chg   !< The change in column potential energy from applying
                                          !! Kddt_h at the present interface, in J m-2.
  real, optional, intent(out) :: dPEc_dKd !< The partial derivative of PE_chg with Kddt_h,
                                          !! in units of J m-2 H-1.
  real, optional, intent(out) :: dPE_max  !< The maximum change in column potential energy that could
                                          !! be realizedd by applying a huge value of Kddt_h at the
                                          !! present interface, in J m-2.
  real, optional, intent(out) :: dPEc_dKd_0 !< The partial derivative of PE_chg with Kddt_h in the
                                            !! limit where Kddt_h = 0, in J m-2 H-1.
  real, optional, intent(out) :: ColHt_cor  !< The correction to PE_chg that is made due to a net
                                            !! change in the column height, in J m-2.

  real :: hps ! The sum of the two effective pivot thicknesses, in H.
  real :: bdt1 ! A product of the two pivot thicknesses plus a diffusive term, in H2.
  real :: dT_c ! The core term in the expressions for the temperature changes, in K H2.
  real :: dS_c ! The core term in the expressions for the salinity changes, in psu H2.
  real :: PEc_core ! The diffusivity-independent core term in the expressions
                   ! for the potential energy changes, J m-3.
  real :: ColHt_core ! The diffusivity-independent core term in the expressions
                     ! for the column height changes, J m-3.
  real :: ColHt_chg  ! The change in the column height, in m.
  real :: y1   ! A local temporary term, in units of H-3 or H-4 in various contexts.

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
    y1 = dKddt_h / (bdt1 * (bdt1 + dKddt_h * hps))
    PE_chg = PEc_core * y1
    ColHt_chg = ColHt_core * y1
    if (ColHt_chg < 0.0) PE_chg = PE_chg - pres * ColHt_chg
    if (present(ColHt_cor)) ColHt_cor = -pres * min(ColHt_chg, 0.0)
  else if (present(ColHt_cor)) then
    y1 = dKddt_h / (bdt1 * (bdt1 + dKddt_h * hps))
    ColHt_cor = -pres * min(ColHt_core * y1, 0.0)
  endif

  if (present(dPEc_dKd)) then
    ! Find the derivative of the potential energy change with dKddt_h.
    y1 = 1.0 / (bdt1 + dKddt_h * hps)**2
    dPEc_dKd = PEc_core * y1
    ColHt_chg = ColHt_core * y1
    if (ColHt_chg < 0.0) dPEc_dKd = dPEc_dKd - pres * ColHt_chg
  endif

  if (present(dPE_max)) then
    ! This expression is the limit of PE_chg for infinite dKddt_h.
    y1 = 1.0 / (bdt1 * hps)
    dPE_max = PEc_core * y1
    ColHt_chg = ColHt_core * y1
    if (ColHt_chg < 0.0) dPE_max = dPE_max - pres * ColHt_chg
  endif

  if (present(dPEc_dKd_0)) then
    ! This expression is the limit of dPEc_dKd for dKddt_h = 0.
    y1 = 1.0 / bdt1**2
    dPEc_dKd_0 = PEc_core * y1
    ColHt_chg = ColHt_core * y1
    if (ColHt_chg < 0.0) dPEc_dKd_0 = dPEc_dKd_0 - pres * ColHt_chg
  endif

end subroutine find_PE_chg

!> This subroutine calculates the change in potential energy and or derivatives
!! for several changes in an interfaces's diapycnal diffusivity times a timestep
!! using the original form used in the first version of ePBL.
subroutine find_PE_chg_orig(Kddt_h, h_k, b_den_1, dTe_term, dSe_term, &
                       dT_km1_t2, dS_km1_t2, dT_to_dPE_k, dS_to_dPE_k, &
                       dT_to_dPEa, dS_to_dPEa, pres, dT_to_dColHt_k, &
                       dS_to_dColHt_k, dT_to_dColHta, dS_to_dColHta, &
                       PE_chg, dPEc_dKd, dPE_max, dPEc_dKd_0)
  real, intent(in)  :: Kddt_h   !< The diffusivity at an interface times the time step and
                                !! divided by the average of the thicknesses around the
                                !! interface, in units of H (m or kg-2).
  real, intent(in)  :: h_k      !< The thickness of the layer below the interface, in H.
  real, intent(in)  :: b_den_1  !< The first term in the denominator of the pivot
                                !! for the tridiagonal solver, given by h_k plus a term that
                                !! is a fraction (determined from the tridiagonal solver) of
                                !! Kddt_h for the interface above, in H.
  real, intent(in)  :: dTe_term !< A diffusivity-independent term related to the
                                !! temperature change in the layer below the interface, in K H.
  real, intent(in)  :: dSe_term !< A diffusivity-independent term related to the
                                !! salinity change in the layer below the interface, in ppt H.
  real, intent(in)  :: dT_km1_t2 !< A diffusivity-independent term related to the
                                 !! temperature change in the layer above the interface, in K.
  real, intent(in)  :: dS_km1_t2 !< A diffusivity-independent term related to the
                                 !! salinity change in the layer above the interface, in ppt.
  real, intent(in)  :: pres      !< The hydrostatic interface pressure, which is used to relate
                                 !! the changes in column thickness to the energy that is radiated
                                 !! as gravity waves and unavailable to drive mixing, in Pa.
  real, intent(in)  :: dT_to_dPE_k !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! potential energy, including all implicit diffusive changes
                                 !! in the temperatures of all the layers below, in J m-2 K-1.
  real, intent(in)  :: dS_to_dPE_k !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! potential energy, including all implicit diffusive changes
                                 !! in the salinities of all the layers below, in J m-2 ppt-1.
  real, intent(in)  :: dT_to_dPEa !< A factor (pres_lay*mass_lay*dSpec_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! potential energy, including all implicit diffusive changes
                                 !! in the temperatures of all the layers above, in J m-2 K-1.
  real, intent(in)  :: dS_to_dPEa !< A factor (pres_lay*mass_lay*dSpec_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! potential energy, including all implicit diffusive changes
                                 !! in the salinities of all the layers above, in J m-2 ppt-1.
  real, intent(in)  :: dT_to_dColHt_k !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the temperatures of all the layers below, in m K-1.
  real, intent(in)  :: dS_to_dColHt_k !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the salinities of all the layers below, in m ppt-1.
  real, intent(in)  :: dT_to_dColHta !< A factor (mass_lay*dSColHtc_vol/dT) relating
                                 !! a layer's temperature change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the temperatures of all the layers above, in m K-1.
  real, intent(in)  :: dS_to_dColHta !< A factor (mass_lay*dSColHtc_vol/dS) relating
                                 !! a layer's salinity change to the change in column
                                 !! height, including all implicit diffusive changes
                                 !! in the salinities of all the layers above, in m ppt-1.

  real, optional, intent(out) :: PE_chg   !< The change in column potential energy from applying
                                          !! Kddt_h at the present interface, in J m-2.
  real, optional, intent(out) :: dPEc_dKd !< The partial derivative of PE_chg with Kddt_h,
                                          !! in units of J m-2 H-1.
  real, optional, intent(out) :: dPE_max  !< The maximum change in column potential energy that could
                                          !! be realizedd by applying a huge value of Kddt_h at the
                                          !! present interface, in J m-2.
  real, optional, intent(out) :: dPEc_dKd_0 !< The partial derivative of PE_chg with Kddt_h in the
                                            !! limit where Kddt_h = 0, in J m-2 H-1.

!   This subroutine determines the total potential energy change due to mixing
! at an interface, including all of the implicit effects of the prescribed
! mixing at interfaces above.  Everything here is derived by careful manipulation
! of the robust tridiagonal solvers used for tracers by MOM6.  The results are
! positive for mixing in a stably stratified environment.
!   The comments describing these arguments are for a downward mixing pass, but
! this routine can also be used for an upward pass with the sense of direction
! reversed.

  real :: b1            ! b1 is used by the tridiagonal solver, in H-1.
  real :: b1Kd          ! Temporary array (nondim.)
  real :: ColHt_chg     ! The change in column thickness in m.
  real :: dColHt_max    ! The change in column thickess for infinite diffusivity, in m.
  real :: dColHt_dKd    ! The partial derivative of column thickess with diffusivity, in s m-1.
  real :: dT_k, dT_km1  ! Temporary arrays in K.
  real :: dS_k, dS_km1  ! Temporary arrays in ppt.
  real :: I_Kr_denom, dKr_dKd   ! Temporary arrays in H-2 and nondim.
  real :: ddT_k_dKd, ddT_km1_dKd ! Temporary arrays in K H-1.
  real :: ddS_k_dKd, ddS_km1_dKd ! Temporary arrays in ppt H-1.

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
    if (ColHt_chg < 0.0) PE_chg = PE_chg - pres * ColHt_chg
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
    if (dColHt_dKd < 0.0) dPEc_dKd = dPEc_dKd - pres * dColHt_dKd
  endif

  if (present(dPE_max)) then
    ! This expression is the limit of PE_chg for infinite Kddt_h.
    dPE_max = (dT_to_dPEa * dT_km1_t2 + dS_to_dPEa * dS_km1_t2) + &
              ((dT_to_dPE_k + dT_to_dPEa) * dTe_term + &
               (dS_to_dPE_k + dS_to_dPEa) * dSe_term) / (b_den_1 + h_k)
    dColHt_max = (dT_to_dColHta * dT_km1_t2 + dS_to_dColHta * dS_km1_t2) + &
              ((dT_to_dColHt_k + dT_to_dColHta) * dTe_term + &
               (dS_to_dColHt_k + dS_to_dColHta) * dSe_term) / (b_den_1 + h_k)
    if (dColHt_max < 0.0) dPE_max = dPE_max - pres*dColHt_max
  endif

  if (present(dPEc_dKd_0)) then
    ! This expression is the limit of dPEc_dKd for Kddt_h = 0.
    dPEc_dKd_0 = (dT_to_dPEa * dT_km1_t2 + dS_to_dPEa * dS_km1_t2) / (b_den_1) + &
                 (dT_to_dPE_k * dTe_term + dS_to_dPE_k * dSe_term) / (h_k*b_den_1)
    dColHt_dKd = (dT_to_dColHta * dT_km1_t2 + dS_to_dColHta * dS_km1_t2) / (b_den_1) + &
                 (dT_to_dColHt_k * dTe_term + dS_to_dColHt_k * dSe_term) / (h_k*b_den_1)
    if (dColHt_dKd < 0.0) dPEc_dKd_0 = dPEc_dKd_0 - pres*dColHt_dKd
  endif

end subroutine find_PE_chg_orig

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

!> Computes wind speed from ustar_air based on COARE 3.5 Cd relationship
subroutine ust_2_u10_coare3p5(USTair,U10,GV)
  real, intent(in)  :: USTair
  type(verticalGrid_type), intent(in) :: GV
  real, intent(out) :: U10
  real, parameter :: vonkar = 0.4
  real, parameter :: nu=1e-6
  real :: z0sm, z0, z0rough, u10a, alpha, CD
  integer :: CT

  ! Uses empirical formula for z0 to convert ustar_air to u10 based on the
  !  COARE 3.5 paper (Edson et al., 2013)
  !alpha=m*U10+b
  !Note in Edson et al. 2013, eq. 13 m is given as 0.017.  However,
  ! m=0.0017 reproduces the curve in their figure 6.

  z0sm = 0.11 * nu / USTair; !Compute z0smooth from ustar guess
  u10 = USTair/sqrt(0.001);  !Guess for u10
  u10a = 1000;

  CT=0
  do while (abs(u10a/u10-1.)>0.001)
    CT=CT+1
    u10a = u10
    alpha = min(0.028,0.0017 * u10 - 0.005)
    z0rough = alpha * USTair**2/GV%g_Earth ! Compute z0rough from ustar guess
    z0=z0sm+z0rough
    CD = ( vonkar / log(10/z0) )**2 ! Compute CD from derived roughness
    u10 = USTair/sqrt(CD);!Compute new u10 from derived CD, while loop
                       ! ends and checks for convergence...CT counter
                       ! makes sure loop doesn't run away if function
                       ! doesn't converge.  This code was produced offline
                       ! and converged rapidly (e.g. 2 cycles)
                       ! for ustar=0.0001:0.0001:10.
    if (CT>20) then
      u10 = USTair/sqrt(0.0015) ! I don't expect to get here, but just
                              !  in case it will output a reasonable value.
      exit
    endif
  enddo
  return
end subroutine ust_2_u10_coare3p5

subroutine get_LA_windsea(ustar, hbl, GV, LA)
! Original description:
! This function returns the enhancement factor, given the 10-meter
! wind (m/s), friction velocity (m/s) and the boundary layer depth (m).
! Update (Jan/25):
! Converted from function to subroutine, now returns Langmuir number.
! Computs 10m wind internally, so only ustar and hbl need passed to
! subroutine.
!
! Qing Li, 160606
! BGR port from CVMix to MOM6 Jan/25/2017
! BGR change output to LA from Efactor
! BGR remove u10 input

! Input
  real, intent(in) :: &
       ! water-side surface friction velocity (m/s)
       ustar, &
       ! boundary layer depth (m)
       hbl
  type(verticalGrid_type), intent(in) :: GV
  real, intent(out) :: LA
! Local variables
  ! parameters
  real, parameter :: &
       ! ratio of U19.5 to U10 (Holthuijsen, 2007)
       u19p5_to_u10 = 1.075, &
       ! ratio of mean frequency to peak frequency for
       ! Pierson-Moskowitz spectrum (Webb, 2011)
       fm_to_fp = 1.296, &
       ! ratio of surface Stokes drift to U10
       us_to_u10 = 0.0162, &
       ! loss ratio of Stokes transport
       r_loss = 0.667
  real :: us, hm0, fm, fp, vstokes, kphil, kstar
  real :: z0, z0i, r1, r2, r3, r4, tmp, us_sl, lasl_sqr_i
  real :: pi, u10
  pi = 4.0*atan(1.0)
  ! Computing u10 based on u_star and COARE 3.5 relationships
  call ust_2_u10_coare3p5(ustar*sqrt(GV%Rho0/1.225),U10,GV)
  if (u10 .gt. 0.0 .and. ustar .gt. 0.0) then
    ! surface Stokes drift
    us = us_to_u10*u10
    !
    ! significant wave height from Pierson-Moskowitz
    ! spectrum (Bouws, 1998)
    hm0 = 0.0246 *u10**2
    !
    ! peak frequency (PM, Bouws, 1998)
    tmp = 2.0 * PI * u19p5_to_u10 * u10
    fp = 0.877 * GV%g_Earth / tmp
    !
    ! mean frequency
    fm = fm_to_fp * fp
    !
    ! total Stokes transport (a factor r_loss is applied to account
    !  for the effect of directional spreading, multidirectional waves
    !  and the use of PM peak frequency and PM significant wave height
    !  on estimating the Stokes transport)
    vstokes = 0.125 * PI * r_loss * fm * hm0**2
    !
    ! the general peak wavenumber for Phillips' spectrum
    ! (Breivik et al., 2016) with correction of directional spreading
    kphil = 0.176 * us / vstokes
    !
    ! surface layer averaged Stokes dirft with Stokes drift profile
    ! estimated from Phillips' spectrum (Breivik et al., 2016)
    ! the directional spreading effect from Webb and Fox-Kemper, 2015
    ! is also included
    kstar = kphil * 2.56
    ! surface layer
    z0 = 0.2 * abs(hbl)
    z0i = 1.0 / z0
    ! term 1 to 4
    r1 = ( 0.151 / kphil * z0i -0.84 ) &
         * ( 1.0 - exp(-2.0 * kphil * z0) )
    r2 = -( 0.84 + 0.0591 / kphil * z0i ) &
         *sqrt( 2.0 * PI * kphil * z0 ) &
         *erfc( sqrt( 2.0 * kphil * z0 ) )
    r3 = ( 0.0632 / kstar * z0i + 0.125 ) &
         * (1.0 - exp(-2.0 * kstar * z0) )
    r4 = ( 0.125 + 0.0946 / kstar * z0i ) &
         *sqrt( 2.0 * PI *kstar * z0) &
         *erfc( sqrt( 2.0 * kstar * z0 ) )
    us_sl = us * (0.715 + r1 + r2 + r3 + r4)
    !
    LA = sqrt(ustar / us_sl)
  else
    LA=1.e8
  endif
endsubroutine Get_LA_windsea

subroutine energetic_PBL_init(Time, G, GV, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(verticalGrid_type), intent(in)    :: GV
  type(param_file_type),   intent(in)    :: param_file
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
  character(len=40)  :: mod = "MOM_energetic_PBL"  ! This module's name.
  real :: omega_frac_dflt
  integer :: isd, ied, jsd, jed
  logical :: use_temperature, use_omega
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (associated(CS)) then
    call MOM_error(WARNING, "mixedlayer_init called with an associated"// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag
  CS%Time => Time

! Set default, read and log parameters
  call log_version(param_file, mod, version, "")

  call get_param(param_file, mod, "MSTAR", CS%mstar, &
                 "The ratio of the friction velocity cubed to the TKE \n"//&
                 "input to the mixed layer.", "units=nondim", default=1.2)
  call get_param(param_file, mod, "MIX_LEN_EXPONENT", CS%MixLenExponent, &
                 "The exponent applied to the ratio of the distance to the MLD \n"//&
                 "and the MLD depth which determines the shape of the mixing length.",&
                 "units=nondim", default=2.0)
  call get_param(param_file, mod, "MSTAR_CAP", CS%mstar_cap, &
                 "Maximum value of mstar allowed in model if non-negative\n"//&
                 "(used if MSTAR_FIXED=false).",&
                 "units=nondim", default=-1.0)
  call get_param(param_file, mod, "MSTAR_SLOPE", CS%mstar_slope, &
                 "The slope of the linear relationship between mstar \n"//&
                 "and the length scale ratio (used if MSTAR_FIXED=false).",&
                 "units=nondim", default=1.0)
  call get_param(param_file, mod, "MSTAR_XINT", CS%mstar_xint, &
                 "The value of the length scale ratio where the mstar \n"//&
                 "is linear above (used if MSTAR_FIXED=false).",&
                 "units=nondim", default=-0.25)
  call get_param(param_file, mod, "MSTAR_AT_XINT", CS%mstar_at_xint, &
                 "The value of mstar at MSTAR_XINT \n"//&
                 "(used if MSTAR_FIXED=false).",&
                 "units=nondim", default=0.13)
  call get_param(param_file, mod, "MSTAR_FIXED", CS%Use_Mstar_Fixed, &
                 "True to use a fixed value of mstar, if false mstar depends \n"//&
                 "on the composite Obhukov length and Ekman length.","units=nondim",&
                 default=.true.)
  call get_param(param_file, mod, "MSTAR_FLATCAP", CS%MSTAR_FLATCAP, &
                 "Set false to use asmptotic cap or defaults to true to use flat cap."&
                 ,"units=nondim",default=.true.)
  call get_param(param_file, mod, "NSTAR", CS%nstar, &
                 "The portion of the buoyant potential energy imparted by \n"//&
                 "surface fluxes that is available to drive entrainment \n"//&
                 "at the base of mixed layer when that energy is positive.", &
                 units="nondim", default=0.2)
  call get_param(param_file, mod, "MKE_TO_TKE_EFFIC", CS%MKE_to_TKE_effic, &
                 "The efficiency with which mean kinetic energy released \n"//&
                 "by mechanically forced entrainment of the mixed layer \n"//&
                 "is converted to turbulent kinetic energy.", units="nondim", &
                 default=0.0)
  call get_param(param_file, mod, "TKE_DECAY", CS%TKE_decay, &
                 "TKE_DECAY relates the vertical rate of decay of the \n"//&
                 "TKE available for mechanical entrainment to the natural \n"//&
                 "Ekman depth.", units="nondim", default=2.5)
!  call get_param(param_file, mod, "HMIX_MIN", CS%Hmix_min, &
!                 "The minimum mixed layer depth if the mixed layer depth \n"//&
!                 "is determined dynamically.", units="m", default=0.0)

  call get_param(param_file, mod, "OMEGA",CS%omega, &
                 "The rotation rate of the earth.", units="s-1", &
                 default=7.2921e-5)
  call get_param(param_file, mod, "ML_USE_OMEGA", use_omega, &
                 "If true, use the absolute rotation rate instead of the \n"//&
                 "vertical component of rotation when setting the decay \n"//&
                 "scale for turbulence.", default=.false., do_not_log=.true.)
  omega_frac_dflt = 0.0
  if (use_omega) then
    call MOM_error(WARNING, "ML_USE_OMEGA is depricated; use ML_OMEGA_FRAC=1.0 instead.")
    omega_frac_dflt = 1.0
  endif
  call get_param(param_file, mod, "ML_OMEGA_FRAC", CS%omega_frac, &
                 "When setting the decay scale for turbulence, use this \n"//&
                 "fraction of the absolute rotation rate blended with the \n"//&
                 "local value of f, as sqrt((1-of)*f^2 + of*4*omega^2).", &
                 units="nondim", default=omega_frac_dflt)
  call get_param(param_file, mod, "WSTAR_USTAR_COEF", CS%wstar_ustar_coef, &
                 "A ratio relating the efficiency with which convectively \n"//&
                 "released energy is converted to a turbulent velocity, \n"//&
                 "relative to mechanically forced TKE. Making this larger \n"//&
                 "increases the BL diffusivity", units="nondim", default=1.0)
  call get_param(param_file, mod, "VSTAR_SCALE_FACTOR", CS%vstar_scale_fac, &
                 "An overall nondimensional scaling factor for v*. \n"//&
                 "Making this larger decreases the PBL diffusivity.", &
                 units="nondim", default=1.0)
  call get_param(param_file, mod, "EKMAN_SCALE_COEF", CS%Ekman_scale_coef, &
                 "A nondimensional scaling factor controlling the inhibition \n"//&
                 "of the diffusive length scale by rotation. Making this larger \n"//&
                 "decreases the PBL diffusivity.", units="nondim", default=1.0)
  call get_param(param_file, mod, "USE_MLD_ITERATION", CS%USE_MLD_ITERATION, &
                 "A logical that specifies whether or not to use the \n"//&
                 "distance to the bottom of the actively turblent boundary \n"//&
                 "layer to help set the EPBL length scale.", default=.false.)
  call get_param(param_file, mod, "ORIG_MLD_ITERATION", CS%ORIG_MLD_ITERATION, &
                 "A logical that specifies whether or not to use the \n"//&
                 "old method for determining MLD depth in iteration, which \n"//&
                 "is limited to resolution.", default=.true.)
  call get_param(param_file, mod, "MLD_ITERATION_GUESS", CS%MLD_ITERATION_GUESS, &
                 "A logical that specifies whether or not to use the \n"//&
                 "previous timestep MLD as a first guess in the MLD iteration.\n"//&
                 "The default is false to facilitate reproducibility.", default=.false.)
  call get_param(param_file, mod, "EPBL_MLD_TOLERANCE", CS%MLD_tol, &
                 "The tolerance for the iteratively determined mixed \n"//&
                 "layer depth.  This is only used with USE_MLD_ITERATION.", &
                 units="meter", default=1.0)
  call get_param(param_file, mod, "EPBL_MIN_MIX_LEN", CS%min_mix_len, &
                 "The minimum mixing length scale that will be used \n"//&
                 "by ePBL.  The default (0) does not set a minimum.", &
                 units="meter", default=0.0)
  call get_param(param_file, mod, "EPBL_ORIGINAL_PE_CALC", CS%orig_PE_calc, &
                 "If true, the ePBL code uses the original form of the \n"//&
                 "potential energy change code.  Otherwise, the newer \n"//&
                 "version that can work with successive increments to the \n"//&
                 "diffusivity in upward or downward passes is used.", default=.true.)
  call get_param(param_file, mod, "EPBL_TRANSITION_SCALE", CS%transLay_scale, &
                 "A scale for the mixing length in the transition layer \n"//&
                 "at the edge of the boundary layer as a fraction of the \n"//&
                 "boundary layer thickness.  The default is 0.1.", &
                 units="nondim", default=0.1)
  if ( CS%USE_MLD_ITERATION .and. abs(CS%transLay_scale-0.5).ge.0.5) then
    call MOM_error(FATAL, "If flag USE_MLD_ITERATION is true, then "//&
                 "EPBL_TRANSITION should be greater than 0 and less than 1.")
  endif
  call get_param(param_file, mod, "N2_DISSIPATION_POS", CS%N2_Dissipation_Scale_Pos, &
                 "A scale for the dissipation of TKE due to stratification \n"//&
                 "in the boundary layer, applied when local stratification \n"//&
                 "is positive.  The default is 0, but should probably be ~0.4.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mod, "N2_DISSIPATION_NEG", CS%N2_Dissipation_Scale_Neg, &
                 "A scale for the dissipation of TKE due to stratification \n"//&
                 "in the boundary layer, applied when local stratification \n"//&
                 "is negative.  The default is 0, but should probably be ~1.", &
                 units="nondim", default=0.0)
   call get_param(param_file, mod, "USE_LA_LI2016", CS%USE_LA_Windsea, &
                 "A logical to use the Li et al. 2016 (submitted) formula to \n"//&
                 " determine the Langmuir number.",&
                 units="nondim", default=.false.)
   call get_param(param_file, mod, "LT_ENHANCE", CS%LT_ENHANCE_FORM, &
                 "Integer for LT enhancement function mode. \n"//&
                 " *Requires USE_LA_LI2016 to be set to True. \n"//&
                 "Options: 0 - no LT enhancement \n"//&
                 "         1 - Van Roekel et al. 2014/Li et al., 2016  \n"//&
                 "         2 - New LES net mixing based derivation (requires coefficients)",&
                 units="nondim", default=0)
   call get_param(param_file, mod, "LT_ENHANCE_COEF", CS%LT_ENHANCE_COEF, &
                 "Coefficient for Langmuir enhancement if LT_ENHANCE = 2",&
                 units="nondim", default=1.57)
   call get_param(param_file, mod, "LT_ENHANCE_EXP", CS%LT_ENHANCE_EXP, &
                 "Exponent for Langmuir enhancement if LT_ENHANCE = 2",&
                 units="nondim", default=-1.5)
   if (CS%LT_ENHANCE_FORM.gt.0 .and. (.not.CS%USE_LA_Windsea)) then
      call MOM_error(FATAL, "If flag USE_LA_LI2016 is false, then "//&
                 " LT_ENHANCE must be 0.")
  endif
  ! This gives a minimum decay scale that is typically much less than Angstrom.
  CS%ustar_min = 2e-4*CS%omega*(GV%Angstrom_z + GV%H_to_m*GV%H_subroundoff)
  call log_param(param_file, mod, "EPBL_USTAR_MIN", CS%ustar_min, &
                 "The (tiny) minimum friction velocity used within the \n"//&
                 "ePBL code, derived from OMEGA and ANGSTROM.", units="meter second-1")

  CS%id_ML_depth = register_diag_field('ocean_model', 'ePBL_h_ML', diag%axesT1, &
      Time, 'Surface mixed layer depth', 'meter')
  CS%id_TKE_wind = register_diag_field('ocean_model', 'ePBL_TKE_wind', diag%axesT1, &
      Time, 'Wind-stirring source of mixed layer TKE', 'meter3 second-3')
  CS%id_TKE_MKE = register_diag_field('ocean_model', 'ePBL_TKE_MKE', diag%axesT1, &
      Time, 'Mean kinetic energy source of mixed layer TKE', 'meter3 second-3')
  CS%id_TKE_conv = register_diag_field('ocean_model', 'ePBL_TKE_conv', diag%axesT1, &
      Time, 'Convective source of mixed layer TKE', 'meter3 second-3')
  CS%id_TKE_forcing = register_diag_field('ocean_model', 'ePBL_TKE_forcing', diag%axesT1, &
      Time, 'TKE consumed by mixing surface forcing or penetrative shortwave radation'//&
            ' through model layers', 'meter3 second-3')
  CS%id_TKE_mixing = register_diag_field('ocean_model', 'ePBL_TKE_mixing', diag%axesT1, &
      Time, 'TKE consumed by mixing that deepens the mixed layer', 'meter3 second-3')
  CS%id_TKE_mech_decay = register_diag_field('ocean_model', 'ePBL_TKE_mech_decay', diag%axesT1, &
      Time, 'Mechanical energy decay sink of mixed layer TKE', 'meter3 second-3')
  CS%id_TKE_conv_decay = register_diag_field('ocean_model', 'ePBL_TKE_conv_decay', diag%axesT1, &
      Time, 'Convective energy decay sink of mixed layer TKE', 'meter3 second-3')
  CS%id_Hsfc_used = register_diag_field('ocean_model', 'ePBL_Hs_used', diag%axesT1, &
      Time, 'Surface region thickness that is used', 'meter')
  CS%id_Mixing_Length = register_diag_field('ocean_model', 'Mixing_Length', diag%axesTi, &
      Time, 'Mixing Length that is used', 'meter')
  CS%id_Velocity_Scale = register_diag_field('ocean_model', 'Velocity_Scale', diag%axesTi, &
      Time, 'Velocity Scale that is used.', 'meter second-1')
  CS%id_LT_enhancement = register_diag_field('ocean_model', 'LT_Enhancement', diag%axesT1, &
      Time, 'LT enhancement that is used.', 'non-dim')
  CS%id_MSTAR_mix = register_diag_field('ocean_model', 'MSTAR', diag%axesT1, &
      Time, 'MSTAR that is used.', 'non-dim')
  CS%id_OSBL = register_diag_field('ocean_model', 'ePBL_OSBL', diag%axesT1, &
      Time, 'Boundary layer depth from the iteration.', 'meter')


  call get_param(param_file, mod, "ENABLE_THERMODYNAMICS", use_temperature, &
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
  call safe_alloc_alloc(CS%Enhance_V, isd, ied, jsd, jed)
  call safe_alloc_alloc(CS%MSTAR_MIX, isd, ied, jsd, jed)


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

end subroutine energetic_PBL_init

subroutine energetic_PBL_end(CS)
  type(energetic_PBL_CS), pointer :: CS

  if (.not.associated(CS)) return

  if (allocated(CS%ML_depth))            deallocate(CS%ML_depth)
  if (allocated(CS%ML_depth2))           deallocate(CS%ML_depth2)
  if (allocated(CS%Enhance_V))           deallocate(CS%Enhance_V)
  if (allocated(CS%MSTAR_MIX))           deallocate(CS%MSTAR_MIX)
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

end module MOM_energetic_PBL
