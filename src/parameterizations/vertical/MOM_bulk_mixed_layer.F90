!> Build mixed layer parameterization
module MOM_bulk_mixed_layer

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : time_type, diag_ctrl, diag_update_remap_grids
use MOM_domains,       only : create_group_pass, do_group_pass, group_pass_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,  only : extractFluxes1d, forcing
use MOM_grid,          only : ocean_grid_type
use MOM_opacity,       only : absorbRemainingSW, optics_type, extract_optics_slice
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>

public bulkmixedlayer, bulkmixedlayer_init

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure with parameters for the MOM_bulk_mixed_layer module
type, public :: bulkmixedlayer_CS ; private
  integer :: nkml            !< The number of layers in the mixed layer.
  integer :: nkbl            !< The number of buffer layers.
  integer :: nsw             !< The number of bands of penetrating shortwave radiation.
  real    :: mstar           !< The ratio of the friction velocity cubed to the
                             !! TKE input to the mixed layer, nondimensional.
  real    :: nstar           !< The fraction of the TKE input to the mixed layer
                             !! available to drive entrainment [nondim].
  real    :: nstar2          !< The fraction of potential energy released by
                             !! convective adjustment that drives entrainment [nondim].
  logical :: absorb_all_SW   !< If true, all shortwave radiation is absorbed by the
                             !! ocean, instead of passing through to the bottom mud.
  real    :: TKE_decay       !< The ratio of the natural Ekman depth to the TKE
                             !! decay scale, nondimensional.
  real    :: bulk_Ri_ML      !< The efficiency with which mean kinetic energy
                             !! released by mechanically forced entrainment of
                             !! the mixed layer is converted to TKE [nondim].
  real    :: bulk_Ri_convective !< The efficiency with which convectively
                             !! released mean kinetic energy becomes TKE [nondim].
  real    :: Hmix_min        !< The minimum mixed layer thickness [H ~> m or kg m-2].
  real    :: H_limit_fluxes  !< When the total ocean depth is less than this
                             !! value [H ~> m or kg m-2], scale away all surface forcing to
                             !! avoid boiling the ocean.
  real    :: ustar_min       !< A minimum value of ustar to avoid numerical problems [Z T-1 ~> m s-1].
                             !! If the value is small enough, this should not affect the solution.
  real    :: omega           !<   The Earth's rotation rate [T-1 ~> s-1].
  real    :: dT_dS_wt        !<   When forced to extrapolate T & S to match the
                             !! layer densities, this factor (in degC / ppt) is
                             !! combined with the derivatives of density with T & S
                             !! to determines what direction is orthogonal to
                             !! density contours.  It should be a typical value of
                             !! (dR/dS) / (dR/dT) in oceanic profiles.
                             !! 6 degC ppt-1 might be reasonable.
  real    :: Hbuffer_min     !< The minimum buffer layer thickness when the mixed layer
                             !! is very large [H ~> m or kg m-2].
  real    :: Hbuffer_rel_min !< The minimum buffer layer thickness relative to the combined
                             !! mixed and buffer layer thicknesses when they are thin [nondim]
  real    :: BL_detrain_time !< A timescale that characterizes buffer layer detrainment
                             !! events [T ~> s].
  real    :: BL_extrap_lim   !< A limit on the density range over which
                             !! extrapolation can occur when detraining from the
                             !! buffer layers, relative to the density range
                             !! within the mixed and buffer layers, when the
                             !! detrainment is going into the lightest interior
                             !! layer  [nondim].
  real :: BL_split_rho_tol   !< The fractional tolerance for matching layer target densities
                             !! when splitting layers to deal with massive interior layers
                             !! that are lighter than one of the mixed or buffer layers [nondim].
  logical :: ML_resort       !<   If true, resort the layers by density, rather than
                             !! doing convective adjustment.
  integer :: ML_presort_nz_conv_adj !< If ML_resort is true, do convective
                             !! adjustment on this many layers (starting from the
                             !! top) before sorting the remaining layers.
  real    :: omega_frac      !<   When setting the decay scale for turbulence, use
                             !! this fraction of the absolute rotation rate blended
                             !! with the local value of f, as sqrt((1-of)*f^2 + of*4*omega^2).
  logical :: correct_absorption !< If true, the depth at which penetrating
                             !! shortwave radiation is absorbed is corrected by
                             !! moving some of the heating upward in the water
                             !! column.  The default is false.
  logical :: Resolve_Ekman   !<   If true, the nkml layers in the mixed layer are
                             !! chosen to optimally represent the impact of the
                             !! Ekman transport on the mixed layer TKE budget.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  logical :: TKE_diagnostics = .false. !< If true, calculate extensive diagnostics of the TKE budget
  logical :: do_rivermix = .false. !< Provide additional TKE to mix river runoff
                             !! at the river mouths to rivermix_depth
  real    :: rivermix_depth = 0.0  !< The depth of mixing if do_rivermix is true [Z ~> m].
  logical :: limit_det       !< If true, limit the extent of buffer layer
                             !! detrainment to be consistent with neighbors.
  real    :: lim_det_dH_sfc  !< The fractional limit in the change between grid
                             !! points of the surface region (mixed & buffer
                             !! layer) thickness [nondim].  0.5 by default.
  real    :: lim_det_dH_bathy !< The fraction of the total depth by which the
                             !! thickness of the surface region (mixed & buffer
                             !! layer) is allowed to change between grid points.
                             !! Nondimensional, 0.2 by default.
  logical :: use_river_heat_content !< If true, use the fluxes%runoff_Hflx field
                             !! to set the heat carried by runoff, instead of
                             !! using SST for temperature of liq_runoff
  logical :: use_calving_heat_content !< Use SST for temperature of froz_runoff
  logical :: salt_reject_below_ML !< It true, add salt below mixed layer (layer mode only)
  logical :: convect_mom_bug !< If true, use code with a bug that causes a loss of momentum
                             !! conservation during mixedlayer convection.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                             !! timing of diagnostic output.
  real    :: Allowed_T_chg   !< The amount by which temperature is allowed
                             !! to exceed previous values during detrainment, K.
  real    :: Allowed_S_chg   !< The amount by which salinity is allowed
                             !! to exceed previous values during detrainment, ppt.

  ! These are terms in the mixed layer TKE budget, all in [Z L2 T-3 ~> m3 s-3] except as noted.
  real, allocatable, dimension(:,:) :: &
    ML_depth, &        !< The mixed layer depth [H ~> m or kg m-2].
    diag_TKE_wind, &   !< The wind source of TKE.
    diag_TKE_RiBulk, & !< The resolved KE source of TKE.
    diag_TKE_conv, &   !< The convective source of TKE.
    diag_TKE_pen_SW, & !< The TKE sink required to mix penetrating shortwave heating.
    diag_TKE_mech_decay, & !< The decay of mechanical TKE.
    diag_TKE_conv_decay, & !< The decay of convective TKE.
    diag_TKE_mixing, & !< The work done by TKE to deepen  the mixed layer.
    diag_TKE_conv_s2, & !< The convective source of TKE due to to mixing in sigma2.
    diag_PE_detrain, & !< The spurious source of potential energy due to mixed layer
                       !! detrainment [R Z L2 T-3 ~> W m-2].
    diag_PE_detrain2   !< The spurious source of potential energy due to mixed layer only
                       !! detrainment [R Z L2 T-3 ~> W m-2].
  logical :: allow_clocks_in_omp_loops  !< If true, clocks can be called from inside loops that can
                                        !! be threaded. To run with multiple threads, set to False.
  type(group_pass_type) :: pass_h_sum_hmbl_prev !< For group halo pass

  !>@{ Diagnostic IDs
  integer :: id_ML_depth = -1, id_TKE_wind = -1, id_TKE_mixing = -1
  integer :: id_TKE_RiBulk = -1, id_TKE_conv = -1, id_TKE_pen_SW = -1
  integer :: id_TKE_mech_decay = -1, id_TKE_conv_decay = -1, id_TKE_conv_s2 = -1
  integer :: id_PE_detrain = -1, id_PE_detrain2 = -1, id_h_mismatch = -1
  integer :: id_Hsfc_used = -1, id_Hsfc_max = -1, id_Hsfc_min = -1
  !!@}
end type bulkmixedlayer_CS

!>@{ CPU clock IDs
integer :: id_clock_detrain=0, id_clock_mech=0, id_clock_conv=0, id_clock_adjustment=0
integer :: id_clock_EOS=0, id_clock_resort=0, id_clock_pass=0
!!@}

contains

!>    This subroutine partially steps the bulk mixed layer model.
!!  The following processes are executed, in the order listed.
!!    1. Undergo convective adjustment into mixed layer.
!!    2. Apply surface heating and cooling.
!!    3. Starting from the top, entrain whatever fluid the TKE budget
!!       permits.  Penetrating shortwave radiation is also applied at
!!      this point.
!!    4. If there is any unentrained fluid that was formerly in the
!!      mixed layer, detrain this fluid into the buffer layer.  This
!!     is equivalent to the mixed layer detraining to the Monin-
!!       Obukhov depth.
!!    5. Divide the fluid in the mixed layer evenly into CS%nkml pieces.
!!    6. Split the buffer layer if appropriate.
!! Layers 1 to nkml are the mixed layer, nkml+1 to nkml+nkbl are the
!! buffer layers. The results of this subroutine are mathematically
!! identical if there are multiple pieces of the mixed layer with
!! the same density or if there is just a single layer. There is no
!! stability limit on the time step.
!!
!!   The key parameters for the mixed layer are found in the control structure.
!! These include mstar, nstar, nstar2, pen_SW_frac, pen_SW_scale, and TKE_decay.
!!   For the Oberhuber (1993) mixed layer, the values of these are:
!!      pen_SW_frac = 0.42, pen_SW_scale = 15.0 m, mstar = 1.25,
!!      nstar = 1, TKE_decay = 2.5, conv_decay = 0.5
!!  TKE_decay is 1/kappa in eq. 28 of Oberhuber (1993), while conv_decay is 1/mu.
!!  Conv_decay has been eliminated in favor of the well-calibrated form for the
!!  efficiency of penetrating convection from Wang (2003).
!!    For a traditional Kraus-Turner mixed layer, the values are:
!!      pen_SW_frac = 0.0, pen_SW_scale = 0.0 m, mstar = 1.25,
!!      nstar = 0.4, TKE_decay = 0.0, conv_decay = 0.0
subroutine bulkmixedlayer(h_3d, u_3d, v_3d, tv, fluxes, dt, ea, eb, G, GV, US, CS, &
                          optics, Hml, aggregate_FW_forcing, dt_diag, last_call)
  type(ocean_grid_type),      intent(inout) :: G      !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h_3d   !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in)    :: u_3d   !< Zonal velocities interpolated to h points
                                                      !! [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in)    :: v_3d   !< Zonal velocities interpolated to h points
                                                      !! [L T-1 ~> m s-1].
  type(thermo_var_ptrs),      intent(inout) :: tv     !< A structure containing pointers to any
                                                      !! available thermodynamic fields. Absent
                                                      !! fields have NULL ptrs.
  type(forcing),              intent(inout) :: fluxes !< A structure containing pointers to any
                                                      !! possible forcing fields.  Unused fields
                                                      !! have NULL ptrs.
  real,                       intent(in)    :: dt     !< Time increment [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: ea     !< The amount of fluid moved downward into a
                                                      !! layer; this should be increased due to
                                                      !! mixed layer detrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: eb     !< The amount of fluid moved upward into a
                                                      !! layer; this should be increased due to
                                                      !! mixed layer entrainment [H ~> m or kg m-2].
  type(bulkmixedlayer_CS),    pointer       :: CS     !< The control structure returned by a
                                                      !! previous call to mixedlayer_init.
  type(optics_type),          pointer       :: optics !< The structure containing the inverse of the
                                                      !! vertical absorption decay scale for
                                                      !! penetrating shortwave radiation [m-1].
  real, dimension(:,:),       pointer       :: Hml    !< Active mixed layer depth [m].
  logical,                    intent(in)    :: aggregate_FW_forcing !< If true, the net incoming and
                                                     !! outgoing surface freshwater fluxes are
                                                     !! combined before being applied, instead of
                                                     !! being applied separately.
  real,             optional, intent(in)    :: dt_diag  !< The diagnostic time step,
                                                      !! which may be less than dt if there are
                                                      !! two callse to mixedlayer [T ~> s].
  logical,          optional, intent(in)    :: last_call !< if true, this is the last call
                                                      !! to mixedlayer in the current time step, so
                                                      !! diagnostics will be written. The default is
                                                      !! .true.

  ! Local variables
  real, dimension(SZI_(G),SZK_(GV)) :: &
    eaml, &     !   The amount of fluid moved downward into a layer due to mixed
                ! layer detrainment [H ~> m or kg m-2]. (I.e. entrainment from above.)
    ebml        !   The amount of fluid moved upward into a layer due to mixed
                ! layer detrainment [H ~> m or kg m-2]. (I.e. entrainment from below.)

  ! If there is resorting, the vertical coordinate for these variables is the
  ! new, sorted index space.  Here layer 0 is an initially massless layer that
  ! will be used to hold the new mixed layer properties.
  real, dimension(SZI_(G),SZK0_(GV)) :: &
    h, &        !   The layer thickness [H ~> m or kg m-2].
    T, &        !   The layer temperatures [degC].
    S, &        !   The layer salinities [ppt].
    R0, &       !   The potential density referenced to the surface [R ~> kg m-3].
    Rcv         !   The coordinate variable potential density [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)) :: &
    u, &        !   The zonal velocity [L T-1 ~> m s-1].
    v, &        !   The meridional velocity [L T-1 ~> m s-1].
    h_orig, &   !   The original thickness [H ~> m or kg m-2].
    d_eb, &     !   The downward increase across a layer in the entrainment from
                ! below [H ~> m or kg m-2].  The sign convention is that positive values of
                ! d_eb correspond to a gain in mass by a layer by upward motion.
    d_ea, &     !   The upward increase across a layer in the entrainment from
                ! above [H ~> m or kg m-2].  The sign convention is that positive values of
                ! d_ea mean a net gain in mass by a layer from downward motion.
    eps         ! The (small) thickness that must remain in a layer [H ~> m or kg m-2].
  integer, dimension(SZI_(G),SZK_(GV)) :: &
    ksort       !   The sorted k-index that each original layer goes to.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    h_miss      !   The summed absolute mismatch [Z ~> m].
  real, dimension(SZI_(G)) :: &
    TKE, &      !   The turbulent kinetic energy available for mixing over a
                ! time step [Z L2 T-2 ~> m3 s-2].
    Conv_En, &  !   The turbulent kinetic energy source due to mixing down to
                ! the depth of free convection [Z L2 T-2 ~> m3 s-2].
    htot, &     !   The total depth of the layers being considered for
                ! entrainment [H ~> m or kg m-2].
    R0_tot, &   !   The integrated potential density referenced to the surface
                ! of the layers which are fully entrained [H R ~> kg m-2 or kg2 m-5].
    Rcv_tot, &  !   The integrated coordinate value potential density of the
                ! layers that are fully entrained [H R ~> kg m-2 or kg2 m-5].
    Ttot, &     !   The integrated temperature of layers which are fully
                ! entrained [degC H ~> degC m or degC kg m-2].
    Stot, &     !   The integrated salt of layers which are fully entrained
                ! [H ppt ~> m ppt or ppt kg m-2].
    uhtot, &    !   The depth integrated zonal and meridional velocities in the
    vhtot, &    ! mixed layer [H L T-1 ~> m2 s-1 or kg m-1 s-1].

    netMassInOut, &  ! The net mass flux (if non-Boussinsq) or volume flux (if
                     ! Boussinesq - i.e. the fresh water flux (P+R-E)) into the
                     ! ocean over a time step [H ~> m or kg m-2].
    NetMassOut,   &  ! The mass flux (if non-Boussinesq) or volume flux (if Boussinesq)
                     ! over a time step from evaporating fresh water [H ~> m or kg m-2]
    Net_heat, & !   The net heating at the surface over a time step [degC H ~> degC m or degC kg m-2].
                ! Any penetrating shortwave radiation is not included in Net_heat.
    Net_salt, & ! The surface salt flux into the ocean over a time step, ppt H.
    Idecay_len_TKE, &  ! The inverse of a turbulence decay length scale [H-1 ~> m-1 or m2 kg-1].
    p_ref, &    !   Reference pressure for the potential density governing mixed
                ! layer dynamics, almost always 0 (or 1e5) Pa.
    p_ref_cv, & !   Reference pressure for the potential density which defines
                ! the coordinate variable, set to P_Ref [Pa].
    dR0_dT, &   !   Partial derivative of the mixed layer potential density with
                ! temperature [R degC-1 ~> kg m-3 degC-1].
    dRcv_dT, &  !   Partial derivative of the coordinate variable potential
                ! density in the mixed layer with temperature [R degC-1 ~> kg m-3 degC-1].
    dR0_dS, &   !   Partial derivative of the mixed layer potential density with
                ! salinity [R ppt-1 ~> kg m-3 ppt-1].
    dRcv_dS, &  !   Partial derivative of the coordinate variable potential
                ! density in the mixed layer with salinity [R ppt-1 ~> kg m-3 ppt-1].
    TKE_river   ! The source of turbulent kinetic energy available for mixing
                ! at rivermouths [Z L2 T-3 ~> m3 s-3].

  real, dimension(max(CS%nsw,1),SZI_(G)) :: &
    Pen_SW_bnd  !   The penetrating fraction of the shortwave heating integrated
                ! over a time step in each band [degC H ~> degC m or degC kg m-2].
  real, dimension(max(CS%nsw,1),SZI_(G),SZK_(GV)) :: &
    opacity_band ! The opacity in each band [H-1 ~> m-1 or m2 kg-1]. The indicies are band, i, k.

  real :: cMKE(2,SZI_(G)) ! Coefficients of HpE and HpE^2 used in calculating the
                          ! denominator of MKE_rate; the two elements have differing
                          ! units of [H-1 ~> m-1 or m2 kg-1] and [H-2 ~> m-2 or m4 kg-2].
  real :: Irho0         ! 1.0 / rho_0 [R-1 ~> m3 kg-1]
  real :: Inkml, Inkmlm1!  1.0 / REAL(nkml) and  1.0 / REAL(nkml-1)
  real :: Ih            !   The inverse of a thickness [H-1 ~> m-1 or m2 kg-1].
  real :: Idt_diag      !   The inverse of the timestep used for diagnostics [T-1 ~> s-1].
  real :: RmixConst

  real, dimension(SZI_(G)) :: &
    dKE_FC, &   !   The change in mean kinetic energy due to free convection
                ! [Z L2 T-2 ~> m3 s-2].
    h_CA        !   The depth to which convective adjustment has gone [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)) :: &
    dKE_CA, &   !   The change in mean kinetic energy due to convective
                ! adjustment [Z L2 T-2 ~> m3 s-2].
    cTKE        !   The turbulent kinetic energy source due to convective
                ! adjustment [Z L2 T-2 ~> m3 s-2].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    Hsfc_max, & ! The thickness of the surface region (mixed and buffer layers)
                ! after entrainment but before any buffer layer detrainment [Z ~> m].
    Hsfc_used, & ! The thickness of the surface region after buffer layer
                ! detrainment [Z ~> m].
    Hsfc_min, & ! The minimum thickness of the surface region based on the
                ! new mixed layer depth and the previous thickness of the
                ! neighboring water columns [Z ~> m].
    h_sum, &    ! The total thickness of the water column [H ~> m or kg m-2].
    hmbl_prev   ! The previous thickness of the mixed and buffer layers [H ~> m or kg m-2].
  real, dimension(SZI_(G)) :: &
    Hsfc, &     !   The thickness of the surface region (mixed and buffer
                ! layers before detrainment in to the interior [H ~> m or kg m-2].
    max_BL_det  !   If non-negative, the maximum amount of entrainment from
                ! the buffer layers that will be allowed this time step [H ~> m or kg m-2].
  real :: dHsfc, dHD ! Local copies of nondimensional parameters.
  real :: H_nbr ! A minimum thickness based on neighboring thicknesses [H ~> m or kg m-2].

  real :: absf_x_H  ! The absolute value of f times the mixed layer thickness [Z T-1 ~> m s-1].
  real :: kU_star   ! Ustar times the Von Karmen constant [Z T-1 ~> m s-1].
  real :: dt__diag  ! A recaled copy of dt_diag (if present) or dt [T ~> s].
  logical :: write_diags  ! If true, write out diagnostics with this step.
  logical :: reset_diags  ! If true, zero out the accumulated diagnostics.
  integer :: i, j, k, is, ie, js, je, nz, nkmb, n
  integer :: nsw    ! The number of bands of penetrating shortwave radiation.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_mixed_layer: "//&
         "Module must be initialized before it is used.")
  if (GV%nkml < 1) return

  if (.not. associated(tv%eqn_of_state)) call MOM_error(FATAL, &
      "MOM_mixed_layer: Temperature, salinity and an equation of state "//&
      "must now be used.")
  if (.NOT. associated(fluxes%ustar)) call MOM_error(FATAL, &
      "MOM_mixed_layer: No surface TKE fluxes (ustar) defined in mixedlayer!")

  nkmb = CS%nkml+CS%nkbl
  Inkml = 1.0 / REAL(CS%nkml)
  if (CS%nkml > 1) Inkmlm1 = 1.0 / REAL(CS%nkml-1)

  Irho0 = 1.0 / (GV%Rho0)
  dt__diag = dt ; if (present(dt_diag)) dt__diag = dt_diag
  Idt_diag = 1.0 / (dt__diag)
  write_diags = .true. ; if (present(last_call)) write_diags = last_call

  p_ref(:) = 0.0 ; p_ref_cv(:) = tv%P_Ref

  nsw = CS%nsw

  if (CS%limit_det .or. (CS%id_Hsfc_min > 0)) then
    !$OMP parallel do default(shared)
    do j=js-1,je+1 ; do i=is-1,ie+1
      h_sum(i,j) = 0.0 ; hmbl_prev(i,j) = 0.0
    enddo ; enddo
    !$OMP parallel do default(shared)
    do j=js-1,je+1
      do k=1,nkmb ; do i=is-1,ie+1
        h_sum(i,j) = h_sum(i,j) + h_3d(i,j,k)
        hmbl_prev(i,j) = hmbl_prev(i,j) + h_3d(i,j,k)
      enddo ; enddo
      do k=nkmb+1,nz ; do i=is-1,ie+1
        h_sum(i,j) = h_sum(i,j) + h_3d(i,j,k)
      enddo ; enddo
    enddo

    call cpu_clock_begin(id_clock_pass)
    call create_group_pass(CS%pass_h_sum_hmbl_prev, h_sum,G%Domain)
    call create_group_pass(CS%pass_h_sum_hmbl_prev, hmbl_prev,G%Domain)
    call do_group_pass(CS%pass_h_sum_hmbl_prev, G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  ! Determine whether to zero out diagnostics before accumulation.
  reset_diags = .true.
  if (present(dt_diag) .and. write_diags .and. (dt__diag > dt)) &
    reset_diags = .false.  ! This is the second call to mixedlayer.

  if (reset_diags) then
    if (CS%TKE_diagnostics) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        CS%diag_TKE_wind(i,j) = 0.0 ; CS%diag_TKE_RiBulk(i,j) = 0.0
        CS%diag_TKE_conv(i,j) = 0.0 ; CS%diag_TKE_pen_SW(i,j) = 0.0
        CS%diag_TKE_mixing(i,j) = 0.0 ; CS%diag_TKE_mech_decay(i,j) = 0.0
        CS%diag_TKE_conv_decay(i,j) = 0.0 ; CS%diag_TKE_conv_s2(i,j) = 0.0
      enddo ; enddo
    endif
    if (allocated(CS%diag_PE_detrain)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        CS%diag_PE_detrain(i,j) = 0.0
      enddo ; enddo
    endif
    if (allocated(CS%diag_PE_detrain2)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        CS%diag_PE_detrain2(i,j) = 0.0
      enddo ; enddo
    endif
  endif

  if (CS%ML_resort) then
    do i=is,ie ; h_CA(i) = 0.0 ; enddo
    do k=1,nz ; do i=is,ie ; dKE_CA(i,k) = 0.0 ; cTKE(i,k) = 0.0 ; enddo ; enddo
  endif
  max_BL_det(:) = -1

  !$OMP parallel default(shared) firstprivate(dKE_CA,cTKE,h_CA,max_BL_det,p_ref,p_ref_cv) &
  !$OMP                 private(h,u,v,h_orig,eps,T,S,opacity_band,d_ea,d_eb,R0,Rcv,ksort, &
  !$OMP                         dR0_dT,dR0_dS,dRcv_dT,dRcv_dS,htot,Ttot,Stot,TKE,Conv_en, &
  !$OMP                         RmixConst,TKE_river,Pen_SW_bnd,netMassInOut,NetMassOut,   &
  !$OMP                         Net_heat,Net_salt,uhtot,vhtot,R0_tot,Rcv_tot,dKE_FC,      &
  !$OMP                         Idecay_len_TKE,cMKE,Hsfc,dHsfc,dHD,H_nbr,kU_Star,         &
  !$OMP                         absf_x_H,ebml,eaml)
  !$OMP do
  do j=js,je
    ! Copy the thicknesses and other fields to 2-d arrays.
    do k=1,nz ; do i=is,ie
      h(i,k) = h_3d(i,j,k) ; u(i,k) = u_3d(i,j,k) ; v(i,k) = v_3d(i,j,k)
      h_orig(i,k) = h_3d(i,j,k)
      eps(i,k) = 0.0 ; if (k > nkmb) eps(i,k) = GV%Angstrom_H
      T(i,k) = tv%T(i,j,k) ; S(i,k) = tv%S(i,j,k)
    enddo ; enddo
    if (nsw>0) call extract_optics_slice(optics, j, G, GV, opacity=opacity_band, opacity_scale=GV%H_to_m)

    do k=1,nz ; do i=is,ie
      d_ea(i,k) = 0.0 ; d_eb(i,k) = 0.0
    enddo ; enddo

    if (id_clock_EOS>0) call cpu_clock_begin(id_clock_EOS)
    ! Calculate an estimate of the mid-mixed layer pressure [Pa]
    do i=is,ie ; p_ref(i) = 0.0 ; enddo
    do k=1,CS%nkml ; do i=is,ie
      p_ref(i) = p_ref(i) + 0.5*GV%H_to_Pa*h(i,k)
    enddo ; enddo
    call calculate_density_derivs(T(:,1), S(:,1), p_ref, dR0_dT, dR0_dS, &
                                  is, ie-is+1, tv%eqn_of_state, scale=US%kg_m3_to_R)
    call calculate_density_derivs(T(:,1), S(:,1), p_ref_cv, dRcv_dT, dRcv_dS, &
                                  is, ie-is+1, tv%eqn_of_state, scale=US%kg_m3_to_R)
    do k=1,nz
      call calculate_density(T(:,k), S(:,k), p_ref, R0(:,k), is, ie-is+1, &
                             tv%eqn_of_state, scale=US%kg_m3_to_R)
      call calculate_density(T(:,k), S(:,k), p_ref_cv, Rcv(:,k), is, &
                             ie-is+1, tv%eqn_of_state, scale=US%kg_m3_to_R)
    enddo
    if (id_clock_EOS>0) call cpu_clock_end(id_clock_EOS)

    if (CS%ML_resort) then
      if (id_clock_resort>0) call cpu_clock_begin(id_clock_resort)
      if (CS%ML_presort_nz_conv_adj > 0) &
        call convective_adjustment(h(:,1:), u, v, R0(:,1:), Rcv(:,1:), T(:,1:), &
                                   S(:,1:), eps, d_eb, dKE_CA, cTKE, j, G, GV, US, CS, &
                                   CS%ML_presort_nz_conv_adj)

      call sort_ML(h(:,1:), R0(:,1:), eps, G, GV, CS, ksort)
      if (id_clock_resort>0) call cpu_clock_end(id_clock_resort)
    else
      do k=1,nz ; do i=is,ie ; ksort(i,k) = k ; enddo ; enddo

      if (id_clock_adjustment>0) call cpu_clock_begin(id_clock_adjustment)
      !  Undergo instantaneous entrainment into the buffer layers and mixed layers
      ! to remove hydrostatic instabilities.  Any water that is lighter than
      ! currently in the mixed or buffer layer is entrained.
      call convective_adjustment(h(:,1:), u, v, R0(:,1:), Rcv(:,1:), T(:,1:), &
                                 S(:,1:), eps, d_eb, dKE_CA, cTKE, j, G, GV, US, CS)
      do i=is,ie ; h_CA(i) = h(i,1) ; enddo

      if (id_clock_adjustment>0) call cpu_clock_end(id_clock_adjustment)
    endif

    if (associated(fluxes%lrunoff) .and. CS%do_rivermix) then

      ! Here we add an additional source of TKE to the mixed layer where river
      ! is present to simulate unresolved estuaries. The TKE input is diagnosed
      ! as follows:
      !   TKE_river[m3 s-3] = 0.5*rivermix_depth*g*Irho0*drho_ds*
      !                       River*(Samb - Sriver) = CS%mstar*U_star^3
      ! where River is in units of [m s-1].
      ! Samb = Ambient salinity at the mouth of the estuary
      ! rivermix_depth =  The prescribed depth over which to mix river inflow
      ! drho_ds = The gradient of density wrt salt at the ambient surface salinity.
      ! Sriver = 0 (i.e. rivers are assumed to be pure freshwater)
      RmixConst = 0.5*CS%rivermix_depth * GV%g_Earth * Irho0**2
      do i=is,ie
        TKE_river(i) = max(0.0, RmixConst*dR0_dS(i)* &
                      (fluxes%lrunoff(i,j) + fluxes%frunoff(i,j)) * S(i,1))
      enddo
    else
      do i=is,ie ; TKE_river(i) = 0.0 ; enddo
    endif


    if (id_clock_conv>0) call cpu_clock_begin(id_clock_conv)

    ! The surface forcing is contained in the fluxes type.
    ! We aggregate the thermodynamic forcing for a time step into the following:
    ! netMassInOut = water [H ~> m or kg m-2] added/removed via surface fluxes
    ! netMassOut   = water [H ~> m or kg m-2] removed via evaporating surface fluxes
    ! net_heat     = heat via surface fluxes [degC H ~> degC m or degC kg m-2]
    ! net_salt     = salt via surface fluxes [ppt H ~> dppt m or gSalt m-2]
    ! Pen_SW_bnd   = components to penetrative shortwave radiation
    call extractFluxes1d(G, GV, US, fluxes, optics, nsw, j, dt, &
                  CS%H_limit_fluxes, CS%use_river_heat_content, CS%use_calving_heat_content, &
                  h(:,1:), T(:,1:), netMassInOut, netMassOut, Net_heat, Net_salt, Pen_SW_bnd,&
                  tv, aggregate_FW_forcing)

    ! This subroutine causes the mixed layer to entrain to depth of free convection.
    call mixedlayer_convection(h(:,1:), d_eb, htot, Ttot, Stot, uhtot, vhtot, &
                               R0_tot, Rcv_tot, u, v, T(:,1:), S(:,1:),       &
                               R0(:,1:), Rcv(:,1:), eps,                      &
                               dR0_dT, dRcv_dT, dR0_dS, dRcv_dS,              &
                               netMassInOut, netMassOut, Net_heat, Net_salt,  &
                               nsw, Pen_SW_bnd, opacity_band, Conv_En,        &
                               dKE_FC, j, ksort, G, GV, US, CS, tv, fluxes, dt, &
                               aggregate_FW_forcing)

    if (id_clock_conv>0) call cpu_clock_end(id_clock_conv)

    !   Now the mixed layer undergoes mechanically forced entrainment.
    ! The mixed layer may entrain down to the Monin-Obukhov depth if the
    ! surface is becoming lighter, and is effecti1336vely detraining.

    !    First the TKE at the depth of free convection that is available
    !  to drive mixing is calculated.
    if (id_clock_mech>0) call cpu_clock_begin(id_clock_mech)

    call find_starting_TKE(htot, h_CA, fluxes, Conv_En, cTKE, dKE_FC, dKE_CA, &
                           TKE, TKE_river, Idecay_len_TKE, cMKE, dt, Idt_diag, &
                           j, ksort, G, GV, US, CS)

    ! Here the mechanically driven entrainment occurs.
    call mechanical_entrainment(h(:,1:), d_eb, htot, Ttot, Stot, uhtot, vhtot, &
                                R0_tot, Rcv_tot, u, v, T(:,1:), S(:,1:), R0(:,1:), Rcv(:,1:), eps, dR0_dT, dRcv_dT, &
                                cMKE, Idt_diag, nsw, Pen_SW_bnd, opacity_band, TKE, &
                                Idecay_len_TKE, j, ksort, G, GV, US, CS)

    call absorbRemainingSW(G, GV, US, h(:,1:), opacity_band, nsw, optics, j, dt, &
                           CS%H_limit_fluxes, CS%correct_absorption, CS%absorb_all_SW, &
                           T(:,1:), Pen_SW_bnd, eps, ksort, htot, Ttot)

    if (CS%TKE_diagnostics) then ; do i=is,ie
      CS%diag_TKE_mech_decay(i,j) = CS%diag_TKE_mech_decay(i,j) - Idt_diag * TKE(i)
    enddo ; endif
    if (id_clock_mech>0) call cpu_clock_end(id_clock_mech)

    ! Calculate the homogeneous mixed layer properties and store them in layer 0.
    do i=is,ie ; if (htot(i) > 0.0) then
      Ih = 1.0 / htot(i)
      R0(i,0) = R0_tot(i) * Ih ; Rcv(i,0) = Rcv_tot(i) * Ih
      T(i,0) = Ttot(i) * Ih ; S(i,0) = Stot(i) * Ih
      h(i,0) = htot(i)
    else ! This may not ever be needed?
      T(i,0) = T(i,1) ; S(i,0) = S(i,1) ; R0(i,0) = R0(i,1) ; Rcv(i,0) = Rcv(i,1)
      h(i,0) = htot(i)
    endif ; enddo
    if (write_diags .and. allocated(CS%ML_depth)) then ; do i=is,ie
      CS%ML_depth(i,j) = h(i,0) * GV%H_to_m  ! Rescale the diagnostic.
    enddo ; endif
    if (associated(Hml)) then ; do i=is,ie
      Hml(i,j) = G%mask2dT(i,j) * (h(i,0) * GV%H_to_m) ! Rescale the diagnostic for output.
    enddo ; endif

! At this point, return water to the original layers, but constrained to
! still be sorted.  After this point, all the water that is in massive
! interior layers will be denser than water remaining in the mixed- and
! buffer-layers.  To achieve this, some of these variable density layers
! might be split between two isopycnal layers that are denser than new
! mixed layer or any remaining water from the old mixed- or buffer-layers.
! Alternately, if there are fewer than nkbl of the old buffer or mixed layers
! with any mass, relatively light interior layers might be transferred to
! these unused layers (but not currently in the code).

    if (CS%ML_resort) then
      if (id_clock_resort>0) call cpu_clock_begin(id_clock_resort)
      call resort_ML(h(:,0:), T(:,0:), S(:,0:), R0(:,0:), Rcv(:,0:), GV%Rlay(:), eps, &
                     d_ea, d_eb, ksort, G, GV, CS, dR0_dT, dR0_dS, dRcv_dT, dRcv_dS)
      if (id_clock_resort>0) call cpu_clock_end(id_clock_resort)
    endif

    if (CS%limit_det .or. (CS%id_Hsfc_max > 0) .or. (CS%id_Hsfc_min > 0)) then
      do i=is,ie ; Hsfc(i) = h(i,0) ; enddo
      do k=1,nkmb ; do i=is,ie ; Hsfc(i) = Hsfc(i) + h(i,k) ; enddo ; enddo

      if (CS%limit_det .or. (CS%id_Hsfc_min > 0)) then
        dHsfc = CS%lim_det_dH_sfc ; dHD = CS%lim_det_dH_bathy
        do i=is,ie
          H_nbr = min(dHsfc*max(hmbl_prev(i-1,j), hmbl_prev(i+1,j), &
                                hmbl_prev(i,j-1), hmbl_prev(i,j+1)), &
                      max(hmbl_prev(i-1,j) - dHD*min(h_sum(i,j),h_sum(i-1,j)), &
                          hmbl_prev(i+1,j) - dHD*min(h_sum(i,j),h_sum(i+1,j)), &
                          hmbl_prev(i,j-1) - dHD*min(h_sum(i,j),h_sum(i,j-1)), &
                          hmbl_prev(i,j+1) - dHD*min(h_sum(i,j),h_sum(i,j+1))) )

          Hsfc_min(i,j) = GV%H_to_Z * max(h(i,0), min(Hsfc(i), H_nbr))

          if (CS%limit_det) max_BL_det(i) = max(0.0, Hsfc(i)-H_nbr)
        enddo
      endif

      if (CS%id_Hsfc_max > 0) then ; do i=is,ie
        Hsfc_max(i,j) = GV%H_to_Z * Hsfc(i)
      enddo ; endif
    endif

! Move water left in the former mixed layer into the buffer layer and
! from the buffer layer into the interior.  These steps might best be
! treated in conjuction.
    if (id_clock_detrain>0) call cpu_clock_begin(id_clock_detrain)
    if (CS%nkbl == 1) then
      call mixedlayer_detrain_1(h(:,0:), T(:,0:), S(:,0:), R0(:,0:), Rcv(:,0:), &
                                GV%Rlay(:), dt, dt__diag, d_ea, d_eb, j, G, GV, US, CS, &
                                dRcv_dT, dRcv_dS, max_BL_det)
    elseif (CS%nkbl == 2) then
      call mixedlayer_detrain_2(h(:,0:), T(:,0:), S(:,0:), R0(:,0:), Rcv(:,0:), &
                                GV%Rlay(:), dt, dt__diag, d_ea, j, G, GV, US, CS, &
                                dR0_dT, dR0_dS, dRcv_dT, dRcv_dS, max_BL_det)
    else ! CS%nkbl not = 1 or 2
      ! This code only works with 1 or 2 buffer layers.
      call MOM_error(FATAL, "MOM_mixed_layer: CS%nkbl must be 1 or 2 for now.")
    endif
    if (id_clock_detrain>0) call cpu_clock_end(id_clock_detrain)


    if (CS%id_Hsfc_used > 0) then
      do i=is,ie ; Hsfc_used(i,j) = GV%H_to_Z * h(i,0) ; enddo
      do k=CS%nkml+1,nkmb ; do i=is,ie
        Hsfc_used(i,j) = Hsfc_used(i,j) + GV%H_to_Z * h(i,k)
      enddo ; enddo
    endif

! Now set the properties of the layers in the mixed layer in the original
! 3-d variables.
    if (CS%Resolve_Ekman .and. (CS%nkml>1)) then
      ! The thickness of the topmost piece of the mixed layer is given by
      ! h_1 = H / (3 + sqrt(|f|*H^2/2*nu_max)), which asymptotes to the Ekman
      ! layer depth and 1/3 of the mixed layer depth.  This curve has been
      ! determined to maximize the impact of the Ekman transport in the mixed
      ! layer TKE budget with nkml=2.  With nkml=3, this should also be used,
      ! as the third piece will then optimally describe mixed layer
      ! restratification.  For nkml>=4 the whole strategy should be revisited.
      do i=is,ie
        kU_star = 0.41*fluxes%ustar(i,j) ! Maybe could be replaced with u*+w*?
        if (associated(fluxes%ustar_shelf) .and. &
            associated(fluxes%frac_shelf_h)) then
          if (fluxes%frac_shelf_h(i,j) > 0.0) &
            kU_star = (1.0 - fluxes%frac_shelf_h(i,j)) * kU_star + &
                      fluxes%frac_shelf_h(i,j) * (0.41*fluxes%ustar_shelf(i,j))
        endif
        absf_x_H = 0.25 * GV%H_to_Z * h(i,0) * &
            ((abs(G%CoriolisBu(I,J)) + abs(G%CoriolisBu(I-1,J-1))) + &
             (abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I-1,J))))
        ! If the mixed layer vertical viscosity specification is changed in
        ! MOM_vert_friction.F90, this line will have to be modified accordingly.
        h_3d(i,j,1) = h(i,0) / (3.0 + sqrt(absf_x_H*(absf_x_H + 2.0*kU_star) / kU_star**2))
        do k=2,CS%nkml
          ! The other layers are evenly distributed through the mixed layer.
          h_3d(i,j,k) = (h(i,0)-h_3d(i,j,1)) * Inkmlm1
          d_ea(i,k) = d_ea(i,k) + h_3d(i,j,k)
          d_ea(i,1) = d_ea(i,1) - h_3d(i,j,k)
        enddo
      enddo
    else
      do i=is,ie
        h_3d(i,j,1) = h(i,0) * Inkml
      enddo
      do k=2,CS%nkml ; do i=is,ie
        h_3d(i,j,k) = h(i,0) * Inkml
        d_ea(i,k) = d_ea(i,k) + h_3d(i,j,k)
        d_ea(i,1) = d_ea(i,1) - h_3d(i,j,k)
      enddo ; enddo
    endif
    do i=is,ie ; h(i,0) = 0.0 ; enddo
    do k=1,CS%nkml ; do i=is,ie
      tv%T(i,j,k) = T(i,0) ; tv%S(i,j,k) = S(i,0)
    enddo ; enddo

    ! These sum needs to be done in the original layer space.

    ! The treatment of layer 1 is atypical because evaporation shows up as
    ! negative ea(i,1), and because all precipitation goes straight into layer 1.
    ! The code is ordered so that any roundoff errors in ea are lost the surface.
!    do i=is,ie ; eaml(i,1) = 0.0 ; enddo
!    do k=2,nz ; do i=is,ie ; eaml(i,k) = eaml(i,k-1) - d_ea(i,k-1) ; enddo ; enddo
!    do i=is,ie ; eaml(i,1) = netMassInOut(i) ; enddo


    do i=is,ie
! eaml(i,nz) is derived from  h(i,nz) - h_orig(i,nz) = eaml(i,nz) - ebml(i,nz-1)
      ebml(i,nz) = 0.0
      eaml(i,nz) = (h(i,nz) - h_orig(i,nz)) - d_eb(i,nz)
    enddo
    do k=nz-1,1,-1 ; do i=is,ie
      ebml(i,k) = ebml(i,k+1) - d_eb(i,k+1)
      eaml(i,k) = eaml(i,k+1) + d_ea(i,k)
    enddo ; enddo
    do i=is,ie ; eaml(i,1) = netMassInOut(i) ; enddo

    ! Copy the interior thicknesses and other fields back to the 3-d arrays.
    do k=CS%nkml+1,nz ; do i=is,ie
      h_3d(i,j,k) = h(i,k); tv%T(i,j,k) = T(i,k) ; tv%S(i,j,k) = S(i,k)
    enddo ; enddo

    do k=1,nz ; do i=is,ie
      ea(i,j,k) = ea(i,j,k) + eaml(i,k)
      eb(i,j,k) = eb(i,j,k) + ebml(i,k)
    enddo ; enddo

    if (CS%id_h_mismatch > 0) then
      do i=is,ie
        h_miss(i,j) = GV%H_to_Z * abs(h_3d(i,j,1) - (h_orig(i,1) + &
          (eaml(i,1) + (ebml(i,1) - eaml(i,1+1)))))
      enddo
      do k=2,nz-1 ; do i=is,ie
        h_miss(i,j) = h_miss(i,j) + GV%H_to_Z * abs(h_3d(i,j,k) - (h_orig(i,k) + &
          ((eaml(i,k) - ebml(i,k-1)) + (ebml(i,k) - eaml(i,k+1)))))
      enddo ; enddo
      do i=is,ie
        h_miss(i,j) = h_miss(i,j) + GV%H_to_Z * abs(h_3d(i,j,nz) - (h_orig(i,nz) + &
          ((eaml(i,nz) - ebml(i,nz-1)) + ebml(i,nz))))
      enddo
    endif

  enddo ! j loop
  !$OMP end parallel

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  ! This needs to happen after the H update and before the next post_data.
  call diag_update_remap_grids(CS%diag)


  if (write_diags) then
    if (CS%id_ML_depth > 0) &
      call post_data(CS%id_ML_depth, CS%ML_depth, CS%diag)
    if (CS%id_TKE_wind > 0) &
      call post_data(CS%id_TKE_wind, CS%diag_TKE_wind, CS%diag)
    if (CS%id_TKE_RiBulk > 0) &
      call post_data(CS%id_TKE_RiBulk, CS%diag_TKE_RiBulk, CS%diag)
    if (CS%id_TKE_conv > 0) &
      call post_data(CS%id_TKE_conv, CS%diag_TKE_conv, CS%diag)
    if (CS%id_TKE_pen_SW > 0) &
      call post_data(CS%id_TKE_pen_SW, CS%diag_TKE_pen_SW, CS%diag)
    if (CS%id_TKE_mixing > 0) &
      call post_data(CS%id_TKE_mixing, CS%diag_TKE_mixing, CS%diag)
    if (CS%id_TKE_mech_decay > 0) &
      call post_data(CS%id_TKE_mech_decay, CS%diag_TKE_mech_decay, CS%diag)
    if (CS%id_TKE_conv_decay > 0) &
      call post_data(CS%id_TKE_conv_decay, CS%diag_TKE_conv_decay, CS%diag)
    if (CS%id_TKE_conv_s2 > 0) &
      call post_data(CS%id_TKE_conv_s2, CS%diag_TKE_conv_s2, CS%diag)
    if (CS%id_PE_detrain > 0) &
      call post_data(CS%id_PE_detrain, CS%diag_PE_detrain, CS%diag)
    if (CS%id_PE_detrain2 > 0) &
      call post_data(CS%id_PE_detrain2, CS%diag_PE_detrain2, CS%diag)
    if (CS%id_h_mismatch > 0) &
      call post_data(CS%id_h_mismatch, h_miss, CS%diag)
    if (CS%id_Hsfc_used > 0) &
      call post_data(CS%id_Hsfc_used, Hsfc_used, CS%diag)
    if (CS%id_Hsfc_max > 0) &
      call post_data(CS%id_Hsfc_max, Hsfc_max, CS%diag)
    if (CS%id_Hsfc_min > 0) &
      call post_data(CS%id_Hsfc_min, Hsfc_min, CS%diag)
  endif

end subroutine bulkmixedlayer

!>   This subroutine does instantaneous convective entrainment into the buffer
!! layers and mixed layers to remove hydrostatic instabilities.  Any water that
!! is lighter than currently in the mixed- or buffer- layer is entrained.
subroutine convective_adjustment(h, u, v, R0, Rcv, T, S, eps, d_eb, &
                                 dKE_CA, cTKE, j, G, GV, US, CS, nz_conv)
  type(ocean_grid_type),             intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),           intent(in)    :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZK_(GV)), intent(inout) :: h    !< Layer thickness [H ~> m or kg m-2].
                                                           !! The units of h are referred to as H below.
  real, dimension(SZI_(G),SZK_(GV)), intent(inout) :: u    !< Zonal velocities interpolated to h
                                                           !! points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK_(GV)), intent(inout) :: v    !< Zonal velocities interpolated to h
                                                           !! points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK_(GV)), intent(inout) :: T    !< Layer temperatures [degC].
  real, dimension(SZI_(G),SZK_(GV)), intent(inout) :: S    !< Layer salinities [ppt].
  real, dimension(SZI_(G),SZK_(GV)), intent(inout) :: R0   !< Potential density referenced to
                                                           !! surface pressure [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), intent(inout) :: Rcv  !< The coordinate defining potential
                                                           !! density [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), intent(inout) :: d_eb !< The downward increase across a layer
                                                           !! in the entrainment from below [H ~> m or kg m-2].
                                                           !! Positive values go with mass gain by
                                                           !! a layer.
  real, dimension(SZI_(G),SZK_(GV)), intent(in)    :: eps  !< The negligibly small amount of water
                                                           !! that will be left in each layer [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)), intent(out)   :: dKE_CA !< The vertically integrated change in
                                                           !! kinetic energy due to convective
                                                           !! adjustment [Z L2 T-2 ~> m3 s-2].
  real, dimension(SZI_(G),SZK_(GV)), intent(out)   :: cTKE !< The buoyant turbulent kinetic energy
                                                           !! source due to convective adjustment
                                                           !! [Z L2 T-2 ~> m3 s-2].
  integer,                           intent(in)    :: j    !< The j-index to work on.
  type(unit_scale_type),             intent(in)    :: US   !< A dimensional unit scaling type
  type(bulkmixedlayer_CS),           pointer       :: CS   !< The control structure for this module.
  integer,                 optional, intent(in)    :: nz_conv !< If present, the number of layers
                                                           !! over which to do convective adjustment
                                                           !! (perhaps CS%nkml).

!   This subroutine does instantaneous convective entrainment into the buffer
! layers and mixed layers to remove hydrostatic instabilities.  Any water that
! is lighter than currently in the mixed- or buffer- layer is entrained.

  ! Local variables
  real, dimension(SZI_(G)) :: &
    htot, &     !   The total depth of the layers being considered for
                ! entrainment [H ~> m or kg m-2].
    R0_tot, &   !   The integrated potential density referenced to the surface
                ! of the layers which are fully entrained [H R ~> kg m-2 or kg2 m-5].
    Rcv_tot, &  !   The integrated coordinate value potential density of the
                ! layers that are fully entrained [H R ~> kg m-2 or kg2 m-5].
    Ttot, &     !   The integrated temperature of layers which are fully
                ! entrained [degC H ~> degC m or degC kg m-2].
    Stot, &     !   The integrated salt of layers which are fully entrained
                ! [H ppt ~> m ppt or ppt kg m-2].
    uhtot, &    !   The depth integrated zonal and meridional velocities in
    vhtot, &    ! the mixed layer [H L T-1 ~> m2 s-1 or kg m-1 s-1].
    KE_orig, &  !   The total mean kinetic energy in the mixed layer before
                ! convection, [H L2 T-2 ~> H m2 s-2].
    h_orig_k1   !   The depth of layer k1 before convective adjustment [H ~> m or kg m-2].
  real :: h_ent !   The thickness from a layer that is entrained [H ~> m or kg m-2].
  real :: Ih    !   The inverse of a thickness [H-1 ~> m-1 or m2 kg-1].
  real :: g_H2_2Rho0  !   Half the gravitational acceleration times the square of
                      ! the conversion from H to Z divided by the mean density,
                      ! in [L2 Z T-3 H-2 R-1 ~> m4 s-3 kg-1 or m10 s-3 kg-3].
  integer :: is, ie, nz, i, k, k1, nzc, nkmb

  is = G%isc ; ie = G%iec ; nz = GV%ke
  g_H2_2Rho0 = (GV%g_Earth * GV%H_to_Z**2) / (2.0 * GV%Rho0)
  nzc = nz ; if (present(nz_conv)) nzc = nz_conv
  nkmb = CS%nkml+CS%nkbl

!   Undergo instantaneous entrainment into the buffer layers and mixed layers
! to remove hydrostatic instabilities.  Any water that is lighter than currently
! in the layer is entrained.
  do k1=min(nzc-1,nkmb),1,-1
    do i=is,ie
      h_orig_k1(i) = h(i,k1)
      KE_orig(i) = 0.5*h(i,k1)*(u(i,k1)**2 + v(i,k1)**2)
      uhtot(i) = h(i,k1)*u(i,k1) ; vhtot(i) = h(i,k1)*v(i,k1)
      R0_tot(i) = R0(i,k1) * h(i,k1)
      cTKE(i,k1) = 0.0 ; dKE_CA(i,k1) = 0.0

      Rcv_tot(i) = Rcv(i,k1) * h(i,k1)
      Ttot(i) = T(i,k1) * h(i,k1) ; Stot(i) = S(i,k1) * h(i,k1)
    enddo
    do k=k1+1,nzc
      do i=is,ie
        if ((h(i,k) > eps(i,k)) .and. (R0_tot(i) > h(i,k1)*R0(i,k))) then
          h_ent = h(i,k)-eps(i,k)
          cTKE(i,k1) = cTKE(i,k1) + h_ent * g_H2_2Rho0 * &
                   (R0_tot(i) - h(i,k1)*R0(i,k)) * CS%nstar2
          if (k < nkmb) then
            cTKE(i,k1) = cTKE(i,k1) + cTKE(i,k)
            dKE_CA(i,k1) = dKE_CA(i,k1) + dKE_CA(i,k)
          endif
          R0_tot(i) = R0_tot(i) + h_ent * R0(i,k)
          KE_orig(i) = KE_orig(i) + 0.5*h_ent* &
              (u(i,k)*u(i,k) + v(i,k)*v(i,k))
          uhtot(i) = uhtot(i) + h_ent*u(i,k)
          vhtot(i) = vhtot(i) + h_ent*v(i,k)

          Rcv_tot(i) = Rcv_tot(i) + h_ent * Rcv(i,k)
          Ttot(i) = Ttot(i) + h_ent * T(i,k)
          Stot(i) = Stot(i) + h_ent * S(i,k)
          h(i,k1) = h(i,k1) + h_ent ; h(i,k) = eps(i,k)

          d_eb(i,k) = d_eb(i,k) - h_ent
          d_eb(i,k1) = d_eb(i,k1) + h_ent
        endif
      enddo
    enddo
! Determine the temperature, salinity, and velocities of the mixed or buffer
! layer in question, if it has entrained.
    do i=is,ie ; if (h(i,k1) > h_orig_k1(i)) then
      Ih = 1.0 / h(i,k1)
      R0(i,k1) = R0_tot(i) * Ih
      u(i,k1) = uhtot(i) * Ih ; v(i,k1) = vhtot(i) * Ih
      dKE_CA(i,k1) = dKE_CA(i,k1) + GV%H_to_Z * (CS%bulk_Ri_convective * &
           (KE_orig(i) - 0.5*h(i,k1)*(u(i,k1)**2 + v(i,k1)**2)))
      Rcv(i,k1) = Rcv_tot(i) * Ih
      T(i,k1) = Ttot(i) * Ih ; S(i,k1) = Stot(i) * Ih
    endif ; enddo
  enddo
! If lower mixed or buffer layers are massless, give them the properties of the
! layer above.
  do k=2,min(nzc,nkmb) ; do i=is,ie ; if (h(i,k) == 0.0) then
    R0(i,k) = R0(i,k-1)
    Rcv(i,k) = Rcv(i,k-1) ; T(i,k) = T(i,k-1) ; S(i,k) = S(i,k-1)
  endif ; enddo ; enddo

end subroutine convective_adjustment

!>   This subroutine causes the mixed layer to entrain to the depth of free
!! convection.  The depth of free convection is the shallowest depth at which the
!! fluid is denser than the average of the fluid above.
subroutine mixedlayer_convection(h, d_eb, htot, Ttot, Stot, uhtot, vhtot,      &
                                 R0_tot, Rcv_tot, u, v, T, S, R0, Rcv, eps,    &
                                 dR0_dT, dRcv_dT, dR0_dS, dRcv_dS,             &
                                 netMassInOut, netMassOut, Net_heat, Net_salt, &
                                 nsw, Pen_SW_bnd, opacity_band, Conv_En,       &
                                 dKE_FC, j, ksort, G, GV, US, CS, tv, fluxes, dt,      &
                                 aggregate_FW_forcing)
  type(ocean_grid_type),    intent(in)    :: G     !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV    !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(inout) :: h     !< Layer thickness [H ~> m or kg m-2].
                                                   !! The units of h are referred to as H below.
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(inout) :: d_eb  !< The downward increase across a layer in the
                                                   !! layer in the entrainment from below [H ~> m or kg m-2].
                                                   !! Positive values go with mass gain by a layer.
  real, dimension(SZI_(G)), intent(out)   :: htot  !< The accumulated mixed layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G)), intent(out)   :: Ttot  !< The depth integrated mixed layer temperature
                                                   !! [degC H ~> degC m or degC kg m-2].
  real, dimension(SZI_(G)), intent(out)   :: Stot  !< The depth integrated mixed layer salinity
                                                   !! [ppt H ~> ppt m or ppt kg m-2].
  real, dimension(SZI_(G)), intent(out)   :: uhtot !< The depth integrated mixed layer zonal
                                                   !! velocity [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G)), intent(out)   :: vhtot !< The integrated mixed layer meridional
                                                   !! velocity [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G)), intent(out)   :: R0_tot !< The integrated mixed layer potential density referenced
                                                   !! to 0 pressure [H R ~> kg m-2 or kg2 m-5].
  real, dimension(SZI_(G)), intent(out)   :: Rcv_tot !< The integrated mixed layer coordinate
                                                   !! variable potential density [H R ~> kg m-2 or kg2 m-5].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: u     !< Zonal velocities interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: v     !< Zonal velocities interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: T     !< Layer temperatures [degC].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: S     !< Layer salinities [ppt].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: R0    !< Potential density referenced to
                                                   !! surface pressure [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: Rcv   !< The coordinate defining potential
                                                   !! density [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: eps   !< The negligibly small amount of water
                                                   !! that will be left in each layer [H ~> m or kg m-2].
  real, dimension(SZI_(G)), intent(in)    :: dR0_dT  !< The partial derivative of R0 with respect to
                                                   !! temperature [R degC-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)), intent(in)    :: dRcv_dT !< The partial derivative of Rcv with respect to
                                                   !! temperature [R degC-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)), intent(in)    :: dR0_dS  !< The partial derivative of R0 with respect to
                                                   !! salinity [R ppt-1 ~> kg m-3 ppt-1].
  real, dimension(SZI_(G)), intent(in)    :: dRcv_dS !< The partial derivative of Rcv with respect to
                                                   !! salinity [R ppt-1 ~> kg m-3 ppt-1].
  real, dimension(SZI_(G)), intent(in)    :: netMassInOut !< The net mass flux (if non-Boussinesq)
                                                   !! or volume flux (if Boussinesq) into the ocean
                                                   !! within a time step [H ~> m or kg m-2]. (I.e. P+R-E.)
  real, dimension(SZI_(G)), intent(in)    :: netMassOut !< The mass or volume flux out of the ocean
                                                   !! within a time step [H ~> m or kg m-2].
  real, dimension(SZI_(G)), intent(in)    :: Net_heat !< The net heating at the surface over a time
                                                   !! step [degC H ~> degC m or degC kg m-2].  Any penetrating
                                                   !! shortwave radiation is not included in Net_heat.
  real, dimension(SZI_(G)), intent(in)    :: Net_salt !< The net surface salt flux into the ocean
                                                   !! over a time step [ppt H ~> ppt m or ppt kg m-2].
  integer,                  intent(in)    :: nsw   !< The number of bands of penetrating
                                                   !! shortwave radiation.
  real, dimension(max(nsw,1),SZI_(G)), intent(inout) :: Pen_SW_bnd !< The penetrating shortwave
                                                   !! heating at the sea surface in each penetrating
                                                   !! band [degC H ~> degC m or degC kg m-2].
  real, dimension(max(nsw,1),SZI_(G),SZK_(GV)), intent(in) :: opacity_band !< The opacity in each band of
                                                   !! penetrating shortwave radiation [H-1 ~> m-1 or m2 kg-1].
  real, dimension(SZI_(G)), intent(out)   :: Conv_En !< The buoyant turbulent kinetic energy source
                                                   !! due to free convection [Z L2 T-2 ~> m3 s-2].
  real, dimension(SZI_(G)), intent(out)   :: dKE_FC !< The vertically integrated change in kinetic
                                                   !! energy due to free convection [Z L2 T-2 ~> m3 s-2].
  integer,                  intent(in)    :: j     !< The j-index to work on.
  integer, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: ksort !< The density-sorted k-indices.
  type(unit_scale_type),    intent(in)    :: US    !< A dimensional unit scaling type
  type(bulkmixedlayer_CS),  pointer       :: CS    !< The control structure for this module.
  type(thermo_var_ptrs),    intent(inout) :: tv    !< A structure containing pointers to any
                                                   !! available thermodynamic fields. Absent
                                                   !! fields have NULL ptrs.
  type(forcing),            intent(inout) :: fluxes  !< A structure containing pointers to any
                                                   !! possible forcing fields.  Unused fields
                                                   !! have NULL ptrs.
  real,                     intent(in)    :: dt    !< Time increment [T ~> s].
  logical,                  intent(in)    :: aggregate_FW_forcing !< If true, the net incoming and
                                                   !! outgoing surface freshwater fluxes are
                                                   !! combined before being applied, instead of
                                                   !! being applied separately.

!   This subroutine causes the mixed layer to entrain to the depth of free
! convection.  The depth of free convection is the shallowest depth at which the
! fluid is denser than the average of the fluid above.

  ! Local variables
  real, dimension(SZI_(G)) :: &
    massOutRem, &      !   Evaporation that remains to be supplied [H ~> m or kg m-2].
    netMassIn          ! mass entering through ocean surface [H ~> m or kg m-2]
  real :: SW_trans     !   The fraction of shortwave radiation
                       ! that is not absorbed in a layer [nondim].
  real :: Pen_absorbed !   The amount of penetrative shortwave radiation
                       ! that is absorbed in a layer [degC H ~> degC m or degC kg m-2].
  real :: h_avail      !   The thickness in a layer available for
                       ! entrainment [H ~> m or kg m-2].
  real :: h_ent        !   The thickness from a layer that is entrained [H ~> m or kg m-2].
  real :: T_precip     !   The temperature of the precipitation [degC].
  real :: C1_3, C1_6   !  1/3 and 1/6.
  real :: En_fn, Frac, x1 !  Nondimensional temporary variables.
  real :: dr, dr0      ! Temporary variables [R H ~> kg m-2 or kg2 m-5].
  real :: dr_ent, dr_comp ! Temporary variables [R H ~> kg m-2 or kg2 m-5].
  real :: dr_dh        ! The partial derivative of dr_ent with h_ent [R ~> kg m-3].
  real :: h_min, h_max !   The minimum, maximum, and previous estimates for
  real :: h_prev       ! h_ent [H ~> m or kg m-2].
  real :: h_evap       !   The thickness that is evaporated [H ~> m or kg m-2].
  real :: dh_Newt      !   The Newton's method estimate of the change in
                       ! h_ent between iterations [H ~> m or kg m-2].
  real :: g_H2_2Rho0   !   Half the gravitational acceleration times the square of
                       ! the conversion from H to Z divided by the mean density,
                       ! [L2 Z T-3 H-2 R-1 ~> m4 s-3 kg-1 or m10 s-3 kg-3].
  real :: Angstrom     !   The minimum layer thickness [H ~> m or kg m-2].
  real :: opacity      !   The opacity converted to inverse thickness units [H-1 ~> m-1 or m2 kg-1]
  real :: sum_Pen_En   !   The potential energy change due to penetrating
                       ! shortwave radiation, integrated over a layer
                       ! [H R ~> kg m-2 or kg2 m-5].
  real :: Idt          ! 1.0/dt [T-1 ~> s-1]
  real :: netHeatOut   ! accumulated heat content of mass leaving ocean
  integer :: is, ie, nz, i, k, ks, itt, n
  real, dimension(max(nsw,1)) :: &
    C2, &              ! Temporary variable R H-1 ~> kg m-4 or m-1].
    r_SW_top           ! Temporary variables [H R ~> kg m-2 or kg2 m-5].

  Angstrom = GV%Angstrom_H
  C1_3 = 1.0/3.0 ; C1_6 = 1.0/6.0
  g_H2_2Rho0 = (GV%g_Earth * GV%H_to_Z**2) / (2.0 * GV%Rho0)
  Idt = 1.0 / dt
  is = G%isc ; ie = G%iec ; nz = GV%ke

  do i=is,ie ; if (ksort(i,1) > 0) then
    k = ksort(i,1)

    if (aggregate_FW_forcing) then
      massOutRem(i) = 0.0
      if (netMassInOut(i) < 0.0) massOutRem(i) = -netMassInOut(i)
      netMassIn(i) = netMassInOut(i) + massOutRem(i)
    else
      massOutRem(i) = -netMassOut(i)
      netMassIn(i)  = netMassInOut(i) - netMassOut(i)
    endif

    ! htot is an Angstrom (taken from layer 1) plus any net precipitation.
    h_ent     = max(min(Angstrom,h(i,k)-eps(i,k)),0.0)
    htot(i)   = h_ent + netMassIn(i)
    h(i,k)    = h(i,k) - h_ent
    d_eb(i,k) = d_eb(i,k) - h_ent

    Pen_absorbed = 0.0
    do n=1,nsw ; if (Pen_SW_bnd(n,i) > 0.0) then
      SW_trans        = exp(-htot(i)*opacity_band(n,i,k))
      Pen_absorbed    = Pen_absorbed + Pen_SW_bnd(n,i) * (1.0-SW_trans)
      Pen_SW_bnd(n,i) = Pen_SW_bnd(n,i) * SW_trans
    endif ; enddo

    ! Precipitation is assumed to have the same temperature and velocity
    ! as layer 1.  Because layer 1 might not be the topmost layer, this
    ! involves multiple terms.
    T_precip = T(i,1)
    Ttot(i) = (Net_heat(i) + (netMassIn(i) * T_precip + h_ent * T(i,k))) + &
              Pen_absorbed
    ! Net_heat contains both heat fluxes and the heat content of mass fluxes.
 !! Ttot(i) = netMassIn(i) * T_precip + h_ent * T(i,k)
 !! Ttot(i) = Net_heat(i) + Ttot(i)
 !! Ttot(i) = Ttot(i) + Pen_absorbed
  ! smg:
  ! Ttot(i)   = (Net_heat(i) + (h_ent * T(i,k))) + Pen_absorbed
    Stot(i)   = h_ent*S(i,k) + Net_salt(i)
    uhtot(i)  = u(i,1)*netMassIn(i) + u(i,k)*h_ent
    vhtot(i)  = v(i,1)*netMassIn(i) + v(i,k)*h_ent
    R0_tot(i) = (h_ent*R0(i,k) + netMassIn(i)*R0(i,1)) + &
!                   dR0_dT(i)*netMassIn(i)*(T_precip - T(i,1)) + &
                (dR0_dT(i)*(Net_heat(i) + Pen_absorbed) - &
                 dR0_dS(i) * (netMassIn(i) * S(i,1) - Net_salt(i)))
    Rcv_tot(i) = (h_ent*Rcv(i,k) + netMassIn(i)*Rcv(i,1)) + &
!                    dRcv_dT(i)*netMassIn(i)*(T_precip - T(i,1)) + &
                 (dRcv_dT(i)*(Net_heat(i) + Pen_absorbed) - &
                  dRcv_dS(i) * (netMassIn(i) * S(i,1) - Net_salt(i)))
    Conv_En(i) = 0.0 ; dKE_FC(i) = 0.0
    if (associated(fluxes%heat_content_massin)) &
      fluxes%heat_content_massin(i,j) = fluxes%heat_content_massin(i,j) + &
                         T_precip * netMassIn(i) * GV%H_to_RZ * fluxes%C_p * Idt
    if (associated(tv%TempxPmE)) tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + &
                         T_precip * netMassIn(i) * GV%H_to_RZ
  else  ! This is a massless column, but zero out the summed variables anyway for safety.
    htot(i) = 0.0 ; Ttot(i) = 0.0 ; Stot(i) = 0.0 ; R0_tot(i) = 0.0 ; Rcv_tot = 0.0
    uhtot(i) = 0.0 ; vhtot(i) = 0.0 ; Conv_En(i) = 0.0 ; dKE_FC(i) = 0.0
  endif ; enddo

  ! Now do netMassOut case in this block.
  ! At this point htot contains an Angstrom of fluid from layer 0 plus netMassIn.
  do ks=1,nz
    do i=is,ie ; if (ksort(i,ks) > 0) then
      k = ksort(i,ks)

      if ((htot(i) < Angstrom) .and. (h(i,k) > eps(i,k))) then
        ! If less than an Angstrom was available from the layers above plus
        ! any precipitation, add more fluid from this layer.
        h_ent = min(Angstrom-htot(i), h(i,k)-eps(i,k))
        htot(i) = htot(i) + h_ent
        h(i,k) = h(i,k) - h_ent
        d_eb(i,k) = d_eb(i,k) - h_ent

        R0_tot(i) = R0_tot(i) + h_ent*R0(i,k)
        uhtot(i) = uhtot(i) + h_ent*u(i,k)
        vhtot(i) = vhtot(i) + h_ent*v(i,k)

        Rcv_tot(i) = Rcv_tot(i) + h_ent*Rcv(i,k)
        Ttot(i) = Ttot(i) + h_ent*T(i,k)
        Stot(i) = Stot(i) + h_ent*S(i,k)
      endif

      ! Water is removed from the topmost layers with any mass.
      ! We may lose layers if they are thin enough.
      ! The salt that is left behind goes into Stot.
      if ((massOutRem(i) > 0.0) .and. (h(i,k) > eps(i,k))) then
        if (massOutRem(i) > (h(i,k) - eps(i,k))) then
          h_evap = h(i,k) - eps(i,k)
          h(i,k) = eps(i,k)
          massOutRem(i) = massOutRem(i) - h_evap
        else
          h_evap = massOutRem(i)
          h(i,k) = h(i,k) - h_evap
          massOutRem(i) = 0.0
        endif

        Stot(i) = Stot(i) + h_evap*S(i,k)
        R0_tot(i) = R0_tot(i) + dR0_dS(i)*h_evap*S(i,k)
        Rcv_tot(i) = Rcv_tot(i) + dRcv_dS(i)*h_evap*S(i,k)
        d_eb(i,k) = d_eb(i,k) - h_evap

        ! smg: when resolve the A=B code, we will set
        ! heat_content_massout = heat_content_massout - T(i,k)*h_evap*GV%H_to_RZ*fluxes%C_p*Idt
        ! by uncommenting the lines here.
        ! we will also then completely remove TempXpme from the model.
        if (associated(fluxes%heat_content_massout)) &
          fluxes%heat_content_massout(i,j) = fluxes%heat_content_massout(i,j) - &
                                      T(i,k)*h_evap*GV%H_to_RZ * fluxes%C_p * Idt
        if (associated(tv%TempxPmE)) tv%TempxPmE(i,j) = tv%TempxPmE(i,j) - &
                                      T(i,k)*h_evap*GV%H_to_RZ

      endif

      ! The following section calculates how much fluid will be entrained.
      h_avail = h(i,k) - eps(i,k)
      if (h_avail > 0.0) then
        dr = R0_tot(i) - htot(i)*R0(i,k)
        h_ent = 0.0

        dr0 = dr
        do n=1,nsw ; if (Pen_SW_bnd(n,i) > 0.0) then
          dr0 = dr0 - (dR0_dT(i)*Pen_SW_bnd(n,i)) * &
                      opacity_band(n,i,k)*htot(i)
        endif ; enddo

        ! Some entrainment will occur from this layer.
        if (dr0 > 0.0) then
          dr_comp = dr
          do n=1,nsw ; if (Pen_SW_bnd(n,i) > 0.0) then
            !   Compare the density at the bottom of a layer with the
            ! density averaged over the mixed layer and that layer.
            opacity = opacity_band(n,i,k)
            SW_trans = exp(-h_avail*opacity)
            dr_comp = dr_comp + (dR0_dT(i)*Pen_SW_bnd(n,i)) * &
                ((1.0 - SW_trans) - opacity*(htot(i)+h_avail)*SW_trans)
          endif ; enddo
          if (dr_comp >= 0.0) then
            ! The entire layer is entrained.
            h_ent = h_avail
          else
            !  The layer is partially entrained.   Iterate to determine how much
            !  entrainment occurs.  Solve for the h_ent at which dr_ent = 0.

            ! Instead of assuming that the curve is linear between the two end
            ! points, assume that the change is concentrated near small values
            ! of entrainment.  On average, this saves about 1 iteration.
            Frac = dr0 / (dr0 - dr_comp)
            h_ent = h_avail * Frac*Frac
            h_min = 0.0 ; h_max = h_avail

            do n=1,nsw
              r_SW_top(n) = dR0_dT(i) * Pen_SW_bnd(n,i)
              C2(n) = r_SW_top(n) * opacity_band(n,i,k)**2
            enddo
            do itt=1,10
              dr_ent = dr ; dr_dh = 0.0
              do n=1,nsw
                opacity = opacity_band(n,i,k)
                SW_trans = exp(-h_ent*opacity)
                dr_ent = dr_ent + r_SW_top(n) * ((1.0 - SW_trans) - &
                           opacity*(htot(i)+h_ent)*SW_trans)
                dr_dh = dr_dh + C2(n) * (htot(i)+h_ent) * SW_trans
              enddo

              if (dr_ent > 0.0) then
                h_min = h_ent
              else
                h_max = h_ent
              endif

              dh_Newt = -dr_ent / dr_dh
              h_prev = h_ent ; h_ent = h_prev+dh_Newt
              if (h_ent > h_max) then
                h_ent = 0.5*(h_prev+h_max)
              elseif (h_ent < h_min) then
                h_ent = 0.5*(h_prev+h_min)
              endif

              if (ABS(dh_Newt) < 0.2*Angstrom) exit
            enddo

          endif

          !  Now that the amount of entrainment (h_ent) has been determined,
          !  calculate changes in various terms.
          sum_Pen_En = 0.0 ; Pen_absorbed = 0.0
          do n=1,nsw ; if (Pen_SW_bnd(n,i) > 0.0) then
            opacity = opacity_band(n,i,k)
            SW_trans = exp(-h_ent*opacity)

            x1 = h_ent*opacity
            if (x1 < 2.0e-5) then
              En_fn = (opacity*htot(i)*(1.0 - 0.5*(x1 - C1_3*x1)) + &
                       x1*x1*C1_6)
            else
              En_fn = ((opacity*htot(i) + 2.0) * &
                       ((1.0-SW_trans) / x1) - 1.0 + SW_trans)
            endif
            sum_Pen_En = sum_Pen_En - (dR0_dT(i)*Pen_SW_bnd(n,i)) * En_fn

            Pen_absorbed = Pen_absorbed + Pen_SW_bnd(n,i) * (1.0 - SW_trans)
            Pen_SW_bnd(n,i) = Pen_SW_bnd(n,i) * SW_trans
          endif ; enddo

          Conv_En(i) = Conv_En(i) + g_H2_2Rho0 * h_ent * &
                       ( (R0_tot(i) - R0(i,k)*htot(i)) + sum_Pen_En )

          R0_tot(i) = R0_tot(i) + (h_ent * R0(i,k) + Pen_absorbed*dR0_dT(i))
          Stot(i) = Stot(i) + h_ent * S(i,k)
          Ttot(i) = Ttot(i) + (h_ent * T(i,k) + Pen_absorbed)
          Rcv_tot(i) = Rcv_tot(i) + (h_ent * Rcv(i,k) + Pen_absorbed*dRcv_dT(i))
        endif ! dr0 > 0.0

        if (h_ent > 0.0) then
          if (htot(i) > 0.0) &
            dKE_FC(i) = dKE_FC(i) + CS%bulk_Ri_convective * 0.5 * &
              ((GV%H_to_Z*h_ent) / (htot(i)*(h_ent+htot(i)))) * &
              ((uhtot(i)-u(i,k)*htot(i))**2 + (vhtot(i)-v(i,k)*htot(i))**2)

          htot(i)  = htot(i)  + h_ent
          h(i,k) = h(i,k) - h_ent
          d_eb(i,k) = d_eb(i,k) - h_ent
          if (CS%convect_mom_bug) then
            uhtot(i) = u(i,k)*h_ent ; vhtot(i) = v(i,k)*h_ent
          else
            uhtot(i) = uhtot(i) + h_ent*u(i,k) ; vhtot(i) = vhtot(i) + h_ent*v(i,k)
          endif
        endif


      endif ! h_avail>0
    endif ; enddo ! i loop
  enddo ! k loop

end subroutine mixedlayer_convection

!>   This subroutine determines the TKE available at the depth of free
!! convection to drive mechanical entrainment.
subroutine find_starting_TKE(htot, h_CA, fluxes, Conv_En, cTKE, dKE_FC, dKE_CA, &
                             TKE, TKE_river, Idecay_len_TKE, cMKE, dt, Idt_diag, &
                             j, ksort, G, GV, US, CS)
  type(ocean_grid_type),      intent(in)    :: G       !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV      !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US      !< A dimensional unit scaling type
  real, dimension(SZI_(G)),   intent(in)    :: htot    !< The accumulated mixed layer thickness
                                                       !! [H ~> m or kg m-2]
  real, dimension(SZI_(G)),   intent(in)    :: h_CA    !< The mixed layer depth after convective
                                                       !! adjustment [H ~> m or kg m-2].
  type(forcing),              intent(in)    :: fluxes  !< A structure containing pointers to any
                                                       !! possible forcing fields.  Unused fields
                                                       !! have NULL ptrs.
  real, dimension(SZI_(G)),   intent(inout) :: Conv_En !< The buoyant turbulent kinetic energy source
                                                       !! due to free convection [Z L2 T-2 ~> m3 s-2].
  real, dimension(SZI_(G)),   intent(in)    :: dKE_FC  !< The vertically integrated change in
                                                       !! kinetic energy due to free convection
                                                       !! [Z L2 T-2 ~> m3 s-2].
  real, dimension(SZI_(G),SZK_(GV)), &
                              intent(in)    :: cTKE    !< The buoyant turbulent kinetic energy
                                                       !! source due to convective adjustment
                                                       !! [Z L2 T-2 ~> m3 s-2].
  real, dimension(SZI_(G),SZK_(GV)), &
                              intent(in)    :: dKE_CA  !< The vertically integrated change in
                                                       !! kinetic energy due to convective
                                                       !! adjustment [Z L2 T-2 ~> m3 s-2].
  real, dimension(SZI_(G)),   intent(out)   :: TKE     !< The turbulent kinetic energy available for
                                                       !! mixing over a time step [Z m2 T-2 ~> m3 s-2].
  real, dimension(SZI_(G)),   intent(out)   :: Idecay_len_TKE !< The inverse of the vertical decay
                                                       !! scale for TKE [H-1 ~> m-1 or m2 kg-1].
  real, dimension(SZI_(G)),   intent(in)    :: TKE_river !< The source of turbulent kinetic energy
                                                       !! available for driving mixing at river mouths
                                                       !! [Z L2 T-3 ~> m3 s-3].
  real, dimension(2,SZI_(G)), intent(out)   :: cMKE    !< Coefficients of HpE and HpE^2 in
                                                       !! calculating the denominator of MKE_rate,
                                                       !! [H-1 ~> m-1 or m2 kg-1] and [H-2 ~> m-2 or m4 kg-2].
  real,                       intent(in)    :: dt      !< The time step [T ~> s].
  real,                       intent(in)    :: Idt_diag !< The inverse of the accumulated diagnostic
                                                       !! time interval [T-1 ~> s-1].
  integer,                    intent(in)    :: j       !< The j-index to work on.
  integer, dimension(SZI_(G),SZK_(GV)), &
                              intent(in)    :: ksort   !< The density-sorted k-indicies.
  type(bulkmixedlayer_CS),    pointer       :: CS      !< The control structure for this module.

!   This subroutine determines the TKE available at the depth of free
! convection to drive mechanical entrainment.

  ! Local variables
  real :: dKE_conv  ! The change in mean kinetic energy due to all convection [Z L2 T-2 ~> m3 s-2].
  real :: nstar_FC  ! The effective efficiency with which the energy released by
                    ! free convection is converted to TKE, often ~0.2 [nondim].
  real :: nstar_CA  ! The effective efficiency with which the energy released by
                    ! convective adjustment is converted to TKE, often ~0.2 [nondim].
  real :: TKE_CA    ! The potential energy released by convective adjustment if
                    ! that release is positive [Z L2 T-2 ~> m3 s-2].
  real :: MKE_rate_CA ! MKE_rate for convective adjustment [nondim], 0 to 1.
  real :: MKE_rate_FC ! MKE_rate for free convection [nondim], 0 to 1.
  real :: totEn_Z   ! The total potential energy released by convection, [Z3 T-2 ~> m3 s-2].
  real :: Ih        ! The inverse of a thickness [H-1 ~> m-1 or m2 kg-1].
  real :: exp_kh    ! The nondimensional decay of TKE across a layer [nondim].
  real :: absf      ! The absolute value of f averaged to thickness points [T-1 ~> s-1].
  real :: U_star    ! The friction velocity [Z T-1 ~> m s-1].
  real :: absf_Ustar  ! The absolute value of f divided by U_star [Z-1 ~> m-1].
  real :: wind_TKE_src ! The surface wind source of TKE [Z L2 T-3 ~> m3 s-3].
  real :: diag_wt   ! The ratio of the current timestep to the diagnostic
                    ! timestep (which may include 2 calls) [nondim].
  integer :: is, ie, nz, i

  is = G%isc ; ie = G%iec ; nz = GV%ke
  diag_wt = dt * Idt_diag

  if (CS%omega_frac >= 1.0) absf = 2.0*CS%omega
  do i=is,ie
    U_star = fluxes%ustar(i,j)
    if (associated(fluxes%ustar_shelf) .and. associated(fluxes%frac_shelf_h)) then
      if (fluxes%frac_shelf_h(i,j) > 0.0) &
        U_star = (1.0 - fluxes%frac_shelf_h(i,j)) * U_star + &
                  fluxes%frac_shelf_h(i,j) * fluxes%ustar_shelf(i,j)
    endif

    if (U_star < CS%ustar_min) U_star = CS%ustar_min
    if (CS%omega_frac < 1.0) then
      absf = 0.25*((abs(G%CoriolisBu(I,J)) + abs(G%CoriolisBu(I-1,J-1))) + &
                   (abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I-1,J))))
      if (CS%omega_frac > 0.0) &
        absf = sqrt(CS%omega_frac*4.0*CS%omega**2 + (1.0-CS%omega_frac)*absf**2)
    endif
    absf_Ustar = absf / U_star
    Idecay_len_TKE(i) = (absf_Ustar * CS%TKE_decay) * GV%H_to_Z

!    The first number in the denominator could be anywhere up to 16.  The
!  value of 3 was chosen to minimize the time-step dependence of the amount
!  of shear-driven mixing in 10 days of a 1-degree global model, emphasizing
!  the equatorial areas.  Although it is not cast as a parameter, it should
!  be considered an empirical parameter, and it might depend strongly on the
!  number of sublayers in the mixed layer and their locations.
!  The 0.41 is VonKarman's constant.  This equation assumes that small & large
!  scales contribute to mixed layer deepening at similar rates, even though
!  small scales are dissipated more rapidly (implying they are less efficient).
!     Ih = 1.0/(16.0*0.41*U_star*dt)
    Ih = GV%H_to_Z/(3.0*0.41*U_star*dt)
    cMKE(1,i) = 4.0 * Ih ; cMKE(2,i) = (absf_Ustar*GV%H_to_Z) * Ih

    if (Idecay_len_TKE(i) > 0.0) then
      exp_kh = exp(-htot(i)*Idecay_len_TKE(i))
    else
      exp_kh = 1.0
    endif

! Here nstar is a function of the natural Rossby number  0.2/(1+0.2/Ro), based
! on a curve fit from the data of Wang (GRL, 2003).
! Note:         Ro = 1.0/sqrt(0.5 * dt * (absf*htot(i))**3 / totEn)
    if (Conv_En(i) < 0.0) Conv_En(i) = 0.0
    if (cTKE(i,1) > 0.0) then ; TKE_CA = cTKE(i,1) ; else ; TKE_CA = 0.0 ; endif
    if ((htot(i) >= h_CA(i)) .or. (TKE_CA == 0.0)) then
      totEn_Z = US%L_to_Z**2 * (Conv_En(i) + TKE_CA)

      if (totEn_Z > 0.0) then
        nstar_FC = CS%nstar * totEn_Z / (totEn_Z + 0.2 * &
                        sqrt(0.5 * dt * (absf*(htot(i)*GV%H_to_Z))**3 * totEn_Z))
      else
        nstar_FC = CS%nstar
      endif
      nstar_CA = nstar_FC
    else
      ! This reconstructs the Buoyancy flux within the topmost htot of water.
      if (Conv_En(i) > 0.0) then
        totEn_Z = US%L_to_Z**2 * (Conv_En(i) + TKE_CA * (htot(i) / h_CA(i)) )
        nstar_FC = CS%nstar * totEn_Z / (totEn_Z + 0.2 * &
                        sqrt(0.5 * dt * (absf*(htot(i)*GV%H_to_Z))**3 * totEn_Z))
      else
        nstar_FC = CS%nstar
      endif

      totEn_Z = US%L_to_Z**2 * (Conv_En(i) + TKE_CA)
      if (TKE_CA > 0.0) then
        nstar_CA = CS%nstar * totEn_Z / (totEn_Z + 0.2 * &
                        sqrt(0.5 * dt * (absf*(h_CA(i)*GV%H_to_Z))**3 * totEn_Z))
      else
        nstar_CA = CS%nstar
      endif
    endif

    if (dKE_FC(i) + dKE_CA(i,1) > 0.0) then
      if (htot(i) >= h_CA(i)) then
        MKE_rate_FC = 1.0 / (1.0 + htot(i)*(cMKE(1,i) + cMKE(2,i)*htot(i)) )
        MKE_rate_CA = MKE_rate_FC
      else
        MKE_rate_FC = 1.0 / (1.0 + htot(i)*(cMKE(1,i) + cMKE(2,i)*htot(i)) )
        MKE_rate_CA = 1.0 / (1.0 + h_CA(i)*(cMKE(1,i) + cMKE(2,i)*h_CA(i)) )
      endif
    else
      ! This branch just saves unnecessary calculations.
      MKE_rate_FC = 1.0 ; MKE_rate_CA = 1.0
    endif

    dKE_conv = dKE_CA(i,1) * MKE_rate_CA + dKE_FC(i) * MKE_rate_FC
! At this point, it is assumed that cTKE is positive and stored in TKE_CA!
! Note: Removed factor of 2 in u*^3 terms.
    TKE(i) = (dt*CS%mstar)*((US%Z_to_L**2*(U_star*U_Star*U_Star))*exp_kh) + &
             (exp_kh * dKE_conv + nstar_FC*Conv_En(i) + nstar_CA * TKE_CA)

    if (CS%do_rivermix) then ! Add additional TKE at river mouths
      TKE(i) = TKE(i) + TKE_river(i)*dt*exp_kh
    endif

    if (CS%TKE_diagnostics) then
      wind_TKE_src = CS%mstar*(US%Z_to_L**2*U_star*U_Star*U_Star) * diag_wt
      CS%diag_TKE_wind(i,j) = CS%diag_TKE_wind(i,j) + &
          ( wind_TKE_src + TKE_river(i) * diag_wt )
      CS%diag_TKE_RiBulk(i,j) = CS%diag_TKE_RiBulk(i,j) + dKE_conv*Idt_diag
      CS%diag_TKE_mech_decay(i,j) = CS%diag_TKE_mech_decay(i,j) + &
          (exp_kh-1.0)*(wind_TKE_src + dKE_conv*Idt_diag)
      CS%diag_TKE_conv(i,j) = CS%diag_TKE_conv(i,j) + &
          Idt_diag * (nstar_FC*Conv_En(i) + nstar_CA*TKE_CA)
      CS%diag_TKE_conv_decay(i,j) = CS%diag_TKE_conv_decay(i,j) + &
          Idt_diag * ((CS%nstar-nstar_FC)*Conv_En(i) + (CS%nstar-nstar_CA)*TKE_CA)
      CS%diag_TKE_conv_s2(i,j) = CS%diag_TKE_conv_s2(i,j) + &
          Idt_diag * (cTKE(i,1)-TKE_CA)
    endif
  enddo

end subroutine find_starting_TKE

!> This subroutine calculates mechanically driven entrainment.
subroutine mechanical_entrainment(h, d_eb, htot, Ttot, Stot, uhtot, vhtot, &
                                  R0_tot, Rcv_tot, u, v, T, S, R0, Rcv, eps, &
                                  dR0_dT, dRcv_dT, cMKE, Idt_diag, nsw, &
                                  Pen_SW_bnd, opacity_band, TKE, &
                                  Idecay_len_TKE, j, ksort, G, GV, US, CS)
  type(ocean_grid_type),    intent(in)    :: G     !< The ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV    !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)    :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(inout) :: h     !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(inout) :: d_eb  !< The downward increase across a layer in the
                                                   !! layer in the entrainment from below [H ~> m or kg m-2].
                                                   !! Positive values go with mass gain by a layer.
  real, dimension(SZI_(G)), intent(inout) :: htot  !< The accumlated mixed layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G)), intent(inout) :: Ttot  !< The depth integrated mixed layer temperature
                                                   !! [degC H ~> degC m or degC kg m-2].
  real, dimension(SZI_(G)), intent(inout) :: Stot  !< The depth integrated mixed layer salinity
                                                   !! [ppt H ~> ppt m or ppt kg m-2].
  real, dimension(SZI_(G)), intent(inout) :: uhtot !< The depth integrated mixed layer zonal
                                                   !! velocity [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G)), intent(inout) :: vhtot !< The integrated mixed layer meridional
                                                   !! velocity [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G)), intent(inout) :: R0_tot !< The integrated mixed layer potential density
                                                   !! referenced to 0 pressure [H R ~> kg m-2 or kg2 m-5].
  real, dimension(SZI_(G)), intent(inout) :: Rcv_tot !< The integrated mixed layer coordinate variable
                                                   !! potential density [H R ~> kg m-2 or kg2 m-5].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: u     !< Zonal velocities interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: v     !< Zonal velocities interpolated to h points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: T     !< Layer temperatures [degC].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: S     !< Layer salinities [ppt].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: R0    !< Potential density referenced to
                                                   !! surface pressure [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: Rcv   !< The coordinate defining potential
                                                   !! density [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: eps   !< The negligibly small amount of water
                                                   !! that will be left in each layer [H ~> m or kg m-2].
  real, dimension(SZI_(G)), intent(in)    :: dR0_dT  !< The partial derivative of R0 with respect to
                                                   !! temperature [R degC-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)), intent(in)    :: dRcv_dT !< The partial derivative of Rcv with respect to
                                                   !! temperature [R degC-1 ~> kg m-3 degC-1].
  real, dimension(2,SZI_(G)), intent(in)  :: cMKE  !< Coefficients of HpE and HpE^2 used in calculating the
                                                   !! denominator of MKE_rate; the two elements have differing
                                                   !! units of [H-1 ~> m-1 or m2 kg-1] and [H-2 ~> m-2 or m4 kg-2].
  real,                     intent(in)    :: Idt_diag !< The inverse of the accumulated diagnostic
                                                   !! time interval [T-1 ~> s-1].
  integer,                  intent(in)    :: nsw   !< The number of bands of penetrating
                                                   !! shortwave radiation.
  real, dimension(max(nsw,1),SZI_(G)), intent(inout) :: Pen_SW_bnd !< The penetrating shortwave
                                                   !! heating at the sea surface in each penetrating
                                                   !! band [degC H ~> degC m or degC kg m-2].
  real, dimension(max(nsw,1),SZI_(G),SZK_(GV)), intent(in) :: opacity_band !< The opacity in each band of
                                                   !! penetrating shortwave radiation [H-1 ~> m-1 or m2 kg-1].
  real, dimension(SZI_(G)), intent(inout) :: TKE   !< The turbulent kinetic energy
                                                   !! available for mixing over a time
                                                   !! step [Z m2 T-2 ~> m3 s-2].
  real, dimension(SZI_(G)), intent(inout) :: Idecay_len_TKE !< The vertical TKE decay rate [H-1 ~> m-1 or m2 kg-1].
  integer,                  intent(in)    :: j     !< The j-index to work on.
  integer, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)    :: ksort !< The density-sorted k-indicies.
  type(bulkmixedlayer_CS),  pointer       :: CS    !< The control structure for this module.

! This subroutine calculates mechanically driven entrainment.

  ! Local variables
  real :: SW_trans  !   The fraction of shortwave radiation that is not
                    ! absorbed in a layer, nondimensional.
  real :: Pen_absorbed  !   The amount of penetrative shortwave radiation
                        ! that is absorbed in a layer [degC H ~> degC m or degC kg m-2].
  real :: h_avail   ! The thickness in a layer available for entrainment [H ~> m or kg m-2].
  real :: h_ent     ! The thickness from a layer that is entrained [H ~> m or kg m-2].
  real :: h_min, h_max ! Limits on the solution for h_ent [H ~> m or kg m-2].
  real :: dh_Newt      !   The Newton's method estimate of the change in
                       ! h_ent between iterations [H ~> m or kg m-2].
  real :: MKE_rate  !   The fraction of the energy in resolved shears
                    ! within the mixed layer that will be eliminated
                    ! within a timestep, nondim, 0 to 1.
  real :: HpE       !   The current thickness plus entrainment [H ~> m or kg m-2].
  real :: g_H_2Rho0   !   Half the gravitational acceleration times the
                      ! conversion from H to m divided by the mean density,
                      ! in [L2 T-2 H-1 R-1 ~> m4 s-2 kg-1 or m7 s-2 kg-2].
  real :: TKE_full_ent  ! The TKE remaining if a layer is fully entrained
                        ! [Z L2 T-2 ~> m3 s-2].
  real :: dRL       ! Work required to mix water from the next layer
                    ! across the mixed layer [L2 T-2 ~> L2 s-2].
  real :: Pen_En_Contrib  ! Penetrating SW contributions to the changes in
                          ! TKE, divided by layer thickness in m [L2 T2 ~> m2 s-2].
  real :: Cpen1     ! A temporary variable [L2 T-2 ~> m2 s-2].
  real :: dMKE      ! A temporary variable related to the release of mean
                    ! kinetic energy [H Z L2 T-2 ~> m4 s-2 or kg m s-2]
  real :: TKE_ent   ! The TKE that remains if h_ent were entrained [Z L2 T-2 ~> m3 s-2].
  real :: TKE_ent1  ! The TKE that would remain, without considering the
                    ! release of mean kinetic energy [Z L2 T-2 ~> m3 s-2].
  real :: dTKE_dh   ! The partial derivative of TKE with h_ent [Z L2 T-2 H-1 ~> m2 s-2 or m5 s-2 kg-1].
  real :: Pen_dTKE_dh_Contrib ! The penetrating shortwave contribution to
                    ! dTKE_dh [L2 T-2 ~> m2 s-2].
  real :: EF4_val   ! The result of EF4() (see later) [H-1 ~> m-1 or m2 kg-1].
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: dEF4_dh   ! The partial derivative of EF4 with h [H-2 ~> m-2 or m4 kg-2].
  real :: Pen_En1   ! A nondimensional temporary variable.
  real :: kh, exp_kh  ! Nondimensional temporary variables related to the
  real :: f1_kh       ! fractional decay of TKE across a layer.
  real :: x1, e_x1      !   Nondimensional temporary variables related to
  real :: f1_x1, f2_x1  ! the relative decay of TKE and SW radiation across
  real :: f3_x1         ! a layer, and exponential-related functions of x1.
  real :: E_HxHpE   ! Entrainment divided by the product of the new and old
                    ! thicknesses [H-1 ~> m-1 or m2 kg-1].
  real :: Hmix_min  ! The minimum mixed layer depth [H ~> m or kg m-2].
  real :: opacity
  real :: C1_3, C1_6, C1_24   !  1/3, 1/6, and 1/24.
  integer :: is, ie, nz, i, k, ks, itt, n

  C1_3 = 1.0/3.0 ; C1_6 = 1.0/6.0 ; C1_24 = 1.0/24.0
  g_H_2Rho0 = (GV%g_Earth * GV%H_to_Z) / (2.0 * GV%Rho0)
  Hmix_min = CS%Hmix_min
  h_neglect = GV%H_subroundoff
  is = G%isc ; ie = G%iec ; nz = GV%ke

  do ks=1,nz

    do i=is,ie ; if (ksort(i,ks) > 0) then
      k = ksort(i,ks)

      h_avail = h(i,k) - eps(i,k)
      if ((h_avail > 0.) .and. ((TKE(i) > 0.) .or. (htot(i) < Hmix_min))) then
        dRL = g_H_2Rho0 * (R0(i,k)*htot(i) - R0_tot(i) )
        dMKE = (GV%H_to_Z * CS%bulk_Ri_ML) * 0.5 * &
            ((uhtot(i)-u(i,k)*htot(i))**2 + (vhtot(i)-v(i,k)*htot(i))**2)

! Find the TKE that would remain if the entire layer were entrained.
        kh = Idecay_len_TKE(i)*h_avail ; exp_kh = exp(-kh)
        if (kh >= 2.0e-5) then ; f1_kh = (1.0-exp_kh) / kh
        else ; f1_kh = (1.0 - kh*(0.5 - C1_6*kh)) ; endif

        Pen_En_Contrib = 0.0
        do n=1,nsw ; if (Pen_SW_bnd(n,i) > 0.0) then
          opacity = opacity_band(n,i,k)
! Two different forms are used here to make sure that only negative
! values are taken into exponentials to avoid excessively large
! numbers.  They are, of course, mathematically identical.
          if (Idecay_len_TKE(i) > opacity) then
            x1 = (Idecay_len_TKE(i) - opacity) * h_avail
            if (x1 >= 2.0e-5) then
              e_x1 = exp(-x1) ; f1_x1 = ((1.0-e_x1)/(x1))
              f3_x1 = ((e_x1-(1.0-x1))/(x1*x1))
            else
              f1_x1 = (1.0 - x1*(0.5 - C1_6*x1))
              f3_x1 = (0.5 - x1*(C1_6 - C1_24*x1))
            endif

            Pen_En1 = exp(-opacity*h_avail) * &
               ((1.0+opacity*htot(i))*f1_x1 + opacity*h_avail*f3_x1)
          else
            x1 = (opacity - Idecay_len_TKE(i)) * h_avail
            if (x1 >= 2.0e-5) then
              e_x1 = exp(-x1) ; f1_x1 = ((1.0-e_x1)/(x1))
              f2_x1 = ((1.0-(1.0+x1)*e_x1)/(x1*x1))
            else
              f1_x1 = (1.0 - x1*(0.5 - C1_6*x1))
              f2_x1 = (0.5 - x1*(C1_3 - 0.125*x1))
            endif

            Pen_En1 = exp_kh * ((1.0+opacity*htot(i))*f1_x1 + &
                                 opacity*h_avail*f2_x1)
          endif
          Pen_En_Contrib = Pen_En_Contrib + &
            (g_H_2Rho0*dR0_dT(i)*Pen_SW_bnd(n,i)) * (Pen_En1 - f1_kh)
        endif ; enddo

        HpE = htot(i)+h_avail
        MKE_rate = 1.0/(1.0 + (cMKE(1,i)*HpE + cMKE(2,i)*HpE**2))
        EF4_val = EF4(htot(i)+h_neglect,h_avail,Idecay_len_TKE(i))
        TKE_full_ent = (exp_kh*TKE(i) - (h_avail*GV%H_to_Z)*(dRL*f1_kh + Pen_En_Contrib)) + &
            MKE_rate*dMKE*EF4_val
        if ((TKE_full_ent >= 0.0) .or. (h_avail+htot(i) <= Hmix_min)) then
          ! The layer will be fully entrained.
          h_ent = h_avail

          if (CS%TKE_diagnostics) then
            E_HxHpE = h_ent / ((htot(i)+h_neglect)*(htot(i)+h_ent+h_neglect))
            CS%diag_TKE_mech_decay(i,j) = CS%diag_TKE_mech_decay(i,j) + &
                Idt_diag * ((exp_kh-1.0)* TKE(i) + (h_ent*GV%H_to_Z)*dRL*(1.0-f1_kh) + &
                            MKE_rate*dMKE*(EF4_val-E_HxHpE))
            CS%diag_TKE_mixing(i,j) = CS%diag_TKE_mixing(i,j) - &
                Idt_diag*(GV%H_to_Z*h_ent)*dRL
            CS%diag_TKE_pen_SW(i,j) = CS%diag_TKE_pen_SW(i,j) - &
                Idt_diag*(GV%H_to_Z*h_ent)*Pen_En_Contrib
            CS%diag_TKE_RiBulk(i,j) = CS%diag_TKE_RiBulk(i,j) + &
                Idt_diag*MKE_rate*dMKE*E_HxHpE
          endif

          TKE(i) = TKE_full_ent
          !### The minimum TKE value in this line may be problematically small.
          if (TKE(i) <= 0.0) TKE(i) = 1.0e-150*US%m_to_Z*US%m_s_to_L_T**2
        else
! The layer is only partially entrained.  The amount that will be
! entrained is determined iteratively.  No further layers will be
! entrained.
          h_min = 0.0 ; h_max = h_avail
          if (TKE(i) <= 0.0) then
            h_ent = 0.0
          else
            h_ent = h_avail * TKE(i) / (TKE(i) - TKE_full_ent)

            do itt=1,15
              ! Evaluate the TKE that would remain if h_ent were entrained.

              kh = Idecay_len_TKE(i)*h_ent ; exp_kh = exp(-kh)
              if (kh >= 2.0e-5) then
                f1_kh = (1.0-exp_kh) / kh
              else
                f1_kh = (1.0 - kh*(0.5 - C1_6*kh))
              endif


              Pen_En_Contrib = 0.0 ; Pen_dTKE_dh_Contrib = 0.0
              do n=1,nsw ; if (Pen_SW_bnd(n,i) > 0.0) then
                ! Two different forms are used here to make sure that only negative
                ! values are taken into exponentials to avoid excessively large
                ! numbers.  They are, of course, mathematically identical.
                opacity = opacity_band(n,i,k)
                SW_trans = exp(-h_ent*opacity)
                if (Idecay_len_TKE(i) > opacity) then
                  x1 = (Idecay_len_TKE(i) - opacity) * h_ent
                  if (x1 >= 2.0e-5) then
                    e_x1 = exp(-x1) ; f1_x1 = ((1.0-e_x1)/(x1))
                    f3_x1 = ((e_x1-(1.0-x1))/(x1*x1))
                  else
                    f1_x1 = (1.0 - x1*(0.5 - C1_6*x1))
                    f3_x1 = (0.5 - x1*(C1_6 - C1_24*x1))
                  endif
                  Pen_En1 = SW_trans * ((1.0+opacity*htot(i))*f1_x1 + &
                                          opacity*h_ent*f3_x1)
                else
                  x1 = (opacity - Idecay_len_TKE(i)) * h_ent
                  if (x1 >= 2.0e-5) then
                    e_x1 = exp(-x1) ; f1_x1 = ((1.0-e_x1)/(x1))
                    f2_x1 = ((1.0-(1.0+x1)*e_x1)/(x1*x1))
                  else
                    f1_x1 = (1.0 - x1*(0.5 - C1_6*x1))
                    f2_x1 = (0.5 - x1*(C1_3 - 0.125*x1))
                  endif

                  Pen_En1 = exp_kh * ((1.0+opacity*htot(i))*f1_x1 + &
                                        opacity*h_ent*f2_x1)
                endif
                Cpen1 = g_H_2Rho0*dR0_dT(i)*Pen_SW_bnd(n,i)
                Pen_En_Contrib = Pen_En_Contrib + Cpen1*(Pen_En1 - f1_kh)
                Pen_dTKE_dh_Contrib = Pen_dTKE_dh_Contrib + &
                           Cpen1*((1.0-SW_trans) - opacity*(htot(i) + h_ent)*SW_trans)
              endif ; enddo ! (Pen_SW_bnd(n,i) > 0.0)

              TKE_ent1 = exp_kh* TKE(i) - (h_ent*GV%H_to_Z)*(dRL*f1_kh + Pen_En_Contrib)
              EF4_val = EF4(htot(i)+h_neglect,h_ent,Idecay_len_TKE(i),dEF4_dh)
              HpE = htot(i)+h_ent
              MKE_rate = 1.0/(1.0 + (cMKE(1,i)*HpE + cMKE(2,i)*HpE**2))
              TKE_ent = TKE_ent1 + dMKE*EF4_val*MKE_rate
              ! TKE_ent is the TKE that would remain if h_ent were entrained.

              dTKE_dh = ((-Idecay_len_TKE(i)*TKE_ent1 - dRL*GV%H_to_Z) + &
                         Pen_dTKE_dh_Contrib*GV%H_to_Z) + dMKE * MKE_rate* &
                        (dEF4_dh - EF4_val*MKE_rate*(cMKE(1,i)+2.0*cMKE(2,i)*HpE))
              !  dh_Newt = -TKE_ent / dTKE_dh
              ! Bisect if the Newton's method prediction is outside of the bounded range.
              if (TKE_ent > 0.0) then
                if ((h_max-h_ent)*(-dTKE_dh) > TKE_ent) then
                  dh_Newt = -TKE_ent / dTKE_dh
                else
                  dh_Newt = 0.5*(h_max-h_ent)
                endif
                h_min = h_ent
              else
                if ((h_min-h_ent)*(-dTKE_dh) < TKE_ent) then
                  dh_Newt = -TKE_ent / dTKE_dh
                else
                  dh_Newt = 0.5*(h_min-h_ent)
                endif
                h_max = h_ent
              endif
              h_ent = h_ent + dh_Newt

              if (ABS(dh_Newt) < 0.2*GV%Angstrom_H) exit
            enddo
          endif

          if (h_ent < Hmix_min-htot(i)) h_ent = Hmix_min - htot(i)

          if (CS%TKE_diagnostics) then
            HpE = htot(i)+h_ent
            MKE_rate = 1.0/(1.0 + cMKE(1,i)*HpE + cMKE(2,i)*HpE**2)
            EF4_val = EF4(htot(i)+h_neglect,h_ent,Idecay_len_TKE(i))

            E_HxHpE = h_ent / ((htot(i)+h_neglect)*(HpE+h_neglect))
            CS%diag_TKE_mech_decay(i,j) = CS%diag_TKE_mech_decay(i,j) + &
                Idt_diag * ((exp_kh-1.0)* TKE(i) + (h_ent*GV%H_to_Z)*dRL*(1.0-f1_kh) + &
                             dMKE*MKE_rate*(EF4_val-E_HxHpE))
            CS%diag_TKE_mixing(i,j) = CS%diag_TKE_mixing(i,j) - &
                Idt_diag*(h_ent*GV%H_to_Z)*dRL
            CS%diag_TKE_pen_SW(i,j) = CS%diag_TKE_pen_SW(i,j) - &
                Idt_diag*(h_ent*GV%H_to_Z)*Pen_En_Contrib
            CS%diag_TKE_RiBulk(i,j) = CS%diag_TKE_RiBulk(i,j) + &
                Idt_diag*dMKE*MKE_rate*E_HxHpE
          endif

          TKE(i) = 0.0
        endif ! TKE_full_ent > 0.0

        Pen_absorbed = 0.0
        do n=1,nsw ; if (Pen_SW_bnd(n,i) > 0.0) then
          SW_trans = exp(-h_ent*opacity_band(n,i,k))
          Pen_absorbed = Pen_absorbed + Pen_SW_bnd(n,i) * (1.0 - SW_trans)
          Pen_SW_bnd(n,i) = Pen_SW_bnd(n,i) * SW_trans
        endif ; enddo

        htot(i)   = htot(i)   + h_ent
        R0_tot(i) = R0_tot(i) + (h_ent * R0(i,k) + Pen_absorbed*dR0_dT(i))
        h(i,k)    = h(i,k)    - h_ent
        d_eb(i,k) = d_eb(i,k) - h_ent

        Stot(i)    = Stot(i)    + h_ent * S(i,k)
        Ttot(i)    = Ttot(i)    + (h_ent * T(i,k) + Pen_absorbed)
        Rcv_tot(i) = Rcv_tot(i) + (h_ent*Rcv(i,k) + Pen_absorbed*dRcv_dT(i))

        uhtot(i) = uhtot(i) + u(i,k)*h_ent
        vhtot(i) = vhtot(i) + v(i,k)*h_ent
      endif ! h_avail > 0.0 .AND TKE(i) > 0.0

    endif ; enddo ! i loop
  enddo ! k loop

end subroutine mechanical_entrainment

!> This subroutine generates an array of indices that are sorted by layer
!! density.
subroutine sort_ML(h, R0, eps, G, GV, CS, ksort)
  type(ocean_grid_type),                intent(in)  :: G     !< The ocean's grid structure.
  type(verticalGrid_type),              intent(in)  :: GV    !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZK_(GV)),    intent(in)  :: h     !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)),    intent(in)  :: R0    !< The potential density used to sort
                                                             !! the layers [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)),    intent(in)  :: eps   !< The (small) thickness that must
                                                             !! remain in each layer [H ~> m or kg m-2].
  type(bulkmixedlayer_CS),              pointer     :: CS    !< The control structure returned by a
                                                             !! previous call to mixedlayer_init.
  integer, dimension(SZI_(G),SZK_(GV)), intent(out) :: ksort !< The k-index to use in the sort.

  ! Local variables
  real :: R0sort(SZI_(G),SZK_(GV))
  integer :: nsort(SZI_(G))
  logical :: done_sorting(SZI_(G))
  integer :: i, k, ks, is, ie, nz, nkmb

  is = G%isc ; ie = G%iec ; nz = GV%ke
  nkmb = CS%nkml+CS%nkbl

  !   Come up with a sorted index list of layers with increasing R0.
  ! Assume that the layers below nkmb are already stably stratified.
  ! Only layers that are thicker than eps are in the list.  Extra elements
  ! have an index of -1.

  !   This is coded using straight insertion, on the assumption that the
  ! layers are usually in the right order (or massless) anyway.

  do k=1,nz ; do i=is,ie ; ksort(i,k) = -1 ; enddo ; enddo

  do i=is,ie ; nsort(i) = 0 ; done_sorting(i) = .false. ; enddo
  do k=1,nz ; do i=is,ie ; if (h(i,k) > eps(i,k)) then
    if (done_sorting(i)) then ; ks = nsort(i) ; else
      do ks=nsort(i),1,-1
        if (R0(i,k) >= R0sort(i,ks)) exit
        R0sort(i,ks+1) = R0sort(i,ks) ; ksort(i,ks+1) = ksort(i,ks)
      enddo
      if ((k > nkmb) .and. (ks == nsort(i))) done_sorting(i) = .true.
    endif

    ksort(i,ks+1) = k
    R0sort(i,ks+1) = R0(i,k)
    nsort(i) = nsort(i) + 1
  endif ; enddo ; enddo

end subroutine sort_ML

!>   This subroutine actually moves properties between layers to achieve a
!! resorted state, with all of the resorted water either moved into the correct
!! interior layers or in the top nkmb layers.
subroutine resort_ML(h, T, S, R0, Rcv, RcvTgt, eps, d_ea, d_eb, ksort, G, GV, CS, &
                     dR0_dT, dR0_dS, dRcv_dT, dRcv_dS)
  type(ocean_grid_type),                intent(in)    :: G       !< The ocean's grid structure.
  type(verticalGrid_type),              intent(in)    :: GV      !< The ocean's vertical grid
                                                                 !! structure.
  real, dimension(SZI_(G),SZK0_(GV)),   intent(inout) :: h       !< Layer thickness [H ~> m or kg m-2].
                                                                 !! Layer 0 is the new mixed layer.
  real, dimension(SZI_(G),SZK0_(GV)),   intent(inout) :: T       !< Layer temperatures [degC].
  real, dimension(SZI_(G),SZK0_(GV)),   intent(inout) :: S       !< Layer salinities [ppt].
  real, dimension(SZI_(G),SZK0_(GV)),   intent(inout) :: R0      !< Potential density referenced to
                                                                 !! surface pressure [R ~> kg m-3].
  real, dimension(SZI_(G),SZK0_(GV)),   intent(inout) :: Rcv     !< The coordinate defining
                                                                 !! potential density [R ~> kg m-3].
  real, dimension(SZK_(GV)),            intent(in)    :: RcvTgt  !< The target value of Rcv for each
                                                                 !! layer [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)),    intent(inout) :: eps     !< The (small) thickness that must
                                                                 !! remain in each layer [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZK_(GV)),    intent(inout) :: d_ea    !< The upward increase across a
                                                                 !! layer in the entrainment from
                                                                 !! above [H ~> m or kg m-2].
                                                                 !!  Positive d_ea goes with layer
                                                                 !! thickness increases.
  real, dimension(SZI_(G),SZK_(GV)),    intent(inout) :: d_eb    !< The downward increase across a
                                                                 !! layer in the entrainment from
                                                                 !! below [H ~> m or kg m-2]. Positive values go
                                                                 !! with mass gain by a layer.
  integer, dimension(SZI_(G),SZK_(GV)), intent(in)    :: ksort   !< The density-sorted k-indicies.
  type(bulkmixedlayer_CS),              pointer       :: CS      !< The control structure for this
                                                                 !! module.
  real, dimension(SZI_(G)),             intent(in)    :: dR0_dT  !< The partial derivative of
                                                                 !! potential density referenced
                                                                 !! to the surface with potential
                                                                 !! temperature [R degC-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)),             intent(in)    :: dR0_dS  !< The partial derivative of
                                                                 !! cpotential density referenced
                                                                 !! to the surface with salinity,
                                                                 !! [R ppt-1 ~> kg m-3 ppt-1].
  real, dimension(SZI_(G)),             intent(in)    :: dRcv_dT !< The partial derivative of
                                                                 !! coordinate defining potential
                                                                 !! density with potential
                                                                 !! temperature [R degC-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)),             intent(in)    :: dRcv_dS !< The partial derivative of
                                                                 !! coordinate defining potential
                                                                 !! density with salinity,
                                                                 !! [R ppt-1 ~> kg m-3 ppt-1].

!   If there are no massive light layers above the deepest of the mixed- and
! buffer layers, do nothing (except perhaps to reshuffle these layers).
!   If there are nkbl or fewer layers above the deepest mixed- or buffer-
! layers, move them (in sorted order) into the buffer layers, even if they
! were previously interior layers.
!   If there are interior layers that are intermediate in density (both in-situ
! and the coordinate density (sigma-2)) between the newly forming mixed layer
! and a residual buffer- or mixed layer, and the number of massive layers above
! the deepest massive buffer or mixed layer is greater than nkbl, then split
! those buffer layers into peices that match the target density of the two
! nearest interior layers.
!   Otherwise, if there are more than nkbl+1 remaining massive layers

  ! Local variables
  real    :: h_move, h_tgt_old, I_hnew
  real    :: dT_dS_wt2, dT_dR, dS_dR, I_denom
  real    :: Rcv_int
  real    :: T_up, S_up, R0_up, I_hup, h_to_up
  real    :: T_dn, S_dn, R0_dn, I_hdn, h_to_dn
  real    :: wt_dn
  real    :: dR1, dR2
  real    :: dPE, hmin, min_dPE, min_hmin
  real, dimension(SZK_(GV)) :: &
    h_tmp, R0_tmp, T_tmp, S_tmp, Rcv_tmp
  integer :: ks_min
  logical :: sorted, leave_in_layer
  integer :: ks_deep(SZI_(G)), k_count(SZI_(G)), ks2_reverse(SZI_(G), SZK_(GV))
  integer :: ks2(SZK_(GV))
  integer :: i, k, ks, is, ie, nz, k1, k2, k_tgt, k_src, k_int_top
  integer :: nks, nkmb, num_interior, top_interior_ks

  is = G%isc ; ie = G%iec ; nz = GV%ke
  nkmb = CS%nkml+CS%nkbl

  dT_dS_wt2 = CS%dT_dS_wt**2

! Find out how many massive layers are above the deepest buffer or mixed layer.
  do i=is,ie ; ks_deep(i) = -1 ; k_count(i) = 0 ; enddo
  do ks=nz,1,-1 ; do i=is,ie ; if (ksort(i,ks) > 0) then
    k = ksort(i,ks)

    if (h(i,k) > eps(i,k)) then
      if (ks_deep(i) == -1) then
        if (k <= nkmb) then
          ks_deep(i) = ks ; k_count(i) = k_count(i) + 1
          ks2_reverse(i,k_count(i)) = k
        endif
      else
        k_count(i) = k_count(i) + 1
        ks2_reverse(i,k_count(i)) = k
      endif
    endif
  endif ; enddo ; enddo

  do i=is,ie ; if (k_count(i) > 1) then
    ! This column might need to be reshuffled.
    nks = k_count(i)

    ! Put ks2 in the right order and check whether reshuffling is needed.
    sorted = .true.
    ks2(nks) = ks2_reverse(i,1)
    do ks=nks-1,1,-1
      ks2(ks) = ks2_reverse(i,1+nks-ks)
      if (ks2(ks) > ks2(ks+1)) sorted = .false.
    enddo

    ! Go to the next column of no reshuffling is needed.
    if (sorted) cycle

    ! Find out how many interior layers are being reshuffled.  If none,
    ! then this is a simple swapping procedure.
    num_interior = 0 ; top_interior_ks = 0
    do ks=nks,1,-1 ; if (ks2(ks) > nkmb) then
      num_interior = num_interior+1 ; top_interior_ks = ks
    endif ; enddo

    if (num_interior >= 1) then
      ! Find the lightest interior layer with a target coordinate density
      ! greater than the newly forming mixed layer.
      do k=nkmb+1,nz ; if (Rcv(i,0) < RcvTgt(k)) exit ; enddo
      k_int_top = k ; Rcv_int = RcvTgt(k)

      I_denom = 1.0 / (dRcv_dS(i)**2 + dT_dS_wt2*dRcv_dT(i)**2)
      dT_dR = dT_dS_wt2*dRcv_dT(i) * I_denom
      dS_dR = dRcv_dS(i) * I_denom


      ! Examine whether layers can be split out of existence.  Stop when there
      ! is a layer that cannot be handled this way, or when the topmost formerly
      ! interior layer has been dealt with.
      do ks = nks,top_interior_ks,-1
        k = ks2(ks)
        leave_in_layer = .false.
        if ((k > nkmb) .and. (Rcv(i,k) <= RcvTgt(k))) then
          if (RcvTgt(k)-Rcv(i,k) < CS%BL_split_rho_tol*(RcvTgt(k) - RcvTgt(k-1))) &
            leave_in_layer = .true.
        elseif (k > nkmb) then
          if (Rcv(i,k)-RcvTgt(k) < CS%BL_split_rho_tol*(RcvTgt(k+1) - RcvTgt(k))) &
            leave_in_layer = .true.
        endif

        if (leave_in_layer) then
          ! Just drop this layer from the sorted list.
          nks = nks-1
        elseif (Rcv(i,k) < Rcv_int) then
          ! There is no interior layer with a target density that is intermediate
          ! between this layer and the mixed layer.
          exit
        else
          ! Try splitting the layer between two interior isopycnal layers.
          ! Find the target densities that bracket this layer.
          do k2=k_int_top+1,nz ; if (Rcv(i,k) < RcvTgt(k2)) exit ; enddo
          if (k2>nz) exit

          ! This layer is bracketed in density between layers k2-1 and k2.

          dR1 = (RcvTgt(k2-1) - Rcv(i,k)) ; dR2 = (RcvTgt(k2) - Rcv(i,k))
          T_up = T(i,k) + dT_dR * dR1
          S_up = S(i,k) + dS_dR * dR1
          T_dn = T(i,k) + dT_dR * dR2
          S_dn = S(i,k) + dS_dR * dR2

          R0_up = R0(i,k) + (dT_dR*dR0_dT(i) + dS_dR*dR0_dS(i)) * dR1
          R0_dn = R0(i,k) + (dT_dR*dR0_dT(i) + dS_dR*dR0_dS(i)) * dR2

          ! Make sure the new properties are acceptable.
          if ((R0_up > R0(i,0)) .or. (R0_dn > R0(i,0))) &
            ! Avoid creating obviously unstable profiles.
            exit

          wt_dn = (Rcv(i,k) - RcvTgt(k2-1)) / (RcvTgt(k2) - RcvTgt(k2-1))
          h_to_up = (h(i,k)-eps(i,k)) * (1.0 - wt_dn)
          h_to_dn = (h(i,k)-eps(i,k)) * wt_dn

          I_hup = 1.0 / (h(i,k2-1) + h_to_up)
          I_hdn = 1.0 / (h(i,k2) + h_to_dn)
          R0(i,k2-1) = (R0(i,k2)*h(i,k2-1) + R0_up*h_to_up) * I_hup
          R0(i,k2) = (R0(i,k2)*h(i,k2) + R0_dn*h_to_dn) * I_hdn

          T(i,k2-1) = (T(i,k2)*h(i,k2-1) + T_up*h_to_up) * I_hup
          T(i,k2) = (T(i,k2)*h(i,k2) + T_dn*h_to_dn) * I_hdn
          S(i,k2-1) = (S(i,k2)*h(i,k2-1) + S_up*h_to_up) * I_hup
          S(i,k2) = (S(i,k2)*h(i,k2) + S_dn*h_to_dn) * I_hdn
          Rcv(i,k2-1) = (Rcv(i,k2)*h(i,k2-1) + RcvTgt(k2-1)*h_to_up) * I_hup
          Rcv(i,k2) = (Rcv(i,k2)*h(i,k2) + RcvTgt(k2)*h_to_dn) * I_hdn

          h(i,k) = eps(i,k)
          h(i,k2) = h(i,k2) + h_to_dn
          h(i,k2-1) = h(i,k2-1) + h_to_up

          if (k > k2-1) then
            d_eb(i,k) = d_eb(i,k) - h_to_up
            d_eb(i,k2-1) = d_eb(i,k2-1) + h_to_up
          elseif (k < k2-1) then
            d_ea(i,k) = d_ea(i,k) - h_to_up
            d_ea(i,k2-1) = d_ea(i,k2-1) + h_to_up
          endif
          if (k > k2) then
            d_eb(i,k) = d_eb(i,k) - h_to_dn
            d_eb(i,k2) = d_eb(i,k2) + h_to_dn
          elseif (k < k2) then
            d_ea(i,k) = d_ea(i,k) - h_to_dn
            d_ea(i,k2) = d_ea(i,k2) + h_to_dn
          endif
          nks = nks-1
        endif
      enddo
    endif

    do while (nks > nkmb)
      ! Having already tried to move surface layers into the interior, there
      ! are still too many layers, and layers must be merged until nks=nkmb.
      ! Examine every merger of a layer with its neighbors, and merge the ones
      ! that increase the potential energy the least.  If there are layers
      ! with (apparently?) negative potential energy change, choose the one
      ! with the smallest total thickness.  Repeat until nkmb layers remain.
      ! Choose the smaller value for the remaining index for convenience.

      ks_min = -1 ; min_dPE = 1.0 ; min_hmin = 0.0
      do ks=1,nks-1
        k1 = ks2(ks) ; k2 = ks2(ks+1)
        dPE = max(0.0, (R0(i,k2)-R0(i,k1)) * h(i,k1) * h(i,k2))
        hmin = min(h(i,k1)-eps(i,k1), h(i,k2)-eps(i,k2))
        if ((ks_min < 0) .or. (dPE < min_dPE) .or. &
            ((dPE <= 0.0) .and. (hmin < min_hmin))) then
           ks_min = ks ; min_dPE = dPE ; min_hmin = hmin
        endif
      enddo

      ! Now merge the two layers that do the least damage.
      k1 = ks2(ks_min) ; k2 = ks2(ks_min+1)
      if (k1 < k2) then ; k_tgt = k1 ; k_src = k2
      else ; k_tgt = k2 ; k_src = k1 ; ks2(ks_min) = k2 ; endif

      h_tgt_old = h(i,k_tgt)
      h_move = h(i,k_src)-eps(i,k_src)
      h(i,k_src) = eps(i,k_src)
      h(i,k_tgt) = h(i,k_tgt) + h_move
      I_hnew = 1.0 / (h(i,k_tgt))
      R0(i,k_tgt) = (R0(i,k_tgt)*h_tgt_old + R0(i,k_src)*h_move) * I_hnew

      T(i,k_tgt) = (T(i,k_tgt)*h_tgt_old + T(i,k_src)*h_move) * I_hnew
      S(i,k_tgt) = (S(i,k_tgt)*h_tgt_old + S(i,k_src)*h_move) * I_hnew
      Rcv(i,k_tgt) = (Rcv(i,k_tgt)*h_tgt_old + Rcv(i,k_src)*h_move) * I_hnew

      d_eb(i,k_src) = d_eb(i,k_src) - h_move
      d_eb(i,k_tgt) = d_eb(i,k_tgt) + h_move

      ! Remove the newly missing layer from the sorted list.
      do ks=ks_min+1,nks ; ks2(ks) = ks2(ks+1) ; enddo
      nks = nks-1
    enddo

    !   Check again whether the layers are sorted, and go on to the next column
    ! if they are.
    sorted = .true.
    do ks=1,nks-1 ; if (ks2(ks) > ks2(ks+1)) sorted = .false. ; enddo
    if (sorted) cycle

    if (nks > 1) then
      ! At this point, actually swap the properties of the layers, and place
      ! the remaining layers in order starting with nkmb.

      ! Save all the properties of the nkmb layers that might be replaced.
      do k=1,nkmb
        h_tmp(k) = h(i,k) ; R0_tmp(k) = R0(i,k)
        T_tmp(k) = T(i,k) ; S_tmp(k) = S(i,k) ; Rcv_tmp(k) = Rcv(i,k)

        h(i,k) = 0.0
      enddo

      do ks=nks,1,-1
        k_tgt = nkmb - nks + ks ; k_src = ks2(ks)
        if (k_tgt == k_src) then
          h(i,k_tgt) = h_tmp(k_tgt)  ! This layer doesn't move, so put the water back.
          cycle
        endif

        ! Note below that eps=0 for k<=nkmb.
        if (k_src > nkmb) then
          h_move = h(i,k_src)-eps(i,k_src)
          h(i,k_src) = eps(i,k_src)
          h(i,k_tgt) = h_move
          R0(i,k_tgt) = R0(i,k_src)

          T(i,k_tgt) = T(i,k_src) ; S(i,k_tgt) = S(i,k_src)
          Rcv(i,k_tgt) = Rcv(i,k_src)

          d_eb(i,k_src) = d_eb(i,k_src) - h_move
          d_eb(i,k_tgt) = d_eb(i,k_tgt) + h_move
        else
          h(i,k_tgt) = h_tmp(k_src)
          R0(i,k_tgt) = R0_tmp(k_src)

          T(i,k_tgt) = T_tmp(k_src) ; S(i,k_tgt) = S_tmp(k_src)
          Rcv(i,k_tgt) = Rcv_tmp(k_src)

          if (k_src > k_tgt) then
            d_eb(i,k_src) = d_eb(i,k_src) - h_tmp(k_src)
            d_eb(i,k_tgt) = d_eb(i,k_tgt) + h_tmp(k_src)
          else
            d_ea(i,k_src) = d_ea(i,k_src) - h_tmp(k_src)
            d_ea(i,k_tgt) = d_ea(i,k_tgt) + h_tmp(k_src)
          endif
        endif
      enddo
    endif

  endif ; enddo

end subroutine resort_ML

!> This subroutine moves any water left in the former mixed layers into the
!! two buffer layers and may also move buffer layer water into the interior
!! isopycnal layers.
subroutine mixedlayer_detrain_2(h, T, S, R0, Rcv, RcvTgt, dt, dt_diag, d_ea, j, G, GV, US, CS, &
                                dR0_dT, dR0_dS, dRcv_dT, dRcv_dS, max_BL_det)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: h    !< Layer thickness [H ~> m or kg m-2].
                                                            !!  Layer 0 is the new mixed layer.
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: T    !< Potential temperature [degC].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: S    !< Salinity [ppt].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: R0   !< Potential density referenced to
                                                            !! surface pressure [R ~> kg m-3].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: Rcv  !< The coordinate defining potential
                                                            !! density [R ~> kg m-3].
  real, dimension(SZK_(GV)),          intent(in)    :: RcvTgt  !< The target value of Rcv for each
                                                            !! layer [R ~> kg m-3].
  real,                               intent(in)    :: dt   !< Time increment [T ~> s].
  real,                               intent(in)    :: dt_diag !< The diagnostic time step [T ~> s].
  real, dimension(SZI_(G),SZK_(GV)),  intent(inout) :: d_ea !< The upward increase across a layer in
                                                            !! the entrainment from above
                                                            !! [H ~> m or kg m-2]. Positive d_ea
                                                            !! goes with layer thickness increases.
  integer,                            intent(in)    :: j    !< The meridional row to work on.
  type(unit_scale_type),              intent(in)    :: US   !< A dimensional unit scaling type
  type(bulkmixedlayer_CS),            pointer       :: CS   !< The control structure returned by a
                                                            !! previous call to mixedlayer_init.
  real, dimension(SZI_(G)),           intent(in)    :: dR0_dT  !< The partial derivative of
                                                            !! potential density referenced to the
                                                            !! surface with potential temperature,
                                                            !! [R degC-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)),           intent(in)    :: dR0_dS  !< The partial derivative of
                                                            !! cpotential density referenced to the
                                                            !! surface with salinity
                                                            !! [R ppt-1 ~> kg m-3 ppt-1].
  real, dimension(SZI_(G)),           intent(in)    :: dRcv_dT !< The partial derivative of
                                                            !! coordinate defining potential density
                                                            !! with potential temperature,
                                                            !! [R degC-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)),           intent(in)    :: dRcv_dS !< The partial derivative of
                                                            !! coordinate defining potential density
                                                            !! with salinity [R ppt-1 ~> kg m-3 ppt-1].
  real, dimension(SZI_(G)),           intent(in)    :: max_BL_det !< If non-negative, the maximum
                                                            !! detrainment permitted from the buffer
                                                            !! layers [H ~> m or kg m-2].

! This subroutine moves any water left in the former mixed layers into the
! two buffer layers and may also move buffer layer water into the interior
! isopycnal layers.

  ! Local variables
  real :: h_to_bl                 ! The total thickness detrained to the buffer
                                  ! layers [H ~> m or kg m-2].
  real :: R0_to_bl                ! The depth integrated amount of R0 that is detrained to the
                                  ! buffer layer [H R ~> kg m-2 or kg2 m-5]
  real :: Rcv_to_bl               ! The depth integrated amount of Rcv that is detrained to the
                                  ! buffer layer [H R ~> kg m-2 or kg2 m-5]
  real :: T_to_bl                 ! The depth integrated amount of T that is detrained to the
                                  ! buffer layer [degC H ~> degC m or degC kg m-2]
  real :: S_to_bl                 ! The depth integrated amount of S that is detrained to the
                                  ! buffer layer [ppt H ~> ppt m or ppt kg m-2]
  real :: h_min_bl                ! The minimum buffer layer thickness [H ~> m or kg m-2].

  real :: h1, h2                  ! Scalar variables holding the values of
                                  ! h(i,CS%nkml+1) and h(i,CS%nkml+2) [H ~> m or kg m-2].
  real :: h1_avail                ! The thickness of the upper buffer layer
                                  ! available to move into the lower buffer
                                  ! layer [H ~> m or kg m-2].
  real :: stays                   ! stays is the thickness of the upper buffer
                                  ! layer that remains there [H ~> m or kg m-2].
  real :: stays_min, stays_max    ! The minimum and maximum permitted values of
                                  ! stays [H ~> m or kg m-2].

  logical :: mergeable_bl         ! If true, it is an option to combine the two
                                  ! buffer layers and create water that matches
                                  ! the target density of an interior layer.
  real :: stays_merge             ! If the two buffer layers can be combined
                                  ! stays_merge is the thickness of the upper
                                  ! layer that remains [H ~> m or kg m-2].
  real :: stays_min_merge         ! The minimum allowed value of stays_merge [H ~> m or kg m-2].

  real :: dR0_2dz, dRcv_2dz       ! Half the vertical gradients of R0 and Rcv [R H-1 ~> kg m-4 or m-1]
!  real :: dT_2dz, dS_2dz         ! Half the vertical gradients of T and S, in degC H-1, and ppt H-1.
  real :: scale_slope             ! A nondimensional number < 1 used to scale down
                                  ! the slope within the upper buffer layer when
                                  ! water MUST be detrained to the lower layer.

  real :: dPE_extrap              ! The potential energy change due to dispersive
                                  ! advection or mixing layers, divided by
                                  ! rho_0*g [H2 ~> m2 or kg2 m-4].
  real :: dPE_det, dPE_merge      ! The energy required to mix the detrained water
                                  ! into the buffer layer or the merge the two
                                  ! buffer layers [R H2 L2 Z-1 T-2 ~> J m-2 or J kg2 m-8].

  real :: h_from_ml               ! The amount of additional water that must be
                                  ! drawn from the mixed layer [H ~> m or kg m-2].
  real :: h_det_h2                ! The amount of detrained water and mixed layer
                                  ! water that will go directly into the lower
                                  ! buffer layer [H ~> m or kg m-2].
  real :: h_det_to_h2, h_ml_to_h2 ! All of the variables hA_to_hB are the thickness fluxes
  real :: h_det_to_h1, h_ml_to_h1 ! from one layer to another [H ~> m or kg m-2],
  real :: h1_to_h2, h1_to_k0      ! with h_det the detrained water, h_ml
  real :: h2_to_k1, h2_to_k1_rem  ! the actively mixed layer, h1 and h2 the upper
                                  ! and lower buffer layers, and k0 and k1 the
                                  ! interior layers that are just lighter and
                                  ! just denser than the lower buffer layer.

  real :: R0_det, T_det, S_det    ! Detrained values of R0 [R ~> kg m-3], T [degC], and S [ppt].
  real :: Rcv_stays, R0_stays     ! Values of Rcv and R0 that stay in a layer.
  real :: T_stays, S_stays        ! Values of T and S that stay in a layer.
  real :: dSpice_det, dSpice_stays! The spiciness difference between an original
                                  ! buffer layer and the water that moves into
                                  ! an interior layer or that stays in that
                                  ! layer [R ~> kg m-3].
  real :: dSpice_lim, dSpice_lim2 ! Limits to the spiciness difference between
                                  ! the lower buffer layer and the water that
                                  ! moves into an interior layer [R ~> kg m-3].
  real :: dSpice_2dz              ! The vertical gradient of spiciness used for
                                  ! advection [R H-1 ~> kg m-4 or m-1].

  real :: dPE_ratio               ! Multiplier of dPE_det at which merging is
                                  ! permitted - here (detrainment_per_day/dt)*30
                                  ! days?
  real :: num_events              ! The number of detrainment events over which
                                  ! to prefer merging the buffer layers.
  real :: dPE_time_ratio          ! Larger of 1 and the detrainment timescale over dt [nondim].
  real :: dT_dS_gauge, dS_dT_gauge ! The relative scales of temperature and
                                  ! salinity changes in defining spiciness, in
                                  ! [degC ppt-1] and [ppt degC-1].
  real :: I_denom                 ! A work variable with units of [ppt2 R-2 ~> ppt2 m6 kg-2].

  real :: g_2                     ! 1/2 g_Earth [L2 Z-1 T-2 ~> m s-2].
  real :: Rho0xG                  ! Rho0 times G_Earth [R L2 Z-1 T-2 ~> kg m-2 s-2].
  real :: I2Rho0                  ! 1 / (2 Rho0) [R-1 ~> m3 kg-1].
  real :: Idt_H2                  ! The square of the conversion from thickness to Z
                                  ! divided by the time step [Z2 H-2 T-1 ~> s-1 or m6 kg-2 s-1].
  logical :: stable_Rcv           ! If true, the buffer layers are stable with
                                  ! respect to the coordinate potential density.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].

  real :: s1en                    ! A work variable [H2 L2 kg m-1 T-3 ~> kg m3 s-3 or kg3 m-3 s-3].
  real :: s1, s2, bh0             ! Work variables [H ~> m or kg m-2].
  real :: s3sq                    ! A work variable [H2 ~> m2 or kg2 m-4].
  real :: I_ya, b1                ! Nondimensional work variables.
  real :: Ih, Ihdet, Ih1f, Ih2f   ! Assorted inverse thickness work variables,
  real :: Ihk0, Ihk1, Ih12        ! all in [H-1 ~> m-1 or m2 kg-1].
  real :: dR1, dR2, dR2b, dRk1    ! Assorted density difference work variables,
  real :: dR0, dR21, dRcv         ! all in [R ~> kg m-3].
  real :: dRcv_stays, dRcv_det, dRcv_lim
  real :: Angstrom                ! The minumum layer thickness [H ~> m or kg m-2].

  real :: h2_to_k1_lim, T_new, S_new, T_max, T_min, S_max, S_min
  character(len=200) :: mesg

  integer :: i, k, k0, k1, is, ie, nz, kb1, kb2, nkmb
  is = G%isc ; ie = G%iec ; nz = GV%ke
  kb1 = CS%nkml+1; kb2 = CS%nkml+2
  nkmb = CS%nkml+CS%nkbl
  h_neglect = GV%H_subroundoff
  g_2 = 0.5 * GV%g_Earth
  Rho0xG = GV%Rho0 * GV%g_Earth
  Idt_H2 = GV%H_to_Z**2 / dt_diag
  I2Rho0 = 0.5 / (GV%Rho0)
  Angstrom = GV%Angstrom_H

  ! This is hard coding of arbitrary and dimensional numbers.
  dT_dS_gauge = CS%dT_dS_wt ; dS_dT_gauge = 1.0 / dT_dS_gauge
  num_events = 10.0

  if (CS%nkbl /= 2) call MOM_error(FATAL, "MOM_mixed_layer"// &
                        "CS%nkbl must be 2 in mixedlayer_detrain_2.")

  if (dt < CS%BL_detrain_time) then ; dPE_time_ratio = CS%BL_detrain_time / (dt)
  else ; dPE_time_ratio = 1.0 ; endif

  do i=is,ie

  ! Determine all of the properties being detrained from the mixed layer.

  ! As coded this has the k and i loop orders switched, but k is CS%nkml is
  ! often just 1 or 2, so this seems like it should not be a problem, especially
  ! since it means that a number of variables can now be scalars, not arrays.
    h_to_bl = 0.0 ; R0_to_bl = 0.0
    Rcv_to_bl = 0.0 ; T_to_bl = 0.0 ; S_to_bl = 0.0

    do k=1,CS%nkml ; if (h(i,k) > 0.0) then
      h_to_bl = h_to_bl + h(i,k)
      R0_to_bl = R0_to_bl + R0(i,k)*h(i,k)

      Rcv_to_bl = Rcv_to_bl + Rcv(i,k)*h(i,k)
      T_to_bl = T_to_bl + T(i,k)*h(i,k)
      S_to_bl = S_to_bl + S(i,k)*h(i,k)

      d_ea(i,k) = d_ea(i,k) - h(i,k)
      h(i,k) = 0.0
    endif ; enddo
    if (h_to_bl > 0.0) then ; R0_det = R0_to_bl / h_to_bl
    else ; R0_det = R0(i,0) ; endif

    ! This code does both downward detrainment from both the mixed layer and the
    ! buffer layers.
    !   Several considerations apply in detraining water into the interior:
    ! (1) Water only moves into the interior from the deeper buffer layer,
    !     so the deeper buffer layer must have some mass.
    ! (2) The upper buffer layer must have some mass so the extrapolation of
    !     density is meaningful (i.e. there is not detrainment from the buffer
    !     layers when there is strong mixed layer entrainment).
    ! (3) The lower buffer layer density extrapolated to its base with a
    !     linear fit between the two layers must exceed the density of the
    !     next denser interior layer.
    ! (4) The average extroplated coordinate density that is moved into the
    !     isopycnal interior matches the target value for that layer.
    ! (5) The potential energy change is calculated and might be used later
    !     to allow the upper buffer layer to mix more into the lower buffer
    !     layer.

    ! Determine whether more must be detrained from the mixed layer to keep a
    ! minimal amount of mass in the buffer layers.  In this case the 5% of the
    ! mixed layer thickness is hard-coded, but probably shouldn't be!
    h_min_bl = MIN(CS%Hbuffer_min, CS%Hbuffer_rel_min*h(i,0))

    stable_Rcv = .true.
    if (((R0(i,kb2)-R0(i,kb1)) * (Rcv(i,kb2)-Rcv(i,kb1)) <= 0.0)) &
      stable_Rcv = .false.

    h1 = h(i,kb1) ; h2 = h(i,kb2)

    h2_to_k1_rem = (h1 + h2) + h_to_bl
    if ((max_BL_det(i) >= 0.0) .and. (h2_to_k1_rem > max_BL_det(i))) &
      h2_to_k1_rem = max_BL_det(i)


    if ((h2 == 0.0) .and. (h1 > 0.0)) then
      ! The lower buffer layer has been eliminated either by convective
      ! adjustment or entrainment from the interior, and its current properties
      ! are not meaningful, but may later be used to determine the properties of
      ! waters moving into the lower buffer layer.  So the properties of the
      ! lower buffer layer are set to be between those of the upper buffer layer
      ! and the next denser interior layer, measured by R0.  This probably does
      ! not happen very often, so I am not too worried about the inefficiency of
      ! the following loop.
      do k1=kb2+1,nz ; if (h(i,k1) > 2.0*Angstrom) exit ; enddo

      R0(i,kb2) = R0(i,kb1)

      Rcv(i,kb2)=Rcv(i,kb1) ; T(i,kb2)=T(i,kb1) ; S(i,kb2)=S(i,kb1)


      if (k1 <= nz) then ; if (R0(i,k1) >= R0(i,kb1)) then
        R0(i,kb2) = 0.5*(R0(i,kb1)+R0(i,k1))

        Rcv(i,kb2) = 0.5*(Rcv(i,kb1)+Rcv(i,k1))
        T(i,kb2) = 0.5*(T(i,kb1)+T(i,k1))
        S(i,kb2) = 0.5*(S(i,kb1)+S(i,k1))
      endif ; endif
    endif ! (h2 = 0 && h1 > 0)

    dPE_extrap = 0.0 ; dPE_merge = 0.0
    mergeable_bl = .false.
    if ((h1 > 0.0) .and. (h2 > 0.0) .and. (h_to_bl > 0.0) .and. &
        (stable_Rcv)) then
      ! Check whether it is permissible for the buffer layers to detrain
      ! into the interior isopycnal layers.

      ! Determine the layer that has the lightest target density that is
      ! denser than the lowermost buffer layer.
      do k1=kb2+1,nz ; if (RcvTgt(k1) >= Rcv(i,kb2)) exit ; enddo ; k0 = k1-1
      dR1 = RcvTgt(k0)-Rcv(i,kb1) ; dR2 = Rcv(i,kb2)-RcvTgt(k0)

      ! Use an energy-balanced combination of downwind advection into the next
      ! denser interior layer and upwind advection from the upper buffer layer
      ! into the lower one, each with an energy change that equals that required
      ! to mix the detrained water with the upper buffer layer.
      h1_avail = h1 - MAX(0.0,h_min_bl-h_to_bl)
      if ((k1<=nz) .and. (h2 > h_min_bl) .and. (h1_avail > 0.0) .and. &
          (R0(i,kb1) < R0(i,kb2)) .and. (h_to_bl*R0(i,kb1) > R0_to_bl)) then
        dRk1 = (RcvTgt(k1) - Rcv(i,kb2)) * (R0(i,kb2) - R0(i,kb1)) / &
                                           (Rcv(i,kb2) - Rcv(i,kb1))
        b1 = dRk1 / (R0(i,kb2) - R0(i,kb1))
        ! b1 = RcvTgt(k1) - Rcv(i,kb2)) / (Rcv(i,kb2) - Rcv(i,kb1))

        ! Apply several limits to the detrainment.
        ! Entrain less than the mass in h2, and keep the base of the buffer
        ! layers from becoming shallower than any neighbors.
        h2_to_k1 = min(h2 - h_min_bl, h2_to_k1_rem)
        ! Balance downwind advection of density into the layer below the
        ! buffer layers with upwind advection from the layer above.
        if (h2_to_k1*(h1_avail + b1*(h1_avail + h2)) > h2*h1_avail) &
          h2_to_k1 = (h2*h1_avail) / (h1_avail + b1*(h1_avail + h2))
        if (h2_to_k1*(dRk1 * h2) > (h_to_bl*R0(i,kb1) - R0_to_bl) * h1) &
          h2_to_k1 = (h_to_bl*R0(i,kb1) - R0_to_bl) * h1 / (dRk1 * h2)

        if ((k1==kb2+1) .and. (CS%BL_extrap_lim > 0.)) then
          ! Simply do not detrain very light water into the lightest isopycnal
          ! coordinate layers if the density jump is too large.
          dRcv_lim = Rcv(i,kb2)-Rcv(i,0)
          do k=1,kb2 ; dRcv_lim = max(dRcv_lim, Rcv(i,kb2)-Rcv(i,k)) ; enddo
          dRcv_lim = CS%BL_extrap_lim*dRcv_lim
          if ((RcvTgt(k1) - Rcv(i,kb2)) >= dRcv_lim) then
            h2_to_k1 = 0.0
          elseif ((RcvTgt(k1) - Rcv(i,kb2)) > 0.5*dRcv_lim) then
            h2_to_k1 = h2_to_k1 * (2.0 - 2.0*((RcvTgt(k1) - Rcv(i,kb2)) / dRcv_lim))
          endif
        endif

        dRcv = (RcvTgt(k1) - Rcv(i,kb2))

        ! Use 2nd order upwind advection of spiciness, limited by the values
        ! in deeper thick layers to determine the detrained temperature and
        ! salinity.
        dSpice_det = (dS_dT_gauge*dRcv_dS(i)*(T(i,kb2)-T(i,kb1)) - &
                      dT_dS_gauge*dRcv_dT(i)*(S(i,kb2)-S(i,kb1))) * &
                      (h2 - h2_to_k1) / (h1 + h2)
        dSpice_lim = 0.0
        if (h(i,k1) > 10.0*Angstrom) then
          dSpice_lim = dS_dT_gauge*dRcv_dS(i)*(T(i,k1)-T(i,kb2)) - &
                       dT_dS_gauge*dRcv_dT(i)*(S(i,k1)-S(i,kb2))
          if (dSpice_det*dSpice_lim <= 0.0) dSpice_lim = 0.0
        endif
        if (k1<nz) then ; if (h(i,k1+1) > 10.0*Angstrom) then
          dSpice_lim2 = dS_dT_gauge*dRcv_dS(i)*(T(i,k1+1)-T(i,kb2)) - &
                        dT_dS_gauge*dRcv_dT(i)*(S(i,k1+1)-S(i,kb2))
          if ((dSpice_det*dSpice_lim2 > 0.0) .and. &
              (abs(dSpice_lim2) > abs(dSpice_lim))) dSpice_lim = dSpice_lim2
        endif ; endif
        if (abs(dSpice_det) > abs(dSpice_lim)) dSpice_det = dSpice_lim

        I_denom = 1.0 / (dRcv_dS(i)**2 + (dT_dS_gauge*dRcv_dT(i))**2)
        T_det = T(i,kb2) + dT_dS_gauge * I_denom * &
            (dT_dS_gauge * dRcv_dT(i) * dRcv + dRcv_dS(i) * dSpice_det)
        S_det = S(i,kb2) + I_denom * &
            (dRcv_dS(i) * dRcv - dT_dS_gauge * dRcv_dT(i) * dSpice_det)
        ! The detrained values of R0 are based on changes in T and S.
        R0_det = R0(i,kb2) + (T_det-T(i,kb2)) * dR0_dT(i) + &
                             (S_det-S(i,kb2)) * dR0_dS(i)

        if (CS%BL_extrap_lim >= 0.) then
          ! Only do this detrainment if the new layer's temperature and salinity
          ! are not too far outside of the range of previous values.
          if (h(i,k1) > 10.0*Angstrom) then
            T_min = min(T(i,kb1), T(i,kb2), T(i,k1)) - CS%Allowed_T_chg
            T_max = max(T(i,kb1), T(i,kb2), T(i,k1)) + CS%Allowed_T_chg
            S_min = min(S(i,kb1), S(i,kb2), S(i,k1)) - CS%Allowed_S_chg
            S_max = max(S(i,kb1), S(i,kb2), S(i,k1)) + CS%Allowed_S_chg
          else
            T_min = min(T(i,kb1), T(i,kb2)) - CS%Allowed_T_chg
            T_max = max(T(i,kb1), T(i,kb2)) + CS%Allowed_T_chg
            S_min = min(S(i,kb1), S(i,kb2)) - CS%Allowed_S_chg
            S_max = max(S(i,kb1), S(i,kb2)) + CS%Allowed_S_chg
          endif
          Ihk1 = 1.0 / (h(i,k1) + h2_to_k1)
          T_new = (h(i,k1)*T(i,k1) + h2_to_k1*T_det) * Ihk1
          S_new = (h(i,k1)*S(i,k1) + h2_to_k1*S_det) * Ihk1
          ! A less restrictive limit might be used here.
          if ((T_new < T_min) .or. (T_new > T_max) .or. &
              (S_new < S_min) .or. (S_new > S_max)) &
            h2_to_k1 = 0.0
        endif

        h1_to_h2 = b1*h2*h2_to_k1 / (h2 - (1.0+b1)*h2_to_k1)

        Ihk1 = 1.0 / (h(i,k1) + h_neglect + h2_to_k1)
        Ih2f = 1.0 / ((h(i,kb2) - h2_to_k1) + h1_to_h2)

        Rcv(i,kb2) = ((h(i,kb2)*Rcv(i,kb2) - h2_to_k1*RcvTgt(k1)) + &
                      h1_to_h2*Rcv(i,kb1))*Ih2f
        Rcv(i,k1) = ((h(i,k1)+h_neglect)*Rcv(i,k1) + h2_to_k1*RcvTgt(k1)) * Ihk1

        T(i,kb2) = ((h(i,kb2)*T(i,kb2) - h2_to_k1*T_det) + &
                    h1_to_h2*T(i,kb1)) * Ih2f
        T(i,k1) = ((h(i,k1)+h_neglect)*T(i,k1) + h2_to_k1*T_det) * Ihk1

        S(i,kb2) = ((h(i,kb2)*S(i,kb2) - h2_to_k1*S_det) + &
                    h1_to_h2*S(i,kb1)) * Ih2f
        S(i,k1) = ((h(i,k1)+h_neglect)*S(i,k1) + h2_to_k1*S_det) * Ihk1

        ! Changes in R0 are based on changes in T and S.
        R0(i,kb2) = ((h(i,kb2)*R0(i,kb2) - h2_to_k1*R0_det) + &
                     h1_to_h2*R0(i,kb1)) * Ih2f
        R0(i,k1) = ((h(i,k1)+h_neglect)*R0(i,k1) + h2_to_k1*R0_det) * Ihk1

        h(i,kb1) = h(i,kb1) - h1_to_h2 ; h1 = h(i,kb1)
        h(i,kb2) = (h(i,kb2) - h2_to_k1) + h1_to_h2 ; h2 = h(i,kb2)
        h(i,k1) = h(i,k1) + h2_to_k1

        d_ea(i,kb1) = d_ea(i,kb1) - h1_to_h2
        d_ea(i,kb2) = (d_ea(i,kb2) - h2_to_k1) + h1_to_h2
        d_ea(i,k1) = d_ea(i,k1) + h2_to_k1
        h2_to_k1_rem = max(h2_to_k1_rem - h2_to_k1, 0.0)

        !   The lower buffer layer has become lighter - it may be necessary to
        ! adjust k1 lighter.
        if ((k1>kb2+1) .and. (RcvTgt(k1-1) >= Rcv(i,kb2))) then
          do k1=k1,kb2+1,-1 ; if (RcvTgt(k1-1) < Rcv(i,kb2)) exit ; enddo
        endif
      endif

      k0 = k1-1
      dR1 = RcvTgt(k0)-Rcv(i,kb1) ; dR2 = Rcv(i,kb2)-RcvTgt(k0)

      if ((k0>kb2) .and. (dR1 > 0.0) .and. (h1 > h_min_bl) .and. &
          (h2*dR2 < h1*dR1) .and. (R0(i,kb2) > R0(i,kb1))) then
        ! An interior isopycnal layer (k0) is intermediate in density between
        ! the two buffer layers, and there can be detrainment. The entire
        ! lower buffer layer is combined with a portion of the upper buffer
        ! layer to match the target density of layer k0.
        stays_merge = 2.0*(h1+h2)*(h1*dR1 - h2*dR2) / &
                     ((dR1+dR2)*h1 + dR1*(h1+h2) + &
                      sqrt((dR2*h1-dR1*h2)**2 + 4*(h1+h2)*h2*(dR1+dR2)*dR2))

        stays_min_merge = MAX(h_min_bl, 2.0*h_min_bl - h_to_bl, &
                  h1 - (h1+h2)*(R0(i,kb1) - R0_det) / (R0(i,kb2) - R0(i,kb1)))
        if ((stays_merge > stays_min_merge) .and. &
            (stays_merge + h2_to_k1_rem >= h1 + h2)) then
          mergeable_bl = .true.
          dPE_merge = g_2*(R0(i,kb2)-R0(i,kb1))*(h1-stays_merge)*(h2-stays_merge)
        endif
      endif

      if ((k1<=nz).and.(.not.mergeable_bl)) then
        ! Check whether linear extrapolation of density (i.e. 2nd order upwind
        ! advection) will allow some of the lower buffer layer to detrain into
        ! the next denser interior layer (k1).
        dR2b = RcvTgt(k1)-Rcv(i,kb2) ; dR21 = Rcv(i,kb2) - Rcv(i,kb1)
        if (dR2b*(h1+h2) < h2*dR21) then
          ! Some of layer kb2 is denser than k1.
          h2_to_k1 = min(h2 - (h1+h2) * dR2b / dR21, h2_to_k1_rem)

          if (h2 > h2_to_k1) then
            dRcv = (RcvTgt(k1) - Rcv(i,kb2))

            ! Use 2nd order upwind advection of spiciness, limited by the values
            ! in deeper thick layers to determine the detrained temperature and
            ! salinity.
            dSpice_det = (dS_dT_gauge*dRcv_dS(i)*(T(i,kb2)-T(i,kb1)) - &
                          dT_dS_gauge*dRcv_dT(i)*(S(i,kb2)-S(i,kb1))) * &
                          (h2 - h2_to_k1) / (h1 + h2)
            dSpice_lim = 0.0
            if (h(i,k1) > 10.0*Angstrom) then
              dSpice_lim = dS_dT_gauge*dRcv_dS(i)*(T(i,k1)-T(i,kb2)) - &
                           dT_dS_gauge*dRcv_dT(i)*(S(i,k1)-S(i,kb2))
              if (dSpice_det*dSpice_lim <= 0.0) dSpice_lim = 0.0
            endif
            if (k1<nz) then; if (h(i,k1+1) > 10.0*Angstrom) then
              dSpice_lim2 = dS_dT_gauge*dRcv_dS(i)*(T(i,k1+1)-T(i,kb2)) - &
                            dT_dS_gauge*dRcv_dT(i)*(S(i,k1+1)-S(i,kb2))
              if ((dSpice_det*dSpice_lim2 > 0.0) .and. &
                  (abs(dSpice_lim2) > abs(dSpice_lim))) dSpice_lim = dSpice_lim2
            endif; endif
            if (abs(dSpice_det) > abs(dSpice_lim)) dSpice_det = dSpice_lim

            I_denom = 1.0 / (dRcv_dS(i)**2 + (dT_dS_gauge*dRcv_dT(i))**2)
            T_det = T(i,kb2) + dT_dS_gauge * I_denom * &
                (dT_dS_gauge * dRcv_dT(i) * dRcv + dRcv_dS(i) * dSpice_det)
            S_det = S(i,kb2) + I_denom * &
                (dRcv_dS(i) * dRcv - dT_dS_gauge * dRcv_dT(i) * dSpice_det)
            ! The detrained values of R0 are based on changes in T and S.
            R0_det = R0(i,kb2) + (T_det-T(i,kb2)) * dR0_dT(i) + &
                                 (S_det-S(i,kb2)) * dR0_dS(i)

            ! Now that the properties of the detrained water are known,
            ! potentially limit the amount of water that is detrained to
            ! avoid creating unphysical properties in the remaining water.
            Ih2f = 1.0 / (h2 - h2_to_k1)

            T_min = min(T(i,kb2), T(i,kb1)) - CS%Allowed_T_chg
            T_max = max(T(i,kb2), T(i,kb1)) + CS%Allowed_T_chg
            T_new = (h2*T(i,kb2) - h2_to_k1*T_det)*Ih2f
            if (T_new < T_min) then
              h2_to_k1_lim = h2 * (T(i,kb2) - T_min) / (T_det - T_min)
!               write(mesg,'("Low temperature limits det to ", &
!                    & 1pe12.5, " from ", 1pe12.5, " at ", 1pg11.4,"E, ",1pg11.4,"N. T=", &
!                    & 5(1pe12.5))') &
!                    h2_to_k1_lim, h2_to_k1, G%geoLonT(i,j), G%geoLatT(i,j), &
!                    T_new, T(i,kb2), T(i,kb1), T_det, T_new-T_min
!               call MOM_error(WARNING, mesg)
              h2_to_k1 = h2_to_k1_lim
              Ih2f = 1.0 / (h2 - h2_to_k1)
            elseif (T_new > T_max) then
              h2_to_k1_lim = h2 * (T(i,kb2) - T_max) / (T_det - T_max)
!               write(mesg,'("High temperature limits det to ", &
!                    & 1pe12.5, " from ", 1pe12.5, " at ", 1pg11.4,"E, ",1pg11.4,"N. T=", &
!                    & 5(1pe12.5))') &
!                    h2_to_k1_lim, h2_to_k1, G%geoLonT(i,j), G%geoLatT(i,j), &
!                    T_new, T(i,kb2), T(i,kb1), T_det, T_new-T_max
!               call MOM_error(WARNING, mesg)
              h2_to_k1 = h2_to_k1_lim
              Ih2f = 1.0 / (h2 - h2_to_k1)
            endif
            S_min = max(min(S(i,kb2), S(i,kb1)) - CS%Allowed_S_chg, 0.0)
            S_max = max(S(i,kb2), S(i,kb1)) + CS%Allowed_S_chg
            S_new = (h2*S(i,kb2) - h2_to_k1*S_det)*Ih2f
            if (S_new < S_min) then
              h2_to_k1_lim = h2 * (S(i,kb2) - S_min) / (S_det - S_min)
!               write(mesg,'("Low salinity limits det to ", &
!                    & 1pe12.5, " from ", 1pe12.5, " at ", 1pg11.4,"E, ",1pg11.4,"N. S=", &
!                    & 5(1pe12.5))') &
!                    h2_to_k1_lim, h2_to_k1, G%geoLonT(i,j), G%geoLatT(i,j), &
!                    S_new, S(i,kb2), S(i,kb1), S_det, S_new-S_min
!               call MOM_error(WARNING, mesg)
              h2_to_k1 = h2_to_k1_lim
              Ih2f = 1.0 / (h2 - h2_to_k1)
            elseif (S_new > S_max) then
              h2_to_k1_lim = h2 * (S(i,kb2) - S_max) / (S_det - S_max)
!               write(mesg,'("High salinity limits det to ", &
!                    & 1pe12.5, " from ", 1pe12.5, " at ", 1pg11.4,"E, ",1pg11.4,"N. S=", &
!                    & 5(1pe12.5))') &
!                    h2_to_k1_lim, h2_to_k1, G%geoLonT(i,j), G%geoLatT(i,j), &
!                    S_new, S(i,kb2), S(i,kb1), S_det, S_new-S_max
!               call MOM_error(WARNING, mesg)
              h2_to_k1 = h2_to_k1_lim
              Ih2f = 1.0 / (h2 - h2_to_k1)
            endif

            Ihk1 = 1.0 / (h(i,k1) + h_neglect + h2_to_k1)
            Rcv(i,k1) = ((h(i,k1)+h_neglect)*Rcv(i,k1) + h2_to_k1*RcvTgt(k1)) * Ihk1
            Rcv(i,kb2) = Rcv(i,kb2) - h2_to_k1*dRcv*Ih2f

            T(i,kb2) = (h2*T(i,kb2) - h2_to_k1*T_det)*Ih2f
            T(i,k1) = ((h(i,k1)+h_neglect)*T(i,k1) + h2_to_k1*T_det) * Ihk1

            S(i,kb2) = (h2*S(i,kb2) - h2_to_k1*S_det) * Ih2f
            S(i,k1) = ((h(i,k1)+h_neglect)*S(i,k1) + h2_to_k1*S_det) * Ihk1

            ! Changes in R0 are based on changes in T and S.
            R0(i,kb2) = (h2*R0(i,kb2) - h2_to_k1*R0_det) * Ih2f
            R0(i,k1) = ((h(i,k1)+h_neglect)*R0(i,k1) + h2_to_k1*R0_det) * Ihk1
          else
            ! h2==h2_to_k1 can happen if dR2b = 0 exactly, but this is very
            ! unlikely.  In this case the entirety of layer kb2 is detrained.
            h2_to_k1 = h2  ! These 2 lines are probably unnecessary.
            Ihk1 = 1.0 / (h(i,k1) + h2)

            Rcv(i,k1) = (h(i,k1)*Rcv(i,k1) + h2*Rcv(i,kb2)) * Ihk1
            T(i,k1) = (h(i,k1)*T(i,k1) + h2*T(i,kb2)) * Ihk1
            S(i,k1) = (h(i,k1)*S(i,k1) + h2*S(i,kb2)) * Ihk1
            R0(i,k1) = (h(i,k1)*R0(i,k1) + h2*R0(i,kb2)) * Ihk1
          endif

          h(i,k1) = h(i,k1) + h2_to_k1
          h(i,kb2) = h(i,kb2) - h2_to_k1 ; h2 = h(i,kb2)
          ! dPE_extrap should be positive here.
          dPE_extrap = I2Rho0*(R0_det-R0(i,kb2))*h2_to_k1*h2

          d_ea(i,kb2) = d_ea(i,kb2) - h2_to_k1
          d_ea(i,k1) = d_ea(i,k1) + h2_to_k1
          h2_to_k1_rem = max(h2_to_k1_rem - h2_to_k1, 0.0)
        endif
      endif ! Detrainment by extrapolation.

    endif ! Detrainment to the interior at all.

    ! Does some of the detrained water go into the lower buffer layer?
    h_det_h2 = MAX(h_min_bl-(h1+h2), 0.0)
    if (h_det_h2 > 0.0) then
      ! Detrained water will go into both upper and lower buffer layers.
      ! h(kb2) will be h_min_bl, but h(kb1) may be larger if there was already
      ! ample detrainment; all water in layer kb1 moves into layer kb2.

      ! Determine the fluxes between the various layers.
      h_det_to_h2 = MIN(h_to_bl, h_det_h2)
      h_ml_to_h2 = h_det_h2 - h_det_to_h2
      h_det_to_h1 = h_to_bl - h_det_to_h2
      h_ml_to_h1 = MAX(h_min_bl-h_det_to_h1,0.0)

      Ih = 1.0/h_min_bl
      Ihdet = 0.0 ; if (h_to_bl > 0.0) Ihdet = 1.0 / h_to_bl
      Ih1f = 1.0 / (h_det_to_h1 + h_ml_to_h1)

      R0(i,kb2) = ((h2*R0(i,kb2) + h1*R0(i,kb1)) + &
                   (h_det_to_h2*R0_to_bl*Ihdet + h_ml_to_h2*R0(i,0))) * Ih
      R0(i,kb1) = (h_det_to_h1*R0_to_bl*Ihdet + h_ml_to_h1*R0(i,0)) * Ih1f

      Rcv(i,kb2) = ((h2*Rcv(i,kb2) + h1*Rcv(i,kb1)) + &
                    (h_det_to_h2*Rcv_to_bl*Ihdet + h_ml_to_h2*Rcv(i,0))) * Ih
      Rcv(i,kb1) = (h_det_to_h1*Rcv_to_bl*Ihdet + h_ml_to_h1*Rcv(i,0)) * Ih1f

      T(i,kb2) = ((h2*T(i,kb2) + h1*T(i,kb1)) + &
                  (h_det_to_h2*T_to_bl*Ihdet + h_ml_to_h2*T(i,0))) * Ih
      T(i,kb1) = (h_det_to_h1*T_to_bl*Ihdet + h_ml_to_h1*T(i,0)) * Ih1f

      S(i,kb2) = ((h2*S(i,kb2) + h1*S(i,kb1)) + &
                  (h_det_to_h2*S_to_bl*Ihdet + h_ml_to_h2*S(i,0))) * Ih
      S(i,kb1) = (h_det_to_h1*S_to_bl*Ihdet + h_ml_to_h1*S(i,0)) * Ih1f

      ! Recall that h1 = h(i,kb1) & h2 = h(i,kb2).
      d_ea(i,1) = d_ea(i,1) - (h_ml_to_h1 + h_ml_to_h2)
      d_ea(i,kb1) = d_ea(i,kb1) + ((h_det_to_h1 + h_ml_to_h1) - h1)
      d_ea(i,kb2) = d_ea(i,kb2) + (h_min_bl - h2)

      h(i,kb1) = h_det_to_h1 + h_ml_to_h1 ; h(i,kb2) = h_min_bl
      h(i,0) = h(i,0) - (h_ml_to_h1 + h_ml_to_h2)


      if (allocated(CS%diag_PE_detrain) .or. allocated(CS%diag_PE_detrain2)) then
        R0_det = R0_to_bl*Ihdet
        s1en = g_2 * Idt_H2 * ( ((R0(i,kb2)-R0(i,kb1))*h1*h2 + &
            h_det_to_h2*( (R0(i,kb1)-R0_det)*h1 + (R0(i,kb2)-R0_det)*h2 ) + &
            h_ml_to_h2*( (R0(i,kb2)-R0(i,0))*h2 + (R0(i,kb1)-R0(i,0))*h1 + &
                         (R0_det-R0(i,0))*h_det_to_h2 ) + &
            h_det_to_h1*h_ml_to_h1*(R0_det-R0(i,0))) - 2.0*GV%Rho0*dPE_extrap )

        if (allocated(CS%diag_PE_detrain)) &
          CS%diag_PE_detrain(i,j) = CS%diag_PE_detrain(i,j) + s1en

        if (allocated(CS%diag_PE_detrain2)) CS%diag_PE_detrain2(i,j) = &
            CS%diag_PE_detrain2(i,j) + s1en + Idt_H2*Rho0xG*dPE_extrap
      endif

    elseif ((h_to_bl > 0.0) .or. (h1 < h_min_bl) .or. (h2 < h_min_bl)) then
    ! Determine how much of the upper buffer layer will be moved into
    ! the lower buffer layer and the properties with which it is moving.
    ! This implementation assumes a 2nd-order upwind advection of density
    ! from the uppermost buffer layer into the next one down.
      h_from_ml = h_min_bl + MAX(h_min_bl-h2,0.0) - h1 - h_to_bl
      if (h_from_ml > 0.0) then
        ! Some water needs to be moved from the mixed layer so that the upper
        ! (and perhaps lower) buffer layers exceed their minimum thicknesses.
        dPE_extrap = dPE_extrap - I2Rho0*h_from_ml*(R0_to_bl - R0(i,0)*h_to_bl)
        R0_to_bl = R0_to_bl + h_from_ml*R0(i,0)
        Rcv_to_bl = Rcv_to_bl + h_from_ml*Rcv(i,0)
        T_to_bl = T_to_bl + h_from_ml*T(i,0)
        S_to_bl = S_to_bl + h_from_ml*S(i,0)

        h_to_bl = h_to_bl + h_from_ml
        h(i,0) = h(i,0) - h_from_ml
        d_ea(i,1) = d_ea(i,1) - h_from_ml
      endif

      ! The absolute value should be unnecessary and 1e9 is just a large number.
      b1 = 1.0e9
      if (R0(i,kb2) - R0(i,kb1) > 1.0e-9*abs(R0(i,kb1) - R0_det)) &
        b1 = abs(R0(i,kb1) - R0_det) / (R0(i,kb2) - R0(i,kb1))
      stays_min = MAX((1.0-b1)*h1 - b1*h2, 0.0, h_min_bl - h_to_bl)
      stays_max = h1 - MAX(h_min_bl-h2,0.0)

      scale_slope = 1.0
      if (stays_max <= stays_min) then
        stays = stays_max
        mergeable_bl = .false.
        if (stays_max < h1) scale_slope = (h1 - stays_min) / (h1 - stays_max)
      else
        ! There are numerous temporary variables used here that should not be
        ! used outside of this "else" branch: s1, s2, s3sq, I_ya, bh0
        bh0 = b1*h_to_bl
        I_ya =  (h1 + h2) / ((h1 + h2) + h_to_bl)
        ! s1 is the amount staying that minimizes the PE increase.
        s1 = 0.5*(h1 + (h2 - bh0) * I_ya) ; s2 = h1 - s1

        if (s2 < 0.0) then
          ! The energy released by detrainment from the lower buffer layer can be
          ! used to mix water from the upper buffer layer into the lower one.
          s3sq = I_ya*MAX(bh0*h1-dPE_extrap, 0.0)
        else
          s3sq = I_ya*(bh0*h1-MIN(dPE_extrap,0.0))
        endif

        if (s3sq == 0.0) then
          ! There is a simple, exact solution to the quadratic equation, namely:
          stays = h1 ! This will revert to stays_max later.
        elseif (s2*s2 <= s3sq) then
          ! There is no solution with 0 PE change - use the minimum energy input.
          stays = s1
        else
          ! The following choose the solutions that are continuous with all water
          ! staying in the upper buffer layer when there is no detrainment,
          ! namely the + root when s2>0 and the - root otherwise. They also
          ! carefully avoid differencing large numbers, using s2 = (h1-s).
          if (bh0 <= 0.0) then ; stays = h1
          elseif (s2 > 0.0) then
  !         stays = s + sqrt(s2*s2 - s3sq) ! Note that s2 = h1-s
            if (s1 >= stays_max) then ; stays = stays_max
            elseif (s1 >= 0.0) then ; stays = s1 + sqrt(s2*s2 - s3sq)
            else ; stays = (h1*(s2-s1) - s3sq) / (-s1 + sqrt(s2*s2 - s3sq))
            endif
          else
  !         stays = s - sqrt(s2*s2 - s3sq) ! Note that s2 = h1-s & stays_min >= 0
            if (s1 <= stays_min) then ; stays = stays_min
            else ; stays = (h1*(s1-s2) + s3sq) / (s1 + sqrt(s2*s2 - s3sq))
            endif
          endif
        endif

        ! Limit the amount that stays so that the motion of water is from the
        ! upper buffer layer into the lower, but no more than is in the upper
        ! layer, and the water left in the upper layer is no lighter than the
        ! detrained water.
        if (stays >= stays_max) then ; stays = stays_max
        elseif (stays < stays_min) then ; stays = stays_min
        endif
      endif

      dPE_det = g_2*((R0(i,kb1)*h_to_bl - R0_to_bl)*stays + &
                     (R0(i,kb2)-R0(i,kb1)) * (h1-stays) * &
                     (h2 - scale_slope*stays*((h1+h2)+h_to_bl)/(h1+h2)) ) - &
                Rho0xG*dPE_extrap

      if (dPE_time_ratio*h_to_bl > h_to_bl+h(i,0)) then
        dPE_ratio = (h_to_bl+h(i,0)) / h_to_bl
      else
        dPE_ratio = dPE_time_ratio
      endif

      if ((mergeable_bl) .and. (num_events*dPE_ratio*dPE_det > dPE_merge)) then
        ! It is energetically preferable to merge the two buffer layers, detrain
        ! them into interior layer (k0), move the remaining upper buffer layer
        ! water into the lower buffer layer, and detrain undiluted into the
        ! upper buffer layer.
        h1_to_k0 = (h1-stays_merge)
        stays = MAX(h_min_bl-h_to_bl,0.0)
        h1_to_h2 = stays_merge - stays

        Ihk0 = 1.0 / ((h1_to_k0 + h2) + h(i,k0))
        Ih1f = 1.0 / (h_to_bl + stays); Ih2f = 1.0 / h1_to_h2
        Ih12 = 1.0 / (h1 + h2)

        dRcv_2dz = (Rcv(i,kb1) - Rcv(i,kb2)) * Ih12
        dRcv_stays = dRcv_2dz*(h1_to_k0 + h1_to_h2)
        dRcv_det = - dRcv_2dz*(stays + h1_to_h2)
        Rcv(i,k0) = ((h1_to_k0*(Rcv(i,kb1) + dRcv_det) + &
                      h2*Rcv(i,kb2)) + h(i,k0)*Rcv(i,k0)) * Ihk0
        Rcv(i,kb2) = Rcv(i,kb1) + dRcv_2dz*(h1_to_k0-stays)
        Rcv(i,kb1) = (Rcv_to_bl + stays*(Rcv(i,kb1) + dRcv_stays)) * Ih1f

        ! Use 2nd order upwind advection of spiciness, limited by the value in
        ! the water from the mixed layer to determine the temperature and
        ! salinity of the water that stays in the buffer layers.
        I_denom = 1.0 / (dRcv_dS(i)**2 + (dT_dS_gauge*dRcv_dT(i))**2)
        dSpice_2dz = (dS_dT_gauge*dRcv_dS(i)*(T(i,kb1)-T(i,kb2)) - &
                      dT_dS_gauge*dRcv_dT(i)*(S(i,kb1)-S(i,kb2))) * Ih12
        dSpice_lim = (dS_dT_gauge*dR0_dS(i)*(T_to_bl-T(i,kb1)*h_to_bl) - &
                      dT_dS_gauge*dR0_dT(i)*(S_to_bl-S(i,kb1)*h_to_bl)) / h_to_bl
        if (dSpice_lim * dSpice_2dz <= 0.0) dSpice_2dz = 0.0

        if (stays > 0.0) then
        ! Limit the spiciness of the water that stays in the upper buffer layer.
          if (abs(dSpice_lim) < abs(dSpice_2dz*(h1_to_k0 + h1_to_h2))) &
            dSpice_2dz = dSpice_lim/(h1_to_k0 + h1_to_h2)

          dSpice_stays = dSpice_2dz*(h1_to_k0 + h1_to_h2)
          T_stays = T(i,kb1) + dT_dS_gauge * I_denom * &
              (dT_dS_gauge * dRcv_dT(i) * dRcv_stays + dRcv_dS(i) * dSpice_stays)
          S_stays = S(i,kb1) + I_denom * &
              (dRcv_dS(i) * dRcv_stays - dT_dS_gauge * dRcv_dT(i) * dSpice_stays)
          ! The values of R0 are based on changes in T and S.
          R0_stays = R0(i,kb1) + (T_stays-T(i,kb1)) * dR0_dT(i) + &
                                 (S_stays-S(i,kb1)) * dR0_dS(i)
        else
          ! Limit the spiciness of the water that moves into the lower buffer layer.
          if (abs(dSpice_lim) < abs(dSpice_2dz*h1_to_k0)) &
            dSpice_2dz = dSpice_lim/h1_to_k0
          ! These will be multiplied by 0 later.
          T_stays = 0.0 ; S_stays = 0.0 ; R0_stays = 0.0
        endif

        dSpice_det = - dSpice_2dz*(stays + h1_to_h2)
        T_det = T(i,kb1) + dT_dS_gauge * I_denom * &
            (dT_dS_gauge * dRcv_dT(i) * dRcv_det + dRcv_dS(i) * dSpice_det)
        S_det = S(i,kb1) + I_denom * &
            (dRcv_dS(i) * dRcv_det - dT_dS_gauge * dRcv_dT(i) * dSpice_det)
        ! The values of R0 are based on changes in T and S.
        R0_det = R0(i,kb1) + (T_det-T(i,kb1)) * dR0_dT(i) + &
                             (S_det-S(i,kb1)) * dR0_dS(i)

        T(i,k0) = ((h1_to_k0*T_det + h2*T(i,kb2)) + h(i,k0)*T(i,k0)) * Ihk0
        T(i,kb2) = (h1*T(i,kb1) - stays*T_stays - h1_to_k0*T_det) * Ih2f
        T(i,kb1) = (T_to_bl + stays*T_stays) * Ih1f

        S(i,k0) = ((h1_to_k0*S_det + h2*S(i,kb2)) + h(i,k0)*S(i,k0)) * Ihk0
        S(i,kb2) = (h1*S(i,kb1) - stays*S_stays - h1_to_k0*S_det) * Ih2f
        S(i,kb1) = (S_to_bl + stays*S_stays) * Ih1f

        R0(i,k0) = ((h1_to_k0*R0_det + h2*R0(i,kb2)) + h(i,k0)*R0(i,k0)) * Ihk0
        R0(i,kb2) = (h1*R0(i,kb1) - stays*R0_stays - h1_to_k0*R0_det) * Ih2f
        R0(i,kb1) = (R0_to_bl + stays*R0_stays) * Ih1f

!        ! The following is 2nd-order upwind advection without limiters.
!        dT_2dz = (T(i,kb1) - T(i,kb2)) * Ih12
!        T(i,k0) = (h1_to_k0*(T(i,kb1) - dT_2dz*(stays+h1_to_h2)) + &
!                     h2*T(i,kb2) + h(i,k0)*T(i,k0)) * Ihk0
!        T(i,kb2) = T(i,kb1) + dT_2dz*(h1_to_k0-stays)
!        T(i,kb1) = (T_to_bl + stays*(T(i,kb1) + &
!                      dT_2dz*(h1_to_k0 + h1_to_h2))) * Ih1f
!        dS_2dz = (S(i,kb1) - S(i,kb2)) * Ih12
!        S(i,k0) = (h1_to_k0*(S(i,kb1) - dS_2dz*(stays+h1_to_h2)) + &
!                     h2*S(i,kb2) + h(i,k0)*S(i,k0)) * Ihk0
!        S(i,kb2) = S(i,kb1) + dS_2dz*(h1_to_k0-stays)
!        S(i,kb1) = (S_to_bl + stays*(S(i,kb1) + &
!                      dS_2dz*(h1_to_k0 + h1_to_h2))) * Ih1f
!        dR0_2dz = (R0(i,kb1) - R0(i,kb2)) * Ih12
!        R0(i,k0) = (h1_to_k0*(R0(i,kb1) - dR0_2dz*(stays+h1_to_h2)) + &
!                    h2*R0(i,kb2) + h(i,k0)*R0(i,k0)) * Ihk0
!        R0(i,kb2) = R0(i,kb1) + dR0_2dz*(h1_to_k0-stays)
!        R0(i,kb1) = (R0_to_bl + stays*(R0(i,kb1) + &
!                     dR0_2dz*(h1_to_k0 + h1_to_h2))) * Ih1f

        d_ea(i,kb1) = (d_ea(i,kb1) + h_to_bl) + (stays - h1)
        d_ea(i,kb2) = d_ea(i,kb2) + (h1_to_h2 - h2)
        d_ea(i,k0) = d_ea(i,k0) + (h1_to_k0 + h2)

        h(i,kb1) = stays + h_to_bl
        h(i,kb2) = h1_to_h2
        h(i,k0) = h(i,k0) + (h1_to_k0 + h2)
        if (allocated(CS%diag_PE_detrain)) &
          CS%diag_PE_detrain(i,j) = CS%diag_PE_detrain(i,j) + Idt_H2*dPE_merge
        if (allocated(CS%diag_PE_detrain2)) CS%diag_PE_detrain2(i,j) = &
             CS%diag_PE_detrain2(i,j) + Idt_H2*(dPE_det + Rho0xG*dPE_extrap)
      else ! Not mergeable_bl.
        ! There is no further detrainment from the buffer layers, and the
        ! upper buffer layer water is distributed optimally between the
        ! upper and lower buffer layer.
        h1_to_h2 = h1 - stays
        Ih1f = 1.0 / (h_to_bl + stays) ; Ih2f = 1.0 / (h2 + h1_to_h2)
        Ih = 1.0 / (h1 + h2)
        dR0_2dz = (R0(i,kb1) - R0(i,kb2)) * Ih
        R0(i,kb2) = (h2*R0(i,kb2) + h1_to_h2*(R0(i,kb1) - &
                     scale_slope*dR0_2dz*stays)) * Ih2f
        R0(i,kb1) = (R0_to_bl + stays*(R0(i,kb1) + &
                        scale_slope*dR0_2dz*h1_to_h2)) * Ih1f

        ! Use 2nd order upwind advection of spiciness, limited by the value
        ! in the detrained water to determine the detrained temperature and
        ! salinity.
        dR0 = scale_slope*dR0_2dz*h1_to_h2
        dSpice_stays = (dS_dT_gauge*dR0_dS(i)*(T(i,kb1)-T(i,kb2)) - &
                        dT_dS_gauge*dR0_dT(i)*(S(i,kb1)-S(i,kb2))) * &
                        scale_slope*h1_to_h2 * Ih
        if (h_to_bl > 0.0) then
          dSpice_lim = (dS_dT_gauge*dR0_dS(i)*(T_to_bl-T(i,kb1)*h_to_bl) - &
                        dT_dS_gauge*dR0_dT(i)*(S_to_bl-S(i,kb1)*h_to_bl)) /&
                        h_to_bl
        else
          dSpice_lim = dS_dT_gauge*dR0_dS(i)*(T(i,0)-T(i,kb1)) - &
                       dT_dS_gauge*dR0_dT(i)*(S(i,0)-S(i,kb1))
        endif
        if (dSpice_stays*dSpice_lim <= 0.0) then
          dSpice_stays = 0.0
        elseif (abs(dSpice_stays) > abs(dSpice_lim)) then
          dSpice_stays = dSpice_lim
        endif
        I_denom = 1.0 / (dR0_dS(i)**2 + (dT_dS_gauge*dR0_dT(i))**2)
        T_stays = T(i,kb1) + dT_dS_gauge * I_denom * &
            (dT_dS_gauge * dR0_dT(i) * dR0 + dR0_dS(i) * dSpice_stays)
        S_stays = S(i,kb1) + I_denom * &
            (dR0_dS(i) * dR0 - dT_dS_gauge * dR0_dT(i) * dSpice_stays)
        ! The detrained values of Rcv are based on changes in T and S.
        Rcv_stays = Rcv(i,kb1) + (T_stays-T(i,kb1)) * dRcv_dT(i) + &
                                 (S_stays-S(i,kb1)) * dRcv_dS(i)

        T(i,kb2) = (h2*T(i,kb2) + h1*T(i,kb1) - T_stays*stays) * Ih2f
        T(i,kb1) = (T_to_bl + stays*T_stays) * Ih1f
        S(i,kb2) = (h2*S(i,kb2) + h1*S(i,kb1) - S_stays*stays) * Ih2f
        S(i,kb1) = (S_to_bl + stays*S_stays) * Ih1f
        Rcv(i,kb2) = (h2*Rcv(i,kb2) + h1*Rcv(i,kb1) - Rcv_stays*stays) * Ih2f
        Rcv(i,kb1) = (Rcv_to_bl + stays*Rcv_stays) * Ih1f

!        ! The following is 2nd-order upwind advection without limiters.
!        dRcv_2dz = (Rcv(i,kb1) - Rcv(i,kb2)) * Ih
!        dRcv = scale_slope*dRcv_2dz*h1_to_h2
!        Rcv(i,kb2) = (h2*Rcv(i,kb2) + h1_to_h2*(Rcv(i,kb1) - &
!                      scale_slope*dRcv_2dz*stays)) * Ih2f
!        Rcv(i,kb1) = (Rcv_to_bl + stays*(Rcv(i,kb1) + dRcv)) * Ih1f
!        dT_2dz = (T(i,kb1) - T(i,kb2)) * Ih
!        T(i,kb2) = (h2*T(i,kb2) + h1_to_h2*(T(i,kb1) - &
!                      scale_slope*dT_2dz*stays)) * Ih2f
!        T(i,kb1) = (T_to_bl + stays*(T(i,kb1) + &
!                      scale_slope*dT_2dz*h1_to_h2)) * Ih1f
!        dS_2dz = (S(i,kb1) - S(i,kb2)) * Ih
!        S(i,kb2) = (h2*S(i,kb2) + h1_to_h2*(S(i,kb1) - &
!                      scale_slope*dS_2dz*stays)) * Ih2f
!        S(i,kb1) = (S_to_bl + stays*(S(i,kb1) + &
!                      scale_slope*dS_2dz*h1_to_h2)) * Ih1f

        d_ea(i,kb1) = d_ea(i,kb1) + ((stays - h1) + h_to_bl)
        d_ea(i,kb2) = d_ea(i,kb2) + h1_to_h2

        h(i,kb1) = stays + h_to_bl
        h(i,kb2) = h(i,kb2) + h1_to_h2

        if (allocated(CS%diag_PE_detrain)) &
          CS%diag_PE_detrain(i,j) = CS%diag_PE_detrain(i,j) + Idt_H2*dPE_det
        if (allocated(CS%diag_PE_detrain2)) CS%diag_PE_detrain2(i,j) = &
          CS%diag_PE_detrain2(i,j) + Idt_H2*(dPE_det + Rho0xG*dPE_extrap)
      endif
    endif ! End of detrainment...

  enddo ! i loop

end subroutine mixedlayer_detrain_2

!> This subroutine moves any water left in the former mixed layers into the
!! single buffer layers and may also move buffer layer water into the interior
!! isopycnal layers.
subroutine mixedlayer_detrain_1(h, T, S, R0, Rcv, RcvTgt, dt, dt_diag, d_ea, d_eb, &
                                j, G, GV, US, CS, dRcv_dT, dRcv_dS, max_BL_det)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: h    !< Layer thickness [H ~> m or kg m-2].
                                                            !! Layer 0 is the new mixed layer.
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: T    !< Potential temperature [degC].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: S    !< Salinity [ppt].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: R0   !< Potential density referenced to
                                                            !! surface pressure [R ~> kg m-3].
  real, dimension(SZI_(G),SZK0_(GV)), intent(inout) :: Rcv  !< The coordinate defining potential
                                                            !! density [R ~> kg m-3].
  real, dimension(SZK_(GV)),          intent(in)    :: RcvTgt !< The target value of Rcv for each
                                                            !! layer [R ~> kg m-3].
  real,                               intent(in)    :: dt   !< Time increment [T ~> s].
  real,                               intent(in)    :: dt_diag !< The accumulated time interval for
                                                            !! diagnostics [T ~> s].
  real, dimension(SZI_(G),SZK_(GV)),  intent(inout) :: d_ea !< The upward increase across a layer in
                                                            !! the entrainment from above
                                                            !! [H ~> m or kg m-2]. Positive d_ea
                                                            !! goes with layer thickness increases.
  real, dimension(SZI_(G),SZK_(GV)),  intent(inout) :: d_eb !< The downward increase across a layer
                                                            !! in the entrainment from below [H ~> m or kg m-2].
                                                            !! Positive values go with mass gain by
                                                            !! a layer.
  integer,                            intent(in)    :: j    !< The meridional row to work on.
  type(unit_scale_type),              intent(in)    :: US   !< A dimensional unit scaling type
  type(bulkmixedlayer_CS),            pointer       :: CS   !< The control structure returned by a
                                                            !! previous call to mixedlayer_init.
  real, dimension(SZI_(G)),           intent(in)    :: dRcv_dT !< The partial derivative of
                                                            !! coordinate defining potential density
                                                            !! with potential temperature
                                                            !! [R degC-1 ~> kg m-3 degC-1].
  real, dimension(SZI_(G)),           intent(in)    :: dRcv_dS    !< The partial derivative of
                                                            !! coordinate defining potential density
                                                            !! with salinity [R ppt-1 ~> kg m-3 ppt-1].
  real, dimension(SZI_(G)),           intent(in)    :: max_BL_det !< If non-negative, the maximum
                                                            !! detrainment permitted from the buffer
                                                            !! layers [H ~> m or kg m-2].

  ! Local variables
  real :: Ih                  ! The inverse of a thickness [H-1 ~> m-1 or m2 kg-1].
  real :: h_ent               ! The thickness from a layer that is
                              ! entrained [H ~> m or kg m-2].
  real :: max_det_rem(SZI_(G)) ! Remaining permitted detrainment [H ~> m or kg m-2].
  real :: detrain(SZI_(G))    ! The thickness of fluid to detrain
                              ! from the mixed layer [H ~> m or kg m-2].
  real :: dT_dR, dS_dR, dRml, dR0_dRcv, dT_dS_wt2
  real :: I_denom             ! A work variable [ppt2 R-2 ~> ppt2 m6 kg-2].
  real :: Sdown, Tdown
  real :: dt_Time             ! The timestep divided by the detrainment timescale [nondim].
  real :: g_H2_2Rho0dt        ! Half the gravitational acceleration times the square of the
                              ! conversion from H to m divided by the mean density times the time
                              ! step [L2 Z T-3 H-2 R-1 ~> m4 s-3 kg-1 or m10 s-3 kg-3].
  real :: g_H2_2dt            ! Half the gravitational acceleration times the square of the
                              ! conversion from H to Z divided by the diagnostic time step
                              ! [L2 Z H-2 T-3 ~> m s-3 or m7 kg-2 s-3].

  logical :: splittable_BL(SZI_(G)), orthogonal_extrap
  real :: x1

  integer :: i, is, ie, k, k1, nkmb, nz
  is = G%isc ; ie = G%iec ; nz = GV%ke
  nkmb = CS%nkml+CS%nkbl
  if (CS%nkbl /= 1) call MOM_error(FATAL,"MOM_mixed_layer: "// &
                        "CS%nkbl must be 1 in mixedlayer_detrain_1.")

  dt_Time = dt / CS%BL_detrain_time
  g_H2_2Rho0dt = (GV%g_Earth * GV%H_to_Z**2) / (2.0 * GV%Rho0 * dt_diag)
  g_H2_2dt = (GV%g_Earth * GV%H_to_Z**2) / (2.0 * dt_diag)

  ! Move detrained water into the buffer layer.
  do k=1,CS%nkml
    do i=is,ie ; if (h(i,k) > 0.0) then
      Ih = 1.0 / (h(i,nkmb) + h(i,k))
      if (CS%TKE_diagnostics) &
        CS%diag_TKE_conv_s2(i,j) = CS%diag_TKE_conv_s2(i,j) + &
            g_H2_2Rho0dt * h(i,k) * h(i,nkmb) * (R0(i,nkmb) - R0(i,k))
      if (allocated(CS%diag_PE_detrain)) &
        CS%diag_PE_detrain(i,j) = CS%diag_PE_detrain(i,j) + &
            g_H2_2dt * h(i,k) * h(i,nkmb) * (R0(i,nkmb) - R0(i,k))
      if (allocated(CS%diag_PE_detrain2)) &
        CS%diag_PE_detrain2(i,j) = CS%diag_PE_detrain2(i,j) + &
            g_H2_2dt * h(i,k) * h(i,nkmb) * (R0(i,nkmb) - R0(i,k))

      R0(i,nkmb) = (R0(i,nkmb)*h(i,nkmb) + R0(i,k)*h(i,k)) * Ih
      Rcv(i,nkmb) = (Rcv(i,nkmb)*h(i,nkmb) + Rcv(i,k)*h(i,k)) * Ih
      T(i,nkmb) = (T(i,nkmb)*h(i,nkmb) + T(i,k)*h(i,k)) * Ih
      S(i,nkmb) = (S(i,nkmb)*h(i,nkmb) + S(i,k)*h(i,k)) * Ih

      d_ea(i,k) = d_ea(i,k) - h(i,k)
      d_ea(i,nkmb) = d_ea(i,nkmb) + h(i,k)
      h(i,nkmb) = h(i,nkmb) + h(i,k)
      h(i,k) = 0.0
    endif ; enddo
  enddo

  do i=is,ie
    max_det_rem(i) = 10.0 * h(i,nkmb)
    if (max_BL_det(i) >= 0.0) max_det_rem(i) = max_BL_det(i)
  enddo

!   If the mixed layer was denser than the densest interior layer,
! but is now lighter than this layer, leaving a buffer layer that
! is denser than this layer, there are problems.  This should prob-
! ably be considered a case of an inadequate choice of resolution in
! density space and should be avoided.  To make the model run sens-
! ibly in this case, it will make the mixed layer denser while making
! the buffer layer the density of the densest interior layer (pro-
! vided that the this will not make the mixed layer denser than the
! interior layer).  Otherwise, make the mixed layer the same density
! as the densest interior layer and lighten the buffer layer with
! the released buoyancy.  With multiple buffer layers, much more
! graceful options are available.
  do i=is,ie ; if (h(i,nkmb) > 0.0) then
    if ((R0(i,0)<R0(i,nz)) .and. (R0(i,nz)<R0(i,nkmb))) then
      if ((R0(i,nz)-R0(i,0))*h(i,0) > &
          (R0(i,nkmb)-R0(i,nz))*h(i,nkmb)) then
        detrain(i) = (R0(i,nkmb)-R0(i,nz))*h(i,nkmb) / (R0(i,nkmb)-R0(i,0))
      else
        detrain(i) = (R0(i,nz)-R0(i,0))*h(i,0) / (R0(i,nkmb)-R0(i,0))
      endif

      d_eb(i,CS%nkml) = d_eb(i,CS%nkml) + detrain(i)
      d_ea(i,CS%nkml) = d_ea(i,CS%nkml) - detrain(i)
      d_eb(i,nkmb) = d_eb(i,nkmb) - detrain(i)
      d_ea(i,nkmb) = d_ea(i,nkmb) + detrain(i)

      if (allocated(CS%diag_PE_detrain)) CS%diag_PE_detrain(i,j) = &
        CS%diag_PE_detrain(i,j) + g_H2_2dt * detrain(i)* &
                     (h(i,0) + h(i,nkmb)) * (R0(i,nkmb) - R0(i,0))
      x1 = R0(i,0)
      R0(i,0) = R0(i,0) - detrain(i)*(R0(i,0)-R0(i,nkmb)) / h(i,0)
      R0(i,nkmb) = R0(i,nkmb) - detrain(i)*(R0(i,nkmb)-x1) / h(i,nkmb)
      x1 = Rcv(i,0)
      Rcv(i,0) = Rcv(i,0) - detrain(i)*(Rcv(i,0)-Rcv(i,nkmb)) / h(i,0)
      Rcv(i,nkmb) = Rcv(i,nkmb) - detrain(i)*(Rcv(i,nkmb)-x1) / h(i,nkmb)
      x1 = T(i,0)
      T(i,0) = T(i,0) - detrain(i)*(T(i,0)-T(i,nkmb)) / h(i,0)
      T(i,nkmb) = T(i,nkmb) - detrain(i)*(T(i,nkmb)-x1) / h(i,nkmb)
      x1 = S(i,0)
      S(i,0) = S(i,0) - detrain(i)*(S(i,0)-S(i,nkmb)) / h(i,0)
      S(i,nkmb) = S(i,nkmb) - detrain(i)*(S(i,nkmb)-x1) / h(i,nkmb)

    endif
  endif ; enddo

  ! Move water out of the buffer layer, if convenient.
!   Split the buffer layer if possible, and replace the buffer layer
! with a small amount of fluid from the mixed layer.
! This is the exponential-in-time splitting, circa 2005.
  do i=is,ie
    if (h(i,nkmb) > 0.0) then ; splittable_BL(i) = .true.
    else ; splittable_BL(i) = .false. ; endif
  enddo

  dT_dS_wt2 = CS%dT_dS_wt**2

  do k=nz-1,nkmb+1,-1 ; do i=is,ie
    if (splittable_BL(i)) then
      if (RcvTgt(k)<=Rcv(i,nkmb)) then
! Estimate dR/drho, dTheta/dR, and dS/dR, where R is the coordinate variable
! and rho is in-situ (or surface) potential density.
! There is no "right" way to do this, so this keeps things reasonable, if
! slightly arbitrary.
        splittable_BL(i) = .false.

        k1 = k+1 ; orthogonal_extrap = .false.
        ! Here we try to find a massive layer to use for interpolating the
        ! temperature and salinity.  If none is available a pseudo-orthogonal
        ! extrapolation is used.  The 10.0 and 0.9 in the following are
        ! arbitrary but probably about right.
        if ((h(i,k+1) < 10.0*GV%Angstrom_H) .or. &
            ((RcvTgt(k+1)-Rcv(i,nkmb)) >= 0.9*(Rcv(i,k1) - Rcv(i,0)))) then
          if (k>=nz-1) then ; orthogonal_extrap = .true.
          elseif ((h(i,k+2) <= 10.0*GV%Angstrom_H) .and. &
              ((RcvTgt(k+1)-Rcv(i,nkmb)) < 0.9*(Rcv(i,k+2)-Rcv(i,0)))) then
            k1 = k+2
          else ; orthogonal_extrap = .true. ; endif
        endif

        if ((R0(i,0) >= R0(i,k1)) .or. (Rcv(i,0) >= Rcv(i,nkmb))) cycle
          ! In this case there is an inversion of in-situ density relative to
          ! the coordinate variable.  Do not detrain from the buffer layer.

        if (orthogonal_extrap) then
          ! 36 here is a typical oceanic value of (dR/dS) / (dR/dT) - it says
          ! that the relative weights of T & S changes is a plausible 6:1.
          ! Also, this was coded on Athena's 6th birthday!
          I_denom = 1.0 / (dRcv_dS(i)**2 + dT_dS_wt2*dRcv_dT(i)**2)
          dT_dR = dT_dS_wt2*dRcv_dT(i) * I_denom
          dS_dR = dRcv_dS(i) * I_denom
        else
          dT_dR = (T(i,0) - T(i,k1)) / (Rcv(i,0) - Rcv(i,k1))
          dS_dR = (S(i,0) - S(i,k1)) / (Rcv(i,0) - Rcv(i,k1))
        endif
        dRml = dt_Time * (R0(i,nkmb) - R0(i,0)) * &
               (Rcv(i,0) - Rcv(i,k1)) / (R0(i,0) - R0(i,k1))
        ! Once again, there is an apparent density inversion in Rcv.
        if (dRml < 0.0) cycle
        dR0_dRcv = (R0(i,0) - R0(i,k1)) / (Rcv(i,0) - Rcv(i,k1))

        if ((Rcv(i,nkmb) - dRml < RcvTgt(k)) .and. (max_det_rem(i) > h(i,nkmb))) then
          ! In this case, the buffer layer is split into two isopycnal layers.
          detrain(i) = h(i,nkmb)*(Rcv(i,nkmb) - RcvTgt(k)) / &
                                  (RcvTgt(k+1) - RcvTgt(k))

          if (allocated(CS%diag_PE_detrain)) CS%diag_PE_detrain(i,j) = &
            CS%diag_PE_detrain(i,j) - g_H2_2dt * detrain(i) * &
                 (h(i,nkmb)-detrain(i)) * (RcvTgt(k+1) - RcvTgt(k)) * dR0_dRcv

          Tdown = detrain(i) * (T(i,nkmb) + dT_dR*(RcvTgt(k+1)-Rcv(i,nkmb)))
          T(i,k) = (h(i,k) * T(i,k) + &
                        (h(i,nkmb) * T(i,nkmb) - Tdown)) / &
                       (h(i,k) + (h(i,nkmb) - detrain(i)))
          T(i,k+1) = (h(i,k+1) * T(i,k+1) + Tdown)/ &
                          (h(i,k+1) + detrain(i))
          T(i,nkmb) = T(i,0)
          Sdown = detrain(i) * (S(i,nkmb) + dS_dR*(RcvTgt(k+1)-Rcv(i,nkmb)))
          S(i,k) = (h(i,k) * S(i,k) + &
                      (h(i,nkmb) * S(i,nkmb) - Sdown)) / &
                      (h(i,k) + (h(i,nkmb) - detrain(i)))
          S(i,k+1) = (h(i,k+1) * S(i,k+1) + Sdown)/ &
                         (h(i,k+1) + detrain(i))
          S(i,nkmb) = S(i,0)
          Rcv(i,nkmb) = Rcv(i,0)

          d_ea(i,k+1) = d_ea(i,k+1) + detrain(i)
          d_ea(i,k) = d_ea(i,k) + (h(i,nkmb) - detrain(i))
          d_ea(i,nkmb) = d_ea(i,nkmb) - h(i,nkmb)

          h(i,k+1) = h(i,k+1) + detrain(i)
          h(i,k) = h(i,k) + h(i,nkmb) - detrain(i)
          h(i,nkmb) = 0.0
        else
          ! Here only part of the buffer layer is moved into the interior.
          detrain(i) = h(i,nkmb) * dRml / (RcvTgt(k+1) - Rcv(i,nkmb) + dRml)
          if (detrain(i) > max_det_rem(i)) detrain(i) = max_det_rem(i)
          Ih = 1.0 / (h(i,k+1) + detrain(i))

          Tdown = (T(i,nkmb) + dT_dR*(RcvTgt(k+1)-Rcv(i,nkmb)))
          T(i,nkmb) = T(i,nkmb) - dT_dR * dRml
          T(i,k+1) = (h(i,k+1) * T(i,k+1) + detrain(i) * Tdown) * Ih
          Sdown = (S(i,nkmb) + dS_dR*(RcvTgt(k+1)-Rcv(i,nkmb)))
!  The following two expressions updating S(nkmb) are mathematically identical.
!            S(i,nkmb) = (h(i,nkmb) * S(i,nkmb) - detrain(i) * Sdown) / &
!                           (h(i,nkmb) - detrain(i))
          S(i,nkmb) = S(i,nkmb) - dS_dR * dRml
          S(i,k+1) = (h(i,k+1) * S(i,k+1) + detrain(i) * Sdown) * Ih

          d_ea(i,k+1) = d_ea(i,k+1) + detrain(i)
          d_ea(i,nkmb) = d_ea(i,nkmb) - detrain(i)

          h(i,k+1) = h(i,k+1) + detrain(i)
          h(i,nkmb) = h(i,nkmb) - detrain(i)

          if (allocated(CS%diag_PE_detrain)) CS%diag_PE_detrain(i,j) = &
            CS%diag_PE_detrain(i,j) - g_H2_2dt * detrain(i) * dR0_dRcv * &
                 (h(i,nkmb)-detrain(i)) * (RcvTgt(k+1) - Rcv(i,nkmb) + dRml)
        endif
      endif ! RcvTgt(k)<=Rcv(i,nkmb)
    endif ! splittable_BL
  enddo ; enddo ! i & k loops

!   The numerical behavior of the buffer layer is dramatically improved
! if it is always at least a small fraction (say 10%) of the thickness
! of the mixed layer.  As the physical distinction between the mixed
! and buffer layers is vague anyway, this seems hard to argue against.
  do i=is,ie
    if (h(i,nkmb) < 0.1*h(i,0)) then
      h_ent =  0.1*h(i,0) - h(i,nkmb)
      Ih = 10.0/h(i,0)
      T(i,nkmb) = (h(i,nkmb)*T(i,nkmb) + h_ent*T(i,0)) * Ih
      S(i,nkmb) = (h(i,nkmb)*S(i,nkmb) + h_ent*S(i,0)) * Ih

      d_ea(i,1) = d_ea(i,1) - h_ent
      d_ea(i,nkmb) = d_ea(i,nkmb) + h_ent

      h(i,0) = h(i,0) - h_ent
      h(i,nkmb) = h(i,nkmb) + h_ent
    endif
  enddo

end subroutine mixedlayer_detrain_1

!> This subroutine initializes the MOM bulk mixed layer module.
subroutine bulkmixedlayer_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< The model's clock with the current time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(bulkmixedlayer_CS), pointer       :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module.
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
  character(len=40)  :: mdl = "MOM_mixed_layer"  ! This module's name.
  real :: BL_detrain_time_dflt ! The default value for BUFFER_LAY_DETRAIN_TIME [s]
  real :: omega_frac_dflt, ustar_min_dflt, Hmix_min_m
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

  if (GV%nkml < 1) return

! Set default, read and log parameters
  call log_version(param_file, mdl, version, "")

  CS%nkml = GV%nkml
  call log_param(param_file, mdl, "NKML", CS%nkml, &
                 "The number of sublayers within the mixed layer if "//&
                 "BULKMIXEDLAYER is true.", units="nondim", default=2)
  CS%nkbl = GV%nk_rho_varies - GV%nkml
  call log_param(param_file, mdl, "NKBL", CS%nkbl, &
                 "The number of variable density buffer layers if "//&
                 "BULKMIXEDLAYER is true.", units="nondim", default=2)
  call get_param(param_file, mdl, "MSTAR", CS%mstar, &
                 "The ratio of the friction velocity cubed to the TKE "//&
                 "input to the mixed layer.", "units=nondim", default=1.2)
  call get_param(param_file, mdl, "NSTAR", CS%nstar, &
                 "The portion of the buoyant potential energy imparted by "//&
                 "surface fluxes that is available to drive entrainment "//&
                 "at the base of mixed layer when that energy is positive.", &
                 units="nondim", default=0.15)
  call get_param(param_file, mdl, "BULK_RI_ML", CS%bulk_Ri_ML, &
                 "The efficiency with which mean kinetic energy released "//&
                 "by mechanically forced entrainment of the mixed layer "//&
                 "is converted to turbulent kinetic energy.", units="nondim",&
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "ABSORB_ALL_SW", CS%absorb_all_sw, &
                 "If true,  all shortwave radiation is absorbed by the "//&
                 "ocean, instead of passing through to the bottom mud.", &
                 default=.false.)
  call get_param(param_file, mdl, "TKE_DECAY", CS%TKE_decay, &
                 "TKE_DECAY relates the vertical rate of decay of the "//&
                 "TKE available for mechanical entrainment to the natural "//&
                 "Ekman depth.", units="nondim", default=2.5)
  call get_param(param_file, mdl, "NSTAR2", CS%nstar2, &
                 "The portion of any potential energy released by "//&
                 "convective adjustment that is available to drive "//&
                 "entrainment at the base of mixed layer. By default "//&
                 "NSTAR2=NSTAR.", units="nondim", default=CS%nstar)
  call get_param(param_file, mdl, "BULK_RI_CONVECTIVE", CS%bulk_Ri_convective, &
                 "The efficiency with which convectively released mean "//&
                 "kinetic energy is converted to turbulent kinetic "//&
                 "energy.  By default BULK_RI_CONVECTIVE=BULK_RI_ML.", &
                 units="nondim", default=CS%bulk_Ri_ML)
  call get_param(param_file, mdl, "HMIX_MIN", CS%Hmix_min, &
                 "The minimum mixed layer depth if the mixed layer depth "//&
                 "is determined dynamically.", units="m", default=0.0, scale=GV%m_to_H, &
                 unscaled=Hmix_min_m)

  call get_param(param_file, mdl, "LIMIT_BUFFER_DETRAIN", CS%limit_det, &
                 "If true, limit the detrainment from the buffer layers "//&
                 "to not be too different from the neighbors.", default=.false.)
  call get_param(param_file, mdl, "ALLOWED_DETRAIN_TEMP_CHG", CS%Allowed_T_chg, &
                 "The amount by which temperature is allowed to exceed "//&
                 "previous values during detrainment.", units="K", default=0.5)
  call get_param(param_file, mdl, "ALLOWED_DETRAIN_SALT_CHG", CS%Allowed_S_chg, &
                 "The amount by which salinity is allowed to exceed "//&
                 "previous values during detrainment.", units="PSU", default=0.1)
  call get_param(param_file, mdl, "ML_DT_DS_WEIGHT", CS%dT_dS_wt, &
                 "When forced to extrapolate T & S to match the layer "//&
                 "densities, this factor (in deg C / PSU) is combined "//&
                 "with the derivatives of density with T & S to determine "//&
                 "what direction is orthogonal to density contours. It "//&
                 "should be a typical value of (dR/dS) / (dR/dT) in "//&
                 "oceanic profiles.", units="degC PSU-1", default=6.0)
  call get_param(param_file, mdl, "BUFFER_LAYER_EXTRAP_LIMIT", CS%BL_extrap_lim, &
                 "A limit on the density range over which extrapolation "//&
                 "can occur when detraining from the buffer layers, "//&
                 "relative to the density range within the mixed and "//&
                 "buffer layers, when the detrainment is going into the "//&
                 "lightest interior layer, nondimensional, or a negative "//&
                 "value not to apply this limit.", units="nondim", default = -1.0)
  call get_param(param_file, mdl, "BUFFER_LAYER_HMIN_THICK", CS%Hbuffer_min, &
                 "The minimum buffer layer thickness when the mixed layer is very thick.", &
                 units="m", default=5.0, scale=GV%m_to_H)
  call get_param(param_file, mdl, "BUFFER_LAYER_HMIN_REL", CS%Hbuffer_rel_min, &
                 "The minimum buffer layer thickness relative to the combined mixed "//&
                 "land buffer ayer thicknesses when they are thin.", &
                 units="nondim", default=0.1/CS%nkbl)
  BL_detrain_time_dflt = 4.0*3600.0 ; if (CS%nkbl==1) BL_detrain_time_dflt = 86400.0*30.0
  call get_param(param_file, mdl, "BUFFER_LAY_DETRAIN_TIME", CS%BL_detrain_time, &
                 "A timescale that characterizes buffer layer detrainment events.", &
                 units="s", default=BL_detrain_time_dflt, scale=US%s_to_T)
  call get_param(param_file, mdl, "BUFFER_SPLIT_RHO_TOL", CS%BL_split_rho_tol, &
                 "The fractional tolerance for matching layer target densities when splitting "//&
                 "layers to deal with massive interior layers that are lighter than one of the "//&
                 "mixed or buffer layers.", units="nondim", default=0.1)

  call get_param(param_file, mdl, "DEPTH_LIMIT_FLUXES", CS%H_limit_fluxes, &
                 "The surface fluxes are scaled away when the total ocean "//&
                 "depth is less than DEPTH_LIMIT_FLUXES.", &
                 units="m", default=0.1*Hmix_min_m, scale=GV%m_to_H)
  call get_param(param_file, mdl, "OMEGA", CS%omega, &
                 "The rotation rate of the earth.", &
                 default=7.2921e-5, units="s-1", scale=US%T_to_s)
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
  call get_param(param_file, mdl, "ML_RESORT", CS%ML_resort, &
                 "If true, resort the topmost layers by potential density "//&
                 "before the mixed layer calculations.", default=.false.)
  if (CS%ML_resort) &
    call get_param(param_file, mdl, "ML_PRESORT_NK_CONV_ADJ", CS%ML_presort_nz_conv_adj, &
                 "Convectively mix the first ML_PRESORT_NK_CONV_ADJ "//&
                 "layers before sorting when ML_RESORT is true.", &
                 units="nondim", default=0, fail_if_missing=.true.) ! Fail added by AJA.
  ! This gives a minimum decay scale that is typically much less than Angstrom.
  ustar_min_dflt = 2e-4*US%s_to_T*CS%omega*(GV%Angstrom_m + GV%H_to_m*GV%H_subroundoff)
  call get_param(param_file, mdl, "BML_USTAR_MIN", CS%ustar_min, &
                 "The minimum value of ustar that should be used by the "//&
                 "bulk mixed layer model in setting vertical TKE decay "//&
                 "scales. This must be greater than 0.", units="m s-1", &
                 default=ustar_min_dflt, scale=US%m_to_Z*US%T_to_s)
  if (CS%ustar_min<=0.0) call MOM_error(FATAL, "BML_USTAR_MIN must be positive.")

  call get_param(param_file, mdl, "RESOLVE_EKMAN", CS%Resolve_Ekman, &
                 "If true, the NKML>1 layers in the mixed layer are "//&
                 "chosen to optimally represent the impact of the Ekman "//&
                 "transport on the mixed layer TKE budget.  Otherwise, "//&
                 "the sublayers are distributed uniformly through the "//&
                 "mixed layer.", default=.false.)
  call get_param(param_file, mdl, "CORRECT_ABSORPTION_DEPTH", CS%correct_absorption, &
                 "If true, the average depth at which penetrating shortwave "//&
                 "radiation is absorbed is adjusted to match the average "//&
                 "heating depth of an exponential profile by moving some "//&
                 "of the heating upward in the water column.", default=.false.)
  call get_param(param_file, mdl, "DO_RIVERMIX", CS%do_rivermix, &
                 "If true, apply additional mixing wherever there is "//&
                 "runoff, so that it is mixed down to RIVERMIX_DEPTH, "//&
                 "if the ocean is that deep.", default=.false.)
  if (CS%do_rivermix) &
    call get_param(param_file, mdl, "RIVERMIX_DEPTH", CS%rivermix_depth, &
                 "The depth to which rivers are mixed if DO_RIVERMIX is "//&
                 "defined.", units="m", default=0.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "USE_RIVER_HEAT_CONTENT", CS%use_river_heat_content, &
                 "If true, use the fluxes%runoff_Hflx field to set the "//&
                 "heat carried by runoff, instead of using SST*CP*liq_runoff.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_CALVING_HEAT_CONTENT", CS%use_calving_heat_content, &
                 "If true, use the fluxes%calving_Hflx field to set the "//&
                 "heat carried by runoff, instead of using SST*CP*froz_runoff.", &
                 default=.false.)
  call get_param(param_file, mdl, "BULKML_CONV_MOMENTUM_BUG", CS%convect_mom_bug, &
                 "If true, use code with a bug that causes a loss of momentum conservation "//&
                 "during mixedlayer convection.", default=.true.)

  call get_param(param_file, mdl, "ALLOW_CLOCKS_IN_OMP_LOOPS", &
                 CS%allow_clocks_in_omp_loops, &
                 "If true, clocks can be called from inside loops that can "//&
                 "be threaded. To run with multiple threads, set to False.", &
                 default=.true.)

  CS%id_ML_depth = register_diag_field('ocean_model', 'h_ML', diag%axesT1, &
      Time, 'Surface mixed layer depth', 'm')
  CS%id_TKE_wind = register_diag_field('ocean_model', 'TKE_wind', diag%axesT1, &
      Time, 'Wind-stirring source of mixed layer TKE', &
      'm3 s-3', conversion=US%Z_to_m*(US%L_to_m**2)*(US%s_to_T**3))
  CS%id_TKE_RiBulk = register_diag_field('ocean_model', 'TKE_RiBulk', diag%axesT1, &
      Time, 'Mean kinetic energy source of mixed layer TKE', &
      'm3 s-3', conversion=US%Z_to_m*(US%L_to_m**2)*(US%s_to_T**3))
  CS%id_TKE_conv = register_diag_field('ocean_model', 'TKE_conv', diag%axesT1, &
      Time, 'Convective source of mixed layer TKE', &
      'm3 s-3', conversion=US%Z_to_m*(US%L_to_m**2)*(US%s_to_T**3))
  CS%id_TKE_pen_SW = register_diag_field('ocean_model', 'TKE_pen_SW', diag%axesT1, &
      Time, 'TKE consumed by mixing penetrative shortwave radation through the mixed layer', &
      'm3 s-3', conversion=US%Z_to_m*(US%L_to_m**2)*(US%s_to_T**3))
  CS%id_TKE_mixing = register_diag_field('ocean_model', 'TKE_mixing', diag%axesT1, &
      Time, 'TKE consumed by mixing that deepens the mixed layer', &
      'm3 s-3', conversion=US%Z_to_m*(US%L_to_m**2)*(US%s_to_T**3))
  CS%id_TKE_mech_decay = register_diag_field('ocean_model', 'TKE_mech_decay', diag%axesT1, &
      Time, 'Mechanical energy decay sink of mixed layer TKE', &
      'm3 s-3', conversion=US%Z_to_m*(US%L_to_m**2)*(US%s_to_T**3))
  CS%id_TKE_conv_decay = register_diag_field('ocean_model', 'TKE_conv_decay', diag%axesT1, &
      Time, 'Convective energy decay sink of mixed layer TKE', &
      'm3 s-3', conversion=US%Z_to_m*(US%L_to_m**2)*(US%s_to_T**3))
  CS%id_TKE_conv_s2 = register_diag_field('ocean_model', 'TKE_conv_s2', diag%axesT1, &
      Time, 'Spurious source of mixed layer TKE from sigma2', &
      'm3 s-3', conversion=US%Z_to_m*(US%L_to_m**2)*(US%s_to_T**3))
  CS%id_PE_detrain = register_diag_field('ocean_model', 'PE_detrain', diag%axesT1, &
      Time, 'Spurious source of potential energy from mixed layer detrainment', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  CS%id_PE_detrain2 = register_diag_field('ocean_model', 'PE_detrain2', diag%axesT1, &
      Time, 'Spurious source of potential energy from mixed layer only detrainment', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  CS%id_h_mismatch = register_diag_field('ocean_model', 'h_miss_ML', diag%axesT1, &
      Time, 'Summed absolute mismatch in entrainment terms', 'm', conversion=US%Z_to_m)
  CS%id_Hsfc_used = register_diag_field('ocean_model', 'Hs_used', diag%axesT1, &
      Time, 'Surface region thickness that is used', 'm', conversion=US%Z_to_m)
  CS%id_Hsfc_max = register_diag_field('ocean_model', 'Hs_max', diag%axesT1, &
      Time, 'Maximum surface region thickness', 'm', conversion=US%Z_to_m)
  CS%id_Hsfc_min = register_diag_field('ocean_model', 'Hs_min', diag%axesT1, &
      Time, 'Minimum surface region thickness', 'm', conversion=US%Z_to_m)
 !CS%lim_det_dH_sfc = 0.5 ; CS%lim_det_dH_bathy = 0.2 ! Technically these should not get used if limit_det is false?
  if (CS%limit_det .or. (CS%id_Hsfc_min > 0)) then
    call get_param(param_file, mdl, "LIMIT_BUFFER_DET_DH_SFC", CS%lim_det_dH_sfc, &
                 "The fractional limit in the change between grid points "//&
                 "of the surface region (mixed & buffer layer) thickness.", &
                 units="nondim", default=0.5)
    call get_param(param_file, mdl, "LIMIT_BUFFER_DET_DH_BATHY", CS%lim_det_dH_bathy, &
                 "The fraction of the total depth by which the thickness "//&
                 "of the surface region (mixed & buffer layer) is allowed "//&
                 "to change between grid points.", units="nondim", default=0.2)
  endif

  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", use_temperature, &
                 "If true, temperature and salinity are used as state "//&
                 "variables.", default=.true.)
  CS%nsw = 0
  if (use_temperature) then
    call get_param(param_file, mdl, "PEN_SW_NBANDS", CS%nsw, default=1)
  endif


  if (max(CS%id_TKE_wind, CS%id_TKE_RiBulk, CS%id_TKE_conv, CS%id_TKE_mixing, &
          CS%id_TKE_pen_SW, CS%id_TKE_mech_decay, CS%id_TKE_conv_decay) > 0) then
    call safe_alloc_alloc(CS%diag_TKE_wind, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%diag_TKE_RiBulk, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%diag_TKE_conv, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%diag_TKE_pen_SW, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%diag_TKE_mixing, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%diag_TKE_mech_decay, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%diag_TKE_conv_decay, isd, ied, jsd, jed)
    call safe_alloc_alloc(CS%diag_TKE_conv_s2, isd, ied, jsd, jed)

    CS%TKE_diagnostics = .true.
  endif
  if (CS%id_PE_detrain > 0) call safe_alloc_alloc(CS%diag_PE_detrain, isd, ied, jsd, jed)
  if (CS%id_PE_detrain2 > 0) call safe_alloc_alloc(CS%diag_PE_detrain2, isd, ied, jsd, jed)
  if (CS%id_ML_depth > 0) call safe_alloc_alloc(CS%ML_depth, isd, ied, jsd, jed)

  if (CS%allow_clocks_in_omp_loops) then
    id_clock_detrain = cpu_clock_id('(Ocean mixed layer detrain)', grain=CLOCK_ROUTINE)
    id_clock_mech = cpu_clock_id('(Ocean mixed layer mechanical entrainment)', grain=CLOCK_ROUTINE)
    id_clock_conv = cpu_clock_id('(Ocean mixed layer convection)', grain=CLOCK_ROUTINE)
    if (CS%ML_resort) then
      id_clock_resort = cpu_clock_id('(Ocean mixed layer resorting)', grain=CLOCK_ROUTINE)
    else
      id_clock_adjustment = cpu_clock_id('(Ocean mixed layer convective adjustment)', grain=CLOCK_ROUTINE)
    endif
    id_clock_EOS = cpu_clock_id('(Ocean mixed layer EOS)', grain=CLOCK_ROUTINE)
  endif

  if (CS%limit_det .or. (CS%id_Hsfc_min > 0)) &
      id_clock_pass = cpu_clock_id('(Ocean mixed layer halo updates)', grain=CLOCK_ROUTINE)


!  if (CS%limit_det) then
!  endif

end subroutine bulkmixedlayer_init

!> This subroutine returns an approximation to the integral
!!   R = exp(-L*(H+E)) integral(LH to L(H+E)) L/(1-(1+x)exp(-x)) dx.
!! The approximation to the integrand is good to within -2% at x~.3
!! and +25% at x~3.5, but the exponential deemphasizes the importance of
!! large x.  When L=0, EF4 returns E/((Ht+E)*Ht).
function EF4(Ht, En, I_L, dR_de)
  real,           intent(in)    :: Ht  !< Total thickness [H ~> m or kg m-2].
  real,           intent(in)    :: En  !< Entrainment [H ~> m or kg m-2].
  real,           intent(in)    :: I_L !< The e-folding scale [H-1 ~> m-1 or m2 kg-1]
  real, optional, intent(inout) :: dR_de !< The partial derivative of the result R with E [H-2 ~> m-2 or m4 kg-2].
  real :: EF4 !< The integral [H-1 ~> m-1 or m2 kg-1].

  ! Local variables
  real :: exp_LHpE ! A nondimensional exponential decay.
  real :: I_HpE    ! An inverse thickness plus entrainment [H-1 ~> m-1 or m2 kg-1].
  real :: Res      ! The result of the integral above [H-1 ~> m-1 or m2 kg-1].

  exp_LHpE = exp(-I_L*(En+Ht))
  I_HpE = 1.0/(Ht+En)
  Res = exp_LHpE * (En*I_HpE/Ht - 0.5*I_L*log(Ht*I_HpE) + 0.5*I_L*I_L*En)
  if (PRESENT(dR_de)) &
    dR_de = -I_L*Res + exp_LHpE*(I_HpE*I_HpE + 0.5*I_L*I_HpE + 0.5*I_L*I_L)
  EF4 = Res

end function EF4

!> \namespace mom_bulk_mixed_layer
!!
!! By Robert Hallberg, 1997 - 2005.
!!
!!   This file contains the subroutine (bulkmixedlayer) that
!! implements a Kraus-Turner-like bulk mixed layer, based on the work
!! of various people, as described in the review paper by Niiler and
!! Kraus (1979), with particular attention to the form proposed by
!! Oberhuber (JPO, 1993, 808-829), with an extension to a refied bulk
!! mixed layer as described in Hallberg (Aha Huliko'a, 2003).  The
!! physical processes portrayed in this subroutine include convective
!! adjustment and mixed layer entrainment and detrainment.
!! Penetrating shortwave radiation and an exponential decay of TKE
!! fluxes are also supported by this subroutine.  Several constants
!! can alternately be set to give a traditional Kraus-Turner mixed
!! layer scheme, although that is not the preferred option.  The
!! physical processes and arguments are described in detail below.

end module MOM_bulk_mixed_layer
