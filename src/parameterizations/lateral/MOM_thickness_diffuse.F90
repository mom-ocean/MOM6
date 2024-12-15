!> Isopycnal height diffusion (or Gent McWilliams diffusion)
module MOM_thickness_diffuse

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,             only : hchksum, uvchksum
use MOM_diag_mediator,         only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator,         only : register_diag_field, safe_alloc_ptr, time_type
use MOM_diag_mediator,         only : diag_update_remap_grids
use MOM_domains,               only : pass_var, CORNER, pass_vector
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_EOS,                   only : calculate_density, calculate_density_derivs, EOS_domain
use MOM_EOS,                   only : calculate_density_second_derivs
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_io,                    only : MOM_read_data, slasher
use MOM_interface_heights,     only : find_eta, thickness_to_dz
use MOM_isopycnal_slopes,      only : vert_fill_TS
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE_types,            only : MEKE_type
use MOM_stochastics,           only : stochastic_CS
use MOM_unit_scaling,          only : unit_scale_type
use MOM_variables,             only : thermo_var_ptrs, cont_diag_ptrs
use MOM_verticalGrid,          only : verticalGrid_type
implicit none ; private

#include <MOM_memory.h>

public thickness_diffuse, thickness_diffuse_init, thickness_diffuse_end
public thickness_diffuse_get_KH

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for thickness_diffuse
type, public :: thickness_diffuse_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real    :: Khth                !< Background isopycnal depth diffusivity [L2 T-1 ~> m2 s-1]
  real    :: Khth_Slope_Cff      !< Slope dependence coefficient of Khth [nondim]
  real    :: max_Khth_CFL        !< Maximum value of the diffusive CFL for isopycnal height diffusion [nondim]
  real    :: Khth_Min            !< Minimum value of Khth [L2 T-1 ~> m2 s-1]
  real    :: Khth_Max            !< Maximum value of Khth [L2 T-1 ~> m2 s-1], or 0 for no max
  real    :: Kh_eta_bg           !< Background isopycnal height diffusivity [L2 T-1 ~> m2 s-1]
  real    :: Kh_eta_vel          !< Velocity scale that is multiplied by the grid spacing to give
                                 !! the isopycnal height diffusivity [L T-1 ~> m s-1]
  real    :: slope_max           !< Slopes steeper than slope_max are limited in some way [Z L-1 ~> nondim]
  real    :: kappa_smooth        !< Vertical diffusivity used to interpolate more sensible values
                                 !! of T & S into thin layers [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  logical :: thickness_diffuse   !< If true, interfaces heights are diffused.
  logical :: full_depth_khth_min !< If true, KHTH_MIN is enforced throughout the whole water column.
                                 !! Otherwise, KHTH_MIN is only enforced at the surface. This parameter
                                 !! is only available when KHTH_USE_EBT_STRUCT=True and KHTH_MIN>0.
  logical :: use_FGNV_streamfn   !< If true, use the streamfunction formulation of
                                 !! Ferrari et al., 2010, which effectively emphasizes
                                 !! graver vertical modes by smoothing in the vertical.
  real    :: FGNV_scale          !< A coefficient scaling the vertical smoothing term in the
                                 !! Ferrari et al., 2010, streamfunction formulation [nondim].
  real    :: FGNV_c_min          !< A minimum wave speed used in the Ferrari et al., 2010,
                                 !! streamfunction formulation [L T-1 ~> m s-1].
  real    :: N2_floor            !< A floor for squared buoyancy frequency in the Ferrari et al., 2010,
                                 !! streamfunction formulation divided by aspect ratio rescaling factors
                                 !! [L2 Z-2 T-2 ~> s-2].
  logical :: detangle_interfaces !< If true, add 3-d structured interface height
                                 !! diffusivities to horizontally smooth jagged layers.
  real    :: detangle_time       !< If detangle_interfaces is true, this is the
                                 !! timescale over which maximally jagged grid-scale
                                 !! thickness variations are suppressed [T ~> s].  This must be
                                 !! longer than DT, or 0 (the default) to use DT.
  integer :: nkml                !< number of layers within mixed layer
  logical :: debug               !< write verbose checksums for debugging purposes
  logical :: use_GME_thickness_diffuse !< If true, passes GM coefficients to MOM_hor_visc for use
                                 !! with GME closure.
  logical :: MEKE_GEOMETRIC      !< If true, uses the GM coefficient formulation from the GEOMETRIC
                                 !! framework (Marshall et al., 2012)
  real    :: MEKE_GEOMETRIC_alpha!< The nondimensional coefficient governing the efficiency of
                                 !! the GEOMETRIC isopycnal height diffusion [nondim]
  real    :: MEKE_GEOMETRIC_epsilon !< Minimum Eady growth rate for the GEOMETRIC thickness
                                 !! diffusivity [T-1 ~> s-1].
  integer :: MEKE_GEOM_answer_date  !< The vintage of the expressions in the MEKE_GEOMETRIC
                                 !! calculation.  Values below 20190101 recover the answers from the
                                 !! original implementation, while higher values use expressions that
                                 !! satisfy rotational symmetry.
  logical :: Use_KH_in_MEKE      !< If true, uses the isopycnal height diffusivity calculated here to diffuse MEKE.
  real    :: MEKE_min_depth_diff !< The minimum total depth over which to average the diffusivity
                                 !! used for MEKE [H ~> m or kg m-2].  When the total depth is less
                                 !! than this, the diffusivity is scaled away.
  logical :: GM_src_alt          !< If true, use the GM energy conversion form S^2*N^2*kappa rather
                                 !! than the streamfunction for the GM source term for MEKE.
  integer :: MEKE_src_answer_date  !< The vintage of the expressions in the GM energy conversion
                                 !! calculation when MEKE_GM_SRC_ALT is true.  Values below 20240601
                                 !! recover the answers from the original implementation, while higher
                                 !! values use expressions that satisfy rotational symmetry.
  logical :: MEKE_src_slope_bug  !< If true, use a bug that limits the positive values, but not the
                                 !! negative values, of the slopes used when MEKE_GM_SRC_ALT is true.
                                 !! When this is true, it breaks rotational symmetry.
  logical :: use_GM_work_bug     !< If true, use the incorrect sign for the
                                 !! top-level work tendency on the top layer.
  real :: Stanley_det_coeff      !< The coefficient correlating SGS temperature variance with the mean
                                 !! temperature gradient in the deterministic part of the Stanley parameterization.
                                 !! Negative values disable the scheme. [nondim]
  logical :: read_khth           !< If true, read a file containing the spatially varying horizontal
                                 !! isopycnal height diffusivity
  logical :: use_stanley_gm      !< If true, also use the Stanley parameterization in MOM_thickness_diffuse

  type(diag_ctrl), pointer :: diag => NULL() !< structure used to regulate timing of diagnostics
  real, allocatable :: GMwork(:,:)        !< Work by isopycnal height diffusion [R Z L2 T-3 ~> W m-2]
  real, allocatable :: diagSlopeX(:,:,:)  !< Diagnostic: zonal neutral slope [Z L-1 ~> nondim]
  real, allocatable :: diagSlopeY(:,:,:)  !< Diagnostic: zonal neutral slope [Z L-1 ~> nondim]

  real, allocatable :: Kh_eta_u(:,:)    !< Isopycnal height diffusivities at u points [L2 T-1 ~> m2 s-1]
  real, allocatable :: Kh_eta_v(:,:)    !< Isopycnal height diffusivities in v points [L2 T-1 ~> m2 s-1]

  real, allocatable :: KH_u_GME(:,:,:)  !< Isopycnal height diffusivities in u-columns [L2 T-1 ~> m2 s-1]
  real, allocatable :: KH_v_GME(:,:,:)  !< Isopycnal height diffusivities in v-columns [L2 T-1 ~> m2 s-1]
  real, allocatable :: khth2d(:,:)      !< 2D isopycnal height diffusivity at h-points [L2 T-1 ~> m2 s-1]

  !>@{
  !! Diagnostic identifier
  integer :: id_uhGM    = -1, id_vhGM    = -1, id_GMwork = -1
  integer :: id_KH_u    = -1, id_KH_v    = -1, id_KH_t   = -1
  integer :: id_KH_u1   = -1, id_KH_v1   = -1, id_KH_t1  = -1
  integer :: id_slope_x = -1, id_slope_y = -1
  integer :: id_sfn_unlim_x = -1, id_sfn_unlim_y = -1, id_sfn_x = -1, id_sfn_y = -1
  !>@}
end type thickness_diffuse_CS

contains

!> Calculates isopycnal height diffusion coefficients and applies isopycnal height diffusion
!! by modifying to the layer thicknesses, h. Diffusivities are limited to ensure stability.
!! Also returns along-layer mass fluxes used in the continuity equation.
subroutine thickness_diffuse(h, uhtr, vhtr, tv, dt, G, GV, US, MEKE, VarMix, CDp, CS, STOCH)
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                      !! [L2 H ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                      !! [L2 H ~> m3 or kg]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamics structure
  real,                                       intent(in)    :: dt     !< Time increment [T ~> s]
  type(MEKE_type),                            intent(inout) :: MEKE   !< MEKE fields
  type(VarMix_CS), target,                    intent(in)    :: VarMix !< Variable mixing coefficients
  type(cont_diag_ptrs),                       intent(inout) :: CDp    !< Diagnostics for the continuity equation
  type(thickness_diffuse_CS),                 intent(inout) :: CS     !< Control structure for thickness_diffuse
  type(stochastic_CS),                        intent(inout) :: STOCH !< Stochastic control structure
  ! Local variables
  real :: e(SZI_(G),SZJ_(G),SZK_(GV)+1) ! heights of interfaces, relative to mean
                                         ! sea level [Z ~> m], positive up.
  real :: uhD(SZIB_(G),SZJ_(G),SZK_(GV)) ! Diffusive u*h fluxes [L2 H T-1 ~> m3 s-1 or kg s-1]
  real :: vhD(SZI_(G),SZJB_(G),SZK_(GV)) ! Diffusive v*h fluxes [L2 H T-1 ~> m3 s-1 or kg s-1]

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: &
    KH_u, &       ! Isopycnal height diffusivities in u-columns [L2 T-1 ~> m2 s-1]
    int_slope_u   ! A nondimensional ratio from 0 to 1 that gives the relative
                  ! weighting of the interface slopes to that calculated also
                  ! using density gradients at u points.  The physically correct
                  ! slopes occur at 0, while 1 is used for numerical closures [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: &
    KH_v, &       ! Isopycnal height diffusivities in v-columns [L2 T-1 ~> m2 s-1]
    int_slope_v   ! A nondimensional ratio from 0 to 1 that gives the relative
                  ! weighting of the interface slopes to that calculated also
                  ! using density gradients at v points.  The physically correct
                  ! slopes occur at 0, while 1 is used for numerical closures [nondim].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    KH_t          ! diagnosed diffusivity at tracer points [L2 T-1 ~> m2 s-1]

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    KH_u_CFL      ! The maximum stable isopycnal height diffusivity at u grid points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    KH_v_CFL      ! The maximum stable isopycnal height diffusivity at v grid points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    htot          ! The sum of the total layer thicknesses [H ~> m or kg m-2]
  real :: Khth_Loc_u(SZIB_(G),SZJ_(G)) ! The isopycnal height diffusivity at u points [L2 T-1 ~> m2 s-1]
  real :: Khth_Loc_v(SZI_(G),SZJB_(G)) ! The isopycnal height diffusivity at v points [L2 T-1 ~> m2 s-1]
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
  real, dimension(:,:), pointer :: cg1 => null() !< Wave speed [L T-1 ~> m s-1]
  real :: hu(SZI_(G),SZJ_(G))       ! A thickness-based mask at u points, used for diagnostics [nondim]
  real :: hv(SZI_(G),SZJ_(G))       ! A thickness-based mask at v points, used for diagnostics [nondim]
  real :: KH_u_lay(SZI_(G),SZJ_(G)) ! Diagnostic of isopycnal height diffusivities at u-points averaged
                                    ! to layer centers [L2 T-1 ~> m2 s-1]
  real :: KH_v_lay(SZI_(G),SZJ_(G)) ! Diagnostic of isopycnal height diffusivities at v-points averaged
                                    ! to layer centers [L2 T-1 ~> m2 s-1]
  logical :: use_VarMix, Resoln_scaled, Depth_scaled, use_stored_slopes, khth_use_vert_struct, use_Visbeck
  logical :: use_QG_Leith
  integer :: i, j, k, is, ie, js, je, nz

  if (.not. CS%initialized) call MOM_error(FATAL, "MOM_thickness_diffuse: "//&
         "Module must be initialized before it is used.")

  if ((.not.CS%thickness_diffuse) &
      .or. .not. (CS%Khth > 0.0 .or. CS%read_khth &
      .or. VarMix%use_variable_mixing)) return

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  h_neglect = GV%H_subroundoff

  if (allocated(MEKE%GM_src)) then
    do j=js,je ; do i=is,ie ; MEKE%GM_src(i,j) = 0. ; enddo ; enddo
  endif

  use_VarMix = .false. ; Resoln_scaled = .false. ; use_stored_slopes = .false.
  khth_use_vert_struct = .false. ; use_Visbeck = .false. ; use_QG_Leith = .false.
  Depth_scaled = .false.

  if (VarMix%use_variable_mixing) then
    use_VarMix = VarMix%use_variable_mixing .and. (CS%KHTH_Slope_Cff > 0.)
    Resoln_scaled = VarMix%Resoln_scaled_KhTh
    Depth_scaled = VarMix%Depth_scaled_KhTh
    use_stored_slopes = VarMix%use_stored_slopes
    khth_use_vert_struct = allocated(VarMix%khth_struct)
    use_Visbeck = VarMix%use_Visbeck
    use_QG_Leith = VarMix%use_QG_Leith_GM
    if (allocated(VarMix%cg1)) cg1 => VarMix%cg1
  else
    cg1 => null()
  endif


  !$OMP parallel do default(shared)
  do j=js,je ; do I=is-1,ie
    KH_u_CFL(I,j) = (0.25*CS%max_Khth_CFL) /  &
      (dt * ((G%IdxCu(I,j)*G%IdxCu(I,j)) + (G%IdyCu(I,j)*G%IdyCu(I,j))))
  enddo ; enddo
  !$OMP parallel do default(shared)
  do J=js-1,je ; do i=is,ie
    KH_v_CFL(i,J) = (0.25*CS%max_Khth_CFL) / &
      (dt * ((G%IdxCv(i,J)*G%IdxCv(i,J)) + (G%IdyCv(i,J)*G%IdyCv(i,J))))
  enddo ; enddo

  ! Calculates interface heights, e, in [Z ~> m].
  call find_eta(h, tv, G, GV, US, e, halo_size=1)

  ! Set the diffusivities.
  !$OMP parallel default(shared)
  if (.not. CS%read_khth) then
    !$OMP do
    do j=js,je ; do I=is-1,ie
      Khth_loc_u(I,j) = CS%Khth
    enddo ; enddo
  else ! use 2d KHTH that was read in from file
    !$OMP do
    do j=js,je ; do I=is-1,ie
      Khth_loc_u(I,j) = 0.5 * (CS%khth2d(i,j) + CS%khth2d(i+1,j))
    enddo ; enddo
  endif

  if (use_VarMix) then
    if (use_Visbeck) then
      !$OMP do
      do j=js,je ; do I=is-1,ie
        Khth_loc_u(I,j) = Khth_loc_u(I,j) + &
          CS%KHTH_Slope_Cff*VarMix%L2u(I,j) * VarMix%SN_u(I,j)
      enddo ; enddo
    endif
  endif

  if (allocated(MEKE%Kh)) then
    if (CS%MEKE_GEOMETRIC) then
      !$OMP do
      do j=js,je ; do I=is-1,ie
        Khth_loc_u(I,j) = Khth_loc_u(I,j) + G%OBCmaskCu(I,j) * CS%MEKE_GEOMETRIC_alpha * &
                          0.5*(MEKE%MEKE(i,j)+MEKE%MEKE(i+1,j)) / &
                          (VarMix%SN_u(I,j) + CS%MEKE_GEOMETRIC_epsilon)
      enddo ; enddo
    else
      do j=js,je ; do I=is-1,ie
        Khth_loc_u(I,j) = Khth_loc_u(I,j) + MEKE%KhTh_fac*sqrt(MEKE%Kh(i,j)*MEKE%Kh(i+1,j))
      enddo ; enddo
    endif
  endif

  if (Resoln_scaled) then
    !$OMP do
    do j=js,je ; do I=is-1,ie
      Khth_loc_u(I,j) = Khth_loc_u(I,j) * VarMix%Res_fn_u(I,j)
    enddo ; enddo
  endif

  if (Depth_scaled) then
    !$OMP do
    do j=js,je ; do I=is-1,ie
      Khth_loc_u(I,j) = Khth_loc_u(I,j) * VarMix%Depth_fn_u(I,j)
    enddo ; enddo
  endif

  if (CS%Khth_Max > 0) then
    !$OMP do
    do j=js,je ; do I=is-1,ie
      Khth_loc_u(I,j) = max(CS%Khth_Min, min(Khth_loc_u(I,j), CS%Khth_Max))
    enddo ; enddo
  else
    !$OMP do
    do j=js,je ; do I=is-1,ie
      Khth_loc_u(I,j) = max(CS%Khth_Min, Khth_loc_u(I,j))
    enddo ; enddo
  endif
  !$OMP do
  do j=js,je ; do I=is-1,ie
    KH_u(I,j,1) = min(KH_u_CFL(I,j), Khth_loc_u(I,j))
  enddo ; enddo

  if (khth_use_vert_struct) then
    if (CS%full_depth_khth_min) then
      !$OMP do
      do K=2,nz+1 ; do j=js,je ; do I=is-1,ie
        KH_u(I,j,K) = KH_u(I,j,1) * 0.5 * ( VarMix%khth_struct(i,j,k-1) + VarMix%khth_struct(i+1,j,k-1) )
        KH_u(I,j,K) = max(KH_u(I,j,K), CS%Khth_Min)
      enddo ; enddo ; enddo
    else
      !$OMP do
      do K=2,nz+1 ; do j=js,je ; do I=is-1,ie
        KH_u(I,j,K) = KH_u(I,j,1) * 0.5 * ( VarMix%khth_struct(i,j,k-1) + VarMix%khth_struct(i+1,j,k-1) )
      enddo ; enddo ; enddo
    endif
  else
    !$OMP do
    do K=2,nz+1 ; do j=js,je ; do I=is-1,ie
      KH_u(I,j,K) = KH_u(I,j,1)
    enddo ; enddo ; enddo
  endif

  if (use_VarMix) then
    if (use_QG_Leith) then
      !$OMP do
      do k=1,nz ; do j=js,je ; do I=is-1,ie
        KH_u(I,j,k) = VarMix%KH_u_QG(I,j,k)
      enddo ; enddo ; enddo
    endif
  endif

  if (CS%use_GME_thickness_diffuse) then
    !$OMP do
    do k=1,nz+1 ; do j=js,je ; do I=is-1,ie
      CS%KH_u_GME(I,j,k) = KH_u(I,j,k)
    enddo ; enddo ; enddo
  endif

  if (.not. CS%read_khth) then
   !$OMP do
    do J=js-1,je ; do i=is,ie
      Khth_loc_v(i,J) = CS%Khth
    enddo ; enddo
  else ! read KHTH from file
   !$OMP do
    do J=js-1,je ; do i=is,ie
      Khth_loc_v(i,J) = 0.5 * (CS%khth2d(i,j) + CS%khth2d(i,j+1))
    enddo ; enddo
  endif

  if (use_VarMix) then
    if (use_Visbeck) then
      !$OMP do
      do J=js-1,je ; do i=is,ie
        Khth_loc_v(i,J) = Khth_loc_v(i,J) + CS%KHTH_Slope_Cff*VarMix%L2v(i,J)*VarMix%SN_v(i,J)
      enddo ; enddo
    endif
  endif
  if (allocated(MEKE%Kh)) then
    if (CS%MEKE_GEOMETRIC) then
      !$OMP do
      do J=js-1,je ; do i=is,ie
        Khth_loc_v(i,J) = Khth_loc_v(i,J) + G%OBCmaskCv(i,J) * CS%MEKE_GEOMETRIC_alpha * &
                        0.5*(MEKE%MEKE(i,j)+MEKE%MEKE(i,j+1)) / &
                        (VarMix%SN_v(i,J) + CS%MEKE_GEOMETRIC_epsilon)
      enddo ; enddo
    else
      do J=js-1,je ; do i=is,ie
        Khth_loc_v(i,J) = Khth_loc_v(i,J) + MEKE%KhTh_fac*sqrt(MEKE%Kh(i,j)*MEKE%Kh(i,j+1))
      enddo ; enddo
    endif
  endif

  if (Resoln_scaled) then
    !$OMP do
    do J=js-1,je ; do i=is,ie
      Khth_loc_v(i,J) = Khth_loc_v(i,J) * VarMix%Res_fn_v(i,J)
    enddo ; enddo
  endif

  if (Depth_scaled) then
    !$OMP do
    do J=js-1,je ; do i=is,ie
      Khth_loc_v(i,J) = Khth_loc_v(i,J) * VarMix%Depth_fn_v(i,J)
    enddo ; enddo
  endif

  if (CS%Khth_Max > 0) then
    !$OMP do
    do J=js-1,je ; do i=is,ie
      Khth_loc_v(i,J) = max(CS%Khth_Min, min(Khth_loc_v(i,J), CS%Khth_Max))
    enddo ; enddo
  else
    !$OMP do
    do J=js-1,je ; do i=is,ie
      Khth_loc_v(i,J) = max(CS%Khth_Min, Khth_loc_v(i,J))
    enddo ; enddo
  endif

  if (CS%max_Khth_CFL > 0.0) then
    !$OMP do
    do J=js-1,je ; do i=is,ie
      KH_v(i,J,1) = min(KH_v_CFL(i,J), Khth_loc_v(i,J))
    enddo ; enddo
  endif

  if (khth_use_vert_struct) then
      if (CS%full_depth_khth_min) then
      !$OMP do
      do K=2,nz+1 ; do J=js-1,je ; do i=is,ie
        KH_v(i,J,K) = KH_v(i,J,1) * 0.5 * ( VarMix%khth_struct(i,j,k-1) + VarMix%khth_struct(i,j+1,k-1) )
        KH_v(i,J,K) = max(KH_v(i,J,K), CS%Khth_Min)
      enddo ; enddo ; enddo
    else
      !$OMP do
      do K=2,nz+1 ; do J=js-1,je ; do i=is,ie
        KH_v(i,J,K) = KH_v(i,J,1) * 0.5 * ( VarMix%khth_struct(i,j,k-1) + VarMix%khth_struct(i,j+1,k-1) )
      enddo ; enddo ; enddo
    endif
  else
    !$OMP do
    do K=2,nz+1 ; do J=js-1,je ; do i=is,ie
      KH_v(i,J,K) = KH_v(i,J,1)
    enddo ; enddo ; enddo
  endif

  if (use_VarMix) then
    if (use_QG_Leith) then
      !$OMP do
      do k=1,nz ; do J=js-1,je ; do i=is,ie
        KH_v(i,J,k) = VarMix%KH_v_QG(i,J,k)
      enddo ; enddo ; enddo
    endif
  endif

  if (CS%use_GME_thickness_diffuse) then
    !$OMP do
    do k=1,nz+1 ; do J=js-1,je ; do i=is,ie
      CS%KH_v_GME(i,J,k) = KH_v(i,J,k)
    enddo ; enddo ; enddo
  endif

  if (allocated(MEKE%Kh)) then
    if (CS%MEKE_GEOMETRIC) then
      if (CS%MEKE_GEOM_answer_date < 20190101) then
        !$OMP do
        do j=js,je ; do i=is,ie
          ! This does not give bitwise rotational symmetry.
          MEKE%Kh(i,j) = CS%MEKE_GEOMETRIC_alpha * MEKE%MEKE(i,j) / &
                         (0.25*(VarMix%SN_u(I,j)+VarMix%SN_u(I-1,j) + &
                                VarMix%SN_v(i,J)+VarMix%SN_v(i,J-1)) + &
                          CS%MEKE_GEOMETRIC_epsilon)
        enddo ; enddo
      else
        !$OMP do
        do j=js,je ; do i=is,ie
          ! With the additional parentheses this gives bitwise rotational symmetry.
          MEKE%Kh(i,j) = CS%MEKE_GEOMETRIC_alpha * MEKE%MEKE(i,j) / &
                         (0.25*((VarMix%SN_u(I,j)+VarMix%SN_u(I-1,j)) + &
                                (VarMix%SN_v(i,J)+VarMix%SN_v(i,J-1))) + &
                          CS%MEKE_GEOMETRIC_epsilon)
        enddo ; enddo
      endif
    endif
  endif

  !$OMP do
  do K=1,nz+1 ; do j=js,je ; do I=is-1,ie ; int_slope_u(I,j,K) = 0.0 ; enddo ; enddo ; enddo
  !$OMP do
  do K=1,nz+1 ; do J=js-1,je ; do i=is,ie ; int_slope_v(i,J,K) = 0.0 ; enddo ; enddo ; enddo
  !$OMP end parallel

  if (CS%detangle_interfaces) then
    call add_detangling_Kh(h, e, Kh_u, Kh_v, KH_u_CFL, KH_v_CFL, tv, dt, G, GV, US, &
                           CS, int_slope_u, int_slope_v)
  endif

  if ((CS%Kh_eta_bg > 0.0) .or. (CS%Kh_eta_vel > 0.0)) then
    call add_interface_Kh(G, GV, US, CS, Kh_u, Kh_v, KH_u_CFL, KH_v_CFL, int_slope_u, int_slope_v)
  endif

  if (CS%debug) then
    call uvchksum("Kh_[uv]", Kh_u, Kh_v, G%HI, haloshift=0, &
                  unscale=(US%L_to_m**2)*US%s_to_T, scalar_pair=.true.)
    call uvchksum("Kh_[uv]_CFL", Kh_u_CFL, Kh_v_CFL, G%HI, haloshift=0, &
                  unscale=(US%L_to_m**2)*US%s_to_T, scalar_pair=.true.)
    if (Resoln_scaled) then
      call uvchksum("Res_fn_[uv]", VarMix%Res_fn_u, VarMix%Res_fn_v, G%HI, haloshift=0, &
                    unscale=1.0, scalar_pair=.true.)
    endif
    call uvchksum("int_slope_[uv]", int_slope_u, int_slope_v, G%HI, haloshift=0)
    call hchksum(h, "thickness_diffuse_1 h", G%HI, haloshift=1, unscale=GV%H_to_m)
    call hchksum(e, "thickness_diffuse_1 e", G%HI, haloshift=1, unscale=US%Z_to_m)
    if (use_stored_slopes) then
      call uvchksum("VarMix%slope_[xy]", VarMix%slope_x, VarMix%slope_y, &
                    G%HI, haloshift=0, unscale=US%Z_to_L)
    endif
    if (associated(tv%eqn_of_state)) then
      call hchksum(tv%T, "thickness_diffuse T", G%HI, haloshift=1, unscale=US%C_to_degC)
      call hchksum(tv%S, "thickness_diffuse S", G%HI, haloshift=1, unscale=US%S_to_ppt)
    endif
  endif

  ! Calculate uhD, vhD from h, e, KH_u, KH_v, tv%T/S
  if (STOCH%skeb_use_gm) then
    if (use_stored_slopes) then
      call thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, cg1, dt, G, GV, US, MEKE, CS, &
                                  int_slope_u, int_slope_v, VarMix%slope_x, VarMix%slope_y, &
                                  STOCH=STOCH, VarMix=VarMix)
    else
      call thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, cg1, dt, G, GV, US, MEKE, CS, &
                                  int_slope_u, int_slope_v, STOCH=STOCH, VarMix=VarMix)
    endif
  else
    if (use_stored_slopes) then
      call thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, cg1, dt, G, GV, US, MEKE, CS, &
                                  int_slope_u, int_slope_v, VarMix%slope_x, VarMix%slope_y)
    else
      call thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, cg1, dt, G, GV, US, MEKE, CS, &
                                  int_slope_u, int_slope_v)
    endif
  endif

  if (VarMix%use_variable_mixing) then
    if (allocated(MEKE%Rd_dx_h) .and. allocated(VarMix%Rd_dx_h)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        MEKE%Rd_dx_h(i,j) = VarMix%Rd_dx_h(i,j)
      enddo ; enddo
    endif
  endif

  ! offer diagnostic fields for averaging
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_uhGM > 0)   call post_data(CS%id_uhGM, uhD, CS%diag)
    if (CS%id_vhGM > 0)   call post_data(CS%id_vhGM, vhD, CS%diag)
    if (CS%id_GMwork > 0) call post_data(CS%id_GMwork, CS%GMwork, CS%diag)
    if (CS%id_KH_u > 0)   call post_data(CS%id_KH_u, KH_u, CS%diag)
    if (CS%id_KH_v > 0)   call post_data(CS%id_KH_v, KH_v, CS%diag)
    if (CS%id_KH_u1 > 0)  call post_data(CS%id_KH_u1, KH_u(:,:,1), CS%diag)
    if (CS%id_KH_v1 > 0)  call post_data(CS%id_KH_v1, KH_v(:,:,1), CS%diag)

    ! Diagnose diffusivity at T-cell point.  Do a simple average, rather than a
    ! thickness-weighted average, so that KH_t is depth-independent when KH_u and KH_v
    ! are depth independent.  If a thickness-weighted average were used, the variations
    ! of thickness could give a spurious depth dependence to the diagnosed KH_t.
    if (CS%id_KH_t > 0 .or. CS%id_KH_t1 > 0 .or. CS%Use_KH_in_MEKE) then
      do k=1,nz
        ! thicknesses across u and v faces, converted to 0/1 mask
        ! layer average of the interface diffusivities KH_u and KH_v
        do j=js,je ; do I=is-1,ie
          ! This expression uses harmonic mean thicknesses:
          ! hu(I,j)       = 2.0*h(i,j,k)*h(i+1,j,k) / (h(i,j,k)+h(i+1,j,k)+h_neglect)
          ! This expression is a 0/1 mask based on depths where there are thick layers:
          hu(I,j) = 0.0 ; if (h(i,j,k)*h(i+1,j,k) /= 0.0) hu(I,j) = 1.0
          KH_u_lay(I,j) = 0.5*(KH_u(I,j,k)+KH_u(I,j,k+1))
        enddo ; enddo
        do J=js-1,je ; do i=is,ie
          ! This expression uses harmonic mean thicknesses:
          ! hv(i,J)       = 2.0*h(i,j,k)*h(i,j+1,k)/(h(i,j,k)+h(i,j+1,k)+h_neglect)
          ! This expression is a 0/1 mask based on depths where there are thick layers:
          hv(i,J) = 0.0 ; if (h(i,j,k)*h(i,j+1,k) /= 0.0) hv(i,J) = 1.0
          KH_v_lay(i,J) = 0.5*(KH_v(i,J,k)+KH_v(i,J,k+1))
        enddo ; enddo
        ! diagnose diffusivity at T-points
        do j=js,je ; do i=is,ie
          Kh_t(i,j,k) = (((hu(I-1,j)*KH_u_lay(i-1,j)) + (hu(I,j)*KH_u_lay(I,j))) + &
                         ((hv(i,J-1)*KH_v_lay(i,J-1)) + (hv(i,J)*KH_v_lay(i,J)))) / &
                        ((hu(I-1,j)+hu(I,j)) + (hv(i,J-1)+hv(i,J)) + 1.0e-20)
          ! Use this denominator instead if hu and hv are actual thicknesses rather than a 0/1 mask:
          !              ((hu(I-1,j)+hu(I,j)) + (hv(i,J-1)+hv(i,J)) + h_neglect)
        enddo ; enddo
      enddo

      if (CS%Use_KH_in_MEKE) then
        MEKE%Kh_diff(:,:) = 0.0
        htot(:,:) = 0.0
        do k=1,nz
          do j=js,je ; do i=is,ie
            MEKE%Kh_diff(i,j) = MEKE%Kh_diff(i,j) + Kh_t(i,j,k) * h(i,j,k)
            htot(i,j) = htot(i,j) + h(i,j,k)
          enddo ; enddo
        enddo

        do j=js,je ; do i=is,ie
          MEKE%Kh_diff(i,j) = MEKE%Kh_diff(i,j) / MAX(CS%MEKE_min_depth_diff, htot(i,j))
        enddo ; enddo
      endif

      if (CS%id_KH_t  > 0) call post_data(CS%id_KH_t,  KH_t,        CS%diag)
      if (CS%id_KH_t1 > 0) call post_data(CS%id_KH_t1, KH_t(:,:,1), CS%diag)
    endif

  endif

  !$OMP parallel do default(shared)
  do k=1,nz
    do j=js,je ; do I=is-1,ie
      uhtr(I,j,k) = uhtr(I,j,k) + uhD(I,j,k) * dt
      if (associated(CDp%uhGM)) CDp%uhGM(I,j,k) = uhD(I,j,k)
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      vhtr(i,J,k) = vhtr(i,J,k) + vhD(i,J,k) * dt
      if (associated(CDp%vhGM)) CDp%vhGM(i,J,k) = vhD(i,J,k)
    enddo ; enddo
    do j=js,je ; do i=is,ie
      h(i,j,k) = h(i,j,k) - dt * G%IareaT(i,j) * &
          ((uhD(I,j,k) - uhD(I-1,j,k)) + (vhD(i,J,k) - vhD(i,J-1,k)))
      if (h(i,j,k) < GV%Angstrom_H) h(i,j,k) = GV%Angstrom_H
    enddo ; enddo
  enddo

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  ! This needs to happen after the H update and before the next post_data.
  call diag_update_remap_grids(CS%diag)

  if (CS%debug) then
    call uvchksum("thickness_diffuse [uv]hD", uhD, vhD, &
                  G%HI, haloshift=0, unscale=GV%H_to_m*US%L_to_m**2*US%s_to_T)
    call uvchksum("thickness_diffuse [uv]htr", uhtr, vhtr, &
                  G%HI, haloshift=0, unscale=US%L_to_m**2*GV%H_to_m)
    call hchksum(h, "thickness_diffuse h", G%HI, haloshift=0, unscale=GV%H_to_m)
  endif

end subroutine thickness_diffuse

!> Calculates parameterized layer transports for use in the continuity equation.
!! Fluxes are limited to give positive definite thicknesses.
!! Called by thickness_diffuse().
subroutine thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, cg1, dt, G, GV, US, MEKE, &
                                  CS, int_slope_u, int_slope_v, slope_x, slope_y, STOCH, VarMix)
  type(ocean_grid_type),                        intent(in)  :: G     !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)  :: GV    !< Vertical grid structure
  type(unit_scale_type),                        intent(in)  :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)  :: h     !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(in)  :: e     !< Interface positions [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in)  :: Kh_u  !< Isopycnal height diffusivity
                                                                     !! at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in)  :: Kh_v  !< Isopycnal height diffusivity
                                                                     !! at v points [L2 T-1 ~> m2 s-1]
  type(thermo_var_ptrs),                        intent(in)  :: tv    !< Thermodynamics structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),   intent(out) :: uhD   !< Zonal mass fluxes
                                                                     !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),   intent(out) :: vhD   !< Meridional mass fluxes
                                                                     !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(:,:),                         pointer     :: cg1   !< Wave speed [L T-1 ~> m s-1]
  real,                                         intent(in)  :: dt    !< Time increment [T ~> s]
  type(MEKE_type),                              intent(inout) :: MEKE !< MEKE fields
  type(thickness_diffuse_CS),                   intent(inout) :: CS  !< Control structure for thickness_diffuse
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in)  :: int_slope_u !< Ratio that determine how much of
                                                                     !! the isopycnal slopes are taken directly from
                                                                     !! the interface slopes without consideration of
                                                                     !! density gradients [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in)  :: int_slope_v !< Ratio that determine how much of
                                                                     !! the isopycnal slopes are taken directly from
                                                                     !! the interface slopes without consideration of
                                                                     !! density gradients [nondim].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), optional, intent(in)  :: slope_x !< Isopyc. slope at u [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), optional, intent(in)  :: slope_y !< Isopyc. slope at v [Z L-1 ~> nondim]
  type(stochastic_CS),                       optional, intent(inout)  :: STOCH !< Stochastic control structure
  type(VarMix_CS), target,                      optional, intent(in)  :: VarMix !< Variable mixing coefficents

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    T, &          ! The temperature [C ~> degC], with the values in
                  ! in massless layers filled vertically by diffusion.
    S, &          ! The filled salinity [S ~> ppt], with the values in
                  ! in massless layers filled vertically by diffusion.
    h_avail, &    ! The mass available for diffusion out of each face, divided
                  ! by dt [H L2 T-1 ~> m3 s-1 or kg s-1].
    h_frac        ! The fraction of the mass in the column above the bottom
                  ! interface of a layer that is within a layer [nondim]. 0<h_frac<=1
  real :: dz(SZI_(G),SZJ_(G),SZK_(GV)) ! Height change across layers [Z ~> m]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: &
    Slope_y_PE, &  ! 3D array of neutral slopes at v-points, set equal to Slope (below) [Z L-1 ~> nondim]
    hN2_y_PE       ! Harmonic mean of thicknesses around the interfaces times the buoyancy frequency
                   ! at v-points with unit conversion factors [H L2 Z-2 T-2 ~> m s-2 or kg m-2 s-2],
                   ! used for calculating the potential energy release
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: &
    Slope_x_PE, &  ! 3D array of neutral slopes at u-points, set equal to Slope (below) [Z L-1 ~> nondim]
    hN2_x_PE       ! Harmonic mean of thicknesses around the interfaces times the buoyancy frequency
                   ! at u-points  with unit conversion factors [H L2 Z-2 T-2 ~> m s-2 or kg m-2 s-2],
                   ! used for calculating the potential energy release
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    pres, &       ! The pressure at an interface [R L2 T-2 ~> Pa].
    h_avail_rsum  ! The running sum of h_avail above an interface [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G)) :: &
    drho_dT_u, &  ! The derivative of density with temperature at u points [R C-1 ~> kg m-3 degC-1]
    drho_dS_u     ! The derivative of density with salinity at u points [R S-1 ~> kg m-3 ppt-1].
  real, dimension(SZIB_(G)) :: scrap ! An array to pass to calculate_density_second_derivs()
                  ! with various units that will be ignored [various]
  real, dimension(SZI_(G)) :: &
    drho_dT_v, &  ! The derivative of density with temperature at v points [R C-1 ~> kg m-3 degC-1]
    drho_dS_v, &  ! The derivative of density with salinity at v points [R S-1 ~> kg m-3 ppt-1].
    drho_dT_dT_h, & ! The second derivative of density with temperature at h points [R C-2 ~> kg m-3 degC-2]
    drho_dT_dT_hr ! The second derivative of density with temperature at h (+1) points [R C-2 ~> kg m-3 degC-2]
  real :: uhtot(SZIB_(G),SZJ_(G))  ! The vertical sum of uhD [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: vhtot(SZI_(G),SZJB_(G))  ! The vertical sum of vhD [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G)) :: &
    T_u, &        ! Temperature on the interface at the u-point [C ~> degC].
    S_u, &        ! Salinity on the interface at the u-point [S ~> ppt].
    pres_u        ! Pressure on the interface at the u-point [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G)) :: &
    T_v, &        ! Temperature on the interface at the v-point [C ~> degC].
    S_v, &        ! Salinity on the interface at the v-point [S ~> ppt].
    pres_v, &     ! Pressure on the interface at the v-point [R L2 T-2 ~> Pa].
    T_h, &        ! Temperature on the interface at the h-point [C ~> degC].
    S_h, &        ! Salinity on the interface at the h-point [S ~> ppt].
    pres_h, &     ! Pressure on the interface at the h-point [R L2 T-2 ~> Pa].
    T_hr, &       ! Temperature on the interface at the h (+1) point [C ~> degC].
    S_hr, &       ! Salinity on the interface at the h (+1) point [S ~> ppt].
    pres_hr       ! Pressure on the interface at the h (+1) point [R L2 T-2 ~> Pa].
  real :: Work_u(SZIB_(G),SZJ_(G)) ! The work done by the isopycnal height diffusion
                                   ! integrated over u-point water columns [R Z L4 T-3 ~> W]
  real :: Work_v(SZI_(G),SZJB_(G)) ! The work done by the isopycnal height diffusion
                                   ! integrated over v-point water columns [R Z L4 T-3 ~> W]
  real :: Work_h        ! The work averaged over an h-cell [R Z L2 T-3 ~> W m-2].
  real :: PE_release_h  ! The amount of potential energy released by GM averaged over an h-cell
                        ! [R Z L2 T-3 ~> W m-2].  The calculation equals rho0 * h * S^2 * N^2 * kappa_GM.
  real :: I4dt          ! 1 / 4 dt [T-1 ~> s-1].
  real :: drdiA, drdiB  ! Along layer zonal potential density  gradients in the layers above (A)
                        ! and below (B) the interface times the grid spacing [R ~> kg m-3].
  real :: drdjA, drdjB  ! Along layer meridional potential density  gradients in the layers above (A)
                        ! and below (B) the interface times the grid spacing [R ~> kg m-3].
  real :: drdkL, drdkR  ! Vertical density differences across an interface [R ~> kg m-3].
  real :: drdi_u(SZIB_(G),SZK_(GV)) ! Copy of drdi at u-points [R ~> kg m-3].
  real :: drdj_v(SZI_(G),SZK_(GV)) ! Copy of drdj at v-points [R ~> kg m-3].
  real :: drdkDe_u(SZIB_(G),SZK_(GV)+1) ! Lateral difference of product of drdk and e at u-points
                                        ! [Z R ~> kg m-2].
  real :: drdkDe_v(SZI_(G),SZK_(GV)+1)  ! Lateral difference of product of drdk and e at v-points
                                        ! [Z R ~> kg m-2].
  real :: hg2A, hg2B, hg2L, hg2R ! Squares of geometric mean thicknesses [H2 ~> m2 or kg2 m-4].
  real :: haA, haB, haL, haR     ! Arithmetic mean thicknesses [H ~> m or kg m-2].
  real :: dzg2A, dzg2B  ! Squares of geometric mean vertical layer extents [Z2 ~> m2].
  real :: dzaA, dzaB    ! Arithmetic mean vertical layer extents [Z ~> m].
  real :: dzaL, dzaR    ! Temporary vertical layer extents [Z ~> m]
  real :: wtA, wtB      ! Unnormalized weights of the slopes above and below [H3 ~> m3 or kg3 m-6]
  real :: wtL, wtR      ! Unnormalized weights of the slopes to the left and right [H3 Z ~> m4 or kg3 m-5]
  real :: drdx, drdy    ! Zonal and meridional density gradients [R L-1 ~> kg m-4].
  real :: drdz          ! Vertical density gradient [R Z-1 ~> kg m-4].
  real :: dz_harm       ! Harmonic mean layer vertical extent [Z ~> m].
  real :: c2_dz_u(SZIB_(G),SZK_(GV)+1) ! Wave speed squared divided by dz at u-points [L2 Z-1 T-2 ~> m s-2]
  real :: c2_dz_v(SZI_(G),SZK_(GV)+1)  ! Wave speed squared divided by dz at v-points [L2 Z-1 T-2 ~> m s-2]
  real :: dzN2_u(SZIB_(G),SZK_(GV)+1) ! Vertical extent times N2 at interfaces above u-points times
                        ! rescaling factors from vertical to horizontal distances [L2 Z-1 T-2 ~> m s-2]
  real :: dzN2_v(SZI_(G),SZK_(GV)+1)  ! Vertical extent times N2 at interfaces above v-points times
                        ! rescaling factors from vertical to horizontal distances [L2 Z-1 T-2 ~> m s-2]
  real :: Sfn_est       ! A preliminary estimate (before limiting) of the overturning
                        ! streamfunction [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: Sfn_unlim_u(SZIB_(G),SZK_(GV)+1) ! Volume streamfunction for u-points [Z L2 T-1 ~> m3 s-1]
  real :: Sfn_unlim_v(SZI_(G),SZK_(GV)+1)  ! Volume streamfunction for v-points [Z L2 T-1 ~> m3 s-1]
  real :: slope2_Ratio_u(SZIB_(G),SZK_(GV)+1) ! The ratio of the slope squared to slope_max squared [nondim]
  real :: slope2_Ratio_v(SZI_(G),SZK_(GV)+1)  ! The ratio of the slope squared to slope_max squared [nondim]
  real :: Sfn_in_h      ! The overturning streamfunction [H L2 T-1 ~> m3 s-1 or kg s-1] (note that
                        ! the units are different from other Sfn vars).
  real :: Sfn_safe      ! The streamfunction that goes linearly back to 0 at the surface
                        ! [H L2 T-1 ~> m3 s-1 or kg s-1].  This is a good value to use when the
                        ! slope is so large as to be meaningless, usually due to weak stratification.
  real :: Slope         ! The slope of density surfaces, calculated in a way that is always
                        ! between -1 and 1 after undoing dimensional scaling, [Z L-1 ~> nondim]
  real :: mag_grad2     ! The squared magnitude of the 3-d density gradient [R2 L-2 ~> kg2 m-8].
  real :: I_slope_max2  ! The inverse of slope_max squared [L2 Z-2 ~> nondim].
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: hn_2          ! Half of h_neglect [H ~> m or kg m-2].
  real :: h_neglect2    ! h_neglect^2 [H2 ~> m2 or kg2 m-4].
  real :: dz_neglect    ! A thickness [Z ~> m], that is so small it is usually lost
                        ! in roundoff and can be neglected [Z ~> m].
  real :: dz_neglect2   ! dz_neglect^2 [Z2 ~> m2]
  real :: G_scale       ! The gravitational acceleration times a unit conversion
                        ! factor [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  logical :: use_EOS    ! If true, density is calculated from T & S using an equation of state.
  logical :: find_work  ! If true, find the change in energy due to the fluxes.
  integer :: nk_linear  ! The number of layers over which the streamfunction goes to 0.
  real :: G_rho0        ! g/Rho0 [L2 R-1 Z-1 T-2 ~> m4 kg-1 s-2].
  real :: Rho_avg       ! The in situ density averaged to an interface [R ~> kg m-3]
  real :: N2_floor      ! A floor for N2 to avoid degeneracy in the elliptic solver
                        ! times unit conversion factors [L2 Z-2 T-2 ~> s-2]
  real :: N2_unlim      ! An unlimited estimate of the buoyancy frequency
                        ! times unit conversion factors [L2 Z-2 T-2 ~> s-2]
  real :: Tl(5)         ! copy of T in local stencil [C ~> degC]
  real :: mn_T          ! mean of T in local stencil [C ~> degC]
  real :: mn_T2         ! mean of T**2 in local stencil [C2 ~> degC2]
  real :: hl(5)         ! Copy of local stencil of H [H ~> m]
  real :: r_sm_H        ! Reciprocal of sum of H in local stencil [H-1 ~> m-1]
  real :: Z_to_H        ! A conversion factor from heights to thicknesses, perhaps based on
                        ! a spatially variable local density [H Z-1 ~> nondim or kg m-3]
  real :: Tsgs2(SZI_(G),SZJ_(G),SZK_(GV)) ! Sub-grid temperature variance [C2 ~> degC2]
  real :: diag_sfn_x(SZIB_(G),SZJ_(G),SZK_(GV)+1)       ! Diagnostic of the x-face streamfunction
                                                        ! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: diag_sfn_unlim_x(SZIB_(G),SZJ_(G),SZK_(GV)+1) ! Diagnostic of the x-face streamfunction before
                                                        ! applying limiters [Z L2 T-1 ~> m3 s-1]
  real :: diag_sfn_y(SZI_(G),SZJB_(G),SZK_(GV)+1)       ! Diagnostic of the y-face streamfunction
                                                        ! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: diag_sfn_unlim_y(SZI_(G),SZJB_(G),SZK_(GV)+1) ! Diagnostic of the y-face streamfunction before
                                                        ! applying limiters [Z L2 T-1 ~> m3 s-1]
                                                        ! applying limiters [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, allocatable :: skeb_gm_work(:,:)                ! Temp array to hold GM work for SKEB
  real, allocatable :: skeb_ebt_norm2(:,:)              ! Used to normalize EBT for SKEB

  logical :: present_slope_x, present_slope_y, calc_derivatives
  integer, dimension(2) :: EOSdom_u  ! The shifted I-computational domain to use for equation of
                                     ! state calculations at u-points.
  integer, dimension(2) :: EOSdom_v  ! The shifted i-computational domain to use for equation of
                                     ! state calculations at v-points.
  integer, dimension(2) :: EOSdom_h1 ! The shifted i-computational domain to use for equation of
                                     ! state calculations at h points with 1 extra halo point
  logical :: use_stanley, skeb_use_gm
  integer :: is, ie, js, je, nz, IsdB, halo
  integer :: i, j, k
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke ; IsdB = G%IsdB

  I4dt = 0.25 / dt
  I_slope_max2 = 1.0 / (CS%slope_max**2)

  h_neglect = GV%H_subroundoff ; h_neglect2 = h_neglect**2 ; hn_2 = 0.5*h_neglect
  dz_neglect = GV%dZ_subroundoff ; dz_neglect2 = dz_neglect**2
  if (GV%Boussinesq) G_rho0 = GV%g_Earth / GV%Rho0
  N2_floor = CS%N2_floor

  use_EOS = associated(tv%eqn_of_state)
  present_slope_x = PRESENT(slope_x)
  present_slope_y = PRESENT(slope_y)

  use_stanley = CS%use_stanley_gm

  skeb_use_gm = .false.
  if (present(STOCH)) skeb_use_gm = STOCH%skeb_use_gm
  if (skeb_use_gm) then
    allocate(skeb_gm_work(is:ie,js:je), source=0.)
    allocate(skeb_ebt_norm2(is:ie,js:je), source=0.)
  endif

  nk_linear = max(GV%nkml, 1)

  Slope_x_PE(:,:,:) = 0.0
  Slope_y_PE(:,:,:) = 0.0
  hN2_x_PE(:,:,:) = 0.0
  hN2_y_PE(:,:,:) = 0.0

  find_work = allocated(MEKE%GM_src)
  find_work = (allocated(CS%GMwork) .or. find_work)
  find_work = (skeb_use_gm .or. find_work)

  if (use_EOS) then
    halo = 1 ! Default halo to fill is 1
    call vert_fill_TS(h, tv%T, tv%S, CS%kappa_smooth*dt, T, S, G, GV, US, halo, larger_h_denom=.true.)
  endif

  ! Rescale the thicknesses, perhaps using the specific volume.
  call thickness_to_dz(h, tv, dz, G, GV, US, halo_size=1)

  if (CS%use_FGNV_streamfn .and. .not. associated(cg1)) call MOM_error(FATAL, &
       "cg1 must be associated when using FGNV streamfunction.")

  !$OMP parallel default(shared) private(hl,r_sm_H,Tl,mn_T,mn_T2)
  ! Find the maximum and minimum permitted streamfunction.
  !$OMP do
  do j=js-1,je+1 ; do i=is-1,ie+1
    h_avail_rsum(i,j,1) = 0.0
    pres(i,j,1) = 0.0
    if (associated(tv%p_surf)) then ; pres(i,j,1) = tv%p_surf(i,j) ; endif

    h_avail(i,j,1) = max(I4dt*G%areaT(i,j)*(h(i,j,1)-GV%Angstrom_H),0.0)
    h_avail_rsum(i,j,2) = h_avail(i,j,1)
    h_frac(i,j,1) = 1.0
    pres(i,j,2) = pres(i,j,1) + (GV%g_Earth*GV%H_to_RZ) * h(i,j,1)
  enddo ; enddo
  do j=js-1,je+1
    do k=2,nz ; do i=is-1,ie+1
      h_avail(i,j,k) = max(I4dt*G%areaT(i,j)*(h(i,j,k)-GV%Angstrom_H),0.0)
      h_avail_rsum(i,j,k+1) = h_avail_rsum(i,j,k) + h_avail(i,j,k)
      h_frac(i,j,k) = 0.0 ; if (h_avail(i,j,k) > 0.0) &
        h_frac(i,j,k) = h_avail(i,j,k) / h_avail_rsum(i,j,k+1)
      pres(i,j,K+1) = pres(i,j,K) + (GV%g_Earth*GV%H_to_RZ) * h(i,j,k)
    enddo ; enddo
  enddo
  !$OMP do
  do j=js,je ; do I=is-1,ie
    uhtot(I,j) = 0.0 ; Work_u(I,j) = 0.0
  enddo ; enddo
  !$OMP do
  do J=js-1,je ; do i=is,ie
    vhtot(i,J) = 0.0 ; Work_v(i,J) = 0.0
  enddo ; enddo
  !$OMP end parallel

  if (CS%id_sfn_x > 0) then ; diag_sfn_x(:,:,1) = 0.0 ; diag_sfn_x(:,:,nz+1) = 0.0 ; endif
  if (CS%id_sfn_y > 0) then ; diag_sfn_y(:,:,1) = 0.0 ; diag_sfn_y(:,:,nz+1) = 0.0 ; endif
  if (CS%id_sfn_unlim_x > 0) then ; diag_sfn_unlim_x(:,:,1) = 0.0 ; diag_sfn_unlim_x(:,:,nz+1) = 0.0 ; endif
  if (CS%id_sfn_unlim_y > 0) then ; diag_sfn_unlim_y(:,:,1) = 0.0 ; diag_sfn_unlim_y(:,:,nz+1) = 0.0 ; endif

  EOSdom_u(1) = (is-1) - (G%IsdB-1) ; EOSdom_u(2) = ie - (G%IsdB-1)
  EOSdom_v(:) = EOS_domain(G%HI)
  EOSdom_h1(:) = EOS_domain(G%HI, halo=1)

  !$OMP parallel do default(none) shared(nz,is,ie,js,je,find_work,use_EOS,G,GV,US,pres,T,S, &
  !$OMP                                  nk_linear,IsdB,tv,h,h_neglect,e,dz,dz_neglect,dz_neglect2, &
  !$OMP                                  h_neglect2,hn_2,I_slope_max2,int_slope_u,KH_u,uhtot, &
  !$OMP                                  h_frac,h_avail_rsum,uhD,h_avail,Work_u,CS,slope_x,cg1, &
  !$OMP                                  diag_sfn_x,diag_sfn_unlim_x,N2_floor,EOSdom_u,EOSdom_h1, &
  !$OMP                                  use_stanley,Tsgs2,present_slope_x,G_rho0,Slope_x_PE,hN2_x_PE) &
  !$OMP                          private(drdiA,drdiB,drdkL,drdkR,pres_u,T_u,S_u,G_scale, &
  !$OMP                                  drho_dT_u,drho_dS_u,hg2A,hg2B,hg2L,hg2R,haA, &
  !$OMP                                  drho_dT_dT_h,scrap,pres_h,T_h,S_h,N2_unlim,  &
  !$OMP                                  haB,haL,haR,dzaL,dzaR,wtA,wtB,wtL,wtR,drdz,  &
  !$OMP                                  dzg2A,dzg2B,dzaA,dzaB,dz_harm,Z_to_H, &
  !$OMP                                  drdx,mag_grad2,Slope,slope2_Ratio_u,dzN2_u,  &
  !$OMP                                  Sfn_unlim_u,Rho_avg,drdi_u,drdkDe_u,c2_dz_u, &
  !$OMP                                  Sfn_safe,Sfn_est,Sfn_in_h,calc_derivatives)
  do j=js,je
    do I=is-1,ie ; dzN2_u(I,1) = 0. ; dzN2_u(I,nz+1) = 0. ; enddo
    do K=nz,2,-1
      if (find_work .and. .not.(use_EOS)) then
        drdiA = 0.0 ; drdiB = 0.0
        drdkL = GV%Rlay(k) - GV%Rlay(k-1) ; drdkR = drdkL
      endif

      calc_derivatives = use_EOS .and. (k >= nk_linear) .and. &
         (find_work .or. .not. present_slope_x .or. CS%use_FGNV_streamfn .or. use_stanley)

      ! Calculate the zonal fluxes and gradients.
      if (calc_derivatives) then
        do I=is-1,ie
          pres_u(I) = 0.5*(pres(i,j,K) + pres(i+1,j,K))
          T_u(I) = 0.25*((T(i,j,k) + T(i+1,j,k)) + (T(i,j,k-1) + T(i+1,j,k-1)))
          S_u(I) = 0.25*((S(i,j,k) + S(i+1,j,k)) + (S(i,j,k-1) + S(i+1,j,k-1)))
        enddo
        call calculate_density_derivs(T_u, S_u, pres_u, drho_dT_u, drho_dS_u, &
                                      tv%eqn_of_state, EOSdom_u)
      endif
      if (use_stanley) then
        do i=is-1,ie+1
          pres_h(i) = pres(i,j,K)
          T_h(i) = 0.5*(T(i,j,k) + T(i,j,k-1))
          S_h(i) = 0.5*(S(i,j,k) + S(i,j,k-1))
        enddo

        ! The second line below would correspond to arguments
        !            drho_dS_dS, drho_dS_dT, drho_dT_dT, drho_dS_dP, drho_dT_dP, &
        call calculate_density_second_derivs(T_h, S_h, pres_h, &
                     scrap, scrap, drho_dT_dT_h, scrap, scrap, &
                     tv%eqn_of_state, EOSdom_h1)
      endif

      do I=is-1,ie
        if (calc_derivatives) then
          ! Estimate the horizontal density gradients along layers.
          drdiA = drho_dT_u(I) * (T(i+1,j,k-1)-T(i,j,k-1)) + &
                  drho_dS_u(I) * (S(i+1,j,k-1)-S(i,j,k-1))
          drdiB = drho_dT_u(I) * (T(i+1,j,k)-T(i,j,k)) + &
                  drho_dS_u(I) * (S(i+1,j,k)-S(i,j,k))

          ! Estimate the vertical density gradients times the grid spacing.
          drdkL = (drho_dT_u(I) * (T(i,j,k)-T(i,j,k-1)) + &
                   drho_dS_u(I) * (S(i,j,k)-S(i,j,k-1)))
          drdkR = (drho_dT_u(I) * (T(i+1,j,k)-T(i+1,j,k-1)) + &
                   drho_dS_u(I) * (S(i+1,j,k)-S(i+1,j,k-1)))
          drdkDe_u(I,K) = (drdkR * e(i+1,j,K)) - (drdkL * e(i,j,K))
        elseif (find_work) then ! This is used in pure stacked SW mode
          drdkDe_u(I,K) = (drdkR * e(i+1,j,K)) - (drdkL * e(i,j,K))
        endif
        if (use_stanley) then
          ! Correction to the horizontal density gradient due to nonlinearity in
          ! the EOS rectifying SGS temperature anomalies
          drdiA = drdiA + 0.5 * ((drho_dT_dT_h(i+1) * tv%varT(i+1,j,k-1)) - &
                                (drho_dT_dT_h(i) * tv%varT(i,j,k-1)) )
          drdiB = drdiB + 0.5 * ((drho_dT_dT_h(i+1) * tv%varT(i+1,j,k)) - &
                                (drho_dT_dT_h(i) * tv%varT(i,j,k)) )
        endif
        if (find_work) drdi_u(I,k) = drdiB

        if (k > nk_linear) then
          if (use_EOS) then
            if (CS%use_FGNV_streamfn .or. find_work .or. .not.present_slope_x) then
              hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
              hg2R = h(i+1,j,k-1)*h(i+1,j,k) + h_neglect2
              haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
              haR = 0.5*(h(i+1,j,k-1) + h(i+1,j,k)) + h_neglect
              if (GV%Boussinesq) then
                dzaL = haL * GV%H_to_Z ; dzaR = haR * GV%H_to_Z
              elseif (GV%semi_Boussinesq) then
                dzaL = 0.5*(e(i,j,K-1) - e(i,j,K+1)) + dz_neglect
                dzaR = 0.5*(e(i+1,j,K-1) - e(i+1,j,K+1)) + dz_neglect
              else
                dzaL = 0.5*(dz(i,j,k-1) + dz(i,j,k)) + dz_neglect
                dzaR = 0.5*(dz(i+1,j,k-1) + dz(i+1,j,k)) + dz_neglect
              endif
              ! Use the harmonic mean thicknesses to weight the horizontal gradients.
              ! These unnormalized weights have been rearranged to minimize divisions.
              wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

              drdz = ((wtL * drdkL) + (wtR * drdkR)) / ((dzaL*wtL) + (dzaR*wtR))
              ! The expression for drdz above is mathematically equivalent to:
              !   drdz = ((hg2L/haL) * drdkL/dzaL + (hg2R/haR) * drdkR/dzaR) / &
              !          ((hg2L/haL) + (hg2R/haR))
              hg2A = h(i,j,k-1)*h(i+1,j,k-1) + h_neglect2
              hg2B = h(i,j,k)*h(i+1,j,k) + h_neglect2
              haA = 0.5*(h(i,j,k-1) + h(i+1,j,k-1)) + h_neglect
              haB = 0.5*(h(i,j,k) + h(i+1,j,k)) + h_neglect

              if (GV%Boussinesq) then
                N2_unlim = drdz*G_rho0
              else
                N2_unlim = (GV%g_Earth*GV%RZ_to_H) * &
                           (((wtL * drdkL) + (wtR * drdkR)) / ((haL*wtL) + (haR*wtR)))
              endif

              dzg2A = dz(i,j,k-1)*dz(i+1,j,k-1) + dz_neglect2
              dzg2B = dz(i,j,k)*dz(i+1,j,k) + dz_neglect2
              dzaA = 0.5*(dz(i,j,k-1) + dz(i+1,j,k-1)) + dz_neglect
              dzaB = 0.5*(dz(i,j,k) + dz(i+1,j,k)) + dz_neglect
              ! dzN2_u is used with the FGNV streamfunction formulation
              dzN2_u(I,K) = (0.5 * ( dzg2A / dzaA + dzg2B / dzaB )) * max(N2_unlim, N2_floor)
              if (find_work .and. CS%GM_src_alt) &
                hN2_x_PE(I,j,k) = (0.5 * ( hg2A / haA + hg2B / haB )) * max(N2_unlim, N2_floor)
            endif

            if (present_slope_x) then
              Slope = slope_x(I,j,k)
              slope2_Ratio_u(I,K) = Slope**2 * I_slope_max2
            else
              ! Use the harmonic mean thicknesses to weight the horizontal gradients.
              ! These unnormalized weights have been rearranged to minimize divisions.
              wtA = hg2A*haB ; wtB = hg2B*haA
              ! This is the gradient of density along geopotentials.
              drdx = ((wtA * drdiA + wtB * drdiB) / (wtA + wtB) - &
                      drdz * (e(i,j,K)-e(i+1,j,K))) * G%IdxCu(I,j)

              ! This estimate of slope is accurate for small slopes, but bounded
              ! to be between -1 and 1.
              mag_grad2 = (US%Z_to_L*drdx)**2 + drdz**2
              if (mag_grad2 > 0.0) then
                Slope = drdx / sqrt(mag_grad2)
                slope2_Ratio_u(I,K) = Slope**2 * I_slope_max2
              else ! Just in case mag_grad2 = 0 ever.
                Slope = 0.0
                slope2_Ratio_u(I,K) = 1.0e20  ! Force the use of the safe streamfunction.
              endif
            endif

            ! Adjust real slope by weights that bias towards slope of interfaces
            ! that ignore density gradients along layers.
            Slope = (1.0 - int_slope_u(I,j,K)) * Slope + &
                    int_slope_u(I,j,K) * ((e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j))
            slope2_Ratio_u(I,K) = (1.0 - int_slope_u(I,j,K)) * slope2_Ratio_u(I,K)

            if (CS%MEKE_src_slope_bug) then
              Slope_x_PE(I,j,k) = MIN(Slope, CS%slope_max)
            else
              Slope_x_PE(I,j,k) = Slope
              if (Slope > CS%slope_max) Slope_x_PE(I,j,k) = CS%slope_max
              if (Slope < -CS%slope_max) Slope_x_PE(I,j,k) = -CS%slope_max
            endif
            if (CS%id_slope_x > 0) CS%diagSlopeX(I,j,k) = Slope

            ! Estimate the streamfunction at each interface [H L2 T-1 ~> m3 s-1 or kg s-1].
            Sfn_unlim_u(I,K) = -(KH_u(I,j,K)*G%dy_Cu(I,j))*Slope

            ! Avoid moving dense water upslope from below the level of
            ! the bottom on the receiving side.
            if (Sfn_unlim_u(I,K) > 0.0) then ! The flow below this interface is positive.
              if (e(i,j,K) < e(i+1,j,nz+1)) then
                Sfn_unlim_u(I,K) = 0.0 ! This is not uhtot, because it may compensate for
                                ! deeper flow in very unusual cases.
              elseif (e(i+1,j,nz+1) > e(i,j,K+1)) then
                ! Scale the transport with the fraction of the donor layer above
                ! the bottom on the receiving side.
                Sfn_unlim_u(I,K) = Sfn_unlim_u(I,K) * ((e(i,j,K) - e(i+1,j,nz+1)) / &
                                         ((e(i,j,K) - e(i,j,K+1)) + dz_neglect))
              endif
            else
              if (e(i+1,j,K) < e(i,j,nz+1)) then ; Sfn_unlim_u(I,K) = 0.0
              elseif (e(i,j,nz+1) > e(i+1,j,K+1)) then
                Sfn_unlim_u(I,K) = Sfn_unlim_u(I,K) * ((e(i+1,j,K) - e(i,j,nz+1)) / &
                                       ((e(i+1,j,K) - e(i+1,j,K+1)) + dz_neglect))
              endif
            endif

          else ! .not. use_EOS
            if (present_slope_x) then
              Slope = slope_x(I,j,k)
            else
              Slope = ((e(i+1,j,K)-e(i,j,K))*G%IdxCu(I,j)) * G%OBCmaskCu(I,j)
            endif
            if (CS%id_slope_x > 0) CS%diagSlopeX(I,j,k) = Slope
            Sfn_unlim_u(I,K) = -(KH_u(I,j,K)*G%dy_Cu(I,j))*Slope
            dzN2_u(I,K) = GV%g_prime(K)
          endif ! if (use_EOS)
        else ! if (k > nk_linear)
          dzN2_u(I,K) = N2_floor * dz_neglect
          Sfn_unlim_u(I,K) = 0.
        endif ! if (k > nk_linear)
        if (CS%id_sfn_unlim_x>0) diag_sfn_unlim_x(I,j,K) = Sfn_unlim_u(I,K)
      enddo ! i-loop
    enddo ! k-loop

    if (CS%use_FGNV_streamfn) then
      do k=1,nz ; do I=is-1,ie ; if (G%OBCmaskCu(I,j)>0.) then
        dz_harm = max( dz_neglect, &
              2. * dz(i,j,k) * dz(i+1,j,k) / ( ( dz(i,j,k) + dz(i+1,j,k) ) + dz_neglect ) )
        c2_dz_u(I,k) = CS%FGNV_scale * ( 0.5*( cg1(i,j) + cg1(i+1,j) ) )**2 / dz_harm
      endif ; enddo ; enddo

      ! Solve an elliptic equation for the streamfunction following Ferrari et al., 2010.
      do I=is-1,ie
        if (G%OBCmaskCu(I,j)>0.) then
          do K=2,nz
            Sfn_unlim_u(I,K) = (1. + CS%FGNV_scale) * Sfn_unlim_u(I,K)
          enddo
          call streamfn_solver(nz, c2_dz_u(I,:), dzN2_u(I,:), Sfn_unlim_u(I,:))
        else
          do K=2,nz
            Sfn_unlim_u(I,K) = 0.
          enddo
        endif
      enddo
    endif

    do K=nz,2,-1
      do I=is-1,ie

        if (allocated(tv%SpV_avg) .and. (find_work .or. (k > nk_linear)) ) then
          Rho_avg = ( ((h(i,j,k) + h(i,j,k-1)) + (h(i+1,j,k) + h(i+1,j,k-1))) + 4.0*hn_2 ) / &
                ( (((h(i,j,k)+hn_2) * tv%SpV_avg(i,j,k))   + ((h(i,j,k-1)+hn_2) * tv%SpV_avg(i,j,k-1))) + &
                  (((h(i+1,j,k)+hn_2)*tv%SpV_avg(i+1,j,k)) + ((h(i+1,j,k-1)+hn_2)*tv%SpV_avg(i+1,j,k-1))) )
          ! Use an average density to convert the volume streamfunction estimate into a mass streamfunction.
          Z_to_H = GV%RZ_to_H*Rho_avg
        else
          Z_to_H = GV%Z_to_H
        endif

        if (k > nk_linear) then
          if (use_EOS) then

            if (uhtot(I,j) <= 0.0) then
              ! The transport that must balance the transport below is positive.
              Sfn_safe = uhtot(I,j) * (1.0 - h_frac(i,j,k))
            else !  (uhtot(I,j) > 0.0)
              Sfn_safe = uhtot(I,j) * (1.0 - h_frac(i+1,j,k))
            endif

            ! Determine the actual streamfunction at each interface.
            Sfn_est = (Z_to_H*Sfn_unlim_u(I,K) + slope2_Ratio_u(I,K)*Sfn_safe) / (1.0 + slope2_Ratio_u(I,K))
          else  ! When use_EOS is false, the layers are constant density.
            Sfn_est = Z_to_H*Sfn_unlim_u(I,K)
          endif

          ! Make sure that there is enough mass above to allow the streamfunction
          ! to satisfy the boundary condition of 0 at the surface.
          Sfn_in_H = min(max(Sfn_est, -h_avail_rsum(i,j,K)), h_avail_rsum(i+1,j,K))

          ! The actual transport is limited by the mass available in the two
          ! neighboring grid cells.
          uhD(I,j,k) = max(min((Sfn_in_H - uhtot(I,j)), h_avail(i,j,k)), &
                           -h_avail(i+1,j,k))

          if (CS%id_sfn_x>0) diag_sfn_x(I,j,K) = diag_sfn_x(I,j,K+1) + uhD(I,j,k)
!         sfn_x(I,j,K) = max(min(Sfn_in_h, uhtot(I,j)+h_avail(i,j,k)), &
!                            uhtot(I,j)-h_avail(i+1,j,K))
!         sfn_slope_x(I,j,K) = max(uhtot(I,j)-h_avail(i+1,j,k), &
!                                  min(uhtot(I,j)+h_avail(i,j,k), &
!               min(h_avail_rsum(i+1,j,K), max(-h_avail_rsum(i,j,K), &
!               (KH_u(I,j,K)*G%dy_Cu(I,j)) * ((e(i,j,K)-e(i+1,j,K))*G%IdxCu(I,j)) )) ))
        else ! k <= nk_linear
          ! Balance the deeper flow with a return flow uniformly distributed
          ! though the remaining near-surface layers.  This is the same as
          ! using Sfn_safe above.  There is no need to apply the limiters in
          ! this case.
          if (uhtot(I,j) <= 0.0) then
            uhD(I,j,k) = -uhtot(I,j) * h_frac(i,j,k)
          else !  (uhtot(I,j) > 0.0)
            uhD(I,j,k) = -uhtot(I,j) * h_frac(i+1,j,k)
          endif

          if (CS%id_sfn_x>0) diag_sfn_x(I,j,K) = diag_sfn_x(I,j,K+1) + uhD(I,j,k)
!         if (sfn_slope_x(I,j,K+1) <= 0.0) then
!           sfn_slope_x(I,j,K) = sfn_slope_x(I,j,K+1) * (1.0 - h_frac(i,j,k))
!         else
!           sfn_slope_x(I,j,K) = sfn_slope_x(I,j,K+1) * (1.0 - h_frac(i+1,j,k))
!         endif

        endif

        uhtot(I,j) = uhtot(I,j) + uhD(I,j,k)

        if (find_work) then
          !   This is the energy tendency based on the original profiles, and does
          ! not include any nonlinear terms due to a finite time step (which would
          ! involve interactions between the fluxes through the different faces.
          !   A second order centered estimate is used for the density transferred
          ! between water columns.

          if (allocated(tv%SpV_avg)) then
            G_scale = GV%H_to_RZ * GV%g_Earth / Rho_avg
          else
            G_scale = GV%g_Earth * GV%H_to_Z
          endif

          Work_u(I,j) = Work_u(I,j) + G_scale * &
            ( uhtot(I,j) * drdkDe_u(I,K) - &
              (uhD(I,j,k) * drdi_u(I,k)) * 0.25 * &
              ((e(i,j,K) + e(i,j,K+1)) + (e(i+1,j,K) + e(i+1,j,K+1))) )
        endif

      enddo
    enddo ! end of k-loop
  enddo ! end of j-loop

  ! Calculate the meridional fluxes and gradients.

  !$OMP parallel do default(none) shared(nz,is,ie,js,je,find_work,use_EOS,G,GV,US,pres,T,S,dz, &
  !$OMP                                  nk_linear,IsdB,tv,h,h_neglect,e,dz_neglect,dz_neglect2, &
  !$OMP                                  h_neglect2,int_slope_v,KH_v,vhtot,h_frac,h_avail_rsum, &
  !$OMP                                  I_slope_max2,vhD,h_avail,Work_v,CS,slope_y,cg1,hn_2,&
  !$OMP                                  diag_sfn_y,diag_sfn_unlim_y,N2_floor,EOSdom_v,use_stanley,&
  !$OMP                                  Tsgs2, present_slope_y,G_rho0,Slope_y_PE,hN2_y_PE)  &
  !$OMP                          private(drdjA,drdjB,drdkL,drdkR,pres_v,T_v,S_v,S_h,S_hr,    &
  !$OMP                                  drho_dT_v,drho_dS_v,hg2A,hg2B,hg2L,hg2R,haA,G_scale, &
  !$OMP                                  drho_dT_dT_h,drho_dT_dT_hr,scrap,pres_h,T_h,T_hr,   &
  !$OMP                                  haB,haL,haR,dzaL,dzaR,wtA,wtB,wtL,wtR,drdz,pres_hr, &
  !$OMP                                  dzg2A,dzg2B,dzaA,dzaB,dz_harm,Z_to_H, &
  !$OMP                                  drdy,mag_grad2,Slope,slope2_Ratio_v,dzN2_v,N2_unlim, &
  !$OMP                                  Sfn_unlim_v,Rho_avg,drdj_v,drdkDe_v,c2_dz_v, &
  !$OMP                                  Sfn_safe,Sfn_est,Sfn_in_h,calc_derivatives)
  do J=js-1,je
    do K=nz,2,-1
      if (find_work .and. .not.(use_EOS)) then
        drdjA = 0.0 ; drdjB = 0.0
        drdkL = GV%Rlay(k) - GV%Rlay(k-1) ; drdkR = drdkL
      endif

      calc_derivatives = use_EOS .and. (k >= nk_linear) .and. &
         (find_work .or. .not. present_slope_y .or. CS%use_FGNV_streamfn .or. use_stanley)

      if (calc_derivatives) then
        do i=is,ie
          pres_v(i) = 0.5*(pres(i,j,K) + pres(i,j+1,K))
          T_v(i) = 0.25*((T(i,j,k) + T(i,j+1,k)) + (T(i,j,k-1) + T(i,j+1,k-1)))
          S_v(i) = 0.25*((S(i,j,k) + S(i,j+1,k)) + (S(i,j,k-1) + S(i,j+1,k-1)))
        enddo
        call calculate_density_derivs(T_v, S_v, pres_v, drho_dT_v, drho_dS_v, &
                                      tv%eqn_of_state, EOSdom_v)
      endif
      if (use_stanley) then
        do i=is,ie
          pres_h(i) = pres(i,j,K)
          T_h(i) = 0.5*(T(i,j,k) + T(i,j,k-1))
          S_h(i) = 0.5*(S(i,j,k) + S(i,j,k-1))

          pres_hr(i) = pres(i,j+1,K)
          T_hr(i) = 0.5*(T(i,j+1,k) + T(i,j+1,k-1))
          S_hr(i) = 0.5*(S(i,j+1,k) + S(i,j+1,k-1))
        enddo

        ! The second line below would correspond to arguments
        !            drho_dS_dS, drho_dS_dT, drho_dT_dT, drho_dS_dP, drho_dT_dP, &
        call calculate_density_second_derivs(T_h, S_h, pres_h, &
                     scrap, scrap, drho_dT_dT_h, scrap, scrap, &
                     tv%eqn_of_state, EOSdom_v)
        call calculate_density_second_derivs(T_hr, S_hr, pres_hr, &
                     scrap, scrap, drho_dT_dT_hr, scrap, scrap, &
                     tv%eqn_of_state, EOSdom_v)
      endif
      do i=is,ie
        if (calc_derivatives) then
          ! Estimate the horizontal density gradients along layers.
          drdjA = drho_dT_v(i) * (T(i,j+1,k-1)-T(i,j,k-1)) + &
                  drho_dS_v(i) * (S(i,j+1,k-1)-S(i,j,k-1))
          drdjB = drho_dT_v(i) * (T(i,j+1,k)-T(i,j,k)) + &
                  drho_dS_v(i) * (S(i,j+1,k)-S(i,j,k))

          ! Estimate the vertical density gradients times the grid spacing.
          drdkL = (drho_dT_v(i) * (T(i,j,k)-T(i,j,k-1)) + &
                   drho_dS_v(i) * (S(i,j,k)-S(i,j,k-1)))
          drdkR = (drho_dT_v(i) * (T(i,j+1,k)-T(i,j+1,k-1)) + &
                   drho_dS_v(i) * (S(i,j+1,k)-S(i,j+1,k-1)))
          drdkDe_v(i,K) =  (drdkR * e(i,j+1,K)) - (drdkL * e(i,j,K))
        elseif (find_work) then ! This is used in pure stacked SW mode
          drdkDe_v(i,K) =  (drdkR * e(i,j+1,K)) - (drdkL * e(i,j,K))
        endif
        if (use_stanley) then
          ! Correction to the horizontal density gradient due to nonlinearity in
          ! the EOS rectifying SGS temperature anomalies
          drdjA = drdjA + 0.5 * ((drho_dT_dT_hr(i) * tv%varT(i,j+1,k-1)) - &
                                (drho_dT_dT_h(i) * tv%varT(i,j,k-1)) )
          drdjB = drdjB + 0.5 * ((drho_dT_dT_hr(i) * tv%varT(i,j+1,k)) - &
                                (drho_dT_dT_h(i) * tv%varT(i,j,k)) )
        endif

        if (find_work) drdj_v(i,k) = drdjB

        if (k > nk_linear) then
          if (use_EOS) then
            if (CS%use_FGNV_streamfn .or. find_work .or. .not. present_slope_y) then
              hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
              hg2R = h(i,j+1,k-1)*h(i,j+1,k) + h_neglect2
              haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
              haR = 0.5*(h(i,j+1,k-1) + h(i,j+1,k)) + h_neglect

              if (GV%Boussinesq) then
                dzaL = haL * GV%H_to_Z ; dzaR = haR * GV%H_to_Z
              elseif (GV%semi_Boussinesq) then
                dzaL = 0.5*(e(i,j,K-1) - e(i,j,K+1)) + dz_neglect
                dzaR = 0.5*(e(i,j+1,K-1) - e(i,j+1,K+1)) + dz_neglect
              else
                dzaL = 0.5*(dz(i,j,k-1) + dz(i,j,k)) + dz_neglect
                dzaR = 0.5*(dz(i,j+1,k-1) + dz(i,j+1,k)) + dz_neglect
              endif
              ! Use the harmonic mean thicknesses to weight the horizontal gradients.
              ! These unnormalized weights have been rearranged to minimize divisions.
              wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

              drdz = ((wtL * drdkL) + (wtR * drdkR)) / ((dzaL*wtL) + (dzaR*wtR))
              ! The expression for drdz above is mathematically equivalent to:
              !   drdz = ((hg2L/haL) * drdkL/dzaL + (hg2R/haR) * drdkR/dzaR) / &
              !          ((hg2L/haL) + (hg2R/haR))
              hg2A = h(i,j,k-1)*h(i,j+1,k-1) + h_neglect2
              hg2B = h(i,j,k)*h(i,j+1,k) + h_neglect2
              haA = 0.5*(h(i,j,k-1) + h(i,j+1,k-1)) + h_neglect
              haB = 0.5*(h(i,j,k) + h(i,j+1,k)) + h_neglect

              if (GV%Boussinesq) then
                N2_unlim = drdz*G_rho0
              else
                N2_unlim = (GV%g_Earth*GV%RZ_to_H) * &
                           (((wtL * drdkL) + (wtR * drdkR)) / ((haL*wtL) + (haR*wtR)))
              endif

              dzg2A = dz(i,j,k-1)*dz(i,j+1,k-1) + dz_neglect2
              dzg2B = dz(i,j,k)*dz(i,j+1,k) + dz_neglect2
              dzaA = 0.5*(dz(i,j,k-1) + dz(i,j+1,k-1)) + dz_neglect
              dzaB = 0.5*(dz(i,j,k) + dz(i,j+1,k)) + dz_neglect

              ! dzN2_v is used with the FGNV streamfunction formulation
              dzN2_v(i,K) = (0.5*( dzg2A / dzaA + dzg2B / dzaB )) * max(N2_unlim, N2_floor)
              if (find_work .and. CS%GM_src_alt) &
                hN2_y_PE(i,J,k) = (0.5*( hg2A / haA + hg2B / haB )) * max(N2_unlim, N2_floor)
            endif
            if (present_slope_y) then
              Slope = slope_y(i,J,k)
              slope2_Ratio_v(i,K) = Slope**2 * I_slope_max2
            else
              ! Use the harmonic mean thicknesses to weight the horizontal gradients.
              ! These unnormalized weights have been rearranged to minimize divisions.
              wtA = hg2A*haB ; wtB = hg2B*haA
              ! This is the gradient of density along geopotentials.
              drdy = ((wtA * drdjA + wtB * drdjB) / (wtA + wtB) - &
                      drdz * (e(i,j,K)-e(i,j+1,K))) * G%IdyCv(i,J)

              ! This estimate of slope is accurate for small slopes, but bounded
              ! to be between -1 and 1.
              mag_grad2 = (US%Z_to_L*drdy)**2 + drdz**2
              if (mag_grad2 > 0.0) then
                Slope = drdy / sqrt(mag_grad2)
                slope2_Ratio_v(i,K) = Slope**2 * I_slope_max2
              else ! Just in case mag_grad2 = 0 ever.
                Slope = 0.0
                slope2_Ratio_v(i,K) = 1.0e20  ! Force the use of the safe streamfunction.
              endif
            endif

            ! Adjust real slope by weights that bias towards slope of interfaces
            ! that ignore density gradients along layers.
            Slope = (1.0 - int_slope_v(i,J,K)) * Slope + &
                    int_slope_v(i,J,K) * ((e(i,j+1,K)-e(i,j,K)) * G%IdyCv(i,J))
            slope2_Ratio_v(i,K) = (1.0 - int_slope_v(i,J,K)) * slope2_Ratio_v(i,K)

            if (CS%MEKE_src_slope_bug) then
              Slope_y_PE(i,J,k) = MIN(Slope, CS%slope_max)
            else
              Slope_y_PE(i,J,k) = Slope
              if (Slope > CS%slope_max) Slope_y_PE(i,J,k) = CS%slope_max
              if (Slope < -CS%slope_max) Slope_y_PE(i,J,k) = -CS%slope_max
            endif
            if (CS%id_slope_y > 0) CS%diagSlopeY(I,j,k) = Slope

            Sfn_unlim_v(i,K) = -((KH_v(i,J,K)*G%dx_Cv(i,J))*Slope)

            ! Avoid moving dense water upslope from below the level of
            ! the bottom on the receiving side.
            if (Sfn_unlim_v(i,K) > 0.0) then ! The flow below this interface is positive.
              if (e(i,j,K) < e(i,j+1,nz+1)) then
                Sfn_unlim_v(i,K) = 0.0 ! This is not vhtot, because it may compensate for
                                ! deeper flow in very unusual cases.
              elseif (e(i,j+1,nz+1) > e(i,j,K+1)) then
                ! Scale the transport with the fraction of the donor layer above
                ! the bottom on the receiving side.
                Sfn_unlim_v(i,K) = Sfn_unlim_v(i,K) * ((e(i,j,K) - e(i,j+1,nz+1)) / &
                                         ((e(i,j,K) - e(i,j,K+1)) + dz_neglect))
              endif
            else
              if (e(i,j+1,K) < e(i,j,nz+1)) then ; Sfn_unlim_v(i,K) = 0.0
              elseif (e(i,j,nz+1) > e(i,j+1,K+1)) then
                Sfn_unlim_v(i,K) = Sfn_unlim_v(i,K) * ((e(i,j+1,K) - e(i,j,nz+1)) / &
                                       ((e(i,j+1,K) - e(i,j+1,K+1)) + dz_neglect))
              endif
            endif

          else ! .not. use_EOS
            if (present_slope_y) then
              Slope = slope_y(i,J,k)
            else
              Slope = ((e(i,j+1,K)-e(i,j,K))*G%IdyCv(i,J)) * G%OBCmaskCv(i,J)
            endif
            if (CS%id_slope_y > 0) CS%diagSlopeY(I,j,k) = Slope
            Sfn_unlim_v(i,K) = -((KH_v(i,J,K)*G%dx_Cv(i,J))*Slope)
            dzN2_v(i,K) = GV%g_prime(K)
          endif ! if (use_EOS)
        else ! if (k > nk_linear)
          dzN2_v(i,K) = N2_floor * dz_neglect
          Sfn_unlim_v(i,K) = 0.
        endif ! if (k > nk_linear)
        if (CS%id_sfn_unlim_y>0) diag_sfn_unlim_y(i,J,K) = Sfn_unlim_v(i,K)
      enddo ! i-loop
    enddo ! k-loop

    if (CS%use_FGNV_streamfn) then
      do k=1,nz ; do i=is,ie ; if (G%OBCmaskCv(i,J)>0.) then
        dz_harm = max( dz_neglect, &
              2. * dz(i,j,k) * dz(i,j+1,k) / ( ( dz(i,j,k) + dz(i,j+1,k) ) + dz_neglect ) )
        c2_dz_v(i,k) = CS%FGNV_scale * ( 0.5*( cg1(i,j) + cg1(i,j+1) ) )**2 / dz_harm
      endif ; enddo ; enddo

      ! Solve an elliptic equation for the streamfunction following Ferrari et al., 2010.
      do i=is,ie
        if (G%OBCmaskCv(i,J)>0.) then
          do K=2,nz
            Sfn_unlim_v(i,K) = (1. + CS%FGNV_scale) * Sfn_unlim_v(i,K)
          enddo
          call streamfn_solver(nz, c2_dz_v(i,:), dzN2_v(i,:), Sfn_unlim_v(i,:))
        else
          do K=2,nz
            Sfn_unlim_v(i,K) = 0.
          enddo
        endif
      enddo
    endif

    do K=nz,2,-1
      do i=is,ie
        if (allocated(tv%SpV_avg) .and. (find_work .or. (k > nk_linear)) ) then
          Rho_avg = ( ((h(i,j,k) + h(i,j,k-1)) + (h(i,j+1,k) + h(i,j+1,k-1))) + 4.0*hn_2 ) / &
              ( (((h(i,j,k)+hn_2) * tv%SpV_avg(i,j,k))   + ((h(i,j,k-1)+hn_2) * tv%SpV_avg(i,j,k-1))) + &
                (((h(i,j+1,k)+hn_2)*tv%SpV_avg(i,j+1,k)) + ((h(i,j+1,k-1)+hn_2)*tv%SpV_avg(i,j+1,k-1))) )
          ! Use an average density to convert the volume streamfunction estimate into a mass streamfunction.
          Z_to_H = GV%RZ_to_H*Rho_avg
        else
          Z_to_H = GV%Z_to_H
        endif

        if (k > nk_linear) then
          if (use_EOS) then

            if (vhtot(i,J) <= 0.0) then
              ! The transport that must balance the transport below is positive.
              Sfn_safe = vhtot(i,J) * (1.0 - h_frac(i,j,k))
            else !  (vhtot(I,j) > 0.0)
              Sfn_safe = vhtot(i,J) * (1.0 - h_frac(i,j+1,k))
            endif

            ! Find the actual streamfunction at each interface.
            Sfn_est = (Z_to_H*Sfn_unlim_v(i,K) + slope2_Ratio_v(i,K)*Sfn_safe) / (1.0 + slope2_Ratio_v(i,K))
          else  ! When use_EOS is false, the layers are constant density.
            Sfn_est = Z_to_H*Sfn_unlim_v(i,K)
          endif

          ! Make sure that there is enough mass above to allow the streamfunction
          ! to satisfy the boundary condition of 0 at the surface.
          Sfn_in_H = min(max(Sfn_est, -h_avail_rsum(i,j,K)), h_avail_rsum(i,j+1,K))

          ! The actual transport is limited by the mass available in the two
          ! neighboring grid cells.
          vhD(i,J,k) = max(min((Sfn_in_H - vhtot(i,J)), h_avail(i,j,k)), -h_avail(i,j+1,k))

          if (CS%id_sfn_y>0) diag_sfn_y(i,J,K) = diag_sfn_y(i,J,K+1) + vhD(i,J,k)
!         sfn_y(i,J,K) = max(min(Sfn_in_h, vhtot(i,J)+h_avail(i,j,k)), &
!                            vhtot(i,J)-h_avail(i,j+1,k))
!         sfn_slope_y(i,J,K) = max(vhtot(i,J)-h_avail(i,j+1,k), &
!                                  min(vhtot(i,J)+h_avail(i,j,k), &
!               min(h_avail_rsum(i,j+1,K), max(-h_avail_rsum(i,j,K), &
!               (KH_v(i,J,K)*G%dx_Cv(i,J)) * ((e(i,j,K)-e(i,j+1,K))*G%IdyCv(i,J)) )) ))
        else ! k <= nk_linear
          ! Balance the deeper flow with a return flow uniformly distributed
          ! though the remaining near-surface layers.  This is the same as
          ! using Sfn_safe above.  There is no need to apply the limiters in
          ! this case.
          if (vhtot(i,J) <= 0.0) then
            vhD(i,J,k) = -vhtot(i,J) * h_frac(i,j,k)
          else !  (vhtot(i,J) > 0.0)
            vhD(i,J,k) = -vhtot(i,J) * h_frac(i,j+1,k)
          endif

          if (CS%id_sfn_y>0) diag_sfn_y(i,J,K) = diag_sfn_y(i,J,K+1) + vhD(i,J,k)
!         if (sfn_slope_y(i,J,K+1) <= 0.0) then
!           sfn_slope_y(i,J,K) = sfn_slope_y(i,J,K+1) * (1.0 - h_frac(i,j,k))
!         else
!           sfn_slope_y(i,J,K) = sfn_slope_y(i,J,K+1) * (1.0 - h_frac(i,j+1,k))
!         endif
        endif

        vhtot(i,J) = vhtot(i,J)  + vhD(i,J,k)

        if (find_work) then
          !   This is the energy tendency based on the original profiles, and does
          ! not include any nonlinear terms due to a finite time step (which would
          ! involve interactions between the fluxes through the different faces.
          !   A second order centered estimate is used for the density transferred
          ! between water columns.

          if (allocated(tv%SpV_avg)) then
            G_scale = GV%H_to_RZ * GV%g_Earth / Rho_avg
          else
            G_scale = GV%g_Earth * GV%H_to_Z
          endif

          Work_v(i,J) = Work_v(i,J) + G_scale * &
            ( vhtot(i,J) * drdkDe_v(i,K) - &
             (vhD(i,J,k) * drdj_v(i,k)) * 0.25 * &
             ((e(i,j,K) + e(i,j,K+1)) + (e(i,j+1,K) + e(i,j+1,K+1))) )
        endif

      enddo
    enddo ! end of k-loop
  enddo ! end of j-loop

  ! In layer 1, enforce the boundary conditions that Sfn(z=0) = 0.0
  if (.not.find_work .or. .not.(use_EOS)) then
    do j=js,je ; do I=is-1,ie ; uhD(I,j,1) = -uhtot(I,j) ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; vhD(i,J,1) = -vhtot(i,J) ; enddo ; enddo
  else
    EOSdom_u(1) = (is-1) - (G%IsdB-1) ; EOSdom_u(2) = ie - (G%IsdB-1)
    !$OMP parallel do default(shared) private(pres_u,T_u,S_u,drho_dT_u,drho_dS_u,drdiB,G_scale)
    do j=js,je
      if (use_EOS) then
        do I=is-1,ie
          pres_u(I) = 0.5*(pres(i,j,1) + pres(i+1,j,1))
          T_u(I) = 0.5*(T(i,j,1) + T(i+1,j,1))
          S_u(I) = 0.5*(S(i,j,1) + S(i+1,j,1))
        enddo
        call calculate_density_derivs(T_u, S_u, pres_u, drho_dT_u, drho_dS_u, &
                                      tv%eqn_of_state, EOSdom_u )
      endif
      do I=is-1,ie
        uhD(I,j,1) = -uhtot(I,j)

        G_scale = GV%g_Earth * GV%H_to_Z
        if (use_EOS) then
          drdiB = drho_dT_u(I) * (T(i+1,j,1)-T(i,j,1)) + &
                  drho_dS_u(I) * (S(i+1,j,1)-S(i,j,1))
          if (allocated(tv%SpV_avg)) then
            G_scale = GV%H_to_RZ * GV%g_Earth * &
                ( ( ((h(i,j,1)+hn_2) * tv%SpV_avg(i,j,1)) + ((h(i+1,j,1)+hn_2) * tv%SpV_avg(i+1,j,1)) ) / &
                  ( (h(i,j,1) + h(i+1,j,1)) + 2.0*hn_2 ) )
          endif
        endif
        if (CS%use_GM_work_bug) then
          Work_u(I,j) = Work_u(I,j) + G_scale * &
              ( (uhD(I,j,1) * drdiB) * 0.25 * &
                ((e(i,j,1) + e(i,j,2)) + (e(i+1,j,1) + e(i+1,j,2))) )
        else
          Work_u(I,j) = Work_u(I,j) - G_scale * &
              ( (uhD(I,j,1) * drdiB) * 0.25 * &
                ((e(i,j,1) + e(i,j,2)) + (e(i+1,j,1) + e(i+1,j,2))) )
        endif
      enddo
    enddo

    EOSdom_v(:) = EOS_domain(G%HI)
    !$OMP parallel do default(shared) private(pres_v,T_v,S_v,drho_dT_v,drho_dS_v,drdjB,G_scale)
    do J=js-1,je
      if (use_EOS) then
        do i=is,ie
          pres_v(i) = 0.5*(pres(i,j,1) + pres(i,j+1,1))
          T_v(i) = 0.5*(T(i,j,1) + T(i,j+1,1))
          S_v(i) = 0.5*(S(i,j,1) + S(i,j+1,1))
        enddo
        call calculate_density_derivs(T_v, S_v, pres_v, drho_dT_v, drho_dS_v, &
                                      tv%eqn_of_state, EOSdom_v)
      endif
      do i=is,ie
        vhD(i,J,1) = -vhtot(i,J)

        G_scale = GV%g_Earth * GV%H_to_Z
        if (use_EOS) then
          drdjB = drho_dT_v(i) * (T(i,j+1,1)-T(i,j,1)) + &
                  drho_dS_v(i) * (S(i,j+1,1)-S(i,j,1))
          if (allocated(tv%SpV_avg)) then
            G_scale = GV%H_to_RZ * GV%g_Earth * &
                ( ( ((h(i,j,1)+hn_2) * tv%SpV_avg(i,j,1)) + ((h(i,j+1,1)+hn_2) * tv%SpV_avg(i,j+1,1)) ) / &
                  ( (h(i,j,1) + h(i,j+1,1)) + 2.0*hn_2 ) )
          endif
        endif
        Work_v(i,J) = Work_v(i,J) - G_scale * &
            ( (vhD(i,J,1) * drdjB) * 0.25 * &
              ((e(i,j,1) + e(i,j,2)) + (e(i,j+1,1) + e(i,j+1,2))) )
      enddo
    enddo
  endif

  if (find_work) then ; do j=js,je ; do i=is,ie
    ! Note that the units of Work_v and Work_u are [R Z L4 T-3 ~> W], while Work_h is in [R Z L2 T-3 ~> W m-2].
    Work_h = 0.5 * G%IareaT(i,j) * &
      ((Work_u(I-1,j) + Work_u(I,j)) + (Work_v(i,J-1) + Work_v(i,J)))
    if (allocated(CS%GMwork)) CS%GMwork(i,j) = Work_h
    if (.not. CS%GM_src_alt) then ; if (allocated(MEKE%GM_src)) then
      MEKE%GM_src(i,j) = MEKE%GM_src(i,j) + Work_h
    endif ; endif
    if (skeb_use_gm) then
      skeb_gm_work(i,j)   = STOCH%skeb_gm_coef * Work_h
      skeb_ebt_norm2(i,j) = 0.0
      do k=1,nz
        skeb_ebt_norm2(i,j) = skeb_ebt_norm2(i,j) + h(i,j,k) * VarMix%ebt_struct(i,j,k)**2
      enddo
      skeb_ebt_norm2(i,j) = GV%H_to_RZ * (skeb_ebt_norm2(i,j) + h_neglect)
    endif
  enddo ; enddo ; endif

  if (skeb_use_gm) then
    ! This block spreads the GM work down through the column using the ebt vertical structure, squared.
    ! Note the sign convention.
    do k=1,nz ; do j=js,je ; do i=is,ie
      STOCH%skeb_diss(i,j,k) = STOCH%skeb_diss(i,j,k) - skeb_gm_work(i,j) * &
                               VarMix%ebt_struct(i,j,k)**2 / skeb_ebt_norm2(i,j)
    enddo ; enddo ; enddo
  endif

  if (find_work .and. CS%GM_src_alt) then ; if (allocated(MEKE%GM_src)) then
    if (CS%MEKE_src_answer_date >= 20240601) then
      do j=js,je ; do i=is,ie ; do k=nz,1,-1
        PE_release_h = -0.25 * GV%H_to_RZ * &
                         ( ((KH_u(I,j,k)*(Slope_x_PE(I,j,k)**2) * hN2_x_PE(I,j,k)) + &
                            (Kh_u(I-1,j,k)*(Slope_x_PE(I-1,j,k)**2) * hN2_x_PE(I-1,j,k))) + &
                           ((Kh_v(i,J,k)*(Slope_y_PE(i,J,k)**2) * hN2_y_PE(i,J,k)) + &
                            (Kh_v(i,J-1,k)*(Slope_y_PE(i,J-1,k)**2) * hN2_y_PE(i,J-1,k))) )
        MEKE%GM_src(i,j) = MEKE%GM_src(i,j) + PE_release_h
      enddo ; enddo ; enddo
    else
      do j=js,je ; do i=is,ie ; do k=nz,1,-1
        PE_release_h = -0.25 * GV%H_to_RZ * &
                           ((KH_u(I,j,k)*(Slope_x_PE(I,j,k)**2) * hN2_x_PE(I,j,k)) + &
                            (Kh_u(I-1,j,k)*(Slope_x_PE(I-1,j,k)**2) * hN2_x_PE(I-1,j,k)) + &
                            (Kh_v(i,J,k)*(Slope_y_PE(i,J,k)**2) * hN2_y_PE(i,J,k)) + &
                            (Kh_v(i,J-1,k)*(Slope_y_PE(i,J-1,k)**2) * hN2_y_PE(i,J-1,k)))
        MEKE%GM_src(i,j) = MEKE%GM_src(i,j) + PE_release_h
      enddo ; enddo ; enddo
    endif

    if (CS%debug) then
      call hchksum(MEKE%GM_src, 'MEKE%GM_src', G%HI, unscale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
      call uvchksum("KH_[uv]", Kh_u, Kh_v, G%HI, unscale=US%L_to_m**2*US%s_to_T, &
                    scalar_pair=.true.)
      call uvchksum("Slope_[xy]_PE", Slope_x_PE, Slope_y_PE, G%HI, unscale=US%Z_to_L)
      call uvchksum("hN2_[xy]_PE", hN2_x_PE, hN2_y_PE, G%HI, unscale=GV%H_to_mks*US%L_to_Z**2*US%s_to_T**2, &
                    scalar_pair=.true.)
    endif
  endif ; endif

  if (CS%id_slope_x > 0) call post_data(CS%id_slope_x, CS%diagSlopeX, CS%diag)
  if (CS%id_slope_y > 0) call post_data(CS%id_slope_y, CS%diagSlopeY, CS%diag)
  if (CS%id_sfn_x > 0) call post_data(CS%id_sfn_x, diag_sfn_x, CS%diag)
  if (CS%id_sfn_y > 0) call post_data(CS%id_sfn_y, diag_sfn_y, CS%diag)
  if (CS%id_sfn_unlim_x > 0) call post_data(CS%id_sfn_unlim_x, diag_sfn_unlim_x, CS%diag)
  if (CS%id_sfn_unlim_y > 0) call post_data(CS%id_sfn_unlim_y, diag_sfn_unlim_y, CS%diag)

end subroutine thickness_diffuse_full

!> Tridiagonal solver for streamfunction at interfaces
subroutine streamfn_solver(nk, c2_h, hN2, sfn)
  integer,               intent(in)    :: nk   !< Number of layers
  real, dimension(nk),   intent(in)    :: c2_h !< Wave speed squared over thickness in layers, rescaled to
                                               !! [H L2 Z-2 T-2 ~> m s-2 or kg m-2 s-2]
  real, dimension(nk+1), intent(in)    :: hN2  !< Thickness times N2 at interfaces times rescaling factors
                                               !! [H L2 Z-2 T-2 ~> m s-2 or kg m-2 s-2]
  real, dimension(nk+1), intent(inout) :: sfn  !< Streamfunction [H L2 T-1 ~> m3 s-1 or kg s-1] or arbitrary units
                                               !! On entry, equals diffusivity times slope.
                                               !! On exit, equals the streamfunction.
  ! Local variables
  real :: c1(nk)  ! The dependence of the final streamfunction on the values below [nondim]
  real :: d1      ! The complement of c1(k) (i.e., 1 - c1(k)) [nondim]
  real :: b_denom ! A term in the denominator of beta [H L2 Z-2 T-2 ~> m s-2 or kg m-2 s-2]
  real :: beta    ! The normalization for the pivot [Z2 T2 H-1 L-2 ~> s2 m-1 or m2 s2 kg-1]
  integer :: k

  sfn(1) = 0.
  b_denom = hN2(2) + c2_h(1)
  beta = 1.0 / ( b_denom + c2_h(2) )
  d1 = beta * b_denom
  sfn(2) = ( beta * hN2(2) )*sfn(2)
  do K=3,nk
    c1(k-1) = beta * c2_h(k-1)
    b_denom = hN2(K) + d1*c2_h(k-1)
    beta = 1.0 / (b_denom + c2_h(k))
    d1 = beta * b_denom
    sfn(K) = beta * (hN2(K)*sfn(K) + c2_h(k-1)*sfn(K-1))
  enddo
  c1(nk) = beta * c2_h(nk)
  sfn(nk+1) = 0.
  do K=nk,2,-1
    sfn(K) = sfn(K) + c1(k)*sfn(K+1)
  enddo

end subroutine streamfn_solver

!> Add a diffusivity that acts on the isopycnal heights, regardless of the densities
subroutine add_interface_Kh(G, GV, US, CS, Kh_u, Kh_v, Kh_u_CFL, Kh_v_CFL, int_slope_u, int_slope_v)
  type(ocean_grid_type),                        intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV   !< Vertical grid structure
  type(unit_scale_type),                        intent(in)    :: US   !< A dimensional unit scaling type
  type(thickness_diffuse_CS),                   intent(in)    :: CS   !< Control structure for thickness_diffuse
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kh_u !< Isopycnal height diffusivity
                                                                      !! at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: Kh_v !< Isopycnal height diffusivity
                                                                      !! at v points [L2 T-1 ~> m2 s-1]
  real, dimension(SZIB_(G),SZJ_(G)),            intent(in)    :: Kh_u_CFL !< Maximum stable isopycnal height
                                                                      !! diffusivity at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G)),            intent(in)    :: Kh_v_CFL !< Maximum stable isopycnal height
                                                                      !! diffusivity at v points [L2 T-1 ~> m2 s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: int_slope_u !< Ratio that determine how much of
                                                                      !! the isopycnal slopes are taken directly from
                                                                      !! the interface slopes without consideration
                                                                      !! of density gradients [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: int_slope_v !< Ratio that determine how much of
                                                                      !! the isopycnal slopes are taken directly from
                                                                      !! the interface slopes without consideration
                                                                      !! of density gradients [nondim].

  ! Local variables
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  do k=1,nz+1 ; do j=js,je ; do I=is-1,ie ; if (CS%Kh_eta_u(I,j) > 0.0) then
    int_slope_u(I,j,K) = (int_slope_u(I,j,K)*Kh_u(I,j,K) + CS%Kh_eta_u(I,j)) / &
                         (Kh_u(I,j,K) + CS%Kh_eta_u(I,j))
    Kh_u(I,j,K) = min(Kh_u(I,j,K) + CS%Kh_eta_u(I,j), Kh_u_CFL(I,j))
  endif ; enddo ; enddo ; enddo

  do k=1,nz+1 ; do J=js-1,je ; do i=is,ie ; if (CS%Kh_eta_v(i,J) > 0.0) then
    int_slope_v(i,J,K) = (int_slope_v(i,J,K)*Kh_v(i,J,K) + CS%Kh_eta_v(i,J)) / &
                         (Kh_v(i,J,K) + CS%Kh_eta_v(i,J))
    Kh_v(i,J,K) = min(Kh_v(i,J,K) + CS%Kh_eta_v(i,J), Kh_v_CFL(i,J))
  endif ; enddo ; enddo ; enddo

end subroutine add_interface_Kh

!> Modifies isopycnal height diffusivities to untangle layer structures
subroutine add_detangling_Kh(h, e, Kh_u, Kh_v, KH_u_CFL, KH_v_CFL, tv, dt, G, GV, US, CS, &
                             int_slope_u, int_slope_v)
  type(ocean_grid_type),                        intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV   !< Vertical grid structure
  type(unit_scale_type),                        intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)    :: h    !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(in)    :: e    !< Interface positions [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kh_u !< Isopycnal height diffusivity
                                                                      !! at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: Kh_v !< Isopycnal height diffusivity
                                                                      !! at v points [L2 T-1 ~> m2 s-1]
  real, dimension(SZIB_(G),SZJ_(G)),            intent(in)    :: Kh_u_CFL !< Maximum stable isopycnal height
                                                                      !! diffusivity at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G)),            intent(in)    :: Kh_v_CFL !< Maximum stable isopycnal height
                                                                      !! diffusivity at v points [L2 T-1 ~> m2 s-1]
  type(thermo_var_ptrs),                        intent(in)    :: tv   !< Thermodynamics structure
  real,                                         intent(in)    :: dt   !< Time increment [T ~> s]
  type(thickness_diffuse_CS),                   intent(in)    :: CS   !< Control structure for thickness_diffuse
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: int_slope_u !< Ratio that determine how much of
                                                                      !! the isopycnal slopes are taken directly from
                                                                      !! the interface slopes without consideration
                                                                      !! of density gradients [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: int_slope_v !< Ratio that determine how much of
                                                                      !! the isopycnal slopes are taken directly from
                                                                      !! the interface slopes without consideration
                                                                      !! of density gradients [nondim].
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    de_top     ! The distances between the top of a layer and the top of the
               ! region where the detangling is applied [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    Kh_lay_u   ! The tentative isopycnal height diffusivity for each layer at
               ! u points [L2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    Kh_lay_v   ! The tentative isopycnal height diffusivity for each layer at
               ! v points [L2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    de_bot     ! The distances from the bottom of the region where the
               ! detangling is applied [H ~> m or kg m-2].
  real :: h1, h2    ! The thinner and thicker surrounding thicknesses [H ~> m or kg m-2],
                    ! with the thinner modified near the boundaries to mask out
                    ! thickness variations due to topography, etc.
  real :: jag_Rat   ! The nondimensional jaggedness ratio for a layer, going
                    ! from 0 (smooth) to 1 (jagged) [nondim].  This is the difference
                    ! between the arithmetic and harmonic mean thicknesses
                    ! normalized by the arithmetic mean thickness.
  real :: Kh_scale  ! A ratio by which Kh_u_CFL is scaled for maximally jagged
                    ! layers [nondim].
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].

  real :: I_sl      ! The absolute value of the larger in magnitude of the slopes
                    ! above and below [L Z-1 ~> nondim].
  real :: Rsl       ! The ratio of the smaller magnitude slope to the larger
                    ! magnitude one [nondim]. 0 <= Rsl <1.
  real :: IRsl      ! The (limited) inverse of Rsl [nondim]. 1 < IRsl <= 1e9.
  real :: dH        ! The thickness gradient divided by the damping timescale
                    ! and the ratio of the face length to the adjacent cell
                    ! areas for comparability with the diffusivities [L Z T-1 ~> m2 s-1].
  real :: adH       ! The absolute value of dH [L Z T-1 ~> m2 s-1].
  real :: sign      ! 1 or -1, with the same sign as the layer thickness gradient [nondim].
  real :: sl_K      ! The sign-corrected slope of the interface above [Z L-1 ~> nondim].
  real :: sl_Kp1    ! The sign-corrected slope of the interface below [Z L-1 ~> nondim].
  real :: I_sl_K    ! The (limited) inverse of sl_K [L Z-1 ~> nondim].
  real :: I_sl_Kp1  ! The (limited) inverse of sl_Kp1 [L Z-1 ~> nondim].
  real :: I_4t      ! A quarter of a flux scaling factor divided by
                    ! the damping timescale [T-1 ~> s-1].
  real :: Fn_R      ! A function of Rsl, such that Rsl < Fn_R < 1 [nondim]
  real :: Idx_eff   ! The effective inverse x-grid spacing at a u-point [L-1 ~> m-1]
  real :: Idy_eff   ! The effective inverse y-grid spacing at a v-point [L-1 ~> m-1]
  real :: slope_sq  ! The sum of the squared slopes above and below a layer [Z2 L-2 ~> nondim]
  real :: Kh_max    ! A local ceiling on the diffusivity [L2 T-1 ~> m2 s-1].
  real :: wt1, wt2  ! Nondimensional weights [nondim].
  !   Variables used only in testing code.
  ! real, dimension(SZK_(GV)) :: uh_here ! The transport in a layer [Z L2 T-1 ~> m3 s-1]
  ! real, dimension(SZK_(GV)+1) :: Sfn ! The streamfunction at an interface [Z L T-1 ~> m2 s-1]
  real :: dKh       ! An increment in the diffusivity [L2 T-1 ~> m2 s-1].

  real, dimension(SZIB_(G),SZK_(GV)+1) :: &
    Kh_bg, &        ! The background (floor) value of Kh [L2 T-1 ~> m2 s-1].
    Kh, &           ! The tentative value of Kh [L2 T-1 ~> m2 s-1].
    Kh_detangle, &  ! The detangling diffusivity that could be used [L2 T-1 ~> m2 s-1].
    Kh_min_max_p, & ! The smallest ceiling that can be placed on Kh(I,K)
                    ! based on the value of Kh(I,K+1) [L2 T-1 ~> m2 s-1].
    Kh_min_max_m, & ! The smallest ceiling that can be placed on Kh(I,K)
                    ! based on the value of Kh(I,K-1) [L2 T-1 ~> m2 s-1].
    ! The following are variables that define the relationships between
    ! successive values of Kh.
    ! Search for Kh that satisfy...
    !    Kh(I,K) >= Kh_min_m(I,K)*Kh(I,K-1) + Kh0_min_m(I,K)
    !    Kh(I,K) >= Kh_min_p(I,K)*Kh(I,K+1) + Kh0_min_p(I,K)
    !    Kh(I,K) <= Kh_max_m(I,K)*Kh(I,K-1) + Kh0_max_m(I,K)
    !    Kh(I,K) <= Kh_max_p(I,K)*Kh(I,K+1) + Kh0_max_p(I,K)
    Kh_min_m , &   ! See above [nondim].
    Kh0_min_m , &  ! See above [L2 T-1 ~> m2 s-1].
    Kh_max_m , &   ! See above [nondim].
    Kh0_max_m, &   ! See above [L2 T-1 ~> m2 s-1].
    Kh_min_p , &   ! See above [nondim].
    Kh0_min_p , &  ! See above [L2 T-1 ~> m2 s-1].
    Kh_max_p , &   ! See above [nondim].
    Kh0_max_p      ! See above [L2 T-1 ~> m2 s-1].
  real, dimension(SZIB_(G)) :: &
    Kh_max_max  ! The maximum diffusivity permitted in a column [L2 T-1 ~> m2 s-1]
  logical, dimension(SZIB_(G)) :: &
    do_i        ! If true, work on a column.
  integer :: i, j, k, n, ish, jsh, is, ie, js, je, nz, k_top
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  k_top = GV%nk_rho_varies + 1
  h_neglect = GV%H_subroundoff
  !   The 0.5 is because we are not using uniform weightings, but are
  ! distributing the diffusivities more effectively (with wt1 & wt2), but this
  ! means that the additions to a single interface can be up to twice as large.
  Kh_scale = 0.5
  if (CS%detangle_time > dt) Kh_scale = 0.5 * dt / CS%detangle_time

  do j=js-1,je+1 ; do i=is-1,ie+1
    de_top(i,j,k_top) = 0.0 ; de_bot(i,j) = 0.0
  enddo ; enddo
  do k=k_top+1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
    de_top(i,j,k) = de_top(i,j,k-1) + h(i,j,k-1)
  enddo ; enddo ; enddo

  do j=js,je ; do I=is-1,ie
    Kh_lay_u(I,j,nz) = 0.0 ; Kh_lay_u(I,j,k_top) = 0.0
  enddo ; enddo
  do J=js-1,je ; do i=is,ie
    Kh_lay_v(i,J,nz) = 0.0 ; Kh_lay_v(i,J,k_top) = 0.0
  enddo ; enddo

  do k=nz-1,k_top+1,-1
    ! Find the diffusivities associated with each layer.
    do j=js-1,je+1 ; do i=is-1,ie+1
      de_bot(i,j) = de_bot(i,j) + h(i,j,k+1)
    enddo ; enddo

    do j=js,je ; do I=is-1,ie ; if (G%OBCmaskCu(I,j) > 0.0) then
      if (h(i,j,k) > h(i+1,j,k)) then
        h2 = h(i,j,k)
        h1 = max( h(i+1,j,k), h2 - min(de_bot(i+1,j), de_top(i+1,j,k)) )
      else
        h2 = h(i+1,j,k)
        h1 = max( h(i,j,k), h2 - min(de_bot(i,j), de_top(i,j,k)) )
      endif
      jag_Rat = (h2 - h1)**2 / (h2 + h1 + h_neglect)**2
      KH_lay_u(I,j,k) = (Kh_scale * KH_u_CFL(I,j)) * jag_Rat**2
    endif ; enddo ; enddo

    do J=js-1,je ; do i=is,ie ; if (G%OBCmaskCv(i,J) > 0.0) then
      if (h(i,j,k) > h(i,j+1,k)) then
        h2 = h(i,j,k)
        h1 = max( h(i,j+1,k), h2 - min(de_bot(i,j+1), de_top(i,j+1,k)) )
      else
        h2 = h(i,j+1,k)
        h1 = max( h(i,j,k), h2 - min(de_bot(i,j), de_top(i,j,k)) )
      endif
      jag_Rat = (h2 - h1)**2 / (h2 + h1 + h_neglect)**2
      KH_lay_v(i,J,k) = (Kh_scale * KH_v_CFL(i,J)) * jag_Rat**2
    endif ; enddo ; enddo
  enddo

  ! Limit the diffusivities

  I_4t = Kh_scale / (4.0 * dt)

  do n=1,2
    if (n==1) then ; jsh = js ; ish = is-1
    else ; jsh = js-1 ; ish = is ; endif

    do j=jsh,je

      ! First, populate the diffusivities
      if (n==1) then ! This is a u-column.
        do i=ish,ie
          do_i(I) = (G%OBCmaskCu(I,j) > 0.0)
          Kh_Max_max(I) = KH_u_CFL(I,j)
        enddo
        do K=1,nz+1 ; do i=ish,ie
          Kh_bg(I,K) = KH_u(I,j,K) ; Kh(I,K) = Kh_bg(I,K)
          Kh_min_max_p(I,K) = Kh_bg(I,K) ; Kh_min_max_m(I,K) = Kh_bg(I,K)
          Kh_detangle(I,K) = 0.0
        enddo ; enddo
      else ! This is a v-column.
        do i=ish,ie
          do_i(i) = (G%OBCmaskCv(i,J) > 0.0) ; Kh_Max_max(I) = KH_v_CFL(i,J)
        enddo
        do K=1,nz+1 ; do i=ish,ie
          Kh_bg(I,K) = KH_v(I,j,K) ; Kh(I,K) = Kh_bg(I,K)
          Kh_min_max_p(I,K) = Kh_bg(I,K) ; Kh_min_max_m(I,K) = Kh_bg(I,K)
          Kh_detangle(I,K) = 0.0
        enddo ; enddo
      endif

      ! Determine the limits on the diffusivities.
      do k=k_top,nz ; do i=ish,ie ; if (do_i(i)) then
        if (n==1) then ! This is a u-column.
          dH = 0.0
          Idx_eff = ((G%IareaT(i+1,j) + G%IareaT(i,j)) * G%dy_Cu(I,j))
          !   This expression uses differences in e in place of h for better
          ! consistency with the slopes.
          if (Idx_eff > 0.0) &
            dH = I_4t * ((e(i+1,j,K) - e(i+1,j,K+1)) - &
                         (e(i,j,K) - e(i,j,K+1))) / Idx_eff
           ! dH = I_4t * (h(i+1,j,k) - h(i,j,k)) / Idx_eff

          adH = abs(dH)
          sign = 1.0 ; if (dH < 0) sign = -1.0
          sl_K = sign * (e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j)
          sl_Kp1 = sign * (e(i+1,j,K+1)-e(i,j,K+1)) * G%IdxCu(I,j)

          ! Add the incremental diffusivities to the surrounding interfaces.
          ! Adding more to the more steeply sloping layers (as below) makes
          ! the diffusivities more than twice as effective.
          slope_sq = (sl_K**2 + sl_Kp1**2)
          wt1 = 0.5 ; wt2 = 0.5
          if (slope_sq > 0.0) then
            wt1 = sl_K**2 / slope_sq ; wt2 = sl_Kp1**2 / slope_sq
          endif
          Kh_detangle(I,K) = Kh_detangle(I,K) + wt1*KH_lay_u(I,j,k)
          Kh_detangle(I,K+1) = Kh_detangle(I,K+1) + wt2*KH_lay_u(I,j,k)
        else ! This is a v-column.
          dH = 0.0
          Idy_eff = ((G%IareaT(i,j+1) + G%IareaT(i,j)) * G%dx_Cv(I,j))
          if (Idy_eff > 0.0) &
            dH = I_4t * ((e(i,j+1,K) - e(i,j+1,K+1)) - &
                         (e(i,j,K) - e(i,j,K+1))) / Idy_eff
           ! dH = I_4t * (h(i,j+1,k) - h(i,j,k)) / Idy_eff

          adH = abs(dH)
          sign = 1.0 ; if (dH < 0) sign = -1.0
          sl_K = sign * (e(i,j+1,K)-e(i,j,K)) * G%IdyCv(i,J)
          sl_Kp1 = sign * (e(i,j+1,K+1)-e(i,j,K+1)) * G%IdyCv(i,J)

          ! Add the incremental diffusivities to the surrounding interfaces.
          ! Adding more to the more steeply sloping layers (as below) makes
          ! the diffusivities more than twice as effective.
          slope_sq = (sl_K**2 + sl_Kp1**2)
          wt1 = 0.5 ; wt2 = 0.5
          if (slope_sq > 0.0) then
            wt1 = sl_K**2 / slope_sq ; wt2 = sl_Kp1**2 / slope_sq
          endif
          Kh_detangle(I,K) = Kh_detangle(I,K) + wt1*KH_lay_v(i,J,k)
          Kh_detangle(I,K+1) = Kh_detangle(I,K+1) + wt2*KH_lay_v(i,J,k)
        endif

        if (adH == 0.0) then
          Kh_min_m(I,K+1) = 1.0 ; Kh0_min_m(I,K+1) = 0.0
          Kh_max_m(I,K+1) = 1.0 ; Kh0_max_m(I,K+1) = 0.0
          Kh_min_p(I,K) = 1.0 ; Kh0_min_p(I,K) = 0.0
          Kh_max_p(I,K) = 1.0 ; Kh0_max_p(I,K) = 0.0
        elseif (adH > 0.0) then
          if (sl_K <= sl_Kp1) then
            ! This case should only arise from nonlinearities in the equation of state.
            ! Treat it as though dedx(K) = dedx(K+1) & dH = 0.
            Kh_min_m(I,K+1) = 1.0 ; Kh0_min_m(I,K+1) = 0.0
            Kh_max_m(I,K+1) = 1.0 ; Kh0_max_m(I,K+1) = 0.0
            Kh_min_p(I,K) = 1.0 ; Kh0_min_p(I,K) = 0.0
            Kh_max_p(I,K) = 1.0 ; Kh0_max_p(I,K) = 0.0
          elseif (sl_K <= 0.0) then   ! Both slopes are opposite to dH
            I_sl = -1.0 / sl_Kp1
            Rsl = -sl_K * I_sl                            ! 0 <= Rsl < 1
            IRsl = 1e9 ; if (Rsl > 1e-9) IRsl = 1.0/Rsl   ! 1 < IRsl <= 1e9

            Fn_R = Rsl
            if (Kh_max_max(I) > 0) &
              Fn_R = min(sqrt(Rsl), Rsl + (adH * I_sl) / (Kh_Max_max(I)))

            Kh_min_m(I,K+1) = Fn_R ; Kh0_min_m(I,K+1) = 0.0
            Kh_max_m(I,K+1) = Rsl ; Kh0_max_m(I,K+1) = adH * I_sl
            Kh_min_p(I,K) = IRsl ; Kh0_min_p(I,K) = -adH * (I_sl*IRsl)
            Kh_max_p(I,K) = 1.0/(Fn_R + 1.0e-30) ; Kh0_max_p(I,K) = 0.0
          elseif (sl_Kp1 < 0.0) then  ! Opposite (nonzero) signs of slopes.
            I_sl_K = 1e18*US%Z_to_L ; if (sl_K > 1e-18*US%L_to_Z) I_sl_K = 1.0 / sl_K
            I_sl_Kp1 = 1e18*US%Z_to_L ; if (-sl_Kp1 > 1e-18*US%L_to_Z) I_sl_Kp1 = -1.0 / sl_Kp1

            Kh_min_m(I,K+1) = 0.0 ; Kh0_min_m(I,K+1) = 0.0
            Kh_max_m(I,K+1) = - sl_K*I_sl_Kp1 ; Kh0_max_m(I,K+1) = adH*I_sl_Kp1
            Kh_min_p(I,K) = 0.0 ; Kh0_min_p(I,K) = 0.0
            Kh_max_p(I,K) = sl_Kp1*I_sl_K ; Kh0_max_p(I,K) = adH*I_sl_K

            ! This limit does not use the slope weighting so that potentially
            ! sharp gradients in diffusivities are not forced to occur.
            Kh_Max = adH / (sl_K - sl_Kp1)
            Kh_min_max_p(I,K) = max(Kh_min_max_p(I,K), Kh_Max)
            Kh_min_max_m(I,K+1) = max(Kh_min_max_m(I,K+1), Kh_Max)
          else ! Both slopes are of the same sign as dH.
            I_sl = 1.0 / sl_K
            Rsl = sl_Kp1 * I_sl                           ! 0 <= Rsl < 1
            IRsl = 1e9 ; if (Rsl > 1e-9) IRsl = 1.0/Rsl   ! 1 < IRsl <= 1e9

            ! Rsl <= Fn_R <= 1
            Fn_R = Rsl
            if (Kh_max_max(I) > 0) &
              Fn_R = min(sqrt(Rsl), Rsl + (adH * I_sl) / Kh_Max_max(I))

            Kh_min_m(I,K+1) = IRsl ; Kh0_min_m(I,K+1) = -adH * (I_sl*IRsl)
            Kh_max_m(I,K+1) = 1.0/(Fn_R + 1.0e-30) ; Kh0_max_m(I,K+1) = 0.0
            Kh_min_p(I,K) = Fn_R ; Kh0_min_p(I,K) = 0.0
            Kh_max_p(I,K) = Rsl ; Kh0_max_p(I,K) = adH * I_sl
          endif
        endif
      endif ; enddo ; enddo ! I-loop & k-loop

      do k=k_top,nz+1,nz+1-k_top ; do i=ish,ie ; if (do_i(i)) then
        ! The diffusivities at k_top and nz+1 are both fixed.
        Kh_min_m(I,k) = 0.0 ; Kh0_min_m(I,k) = 0.0
        Kh_max_m(I,k) = 0.0 ; Kh0_max_m(I,k) = 0.0
        Kh_min_p(I,k) = 0.0 ; Kh0_min_p(I,k) = 0.0
        Kh_max_p(I,k) = 0.0 ; Kh0_max_p(I,k) = 0.0
        Kh_min_max_p(I,K) = Kh_bg(I,K)
        Kh_min_max_m(I,K) = Kh_bg(I,K)
      endif ; enddo ; enddo ! I-loop and k_top/nz+1 loop

      ! Search for Kh that satisfy...
      !    Kh(I,K) >= Kh_min_m(I,K)*Kh(I,K-1) + Kh0_min_m(I,K)
      !    Kh(I,K) >= Kh_min_p(I,K)*Kh(I,K+1) + Kh0_min_p(I,K)
      !    Kh(I,K) <= Kh_max_m(I,K)*Kh(I,K-1) + Kh0_max_m(I,K)
      !    Kh(I,K) <= Kh_max_p(I,K)*Kh(I,K+1) + Kh0_max_p(I,K)

      ! Increase the diffusivities to satisfy the min constraints.
      ! All non-zero min constraints on one diffusivity are max constraints on another.
      do K=k_top+1,nz ; do i=ish,ie ; if (do_i(i)) then
        Kh(I,K) = max(Kh_bg(I,K), Kh_detangle(I,K), &
                      min(Kh_min_m(I,K)*Kh(I,K-1) + Kh0_min_m(I,K), Kh(I,K-1)))

        if (Kh0_max_m(I,K) > Kh_bg(I,K)) Kh(I,K) = min(Kh(I,K), Kh0_max_m(I,K))
        if (Kh0_max_p(I,K) > Kh_bg(I,K)) Kh(I,K) = min(Kh(I,K), Kh0_max_p(I,K))
      endif ; enddo ; enddo ! I-loop & k-loop
      ! This is still true... do i=ish,ie ; Kh(I,nz+1) = Kh_bg(I,nz+1) ; enddo
      do K=nz,k_top+1,-1 ; do i=ish,ie ; if (do_i(i)) then
        Kh(I,k) = max(Kh(I,K), min(Kh_min_p(I,K)*Kh(I,K+1) + Kh0_min_p(I,K), Kh(I,K+1)))

        Kh_Max = max(Kh_min_max_p(I,K), Kh_max_p(I,K)*Kh(I,K+1) + Kh0_max_p(I,K))
        Kh(I,k) = min(Kh(I,k), Kh_Max)
      endif ; enddo ; enddo ! I-loop & k-loop
      !  All non-zero min constraints on one diffusivity are max constraints on
      ! another layer, so the min constraints can now be discounted.

      ! Decrease the diffusivities to satisfy the max constraints.
        do K=k_top+1,nz ; do i=ish,ie ; if (do_i(i)) then
          Kh_Max = max(Kh_min_max_m(I,K), Kh_max_m(I,K)*Kh(I,K-1) + Kh0_max_m(I,K))
          if (Kh(I,k) > Kh_Max) Kh(I,k) = Kh_Max
        endif ; enddo ; enddo  ! i- and K-loops

      ! This code tests the solutions...
!     do i=ish,ie
!       Sfn(:) = 0.0 ; uh_here(:) = 0.0
!       do K=k_top,nz
!         if ((Kh(i,K) > Kh_bg(i,K)) .or. (Kh(i,K+1) > Kh_bg(i,K+1))) then
!           if (n==1) then ! u-point.
!             if ((h(i+1,j,k) - h(i,j,k)) * &
!                 ((e(i+1,j,K)-e(i+1,j,K+1)) - (e(i,j,K)-e(i,j,K+1))) > 0.0) then
!               Sfn(K) = -Kh(i,K) * (e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j)
!               Sfn(K+1) = -Kh(i,K+1) * (e(i+1,j,K+1)-e(i,j,K+1)) * G%IdxCu(I,j)
!               uh_here(k) = (Sfn(K) - Sfn(K+1))*G%dy_Cu(I,j)
!               if (abs(uh_here(k)) * min(G%IareaT(i,j), G%IareaT(i+1,j)) > &
!                   (1e-10*GV%m_to_H)) then
!                 if (uh_here(k) * (h(i+1,j,k) - h(i,j,k)) > 0.0) then
!                   call MOM_error(WARNING, "Corrective u-transport is up the thickness gradient.", .true.)
!                 endif
!                 if (((h(i,j,k) - 4.0*dt*G%IareaT(i,j)*uh_here(k)) - &
!                      (h(i+1,j,k) + 4.0*dt*G%IareaT(i+1,j)*uh_here(k))) * &
!                     (h(i,j,k) - h(i+1,j,k)) < 0.0) then
!                   call MOM_error(WARNING, "Corrective u-transport is too large.", .true.)
!                 endif
!               endif
!             endif
!           else ! v-point
!             if ((h(i,j+1,k) - h(i,j,k)) * &
!                 ((e(i,j+1,K)-e(i,j+1,K+1)) - (e(i,j,K)-e(i,j,K+1))) > 0.0) then
!               Sfn(K) = -Kh(i,K) * (e(i,j+1,K)-e(i,j,K)) * G%IdyCv(i,J)
!               Sfn(K+1) = -Kh(i,K+1) * (e(i,j+1,K+1)-e(i,j,K+1)) * G%IdyCv(i,J)
!               uh_here(k) = (Sfn(K) - Sfn(K+1))*G%dx_Cv(i,J)
!               if (abs(uh_here(K)) * min(G%IareaT(i,j), G%IareaT(i,j+1)) > &
!                   (1e-10*GV%m_to_H)) then
!                 if (uh_here(K) * (h(i,j+1,k) - h(i,j,k)) > 0.0) then
!                   call MOM_error(WARNING, &
!                          "Corrective v-transport is up the thickness gradient.", .true.)
!                 endif
!                 if (((h(i,j,k) - 4.0*dt*G%IareaT(i,j)*uh_here(K)) - &
!                      (h(i,j+1,k) + 4.0*dt*G%IareaT(i,j+1)*uh_here(K))) * &
!                     (h(i,j,k) - h(i,j+1,k)) < 0.0) then
!                   call MOM_error(WARNING, &
!                          "Corrective v-transport is too large.", .true.)
!                 endif
!               endif
!             endif
!           endif ! u- or v- selection.
!          !  de_dx(I,K) = (e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j)
!         endif
!       enddo
!     enddo

      if (n==1) then ! This is a u-column.
        do K=k_top+1,nz ; do i=ish,ie
          if (Kh(I,K) > KH_u(I,j,K)) then
            dKh = (Kh(I,K) - KH_u(I,j,K))
            int_slope_u(I,j,K) = dKh / Kh(I,K)
            KH_u(I,j,K) = Kh(I,K)
          endif
        enddo ; enddo
      else ! This is a v-column.
        do K=k_top+1,nz ; do i=ish,ie
          if (Kh(i,K) > KH_v(i,J,K)) then
            dKh = Kh(i,K) - KH_v(i,J,K)
            int_slope_v(i,J,K) = dKh / Kh(i,K)
            KH_v(i,J,K) = Kh(i,K)
          endif
        enddo ; enddo
      endif

    enddo ! j-loop
  enddo  ! n-loop over u- and v- directions.

end subroutine add_detangling_Kh

!> Initialize the isopycnal height diffusion module and its control structure
subroutine thickness_diffuse_init(Time, G, GV, US, param_file, diag, CDp, CS)
  type(time_type),         intent(in) :: Time    !< Current model time
  type(ocean_grid_type),   intent(in) :: G       !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV      !< Vertical grid structure
  type(unit_scale_type),   intent(in) :: US      !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(cont_diag_ptrs),    intent(inout) :: CDp  !< Continuity equation diagnostics
  type(thickness_diffuse_CS), intent(inout) :: CS !< Control structure for thickness_diffuse

  ! Local variables
  character(len=40)  :: mdl = "MOM_thickness_diffuse" ! This module's name.
  character(len=200) :: khth_file, inputdir, khth_varname
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  real :: grid_sp      ! The local grid spacing [L ~> m]
  real :: omega        ! The Earth's rotation rate [T-1 ~> s-1]
  real :: strat_floor  ! A floor for buoyancy frequency in the Ferrari et al. 2010,
                       ! streamfunction formulation, expressed as a fraction of planetary
                       ! rotation divided by an aspect ratio rescaling factor [L Z-1 ~> nondim]
  real :: Stanley_coeff ! Coefficient relating the temperature gradient and sub-gridscale
                        ! temperature variance [nondim]
  logical :: khth_use_ebt_struct ! If true, uses the equivalent barotropic structure
                                 ! as the vertical structure of thickness diffusivity.
                                 ! Used to determine if FULL_DEPTH_KHTH_MIN should be
                                 ! available.
  logical :: use_meke = .false. ! If true, use the MEKE formulation for the thickness diffusivity.
  integer :: default_answer_date ! The default setting for the various ANSWER_DATE flags.
  integer :: i, j

  CS%initialized = .true.
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "THICKNESSDIFFUSE", CS%thickness_diffuse, &
                 "If true, interface heights are diffused with a "//&
                 "coefficient of KHTH.", default=.false.)
  call get_param(param_file, mdl, "KHTH", CS%Khth, &
                 "The background horizontal thickness diffusivity.", &
                 default=0.0, units="m2 s-1", scale=US%m_to_L**2*US%T_to_s)
  call get_param(param_file, mdl, "READ_KHTH", CS%read_khth, &
                 "If true, read a file (given by KHTH_FILE) containing the "//&
                 "spatially varying horizontal isopycnal height diffusivity.", &
                 default=.false.)
  if (CS%read_khth) then
    if (CS%Khth > 0) then
        call MOM_error(FATAL, "thickness_diffuse_init: KHTH > 0 is not "// &
              "compatible with READ_KHTH = TRUE. ")
    endif
    call get_param(param_file, mdl, "INPUTDIR", inputdir, &
                 "The directory in which all input files are found.", &
                 default=".", do_not_log=.true.)
    inputdir = slasher(inputdir)
    call get_param(param_file, mdl, "KHTH_FILE", khth_file, &
                 "The file containing the spatially varying horizontal "//&
                 "isopycnal height diffusivity.", default="khth.nc")
    call get_param(param_file, mdl, "KHTH_VARIABLE", khth_varname, &
                 "The name of the isopycnal height diffusivity variable to read "//&
                 "from KHTH_FILE.", &
                 default="khth")
    khth_file = trim(inputdir) // trim(khth_file)

    allocate(CS%khth2d(G%isd:G%ied, G%jsd:G%jed), source=0.0)
    call MOM_read_data(khth_file, khth_varname, CS%khth2d(:,:), G%domain, scale=US%m_to_L**2*US%T_to_s)
    call pass_var(CS%khth2d, G%domain)
  endif
  call get_param(param_file, mdl, "KHTH_SLOPE_CFF", CS%KHTH_Slope_Cff, &
                 "The nondimensional coefficient in the Visbeck formula for "//&
                 "the interface depth diffusivity", units="nondim", default=0.0)
  call get_param(param_file, mdl, "KHTH_MIN", CS%KHTH_Min, &
                 "The minimum horizontal thickness diffusivity.", &
                 default=0.0, units="m2 s-1", scale=US%m_to_L**2*US%T_to_s)
  call get_param(param_file, mdl, "KHTH_USE_EBT_STRUCT", khth_use_ebt_struct, &
                 "If true, uses the equivalent barotropic structure "//&
                 "as the vertical structure of thickness diffusivity.",&
                 default=.false., do_not_log=.true.)
  if (khth_use_ebt_struct .and. CS%KHTH_Min>0.0) then
    call get_param(param_file, mdl, "FULL_DEPTH_KHTH_MIN", CS%full_depth_khth_min, &
                   "If true, KHTH_MIN is enforced throughout the whole water column. "//&
                   "Otherwise, KHTH_MIN is only enforced at the surface. This parameter "//&
                   "is only available when KHTH_USE_EBT_STRUCT=True and KHTH_MIN>0.",      &
                   default=.false.)
  endif
  call get_param(param_file, mdl, "KHTH_MAX", CS%KHTH_Max, &
                 "The maximum horizontal thickness diffusivity.", &
                 default=0.0, units="m2 s-1", scale=US%m_to_L**2*US%T_to_s)
  call get_param(param_file, mdl, "KHTH_MAX_CFL", CS%max_Khth_CFL, &
                 "The maximum value of the local diffusive CFL ratio that "//&
                 "is permitted for the thickness diffusivity. 1.0 is the "//&
                 "marginally unstable value in a pure layered model, but "//&
                 "much smaller numbers (e.g. 0.1) seem to work better for "//&
                 "ALE-based models.", units="nondimensional", default=0.8)

  call get_param(param_file, mdl, "KH_ETA_CONST", CS%Kh_eta_bg, &
                 "The background horizontal diffusivity of the interface heights (without "//&
                 "considering the layer density structure).  If diffusive CFL limits are "//&
                 "encountered, the diffusivities of the isopycnals and the interfaces heights "//&
                 "are scaled back proportionately.", &
                 default=0.0, units="m2 s-1", scale=US%m_to_L**2*US%T_to_s)
  call get_param(param_file, mdl, "KH_ETA_VEL_SCALE", CS%Kh_eta_vel, &
                 "A velocity scale that is multiplied by the grid spacing to give a contribution "//&
                 "to the horizontal diffusivity of the interface heights (without considering "//&
                 "the layer density structure).", &
                 default=0.0, units="m s-1", scale=US%m_to_L*US%T_to_s)

  if ((CS%Kh_eta_bg > 0.0) .or. (CS%Kh_eta_vel > 0.0)) then
    allocate(CS%Kh_eta_u(G%IsdB:G%IedB, G%jsd:G%jed), source=0.)
    allocate(CS%Kh_eta_v(G%isd:G%ied, G%JsdB:G%JedB), source=0.)
    do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
      grid_sp = sqrt((2.0*G%dxCu(I,j)**2 * G%dyCu(I,j)**2) / ((G%dxCu(I,j)**2) + (G%dyCu(I,j)**2)))
      CS%Kh_eta_u(I,j) = G%OBCmaskCu(I,j) * MAX(0.0, CS%Kh_eta_bg + CS%Kh_eta_vel * grid_sp)
    enddo ; enddo
    do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
      grid_sp = sqrt((2.0*G%dxCv(i,J)**2 * G%dyCv(i,J)**2) / ((G%dxCv(i,J)**2) + (G%dyCv(i,J)**2)))
      CS%Kh_eta_v(i,J) = G%OBCmaskCv(i,J) * MAX(0.0, CS%Kh_eta_bg + CS%Kh_eta_vel * grid_sp)
    enddo ; enddo
  endif

  if (CS%max_Khth_CFL < 0.0) CS%max_Khth_CFL = 0.0
  call get_param(param_file, mdl, "DETANGLE_INTERFACES", CS%detangle_interfaces, &
                 "If defined add 3-d structured enhanced interface height "//&
                 "diffusivities to horizontally smooth jagged layers.", &
                 default=.false.)
  CS%detangle_time = 0.0
  if (CS%detangle_interfaces) &
    call get_param(param_file, mdl, "DETANGLE_TIMESCALE", CS%detangle_time, &
                 "A timescale over which maximally jagged grid-scale "//&
                 "thickness variations are suppressed.  This must be "//&
                 "longer than DT, or 0 to use DT.", units="s", default=0.0, scale=US%s_to_T)
  call get_param(param_file, mdl, "KHTH_SLOPE_MAX", CS%slope_max, &
                 "A slope beyond which the calculated isopycnal slope is "//&
                 "not reliable and is scaled away.", units="nondim", default=0.01, scale=US%L_to_Z)
  call get_param(param_file, mdl, "KD_SMOOTH", CS%kappa_smooth, &
                 "A diapycnal diffusivity that is used to interpolate "//&
                 "more sensible values of T & S into thin layers.", &
                 units="m2 s-1", default=1.0e-6, scale=GV%m2_s_to_HZ_T)
  call get_param(param_file, mdl, "KHTH_USE_FGNV_STREAMFUNCTION", CS%use_FGNV_streamfn, &
                 "If true, use the streamfunction formulation of "//&
                 "Ferrari et al., 2010, which effectively emphasizes "//&
                 "graver vertical modes by smoothing in the vertical.",  &
                 default=.false.)
  call get_param(param_file, mdl, "FGNV_FILTER_SCALE", CS%FGNV_scale, &
                 "A coefficient scaling the vertical smoothing term in the "//&
                 "Ferrari et al., 2010, streamfunction formulation.", &
                 units="nondim", default=1., do_not_log=.not.CS%use_FGNV_streamfn)
  call get_param(param_file, mdl, "FGNV_C_MIN", CS%FGNV_c_min, &
                 "A minium wave speed used in the Ferrari et al., 2010, "//&
                 "streamfunction formulation.", &
                 default=0., units="m s-1", scale=US%m_s_to_L_T, do_not_log=.not.CS%use_FGNV_streamfn)
  call get_param(param_file, mdl, "FGNV_STRAT_FLOOR", strat_floor, &
                 "A floor for Brunt-Vasaila frequency in the Ferrari et al., 2010, "//&
                 "streamfunction formulation, expressed as a fraction of planetary "//&
                 "rotation, OMEGA. This should be tiny but non-zero to avoid degeneracy.", &
                 default=1.e-15, units="nondim", scale=US%Z_to_L, do_not_log=.not.CS%use_FGNV_streamfn)
  call get_param(param_file, mdl, "USE_STANLEY_GM", CS%use_stanley_gm, &
                 "If true, turn on Stanley SGS T variance parameterization "// &
                 "in GM code.", default=.false.)
  if (CS%use_stanley_gm) then
    call get_param(param_file, mdl, "STANLEY_COEFF", Stanley_coeff, &
                 "Coefficient correlating the temperature gradient and SGS T variance.", &
                 units="nondim", default=-1.0, do_not_log=.true.)
    if (Stanley_coeff < 0.0) call MOM_error(FATAL, &
                 "STANLEY_COEFF must be set >= 0 if USE_STANLEY_GM is true.")
  endif
  call get_param(param_file, mdl, "OMEGA", omega, &
                 "The rotation rate of the earth.", &
                 default=7.2921e-5, units="s-1", scale=US%T_to_s, do_not_log=.not.CS%use_FGNV_streamfn)
  CS%N2_floor = 0.
  if (CS%use_FGNV_streamfn) CS%N2_floor = (strat_floor*omega)**2
  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)

  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231, do_not_log=.true.)

  call get_param(param_file, mdl, "MEKE_GM_SRC_ALT", CS%GM_src_alt, &
                 "If true, use the GM energy conversion form S^2*N^2*kappa rather "//&
                 "than the streamfunction for the GM source term.", default=.false.)
  call get_param(param_file, mdl, "MEKE_GM_SRC_ANSWER_DATE", CS%MEKE_src_answer_date, &
                 "The vintage of the expressions in the GM energy conversion calculation when "//&
                 "MEKE_GM_SRC_ALT is true.  Values below 20240601 recover the answers from the "//&
                 "original implementation, while higher values use expressions that satisfy "//&
                 "rotational symmetry.", &
                 default=default_answer_date, do_not_log=.not.CS%GM_src_alt)
  call get_param(param_file, mdl, "MEKE_GM_SRC_ALT_SLOPE_BUG", CS%MEKE_src_slope_bug, &
                 "If true, use a bug that limits the positive values, but not the negative values, "//&
                 "of the slopes used when MEKE_GM_SRC_ALT is true.  When this is true, it breaks "//&
                 "all of the symmetry rules that MOM6 is supposed to obey.", &
                 default=.false., do_not_log=.not.CS%GM_src_alt)

  call get_param(param_file, mdl, "MEKE_GEOMETRIC", CS%MEKE_GEOMETRIC, &
                 "If true, uses the GM coefficient formulation from the GEOMETRIC "//&
                 "framework (Marshall et al., 2012).", default=.false.)
  if (CS%MEKE_GEOMETRIC) then
    call get_param(param_file, mdl, "MEKE_GEOMETRIC_EPSILON", CS%MEKE_GEOMETRIC_epsilon, &
                 "Minimum Eady growth rate used in the calculation of GEOMETRIC "//&
                 "thickness diffusivity.", units="s-1", default=1.0e-7, scale=US%T_to_s)
    call get_param(param_file, mdl, "MEKE_GEOMETRIC_ALPHA", CS%MEKE_GEOMETRIC_alpha, &
                 "The nondimensional coefficient governing the efficiency of the GEOMETRIC "//&
                 "thickness diffusion.", units="nondim", default=0.05)

    call get_param(param_file, mdl, "MEKE_GEOMETRIC_ANSWER_DATE", CS%MEKE_GEOM_answer_date, &
                 "The vintage of the expressions in the MEKE_GEOMETRIC calculation.  "//&
                 "Values below 20190101 recover the answers from the original implementation, "//&
                 "while higher values use expressions that satisfy rotational symmetry.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
    if (.not.GV%Boussinesq) CS%MEKE_GEOM_answer_date = max(CS%MEKE_GEOM_answer_date, 20230701)
  endif

  call get_param(param_file, mdl, "USE_MEKE", use_meke, default=.false., do_not_log=.true.)
  if (use_meke) then
    call get_param(param_file, mdl, "USE_KH_IN_MEKE", CS%Use_KH_in_MEKE, &
                   "If true, uses the thickness diffusivity calculated here to diffuse MEKE.", &
                   default=.false.)
    call get_param(param_file, mdl, "MEKE_MIN_DEPTH_DIFF", CS%MEKE_min_depth_diff, &
                   "The minimum total depth over which to average the diffusivity used for MEKE.  "//&
                   "When the total depth is less than this, the diffusivity is scaled away.", &
                   units="m", default=1.0, scale=GV%m_to_H, do_not_log=.not.CS%Use_KH_in_MEKE)
  endif

  call get_param(param_file, mdl, "USE_GME", CS%use_GME_thickness_diffuse, &
                 "If true, use the GM+E backscatter scheme in association "//&
                 "with the Gent and McWilliams parameterization.", default=.false.)

  call get_param(param_file, mdl, "USE_GM_WORK_BUG", CS%use_GM_work_bug, &
                 "If true, compute the top-layer work tendency on the u-grid "//&
                 "with the incorrect sign, for legacy reproducibility.", &
                 default=.false.)

  if (CS%use_GME_thickness_diffuse) then
    allocate(CS%KH_u_GME(G%IsdB:G%IedB, G%jsd:G%jed, GV%ke+1), source=0.)
    allocate(CS%KH_v_GME(G%isd:G%ied, G%JsdB:G%JedB, GV%ke+1), source=0.)
  endif

  CS%id_uhGM = register_diag_field('ocean_model', 'uhGM', diag%axesCuL, Time, &
           'Time Mean Diffusive Zonal Thickness Flux', &
           'kg s-1', conversion=GV%H_to_kg_m2*US%L_to_m**2*US%s_to_T, &
           y_cell_method='sum', v_extensive=.true.)
  if (CS%id_uhGM > 0) call safe_alloc_ptr(CDp%uhGM,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  CS%id_vhGM = register_diag_field('ocean_model', 'vhGM', diag%axesCvL, Time, &
           'Time Mean Diffusive Meridional Thickness Flux', &
           'kg s-1', conversion=GV%H_to_kg_m2*US%L_to_m**2*US%s_to_T, &
           x_cell_method='sum', v_extensive=.true.)
  if (CS%id_vhGM > 0) call safe_alloc_ptr(CDp%vhGM,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)

  CS%id_GMwork = register_diag_field('ocean_model', 'GMwork', diag%axesT1, Time, &
          'Integrated Tendency of Ocean Mesoscale Eddy KE from Parameterized Eddy Advection', &
          'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2, cmor_field_name='tnkebto', &
          cmor_long_name='Integrated Tendency of Ocean Mesoscale Eddy KE from Parameterized Eddy Advection', &
          cmor_standard_name='tendency_of_ocean_eddy_kinetic_energy_content_due_to_parameterized_eddy_advection')
  if (CS%id_GMwork > 0) &
    allocate(CS%GMwork(G%isd:G%ied,G%jsd:G%jed), source=0.)

  CS%id_KH_u = register_diag_field('ocean_model', 'KHTH_u', diag%axesCui, Time, &
           'Parameterized mesoscale eddy advection diffusivity at U-point', &
           'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  CS%id_KH_v = register_diag_field('ocean_model', 'KHTH_v', diag%axesCvi, Time, &
           'Parameterized mesoscale eddy advection diffusivity at V-point', &
           'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  CS%id_KH_t = register_diag_field('ocean_model', 'KHTH_t', diag%axesTL, Time, &
          'Ocean Tracer Diffusivity due to Parameterized Mesoscale Advection', &
          'm2 s-1', conversion=US%L_to_m**2*US%s_to_T, &
          cmor_field_name='diftrblo', &
          cmor_long_name='Ocean Tracer Diffusivity due to Parameterized Mesoscale Advection', &
          cmor_standard_name='ocean_tracer_diffusivity_due_to_parameterized_mesoscale_advection')

  CS%id_KH_u1 = register_diag_field('ocean_model', 'KHTH_u1', diag%axesCu1, Time,         &
           'Parameterized mesoscale eddy advection diffusivity at U-points (2-D)', &
           'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  CS%id_KH_v1 = register_diag_field('ocean_model', 'KHTH_v1', diag%axesCv1, Time,         &
           'Parameterized mesoscale eddy advection diffusivity at V-points (2-D)', &
           'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  CS%id_KH_t1 = register_diag_field('ocean_model', 'KHTH_t1', diag%axesT1, Time, &
           'Parameterized mesoscale eddy advection diffusivity at T-points (2-D)', &
           'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)

  CS%id_slope_x =  register_diag_field('ocean_model', 'neutral_slope_x', diag%axesCui, Time, &
           'Zonal slope of neutral surface', 'nondim', conversion=US%Z_to_L)
  if (CS%id_slope_x > 0) &
    allocate(CS%diagSlopeX(G%IsdB:G%IedB,G%jsd:G%jed,GV%ke+1), source=0.)

  CS%id_slope_y =  register_diag_field('ocean_model', 'neutral_slope_y', diag%axesCvi, Time, &
           'Meridional slope of neutral surface', 'nondim', conversion=US%Z_to_L)
  if (CS%id_slope_y > 0) &
    allocate(CS%diagSlopeY(G%isd:G%ied,G%JsdB:G%JedB,GV%ke+1), source=0.)

  CS%id_sfn_x =  register_diag_field('ocean_model', 'GM_sfn_x', diag%axesCui, Time, &
           'Parameterized Zonal Overturning Streamfunction', &
           'm3 s-1', conversion=GV%H_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_sfn_y =  register_diag_field('ocean_model', 'GM_sfn_y', diag%axesCvi, Time, &
           'Parameterized Meridional Overturning Streamfunction', &
           'm3 s-1', conversion=GV%H_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_sfn_unlim_x =  register_diag_field('ocean_model', 'GM_sfn_unlim_x', diag%axesCui, Time, &
           'Parameterized Zonal Overturning Streamfunction before limiting/smoothing', &
           'm3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_sfn_unlim_y =  register_diag_field('ocean_model', 'GM_sfn_unlim_y', diag%axesCvi, Time, &
           'Parameterized Meridional Overturning Streamfunction before limiting/smoothing', &
           'm3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)

end subroutine thickness_diffuse_init

!> Copies KH_u_GME and KH_v_GME from private type into arrays provided as arguments
subroutine thickness_diffuse_get_KH(CS, KH_u_GME, KH_v_GME, G, GV)
  type(thickness_diffuse_CS),          intent(in)  :: CS   !< Control structure for this module
  type(ocean_grid_type),               intent(in)  :: G    !< Grid structure
  type(verticalGrid_type),             intent(in)  :: GV   !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: KH_u_GME !< Isopycnal height
                                                   !! diffusivities at u-faces [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: KH_v_GME !< Isopycnal height
                                                   !! diffusivities at v-faces [L2 T-1 ~> m2 s-1]
  ! Local variables
  integer :: i,j,k

  do k=1,GV%ke+1 ; do j = G%jsc, G%jec ; do I = G%isc-1, G%iec
    KH_u_GME(I,j,k) = CS%KH_u_GME(I,j,k)
  enddo ; enddo ; enddo

  do k=1,GV%ke+1 ; do J = G%jsc-1, G%jec ; do i = G%isc, G%iec
    KH_v_GME(i,J,k) = CS%KH_v_GME(i,J,k)
  enddo ; enddo ; enddo

end subroutine thickness_diffuse_get_KH

!> Deallocate the thickness_diffus3 control structure
subroutine thickness_diffuse_end(CS, CDp)
  type(thickness_diffuse_CS), intent(inout) :: CS !< Control structure for thickness_diffuse
  type(cont_diag_ptrs), intent(inout) :: CDp      !< Continuity diagnostic control structure

  if (CS%id_slope_x > 0) deallocate(CS%diagSlopeX)
  if (CS%id_slope_y > 0) deallocate(CS%diagSlopeY)

  if (CS%id_GMwork > 0) deallocate(CS%GMwork)

  ! NOTE: [uv]hGM may be allocated either here or the diagnostic module
  if (associated(CDp%uhGM)) deallocate(CDp%uhGM)
  if (associated(CDp%vhGM)) deallocate(CDp%vhGM)

  if (CS%use_GME_thickness_diffuse) then
    deallocate(CS%KH_u_GME)
    deallocate(CS%KH_v_GME)
  endif

  if (allocated(CS%khth2d)) deallocate(CS%khth2d)
end subroutine thickness_diffuse_end

!> \namespace mom_thickness_diffuse
!!
!! \section section_gm Isopycnal height diffusion (aka Gent-McWilliams)
!!
!! Isopycnal height diffusion is implemented via along-layer mass fluxes
!! \f[
!! h^\dagger \leftarrow h^n - \Delta t \nabla \cdot ( \vec{uh}^* )
!! \f]
!! where the mass fluxes are cast as the difference in vector streamfunction
!!
!! \f[
!! \vec{uh}^* = \delta_k \vec{\psi} .
!! \f]
!!
!! The GM implementation of isopycnal height diffusion made the streamfunction proportional
!! to the potential density slope
!! \f[
!! \vec{\psi} = - \kappa_h \frac{\nabla_z \rho}{\partial_z \rho}
!! = \frac{g\kappa_h}{\rho_o} \frac{\nabla \rho}{N^2} = -\kappa_h \frac{M^2}{N^2}
!! \f]
!! but for robustness the scheme is implemented as
!! \f[
!! \vec{\psi} = -\kappa_h \frac{M^2}{\sqrt{N^4 + M^4}}
!! \f]
!! since the quantity \f$\frac{M^2}{\sqrt{N^4 + M^4}}\f$ is bounded between $-1$ and $1$ and does not change sign
!! if \f$N^2<0\f$.
!!
!! Optionally, the method of \cite ferrari2010, can be used to obtain the streamfunction which solves the
!! vertically elliptic equation:
!! \f[
!! \gamma_F \partial_z c^2 \partial_z \psi - N_*^2 \psi  = -( 1 + \gamma_F ) \kappa_h N_*^2 \frac{M^2}{\sqrt{N^4+M^4}}
!! \f]
!! which recovers the previous streamfunction relation in the limit that \f$ c \rightarrow 0 \f$.
!! Here, \f$c=\max(c_{min},c_g)\f$ is the maximum of either \f$c_{min}\f$ and either the first baroclinic mode
!! wave-speed or the equivalent barotropic mode wave-speed.
!! \f$N_*^2 = \max(N^2,0)\f$ is a non-negative form of the square of the buoyancy frequency.
!! The parameter \f$\gamma_F\f$ is used to reduce the vertical smoothing length scale.
!! \f[
!! \kappa_h = \left( \kappa_o + \alpha_{s} L_{s}^2 < S N > + \alpha_{M} \kappa_{M} \right) r(\Delta x,L_d)
!! \f]
!! where \f$ S \f$ is the isoneutral slope magnitude, \f$ N \f$ is the buoyancy frequency,
!! \f$\kappa_{M}\f$ is the diffusivity calculated by the MEKE parameterization (mom_meke module) and
!! \f$ r(\Delta x,L_d) \f$ is a function of the local resolution (ratio of grid-spacing, \f$\Delta x\f$,
!! to deformation radius, \f$L_d\f$). The length \f$L_s\f$ is provided by the mom_lateral_mixing_coeffs module
!! (enabled with <code>USE_VARIABLE_MIXING=True</code> and the term \f$<SN>\f$ is the vertical average slope
!! times the buoyancy frequency prescribed by \cite visbeck1996.
!!
!! The result of the above expression is subsequently bounded by minimum and maximum values, including an upper
!! diffusivity consistent with numerical stability (\f$ \kappa_{cfl} \f$ is calculated internally).
!! \f[
!! \kappa_h \leftarrow \min{\left( \kappa_{max}, \kappa_{cfl}, \max{\left( \kappa_{min}, \kappa_h \right)} \right)}
!!                      f(c_g,z)
!! \f]
!!
!! where \f$f(c_g,z)\f$ is a vertical structure function.
!! \f$f(c_g,z)\f$ is calculated in module mom_lateral_mixing_coeffs.
!! If <code>KHTH_USE_EBT_STRUCT=True</code> then \f$f(c_g,z)\f$ is set to look like the equivalent barotropic
!! modal velocity structure. Otherwise \f$f(c_g,z)=1\f$ and the diffusivity is independent of depth.
!!
!! In order to calculate meaningful slopes in vanished layers, temporary copies of the thermodynamic variables
!! are passed through a vertical smoother, function vert_fill_ts():
!! \f{eqnarray*}{
!! \left[ 1 + \Delta t \kappa_{smth} \frac{\partial^2}{\partial_z^2} \right] \theta & \leftarrow & \theta \\
!! \left[ 1 + \Delta t \kappa_{smth} \frac{\partial^2}{\partial_z^2} \right] s & \leftarrow & s
!! \f}
!!
!! \subsection section_khth_module_parameters Module mom_thickness_diffuse parameters
!!
!! | Symbol                | Module parameter |
!! | ------                | --------------- |
!! | -                     | <code>THICKNESSDIFFUSE</code> |
!! | \f$ \kappa_o \f$      | <code>KHTH</code> |
!! | \f$ \alpha_{s} \f$    | <code>KHTH_SLOPE_CFF</code> |
!! | \f$ \kappa_{min} \f$  | <code>KHTH_MIN</code> |
!! | \f$ \kappa_{max} \f$  | <code>KHTH_MAX</code> |
!! | -                     | <code>KHTH_MAX_CFL</code> |
!! | \f$ \kappa_{smth} \f$ | <code>KD_SMOOTH</code> |
!! | \f$ \alpha_{M} \f$    | <code>MEKE_KHTH_FAC</code> (from mom_meke module) |
!! | -                     | <code>KHTH_USE_EBT_STRUCT</code> (from mom_lateral_mixing_coeffs module) |
!! | -                     | <code>KHTH_USE_FGNV_STREAMFUNCTION</code> |
!! | \f$ \gamma_F \f$      | <code>FGNV_FILTER_SCALE</code> |
!! | \f$ c_{min} \f$       | <code>FGNV_C_MIN</code> |
!!
!! \subsection section_khth_module_reference References
!!
!! Ferrari, R., S.M. Griffies, A.J.G. Nurser and G.K. Vallis, 2010:
!! A boundary-value problem for the parameterized mesoscale eddy transport.
!! Ocean Modelling, 32, 143-156. http://doi.org/10.1016/j.ocemod.2010.01.004
!!
!! Visbeck, M., J.C. Marshall, H. Jones, 1996:
!! Dynamics of isolated convective regions in the ocean. J. Phys. Oceangr., 26, 1721-1734.
!! http://dx.doi.org/10.1175/1520-0485(1996)026%3C1721:DOICRI%3E2.0.CO;2

end module MOM_thickness_diffuse
