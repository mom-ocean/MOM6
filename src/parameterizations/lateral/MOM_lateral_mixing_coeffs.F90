!> Variable mixing coefficients
module MOM_lateral_mixing_coeffs

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,         only : hchksum, uvchksum
use MOM_error_handler,     only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_diag_mediator,     only : register_diag_field, safe_alloc_ptr, post_data
use MOM_diag_mediator,     only : diag_ctrl, time_type, query_averaging_enabled
use MOM_domains,           only : create_group_pass, do_group_pass
use MOM_domains,           only : group_pass_type, pass_var, pass_vector
use MOM_file_parser,       only : get_param, log_version, param_file_type
use MOM_interface_heights, only : find_eta, thickness_to_dz
use MOM_isopycnal_slopes,  only : calc_isoneutral_slopes
use MOM_grid,              only : ocean_grid_type
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : thermo_var_ptrs
use MOM_verticalGrid,      only : verticalGrid_type
use MOM_wave_speed,        only : wave_speed, wave_speed_CS, wave_speed_init
use MOM_open_boundary,     only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary,     only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_MEKE_types,        only : MEKE_type

implicit none ; private

#include <MOM_memory.h>

!> Variable mixing coefficients
type, public :: VarMix_CS
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: use_variable_mixing  !< If true, use the variable mixing.
  logical :: Resoln_scaling_used  !< If true, a resolution function is used somewhere to scale
                                  !! away one of the viscosities or diffusivities when the
                                  !! deformation radius is well resolved.
  logical :: Resoln_scaled_Kh     !< If true, scale away the Laplacian viscosity
                                  !! when the deformation radius is well resolved.
  logical :: Resoln_scaled_KhTh   !< If true, scale away the thickness diffusivity
                                  !! when the deformation radius is well resolved.
  logical :: Depth_scaled_KhTh    !< If true, KHTH is scaled away when the depth is
                                  !! shallower than a reference depth.
  logical :: Resoln_scaled_KhTr   !< If true, scale away the tracer diffusivity
                                  !! when the deformation radius is well resolved.
  logical :: interpolate_Res_fn   !< If true, interpolate the resolution function
                                  !! to the velocity points from the thickness
                                  !! points; otherwise interpolate the wave
                                  !! speed and calculate the resolution function
                                  !! independently at each point.
  logical :: use_stored_slopes    !< If true, stores isopycnal slopes in this structure.
  logical :: Resoln_use_ebt       !< If true, uses the equivalent barotropic wave speed instead
                                  !! of first baroclinic wave for calculating the resolution fn.
  logical :: khth_use_ebt_struct  !< If true, uses the equivalent barotropic structure
                                  !! as the vertical structure of thickness diffusivity.
  logical :: kdgl90_use_ebt_struct  !< If true, uses the equivalent barotropic structure
                                  !! as the vertical structure of diffusivity in the GL90 scheme.
  logical :: kdgl90_use_sqg_struct  !< If true, uses the surface quasigeostrophic structure
                                  !! as the vertical structure of diffusivity in the GL90 scheme.
  logical :: khth_use_sqg_struct  !< If true, uses the surface quasigeostrophic structure
                                  !! as the vertical structure of thickness diffusivity.
  logical :: khtr_use_ebt_struct  !< If true, uses the equivalent barotropic structure
                                  !! as the vertical structure of tracer diffusivity.
  logical :: khtr_use_sqg_struct  !< If true, uses the surface quasigeostrophic structure
                                  !! as the vertical structure of tracer diffusivity.
  logical :: calculate_cg1        !< If true, calls wave_speed() to calculate the first
                                  !! baroclinic wave speed and populate CS%cg1.
                                  !! This parameter is set depending on other parameters.
  logical :: calculate_Rd_dx      !< If true, calculates Rd/dx and populate CS%Rd_dx_h.
                                  !! This parameter is set depending on other parameters.
  logical :: calculate_res_fns    !< If true, calculate all the resolution factors.
                                  !! This parameter is set depending on other parameters.
  logical :: calculate_depth_fns !< If true, calculate all the depth factors.
                                  !! This parameter is set depending on other parameters.
  logical :: calculate_Eady_growth_rate !< If true, calculate all the Eady growth rates.
                                  !! This parameter is set depending on other parameters.
  logical :: use_stanley_iso      !< If true, use Stanley parameterization in MOM_isopycnal_slopes
  logical :: use_simpler_Eady_growth_rate !< If true, use a simpler method to calculate the
                                  !! Eady growth rate that avoids division by layer thickness.
                                  !! This parameter is set depending on other parameters.
  logical :: full_depth_Eady_growth_rate !< If true, calculate the Eady growth rate based on an
                                  !! average that includes contributions from sea-level changes
                                  !! in its denominator, rather than just the nominal depth of
                                  !! the bathymetry.  This only applies when using the model
                                  !! interface heights as a proxy for isopycnal slopes.
  logical :: OBC_friendly         !< If true, use only interior data for thickness weighting and
                                  !! to calculate stratification and other fields at open boundary
                                  !!  condition faces.
  real :: cropping_distance       !< Distance from surface or bottom to filter out outcropped or
                                  !! incropped interfaces for the Eady growth rate calc [Z ~> m]
  real :: h_min_N2                !< The minimum vertical distance to use in the denominator of the
                                  !! buoyancy frequency used in the slope calculation [H ~> m or kg m-2]

  real, allocatable :: SN_u(:,:)      !< S*N at u-points [T-1 ~> s-1]
  real, allocatable :: SN_v(:,:)      !< S*N at v-points [T-1 ~> s-1]
  real, allocatable :: L2u(:,:)       !< Length scale^2 at u-points [L2 ~> m2]
  real, allocatable :: L2v(:,:)       !< Length scale^2 at v-points [L2 ~> m2]
  real, allocatable :: cg1(:,:)       !< The first baroclinic gravity wave speed [L T-1 ~> m s-1].
  real, allocatable :: Res_fn_h(:,:)  !< Non-dimensional function of the ratio the first baroclinic
                                      !! deformation radius to the grid spacing at h points [nondim].
  real, allocatable :: Res_fn_q(:,:)  !< Non-dimensional function of the ratio the first baroclinic
                                      !! deformation radius to the grid spacing at q points [nondim].
  real, allocatable :: Res_fn_u(:,:)  !< Non-dimensional function of the ratio the first baroclinic
                                      !! deformation radius to the grid spacing at u points [nondim].
  real, allocatable :: Res_fn_v(:,:)  !< Non-dimensional function of the ratio the first baroclinic
                                      !! deformation radius to the grid spacing at v points [nondim].
  real, allocatable :: Depth_fn_u(:,:) !< Non-dimensional function of the ratio of the depth to
                                      !! a reference depth (maximum 1) at u points [nondim]
  real, allocatable :: Depth_fn_v(:,:) !< Non-dimensional function of the ratio of the depth to
                                      !! a reference depth (maximum 1) at v points [nondim]
  real, allocatable :: beta_dx2_h(:,:) !< The magnitude of the gradient of the Coriolis parameter
                                      !! times the grid spacing squared at h points [L T-1 ~> m s-1].
  real, allocatable :: beta_dx2_q(:,:) !< The magnitude of the gradient of the Coriolis parameter
                                      !! times the grid spacing squared at q points [L T-1 ~> m s-1].
  real, allocatable :: beta_dx2_u(:,:) !< The magnitude of the gradient of the Coriolis parameter
                                      !! times the grid spacing squared at u points [L T-1 ~> m s-1].
  real, allocatable :: beta_dx2_v(:,:) !< The magnitude of the gradient of the Coriolis parameter
                                      !! times the grid spacing squared at v points [L T-1 ~> m s-1].
  real, allocatable :: f2_dx2_h(:,:)  !< The Coriolis parameter squared times the grid
                                      !! spacing squared at h [L2 T-2 ~> m2 s-2].
  real, allocatable :: f2_dx2_q(:,:)  !< The Coriolis parameter squared times the grid
                                      !! spacing squared at q [L2 T-2 ~> m2 s-2].
  real, allocatable :: f2_dx2_u(:,:)  !< The Coriolis parameter squared times the grid
                                      !! spacing squared at u [L2 T-2 ~> m2 s-2].
  real, allocatable :: f2_dx2_v(:,:)  !< The Coriolis parameter squared times the grid
                                      !! spacing squared at v [L2 T-2 ~> m2 s-2].
  real, allocatable :: Rd_dx_h(:,:)   !< Deformation radius over grid spacing [nondim]

  real, allocatable :: slope_x(:,:,:)     !< Zonal isopycnal slope [Z L-1 ~> nondim]
  real, allocatable :: slope_y(:,:,:)     !< Meridional isopycnal slope [Z L-1 ~> nondim]
  real, allocatable :: ebt_struct(:,:,:)  !< EBT vertical structure to scale diffusivities with [nondim]
  real, allocatable :: sqg_struct(:,:,:)  !< SQG vertical structure to scale diffusivities with [nondim]
  real, allocatable :: BS_struct(:,:,:) !< Vertical structure function used in backscatter [nondim]
  real, allocatable :: khth_struct(:,:,:) !< Vertical structure function used in thickness diffusivity [nondim]
  real, allocatable :: khtr_struct(:,:,:) !< Vertical structure function used in tracer diffusivity [nondim]
  real, allocatable :: kdgl90_struct(:,:,:) !< Vertical structure function used in GL90 diffusivity [nondim]
  real :: BS_EBT_power                !< Power to raise EBT vertical structure to. Default 0.0.
  real :: sqg_expo     !< Exponent for SQG vertical structure [nondim]. Default 1.0
  logical :: BS_use_sqg_struct   !< If true, use sqg_stuct for backscatter vertical structure.


  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    Laplac3_const_u       !< Laplacian metric-dependent constants [L3 ~> m3]

  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    Laplac3_const_v       !< Laplacian metric-dependent constants [L3 ~> m3]

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    KH_u_QG               !< QG Leith GM coefficient at u-points [L2 T-1 ~> m2 s-1]

  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    KH_v_QG               !< QG Leith GM coefficient at v-points [L2 T-1 ~> m2 s-1]

  ! Parameters
  logical :: use_Visbeck  !< Use Visbeck formulation for thickness diffusivity
  integer :: VarMix_Ktop  !< Top layer to start downward integrals
  real :: Visbeck_L_scale !< Fixed length scale in Visbeck formula [L ~> m], or if negative a scaling
                          !! factor [nondim] relating this length scale squared to the cell area
  real :: Eady_GR_D_scale !< Depth over which to average SN [Z ~> m]
  real :: Res_coef_khth   !< A coefficient [nondim] that determines the function
                          !! of resolution, used for thickness and tracer mixing, as:
                          !!  F = 1 / (1 + (Res_coef_khth*Ld/dx)^Res_fn_power)
  real :: Res_coef_visc   !< A coefficient [nondim] that determines the function
                          !! of resolution, used for lateral viscosity, as:
                          !!  F = 1 / (1 + (Res_coef_visc*Ld/dx)^Res_fn_power)
  real :: depth_scaled_khth_h0 !< The depth above which KHTH is linearly scaled away [Z ~> m]
  real :: depth_scaled_khth_exp !< The exponent used in the depth dependent scaling function for KHTH [nondim]
  real :: kappa_smooth    !< A diffusivity for smoothing T/S in vanished layers [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  integer :: Res_fn_power_khth !< The power of dx/Ld in the KhTh resolution function.  Any
                               !! positive integer power may be used, but even powers
                               !! and especially 2 are coded to be more efficient.
  integer :: Res_fn_power_visc !< The power of dx/Ld in the Kh resolution function.  Any
                               !! positive integer power may be used, but even powers
                               !! and especially 2 are coded to be more efficient.
  real :: Visbeck_S_max   !< Upper bound on slope used in Eady growth rate [Z L-1 ~> nondim].

  ! Leith parameters
  logical :: use_QG_Leith_GM      !< If true, uses the QG Leith viscosity as the GM coefficient
  logical :: use_beta_in_QG_Leith !< If true, includes the beta term in the QG Leith GM coefficient

  ! Diagnostics
  !>@{
  !! Diagnostic identifier
  integer :: id_SN_u=-1, id_SN_v=-1, id_L2u=-1, id_L2v=-1, id_Res_fn = -1
  integer :: id_N2_u=-1, id_N2_v=-1, id_S2_u=-1, id_S2_v=-1
  integer :: id_dzu=-1, id_dzv=-1, id_dzSxN=-1, id_dzSyN=-1
  integer :: id_Rd_dx=-1, id_KH_u_QG = -1, id_KH_v_QG = -1
  integer :: id_sqg_struct=-1, id_BS_struct=-1, id_khth_struct=-1, id_khtr_struct=-1
  integer :: id_kdgl90_struct=-1
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.
  !>@}

  type(wave_speed_CS) :: wave_speed !< Wave speed control structure
  type(group_pass_type) :: pass_cg1 !< For group halo pass
  logical :: debug      !< If true, write out checksums of data for debugging
end type VarMix_CS

public VarMix_init, VarMix_end, calc_slope_functions, calc_resoln_function
public calc_QG_slopes, calc_QG_Leith_viscosity, calc_depth_function, calc_sqg_struct

contains

!> Calculates the non-dimensional depth functions.
subroutine calc_depth_function(G, CS)
  type(ocean_grid_type),  intent(in)    :: G  !< Ocean grid structure
  type(VarMix_CS),        intent(inout) :: CS !< Variable mixing control structure

  ! Local variables
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: i, j
  real    :: H0   ! The depth above which KHTH is linearly scaled away [Z ~> m]
  real    :: expo ! exponent used in the depth dependent scaling [nondim]
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not. CS%initialized) call MOM_error(FATAL, "calc_depth_function: "// &
         "Module must be initialized before it is used.")

  if (.not. CS%calculate_depth_fns) return
  if (.not. allocated(CS%Depth_fn_u)) call MOM_error(FATAL, &
    "calc_depth_function: %Depth_fn_u is not associated with Depth_scaled_KhTh.")
  if (.not. allocated(CS%Depth_fn_v)) call MOM_error(FATAL, &
    "calc_depth_function: %Depth_fn_v is not associated with Depth_scaled_KhTh.")

  ! For efficiency, the reciprocal of H0 should be used instead.
  H0 = CS%depth_scaled_khth_h0
  expo = CS%depth_scaled_khth_exp
!$OMP do
  do j=js,je ; do I=is-1,Ieq
    CS%Depth_fn_u(I,j) = (MIN(1.0, (0.5*(G%bathyT(i,j) + G%bathyT(i+1,j)) + G%Z_ref)/H0))**expo
  enddo ; enddo
!$OMP do
  do J=js-1,Jeq ; do i=is,ie
    CS%Depth_fn_v(i,J) = (MIN(1.0, (0.5*(G%bathyT(i,j) + G%bathyT(i,j+1)) + G%Z_ref)/H0))**expo
  enddo ; enddo

end subroutine calc_depth_function

!> Calculates and stores the non-dimensional resolution functions
subroutine calc_resoln_function(h, tv, G, GV, US, CS, MEKE, OBC, dt)
  type(ocean_grid_type),                     intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(in)    :: tv !< Thermodynamic variables
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  type(VarMix_CS),                           intent(inout) :: CS !< Variable mixing control structure
  type(MEKE_type),                           intent(in)    :: MEKE !< MEKE struct
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure
  real,                                      intent(in)    :: dt !< Time increment [T ~> s]

  ! Local variables
  ! Depending on the power-function being used, dimensional rescaling may be limited, so some
  ! of the following variables have units that depend on that power.
  real :: cg1_q  ! The gravity wave speed interpolated to q points [L T-1 ~> m s-1] or [m s-1].
  real :: cg1_u  ! The gravity wave speed interpolated to u points [L T-1 ~> m s-1] or [m s-1].
  real :: cg1_v  ! The gravity wave speed interpolated to v points [L T-1 ~> m s-1] or [m s-1].
  real :: dx_term ! A term in the denominator [L2 T-2 ~> m2 s-2] or [m2 s-2]
  integer :: power_2
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not. CS%initialized) call MOM_error(FATAL, "calc_resoln_function: "// &
         "Module must be initialized before it is used.")

  if (CS%calculate_cg1) then
    if (.not. allocated(CS%cg1)) call MOM_error(FATAL, &
      "calc_resoln_function: %cg1 is not associated with Resoln_scaled_Kh.")
    if (CS%khth_use_ebt_struct .or. CS%kdgl90_use_ebt_struct &
     .or. CS%khtr_use_ebt_struct .or. CS%BS_EBT_power>0.) then
      if (.not. allocated(CS%ebt_struct)) call MOM_error(FATAL, &
        "calc_resoln_function: %ebt_struct is not associated with RESOLN_USE_EBT.")
      if (CS%Resoln_use_ebt) then
        ! Both resolution fn and vertical structure are using EBT
        call wave_speed(h, tv, G, GV, US, CS%cg1, CS%wave_speed, modal_structure=CS%ebt_struct)
      else
        ! Use EBT to get vertical structure first and then re-calculate cg1 using first baroclinic mode
        call wave_speed(h, tv, G, GV, US, CS%cg1, CS%wave_speed, modal_structure=CS%ebt_struct, &
                        use_ebt_mode=.true.)
        call wave_speed(h, tv, G, GV, US, CS%cg1, CS%wave_speed)
      endif
      call pass_var(CS%ebt_struct, G%Domain)
    else
      call wave_speed(h, tv, G, GV, US, CS%cg1, CS%wave_speed)
    endif

    call create_group_pass(CS%pass_cg1, CS%cg1, G%Domain)
    call do_group_pass(CS%pass_cg1, G%Domain)
  endif
  if (CS%BS_use_sqg_struct .or. CS%khth_use_sqg_struct .or. CS%khtr_use_sqg_struct &
      .or. CS%kdgl90_use_sqg_struct .or. CS%id_sqg_struct>0) then
    call calc_sqg_struct(h, tv, G, GV, US, CS, dt, MEKE, OBC)
    call pass_var(CS%sqg_struct, G%Domain)
  endif

  if (CS%BS_EBT_power>0.) then
    do k=1,nz ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
      CS%BS_struct(i,j,k) = CS%ebt_struct(i,j,k)**CS%BS_EBT_power
    enddo ; enddo ; enddo
  elseif (CS%BS_use_sqg_struct) then
    do k=1,nz ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
      CS%BS_struct(i,j,k) = CS%sqg_struct(i,j,k)
    enddo ; enddo ; enddo
  endif

  if (CS%khth_use_ebt_struct) then
    do k=1,nz ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
      CS%khth_struct(i,j,k) = CS%ebt_struct(i,j,k)
    enddo ; enddo ; enddo
  elseif (CS%khth_use_sqg_struct) then
    do k=1,nz ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
      CS%khth_struct(i,j,k) = CS%sqg_struct(i,j,k)
    enddo ; enddo ; enddo
  endif

  if (CS%khtr_use_ebt_struct) then
    do k=1,nz ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
      CS%khtr_struct(i,j,k) = CS%ebt_struct(i,j,k)
    enddo ; enddo ; enddo
  elseif (CS%khtr_use_sqg_struct) then
    do k=1,nz ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
      CS%khtr_struct(i,j,k) = CS%sqg_struct(i,j,k)
    enddo ; enddo ; enddo
  endif

  if (CS%kdgl90_use_ebt_struct) then
    do k=1,nz ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
      CS%kdgl90_struct(i,j,k) = CS%ebt_struct(i,j,k)
    enddo ; enddo ; enddo
  elseif (CS%kdgl90_use_sqg_struct) then
    do k=1,nz ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
      CS%kdgl90_struct(i,j,k) = CS%sqg_struct(i,j,k)
    enddo ; enddo ; enddo
  endif

  ! Calculate and store the ratio between deformation radius and grid-spacing
  ! at h-points [nondim].
  if (CS%calculate_rd_dx) then
    if (.not. allocated(CS%Rd_dx_h)) call MOM_error(FATAL, &
      "calc_resoln_function: %Rd_dx_h is not associated with calculate_rd_dx.")
    !$OMP parallel do default(shared)
    do j=js-1,je+1 ; do i=is-1,ie+1
      CS%Rd_dx_h(i,j) = CS%cg1(i,j) / &
            (sqrt(CS%f2_dx2_h(i,j) + CS%cg1(i,j)*CS%beta_dx2_h(i,j)))
    enddo ; enddo
    if (query_averaging_enabled(CS%diag)) then
      if (CS%id_Rd_dx > 0) call post_data(CS%id_Rd_dx, CS%Rd_dx_h, CS%diag)
    endif
  endif

  if (.not. CS%calculate_res_fns) return

  if (.not. allocated(CS%Res_fn_h)) call MOM_error(FATAL, &
    "calc_resoln_function: %Res_fn_h is not associated with Resoln_scaled_Kh.")
  if (.not. allocated(CS%Res_fn_q)) call MOM_error(FATAL, &
    "calc_resoln_function: %Res_fn_q is not associated with Resoln_scaled_Kh.")
  if (.not. allocated(CS%Res_fn_u)) call MOM_error(FATAL, &
    "calc_resoln_function: %Res_fn_u is not associated with Resoln_scaled_Kh.")
  if (.not. allocated(CS%Res_fn_v)) call MOM_error(FATAL, &
    "calc_resoln_function: %Res_fn_v is not associated with Resoln_scaled_Kh.")
  if (.not. allocated(CS%f2_dx2_h)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_h is not associated with Resoln_scaled_Kh.")
  if (.not. allocated(CS%f2_dx2_q)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_q is not associated with Resoln_scaled_Kh.")
  if (.not. allocated(CS%f2_dx2_u)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_u is not associated with Resoln_scaled_Kh.")
  if (.not. allocated(CS%f2_dx2_v)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_v is not associated with Resoln_scaled_Kh.")
  if (.not. allocated(CS%beta_dx2_h)) call MOM_error(FATAL, &
    "calc_resoln_function: %beta_dx2_h is not associated with Resoln_scaled_Kh.")
  if (.not. allocated(CS%beta_dx2_q)) call MOM_error(FATAL, &
    "calc_resoln_function: %beta_dx2_q is not associated with Resoln_scaled_Kh.")
  if (.not. allocated(CS%beta_dx2_u)) call MOM_error(FATAL, &
    "calc_resoln_function: %beta_dx2_u is not associated with Resoln_scaled_Kh.")
  if (.not. allocated(CS%beta_dx2_v)) call MOM_error(FATAL, &
    "calc_resoln_function: %beta_dx2_v is not associated with Resoln_scaled_Kh.")

  !   Do this calculation on the extent used in MOM_hor_visc.F90, and
  ! MOM_tracer.F90 so that no halo update is needed.

!$OMP parallel default(none) shared(is,ie,js,je,Ieq,Jeq,CS,US) &
!$OMP                       private(dx_term,cg1_q,power_2,cg1_u,cg1_v)
  if (CS%Res_fn_power_visc >= 100) then
!$OMP do
    do j=js-1,je+1 ; do i=is-1,ie+1
      dx_term = CS%f2_dx2_h(i,j) + CS%cg1(i,j)*CS%beta_dx2_h(i,j)
      if ((CS%Res_coef_visc * CS%cg1(i,j))**2 > dx_term) then
        CS%Res_fn_h(i,j) = 0.0
      else
        CS%Res_fn_h(i,j) = 1.0
      endif
    enddo ; enddo
!$OMP do
    do J=js-1,Jeq ; do I=is-1,Ieq
      cg1_q = 0.25 * ((CS%cg1(i,j) + CS%cg1(i+1,j+1)) + (CS%cg1(i+1,j) + CS%cg1(i,j+1)))
      dx_term = CS%f2_dx2_q(I,J) +  cg1_q * CS%beta_dx2_q(I,J)
      if ((CS%Res_coef_visc * cg1_q)**2 > dx_term) then
        CS%Res_fn_q(I,J) = 0.0
      else
        CS%Res_fn_q(I,J) = 1.0
      endif
    enddo ; enddo
  elseif (CS%Res_fn_power_visc == 2) then
!$OMP do
    do j=js-1,je+1 ; do i=is-1,ie+1
      dx_term = CS%f2_dx2_h(i,j) + CS%cg1(i,j)*CS%beta_dx2_h(i,j)
      CS%Res_fn_h(i,j) = dx_term / (dx_term + (CS%Res_coef_visc * CS%cg1(i,j))**2)
    enddo ; enddo
!$OMP do
    do J=js-1,Jeq ; do I=is-1,Ieq
      cg1_q = 0.25 * ((CS%cg1(i,j) + CS%cg1(i+1,j+1)) + (CS%cg1(i+1,j) + CS%cg1(i,j+1)))
      dx_term = CS%f2_dx2_q(I,J) +  cg1_q * CS%beta_dx2_q(I,J)
      CS%Res_fn_q(I,J) = dx_term / (dx_term + (CS%Res_coef_visc * cg1_q)**2)
    enddo ; enddo
  elseif (mod(CS%Res_fn_power_visc, 2) == 0) then
    power_2 = CS%Res_fn_power_visc / 2
!$OMP do
    do j=js-1,je+1 ; do i=is-1,ie+1
      dx_term = (US%L_T_to_m_s**2*(CS%f2_dx2_h(i,j) + CS%cg1(i,j)*CS%beta_dx2_h(i,j)))**power_2
      CS%Res_fn_h(i,j) = dx_term / &
          (dx_term + (CS%Res_coef_visc * US%L_T_to_m_s*CS%cg1(i,j))**CS%Res_fn_power_visc)
    enddo ; enddo
!$OMP do
    do J=js-1,Jeq ; do I=is-1,Ieq
      cg1_q = 0.25 * ((CS%cg1(i,j) + CS%cg1(i+1,j+1)) + (CS%cg1(i+1,j) + CS%cg1(i,j+1)))
      dx_term = (US%L_T_to_m_s**2*(CS%f2_dx2_q(I,J) + cg1_q * CS%beta_dx2_q(I,J)))**power_2
      CS%Res_fn_q(I,J) = dx_term / &
          (dx_term + (CS%Res_coef_visc * US%L_T_to_m_s*cg1_q)**CS%Res_fn_power_visc)
    enddo ; enddo
  else
!$OMP do
    do j=js-1,je+1 ; do i=is-1,ie+1
      dx_term = (US%L_T_to_m_s*sqrt(CS%f2_dx2_h(i,j) + &
                                    CS%cg1(i,j)*CS%beta_dx2_h(i,j)))**CS%Res_fn_power_visc
      CS%Res_fn_h(i,j) = dx_term / &
         (dx_term + (CS%Res_coef_visc * US%L_T_to_m_s*CS%cg1(i,j))**CS%Res_fn_power_visc)
    enddo ; enddo
!$OMP do
    do J=js-1,Jeq ; do I=is-1,Ieq
      cg1_q = 0.25 * ((CS%cg1(i,j) + CS%cg1(i+1,j+1)) + (CS%cg1(i+1,j) + CS%cg1(i,j+1)))
      dx_term = (US%L_T_to_m_s*sqrt(CS%f2_dx2_q(I,J) + &
                                    cg1_q * CS%beta_dx2_q(I,J)))**CS%Res_fn_power_visc
      CS%Res_fn_q(I,J) = dx_term / &
          (dx_term + (CS%Res_coef_visc * US%L_T_to_m_s*cg1_q)**CS%Res_fn_power_visc)
    enddo ; enddo
  endif

  if (CS%interpolate_Res_fn) then
    do j=js,je ; do I=is-1,Ieq
      CS%Res_fn_u(I,j) = 0.5*(CS%Res_fn_h(i,j) + CS%Res_fn_h(i+1,j))
    enddo ; enddo
    do J=js-1,Jeq ; do i=is,ie
      CS%Res_fn_v(i,J) = 0.5*(CS%Res_fn_h(i,j) + CS%Res_fn_h(i,j+1))
    enddo ; enddo
  else ! .not.CS%interpolate_Res_fn
    if (CS%Res_fn_power_khth >= 100) then
!$OMP do
      do j=js,je ; do I=is-1,Ieq
        cg1_u = 0.5 * (CS%cg1(i,j) + CS%cg1(i+1,j))
        dx_term = CS%f2_dx2_u(I,j) + cg1_u * CS%beta_dx2_u(I,j)
        if ((CS%Res_coef_khth * cg1_u)**2 > dx_term) then
          CS%Res_fn_u(I,j) = 0.0
        else
          CS%Res_fn_u(I,j) = 1.0
        endif
      enddo ; enddo
!$OMP do
      do J=js-1,Jeq ; do i=is,ie
        cg1_v = 0.5 * (CS%cg1(i,j) + CS%cg1(i,j+1))
        dx_term = CS%f2_dx2_v(i,J) + cg1_v * CS%beta_dx2_v(i,J)
        if ((CS%Res_coef_khth * cg1_v)**2 > dx_term) then
          CS%Res_fn_v(i,J) = 0.0
        else
          CS%Res_fn_v(i,J) = 1.0
        endif
      enddo ; enddo
    elseif (CS%Res_fn_power_khth == 2) then
!$OMP do
      do j=js,je ; do I=is-1,Ieq
        cg1_u = 0.5 * (CS%cg1(i,j) + CS%cg1(i+1,j))
        dx_term = CS%f2_dx2_u(I,j) + cg1_u * CS%beta_dx2_u(I,j)
        CS%Res_fn_u(I,j) = dx_term / (dx_term + (CS%Res_coef_khth * cg1_u)**2)
      enddo ; enddo
!$OMP do
      do J=js-1,Jeq ; do i=is,ie
        cg1_v = 0.5 * (CS%cg1(i,j) + CS%cg1(i,j+1))
        dx_term = CS%f2_dx2_v(i,J) + cg1_v * CS%beta_dx2_v(i,J)
        CS%Res_fn_v(i,J) = dx_term / (dx_term + (CS%Res_coef_khth * cg1_v)**2)
      enddo ; enddo
    elseif (mod(CS%Res_fn_power_khth, 2) == 0) then
      power_2 = CS%Res_fn_power_khth / 2
!$OMP do
      do j=js,je ; do I=is-1,Ieq
        cg1_u = 0.5 * (CS%cg1(i,j) + CS%cg1(i+1,j))
        dx_term = (US%L_T_to_m_s**2 * (CS%f2_dx2_u(I,j) + cg1_u * CS%beta_dx2_u(I,j)))**power_2
        CS%Res_fn_u(I,j) = dx_term / &
            (dx_term + (CS%Res_coef_khth * US%L_T_to_m_s*cg1_u)**CS%Res_fn_power_khth)
      enddo ; enddo
!$OMP do
      do J=js-1,Jeq ; do i=is,ie
        cg1_v = 0.5 * (CS%cg1(i,j) + CS%cg1(i,j+1))
        dx_term = (US%L_T_to_m_s**2 * (CS%f2_dx2_v(i,J) + cg1_v * CS%beta_dx2_v(i,J)))**power_2
        CS%Res_fn_v(i,J) = dx_term / &
            (dx_term + (CS%Res_coef_khth * US%L_T_to_m_s*cg1_v)**CS%Res_fn_power_khth)
      enddo ; enddo
    else
!$OMP do
      do j=js,je ; do I=is-1,Ieq
        cg1_u = 0.5 * (CS%cg1(i,j) + CS%cg1(i+1,j))
        dx_term = (US%L_T_to_m_s*sqrt(CS%f2_dx2_u(I,j) + &
                                      cg1_u * CS%beta_dx2_u(I,j)))**CS%Res_fn_power_khth
        CS%Res_fn_u(I,j) = dx_term / &
            (dx_term + (CS%Res_coef_khth * US%L_T_to_m_s*cg1_u)**CS%Res_fn_power_khth)
      enddo ; enddo
!$OMP do
      do J=js-1,Jeq ; do i=is,ie
        cg1_v = 0.5 * (CS%cg1(i,j) + CS%cg1(i,j+1))
        dx_term = (US%L_T_to_m_s*sqrt(CS%f2_dx2_v(i,J) + &
                                      cg1_v * CS%beta_dx2_v(i,J)))**CS%Res_fn_power_khth
        CS%Res_fn_v(i,J) = dx_term / &
            (dx_term + (CS%Res_coef_khth * US%L_T_to_m_s*cg1_v)**CS%Res_fn_power_khth)
      enddo ; enddo
    endif
  endif
!$OMP end parallel

  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_Res_fn > 0) call post_data(CS%id_Res_fn, CS%Res_fn_h, CS%diag)
    if (CS%id_BS_struct > 0) call post_data(CS%id_BS_struct, CS%BS_struct, CS%diag)
    if (CS%id_khth_struct > 0) call post_data(CS%id_khth_struct, CS%khth_struct, CS%diag)
    if (CS%id_khtr_struct > 0) call post_data(CS%id_khtr_struct, CS%khtr_struct, CS%diag)
    if (CS%id_kdgl90_struct > 0) call post_data(CS%id_kdgl90_struct, CS%kdgl90_struct, CS%diag)
  endif

  if (CS%debug) then
    call hchksum(CS%cg1, "calc_resoln_fn cg1", G%HI, haloshift=1, unscale=US%L_T_to_m_s)
    call uvchksum("Res_fn_[uv]", CS%Res_fn_u, CS%Res_fn_v, G%HI, haloshift=0, &
                  unscale=1.0, scalar_pair=.true.)
  endif

end subroutine calc_resoln_function

!> Calculates and stores functions of SQG mode
subroutine calc_sqg_struct(h, tv, G, GV, US, CS, dt, MEKE, OBC)
  type(ocean_grid_type),                     intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(in)    :: tv !<Thermodynamic variables
   real,                                     intent(in)    :: dt !< Time increment [T ~> s]
  type(VarMix_CS),                           intent(inout) :: CS !< Variable mixing control struct
  type(MEKE_type),                           intent(in)    :: MEKE !< MEKE struct
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure

  ! Local variables
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)+1) :: e    ! The interface heights relative to mean sea level [Z ~> m]
  real, dimension(SZIB_(G), SZJ_(G),SZK_(GV)+1) :: N2_u ! Square of buoyancy frequency at u-points [L2 Z-2 T-2 ~> s-2]
  real, dimension(SZI_(G), SZJB_(G),SZK_(GV)+1) :: N2_v ! Square of buoyancy frequency at v-points [L2 Z-2 T-2 ~> s-2]
  real, dimension(SZIB_(G), SZJ_(G),SZK_(GV)+1) :: dzu ! Z-thickness at u-points [Z ~> m]
  real, dimension(SZI_(G), SZJB_(G),SZK_(GV)+1) :: dzv ! Z-thickness at v-points [Z ~> m]
  real, dimension(SZIB_(G), SZJ_(G),SZK_(GV)+1) :: dzSxN ! |Sx| N times dz at u-points [Z T-1 ~> m s-1]
  real, dimension(SZI_(G), SZJB_(G),SZK_(GV)+1) :: dzSyN ! |Sy| N times dz at v-points [Z T-1 ~> m s-1]
  real, dimension(SZI_(G), SZJ_(G)) :: f  ! Absolute value of the Coriolis parameter at h point [T-1 ~> s-1]
  real :: N2             ! Positive buoyancy frequency square or zero [L2 Z-2 T-2 ~> s-2]
  real :: dzc            ! Spacing between two adjacent layers in stretched vertical coordinate [Z ~> m]
  real :: f_subround     ! The minimal resolved value of Coriolis parameter to prevent division by zero [T-1 ~> s-1]
  real, dimension(SZI_(G), SZJ_(G)) :: Le  ! Eddy length scale [L ~> m]
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  f_subround = 1.0e-40 * US%s_to_T

  if (.not. CS%initialized) call MOM_error(FATAL, "MOM_lateral_mixing_coeffs.F90, calc_slope_functions: "//&
         "Module must be initialized before it is used.")

  call find_eta(h, tv, G, GV, US, e, halo_size=2)
  call calc_isoneutral_slopes(G, GV, US, h, e, tv, dt*CS%kappa_smooth, CS%use_stanley_iso, &
                                  CS%slope_x, CS%slope_y, N2_u=N2_u, N2_v=N2_v, dzu=dzu, dzv=dzv, &
                                  dzSxN=dzSxN, dzSyN=dzSyN, halo=1, OBC=OBC, OBC_N2=CS%OBC_friendly)

  if (CS%sqg_expo<=0.) then
    CS%sqg_struct(:,:,:) = 1.
  else
    do j=js,je ; do i=is,ie
      CS%sqg_struct(i,j,1) = 1.0
    enddo ; enddo
    if (allocated(MEKE%Le)) then
      do j=js,je ; do i=is,ie
        Le(i,j) = MEKE%Le(i,j)
        f(i,j) = max(0.25 * abs((G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J-1)) + &
                         (G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J-1))), f_subround)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        Le(i,j) = sqrt(G%areaT(i,j))
        f(i,j) = max(0.25 * abs((G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J-1)) + &
                         (G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J-1))), f_subround)
      enddo ; enddo
    endif
    do k=2,nz ; do j=js,je ; do i=is,ie
      N2 = max(0.25 * ((N2_u(I-1,j,k) + N2_u(I,j,k)) + (N2_v(i,J-1,k) + N2_v(i,J,k))), 0.0)
      dzc = 0.25 * ((dzu(I-1,j,k) + dzu(I,j,k)) + (dzv(i,J-1,k) + dzv(i,J,k)))
      CS%sqg_struct(i,j,k) = CS%sqg_struct(i,j,k-1) * &
              exp(-CS%sqg_expo * (dzc * sqrt(N2)/(f(i,j) * Le(i,j))))
    enddo ; enddo ; enddo
  endif


  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_sqg_struct > 0) call post_data(CS%id_sqg_struct, CS%sqg_struct, CS%diag)
    if (CS%id_N2_u > 0) call post_data(CS%id_N2_u, N2_u, CS%diag)
    if (CS%id_N2_v > 0) call post_data(CS%id_N2_v, N2_v, CS%diag)
  endif

end subroutine calc_sqg_struct

!> Calculates and stores functions of isopycnal slopes, e.g. Sx, Sy, S*N, mostly used in the Visbeck et al.
!! style scaling of diffusivity
subroutine calc_slope_functions(h, tv, dt, G, GV, US, CS, OBC)
  type(ocean_grid_type),                     intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(in)    :: tv !< Thermodynamic variables
  real,                                      intent(in)    :: dt !< Time increment [T ~> s]
  type(VarMix_CS),                           intent(inout) :: CS !< Variable mixing control structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1)  :: e    ! The interface heights relative to mean sea level [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: N2_u ! Square of buoyancy frequency at u-points [L2 Z-2 T-2 ~> s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: N2_v ! Square of buoyancy frequency at v-points [L2 Z-2 T-2 ~> s-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: dzu  ! Z-thickness at u-points [Z ~> m]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: dzv  ! Z-thickness at v-points [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: dzSxN ! |Sx| N times dz at u-points [Z T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: dzSyN ! |Sy| N times dz at v-points [Z T-1 ~> m s-1]

  if (.not. CS%initialized) call MOM_error(FATAL, "MOM_lateral_mixing_coeffs.F90, calc_slope_functions: "//&
         "Module must be initialized before it is used.")

  if (CS%calculate_Eady_growth_rate) then
    call find_eta(h, tv, G, GV, US, e, halo_size=2)
    if (CS%use_simpler_Eady_growth_rate) then
      call calc_isoneutral_slopes(G, GV, US, h, e, tv, dt*CS%kappa_smooth, CS%use_stanley_iso, &
                                  CS%slope_x, CS%slope_y, N2_u=N2_u, N2_v=N2_v, dzu=dzu, dzv=dzv, &
                                  dzSxN=dzSxN, dzSyN=dzSyN, halo=1, OBC=OBC, OBC_N2=CS%OBC_friendly)
      call calc_Eady_growth_rate_2D(CS, G, GV, US, h, e, dzu, dzv, dzSxN, dzSyN, CS%SN_u, CS%SN_v)
    elseif (CS%use_stored_slopes) then
      call calc_isoneutral_slopes(G, GV, US, h, e, tv, dt*CS%kappa_smooth, CS%use_stanley_iso, &
                                  CS%slope_x, CS%slope_y, N2_u=N2_u, N2_v=N2_v, halo=1, OBC=OBC, &
                                  OBC_N2=CS%OBC_friendly)
      call calc_Visbeck_coeffs_old(h, CS%slope_x, CS%slope_y, N2_u, N2_v, G, GV, US, CS, OBC)
    else
      call calc_slope_functions_using_just_e(h, G, GV, US, CS, e)
    endif
  endif

  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_dzu > 0) call post_data(CS%id_dzu, dzu, CS%diag)
    if (CS%id_dzv > 0) call post_data(CS%id_dzv, dzv, CS%diag)
    if (CS%id_dzSxN > 0) call post_data(CS%id_dzSxN, dzSxN, CS%diag)
    if (CS%id_dzSyN > 0) call post_data(CS%id_dzSyN, dzSyN, CS%diag)
    if (CS%id_SN_u > 0) call post_data(CS%id_SN_u, CS%SN_u, CS%diag)
    if (CS%id_SN_v > 0) call post_data(CS%id_SN_v, CS%SN_v, CS%diag)
    if (CS%id_L2u > 0)  call post_data(CS%id_L2u, CS%L2u, CS%diag)
    if (CS%id_L2v > 0)  call post_data(CS%id_L2v, CS%L2v, CS%diag)
    if (CS%id_N2_u > 0) call post_data(CS%id_N2_u, N2_u, CS%diag)
    if (CS%id_N2_v > 0) call post_data(CS%id_N2_v, N2_v, CS%diag)
  endif

end subroutine calc_slope_functions

!> Calculates factors used when setting diffusivity coefficients similar to Visbeck et al., 1997.
!! This is on older implementation that is susceptible to large values of Eady growth rate
!! for incropping layers.
subroutine calc_Visbeck_coeffs_old(h, slope_x, slope_y, N2_u, N2_v, G, GV, US, CS, OBC)
  type(ocean_grid_type),                        intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: slope_x !< Zonal isoneutral slope [Z L-1 ~> nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: N2_u    !< Buoyancy (Brunt-Vaisala) frequency
                                                                         !! at u-points [L2 Z-2 T-2 ~> s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in)    :: slope_y !< Meridional isoneutral slope
                                                                         !! [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in)    :: N2_v    !< Buoyancy (Brunt-Vaisala) frequency
                                                                         !! at v-points [L2 Z-2 T-2 ~> s-2]
  type(unit_scale_type),                        intent(in)    :: US !< A dimensional unit scaling type
  type(VarMix_CS),                              intent(inout) :: CS !< Variable mixing control structure
  type(ocean_OBC_type),                         pointer       :: OBC  !< Open boundaries control structure.

  ! Local variables
  real :: S2            ! Interface slope squared [Z2 L-2 ~> nondim]
  real :: N2            ! Positive buoyancy frequency or zero [L2 Z-2 T-2 ~> s-2]
  real :: Hup, Hdn      ! Thickness from above, below [H ~> m or kg m-2]
  real :: H_geom        ! The geometric mean of Hup and Hdn [H ~> m or kg m-2].
  real :: S2max         ! An upper bound on the squared slopes [Z2 L-2 ~> nondim]
  real :: wNE, wSE, wSW, wNW ! Weights of adjacent points [nondim]
  real :: H_u(SZIB_(G)), H_v(SZI_(G)) ! Layer thicknesses at u- and v-points [H ~> m or kg m-2]

  ! Note that at some points in the code S2_u and S2_v hold the running depth
  ! integrals of the squared slope [H ~> m or kg m-2] before the average is taken.
  real :: S2_u(SZIB_(G),SZJ_(G)) ! At first the thickness-weighted depth integral of the squared
                                 ! slope [H Z2 L-2 ~> m or kg m-2] and then the average of the
                                 ! squared slope [Z2 L-2 ~> nondim] at u points.
  real :: S2_v(SZI_(G),SZJB_(G)) ! At first the thickness-weighted depth integral of the squared
                                 ! slope [H Z2 L-2 ~> m or kg m-2] and then the average of the
                                 ! squared slope [Z2 L-2 ~> nondim] at v points.
  integer :: OBC_dir_u(SZIB_(G),SZJ_(G))  ! An integer indicating where there are u OBCs: +1 for
                                 ! eastern OBCs, -1 for western OBCs and 0 at points with no OBCs.
  integer :: OBC_dir_v(SZI_(G),SZJB_(G))  ! An integer indicating where there are v OBCs: +1 for
                                 ! northern OBCs, -1 for southern OBCs and 0 at points with no OBCs.
  real :: h4_u(SZIB_(G),SZJ_(G),SZK_(GV)+1)  ! The product of the 4 thicknesses surrounding a u-point
                                 ! interface or the inward equivalent with OBCs [H4 ~> m4 or kg2 m-4]
  real :: h4_v(SZI_(G),SZJB_(G),SZK_(GV)+1)  ! The product of the 4 thicknesses surrounding a v-point
                                 ! interface or the inward equivalent with OBCs [H4 ~> m4 or kg2 m-4]
  integer :: i, j, k, is, ie, js, je, nz

  if (.not. CS%initialized) call MOM_error(FATAL, "calc_Visbeck_coeffs_old: "// &
         "Module must be initialized before it is used.")

  if (.not. CS%calculate_Eady_growth_rate) return
  if (.not. allocated(CS%SN_u)) call MOM_error(FATAL, "calc_slope_function:"// &
         "%SN_u is not associated with use_variable_mixing.")
  if (.not. allocated(CS%SN_v)) call MOM_error(FATAL, "calc_slope_function:R"// &
         "%SN_v is not associated with use_variable_mixing.")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  S2max = CS%Visbeck_S_max**2

  CS%SN_u(:,:) = 0.0
  CS%SN_v(:,:) = 0.0

  ! These settings apply where there are not open boundary conditions.
  OBC_dir_u(:,:) = 0 ; OBC_dir_v(:,:) = 0

  if (associated(OBC).and. CS%OBC_friendly) then
   ! Store the direction of any OBC faces.
   !$OMP parallel do default(shared)
    do j=js-1,je+1 ; do I=is-1,ie ; if (OBC%segnum_u(I,j) /= 0) then
      if (OBC%segnum_u(I,j) > 0) OBC_dir_u(I,j) = 1   !  OBC_DIRECTION_E
      if (OBC%segnum_u(I,j) < 0) OBC_dir_u(I,j) = -1  !  OBC_DIRECTION_W
    endif ; enddo ; enddo
   !$OMP parallel do default(shared)
    do J=js-1,je ; do i=is-1,ie+1 ; if (OBC%segnum_v(i,J) /= 0) then
      if (OBC%segnum_v(i,J) > 0) OBC_dir_v(i,J) = 1   ! OBC_DIRECTION_N
      if (OBC%segnum_v(i,J) < 0) OBC_dir_v(i,J) = -1  !  OBC_DIRECTION_S
    endif ; enddo ; enddo

    ! Use the masked product of the 4 (or 2) thicknesses around a velocity-point interface for weights.
    !$OMP parallel do default(shared)
    do K=2,nz
      do j=js-1,je+1 ; do I=is-1,ie
        if (OBC_dir_u(I,j) == 0) then
          h4_u(I,j,K) = G%mask2dCu(I,j) * ( (h(i,j,k)*h(i+1,j,k)) * (h(i,j,k-1)*h(i+1,j,k-1)) )
        elseif (OBC_dir_u(I,j) == 1) then  ! OBC_DIRECTION_E
          h4_u(I,j,K) = G%mask2dCu(I,j) * ( (h(i,j,k)**2) * (h(i,j,k-1)**2) )
        elseif (OBC_dir_u(I,j) == -1) then  ! OBC_DIRECTION_W
          h4_u(I,j,K) = G%mask2dCu(I,j) * ( (h(i+1,j,k)**2) * (h(i+1,j,k-1)**2) )
        endif
      enddo ; enddo
      do J=js-1,je ; do i=is-1,ie+1
        if (OBC_dir_v(i,J) == 0) then
          h4_v(i,J,K) = G%mask2dCv(i,J) * ( (h(i,j,k)*h(i,j+1,k)) * (h(i,j,k-1)*h(i,j+1,k-1)) )
        elseif (OBC_dir_v(i,J) == 1) then  ! OBC_DIRECTION_N
          h4_v(i,J,K) = G%mask2dCv(i,J) * ( (h(i,j,k)**2) * (h(i,j,k-1)**2) )
        elseif (OBC_dir_v(i,J) == -1) then  ! OBC_DIRECTION_S
          h4_v(i,J,K) = G%mask2dCv(i,J) * ( (h(i,j+1,k)**2) * (h(i,j+1,k-1)**2) )
        endif
      enddo ; enddo
    enddo
  else  ! The land mask is sufficient and there are no special considerations taken at OBC points.
    ! Use the masked product of the 4 thicknesses around a velocity-point interface for weights.
    !$OMP parallel do default(shared)
    do K=2,nz
      do j=js-1,je+1 ; do I=is-1,ie
        h4_u(I,j,K) = G%mask2dCu(I,j) * ( (h(i,j,k)*h(i+1,j,k)) * (h(i,j,k-1)*h(i+1,j,k-1)) )
      enddo ; enddo
      do J=js-1,je ; do i=is-1,ie+1
        h4_v(i,J,K) = G%mask2dCv(i,J) * ( (h(i,j,k)*h(i,j+1,k)) * (h(i,j,k-1)*h(i,j+1,k-1)) )
      enddo ; enddo
    enddo
  endif

  ! To set the length scale based on the deformation radius, use wave_speed to
  ! calculate the first-mode gravity wave speed and then blend the equatorial
  ! and midlatitude deformation radii, using calc_resoln_function as a template.

  !$OMP parallel do default(shared) private(S2,H_u,Hdn,Hup,H_geom,N2,wNE,wSE,wSW,wNW)
  do j=js,je
    do I=is-1,ie
      CS%SN_u(I,j) = 0. ; H_u(I) = 0. ; S2_u(I,j) = 0.
    enddo
    do K=2,nz ; do I=is-1,ie
      Hdn = sqrt( h(i,j,k) * h(i+1,j,k) )
      Hup = sqrt( h(i,j,k-1) * h(i+1,j,k-1) )
      H_geom = sqrt( Hdn * Hup )
     !H_geom = H_geom * sqrt(N2) ! WKB-ish
     !H_geom = H_geom * N2       ! WKB-ish
      wSE = h4_v(i+1,J-1,K)
      wNW = h4_v(i,J,K)
      wNE = h4_v(i+1,J,K)
      wSW = h4_v(i,J-1,K)
      if (OBC_dir_u(I,j) == 1) then  ! OBC_DIRECTION_E
        wSE = 0.0 ; wNE = 0.0
        H_geom = sqrt( h(i,j,k) * h(i,j,k-1) )
      elseif (OBC_dir_u(I,j) == -1) then  ! OBC_DIRECTION_W
        wSW = 0.0 ; wNW = 0.0
        H_geom = sqrt( h(i+1,j,k) * h(i+1,j,k-1) )
      endif
      S2 =  slope_x(I,j,K)**2 + &
              (((wNW*slope_y(i,J,K)**2) + (wSE*slope_y(i+1,J-1,K)**2)) + &
               ((wNE*slope_y(i+1,J,K)**2) + (wSW*slope_y(i,J-1,K)**2)) ) / &
              ( ((wSE+wNW) + (wNE+wSW)) + GV%H_subroundoff**4 )
      if (S2max>0.) S2 = S2 * S2max / (S2 + S2max) ! Limit S2

      N2 = max(0., N2_u(I,j,k))
      CS%SN_u(I,j) = CS%SN_u(I,j) + sqrt( S2*N2 )*H_geom
      S2_u(I,j) = S2_u(I,j) + S2*H_geom
      H_u(I) = H_u(I) + H_geom
    enddo ; enddo
    do I=is-1,ie
      if (H_u(I)>0.) then
        CS%SN_u(I,j) = G%OBCmaskCu(I,j) * CS%SN_u(I,j) / H_u(I)
        S2_u(I,j) =  G%OBCmaskCu(I,j) * S2_u(I,j) / H_u(I)
      else
        CS%SN_u(I,j) = 0.
      endif
    enddo
  enddo

  !$OMP parallel do default(shared) private(S2,H_v,Hdn,Hup,H_geom,N2,wNE,wSE,wSW,wNW)
  do J=js-1,je
    do i=is,ie
      CS%SN_v(i,J) = 0. ; H_v(i) = 0. ; S2_v(i,J) = 0.
    enddo
    do K=2,nz ; do i=is,ie
      Hdn = sqrt( h(i,j,k) * h(i,j+1,k) )
      Hup = sqrt( h(i,j,k-1) * h(i,j+1,k-1) )
      H_geom = sqrt( Hdn * Hup )
     !H_geom = H_geom * sqrt(N2) ! WKB-ish
     !H_geom = H_geom * N2       ! WKB-ish
      wSE = h4_u(I,j,K)
      wNW = h4_u(I-1,j+1,K)
      wNE = h4_u(I,j+1,K)
      wSW = h4_u(I-1,j,K)
      if (OBC_dir_v(i,J) == 1) then  ! OBC_DIRECTION_N
        wNW = 0.0 ; wNE = 0.0
        H_geom = sqrt( h(i,j,k) *  h(i,j,k-1) )
      elseif (OBC_dir_v(i,J) == -1) then  ! OBC_DIRECTION_S
        wSW = 0.0 ; wSE = 0.0
        H_geom = sqrt( h(i,j+1,k) * h(i,j+1,k-1) )
      endif
      S2 = slope_y(i,J,K)**2 + &
             (((wSE*slope_x(I,j,K)**2) + (wNW*slope_x(I-1,j+1,K)**2)) + &
              ((wNE*slope_x(I,j+1,K)**2) + (wSW*slope_x(I-1,j,K)**2)) ) / &
             ( ((wSE+wNW) + (wNE+wSW)) + GV%H_subroundoff**4 )
      if (S2max>0.) S2 = S2 * S2max / (S2 + S2max) ! Limit S2

      N2 = max(0., N2_v(i,J,K))
      CS%SN_v(i,J) = CS%SN_v(i,J) + sqrt( S2*N2 )*H_geom
      S2_v(i,J) = S2_v(i,J) + S2*H_geom
      H_v(i) = H_v(i) + H_geom
    enddo ; enddo
    do i=is,ie
      if (H_v(i)>0.) then
        CS%SN_v(i,J) = G%OBCmaskCv(i,J) * CS%SN_v(i,J) / H_v(i)
        S2_v(i,J) = G%OBCmaskCv(i,J) * S2_v(i,J) / H_v(i)
      else
        CS%SN_v(i,J) = 0.
      endif
    enddo
  enddo

  ! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_S2_u > 0) call post_data(CS%id_S2_u, S2_u, CS%diag)
    if (CS%id_S2_v > 0) call post_data(CS%id_S2_v, S2_v, CS%diag)
  endif

  if (CS%debug) then
    call uvchksum("calc_Visbeck_coeffs_old slope_[xy]", slope_x, slope_y, G%HI, &
                  unscale=US%Z_to_L, haloshift=1)
    ! call uvchksum("calc_Visbeck_coeffs_old S2_[uv]", S2_u, S2_v, G%HI, &
    !               unscale=US%Z_to_L**2, scalar_pair=.true.)
    call uvchksum("calc_Visbeck_coeffs_old N2_u, N2_v", N2_u, N2_v, G%HI, &
                  unscale=US%L_to_Z**2*US%s_to_T**2, scalar_pair=.true.)
    call uvchksum("calc_Visbeck_coeffs_old SN_[uv]", CS%SN_u, CS%SN_v, G%HI, &
                  unscale=US%s_to_T, scalar_pair=.true.)
  endif

end subroutine calc_Visbeck_coeffs_old

!> Calculates the Eady growth rate (2D fields) for use in MEKE and the Visbeck schemes
subroutine calc_Eady_growth_rate_2D(CS, G, GV, US, h, e, dzu, dzv, dzSxN, dzSyN, SN_u, SN_v)
  type(VarMix_CS),                              intent(inout) :: CS !< Variable mixing coefficients
  type(ocean_grid_type),                        intent(in) :: G   !< Ocean grid structure
  type(verticalGrid_type),                      intent(in) :: GV  !< Vertical grid structure
  type(unit_scale_type),                        intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: h   !< Interface height [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(in) :: e   !< Interface height [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: dzu !< dz at u-points [Z ~> m]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in) :: dzv !< dz at v-points [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: dzSxN !< dz Sx N at u-points [Z T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in) :: dzSyN !< dz Sy N at v-points [Z T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: SN_u !< SN at u-points [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: SN_v !< SN at v-points [T-1 ~> s-1]
  ! Local variables
  real :: D_scale ! The depth over which to average SN [Z ~> m]
  real :: dnew ! Depth of bottom of layer [Z ~> m]
  real :: dz ! Limited thickness of this layer [Z ~> m]
  real :: weight ! Fraction of this layer that contributes to integral [nondim]
  real :: sum_dz(SZI_(G)) ! Cumulative sum of z-thicknesses [Z ~> m]
  real :: vint_SN(SZIB_(G)) ! Cumulative integral of SN [Z T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G)) :: SN_cpy !< SN at u-points [T-1 ~> s-1]
  real :: dz_neglect ! A negligibly small distance to avoid division by zero [Z ~> m]
  real :: r_crp_dist ! The inverse of the distance over which to scale the cropping [Z-1 ~> m-1]
  real :: dB, dT ! Elevation variables used when cropping [Z ~> m]
  integer :: i, j, k, l_seg
  logical :: crop

  dz_neglect = GV%dZ_subroundoff
  D_scale = CS%Eady_GR_D_scale
  if (D_scale<=0.) D_scale = 64.*GV%max_depth ! 0 means use full depth so choose something big
  r_crp_dist = 1. / max( dz_neglect, CS%cropping_distance )
  crop = CS%cropping_distance>=0. ! Only filter out in-/out-cropped interface is parameter if non-negative

  if (CS%debug) then
    call uvchksum("calc_Eady_growth_rate_2D dz[uv]", dzu, dzv, G%HI, unscale=US%Z_to_m, scalar_pair=.true.)
    call uvchksum("calc_Eady_growth_rate_2D dzS2N2[uv]", dzSxN, dzSyN, G%HI, &
                  unscale=US%Z_to_m*US%s_to_T, scalar_pair=.true.)
  endif

  !$OMP parallel do default(shared)
  do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
    CS%SN_u(i,j) = 0.0
    CS%SN_v(i,j) = 0.0
  enddo ; enddo

  !$OMP parallel do default(shared) private(dnew,dz,weight,l_seg,vint_SN,sum_dz,dT,dB)
  do j=G%jsc-1,G%jec+1
    do I=G%isc-1,G%iec
      vint_SN(I) = 0.
      sum_dz(I) = dz_neglect
    enddo
    if (crop) then
      do K=2,GV%ke ; do I=G%isc-1,G%iec
        dnew = sum_dz(I) + dzu(I,j,K) ! This is where the bottom of the layer is
        dnew = min(dnew, D_scale) ! This limits the depth to D_scale
        dz = max(0., dnew - sum_dz(I)) ! This is the part of the layer to be included in the integral.
                                       ! When D_scale>dnew, dz=dzu (+roundoff error).
                                       ! When sum_dz<D_scale<dnew, 0<dz<dzu.
                                       ! When D_scale<sum_dz, dz=0.
        weight = dz / ( dzu(I,j,K) + dz_neglect ) ! Fraction of this layer to include
        dT = min( e(i,j,1), e(i+1,j,1) ) ! Deepest sea surface
        dB = max( e(i,j,K), e(i+1,j,K) ) ! Shallowest interface
        weight = weight * min( max( 0., (dT-dB)*r_crp_dist ), 1. )
        dT = min( e(i,j,K), e(i+1,j,K) ) ! Deepest interface
        dB = max( e(i,j,GV%ke+1), e(i+1,j,GV%ke+1) ) ! Shallowest topography
        weight = weight * min( max( 0., (dT-dB)*r_crp_dist ), 1. )
        vint_SN(I) = vint_SN(I) + weight * dzSxN(I,j,K)
        sum_dz(I) = sum_dz(I) + weight * dzu(I,j,K)
      enddo ; enddo
    else
      do K=2,GV%ke ; do I=G%isc-1,G%iec
        dnew = sum_dz(I) + dzu(I,j,K) ! This is where the bottom of the layer is
        dnew = min(dnew, D_scale) ! This limits the depth to D_scale
        dz = max(0., dnew - sum_dz(I)) ! This is the part of the layer to be included in the integral.
                                       ! When D_scale>dnew, dz=dzu (+roundoff error).
                                       ! When sum_dz<D_scale<dnew, 0<dz<dzu.
                                       ! When D_scale<sum_dz, dz=0.
        weight = dz / ( dzu(I,j,K) + dz_neglect ) ! Fraction of this layer to include
        vint_SN(I) = vint_SN(I) + weight * dzSxN(I,j,K)
        sum_dz(I) = sum_dz(I) + weight * dzu(I,j,K)
      enddo ; enddo
    endif
    do I=G%isc-1,G%iec
      CS%SN_u(I,j) = G%OBCmaskCu(I,j) * ( vint_SN(I) / sum_dz(I) )
      SN_cpy(I,j) = G%OBCmaskCu(I,j) * ( vint_SN(I) / sum_dz(I) )
    enddo
  enddo

  !$OMP parallel do default(shared) private(dnew,dz,weight,l_seg,vint_SN,sum_dz,dT,dB)
  do J=G%jsc-1,G%jec
    do i=G%isc-1,G%iec+1
      vint_SN(i) = 0.
      sum_dz(i) = dz_neglect
    enddo
    if (crop) then
      do K=2,GV%ke ; do i=G%isc-1,G%iec+1
        dnew = sum_dz(i) + dzv(i,J,K) ! This is where the bottom of the layer is
        dnew = min(dnew, D_scale) ! This limits the depth to D_scale
        dz = max(0., dnew - sum_dz(i)) ! This is the part of the layer to be included in the integral.
                                       ! When D_scale>dnew, dz=dzu (+roundoff error).
                                       ! When sum_dz<D_scale<dnew, 0<dz<dzu.
                                       ! When D_scale<sum_dz, dz=0.
        weight = dz / ( dzv(i,J,K) + dz_neglect ) ! Fraction of this layer to include
        dT = min( e(i,j,1), e(i,j+1,1) ) ! Deepest sea surface
        dB = max( e(i,j,K), e(i,j+1,K) ) ! Shallowest interface
        weight = weight * min( max( 0., (dT-dB)*r_crp_dist ), 1. )
        dT = min( e(i,j,K), e(i,j+1,K) )! Deepest interface
        dB = max( e(i,j,GV%ke+1), e(i,j+1,GV%ke+1) ) ! Shallowest topography
        weight = weight * min( max( 0., (dT-dB)*r_crp_dist ), 1. )
        vint_SN(I) = vint_SN(I) + weight**2 * dzSyN(i,J,K)
        sum_dz(i) = sum_dz(i) + weight * dzv(i,J,K)
      enddo ; enddo
    else
      do K=2,GV%ke ; do i=G%isc-1,G%iec+1
        dnew = sum_dz(i) + dzv(i,J,K) ! This is where the bottom of the layer is
        dnew = min(dnew, D_scale) ! This limits the depth to D_scale
        dz = max(0., dnew - sum_dz(i)) ! This is the part of the layer to be included in the integral.
                                       ! When D_scale>dnew, dz=dzu (+roundoff error).
                                       ! When sum_dz<D_scale<dnew, 0<dz<dzu.
                                       ! When D_scale<sum_dz, dz=0.
        weight = dz / ( dzv(i,J,K) + dz_neglect ) ! Fraction of this layer to include
        vint_SN(I) = vint_SN(I) + weight**2 * dzSyN(i,J,K)
        sum_dz(i) = sum_dz(i) + weight * dzv(i,J,K)
      enddo ; enddo
    endif
    do i=G%isc-1,G%iec+1
      CS%SN_v(i,J) = G%OBCmaskCv(i,J) * ( vint_SN(i) / sum_dz(i) )
    enddo
  enddo

  do j=G%jsc,G%jec
    do I=G%isc-1,G%iec
      CS%SN_u(I,j) = sqrt( SN_cpy(I,j)**2 &
                         + 0.25*( ((CS%SN_v(i,J)**2) + (CS%SN_v(i+1,J-1)**2)) &
                                + ((CS%SN_v(i+1,J)**2) + (CS%SN_v(i,J-1)**2)) ) )
    enddo
  enddo
  do J=G%jsc-1,G%jec
    do i=G%isc,G%iec
      CS%SN_v(i,J) = sqrt( CS%SN_v(i,J)**2 &
                         + 0.25*( ((SN_cpy(I,j)**2) + (SN_cpy(I-1,j+1)**2)) &
                                + ((SN_cpy(I,j+1)**2) + (SN_cpy(I-1,j)**2)) ) )
    enddo
  enddo

  if (CS%debug) then
    call uvchksum("calc_Eady_growth_rate_2D SN_[uv]", CS%SN_u, CS%SN_v, G%HI, &
                  unscale=US%s_to_T, scalar_pair=.true.)
  endif

end subroutine calc_Eady_growth_rate_2D

!> The original calc_slope_function() that calculated slopes using
!! interface positions only, not accounting for density variations.
subroutine calc_slope_functions_using_just_e(h, G, GV, US, CS, e)
  type(ocean_grid_type),                       intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                     intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  type(unit_scale_type),                       intent(in)    :: US !< A dimensional unit scaling type
  type(VarMix_CS),                             intent(inout) :: CS !< Variable mixing control structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: e  !< Interface position [Z ~> m]
  ! type(thermo_var_ptrs),                     intent(in)    :: tv !< Thermodynamic variables
  ! Local variables
  real :: E_x(SZIB_(G),SZJ_(G))  ! X-slope of interface at u points [Z L-1 ~> nondim] (for diagnostics)
  real :: E_y(SZI_(G),SZJB_(G))  ! Y-slope of interface at v points [Z L-1 ~> nondim] (for diagnostics)
  real :: dz_tot(SZI_(G),SZJ_(G)) ! The total thickness of the water columns [Z ~> m]
  ! real :: dz(SZI_(G),SZJ_(G),SZK_(GV)) ! The vertical distance across each layer [Z ~> m]
  real :: H_cutoff      ! Local estimate of a minimum thickness for masking [H ~> m or kg m-2]
  real :: dZ_cutoff     ! A minimum water column depth for masking [H ~> m or kg m-2]
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: S2            ! Interface slope squared [Z2 L-2 ~> nondim]
  real :: N2            ! Brunt-Vaisala frequency squared [L2 Z-2 T-2 ~> s-2]
  real :: Hup, Hdn      ! Thickness from above, below [H ~> m or kg m-2]
  real :: H_geom        ! The geometric mean of Hup*Hdn [H ~> m or kg m-2].
  real :: S2N2_u_local(SZIB_(G),SZJ_(G),SZK_(GV)) ! The depth integral of the slope times
                        ! the buoyancy frequency squared at u-points [Z T-2 ~> m s-2]
  real :: S2N2_v_local(SZI_(G),SZJB_(G),SZK_(GV)) ! The depth integral of the slope times
                        ! the buoyancy frequency squared at v-points [Z T-2 ~> m s-2]
  logical :: use_dztot  ! If true, use the total water column thickness rather than the
                        ! bathymetric depth for certain calculations.
  integer :: is, ie, js, je, nz
  integer :: i, j, k
  integer :: l_seg

  if (.not. CS%initialized) call MOM_error(FATAL, "calc_slope_functions_using_just_e: "// &
         "Module must be initialized before it is used.")

  if (.not. CS%calculate_Eady_growth_rate) return
  if (.not. allocated(CS%SN_u)) call MOM_error(FATAL, "calc_slope_function:"// &
         "%SN_u is not associated with use_variable_mixing.")
  if (.not. allocated(CS%SN_v)) call MOM_error(FATAL, "calc_slope_function:"// &
         "%SN_v is not associated with use_variable_mixing.")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  h_neglect = GV%H_subroundoff
  H_cutoff = real(2*nz) * (GV%Angstrom_H + h_neglect)
  dZ_cutoff = real(2*nz) * (GV%Angstrom_Z + GV%dz_subroundoff)

  use_dztot = CS%full_depth_Eady_growth_rate ! .or. .not.(GV%Boussinesq or GV%semi_Boussinesq)

  if (use_dztot) then
    !$OMP parallel do default(shared)
    do j=js-1,je+1 ; do i=is-1,ie+1
      dz_tot(i,j) = e(i,j,1) - e(i,j,nz+1)
    enddo ; enddo
    ! The following mathematically equivalent expression is more expensive but is less
    ! sensitive to roundoff for large Z_ref:
    ! call thickness_to_dz(h, tv, dz, G, GV, US, halo_size=1)
    ! do j=js-1,je+1
    !   do i=is-1,ie+1 ; dz_tot(i,j) = 0.0 ; enddo
    !   do k=1,nz ; do i=is-1,ie+1
    !     dz_tot(i,j) = dz_tot(i,j) + dz(i,j,k)
    !   enddo ; enddo
    ! enddo
  endif

  ! To set the length scale based on the deformation radius, use wave_speed to
  ! calculate the first-mode gravity wave speed and then blend the equatorial
  ! and midlatitude deformation radii, using calc_resoln_function as a template.

  !$OMP parallel do default(shared) private(E_x,E_y,S2,Hdn,Hup,H_geom,N2)
  do k=nz,CS%VarMix_Ktop,-1

    ! Calculate the interface slopes E_x and E_y and u- and v- points respectively
    do j=js-1,je+1 ; do I=is-1,ie
      E_x(I,j) = (e(i+1,j,K)-e(i,j,K))*G%IdxCu(I,j)
      ! Mask slopes where interface intersects topography
      if (min(h(i,j,k),h(i+1,j,k)) < H_cutoff) E_x(I,j) = 0.
    enddo ; enddo
    do J=js-1,je ; do i=is-1,ie+1
      E_y(i,J) = (e(i,j+1,K)-e(i,j,K))*G%IdyCv(i,J)
      ! Mask slopes where interface intersects topography
      if (min(h(i,j,k),h(i,j+1,k)) < H_cutoff) E_y(i,J) = 0.
    enddo ; enddo

    ! Calculate N*S*h from this layer and add to the sum
    do j=js,je ; do I=is-1,ie
      S2 = ( E_x(I,j)**2  + 0.25*( &
            ((E_y(i,J)**2) + (E_y(i+1,J-1)**2)) + ((E_y(i+1,J)**2) + (E_y(i,J-1)**2)) ) )
      if (min(h(i,j,k-1), h(i+1,j,k-1), h(i,j,k), h(i+1,j,k)) < H_cutoff) S2 = 0.0

      Hdn = 2.*h(i,j,k)*h(i,j,k-1) / (h(i,j,k) + h(i,j,k-1) + h_neglect)
      Hup = 2.*h(i+1,j,k)*h(i+1,j,k-1) / (h(i+1,j,k) + h(i+1,j,k-1) + h_neglect)
      H_geom = sqrt(Hdn*Hup)
      ! N2 = GV%g_prime(k) / (GV%H_to_Z * max(Hdn, Hup, CS%h_min_N2))
      S2N2_u_local(I,j,k) = (H_geom * S2) * (GV%g_prime(k) / max(Hdn, Hup, CS%h_min_N2) )
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      S2 = ( E_y(i,J)**2  + 0.25*( &
            ((E_x(I,j)**2) + (E_x(I-1,j+1)**2)) + ((E_x(I,j+1)**2) + (E_x(I-1,j)**2)) ) )
      if (min(h(i,j,k-1), h(i,j+1,k-1), h(i,j,k), h(i,j+1,k)) < H_cutoff) S2 = 0.0

      Hdn = 2.*h(i,j,k)*h(i,j,k-1) / (h(i,j,k) + h(i,j,k-1) + h_neglect)
      Hup = 2.*h(i,j+1,k)*h(i,j+1,k-1) / (h(i,j+1,k) + h(i,j+1,k-1) + h_neglect)
      H_geom = sqrt(Hdn*Hup)
      ! N2 = GV%g_prime(k) / (GV%H_to_Z * max(Hdn, Hup, CS%h_min_N2))
      S2N2_v_local(i,J,k) = (H_geom * S2) * (GV%g_prime(k) / (max(Hdn, Hup, CS%h_min_N2)))
    enddo ; enddo

  enddo ! k

  !$OMP parallel do default(shared)
  do j=js,je
    do I=is-1,ie ; CS%SN_u(I,j) = 0.0 ; enddo
    do k=nz,CS%VarMix_Ktop,-1 ; do I=is-1,ie
      CS%SN_u(I,j) = CS%SN_u(I,j) + S2N2_u_local(I,j,k)
    enddo ; enddo
    ! SN above contains S^2*N^2*H, convert to vertical average of S*N

    if (use_dztot) then
      do I=is-1,ie
        CS%SN_u(I,j) = G%OBCmaskCu(I,j) * sqrt( CS%SN_u(I,j) / &
                                                max(dz_tot(i,j), dz_tot(i+1,j), GV%dz_subroundoff) )
      enddo
    else
      do I=is-1,ie
        if ( min(G%bathyT(i,j), G%bathyT(i+1,j)) + G%Z_ref > dZ_cutoff ) then
          CS%SN_u(I,j) = G%OBCmaskCu(I,j) * sqrt( CS%SN_u(I,j) / &
                                                  (max(G%bathyT(i,j), G%bathyT(i+1,j)) + G%Z_ref) )
        else
          CS%SN_u(I,j) = 0.0
        endif
      enddo
    endif
  enddo
  !$OMP parallel do default(shared)
  do J=js-1,je
    do i=is,ie ; CS%SN_v(i,J) = 0.0 ; enddo
    do k=nz,CS%VarMix_Ktop,-1 ; do i=is,ie
      CS%SN_v(i,J) = CS%SN_v(i,J) + S2N2_v_local(i,J,k)
    enddo ; enddo
    if (use_dztot) then
      do i=is,ie
        CS%SN_v(i,J) = G%OBCmaskCv(i,J) * sqrt( CS%SN_v(i,J) / &
                                                max(dz_tot(i,j), dz_tot(i,j+1), GV%dz_subroundoff) )
      enddo
    else
      do i=is,ie
        ! There is a primordial horizontal indexing bug on the following line from the previous
        ! versions of the code.  This comment should be deleted by the end of 2024.
        ! if ( min(G%bathyT(i,j), G%bathyT(i+1,j)) + G%Z_ref > dZ_cutoff ) then
        if ( min(G%bathyT(i,j), G%bathyT(i,j+1)) + G%Z_ref > dZ_cutoff ) then
          CS%SN_v(i,J) = G%OBCmaskCv(i,J) * sqrt( CS%SN_v(i,J) / &
                                                  (max(G%bathyT(i,j), G%bathyT(i,j+1)) + G%Z_ref) )
        else
          CS%SN_v(i,J) = 0.0
        endif
      enddo
    endif
  enddo

end subroutine calc_slope_functions_using_just_e


!> Calculates and returns isopycnal slopes with wider halos for use in finding QG viscosity.
subroutine calc_QG_slopes(h, tv, dt, G, GV, US, slope_x, slope_y, CS, OBC)
  type(ocean_grid_type),                        intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),                        intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                        intent(in)    :: tv !< Thermodynamic variables
  real,                                         intent(in)    :: dt !< Time increment [T ~> s]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: slope_x !< Isopycnal slope in i-dir [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: slope_y !< Isopycnal slope in j-dir [Z L-1 ~> nondim]
  type(VarMix_CS),                              intent(in)    :: CS !< Variable mixing control structure
  type(ocean_OBC_type),                         pointer       :: OBC !< Open boundaries control structure
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1)  :: e    ! The interface heights relative to mean sea level [Z ~> m]

  if (.not. CS%initialized) call MOM_error(FATAL, "MOM_lateral_mixing_coeffs.F90, calc_QG_slopes: "//&
         "Module must be initialized before it is used.")

  call find_eta(h, tv, G, GV, US, e, halo_size=3)
  call calc_isoneutral_slopes(G, GV, US, h, e, tv, dt*CS%kappa_smooth, CS%use_stanley_iso, &
                              slope_x, slope_y, halo=2, OBC=OBC, OBC_N2=CS%OBC_friendly)

end subroutine calc_QG_slopes

!> Calculates the Leith Laplacian and bi-harmonic viscosity coefficients
subroutine calc_QG_Leith_viscosity(CS, G, GV, US, h, dz, k, div_xx_dx, div_xx_dy, slope_x, slope_y, &
                                   vort_xy_dx, vort_xy_dy)
  type(VarMix_CS),                           intent(inout) :: CS !< Variable mixing coefficients
  type(ocean_grid_type),                     intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: dz !< Layer vertical extents [Z ~> m]
  integer,                                   intent(in)    :: k  !< Layer for which to calculate vorticity magnitude
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: div_xx_dx  !< x-derivative of horizontal divergence
                                                                 !! (d/dx(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
  real, dimension(SZI_(G),SZJB_(G)),         intent(in)    :: div_xx_dy  !< y-derivative of horizontal divergence
                                                                 !! (d/dy(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: slope_x !< Isopycnal slope in i-dir [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: slope_y !< Isopycnal slope in j-dir [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G)),         intent(inout) :: vort_xy_dx !< x-derivative of vertical vorticity
                                                                 !! (d/dx(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
  real, dimension(SZIB_(G),SZJ_(G)),         intent(inout) :: vort_xy_dy !< y-derivative of vertical vorticity
                                                                 !! (d/dy(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
  ! Local variables
  real, dimension(SZI_(G),SZJB_(G)) :: &
    dslopey_dz, & ! z-derivative of y-slope at v-points [L-1 ~> m-1]
    h_at_v,     & ! Thickness at v-points [H ~> m or kg m-2]
    beta_v,     & ! Beta at v-points [T-1 L-1 ~> s-1 m-1]
    grad_vort_mag_v, & ! Magnitude of vorticity gradient at v-points [T-1 L-1 ~> s-1 m-1]
    grad_div_mag_v     ! Magnitude of divergence gradient at v-points [T-1 L-1 ~> s-1 m-1]

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    dslopex_dz, & ! z-derivative of x-slope at u-points [L-1 ~> m-1]
    h_at_u,     & ! Thickness at u-points [H ~> m or kg m-2]
    beta_u,     & ! Beta at u-points [T-1 L-1 ~> s-1 m-1]
    grad_vort_mag_u, & ! Magnitude of vorticity gradient at u-points [T-1 L-1 ~> s-1 m-1]
    grad_div_mag_u     ! Magnitude of divergence gradient at u-points [T-1 L-1 ~> s-1 m-1]
  real :: h_at_slope_above ! The thickness above [H ~> m or kg m-2]
  real :: h_at_slope_below ! The thickness below [H ~> m or kg m-2]
  real :: Ih ! The inverse of a combination of thicknesses [H-1 ~> m-1 or m2 kg-1]
  real :: f  ! A copy of the Coriolis parameter [T-1 ~> s-1]
  real :: Z_to_H  ! A local copy of depth to thickness conversion factors or the inverse of the
                  ! mass-weighted average specific volumes around an interface [H Z-1 ~> nondim or kg m-3]
  real :: inv_PI3 ! The inverse of pi cubed [nondim]
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nz = GV%ke

  inv_PI3 = 1.0 / ((4.0*atan(1.0))**3)
  Z_to_H = GV%Z_to_H  ! This will be replaced with a varying value in non-Boussinesq mode.

  if ((k > 1) .and. (k < nz)) then

    do j=js-2,je+2 ; do I=is-2,ie+1
      h_at_slope_above = 2. * ( h(i,j,k-1) * h(i+1,j,k-1) ) * ( h(i,j,k) * h(i+1,j,k) ) / &
                         ( ( h(i,j,k-1) * h(i+1,j,k-1) ) * ( h(i,j,k) + h(i+1,j,k) ) &
                         + ( h(i,j,k) * h(i+1,j,k) ) * ( h(i,j,k-1) + h(i+1,j,k-1) ) + GV%H_subroundoff**3 )
      h_at_slope_below = 2. * ( h(i,j,k) * h(i+1,j,k) ) * ( h(i,j,k+1) * h(i+1,j,k+1) ) / &
                         ( ( h(i,j,k) * h(i+1,j,k) ) * ( h(i,j,k+1) + h(i+1,j,k+1) ) &
                         + ( h(i,j,k+1) * h(i+1,j,k+1) ) * ( h(i,j,k) + h(i+1,j,k) ) + GV%H_subroundoff**3 )
      Ih = 1./ ( h_at_slope_above + h_at_slope_below + GV%H_subroundoff )
      if (.not.GV%Boussinesq) &
        Z_to_H = ( (h(i,j,k-1) + h(i+1,j,k-1)) + (h(i,j,k) + h(i+1,j,k)) ) / &
                 ( (dz(i,j,k-1) + dz(i+1,j,k-1)) + (dz(i,j,k) + dz(i+1,j,k)) + GV%dZ_subroundoff)
      dslopex_dz(I,j) = 2. * ( slope_x(I,j,k) - slope_x(I,j,k+1) ) * (Z_to_H * Ih)
      h_at_u(I,j) = 2. * ( h_at_slope_above * h_at_slope_below ) * Ih
    enddo ; enddo

    do J=js-2,je+1 ; do i=is-2,ie+2
      h_at_slope_above = 2. * ( h(i,j,k-1) * h(i,j+1,k-1) ) * ( h(i,j,k) * h(i,j+1,k) ) / &
                         ( ( h(i,j,k-1) * h(i,j+1,k-1) ) * ( h(i,j,k) + h(i,j+1,k) ) &
                         + ( h(i,j,k) * h(i,j+1,k) ) * ( h(i,j,k-1) + h(i,j+1,k-1) ) + GV%H_subroundoff**3 )
      h_at_slope_below = 2. * ( h(i,j,k) * h(i,j+1,k) ) * ( h(i,j,k+1) * h(i,j+1,k+1) ) / &
                         ( ( h(i,j,k) * h(i,j+1,k) ) * ( h(i,j,k+1) + h(i,j+1,k+1) ) &
                         + ( h(i,j,k+1) * h(i,j+1,k+1) ) * ( h(i,j,k) + h(i,j+1,k) ) + GV%H_subroundoff**3 )
      Ih = 1./ ( h_at_slope_above + h_at_slope_below + GV%H_subroundoff )
      if (.not.GV%Boussinesq) &
        Z_to_H = ( (h(i,j,k-1) + h(i,j+1,k-1)) + (h(i,j,k) + h(i,j+1,k)) ) / &
                 ( (dz(i,j,k-1) + dz(i,j+1,k-1)) + (dz(i,j,k) + dz(i,j+1,k)) + GV%dZ_subroundoff)
      dslopey_dz(i,J) = 2. * ( slope_y(i,J,k) - slope_y(i,J,k+1) ) * (Z_to_H * Ih)
      h_at_v(i,J) = 2. * ( h_at_slope_above * h_at_slope_below ) * Ih
    enddo ; enddo

    do J=js-2,je+1 ; do i=is-1,ie+1
      f = 0.5 * ( G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J) )
      vort_xy_dx(i,J) = vort_xy_dx(i,J) - f * &
            ( ( (h_at_u(I,j) * dslopex_dz(I,j)) + (h_at_u(I-1,j+1) * dslopex_dz(I-1,j+1)) ) &
            + ( (h_at_u(I-1,j) * dslopex_dz(I-1,j)) + (h_at_u(I,j+1) * dslopex_dz(I,j+1)) ) ) / &
              ( ( h_at_u(I,j) + h_at_u(I-1,j+1) ) + ( h_at_u(I-1,j) + h_at_u(I,j+1) ) + GV%H_subroundoff)
    enddo ; enddo

    do j=js-1,je+1 ; do I=is-2,ie+1
      f = 0.5 * ( G%CoriolisBu(I,J) + G%CoriolisBu(I,J-1) )
      vort_xy_dy(I,j) = vort_xy_dy(I,j) - f * &
            ( ( (h_at_v(i,J) * dslopey_dz(i,J)) + (h_at_v(i+1,J-1) * dslopey_dz(i+1,J-1)) ) &
            + ( (h_at_v(i,J-1) * dslopey_dz(i,J-1)) + (h_at_v(i+1,J) * dslopey_dz(i+1,J)) ) ) / &
              ( ( h_at_v(i,J) + h_at_v(i+1,J-1) ) + ( h_at_v(i,J-1) + h_at_v(i+1,J) ) + GV%H_subroundoff)
    enddo ; enddo
  endif ! k > 1

  if (CS%use_QG_Leith_GM) then

    do j=js,je ; do I=is-1,Ieq
      grad_vort_mag_u(I,j) = SQRT(vort_xy_dy(I,j)**2  + (0.25*((vort_xy_dx(i,J) + vort_xy_dx(i+1,J-1)) &
                                                             + (vort_xy_dx(i+1,J) + vort_xy_dx(i,J-1))))**2)
      grad_div_mag_u(I,j) = SQRT(div_xx_dx(I,j)**2  + (0.25*((div_xx_dy(i,J) + div_xx_dy(i+1,J-1)) &
                                                           + (div_xx_dy(i+1,J) + div_xx_dy(i,J-1))))**2)
      if (CS%use_beta_in_QG_Leith) then
        beta_u(I,j) = sqrt((0.5*(G%dF_dx(i,j)+G%dF_dx(i+1,j))**2) + &
                           (0.5*(G%dF_dy(i,j)+G%dF_dy(i+1,j))**2))
        CS%KH_u_QG(I,j,k) = MIN(grad_vort_mag_u(I,j) + grad_div_mag_u(I,j), 3.0*beta_u(I,j)) * &
                            CS%Laplac3_const_u(I,j) * inv_PI3
      else
        CS%KH_u_QG(I,j,k) = (grad_vort_mag_u(I,j) + grad_div_mag_u(I,j)) * &
                            CS%Laplac3_const_u(I,j) * inv_PI3
      endif
    enddo ; enddo

    do J=js-1,Jeq ; do i=is,ie
      grad_vort_mag_v(i,J) = SQRT(vort_xy_dx(i,J)**2  + (0.25*((vort_xy_dy(I,j) + vort_xy_dy(I-1,j+1)) &
                                                             + (vort_xy_dy(I,j+1) + vort_xy_dy(I-1,j))))**2)
      grad_div_mag_v(i,J) = SQRT(div_xx_dy(i,J)**2  + (0.25*((div_xx_dx(I,j) + div_xx_dx(I-1,j+1)) &
                                                           + (div_xx_dx(I,j+1) + div_xx_dx(I-1,j))))**2)
      if (CS%use_beta_in_QG_Leith) then
        beta_v(i,J) = sqrt((0.5*(G%dF_dx(i,j)+G%dF_dx(i,j+1))**2) + &
                           (0.5*(G%dF_dy(i,j)+G%dF_dy(i,j+1))**2))
        CS%KH_v_QG(i,J,k) = MIN(grad_vort_mag_v(i,J) + grad_div_mag_v(i,J), 3.0*beta_v(i,J)) * &
                            CS%Laplac3_const_v(i,J) * inv_PI3
      else
        CS%KH_v_QG(i,J,k) = (grad_vort_mag_v(i,J) + grad_div_mag_v(i,J)) * &
                            CS%Laplac3_const_v(i,J) * inv_PI3
      endif
    enddo ; enddo
    ! post diagnostics

    if (k==nz) then
      if (CS%id_KH_v_QG > 0)  call post_data(CS%id_KH_v_QG, CS%KH_v_QG, CS%diag)
      if (CS%id_KH_u_QG > 0)  call post_data(CS%id_KH_u_QG, CS%KH_u_QG, CS%diag)
    endif
  endif

end subroutine calc_QG_Leith_viscosity

!> Initializes the variables mixing coefficients container
subroutine VarMix_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type),            intent(in) :: Time !< Current model time
  type(ocean_grid_type),      intent(in) :: G    !< Ocean grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(VarMix_CS),         intent(inout) :: CS   !< Variable mixing coefficients

  ! Local variables
  real :: KhTr_Slope_Cff ! The nondimensional coefficient in the Visbeck formula
                         ! for the epipycnal tracer diffusivity [nondim]
  real :: KhTh_Slope_Cff ! The nondimensional coefficient in the Visbeck formula
                         ! for the interface depth diffusivity [nondim]
  real :: oneOrTwo ! A variable that may be 1 or 2, depending on which form
                   ! of the equatorial deformation radius us used [nondim]
  real :: N2_filter_depth  ! A depth below which stratification is treated as monotonic when
                           ! calculating the first-mode wave speed [H ~> m or kg m-2]
  real :: KhTr_passivity_coeff ! Coefficient setting the ratio between along-isopycnal tracer
                               ! mixing and interface height mixing [nondim]
  real :: absurdly_small_freq  ! A miniscule frequency that is used to avoid division by 0 [T-1 ~> s-1].  The
             ! default value is roughly (pi / (the age of the universe)).
  logical :: Gill_equatorial_Ld, use_FGNV_streamfn, use_MEKE, in_use
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  integer :: remap_answer_date    ! The vintage of the order of arithmetic and expressions to use
                                  ! for remapping.  Values below 20190101 recover the remapping
                                  ! answers from 2018, while higher values use more robust
                                  ! forms of the same remapping expressions.
  real :: MLE_front_length        ! The frontal-length scale used to calculate the upscaling of
                                  ! buoyancy gradients in boundary layer parameterizations [L ~> m]
  real :: Leith_Lap_const      ! The non-dimensional coefficient in the Leith viscosity [nondim]
  real :: grid_sp_u2, grid_sp_v2 ! Intermediate quantities for Leith metrics [L2 ~> m2]
  real :: grid_sp_u3, grid_sp_v3 ! Intermediate quantities for Leith metrics [L3 ~> m3]
  real :: wave_speed_min      ! A floor in the first mode speed below which 0 is returned [L T-1 ~> m s-1]
  real :: wave_speed_tol      ! The fractional tolerance for finding the wave speeds [nondim]
  logical :: Resoln_scaled_MEKE_visc ! If true, the viscosity contribution from MEKE is
                                  ! scaled by the resolution function.
  logical :: better_speed_est ! If true, use a more robust estimate of the first
                              ! mode wave speed as the starting point for iterations.
  real :: Stanley_coeff    ! Coefficient relating the temperature gradient and sub-gridscale
                           ! temperature variance [nondim]
  logical :: om4_remap_via_sub_cells ! Use the OM4-era remap_via_sub_cells for calculating the EBT structure
  logical :: enable_bugs   ! If true, the defaults for recently added bug-fix flags are set to
                           ! recreate the bugs, or if false bugs are only used if actively selected.
  logical :: mixing_coefs_OBC_bug ! If false, use only interior data for thickness weighting in
                           ! lateral mixing coefficient calculations and to calculate stratification
                           ! and other fields at open boundary condition faces.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_lateral_mixing_coeffs" ! This module's name.
  integer :: number_of_OBC_segments
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, i, j
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%initialized = .true.
  in_use = .false. ! Set to true to avoid deallocating
  CS%diag => diag ! Diagnostics pointer
  CS%calculate_cg1 = .false.
  CS%calculate_Rd_dx = .false.
  CS%calculate_res_fns = .false.
  CS%use_simpler_Eady_growth_rate  = .false.
  CS%full_depth_Eady_growth_rate = .false.
  CS%calculate_depth_fns = .false.
  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "USE_VARIABLE_MIXING", CS%use_variable_mixing,&
                 "If true, the variable mixing code will be called.  This "//&
                 "allows diagnostics to be created even if the scheme is "//&
                 "not used.  If KHTR_SLOPE_CFF>0 or  KhTh_Slope_Cff>0, "//&
                 "this is set to true regardless of what is in the "//&
                 "parameter file.", default=.false.)
  call get_param(param_file, mdl, "USE_VISBECK", CS%use_Visbeck,&
                 "If true, use the Visbeck et al. (1997) formulation for \n"//&
                 "thickness diffusivity.", default=.false.)
  call get_param(param_file, mdl, "RESOLN_SCALED_KH", CS%Resoln_scaled_Kh, &
                 "If true, the Laplacian lateral viscosity is scaled away "//&
                 "when the first baroclinic deformation radius is well "//&
                 "resolved.", default=.false.)
  call get_param(param_file, mdl, "DEPTH_SCALED_KHTH", CS%Depth_scaled_KhTh, &
                 "If true, KHTH is scaled away when the depth is shallower"//&
                 "than a reference depth: KHTH = MIN(1,H/H0)**N * KHTH, "//&
                 "where H0 is a reference depth, controlled via DEPTH_SCALED_KHTH_H0, "//&
                 "and the exponent (N) is controlled via DEPTH_SCALED_KHTH_EXP.",&
                 default=.false.)
  call get_param(param_file, mdl, "RESOLN_SCALED_KHTH", CS%Resoln_scaled_KhTh, &
                 "If true, the interface depth diffusivity is scaled away "//&
                 "when the first baroclinic deformation radius is well "//&
                 "resolved.", default=.false.)
  call get_param(param_file, mdl, "RESOLN_SCALED_KHTR", CS%Resoln_scaled_KhTr, &
                 "If true, the epipycnal tracer diffusivity is scaled "//&
                 "away when the first baroclinic deformation radius is "//&
                 "well resolved.", default=.false.)
  call get_param(param_file, mdl, "USE_MEKE", use_MEKE, &
                 default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "RES_SCALE_MEKE_VISC", Resoln_scaled_MEKE_visc, &
                 "If true, the viscosity contribution from MEKE is scaled by "//&
                 "the resolution function.", default=.false., do_not_log=.true.) ! Logged elsewhere.
  if (.not.use_MEKE) Resoln_scaled_MEKE_visc = .false.
  call get_param(param_file, mdl, "RESOLN_USE_EBT", CS%Resoln_use_ebt, &
                 "If true, uses the equivalent barotropic wave speed instead "//&
                 "of first baroclinic wave for calculating the resolution fn.",&
                 default=.false.)
  call get_param(param_file, mdl, "BACKSCAT_EBT_POWER", CS%BS_EBT_power, &
                 "Power to raise EBT vertical structure to when backscatter "// &
                 "has vertical structure.", units="nondim", default=0.0)
  call get_param(param_file, mdl, "BS_USE_SQG_STRUCT", CS%BS_use_sqg_struct, &
                 "If true, the SQG vertical structure is used for backscatter "//&
                 "on the condition that BS_EBT_power=0", &
                 default=.false.)
  call get_param(param_file, mdl, "SQG_EXPO", CS%sqg_expo, &
                 "Nondimensional exponent coeffecient of the SQG mode "// &
                 "that is used for the vertical struture of diffusivities.", units="nondim", default=1.0)
  call get_param(param_file, mdl, "KHTH_USE_EBT_STRUCT", CS%khth_use_ebt_struct, &
                 "If true, uses the equivalent barotropic structure "//&
                 "as the vertical structure of thickness diffusivity.",&
                 default=.false.)
  call get_param(param_file, mdl, "KHTH_USE_SQG_STRUCT", CS%khth_use_sqg_struct, &
                 "If true, uses the surface quasigeostrophic structure "//&
                 "as the vertical structure of thickness diffusivity.",&
                 default=.false.)
  call get_param(param_file, mdl, "KHTR_USE_EBT_STRUCT", CS%khtr_use_ebt_struct, &
                 "If true, uses the equivalent barotropic structure "//&
                 "as the vertical structure of tracer diffusivity.",&
                 default=.false.)
  call get_param(param_file, mdl, "KHTR_USE_SQG_STRUCT", CS%khtr_use_sqg_struct, &
                 "If true, uses the surface quasigeostrophic structure "//&
                 "as the vertical structure of tracer diffusivity.",&
                 default=.false.)
  call get_param(param_file, mdl, "KD_GL90_USE_EBT_STRUCT", CS%kdgl90_use_ebt_struct, &
                 "If true, uses the equivalent barotropic structure "//&
                 "as the vertical structure of diffusivity in the GL90 scheme.",&
                 default=.false.)
  call get_param(param_file, mdl, "KD_GL90_USE_SQG_STRUCT", CS%kdgl90_use_sqg_struct, &
                 "If true, uses the equivalent barotropic structure "//&
                 "as the vertical structure of diffusivity in the GL90 scheme.",&
                 default=.false.)
  call get_param(param_file, mdl, "KHTH_SLOPE_CFF", KhTh_Slope_Cff, &
                 "The nondimensional coefficient in the Visbeck formula "//&
                 "for the interface depth diffusivity", units="nondim", default=0.0)
  call get_param(param_file, mdl, "KHTR_SLOPE_CFF", KhTr_Slope_Cff, &
                 "The nondimensional coefficient in the Visbeck formula "//&
                 "for the epipycnal tracer diffusivity", units="nondim", default=0.0)
  call get_param(param_file, mdl, "USE_STORED_SLOPES", CS%use_stored_slopes,&
                 "If true, the isopycnal slopes are calculated once and "//&
                 "stored for re-use. This uses more memory but avoids calling "//&
                 "the equation of state more times than should be necessary.", &
                 default=.false.)
  call get_param(param_file, mdl, "VERY_SMALL_FREQUENCY", absurdly_small_freq, &
                 "A miniscule frequency that is used to avoid division by 0.  The default "//&
                 "value is roughly (pi / (the age of the universe)).", &
                 default=1.0e-17, units="s-1", scale=US%T_to_s)
  call get_param(param_file, mdl, "KHTH_USE_FGNV_STREAMFUNCTION", use_FGNV_streamfn, &
                 default=.false., do_not_log=.true.)
  CS%calculate_cg1 = CS%calculate_cg1 .or. use_FGNV_streamfn .or. CS%khth_use_ebt_struct &
                     .or. CS%kdgl90_use_ebt_struct .or. CS%BS_EBT_power>0.
  CS%calculate_Rd_dx = CS%calculate_Rd_dx .or. use_MEKE
  ! Indicate whether to calculate the Eady growth rate
  CS%calculate_Eady_growth_rate = use_MEKE .or. (KhTr_Slope_Cff>0.) .or. (KhTh_Slope_Cff>0.)
  call get_param(param_file, mdl, "KHTR_PASSIVITY_COEFF", KhTr_passivity_coeff, &
                 units="nondim", default=0., do_not_log=.true.)
  CS%calculate_Rd_dx = CS%calculate_Rd_dx .or. (KhTr_passivity_coeff>0.)
  call get_param(param_file, mdl, "MLE_FRONT_LENGTH", MLE_front_length, &
                 units="m", default=0.0, scale=US%m_to_L, do_not_log=.true.)
  CS%calculate_Rd_dx = CS%calculate_Rd_dx .or. (MLE_front_length>0.)

  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false., do_not_log=.true.)

  call get_param(param_file, mdl, "USE_STANLEY_ISO", CS%use_stanley_iso, &
                 "If true, turn on Stanley SGS T variance parameterization "// &
                 "in isopycnal slope code.", default=.false.)
  if (CS%use_stanley_iso) then
    call get_param(param_file, mdl, "STANLEY_COEFF", Stanley_coeff, &
                 "Coefficient correlating the temperature gradient and SGS T variance.", &
                 units="nondim", default=-1.0, do_not_log=.true.)
    if (Stanley_coeff < 0.0) call MOM_error(FATAL, &
                 "STANLEY_COEFF must be set >= 0 if USE_STANLEY_ISO is true.")
  endif
  call get_param(param_file, mdl, "OBC_NUMBER_OF_SEGMENTS", number_of_OBC_segments, &
                 default=0, do_not_log=.true.)
  call get_param(param_file, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
  call get_param(param_file, mdl, "MIXING_COEFS_OBC_BUG", mixing_coefs_OBC_bug, &
                 "If false, use only interior data for thickness weighting in lateral mixing "//&
                 "coefficient calculations and to calculate stratification and other fields at "//&
                 "open boundary condition faces.", &
                 default=enable_bugs, do_not_log=(number_of_OBC_segments<=0))
  CS%OBC_friendly = .not. MIXING_COEFS_OBC_BUG

  if (CS%Resoln_use_ebt .or. CS%khth_use_ebt_struct .or. CS%kdgl90_use_ebt_struct &
      .or. CS%BS_EBT_power>0. .or. CS%khtr_use_ebt_struct) then
    in_use = .true.
    call get_param(param_file, mdl, "RESOLN_N2_FILTER_DEPTH", N2_filter_depth, &
                 "The depth below which N2 is monotonized to avoid stratification "//&
                 "artifacts from altering the equivalent barotropic mode structure.  "//&
                 "This monotonzization is disabled if this parameter is negative.", &
                 units="m", default=-1.0, scale=GV%m_to_H)
    allocate(CS%ebt_struct(isd:ied,jsd:jed,GV%ke), source=0.0)
  endif


  if (CS%BS_EBT_power>0. .and. CS%BS_use_sqg_struct) then
    call MOM_error(FATAL, &
                   "calc_resoln_function: BS_EBT_POWER>0. &
                   & and BS_USE_SQG=True cannot be set together")
  endif

  if (CS%khth_use_ebt_struct .and. CS%khth_use_sqg_struct) then
    call MOM_error(FATAL, &
                   "calc_resoln_function: Only one of KHTH_USE_EBT_STRUCT &
                   & and KHTH_USE_SQG_STRUCT can be true")
  endif

  if (CS%khtr_use_ebt_struct .and. CS%khtr_use_sqg_struct) then
    call MOM_error(FATAL, &
                   "calc_resoln_function: Only one of KHTR_USE_EBT_STRUCT &
                   & and KHTR_USE_SQG_STRUCT can be true")
  endif

  if (CS%kdgl90_use_ebt_struct .and. CS%kdgl90_use_sqg_struct) then
    call MOM_error(FATAL, &
                   "calc_resoln_function: Only one of KD_GL90_USE_EBT_STRUCT &
                   & and KD_GL90_USE_SQG_STRUCT can be true")
  endif

  if (CS%BS_EBT_power>0. .or. CS%BS_use_sqg_struct) then
    allocate(CS%BS_struct(isd:ied,jsd:jed,GV%ke), source=0.0)
  endif

  if (CS%khth_use_ebt_struct .or. CS%khth_use_sqg_struct) then
    allocate(CS%khth_struct(isd:ied, jsd:jed, gv%ke), source=0.0)
  endif

  if (CS%khtr_use_ebt_struct .or. CS%khtr_use_sqg_struct) then
    allocate(CS%khtr_struct(isd:ied, jsd:jed, gv%ke), source=0.0)
  endif

  if (CS%kdgl90_use_ebt_struct .or. CS%kdgl90_use_sqg_struct) then
    allocate(CS%kdgl90_struct(isd:ied, jsd:jed, gv%ke), source=0.0)
  endif

  if (CS%use_stored_slopes) then
    if (KhTr_Slope_Cff>0. .or. KhTh_Slope_Cff>0.) then
      call get_param(param_file, mdl, "VISBECK_MAX_SLOPE", CS%Visbeck_S_max, &
            "If non-zero, is an upper bound on slopes used in the "//&
            "Visbeck formula for diffusivity. This does not affect the "//&
            "isopycnal slope calculation used within thickness diffusion.",  &
            units="nondim", default=0.0, scale=US%L_to_Z)
    else
      CS%Visbeck_S_max = 0.
    endif
  endif

  if (CS%use_stored_slopes .or. CS%sqg_expo>0.0) then
    ! CS%calculate_Eady_growth_rate=.true.
    in_use = .true.
    allocate(CS%slope_x(IsdB:IedB,jsd:jed,GV%ke+1), source=0.0)
    allocate(CS%slope_y(isd:ied,JsdB:JedB,GV%ke+1), source=0.0)
    call get_param(param_file, mdl, "KD_SMOOTH", CS%kappa_smooth, &
                 "A diapycnal diffusivity that is used to interpolate "//&
                 "more sensible values of T & S into thin layers.", &
                 units="m2 s-1", default=1.0e-6, scale=GV%m2_s_to_HZ_T)
  endif

  if (CS%calculate_Eady_growth_rate) then
    in_use = .true.
    allocate(CS%SN_u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%SN_v(isd:ied,JsdB:JedB), source=0.0)
    CS%id_SN_u = register_diag_field('ocean_model', 'SN_u', diag%axesCu1, Time, &
       'Inverse eddy time-scale, S*N, at u-points', 's-1', conversion=US%s_to_T)
    CS%id_SN_v = register_diag_field('ocean_model', 'SN_v', diag%axesCv1, Time, &
       'Inverse eddy time-scale, S*N, at v-points', 's-1', conversion=US%s_to_T)
    call get_param(param_file, mdl, "USE_SIMPLER_EADY_GROWTH_RATE", CS%use_simpler_Eady_growth_rate, &
                   "If true, use a simpler method to calculate the Eady growth rate "//&
                   "that avoids division by layer thickness. Recommended.", default=.false.)
    if (CS%use_simpler_Eady_growth_rate) then
      if (.not. CS%use_stored_slopes) call MOM_error(FATAL, &
           "MOM_lateral_mixing_coeffs.F90, VarMix_init:"//&
           "When USE_SIMPLER_EADY_GROWTH_RATE=True, USE_STORED_SLOPES must also be True.")
      call get_param(param_file, mdl, "EADY_GROWTH_RATE_D_SCALE", CS%Eady_GR_D_scale, &
                     "The depth from surface over which to average SN when calculating "//&
                     "a 2D Eady growth rate. Zero mean use full depth.", &
                      units="m", default=0., scale=US%m_to_Z)
      call get_param(param_file, mdl, "EADY_GROWTH_RATE_CROPPING_DISTANCE", CS%cropping_distance, &
                     "Distance from surface or bottom to filter out outcropped or "//&
                     "incropped interfaces for the Eady growth rate calc. "//&
                     "Negative values disables cropping.", units="m", default=0., scale=US%m_to_Z)
    else
      call get_param(param_file, mdl, "VARMIX_KTOP", CS%VarMix_Ktop, &
                     "The layer number at which to start vertical integration "//&
                     "of S*N for purposes of finding the Eady growth rate.", &
                     units="nondim", default=2)
      call get_param(param_file, mdl, "MIN_DZ_FOR_SLOPE_N2", CS%h_min_N2, &
                     "The minimum vertical distance to use in the denominator of the "//&
                     "bouyancy frequency used in the slope calculation.", &
                     units="m", default=1.0, scale=GV%m_to_H, do_not_log=CS%use_stored_slopes)

      call get_param(param_file, mdl, "FULL_DEPTH_EADY_GROWTH_RATE", CS%full_depth_Eady_growth_rate, &
                   "If true, calculate the Eady growth rate based on average slope times "//&
                   "stratification that includes contributions from sea-level changes "//&
                   "in its denominator, rather than just the nominal depth of the bathymetry.  "//&
                   "This only applies when using the model interface heights as a proxy for "//&
                   "isopycnal slopes.", default=.not.(GV%Boussinesq.or.GV%semi_Boussinesq), &
                   do_not_log=CS%use_stored_slopes)
    endif
  endif

  if (KhTr_Slope_Cff>0. .or. KhTh_Slope_Cff>0.) then
    in_use = .true.
    call get_param(param_file, mdl, "VISBECK_L_SCALE", CS%Visbeck_L_scale, &
                 "The fixed length scale in the Visbeck formula, or if negative a nondimensional "//&
                 "scaling factor relating this length scale squared to the cell areas.", &
                 units="m or nondim", default=0.0, scale=US%m_to_L)
    allocate(CS%L2u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%L2v(isd:ied,JsdB:JedB), source=0.0)
    if (CS%Visbeck_L_scale<0) then
      ! Undo the rescaling of CS%Visbeck_L_scale.
      do j=js,je ; do I=is-1,Ieq
        CS%L2u(I,j) = (US%L_to_m*CS%Visbeck_L_scale)**2 * G%areaCu(I,j)
      enddo ; enddo
      do J=js-1,Jeq ; do i=is,ie
        CS%L2v(i,J) = (US%L_to_m*CS%Visbeck_L_scale)**2 * G%areaCv(i,J)
      enddo ; enddo
    else
      CS%L2u(:,:) = CS%Visbeck_L_scale**2
      CS%L2v(:,:) = CS%Visbeck_L_scale**2
    endif

    CS%id_L2u = register_diag_field('ocean_model', 'L2u', diag%axesCu1, Time, &
       'Length scale squared for mixing coefficient, at u-points', &
       'm2', conversion=US%L_to_m**2)
    CS%id_L2v = register_diag_field('ocean_model', 'L2v', diag%axesCv1, Time, &
       'Length scale squared for mixing coefficient, at v-points', &
       'm2', conversion=US%L_to_m**2)
  endif

  CS%id_sqg_struct = register_diag_field('ocean_model', 'sqg_struct', diag%axesTl, Time, &
            'Vertical structure of SQG mode', 'nondim')
  if (CS%BS_use_sqg_struct .or. CS%khth_use_sqg_struct .or. CS%khtr_use_sqg_struct &
      .or. CS%kdgl90_use_sqg_struct .or. CS%id_sqg_struct>0) then
    allocate(CS%sqg_struct(isd:ied,jsd:jed,GV%ke), source=0.0)
  endif

  if (CS%BS_EBT_power>0. .or. CS%BS_use_sqg_struct) then
    CS%id_BS_struct = register_diag_field('ocean_model', 'BS_struct', diag%axesTl, Time, &
              'Vertical structure of backscatter', 'nondim')
  endif
  if (CS%khth_use_ebt_struct .or. CS%khth_use_sqg_struct) then
    CS%id_khth_struct = register_diag_field('ocean_model', 'khth_struct', diag%axesTl, Time, &
            'Vertical structure of thickness diffusivity', 'nondim')
  endif
  if (CS%khtr_use_ebt_struct .or. CS%khtr_use_sqg_struct) then
    CS%id_khtr_struct = register_diag_field('ocean_model', 'khtr_struct', diag%axesTl, Time, &
            'Vertical structure of tracer diffusivity', 'nondim')
  endif
  if (CS%kdgl90_use_ebt_struct .or. CS%kdgl90_use_sqg_struct) then
    CS%id_kdgl90_struct = register_diag_field('ocean_model', 'kdgl90_struct', diag%axesTl, Time, &
            'Vertical structure of GL90 diffusivity', 'nondim')
  endif

  if ((CS%calculate_Eady_growth_rate .and. CS%use_stored_slopes) ) then
    CS%id_N2_u = register_diag_field('ocean_model', 'N2_u', diag%axesCui, Time, &
         'Square of Brunt-Vaisala frequency, N^2, at u-points, as used in Visbeck et al.', &
         's-2', conversion=(US%L_to_Z*US%s_to_T)**2)
    CS%id_N2_v = register_diag_field('ocean_model', 'N2_v', diag%axesCvi, Time, &
         'Square of Brunt-Vaisala frequency, N^2, at v-points, as used in Visbeck et al.', &
         's-2', conversion=(US%L_to_Z*US%s_to_T)**2)
  endif
  if (CS%use_simpler_Eady_growth_rate) then
    CS%id_dzu = register_diag_field('ocean_model', 'dzu_Visbeck', diag%axesCui, Time, &
         'dz at u-points, used in calculating Eady growth rate in Visbeck et al..', &
         'm', conversion=US%Z_to_m)
    CS%id_dzv = register_diag_field('ocean_model', 'dzv_Visbeck', diag%axesCvi, Time, &
         'dz at v-points, used in calculating Eady growth rate in Visbeck et al..', &
         'm', conversion=US%Z_to_m)
    CS%id_dzSxN = register_diag_field('ocean_model', 'dzSxN', diag%axesCui, Time, &
         'dz * |slope_x| * N, used in calculating Eady growth rate in '//&
         'Visbeck et al..', 'm s-1', conversion=US%Z_to_m*US%s_to_T)
    CS%id_dzSyN = register_diag_field('ocean_model', 'dzSyN', diag%axesCvi, Time, &
         'dz * |slope_y| * N, used in calculating Eady growth rate in '//&
         'Visbeck et al..', 'm s-1', conversion=US%Z_to_m*US%s_to_T)
  endif
  if (CS%use_stored_slopes) then
    CS%id_S2_u = register_diag_field('ocean_model', 'S2_u', diag%axesCu1, Time, &
         'Depth average square of slope magnitude, S^2, at u-points, as used in Visbeck et al.', &
         'nondim', conversion=US%Z_to_L**2)
    CS%id_S2_v = register_diag_field('ocean_model', 'S2_v', diag%axesCv1, Time, &
         'Depth average square of slope magnitude, S^2, at v-points, as used in Visbeck et al.', &
         'nondim', conversion=US%Z_to_L**2)
  endif

  oneOrTwo = 1.0
  CS%Resoln_scaling_used = CS%Resoln_scaled_Kh .or. CS%Resoln_scaled_KhTh .or. &
                           CS%Resoln_scaled_KhTr .or. Resoln_scaled_MEKE_visc
  if (CS%Resoln_scaling_used) then
    CS%calculate_Rd_dx = .true.
    CS%calculate_res_fns = .true.
    allocate(CS%Res_fn_h(isd:ied,jsd:jed), source=0.0)
    allocate(CS%Res_fn_q(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%Res_fn_u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%Res_fn_v(isd:ied,JsdB:JedB), source=0.0)
    allocate(CS%beta_dx2_q(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%beta_dx2_u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%beta_dx2_v(isd:ied,JsdB:JedB), source=0.0)
    allocate(CS%f2_dx2_q(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%f2_dx2_u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%f2_dx2_v(isd:ied,JsdB:JedB), source=0.0)

    CS%id_Res_fn = register_diag_field('ocean_model', 'Res_fn', diag%axesT1, Time, &
       'Resolution function for scaling diffusivities', 'nondim')

    call get_param(param_file, mdl, "KH_RES_SCALE_COEF", CS%Res_coef_khth, &
                 "A coefficient that determines how KhTh is scaled away if "//&
                 "RESOLN_SCALED_... is true, as "//&
                 "F = 1 / (1 + (KH_RES_SCALE_COEF*Rd/dx)^KH_RES_FN_POWER).", &
                 units="nondim", default=1.0)
    call get_param(param_file, mdl, "KH_RES_FN_POWER", CS%Res_fn_power_khth, &
                 "The power of dx/Ld in the Kh resolution function.  Any "//&
                 "positive integer may be used, although even integers "//&
                 "are more efficient to calculate.  Setting this greater "//&
                 "than 100 results in a step-function being used.", &
                 default=2)
    call get_param(param_file, mdl, "VISC_RES_SCALE_COEF", CS%Res_coef_visc, &
                 "A coefficient that determines how Kh is scaled away if "//&
                 "RESOLN_SCALED_... is true, as "//&
                 "F = 1 / (1 + (KH_RES_SCALE_COEF*Rd/dx)^KH_RES_FN_POWER). "//&
                 "This function affects lateral viscosity, Kh, and not KhTh.", &
                 units="nondim", default=CS%Res_coef_khth)
    call get_param(param_file, mdl, "VISC_RES_FN_POWER", CS%Res_fn_power_visc, &
                 "The power of dx/Ld in the Kh resolution function.  Any "//&
                 "positive integer may be used, although even integers "//&
                 "are more efficient to calculate.  Setting this greater "//&
                 "than 100 results in a step-function being used. "//&
                 "This function affects lateral viscosity, Kh, and not KhTh.", &
                 default=CS%Res_fn_power_khth)
    call get_param(param_file, mdl, "INTERPOLATE_RES_FN", CS%interpolate_Res_fn, &
                 "If true, interpolate the resolution function to the "//&
                 "velocity points from the thickness points; otherwise "//&
                 "interpolate the wave speed and calculate the resolution "//&
                 "function independently at each point.", default=.false.)
    if (CS%interpolate_Res_fn) then
      if (CS%Res_coef_visc /= CS%Res_coef_khth) call MOM_error(FATAL, &
           "MOM_lateral_mixing_coeffs.F90, VarMix_init:"//&
           "When INTERPOLATE_RES_FN=True, VISC_RES_FN_POWER must equal KH_RES_SCALE_COEF.")
      if (CS%Res_fn_power_visc /= CS%Res_fn_power_khth) call MOM_error(FATAL, &
           "MOM_lateral_mixing_coeffs.F90, VarMix_init:"//&
           "When INTERPOLATE_RES_FN=True, VISC_RES_FN_POWER must equal KH_RES_FN_POWER.")
    endif
    call get_param(param_file, mdl, "GILL_EQUATORIAL_LD", Gill_equatorial_Ld, &
                 "If true, uses Gill's definition of the baroclinic "//&
                 "equatorial deformation radius, otherwise, if false, use "//&
                 "Pedlosky's definition. These definitions differ by a factor "//&
                 "of 2 in front of the beta term in the denominator. Gill's "//&
                 "is the more appropriate definition.", default=.true.)
    if (Gill_equatorial_Ld) then
      oneOrTwo = 2.0
    endif

    do J=js-1,Jeq ; do I=is-1,Ieq
      CS%f2_dx2_q(I,J) = ((G%dxBu(I,J)**2) + (G%dyBu(I,J)**2)) * &
                         max(G%Coriolis2Bu(I,J), absurdly_small_freq**2)
      CS%beta_dx2_q(I,J) = oneOrTwo * ((G%dxBu(I,J)**2) + (G%dyBu(I,J)**2)) * (sqrt(0.5 * &
          ( ((((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2) + &
             (((G%CoriolisBu(I+1,J)-G%CoriolisBu(I,J)) * G%IdxCv(i+1,J))**2)) + &
            ((((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2) + &
             (((G%CoriolisBu(I,J+1)-G%CoriolisBu(I,J)) * G%IdyCu(I,j+1))**2)) ) ))
    enddo ; enddo

    do j=js,je ; do I=is-1,Ieq
      CS%f2_dx2_u(I,j) = ((G%dxCu(I,j)**2) + (G%dyCu(I,j)**2)) * &
          max(0.5* (G%Coriolis2Bu(I,J)+G%Coriolis2Bu(I,J-1)), absurdly_small_freq**2)
      CS%beta_dx2_u(I,j) = oneOrTwo * ((G%dxCu(I,j)**2) + (G%dyCu(I,j)**2)) * (sqrt( &
          ((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
          0.25*( ((((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2) + &
                  (((G%CoriolisBu(I+1,J)-G%CoriolisBu(I,J)) * G%IdxCv(i+1,J))**2)) + &
                 ((((G%CoriolisBu(I+1,J-1)-G%CoriolisBu(I,J-1)) * G%IdxCv(i+1,J-1))**2) + &
                  (((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2)) ) ))
    enddo ; enddo

    do J=js-1,Jeq ; do i=is,ie
      CS%f2_dx2_v(i,J) = ((G%dxCv(i,J)**2) + (G%dyCv(i,J)**2)) * &
          max(0.5*(G%Coriolis2Bu(I,J)+G%Coriolis2Bu(I-1,J)), absurdly_small_freq**2)
      CS%beta_dx2_v(i,J) = oneOrTwo * ((G%dxCv(i,J)**2) + (G%dyCv(i,J)**2)) * (sqrt( &
          ((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
          0.25*( ((((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2) + &
                  (((G%CoriolisBu(I-1,J+1)-G%CoriolisBu(I-1,J)) * G%IdyCu(I-1,j+1))**2)) + &
                 ((((G%CoriolisBu(I,J+1)-G%CoriolisBu(I,J)) * G%IdyCu(I,j+1))**2) + &
                  (((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2)) ) ))
    enddo ; enddo

  endif

  if (CS%Depth_scaled_KhTh) then
    CS%calculate_depth_fns = .true.
    allocate(CS%Depth_fn_u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%Depth_fn_v(isd:ied,JsdB:JedB), source=0.0)
    call get_param(param_file, mdl, "DEPTH_SCALED_KHTH_H0", CS%depth_scaled_khth_h0, &
                   "The depth above which KHTH is scaled away.", &
                   units="m", scale=US%m_to_Z, default=1000.)
    call get_param(param_file, mdl, "DEPTH_SCALED_KHTH_EXP", CS%depth_scaled_khth_exp, &
                   "The exponent used in the depth dependent scaling function for KHTH.", &
                   units="nondim", default=3.0)
  endif

  ! Resolution %Rd_dx_h
  CS%id_Rd_dx = register_diag_field('ocean_model', 'Rd_dx', diag%axesT1, Time, &
       'Ratio between deformation radius and grid spacing', 'm m-1')
  CS%calculate_Rd_dx = CS%calculate_Rd_dx .or. (CS%id_Rd_dx>0)

  if (CS%calculate_Rd_dx) then
    CS%calculate_cg1 = .true. ! We will need %cg1
    allocate(CS%Rd_dx_h(isd:ied,jsd:jed), source=0.0)
    allocate(CS%beta_dx2_h(isd:ied,jsd:jed), source=0.0)
    allocate(CS%f2_dx2_h(isd:ied,jsd:jed), source=0.0)
    do j=js-1,je+1 ; do i=is-1,ie+1
      CS%f2_dx2_h(i,j) = ((G%dxT(i,j)**2) + (G%dyT(i,j)**2)) * &
          max(0.25 * ((G%Coriolis2Bu(I,J) + G%Coriolis2Bu(I-1,J-1)) + &
                      (G%Coriolis2Bu(I-1,J) + G%Coriolis2Bu(I,J-1))), &
              absurdly_small_freq**2)
      CS%beta_dx2_h(i,j) = oneOrTwo * ((G%dxT(i,j)**2) + (G%dyT(i,j)**2)) * (sqrt(0.5 * &
          ( ((((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2) + &
             (((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2)) + &
            ((((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2) + &
             (((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2)) ) ))
    enddo ; enddo
  endif

  if (CS%calculate_cg1) then
    in_use = .true.
    allocate(CS%cg1(isd:ied,jsd:jed), source=0.0)
    call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
    call get_param(param_file, mdl, "REMAPPING_ANSWER_DATE", remap_answer_date, &
                 "The vintage of the expressions and order of arithmetic to use for remapping.  "//&
                 "Values below 20190101 result in the use of older, less accurate expressions "//&
                 "that were in use at the end of 2018.  Higher values result in the use of more "//&
                 "robust and accurate forms of mathematically equivalent expressions.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) remap_answer_date = max(remap_answer_date, 20230701)

    call get_param(param_file, mdl, "INTERNAL_WAVE_SPEED_TOL", wave_speed_tol, &
                 "The fractional tolerance for finding the wave speeds.", &
                 units="nondim", default=0.001)
    !### Set defaults so that wave_speed_min*wave_speed_tol >= 1e-9 m s-1
    call get_param(param_file, mdl, "INTERNAL_WAVE_SPEED_MIN", wave_speed_min, &
                 "A floor in the first mode speed below which 0 used instead.", &
                 units="m s-1", default=0.0, scale=US%m_s_to_L_T)
    call get_param(param_file, mdl, "INTERNAL_WAVE_SPEED_BETTER_EST", better_speed_est, &
                 "If true, use a more robust estimate of the first mode wave speed as the "//&
                 "starting point for iterations.", default=.true.)
    call get_param(param_file, mdl, "REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                   do_not_log=.true., default=.true.)
    call get_param(param_file, mdl, "EBT_REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 "If true, use the OM4 remapping-via-subcells algorithm for calculating EBT structure. "//&
                 "See REMAPPING_USE_OM4_SUBCELLS for details. "//&
                 "We recommend setting this option to false.", default=om4_remap_via_sub_cells)
    call wave_speed_init(CS%wave_speed, GV, use_ebt_mode=CS%Resoln_use_ebt, &
                         mono_N2_depth=N2_filter_depth, remap_answer_date=remap_answer_date, &
                         better_speed_est=better_speed_est, min_speed=wave_speed_min, &
                         om4_remap_via_sub_cells=om4_remap_via_sub_cells, wave_speed_tol=wave_speed_tol)
  endif

  ! Leith parameters
  call get_param(param_file, mdl, "USE_QG_LEITH_GM", CS%use_QG_Leith_GM, &
               "If true, use the QG Leith viscosity as the GM coefficient.", &
               default=.false.)

  if (CS%Use_QG_Leith_GM) then
    call get_param(param_file, mdl, "LEITH_LAP_CONST", Leith_Lap_const, &
               "The nondimensional Laplacian Leith constant, \n"//&
               "often set to 1.0", units="nondim", default=0.0)

    call get_param(param_file, mdl, "USE_BETA_IN_LEITH", CS%use_beta_in_QG_Leith, &
               "If true, include the beta term in the Leith nonlinear eddy viscosity.", &
               default=.true.)

    ALLOC_(CS%Laplac3_const_u(IsdB:IedB,jsd:jed)) ; CS%Laplac3_const_u(:,:) = 0.0
    ALLOC_(CS%Laplac3_const_v(isd:ied,JsdB:JedB)) ; CS%Laplac3_const_v(:,:) = 0.0
    ALLOC_(CS%KH_u_QG(IsdB:IedB,jsd:jed,GV%ke)) ; CS%KH_u_QG(:,:,:) = 0.0
    ALLOC_(CS%KH_v_QG(isd:ied,JsdB:JedB,GV%ke)) ; CS%KH_v_QG(:,:,:) = 0.0
    ! register diagnostics

    CS%id_KH_u_QG = register_diag_field('ocean_model', 'KH_u_QG', diag%axesCuL, Time, &
       'Horizontal viscosity from Leith QG, at u-points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
    CS%id_KH_v_QG = register_diag_field('ocean_model', 'KH_v_QG', diag%axesCvL, Time, &
       'Horizontal viscosity from Leith QG, at v-points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)

    do j=Jsq,Jeq+1 ; do I=is-1,Ieq
      ! Static factors in the Leith schemes
      grid_sp_u2 = G%dyCu(I,j)*G%dxCu(I,j)
      grid_sp_u3 = grid_sp_u2*sqrt(grid_sp_u2)
      CS%Laplac3_const_u(I,j) = Leith_Lap_const * grid_sp_u3
    enddo ; enddo
    do j=js-1,Jeq ; do I=Isq,Ieq+1
      ! Static factors in the Leith schemes
      grid_sp_v2 = G%dyCv(i,J)*G%dxCv(i,J)
      grid_sp_v3 = grid_sp_v2*sqrt(grid_sp_v2)
      CS%Laplac3_const_v(i,J) = Leith_Lap_const * grid_sp_v3
    enddo ; enddo

    if (.not. CS%use_stored_slopes) call MOM_error(FATAL, &
           "MOM_lateral_mixing_coeffs.F90, VarMix_init:"//&
           "USE_STORED_SLOPES must be True when using QG Leith.")
  endif

  ! Re-enable variable mixing if one of the schemes was enabled
  CS%use_variable_mixing = in_use .or. CS%use_variable_mixing
end subroutine VarMix_init

!> Destructor for VarMix control structure
subroutine VarMix_end(CS)
  type(VarMix_CS), intent(inout) :: CS

  if (CS%Resoln_use_ebt .or. CS%khth_use_ebt_struct .or. CS%kdgl90_use_ebt_struct &
      .or. CS%BS_EBT_power>0. .or. CS%khtr_use_ebt_struct) deallocate(CS%ebt_struct)
  if (allocated(CS%sqg_struct)) deallocate(CS%sqg_struct)
  if (allocated(CS%BS_struct)) deallocate(CS%BS_struct)
  if (CS%khth_use_ebt_struct .or. CS%khth_use_sqg_struct) deallocate(CS%khth_struct)
  if (CS%khtr_use_ebt_struct .or. CS%khtr_use_sqg_struct) deallocate(CS%khtr_struct)
  if (CS%kdgl90_use_ebt_struct .or. CS%kdgl90_use_sqg_struct) deallocate(CS%kdgl90_struct)

  if (CS%use_stored_slopes .or. CS%sqg_expo>0.0) then
    deallocate(CS%slope_x)
    deallocate(CS%slope_y)
  endif

  if (CS%calculate_Eady_growth_rate) then
    deallocate(CS%SN_u)
    deallocate(CS%SN_v)
  endif

  if (allocated(CS%L2u)) deallocate(CS%L2u)
  if (allocated(CS%L2v)) deallocate(CS%L2v)

  if (CS%Resoln_scaling_used) then
    deallocate(CS%Res_fn_h)
    deallocate(CS%Res_fn_q)
    deallocate(CS%Res_fn_u)
    deallocate(CS%Res_fn_v)
    deallocate(CS%beta_dx2_q)
    deallocate(CS%beta_dx2_u)
    deallocate(CS%beta_dx2_v)
    deallocate(CS%f2_dx2_q)
    deallocate(CS%f2_dx2_u)
    deallocate(CS%f2_dx2_v)
  endif

  if (CS%Depth_scaled_KhTh) then
    deallocate(CS%Depth_fn_u)
    deallocate(CS%Depth_fn_v)
  endif

  if (CS%calculate_Rd_dx) then
    deallocate(CS%Rd_dx_h)
    deallocate(CS%beta_dx2_h)
    deallocate(CS%f2_dx2_h)
  endif

  if (CS%calculate_cg1) then
    deallocate(CS%cg1)
  endif

  if (CS%Use_QG_Leith_GM) then
    DEALLOC_(CS%Laplac3_const_u)
    DEALLOC_(CS%Laplac3_const_v)
    DEALLOC_(CS%KH_u_QG)
    DEALLOC_(CS%KH_v_QG)
  endif
end subroutine VarMix_end

!> \namespace mom_lateral_mixing_coeffs
!!
!! This module provides a container for various factors used in prescribing diffusivities, that are
!! a function of the state (in particular the stratification and isoneutral slopes).
!!
!! \section section_Resolution_Function The resolution function
!!
!! The resolution function is expressed in terms of the ratio of grid-spacing to deformation radius.
!! The square of the resolution parameter is
!!
!! \f[
!! R^2 = \frac{L_d^2}{\Delta^2} = \frac{ c_g^2 }{ f^2 \Delta^2 + c_g \beta \Delta^2 }
!! \f]
!!
!! where the grid spacing is calculated as
!!
!! \f[
!! \Delta^2 = \Delta x^2 + \Delta y^2 .
!! \f]
!!
!! \todo Check this reference to Bob on/off paper.
!! The resolution function used in scaling diffusivities (\cite hallberg2013) is
!!
!! \f[
!! r(\Delta,L_d) = \frac{1}{1+(\alpha R)^p}
!! \f]
!!
!! The resolution function can be applied independently to thickness diffusion \(module mom_thickness_diffuse\),
!! tracer diffusion \(mom_tracer_hordiff\) lateral viscosity \(mom_hor_visc\).
!!
!! Robert Hallberg, 2013: Using a resolution function to regulate parameterizations of oceanic mesoscale eddy effects.
!! Ocean Modelling, 71, pp 92-103.  http://dx.doi.org/10.1016/j.ocemod.2013.08.007
!!
!! | Symbol                | Module parameter |
!! | ------                | --------------- |
!! | -                     | <code>USE_VARIABLE_MIXING</code> |
!! | -                     | <code>RESOLN_SCALED_KH</code> |
!! | -                     | <code>RESOLN_SCALED_KHTH</code> |
!! | -                     | <code>RESOLN_SCALED_KHTR</code> |
!! | \f$ \alpha \f$        | <code>KH_RES_SCALE_COEF</code> (for thickness and tracer diffusivity) |
!! | \f$ p \f$             | <code>KH_RES_FN_POWER</code> (for thickness and tracer diffusivity) |
!! | \f$ \alpha \f$        | <code>VISC_RES_SCALE_COEF</code> (for lateral viscosity) |
!! | \f$ p \f$             | <code>VISC_RES_FN_POWER</code> (for lateral viscosity) |
!! | -                     | <code>GILL_EQUATORIAL_LD</code> |
!!
!!
!!
!! \section section_Vicbeck Visbeck diffusivity
!!
!! This module also calculates factors used in setting the thickness diffusivity similar to a Visbeck et al., 1997,
!! scheme.  The factors are combined in mom_thickness_diffuse::thickness_diffuse but calculated in this module.
!!
!! \f[
!! \kappa_h = \alpha_s L_s^2 S N
!! \f]
!!
!! where \f$S\f$ is the magnitude of the isoneutral slope and \f$N\f$ is the Brunt-Vaisala frequency.
!!
!! Visbeck, Marshall, Haine and Spall, 1997: Specification of Eddy Transfer Coefficients in Coarse-Resolution
!! Ocean Circulation Models. J. Phys. Oceanogr. http://dx.doi.org/10.1175/1520-0485(1997)027%3C0381:SOETCI%3E2.0.CO;2
!!
!! | Symbol                | Module parameter |
!! | ------                | --------------- |
!! | -                     | <code>USE_VARIABLE_MIXING</code> |
!! | \f$ \alpha_s \f$      | <code>KHTH_SLOPE_CFF</code> (for mom_thickness_diffuse module)|
!! | \f$ \alpha_s \f$      | <code>KHTR_SLOPE_CFF</code> (for mom_tracer_hordiff module)|
!! | \f$ L_{s} \f$         | <code>VISBECK_L_SCALE</code> |
!! | \f$ S_{max} \f$       | <code>VISBECK_MAX_SLOPE</code> |
!!
!!
!! \section section_vertical_structure_khth Vertical structure function for KhTh
!!
!! The thickness diffusivity can be prescribed a vertical distribution with the shape of the equivalent barotropic
!! velocity mode.  The structure function is stored in the control structure for this module (varmix_cs) but is
!! calculated using subroutines in mom_wave_speed.
!!
!! | Symbol                | Module parameter |
!! | ------                | --------------- |
!! | -                     | <code>KHTH_USE_EBT_STRUCT</code> |

end module MOM_lateral_mixing_coeffs
