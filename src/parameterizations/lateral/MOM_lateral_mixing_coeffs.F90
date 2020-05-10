!> Variable mixing coefficients
module MOM_lateral_mixing_coeffs

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,     only : hchksum, uvchksum
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, post_data
use MOM_diag_mediator, only : diag_ctrl, time_type, query_averaging_enabled
use MOM_domains,       only : create_group_pass, do_group_pass
use MOM_domains,       only : group_pass_type, pass_var, pass_vector
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_interface_heights, only : find_eta
use MOM_isopycnal_slopes, only : calc_isoneutral_slopes
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_wave_speed, only : wave_speed, wave_speed_CS, wave_speed_init

implicit none ; private

#include <MOM_memory.h>

!> Variable mixing coefficients
type, public :: VarMix_CS
  logical :: use_variable_mixing  !< If true, use the variable mixing.
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
  logical :: calculate_cg1        !< If true, calls wave_speed() to calculate the first
                                  !! baroclinic wave speed and populate CS%cg1.
                                  !! This parameter is set depending on other parameters.
  logical :: calculate_Rd_dx      !< If true, calculates Rd/dx and populate CS%Rd_dx_h.
                                  !! This parameter is set depending on other parameters.
  logical :: calculate_res_fns    !< If true, calculate all the resolution factors.
                                  !! This parameter is set depending on other parameters.
  logical :: calculate_depth_fns !< If true, calculate all the depth factors.
                                  !! This parameter is set depending on other parameters.
  logical :: calculate_Eady_growth_rate !< If true, calculate all the Eady growth rate.
                                  !! This parameter is set depending on other parameters.
  real, dimension(:,:), pointer :: &
    SN_u => NULL(), &     !< S*N at u-points [T-1 ~> s-1]
    SN_v => NULL(), &     !< S*N at v-points [T-1 ~> s-1]
    L2u => NULL(), &      !< Length scale^2 at u-points [L2 ~> m2]
    L2v => NULL(), &      !< Length scale^2 at v-points [L2 ~> m2]
    cg1 => NULL(), &      !< The first baroclinic gravity wave speed [L T-1 ~> m s-1].
    Res_fn_h => NULL(), & !< Non-dimensional function of the ratio the first baroclinic
                          !! deformation radius to the grid spacing at h points [nondim].
    Res_fn_q => NULL(), & !< Non-dimensional function of the ratio the first baroclinic
                          !! deformation radius to the grid spacing at q points [nondim].
    Res_fn_u => NULL(), & !< Non-dimensional function of the ratio the first baroclinic
                          !! deformation radius to the grid spacing at u points [nondim].
    Res_fn_v => NULL(), & !< Non-dimensional function of the ratio the first baroclinic
                          !! deformation radius to the grid spacing at v points [nondim].
    Depth_fn_u => NULL(), & !< Non-dimensional function of the ratio of the depth to
                            !! a reference depth (maximum 1) at u points [nondim]
    Depth_fn_v => NULL(), & !< Non-dimensional function of the ratio of the depth to
                            !! a reference depth (maximum 1) at v points [nondim]
    beta_dx2_h => NULL(), & !< The magnitude of the gradient of the Coriolis parameter
                            !! times the grid spacing squared at h points [L T-1 ~> m s-1].
    beta_dx2_q => NULL(), & !< The magnitude of the gradient of the Coriolis parameter
                            !! times the grid spacing squared at q points [L T-1 ~> m s-1].
    beta_dx2_u => NULL(), & !< The magnitude of the gradient of the Coriolis parameter
                            !! times the grid spacing squared at u points [L T-1 ~> m s-1].
    beta_dx2_v => NULL(), & !< The magnitude of the gradient of the Coriolis parameter
                            !! times the grid spacing squared at v points [L T-1 ~> m s-1].
    f2_dx2_h => NULL(), & !< The Coriolis parameter squared times the grid
                          !! spacing squared at h [L2 T-2 ~> m2 s-2].
    f2_dx2_q => NULL(), & !< The Coriolis parameter squared times the grid
                          !! spacing squared at q [L2 T-2 ~> m2 s-2].
    f2_dx2_u => NULL(), & !< The Coriolis parameter squared times the grid
                          !! spacing squared at u [L2 T-2 ~> m2 s-2].
    f2_dx2_v => NULL(), & !< The Coriolis parameter squared times the grid
                          !! spacing squared at v [L2 T-2 ~> m2 s-2].
    Rd_dx_h => NULL()     !< Deformation radius over grid spacing [nondim]

  real, dimension(:,:,:), pointer :: &
    slope_x => NULL(), &  !< Zonal isopycnal slope [nondim]
    slope_y => NULL(), &  !< Meridional isopycnal slope [nondim]
    ebt_struct => NULL()  !< Vertical structure function to scale diffusivities with [nondim]
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
  real :: Visbeck_L_scale !< Fixed length scale in Visbeck formula
  real :: Res_coef_khth   !< A non-dimensional number that determines the function
                          !! of resolution, used for thickness and tracer mixing, as:
                          !!  F = 1 / (1 + (Res_coef_khth*Ld/dx)^Res_fn_power)
  real :: Res_coef_visc   !< A non-dimensional number that determines the function
                          !! of resolution, used for lateral viscosity, as:
                          !!  F = 1 / (1 + (Res_coef_visc*Ld/dx)^Res_fn_power)
  real :: depth_scaled_khth_h0 !< The depth above which KHTH is linearly scaled away [Z ~> m]
  real :: depth_scaled_khth_exp !< The exponent used in the depth dependent scaling function for KHTH [nondim]
  real :: kappa_smooth    !< A diffusivity for smoothing T/S in vanished layers [Z2 T-1 ~> m2 s-1]
  integer :: Res_fn_power_khth !< The power of dx/Ld in the KhTh resolution function.  Any
                               !! positive integer power may be used, but even powers
                               !! and especially 2 are coded to be more efficient.
  integer :: Res_fn_power_visc !< The power of dx/Ld in the Kh resolution function.  Any
                               !! positive integer power may be used, but even powers
                               !! and especially 2 are coded to be more efficient.
  real :: Visbeck_S_max   !< Upper bound on slope used in Eady growth rate [nondim].

  ! Leith parameters
  logical :: use_QG_Leith_GM      !< If true, uses the QG Leith viscosity as the GM coefficient
  logical :: use_beta_in_QG_Leith !< If true, includes the beta term in the QG Leith GM coefficient

  ! Diagnostics
  !>@{
  !! Diagnostic identifier
  integer :: id_SN_u=-1, id_SN_v=-1, id_L2u=-1, id_L2v=-1, id_Res_fn = -1
  integer :: id_N2_u=-1, id_N2_v=-1, id_S2_u=-1, id_S2_v=-1
  integer :: id_Rd_dx=-1, id_KH_u_QG = -1, id_KH_v_QG = -1
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.
  !>@}

  type(wave_speed_CS), pointer :: wave_speed_CSp => NULL() !< Wave speed control structure
  type(group_pass_type) :: pass_cg1 !< For group halo pass
  logical :: debug      !< If true, write out checksums of data for debugging
end type VarMix_CS

public VarMix_init, calc_slope_functions, calc_resoln_function
public calc_QG_Leith_viscosity, calc_depth_function

contains

!> Calculates the non-dimensional depth functions.
subroutine calc_depth_function(G, CS)
  type(ocean_grid_type),                    intent(in) :: G  !< Ocean grid structure
  type(VarMix_CS),                          pointer       :: CS !< Variable mixing coefficients

  ! Local variables
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: i, j
  real    :: H0 ! local variable for reference depth
  real    :: expo ! exponent used in the depth dependent scaling
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not. associated(CS)) call MOM_error(FATAL, "calc_depth_function:"// &
         "Module must be initialized before it is used.")
  if (.not. CS%calculate_depth_fns) return
  if (.not. associated(CS%Depth_fn_u)) call MOM_error(FATAL, &
    "calc_depth_function: %Depth_fn_u is not associated with Depth_scaled_KhTh.")
  if (.not. associated(CS%Depth_fn_v)) call MOM_error(FATAL, &
    "calc_depth_function: %Depth_fn_v is not associated with Depth_scaled_KhTh.")

  H0 = CS%depth_scaled_khth_h0
  expo = CS%depth_scaled_khth_exp
!$OMP do
  do j=js,je ; do I=is-1,Ieq
    CS%Depth_fn_u(I,j) = (MIN(1.0, 0.5*(G%bathyT(i,j) + G%bathyT(i+1,j))/H0))**expo
  enddo ; enddo
!$OMP do
  do J=js-1,Jeq ; do i=is,ie
    CS%Depth_fn_v(i,J) = (MIN(1.0, 0.5*(G%bathyT(i,j) + G%bathyT(i,j+1))/H0))**expo
  enddo ; enddo

end subroutine calc_depth_function

!> Calculates and stores the non-dimensional resolution functions
subroutine calc_resoln_function(h, tv, G, GV, US, CS)
  type(ocean_grid_type),                    intent(inout) :: G  !< Ocean grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                    intent(in)    :: tv !< Thermodynamic variables
  type(verticalGrid_type),                  intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),                    intent(in)    :: US !< A dimensional unit scaling type
  type(VarMix_CS),                          pointer       :: CS !< Variable mixing coefficients

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
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not. associated(CS)) call MOM_error(FATAL, "calc_resoln_function:"// &
         "Module must be initialized before it is used.")
  if (CS%calculate_cg1) then
    if (.not. associated(CS%cg1)) call MOM_error(FATAL, &
      "calc_resoln_function: %cg1 is not associated with Resoln_scaled_Kh.")
    if (CS%khth_use_ebt_struct) then
      if (.not. associated(CS%ebt_struct)) call MOM_error(FATAL, &
        "calc_resoln_function: %ebt_struct is not associated with RESOLN_USE_EBT.")
      if (CS%Resoln_use_ebt) then
        ! Both resolution fn and vertical structure are using EBT
        call wave_speed(h, tv, G, GV, US, CS%cg1, CS%wave_speed_CSp, modal_structure=CS%ebt_struct)
      else
        ! Use EBT to get vertical structure first and then re-calculate cg1 using first baroclinic mode
        call wave_speed(h, tv, G, GV, US, CS%cg1, CS%wave_speed_CSp, modal_structure=CS%ebt_struct, &
                        use_ebt_mode=.true.)
        call wave_speed(h, tv, G, GV, US, CS%cg1, CS%wave_speed_CSp)
      endif
      call pass_var(CS%ebt_struct, G%Domain)
    else
      call wave_speed(h, tv, G, GV, US, CS%cg1, CS%wave_speed_CSp)
    endif

    call create_group_pass(CS%pass_cg1, CS%cg1, G%Domain)
    call do_group_pass(CS%pass_cg1, G%Domain)
  endif

  ! Calculate and store the ratio between deformation radius and grid-spacing
  ! at h-points [nondim].
  if (CS%calculate_rd_dx) then
    if (.not. associated(CS%Rd_dx_h)) call MOM_error(FATAL, &
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

  if (.not. associated(CS%Res_fn_h)) call MOM_error(FATAL, &
    "calc_resoln_function: %Res_fn_h is not associated with Resoln_scaled_Kh.")
  if (.not. associated(CS%Res_fn_q)) call MOM_error(FATAL, &
    "calc_resoln_function: %Res_fn_q is not associated with Resoln_scaled_Kh.")
  if (.not. associated(CS%Res_fn_u)) call MOM_error(FATAL, &
    "calc_resoln_function: %Res_fn_u is not associated with Resoln_scaled_Kh.")
  if (.not. associated(CS%Res_fn_v)) call MOM_error(FATAL, &
    "calc_resoln_function: %Res_fn_v is not associated with Resoln_scaled_Kh.")
  if (.not. associated(CS%f2_dx2_h)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_h is not associated with Resoln_scaled_Kh.")
  if (.not. associated(CS%f2_dx2_q)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_q is not associated with Resoln_scaled_Kh.")
  if (.not. associated(CS%f2_dx2_u)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_u is not associated with Resoln_scaled_Kh.")
  if (.not. associated(CS%f2_dx2_v)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_v is not associated with Resoln_scaled_Kh.")
  if (.not. associated(CS%beta_dx2_h)) call MOM_error(FATAL, &
    "calc_resoln_function: %beta_dx2_h is not associated with Resoln_scaled_Kh.")
  if (.not. associated(CS%beta_dx2_q)) call MOM_error(FATAL, &
    "calc_resoln_function: %beta_dx2_q is not associated with Resoln_scaled_Kh.")
  if (.not. associated(CS%beta_dx2_u)) call MOM_error(FATAL, &
    "calc_resoln_function: %beta_dx2_u is not associated with Resoln_scaled_Kh.")
  if (.not. associated(CS%beta_dx2_v)) call MOM_error(FATAL, &
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
  endif

end subroutine calc_resoln_function

!> Calculates and stores functions of isopycnal slopes, e.g. Sx, Sy, S*N, mostly used in the Visbeck et al.
!! style scaling of diffusivity
subroutine calc_slope_functions(h, tv, dt, G, GV, US, CS)
  type(ocean_grid_type),                    intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                  intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),                    intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                    intent(in)    :: tv !< Thermodynamic variables
  real,                                     intent(in)    :: dt !< Time increment [T ~> s]
  type(VarMix_CS),                          pointer       :: CS !< Variable mixing coefficients
  ! Local variables
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)+1) :: &
    e             ! The interface heights relative to mean sea level [Z ~> m].
  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)+1) :: N2_u ! Square of Brunt-Vaisala freq at u-points [T-2 ~> s-2]
  real, dimension(SZI_(G), SZJB_(G), SZK_(G)+1) :: N2_v ! Square of Brunt-Vaisala freq at v-points [T-2 ~> s-2]

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_lateral_mixing_coeffs.F90, calc_slope_functions:"//&
         "Module must be initialized before it is used.")

  if (CS%calculate_Eady_growth_rate) then
    call find_eta(h, tv, G, GV, US, e, halo_size=2)
    if (CS%use_stored_slopes) then
      call calc_isoneutral_slopes(G, GV, US, h, e, tv, dt*CS%kappa_smooth, &
                                  CS%slope_x, CS%slope_y, N2_u, N2_v, 1)
      call calc_Visbeck_coeffs(h, CS%slope_x, CS%slope_y, N2_u, N2_v, G, GV, US, CS)
!     call calc_slope_functions_using_just_e(h, G, CS, e, .false.)
    else
      !call calc_isoneutral_slopes(G, GV, h, e, tv, dt*CS%kappa_smooth, CS%slope_x, CS%slope_y)
      call calc_slope_functions_using_just_e(h, G, GV, US, CS, e, .true.)
    endif
  endif

  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_SN_u > 0) call post_data(CS%id_SN_u, CS%SN_u, CS%diag)
    if (CS%id_SN_v > 0) call post_data(CS%id_SN_v, CS%SN_v, CS%diag)
    if (CS%id_L2u > 0)  call post_data(CS%id_L2u, CS%L2u, CS%diag)
    if (CS%id_L2v > 0)  call post_data(CS%id_L2v, CS%L2v, CS%diag)
    if (CS%calculate_Eady_growth_rate .and. CS%use_stored_slopes) then
      if (CS%id_N2_u > 0) call post_data(CS%id_N2_u, N2_u, CS%diag)
      if (CS%id_N2_v > 0) call post_data(CS%id_N2_v, N2_v, CS%diag)
    endif
  endif

end subroutine calc_slope_functions

!> Calculates factors used when setting diffusivity coefficients similar to Visbeck et al.
subroutine calc_Visbeck_coeffs(h, slope_x, slope_y, N2_u, N2_v, G, GV, US, CS)
  type(ocean_grid_type),                       intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                     intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),    intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1), intent(in)    :: slope_x !< Zonal isoneutral slope
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1), intent(in)    :: N2_u    !< Buoyancy (Brunt-Vaisala) frequency
                                                                        !! at u-points [T-2 ~> s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1), intent(in)    :: slope_y !< Meridional isoneutral slope
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1), intent(in)    :: N2_v    !< Buoyancy (Brunt-Vaisala) frequency
                                                                        !! at v-points [T-2 ~> s-2]
  type(unit_scale_type),                       intent(in)    :: US !< A dimensional unit scaling type
  type(VarMix_CS),                             pointer       :: CS !< Variable mixing coefficients

  ! Local variables
  real :: S2            ! Interface slope squared [nondim]
  real :: N2            ! Positive buoyancy frequency or zero [T-2 ~> s-2]
  real :: Hup, Hdn      ! Thickness from above, below [H ~> m or kg m-2]
  real :: H_geom        ! The geometric mean of Hup*Hdn [H ~> m or kg m-2].
  integer :: is, ie, js, je, nz
  integer :: i, j, k, kb_max
  real :: S2max, wNE, wSE, wSW, wNW
  real :: H_u(SZIB_(G)), H_v(SZI_(G))
  real :: S2_u(SZIB_(G), SZJ_(G))
  real :: S2_v(SZI_(G), SZJB_(G))

  if (.not. associated(CS)) call MOM_error(FATAL, "calc_slope_function:"// &
         "Module must be initialized before it is used.")
  if (.not. CS%calculate_Eady_growth_rate) return
  if (.not. associated(CS%SN_u)) call MOM_error(FATAL, "calc_slope_function:"// &
         "%SN_u is not associated with use_variable_mixing.")
  if (.not. associated(CS%SN_v)) call MOM_error(FATAL, "calc_slope_function:"// &
         "%SN_v is not associated with use_variable_mixing.")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  S2max = CS%Visbeck_S_max**2

  !$OMP parallel do default(shared)
  do j=js-1,je+1 ; do i=is-1,ie+1
    CS%SN_u(i,j) = 0.0
    CS%SN_v(i,j) = 0.0
  enddo ; enddo

  ! To set the length scale based on the deformation radius, use wave_speed to
  ! calculate the first-mode gravity wave speed and then blend the equatorial
  ! and midlatitude deformation radii, using calc_resoln_function as a template.

  !$OMP parallel do default(shared) private(S2,H_u,Hdn,Hup,H_geom,N2,wNE,wSE,wSW,wNW)
  do j = js,je
    do I=is-1,ie
      CS%SN_u(I,j) = 0. ; H_u(I) = 0. ; S2_u(I,j) = 0.
    enddo
    do K=2,nz ; do I=is-1,ie
      Hdn = sqrt( h(i,j,k) * h(i+1,j,k) )
      Hup = sqrt( h(i,j,k-1) * h(i+1,j,k-1) )
      H_geom = sqrt( Hdn * Hup )
     !H_geom = H_geom * sqrt(N2) ! WKB-ish
     !H_geom = H_geom * N2       ! WKB-ish
      wSE = G%mask2dCv(i+1,J-1) * ( (h(i+1,j,k)*h(i+1,j-1,k)) * (h(i+1,j,k-1)*h(i+1,j-1,k-1)) )
      wNW = G%mask2dCv(i  ,J  ) * ( (h(i  ,j,k)*h(i  ,j+1,k)) * (h(i  ,j,k-1)*h(i  ,j+1,k-1)) )
      wNE = G%mask2dCv(i+1,J  ) * ( (h(i+1,j,k)*h(i+1,j+1,k)) * (h(i+1,j,k-1)*h(i+1,j+1,k-1)) )
      wSW = G%mask2dCv(i  ,J-1) * ( (h(i  ,j,k)*h(i  ,j-1,k)) * (h(i  ,j,k-1)*h(i  ,j-1,k-1)) )
      S2 =  slope_x(I,j,K)**2 + &
              ((wNW*slope_y(i,J,K)**2 + wSE*slope_y(i+1,J-1,K)**2) + &
               (wNE*slope_y(i+1,J,K)**2 + wSW*slope_y(i,J-1,K)**2) ) / &
              ( ((wSE+wNW) + (wNE+wSW)) + GV%H_subroundoff**4 )
      if (S2max>0.) S2 = S2 * S2max / (S2 + S2max) ! Limit S2

      N2 = max(0., N2_u(I,j,k))
      CS%SN_u(I,j) = CS%SN_u(I,j) + sqrt( S2*N2 )*H_geom
      S2_u(I,j) = S2_u(I,j) + S2*H_geom
      H_u(I) = H_u(I) + H_geom
    enddo ; enddo
    do I=is-1,ie
      if (H_u(I)>0.) then
        CS%SN_u(I,j) = G%mask2dCu(I,j) * CS%SN_u(I,j) / H_u(I)
        S2_u(I,j) =  G%mask2dCu(I,j) * S2_u(I,j) / H_u(I)
      else
        CS%SN_u(I,j) = 0.
      endif
    enddo
  enddo

  !$OMP parallel do default(shared) private(S2,H_v,Hdn,Hup,H_geom,N2,wNE,wSE,wSW,wNW)
  do J = js-1,je
    do i=is,ie
      CS%SN_v(i,J) = 0. ; H_v(i) = 0. ; S2_v(i,J) = 0.
    enddo
    do K=2,nz ; do i=is,ie
      Hdn = sqrt( h(i,j,k) * h(i,j+1,k) )
      Hup = sqrt( h(i,j,k-1) * h(i,j+1,k-1) )
      H_geom = sqrt( Hdn * Hup )
     !H_geom = H_geom * sqrt(N2) ! WKB-ish
     !H_geom = H_geom * N2       ! WKB-ish
      wSE = G%mask2dCu(I,j)     * ( (h(i,j  ,k)*h(i+1,j  ,k)) * (h(i,j  ,k-1)*h(i+1,j  ,k-1)) )
      wNW = G%mask2dCu(I-1,j+1) * ( (h(i,j+1,k)*h(i-1,j+1,k)) * (h(i,j+1,k-1)*h(i-1,j+1,k-1)) )
      wNE = G%mask2dCu(I,j+1)   * ( (h(i,j+1,k)*h(i+1,j+1,k)) * (h(i,j+1,k-1)*h(i+1,j+1,k-1)) )
      wSW = G%mask2dCu(I-1,j)   * ( (h(i,j  ,k)*h(i-1,j  ,k)) * (h(i,j  ,k-1)*h(i-1,j  ,k-1)) )
      S2 = slope_y(i,J,K)**2 + &
             ((wSE*slope_x(I,j,K)**2 + wNW*slope_x(I-1,j+1,K)**2) + &
              (wNE*slope_x(I,j+1,K)**2 + wSW*slope_x(I-1,j,K)**2) ) / &
             ( ((wSE+wNW) + (wNE+wSW)) + GV%H_subroundoff**4 )
      if (S2max>0.) S2 = S2 * S2max / (S2 + S2max) ! Limit S2

      N2 = max(0., N2_v(i,J,K))
      CS%SN_v(i,J) = CS%SN_v(i,J) + sqrt( S2*N2 )*H_geom
      S2_v(i,J) = S2_v(i,J) + S2*H_geom
      H_v(i) = H_v(i) + H_geom
    enddo ; enddo
    do i=is,ie
      if (H_v(i)>0.) then
        CS%SN_v(i,J) = G%mask2dCv(i,J) * CS%SN_v(i,J) / H_v(i)
        S2_v(i,J) = G%mask2dCv(i,J) * S2_v(i,J) / H_v(i)
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
    call uvchksum("calc_Visbeck_coeffs slope_[xy]", slope_x, slope_y, G%HI, haloshift=1)
    call uvchksum("calc_Visbeck_coeffs N2_u, N2_v", N2_u, N2_v, G%HI, &
                  scale=US%s_to_T**2, scalar_pair=.true.)
    call uvchksum("calc_Visbeck_coeffs SN_[uv]", CS%SN_u, CS%SN_v, G%HI, &
                  scale=US%s_to_T, scalar_pair=.true.)
  endif

end subroutine calc_Visbeck_coeffs

!> The original calc_slope_function() that calculated slopes using
!! interface positions only, not accounting for density variations.
subroutine calc_slope_functions_using_just_e(h, G, GV, US, CS, e, calculate_slopes)
  type(ocean_grid_type),                      intent(inout) :: G  !< Ocean grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  type(verticalGrid_type),                    intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US !< A dimensional unit scaling type
  type(VarMix_CS),                            pointer       :: CS !< Variable mixing coefficients
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(in)    :: e  !< Interface position [Z ~> m]
  logical,                                    intent(in)    :: calculate_slopes !< If true, calculate slopes internally
                                                                  !! otherwise use slopes stored in CS
  ! Local variables
  real :: E_x(SZIB_(G), SZJ_(G))  ! X-slope of interface at u points [nondim] (for diagnostics)
  real :: E_y(SZI_(G), SZJB_(G))  ! Y-slope of interface at v points [nondim] (for diagnostics)
  real :: H_cutoff      ! Local estimate of a minimum thickness for masking [H ~> m or kg m-2]
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: S2            ! Interface slope squared [nondim]
  real :: N2            ! Brunt-Vaisala frequency squared [T-2 ~> s-2]
  real :: Hup, Hdn      ! Thickness from above, below [H ~> m or kg m-2]
  real :: H_geom        ! The geometric mean of Hup*Hdn [H ~> m or kg m-2].
  real :: one_meter     ! One meter in thickness units [H ~> m or kg m-2].
  integer :: is, ie, js, je, nz
  integer :: i, j, k, kb_max
  real    :: S2N2_u_local(SZIB_(G), SZJ_(G),SZK_(G))
  real    :: S2N2_v_local(SZI_(G), SZJB_(G),SZK_(G))

  if (.not. associated(CS)) call MOM_error(FATAL, "calc_slope_function:"// &
         "Module must be initialized before it is used.")
  if (.not. CS%calculate_Eady_growth_rate) return
  if (.not. associated(CS%SN_u)) call MOM_error(FATAL, "calc_slope_function:"// &
         "%SN_u is not associated with use_variable_mixing.")
  if (.not. associated(CS%SN_v)) call MOM_error(FATAL, "calc_slope_function:"// &
         "%SN_v is not associated with use_variable_mixing.")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  one_meter = 1.0 * GV%m_to_H
  h_neglect = GV%H_subroundoff
  H_cutoff = real(2*nz) * (GV%Angstrom_H + h_neglect)

  ! To set the length scale based on the deformation radius, use wave_speed to
  ! calculate the first-mode gravity wave speed and then blend the equatorial
  ! and midlatitude deformation radii, using calc_resoln_function as a template.

  !$OMP parallel do default(shared) private(E_x,E_y,S2,Hdn,Hup,H_geom,N2)
  do k=nz,CS%VarMix_Ktop,-1

    if (calculate_slopes) then
      ! Calculate the interface slopes E_x and E_y and u- and v- points respectively
      do j=js-1,je+1 ; do I=is-1,ie
        E_x(I,j) = US%Z_to_L*(e(i+1,j,K)-e(i,j,K))*G%IdxCu(I,j)
        ! Mask slopes where interface intersects topography
        if (min(h(I,j,k),h(I+1,j,k)) < H_cutoff) E_x(I,j) = 0.
      enddo ; enddo
      do J=js-1,je ; do i=is-1,ie+1
        E_y(i,J) = US%Z_to_L*(e(i,j+1,K)-e(i,j,K))*G%IdyCv(i,J)
        ! Mask slopes where interface intersects topography
        if (min(h(i,J,k),h(i,J+1,k)) < H_cutoff) E_y(I,j) = 0.
      enddo ; enddo
    else
      do j=js-1,je+1 ; do I=is-1,ie
        E_x(I,j) = CS%slope_x(I,j,k)
        if (min(h(I,j,k),h(I+1,j,k)) < H_cutoff) E_x(I,j) = 0.
      enddo ; enddo
      do j=js-1,je ; do I=is-1,ie+1
        E_y(i,J) = CS%slope_y(i,J,k)
        if (min(h(i,J,k),h(i,J+1,k)) < H_cutoff) E_y(I,j) = 0.
      enddo ; enddo
    endif

    ! Calculate N*S*h from this layer and add to the sum
    do j=js,je ; do I=is-1,ie
      S2 = ( E_x(I,j)**2  + 0.25*( &
            (E_y(I,j)**2+E_y(I+1,j-1)**2)+(E_y(I+1,j)**2+E_y(I,j-1)**2) ) )
      Hdn = 2.*h(i,j,k)*h(i,j,k-1) / (h(i,j,k) + h(i,j,k-1) + h_neglect)
      Hup = 2.*h(i+1,j,k)*h(i+1,j,k-1) / (h(i+1,j,k) + h(i+1,j,k-1) + h_neglect)
      H_geom = sqrt(Hdn*Hup)
      N2 = GV%g_prime(k)*US%L_to_Z**2 / (GV%H_to_Z * max(Hdn,Hup,one_meter))
      if (min(h(i,j,k-1), h(i+1,j,k-1), h(i,j,k), h(i+1,j,k)) < H_cutoff) &
        S2 = 0.0
      S2N2_u_local(I,j,k) = (H_geom * GV%H_to_Z) * S2 * N2
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      S2 = ( E_y(i,J)**2  + 0.25*( &
            (E_x(i,J)**2+E_x(i-1,J+1)**2)+(E_x(i,J+1)**2+E_x(i-1,J)**2) ) )
      Hdn = 2.*h(i,j,k)*h(i,j,k-1) / (h(i,j,k) + h(i,j,k-1) + h_neglect)
      Hup = 2.*h(i,j+1,k)*h(i,j+1,k-1) / (h(i,j+1,k) + h(i,j+1,k-1) + h_neglect)
      H_geom = sqrt(Hdn*Hup)
      N2 = GV%g_prime(k)*US%L_to_Z**2 / (GV%H_to_Z * max(Hdn,Hup,one_meter))
      if (min(h(i,j,k-1), h(i,j+1,k-1), h(i,j,k), h(i,j+1,k)) < H_cutoff) &
        S2 = 0.0
      S2N2_v_local(i,J,k) = (H_geom * GV%H_to_Z) * S2 * N2
    enddo ; enddo

  enddo ! k
  !$OMP parallel do default(shared)
  do j=js,je
    do I=is-1,ie ; CS%SN_u(I,j) = 0.0 ; enddo
    do k=nz,CS%VarMix_Ktop,-1 ; do I=is-1,ie
      CS%SN_u(I,j) = CS%SN_u(I,j) + S2N2_u_local(I,j,k)
    enddo ; enddo
    ! SN above contains S^2*N^2*H, convert to vertical average of S*N
    do I=is-1,ie
      !SN_u(I,j) = sqrt( SN_u(I,j) / ( max(G%bathyT(I,j), G%bathyT(I+1,j)) + GV%Angstrom_Z ) ))
      !The code below behaves better than the line above. Not sure why? AJA
      if ( min(G%bathyT(I,j), G%bathyT(I+1,j)) > H_cutoff*GV%H_to_Z ) then
        CS%SN_u(I,j) = G%mask2dCu(I,j) * sqrt( CS%SN_u(I,j) / &
                                               (max(G%bathyT(I,j), G%bathyT(I+1,j))) )
      else
        CS%SN_u(I,j) = 0.0
      endif
    enddo
  enddo
  !$OMP parallel do default(shared)
  do J=js-1,je
    do i=is,ie ; CS%SN_v(i,J) = 0.0 ; enddo
    do k=nz,CS%VarMix_Ktop,-1 ; do i=is,ie
      CS%SN_v(i,J) = CS%SN_v(i,J) + S2N2_v_local(i,J,k)
    enddo ; enddo
    do i=is,ie
      !SN_v(i,J) = sqrt( SN_v(i,J) / ( max(G%bathyT(i,J), G%bathyT(i,J+1)) + GV%Angstrom_Z ) ))
      !The code below behaves better than the line above. Not sure why? AJA
      if ( min(G%bathyT(I,j), G%bathyT(I+1,j)) > H_cutoff*GV%H_to_Z ) then
        CS%SN_v(i,J) = G%mask2dCv(i,J) * sqrt( CS%SN_v(i,J) / &
                                               (max(G%bathyT(i,J), G%bathyT(i,J+1))) )
      else
        CS%SN_v(I,j) = 0.0
      endif
    enddo
  enddo

end subroutine calc_slope_functions_using_just_e

!> Calculates the Leith Laplacian and bi-harmonic viscosity coefficients
subroutine calc_QG_Leith_viscosity(CS, G, GV, US, h, k, div_xx_dx, div_xx_dy, vort_xy_dx, vort_xy_dy)
  type(VarMix_CS),                           pointer     :: CS !< Variable mixing coefficients
  type(ocean_grid_type),                     intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),                     intent(in)  :: US   !< A dimensional unit scaling type
! real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)  :: u  !< Zonal flow [L T-1 ~> m s-1]
! real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)  :: v  !< Meridional flow [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h !< Layer thickness [H ~> m or kg m-2]
  integer,                                   intent(in)  :: k  !< Layer for which to calculate vorticity magnitude
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)  :: div_xx_dx  !< x-derivative of horizontal divergence
                                                                 !! (d/dx(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
  real, dimension(SZI_(G),SZJB_(G)),         intent(in)  :: div_xx_dy  !< y-derivative of horizontal divergence
                                                                 !! (d/dy(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
  real, dimension(SZI_(G),SZJB_(G)),         intent(inout) :: vort_xy_dx !< x-derivative of vertical vorticity
                                                                 !! (d/dx(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
  real, dimension(SZIB_(G),SZJ_(G)),         intent(inout) :: vort_xy_dy !< y-derivative of vertical vorticity
                                                                 !! (d/dy(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
!  real, dimension(SZI_(G),SZJ_(G)),          intent(out) :: Leith_Kh_h !< Leith Laplacian viscosity
                                                                 !! at h-points [L2 T-1 ~> m2 s-1]
!  real, dimension(SZIB_(G),SZJB_(G)),        intent(out) :: Leith_Kh_q !< Leith Laplacian viscosity
                                                                 !! at q-points [L2 T-1 ~> m2 s-1]
!  real, dimension(SZI_(G),SZJ_(G)),          intent(out) :: Leith_Ah_h !< Leith bi-harmonic viscosity
                                                                 !! at h-points [L4 T-1 ~> m4 s-1]
!  real, dimension(SZIB_(G),SZJB_(G)),        intent(out) :: Leith_Ah_q !< Leith bi-harmonic viscosity
                                                                 !! at q-points [L4 T-1 ~> m4 s-1]

  ! Local variables
  real, dimension(SZI_(G),SZJB_(G)) :: &
    dslopey_dz, & ! z-derivative of y-slope at v-points [Z-1 ~> m-1]
    h_at_v,     & ! Thickness at v-points [H ~> m or kg m-2]
    beta_v,     & ! Beta at v-points [T-1 L-1 ~> s-1 m-1]
    grad_vort_mag_v, & ! Magnitude of vorticity gradient at v-points [T-1 L-1 ~> s-1 m-1]
    grad_div_mag_v     ! Magnitude of divergence gradient at v-points [T-1 L-1 ~> s-1 m-1]

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    dslopex_dz, & ! z-derivative of x-slope at u-points [Z-1 ~> m-1]
    h_at_u,     & ! Thickness at u-points [H ~> m or kg m-2]
    beta_u,     & ! Beta at u-points [T-1 L-1 ~> s-1 m-1]
    grad_vort_mag_u, & ! Magnitude of vorticity gradient at u-points [T-1 L-1 ~> s-1 m-1]
    grad_div_mag_u     ! Magnitude of divergence gradient at u-points [T-1 L-1 ~> s-1 m-1]
  real :: h_at_slope_above ! The thickness above [H ~> m or kg m-2]
  real :: h_at_slope_below ! The thickness below [H ~> m or kg m-2]
  real :: Ih ! The inverse of a combination of thicknesses [H-1 ~> m-1 or m2 kg-1]
  real :: f  ! A copy of the Coriolis parameter [T-1 ~> s-1]
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq,nz
  real :: inv_PI3

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nz = G%ke

  inv_PI3 = 1.0/((4.0*atan(1.0))**3)

  !### I believe this halo update to be unnecessary. -RWH
  call pass_var(h, G%Domain)

  if ((k > 1) .and. (k < nz)) then

  ! Add in stretching term for the QG Leith vsicosity
!  if (CS%use_QG_Leith) then

    !### do j=js-1,je+1 ; do I=is-2,Ieq+1
    do j=js-2,Jeq+2 ; do I=is-2,Ieq+1
      h_at_slope_above = 2. * ( h(i,j,k-1) * h(i+1,j,k-1) ) * ( h(i,j,k) * h(i+1,j,k) ) / &
                         ( ( h(i,j,k-1) * h(i+1,j,k-1) ) * ( h(i,j,k) + h(i+1,j,k) ) &
                         + ( h(i,j,k) * h(i+1,j,k) ) * ( h(i,j,k-1) + h(i+1,j,k-1) ) + GV%H_subroundoff )
      h_at_slope_below = 2. * ( h(i,j,k) * h(i+1,j,k) ) * ( h(i,j,k+1) * h(i+1,j,k+1) ) / &
                         ( ( h(i,j,k) * h(i+1,j,k) ) * ( h(i,j,k+1) + h(i+1,j,k+1) ) &
                         + ( h(i,j,k+1) * h(i+1,j,k+1) ) * ( h(i,j,k) + h(i+1,j,k) ) + GV%H_subroundoff )
      Ih = 1./ ( ( h_at_slope_above + h_at_slope_below + GV%H_subroundoff ) * GV%H_to_Z )
      dslopex_dz(I,j) = 2. * ( CS%slope_x(i,j,k) - CS%slope_x(i,j,k+1) ) * Ih
      h_at_u(I,j) = 2. * ( h_at_slope_above * h_at_slope_below ) * Ih
    enddo ; enddo

    !### do J=js-2,Jeq+1 ; do i=is-1,ie+1
    do J=js-2,Jeq+1 ; do i=is-2,Ieq+2
      h_at_slope_above = 2. * ( h(i,j,k-1) * h(i,j+1,k-1) ) * ( h(i,j,k) * h(i,j+1,k) ) / &
                         ( ( h(i,j,k-1) * h(i,j+1,k-1) ) * ( h(i,j,k) + h(i,j+1,k) ) &
                         + ( h(i,j,k) * h(i,j+1,k) ) * ( h(i,j,k-1) + h(i,j+1,k-1) ) + GV%H_subroundoff )
      h_at_slope_below = 2. * ( h(i,j,k) * h(i,j+1,k) ) * ( h(i,j,k+1) * h(i,j+1,k+1) ) / &
                         ( ( h(i,j,k) * h(i,j+1,k) ) * ( h(i,j,k+1) + h(i,j+1,k+1) ) &
                         + ( h(i,j,k+1) * h(i,j+1,k+1) ) * ( h(i,j,k) + h(i,j+1,k) ) + GV%H_subroundoff )
      Ih = 1./ ( ( h_at_slope_above + h_at_slope_below + GV%H_subroundoff ) * GV%H_to_Z )
      dslopey_dz(i,J) = 2. * ( CS%slope_y(i,j,k) - CS%slope_y(i,j,k+1) ) * Ih
      h_at_v(i,J) = 2. * ( h_at_slope_above * h_at_slope_below ) * Ih
    enddo ; enddo

    !### do J=js-1,je ; do i=is-1,Ieq+1
    do J=js-2,Jeq+1 ; do i=is-1,Ieq+1
      f = 0.5 * ( G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J) )
      vort_xy_dx(i,J) = vort_xy_dx(i,J) - f * US%L_to_Z * &
            ( ( h_at_u(I,j) * dslopex_dz(I,j) + h_at_u(I-1,j+1) * dslopex_dz(I-1,j+1) ) &
            + ( h_at_u(I-1,j) * dslopex_dz(I-1,j) + h_at_u(I,j+1) * dslopex_dz(I,j+1) ) ) / &
              ( ( h_at_u(I,j) + h_at_u(I-1,j+1) ) + ( h_at_u(I-1,j) + h_at_u(I,j+1) ) + GV%H_subroundoff)
    enddo ; enddo

    !### do j=js-1,Jeq+1 ; do I=is-1,ie
    do j=js-1,Jeq+1 ; do I=is-2,Ieq+1
      f = 0.5 * ( G%CoriolisBu(I,J) + G%CoriolisBu(I,J-1) )
      !### I think that this should be vort_xy_dy(I,j) = vort_xy_dy(I,j) - f * &
      vort_xy_dy(I,j) = vort_xy_dx(I,j) - f * US%L_to_Z * &
            ( ( h_at_v(i,J) * dslopey_dz(i,J) + h_at_v(i+1,J-1) * dslopey_dz(i+1,J-1) ) &
            + ( h_at_v(i,J-1) * dslopey_dz(i,J-1) + h_at_v(i+1,J) * dslopey_dz(i+1,J) ) ) / &
              ( ( h_at_v(i,J) + h_at_v(i+1,J-1) ) + ( h_at_v(i,J-1) + h_at_v(i+1,J) ) + GV%H_subroundoff)
    enddo ; enddo
  endif ! k > 1

  !### I believe this halo update to be unnecessary. -RWH
  call pass_vector(vort_xy_dy,vort_xy_dx,G%Domain)

  if (CS%use_QG_Leith_GM) then

    do j=js,je ; do I=is-1,Ieq
      !### These expressions are not rotationally symmetric.  Add parentheses and regroup, as in:
    ! grad_vort_mag_u(I,j) = SQRT(vort_xy_dy(I,j)**2 + (0.25*((vort_xy_dx(i,J) + vort_xy_dx(i+1,J-1)) +
    !                                                         (vort_xy_dx(i+1,J) + vort_xy_dx(i,J-1))))**2 )
      grad_vort_mag_u(I,j) = SQRT(vort_xy_dy(I,j)**2  + (0.25*(vort_xy_dx(i,J) + vort_xy_dx(i+1,J) &
                                                             + vort_xy_dx(i,J-1) + vort_xy_dx(i+1,J-1)))**2)
      grad_div_mag_u(I,j) = SQRT(div_xx_dx(I,j)**2  + (0.25*(div_xx_dy(i,J) + div_xx_dy(i+1,J) &
                                                           + div_xx_dy(i,J-1) + div_xx_dy(i+1,J-1)))**2)
      if (CS%use_beta_in_QG_Leith) then
        beta_u(I,j) = sqrt( (0.5*(G%dF_dx(i,j)+G%dF_dx(i+1,j))**2) + &
                                      (0.5*(G%dF_dy(i,j)+G%dF_dy(i+1,j))**2) )
        CS%KH_u_QG(I,j,k) = MIN(grad_vort_mag_u(I,j) + grad_div_mag_u(I,j), 3.0*beta_u(I,j)) * &
                            CS%Laplac3_const_u(I,j) * inv_PI3
      else
        CS%KH_u_QG(I,j,k) = (grad_vort_mag_u(I,j) + grad_div_mag_u(I,j)) * &
                            CS%Laplac3_const_u(I,j) * inv_PI3
      endif
    enddo ; enddo

    do J=js-1,Jeq ; do i=is,ie
      !### These expressions are not rotationally symmetric.  Add parentheses and regroup.
      grad_vort_mag_v(i,J) = SQRT(vort_xy_dx(i,J)**2  + (0.25*(vort_xy_dy(I,j) + vort_xy_dy(I-1,j) &
                                                             + vort_xy_dy(I,j+1) + vort_xy_dy(I-1,j+1)))**2)
      grad_div_mag_v(i,J) = SQRT(div_xx_dy(i,J)**2  + (0.25*(div_xx_dx(I,j) + div_xx_dx(I-1,j) &
                                                           + div_xx_dx(I,j+1) + div_xx_dx(I-1,j+1)))**2)
      if (CS%use_beta_in_QG_Leith) then
        beta_v(i,J) = sqrt( (0.5*(G%dF_dx(i,j)+G%dF_dx(i,j+1))**2) + &
                            (0.5*(G%dF_dy(i,j)+G%dF_dy(i,j+1))**2) )
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
  type(VarMix_CS),               pointer :: CS   !< Variable mixing coefficients
  ! Local variables
  real :: KhTr_Slope_Cff, KhTh_Slope_Cff, oneOrTwo
  real :: N2_filter_depth  ! A depth below which stratification is treated as monotonic when
                           ! calculating the first-mode wave speed [Z ~> m]
  real :: KhTr_passivity_coeff
  real :: absurdly_small_freq  ! A miniscule frequency that is used to avoid division by 0 [T-1 ~> s-1].  The
             ! default value is roughly (pi / (the age of the universe)).
  logical :: Gill_equatorial_Ld, use_FGNV_streamfn, use_MEKE, in_use
  logical :: default_2018_answers, remap_answers_2018
  real :: MLE_front_length
  real :: Leith_Lap_const      ! The non-dimensional coefficient in the Leith viscosity
  real :: grid_sp_u2, grid_sp_v2 ! Intermediate quantities for Leith metrics [L2 ~> m2]
  real :: grid_sp_u3, grid_sp_v3 ! Intermediate quantities for Leith metrics [L3 ~> m3]
  real :: wave_speed_min      ! A floor in the first mode speed below which 0 is returned [L T-1 ~> m s-1]
  real :: wave_speed_tol      ! The fractional tolerance for finding the wave speeds [nondim]
  logical :: better_speed_est ! If true, use a more robust estimate of the first
                              ! mode wave speed as the starting point for iterations.
! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_lateral_mixing_coeffs" ! This module's name.
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, i, j
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "VarMix_init called with an associated "// &
                             "control structure.")
    return
  endif

  allocate(CS)
  in_use = .false. ! Set to true to avoid deallocating
  CS%diag => diag ! Diagnostics pointer
  CS%calculate_cg1 = .false.
  CS%calculate_Rd_dx = .false.
  CS%calculate_res_fns = .false.
  CS%calculate_Eady_growth_rate = .false.
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
  call get_param(param_file, mdl, "RESOLN_USE_EBT", CS%Resoln_use_ebt, &
                 "If true, uses the equivalent barotropic wave speed instead "//&
                 "of first baroclinic wave for calculating the resolution fn.",&
                 default=.false.)
  call get_param(param_file, mdl, "KHTH_USE_EBT_STRUCT", CS%khth_use_ebt_struct, &
                 "If true, uses the equivalent barotropic structure "//&
                 "as the vertical structure of thickness diffusivity.",&
                 default=.false.)
  call get_param(param_file, mdl, "KHTH_SLOPE_CFF", KhTh_Slope_Cff, &
                 "The nondimensional coefficient in the Visbeck formula "//&
                 "for the interface depth diffusivity", units="nondim", &
                 default=0.0)
  call get_param(param_file, mdl, "KHTR_SLOPE_CFF", KhTr_Slope_Cff, &
                 "The nondimensional coefficient in the Visbeck formula "//&
                 "for the epipycnal tracer diffusivity", units="nondim", &
                 default=0.0)
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
  CS%calculate_cg1 = CS%calculate_cg1 .or. use_FGNV_streamfn
  call get_param(param_file, mdl, "USE_MEKE", use_MEKE, &
                 default=.false., do_not_log=.true.)
  CS%calculate_Rd_dx = CS%calculate_Rd_dx .or. use_MEKE
  CS%calculate_Eady_growth_rate = CS%calculate_Eady_growth_rate .or. use_MEKE
  call get_param(param_file, mdl, "KHTR_PASSIVITY_COEFF", KhTr_passivity_coeff, &
                 default=0., do_not_log=.true.)
  CS%calculate_Rd_dx = CS%calculate_Rd_dx .or. (KhTr_passivity_coeff>0.)
  call get_param(param_file, mdl, "MLE_FRONT_LENGTH", MLE_front_length, &
                 default=0., do_not_log=.true.)
  CS%calculate_Rd_dx = CS%calculate_Rd_dx .or. (MLE_front_length>0.)

  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false., do_not_log=.true.)


  if (CS%Resoln_use_ebt .or. CS%khth_use_ebt_struct) then
    in_use = .true.
    call get_param(param_file, mdl, "RESOLN_N2_FILTER_DEPTH", N2_filter_depth, &
                 "The depth below which N2 is monotonized to avoid stratification "//&
                 "artifacts from altering the equivalent barotropic mode structure.",&
                 units="m", default=2000., scale=US%m_to_Z)
    allocate(CS%ebt_struct(isd:ied,jsd:jed,G%ke)) ; CS%ebt_struct(:,:,:) = 0.0
  endif

  if (KhTr_Slope_Cff>0. .or. KhTh_Slope_Cff>0.) then
    CS%calculate_Eady_growth_rate = .true.
    call get_param(param_file, mdl, "VISBECK_MAX_SLOPE", CS%Visbeck_S_max, &
          "If non-zero, is an upper bound on slopes used in the "//&
          "Visbeck formula for diffusivity. This does not affect the "//&
          "isopycnal slope calculation used within thickness diffusion.",  &
          units="nondim", default=0.0)
  endif

  if (CS%use_stored_slopes) then
    in_use = .true.
    allocate(CS%slope_x(IsdB:IedB,jsd:jed,G%ke+1)) ; CS%slope_x(:,:,:) = 0.0
    allocate(CS%slope_y(isd:ied,JsdB:JedB,G%ke+1)) ; CS%slope_y(:,:,:) = 0.0
    call get_param(param_file, mdl, "KD_SMOOTH", CS%kappa_smooth, &
                 "A diapycnal diffusivity that is used to interpolate "//&
                 "more sensible values of T & S into thin layers.", &
                 units="m2 s-1", default=1.0e-6, scale=US%m_to_Z**2*US%T_to_s)
  endif

  if (CS%calculate_Eady_growth_rate) then
    in_use = .true.
    allocate(CS%SN_u(IsdB:IedB,jsd:jed)) ; CS%SN_u(:,:) = 0.0
    allocate(CS%SN_v(isd:ied,JsdB:JedB)) ; CS%SN_v(:,:) = 0.0
    CS%id_SN_u = register_diag_field('ocean_model', 'SN_u', diag%axesCu1, Time, &
       'Inverse eddy time-scale, S*N, at u-points', 's-1', conversion=US%s_to_T)
    CS%id_SN_v = register_diag_field('ocean_model', 'SN_v', diag%axesCv1, Time, &
       'Inverse eddy time-scale, S*N, at v-points', 's-1', conversion=US%s_to_T)
    call get_param(param_file, mdl, "VARMIX_KTOP", CS%VarMix_Ktop, &
                 "The layer number at which to start vertical integration "//&
                 "of S*N for purposes of finding the Eady growth rate.", &
                 units="nondim", default=2)
  endif

  if (KhTr_Slope_Cff>0. .or. KhTh_Slope_Cff>0.) then
    in_use = .true.
    call get_param(param_file, mdl, "VISBECK_L_SCALE", CS%Visbeck_L_scale, &
                 "The fixed length scale in the Visbeck formula.", units="m", &
                 default=0.0)
    allocate(CS%L2u(IsdB:IedB,jsd:jed)) ; CS%L2u(:,:) = 0.0
    allocate(CS%L2v(isd:ied,JsdB:JedB)) ; CS%L2v(:,:) = 0.0
    if (CS%Visbeck_L_scale<0) then
      do j=js,je ; do I=is-1,Ieq
        CS%L2u(I,j) = CS%Visbeck_L_scale**2 * G%areaCu(I,j)
      enddo; enddo
      do J=js-1,Jeq ; do i=is,ie
        CS%L2v(i,J) = CS%Visbeck_L_scale**2 * G%areaCv(i,J)
      enddo; enddo
    else
      CS%L2u(:,:) = US%m_to_L**2*CS%Visbeck_L_scale**2
      CS%L2v(:,:) = US%m_to_L**2*CS%Visbeck_L_scale**2
    endif

    CS%id_L2u = register_diag_field('ocean_model', 'L2u', diag%axesCu1, Time, &
       'Length scale squared for mixing coefficient, at u-points', &
       'm2', conversion=US%L_to_m**2)
    CS%id_L2v = register_diag_field('ocean_model', 'L2v', diag%axesCv1, Time, &
       'Length scale squared for mixing coefficient, at v-points', &
       'm2', conversion=US%L_to_m**2)
  endif

  if (CS%calculate_Eady_growth_rate .and. CS%use_stored_slopes) then
    CS%id_N2_u = register_diag_field('ocean_model', 'N2_u', diag%axesCui, Time, &
         'Square of Brunt-Vaisala frequency, N^2, at u-points, as used in Visbeck et al.', &
         's-2', conversion=US%s_to_T**2)
    CS%id_N2_v = register_diag_field('ocean_model', 'N2_v', diag%axesCvi, Time, &
         'Square of Brunt-Vaisala frequency, N^2, at v-points, as used in Visbeck et al.', &
         's-2', conversion=US%s_to_T**2)
  endif
  if (CS%use_stored_slopes) then
    CS%id_S2_u = register_diag_field('ocean_model', 'S2_u', diag%axesCu1, Time, &
         'Depth average square of slope magnitude, S^2, at u-points, as used in Visbeck et al.', 'nondim')
    CS%id_S2_v = register_diag_field('ocean_model', 'S2_v', diag%axesCv1, Time, &
         'Depth average square of slope magnitude, S^2, at v-points, as used in Visbeck et al.', 'nondim')
  endif

  oneOrTwo = 1.0
  if (CS%Resoln_scaled_Kh .or. CS%Resoln_scaled_KhTh .or. CS%Resoln_scaled_KhTr) then
    CS%calculate_Rd_dx = .true.
    CS%calculate_res_fns = .true.
    allocate(CS%Res_fn_h(isd:ied,jsd:jed))       ; CS%Res_fn_h(:,:) = 0.0
    allocate(CS%Res_fn_q(IsdB:IedB,JsdB:JedB))   ; CS%Res_fn_q(:,:) = 0.0
    allocate(CS%Res_fn_u(IsdB:IedB,jsd:jed))     ; CS%Res_fn_u(:,:) = 0.0
    allocate(CS%Res_fn_v(isd:ied,JsdB:JedB))     ; CS%Res_fn_v(:,:) = 0.0
    allocate(CS%beta_dx2_q(IsdB:IedB,JsdB:JedB)) ; CS%beta_dx2_q(:,:) = 0.0
    allocate(CS%beta_dx2_u(IsdB:IedB,jsd:jed))   ; CS%beta_dx2_u(:,:) = 0.0
    allocate(CS%beta_dx2_v(isd:ied,JsdB:JedB))   ; CS%beta_dx2_v(:,:) = 0.0
    allocate(CS%f2_dx2_q(IsdB:IedB,JsdB:JedB))   ; CS%f2_dx2_q(:,:) = 0.0
    allocate(CS%f2_dx2_u(IsdB:IedB,jsd:jed))     ; CS%f2_dx2_u(:,:) = 0.0
    allocate(CS%f2_dx2_v(isd:ied,JsdB:JedB))     ; CS%f2_dx2_v(:,:) = 0.0

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
                 units="nondim", default=2)
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
                 units="nondim", default=CS%Res_fn_power_khth)
    call get_param(param_file, mdl, "INTERPOLATE_RES_FN", CS%interpolate_Res_fn, &
                 "If true, interpolate the resolution function to the "//&
                 "velocity points from the thickness points; otherwise "//&
                 "interpolate the wave speed and calculate the resolution "//&
                 "function independently at each point.", default=.true.)
    if (CS%interpolate_Res_fn) then
      if (CS%Res_coef_visc /= CS%Res_coef_khth) call MOM_error(FATAL, &
           "MOM_lateral_mixing_coeffs.F90, VarMix_init:"//&
           "When INTERPOLATE_RES_FN=True, VISC_RES_FN_POWER must equal KH_RES_SCALE_COEF.")
      if (CS%Res_fn_power_visc /= CS%Res_fn_power_khth) call MOM_error(FATAL, &
           "MOM_lateral_mixing_coeffs.F90, VarMix_init:"//&
           "When INTERPOLATE_RES_FN=True, VISC_RES_FN_POWER must equal KH_RES_FN_POWER.")
    endif
    !### Change the default of GILL_EQUATORIAL_LD to True.
    call get_param(param_file, mdl, "GILL_EQUATORIAL_LD", Gill_equatorial_Ld, &
                 "If true, uses Gill's definition of the baroclinic "//&
                 "equatorial deformation radius, otherwise, if false, use "//&
                 "Pedlosky's definition. These definitions differ by a factor "//&
                 "of 2 in front of the beta term in the denominator. Gill's "//&
                 "is the more appropriate definition.", default=.false.)
    if (Gill_equatorial_Ld) then
      oneOrTwo = 2.0
    endif

    do J=js-1,Jeq ; do I=is-1,Ieq
      CS%f2_dx2_q(I,J) = (G%dxBu(I,J)**2 + G%dyBu(I,J)**2) * &
                         max(G%CoriolisBu(I,J)**2, absurdly_small_freq**2)
      CS%beta_dx2_q(I,J) = oneOrTwo * ((G%dxBu(I,J))**2 + (G%dyBu(I,J))**2) * (sqrt(0.5 * &
          ( (((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
             ((G%CoriolisBu(I+1,J)-G%CoriolisBu(I,J)) * G%IdxCv(i+1,J))**2) + &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
             ((G%CoriolisBu(I,J+1)-G%CoriolisBu(I,J)) * G%IdyCu(I,j+1))**2) ) ))
    enddo ; enddo

    do j=js,je ; do I=is-1,Ieq
      CS%f2_dx2_u(I,j) = (G%dxCu(I,j)**2 + G%dyCu(I,j)**2) * &
          max(0.5* (G%CoriolisBu(I,J)**2+G%CoriolisBu(I,J-1)**2), absurdly_small_freq**2)
      CS%beta_dx2_u(I,j) = oneOrTwo * ((G%dxCu(I,j))**2 + (G%dyCu(I,j))**2) * (sqrt( &
          0.25*( (((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2 + &
                  ((G%CoriolisBu(I+1,J)-G%CoriolisBu(I,J)) * G%IdxCv(i+1,J))**2) + &
                 (((G%CoriolisBu(I+1,J-1)-G%CoriolisBu(I,J-1)) * G%IdxCv(i+1,J-1))**2 + &
                  ((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2) ) + &
                  ((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 ))
    enddo ; enddo

    do J=js-1,Jeq ; do i=is,ie
      CS%f2_dx2_v(i,J) = ((G%dxCv(i,J))**2 + (G%dyCv(i,J))**2) * &
          max(0.5*(G%CoriolisBu(I,J)**2+G%CoriolisBu(I-1,J)**2), absurdly_small_freq**2)
      CS%beta_dx2_v(i,J) = oneOrTwo * ((G%dxCv(i,J))**2 + (G%dyCv(i,J))**2) * (sqrt( &
          ((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
          0.25*( (((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
                  ((G%CoriolisBu(I-1,J+1)-G%CoriolisBu(I-1,J)) * G%IdyCu(I-1,j+1))**2) + &
                 (((G%CoriolisBu(I,J+1)-G%CoriolisBu(I,J)) * G%IdyCu(I,j+1))**2 + &
                  ((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2) ) ))
    enddo ; enddo

  endif

  if (CS%Depth_scaled_KhTh) then
    CS%calculate_depth_fns = .true.
    allocate(CS%Depth_fn_u(IsdB:IedB,jsd:jed))     ; CS%Depth_fn_u(:,:) = 0.0
    allocate(CS%Depth_fn_v(isd:ied,JsdB:JedB))     ; CS%Depth_fn_v(:,:) = 0.0
    call get_param(param_file, mdl, "DEPTH_SCALED_KHTH_H0", CS%depth_scaled_khth_h0, &
    "The depth above which KHTH is scaled away.",&
    units="m", default=1000.)
    call get_param(param_file, mdl, "DEPTH_SCALED_KHTH_EXP", CS%depth_scaled_khth_exp, &
    "The exponent used in the depth dependent scaling function for KHTH.",&
    units="nondim", default=3.0)
  endif

  ! Resolution %Rd_dx_h
  CS%id_Rd_dx = register_diag_field('ocean_model', 'Rd_dx', diag%axesT1, Time, &
       'Ratio between deformation radius and grid spacing', 'm m-1')
  CS%calculate_Rd_dx = CS%calculate_Rd_dx .or. (CS%id_Rd_dx>0)

  if (CS%calculate_Rd_dx) then
    CS%calculate_cg1 = .true. ! We will need %cg1
    allocate(CS%Rd_dx_h(isd:ied,jsd:jed))   ; CS%Rd_dx_h(:,:) = 0.0
    allocate(CS%beta_dx2_h(isd:ied,jsd:jed)); CS%beta_dx2_h(:,:) = 0.0
    allocate(CS%f2_dx2_h(isd:ied,jsd:jed))  ; CS%f2_dx2_h(:,:) = 0.0
    do j=js-1,je+1 ; do i=is-1,ie+1
      CS%f2_dx2_h(i,j) = (G%dxT(i,j)**2 + G%dyT(i,j)**2) * &
          max(0.25 * ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                      (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2)), &
              absurdly_small_freq**2)
      CS%beta_dx2_h(i,j) = oneOrTwo * ((G%dxT(i,j))**2 + (G%dyT(i,j))**2) * (sqrt(0.5 * &
          ( (((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
             ((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2) + &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
             ((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2) ) ))
    enddo ; enddo
  endif

  if (CS%calculate_cg1) then
    in_use = .true.
    allocate(CS%cg1(isd:ied,jsd:jed)) ; CS%cg1(:,:) = 0.0
    call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.true.)
    call get_param(param_file, mdl, "REMAPPING_2018_ANSWERS", remap_answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
    call get_param(param_file, mdl, "INTERNAL_WAVE_SPEED_TOL", wave_speed_tol, &
                 "The fractional tolerance for finding the wave speeds.", &
                 units="nondim", default=0.001)
    !### Set defaults so that wave_speed_min*wave_speed_tol >= 1e-9 m s-1
    call get_param(param_file, mdl, "INTERNAL_WAVE_SPEED_MIN", wave_speed_min, &
                 "A floor in the first mode speed below which 0 used instead.", &
                 units="m s-1", default=0.0, scale=US%m_s_to_L_T)
    call get_param(param_file, mdl, "INTERNAL_WAVE_SPEED_BETTER_EST", better_speed_est, &
                 "If true, use a more robust estimate of the first mode wave speed as the "//&
                 "starting point for iterations.", default=.false.) !### Change the default.
    call wave_speed_init(CS%wave_speed_CSp, use_ebt_mode=CS%Resoln_use_ebt, &
                         mono_N2_depth=N2_filter_depth, remap_answers_2018=remap_answers_2018, &
                         better_speed_est=better_speed_est, min_speed=wave_speed_min, &
                         wave_speed_tol=wave_speed_tol)
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
    ALLOC_(CS%KH_u_QG(IsdB:IedB,jsd:jed,G%ke)) ; CS%KH_u_QG(:,:,:) = 0.0
    ALLOC_(CS%KH_v_QG(isd:ied,JsdB:JedB,G%ke)) ; CS%KH_v_QG(:,:,:) = 0.0
    ! register diagnostics

    CS%id_KH_u_QG = register_diag_field('ocean_model', 'KH_u_QG', diag%axesCuL, Time, &
       'Horizontal viscosity from Leith QG, at u-points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
    CS%id_KH_v_QG = register_diag_field('ocean_model', 'KH_v_QG', diag%axesCvL, Time, &
       'Horizontal viscosity from Leith QG, at v-points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)

    do j=Jsq,Jeq+1 ; do I=is-1,Ieq
      ! Static factors in the Leith schemes
      grid_sp_u2 = G%dyCu(I,j)*G%dxCu(I,j)
      grid_sp_u3 = sqrt(grid_sp_u2)
      CS%Laplac3_const_u(I,j) = Leith_Lap_const * grid_sp_u3
    enddo ; enddo
    do j=js-1,Jeq ; do I=Isq,Ieq+1
      ! Static factors in the Leith schemes
      !### The second factor here is wrong.  It should be G%dxCv(i,J).
      grid_sp_v2 = G%dyCv(i,J)*G%dxCu(i,J)
      grid_sp_v3 = grid_sp_v2*sqrt(grid_sp_v2)
      CS%Laplac3_const_v(i,J) = Leith_Lap_const * grid_sp_v3
    enddo ; enddo

    if (.not. CS%use_stored_slopes) call MOM_error(FATAL, &
           "MOM_lateral_mixing_coeffs.F90, VarMix_init:"//&
           "USE_STORED_SLOPES must be True when using QG Leith.")
  endif

  ! If nothing is being stored in this class then deallocate
  if (in_use) then
    CS%use_variable_mixing = .true.
  else
    deallocate(CS)
    return
  endif

end subroutine VarMix_init

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
!! The resolution function used in scaling diffusivities (Hallberg, 2010) is
!!
!! \f[
!! r(\Delta,L_d) = \frac{1}{1+(\alpha R)^p}
!! \f]
!!
!! The resolution function can be applied independently to thickness diffusion (module mom_thickness_diffuse),
!! tracer diffusion (mom_tracer_hordiff) lateral viscosity (mom_hor_visc).
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
!! scheme.  The factors are combined in mom_thickness_diffuse::thickness_diffuse() but calculated in this module.
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
!! velocity mode.  The structure function is stored in the control structure for thie module (varmix_cs) but is
!! calculated using subroutines in mom_wave_speed.
!!
!! | Symbol                | Module parameter |
!! | ------                | --------------- |
!! | -                     | <code>KHTH_USE_EBT_STRUCT</code> |

end module MOM_lateral_mixing_coeffs
