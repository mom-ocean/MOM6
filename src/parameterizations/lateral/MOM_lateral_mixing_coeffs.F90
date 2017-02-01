!> Variable mixing coefficients
module MOM_lateral_mixing_coeffs

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, post_data
use MOM_diag_mediator, only : diag_ctrl, time_type, query_averaging_enabled
use MOM_domains,       only : create_group_pass, do_group_pass
use MOM_domains,       only : group_pass_type, pass_var
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_interface_heights, only : find_eta
use MOM_isopycnal_slopes, only : calc_isoneutral_slopes
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_wave_speed, only : wave_speed, wave_speed_CS, wave_speed_init

implicit none ; private

#include <MOM_memory.h>

!> Variable mixing coefficients
type, public :: VarMix_CS ;
  logical :: use_variable_mixing  !< If true, use the variable mixing.
  logical :: Resoln_scaled_Kh     !< If true, scale away the Laplacian viscosity
                                  !! when the deformation radius is well resolved.
  logical :: Resoln_scaled_KhTh   !< If true, scale away the thickness diffusivity
                                  !! when the deformation radius is well resolved.
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
  real, dimension(:,:), pointer :: &
    SN_u => NULL(), &   !< S*N at u-points (s^-1)
    SN_v => NULL(), &  !< S*N at v-points (s^-1)
    L2u => NULL(), &   !< Length scale^2 at u-points (m^2)
    L2v => NULL(), &   !< Length scale^2 at v-points (m^2)
    cg1 => NULL(), &   !< The first baroclinic gravity wave speed in m s-1.
    Res_fn_h => NULL(), & !< Non-dimensional function of the ratio the first baroclinic
                          !! deformation radius to the grid spacing at h points.
    Res_fn_q => NULL(), & !< Non-dimensional function of the ratio the first baroclinic
                          !! deformation radius to the grid spacing at q points.
    Res_fn_u => NULL(), & !< Non-dimensional function of the ratio the first baroclinic
                          !! deformation radius to the grid spacing at u points.
    Res_fn_v => NULL(), & !< Non-dimensional function of the ratio the first baroclinic
                          !! deformation radius to the grid spacing at v points.
    beta_dx2_h => NULL(), & !< The magnitude of the gradient of the Coriolis parameter
                            !! times the grid spacing squared at h points.
    beta_dx2_q => NULL(), & !< The magnitude of the gradient of the Coriolis parameter
                            !! times the grid spacing squared at q points.
    beta_dx2_u => NULL(), & !< The magnitude of the gradient of the Coriolis parameter
                            !! times the grid spacing squared at u points.
    beta_dx2_v => NULL(), & !< The magnitude of the gradient of the Coriolis parameter
                            !! times the grid spacing squared at v points.
    f2_dx2_h => NULL(), & !< The Coriolis parameter squared times the grid
                          !! spacing squared at h, in m2 s-2.
    f2_dx2_q => NULL(), & !< The Coriolis parameter squared times the grid
                          !! spacing squared at q, in m2 s-2.
    f2_dx2_u => NULL(), & !< The Coriolis parameter squared times the grid
                          !! spacing squared at u, in m2 s-2.
    f2_dx2_v => NULL(), & !< The Coriolis parameter squared times the grid
                          !! spacing squared at v, in m2 s-2.
    Rd_dx_h => NULL()     !< Deformation radius over grid spacing (non-dim.)

  real, dimension(:,:,:), pointer :: &
    slope_x => NULL(), &  !< Zonal isopycnal slope (non-dimensional)
    slope_y => NULL(), &  !< Meridional isopycnal slope (non-dimensional)
    ebt_struct => NULL()  !< Vertical structure function to scale diffusivities with (non-dim)

  ! Parameters
  integer :: VarMix_Ktop  !< Top layer to start downward integrals
  real :: Visbeck_L_scale !< Fixed length scale in Visbeck formula
  real :: Res_coef_khth   !< A non-dimensional number that determines the function
                          !! of resolution, used for thickness and tracer mixing, as:
                          !!  F = 1 / (1 + (Res_coef_khth*Ld/dx)^Res_fn_power)
  real :: Res_coef_visc   !< A non-dimensional number that determines the function
                          !! of resolution, used for lateral viscosity, as:
                          !!  F = 1 / (1 + (Res_coef_visc*Ld/dx)^Res_fn_power)
  real :: kappa_smooth    !< A diffusivity for smoothing T/S in vanished layers (m2/s)
  integer :: Res_fn_power_khth !< The power of dx/Ld in the KhTh resolution function.  Any
                               !! positive integer power may be used, but even powers
                               !! and especially 2 are coded to be more efficient.
  integer :: Res_fn_power_visc !< The power of dx/Ld in the Kh resolution function.  Any
                               !! positive integer power may be used, but even powers
                               !! and especially 2 are coded to be more efficient.
  real :: Visbeck_S_max   !< Upper bound on slope used in Eady growth rate (nondim).

  ! Diagnostics
  !>@{
  !! Diagnostic identifier
  integer :: id_SN_u=-1, id_SN_v=-1, id_L2u=-1, id_L2v=-1, id_Res_fn = -1
  integer :: id_N2_u=-1, id_N2_v=-1, id_S2_u=-1, id_S2_v=-1
  integer :: id_Rd_dx=-1
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.
  !>@}

  type(wave_speed_CS), pointer :: wave_speed_CSp => NULL() !< Wave speed control structure
  type(group_pass_type) :: pass_cg1 !< For group halo pass

end type VarMix_CS

public VarMix_init, calc_slope_functions, calc_resoln_function

contains

!> Calculates and stores the non-dimensional resolution functions
subroutine calc_resoln_function(h, tv, G, GV, CS)
  type(ocean_grid_type),                    intent(inout) :: G  !< Ocean grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h  !< Layer thickness (m or kg/m2)
  type(thermo_var_ptrs),                    intent(in)    :: tv !< Thermodynamic variables
  type(verticalGrid_type),                  intent(in)    :: GV !< Vertical grid structure
  type(VarMix_CS),                          pointer       :: CS !< Variable mixing coefficients
  ! Local variables
  real :: cg1_q  ! The gravity wave speed interpolated to q points, in m s-1.
  real :: cg1_u  ! The gravity wave speed interpolated to u points, in m s-1.
  real :: cg1_v  ! The gravity wave speed interpolated to v points, in m s-1.
  real :: dx_term
  integer :: power_2
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not. ASSOCIATED(CS)) call MOM_error(FATAL, "calc_resoln_function:"// &
         "Module must be initialized before it is used.")
  if (.not. (CS%Resoln_scaled_Kh .or. CS%Resoln_scaled_KhTh .or. &
             CS%Resoln_scaled_KhTr)) return
  if (.not. ASSOCIATED(CS%cg1)) call MOM_error(FATAL, &
    "calc_resoln_function: %cg1 is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%Res_fn_h)) call MOM_error(FATAL, &
    "calc_resoln_function: %Res_fn_h is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%Res_fn_q)) call MOM_error(FATAL, &
    "calc_resoln_function: %Res_fn_q is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%Res_fn_u)) call MOM_error(FATAL, &
    "calc_resoln_function: %Res_fn_u is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%Res_fn_v)) call MOM_error(FATAL, &
    "calc_resoln_function: %Res_fn_v is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%f2_dx2_h)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_h is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%f2_dx2_q)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_q is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%f2_dx2_u)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_u is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%f2_dx2_v)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_v is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%beta_dx2_h)) call MOM_error(FATAL, &
    "calc_resoln_function: %beta_dx2_h is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%beta_dx2_q)) call MOM_error(FATAL, &
    "calc_resoln_function: %beta_dx2_q is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%beta_dx2_u)) call MOM_error(FATAL, &
    "calc_resoln_function: %beta_dx2_u is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%beta_dx2_v)) call MOM_error(FATAL, &
    "calc_resoln_function: %beta_dx2_v is not associated with Resoln_scaled_Kh.")

  if (CS%khth_use_ebt_struct) then
    if (CS%Resoln_use_ebt) then
      ! Both resolution fn and vertical structure are using EBT
      call wave_speed(h, tv, G, GV, CS%cg1, CS%wave_speed_CSp, modal_structure=CS%ebt_struct)
    else
      ! Use EBT to get vertical structure first and then re-calculate cg1 using first baroclinic mode
      call wave_speed(h, tv, G, GV, CS%cg1, CS%wave_speed_CSp, modal_structure=CS%ebt_struct, use_ebt_mode=.true.)
      call wave_speed(h, tv, G, GV, CS%cg1, CS%wave_speed_CSp)
    endif
    call pass_var(CS%ebt_struct, G%Domain)
  else
    call wave_speed(h, tv, G, GV, CS%cg1, CS%wave_speed_CSp)
  endif

  call create_group_pass(CS%pass_cg1, CS%cg1, G%Domain)
  call do_group_pass(CS%pass_cg1, G%Domain)

  !   Do this calculation on the extent used in MOM_hor_visc.F90, and
  ! MOM_tracer.F90 so that no halo update is needed.

!$OMP parallel default(none) shared(is,ie,js,je,Ieq,Jeq,CS) &
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
      cg1_q = 0.25 * ((CS%cg1(i,j) + CS%cg1(i+1,j+1)) + &
                      (CS%cg1(i+1,j) + CS%cg1(i,j+1)))
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
      cg1_q = 0.25 * ((CS%cg1(i,j) + CS%cg1(i+1,j+1)) + &
                      (CS%cg1(i+1,j) + CS%cg1(i,j+1)))
      dx_term = CS%f2_dx2_q(I,J) +  cg1_q * CS%beta_dx2_q(I,J)
      CS%Res_fn_q(I,J) = dx_term / (dx_term + (CS%Res_coef_visc * cg1_q)**2)
    enddo ; enddo
  elseif (mod(CS%Res_fn_power_visc, 2) == 0) then
    power_2 = CS%Res_fn_power_visc / 2
!$OMP do
    do j=js-1,je+1 ; do i=is-1,ie+1
      dx_term = (CS%f2_dx2_h(i,j) + CS%cg1(i,j)*CS%beta_dx2_h(i,j))**power_2
      CS%Res_fn_h(i,j) = dx_term / &
          (dx_term + (CS%Res_coef_visc * CS%cg1(i,j))**CS%Res_fn_power_visc)
    enddo ; enddo
!$OMP do
    do J=js-1,Jeq ; do I=is-1,Ieq
      cg1_q = 0.25 * ((CS%cg1(i,j) + CS%cg1(i+1,j+1)) + &
                      (CS%cg1(i+1,j) + CS%cg1(i,j+1)))
      dx_term = (CS%f2_dx2_q(I,J) +  cg1_q * CS%beta_dx2_q(I,J))**power_2
      CS%Res_fn_q(I,J) = dx_term / &
          (dx_term + (CS%Res_coef_visc * cg1_q)**CS%Res_fn_power_visc)
    enddo ; enddo
  else
!$OMP do
    do j=js-1,je+1 ; do i=is-1,ie+1
      dx_term = (sqrt(CS%f2_dx2_h(i,j) + &
                      CS%cg1(i,j)*CS%beta_dx2_h(i,j)))**CS%Res_fn_power_visc
      CS%Res_fn_h(i,j) = dx_term / &
         (dx_term + (CS%Res_coef_visc * CS%cg1(i,j))**CS%Res_fn_power_visc)
    enddo ; enddo
!$OMP do
    do J=js-1,Jeq ; do I=is-1,Ieq
      cg1_q = 0.25 * ((CS%cg1(i,j) + CS%cg1(i+1,j+1)) + &
                      (CS%cg1(i+1,j) + CS%cg1(i,j+1)))
      dx_term = (sqrt(CS%f2_dx2_q(I,J) + &
                      cg1_q * CS%beta_dx2_q(I,J)))**CS%Res_fn_power_visc
      CS%Res_fn_q(I,J) = dx_term / &
          (dx_term + (CS%Res_coef_visc * cg1_q)**CS%Res_fn_power_visc)
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
        dx_term = (CS%f2_dx2_u(I,j) + cg1_u * CS%beta_dx2_u(I,j))**power_2
        CS%Res_fn_u(I,j) = dx_term / &
            (dx_term + (CS%Res_coef_khth * cg1_u)**CS%Res_fn_power_khth)
      enddo ; enddo
!$OMP do
      do J=js-1,Jeq ; do i=is,ie
        cg1_v = 0.5 * (CS%cg1(i,j) + CS%cg1(i,j+1))
        dx_term = (CS%f2_dx2_v(i,J) + cg1_v * CS%beta_dx2_v(i,J))**power_2
        CS%Res_fn_v(i,J) = dx_term / &
            (dx_term + (CS%Res_coef_khth * cg1_v)**CS%Res_fn_power_khth)
      enddo ; enddo
    else
!$OMP do
      do j=js,je ; do I=is-1,Ieq
        cg1_u = 0.5 * (CS%cg1(i,j) + CS%cg1(i+1,j))
        dx_term = (sqrt(CS%f2_dx2_u(I,j) + &
                        cg1_u * CS%beta_dx2_u(I,j)))**CS%Res_fn_power_khth
        CS%Res_fn_u(I,j) = dx_term / &
            (dx_term + (CS%Res_coef_khth * cg1_u)**CS%Res_fn_power_khth)
      enddo ; enddo
!$OMP do
      do J=js-1,Jeq ; do i=is,ie
        cg1_v = 0.5 * (CS%cg1(i,j) + CS%cg1(i,j+1))
        dx_term = (sqrt(CS%f2_dx2_v(i,J) + &
                        cg1_v * CS%beta_dx2_v(i,J)))**CS%Res_fn_power_khth
        CS%Res_fn_v(i,J) = dx_term / &
            (dx_term + (CS%Res_coef_khth * cg1_v)**CS%Res_fn_power_khth)
      enddo ; enddo
    endif
  endif

  ! Calculate and store the ratio between deformation radius and grid-spacing
  ! at h-points (non-dimensional).
!$OMP do
  do j=js-1,je+1 ; do i=is-1,ie+1
    CS%Rd_dx_h(i,j) = CS%cg1(i,j) / &
          (sqrt(CS%f2_dx2_h(i,j) + CS%cg1(i,j)*CS%beta_dx2_h(i,j)))
  enddo ; enddo
!$OMP end parallel

  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_Res_fn > 0) call post_data(CS%id_Res_fn, CS%Res_fn_h, CS%diag)
    if (CS%id_Rd_dx > 0) call post_data(CS%id_Rd_dx, CS%Rd_dx_h, CS%diag)
  endif

end subroutine calc_resoln_function

!> Calculates and stores functions of isopycnal slopes, e.g. Sx, Sy, S*N, mostly used in the Visbeck et al.
!! style scaling of diffusivity
subroutine calc_slope_functions(h, tv, dt, G, GV, CS)
  type(ocean_grid_type),                    intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                  intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: h  !< Layer thickness (m or kg/m2)
  type(thermo_var_ptrs),                    intent(in)    :: tv !< Thermodynamic variables
  real,                                     intent(in)    :: dt !< Time increment (s)
  type(VarMix_CS),                          pointer       :: CS !< Variable mixing coefficients
  ! Local variables
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)+1) :: &
    e             ! The interface heights relative to mean sea level, in m.
  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)+1) :: N2_u ! Square of Brunt-Vaisala freq at u-points
  real, dimension(SZI_(G), SZJB_(G), SZK_(G)+1) :: N2_v ! Square of Brunt-Vaisala freq at u-points

  if (.not. ASSOCIATED(CS)) call MOM_error(FATAL, "MOM_lateral_mixing_coeffs.F90, calc_slope_functions:"//&
         "Module must be initialized before it is used.")

  call find_eta(h, tv, GV%g_Earth, G, GV, e, halo_size=2)
  if (CS%use_variable_mixing) then
    if (CS%use_stored_slopes) then
      call calc_isoneutral_slopes(G, GV, h, e, tv, dt*CS%kappa_smooth, &
                                  CS%slope_x, CS%slope_y, N2_u, N2_v, 1)
      call calc_Visbeck_coeffs(h, e, CS%slope_x, CS%slope_y, N2_u, N2_v, G, GV, CS)
!     call calc_slope_functions_using_just_e(h, G, CS, e, .false.)
    else
      !call calc_isoneutral_slopes(G, GV, h, e, tv, dt*CS%kappa_smooth, CS%slope_x, CS%slope_y)
      call calc_slope_functions_using_just_e(h, G, GV, CS, e, .true.)
    endif
  endif

  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_SN_u > 0) call post_data(CS%id_SN_u, CS%SN_u, CS%diag)
    if (CS%id_SN_v > 0) call post_data(CS%id_SN_v, CS%SN_v, CS%diag)
    if (CS%id_L2u > 0) call post_data(CS%id_L2u, CS%L2u, CS%diag)
    if (CS%id_L2v > 0) call post_data(CS%id_L2v, CS%L2v, CS%diag)
    if (CS%use_stored_slopes) then
      if (CS%id_N2_u > 0) call post_data(CS%id_N2_u, N2_u, CS%diag)
      if (CS%id_N2_v > 0) call post_data(CS%id_N2_v, N2_v, CS%diag)
    endif
  endif

end subroutine calc_slope_functions

!> Calculates factors used when setting diffusivity coefficients similar to Visbeck et al.
subroutine calc_Visbeck_coeffs(h, e, slope_x, slope_y, N2_u, N2_v, G, GV, CS)
  type(ocean_grid_type),                       intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type),                     intent(in)    :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),    intent(in)    :: h  !< Layer thickness (m or kg/m2)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1),  intent(in)    :: e  !< Interface position (m)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1), intent(in)    :: slope_x !< Zonal isoneutral slope
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1), intent(in)    :: N2_u    !< Brunt-Vaisala frequency at u-points (1/s2)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1), intent(in)    :: slope_y !< Meridional isoneutral slope
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1), intent(in)    :: N2_v    !< Brunt-Vaisala frequency at v-points (1/s2)
  type(VarMix_CS),                             intent(inout) :: CS !< Variable mixing coefficients
  ! Local variables
  real :: E_x(SZIB_(G), SZJ_(G))  ! X-slope of interface at u points (for diagnostics)
  real :: E_y(SZI_(G), SZJB_(G))  ! Y-slope of interface at u points (for diagnostics)
  real :: Khth_Loc      ! Locally calculated thickness mixing coefficient (m2/s)
  real :: S2            ! Interface slope squared (non-dim)
  real :: N2            ! Brunt-Vaisala frequency (1/s)
  real :: Hup, Hdn      ! Thickness from above, below (m or kg m-2)
  real :: H_geom        ! The geometric mean of Hup*Hdn, in m or kg m-2.
  integer :: is, ie, js, je, nz
  integer :: i, j, k, kb_max
  real :: S2max, wNE, wSE, wSW, wNW
  real :: SN_u_local(SZIB_(G), SZJ_(G),SZK_(G))
  real :: SN_v_local(SZI_(G), SZJB_(G),SZK_(G))
  real :: H_u(SZIB_(G)), H_v(SZI_(G))
  real :: S2_u(SZIB_(G), SZJ_(G))
  real :: S2_v(SZI_(G), SZJB_(G))

  if (LOC(CS)==0) call MOM_error(FATAL, "calc_slope_function:"// &
         "Module must be initialized before it is used.")
  if (.not. CS%use_variable_mixing) return
  if (.not. ASSOCIATED(CS%SN_u)) call MOM_error(FATAL, "calc_slope_function:"// &
         "%SN_u is not associated with use_variable_mixing.")
  if (.not. ASSOCIATED(CS%SN_v)) call MOM_error(FATAL, "calc_slope_function:"// &
         "%SN_v is not associated with use_variable_mixing.")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  S2max = CS%Visbeck_S_max**2

!$OMP parallel default(none) shared(is,ie,js,je,CS,nz,e,G,GV,h, &
!$OMP                               S2_u,S2_v,slope_x,slope_y,   &
!$OMP                               SN_u_local,SN_v_local,N2_u,N2_v, S2max)   &
!$OMP                       private(E_x,E_y,S2,H_u,H_v,Hdn,Hup,H_geom,N2, &
!$OMP                       wNE, wSE, wSW, wNW)
!$OMP do
  do j=js-1,je+1 ; do i=is-1,ie+1
    CS%SN_u(i,j) = 0.0
    CS%SN_v(i,j) = 0.0
  enddo ; enddo


  ! To set the length scale based on the deformation radius, use wave_speed to
  ! calculate the first-mode gravity wave speed and then blend the equatorial
  ! and midlatitude deformation radii, using calc_resoln_function as a template.

  ! Set the length scale at u-points.
!$OMP do
  do j=js,je ; do I=is-1,ie
    CS%L2u(I,j) = CS%Visbeck_L_scale**2
  enddo ; enddo
  ! Set length scale at v-points
!$OMP do
  do J=js-1,je ; do i=is,ie
    CS%L2v(i,J) = CS%Visbeck_L_scale**2
  enddo ; enddo

!$OMP do
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
      wSE = h(i+1,j,k)*h(i+1,j-1,k) * h(i+1,j,k)*h(i+1,j-1,k-1)
      wNW = h(i  ,j,k)*h(i  ,j+1,k) * h(i  ,j,k)*h(i  ,j+1,k-1)
      wNE = h(i+1,j,k)*h(i+1,j+1,k) * h(i+1,j,k)*h(i+1,j+1,k-1)
      wSW = h(i  ,j,k)*h(i  ,j-1,k) * h(i  ,j,k)*h(i  ,j-1,k-1)
      S2 =  slope_x(I,j,K)**2  + ( &
           (wNW*slope_y(i,J,K)**2+wSE*slope_y(i+1,J-1,K)**2)     &
          +(wNE*slope_y(i+1,J,K)**2+wSW*slope_y(i,J-1,K)**2) ) / &
           ( ((wSE+wNW) + (wNE+wSW)) + GV%H_subroundoff**2 ) !### This should be **4 for consistent units.
      if (S2max>0.) S2 = S2 * S2max / (S2 + S2max) ! Limit S2
      N2 = max(0., N2_u(I,j,k))
      CS%SN_u(I,j) = CS%SN_u(I,j) + sqrt( S2*N2 )*H_geom
      S2_u(I,j) = S2_u(I,j) + S2*H_geom
      H_u(I) = H_u(I) + H_geom
    enddo ; enddo
    do I=is-1,ie
      if (H_u(I)>0.) then
        CS%SN_u(I,j) = CS%SN_u(I,j) / H_u(I)
        S2_u(I,j) = S2_u(I,j) / H_u(I)
      else
        CS%SN_u(I,j) = 0.
      endif
    enddo
  enddo

!$OMP do
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
      wSE = h(i,j  ,k)*h(i+1,j  ,k) * h(i,j  ,k)*h(i+1,j  ,k-1)
      wNW = h(i,j+1,k)*h(i-1,j+1,k) * h(i,j+1,k)*h(i-1,j+1,k-1)
      wNE = h(i,j+1,k)*h(i+1,j+1,k) * h(i,j+1,k)*h(i+1,j+1,k-1)
      wSW = h(i,j  ,k)*h(i-1,j  ,k) * h(i,j  ,k)*h(i-1,j  ,k-1)
      S2 =  slope_y(i,J,K)**2  + ( &
           (wSE*slope_x(I,j,K)**2+wNW*slope_x(I-1,j+1,K)**2)     &
          +(wNE*slope_x(I,j+1,K)**2+wSW*slope_x(I-1,j,K)**2) ) / &
           ( ((wSE+wNW) + (wNE+wSW)) + GV%H_subroundoff**2 ) !### This should be **4 for consistent units.
      if (S2max>0.) S2 = S2 * S2max / (S2 + S2max) ! Limit S2
      N2 = max(0., N2_v(i,J,K))
      CS%SN_v(i,J) = CS%SN_v(i,J) + sqrt( S2*N2 )*H_geom
      S2_v(i,J) = S2_v(i,J) + S2*H_geom
      H_v(i) = H_v(i) + H_geom
    enddo ; enddo
    do i=is,ie
      if (H_v(i)>0.) then
        CS%SN_v(i,J) = CS%SN_v(i,J) / H_v(i)
        S2_v(i,J) = S2_v(i,J) / H_v(i)
      else
        CS%SN_v(i,J) = 0.
      endif
    enddo
  enddo

!$OMP end parallel

! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_S2_u > 0) call post_data(CS%id_S2_u, S2_u, CS%diag)
    if (CS%id_S2_v > 0) call post_data(CS%id_S2_v, S2_v, CS%diag)
  endif

end subroutine calc_Visbeck_coeffs

!> The original calc_slope_function() that calculated slopes using
!! interface positions only, not accounting for density variations.
subroutine calc_slope_functions_using_just_e(h, G, GV, CS, e, calculate_slopes)
  type(ocean_grid_type),                      intent(inout) :: G  !< Ocean grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(inout) :: h  !< Layer thickness (m or kg/m2)
  type(verticalGrid_type),                    intent(in)    :: GV !< Vertical grid structure
  type(VarMix_CS),                            pointer       :: CS !< Variable mixing coefficients
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(in)    :: e  !< Interface position (m)
  logical,                                    intent(in)    :: calculate_slopes !< If true, calculate slopes internally
                                                                  !! otherwise use slopes stored in CS
  ! Local variables
  real :: E_x(SZIB_(G), SZJ_(G))  ! X-slope of interface at u points (for diagnostics)
  real :: E_y(SZI_(G), SZJB_(G))  ! Y-slope of interface at u points (for diagnostics)
  real :: Khth_Loc      ! Locally calculated thickness mixing coefficient (m2/s)
  real :: H_cutoff      ! Local estimate of a minimum thickness for masking (m)
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected, in H.
  real :: S2            ! Interface slope squared (non-dim)
  real :: N2            ! Brunt-Vaisala frequency (1/s)
  real :: Hup, Hdn      ! Thickness from above, below (m or kg m-2)
  real :: H_geom        ! The geometric mean of Hup*Hdn, in m or kg m-2.
  real :: one_meter     ! One meter in thickness units of m or kg m-2.
  integer :: is, ie, js, je, nz
  integer :: i, j, k, kb_max
  real    :: SN_u_local(SZIB_(G), SZJ_(G),SZK_(G))
  real    :: SN_v_local(SZI_(G), SZJB_(G),SZK_(G))

  if (.not. ASSOCIATED(CS)) call MOM_error(FATAL, "calc_slope_function:"// &
         "Module must be initialized before it is used.")
  if (.not. CS%use_variable_mixing) return
  if (.not. ASSOCIATED(CS%SN_u)) call MOM_error(FATAL, "calc_slope_function:"// &
         "%SN_u is not associated with use_variable_mixing.")
  if (.not. ASSOCIATED(CS%SN_v)) call MOM_error(FATAL, "calc_slope_function:"// &
         "%SN_v is not associated with use_variable_mixing.")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  one_meter = 1.0 * GV%m_to_H
  h_neglect = GV%H_subroundoff
  H_cutoff = real(2*nz) * (GV%Angstrom + h_neglect)

!$OMP parallel default(none) shared(is,ie,js,je,CS,nz,e,G,GV,h,H_cutoff,h_neglect, &
!$OMP                               one_meter,SN_u_local,SN_v_local,calculate_slopes)   &
!$OMP                       private(E_x,E_y,S2,Hdn,Hup,H_geom,N2)
!$OMP do
  do j=js-1,je+1 ; do i=is-1,ie+1
    CS%SN_u(i,j) = 0.0
    CS%SN_v(i,j) = 0.0
  enddo ; enddo

  ! To set the length scale based on the deformation radius, use wave_speed to
  ! calculate the first-mode gravity wave speed and then blend the equatorial
  ! and midlatitude deformation radii, using calc_resoln_function as a template.

  ! Set the length scale at u-points.
!$OMP do
  do j=js,je ; do I=is-1,ie
    CS%L2u(I,j) = CS%Visbeck_L_scale**2
  enddo ; enddo
  ! Set length scale at v-points
!$OMP do
  do J=js-1,je ; do i=is,ie
    CS%L2v(i,J) = CS%Visbeck_L_scale**2
  enddo ; enddo
!$OMP do
  do k=nz,CS%VarMix_Ktop,-1

    if (calculate_slopes) then
      ! Calculate the interface slopes E_x and E_y and u- and v- points respectively
      do j=js-1,je+1 ; do I=is-1,ie
        E_x(I,j) = (e(i+1,j,K)-e(i,j,K))*G%IdxCu(I,j)
        ! Mask slopes where interface intersects topography
        if (min(h(I,j,k),h(I+1,j,k)) < H_cutoff) E_x(I,j) = 0.
      enddo ; enddo
      do J=js-1,je ; do i=is-1,ie+1
        E_y(i,J) = (e(i,j+1,K)-e(i,j,K))*G%IdyCv(i,J)
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
      N2 = GV%g_prime(k) / (GV%H_to_m * max(Hdn,Hup,one_meter))
      if (min(h(i,j,k-1), h(i+1,j,k-1), h(i,j,k), h(i+1,j,k)) < H_cutoff) &
        S2 = 0.0
      SN_u_local(I,j,k) = (H_geom * GV%H_to_m) * S2 * N2
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      S2 = ( E_y(i,J)**2  + 0.25*( &
            (E_x(i,J)**2+E_x(i-1,J+1)**2)+(E_x(i,J+1)**2+E_x(i-1,J)**2) ) )
      Hdn = 2.*h(i,j,k)*h(i,j,k-1) / (h(i,j,k) + h(i,j,k-1) + h_neglect)
      Hup = 2.*h(i,j+1,k)*h(i,j+1,k-1) / (h(i,j+1,k) + h(i,j+1,k-1) + h_neglect)
      H_geom = sqrt(Hdn*Hup)
      N2 = GV%g_prime(k) / (GV%H_to_m * max(Hdn,Hup,one_meter))
      if (min(h(i,j,k-1), h(i,j+1,k-1), h(i,j,k), h(i,j+1,k)) < H_cutoff) &
        S2 = 0.0
      SN_v_local(i,J,k) = (H_geom * GV%H_to_m) * S2 * N2
    enddo ; enddo

  enddo ! k
!$OMP do
  do j = js,je;
    do k=nz,CS%VarMix_Ktop,-1 ; do I=is-1,ie
      CS%SN_u(I,j) = CS%SN_u(I,j) + SN_u_local(I,j,k)
    enddo ; enddo
    ! SN above contains S^2*N^2*H, convert to vertical average of S*N
    do I=is-1,ie
      !SN_u(I,j) = sqrt( SN_u(I,j) / ( max(G%bathyT(I,j), G%bathyT(I+1,j)) + GV%Angstrom ) )
      !The code below behaves better than the line above. Not sure why? AJA
      if ( min(G%bathyT(I,j), G%bathyT(I+1,j)) > H_cutoff ) then
        CS%SN_u(I,j) = sqrt( CS%SN_u(I,j) / max(G%bathyT(I,j), G%bathyT(I+1,j)) )
      else
        CS%SN_u(I,j) = 0.0
      endif
    enddo
  enddo
!$OMP do
  do J=js-1,je
    do k=nz,CS%VarMix_Ktop,-1 ; do I=is,ie
      CS%SN_v(i,J) = CS%SN_v(i,J) + SN_v_local(i,J,k)
    enddo ; enddo
    do i=is,ie
      !SN_v(i,J) = sqrt( SN_v(i,J) / ( max(G%bathyT(i,J), G%bathyT(i,J+1)) + GV%Angstrom ) )
      !The code below behaves better than the line above. Not sure why? AJA
      if ( min(G%bathyT(I,j), G%bathyT(I+1,j)) > H_cutoff ) then
        CS%SN_v(i,J) = sqrt( CS%SN_v(i,J) / max(G%bathyT(i,J), G%bathyT(i,J+1)) )
      else
        CS%SN_v(I,j) = 0.0
      endif
    enddo
  enddo
!$OMP end parallel

end subroutine calc_slope_functions_using_just_e

!> Initializes the variables mixing coefficients container
subroutine VarMix_init(Time, G, param_file, diag, CS)
  type(time_type),            intent(in) :: Time !< Current model time
  type(ocean_grid_type),      intent(in) :: G    !< Ocean grid structure
  type(param_file_type),      intent(in) :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(VarMix_CS),               pointer :: CS   !< Variable mixing coefficients
  ! Local variables
  real :: KhTr_Slope_Cff, KhTh_Slope_Cff, oneOrTwo, N2_filter_depth
  real, parameter :: absurdly_small_freq2 = 1e-34  ! A miniscule frequency
             ! squared that is used to avoid division by 0, in s-2.  This
             ! value is roughly (pi / (the age of the universe) )^2.
  logical :: use_variable_mixing, Gill_equatorial_Ld, use_stored_slopes
  logical :: Resoln_scaled_Kh, Resoln_scaled_KhTh, Resoln_scaled_KhTr
  logical :: Resoln_use_ebt, khth_use_ebt_struct
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_lateral_mixing_coeffs" ! This module's name.
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

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  !   This first set of parameters are read into local variables first, in case
  ! the control structure should not be allocated.
  call get_param(param_file, mod, "USE_VARIABLE_MIXING", use_variable_mixing,&
                 "If true, the variable mixing code will be called.  This \n"//&
                 "allows diagnostics to be created even if the scheme is \n"//&
                 "not used.  If KHTR_SLOPE_CFF>0 or  KhTh_Slope_Cff>0, \n"//&
                 "this is set to true regardless of what is in the \n"//&
                 "parameter file.", default=.false.)
  call get_param(param_file, mod, "RESOLN_SCALED_KH", Resoln_scaled_Kh, &
                 "If true, the Laplacian lateral viscosity is scaled away \n"//&
                 "when the first baroclinic deformation radius is well \n"//&
                 "resolved.", default=.false.)
  call get_param(param_file, mod, "RESOLN_SCALED_KHTH", Resoln_scaled_KhTh, &
                 "If true, the interface depth diffusivity is scaled away \n"//&
                 "when the first baroclinic deformation radius is well \n"//&
                 "resolved.", default=.false.)
  call get_param(param_file, mod, "RESOLN_SCALED_KHTR", Resoln_scaled_KhTr, &
                 "If true, the epipycnal tracer diffusivity is scaled \n"//&
                 "away when the first baroclinic deformation radius is \n"//&
                 "well resolved.", default=.false.)
  call get_param(param_file, mod, "RESOLN_USE_EBT", Resoln_use_ebt, &
                 "If true, uses the equivalent barotropic wave speed instead\n"//&
                 "of first baroclinic wave for calculating the resolution fn.",&
                 default=.false.)
  call get_param(param_file, mod, "KHTH_USE_EBT_STRUCT", khth_use_ebt_struct, &
                 "If true, uses the equivalent barotropic structure\n"//&
                 "as the vertical structure of thickness diffusivity.",&
                 default=.false.)
  call get_param(param_file, mod, "KHTH_SLOPE_CFF", KhTh_Slope_Cff, &
                 "The nondimensional coefficient in the Visbeck formula \n"//&
                 "for the interface depth diffusivity", units="nondim", &
                 default=0.0)
  call get_param(param_file, mod, "KHTR_SLOPE_CFF", KhTr_Slope_Cff, &
                 "The nondimensional coefficient in the Visbeck formula \n"//&
                 "for the epipycnal tracer diffusivity", units="nondim", &
                 default=0.0)
  call get_param(param_file, mod, "USE_STORED_SLOPES", use_stored_slopes,&
                 "If true, the isopycnal slopes are calculated once and\n"//&
                 "stored for re-use. This uses more memory but avoids calling\n"//&
                 "the equation of state more times than should be necessary.", &
                 default=.false.)
  if (KhTr_Slope_Cff>0. .or. KhTh_Slope_Cff>0.) use_variable_mixing = .true.

  if (use_variable_mixing .or. Resoln_scaled_Kh .or. Resoln_scaled_KhTh .or. &
      Resoln_scaled_KhTr .or. use_stored_slopes .or. khth_use_ebt_struct) then
    allocate(CS)
    CS%diag => diag ! Diagnostics pointer
    CS%Resoln_scaled_Kh = Resoln_scaled_Kh
    CS%Resoln_scaled_KhTh = Resoln_scaled_KhTh
    CS%Resoln_scaled_KhTr = Resoln_scaled_KhTr
    CS%Resoln_use_ebt = Resoln_use_ebt
    CS%khth_use_ebt_struct = khth_use_ebt_struct
    CS%use_variable_mixing = use_variable_mixing
    CS%use_stored_slopes = use_stored_slopes
  else
    return
  endif
  if (Resoln_use_ebt .or. khth_use_ebt_struct) then
    call get_param(param_file, mod, "RESOLN_N2_FILTER_DEPTH", N2_filter_depth, &
                 "The depth below which N2 is monotonized to avoid stratification\n"//&
                 "artifacts from altering the equivalent barotropic mode structure.",&
                 units='m', default=2000.)
  endif
  if (khth_use_ebt_struct) then
    allocate(CS%ebt_struct(isd:ied,JsdB:JedB,G%ke)) ; CS%ebt_struct(:,:,:) = 0.0
  endif
  if (use_variable_mixing) then
    call get_param(param_file, mod, "VISBECK_MAX_SLOPE", CS%Visbeck_S_max, &
          "If non-zero, is an upper bound on slopes used in the\n"//       &
          "Visbeck formula for diffusivity. This does not affect the\n"//  &
          "isopycnal slope calculation used within thickness diffusion.",  &
          units="nondim", default=0.0)
  endif

! Allocate CS and memory
  if (CS%use_stored_slopes) then
    allocate(CS%slope_x(IsdB:IedB,jsd:jed,G%ke+1)) ; CS%slope_x(:,:,:) = 0.0
    allocate(CS%slope_y(isd:ied,JsdB:JedB,G%ke+1)) ; CS%slope_y(:,:,:) = 0.0
    call get_param(param_file, mod, "KD_SMOOTH", CS%kappa_smooth, &
                 "A diapycnal diffusivity that is used to interpolate \n"//&
                 "more sensible values of T & S into thin layers.", &
                 default=1.0e-6)
  endif

  if (CS%use_variable_mixing) then
    allocate(CS%SN_u(IsdB:IedB,jsd:jed)) ; CS%SN_u(:,:) = 0.0
    allocate(CS%SN_v(isd:ied,JsdB:JedB)) ; CS%SN_v(:,:) = 0.0
    allocate(CS%L2u(IsdB:IedB,jsd:jed)) ; CS%L2u(:,:) = 0.0
    allocate(CS%L2v(isd:ied,JsdB:JedB)) ; CS%L2v(:,:) = 0.0
    call MOM_mesg("VarMix_init: memory allocated for use_variable_mixing", 5)

  ! More run-time parameters
    call get_param(param_file, mod, "VARMIX_KTOP", CS%VarMix_Ktop, &
                 "The layer number at which to start vertical integration \n"//&
                 "of S*N for purposes of finding the Eady growth rate.", &
                 units="nondim", default=2)
    call get_param(param_file, mod, "VISBECK_L_SCALE", CS%Visbeck_L_scale, &
                 "The fixed length scale in the Visbeck formula.", units="m", &
                 default=0.0)

  ! Register fields for output from this module.
    CS%id_SN_u = register_diag_field('ocean_model', 'SN_u', diag%axesCu1, Time, &
       'Inverse eddy time-scale, S*N, at u-points', 's^-1')
    CS%id_SN_v = register_diag_field('ocean_model', 'SN_v', diag%axesCv1, Time, &
       'Inverse eddy time-scale, S*N, at v-points', 's^-1')
    CS%id_L2u = register_diag_field('ocean_model', 'L2u', diag%axesCu1, Time, &
       'Length scale squared for mixing coefficient, at u-points', 'm^2')
    CS%id_L2v = register_diag_field('ocean_model', 'L2v', diag%axesCv1, Time, &
       'Length scale squared for mixing coefficient, at v-points', 'm^2')

    if (CS%use_stored_slopes) then
      CS%id_N2_u = register_diag_field('ocean_model', 'N2_u', diag%axesCui, Time, &
         'Square of Brunt-Vaisala frequency, N^2, at u-points, as used in Visbeck et al.', 's^-2')
      CS%id_N2_v = register_diag_field('ocean_model', 'N2_v', diag%axesCvi, Time, &
         'Square of Brunt-Vaisala frequency, N^2, at v-points, as used in Visbeck et al.', 's^-2')
      CS%id_S2_u = register_diag_field('ocean_model', 'S2_u', diag%axesCu1, Time, &
         'Depth average square of slope magnitude, S^2, at u-points, as used in Visbeck et al.', 's^-2')
      CS%id_S2_v = register_diag_field('ocean_model', 'S2_v', diag%axesCv1, Time, &
         'Depth average square of slope magnitude, S^2, at v-points, as used in Visbeck et al.', 's^-2')
    endif
  endif

  if (CS%Resoln_scaled_Kh .or. Resoln_scaled_KhTh .or. Resoln_scaled_KhTr) then
    call wave_speed_init(CS%wave_speed_CSp, use_ebt_mode=Resoln_use_ebt, mono_N2_depth=N2_filter_depth)

    ! Allocate and initialize various arrays.
    allocate(CS%Res_fn_h(isd:ied,jsd:jed))       ; CS%Res_fn_h(:,:) = 0.0
    allocate(CS%Res_fn_q(IsdB:IedB,JsdB:JedB))   ; CS%Res_fn_q(:,:) = 0.0
    allocate(CS%Res_fn_u(IsdB:IedB,jsd:jed))     ; CS%Res_fn_u(:,:) = 0.0
    allocate(CS%Res_fn_v(isd:ied,JsdB:JedB))     ; CS%Res_fn_v(:,:) = 0.0
    allocate(CS%cg1(isd:ied,jsd:jed))            ; CS%cg1(:,:) = 0.0
    allocate(CS%beta_dx2_h(isd:ied,jsd:jed))     ; CS%beta_dx2_h(:,:) = 0.0
    allocate(CS%beta_dx2_q(IsdB:IedB,JsdB:JedB)) ; CS%beta_dx2_q(:,:) = 0.0
    allocate(CS%beta_dx2_u(IsdB:IedB,jsd:jed))   ; CS%beta_dx2_u(:,:) = 0.0
    allocate(CS%beta_dx2_v(isd:ied,JsdB:JedB))   ; CS%beta_dx2_v(:,:) = 0.0
    allocate(CS%f2_dx2_h(isd:ied,jsd:jed))       ; CS%f2_dx2_h(:,:) = 0.0
    allocate(CS%f2_dx2_q(IsdB:IedB,JsdB:JedB))   ; CS%f2_dx2_q(:,:) = 0.0
    allocate(CS%f2_dx2_u(IsdB:IedB,jsd:jed))     ; CS%f2_dx2_u(:,:) = 0.0
    allocate(CS%f2_dx2_v(isd:ied,JsdB:JedB))     ; CS%f2_dx2_v(:,:) = 0.0
    allocate(CS%Rd_dx_h(isd:ied,jsd:jed))        ; CS%Rd_dx_h(:,:) = 0.0

    CS%id_Res_fn = register_diag_field('ocean_model', 'Res_fn', diag%axesT1, Time, &
       'Resolution function for scaling diffusivities', 'Nondim')
    CS%id_Rd_dx = register_diag_field('ocean_model', 'Rd_dx', diag%axesT1, Time, &
       'Ratio between deformation radius and grid spacing', 'Nondim')

    call get_param(param_file, mod, "KH_RES_SCALE_COEF", CS%Res_coef_khth, &
                 "A coefficient that determines how KhTh is scaled away if \n"//&
                 "RESOLN_SCALED_... is true, as \n"//&
                 "F = 1 / (1 + (KH_RES_SCALE_COEF*Rd/dx)^KH_RES_FN_POWER).", &
                 units="nondim", default=1.0)
    call get_param(param_file, mod, "KH_RES_FN_POWER", CS%Res_fn_power_khth, &
                 "The power of dx/Ld in the Kh resolution function.  Any \n"//&
                 "positive integer may be used, although even integers \n"//&
                 "are more efficient to calculate.  Setting this greater \n"//&
                 "than 100 results in a step-function being used.", &
                 units="nondim", default=2)
    call get_param(param_file, mod, "VISC_RES_SCALE_COEF", CS%Res_coef_visc, &
                 "A coefficient that determines how Kh is scaled away if \n"//&
                 "RESOLN_SCALED_... is true, as \n"//&
                 "F = 1 / (1 + (KH_RES_SCALE_COEF*Rd/dx)^KH_RES_FN_POWER).\n"//&
                 "This function affects lateral viscosity, Kh, and not KhTh.", &
                 units="nondim", default=CS%Res_coef_khth)
    call get_param(param_file, mod, "VISC_RES_FN_POWER", CS%Res_fn_power_visc, &
                 "The power of dx/Ld in the Kh resolution function.  Any \n"//&
                 "positive integer may be used, although even integers \n"//&
                 "are more efficient to calculate.  Setting this greater \n"//&
                 "than 100 results in a step-function being used.\n"//&
                 "This function affects lateral viscosity, Kh, and not KhTh.", &
                 units="nondim", default=CS%Res_fn_power_khth)
    call get_param(param_file, mod, "INTERPOLATE_RES_FN", CS%interpolate_Res_fn, &
                 "If true, interpolate the resolution function to the \n"//&
                 "velocity points from the thickness points; otherwise \n"//&
                 "interpolate the wave speed and calculate the resolution \n"//&
                 "function independently at each point.", default=.true.)
    if (CS%interpolate_Res_fn) then
      if (CS%Res_coef_visc .ne. CS%Res_coef_khth) call MOM_error(FATAL, &
           "MOM_lateral_mixing_coeffs.F90, VarMix_init:"//&
           "When INTERPOLATE_RES_FN=True, VISC_RES_FN_POWER must equal KH_RES_SCALE_COEF.")
      if (CS%Res_fn_power_visc .ne. CS%Res_fn_power_khth) call MOM_error(FATAL, &
           "MOM_lateral_mixing_coeffs.F90, VarMix_init:"//&
           "When INTERPOLATE_RES_FN=True, VISC_RES_FN_POWER must equal KH_RES_FN_POWER.")
    endif
    call get_param(param_file, mod, "GILL_EQUATORIAL_LD", Gill_equatorial_Ld, &
                 "If true, uses Gill's definition of the baroclinic\n"//&
                 "equatorial deformation radius, otherwise, if false, use\n"//&
                 "Pedlosky's definition. These definitions differ by a factor\n"//&
                 "of 2 infront of the beta term in the denominator. Gill's"//&
                 "is the more appropriate definition.\n", default=.false.)

    ! Pre-calculate several static expressions for later use.
    if (Gill_equatorial_Ld) then; oneOrTwo = 2.0
      else; oneOrTwo = 1.0; endif

    do j=js-1,je+1 ; do i=is-1,ie+1
      CS%f2_dx2_h(i,j) = (G%dxT(i,j)**2 + G%dyT(i,j)**2) * &
          max(0.25 * ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                      (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2)), &
              absurdly_small_freq2)
      CS%beta_dx2_h(i,j) = oneOrTwo * (G%dxT(i,j)**2 + G%dyT(i,j)**2) * (sqrt(0.5 * &
          ( (((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
             ((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2) + &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
             ((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2) ) ))
    enddo ; enddo

    do J=js-1,Jeq ; do I=is-1,Ieq
      CS%f2_dx2_q(I,J) = (G%dxBu(I,J)**2 + G%dyBu(I,J)**2) * &
                         max(G%CoriolisBu(I,J)**2, absurdly_small_freq2)
      CS%beta_dx2_q(I,J) = oneOrTwo * (G%dxBu(I,J)**2 + G%dyBu(I,J)**2) * (sqrt(0.5 * &
          ( (((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
             ((G%CoriolisBu(I+1,J)-G%CoriolisBu(I,J)) * G%IdxCv(i+1,J))**2) + &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
             ((G%CoriolisBu(I,J+1)-G%CoriolisBu(I,J)) * G%IdyCu(I,j+1))**2) ) ))
    enddo ; enddo

    do j=js,je ; do I=is-1,Ieq
      CS%f2_dx2_u(I,j) = (G%dxCu(I,j)**2 + G%dyCu(I,j)**2) * &
          max(0.5*(G%CoriolisBu(I,J)**2+G%CoriolisBu(I,J-1)**2), absurdly_small_freq2)
      CS%beta_dx2_u(I,j) = oneOrTwo * (G%dxCu(I,j)**2 + G%dyCu(I,j)**2) * (sqrt( &
          0.25*( (((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2 + &
                  ((G%CoriolisBu(I+1,J)-G%CoriolisBu(I,J)) * G%IdxCv(i+1,J))**2) + &
                 (((G%CoriolisBu(I+1,J-1)-G%CoriolisBu(I,J-1)) * G%IdxCv(i+1,J-1))**2 + &
                  ((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2) ) + &
                  ((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 ))
    enddo ; enddo

    do J=js-1,Jeq ; do i=is,ie
      CS%f2_dx2_v(i,J) = (G%dxCv(i,J)**2 + G%dyCv(i,J)**2) * &
          max(0.5*(G%CoriolisBu(I,J)**2+G%CoriolisBu(I-1,J)**2), absurdly_small_freq2)
      CS%beta_dx2_v(i,J) = oneOrTwo * (G%dxCv(i,J)**2 + G%dyCv(i,J)**2) * (sqrt( &
          ((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
          0.25*( (((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
                  ((G%CoriolisBu(I-1,J+1)-G%CoriolisBu(I-1,J)) * G%IdyCu(I-1,j+1))**2) + &
                 (((G%CoriolisBu(I,J+1)-G%CoriolisBu(I,J)) * G%IdyCu(I,j+1))**2 + &
                  ((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2) ) ))
    enddo ; enddo

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
!! The resolution function can be applied independently to thickness diffusion (module mom_thickness_diffuse), tracer diffusion (mom_tracer_hordiff)
!! lateral viscosity (mom_hor_visc).
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
!! This module also calculates factors used in setting the thickness diffusivity similar to a Visbeck et al., 1997, scheme.
!! The factors are combined in mom_thickness_diffuse::thickness_diffuse() but calculated in this module.
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
!! The thickness diffusivity can be prescribed a vertical distribution with the shape of the equivalent barotropic velocity mode.
!! The structure function is stored in the control structure for thie module (varmix_cs) but is calculated use subroutines in
!! mom_wave_speed.
!!
!! | Symbol                | Module parameter |
!! | ------                | --------------- |
!! | -                     | <code>KHTH_USE_EBT_STRUCT</code> |

end module MOM_lateral_mixing_coeffs
