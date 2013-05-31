module MOM_lateral_mixing_coeffs
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
!*    This program contains the subroutine that calculates the         *
!*  effects of horizontal viscosity, including parameterizations of    *
!*  the value of the viscosity itself. mesoscale_EKE calculates        *
!*  the evolution of sub-grid scale mesoscale EKE.                     *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, post_data
use MOM_diag_mediator, only : diag_ptrs, time_type, query_averaging_enabled
use MOM_domains, only : pass_var, pass_vector, CGRID_NE, To_All, Scalar_Pair
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_interface_heights, only : find_eta
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_wave_speed, only : wave_speed_init, wave_speed, wave_speed_CS

implicit none ; private

#include <MOM_memory.h>

type, public :: VarMix_CS ;
  logical :: use_variable_mixing  ! If true, use the variable mixing.
  logical :: Resoln_scaled_Kh     ! If true, scale away the Laplacian viscosity
                                  ! when the deformation radius is well resolved.
  logical :: Resoln_scaled_KhTh   ! If true, scale away the thickness diffusivity
                                  ! when the deformation radius is well resolved.
  logical :: Resoln_scaled_KhTr   ! If true, scale away the tracer diffusivity
                                  ! when the deformation radius is well resolved.
  real, dimension(:,:), pointer :: &
    SN_u => NULL(), &   ! S*N at u-points (s^-1)
    SN_v => NULL(), &  ! S*N at v-points (s^-1)
    L2u => NULL(), &   ! Length scale^2 at u-points (m^2)
    L2v => NULL(), &   ! Length scale^2 at v-points (m^2)
    cg1 => NULL(), &   ! The first baroclinic gravity wave speed in m s-1.
    Res_fn_h => NULL(), & ! Res_fn_h and Res_fn_q are nondimensional functions
    Res_fn_q => NULL(), & ! of the ratio the first baroclinic deformation
                          ! radius to the grid spacing at h and q points,
                          ! respectively. These can be used to scale away
                          ! horizontal viscosities or thickness diffusivities
                          ! when the deformation radius is well resolved.
    beta_dx2_h => NULL(), &  ! The magnitude of the gradient of the Coriolis
    beta_dx2_q => NULL(), &  ! parameter times the grid spacing squared at
                             ! h and q points, in m s-1.
    f2_dx2_h => NULL(), & ! The Coriolis parameter squared times the grid
    f2_dx2_q => NULL(), & ! spacing squared at h and q points, in m2 s-2.
    Rd_dx_h => NULL()     ! Deformation radius over grid spacing (non-dim.)
  ! Parameters
  integer :: VarMix_Ktop  ! Top layer to start downward integrals
  real :: Visbeck_L_scale ! Fixed length scale in Visbeck formula
  real :: Res_coef        ! A nondimensional number that determines the function
                          ! of resolution as:
                          !  F = 1 / (1 + (Res_coef*Ld/dx)^Res_fn_power)
                          
  integer :: Res_fn_power ! The power of dx/Ld in the resolution function.  Any
                          ! positive integer power may be used, but even powers
                          ! and especially 2 are coded to be more efficient.

  ! Diagnostics
  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
  type(wave_speed_CS), pointer :: wave_speed_CSp => NULL()
  integer :: id_SN_u=-1, id_SN_v=-1, id_L2u=-1, id_L2v=-1, id_Res_fn = -1
  integer :: id_Rd_dx=-1
end type VarMix_CS

public VarMix_init, calc_slope_function, calc_resoln_function

interface calc_slope_function
  module procedure calc_slope_function_, calc_slope_function_need_e
end interface calc_slope_function

contains

subroutine calc_resoln_function(h, tv, G, CS)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h
  type(thermo_var_ptrs),                 intent(in)    :: tv
  type(ocean_grid_type),                 intent(inout) :: G
  type(VarMix_CS),                       pointer       :: CS
!    This subroutine determines a function of the ratio of the grid
! spacing to the deformation radius that is used to scale horizontal
! viscosities or diffusivities.
!   In the current implementation, the function is F(R) = 1/(1+R) where
! R = cg^2/(f^2*(dx^2+dy^2))
! Arguments: h - Layer thickness, in m or kg m-2.
!  (in)      tv - A structure containing the thermodynamic variables.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 VarMix_init.

  real :: cg1_q  ! The gravity wave speed interpolated to q points, in m s-1.
  real :: dx_term
  integer :: mod_power_2, power_2
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
  if (.not. ASSOCIATED(CS%f2_dx2_h)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_h is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%beta_dx2_h)) call MOM_error(FATAL, &
    "calc_resoln_function: %beta_dx2_h is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%f2_dx2_q)) call MOM_error(FATAL, &
    "calc_resoln_function: %f2_dx2_q is not associated with Resoln_scaled_Kh.")
  if (.not. ASSOCIATED(CS%beta_dx2_q)) call MOM_error(FATAL, &
    "calc_resoln_function: %beta_dx2_q is not associated with Resoln_scaled_Kh.")

  call wave_speed(h, tv, G, CS%cg1, CS%wave_speed_CSp)

  call pass_var(CS%cg1, G%Domain)

  !   Do this calculation on the extent used in MOM_hor_visc.F90, and
  ! MOM_tracer.F90 so that no halo update is needed.
  mod_power_2 = mod(CS%Res_fn_power, 2)
  if (CS%Res_fn_power >= 100) then
    do j=js-1,je+1 ; do i=is-1,ie+1
      dx_term = CS%f2_dx2_h(i,j) + CS%cg1(i,j)*CS%beta_dx2_h(i,j)
      if ((CS%Res_coef * CS%cg1(i,j))**2 > dx_term) then
        CS%Res_fn_h(i,j) = 0.0
      else
        CS%Res_fn_h(i,j) = 1.0
      endif 
    enddo ; enddo

    do J=js-1,Jeq ; do I=is-1,Ieq
      cg1_q = 0.25 * ((CS%cg1(i,j) + CS%cg1(i+1,j+1)) + &
                      (CS%cg1(i+1,j) + CS%cg1(i,j+1)))
      dx_term = CS%f2_dx2_q(I,J) +  cg1_q * CS%beta_dx2_q(I,J)
      if ((CS%Res_coef * cg1_q)**2 > dx_term) then
        CS%Res_fn_q(I,J) = 0.0
      else
        CS%Res_fn_q(I,J) = 1.0
      endif
    enddo ; enddo
  elseif (CS%Res_fn_power == 2) then
    do j=js-1,je+1 ; do i=is-1,ie+1
      dx_term = CS%f2_dx2_h(i,j) + CS%cg1(i,j)*CS%beta_dx2_h(i,j)
      CS%Res_fn_h(i,j) = dx_term / (dx_term + (CS%Res_coef * CS%cg1(i,j))**2)
    enddo ; enddo

    do J=js-1,Jeq ; do I=is-1,Ieq
      cg1_q = 0.25 * ((CS%cg1(i,j) + CS%cg1(i+1,j+1)) + &
                      (CS%cg1(i+1,j) + CS%cg1(i,j+1)))
      dx_term = CS%f2_dx2_q(I,J) +  cg1_q * CS%beta_dx2_q(I,J)
      CS%Res_fn_q(I,J) = dx_term / (dx_term + (CS%Res_coef * cg1_q)**2)
    enddo ; enddo
  elseif (mod_power_2 == 0) then
    power_2 = CS%Res_fn_power / 2
    do j=js-1,je+1 ; do i=is-1,ie+1
      dx_term = (CS%f2_dx2_h(i,j) + CS%cg1(i,j)*CS%beta_dx2_h(i,j))**power_2
      CS%Res_fn_h(i,j) = dx_term / &
          (dx_term + (CS%Res_coef * CS%cg1(i,j))**CS%Res_fn_power)
    enddo ; enddo

    do J=js-1,Jeq ; do I=is-1,Ieq
      cg1_q = 0.25 * ((CS%cg1(i,j) + CS%cg1(i+1,j+1)) + &
                      (CS%cg1(i+1,j) + CS%cg1(i,j+1)))
      dx_term = (CS%f2_dx2_q(I,J) +  cg1_q * CS%beta_dx2_q(I,J))**power_2
      CS%Res_fn_q(I,J) = dx_term / &
          (dx_term + (CS%Res_coef * cg1_q)**CS%Res_fn_power)
    enddo ; enddo
  else
    do j=js-1,je+1 ; do i=is-1,ie+1
      dx_term = (sqrt(CS%f2_dx2_h(i,j) + &
                      CS%cg1(i,j)*CS%beta_dx2_h(i,j)))**CS%Res_fn_power
      CS%Res_fn_h(i,j) = dx_term / &
         (dx_term + (CS%Res_coef * CS%cg1(i,j))**CS%Res_fn_power)
    enddo ; enddo


    do J=js-1,Jeq ; do I=is-1,Ieq
      cg1_q = 0.25 * ((CS%cg1(i,j) + CS%cg1(i+1,j+1)) + &
                      (CS%cg1(i+1,j) + CS%cg1(i,j+1)))
      dx_term = (sqrt(CS%f2_dx2_q(I,J) + &
                      cg1_q * CS%beta_dx2_q(I,J)))**CS%Res_fn_power
      CS%Res_fn_q(I,J) = dx_term / &
          (dx_term + (CS%Res_coef * cg1_q)**CS%Res_fn_power)
    enddo ; enddo
  endif

  ! Calculate and store the ratio between deformation radius and grid-spacing
  ! at h-points (non-dimensional).
  do j=js-1,je+1 ; do i=is-1,ie+1
    CS%Rd_dx_h(i,j) = CS%cg1(i,j) / &
          (sqrt(CS%f2_dx2_h(i,j) + CS%cg1(i,j)*CS%beta_dx2_h(i,j)))
  enddo ; enddo

  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_Res_fn > 0) call post_data(CS%id_Res_fn, CS%Res_fn_h, CS%diag)
    if (CS%id_Rd_dx > 0) call post_data(CS%id_Rd_dx, CS%Rd_dx_h, CS%diag)
  endif

end subroutine calc_resoln_function

subroutine calc_slope_function_need_e(h, tv, G, CS)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: h
  type(thermo_var_ptrs),                 intent(in)    :: tv
  type(ocean_grid_type),                 intent(inout) :: G
  type(VarMix_CS),                       pointer       :: CS
!    This subroutine calls for the calculation of the interface heights, and
!  then calls for the calculation of the slope function S*N for the Visbeck et
!  al. style scaling for the various horizontal diffusivities.
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)+1) :: &
    e             ! The interface heights relative to mean sea level, in m.

  call find_eta(h, tv, G%g_Earth, G, e, halo_size=1)
  
  call calc_slope_function_(h, tv, G, CS, e)

end subroutine calc_slope_function_need_e

subroutine calc_slope_function_(h, tv, G, CS, e)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: h
  type(thermo_var_ptrs),                 intent(in)    :: tv
  type(ocean_grid_type),                 intent(inout) :: G
  type(VarMix_CS),                       pointer       :: CS
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_), intent(in) :: e
!    This subroutine calculates the slope function S*N for the Visbeck et
!  al. style scaling for the various horizontal diffusivities.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      G - The ocean's grid structure.
!  (inout)   CS - The control structure returned by a previous call to
!                 thickness_diffuse_init.
!  (in)      e - The interface heights, in m.
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

  if (.not. ASSOCIATED(CS)) call MOM_error(FATAL, "calc_slope_function:"// &
         "Module must be initialized before it is used.")
  if (.not. CS%use_variable_mixing) return
  if (.not. ASSOCIATED(CS%SN_u)) call MOM_error(FATAL, "calc_slope_function:"// &
         "%SN_u is not associated with use_variable_mixing.")
  if (.not. ASSOCIATED(CS%SN_v)) call MOM_error(FATAL, "calc_slope_function:"// &
         "%SN_v is not associated with use_variable_mixing.")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  one_meter = 1.0 * G%m_to_H
  h_neglect = G%H_subroundoff
  H_cutoff = real(2*nz) * (G%Angstrom + h_neglect)

  do j=js-1,je+1 ; do i=is-1,ie+1
    CS%SN_u(i,j) = 0.0
    CS%SN_v(i,j) = 0.0
  enddo ; enddo

  ! To set the length scale based on the deformation radius, use wave_speed to
  ! calculate the first-mode gravity wave speed and then blend the equatorial
  ! and midlatitude deformation radii, using calc_resoln_function as a template.

  ! Set the length scale at u-points.
  do j=js,je ; do I=is-1,ie
    CS%L2u(I,j) = CS%Visbeck_L_scale**2
  enddo ; enddo
  ! Set length scale at v-points
  do J=js-1,je ; do i=is,ie
    CS%L2v(i,J) = CS%Visbeck_L_scale**2
  enddo ; enddo

  do k=nz,CS%VarMix_Ktop,-1

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

    ! Calculate N*S*h from this layer and add to the sum
    do j=js,je ; do I=is-1,ie
      S2 =  ( E_x(I,j)**2  + 0.25*( &
            (E_y(I,j)**2+E_y(I+1,j-1)**2)+(E_y(I+1,j)**2+E_y(I,j-1)**2) ) )
      Hdn = 2.*h(i,j,k)*h(i,j,k-1) / (h(i,j,k) + h(i,j,k-1) + h_neglect)
      Hup = 2.*h(i+1,j,k)*h(i+1,j,k-1) / (h(i+1,j,k) + h(i+1,j,k-1) + h_neglect)
      H_geom = sqrt(Hdn*Hup)
      N2 = G%g_prime(k) / (G%H_to_m * max(Hdn,Hup,one_meter))
      if (min(h(i,j,k-1), h(i+1,j,k-1), h(i,j,k), h(i+1,j,k)) < H_cutoff) &
        S2 = 0.0
      CS%SN_u(I,j) = CS%SN_u(I,j) + (H_geom * G%H_to_m) * S2 * N2
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      S2 = ( E_y(i,J)**2  + 0.25*( &
            (E_x(i,J)**2+E_x(i-1,J+1)**2)+(E_x(i,J+1)**2+E_x(i-1,J)**2) ) )
      Hdn = 2.*h(i,j,k)*h(i,j,k-1) / (h(i,j,k) + h(i,j,k-1) + h_neglect)
      Hup = 2.*h(i,j+1,k)*h(i,j+1,k-1) / (h(i,j+1,k) + h(i,j+1,k-1) + h_neglect)
      H_geom = sqrt(Hdn*Hup)
      N2 = G%g_prime(k) / (G%H_to_m * max(Hdn,Hup,one_meter))
      if (min(h(i,j,k-1), h(i,j+1,k-1), h(i,j,k), h(i,j+1,k)) < H_cutoff) &
        S2 = 0.0
      CS%SN_v(i,J) = CS%SN_v(i,J) + (H_geom * G%H_to_m) * S2 * N2
    enddo ; enddo

  enddo ! k

  ! SN above contains S^2*N^2*H, convert to vertical average of S*N
  do j=js,je ; do I=is-1,ie
   !SN_u(I,j) = sqrt( SN_u(I,j) / ( max(G%bathyT(I,j), G%bathyT(I+1,j)) + G%Angstrom ) )
   !The code below behaves better than the line above. Not sure why? AJA
    if ( min(G%bathyT(I,j), G%bathyT(I+1,j)) > H_cutoff ) then
      CS%SN_u(I,j) = sqrt( CS%SN_u(I,j) / max(G%bathyT(I,j), G%bathyT(I+1,j)) )
    else
      CS%SN_u(I,j) = 0.0
    endif
  enddo ; enddo
  do J=js-1,je ; do i=is,ie
   !SN_v(i,J) = sqrt( SN_v(i,J) / ( max(G%bathyT(i,J), G%bathyT(i,J+1)) + G%Angstrom ) )
   !The code below behaves better than the line above. Not sure why? AJA
    if ( min(G%bathyT(I,j), G%bathyT(I+1,j)) > H_cutoff ) then
      CS%SN_v(i,J) = sqrt( CS%SN_v(i,J) / max(G%bathyT(i,J), G%bathyT(i,J+1)) )
    else
      CS%SN_v(I,j) = 0.0
    endif
  enddo ; enddo

! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_SN_u > 0) call post_data(CS%id_SN_u, CS%SN_u, CS%diag)
    if (CS%id_SN_v > 0) call post_data(CS%id_SN_v, CS%SN_v, CS%diag)
    if (CS%id_L2u > 0) call post_data(CS%id_L2u, CS%L2u, CS%diag)
    if (CS%id_L2v > 0) call post_data(CS%id_L2v, CS%L2v, CS%diag)
  endif

end subroutine calc_slope_function_

subroutine VarMix_init(Time, G, param_file, diag, CS)
  type(time_type),            intent(in) :: Time
  type(ocean_grid_type),      intent(in) :: G
  type(param_file_type),      intent(in) :: param_file
  type(diag_ptrs), target, intent(inout) :: diag
  type(VarMix_CS),               pointer :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
  real :: KhTr_Slope_Cff, KhTh_Slope_Cff
  real, parameter :: absurdly_small_freq2 = 1e-34  ! A miniscule frequency
             ! squared that is used to avoid division by 0, in s-2.  This
             ! value is roughly (pi / (the age of the universe) )^2.
  logical :: use_variable_mixing
  logical :: Resoln_scaled_Kh, Resoln_scaled_KhTh, Resoln_scaled_KhTr
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
  call get_param(param_file, mod, "KHTH_SLOPE_CFF", KhTh_Slope_Cff, &
                 "The nondimensional coefficient in the Visbeck formula \n"//&
                 "for the interface depth diffusivity", units="nondim", &
                 default=0.0)
  call get_param(param_file, mod, "KHTR_SLOPE_CFF", KhTr_Slope_Cff, &
                 "The nondimensional coefficient in the Visbeck formula \n"//&
                 "for the epipycnal tracer diffusivity", units="nondim", &
                 default=0.0)

  if (KhTr_Slope_Cff>0. .or. KhTh_Slope_Cff>0.) use_variable_mixing = .true.

  if (use_variable_mixing .or. Resoln_scaled_Kh .or. Resoln_scaled_KhTh .or. &
      Resoln_scaled_KhTr) then
    allocate(CS)
    CS%diag => diag ! Diagnostics pointer
    CS%Resoln_scaled_Kh = Resoln_scaled_Kh
    CS%Resoln_scaled_KhTh = Resoln_scaled_KhTh
    CS%Resoln_scaled_KhTr = Resoln_scaled_KhTr
    CS%use_variable_mixing = use_variable_mixing
  else
    return
  endif

! Allocate CS and memory
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
    CS%id_SN_u = register_diag_field('ocean_model', 'SN_u', G%axesCu1, Time, &
       'Inverse eddy time-scale, S*N, at u-points', 's^-1')
    CS%id_SN_v = register_diag_field('ocean_model', 'SN_v', G%axesCv1, Time, &
       'Inverse eddy time-scale, S*N, at v-points', 's^-1')
    CS%id_L2u = register_diag_field('ocean_model', 'L2u', G%axesCu1, Time, &
       'Length scale squared for mixing coefficient, at u-points', 'm^2')
    CS%id_L2v = register_diag_field('ocean_model', 'L2v', G%axesCv1, Time, &
       'Length scale squared for mixing coefficient, at v-points', 'm^2')
  endif

  if (CS%Resoln_scaled_Kh .or. Resoln_scaled_KhTh .or. Resoln_scaled_KhTr) then
    call wave_speed_init(Time, G, param_file, diag, CS%wave_speed_CSp)

    ! Allocate and initialize various arrays.
    allocate(CS%Res_fn_h(isd:ied,jsd:jed))       ; CS%Res_fn_h(:,:) = 0.0
    allocate(CS%Res_fn_q(IsdB:IedB,JsdB:JedB))   ; CS%Res_fn_q(:,:) = 0.0
    allocate(CS%cg1(isd:ied,jsd:jed))            ; CS%cg1(:,:) = 0.0
    allocate(CS%beta_dx2_h(isd:ied,jsd:jed))     ; CS%beta_dx2_h(:,:) = 0.0
    allocate(CS%beta_dx2_q(IsdB:IedB,JsdB:JedB)) ; CS%beta_dx2_q(:,:) = 0.0
    allocate(CS%f2_dx2_h(isd:ied,jsd:jed))       ; CS%f2_dx2_h(:,:) = 0.0
    allocate(CS%f2_dx2_q(IsdB:IedB,JsdB:JedB))   ; CS%f2_dx2_q(:,:) = 0.0
    allocate(CS%Rd_dx_h(isd:ied,jsd:jed))        ; CS%Rd_dx_h(:,:) = 0.0


    CS%id_Res_fn = register_diag_field('ocean_model', 'Res_fn', G%axesT1, Time, &
       'Resolution function for scaling diffusivities', 'Nondim')
    CS%id_Rd_dx = register_diag_field('ocean_model', 'Rd_dx', G%axesT1, Time, &
       'Ratio between deformation radius and grid spacing', 'Nondim')

    call get_param(param_file, mod, "KH_RES_SCALE_COEF", CS%Res_coef, &
                 "A coefficient that determines how Kh is scaled away if \n"//&
                 "RESOLN_SCALED_... is true, as \n"//&
                 "F = 1 / (1 + (KH_RES_SCALE_COEF*Rd/dx)^KH_RES_FN_POWER).", &
                 units="nondim", default=1.0)
    call get_param(param_file, mod, "KH_RES_FN_POWER", CS%Res_fn_power, &
                 "The power of dx/Ld in the resolution function.  Any \n"//&
                 "positive integer may be used, although even integers \n"//&
                 "are more efficient to calculate.  Setting this greater \n"//&
                 "than 100 results in a step-function being used.", &
                 units="nondim", default=2)

    ! Pre-calculate several static expressions for later use.
    do j=js-1,je+1 ; do i=is-1,ie+1
      CS%f2_dx2_h(i,j) = (G%dxT(i,j)**2 + G%dyT(i,j)**2) * &
          max(0.25 * ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                      (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2)), &
              absurdly_small_freq2)
      CS%beta_dx2_h(i,j) = (G%dxT(i,j)**2 + G%dyT(i,j)**2) * (sqrt(0.5 * &
          ( (((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
             ((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2) + &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
             ((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2) ) ))
    enddo ; enddo

    do J=js-1,Jeq ; do I=is-1,Ieq
      CS%f2_dx2_q(I,J) = (G%dxBu(i,j)**2 + G%dyBu(i,j)**2) * &
                         max(G%CoriolisBu(I,J)**2, absurdly_small_freq2)
      CS%beta_dx2_q(I,J) = (G%dxBu(i,j)**2 + G%dyBu(i,j)**2) * (sqrt(0.5 * &
          ( (((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
             ((G%CoriolisBu(I+1,J)-G%CoriolisBu(I,J)) * G%IdxCv(i+1,J))**2) + &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
             ((G%CoriolisBu(I,J+1)-G%CoriolisBu(I,J)) * G%IdyCu(I,j+1))**2) ) ))
    enddo ; enddo

  endif

end subroutine VarMix_init

end module MOM_lateral_mixing_coeffs
