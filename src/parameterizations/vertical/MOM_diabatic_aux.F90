module MOM_diabatic_aux
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
!*  By Robert Hallberg, April 1994 - July 2000                         *
!*     Alistair Adcroft, and Stephen Griffies                          *
!*                                                                     *
!*    This program contains the subroutine that, along with the        *
!*  subroutines that it calls, implements diapycnal mass and momentum  *
!*  fluxes and a bulk mixed layer.  The diapycnal diffusion can be     *
!*  used without the bulk mixed layer.                                 *
!*                                                                     *
!*    diabatic first determines the (diffusive) diapycnal mass fluxes  *
!*  based on the convergence of the buoyancy fluxes within each layer. *
!*  The dual-stream entrainment scheme of MacDougall and Dewar (JPO,   *
!*  1997) is used for combined diapycnal advection and diffusion,      *
!*  calculated implicitly and potentially with the Richardson number   *
!*  dependent mixing, as described by Hallberg (MWR, 2000). Diapycnal  *
!*  advection is fundamentally the residual of diapycnal diffusion,    *
!*  so the fully implicit upwind differencing scheme that is used is   *
!*  entirely appropriate.  The downward buoyancy flux in each layer    *
!*  is determined from an implicit calculation based on the previously *
!*  calculated flux of the layer above and an estimated flux in the    *
!*  layer below.  This flux is subject to the following conditions:    *
!*  (1) the flux in the top and bottom layers are set by the boundary  *
!*  conditions, and (2) no layer may be driven below an Angstrom thick-*
!*  ness.  If there is a bulk mixed layer, the buffer layer is treat-  *
!*  ed as a fixed density layer with vanishingly small diffusivity.    *
!*                                                                     *
!*    diabatic takes 5 arguments:  the two velocities (u and v), the   *
!*  thicknesses (h), a structure containing the forcing fields, and    *
!*  the length of time over which to act (dt).  The velocities and     *
!*  thickness are taken as inputs and modified within the subroutine.  *
!*  There is no limit on the time step.                                *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, T, S, buoy, ustar, ea, eb, etc.       *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,     only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_EOS,           only : calculate_density, calculate_TFreeze
use MOM_EOS,           only : calculate_specific_vol_derivs
use MOM_error_handler, only : MOM_error, FATAL, WARNING, callTree_showQuery
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_forcing_type,  only : forcing, extractFluxes1d, forcing_SinglePointPrint
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : vardesc
use MOM_shortwave_abs, only : absorbRemainingSW, optics_type
use MOM_variables,     only : thermo_var_ptrs, vertvisc_type! , accel_diag_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public diabatic_aux_init, diabatic_aux_end
public make_frazil, adjust_salt, insert_brine, differential_diffuse_T_S, triDiagTS
public find_uv_at_h, diagnoseMLDbyDensityDifference, applyBoundaryFluxesInOut

!> Control structure for diabatic_aux
type, public :: diabatic_aux_CS ; private
  logical :: do_rivermix = .false. !< Provide additional TKE to mix river runoff
                                   !! at the river mouths to "rivermix_depth" meters
  real    :: rivermix_depth = 0.0  !< The depth to which rivers are mixed if
                                   !! do_rivermix = T, in m.
  real, public    :: minimum_forcing_depth = 0.001 !< The smallest depth over which forcing is
                                   !! applied, in m.
  real, public    :: evap_CFL_limit = 0.8  !< The largest fraction of a layer that can be
                                   !! evaporated in one time-step (non-dim).

  logical :: reclaim_frazil  !<   If true, try to use any frazil heat deficit to
                             !! to cool the topmost layer down to the freezing
                             !! point.  The default is false.
  logical :: pressure_dependent_frazil  !< If true, use a pressure dependent
                             !! freezing temperature when making frazil.  The
                             !! default is false, which will be faster but is
                             !! inappropriate with ice-shelf cavities.
  logical :: use_river_heat_content !< If true, assumes that ice-ocean boundary
                             !! has provided a river heat content. Otherwise, runoff
                             !! is added with a temperature of the local SST.
  logical :: use_calving_heat_content !< If true, assumes that ice-ocean boundary
                             !! has provided a calving heat content. Otherwise, calving
                             !! is added with a temperature of the local SST.

  type(diag_ctrl), pointer :: diag !< Structure used to regulate timing of diagnostic output

  ! Diagnostic handles
  integer :: id_createdH       = -1
  integer :: id_brine_lay      = -1
  integer :: id_penSW_diag     = -1 !< Penetrative shortwave heating (flux convergence) diagnostic
  integer :: id_penSWflux_diag = -1 !< Penetrative shortwave flux diagnostic
  integer :: id_nonpenSW_diag  = -1 !< Non-penetrative shortwave heating diagnostic

  ! Optional diagnostic arrays
  real, allocatable, dimension(:,:)   :: createdH       !< The amount of volume added in order to avoid grounding (m/s)
  real, allocatable, dimension(:,:,:) :: penSW_diag     !< Heating in a layer from convergence of penetrative SW (W/m2)
  real, allocatable, dimension(:,:,:) :: penSWflux_diag !< Penetrative SW flux at base of grid layer (W/m2)
  real, allocatable, dimension(:,:)   :: nonpenSW_diag  !< Non-downwelling SW radiation (W/m2) at ocean surface

end type diabatic_aux_CS

integer :: id_clock_uv_at_h, id_clock_frazil

contains

subroutine make_frazil(h, tv, G, GV, CS, p_surf)
  type(ocean_grid_type),                 intent(in)    :: G
  type(verticalGrid_type),               intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
  type(thermo_var_ptrs),                 intent(inout) :: tv
  type(diabatic_aux_CS),                 intent(in)    :: CS
  real, dimension(SZI_(G),SZJ_(G)), intent(in), optional :: p_surf

!   Frazil formation keeps the temperature above the freezing point.
! This subroutine warms any water that is colder than the (currently
! surface) freezing point up to the freezing point and accumulates
! the required heat (in J m-2) in tv%frazil.
!   The expression, below, for the freezing point of sea water comes
! from Millero (1978) via Appendix A of Gill, 1982.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (in/out)  tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 diabatic_driver_init.
  real, dimension(SZI_(G)) :: &
    fraz_col, & ! The accumulated heat requirement due to frazil, in J.
    T_freeze, & ! The freezing potential temperature at the current salinity, C.
    ps          ! pressure
  real, dimension(SZI_(G),SZK_(G)) :: &
    pressure    ! The pressure at the middle of each layer in Pa.
  real :: hc    ! A layer's heat capacity in J m-2 K-1.
  logical :: T_fr_set  ! True if the freezing point has been calculated for a
                       ! row of points.
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call cpu_clock_begin(id_clock_frazil)

  if (.not.CS%pressure_dependent_frazil) then
    do k=1,nz ; do i=is,ie ; pressure(i,k) = 0.0 ; enddo ; enddo
  endif
!$OMP parallel do default(none) shared(is,ie,js,je,CS,G,GV,h,nz,tv,p_surf) &
!$OMP                           private(fraz_col,T_fr_set,T_freeze,hc,ps)  &
!$OMP                      firstprivate(pressure)
  do j=js,je
     ps(:) = 0.0
     if (PRESENT(p_surf)) then
       ps(:) = p_surf(:,j)
     endif

    do i=is,ie ; fraz_col(:) = 0.0 ; enddo

    if (CS%pressure_dependent_frazil) then
      do i=is,ie
        pressure(i,1) = ps(i) + (0.5*GV%H_to_Pa)*h(i,j,1)
      enddo
      do k=2,nz ; do i=is,ie
        pressure(i,k) = pressure(i,k-1) + &
          (0.5*GV%H_to_Pa) * (h(i,j,k) + h(i,j,k-1))
      enddo ; enddo
    endif

    if (CS%reclaim_frazil) then
      T_fr_set = .false.
      do i=is,ie ; if (tv%frazil(i,j) > 0.0) then
        if (.not.T_fr_set) then
          call calculate_TFreeze(tv%S(i:,j,1), pressure(i:,1), T_freeze(i:), &
                                 1, ie-i+1, tv%eqn_of_state)
          T_fr_set = .true.
        endif

        if (tv%T(i,j,1) > T_freeze(i)) then
    ! If frazil had previously been formed, but the surface temperature is now
    ! above freezing, cool the surface layer with the frazil heat deficit.
          hc = (tv%C_p*GV%H_to_kg_m2) * h(i,j,1)
          if (tv%frazil(i,j) - hc * (tv%T(i,j,1) - T_freeze(i)) <= 0.0) then
            tv%T(i,j,1) = tv%T(i,j,1) - tv%frazil(i,j)/hc
            tv%frazil(i,j) = 0.0
          else
            tv%frazil(i,j) = tv%frazil(i,j) - hc * (tv%T(i,j,1) - T_freeze(i))
            tv%T(i,j,1) = T_freeze(i)
          endif
        endif
      endif ; enddo
    endif

    do k=nz,1,-1
      T_fr_set = .false.
      do i=is,ie
        if ((G%mask2dT(i,j) > 0.0) .and. &
            ((tv%T(i,j,k) < 0.0) .or. (fraz_col(i) > 0.0))) then
          if (.not.T_fr_set) then
            call calculate_TFreeze(tv%S(i:,j,k), pressure(i:,k), T_freeze(i:), &
                                   1, ie-i+1, tv%eqn_of_state)
            T_fr_set = .true.
          endif

          hc = (tv%C_p*GV%H_to_kg_m2) * h(i,j,k)
          if (h(i,j,k) <= 10.0*GV%Angstrom) then
            ! Very thin layers should not be cooled by the frazil flux.
            if (tv%T(i,j,k) < T_freeze(i)) then
              fraz_col(i) = fraz_col(i) + hc * (T_freeze(i) - tv%T(i,j,k))
              tv%T(i,j,k) = T_freeze(i)
            endif
          else
            if (fraz_col(i) + hc * (T_freeze(i) - tv%T(i,j,k)) <= 0.0) then
              tv%T(i,j,k) = tv%T(i,j,k) - fraz_col(i)/hc
              fraz_col(i) = 0.0
            else
              fraz_col(i) = fraz_col(i) + hc * (T_freeze(i) - tv%T(i,j,k))
              tv%T(i,j,k) = T_freeze(i)
            endif
          endif
        endif
      enddo
    enddo
    do i=is,ie
      tv%frazil(i,j) = tv%frazil(i,j) + fraz_col(i)
    enddo
  enddo
  call cpu_clock_end(id_clock_frazil)

end subroutine make_frazil

subroutine differential_diffuse_T_S(h, tv, visc, dt, G, GV)
  type(ocean_grid_type),                 intent(in)    :: G
  type(verticalGrid_type),               intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
  type(thermo_var_ptrs),                 intent(inout) :: tv
  type(vertvisc_type),                   intent(in)    :: visc
  real,                                  intent(in)    :: dt

! This subroutine applies double diffusion to T & S, assuming no diapycal mass
! fluxes, using a simple triadiagonal solver.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      visc - A structure containing vertical viscosities, bottom boundary
!                   layer properies, and related fields.
!  (in)      dt - Time increment, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.

  real, dimension(SZI_(G)) :: &
    b1_T, b1_S, &  !  Variables used by the tridiagonal solvers of T & S, in H.
    d1_T, d1_S     !  Variables used by the tridiagonal solvers, nondim.
  real, dimension(SZI_(G),SZK_(G)) :: &
    c1_T, c1_S     !  Variables used by the tridiagonal solvers, in m or kg m-2.
  real, dimension(SZI_(G),SZK_(G)+1) :: &
    mix_T, mix_S   !  Mixing distances in both directions across each
                   !  interface, in m or kg m-2.
  real :: h_tr         ! h_tr is h at tracer points with a tiny thickness
                       ! added to ensure positive definiteness, in m or kg m-2.
  real :: h_neglect    ! A thickness that is so small it is usually lost
                       ! in roundoff and can be neglected, in m or kg m-2.
  real :: I_h_int      ! The inverse of the thickness associated with an
                       ! interface, in m-1 or m2 kg-1.
  real :: b_denom_T    ! The first term in the denominators for the expressions
  real :: b_denom_S    ! for b1_T and b1_S, both in m or kg m-2.

  integer :: i, j, k, is, ie, js, je, nz
  real, pointer :: T(:,:,:), S(:,:,:), Kd_T(:,:,:), Kd_S(:,:,:)
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  h_neglect = GV%H_subroundoff

  if (.not.associated(tv%T)) call MOM_error(FATAL, &
      "differential_diffuse_T_S: Called with an unassociated tv%T")
  if (.not.associated(tv%S)) call MOM_error(FATAL, &
      "differential_diffuse_T_S: Called with an unassociated tv%S")
  if (.not.associated(visc%Kd_extra_T)) call MOM_error(FATAL, &
      "differential_diffuse_T_S: Called with an unassociated visc%Kd_extra_T")
  if (.not.associated(visc%Kd_extra_S)) call MOM_error(FATAL, &
      "differential_diffuse_T_S: Called with an unassociated visc%Kd_extra_S")

  T => tv%T ; S => tv%S
  Kd_T => visc%Kd_extra_T ; Kd_S => visc%Kd_extra_S
!$OMP parallel do default(none) shared(is,ie,js,je,h,h_neglect,dt,Kd_T,Kd_S,G,GV,T,S,nz) &
!$OMP                          private(I_h_int,mix_T,mix_S,h_tr,b1_T,b1_S, &
!$OMP                                  d1_T,d1_S,c1_T,c1_S,b_denom_T,b_denom_S)
  do j=js,je
    do i=is,ie
      I_h_int = 1.0 / (0.5 * (h(i,j,1) + h(i,j,2)) + h_neglect)
      mix_T(i,2) = ((dt * Kd_T(i,j,2)) * GV%m_to_H**2) * I_h_int
      mix_S(i,2) = ((dt * Kd_S(i,j,2)) * GV%m_to_H**2) * I_h_int

      h_tr = h(i,j,1) + h_neglect
      b1_T(i) = 1.0 / (h_tr + mix_T(i,2))
      b1_S(i) = 1.0 / (h_tr + mix_S(i,2))
      d1_T(i) = h_tr * b1_T(i)
      d1_S(i) = h_tr * b1_S(i)
      T(i,j,1) = (b1_T(i)*h_tr)*T(i,j,1)
      S(i,j,1) = (b1_S(i)*h_tr)*S(i,j,1)
    enddo
    do k=2,nz-1 ; do i=is,ie
      ! Calculate the mixing across the interface below this layer.
      I_h_int = 1.0 / (0.5 * (h(i,j,k) + h(i,j,k+1)) + h_neglect)
      mix_T(i,K+1) = ((dt * Kd_T(i,j,K+1)) * GV%m_to_H**2) * I_h_int
      mix_S(i,K+1) = ((dt * Kd_S(i,j,K+1)) * GV%m_to_H**2) * I_h_int

      c1_T(i,k) = mix_T(i,K) * b1_T(i)
      c1_S(i,k) = mix_S(i,K) * b1_S(i)

      h_tr = h(i,j,k) + h_neglect
      b_denom_T = h_tr + d1_T(i)*mix_T(i,K)
      b_denom_S = h_tr + d1_S(i)*mix_S(i,K)
      b1_T(i) = 1.0 / (b_denom_T + mix_T(i,K+1))
      b1_S(i) = 1.0 / (b_denom_S + mix_S(i,K+1))
      d1_T(i) = b_denom_T * b1_T(i)
      d1_S(i) = b_denom_S * b1_S(i)

      T(i,j,k) = b1_T(i) * (h_tr*T(i,j,k) + mix_T(i,K)*T(i,j,k-1))
      S(i,j,k) = b1_S(i) * (h_tr*S(i,j,k) + mix_S(i,K)*S(i,j,k-1))
    enddo ; enddo
    do i=is,ie
      c1_T(i,nz) = mix_T(i,nz) * b1_T(i)
      c1_S(i,nz) = mix_S(i,nz) * b1_S(i)

      h_tr = h(i,j,nz) + h_neglect
      b1_T(i) = 1.0 / (h_tr + d1_T(i)*mix_T(i,nz))
      b1_S(i) = 1.0 / (h_tr + d1_S(i)*mix_S(i,nz))

      T(i,j,nz) = b1_T(i) * (h_tr*T(i,j,nz) + mix_T(i,nz)*T(i,j,nz-1))
      S(i,j,nz) = b1_S(i) * (h_tr*S(i,j,nz) + mix_S(i,nz)*S(i,j,nz-1))
    enddo
    do k=nz-1,1,-1 ; do i=is,ie
      T(i,j,k) = T(i,j,k) + c1_T(i,k+1)*T(i,j,k+1)
      S(i,j,k) = S(i,j,k) + c1_S(i,k+1)*S(i,j,k+1)
    enddo ; enddo
  enddo

end subroutine differential_diffuse_T_S

subroutine adjust_salt(h, tv, G, GV, CS)
  type(ocean_grid_type),                 intent(in)    :: G
  type(verticalGrid_type),               intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
  type(thermo_var_ptrs),                 intent(inout) :: tv
  type(diabatic_aux_CS),                 intent(in)    :: CS

!  Keep salinity from falling below a small but positive threshold
!  This occurs when the ice model attempts to extract more salt then
!  is actually available to it from the ocean.

! Arguments: h - Layer thickness, in m.
!  (in/out)  tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 diabatic_driver_init.
  real :: salt_add_col(SZI_(G),SZJ_(G)) ! The accumulated salt requirement
  real :: S_min      ! The minimum salinity
  real :: mc         ! A layer's mass kg  m-2 .
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

!  call cpu_clock_begin(id_clock_adjust_salt)

  S_min = 0.01

  salt_add_col(:,:) = 0.0

!$OMP parallel do default(none) shared(is,ie,js,je,nz,G,GV,tv,h,salt_add_col, S_min) &
!$OMP                          private(mc)
  do j=js,je
    do k=nz,1,-1 ; do i=is,ie
      if ((G%mask2dT(i,j) > 0.0) .and. &
           ((tv%S(i,j,k) < S_min) .or. (salt_add_col(i,j) > 0.0))) then
        mc = GV%H_to_kg_m2 * h(i,j,k)
        if (h(i,j,k) <= 10.0*GV%Angstrom) then
          ! Very thin layers should not be adjusted by the salt flux
          if (tv%S(i,j,k) < S_min) then
            salt_add_col(i,j) = salt_add_col(i,j) +  mc * (S_min - tv%S(i,j,k))
            tv%S(i,j,k) = S_min
          endif
        else
          if (salt_add_col(i,j) + mc * (S_min - tv%S(i,j,k)) <= 0.0) then
            tv%S(i,j,k) = tv%S(i,j,k) - salt_add_col(i,j)/mc
            salt_add_col(i,j) = 0.0
          else
            salt_add_col(i,j) = salt_add_col(i,j) + mc * (S_min - tv%S(i,j,k))
            tv%S(i,j,k) = S_min
          endif
        endif
      endif
    enddo ; enddo
    do i=is,ie
      tv%salt_deficit(i,j) = tv%salt_deficit(i,j) + salt_add_col(i,j)
    enddo
  enddo
!  call cpu_clock_end(id_clock_adjust_salt)

end subroutine adjust_salt

subroutine insert_brine(h, tv, G, GV, fluxes, nkmb, CS, dt, id_brine_lay)
  type(ocean_grid_type),                 intent(in)    :: G
  type(verticalGrid_type),               intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
  type(thermo_var_ptrs),                 intent(inout) :: tv
  type(forcing),                         intent(in)    :: fluxes
  integer,                               intent(in)    :: nkmb
  type(diabatic_aux_CS),                 intent(in)    :: CS
  real,                                  intent(in)    :: dt
  integer,                               intent(in)    :: id_brine_lay

! Insert salt from brine rejection into the first layer below
! the mixed layer which both contains mass and in which the
! change in layer density remains stable after the addition
! of salt via brine rejection.

! Arguments: h - Layer thickness, in m.
!  (in/out)  tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      fluxes = A structure containing pointers to any possible
!                     forcing fields; unused fields have NULL ptrs.
!  (in)      nkmb - The number of layers in the mixed and buffer layers.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 diabatic_driver_init.

  real :: salt(SZI_(G)) ! The amount of salt rejected from
                        !  sea ice. [grams]
  real :: dzbr(SZI_(G)) ! cumulative depth over which brine is distributed
  real :: inject_layer(SZI_(G),SZJ_(G)) ! diagnostic

  real :: p_ref_cv(SZI_(G))
  real :: T(SZI_(G),SZK_(G))
  real :: S(SZI_(G),SZK_(G))
  real :: h_2d(SZI_(G),SZK_(G))
  real :: Rcv(SZI_(G),SZK_(G))
  real :: mc  ! A layer's mass in kg m-2 .
  real :: s_new,R_new,t0,scale, cdz
  integer :: i, j, k, is, ie, js, je, nz, ks

  real, parameter :: brine_dz = 1.0  ! minumum thickness over which to distribute brine
  real, parameter :: s_max = 45.0    ! salinity bound

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.ASSOCIATED(fluxes%salt_flux)) return

  p_ref_cv(:)  = tv%P_ref

  inject_layer = nz

  do j=js,je

    salt(:)=0.0 ; dzbr(:)=0.0

    do i=is,ie ; if (G%mask2dT(i,j) > 0.) then
      salt(i) = dt * (1000. * fluxes%salt_flux(i,j))
    endif ; enddo

    do k=1,nz
      do i=is,ie
        T(i,k)=tv%T(i,j,k); S(i,k)=tv%S(i,j,k); h_2d(i,k)=h(i,j,k)
      enddo

      call calculate_density(T(:,k), S(:,k), p_ref_cv, Rcv(:,k), is, &
                             ie-is+1, tv%eqn_of_state)
    enddo

    ! First, try to find an interior layer where inserting all the salt
    ! will not cause the layer to become statically unstable.
    ! Bias towards deeper layers.

    do k=nkmb+1,nz-1 ; do i=is,ie
      if ((G%mask2dT(i,j) > 0.0) .and. dzbr(i) < brine_dz .and. salt(i) > 0.) then
        mc = GV%H_to_kg_m2 * h_2d(i,k)
        s_new = S(i,k) + salt(i)/mc
        t0 = T(i,k)
        call calculate_density(t0,s_new,tv%P_Ref,R_new,tv%eqn_of_state)
        if (R_new < 0.5*(Rcv(i,k)+Rcv(i,k+1)) .and. s_new<s_max) then
          dzbr(i)=dzbr(i)+h_2d(i,k)
          inject_layer(i,j) = min(inject_layer(i,j),real(k))
        endif
      endif
    enddo ; enddo

    ! Then try to insert into buffer layers if they exist
    do k=nkmb,GV%nkml+1,-1 ; do i=is,ie
      if ((G%mask2dT(i,j) > 0.0) .and. dzbr(i) < brine_dz .and. salt(i) > 0.) then
        mc = GV%H_to_kg_m2 * h_2d(i,k)
        dzbr(i)=dzbr(i)+h_2d(i,k)
        inject_layer(i,j) = min(inject_layer(i,j),real(k))
      endif
    enddo ; enddo

    ! finally if unable to find a layer to insert, then place in mixed layer

    do k=1,GV%nkml ; do i=is,ie
      if ((G%mask2dT(i,j) > 0.0) .and. dzbr(i) < brine_dz .and. salt(i) > 0.) then
        mc = GV%H_to_kg_m2 * h_2d(i,k)
        dzbr(i)=dzbr(i)+h_2d(i,k)
        inject_layer(i,j) = min(inject_layer(i,j),real(k))
      endif
    enddo ; enddo


    do i=is,ie
      if ((G%mask2dT(i,j) > 0.0) .and. salt(i) > 0.) then
    !   if (dzbr(i)< brine_dz) call MOM_error(FATAL,"insert_brine: failed")
        ks=inject_layer(i,j)
        cdz=0.0
        do k=ks,nz
          mc = GV%H_to_kg_m2 * h_2d(i,k)
          scale = h_2d(i,k)/dzbr(i)
          cdz=cdz+h_2d(i,k)
          if (cdz > 1.0) exit
          tv%S(i,j,k) = tv%S(i,j,k) + scale*salt(i)/mc
        enddo
      endif
    enddo

  enddo

  if (CS%id_brine_lay > 0) call post_data(CS%id_brine_lay,inject_layer,CS%diag)

end subroutine insert_brine

subroutine triDiagTS(G, GV, is, ie, js, je, hold, ea, eb, T, S)
! Simple tri-diagnonal solver for T and S
! "Simple" means it only uses arrays hold, ea and eb
  ! Arguments
  type(ocean_grid_type),                    intent(in)    :: G
  type(verticalGrid_type),                  intent(in)    :: GV
  integer,                                  intent(in)    :: is, ie, js, je
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: hold, ea, eb
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: T, S
  ! Local variables
  real :: b1(SZIB_(G)), d1(SZIB_(G)) ! b1, c1, and d1 are variables used by the
  real :: c1(SZIB_(G),SZK_(G))       ! tridiagonal solver.
  real :: h_tr, b_denom_1
  integer :: i, j, k
!$OMP parallel do default(none) shared(is,ie,js,je,G,GV,hold,eb,T,S,ea) &
!$OMP                          private(h_tr,b1,d1,c1,b_denom_1)
  do j=js,je
    do i=is,ie
      h_tr = hold(i,j,1) + GV%H_subroundoff
      b1(i) = 1.0 / (h_tr + eb(i,j,1))
      d1(i) = h_tr * b1(i)
      T(i,j,1) = (b1(i)*h_tr)*T(i,j,1)
      S(i,j,1) = (b1(i)*h_tr)*S(i,j,1)
    enddo
    do k=2,G%ke ; do i=is,ie
      c1(i,k) = eb(i,j,k-1) * b1(i)
      h_tr = hold(i,j,k) + GV%H_subroundoff
      b_denom_1 = h_tr + d1(i)*ea(i,j,k)
      b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
      d1(i) = b_denom_1 * b1(i)
      T(i,j,k) = b1(i) * (h_tr*T(i,j,k) + ea(i,j,k)*T(i,j,k-1))
      S(i,j,k) = b1(i) * (h_tr*S(i,j,k) + ea(i,j,k)*S(i,j,k-1))
    enddo ; enddo
    do k=G%ke-1,1,-1 ; do i=is,ie
      T(i,j,k) = T(i,j,k) + c1(i,k+1)*T(i,j,k+1)
      S(i,j,k) = S(i,j,k) + c1(i,k+1)*S(i,j,k+1)
    enddo ; enddo
  enddo
end subroutine triDiagTS


subroutine find_uv_at_h(u, v, h, u_h, v_h, G, GV, ea, eb)
  type(ocean_grid_type),                     intent(in)  :: G
  type(verticalGrid_type),                   intent(in)  :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)  :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)  :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)  :: h
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(out) :: u_h, v_h
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in), optional  :: ea, eb
!   This subroutine calculates u_h and v_h (velocities at thickness
! points), optionally using the entrainments (in m) passed in as arguments.

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (out)     u_h - The zonal velocity at thickness points after
!                  entrainment, in m s-1.
!  (out)     v_h - The meridional velocity at thickness points after
!                  entrainment, in m s-1.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in, opt) ea - The amount of fluid entrained from the layer above within
!                 this time step, in units of m or kg m-2.  Omitting ea is the
!                 same as setting it to 0.
!  (in, opt) eb - The amount of fluid entrained from the layer below within
!                 this time step, in units of m or kg m-2.  Omitting eb is the
!                 same as setting it to 0.  ea and eb must either be both
!                 present or both absent.

  real :: b_denom_1    ! The first term in the denominator of b1 in m or kg m-2.
  real :: h_neglect    ! A thickness that is so small it is usually lost
                       ! in roundoff and can be neglected, in m or kg m-2.
  real :: b1(SZI_(G)), d1(SZI_(G)), c1(SZI_(G),SZK_(G))
  real :: a_n(SZI_(G)), a_s(SZI_(G))  ! Fractional weights of the neighboring
  real :: a_e(SZI_(G)), a_w(SZI_(G))  ! velocity points, ~1/2 in the open
                                      ! ocean, nondimensional.
  real :: s, Idenom
  logical :: mix_vertically
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  call cpu_clock_begin(id_clock_uv_at_h)
  h_neglect = GV%H_subroundoff

  mix_vertically = present(ea)
  if (present(ea) .neqv. present(eb)) call MOM_error(FATAL, &
      "find_uv_at_h: Either both ea and eb or neither one must be present "// &
      "in call to find_uv_at_h.")
!$OMP parallel do default(none) shared(is,ie,js,je,G,GV,mix_vertically,h,h_neglect, &
!$OMP                                  eb,u_h,u,v_h,v,nz,ea)                     &
!$OMP                          private(s,Idenom,a_w,a_e,a_s,a_n,b_denom_1,b1,d1,c1)
  do j=js,je
    do i=is,ie
      s = G%areaCu(I-1,j)+G%areaCu(I,j)
      if (s>0.0) then
        Idenom = sqrt(0.5*G%IareaT(i,j)/s)
        a_w(i) = G%areaCu(I-1,j)*Idenom
        a_e(i) = G%areaCu(I,j)*Idenom
      else
        a_w(i) = 0.0 ; a_e(i) = 0.0
      endif

      s = G%areaCv(i,J-1)+G%areaCv(i,J)
      if (s>0.0) then
        Idenom = sqrt(0.5*G%IareaT(i,j)/s)
        a_s(i) = G%areaCv(i,J-1)*Idenom
        a_n(i) = G%areaCv(i,J)*Idenom
      else
        a_s(i) = 0.0 ; a_n(i) = 0.0
      endif
    enddo

    if (mix_vertically) then
      do i=is,ie
        b_denom_1 = h(i,j,1) + h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,1))
        d1(i) = b_denom_1 * b1(i)
        u_h(i,j,1) = (h(i,j,1)*b1(i)) * (a_e(i)*u(I,j,1) + a_w(i)*u(I-1,j,1))
        v_h(i,j,1) = (h(i,j,1)*b1(i)) * (a_n(i)*v(i,J,1) + a_s(i)*v(i,J-1,1))
      enddo
      do k=2,nz ; do i=is,ie
        c1(i,k) = eb(i,j,k-1) * b1(i)
        b_denom_1 = h(i,j,k) + d1(i)*ea(i,j,k) + h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
        d1(i) = b_denom_1 * b1(i)
        u_h(i,j,k) = (h(i,j,k) * (a_e(i)*u(I,j,k) + a_w(i)*u(I-1,j,k)) + &
                      ea(i,j,k)*u_h(i,j,k-1))*b1(i)
        v_h(i,j,k) = (h(i,j,k) * (a_n(i)*v(i,J,k) + a_s(i)*v(i,J-1,k)) + &
                      ea(i,j,k)*v_h(i,j,k-1))*b1(i)
      enddo ; enddo
      do k=nz-1,1,-1 ; do i=is,ie
        u_h(i,j,k) = u_h(i,j,k) + c1(i,k+1)*u_h(i,j,k+1)
        v_h(i,j,k) = v_h(i,j,k) + c1(i,k+1)*v_h(i,j,k+1)
      enddo ; enddo
    else
      do k=1,nz ; do i=is,ie
        u_h(i,j,k) = a_e(i)*u(I,j,k) + a_w(i)*u(I-1,j,k)
        v_h(i,j,k) = a_n(i)*v(i,J,k) + a_s(i)*v(i,J-1,k)
      enddo ; enddo
    endif
  enddo

  call cpu_clock_end(id_clock_uv_at_h)
end subroutine find_uv_at_h


!> Diagnose a mixed layer depth (MLD) determined by a given density difference with the surface.
!> This routine is appropriate in MOM_diabatic_driver due to its position within the time stepping.
subroutine diagnoseMLDbyDensityDifference(id_MLD, h, tv, densityDiff, G, GV, diagPtr, id_N2subML, id_MLDsq)
  type(ocean_grid_type),                 intent(in) :: G           !< Grid type
  type(verticalGrid_type),               intent(in) :: GV          !< ocean vertical grid structure
  integer,                               intent(in) :: id_MLD      !< Handle (ID) of MLD diagnostic
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h        !< Layer thickness
  type(thermo_var_ptrs),                 intent(in) :: tv          !< Thermodynamics type
  real,                                  intent(in) :: densityDiff !< Density difference to determine MLD (kg/m3)
  type(diag_ctrl),                       pointer    :: diagPtr     !< Diagnostics structure
  integer,                     optional, intent(in) :: id_N2subML  !< Optional handle (ID) of subML stratification
  integer,                     optional, intent(in) :: id_MLDsq    !< Optional handle (ID) of squared MLD

  ! Local variables
  real, dimension(SZI_(G))          :: rhoSurf, deltaRhoAtKm1, deltaRhoAtK, dK, dKm1, pRef_MLD
  real, dimension(SZI_(G))          :: rhoAtK, rho1, d1, pRef_N2 ! Used for N2
  real, dimension(SZI_(G), SZJ_(G)) :: MLD ! Diagnosed mixed layer depth
  real, dimension(SZI_(G), SZJ_(G)) :: subMLN2 ! Diagnosed stratification below ML
  real, dimension(SZI_(G), SZJ_(G)) :: MLD2 ! Diagnosed MLD^2
  real, parameter                   :: dz_subML = 50. ! Depth below ML over which to diagnose stratification (m)
  integer :: i, j, is, ie, js, je, k, nz, id_N2, id_SQ
  real :: aFac, ddRho

  id_N2 = -1
  if (PRESENT(id_N2subML)) id_N2 = id_N2subML

  id_SQ = -1
  if (PRESENT(id_N2subML)) id_SQ = id_MLDsq

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  pRef_MLD(:) = 0. ; pRef_N2(:) = 0.
  do j = js, je
    dK(:) = 0.5 * h(:,j,1) * GV%H_to_m ! Depth of center of surface layer
    call calculate_density(tv%T(:,j,1), tv%S(:,j,1), pRef_MLD, rhoSurf, is, ie-is+1, tv%eqn_of_state)
    deltaRhoAtK(:) = 0.
    MLD(:,j) = 0.
    if (id_N2>0) then
      subMLN2(:,j) = 0.
      rho1(:) = 0.
      d1(:) = 0.
      pRef_N2(:) = GV%g_Earth * GV%Rho0 * h(:,j,1) * GV%H_to_m ! Boussinesq approximation!!!! ?????
    endif
    do k = 2, nz
      dKm1(:) = dK(:) ! Depth of center of layer K-1
      dK(:) = dK(:) + 0.5 * ( h(:,j,k) + h(:,j,k-1) ) * GV%H_to_m ! Depth of center of layer K

      ! Stratification, N2, immediately below the mixed layer, averaged over at least 50 m.
      if (id_N2>0) then
        pRef_N2(:) = pRef_N2(:) + GV%g_Earth * GV%Rho0 * h(:,j,k) * GV%H_to_m ! Boussinesq approximation!!!! ?????
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pRef_N2, rhoAtK, is, ie-is+1, tv%eqn_of_state)
        do i = is, ie
          if (MLD(i,j)>0. .and. subMLN2(i,j)==0.) then ! This block is below the mixed layer
            if (d1(i)==0.) then ! Record the density, depth and pressure, immediately below the ML
              rho1(i) = rhoAtK(i)
              d1(i) = dK(i)
              ! Use pressure at the bottom of the upper layer used in calculating d/dz rho
              pRef_N2(i) = pRef_N2(i) + GV%g_Earth * GV%Rho0 * h(i,j,k) * GV%H_to_m ! Boussinesq approximation!!!! ?????
            endif
            if (d1(i)>0. .and. dK(i)-d1(i)>=dz_subML) then
              subMLN2(i,j) = GV%g_Earth/ GV%Rho0 * (rho1(i)-rhoAtK(i)) / (d1(i) - dK(i))
            endif
          endif
        enddo ! i-loop
      endif ! id_N2>0

      ! Mixed-layer depth, using sigma-0 (surface reference pressure)
      deltaRhoAtKm1(:) = deltaRhoAtK(:) ! Store value from previous iteration of K
      call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pRef_MLD, deltaRhoAtK, is, ie-is+1, tv%eqn_of_state)
      do i = is, ie
        deltaRhoAtK(i) = deltaRhoAtK(i) - rhoSurf(i) ! Density difference between layer K and surface
        ddRho = deltaRhoAtK(i) - deltaRhoAtKm1(i)
        if ((MLD(i,j)==0.) .and. (ddRho>0.) .and. &
            (deltaRhoAtKm1(i)<densityDiff) .and. (deltaRhoAtK(i)>=densityDiff)) then
          aFac = ( densityDiff - deltaRhoAtKm1(i) ) / ddRho
          MLD(i,j) = dK(i) * aFac + dKm1(i) * (1. - aFac)
        endif
      enddo ! i-loop
      if (id_SQ > 0) MLD2(is:ie,j) = MLD(is:ie,j)**2
    enddo ! k-loop
    do i = is, ie
      if ((MLD(i,j)==0.) .and. (deltaRhoAtK(i)<densityDiff)) MLD(i,j) = dK(i) ! Assume mixing to the bottom
   !  if (id_N2>0 .and. subMLN2(i,j)==0. .and. d1(i)>0. .and. dK(i)-d1(i)>0.) then
   !    ! Use what ever stratification we can, measured over what ever distance is available
   !    subMLN2(i,j) = GV%g_Earth/ GV%Rho0 * (rho1(i)-rhoAtK(i)) / (d1(i) - dK(i))
   !  endif
    enddo
  enddo ! j-loop

  if (id_MLD > 0) call post_data(id_MLD, MLD, diagPtr)
  if (id_N2 > 0)  call post_data(id_N2, subMLN2 , diagPtr)
  if (id_SQ > 0)  call post_data(id_SQ, MLD2, diagPtr)

end subroutine diagnoseMLDbyDensityDifference

!> Update the thickness, temperature, and salinity due to thermodynamic
!! boundary forcing (contained in fluxes type) applied to h, tv%T and tv%S,
!! and calculate the TKE implications of this heating.
subroutine applyBoundaryFluxesInOut(CS, G, GV, dt, fluxes, optics, h, tv, &
                                    aggregate_FW_forcing, cTKE, dSV_dT, dSV_dS)
  type(diabatic_aux_CS),                 pointer       :: CS !< Control structure for diabatic_aux
  type(ocean_grid_type),                 intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),               intent(in)    :: GV        !< ocean vertical grid structure
  real,                                  intent(in)    :: dt !< Time-step over which forcing is applied (s)
  type(forcing),                         intent(inout) :: fluxes !< Surface fluxes container
  type(optics_type),                     pointer       :: optics !< Optical properties container
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: h  !< Layer thickness in H units
  type(thermo_var_ptrs),                 intent(inout) :: tv !< Thermodynamics container
  !> If False, treat in/out fluxes separately.
  logical,                               intent(in)    :: aggregate_FW_forcing
  !> Turbulent kinetic energy requirement to mix forcing through each layer, in W m-2
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), optional, intent(out) :: cTKE
  !> Partial derivative of specific volume with potential temperature, in m3 kg-1 K-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), optional, intent(out) :: dSV_dT
  !> Partial derivative of specific a volume with potential salinity, in m3 kg-1 / (g kg-1).
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), optional, intent(out) :: dSV_dS

  ! Local variables
  integer, parameter :: maxGroundings = 5
  integer :: numberOfGroundings, iGround(maxGroundings), jGround(maxGroundings)
  real :: H_limit_fluxes, IforcingDepthScale, Idt
  real :: dThickness, dTemp, dSalt
  real :: fractionOfForcing, hOld, Ithickness
  real :: RivermixConst  ! A constant used in implementing river mixing, in Pa s.
  real, dimension(SZI_(G)) :: &
    d_pres,       &  ! pressure change across a layer (Pa)
    p_lay,        &  ! average pressure in a layer (Pa)
    pres,         &  ! pressure at an interface (Pa)
    netMassInOut, &  ! surface water fluxes (H units) over time step
    netMassIn,    &  ! mass entering ocean surface (H units) over a time step
    netMassOut,   &  ! mass leaving ocean surface (H units) over a time step
    netHeat,      &  ! heat (degC * H) via surface fluxes, excluding
                     ! Pen_SW_bnd and netMassOut
    netSalt,      &  ! surface salt flux ( g(salt)/m2 for non-Bouss and ppt*H for Bouss )
    nonpenSW         ! non-downwelling SW, which is absorbed at ocean surface
  real, dimension(SZI_(G), SZK_(G))                     :: h2d, T2d
  real, dimension(SZI_(G), SZK_(G))                     :: pen_TKE_2d, dSV_dT_2d
  real, dimension(max(optics%nbands,1),SZI_(G))         :: Pen_SW_bnd
  real, dimension(max(optics%nbands,1),SZI_(G),SZK_(G)) :: opacityBand
  real                                                  :: hGrounding(maxGroundings)
  real    :: Temp_in, Salin_in
  real    :: I_G_Earth, g_Hconv2
  logical :: calculate_energetics
  integer :: i, j, is, ie, js, je, k, nz, n, nsw
  character(len=45) :: mesg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! Only apply forcing if fluxes%sw is associated.
  if (.not.ASSOCIATED(fluxes%sw)) return

#define _OLD_ALG_
  nsw = optics%nbands
  Idt = 1.0/dt

  calculate_energetics = (present(cTKE) .and. present(dSV_dT) .and. present(dSV_dS))
  I_G_Earth = 1.0 / GV%g_Earth
  g_Hconv2 = GV%g_Earth * GV%H_to_kg_m2**2

  if (present(cTKE)) cTKE(:,:,:) = 0.0

  ! H_limit_fluxes is used by extractFluxes1d to scale down fluxes if the total
  ! depth of the ocean is vanishing. It does not (yet) handle a value of zero.
  ! To accommodate vanishing upper layers, we need to allow for an instantaneous
  ! distribution of forcing over some finite vertical extent. The bulk mixed layer
  ! code handles this issue properly.
  H_limit_fluxes = max(GV%Angstrom, 1.E-30*GV%m_to_H)

  ! diagnostic to see if need to create mass to avoid grounding
  if (CS%id_createdH>0) CS%createdH(:,:) = 0.
  numberOfGroundings = 0

!$OMP parallel do default(none) shared(is,ie,js,je,nz,h,tv,nsw,G,GV,optics,fluxes,dt,    &
!$OMP                                  H_limit_fluxes,IforcingDepthScale,                &
!$OMP                                  numberOfGroundings,iGround,jGround,nonPenSW,      &
!$OMP                                  hGrounding,CS,Idt,aggregate_FW_forcing,           &
!$OMP                                  calculate_energetics,dSV_dT,dSV_dS,cTKE,g_Hconv2) &
!$OMP                          private(opacityBand,h2d,T2d,netMassInOut,netMassOut,      &
!$OMP                                  netHeat,netSalt,Pen_SW_bnd,fractionOfForcing,     &
!$OMP                                  dThickness,dTemp,dSalt,hOld,Ithickness,           &
!$OMP                                  netMassIn,pres,d_pres,p_lay,dSV_dT_2d,            &
!$OMP                                  pen_TKE_2d,Temp_in,Salin_in,RivermixConst)

  ! Work in vertical slices for efficiency
  do j=js,je

    ! Copy state into 2D-slice arrays
    do k=1,nz ; do i=is,ie
      h2d(i,k) = h(i,j,k)
      T2d(i,k) = tv%T(i,j,k)
      do n=1,nsw
        opacityBand(n,i,k) = (1.0 / GV%m_to_H)*optics%opacity_band(n,i,j,k)
      enddo
    enddo ; enddo

    if (calculate_energetics) then
      ! The partial derivatives of specific volume with temperature and
      ! salinity need to be precalculated to avoid having heating of
      ! tiny layers give nonsensical values.
      do i=is,ie ; pres(i) = 0.0 ; enddo ! Add surface pressure?
      do k=1,nz
        do i=is,ie
          d_pres(i) = GV%g_Earth * GV%H_to_kg_m2 * h2d(i,k)
          p_lay(i) = pres(i) + 0.5*d_pres(i)
          pres(i) = pres(i) + d_pres(i)
        enddo
        call calculate_specific_vol_derivs(T2d(:,k), tv%S(:,j,k), p_lay(:),&
                 dSV_dT(:,j,k), dSV_dS(:,j,k), is, ie-is+1, tv%eqn_of_state)
        do i=is,ie ; dSV_dT_2d(i,k) = dSV_dT(i,j,k) ; enddo
!        do i=is,ie
!          dT_to_dPE(i,k) = I_G_Earth * d_pres(i) * p_lay(i) * dSV_dT(i,j,k)
!          dS_to_dPE(i,k) = I_G_Earth * d_pres(i) * p_lay(i) * dSV_dS(i,j,k)
!        enddo
      enddo
      pen_TKE_2d(:,:) = 0.0
    endif

    ! The surface forcing is contained in the fluxes type.
    ! We aggregate the thermodynamic forcing for a time step into the following:
    ! netMassInOut = surface water fluxes (H units) over time step
    !              = lprec + fprec + vprec + evap + lrunoff + frunoff
    !                note that lprec generally has sea ice melt/form included.
    ! netMassOut   = net mass leaving ocean surface (H units) over a time step.
    !                netMassOut < 0 means mass leaves ocean.
    ! netHeat      = heat (degC * H) via surface fluxes, excluding the part
    !                contained in Pen_SW_bnd; and excluding heat_content of netMassOut < 0.
    ! netSalt      = surface salt fluxes ( g(salt)/m2 for non-Bouss and ppt*H for Bouss )
    ! Pen_SW_bnd   = components to penetrative shortwave radiation split according to bands.
    !                This field provides that portion of SW from atmosphere that in fact
    !                enters to the ocean and participates in pentrative SW heating.
    ! nonpenSW     = non-downwelling SW flux, which is absorbed in ocean surface
    !                (in tandem w/ LW,SENS,LAT); saved only for diagnostic purposes.
    call extractFluxes1d(G, GV, fluxes, optics, nsw, j, dt,                               &
                  H_limit_fluxes, CS%use_river_heat_content, CS%use_calving_heat_content, &
                  h2d, T2d, netMassInOut, netMassOut, netHeat, netSalt,                   &
                  Pen_SW_bnd, tv, aggregate_FW_forcing, nonpenSW)

    ! ea is for passive tracers
    do i=is,ie
    !  ea(i,j,1) = netMassInOut(i)
      if (aggregate_FW_forcing) then
        netMassOut(i) = netMassInOut(i)
        netMassIn(i) = 0.
      else
        netMassIn(i) = netMassInOut(i) - netMassOut(i)
      endif
      if (G%mask2dT(i,j)>0.0) then
        fluxes%netMassOut(i,j) = netMassOut(i)
        fluxes%netMassIn(i,j) = netMassIn(i)
      else
        fluxes%netMassOut(i,j) = 0.0
        fluxes%netMassIn(i,j) = 0.0
      endif
    enddo

    ! Apply the surface boundary fluxes in three steps:
    ! A/ update mass, temp, and salinity due to all terms except mass leaving
    !    ocean (and corresponding outward heat content), and ignoring penetrative SW.
    ! B/ update mass, salt, temp from mass leaving ocean.
    ! C/ update temp due to penetrative SW
    do i=is,ie
      if (G%mask2dT(i,j)>0.) then

        ! A/ Update mass, temp, and salinity due to incoming mass flux.
        do k=1,1

          ! Change in state due to forcing
          dThickness = netMassIn(i) ! Since we are adding mass, we can use all of it
          dTemp = 0.
          dSalt = 0.

          ! Update the forcing by the part to be consumed within the present k-layer.
          ! If fractionOfForcing = 1, then updated netMassIn, netHeat, and netSalt vanish.
          netMassIn(i) = netMassIn(i) - dThickness
          ! This line accounts for the temperature of the mass exchange
          Temp_in = T2d(i,k)
          Salin_in = 0.0
          dTemp = dTemp + dThickness*Temp_in

          ! Diagnostics of heat content associated with mass fluxes
          if (ASSOCIATED(fluxes%heat_content_massin))                             &
            fluxes%heat_content_massin(i,j) = fluxes%heat_content_massin(i,j) +   &
                         T2d(i,k) * max(0.,dThickness) * GV%H_to_kg_m2 * fluxes%C_p * Idt
          if (ASSOCIATED(fluxes%heat_content_massout))                            &
            fluxes%heat_content_massout(i,j) = fluxes%heat_content_massout(i,j) + &
                         T2d(i,k) * min(0.,dThickness) * GV%H_to_kg_m2 * fluxes%C_p * Idt
          if (ASSOCIATED(tv%TempxPmE)) tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + &
                         T2d(i,k) * dThickness * GV%H_to_kg_m2

          ! Determine the energetics of river mixing before updating the state.
          if (calculate_energetics .and. associated(fluxes%lrunoff) .and. CS%do_rivermix) then
            ! Here we add an additional source of TKE to the mixed layer where river
            ! is present to simulate unresolved estuaries. The TKE input is diagnosed
            ! as follows:
            !   TKE_river[m3 s-3] = 0.5*rivermix_depth*g*(1/rho)*drho_ds*
            !                       River*(Samb - Sriver) = CS%mstar*U_star^3
            ! where River is in units of m s-1.
            ! Samb = Ambient salinity at the mouth of the estuary
            ! rivermix_depth =  The prescribed depth over which to mix river inflow
            ! drho_ds = The gradient of density wrt salt at the ambient surface salinity.
            ! Sriver = 0 (i.e. rivers are assumed to be pure freshwater)
            RivermixConst = -0.5*(CS%rivermix_depth*dt)*GV%m_to_H*GV%H_to_Pa

            cTKE(i,j,k) = cTKE(i,j,k) + max(0.0, RivermixConst*dSV_dS(i,j,1) * &
                  (fluxes%lrunoff(i,j) + fluxes%frunoff(i,j)) * tv%S(i,j,1))
          endif

          ! Update state
          hOld     = h2d(i,k)               ! Keep original thickness in hand
          h2d(i,k) = h2d(i,k) + dThickness  ! New thickness
          if (h2d(i,k) > 0.0) then
            if (calculate_energetics .and. (dThickness > 0.)) then
              ! Calculate the energy required to mix the newly added water over
              ! the topmost grid cell.  ###CHECK THE SIGNS!!!
              cTKE(i,j,k) = cTKE(i,j,k) + 0.5*g_Hconv2*(hOld*dThickness) * &
                 ((T2d(i,k) - Temp_in) * dSV_dT(i,j,k) + (tv%S(i,j,k) - Salin_in) * dSV_dS(i,j,k))
            endif
            Ithickness  = 1.0/h2d(i,k)      ! Inverse new thickness
            ! The "if"s below avoid changing T/S by roundoff unnecessarily
            if (dThickness /= 0. .or. dTemp /= 0.) T2d(i,k)    = (hOld*T2d(i,k)    + dTemp)*Ithickness
            if (dThickness /= 0. .or. dSalt /= 0.) tv%S(i,j,k) = (hOld*tv%S(i,j,k) + dSalt)*Ithickness

          endif

        enddo ! k=1,1

        ! B/ Update mass, salt, temp from mass leaving ocean and other fluxes of heat and salt.
        do k=1,nz

          ! Place forcing into this layer if this layer has nontrivial thickness.
          ! For layers thin relative to 1/IforcingDepthScale, then distribute
          ! forcing into deeper layers.
          IforcingDepthScale = 1. / max(GV%H_subroundoff, CS%minimum_forcing_depth*GV%m_to_H - netMassOut(i) )
          ! fractionOfForcing = 1.0, unless h2d is less than IforcingDepthScale.
          fractionOfForcing = min(1.0, h2d(i,k)*IforcingDepthScale)

          ! In the case with (-1)*netMassOut*fractionOfForcing greater than cfl*h, we
          ! limit the forcing applied to this cell, leaving the remaining forcing to
          ! be distributed downwards.
          if (-fractionOfForcing*netMassOut(i) > CS%evap_CFL_limit*h2d(i,k)) then
            fractionOfForcing = -CS%evap_CFL_limit*h2d(i,k)/netMassOut(i)
          endif

          ! Change in state due to forcing

          dThickness = max( fractionOfForcing*netMassOut(i), -h2d(i,k) )
          dTemp      = fractionOfForcing*netHeat(i)
          !   ### The 0.9999 here should become a run-time parameter?
          dSalt = max( fractionOfForcing*netSalt(i), -0.9999*h2d(i,k)*tv%S(i,j,k))

          ! Update the forcing by the part to be consumed within the present k-layer.
          ! If fractionOfForcing = 1, then new netMassOut vanishes.
          netMassOut(i) = netMassOut(i) - dThickness
          netHeat(i) = netHeat(i) - dTemp
          netSalt(i) = netSalt(i) - dSalt

          ! This line accounts for the temperature of the mass exchange
          dTemp = dTemp + dThickness*T2d(i,k)

          ! Diagnostics of heat content associated with mass fluxes
          if (ASSOCIATED(fluxes%heat_content_massin))                             &
            fluxes%heat_content_massin(i,j) = fluxes%heat_content_massin(i,j) +   &
                         tv%T(i,j,k) * max(0.,dThickness) * GV%H_to_kg_m2 * fluxes%C_p * Idt
          if (ASSOCIATED(fluxes%heat_content_massout))                            &
            fluxes%heat_content_massout(i,j) = fluxes%heat_content_massout(i,j) + &
                         tv%T(i,j,k) * min(0.,dThickness) * GV%H_to_kg_m2 * fluxes%C_p * Idt
          if (ASSOCIATED(tv%TempxPmE)) tv%TempxPmE(i,j) = tv%TempxPmE(i,j) + &
                         tv%T(i,j,k) * dThickness * GV%H_to_kg_m2
!NOTE tv%T should be T2d

          ! Update state by the appropriate increment.
          hOld     = h2d(i,k)               ! Keep original thickness in hand
          h2d(i,k) = h2d(i,k) + dThickness  ! New thickness

          if (h2d(i,k) > 0.) then
            if (calculate_energetics) then
              ! Calculate the energy required to mix the newly added water over
              ! the topmost grid cell, assuming that the fluxes of heat and salt
              ! and rejected brine are initially applied in vanishingly thin
              ! layers at the top of the layer before being mixed throughout
              ! the layer.  Note that dThickness is always <= 0. ###CHECK THE SIGNS!!!
              cTKE(i,j,k) = cTKE(i,j,k) - (0.5*h2d(i,k)*g_Hconv2) * &
                 ((dTemp - dthickness*T2d(i,k)) * dSV_dT(i,j,k) + &
                  (dSalt - dthickness*tv%S(i,j,k)) * dSV_dS(i,j,k))
            endif
            Ithickness  = 1.0/h2d(i,k) ! Inverse of new thickness
            T2d(i,k)    = (hOld*T2d(i,k) + dTemp)*Ithickness
            tv%S(i,j,k) = (hOld*tv%S(i,j,k) + dSalt)*Ithickness
          elseif (h2d(i,k) < 0.0) then ! h2d==0 is a special limit that needs no extra handling
            call forcing_SinglePointPrint(fluxes,G,i,j,'applyBoundaryFluxesInOut (h<0)')
            write(0,*) 'applyBoundaryFluxesInOut(): lon,lat=',G%geoLonT(i,j),G%geoLatT(i,j)
            write(0,*) 'applyBoundaryFluxesInOut(): netT,netS,netH=',netHeat(i),netSalt(i),netMassInOut(i)
            write(0,*) 'applyBoundaryFluxesInOut(): dT,dS,dH=',dTemp,dSalt,dThickness
            write(0,*) 'applyBoundaryFluxesInOut(): h(n),h(n+1),k=',hOld,h2d(i,k),k
            call MOM_error(FATAL, "MOM_diabatic_driver.F90, applyBoundaryFluxesInOut(): "//&
                           "Complete mass loss in column!")
          endif

        enddo ! k

      ! Check if trying to apply fluxes over land points
      elseif((abs(netHeat(i))+abs(netSalt(i))+abs(netMassIn(i))+abs(netMassOut(i)))>0.) then
        call forcing_SinglePointPrint(fluxes,G,i,j,'applyBoundaryFluxesInOut (land)')
        write(0,*) 'applyBoundaryFluxesInOut(): lon,lat=',G%geoLonT(i,j),G%geoLatT(i,j)
        write(0,*) 'applyBoundaryFluxesInOut(): netHeat,netSalt,netMassIn,netMassOut=',&
                   netHeat(i),netSalt(i),netMassIn(i),netMassOut(i)
        call MOM_error(FATAL, "MOM_diabatic_driver.F90, applyBoundaryFluxesInOut(): "//&
                              "Mass loss over land?")
      endif

      ! If anything remains after the k-loop, then we have grounded out, which is a problem.
      if (netMassIn(i)+netMassOut(i) /= 0.0) then
!$OMP critical
        numberOfGroundings = numberOfGroundings +1
        if (numberOfGroundings<=maxGroundings) then
          iGround(numberOfGroundings) = i ! Record i,j location of event for
          jGround(numberOfGroundings) = j ! warning message
          hGrounding(numberOfGroundings) = netMassIn(i)+netMassOut(i)
        endif
!$OMP end critical
        if (CS%id_createdH>0) CS%createdH(i,j) = CS%createdH(i,j) - (netMassIn(i)+netMassOut(i))/dt
      endif

    enddo ! i

    ! Step C/ in the application of fluxes
    ! Heat by the convergence of penetrating SW.
    ! SW penetrative heating uses the updated thickness from above.

    ! Save temperature before increment with SW heating
    ! and initialize CS%penSWflux_diag to zero.
    if(CS%id_penSW_diag > 0 .or. CS%id_penSWflux_diag > 0) then
      do k=1,nz ; do i=is,ie
        CS%penSW_diag(i,j,k)     = T2d(i,k)
        CS%penSWflux_diag(i,j,k) = 0.0
      enddo ; enddo
      k=nz+1 ; do i=is,ie
         CS%penSWflux_diag(i,j,k) = 0.0
      enddo
    endif

    if (calculate_energetics) then
      call absorbRemainingSW(G, GV, h2d, opacityBand, nsw, j, dt, H_limit_fluxes, &
                             .false., .true., T2d, Pen_SW_bnd, TKE=pen_TKE_2d, dSV_dT=dSV_dT_2d)
      k = 1 ! For setting break-points.
      do k=1,nz ; do i=is,ie
        cTKE(i,j,k) = cTKE(i,j,k) + pen_TKE_2d(i,k)
      enddo ; enddo
    else
      call absorbRemainingSW(G, GV, h2d, opacityBand, nsw, j, dt, H_limit_fluxes, &
                             .false., .true., T2d, Pen_SW_bnd)
    endif


    ! Step D/ copy updated thickness and temperature
    ! 2d slice now back into model state.
    do k=1,nz ; do i=is,ie
      h(i,j,k)    = h2d(i,k)
      tv%T(i,j,k) = T2d(i,k)
    enddo ; enddo

    ! Diagnose heating (W/m2) applied to a grid cell from SW penetration
    ! Also diagnose the penetrative SW heat flux at base of layer.
    if(CS%id_penSW_diag > 0 .or. CS%id_penSWflux_diag > 0) then

      ! convergence of SW into a layer
      do k=1,nz ; do i=is,ie
        CS%penSW_diag(i,j,k) = (T2d(i,k)-CS%penSW_diag(i,j,k))*h(i,j,k) * Idt * tv%C_p * GV%H_to_kg_m2
      enddo ; enddo

      ! Perform a cumulative sum upwards from bottom to
      ! diagnose penetrative SW flux at base of tracer cell.
      ! CS%penSWflux_diag(i,j,k=1)    is penetrative shortwave at top of ocean.
      ! CS%penSWflux_diag(i,j,k=kbot+1) is zero, since assume no SW penetrates rock.
      ! CS%penSWflux_diag = rsdo  and CS%penSW_diag = rsdoabsorb
      ! rsdoabsorb(k) = rsdo(k) - rsdo(k+1), so that rsdo(k) = rsdo(k+1) + rsdoabsorb(k)
      if(CS%id_penSWflux_diag > 0) then
        do k=nz,1,-1 ; do i=is,ie
          CS%penSWflux_diag(i,j,k) = CS%penSW_diag(i,j,k) + CS%penSWflux_diag(i,j,k+1)
        enddo ; enddo
      endif

    endif

    ! Fill CS%nonpenSW_diag
    if(CS%id_nonpenSW_diag > 0) then
      do i=is,ie
        CS%nonpenSW_diag(i,j) = nonpenSW(i)
      enddo
    endif

  enddo ! j-loop finish

  ! Post the diagnostics
  if (CS%id_createdH       > 0) call post_data(CS%id_createdH      , CS%createdH      , CS%diag)
  if (CS%id_penSW_diag     > 0) call post_data(CS%id_penSW_diag    , CS%penSW_diag    , CS%diag)
  if (CS%id_penSWflux_diag > 0) call post_data(CS%id_penSWflux_diag, CS%penSWflux_diag, CS%diag)
  if (CS%id_nonpenSW_diag  > 0) call post_data(CS%id_nonpenSW_diag , CS%nonpenSW_diag , CS%diag)

  if (numberOfGroundings>0) then
    do i = 1, min(numberOfGroundings, maxGroundings)
      call forcing_SinglePointPrint(fluxes,G,iGround(i),jGround(i),'applyBoundaryFluxesInOut (grounding)')
      write(mesg(1:45),'(3es15.3)') G%geoLonT( iGround(i), jGround(i) ), &
                             G%geoLatT( iGround(i), jGround(i)) , hGrounding(i)
      call MOM_error(WARNING, "MOM_diabatic_driver.F90, applyBoundaryFluxesInOut(): "//&
                              "Mass created. x,y,dh= "//trim(mesg), all_print=.true.)
    enddo

    if (numberOfGroundings - maxGroundings > 0) then
      write(mesg, '(i4)') numberOfGroundings - maxGroundings
      call MOM_error(WARNING, "MOM_diabatic_driver:F90, applyBoundaryFluxesInOut(): "//&
                              trim(mesg) // " groundings remaining")
    endif
  endif

end subroutine applyBoundaryFluxesInOut


subroutine diabatic_aux_init(Time, G, GV, param_file, diag, CS, useALEalgorithm, use_ePBL)
  type(time_type),         intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(verticalGrid_type),               intent(in)    :: GV
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ctrl), target, intent(inout) :: diag
  type(diabatic_aux_CS),   pointer       :: CS
  logical,                 intent(in)    :: useALEalgorithm
  logical,                 intent(in)    :: use_ePBL

! Arguments:
!  (in)     Time       = current model time
!  (in)     G          = ocean grid structure
!  (in)     GV - The ocean's vertical grid structure.
!  (in)     param_file = structure indicating the open file to parse for parameter values
!  (in)     diag       = structure used to regulate diagnostic output
!  (in/out) CS         = pointer set to point to the control structure for this module
!  (in)     use_ePBL   = If true, use the implicit energetics planetary boundary
!                        layer scheme to determine the diffusivity in the
!                        surface boundary layer.
  type(vardesc) :: vd

! This "include" declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod  = "MOM_diabatic_aux" ! This module's name.
  character(len=48)  :: thickness_units
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz, nbands
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "diabatic_aux_init called with an "// &
                            "associated control structure.")
    return
  else
    allocate(CS)
   endif

  CS%diag => diag

! Set default, read and log parameters
  call log_version(param_file, mod, version, &
                   "The following parameters are used for auxiliary diabatic processes.")

  call get_param(param_file, mod, "RECLAIM_FRAZIL", CS%reclaim_frazil, &
                 "If true, try to use any frazil heat deficit to cool any\n"//&
                 "overlying layers down to the freezing point, thereby \n"//&
                 "avoiding the creation of thin ice when the SST is above \n"//&
                 "the freezing point.", default=.true.)
  call get_param(param_file, mod, "PRESSURE_DEPENDENT_FRAZIL", &
                                CS%pressure_dependent_frazil, &
                 "If true, use a pressure dependent freezing temperature \n"//&
                 "when making frazil. The default is false, which will be \n"//&
                 "faster but is inappropriate with ice-shelf cavities.", &
                 default=.false.)
  call get_param(param_file, mod, "MINIMUM_FORCING_DEPTH", CS%minimum_forcing_depth, &
                 "The smallest depth over which forcing can be applied. This\n"//&
                 "only takes effect when near-surface layers become thin\n"//&
                 "relative to this scale, in which case the forcing tendencies\n"//&
                 "scaled down by distributing the forcing over this depth scale.", &
                 units="m", default=0.001)
  call get_param(param_file, mod, "EVAP_CFL_LIMIT", CS%evap_CFL_limit, &
                 "The largest fraction of a layer than can be lost to forcing\n"//&
                 "(e.g. evaporation, sea-ice formation) in one time-step. The unused\n"//&
                 "mass loss is passed down through the column.", &
                 units="nondim", default=0.8)

  if (use_ePBL) then
    call get_param(param_file, mod, "DO_RIVERMIX", CS%do_rivermix, &
                 "If true, apply additional mixing whereever there is \n"//&
                 "runoff, so that it is mixed down to RIVERMIX_DEPTH \n"//&
                 "if the ocean is that deep.", default=.false.)
    if (CS%do_rivermix) &
      call get_param(param_file, mod, "RIVERMIX_DEPTH", CS%rivermix_depth, &
                 "The depth to which rivers are mixed if DO_RIVERMIX is \n"//&
                 "defined.", units="m", default=0.0)
  else ; CS%do_rivermix = .false. ; CS%rivermix_depth = 0.0 ; endif
  if (GV%nkml == 0) then
    call get_param(param_file, mod, "USE_RIVER_HEAT_CONTENT", CS%use_river_heat_content, &
                   "If true, use the fluxes%runoff_Hflx field to set the \n"//&
                   "heat carried by runoff, instead of using SST*CP*liq_runoff.", &
                   default=.false.)
    call get_param(param_file, mod, "USE_CALVING_HEAT_CONTENT", CS%use_calving_heat_content, &
                   "If true, use the fluxes%calving_Hflx field to set the \n"//&
                   "heat carried by runoff, instead of using SST*CP*froz_runoff.", &
                   default=.false.)
  else
    CS%use_river_heat_content = .false.
    CS%use_calving_heat_content = .false.
  endif

  if (useALEalgorithm) then
    CS%id_createdH = register_diag_field('ocean_model',"created_H",diag%axesT1, &
        Time, "The volume flux added to stop the ocean from drying out and becoming negative in depth", &
        "meter second-1")
    if (CS%id_createdH>0) allocate(CS%createdH(isd:ied,jsd:jed))

    ! diagnostic for heating of a grid cell from convergence of SW heat into the cell
    CS%id_penSW_diag = register_diag_field('ocean_model', 'rsdoabsorb',                     &
          diag%axesTL, Time, 'Convergence of Penetrative Shortwave Flux in Sea Water Layer',&
          'Watt meter-2', standard_name='net_rate_of_absorption_of_shortwave_energy_in_ocean_layer')

    ! diagnostic for penetrative SW heat flux at top interface of tracer cell (nz+1 interfaces)
    ! k=1 gives penetrative SW at surface; SW(k=nz+1)=0 (no penetration through rock).
    CS%id_penSWflux_diag = register_diag_field('ocean_model', 'rsdo',                               &
          diag%axesTi, Time, 'Downwelling Shortwave Flux in Sea Water at Grid Cell Upper Interface',&
          'Watt meter-2', standard_name='downwelling_shortwave_flux_in_sea_water')

    ! need both arrays for the SW diagnostics (one for flux, one for convergence)
    if (CS%id_penSW_diag>0 .or. CS%id_penSWflux_diag>0) then
       allocate(CS%penSW_diag(isd:ied,jsd:jed,nz))
       CS%penSW_diag(:,:,:) = 0.0
       allocate(CS%penSWflux_diag(isd:ied,jsd:jed,nz+1))
       CS%penSWflux_diag(:,:,:) = 0.0
    endif

    ! diagnostic for non-downwelling SW radiation (i.e., SW absorbed at ocean surface)
    CS%id_nonpenSW_diag = register_diag_field('ocean_model', 'nonpenSW',                       &
          diag%axesT1, Time,                                                                   &
          'Non-downwelling SW radiation (i.e., SW absorbed in ocean surface with LW,SENS,LAT)',&
          'Watt meter-2', standard_name='nondownwelling_shortwave_flux_in_sea_water')
    if (CS%id_nonpenSW_diag > 0) then
       allocate(CS%nonpenSW_diag(isd:ied,jsd:jed))
       CS%nonpenSW_diag(:,:) = 0.0
    endif
  endif

  id_clock_uv_at_h = cpu_clock_id('(Ocean find_uv_at_h)', grain=CLOCK_ROUTINE)
  id_clock_frazil  = cpu_clock_id('(Ocean frazil)', grain=CLOCK_ROUTINE)

end subroutine diabatic_aux_init


subroutine diabatic_aux_end(CS)
  type(diabatic_aux_CS), pointer :: CS

  if (.not.associated(CS)) return

  if (CS%id_createdH       >0) deallocate(CS%createdH)
  if (CS%id_penSW_diag     >0) deallocate(CS%penSW_diag)
  if (CS%id_penSWflux_diag >0) deallocate(CS%penSWflux_diag)
  if (CS%id_nonpenSW_diag  >0) deallocate(CS%nonpenSW_diag)

  if (associated(CS)) deallocate(CS)

end subroutine diabatic_aux_end

end module MOM_diabatic_aux
