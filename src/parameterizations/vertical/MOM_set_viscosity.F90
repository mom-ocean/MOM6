module MOM_set_visc
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
!*  By Robert Hallberg, April 1994 - October 2006                      *
!*  Quadratic Bottom Drag by James Stephens and R. Hallberg.           *
!*                                                                     *
!*    This file contains the subroutine that calculates various values *
!*  related to the bottom boundary layer, such as the viscosity and    *
!*  thickness of the BBL (set_viscous_BBL).  This would also be the    *
!*  module in which other viscous quantities that are flow-independent *
!*  might be set.  This information is transmitted to other modules    *
!*  via a vertvisc type structure.                                     *
!*                                                                     *
!*    The same code is used for the two velocity components, by        *
!*  indirectly referencing the velocities and defining a handful of    *
!*  direction-specific defined variables.                              *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, frhatv, tauy                          *
!*    j    x ^ x ^ x   At >:  u, frhatu, taux                          *
!*    j    > o > o >   At o:  h                                        *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_debugging, only : uchksum, vchksum
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_kappa_shear, only : kappa_shear_is_used
use MOM_CVMix_shear, only : CVMix_shear_is_used
use MOM_io, only : vardesc, var_desc
use MOM_restart, only : register_restart_field, MOM_restart_CS
use MOM_variables, only : thermo_var_ptrs
use MOM_variables, only : vertvisc_type
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>

public set_viscous_BBL, set_viscous_ML, set_visc_init, set_visc_end
public set_visc_register_restarts

type, public :: set_visc_CS ; private
  real    :: Hbbl           ! The static bottom boundary layer thickness, in
                            ! the same units as thickness (m or kg m-2).
  real    :: cdrag          ! The quadratic drag coefficient.
  real    :: c_Smag         ! The Laplacian Smagorinsky coefficient for
                            ! calculating the drag in channels.
  real    :: drag_bg_vel    ! An assumed unresolved background velocity for
                            ! calculating the bottom drag, in m s-1.
  real    :: BBL_thick_min  ! The minimum bottom boundary layer thickness in
                            ! the same units as thickness (m or kg m-2).
                            ! This might be Kv / (cdrag * drag_bg_vel) to give
                            ! Kv as the minimum near-bottom viscosity.
  real    :: Htbl_shelf     ! A nominal thickness of the surface boundary layer
                            ! for use in calculating the near-surface velocity,
                            ! in units of m.
  real    :: Htbl_shelf_min ! The minimum surface boundary layer thickness in m.
  real    :: KV_BBL_min     ! The minimum viscosities in the bottom and top
  real    :: KV_TBL_min     ! boundary layers, both in m2 s-1.

  logical :: bottomdraglaw  ! If true, the  bottom stress is calculated with a
                            ! drag law c_drag*|u|*u. The velocity magnitude
                            ! may be an assumed value or it may be based on the
                            ! actual velocity in the bottommost HBBL, depending
                            ! on whether linear_drag is true.
  logical :: BBL_use_EOS    ! If true, use the equation of state in determining
                            ! the properties of the bottom boundary layer.
  logical :: linear_drag    ! If true, the drag law is cdrag*DRAG_BG_VEL*u.
  logical :: Channel_drag   ! If true, the drag is exerted directly on each
                            ! layer according to what fraction of the bottom
                            ! they overlie.
  logical :: RiNo_mix       ! If true, use Richardson number dependent mixing.
  logical :: dynamic_viscous_ML  ! If true, use a bulk Richardson number criterion to
                            ! determine the mixed layer thickness for viscosity.
  real    :: bulk_Ri_ML     ! The bulk mixed layer used to determine the
                            ! thickness of the viscous mixed layer.  Nondim.
  real    :: omega          !   The Earth's rotation rate, in s-1.
  real    :: ustar_min      ! A minimum value of ustar to avoid numerical
                            ! problems, in m s-1.  If the value is small enough,
                            ! this should not affect the solution.
  real    :: TKE_decay      ! The ratio of the natural Ekman depth to the TKE
                            ! decay scale, nondimensional.
  real    :: omega_frac     !   When setting the decay scale for turbulence, use
                            ! this fraction of the absolute rotation rate blended
                            ! with the local value of f, as sqrt((1-of)*f^2 + of*4*omega^2).
  logical :: debug          ! If true, write verbose checksums for debugging purposes.
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                            ! timing of diagnostic output.
  integer :: id_bbl_thick_u = -1, id_kv_bbl_u = -1
  integer :: id_bbl_thick_v = -1, id_kv_bbl_v = -1
  integer :: id_Ray_u = -1, id_Ray_v = -1
  integer :: id_nkml_visc_u = -1, id_nkml_visc_v = -1

end type set_visc_CS

contains

subroutine set_viscous_BBL(u, v, h, tv, visc, G, GV, CS)
  type(ocean_grid_type),                  intent(inout) :: G
  type(verticalGrid_type),                intent(in)    :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in) :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in) :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h
  type(thermo_var_ptrs),                     intent(in) :: tv
  type(vertvisc_type),                    intent(inout) :: visc
  type(set_visc_CS),                         pointer    :: CS
!   The following subroutine calculates the thickness of the bottom
! boundary layer and the viscosity within that layer.  A drag law is
! used, either linearized about an assumed bottom velocity or using
! the actual near-bottom velocities combined with an assumed
! unresolved velocity.  The bottom boundary layer thickness is
! limited by a combination of stratification and rotation, as in the
! paper of Killworth and Edwards, JPO 1999.  It is not necessary to
! calculate the thickness and viscosity every time step; instead
! previous values may be used.
!
! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m or kg m-2.  In the comments below,
!                the units of h are denoted as H.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (out)     visc - A structure containing vertical viscosities and related
!                   fields.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 vertvisc_init.
  real, dimension(SZIB_(G)) :: &
    ustar, &    !   The bottom friction velocity, in m s-1.
    T_EOS, &    !   The temperature used to calculate the partial derivatives
                ! of density with T and S, in deg C.
    S_EOS, &    !   The salinity used to calculate the partial derivatives
                ! of density with T and S, in PSU.
    dR_dT, &    !   Partial derivative of the density in the bottom boundary
                ! layer with temperature, in units of kg m-3 K-1.
    dR_dS, &    !   Partial derivative of the density in the bottom boundary
                ! layer with salinity, in units of kg m-3 psu-1.
    press       !   The pressure at which dR_dT and dR_dS are evaluated, in Pa.
  real :: htot             ! Sum of the layer thicknesses up to some
                           ! point, in H (i.e., m or kg m-2).
  real :: htot_vel         ! Sum of the layer thicknesses up to some
                           ! point, in H (i.e., m or kg m-2).

  real :: Rhtot            ! Running sum of thicknesses times the
                           ! layer potential densities in H kg m-3.
  real :: h_at_vel(SZIB_(G),SZK_(G))! Layer thickness at a velocity point,
                           ! using an upwind-biased second order
                           ! accurate estimate based on the previous
                           ! velocity direction, in H.
  real :: h_vel            ! The arithmetic mean thickness at a velocity point,
                           ! in H.
  real :: ustarsq          ! 400 times the square of ustar, times
                           ! Rho0 divided by G_Earth and the conversion
                           ! from m to thickness units, in kg m-2 or kg2 m-5.
  real :: cdrag_sqrt       ! Square root of the drag coefficient, nd.
  real :: oldfn            ! The integrated energy required to
                           ! entrain up to the bottom of the layer,
                           ! divided by G_Earth, in H kg m-3.
  real :: Dfn              ! The increment in oldfn for entraining
                           ! the layer, in H kg m-3.
  real :: Dh               ! The increment in layer thickness from
                           ! the present layer, in H.
  real :: bbl_thick        ! The thickness of the bottom boundary layer in m.
  real :: Rl               ! The potential density of the current layer and the
  real :: Rla              ! layer above, in kg m-3.
  real :: TLay, Tabove     ! The temperature and salinity of the current layer
  real :: SLay, Sabove     ! and the layer above, in deg C and PSU.
  real :: C2f              ! C2f = 2*f at velocity points.

  real :: U_bg_sq          ! The square of an assumed background
                           ! velocity, for calculating the mean
                           ! magnitude near the bottom for use in the
                           ! quadratic bottom drag, in m2.
  real :: hwtot            ! Sum of the thicknesses used to calculate
                           ! the near-bottom velocity magnitude, in H.
  real :: hutot            ! Running sum of thicknesses times the
                           ! velocity magnitudes, in H m s-1.
  real :: Thtot            ! Running sum of thickness times temperature, in H C.
  real :: Shtot            ! Running sum of thickness times salinity, in H psu.
  real :: hweight          ! The thickness of a layer that is within Hbbl
                           ! of the bottom, in H.
  real :: v_at_u, u_at_v   ! v at a u point or vice versa, m s-1.
  real :: Rho0x400_G       ! 400*Rho0/G_Earth, in kg s2 m-4.  The 400 is a
                           ! constant proposed by Killworth and Edwards, 1999.
  real, dimension(SZI_(G),SZJ_(G),max(GV%nk_rho_varies,1)) :: &
    Rml                    ! The mixed layer coordinate density, in kg m-3.
  real :: p_ref(SZI_(G))   !   The pressure used to calculate the coordinate
                           ! density, in Pa (usually set to 2e7 Pa = 2000 dbar).

  ! The units H in the following are thickness units - typically m or kg m-2.
  real :: D_vel            ! The bottom depth at a velocity point, in H.
  real :: Dp, Dm           ! The depths at the edges of a velocity cell, in H.
  real :: a                ! a is the curvature of the bottom depth across a
                           ! cell, times the cell width squared, in H.
  real :: a_3, a_12, C24_a ! a/3, a/12, and 24/a, in H, H, and H-1.
  real :: slope            ! The absolute value of the bottom depth slope across
                           ! a cell times the cell width, in H.
  real :: apb_4a, ax2_3apb ! Various nondimensional ratios of a and slope.
  real :: a2x48_apb3, Iapb, Ibma_2 ! Combinations of a and slope with units of H-1.
  ! All of the following "volumes" have units of meters as they are normalized
  ! by the full horizontal area of a velocity cell.
  real :: Vol_open         ! The cell volume above which it is open, in H.
  real :: Vol_direct       ! With less than Vol_direct (in H), there is a direct
                           ! solution of a cubic equation for L.
  real :: Vol_2_reg        ! The cell volume above which there are two separate
                           ! open areas that must be integrated, in H.
  real :: vol              ! The volume below the interface whose normalized
                           ! width is being sought, in H.
  real :: vol_below        ! The volume below the interface below the one that
                           ! is currently under consideration, in H.
  real :: Vol_err          ! The error in the volume with the latest estimate of
                           ! L, or the error for the interface below, in H.
  real :: Vol_quit         ! The volume error below which to quit iterating, in H.
  real :: Vol_tol          ! A volume error tolerance, in H.
  real :: L(SZK_(G)+1)     ! The fraction of the full cell width that is open at
                           ! the depth of each interface, nondimensional.
  real :: L_direct         ! The value of L above volume Vol_direct, nondim.
  real :: L_max, L_min     ! Upper and lower bounds on the correct value for L.
  real :: Vol_err_max      ! The volume errors for the upper and lower bounds on
  real :: Vol_err_min      ! the correct value for L, in H.
  real :: Vol_0            ! A deeper volume with known width L0, in H.
  real :: L0               ! The value of L above volume Vol_0, nondim.
  real :: dVol             ! vol - Vol_0, in H.
  real :: dV_dL2           ! The partial derivative of volume with L squared
                           ! evaluated at L=L0, in H.
  real :: h_neglect        ! A thickness that is so small it is usually lost
                           ! in roundoff and can be neglected, in H.
  real :: ustH             ! ustar converted to units of H s-1.
  real :: root             ! A temporary variable with units of H s-1.
  real :: H_to_m, m_to_H   ! Local copies of unit conversion factors.

  real :: Cell_width       ! The transverse width of the velocity cell, in m.
  real :: Rayleigh         ! A nondimensional value that is multiplied by the
                           ! layer's velocity magnitude to give the Rayleigh
                           ! drag velocity.
  real :: gam              ! The ratio of the change in the open interface width
                           ! to the open interface width atop a cell, nondim.
  real :: BBL_frac         ! The fraction of a layer's drag that goes into the
                           ! viscous bottom boundary layer, nondim.
  real :: BBL_visc_frac    ! The fraction of all the drag that is expressed as
                           ! a viscous bottom boundary layer, nondim.
  real, parameter :: C1_3 = 1.0/3.0, C1_6 = 1.0/6.0, C1_12 = 1.0/12.0
  real :: C2pi_3           ! An irrational constant, 2/3 pi.
  real :: tmp              ! A temporary variable.
  logical :: use_BBL_EOS, do_i(SZIB_(G))
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, m, K2, nkmb, nkml
  integer :: itt, maxitt=20
  real :: tmp_val_m1_to_p1
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nkmb = GV%nk_rho_varies ; nkml = GV%nkml
  h_neglect = GV%H_subroundoff
  Rho0x400_G = 400.0*(GV%Rho0/GV%g_Earth)*GV%m_to_H
  Vol_quit = 0.9*GV%Angstrom + h_neglect
  H_to_m = GV%H_to_m ; m_to_H = GV%m_to_H
  C2pi_3 = 8.0*atan(1.0)/3.0

  if (.not.associated(CS)) call MOM_error(FATAL,"MOM_vert_friction(BBL): "//&
         "Module must be initialized before it is used.")
  if (.not.CS%bottomdraglaw) return

  use_BBL_EOS = associated(tv%eqn_of_state) .and. CS%BBL_use_EOS

  U_bg_sq = CS%drag_bg_vel * CS%drag_bg_vel
  cdrag_sqrt=sqrt(CS%cdrag)
  K2 = max(nkmb+1, 2)

!  With a linear drag law, the friction velocity is already known.
!  if (CS%linear_drag) ustar(:) = cdrag_sqrt*CS%drag_bg_vel

  if ((nkml>0) .and. .not.use_BBL_EOS) then
    do i=G%IscB,G%IecB+1 ; p_ref(i) = tv%P_ref ; enddo
!$OMP parallel do default(none) shared(Jsq,Jeq,Isq,Ieq,nkmb,tv,p_ref,Rml)
    do j=Jsq,Jeq+1 ; do k=1,nkmb
      call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_ref, &
                      Rml(:,j,k), Isq, Ieq-Isq+2, tv%eqn_of_state)
    enddo ; enddo
  endif

!$OMP parallel do default(none) shared(u, v, h, tv, visc, G, GV, CS, Rml, is, ie, js, je,  &
!$OMP                                  nz, Isq, Ieq, Jsq, Jeq, nkmb, h_neglect, Rho0x400_G,&
!$OMP                                  C2pi_3, U_bg_sq, cdrag_sqrt,K2,use_BBL_EOS,         &
!$OMP                                  maxitt,nkml,m_to_H,H_to_m,Vol_quit)                 &
!$OMP                          private(do_i,h_at_vel,htot_vel,hwtot,hutot,Thtot,Shtot,     &
!$OMP                                  hweight,v_at_u,u_at_v,ustar,T_EOS,S_EOS,press,      &
!$OMP                                  dR_dT, dR_dS,ustarsq,htot,TLay,SLay,Tabove,Sabove,  &
!$OMP                                  oldfn,Dfn,Dh,Rhtot,Rla,Rl,C2f,ustH,root,bbl_thick,  &
!$OMP                                  D_vel,tmp,Dp,Dm,a_3,a,a_12,slope,Vol_open,Vol_2_reg,&
!$OMP                                  C24_a,apb_4a,Iapb,a2x48_apb3,ax2_3apb,Vol_direct,   &
!$OMP                                  L_direct,Ibma_2,L,vol,vol_below,Vol_err,            &
!$OMP                                  BBL_visc_frac,h_vel,L0,Vol_0,dV_dL2,dVol,L_max,     &
!$OMP                                  L_min,Vol_err_min,Vol_err_max,BBL_frac,Cell_width,  &
!$OMP                                  gam,Rayleigh, Vol_tol, tmp_val_m1_to_p1)
  do j=G%JscB,G%JecB ; do m=1,2

    if (m==1) then
      if (j<G%Jsc) cycle
      is = G%IscB ; ie = G%IecB
      do i=is,ie
        do_i(i) = .false.
        if (G%mask2dCu(I,j) > 0) do_i(i) = .true.
      enddo
    else
      is = G%isc ; ie = G%iec
      do i=is,ie
        do_i(i) = .false.
        if (G%mask2dCv(i,J) > 0) do_i(i) = .true.
      enddo
    endif

    if (m==1) then
      do k=1,nz ; do i=is,ie ; if (do_i(i)) then
        if (u(I,j,k) *(h(i+1,j,k) - h(i,j,k)) >= 0) then
          h_at_vel(i,k) = 2.0*h(i,j,k)*h(i+1,j,k) / &
                          (h(i,j,k) + h(i+1,j,k) + h_neglect)
        else
          h_at_vel(i,k) =  0.5 * (h(i,j,k) + h(i+1,j,k))
        endif
      endif ; enddo ; enddo
    else
      do k=1,nz ; do i=is,ie ; if (do_i(i)) then
        if (v(i,J,k) * (h(i,j+1,k) - h(i,j,k)) >= 0) then
          h_at_vel(i,k) = 2.0*h(i,j,k)*h(i,j+1,k) / &
                          (h(i,j,k) + h(i,j+1,k) + h_neglect)
        else
          h_at_vel(i,k) = 0.5 * (h(i,j,k) + h(i,j+1,k))
        endif
      endif ; enddo ; enddo
    endif

    if (use_BBL_EOS .or. .not.CS%linear_drag) then
      do i=is,ie ; if (do_i(i)) then
!   This block of code calculates the mean velocity magnitude over
! the bottommost CS%Hbbl of the water column for determining
! the quadratic bottom drag.
        htot_vel = 0.0 ; hwtot = 0.0 ; hutot = 0.0
        Thtot = 0.0 ; Shtot = 0.0
        do k=nz,1,-1

          if (htot_vel>=CS%Hbbl) exit ! terminate the k loop

          hweight = MIN(CS%Hbbl - htot_vel, h_at_vel(i,k))
          if (hweight < 1.5*GV%Angstrom + h_neglect) cycle

          htot_vel  = htot_vel + h_at_vel(i,k)
          hwtot = hwtot + hweight

          if ((.not.CS%linear_drag) .and. (hweight >= 0.0)) then ; if (m==1) then
            v_at_u = set_v_at_u(v, h, G, i, j, k)
            hutot = hutot + hweight * sqrt(u(I,j,k)*u(I,j,k) + &
                                           v_at_u*v_at_u + U_bg_sq)
          else
            u_at_v = set_u_at_v(u, h, G, i, j, k)
            hutot = hutot + hweight * sqrt(v(i,J,k)*v(i,J,k) + &
                                           u_at_v*u_at_v + U_bg_sq)
          endif ; endif

          if (use_BBL_EOS .and. (hweight >= 0.0)) then ; if (m==1) then
            Thtot = Thtot + hweight * 0.5 * (tv%T(i,j,k) + tv%T(i+1,j,k))
            Shtot = Shtot + hweight * 0.5 * (tv%S(i,j,k) + tv%S(i+1,j,k))
          else
            Thtot = Thtot + hweight * 0.5 * (tv%T(i,j,k) + tv%T(i,j+1,k))
            Shtot = Shtot + hweight * 0.5 * (tv%S(i,j,k) + tv%S(i,j+1,k))
          endif ; endif
        enddo ! end of k loop

        if (.not.CS%linear_drag .and. (hwtot > 0.0)) then
          ustar(i) = cdrag_sqrt*hutot/hwtot
        else
          ustar(i) = cdrag_sqrt*CS%drag_bg_vel
        endif

        if (use_BBL_EOS) then ; if (hwtot > 0.0) then
          T_EOS(i) = Thtot/hwtot ; S_EOS(i) = Shtot/hwtot
        else
          T_EOS(i) = 0.0 ; S_EOS(i) = 0.0
        endif ; endif
      endif ; enddo
    else
      do i=is,ie ; ustar(i) = cdrag_sqrt*CS%drag_bg_vel ; enddo
    endif ! Not linear_drag

    if (use_BBL_EOS) then
      do i=is,ie
        press(i) = 0.0 ! or = fluxes%p_surf(i,j)
        if (.not.do_i(i)) then ; T_EOS(i) = 0.0 ; S_EOS(i) = 0.0 ; endif
      enddo
      if (m==1) then ; do k=1,nz ; do i=is,ie
        press(I) = press(I) + GV%H_to_Pa * 0.5 * (h(i,j,k) + h(i+1,j,k))
      enddo ; enddo ; else ; do k=1,nz ; do i=is,ie
        press(i) = press(i) + GV%H_to_Pa * 0.5 * (h(i,j,k) + h(i,j+1,k))
      enddo ; enddo ; endif
      call calculate_density_derivs(T_EOS, S_EOS, press, dR_dT, dR_dS, &
                                    is-G%IsdB+1, ie-is+1, tv%eqn_of_state)
    endif

    do i=is,ie ; if (do_i(i)) then
!  The 400.0 in this expression is the square of a constant proposed
!  by Killworth and Edwards, 1999, in equation (2.20).
      ustarsq = Rho0x400_G * ustar(i)**2
      htot = 0.0

!   This block of code calculates the thickness of a stratification
! limited bottom boundary layer, using the prescription from
! Killworth and Edwards, 1999, as described in Stephens and Hallberg
! 2000 (unpublished and lost manuscript).
      if (use_BBL_EOS) then
        Thtot = 0.0 ; Shtot = 0.0 ; oldfn = 0.0
        do k=nz,2,-1
          if (h_at_vel(i,k) <= 0.0) cycle

          if (m==1) then ! Work on a u-point
            ! Perhaps these should be thickness weighted?
            TLay = 0.5 * (tv%T(i,j,k) + tv%T(i+1,j,k))
            SLay = 0.5 * (tv%S(i,j,k) + tv%S(i+1,j,k))
            Tabove = 0.5 * (tv%T(i,j,k-1) + tv%T(i+1,j,k-1))
            Sabove = 0.5 * (tv%S(i,j,k-1) + tv%S(i+1,j,k-1))
          else ! Work on a v-point
            TLay = 0.5 * (tv%T(i,j,k) + tv%T(i,j+1,k))
            SLay = 0.5 * (tv%S(i,j,k) + tv%S(i,j+1,k))
            Tabove = 0.5 * (tv%T(i,j,k-1) + tv%T(i,j+1,k-1))
            Sabove = 0.5 * (tv%S(i,j,k-1) + tv%S(i,j+1,k-1))
          endif
          oldfn = dR_dT(i)*(Thtot - TLay*htot) + dR_dS(i)*(Shtot - SLay*htot)
          if (oldfn >= ustarsq) exit

          Dfn = (dR_dT(i)*(TLay - Tabove) + dR_dS(i)*(SLay-Sabove)) * &
                (h_at_vel(i,k)+htot)

          if ((oldfn + Dfn) <= ustarsq) then
            Dh = h_at_vel(i,k)
          else
            Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn)/Dfn)
          endif

          htot = htot + Dh
          Thtot = Thtot + TLay*Dh ; Shtot = Shtot + SLay*Dh
        enddo
        if ((oldfn < ustarsq) .and. h_at_vel(i,1) > 0.0) then
          ! Layer 1 might be part of the BBL.
          if (m==1) then ! Work on a u-point
            TLay = 0.5 * (tv%T(i,j,1) + tv%T(i+1,j,1))
            SLay = 0.5 * (tv%S(i,j,1) + tv%S(i+1,j,1))
          else ! Work on a v-point
            TLay = 0.5 * (tv%T(i,j,1) +tv%T(i,j+1,1))
            SLay = 0.5 * (tv%S(i,j,1) +tv%S(i,j+1,1))
          endif
          if (dR_dT(i)*(Thtot-TLay*htot) + dR_dS(i)*(Shtot-SLay*htot) < ustarsq) &
            htot = htot + h_at_vel(i,1)
        endif ! Examination of layer 1.
      else  ! Use Rlay and/or the coordinate density as density variables.
        Rhtot = 0.0
        do k=nz,K2,-1
          oldfn = Rhtot - GV%Rlay(k)*htot
          Dfn = (GV%Rlay(k) - GV%Rlay(k-1))*(h_at_vel(i,k)+htot)

          if (oldfn >= ustarsq) then
            cycle
          else if ((oldfn + Dfn) <= ustarsq) then
            Dh = h_at_vel(i,k)
          else
            Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn)/Dfn)
          endif

          htot = htot + Dh
          Rhtot = Rhtot + GV%Rlay(k)*Dh
        enddo
        if (nkml>0) then
          if (m==1) then
            Rla = 0.5*(Rml(i,j,nkmb) + Rml(i+1,j,nkmb))
          else
            Rla = 0.5*(Rml(i,j,nkmb) + Rml(i,j+1,nkmb))
          endif
          do k=nkmb,2,-1
            Rl = Rla
            if (m==1) then
              Rla = 0.5*(Rml(i,j,k-1) + Rml(i+1,j,k-1))
            else
              Rla = 0.5*(Rml(i,j,k-1) + Rml(i,j+1,k-1))
            endif

            oldfn = Rhtot - Rl*htot
            Dfn = (Rl - Rla)*(h_at_vel(i,k)+htot)

            if (oldfn >= ustarsq) then
              cycle
            else if ((oldfn + Dfn) <= ustarsq) then
              Dh = h_at_vel(i,k)
            else
              Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn)/Dfn)
            endif

            htot = htot + Dh
            Rhtot = Rhtot + Rl*Dh
          enddo
          if (Rhtot - Rla*htot < ustarsq) htot = htot + h_at_vel(i,1)
        else
          if (Rhtot - GV%Rlay(1)*htot < ustarsq) htot = htot + h_at_vel(i,1)
        endif
      endif ! use_BBL_EOS

! The Coriolis limit is 0.5*ustar/f. The buoyancy limit here is htot.
! The  bottom boundary layer thickness is found by solving the same
! equation as in Killworth and Edwards:    (h/h_f)^2 + h/h_N = 1.

      if (m==1) then ; C2f = (G%CoriolisBu(I,J-1)+G%CoriolisBu(I,J))
      else ; C2f = (G%CoriolisBu(I-1,J)+G%CoriolisBu(I,J)) ; endif

      if (CS%cdrag * U_bg_sq <= 0.0) then
        ! This avoids NaNs and overflows, and could be used in all cases,
        ! but is not bitwise identical to the current code.
        ustH = ustar(i)*m_to_H ; root = sqrt(0.25*ustH**2 + (htot*C2f)**2)
        if (htot*ustH <= (CS%BBL_thick_min+h_neglect) * (0.5*ustH + root)) then
          bbl_thick = CS%BBL_thick_min
        else
          bbl_thick = (htot * ustH) / (0.5*ustH + root)
        endif
      else
        bbl_thick = htot / (0.5 + sqrt(0.25 + htot*htot*C2f*C2f/ &
          ((ustar(i)*ustar(i)) * (m_to_H**2) )))

        if (bbl_thick < CS%BBL_thick_min) bbl_thick = CS%BBL_thick_min
      endif
! If there is Richardson number dependent mixing, that determines
! the vertical extent of the bottom boundary layer, and there is no
! need to set that scale here.  In fact, viscously reducing the
! shears over an excessively large region reduces the efficacy of
! the Richardson number dependent mixing.
      if ((bbl_thick > 0.5*CS%Hbbl) .and. (CS%RiNo_mix)) bbl_thick = 0.5*CS%Hbbl

      if (CS%Channel_drag) then
        ! The drag within the bottommost bbl_thick is applied as a part of
        ! an enhanced bottom viscosity, while above this the drag is applied
        ! directly to the layers in question as a Rayleigh drag term.
        if (m==1) then
          D_vel = 0.5*(G%bathyT(i,j) + G%bathyT(i+1,j))
          tmp = G%mask2dCu(I,j+1) * 0.5*(G%bathyT(i,j+1) + G%bathyT(i+1,j+1))
          Dp = 2.0 * D_vel * tmp / (D_vel + tmp)
          tmp = G%mask2dCu(I,j-1) * 0.5*(G%bathyT(i,j-1) + G%bathyT(i+1,j-1))
          Dm = 2.0 * D_vel * tmp / (D_vel + tmp)
        else
          D_vel = 0.5*(G%bathyT(i,j) + G%bathyT(i,j+1))
          tmp = G%mask2dCv(i+1,J) * 0.5*(G%bathyT(i+1,j) + G%bathyT(i+1,j+1))
          Dp = 2.0 * D_vel * tmp / (D_vel + tmp)
          tmp = G%mask2dCv(i,J-1) * 0.5*(G%bathyT(i-1,j) + G%bathyT(i-1,j+1))
          Dm = 2.0 * D_vel * tmp / (D_vel + tmp)
        endif
        if (Dm > Dp) then ; tmp = Dp ; Dp = Dm ; Dm = tmp ; endif

        ! Convert the D's to the units of thickness.
        Dp = m_to_H*Dp ; Dm = m_to_H*Dm ; D_vel = m_to_H*D_vel

        a_3 = (Dp + Dm - 2.0*D_vel) ; a = 3.0*a_3 ; a_12 = 0.25*a_3
        slope = Dp - Dm
        ! If the curvature is small enough, there is no reason not to assume
        ! a uniformly sloping or flat bottom.
        if (abs(a) < 1e-2*(slope + CS%BBL_thick_min)) a = 0.0
        ! Each cell extends from x=-1/2 to 1/2, and has a topography
        ! given by D(x) = a*x^2 + b*x + D - a/12.

        ! Calculate the volume above which the entire cell is open and the
        ! other volumes at which the equation that is solved for L changes.
        if (a > 0.0) then
          if (slope >= a) then
            Vol_open = D_vel - Dm ; Vol_2_reg = Vol_open
          else
            tmp = slope/a
            Vol_open = 0.25*slope*tmp + C1_12*a
            Vol_2_reg = 0.5*tmp**2 * (a - C1_3*slope)
          endif
          ! Define some combinations of a & b for later use.
          C24_a = 24.0/a ; Iapb = 1.0/(a+slope)
          apb_4a = (slope+a)/(4.0*a) ; a2x48_apb3 = (48.0*(a*a))*(Iapb**3)
          ax2_3apb = 2.0*C1_3*a*Iapb
        elseif (a == 0.0) then
          Vol_open = 0.5*slope
          if (slope > 0) Iapb = 1.0/slope
        else ! a < 0.0
          Vol_open = D_vel - Dm
          if (slope >= -a) then
            Iapb = 1.0e30 ; if (slope+a /= 0.0) Iapb = 1.0/(a+slope)
            Vol_direct = 0.0 ; L_direct = 0.0 ; C24_a = 0.0
          else
            C24_a = 24.0/a ; Iapb = 1.0/(a+slope)
            L_direct = 1.0 + slope/a ! L_direct < 1 because a < 0
            Vol_direct = -C1_6*a*L_direct**3
          endif
          Ibma_2 = 2.0 / (slope - a)
        endif

        L(nz+1) = 0.0 ; vol = 0.0 ; Vol_err = 0.0 ; BBL_visc_frac = 0.0
        ! Determine the normalized open length at each interface.
        do K=nz,1,-1
          vol_below = vol
          if (m==1) then ; h_vel = 0.5*(h(i,j,k) + h(i+1,j,k))
          else ; h_vel = 0.5*(h(i,j,k) + h(i,j+1,k)) ; endif
          vol = vol + h_vel
          h_vel = h_vel + h_neglect

          if (vol >= Vol_open) then ; L(K) = 1.0
          elseif (a == 0) then ! The bottom has no curvature.
            L(K) = sqrt(2.0*vol*Iapb)
          elseif (a > 0) then
            ! There may be a minimum depth, and there are
            ! analytic expressions for L for all cases.
            if (vol < Vol_2_reg) then
              ! In this case, there is a contiguous open region and
              !   vol = 0.5*L^2*(slope + a/3*(3-4L)).
               if (a2x48_apb3*vol < 1e-8) then ! Could be 1e-7?
                ! There is a very good approximation here for massless layers.
                L0 = sqrt(2.0*vol*Iapb) ; L(K) = L0*(1.0 + ax2_3apb*L0)
              else
                L(K) = apb_4a * (1.0 - &
                         2.0 * cos(C1_3*acos(a2x48_apb3*vol - 1.0) - C2pi_3))
              endif
              ! To check the answers.
              ! Vol_err = 0.5*(L(K)*L(K))*(slope + a_3*(3.0-4.0*L(K))) - vol
            else ! There are two separate open regions.
              !   vol = slope^2/4a + a/12 - (a/12)*(1-L)^2*(1+2L)
              ! At the deepest volume, L = slope/a, at the top L = 1.
              !L(K) = 0.5 - cos(C1_3*acos(1.0 - C24_a*(Vol_open - vol)) - C2pi_3)
              tmp_val_m1_to_p1 = 1.0 - C24_a*(Vol_open - vol)
              tmp_val_m1_to_p1 = max(-1., min(1., tmp_val_m1_to_p1))
              L(K) = 0.5 - cos(C1_3*acos(tmp_val_m1_to_p1) - C2pi_3)
              ! To check the answers.
              ! Vol_err = Vol_open - a_12*(1.0+2.0*L(K)) * (1.0-L(K))**2 - vol
            endif
          else ! a < 0.
            if (vol <= Vol_direct) then
              ! Both edges of the cell are bounded by walls.
              L(K) = (-0.25*C24_a*vol)**C1_3
            else
              ! x_R is at 1/2 but x_L is in the interior & L is found by solving
              !   vol = 0.5*L^2*(slope + a/3*(3-4L))

              !  Vol_err = 0.5*(L(K+1)*L(K+1))*(slope + a_3*(3.0-4.0*L(K+1))) - vol_below
              ! Change to ...
              !   if (min(Vol_below + Vol_err, vol) <= Vol_direct) then ?
              if (vol_below + Vol_err <= Vol_direct) then
                L0 = L_direct ; Vol_0 = Vol_direct
              else
                L0 = L(K+1) ; Vol_0 = Vol_below + Vol_err
                ! Change to   Vol_0 = min(Vol_below + Vol_err, vol) ?
              endif

              !   Try a relatively simple solution that usually works well
              ! for massless layers.
              dV_dL2 = 0.5*(slope+a) - a*L0 ; dVol = (vol-Vol_0)
           !  dV_dL2 = 0.5*(slope+a) - a*L0 ; dVol = max(vol-Vol_0, 0.0)

           !### The following code is more robust when GV%Angstrom=0, but it
           !### changes answers.
           !   Vol_tol = max(0.5*GV%Angstrom + GV%H_subroundoff, 1e-14*vol)
           !   Vol_quit = max(0.9*GV%Angstrom + GV%H_subroundoff, 1e-14*vol)

           !   if (dVol <= 0.0) then
           !     L(K) = L0
           !     Vol_err = 0.5*(L(K)*L(K))*(slope + a_3*(3.0-4.0*L(K))) - vol
           !   elseif (a*a*dVol**3 < Vol_tol*dV_dL2**2 * &
           !                     (dV_dL2*Vol_tol - 2.0*a*L0*dVol)) then
              if (a*a*dVol**3 < GV%Angstrom*dV_dL2**2 * &
                                (0.25*dV_dL2*GV%Angstrom - a*L0*dVol)) then
                ! One iteration of Newton's method should give an estimate
                ! that is accurate to within Vol_tol.
                L(K) = sqrt(L0*L0 + dVol / dV_dL2)
                Vol_err = 0.5*(L(K)*L(K))*(slope + a_3*(3.0-4.0*L(K))) - vol
              else
                if (dV_dL2*(1.0-L0*L0) < dVol + &
                    dV_dL2 * (Vol_open - Vol)*Ibma_2) then
                  L_max = sqrt(1.0 - (Vol_open - Vol)*Ibma_2)
                else
                  L_max = sqrt(L0*L0 + dVol / dV_dL2)
                endif
                L_min = sqrt(L0*L0 + dVol / (0.5*(slope+a) - a*L_max))

                Vol_err_min = 0.5*(L_min**2)*(slope + a_3*(3.0-4.0*L_min)) - vol
                Vol_err_max = 0.5*(L_max**2)*(slope + a_3*(3.0-4.0*L_max)) - vol
           !    if ((abs(Vol_err_min) <= Vol_quit) .or. (Vol_err_min >= Vol_err_max)) then
                if (abs(Vol_err_min) <= Vol_quit) then
                  L(K) = L_min ; Vol_err = Vol_err_min
                else
                  L(K) = sqrt((L_min**2*Vol_err_max - L_max**2*Vol_err_min) / &
                              (Vol_err_max - Vol_err_min))
                  do itt=1,maxitt
                    Vol_err = 0.5*(L(K)*L(K))*(slope + a_3*(3.0-4.0*L(K))) - vol
                    if (abs(Vol_err) <= Vol_quit) exit
                    ! Take a Newton's method iteration. This equation has proven
                    ! robust enough not to need bracketing.
                    L(K) = L(K) - Vol_err / (L(K)* (slope + a - 2.0*a*L(K)))
                    ! This would be a Newton's method iteration for L^2:
                    !   L(K) = sqrt(L(K)*L(K) - Vol_err / (0.5*(slope+a) - a*L(K)))
                  enddo
                endif ! end of iterative solver
              endif ! end of 1-boundary alternatives.
            endif ! end of a<0 cases.
          endif

          ! Determine the drag contributing to the bottom boundary layer
          ! and the Raleigh drag that acting on each layer.
          if (L(K) > L(K+1)) then
            if (vol_below < bbl_thick) then
              BBL_frac = (1.0-vol_below/bbl_thick)**2
              BBL_visc_frac = BBL_visc_frac + BBL_frac*(L(K) - L(K+1))
            else
              BBL_frac = 0.0
            endif

            if (m==1) then ; Cell_width = G%dy_Cu(I,j)
            else ; Cell_width = G%dx_Cv(i,J) ; endif
            gam = 1.0 - L(K+1)/L(K)
            Rayleigh = CS%cdrag * (L(K)-L(K+1)) * (1.0-BBL_frac) * &
                (12.0*CS%c_Smag*h_vel) /  (12.0*CS%c_Smag*h_vel + m_to_H * &
                 CS%cdrag * gam*(1.0-gam)*(1.0-1.5*gam) * L(K)**2 * Cell_width)
          else ! This layer feels no drag.
            Rayleigh = 0.0
          endif

          if (m==1) then
            if (Rayleigh > 0.0) then
              v_at_u = set_v_at_u(v, h, G, i, j, k)
              visc%Ray_u(I,j,k) = Rayleigh*sqrt(u(I,j,k)*u(I,j,k) + &
                                                v_at_u*v_at_u + U_bg_sq)
            else ; visc%Ray_u(I,j,k) = 0.0 ; endif
          else
            if (Rayleigh > 0.0) then
              u_at_v = set_u_at_v(u, h, G, i, j, k)
              visc%Ray_v(i,J,k) = Rayleigh*sqrt(v(i,J,k)*v(i,J,k) + &
                                                u_at_v*u_at_v + U_bg_sq)
            else ; visc%Ray_v(i,J,k) = 0.0 ; endif
          endif

        enddo ! k loop to determine L(K).

        bbl_thick = bbl_thick * H_to_m
        if (m==1) then
          visc%kv_bbl_u(I,j) = max(CS%KV_BBL_min, &
                                   cdrag_sqrt*ustar(i)*bbl_thick*BBL_visc_frac)
          visc%bbl_thick_u(I,j) = bbl_thick
        else
          visc%kv_bbl_v(i,J) = max(CS%KV_BBL_min, &
                                   cdrag_sqrt*ustar(i)*bbl_thick*BBL_visc_frac)
          visc%bbl_thick_v(i,J) = bbl_thick
        endif

      else ! Not Channel_drag.
!   Here the near-bottom viscosity is set to a value which will give
! the correct stress when the shear occurs over bbl_thick.
        bbl_thick = bbl_thick * H_to_m
        if (m==1) then
          visc%kv_bbl_u(I,j) = max(CS%KV_BBL_min, cdrag_sqrt*ustar(i)*bbl_thick)
          visc%bbl_thick_u(I,j) = bbl_thick
        else
          visc%kv_bbl_v(i,J) = max(CS%KV_BBL_min, cdrag_sqrt*ustar(i)*bbl_thick)
          visc%bbl_thick_v(i,J) = bbl_thick
        endif
      endif
    endif ; enddo ! end of i loop
  enddo ; enddo ! end of m & j loops

! Offer diagnostics for averaging
  if (CS%id_bbl_thick_u > 0) &
    call post_data(CS%id_bbl_thick_u, visc%bbl_thick_u, CS%diag)
  if (CS%id_kv_bbl_u > 0) &
    call post_data(CS%id_kv_bbl_u, visc%kv_bbl_u, CS%diag)
  if (CS%id_bbl_thick_v > 0) &
    call post_data(CS%id_bbl_thick_v, visc%bbl_thick_v, CS%diag)
  if (CS%id_kv_bbl_v > 0) &
    call post_data(CS%id_kv_bbl_v, visc%kv_bbl_v, CS%diag)
  if (CS%id_Ray_u > 0) &
    call post_data(CS%id_Ray_u, visc%Ray_u, CS%diag)
  if (CS%id_Ray_v > 0) &
    call post_data(CS%id_Ray_v, visc%Ray_v, CS%diag)

  if (CS%debug) then
    if (associated(visc%Ray_u)) call uchksum(visc%Ray_u,"Ray u",G%HI,haloshift=0)
    if (associated(visc%Ray_v)) call vchksum(visc%Ray_v,"Ray v",G%HI,haloshift=0)
    if (associated(visc%kv_bbl_u)) call uchksum(visc%kv_bbl_u,"kv_bbl_u",G%HI,haloshift=0)
    if (associated(visc%kv_bbl_v)) call vchksum(visc%kv_bbl_v,"kv_bbl_v",G%HI,haloshift=0)
    if (associated(visc%bbl_thick_u)) call uchksum(visc%bbl_thick_u,"bbl_thick_u",G%HI,haloshift=0)
    if (associated(visc%bbl_thick_v)) call vchksum(visc%bbl_thick_v,"bbl_thick_v",G%HI,haloshift=0)
  endif

end subroutine set_viscous_BBL

function set_v_at_u(v, h, G, i, j, k)
  type(ocean_grid_type),                     intent(in) :: G
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in) :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h
  integer,                                   intent(in) :: i, j, k
  real                                                  :: set_v_at_u
  ! This subroutine finds a thickness-weighted value of v at the u-points.
  real :: hwt(4)           ! Masked weights used to average v onto u, in H.
  real :: hwt_tot          ! The sum of the masked thicknesses, in H.

  hwt(1) = (h(i,j-1,k) + h(i,j,k)) * G%mask2dCv(i,J-1)
  hwt(2) = (h(i+1,j-1,k) + h(i+1,j,k)) * G%mask2dCv(i+1,J-1)
  hwt(3) = (h(i,j,k) + h(i,j+1,k)) * G%mask2dCv(i,J)
  hwt(4) = (h(i+1,j,k) + h(i+1,j+1,k)) * G%mask2dCv(i+1,J)
  hwt_tot = (hwt(1) + hwt(4)) + (hwt(2) + hwt(3))
  set_v_at_u = 0.0
  if (hwt_tot > 0.0) set_v_at_u = &
          ((hwt(3) * v(i,J,k) + hwt(2) * v(i+1,J-1,k)) + &
           (hwt(4) * v(i+1,J,k) + hwt(1) * v(i,J-1,k))) / hwt_tot
end function set_v_at_u

function set_u_at_v(u, h, G, i, j, k)
  type(ocean_grid_type),                  intent(in) :: G
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in) :: u
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h
  integer,                                intent(in) :: i, j, k
  real                                               :: set_u_at_v
  ! This subroutine finds a thickness-weighted value of u at the v-points.
  real :: hwt(4)           ! Masked weights used to average u onto v, in H.
  real :: hwt_tot          ! The sum of the masked thicknesses, in H.

  hwt(1) = (h(i-1,j,k) + h(i,j,k)) * G%mask2dCu(I-1,j)
  hwt(2) = (h(i,j,k) + h(i+1,j,k)) * G%mask2dCu(I,j)
  hwt(3) = (h(i-1,j+1,k) + h(i,j+1,k)) * G%mask2dCu(I-1,j+1)
  hwt(4) = (h(i,j+1,k) + h(i+1,j+1,k)) * G%mask2dCu(I,j+1)
  hwt_tot = (hwt(1) + hwt(4)) + (hwt(2) + hwt(3))
  set_u_at_v = 0.0
  if (hwt_tot > 0.0) set_u_at_v = &
          ((hwt(2) * u(I,j,k) + hwt(3) * u(I-1,j+1,k)) + &
           (hwt(1) * u(I-1,j,k) + hwt(4) * u(I,j+1,k))) / hwt_tot
end function set_u_at_v

subroutine set_viscous_ML(u, v, h, tv, fluxes, visc, dt, G, GV, CS)
  type(ocean_grid_type),               intent(inout) :: G
  type(verticalGrid_type),             intent(in)    :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in) :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in) :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h
  type(thermo_var_ptrs),               intent(in)    :: tv
  type(forcing),                       intent(in)    :: fluxes
  type(vertvisc_type),                 intent(inout) :: visc
  real,                                intent(in)    :: dt
  type(set_visc_CS),                   pointer       :: CS
!   The following subroutine calculates the thickness of the surface boundary
! layer for applying an elevated viscosity.  A bulk Richardson criterion or
! the thickness of the topmost NKML layers (with a bulk mixed layer) are
! currently used.  The thicknesses are given in terms of fractional layers, so
! that this thickness will move as the thickness of the topmost layers change.
!
! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m or kg m-2.  In the comments below,
!                the units of h are denoted as H.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (out)     visc - A structure containing vertical viscosities and related
!                   fields.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 vertvisc_init.

  real, dimension(SZIB_(G)) :: &
    htot, &     !   The total depth of the layers being that are within the
                ! surface mixed layer, in H.
    Thtot, &    !   The integrated temperature of layers that are within the
                ! surface mixed layer, in H degC.
    Shtot, &    !   The integrated salt of layers that are within the
                ! surface mixed layer H PSU.
    Rhtot, &    !   The integrated density of layers that are within the
                ! surface mixed layer,  in H kg m-3.  Rhtot is only used if no
                ! equation of state is used.
    uhtot, &    !   The depth integrated zonal and meridional velocities within
    vhtot, &    ! the surface mixed layer, in H m s-1.
    Idecay_len_TKE, &  ! The inverse of a turbulence decay length scale, in H-1.
    dR_dT, &    !   Partial derivative of the density at the base of layer nkml
                ! (roughly the base of the mixed layer) with temperature, in
                ! units of kg m-3 K-1.
    dR_dS, &    !   Partial derivative of the density at the base of layer nkml
                ! (roughly the base of the mixed layer) with salinity, in units
                ! of kg m-3 psu-1.
    ustar, &    !   The surface friction velocity under ice shelves, in m s-1.
    press, &    ! The pressure at which dR_dT and dR_dS are evaluated, in Pa.
    T_EOS, &    ! T_EOS and S_EOS are the potential temperature and salnity at which dR_dT and dR_dS
    S_EOS       ! which dR_dT and dR_dS are evaluated, in degC and PSU.
  real :: h_at_vel(SZIB_(G),SZK_(G))! Layer thickness at velocity points,
                ! using an upwind-biased second order accurate estimate based
                ! on the previous velocity direction, in H.
  integer :: k_massive(SZIB_(G)) ! The k-index of the deepest layer yet found
                ! that has more than h_tiny thickness and will be in the
                ! viscous mixed layer.
  real :: Uh2   ! The squared magnitude of the difference between the velocity
                ! integrated through the mixed layer and the velocity of the
                ! interior layer layer times the depth of the the mixed layer,
                ! in H2 m2 s-2.
  real :: htot_vel  ! Sum of the layer thicknesses up to some
                    ! point, in H (i.e., m or kg m-2).
  real :: hwtot     ! Sum of the thicknesses used to calculate
                    ! the near-bottom velocity magnitude, in H.
  real :: hutot     ! Running sum of thicknesses times the
                    ! velocity magnitudes, in H m s-1.
  real :: hweight   ! The thickness of a layer that is within Hbbl
                    ! of the bottom, in H.

  real :: hlay      ! The layer thickness at velocity points, in H.
  real :: I_2hlay   ! 1 / 2*hlay, in H-1.
  real :: T_lay     ! The layer temperature at velocity points, in deg C.
  real :: S_lay     ! The layer salinity at velocity points, in PSU.
  real :: Rlay      ! The layer potential density at velocity points, in kg m-3.
  real :: Rlb       ! The potential density of the layer below, in kg m-3.
  real :: v_at_u    ! The meridonal velocity at a zonal velocity point in m s-1.
  real :: u_at_v    ! The zonal velocity at a meridonal velocity point in m s-1.
  real :: gHprime   ! The mixed-layer internal gravity wave speed squared, based
                    ! on the mixed layer thickness and density difference across
                    ! the base of the mixed layer, in m2 s-2.
  real :: RiBulk    ! The bulk Richardson number below which water is in the
                    ! viscous mixed layer, including reduction for turbulent
                    ! decay. Nondimensional.
  real :: dt_Rho0   ! The time step divided by the conversion from the layer
                    ! thickness to layer mass, in s H m2 kg-1.
  real :: g_H_Rho0  !   The gravitational acceleration times the conversion from
                    ! H to m divided by the mean density, in m5 s-2 H-1 kg-1.
  real :: ustarsq     ! 400 times the square of ustar, times
                      ! Rho0 divided by G_Earth and the conversion
                      ! from m to thickness units, in kg m-2 or kg2 m-5.
  real :: cdrag_sqrt  ! Square root of the drag coefficient, nd.
  real :: oldfn       ! The integrated energy required to
                      ! entrain up to the bottom of the layer,
                      ! divided by G_Earth, in H kg m-3.
  real :: Dfn         ! The increment in oldfn for entraining
                      ! the layer, in H kg m-3.
  real :: Dh          ! The increment in layer thickness from
                      ! the present layer, in H.
  real :: U_bg_sq   ! The square of an assumed background velocity, for
                    ! calculating the mean magnitude near the top for use in
                    ! the quadratic surface drag, in m2.
  real :: h_tiny    ! A very small thickness, in H. Layers that are less than
                    ! h_tiny can not be the deepest in the viscous mixed layer.
  real :: absf      ! The absolute value of f averaged to velocity points, s-1.
  real :: U_star    ! The friction velocity at velocity points, in m s-1.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in H.
  real :: Rho0x400_G  ! 400*Rho0/G_Earth, in kg s2 m-4.  The 400 is a
                      ! constant proposed by Killworth and Edwards, 1999.
  real :: H_to_m, m_to_H   ! Local copies of unit conversion factors.
  real :: ustar1    ! ustar in units of H/s
  real :: h2f2      ! (h*2*f)^2
  logical :: use_EOS, do_any, do_any_shelf, do_i(SZIB_(G))
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, K2, nkmb, nkml

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nkmb = GV%nk_rho_varies ; nkml = GV%nkml

  if (.not.associated(CS)) call MOM_error(FATAL,"MOM_vert_friction(visc_ML): "//&
         "Module must be initialized before it is used.")
  if (.not.(CS%dynamic_viscous_ML .or. associated(fluxes%frac_shelf_h))) return

  Rho0x400_G = 400.0*(GV%Rho0/GV%g_Earth)*GV%m_to_H
  U_bg_sq = CS%drag_bg_vel * CS%drag_bg_vel
  cdrag_sqrt=sqrt(CS%cdrag)

  use_EOS = associated(tv%eqn_of_state)
  dt_Rho0 = dt/GV%H_to_kg_m2
  h_neglect = GV%H_subroundoff
  h_tiny = 2.0*GV%Angstrom + h_neglect
  g_H_Rho0 = (GV%g_Earth * GV%H_to_m) / GV%Rho0
  H_to_m = GV%H_to_m ; m_to_H = GV%m_to_H

  if (associated(fluxes%frac_shelf_h)) then
    ! This configuration has ice shelves, and the appropriate variables need to
    ! be allocated.
    if (.not.associated(fluxes%frac_shelf_u)) call MOM_error(FATAL, &
      "set_viscous_ML: fluxes%frac_shelf_h is associated, but " // &
      "fluxes%frac_shelf_u is not.")
    if (.not.associated(fluxes%frac_shelf_v)) call MOM_error(FATAL, &
      "set_viscous_ML: fluxes%frac_shelf_h is associated, but " // &
      "fluxes%frac_shelf_v is not.")
    if (.not.associated(visc%taux_shelf)) then
      allocate(visc%taux_shelf(G%IsdB:G%IedB,G%jsd:G%jed))
      visc%taux_shelf(:,:) = 0.0
    endif
    if (.not.associated(visc%tauy_shelf)) then
      allocate(visc%tauy_shelf(G%isd:G%ied,G%JsdB:G%JedB))
      visc%tauy_shelf(:,:) = 0.0
    endif
    if (.not.associated(visc%tbl_thick_shelf_u)) then
      allocate(visc%tbl_thick_shelf_u(G%IsdB:G%IedB,G%jsd:G%jed))
      visc%tbl_thick_shelf_u(:,:) = 0.0
    endif
    if (.not.associated(visc%tbl_thick_shelf_v)) then
      allocate(visc%tbl_thick_shelf_v(G%isd:G%ied,G%JsdB:G%JedB))
      visc%tbl_thick_shelf_v(:,:) = 0.0
    endif
    if (.not.associated(visc%kv_tbl_shelf_u)) then
      allocate(visc%kv_tbl_shelf_u(G%IsdB:G%IedB,G%jsd:G%jed))
      visc%kv_tbl_shelf_u(:,:) = 0.0
    endif
    if (.not.associated(visc%kv_tbl_shelf_v)) then
      allocate(visc%kv_tbl_shelf_v(G%isd:G%ied,G%JsdB:G%JedB))
      visc%kv_tbl_shelf_v(:,:) = 0.0
    endif

    !  With a linear drag law, the friction velocity is already known.
!    if (CS%linear_drag) ustar(:) = cdrag_sqrt*CS%drag_bg_vel
  endif

!$OMP parallel do default(none) shared(u, v, h, tv, fluxes, visc, dt, G, GV, CS, use_EOS, &
!$OMP                                  dt_Rho0, h_neglect, h_tiny, g_H_Rho0,js,je,        &
!$OMP                                  H_to_m, m_to_H, Isq, Ieq, nz, U_bg_sq,             &
!$OMP                                  cdrag_sqrt,Rho0x400_G,nkml) &
!$OMP                          private(do_any,htot,do_i,k_massive,Thtot,uhtot,vhtot,U_Star, &
!$OMP                                  Idecay_len_TKE,press,k2,I_2hlay,T_EOS,S_EOS,dR_dT,   &
!$OMP                                  dR_dS,hlay,v_at_u,Uh2,T_lay,S_lay,gHprime,           &
!$OMP                                  RiBulk,Shtot,Rhtot,absf,do_any_shelf,                &
!$OMP                                  h_at_vel,ustar,htot_vel,hwtot,hutot,hweight,ustarsq, &
!$OMP                                  oldfn,Dfn,Dh,Rlay,Rlb,h2f2,ustar1)
  do j=js,je  ! u-point loop
    if (CS%dynamic_viscous_ML) then
      do_any = .false.
      do I=Isq,Ieq
        htot(I) = 0.0
        if (G%mask2dCu(I,j) < 0.5) then
          do_i(I) = .false. ; visc%nkml_visc_u(I,j) = nkml
        else
          do_i(I) = .true. ; do_any = .true.
          k_massive(I) = nkml
          Thtot(I) = 0.0 ; Shtot(I) = 0.0 ; Rhtot(i) = 0.0
          uhtot(I) = dt_Rho0 * fluxes%taux(I,j)
          vhtot(I) = 0.25 * dt_Rho0 * ((fluxes%tauy(i,J) + fluxes%tauy(i+1,J-1)) + &
                                       (fluxes%tauy(i,J-1) + fluxes%tauy(i+1,J)))

          if (CS%omega_frac >= 1.0) then ; absf = 2.0*CS%omega ; else
            absf = 0.5*(abs(G%CoriolisBu(I,J)) + abs(G%CoriolisBu(I,J-1)))
            if (CS%omega_frac > 0.0) &
              absf = sqrt(CS%omega_frac*4.0*CS%omega**2 + (1.0-CS%omega_frac)*absf**2)
          endif
          U_Star = max(CS%ustar_min, 0.5 * (fluxes%ustar(i,j) + fluxes%ustar(i+1,j)))
          Idecay_len_TKE(I) = ((absf / U_Star) * CS%TKE_decay) * H_to_m
        endif
      enddo

      if (do_any) then ; do k=1,nz
        if (k > nkml) then
          do_any = .false.
          if (use_EOS .and. (k==nkml+1)) then
            ! Find dRho/dT and dRho_dS.
            do I=Isq,Ieq
              press(I) = GV%H_to_Pa * htot(I)
              k2 = max(1,nkml)
              I_2hlay = 1.0 / (h(i,j,k2) + h(i+1,j,k2) + h_neglect)
              T_EOS(I) = (h(i,j,k2)*tv%T(i,j,k2) + h(i+1,j,k2)*tv%T(i+1,j,k2)) * I_2hlay
              S_EOS(I) = (h(i,j,k2)*tv%S(i,j,k2) + h(i+1,j,k2)*tv%S(i+1,j,k2)) * I_2hlay
            enddo
            call calculate_density_derivs(T_EOS, S_EOS, press, dR_dT, dR_dS, &
                                          Isq-G%IsdB+1, Ieq-Isq+1, tv%eqn_of_state)
          endif

          do I=Isq,Ieq ; if (do_i(I)) then

            hlay = 0.5*(h(i,j,k) + h(i+1,j,k))
            if (hlay > h_tiny) then ! Only consider non-vanished layers.
              I_2hlay = 1.0 / (h(i,j,k) + h(i+1,j,k))
              v_at_u = 0.5 * (h(i,j,k)   * (v(i,J,k) + v(i,J-1,k)) + &
                              h(i+1,j,k) * (v(i+1,J,k) + v(i+1,J-1,k))) * I_2hlay
              Uh2 = (uhtot(I) - htot(I)*u(I,j,k))**2 + (vhtot(I) - htot(I)*v_at_u)**2

              if (use_EOS) then
                T_lay = (h(i,j,k)*tv%T(i,j,k) + h(i+1,j,k)*tv%T(i+1,j,k)) * I_2hlay
                S_lay = (h(i,j,k)*tv%S(i,j,k) + h(i+1,j,k)*tv%S(i+1,j,k)) * I_2hlay
                gHprime = g_H_Rho0 * (dR_dT(I) * (T_lay*htot(I) - Thtot(I)) + &
                                      dR_dS(I) * (S_lay*htot(I) - Shtot(I)))
              else
                gHprime = g_H_Rho0 * (GV%Rlay(k)*htot(I) - Rhtot(I))
              endif

              if (gHprime > 0.0) then
                RiBulk = CS%bulk_Ri_ML * exp(-htot(I) * Idecay_len_TKE(I))
                if (RiBulk * Uh2 <= (htot(I)**2) * gHprime) then
                  visc%nkml_visc_u(I,j) = real(k_massive(I))
                  do_i(I) = .false.
                elseif (RiBulk * Uh2 <= (htot(I) + hlay)**2 * gHprime) then
                  visc%nkml_visc_u(I,j) = real(k-1) + &
                    ( sqrt(RiBulk * Uh2 / gHprime) - htot(I) ) / hlay
                  do_i(I) = .false.
                endif
              endif
              k_massive(I) = k
            endif ! hlay > h_tiny

            if (do_i(I)) do_any = .true.
          endif ; enddo

          if (.not.do_any) exit ! All columns are done.
        endif

        do I=Isq,Ieq ; if (do_i(I)) then
          htot(I) = htot(I) + 0.5 * (h(i,j,k) + h(i+1,j,k))
          uhtot(I) = uhtot(I) + 0.5 * (h(i,j,k) + h(i+1,j,k)) * u(I,j,k)
          vhtot(I) = vhtot(I) + 0.25 * (h(i,j,k) * (v(i,J,k) + v(i,J-1,k)) + &
                                        h(i+1,j,k) * (v(i+1,J,k) + v(i+1,J-1,k)))
          if (use_EOS) then
            Thtot(I) = Thtot(I) + 0.5 * (h(i,j,k)*tv%T(i,j,k) + h(i+1,j,k)*tv%T(i+1,j,k))
            Shtot(I) = Shtot(I) + 0.5 * (h(i,j,k)*tv%S(i,j,k) + h(i+1,j,k)*tv%S(i+1,j,k))
          else
            Rhtot(i) = Rhtot(i) + 0.5 * (h(i,j,k) + h(i+1,j,k)) * GV%Rlay(k)
          endif
        endif ; enddo
      enddo ; endif

      if (do_any) then ; do I=Isq,Ieq ; if (do_i(I)) then
        visc%nkml_visc_u(I,j) = k_massive(I)
      endif ; enddo ; endif
    endif ! dynamic_viscous_ML

    do_any_shelf = .false.
    if (associated(fluxes%frac_shelf_h)) then
      do I=Isq,Ieq
        if (fluxes%frac_shelf_u(I,j)*G%mask2dCu(I,j) == 0.0) then
          do_i(I) = .false.
          visc%tbl_thick_shelf_u(I,j) = 0.0 ; visc%kv_tbl_shelf_u(I,j) = 0.0
        else
          do_i(I) = .true. ; do_any_shelf = .true.
        endif
      enddo
    endif

    if (do_any_shelf) then
      do k=1,nz ; do I=Isq,Ieq ; if (do_i(I)) then
        if (u(I,j,k) *(h(i+1,j,k) - h(i,j,k)) >= 0) then
          h_at_vel(i,k) = 2.0*h(i,j,k)*h(i+1,j,k) / &
                          (h(i,j,k) + h(i+1,j,k) + h_neglect)
        else
          h_at_vel(i,k) =  0.5 * (h(i,j,k) + h(i+1,j,k))
        endif
      else
        h_at_vel(I,k) = 0.0 ; ustar(I) = 0.0
      endif ; enddo ; enddo

      do I=Isq,Ieq ; if (do_i(I)) then
        htot_vel = 0.0 ; hwtot = 0.0 ; hutot = 0.0
        Thtot(I) = 0.0 ; Shtot(I) = 0.0
        if (use_EOS .or. .not.CS%linear_drag) then ; do k=1,nz
          if (htot_vel>=CS%Htbl_shelf) exit ! terminate the k loop
          hweight = MIN(CS%Htbl_shelf - htot_vel, h_at_vel(i,k))
          if (hweight <= 1.5*GV%Angstrom + h_neglect) cycle

          htot_vel  = htot_vel + h_at_vel(i,k)
          hwtot = hwtot + hweight

          if (.not.CS%linear_drag) then
            v_at_u = set_v_at_u(v, h, G, i, j, k)
            hutot = hutot + hweight * sqrt(u(I,j,k)**2 + &
                                           v_at_u**2 + U_bg_sq)
          endif
          if (use_EOS) then
            Thtot(I) = Thtot(I) + hweight * 0.5 * (tv%T(i,j,k) + tv%T(i+1,j,k))
            Shtot(I) = Shtot(I) + hweight * 0.5 * (tv%S(i,j,k) + tv%S(i+1,j,k))
          endif
        enddo ; endif

        if ((.not.CS%linear_drag) .and. (hwtot > 0.0)) then
          ustar(I) = cdrag_sqrt*hutot/hwtot
        else
          ustar(I) = cdrag_sqrt*CS%drag_bg_vel
        endif

        if (use_EOS) then ; if (hwtot > 0.0) then
          T_EOS(I) = Thtot(I)/hwtot ; S_EOS(I) = Shtot(I)/hwtot
        else
          T_EOS(I) = 0.0 ; S_EOS(I) = 0.0
        endif ; endif
      endif ; enddo ! I-loop

      if (use_EOS) then
        call calculate_density_derivs(T_EOS, S_EOS, fluxes%p_surf(:,j), &
                     dR_dT, dR_dS, Isq-G%IsdB+1, Ieq-Isq+1, tv%eqn_of_state)
      endif

      do I=Isq,Ieq ; if (do_i(I)) then
  !  The 400.0 in this expression is the square of a constant proposed
  !  by Killworth and Edwards, 1999, in equation (2.20).
        ustarsq = Rho0x400_G * ustar(i)**2
        htot(i) = 0.0
        if (use_EOS) then
          Thtot(i) = 0.0 ; Shtot(i) = 0.0
          do k=1,nz-1
            if (h_at_vel(i,k) <= 0.0) cycle
            T_Lay = 0.5 * (tv%T(i,j,k) + tv%T(i+1,j,k))
            S_Lay = 0.5 * (tv%S(i,j,k) + tv%S(i+1,j,k))
            oldfn = dR_dT(i)*(T_Lay*htot(i) - Thtot(i)) + dR_dS(i)*(S_Lay*htot(i) - Shtot(i))
            if (oldfn >= ustarsq) exit

            Dfn = (dR_dT(i)*(0.5*(tv%T(i,j,k+1)+tv%T(i+1,j,k+1)) - T_Lay) + &
                   dR_dS(i)*(0.5*(tv%S(i,j,k+1)+tv%S(i+1,j,k+1)) - S_Lay)) * &
                  (h_at_vel(i,k)+htot(i))
            if ((oldfn + Dfn) <= ustarsq) then
              Dh = h_at_vel(i,k)
            else
              Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn)/Dfn)
            endif

            htot(i) = htot(i) + Dh
            Thtot(i) = Thtot(i) + T_Lay*Dh ; Shtot(i) = Shtot(i) + S_Lay*Dh
          enddo
          if ((oldfn < ustarsq) .and. (h_at_vel(i,nz) > 0.0)) then
            T_Lay = 0.5*(tv%T(i,j,nz) + tv%T(i+1,j,nz))
            S_Lay = 0.5*(tv%S(i,j,nz) + tv%S(i+1,j,nz))
            if (dR_dT(i)*(T_Lay*htot(i) - Thtot(i)) + &
                dR_dS(i)*(S_Lay*htot(i) - Shtot(i)) < ustarsq) &
              htot(i) = htot(i) + h_at_vel(i,nz)
          endif ! Examination of layer nz.
        else  ! Use Rlay as the density variable.
          Rhtot = 0.0
          do k=1,nz-1
            Rlay = GV%Rlay(k) ; Rlb = GV%Rlay(k+1)

            oldfn = Rlay*htot(i) - Rhtot(i)
            if (oldfn >= ustarsq) exit

            Dfn = (Rlb - Rlay)*(h_at_vel(i,k)+htot(i))
            if ((oldfn + Dfn) <= ustarsq) then
              Dh = h_at_vel(i,k)
            else
              Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn)/Dfn)
            endif

            htot(i) = htot(i) + Dh
            Rhtot(i) = Rhtot(i) + Rlay*Dh
          enddo
          if (GV%Rlay(nz)*htot(i) - Rhtot(i) < ustarsq) &
            htot(i) = htot(i) + h_at_vel(i,nz)
        endif ! use_EOS

       !visc%tbl_thick_shelf_u(I,j) = max(CS%Htbl_shelf_min, &
       !    htot(I) / (0.5 + sqrt(0.25 + &
       !                 (htot(i)*(G%CoriolisBu(I,J-1)+G%CoriolisBu(I,J)))**2 / &
       !                 (ustar(i)*m_to_H)**2 )) )
        ustar1 = ustar(i)*m_to_H
        h2f2 = (htot(i)*(G%CoriolisBu(I,J-1)+G%CoriolisBu(I,J)) + h_neglect*CS%Omega)**2
        visc%tbl_thick_shelf_u(I,j) = max(CS%Htbl_shelf_min, &
            ( htot(I)*ustar1 ) / ( 0.5*ustar1 + sqrt((0.5*ustar1)**2 + h2f2 ) ) )
        visc%kv_tbl_shelf_u(I,j) = max(CS%KV_TBL_min, &
                       cdrag_sqrt*ustar(I)*visc%tbl_thick_shelf_u(I,j))
      endif ; enddo ! I-loop
    endif ! do_any_shelf

  enddo ! j-loop at u-points

!$OMP parallel do default(none) shared(u, v, h, tv, fluxes, visc, dt, G, GV, CS, use_EOS,&
!$OMP                                  dt_Rho0, h_neglect, h_tiny, g_H_Rho0,is,ie,       &
!$OMP                                  Jsq,Jeq,nz,U_bg_sq,cdrag_sqrt,Rho0x400_G,nkml,    &
!$OMP                                  m_to_H,H_to_m) &
!$OMP                          private(do_any,htot,do_i,k_massive,Thtot,vhtot,uhtot,absf,&
!$OMP                                  U_Star,Idecay_len_TKE,press,k2,I_2hlay,T_EOS,     &
!$OMP                                  S_EOS,dR_dT, dR_dS,hlay,u_at_v,Uh2,               &
!$OMP                                  T_lay,S_lay,gHprime,RiBulk,do_any_shelf,          &
!$OMP                                  Shtot,Rhtot,ustar,h_at_vel,htot_vel,hwtot,        &
!$OMP                                  hutot,hweight,ustarsq,oldfn,Dh,Rlay,Rlb,Dfn,      &
!$OMP                                  h2f2,ustar1)
  do J=Jsq,Jeq  ! v-point loop
    if (CS%dynamic_viscous_ML) then
      do_any = .false.
      do i=is,ie
        htot(i) = 0.0
        if (G%mask2dCv(i,J) < 0.5) then
          do_i(i) = .false. ; visc%nkml_visc_v(i,J) = nkml
        else
          do_i(i) = .true. ; do_any = .true.
          k_massive(i) = nkml
          Thtot(i) = 0.0 ; Shtot(i) = 0.0 ; Rhtot(i) = 0.0
          vhtot(i) = dt_Rho0 * fluxes%tauy(i,J)
          uhtot(i) = 0.25 * dt_Rho0 * ((fluxes%taux(I,j) + fluxes%tauy(I-1,j+1)) + &
                                       (fluxes%taux(I-1,j) + fluxes%tauy(I,j+1)))

         if (CS%omega_frac >= 1.0) then ; absf = 2.0*CS%omega ; else
           absf = 0.5*(abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))
           if (CS%omega_frac > 0.0) &
             absf = sqrt(CS%omega_frac*4.0*CS%omega**2 + (1.0-CS%omega_frac)*absf**2)
         endif

         U_Star = max(CS%ustar_min, 0.5 * (fluxes%ustar(i,j) + fluxes%ustar(i,j+1)))
         Idecay_len_TKE(i) = ((absf / U_Star) * CS%TKE_decay) * H_to_m

        endif
      enddo

      if (do_any) then ; do k=1,nz
        if (k > nkml) then
          do_any = .false.
          if (use_EOS .and. (k==nkml+1)) then
            ! Find dRho/dT and dRho_dS.
            do i=is,ie
              press(i) = GV%H_to_Pa * htot(i)
              k2 = max(1,nkml)
              I_2hlay = 1.0 / (h(i,j,k2) + h(i,j+1,k2) + h_neglect)
              T_EOS(i) = (h(i,j,k2)*tv%T(i,j,k2) + h(i,j+1,k2)*tv%T(i,j+1,k2)) * I_2hlay
              S_EOS(i) = (h(i,j,k2)*tv%S(i,j,k2) + h(i,j+1,k2)*tv%S(i,j+1,k2)) * I_2hlay
            enddo
            call calculate_density_derivs(T_EOS, S_EOS, press, dR_dT, dR_dS, &
                                          is-G%IsdB+1, ie-is+1, tv%eqn_of_state)
          endif

          do i=is,ie ; if (do_i(i)) then

            hlay = 0.5*(h(i,j,k) + h(i,j+1,k))
            if (hlay > h_tiny) then ! Only consider non-vanished layers.
              I_2hlay = 1.0 / (h(i,j,k) + h(i,j+1,k))
              u_at_v = 0.5 * (h(i,j,k)   * (u(I-1,j,k)   + u(I,j,k)) + &
                              h(i,j+1,k) * (u(I-1,j+1,k) + u(I,j+1,k))) * I_2hlay
              Uh2 = (uhtot(I) - htot(I)*u_at_v)**2 + (vhtot(I) - htot(I)*v(i,J,k))**2

              if (use_EOS) then
                T_lay = (h(i,j,k)*tv%T(i,j,k) + h(i,j+1,k)*tv%T(i,j+1,k)) * I_2hlay
                S_lay = (h(i,j,k)*tv%S(i,j,k) + h(i,j+1,k)*tv%S(i,j+1,k)) * I_2hlay
                gHprime = g_H_Rho0 * (dR_dT(i) * (T_lay*htot(i) - Thtot(i)) + &
                                      dR_dS(i) * (S_lay*htot(i) - Shtot(i)))
              else
                gHprime = g_H_Rho0 * (GV%Rlay(k)*htot(i) - Rhtot(i))
              endif

              if (gHprime > 0.0) then
                RiBulk = CS%bulk_Ri_ML * exp(-htot(i) * Idecay_len_TKE(i))
                if (RiBulk * Uh2 <= htot(i)**2 * gHprime) then
                  visc%nkml_visc_v(i,J) = real(k_massive(i))
                  do_i(i) = .false.
                elseif (RiBulk * Uh2 <= (htot(i) + hlay)**2 * gHprime) then
                  visc%nkml_visc_v(i,J) = real(k-1) + &
                    ( sqrt(RiBulk * Uh2 / gHprime) - htot(i) ) / hlay
                  do_i(i) = .false.
                endif
              endif
              k_massive(i) = k
            endif ! hlay > h_tiny

            if (do_i(i)) do_any = .true.
          endif ; enddo

          if (.not.do_any) exit ! All columns are done.
        endif

        do i=is,ie ; if (do_i(i)) then
          htot(i) = htot(i) + 0.5 * (h(i,J,k) + h(i,j+1,k))
          vhtot(i) = vhtot(i) + 0.5 * (h(i,j,k) + h(i,j+1,k)) * v(i,J,k)
          uhtot(i) = uhtot(i) + 0.25 * (h(i,j,k) * (u(I-1,j,k) + u(I,j,k)) + &
                                        h(i,j+1,k) * (u(I-1,j+1,k) + u(I,j+1,k)))
          if (use_EOS) then
            Thtot(i) = Thtot(i) + 0.5 * (h(i,j,k)*tv%T(i,j,k) + h(i,j+1,k)*tv%T(i,j+1,k))
            Shtot(i) = Shtot(i) + 0.5 * (h(i,j,k)*tv%S(i,j,k) + h(i,j+1,k)*tv%S(i,j+1,k))
          else
            Rhtot(i) = Rhtot(i) + 0.5 * (h(i,j,k) + h(i,j+1,k)) * GV%Rlay(k)
          endif
        endif ; enddo
      enddo ; endif

      if (do_any) then ; do i=is,ie ; if (do_i(i)) then
        visc%nkml_visc_v(i,J) = k_massive(i)
      endif ; enddo ; endif
    endif ! dynamic_viscous_ML

    do_any_shelf = .false.
    if (associated(fluxes%frac_shelf_h)) then
      do I=Is,Ie
        if (fluxes%frac_shelf_v(i,J)*G%mask2dCv(i,J) == 0.0) then
          do_i(I) = .false.
          visc%tbl_thick_shelf_v(i,J) = 0.0 ; visc%kv_tbl_shelf_v(i,J) = 0.0
        else
          do_i(I) = .true. ; do_any_shelf = .true.
        endif
      enddo
    endif

    if (do_any_shelf) then
      do k=1,nz ; do i=is,ie ; if (do_i(i)) then
        if (v(i,J,k) * (h(i,j+1,k) - h(i,j,k)) >= 0) then
          h_at_vel(i,k) = 2.0*h(i,j,k)*h(i,j+1,k) / &
                          (h(i,j,k) + h(i,j+1,k) + h_neglect)
        else
          h_at_vel(i,k) =  0.5 * (h(i,j,k) + h(i,j+1,k))
        endif
      else
        h_at_vel(I,k) = 0.0 ; ustar(i) = 0.0
      endif ; enddo ; enddo

      do i=is,ie ; if (do_i(i)) then
        htot_vel = 0.0 ; hwtot = 0.0 ; hutot = 0.0
        Thtot(i) = 0.0 ; Shtot(i) = 0.0
        if (use_EOS .or. .not.CS%linear_drag) then ; do k=1,nz
          if (htot_vel>=CS%Htbl_shelf) exit ! terminate the k loop
          hweight = MIN(CS%Htbl_shelf - htot_vel, h_at_vel(i,k))
          if (hweight <= 1.5*GV%Angstrom + h_neglect) cycle

          htot_vel  = htot_vel + h_at_vel(i,k)
          hwtot = hwtot + hweight

          if (.not.CS%linear_drag) then
            u_at_v = set_u_at_v(u, h, G, i, J, k)
            hutot = hutot + hweight * sqrt(v(i,J,k)**2 + &
                                           u_at_v**2 + U_bg_sq)
          endif
          if (use_EOS) then
            Thtot(i) = Thtot(i) + hweight * 0.5 * (tv%T(i,j,k) + tv%T(i,j+1,k))
            Shtot(i) = Shtot(i) + hweight * 0.5 * (tv%S(i,j,k) + tv%S(i,j+1,k))
          endif
        enddo ; endif

        if (.not.CS%linear_drag) then ; if (hwtot > 0.0) then
          ustar(i) = cdrag_sqrt*hutot/hwtot
        else
          ustar(i) = cdrag_sqrt*CS%drag_bg_vel
        endif ; endif

        if (use_EOS) then ; if (hwtot > 0.0) then
          T_EOS(i) = Thtot(i)/hwtot ; S_EOS(i) = Shtot(i)/hwtot
        else
          T_EOS(i) = 0.0 ; S_EOS(i) = 0.0
        endif ; endif
      endif ; enddo ! I-loop

      if (use_EOS) then
        call calculate_density_derivs(T_EOS, S_EOS, fluxes%p_surf(:,j), &
                     dR_dT, dR_dS, is-G%IsdB+1, ie-is+1, tv%eqn_of_state)
      endif

      do i=is,ie ; if (do_i(i)) then
  !  The 400.0 in this expression is the square of a constant proposed
  !  by Killworth and Edwards, 1999, in equation (2.20).
        ustarsq = Rho0x400_G * ustar(i)**2
        htot(i) = 0.0
        if (use_EOS) then
          Thtot(i) = 0.0 ; Shtot(i) = 0.0
          do k=1,nz-1
            if (h_at_vel(i,k) <= 0.0) cycle
            T_Lay = 0.5 * (tv%T(i,j,k) + tv%T(i,j+1,k))
            S_Lay = 0.5 * (tv%S(i,j,k) + tv%S(i,j+1,k))
            oldfn = dR_dT(i)*(T_Lay*htot(i) - Thtot(i)) + dR_dS(i)*(S_Lay*htot(i) - Shtot(i))
            if (oldfn >= ustarsq) exit

            Dfn = (dR_dT(i)*(0.5*(tv%T(i,j,k+1)+tv%T(i,j+1,k+1)) - T_Lay) + &
                   dR_dS(i)*(0.5*(tv%S(i,j,k+1)+tv%S(i,j+1,k+1)) - S_Lay)) * &
                  (h_at_vel(i,k)+htot(i))
            if ((oldfn + Dfn) <= ustarsq) then
              Dh = h_at_vel(i,k)
            else
              Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn)/Dfn)
            endif

            htot(i) = htot(i) + Dh
            Thtot(i) = Thtot(i) + T_Lay*Dh ; Shtot(i) = Shtot(i) + S_Lay*Dh
          enddo
          if ((oldfn < ustarsq) .and. (h_at_vel(i,nz) > 0.0)) then
            T_Lay = 0.5*(tv%T(i,j,nz) + tv%T(i,j+1,nz))
            S_Lay = 0.5*(tv%S(i,j,nz) + tv%S(i,j+1,nz))
            if (dR_dT(i)*(T_Lay*htot(i) - Thtot(i)) + &
                dR_dS(i)*(S_Lay*htot(i) - Shtot(i)) < ustarsq) &
              htot(i) = htot(i) + h_at_vel(i,nz)
          endif ! Examination of layer nz.
        else  ! Use Rlay as the density variable.
          Rhtot = 0.0
          do k=1,nz-1
            Rlay = GV%Rlay(k) ; Rlb = GV%Rlay(k+1)

            oldfn = Rlay*htot(i) - Rhtot(i)
            if (oldfn >= ustarsq) exit

            Dfn = (Rlb - Rlay)*(h_at_vel(i,k)+htot(i))
            if ((oldfn + Dfn) <= ustarsq) then
              Dh = h_at_vel(i,k)
            else
              Dh = h_at_vel(i,k) * sqrt((ustarsq-oldfn)/Dfn)
            endif

            htot(i) = htot(i) + Dh
            Rhtot = Rhtot + Rlay*Dh
          enddo
          if (GV%Rlay(nz)*htot(i) - Rhtot(i) < ustarsq) &
            htot(i) = htot(i) + h_at_vel(i,nz)
        endif ! use_EOS

       !visc%tbl_thick_shelf_v(i,J) = max(CS%Htbl_shelf_min, &
       !    htot(i) / (0.5 + sqrt(0.25 + &
       !        (htot(i)*(G%CoriolisBu(I-1,J)+G%CoriolisBu(I,J)))**2 / &
       !        (ustar(i)*m_to_H)**2 )) )
        ustar1 = ustar(i)*m_to_H
        h2f2 = (htot(i)*(G%CoriolisBu(I-1,J)+G%CoriolisBu(I,J)) + h_neglect*CS%Omega)**2
        visc%tbl_thick_shelf_v(i,J) = max(CS%Htbl_shelf_min, &
            ( htot(i)*ustar1 ) / ( 0.5*ustar1 + sqrt((0.5*ustar1)**2 + h2f2 ) ) )
        visc%kv_tbl_shelf_v(i,J) = max(CS%KV_TBL_min, &
                       cdrag_sqrt*ustar(i)*visc%tbl_thick_shelf_v(i,J))
      endif ; enddo ! i-loop
    endif ! do_any_shelf

  enddo ! J-loop at v-points

  if (CS%debug) then
    if (associated(visc%nkml_visc_u)) &
      call uchksum(visc%nkml_visc_u,"nkml_visc_u",G%HI,haloshift=0)
    if (associated(visc%nkml_visc_v)) &
      call vchksum(visc%nkml_visc_v,"nkml_visc_v",G%HI,haloshift=0)
  endif
  if (CS%id_nkml_visc_u > 0) &
    call post_data(CS%id_nkml_visc_u, visc%nkml_visc_u, CS%diag)
  if (CS%id_nkml_visc_v > 0) &
    call post_data(CS%id_nkml_visc_v, visc%nkml_visc_v, CS%diag)

end subroutine set_viscous_ML

subroutine set_visc_register_restarts(HI, GV, param_file, visc, restart_CS)
  type(hor_index_type),    intent(in)    :: HI
  type(verticalGrid_type), intent(in)    :: GV
  type(param_file_type),   intent(in)    :: param_file
  type(vertvisc_type),     intent(inout) :: visc
  type(MOM_restart_CS),    pointer       :: restart_CS
!   This subroutine is used to register any fields associated with the
! vertvisc_type.
! Arguments: HI - A horizontal index type structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (out)     visc - A structure containing vertical viscosities and related
!                   fields.  Allocated here.
!  (in)      restart_CS - A pointer to the restart control structure.
  type(vardesc) :: vd
  logical :: use_kappa_shear, adiabatic, useKPP, useEPBL
  logical :: use_CVMix, MLE_use_PBL_MLD
  integer :: isd, ied, jsd, jed, nz
  character(len=40)  :: mod = "MOM_set_visc"  ! This module's name.
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  call get_param(param_file, mod, "ADIABATIC", adiabatic, default=.false., &
                 do_not_log=.true.)
  use_kappa_shear = .false. ; use_CVMix = .false. ;
  useKPP = .false. ; useEPBL = .false.
  if (.not.adiabatic) then
    use_kappa_shear = kappa_shear_is_used(param_file)
    use_CVMix = CVMix_shear_is_used(param_file)
    call get_param(param_file, mod, "USE_KPP", useKPP, &
                 "If true, turns on the [CVmix] KPP scheme of Large et al., 1984,\n"// &
                 "to calculate diffusivities and non-local transport in the OBL.", &
                 default=.false., do_not_log=.true.)
    call get_param(param_file, mod, "ENERGETICS_SFC_PBL", useEPBL, &
                 "If true, use an implied energetics planetary boundary \n"//&
                 "layer scheme to determine the diffusivity and viscosity \n"//&
                 "in the surface boundary layer.", default=.false., do_not_log=.true.)
  endif

  if (use_kappa_shear .or. useKPP .or. useEPBL .or. use_CVMix) then
    allocate(visc%Kd_turb(isd:ied,jsd:jed,nz+1)) ; visc%Kd_turb(:,:,:) = 0.0
    allocate(visc%TKE_turb(isd:ied,jsd:jed,nz+1)) ; visc%TKE_turb(:,:,:) = 0.0
    allocate(visc%Kv_turb(isd:ied,jsd:jed,nz+1)) ; visc%Kv_turb(:,:,:) = 0.0

    vd = var_desc("Kd_turb","m2 s-1","Turbulent diffusivity at interfaces", &
                  hor_grid='h', z_grid='i')
    call register_restart_field(visc%Kd_turb, vd, .false., restart_CS)

    vd = var_desc("TKE_turb","m2 s-2","Turbulent kinetic energy per unit mass at interfaces", &
                  hor_grid='h', z_grid='i')
    call register_restart_field(visc%TKE_turb, vd, .false., restart_CS)
    vd = var_desc("Kv_turb","m2 s-1","Turbulent viscosity at interfaces", &
                  hor_grid='h', z_grid='i')
    call register_restart_field(visc%Kv_turb, vd, .false., restart_CS)
  endif

  ! visc%MLD is used to communicate the state of the (e)PBL to the rest of the model
  call get_param(param_file, mod, "MLE_USE_PBL_MLD", MLE_use_PBL_MLD, &
                 default=.false., do_not_log=.true.)
  if (MLE_use_PBL_MLD) then
    allocate(visc%MLD(isd:ied,jsd:jed)) ; visc%MLD(:,:) = 0.0
    vd = var_desc("MLD","m","Instantaneous active mixing layer depth", &
                  hor_grid='h', z_grid='1')
    call register_restart_field(visc%MLD, vd, .false., restart_CS)
  endif

end subroutine set_visc_register_restarts

subroutine set_visc_init(Time, G, GV, param_file, diag, visc, CS)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(verticalGrid_type), intent(in)    :: GV
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ctrl), target, intent(inout) :: diag
  type(vertvisc_type),     intent(inout) :: visc
  type(set_visc_CS),       pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (out)     visc - A structure containing vertical viscosities and related
!                   fields.  Allocated here.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
  real    :: Csmag_chan_dflt, smag_const1, TKE_decay_dflt, bulk_Ri_ML_dflt
  real    :: Kv_background
  real    :: omega_frac_dflt
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  logical :: use_kappa_shear, adiabatic, differential_diffusion, use_omega
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_set_visc"  ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "set_visc_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%diag => diag

! Set default, read and log parameters
  call log_version(param_file, mod, version, "")
  CS%RiNo_mix = .false.
  use_kappa_shear = .false. ; differential_diffusion = .false. !; adiabatic = .false.  ! Needed? -AJA
  call get_param(param_file, mod, "BOTTOMDRAGLAW", CS%bottomdraglaw, &
                 "If true, the bottom stress is calculated with a drag \n"//&
                 "law of the form c_drag*|u|*u. The velocity magnitude \n"//&
                 "may be an assumed value or it may be based on the \n"//&
                 "actual velocity in the bottommost HBBL, depending on \n"//&
                 "LINEAR_DRAG.", default=.true.)
  call get_param(param_file, mod, "CHANNEL_DRAG", CS%Channel_drag, &
                 "If true, the bottom drag is exerted directly on each \n"//&
                 "layer proportional to the fraction of the bottom it \n"//&
                 "overlies.", default=.false.)
  call get_param(param_file, mod, "LINEAR_DRAG", CS%linear_drag, &
                 "If LINEAR_DRAG and BOTTOMDRAGLAW are defined the drag \n"//&
                 "law is cdrag*DRAG_BG_VEL*u.", default=.false.)
  call get_param(param_file, mod, "ADIABATIC", adiabatic, default=.false., &
                 do_not_log=.true.)
  if (adiabatic) then
    call log_param(param_file, mod, "ADIABATIC",adiabatic, &
                 "There are no diapycnal mass fluxes if ADIABATIC is \n"//&
                 "true. This assumes that KD = KDML = 0.0 and that \n"//&
                 "there is no buoyancy forcing, but makes the model \n"//&
                 "faster by eliminating subroutine calls.", default=.false.)
  endif

  if (.not.adiabatic) then
    use_kappa_shear = kappa_shear_is_used(param_file)
    CS%RiNo_mix = use_kappa_shear
    call get_param(param_file, mod, "DOUBLE_DIFFUSION", differential_diffusion, &
                 "If true, increase diffusivitives for temperature or salt \n"//&
                 "based on double-diffusive paramaterization from MOM4/KPP.", &
                 default=.false.)
  endif
  call get_param(param_file, mod, "PRANDTL_TURB", visc%Prandtl_turb, &
                 "The turbulent Prandtl number applied to shear \n"//&
                 "instability.", units="nondim", default=1.0)
  call get_param(param_file, mod, "DEBUG", CS%debug, default=.false.)

  call get_param(param_file, mod, "DYNAMIC_VISCOUS_ML", CS%dynamic_viscous_ML, &
                 "If true, use a bulk Richardson number criterion to \n"//&
                 "determine the mixed layer thickness for viscosity.", &
                 default=.false.)
  if (CS%dynamic_viscous_ML) then
    call get_param(param_file, mod, "BULK_RI_ML", bulk_Ri_ML_dflt, default=0.0)
    call get_param(param_file, mod, "BULK_RI_ML_VISC", CS%bulk_Ri_ML, &
                 "The efficiency with which mean kinetic energy released \n"//&
                 "by mechanically forced entrainment of the mixed layer \n"//&
                 "is converted to turbulent kinetic energy.  By default, \n"//&
                 "BULK_RI_ML_VISC = BULK_RI_ML or 0.", units="nondim", &
                 default=bulk_Ri_ML_dflt)
    call get_param(param_file, mod, "TKE_DECAY", TKE_decay_dflt, default=0.0)
    call get_param(param_file, mod, "TKE_DECAY_VISC", CS%TKE_decay, &
                 "TKE_DECAY_VISC relates the vertical rate of decay of \n"//&
                 "the TKE available for mechanical entrainment to the \n"//&
                 "natural Ekman depth for use in calculating the dynamic \n"//&
                 "mixed layer viscosity.  By default, \n"//&
                 "TKE_DECAY_VISC = TKE_DECAY or 0.", units="nondim", &
                 default=TKE_decay_dflt)
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
    call get_param(param_file, mod, "OMEGA", CS%omega, &
                 "The rotation rate of the earth.", units="s-1", &
                 default=7.2921e-5)
    ! This give a minimum decay scale that is typically much less than Angstrom.
    CS%ustar_min = 2e-4*CS%omega*(GV%Angstrom_z + GV%H_to_m*GV%H_subroundoff)
  else
    call get_param(param_file, mod, "OMEGA", CS%omega, &
                 "The rotation rate of the earth.", units="s-1", &
                 default=7.2921e-5)
  endif

  call get_param(param_file, mod, "HBBL", CS%Hbbl, &
                 "The thickness of a bottom boundary layer with a \n"//&
                 "viscosity of KVBBL if BOTTOMDRAGLAW is not defined, or \n"//&
                 "the thickness over which near-bottom velocities are \n"//&
                 "averaged for the drag law if BOTTOMDRAGLAW is defined \n"//&
                 "but LINEAR_DRAG is not.", units="m", fail_if_missing=.true.)
  if (CS%bottomdraglaw) then
    call get_param(param_file, mod, "CDRAG", CS%cdrag, &
                 "CDRAG is the drag coefficient relating the magnitude of \n"//&
                 "the velocity field to the bottom stress. CDRAG is only \n"//&
                 "used if BOTTOMDRAGLAW is defined.", units="nondim", &
                 default=0.003)
    call get_param(param_file, mod, "DRAG_BG_VEL", CS%drag_bg_vel, &
                 "DRAG_BG_VEL is either the assumed bottom velocity (with \n"//&
                 "LINEAR_DRAG) or an unresolved  velocity that is \n"//&
                 "combined with the resolved velocity to estimate the \n"//&
                 "velocity magnitude.  DRAG_BG_VEL is only used when \n"//&
                 "BOTTOMDRAGLAW is defined.", units="m s-1", default=0.0)
    call get_param(param_file, mod, "BBL_USE_EOS", CS%BBL_use_EOS, &
                 "If true, use the equation of state in determining the \n"//&
                 "properties of the bottom boundary layer.  Otherwise use \n"//&
                 "the layer target potential densities.", default=.false.)
  endif
  call get_param(param_file, mod, "BBL_THICK_MIN", CS%BBL_thick_min, &
                 "The minimum bottom boundary layer thickness that can be \n"//&
                 "used with BOTTOMDRAGLAW. This might be \n"//&
                 "Kv / (cdrag * drag_bg_vel) to give Kv as the minimum \n"//&
                 "near-bottom viscosity.", units="m", default=0.0)
  call get_param(param_file, mod, "HTBL_SHELF_MIN", CS%Htbl_shelf_min, &
                 "The minimum top boundary layer thickness that can be \n"//&
                 "used with BOTTOMDRAGLAW. This might be \n"//&
                 "Kv / (cdrag * drag_bg_vel) to give Kv as the minimum \n"//&
                 "near-top viscosity.", units="m", default=CS%BBL_thick_min)
  call get_param(param_file, mod, "HTBL_SHELF", CS%Htbl_shelf, &
                 "The thickness over which near-surface velocities are \n"//&
                 "averaged for the drag law under an ice shelf.  By \n"//&
                 "default this is the same as HBBL", units="m", default=CS%Hbbl)

  call get_param(param_file, mod, "KV", Kv_background, &
                 "The background kinematic viscosity in the interior. \n"//&
                 "The molecular value, ~1e-6 m2 s-1, may be used.", &
                 units="m2 s-1", fail_if_missing=.true.)
  call get_param(param_file, mod, "KV_BBL_MIN", CS%KV_BBL_min, &
                 "The minimum viscosities in the bottom boundary layer.", &
                 units="m2 s-1", default=Kv_background)
  call get_param(param_file, mod, "KV_TBL_MIN", CS%KV_TBL_min, &
                 "The minimum viscosities in the top boundary layer.", &
                 units="m2 s-1", default=Kv_background)

  if (CS%Channel_drag) then
    call get_param(param_file, mod, "SMAG_LAP_CONST", smag_const1, default=-1.0)

    cSmag_chan_dflt = 0.15
    if (smag_const1 >= 0.0) cSmag_chan_dflt = smag_const1

    call get_param(param_file, mod, "SMAG_CONST_CHANNEL", CS%c_Smag, &
                 "The nondimensional Laplacian Smagorinsky constant used \n"//&
                 "in calculating the channel drag if it is enabled.  The \n"//&
                 "default is to use the same value as SMAG_LAP_CONST if \n"//&
                 "it is defined, or 0.15 if it is not. The value used is \n"//&
                 "also 0.15 if the specified value is negative.", &
                 units="nondim", default=cSmag_chan_dflt)
    if (CS%c_Smag < 0.0) CS%c_Smag = 0.15
  endif

  if (CS%bottomdraglaw) then
    allocate(visc%bbl_thick_u(IsdB:IedB,jsd:jed)) ; visc%bbl_thick_u = 0.0
    allocate(visc%kv_bbl_u(IsdB:IedB,jsd:jed)) ; visc%kv_bbl_u = 0.0
    allocate(visc%bbl_thick_v(isd:ied,JsdB:JedB)) ; visc%bbl_thick_v = 0.0
    allocate(visc%kv_bbl_v(isd:ied,JsdB:JedB)) ; visc%kv_bbl_v = 0.0
    allocate(visc%ustar_bbl(isd:ied,jsd:jed)) ; visc%ustar_bbl = 0.0
    allocate(visc%TKE_bbl(isd:ied,jsd:jed)) ; visc%TKE_bbl = 0.0

    CS%id_bbl_thick_u = register_diag_field('ocean_model', 'bbl_thick_u', &
       diag%axesCu1, Time, 'BBL thickness at u points', 'meter')
    CS%id_kv_bbl_u = register_diag_field('ocean_model', 'kv_bbl_u', diag%axesCu1, &
       Time, 'BBL viscosity at u points', 'meter2 second-1')
    CS%id_bbl_thick_v = register_diag_field('ocean_model', 'bbl_thick_v', &
       diag%axesCv1, Time, 'BBL thickness at v points', 'meter')
    CS%id_kv_bbl_v = register_diag_field('ocean_model', 'kv_bbl_v', diag%axesCv1, &
       Time, 'BBL viscosity at v points', 'meter2 second-1')
  endif
  if (CS%Channel_drag) then
    allocate(visc%Ray_u(IsdB:IedB,jsd:jed,nz)) ; visc%Ray_u = 0.0
    allocate(visc%Ray_v(isd:ied,JsdB:JedB,nz)) ; visc%Ray_v = 0.0
    CS%id_Ray_u = register_diag_field('ocean_model', 'Rayleigh_u', diag%axesCuL, &
       Time, 'Rayleigh drag velocity at u points', 'meter second-1')
    CS%id_Ray_v = register_diag_field('ocean_model', 'Rayleigh_v', diag%axesCvL, &
       Time, 'Rayleigh drag velocity at v points', 'meter second-1')
  endif

  if (differential_diffusion) then
    allocate(visc%Kd_extra_T(isd:ied,jsd:jed,nz+1)) ; visc%Kd_extra_T = 0.0
    allocate(visc%Kd_extra_S(isd:ied,jsd:jed,nz+1)) ; visc%Kd_extra_S = 0.0
  endif

  if (CS%dynamic_viscous_ML) then
    allocate(visc%nkml_visc_u(IsdB:IedB,jsd:jed)) ; visc%nkml_visc_u = 0.0
    allocate(visc%nkml_visc_v(isd:ied,JsdB:JedB)) ; visc%nkml_visc_v = 0.0
    CS%id_nkml_visc_u = register_diag_field('ocean_model', 'nkml_visc_u', &
       diag%axesCu1, Time, 'Number of layers in viscous mixed layer at u points', 'meter')
    CS%id_nkml_visc_v = register_diag_field('ocean_model', 'nkml_visc_v', &
       diag%axesCv1, Time, 'Number of layers in viscous mixed layer at v points', 'meter')
  endif

  CS%Hbbl = CS%Hbbl * GV%m_to_H
  CS%BBL_thick_min = CS%BBL_thick_min * GV%m_to_H

end subroutine set_visc_init

subroutine set_visc_end(visc, CS)
  type(vertvisc_type), intent(inout) :: visc
  type(set_visc_CS),   pointer       :: CS
  if (CS%bottomdraglaw) then
    deallocate(visc%bbl_thick_u) ; deallocate(visc%bbl_thick_v)
    deallocate(visc%kv_bbl_u) ; deallocate(visc%kv_bbl_v)
  endif
  if (CS%Channel_drag) then
    deallocate(visc%Ray_u) ; deallocate(visc%Ray_v)
  endif
  if (CS%dynamic_viscous_ML) then
    deallocate(visc%nkml_visc_u) ; deallocate(visc%nkml_visc_v)
  endif
  if (associated(visc%Kd_turb)) deallocate(visc%Kd_turb)
  if (associated(visc%TKE_turb)) deallocate(visc%TKE_turb)
  if (associated(visc%Kv_turb)) deallocate(visc%Kv_turb)
  if (associated(visc%ustar_bbl)) deallocate(visc%ustar_bbl)
  if (associated(visc%TKE_bbl)) deallocate(visc%TKE_bbl)
  if (associated(visc%taux_shelf)) deallocate(visc%taux_shelf)
  if (associated(visc%tauy_shelf)) deallocate(visc%tauy_shelf)
  if (associated(visc%tbl_thick_shelf_u)) deallocate(visc%tbl_thick_shelf_u)
  if (associated(visc%tbl_thick_shelf_v)) deallocate(visc%tbl_thick_shelf_v)
  if (associated(visc%kv_tbl_shelf_u)) deallocate(visc%kv_tbl_shelf_u)
  if (associated(visc%kv_tbl_shelf_v)) deallocate(visc%kv_tbl_shelf_v)

  deallocate(CS)
end subroutine set_visc_end

end module MOM_set_visc
