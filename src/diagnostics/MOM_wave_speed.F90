module MOM_wave_speed
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
!
!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, May 2008                                       *
!*                                                                     *
!*    The subroutine in this module calculates the first baroclinic    *
!*  mode internal wave speed.                                          *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, vh, vav                               *
!*    j    x ^ x ^ x   At >:  u, uh, uav                               *
!*    j    > o > o >   At o:  h                                        *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ptrs
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : read_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_EOS, only : calculate_density, calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>

public wave_speed, wave_speed_init

type, public :: wave_speed_CS ; private
  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                             ! ocean diagnostic fields and control variables.
end type wave_speed_CS

contains

subroutine wave_speed(h, tv, G, cg1, CS, full_halos)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)  :: h
  type(thermo_var_ptrs),                 intent(in)  :: tv
  real, dimension(NIMEM_,NJMEM_),        intent(out) :: cg1
  type(ocean_grid_type),                 intent(in)  :: G
  type(wave_speed_CS), optional,         pointer     :: CS
  logical,             optional,         intent(in)  :: full_halos
!    This subroutine determines the first mode internal wave speed.
! Arguments: h - Layer thickness, in m or kg m-2.
!  (in)      tv - A structure containing the thermobaric variables.
!  (out)     cg1 - The first mode internal gravity wave speed, in m s-1.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 wave_speed_init.
!  (in,opt)  full_halos - If true, do the calculation over the entire
!                         computational domain.

!   This subroutine solves for the first baroclinic mode wave speed.  (It could
! solve for all the wave speeds, but the iterative approach taken here means
! that this is not particularly efficient.)

!   If e(k) is the perturbation interface height, this means solving for the
! smallest eigenvalue (lam = 1/c^2) of the system
!
!   -Igu(k)*e(k-1) + (Igu(k)+Igl(k)-lam)*e(k) - Igl(k)*e(k+1) = 0.0
!
! with rigid lid boundary conditions e(1) = e(nz+1) = 0.0 giving
!
!   (Igu(2)+Igl(2)-lam)*e(2) - Igl(2)*e(3) = 0.0
!   -Igu(nz)*e(nz-1) + (Igu(nz)+Igl(nz)-lam)*e(nz) = 0.0
!
! Here
!   Igl(k) = 1.0/(gprime(k)*H(k)) ; Igu(k) = 1.0/(gprime(k)*H(k-1))
!
!   Alternately, these same eigenvalues can be found from the second smallest
! eigenvalue of the Montgomery potential (M(k)) calculation:
!
!   -Igl(k)*M(k-1) + (Igl(k)+Igu(k+1)-lam)*M(k) - Igu(k+1)*M(k+1) = 0.0
!
! with rigid lid and flat bottom boundary conditions
!
!   (Igu(2)-lam)*M(1) - Igu(2)*M(2) = 0.0
!   -Igl(nz)*M(nz-1) + (Igl(nz)-lam)*M(nz) = 0.0
!
! Note that the barotropic mode has been eliminated from the rigid lid
! interface height equations, hence the matrix is one row smaller.  Without
! the rigid lid, the top boundary condition is simpler to implement with
! the M equations.

  real, dimension(SZK_(G)+1) :: &
    dRho_dT, dRho_dS, &
    pres, T_int, S_int, &
    gprime        ! The reduced gravity across each interface, in m s-2.
  real, dimension(SZK_(G)) :: &
    Igl, Igu      ! The inverse of the reduced gravity across an interface times
                  ! the thickness of the layer below (Igl) or above (Igu) it,
                  ! in units of s2 m-2.
  real, dimension(SZK_(G),SZI_(G)) :: &
    Hf, Tf, Sf, Rf
  real, dimension(SZK_(G)) :: &
    Hc, Tc, Sc, Rc, &
    det, ddet
  real :: lam, dlam, lam0
  real :: min_h_frac
  real :: H_to_pres
  real, dimension(SZI_(G)) :: &
    htot, hmin, &  ! Thicknesses in m.
    H_here, HxT_here, HxS_here, HxR_here
  real :: speed2_tot
  real :: I_Hnew, drxh_sum
  real, parameter :: tol1  = 0.0001, tol2 = 0.001
  real, pointer, dimension(:,:,:) :: T, S
  real :: g_Rho0  ! G_Earth/Rho0 in m4 s-2 kg-1.
  real :: rescale, I_rescale
  integer :: kf(SZI_(G))
  integer, parameter :: max_itt = 10
  real :: lam_it(max_itt), det_it(max_itt), ddet_it(max_itt)
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  integer :: kc
  integer :: i, j, k, k2, itt, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (present(CS)) then
    if (.not. associated(CS)) call MOM_error(FATAL, "MOM_wave_speed: "// &
           "Module must be initialized before it is used.")
  endif
  if (present(full_halos)) then ; if (full_halos) then
    is = G%isd ; ie = G%ied ; js = G%jsd ; je = G%jed
  endif ; endif

  S => tv%S ; T => tv%T
  g_Rho0 = G%g_Earth/G%Rho0
  use_EOS = associated(tv%eqn_of_state)

  H_to_pres = G%g_Earth * G%Rho0
  rescale = 1024.0**4 ; I_rescale = 1.0/rescale

  min_h_frac = tol1 / real(nz)

  do j=js,je
    !   First merge very thin layers with the one above (or below if they are
    ! at the top).  This also transposes the row order so that columns can
    ! be worked upon one at a time.
    do i=is,ie ; htot(i) = 0.0 ; enddo
    do k=1,nz ; do i=is,ie ; htot(i) = htot(i) + h(i,j,k)*G%H_to_m ; enddo ; enddo

    do i=is,ie
      hmin(i) = htot(i)*min_h_frac ; kf(i) = 1 ; H_here(i) = 0.0
      HxT_here(i) = 0.0 ; HxS_here(i) = 0.0 ; HxR_here(i) = 0.0
    enddo
    if (use_EOS) then
      do k=1,nz ; do i=is,ie
        if ((H_here(i) > hmin(i)) .and. (h(i,j,k)*G%H_to_m > hmin(i))) then
          Hf(kf(i),i) = H_here(i)
          Tf(kf(i),i) = HxT_here(i) / H_here(i)
          Sf(kf(i),i) = HxS_here(i) / H_here(i)
          kf(i) = kf(i) + 1

          ! Start a new layer
          H_here(i) = h(i,j,k)*G%H_to_m
          HxT_here(i) = (h(i,j,k)*G%H_to_m)*T(i,j,k)
          HxS_here(i) = (h(i,j,k)*G%H_to_m)*S(i,j,k)
        else
          H_here(i) = H_here(i) + h(i,j,k)*G%H_to_m
          HxT_here(i) = HxT_here(i) + (h(i,j,k)*G%H_to_m)*T(i,j,k)
          HxS_here(i) = HxS_here(i) + (h(i,j,k)*G%H_to_m)*S(i,j,k)
        endif
      enddo ; enddo
      do i=is,ie ; if (H_here(i) > 0.0) then
        Hf(kf(i),i) = H_here(i)
        Tf(kf(i),i) = HxT_here(i) / H_here(i)
        Sf(kf(i),i) = HxS_here(i) / H_here(i)
      endif ; enddo
    else
      do k=1,nz ; do i=is,ie
        if ((H_here(i) > hmin(i)) .and. (h(i,j,k)*G%H_to_m > hmin(i))) then
          Hf(kf(i),i) = H_here(i) ; Rf(kf(i),i) = HxR_here(i) / H_here(i)
          kf(i) = kf(i) + 1

          ! Start a new layer
          H_here(i) = h(i,j,k)*G%H_to_m
          HxR_here(i) = (h(i,j,k)*G%H_to_m)*G%Rlay(k)
        else
          H_here(i) = H_here(i) + h(i,j,k)*G%H_to_m
          HxR_here(i) = HxR_here(i) + (h(i,j,k)*G%H_to_m)*G%Rlay(k)
        endif
      enddo ; enddo
      do i=is,ie ; if (H_here(i) > 0.0) then
        Hf(kf(i),i) = H_here(i) ; Rf(kf(i),i) = HxR_here(i) / H_here(i)
      endif ; enddo
    endif

    ! From this point, we can work on individual columns without causing memory
    ! to have page faults.
    do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
      if (use_EOS) then
        pres(1) = 0.0
        do k=2,kf(i)
          pres(k) = pres(k-1) + H_to_pres*Hf(k-1,i)
          T_int(k) = 0.5*(Tf(k,i)+Tf(k-1,i))
          S_int(k) = 0.5*(Sf(k,i)+Sf(k-1,i))
        enddo
        call calculate_density_derivs(T_int, S_int, pres, drho_dT, drho_dS, 2, &
                                      kf(i)-1, tv%eqn_of_state)

        ! Sum the reduced gravities to find out how small a density difference
        ! is negligibly small.
        drxh_sum = 0.0
        do k=2,kf(i)
          drxh_sum = drxh_sum + 0.5*(Hf(k-1,i)+Hf(k,i)) * &
              max(0.0,drho_dT(k)*(Tf(k,i)-Tf(k-1,i)) + &
                      drho_dS(k)*(Sf(k,i)-Sf(k-1,i)))
        enddo
      else
        drxh_sum = 0.0
        do k=2,kf(i)
          drxh_sum = drxh_sum + 0.5*(Hf(k-1,i)+Hf(k,i)) * &
                            max(0.0,Rf(k,i)-Rf(k-1,i))
        enddo
      endif
  !   Find gprime across each internal interface, taking care of convective
  ! instabilities by merging layers.
      if (drxh_sum <= 0.0) then
        cg1(i,j) = 0.0
      else
        ! Merge layers to eliminate convective instabilities or exceedingly
        ! small reduced gravities.
        if (use_EOS) then
          kc = 1
          Hc(1) = Hf(1,i) ; Tc(1) = Tf(1,i) ; Sc(1) = Sf(1,i)
          do k=2,kf(i)
            if ((drho_dT(k)*(Tf(k,i)-Tc(kc)) + drho_dS(k)*(Sf(k,i)-Sc(kc))) * &
                (Hc(kc) + Hf(k,i)) < 2.0 * tol2*drxh_sum) then
              ! Merge this layer with the one above and backtrack.
              I_Hnew = 1.0 / (Hc(kc) + Hf(k,i))
              Tc(kc) = (Hc(kc)*Tc(kc) + Hf(k,i)*Tf(k,i)) * I_Hnew
              Sc(kc) = (Hc(kc)*Sc(kc) + Hf(k,i)*Sf(k,i)) * I_Hnew
              Hc(kc) = (Hc(kc) + Hf(k,i))
              ! Backtrack to remove any convective instabilities above...  Note
              ! that the tolerance is a factor of two larger, to avoid limit how
              ! far back we go.
              do k2=kc,2,-1
                if ((drho_dT(k2)*(Tc(k2)-Tc(k2-1)) + drho_dS(k2)*(Sc(k2)-Sc(k2-1))) * &
                    (Hc(k2) + Hc(k2-1)) < tol2*drxh_sum) then
                  ! Merge the two bottommost layers.  At this point kc = k2.
                  I_Hnew = 1.0 / (Hc(kc) + Hc(kc-1))
                  Tc(kc-1) = (Hc(kc)*Tc(kc) + Hc(kc-1)*Tc(kc-1)) * I_Hnew
                  Sc(kc-1) = (Hc(kc)*Sc(kc) + Hc(kc-1)*Sc(kc-1)) * I_Hnew
                  Hc(kc-1) = (Hc(kc) + Hc(kc-1))
                  kc = kc - 1
                else ; exit ; endif
              enddo
            else
              ! Add a new layer to the column.
              kc = kc + 1
              drho_dS(kc) = drho_dS(k) ; drho_dT(kc) = drho_dT(k)
              Tc(kc) = Tf(k,i) ; Sc(kc) = Sf(k,i) ; Hc(kc) = Hf(k,i) 
            endif
          enddo
          ! At this point there are kc layers and the gprimes should be positive.
          do k=2,kc ! Revisit this if non-Boussinesq.
            gprime(k) = g_Rho0 * (drho_dT(k)*(Tc(k)-Tc(k-1)) + &
                                  drho_dS(k)*(Sc(k)-Sc(k-1)))
          enddo
        else  ! .not.use_EOS
          ! Do the same with density directly...
          kc = 1
          Hc(1) = Hf(1,i) ; Rc(1) = Rf(1,i)
          do k=2,kf(i)
            if ((Rf(k,i) - Rc(kc)) * (Hc(kc) + Hf(k,i)) < 2.0*tol2*drxh_sum) then
              ! Merge this layer with the one above and backtrack.
              Rc(kc) = (Hc(kc)*Rc(kc) + Hf(k,i)*Rf(k,i)) / (Hc(kc) + Hf(k,i))
              Hc(kc) = (Hc(kc) + Hf(k,i))
              ! Backtrack to remove any convective instabilities above...  Note
              ! that the tolerance is a factor of two larger, to avoid limit how
              ! far back we go.
              do k2=kc,2,-1
                if ((Rc(k2)-Rc(k2-1)) * (Hc(k2)+Hc(k2-1)) < tol2*drxh_sum) then
                  ! Merge the two bottommost layers.  At this point kc = k2.
                  Rc(kc-1) = (Hc(kc)*Rc(kc) + Hc(kc-1)*Rc(kc-1)) / (Hc(kc) + Hc(kc-1))
                  Hc(kc-1) = (Hc(kc) + Hc(kc-1))
                  kc = kc - 1
                else ; exit ; endif
              enddo
            else
              ! Add a new layer to the column.
              kc = kc + 1
              Rc(kc) = Rf(k,i) ; Hc(kc) = Hf(k,i) 
            endif
          enddo
          ! At this point there are kc layers and the gprimes should be positive.
          do k=2,kc ! Revisit this if non-Boussinesq.
            gprime(k) = g_Rho0 * (Rc(k)-Rc(k-1))
          enddo
        endif  ! use_EOS
        
        !   Sum the contributions from all of the interfaces to give an over-estimate
        ! of the first-mode wave speed.
        if (kc >= 2) then
          speed2_tot = 0.0
          do k=2,kc
            Igl(k) = 1.0/(gprime(k)*Hc(k)) ; Igu(k) = 1.0/(gprime(k)*Hc(k-1))
            speed2_tot = speed2_tot + gprime(k)*(Hc(k-1)+Hc(k))
          enddo

          ! Overestimate the speed to start with.
          lam0 = 1.0 / speed2_tot ; lam = lam0
          ! Find the determinant and its derivative with lam.
          do itt=1,max_itt
            lam_it(itt) = lam
            det(1) = 1.0 ; ddet(1) = 0.0
            det(2) = (Igu(2)+Igl(2)-lam) ; ddet(2) = -1.0
            do k=3,kc
              det(k) = (Igu(k)+Igl(k)-lam)*det(k-1) - (Igu(k)*Igl(k-1))*det(k-2)
              ddet(k) = (Igu(k)+Igl(k)-lam)*ddet(k-1) - (Igu(k)*Igl(k-1))*ddet(k-2) - &
                        det(k-1)
                        
              ! Rescale det & ddet if det is getting too large.
              if (abs(det(k)) > rescale) then
                det(k) = I_rescale*det(k) ; det(k-1) = I_rescale*det(k-1)
                ddet(k) = I_rescale*ddet(k) ; ddet(k-1) = I_rescale*ddet(k-1)
              endif
            enddo
            ! Use Newton's method iteration to find a new estimate of lam.
            det_it(itt) = det(kc) ; ddet_it(itt) = ddet(kc)

            if ((ddet(kc) >= 0.0) .or. (-det(kc) > -0.5*lam*ddet(kc))) then
              ! lam was not an under-estimate, as intended, so Newton's method
              ! may not be reliable; lam must be reduced, but not by more
              ! than half.
              lam = 0.5 * lam
            else  ! Newton's method is OK.
              dlam = - det(kc) / ddet(kc)
              lam = lam + dlam
              if (abs(dlam) < tol2*lam) exit
            endif
          enddo

          cg1(i,j) = 0.0
          if (lam > 0.0) cg1(i,j) = 1.0 / sqrt(lam)
        else
          cg1(i,j) = 0.0
        endif
       
      endif ! cg1 /= 0.0
    else
      cg1(i,j) = 0.0 ! This is a land point.
    endif ; enddo ! i-loop
  enddo ! j-loop

end subroutine wave_speed

subroutine wave_speed_init(Time, G, param_file, diag, CS)
  type(time_type),             intent(in)    :: Time
  type(ocean_grid_type),       intent(in)    :: G
  type(param_file_type),       intent(in)    :: param_file
  type(diag_ptrs), target,     intent(inout) :: diag
  type(wave_speed_CS),         pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                  for this module
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'
  character(len=40)  :: mod = "MOM_wave_speed"  ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "wave_speed_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag

  ! Write all relevant parameters to the model log.
  call log_version(param_file, mod, version, tagname, "")

end subroutine wave_speed_init

end module MOM_wave_speed
