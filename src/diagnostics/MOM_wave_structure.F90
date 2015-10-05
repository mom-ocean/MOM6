module MOM_wave_structure
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
!*  By Ben Mater, September, 2015                                      *
!*                                                                     *
!*    The subroutine in this module calculates the vertical structure  *
!*    functions of the first baroclinic mode internal wave speed.      *
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

use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : log_version, param_file_type
use MOM_grid,       only : ocean_grid_type
use MOM_variables,  only : thermo_var_ptrs

implicit none ; private

#include <MOM_memory.h>

public wave_structure, wave_structure_init

type, public :: wave_structure_CS ; !private
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                                   ! timing of diagnostic output.
  real, allocatable, dimension(:,:,:) :: w_strct
                                   ! Vertical structure of vertical velocity (normalized), in m s-1.
  real, allocatable, dimension(:,:,:) :: u_strct
                                   ! Vertical structure of horizontal velocity (normalized), in m s-1.
  real, allocatable, dimension(:,:,:) :: W_profile
                                   ! Vertical profile of w_hat(z), where 
                                   ! w(x,y,z,t) = w_hat(z)*exp(i(kx+ly-freq*t)) is the full time-
                                   ! varying vertical velocity with w_hat(z) = W0*w_strct(z), in m s-1.
  real, allocatable, dimension(:,:,:) :: Uavg_profile
                                   ! Vertical profile of the magnitude of horizontal velocity,
                                   ! (u^2+v^2)^0.5, averaged over a period, in m s-1.
  real, allocatable, dimension(:,:,:) :: z_depths
                                   ! Depths of layer interfaces, in m.
  real, allocatable, dimension(:,:,:) :: N2
                                   ! Squared buoyancy frequency at each interface
  integer, allocatable, dimension(:,:):: num_intfaces
                                   ! Number of layer interfaces (including surface and bottom)
                                   
end type wave_structure_CS

contains

subroutine wave_structure(h, tv, G, cn, freq, CS, En, full_halos)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),    intent(in)  :: h
  type(thermo_var_ptrs),                    intent(in)  :: tv
  type(ocean_grid_type),                    intent(in)  :: G
  real, dimension(NIMEM_,NJMEM_),           intent(in)  :: cn
  real,                                     intent(in)  :: freq
  type(wave_structure_CS),                  pointer     :: CS
  real, dimension(NIMEM_,NJMEM_), optional, intent(in)  :: En
  logical,optional,                         intent(in)  :: full_halos
  
!    This subroutine determines the first mode internal wave velocity structure.
! Arguments: h - Layer thickness, in m or kg m-2.
!  (in)      tv - A structure containing the thermobaric variables.
!  (in)      G - The ocean's grid structure.
!  (in)      cn - The (non-rotational) mode internal gravity wave speed, in m s-1.
!  (in)      freq - intrinsic wave frequency, in s-1
!  (in)      CS - The control structure returned by a previous call to
!                 wave_structure_init.
!  (in,opt)  En - Internal wave energy density, in Jm-2
!  (in,opt)  full_halos - If true, do the calculation over the entire
!                         computational domain.
!
! This subroutine solves for the eigen vector [vertical structure, e(k)] associated with 
! the first baroclinic mode speed [i.e., smallest eigen value (lam = 1/c^2)] of the 
! system d2e/dz2 = -(N2/cn2)e, or (A-lam*I)e = 0, where A = -(1/N2)(d2/dz2), lam = 1/c^2,
! and I is the identity matrix. 2nd order discretization in the vertical lets this system
! be represented as 
!
!   -Igu(k)*e(k-1) + (Igu(k)+Igl(k)-lam)*e(k) - Igl(k)*e(k+1) = 0.0
!
! with rigid lid boundary conditions e(1) = e(nz+1) = 0.0 giving
!
!   (Igu(2)+Igl(2)-lam)*e(2) - Igl(2)*e(3) = 0.0
!   -Igu(nz)*e(nz-1) + (Igu(nz)+Igl(nz)-lam)*e(nz) = 0.0
!
! where, upon noting N2 = reduced gravity/layer thickness, we get
!    Igl(k) = 1.0/(gprime(k)*H(k)) ; Igu(k) = 1.0/(gprime(k)*H(k-1))
!
! The eigen value for this system is approximated using "wave_speed." This subroutine uses 
! these eigen values (mode speeds) to estimate the corresponding eigen vectors (velocity
! structure) using the "inverse iteration with shift" method. The algorithm is 
! 
!   Pick a starting vector reasonably close to mode structure and with unit magnitude, b_guess
!   For n=1,2,3,...
!     Solve (A-lam*I)e = e_guess for e
!     Set e_guess=e/|e| and repeat, with each iteration refining the estimate of e

  real, dimension(SZK_(G)+1) :: &
    dRho_dT, dRho_dS, &
    pres, T_int, S_int, &
    gprime        ! The reduced gravity across each interface, in m s-2.
  real, dimension(SZK_(G)) :: &
    Igl, Igu      ! The inverse of the reduced gravity across an interface times
                  ! the thickness of the layer below (Igl) or above (Igu) it,
                  ! in units of s2 m-2.
  real, dimension(SZK_(G),SZI_(G)) :: &
    Hf, Tf, Sf, Rf, htot
  real, dimension(SZK_(G)) :: &
    Hc, Tc, Sc, Rc, &
    det, ddet
  real :: lam, dlam, lam0
  real :: min_h_frac
  real :: H_to_pres
  real, dimension(SZI_(G)) :: &
    hmin, &  ! Thicknesses in m.
    H_here, HxT_here, HxS_here, HxR_here
  real :: speed2_tot
  real :: I_Hnew, drxh_sum
  real, parameter :: tol1  = 0.0001, tol2 = 0.001
  real, pointer, dimension(:,:,:) :: T, S
  real :: g_Rho0  ! G_Earth/Rho0 in m4 s-2 kg-1.
  real :: rescale, I_rescale
  integer :: kf(SZI_(G))
  integer, parameter :: max_itt = 10
  real, parameter :: cg_subRO = 1e-100
  real, parameter :: a_int = 0.5 ! value of normalized integral: \int(w_strct^2)dz = a_int
  real            :: I_a_int     ! inverse of a_int
  real            :: f2          ! squared Coriolis frequency
  real            :: Kmag2       ! magnitude of horizontal wave number squared
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  real, allocatable, dimension(:) :: w_strct, u_strct, W_profile, Uavg_profile, z_int, N2 
  real, allocatable, dimension(:) :: w_strct2, u_strct2 ! squared values
  real, allocatable, dimension(:) :: dz      ! thicknesses of merged layers (same as Hc I hope)
  real, allocatable, dimension(:) :: dWdz_profile ! profile of dW/dz
  real                            :: w2avg   ! average of squared vertical velocity structure funtion
  real                            :: int_dwdz2, int_w2, int_N2w2, KE_term, PE_term, W0
                                             ! terms in vertically averaged energy equation
  real, allocatable, dimension(:) :: lam_z   ! product of eigen value and gprime(k)
  real, allocatable, dimension(:) :: a_diag, b_diag, c_diag
                                             ! diagonals of tridiagonal matrix
  real, allocatable, dimension(:) :: e_guess ! guess at eigen vector with unit amplitde (for TDMA)
  real, allocatable, dimension(:) :: e_itt   ! improved guess at eigen vector (from TDMA)
  real    :: Pi
  integer :: kc
  integer :: i, j, k, k2, itt, is, ie, js, je, nz, row

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  I_a_int = 1/a_int
  
  !if (present(CS)) then
    if (.not. associated(CS)) call MOM_error(FATAL, "MOM_wave_structure: "// &
           "Module must be initialized before it is used.")
  !endif
  if (present(full_halos)) then ; if (full_halos) then
    is = G%isd ; ie = G%ied ; js = G%jsd ; je = G%jed
  endif ; endif
  
  Pi = (4.0*atan(1.0))
  
  S => tv%S ; T => tv%T
  g_Rho0 = G%g_Earth/G%Rho0
  use_EOS = associated(tv%eqn_of_state)

  H_to_pres = G%g_Earth * G%Rho0
  rescale = 1024.0**4 ; I_rescale = 1.0/rescale

  min_h_frac = tol1 / real(nz)
!$OMP parallel do default(none) shared(is,ie,js,je,nz,h,G,min_h_frac,use_EOS,T,S,      &
!$OMP                                  H_to_pres,tv,cg1,g_Rho0,rescale,I_rescale)      &
!$OMP                          private(htot,hmin,kf,H_here,HxT_here,HxS_here,HxR_here, &
!$OMP                                  Hf,Tf,Sf,Rf,pres,T_int,S_int,drho_dT,           &
!$OMP                                  drho_dS,drxh_sum,kc,Hc,Tc,Sc,I_Hnew,gprime,     &
!$OMP                                  Rc,speed2_tot,Igl,Igu,lam0,lam,lam_it,dlam,     &
!$OMP                                  det,ddet,det_it,ddet_it)

  do j=js,je
    !   First merge very thin layers with the one above (or below if they are
    ! at the top).  This also transposes the row order so that columns can
    ! be worked upon one at a time.
    do i=is,ie ; htot(i,j) = 0.0 ; enddo
    do k=1,nz ; do i=is,ie ; htot(i,j) = htot(i,j) + h(i,j,k)*G%H_to_m ; enddo ; enddo

    do i=is,ie
      hmin(i) = htot(i,j)*min_h_frac ; kf(i) = 1 ; H_here(i) = 0.0
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
    endif ! use_EOS?

    ! From this point, we can work on individual columns without causing memory
    ! to have page faults.
    do i=is,ie
      if (G%mask2dT(i,j) > 0.5) then

        lam = 1/(cn(i,j)**2)
        
        ! Calculate drxh_sum
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
        endif ! use_EOS?
        
        !   Find gprime across each internal interface, taking care of convective
        ! instabilities by merging layers.
        if (drxh_sum >= 0.0) then
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
          endif  ! use_EOS?
          
          ! Construct and solve tridiagonal system
          if (kc >= 2) then
            z_int(1) = 0.0
            ! Calculate Igu, Igl, depth, and N2 at each interior interface
            do k=2,kc
              Igl(k) = 1.0/(gprime(k)*Hc(k)) ; Igu(k) = 1.0/(gprime(k)*Hc(k-1))
              z_int(k) = z_int(k-1) + Hc(k-1)
              N2(k) = gprime(k)/(0.5*(Hc(k)+Hc(k-1)))
              lam_z(k) = lam*gprime(k)
            enddo
            z_int(kc+1) = z_int(kc)+Hc(kc)
            if (abs(z_int(kc+1)-htot(i,j)) > 1.e-30) then
              call MOM_error(WARNING, "wave_structure: mismatch in total depths")
            endif

            ! Populate interior rows of tridiagonal matrix; must multiply through by
            ! gprime to get tridiagonal matrix to the symmetrical form:
            ! [-1/H(k-1)]e(k-1) + [1/H(k-1)+1/H(k)-lam_z]e(k) + [-1/H(k)]e(k+1) = 0,
            ! where lam_z = lam*gprime is now a function of depth
            do k=3,kc-1
              row = k-1
              a_diag(row) = gprime(k)*(-Igu(k))
              b_diag(row) = gprime(k)*(Igu(k)+Igl(k)) - lam_z(k)
              c_diag(row) = gprime(k)*(-Igl(k))
            enddo
            ! Populate top row of tridiagonal matrix
            k=2 ; row = k-1
            a_diag(row) = 0.0
            b_diag(row) = gprime(k)*(Igu(k)+Igl(k)) - lam_z(k)
            c_diag(row) = gprime(k)*(-Igl(k))
            ! Populate bottom row of tridiagonal matrix
            k=kc ; row = k-1
            a_diag(row) = gprime(k)*(-Igu(k))
            b_diag(row) = gprime(k)*(Igu(k)+Igl(k)) - lam_z(k)
            c_diag(row) = 0.0
            
            ! Guess a vector shape to start with
            e_guess = cos(z_int(2:kc)/htot(i,j)*Pi)
            e_guess = e_guess/sqrt(sum(e_guess**2))
            
            ! Perform inverse iteration with tri-diag solver
            do itt=1,max_itt
              call tridiag_solver(a_diag,b_diag,c_diag,lam_z,e_guess,'TDMA_T',e_itt)
              e_guess = e_itt/sqrt(sum(e_itt**2))                           
            enddo ! itt-loop
            w_strct(2:kc) = e_guess
            w_strct(1)    = 0.0    ! rigid lid at surface
            w_strct(kc+1) = 0.0 ! zero-flux at bottom
            
            ! Normalize vertical structure function of w such that
            ! \int(w_strct)^2dz = a_int (a_int could be any value, e.g., 0.5)
            nz = kc+1 ! number of layer interfaces (including surface and bottom)
            w2avg = 0.0
            do k=1,nz-1
              dz(K) = Hc(K)
              w2avg = w2avg + 0.5*(w_strct(k)+w_strct(k+1))*dz(K)
            enddo
            w2avg = w2avg/htot(i,j)
            w_strct = w_strct/sqrt(htot(i,j)*w2avg*I_a_int)
            
            ! Calculate vertical structure function of u (i.e. dw/dz)
            do k=2,nz-1
              u_strct(k) = (w_strct(k-1) - w_strct(k)  )/dz(K-1) + &
                           (w_strct(k)   - w_strct(k+1))/dz(K)
            enddo
            u_strct(1)   = (w_strct(1)   -  w_strct(2) )/dz(1)
            u_strct(nz)  = (w_strct(nz-1)-  w_strct(nz))/dz(nz-1)
            
            ! Calculate wavenumber magnitude
            f2 = 0.25*((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                (G%CoriolisBu(I,J-1)**2 + G%CoriolisBu(I-1,J)**2))
            Kmag2 = (freq**2 - f2) / (cn(i,j)**2 + cg_subRO**2)
            
            ! Calculate terms in vertically integrated energy equation
            int_dwdz2 = 0.0 ; int_w2 = 0.0 ; int_N2w2 = 0.0
            u_strct2 = u_strct(1:nz)**2
            w_strct2 = w_strct(1:nz)**2
            ! vertical integration with Trapezoidal rule
            do k=1,nz-1
              int_dwdz2 = int_dwdz2 + 0.5*(u_strct2(k)+u_strct2(k+1))*dz(K)
              int_w2    = int_w2    + 0.5*(w_strct2(k)+w_strct2(k+1))*dz(K)
              int_N2w2  = int_N2w2  + 0.5*(w_strct2(k)*N2(k)+w_strct2(k+1)*N2(k+1))*dz(K)
            enddo
            
            KE_term = 0.25*G%Rho0*( (1+f2/freq**2)/Kmag2*int_dwdz2 + int_w2 )
            PE_term = 0.25*G%Rho0*( int_N2w2/freq**2 )
            
            ! Back-calculate amplitude from energy equation
            W0 = sqrt( En(i,j)/(KE_term + PE_term) )
            
            ! Calculate actual vertical velocity profile and derivative
            W_profile    = W0*w_strct
            dWdz_profile = W0*u_strct
            
            ! Calculate average magnitude of actual horizontal velocity over a period
            Uavg_profile = abs(dWdz_profile) * sqrt((1+f2/freq**2)/(2.0*Kmag2))
             
            ! Store values in control structure
            CS%w_strct(i,j,1:kc+1)     = w_strct
            CS%u_strct(i,j,1:kc+1)     = u_strct
            CS%W_profile(i,j,1:kc+1)   = W_profile
            CS%Uavg_profile(i,j,1:kc+1)= Uavg_profile
            CS%z_depths(i,j,1:kc+1)    = z_int
            CS%N2(i,j,1:kc+1)          = N2
            CS%num_intfaces(i,j)       = nz
          endif  ! kc >= 2?
        endif ! drxh_sum >= 0?
      endif ! mask2dT > 0.5?
    enddo ! i-loop
  enddo ! j-loop
  
end subroutine wave_structure

subroutine tridiag_solver(a,b,c,h,y,method,x)
  real, dimension(:), intent(in)  :: a, b, c, h, y
  character(len=*), intent(in)    :: method
  real, dimension(:), intent(out) :: x

!    This subroutine solves a tri-diagonal system Ax=y using either the standard
! Thomas algorithim (TDMA_T) or its more stable variant that invokes the
! "Hallberg substitution" (TDMA_H).  
!
! Arguments: 
!  (in)      a - lower diagonal with first entry equal to zero
!  (in)      b - middle diagonal
!  (in)      c - upper diagonal with last entry equal to zero
!  (in)      h - vector of values that have already been added to b; used for 
!                systems of the form (e.g. average layer thickness in vertical diffusion case):
!                [ -alpha(k-1/2) ]                       * e(k-1) +
!                [  alpha(k-1/2) + alpha(k+1/2) + h(k) ] * e(k)   +
!                [ -alpha(k+1/2) ]                       * e(k+1) = y(k)
!                where a(k)=[-alpha(k-1/2)], b(k)=[alpha(k-1/2)+alpha(k+1/2) + h(k)],
!                and c(k)=[-alpha(k+1/2)].
!  (in)      y - vector of known values on right hand side
!  (out)     x - vector of unknown values to solve for

  integer :: nrow                         ! number of rows in A matrix
  real, allocatable, dimension(:) :: c_prime, y_prime, q, alpha
                                          ! intermediate values for solvers
  real    :: Q_prime, beta                ! intermediate values for solver
  integer :: k                            ! row (e.g. interface) index
  
  nrow = size(y)
  
  if (method == 'TDMA_T') then
    ! Standard Thomas algoritim (4th variant)
    c_prime(:) = 0.0       ; y_prime(:) = 0.0
    c_prime(1) = c(1)/b(1) ; y_prime(1) = y(1)/b(1)
    ! Forward sweep
    do k=2,nrow-1
      c_prime(k) = c(k)/(b(k)-a(k)*c_prime(k-1))
    enddo
    do k=2,nrow
      y_prime(k) = (y(k)-a(k)*y_prime(k-1))/(b(k)-a(k)*c_prime(k-1))
    enddo
    x(nrow) = y_prime(nrow)
    ! Backward sweep
    do k=nrow-1,-1,1
      x(k) = y_prime(k)-c_prime(k)*x(k+1)
    enddo
  elseif (method == 'TDMA_H') then
    ! Thomas algoritim (4th variant) w/ Hallberg substitution.
    ! For a layered system where k is at interfaces, alpha{k+1/2} refers to
    ! some property (e.g. inverse thickness for mode-structure problem) of the 
    ! layer below and alpha{k-1/2} refers to the layer above. 
    ! Here, alpha(k)=alpha{k+1/2} and alpha(k-1)=alpha{k-1/2}.
    ! Formulation requires symmetry - check for it
    do k=1,nrow-1
      if (abs(a(k+1)-c(k)) > 1.e-30) then
        call MOM_error(WARNING, "tridiag_solver: matrix not symmetric; need symmetry when invoking TDMA_H")
      endif
    enddo
    alpha = -c
    ! Alpha of the bottom-most layer is not necessarily zero. Therefore,
    ! back out the value from the provided b(nrow and h(nrow) values
    alpha(nrow)   = b(nrow)-h(nrow)-alpha(nrow-1)
    ! Prime other variables
    beta       = 1/b(1)
    y_prime(:) = 0.0       ; q(:) = 0.0
    y_prime(1) = beta*y(1) ; q(1) = beta*alpha(1)
    Q_prime    = 1-q(1)
    ! Forward sweep
    do k=2,nrow-1
      beta = 1/(h(k)+alpha(k-1)*Q_prime+alpha(k))
      q(k) = beta*alpha(k)
      y_prime(k) = beta*(y(k)+alpha(k-1)*y_prime(k-1))
      Q_prime = beta*(h(k)+alpha(k-1)*Q_prime)
    enddo
    beta = 1/(h(nrow)+alpha(nrow-1)*Q_prime+alpha(nrow))
    y_prime(nrow) = beta*(y(nrow)+alpha(nrow-1)*y_prime(nrow-1))
    x(nrow) = y_prime(nrow)
    ! Backward sweep
    do k=nrow-1,-1,1
      x(k) = y_prime(k)+q(k)*x(k+1)
    enddo
  endif
  
end subroutine tridiag_solver

subroutine wave_structure_init(Time, G, param_file, diag, CS)
  type(time_type),             intent(in)    :: Time
  type(ocean_grid_type),       intent(in)    :: G
  type(param_file_type),       intent(in)    :: param_file
  type(diag_ctrl), target,     intent(inout) :: diag
  type(wave_structure_CS),     pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                  for this module
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_wave_structure"  ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "wave_structure_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag

  ! Write all relevant parameters to the model log.
  call log_version(param_file, mod, version, "")

end subroutine wave_structure_init

end module MOM_wave_structure
