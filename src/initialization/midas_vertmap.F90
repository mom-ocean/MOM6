!> Routines for initialization callable from MOM6 or Python (MIDAS)
module MIDAS_vertmap

! This file is part of MOM6. See LICENSE.md for the license.

! If calling from MOM6, use MOM6 interfaces for EOS functions
#ifndef PY_SOLO
use MOM_EOS, only : EOS_type, calculate_density,calculate_density_derivs

implicit none ; private

public tracer_z_init, determine_temperature, fill_boundaries
public find_interfaces, meshgrid
#endif

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Fill grid edges
interface fill_boundaries
   module procedure fill_boundaries_real
   module procedure fill_boundaries_int
end interface

! real, parameter :: epsln=1.e-10 !< A hard-wired constant!
                                  !! \todo Get rid of this constant

contains

#ifdef PY_SOLO
!> Calculate seawater equation of state, given T[degC], S[PSU], and p[Pa]
!! Returns density [kg m-3]
!!
!! These EOS routines are needed only for the stand-alone version of the code
!! The subroutines in this file implement the equation of state for
!! sea water using the formulae given by  Wright, 1997, J. Atmos.
!! Ocean. Tech., 14, 735-740.
function wright_eos_2d(T,S,p) result(rho)
  real(kind=8), dimension(:,:), intent(in) :: T,S !< temperature [degC] and Salinity [psu]
  real, intent(in) :: p  !< pressure [Pa]
  real(kind=8), dimension(size(T,1),size(T,2)) :: rho !< potential density [kg m-3]
  ! Local variables
  real(kind=8) :: a0,a1,a2,b0,b1,b2,b3,b4,b5,c0,c1,c2,c3,c4,c5
  real(kind=8) :: al0,lam,p0,I_denom
  integer :: i,k

  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4
  b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422
  c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464

  do k=1,size(T,2)
    do i=1,size(T,1)
      al0 = a0 + a1*T(i,k) +a2*S(i,k)
      p0  = b0 + b4*S(i,k) + T(i,k) * (b1 + T(i,k)*(b2 + &
           b3*T(i,k)) + b5*S(i,k))
      lam = c0 +c4*S(i,k) + T(i,k) * (c1 + T(i,k)*(c2 + &
           c3*T(i,k)) + c5*S(i,k))
      I_denom = 1.0 / (lam + al0*(p+p0))
      rho(i,k) = (p + p0) * I_denom
    enddo
  enddo

  return
end function wright_eos_2d

!> Calculate seawater thermal expansion coefficient given T[degC],S[PSU],p[Pa]
!! Returns density [kg m-3 degC-1]
!!
!! The subroutines in this file implement the equation of state for
!! sea water using the formulae given by  Wright, 1997, J. Atmos.
!! Ocean. Tech., 14, 735-740.
function alpha_wright_eos_2d(T,S,p) result(drho_dT)
  real(kind=8), dimension(:,:), intent(in) :: T,S !< temperature [degC] and Salinity [psu]
  real, intent(in) :: p !< pressure [Pa]
  real(kind=8), dimension(size(T,1),size(T,2)) :: drho_dT !< partial derivative of density with
                                                          !! respect to temperature [kg m-3 degC-1]
  ! Local variables
  real(kind=8) :: a0,a1,a2,b0,b1,b2,b3,b4,b5,c0,c1,c2,c3,c4,c5
  real(kind=8) :: al0,lam,p0,I_denom,I_denom2
  integer :: i,k

  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4
  b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422
  c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464

  do k=1,size(T,2)
    do i=1,size(T,1)
      al0 = a0 + a1*T(i,k) +a2*S(i,k)
      p0  = b0 + b4*S(i,k) + T(i,k) * (b1 + T(i,k)*(b2 + &
           b3*T(i,k)) + b5*S(i,k))
      lam = c0 +c4*S(i,k) + T(i,k) * (c1 + T(i,k)*(c2 +  &
           c3*T(i,k)) + c5*S(i,k))
      I_denom = 1.0 / (lam + al0*(p+p0))
      I_denom2 = I_denom*I_denom
      drho_dT(i,k) = I_denom2*(lam*(b1+T(i,k)*(2*b2 + &
           3*b3*T(i,k)) + b5*S(i,k)) - (p+p0)*((p+p0)*a1 + &
           (c1+T(i,k)*(2*c2 + 3*c3*T(i,k)) + c5*S(i,k))))
    enddo
  enddo

  return
end function alpha_wright_eos_2d

!> Calculate seawater haline expansion coefficient given T[degC],S[PSU],p[Pa]
!! Returns density [kg m-3 PSU-1]
!!
!! The subroutines in this file implement the equation of state for
!! sea water using the formulae given by  Wright, 1997, J. Atmos.
!! Ocean. Tech., 14, 735-740.
function beta_wright_eos_2d(T,S,p) result(drho_dS)
  real(kind=8), dimension(:,:), intent(in) :: T,S !< temperature [degC] and salinity [psu]
  real, intent(in) :: p !< pressure [Pa]
  real(kind=8), dimension(size(T,1),size(T,2)) :: drho_dS !< partial derivative of density with
                                                          !! respect to salinity [kg m-3 PSU-1]
  ! Local variables
  real(kind=8) :: a0,a1,a2,b0,b1,b2,b3,b4,b5,c0,c1,c2,c3,c4,c5
  real(kind=8) :: al0,lam,p0,I_denom,I_denom2
  integer :: i,k

  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4
  b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422
  c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464

  do k=1,size(T,2)
    do i=1,size(T,1)
      al0 = a0 + a1*T(i,k) +a2*S(i,k)
      p0  = b0 + b4*S(i,k) + T(i,k) * (b1 + T(i,k)*(b2 + &
           b3*T(i,k)) + b5*S(i,k))
      lam = c0 +c4*S(i,k) + T(i,k) * (c1 + T(i,k)*(c2 + &
           c3*T(i,k)) + c5*S(i,k))
      I_denom = 1.0 / (lam + al0*(p+p0))
      I_denom2 = I_denom*I_denom
      drho_dS(i,k) = I_denom2*(lam*(b4+b5*T(i,k)) - &
           (p+p0)*((p+p0)*a2 + (c4+c5*T(i,k))))
    enddo
  enddo

  return
end function beta_wright_eos_2d
#endif

!> Layer model routine for remapping tracers
function tracer_z_init(tr_in, z_edges, e, nkml, nkbl, land_fill, wet, nlay, nlevs, &
                       debug, i_debug, j_debug, eps_z) result(tr)
  real, dimension(:,:,:),           intent(in) :: tr_in !< The z-space array of tracer concentrations that is read in.
  real, dimension(size(tr_in,3)+1), intent(in) :: z_edges !< The depths of the cell edges in the input z* data
                                                       !! [Z ~> m or m]
  integer,                          intent(in) :: nlay !< The number of vertical layers in the target grid
  real, dimension(size(tr_in,1),size(tr_in,2),nlay+1), &
                                    intent(in) :: e !< The depths of the target layer interfaces [Z ~> m or m]
  integer,                          intent(in) :: nkml !< The number of mixed layers
  integer,                          intent(in) :: nkbl !< The number of buffer layers
  real,                             intent(in) :: land_fill !< fill in data over land (1)
  real, dimension(size(tr_in,1),size(tr_in,2)), &
                                    intent(in) :: wet !< The wet mask for the source data (valid points)
  real, dimension(size(tr_in,1),size(tr_in,2)), &
                          optional, intent(in) :: nlevs !< The number of input levels with valid data
  logical,                optional, intent(in) :: debug !< optional debug flag
  integer,                optional, intent(in) :: i_debug !< i-index of point for debugging
  integer,                optional, intent(in) :: j_debug !< j-index of point for debugging
  real,                   optional, intent(in) :: eps_z   !< A negligibly small layer thickness [Z ~> m or m].
  real, dimension(size(tr_in,1),size(tr_in,2),nlay) :: tr !< tracers in layer space

  ! Local variables
  real, dimension(size(tr_in,3)) :: tr_1d !< a copy of the input tracer concentrations in a column.
  real, dimension(nlay+1) :: e_1d ! A 1-d column of intreface heights, in the same units as e.
  real, dimension(nlay) :: tr_    ! A 1-d column of tracer concentrations
  integer, dimension(size(tr_in,1),size(tr_in,2)) :: nlevs_data !< number of valid levels in the input dataset
  integer :: n,i,j,k,l,nx,ny,nz,nt,kz
  integer :: k_top,k_bot,k_bot_prev,kk,kstart
  real    :: sl_tr    ! The tracer concentration slope times the layer thickness, in tracer units.
  real    :: epsln_Z  ! A negligibly thin layer thickness [Z ~> m].
  real, dimension(size(tr_in,3)) :: wt !< The fractional weight for each layer in the range between z1 and z2
  real, dimension(size(tr_in,3)) :: z1, z2 ! z1 and z2 are the fractional depths of the top and bottom
                                  ! limits of the part of a z-cell that contributes to a layer, relative
                                  ! to the cell center and normalized by the cell thickness [nondim].
                                  ! Note that -1/2 <= z1 <= z2 <= 1/2.

  logical :: debug_msg, debug_, debug_pt

  nx = size(tr_in,1); ny=size(tr_in,2); nz = size(tr_in,3)

  nlevs_data = size(tr_in,3)
  if (PRESENT(nlevs)) nlevs_data = anint(nlevs)
  epsln_Z = 1.0e-10 ; if (PRESENT(eps_z)) epsln_Z = eps_z

  debug_=.false. ; if (PRESENT(debug)) debug_ = debug
  debug_msg = debug_
  debug_pt = debug_ ; if (PRESENT(i_debug) .and. PRESENT(j_debug)) debug_pt = debug_

  do j=1,ny
    i_loop: do i=1,nx
      if (nlevs_data(i,j) == 0 .or. wet(i,j) == 0.) then
        tr(i,j,:) = land_fill
        cycle i_loop
      endif

      do k=1,nz
        tr_1d(k) = tr_in(i,j,k)
      enddo

      do k=1,nlay+1
        e_1d(k) = e(i,j,k)
      enddo
      k_bot = 1 ; k_bot_prev = -1
      do k=1,nlay
        if (e_1d(k+1) > z_edges(1)) then
          tr(i,j,k) = tr_1d(1)
        elseif (e_1d(k) < z_edges(nlevs_data(i,j)+1)) then
          if (debug_msg) then
            print *,'*** WARNING : Found interface below valid range of z data '
            print *,'(i,j,z_bottom,interface)= ',&
                 i,j,z_edges(nlevs_data(i,j)+1),e_1d(k)
            print *,'z_edges= ',z_edges
            print *,'e=',e_1d
            print *,'*** I will extrapolate below using the bottom-most valid values'
            debug_msg = .false.
          endif
          tr(i,j,k) = tr_1d(nlevs_data(i,j))

        else
          kstart=k_bot
          call find_overlap(z_edges, e_1d(k), e_1d(k+1), nlevs_data(i,j), &
                            kstart, k_top, k_bot, wt, z1, z2)

          if (debug_pt) then ; if ((i == i_debug) .and. (j == j_debug)) then
            print *,'0001 k,k_top,k_bot,sum(wt),sum(z2-z1) = ',k,k_top,k_bot,sum(wt),sum(z2-z1)
          endif ; endif
          kz = k_top
          sl_tr=0.0; ! cur_tr=0.0
          if (kz /= k_bot_prev) then
            ! Calculate the intra-cell profile.
            if ((kz < nlevs_data(i,j)) .and. (kz > 1)) then
               sl_tr = find_limited_slope(tr_1d, z_edges, kz)
            endif
          endif
          if (kz > nlevs_data(i,j)) kz = nlevs_data(i,j)
          ! This is the piecewise linear form.
          tr(i,j,k) = wt(kz) * (tr_1d(kz) + 0.5*sl_tr*(z2(kz) + z1(kz)))
          ! For the piecewise parabolic form add the following...
          !     + C1_3*wt(kz) * cur_tr*(z2(kz)**2 + z2(kz)*z1(kz) + z1(kz)**2))
          if (debug_pt) then ; if ((i == i_debug) .and. (j == j_debug)) then
            print *,'0002 k,k_top,k_bot,k_bot_prev,sl_tr = ',k,k_top,k_bot,k_bot_prev,sl_tr
          endif ; endif

          do kz=k_top+1,k_bot-1
            tr(i,j,k) = tr(i,j,k) + wt(kz)*tr_1d(kz)
          enddo

          if (debug_pt) then ; if ((i == i_debug) .and. (j == j_debug)) then
            print *,'0003 k,tr = ',k,tr(i,j,k)
          endif ; endif

          if (k_bot > k_top) then
            kz = k_bot
            ! Calculate the intra-cell profile.
            sl_tr = 0.0 ! ; cur_tr = 0.0
            if ((kz < nlevs_data(i,j)) .and. (kz > 1)) then
               sl_tr = find_limited_slope(tr_1d, z_edges, kz)
            endif
            ! This is the piecewise linear form.
            tr(i,j,k) = tr(i,j,k) + wt(kz) * &
                 (tr_1d(kz) + 0.5*sl_tr*(z2(kz) + z1(kz)))
            ! For the piecewise parabolic form add the following...
            !     + C1_3*cur_tr*(z2(kz)**2 + z2(kz)*z1(kz) + z1(kz)**2))

            if (debug_pt) then ; if ((i == i_debug) .and. (j == j_debug)) then
              print *,'0004 k,kz,nlevs,sl_tr,tr = ',k,kz,nlevs_data(i,j),sl_tr,tr(i,j,k)
              print *,'0005 k,kz,tr(kz),tr(kz-1),tr(kz+1) = ',k,kz,tr_1d(kz),tr_1d(kz-1),tr_1d(kz+1),z_edges(kz+2)
            endif ; endif

          endif
          k_bot_prev = k_bot

        endif
      enddo ! k-loop

      do k=2,nlay  ! simply fill vanished layers with adjacent value
        if (e_1d(k)-e_1d(k+1) <= epsln_Z) tr(i,j,k)=tr(i,j,k-1)
      enddo

    enddo i_loop
  enddo

end function tracer_z_init

!> Return the index where to insert item x in list a, assuming a is sorted.
!! The return values [i] is such that all e in a[:i-1] have e <= x, and all e in
!! a[i:] have e > x. So if x already appears in the list, will
!! insert just after the rightmost x already there.
!! Optional args lo (default 1) and hi (default len(a)) bound the
!! slice of a to be searched.
function bisect_fast(a, x, lo, hi) result(bi_r)
  real, dimension(:,:), intent(in) :: a !< Sorted list
  real, dimension(:), intent(in) :: x !< Item to be inserted
  integer, dimension(size(a,1)), optional, intent(in) :: lo !< Lower bracket of optional range to search
  integer, dimension(size(a,1)), optional, intent(in) :: hi !< Upper bracket of optional range to search
  integer, dimension(size(a,1),size(x,1))  :: bi_r

  integer :: mid,num_x,num_a,i
  integer, dimension(size(a,1))  :: lo_,hi_,lo0,hi0
  integer :: nprofs,j

  lo_=1;hi_=size(a,2);num_x=size(x,1);bi_r=-1;nprofs=size(a,1)

  if (PRESENT(lo)) then
     where (lo>0) lo_=lo
  endif
  if (PRESENT(hi)) then
     where (hi>0) hi_=hi
  endif

  lo0=lo_;hi0=hi_

  do j=1,nprofs
    do i=1,num_x
      lo_=lo0;hi_=hi0
      do while (lo_(j) < hi_(j))
        mid = (lo_(j)+hi_(j))/2
        if (x(i) < a(j,mid)) then
           hi_(j) = mid
        else
           lo_(j) = mid+1
        endif
      enddo
      bi_r(j,i)=lo_(j)
    enddo
  enddo


  return

end function bisect_fast

#ifdef PY_SOLO
!  Only for stand-alone python

!> This subroutine determines the potential temperature and salinity that
!! is consistent with the target density using provided initial guess
subroutine determine_temperature(temp, salt, R, p_ref, niter, land_fill, h, k_start)
  real(kind=8), dimension(:,:,:), intent(inout) :: temp !< potential temperature [degC]
  real(kind=8), dimension(:,:,:), intent(inout) :: salt !< salinity [PSU]
  real(kind=8), dimension(size(temp,3)), intent(in) :: R !< desired potential density [kg m-3].
  real, intent(in) :: p_ref !< reference pressure [Pa].
  integer, intent(in) :: niter !< maximum number of iterations
  integer, intent(in) :: k_start !< starting index (i.e. below the buffer layer)
  real, intent(in) :: land_fill !< land fill value
  real(kind=8), dimension(:,:,:), intent(in) :: h !< layer thickness . Do not iterate for massless layers

  ! Local variables
  real, parameter :: T_max = 35.0, T_min = -2.0
#else
!> This subroutine determines the potential temperature and salinity that
!! is consistent with the target density using provided initial guess
subroutine determine_temperature(temp, salt, R, p_ref, niter, land_fill, h, k_start, eos)
  real, dimension(:,:,:),        intent(inout) :: temp !< potential temperature [degC]
  real, dimension(:,:,:),        intent(inout) :: salt !< salinity [PSU]
  real, dimension(size(temp,3)), intent(in)    :: R !< desired potential density [kg m-3].
  real,                          intent(in)    :: p_ref !< reference pressure [Pa].
  integer,                       intent(in)    :: niter !< maximum number of iterations
  integer,                       intent(in)    :: k_start !< starting index (i.e. below the buffer layer)
  real,                          intent(in)    :: land_fill !< land fill value
  real, dimension(:,:,:),        intent(in)    :: h   !< layer thickness, used only to avoid working on massless layers
  type(eos_type),                pointer       :: eos !< seawater equation of state control structure

  real, parameter :: T_max = 31.0, T_min = -2.0
#endif
  ! Local variables (All of which need documentation!)
  real(kind=8), dimension(size(temp,1),size(temp,3)) :: T, S, dT, dS, rho, hin
  real(kind=8), dimension(size(temp,1),size(temp,3)) :: drho_dT, drho_dS
  real(kind=8), dimension(size(temp,1)) :: press
  integer :: nx, ny, nz, nt, i, j, k, n, itt
  real    :: dT_dS
  logical :: adjust_salt, old_fit
  real, parameter :: S_min = 0.5, S_max=65.0
  real, parameter :: tol=1.e-4, max_t_adj=1.0, max_s_adj = 0.5

  old_fit = .true.   ! reproduces siena behavior
  ! will switch to the newer method which simultaneously adjusts
  ! temp and salt based on the ratio of the thermal and haline coefficients.

  nx=size(temp,1) ; ny=size(temp,2) ; nz=size(temp,3)

  press(:) = p_ref

  do j=1,ny
    dS(:,:) = 0. ! Needs to be zero everywhere since there is a maxval(abs(dS)) later...
    T=temp(:,j,:)
    S=salt(:,j,:)
    hin=h(:,j,:)
    dT=0.0
    adjust_salt = .true.
    iter_loop: do itt = 1,niter
#ifdef PY_SOLO
      rho=wright_eos_2d(T,S,p_ref)
      drho_dT=alpha_wright_eos_2d(T,S,p_ref)
#else
      do k=1, nz
        call calculate_density(T(:,k), S(:,k), press, rho(:,k), 1, nx, eos)
        call calculate_density_derivs(T(:,k), S(:,k), press, drho_dT(:,k), drho_dS(:,k), 1, nx, eos)
      enddo
#endif
      do k=k_start,nz ; do i=1,nx

!       if (abs(rho(i,k)-R(k))>tol .and. hin(i,k)>epsln .and. abs(T(i,k)-land_fill) < epsln) then
        if (abs(rho(i,k)-R(k))>tol) then
          if (old_fit) then
             dT(i,k) = max(min((R(k)-rho(i,k)) / drho_dT(i,k), max_t_adj), -max_t_adj)
             T(i,k) = max(min(T(i,k)+dT(i,k), T_max), T_min)
          else
             dT_dS = 10.0 - min(-drho_dT(i,k)/drho_dS(i,k),10.)
             !### RWH: Based on the dimensions alone, the expression above should be:
             ! dT_dS = 10.0 - min(-drho_dS(i,k)/drho_dT(i,k),10.)
             dS(i,k) = (R(k)-rho(i,k)) / (drho_dS(i,k) - drho_dT(i,k)*dT_dS )
             dT(i,k) = -dT_dS*dS(i,k)
           ! dT(i,k) = max(min(dT(i,k), max_t_adj), -max_t_adj)
             T(i,k) = max(min(T(i,k)+dT(i,k), T_max), T_min)
             S(i,k) = max(min(S(i,k)+dS(i,k), S_max), S_min)
          endif
        endif
      enddo ; enddo
      if (maxval(abs(dT)) < tol) then
         adjust_salt = .false.
         exit iter_loop
      endif
    enddo iter_loop

    if (adjust_salt .and. old_fit) then ; do itt = 1,niter
#ifdef PY_SOLO
      rho = wright_eos_2d(T,S,p_ref)
      drho_dS = beta_wright_eos_2d(T,S,p_ref)
#else
      do k=1, nz
        call calculate_density(T(:,k),S(:,k),press,rho(:,k),1,nx,eos)
        call calculate_density_derivs(T(:,k),S(:,k),press,drho_dT(:,k),drho_dS(:,k),1,nx,eos)
      enddo
#endif
      do k=k_start,nz ; do i=1,nx
!       if (abs(rho(i,k)-R(k))>tol .and. hin(i,k)>epsln .and. abs(T(i,k)-land_fill) < epsln ) then
        if (abs(rho(i,k)-R(k)) > tol) then
          dS(i,k) = max(min((R(k)-rho(i,k)) / drho_dS(i,k), max_s_adj), -max_s_adj)
          S(i,k) = max(min(S(i,k)+dS(i,k), S_max), S_min)
        endif
      enddo ; enddo
      if (maxval(abs(dS)) < tol) exit
    enddo ; endif

    temp(:,j,:)=T(:,:)
    salt(:,j,:)=S(:,:)
  enddo

end subroutine determine_temperature

!> This subroutine determines the layers bounded by interfaces e that overlap
!! with the depth range between Z_top and Z_bot, and also the fractional weights
!! of each layer. It also calculates the normalized relative depths of the range
!! of each layer that overlaps that depth range.
!! Note that by convention, e decreases with increasing k and Z_top > Z_bot.
subroutine find_overlap(e, Z_top, Z_bot, k_max, k_start, k_top, k_bot, wt, z1, z2)
  real, dimension(:), intent(in)  :: e  !< The interface positions, [Z ~> m] or other units.
  real,               intent(in)  :: Z_top !< The top of the range being mapped to, [Z ~> m] or other units.
  real,               intent(in)  :: Z_bot !< The bottom of the range being mapped to, [Z ~> m] or other units.
  integer,            intent(in)  :: k_max !< The number of valid layers.
  integer,            intent(in)  :: k_start !< The layer at which to start searching.
  integer,            intent(out) :: k_top !< The index of the top layer that overlap with the depth range.
  integer,            intent(out) :: k_bot !< The index of the bottom layer that overlap with the depth range.
  real, dimension(:), intent(out) :: wt !< The relative weights of each layer from k_top to k_bot [nondim].
  real, dimension(:), intent(out) :: z1 !< Depth of the top limit of layer that contributes to a level [nondim].
  real, dimension(:), intent(out) :: z2 !< Depth of the bottom limit of layer that contributes to a level [nondim].

  ! Local variables
  real :: Ih, e_c, tot_wt, I_totwt
  integer :: k

  wt(:)=0.0 ; z1(:)=0.0 ; z2(:)=0.0
  k_top = k_start ; k_bot = k_start ; wt(1) = 1.0 ; z1(1) = -0.5 ; z2(1) = 0.5

  do k=k_start,k_max ; if (e(K+1) < Z_top) exit ; enddo
  k_top = k

  if (k>k_max) return

  ! Determine the fractional weights of each layer.
  ! Note that by convention, e and Z_int decrease with increasing k.
  if (e(K+1) <= Z_bot) then
    wt(k) = 1.0 ; k_bot = k
    Ih = 0.0 ; if (e(K) /= e(K+1)) Ih = 1.0 / (e(K)-e(K+1))
    e_c = 0.5*(e(K)+e(K+1))
    z1(k) = (e_c - MIN(e(K), Z_top)) * Ih
    z2(k) = (e_c - Z_bot) * Ih
  else
    wt(k) = MIN(e(K),Z_top) - e(K+1) ; tot_wt = wt(k) ! These are always > 0.
    ! Ih = 0.0 ; if (e(K) /= e(K+1)) Ih = 1.0 / (e(K)-e(K+1))
    if (e(K) /= e(K+1)) then
      z1(k) = (0.5*(e(K)+e(K+1)) - MIN(e(K), Z_top)) / (e(K)-e(K+1))
    else ; z1(k) = -0.5 ; endif
    z2(k) = 0.5
    k_bot = k_max
    do k=k_top+1,k_max
      if (e(K+1) <= Z_bot) then
        k_bot = k
        wt(k) = e(K) - Z_bot ; z1(k) = -0.5
        if (e(K) /= e(K+1)) then
          z2(k) = (0.5*(e(K)+e(K+1)) - Z_bot) / (e(K)-e(K+1))
        else ; z2(k) = 0.5 ; endif
      else
        wt(k) = e(K) - e(K+1) ; z1(k) = -0.5 ; z2(k) = 0.5
      endif
      tot_wt = tot_wt + wt(k) ! wt(k) is always > 0.
      if (k>=k_bot) exit
    enddo

    I_totwt = 0.0 ; if (tot_wt > 0.0) I_totwt = 1.0 / tot_wt
    do k=k_top,k_bot ; wt(k) = I_totwt*wt(k) ; enddo
  endif

end subroutine find_overlap

!> This subroutine determines a limited slope for val to be advected with
!! a piecewise limited scheme.
function find_limited_slope(val, e, k) result(slope)
  real, dimension(:), intent(in) :: val !< An column the values that are being interpolated.
  real, dimension(:), intent(in) :: e   !< A column's interface heights [Z ~> m] or other units.
  integer,            intent(in) :: k   !< The layer whose slope is being determined.
  real :: slope !< The normalized slope in the intracell distribution of val.
  ! Local variables
  real :: amn, cmn
  real :: d1, d2

  if ((val(k)-val(k-1)) * (val(k)-val(k+1)) >= 0.0) then
    slope = 0.0 ! ; curvature = 0.0
  else
    d1 = 0.5*(e(K-1)-e(K+1)) ; d2 = 0.5*(e(K)-e(K+2))
    if (d1*d2 > 0.0) then
      slope = ((d1**2)*(val(k+1) - val(k)) + (d2**2)*(val(k) - val(k-1))) * &
              (e(K) - e(K+1)) / (d1*d2*(d1+d2))
      ! slope = 0.5*(val(k+1) - val(k-1))
      ! This is S.J. Lin's form of the PLM limiter.
      amn = min(abs(slope), 2.0*(max(val(k-1), val(k), val(k+1)) - val(k)))
      cmn = 2.0*(val(k) - min(val(k-1), val(k), val(k+1)))
      slope = sign(1.0, slope) * min(amn, cmn)

      ! min(abs(slope), 2.0*(max(val(k-1),val(k),val(k+1)) - val(k)), &
      !                 2.0*(val(k) - min(val(k-1),val(k),val(k+1))))
      ! curvature = 0.0
    else
      slope = 0.0 ! ; curvature = 0.0
    endif
  endif

end function find_limited_slope

!> Find interface positions corresponding to density profile
function find_interfaces(rho, zin, Rb, depth, nlevs, nkml, nkbl, hml, debug, eps_z, eps_rho) result(zi)
  real, dimension(:,:,:), &
                      intent(in) :: rho   !< potential density in z-space [kg m-3 or R ~> kg m-3]
  real, dimension(size(rho,3)), &
                      intent(in) :: zin   !< Input data levels [m or Z ~> m].
  real, dimension(:), intent(in) :: Rb    !< target interface densities [kg m-3 or R ~> kg m-3]
  real, dimension(size(rho,1),size(rho,2)), &
                      intent(in) :: depth !< ocean depth [Z ~> m].
  real, dimension(size(rho,1),size(rho,2)), &
            optional, intent(in) :: nlevs !< number of valid points in each column
  logical,  optional, intent(in) :: debug !< optional debug flag
  integer,  optional, intent(in) :: nkml  !< number of mixed layer pieces
  integer,  optional, intent(in) :: nkbl  !< number of buffer layer pieces
  real,     optional, intent(in) :: hml   !< mixed layer depth [Z ~> m].
  real,     optional, intent(in) :: eps_z !< A negligibly small layer thickness [m or Z ~> m].
  real,     optional, intent(in) :: eps_rho !< A negligibly small density difference [kg m-3 or R ~> kg m-3].
  real, dimension(size(rho,1),size(rho,2),size(Rb,1)) :: zi !< The returned interface, in the same units az zin.

  ! Local variables
  real, dimension(size(rho,1),size(rho,3)) :: rho_ ! A slice of densities [R ~> kg m-3]
  real, dimension(size(rho,1)) :: depth_
  logical :: unstable
  integer :: dir
  integer, dimension(size(rho,1),size(Rb,1)) :: ki_
  real, dimension(size(rho,1),size(Rb,1)) :: zi_
  integer, dimension(size(rho,1),size(rho,2)) :: nlevs_data
  integer, dimension(size(rho,1)) :: lo, hi
  real :: slope,rsm,drhodz,hml_
  integer :: n,i,j,k,l,nx,ny,nz,nt
  integer :: nlay,kk,nkml_,nkbl_
  logical :: debug_ = .false.
  real    :: epsln_Z    ! A negligibly thin layer thickness [m or Z ~> m].
  real    :: epsln_rho  ! A negligibly small density change [kg m-3 or R ~> kg m-3].
  real, parameter :: zoff=0.999

  nlay=size(Rb)-1

  zi(:,:,:) = 0.0

  if (PRESENT(debug)) debug_=debug

  nx = size(rho,1); ny=size(rho,2); nz = size(rho,3)
  nlevs_data(:,:) = size(rho,3)

  nkml_ = 0 ;  if (PRESENT(nkml)) nkml_ = max(0, nkml)
  nkbl_ = 0 ;  if (PRESENT(nkbl)) nkbl_ = max(0, nkbl)
  hml_ = 0.0 ; if (PRESENT(hml)) hml_ = hml
  epsln_Z = 1.0e-10 ; if (PRESENT(eps_z)) epsln_Z = eps_z
  epsln_rho = 1.0e-10 ; if (PRESENT(eps_rho)) epsln_rho = eps_rho

  if (PRESENT(nlevs)) then
    nlevs_data(:,:) = nlevs(:,:)
  endif

  do j=1,ny
    rho_(:,:) = rho(:,j,:)
    i_loop: do i=1,nx
      if (debug_) then
        print *,'looking for interfaces, i,j,nlevs= ',i,j,nlevs_data(i,j)
        print *,'initial density profile= ', rho_(i,:)
      endif
      unstable=.true.
      dir=1
      do while (unstable)
        unstable=.false.
        if (dir == 1) then
          do k=2,nlevs_data(i,j)-1
            if (rho_(i,k) - rho_(i,k-1) < 0.0 ) then
              if (k == 2) then
                rho_(i,k-1) = rho_(i,k)-epsln_rho
              else
                drhodz = (rho_(i,k+1)-rho_(i,k-1)) / (zin(k+1)-zin(k-1))
                if (drhodz < 0.0) unstable=.true.
                rho_(i,k) = rho_(i,k-1) + drhodz*zoff*(zin(k)-zin(k-1))
              endif
            endif
          enddo
          dir = -1*dir
        else
          do k=nlevs_data(i,j)-1,2,-1
            if (rho_(i,k+1) - rho_(i,k) < 0.0) then
              if (k == nlevs_data(i,j)-1) then
                rho_(i,k+1) = rho_(i,k-1)+epsln_rho
              else
                drhodz = (rho_(i,k+1)-rho_(i,k-1))/(zin(k+1)-zin(k-1))
                if (drhodz  < 0.0) unstable=.true.
                rho_(i,k) = rho_(i,k+1)-drhodz*(zin(k+1)-zin(k))
              endif
            endif
          enddo
          dir = -1*dir
        endif
      enddo
      if (debug_) then
        print *,'final density profile= ', rho_(i,:)
      endif
    enddo i_loop

    ki_(:,:) = 0
    zi_(:,:) = 0.0
    depth_(:) = -1.0*depth(:,j)
    lo(:) = 1
    hi(:) = nlevs_data(:,j)
    ki_ = bisect_fast(rho_, Rb, lo, hi)
    ki_(:,:) = max(1, ki_(:,:)-1)
    do i=1,nx
      do l=2,nlay
        slope = (zin(ki_(i,l)+1) - zin(ki_(i,l))) / max(rho_(i,ki_(i,l)+1) - rho_(i,ki_(i,l)),epsln_rho)
        zi_(i,l) = -1.0*(zin(ki_(i,l)) + slope*(Rb(l)-rho_(i,ki_(i,l))))
        zi_(i,l) = max(zi_(i,l), depth_(i))
        zi_(i,l) = min(zi_(i,l), -1.0*hml_)
      enddo
      zi_(i,nlay+1) = depth_(i)
      do l=2,nkml_+1
        zi_(i,l) = max(hml_*((1.0-real(l))/real(nkml_)), depth_(i))
      enddo
      do l=nlay,nkml_+2,-1
        if (zi_(i,l) < zi_(i,l+1) + epsln_Z) zi_(i,l) = zi_(i,l+1) + epsln_Z
        if (zi_(i,l) > -1.0*hml_)  zi_(i,l) = max(-1.0*hml_, depth_(i))
      enddo
    enddo
    zi(:,j,:) = zi_(:,:)
  enddo

end function find_interfaces

!> Create a 2d-mesh of grid coordinates from 1-d arrays
subroutine meshgrid(x,y,x_T,y_T)
  real, dimension(:), intent(in) :: x !< input x coordinates
  real, dimension(:), intent(in) :: y !< input y coordinates
  real, dimension(size(x,1),size(y,1)), intent(inout) :: x_T !< output 2-d version
  real, dimension(size(x,1),size(y,1)), intent(inout) :: y_T !< output 2-d version

  integer :: ni,nj,i,j

  ni=size(x,1);nj=size(y,1)

  do j=1,nj
    x_T(:,j)=x(:)
  enddo

  do i=1,ni
    y_T(i,:)=y(:)
  enddo

  return

end subroutine meshgrid

!> Solve del2 (zi) = 0 using successive iterations
!! with a 5 point stencil. Only points fill==1 are
!! modified. Except where bad==1, information propagates
!! isotropically in index space.  The resulting solution
!! in each region is an approximation to del2(zi)=0 subject to
!! boundary conditions along the valid points curve bounding this region.
subroutine smooth_heights(zi,fill,bad,sor,niter,cyclic_x, tripolar_n)
  real, dimension(:,:), intent(inout) :: zi !< interface positions [m] or arbitrary
  integer, dimension(size(zi,1),size(zi,2)), intent(in) :: fill !< points to be smoothed
  integer, dimension(size(zi,1),size(zi,2)), intent(in) :: bad !< ignore these points
  real, intent(in)  :: sor !< successive over-relaxation coefficient (typically 0.6)
  integer, intent(in) :: niter !< maximum number of iterations
  logical, intent(in) :: cyclic_x !< input grid cyclic condition in the zonal direction
  logical, intent(in) :: tripolar_n !< tripolar Arctic fold flag

  integer :: i,j,k,n
  integer :: ni,nj

  real, dimension(size(zi,1),size(zi,2)) :: res, m
  integer, dimension(size(zi,1),size(zi,2),4) :: B
  real, dimension(0:size(zi,1)+1,0:size(zi,2)+1) :: mp
  integer, dimension(0:size(zi,1)+1,0:size(zi,2)+1) :: nm

  real :: Isum, bsum

  ni=size(zi,1); nj=size(zi,2)


  mp=fill_boundaries(zi,cyclic_x,tripolar_n)

  B(:,:,:)=0.0
  nm=fill_boundaries(bad,cyclic_x,tripolar_n)

  do j=1,nj
    do i=1,ni
      if (fill(i,j) == 1) then
         B(i,j,1)=1-nm(i+1,j);B(i,j,2)=1-nm(i-1,j)
         B(i,j,3)=1-nm(i,j+1);B(i,j,4)=1-nm(i,j-1)
      endif
    enddo
  enddo

  do n=1,niter
    do j=1,nj
      do i=1,ni
        if (fill(i,j) == 1) then
           bsum = real(B(i,j,1)+B(i,j,2)+B(i,j,3)+B(i,j,4))
           Isum = 1.0/bsum
           res(i,j)=Isum*(B(i,j,1)*mp(i+1,j)+B(i,j,2)*mp(i-1,j)+&
                B(i,j,3)*mp(i,j+1)+B(i,j,4)*mp(i,j-1)) - mp(i,j)
        endif
      enddo
    enddo
    res(:,:)=res(:,:)*sor

    do j=1,nj
      do i=1,ni
        mp(i,j)=mp(i,j)+res(i,j)
      enddo
    enddo

    zi(:,:)=mp(1:ni,1:nj)
    mp = fill_boundaries(zi,cyclic_x,tripolar_n)
  enddo

  return

end subroutine smooth_heights

!> Fill grid edges
function fill_boundaries_int(m,cyclic_x,tripolar_n) result(mp)
  integer, dimension(:,:), intent(in) :: m !< input array
  logical, intent(in) :: cyclic_x !< zonal cyclic condition
  logical, intent(in) :: tripolar_n  !< northern fold condition
  integer, dimension(0:size(m,1)+1,0:size(m,2)+1) :: mp !< output filled array
  ! Local variables
  real, dimension(size(m,1),size(m,2)) :: m_real
  real, dimension(0:size(m,1)+1,0:size(m,2)+1) :: mp_real

  m_real = real(m)

  mp_real = fill_boundaries_real(m_real,cyclic_x,tripolar_n)

  mp = int(mp_real)

  return

end function fill_boundaries_int

!> fill grid edges
function fill_boundaries_real(m,cyclic_x,tripolar_n) result(mp)
  real, dimension(:,:), intent(in) :: m !< input array
  logical, intent(in) :: cyclic_x !< zonal cyclic condition
  logical, intent(in) :: tripolar_n  !< northern fold condition
  real, dimension(0:size(m,1)+1,0:size(m,2)+1) :: mp !< output filled array

  integer :: ni,nj,i,j

  ni=size(m,1); nj=size(m,2)

  mp(1:ni,1:nj)=m(:,:)

  if (cyclic_x) then
     mp(0,1:nj)=m(ni,1:nj)
     mp(ni+1,1:nj)=m(1,1:nj)
  else
     mp(0,1:nj)=m(1,1:nj)
     mp(ni+1,1:nj)=m(ni,1:nj)
  endif

  mp(1:ni,0)=m(1:ni,1)
  if (tripolar_n) then
     do i=1,ni
       mp(i,nj+1)=m(ni-i+1,nj)
     enddo
  else
     mp(1:ni,nj+1)=m(1:ni,nj)
  endif

  return

end function fill_boundaries_real

end module MIDAS_vertmap
