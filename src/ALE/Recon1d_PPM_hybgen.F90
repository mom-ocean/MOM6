!> Piecewise Parabolic Method 1D reconstruction following Colella and Woodward, 1984
!!
!! This implementation of PPM follows Colella and Woodward, 1984 \cite colella1984, with
!! cells resorting to PCM for extrema including first and last cells in column. The algorithm was
!! first ported from Hycom as hybgen_ppm_coefs() in the mom_hybgen_remap module. This module is
!! a refactor to facilitate more complete testing and evaluation.
!!
!! The mom_hybgen_remap.hybgen_ppm_coefs() function (reached with "PPM_HYGEN"),
!! regrid_edge_values.edge_values_explicit_h4cw() function followed by ppm_functions.ppm_reconstruction()
!! (reached with "PPM_CW"), are equivalent. Similarly recon1d_ppm_hybgen (this implementation) is equivalent also.
module Recon1d_PPM_hybgen

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_type, only : testing
use Recon1d_PPM_CW, only : PPM_CW

implicit none ; private

public PPM_hybgen, testing

!> PPM reconstruction following White and Adcroft, 2008
!!
!! Implemented by extending recon1d_ppm_cwk.
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_ppm_cw.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_ppm_cw.average()
!! - f()                    -> recon1d_ppm_cw.f()
!! - dfdx()                 -> recon1d_ppm_cw.dfdx()
!! - check_reconstruction()    *locally defined
!! - unit_tests()           -> recon1d_ppm_cw.unit_tests()
!! - destroy()              -> recon1d_ppm_cw.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (PPM_CW) :: PPM_hybgen

contains
  !> Implementation of the PPM_hybgen reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of check reconstruction for the PPM_hybgen reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PPM_hybgen reconstruction
  procedure :: unit_tests => unit_tests

end type PPM_hybgen

contains

!> Calculate a 1D PPM_hybgen reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(PPM_hybgen), intent(inout) :: this !< This reconstruction
  real,              intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,              intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  real :: h0, h1, h2, h3 ! Cell thickness h(k-2), h(k-1), h(k), h(k+1) in K loop [H]
  real :: h01_h112, h23_h122 ! Approximately 2/3 [nondim]
  real :: h112, h122 ! Approximately 3 h [H]
  real :: ddh ! Approximately 0 [nondim]
  real :: I_h12, I_h01, I_h0123 ! Reciprocals of d12 and sum(h) [H-1]
  real :: dul, dur ! Left and right cell PLM slopes [A]
  real :: u0, u1, u2 ! Far left, left, and right cell values [A]
  real :: edge ! Edge value between cell k-1 and k [A]
  real :: u_min, u_max ! Minimum and maximum value across edge [A]
  real :: a6 ! Colella and Woodward curvature [A]
  real :: du, duc ! Difference between edges across cell [A]
  real :: slp(this%n) ! PLM slope [A]
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
                                    ! differences across the cell [A]
  real :: slope_x_h             ! retained PLM slope times  half grid step [A]
  real :: edge_l, edge_r        ! Edge values (left and right) [A]
  real :: expr1, expr2          ! Temporary expressions [A2]
  real :: u0_avg                ! avg value at given edge [A]
  integer :: k, n, km1, kp1

  n = this%n

  ! First populate the PLM reconstructions
  slp(1) = 0.
  do k = 2, n-1
    h0 = max( this%h_neglect, h(k-1) )
    h1 = max( this%h_neglect, h(k) )
    h2 = max( this%h_neglect, h(k+1) )
    dul = u(k) - u(k-1)
    dur = u(k+1) - u(k)
    h112 = ( 2.0 * h0 + h1 )
    h122 = ( h1 + 2.0 * h2 )
    I_h01 = 1. / ( h0 + h1 )
    I_h12 = 1. / ( h1 + h2 )
    h01_h112 = ( 2.0 * h0 + h1 ) / ( h0 + h1 ) ! When uniform -> 3/2
    h23_h122 = ( 2.0 * h2 + h1 ) / ( h2 + h1 ) ! When uniform -> 3/2
    if ( dul * dur > 0.) then
      du = ( h1 / ( h1 + ( h0 + h2 ) ) ) * ( h112 * dur * I_h12 + h122 * dul * I_h01 )
      slp(k) = sign( min( abs(2.0 * dul), abs(du), abs(2.0 * dur) ), du)
    else
      slp(k) = 0.
    endif
  enddo
  slp(n) = 0.

  this%ul(1) = u(1) ! PCM
  this%ur(1) = u(1) ! PCM
  this%ul(2) = u(1) ! PCM
  do K = 3, n-1 ! K=3 is interface between cells 2 and 3
    h0 = max( this%h_neglect, h(k-2) )
    h1 = max( this%h_neglect, h(k-1) )
    h2 = max( this%h_neglect, h(k) )
    h3 = max( this%h_neglect, h(k+1) )
    h01_h112 = ( h0 + h1 ) / ( 2. * h1 + h2 )     ! When uniform -> 2/3
    h23_h122 = ( h2 + h3 ) / ( h1 + 2. * h2 )     ! When uniform -> 2/3
    ddh = h01_h112 - h23_h122                     ! When uniform -> 0
    I_h12 = 1.0 / ( h1 + h2 )                     ! When uniform -> 1/(2h)
    I_h0123 = 1.0 / ( ( h0 + h1 ) + ( h2 + h3 ) ) ! When uniform -> 1/(4h)
    dul = slp(k-1)
    dur = slp(k)
    u1 = u(k-1)
    u2 = u(k)
    edge = I_h12 * ( h2 * u1 + h1 * u2 ) &                              ! 1/2 u1 + 1/2 u2
         + I_h0123 * ( 2.0 * h1 * h2 * I_h12 * ( u2 - u1 ) * ddh &      ! 0
                     + ( h2 * dul * h23_h122 - h1 * dur * h01_h112 ) )  ! 1/6 dul - 1/6 dur
    this%ur(k-1) = edge
    this%ul(k) = edge
  enddo
  this%ur(n-1) = u(n) ! PCM
  this%ur(n) = u(n) ! PCM
  this%ul(n) = u(n) ! PCM

  do K = 2, n ! K=2 is interface between cells 1 and 2
    u0 = u(k-1)
    u1 = u(k)
    u2 = u(k+1)
    a6 = 3.0 * ( ( u1 - this%ul(k) ) + ( u1 - this%ur(k) ) )
    a6 = 6.0 * u1 - 3.0 * ( this%ul(k) + this%ur(k) )
    du = this%ur(k) - this%ul(k)
    if ( ( u2 - u1 ) * ( u1 - u0 ) <- 0.0 ) then ! Large scale extrema
      this%ul(k) = u1
      this%ur(k) = u1
    elseif ( du * a6 > du * du ) then ! Extrema on right
      edge = 3.0 * u1 - 2.0 * this%ur(k) ! Subject to round off
    ! u_min = min( u0, u1 )
    ! u_max = max( u0, u1 )
    ! edge = max( min( edge, u_max), u_min )
      this%ul(k) = edge
    elseif ( du * a6 < - du * du ) then ! Extrema on left
      edge = 3.0 * u1 - 2.0 * this%ul(k) ! Subject to round off
    ! u_min = min( u1, u2 )
    ! u_max = max( u1, u2 )
    ! edge = max( min( edge, u_max), u_min )
      this%ur(k) = edge
    endif
  enddo

  ! ### Note that the PPM_HYBGEM option calculated the CW PPM coefficients and then
  ! invoked the OM4-era limiters afterwards, effectively doing the limiters twice.
  ! This second pass does change answers!

  ! Loop on cells to bound edge value
  do k = 1, n

    ! For the sake of bounding boundary edge values, the left neighbor of the left boundary cell
    ! is assumed to be the same as the left boundary cell and the right neighbor of the right
    ! boundary cell is assumed to be the same as the right boundary cell. This effectively makes
    ! boundary cells look like extrema.
    km1 = max(1,k-1) ; kp1 = min(k+1,N)

    slope_x_h = 0.0
    sigma_l = ( u(k) - u(km1) )
    if ( (h(km1) + h(kp1)) + 2.0*h(k) > 0. ) then
      sigma_c = ( u(kp1) - u(km1) ) * ( h(k) / ((h(km1) + h(kp1)) + 2.0*h(k)) )
    else
      sigma_c = 0.
    endif
    sigma_r = ( u(kp1) - u(k) )

    ! The limiter is used in the local coordinate system to each cell, so for convenience store
    ! the slope times a half grid spacing.  (See White and Adcroft JCP 2008 Eqs 19 and 20)
    if ( (sigma_l * sigma_r) > 0.0 ) &
      slope_x_h = sign( min(abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )

    ! Limit the edge values
    if ( (u(km1)-this%ul(k)) * (this%ul(k)-u(k)) < 0.0 ) then
      this%ul(k) = u(k) - sign( min( abs(slope_x_h), abs(this%ul(k)-u(k)) ), slope_x_h )
    endif

    if ( (u(kp1)-this%ur(k)) * (this%ur(k)-u(k)) < 0.0 ) then
      this%ur(k) = u(k) + sign( min( abs(slope_x_h), abs(this%ur(k)-u(k)) ), slope_x_h )
    endif

    ! Finally bound by neighboring cell means in case of roundoff
    this%ul(k) = max( min( this%ul(k), max(u(km1), u(k)) ), min(u(km1), u(k)) )
    this%ur(k) = max( min( this%ur(k), max(u(kp1), u(k)) ), min(u(kp1), u(k)) )

  enddo ! loop on interior edges

  do k = 1, n-1
    if ( (this%ul(k+1) - this%ur(k)) * (u(k+1) - u(k)) < 0.0 ) then
      u0_avg = 0.5 * ( this%ur(k) + this%ul(k+1) )
      u0_avg = max( min( u0_avg, max(u(k), u(k+1)) ), min(u(k), u(k+1)) )
      this%ur(k) = u0_avg
      this%ul(k+1) = u0_avg
    endif
  enddo ! end loop on interior edges

  ! Loop on interior cells to apply the standard
  ! PPM limiter (Colella & Woodward, JCP 84)
  do k = 2, n-1

    ! Get cell averages
    u0 = u(k-1)
    u1 = u(k)
    u2 = u(k+1)

    edge_l = this%ul(k)
    edge_r = this%ur(k)

    if ( (u2 - u1)*(u1 - u0) <= 0.0) then
      ! Flatten extremum
      edge_l = u1
      edge_r = u1
    else
      expr1 = 3.0 * (edge_r - edge_l) * ( (u1 - edge_l) + (u1 - edge_r))
      expr2 = (edge_r - edge_l) * (edge_r - edge_l)
      if ( expr1 > expr2 ) then
        ! Place extremum at right edge of cell by adjusting left edge value
        edge_l = u1 + 2.0 * ( u1 - edge_r )
        edge_l = max( min( edge_l, max(u0, u1) ), min(u0, u1) ) ! In case of round off
      elseif ( expr1 < -expr2 ) then
        ! Place extremum at left edge of cell by adjusting right edge value
        edge_r = u1 + 2.0 * ( u1 - edge_l )
        edge_r = max( min( edge_r, max(u2, u1) ), min(u2, u1) ) ! In case of round off
      endif
    endif
    ! This checks that the difference in edge values is representable
    ! and avoids overshoot problems due to round off.
    !### The 1.e-60 needs to have units of [A], so this dimensionally inconsistent.
    if ( abs( edge_r - edge_l )<max(1.e-60,epsilon(u1)*abs(u1)) ) then
      edge_l = u1
      edge_r = u1
    endif

    this%ul(k) = edge_l
    this%ur(k) = edge_r

  enddo ! end loop on interior cells


  ! After the limiter, are ur and ul bounded???? -AJA

  ! Store mean
  do k = 1, n
    this%u_mean(k) = u(k)
  enddo

end subroutine reconstruct

!> Checks the PPM_hybgen reconstruction for consistency
logical function check_reconstruction(this, h, u)
  class(PPM_hybgen), intent(in) :: this !< This reconstruction
  real,              intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,              intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: k

  check_reconstruction = .false.

  ! Simply checks the internal copy of "u" is exactly equal to "u"
  do k = 1, this%n
    if ( abs( this%u_mean(k) - u(k) ) > 0. ) check_reconstruction = .true.
  enddo

  ! If (u - ul) has the opposite sign from (ur - u), then this cell has an interior extremum
  do k = 1, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ur(k) - this%u_mean(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! The following consistency checks would fail for this implementation of PPM CW,
  ! due to round off in the final limiter violating the monotonicity of edge values,
  ! but actually passes due to the second pass of the limiters with explicit bounding.
  ! i.e. This implementation cheats!

  ! Check bounding of right edges, w.r.t. the cell means
  do K = 1, this%n-1
    if ( ( this%ur(k) - this%u_mean(k) ) * ( this%u_mean(k+1) - this%ur(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check bounding of left edges, w.r.t. the cell means
  do K = 2, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ul(k) - this%u_mean(k-1) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check bounding of right edges, w.r.t. this cell mean and the next cell left edge
  do K = 1, this%n-1
    if ( ( this%ur(k) - this%u_mean(k) ) * ( this%ul(k+1) - this%ur(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check bounding of left edges, w.r.t. this cell mean and the previous cell right edge
  do K = 2, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ul(k) - this%ur(k-1) ) < 0. ) check_reconstruction = .true.
  enddo

end function check_reconstruction

!> Runs PPM_hybgen reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(PPM_hybgen), intent(inout) :: this    !< This reconstruction
  logical,           intent(in)    :: verbose !< True, if verbose
  integer,           intent(in)    :: stdout  !< I/O channel for stdout
  integer,           intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables
  real, allocatable :: ul(:), ur(:), um(:) ! test values [A]
  real, allocatable :: ull(:), urr(:) ! test values [A]
  type(testing) :: test ! convenience functions
  integer :: k

  call test%set( stdout=stdout ) ! Sets the stdout channel in test
  call test%set( stderr=stderr ) ! Sets the stderr channel in test
  call test%set( verbose=verbose ) ! Sets the verbosity flag in test

  if (verbose) write(stdout,'(a)') 'PPM_hybgen:unit_tests testing with linear fn'

  call this%init(5)
  call test%test( this%n /= 5, 'Setting number of levels')
  allocate( um(5), ul(5), ur(5), ull(5), urr(5) )

  ! Straight line, f(x) = x , or  f(K) = 2*K
  call this%reconstruct( (/2.,2.,2.,2.,2./), (/1.,4.,7.,10.,13./) )
  call test%real_arr(5, this%u_mean, (/1.,4.,7.,10.,13./), 'Setting cell values')
  !   Without PLM extrapolation we get l(2)=2 and r(4)=12 due to PLM=0 in boundary cells. -AJA
  call test%real_arr(5, this%ul, (/1.,1.,5.5,8.5,13./), 'Left edge values')
  call test%real_arr(5, this%ur, (/1.,5.5,8.5,13.,13./), 'Right edge values')

  do k = 1, 5
    ul(k) = this%f(k, 0.)
    um(k) = this%f(k, 0.5)
    ur(k) = this%f(k, 1.)
  enddo
  call test%real_arr(5, ul, this%ul, 'Evaluation on left edge')
  call test%real_arr(5, um, (/1.,4.375,7.,9.625,13./), 'Evaluation in center')
  call test%real_arr(5, ur, this%ur, 'Evaluation on right edge')

  do k = 1, 5
    ul(k) = this%dfdx(k, 0.)
    um(k) = this%dfdx(k, 0.5)
    ur(k) = this%dfdx(k, 1.)
  enddo
  ! Most of these values are affected by the PLM boundary cells
  call test%real_arr(5, ul, (/0.,0.,3.,9.,0./), 'dfdx on left edge')
  call test%real_arr(5, um, (/0.,4.5,3.,4.5,0./), 'dfdx in center')
  call test%real_arr(5, ur, (/0.,9.,3.,0.,0./), 'dfdx on right edge')

  do k = 1, 5
    um(k) = this%average(k, 0.5, 0.75) ! Average from x=0.25 to 0.75 in each cell
  enddo
  ! Most of these values are affected by the PLM boundary cells
  call test%real_arr(5, um, (/1.,4.84375,7.375,10.28125,13./), 'Return interval average')

  if (verbose) write(stdout,'(a)') 'PPM_hybgen:unit_tests testing with parabola'

  ! x = 2 i   i=0 at origin
  ! f(x) = 3/4 x^2    = (2 i)^2
  ! f[i] = 3/4 ( 2 i - 1 )^2 on centers
  ! f[I] = 3/4 ( 2 I )^2 on edges
  ! f[i] = 1/8 [ x^3 ] for means
  ! edges:        0,  1, 12, 27, 48, 75
  ! means:          1,  7, 19, 37, 61
  ! cengters:      0.75, 6.75, 18.75, 36.75, 60.75
  call this%reconstruct( (/2.,2.,2.,2.,2./), (/1.,7.,19.,37.,61./) )
  do k = 1, 5
    ul(k) = this%f(k, 0.)
    um(k) = this%f(k, 0.5)
    ur(k) = this%f(k, 1.)
  enddo
  call test%real_arr(5, ul, (/1.,1.,12.,27.,61./), 'Return left edge')
  call test%real_arr(5, um, (/1.,7.25,18.75,34.5,61./), 'Return center')
  call test%real_arr(5, ur, (/1.,12.,27.,57.,61./), 'Return right edge')

  ! x = 3 i   i=0 at origin
  ! f(x) = x^2 / 3   = 3 i^2
  ! f[i] = [ ( 3 i )^3 - ( 3 i - 3 )^3 ]    i=1,2,3,4,5
  ! means:   1, 7, 19, 37, 61
  ! edges:  0, 3, 12, 27, 48, 75
  call this%reconstruct( (/3.,3.,3.,3.,3./), (/1.,7.,19.,37.,61./) )
  do k = 1, 5
    ul(k) = this%f(k, 0.)
    um(k) = this%f(k, 0.5)
    ur(k) = this%f(k, 1.)
  enddo
  call test%real_arr(5, ul, (/1.,1.,12.,27.,61./), 'Return left edge')
  call test%real_arr(5, ur, (/1.,12.,27.,57.,61./), 'Return right edge')

  call this%destroy()
  deallocate( um, ul, ur, ull, urr )

  unit_tests = test%summarize('PPM_hybgen:unit_tests')

end function unit_tests

!> \namespace recon1d_ppm_hybgen
!!

end module Recon1d_PPM_hybgen
