!> Piecewise Parabolic Method 1D reconstruction with h4 interpolation for edges (2018 version)
!!
!! This implementation of PPM follows White and Adcroft 2008 \cite white2008, with cells
!! resorting to PCM for extrema including first and last cells in column.
!! This scheme differs from Colella and Woodward, 1984 \cite colella1984, in the method
!! of first estimating the fourth-order accurate edge values.
!! This uses numerical expressions that predate a 2019 refactoring.
!! The first and last cells are always limited to PCM.
module Recon1d_PPM_H4_2018

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_PPM_H4_2019, only : PPM_H4_2019, testing
use regrid_edge_values, only : bound_edge_values, check_discontinuous_edge_values
use regrid_solvers, only :  solve_linear_system

implicit none ; private

public PPM_H4_2018, testing

!> PPM reconstruction following White and Adcroft, 2008
!!
!! Implemented by extending recon1d_ppm_h4_2019.
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_ppm_h4_2019.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_ppm_h4_2019.average()
!! - f()                    -> recon1d_ppm_h4_2019.f()
!! - dfdx()                 -> recon1d_ppm_h4_2019.dfdx()
!! - check_reconstruction() -> recon1d_ppm_h4_2019.check_reconstruction()
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_ppm_h4_2019.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> recon1d_ppm_h4_2019.init()
!! - reconstruct_parent()   -> recon1d_ppm_h4_2019.reconstruct()
type, extends (PPM_H4_2019) :: PPM_H4_2018

contains
  !> Implementation of the PPM_H4_2018 reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of unit tests for the PPM_H4_2018 reconstruction
  procedure :: unit_tests => unit_tests

end type PPM_H4_2018

contains

!> Calculate a 1D PPM_H4_2018 reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(PPM_H4_2018), intent(inout) :: this !< This reconstruction
  real,               intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,               intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  real :: slp ! The PLM slopes (difference across cell) [A]
  real :: u_l, u_r, u_c ! Left, right, and center values [A]
  real :: h_l, h_c, h_r ! Thickness of left, center and right cells [H]
  real :: h0, h1, h2, h3        ! temporary thicknesses [H]
  real :: h_min                 ! A minimal cell width [H]
  real :: f1                    ! An auxiliary variable [H]
  real :: f2                    ! An auxiliary variable [A H]
  real :: f3                    ! An auxiliary variable [H-1]
  real :: et1, et2, et3         ! terms the expression for edge values [A H]
  real :: dx                    ! Difference of successive values of x [H]
  real :: f                     ! value of polynomial at x in arbitrary units [A]
  real :: edge_l, edge_r        ! Edge values (left and right) [A]
  real :: expr1, expr2          ! Temporary expressions [A2]
  real, parameter :: hMinFrac = 1.e-5  !< A minimum fraction for min(h)/sum(h) [nondim]
  real, dimension(5) :: x ! Coordinate system with 0 at edges [H]
  real :: edge_values(this%n,2) ! Edge values [A]
  real :: ppoly_coef(this%n,3) ! Polynomial coefficients [A]
  real, dimension(4,4) :: A ! Differences in successive positions raised to various powers,
                            ! in units that vary with the second (j) index as [H^j]
  real, dimension(4)   :: B ! The right hand side of the system to solve for C [A H]
  real, dimension(4)   :: C ! The coefficients of a fit polynomial in units that vary
                            ! with the index (j) as [A H^(j-1)]
  integer :: k, n, j

  n = this%n

  ! Loop on interior cells
  do K = 3, n-1

    h0 = h(k-2)
    h1 = h(k-1)
    h2 = h(k)
    h3 = h(k+1)

    ! Avoid singularities when consecutive pairs of h vanish
    if (h0+h1==0.0 .or. h1+h2==0.0 .or. h2+h3==0.0) then
      h_min = hMinFrac*max( this%h_neglect, h0+h1+h2+h3 )
      h0 = max( h_min, h(k-2) )
      h1 = max( h_min, h(k-1) )
      h2 = max( h_min, h(k) )
      h3 = max( h_min, h(k+1) )
    endif

    f1 = (h0+h1) * (h2+h3) / (h1+h2)
    f2 = h2 * u(k-1) + h1 * u(k)
    f3 = 1.0 / (h0+h1+h2) + 1.0 / (h1+h2+h3)
    et1 = f1 * f2 * f3
    et2 = ( h2 * (h2+h3) / ( (h0+h1+h2)*(h0+h1) ) ) * &
          ((h0+2.0*h1) * u(k-1) - h1 * u(k-2))
    et3 = ( h1 * (h0+h1) / ( (h1+h2+h3)*(h2+h3) ) ) * &
          ((2.0*h2+h3) * u(k) - h2 * u(k+1))
    edge_values(k,1) = (et1 + et2 + et3) / ( h0 + h1 + h2 + h3)
    edge_values(k-1,2) = edge_values(k,1)

  enddo ! end loop on interior cells

  ! Determine first two edge values
  h_min = max( this%h_neglect, hMinFrac*sum(h(1:4)) )
  x(1) = 0.0
  do k = 1,4
    dx = max(h_min, h(k) )
    x(k+1) = x(k) + dx
    do j = 1,4 ; A(k,j) = ( (x(k+1)**j) - (x(k)**j) ) / real(j) ; enddo
    B(k) = u(k) * dx
  enddo

  call solve_linear_system( A, B, C, 4 )

  ! Set the edge values of the first cell
  f = 0.0
  do k = 1, 4
    f = f + C(k) * ( x(1)**(k-1) )
  enddo
  edge_values(1,1) = f
  f = 0.0
  do k = 1, 4
    f = f + C(k) * ( x(2)**(k-1) )
  enddo
  edge_values(1,2) = f
  edge_values(2,1) = edge_values(1,2)

  ! Determine two edge values of the last cell
  h_min = max( this%h_neglect, hMinFrac*sum(h(n-3:n)) )
  x(1) = 0.0
  do k = 1,4
    dx = max(h_min, h(n-4+k) )
    x(k+1) = x(k) + dx
    do j = 1,4 ; A(k,j) = ( (x(k+1)**j) - (x(k)**j) ) / real(j) ; enddo
    B(k) = u(n-4+k) * dx
  enddo

  call solve_linear_system( A, B, C, 4 )

  ! Set the last and second to last edge values
  f = 0.0
  do k = 1, 4
    f = f + C(k) * ( x(5)**(k-1) )
  enddo
  edge_values(n,2) = f
  f = 0.0
  do k = 1, 4
    f = f + C(k) * ( x(4)**(k-1) )
  enddo
  edge_values(n,1) = f
  edge_values(n-1,2) = edge_values(n,1)

  ! Bound edge values
  call bound_edge_values( n, h, u, edge_values, this%h_neglect, answer_date=20180101 )

  ! Make discontinuous edge values monotonic
  call check_discontinuous_edge_values( n, u, edge_values )

  ! Loop on interior cells to apply the standard
  ! PPM limiter (Colella & Woodward, JCP 84)
  do k = 2,n-1

    ! Get cell averages
    u_l = u(k-1)
    u_c = u(k)
    u_r = u(k+1)

    edge_l = edge_values(k,1)
    edge_r = edge_values(k,2)

    if ( (u_r - u_c)*(u_c - u_l) <= 0.0) then
      ! Flatten extremum
      edge_l = u_c
      edge_r = u_c
    else
      expr1 = 3.0 * (edge_r - edge_l) * ( (u_c - edge_l) + (u_c - edge_r))
      expr2 = (edge_r - edge_l) * (edge_r - edge_l)
      if ( expr1 > expr2 ) then
        ! Place extremum at right edge of cell by adjusting left edge value
        edge_l = u_c + 2.0 * ( u_c - edge_r )
        edge_l = max( min( edge_l, max(u_l, u_c) ), min(u_l, u_c) ) ! In case of round off
      elseif ( expr1 < -expr2 ) then
        ! Place extremum at left edge of cell by adjusting right edge value
        edge_r = u_c + 2.0 * ( u_c - edge_l )
        edge_r = max( min( edge_r, max(u_r, u_c) ), min(u_r, u_c) ) ! In case of round off
      endif
    endif
    ! This checks that the difference in edge values is representable
    ! and avoids overshoot problems due to round off.
    !### The 1.e-60 needs to have units of [A], so this dimensionally inconsistent.
    if ( abs( edge_r - edge_l )<max(1.e-60,epsilon(u_c)*abs(u_c)) ) then
      edge_l = u_c
      edge_r = u_c
    endif

    edge_values(k,1) = edge_l
    edge_values(k,2) = edge_r

  enddo ! end loop on interior cells

  ! PCM within boundary cells
  edge_values(1,:) = u(1)
  edge_values(n,:) = u(n)

  ! Store reconstruction
  do k = 1, n
    this%u_mean(k) = u(k)
    this%ul(k) = edge_values(k,1)
    this%ur(k) = edge_values(k,2)
  enddo

end subroutine reconstruct

!> Runs PPM_H4_2018 reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(PPM_H4_2018), intent(inout) :: this    !< This reconstruction
  logical,            intent(in)    :: verbose !< True, if verbose
  integer,            intent(in)    :: stdout  !< I/O channel for stdout
  integer,            intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables
  real, allocatable :: ul(:), ur(:), um(:) ! test values [A]
  real, allocatable :: ull(:), urr(:) ! test values [A]
  type(testing) :: test ! convenience functions
  integer :: k

  call test%set( stdout=stdout ) ! Sets the stdout channel in test
  call test%set( stderr=stderr ) ! Sets the stderr channel in test
  call test%set( verbose=verbose ) ! Sets the verbosity flag in test

  if (verbose) write(stdout,'(a)') 'PPM_H4_2018:unit_tests testing with linear fn'

  call this%init(5)
  call test%test( this%n /= 5, 'Setting number of levels')
  allocate( um(5), ul(5), ur(5), ull(5), urr(5) )

  ! Straight line, f(x) = x , or  f(K) = 2*K
  call this%reconstruct( (/2.,2.,2.,2.,2./), (/1.,3.,5.,7.,9./) )
  call test%real_arr(5, this%u_mean, (/1.,3.,5.,7.,9./), 'Setting cell values')
  call test%real_arr(5, this%ul, (/1.,2.,4.,6.,9./), 'Left edge values', robits=2)
  call test%real_arr(5, this%ur, (/1.,4.,6.,8.,9./), 'Right edge values', robits=1)
  do k = 1, 5
    um(k) = this%u_mean(k)
  enddo
  call test%real_arr(5, um, (/1.,3.,5.,7.,9./), 'Return cell mean')

  do k = 1, 5
    ul(k) = this%f(k, 0.)
    um(k) = this%f(k, 0.5)
    ur(k) = this%f(k, 1.)
  enddo
  call test%real_arr(5, ul, this%ul, 'Evaluation on left edge')
  call test%real_arr(5, um, (/1.,3.,5.,7.,9./), 'Evaluation in center')
  call test%real_arr(5, ur, this%ur, 'Evaluation on right edge')

  do k = 1, 5
    ul(k) = this%dfdx(k, 0.)
    um(k) = this%dfdx(k, 0.5)
    ur(k) = this%dfdx(k, 1.)
  enddo
  call test%real_arr(5, ul, (/0.,2.,2.,2.,0./), 'dfdx on left edge', robits=4)
  call test%real_arr(5, um, (/0.,2.,2.,2.,0./), 'dfdx in center', robits=2)
  call test%real_arr(5, ur, (/0.,2.,2.,2.,0./), 'dfdx on right edge', robits=6)

  do k = 1, 5
    um(k) = this%average(k, 0.5, 0.75) ! Average from x=0.25 to 0.75 in each cell
  enddo
  call test%real_arr(5, um, (/1.,3.25,5.25,7.25,9./), 'Return interval average')

  if (verbose) write(stdout,'(a)') 'PPM_H4_2018:unit_tests testing with parabola'

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
  call test%real_arr(5, ul, (/1.,3.,12.,27.,61./), 'Return left edge', robits=2)
  call test%real_arr(5, ur, (/1.,12.,27.,48.,61./), 'Return right edge', robits=1)

  call this%destroy()
  deallocate( um, ul, ur, ull, urr )

  unit_tests = test%summarize('PPM_H4_2018:unit_tests')

end function unit_tests

!> \namespace recon1d_ppm_h4_2018
!!

end module Recon1d_PPM_H4_2018
