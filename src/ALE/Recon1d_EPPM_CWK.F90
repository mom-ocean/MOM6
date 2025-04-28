!> Piecewise Parabolic Method 1D reconstruction in model index space with linear
!! extrapolation for first and last cells
!!
!! This implementation of PPM follows Colella and Woodward, 1984, using uniform thickness
!! and with cells resorting to PCM for local extrema. First and last cells use a PLM
!! representation with slope set by matching the edge of the first interior cell.
module Recon1d_EPPM_CWK

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_type, only : Recon1d, testing
use Recon1d_PPM_CWK, only : PPM_CWK

implicit none ; private

public EPPM_CWK, testing

!> PPM reconstruction in index space (no grid dependence) with linear extrapolation
!! for first and last cells.
!!
!! Implemented by extending recon1d_ppm_cwk.
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_ppm_cwk.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_ppm_cwk.average()
!! - f()                    -> recon1d_ppm_cwk.f()
!! - dfdx()                 -> recon1d_ppm_cwk.dfdx()
!! - check_reconstruction() -> recon1d_ppm_cwk.check_reconstruction()
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_ppm_cwk.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> recon1d_ppm_cwk.init()
!! - reconstruct_parent()   -> recon1d_ppm_cwk.reconstruct()
type, extends (PPM_CWK) :: EPPM_CWK

contains
  !> Implementation of the EPPM_CWK reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of unit tests for the EPPM_CWK reconstruction
  procedure :: unit_tests => unit_tests

end type EPPM_CWK

contains

!> Calculate a 1D EPPM_CWK reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(EPPM_CWK), intent(inout) :: this !< This reconstruction
  real,            intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,            intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  real :: dul, dur ! Left and right cell PLM slopes [A]
  real :: u0, u1, u2 ! Far left, left, and right cell values [A]
  real :: edge ! Edge value between cell k-1 and k [A]
  real :: u_min, u_max ! Minimum and maximum value across edge [A]
  real :: a6 ! Colella and Woodward curvature [A]
  real :: du ! Difference between edges across cell [A]
  real :: slp(this%n) ! PLM slope [A]
  real, parameter :: one_sixth = 1. / 6. ! 1/6 [nondim]
  integer :: k, n

  n = this%n

  call this%reconstruct_parent( h, u )

  ! Extrapolate in first cell
  this%ur(1) = this%ul(2) ! Assume ur=ul on right edge
  this%ul(1) = u(1) + ( u(1) - this%ur(1) ) ! Linearly extrapolat across cell

  ! Extrapolate in last cell
  this%ul(n) = this%ur(n-1) ! Assume ul=ur on left edge
  this%ur(n) = u(n) + ( u(n) - this%ul(n) ) ! Linearly extrapolat across cell

end subroutine reconstruct

!> Runs EPPM_CWK reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(EPPM_CWK), intent(inout) :: this    !< This reconstruction
  logical,         intent(in)    :: verbose !< True, if verbose
  integer,         intent(in)    :: stdout  !< I/O channel for stdout
  integer,         intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables
  real, allocatable :: ul(:), ur(:), um(:) ! test values [A]
  real, allocatable :: ull(:), urr(:) ! test values [A]
  type(testing) :: test ! convenience functions
  integer :: k

  call test%set( stdout=stdout ) ! Sets the stdout channel in test
  call test%set( stderr=stderr ) ! Sets the stderr channel in test
  call test%set( verbose=verbose ) ! Sets the verbosity flag in test

  if (verbose) write(stdout,'(a)') 'EPPM_CWK:unit_tests testing with linear fn'

  call this%init(5)
  call test%test( this%n /= 5, 'Setting number of levels')
  allocate( um(5), ul(5), ur(5), ull(5), urr(5) )

  ! Straight line, f(x) = x , or  f(K) = 2*K
  call this%reconstruct( (/2.,2.,2.,2.,2./), (/1.,4.,7.,10.,13./) )
  call test%real_arr(5, this%u_mean, (/1.,4.,7.,10.,13./), 'Setting cell values')
  call test%real_arr(5, this%ul, (/-0.5,2.5,5.5,8.5,11.5/), 'Left edge values')
  call test%real_arr(5, this%ur, (/2.5,5.5,8.5,11.5,14.5/), 'Right edge values')

  do k = 1, 5
    ul(k) = this%f(k, 0.)
    um(k) = this%f(k, 0.5)
    ur(k) = this%f(k, 1.)
  enddo
  call test%real_arr(5, ul, this%ul, 'Evaluation on left edge')
  call test%real_arr(5, um, (/1.,4.,7.,10.,13./), 'Evaluation in center')
  call test%real_arr(5, ur, this%ur, 'Evaluation on right edge')

  do k = 1, 5
    ul(k) = this%dfdx(k, 0.)
    um(k) = this%dfdx(k, 0.5)
    ur(k) = this%dfdx(k, 1.)
  enddo
  ! Most of these values are affected by the PLM boundary cells
  call test%real_arr(5, ul, (/3.,3.,3.,3.,3./), 'dfdx on left edge')
  call test%real_arr(5, um, (/3.,3.,3.,3.,3./), 'dfdx in center')
  call test%real_arr(5, ur, (/3.,3.,3.,3.,3./), 'dfdx on right edge')

  do k = 1, 5
    um(k) = this%average(k, 0.5, 0.75) ! Average from x=0.25 to 0.75 in each cell
  enddo
  ! Most of these values are affected by the PLM boundary cells
  call test%real_arr(5, um, (/1.375,4.375,7.375,10.375,13.375/), 'Return interval average')

  if (verbose) write(stdout,'(a)') 'EPPM_CWK:unit_tests testing with parabola'

  ! x = 2 i   i=0 at origin
  ! f(x) = 3/4 x^2    = (2 i)^2
  ! f[i] = 3/4 ( 2 i - 1 )^2 on centers
  ! f[I] = 3/4 ( 2 I )^2 on edges
  ! f[i] = 1/8 [ x^3 ] for means
  ! edges:        0,  1, 12, 27, 48, 75
  ! means:          1,  7, 19, 37, 61
  ! centers:      0.75, 6.75, 18.75, 36.75, 60.75
  call this%reconstruct( (/2.,2.,2.,2.,2./), (/1.,7.,19.,37.,61./) )
  do k = 1, 5
    ul(k) = this%f(k, 0.)
    um(k) = this%f(k, 0.5)
    ur(k) = this%f(k, 1.)
  enddo
  call test%real_arr(5, ul, (/-1.,3.,12.,27.,48./), 'Return left edge')
  call test%real_arr(5, um, (/1.,6.75,18.75,36.75,61./), 'Return center')
  call test%real_arr(5, ur, (/3.,12.,27.,48.,74./), 'Return right edge')

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
  call test%real_arr(5, ul, (/-1.,3.,12.,27.,48./), 'Return left edge')
  call test%real_arr(5, um, (/1.,6.75,18.75,36.75,61./), 'Return center')
  call test%real_arr(5, ur, (/3.,12.,27.,48.,74./), 'Return right edge')

  call this%destroy()
  deallocate( um, ul, ur, ull, urr )

  unit_tests = test%summarize('EPPM_CWK:unit_tests')

end function unit_tests

!> \namespace recon1d_eppm_cwk
!!

end module Recon1d_EPPM_CWK
