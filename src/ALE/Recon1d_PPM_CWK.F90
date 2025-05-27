!> Piecewise Parabolic Method 1D reconstruction in model index space
!!
!! This implementation of PPM follows Colella and Woodward, 1984, using uniform thickness
!! and with cells resorting to PCM for local extrema including the first and last cells.
!!
!! "Fourth order" estimates of edge values use PLM also calculated in index space
!! (i.e. with no grid dependence). First and last PLM slopes are extrapolated.
!! Limiting follows Colella and Woodward thereafter. The high accuracy of this scheme is
!! realized only when the grid-spacing is exactly uniform. This scheme deviates from CW84
!! when the grid spacing is variable.
module Recon1d_PPM_CWK

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_type, only : Recon1d, testing
use Recon1d_PLM_CWK, only : PLM_CWK

implicit none ; private

public PPM_CWK, testing

!> PPM reconstruction in index space (no grid dependence).
!!
!! The source for the methods ultimately used by this class are:
!! - init()                    *locally defined
!! - reconstruct()             *locally defined
!! - average()                 *locally defined
!! - f()                       *locally defined
!! - dfdx()                    *locally defined
!! - check_reconstruction()    *locally defined
!! - unit_tests()              *locally defined
!! - destroy()                 *locally defined
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (Recon1d) :: PPM_CWK

  real, allocatable :: ul(:) !< Left edge value [A]
  real, allocatable :: ur(:) !< Right edge value [A]
  type(PLM_CWK) :: PLM !< The PLM reconstruction used to estimate edge values

contains
  !> Implementation of the PPM_CWK initialization
  procedure :: init => init
  !> Implementation of the PPM_CWK reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of the PPM_CWK average over an interval [A]
  procedure :: average => average
  !> Implementation of evaluating the PPM_CWK reconstruction at a point [A]
  procedure :: f => f
  !> Implementation of the derivative of the PPM_CWK reconstruction at a point [A]
  procedure :: dfdx => dfdx
  !> Implementation of deallocation for PPM_CWK
  procedure :: destroy => destroy
  !> Implementation of check reconstruction for the PPM_CWK reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PPM_CWK reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to init()
  procedure :: init_parent => init
  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type PPM_CWK

contains

!> Initialize a 1D PPM_CWK reconstruction for n cells
subroutine init(this, n, h_neglect, check)
  class(PPM_CWK),     intent(out) :: this      !< This reconstruction
  integer,            intent(in)  :: n         !< Number of cells in this column
  real, optional,     intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  logical, optional,  intent(in)  :: check     !< If true, enable some consistency checking

  this%n = n

  allocate( this%u_mean(n) )
  allocate( this%ul(n) )
  allocate( this%ur(n) )

  ! This incurs an extra store of u_mean but by using PCM_CW
  ! we avoid duplicating and testing more code
  call this%PLM%init( n, h_neglect=h_neglect, check=check )

  this%h_neglect = tiny( this%u_mean(1) )
  if (present(h_neglect)) this%h_neglect = h_neglect
  this%check = .false.
  if (present(check)) this%check = check

end subroutine init

!> Calculate a 1D PPM_CWK reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(PPM_CWK), intent(inout) :: this !< This reconstruction
  real,           intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in)    :: u(*) !< Cell mean values [A]
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

  ! First populate the PLM (k-space) reconstructions
  call this%PLM%reconstruct( h, u )
  do k = 1, n
    slp(k) = this%PLM%ur(k) - this%PLM%ul(k)
  enddo
  ! Extrapolate from interior for boundary PLM slopes
  ! Note: this is not conventional but helps retain accuracy near top/bottom
  ! boudaries and reduces the adverse influence of the boudnaries int he interior
  ! reconstructions. The final PPM reconstruction is still bounded to PCM.
  slp(1) = 2.0 * ( this%PLM%ul(2) - u(1) )
  slp(n) = 2.0 * ( u(n) - this%PLM%ur(n-1) )

  do K = 2, n ! K=2 is interface between cells 1 and 2
    dul = slp(k-1)
    dur = slp(k)
    u2 = u(k)
    u1 = u(k-1)
    edge = 0.5 * ( u1 + u2 ) + one_sixth * ( dul - dur )  ! Eq. 1.6 with uniform h
    u_min = min( u1, u2 )
    u_max = max( u1, u2 )
    edge = max( min( edge, u_max), u_min ) ! Unclear if we need this bounding in the interior
    this%ur(k-1) = edge
    this%ul(k) = edge
  enddo
  this%ul(1) = u(1) ! PCM
  this%ur(1) = u(1) ! PCM
  this%ur(n) = u(n) ! PCM
  this%ul(n) = u(n) ! PCM

  do K = 2, n ! K=2 is interface between cells 1 and 2
    u0 = u(k-1)
    u1 = u(k)
    u2 = u(k+1)
    a6 = 3.0 * ( ( u1 - this%ul(k) ) + ( u1 - this%ur(k) ) )
    du = this%ur(k) - this%ul(k)
    if ( ( u2 - u1 ) * ( u1 - u0 ) <- 0.0 ) then ! Large scale extrema
      this%ul(k) = u1
      this%ur(k) = u1
    elseif ( du * a6 > du * du ) then ! Extrema on right
      edge = u1 + 2.0 * ( u1 - this%ur(k) )
    ! u_min = min( u0, u1 )
    ! u_max = max( u0, u1 )
    ! edge = max( min( edge, u_max), u_min )
      this%ul(k) = edge
    elseif ( du * a6 < - du * du ) then ! Extrema on left
      edge = u1 + 2.0 * ( u1 - this%ul(k) )
    ! u_min = min( u1, u2 )
    ! u_max = max( u1, u2 )
    ! edge = max( min( edge, u_max), u_min )
      this%ur(k) = edge
    endif
  enddo

  ! After the limiter, are ur and ul bounded???? -AJA

  ! Store mean
  do k = 1, n
    this%u_mean(k) = u(k)
  enddo

end subroutine reconstruct

!> Value of PPM_CWK reconstruction at a point in cell k [A]
real function f(this, k, x)
  class(PPM_CWK), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: x    !< Non-dimensional position within element [nondim]
  real :: xc ! Bounded version of x [nondim]
  real :: du ! Difference across cell [A]
  real :: a6 ! Collela and Woordward curvature parameter [A]
  real :: u_a, u_b ! Two estimate of f [A]
  real :: lmx ! 1 - x [nondim]
  real :: wb ! Weight based on x [nondim]

  du = this%ur(k) - this%ul(k)
  a6 = 3.0 * ( ( this%u_mean(k) - this%ul(k) ) + ( this%u_mean(k) - this%ur(k) ) )
  xc = max( 0., min( 1., x ) )
  lmx = 1.0 - xc

  ! This expression for u_a can overshoot u_r but is good for x<<1
  u_a = this%ul(k) + xc * ( du + a6 * lmx )
  ! This expression for u_b can overshoot u_l but is good for 1-x<<1
  u_b = this%ur(k) + lmx * ( - du + a6 * xc )

  ! Since u_a and u_b are both side-bounded, using weights=0 or 1 will preserve uniformity
  wb = 0.5 + sign(0.5, xc - 0.5 ) ! = 1 @ x=0, = 0 @ x=1
  f = ( ( 1. - wb ) * u_a ) + ( wb * u_b )

end function f

!> Derivative of PPM_CWK reconstruction at a point in cell k [A]
real function dfdx(this, k, x)
  class(PPM_CWK), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: x    !< Non-dimensional position within element [nondim]
  real :: xc ! Bounded version of x [nondim]
  real :: du ! Difference across cell [A]
  real :: a6 ! Collela and Woordward curvature parameter [A]

  du = this%ur(k) - this%ul(k)
  a6 = 3.0 * ( ( this%u_mean(k) - this%ul(k) ) + ( this%u_mean(k) - this%ur(k) ) )
  xc = max( 0., min( 1., x ) )

  dfdx = du + a6 * ( 2.0 * xc - 1.0 )

end function dfdx

!> Average between xa and xb for cell k of a 1D PPM reconstruction [A]
real function average(this, k, xa, xb)
  class(PPM_CWK), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real,           intent(in) :: xb   !< End of averaging interval on element (0 to 1)
  real :: xapxb                      ! A sum of fracional positions [nondim]
  real :: mx, Ya, Yb, my             ! Various fractional positions [nondim]
  real :: u_a, u_b                   ! Values at xa and xb [A]
  real :: xa2pxb2,  xa2b2ab, Ya2b2ab ! Sums of squared fractional positions [nondim]
  real :: a_L, a_R, u_c, a_c         ! Values of the polynomial at various locations [A]

  mx = 0.5 * ( xa + xb )
  a_L = this%ul(k)
  a_R = this%ur(k)
  u_c = this%u_mean(k)
  a_c = 0.5 * ( ( u_c - a_L ) + ( u_c - a_R ) ) ! a_6 / 6
  if (mx<0.5) then
    ! This integration of the PPM reconstruction is expressed in distances from the left edge
    xa2b2ab = (xa * xa + xb * xb) + xa * xb
    average = a_L + ( ( a_R - a_L ) * mx &
                    + a_c * ( 3. * ( xb + xa ) - 2. * xa2b2ab ) )
  else
    ! This integration of the PPM reconstruction is expressed in distances from the right edge
    Ya = 1. - xa
    Yb = 1. - xb
    my = 0.5 * ( Ya + Yb )
    Ya2b2ab = (Ya * Ya + Yb * Yb) + Ya * Yb
    average = a_R + ( ( a_L - a_R ) * my &
                    + a_c * ( 3. * ( Yb + Ya ) - 2. * Ya2b2ab ) )
  endif

end function average

!> Deallocate the PPM_CWK reconstruction
subroutine destroy(this)
  class(PPM_CWK), intent(inout) :: this !< This reconstruction

  deallocate( this%u_mean, this%ul, this%ur )

end subroutine destroy

!> Checks the PPM_CWK reconstruction for consistency
logical function check_reconstruction(this, h, u)
  class(PPM_CWK), intent(in) :: this !< This reconstruction
  real,           intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in) :: u(*) !< Cell mean values [A]
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

!> Runs PPM_CWK reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(PPM_CWK), intent(inout) :: this    !< This reconstruction
  logical,        intent(in)    :: verbose !< True, if verbose
  integer,        intent(in)    :: stdout  !< I/O channel for stdout
  integer,        intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables
  real, allocatable :: ul(:), ur(:), um(:) ! test values [A]
  real, allocatable :: ull(:), urr(:) ! test values [A]
  type(testing) :: test ! convenience functions
  integer :: k

  call test%set( stdout=stdout ) ! Sets the stdout channel in test
  call test%set( stderr=stderr ) ! Sets the stderr channel in test
  call test%set( verbose=verbose ) ! Sets the verbosity flag in test

  if (verbose) write(stdout,'(a)') 'PPM_CWK:unit_tests testing with linear fn'

  call this%init(5)
  call test%test( this%n /= 5, 'Setting number of levels')
  allocate( um(5), ul(5), ur(5), ull(5), urr(5) )

  ! Straight line, f(x) = x , or  f(K) = 2*K
  call this%reconstruct( (/2.,2.,2.,2.,2./), (/1.,4.,7.,10.,13./) )
  call test%real_arr(5, this%u_mean, (/1.,4.,7.,10.,13./), 'Setting cell values')
  !   Without PLM extrapolation we get l(2)=2 and r(4)=12 due to PLM=0 in boundary cells. -AJA
  call test%real_arr(5, this%ul, (/1.,2.5,5.5,8.5,13./), 'Left edge values')
  call test%real_arr(5, this%ur, (/1.,5.5,8.5,11.5,13./), 'Right edge values')

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
  call test%real_arr(5, ul, (/0.,3.,3.,3.,0./), 'dfdx on left edge')
  call test%real_arr(5, um, (/0.,3.,3.,3.,0./), 'dfdx in center')
  call test%real_arr(5, ur, (/0.,3.,3.,3.,0./), 'dfdx on right edge')

  do k = 1, 5
    um(k) = this%average(k, 0.5, 0.75) ! Average from x=0.25 to 0.75 in each cell
  enddo
  ! Most of these values are affected by the PLM boundary cells
  call test%real_arr(5, um, (/1.,4.375,7.375,10.375,13./), 'Return interval average')

  if (verbose) write(stdout,'(a)') 'PPM_CWK:unit_tests testing with parabola'

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
  call test%real_arr(5, ul, (/1.,3.,12.,27.,61./), 'Return left edge')
  call test%real_arr(5, um, (/1.,6.75,18.75,36.75,61./), 'Return center')
  call test%real_arr(5, ur, (/1.,12.,27.,48.,61./), 'Return right edge')

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
  call test%real_arr(5, ul, (/1.,3.,12.,27.,61./), 'Return left edge')
  call test%real_arr(5, um, (/1.,6.75,18.75,36.75,61./), 'Return center')
  call test%real_arr(5, ur, (/1.,12.,27.,48.,61./), 'Return right edge')

  call this%destroy()
  deallocate( um, ul, ur, ull, urr )

  unit_tests = test%summarize('PPM_CWK:unit_tests')

end function unit_tests

!> \namespace recon1d_ppm_cwk
!!

end module Recon1d_PPM_CWK
