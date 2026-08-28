! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

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
!! - average()                 *locally defined but calls recon1d_ppm_cw.average()
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
  !> Implementation of the PPM_hybgen average over an interval [A]
  procedure :: average => average
  !> Implementation of check reconstruction for the PPM_hybgen reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PPM_hybgen reconstruction
  procedure :: unit_tests => unit_tests

end type PPM_hybgen

contains

!> Calculate a 1D PPM_hybgen reconstruction based on h(:) and u(:)
!!
!! First pass: hybgen_ppm_coefs() computes initial edge estimates with CW monotonicity.
!! Second pass: applies OM4-era bound_edge_values() and check_discontinuous_edge_values(),
!! then the standard CW PPM limiter (post-2018 expressions, answer_date=99991231).
!! This reproduces bit-for-bit the behavior of the old-style PPM_HYBGEN scheme.
subroutine reconstruct(this, h, u)
  class(PPM_hybgen), intent(inout) :: this !< This reconstruction
  real,              intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,              intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: k, n
  real :: ppoly_e(this%n, 2) ! PPM edge values [A]
  real :: u_l, u_c, u_r      ! Left, center, right cell averages [A]
  real :: edge_l, edge_r     ! Left and right edge values [A]
  real :: expr1, expr2       ! Temporary expressions [A2]

  n = this%n

  ! First pass: compute initial edge estimates using the hybgen algorithm with CW monotonicity
  call hybgen_ppm_coefs(u, h, ppoly_e, n, this%h_neglect)

  ! Second pass: apply OM4-era PPM limiters (post-2018 answers via answer_date=99991231)
  call bound_edge_values(n, h, u, ppoly_e, this%h_neglect, answer_date=99991231)
  call check_discontinuous_edge_values(n, u, ppoly_e)

  ! Apply the standard CW PPM limiter (Colella & Woodward, JCP 84) on interior cells
  do k = 2, n-1
    u_l = u(k-1) ; u_c = u(k) ; u_r = u(k+1)
    edge_l = ppoly_e(k,1) ; edge_r = ppoly_e(k,2)
    if ( (u_r - u_c)*(u_c - u_l) <= 0.0 ) then
      edge_l = u_c ; edge_r = u_c
    else
      expr1 = 3.0 * (edge_r - edge_l) * ( (u_c - edge_l) + (u_c - edge_r) )
      expr2 = (edge_r - edge_l) * (edge_r - edge_l)
      if ( expr1 > expr2 ) then
        edge_l = u_c + 2.0 * ( u_c - edge_r )
        edge_l = max( min( edge_l, max(u_l, u_c) ), min(u_l, u_c) )
      elseif ( expr1 < -expr2 ) then
        edge_r = u_c + 2.0 * ( u_c - edge_l )
        edge_r = max( min( edge_r, max(u_r, u_c) ), min(u_r, u_c) )
      endif
    endif
    !### The 1.e-60 needs to have units of [A], so this is dimensionally inconsistent.
    if ( abs( edge_r - edge_l ) < max(1.e-60, epsilon(u_c)*abs(u_c)) ) then
      edge_l = u_c ; edge_r = u_c
    endif
    ppoly_e(k,1) = edge_l ; ppoly_e(k,2) = edge_r
  enddo
  ! Boundary cells are PCM
  ppoly_e(1,:) = u(1) ; ppoly_e(n,:) = u(n)

  do k = 1, n
    this%ul(k) = ppoly_e(k, 1)
    this%ur(k) = ppoly_e(k, 2)
    this%u_mean(k) = u(k)
  enddo

end subroutine reconstruct

!> Average between xa and xb for cell k of a PPM_hybgen reconstruction [A]
!!
!! Calls the parent PPM_CW average and then clamps the result to [min(ul,ur), max(ul,ur)].
!! This replicates the force_bounds_in_subcell behavior of the equivalent old-style PPM_HYBGEN
!! scheme.
real function average(this, k, xa, xb)
  class(PPM_hybgen), intent(in) :: this !< This reconstruction
  integer,           intent(in) :: k    !< Cell number
  real,              intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real,              intent(in) :: xb   !< End of averaging interval on element (0 to 1)
  real :: u_lo, u_hi ! Bounds on the sub-cell average given by the edge values [A]

  average = this%PPM_CW%average(k, xa, xb)
  u_lo = min(this%ul(k), this%ur(k))
  u_hi = max(this%ul(k), this%ur(k))
  average = max(u_lo, min(u_hi, average))

end function average

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

! ============================================================================
! Private subroutines copied from phased-out modules to avoid dependencies.
! These reproduce bit-for-bit the results of the original functions they replace.
! ============================================================================

!> Set up edge values for PPM reconstruction using the hybgen (HYCOM) algorithm.
!!
!! Copied from MOM_hybgen_remap.hybgen_ppm_coefs().
!! Original code by Tim Campbell (MSU, 2002) and Alan Wallcraft (NRL, 2007).
subroutine hybgen_ppm_coefs(s, h_src, edges, nk, thin, PCM_lay)
  integer, intent(in)  :: nk        !< The number of input layers
  real,    intent(in)  :: s(nk)     !< The input scalar fields [A]
  real,    intent(in)  :: h_src(nk) !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: edges(nk,2) !< The PPM interpolation edge values [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(nk) !< If true for a layer, use PCM remapping

  real :: dp(nk) ! Input grid layer thicknesses, but with a minimum thickness given by thin [H ~> m or kg m-2]
  logical :: PCM_layer(nk) ! True for layers that should use PCM remapping
  real :: da        ! Difference between the unlimited scalar edge value estimates [A]
  real :: a6        ! Scalar field differences that are proportional to the curvature [A]
  real :: slk, srk  ! Differences between adjacent cell averages of scalars [A]
  real :: sck       ! Scalar differences across a cell [A]
  real :: as(nk)    ! Scalar field difference across each cell [A]
  real :: al(nk), ar(nk)   ! Scalar field at the left and right edges of a cell [A]
  real :: h112(nk+1), h122(nk+1)  ! Combinations of thicknesses [H ~> m or kg m-2]
  real :: I_h12(nk+1) ! Inverses of combinations of thicknesses [H-1 ~> m-1 or m2 kg-1]
  real :: h2_h123(nk)  ! A ratio of a layer thickness to the sum of 3 adjacent thicknesses [nondim]
  real :: I_h0123(nk)     ! Inverse of the sum of 4 adjacent thicknesses [H-1 ~> m-1 or m2 kg-1]
  real :: h01_h112(nk+1) ! A ratio of sums of adjacent thicknesses [nondim]
  real :: h23_h122(nk+1) ! A ratio of sums of adjacent thicknesses [nondim]
  integer :: k

  do k=1,nk ; dp(k) = max(h_src(k), thin) ; enddo

  if (present(PCM_lay)) then
    do k=1,nk ; PCM_layer(k) = (PCM_lay(k) .or. dp(k) <= thin) ; enddo
  else
    do k=1,nk ; PCM_layer(k) = (dp(k) <= thin) ; enddo
  endif

  do k=2,nk
    h112(K) = 2.*dp(k-1) + dp(k)
    h122(K) = dp(k-1) + 2.*dp(k)
    I_h12(K) = 1.0 / (dp(k-1) + dp(k))
  enddo
  do k=2,nk-1
    h2_h123(k) = dp(k) / (dp(k) + (dp(k-1)+dp(k+1)))
  enddo
  do K=3,nk-1
    I_h0123(K) = 1.0 / ((dp(k-2) + dp(k-1)) + (dp(k) + dp(k+1)))
    h01_h112(K) = (dp(k-2) + dp(k-1)) / (2.0*dp(k-1) + dp(k))
    h23_h122(K) = (dp(k) + dp(k+1))   / (dp(k-1) + 2.0*dp(k))
  enddo

    as(1) = 0.
    do k=2,nk-1
      if (PCM_layer(k)) then
        as(k) = 0.0
      else
        slk = s(k)-s(k-1)
        srk = s(k+1)-s(k)
        if (slk*srk > 0.) then
          sck = h2_h123(k)*( h112(K)*srk*I_h12(K+1) + h122(K+1)*slk*I_h12(K) )
          as(k) = sign(min(abs(2.0*slk), abs(sck), abs(2.0*srk)), sck)
        else
          as(k) = 0.
        endif
      endif
    enddo
    as(nk) = 0.
    al(1) = s(1)
    ar(1) = s(1)
    al(2) = s(1)
    do K=3,nk-1
      al(k) = (dp(k)*s(k-1) + dp(k-1)*s(k)) * I_h12(K) &
            + I_h0123(K)*( 2.*dp(k)*dp(k-1)*I_h12(K)*(s(k)-s(k-1)) * &
                           ( h01_h112(K) - h23_h122(K) ) &
                    + (dp(k)*as(k-1)*h23_h122(K) - dp(k-1)*as(k)*h01_h112(K)) )
      ar(k-1) = al(k)
    enddo
    ar(nk-1) = s(nk)
    al(nk)  = s(nk)
    ar(nk)  = s(nk)
    do k=2,nk-1
      if ((PCM_layer(k)) .or. ((s(k+1)-s(k))*(s(k)-s(k-1)) <= 0.)) then
        al(k) = s(k)
        ar(k) = s(k)
      else
        da = ar(k)-al(k)
        a6 = 6.0*s(k) - 3.0*(al(k)+ar(k))
        if (da*a6 > da*da) then
          al(k) = 3.0*s(k) - 2.0*ar(k)
        elseif (da*a6 < -da*da) then
          ar(k) = 3.0*s(k) - 2.0*al(k)
        endif
      endif
    enddo
    do k=1,nk
      edges(k,1) = al(k)
      edges(k,2) = ar(k)
    enddo

end subroutine hybgen_ppm_coefs

!> Bound edge values by the averages of the neighboring cells.
!!
!! Copied from regrid_edge_values.bound_edge_values().
subroutine bound_edge_values(N, h, u, edge_val, h_neglect, answer_date)
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< Cell widths [H]
  real, dimension(N),   intent(in)    :: u !< Cell averages [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Edge values [A]
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]
  integer,    optional, intent(in)    :: answer_date !< The vintage of the expressions to use

  real    :: sigma_l, sigma_c, sigma_r
  real    :: slope_x_h
  logical :: use_2018_answers
  integer :: k, km1, kp1

  use_2018_answers = .true. ; if (present(answer_date)) use_2018_answers = (answer_date < 20190101)

  do k = 1,N
    km1 = max(1,k-1) ; kp1 = min(k+1,N)
    slope_x_h = 0.0
    if (use_2018_answers) then
      sigma_l = 2.0 * ( u(k) - u(km1) ) / ( h(k) + h_neglect )
      sigma_c = 2.0 * ( u(kp1) - u(km1) ) / ( h(km1) + 2.0*h(k) + h(kp1) + h_neglect )
      sigma_r = 2.0 * ( u(kp1) - u(k) ) / ( h(k) + h_neglect )
      if ( (sigma_l * sigma_r) > 0.0 ) &
        slope_x_h = 0.5 * h(k) * sign( min(abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )
    elseif ( ((h(km1) + h(kp1)) + 2.0*h(k)) > 0.0 ) then
      sigma_l = ( u(k) - u(km1) )
      sigma_c = ( u(kp1) - u(km1) ) * ( h(k) / ((h(km1) + h(kp1)) + 2.0*h(k)) )
      sigma_r = ( u(kp1) - u(k) )
      if ( (sigma_l * sigma_r) > 0.0 ) &
        slope_x_h = sign( min(abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )
    endif
    if ( (u(km1)-edge_val(k,1)) * (edge_val(k,1)-u(k)) < 0.0 ) then
      edge_val(k,1) = u(k) - sign( min( abs(slope_x_h), abs(edge_val(k,1)-u(k)) ), slope_x_h )
    endif
    if ( (u(kp1)-edge_val(k,2)) * (edge_val(k,2)-u(k)) < 0.0 ) then
      edge_val(k,2) = u(k) + sign( min( abs(slope_x_h), abs(edge_val(k,2)-u(k)) ), slope_x_h )
    endif
    edge_val(k,1) = max( min( edge_val(k,1), max(u(km1), u(k)) ), min(u(km1), u(k)) )
    edge_val(k,2) = max( min( edge_val(k,2), max(u(kp1), u(k)) ), min(u(kp1), u(k)) )
  enddo

end subroutine bound_edge_values

!> Replace discontinuous edge values with their average when not monotonic.
!!
!! Copied from regrid_edge_values.check_discontinuous_edge_values().
subroutine check_discontinuous_edge_values(N, u, edge_val)
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: u !< Cell averages [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Edge values [A]

  integer :: k
  real    :: u0_avg

  do k = 1,N-1
    if ( (edge_val(k+1,1) - edge_val(k,2)) * (u(k+1) - u(k)) < 0.0 ) then
      u0_avg = 0.5 * ( edge_val(k,2) + edge_val(k+1,1) )
      u0_avg = max( min( u0_avg, max(u(k), u(k+1)) ), min(u(k), u(k+1)) )
      edge_val(k,2) = u0_avg
      edge_val(k+1,1) = u0_avg
    endif
  enddo

end subroutine check_discontinuous_edge_values

end module Recon1d_PPM_hybgen
