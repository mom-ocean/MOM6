!> Extrapolated-Monotonized Piecewise Linear Method 1D reconstruction
!!
!! This extends MPLM_WA, following White and Adcroft, 2008 \cite white2008, by extrapolating for the slopes of the
!! first and last cells. This extrapolation is used by White et al., 2009, during grid-generation.
module Recon1d_EMPLM_WA

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_MPLM_WA, only : MPLM_WA, testing

implicit none ; private

public EMPLM_WA

!> Extraplated Monotonic PLM reconstruction of White and Adcroft, 2008
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_mplm_wa           -> recon1d_plm_cw.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_mplm_wa           -> recon1d_plm_cw.average()
!! - f()                    -> recon1d_mplm_wa           -> recon1d_plm_cw.f()
!! - dfdx()                 -> recon1d_mplm_wa           -> recon1d_plm_cw.dfdx()
!! - check_reconstruction() -> recon1d_mplm_wa.check_reconstruction()
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_mplm_wa           -> recon1d_plm_cw.destroy()
!! - remap_to_sub_grid()    -> recon1d_mplm_wa           -> recon1d_plm_cw        -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> recon1d_mplm_wa           -> recon1d_plm_cw.init()
!! - reconstruct_parent()   -> recon1d_mplm_wa.reconstruct()
type, extends (MPLM_WA) :: EMPLM_WA

contains
  !> Implementation of the EMPLM_WA reconstruction with boundary extrapolation
  procedure :: reconstruct => reconstruct
  !> Implementation of unit tests for the EMPLM_WA reconstruction
  procedure :: unit_tests => unit_tests

end type EMPLM_WA

contains

!> Calculate a 1D PLM reconstruction based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(EMPLM_WA), intent(inout) :: this !< This reconstruction
  real,            intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,            intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: n
  real :: slope ! Difference of u across cell [A]
  real :: edge_h2 ! Edge value found by linear interpolation [A]
  real :: slope_h2 ! Twice the difference between cell center and 2nd order edge value [A]
  real :: slope_e ! Twice the difference between cell center and neighbor edge value [A]
  real :: hn, hc ! Neighbor and central cell thicknesses adjusted by h_neglect [H]
  real :: u_min, u_max ! Working values for bounding edge values [A]

  ! Use parent (MPLM_WA) reconstruction
  call this%reconstruct_parent(h, u)

  ! Fix reconstruction for first cell
  ! Avoid division by zero for vanished cells
  hn = h(2) + this%h_neglect
  hc = h(1) + this%h_neglect
  edge_h2 = ( u(2) * hc + u(1) * hn ) / ( hn + hc )
  slope_h2 = 2.0 * ( edge_h2 - u(1) )
  slope_e = 2.0 * ( this%ul(2) - u(1) )
  slope = sign( min( abs(slope_h2), abs(slope_e) ), u(2) - u(1) )
  edge_h2 = u(1) + 0.5 * slope
  u_min = min( this%ul(2), u(1) )
  u_max = max( this%ul(2), u(1) )
  this%ur(1) = max( u_min, min( u_max, edge_h2 ) )
  this%ul(1) = u(1) - 0.5 * slope
! slope = - PLM_extrapolate_slope( h(2), h(1), this%h_neglect, this%ul(2), u(1) )
! this%ul(1) = u(1) - 0.5 * slope
! this%ur(1) = u(1) + 0.5 * slope

  ! Fix reconstruction for last cell
  n = this%n
  ! Avoid division by zero for vanished cells
  hn = h(n-1) + this%h_neglect
  hc = h(n) + this%h_neglect
  edge_h2 = ( u(n-1) * hc + u(n) * hn ) / ( hn + hc )
  slope_h2 = 2.0 * ( u(n) - edge_h2 )
  slope_e = 2.0 * ( u(n) - this%ur(n-1) )
  slope = sign( min( abs(slope_h2), abs(slope_e) ), u(n) - u(n-1) )
  edge_h2 = u(n) - 0.5 * slope
  u_min = min( this%ur(n-1), u(n) )
  u_max = max( this%ur(n-1), u(n) )
  this%ul(n) = max( u_min, min( u_max, edge_h2 ) )
  this%ur(n) = u(n) + 0.5 * slope
! slope = PLM_extrapolate_slope( h(n-1), h(n), this%h_neglect, this%ur(n-1), u(n) )
! this%ul(n) = u(n) - 0.5 * slope
! this%ur(n) = u(n) + 0.5 * slope

end subroutine reconstruct

!> Returns a PLM slope using h2 extrapolation from a cell to the left, in the same
!! arbitrary units as the input values [A].
!! Use the negative to extrapolate from the cell to the right.
real elemental pure function PLM_extrapolate_slope(h_l, h_c, h_neglect, u_l, u_c)
  real, intent(in) :: h_l !< Thickness of left cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_c !< Thickness of center cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_neglect !< A negligible thickness [H]
  real, intent(in) :: u_l !< Value of left cell in arbitrary units [A]
  real, intent(in) :: u_c !< Value of center cell in arbitrary units [A]
  ! Local variables
  real :: left_edge ! Left edge value [A]
  real :: hl, hc ! Left and central cell thicknesses [H]

  ! Avoid division by zero for vanished cells
  hl = h_l + h_neglect
  hc = h_c + h_neglect

  ! The h2 scheme is used to compute the left edge value
  left_edge = (u_l*hc + u_c*hl) / (hl + hc)

  PLM_extrapolate_slope = 2.0 * ( u_c - left_edge )

end function PLM_extrapolate_slope

!> Runs PLM reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(EMPLM_WA), intent(inout) :: this    !< This reconstruction
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

  call this%init(3)
  call test%test( this%n /= 3, 'Setting number of levels')
  allocate( um(3), ul(3), ur(3), ull(3), urr(3) )

  call this%reconstruct( (/2.,2.,2./), (/1.,3.,5./) )
  call test%real_arr(3, this%u_mean, (/1.,3.,5./), 'Setting cell values')

  do k = 1, 3
    ul(k) = this%f(k, 0.)
    um(k) = this%f(k, 0.5)
    ur(k) = this%f(k, 1.)
  enddo
  call test%real_arr(3, ul, (/0.,2.,4./), 'Evaluation on left edge')
  call test%real_arr(3, um, (/1.,3.,5./), 'Evaluation in center')
  call test%real_arr(3, ur, (/2.,4.,6./), 'Evaluation on right edge')

  do k = 1, 3
    ul(k) = this%dfdx(k, 0.)
    um(k) = this%dfdx(k, 0.5)
    ur(k) = this%dfdx(k, 1.)
  enddo
  call test%real_arr(3, ul, (/2.,2.,2./), 'dfdx on left edge')
  call test%real_arr(3, um, (/2.,2.,2./), 'dfdx in center')
  call test%real_arr(3, ur, (/2.,2.,2./), 'dfdx on right edge')

  do k = 1, 3
    um(k) = this%average(k, 0.5, 0.75) ! Average from x=0.25 to 0.75 in each cell
  enddo
  call test%real_arr(3, um, (/1.25,3.25,5.25/), 'Return interval average')

  unit_tests = test%summarize('EMPLM_WA:unit_tests')

end function unit_tests

!> \namespace recon1d_emplm_wa
!!

end module Recon1d_EMPLM_WA
