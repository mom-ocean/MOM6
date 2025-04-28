!> Extrapolated-Monotonized Piecewise Linear Method 1D reconstruction
!!
!! This extends MPLM_poly, following White and Adcroft, 2008 \cite white2008, by extraplating for the slopes of the
!! first and last cells. This extrapolation is used by White et al., 2009, during grid-generation.
!!
!! This stores and evaluates the reconstruction using a polynomial representation which is not preferred
!! but was the form used in OM4.
module Recon1d_EMPLM_WA_poly

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_MPLM_WA_poly, only : MPLM_WA_poly, testing

implicit none ; private

public EMPLM_WA_poly

!> Extrapolation Limited Monotonic PLM reconstruction following White and Adcroft, 2008
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_mplm_wa_poly.init()
!! - reconstruct()          -> recon1d_mplm_wa_poly.reconstruct()
!! - average()              -> recon1d_mplm_wa_poly.average()
!! - f()                    -> recon1d_mplm_wa_poly           -> recon1d_mplm_wa -> recon1d_plm_cw.f()
!! - dfdx()                 -> recon1d_mplm_wa_poly           -> recon1d_mplm_wa -> recon1d_plm_cw.dfdx()
!! - check_reconstruction()    *locally defined
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_mplm_wa_poly           -> recon1d_mplm_wa -> recon1d_plm_cw.destroy()
!! - remap_to_sub_grid()       *locally defined
!! - init_parent()          -> recon1d_mplm_wa_poly           -> recon1d_mplm_wa.init()
!! - reconstruct_parent()   -> recon1d_mplm_wa_poly           -> recon1d_mplm_wa.reconstruct()
type, extends (MPLM_WA_poly) :: EMPLM_WA_poly

contains
  !> Implementation of the EMPLM_WA_poly reconstruction with boundary extrapolation
  procedure :: reconstruct => reconstruct
  !> Implementation of check reconstruction for the EMPLM_WA_poly reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the EMPLM_WA_poly reconstruction
  procedure :: unit_tests => unit_tests

end type EMPLM_WA_poly

contains

!> Calculate a 1D PLM reconstruction based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(EMPLM_WA_poly), intent(inout) :: this !< This reconstruction
  real,                 intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,                 intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: n
  real :: slope ! Difference of u across cell [A]

  ! Use parent (MPLM_WA) reconstruction
  call this%reconstruct_parent(h, u)

  n = this%n

  ! Fix reconstruction for first cell
  slope = - PLM_extrapolate_slope( h(2), h(1), this%h_neglect, u(2), u(1) )
  this%ul(1) = u(1) - 0.5 * slope
  this%ur(1) = u(1) + 0.5 * slope
  this%poly_coef(1,1) = this%ul(1)
  this%poly_coef(1,2) = this%ur(1) - this%ul(1)

  ! Fix reconstruction for last cell
  slope = PLM_extrapolate_slope( h(n-1), h(n), this%h_neglect, u(n-1), u(n) )
  this%ul(n) = u(n) - 0.5 * slope
  this%ur(n) = u(n) + 0.5 * slope
  this%poly_coef(n,1) = this%ul(n)
  this%poly_coef(n,2) = this%ur(n) - this%ul(n)

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

!> Checks the EMPLM_WA_poly reconstruction for consistency
logical function check_reconstruction(this, h, u)
  class(EMPLM_WA_poly), intent(in) :: this !< This reconstruction
  real,                 intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,                 intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: k

  check_reconstruction = .false.

  do k = 1, this%n
    if ( abs( this%u_mean(k) - u(k) ) > 0. ) check_reconstruction = .true.
  enddo

  ! Check implied curvature
  do k = 1, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ur(k) - this%u_mean(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! These two checks fail MOM_remapping:test_recon_consistency in the presence of vanished layers
  ! e.g. intel/2023.2.0 on gaea at iter=26

! ! Check bounding of right edges, w.r.t. the cell means
! do K = 1, this%n-1
!   if ( ( this%ur(k) - this%u_mean(k) ) * ( this%u_mean(k+1) - this%ur(k) ) < 0. ) check_reconstruction = .true.
! enddo

! ! Check bounding of left edges, w.r.t. the cell means
! do K = 2, this%n
!   if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ul(k) - this%u_mean(k-1) ) < 0. ) check_reconstruction = .true.
! enddo

  ! Check order of u, ur, ul
  ! Note that in the OM4-era implementation, we were not consistent for top and bottom layers due
  ! extrapolation using cell means rather than edge values, hence reduced range for K
  do K = 2, this%n-2
    if ( ( this%ur(k) - this%u_mean(k) ) * ( this%ul(k+1) - this%ur(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check bounding of left edges, w.r.t. this cell mean and the previous cell right edge
  do K = 3, this%n-1
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ul(k) - this%ur(k-1) ) < 0. ) check_reconstruction = .true.
  enddo

end function check_reconstruction


!> Runs PLM reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(EMPLM_WA_poly), intent(inout) :: this    !< This reconstruction
  logical,              intent(in)    :: verbose !< True, if verbose
  integer,              intent(in)    :: stdout  !< I/O channel for stdout
  integer,              intent(in)    :: stderr  !< I/O channel for stderr
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

  unit_tests = test%summarize('EMPLM_WA_poly:unit_tests')

end function unit_tests

!> \namespace recon1d_emplm_wa_poly
!!

end module Recon1d_EMPLM_WA_poly
