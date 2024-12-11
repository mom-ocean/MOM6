!> Piecewise Linear Method 1D reconstruction in index space and boundary extrapolation
!!
!! This implementation of PLM follows Colella and Woodward, 1984 \cite colella1984, except for assuming
!! uniform resolution so that the method is independent of grid spacing. The cell-wise reconstructions
!! are limited so that the edge values (which are also the extrema in a cell) are bounded by the neighbors.
!! The slope of the first and last cells are set so that the first interior edge values match the interior
!! cell (i.e. extrapolates from the interior).
module Recon1d_EMPLM_CWK

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_type, only : testing
use Recon1d_MPLM_CWK, only : MPLM_CWK

implicit none ; private

public EMPLM_CWK, testing

!> PLM reconstruction following Colella and Woodward, 1984
!!
!! Implemented by extending recon1d_mplm_cwk.
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_mplm_cwk -> recon1d_plm_cw.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_mplm_cwk -> recon1d_plm_cw.average()
!! - f()                    -> recon1d_mplm_cwk -> recon1d_plm_cw.f()
!! - dfdx()                 -> recon1d_mplm_cwk -> recon1d_plm_cw.dfdx()
!! - check_reconstruction() -> recon1d_mplm_cwk.check_reconstruction()
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_mplm_cwk -> recon1d_plm_cw.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> recon1d_mplm_cwk.reconstruct()
type, extends (MPLM_CWK) :: EMPLM_CWK

contains
  !> Implementation of the EMPLM_CWK reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of unit tests for the EMPLM_CWK reconstruction
  procedure :: unit_tests => unit_tests

end type EMPLM_CWK

contains

!> Calculate a 1D PLM reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(EMPLM_CWK), intent(inout) :: this !< This reconstruction
  real,             intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,             intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  real :: slp ! The PLM slopes (difference across cell) [A]
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
                                    ! differences across the cell [A]
  real :: u_min, u_max ! Minimum and maximum value across cell [A]
  real :: u_l, u_r, u_c ! Left, right, and center values [A]
  real :: u_e(this%n+1) ! Average of edge values [A]
  integer :: k, n

  n = this%n

  call this%reconstruct_parent(h, u)

  this%ur(1) = this%ul(2)
  this%ul(1) = u(1) + ( u(1) - this%ur(1) )

  this%ul(n) = this%ur(n-1)
  this%ur(n) = u(n) + ( u(n) - this%ul(n) )

end subroutine reconstruct

!> Runs PLM reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(EMPLM_CWK), intent(inout) :: this    !< This reconstruction
  logical,          intent(in)    :: verbose !< True, if verbose
  integer,          intent(in)    :: stdout  !< I/O channel for stdout
  integer,          intent(in)    :: stderr  !< I/O channel for stderr
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
    um(k) = this%average(k, 0.25, 0.75) ! Average from x=0.25 to 0.75 in each cell
  enddo
  call test%real_arr(3, um, (/1.,3.,5./), 'Return interval average')

  call this%destroy()
  deallocate( um, ul, ur, ull, urr )

  allocate( um(4), ul(4), ur(4) )
  call this%init(4)

  ! These values lead to non-monotonic reconstuctions which are
  ! valid for transport problems but not always appropriate for
  ! remapping to arbitrary resolution grids.
  ! The O(h^2) slopes are -, 2, 2, - and the limited
  ! slopes are 0, 1, 1, 0 so the everywhere the reconstructions
  ! are bounded by neighbors but ur(2) and ul(3) are out-of-order.
  call this%reconstruct( (/1.,1.,1.,1./), (/0.,3.,4.,7./) )
  do k = 1, 4
    ul(k) = this%f(k, 0.)
    ur(k) = this%f(k, 1.)
  enddo
  call test%real_arr(4, ul, (/-2.5,2.5,3.5,4.5/), 'Evaluation on left edge')
  call test%real_arr(4, ur, (/2.5,3.5,4.5,9.5/), 'Evaluation on right edge')

  deallocate( um, ul, ur )

  unit_tests = test%summarize('EMPLM_CWK:unit_tests')

end function unit_tests

!> \namespace recon1d_emplm_cwk
!!

end module Recon1d_EMPLM_CWK
