!> Piecewise Linear Method 1D reconstruction in index space
!!
!! This implementation of PLM follows Colella and Woodward, 1984 \cite colella1984, except for assuming
!! uniform resolution so that the method is independent of grid spacing. The cell-wise reconstructions
!! are limited so that the edge values (which are also the extrema in a cell) are bounded by the neighbors.
!! The first and last cells are always limited to PCM.
module Recon1d_MPLM_CWK

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_type, only : testing
use Recon1d_PLM_CWK, only : PLM_CWK

implicit none ; private

public MPLM_CWK, testing

!> PLM reconstruction following Colella and Woodward, 1984
!!
!! Implemented by extending recon1d_plm_cwk.
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_plm_cwk -> recon1d_plm_cw.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_plm_cwk -> recon1d_plm_cw.average()
!! - f()                    -> recon1d_plm_cwk -> recon1d_plm_cw.f()
!! - dfdx()                 -> recon1d_plm_cwk -> recon1d_plm_cw.dfdx()
!! - check_reconstruction()    *locally defined
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_plm_cwk -> recon1d_plm_cw.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (PLM_CWK) :: MPLM_CWK

contains
  !> Implementation of the MPLM_CWK reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of check reconstruction for the MPLM_CWK reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the MPLM_CWK reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct
end type MPLM_CWK

contains

!> Calculate a 1D PLM reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(MPLM_CWK), intent(inout) :: this !< This reconstruction
  real,            intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,            intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  real :: slp ! The PLM slopes (difference across cell) [A]
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
                                    ! differences across the cell [A]
  real :: u_min, u_max ! Minimum and maximum value across cell [A]
  real :: u_l, u_r, u_c ! Left, right, and center values [A]
  real :: u_e(this%n+1) ! Average of edge values [A]
  integer :: k, n

  n = this%n

  ! Loop over all cells
  do k = 1, n
    this%u_mean(k) = u(k)
  enddo

  ! Boundary cells use PCM
  this%ul(1) = u(1)
  this%ur(1) = u(1)

  ! Loop over interior cells
  do k = 2, n-1
    u_l = u(k-1)
    u_c = u(k)
    u_r = u(k+1)

    ! Side differences
    sigma_r = u_r - u_c
    sigma_l = u_c - u_l

    ! This is the second order slope given by equation 1.7 of
    ! Piecewise Parabolic Method, Colella and Woodward (1984),
    ! http://dx.doi.org/10.1016/0021-991(84)90143-8.
    ! For uniform resolution it simplifies to ( u_r - u_l )/2 .
    sigma_c = 0.5 * ( u_r - u_l )

    ! Limit slope so that reconstructions are bounded by neighbors
    u_min = min( u_l, u_c, u_r )
    u_max = max( u_l, u_c, u_r )

    if ( (sigma_l * sigma_r) > 0.0 ) then
      ! This limits the slope so that the edge values are bounded by the two cell averages spanning the edge
      slp = sign( min( abs(sigma_c), 2.*min( u_c - u_min, u_max - u_c ) ), sigma_c )
    else
      ! Extrema in the mean values require a PCM reconstruction
      slp = 0.0
    endif

    ! Left edge
    u_min = min( u_c, u_l )
    u_max = max( u_c, u_l )
    u_l = u_c - 0.5 * slp
    this%ul(k) = max( min( u_l, u_max), u_min )

    ! Right edge
    u_min = min( u_c, u_r )
    u_max = max( u_c, u_r )
    u_r = u_c + 0.5 * slp
    this%ur(k) = max( min( u_r, u_max), u_min )
  enddo

  ! Boundary cells use PCM
  this%ul(n) = u(n)
  this%ur(n) = u(n)

  ! Average edge values
  u_e(1) = this%ul(1)
  do K = 2, n
    u_e(K) = 0.5 * ( this%ur(k-1) + this%ul(k) )
  enddo
  u_e(n+1) = this%ur(n)

  ! Loop over interior cells, redo PLM slope limiting using average edge as neighbor cell values
  do k = 2, n-1
    u_l = u_e(k)
    u_c = u(k)
    u_r = u_e(k+1)

    ! Side differences
    sigma_r = u_r - u_c
    sigma_l = u_c - u_l

    ! This is the second order slope given by equation 1.7 of
    ! Piecewise Parabolic Method, Colella and Woodward (1984),
    ! http://dx.doi.org/10.1016/0021-991(84)90143-8.
    ! For uniform resolution it simplifies to ( u_r - u_l )/2 .
    sigma_c = this%ur(k) - this%ul(k)

    ! Limit slope so that reconstructions are bounded by neighbors
    u_min = min( u_l, u_c, u_r )
    u_max = max( u_l, u_c, u_r )

    if ( (sigma_l * sigma_r) > 0.0 ) then
      ! This limits the slope so that the edge values are bounded by the two cell averages spanning the edge
      slp = sign( min( abs(sigma_c), 2.*min( u_c - u_min, u_max - u_c ) ), sigma_c )
    else
      ! Extrema in the mean values require a PCM reconstruction
      slp = 0.0
    endif

    ! Left edge
    u_min = min( u_c, u_l )
    u_max = max( u_c, u_l )
    u_l = u_c - 0.5 * slp
    this%ul(k) = max( min( u_l, u_max), u_min )

    ! Right edge
    u_min = min( u_c, u_r )
    u_max = max( u_c, u_r )
    u_r = u_c + 0.5 * slp
    this%ur(k) = max( min( u_r, u_max), u_min )
  enddo

end subroutine reconstruct

!> Checks the MPLM_CWK reconstruction for consistency
logical function check_reconstruction(this, h, u)
  class(MPLM_CWK), intent(in) :: this !< This reconstruction
  real,            intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,            intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: k

  check_reconstruction = .false.

  do k = 1, this%n
    if ( abs( this%u_mean(k) - u(k) ) > 0. ) check_reconstruction = .true.
  enddo

  ! Check the cell reconstruction is monotonic within each cell (it should be as a straight line)
  do k = 1, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ur(k) - this%u_mean(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check the cell is a straight line (to within machine precision)
  do k = 1, this%n
    if ( abs(2. * this%u_mean(k) - ( this%ul(k) + this%ur(k) )) > epsilon(this%u_mean(1)) * &
         max(abs(2. * this%u_mean(k)), abs(this%ul(k)), abs(this%ur(k))) ) check_reconstruction = .true.
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

!> Runs PLM reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(MPLM_CWK), intent(inout) :: this    !< This reconstruction
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
  call test%real_arr(3, ul, (/1.,2.,5./), 'Evaluation on left edge')
  call test%real_arr(3, um, (/1.,3.,5./), 'Evaluation in center')
  call test%real_arr(3, ur, (/1.,4.,5./), 'Evaluation on right edge')

  do k = 1, 3
    ul(k) = this%dfdx(k, 0.)
    um(k) = this%dfdx(k, 0.5)
    ur(k) = this%dfdx(k, 1.)
  enddo
  call test%real_arr(3, ul, (/0.,2.,0./), 'dfdx on left edge')
  call test%real_arr(3, um, (/0.,2.,0./), 'dfdx in center')
  call test%real_arr(3, ur, (/0.,2.,0./), 'dfdx on right edge')

  do k = 1, 3
    um(k) = this%average(k, 0.5, 0.75) ! Average from x=0.25 to 0.75 in each cell
  enddo
  call test%real_arr(3, um, (/1.,3.25,5./), 'Return interval average')

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
  call test%real_arr(4, ul, (/0.,2.5,3.5,7./), 'Evaluation on left edge')
  call test%real_arr(4, ur, (/0.,3.5,4.5,7./), 'Evaluation on right edge')

  deallocate( um, ul, ur )

  unit_tests = test%summarize('MPLM_CWK:unit_tests')

end function unit_tests

!> \namespace recon1d_mplm_cwk
!!

end module Recon1d_MPLM_CWK
