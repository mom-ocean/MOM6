!> Piecewise Linear Method 1D reconstruction
!!
!! This implementation of PLM follows Colella and Woodward, 1984, except for assuming
!! uniform cell thicknesses. Cells resort to PCM for extrema including first and last cells in column.
!! The cell-wise reconstructions are limited so that the edge values (which are also the
!! extrema in a cell) are bounded by the neighbor cell means. However, this does not yield
!! monotonic profiles for the whole column.
!!
!! Note that internally the edge values, rather than the PLM slope, are stored to ensure
!! resulting calculations are properly bounded.
module Recon1d_PLM_CWK

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_type, only : testing
use Recon1d_PLM_CW, only : PLM_CW

implicit none ; private

public PLM_CWK, testing

!> PLM reconstruction following Colella and Woodward, 1984
!!
!! Implemented by extending recon1d_plm_cw.
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_plm_cw.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_plm_cw.average()
!! - f()                    -> recon1d_plm_cw.f()
!! - dfdx()                 -> recon1d_plm_cw.dfdx()
!! - check_reconstruction() -> recon1d_plm_cw.check_reconstruction()
!! - unit_tests()           -> recon1d_plm_cw.unit_tests()
!! - destroy()              -> recon1d_plm_cw.destroy()
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (PLM_CW) :: PLM_CWK

contains
  !> Implementation of the PLM_CWK reconstruction
  procedure :: reconstruct => reconstruct

end type PLM_CWK

contains

!> Calculate a 1D PLM reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(PLM_CWK), intent(inout) :: this !< This reconstruction
  real,           intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  real :: slp ! The PLM slopes (difference across cell) [A]
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
                                    ! differences across the cell [A]
  real :: u_min, u_max ! Minimum and maximum value across cell [A]
  real :: u_l, u_r, u_c ! Left, right, and center values [A]
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
    ! but for uniform resolution.
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

end subroutine reconstruct

!> \namespace recon1d_plm_cwk
!!

end module Recon1d_PLM_CWK
