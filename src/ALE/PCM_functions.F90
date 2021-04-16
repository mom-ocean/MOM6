!> Piecewise constant reconstruction functions
module PCM_functions

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public PCM_reconstruction

contains

!> Reconstruction by constant polynomials within each cell. There is nothing to
!! do but this routine is provided to ensure a homogeneous interface
!! throughout the regridding toolbox.
!!
!! It is assumed that the dimension of 'u' is equal to the number of cells
!! defining 'grid' and 'ppoly'. No consistency check is performed.
subroutine PCM_reconstruction( N, u, edge_values, ppoly_coef )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: u !< cell averages
  real, dimension(:,:), intent(inout) :: edge_values !< Edge value of polynomial,
                                           !! with the same units as u.
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Coefficients of polynomial,
                                           !! with the same units as u.

  ! Local variables
  integer :: k

  ! The coefficients of the piecewise constant polynomial are simply
  ! the cell averages.
  ppoly_coef(:,1) = u(:)

  ! The edge values are equal to the cell average
  do k = 1,N
    edge_values(k,:) = u(k)
  enddo

end subroutine PCM_reconstruction

!> \namespace PCM_functions
!!
!! Date of creation: 2008.06.06
!! L. White
!!
!! This module contains routines that handle one-dimensionnal finite volume
!! reconstruction using the piecewise constant method (PCM).

end module PCM_functions
