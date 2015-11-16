module PCM_functions
!==============================================================================
!
! This file is part of MOM.
!
! Date of creation: 2008.06.06
! L. White
!
! This module contains routines that handle one-dimensionnal finite volume
! reconstruction using the piecewise constant method (PCM).
!
!==============================================================================

implicit none ; private

public PCM_reconstruction

contains

!------------------------------------------------------------------------------
! pcm_reconstruction
!------------------------------------------------------------------------------
subroutine PCM_reconstruction( N, u, ppoly_E, ppoly_coefficients )
!------------------------------------------------------------------------------
! Reconstruction by constant polynomials within each cell. There is nothing to
! do but this routine is provided to ensure a homogeneous interface
! throughout the regridding toolbox.
!
! N:     number of cells in grid
! h:     thicknesses of grid cells
! u:     cell averages to use in constructing piecewise polynomials
! ppoly_E : edge values of piecewise polynomials
! ppoly_coefficients : coefficients of piecewise polynomials
!
! It is assumed that the dimension of 'u' is equal to the number of cells
! defining 'grid' and 'ppoly'. No consistency check is performed.
!------------------------------------------------------------------------------

  ! Arguments
  integer,              intent(in)    :: N ! Number of cells
  real, dimension(:),   intent(in)    :: u ! cell averages
  real, dimension(:,:), intent(inout) :: ppoly_E            !Edge value of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_coefficients !Coefficients of polynomial

  ! Local variables
  integer :: k

  ! The coefficients of the piecewise constant polynomial are simply
  ! the cell averages.
  ppoly_coefficients(:,1) = u(:)

  ! The edge values are equal to the cell average
  do k = 1,N
    ppoly_E(k,:) = u(k)
  end do

end subroutine PCM_reconstruction

end module PCM_functions
