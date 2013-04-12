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
use regrid_grid1d_class, only : grid1D_t
use regrid_ppoly_class, only : ppoly_t

implicit none ; private

public PCM_reconstruction

contains

!------------------------------------------------------------------------------
! pcm_reconstruction
!------------------------------------------------------------------------------
subroutine PCM_reconstruction( grid, u, ppoly )
!------------------------------------------------------------------------------
! Reconstruction by constant polynomials within each cell. There is nothing to
! do but this routine is provided to ensure a homogeneous interface
! throughout the regridding toolbox.
!
! grid:  one-dimensional grid (properly initialized)
! ppoly: piecewise constant polynomial to be reconstructed (properly
!        initialized)
! u:     cell averages
!
! It is assumed that the dimension of 'u' is equal to the number of cells
! defining 'grid' and 'ppoly'. No consistency check is performed.
!------------------------------------------------------------------------------

  ! Arguments
  type(grid1D_t), intent(in)     :: grid
  real, dimension(:), intent(in) :: u
  type(ppoly_t), intent(inout)   :: ppoly

  ! Local variables
  integer   :: k

  ! The coefficients of the piecewise constant polynomial are simply
  ! the cell averages.
  ppoly%coefficients(:,1) = u

  ! The edge values are equal to the cell average
  do k = 1,grid%nb_cells
    ppoly%E(k,:) = u(k)
  end do

end subroutine PCM_reconstruction

end module PCM_functions
