module regrid_ppoly_class

implicit none ; private

public :: ppoly_t, ppoly_init, ppoly_destroy

! -----------------------------------------------------------------------------
! Definition of the piecewise polynomial structure
! -----------------------------------------------------------------------------
type ppoly_t
  ! Number of piecewise polynomials (i.e., number of grid cells)
  integer                           :: nb_cells
  
  ! 'E' and 'S' are arrays of edge values (E) and edge slopes (S). There
  ! are two edge values and slopes per cell
  real, dimension(:,:), allocatable :: E
  real, dimension(:,:), allocatable :: S

  ! This array holds the coefficients of each polynomial in each cell, 
  ! expressed in terms of the normalized coordinate xi \in [0,1]. The size 
  ! of this array is nb_cells x (degree + 1), where degree is the degree
  ! of the polynomial used for the reconstruction. E.g., a polynomial of
  ! degree 1 needs two coefficients. Note that polynomials are expressed
  ! as follows: P(\xi) = c_0 + c_1 \xi + c_2 \xi^2 + ...
  real, dimension(:,:), allocatable :: coefficients
  
end type ppoly_t 

contains

!------------------------------------------------------------------------------
! ppoly_init
! -----------------------------------------------------------------------------
subroutine ppoly_init( ppoly, nb_cells, degree )
!------------------------------------------------------------------------------
! Initialization (memory allocation) for a piecewise polynomial ppoly.
!------------------------------------------------------------------------------

  type(ppoly_t), intent(inout) :: ppoly
  integer, intent(in)          :: nb_cells
  integer, intent(in)          :: degree

  allocate( ppoly%E(nb_cells,2) )
  allocate( ppoly%S(nb_cells,2) )
  allocate( ppoly%coefficients(nb_cells,degree+1) )
  ppoly%nb_cells = nb_cells

end subroutine ppoly_init

!------------------------------------------------------------------------------
! ppoly_destroy
! -----------------------------------------------------------------------------
subroutine ppoly_destroy( ppoly )
!------------------------------------------------------------------------------
! Reclaim previously allocated memory for a piecewise polynomial ppoly.
!------------------------------------------------------------------------------

  type(ppoly_t), intent(inout) :: ppoly

  deallocate( ppoly%E )
  deallocate( ppoly%S )
  deallocate( ppoly%coefficients )

end subroutine ppoly_destroy

end module regrid_ppoly_class
