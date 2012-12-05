module regrid_grid1d_class

implicit none ; private

public :: grid1d_t, grid1d_init, grid1d_destroy

! -----------------------------------------------------------------------------
! Definition of the one-dimensional grid structure
! -----------------------------------------------------------------------------
type grid1d_t
  ! Number of grid cells
  integer                           :: nb_cells;
  
  ! Cell widths (size = nb_cells)
  real, dimension(:), allocatable   :: h;
  
  ! Coordinates of nodes (size = nb_cells+1)
  real, dimension(:), allocatable   :: x;   

end type grid1d_t 

contains

!------------------------------------------------------------------------------
! grid_init
! -----------------------------------------------------------------------------
subroutine grid1d_init ( grid, nb_cells )
!------------------------------------------------------------------------------
! Initialization (memory allocation) of a grid
!------------------------------------------------------------------------------

  type(grid1d_t), intent(inout)  :: grid;
  integer, intent(in)            :: nb_cells;

  ! Set the number of cells
  grid%nb_cells = nb_cells;

  ! Memory allocation
  allocate ( grid%h(nb_cells) );
  allocate ( grid%x(nb_cells+1) );

  ! Set all entries to zero
  grid%h = 0.0
  grid%x = 0.0

end subroutine grid1d_init

!------------------------------------------------------------------------------
! grid_destroy
! -----------------------------------------------------------------------------
subroutine grid1d_destroy ( grid )
!------------------------------------------------------------------------------
! Reclaim previously allocated memory for a grid
!------------------------------------------------------------------------------

  type(grid1d_t), intent(inout) :: grid;

  deallocate ( grid%h );
  deallocate ( grid%x );

end subroutine grid1d_destroy

end module regrid_grid1d_class
