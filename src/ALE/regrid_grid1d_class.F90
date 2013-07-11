module regrid_grid1d_class

implicit none ; private

public :: grid1D_t, grid1Dconstruct, grid1Ddestroy, grid1Dprint

! -----------------------------------------------------------------------------
! Definition of the one-dimensional grid structure
! -----------------------------------------------------------------------------
type grid1D_t
  ! Number of grid cells
  integer                           :: nb_cells

  ! Cell widths (size = nb_cells)
  real, dimension(:), allocatable   :: h

  ! Coordinates of nodes (size = nb_cells+1)
  real, dimension(:), allocatable   :: x

end type grid1D_t

contains

!------------------------------------------------------------------------------
! grid1Dconstruct
! -----------------------------------------------------------------------------
subroutine grid1Dconstruct( grid, nb_cells )
!------------------------------------------------------------------------------
! Initialization (memory allocation) of a grid
!------------------------------------------------------------------------------

  type(grid1D_t), intent(inout)  :: grid
  integer, intent(in)            :: nb_cells

  ! Set the number of cells
  grid%nb_cells = nb_cells

  ! Memory allocation
  allocate( grid%h(nb_cells) )
  allocate( grid%x(nb_cells+1) )

  ! Set all entries to zero
  grid%h(:) = 0.0
  grid%x(:) = 0.0

end subroutine grid1Dconstruct

!------------------------------------------------------------------------------
! grid1Dprint
! -----------------------------------------------------------------------------
subroutine grid1Dprint( ioChannel, grid )
!------------------------------------------------------------------------------
! Initialization (memory allocation) of a grid
!------------------------------------------------------------------------------

  integer, intent(in) :: ioChannel
  type(grid1D_t), intent(in)  :: grid
  integer :: k

  write(ioChannel,'(a4,a12,a12)') 'k','x','h'
  write(ioChannel,'(i4,es12.4)') 1, grid%x(1)
  do k = 1, grid%nb_cells
    write(ioChannel,'(i4,12x,es12.4)') k, grid%h(k)
    write(ioChannel,'(i4,es12.4)') k+1, grid%x(k+1)
  enddo

end subroutine grid1Dprint

!------------------------------------------------------------------------------
! grid1Ddestroy
! -----------------------------------------------------------------------------
subroutine grid1Ddestroy( grid )
!------------------------------------------------------------------------------
! Reclaim previously allocated memory for a grid
!------------------------------------------------------------------------------

  type(grid1D_t), intent(inout) :: grid

  deallocate( grid%h )
  deallocate( grid%x )

end subroutine grid1Ddestroy

end module regrid_grid1d_class
