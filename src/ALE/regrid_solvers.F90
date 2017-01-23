module regrid_solvers
!==============================================================================
!
! This file is part of MOM.
!
! Date of creation: 2008.06.12
! L. White
!
! This module contains solvers of linear systems.
! These routines could (should ?) be replaced later by more efficient ones.
!
!
!==============================================================================

use MOM_error_handler, only : MOM_error, FATAL

implicit none ; private

public :: solve_linear_system, solve_tridiagonal_system

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

! -----------------------------------------------------------------------------
! Solve the linear system AX = B
! -----------------------------------------------------------------------------
subroutine solve_linear_system( A, B, X, system_size )
! -----------------------------------------------------------------------------
! This routine uses Gauss's algorithm to transform the system's original
! matrix into an upper triangular matrix. Back substitution yields the answer.
! The matrix A must be square and its size must be that of the vectors B and X.
! -----------------------------------------------------------------------------

  ! Arguments
  real, dimension(:,:), intent(inout)   :: A
  real, dimension(:), intent(inout)     :: B
  real, dimension(:), intent(inout)     :: X
  integer                               :: system_size

  ! Local variables
  integer               :: i, j, k
  real, parameter       :: eps = 0.0        ! Minimum pivot magnitude allowed
  real                  :: factor
  real                  :: pivot
  real                  :: swap_a, swap_b
  logical               :: found_pivot      ! boolean indicating whether
                                            ! a pivot has been found
  ! Loop on rows
  do i = 1,system_size-1

    found_pivot = .false.

    ! Start to look for a pivot in row i. If the pivot
    ! in row i -- which is the current row -- is not valid,
    ! we keep looking for a valid pivot by searching the
    ! entries of column i in rows below row i. Once a valid
    ! pivot is found (say in row k), rows i and k are swaped.
    k = i
    do while ( ( .NOT. found_pivot ) .AND. ( k .LE. system_size ) )

        if ( abs( A(k,i) ) .GT. eps ) then  ! a valid pivot is found
          found_pivot = .true.
        else                                ! Go to the next row to see
                                            ! if there is a valid pivot there
          k = k + 1
        end if

    end do ! end loop to find pivot

    ! If no pivot could be found, the system is singular and we need
    ! to end the execution
    if ( .NOT. found_pivot ) then
      write(0,*) ' A=',A
      call MOM_error( FATAL, 'The linear system is singular !' )
    end if

    ! If the pivot is in a row that is different than row i, that is if
    ! k is different than i, we need to swap those two rows
    if ( k .NE. i ) then
      do j = 1,system_size
        swap_a = A(i,j)
        A(i,j) = A(k,j)
        A(k,j) = swap_a
      end do
      swap_b = B(i)
      B(i) = B(k)
      B(k) = swap_b
    end if

    ! Transform pivot to 1 by dividing the entire row
    ! (right-hand side included) by the pivot
    pivot = A(i,i)
    do j = i,system_size
      A(i,j) = A(i,j) / pivot
    end do
    B(i) = B(i) / pivot

    ! #INV: At this point, A(i,i) is a suitable pivot and it is equal to 1

    ! Put zeros in column for all rows below that containing
    ! pivot (which is row i)
    do k = (i+1),system_size    ! k is the row index
      factor = A(k,i)
      do j = (i+1),system_size      ! j is the column index
        A(k,j) = A(k,j) - factor * A(i,j)
      end do
      B(k) = B(k) - factor * B(i)
    end do

  end do ! end loop on i


  ! Solve system by back substituting
  X(system_size) = B(system_size) / A(system_size,system_size)
  do i = system_size-1,1,-1 ! loop on rows, starting from second to last row
    X(i) = B(i)
    do j = (i+1),system_size
      X(i) = X(i) - A(i,j) * X(j)
    end do
    X(i) = X(i) / A(i,i)
  end do

end subroutine solve_linear_system


! -----------------------------------------------------------------------------
! Solve the tridiagonal system AX = B
! -----------------------------------------------------------------------------
subroutine solve_tridiagonal_system( Al, Ad, Au, B, X, system_size )
! -----------------------------------------------------------------------------
! This routine uses Thomas's algorithm to solve the tridiagonal system AX = B.
! (A is made up of lower, middle and upper diagonals)
! -----------------------------------------------------------------------------
  ! Arguments
  real, dimension(:), intent(inout) :: Al, Ad, Au   ! lo., mid. and up. diagonals
  real, dimension(:), intent(inout) :: B            ! system right-hand side
  real, dimension(:), intent(inout) :: X            ! solution vector
  integer, intent(in)               :: system_size

  ! Local variables
  integer                               :: k        ! Loop index
  integer                               :: N        ! system size

  N = system_size

  ! Factorization
  do k = 1,N-1
    Al(k+1) = Al(k+1) / Ad(k)
    Ad(k+1) = Ad(k+1) - Al(k+1) * Au(k)
  end do

  ! Forward sweep
  do k = 2,N
    B(k) = B(k) - Al(k) * B(k-1)
  end do

  ! Backward sweep
  X(N) = B(N) / Ad(N)
  do k = N-1,1,-1
    X(k) = ( B(k) - Au(k)*X(k+1) ) / Ad(k)
  end do

end subroutine solve_tridiagonal_system

end module regrid_solvers
