!> Solvers of linear systems.
module regrid_solvers

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL

implicit none ; private

public :: solve_linear_system, linear_solver, solve_tridiagonal_system, solve_diag_dominant_tridiag

contains

!> Solve the linear system AX = R by Gaussian elimination
!!
!! This routine uses Gauss's algorithm to transform the system's original
!! matrix into an upper triangular matrix. Back substitution yields the answer.
!! The matrix A must be square, with the first index varing down the column.
subroutine solve_linear_system( A, R, X, N, answers_2018 )
  integer,              intent(in)    :: N  !< The size of the system
  real, dimension(N,N), intent(inout) :: A  !< The matrix being inverted [nondim]
  real, dimension(N),   intent(inout) :: R  !< system right-hand side [A]
  real, dimension(N),   intent(inout) :: X  !< solution vector [A]
  logical,    optional, intent(in)    :: answers_2018 !< If true or absent use older, less efficient expressions.
  ! Local variables
  real, parameter       :: eps = 0.0        ! Minimum pivot magnitude allowed
  real    :: factor       ! The factor that eliminates the leading nonzero element in a row.
  real    :: pivot, I_pivot ! The pivot value and its reciprocal [nondim]
  real    :: swap_a, swap_b
  logical :: found_pivot  ! If true, a pivot has been found
  logical :: old_answers  ! If true, use expressions that give the original (2008 through 2018) MOM6 answers
  integer :: i, j, k

  old_answers = .true. ; if (present(answers_2018)) old_answers = answers_2018

  ! Loop on rows to transform the problem into multiplication by an upper-right matrix.
  do i = 1,N-1


    ! Start to look for a pivot in the current row, i.  If the pivot in row i is not valid,
    ! keep looking for a valid pivot by searching the entries of column i in rows below row i.
    ! Once a valid pivot is found (say in row k), rows i and k are swaped.
    found_pivot = .false.
    k = i
    do while ( ( .NOT. found_pivot ) .AND. ( k <= N ) )
      if ( abs( A(k,i) ) > eps ) then  ! A valid pivot has been found
        found_pivot = .true.
      else                             ! Seek a valid pivot in the next row
        k = k + 1
      endif
    enddo ! end loop to find pivot

    ! If no pivot could be found, the system is singular.
    if ( .NOT. found_pivot ) then
      write(0,*) ' A=',A
      call MOM_error( FATAL, 'The linear system is singular !' )
    endif

    ! If the pivot is in a row that is different than row i, that is if
    ! k is different than i, we need to swap those two rows
    if ( k /= i ) then
      do j = 1,N
        swap_a = A(i,j) ; A(i,j) = A(k,j) ; A(k,j) = swap_a
      enddo
      swap_b = R(i) ; R(i) = R(k) ; R(k) = swap_b
    endif

    ! Transform pivot to 1 by dividing the entire row (right-hand side included) by the pivot
    if (old_answers) then
      pivot = A(i,i)
      do j = i,N ; A(i,j) = A(i,j) / pivot ; enddo
      R(i) = R(i) / pivot
    else
      I_pivot = 1.0 / A(i,i)
      A(i,i) = 1.0
      do j = i+1,N ; A(i,j) = A(i,j) * I_pivot ; enddo
      R(i) = R(i) * I_pivot
    endif

    ! #INV: At this point, A(i,i) is a suitable pivot and it is equal to 1

    ! Put zeros in column for all rows below that contain the pivot (which is row i)
    do k = i+1,N    ! k is the row index
      factor = A(k,i)
      ! A(k,i) = 0.0  ! These elements are not used again, so this line can be skipped for speed.
      do j = i+1,N  ! j is the column index
        A(k,j) = A(k,j) - factor * A(i,j)
      enddo
      R(k) = R(k) - factor * R(i)
    enddo

  enddo ! end loop on i

  ! Solve system by back substituting in what is now an upper-right matrix.
  X(N) = R(N) / A(N,N)  ! The last row is now trivially solved.
  do i = N-1,1,-1 ! loop on rows, starting from second to last row
    X(i) = R(i)
    do j = i+1,N
      X(i) = X(i) - A(i,j) * X(j)
    enddo
    if (old_answers) X(i) = X(i) / A(i,i)
  enddo

end subroutine solve_linear_system

!> Solve the linear system AX = R by Gaussian elimination
!!
!! This routine uses Gauss's algorithm to transform the system's original
!! matrix into an upper triangular matrix. Back substitution then yields the answer.
!! The matrix A must be square, with the first index varing along the row.
subroutine linear_solver( N, A, R, X )
  integer,              intent(in)    :: N  !< The size of the system
  real, dimension(N,N), intent(inout) :: A  !< The matrix being inverted [nondim]
  real, dimension(N),   intent(inout) :: R  !< system right-hand side [A]
  real, dimension(N),   intent(inout) :: X  !< solution vector [A]

  ! Local variables
  real, parameter :: eps = 0.0   ! Minimum pivot magnitude allowed
  real    :: factor       ! The factor that eliminates the leading nonzero element in a row.
  real    :: I_pivot      ! The reciprocal of the pivot value [inverse of the input units of a row of A]
  real    :: swap
  logical :: found_pivot  ! If true, a pivot has been found
  integer :: i, j, k

  ! Loop on rows to transform the problem into multiplication by an upper-right matrix.
  do i=1,N-1
    ! Seek a pivot for column i starting in row i, and continuing into the remaining rows.  If the
    ! pivot is in a row other than i, swap them.  If no valid pivot is found, i = N+1 after this loop.
    do k=i,N ; if ( abs( A(i,k) ) > eps ) exit ; enddo ! end loop to find pivot
    if ( k > N ) then  ! No pivot could be found and the system is singular.
      write(0,*) ' A=',A
      call MOM_error( FATAL, 'The linear system is singular !' )
    endif

    ! If the pivot is in a row that is different than row i, swap those two rows, noting that both
    ! rows start with i-1 zero values.
    if ( k /= i ) then
      do j=i,N ; swap = A(j,i) ; A(j,i) = A(j,k) ; A(j,k) = swap ; enddo
      swap = R(i) ; R(i) = R(k) ; R(k) = swap
    endif

    ! Transform the pivot to 1 by dividing the entire row (right-hand side included) by the pivot
    I_pivot = 1.0 / A(i,i)
    A(i,i) = 1.0
    do j=i+1,N ; A(j,i) = A(j,i) * I_pivot ; enddo
    R(i) = R(i) * I_pivot

    ! Put zeros in column for all rows below that contain the pivot (which is row i)
    do k=i+1,N    ! k is the row index
      factor = A(i,k)
      ! A(i,k) = 0.0  ! These elements are not used again, so this line can be skipped for speed.
      do j=i+1,N ; A(j,k) = A(j,k) - factor * A(j,i) ; enddo
      R(k) = R(k) - factor * R(i)
    enddo

  enddo ! end loop on i

  ! Solve the system by back substituting into what is now an upper-right matrix.
  X(N) = R(N) / A(N,N)  ! The last row is now trivially solved.
  do i=N-1,1,-1 ! loop on rows, starting from second to last row
    X(i) = R(i)
    do j=i+1,N ; X(i) = X(i) - A(j,i) * X(j) ; enddo
  enddo

end subroutine linear_solver


!> Solve the tridiagonal system AX = R
!!
!! This routine uses Thomas's algorithm to solve the tridiagonal system AX = R.
!! (A is made up of lower, middle and upper diagonals)
subroutine solve_tridiagonal_system( Al, Ad, Au, R, X, N, answers_2018 )
  integer,            intent(in)  :: N   !< The size of the system
  real, dimension(N), intent(in)  :: Ad  !< Matrix center diagonal
  real, dimension(N), intent(in)  :: Al  !< Matrix lower diagonal
  real, dimension(N), intent(in)  :: Au  !< Matrix upper diagonal
  real, dimension(N), intent(in)  :: R   !< system right-hand side
  real, dimension(N), intent(out) :: X   !< solution vector
  logical,  optional, intent(in)  :: answers_2018 !< If true use older, less acccurate expressions.
  ! Local variables
  real, dimension(N) :: pivot, Al_piv
  real, dimension(N) :: c1       ! Au / pivot for the backward sweep
  real    :: I_pivot  ! The inverse of the most recent pivot
  integer :: k        ! Loop index
  logical :: old_answers  ! If true, use expressions that give the original (2008 through 2018) MOM6 answers

  old_answers = .true. ; if (present(answers_2018)) old_answers = answers_2018

  if (old_answers) then
    ! This version gives the same answers as the original (2008 through 2018) MOM6 code
    ! Factorization and forward sweep
    pivot(1) = Ad(1)
    X(1) = R(1)
    do k = 2,N
      Al_piv(k) = Al(k) / pivot(k-1)
      pivot(k) = Ad(k) - Al_piv(k) * Au(k-1)
      X(k) = R(k) - Al_piv(k) * X(k-1)
    enddo

    ! Backward sweep
    X(N) = R(N) / pivot(N)  ! This should be X(N) / pivot(N), but is OK if Al(N) = 0.
    do k = N-1,1,-1
      X(k) = ( X(k) - Au(k)*X(k+1) ) / pivot(k)
    enddo
  else
    ! This is a more typical implementation of a tridiagonal solver than the one above.
    ! It is mathematically equivalent but differs at roundoff, which can cascade up to larger values.

    ! Factorization and forward sweep
    I_pivot = 1.0 / Ad(1)
    X(1) = R(1) * I_pivot
    do k = 2,N
      c1(K-1) = Au(k-1) * I_pivot
      I_pivot = 1.0 / (Ad(k) - Al(k) * c1(K-1))
      X(k) = (R(k) - Al(k) * X(k-1)) * I_pivot
    enddo
    ! Backward sweep
    do k = N-1,1,-1
      X(k) = X(k) - c1(K) * X(k+1)
    enddo

  endif

end subroutine solve_tridiagonal_system


!> Solve the tridiagonal system AX = R
!!
!! This routine uses a variant of Thomas's algorithm to solve the tridiagonal system AX = R, in
!! a form that is guaranteed to avoid dividing by a zero pivot.  The matrix A is made up of
!! lower (Al) and upper diagonals (Au) and a central diagonal Ad = Ac+Al+Au, where
!! Al, Au, and Ac are all positive (or negative) definite.  However when Ac is smaller than
!! roundoff compared with (Al+Au), the answers are prone to inaccuracy.
subroutine solve_diag_dominant_tridiag( Al, Ac, Au, R, X, N )
  integer,            intent(in)  :: N   !< The size of the system
  real, dimension(N), intent(in)  :: Ac  !< Matrix center diagonal offset from Al + Au
  real, dimension(N), intent(in)  :: Al  !< Matrix lower diagonal
  real, dimension(N), intent(in)  :: Au  !< Matrix upper diagonal
  real, dimension(N), intent(in)  :: R   !< system right-hand side
  real, dimension(N), intent(out) :: X   !< solution vector
  ! Local variables
  real, dimension(N) :: c1       ! Au / pivot for the backward sweep
  real               :: d1       ! The next value of 1.0 - c1
  real               :: I_pivot  ! The inverse of the most recent pivot
  real               :: denom_t1 ! The first term in the denominator of the inverse of the pivot.
  integer            :: k        ! Loop index

  ! Factorization and forward sweep, in a form that will never give a division by a
  ! zero pivot for positive definite Ac, Al, and Au.
  I_pivot = 1.0 / (Ac(1) + Au(1))
  d1 = Ac(1) * I_pivot
  c1(1) = Au(1) * I_pivot
  X(1) = R(1) * I_pivot
  do k=2,N-1
    denom_t1 = Ac(k) + d1 * Al(k)
    I_pivot = 1.0 / (denom_t1 + Au(k))
    d1 = denom_t1 * I_pivot
    c1(k) = Au(k) * I_pivot
    X(k) = (R(k) - Al(k) * X(k-1)) * I_pivot
  enddo
  I_pivot = 1.0 / (Ac(N) + d1 * Al(N))
  X(N) = (R(N) - Al(N) * X(N-1)) * I_pivot
  ! Backward sweep
  do k=N-1,1,-1
    X(k) = X(k) - c1(k) * X(k+1)
  enddo

end subroutine solve_diag_dominant_tridiag


!> \namespace regrid_solvers
!!
!! Date of creation: 2008.06.12
!! L. White
!!
!! This module contains solvers of linear systems.
!! These routines have now been updated for greater efficiency, especially in special cases.

end module regrid_solvers
