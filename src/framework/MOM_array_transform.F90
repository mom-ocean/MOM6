!> Module for supporting the rotation of a field's index map.
!! The implementation of each angle is described below.
!!
!! +90deg: B(i,j) = A(n-j,i)
!!                = transpose, then row reverse
!! 180deg: B(i,j) = A(m-i,n-j)
!!                = row reversal + column reversal
!! -90deg: B(i,j) = A(j,m-i)
!!                = row reverse, then transpose
!!
!! 90 degree rotations change the shape of the field, and are handled
!! separately from 180 degree rotations.

module MOM_array_transform

implicit none

private
public rotate_array
public rotate_array_pair
public rotate_vector
public allocate_rotated_array


!> Rotate the elements of an array to the rotated set of indices.
!! Rotation is applied across the first and second axes of the array.
interface rotate_array
  module procedure rotate_array_real_2d
  module procedure rotate_array_real_3d
  module procedure rotate_array_real_4d
  module procedure rotate_array_integer
  module procedure rotate_array_logical
end interface rotate_array


!> Rotate a pair of arrays which map to a rotated set of indices.
!! Rotation is applied across the first and second axes of the array.
!! This rotation should be applied when one field is mapped onto the other.
!! For example, a tracer indexed along u or v face points will map from one
!! to the other after a quarter turn, and back onto itself after a half turn.
interface rotate_array_pair
  module procedure rotate_array_pair_real_2d
  module procedure rotate_array_pair_real_3d
  module procedure rotate_array_pair_integer
end interface rotate_array_pair


!> Rotate an array pair representing the components of a vector.
!! Rotation is applied across the first and second axes of the array.
!! This rotation should be applied when the fields satisfy vector
!! transformation rules.  For example, the u and v components of a velocity
!! will map from one to the other for quarter turns, with a sign change in one
!! component.  A half turn will map elements onto themselves with sign changes
!! in both components.
interface rotate_vector
  module procedure rotate_vector_real_2d
  module procedure rotate_vector_real_3d
  module procedure rotate_vector_real_4d
end interface rotate_vector


!> Allocate an array based on the rotated index map of an unrotated reference
!! array.
interface allocate_rotated_array
  module procedure allocate_rotated_array_real_2d
  module procedure allocate_rotated_array_real_3d
  module procedure allocate_rotated_array_real_4d
  module procedure allocate_rotated_array_integer
end interface allocate_rotated_array

contains

!> Rotate the elements of a 2d real array along first and second axes.
subroutine rotate_array_real_2d(A_in, turns, A)
  real, intent(in) :: A_in(:,:) !< Unrotated array
  integer, intent(in) :: turns  !< Number of quarter turns
  real, intent(out) :: A(:,:)   !< Rotated array

  integer :: m, n

  m = size(A_in, 1)
  n = size(A_in, 2)

  select case (modulo(turns, 4))
    case(0)
      A(:,:) = A_in(:,:)
    case(1)
      A(:,:) = transpose(A_in)
      A(:,:) = A(n:1:-1, :)
    case(2)
      A(:,:) = A_in(m:1:-1, n:1:-1)
    case(3)
      A(:,:) = transpose(A_in(m:1:-1, :))
  end select
end subroutine rotate_array_real_2d


!> Rotate the elements of a 3d real array along first and second axes.
subroutine rotate_array_real_3d(A_in, turns, A)
  real, intent(in) :: A_in(:,:,:) !< Unrotated array
  integer, intent(in) :: turns    !< Number of quarter turns
  real, intent(out) :: A(:,:,:)   !< Rotated array

  integer :: k

  do k = 1, size(A_in, 3)
    call rotate_array(A_in(:,:,k), turns, A(:,:,k))
  enddo
end subroutine rotate_array_real_3d


!> Rotate the elements of a 4d real array along first and second axes.
subroutine rotate_array_real_4d(A_in, turns, A)
  real, intent(in) :: A_in(:,:,:,:) !< Unrotated array
  integer, intent(in) :: turns      !< Number of quarter turns
  real, intent(out) :: A(:,:,:,:)   !< Rotated array

  integer :: n

  do n = 1, size(A_in, 4)
    call rotate_array(A_in(:,:,:,n), turns, A(:,:,:,n))
  enddo
end subroutine rotate_array_real_4d


!> Rotate the elements of a 2d integer array along first and second axes.
subroutine rotate_array_integer(A_in, turns, A)
  integer, intent(in) :: A_in(:,:)  !< Unrotated array
  integer, intent(in) :: turns      !< Number of quarter turns
  integer, intent(out) :: A(:,:)    !< Rotated array

  integer :: m, n

  m = size(A_in, 1)
  n = size(A_in, 2)

  select case (modulo(turns, 4))
    case(0)
      A(:,:) = A_in(:,:)
    case(1)
      A(:,:) = transpose(A_in)
      A(:,:) = A(n:1:-1, :)
    case(2)
      A(:,:) = A_in(m:1:-1, n:1:-1)
    case(3)
      A(:,:) = transpose(A_in(m:1:-1, :))
  end select
end subroutine rotate_array_integer


!> Rotate the elements of a 2d logical array along first and second axes.
subroutine rotate_array_logical(A_in, turns, A)
  logical, intent(in) :: A_in(:,:)  !< Unrotated array
  integer, intent(in) :: turns      !< Number of quarter turns
  logical, intent(out) :: A(:,:)    !< Rotated array

  integer :: m, n

  m = size(A_in, 1)
  n = size(A_in, 2)

  select case (modulo(turns, 4))
    case(0)
      A(:,:) = A_in(:,:)
    case(1)
      A(:,:) = transpose(A_in)
      A(:,:) = A(n:1:-1, :)
    case(2)
      A(:,:) = A_in(m:1:-1, n:1:-1)
    case(3)
      A(:,:) = transpose(A_in(m:1:-1, :))
  end select
end subroutine rotate_array_logical


!> Rotate the elements of a 2d real array pair along first and second axes.
subroutine rotate_array_pair_real_2d(A_in, B_in, turns, A, B)
  real, intent(in) :: A_in(:,:)   !< Unrotated scalar array pair
  real, intent(in) :: B_in(:,:)   !< Unrotated scalar array pair
  integer, intent(in) :: turns    !< Number of quarter turns
  real, intent(out) :: A(:,:)     !< Rotated scalar array pair
  real, intent(out) :: B(:,:)     !< Rotated scalar array pair

  if (modulo(turns, 2) /= 0) then
    call rotate_array(B_in, turns, A)
    call rotate_array(A_in, turns, B)
  else
    call rotate_array(A_in, turns, A)
    call rotate_array(B_in, turns, B)
  endif
end subroutine rotate_array_pair_real_2d


!> Rotate the elements of a 3d real array pair along first and second axes.
subroutine rotate_array_pair_real_3d(A_in, B_in, turns, A, B)
  real, intent(in) :: A_in(:,:,:)   !< Unrotated scalar array pair
  real, intent(in) :: B_in(:,:,:)   !< Unrotated scalar array pair
  integer, intent(in) :: turns      !< Number of quarter turns
  real, intent(out) :: A(:,:,:)     !< Rotated scalar array pair
  real, intent(out) :: B(:,:,:)     !< Rotated scalar array pair

  integer :: k

  do k = 1, size(A_in, 3)
    call rotate_array_pair(A_in(:,:,k), B_in(:,:,k), turns, &
        A(:,:,k), B(:,:,k))
  enddo
end subroutine rotate_array_pair_real_3d


!> Rotate the elements of a 4d real array pair along first and second axes.
subroutine rotate_array_pair_integer(A_in, B_in, turns, A, B)
  integer, intent(in) :: A_in(:,:)  !< Unrotated scalar array pair
  integer, intent(in) :: B_in(:,:)  !< Unrotated scalar array pair
  integer, intent(in) :: turns      !< Number of quarter turns
  integer, intent(out) :: A(:,:)    !< Rotated scalar array pair
  integer, intent(out) :: B(:,:)    !< Rotated scalar array pair

  if (modulo(turns, 2) /= 0) then
    call rotate_array(B_in, turns, A)
    call rotate_array(A_in, turns, B)
  else
    call rotate_array(A_in, turns, A)
    call rotate_array(B_in, turns, B)
  endif
end subroutine rotate_array_pair_integer


!> Rotate the elements of a 2d real vector along first and second axes.
subroutine rotate_vector_real_2d(A_in, B_in, turns, A, B)
  real, intent(in) :: A_in(:,:) !< First component of unrotated vector
  real, intent(in) :: B_in(:,:) !< Second component of unrotated vector
  integer, intent(in) :: turns  !< Number of quarter turns
  real, intent(out) :: A(:,:)   !< First component of rotated vector
  real, intent(out) :: B(:,:)   !< Second component of unrotated vector

  call rotate_array_pair(A_in, B_in, turns, A, B)

  if (modulo(turns, 4) == 1 .or. modulo(turns, 4) == 2) &
    A(:,:) = -A(:,:)

  if (modulo(turns, 4) == 2 .or. modulo(turns, 4) == 3) &
    B(:,:) = -B(:,:)
end subroutine rotate_vector_real_2d


!> Rotate the elements of a 3d real vector along first and second axes.
subroutine rotate_vector_real_3d(A_in, B_in, turns, A, B)
  real, intent(in) :: A_in(:,:,:) !< First component of unrotated vector
  real, intent(in) :: B_in(:,:,:) !< Second component of unrotated vector
  integer, intent(in) :: turns    !< Number of quarter turns
  real, intent(out) :: A(:,:,:)   !< First component of rotated vector
  real, intent(out) :: B(:,:,:)   !< Second component of unrotated vector

  integer :: k

  do k = 1, size(A_in, 3)
    call rotate_vector(A_in(:,:,k), B_in(:,:,k), turns, A(:,:,k), B(:,:,k))
  enddo
end subroutine rotate_vector_real_3d


!> Rotate the elements of a 4d real vector along first and second axes.
subroutine rotate_vector_real_4d(A_in, B_in, turns, A, B)
  real, intent(in) :: A_in(:,:,:,:) !< First component of unrotated vector
  real, intent(in) :: B_in(:,:,:,:) !< Second component of unrotated vector
  integer, intent(in) :: turns      !< Number of quarter turns
  real, intent(out) :: A(:,:,:,:)   !< First component of rotated vector
  real, intent(out) :: B(:,:,:,:)   !< Second component of unrotated vector

  integer :: n

  do n = 1, size(A_in, 4)
    call rotate_vector(A_in(:,:,:,n), B_in(:,:,:,n), turns, &
        A(:,:,:,n), B(:,:,:,n))
  enddo
end subroutine rotate_vector_real_4d


!> Allocate a 2d real array on the rotated index map of a reference array.
subroutine allocate_rotated_array_real_2d(A_in, lb, turns, A)
  ! NOTE: lb must be declared before A_in
  integer, intent(in) :: lb(2)                !< Lower index bounds of A_in
  real, intent(in) :: A_in(lb(1):, lb(2):)    !< Reference array
  integer, intent(in) :: turns                !< Number of quarter turns
  real, allocatable, intent(inout) :: A(:,:)  !< Array on rotated index

  integer :: ub(2)

  ub(:) = ubound(A_in)

  if (modulo(turns, 2) /= 0) then
    allocate(A(lb(2):ub(2), lb(1):ub(1)))
  else
    allocate(A(lb(1):ub(1), lb(2):ub(2)))
  endif
end subroutine allocate_rotated_array_real_2d


!> Allocate a 3d real array on the rotated index map of a reference array.
subroutine allocate_rotated_array_real_3d(A_in, lb, turns, A)
  ! NOTE: lb must be declared before A_in
  integer, intent(in) :: lb(3)                    !< Lower index bounds of A_in
  real, intent(in) :: A_in(lb(1):, lb(2):, lb(3):)  !< Reference array
  integer, intent(in) :: turns                    !< Number of quarter turns
  real, allocatable, intent(inout) :: A(:,:,:)    !< Array on rotated index

  integer :: ub(3)

  ub(:) = ubound(A_in)

  if (modulo(turns, 2) /= 0) then
    allocate(A(lb(2):ub(2), lb(1):ub(1), lb(3):ub(3)))
  else
    allocate(A(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))
  endif
end subroutine allocate_rotated_array_real_3d


!> Allocate a 4d real array on the rotated index map of a reference array.
subroutine allocate_rotated_array_real_4d(A_in, lb, turns, A)
  ! NOTE: lb must be declared before A_in
  integer, intent(in) :: lb(4)                    !< Lower index bounds of A_in
  real, intent(in) :: A_in(lb(1):,lb(2):,lb(3):,lb(4):) !< Reference array
  integer, intent(in) :: turns                    !< Number of quarter turns
  real, allocatable, intent(inout) :: A(:,:,:,:)  !< Array on rotated index

  integer:: ub(4)

  ub(:) = ubound(A_in)

  if (modulo(turns, 2) /= 0) then
    allocate(A(lb(2):ub(2), lb(1):ub(1), lb(3):ub(3), lb(4):ub(4)))
  else
    allocate(A(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))
  endif
end subroutine allocate_rotated_array_real_4d


!> Allocate a 2d integer array on the rotated index map of a reference array.
subroutine allocate_rotated_array_integer(A_in, lb, turns, A)
  integer, intent(in) :: lb(2)                  !< Lower index bounds of A_in
  integer, intent(in) :: A_in(lb(1):,lb(2):)    !< Reference array
  integer, intent(in) :: turns                  !< Number of quarter turns
  integer, allocatable, intent(inout) :: A(:,:) !< Array on rotated index

  integer :: ub(2)

  ub(:) = ubound(A_in)

  if (modulo(turns, 2) /= 0) then
    allocate(A(lb(2):ub(2), lb(1):ub(1)))
  else
    allocate(A(lb(1):ub(1), lb(2):ub(2)))
  endif
end subroutine allocate_rotated_array_integer

end module MOM_array_transform
