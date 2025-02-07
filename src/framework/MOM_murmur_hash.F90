!> MurmurHash is a non-cryptographic hash function developed by Austin Appleby.
!!
!! This module provides an implementation of the 32-bit MurmurHash3 algorithm.
!! It is used in MOM6 to generate unique hashes of field arrays.  The hash is
!! sensitive to order of elements and can detect changes that would otherwise
!! be missed by the mean/min/max/bitcount tests.
!!
!! Sensitivity to order means that it must be used with care for tests such as
!! processor layout.
!!
!! This implementation assumes data sizes of either 32 or 64 bits.  It cannot
!! be used for smaller types such as strings.
!!
!! https://github.com/aappleby/smhasher
module MOM_murmur_hash

use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64

implicit none ; private

public :: murmur_hash

!> Return the murmur3 hash of an array.
interface murmur_hash
  procedure murmurhash3_i32
  procedure murmurhash3_i64
  procedure murmurhash3_r32
  procedure murmurhash3_r32_1d
  procedure murmurhash3_r32_2d
  procedure murmurhash3_r32_3d
  procedure murmurhash3_r32_4d
  procedure murmurhash3_r64
  procedure murmurhash3_r64_1d
  procedure murmurhash3_r64_2d
  procedure murmurhash3_r64_3d
  procedure murmurhash3_r64_4d
end interface murmur_hash

contains

!> Return the murmur3 hash for a 32-bit integer array.
function murmurhash3_i32(key, seed) result(hash)
  integer(int32), intent(in) :: key(:)
    !< Input array
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array

  integer(int32), parameter :: c1 = int(z'cc9e2d51', kind=int32)
  integer(int32), parameter :: c2 = int(z'1b873593', kind=int32)
  integer(int32), parameter :: c3 = int(z'e6546b64', kind=int32)

  integer(int32), parameter :: c4 = int(z'85ebca6b', kind=int32)
  integer(int32), parameter :: c5 = int(z'c2b2ae35', kind=int32)

  integer :: i
  integer(int32) :: k

  hash = 0
  if (present(seed)) hash = seed

  do i = 1, size(key)
    k = key(i)
    k = k * c1
    k = ishftc(k, 15)
    k = k * c2

    hash = ieor(hash, k)
    hash = ishftc(hash, 13)
    hash = 5 * hash + c3
  enddo

  ! NOTE: This is the point where the algorithm would handle trailing bytes.
  ! Since our arrays are comprised of 4 or 8 byte elements, we skip this part.

  hash = ieor(hash, 4*size(key))

  hash = ieor(hash, ishft(hash, -16))
  hash = hash * c4
  hash = ieor(hash, ishft(hash, -13))
  hash = hash * c5
  hash = ieor(hash, ishft(hash, -16))
end function murmurhash3_i32


!> Return the murmur3 hash for a 64-bit integer array.
function murmurhash3_i64(key, seed) result(hash)
  integer(int64), intent(in) :: key(:)
    !< Input array
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array

  integer(int32) :: ikey(2*size(key))

  hash = murmur_hash(transfer(key, ikey), seed=seed)
end function murmurhash3_i64


!> Return the murmur3 hash for a 32-bit real array.
function murmurhash3_r32(key, seed) result(hash)
  real(real32), intent(in) :: key
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array

  integer(int32) :: ikey(1)

  hash = murmur_hash(transfer(key, ikey), seed=seed)
end function murmurhash3_r32


!> Return the murmur3 hash for a 32-bit real array.
function murmurhash3_r32_1d(key, seed) result(hash)
  real(real32), intent(in) :: key(:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array

  integer(int32) :: ikey(size(key))

  hash = murmur_hash(transfer(key, ikey), seed=seed)
end function murmurhash3_r32_1d


!> Return the murmur3 hash for a 32-bit real 2D array.
function murmurhash3_r32_2d(key, seed) result(hash)
  real(real32), intent(in) :: key(:,:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array

  integer(int32) :: ikey(size(key))

  hash = murmur_hash(transfer(key, ikey), seed=seed)
end function murmurhash3_r32_2d


!> Return the murmur3 hash for a 32-bit real 3D array.
function murmurhash3_r32_3d(key, seed) result(hash)
  real(real32), intent(in) :: key(:,:,:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array

  integer(int32) :: ikey(size(key))

  hash = murmur_hash(transfer(key, ikey), seed=seed)
end function murmurhash3_r32_3d


!> Return the murmur3 hash for a 32-bit real 4D array.
function murmurhash3_r32_4d(key, seed) result(hash)
  real(real32), intent(in) :: key(:,:,:,:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array

  integer(int32) :: ikey(size(key))

  hash = murmur_hash(transfer(key, ikey), seed=seed)
end function murmurhash3_r32_4d


!> Return the murmur3 hash for a 64-bit real array.
function murmurhash3_r64(key, seed) result(hash)
  real(real64), intent(in) :: key
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array

  integer(int32) :: ikey(2)

  hash = murmur_hash(transfer(key, ikey), seed=seed)
end function murmurhash3_r64


!> Return the murmur3 hash for a 64-bit real array.
function murmurhash3_r64_1d(key, seed) result(hash)
  real(real64), intent(in) :: key(:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array

  integer(int32) :: ikey(2*size(key))

  hash = murmur_hash(transfer(key, ikey), seed=seed)
end function murmurhash3_r64_1d


!> Return the murmur3 hash for a 64-bit real 2D array.
function murmurhash3_r64_2d(key, seed) result(hash)
  real(real64), intent(in) :: key(:,:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array

  integer(int32) :: ikey(2*size(key))

  hash = murmur_hash(transfer(key, ikey), seed=seed)
end function murmurhash3_r64_2d


!> Return the murmur3 hash for a 64-bit real 3D array.
function murmurhash3_r64_3d(key, seed) result(hash)
  real(real64), intent(in) :: key(:,:,:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array

  integer(int32) :: ikey(2*size(key))

  hash = murmur_hash(transfer(key, ikey), seed=seed)
end function murmurhash3_r64_3d


!> Return the murmur3 hash for a 64-bit real 4D array.
function murmurhash3_r64_4d(key, seed) result(hash)
  real(real64), intent(in) :: key(:,:,:,:)
    !< Input array [arbitrary]
  integer(int32), intent(in), optional :: seed
    !< Hash seed
  integer(int32) :: hash
    !< Murmur hash of array

  integer(int32) :: ikey(2*size(key))

  hash = murmur_hash(transfer(key, ikey), seed=seed)
end function murmurhash3_r64_4d

end module MOM_murmur_hash
