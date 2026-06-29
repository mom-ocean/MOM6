! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Interfaces to non-domain-oriented communication subroutines, including the
!! MOM6 reproducing sums facility
module MOM_coms

use, intrinsic :: iso_fortran_env, only : int64
use MOM_coms_infra,    only : PE_here, root_PE, num_PEs, set_rootPE, Set_PElist, Get_PElist
use MOM_coms_infra,    only : broadcast, field_chksum, MOM_infra_init, MOM_infra_end
use MOM_coms_infra,    only : sum_across_PEs, max_across_PEs, min_across_PEs
use MOM_coms_infra,    only : all_across_PEs, any_across_PEs
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_coms_infra,    only : sync_PEs

implicit none ; private

public :: PE_here, root_PE, num_PEs, MOM_infra_init, MOM_infra_end
public :: sync_PEs
public :: broadcast, sum_across_PEs, min_across_PEs, max_across_PEs, field_chksum
public :: all_across_PEs, any_across_PEs
public :: set_PElist, Get_PElist, Set_rootPE
public :: reproducing_sum, reproducing_sum_EFP, EFP_sum_across_PEs, EFP_list_sum_across_PEs
public :: EFP_plus, EFP_minus, EFP_to_real, real_to_EFP, EFP_real_diff
public :: operator(+), operator(-), assignment(=)
public :: query_EFP_overflow_error, reset_EFP_overflow_error

integer, parameter :: accum_width = digits(1_int64)
  !< Accumulator width; total available bits for summation (excluding sign bit)
integer, parameter :: prec_width = 46
  !< Precision width; total bits for computed results
integer, parameter :: guard_width = accum_width - prec_width
  !< Number of guard bits reserved for carry overflow

! A sum of N points does N - 1 additions, which at most adds N - 1 carry bits.
! For G guard bits, the maximum value is 2**G - 1.  A summation of N values
! therefore requires that N - 1 <= 2**G - 1, or simply N <= 2**G.

integer, parameter :: max_summands = 2**guard_width
  !< Maximum number of summable points that can guarantee no carry overflow.
  !! Assumes that guard_bits is less than number of bits in a default integer.

integer(kind=int64), parameter :: prec = (2_int64)**prec_width
  !< EPF upper bound (exclusive).  For each EPF bin e(i), 0 <= e(i) < prec.

real, parameter :: r_prec = 2.**prec_width
  !< Real-value of prec [nondim]
real, parameter :: I_prec = 2.**(-prec_width)
  !< Inverse real-value of prec [nondim]

integer, parameter :: efp_digits = 6
  !< The number of base `prec` digits used to represent an EFP value.
real, parameter, dimension(efp_digits) :: &
    pr = [r_prec**2, r_prec, 1., r_prec**(-1), r_prec**(-2), r_prec**(-3)]
  !< An array of the real precision of each of the integers in arbitrary
  !! units [a]
real, parameter, dimension(efp_digits) :: &
    I_pr = [r_prec**(-2), r_prec**(-1), 1., r_prec, r_prec**2, r_prec**3]
  !< An array of the inverse of the real precision of each of the integers in
  !! arbitrary units [a-1]
real, parameter :: max_efp_float = pr(1) * real(huge(1_int64))
  !< The largest float with an EFP representation in arbitrary units [a].
  !! NOTE: Only the first bin can exceed precision, but is bounded by the
  !! largest signed integer.

logical :: overflow_error = .false.
  !< This becomes true if an overflow is encountered.
logical :: NaN_error = .false.
  !< This becomes true if a NaN is encountered.
logical :: debug = .false.
  !< Making this true enables debugging output.

! This module provides interfaces to the non-domain-oriented communication subroutines.

!> Find an accurate and order-invariant sum of a distributed 2d or 3d field, in some cases after
!! undoing the scaling of the input array and restoring that scaling in the returned value
interface reproducing_sum
  module procedure reproducing_sum_2d, reproducing_sum_3d
end interface reproducing_sum

!> Find an accurate and order-invariant sum of a distributed 2d field, returning the result
!! in the form of an extended fixed point value that can be converted back with EFP_to_real.
interface reproducing_sum_EFP
  module procedure reproducing_EFP_sum_2d
end interface reproducing_sum_EFP

!> Sum a value or 1-d array of values across processors, returning the sums in place
interface EFP_sum_across_PEs
  module procedure EFP_list_sum_across_PEs, EFP_val_sum_across_PEs
end interface EFP_sum_across_PEs

!> The Extended Fixed Point (EFP) type provides a public interface for doing sums
!! and taking differences with this type.
!!
!! The use of this type is documented in
!!   Hallberg, R. & A. Adcroft, 2014: An Order-invariant Real-to-Integer Conversion Sum.
!!   Parallel Computing, 40(5-6), doi:10.1016/j.parco.2014.04.007.
type, public :: EFP_type ; private
  integer(kind=int64), dimension(efp_digits) :: v !< The value in this type
end type EFP_type

!> Add two extended-fixed-point numbers
interface operator (+) ; module procedure EFP_plus  ; end interface
!> Subtract one extended-fixed-point number from another
interface operator (-) ; module procedure EFP_minus ; end interface
!> Copy the value of one extended-fixed-point number into another
interface assignment(=); module procedure EFP_assign ; end interface

contains

!> This subroutine uses a conversion to an integer representation of real numbers to give an
!! order-invariant sum of distributed 2-D arrays that reproduces across domain decomposition, with
!! the result returned as an extended fixed point type that can be converted back to a real number
!! using EFP_to_real.  This technique is described in Hallberg & Adcroft, 2014, Parallel Computing,
!! doi:10.1016/j.parco.2014.04.007.
function reproducing_EFP_sum_2d(array, isr, ier, jsr, jer, overflow_check, err, only_on_PE, unscale) result(EFP_sum)
  real, dimension(:,:),     intent(in)  :: array   !< The array to be summed in arbitrary units [a], or in
                                                   !! arbitrary scaled units [A ~> a] if unscale is present
  integer,        optional, intent(in)  :: isr     !< The starting i-index of the sum, noting
                                                   !! that the array indices starts at 1
  integer,        optional, intent(in)  :: ier     !< The ending i-index of the sum, noting
                                                   !! that the array indices starts at 1
  integer,        optional, intent(in)  :: jsr     !< The starting j-index of the sum, noting
                                                   !! that the array indices starts at 1
  integer,        optional, intent(in)  :: jer     !< The ending j-index of the sum, noting
                                                   !! that the array indices starts at 1
  logical,        optional, intent(in)  :: overflow_check !< If present and false, disable
                                                   !! checking for overflows in incremental results.
                                                   !! This can speed up calculations if the number
                                                   !! of values being summed is small enough
  integer,        optional, intent(out) :: err     !< If present, return an error code instead of
                                                   !! triggering any fatal errors directly from
                                                   !! this routine.
  logical,        optional, intent(in)  :: only_on_PE !< If present and true, do not do the sum
                                                   !! across processors, only reporting the local sum
  real,           optional, intent(in)  :: unscale !< A factor that is used to undo scaling of array before it is
                                                   !! summed, often to compensate for the scaling in [a A-1 ~> 1]
  type(EFP_type)                        :: EFP_sum !< The result in extended fixed point format

  !   This subroutine uses a conversion to an integer representation
  ! of real numbers to give order-invariant sums that will reproduce
  ! across PE count.  This idea comes from R. Hallberg and A. Adcroft.

  integer(kind=int64), dimension(efp_digits)  :: ints_sum
  integer(kind=int64) :: ival, prec_error
  real    :: rs ! The remaining value to add, in arbitrary units [a]
  real    :: max_mag_term ! A running maximum magnitude of the values in arbitrary units [a]
  real    :: descale    ! A local copy of unscale if it is present [a A-1 ~> 1] or 1
  logical :: over_check, do_sum_across_PEs, do_unscale
  character(len=256) :: mesg
  integer :: i, j, n, is, ie, js, je, sgn

  if (num_PEs() > max_summands) call MOM_error(FATAL, &
    "reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/num_PEs.")

  prec_error = huge(1_int64) / num_PEs()

  is = 1 ; ie = size(array,1) ; js = 1 ; je = size(array,2)
  if (present(isr)) then
    if (isr < is) call MOM_error(FATAL, "Value of isr too small in reproducing_EFP_sum_2d.")
    is = isr
  endif
  if (present(ier)) then
    if (ier > ie) call MOM_error(FATAL, "Value of ier too large in reproducing_EFP_sum_2d.")
    ie = ier
  endif
  if (present(jsr)) then
    if (jsr < js) call MOM_error(FATAL, "Value of jsr too small in reproducing_EFP_sum_2d.")
    js = jsr
  endif
  if (present(jer)) then
    if (jer > je) call MOM_error(FATAL, "Value of jer too large in reproducing_EFP_sum_2d.")
    je = jer
  endif

  over_check = .true. ; if (present(overflow_check)) over_check = overflow_check
  do_sum_across_PEs = .true. ; if (present(only_on_PE)) do_sum_across_PEs = .not.only_on_PE
  do_unscale = .false. ; if (present(unscale)) do_unscale = (unscale /= 1.0)
  descale = 1.0 ; if (do_unscale) descale = unscale

  overflow_error = .false. ; NaN_error = .false. ; max_mag_term = 0.0

  ints_sum(:) = 0
  if (over_check) then
    call increment_block_ints(array, is, ie, js, je, descale, ints_sum, &
        max_mag_term, prec_error)
  else
    do j=js,je ; do i=is,ie
      sgn = 1 ; if (array(i,j)<0.0) sgn = -1
      rs = abs(descale*array(i,j))
      do n=1,efp_digits
        ival = int(rs*I_pr(n), kind=int64)
        rs = rs - ival*pr(n)
        ints_sum(n) = ints_sum(n) + sgn*ival
      enddo
    enddo ; enddo
    call carry_overflow(ints_sum, prec_error)
  endif

  if (present(err)) then
    err = 0
    if (overflow_error) &
      err = err+2
    if (NaN_error) &
      err = err+4
    if (err > 0) then ; do n=1,efp_digits ; ints_sum(n) = 0 ; enddo ; endif
  else
    if (NaN_error) then
      call MOM_error(FATAL, "NaN in input field of reproducing_EFP_sum(_2d).")
    endif
    if (abs(max_mag_term) >= prec_error*pr(1)) then
      write(mesg, '(ES13.5)') max_mag_term
      call MOM_error(FATAL,"Overflow in reproducing_EFP_sum(_2d) conversion of "//trim(mesg))
    endif
    if (overflow_error) then
      call MOM_error(FATAL, "Overflow in reproducing_EFP_sum(_2d).")
    endif
  endif

  if (do_sum_across_PEs) call sum_across_PEs(ints_sum, efp_digits)

  call regularize_ints(ints_sum)

  EFP_sum%v(:) = ints_sum(:)

end function reproducing_EFP_sum_2d


!> This subroutine uses a conversion to an integer representation of real numbers to give an
!! order-invariant sum of distributed 2-D arrays that reproduces across domain decomposition.
!! This technique is described in Hallberg & Adcroft, 2014, Parallel Computing,
!! doi:10.1016/j.parco.2014.04.007.
function reproducing_sum_2d(array, isr, ier, jsr, jer, EFP_sum, reproducing, &
                            overflow_check, err, only_on_PE, unscale) result(sum)
  real, dimension(:,:),     intent(in)  :: array   !< The array to be summed in arbitrary units [a], or in
                                                   !! arbitrary scaled units [A ~> a] if unscale is present
  integer,        optional, intent(in)  :: isr     !< The starting i-index of the sum, noting
                                                   !! that the array indices starts at 1
  integer,        optional, intent(in)  :: ier     !< The ending i-index of the sum, noting
                                                   !! that the array indices starts at 1
  integer,        optional, intent(in)  :: jsr     !< The starting j-index of the sum, noting
                                                   !! that the array indices starts at 1
  integer,        optional, intent(in)  :: jer     !< The ending j-index of the sum, noting
                                                   !! that the array indices starts at 1
  type(EFP_type), optional, intent(out) :: EFP_sum !< The result in extended fixed point format
  logical,        optional, intent(in)  :: reproducing !< If present and false, do the sum
                                                !! using the naive non-reproducing approach
  logical,        optional, intent(in)  :: overflow_check !< If present and false, disable
                                                !! checking for overflows in incremental results.
                                                !! This can speed up calculations if the number
                                                !! of values being summed is small enough
  integer,        optional, intent(out) :: err  !< If present, return an error code instead of
                                                !! triggering any fatal errors directly from
                                                !! this routine.
  logical,        optional, intent(in)  :: only_on_PE !< If present and true, do not do the sum
                                                !! across processors, only reporting the local sum
  real,           optional, intent(in)  :: unscale !< A factor that is used to undo scaling of array before it is
                                                   !! summed, often to compensate for the scaling in [a A-1 ~> 1]
  real                                  :: sum     !< The sum of the values in array in the same
                                                   !! arbitrary units as array [a] or [A ~> a]

  ! Local variables
  integer(kind=int64), dimension(efp_digits)  :: ints_sum
  integer(kind=int64) :: prec_error
  real    :: rsum(1)    ! The running sum, in arbitrary units [a]
  real    :: descale    ! A local copy of unscale if it is present [a A-1 ~> 1] or 1
  real    :: I_unscale  ! The reciprocal of unscale [A a-1 ~> 1]
  logical :: repro, do_sum_across_PEs, do_unscale
  character(len=256) :: mesg
  type(EFP_type) :: EFP_val ! An extended fixed point version of the sum
  integer :: i, j, is, ie, js, je

  if (num_PEs() > max_summands) call MOM_error(FATAL, &
    "reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/num_PEs.")

  prec_error = huge(1_int64) / num_PEs()

  is = 1 ; ie = size(array,1) ; js = 1 ; je = size(array,2)
  if (present(isr)) then
    if (isr < is) call MOM_error(FATAL, "Value of isr too small in reproducing_sum_2d.")
    is = isr
  endif
  if (present(ier)) then
    if (ier > ie) call MOM_error(FATAL, "Value of ier too large in reproducing_sum_2d.")
    ie = ier
  endif
  if (present(jsr)) then
    if (jsr < js) call MOM_error(FATAL, "Value of jsr too small in reproducing_sum_2d.")
    js = jsr
  endif
  if (present(jer)) then
    if (jer > je) call MOM_error(FATAL, "Value of jer too large in reproducing_sum_2d.")
    je = jer
  endif

  repro = .true. ; if (present(reproducing)) repro = reproducing
  do_sum_across_PEs = .true. ; if (present(only_on_PE)) do_sum_across_PEs = .not.only_on_PE
  do_unscale = .false. ; if (present(unscale)) do_unscale = (unscale /= 1.0)
  descale = 1.0 ;  I_unscale = 1.0
  if (do_unscale) then
    descale = unscale
    if (abs(unscale) > 0.0) I_unscale = 1.0 / unscale
  endif

  if (repro) then
    EFP_val = reproducing_EFP_sum_2d(array, isr, ier, jsr, jer, overflow_check, err, only_on_PE, unscale)
    sum = ints_to_real(EFP_val%v) * I_unscale
    if (present(EFP_sum)) EFP_sum = EFP_val
    if (debug) ints_sum(:) = EFP_sum%v(:)
  else
    rsum(1) = 0.0
    do j=js,je ; do i=is,ie
      rsum(1) = rsum(1) + descale*array(i,j)
    enddo ; enddo
    if (do_sum_across_PEs) call sum_across_PEs(rsum,1)
    sum = rsum(1) * I_unscale

    if (present(err)) then ; err = 0 ; endif

    if (debug .or. present(EFP_sum)) then
      overflow_error = .false.
      ints_sum = real_to_ints(sum, prec_error, overflow_error)
      if (overflow_error) then
        if (present(err)) then
          err = err + 2
        else
          write(mesg, '(ES13.5)') sum
          call MOM_error(FATAL,"Repro_sum_2d: Overflow in real_to_ints conversion of "//trim(mesg))
        endif
      endif
    endif
    if (present(EFP_sum)) EFP_sum%v(:) = ints_sum(:)
  endif

  if (debug) then
    write(mesg,'("2d RS: ", ES24.16, 6 Z17.16)') sum*descale, ints_sum(1:efp_digits)
    call MOM_mesg(mesg, 3)
  endif

end function reproducing_sum_2d

!> This subroutine uses a conversion to an integer representation of real numbers to give an
!! order-invariant sum of distributed 3-D arrays that reproduces across domain decomposition.
!! This technique is described in Hallberg & Adcroft, 2014, Parallel Computing,
!! doi:10.1016/j.parco.2014.04.007.
function reproducing_sum_3d(array, isr, ier, jsr, jer, sums, EFP_sum, EFP_lay_sums, err, only_on_PE, unscale) &
                            result(sum)
  real, dimension(:,:,:),       intent(in)  :: array   !< The array to be summed in arbitrary units [a], or in
                                                       !! arbitrary scaled units [A ~> a] if unscale is present
  integer,            optional, intent(in)  :: isr     !< The starting i-index of the sum, noting
                                                       !! that the array indices starts at 1
  integer,            optional, intent(in)  :: ier     !< The ending i-index of the sum, noting
                                                       !! that the array indices starts at 1
  integer,            optional, intent(in)  :: jsr     !< The starting j-index of the sum, noting
                                                       !! that the array indices starts at 1
  integer,            optional, intent(in)  :: jer     !< The ending j-index of the sum, noting
                                                       !! that the array indices starts at 1
  real, dimension(:), optional, intent(out) :: sums    !< The sums by vertical layer in the same
                                                       !! arbitrary units as array [a] or [A ~> a]
  type(EFP_type),     optional, intent(out) :: EFP_sum !< The result in extended fixed point format
  type(EFP_type), dimension(:), &
                      optional, intent(out) :: EFP_lay_sums !< The sums by vertical layer in EFP format
  integer,            optional, intent(out) :: err     !< If present, return an error code instead of
                                                       !! triggering any fatal errors directly from
                                                       !! this routine.
  logical,            optional, intent(in)  :: only_on_PE !< If present and true, do not do the sum
                                                       !! across processors, only reporting the local sum
  real,               optional, intent(in)  :: unscale !< A factor that is used to undo scaling of array before it is
                                                       !! summed, often to compensate for the scaling in [a A-1 ~> 1]
  real                                      :: sum     !< The sum of the values in array in the same
                                                       !! arbitrary units as array [a] or [A ~> a]

  ! Local variables
  real    :: val ! The real number that is extracted in arbitrary units [a]
  real    :: max_mag_term ! A running maximum magnitude of the val's in arbitrary units [a]
  real    :: descale    ! A local copy of unscale if it is present [a A-1 ~> 1] or 1
  real    :: I_unscale  ! The Adcroft reciprocal of unscale [A a-1 ~> 1]
  integer(kind=int64), dimension(efp_digits)  :: ints_sum
  integer(kind=int64), dimension(efp_digits,size(array,3))  :: ints_sums
  integer(kind=int64) :: prec_error
  character(len=256) :: mesg
  logical :: do_sum_across_PEs, do_unscale
  integer :: i, j, k, is, ie, js, je, ke, isz, jsz, n

  if (num_PEs() > max_summands) call MOM_error(FATAL, &
    "reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/num_PEs.")

  prec_error = huge(1_int64) / num_PEs()
  max_mag_term = 0.0

  is = 1 ; ie = size(array,1) ; js = 1 ; je = size(array,2) ; ke = size(array,3)
  if (present(isr)) then
    if (isr < is) call MOM_error(FATAL, "Value of isr too small in reproducing_sum(_3d).")
    is = isr
  endif
  if (present(ier)) then
    if (ier > ie) call MOM_error(FATAL, "Value of ier too large in reproducing_sum(_3d).")
    ie = ier
  endif
  if (present(jsr)) then
    if (jsr < js) call MOM_error(FATAL, "Value of jsr too small in reproducing_sum(_3d).")
    js = jsr
  endif
  if (present(jer)) then
    if (jer > je) call MOM_error(FATAL, "Value of jer too large in reproducing_sum(_3d).")
    je = jer
  endif
  jsz = je+1-js ; isz = ie+1-is

  do_sum_across_PEs = .true. ; if (present(only_on_PE)) do_sum_across_PEs = .not.only_on_PE
  do_unscale = .false. ; if (present(unscale)) do_unscale = (unscale /= 1.0)
  descale = 1.0 ; if (do_unscale) descale = unscale

  if (present(sums) .or. present(EFP_lay_sums)) then
    if (present(sums)) then ; if (size(sums) < ke) then
      call MOM_error(FATAL, "Sums is smaller than the vertical extent of array in reproducing_sum(_3d).")
    endif ; endif
    if (present(EFP_lay_sums)) then ; if (size(EFP_lay_sums) < ke) then
      call MOM_error(FATAL, "Sums is smaller than the vertical extent of array in reproducing_sum(_3d).")
    endif ; endif

    overflow_error = .false. ; NaN_error = .false. ; max_mag_term = 0.0

    ints_sums(:,:) = 0
    do k=1,ke
      call increment_block_ints(array(:,:,k), is, ie, js, je, descale, &
          ints_sums(:,k), max_mag_term, prec_error)
    enddo

    if (present(err)) then
      err = 0
      if (abs(max_mag_term) >= prec_error*pr(1)) err = err+1
      if (overflow_error) err = err+2
      if (NaN_error) err = err+2
      if (err > 0) then ; do k=1,ke ; do n=1,efp_digits ; ints_sums(n,k) = 0 ; enddo ; enddo ; endif
    else
      if (NaN_error) call MOM_error(FATAL, "NaN in input field of reproducing_sum(_3d).")
      if (abs(max_mag_term) >= prec_error*pr(1)) then
        write(mesg, '(ES13.5)') max_mag_term
        call MOM_error(FATAL,"Overflow in reproducing_sum(_3d) conversion of "//trim(mesg))
      endif
      if (overflow_error) call MOM_error(FATAL, "Overflow in reproducing_sum(_3d).")
    endif

    if (do_sum_across_PEs) call sum_across_PEs(ints_sums(:,1:ke), efp_digits*ke)

    sum = 0.0
    do k=1,ke
      call regularize_ints(ints_sums(:,k))
      val = ints_to_real(ints_sums(:,k))
      if (present(sums)) sums(k) = val
      sum = sum + val
    enddo
    if (present(EFP_lay_sums)) then ; do k=1,ke
      EFP_lay_sums(k)%v(:) = ints_sums(:,k)
    enddo ; endif

    if (present(EFP_sum)) then
      EFP_sum%v(:) = 0
      do k=1,ke ; call increment_ints(EFP_sum%v(:), ints_sums(:,k)) ; enddo
    endif

    if (debug) then
      do n=1,efp_digits ; ints_sum(n) = 0 ; enddo
      do k=1,ke ; do n=1,efp_digits ; ints_sum(n) = ints_sum(n) + ints_sums(n,k) ; enddo ; enddo
      write(mesg,'("3D RS: ", ES24.16, 6 Z17.16)') sum, ints_sum(1:efp_digits)
      call MOM_mesg(mesg, 3)
    endif
  else
    overflow_error = .false. ; NaN_error = .false. ; max_mag_term = 0.0

    ints_sum(:) = 0
    do k=1,ke
      call increment_block_ints(array(:,:,k), is, ie, js, je, descale, &
          ints_sum, max_mag_term, prec_error)
    enddo

    if (present(err)) then
      err = 0
      if (abs(max_mag_term) >= prec_error*pr(1)) err = err+1
      if (overflow_error) err = err+2
      if (NaN_error) err = err+2
      if (err > 0) then ; do n=1,efp_digits ; ints_sum(n) = 0 ; enddo ; endif
    else
      if (NaN_error) call MOM_error(FATAL, "NaN in input field of reproducing_sum(_3d).")
      if (abs(max_mag_term) >= prec_error*pr(1)) then
        write(mesg, '(ES13.5)') max_mag_term
        call MOM_error(FATAL,"Overflow in reproducing_sum(_3d) conversion of "//trim(mesg))
      endif
      if (overflow_error) call MOM_error(FATAL, "Overflow in reproducing_sum(_3d).")
    endif

    if (do_sum_across_PEs) call sum_across_PEs(ints_sum, efp_digits)

    call regularize_ints(ints_sum)
    sum = ints_to_real(ints_sum)

    if (present(EFP_sum)) EFP_sum%v(:) = ints_sum(:)

    if (debug) then
      write(mesg,'("3d RS: ", ES24.16, 6 Z17.16)') sum, ints_sum(1:efp_digits)
      call MOM_mesg(mesg, 3)
    endif
  endif

  if (do_unscale) then
    ! Revise the sum to restore the scaling of input array before it is returned
    I_unscale = 0.0 ; if (abs(unscale) > 0.0) I_unscale = 1.0 / unscale
    sum = sum * I_unscale
    if (present(sums)) then
      do k=1,ke ; sums(k) = sums(k) * I_unscale ; enddo
    endif
  endif

end function reproducing_sum_3d

!> Convert a real number into the array of integers constitute its extended-fixed-point representation
function real_to_ints(r, prec_error, overflow) result(ints)
  real,                      intent(in) :: r  !< The real number being converted in arbitrary units [a]
  integer(kind=int64), optional, intent(in) :: prec_error  !< The PE-count dependent precision of the
                                              !! integers that is safe from overflows during global
                                              !! sums.  This will be larger than the compile-time
                                              !! precision parameter, and is used to detect overflows.
  logical,         optional, intent(inout) :: overflow !< Returns true if the conversion is being
                                              !! done on a value that is too large to be represented
  integer(kind=int64), dimension(efp_digits)  :: ints

  !   This subroutine converts a real number to an equivalent representation
  ! using several long integers.

  ! Local variables
  real :: rs  ! The remaining value to add, in arbitrary units [a]
  character(len=80) :: mesg
  integer(kind=int64) :: ival, prec_err
  integer :: sgn, i

  prec_err = prec ; if (present(prec_error)) prec_err = prec_error
  ints(:) = 0
  if ((r >= 1e30) .eqv. (r < 1e30)) then ; NaN_error = .true. ; return ; endif

  sgn = 1 ; if (r<0.0) sgn = -1
  rs = abs(r)

  if (present(overflow)) then
    if (.not.(rs < prec_err*pr(1))) overflow = .true.
    if ((r >= 1e30) .eqv. (r < 1e30)) overflow = .true.
  elseif (.not.(rs < prec_err*pr(1))) then
    write(mesg, '(ES13.5)') r
    call MOM_error(FATAL,"Overflow in real_to_ints conversion of "//trim(mesg))
  endif

  do i=1,efp_digits
    ival = int(rs*I_pr(i), kind=int64)
    rs = rs - ival*pr(i)
    ints(i) = sgn*ival
  enddo

end function real_to_ints

!> Convert the array of integers that constitute an extended-fixed-point
!! representation into a real number
function ints_to_real(ints) result(r)
  integer(kind=int64), dimension(efp_digits), intent(in) :: ints !< The array of EFP integers
  real :: r  ! The real number that is extracted in arbitrary units [a]
  ! This subroutine reverses the conversion in real_to_ints.

  integer :: i

  r = 0.0
  do i=1,efp_digits ; r = r + pr(i)*ints(i) ; enddo
end function ints_to_real

!> Increment an array of integers that constitutes an extended-fixed-point
!! representation with a another EFP number
subroutine increment_ints(int_sum, int2, prec_error)
  integer(kind=int64), dimension(efp_digits), intent(inout) :: int_sum !< The array of EFP integers being incremented
  integer(kind=int64), dimension(efp_digits), intent(in)    :: int2    !< The array of EFP integers being added
  integer(kind=int64), optional,      intent(in)    :: prec_error !< The PE-count dependent precision of the
                                              !! integers that is safe from overflows during global
                                              !! sums.  This will be larger than the compile-time
                                              !! precision parameter, and is used to detect overflows.

  ! This subroutine increments a number with another, both using the integer
  ! representation in real_to_ints.
  integer :: i

  do i=efp_digits,2,-1
    int_sum(i) = int_sum(i) + int2(i)
    ! Carry the local overflow.
    if (int_sum(i) > prec) then
      int_sum(i) = int_sum(i) - prec
      int_sum(i-1) = int_sum(i-1) + 1
    elseif (int_sum(i) < -prec) then
      int_sum(i) = int_sum(i) + prec
      int_sum(i-1) = int_sum(i-1) - 1
    endif
  enddo
  int_sum(1) = int_sum(1) + int2(1)
  if (present(prec_error)) then
    if (abs(int_sum(1)) > prec_error) overflow_error = .true.
  else
    if (abs(int_sum(1)) > prec) overflow_error = .true.
  endif

end subroutine increment_ints


!> Sum the elements of an array in EFP form and append the result to an
!! existing EFP array.
subroutine increment_block_ints(array, is, ie, js, je, descale, ints_sum, &
    max_mag_term, prec_error)
  real, intent(in) :: array(:,:)
    !< The field being added, in arbitrary units [A ~> a]
  integer, intent(in) :: is
    !< Start i-index of the summed domain
  integer, intent(in) :: ie
    !< End i-index of the summed domain
  integer, intent(in) :: js
    !< Start j-index of the summed domain
  integer, intent(in) :: je
    !< End j-index of the summed domain
  real, intent(in) :: descale
    !< Factor to descale array to physical value [a A-1 ~> 1]
  integer(kind=int64), intent(inout) :: ints_sum(efp_digits)
    !< The array of EFP integers being incremented
  real, intent(inout) :: max_mag_term
    !< A running maximum magnitude of the r's, in arbitrary units [a]
  integer(kind=int64), intent(in) :: prec_error
    !< The maximum resolvable value for a given number of PEs

  integer :: i, j, ib, jb, ibs, ibe, jbs, jbe
    ! Loop indices
  integer :: b
    ! Block counter
  integer :: ni, nj
    ! Array summation domain size along each axis
  integer :: isize_max
    ! Largest block size in i.  Typically equal to ni
  integer :: jsize
    ! Number of j-rows per block.
  integer :: nblocks, niblocks, njblocks
    ! Number of total blocks, and number of blocks in i and j
  integer(kind=int64) :: e(efp_digits)
    ! The EPF representation of each array element
  integer(kind=int64) :: block_sum(efp_digits), array_sum(efp_digits)
    ! The cumulant per-block and total array EFP sums
  real :: r, rmag
    ! Local array element value and its magnitude [a]
  real :: max_pos, max_neg, block_max_pos, block_max_neg
    ! Largest positive and negative values (whole array and per-block) used to
    ! find the largest maximum magnitude of array in a thread-safe manner [a]
  integer :: inan, iovf, lnan, lovf
    ! Thread-safe tracking of NaN and overflow state
  integer :: max_sum_count
    ! The total number of local sum operations that ensures no carry overflow

  max_pos = max(0., max_mag_term)
  max_neg = max(0., -max_mag_term)
  inan = 0 ; iovf = 0

  ! Reduce the maximum number of summations to account for the cumulant
  ! summations of array_sum and ints_sum.
  max_sum_count = max_summands - 2

  ! Get the compute domain size
  ni = ie - is + 1
  nj = je - js + 1

  ! Partition in i so that the widest i-slice fits within max_sum_count.
  niblocks = (ni + max_sum_count - 1) / max_sum_count
         ! = ⌈ni / max_sum_count⌉

  ! NOTE: niblocks is typically one, since default max_sum_count is ~130k.

  ! For a balanced i-partition, the number of i-points per block is either
  !   ⌊ni / niblocks⌋ or ⌈ni / niblocks⌉.  Use the upper bound to find jsize.

  isize_max = (ni + niblocks - 1) / niblocks
          ! = ⌈ni / niblocks⌉

  ! Set jsize so that the widest i-slice times the number of j-rows does not
  !   exceed max_sum_count.
  jsize = max_sum_count / isize_max
      ! = ⌊max_sum_count / isize_max⌋

  ! Choose enough j-blocks so that no j-block has more than jsize rows.
  njblocks = (nj + jsize - 1) / jsize
         ! = ⌈nj / jsize⌉

  nblocks = niblocks * njblocks

  ! Abort if the number of blocks also exceeds the carry-bit summation limit.
  ! For default settings, this would be over 17 billion points per PE.
  if (nblocks > max_sum_count) call MOM_error(FATAL, &
      "reproducing sum: Number of blocks exceeds summmation carry limit.")

  array_sum(:) = 0

  do jb=1,njblocks ; do ib=1,niblocks
    ! Use evenly distributed blocks, either ⌊n / nblocks⌋ or ⌈n / nblocks⌉.
    jbs = js + ((jb - 1) * nj) / njblocks
    jbe = js + (jb * nj) / njblocks - 1

    ibs = is + ((ib - 1) * ni) / niblocks
    ibe = is + (ib * ni) / niblocks - 1

    block_sum(:) = 0
    block_max_pos = 0.
    block_max_neg = 0.

    ! Compute the sum of each block
    do j=jbs,jbe ; do i=ibs,ibe

      ! Convert array(i,j) to EFP form
      r = descale * array(i,j)
      call efp_decompose(r, e, rmag, lnan, lovf)

      ! Verify that the conversion was completed
      inan = max(inan, lnan)
      iovf = max(iovf, lovf)

      if (r >= 0.) then
        if (rmag > block_max_pos) block_max_pos = rmag
      else
        if (rmag > block_max_neg) block_max_neg = rmag
      endif

      ! Add the EFP result (including potential carry bits)
      block_sum(:) = block_sum(:) + e(:)
    enddo ; enddo

    array_sum(:) = array_sum(:) + block_sum(:)

    ! Redistribute carry bits across bins
    ! For the final pass (or single pass) this is handled by ints_sum.
    b = (jb - 1) * niblocks + ib
    if (b < nblocks) call carry_overflow(array_sum, prec_error)

    ! Update maximum magnitudes
    max_pos = max(max_pos, block_max_pos)
    max_neg = max(max_neg, block_max_neg)
  enddo ; enddo

  ! Finally, apply the cumulant result
  ints_sum(:) = ints_sum(:) + array_sum(:)

  ! Redistribute carry bits to normalize the final result.
  call carry_overflow(ints_sum, prec_error)

  ! Extract the maximum value while preserving sign (NOTE: ties to positive)
  if (max_pos >= max_neg) then
    max_mag_term = max_pos
  else
    max_mag_term = -max_neg
  endif

  ! Transfer error/warning signals to module flags
  if (inan /= 0) NaN_error = .true.
  if (iovf /= 0) overflow_error = .true.
end subroutine increment_block_ints


!> Decompose one real into its 6 signed EFP bin contributions.  NaNs and
!! overflows are reported by flags, rather than the module-level error
!! logicals, so that the routine is free of side effects.
pure subroutine efp_decompose(r, e, rmag, is_nan, is_ovf)
  real, intent(in)  :: r
    !< The real number being decomposed [a]
  integer(kind=int64), intent(out) :: e(efp_digits)
    !< Signed contribution to EFP bins
  real, intent(out) :: rmag
    !< Equals abs(r), or 0 if r is NaN/Inf [a]
  integer, intent(out) :: is_nan
    !< Equals 1 if r is a NaN or Inf, else 0
  integer, intent(out) :: is_ovf
    !< Equals 1 if abs(r) has no EFP representation, else 0

  real :: rs
    ! The remaining value to add, in arbitrary units [a]
  integer(kind=int64) :: ival
  integer :: sgn
  integer :: n

  e(:) = 0
  rmag = 0.0 ; is_nan = 0 ; is_ovf = 0

  if ((r >= 1e30) .eqv. (r < 1e30)) then
    is_nan = 1
    return
  endif

  sgn = 1
  if (r < 0.0) sgn = -1

  rs = abs(r) ; rmag = rs

  ! Abort if the number has no EFP representation
  if (rs > max_efp_float) then
    is_ovf = 1
    return
  endif

  do n=1,efp_digits
    ival = int(rs * I_pr(n), kind=int64)
    rs = rs - ival * pr(n)
    e(n) = sgn * ival
  enddo
end subroutine efp_decompose


!> This subroutine handles carrying of the overflow.
subroutine carry_overflow(int_sum, prec_error)
  integer(kind=int64), dimension(efp_digits), intent(inout) :: int_sum  !< The array of EFP integers being
                                              !! modified by carries, but without changing value.
  integer(kind=int64),                intent(in)    :: prec_error  !< The PE-count dependent precision of the
                                              !! integers that is safe from overflows during global
                                              !! sums.  This will be larger than the compile-time
                                              !! precision parameter, and is used to detect overflows.

  ! This subroutine handles carrying of the overflow.
  integer :: i, num_carry

  do i=efp_digits,2,-1 ; if (abs(int_sum(i)) >= prec) then
    num_carry = int(int_sum(i) * I_prec)
    int_sum(i) = int_sum(i) - num_carry*prec
    int_sum(i-1) = int_sum(i-1) + num_carry
  endif ; enddo
  if (abs(int_sum(1)) > prec_error) then
    overflow_error = .true.
  endif

end subroutine carry_overflow

!> This subroutine carries the overflow, and then makes sure that
!! all integers are of the same sign as the overall value.
subroutine regularize_ints(int_sum)
  integer(kind=int64), dimension(efp_digits), &
    intent(inout) :: int_sum !< The array of integers being modified to take a
                             !! regular form with all integers of the same sign,
                             !! but without changing value.

  ! This subroutine carries the overflow, and then makes sure that
  ! all integers are of the same sign as the overall value.
  logical :: positive
  integer :: i, num_carry

  do i=efp_digits,2,-1 ; if (abs(int_sum(i)) >= prec) then
    num_carry = int(int_sum(i) * I_prec)
    int_sum(i) = int_sum(i) - num_carry*prec
    int_sum(i-1) = int_sum(i-1) + num_carry
  endif ; enddo

  ! Determine the sign of the final number.
  positive = .true.
  do i=1,efp_digits
    if (abs(int_sum(i)) > 0) then
      if (int_sum(i) < 0) positive = .false.
      exit
    endif
  enddo

  if (positive) then
    do i=efp_digits,2,-1 ; if (int_sum(i) < 0) then
      int_sum(i) = int_sum(i) + prec
      int_sum(i-1) = int_sum(i-1) - 1
    endif ; enddo
  else
    do i=efp_digits,2,-1 ; if (int_sum(i) > 0) then
      int_sum(i) = int_sum(i) - prec
      int_sum(i-1) = int_sum(i-1) + 1
    endif ; enddo
  endif

end subroutine regularize_ints

!> Returns the status of the module's error flag
function query_EFP_overflow_error()
  logical :: query_EFP_overflow_error
  query_EFP_overflow_error = overflow_error
end function query_EFP_overflow_error

!> Reset the module's error flag to false
subroutine reset_EFP_overflow_error()
  overflow_error = .false.
end subroutine reset_EFP_overflow_error

!> Add two extended-fixed-point numbers
function EFP_plus(EFP1, EFP2)
  type(EFP_type)             :: EFP_plus !< The result in extended fixed point format
  type(EFP_type), intent(in) :: EFP1 !< The first extended fixed point number
  type(EFP_type), intent(in) :: EFP2 !< The second extended fixed point number

  EFP_plus = EFP1

  call increment_ints(EFP_plus%v(:), EFP2%v(:))
end function EFP_plus

!> Subract one extended-fixed-point number from another
function EFP_minus(EFP1, EFP2)
  type(EFP_type)             :: EFP_minus !< The result in extended fixed point format
  type(EFP_type), intent(in) :: EFP1 !< The first extended fixed point number
  type(EFP_type), intent(in) :: EFP2 !< The extended fixed point number being
                        !! subtracted from the first extended fixed point number
  integer :: i

  do i=1,efp_digits ; EFP_minus%v(i) = -1*EFP2%v(i) ; enddo

  call increment_ints(EFP_minus%v(:), EFP1%v(:))
end function EFP_minus

!> Copy one extended-fixed-point number into another
subroutine EFP_assign(EFP1, EFP2)
  type(EFP_type), intent(out) :: EFP1 !< The recipient extended fixed point number
  type(EFP_type), intent(in)  :: EFP2 !< The source extended fixed point number
  integer i
  ! This subroutine assigns all components of the extended fixed point type
  ! variable on the RHS (EFP2) to the components of the variable on the LHS
  ! (EFP1).

  do i=1,efp_digits ; EFP1%v(i) = EFP2%v(i) ; enddo
end subroutine EFP_assign

!> Return the real number that an extended-fixed-point number corresponds with
function EFP_to_real(EFP1)
  type(EFP_type), intent(inout) :: EFP1 !< The extended fixed point number being converted
  real :: EFP_to_real  !< The real version of the number in arbitrary units [a]

  call regularize_ints(EFP1%v)
  EFP_to_real = ints_to_real(EFP1%v)
end function EFP_to_real

!> Take the difference between two extended-fixed-point numbers (EFP1 - EFP2)
!! and return the result as a real number
function EFP_real_diff(EFP1, EFP2)
  type(EFP_type), intent(in) :: EFP1  !< The first extended fixed point number
  type(EFP_type), intent(in) :: EFP2  !< The extended fixed point number being
                        !! subtracted from the first extended fixed point number
  real :: EFP_real_diff !< The real result in arbitrary units [a]

  type(EFP_type)             :: EFP_diff

  EFP_diff = EFP1 - EFP2
  EFP_real_diff = EFP_to_real(EFP_diff)

end function EFP_real_diff

!> Return the extended-fixed-point number that a real number corresponds with
function real_to_EFP(val, overflow)
  real,              intent(in)    :: val !< The real number being converted in arbitrary units [a]
  logical, optional, intent(inout) :: overflow !< Returns true if the conversion is being
                                          !! done on a value that is too large to be represented
  type(EFP_type) :: real_to_EFP

  logical :: over
  character(len=80) :: mesg

  if (present(overflow)) then
    real_to_EFP%v(:) = real_to_ints(val, overflow=overflow)
  else
    over = .false.
    real_to_EFP%v(:) = real_to_ints(val, overflow=over)
    if (over) then
      write(mesg, '(ES13.5)') val
      call MOM_error(FATAL,"Overflow in real_to_EFP conversion of "//trim(mesg))
    endif
  endif

end function real_to_EFP

!>   This subroutine does a sum across PEs of a list of EFP variables,
!! returning the sums in place, with all overflows carried.
subroutine EFP_list_sum_across_PEs(EFPs, nval, errors)
  type(EFP_type), dimension(:), &
              intent(inout) :: EFPs   !< The list of extended fixed point numbers
                                      !! being summed across PEs.
  integer,    intent(in)    :: nval   !< The number of values being summed.
  logical, dimension(:), &
           optional, intent(out)   :: errors !< A list of error flags for each sum

  !   This subroutine does a sum across PEs of a list of EFP variables,
  ! returning the sums in place, with all overflows carried.

  integer(kind=int64), dimension(efp_digits,nval) :: ints
  integer(kind=int64) :: prec_error
  logical :: error_found
  character(len=256) :: mesg
  integer :: i, n

  if (num_PEs() > max_summands) call MOM_error(FATAL, &
    "reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/num_PEs.")

  prec_error = huge(1_int64) / num_PEs()

  ! overflow_error is an overflow error flag for the whole module.
  overflow_error = .false. ; error_found = .false.

  do i=1,nval ; do n=1,efp_digits ; ints(n,i) = EFPs(i)%v(n) ; enddo ; enddo

  call sum_across_PEs(ints(:,:), efp_digits*nval)

  if (present(errors)) errors(:) = .false.
  do i=1,nval
    overflow_error = .false.
    call carry_overflow(ints(:,i), prec_error)
    do n=1,efp_digits ; EFPs(i)%v(n) = ints(n,i) ; enddo
    if (present(errors)) errors(i) = overflow_error
    if (overflow_error) then
      write (mesg,'("EFP_list_sum_across_PEs error at ",i0," val was ",ES12.6, ", prec_error = ",ES12.6)') &
             i, EFP_to_real(EFPs(i)), real(prec_error)
      call MOM_error(WARNING, mesg)
    endif
    error_found = error_found .or. overflow_error
  enddo
  if (error_found .and. .not.(present(errors))) then
    call MOM_error(FATAL, "Overflow in EFP_list_sum_across_PEs.")
  endif

end subroutine EFP_list_sum_across_PEs

!>   This subroutine does a sum across PEs of an EFP variable,
!! returning the sums in place, with all overflows carried.
subroutine EFP_val_sum_across_PEs(EFP, error)
  type(EFP_type),  intent(inout) :: EFP   !< The extended fixed point numbers
                                          !! being summed across PEs.
  logical, optional, intent(out) :: error !< An error flag for this sum

  !   This subroutine does a sum across PEs of a list of EFP variables,
  ! returning the sums in place, with all overflows carried.

  integer(kind=int64), dimension(efp_digits) :: ints
  integer(kind=int64) :: prec_error
  logical :: error_found
  character(len=256) :: mesg
  integer :: n

  if (num_PEs() > max_summands) call MOM_error(FATAL, &
    "reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/num_PEs.")

  prec_error = huge(1_int64) / num_PEs()

  ! overflow_error is an overflow error flag for the whole module.
  overflow_error = .false. ; error_found = .false.

  do n=1,efp_digits ; ints(n) = EFP%v(n) ; enddo

  call sum_across_PEs(ints(:), efp_digits)

  if (present(error)) error = .false.

  overflow_error = .false.
  call carry_overflow(ints(:), prec_error)
  do n=1,efp_digits ; EFP%v(n) = ints(n) ; enddo
  if (present(error)) error = overflow_error
  if (overflow_error) then
    write (mesg,'("EFP_val_sum_across_PEs error val was ",ES12.6, ", prec_error = ",ES12.6)') &
           EFP_to_real(EFP), real(prec_error)
    call MOM_error(WARNING, mesg)
  endif
  error_found = error_found .or. overflow_error

  if (error_found .and. .not.(present(error))) then
    call MOM_error(FATAL, "Overflow in EFP_val_sum_across_PEs.")
  endif

end subroutine EFP_val_sum_across_PEs

end module MOM_coms
