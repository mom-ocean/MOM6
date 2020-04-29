!> Interfaces to non-domain-oriented communication subroutines, including the
!! MOM6 reproducing sums facility
module MOM_coms

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use fms_mod, only : fms_end, MOM_infra_init => fms_init
use memutils_mod, only : print_memuse_stats
use mpp_mod, only : PE_here => mpp_pe, root_PE => mpp_root_pe, num_PEs => mpp_npes
use mpp_mod, only : Set_PElist => mpp_set_current_pelist, Get_PElist => mpp_get_current_pelist
use mpp_mod, only : broadcast => mpp_broadcast
use mpp_mod, only : sum_across_PEs => mpp_sum, max_across_PEs => mpp_max, min_across_PEs => mpp_min

implicit none ; private

public :: PE_here, root_PE, num_PEs, MOM_infra_init, MOM_infra_end
public :: broadcast, sum_across_PEs, min_across_PEs, max_across_PEs
public :: reproducing_sum, reproducing_sum_EFP, EFP_sum_across_PEs, EFP_list_sum_across_PEs
public :: EFP_plus, EFP_minus, EFP_to_real, real_to_EFP, EFP_real_diff
public :: operator(+), operator(-), assignment(=)
public :: query_EFP_overflow_error, reset_EFP_overflow_error
public :: Set_PElist, Get_PElist

! This module provides interfaces to the non-domain-oriented communication subroutines.

integer(kind=8), parameter :: prec=2_8**46 !< The precision of each integer.
real, parameter :: r_prec=2.0**46  !< A real version of prec.
real, parameter :: I_prec=1.0/(2.0**46) !< The inverse of prec.
integer, parameter :: max_count_prec=2**(63-46)-1
                              !< The number of values that can be added together
                              !! with the current value of prec before there will
                              !! be roundoff problems.

integer, parameter :: ni=6    !< The number of long integers to use to represent
                              !< a real number.
real, parameter, dimension(ni) :: &
  pr = (/ r_prec**2, r_prec, 1.0, 1.0/r_prec, 1.0/r_prec**2, 1.0/r_prec**3 /)
    !< An array of the real precision of each of the integers
real, parameter, dimension(ni) :: &
  I_pr = (/ 1.0/r_prec**2, 1.0/r_prec, 1.0, r_prec, r_prec**2, r_prec**3 /)
    !< An array of the inverse of the real precision of each of the integers
real, parameter :: max_efp_float = pr(1) * (2.**63 - 1.)
                              !< The largest float with an EFP representation.
                              !! NOTE: Only the first bin can exceed precision,
                              !! but is bounded by the largest signed integer.

logical :: overflow_error = .false. !< This becomes true if an overflow is encountered.
logical :: NaN_error = .false.      !< This becomes true if a NaN is encountered.
logical :: debug = .false.          !< Making this true enables debugging output.

!> Find an accurate and order-invariant sum of a distributed 2d or 3d field
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
  integer(kind=8), dimension(ni) :: v !< The value in this type
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
function reproducing_EFP_sum_2d(array, isr, ier, jsr, jer, overflow_check, err, only_on_PE) result(EFP_sum)
  real, dimension(:,:),     intent(in)  :: array   !< The array to be summed
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
  integer,        optional, intent(out) :: err  !< If present, return an error code instead of
                                                !! triggering any fatal errors directly from
                                                !! this routine.
  logical,        optional, intent(in)  :: only_on_PE !< If present and true, do not do the sum
                                                !! across processors, only reporting the local sum
  type(EFP_type)                        :: EFP_sum  !< The result in extended fixed point format

  !   This subroutine uses a conversion to an integer representation
  ! of real numbers to give order-invariant sums that will reproduce
  ! across PE count.  This idea comes from R. Hallberg and A. Adcroft.

  integer(kind=8), dimension(ni)  :: ints_sum
  integer(kind=8) :: ival, prec_error
  real    :: rs
  real    :: max_mag_term
  logical :: over_check, do_sum_across_PEs
  character(len=256) :: mesg
  integer :: i, j, n, is, ie, js, je, sgn

  if (num_PEs() > max_count_prec) call MOM_error(FATAL, &
    "reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/num_PEs.")

  prec_error = (2_8**62 + (2_8**62 - 1)) / num_PEs()

  is = 1 ; ie = size(array,1) ; js = 1 ; je = size(array,2 )
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

  overflow_error = .false. ; NaN_error = .false. ; max_mag_term = 0.0
  ints_sum(:) = 0
  if (over_check) then
    if ((je+1-js)*(ie+1-is) < max_count_prec) then
      do j=js,je ; do i=is,ie
        call increment_ints_faster(ints_sum, array(i,j), max_mag_term)
      enddo ; enddo
      call carry_overflow(ints_sum, prec_error)
    elseif ((ie+1-is) < max_count_prec) then
      do j=js,je
        do i=is,ie
          call increment_ints_faster(ints_sum, array(i,j), max_mag_term)
        enddo
        call carry_overflow(ints_sum, prec_error)
      enddo
    else
      do j=js,je ; do i=is,ie
        call increment_ints(ints_sum, real_to_ints(array(i,j), prec_error), &
                            prec_error)
      enddo ; enddo
    endif
  else
    do j=js,je ; do i=is,ie
      sgn = 1 ; if (array(i,j)<0.0) sgn = -1
      rs = abs(array(i,j))
      do n=1,ni
        ival = int(rs*I_pr(n), 8)
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
    if (err > 0) then ; do n=1,ni ; ints_sum(n) = 0 ; enddo ; endif
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

  if (do_sum_across_PEs) call sum_across_PEs(ints_sum, ni)

  call regularize_ints(ints_sum)

  EFP_sum%v(:) = ints_sum(:)

end function reproducing_EFP_sum_2d

!> This subroutine uses a conversion to an integer representation of real numbers to give an
!! order-invariant sum of distributed 2-D arrays that reproduces across domain decomposition.
!! This technique is described in Hallberg & Adcroft, 2014, Parallel Computing,
!! doi:10.1016/j.parco.2014.04.007.
function reproducing_sum_2d(array, isr, ier, jsr, jer, EFP_sum, reproducing, &
                            overflow_check, err, only_on_PE) result(sum)
  real, dimension(:,:),     intent(in)  :: array   !< The array to be summed
  integer,        optional, intent(in)  :: isr     !< The starting i-index of the sum, noting
                                                   !! that the array indices starts at 1
  integer,        optional, intent(in)  :: ier     !< The ending i-index of the sum, noting
                                                   !! that the array indices starts at 1
  integer,        optional, intent(in)  :: jsr     !< The starting j-index of the sum, noting
                                                   !! that the array indices starts at 1
  integer,        optional, intent(in)  :: jer     !< The ending j-index of the sum, noting
                                                   !! that the array indices starts at 1
  type(EFP_type), optional, intent(out) :: EFP_sum  !< The result in extended fixed point format
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
  real                                  :: sum  !< Result

  !   This subroutine uses a conversion to an integer representation
  ! of real numbers to give order-invariant sums that will reproduce
  ! across PE count.  This idea comes from R. Hallberg and A. Adcroft.

  integer(kind=8), dimension(ni)  :: ints_sum
  integer(kind=8) :: prec_error
  real    :: rsum(1), rs
  logical :: repro, do_sum_across_PEs
  character(len=256) :: mesg
  type(EFP_type) :: EFP_val ! An extended fixed point version of the sum
  integer :: i, j, n, is, ie, js, je

  if (num_PEs() > max_count_prec) call MOM_error(FATAL, &
    "reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/num_PEs.")

  prec_error = (2_8**62 + (2_8**62 - 1)) / num_PEs()

  is = 1 ; ie = size(array,1) ; js = 1 ; je = size(array,2 )
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

  if (repro) then
    EFP_val = reproducing_EFP_sum_2d(array, isr, ier, jsr, jer, overflow_check, err, only_on_PE)
    sum = ints_to_real(EFP_val%v)
    if (present(EFP_sum)) EFP_sum = EFP_val
    if (debug) ints_sum(:) = EFP_sum%v(:)
  else
    rsum(1) = 0.0
    do j=js,je ; do i=is,ie
      rsum(1) = rsum(1) + array(i,j)
    enddo ; enddo
    if (do_sum_across_PEs) call sum_across_PEs(rsum,1)
    sum = rsum(1)

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
    write(mesg,'("2d RS: ", ES24.16, 6 Z17.16)') sum, ints_sum(1:ni)
    call MOM_mesg(mesg, 3)
  endif

end function reproducing_sum_2d

!> This subroutine uses a conversion to an integer representation of real numbers to give an
!! order-invariant sum of distributed 3-D arrays that reproduces across domain decomposition.
!! This technique is described in Hallberg & Adcroft, 2014, Parallel Computing,
!! doi:10.1016/j.parco.2014.04.007.
function reproducing_sum_3d(array, isr, ier, jsr, jer, sums, EFP_sum, EFP_lay_sums, err, only_on_PE) &
                            result(sum)
  real, dimension(:,:,:),       intent(in)  :: array   !< The array to be summed
  integer,            optional, intent(in)  :: isr     !< The starting i-index of the sum, noting
                                                       !! that the array indices starts at 1
  integer,            optional, intent(in)  :: ier     !< The ending i-index of the sum, noting
                                                       !! that the array indices starts at 1
  integer,            optional, intent(in)  :: jsr     !< The starting j-index of the sum, noting
                                                       !! that the array indices starts at 1
  integer,            optional, intent(in)  :: jer     !< The ending j-index of the sum, noting
                                                       !! that the array indices starts at 1
  real, dimension(:), optional, intent(out) :: sums    !< The sums by vertical layer
  type(EFP_type),     optional, intent(out) :: EFP_sum !< The result in extended fixed point format
  type(EFP_type), dimension(:), &
                      optional, intent(out) :: EFP_lay_sums !< The sums by vertical layer in EFP format
  integer,            optional, intent(out) :: err  !< If present, return an error code instead of
                                                    !! triggering any fatal errors directly from
                                                    !! this routine.
  logical,            optional, intent(in)  :: only_on_PE !< If present and true, do not do the sum
                                                    !! across processors, only reporting the local sum
  real                                      :: sum  !< Result

  !   This subroutine uses a conversion to an integer representation
  ! of real numbers to give order-invariant sums that will reproduce
  ! across PE count.  This idea comes from R. Hallberg and A. Adcroft.

  real    :: val, max_mag_term
  integer(kind=8), dimension(ni)  :: ints_sum
  integer(kind=8), dimension(ni,size(array,3))  :: ints_sums
  integer(kind=8) :: prec_error
  character(len=256) :: mesg
  logical :: do_sum_across_PEs
  integer :: i, j, k, is, ie, js, je, ke, isz, jsz, n

  if (num_PEs() > max_count_prec) call MOM_error(FATAL, &
    "reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/num_PEs.")

  prec_error = (2_8**62 + (2_8**62 - 1)) / num_PEs()
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
  jsz = je+1-js; isz = ie+1-is

  do_sum_across_PEs = .true. ; if (present(only_on_PE)) do_sum_across_PEs = .not.only_on_PE

  if (present(sums) .or. present(EFP_lay_sums)) then
    if (present(sums)) then ; if (size(sums) < ke) then
      call MOM_error(FATAL, "Sums is smaller than the vertical extent of array in reproducing_sum(_3d).")
    endif ; endif
    if (present(EFP_lay_sums)) then ; if (size(EFP_lay_sums) < ke) then
      call MOM_error(FATAL, "Sums is smaller than the vertical extent of array in reproducing_sum(_3d).")
    endif ; endif
    ints_sums(:,:) = 0
    overflow_error = .false. ; NaN_error = .false. ; max_mag_term = 0.0
    if (jsz*isz < max_count_prec) then
      do k=1,ke
        do j=js,je ; do i=is,ie
          call increment_ints_faster(ints_sums(:,k), array(i,j,k), max_mag_term)
        enddo ; enddo
        call carry_overflow(ints_sums(:,k), prec_error)
      enddo
    elseif (isz < max_count_prec) then
      do k=1,ke ; do j=js,je
        do i=is,ie
          call increment_ints_faster(ints_sums(:,k), array(i,j,k), max_mag_term)
        enddo
        call carry_overflow(ints_sums(:,k), prec_error)
      enddo ; enddo
    else
      do k=1,ke ; do j=js,je ; do i=is,ie
        call increment_ints(ints_sums(:,k), &
                            real_to_ints(array(i,j,k), prec_error), prec_error)
      enddo ; enddo ; enddo
    endif
    if (present(err)) then
      err = 0
      if (abs(max_mag_term) >= prec_error*pr(1)) err = err+1
      if (overflow_error) err = err+2
      if (NaN_error) err = err+2
      if (err > 0) then ; do k=1,ke ; do n=1,ni ; ints_sums(n,k) = 0 ; enddo ; enddo ; endif
    else
      if (NaN_error) call MOM_error(FATAL, "NaN in input field of reproducing_sum(_3d).")
      if (abs(max_mag_term) >= prec_error*pr(1)) then
        write(mesg, '(ES13.5)') max_mag_term
        call MOM_error(FATAL,"Overflow in reproducing_sum(_3d) conversion of "//trim(mesg))
      endif
      if (overflow_error) call MOM_error(FATAL, "Overflow in reproducing_sum(_3d).")
    endif

    if (do_sum_across_PEs) call sum_across_PEs(ints_sums(:,1:ke), ni*ke)

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
      do n=1,ni ; ints_sum(n) = 0 ; enddo
      do k=1,ke ; do n=1,ni ; ints_sum(n) = ints_sum(n) + ints_sums(n,k) ; enddo ; enddo
      write(mesg,'("3D RS: ", ES24.16, 6 Z17.16)') sum, ints_sum(1:ni)
      call MOM_mesg(mesg, 3)
    endif
  else
    ints_sum(:) = 0
    overflow_error = .false. ; NaN_error = .false. ; max_mag_term = 0.0
    if (jsz*isz < max_count_prec) then
      do k=1,ke
        do j=js,je ; do i=is,ie
          call increment_ints_faster(ints_sum, array(i,j,k), max_mag_term)
        enddo ; enddo
        call carry_overflow(ints_sum, prec_error)
      enddo
    elseif (isz < max_count_prec) then
      do k=1,ke ; do j=js,je
        do i=is,ie
          call increment_ints_faster(ints_sum, array(i,j,k), max_mag_term)
        enddo
        call carry_overflow(ints_sum, prec_error)
      enddo ; enddo
    else
      do k=1,ke ; do j=js,je ; do i=is,ie
        call increment_ints(ints_sum, real_to_ints(array(i,j,k), prec_error), &
                            prec_error)
      enddo ; enddo ; enddo
    endif
    if (present(err)) then
      err = 0
      if (abs(max_mag_term) >= prec_error*pr(1)) err = err+1
      if (overflow_error) err = err+2
      if (NaN_error) err = err+2
      if (err > 0) then ; do n=1,ni ; ints_sum(n) = 0 ; enddo ; endif
    else
      if (NaN_error) call MOM_error(FATAL, "NaN in input field of reproducing_sum(_3d).")
      if (abs(max_mag_term) >= prec_error*pr(1)) then
        write(mesg, '(ES13.5)') max_mag_term
        call MOM_error(FATAL,"Overflow in reproducing_sum(_3d) conversion of "//trim(mesg))
      endif
      if (overflow_error) call MOM_error(FATAL, "Overflow in reproducing_sum(_3d).")
    endif

    if (do_sum_across_PEs) call sum_across_PEs(ints_sum, ni)

    call regularize_ints(ints_sum)
    sum = ints_to_real(ints_sum)

    if (present(EFP_sum)) EFP_sum%v(:) = ints_sum(:)

    if (debug) then
      write(mesg,'("3d RS: ", ES24.16, 6 Z17.16)') sum, ints_sum(1:ni)
      call MOM_mesg(mesg, 3)
    endif
  endif

end function reproducing_sum_3d

!> Convert a real number into the array of integers constitute its extended-fixed-point representation
function real_to_ints(r, prec_error, overflow) result(ints)
  real,                      intent(in) :: r  !< The real number being converted
  integer(kind=8), optional, intent(in) :: prec_error  !< The PE-count dependent precision of the
                                              !! integers that is safe from overflows during global
                                              !! sums.  This will be larger than the compile-time
                                              !! precision parameter, and is used to detect overflows.
  logical,         optional, intent(inout) :: overflow !< Returns true if the conversion is being
                                              !! done on a value that is too large to be represented
  integer(kind=8), dimension(ni)  :: ints
  !   This subroutine converts a real number to an equivalent representation
  ! using several long integers.

  real :: rs
  character(len=80) :: mesg
  integer(kind=8) :: ival, prec_err
  integer :: sgn, i

  prec_err = prec ; if (present(prec_error)) prec_err = prec_error
  ints(:) = 0_8
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

  do i=1,ni
    ival = int(rs*I_pr(i), 8)
    rs = rs - ival*pr(i)
    ints(i) = sgn*ival
  enddo

end function real_to_ints

!> Convert the array of integers that constitute an extended-fixed-point
!! representation into a real number
function ints_to_real(ints) result(r)
  integer(kind=8), dimension(ni), intent(in) :: ints !< The array of EFP integers
  real :: r
  ! This subroutine reverses the conversion in real_to_ints.

  integer :: i

  r = 0.0
  do i=1,ni ; r = r + pr(i)*ints(i) ; enddo
end function ints_to_real

!> Increment an array of integers that constitutes an extended-fixed-point
!! representation with a another EFP number
subroutine increment_ints(int_sum, int2, prec_error)
  integer(kind=8), dimension(ni), intent(inout) :: int_sum !< The array of EFP integers being incremented
  integer(kind=8), dimension(ni), intent(in)    :: int2    !< The array of EFP integers being added
  integer(kind=8), optional,      intent(in)    :: prec_error !< The PE-count dependent precision of the
                                              !! integers that is safe from overflows during global
                                              !! sums.  This will be larger than the compile-time
                                              !! precision parameter, and is used to detect overflows.

  ! This subroutine increments a number with another, both using the integer
  ! representation in real_to_ints.
  integer :: i

  do i=ni,2,-1
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

!> Increment an EFP number with a real number without doing any carrying of
!! of overflows and using only minimal error checking.
subroutine increment_ints_faster(int_sum, r, max_mag_term)
  integer(kind=8), dimension(ni), intent(inout) :: int_sum  !< The array of EFP integers being incremented
  real,                           intent(in)    :: r        !< The real number being added.
  real,                           intent(inout) :: max_mag_term !< A running maximum magnitude of the r's.

  ! This subroutine increments a number with another, both using the integer
  ! representation in real_to_ints, but without doing any carrying of overflow.
  ! The entire operation is embedded in a single call for greater speed.
  real :: rs
  integer(kind=8) :: ival
  integer :: sgn, i

  if ((r >= 1e30) .eqv. (r < 1e30)) then ; NaN_error = .true. ; return ; endif
  sgn = 1 ; if (r<0.0) sgn = -1
  rs = abs(r)
  if (rs > abs(max_mag_term)) max_mag_term = r

  ! Abort if the number has no EFP representation
  if (rs > max_efp_float) then
    overflow_error = .true.
    return
  endif

  do i=1,ni
    ival = int(rs*I_pr(i), 8)
    rs = rs - ival*pr(i)
    int_sum(i) = int_sum(i) + sgn*ival
  enddo

end subroutine increment_ints_faster

!> This subroutine handles carrying of the overflow.
subroutine carry_overflow(int_sum, prec_error)
  integer(kind=8), dimension(ni), intent(inout) :: int_sum  !< The array of EFP integers being
                                              !! modified by carries, but without changing value.
  integer(kind=8),                intent(in)    :: prec_error  !< The PE-count dependent precision of the
                                              !! integers that is safe from overflows during global
                                              !! sums.  This will be larger than the compile-time
                                              !! precision parameter, and is used to detect overflows.

  ! This subroutine handles carrying of the overflow.
  integer :: i, num_carry

  do i=ni,2,-1 ; if (abs(int_sum(i)) >= prec) then
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
  integer(kind=8), dimension(ni), &
    intent(inout) :: int_sum !< The array of integers being modified to take a
                             !! regular form with all integers of the same sign,
                             !! but without changing value.

  ! This subroutine carries the overflow, and then makes sure that
  ! all integers are of the same sign as the overall value.
  logical :: positive
  integer :: i, num_carry

  do i=ni,2,-1 ; if (abs(int_sum(i)) >= prec) then
    num_carry = int(int_sum(i) * I_prec)
    int_sum(i) = int_sum(i) - num_carry*prec
    int_sum(i-1) = int_sum(i-1) + num_carry
  endif ; enddo

  ! Determine the sign of the final number.
  positive = .true.
  do i=1,ni
    if (abs(int_sum(i)) > 0) then
      if (int_sum(i) < 0) positive = .false.
      exit
    endif
  enddo

  if (positive) then
    do i=ni,2,-1 ; if (int_sum(i) < 0) then
      int_sum(i) = int_sum(i) + prec
      int_sum(i-1) = int_sum(i-1) - 1
    endif ; enddo
  else
    do i=ni,2,-1 ; if (int_sum(i) > 0) then
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

  do i=1,ni ; EFP_minus%v(i) = -1*EFP2%v(i) ; enddo

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

  do i=1,ni ; EFP1%v(i) = EFP2%v(i) ; enddo
end subroutine EFP_assign

!> Return the real number that an extended-fixed-point number corresponds with
function EFP_to_real(EFP1)
  type(EFP_type), intent(inout) :: EFP1 !< The extended fixed point number being converted
  real :: EFP_to_real

  call regularize_ints(EFP1%v)
  EFP_to_real = ints_to_real(EFP1%v)
end function EFP_to_real

!> Take the difference between two extended-fixed-point numbers (EFP1 - EFP2)
!! and return the result as a real number
function EFP_real_diff(EFP1, EFP2)
  type(EFP_type), intent(in) :: EFP1  !< The first extended fixed point number
  type(EFP_type), intent(in) :: EFP2  !< The extended fixed point number being
                        !! subtracted from the first extended fixed point number
  real :: EFP_real_diff !< The real result

  type(EFP_type)             :: EFP_diff

  EFP_diff = EFP1 - EFP2
  EFP_real_diff = EFP_to_real(EFP_diff)

end function EFP_real_diff

!> Return the extended-fixed-point number that a real number corresponds with
function real_to_EFP(val, overflow)
  real,              intent(in)    :: val !< The real number being converted
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

  integer(kind=8), dimension(ni,nval) :: ints
  integer(kind=8) :: prec_error
  logical :: error_found
  character(len=256) :: mesg
  integer :: i, n

  if (num_PEs() > max_count_prec) call MOM_error(FATAL, &
    "reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/num_PEs.")

  prec_error = (2_8**62 + (2_8**62 - 1)) / num_PEs()
  ! overflow_error is an overflow error flag for the whole module.
  overflow_error = .false. ; error_found = .false.

  do i=1,nval ; do n=1,ni ; ints(n,i) = EFPs(i)%v(n) ; enddo ; enddo

  call sum_across_PEs(ints(:,:), ni*nval)

  if (present(errors)) errors(:) = .false.
  do i=1,nval
    overflow_error = .false.
    call carry_overflow(ints(:,i), prec_error)
    do n=1,ni ; EFPs(i)%v(n) = ints(n,i) ; enddo
    if (present(errors)) errors(i) = overflow_error
    if (overflow_error) then
      write (mesg,'("EFP_list_sum_across_PEs error at ",i6," val was ",ES12.6, ", prec_error = ",ES12.6)') &
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

  integer(kind=8), dimension(ni) :: ints
  integer(kind=8) :: prec_error
  logical :: error_found
  character(len=256) :: mesg
  integer :: n

  if (num_PEs() > max_count_prec) call MOM_error(FATAL, &
    "reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/num_PEs.")

  prec_error = (2_8**62 + (2_8**62 - 1)) / num_PEs()
  ! overflow_error is an overflow error flag for the whole module.
  overflow_error = .false. ; error_found = .false.

  do n=1,ni ; ints(n) = EFP%v(n) ; enddo

  call sum_across_PEs(ints(:), ni)

  if (present(error)) error = .false.

  overflow_error = .false.
  call carry_overflow(ints(:), prec_error)
  do n=1,ni ; EFP%v(n) = ints(n) ; enddo
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


!> This subroutine carries out all of the calls required to close out the infrastructure cleanly.
!! This should only be called in ocean-only runs, as the coupler takes care of this in coupled runs.
subroutine MOM_infra_end
  call print_memuse_stats( 'Memory HiWaterMark', always=.TRUE. )
  call fms_end
end subroutine MOM_infra_end

end module MOM_coms
