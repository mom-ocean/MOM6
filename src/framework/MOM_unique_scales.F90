!> This module provides tools that can be used to check the uniqueness of the dimensional
!! scaling factors used by the MOM6 ocean model or other models
module MOM_unique_scales

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, assert, MOM_get_verbosity

implicit none ; private

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

public check_scaling_uniqueness, scales_to_powers

contains

!> This subroutine does a checks whether the provided dimensional scaling factors give a unique
!! overall scaling for each of the combinations of units in description, and suggests a better
!! combination if it is not unique.  However, this subroutine does nothing if the verbosity level
!! for this run is below 3.
subroutine check_scaling_uniqueness(component, descs, weights, key, scales, max_powers)
  character(len=*),  intent(in) :: component  !< The name of the component (e.g., MOM6) to use in messages
  character(len=*),  intent(in) :: descs(:)   !< The descriptions for each combination of units
  integer,           intent(in) :: weights(:) !< A list of the weights for each described combination
  character(len=*),  intent(in) :: key(:)     !< The key for the unit scaling
  integer,           intent(in) :: scales(:)  !< The powers of 2 that give the scaling for each unit in key
  integer, optional, intent(in) :: max_powers !< The maximum range of powers of 2 to search for
                                              !! suggestions of better scaling factors, or 0 to avoid
                                              !! suggesting improved factors.

  ! Local variables
  integer, dimension(size(key)) :: next_scales, prev_scales, better_scales
  character(len=512) :: mesg
  character(len=64) :: msg_frag
  integer, dimension(size(key), size(weights)) :: list
  integer :: verbosity
  logical :: same_key
  integer :: orig_cost, test_cost, better_cost, prev_cost ! Various squared-weight mismatch costs.
  integer :: better_dp ! The absolute change in powers with the better estimate.
  integer :: ndims, ns, m, n, i, p, itt, max_itt, max_pow

  call assert((size(scales) == size(key)), "check_scaling_factors: Mismatched scales and key sizes.")
  call assert((size(descs) == size(weights)), "check_scaling_factors: Mismatched descs and weights.")

  verbosity = MOM_get_verbosity()
  ! Skip the rest of this routine if it would not write anything out.
  if (verbosity < 3) return

  ndims = size(key)
  ns = size(weights)
  max_pow = 0 ; if (present(max_powers)) max_pow = max_powers

  list(:,:) = 0
  do n=1,ns
    call encode_dim_powers(descs(n), key, list(:,n))
  enddo

  if (verbosity >= 7) then
    write(mesg, '(I8)') ns
    call MOM_mesg(trim(component)//": Extracted "//trim(adjustl(mesg))//" unit combinations from the list.")
    mesg = "Dim Key:  ["
    do i=1,ndims ; mesg = trim(mesg)//"  "//trim(key(i)) ; enddo
    mesg = trim(mesg)//"]:"
    call MOM_mesg(mesg)
    do n=1,ns
      call MOM_mesg(trim(component)//": Extracted ["//trim(int_array_msg(list(:,n)))//"] from "//trim(descs(n)))
    enddo

    do n=1,ns ; do m=1,n-1
      same_key = .true.
      do i=1,ndims ; if (list(i,n) /= list(i,m)) same_key = .false. ; enddo
      if (same_key) then
        call MOM_mesg(trim(component)//": The same powers occur for "//&
                      trim(descs(n))//" and "//trim(descs(m))//"." )
      endif
    enddo ; enddo
  endif

  orig_cost = non_unique_scales(scales, list, descs, weights, silent=(verbosity<4))

  max_itt = 3*ndims  ! Do up to 3 iterations for each rescalable dimension.
  if (orig_cost /= 0) then
    call MOM_mesg(trim(component)//": The dimensional scaling factors are not unique.")
    prev_cost = orig_cost
    prev_scales(:) = scales(:)
    do itt=1,max_itt
      ! Iterate to find a better solution.
      better_scales(:) = prev_scales(:)
      better_cost = prev_cost
      better_dp = 0
      do i=1,ndims
        if (scales(i) == 0) cycle  ! DO not optimize unscaled dimensions.
        next_scales(:) = prev_scales(:)
        do p=-max_pow,max_pow
          if ((p==0) .or. (p==prev_scales(i))) cycle
          next_scales(i) = p
          test_cost = non_unique_scales(next_scales, list, descs, weights, silent=.true.)
          if ((test_cost < better_cost) .or. &
              ((test_cost == better_cost) .and. (abs(p-prev_scales(i)) < better_dp))) then
            ! This is a better scaling or has the same weighted mismatches but smaller
            ! changes in rescaling factors, so it could be the next guess.
            better_scales(:) = next_scales(:)
            better_cost = test_cost
            better_dp = abs(p - prev_scales(i))
          endif
        enddo
      enddo
      if (better_cost < prev_cost) then
        ! Store the new best guess and try again.
        prev_scales(:) = better_scales(:)
        prev_cost = better_cost
      else ! No further optimization is possible.
        exit
      endif
      if (better_cost == 0) exit
      if (verbosity >= 7) then
        write(mesg, '("Iteration ",I2," scaling cost reduced from ",I8," with original scales to ", I8)') &
                    itt, orig_cost, better_cost
        call MOM_mesg(trim(component)//": "//trim(mesg)//" with revised scaling factors.")
      endif
    enddo
    if (prev_cost < orig_cost) then
      test_cost = non_unique_scales(prev_scales, list, descs, weights, silent=(verbosity<4))
      mesg = trim(component)//": Suggested improved scales: "
      do i=1,ndims ; if ((prev_scales(i) /= scales(i)) .and. (scales(i) /= 0)) then
        write(msg_frag, '(I3)') prev_scales(i)
        mesg = trim(mesg)//" "//trim(key(i))//"_RESCALE_POWER = "//trim(adjustl(msg_frag))
      endif ; enddo
      call MOM_mesg(mesg)

      write(mesg, '(I8)') orig_cost
      write(msg_frag, '(I8)') test_cost
      mesg = trim(component)//": Scaling overlaps reduced from "//trim(adjustl(mesg))//&
             " with original scales to "//trim(adjustl(msg_frag))//" with suggested scales."
      call MOM_mesg(mesg)
    endif

  endif

end subroutine check_scaling_uniqueness

!> Convert a unit scaling descriptor into an array of the dimensions of powers given in the key
subroutine encode_dim_powers(scaling, key, dim_powers)

  character(len=*),               intent(in)  :: scaling   !< The unit description that will be converted
  character(len=*), dimension(:), intent(in)  :: key(:)    !< The key for the unit scaling
  integer, dimension(size(key)),  intent(out) :: dim_powers !< The dimensions in scaling of each
                                                           !! element of they key.

  ! Local variables
  character(len=:), allocatable :: actstr ! The full active remaining string to be parsed.
  character(len=:), allocatable :: fragment ! The space-delimited fragment being parsed.
  character(len=:), allocatable :: dimnm  ! The probable dimension name
  character(len=11) :: numbers ! The list of characters that could make up the exponent.
  ! character(len=128) :: mesg
  integer :: istart, iend, ieq, ief, ipow  ! Positions in strings.
  integer :: dp   ! The power for this dimension.
  integer :: ndim ! The number of dimensional scaling factors to consider.
  integer :: n

  dim_powers(:) = 0

  iend = index(scaling, "~>") - 1
  if (iend < 1) return

  ! Parse the key.
  ndim = size(key)
  numbers = "-0123456789"

  ! Strip away any leading square brace.
  istart = index(scaling(:iend), "[") + 1
  ! If there is an "=" in the string, start after this.
  ieq = index(scaling(istart:iend), "=", back=.true.)
  if (ieq > 0) istart = istart + ieq

  ! Set up the active string to work on.
  actstr = trim(adjustl(scaling(istart:iend)))
  do  ! Loop over each of the elements in the unit scaling descriptor.
    if (len_trim(actstr) == 0) exit
    ief = index(actstr, " ") - 1
    if (ief <= 0) ief = len_trim(actstr)
    fragment = actstr(:ief)
    ipow = scan(fragment, "-")
    if (ipow == 0) ipow = scan(fragment, numbers)

    if (ipow == 0) then ! There is no exponent
      dimnm = fragment
      dp = 1
      ! call MOM_mesg("Parsing powerless fragment "//trim(fragment)//" from "//trim(scaling))
    else
      if (verify(fragment(ipow:), numbers) == 0) then
        read(fragment(ipow:),*) dp
        dimnm = fragment(:ipow-1)
        ! write(mesg, '(I3)') dp
        ! call MOM_mesg("Parsed fragment "//trim(fragment)//" from "//trim(scaling)//&
        !               " as "//trim(dimnm)//trim(adjustl(mesg)))
      else
        dimnm = fragment
        dp = 1
        ! call MOM_mesg("Unparsed fragment "//trim(fragment)//" from "//trim(scaling))
      endif
    endif

    do n=1,ndim
      if (trim(dimnm) == trim(key(n))) then
        dim_powers(n) = dim_powers(n) + dp
        exit
      endif
    enddo

    ! Remove the leading fragment that has been parsed from actstr
    actstr = trim(adjustl(actstr(ief+1:)))
  enddo

end subroutine encode_dim_powers

!> Find the integer power of two that describe each of the scaling factors, or return 0 for
!! scaling factors that are not exceptionally close to an integer power of 2.
subroutine scales_to_powers(scale, pow2)
  real,    intent(in)  :: scale(:)  !< The scaling factor for each dimension
  integer, intent(out) :: pow2(:)   !< The exact powers of 2 for each scale, or 0 for non-exact powers of 2.

  real :: log2_sc        ! The log base 2 of an element of scale
  integer :: n, ndim

  ndim = size(scale)

  ! Find the integer power of two for the scaling factors, but skip the analysis of any factors
  ! that are not close enough to being integer powers of 2.
  do n=1,ndim
    if (abs(scale(n)) > 0.0) then
      log2_sc = log(abs(scale(n))) / log(2.0)
    else
      log2_sc = 0.0
    endif
    if (abs(log2_sc - nint(log2_sc)) < 1.0e-6) then
      ! This is close to an integer power of two.
      pow2(n) = nint(log2_sc)
    else
      ! This is not being scaled by an integer power of 2, so return 0.
      pow2(n) = 0
    endif
  enddo

end subroutine scales_to_powers

!> Determine from the list of scaling factors and the unit combinations that are in use whether
!! all these combinations scale uniquely.
integer function non_unique_scales(scales, list, descs, weights, silent)
  integer,           intent(in) :: scales(:)  !< The power of 2 that gives the scaling factor for each dimension
  integer,           intent(in) :: list(:,:)  !< A list of the integers for each scaling
  character(len=*),  intent(in) :: descs(:)   !< The unit descriptions that have been converted
  integer,           intent(in) :: weights(:) !< A list of the weights for each scaling
  logical, optional, intent(in) :: silent     !< If present and true, do not write any output.

  ! Local variables
  integer, dimension(size(weights)) :: res_pow  ! The net rescaling power for each combination.
  integer, dimension(size(weights)) :: wt_merge ! The merged weights of scaling factors with common powers
                                                ! for the dimensions being tested.
  logical :: same_key, same_scales, verbose
  character(len=256) :: mesg
  integer :: nonzero_count  ! The number of non-zero scaling factors
  integer :: ndim           ! The number of dimensional scaling factors to work with
  integer :: i, n, m, ns

  verbose = .true. ; if (present(silent)) verbose = .not.silent

  ndim = size(scales)
  ns = size(descs)
  call assert((size(scales) == size(list, 1)), "non_unique_scales: Mismatched scales and list sizes.")
  call assert((size(descs) == size(list, 2)), "non_unique_scales: Mismatched descs and list sizes.")
  call assert((size(descs) == size(weights)), "non_unique_scales: Mismatched descs and weights.")

  ! Return .true. if all scaling powers are 0, or there is only one scaling factor in use.
  nonzero_count = 0 ; do n=1,ndim ; if (scales(n) /= 0) nonzero_count = nonzero_count + 1 ; enddo
  if (nonzero_count <= 1) return

  ! Figure out which unit combinations are unique for the set of dimensions and scaling factors
  ! that are being tested, and combine the weights for scaling factors.
  wt_merge(:) = weights(:)
  do n=1,ns ; do m=1,n-1
    same_key = .true.
    same_scales = .true.
    do i=1,ndim
      if (list(i,n) /= list(i,m)) same_key = .false.
      if ((scales(i) /= 0) .and. (list(i,n) /= list(i,m))) same_scales = .false.
    enddo
    if (same_key .or. same_scales) then
      if (wt_merge(n) > wt_merge(m)) then
        wt_merge(n) = wt_merge(n) + wt_merge(m)
        wt_merge(m) = 0
      else
        wt_merge(m) = wt_merge(m) + wt_merge(n)
        wt_merge(n) = 0
      endif
    endif
    if (wt_merge(n) == 0) exit ! Go to the next value of n.
  enddo ; enddo

  do n=1,ns
    res_pow(n) = 0
    do i=1,ndim
      res_pow(n) = res_pow(n) + scales(i) * list(i,n)
    enddo
  enddo

  ! Determine the weighted cost of non-unique scaling factors.
  non_unique_scales = 0
  do n=1,ns ; if (wt_merge(n) > 0) then ; do m=1,n-1 ; if (wt_merge(m) > 0) then
    if (res_pow(n) == res_pow(m)) then
      ! Use the product of the weights as the cost, as this should be vaguely proportional to
      ! the likelihood that these factors would be combined in an expression.
      non_unique_scales = min(non_unique_scales + wt_merge(n) * wt_merge(m), 99999999)
      if (verbose) then
        write(mesg, '(I8)') res_pow(n)
        call MOM_mesg("The factors "//trim(descs(n))//" and "//trim(descs(m))//" both scale to "//&
                      trim(adjustl(mesg))//" for the given powers.")

        ! call MOM_mesg("Powers ["//trim(int_array_msg(list(:,n)))//"] and ["//&
        !                    trim(int_array_msg(list(:,m)))//"] with rescaling by ["//&
        !                    trim(int_array_msg(scales))//"]")
      endif
    endif
  endif ; enddo ; endif ; enddo

end function non_unique_scales

!> Return a string the elements of an array of integers
function int_array_msg(array)
  integer,  intent(in) :: array(:)  !< The array whose values are to be written.
  character(len=16*size(array)) :: int_array_msg

  character(len=12) :: msg_frag
  integer :: i, ni
  ni = size(array)

  int_array_msg = ""
  if (ni < 1) return

  do i=1,ni
    write(msg_frag, '(I8)') array(i)
    msg_frag = adjustl(msg_frag)
    if (i == 1) then
      int_array_msg = trim(msg_frag)
    else
      int_array_msg = trim(int_array_msg)//" "//trim(msg_frag)
    endif
  enddo
end function int_array_msg

end module MOM_unique_scales
