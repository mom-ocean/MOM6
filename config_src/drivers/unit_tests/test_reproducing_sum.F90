program test_reproducing_sum

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms, only : PE_here, root_PE, num_PEs, reproducing_sum
use MOM_coms, only : sum_across_PEs, max_across_PEs, max_count_prec
use MOM_domains, only : MOM_domain_type, create_MOM_domain, MOM_infra_init, MOM_infra_end
use MOM_domains, only : MOM_define_layout
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, MOM_set_verbosity
use MOM_hor_index, only : hor_index_type, hor_index_init

  implicit none

  type(MOM_domain_type), pointer :: Domain => NULL() ! Ocean model domain
  type(hor_index_type) :: HI ! A hor_index_type for array extents
  real, allocatable :: array(:,:) ! An array with values to sum over [A]
  real :: tot_R, tot_std, tot_fastR ! Sums via different methods [A]
  real :: error_bound, likely_error ! Errors via different methods [A]
  character(len=200) :: mesg ! String for messages
  integer :: n_repeat ! Number of times to repeat the sum call
  integer :: n ! Loop counter
  integer :: io_unit ! i/o unit for creating input.nml (sigh)
  integer :: n_global(2) ! Global i-, j- dimensions of domain (h-points)
  integer :: layout(2)   ! PE count in i-, j- directions
  integer :: PEs_used    ! Number of PEs available to executable
  logical :: tests_failed ! True if a fail is encountered
  integer :: i, j, ig, jg ! Spatial indices

  ! FMS requires the file "input.nml" to exist ...
  open(newunit=io_unit, file="input.nml", status="replace", action="write")
  close(io_unit) ! ... but an empty input.nml is sufficient

  call MOM_infra_init()

  n_repeat = 100

  ! Optionally use command-line to change size of the problem
  ! Usage:  ./executable [tile-size] [number-of-calls]
  n = command_argument_count()
  if (n==2) then
    call get_command_argument(1, mesg)
    read(mesg,*) n_global(1)
    n_global(2) = n_global(1)
    call get_command_argument(2, mesg)
    read(mesg,*) n_repeat
  elseif (n==1) then
    call get_command_argument(1, mesg)
    read(mesg,*) n_global(1)
    n_global(2) = n_global(1)
  else
    n_global = (/200, 300/) ! Fallback value if no argument provided
  endif

  tests_failed = .false.
  call MOM_set_verbosity(2)

  ! Setup distributed domain
  PEs_used = num_PEs()
  call MOM_define_layout(n_global, PEs_used, layout)
  call create_MOM_domain(Domain, n_global, (/2,2/), (/.false.,.false./), .false., layout)
  call hor_index_init(Domain, HI)

  allocate( array(HI%isd:HI%ied,HI%jsd:HI%jed), source=0. )

  ! Set up an array of values to sum
  call generate_array_of_values(array, HI, n_global)

  ! This estimates the maximum possible accumulated round off error, and likely error
  ! from a random walk of round off errors
  error_bound = 0.
  tot_std = 0.
  do j = HI%jsc, HI%jec ; do i = HI%isc, HI%iec
    ! Actual round off error for adding tot_std + array(i,j)
    error_bound = error_bound + max( abs(tot_std), abs(array(i,j)) ) * epsilon(error_bound)
    tot_std = tot_std + array(i,j)
  enddo ; enddo
  call sum_across_PEs( error_bound )
  call sum_across_PEs( tot_std )
  N = n_global(1) * n_global(2)
  likely_error = tot_std * epsilon(tot_std) * sqrt( real( N ) )
  if (likely_error > error_bound) call MOM_error(FATAL, 'Something went wrong in error estimate!')

  tot_std = reproducing_sum(array, HI%isc, HI%iec, HI%jsc, HI%jec, reproducing=.false.)
  tot_R = reproducing_sum(array, HI%isc, HI%iec, HI%jsc, HI%jec)
  tot_fastR = reproducing_sum(array, HI%isc, HI%iec, HI%jsc, HI%jec, overflow_check=.false.)

  ! tot_std and tot_R should differ only by round off, if at all
  if (abs(tot_std - tot_R) > likely_error) then
    write(mesg,'("Mismatch between standard and reproducing sum.",4ES13.5)') &
       tot_std, tot_R, tot_std - tot_R, ( tot_std - tot_R ) / tot_R
    call MOM_mesg(mesg)
    tests_failed = tests_failed .or. .true.
  endif
  ! tot_fastR and tot_R should be identical unless too many values are summed
  if (abs(tot_fastR - tot_R) > 0.) then
    if (n < max_count_prec) then
      write(mesg,'("Mismatch between reproducing and fast reproducing sums.",4ES13.5)') &
         tot_fastR, tot_R, tot_fastR - tot_R, ( tot_fastR - tot_R ) / tot_R
      tests_failed = tests_failed .or. .true.
    else
      write(mesg,'("Too many values were summed for the fast reproducing sum to work.")')
    endif
    call MOM_mesg(mesg)
  endif

  ! Now check the reproducing sums give the exact answer for known sets of values

  ! Fill array with values 1, 2, ..., Ni*Nj  whose sum is N ( N + 1 ) / 2 where N + Ni*Nj
  do j = HI%jsc, HI%jec ; do i = HI%isc, HI%iec
    jg = j + HI%jdg_offset - 1 ! 0 .. Nj-1
    ig = i + HI%idg_offset - 1 ! 0 .. Ni-1
    array(i,j) = 1 + ig + n_global(1) * jg
  enddo ; enddo
  tot_std = 0.5 * real(N) * real(N + 1) ! tot_std will contain analytic solution
  tot_R = reproducing_sum(array, HI%isc, HI%iec, HI%jsc, HI%jec)
  if (abs(tot_R - tot_std) > 0.) then
    write(mesg,'("Sum_k=1^N k != N(N+1)/2",2ES13.5)') tot_R, tot_std
    call MOM_mesg(mesg)
    tests_failed = tests_failed .or. .true.
  endif

  ! Change the order of values in the arrya to check the sum is truly order invariant
  do i = 1, n_repeat
    call randomly_swap_elements(HI, array)
    tot_R = reproducing_sum(array, HI%isc, HI%iec, HI%jsc, HI%jec)
    if (abs(tot_R - tot_std) > 0.) then
      write(mesg,'("Reordered list changed sum",2ES13.5)') tot_R, tot_std
      call MOM_mesg(mesg)
      tests_failed = tests_failed .or. .true.
    endif
  enddo

  call random_number( array ) ! This will also fill the halos but they will be ignored
  tot_std = reproducing_sum(array, HI%isc, HI%iec, HI%jsc, HI%jec) ! Use this as the true value
  ! Change the order of values in the arrya to check the sum is truly order invariant
  do i = 1, n_repeat
    call randomly_swap_elements(HI, array)
    tot_R = reproducing_sum(array, HI%isc, HI%iec, HI%jsc, HI%jec)
    if (abs(tot_R - tot_std) > 0.) then
      write(mesg,'("Reordered list of random numbers changed sum",2ES13.5)') tot_R, tot_std
      call MOM_mesg(mesg)
      tests_failed = tests_failed .or. .true.
    endif
  enddo

  ! Cleanup the "input.nml" file created to boot FMS
  if (PE_here() == root_PE()) then ! Can only delete the file once (i.e. on root PE)
    open(newunit=io_unit, file="input.nml", status="replace", action="write")
    close(io_unit, status="delete") ! we could leave this in place but that would be untidy
  endif

  call MOM_infra_end
  if (tests_failed) stop 1

contains

!> Randomly swap elements within the computational domain of an array
subroutine randomly_swap_elements(HI, array)
  type(hor_index_type), intent(in)  :: HI !< The horizontal index type
  real, intent(inout) :: array(HI%isd:HI%ied,HI%jsd:HI%jed) !< Array of values to play with [A]
  ! Local variables
  integer :: n_swaps !< Number of swaps to perform
  integer :: i0, j0, i1, j1, iter ! Indices and counter
  real :: r(4) ! Random numbers [nondim]
  real :: v ! Value being swapped

  n_swaps = ( HI%iec - HI%isc ) * ( HI%jec - HI%jsc )
  do iter = 1, n_swaps
    do
      call random_number( r ) ! Random numbers 0..1
      i0 = HI%isc + int( r(1) * real( HI%iec - HI%isc ) )
      j0 = HI%jsc + int( r(2) * real( HI%jec - HI%jsc ) )
      i1 = HI%isc + int( r(3) * real( HI%iec - HI%isc ) )
      j1 = HI%jsc + int( r(4) * real( HI%jec - HI%jsc ) )
      if (i0 /= i1 .and. j0 /= j1) exit ! Repeat dice roll if points are the same
    enddo
    v = array(i0,j0)
    array(i0,j0) = array(i1,j1)
    array(i1,j1) = v
  enddo
end subroutine randomly_swap_elements

!> Generate some "spatial" data, reminiscent of benchmark topography
subroutine generate_array_of_values(D, HI, n_global)
  type(hor_index_type), intent(in)  :: HI !< The horizontal index type
  real, intent(out) :: D(HI%isd:HI%ied,HI%jsd:HI%jed) !< Ocean bottom depth in [m]
  integer, intent(in) :: n_global(2) !< Global i-, j- dimensions of domain (h-points)
  ! Local variables
  real :: PI ! 3.1415926... calculated as 4*atan(1) [nondim]
  real :: x ! A fractional position in the x-direction [nondim]
  real :: y ! A fractional position in the y-direction [nondim]
  integer :: i, j ! Loop indices

  PI = 4.0*atan(1.0)

  !  Calculate the depth of the bottom.
  do concurrent( j=HI%jsc:HI%jec, i=HI%isc:HI%iec )
    x = real( i + HI%idg_offset ) / real( n_global(1) )
    y = real( j + HI%idg_offset ) / real( n_global(2) )
    D(i,j) = -3000.0  * ( y*(1.0 + 0.6*cos(4.0*PI*x)) &
                          + 0.75*exp(-6.0*y) &
                          + 0.05*cos(10.0*PI*x) - 0.7 )
    if (D(i,j) > 3000.0) D(i,j) = 3000.0
    if (D(i,j) < 1.) D(i,j) = 0.
  enddo

end subroutine generate_array_of_values

end program test_reproducing_sum
