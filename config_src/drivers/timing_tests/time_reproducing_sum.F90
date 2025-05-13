program time_reproducing_sum

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms, only : PE_here, root_PE, num_PEs, reproducing_sum
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_domains, only : MOM_domain_type, create_MOM_domain, MOM_infra_init, MOM_infra_end
use MOM_domains, only : MOM_define_layout
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, MOM_set_verbosity
use MOM_hor_index, only : hor_index_type, hor_index_init

  implicit none

  type(MOM_domain_type), pointer :: Domain => NULL() ! Ocean model domain
  type(hor_index_type) :: HI ! A hor_index_type for array extents
  real, allocatable, dimension(:) :: depth_tot_R, depth_tot_std, depth_tot_fastR ! Various sums of depths [m]
  real, allocatable :: array(:,:) ! An array with values to sum over [m]
  character(len=200) :: mesg ! String for messages
  integer :: num_sums ! Number of times to repeat the sum call
  integer :: n ! Loop counter
  integer :: io_unit ! i/o unit for creating input.nml (sigh)
  integer :: reproClock, fastreproClock, stdClock, initClock ! Clocks for each sum
  integer :: n_global(2) ! Global i-, j- dimensions of domain (h-points)
  integer :: layout(2)   ! PE count in i-, j- directions
  integer :: PEs_used    ! Number of PEs available to executable

  ! FMS requires the file "input.nml" to exist ...
  open(newunit=io_unit, file="input.nml", status="replace", action="write")
  close(io_unit) ! ... but an empty input.nml is sufficient

  call MOM_infra_init()

  ! These clocks are on the global pelist.
  initClock = cpu_clock_id( 'Initialization' )
  stdClock = cpu_clock_id( 'Standard Sums' )
  reproClock = cpu_clock_id( 'Reproducing Sums' )
  fastreproClock = cpu_clock_id( 'Fast Reproducing Sums' )
  num_sums = 100

  call cpu_clock_begin(initClock)
  ! Optionally use command-line to change size of the problem
  ! Usage:  ./executable [tile-size] [number-of-calls]
  n = command_argument_count()
  if (n==2) then
    call get_command_argument(1, mesg)
    read(mesg,*) n_global(1)
    n_global(2) = n_global(1)
    call get_command_argument(2, mesg)
    read(mesg,*) num_sums
  elseif (n==1) then
    call get_command_argument(1, mesg)
    read(mesg,*) n_global(1)
    n_global(2) = n_global(1)
  else
    n_global = (/500, 300/) ! Fallback value if no argument provided
  endif

  call MOM_mesg('======== Unit test being driven by MOM_sum_driver ========', 2)
  call MOM_set_verbosity(2)

  ! Setup distributed domain
  PEs_used = num_PEs()
  call MOM_define_layout(n_global, PEs_used, layout)
  call create_MOM_domain(Domain, n_global, (/2,2/), (/.false.,.false./), .false., layout)
  call hor_index_init(Domain, HI)

  allocate( array(HI%isd:HI%ied,HI%jsd:HI%jed), source=0. )
  allocate( depth_tot_std(num_sums), source=0. )
  allocate( depth_tot_R(num_sums), source=0. )
  allocate( depth_tot_fastR(num_sums), source=0. )

  ! Set up an array of values to sum
  call generate_array_of_values(array, HI, n_global)

  call cpu_clock_end(initClock) !end initialization
  call MOM_mesg("Done with initialization.", 5)

  call MOM_mesg('==== Standard Non-reproducing Sum ===', 2)
  do n=1,num_sums
    call cpu_clock_begin(stdClock)
    depth_tot_std(n) = reproducing_sum(array, HI%isc, HI%iec, HI%jsc, HI%jec, reproducing=.false.)
    call cpu_clock_end(stdClock)
  enddo

  call MOM_mesg('==== Reproducing Fixed Point Sum ===', 2)
  do n=1,num_sums
    call cpu_clock_begin(reproClock)
    depth_tot_R(n) = reproducing_sum(array, HI%isc, HI%iec, HI%jsc, HI%jec)
    call cpu_clock_end(reproClock)
  enddo

  call MOM_mesg('==== No Error Handling Reproducing Fixed Point Sum ===', 2)
  do n=1,num_sums
    call cpu_clock_begin(fastreproClock)
    depth_tot_fastR(n) = reproducing_sum(array, HI%isc, HI%iec, HI%jsc, HI%jec, overflow_check=.false.)
    call cpu_clock_end(fastreproClock)
  enddo

  ! Cleanup the "input.nml" file created to boot FMS
  if (PE_here() == root_PE()) then ! Can only delete the file once (i.e. on root PE)
    open(newunit=io_unit, file="input.nml", status="replace", action="write")
    close(io_unit, status="delete") ! we could leave this in place but that would be untidy
  endif

  call MOM_infra_end

contains

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

end program time_reproducing_sum
