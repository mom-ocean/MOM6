!> Streaming band-pass filter for detecting the instantaneous tidal signals in the simulation

module MOM_streaming_filter

use MOM_error_handler, only : MOM_mesg, MOM_error, NOTE, FATAL
use MOM_file_parser,   only : get_param, param_file_type
use MOM_hor_index,     only : hor_index_type
use MOM_io,            only : axis_info, set_axis_info
use MOM_restart,       only : register_restart_field, query_initialized, MOM_restart_CS
use MOM_tidal_forcing, only : tidal_frequency
use MOM_time_manager,  only : time_type, time_type_to_real
use MOM_unit_scaling,  only : unit_scale_type

implicit none ; private

public Filt_register, Filt_init, Filt_accum

#include <MOM_memory.h>

!> Control structure for the MOM_streaming_filter module
type, public :: Filter_CS ; private
  integer :: nf                        !< Number of filters to be used in the simulation
  !>@{ Lower and upper bounds of input data
  integer :: is, ie, js, je
  !>@}
  character(len=8) :: key              !< Identifier of the variable to be filtered
  character(len=2), allocatable, dimension(:) :: filter_names !< Names of filters
  real, allocatable, dimension(:)      :: filter_omega !< Target frequencies of filters [rad T-1 ~> rad s-1]
  real, allocatable, dimension(:)      :: filter_alpha !< Bandwidth parameters of filters [nondim]
  real, allocatable, dimension(:,:,:)  :: s1, &        !< A dummy variable for solving the system of ODEs [A]
                                          u1           !< Filtered data, representing the narrow-band signal
                                                       !< oscillating around the target frequency [A]
  real :: old_time = -1.0              !< The time of the previous accumulating step [T ~> s]
end type Filter_CS

contains

!> This subroutine registers the filter variables given the number of filters and the grid
subroutine Filt_register(nf, key, grid, HI, CS, restart_CS)
  integer,               intent(in)    :: nf           !< Number of filters to be used in the simulation
  character(len=*),      intent(in)    :: key          !< Identifier of the variable to be filtered
  character(len=*),      intent(in)    :: grid         !< Horizontal grid location: "h", "u", or "v"
  type(hor_index_type),  intent(in)    :: HI           !< Horizontal index type structure
  type(Filter_CS),       intent(out)   :: CS           !< Control structure of MOM_streaming_filter
  type(MOM_restart_CS),  intent(inout) :: restart_CS   !< MOM restart control structure

  ! Local variables
  type(axis_info) :: filter_axis(1)
  real, dimension(:), allocatable :: n_filters         !< Labels of filters [nondim]
  integer :: c

  CS%nf  = nf
  CS%key = key

  select case (trim(grid))
    case ('h')
      CS%is = HI%isd  ; CS%ie = HI%ied  ; CS%js = HI%jsd  ; CS%je = HI%jed
    case ('u')
      CS%is = HI%IsdB ; CS%ie = HI%IedB ; CS%js = HI%jsd  ; CS%je = HI%jed
    case ('v')
      CS%is = HI%isd  ; CS%ie = HI%ied  ; CS%js = HI%JsdB ; CS%je = HI%JedB
    case default
      call MOM_error(FATAL, "MOM_streaming_filter: horizontal grid not supported")
  end select

  allocate(CS%s1(CS%is:CS%ie, CS%js:CS%je, nf), source=0.0)
  allocate(CS%u1(CS%is:CS%ie, CS%js:CS%je, nf), source=0.0)

  ! Register restarts for s1 and u1
  allocate(n_filters(nf))

  do c=1,nf ; n_filters(c) = c ; enddo

  call set_axis_info(filter_axis(1), "n_filters", "", "number of filters", nf, n_filters, "N", 1)

  call register_restart_field(CS%s1(:,:,:), "Filter_"//trim(key)//"_s1", .false., restart_CS, &
                              longname="Dummy variable for streaming band-pass filter", &
                              hor_grid=trim(grid), z_grid="1", t_grid="s", extra_axes=filter_axis)
  call register_restart_field(CS%u1(:,:,:), "Filter_"//trim(key)//"_u1", .false., restart_CS, &
                              longname="Output of streaming band-pass filter", &
                              hor_grid=trim(grid), z_grid="1", t_grid="s", extra_axes=filter_axis)

end subroutine Filt_register

!> This subroutine initializes the filters
subroutine Filt_init(param_file, US, CS, restart_CS)
  type(param_file_type), intent(in)    :: param_file   !< A structure to parse for run-time parameters
  type(unit_scale_type), intent(in)    :: US           !< A dimensional unit scaling type
  type(Filter_CS),       intent(inout) :: CS           !< Control structure of MOM_streaming_filter
  type(MOM_restart_CS),  intent(in)    :: restart_CS   !< MOM restart control structure

  ! Local variables
  character(len=40)  :: mdl = "MOM_streaming_filter"   !< This module's name
  character(len=50)  :: filter_name_str                !< List of filters to be registered
  character(len=200) :: mesg
  integer :: c

  call get_param(param_file, mdl, "FILTER_NAMES", filter_name_str, &
                 "Names of streaming band-pass filters to be used in the simulation.", &
                 fail_if_missing=.true.)
  allocate(CS%filter_names(CS%nf))
  allocate(CS%filter_omega(CS%nf))
  allocate(CS%filter_alpha(CS%nf))
  read(filter_name_str, *) CS%filter_names

  do c=1,CS%nf
    ! If filter_name_str consists of tidal constituents, use tidal frequencies.
    call get_param(param_file, mdl, "FILTER_"//trim(CS%filter_names(c))//"_OMEGA", &
                   CS%filter_omega(c), "Target frequency of the "//trim(CS%filter_names(c))//&
                   " filter. This is used if USE_FILTER is true and "//trim(CS%filter_names(c))//&
                   " is in FILTER_NAMES.", units="rad s-1", scale=US%T_to_s, default=0.0)
    call get_param(param_file, mdl, "FILTER_"//trim(CS%filter_names(c))//"_ALPHA", &
                   CS%filter_alpha(c), "Bandwidth parameter of the "//trim(CS%filter_names(c))//&
                   " filter. Must be positive.", units="nondim", fail_if_missing=.true.)

    if (CS%filter_omega(c)<=0.0) CS%filter_omega(c) = tidal_frequency(trim(CS%filter_names(c)))
    if (CS%filter_alpha(c)<=0.0) call MOM_error(FATAL, "MOM_streaming_filter: bandwidth <= 0")

    write(mesg,*) "MOM_streaming_filter: ", trim(CS%filter_names(c)), &
                  " filter registered, target frequency = ", CS%filter_omega(c), &
                  ", bandwidth = ", CS%filter_alpha(c)
    call MOM_error(NOTE, trim(mesg))
  enddo

  if (query_initialized(CS%s1, "Filter_"//trim(CS%key)//"_s1", restart_CS)) then
    write(mesg,*) "MOM_streaming_filter: Dummy variable for filter ", trim(CS%key), &
                  " found in restart files."
  else
    write(mesg,*) "MOM_streaming_filter: Dummy variable for filter ", trim(CS%key), &
                  " not found in restart files. The filter will spin up from zeros."
  endif
  call MOM_error(NOTE, trim(mesg))

  if (query_initialized(CS%u1, "Filter_"//trim(CS%key)//"_u1", restart_CS)) then
    write(mesg,*) "MOM_streaming_filter: Output of filter ", trim(CS%key), &
                  " found in restart files."
  else
    write(mesg,*) "MOM_streaming_filter: Output of filter ", trim(CS%key), &
                  " not found in restart files. The filter will spin up from zeros."
  endif
  call MOM_error(NOTE, trim(mesg))

end subroutine Filt_init

!> This subroutine timesteps the filter equations. Here, u is the broadband input signal from the model,
!! and u1 is the filtered, narrowband output signal, obtained from the solution of the filter equations.
subroutine Filt_accum(u, u1, Time, US, CS)
  real, dimension(:,:,:), pointer, intent(out)   :: u1   !< Output of the filter [A]
  type(time_type),                 intent(in)    :: Time !< The current model time
  type(unit_scale_type),           intent(in)    :: US   !< A dimensional unit scaling type
  type(Filter_CS),        target,  intent(inout) :: CS   !< Control structure of MOM_streaming_filter
  real, dimension(CS%is:CS%ie,CS%js:CS%je), intent(in) :: u !< Input into the filter [A]

  ! Local variables
  real    :: now, &              !< The current model time [T ~> s]
             dt, &               !< Time step size for the filter equations [T ~> s]
             c1, c2              !< Coefficients for the filter equations [nondim]
  integer :: i, j, k

  now = US%s_to_T * time_type_to_real(Time)

  ! Initialize CS%old_time at the first time step
  if (CS%old_time<0.0) CS%old_time = now

  ! Timestep the filter equations only if we are in a new time step
  if (CS%old_time<now) then
    dt = now - CS%old_time
    CS%old_time = now

    do k=1,CS%nf
      c1 = CS%filter_omega(k) * dt
      c2 = 1.0 - CS%filter_alpha(k) * c1

      do j=CS%js,CS%je ; do i=CS%is,CS%ie
        CS%s1(i,j,k) =  c1 *  CS%u1(i,j,k) + CS%s1(i,j,k)
        CS%u1(i,j,k) = -c1 * (CS%s1(i,j,k) - CS%filter_alpha(k) * u(i,j)) + c2 * CS%u1(i,j,k)
      enddo; enddo
    enddo ! k=1,CS%nf
  endif ! (CS%old_time<now)

  u1 => CS%u1

end subroutine Filt_accum

!> \namespace mom_streaming_filter
!!
!! By Chengzhu Xu (chengzhu.xu@oregonstate.edu) and Edward D. Zaron
!!
!! The algorithm detects the instantaneous, narrowband tidal signals (u1) from the broadband
!! model output (u) by solving a set of coupled ODEs (the filter equations) at each time step.
!! In the filter equations, u1 is approximately the part of the signal that oscillates at the
!! filter's target frequency, and s1 is approximately the imaginary complement in time of u1.
!!
!! Major revision on Dec 9, 2024: The filters are no longer hard-coded. Instead, multiple filters
!! with tidal frequencies or arbitrary frequencies as their target frequencies can be turned on.
!! The filter names are specified in MOM_input and must consist of two letters/numbers. If the
!! name of a filter is the same as the name of a tidal constituent, then the corresponding tidal
!! frequency will be used as its target frequency. Otherwise, the user must specify the target
!! frequency. In either case, the target frequency is specified by "FILTER_${FILTER_NAME}_OMEGA".
!!
!! The restarting capability has also been implemented. Because the filtering is a point-wise
!! operation, all variables are considered as fields, even if they are velocity components.
!!
!! Xu, C., & Zaron, E. D. (2024). Detecting instantaneous tidal signals in ocean models utilizing
!! streaming band-pass filters. Journal of Advances in Modeling Earth Systems, 16, e2024MS004319.
!! https://doi.org/10.1029/2024MS004319

end module MOM_streaming_filter

