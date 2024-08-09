program time_MOM_EOS

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_EOS, only : EOS_type
use MOM_EOS, only : EOS_manual_init
use MOM_EOS, only : calculate_density, calculate_spec_vol
use MOM_EOS, only : list_of_eos, get_EOS_name

implicit none

! This macro is used to write out timings of a single test rather than conduct
! a suite of tests. It is not meant for general consumption.
#undef PDF_ONLY

integer, parameter :: n_fns = 4
character(len=40) :: fn_labels(n_fns)

! Testing parameters:
!  nic is number of elements to compute density for (array size), per call
!  halo is data on either end of the array that should not be used
!  nits is how many times to repeat the call between turning the timer on/off
!       to overcome limited resolution of the timer
!  nsamp repeats the timing to collect statistics on the measurement
#ifdef PDF_ONLY
integer, parameter :: nic=26, halo=4, nits=10000, nsamp=400
#else
integer, parameter :: nic=23, halo=4, nits=1000, nsamp=400
#endif

real :: times(nsamp) ! CPU times for observing the PDF [seconds]

! Arrays to hold timings in [seconds]:
!  first axis corresponds to the form of EOS
!  second axis corresponds to the function being timed
real, dimension(:,:), allocatable :: timings, tmean, tstd, tmin, tmax
integer :: n_eos, i, j

n_eos = size(list_of_eos)
allocate( timings(n_eos,n_fns), tmean(n_eos,n_fns) )
allocate( tstd(n_eos,n_fns), tmin(n_eos,n_fns), tmax(n_eos,n_fns) )

fn_labels(1) = 'calculate_density_scalar()'
fn_labels(2) = 'calculate_density_array()'
fn_labels(3) = 'calculate_spec_vol_scalar()'
fn_labels(4) = 'calculate_spec_vol_array()'

tmean(:,:) = 0.
tstd(:,:) = 0.
tmin(:,:) = 1.e9
tmax(:,:) = 0.
do i = 1, nsamp
#ifdef PDF_ONLY
  call run_one(list_of_EOS, nic, halo, nits, times(i))
#else
  call run_suite(list_of_EOS, nic, halo, nits, timings)
  tmean(:,:) = tmean(:,:) + timings(:,:)
  tstd(:,:) = tstd(:,:) + timings(:,:)**2 ! tstd contains sum or squares here
  tmin(:,:) = min( tmin(:,:), timings(:,:) )
  tmax(:,:) = max( tmax(:,:), timings(:,:) )
#endif
enddo
tmean(:,:) = tmean(:,:) / real(nsamp)
tstd(:,:) = tstd(:,:) / real(nsamp) ! convert to mean of squares
tstd(:,:) = tstd(:,:) - tmean(:,:)**2 ! convert to variance
tstd(:,:) = sqrt( tstd(:,:) * ( real(nsamp) / real(nsamp-1) ) ) ! Standard deviation

#ifdef PDF_ONLY
open(newunit=i, file='times.txt', status='replace', action='write')
write(i,'(1pE9.3)') times(:)
close(i)
#else

! Display results in YAML
write(*,'(a)') "{"
do i = 1, n_eos
  do j = 1, n_fns
    write(*,"(2x,5a)") '"MOM_EOS_', trim(get_EOS_name(list_of_EOS(i))), &
                       ' ', trim(fn_labels(j)), '": {'
    write(*,"(4x,a,1pe11.4,',')") '"min": ',tmin(i,j)
    write(*,"(4x,a,1pe11.4,',')") '"mean":',tmean(i,j)
    write(*,"(4x,a,1pe11.4,',')") '"std": ',tstd(i,j)
    write(*,"(4x,a,i7,',')") '"n_samples": ',nsamp
    if (i*j.ne.n_eos*n_fns) then
      write(*,"(4x,a,1pe11.4,'},')") '"max": ',tmax(i,j)
    else
      write(*,"(4x,a,1pe11.4,'}')") '"max": ',tmax(i,j)
    endif
  enddo
enddo
write(*,'(a)') "}"
#endif

contains

subroutine run_suite(EOS_list, nic, halo, nits, timings)
  integer, intent(in)  :: EOS_list(n_eos) !< IDs of EOS forms to loop over
  integer, intent(in)  :: nic          !< Width of computational domain
  integer, intent(in)  :: halo         !< Width of halo to add on either end
  integer, intent(in)  :: nits         !< Number of calls to sample
                                       !! (large enough that the CPU timers can resolve
                                       !! the loop)
  real,    intent(out) :: timings(n_eos,n_fns) !< The average time taken for nits calls [seconds]
                                       !! First index corresponds to EOS
                                       !! Second index: 1 = scalar args,
                                       !! 2 = array args without halo,
                                       !! 3 = array args with halo and "dom".
  type(EOS_type) :: EOS
  integer :: e, i, dom(2)
  real :: start, finish  ! CPU times [seconds]
  real :: T  ! A potential or conservative temperature [degC]
  real :: S  ! A practical salinity or absolute salinity [ppt]
  real :: P  ! A pressure [Pa]
  real :: rho ! A density [kg m-3] or specific volume [m3 kg-1]
  real, dimension(nic+2*halo) :: T1, S1, P1, rho1

  T = 10.
  S = 35.
  P = 2000.e4

  ! Time the scalar interface
  do e = 1, n_eos
    call EOS_manual_init(EOS, form_of_EOS=EOS_list(e), &
                         Rho_T0_S0=1030., dRho_dT=0.2, dRho_dS=-0.7)

    call cpu_time(start)
    do i = 1, nits*nic ! Calling nic* to make similar cost to array call
      call calculate_density(T, S, P, rho, EOS)
    enddo
    call cpu_time(finish)
    timings(e,1) = (finish - start) / real(nits)

    call cpu_time(start)
    do i = 1, nits*nic ! Calling nic* to make similar cost to array call
      call calculate_spec_vol(T, S, P, rho, EOS)
    enddo
    call cpu_time(finish)
    timings(e,2) = (finish - start) / real(nits)

  enddo

  ! Time the "dom" interface, 1D array + halos
  T1(:) = T
  S1(:) = S
  P1(:) = P
  dom(:) = [1+halo,nic+halo]

  do e = 1, n_eos
    call EOS_manual_init(EOS, form_of_EOS=EOS_list(e), &
                         Rho_T0_S0=1030., dRho_dT=0.2, dRho_dS=-0.7)

    call cpu_time(start)
    do i = 1, nits
      call calculate_density(T1, S1, P1, rho1, EOS, dom)
    enddo
    call cpu_time(finish)
    timings(e,3) = (finish - start) / real(nits)

    call cpu_time(start)
    do i = 1, nits
      call calculate_spec_vol(T1, S1, P1, rho1, EOS, dom)
    enddo
    call cpu_time(finish)
    timings(e,4) = (finish - start) / real(nits)

  enddo

end subroutine run_suite

!> Return timing for just one fixed call to explore the PDF
subroutine run_one(EOS_list, nic, halo, nits, timing)
  integer, intent(in)  :: EOS_list(n_eos) !< IDs of EOS forms to loop over
  integer, intent(in)  :: nic          !< Width of computational domain
  integer, intent(in)  :: halo         !< Width of halo to add on either end
  integer, intent(in)  :: nits         !< Number of calls to sample
                                       !! (large enough that the CPU timers can resolve
                                       !! the loop)
  real,    intent(out) :: timing       !< The average time taken for nits calls [seconds]
                                       !! First index corresponds to EOS
                                       !! Second index: 1 = scalar args,
                                       !! 2 = array args without halo,
                                       !! 3 = array args with halo and "dom".
  type(EOS_type) :: EOS
  integer :: i, dom(2)
  real :: start, finish  ! CPU times [seconds]
  real, dimension(nic+2*halo) :: T1   ! Potential or conservative temperatures [degC]
  real, dimension(nic+2*halo) :: S1   ! A practical salinities or absolute salinities [ppt]
  real, dimension(nic+2*halo) :: P1   ! Pressures [Pa]
  real, dimension(nic+2*halo) :: rho1 ! Densities [kg m-3] or specific volumes [m3 kg-1]

  ! Time the scalar interface
  call EOS_manual_init(EOS, form_of_EOS=EOS_list(5), &
                       Rho_T0_S0=1030., dRho_dT=0.2, dRho_dS=-0.7)

  ! Time the "dom" interface, 1D array + halos
  T1(:) = 10.
  S1(:) = 35.
  P1(:) = 2000.e4
  dom(:) = [1+halo,nic+halo]

  call EOS_manual_init(EOS, form_of_EOS=EOS_list(5), &
                       Rho_T0_S0=1030., dRho_dT=0.2, dRho_dS=-0.7)

  call cpu_time(start)
  do i = 1, nits
    call calculate_density(T1, S1, P1, rho1, EOS, dom)
  enddo
  call cpu_time(finish)
  timing = (finish-start)/real(nits)

end subroutine run_one

end program time_MOM_EOS
