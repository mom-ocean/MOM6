program time_MOM_remapping

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_remapping, only : remapping_CS
use MOM_remapping, only : initialize_remapping
use MOM_remapping, only : remapping_core_h

implicit none

type(remapping_CS) :: CS
integer, parameter :: nk=75, nij=20*20, nits=10, nsamp=100, nschemes = 2
character(len=10) :: scheme_labels(nschemes)
real, dimension(nschemes) :: timings ! Time for nits of nij calls for each scheme [s]
real, dimension(nschemes) :: tmean ! Mean time for a call [s]
real, dimension(nschemes) :: tstd ! Standard deviation of time for a call [s]
real, dimension(nschemes) :: tmin ! Shortest time for a call [s]
real, dimension(nschemes) :: tmax ! Longest time for a call [s]
real, dimension(:,:), allocatable :: u0, u1 ! Source/target values [arbitrary but same units as each other]
real, dimension(:,:), allocatable :: h0, h1 ! Source target thicknesses [0..1] [nondim]
real :: start, finish ! Times [s]
real :: h_neglect    ! A negligible thickness [nondim]
real :: h0sum, h1sum ! Totals of h0 and h1 [nondim]
integer :: ij, k, isamp, iter, ischeme ! Indices and counters
integer :: seed_size ! Number of integers used by seed
integer, allocatable :: seed(:) ! Random number seed

! Set seed for random numbers
call random_seed(size=seed_size)
allocate( seed(seed_Size) )
seed(:) = 102030405
call random_seed(put=seed)

scheme_labels(1) = 'PCM'
scheme_labels(2) = 'PLM'

! Set up some test data (note: using k,i indexing rather than i,k)
allocate( u0(nk,nij), h0(nk,nij), u1(nk,nij), h1(nk,nij) )
call random_number(u0) ! In range 0-1
call random_number(h0) ! In range 0-1
call random_number(h1) ! In range 0-1
do ij = 1, nij
  h0(:,ij) = max(0., h0(:,ij) - 0.05) ! Make 5% of values equal to zero
  h1(:,ij) = max(0., h1(:,ij) - 0.05) ! Make 5% of values equal to zero
  h0sum = h0(1,ij)
  h1sum = h1(1,ij)
  do k = 2, nk
    h0sum = h0sum + h0(k,ij)
    h1sum = h1sum + h1(k,ij)
  enddo
  h0(:,ij) = h0(:,ij) / h0sum
  h1(:,ij) = h1(:,ij) / h1sum
enddo
h_neglect = 1.0-30

! Loop over many samples of timing loop to collect statistics
tmean(:) = 0.
tstd(:) = 0.
tmin(:) = 1.e9
tmax(:) = 0.
do isamp = 1, nsamp
  ! Time reconstruction + remapping
  do ischeme = 1, nschemes
    call initialize_remapping(CS, remapping_scheme=trim(scheme_labels(ischeme)), &
                 h_neglect=h_neglect, h_neglect_edge=h_neglect)
    call cpu_time(start)
    do iter = 1, nits ! Make many passes to reduce sampling error
      do ij = 1, nij ! Calling nij times to make similar to cost in MOM_ALE()
        call remapping_core_h(CS, nk, h0(:,ij), u0(:,ij), nk, h1(:,ij), u1(:,ij))
      enddo
    enddo
    call cpu_time(finish)
    timings(ischeme) = (finish-start)/real(nits*nij) ! Average time per call
  enddo
  tmean(:) = tmean(:) + timings(:)
  tstd(:) = tstd(:) + timings(:)**2 ! tstd contains sum of squares here
  tmin(:) = min( tmin(:), timings(:) )
  tmax(:) = max( tmax(:), timings(:) )
enddo
tmean(:) = tmean(:) / real(nsamp) ! convert to mean
tstd(:) = tstd(:) / real(nsamp) ! convert to mean of squares
tstd(:) = tstd(:) - tmean(:)**2  ! convert to variance
tstd(:) = sqrt( tstd(:) * real(nsamp) / real(nsamp-1) ) ! convert to standard deviation


! Display results in YAML
write(*,'(a)') "{"
do ischeme = 1, nschemes
  write(*,"(2x,5a)") '"MOM_remapping remapping_core_h(remapping_scheme=', &
                     trim(scheme_labels(ischeme)), ')": {'
  write(*,"(4x,a,1pe11.4,',')") '"min": ',tmin(ischeme)
  write(*,"(4x,a,1pe11.4,',')") '"mean":',tmean(ischeme)
  write(*,"(4x,a,1pe11.4,',')") '"std": ',tstd(ischeme)
  write(*,"(4x,a,i7,',')") '"n_samples": ',nsamp
  if (ischeme.ne.nschemes) then
    write(*,"(4x,a,1pe11.4,'},')") '"max": ',tmax(ischeme)
  else
    write(*,"(4x,a,1pe11.4,'}')") '"max": ',tmax(ischeme)
  endif
enddo
write(*,'(a)') "}"

end program time_MOM_remapping
