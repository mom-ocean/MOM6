!> Provides gridded random number capability
module MOM_random

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_hor_index,    only : hor_index_type
use MOM_time_manager, only : time_type, set_date, get_date

use iso_fortran_env,  only : stdout=>output_unit, stderr=>error_unit

implicit none ; private

public :: random_0d_constructor
public :: random_01
public :: random_01_CB
public :: random_norm
public :: random_2d_constructor
public :: random_2d_01
public :: random_2d_norm
public :: random_unit_tests

! Private period parameters for the Mersenne Twister
integer, parameter :: blockSize = 624,           & !< Size of the state vector
                      M         = 397,           & !< Pivot element in state vector
                      MATRIX_A  = -1727483681,   & !< constant vector a         (0x9908b0dfUL)
                      UMASK     = -2147483648_8, & !< most significant w-r bits (0x80000000UL)
                      LMASK     =  2147483647      !< least significant r bits  (0x7fffffffUL)
! Private tempering parameters for the Mersenne Twister
integer, parameter :: TMASKB= -1658038656, & !< (0x9d2c5680UL)
                      TMASKC= -272236544     !< (0xefc60000UL)

!> A private type used by the Mersenne Twistor
type randomNumberSequence
  integer                            :: currentElement !< Index into state vector
  integer, dimension(0:blockSize -1) :: state          !< State vector
end type randomNumberSequence

!> Container for pseudo-random number generators
type, public :: PRNG ; private

  !> Scalar random number generator for whole model
  type(randomNumberSequence) :: stream0d

  !> Random number generator for each cell on horizontal grid
  type(randomNumberSequence), dimension(:,:), allocatable :: stream2d

end type PRNG

contains

!> Returns a random number between 0 and 1
real function random_01(CS)
  type(PRNG), intent(inout) :: CS !< Container for pseudo-random number generators

  random_01 = getRandomReal(CS%stream0d)

end function random_01

!> Returns a random number between 0 and 1
!! See https://arxiv.org/abs/2004.06278. Not an exact reproduction of "squares" because Fortran
!! doesn't have a uint64 type, and not all compilers provide integers with > 64 bits...
real function random_01_CB(ctr, key)
  use iso_fortran_env, only : int64
  integer, intent(in)  :: ctr !< ctr should be incremented each time you call the function
  integer, intent(in)  :: key !< key is like a seed: use a different key for each random stream
  integer(kind=int64) :: x, y, z ! Follows "Squares" naming convention

  x = (ctr + 1) * (key + 65536) ! 65536 added because keys below that don't work.
  y = (ctr + 1) * (key + 65536)
  z = y + (key + 65536)
  x = x*x + y
  x = ior(ishft(x,32),ishft(x,-32))
  x = x*x + z
  x = ior(ishft(x,32),ishft(x,-32))
  x = x*x + y
  x = ior(ishft(x,32),ishft(x,-32))
  x = x*x + z
  random_01_CB = .5*(1. + .5*real(int(ishft(x,-32)))/real(2**30))

end function

!> Returns an approximately normally distributed random number with mean 0 and variance 1
real function random_norm(CS)
  type(PRNG), intent(inout) :: CS !< Container for pseudo-random number generators
  ! Local variables
  integer :: i

  random_norm = getRandomReal(CS%stream0d) - 0.5
  do i = 1,11
    random_norm = random_norm + ( getRandomReal(CS%stream0d) - 0.5 )
  enddo

end function random_norm

!> Generates random numbers between 0 and 1 for each cell of the model grid
subroutine random_2d_01(CS, HI, rand)
  type(PRNG),           intent(inout) :: CS !< Container for pseudo-random number generators
  type(hor_index_type), intent(in)    :: HI !< Horizontal index structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), intent(out) :: rand !< Random numbers between 0 and 1
  ! Local variables
  integer :: i,j

  do j = HI%jsd,HI%jed
    do i = HI%isd,HI%ied
      rand(i,j) = getRandomReal( CS%stream2d(i,j) )
    enddo
  enddo

end subroutine random_2d_01

!> Returns an approximately normally distributed random number with mean 0 and variance 1
!! for each cell of the model grid
subroutine random_2d_norm(CS, HI, rand)
  type(PRNG),           intent(inout) :: CS !< Container for pseudo-random number generators
  type(hor_index_type), intent(in)    :: HI !< Horizontal index structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), intent(out) :: rand !< Random numbers between 0 and 1
  ! Local variables
  integer :: i,j,n

  do j = HI%jsd,HI%jed
    do i = HI%isd,HI%ied
      rand(i,j) = getRandomReal( CS%stream2d(i,j) ) - 0.5
    enddo
    do n = 1,11
      do i = HI%isd,HI%ied
        rand(i,j) = rand(i,j) + ( getRandomReal( CS%stream2d(i,j) ) - 0.5 )
      enddo
    enddo
  enddo

end subroutine random_2d_norm

!> Constructor for scalar PRNG. Can be used to reset the sequence.
subroutine random_0d_constructor(CS, Time, seed)
  type(PRNG),      intent(inout) :: CS   !< Container for pseudo-random number generators
  type(time_type), intent(in)    :: Time !< Current model time
  integer,         intent(in)    :: seed !< Seed for PRNG
  ! Local variables
  integer :: tseed

  tseed = seed_from_time(Time)
  tseed = ieor(tseed, seed)
  CS%stream0d = new_RandomNumberSequence(tseed)

end subroutine random_0d_constructor

!> Constructor for gridded PRNG. Can be used to reset the sequence.
subroutine random_2d_constructor(CS, HI, Time, seed)
  type(PRNG),           intent(inout) :: CS   !< Container for pseudo-random number generators
  type(hor_index_type), intent(in)    :: HI   !< Horizontal index structure
  type(time_type),      intent(in)    :: Time !< Current model time
  integer,              intent(in)    :: seed !< Seed for PRNG
  ! Local variables
  integer :: i,j,sseed,tseed

  if (.not. allocated(CS%stream2d)) allocate( CS%stream2d(HI%isd:HI%ied,HI%jsd:HI%jed) )

  tseed = seed_from_time(Time)
  tseed = ieor(tseed*9007, seed)
  do j = HI%jsd,HI%jed
    do i = HI%isd,HI%ied
      sseed = seed_from_index(HI, i, j)
      sseed = ieor(tseed, sseed*7993)
      CS%stream2d(i,j) = new_RandomNumberSequence(sseed)
    enddo
  enddo

end subroutine random_2d_constructor

!> Return a seed derived as hash of values in Time
integer function seed_from_time(Time)
  type(time_type), intent(in)    :: Time !< Current model time
  ! Local variables
  integer :: yr,mo,dy,hr,mn,sc,s1,s2

  call get_date(Time,yr,mo,dy,hr,mn,sc)
  s1 = sc + 61*(mn + 61*hr) + 379 ! Range 379 .. 89620
  ! Fun fact: 2147483647 is the eighth Mersenne prime.
  ! This is not the reason for using 2147483647 here. It is the
  ! largest integer of kind=4.
  s2 = modulo(dy + 32*(mo + 13*yr), 2147483647_4) ! Range 0 .. 2147483646
  seed_from_time = ieor(s1*4111, s2)

end function seed_from_time

!> Create seed from position index
integer function seed_from_index(HI, i, j)
  type(hor_index_type), intent(in) :: HI !< Horizontal index structure
  integer,              intent(in) :: i !< i-index (of h-cell)
  integer,              intent(in) :: j !< j-index (of h-cell)
  ! Local variables
  integer :: ig, jg, ni, nj, ij

  ni = HI%niglobal
  nj = HI%njglobal
  ! Periodicity is assumed here but does not break non-periodic models
  ig = mod(HI%idg_offset + i - 1 + ni, ni)+1
  jg = max(HI%jdg_offset + j, 0)
  if (jg>nj) then ! Tri-polar hard-coded until we put needed info in HI **TODO**
    jg = 2*nj+1-jg
    ig = ni+1-ig
  endif
  seed_from_index = ig + ni*(jg-1)

end function seed_from_index

!> Destructor for PRNG
subroutine random_destruct(CS)
  type(PRNG), pointer :: CS !< Container for pseudo-random number generators

  if (allocated(CS%stream2d)) deallocate(CS%stream2d)
  !deallocate(CS)
end subroutine random_destruct

!> Return an initialized twister using seed
!!
!! Code was based on initialize_scaler() from the FMS implementation of the Mersenne Twistor
function new_RandomNumberSequence(seed) result(twister)
  integer, intent(in) :: seed !< Seed to initialize twister
  type(randomNumberSequence) :: twister !< The Mersenne Twister container
  ! Local variables
  integer :: i

  twister%state(0) = iand(seed, -1)
  do i = 1,  blockSize - 1 ! ubound(twister%state)
     twister%state(i) = 1812433253 * ieor(twister%state(i-1), &
                                          ishft(twister%state(i-1), -30)) + i
     twister%state(i) = iand(twister%state(i), -1) ! for >32 bit machines
  end do
  twister%currentElement = blockSize
end function new_RandomNumberSequence

!> Return a random integer on interval [0,0xffffffff]
!!
!! Code was based on getRandomInt() from the FMS implementation of the Mersenne Twistor
integer function getRandomInt(twister)
  type(randomNumberSequence), intent(inout) :: twister !< The Mersenne Twister container

  if(twister%currentElement >= blockSize) call nextState(twister)
  getRandomInt = temper(twister%state(twister%currentElement))
  twister%currentElement = twister%currentElement + 1

end function getRandomInt

!> Return a random real number on interval [0,1]
!!
!! Code was based on getRandomReal() from the FMS implementation of the Mersenne Twistor
double precision function getRandomReal(twister)
  type(randomNumberSequence), intent(inout) :: twister
  ! Local variables
  integer :: localInt

  localInt = getRandomInt(twister)
  if(localInt < 0) then
    getRandomReal = dble(localInt + 2.0d0**32)/(2.0d0**32 - 1.0d0)
  else
    getRandomReal = dble(localInt            )/(2.0d0**32 - 1.0d0)
  end if
end function getRandomReal

!> Merge bits of u and v
integer function mixbits(u, v)
  integer, intent(in) :: u !< An integer
  integer, intent(in) :: v !< An integer

  mixbits = ior(iand(u, UMASK), iand(v, LMASK))
end function mixbits

!> Twist bits of u and v
integer function twist(u, v)
  integer, intent(in) :: u !< An integer
  integer, intent(in) :: v !< An integer
  ! Local variable
  integer, parameter, dimension(0:1) :: t_matrix = (/ 0, MATRIX_A /)

  twist = ieor(ishft(mixbits(u, v), -1), t_matrix(iand(v, 1)))
  twist = ieor(ishft(mixbits(u, v), -1), t_matrix(iand(v, 1)))
end function twist

!> Update internal state of twister to the next state in the sequence
subroutine nextState(twister)
  type(randomNumberSequence), intent(inout) :: twister !< Container for the Mersenne Twister
  ! Local variables
  integer :: k

  do k = 0, blockSize - M - 1
    twister%state(k) = ieor(twister%state(k + M), &
                            twist(twister%state(k), twister%state(k + 1)))
  end do
  do k = blockSize - M, blockSize - 2
    twister%state(k) = ieor(twister%state(k + M - blockSize), &
                            twist(twister%state(k), twister%state(k + 1)))
  end do
  twister%state(blockSize - 1) = ieor(twister%state(M - 1), &
                                      twist(twister%state(blockSize - 1), twister%state(0)))
  twister%currentElement = 0
end subroutine nextState

!> Tempering of bits in y
elemental integer function temper(y)
  integer, intent(in) :: y !< An integer
  ! Local variables
  integer :: x

  x      = ieor(y, ishft(y, -11))
  x      = ieor(x, iand(ishft(x,  7), TMASKB))
  x      = ieor(x, iand(ishft(x, 15), TMASKC))
  temper = ieor(x, ishft(x, -18))
end function temper

!> Runs some statistical tests on the PRNG
logical function random_unit_tests(verbose)
  logical :: verbose !< True if results should be written to stdout
  ! Local variables
  type(PRNG) :: test_rng ! Generator
  type(time_type) :: Time ! Model time
  real :: r1, r2, r3 ! Some random numbers and re-used work variables
  real :: mean, var, ar1, std ! Some statistics
  integer :: stdunit ! For messages
  integer, parameter :: n_samples = 800
  integer :: i, j, ni, nj
  ! Fake being on a decomposed domain
  type(hor_index_type), pointer :: HI => null() !< Not the real HI
  real, dimension(:,:), allocatable :: r2d ! Random numbers

  ! Fake a decomposed domain
  ni = 6
  nj = 9
  allocate(HI)
  HI%isd = 0
  HI%ied = ni+1
  HI%jsd = 0
  HI%jed = nj+1
  HI%niglobal = ni
  HI%njglobal = nj
  HI%idg_offset = 0
  HI%jdg_offset = 0

  random_unit_tests = .false.
  stdunit = stdout
  write(stdunit,'(1x,a)') '==== MOM_random: random_unit_tests ======================='

  if (verbose) write(stdunit,'(1x,"random: ",a)') '-- Time-based seeds ---------------------'
  ! Check time-based seed generation
  Time = set_date(1903, 11, 21, 13, 47, 29)
  i = seed_from_time(Time)
  random_unit_tests = random_unit_tests .or. &
      test_fn(verbose, i==212584341, 'time seed 1903/11/21 13:47:29', ivalue=i)
  Time = set_date(1903, 11, 22, 13, 47, 29)
  i = seed_from_time(Time)
  random_unit_tests = random_unit_tests .or.&
       test_fn(verbose, i==212584342, 'time seed 1903/11/22 13:47:29', ivalue=i)
  Time = set_date(1903, 11, 21, 13, 47, 30)
  i = seed_from_time(Time)
  random_unit_tests = random_unit_tests .or.&
       test_fn(verbose, i==212596634, 'time seed 1903/11/21 13:47:30', ivalue=i)

  if (verbose) write(stdunit,'(1x,"random: ",a)') '-- PRNG tests ---------------------------'
  ! Generate a random number, r1
  call random_0d_constructor(test_rng, Time, 1)
  r1 = random_01(test_rng)
  random_unit_tests = random_unit_tests .or. &
       test_fn(verbose, abs(r1-4.75310122e-2)<1.e-9, 'first call', r1)

  ! Check that we get a different number, r2, on a second call
  r2 = random_01(test_rng)
  random_unit_tests = random_unit_tests .or. &
       test_fn(verbose, abs(r2-2.71289742e-1)<1.e-9, 'consecutive test', r2)

  ! Check that we can reproduce r1 by resetting the seed
  call random_0d_constructor(test_rng, Time, 1)
  r2 = random_01(test_rng)
  random_unit_tests = random_unit_tests .or. &
       test_fn(verbose, abs(r2-r1)==0., 'reproduce test', r2)

  ! Check that we get a different number, r2, with a different seed but same date
  call random_0d_constructor(test_rng, Time, 2)
  r2 = random_01(test_rng)
  random_unit_tests = random_unit_tests .or. &
       test_fn(verbose, abs(r2-7.15508473e-1)<1.e-9, 'different seed test', r2)

  ! Check that we get a different number, r2, for a different date but same seed
  Time = set_date(1903, 11, 21, 13, 0, 29)
  call random_0d_constructor(test_rng, Time, 1)
  r2 = random_01(test_rng)
  random_unit_tests = random_unit_tests .or. &
       test_fn(verbose, abs(r2-9.56667163e-1)<1.e-9, 'different date test', r2)

  if (verbose) write(stdunit,'(1x,"random: ",a)') '-- index-based seeds --------------------'
  ! Check index-based seed
  i = seed_from_index(HI,1,1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, i==1, 'seed from index (1,1)', ivalue=i)
  j = seed_from_index(HI,ni+1,1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, j==i, 'seed from index (n+1,1)', ivalue=j)
  i = seed_from_index(HI,ni,1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, i==6, 'seed from index (n,1)', ivalue=i)
  j = seed_from_index(HI,0,1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, j==i, 'seed from index (0,1)', ivalue=j)
  i = seed_from_index(HI,1,nj)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, i==49, 'seed from index (1,n)', ivalue=i)
  j = seed_from_index(HI,ni,nj+1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, j==i, 'seed from index (n,n+1)', ivalue=j)
  i = seed_from_index(HI,ni,nj)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, i==54, 'seed from index (n,n)', ivalue=i)
  j = seed_from_index(HI,1,nj+1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, j==i, 'seed from index (1,n+1)', ivalue=j)

  if (.not.random_unit_tests) write(stdunit,'(1x,a)') 'Passed unit tests'
  ! The rest of these are not unit tests but statistical tests and as such
  ! could fail for different sample sizes but happen to pass here.

  ! Check statistics of large samples for uniform generator
  mean = 0. ; var = 0. ; ar1 = 0. ; r2 = 0.
  do i = 1, n_samples
    r1 = random_01(test_rng) - 0.5
    mean = mean + r1
    var = var + r1**2
    ar1 = ar1 + r1*r2
    r2 = r1 ! Keep copy of last value
  enddo
  mean = mean / real(n_samples) ! Expected mean is 0
  var = var / real(n_samples) ! Expected variance is 1/12
  ar1 = ar1 / real(n_samples-1) ! Autocovariance
  std = sqrt(var) ! Expected std is sqrt(1/12)
  r2 = mean*sqrt(real(12*n_samples)) ! Normalized error in mean
  r3 = std*sqrt(12.) ! Normalized standard deviation
  r1 = ( ar1 * sqrt(real(n_samples-1)) ) / var
  if (verbose) then
    write(stdunit,'(1x,"random: ",a)') '-- Uniform -0.5 .. 0.5 generator --------'
    write(stdunit,'(1x,"random: ",a,f12.9)') 'mean =',mean,'std =',std,'AR1 =',ar1
    write(stdunit,'(1x,"random: ",a,f12.9)') 'norm. mean =',r2, &
                                             'norm. std =',r3,'norm. AR1 =',r1
  endif
  random_unit_tests = random_unit_tests .or. &
                      test_fn(verbose, abs(r2)<2., &
                              'n>>1, mean within 2 sigma [uniform]', r2)
  random_unit_tests = random_unit_tests .or. &
                      test_fn(verbose, abs(r3-1.)<1./sqrt(real(n_samples)), &
                              'n>>1, std ~ 1/sqrt(12) [uniform]', r3-1.)
  random_unit_tests = random_unit_tests .or. &
                      test_fn(verbose, abs(r1)<2., &
                              'n>>1, AR1 < std/sqrt(n) [uniform]', r1)

  ! Check statistics of large samples for normal generator
  mean = 0. ; var = 0. ; ar1 = 0. ; r2 = 0.
  do i = 1, n_samples
    r1 = random_norm(test_rng)
    mean = mean + r1
    var = var + r1**2
    ar1 = ar1 + r1*r2
    r2 = r1 ! Keep copy of last value for AR calculation
  enddo
  mean = mean / real(n_samples)
  var = var / real(n_samples)
  ar1 = ar1 / real(n_samples)
  std = sqrt(var)
  r3 = 1./sqrt(real(n_samples)) ! Standard error of mean
  r2 = mean*sqrt(real(n_samples)) ! Normalized error in mean
  r3 = std ! Normalized standard deviation
  r1 = ( ar1 * sqrt(real(n_samples-1)) ) / var
  if (verbose) then
    write(stdunit,'(1x,"random: ",a)') '-- Normal distribution generator --------'
    write(stdunit,'(1x,"random: ",a,f12.9)') 'mean =',mean,'std =',std,'AR1 =',ar1
    write(stdunit,'(1x,"random: ",a,f12.9)') 'norm. error in mean =',r2, &
                                             'norm. standard deviation =',r3,'norm. AR1 =',r1
  endif
  random_unit_tests = random_unit_tests .or. &
                      test_fn(verbose, abs(r2)<2., &
                              'n>>1, mean within 2 sigma [norm]', r2)
  random_unit_tests = random_unit_tests .or. &
                      test_fn(verbose, abs(r3-1.)<1./sqrt(real(n_samples)), &
                              'n>>1, std ~ 1 [norm]', r3-1.)
  random_unit_tests = random_unit_tests .or. &
                      test_fn(verbose, abs(r1)<2., &
                              'n>>1, AR1 < std/sqrt(n) [norm]', r1)

  if (verbose) write(stdunit,'(1x,"random: ",a)') '-- 2d PRNG ------------------------------'
  ! Check 2d random number generator 0..1
  allocate( r2d(HI%isd:HI%ied,HI%jsd:HI%jed) )
  call random_2d_constructor(test_rng, HI, Time, 123)
  r2d(:,:) = -999. ! Use -9. to detect unset values
  call random_2d_01(test_rng, HI, r2d)
  if (any(abs(r2d(:,:)+999.)<=0.)) random_unit_tests=.true.
  r1 = minval(r2d)
  r2 = maxval(r2d)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, r1>=0., '2d all set', r1)
  random_unit_tests = random_unit_tests .or. test_fn(verbose, r2<=1., '2d all valid', r2)
  mean = sum( r2d(1:ni,1:nj) - 0.5 )/real(ni*nj)
  var = sum( (r2d(1:ni,1:nj) - 0.5 - mean)**2 )/real(ni*nj)
  std = sqrt(var)
  r3 = 1./sqrt(real(12*ni*nj)) ! Standard error of mean
  r2 = mean*sqrt(real(12*ni*nj)) ! Normalized error in mean
  r3 = std*sqrt(12.) ! Normalized standard deviation
  if (verbose) then
    write(stdunit,'(1x,"random: ",a)') '2D uniform 0..1 generator'
    write(stdunit,'(1x,"random: ",a,f12.9)') 'mean =',mean,'std =',std
    write(stdunit,'(1x,"random: ",a,f12.9)') 'norm. error in mean =',r2
    write(stdunit,'(1x,"random: ",a,f12.9)') 'norm. standard deviation =',r3
  endif
  random_unit_tests = random_unit_tests .or. &
                      test_fn(verbose, abs(r2)<2., &
                              '2d, mean within 2 sigma [uniform]', r2)
  random_unit_tests = random_unit_tests .or. &
                      test_fn(verbose, abs(r3-1.)<1./sqrt(real(ni*nj)), &
                              '2d, std ~ 1/sqrt(12) [uniform]', r3-1.)
  if (verbose) then
    write(stdunit,'(1x,"random:")')
    write(stdunit,'(1x,"random:",8f8.5)') r2d
    write(stdunit,'(1x,"random:")')
  endif

  ! Check 2d normal random number generator
  call random_2d_norm(test_rng, HI, r2d)
  mean = sum( r2d(1:ni,1:nj) )/real(ni*nj)
  var = sum( r2d(1:ni,1:nj)**2 )/real(ni*nj)
  std = sqrt(var)
  r3 = 1./sqrt(real(ni*nj)) ! Standard error of mean
  r2 = mean*sqrt(real(ni*nj)) ! Normalized error in mean
  r3 = std ! Normalized standard deviation
  if (verbose) then
    write(stdunit,'(1x,"random: ",a)') '2D normal generator'
    write(stdunit,'(1x,"random: ",a,f12.9)') 'mean =',mean,'std =',std
    write(stdunit,'(1x,"random: ",a,f12.9)') 'norm. error in mean =',r2
    write(stdunit,'(1x,"random: ",a,f12.9)') 'norm. standard deviation =',r3
  endif
  random_unit_tests = random_unit_tests .or. &
                      test_fn(verbose, abs(r2)<2., &
                              '2d, mean within 2 sigma [norm]', r2)
  random_unit_tests = random_unit_tests .or. &
                      test_fn(verbose, abs(r3-1.)<1./sqrt(real(ni*nj)), &
                              '2d, std ~ 1/sqrt(12) [norm]', r3-1.)

  ! Clean up
  deallocate(r2d)
  deallocate(HI)

  if (.not.random_unit_tests) write(stdunit,'(1x,a)') 'Passed statistical tests'

end function random_unit_tests

!> Convenience function for reporting result of test
logical function test_fn(verbose, good, label, rvalue, ivalue)
  logical,          intent(in) :: verbose !< Verbosity
  logical,          intent(in) :: good !< True if pass, false otherwise
  character(len=*), intent(in) :: label !< Label for messages
  real,             intent(in) :: rvalue !< Result of calculation
  integer,          intent(in) :: ivalue !< Result of calculation
  optional :: rvalue, ivalue

  if (present(ivalue)) then
    if (.not. good) then
      write(stdout,'(1x,a,i10,1x,a,a)') 'random: result =',ivalue,label,' <------- FAIL!'
      write(stderr,'(1x,a,i10,1x,a,a)') 'random: result =',ivalue,label,' <------- FAIL!'
    elseif (verbose) then
      write(stdout,'(1x,a,i10,1x,a)') 'random: result =',ivalue,label
    endif
  else
    if (.not. good) then
      write(stdout,'(1x,a,1pe15.8,1x,a,a)') 'random: result =',rvalue,label,' <------- FAIL!'
      write(stderr,'(1x,a,1pe15.8,1x,a,a)') 'random: result =',rvalue,label,' <------- FAIL!'
    elseif (verbose) then
      write(stdout,'(1x,a,1pe15.8,1x,a)') 'random: result =',rvalue,label
    endif
  endif
  test_fn = .not. good

end function test_fn

end module MOM_random

!> \namespace mom_random
!!
!! Provides MOM6 implementation of the Mersenne Twistor, copied from the FMS implementation
!! which was originally written by Robert Pincus (Robert.Pincus@colorado.edu).
!! We once used the FMS implementation directly but since random numers do not need to be
!! infrastructure specific, and because MOM6 should be infrastructure agnostic, we have copied
!! the parts of MT that we used here.
!!
!! Example usage:
!! \code
!! type(PRNG) :: rng
!! real :: rn
!! call random_0d_constructor(rng, Time, seed) ! Call this each time-step
!! rn = random_01(rng)
!! rn = random_norm(rng)
!!
!! type(PRNG) :: rng
!! real, dimension(:,:) :: rn2d
!! call random_2d_constructor(rng, HI, Time, seed) ! Call this each time-step
!! call random_2d_01(rng, HI, rn2d)
!! call random_2d_norm(rng, HI, rn2d)
!!
!! Note: reproducibility across restarts is implemented by using time-derived
!! seeds to pass to the Mersenne twister. It is therefore important that any
!! PRNG type be re-initialized each time-step.
!! \endcode
