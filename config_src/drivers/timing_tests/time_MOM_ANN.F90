! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

program time_MOM_ANN

use MOM_ANN, only : ANN_CS
use MOM_ANN, only : ANN_allocate, ANN_apply, ANN_end
use MOM_ANN, only : ANN_apply_vector_orig, ANN_apply_vector_oi
use MOM_ANN, only : ANN_apply_array_sio
use MOM_ANN, only : ANN_random

implicit none

! Command line options
integer :: nargs ! Number of command line arguments
character(len=12) :: cmd_ln_arg !< Command line argument (if any)

! ANN parameters
integer :: nlayers ! Number of layers
integer :: nin ! Number of inputs
integer :: layer_width ! Width of hidden layers
integer :: nout ! Number of outputs
! Timing parameters
integer :: nsamp ! Number of measurements
integer :: nits ! Number of calls to time
integer :: nxy ! Spatial dimension

nlayers = 7 ; nin = 4 ; layer_width = 16 ; nout = 1 ! Deep network
!nlayers = 4 ; nin = 4 ; layer_width = 48 ; nout = 1 ! Shallow-wide network
!nlayers = 3 ; nin = 4 ; layer_width = 20 ; nout = 1 ! Small network

nsamp = 100
nits = 20000
!nits = 300000 ! Needed for robust measurements on small networks
nxy = 100 ! larger array
!nxy = 10 ! small array

! Optionally grab ANN and timing parameters from the command line
nargs = command_argument_count()
if (nargs==7) then
  call get_command_argument(1, cmd_ln_arg)
  read(cmd_ln_arg,*) nlayers
  call get_command_argument(2, cmd_ln_arg)
  read(cmd_ln_arg,*) nin
  call get_command_argument(3, cmd_ln_arg)
  read(cmd_ln_arg,*) layer_width
  call get_command_argument(4, cmd_ln_arg)
  read(cmd_ln_arg,*) nout
  call get_command_argument(5, cmd_ln_arg)
  read(cmd_ln_arg,*) nsamp
  call get_command_argument(6, cmd_ln_arg)
  read(cmd_ln_arg,*) nits
  call get_command_argument(7, cmd_ln_arg)
  read(cmd_ln_arg,*) nxy
endif

! Fastest variants on Intel Xeon W-2223 CPU @ 3.60GHz (gfortran-13.2 -O3)
!                   | vector(nxy=1)  |   nxy = 10   |   nxy = 100
! ----------------------------------------------------------------------------
! Small ANN         |   vector_oi    |   array_soi  |   array_sio
! Shallow-wide ANN  |   vector_oi    |   array_ois  |   array_sio
! Deep ANN          |   vector_oi    |   array_ois  |   array_sio

write(*,'(a)') "{"

call time_ANN(nlayers, nin, layer_width, nout, nsamp, nits, nxy, &
              0, "MOM_ANN:ANN_apply(vector)")
write(*,"(',')")
call time_ANN(nlayers, nin, layer_width, nout, nsamp, nits, nxy, &
              1, "MOM_ANN:ANN_apply_vector_orig(array)")
write(*,"(',')")
call time_ANN(nlayers, nin, layer_width, nout, nsamp, nits, nxy, &
              2, "MOM_ANN:ANN_apply_vector_oi(array)")
write(*,"(',')")
call time_ANN(nlayers, nin, layer_width, nout, nsamp, nits, nxy, &
              12, "MOM_ANN:ANN_apply_array_sio(array)")
write(*,"()")

write(*,'(a)') "}"

contains

!> Time ANN inference.
!!
!! Times are measured over the "nits effective calls" and appropriately scaled to the
!! time per call per single vector of input features. For array inputs, the number of
!! actual calls is reduced by the size of the array.  The timing measurement is repeated
!! "nsamp" times, to check the statistics of the timing measurement.
subroutine time_ANN(nlayers, nin, width, nout, nsamp, nits, nxy, impl, label)
  integer,          intent(in)  :: nlayers !< Number of layers
  integer,          intent(in)  :: nin     !< Number of inputs
  integer,          intent(in)  :: width   !< Width of hidden layers
  integer,          intent(in)  :: nout    !< Number of outputs
  integer,          intent(in)  :: nsamp   !< Number of measurements
  integer,          intent(in)  :: nits    !< Number of calls to time
  integer,          intent(in)  :: nxy     !< Spatial dimension
  integer,          intent(in)  :: impl    !< Implementation to time
  character(len=*), intent(in)  :: label   !< Label for YAML output
  ! Local variables
  type(ANN_CS) :: ANN ! ANN
  integer :: widths(nlayers) ! Width of each layer
  real :: x_s(nin) ! Inputs (just features) [nondim]
  real :: y_s(nin) ! Outputs (just features) [nondim]
  real :: x_fs(nin,nxy) ! Inputs (feature, space) [nondim]
  real :: y_fs(nin,nxy) ! Outputs (feature, space) [nondim]
  real :: x_sf(nin,nxy) ! Inputs (space, feature) [nondim]
  real :: y_sf(nin,nxy) ! Outputs (space, feature) [nondim]
  integer :: iter, samp ! Loop counters
  integer :: ij ! Horizontal loop index
  real :: start, finish, timing ! CPU times [s]
  real :: tmin, tmax, tmean, tstd ! Min, max, mean, and standard deviation, of CPU times [s]
  integer :: asamp ! Actual samples of timings
  integer :: aits ! Actual iterations
  real :: words_per_sec ! Operations per sec estimated from parameters [# s-1]

  widths(:) = width
  widths(1) = nin
  widths(nlayers) = nout

  call ANN_random(ANN, nlayers, widths)
  call random_number(x_fs)
  call random_number(x_sf)


  tmin = 1e9
  tmax = 0.
  tmean = 0.
  tstd = 0.
  asamp = nits ! Most cases below use this
  aits = nits / nxy ! Most cases below use this

  do samp = 1, nsamp
    select case (impl)
      case (0)
        aits = nits
        call cpu_time(start)
        do iter = 1, nits ! Make many passes to reduce sampling error
          call ANN_apply(x_s, y_s, ANN)
        enddo
        call cpu_time(finish)
      case (1)
        call cpu_time(start)
        do iter = 1, aits ! Make many passes to reduce sampling error
          do ij = 1, nxy
            call ANN_apply_vector_orig(x_fs(:,ij), y_fs(:,ij), ANN)
          enddo
        enddo
        call cpu_time(finish)
      case (2)
        call cpu_time(start)
        do iter = 1, aits ! Make many passes to reduce sampling error
          do ij = 1, nxy
            call ANN_apply_vector_oi(x_fs(:,ij), y_fs(:,ij), ANN)
          enddo
        enddo
        call cpu_time(finish)
      case (12)
        call cpu_time(start)
        do iter = 1, aits ! Make many passes to reduce sampling error
          call ANN_apply_array_sio(nxy, x_sf(:,:), y_sf(:,:), ANN)
        enddo
        call cpu_time(finish)
        asamp = nsamp * aits ! Account for working on whole arrays
    end select

    timing = ( finish - start ) / real(nits) ! Average time per call

    tmin = min( tmin, timing )
    tmax = max( tmax, timing )
    tmean = tmean + timing
    tstd = tstd + timing**2
  enddo

  tmean = tmean / real(nsamp)
  tstd = tstd / real(nsamp) ! convert to mean of squares
  tstd = tstd - tmean**2  ! convert to variance
  tstd = sqrt( tstd * real(nsamp) / real(nsamp-1) ) ! convert to standard deviation
  words_per_sec = ANN%parameters / ( tmean * 1024 * 1024 )

  write(*,"(2x,3a)") '"', trim(label), '": {'
  write(*,"(4x,a,1pe11.4,',')") '"min": ', tmin
  write(*,"(4x,a,1pe11.4,',')") '"mean":', tmean
  write(*,"(4x,a,1pe11.4,',')") '"std": ', tstd
  write(*,"(4x,a,i0,',')") '"n_samples": ', asamp
  write(*,"(4x,a,1pe11.4,',')") '"max": ', tmax
  write(*,"(4x,a,1pe11.4,'}')", advance="no") '"MBps": ', words_per_sec

end subroutine time_ANN

end program time_MOM_ANN
