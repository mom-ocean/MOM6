program MOM_main
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*                  The Modular Ocean Model                            *
!*                               MOM                                   *
!*                                                                     *
!*  By Robert Hallberg                                                 *
!*                                                                     *
!*    This file is a simple driver for unit testing the distributed    *
!*  sums code.                                                         *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

  use MOM_coms, only : sum_across_PEs, PE_here, root_PE, num_PEs, reproducing_sum
  use MOM_coms, only : EFP_type, operator(+), operator(-), assignment(=), EFP_to_real, real_to_EFP
  use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
  use MOM_cpu_clock, only : CLOCK_COMPONENT
!  use MOM_diag_mediator, only : diag_mediator_end, diag_mediator_init
!  use MOM_diag_mediator, only : diag_mediator_close_registration
  use MOM_domains, only : MOM_domains_init, MOM_infra_init, MOM_infra_end
  use MOM_error_handler, only : MOM_error, MOM_mesg, WARNING, FATAL, is_root_pe
  use MOM_error_handler, only : MOM_set_verbosity
  use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
  use MOM_file_parser, only : open_param_file, close_param_file
  use MOM_grid, only : MOM_grid_init, ocean_grid_type
  use MOM_grid_initialize, only : set_grid_metrics
  use MOM_io, only : MOM_io_init, file_exists, open_file, close_file
  use MOM_io, only : check_nml_error, io_infra_init, io_infra_end
  use MOM_io, only : APPEND_FILE, ASCII_FILE, READONLY_FILE, SINGLE_FILE

  implicit none

#include <MOM_memory.h>

  type(ocean_grid_type) :: grid ! A structure containing metrics and grid info.

  type(param_file_type) :: param_file ! The structure indicating the file(s)
                                ! containing all run-time parameters.
  real    :: max_depth
  integer :: verbosity
  integer :: num_sums
  integer :: n, i, j, is, ie, js, je, nz
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  integer :: unit, io_status, ierr
  logical :: unit_in_use

  real, allocatable, dimension(:) :: &
    depth_tot_R, depth_tot_std, depth_tot_fastR
  integer :: reproClock, fastreproClock, stdClock, initClock

  !-----------------------------------------------------------------------

  character(len=4), parameter :: vers_num = 'v2.0'
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_main (MOM_sum_driver)" ! This module's name.
  character(len=200) :: mesg

  !#######################################################################

  call MOM_infra_init() ; call io_infra_init()

  ! These clocks are on the global pelist.
  initClock = cpu_clock_id( 'Initialization' )
  reproClock = cpu_clock_id( 'Reproducing Sums' )
  fastreproClock = cpu_clock_id( 'Fast Reproducing Sums' )
  stdClock = cpu_clock_id( 'Standard Sums' )

  call cpu_clock_begin(initClock)

  call MOM_mesg('======== Unit test being driven by MOM_sum_driver ========', 2)

  call open_param_file("./MOM_input", param_file)

  verbosity = 2 ; call read_param(param_file, "VERBOSITY", verbosity)
  call MOM_set_verbosity(verbosity)

  call MOM_domains_init(grid%domain, param_file)

  call MOM_io_init(param_file)
!  call diag_mediator_init(param_file)
  call MOM_grid_init(grid, param_file)

  is = grid%isc ; ie = grid%iec ; js = grid%jsc ; je = grid%jec ; nz = grid%ke
  isd = grid%isd ; ied = grid%ied ; jsd = grid%jsd ; jed = grid%jed
  IsdB = grid%IsdB ; IedB = grid%IedB ; JsdB = grid%JsdB ; JedB = grid%JedB

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, "MOM", version, "")
  call get_param(param_file, "MOM", "VERBOSITY", verbosity,  &
                 "Integer controlling level of messaging\n" // &
                 "\t0 = Only FATAL messages\n" // &
                 "\t2 = Only FATAL, WARNING, NOTE [default]\n" // &
                 "\t9 = All)", default=2)
  call get_param(param_file, "MOM", "NUMBER_OF_SUMS", num_sums, &
                 "The number of times to do the global sums.", default=1)

  allocate(depth_tot_R(num_sums))     ; depth_tot_R(:) = 0.0
  allocate(depth_tot_std(num_sums))   ; depth_tot_std(:) = 0.0
  allocate(depth_tot_fastR(num_sums)) ; depth_tot_fastR(:) = 0.0

! Set up the parameters of the physical domain (i.e. the grid), G
  call set_grid_metrics(grid, param_file)

! Set up the bottom depth, grid%bathyT either analytically or from file
  call get_param(param_file, "MOM", "MAXIMUM_DEPTH", max_depth, &
                 "The maximum depth of the ocean.", units="m", default=4000.0)
  call benchmark_init_topog_local(grid%bathyT, grid, param_file, max_depth)

  ! Close the param_file.  No further parsing of input is possible after this.
  call close_param_file(param_file)

  call cpu_clock_end(initClock) !end initialization
  call MOM_mesg("Done with initialization.", 5)

  call MOM_mesg('==== Reproducing Fixed Point Sum ===', 2)

  call cpu_clock_begin(reproClock)
  do n=1,num_sums
    depth_tot_R(n) = reproducing_sum(grid%bathyT, is, ie, js, je)
  enddo
  call cpu_clock_end(reproClock)

  call MOM_mesg('==== Standard Non-reproducing Sum ===', 2)

  call cpu_clock_begin(stdClock)
!  do n=1,num_sums
!    do j=js,je ; do i=is,ie
!      depth_tot_std(n) = depth_tot_std(n) + grid%bathyT(i,j)
!    enddo ; enddo
!    call sum_across_PEs(depth_tot_std(n:),1)
!  enddo
  do n=1,num_sums
    depth_tot_fastR(n) = reproducing_sum(grid%bathyT, is, ie, js, je, reproducing=.false.)
  enddo
  call cpu_clock_end(stdClock)

  call MOM_mesg('==== No Error Handling Reproducing Fixed Point Sum ===', 2)

  call cpu_clock_begin(fastreproClock)
  do n=1,num_sums
    depth_tot_fastR(n) = reproducing_sum(grid%bathyT, is, ie, js, je, overflow_check=.false.)
  enddo
  call cpu_clock_end(fastreproClock)

  do n=1,num_sums
    if ((depth_tot_std(n) - depth_tot_R(n)) > 1e-15*depth_tot_R(n)) then
      write(mesg,'("Mismatch between standard and reproducing sum.",2ES13.5)') &
         depth_tot_std(n) - depth_tot_R(n),  depth_tot_R(n)
      call MOM_mesg(mesg) ; exit
    endif
    if ((depth_tot_fastR(n) - depth_tot_R(n)) > 1e-15*depth_tot_R(n)) then
      write(mesg,'("Mismatch between reproducing and fast reproducing sums.",2ES13.5)') &
         depth_tot_fastR(n) - depth_tot_R(n),  depth_tot_R(n)
      call MOM_mesg(mesg) ; exit
!       call MOM_mesg("Mismatch between reproducing and fast reproducing sums.")
    endif
  enddo

  call io_infra_end ; call MOM_infra_end

contains

subroutine benchmark_init_topog_local(D, G, param_file, max_depth)
  type(ocean_grid_type), intent(in)             :: G
  real, intent(out), dimension(SZI_(G),SZJ_(G)) :: D
  type(param_file_type), intent(in)             :: param_file
  real,                  intent(in)             :: max_depth
! Arguments: D          - the bottom depth in m. Intent out.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets up the benchmark test case topography
  real :: min_depth            ! The minimum and maximum depths in m.
  real :: PI                   ! 3.1415926... calculated as 4*atan(1)
  real :: D0                   ! A constant to make the maximum     !
                               ! basin depth MAXIMUM_DEPTH.         !
  real :: x, y
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "benchmark_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call MOM_mesg("  benchmark_initialization.F90, benchmark_initialize_topography: setting topography", 5)

  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)

  PI = 4.0*atan(1.0)
  D0 = max_depth / 0.5;

!  Calculate the depth of the bottom.
  do i=is,ie ; do j=js,je
    x=(G%geoLonT(i,j)-G%west_lon)/G%len_lon
    y=(G%geoLatT(i,j)-G%south_lat)/G%len_lat
!  This sets topography that has a reentrant channel to the south.
    D(i,j) = -D0 * ( y*(1.0 + 0.6*cos(4.0*PI*x)) &
                   + 0.75*exp(-6.0*y) &
                   + 0.05*cos(10.0*PI*x) - 0.7 )
    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.
  enddo ; enddo

end subroutine benchmark_init_topog_local

end program MOM_main
