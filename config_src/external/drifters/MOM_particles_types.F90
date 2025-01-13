!> Dummy data structures and methods for drifters package
module particles_types_mod

! This file is part of MOM6. See LICENSE.md for the license.

use, intrinsic :: iso_fortran_env, only : int64
use MOM_grid, only : ocean_grid_type
use MOM_domains, only: domain2D

implicit none ; private

!> Container for gridded fields
type, public :: particles_gridded
  type(domain2D), pointer :: domain !< MPP parallel domain
  integer :: halo !< Nominal halo width
  integer :: isc !< Start i-index of computational domain
  integer :: iec !< End i-index of computational domain
  integer :: jsc !< Start j-index of computational domain
  integer :: jec !< End j-index of computational domain
  integer :: isd !< Start i-index of data domain
  integer :: ied !< End i-index of data domain
  integer :: jsd !< Start j-index of data domain
  integer :: jed !< End j-index of data domain
  integer :: isg !< Start i-index of global domain
  integer :: ieg !< End i-index of global domain
  integer :: jsg !< Start j-index of global domain
  integer :: jeg !< End j-index of global domain
  integer :: is_offset=0 !< add to i to recover global i-index
  integer :: js_offset=0 !< add to j to recover global j-index
  integer :: my_pe !< MPI PE index
  integer :: pe_N !< MPI PE index of PE to the north
  integer :: pe_S !< MPI PE index of PE to the south
  integer :: pe_E !< MPI PE index of PE to the east
  integer :: pe_W !< MPI PE index of PE to the west
  logical :: grid_is_latlon !< Flag to say whether the coordinate is in lat-lon degrees, or meters
  logical :: grid_is_regular !< Flag to say whether point in cell can be found assuming regular Cartesian grid
  real :: Lx !< Length of the domain in x direction
  real, dimension(:,:), allocatable :: lon !< Longitude of cell corners (degree E)
  real, dimension(:,:), allocatable :: lat !< Latitude of cell corners (degree N)
  real, dimension(:,:), allocatable :: lonc !< Longitude of cell centers (degree E)
  real, dimension(:,:), allocatable :: latc !< Latitude of cell centers (degree N)
  real, dimension(:,:), allocatable :: dx !< Length of cell edge (m)
  real, dimension(:,:), allocatable :: dy !< Length of cell edge (m)
  real, dimension(:,:), allocatable :: area !< Area of cell (m^2)
  real, dimension(:,:), allocatable :: msk !< Ocean-land mask (1=ocean)
  real, dimension(:,:), allocatable :: cos !< Cosine from rotation matrix to lat-lon coords
  real, dimension(:,:), allocatable :: sin !< Sine from rotation matrix to lat-lon coords
  real, dimension(:,:), allocatable :: ocean_depth !< Depth of ocean (m)
  real, dimension(:,:), allocatable :: uo !< Ocean zonal flow (m/s)
  real, dimension(:,:), allocatable :: vo !< Ocean meridional flow (m/s)
  real, dimension(:,:), allocatable :: tmp !< Temporary work space
  real, dimension(:,:), allocatable :: tmpc !< Temporary work space
  real, dimension(:,:), allocatable :: parity_x !< X component of vector point from i,j to i+1,j+1
  real, dimension(:,:), allocatable :: parity_y !< Y component of vector point from i,j to i+1,j+1
  ! (For detecting tri-polar fold)
  integer, dimension(:,:), allocatable :: particle_counter_grd !< Counts particles created for naming purposes
  !>@{
  !! Diagnostic handle
  integer :: id_uo=-1, id_vo=-1, id_unused=-1
  integer :: id_count=-1, id_chksum=-1
  !>@}

end type particles_gridded


!>xyt is a data structure containing particle position and velocity fields.
type, public :: xyt
  real :: lon !< Longitude of particle (degree N or unit of grid coordinate)
  real :: lat !< Latitude of particle (degree N or unit of grid coordinate)
  real :: day      !< Day of this record (days)
  real :: lat_old  !< Previous latitude
  real :: lon_old   !< Previous longitude
  real :: uvel       !< Zonal velocity of particle (m/s)
  real :: vvel       !< Meridional velocity of particle (m/s)
  real :: uvel_old  !< Previous zonal velocity component (m/s)
  real :: vvel_old !< Previous meridional velocity component (m/s)
  integer :: year  !< Year of this record
  integer :: particle_num  !< Current particle number
  integer(kind=int64) :: id = -1 !< Particle Identifier
  type(xyt), pointer :: next=>null()  !< Pointer to the next position in the list
end type xyt

!>particle types are data structures describing a tracked particle
type, public :: particle
  type(particle), pointer :: prev=>null() !< Previous link in list
  type(particle), pointer :: next=>null() !< Next link in list
! State variables (specific to the particles, needed for restarts)
  real :: lon !< Longitude of particle (degree N or unit of grid coordinate)
  real :: lat !< Latitude of particle (degree E or unit of grid coordinate)
  real :: depth !< Depth of particle
  real :: uvel !< Zonal velocity of particle (m/s)
  real :: vvel !< Meridional velocity of particle (m/s)
  real :: lon_old !< previous lon (degrees)
  real :: lat_old !< previous lat (degrees)
  real :: uvel_old  !< previous uvel
  real :: vvel_old  !< previous vvel
  real :: start_lon !< starting longitude where particle was created
  real :: start_lat !< starting latitude where particle was created
  real :: start_day       !< origination position (degrees) and day
  integer :: start_year                         !< origination year
  real :: halo_part  !< equal to zero for particles on the computational domain, and 1 for particles on the halo
  integer(kind=int64) :: id                      !< particle identifier
  integer(kind=int64) :: drifter_num             !< particle identifier
  integer :: ine                           !< nearest i-index in NE direction (for convenience)
  integer :: jne                           !< nearest j-index in NE direction (for convenience)
  real :: xi                               !< non-dimensional x-coordinate within current cell (0..1)
  real :: yj                                !< non-dimensional y-coordinate within current cell (0..1)
  real :: uo                                !< zonal ocean velocity
  real :: vo                                !< meridional ocean velocity
                                                !< by the particle (m/s)
  type(xyt), pointer :: trajectory=>null() !< Trajectory for this particle
end type particle


!>A buffer structure for message passing
type, public :: buffer
  integer :: size=0 !< Size of buffer
  real, dimension(:,:), pointer :: data !< Buffer memory
end type buffer

!> A wrapper for the particle linked list (since an array of pointers is not allowed)
type, public :: linked_list
  type(particle), pointer :: first=>null() !< Pointer to the beginning of a linked list of parts
end type linked_list


!> A grand data structure for the particles in the local MOM domain
type, public :: particles !; private
  type(particles_gridded) :: grd !< Container with all gridded data
  type(linked_list), dimension(:,:), allocatable :: list !< Linked list of particles
  type(xyt), pointer :: trajectories=>null() !< A linked list for detached segments of trajectories
  real :: dt !< Time-step between particle calls
  integer :: current_year !< Current year (years)
  real :: current_yearday !< Current year-day, 1.00-365.99, (days)
  integer :: traj_sample_hrs !< Period between sampling for trajectories (hours)
  integer :: traj_write_hrs !< Period between writing of trajectories (hours)
  integer :: verbose_hrs !< Period between terminal status reports (hours)
  !>@{
  !! Handles for clocks
  integer :: clock, clock_mom, clock_the, clock_int, clock_cal, clock_com, clock_ini, clock_ior, clock_iow, clock_dia
  integer :: clock_trw, clock_trp
  !>@}
  logical :: restarted=.false. !< Indicate whether we read state from a restart or not
  logical :: Runge_not_Verlet=.True. !< True=Runge-Kutta, False=Verlet.
  logical :: ignore_missing_restart_parts=.False. !< True allows the model to ignore particles missing in the restart.
  logical :: halo_debugging=.False. !< Use for debugging halos (remove when its working)
  logical :: save_short_traj=.false. !< True saves only lon,lat,time,id in particle_trajectory.nc
  logical :: ignore_traj=.False. !< If true, then model does not write trajectory data at all
  logical :: use_new_predictive_corrective =.False. !< Flag to use Bob's predictive corrective particle scheme
  !Added by Alon
  integer(kind=int64) :: debug_particle_with_id = -1 !< If positive, monitors a part with this id
  type(buffer), pointer :: obuffer_n=>null() !< Buffer for outgoing parts to the north
  type(buffer), pointer :: ibuffer_n=>null() !< Buffer for incoming parts from the north
  type(buffer), pointer :: obuffer_s=>null() !< Buffer for outgoing parts to the south
  type(buffer), pointer :: ibuffer_s=>null() !< Buffer for incoming parts from the south
  type(buffer), pointer :: obuffer_e=>null() !< Buffer for outgoing parts to the east
  type(buffer), pointer :: ibuffer_e=>null() !< Buffer for incoming parts from the east
  type(buffer), pointer :: obuffer_w=>null() !< Buffer for outgoing parts to the west
  type(buffer), pointer :: ibuffer_w=>null() !< Buffer for incoming parts from the west
  type(buffer), pointer :: obuffer_io=>null() !< Buffer for outgoing parts during i/o
  type(buffer), pointer :: ibuffer_io=>null() !< Buffer for incoming parts during i/o
end type particles


end module particles_types_mod
