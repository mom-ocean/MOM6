!> Dummy data structures and methods for drifters package
module particles_types_mod


!>xyt is a data structure containing particle position and velocity fields.
type :: xyt
  real :: lon, lat, day      !< Current position (degrees) and day
  real :: lat_old, lon_old   !< Previous position (degrees)
  real :: uvel, vvel         !< Current velocity components (m/s)
  real :: uvel_old, vvel_old !< Previous velocity components (m/s)
  integer :: year, particle_num  !< Current year and particle number
  integer(kind=8) :: id = -1 !< Particle Identifier
  type(xyt), pointer :: next=>null()  !< Pointer to the next position in the list
end type xyt


!>A buffer structure for message passing
type :: buffer
  integer :: size=0
  real, dimension(:,:), pointer :: data
end type buffer

!> A wrapper for the particle linked list (since an array of pointers is not allowed)
type :: linked_list
  type(particle), pointer :: first=>null() !< Pointer to the beginning of a linked list of parts
end type linked_list


type :: particles !; private
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
  logical :: use_new_predictive_corrective =.False. !< Flag to use Bob's predictive corrective particle scheme- Added by Alon
  integer(kind=8) :: debug_particle_with_id = -1 !< If positive, monitors a part with this id
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
