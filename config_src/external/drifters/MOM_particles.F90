!> A set of dummy interfaces for compiling the MOM6 drifters code
module MOM_particles_mod

use MOM_time_manager, only : time_type, get_date, operator(-)
use MOM_variables,    only : thermo_var_ptrs


use particles_types_mod, only: particles, particles_gridded

public particles_run

contains

subroutine particles_init(parts, Grid, Time, dt, u, v)
 type(particles), pointer, intent(out) :: parts
 type(ocean_grid_type), target, intent(in) :: Grid !< Grid type from parent model
 type(time_type), intent(in) :: Time !< Time type from parent model
 real, intent(in)            :: dt !< particle timestep in seconds
 real, dimension(:,:,:),intent(in)      :: u, v !< Horizontal velocity fields

end subroutine particles_init


subroutine particles_run(parts, time, uo, vo, ho, tv, stagger)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  type(time_type), intent(in) :: time !< Model time
  real, dimension(:,:,:),intent(in) :: uo !< Ocean zonal velocity (m/s)
  real, dimension(:,:,:),intent(in) :: vo !< Ocean meridional velocity (m/s)
  real, dimension(:,:,:),intent(in) :: ho !< Ocean layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs), intent(in) :: tv !< structure containing pointers to available thermodynamic fields
  integer, optional, intent(in) :: stagger

end subroutine particles_run

subroutine particles_save_restart(parts,temp,salt)
! Arguments
type(particles), pointer :: parts
real,dimension(:,:,:),optional,intent(in) :: temp, salt

end subroutine particles_save_restart

end module MOM_particles_mod
