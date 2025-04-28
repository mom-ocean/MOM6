!> A set of dummy interfaces for compiling the MOM6 drifters code
module MOM_particles_mod

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_grid,         only : ocean_grid_type
use MOM_time_manager, only : time_type, get_date, operator(-)
use MOM_variables,    only : thermo_var_ptrs
use particles_types_mod, only : particles, particles_gridded

implicit none ; private

public particles, particles_run, particles_init, particles_save_restart, particles_end
public particles_to_k_space, particles_to_z_space

contains

!> Initializes particles container "parts"
subroutine particles_init(parts, Grid, Time, dt, u, v, h)
  ! Arguments
  type(particles), pointer, intent(out) :: parts !< Container for all types and memory
  type(ocean_grid_type), target, intent(in) :: Grid !< Grid type from parent model
  type(time_type), intent(in) :: Time !< Time type from parent model
  real, intent(in)            :: dt !< particle timestep in seconds [T ~> s]
  real, dimension(:,:,:),intent(in)      :: u !< Zonal velocity field [L T-1 ~> m s-1]
  real, dimension(:,:,:),intent(in)      :: v !< Meridional velocity field [L T-1 ~> m s-1]
  real, dimension(:,:,:),intent(in)      :: h !< Thickness of each layer [H ~> m or kg m-2]
end subroutine particles_init

!> The main driver the steps updates particles
subroutine particles_run(parts, time, uo, vo, ho, tv, dt_adv, use_uh)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  type(time_type), intent(in) :: time !< Model time
  real, dimension(:,:,:), intent(in) :: uo !< If use_uh is false, ocean zonal velocity [L T-1 ~>m s-1].
                                           !! If use_uh is true, accumulated zonal thickness fluxes
                                           !! that are used to advect tracers [H L2 ~> m3 or kg]
  real, dimension(:,:,:), intent(in) :: vo !< If use_uh is false, ocean meridional velocity [L T-1 ~>m s-1].
                                           !! If use_uh is true, accumulated meridional thickness fluxes
                                           !! that are used to advect tracers [H L2 ~> m3 or kg]
  real, dimension(:,:,:), intent(in) :: ho !< Ocean layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),  intent(in) :: tv !< structure containing pointers to available thermodynamic fields
  real, intent(in) :: dt_adv !< timestep for advecting particles [s]
  logical :: use_uh !< Flag for whether u and v are weighted by thickness

end subroutine particles_run


!>Save particle locations (and sometimes other vars) to restart file
subroutine particles_save_restart(parts, h, directory, time, time_stamped)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  real, dimension(:,:,:),intent(in)      :: h !< Thickness of each layer [H ~> m or kg m-2]
  character(len=*), intent(in) :: directory !< The directory where the restart files are to be written
  type(time_type), intent(in) :: time !< The current model time
  logical, optional, intent(in) :: time_stamped !< If present and true, add time-stamp to the restart file names

end subroutine particles_save_restart

!> Deallocate all memory and disassociated pointer
subroutine particles_end(parts, h, temp, salt)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  real, dimension(:,:,:),intent(in)      :: h !< Thickness of each layer [H ~> m or kg m-2]
  real, dimension(:,:,:), optional, intent(in) :: temp !< Optional container for temperature [C ~> degC]
  real, dimension(:,:,:), optional, intent(in) :: salt !< Optional container for salinity [S ~> ppt]

end subroutine particles_end

subroutine particles_to_k_space(parts, h)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  real, dimension(:,:,:),intent(in)      :: h !< Thickness of layers [H ~> m or kg m-2]

end subroutine particles_to_k_space


subroutine particles_to_z_space(parts, h)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  real, dimension(:,:,:),intent(in)      :: h !< Thickness of layers [H ~> m or kg m-2]

end subroutine particles_to_z_space

end module MOM_particles_mod
