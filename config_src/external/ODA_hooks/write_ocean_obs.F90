!> Dummy interfaces for writing ODA data
module write_ocean_obs_mod


 use ocean_da_types_mod, only : ocean_profile_type
 use MOM_time_manager, only : time_type, get_time, set_date

 implicit none

 private

 public :: open_profile_file, write_profile, close_profile_file, &
            write_ocean_obs_init

contains

!> Open a profile file
integer function open_profile_file(name, nvar, grid_lon, grid_lat,thread,fset)
  character(len=*), intent(in) :: name !< File name
  integer, intent(in), optional :: nvar !< Number of variables
  real, dimension(:), optional, intent(in) :: grid_lon !< Longitude [degreeE]
  real, dimension(:), optional, intent(in) :: grid_lat !< Latitude [degreeN]
  integer, intent(in), optional :: thread !< Thread
  integer, intent(in), optional :: fset !< File set

  open_profile_file=-1
end function open_profile_file

!> Write a profile
subroutine write_profile(unit,profile)
  integer, intent(in) :: unit !< File unit
  type(ocean_profile_type), intent(in) :: profile !< Profile

  return
end subroutine write_profile

!> Close a profile file
subroutine close_profile_file(unit)
  integer, intent(in) :: unit !< File unit

  return
end subroutine close_profile_file

!> Initialize write_ocean_obs module
subroutine write_ocean_obs_init()

  return
end subroutine write_ocean_obs_init

end module write_ocean_obs_mod
