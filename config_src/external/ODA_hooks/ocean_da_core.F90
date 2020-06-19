!> A set of dummy interfaces for compiling the MOM6 DA driver code.
module ocean_da_core_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains a set of dummy interfaces for compiling the MOM6 DA
! driver code. These interfaces are not finalized and will be replaced by supported
! interfaces at some later date.
!
! 3/22/18
! matthew.harrison@noaa.gov
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use mpp_domains_mod, only : domain2d
  use time_manager_mod, only : time_type, set_time, get_date
  ! ODA_tools modules
  use ocean_da_types_mod, only : ocean_profile_type, grid_type


  implicit none
  private
  public :: ocean_da_core_init, open_profile_dataset
  public :: get_profiles, copy_profiles

contains

  !> Initialize ODA
  subroutine ocean_da_core_init(Domain, T_grid, Profiles, model_time)
    type(domain2d), pointer, intent(in) :: Domain !< MOM type for the local domain`
    type(grid_type), pointer, intent(in) :: T_grid !< MOM grid type for the local domain
    type(ocean_profile_type), pointer :: Profiles
    !< This is an unstructured recursive list of profiles
    !< which are either within the localized domain corresponding
    !< to the Domain argument, or the global profile list
    type(time_type), intent(in) :: model_time !< Model time

    Profiles=>NULL()
    return
  end subroutine ocean_da_core_init

  !> Open a profile dataset
  subroutine open_profile_dataset(Profiles, Domain, T_grid, &
                  filename, time_start, time_end, obs_variable, localize)
    type(ocean_profile_type), pointer :: Profiles
    !< This is an unstructured recursive list of profiles
    !< which are either within the localized domain corresponding
    !< to the Domain argument, or the global profile list
    type(domain2d), pointer, intent(in) :: Domain !< MOM type for the local domain
    type(grid_type), pointer, intent(in) :: T_grid !< MOM grid type for the local domain
    character(len=*), intent(in) :: filename !< filename containing profile data
    type(time_type), intent(in) :: time_start, time_end !< start and end times for the analysis
    integer, intent(in), optional :: obs_variable !< If present, then extract corresponding data
    !< from file, otherwise, extract all available data which.
    logical, intent(in), optional :: localize !< Localize the observations to the current computational domain

    return

  end subroutine open_profile_dataset

  !> Get profiles obs relevant to current analysis interval
  subroutine get_profiles(model_time, Profiles, Current_profiles)
    type(time_type), intent(in) :: model_time
    type(ocean_profile_type), pointer :: Profiles
    type(ocean_profile_type), pointer :: Current_profiles

    Profiles=>NULL()
    Current_Profiles=>NULL()

    return
  end subroutine get_profiles

  !> Copy profiles at current analysis step from a linked list to an array
  !! feasible now since the number of localized current profiles is small
  subroutine copy_profiles(Current_profiles, Profiles)
    type(ocean_profile_type), pointer :: Current_profiles
    type(ocean_profile_type), pointer, dimension(:) :: Profiles

    return

  end subroutine copy_profiles

  !> Copy observations
  subroutine copy_obs(obs_in, obs_out)
    type(ocean_profile_type), pointer :: obs_in
    type(ocean_profile_type), pointer :: obs_out

    return
  end subroutine copy_obs

end module ocean_da_core_mod
