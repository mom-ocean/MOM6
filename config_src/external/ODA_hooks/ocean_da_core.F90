!> A set of dummy interfaces for compiling the MOM6 DA driver code.
module ocean_da_core_mod
  ! MOM modules
  use MOM_domains, only : MOM_domain_type, domain2D
  use MOM_time_manager, only : time_type, set_time, get_date
  ! ODA_tools modules
  use ocean_da_types_mod, only : ocean_profile_type, grid_type
  use kdtree, only : kd_root

  implicit none
  private
  public :: ocean_da_core_init
  public :: get_profiles

contains

  !> Initializes the MOM6 DA driver code.
  subroutine ocean_da_core_init(Domain, global_grid, Profiles, model_time)
    type(domain2D), pointer, intent(in) :: Domain !< A MOM domain type
    type(grid_type), pointer, intent(in) :: global_grid !< The global ODA horizontal grid type
    type(ocean_profile_type), pointer :: Profiles  !< This is an unstructured recursive list of profiles
                                          !! which are either within the localized domain corresponding
                                          !! to the Domain argument, or the global profile list (type).
    type(time_type), intent(in) :: model_time !< The current model time type.



    Profiles=>NULL()
    return
  end subroutine ocean_da_core_init


  !> Get profiles obs within the current analysis interval
  subroutine get_profiles(model_time, Profiles, Current_profiles)
    type(time_type), intent(in) :: model_time !< The current analysis time.
    type(ocean_profile_type), pointer :: Profiles !< The full recursive list of profiles.
    type(ocean_profile_type), pointer :: Current_profiles !< A returned list of profiles for the
                                                          !! current analysis step.

    Profiles=>NULL()
    Current_Profiles=>NULL()

    return
  end subroutine get_profiles


end module ocean_da_core_mod
