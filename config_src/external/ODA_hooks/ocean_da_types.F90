!> Dummy aata structures and methods for ocean data assimilation.
module ocean_da_types_mod

  use MOM_time_manager, only : time_type

  implicit none

  private

  !> Example type for ocean ensemble DA state
  type, public :: OCEAN_CONTROL_STRUCT
     integer :: ensemble_size        !< ensemble size
     real, pointer, dimension(:,:,:)   :: SSH=>NULL() !<sea surface height (m) across ensembles
     real, pointer, dimension(:,:,:,:) :: h=>NULL() !<layer thicknesses (m or kg) across ensembles
     real, pointer, dimension(:,:,:,:) :: T=>NULL() !<layer potential temperature (degC) across ensembles
     real, pointer, dimension(:,:,:,:) :: S=>NULL() !<layer salinity (psu or g kg-1) across ensembles
     real, pointer, dimension(:,:,:,:) :: U=>NULL() !<layer zonal velocity (m s-1) across ensembles
     real, pointer, dimension(:,:,:,:) :: V=>NULL() !<layer meridional velocity (m s-1) across ensembles
  end type OCEAN_CONTROL_STRUCT

  !> Example of a profile type
  type, public :: ocean_profile_type
     integer :: inst_type !< A numeric code indicating the type of instrument (e.g. ARGO drifter, CTD, ...)
     logical :: initialized !< a True value indicates that this profile has been allocated for use
     logical :: colocated !< a True value indicated that the measurements of (num_variables) data are
                          !! co-located in space-time
     integer :: ensemble_size !< size of the ensemble of model states used in association with this profile
     integer :: num_variables !< number of measurement types associated with this profile.
     integer, pointer, dimension(:) :: var_id !< variable ids are defined by the ocean_types module
     integer :: platform !< platform types are defined by platform class (e.g. MOORING, DROP, etc.)
                         !! and instrument type (XBT, CDT, etc.)
     integer :: levels !< number of levels in the current profile
     integer :: basin_mask !< 1:Southern Ocean, 2:Atlantic Ocean, 3:Pacific Ocean,
                           !! 4:Arctic Ocean, 5:Indian Ocean, 6:Mediterranean Sea, 7:Black Sea,
                           !! 8:Hudson Bay, 9:Baltic Sea, 10:Red Sea, 11:Persian Gulf
     integer :: profile_flag !< an overall flag for the profile
     real :: lat      !< latitude [degrees_N]
     real :: lon      !< longitude [degrees_E]
     logical :: accepted !< logical flag to disable a profile
     type(time_type) :: time_window !< The time window associated with this profile [s]
     real, pointer, dimension(:) :: obs_error  !< The observation error by variable
     real  :: loc_dist   !< The impact radius of this observation (m)
     type(ocean_profile_type), pointer :: next=>NULL() !< all profiles are stored as linked list.
     type(ocean_profile_type), pointer :: prev=>NULL() !< previous
     type(ocean_profile_type), pointer :: cnext=>NULL() !< current profiles are stored as linked list.
     type(ocean_profile_type), pointer :: cprev=>NULL() !< previous
     integer :: nbr_xi !< x nearest neighbor model gridpoint for the profile
     integer :: nbr_yi !< y nearest neighbor model gridpoint for the profile
     real :: nbr_dist !< distance to nearest neighbor model gridpoint
     logical :: compute !< profile is within current compute domain
     real, dimension(:,:), pointer :: depth => NULL() !< depth of measurement [m]
     real, dimension(:,:), pointer :: data => NULL() !< data by variable type
     integer, dimension(:,:), pointer :: flag => NULL() !< flag by depth and variable type
     real, dimension(:,:,:), pointer :: forecast => NULL() !< ensemble member first guess
     real, dimension(:,:,:), pointer :: analysis => NULL() !< ensemble member analysis
     type(forward_operator_type), pointer :: obs_def => NULL() !< observation forward operator
     type(time_type) :: time !< profile time type
     real :: i_index !< model longitude indices respectively
     real :: j_index !< model latitude indices respectively
     real, dimension(:,:), pointer :: k_index !< model depth indices
     type(time_type) :: tdiff !< difference between model time and observation time
     character(len=128) :: filename !< a filename
  end type ocean_profile_type

  !>  Example forward operator type.
  type, public :: forward_operator_type
     integer :: num                      !< how many?
     integer, dimension(2) :: state_size !< for
     integer, dimension(:), pointer :: state_var_index !< for flattened data
     integer, dimension(:), pointer :: i_index !< i-dimension index
     integer, dimension(:), pointer :: j_index !< j-dimension index
     real, dimension(:), pointer :: coef !< coefficient
  end type forward_operator_type

  !> Grid type for DA
  type, public :: grid_type
     real, pointer, dimension(:,:) :: x=>NULL()    !< x
     real, pointer, dimension(:,:) :: y=>NULL()    !< y
     real, pointer, dimension(:,:,:) :: z=>NULL()  !< z
     real, pointer, dimension(:,:,:) :: h=>NULL()  !< h
     real, pointer, dimension(:,:) :: basin_mask => NULL() !< basin mask
     real, pointer, dimension(:,:,:) :: mask => NULL()     !< land mask?
     real, pointer, dimension(:,:) :: bathyT => NULL()     !< bathymetry at T points
     logical :: tripolar_N    !< True for tripolar grids
     integer :: ni !< ni
     integer :: nj !< nj
     integer :: nk !< nk
  end type grid_type

end module ocean_da_types_mod
