!> This module contains a set of data structures and interfaces for compiling the MOM6 DA
!! driver code.
module ocean_da_types_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains a set of data structures and interfaces for compiling the MOM6 DA
! driver code. This code is not yet finalized and will be replaced by supported
! software at some later date.
!
! 3/22/18
! matthew.harrison@noaa.gov
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef MAX_LEVS_FILE_
#define MAX_LEVS_FILE_ 50
#endif

#ifndef MAX_LINKS_
#define MAX_LINKS_ 100
#endif

!============================================================
! This module contains type declarations and default values
! for oda modules.
!============================================================

! Contact: Matthew.Harrison@noaa.gov and Feiyu.Lu@noaa.gov

  use time_manager_mod, only : time_type
  !use obs_tools_mod, only : obs_def_type
  !use mpp_domains_mod, only : domain2d

  implicit none

  private

  integer, parameter, public :: MAX_LEVELS_FILE = MAX_LEVS_FILE_ !< Controls record length for optimal storage
  integer, parameter, public :: MAX_LINKS = MAX_LINKS_ !< Maximum number of records per profile for storage for profiles
  integer, parameter, public :: UNKNOWN = 0

  integer, save, public :: TEMP_ID = 1
  integer, save, public :: SALT_ID = 2
  real, parameter, public :: MISSING_VALUE = -1.e10

  !> Type for ocean state in DA space (same decomposition and vertical grid)
  type, public :: OCEAN_CONTROL_STRUCT
     integer :: ensemble_size
     real, pointer, dimension(:,:,:) :: SSH=>NULL() !<sea surface height (m) across ensembles
     real, pointer, dimension(:,:,:,:) :: h=>NULL() !<layer thicknesses (m or kg) across ensembles
     real, pointer, dimension(:,:,:,:) :: T=>NULL() !<layer potential temperature (degC) across ensembles
     real, pointer, dimension(:,:,:,:) :: S=>NULL() !<layer salinity (psu or g kg-1) across ensembles
     real, pointer, dimension(:,:,:,:) :: U=>NULL() !<layer zonal velocity (m s-1) across ensembles
     real, pointer, dimension(:,:,:,:) :: V=>NULL() !<layer meridional velocity (m s-1) across ensembles
     integer, dimension(:), pointer :: id_t=>NULL(), id_s=>NULL() !< diagnostic IDs for temperature and salinity
     integer, dimension(:), pointer :: id_u=>NULL(), id_v=>NULL() !< diagnostic IDs for zonal and meridional velocity
     integer, dimension(:), pointer :: id_ssh=>NULL()  !< diagnostic IDs for SSH
  end type OCEAN_CONTROL_STRUCT

  !> Profile
  type, public :: ocean_profile_type
     integer :: variable !< variable ids are defined by the ocean_types module (e.g. TEMP_ID, SALT_ID)
     integer :: inst_type !< instrument types are defined by platform class
                          !! (e.g. MOORING, DROP, etc.) and instrument type (XBT, CDT, etc.)
     integer :: nvar !< number of observations types associated with the current profile
     real    :: project !< e.g. FGGE, COARE, ACCE, ...
     real    :: probe !< MBT, XBT, drifting buoy
     real    :: ref_inst !< instrument (thermograph, hull sensor, ...)
     integer :: wod_cast_num !< NODC world ocean dataset unique id
     real    :: fix_depth !< adjust profile depths (for XBT drop rate corrections)
     real    :: ocn_vehicle !< ocean vehicle type
     real    :: database_id !< a unique profile id
     integer :: levels !< number of levels in the current profile
     integer :: basin_mask !<1:Southern Ocean, 2:Atlantic Ocean, 3:Pacific Ocean,
                           !! 4:Arctic Ocean, 5:Indian Ocean, 6:Mediterranean Sea, 7:Black Sea,
                           !!  8:Hudson Bay, 9:Baltic Sea, 10:Red Sea, 11:Persian Gulf
     integer :: profile_flag !< an overall flag for the profile
     integer :: profile_flag_s !< an overall flag for the profile salinity
     real :: lat, lon !< latitude and longitude (degrees E and N)
     logical :: accepted !< logical flag to disable a profile
     integer :: nlinks !< number of links used to construct the profile (when reading from disk)
     type(ocean_profile_type), pointer :: next=>NULL() !< all profiles are stored as linked list.
     type(ocean_profile_type), pointer :: prev=>NULL()
     type(ocean_profile_type), pointer :: cnext=>NULL() ! current profiles are stored as linked list.
     type(ocean_profile_type), pointer :: cprev=>NULL()
     integer :: nbr_xi, nbr_yi ! nearest neighbor model gridpoint for the profile
     real :: nbr_dist ! distance to nearest neighbor model gridpoint
     real, dimension(:), pointer :: depth
     real, dimension(:), pointer :: data_t => NULL(), data_s => NULL()
     real, dimension(:), pointer :: data
     !integer, dimension(:), pointer :: flag_t
     !integer, dimension(:), pointer :: flag_s ! level-by-level flags for salinity
     !::sdu:: For now ECDA use flag as a logical, will likely change in future releases.
     logical, dimension(:), pointer :: flag
     real    :: temp_err, salt_err ! measurement error
     !real, dimension(:), pointer :: ms_t ! ms temperature by level
     !real, dimension(:), pointer :: ms_s ! ms salinity by level
     real, dimension(:), pointer :: ms_inv => NULL()
     real, dimension(:), pointer :: ms => NULL()
!     type(obs_def_type), dimension(:), pointer :: obs_def => NULL()
     type(time_type) :: time
     integer         :: yyyy
     integer         :: mmdd
     !type(time_type), pointer :: Model_time ! each profile can be associated
                                             ! with a first-guess field with an associated time and grid
     !type(grid_type), pointer :: Model_grid
     real :: i_index, j_index ! model longitude and latitude indices respectively
     real, dimension(:), pointer :: k_index ! model depth indices
     type(time_type) :: tdiff      ! positive difference between model time and observation time
  end type ocean_profile_type

  !> Grid information for ODA purposes, including arrays of
  !! lat, lon, depth, thickness, basin and land mask
  type, public :: grid_type
     real, pointer, dimension(:,:) :: x=>NULL(), y=>NULL()
     !real, pointer, dimension(:,:) :: x_bound=>NULL(), y_bound=>NULL()
     !real, pointer, dimension(:,:) :: dx=>NULL(), dy=>NULL()
     real, pointer, dimension(:,:,:) :: z=>NULL()
     real, pointer, dimension(:,:,:) :: h=>NULL()
     !real, pointer, dimension(:) :: z_bound=>NULL()
     !real, pointer, dimension(:) :: dz => NULL()
     real, pointer, dimension(:,:) :: basin_mask => NULL()
     real, pointer, dimension(:,:,:) :: mask => NULL()
     !type(domain2d), pointer :: Dom ! FMS domain type
     !logical :: cyclic
     integer :: ni, nj, nk
  end type grid_type

end module ocean_da_types_mod
