module coupler_types_mod

! This file is part of MOM6. See LICENSE.md for the license.

!   This module contains the coupler-type declarations and methods for use in
! ocean-only configurations of MOM6.  It is intended that the version of
! coupler_types_mod that is avialable from FMS will conform to this version with
! the FMS city release after warsaw.

use fms_io_mod,        only: restart_file_type, register_restart_field
use fms_io_mod,        only: query_initialized, restore_state
use time_manager_mod,  only: time_type
use diag_manager_mod,  only: register_diag_field, send_data
use data_override_mod, only: data_override
use mpp_domains_mod,   only: domain2D, mpp_redistribute
use mpp_mod,           only: stdout, mpp_error, FATAL, mpp_chksum

implicit none ; private

public coupler_type_copy, coupler_type_spawn, coupler_type_set_diags
public coupler_type_write_chksums, coupler_type_send_data, coupler_type_data_override
public coupler_type_register_restarts, coupler_type_restore_state
public coupler_type_increment_data, coupler_type_rescale_data
public coupler_type_copy_data, coupler_type_redistribute_data
public coupler_type_destructor, coupler_type_initialized
public coupler_type_extract_data, coupler_type_set_data

public coupler_type_copy_1d_2d
public coupler_type_copy_1d_3d

!
!       3-d fields
!
type, public :: coupler_3d_values_type
  character(len=48)       :: name = ' '  !< The diagnostic name for this array
  real, pointer, contiguous, dimension(:,:,:) :: values => NULL() !< The pointer to the
                                         !! array of values for this field; this
                                         !! should be changed to allocatable
  logical                 :: mean = .true. !< mean
  logical                 :: override = .false. !< override
  integer                 :: id_diag = 0 !< The diagnostic id for this array
  character(len=128)      :: long_name = ' ' !< The diagnostic long_name for this array
  character(len=128)      :: units = ' ' !< The units for this array
  integer                 :: id_rest = 0 !< The id of this array in the restart field
  logical                 :: may_init = .true. !< If true, there is an internal method
                                         !! that can be used to initialize this field
                                         !! if it can not be read from a restart file
end type coupler_3d_values_type

type, public :: coupler_3d_field_type
  character(len=48)                 :: name = ' ' !< name
  integer                           :: num_fields = 0 !< num_fields
  type(coupler_3d_values_type), pointer, dimension(:) :: field => NULL() !< field
  character(len=128)                :: flux_type = ' ' !< flux_type
  character(len=128)                :: implementation = ' ' !< implementation
  real, pointer, dimension(:)       :: param => NULL() !< param
  logical, pointer, dimension(:)    :: flag => NULL() !< flag
  integer                           :: atm_tr_index = 0 !< atm_tr_index
  character(len=128)                :: ice_restart_file = ' ' !< ice_restart_file
  character(len=128)                :: ocean_restart_file = ' ' !< ocean_restart_file
  type(restart_file_type), pointer  :: rest_type => NULL() !< A pointer to the restart_file_type
                                                           !! that is used for this field.
  logical                           :: use_atm_pressure !< use_atm_pressure
  logical                           :: use_10m_wind_speed !< use_10m_wind_speed
  logical                           :: pass_through_ice !< pass_through_ice
  real                              :: mol_wt = 0.0 !< mol_wt
end type coupler_3d_field_type

type, public :: coupler_3d_bc_type
  integer                                            :: num_bcs = 0  !< The number of boundary condition fields
  type(coupler_3d_field_type), dimension(:), pointer :: bc => NULL() !< A pointer to the array of boundary condition fields
  logical    :: set = .false.       !< If true, this type has been initialized
  integer    :: isd, isc, iec, ied  !< The i-direction data and computational domain index ranges for this type
  integer    :: jsd, jsc, jec, jed  !< The j-direction data and computational domain index ranges for this type
  integer    :: ks, ke              !< The k-direction index ranges for this type
end type coupler_3d_bc_type

!
!       2-d fields
!
type, public    :: coupler_2d_values_type
  character(len=48)       :: name = ' '  !< The diagnostic name for this array
  real, pointer, contiguous, dimension(:,:) :: values => NULL() !< The pointer to the
                                         !! array of values for this field; this
                                         !! should be changed to allocatable
  logical                 :: mean = .true. !< mean
  logical                 :: override = .false. !< override
  integer                 :: id_diag = 0 !< The diagnostic id for this array
  character(len=128)      :: long_name = ' ' !< The diagnostic long_name for this array
  character(len=128)      :: units = ' ' !< The units for this array
  integer                 :: id_rest = 0 !< The id of this array in the restart field
  logical                 :: may_init = .true. !< If true, there is an internal method
                                         !! that can be used to initialize this field
                                         !! if it can not be read from a restart file
end type coupler_2d_values_type

type, public    :: coupler_2d_field_type
  character(len=48)                 :: name = ' ' !< name
  integer                           :: num_fields = 0 !< num_fields
  type(coupler_2d_values_type), pointer, dimension(:)   :: field => NULL() !< field
  character(len=128)                :: flux_type = ' ' !< flux_type
  character(len=128)                :: implementation = ' ' !< implementation
  real, pointer, dimension(:)       :: param => NULL() !< param
  logical, pointer, dimension(:)    :: flag => NULL() !< flag
  integer                           :: atm_tr_index = 0 !< atm_tr_index
  character(len=128)                :: ice_restart_file = ' ' !< ice_restart_file
  character(len=128)                :: ocean_restart_file = ' ' !< ocean_restart_file
  type(restart_file_type), pointer  :: rest_type => NULL() !< A pointer to the restart_file_type
                                                          !! that is used for this field.
  logical                           :: use_atm_pressure !< use_atm_pressure
  logical                           :: use_10m_wind_speed !< use_10m_wind_speed
  logical                           :: pass_through_ice !< pass_through_ice
  real                              :: mol_wt = 0.0 !< mol_wt
end type coupler_2d_field_type

type, public    :: coupler_2d_bc_type
  integer                                            :: num_bcs = 0  !< The number of boundary condition fields
  type(coupler_2d_field_type), dimension(:), pointer :: bc => NULL() !< A pointer to the array of boundary condition fields
  logical    :: set = .false.       !< If true, this type has been initialized
  integer    :: isd, isc, iec, ied  !< The i-direction data and computational domain index ranges for this type
  integer    :: jsd, jsc, jec, jed  !< The j-direction data and computational domain index ranges for this type
end type coupler_2d_bc_type

!
!       1-d fields
!
type, public    :: coupler_1d_values_type
  character(len=48)           :: name = ' '  !< The diagnostic name for this array
  real, pointer, dimension(:) :: values => NULL() !< The pointer to the array of values
  logical                     :: mean = .true. !< mean
  logical                     :: override = .false. !< override
  integer                     :: id_diag = 0 !< The diagnostic id for this array
  character(len=128)          :: long_name = ' ' !< The diagnostic long_name for this array
  character(len=128)          :: units = ' ' !< The units for this array
  logical                     :: may_init = .true. !< If true, there is an internal method
                                             !! that can be used to initialize this field
                                             !! if it can not be read from a restart file
end type coupler_1d_values_type

type, public    :: coupler_1d_field_type
  character(len=48)              :: name = ' ' !< name
  integer                        :: num_fields = 0 !< num_fields
  type(coupler_1d_values_type), pointer, dimension(:)   :: field => NULL() !< field
  character(len=128)             :: flux_type = ' ' !< flux_type
  character(len=128)             :: implementation = ' ' !< implementation
  real, pointer, dimension(:)    :: param => NULL() !< param
  logical, pointer, dimension(:) :: flag => NULL() !< flag
  integer                        :: atm_tr_index = 0 !< atm_tr_index
  character(len=128)             :: ice_restart_file = ' ' !< ice_restart_file
  character(len=128)             :: ocean_restart_file = ' ' !< ocean_restart_file
  logical                        :: use_atm_pressure !< use_atm_pressure
  logical                        :: use_10m_wind_speed !< use_10m_wind_speed
  logical                        :: pass_through_ice !< pass_through_ice
  real                           :: mol_wt = 0.0 !< mol_wt
end type coupler_1d_field_type

type, public    :: coupler_1d_bc_type
  integer                                            :: num_bcs = 0  !< The number of boundary condition fields
  type(coupler_1d_field_type), dimension(:), pointer :: bc => NULL() !< A pointer to the array of boundary condition fields
  logical    :: set = .false.       !< If true, this type has been initialized
end type coupler_1d_bc_type

!----------------------------------------------------------------------
!   The following public parameters can help in selecting the sub-elements of a
! coupler type.  There are duplicate values because different boundary
! conditions have different sub-elements.
integer, parameter, public :: ind_pcair = 1 !< The index of the atmospheric concentration
integer, parameter, public :: ind_u10 = 2   !< The index of the 10 m wind speed
integer, parameter, public :: ind_psurf = 3 !< The index of the surface atmospheric pressure
integer, parameter, public :: ind_alpha = 1 !< The index of the solubility array for a tracer
integer, parameter, public :: ind_csurf = 2 !< The index of the ocean surface concentration
integer, parameter, public :: ind_sc_no = 3 !< The index for the Schmidt number for a tracer flux
integer, parameter, public :: ind_flux = 1  !< The index for the tracer flux
integer, parameter, public :: ind_deltap= 2 !< The index for ocean-air gas partial pressure change
integer, parameter, public :: ind_kw = 3    !< The index for the piston velocity
integer, parameter, public :: ind_deposition = 1 !< The index for the atmospheric deposition flux
integer, parameter, public :: ind_runoff = 1 !< The index for a runoff flux

!----------------------------------------------------------------------
!        Interface definitions for overloaded routines
!----------------------------------------------------------------------

!> This is the interface to spawn one coupler_bc_type into another and then
!! register diagnostics associated with the new type.
interface  coupler_type_copy
  module procedure coupler_type_copy_1d_2d, coupler_type_copy_1d_3d
  module procedure coupler_type_copy_2d_2d, coupler_type_copy_2d_3d
  module procedure coupler_type_copy_3d_2d, coupler_type_copy_3d_3d
end interface coupler_type_copy

!> This is the interface to spawn one coupler_bc_type into another.
interface  coupler_type_spawn
  module procedure CT_spawn_1d_2d, CT_spawn_2d_2d, CT_spawn_3d_2d
  module procedure CT_spawn_1d_3d, CT_spawn_2d_3d, CT_spawn_3d_3d
end interface coupler_type_spawn

!> This is the interface to copy the field data from one coupler_bc_type
!! to another of the same rank, size and decomposition.
interface coupler_type_copy_data
  module procedure CT_copy_data_2d, CT_copy_data_3d, CT_copy_data_2d_3d
end interface coupler_type_copy_data

!> This is the interface to redistribute the field data from one coupler_bc_type
!! to another of the same rank and global size, but a different decomposition.
interface coupler_type_redistribute_data
  module procedure CT_redistribute_data_2d, CT_redistribute_data_3d
end interface coupler_type_redistribute_data

!> This is the interface to rescale the field data in a coupler_bc_type.
interface coupler_type_rescale_data
  module procedure CT_rescale_data_2d, CT_rescale_data_3d
end interface coupler_type_rescale_data

!> This is the interface to increment the field data from one coupler_bc_type
!! with the data from another.  Both must have the same horizontal size and
!! decomposition, but a 2d type may be incremented by a 2d or 3d type
interface coupler_type_increment_data
  module procedure CT_increment_data_2d_2d, CT_increment_data_3d_3d, CT_increment_data_2d_3d
end interface coupler_type_increment_data

!> This is the interface to extract a field in a coupler_bc_type into an array.
interface coupler_type_extract_data
  module procedure CT_extract_data_2d, CT_extract_data_3d, CT_extract_data_3d_2d
end interface coupler_type_extract_data

!> This is the interface to set a field in a coupler_bc_type from an array.
interface coupler_type_set_data
  module procedure CT_set_data_2d, CT_set_data_3d, CT_set_data_2d_3d
end interface coupler_type_set_data

!> This is the interface to set diagnostics for the arrays in a coupler_bc_type.
interface coupler_type_set_diags
  module procedure CT_set_diags_2d, CT_set_diags_3d
end interface coupler_type_set_diags

!> This is the interface to write out checksums for the elements of a coupler_bc_type.
interface coupler_type_write_chksums
  module procedure CT_write_chksums_2d, CT_write_chksums_3d
end interface coupler_type_write_chksums

!> This is the interface to write out diagnostics of the arrays in a coupler_bc_type.
interface coupler_type_send_data
  module procedure CT_send_data_2d, CT_send_data_3d
end interface coupler_type_send_data

!> This is the interface to override the values of the arrays in a coupler_bc_type.
interface coupler_type_data_override
  module procedure CT_data_override_2d, CT_data_override_3d
end interface coupler_type_data_override

!> This is the interface to register the fields in a coupler_bc_type to be saved
!! in restart files.
interface coupler_type_register_restarts
  module procedure CT_register_restarts_2d, CT_register_restarts_3d
  module procedure CT_register_restarts_to_file_2d, CT_register_restarts_to_file_3d
end interface coupler_type_register_restarts

!> This is the interface to read in the fields in a coupler_bc_type that have
!! been saved in restart files.
interface coupler_type_restore_state
  module procedure CT_restore_state_2d, CT_restore_state_3d
end interface coupler_type_restore_state

!> This function interface indicates whether a coupler_bc_type has been initialized.
interface coupler_type_initialized
  module procedure CT_initialized_1d, CT_initialized_2d, CT_initialized_3d
end interface coupler_type_initialized

!> This is the interface to deallocate any data associated with a coupler_bc_type.
interface coupler_type_destructor
  module procedure CT_destructor_1d, CT_destructor_2d, CT_destructor_3d
end interface coupler_type_destructor

contains

!#######################################################################
!> \brief Copy fields from one coupler type to another. 1-D to 2-D version for generic coupler_type_copy.
!!
!! Template:
!!
!! ~~~~~~~~~~{.f90}
!!   call coupler_type_copy(var_in, var_out, is, ie, js, je, &
!!        diag_name, axes, time, suffix = 'something')
!! ~~~~~~~~~~
subroutine coupler_type_copy_1d_2d(var_in, var_out, is, ie, js, je,     &
     diag_name, axes, time, suffix)

  type(coupler_1d_bc_type), intent(in)    :: var_in !< variable to copy information from
  type(coupler_2d_bc_type), intent(inout) :: var_out !< variable to copy information to
  integer, intent(in)                     :: is !< lower bound of first dimension
  integer, intent(in)                     :: ie !< upper bound of first dimension
  integer, intent(in)                     :: js !< lower bound of second dimension
  integer, intent(in)                     :: je !< upper bound of second dimension
  character(len=*), intent(in)            :: diag_name !< name for diagnostic file--if blank, then don't register the fields
  integer, dimension(:), intent(in)       :: axes !< array of axes identifiers for diagnostic variable registration
  type(time_type), intent(in)             :: time !< model time variable for registering diagnostic field
  character(len=*), intent(in), optional  :: suffix !< optional suffix to make the name identifier unique

  character(len=256), parameter :: error_header = &
       '==>Error from coupler_types_mod (coupler_type_copy_1d_2d):'
  character(len=400)      :: error_msg
  integer                 :: m, n

  if (var_out%num_bcs > 0) then
    ! It is an error if the number of output fields exceeds zero, because it means this
    ! type has already been populated.
    call mpp_error(FATAL, trim(error_header) // ' Number of output fields exceeds zero')
  endif

  if (var_in%num_bcs >= 0) &
    call CT_spawn_1d_2d(var_in, var_out, (/ is, is, ie, ie /), (/ js, js, je, je /), suffix)

  if ((var_out%num_bcs > 0) .and. (diag_name .ne. ' ')) &
    call CT_set_diags_2d(var_out, diag_name, axes, time)

end subroutine  coupler_type_copy_1d_2d

!#######################################################################
!> \brief Copy fields from one coupler type to another. 1-D to 3-D version for generic coupler_type_copy.
!!
!! Template:
!!
!! ~~~~~~~~~~{.f90}
!!   call coupler_type_copy(var_in, var_out, is, ie, js, je, kd, &
!!        diag_name, axes, time, suffix = 'something')
!! ~~~~~~~~~~
!!
!! \throw FATAL, "Number of output fields is non-zero"
!! \throw FATAL, "var_out%bc already associated"
!! \throw FATAL, "var_out%bc([n])%field already associated"
!! \throw FATAL, "var_out%bc([n])%field([m])%values already associated"
!! \throw FATAL, "axes less than 3 elements"
subroutine coupler_type_copy_1d_3d(var_in, var_out, is, ie, js, je, kd, &
     diag_name, axes, time, suffix)

  type(coupler_1d_bc_type), intent(in)    :: var_in !< variable to copy information from
  type(coupler_3d_bc_type), intent(inout) :: var_out !< variable to copy information to
  integer, intent(in)                     :: is !< lower bound of first dimension
  integer, intent(in)                     :: ie !< upper bound of first dimension
  integer, intent(in)                     :: js !< lower bound of second dimension
  integer, intent(in)                     :: je !< upper bound of second dimension
  integer, intent(in)                     :: kd !< third dimension
  character(len=*), intent(in)            :: diag_name !< name for diagnostic file--if blank, then don't register the fields
  integer, dimension(:), intent(in)       :: axes !< array of axes identifiers for diagnostic variable registration
  type(time_type), intent(in)             :: time !< model time variable for registering diagnostic field
  character(len=*), intent(in), optional  :: suffix !< optional suffix to make the name identifier unique

  character(len=256), parameter :: error_header = &
     '==>Error from coupler_types_mod (coupler_type_copy_1d_3d):'
  character(len=400)      :: error_msg
  integer                 :: m, n


  if (var_out%num_bcs > 0) then
    ! It is an error if the number of output fields exceeds zero, because it means this
    ! type has already been populated.
    call mpp_error(FATAL, trim(error_header) // ' Number of output fields exceeds zero')
  endif

  if (var_in%num_bcs >= 0) &
    call CT_spawn_1d_3d(var_in, var_out,  (/ is, is, ie, ie /), (/ js, js, je, je /), (/1, kd/), suffix)

  if ((var_out%num_bcs > 0) .and. (diag_name .ne. ' ')) &
    call CT_set_diags_3d(var_out, diag_name, axes, time)

end subroutine  coupler_type_copy_1d_3d

!#######################################################################
!> \brief Copy fields from one coupler type to another. 2-D to 2-D version for generic coupler_type_copy.
!!
!! Template:
!!
!! ~~~~~~~~~~{.f90}
!!   call coupler_type_copy(var_in, var_out, is, ie, js, je, &
!!        diag_name, axes, time, suffix = 'something')
!! ~~~~~~~~~~
subroutine coupler_type_copy_2d_2d(var_in, var_out, is, ie, js, je,     &
     diag_name, axes, time, suffix)

  type(coupler_2d_bc_type), intent(in)    :: var_in !< variable to copy information from
  type(coupler_2d_bc_type), intent(inout) :: var_out !< variable to copy information to
  integer, intent(in)                     :: is !< lower bound of first dimension
  integer, intent(in)                     :: ie !< upper bound of first dimension
  integer, intent(in)                     :: js !< lower bound of second dimension
  integer, intent(in)                     :: je !< upper bound of second dimension
  character(len=*), intent(in)            :: diag_name !< name for diagnostic file--if blank, then don't register the fields
  integer, dimension(:), intent(in)       :: axes !< array of axes identifiers for diagnostic variable registration
  type(time_type), intent(in)             :: time !< model time variable for registering diagnostic field
  character(len=*), intent(in), optional  :: suffix !< optional suffix to make the name identifier unique

  character(len=256), parameter :: error_header = &
       '==>Error from coupler_types_mod (coupler_type_copy_2d_2d):'
  character(len=400)      :: error_msg
  integer                 :: m, n

  if (var_out%num_bcs > 0) then
    ! It is an error if the number of output fields exceeds zero, because it means this
    ! type has already been populated.
    call mpp_error(FATAL, trim(error_header) // ' Number of output fields exceeds zero')
  endif

  if (var_in%num_bcs >= 0) &
    call CT_spawn_2d_2d(var_in, var_out, (/ is, is, ie, ie /), (/ js, js, je, je /), suffix)

  if ((var_out%num_bcs > 0) .and. (diag_name .ne. ' ')) &
    call CT_set_diags_2d(var_out, diag_name, axes, time)

end subroutine  coupler_type_copy_2d_2d

!#######################################################################
!> \brief Copy fields from one coupler type to another. 2-D to 3-D version for generic coupler_type_copy.
!!
!! Template:
!!
!! ~~~~~~~~~~{.f90}
!!   call coupler_type_copy(var_in, var_out, is, ie, js, je, kd, &
!!        diag_name, axes, time, suffix = 'something')
!! ~~~~~~~~~~
!!
!! \throw FATAL, "Number of output fields is non-zero"
!! \throw FATAL, "var_out%bc already associated"
!! \throw FATAL, "var_out%bc([n])%field already associated"
!! \throw FATAL, "var_out%bc([n])%field([m])%values already associated"
!! \throw FATAL, "axes less than 3 elements"
subroutine coupler_type_copy_2d_3d(var_in, var_out, is, ie, js, je, kd, &
     diag_name, axes, time, suffix)

  type(coupler_2d_bc_type), intent(in)    :: var_in !< variable to copy information from
  type(coupler_3d_bc_type), intent(inout) :: var_out !< variable to copy information to
  integer, intent(in)                     :: is !< lower bound of first dimension
  integer, intent(in)                     :: ie !< upper bound of first dimension
  integer, intent(in)                     :: js !< lower bound of second dimension
  integer, intent(in)                     :: je !< upper bound of second dimension
  integer, intent(in)                     :: kd !< third dimension
  character(len=*), intent(in)            :: diag_name !< name for diagnostic file--if blank, then don't register the fields
  integer, dimension(:), intent(in)       :: axes !< array of axes identifiers for diagnostic variable registration
  type(time_type), intent(in)             :: time !< model time variable for registering diagnostic field
  character(len=*), intent(in), optional  :: suffix !< optional suffix to make the name identifier unique

  character(len=256), parameter :: error_header = &
     '==>Error from coupler_types_mod (coupler_type_copy_2d_3d):'
  character(len=400)      :: error_msg
  integer                 :: m, n


  if (var_out%num_bcs > 0) then
    ! It is an error if the number of output fields exceeds zero, because it means this
    ! type has already been populated.
    call mpp_error(FATAL, trim(error_header) // ' Number of output fields exceeds zero')
  endif

  if (var_in%num_bcs >= 0) &
    call CT_spawn_2d_3d(var_in, var_out,  (/ is, is, ie, ie /), (/ js, js, je, je /), (/1, kd/), suffix)

  if ((var_out%num_bcs > 0) .and. (diag_name .ne. ' ')) &
    call CT_set_diags_3d(var_out, diag_name, axes, time)

end subroutine  coupler_type_copy_2d_3d

!#######################################################################
!> \brief Copy fields from one coupler type to another. 3-D to 2-D version for generic coupler_type_copy.
!!
!! Template:
!!
!! ~~~~~~~~~~{.f90}
!!   call coupler_type_copy(var_in, var_out, is, ie, js, je, &
!!        diag_name, axes, time, suffix = 'something')
!! ~~~~~~~~~~
subroutine coupler_type_copy_3d_2d(var_in, var_out, is, ie, js, je,     &
     diag_name, axes, time, suffix)

  type(coupler_3d_bc_type), intent(in)    :: var_in !< variable to copy information from
  type(coupler_2d_bc_type), intent(inout) :: var_out !< variable to copy information to
  integer, intent(in)                     :: is !< lower bound of first dimension
  integer, intent(in)                     :: ie !< upper bound of first dimension
  integer, intent(in)                     :: js !< lower bound of second dimension
  integer, intent(in)                     :: je !< upper bound of second dimension
  character(len=*), intent(in)            :: diag_name !< name for diagnostic file--if blank, then don't register the fields
  integer, dimension(:), intent(in)       :: axes !< array of axes identifiers for diagnostic variable registration
  type(time_type), intent(in)             :: time !< model time variable for registering diagnostic field
  character(len=*), intent(in), optional  :: suffix !< optional suffix to make the name identifier unique

  character(len=256), parameter :: error_header = &
       '==>Error from coupler_types_mod (coupler_type_copy_3d_2d):'
  character(len=400)      :: error_msg
  integer                 :: m, n

  if (var_out%num_bcs > 0) then
    ! It is an error if the number of output fields exceeds zero, because it means this
    ! type has already been populated.
    call mpp_error(FATAL, trim(error_header) // ' Number of output fields exceeds zero')
  endif

  if (var_in%num_bcs >= 0) &
    call CT_spawn_3d_2d(var_in, var_out, (/ is, is, ie, ie /), (/ js, js, je, je /), suffix)

  if ((var_out%num_bcs > 0) .and. (diag_name .ne. ' ')) &
    call CT_set_diags_2d(var_out, diag_name, axes, time)

end subroutine  coupler_type_copy_3d_2d

!#######################################################################
!> \brief Copy fields from one coupler type to another. 3-D to 3-D version for generic coupler_type_copy.
!!
!! Template:
!!
!! ~~~~~~~~~~{.f90}
!!   call coupler_type_copy(var_in, var_out, is, ie, js, je, kd, &
!!        diag_name, axes, time, suffix = 'something')
!! ~~~~~~~~~~
!!
!! \throw FATAL, "Number of output fields is non-zero"
!! \throw FATAL, "var_out%bc already associated"
!! \throw FATAL, "var_out%bc([n])%field already associated"
!! \throw FATAL, "var_out%bc([n])%field([m])%values already associated"
!! \throw FATAL, "axes less than 3 elements"
subroutine coupler_type_copy_3d_3d(var_in, var_out, is, ie, js, je, kd, &
     diag_name, axes, time, suffix)

  type(coupler_3d_bc_type), intent(in)    :: var_in !< variable to copy information from
  type(coupler_3d_bc_type), intent(inout) :: var_out !< variable to copy information to
  integer, intent(in)                     :: is !< lower bound of first dimension
  integer, intent(in)                     :: ie !< upper bound of first dimension
  integer, intent(in)                     :: js !< lower bound of second dimension
  integer, intent(in)                     :: je !< upper bound of second dimension
  integer, intent(in)                     :: kd !< third dimension
  character(len=*), intent(in)            :: diag_name !< name for diagnostic file--if blank, then don't register the fields
  integer, dimension(:), intent(in)       :: axes !< array of axes identifiers for diagnostic variable registration
  type(time_type), intent(in)             :: time !< model time variable for registering diagnostic field
  character(len=*), intent(in), optional  :: suffix !< optional suffix to make the name identifier unique

  character(len=256), parameter :: error_header = &
     '==>Error from coupler_types_mod (coupler_type_copy_3d_3d):'
  character(len=400)      :: error_msg
  integer                 :: m, n


  if (var_out%num_bcs > 0) then
    ! It is an error if the number of output fields exceeds zero, because it means this
    ! type has already been populated.
    call mpp_error(FATAL, trim(error_header) // ' Number of output fields exceeds zero')
  endif

  if (var_in%num_bcs >= 0) &
    call CT_spawn_3d_3d(var_in, var_out,  (/ is, is, ie, ie /), (/ js, js, je, je /), (/1, kd/), suffix)

  if ((var_out%num_bcs > 0) .and. (diag_name .ne. ' ')) &
    call CT_set_diags_3d(var_out, diag_name, axes, time)

end subroutine  coupler_type_copy_3d_3d


!#######################################################################
!> \brief Generate one coupler type using another as a template. 1-D to 2-D version for generic coupler_type_spawn.
!!
!! Template:
!!
!! ~~~~~~~~~~{.f90}
!!   call coupler_type_spawn(var_in, var_out, idim, jdim, suffix = 'something')
!! ~~~~~~~~~~
!!
!! \throw FATAL, "Number of output fields is non-zero"
!! \throw FATAL, "var_out%bc already associated"
!! \throw FATAL, "var_out%bc([n])%field already associated"
!! \throw FATAL, "var_out%bc([n])%field([m])%values already associated"
subroutine CT_spawn_1d_2d(var_in, var, idim, jdim, suffix, as_needed)

  type(coupler_1d_bc_type), intent(in)    :: var_in  !< structure from which to copy information
  type(coupler_2d_bc_type), intent(inout) :: var     !< structure into which to copy information
  integer, dimension(4),    intent(in)    :: idim    !< The data and computational domain extents of
                                                     !! the first dimension in a non-decreasing list
  integer, dimension(4),    intent(in)    :: jdim    !< The data and computational domain extents of
                                                     !! the second dimension in a non-decreasing list
  character(len=*), optional, intent(in)  :: suffix  !< optional suffix to make the name identifier unique
  logical,          optional, intent(in)  :: as_needed !< Only do the spawn if the target type (var)
                                                     !! is not set and the parent type (var_in) is set.

  character(len=256), parameter :: error_header = &
       '==>Error from coupler_types_mod (CT_spawn_1d_2d):'
  character(len=400)      :: error_msg
  integer                 :: m, n

  if (present(as_needed)) then ; if (as_needed) then
    if ((var%set) .or. (.not.var_in%set)) return
  endif ; endif

  if (var%set) &
    call mpp_error(FATAL, trim(error_header) // ' The output type has already been initialized.')
  if (.not.var_in%set) &
    call mpp_error(FATAL, trim(error_header) // ' The parent type has not been initialized.')

  var%num_bcs = var_in%num_bcs ; var%set = .true.

  if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
    write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
    call mpp_error(FATAL, trim(error_msg))
  endif
  if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
    write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list  ', jdim
    call mpp_error(FATAL, trim(error_msg))
  endif
  var%isd = idim(1) ; var%isc = idim(2) ; var%iec = idim(3) ; var%ied = idim(4)
  var%jsd = jdim(1) ; var%jsc = jdim(2) ; var%jec = jdim(3) ; var%jed = jdim(4)

  if (var%num_bcs > 0) then
    if (associated(var%bc)) then
      call mpp_error(FATAL, trim(error_header) // ' var%bc already associated')
    endif
    allocate ( var%bc(var%num_bcs) )
    do n = 1, var%num_bcs
      var%bc(n)%name = var_in%bc(n)%name
      var%bc(n)%atm_tr_index = var_in%bc(n)%atm_tr_index
      var%bc(n)%flux_type = var_in%bc(n)%flux_type
      var%bc(n)%implementation = var_in%bc(n)%implementation
      var%bc(n)%ice_restart_file = var_in%bc(n)%ice_restart_file
      var%bc(n)%ocean_restart_file = var_in%bc(n)%ocean_restart_file
      var%bc(n)%use_atm_pressure = var_in%bc(n)%use_atm_pressure
      var%bc(n)%use_10m_wind_speed = var_in%bc(n)%use_10m_wind_speed
      var%bc(n)%pass_through_ice = var_in%bc(n)%pass_through_ice
      var%bc(n)%mol_wt = var_in%bc(n)%mol_wt
      var%bc(n)%num_fields = var_in%bc(n)%num_fields
      if (associated(var%bc(n)%field)) then
        write (error_msg, *) trim(error_header), ' var%bc(', n, ')%field already associated'
        call mpp_error(FATAL, trim(error_msg))
      endif
      allocate ( var%bc(n)%field(var%bc(n)%num_fields) )
      do m = 1, var%bc(n)%num_fields
        if (present(suffix)) then
          var%bc(n)%field(m)%name = trim(var_in%bc(n)%field(m)%name) // trim(suffix)
        else
          var%bc(n)%field(m)%name = var_in%bc(n)%field(m)%name
        endif
        var%bc(n)%field(m)%long_name = var_in%bc(n)%field(m)%long_name
        var%bc(n)%field(m)%units = var_in%bc(n)%field(m)%units
        var%bc(n)%field(m)%may_init = var_in%bc(n)%field(m)%may_init
        var%bc(n)%field(m)%mean = var_in%bc(n)%field(m)%mean
        if (associated(var%bc(n)%field(m)%values)) then
          write (error_msg, *) trim(error_header), ' var%bc(', n, ')%field(', m, ')%values already associated'
          call mpp_error(FATAL, trim(error_msg))
        endif
        ! Note that this may be allocating a zero-sized array, which is legal in Fortran.
        allocate ( var%bc(n)%field(m)%values(var%isd:var%ied,var%jsd:var%jed) )
        var%bc(n)%field(m)%values(:,:) = 0.0
      enddo
    enddo

  endif

end subroutine  CT_spawn_1d_2d

!#######################################################################
!> \brief Generate one coupler type using another as a template. 1-D to 3-D version for generic CT_spawn.
!!
!! Template:
!!
!! ~~~~~~~~~~{.f90}
!!   call coupler_type_spawn(var_in, var, idim, jdim, kdim, suffix = 'something')
!! ~~~~~~~~~~
!!
!! \throw FATAL, "Number of output fields is non-zero"
!! \throw FATAL, "var%bc already associated"
!! \throw FATAL, "var%bc([n])%field already associated"
!! \throw FATAL, "var%bc([n])%field([m])%values already associated"
subroutine CT_spawn_1d_3d(var_in, var, idim, jdim, kdim, suffix, as_needed)

  type(coupler_1d_bc_type), intent(in)    :: var_in  !< structure from which to copy information
  type(coupler_3d_bc_type), intent(inout) :: var     !< structure into which to copy information
  integer, dimension(4),    intent(in)    :: idim    !< The data and computational domain extents of
                                                     !! the first dimension in a non-decreasing list
  integer, dimension(4),    intent(in)    :: jdim    !< The data and computational domain extents of
                                                     !! the second dimension in a non-decreasing list
  integer, dimension(2),    intent(in)    :: kdim    !< The array extents of the third dimension in
                                                     !! a non-decreasing list
  character(len=*), optional, intent(in)  :: suffix  !< optional suffix to make the name identifier unique
  logical,          optional, intent(in)  :: as_needed !< Only do the spawn if the target type (var)
                                                     !! is not set and the parent type (var_in) is set.

  character(len=256), parameter :: error_header = &
     '==>Error from coupler_types_mod (CT_spawn_1d_3d):'
  character(len=400)      :: error_msg
  integer                 :: m, n

  if (present(as_needed)) then ; if (as_needed) then
    if ((var%set) .or. (.not.var_in%set)) return
  endif ; endif

  if (var%set) &
    call mpp_error(FATAL, trim(error_header) // ' The output type has already been initialized.')
  if (.not.var_in%set) &
    call mpp_error(FATAL, trim(error_header) // ' The parent type has not been initialized.')

  var%num_bcs = var_in%num_bcs ; var%set = .true.

  ! Store the array extents that are to be used with this bc_type.
  if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
    write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
    call mpp_error(FATAL, trim(error_msg))
  endif
  if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
    write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list  ', jdim
    call mpp_error(FATAL, trim(error_msg))
  endif
  if (kdim(1) > kdim(2)) then
    write (error_msg, *) trim(error_header), ' Disordered k-dimension index bound list  ', kdim
    call mpp_error(FATAL, trim(error_msg))
  endif
  var%isd = idim(1) ; var%isc = idim(2) ; var%iec = idim(3) ; var%ied = idim(4)
  var%jsd = jdim(1) ; var%jsc = jdim(2) ; var%jec = jdim(3) ; var%jed = jdim(4)
  var%ks  = kdim(1) ; var%ke  = kdim(2)

  if (var%num_bcs > 0) then
    if (associated(var%bc)) then
      call mpp_error(FATAL, trim(error_header) // ' var%bc already associated')
    endif
    allocate ( var%bc(var%num_bcs) )
    do n = 1, var%num_bcs
      var%bc(n)%name = var_in%bc(n)%name
      var%bc(n)%atm_tr_index = var_in%bc(n)%atm_tr_index
      var%bc(n)%flux_type = var_in%bc(n)%flux_type
      var%bc(n)%implementation = var_in%bc(n)%implementation
      var%bc(n)%ice_restart_file = var_in%bc(n)%ice_restart_file
      var%bc(n)%ocean_restart_file = var_in%bc(n)%ocean_restart_file
      var%bc(n)%use_atm_pressure = var_in%bc(n)%use_atm_pressure
      var%bc(n)%use_10m_wind_speed = var_in%bc(n)%use_10m_wind_speed
      var%bc(n)%pass_through_ice = var_in%bc(n)%pass_through_ice
      var%bc(n)%mol_wt = var_in%bc(n)%mol_wt
      var%bc(n)%num_fields = var_in%bc(n)%num_fields
      if (associated(var%bc(n)%field)) then
        write (error_msg, *) trim(error_header), ' var%bc(', n, ')%field already associated'
        call mpp_error(FATAL, trim(error_msg))
      endif
      allocate ( var%bc(n)%field(var%bc(n)%num_fields) )
      do m = 1, var%bc(n)%num_fields
        if (present(suffix)) then
          var%bc(n)%field(m)%name = trim(var_in%bc(n)%field(m)%name) // trim(suffix)
        else
          var%bc(n)%field(m)%name = var_in%bc(n)%field(m)%name
        endif
        var%bc(n)%field(m)%long_name = var_in%bc(n)%field(m)%long_name
        var%bc(n)%field(m)%units = var_in%bc(n)%field(m)%units
        var%bc(n)%field(m)%may_init = var_in%bc(n)%field(m)%may_init
        var%bc(n)%field(m)%mean = var_in%bc(n)%field(m)%mean
        if (associated(var%bc(n)%field(m)%values)) then
          write (error_msg, *) trim(error_header), ' var%bc(', n, ')%field(', m, ')%values already associated'
          call mpp_error(FATAL, trim(error_msg))
        endif
        ! Note that this may be allocating a zero-sized array, which is legal in Fortran.
        allocate ( var%bc(n)%field(m)%values(var%isd:var%ied,var%jsd:var%jed,var%ks:var%ke) )
        var%bc(n)%field(m)%values(:,:,:) = 0.0
      enddo
    enddo

  endif

end subroutine  CT_spawn_1d_3d

!#######################################################################
!> \brief Generate one coupler type using another as a template. 2-D to 2-D version for generic CT_spawn.
!!
!! Template:
!!
!! ~~~~~~~~~~{.f90}
!!   call coupler_type_spawn(var_in, var, idim, jdim, suffix = 'something')
!! ~~~~~~~~~~
!!
!! \throw FATAL, "Number of output fields is non-zero"
!! \throw FATAL, "var%bc already associated"
!! \throw FATAL, "var%bc([n])%field already associated"
!! \throw FATAL, "var%bc([n])%field([m])%values already associated"
subroutine CT_spawn_2d_2d(var_in, var, idim, jdim, suffix, as_needed)

  type(coupler_2d_bc_type), intent(in)    :: var_in  !< structure from which to copy information
  type(coupler_2d_bc_type), intent(inout) :: var     !< structure into which to copy information
  integer, dimension(4),    intent(in)    :: idim    !< The data and computational domain extents of
                                                     !! the first dimension in a non-decreasing list
  integer, dimension(4),    intent(in)    :: jdim    !< The data and computational domain extents of
                                                     !! the second dimension in a non-decreasing list
  character(len=*), optional, intent(in)  :: suffix  !< optional suffix to make the name identifier unique
  logical,          optional, intent(in)  :: as_needed !< Only do the spawn if the target type (var)
                                                     !! is not set and the parent type (var_in) is set.

  character(len=256), parameter :: error_header = &
       '==>Error from coupler_types_mod (CT_spawn_2d_2d):'
  character(len=400)      :: error_msg
  integer                 :: m, n

  if (present(as_needed)) then ; if (as_needed) then
    if ((var%set) .or. (.not.var_in%set)) return
  endif ; endif

  if (var%set) &
    call mpp_error(FATAL, trim(error_header) // ' The output type has already been initialized.')
  if (.not.var_in%set) &
    call mpp_error(FATAL, trim(error_header) // ' The parent type has not been initialized.')

  var%num_bcs = var_in%num_bcs ; var%set = .true.

  if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
    write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
    call mpp_error(FATAL, trim(error_msg))
  endif
  if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
    write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list  ', jdim
    call mpp_error(FATAL, trim(error_msg))
  endif
  var%isd = idim(1) ; var%isc = idim(2) ; var%iec = idim(3) ; var%ied = idim(4)
  var%jsd = jdim(1) ; var%jsc = jdim(2) ; var%jec = jdim(3) ; var%jed = jdim(4)

  if (var%num_bcs > 0) then
    if (associated(var%bc)) then
      call mpp_error(FATAL, trim(error_header) // ' var%bc already associated')
    endif
    allocate ( var%bc(var%num_bcs) )
    do n = 1, var%num_bcs
      var%bc(n)%name = var_in%bc(n)%name
      var%bc(n)%atm_tr_index = var_in%bc(n)%atm_tr_index
      var%bc(n)%flux_type = var_in%bc(n)%flux_type
      var%bc(n)%implementation = var_in%bc(n)%implementation
      var%bc(n)%ice_restart_file = var_in%bc(n)%ice_restart_file
      var%bc(n)%ocean_restart_file = var_in%bc(n)%ocean_restart_file
      var%bc(n)%use_atm_pressure = var_in%bc(n)%use_atm_pressure
      var%bc(n)%use_10m_wind_speed = var_in%bc(n)%use_10m_wind_speed
      var%bc(n)%pass_through_ice = var_in%bc(n)%pass_through_ice
      var%bc(n)%mol_wt = var_in%bc(n)%mol_wt
      var%bc(n)%num_fields = var_in%bc(n)%num_fields
      if (associated(var%bc(n)%field)) then
        write (error_msg, *) trim(error_header), ' var%bc(', n, ')%field already associated'
        call mpp_error(FATAL, trim(error_msg))
      endif
      allocate ( var%bc(n)%field(var%bc(n)%num_fields) )
      do m = 1, var%bc(n)%num_fields
        if (present(suffix)) then
          var%bc(n)%field(m)%name = trim(var_in%bc(n)%field(m)%name) // trim(suffix)
        else
          var%bc(n)%field(m)%name = var_in%bc(n)%field(m)%name
        endif
        var%bc(n)%field(m)%long_name = var_in%bc(n)%field(m)%long_name
        var%bc(n)%field(m)%units = var_in%bc(n)%field(m)%units
        var%bc(n)%field(m)%may_init = var_in%bc(n)%field(m)%may_init
        var%bc(n)%field(m)%mean = var_in%bc(n)%field(m)%mean
        if (associated(var%bc(n)%field(m)%values)) then
          write (error_msg, *) trim(error_header), ' var%bc(', n, ')%field(', m, ')%values already associated'
          call mpp_error(FATAL, trim(error_msg))
        endif
        ! Note that this may be allocating a zero-sized array, which is legal in Fortran.
        allocate ( var%bc(n)%field(m)%values(var%isd:var%ied,var%jsd:var%jed) )
        var%bc(n)%field(m)%values(:,:) = 0.0
      enddo
    enddo

  endif

end subroutine  CT_spawn_2d_2d

!#######################################################################
!> \brief Generate one coupler type using another as a template. 2-D to 3-D version for generic CT_spawn.
!!
!! Template:
!!
!! ~~~~~~~~~~{.f90}
!!   call coupler_type_spawn(var_in, var, idim, jdim, kdim, suffix = 'something')
!! ~~~~~~~~~~
!!
!! \throw FATAL, "Number of output fields is non-zero"
!! \throw FATAL, "var%bc already associated"
!! \throw FATAL, "var%bc([n])%field already associated"
!! \throw FATAL, "var%bc([n])%field([m])%values already associated"
subroutine CT_spawn_2d_3d(var_in, var, idim, jdim, kdim, suffix, as_needed)

  type(coupler_2d_bc_type), intent(in)    :: var_in  !< structure from which to copy information
  type(coupler_3d_bc_type), intent(inout) :: var     !< structure into which to copy information
  integer, dimension(4),    intent(in)    :: idim    !< The data and computational domain extents of
                                                     !! the first dimension in a non-decreasing list
  integer, dimension(4),    intent(in)    :: jdim    !< The data and computational domain extents of
                                                     !! the second dimension in a non-decreasing list
  integer, dimension(2),    intent(in)    :: kdim    !< The array extents of the third dimension in
                                                     !! a non-decreasing list
  character(len=*), optional, intent(in)  :: suffix  !< optional suffix to make the name identifier unique
  logical,          optional, intent(in)  :: as_needed !< Only do the spawn if the target type (var)
                                                     !! is not set and the parent type (var_in) is set.

  character(len=256), parameter :: error_header = &
     '==>Error from coupler_types_mod (CT_spawn_2d_3d):'
  character(len=400)      :: error_msg
  integer                 :: m, n

  if (present(as_needed)) then ; if (as_needed) then
    if ((var%set) .or. (.not.var_in%set)) return
  endif ; endif

  if (var%set) &
    call mpp_error(FATAL, trim(error_header) // ' The output type has already been initialized.')
  if (.not.var_in%set) &
    call mpp_error(FATAL, trim(error_header) // ' The parent type has not been initialized.')

  var%num_bcs = var_in%num_bcs ; var%set = .true.

  ! Store the array extents that are to be used with this bc_type.
  if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
    write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
    call mpp_error(FATAL, trim(error_msg))
  endif
  if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
    write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list  ', jdim
    call mpp_error(FATAL, trim(error_msg))
  endif
  if (kdim(1) > kdim(2)) then
    write (error_msg, *) trim(error_header), ' Disordered k-dimension index bound list  ', kdim
    call mpp_error(FATAL, trim(error_msg))
  endif
  var%isd = idim(1) ; var%isc = idim(2) ; var%iec = idim(3) ; var%ied = idim(4)
  var%jsd = jdim(1) ; var%jsc = jdim(2) ; var%jec = jdim(3) ; var%jed = jdim(4)
  var%ks  = kdim(1) ; var%ke = kdim(2)

  if (var%num_bcs > 0) then
    if (associated(var%bc)) then
      call mpp_error(FATAL, trim(error_header) // ' var%bc already associated')
    endif
    allocate ( var%bc(var%num_bcs) )
    do n = 1, var%num_bcs
      var%bc(n)%name = var_in%bc(n)%name
      var%bc(n)%atm_tr_index = var_in%bc(n)%atm_tr_index
      var%bc(n)%flux_type = var_in%bc(n)%flux_type
      var%bc(n)%implementation = var_in%bc(n)%implementation
      var%bc(n)%ice_restart_file = var_in%bc(n)%ice_restart_file
      var%bc(n)%ocean_restart_file = var_in%bc(n)%ocean_restart_file
      var%bc(n)%use_atm_pressure = var_in%bc(n)%use_atm_pressure
      var%bc(n)%use_10m_wind_speed = var_in%bc(n)%use_10m_wind_speed
      var%bc(n)%pass_through_ice = var_in%bc(n)%pass_through_ice
      var%bc(n)%mol_wt = var_in%bc(n)%mol_wt
      var%bc(n)%num_fields = var_in%bc(n)%num_fields
      if (associated(var%bc(n)%field)) then
        write (error_msg, *) trim(error_header), ' var%bc(', n, ')%field already associated'
        call mpp_error(FATAL, trim(error_msg))
      endif
      allocate ( var%bc(n)%field(var%bc(n)%num_fields) )
      do m = 1, var%bc(n)%num_fields
        if (present(suffix)) then
          var%bc(n)%field(m)%name = trim(var_in%bc(n)%field(m)%name) // trim(suffix)
        else
          var%bc(n)%field(m)%name = var_in%bc(n)%field(m)%name
        endif
        var%bc(n)%field(m)%long_name = var_in%bc(n)%field(m)%long_name
        var%bc(n)%field(m)%units = var_in%bc(n)%field(m)%units
        var%bc(n)%field(m)%may_init = var_in%bc(n)%field(m)%may_init
        var%bc(n)%field(m)%mean = var_in%bc(n)%field(m)%mean
        if (associated(var%bc(n)%field(m)%values)) then
          write (error_msg, *) trim(error_header), ' var%bc(', n, ')%field(', m, ')%values already associated'
          call mpp_error(FATAL, trim(error_msg))
        endif
        ! Note that this may be allocating a zero-sized array, which is legal in Fortran.
        allocate ( var%bc(n)%field(m)%values(var%isd:var%ied,var%jsd:var%jed,var%ks:var%ke) )
        var%bc(n)%field(m)%values(:,:,:) = 0.0
      enddo
    enddo

  endif

end subroutine  CT_spawn_2d_3d

!#######################################################################
!> \brief Generate one coupler type using another as a template. 3-D to 2-D version for generic CT_spawn.
!!
!! Template:
!!
!! ~~~~~~~~~~{.f90}
!!   call coupler_type_spawn(var_in, var, idim, jdim, suffix = 'something')
!! ~~~~~~~~~~
!!
!! \throw FATAL, "Number of output fields is non-zero"
!! \throw FATAL, "var%bc already associated"
!! \throw FATAL, "var%bc([n])%field already associated"
!! \throw FATAL, "var%bc([n])%field([m])%values already associated"
subroutine CT_spawn_3d_2d(var_in, var, idim, jdim, suffix, as_needed)

  type(coupler_3d_bc_type), intent(in)    :: var_in  !< structure from which to copy information
  type(coupler_2d_bc_type), intent(inout) :: var     !< structure into which to copy information
  integer, dimension(4),    intent(in)    :: idim    !< The data and computational domain extents of
                                                     !! the first dimension in a non-decreasing list
  integer, dimension(4),    intent(in)    :: jdim    !< The data and computational domain extents of
                                                     !! the second dimension in a non-decreasing list
  character(len=*), optional, intent(in)  :: suffix  !< optional suffix to make the name identifier unique
  logical,          optional, intent(in)  :: as_needed !< Only do the spawn if the target type (var)
                                                     !! is not set and the parent type (var_in) is set.

  character(len=256), parameter :: error_header = &
       '==>Error from coupler_types_mod (CT_spawn_3d_2d):'
  character(len=400)      :: error_msg
  integer                 :: m, n

  if (present(as_needed)) then ; if (as_needed) then
    if ((var%set) .or. (.not.var_in%set)) return
  endif ; endif

  if (var%set) &
    call mpp_error(FATAL, trim(error_header) // ' The output type has already been initialized.')
  if (.not.var_in%set) &
    call mpp_error(FATAL, trim(error_header) // ' The parent type has not been initialized.')

  var%num_bcs = var_in%num_bcs ; var%set = .true.

  if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
    write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
    call mpp_error(FATAL, trim(error_msg))
  endif
  if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
    write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list  ', jdim
    call mpp_error(FATAL, trim(error_msg))
  endif
  var%isd = idim(1) ; var%isc = idim(2) ; var%iec = idim(3) ; var%ied = idim(4)
  var%jsd = jdim(1) ; var%jsc = jdim(2) ; var%jec = jdim(3) ; var%jed = jdim(4)

  if (var%num_bcs > 0) then
    if (associated(var%bc)) then
      call mpp_error(FATAL, trim(error_header) // ' var%bc already associated')
    endif
    allocate ( var%bc(var%num_bcs) )
    do n = 1, var%num_bcs
      var%bc(n)%name = var_in%bc(n)%name
      var%bc(n)%atm_tr_index = var_in%bc(n)%atm_tr_index
      var%bc(n)%flux_type = var_in%bc(n)%flux_type
      var%bc(n)%implementation = var_in%bc(n)%implementation
      var%bc(n)%ice_restart_file = var_in%bc(n)%ice_restart_file
      var%bc(n)%ocean_restart_file = var_in%bc(n)%ocean_restart_file
      var%bc(n)%use_atm_pressure = var_in%bc(n)%use_atm_pressure
      var%bc(n)%use_10m_wind_speed = var_in%bc(n)%use_10m_wind_speed
      var%bc(n)%pass_through_ice = var_in%bc(n)%pass_through_ice
      var%bc(n)%mol_wt = var_in%bc(n)%mol_wt
      var%bc(n)%num_fields = var_in%bc(n)%num_fields
      if (associated(var%bc(n)%field)) then
        write (error_msg, *) trim(error_header), ' var%bc(', n, ')%field already associated'
        call mpp_error(FATAL, trim(error_msg))
      endif
      allocate ( var%bc(n)%field(var%bc(n)%num_fields) )
      do m = 1, var%bc(n)%num_fields
        if (present(suffix)) then
          var%bc(n)%field(m)%name = trim(var_in%bc(n)%field(m)%name) // trim(suffix)
        else
          var%bc(n)%field(m)%name = var_in%bc(n)%field(m)%name
        endif
        var%bc(n)%field(m)%long_name = var_in%bc(n)%field(m)%long_name
        var%bc(n)%field(m)%units = var_in%bc(n)%field(m)%units
        var%bc(n)%field(m)%may_init = var_in%bc(n)%field(m)%may_init
        var%bc(n)%field(m)%mean = var_in%bc(n)%field(m)%mean
        if (associated(var%bc(n)%field(m)%values)) then
          write (error_msg, *) trim(error_header), ' var%bc(', n, ')%field(', m, ')%values already associated'
          call mpp_error(FATAL, trim(error_msg))
        endif
        ! Note that this may be allocating a zero-sized array, which is legal in Fortran.
        allocate ( var%bc(n)%field(m)%values(var%isd:var%ied,var%jsd:var%jed) )
        var%bc(n)%field(m)%values(:,:) = 0.0
      enddo
    enddo

  endif

end subroutine  CT_spawn_3d_2d

!#######################################################################
!> \brief Generate one coupler type using another as a template. 3-D to 3-D version for generic CT_spawn.
!!
!! Template:
!!
!! ~~~~~~~~~~{.f90}
!!   call coupler_type_spawn(var_in, var, idim, jdim, kdim, suffix = 'something')
!! ~~~~~~~~~~
!!
!! \throw FATAL, "Number of output fields is non-zero"
!! \throw FATAL, "var%bc already associated"
!! \throw FATAL, "var%bc([n])%field already associated"
!! \throw FATAL, "var%bc([n])%field([m])%values already associated"
subroutine CT_spawn_3d_3d(var_in, var, idim, jdim, kdim, suffix, as_needed)

  type(coupler_3d_bc_type), intent(in)    :: var_in  !< structure from which to copy information
  type(coupler_3d_bc_type), intent(inout) :: var     !< structure into which to copy information
  integer, dimension(4),    intent(in)    :: idim    !< The data and computational domain extents of
                                                     !! the first dimension in a non-decreasing list
  integer, dimension(4),    intent(in)    :: jdim    !< The data and computational domain extents of
                                                     !! the second dimension in a non-decreasing list
  integer, dimension(2),    intent(in)    :: kdim    !< The array extents of the third dimension in
                                                     !! a non-decreasing list
  character(len=*), optional, intent(in)  :: suffix  !< optional suffix to make the name identifier unique
  logical,          optional, intent(in)  :: as_needed !< Only do the spawn if the target type (var)
                                                     !! is not set and the parent type (var_in) is set.

  character(len=256), parameter :: error_header = &
     '==>Error from coupler_types_mod (CT_spawn_3d_3d):'
  character(len=400)      :: error_msg
  integer                 :: m, n

  if (present(as_needed)) then ; if (as_needed) then
    if ((var%set) .or. (.not.var_in%set)) return
  endif ; endif

  if (var%set) &
    call mpp_error(FATAL, trim(error_header) // ' The output type has already been initialized.')
  if (.not.var_in%set) &
    call mpp_error(FATAL, trim(error_header) // ' The parent type has not been initialized.')

  var%num_bcs = var_in%num_bcs ; var%set = .true.

  if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
    write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
    call mpp_error(FATAL, trim(error_msg))
  endif
  if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
    write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list  ', jdim
    call mpp_error(FATAL, trim(error_msg))
  endif
  if (kdim(1) > kdim(2)) then
    write (error_msg, *) trim(error_header), ' Disordered k-dimension index bound list  ', kdim
    call mpp_error(FATAL, trim(error_msg))
  endif
  var%isd = idim(1) ; var%isc = idim(2) ; var%iec = idim(3) ; var%ied = idim(4)
  var%jsd = jdim(1) ; var%jsc = jdim(2) ; var%jec = jdim(3) ; var%jed = jdim(4)
  var%ks  = kdim(1) ; var%ke  = kdim(2)

  if (var%num_bcs > 0) then
    if (associated(var%bc)) then
      call mpp_error(FATAL, trim(error_header) // ' var%bc already associated')
    endif
    allocate ( var%bc(var%num_bcs) )
    do n = 1, var%num_bcs
      var%bc(n)%name = var_in%bc(n)%name
      var%bc(n)%atm_tr_index = var_in%bc(n)%atm_tr_index
      var%bc(n)%flux_type = var_in%bc(n)%flux_type
      var%bc(n)%implementation = var_in%bc(n)%implementation
      var%bc(n)%ice_restart_file = var_in%bc(n)%ice_restart_file
      var%bc(n)%ocean_restart_file = var_in%bc(n)%ocean_restart_file
      var%bc(n)%use_atm_pressure = var_in%bc(n)%use_atm_pressure
      var%bc(n)%use_10m_wind_speed = var_in%bc(n)%use_10m_wind_speed
      var%bc(n)%pass_through_ice = var_in%bc(n)%pass_through_ice
      var%bc(n)%mol_wt = var_in%bc(n)%mol_wt
      var%bc(n)%num_fields = var_in%bc(n)%num_fields
      if (associated(var%bc(n)%field)) then
        write (error_msg, *) trim(error_header), ' var%bc(', n, ')%field already associated'
        call mpp_error(FATAL, trim(error_msg))
      endif
      allocate ( var%bc(n)%field(var%bc(n)%num_fields) )
      do m = 1, var%bc(n)%num_fields
        if (present(suffix)) then
          var%bc(n)%field(m)%name = trim(var_in%bc(n)%field(m)%name) // trim(suffix)
        else
          var%bc(n)%field(m)%name = var_in%bc(n)%field(m)%name
        endif
        var%bc(n)%field(m)%long_name = var_in%bc(n)%field(m)%long_name
        var%bc(n)%field(m)%units = var_in%bc(n)%field(m)%units
        var%bc(n)%field(m)%may_init = var_in%bc(n)%field(m)%may_init
        var%bc(n)%field(m)%mean = var_in%bc(n)%field(m)%mean
        if (associated(var%bc(n)%field(m)%values)) then
          write (error_msg, *) trim(error_header), ' var%bc(', n, ')%field(', m, ')%values already associated'
          call mpp_error(FATAL, trim(error_msg))
        endif

        ! Note that this may be allocating a zero-sized array, which is legal in Fortran.
        allocate ( var%bc(n)%field(m)%values(var%isd:var%ied,var%jsd:var%jed,var%ks:var%ke) )
        var%bc(n)%field(m)%values(:,:,:) = 0.0
      enddo
    enddo

  endif

end subroutine  CT_spawn_3d_3d


!> This subroutine does a direct copy of the data in all elements of one
!! coupler_2d_bc_type into another.  Both must have the same array sizes.
subroutine CT_copy_data_2d(var_in, var, halo_size, bc_index, field_index, &
                           exclude_flux_type, only_flux_type, pass_through_ice)
  type(coupler_2d_bc_type),   intent(in)    :: var_in  !< BC_type structure with the data to copy
  type(coupler_2d_bc_type),   intent(inout) :: var !< The recipient BC_type structure
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                       !! that is being copied
  integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                       !! boundary condition that is being copied
  character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types of fluxes to exclude from this copy.
  character(len=*), optional, intent(in)    :: only_flux_type    !< A string describing which types of fluxes to include from this copy.
  logical,          optional, intent(in)    :: pass_through_ice !< If true, only copy BCs whose
                                                       !! value of pass_through ice matches this
  logical :: copy_bc
  integer :: i, j, m, n, n1, n2, halo, i_off, j_off

  if (present(bc_index)) then
    if (bc_index > var_in%num_bcs) &
      call mpp_error(FATAL, "CT_copy_data_2d: bc_index is present and exceeds var_in%num_bcs.")
    if (present(field_index)) then ; if (field_index > var_in%bc(bc_index)%num_fields) &
      call mpp_error(FATAL, "CT_copy_data_2d: field_index is present and exceeds num_fields for" //&
                     trim(var_in%bc(bc_index)%name) )
    endif
  elseif (present(field_index)) then
    call mpp_error(FATAL, "CT_copy_data_2d: bc_index must be present if field_index is present.")
  endif

  halo = 0 ; if (present(halo_size)) halo = halo_size

  n1 = 1 ; n2 = var_in%num_bcs
  if (present(bc_index)) then ; n1 = bc_index ; n2 = bc_index ; endif

  if (n2 >= n1) then
    ! A more consciencious implementation would include a more descriptive error messages.
    if ((var_in%iec-var_in%isc) /= (var%iec-var%isc)) &
      call mpp_error(FATAL, "CT_copy_data_2d: There is an i-direction computional domain size mismatch.")
    if ((var_in%jec-var_in%jsc) /= (var%jec-var%jsc)) &
      call mpp_error(FATAL, "CT_copy_data_2d: There is a j-direction computional domain size mismatch.")
    if ((var_in%isc-var_in%isd < halo) .or. (var_in%ied-var_in%iec < halo)) &
      call mpp_error(FATAL, "CT_copy_data_2d: Excessive i-direction halo size for the input structure.")
    if ((var_in%jsc-var_in%jsd < halo) .or. (var_in%jed-var_in%jec < halo)) &
      call mpp_error(FATAL, "CT_copy_data_2d: Excessive j-direction halo size for the input structure.")
    if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo)) &
      call mpp_error(FATAL, "CT_copy_data_2d: Excessive i-direction halo size for the output structure.")
    if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo)) &
      call mpp_error(FATAL, "CT_copy_data_2d: Excessive j-direction halo size for the output structure.")

    i_off = var_in%isc - var%isc ; j_off = var_in%jsc - var%jsc
  endif

  do n = n1, n2

    copy_bc = .true.
    if (copy_bc .and. present(exclude_flux_type)) &
      copy_bc = .not.(trim(var%bc(n)%flux_type) == trim(exclude_flux_type))
    if (copy_bc .and. present(only_flux_type)) &
      copy_bc = (trim(var%bc(n)%flux_type) == trim(only_flux_type))
    if (copy_bc .and. present(pass_through_ice)) &
      copy_bc = (pass_through_ice .eqv. var%bc(n)%pass_through_ice)
    if (.not.copy_bc) cycle

    do m = 1, var%bc(n)%num_fields
      if (present(field_index)) then ; if (m /= field_index) cycle ; endif
      if ( associated(var%bc(n)%field(m)%values) ) then
        do j=var%jsc-halo,var%jec+halo ; do i=var%isc-halo,var%iec+halo
          var%bc(n)%field(m)%values(i,j) = var_in%bc(n)%field(m)%values(i+i_off,j+j_off)
        enddo ; enddo
      endif
    enddo
  enddo

end subroutine CT_copy_data_2d

!> This subroutine does a direct copy of the data in all elements of one
!! coupler_3d_bc_type into another.  Both types must have the same array sizes.
subroutine CT_copy_data_3d(var_in, var, halo_size, bc_index, field_index, &
                           exclude_flux_type, only_flux_type, pass_through_ice)
  type(coupler_3d_bc_type),   intent(in)    :: var_in  !< BC_type structure with the data to copy
  type(coupler_3d_bc_type),   intent(inout) :: var !< The recipient BC_type structure
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                       !! that is being copied
  integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                       !! boundary condition that is being copied
  character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types of fluxes to exclude from this copy.
  character(len=*), optional, intent(in)    :: only_flux_type    !< A string describing which types of fluxes to include from this copy.
  logical,          optional, intent(in)    :: pass_through_ice !< If true, only copy BCs whose
                                                       !! value of pass_through ice matches this
  logical :: copy_bc
  integer :: i, j, k, m, n, n1, n2, halo, i_off, j_off, k_off

  if (present(bc_index)) then
    if (bc_index > var_in%num_bcs) &
      call mpp_error(FATAL, "CT_copy_data_3d: bc_index is present and exceeds var_in%num_bcs.")
    if (present(field_index)) then ; if (field_index > var_in%bc(bc_index)%num_fields) &
      call mpp_error(FATAL, "CT_copy_data_3d: field_index is present and exceeds num_fields for" //&
                     trim(var_in%bc(bc_index)%name) )
    endif
  elseif (present(field_index)) then
    call mpp_error(FATAL, "CT_copy_data_3d: bc_index must be present if field_index is present.")
  endif

  halo = 0 ; if (present(halo_size)) halo = halo_size

  n1 = 1 ; n2 = var_in%num_bcs
  if (present(bc_index)) then ; n1 = bc_index ; n2 = bc_index ; endif

  if (n2 >= n1) then
    ! A more consciencious implementation would include a more descriptive error messages.
    if ((var_in%iec-var_in%isc) /= (var%iec-var%isc)) &
      call mpp_error(FATAL, "CT_copy_data_3d: There is an i-direction computional domain size mismatch.")
    if ((var_in%jec-var_in%jsc) /= (var%jec-var%jsc)) &
      call mpp_error(FATAL, "CT_copy_data_3d: There is a j-direction computional domain size mismatch.")
    if ((var_in%ke-var_in%ks) /= (var%ke-var%ks)) &
      call mpp_error(FATAL, "CT_copy_data_3d: There is a k-direction computional domain size mismatch.")
    if ((var_in%isc-var_in%isd < halo) .or. (var_in%ied-var_in%iec < halo)) &
      call mpp_error(FATAL, "CT_copy_data_3d: Excessive i-direction halo size for the input structure.")
    if ((var_in%jsc-var_in%jsd < halo) .or. (var_in%jed-var_in%jec < halo)) &
      call mpp_error(FATAL, "CT_copy_data_3d: Excessive j-direction halo size for the input structure.")
    if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo)) &
      call mpp_error(FATAL, "CT_copy_data_3d: Excessive i-direction halo size for the output structure.")
    if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo)) &
      call mpp_error(FATAL, "CT_copy_data_3d: Excessive j-direction halo size for the output structure.")

    i_off = var_in%isc - var%isc ; j_off = var_in%jsc - var%jsc ; k_off = var_in%ks - var%ks
  endif

  do n = n1, n2

    copy_bc = .true.
    if (copy_bc .and. present(exclude_flux_type)) &
      copy_bc = .not.(trim(var%bc(n)%flux_type) == trim(exclude_flux_type))
    if (copy_bc .and. present(only_flux_type)) &
      copy_bc = (trim(var%bc(n)%flux_type) == trim(only_flux_type))
    if (copy_bc .and. present(pass_through_ice)) &
      copy_bc = (pass_through_ice .eqv. var%bc(n)%pass_through_ice)
    if (.not.copy_bc) cycle

    do m = 1, var_in%bc(n)%num_fields
      if (present(field_index)) then ; if (m /= field_index) cycle ; endif
      if ( associated(var%bc(n)%field(m)%values) ) then
        do k=var%ks,var%ke ; do j=var%jsc-halo,var%jec+halo ; do i=var%isc-halo,var%iec+halo
          var%bc(n)%field(m)%values(i,j,k) = var_in%bc(n)%field(m)%values(i+i_off,j+j_off,k+k_off)
        enddo ; enddo ; enddo
      endif
    enddo
  enddo

end subroutine CT_copy_data_3d

!> This subroutine does a direct copy of the data in all elements of a
!! coupler_2d_bc_type into a coupler_3d_bc_type.  Both types must have the same
!! array sizes for their first two dimensions, while the extent of the 3rd dimension
!! that is being filled may be specified via optional arguments.
subroutine CT_copy_data_2d_3d(var_in, var, halo_size, bc_index, field_index, &
                           exclude_flux_type, only_flux_type, pass_through_ice, &
                           ind3_start, ind3_end)
  type(coupler_2d_bc_type),   intent(in)    :: var_in  !< BC_type structure with the data to copy
  type(coupler_3d_bc_type),   intent(inout) :: var !< The recipient BC_type structure
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                       !! that is being copied
  integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                       !! boundary condition that is being copied
  character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types of fluxes to exclude from this copy.
  character(len=*), optional, intent(in)    :: only_flux_type    !< A string describing which types of fluxes to include from this copy.
  logical,          optional, intent(in)    :: pass_through_ice !< If true, only copy BCs whose
                                                       !! value of pass_through ice matches this
  integer,          optional, intent(in)    :: ind3_start  !< The starting value of the 3rd
                                                       !! index of the 3d type to fill in.
  integer,          optional, intent(in)    :: ind3_end    !< The ending value of the 3rd
                                                       !! index of the 3d type to fill in.
  logical :: copy_bc
  integer :: i, j, k, m, n, n1, n2, halo, i_off, j_off, ks, ke

  if (present(bc_index)) then
    if (bc_index > var_in%num_bcs) &
      call mpp_error(FATAL, "CT_copy_data_2d_3d: bc_index is present and exceeds var_in%num_bcs.")
    if (present(field_index)) then ; if (field_index > var_in%bc(bc_index)%num_fields) &
      call mpp_error(FATAL, "CT_copy_data_2d_3d: field_index is present and exceeds num_fields for" //&
                     trim(var_in%bc(bc_index)%name) )
    endif
  elseif (present(field_index)) then
    call mpp_error(FATAL, "CT_copy_data_2d_3d: bc_index must be present if field_index is present.")
  endif

  halo = 0 ; if (present(halo_size)) halo = halo_size

  n1 = 1 ; n2 = var_in%num_bcs
  if (present(bc_index)) then ; n1 = bc_index ; n2 = bc_index ; endif

  if (n2 >= n1) then
    ! A more consciencious implementation would include a more descriptive error messages.
    if ((var_in%iec-var_in%isc) /= (var%iec-var%isc)) &
      call mpp_error(FATAL, "CT_copy_data_2d_3d: There is an i-direction computional domain size mismatch.")
    if ((var_in%jec-var_in%jsc) /= (var%jec-var%jsc)) &
      call mpp_error(FATAL, "CT_copy_data_2d_3d: There is a j-direction computional domain size mismatch.")
    if ((var_in%isc-var_in%isd < halo) .or. (var_in%ied-var_in%iec < halo)) &
      call mpp_error(FATAL, "CT_copy_data_2d_3d: Excessive i-direction halo size for the input structure.")
    if ((var_in%jsc-var_in%jsd < halo) .or. (var_in%jed-var_in%jec < halo)) &
      call mpp_error(FATAL, "CT_copy_data_2d_3d: Excessive j-direction halo size for the input structure.")
    if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo)) &
      call mpp_error(FATAL, "CT_copy_data_2d_3d: Excessive i-direction halo size for the output structure.")
    if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo)) &
      call mpp_error(FATAL, "CT_copy_data_2d_3d: Excessive j-direction halo size for the output structure.")
  endif

  i_off = var_in%isc - var%isc ; j_off = var_in%jsc - var%jsc
  do n = n1, n2

    copy_bc = .true.
    if (copy_bc .and. present(exclude_flux_type)) &
      copy_bc = .not.(trim(var_in%bc(n)%flux_type) == trim(exclude_flux_type))
    if (copy_bc .and. present(only_flux_type)) &
      copy_bc = (trim(var_in%bc(n)%flux_type) == trim(only_flux_type))
    if (copy_bc .and. present(pass_through_ice)) &
      copy_bc = (pass_through_ice .eqv. var_in%bc(n)%pass_through_ice)
    if (.not.copy_bc) cycle

    do m = 1, var_in%bc(n)%num_fields
      if (present(field_index)) then ; if (m /= field_index) cycle ; endif
      if ( associated(var%bc(n)%field(m)%values) ) then
        ks = var%ks ; if (present(ind3_start)) ks = max(ks, ind3_start)
        ke = var%ke ; if (present(ind3_end)) ke = max(ke, ind3_end)
        do k=ks,ke ; do j=var%jsc-halo,var%jec+halo ; do i=var%isc-halo,var%iec+halo
          var%bc(n)%field(m)%values(i,j,k) = var_in%bc(n)%field(m)%values(i+i_off,j+j_off)
        enddo ; enddo ; enddo
      endif
    enddo
  enddo

end subroutine CT_copy_data_2d_3d


!> This subroutine redistributes the data in all elements of one coupler_2d_bc_type
!! into another, which may be on different processors with a different decomposition.
subroutine CT_redistribute_data_2d(var_in, domain_in, var_out, domain_out, complete)
  type(coupler_2d_bc_type), intent(in)    :: var_in     !< BC_type structure with the data to copy (intent in)
  type(domain2D),           intent(in)    :: domain_in  !< The FMS domain for the input structure
  type(coupler_2d_bc_type), intent(inout) :: var_out    !< The recipient BC_type structure (data intent out)
  type(domain2D),           intent(in)    :: domain_out !< The FMS domain for the output structure
  logical,        optional, intent(in)    :: complete   !< If true, complete the updates

  real, pointer, dimension(:,:) :: null_ptr2D => NULL()
  logical :: do_in, do_out, do_complete
  integer :: m, n, fc, fc_in, fc_out

  do_complete = .true. ; if (present(complete)) do_complete = complete

  ! Figure out whether this PE has valid input or output fields or both.
  do_in = var_in%set
  do_out = var_out%set

  fc_in = 0 ; fc_out = 0
  if (do_in) then ; do n = 1, var_in%num_bcs ; do m = 1, var_in%bc(n)%num_fields
    if (associated(var_in%bc(n)%field(m)%values)) fc_in = fc_in + 1
  enddo ; enddo ; endif
  if (fc_in == 0) do_in = .false.
  if (do_out) then ; do n = 1, var_out%num_bcs ; do m = 1, var_out%bc(n)%num_fields
    if (associated(var_out%bc(n)%field(m)%values)) fc_out = fc_out + 1
  enddo ; enddo ; endif
  if (fc_out == 0) do_out = .false.

  if (do_in .and. do_out) then
    if (var_in%num_bcs /= var_out%num_bcs) call mpp_error(FATAL, &
      "Mismatch in num_bcs in CT_copy_data_2d.")
    if (fc_in /= fc_out) call mpp_error(FATAL, &
      "Mismatch in the total number of fields in CT_redistribute_data_2d.")
  endif

  if (.not.(do_in .or. do_out)) return

  fc = 0
  if (do_in .and. do_out) then
    do n = 1, var_in%num_bcs ; do m = 1, var_in%bc(n)%num_fields
      if ( associated(var_in%bc(n)%field(m)%values) .neqv. &
           associated(var_out%bc(n)%field(m)%values) ) &
        call mpp_error(FATAL, &
          "Mismatch in which fields are associated in CT_redistribute_data_2d.")

      if ( associated(var_in%bc(n)%field(m)%values) ) then
        fc = fc + 1
        call mpp_redistribute(domain_in, var_in%bc(n)%field(m)%values, &
                              domain_out, var_out%bc(n)%field(m)%values, &
                              complete=(do_complete.and.(fc==fc_in)) )
      endif
    enddo ; enddo
  elseif (do_in) then
    do n = 1, var_in%num_bcs ; do m = 1, var_in%bc(n)%num_fields
      if ( associated(var_in%bc(n)%field(m)%values) ) then
        fc = fc + 1
        call mpp_redistribute(domain_in, var_in%bc(n)%field(m)%values, &
                              domain_out, null_ptr2D, &
                              complete=(do_complete.and.(fc==fc_in)) )
      endif
    enddo ; enddo
  elseif (do_out) then
    do n = 1, var_out%num_bcs ; do m = 1, var_out%bc(n)%num_fields
      if ( associated(var_out%bc(n)%field(m)%values) ) then
        fc = fc + 1
        call mpp_redistribute(domain_in, null_ptr2D, &
                              domain_out, var_out%bc(n)%field(m)%values, &
                              complete=(do_complete.and.(fc==fc_out)) )
      endif
    enddo ; enddo
  endif

end subroutine CT_redistribute_data_2d

!> This subroutine redistributes the data in all elements of one coupler_2d_bc_type
!! into another, which may be on different processors with a different decomposition.
subroutine CT_redistribute_data_3d(var_in, domain_in, var_out, domain_out, complete)
  type(coupler_3d_bc_type), intent(in)    :: var_in     !< BC_type structure with the data to copy (intent in)
  type(domain2D),           intent(in)    :: domain_in  !< The FMS domain for the input structure
  type(coupler_3d_bc_type), intent(inout) :: var_out    !< The recipient BC_type structure (data intent out)
  type(domain2D),           intent(in)    :: domain_out !< The FMS domain for the output structure
  logical,        optional, intent(in)    :: complete   !< If true, complete the updates

  real, pointer, dimension(:,:,:) :: null_ptr3D => NULL()
  logical :: do_in, do_out, do_complete
  integer :: m, n, fc, fc_in, fc_out

  do_complete = .true. ; if (present(complete)) do_complete = complete

  ! Figure out whether this PE has valid input or output fields or both.
  do_in = var_in%set
  do_out = var_out%set

  fc_in = 0 ; fc_out = 0
  if (do_in) then ; do n = 1, var_in%num_bcs ; do m = 1, var_in%bc(n)%num_fields
    if (associated(var_in%bc(n)%field(m)%values)) fc_in = fc_in + 1
  enddo ; enddo ; endif
  if (fc_in == 0) do_in = .false.
  if (do_out) then ; do n = 1, var_out%num_bcs ; do m = 1, var_out%bc(n)%num_fields
    if (associated(var_out%bc(n)%field(m)%values)) fc_out = fc_out + 1
  enddo ; enddo ; endif
  if (fc_out == 0) do_out = .false.

  if (do_in .and. do_out) then
    if (var_in%num_bcs /= var_out%num_bcs) call mpp_error(FATAL, &
      "Mismatch in num_bcs in CT_copy_data_3d.")
    if (fc_in /= fc_out) call mpp_error(FATAL, &
      "Mismatch in the total number of fields in CT_redistribute_data_3d.")
  endif

  if (.not.(do_in .or. do_out)) return

  fc = 0
  if (do_in .and. do_out) then
    do n = 1, var_in%num_bcs ; do m = 1, var_in%bc(n)%num_fields
      if ( associated(var_in%bc(n)%field(m)%values) .neqv. &
           associated(var_out%bc(n)%field(m)%values) ) &
        call mpp_error(FATAL, &
          "Mismatch in which fields are associated in CT_redistribute_data_3d.")

      if ( associated(var_in%bc(n)%field(m)%values) ) then
        fc = fc + 1
        call mpp_redistribute(domain_in, var_in%bc(n)%field(m)%values, &
                              domain_out, var_out%bc(n)%field(m)%values, &
                              complete=(do_complete.and.(fc==fc_in)) )
      endif
    enddo ; enddo
  elseif (do_in) then
    do n = 1, var_in%num_bcs ; do m = 1, var_in%bc(n)%num_fields
      if ( associated(var_in%bc(n)%field(m)%values) ) then
        fc = fc + 1
        call mpp_redistribute(domain_in, var_in%bc(n)%field(m)%values, &
                              domain_out, null_ptr3D, &
                              complete=(do_complete.and.(fc==fc_in)) )
      endif
    enddo ; enddo
  elseif (do_out) then
    do n = 1, var_out%num_bcs ; do m = 1, var_out%bc(n)%num_fields
      if ( associated(var_out%bc(n)%field(m)%values) ) then
        fc = fc + 1
        call mpp_redistribute(domain_in, null_ptr3D, &
                              domain_out, var_out%bc(n)%field(m)%values, &
                              complete=(do_complete.and.(fc==fc_out)) )
      endif
    enddo ; enddo
  endif

end subroutine CT_redistribute_data_3d


!> This subroutine rescales the fields in the elements of a coupler_2d_bc_type
!! by multiplying by a factor scale.  If scale is 0, this is a direct
!! assignment to 0, so that NaNs will not persist.
subroutine CT_rescale_data_2d(var, scale, halo_size, bc_index, field_index, &
                           exclude_flux_type, only_flux_type, pass_through_ice)
  type(coupler_2d_bc_type),   intent(inout) :: var !< The BC_type structure whose fields are being rescaled
  real,                       intent(in)    :: scale   !< A scaling factor to multiply fields by
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default or
                                                         !! the full arrays if scale is 0.
  integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                       !! that is being copied
  integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                       !! boundary condition that is being copied
  character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types of fluxes to exclude from this copy.
  character(len=*), optional, intent(in)    :: only_flux_type    !< A string describing which types of fluxes to include from this copy.
  logical,          optional, intent(in)    :: pass_through_ice !< If true, only copy BCs whose
                                                       !! value of pass_through ice matches this
  logical :: do_bc
  integer :: i, j, m, n, n1, n2, halo

  if (present(bc_index)) then
    if (bc_index > var%num_bcs) &
      call mpp_error(FATAL, "CT_rescale_data_2d: bc_index is present and exceeds var%num_bcs.")
    if (present(field_index)) then ; if (field_index > var%bc(bc_index)%num_fields) &
      call mpp_error(FATAL, "CT_rescale_data_2d: field_index is present and exceeds num_fields for" //&
                     trim(var%bc(bc_index)%name) )
    endif
  elseif (present(field_index)) then
    call mpp_error(FATAL, "CT_rescale_data_2d: bc_index must be present if field_index is present.")
  endif

  halo = 0 ; if (present(halo_size)) halo = halo_size

  n1 = 1 ; n2 = var%num_bcs
  if (present(bc_index)) then ; n1 = bc_index ; n2 = bc_index ; endif

  if (n2 >= n1) then
    ! A more consciencious implementation would include a more descriptive error messages.
    if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo)) &
      call mpp_error(FATAL, "CT_rescale_data_2d: Excessive i-direction halo size.")
    if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo)) &
      call mpp_error(FATAL, "CT_rescale_data_2d: Excessive j-direction halo size.")
  endif

  do n = n1, n2

    do_bc = .true.
    if (do_bc .and. present(exclude_flux_type)) &
      do_bc = .not.(trim(var%bc(n)%flux_type) == trim(exclude_flux_type))
    if (do_bc .and. present(only_flux_type)) &
      do_bc = (trim(var%bc(n)%flux_type) == trim(only_flux_type))
    if (do_bc .and. present(pass_through_ice)) &
      do_bc = (pass_through_ice .eqv. var%bc(n)%pass_through_ice)
    if (.not.do_bc) cycle

    do m = 1, var%bc(n)%num_fields
      if (present(field_index)) then ; if (m /= field_index) cycle ; endif
      if ( associated(var%bc(n)%field(m)%values) ) then
        if (scale == 0.0) then
          if (present(halo_size)) then
            do j=var%jsc-halo,var%jec+halo ; do i=var%isc-halo,var%iec+halo
              var%bc(n)%field(m)%values(i,j) = 0.0
            enddo ; enddo
          else
            var%bc(n)%field(m)%values(:,:) = 0.0
          endif
        else
          do j=var%jsc-halo,var%jec+halo ; do i=var%isc-halo,var%iec+halo
            var%bc(n)%field(m)%values(i,j) = scale * var%bc(n)%field(m)%values(i,j)
          enddo ; enddo
        endif
      endif
    enddo
  enddo

end subroutine CT_rescale_data_2d

!> This subroutine rescales the fields in the elements of a coupler_3d_bc_type
!! by multiplying by a factor scale.  If scale is 0, this is a direct
!! assignment to 0, so that NaNs will not persist.
subroutine CT_rescale_data_3d(var, scale, halo_size, bc_index, field_index, &
                           exclude_flux_type, only_flux_type, pass_through_ice)
  type(coupler_3d_bc_type),   intent(inout) :: var !< The BC_type structure whose fields are being rescaled
  real,                       intent(in)    :: scale   !< A scaling factor to multiply fields by
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default or
                                                         !! the full arrays if scale is 0.
  integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                       !! that is being copied
  integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                       !! boundary condition that is being copied
  character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types of fluxes to exclude from this copy.
  character(len=*), optional, intent(in)    :: only_flux_type    !< A string describing which types of fluxes to include from this copy.
  logical,          optional, intent(in)    :: pass_through_ice !< If true, only copy BCs whose
                                                       !! value of pass_through ice matches this
  logical :: do_bc
  integer :: i, j, k, m, n, n1, n2, halo

  if (present(bc_index)) then
    if (bc_index > var%num_bcs) &
      call mpp_error(FATAL, "CT_rescale_data_2d: bc_index is present and exceeds var%num_bcs.")
    if (present(field_index)) then ; if (field_index > var%bc(bc_index)%num_fields) &
      call mpp_error(FATAL, "CT_rescale_data_2d: field_index is present and exceeds num_fields for" //&
                     trim(var%bc(bc_index)%name) )
    endif
  elseif (present(field_index)) then
    call mpp_error(FATAL, "CT_rescale_data_2d: bc_index must be present if field_index is present.")
  endif

  halo = 0 ; if (present(halo_size)) halo = halo_size

  n1 = 1 ; n2 = var%num_bcs
  if (present(bc_index)) then ; n1 = bc_index ; n2 = bc_index ; endif

  if (n2 >= n1) then
    ! A more consciencious implementation would include a more descriptive error messages.
    if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo)) &
      call mpp_error(FATAL, "CT_rescale_data_3d: Excessive i-direction halo size.")
    if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo)) &
      call mpp_error(FATAL, "CT_rescale_data_3d: Excessive j-direction halo size.")
  endif

  do n = n1, n2

    do_bc = .true.
    if (do_bc .and. present(exclude_flux_type)) &
      do_bc = .not.(trim(var%bc(n)%flux_type) == trim(exclude_flux_type))
    if (do_bc .and. present(only_flux_type)) &
      do_bc = (trim(var%bc(n)%flux_type) == trim(only_flux_type))
    if (do_bc .and. present(pass_through_ice)) &
      do_bc = (pass_through_ice .eqv. var%bc(n)%pass_through_ice)
    if (.not.do_bc) cycle

    do m = 1, var%bc(n)%num_fields
      if (present(field_index)) then ; if (m /= field_index) cycle ; endif
      if ( associated(var%bc(n)%field(m)%values) ) then
        if (scale == 0.0) then
          if (present(halo_size)) then
            do k=var%ks,var%ke ; do j=var%jsc-halo,var%jec+halo ; do i=var%isc-halo,var%iec+halo
              var%bc(n)%field(m)%values(i,j,k) = 0.0
            enddo ; enddo ; enddo
          else
            var%bc(n)%field(m)%values(:,:,:) = 0.0
          endif
        else
          do k=var%ks,var%ke ; do j=var%jsc-halo,var%jec+halo ; do i=var%isc-halo,var%iec+halo
            var%bc(n)%field(m)%values(i,j,k) = scale * var%bc(n)%field(m)%values(i,j,k)
          enddo ; enddo ; enddo
        endif
      endif
    enddo
  enddo

end subroutine CT_rescale_data_3d


!> This subroutine does a direct increment of the data in all elements of one
!! coupler_2d_bc_type into another.  Both must have the same array sizes.
subroutine CT_increment_data_2d_2d(var_in, var, halo_size, bc_index, field_index, &
                           scale_factor, scale_prev, exclude_flux_type, only_flux_type, pass_through_ice)
  type(coupler_2d_bc_type),   intent(in)    :: var_in  !< BC_type structure with the data to add to the other type
  type(coupler_2d_bc_type),   intent(inout) :: var !< The BC_type structure whose fields are being incremented
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to increment; 0 by default
  integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                       !! that is being copied
  integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                       !! boundary condition that is being copied
  real,             optional, intent(in)    :: scale_factor  !< A scaling factor for the data that is being added
  real,             optional, intent(in)    :: scale_prev    !< A scaling factor for the data that is already here
  character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types of fluxes to exclude from this increment.
  character(len=*), optional, intent(in)    :: only_flux_type    !< A string describing which types of fluxes to include from this increment.
  logical,          optional, intent(in)    :: pass_through_ice !< If true, only increment BCs whose
                                                       !! value of pass_through ice matches this

  real :: scale, sc_prev
  logical :: increment_bc
  integer :: i, j, m, n, n1, n2, halo, i_off, j_off

  scale = 1.0 ; if (present(scale_factor)) scale = scale_factor
  sc_prev = 1.0 ; if (present(scale_prev)) sc_prev = scale_prev

  if (present(bc_index)) then
    if (bc_index > var_in%num_bcs) &
      call mpp_error(FATAL, "CT_increment_data_2d_2d: bc_index is present and exceeds var_in%num_bcs.")
    if (present(field_index)) then ; if (field_index > var_in%bc(bc_index)%num_fields) &
      call mpp_error(FATAL, "CT_increment_data_2d_2d: field_index is present and exceeds num_fields for" //&
                     trim(var_in%bc(bc_index)%name) )
    endif
  elseif (present(field_index)) then
    call mpp_error(FATAL, "CT_increment_data_2d_2d: bc_index must be present if field_index is present.")
  endif

  halo = 0 ; if (present(halo_size)) halo = halo_size

  n1 = 1 ; n2 = var_in%num_bcs
  if (present(bc_index)) then ; n1 = bc_index ; n2 = bc_index ; endif

  if (n2 >= n1) then
    ! A more consciencious implementation would include a more descriptive error messages.
    if ((var_in%iec-var_in%isc) /= (var%iec-var%isc)) &
      call mpp_error(FATAL, "CT_increment_data_2d: There is an i-direction computional domain size mismatch.")
    if ((var_in%jec-var_in%jsc) /= (var%jec-var%jsc)) &
      call mpp_error(FATAL, "CT_increment_data_2d: There is a j-direction computional domain size mismatch.")
    if ((var_in%isc-var_in%isd < halo) .or. (var_in%ied-var_in%iec < halo)) &
      call mpp_error(FATAL, "CT_increment_data_2d: Excessive i-direction halo size for the input structure.")
    if ((var_in%jsc-var_in%jsd < halo) .or. (var_in%jed-var_in%jec < halo)) &
      call mpp_error(FATAL, "CT_increment_data_2d: Excessive j-direction halo size for the input structure.")
    if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo)) &
      call mpp_error(FATAL, "CT_increment_data_2d: Excessive i-direction halo size for the output structure.")
    if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo)) &
      call mpp_error(FATAL, "CT_increment_data_2d: Excessive j-direction halo size for the output structure.")

    i_off = var_in%isc - var%isc ; j_off = var_in%jsc - var%jsc
  endif

  do n = n1, n2

    increment_bc = .true.
    if (increment_bc .and. present(exclude_flux_type)) &
      increment_bc = .not.(trim(var%bc(n)%flux_type) == trim(exclude_flux_type))
    if (increment_bc .and. present(only_flux_type)) &
      increment_bc = (trim(var%bc(n)%flux_type) == trim(only_flux_type))
    if (increment_bc .and. present(pass_through_ice)) &
      increment_bc = (pass_through_ice .eqv. var%bc(n)%pass_through_ice)
    if (.not.increment_bc) cycle

    do m = 1, var_in%bc(n)%num_fields
      if (present(field_index)) then ; if (m /= field_index) cycle ; endif
      if ( associated(var%bc(n)%field(m)%values) ) then
        do j=var%jsc-halo,var%jec+halo ; do i=var%isc-halo,var%iec+halo
          var%bc(n)%field(m)%values(i,j) = sc_prev * var%bc(n)%field(m)%values(i,j) + &
                          scale * var_in%bc(n)%field(m)%values(i+i_off,j+j_off)
        enddo ; enddo
      endif
    enddo
  enddo

end subroutine CT_increment_data_2d_2d


!> This subroutine does a direct increment of the data in all elements of one
!! coupler_3d_bc_type into another.  Both must have the same array sizes.
subroutine CT_increment_data_3d_3d(var_in, var, halo_size, bc_index, field_index, &
                           scale_factor, scale_prev, exclude_flux_type, only_flux_type, pass_through_ice)
  type(coupler_3d_bc_type),   intent(in)    :: var_in  !< BC_type structure with the data to add to the other type
  type(coupler_3d_bc_type),   intent(inout) :: var !< The BC_type structure whose fields are being incremented
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                       !! that is being copied
  integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                       !! boundary condition that is being copied
  real,             optional, intent(in)    :: scale_factor  !< A scaling factor for the data that is being added
  real,             optional, intent(in)    :: scale_prev    !< A scaling factor for the data that is already here
  character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types of fluxes to exclude from this increment.
  character(len=*), optional, intent(in)    :: only_flux_type    !< A string describing which types of fluxes to include from this increment.
  logical,          optional, intent(in)    :: pass_through_ice !< If true, only increment BCs whose
                                                       !! value of pass_through ice matches this

  real :: scale, sc_prev
  logical :: increment_bc
  integer :: i, j, k, m, n, n1, n2, halo, i_off, j_off, k_off

  scale = 1.0 ; if (present(scale_factor)) scale = scale_factor
  sc_prev = 1.0 ; if (present(scale_prev)) sc_prev = scale_prev

  if (present(bc_index)) then
    if (bc_index > var_in%num_bcs) &
      call mpp_error(FATAL, "CT_increment_data_3d_3d: bc_index is present and exceeds var_in%num_bcs.")
    if (present(field_index)) then ; if (field_index > var_in%bc(bc_index)%num_fields) &
      call mpp_error(FATAL, "CT_increment_data_3d_3d: field_index is present and exceeds num_fields for" //&
                     trim(var_in%bc(bc_index)%name) )
    endif
  elseif (present(field_index)) then
    call mpp_error(FATAL, "CT_increment_data_3d_3d: bc_index must be present if field_index is present.")
  endif

  halo = 0 ; if (present(halo_size)) halo = halo_size

  n1 = 1 ; n2 = var_in%num_bcs
  if (present(bc_index)) then ; n1 = bc_index ; n2 = bc_index ; endif

  if (n2 >= n1) then
    ! A more consciencious implementation would include a more descriptive error messages.
    if ((var_in%iec-var_in%isc) /= (var%iec-var%isc)) &
      call mpp_error(FATAL, "CT_increment_data_3d: There is an i-direction computional domain size mismatch.")
    if ((var_in%jec-var_in%jsc) /= (var%jec-var%jsc)) &
      call mpp_error(FATAL, "CT_increment_data_3d: There is a j-direction computional domain size mismatch.")
    if ((var_in%ke-var_in%ks) /= (var%ke-var%ks)) &
      call mpp_error(FATAL, "CT_increment_data_3d: There is a k-direction computional domain size mismatch.")
    if ((var_in%isc-var_in%isd < halo) .or. (var_in%ied-var_in%iec < halo)) &
      call mpp_error(FATAL, "CT_increment_data_3d: Excessive i-direction halo size for the input structure.")
    if ((var_in%jsc-var_in%jsd < halo) .or. (var_in%jed-var_in%jec < halo)) &
      call mpp_error(FATAL, "CT_increment_data_3d: Excessive j-direction halo size for the input structure.")
    if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo)) &
      call mpp_error(FATAL, "CT_increment_data_3d: Excessive i-direction halo size for the output structure.")
    if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo)) &
      call mpp_error(FATAL, "CT_increment_data_3d: Excessive j-direction halo size for the output structure.")

    i_off = var_in%isc - var%isc ; j_off = var_in%jsc - var%jsc ; k_off = var_in%ks - var%ks
  endif

  do n = n1, n2

    increment_bc = .true.
    if (increment_bc .and. present(exclude_flux_type)) &
      increment_bc = .not.(trim(var%bc(n)%flux_type) == trim(exclude_flux_type))
    if (increment_bc .and. present(only_flux_type)) &
      increment_bc = (trim(var%bc(n)%flux_type) == trim(only_flux_type))
    if (increment_bc .and. present(pass_through_ice)) &
      increment_bc = (pass_through_ice .eqv. var%bc(n)%pass_through_ice)
    if (.not.increment_bc) cycle

    do m = 1, var_in%bc(n)%num_fields
      if (present(field_index)) then ; if (m /= field_index) cycle ; endif
      if ( associated(var%bc(n)%field(m)%values) ) then
        do k=var%ks,var%ke ; do j=var%jsc-halo,var%jec+halo ; do i=var%isc-halo,var%iec+halo
          var%bc(n)%field(m)%values(i,j,k) = sc_prev * var%bc(n)%field(m)%values(i,j,k) + &
                       scale * var_in%bc(n)%field(m)%values(i+i_off,j+j_off,k+k_off)
        enddo ; enddo ; enddo
      endif
    enddo
  enddo

end subroutine CT_increment_data_3d_3d

!> This subroutine does increments the data in the elements of a coupler_2d_bc_type
!! with the weighed average of the elements of a coupler_3d_bc_type. Both must have
!! the same horizontal array sizes and the normalized weight array must match the
!! array sizes of the coupler_3d_bc_type.
subroutine CT_increment_data_2d_3d(var_in, weights, var, halo_size, bc_index, field_index, &
                           scale_factor, scale_prev, exclude_flux_type, only_flux_type, pass_through_ice)
  type(coupler_3d_bc_type),   intent(in)    :: var_in  !< BC_type structure with the data to add to the other type
  real, dimension(:,:,:),     intent(in)    :: weights !< An array of normalized weights for the 3d-data to
                                                       !! increment the 2d-data.  There is no renormalization,
                                                       !! so if the weights do not sum to 1 in the 3rd dimension
                                                       !! there may be adverse consequences!
  type(coupler_2d_bc_type),   intent(inout) :: var !< The BC_type structure whose fields are being incremented
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                       !! that is being copied
  integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                       !! boundary condition that is being copied
  real,             optional, intent(in)    :: scale_factor  !< A scaling factor for the data that is being added
  real,             optional, intent(in)    :: scale_prev    !< A scaling factor for the data that is already here
  character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types of fluxes to exclude from this increment.
  character(len=*), optional, intent(in)    :: only_flux_type    !< A string describing which types of fluxes to include from this increment.
  logical,          optional, intent(in)    :: pass_through_ice !< If true, only increment BCs whose
                                                       !! value of pass_through ice matches this

  real :: scale, sc_prev
  logical :: increment_bc
  integer :: i, j, k, m, n, n1, n2, halo
  integer :: io1, jo1, iow, jow, kow  ! Offsets to account for different index conventions.

  scale = 1.0 ; if (present(scale_factor)) scale = scale_factor
  sc_prev = 1.0 ; if (present(scale_prev)) sc_prev = scale_prev

  if (present(bc_index)) then
    if (bc_index > var_in%num_bcs) &
      call mpp_error(FATAL, "CT_increment_data_2d_3d: bc_index is present and exceeds var_in%num_bcs.")
    if (present(field_index)) then ; if (field_index > var_in%bc(bc_index)%num_fields) &
      call mpp_error(FATAL, "CT_increment_data_2d_3d: field_index is present and exceeds num_fields for" //&
                     trim(var_in%bc(bc_index)%name) )
    endif
  elseif (present(field_index)) then
    call mpp_error(FATAL, "CT_increment_data_2d_3d: bc_index must be present if field_index is present.")
  endif

  halo = 0 ; if (present(halo_size)) halo = halo_size

  n1 = 1 ; n2 = var_in%num_bcs
  if (present(bc_index)) then ; n1 = bc_index ; n2 = bc_index ; endif

  if (n2 >= n1) then
    ! A more consciencious implementation would include a more descriptive error messages.
    if ((var_in%iec-var_in%isc) /= (var%iec-var%isc)) &
      call mpp_error(FATAL, "CT_increment_data_2d_3d: There is an i-direction computional domain size mismatch.")
    if ((var_in%jec-var_in%jsc) /= (var%jec-var%jsc)) &
      call mpp_error(FATAL, "CT_increment_data_2d_3d: There is a j-direction computional domain size mismatch.")
    if ((1+var_in%ke-var_in%ks) /= size(weights,3)) &
      call mpp_error(FATAL, "CT_increment_data_2d_3d: There is a k-direction size mismatch with the weights array.")
    if ((var_in%isc-var_in%isd < halo) .or. (var_in%ied-var_in%iec < halo)) &
      call mpp_error(FATAL, "CT_increment_data_2d_3d: Excessive i-direction halo size for the input structure.")
    if ((var_in%jsc-var_in%jsd < halo) .or. (var_in%jed-var_in%jec < halo)) &
      call mpp_error(FATAL, "CT_increment_data_2d_3d: Excessive j-direction halo size for the input structure.")
    if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo)) &
      call mpp_error(FATAL, "CT_increment_data_2d_3d: Excessive i-direction halo size for the output structure.")
    if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo)) &
      call mpp_error(FATAL, "CT_increment_data_2d_3d: Excessive j-direction halo size for the output structure.")

    if ((1+var%iec-var%isc) == size(weights,1)) then
      iow = 1 - var%isc
    elseif ((1+var%ied-var%isd) == size(weights,1)) then
      iow = 1 - var%isd
    elseif ((1+var_in%ied-var_in%isd) == size(weights,1)) then
      iow = 1 + (var_in%isc - var_in%isd) - var%isc
    else
      call mpp_error(FATAL, "CT_increment_data_2d_3d: weights array must be the i-size of a computational or data domain.")
    endif
    if ((1+var%jec-var%jsc) == size(weights,2)) then
      jow = 1 - var%jsc
    elseif ((1+var%jed-var%jsd) == size(weights,2)) then
      jow = 1 - var%jsd
    elseif ((1+var_in%jed-var_in%jsd) == size(weights,2)) then
      jow = 1 + (var_in%jsc - var_in%jsd) - var%jsc
    else
      call mpp_error(FATAL, "CT_increment_data_2d_3d: weights array must be the j-size of a computational or data domain.")
    endif

    io1 = var_in%isc - var%isc ; jo1 = var_in%jsc - var%jsc ; kow = 1 - var_in%ks
  endif

  do n = n1, n2

    increment_bc = .true.
    if (increment_bc .and. present(exclude_flux_type)) &
      increment_bc = .not.(trim(var_in%bc(n)%flux_type) == trim(exclude_flux_type))
    if (increment_bc .and. present(only_flux_type)) &
      increment_bc = (trim(var_in%bc(n)%flux_type) == trim(only_flux_type))
    if (increment_bc .and. present(pass_through_ice)) &
      increment_bc = (pass_through_ice .eqv. var_in%bc(n)%pass_through_ice)
    if (.not.increment_bc) cycle

    do m = 1, var_in%bc(n)%num_fields
      if (present(field_index)) then ; if (m /= field_index) cycle ; endif
      if ( associated(var%bc(n)%field(m)%values) ) then
        do k=var_in%ks,var_in%ke ; do j=var%jsc-halo,var%jec+halo ; do i=var%isc-halo,var%iec+halo
          var%bc(n)%field(m)%values(i,j) = sc_prev * var%bc(n)%field(m)%values(i,j) + &
                     (scale * weights(i+iow,j+jow,k+kow)) * var_in%bc(n)%field(m)%values(i+io1,j+io1,k)
        enddo ; enddo ; enddo
      endif
    enddo
  enddo

end subroutine CT_increment_data_2d_3d


!> This subroutine extracts a single 2-d field from a coupler_2d_bc_type into
!! a two-dimensional array.
subroutine CT_extract_data_2d(var_in, bc_index, field_index, array_out, &
                              scale_factor, halo_size, idim, jdim)
  type(coupler_2d_bc_type),   intent(in)    :: var_in    !< BC_type structure with the data to extract
  integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                         !! that is being copied
  integer,                    intent(in)    :: field_index !< The index of the field in the
                                                         !! boundary condition that is being copied
  real, dimension(1:,1:),     intent(out)   :: array_out !< The recipient array for the field; its size
                                                         !! must match the size of the data being copied
                                                         !! unless idim and jdim are supplied.
  real,             optional, intent(in)    :: scale_factor !< A scaling factor for the data that is being added
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer, dimension(4), optional, intent(in) :: idim    !< The data and computational domain extents of
                                                         !! the first dimension of the output array
                                                         !! in a non-decreasing list
  integer, dimension(4), optional, intent(in) :: jdim    !< The data and computational domain extents of
                                                         !! the second dimension of the output array
                                                         !! in a non-decreasing list
  character(len=256), parameter :: error_header = &
     '==>Error from coupler_types_mod (CT_extract_data_2d):'
  character(len=400)      :: error_msg

  real :: scale
  integer :: i, j, halo, i_off, j_off

  if (bc_index <= 0) then
    array_out(:,:) = 0.0
    return
  endif

  halo = 0 ; if (present(halo_size)) halo = halo_size
  scale = 1.0 ; if (present(scale_factor)) scale = scale_factor

  if ((var_in%isc-var_in%isd < halo) .or. (var_in%ied-var_in%iec < halo)) &
    call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the input structure.")
  if ((var_in%jsc-var_in%jsd < halo) .or. (var_in%jed-var_in%jec < halo)) &
    call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the input structure.")

  if (bc_index > var_in%num_bcs) &
    call mpp_error(FATAL, trim(error_header)//" bc_index exceeds var_in%num_bcs.")
  if (field_index > var_in%bc(bc_index)%num_fields) &
    call mpp_error(FATAL, trim(error_header)//" field_index exceeds num_fields for" //&
                   trim(var_in%bc(bc_index)%name) )

  ! Do error checking on the i-dimension and determine the array offsets.
  if (present(idim)) then
    if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
      write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
      call mpp_error(FATAL, trim(error_msg))
    endif
    if (size(array_out,1) /= (1+idim(4)-idim(1))) then
      write (error_msg, *) trim(error_header), ' The declared i-dimension size of ', &
            (1+idim(4)-idim(1)), ' does not match the actual size of ', size(array_out,1)
      call mpp_error(FATAL, trim(error_msg))
    endif
    if ((var_in%iec-var_in%isc) /= (idim(3)-idim(2))) &
      call mpp_error(FATAL, trim(error_header)//" There is an i-direction computional domain size mismatch.")
    if ((idim(2)-idim(1) < halo) .or. (idim(4)-idim(3) < halo)) &
      call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the output array.")
    if (size(array_out,1) < 2*halo + 1 + var_in%iec - var_in%isc) then
      write (error_msg, *) trim(error_header), ' The target array with i-dimension size ', &
            (1+idim(4)-idim(1)), ' is too small to match the data of size ', &
            (2*halo + 1 + var_in%iec - var_in%isc)
      call mpp_error(FATAL, trim(error_msg))
    endif

    i_off = (1-idim(1)) + (idim(2)-var_in%isc)
  else
    if (size(array_out,1) < 2*halo + 1 + var_in%iec - var_in%isc) then
      write (error_msg, *) trim(error_header), ' The target array with i-dimension size ', &
            size(array_out,1), ' does not match the data of size ', &
            (2*halo + 1 + var_in%iec - var_in%isc)
      call mpp_error(FATAL, trim(error_msg))
    endif
    i_off = 1 - (var_in%isc-halo)
  endif

  ! Do error checking on the j-dimension and determine the array offsets.
  if (present(jdim)) then
    if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
      write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list ', jdim
      call mpp_error(FATAL, trim(error_msg))
    endif
    if (size(array_out,2) /= (1+jdim(4)-jdim(1))) then
      write (error_msg, *) trim(error_header), ' The declared j-dimension size of ', &
            (1+jdim(4)-jdim(1)), ' does not match the actual size of ', size(array_out,2)
      call mpp_error(FATAL, trim(error_msg))
    endif
    if ((var_in%jec-var_in%jsc) /= (jdim(3)-jdim(2))) &
      call mpp_error(FATAL, trim(error_header)//" There is an j-direction computional domain size mismatch.")
    if ((jdim(2)-jdim(1) < halo) .or. (jdim(4)-jdim(3) < halo)) &
      call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the output array.")
    if (size(array_out,2) < 2*halo + 1 + var_in%jec - var_in%jsc) then
      write (error_msg, *) trim(error_header), ' The target array with j-dimension size ', &
            (1+jdim(4)-jdim(1)), ' is too small to match the data of size ', &
            (2*halo + 1 + var_in%jec - var_in%jsc)
      call mpp_error(FATAL, trim(error_msg))
    endif

    j_off = (1-jdim(1)) + (jdim(2)-var_in%jsc)
  else
    if (size(array_out,2) < 2*halo + 1 + var_in%jec - var_in%jsc) then
      write (error_msg, *) trim(error_header), ' The target array with j-dimension size ', &
            size(array_out,2), ' does not match the data of size ', &
            (2*halo + 1 + var_in%jec - var_in%jsc)
      call mpp_error(FATAL, trim(error_msg))
    endif
    j_off = 1 - (var_in%jsc-halo)
  endif

  do j=var_in%jsc-halo,var_in%jec+halo ; do i=var_in%isc-halo,var_in%iec+halo
    array_out(i+i_off,j+j_off) = scale * var_in%bc(bc_index)%field(field_index)%values(i,j)
  enddo ; enddo

end subroutine CT_extract_data_2d

!> This subroutine extracts a single k-level of a 3-d field from a coupler_3d_bc_type
!! into a two-dimensional array.
subroutine CT_extract_data_3d_2d(var_in, bc_index, field_index, k_in, array_out, &
                                 scale_factor, halo_size, idim, jdim)
  type(coupler_3d_bc_type),   intent(in)    :: var_in    !< BC_type structure with the data to extract
  integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                         !! that is being copied
  integer,                    intent(in)    :: field_index !< The index of the field in the
                                                         !! boundary condition that is being copied
  integer,                    intent(in)    :: k_in      !< The k-index to extract
  real, dimension(1:,1:),     intent(out)   :: array_out !< The recipient array for the field; its size
                                                         !! must match the size of the data being copied
                                                         !! unless idim and jdim are supplied.
  real,             optional, intent(in)    :: scale_factor !< A scaling factor for the data that is being added
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer, dimension(4), optional, intent(in) :: idim    !< The data and computational domain extents of
                                                         !! the first dimension of the output array
                                                         !! in a non-decreasing list
  integer, dimension(4), optional, intent(in) :: jdim    !< The data and computational domain extents of
                                                         !! the second dimension of the output array
                                                         !! in a non-decreasing list
  character(len=256), parameter :: error_header = &
     '==>Error from coupler_types_mod (CT_extract_data_3d_2d):'
  character(len=400)      :: error_msg

  real :: scale
  integer :: i, j, k, halo, i_off, j_off

  if (bc_index <= 0) then
    array_out(:,:) = 0.0
    return
  endif

  halo = 0 ; if (present(halo_size)) halo = halo_size
  scale = 1.0 ; if (present(scale_factor)) scale = scale_factor

  if ((var_in%isc-var_in%isd < halo) .or. (var_in%ied-var_in%iec < halo)) &
    call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the input structure.")
  if ((var_in%jsc-var_in%jsd < halo) .or. (var_in%jed-var_in%jec < halo)) &
    call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the input structure.")

  if (bc_index > var_in%num_bcs) &
    call mpp_error(FATAL, trim(error_header)//" bc_index exceeds var_in%num_bcs.")
  if (field_index > var_in%bc(bc_index)%num_fields) &
    call mpp_error(FATAL, trim(error_header)//" field_index exceeds num_fields for" //&
                   trim(var_in%bc(bc_index)%name) )

  ! Do error checking on the i-dimension and determine the array offsets.
  if (present(idim)) then
    if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
      write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
      call mpp_error(FATAL, trim(error_msg))
    endif
    if (size(array_out,1) /= (1+idim(4)-idim(1))) then
      write (error_msg, *) trim(error_header), ' The declared i-dimension size of ', &
            (1+idim(4)-idim(1)), ' does not match the actual size of ', size(array_out,1)
      call mpp_error(FATAL, trim(error_msg))
    endif
    if ((var_in%iec-var_in%isc) /= (idim(3)-idim(2))) &
      call mpp_error(FATAL, trim(error_header)//" There is an i-direction computional domain size mismatch.")
    if ((idim(2)-idim(1) < halo) .or. (idim(4)-idim(3) < halo)) &
      call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the output array.")
    if (size(array_out,1) < 2*halo + 1 + var_in%iec - var_in%isc) then
      write (error_msg, *) trim(error_header), ' The target array with i-dimension size ', &
            (1+idim(4)-idim(1)), ' is too small to match the data of size ', &
            (2*halo + 1 + var_in%iec - var_in%isc)
      call mpp_error(FATAL, trim(error_msg))
    endif

    i_off = (1-idim(1)) + (idim(2)-var_in%isc)
  else
    if (size(array_out,1) < 2*halo + 1 + var_in%iec - var_in%isc) then
      write (error_msg, *) trim(error_header), ' The target array with i-dimension size ', &
            size(array_out,1), ' does not match the data of size ', &
            (2*halo + 1 + var_in%iec - var_in%isc)
      call mpp_error(FATAL, trim(error_msg))
    endif
    i_off = 1 - (var_in%isc-halo)
  endif

  ! Do error checking on the j-dimension and determine the array offsets.
  if (present(jdim)) then
    if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
      write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list ', jdim
      call mpp_error(FATAL, trim(error_msg))
    endif
    if (size(array_out,2) /= (1+jdim(4)-jdim(1))) then
      write (error_msg, *) trim(error_header), ' The declared j-dimension size of ', &
            (1+jdim(4)-jdim(1)), ' does not match the actual size of ', size(array_out,2)
      call mpp_error(FATAL, trim(error_msg))
    endif
    if ((var_in%jec-var_in%jsc) /= (jdim(3)-jdim(2))) &
      call mpp_error(FATAL, trim(error_header)//" There is an j-direction computional domain size mismatch.")
    if ((jdim(2)-jdim(1) < halo) .or. (jdim(4)-jdim(3) < halo)) &
      call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the output array.")
    if (size(array_out,2) < 2*halo + 1 + var_in%jec - var_in%jsc) then
      write (error_msg, *) trim(error_header), ' The target array with j-dimension size ', &
            (1+jdim(4)-jdim(1)), ' is too small to match the data of size ', &
            (2*halo + 1 + var_in%jec - var_in%jsc)
      call mpp_error(FATAL, trim(error_msg))
    endif

    j_off = (1-jdim(1)) + (jdim(2)-var_in%jsc)
  else
    if (size(array_out,2) < 2*halo + 1 + var_in%jec - var_in%jsc) then
      write (error_msg, *) trim(error_header), ' The target array with j-dimension size ', &
            size(array_out,2), ' does not match the data of size ', &
            (2*halo + 1 + var_in%jec - var_in%jsc)
      call mpp_error(FATAL, trim(error_msg))
    endif
    j_off = 1 - (var_in%jsc-halo)
  endif

  if ((k_in > var_in%ke) .or. (k_in < var_in%ks)) then
    write (error_msg, *) trim(error_header), ' The extracted k-index of ', k_in, &
           ' is outside of the valid range of ', var_in%ks, ' to ', var_in%ke
    call mpp_error(FATAL, trim(error_msg))
  endif

  do j=var_in%jsc-halo,var_in%jec+halo ; do i=var_in%isc-halo,var_in%iec+halo
    array_out(i+i_off,j+j_off) = scale * var_in%bc(bc_index)%field(field_index)%values(i,j,k_in)
  enddo ; enddo

end subroutine CT_extract_data_3d_2d

!> This subroutine extracts a single 3-d field from a coupler_3d_bc_type into
!! a three-dimensional array.
subroutine CT_extract_data_3d(var_in, bc_index, field_index, array_out, &
                              scale_factor, halo_size, idim, jdim)
  type(coupler_3d_bc_type),   intent(in)    :: var_in    !< BC_type structure with the data to extract
  integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                         !! that is being copied
  integer,                    intent(in)    :: field_index !< The index of the field in the
                                                         !! boundary condition that is being copied
  real, dimension(1:,1:,1:),  intent(out)   :: array_out !< The recipient array for the field; its size
                                                         !! must match the size of the data being copied
                                                         !! unless idim and jdim are supplied.
  real,             optional, intent(in)    :: scale_factor !< A scaling factor for the data that is being added
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer, dimension(4), optional, intent(in) :: idim    !< The data and computational domain extents of
                                                         !! the first dimension of the output array
                                                         !! in a non-decreasing list
  integer, dimension(4), optional, intent(in) :: jdim    !< The data and computational domain extents of
                                                         !! the second dimension of the output array
                                                         !! in a non-decreasing list
  character(len=256), parameter :: error_header = &
     '==>Error from coupler_types_mod (CT_extract_data_3d):'
  character(len=400)      :: error_msg

  real :: scale
  integer :: i, j, k, halo, i_off, j_off, k_off

  if (bc_index <= 0) then
    array_out(:,:,:) = 0.0
    return
  endif

  halo = 0 ; if (present(halo_size)) halo = halo_size
  scale = 1.0 ; if (present(scale_factor)) scale = scale_factor

  if ((var_in%isc-var_in%isd < halo) .or. (var_in%ied-var_in%iec < halo)) &
    call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the input structure.")
  if ((var_in%jsc-var_in%jsd < halo) .or. (var_in%jed-var_in%jec < halo)) &
    call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the input structure.")

  if (bc_index > var_in%num_bcs) &
    call mpp_error(FATAL, trim(error_header)//" bc_index exceeds var_in%num_bcs.")
  if (field_index > var_in%bc(bc_index)%num_fields) &
    call mpp_error(FATAL, trim(error_header)//" field_index exceeds num_fields for" //&
                   trim(var_in%bc(bc_index)%name) )

  ! Do error checking on the i-dimension and determine the array offsets.
  if (present(idim)) then
    if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
      write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
      call mpp_error(FATAL, trim(error_msg))
    endif
    if (size(array_out,1) /= (1+idim(4)-idim(1))) then
      write (error_msg, *) trim(error_header), ' The declared i-dimension size of ', &
            (1+idim(4)-idim(1)), ' does not match the actual size of ', size(array_out,1)
      call mpp_error(FATAL, trim(error_msg))
    endif
    if ((var_in%iec-var_in%isc) /= (idim(3)-idim(2))) &
      call mpp_error(FATAL, trim(error_header)//" There is an i-direction computional domain size mismatch.")
    if ((idim(2)-idim(1) < halo) .or. (idim(4)-idim(3) < halo)) &
      call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the output array.")
    if (size(array_out,1) < 2*halo + 1 + var_in%iec - var_in%isc) then
      write (error_msg, *) trim(error_header), ' The target array with i-dimension size ', &
            (1+idim(4)-idim(1)), ' is too small to match the data of size ', &
            (2*halo + 1 + var_in%iec - var_in%isc)
      call mpp_error(FATAL, trim(error_msg))
    endif

    i_off = (1-idim(1)) + (idim(2)-var_in%isc)
  else
    if (size(array_out,1) < 2*halo + 1 + var_in%iec - var_in%isc) then
      write (error_msg, *) trim(error_header), ' The target array with i-dimension size ', &
            size(array_out,1), ' does not match the data of size ', &
            (2*halo + 1 + var_in%iec - var_in%isc)
      call mpp_error(FATAL, trim(error_msg))
    endif
    i_off = 1 - (var_in%isc-halo)
  endif

  ! Do error checking on the j-dimension and determine the array offsets.
  if (present(jdim)) then
    if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
      write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list ', jdim
      call mpp_error(FATAL, trim(error_msg))
    endif
    if (size(array_out,2) /= (1+jdim(4)-jdim(1))) then
      write (error_msg, *) trim(error_header), ' The declared j-dimension size of ', &
            (1+jdim(4)-jdim(1)), ' does not match the actual size of ', size(array_out,2)
      call mpp_error(FATAL, trim(error_msg))
    endif
    if ((var_in%jec-var_in%jsc) /= (jdim(3)-jdim(2))) &
      call mpp_error(FATAL, trim(error_header)//" There is an j-direction computional domain size mismatch.")
    if ((jdim(2)-jdim(1) < halo) .or. (jdim(4)-jdim(3) < halo)) &
      call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the output array.")
    if (size(array_out,2) < 2*halo + 1 + var_in%jec - var_in%jsc) then
      write (error_msg, *) trim(error_header), ' The target array with j-dimension size ', &
            (1+jdim(4)-jdim(1)), ' is too small to match the data of size ', &
            (2*halo + 1 + var_in%jec - var_in%jsc)
      call mpp_error(FATAL, trim(error_msg))
    endif

    j_off = (1-jdim(1)) + (jdim(2)-var_in%jsc)
  else
    if (size(array_out,2) < 2*halo + 1 + var_in%jec - var_in%jsc) then
      write (error_msg, *) trim(error_header), ' The target array with j-dimension size ', &
            size(array_out,2), ' does not match the data of size ', &
            (2*halo + 1 + var_in%jec - var_in%jsc)
      call mpp_error(FATAL, trim(error_msg))
    endif
    j_off = 1 - (var_in%jsc-halo)
  endif

  if (size(array_out,3) /= 1 + var_in%ke - var_in%ks) then
    write (error_msg, *) trim(error_header), ' The target array with k-dimension size ', &
          size(array_out,3), ' does not match the data of size ', &
          (1 + var_in%ke - var_in%ks)
    call mpp_error(FATAL, trim(error_msg))
  endif
  k_off = 1 - var_in%ks

  do k=var_in%ks,var_in%ke ; do j=var_in%jsc-halo,var_in%jec+halo ; do i=var_in%isc-halo,var_in%iec+halo
    array_out(i+i_off,j+j_off,k+k_off) = scale * var_in%bc(bc_index)%field(field_index)%values(i,j,k)
  enddo ; enddo ; enddo

end subroutine CT_extract_data_3d


!> This subroutine sets a single 2-d field in a coupler_3d_bc_type from
!! a two-dimensional array.
subroutine CT_set_data_2d(array_in, bc_index, field_index, var, &
                          scale_factor, halo_size, idim, jdim)
  real, dimension(1:,1:),     intent(in)   :: array_in   !< The source array for the field; its size
                                                         !! must match the size of the data being copied
                                                         !! unless idim and jdim are supplied.
  integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                         !! that is being copied
  integer,                    intent(in)    :: field_index !< The index of the field in the
                                                         !! boundary condition that is being copied
  type(coupler_2d_bc_type),   intent(inout) :: var       !< BC_type structure with the data to set
  real,             optional, intent(in)    :: scale_factor !< A scaling factor for the data that is being added
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer, dimension(4), optional, intent(in) :: idim    !< The data and computational domain extents of
                                                         !! the first dimension of the output array
                                                         !! in a non-decreasing list
  integer, dimension(4), optional, intent(in) :: jdim    !< The data and computational domain extents of
                                                         !! the second dimension of the output array
                                                         !! in a non-decreasing list
  character(len=256), parameter :: error_header = &
     '==>Error from coupler_types_mod (CT_set_data_2d):'
  character(len=400)      :: error_msg

  real :: scale
  integer :: i, j, halo, i_off, j_off

  if (bc_index <= 0) return

  halo = 0 ; if (present(halo_size)) halo = halo_size
  scale = 1.0 ; if (present(scale_factor)) scale = scale_factor

  if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo)) &
    call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the input structure.")
  if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo)) &
    call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the input structure.")

  if (bc_index > var%num_bcs) &
    call mpp_error(FATAL, trim(error_header)//" bc_index exceeds var%num_bcs.")
  if (field_index > var%bc(bc_index)%num_fields) &
    call mpp_error(FATAL, trim(error_header)//" field_index exceeds num_fields for" //&
                   trim(var%bc(bc_index)%name) )

  ! Do error checking on the i-dimension and determine the array offsets.
  if (present(idim)) then
    if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
      write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
      call mpp_error(FATAL, trim(error_msg))
    endif
    if (size(array_in,1) /= (1+idim(4)-idim(1))) then
      write (error_msg, *) trim(error_header), ' The declared i-dimension size of ', &
            (1+idim(4)-idim(1)), ' does not match the actual size of ', size(array_in,1)
      call mpp_error(FATAL, trim(error_msg))
    endif
    if ((var%iec-var%isc) /= (idim(3)-idim(2))) &
      call mpp_error(FATAL, trim(error_header)//" There is an i-direction computional domain size mismatch.")
    if ((idim(2)-idim(1) < halo) .or. (idim(4)-idim(3) < halo)) &
      call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the output array.")
    if (size(array_in,1) < 2*halo + 1 + var%iec - var%isc) then
      write (error_msg, *) trim(error_header), ' The target array with i-dimension size ', &
            (1+idim(4)-idim(1)), ' is too small to match the data of size ', &
            (2*halo + 1 + var%iec - var%isc)
      call mpp_error(FATAL, trim(error_msg))
    endif

    i_off = (1-idim(1)) + (idim(2)-var%isc)
  else
    if (size(array_in,1) < 2*halo + 1 + var%iec - var%isc) then
      write (error_msg, *) trim(error_header), ' The target array with i-dimension size ', &
            size(array_in,1), ' does not match the data of size ', &
            (2*halo + 1 + var%iec - var%isc)
      call mpp_error(FATAL, trim(error_msg))
    endif
    i_off = 1 - (var%isc-halo)
  endif

  ! Do error checking on the j-dimension and determine the array offsets.
  if (present(jdim)) then
    if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
      write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list ', jdim
      call mpp_error(FATAL, trim(error_msg))
    endif
    if (size(array_in,2) /= (1+jdim(4)-jdim(1))) then
      write (error_msg, *) trim(error_header), ' The declared j-dimension size of ', &
            (1+jdim(4)-jdim(1)), ' does not match the actual size of ', size(array_in,2)
      call mpp_error(FATAL, trim(error_msg))
    endif
    if ((var%jec-var%jsc) /= (jdim(3)-jdim(2))) &
      call mpp_error(FATAL, trim(error_header)//" There is an j-direction computional domain size mismatch.")
    if ((jdim(2)-jdim(1) < halo) .or. (jdim(4)-jdim(3) < halo)) &
      call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the output array.")
    if (size(array_in,2) < 2*halo + 1 + var%jec - var%jsc) then
      write (error_msg, *) trim(error_header), ' The target array with j-dimension size ', &
            (1+jdim(4)-jdim(1)), ' is too small to match the data of size ', &
            (2*halo + 1 + var%jec - var%jsc)
      call mpp_error(FATAL, trim(error_msg))
    endif

    j_off = (1-jdim(1)) + (jdim(2)-var%jsc)
  else
    if (size(array_in,2) < 2*halo + 1 + var%jec - var%jsc) then
      write (error_msg, *) trim(error_header), ' The target array with j-dimension size ', &
            size(array_in,2), ' does not match the data of size ', &
            (2*halo + 1 + var%jec - var%jsc)
      call mpp_error(FATAL, trim(error_msg))
    endif
    j_off = 1 - (var%jsc-halo)
  endif

  do j=var%jsc-halo,var%jec+halo ; do i=var%isc-halo,var%iec+halo
    var%bc(bc_index)%field(field_index)%values(i,j) = scale * array_in(i+i_off,j+j_off)
  enddo ; enddo

end subroutine CT_set_data_2d

!> This subroutine sets a one k-level of a single 3-d field in a
!! coupler_3d_bc_type from a two-dimensional array.
subroutine CT_set_data_2d_3d(array_in, bc_index, field_index, k_out, var, &
                             scale_factor, halo_size, idim, jdim)
  real, dimension(1:,1:),     intent(in)    :: array_in  !< The source array for the field; its size
                                                         !! must match the size of the data being copied
                                                         !! unless idim and jdim are supplied.
  integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                         !! that is being copied
  integer,                    intent(in)    :: field_index !< The index of the field in the
                                                         !! boundary condition that is being copied
  integer,                    intent(in)    :: k_out     !< The k-index to set
  type(coupler_3d_bc_type),   intent(inout) :: var       !< BC_type structure with the data to be set
  real,             optional, intent(in)    :: scale_factor !< A scaling factor for the data that is being added
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer, dimension(4), optional, intent(in) :: idim    !< The data and computational domain extents of
                                                         !! the first dimension of the output array
                                                         !! in a non-decreasing list
  integer, dimension(4), optional, intent(in) :: jdim    !< The data and computational domain extents of
                                                         !! the second dimension of the output array
                                                         !! in a non-decreasing list
  character(len=256), parameter :: error_header = &
     '==>Error from coupler_types_mod (CT_set_data_3d_2d):'
  character(len=400)      :: error_msg

  real :: scale
  integer :: i, j, halo, i_off, j_off

  if (bc_index <= 0) return

  halo = 0 ; if (present(halo_size)) halo = halo_size
  scale = 1.0 ; if (present(scale_factor)) scale = scale_factor

  if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo)) &
    call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the input structure.")
  if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo)) &
    call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the input structure.")

  if (bc_index > var%num_bcs) &
    call mpp_error(FATAL, trim(error_header)//" bc_index exceeds var%num_bcs.")
  if (field_index > var%bc(bc_index)%num_fields) &
    call mpp_error(FATAL, trim(error_header)//" field_index exceeds num_fields for" //&
                   trim(var%bc(bc_index)%name) )

  ! Do error checking on the i-dimension and determine the array offsets.
  if (present(idim)) then
    if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
      write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
      call mpp_error(FATAL, trim(error_msg))
    endif
    if (size(array_in,1) /= (1+idim(4)-idim(1))) then
      write (error_msg, *) trim(error_header), ' The declared i-dimension size of ', &
            (1+idim(4)-idim(1)), ' does not match the actual size of ', size(array_in,1)
      call mpp_error(FATAL, trim(error_msg))
    endif
    if ((var%iec-var%isc) /= (idim(3)-idim(2))) &
      call mpp_error(FATAL, trim(error_header)//" There is an i-direction computional domain size mismatch.")
    if ((idim(2)-idim(1) < halo) .or. (idim(4)-idim(3) < halo)) &
      call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the output array.")
    if (size(array_in,1) < 2*halo + 1 + var%iec - var%isc) then
      write (error_msg, *) trim(error_header), ' The target array with i-dimension size ', &
            (1+idim(4)-idim(1)), ' is too small to match the data of size ', &
            (2*halo + 1 + var%iec - var%isc)
      call mpp_error(FATAL, trim(error_msg))
    endif

    i_off = (1-idim(1)) + (idim(2)-var%isc)
  else
    if (size(array_in,1) < 2*halo + 1 + var%iec - var%isc) then
      write (error_msg, *) trim(error_header), ' The target array with i-dimension size ', &
            size(array_in,1), ' does not match the data of size ', &
            (2*halo + 1 + var%iec - var%isc)
      call mpp_error(FATAL, trim(error_msg))
    endif
    i_off = 1 - (var%isc-halo)
  endif

  ! Do error checking on the j-dimension and determine the array offsets.
  if (present(jdim)) then
    if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
      write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list ', jdim
      call mpp_error(FATAL, trim(error_msg))
    endif
    if (size(array_in,2) /= (1+jdim(4)-jdim(1))) then
      write (error_msg, *) trim(error_header), ' The declared j-dimension size of ', &
            (1+jdim(4)-jdim(1)), ' does not match the actual size of ', size(array_in,2)
      call mpp_error(FATAL, trim(error_msg))
    endif
    if ((var%jec-var%jsc) /= (jdim(3)-jdim(2))) &
      call mpp_error(FATAL, trim(error_header)//" There is an j-direction computional domain size mismatch.")
    if ((jdim(2)-jdim(1) < halo) .or. (jdim(4)-jdim(3) < halo)) &
      call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the output array.")
    if (size(array_in,2) < 2*halo + 1 + var%jec - var%jsc) then
      write (error_msg, *) trim(error_header), ' The target array with j-dimension size ', &
            (1+jdim(4)-jdim(1)), ' is too small to match the data of size ', &
            (2*halo + 1 + var%jec - var%jsc)
      call mpp_error(FATAL, trim(error_msg))
    endif

    j_off = (1-jdim(1)) + (jdim(2)-var%jsc)
  else
    if (size(array_in,2) < 2*halo + 1 + var%jec - var%jsc) then
      write (error_msg, *) trim(error_header), ' The target array with j-dimension size ', &
            size(array_in,2), ' does not match the data of size ', &
            (2*halo + 1 + var%jec - var%jsc)
      call mpp_error(FATAL, trim(error_msg))
    endif
    j_off = 1 - (var%jsc-halo)
  endif

  if ((k_out > var%ke) .or. (k_out < var%ks)) then
    write (error_msg, *) trim(error_header), ' The seted k-index of ', k_out, &
           ' is outside of the valid range of ', var%ks, ' to ', var%ke
    call mpp_error(FATAL, trim(error_msg))
  endif

  do j=var%jsc-halo,var%jec+halo ; do i=var%isc-halo,var%iec+halo
    var%bc(bc_index)%field(field_index)%values(i,j,k_out) = scale * array_in(i+i_off,j+j_off)
  enddo ; enddo

end subroutine CT_set_data_2d_3d

!> This subroutine sets a single 3-d field in a coupler_3d_bc_type from
!! a three-dimensional array.
subroutine CT_set_data_3d(array_in, bc_index, field_index, var, &
                          scale_factor, halo_size, idim, jdim)
  real, dimension(1:,1:,1:),  intent(in)    :: array_in  !< The source array for the field; its size
                                                         !! must match the size of the data being copied
                                                         !! unless idim and jdim are supplied.
  integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                         !! that is being copied
  integer,                    intent(in)    :: field_index !< The index of the field in the
                                                         !! boundary condition that is being copied
  type(coupler_3d_bc_type),   intent(inout) :: var       !< BC_type structure with the data to be set
  real,             optional, intent(in)    :: scale_factor !< A scaling factor for the data that is being added
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer, dimension(4), optional, intent(in) :: idim    !< The data and computational domain extents of
                                                         !! the first dimension of the output array
                                                         !! in a non-decreasing list
  integer, dimension(4), optional, intent(in) :: jdim    !< The data and computational domain extents of
                                                         !! the second dimension of the output array
                                                         !! in a non-decreasing list
  character(len=256), parameter :: error_header = &
     '==>Error from coupler_types_mod (CT_set_data_3d):'
  character(len=400)      :: error_msg

  real :: scale
  integer :: i, j, k, halo, i_off, j_off, k_off

  if (bc_index <= 0) return

  halo = 0 ; if (present(halo_size)) halo = halo_size
  scale = 1.0 ; if (present(scale_factor)) scale = scale_factor

  if ((var%isc-var%isd < halo) .or. (var%ied-var%iec < halo)) &
    call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the input structure.")
  if ((var%jsc-var%jsd < halo) .or. (var%jed-var%jec < halo)) &
    call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the input structure.")

  if (bc_index > var%num_bcs) &
    call mpp_error(FATAL, trim(error_header)//" bc_index exceeds var%num_bcs.")
  if (field_index > var%bc(bc_index)%num_fields) &
    call mpp_error(FATAL, trim(error_header)//" field_index exceeds num_fields for" //&
                   trim(var%bc(bc_index)%name) )

  ! Do error checking on the i-dimension and determine the array offsets.
  if (present(idim)) then
    if ((idim(1) > idim(2)) .or. (idim(3) > idim(4))) then
      write (error_msg, *) trim(error_header), ' Disordered i-dimension index bound list ', idim
      call mpp_error(FATAL, trim(error_msg))
    endif
    if (size(array_in,1) /= (1+idim(4)-idim(1))) then
      write (error_msg, *) trim(error_header), ' The declared i-dimension size of ', &
            (1+idim(4)-idim(1)), ' does not match the actual size of ', size(array_in,1)
      call mpp_error(FATAL, trim(error_msg))
    endif
    if ((var%iec-var%isc) /= (idim(3)-idim(2))) &
      call mpp_error(FATAL, trim(error_header)//" There is an i-direction computional domain size mismatch.")
    if ((idim(2)-idim(1) < halo) .or. (idim(4)-idim(3) < halo)) &
      call mpp_error(FATAL, trim(error_header)//" Excessive i-direction halo size for the output array.")
    if (size(array_in,1) < 2*halo + 1 + var%iec - var%isc) then
      write (error_msg, *) trim(error_header), ' The target array with i-dimension size ', &
            (1+idim(4)-idim(1)), ' is too small to match the data of size ', &
            (2*halo + 1 + var%iec - var%isc)
      call mpp_error(FATAL, trim(error_msg))
    endif

    i_off = (1-idim(1)) + (idim(2)-var%isc)
  else
    if (size(array_in,1) < 2*halo + 1 + var%iec - var%isc) then
      write (error_msg, *) trim(error_header), ' The target array with i-dimension size ', &
            size(array_in,1), ' does not match the data of size ', &
            (2*halo + 1 + var%iec - var%isc)
      call mpp_error(FATAL, trim(error_msg))
    endif
    i_off = 1 - (var%isc-halo)
  endif

  ! Do error checking on the j-dimension and determine the array offsets.
  if (present(jdim)) then
    if ((jdim(1) > jdim(2)) .or. (jdim(3) > jdim(4))) then
      write (error_msg, *) trim(error_header), ' Disordered j-dimension index bound list ', jdim
      call mpp_error(FATAL, trim(error_msg))
    endif
    if (size(array_in,2) /= (1+jdim(4)-jdim(1))) then
      write (error_msg, *) trim(error_header), ' The declared j-dimension size of ', &
            (1+jdim(4)-jdim(1)), ' does not match the actual size of ', size(array_in,2)
      call mpp_error(FATAL, trim(error_msg))
    endif
    if ((var%jec-var%jsc) /= (jdim(3)-jdim(2))) &
      call mpp_error(FATAL, trim(error_header)//" There is an j-direction computional domain size mismatch.")
    if ((jdim(2)-jdim(1) < halo) .or. (jdim(4)-jdim(3) < halo)) &
      call mpp_error(FATAL, trim(error_header)//" Excessive j-direction halo size for the output array.")
    if (size(array_in,2) < 2*halo + 1 + var%jec - var%jsc) then
      write (error_msg, *) trim(error_header), ' The target array with j-dimension size ', &
            (1+jdim(4)-jdim(1)), ' is too small to match the data of size ', &
            (2*halo + 1 + var%jec - var%jsc)
      call mpp_error(FATAL, trim(error_msg))
    endif

    j_off = (1-jdim(1)) + (jdim(2)-var%jsc)
  else
    if (size(array_in,2) < 2*halo + 1 + var%jec - var%jsc) then
      write (error_msg, *) trim(error_header), ' The target array with j-dimension size ', &
            size(array_in,2), ' does not match the data of size ', &
            (2*halo + 1 + var%jec - var%jsc)
      call mpp_error(FATAL, trim(error_msg))
    endif
    j_off = 1 - (var%jsc-halo)
  endif

  if (size(array_in,3) /= 1 + var%ke - var%ks) then
    write (error_msg, *) trim(error_header), ' The target array with k-dimension size ', &
          size(array_in,3), ' does not match the data of size ', &
          (1 + var%ke - var%ks)
    call mpp_error(FATAL, trim(error_msg))
  endif
  k_off = 1 - var%ks

  do k=var%ks,var%ke ; do j=var%jsc-halo,var%jec+halo ; do i=var%isc-halo,var%iec+halo
    var%bc(bc_index)%field(field_index)%values(i,j,k) = scale * array_in(i+i_off,j+j_off,k+k_off)
  enddo ; enddo ; enddo

end subroutine CT_set_data_3d


!> This routine registers the diagnostics of a coupler_2d_bc_type.
subroutine CT_set_diags_2d(var, diag_name, axes, time)
  type(coupler_2d_bc_type), intent(inout) :: var  !< BC_type structure for which to register diagnostics
  character(len=*),         intent(in)    :: diag_name !< name for diagnostic file--if blank, then don't register the fields
  integer, dimension(:),    intent(in)    :: axes !< array of axes identifiers for diagnostic variable registration
  type(time_type),          intent(in)    :: time !< model time variable for registering diagnostic field

  integer :: m, n

  if (diag_name == ' ') return

  if (size(axes) < 2) then
    call mpp_error(FATAL, '==>Error from coupler_types_mod' //&
             '(coupler_types_set_diags_3d): axes has less than 2 elements')
  endif

  do n = 1, var%num_bcs
    do m = 1, var%bc(n)%num_fields
      var%bc(n)%field(m)%id_diag = register_diag_field(diag_name,                &
           var%bc(n)%field(m)%name, axes(1:2), Time,                             &
           var%bc(n)%field(m)%long_name, var%bc(n)%field(m)%units )
    enddo
  enddo

end subroutine CT_set_diags_2d

!> This routine registers the diagnostics of a coupler_3d_bc_type.
subroutine CT_set_diags_3d(var, diag_name, axes, time)
  type(coupler_3d_bc_type), intent(inout) :: var  !< BC_type structure for which to register diagnostics
  character(len=*),         intent(in)    :: diag_name !< name for diagnostic file--if blank, then don't register the fields
  integer, dimension(:),    intent(in)    :: axes !< array of axes identifiers for diagnostic variable registration
  type(time_type),          intent(in)    :: time !< model time variable for registering diagnostic field

  integer :: m, n

  if (diag_name == ' ') return

  if (size(axes) < 3) then
    call mpp_error(FATAL, '==>Error from coupler_types_mod' //&
             '(coupler_types_set_diags_3d): axes has less than 3 elements')
  endif

  do n = 1, var%num_bcs
    do m = 1, var%bc(n)%num_fields
      var%bc(n)%field(m)%id_diag = register_diag_field(diag_name,                &
           var%bc(n)%field(m)%name, axes(1:3), Time,                             &
           var%bc(n)%field(m)%long_name, var%bc(n)%field(m)%units )
    enddo
  enddo

end subroutine CT_set_diags_3d


!> This subroutine writes out all diagnostics of elements of a coupler_2d_bc_type
subroutine CT_send_data_2d(var, Time)
  type(coupler_2d_bc_type), intent(in) :: var  !< BC_type structure with the diagnostics to write
  type(time_type),          intent(in) :: time !< The current model time

  integer :: m, n
  logical :: used

  do n = 1, var%num_bcs ; do m = 1, var%bc(n)%num_fields
    used = send_data(var%bc(n)%field(m)%id_diag, var%bc(n)%field(m)%values, Time)
  enddo ; enddo

end subroutine CT_send_data_2d

!> This subroutine writes out all diagnostics of elements of a coupler_2d_bc_type
subroutine CT_send_data_3d(var, Time)
  type(coupler_3d_bc_type), intent(in) :: var  !< BC_type structure with the diagnostics to write
  type(time_type),          intent(in) :: time !< The current model time

  integer :: m, n
  logical :: used

  do n = 1, var%num_bcs ; do m = 1, var%bc(n)%num_fields
    used = send_data(var%bc(n)%field(m)%id_diag, var%bc(n)%field(m)%values, Time)
  enddo ; enddo

end subroutine CT_send_data_3d


!> This subroutine registers the fields in a coupler_2d_bc_type to be saved
!! in restart files specified in the field table.
subroutine CT_register_restarts_2d(var, bc_rest_files, num_rest_files, mpp_domain, ocean_restart)
  type(coupler_2d_bc_type), intent(inout) :: var  !< BC_type structure to be registered for restarts
  type(restart_file_type),  dimension(:), pointer :: bc_rest_files !< Structures describing the restart files
  integer,                  intent(out) :: num_rest_files !< The number of restart files to use
  type(domain2D),           intent(in)  :: mpp_domain     !< The FMS domain to use for this registration call
  logical,        optional, intent(in)  :: ocean_restart  !< If true, use the ocean restart file name.

  character(len=80), dimension(max(1,var%num_bcs)) :: rest_file_names
  character(len=80) :: file_nm
  logical :: ocn_rest
  integer :: f, n, m

  ocn_rest = .true. ; if (present(ocean_restart)) ocn_rest = ocean_restart

  ! Determine the number and names of the restart files
  num_rest_files = 0
  do n = 1, var%num_bcs
    if (var%bc(n)%num_fields <= 0) cycle
    file_nm = trim(var%bc(n)%ice_restart_file)
    if (ocn_rest) file_nm = trim(var%bc(n)%ocean_restart_file)
    do f = 1, num_rest_files
      if (trim(file_nm) == trim(rest_file_names(f))) exit
    enddo
    if (f>num_rest_files) then
      num_rest_files = num_rest_files + 1
      rest_file_names(f) = trim(file_nm)
    endif
  enddo

  if (num_rest_files == 0) return

  ! Register the fields with the restart files
  allocate(bc_rest_files(num_rest_files))
  do n = 1, var%num_bcs
    if (var%bc(n)%num_fields <= 0) cycle

    file_nm = trim(var%bc(n)%ice_restart_file)
    if (ocn_rest) file_nm = trim(var%bc(n)%ocean_restart_file)
    do f = 1, num_rest_files
      if (trim(file_nm) == trim(rest_file_names(f))) exit
    enddo

    var%bc(n)%rest_type => bc_rest_files(f)
    do m = 1, var%bc(n)%num_fields
      var%bc(n)%field(m)%id_rest = register_restart_field(bc_rest_files(f), &
              rest_file_names(f), var%bc(n)%field(m)%name, var%bc(n)%field(m)%values, &
              mpp_domain, mandatory=.not.var%bc(n)%field(m)%may_init )
    enddo
  enddo

end subroutine CT_register_restarts_2d

!> This subroutine registers the fields in a coupler_2d_bc_type to be saved
!! in the specified restart file.
subroutine CT_register_restarts_to_file_2d(var, file_name, rest_file, mpp_domain, &
                                           varname_prefix)
  type(coupler_2d_bc_type), intent(inout) :: var  !< BC_type structure to be registered for restarts
  character(len=*),         intent(in)    :: file_name !< The name of the restart file
  type(restart_file_type),  pointer       :: rest_file !< A (possibly associated) structure describing the restart file
  type(domain2D),           intent(in)    :: mpp_domain !< The FMS domain to use for this registration call
  character(len=*), optional, intent(in)  :: varname_prefix !< A prefix for the variable name
                                                  !! in the restart file, intended to allow
                                                  !! multiple BC_type variables to use the
                                                  !! same restart files.

  character(len=128) :: var_name
  integer :: n, m

  ! Register the fields with the restart file
  if (.not.associated(rest_file)) allocate(rest_file)
  do n = 1, var%num_bcs
    if (var%bc(n)%num_fields <= 0) cycle

    var%bc(n)%rest_type => rest_file
    do m = 1, var%bc(n)%num_fields
      var_name = trim(var%bc(n)%field(m)%name)
      if (present(varname_prefix)) var_name = trim(varname_prefix)//trim(var_name)
      var%bc(n)%field(m)%id_rest = register_restart_field(rest_file, &
              file_name, var_name, var%bc(n)%field(m)%values, &
              mpp_domain, mandatory=.not.var%bc(n)%field(m)%may_init )
    enddo
  enddo

end subroutine CT_register_restarts_to_file_2d

!> This subroutine registers the fields in a coupler_3d_bc_type to be saved
!! in restart files specified in the field table.
subroutine CT_register_restarts_3d(var, bc_rest_files, num_rest_files, mpp_domain, ocean_restart)
  type(coupler_3d_bc_type), intent(inout) :: var  !< BC_type structure to be registered for restarts
  type(restart_file_type),  dimension(:), pointer :: bc_rest_files !< Structures describing the restart files
  integer,                  intent(out)   :: num_rest_files !< The number of restart files to use
  type(domain2D),           intent(in)    :: mpp_domain     !< The FMS domain to use for this registration call
  logical,        optional, intent(in)    :: ocean_restart  !< If true, use the ocean restart file name.

  character(len=80), dimension(max(1,var%num_bcs)) :: rest_file_names
  character(len=80) :: file_nm
  logical :: ocn_rest
  integer :: f, n, m, id_restart

  ocn_rest = .true. ; if (present(ocean_restart)) ocn_rest = ocean_restart

  ! Determine the number and names of the restart files
  num_rest_files = 0
  do n = 1, var%num_bcs
    if (var%bc(n)%num_fields <= 0) cycle
    file_nm = trim(var%bc(n)%ice_restart_file)
    if (ocn_rest) file_nm = trim(var%bc(n)%ocean_restart_file)
    do f = 1, num_rest_files
      if (trim(file_nm) == trim(rest_file_names(f))) exit
    enddo
    if (f>num_rest_files) then
      num_rest_files = num_rest_files + 1
      rest_file_names(f) = trim(file_nm)
    endif
  enddo

  if (num_rest_files == 0) return

  ! Register the fields with the restart files
  allocate(bc_rest_files(num_rest_files))
  do n = 1, var%num_bcs
    if (var%bc(n)%num_fields <= 0) cycle
    file_nm = trim(var%bc(n)%ice_restart_file)
    if (ocn_rest) file_nm = trim(var%bc(n)%ocean_restart_file)
    do f = 1, num_rest_files
      if (trim(file_nm) == trim(rest_file_names(f))) exit
    enddo

    var%bc(n)%rest_type => bc_rest_files(f)
    do m = 1, var%bc(n)%num_fields
      var%bc(n)%field(m)%id_rest = register_restart_field(bc_rest_files(f), &
              rest_file_names(f), var%bc(n)%field(m)%name, var%bc(n)%field(m)%values, &
              mpp_domain, mandatory=.not.var%bc(n)%field(m)%may_init )
    enddo
  enddo

end subroutine CT_register_restarts_3d

!> This subroutine registers the fields in a coupler_3d_bc_type to be saved
!! in the specified restart file.
subroutine CT_register_restarts_to_file_3d(var, file_name, rest_file, mpp_domain, &
                                           varname_prefix)
  type(coupler_3d_bc_type), intent(inout) :: var  !< BC_type structure to be registered for restarts
  character(len=*),         intent(in)  :: file_name !< The name of the restart file
  type(restart_file_type),  pointer     :: rest_file !< A (possibly associated) structure describing the restart file
  type(domain2D),           intent(in)  :: mpp_domain     !< The FMS domain to use for this registration call

  character(len=*), optional, intent(in)  :: varname_prefix !< A prefix for the variable name
                                                  !! in the restart file, intended to allow
                                                  !! multiple BC_type variables to use the
                                                  !! same restart files.

  character(len=128) :: var_name
  integer :: n, m

  ! Register the fields with the restart file
  if (.not.associated(rest_file)) allocate(rest_file)
  do n = 1, var%num_bcs
    if (var%bc(n)%num_fields <= 0) cycle

    var%bc(n)%rest_type => rest_file
    do m = 1, var%bc(n)%num_fields
      var_name = trim(var%bc(n)%field(m)%name)
      if (present(varname_prefix)) var_name = trim(varname_prefix)//trim(var_name)
      var%bc(n)%field(m)%id_rest = register_restart_field(rest_file, &
              file_name, var_name, var%bc(n)%field(m)%values, &
              mpp_domain, mandatory=.not.var%bc(n)%field(m)%may_init )
    enddo
  enddo

end subroutine CT_register_restarts_to_file_3d


!> This subroutine reads in the fields in a coupler_2d_bc_type that have
!! been saved in restart files.
subroutine CT_restore_state_2d(var, directory, all_or_nothing, &
                                         all_required, test_by_field)
  type(coupler_2d_bc_type), intent(inout) :: var  !< BC_type structure to restore from restart files
  character(len=*), optional, intent(in)  :: directory !< A directory where the restart files should
                                                  !! be found.  The default for FMS is 'INPUT'.
  logical,        optional, intent(in)    :: all_or_nothing !< If true and there are non-mandatory
                                                  !! restart fields, it is still an error if some
                                                  !! fields are read successfully but others are not.
  logical,        optional, intent(in)    :: all_required !< If true, all fields must be successfully
                                                  !! read from the restart file, even if they were
                                                  !! registered as not mandatory.
  logical,        optional, intent(in)    :: test_by_field !< If true, all or none of the variables
                                                  !! in a single field must be read successfully.

  integer :: n, m, num_fld
  character(len=80) :: unset_varname
  logical :: any_set, all_set, all_var_set, any_var_set, var_set

  any_set = .false. ; all_set = .true. ; num_fld = 0 ; unset_varname = ""

  do n = 1, var%num_bcs
    any_var_set = .false. ; all_var_set = .true.
    do m = 1, var%bc(n)%num_fields
      var_set = .false.
      if (var%bc(n)%field(m)%id_rest > 0) then
        var_set = query_initialized(var%bc(n)%rest_type, var%bc(n)%field(m)%id_rest)
        if (.not.var_set) then
          call restore_state(var%bc(n)%rest_type, var%bc(n)%field(m)%id_rest, &
                             directory=directory, nonfatal_missing_files=.true.)
          var_set = query_initialized(var%bc(n)%rest_type, var%bc(n)%field(m)%id_rest)
        endif
      endif

      if (.not.var_set) unset_varname = trim(var%bc(n)%field(m)%name)
      if (var_set) any_set = .true.
      if (all_set) all_set = var_set
      if (var_set) any_var_set = .true.
      if (all_var_set) all_var_set = var_set
    enddo

    num_fld = num_fld + var%bc(n)%num_fields
    if ((var%bc(n)%num_fields > 0) .and. present(test_by_field)) then
      if (test_by_field .and. (all_var_set .neqv. any_var_set)) call mpp_error(FATAL, &
             "CT_restore_state_2d: test_by_field is true, and "//&
             trim(unset_varname)//" was not read but some other fields in "//&
             trim(trim(var%bc(n)%name))//" were.")
    endif
  enddo

  if ((num_fld > 0) .and. present(all_or_nothing)) then
    if (all_or_nothing .and. (all_set .neqv. any_set)) call mpp_error(FATAL, &
           "CT_restore_state_2d: all_or_nothing is true, and "//&
           trim(unset_varname)//" was not read but some other fields were.")
  endif

  if (present(all_required)) then ; if (all_required .and. .not.all_set) then
    call mpp_error(FATAL, "CT_restore_state_2d: all_required is true, but "//&
           trim(unset_varname)//" was not read from its restart file.")
  endif ; endif

end subroutine CT_restore_state_2d

!> This subroutine reads in the fields in a coupler_3d_bc_type that have
!! been saved in restart files.
subroutine CT_restore_state_3d(var, directory, all_or_nothing, &
                                         all_required, test_by_field)
  type(coupler_3d_bc_type), intent(inout) :: var  !< BC_type structure to restore from restart files
  character(len=*), optional, intent(in)  :: directory !< A directory where the restart files should
                                                  !! be found.  The default for FMS is 'INPUT'.
  logical,        optional, intent(in)    :: all_or_nothing !< If true and there are non-mandatory
                                                  !! restart fields, it is still an error if some
                                                  !! fields are read successfully but others are not.
  logical,        optional, intent(in)    :: all_required !< If true, all fields must be successfully
                                                  !! read from the restart file, even if they were
                                                  !! registered as not mandatory.
  logical,        optional, intent(in)    :: test_by_field !< If true, all or none of the variables
                                                  !! in a single field must be read successfully.

  integer :: n, m, num_fld
  character(len=80) :: unset_varname
  logical :: any_set, all_set, all_var_set, any_var_set, var_set

  any_set = .false. ; all_set = .true. ; num_fld = 0 ; unset_varname = ""

  do n = 1, var%num_bcs
    any_var_set = .false. ; all_var_set = .true.
    do m = 1, var%bc(n)%num_fields
      var_set = .false.
      if (var%bc(n)%field(m)%id_rest > 0) then
        var_set = query_initialized(var%bc(n)%rest_type, var%bc(n)%field(m)%id_rest)
        if (.not.var_set) then
          call restore_state(var%bc(n)%rest_type, var%bc(n)%field(m)%id_rest, &
                             directory=directory, nonfatal_missing_files=.true.)
          var_set = query_initialized(var%bc(n)%rest_type, var%bc(n)%field(m)%id_rest)
        endif
      endif

      if (.not.var_set) unset_varname = trim(var%bc(n)%field(m)%name)

      if (var_set) any_set = .true.
      if (all_set) all_set = var_set
      if (var_set) any_var_set = .true.
      if (all_var_set) all_var_set = var_set
    enddo

    num_fld = num_fld + var%bc(n)%num_fields
    if ((var%bc(n)%num_fields > 0) .and. present(test_by_field)) then
      if (test_by_field .and. (all_var_set .neqv. any_var_set)) call mpp_error(FATAL, &
             "CT_restore_state_3d: test_by_field is true, and "//&
             trim(unset_varname)//" was not read but some other fields in "//&
             trim(trim(var%bc(n)%name))//" were.")
    endif
  enddo

  if ((num_fld > 0) .and. present(all_or_nothing)) then
    if (all_or_nothing .and. (all_set .neqv. any_set)) call mpp_error(FATAL, &
           "CT_restore_state_3d: all_or_nothing is true, and "//&
           trim(unset_varname)//" was not read but some other fields were.")
  endif

  if (present(all_required)) then ; if (all_required .and. .not.all_set) then
    call mpp_error(FATAL, "CT_restore_state_3d: all_required is true, but "//&
           trim(unset_varname)//" was not read from its restart file.")
  endif ; endif

end subroutine CT_restore_state_3d


!> This subroutine potentially overrides the values in a coupler_2d_bc_type
subroutine CT_data_override_2d(gridname, var, Time)
  character(len=3),         intent(in)    :: gridname !< 3-character long model grid ID
  type(coupler_2d_bc_type), intent(inout) :: var  !< BC_type structure to override
  type(time_type),          intent(in)    :: time !< The current model time

  integer :: m, n

  do n = 1, var%num_bcs ; do m = 1, var%bc(n)%num_fields
    call data_override(gridname, var%bc(n)%field(m)%name, var%bc(n)%field(m)%values, Time)
  enddo ; enddo

end subroutine CT_data_override_2d

!> This subroutine potentially overrides the values in a coupler_3d_bc_type
subroutine CT_data_override_3d(gridname, var, Time)
  character(len=3),         intent(in)    :: gridname !< 3-character long model grid ID
  type(coupler_3d_bc_type), intent(inout) :: var  !< BC_type structure to override
  type(time_type),          intent(in)    :: time !< The current model time

  integer :: m, n

  do n = 1, var%num_bcs ; do m = 1, var%bc(n)%num_fields
    call data_override(gridname, var%bc(n)%field(m)%name, var%bc(n)%field(m)%values, Time)
  enddo ; enddo

end subroutine CT_data_override_3d


!> This subroutine writes out checksums for the elements of a coupler_2d_bc_type
subroutine CT_write_chksums_2d(var, outunit, name_lead)
  type(coupler_2d_bc_type),   intent(in) :: var  !< BC_type structure for which to register diagnostics
  integer,                    intent(in) :: outunit !< The index of a open output file
  character(len=*), optional, intent(in) :: name_lead !< An optional prefix for the variable names

  character(len=120) :: var_name
  integer :: m, n

  do n = 1, var%num_bcs ; do m = 1, var%bc(n)%num_fields
    if (present(name_lead)) then
      var_name = trim(name_lead)//trim(var%bc(n)%field(m)%name)
    else
      var_name = trim(var%bc(n)%field(m)%name)
    endif
    write(outunit, '("   CHECKSUM:: ",A40," = ",Z20)') trim(var_name), &
      mpp_chksum(var%bc(n)%field(m)%values(var%isc:var%iec,var%jsc:var%jec) )
  enddo ; enddo

end subroutine CT_write_chksums_2d

!> This subroutine writes out checksums for the elements of a coupler_3d_bc_type
subroutine CT_write_chksums_3d(var, outunit, name_lead)
  type(coupler_3d_bc_type),   intent(in) :: var  !< BC_type structure for which to register diagnostics
  integer,                    intent(in) :: outunit !< The index of a open output file
  character(len=*), optional, intent(in) :: name_lead !< An optional prefix for the variable names

  character(len=120) :: var_name
  integer :: m, n

  do n = 1, var%num_bcs ; do m = 1, var%bc(n)%num_fields
    if (present(name_lead)) then
      var_name = trim(name_lead)//trim(var%bc(n)%field(m)%name)
    else
      var_name = trim(var%bc(n)%field(m)%name)
    endif
    write(outunit, '("   CHECKSUM:: ",A40," = ",Z20)') var_name, &
      mpp_chksum(var%bc(n)%field(m)%values(var%isc:var%iec,var%jsc:var%jec,:) )
  enddo ; enddo

end subroutine CT_write_chksums_3d


!> This function indicates whether a coupler_1d_bc_type has been initialized.
function CT_initialized_1d(var)
  type(coupler_1d_bc_type), intent(in) :: var  !< BC_type structure to be deconstructed
  logical :: CT_initialized_1d  !< The return value, indicating whether this type has been initialized

  CT_initialized_1d = var%set
end function CT_initialized_1d

!> This function indicates whether a coupler_2d_bc_type has been initialized.
function CT_initialized_2d(var)
  type(coupler_2d_bc_type), intent(in) :: var  !< BC_type structure to be deconstructed
  logical :: CT_initialized_2d  !< The return value, indicating whether this type has been initialized

  CT_initialized_2d = var%set
end function CT_initialized_2d

!> This function indicates whether a coupler_3d_bc_type has been initialized.
function CT_initialized_3d(var)
  type(coupler_3d_bc_type), intent(in) :: var  !< BC_type structure to be deconstructed
  logical :: CT_initialized_3d  !< The return value, indicating whether this type has been initialized

  CT_initialized_3d = var%set
end function CT_initialized_3d


!> This subroutine deallocates all data associated with a coupler_1d_bc_type
subroutine CT_destructor_1d(var)
  type(coupler_1d_bc_type), intent(inout) :: var  !< BC_type structure to be deconstructed

  integer :: m, n

  if (var%num_bcs > 0) then
    do n = 1, var%num_bcs
      do m = 1, var%bc(n)%num_fields
        deallocate ( var%bc(n)%field(m)%values )
      enddo
      deallocate ( var%bc(n)%field )
    enddo
    deallocate ( var%bc )
  endif

  var%num_bcs = 0 ; var%set = .false.

end subroutine CT_destructor_1d

!> This subroutine deallocates all data associated with a coupler_2d_bc_type
subroutine CT_destructor_2d(var)
  type(coupler_2d_bc_type), intent(inout) :: var  !< BC_type structure to be deconstructed

  integer :: m, n

  if (var%num_bcs > 0) then
    do n = 1, var%num_bcs
      do m = 1, var%bc(n)%num_fields
        deallocate ( var%bc(n)%field(m)%values )
      enddo
      deallocate ( var%bc(n)%field )
    enddo
    deallocate ( var%bc )
  endif

  var%num_bcs = 0 ; var%set = .false.

end subroutine CT_destructor_2d


!> This subroutine deallocates all data associated with a coupler_3d_bc_type
subroutine CT_destructor_3d(var)
  type(coupler_3d_bc_type), intent(inout) :: var  !< BC_type structure to be deconstructed

  integer :: m, n

  if (var%num_bcs > 0) then
    do n = 1, var%num_bcs
      do m = 1, var%bc(n)%num_fields
        deallocate ( var%bc(n)%field(m)%values )
      enddo
      deallocate ( var%bc(n)%field )
    enddo
    deallocate ( var%bc )
  endif

  var%num_bcs = 0 ; var%set = .false.

end subroutine CT_destructor_3d

end module coupler_types_mod
