!> A simple (very thin) wrapper for the FMS diag_manager routines, with some name changes
module MOM_diag_manager

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_time_manager, only : time_type
use diag_axis_mod, only : diag_axis_init, get_diag_axis_name, EAST, NORTH
use diag_data_mod, only : null_axis_id
use diag_manager_mod, only : diag_manager_init, diag_manager_end
use diag_manager_mod, only : send_data, diag_field_add_attribute, DIAG_FIELD_NOT_FOUND
use diag_manager_mod, only : register_diag_field
use diag_manager_mod, only : register_static_field_fms=>register_static_field
use diag_manager_mod, only : get_diag_field_id_fms=>get_diag_field_id

implicit none ; private

public diag_manager_init, diag_manager_end
public diag_axis_init, get_diag_axis_name, EAST, NORTH, null_axis_id
public send_data, diag_field_add_attribute, DIAG_FIELD_NOT_FOUND
public register_diag_field_fms, register_static_field_fms,  get_diag_field_id_fms

!> A wrapper for register_diag_field_array()
interface register_diag_field_fms
  module procedure register_diag_field_array_fms, register_diag_field_scalar_fms
end interface

contains

!> An integer handle for a diagnostic array returned by register_diag_field()
integer function register_diag_field_array_fms(module_name, field_name, axes, init_time, &
                     long_name, units, missing_value, range, mask_variant, standard_name, &
                     verbose, do_not_log, err_msg, interp_method, tile_count, area, volume)
  character(len=*), intent(in) :: module_name             !< Name of this module, usually "ocean_model" or
                                                          !! "ice_shelf_model"
  character(len=*), intent(in) :: field_name              !< Name of the diagnostic field
  integer,          intent(in) :: axes(:)                 !< Container w/ up to 3 integer handles that
                                                          !! indicates axes for this field
  type(time_type),  intent(in) :: init_time               !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name     !< Long name of a field.
  character(len=*), optional, intent(in) :: units         !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values.
  real,             optional, intent(in) :: range(2)      !< Valid range of a variable (not used in MOM?)
  logical,          optional, intent(in) :: mask_variant  !< If true a logical mask must be provided with
                                                          !! post_data calls (not used in MOM?)
  logical,          optional, intent(in) :: verbose       !< If true, FMS is verbose (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log    !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg       !< String into which an error message might be
                                                          !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not be
                                                          !! interpolated as a scalar
  integer,          optional, intent(in) :: tile_count    !< no clue (not used in MOM?)
  integer,          optional, intent(in) :: area          !< The FMS id of cell area
  integer,          optional, intent(in) :: volume        !< The FMS id of cell volume
  ! Local variables

  register_diag_field_array_fms = register_diag_field(module_name, field_name, axes,   &
             init_time, long_name=long_name, units=units, missing_value=missing_value, &
             mask_variant=mask_variant, standard_name=standard_name,                   &
             verbose=verbose, do_not_log=do_not_log, err_msg=err_msg,                  &
             area=area, volume=volume, interp_method=interp_method)

end function register_diag_field_array_fms

!> An integer handle for a diagnostic scalar array returned by register_diag_field()
integer function register_diag_field_scalar_fms(module_name, field_name, init_time, &
                     long_name, units, missing_value, range, mask_variant, standard_name, &
                     verbose, do_not_log, err_msg, interp_method, tile_count, area, volume)
  character(len=*), intent(in) :: module_name             !< Name of this module, usually "ocean_model"
                                                          !! or "ice_shelf_model"
  character(len=*), intent(in) :: field_name              !< Name of the diagnostic field
  type(time_type),  intent(in) :: init_time               !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name     !< Long name of a field.
  character(len=*), optional, intent(in) :: units         !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values.
  real,             optional, intent(in) :: range(2)      !< Valid range of a variable (not used in MOM?)
  logical,          optional, intent(in) :: mask_variant  !< If true a logical mask must be provided with
                                                          !! post_data calls (not used in MOM?)
  logical,          optional, intent(in) :: verbose       !< If true, FMS is verbose (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log    !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg       !< String into which an error message might
                                                          !! be placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                          !! be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count    !< no clue (not used in MOM?)
  integer,          optional, intent(in) :: area          !< The FMS id of cell area (not used for scalars)
  integer,          optional, intent(in) :: volume        !< The FMS id of cell volume (not used for scalars)
  ! Local variables

  register_diag_field_scalar_fms = register_diag_field(module_name, field_name,        &
             init_time, long_name=long_name, units=units, missing_value=missing_value, &
             standard_name=standard_name, do_not_log=do_not_log, err_msg=err_msg)

end function register_diag_field_scalar_fms

!> \namespace mom_diag_manager
!!
!! This module simply wraps register_diag_field() from FMS's diag_manager_mod.
!! We used to be able to import register_diag_field and rename it to register_diag_field_fms
!! with a simple "use, only : register_diag_field_fms => register_diag_field" but PGI 16.5
!! has a bug that refuses to compile this - earlier versions did work.

end module MOM_diag_manager
