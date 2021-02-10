!> A wrapper for the FMS diag_manager routines. This module should be the
!! only MOM6 module which imports the FMS shared infrastructure for
!! diagnostics. Pass through interfaces are being documented
!! here and renamed in order to clearly identify these APIs as being
!! consistent with the FMS infrastructure (Any future updates to
!! those APIs would be applied here).
module MOM_diag_manager_infra

! This file is part of MOM6. See LICENSE.md for the license.

use diag_axis_mod,    only : fms_axis_init=>diag_axis_init
use diag_axis_mod,    only : fms_get_diag_axis_name => get_diag_axis_name
use diag_axis_mod,    only : EAST, NORTH
use diag_data_mod,    only : null_axis_id
use diag_manager_mod, only : fms_diag_manager_init => diag_manager_init
use diag_manager_mod, only : fms_diag_manager_end => diag_manager_end
use diag_manager_mod, only : send_data_fms => send_data
use diag_manager_mod, only : fms_diag_field_add_attribute => diag_field_add_attribute
use diag_manager_mod, only : DIAG_FIELD_NOT_FOUND
use diag_manager_mod, only : register_diag_field_fms => register_diag_field
use diag_manager_mod, only : register_static_field_fms => register_static_field
use diag_manager_mod, only : get_diag_field_id_fms => get_diag_field_id
use time_manager_mod, only : time_type
use MOM_domain_infra, only : MOM_domain_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING

implicit none ; private

!> transmit data for diagnostic output
interface register_diag_field_infra
  module procedure register_diag_field_infra_scalar
  module procedure register_diag_field_infra_array
end interface register_diag_field_infra

!> transmit data for diagnostic output
interface send_data_infra
  module procedure send_data_infra_0d, send_data_infra_1d
  module procedure send_data_infra_2d, send_data_infra_3d
#ifdef OVERLOAD_R8
  module procedure send_data_infra_2d_r8, send_data_infra_3d_r8
#endif
end interface send_data_infra

!> Add an attribute to a diagnostic field
interface MOM_diag_field_add_attribute
  module procedure MOM_diag_field_add_attribute_scalar_r
  module procedure MOM_diag_field_add_attribute_scalar_i
  module procedure MOM_diag_field_add_attribute_scalar_c
  module procedure MOM_diag_field_add_attribute_r1d
  module procedure MOM_diag_field_add_attribute_i1d
end interface MOM_diag_field_add_attribute


! Public interfaces
public MOM_diag_axis_init
public get_MOM_diag_axis_name
public MOM_diag_manager_init
public MOM_diag_manager_end
public send_data_infra
public MOM_diag_field_add_attribute
public register_diag_field_infra
public register_static_field_infra
public get_MOM_diag_field_id
! Public data
public null_axis_id
public DIAG_FIELD_NOT_FOUND
public EAST, NORTH


contains

!> Initialize a diagnostic axis
integer function MOM_diag_axis_init(name, data, units, cart_name, long_name, MOM_domain, position, &
          & direction, edges, set_name, coarsen, null_axis)
   character(len=*),   intent(in) :: name      !< The name of this axis
   real, dimension(:), intent(in) :: data      !< The array of coordinate values
   character(len=*),   intent(in) :: units     !< The units for the axis data
   character(len=*),   intent(in) :: cart_name !< Cartesian axis ("X", "Y", "Z", "T", or "N" for none)
   character(len=*), &
             optional, intent(in) :: long_name !< The long name of this axis
   type(MOM_domain_type), &
             optional, intent(in) :: MOM_Domain !< A MOM_Domain that describes the decomposition
   integer,  optional, intent(in) :: position  !< This indicates the relative position of this
                                               !! axis.  The default is CENTER, but EAST and NORTH
                                               !! are common options.
   integer,  optional, intent(in) :: direction !< This indicates the direction along which this
                                               !! axis increases: 1 for upward, -1 for downward, or
                                               !! 0 for non-vertical axes (the default)
   integer,  optional, intent(in) :: edges     !< The axis_id of the complementary axis that
                                               !! describes the edges of this axis
   character(len=*), &
             optional, intent(in) :: set_name  !< A name to use for this set of axes.
   integer,  optional, intent(in) :: coarsen   !< An optional degree of coarsening for the grid, 1
                                               !! by default.
   logical,  optional, intent(in) :: null_axis !< If present and true, return the special null axis
                                               !! id for use with scalars.

   integer :: coarsening ! The degree of grid coarsening

   if (present(null_axis)) then ; if (null_axis) then
     ! Return the special null axis id for scalars
     MOM_diag_axis_init = null_axis_id
     return
   endif ; endif

  if (present(MOM_domain)) then
    coarsening = 1 ; if (present(coarsen)) coarsening = coarsen
    if (coarsening == 1) then
      MOM_diag_axis_init = fms_axis_init(name, data, units, cart_name, long_name=long_name, &
              direction=direction, set_name=set_name, edges=edges, &
              domain2=MOM_domain%mpp_domain, domain_position=position)
    elseif (coarsening == 2) then
      MOM_diag_axis_init = fms_axis_init(name, data, units, cart_name, long_name=long_name, &
              direction=direction, set_name=set_name, edges=edges, &
              domain2=MOM_domain%mpp_domain_d2, domain_position=position)
    else
      call MOM_error(FATAL, "diag_axis_init called with an invalid value of coarsen.")
    endif
  else
    if (present(coarsen)) then ; if (coarsen /= 1) then
      call MOM_error(FATAL, "diag_axis_init does not support grid coarsening without a MOM_domain.")
    endif ; endif
    MOM_diag_axis_init = fms_axis_init(name, data, units, cart_name, long_name=long_name, &
            direction=direction, set_name=set_name, edges=edges)
  endif

end function MOM_diag_axis_init

!> Returns the short name of the axis
subroutine get_MOM_diag_axis_name(id, name)
  integer,          intent(in)  :: id   !< The axis numeric id
  character(len=*), intent(out) :: name !< The short name of the axis

  call fms_get_diag_axis_name(id, name)

end subroutine get_MOM_diag_axis_name

!> Return a unique numeric ID field a module/field name combination.
integer function get_MOM_diag_field_id(module_name, field_name)
  character(len=*), intent(in) :: module_name !< A module name string to query.
  character(len=*), intent(in) :: field_name  !< A field name string to query.


  get_MOM_diag_field_id = -1
  get_MOM_diag_field_id = get_diag_field_id_fms(module_name, field_name)

end function get_MOM_diag_field_id

!> Initializes the diagnostic manager
subroutine MOM_diag_manager_init(diag_model_subset, time_init, err_msg)
  integer,               optional, intent(in) :: diag_model_subset !< An optional diagnostic subset
  integer, dimension(6), optional, intent(in) :: time_init !< An optional reference time for diagnostics
                                                           !! The default uses the value contained in the
                                                           !! diag_table. Format is Y-M-D-H-M-S
  character(len=*),     optional, intent(out) :: err_msg   !< Error message.
  call FMS_diag_manager_init(diag_model_subset, time_init, err_msg)

end subroutine MOM_diag_manager_init

!> Close the diagnostic manager
subroutine MOM_diag_manager_end(time)
  type(time_type), intent(in) :: time !< Model time at call to close.

  call FMS_diag_manager_end(time)

end subroutine MOM_diag_manager_end

!> Register a MOM diagnostic field for scalars
integer function register_diag_field_infra_scalar(module_name, field_name, init_time, &
                        long_name, units, missing_value, range, standard_name, do_not_log, &
                        err_msg, area, volume)
  character(len=*),              intent(in) :: module_name !< The name of the associated module
  character(len=*),              intent(in) :: field_name !< The name of the field
  type(time_type),     optional, intent(in) :: init_time !< The registration time
  character(len=*),    optional, intent(in) :: long_name !< A long name for the field
  character(len=*),    optional, intent(in) :: units     !< Field units
  character(len=*),    optional, intent(in) :: standard_name !< A standard name for the field
  real,                optional, intent(in) :: missing_value !< Missing value attribute
  real,  dimension(2), optional, intent(in) :: range     !< A valid range of the field
  logical,             optional, intent(in) :: do_not_log !< if TRUE, field information is not logged
  character(len=*),    optional, intent(out):: err_msg   !< An error message to return
  integer,             optional, intent(in) :: area      !< Diagnostic ID of the field containing the area attribute
  integer,             optional, intent(in) :: volume    !< Diagnostic ID of the field containing the volume attribute

  register_diag_field_infra_scalar = register_diag_field_fms(module_name, field_name, init_time, &
        long_name, units, missing_value, range, standard_name, do_not_log, err_msg, area, volume)

end function register_diag_field_infra_scalar

!> Register a MOM diagnostic field for scalars
integer function register_diag_field_infra_array(module_name, field_name, axes, init_time, &
                        long_name, units, missing_value, range, mask_variant, standard_name, verbose, &
                        do_not_log, err_msg, interp_method, tile_count, area, volume)
  character(len=*),             intent(in) :: module_name !< The name of the associated module
  character(len=*),             intent(in) :: field_name !< The name of the field
  integer, dimension(:),        intent(in) :: axes      !< Diagnostic IDs of axis attributes for the field
  type(time_type),    optional, intent(in) :: init_time !< The registration time
  character(len=*),   optional, intent(in) :: long_name !< A long name for the field
  character(len=*),   optional, intent(in) :: units     !< Units of the field
  real,               optional, intent(in) :: missing_value !< Missing value attribute
  real, dimension(2), optional, intent(in) :: range     !< A valid range of the field
  logical,            optional, intent(in) :: mask_variant !< If true, the field mask is varying in time
  character(len=*),   optional, intent(in) :: standard_name !< A standard name for the field
  logical,            optional, intent(in) :: verbose    !< If true, provide additional log information
  logical,            optional, intent(in) :: do_not_log !< if TRUE, field information is not logged
  character(len=*),   optional, intent(in) :: interp_method !< If 'none' indicates the field should
                                                         !! not be interpolated as a scalar
  integer,            optional, intent(in) :: tile_count !< The tile number for the current PE
  character(len=*),   optional, intent(out):: err_msg   !< An error message to return
  integer,            optional, intent(in) :: area      !< Diagnostic ID of the field containing the area attribute
  integer,            optional, intent(in) :: volume    !< Diagnostic ID of the field containing the volume attribute

  register_diag_field_infra_array = register_diag_field_fms(module_name, field_name, axes, init_time, &
        long_name, units, missing_value, range, mask_variant, standard_name, verbose, do_not_log, &
        err_msg, interp_method, tile_count, area, volume)

end function register_diag_field_infra_array


integer function register_static_field_infra(module_name, field_name, axes, long_name, units, &
                        missing_value, range, mask_variant, standard_name, do_not_log, interp_method, &
                        tile_count, area, volume)
  character(len=*),             intent(in) :: module_name !< The name of the associated module
  character(len=*),             intent(in) :: field_name !< The name of the field
  integer, dimension(:),        intent(in) :: axes      !< Diagnostic IDs of axis attributes for the field
  character(len=*),   optional, intent(in) :: long_name !< A long name for the field
  character(len=*),   optional, intent(in) :: units     !< Units of the field
  real,               optional, intent(in) :: missing_value !< Missing value attribute
  real, dimension(2), optional, intent(in) :: range     !< A valid range of the field
  logical,            optional, intent(in) :: mask_variant !< If true, the field mask is varying in time
  character(len=*),   optional, intent(in) :: standard_name !< A standard name for the field
  logical,            optional, intent(in) :: do_not_log !< if TRUE, field information is not logged
  character(len=*),   optional, intent(in) :: interp_method !< If 'none' indicates the field should
                                                         !! not be interpolated as a scalar
  integer,            optional, intent(in) :: tile_count !< The tile number for the current PE
  integer,            optional, intent(in) :: area      !< Diagnostic ID of the field containing the area attribute
  integer,            optional, intent(in) :: volume    !< Diagnostic ID of the field containing the volume attribute

  register_static_field_infra = register_static_field_fms(module_name, field_name, axes, long_name, units,&
       & missing_value, range, mask_variant, standard_name, dynamic=.false.,do_not_log=do_not_log, &
       interp_method=interp_method,tile_count=tile_count, area=area, volume=volume)
end function register_static_field_infra

!> Returns true if the argument data are successfully passed to a diagnostic manager
!! with the indicated unique reference id, false otherwise.
logical function send_data_infra_0d(diag_field_id, field, time, err_msg)
  integer,                    intent(in)  :: diag_field_id !< The diagnostic manager identifier for this field
  real,                       intent(in)  :: field   !< The value being recorded
  TYPE(time_type),  optional, intent(in)  :: time    !< The time for the current record
  CHARACTER(len=*), optional, intent(out) :: err_msg !< An optional error message

  send_data_infra_0d = send_data_fms(diag_field_id, field, time, err_msg)
end function send_data_infra_0d

!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_infra_1d(diag_field_id, field, is_in, ie_in, time, mask, rmask, weight, err_msg)
  integer,                         intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  real, dimension(:),              intent(in) :: field !< A 1-d array of values being recorded
  integer,               optional, intent(in) :: is_in !< The starting index for the data being recorded
  integer,               optional, intent(in) :: ie_in !< The end index for the data being recorded
  type(time_type),       optional, intent(in) :: time  !< The time for the current record
  logical, dimension(:), optional, intent(in) :: mask  !< An optional rank 1 logical mask
  real, dimension(:),    optional, intent(in) :: rmask !< An optional rank 1 mask array
  real,                  optional, intent(in) :: weight !< A scalar weight factor to apply to the current
                                                       !! record if there is averaging in time
  character(len=*),      optional, intent(out) :: err_msg !< A log indicating the status of the post upon
                                                       !! returning to the calling routine

  send_data_infra_1d = send_data_fms(diag_field_id, field, time, is_in, mask, rmask, ie_in, weight, err_msg)

end function send_data_infra_1d

!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_infra_2d(diag_field_id, field, is_in, ie_in, js_in, je_in, &
                                    time, mask, rmask, weight, err_msg)
  integer,                           intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  real, dimension(:,:),              intent(in) :: field !< A 2-d array of values being recorded
  integer,                 optional, intent(in) :: is_in !< The starting i-index for the data being recorded
  integer,                 optional, intent(in) :: ie_in !< The end i-index for the data being recorded
  integer,                 optional, intent(in) :: js_in !< The starting j-index for the data being recorded
  integer,                 optional, intent(in) :: je_in !< The end j-index for the data being recorded
  type(time_type),         optional, intent(in) :: time  !< The time for the current record
  logical, dimension(:,:), optional, intent(in) :: mask  !< An optional 2-d logical mask
  real, dimension(:,:),    optional, intent(in) :: rmask !< An optional 2-d mask array
  real,                    optional, intent(in) :: weight !< A scalar weight factor to apply to the current
                                                         !! record if there is averaging in time
  character(len=*),        optional, intent(out) :: err_msg !< A log indicating the status of the post upon
                                                         !! returning to the calling routine

  send_data_infra_2d = send_data_fms(diag_field_id, field, time, is_in, js_in, mask, &
                                rmask, ie_in, je_in, weight, err_msg)

end function send_data_infra_2d

!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_infra_3d(diag_field_id, field, is_in, ie_in, js_in, je_in, ks_in, ke_in, &
                                    time, mask, rmask, weight, err_msg)
  integer,                             intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  real, dimension(:,:,:),              intent(in) :: field !< A rank 1 array of floating point values being recorded
  integer,                   optional, intent(in) :: is_in !< The starting i-index for the data being recorded
  integer,                   optional, intent(in) :: ie_in !< The end i-index for the data being recorded
  integer,                   optional, intent(in) :: js_in !< The starting j-index for the data being recorded
  integer,                   optional, intent(in) :: je_in !< The end j-index for the data being recorded
  integer,                   optional, intent(in) :: ks_in !< The starting k-index for the data being recorded
  integer,                   optional, intent(in) :: ke_in !< The end k-index for the data being recorded
  type(time_type),           optional, intent(in) :: time  !< The time for the current record
  logical, dimension(:,:,:), optional, intent(in) :: mask  !< An optional 3-d logical mask
  real, dimension(:,:,:),    optional, intent(in) :: rmask !< An optional 3-d mask array
  real,                      optional, intent(in) :: weight !< A scalar weight factor to apply to the current
                                                           !! record if there is averaging in time
  character(len=*),          optional, intent(out) :: err_msg !< A log indicating the status of the post upon
                                                           !! returning to the calling routine

  send_data_infra_3d = send_data_fms(diag_field_id, field, time, is_in, js_in, ks_in, mask, &
                               rmask, ie_in, je_in, ke_in, weight, err_msg)

end function send_data_infra_3d


#ifdef OVERLOAD_R8
!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_infra_2d_r8(diag_field_id, field, is_in, ie_in, js_in, je_in, &
                                       time, mask, rmask, weight, err_msg)
  integer,                           intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  real(kind=8), dimension(:,:),      intent(in) :: field !< A 2-d array of values being recorded
  integer,                 optional, intent(in) :: is_in !< The starting i-index for the data being recorded
  integer,                 optional, intent(in) :: ie_in !< The end i-index for the data being recorded
  integer,                 optional, intent(in) :: js_in !< The starting j-index for the data being recorded
  integer,                 optional, intent(in) :: je_in !< The end j-index for the data being recorded
  type(time_type),         optional, intent(in) :: time  !< The time for the current record
  logical, dimension(:,:), optional, intent(in) :: mask  !< An optional 2-d logical mask
  real, dimension(:,:),    optional, intent(in) :: rmask !< An optional 2-d mask array
  real,                    optional, intent(in) :: weight !< A scalar weight factor to apply to the current
                                                         !! record if there is averaging in time
  character(len=*),        optional, intent(out) :: err_msg !< A log indicating the status of the post upon
                                                         !! returning to the calling routine

  send_data_infra_2d_r8 = send_data_fms(diag_field_id, field, time, is_in, js_in, mask, &
                                   rmask, ie_in, je_in, weight, err_msg)

end function send_data_infra_2d_r8

!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_infra_3d_r8(diag_field_id, field, is_in, ie_in, js_in, je_in, ks_in, ke_in, &
                                    time, mask, rmask, weight, err_msg)
  integer,                             intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  real(kind=8), dimension(:,:,:),      intent(in) :: field !< A rank 1 array of floating point values being recorded
  integer,                   optional, intent(in) :: is_in !< The starting i-index for the data being recorded
  integer,                   optional, intent(in) :: ie_in !< The end i-index for the data being recorded
  integer,                   optional, intent(in) :: js_in !< The starting j-index for the data being recorded
  integer,                   optional, intent(in) :: je_in !< The end j-index for the data being recorded
  integer,                   optional, intent(in) :: ks_in !< The starting k-index for the data being recorded
  integer,                   optional, intent(in) :: ke_in !< The end k-index for the data being recorded
  type(time_type),           optional, intent(in) :: time  !< The time for the current record
  logical, dimension(:,:,:), optional, intent(in) :: mask  !< An optional 3-d logical mask
  real, dimension(:,:,:),    optional, intent(in) :: rmask !< An optional 3-d mask array
  real,                      optional, intent(in) :: weight !< A scalar weight factor to apply to the current
                                                           !! record if there is averaging in time
  character(len=*),          optional, intent(out) :: err_msg !< A log indicating the status of the post upon
                                                           !! returning to the calling routine

  send_data_infra_3d_r8 = send_data_fms(diag_field_id, field, time, is_in, js_in, ks_in, mask, rmask, &
                                ie_in, je_in, ke_in, weight, err_msg)

end function send_data_infra_3d_r8
#endif

!> Add a real scalar attribute to a diagnostic field
subroutine MOM_diag_field_add_attribute_scalar_r(diag_field_id, att_name, att_value)
  integer,          intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  character(len=*), intent(in) :: att_name  !< The name of the attribute
  real,             intent(in) :: att_value !< A real scalar value

  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end subroutine MOM_diag_field_add_attribute_scalar_r

!> Add an integer attribute to a diagnostic field
subroutine MOM_diag_field_add_attribute_scalar_i(diag_field_id, att_name, att_value)
  integer,          intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  character(len=*), intent(in) :: att_name  !< The name of the attribute
  integer,          intent(in) :: att_value !< An integer scalar value

  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end subroutine MOM_diag_field_add_attribute_scalar_i

!> Add a character string attribute to a diagnostic field
subroutine MOM_diag_field_add_attribute_scalar_c(diag_field_id, att_name, att_value)
  integer,          intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  character(len=*), intent(in) :: att_name  !< The name of the attribute
  character(len=*), intent(in) :: att_value !< A character string value

  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end subroutine MOM_diag_field_add_attribute_scalar_c

!> Add a real list of attributes attribute to a diagnostic field
subroutine MOM_diag_field_add_attribute_r1d(diag_field_id, att_name, att_value)
  integer,            intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  character(len=*),   intent(in) :: att_name  !< The name of the attribute
  real, dimension(:), intent(in) :: att_value !< An array of real values

  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end subroutine MOM_diag_field_add_attribute_r1d

!> Add a integer list of attributes attribute to a diagnostic field
subroutine MOM_diag_field_add_attribute_i1d(diag_field_id, att_name, att_value)
  integer,               intent(in) :: diag_field_id !< The diagnostic manager identifier for this field
  character(len=*),      intent(in) :: att_name  !< The name of the attribute
  integer, dimension(:), intent(in) :: att_value !< An array of integer values

  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end subroutine MOM_diag_field_add_attribute_i1d

end module MOM_diag_manager_infra
