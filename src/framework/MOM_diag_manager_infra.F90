!> A simple (very thin) wrapper for the FMS diag_manager routines, with some name changes
module MOM_diag_manager_infra

! This file is part of MOM6. See LICENSE.md for the license.

use diag_axis_mod,    only : axis_init=>diag_axis_init
use diag_axis_mod,    only : FMS_get_diag_axis_name => get_diag_axis_name
use diag_axis_mod,    only : EAST, NORTH
use diag_data_mod,    only : null_axis_id
use diag_manager_mod, only : FMS_diag_manager_init => diag_manager_init
use diag_manager_mod, only : FMS_diag_manager_end => diag_manager_end
use diag_manager_mod, only : send_data_fms => send_data
use diag_manager_mod, only : FMS_diag_field_add_attribute => diag_field_add_attribute
use diag_manager_mod, only : DIAG_FIELD_NOT_FOUND
use diag_manager_mod, only : register_diag_field_fms => register_diag_field
use diag_manager_mod, only : register_static_field_fms => register_static_field
use diag_manager_mod, only : get_diag_field_id_fms => get_diag_field_id
use time_manager_mod, only : time_type
use MOM_domain_infra, only : MOM_domain_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING
implicit none ; private

!> transmit data for diagnostic output
interface send_data
  module procedure send_data_0d
  module procedure send_data_1d
  module procedure send_data_2d
  module procedure send_data_3d
#ifdef OVERLOAD_R8
  module procedure send_data_2d_r8
  module procedure send_data_3d_r8
#endif
end interface send_data

!> Add an attribute to a diagnostic field
interface diag_field_add_attribute
  module procedure diag_field_add_attribute_scalar_r
  module procedure diag_field_add_attribute_scalar_i
  module procedure diag_field_add_attribute_scalar_c
  module procedure diag_field_add_attribute_r1d
  module procedure diag_field_add_attribute_i1d
end interface diag_field_add_attribute


! Public interfaces
public diag_axis_init
public get_diag_axis_name
public diag_manager_init
public diag_manager_end
public send_data
public diag_field_add_attribute
public register_diag_field_fms
public register_static_field_fms
public get_diag_field_id_fms
! Public data
public null_axis_id
public DIAG_FIELD_NOT_FOUND
public EAST, NORTH


contains

!> Initialize a diagnostic axis
integer function diag_axis_init(name, data, units, cart_name, long_name, MOM_domain, position, &
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
     diag_axis_init = null_axis_id
     return
   endif ; endif

  if (present(MOM_domain)) then
    coarsening = 1 ; if (present(coarsen)) coarsening = coarsen
    if (coarsening == 1) then
      diag_axis_init = axis_init(name, data, units, cart_name, long_name=long_name, &
              direction=direction, set_name=set_name, edges=edges, &
              domain2=MOM_domain%mpp_domain, domain_position=position)
    elseif (coarsening == 2) then
      diag_axis_init = axis_init(name, data, units, cart_name, long_name=long_name, &
              direction=direction, set_name=set_name, edges=edges, &
              domain2=MOM_domain%mpp_domain_d2, domain_position=position)
    else
      call MOM_error(FATAL, "diag_axis_init called with an invalid value of coarsen.")
    endif
  else
    if (present(coarsen)) then ; if (coarsen /= 1) then
      call MOM_error(FATAL, "diag_axis_init does not support grid coarsening without a MOM_domain.")
    endif ; endif
    diag_axis_init = axis_init(name, data, units, cart_name, long_name=long_name, &
            direction=direction, set_name=set_name, edges=edges)
  endif

end function diag_axis_init

!> Returns the short name of the axis
subroutine get_diag_axis_name(id, name)
  integer,  intent(in)          :: id   !< The axis numeric id
  character(len=*), intent(out) :: name !< The short name of the axis

  call FMS_get_diag_axis_name(id, name)

end subroutine get_diag_axis_name

!> Initializes the diagnostic manager
subroutine diag_manager_init(diag_model_subset, time_init, err_msg)
  integer, optional, intent(in) :: diag_model_subset       !< An optional diagnostic subset
  integer, dimension(6), optional, intent(in) :: time_init !< An optional reference time for diagnostics
                                                           !! The default uses the value contained in the
                                                           !! diag_table. Format is Y-M-D-H-M-S
  character(len=*), intent(out), optional :: err_msg        !< Error message.
  call FMS_diag_manager_init(diag_model_subset, time_init, err_msg)

end subroutine diag_manager_init

!> Close the diag manager
subroutine diag_manager_end(time)
  type(time_type), intent(in) :: time !< Model time at call to close.

  call FMS_diag_manager_end(time)

end subroutine diag_manager_end

!> Returns true if the argument data are successfully passed to a diagnostic manager
!! with the indicated unique reference id, false otherwise.
logical function send_data_0d(diag_field_id, field, time, err_msg)
  integer, intent(in) :: diag_field_id !< A unique identifier for this data to the diagnostic manager
  real, intent(in) :: field !< Floating point value being recorded
  TYPE(time_type), intent(in), optional :: time !< Time slice for this record
  CHARACTER(len=*), intent(out), optional :: err_msg !< An optional error message

  send_data_0d= send_data_fms(diag_field_id, field, time, err_msg)
end function send_data_0d

!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_1d(diag_field_id, field, time, is_in, mask, rmask, ie_in, weight, err_msg)
    integer, intent(in) :: diag_field_id !< A unique identifier for this data to the diagnostic manager
    real, dimension(:), intent(in) :: field !< A rank 1 array of floating point values being recorded
    type (time_type), intent(in), optional :: time  !< The time for the current record.
    logical, intent(in), dimension(:), optional :: mask !< An optional rank 1 logical mask.
    real, intent(in), dimension(:), optional :: rmask !< An optional rank 1 mask array
    integer, intent(in), optional :: is_in, ie_in !< An optional range for subsetting the data being recorded.
    real, intent(in), optional :: weight !< An optional scalar weight factor to apply to the current record
                                         !! in the case where data a data reduction in time is being performed.
    character(len=*), intent(out), optional :: err_msg !< A log indicating the status of the post upon
                                                       !! returning to the calling routine.

    send_data_1d= send_data_fms(diag_field_id, field, time, is_in, mask, rmask, ie_in, weight, err_msg)

end function send_data_1d

!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_2d(diag_field_id, field, time, is_in, js_in, mask, rmask, &
                            &  ie_in, je_in, weight, err_msg)
    integer, intent(in) :: diag_field_id !< A unique identifier for this data to the diagnostic manager
    real, dimension(:,:), intent(in) :: field !< A rank 1 array of floating point values being recorded
    type (time_type), intent(in), optional :: time  !< The time for the current record.
    logical, intent(in), dimension(:,:), optional :: mask !< An optional rank 1 logical mask.
    real, intent(in), dimension(:,:), optional :: rmask !< An optional rank 1 mask array
    integer, intent(in), optional :: is_in, ie_in !< An optional range for subsetting the data being recorded.
    integer, intent(in), optional :: js_in, je_in !< An optional range for subsetting the data being recorded.
    real, intent(in), optional :: weight !< An optional scalar weight factor to apply to the current record
                                         !! in the case where data a data reduction in time is being performed.
    character(len=*), intent(out), optional :: err_msg !< A log indicating the status of the post upon
                                                       !! returning to the calling routine.

    send_data_2d= send_data_fms(diag_field_id, field, time, is_in, js_in, mask, &
                                rmask, ie_in, je_in, weight, err_msg)

end function send_data_2d

#ifdef OVERLOAD_R8
!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_2d_r8(diag_field_id, field, time, is_in, js_in, mask, rmask, &
                            &  ie_in, je_in, weight, err_msg)
    integer, intent(in) :: diag_field_id !< A unique identifier for this data to the diagnostic manager
    real(kind=8), dimension(:,:), intent(in) :: field !< A rank 1 array of floating point values being recorded
    type (time_type), intent(in), optional :: time  !< The time for the current record.
    logical, intent(in), dimension(:,:), optional :: mask !< An optional rank 1 logical mask.
    real, intent(in), dimension(:,:), optional :: rmask !< An optional rank 1 mask array
    integer, intent(in), optional :: is_in, ie_in !< An optional range for subsetting the data being recorded.
    integer, intent(in), optional :: js_in, je_in !< An optional range for subsetting the data being recorded.
    real, intent(in), optional :: weight !< An optional scalar weight factor to apply to the current record
                                         !! in the case where data a data reduction in time is being performed.
    character(len=*), intent(out), optional :: err_msg !< A log indicating the status of the post upon
                                                       !! returning to the calling routine.

    send_data_2d_r8 = send_data_fms(diag_field_id, field, time, is_in, js_in, mask, &
                                   rmask, ie_in, je_in, weight, err_msg)

end function send_data_2d_r8
#endif

!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_3d(diag_field_id, field, time, is_in, js_in, ks_in, mask, rmask, &
                            &  ie_in, je_in, ke_in, weight, err_msg)
    integer, intent(in) :: diag_field_id !< A unique identifier for this data to the diagnostic manager
    real, dimension(:,:,:), intent(in) :: field !< A rank 1 array of floating point values being recorded
    type (time_type), intent(in), optional :: time  !< The time for the current record.
    logical, intent(in), dimension(:,:,:), optional :: mask !< An optional rank 1 logical mask.
    real, intent(in), dimension(:,:,:), optional :: rmask !< An optional rank 1 mask array
    integer, intent(in), optional :: is_in, ie_in !< An optional range for subsetting the data being recorded.
    integer, intent(in), optional :: js_in, je_in !< An optional range for subsetting the data being recorded.
    integer, intent(in), optional :: ks_in, ke_in !< An optional range for subsetting the data being recorded.
    real, intent(in), optional :: weight !< An optional scalar weight factor to apply to the current record
                                         !! in the case where data a data reduction in time is being performed.
    character(len=*), intent(out), optional :: err_msg !< A log indicating the status of the post upon
                                                       !! returning to the calling routine.

    send_data_3d = send_data_fms(diag_field_id, field, time, is_in, js_in, ks_in, mask, &
                               rmask, ie_in, je_in, ke_in, weight, err_msg)

end function send_data_3d


#ifdef OVERLOAD_R8
!> Returns true if the argument data are successfully passed to a diagnostic manager
!!  with the indicated unique reference id, false otherwise.
logical function send_data_3d_r8(diag_field_id, field, time, is_in, js_in, ks_in, mask, rmask, &
                            &  ie_in, je_in, ke_in, weight, err_msg)
    integer, intent(in) :: diag_field_id !< A unique identifier for this data to the diagnostic manager
    real(kind=8), dimension(:,:,:), intent(in) :: field !< A rank 1 array of floating point values being recorded
    type (time_type), intent(in), optional :: time  !< The time for the current record.
    logical, intent(in), dimension(:,:,:), optional :: mask !< An optional rank 1 logical mask.
    real, intent(in), dimension(:,:,:), optional :: rmask !< An optional rank 1 mask array
    integer, intent(in), optional :: is_in, ie_in !< An optional range for subsetting the data being recorded.
    integer, intent(in), optional :: js_in, je_in !< An optional range for subsetting the data being recorded.
    integer, intent(in), optional :: ks_in, ke_in !< An optional range for subsetting the data being recorded.
    real, intent(in), optional :: weight !< An optional scalar weight factor to apply to the current record
                                         !! in the case where data a data reduction in time is being performed.
    character(len=*), intent(out), optional :: err_msg !< A log indicating the status of the post upon
                                                       !! returning to the calling routine.

    send_data_3d_r8 = send_data_fms(diag_field_id, field, time, is_in, js_in, ks_in, mask, rmask, &
                                ie_in, je_in, ke_in, weight, err_msg)

end function send_data_3d_r8
#endif

!> Add a real scalar attribute to a diagnostic field
subroutine diag_field_add_attribute_scalar_r(diag_field_id, att_name, att_value)
  integer, intent(in) :: diag_field_id !< A unique numeric field id
  character(len=*), intent(in) :: att_name !< The name of the attribute
  real, intent(in) :: att_value !< A real scalar value

  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end subroutine diag_field_add_attribute_scalar_r

!> Add an integer attribute to a diagnostic field
subroutine diag_field_add_attribute_scalar_i(diag_field_id, att_name, att_value)
  integer, intent(in) :: diag_field_id !< A unique numeric field id
  character(len=*), intent(in) :: att_name !< The name of the attribute
  integer, intent(in) :: att_value !< A real scalar value

  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end subroutine diag_field_add_attribute_scalar_i

!> Add a character string attribute to a diagnostic field
subroutine diag_field_add_attribute_scalar_c(diag_field_id, att_name, att_value)
  integer, intent(in) :: diag_field_id !< A unique numeric field id
  character(len=*), intent(in) :: att_name !< The name of the attribute
  character(len=*), intent(in) :: att_value !< A real scalar value

  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end subroutine diag_field_add_attribute_scalar_c

!> Add a real list of attributes attribute to a diagnostic field
subroutine diag_field_add_attribute_r1d(diag_field_id, att_name, att_value)
  integer, intent(in) :: diag_field_id !< A unique numeric field id
  character(len=*), intent(in) :: att_name !< The name of the attribute
  real, dimension(:), intent(in) :: att_value !< A real scalar value

  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end subroutine diag_field_add_attribute_r1d

!> Add a integer list of attributes attribute to a diagnostic field
subroutine diag_field_add_attribute_i1d(diag_field_id, att_name, att_value)
  integer, intent(in) :: diag_field_id !< A unique numeric field id
  character(len=*), intent(in) :: att_name !< The name of the attribute
  integer, dimension(:), intent(in) :: att_value !< A integer list of values

  call FMS_diag_field_add_attribute(diag_field_id, att_name, att_value)

end subroutine diag_field_add_attribute_i1d


end module MOM_diag_manager_infra
