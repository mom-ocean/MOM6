module FMS_coupler_util

use coupler_types_mod, only : coupler_2d_bc_type

implicit none ; private

public :: extract_coupler_values, set_coupler_values

contains

subroutine extract_coupler_values(BC_struc, BC_index, BC_element, array_out, ilb, jlb, &
                                  is, ie, js, je, conversion)
  real, dimension(ilb:,jlb:),intent(out) :: array_out
  integer,                   intent(in)  :: ilb, jlb
  type(coupler_2d_bc_type),  intent(in)  :: BC_struc
  integer,                   intent(in)  :: BC_index, BC_element
  integer,        optional,  intent(in)  :: is, ie, js, je
  real,           optional,  intent(in)  :: conversion
end subroutine extract_coupler_values

subroutine set_coupler_values(array_in, BC_struc, BC_index, BC_element, ilb, jlb,&
                              is, ie, js, je, conversion)
  real, dimension(ilb:,jlb:), intent(in)  :: array_in
  integer,                  intent(in)    :: ilb, jlb
  type(coupler_2d_bc_type), intent(inout) :: BC_struc
  integer,                  intent(in)    :: BC_index, BC_element
  integer,        optional, intent(in)    :: is, ie, js, je
  real,           optional, intent(in)    :: conversion
end subroutine set_coupler_values

end module FMS_coupler_util
