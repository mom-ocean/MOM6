!> Provides a couple of interfaces to allow more transparent and
!! robust extraction of the various fields in the coupler types.
module coupler_util

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL, WARNING
use coupler_types_mod, only : coupler_2d_bc_type, ind_flux, ind_alpha
use coupler_types_mod, only : ind_csurf

implicit none ; private

public :: extract_coupler_values, set_coupler_values
public :: ind_flux, ind_alpha, ind_csurf

contains

!> Extract an array of values in a coupler bc type
subroutine extract_coupler_values(BC_struc, BC_index, BC_element, array_out, &
                                  is, ie, js, je, conversion)
  type(coupler_2d_bc_type), intent(in)  :: BC_struc !< The type from which the data is being extracted.
  integer,                  intent(in)  :: BC_index !< The boundary condition number being extracted.
  integer,                  intent(in)  :: BC_element !< The element of the boundary condition being extracted.
  real, dimension(:,:),     intent(out) :: array_out !< The array being filled with the input values.
  integer,        optional, intent(in)  :: is !< Start i-index
  integer,        optional, intent(in)  :: ie !< End i-index
  integer,        optional, intent(in)  :: js !< Start j-index
  integer,        optional, intent(in)  :: je !< End j-index
  real,           optional, intent(in)  :: conversion !< A number that every element is multiplied by, to
                                           !! permit sign convention or unit conversion.
  ! Local variables
  real, pointer, dimension(:,:) :: Array_in
  real :: conv
  integer :: i, j, is0, ie0, js0, je0, i_offset, j_offset

  if ((BC_element /= ind_flux) .and. (BC_element /= ind_alpha) .and. &
      (BC_element /= ind_csurf)) then
    call MOM_error(FATAL,"extract_coupler_values: Unrecognized BC_element.")
  endif

  ! These error messages should be made more explicit.
!  if (.not.associated(BC_struc%bc(BC_index))) &
  if (.not.associated(BC_struc%bc)) &
    call MOM_error(FATAL,"extract_coupler_values: " // &
       "The requested boundary condition is not associated.")
!  if (.not.associated(BC_struc%bc(BC_index)%field(BC_element))) &
  if (.not.associated(BC_struc%bc(BC_index)%field)) &
    call MOM_error(FATAL,"extract_coupler_values: " // &
       "The requested boundary condition element is not associated.")
  if (.not.associated(BC_struc%bc(BC_index)%field(BC_element)%values)) &
    call MOM_error(FATAL,"extract_coupler_values: " // &
       "The requested boundary condition value array is not associated.")

  Array_in => BC_struc%bc(BC_index)%field(BC_element)%values

  if (present(is)) then ; is0 = is ; else ; is0 = LBOUND(array_out,1) ; endif
  if (present(ie)) then ; ie0 = ie ; else ; ie0 = UBOUND(array_out,1) ; endif
  if (present(js)) then ; js0 = js ; else ; js0 = LBOUND(array_out,2) ; endif
  if (present(je)) then ; je0 = je ; else ; je0 = UBOUND(array_out,2) ; endif

  conv = 1.0 ; if (present(conversion)) conv = conversion

  if (size(Array_in,1) /= ie0 - is0 + 1) &
    call MOM_error(FATAL,"extract_coupler_values: Mismatch in i-size " // &
                   "between BC array and output array or computational domain.")
  if (size(Array_in,2) /= je0 - js0 + 1) &
    call MOM_error(FATAL,"extract_coupler_values: Mismatch in i-size " // &
                   "between BC array and output array or computational domain.")
  i_offset = lbound(Array_in,1) - is0
  j_offset = lbound(Array_in,2) - js0
  do j=js0,je0 ; do i=is0,ie0
    array_out(i,j) = conv * Array_in(i+i_offset,j+j_offset)
  enddo ; enddo

end subroutine extract_coupler_values

!> Set an array of values in a coupler bc type
subroutine set_coupler_values(array_in, BC_struc, BC_index, BC_element, &
                              is, ie, js, je, conversion)
  real, dimension(:,:),     intent(in)    :: array_in !< The array containing the values to load into the BC.
  type(coupler_2d_bc_type), intent(inout) :: BC_struc !< The type from which the data is being extracted.
  integer,                  intent(in)    :: BC_index !< The boundary condition number being extracted.
  integer,                  intent(in)    :: BC_element !< The element of the boundary condition being extracted.
                                             !! This could be ind_csurf, ind_alpha, ind_flux or ind_deposition.
  integer,        optional, intent(in)    :: is !< Start i-index
  integer,        optional, intent(in)    :: ie !< End i-index
  integer,        optional, intent(in)    :: js !< Start j-index
  integer,        optional, intent(in)    :: je !< End j-index
  real,           optional, intent(in)    :: conversion !< A number that every element is multiplied by, to
                                             !! permit sign convention or unit conversion.
  ! Local variables
  real, pointer, dimension(:,:) :: Array_out
  real :: conv
  integer :: i, j, is0, ie0, js0, je0, i_offset, j_offset

  if ((BC_element /= ind_flux) .and. (BC_element /= ind_alpha) .and. &
      (BC_element /= ind_csurf)) then
    call MOM_error(FATAL,"extract_coupler_values: Unrecognized BC_element.")
  endif

  ! These error messages should be made more explicit.
!  if (.not.associated(BC_struc%bc(BC_index))) &
  if (.not.associated(BC_struc%bc)) &
    call MOM_error(FATAL,"set_coupler_values: " // &
       "The requested boundary condition is not associated.")
!  if (.not.associated(BC_struc%bc(BC_index)%field(BC_element))) &
  if (.not.associated(BC_struc%bc(BC_index)%field)) &
    call MOM_error(FATAL,"set_coupler_values: " // &
       "The requested boundary condition element is not associated.")
  if (.not.associated(BC_struc%bc(BC_index)%field(BC_element)%values)) &
    call MOM_error(FATAL,"set_coupler_values: " // &
       "The requested boundary condition value array is not associated.")

  Array_out => BC_struc%bc(BC_index)%field(BC_element)%values

  if (present(is)) then ; is0 = is ; else ; is0 = LBOUND(array_in,1) ; endif
  if (present(ie)) then ; ie0 = ie ; else ; ie0 = UBOUND(array_in,1) ; endif
  if (present(js)) then ; js0 = js ; else ; js0 = LBOUND(array_in,2) ; endif
  if (present(je)) then ; je0 = je ; else ; je0 = UBOUND(array_in,2) ; endif

  conv = 1.0 ; if (present(conversion)) conv = conversion

  if (size(Array_out,1) /= ie0 - is0 + 1) &
    call MOM_error(FATAL,"extract_coupler_values: Mismatch in i-size " // &
                   "between BC array and input array or computational domain.")
  if (size(Array_out,2) /= je0 - js0 + 1) &
    call MOM_error(FATAL,"extract_coupler_values: Mismatch in i-size " // &
                   "between BC array and input array or computational domain.")
  i_offset = lbound(Array_out,1) - is0
  j_offset = lbound(Array_out,2) - js0
  do j=js0,je0 ; do i=is0,ie0
    Array_out(i+i_offset,j+j_offset) = conv * array_in(i,j)
  enddo ; enddo

end subroutine set_coupler_values

end module coupler_util
