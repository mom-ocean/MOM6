module coupler_types_mod
!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*    This file is a part of MOM.  See MOM.F90 for licensing.          *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**
!
! This file was modified for use in MOM from MOM4 code by Stephen Griffies.
!
! This module contains simplified type declarations for the coupler
! of an ocean-only model. With mom4, this is used with ocean_solo_mod
! as the driver, and with MOM it is used with MOM_driver.
!
! This module contains simplified type declarations for the coupler
! of an ocean-only model. With mom4, this is used with ocean_solo_mod
! as the driver, and with MOM it is used with MOM_driver.

implicit none ; private

type, public :: coupler_2d_bc_type
  integer    :: num_bcs = 0
end type coupler_2d_bc_type
integer, public :: ind_flux=-1, ind_alpha=-2, ind_csurf=-3

end module coupler_types_mod


!============= DUMMY VERSION OF coupler_util ===================================

module coupler_util

use MOM_error_handler, only : MOM_error, FATAL, WARNING
use coupler_types_mod, only : coupler_2d_bc_type, ind_flux, ind_alpha, ind_csurf

implicit none ; private

public :: extract_coupler_values, set_coupler_values
public :: ind_flux, ind_alpha, ind_csurf

contains

subroutine extract_coupler_values(BC_struc, BC_index, BC_element, array_out, &
                                  is, ie, js, je, conversion)
  type(coupler_2d_bc_type), intent(in)  :: BC_struc
  integer,                  intent(in)  :: BC_index, BC_element
  real, dimension(:,:),     intent(inout) :: array_out
  integer,        optional, intent(in)  :: is, ie, js, je
  real,           optional, intent(in)  :: conversion
! Arguments: BC_struc - The type from which the data is being extracted.
!  (in)      BC_index - The boundary condition number being extracted.
!  (in)      BC_element - The element of the boundary condition being extracted.
!            This could be ind_csurf, ind_alpha, ind_flux or ind_deposition.
!  (out)     array_out - The array being filled with the input values.
!  (in, opt) is, ie, js, je - The i- and j- limits of array_out to be filled.
!            These must match the size of the corresponding value array or an
!            error message is issued.
!  (in, opt) conversion - A number that every element is multiplied by, to
!                         permit sign convention or unit conversion.

  call MOM_error(FATAL,"The called version of extract_coupler_values does" // &
                 " nothing but issue this fatal error message.")

end subroutine extract_coupler_values

subroutine set_coupler_values(array_in, BC_struc, BC_index, BC_element, &
                              is, ie, js, je, conversion)
  real, dimension(:,:),     intent(in)    :: array_in
  type(coupler_2d_bc_type), intent(inout) :: BC_struc
  integer,                  intent(in)    :: BC_index, BC_element
  integer,        optional, intent(in)    :: is, ie, js, je
  real,           optional, intent(in)    :: conversion
! Arguments: array_in - The array containing the values to load into the BC.
!  (out)     BC_struc - The type into which the data is being loaded.
!  (in)      BC_index - The boundary condition number being extracted.
!  (in)      BC_element - The element of the boundary condition being extracted.
!            This could be ind_csurf, ind_alpha, ind_flux or ind_deposition.
!  (in, opt) is, ie, js, je - The i- and j- limits of array_out to be filled.
!            These must match the size of the corresponding value array or an
!            error message is issued.
!  (in, opt) conversion - A number that every element is multiplied by, to
!                         permit sign convention or unit conversion.

  call MOM_error(FATAL,"The called version of set_coupler_values does" // &
                 " nothing but issue this fatal error message.")

end subroutine set_coupler_values

end module coupler_util
