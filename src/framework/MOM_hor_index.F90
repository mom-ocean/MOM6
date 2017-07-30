!> Defines the horizontal index type (hor_index_type) used for providing index ranges
module MOM_hor_index

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domains, only : MOM_domain_type, get_domain_extent
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type

implicit none ; private

public :: hor_index_init, assignment(=)

!> Container for horizontal index ranges for data, computational and global domains
type, public :: hor_index_type
  integer :: isc !< The start i-index of cell centers within the computational domain
  integer :: iec !< The end i-index of cell centers within the computational domain
  integer :: jsc !< The start j-index of cell centers within the computational domain
  integer :: jec !< The end j-index of cell centers within the computational domain

  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain

  integer :: isg !< The start i-index of cell centers within the global domain
  integer :: ieg !< The end i-index of cell centers within the global domain
  integer :: jsg !< The start j-index of cell centers within the global domain
  integer :: jeg !< The end j-index of cell centers within the global domain

  integer :: IscB !< The start i-index of cell vertices within the computational domain
  integer :: IecB !< The end i-index of cell vertices within the computational domain
  integer :: JscB !< The start j-index of cell vertices within the computational domain
  integer :: JecB !< The end j-index of cell vertices within the computational domain

  integer :: IsdB !< The start i-index of cell vertices within the data domain
  integer :: IedB !< The end i-index of cell vertices within the data domain
  integer :: JsdB !< The start j-index of cell vertices within the data domain
  integer :: JedB !< The end j-index of cell vertices within the data domain

  integer :: IsgB !< The start i-index of cell vertices within the global domain
  integer :: IegB !< The end i-index of cell vertices within the global domain
  integer :: JsgB !< The start j-index of cell vertices within the global domain
  integer :: JegB !< The end j-index of cell vertices within the global domain

  integer :: idg_offset !< The offset between the corresponding global and local i-indices.
  integer :: jdg_offset !< The offset between the corresponding global and local j-indices.
  logical :: symmetric  !< True if symmetric memory is used.
end type hor_index_type

interface assignment(=); module procedure HIT_assign ; end interface

contains

!> Sets various index values in a hor_index_type.
subroutine hor_index_init(Domain, HI, param_file, local_indexing, index_offset)
  type(MOM_domain_type),  intent(in)    :: Domain     !< The MOM domain from which to extract information.
  type(hor_index_type),   intent(inout) :: HI         !< A horizontal index type to populate with data
  type(param_file_type),  intent(in)    :: param_file !< Parameter file handle
  logical, optional,      intent(in)    :: local_indexing !< If true, all tracer data domains start at 1
  integer, optional,      intent(in)    :: index_offset   !< A fixed additional offset to all indices

! This include declares and sets the variable "version".
#include "version_variable.h"

  ! get_domain_extent ensures that domains start at 1 for compatibility between
  ! static and dynamically allocated arrays.
  call get_domain_extent(Domain, HI%isc, HI%iec, HI%jsc, HI%jec, &
                         HI%isd, HI%ied, HI%jsd, HI%jed, &
                         HI%isg, HI%ieg, HI%jsg, HI%jeg, &
                         HI%idg_offset, HI%jdg_offset, HI%symmetric, &
                         local_indexing=local_indexing)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, "MOM_hor_index", version, &
                   "Sets the horizontal array index types.")

  HI%IscB = HI%isc ; HI%JscB = HI%jsc
  HI%IsdB = HI%isd ; HI%JsdB = HI%jsd
  HI%IsgB = HI%isg ; HI%JsgB = HI%jsg
  if (HI%symmetric) then
    HI%IscB = HI%isc-1 ; HI%JscB = HI%jsc-1
    HI%IsdB = HI%isd-1 ; HI%JsdB = HI%jsd-1
    HI%IsgB = HI%isg-1 ; HI%JsgB = HI%jsg-1
  endif
  HI%IecB = HI%iec ; HI%JecB = HI%jec
  HI%IedB = HI%ied ; HI%JedB = HI%jed
  HI%IegB = HI%ieg ; HI%JegB = HI%jeg

end subroutine hor_index_init

!> HIT_assign copies one hor_index_type into another.  It is accessed via an
!! assignment (=) operator.
subroutine HIT_assign(HI1, HI2)
  type(hor_index_type), intent(out) :: HI1 !< Horizontal index type to copy to
  type(hor_index_type), intent(in)  :: HI2 !< Horizontal index type to copy from
  ! This subroutine copies all components of the horizontal array index type
  ! variable on the RHS (HI2) to the variable on the LHS (HI1).

  HI1%isc = HI2%isc ; HI1%iec = HI2%iec ; HI1%jsc = HI2%jsc ; HI1%jec = HI2%jec
  HI1%isd = HI2%isd ; HI1%ied = HI2%ied ; HI1%jsd = HI2%jsd ; HI1%jed = HI2%jed
  HI1%isg = HI2%isg ; HI1%ieg = HI2%ieg ; HI1%jsg = HI2%jsg ; HI1%jeg = HI2%jeg

  HI1%IscB = HI2%IscB ; HI1%IecB = HI2%IecB ; HI1%JscB = HI2%JscB ; HI1%JecB = HI2%JecB
  HI1%IsdB = HI2%IsdB ; HI1%IedB = HI2%IedB ; HI1%JsdB = HI2%JsdB ; HI1%JedB = HI2%JedB
  HI1%IsgB = HI2%IsgB ; HI1%IegB = HI2%IegB ; HI1%JsgB = HI2%JsgB ; HI1%JegB = HI2%JegB

  HI1%idg_offset = HI2%idg_offset ; HI1%jdg_offset = HI2%jdg_offset
  HI1%symmetric = HI2%symmetric

end subroutine HIT_assign

!> \namespace mom_hor_index
!!
!! The hor_index_type provides the decalarations and loop ranges for almost all data with horizontal extent.
!!
!! Declarations and loop ranges should always be coded with the symmetric memory model in mind.
!! The non-symmetric memory mode will then also work, albeit with a different (less efficient) communication pattern.
!!
!! Using the hor_index_type HI:
!! - declaration of h-point data is of the form `h(HI%%isd:HI%%ied,HI%%jsd:HI%%jed)`;
!! - declaration of q-point data is of the form `q(HI%%IsdB:HI%%IedB,HI%%JsdB:HI%%JedB)`;
!! - declaration of u-point data is of the form `u(HI%%IsdB:HI%%IedB,HI%%jsd:HI%%jed)`;
!! - declaration of v-point data is of the form `v(HI%%isd:HI%%ied,HI%%JsdB:HI%%JedB)`.
!!
!! For more detail explanation of horizontal indexing see \ref Horizontal_indexing.


end module MOM_hor_index
