!> This module provides coupler type interfaces for use by MOM6
module MOM_coupler_types

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_couplertype_infra, only : CT_spawn, CT_initialized, CT_destructor, atmos_ocn_coupler_flux
use MOM_couplertype_infra, only : CT_set_diags, CT_send_data, CT_write_chksums, CT_data_override
use MOM_couplertype_infra, only : CT_copy_data, CT_increment_data, CT_rescale_data
use MOM_couplertype_infra, only : CT_set_data, CT_extract_data, CT_redistribute_data
use MOM_couplertype_infra, only : coupler_1d_bc_type, coupler_2d_bc_type, coupler_3d_bc_type
use MOM_couplertype_infra, only : ind_flux, ind_alpha, ind_csurf
use MOM_domain_infra,      only : domain2D
use MOM_time_manager,      only : time_type

implicit none ; private

public :: coupler_type_spawn, coupler_type_destructor, coupler_type_initialized
public :: coupler_type_set_diags, coupler_type_send_data, coupler_type_write_chksums
public :: set_coupler_type_data, extract_coupler_type_data, coupler_type_redistribute_data
public :: coupler_type_copy_data, coupler_type_increment_data, coupler_type_rescale_data
public :: atmos_ocn_coupler_flux, coupler_type_data_override
public :: ind_flux, ind_alpha, ind_csurf
public :: coupler_1d_bc_type, coupler_2d_bc_type, coupler_3d_bc_type

!> This is the interface to spawn one coupler_bc_type into another.
interface coupler_type_spawn
  module procedure CT_spawn_1d_2d, CT_spawn_1d_3d, CT_spawn_2d_2d, CT_spawn_2d_3d
  module procedure CT_spawn_3d_2d, CT_spawn_3d_3d
end interface coupler_type_spawn

!> This function interface indicates whether a coupler_bc_type has been initialized.
interface coupler_type_initialized
  module procedure CT_initialized_1d, CT_initialized_2d, CT_initialized_3d
end interface coupler_type_initialized

!> This is the interface to deallocate any data associated with a coupler_bc_type.
interface coupler_type_destructor
  module procedure CT_destructor_1d, CT_destructor_2d
end interface coupler_type_destructor

!> Copy all elements of the data in either a coupler_2d_bc_type or a coupler_3d_bc_type into
!! another structure of the same or the other type.  Both must have the same array sizes in common
!! dimensions, while the details of any expansion from 2d to 3d are controlled by arguments.
interface coupler_type_copy_data
  module procedure CT_copy_data_2d, CT_copy_data_3d, CT_copy_data_2d_3d
end interface coupler_type_copy_data

!> Increment data in all elements of one coupler_2d_bc_type or coupler_3d_bc_type with
!! the data from another. Both must have the same array sizes and rank.
interface coupler_type_increment_data
  module procedure CT_increment_data_2d, CT_increment_data_3d, CT_increment_data_2d_3d
end interface coupler_type_increment_data

!> Rescale the fields in the elements of a coupler_2d_bc_type by multiplying by a factor scale.
!! If scale is 0, this is a direct assignment to 0, so that NaNs will not persist.
interface coupler_type_rescale_data
  module procedure CT_rescale_data_2d, CT_rescale_data_3d
end interface coupler_type_rescale_data

!> Redistribute the data in all elements of one coupler_2d_bc_type or coupler_3d_bc_type into
!! another, which may be on different processors with a different decomposition.
interface coupler_type_redistribute_data
  module procedure CT_redistribute_data_2d, CT_redistribute_data_3d
end interface coupler_type_redistribute_data

!> Write out checksums for the elements of a coupler_2d_bc_type or coupler_3d_bc_type
interface coupler_type_write_chksums
  module procedure CT_write_chksums_2d, CT_write_chksums_3d
end interface coupler_type_write_chksums

contains

!> Generate a 2-D coupler type using a 1-D coupler type as a template.
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

  call CT_spawn(var_in, var, idim, jdim, suffix=suffix, as_needed=as_needed)

end subroutine  CT_spawn_1d_2d

!> Generate a 3-D coupler type using a 1-D coupler type as a template.
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

  call CT_spawn(var_in, var, idim, jdim, kdim, suffix=suffix, as_needed=as_needed)

end subroutine  CT_spawn_1d_3d

!> Generate one 2-D coupler type using another 2-D coupler type as a template.
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

  call CT_spawn(var_in, var, idim, jdim, suffix=suffix, as_needed=as_needed)

end subroutine  CT_spawn_2d_2d

!> Generate a 3-D coupler type using a 2-D coupler type as a template.
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

  call CT_spawn(var_in, var, idim, jdim, kdim, suffix=suffix, as_needed=as_needed)

end subroutine  CT_spawn_2d_3d

!> Generate a 2-D coupler type using a 3-D coupler type as a template.
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

  call CT_spawn(var_in, var, idim, jdim, suffix=suffix, as_needed=as_needed)

end subroutine  CT_spawn_3d_2d

!> Generate a 3-D coupler type using another 3-D coupler type as a template.
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

  call CT_spawn(var_in, var, idim, jdim, kdim, suffix=suffix, as_needed=as_needed)

end subroutine  CT_spawn_3d_3d

!> Copy all elements of the data in a coupler_2d_bc_type into another.  Both must have the same array sizes.
subroutine CT_copy_data_2d(var_in, var, halo_size, bc_index, field_index, &
                  exclude_flux_type, only_flux_type, pass_through_ice)
  type(coupler_2d_bc_type),   intent(in)    :: var_in  !< BC_type structure with the data to copy
  type(coupler_2d_bc_type),   intent(inout) :: var     !< The recipient BC_type structure
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                         !! that is being copied
  integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                         !! boundary condition that is being copied
  character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types of fluxes
                                                         !! to exclude from this copy.
  character(len=*), optional, intent(in)    :: only_flux_type    !< A string describing which types of fluxes
                                                         !! to include from this copy.
  logical,          optional, intent(in)    :: pass_through_ice !< If true, only copy BCs whose
                                                         !! value of pass_through ice matches this

  call CT_copy_data(var_in, var, halo_size, bc_index, field_index, &
                    exclude_flux_type, only_flux_type, pass_through_ice)
end subroutine CT_copy_data_2d

!> Copy all elements of the data in a coupler_3d_bc_type into another.  Both must have the same array sizes.
subroutine CT_copy_data_3d(var_in, var, halo_size, bc_index, field_index, &
                  exclude_flux_type, only_flux_type, pass_through_ice)
  type(coupler_3d_bc_type),   intent(in)    :: var_in  !< BC_type structure with the data to copy
  type(coupler_3d_bc_type),   intent(inout) :: var     !< The recipient BC_type structure
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                         !! that is being copied
  integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                         !! boundary condition that is being copied
  character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types of fluxes
                                                         !! to exclude from this copy.
  character(len=*), optional, intent(in)    :: only_flux_type    !< A string describing which types of fluxes
                                                         !! to include from this copy.
  logical,          optional, intent(in)    :: pass_through_ice !< If true, only copy BCs whose
                                                         !! value of pass_through ice matches this

  call CT_copy_data(var_in, var, halo_size, bc_index, field_index, &
                    exclude_flux_type, only_flux_type, pass_through_ice)
end subroutine CT_copy_data_3d

!> Copy all elements of the data in a coupler_2d_bc_type into a coupler_3d_bc_type.
!! Both must have the same array sizes for the first two dimensions, while the extent
!! of the 3rd dimension that is being filled may be specified via optional arguments.
subroutine CT_copy_data_2d_3d(var_in, var, halo_size, bc_index, field_index, &
                  exclude_flux_type, only_flux_type, pass_through_ice, ind3_start, ind3_end)
  type(coupler_2d_bc_type),   intent(in)    :: var_in  !< BC_type structure with the data to copy
  type(coupler_3d_bc_type),   intent(inout) :: var     !< The recipient BC_type structure
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer,          optional, intent(in)    :: bc_index  !< The index of the boundary condition
                                                         !! that is being copied
  integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                         !! boundary condition that is being copied
  character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types of fluxes
                                                         !! to exclude from this copy.
  character(len=*), optional, intent(in)    :: only_flux_type    !< A string describing which types of fluxes
                                                         !! to include from this copy.
  logical,          optional, intent(in)    :: pass_through_ice !< If true, only copy BCs whose
                                                         !! value of pass_through ice matches this
  integer,          optional, intent(in)    :: ind3_start !< The starting value of the 3rd
                                                         !! index of the 3d type to fill in.
  integer,          optional, intent(in)    :: ind3_end  !< The ending value of the 3rd
                                                         !! index of the 3d type to fill in.

  call CT_copy_data(var_in, var, halo_size, bc_index, field_index, &
                    exclude_flux_type, only_flux_type, pass_through_ice, ind3_start, ind3_end)
end subroutine CT_copy_data_2d_3d

!> Increment data in all elements of one coupler_2d_bc_type with the data from another. Both
!! must have the same array sizes.
subroutine CT_increment_data_2d(var_in, var, halo_size, scale_factor, scale_prev)
  type(coupler_2d_bc_type),   intent(in)    :: var_in  !< coupler_type structure with the data to add to the other type
  type(coupler_2d_bc_type),   intent(inout) :: var     !< The coupler_type structure whose fields are being incremented
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to increment; 0 by default
  real,             optional, intent(in)    :: scale_factor  !< A scaling factor for the data that is being added
  real,             optional, intent(in)    :: scale_prev    !< A scaling factor for the data that is already here

  call CT_increment_data(var_in, var, halo_size=halo_size, scale_factor=scale_factor, &
                         scale_prev=scale_prev)

end subroutine CT_increment_data_2d

!> Increment data in all elements of one coupler_3d_bc_type with the data from another. Both
!! must have the same array sizes.
subroutine CT_increment_data_3d(var_in, var, halo_size, scale_factor, scale_prev, exclude_flux_type, only_flux_type)
  type(coupler_3d_bc_type),   intent(in)    :: var_in  !< coupler_type structure with the data to add to the other type
  type(coupler_3d_bc_type),   intent(inout) :: var     !< The coupler_type structure whose fields are being incremented
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to increment; 0 by default
  real,             optional, intent(in)    :: scale_factor  !< A scaling factor for the data that is being added
  real,             optional, intent(in)    :: scale_prev    !< A scaling factor for the data that is already here
  character(len=*), optional, intent(in)    :: exclude_flux_type !< A string describing which types
                                                         !! of fluxes to exclude from this increment.
  character(len=*), optional, intent(in)    :: only_flux_type    !< A string describing which types
                                                         !! of fluxes to include from this increment.

  call CT_increment_data(var_in, var, halo_size=halo_size, scale_factor=scale_factor, &
                         scale_prev=scale_prev, exclude_flux_type=exclude_flux_type, &
                         only_flux_type=only_flux_type)

end subroutine CT_increment_data_3d

!> Increment data in the elements of a coupler_2d_bc_type with weighted averages of elements of a
!! coupler_3d_bc_type
subroutine CT_increment_data_2d_3d(var_in, weights, var, halo_size)
  type(coupler_3d_bc_type),   intent(in)    :: var_in  !< coupler_type structure with the data to add to the other type
  real, dimension(:,:,:),     intent(in)    :: weights !< An array of normalized weights for the 3d-data to
                                                       !! increment the 2d-data.  There is no renormalization,
                                                       !! so if the weights do not sum to 1 in the 3rd dimension
                                                       !! there may be adverse consequences!
  type(coupler_2d_bc_type),   intent(inout) :: var     !< The coupler_type structure whose fields are being incremented
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to increment; 0 by default

  call CT_increment_data(var_in, weights, var, halo_size=halo_size)

end subroutine CT_increment_data_2d_3d

!> Rescales the fields in the elements of a coupler_2d_bc_type by multiplying by a factor scale.
!! If scale is 0, this is a direct assignment to 0, so that NaNs will not persist.
subroutine CT_rescale_data_2d(var, scale)
  type(coupler_2d_bc_type),   intent(inout) :: var   !< The BC_type structure whose fields are being rescaled
  real,                       intent(in)    :: scale !< A scaling factor to multiply fields by

  call CT_rescale_data(var, scale)

end subroutine CT_rescale_data_2d

!> Rescales the fields in the elements of a coupler_3d_bc_type by multiplying by a factor scale.
!! If scale is 0, this is a direct assignment to 0, so that NaNs will not persist.
subroutine CT_rescale_data_3d(var, scale)
  type(coupler_3d_bc_type),   intent(inout) :: var   !< The BC_type structure whose fields are being rescaled
  real,                       intent(in)    :: scale !< A scaling factor to multiply fields by

  call CT_rescale_data(var, scale)

end subroutine CT_rescale_data_3d

!> Redistribute the data in all elements of one coupler_2d_bc_type into another, which may be on
!! different processors with a different decomposition.
subroutine CT_redistribute_data_2d(var_in, domain_in, var_out, domain_out, complete)
  type(coupler_2d_bc_type), intent(in)    :: var_in     !< BC_type structure with the data to copy
  type(domain2D),           intent(in)    :: domain_in  !< The FMS domain for the input structure
  type(coupler_2d_bc_type), intent(inout) :: var_out    !< The recipient BC_type structure (data intent out)
  type(domain2D),           intent(in)    :: domain_out !< The FMS domain for the output structure
  logical,        optional, intent(in)    :: complete   !< If true, complete the updates

  call CT_redistribute_data(var_in, domain_in, var_out, domain_out, complete)
end subroutine CT_redistribute_data_2d

!> Redistribute the data in all elements of one coupler_3d_bc_type into another, which may be on
!! different processors with a different decomposition.
subroutine CT_redistribute_data_3d(var_in, domain_in, var_out, domain_out, complete)
  type(coupler_3d_bc_type), intent(in)    :: var_in     !< BC_type structure with the data to copy
  type(domain2D),           intent(in)    :: domain_in  !< The FMS domain for the input structure
  type(coupler_3d_bc_type), intent(inout) :: var_out    !< The recipient BC_type structure (data intent out)
  type(domain2D),           intent(in)    :: domain_out !< The FMS domain for the output structure
  logical,        optional, intent(in)    :: complete   !< If true, complete the updates

  call CT_redistribute_data(var_in, domain_in, var_out, domain_out, complete)
end subroutine CT_redistribute_data_3d


!> Potentially override the values in a coupler_2d_bc_type
subroutine coupler_type_data_override(gridname, var, time)
  character(len=3),         intent(in)    :: gridname !< 3-character long model grid ID
  type(coupler_2d_bc_type), intent(inout) :: var  !< BC_type structure to override
  type(time_type),          intent(in)    :: time !< The current model time

  call CT_data_override(gridname, var, time)
end subroutine coupler_type_data_override


!> Extract a 2d field from a coupler_2d_bc_type into a two-dimensional array, using a
!! MOM-specific interface.
subroutine extract_coupler_type_data(var_in, bc_index, array_out, scale_factor, &
                                     halo_size, idim, jdim, field_index)
  type(coupler_2d_bc_type),   intent(in)    :: var_in    !< BC_type structure with the data to extract
  integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                         !! that is being copied
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
  integer,          optional, intent(in)    :: field_index !< The index of the field in the boundary
                                                         !! condition that is being copied, or the
                                                         !! surface flux by default.

  if (present(field_index)) then
    call CT_extract_data(var_in, bc_index, field_index, array_out, &
                    scale_factor=scale_factor, halo_size=halo_size, idim=idim, jdim=jdim)
  else
    call CT_extract_data(var_in, bc_index, ind_flux, array_out, &
                    scale_factor=scale_factor, halo_size=halo_size, idim=idim, jdim=jdim)
  endif

end subroutine extract_coupler_type_data

!> Set single 2d field in coupler_2d_bc_type from a two-dimensional array, using a
!! MOM-specific interface.
subroutine set_coupler_type_data(array_in, bc_index, var, solubility, scale_factor, &
                                 halo_size, idim, jdim, field_index)
  real, dimension(1:,1:),     intent(in)   :: array_in   !< The source array for the field; its size
                                                         !! must match the size of the data being copied
                                                         !! unless idim and jdim are supplied.
  integer,                    intent(in)    :: bc_index  !< The index of the boundary condition
                                                         !! that is being copied
  type(coupler_2d_bc_type),   intent(inout) :: var       !< BC_type structure with the data to set
  logical,          optional, intent(in)    :: solubility !< If true and field index is missing, set
                                                         !! the solubility field.  Otherwise set the
                                                         !! surface concentration (the default).
  real,             optional, intent(in)    :: scale_factor !< A scaling factor for the data that is being added
  integer,          optional, intent(in)    :: halo_size !< The extent of the halo to copy; 0 by default
  integer, dimension(4), optional, intent(in) :: idim    !< The data and computational domain extents of
                                                         !! the first dimension of the output array
                                                         !! in a non-decreasing list
  integer, dimension(4), optional, intent(in) :: jdim    !< The data and computational domain extents of
                                                         !! the second dimension of the output array
                                                         !! in a non-decreasing list
  integer,          optional, intent(in)    :: field_index !< The index of the field in the
                                                         !! boundary condition that is being set. The
                                                         !! surface concentration is set by default.

  integer :: subfield ! An integer indicating which field to set.

  subfield = ind_csurf
  if (present(solubility)) then ; if (solubility) subfield = ind_alpha ; endif
  if (present(field_index)) subfield = field_index

  call CT_set_data(array_in, bc_index, subfield, var, &
                   scale_factor=scale_factor, halo_size=halo_size, idim=idim, jdim=jdim)

end subroutine set_coupler_type_data

!> Register the diagnostics of a coupler_2d_bc_type
subroutine coupler_type_set_diags(var, diag_name, axes, time)
  type(coupler_2d_bc_type), intent(inout) :: var  !< BC_type structure for which to register diagnostics
  character(len=*),         intent(in)    :: diag_name !< name for diagnostic file, or blank not to register the fields
  integer, dimension(:),    intent(in)    :: axes !< array of axes identifiers for diagnostic variable registration
  type(time_type),          intent(in)    :: time !< model time variable for registering diagnostic field

  call CT_set_diags(var, diag_name, axes, time)

end subroutine coupler_type_set_diags

!> Write out all diagnostics of elements of a coupler_2d_bc_type
subroutine coupler_type_send_data(var, Time)
  type(coupler_2d_bc_type), intent(in) :: var  !< BC_type structure with the diagnostics to write
  type(time_type),          intent(in) :: time !< The current model time

  call CT_send_data(var, Time)
end subroutine coupler_type_send_data

!> Write out checksums for the elements of a coupler_2d_bc_type
subroutine CT_write_chksums_2d(var, outunit, name_lead)
  type(coupler_2d_bc_type),   intent(in) :: var  !< BC_type structure for which to register diagnostics
  integer,                    intent(in) :: outunit !< The index of a open output file
  character(len=*), optional, intent(in) :: name_lead !< An optional prefix for the variable names

  call CT_write_chksums(var, outunit, name_lead)

end subroutine CT_write_chksums_2d

!> Write out checksums for the elements of a coupler_3d_bc_type
subroutine CT_write_chksums_3d(var, outunit, name_lead)
  type(coupler_3d_bc_type),   intent(in) :: var  !< BC_type structure for which to register diagnostics
  integer,                    intent(in) :: outunit !< The index of a open output file
  character(len=*), optional, intent(in) :: name_lead !< An optional prefix for the variable names

  call CT_write_chksums(var, outunit, name_lead)

end subroutine CT_write_chksums_3d

!> Indicate whether a coupler_1d_bc_type has been initialized.
logical function CT_initialized_1d(var)
  type(coupler_1d_bc_type), intent(in) :: var  !< BC_type structure to be deconstructed

  CT_initialized_1d = CT_initialized(var)
end function CT_initialized_1d

!> Indicate whether a coupler_2d_bc_type has been initialized.
logical function CT_initialized_2d(var)
  type(coupler_2d_bc_type), intent(in) :: var  !< BC_type structure to be deconstructed

  CT_initialized_2d = CT_initialized(var)
end function CT_initialized_2d

!> Indicate whether a coupler_3d_bc_type has been initialized.
logical function CT_initialized_3d(var)
  type(coupler_3d_bc_type), intent(in) :: var  !< BC_type structure to be deconstructed

  CT_initialized_3d = CT_initialized(var)
end function CT_initialized_3d

!> Deallocate all data associated with a coupler_1d_bc_type
subroutine CT_destructor_1d(var)
  type(coupler_1d_bc_type), intent(inout) :: var  !< BC_type structure to be deconstructed

  call CT_destructor(var)

end subroutine CT_destructor_1d

!> Deallocate all data associated with a coupler_2d_bc_type
subroutine CT_destructor_2d(var)
  type(coupler_2d_bc_type), intent(inout) :: var  !< BC_type structure to be deconstructed

  call CT_destructor(var)

end subroutine CT_destructor_2d

end module MOM_coupler_types
