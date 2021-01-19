!> This module wraps the FMS temporal and spatial interpolation routines
module MOM_interp_infra

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domain_infra,    only : MOM_domain_type, domain2d
use MOM_io_infra,        only : axistype
use MOM_time_manager,    only : time_type
use horiz_interp_mod,    only : horiz_interp_new, horiz_interp, horiz_interp_init, horiz_interp_type
use time_interp_external_mod, only : time_interp_external
use time_interp_external_mod, only : init_external_field, time_interp_external_init
use time_interp_external_mod, only : get_external_field_size
use time_interp_external_mod, only : get_external_field_axes, get_external_field_missing

implicit none ; private

public :: time_interp_extern, init_extern_field, time_interp_external_init
public :: get_external_field_info
public :: horiz_interp_type, horiz_interp_init, horiz_interp, horiz_interp_new

!> Read a field based on model time, and rotate to the model domain.
interface time_interp_extern
  module procedure time_interp_extern_0d
  module procedure time_interp_extern_2d
  module procedure time_interp_extern_3d
end interface time_interp_extern

contains

!> Get information about the external fields.
subroutine get_external_field_info(field_id, size, axes, missing)
  integer,                                intent(in)    :: field_id !< The integer index of the external
                                                                    !! field returned from a previous
                                                                    !! call to init_external_field()
  integer,        dimension(4), optional, intent(inout) :: size    !< Dimension sizes for the input data
  type(axistype), dimension(4), optional, intent(inout) :: axes    !< Axis types for the input data
  real,                         optional, intent(inout) :: missing !< Missing value for the input data

  if (present(size)) then
    size(1:4) = get_external_field_size(field_id)
  endif

  if (present(axes)) then
    axes(1:4) = get_external_field_axes(field_id)
  endif

  if (present(missing)) then
    missing = get_external_field_missing(field_id)
  endif

end subroutine get_external_field_info


!> Read a scalar field based on model time.
subroutine time_interp_extern_0d(field_id, time, data_in, verbose)
  integer,           intent(in)    :: field_id !< The integer index of the external field returned
                                               !! from a previous call to init_external_field()
  type(time_type),   intent(in)    :: time     !< The target time for the data
  real,              intent(inout) :: data_in  !< The interpolated value
  logical, optional, intent(in)    :: verbose  !< If true, write verbose output for debugging

  call time_interp_external(field_id, time, data_in, verbose=verbose)
end subroutine time_interp_extern_0d

!> Read a 2d field from an external based on model time, potentially including horizontal
!! interpolation and rotation of the data
subroutine time_interp_extern_2d(field_id, time, data_in, interp, verbose, horz_interp, mask_out)
  integer,              intent(in)    :: field_id !< The integer index of the external field returned
                                                  !! from a previous call to init_external_field()
  type(time_type),      intent(in)    :: time     !< The target time for the data
  real, dimension(:,:), intent(inout) :: data_in  !< The array in which to store the interpolated values
  integer,    optional, intent(in)    :: interp   !< A flag indicating the temporal interpolation method
  logical,    optional, intent(in)    :: verbose  !< If true, write verbose output for debugging
  type(horiz_interp_type), &
              optional, intent(in)    :: horz_interp !< A structure to control horizontal interpolation
  logical, dimension(:,:), &
              optional, intent(out)   :: mask_out !< An array that is true where there is valid data

  call time_interp_external(field_id, time, data_in, interp=interp, verbose=verbose, &
                            horz_interp=horz_interp, mask_out=mask_out)
end subroutine time_interp_extern_2d


!> Read a 3d field based on model time, and rotate to the model grid
subroutine time_interp_extern_3d(field_id, time, data_in, interp, verbose, horz_interp, mask_out)
  integer,                intent(in)    :: field_id !< The integer index of the external field returned
                                                    !! from a previous call to init_external_field()
  type(time_type),        intent(in)    :: time     !< The target time for the data
  real, dimension(:,:,:), intent(inout) :: data_in  !< The array in which to store the interpolated values
  integer,      optional, intent(in)    :: interp   !< A flag indicating the temporal interpolation method
  logical,      optional, intent(in)    :: verbose  !< If true, write verbose output for debugging
  type(horiz_interp_type), &
                optional, intent(in)    :: horz_interp !< A structure to control horizontal interpolation
  logical, dimension(:,:,:), &
                optional, intent(out)   :: mask_out !< An array that is true where there is valid data

  call time_interp_external(field_id, time, data_in, interp=interp, verbose=verbose, &
                            horz_interp=horz_interp, mask_out=mask_out)
end subroutine time_interp_extern_3d

integer function init_extern_field(file, fieldname, MOM_domain, domain, verbose, &
                                   threading, ierr, ignore_axis_atts )

  character(len=*),         intent(in)  :: file  !< The name of the file to read
  character(len=*),         intent(in)  :: fieldname !< The name of the field in the file
  integer,        optional, intent(in)  :: threading !< A flag specifying whether the root PE reads
                                                 !! the data and broadcasts it (SINGLE_FILE) or all
                                                 !! processors read (MULTIPLE, the default).
  logical,        optional, intent(in)  :: verbose !< If true, write verbose output for debugging
  type(domain2d), optional, intent(in)  :: domain !< A domain2d type that describes the decomposition
  type(MOM_domain_type), &
                  optional, intent(in)  :: MOM_Domain !< A MOM_Domain that describes the decomposition
  integer,        optional, intent(out) :: ierr  !< Returns a non-zero error code in case of failure
  logical,        optional, intent(in)  :: ignore_axis_atts !< If present and true, do not issue a
                                                 !! fatal error if the axis Cartesian attribute is
                                                 !! not set to a recognized value.

  if (present(MOM_Domain)) then
    init_extern_field = init_external_field(file, fieldname, domain=MOM_domain%mpp_domain, &
             verbose=verbose, threading=threading, ierr=ierr, ignore_axis_atts=ignore_axis_atts)
  else
    init_extern_field = init_external_field(file, fieldname, domain=domain, &
             verbose=verbose, threading=threading, ierr=ierr, ignore_axis_atts=ignore_axis_atts)
  endif

end function init_extern_field

end module MOM_interp_infra
