!> These interfaces allow for ocean or sea-ice variables to be replaced with data.
module MOM_data_override_infra

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domain_infra,  only : MOM_domain_type, domain2d
use MOM_domain_infra,  only : get_simple_array_i_ind, get_simple_array_j_ind
use MOM_time_manager,  only : time_type
use data_override_mod, only : data_override_init
use data_override_mod, only : data_override
use data_override_mod, only : data_override_unset_domains

implicit none ; private

public :: impose_data_init, impose_data, impose_data_unset_domains

!> Potentially override the values of a field in the model with values from a dataset.
interface impose_data
  module procedure data_override_MD, data_override_2d
end interface

contains

!> Initialize the data override capability and set the domains for the ocean and ice components.
!> There should be a call to impose_data_init before impose_data is called.
subroutine impose_data_init(MOM_domain_in, Ocean_domain_in, Ice_domain_in)
  type (MOM_domain_type), intent(in), optional :: MOM_domain_in
  type (domain2d),        intent(in), optional :: Ocean_domain_in
  type (domain2d),        intent(in), optional :: Ice_domain_in

  if (present(MOM_domain_in)) then
    call data_override_init(Ocean_domain_in=MOM_domain_in%mpp_domain, Ice_domain_in=Ice_domain_in)
  else
    call data_override_init(Ocean_domain_in=Ocean_domain_in, Ice_domain_in=Ice_domain_in)
  endif
end subroutine impose_data_init


!> Potentially override a 2-d field on a MOM6 domain with values from a dataset.
subroutine data_override_MD(domain, fieldname, data_2D, time, scale, override, is_ice)
  type(MOM_domain_type), intent(in)   :: domain   !< MOM domain from which to extract information
  character(len=*),     intent(in)    :: fieldname !< Name of the field to override
  real, dimension(:,:), intent(inout) :: data_2D  !< Data that may be modified by this call.
  type(time_type),      intent(in)    :: time     !< The model time, and the time for the data
  real,       optional, intent(in)    :: scale    !< A scaling factor that an overridden field is
                                                  !! multiplied by before it is returned.  However,
                                                  !! if there is no override, there is no rescaling.
  logical,    optional, intent(out)   :: override !< True if the field has been overridden successfully
  logical,    optional, intent(in)    :: is_ice   !< If present and true, use the ice domain.

  logical :: overridden, is_ocean
  integer :: i, j, is, ie, js, je

  overridden = .false.
  is_ocean = .true. ; if (present(is_ice)) is_ocean = .not.is_ice
  if (is_ocean) then
    call data_override('OCN', fieldname, data_2D, time, override=overridden)
  else
    call data_override('ICE', fieldname, data_2D, time, override=overridden)
  endif

  if (overridden .and. present(scale)) then ; if (scale /= 1.0) then
    ! Rescale data in the computational domain if the data override has occurred.
    call get_simple_array_i_ind(domain, size(data_2D,1), is, ie)
    call get_simple_array_j_ind(domain, size(data_2D,2), js, je)
    do j=js,je ; do i=is,ie
      data_2D(i,j) = scale*data_2D(i,j)
    enddo ; enddo
  endif ; endif

  if (present(override)) override = overridden

end subroutine data_override_MD


!> Potentially override a 2-d field with values from a dataset.
subroutine data_override_2d(gridname, fieldname, data_2D, time, override)
  character(len=3),     intent(in)    :: gridname !< String identifying the model component, in MOM6
                                                  !! and SIS this may be either 'OCN' or 'ICE'
  character(len=*),     intent(in)    :: fieldname !< Name of the field to override
  real, dimension(:,:), intent(inout) :: data_2D  !< Data that may be modified by this call
  type(time_type),      intent(in)    :: time     !< The model time, and the time for the data
  logical,    optional, intent(out)   :: override !< True if the field has been overridden successfully

  call data_override(gridname, fieldname, data_2D, time, override)

end subroutine data_override_2d

!> Unset domains that had previously been set for use by data_override.
subroutine impose_data_unset_domains(unset_Ocean, unset_Ice, must_be_set)
  logical, intent(in), optional :: unset_Ocean !< If present and true, unset the ocean domain for overrides
  logical, intent(in), optional :: unset_Ice   !< If present and true, unset the sea-ice domain for overrides
  logical, intent(in), optional :: must_be_set !< If present and true, it is a fatal error to unset
                                               !! a domain that is not set.

  call data_override_unset_domains(unset_Ocean=unset_Ocean, unset_Ice=unset_Ice, &
                                   must_be_set=must_be_set)
end subroutine impose_data_unset_domains

end module MOM_data_override_infra

!> \namespace MOM_data_override_infra
!!
!! The routines here wrap routines from the FMS module data_override_mod, which potentially replace
!! model values with values read from a data file.
