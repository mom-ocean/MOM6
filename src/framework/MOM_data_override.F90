!> These interfaces allow for ocean or sea-ice variables to be replaced with data.
module MOM_data_override

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_data_override_infra, only : data_override_init => impose_data_init
use MOM_data_override_infra, only : data_override => impose_data
use MOM_data_override_infra, only : data_override_unset_domains => impose_data_unset_domains

implicit none ; private

!> Public functions:
!> mom_data_override_infra:impose_data_init
public :: data_override_init
!> mom_data_override_infra:impose_data
public :: data_override
!> mom_data_override_infra:impose_data_unset_domains
public :: data_override_unset_domains

end module MOM_data_override

!> \namespace MOM_data_override
!!
!! APIs are defined and implemented in MOM_data_override_infra
