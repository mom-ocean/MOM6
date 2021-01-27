!> Manages ensemble member layout information
module MOM_ensemble_manager

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_ensemble_manager_infra, only : ensemble_manager_init
use MOM_ensemble_manager_infra, only : ensemble_pelist_setup
use MOM_ensemble_manager_infra, only : get_ensemble_id
use MOM_ensemble_manager_infra, only : get_ensemble_size
use MOM_ensemble_manager_infra, only : get_ensemble_pelist
use MOM_ensemble_manager_infra, only : get_ensemble_filter_pelist

implicit none ; private

!> Public functions:
!> mom_ensemble_manager_infra:ensemble_manager_init
public :: ensemble_manager_init
!> mom_ensemble_manager_infra:ensemble_pelist_setup
public :: ensemble_pelist_setup
!> mom_ensemble_manager_infra:get_ensemble_id
public :: get_ensemble_id
!> mom_ensemble_manager_infra:get_ensemble_size
public :: get_ensemble_size
!> mom_ensemble_manager_infra:get_ensemble_pelist
public :: get_ensemble_pelist
!> mom_ensemble_manager_infra:get_ensemble_filter_pelist
public :: get_ensemble_filter_pelist




end module MOM_ensemble_manager

!> \namespace mom_ensemble_manager
!!
!! APIs are defined and implemented in MOM_ensemble_manager_infra
