!> A simple (very thin) wrapper for managing ensemble member layout information
module MOM_ensemble_manager_infra

! This file is part of MOM6. See LICENSE.md for the license.

use ensemble_manager_mod, only : FMS_ensemble_manager_init => ensemble_manager_init
use ensemble_manager_mod, only : FMS_ensemble_pelist_setup => ensemble_pelist_setup
use ensemble_manager_mod, only : FMS_get_ensemble_id => get_ensemble_id
use ensemble_manager_mod, only : FMS_get_ensemble_size => get_ensemble_size
use ensemble_manager_mod, only : FMS_get_ensemble_pelist => get_ensemble_pelist
use ensemble_manager_mod, only : FMS_get_ensemble_filter_pelist => get_ensemble_filter_pelist

implicit none ; private

public :: ensemble_manager_init, ensemble_pelist_setup
public :: get_ensemble_id, get_ensemble_size
public :: get_ensemble_pelist, get_ensemble_filter_pelist

contains

!> Initializes the ensemble manager which divides available resources
!! in order to concurrently execute an ensemble of model realizations.
subroutine ensemble_manager_init()

  call FMS_ensemble_manager_init()

end subroutine ensemble_manager_init

!> Create a list of processing elements (PEs) across components
!! associated with the current ensemble member.
subroutine ensemble_pelist_setup(concurrent, atmos_npes, ocean_npes, land_npes, ice_npes, &
                                   Atm_pelist, Ocean_pelist, Land_pelist, Ice_pelist)
  logical,               intent(in)    :: concurrent !< A logical flag, if True, then ocean fast
                                                     !! PEs are run concurrently with
                                                     !! slow PEs within the coupler.
  integer,               intent(in)    :: atmos_npes !< The number of atmospheric (fast) PEs
  integer,               intent(in)    :: ocean_npes !< The number of ocean (slow) PEs
  integer,               intent(in)    :: land_npes  !< The number of land PEs (fast)
  integer,               intent(in)    :: ice_npes   !< The number of ice (fast) PEs
  integer, dimension(:), intent(inout) :: Atm_pelist !< A list of Atm PEs
  integer, dimension(:), intent(inout) :: Ocean_pelist !< A list of Ocean PEs
  integer, dimension(:), intent(inout) :: Land_pelist !< A list of Land PEs
  integer, dimension(:), intent(inout) :: Ice_pelist !< A list of Ice PEs


  call FMS_ensemble_pelist_setup(concurrent, atmos_npes, ocean_npes, land_npes, ice_npes, &
         Atm_pelist, Ocean_pelist, Land_pelist, Ice_pelist)

end subroutine ensemble_pelist_setup

!> Returns the numeric id for the current ensemble member
function get_ensemble_id()
  integer :: get_ensemble_id

  get_ensemble_id = FMS_get_ensemble_id()

end function get_ensemble_id

!> Returns ensemble information as follows,
!! index (1) :: ensemble size
!! index (2) :: Number of PEs per ensemble member
!! index (3) :: Number of ocean PEs per ensemble member
!! index (4) :: Number of atmos PEs per ensemble member
!! index (5) :: Number of land PEs per ensemble member
!! index (6) :: Number of ice PEs per ensemble member
function get_ensemble_size()
  integer, dimension(6) :: get_ensemble_size

  get_ensemble_size = FMS_get_ensemble_size()

end function get_ensemble_size

!> Returns the list of PEs associated with all ensemble members
!! Results are stored in the argument array which must be large
!! enough to contain the list.  If the optional name argument is present,
!! the returned processor list are for a particular component (atmos, ocean ,land, ice)
subroutine get_ensemble_pelist(pelist, name)
  integer,                    intent(inout) :: pelist(:,:) !< A processor list for all ensemble members
  character(len=*), optional, intent(in)    :: name !< An optional component name (atmos, ocean, land, ice)

  call FMS_get_ensemble_pelist(pelist, name)

end subroutine get_ensemble_pelist

!> Returns the list of PEs associated with the named ensemble filter application.
!! Valid component names include ('atmos', 'ocean', 'land', and 'ice')
subroutine get_ensemble_filter_pelist(pelist, name)
  integer,          intent(inout) :: pelist(:) !< A processor list for the ensemble filter
  character(len=*), intent(in)    :: name      !< The component name (atmos, ocean, land, ice)

  call FMS_get_Ensemble_filter_pelist(pelist, name)

end subroutine get_ensemble_filter_pelist

end module MOM_ensemble_manager_infra
