!> A dummy version of atmos_ocean_fluxes_mod module for
!! use when the vastly larger FMS package is not needed.
module atmos_ocean_fluxes_mod

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public :: aof_set_coupler_flux

contains

!> This subroutine duplicates an interface used by the FMS coupler, but only
!! returns a value of -1.  None of the arguments are used for anything.
function aof_set_coupler_flux(name, flux_type, implementation, atm_tr_index,     &
                              param, flag, ice_restart_file, ocean_restart_file, &
                              units, caller, verbosity)  result (coupler_index)

  character(len=*),                intent(in) :: name !< An unused argument
  character(len=*),                intent(in) :: flux_type !< An unused argument
  character(len=*),                intent(in) :: implementation !< An unused argument
  integer,               optional, intent(in) :: atm_tr_index !< An unused argument
  real,    dimension(:), optional, intent(in) :: param !< An unused argument
  logical, dimension(:), optional, intent(in) :: flag !< An unused argument
  character(len=*),      optional, intent(in) :: ice_restart_file !< An unused argument
  character(len=*),      optional, intent(in) :: ocean_restart_file !< An unused argument
  character(len=*),      optional, intent(in) :: units !< An unused argument
  character(len=*),      optional, intent(in) :: caller !< An unused argument
  integer,               optional, intent(in) :: verbosity !< An unused argument

  ! None of these arguments are used for anything.

  integer :: coupler_index
  coupler_index = -1

end function aof_set_coupler_flux

end module atmos_ocean_fluxes_mod
