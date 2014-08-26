module atmos_ocean_fluxes_mod
!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*    This file is a part of MOM.  See MOM.F90 for licensing.          *
!*                                                                     *
!*    This module is the dummy version of atmos_ocean_fluxes_mod for   *
!*    use when the GFDL-FMS package is not needed.                     *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

implicit none ; private

public :: aof_set_coupler_flux

contains

function aof_set_coupler_flux(name, flux_type, implementation, atm_tr_index,     &
                              param, flag, ice_restart_file, ocean_restart_file, &
                              units, caller)  result (coupler_index)

  character(len=*), intent(in)                          :: name
  character(len=*), intent(in)                          :: flux_type
  character(len=*), intent(in)                          :: implementation
  integer,          intent(in), optional                :: atm_tr_index
  real,             intent(in), dimension(:), optional  :: param
  logical,          intent(in), dimension(:), optional  :: flag
  character(len=*), intent(in), optional                :: ice_restart_file
  character(len=*), intent(in), optional                :: ocean_restart_file
  character(len=*), intent(in), optional                :: units
  character(len=*), intent(in), optional                :: caller

  ! None of these arguments are used for anything.

  integer :: coupler_index
  coupler_index = -1

end function aof_set_coupler_flux

end module atmos_ocean_fluxes_mod
