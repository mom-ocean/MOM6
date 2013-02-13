module regrid_consts
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : uppercase
!==============================================================================
!
! This file is part of MOM.
!
! This module contains constants for interpretting input parameters that
! control regridding.
!
!==============================================================================
implicit none ; public

! List of regridding types
integer, parameter :: REGRIDDING_LAYER     = -1     ! Layer mode
integer, parameter :: REGRIDDING_ZSTAR     = 0      ! z* coordinates
integer, parameter :: REGRIDDING_RHO       = 1      ! target interface densities
integer, parameter :: REGRIDDING_SIGMA     = 2      ! sigma coordinates
integer, parameter :: REGRIDDING_ARBITRARY = 3      ! arbitrary coordinates
character(len=5), parameter :: REGRIDDING_LAYER_STRING = "LAYER"
character(len=2), parameter :: REGRIDDING_ZSTAR_STRING = "Z*"
character(len=3), parameter :: REGRIDDING_RHO_STRING   = "RHO"
character(len=5), parameter :: REGRIDDING_SIGMA_STRING = "SIGMA"
character(len=5), parameter :: DEFAULT_COORDINATE_MODE = REGRIDDING_LAYER_STRING

contains

! =============================================================================

function coordinateMode(string)
! Use this function to parse a string parameter specifying the
! coorindate mode and return the appropriate enumerated integer
  integer :: coordinateMode ! 
  character(len=*), intent(in) :: string
  select case ( uppercase(trim(string)) )
    case (REGRIDDING_LAYER_STRING); coordinateMode = REGRIDDING_LAYER
    case (REGRIDDING_ZSTAR_STRING); coordinateMode = REGRIDDING_ZSTAR
    case (REGRIDDING_RHO_STRING);   coordinateMode = REGRIDDING_RHO
    case (REGRIDDING_SIGMA_STRING); coordinateMode = REGRIDDING_SIGMA
    case default ; call MOM_error(FATAL, "coordinateMode: "//&
       "Unrecognized choice of coordinate ("//trim(string)//").")
  end select
end function coordinateMode

end module regrid_consts
