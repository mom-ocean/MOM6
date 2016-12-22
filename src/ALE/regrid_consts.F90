!> Contains constants for interpreting input parameters that control regridding.
module regrid_consts

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_string_functions, only : uppercase

implicit none ; public

integer, parameter :: REGRIDDING_NUM_TYPES  = 2

! List of regridding types. These should be consecutive and starting at 1.
! This allows them to be used as array indices.
integer, parameter :: REGRIDDING_LAYER     = 1      !< Layer mode
integer, parameter :: REGRIDDING_ZSTAR     = 2      !< z* coordinates
integer, parameter :: REGRIDDING_RHO       = 3      !< Target interface densities
integer, parameter :: REGRIDDING_SIGMA     = 4      !< Sigma coordinates
integer, parameter :: REGRIDDING_ARBITRARY = 5      !< Arbitrary coordinates
integer, parameter :: REGRIDDING_HYCOM1    = 6      !< Simple HyCOM coordinates without BBL
integer, parameter :: REGRIDDING_SLIGHT    = 7      !< Stretched coordinates in the
integer, parameter :: REGRIDDING_SIGMA_SHELF_ZSTAR = 8   !< z* coordinates at the bottom, sigma-near the top
                                                    !! lightest water, isopycnal below
character(len=6), parameter :: REGRIDDING_LAYER_STRING = "LAYER"   !< Layer string
character(len=6), parameter :: REGRIDDING_ZSTAR_STRING_OLD = "Z*"  !< z* string (legacy name)
character(len=6), parameter :: REGRIDDING_ZSTAR_STRING = "ZSTAR"   !< z* string
character(len=6), parameter :: REGRIDDING_RHO_STRING   = "RHO"     !< Rho string
character(len=6), parameter :: REGRIDDING_SIGMA_STRING = "SIGMA"   !< Sigma string
character(len=6), parameter :: REGRIDDING_ARBITRARY_STRING = "ARB" !< Arbitrary coordinates
character(len=6), parameter :: REGRIDDING_HYCOM1_STRING = "HYCOM1" !< Hycom string
character(len=6), parameter :: REGRIDDING_SLIGHT_STRING = "SLIGHT" !< Hybrid S-rho string
character(len=17), parameter :: REGRIDDING_SIGMA_SHELF_ZSTAR_STRING = "SIGMA_SHELF_ZSTAR" !< Hybrid z*/sigma
character(len=6), parameter :: DEFAULT_COORDINATE_MODE = REGRIDDING_LAYER_STRING !< Default coordinate mode

integer, dimension(REGRIDDING_NUM_TYPES), parameter :: vertical_coords = &
  (/ REGRIDDING_LAYER, REGRIDDING_ZSTAR /)
 !(/ REGRIDDING_LAYER, REGRIDDING_ZSTAR, REGRIDDING_RHO, &
 !  REGRIDDING_SIGMA, REGRIDDING_ARBITRARY, &
 !  REGRIDDING_HYCOM1, REGRIDDING_SLIGHT /)

character(len=*), dimension(REGRIDDING_NUM_TYPES), parameter :: vertical_coord_strings = &
  (/ REGRIDDING_LAYER_STRING, REGRIDDING_ZSTAR_STRING /)
 !(/ REGRIDDING_LAYER_STRING, REGRIDDING_ZSTAR_STRING, REGRIDDING_RHO_STRING, &
 !  REGRIDDING_SIGMA_STRING, REGRIDDING_ARBITRARY_STRING, &
 !  REGRIDDING_HYCOM1_STRING, REGRIDDING_SLIGHT_STRING /)

interface coordinateUnits
  module procedure coordinateUnitsI
  module procedure coordinateUnitsS
end interface

interface state_dependent
  module procedure state_dependent_char
  module procedure state_dependent_int
end interface

contains

!> Parse a string parameter specifying the coordinate mode and
!! return the appropriate enumerated integer
function coordinateMode(string)
  integer :: coordinateMode !< Enumerated integer indicating coordinate mode
  character(len=*), intent(in) :: string !< String to indicate coordinate mode
  select case ( uppercase(trim(string)) )
    case (trim(REGRIDDING_LAYER_STRING)); coordinateMode = REGRIDDING_LAYER
    case (trim(REGRIDDING_ZSTAR_STRING)); coordinateMode = REGRIDDING_ZSTAR
    case (trim(REGRIDDING_ZSTAR_STRING_OLD)); coordinateMode = REGRIDDING_ZSTAR
    case (trim(REGRIDDING_RHO_STRING));   coordinateMode = REGRIDDING_RHO
    case (trim(REGRIDDING_SIGMA_STRING)); coordinateMode = REGRIDDING_SIGMA
    case (trim(REGRIDDING_HYCOM1_STRING)); coordinateMode = REGRIDDING_HYCOM1
    case (trim(REGRIDDING_SLIGHT_STRING)); coordinateMode = REGRIDDING_SLIGHT
    case (trim(REGRIDDING_ARBITRARY_STRING)); coordinateMode = REGRIDDING_ARBITRARY
    case (trim(REGRIDDING_SIGMA_SHELF_ZSTAR_STRING)); coordinateMode = REGRIDDING_SIGMA_SHELF_ZSTAR
    case default ; call MOM_error(FATAL, "coordinateMode: "//&
       "Unrecognized choice of coordinate ("//trim(string)//").")
  end select
end function coordinateMode

!> Returns a string with the coordinate units associated with the
!! enumerated integer,
function coordinateUnitsI(coordMode)
  character(len=16) :: coordinateUnitsI !< Units of coordinate
  integer, intent(in) :: coordMode !< Coordinate mode
  select case ( coordMode )
    case (REGRIDDING_LAYER); coordinateUnitsI = "kg m^-3"
    case (REGRIDDING_ZSTAR); coordinateUnitsI = "m"
    case (REGRIDDING_SIGMA_SHELF_ZSTAR); coordinateUnitsI = "m"
    case (REGRIDDING_RHO);   coordinateUnitsI = "kg m^-3"
    case (REGRIDDING_SIGMA); coordinateUnitsI = "Non-dimensional"
    case (REGRIDDING_HYCOM1); coordinateUnitsI = "m"
    case (REGRIDDING_SLIGHT); coordinateUnitsI = "m"
    case default ; call MOM_error(FATAL, "coordinateUnts: "//&
       "Unrecognized coordinate mode.")
  end select
end function coordinateUnitsI

!> Returns a string with the coordinate units associated with the
!! string defining the coordinate mode.
function coordinateUnitsS(string)
  character(len=16) :: coordinateUnitsS !< Units of coordinate
  character(len=*), intent(in) :: string !< Coordinate mode
  integer :: coordMode
  coordMode = coordinateMode(string)
  coordinateUnitsS = coordinateUnitsI(coordMode)
end function coordinateUnitsS

!> Returns true if the coordinate is dependent on the state density, returns false otherwise.
logical function state_dependent_char(string)
  character(len=*), intent(in) :: string !< String to indicate coordinate mode

  state_dependent_char = state_dependent_int( coordinateMode(string) )

end function state_dependent_char

!> Returns true if the coordinate is dependent on the state density, returns false otherwise.
logical function state_dependent_int(mode)
  integer, intent(in) :: mode !< Coordinate mode
  select case ( mode )
    case (REGRIDDING_LAYER); state_dependent_int = .true.
    case (REGRIDDING_ZSTAR); state_dependent_int = .false.
    case (REGRIDDING_SIGMA_SHELF_ZSTAR); state_dependent_int = .false.
    case (REGRIDDING_RHO);   state_dependent_int = .true.
    case (REGRIDDING_SIGMA); state_dependent_int = .false.
    case (REGRIDDING_HYCOM1); state_dependent_int = .true.
    case (REGRIDDING_SLIGHT); state_dependent_int = .true.
    case default ; call MOM_error(FATAL, "state_dependent: "//&
       "Unrecognized choice of coordinate.")
  end select
end function state_dependent_int

end module regrid_consts
