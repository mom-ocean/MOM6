!> Provides a few physical constants
module MOM_constants

! This file is part of MOM6. See LICENSE.md for the license.

use constants_mod, only : HLV, HLF

implicit none ; private

!> The constant offset for converting temperatures in Kelvin to Celsius
real, public, parameter :: CELSIUS_KELVIN_OFFSET = 273.15
public :: HLV, HLF

end module MOM_constants
