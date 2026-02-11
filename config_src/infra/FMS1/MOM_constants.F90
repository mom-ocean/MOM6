! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides a few physical constants
module MOM_constants

use constants_mod, only : FMS_HLV => HLV
use constants_mod, only : FMS_HLF => HLF

implicit none ; private

real, public, parameter :: CELSIUS_KELVIN_OFFSET = 273.15
  !< The constant offset for converting temperatures in Kelvin to Celsius [K]
real, public, parameter :: HLV = real(FMS_HLV, kind=kind(1.0))
  !< Latent heat of vaporization [J kg-1]
real, public, parameter :: HLF = real(FMS_HLF, kind=kind(1.0))
  !< Latent heat of fusion [J kg-1]

end module MOM_constants
