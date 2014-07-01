module regrid_defs
!==============================================================================
!
! This file is part of MOM.
!
! Date of creation: 2008.12.15
! L. White
!
! This module contains the parameters and types used for the
! regridding/remapping.
!
!==============================================================================
implicit none ; public

! List of reconstruction schemes for pressure gradient calculation
integer, parameter  :: PRESSURE_RECONSTRUCTION_PLM   = 1
integer, parameter  :: PRESSURE_RECONSTRUCTION_PPM   = 2

end module regrid_defs
