module MOM_verticalGrid

!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type

implicit none ; private

#include <MOM_memory.h>

public verticalGridInit, verticalGridEnd
public setVerticalGridAxes
public get_flux_units, get_thickness_units, get_tr_flux_units

type, public :: verticalGrid_type

  ! Commonly used parameters
  integer :: ke     ! The number of layers/levels in the vertical
  real :: max_depth ! The maximum depth of the ocean in meters.
  real :: g_Earth   ! The gravitational acceleration in m s-2.
  real :: Rho0      !   The density used in the Boussinesq approximation or
                    ! nominal density used to convert depths into mass
                    ! units, in kg m-3.

  ! Vertical coordinate descriptions for diagnostics and I/O
  character(len=40) :: &
    zAxisUnits, & ! The units that vertical coordinates are written in
    zAxisLongName ! Coordinate name to appear in files,
                  ! e.g. "Target Potential Density" or "Height"
  real ALLOCABLE_, dimension(NKMEM_) :: sLayer ! Coordinate values of layer centers
  real ALLOCABLE_, dimension(NK_INTERFACE_) :: sInterface ! Coordinate values on interfaces
  integer :: direction = 1 ! Direction defaults to 1, positive up.

  ! The following variables give information about the vertical grid.
  logical :: Boussinesq     ! If true, make the Boussinesq approximation.
  real :: Angstrom      !   A one-Angstrom thickness in the model's thickness
                        ! units.  (This replaces the old macro EPSILON.)
  real :: Angstrom_z    !   A one-Angstrom thickness in m.
  real :: H_subroundoff !   A thickness that is so small that it can be added to
                        ! a thickness of Angstrom or larger without changing it
                        ! at the bit level, in thickness units.  If Angstrom is
                        ! 0 or exceedingly small, this is negligible compared to
                        ! a thickness of 1e-17 m.
  real ALLOCABLE_, dimension(NK_INTERFACE_) :: &
    g_prime, &          ! The reduced gravity at each interface, in m s-2.
    Rlay                ! The target coordinate value (potential density) in
                        ! in each layer in kg m-3.
  integer :: nkml = 0   ! The number of layers at the top that should be treated
                        ! as parts of a homogenous region.
  integer :: nk_rho_varies = 0 ! The number of layers at the top where the
                        ! density does not track any target density.
  real :: H_to_kg_m2    ! A constant that translates thicknesses from the units
                        ! of thickness to kg m-2.
  real :: kg_m2_to_H    ! A constant that translates thicknesses from kg m-2 to
                        ! the units of thickness.
  real :: m_to_H        ! A constant that translates distances in m to the
                        ! units of thickness.
  real :: H_to_m        ! A constant that translates distances in the units of
                        ! thickness to m.
  real :: H_to_Pa       ! A constant that translates the units of thickness to
                        ! to pressure in Pa.
end type verticalGrid_type

contains

!> Allocates and initializes the model's vertical grid structure.
subroutine verticalGridInit( param_file, GV )
! This routine initializes the verticalGrid_type structure (GV).
! All memory is allocated but not necessarily set to meaningful values until later.
  type(param_file_type),   intent(in) :: param_file ! Parameter file handle/type
  type(verticalGrid_type), pointer    :: GV         ! The container for vertical grid data
! This include declares and sets the variable "version".
#include "version_variable.h"
  integer :: nk
  character(len=16) :: mod = 'MOM_verticalGrid'

  if (associated(GV)) call MOM_error(FATAL, &
     'verticalGridInit: called with an associated GV pointer.')
  allocate(GV)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, &
                   "Parameters providing information about the vertical grid.")
  call get_param(param_file, mod, "G_EARTH", GV%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)
  call get_param(param_file, mod, "RHO_0", GV%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  call get_param(param_file, mod, "BOUSSINESQ", GV%Boussinesq, &
                 "If true, make the Boussinesq approximation.", default=.true.)
  call get_param(param_file, mod, "ANGSTROM", GV%Angstrom_z, &
                 "The minumum layer thickness, usually one-Angstrom.", &
                 units="m", default=1.0e-10)
  if (.not.GV%Boussinesq) then
    call get_param(param_file, mod, "H_TO_KG_M2", GV%H_to_kg_m2,&
                 "A constant that translates thicknesses from the model's \n"//&
                 "internal units of thickness to kg m-2.", units="kg m-2 H-1", &
                 default=1.0)
  else
    call get_param(param_file, mod, "H_TO_M", GV%H_to_m, &
                 "A constant that translates the model's internal \n"//&
                 "units of thickness into m.", units="m H-1", default=1.0)
  endif
#ifdef STATIC_MEMORY_
  ! Here NK_ is a macro, while nk is a variable.
  call get_param(param_file, mod, "NK", nk, &
                 "The number of model layers.", units="nondim", &
                 static_value=NK_)
  if (nk /= NK_) call MOM_error(FATAL, "verticalGridInit: " // &
       "Mismatched number of layers NK_ between MOM_memory.h and param_file")

#else
  call get_param(param_file, mod, "NK", nk, &
                 "The number of model layers.", units="nondim", fail_if_missing=.true.)
#endif
  GV%ke = nk

  if (GV%Boussinesq) then
    GV%H_to_kg_m2 = GV%Rho0 * GV%H_to_m
    GV%kg_m2_to_H = 1.0 / GV%H_to_kg_m2
    GV%m_to_H = 1.0 / GV%H_to_m
    GV%Angstrom = GV%Angstrom_z
  else
    GV%kg_m2_to_H = 1.0 / GV%H_to_kg_m2
    GV%m_to_H = GV%Rho0 * GV%kg_m2_to_H
    GV%H_to_m = GV%H_to_kg_m2 / GV%Rho0
    GV%Angstrom = GV%Angstrom_z*1000.0*GV%kg_m2_to_H
  endif
  GV%H_subroundoff = 1e-20 * max(GV%Angstrom,GV%m_to_H*1e-17)
  GV%H_to_Pa = GV%g_Earth * GV%H_to_kg_m2

! Log derivative values.
  call log_param(param_file, mod, "M to THICKNESS", GV%m_to_H)

  ALLOC_( GV%sInterface(nk+1) )
  ALLOC_( GV%sLayer(nk) )
  ALLOC_( GV%g_prime(nk+1) ) ; GV%g_prime(:) = 0.0
  ! The extent of Rlay should be changed to nk?
  ALLOC_( GV%Rlay(nk+1) )    ; GV%Rlay(:) = 0.0

end subroutine verticalGridInit

!> Returns the model's thickness units, usually m or kg/m^2.
function get_thickness_units(GV)
  character(len=48)                 :: get_thickness_units
  type(verticalGrid_type), intent(in) :: GV
!   This subroutine returns the appropriate units for thicknesses,
! depending on whether the model is Boussinesq or not and the scaling for
! the vertical thickness.

! Arguments: G - The ocean's grid structure.
!  (ret)     get_thickness_units - The model's vertical thickness units.

  if (GV%Boussinesq) then
    get_thickness_units = "meter"
  else
    get_thickness_units = "kilogram meter-2"
  endif
end function get_thickness_units

!> Returns the model's thickness flux units, usually m^3/s or kg/s.
function get_flux_units(GV)
  character(len=48)                 :: get_flux_units
  type(verticalGrid_type), intent(in) :: GV
!   This subroutine returns the appropriate units for thickness fluxes,
! depending on whether the model is Boussinesq or not and the scaling for
! the vertical thickness.

! Arguments: G - The ocean's grid structure.
!  (ret)     get_flux_units - The model's thickness flux units.

  if (GV%Boussinesq) then
    get_flux_units = "meter3 second-1"
  else
    get_flux_units = "kilogram second-1"
  endif
end function get_flux_units

!> Returns the model's tracer flux units.
function get_tr_flux_units(GV, tr_units, tr_vol_conc_units, tr_mass_conc_units)
  character(len=48)                      :: get_tr_flux_units
  type(verticalGrid_type),    intent(in) :: GV
  character(len=*), optional, intent(in) :: tr_units
  character(len=*), optional, intent(in) :: tr_vol_conc_units
  character(len=*), optional, intent(in) :: tr_mass_conc_units
!   This subroutine returns the appropriate units for thicknesses and fluxes,
! depending on whether the model is Boussinesq or not and the scaling for
! the vertical thickness.

! Arguments: G - The ocean's grid structure.
!      One of the following three arguments must be present.
!  (in,opt)  tr_units - Units for a tracer, for example Celsius or PSU.
!  (in,opt)  tr_vol_conc_units - The concentration units per unit volume, for
!                                example if the units are umol m-3,
!                                tr_vol_conc_units would be umol.
!  (in,opt)  tr_mass_conc_units - The concentration units per unit mass of sea
!                                water, for example if the units are mol kg-1,
!                                tr_vol_conc_units would be mol.
!  (ret)     get_tr_flux_units - The model's flux units for a tracer.
  integer :: cnt

  cnt = 0
  if (present(tr_units)) cnt = cnt+1
  if (present(tr_vol_conc_units)) cnt = cnt+1
  if (present(tr_mass_conc_units)) cnt = cnt+1

  if (cnt == 0) call MOM_error(FATAL, "get_tr_flux_units: One of the three "//&
    "arguments tr_units, tr_vol_conc_units, or tr_mass_conc_units "//&
    "must be present.")
  if (cnt > 1) call MOM_error(FATAL, "get_tr_flux_units: Only one of "//&
    "tr_units, tr_vol_conc_units, and tr_mass_conc_units may be present.")
  if (present(tr_units)) then
    if (GV%Boussinesq) then
      get_tr_flux_units = trim(tr_units)//" meter3 second-1"
    else
      get_tr_flux_units = trim(tr_units)//" kilogram second-1"
    endif
  endif
  if (present(tr_vol_conc_units)) then
    if (GV%Boussinesq) then
      get_tr_flux_units = trim(tr_vol_conc_units)//" second-1"
    else
      get_tr_flux_units = trim(tr_vol_conc_units)//" m-3 kg s-1"
    endif
  endif
  if (present(tr_mass_conc_units)) then
    if (GV%Boussinesq) then
      get_tr_flux_units = trim(tr_mass_conc_units)//" kg-1 m3 s-1"
    else
      get_tr_flux_units = trim(tr_mass_conc_units)//" second-1"
    endif
  endif

end function get_tr_flux_units

!> This sets the coordinate data for the "layer mode" of the isopycnal model.
subroutine setVerticalGridAxes( Rlay, GV )
  ! Arguments
  type(verticalGrid_type), intent(inout) :: GV   !< The container for vertical grid data
  real, dimension(GV%ke),  intent(in)    :: Rlay !< The layer target density
  ! Local variables
  integer :: k, nk

  nk = GV%ke

  GV%zAxisLongName = 'Target Potential Density'
  GV%zAxisUnits = 'kg m-3'
  do k=1,nk ; GV%sLayer(k) = Rlay(k) ; enddo
  if (nk > 1) then
    GV%sInterface(1) = 1.5*Rlay(1) - 0.5*Rlay(2)
    do K=2,nk ; GV%sInterface(K) = 0.5*( Rlay(k-1) + Rlay(k) ) ; enddo
    GV%sInterface(nk+1) = 1.5*Rlay(nk) - 0.5*Rlay(nk-1)
  else
    GV%sInterface(1) = 0.0 ; GV%sInterface(nk+1) = 2.0*Rlay(nk)
  endif

end subroutine setVerticalGridAxes

!> Deallocates the model's vertical grid structure.
subroutine verticalGridEnd( GV )
! Arguments: G - The ocean's grid structure.
  type(verticalGrid_type), pointer :: GV

  DEALLOC_(GV%g_prime) ; DEALLOC_(GV%Rlay)
  DEALLOC_( GV%sInterface )
  DEALLOC_( GV%sLayer )
  deallocate( GV )

end subroutine verticalGridEnd

end module MOM_verticalGrid
