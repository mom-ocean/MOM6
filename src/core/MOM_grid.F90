module MOM_grid

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

use MOM_domains, only : MOM_domain_type, get_domain_extent
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type

implicit none ; private

#include <MOM_memory.h>

public MOM_grid_init, MOM_grid_end, set_first_direction
public get_flux_units, get_thickness_units, get_tr_flux_units

type, public :: ocean_grid_type
  type(MOM_domain_type), pointer :: Domain => NULL()
  type(MOM_domain_type), pointer :: Domain_aux => NULL()
  integer :: isc, iec, jsc, jec ! The range of the computational domain indicies
  integer :: isd, ied, jsd, jed ! and data domain indicies at tracer cell centers.
  integer :: isg, ieg, jsg, jeg ! The range of the global domain tracer cell indicies.
  integer :: IscB, IecB, JscB, JecB ! The range of the computational domain indicies
  integer :: IsdB, IedB, JsdB, JedB ! and data domain indicies at tracer cell vertices.
  integer :: IsgB, IegB, JsgB, JegB ! The range of the global domain vertex indicies.
  integer :: isd_global         ! The values of isd and jsd in the global
  integer :: jsd_global         ! (decomposition invariant) index space.
  integer :: ks, ke             ! The range of layer's vertical indicies.
  logical :: symmetric          ! True if symmetric memory is used.
  logical :: nonblocking_updates  ! If true, non-blocking halo updates are
                                  ! allowed.  The default is .false. (for now).
  integer :: first_direction ! An integer that indicates which direction is
                             ! to be updated first in directionally split
                             ! parts of the calculation.  This can be altered
                             ! during the course of the run via calls to
                             ! set_first_direction.

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    mask2dT, &   ! 0 for land points and 1 for ocean points on the h-grid. Nd.
    geoLatT, & ! The geographic latitude at q points in degrees of latitude or m.
    geoLonT, & ! The geographic longitude at q points in degrees of longitude or m.
    dxT, IdxT, & ! dxT is delta x at h points, in m, and IdxT is 1/dxT in m-1.
    dyT, IdyT, & ! dyT is delta y at h points, in m, and IdyT is 1/dyT in m-1.
    areaT, &     ! areaT is the area of an h-cell, in m2.
    IareaT       ! IareaT = 1/areaT, in m-2.

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    mask2dCu, &  ! 0 for boundary points and 1 for ocean points on the u grid.  Nondim.
    geoLatCu, &  ! The geographic latitude at u points in degrees of latitude or m.
    geoLonCu, &  ! The geographic longitude at u points in degrees of longitude or m.
    dxCu, IdxCu, & ! dxCu is delta x at u points, in m, and IdxCu is 1/dxCu in m-1.
    dyCu, IdyCu, & ! dyCu is delta y at u points, in m, and IdyCu is 1/dyCu in m-1.
    dy_Cu, &     ! The unblocked lengths of the u-faces of the h-cell in m.
    dy_Cu_obc, & ! The unblocked lengths of the u-faces of the h-cell in m for OBC.
    IareaCu, &   ! The masked inverse areas of u-grid cells in m2.
    areaCu       ! The areas of the u-grid cells in m2.

  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    mask2dCv, &  ! 0 for boundary points and 1 for ocean points on the v grid.  Nondim.
    geoLatCv, &  ! The geographic latitude at v points in degrees of latitude or m.
    geoLonCv, &  !  The geographic longitude at v points in degrees of longitude or m.
    dxCv, IdxCv, & ! dxCv is delta x at v points, in m, and IdxCv is 1/dxCv in m-1.
    dyCv, IdyCv, & ! dyCv is delta y at v points, in m, and IdyCv is 1/dyCv in m-1.
    dx_Cv, &     ! The unblocked lengths of the v-faces of the h-cell in m.
    dx_Cv_obc, & ! The unblocked lengths of the v-faces of the h-cell in m for OBC.
    IareaCv, &   ! The masked inverse areas of v-grid cells in m2.
    areaCv       ! The areas of the v-grid cells in m2.

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    mask2dBu, &  ! 0 for boundary points and 1 for ocean points on the q grid.  Nondim.
    geoLatBu, &  ! The geographic latitude at q points in degrees of latitude or m.
    geoLonBu, &  ! The geographic longitude at q points in degrees of longitude or m.
    dxBu, IdxBu, & ! dxBu is delta x at q points, in m, and IdxBu is 1/dxBu in m-1.
    dyBu, IdyBu, & ! dyBu is delta y at q points, in m, and IdyBu is 1/dyBu in m-1.
    areaBu, &    ! areaBu is the area of a q-cell, in m2
    IareaBu      ! IareaBu = 1/areaBu in m-2.

  real, pointer, dimension(:) :: &
    gridLatT => NULL(), gridLatB => NULL() ! The latitude of T or B points for
                        ! the purpose of labeling the output axes.
                        ! On many grids these are the same as geoLatT & geoLatBu.
  real, pointer, dimension(:) :: &
    gridLonT => NULL(), gridLonB => NULL() ! The longitude of T or B points for
                        ! the purpose of labeling the output axes.
                        ! On many grids these are the same as geoLonT & geoLonBu.
  character(len=40) :: &
    x_axis_units, &     !   The units that are used in labeling the coordinate
    y_axis_units        ! axes.

  ! These parameters are run-time parameters that are used during some
  ! initialization routines (but not all)
  real :: south_lat,   &! The latitude (or y-coordinate) of the first v-line
          west_lon,    &! The longitude (or x-coordinate) of the first u-line
          len_lat = 0.,&! The latitudinal (or y-coord) extent of physical domain
          len_lon = 0.,&! The longitudinal (or x-coord) extent of physical domain
          Rad_Earth = 6.378e6 ! The radius of the planet in meters.
  real :: max_depth     ! The maximum depth of the ocean in meters.
  character(len=40) :: axis_units = ' '! Units for the horizontal coordinates.

  real :: g_Earth !   The gravitational acceleration in m s-2.
  real :: Rho0    !   The density used in the Boussinesq approximation or
                  ! nominal density used to convert depths into mass
                  ! units, in kg m-3.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    bathyT        ! Ocean bottom depth at tracer points, in m.

  logical :: bathymetry_at_vel  ! If true, there are separate values for the
                  ! basin depths at velocity points.  Otherwise the effects of
                  ! of topography are entirely determined from thickness points.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    Dblock_u, &   ! Topographic depths at u-points at which the flow is blocked
    Dopen_u       ! (Dblock_u) and open at width dy_Cu (Dopen_u), both in m.
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    Dblock_v, &   ! Topographic depths at v-points at which the flow is blocked
    Dopen_v       ! (Dblock_v) and open at width dx_Cv (Dopen_v), both in m.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    CoriolisBu    ! The Coriolis parameter at corner points, in s-1.

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

  ! The following are axis types defined for output.
  integer, dimension(3) :: axesBL, axesTL, axesCuL, axesCvL
  integer, dimension(3) :: axesBi, axesTi, axesCui, axesCvi
  integer, dimension(2) :: axesB1, axesT1, axesCu1, axesCv1
  integer, dimension(1) :: axeszi, axeszL

end type ocean_grid_type

contains

subroutine MOM_grid_init(grid, param_file)
  type(ocean_grid_type), intent(inout) :: grid
  type(param_file_type), intent(in)    :: param_file
! Arguments: grid - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
! This include declares and sets the variable "version".
#include "version_variable.h"
  integer :: isd, ied, jsd, jed, nk, idg_off, jdg_off
  integer :: IsdB, IedB, JsdB, JedB

  ! get_domain_extent ensures that domains start at 1 for compatibility between
  ! static and dynamically allocated arrays.
  call get_domain_extent(grid%Domain, grid%isc, grid%iec, grid%jsc, grid%jec, &
                         grid%isd, grid%ied, grid%jsd, grid%jed, &
                         grid%isg, grid%ieg, grid%jsg, grid%jeg, &
                         idg_off, jdg_off, grid%symmetric)
  grid%isd_global = grid%isd+idg_off ; grid%jsd_global = grid%jsd+jdg_off

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, "MOM_grid", version, &
                   "Parameters providing information about the vertical grid.")
  call get_param(param_file, "MOM", "G_EARTH", grid%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)
  call get_param(param_file, "MOM", "RHO_0", grid%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  call get_param(param_file, "MOM_grid", "FIRST_DIRECTION", grid%first_direction, &
                 "An integer that indicates which direction goes first \n"//&
                 "in parts of the code that use directionally split \n"//&
                 "updates, with even numbers (or 0) used for x- first \n"//&
                 "and odd numbers used for y-first.", default=0)
  call get_param(param_file, "MOM_grid", "BOUSSINESQ", grid%Boussinesq, &
                 "If true, make the Boussinesq approximation.", default=.true.)
  call get_param(param_file, "MOM_grid", "ANGSTROM", grid%Angstrom_z, &
                 "The minumum layer thickness, usually one-Angstrom.", &
                 units="m", default=1.0e-10)
  if (.not.grid%Boussinesq) &
    call get_param(param_file, "MOM_grid", "H_TO_KG_M2", grid%H_to_kg_m2,&
                 "A constant that translates thicknesses from the model's \n"//&
                 "internal units of thickness to kg m-2.", units="kg m-2 H-1", &
                 default=1.0)
  call get_param(param_file, "MOM_grid", "BATHYMETRY_AT_VEL", grid%bathymetry_at_vel, &
                 "If true, there are separate values for the basin depths \n"//&
                 "at velocity points.  Otherwise the effects of of \n"//&
                 "topography are entirely determined from thickness points.", &
                 default=.false.)
#ifdef STATIC_MEMORY_
  ! Here NK_ is a macro, while nk is a variable.
  call get_param(param_file, "MOM_grid", "NK", nk, &
                 "The number of model layers.", units="nondim", default=NK_)
  if (nk /= NK_) call MOM_error(FATAL, "MOM_grid_init: " // &
       "Mismatched number of layers NK_ between MOM_memory.h and param_file")

#else
  call get_param(param_file, "MOM_grid", "NK", nk, &
                 "The number of model layers.", units="nondim", fail_if_missing=.true.)
#endif


  grid%ks = 1 ; grid%ke = nk

  grid%nonblocking_updates = grid%Domain%nonblocking_updates

  grid%IscB = grid%isc ; grid%JscB = grid%jsc
  grid%IsdB = grid%isd ; grid%JsdB = grid%jsd
  grid%IsgB = grid%isg ; grid%JsgB = grid%jsg
  if (grid%symmetric) then
    grid%IscB = grid%isc-1 ; grid%JscB = grid%jsc-1
    grid%IsdB = grid%isd-1 ; grid%JsdB = grid%jsd-1
    grid%IsgB = grid%isg-1 ; grid%JsgB = grid%jsg-1
  endif
  grid%IecB = grid%iec ; grid%JecB = grid%jec
  grid%IedB = grid%ied ; grid%JedB = grid%jed
  grid%IegB = grid%ieg ; grid%JegB = grid%jeg


  isd = grid%isd ; ied = grid%ied ; jsd = grid%jsd ; jed = grid%jed
  IsdB = grid%IsdB ; IedB = grid%IedB ; JsdB = grid%JsdB ; JedB = grid%JedB
  ALLOC_(grid%bathyT(isd:ied, jsd:jed)) ; grid%bathyT(:,:) = grid%Angstrom_z
  ALLOC_(grid%CoriolisBu(IsdB:IedB, JsdB:JedB)) ; grid%CoriolisBu(:,:) = 0.0
  ALLOC_(grid%g_prime(nk+1)) ; grid%g_prime(:) = 0.0
  ALLOC_(grid%Rlay(nk+1))    ; grid%Rlay(:) = 0.0

  if (grid%bathymetry_at_vel) then
    ALLOC_(grid%Dblock_u(IsdB:IedB, jsd:jed)) ; grid%Dblock_u(:,:) = 0.0
    ALLOC_(grid%Dopen_u(IsdB:IedB, jsd:jed))  ; grid%Dopen_u(:,:) = 0.0
    ALLOC_(grid%Dblock_v(isd:ied, JsdB:JedB)) ; grid%Dblock_v(:,:) = 0.0
    ALLOC_(grid%Dopen_v(isd:ied, JsdB:JedB))  ; grid%Dopen_v(:,:) = 0.0
  endif

  if (grid%Boussinesq) then
    grid%H_to_kg_m2 = grid%Rho0
    grid%kg_m2_to_H = 1.0/grid%Rho0
    grid%m_to_H = 1.0
    grid%H_to_m = 1.0
    grid%Angstrom = grid%Angstrom_z
  else
!    grid%H_to_kg_m2 = 1.0
!    grid%kg_m2_to_H = 1.0
!    grid%m_to_H = Rho0
!    grid%H_to_m = 1.0 / Rho0
    grid%kg_m2_to_H = 1.0 / grid%H_to_kg_m2
    grid%m_to_H = grid%Rho0 * grid%kg_m2_to_H
    grid%H_to_m = grid%H_to_kg_m2 / grid%Rho0
    grid%Angstrom = grid%Angstrom_z*1000.0*grid%kg_m2_to_H
  endif
  grid%H_subroundoff = 1e-20 * max(grid%Angstrom,grid%m_to_H*1e-17)
  grid%H_to_Pa = grid%g_Earth * grid%H_to_kg_m2

  allocate(grid%gridLatT(grid%Domain%njglobal+2*grid%Domain%njhalo))
  allocate(grid%gridLatB(grid%Domain%njglobal+2*grid%Domain%njhalo))
  grid%gridLatT(:) = 0.0 ; grid%gridLatB(:) = 0.0
  allocate(grid%gridLonT(grid%Domain%niglobal+2*grid%Domain%nihalo))
  allocate(grid%gridLonB(grid%Domain%niglobal+2*grid%Domain%nihalo))
  grid%gridLonT(:) = 0.0 ; grid%gridLonB(:) = 0.0

! Log derivative values.
  call log_param(param_file, "MOM_grid", "M to THICKNESS", grid%m_to_H)

end subroutine MOM_grid_init

subroutine set_first_direction(grid, y_first)
  type(ocean_grid_type), intent(inout) :: grid
  integer,               intent(in) :: y_first

  grid%first_direction = y_first
end subroutine set_first_direction

function get_thickness_units(grid)
  character(len=48)                 :: get_thickness_units
  type(ocean_grid_type), intent(in) :: grid
!   This subroutine returns the appropriate units for thicknesses,
! depending on whether the model is Boussinesq or not and the scaling for
! the vertical thickness.

! Arguments: grid - The ocean's grid structure.
!  (ret)     get_thickness_units - The model's vertical thickness units.

  if (grid%Boussinesq) then
    get_thickness_units = "meter"
  else
    get_thickness_units = "kilogram meter-2"
  endif
end function get_thickness_units

function get_flux_units(grid)
  character(len=48)                 :: get_flux_units
  type(ocean_grid_type), intent(in) :: grid
!   This subroutine returns the appropriate units for thickness fluxes,
! depending on whether the model is Boussinesq or not and the scaling for
! the vertical thickness.

! Arguments: grid - The ocean's grid structure.
!  (ret)     get_flux_units - The model's thickness flux units.

  if (grid%Boussinesq) then
    get_flux_units = "meter3 second-1"
  else
    get_flux_units = "kilogram second-1"
  endif
end function get_flux_units

function get_tr_flux_units(grid, tr_units, tr_vol_conc_units,tr_mass_conc_units)
  character(len=48)                      :: get_tr_flux_units
  type(ocean_grid_type),      intent(in) :: grid
  character(len=*), optional, intent(in) :: tr_units
  character(len=*), optional, intent(in) :: tr_vol_conc_units
  character(len=*), optional, intent(in) :: tr_mass_conc_units
!   This subroutine returns the appropriate units for thicknesses and fluxes,
! depending on whether the model is Boussinesq or not and the scaling for
! the vertical thickness.

! Arguments: grid - The ocean's grid structure.
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
    if (grid%Boussinesq) then
      get_tr_flux_units = trim(tr_units)//" meter3 second-1"
    else
      get_tr_flux_units = trim(tr_units)//" kilogram second-1"
    endif
  endif
  if (present(tr_vol_conc_units)) then
    if (grid%Boussinesq) then
      get_tr_flux_units = trim(tr_vol_conc_units)//" second-1"
    else
      get_tr_flux_units = trim(tr_vol_conc_units)//" m-3 kg s-1"
    endif
  endif
  if (present(tr_mass_conc_units)) then
    if (grid%Boussinesq) then
      get_tr_flux_units = trim(tr_mass_conc_units)//" kg-1 m3 s-1"
    else
      get_tr_flux_units = trim(tr_mass_conc_units)//" second-1"
    endif
  endif

end function get_tr_flux_units

subroutine MOM_grid_end(grid)
! Arguments: grid - The ocean's grid structure.
  type(ocean_grid_type), intent(inout) :: grid

  DEALLOC_(grid%bathyT)  ; DEALLOC_(grid%CoriolisBu)
  DEALLOC_(grid%g_prime) ; DEALLOC_(grid%Rlay)
  deallocate(grid%gridLonT) ; deallocate(grid%gridLatT)
  deallocate(grid%gridLonB) ; deallocate(grid%gridLatB)
end subroutine MOM_grid_end

end module MOM_grid
