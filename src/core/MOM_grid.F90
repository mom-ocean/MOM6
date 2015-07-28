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

use MOM_domains, only : MOM_domain_type, get_domain_extent, compute_block_extent
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_verticalGrid, only : verticalGrid_type
use MOM_verticalGrid, only : verticalGridInit, verticalGridEnd

implicit none ; private

#include <MOM_memory.h>

public MOM_grid_init, MOM_grid_end, set_first_direction
public get_flux_units, get_thickness_units, get_tr_flux_units
public isPointInCell

type, public :: ocean_block_type
  integer :: isc, iec, jsc, jec     ! The range of the computational indices and 
  integer :: isd, ied, jsd, jed     ! data indices at tracer cell enters for each block.
  integer :: IscB, IecB, JscB, JecB ! The range of the computational indices and 
  integer :: IsdB, IedB, JsdB, JedB ! data indices at tracer cell vertices 
                                    ! for each block.
  integer :: ioff, joff             ! index offset of block indices relative to 
                                    ! domain indices
end type ocean_block_type

type, public :: ocean_grid_type
  type(MOM_domain_type), pointer :: Domain => NULL()
  type(MOM_domain_type), pointer :: Domain_aux => NULL()
  type(verticalGrid_type), pointer :: GV => NULL()
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
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    dF_dx, dF_dy  ! Derivatives of f (Coriolis parameter) at h-points, in s-1 m-1.

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

  ! These variables are global sums that are useful for 1-d diagnostics
  real :: areaT_global  ! Global sum of h-cell area in m2
  real :: IareaT_global ! Global sum of inverse h-cell area (1/areaT_global)
                        ! in m2
  ! These variables are for block strucutre.
  integer                   :: nblocks
  type(ocean_block_type), pointer :: Block(:) => NULL() ! store indices for each block
  integer :: isd_bk, ied_bk, jsd_bk, jed_bk     ! block data domain indices at 
                                                ! tracer cell centers.
  integer :: isdB_bk, iedB_bk, jsdB_bk, jedB_bk ! block data domain indices at 
                                                ! tracer cell vertices.
end type ocean_grid_type

contains

subroutine MOM_grid_init(G, param_file)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
! This include declares and sets the variable "version".
#include "version_variable.h"
  integer :: isd, ied, jsd, jed, nk, idg_off, jdg_off
  integer :: IsdB, IedB, JsdB, JedB
  integer :: niblock, njblock, nihalo, njhalo, nblocks, n, i, j
  integer, allocatable, dimension(:) :: ibegin, iend, jbegin, jend

  ! get_domain_extent ensures that domains start at 1 for compatibility between
  ! static and dynamically allocated arrays.
  call get_domain_extent(G%Domain, G%isc, G%iec, G%jsc, G%jec, &
                         G%isd, G%ied, G%jsd, G%jed, &
                         G%isg, G%ieg, G%jsg, G%jeg, &
                         idg_off, jdg_off, G%symmetric)
  G%isd_global = G%isd+idg_off ; G%jsd_global = G%jsd+jdg_off

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, "MOM_grid", version, &
                   "Parameters providing information about the vertical grid.")
  call get_param(param_file, "MOM", "G_EARTH", G%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)
  call get_param(param_file, "MOM", "RHO_0", G%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  call get_param(param_file, "MOM_grid", "FIRST_DIRECTION", G%first_direction, &
                 "An integer that indicates which direction goes first \n"//&
                 "in parts of the code that use directionally split \n"//&
                 "updates, with even numbers (or 0) used for x- first \n"//&
                 "and odd numbers used for y-first.", default=0)
  call get_param(param_file, "MOM_grid", "BOUSSINESQ", G%Boussinesq, &
                 "If true, make the Boussinesq approximation.", default=.true.)
  call get_param(param_file, "MOM_grid", "ANGSTROM", G%Angstrom_z, &
                 "The minumum layer thickness, usually one-Angstrom.", &
                 units="m", default=1.0e-10)

  if (.not.G%Boussinesq) then
    call get_param(param_file, "MOM_grid", "H_TO_KG_M2", G%H_to_kg_m2,&
                 "A constant that translates thicknesses from the model's \n"//&
                 "internal units of thickness to kg m-2.", units="kg m-2 H-1", &
                 default=1.0)
  else
    call get_param(param_file, "MOM_grid", "H_TO_M", G%H_to_m, default=1.0)
  endif

  call get_param(param_file, "MOM_grid", "BATHYMETRY_AT_VEL", G%bathymetry_at_vel, &
                 "If true, there are separate values for the basin depths \n"//&
                 "at velocity points.  Otherwise the effects of of \n"//&
                 "topography are entirely determined from thickness points.", &
                 default=.false.)
#ifdef STATIC_MEMORY_
  ! Here NK_ is a macro, while nk is a variable.
  call get_param(param_file, "MOM_grid", "NK", nk, &
                 "The number of model layers.", units="nondim", &
                 static_value=NK_)
  if (nk /= NK_) call MOM_error(FATAL, "MOM_grid_init: " // &
       "Mismatched number of layers NK_ between MOM_memory.h and param_file")
  niblock = NIBLOCK_
  njblock = NJBLOCK_
  call log_param(param_file, "MOM_grid", "NIBLOCK", niblock, "The number of blocks "// &
                 "in the x-direction on each processor (for openmp).", default=1, &
                 layoutParam=.true.)
  call log_param(param_file, "MOM_grid", "NJBLOCK", njblock, "The number of blocks "// &
                 "in the y-direction on each processor (for openmp).", default=1, &
                 layoutParam=.true.)
#else
  call get_param(param_file, "MOM_grid", "NK", nk, &
                 "The number of model layers.", units="nondim", fail_if_missing=.true.)
  call get_param(param_file, "MOM_grid", "NIBLOCK", niblock, "The number of blocks "// &
                 "in the x-direction on each processor (for openmp).", default=1, &
                 layoutParam=.true.)
  call get_param(param_file, "MOM_grid", "NJBLOCK", njblock, "The number of blocks "// &
                 "in the y-direction on each processor (for openmp).", default=1, &
                 layoutParam=.true.)
#endif

  G%ks = 1 ; G%ke = nk

  G%nonblocking_updates = G%Domain%nonblocking_updates

  G%IscB = G%isc ; G%JscB = G%jsc
  G%IsdB = G%isd ; G%JsdB = G%jsd
  G%IsgB = G%isg ; G%JsgB = G%jsg
  if (G%symmetric) then
    G%IscB = G%isc-1 ; G%JscB = G%jsc-1
    G%IsdB = G%isd-1 ; G%JsdB = G%jsd-1
    G%IsgB = G%isg-1 ; G%JsgB = G%jsg-1
  endif
  G%IecB = G%iec ; G%JecB = G%jec
  G%IedB = G%ied ; G%JedB = G%jed
  G%IegB = G%ieg ; G%JegB = G%jeg

  call MOM_mesg("  MOM_grid.F90, MOM_grid_init: allocating metrics", 5)
 
  call allocate_metrics(G)

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  ALLOC_(G%bathyT(isd:ied, jsd:jed)) ; G%bathyT(:,:) = G%Angstrom_z
  ALLOC_(G%g_prime(nk+1)) ; G%g_prime(:) = 0.0
  ALLOC_(G%Rlay(nk+1))    ; G%Rlay(:) = 0.0

  if (G%bathymetry_at_vel) then
    ALLOC_(G%Dblock_u(IsdB:IedB, jsd:jed)) ; G%Dblock_u(:,:) = 0.0
    ALLOC_(G%Dopen_u(IsdB:IedB, jsd:jed))  ; G%Dopen_u(:,:) = 0.0
    ALLOC_(G%Dblock_v(isd:ied, JsdB:JedB)) ; G%Dblock_v(:,:) = 0.0
    ALLOC_(G%Dopen_v(isd:ied, JsdB:JedB))  ; G%Dopen_v(:,:) = 0.0
  endif

  if (G%Boussinesq) then
    G%H_to_kg_m2 = G%Rho0 * G%H_to_m
    G%kg_m2_to_H = 1.0 / G%H_to_kg_m2
    G%m_to_H = 1.0 / G%H_to_m
    G%Angstrom = G%Angstrom_z*G%m_to_H
  else
    G%kg_m2_to_H = 1.0 / G%H_to_kg_m2
    G%m_to_H = G%Rho0 * G%kg_m2_to_H
    G%H_to_m = G%H_to_kg_m2 / G%Rho0
    G%Angstrom = G%Angstrom_z*1000.0*G%kg_m2_to_H
  endif
  G%H_subroundoff = 1e-20 * max(G%Angstrom,G%m_to_H*1e-17)
  G%H_to_Pa = G%g_Earth * G%H_to_kg_m2

  allocate(G%gridLatT(G%jsg:G%jeg))
  allocate(G%gridLatB(G%JsgB:G%JegB))
  G%gridLatT(:) = 0.0 ; G%gridLatB(:) = 0.0
  allocate(G%gridLonT(G%isg:G%ieg))
  allocate(G%gridLonB(G%IsgB:G%IegB))
  G%gridLonT(:) = 0.0 ; G%gridLonB(:) = 0.0

! setup block indices.
  nihalo = G%Domain%nihalo
  njhalo = G%Domain%njhalo
  nblocks = niblock * njblock
  if (nblocks < 1) call MOM_error(FATAL, "MOM_grid_init: " // &
       "nblocks(=NI_BLOCK*NJ_BLOCK) must be no less than 1")

  allocate(ibegin(niblock), iend(niblock), jbegin(njblock), jend(njblock))
  call compute_block_extent(G%isc,G%iec,niblock,ibegin,iend)
  call compute_block_extent(G%jsc,G%jec,njblock,jbegin,jend)

  G%nblocks = nblocks
  allocate(G%Block(nblocks))
  do n = 1,nblocks
    i = mod((n-1), niblock) + 1
    j = (n-1)/niblock + 1
    !--- isd and jsd are always 1 for each block
    G%Block(n)%isd = 1 ; G%Block(n)%jsd = 1
    G%Block(n)%isc = G%Block(n)%isd+nihalo 
    G%Block(n)%jsc = G%Block(n)%jsd+njhalo
    G%Block(n)%iec = G%Block(n)%isc + iend(i) - ibegin(i)
    G%Block(n)%jec = G%Block(n)%jsc + jend(j) - jbegin(j)
    G%Block(n)%ied = G%Block(n)%iec + nihalo
    G%Block(n)%jed = G%Block(n)%jec + njhalo
    G%Block(n)%IscB = G%Block(n)%isc; G%Block(n)%IecB = G%Block(n)%iec
    G%Block(n)%JscB = G%Block(n)%jsc; G%Block(n)%JecB = G%Block(n)%jec
    !--- For symmetry domain, the first block will have the extra point.
    if (G%symmetric) then
      if (i==1) G%Block(n)%IscB = G%Block(n)%IscB-1
      if (j==1) G%Block(n)%JscB = G%Block(n)%JscB-1
    endif
    G%Block(n)%IsdB = G%Block(n)%isd; G%Block(n)%IedB = G%Block(n)%ied
    G%Block(n)%JsdB = G%Block(n)%jsd; G%Block(n)%JedB = G%Block(n)%jed
    !--- For symmetry domain, every block will have the extra point.
    if (G%symmetric) then
      G%Block(n)%IsdB = G%Block(n)%IsdB-1
      G%Block(n)%JsdB = G%Block(n)%JsdB-1
    endif
    G%Block(n)%ioff = ibegin(i) - G%Block(n)%isc
    G%Block(n)%joff = jbegin(j) - G%Block(n)%jsc
  enddo

  !-- make sure the last block is the largest.
  do i = 1, niblock-1
    if (iend(i)-ibegin(i) > iend(niblock)-ibegin(niblock) ) call MOM_error(FATAL, &
       "MOM_grid_init: the last block size in x-direction is not the largest")
  enddo
  do j = 1, njblock-1
    if (jend(j)-jbegin(j) > jend(njblock)-jbegin(njblock) ) call MOM_error(FATAL, &
       "MOM_grid_init: the last block size in y-direction is not the largest")
  enddo

  !-- make sure

  !--- define the block memory domain ( maximum data domain size of all blocks )
  G%isd_bk  = G%block(nblocks)%isd  ; G%ied_bk  = G%block(nblocks)%ied
  G%jsd_bk  = G%block(nblocks)%jsd  ; G%jed_bk  = G%block(nblocks)%jed
  G%isdB_bk = G%block(nblocks)%isdB ; G%iedB_bk = G%block(nblocks)%iedB
  G%jsdB_bk = G%block(nblocks)%jsdB ; G%jedB_bk = G%block(nblocks)%jedB

  !-- do some bound check
  if ( G%ied_bk+G%block(nblocks)%ioff > G%ied ) call MOM_error(FATAL, &
       "MOM_grid_init: G%ied_bk+G%block(nblocks)%ioff > G%ied")
  if ( G%jed_bk+G%block(nblocks)%joff > G%jed ) call MOM_error(FATAL, &
       "MOM_grid_init: G%jed_bk+G%block(nblocks)%joff > G%jed")

  !--- For static memory, make sure G%iem_bk - G%ism_bk + 1 = NI_MEM_BK_
  !---                         and  G%jem_bk - G%jsm_bk + 1 = NJ_MEM_BK_
#ifdef STATIC_MEMORY_
  if ( (G%ied_bk-G%isd_bk+1) .NE. NIMEM_BK_ ) call MOM_error(FATAL, &
       "MOM_grid_init:  (G%ied_bk-G%isd_bk+1) .NE. NIMEM_BK_ for static memory ")
  if ( (G%jed_bk-G%jsd_bk+1) .NE. NJMEM_BK_ ) call MOM_error(FATAL, &
       "MOM_grid_init:  (G%jed_bk-G%jsd_bk+1) .NE. NJMEM_BK_ for static memory ")
#endif

! Log derivative values.
  call log_param(param_file, "MOM_grid", "M to THICKNESS", G%m_to_H)

  call verticalGridInit( param_file, G%GV )

end subroutine MOM_grid_init

!> Returns true if the coordinates (x,y) are within the h-cell (i,j)
logical function isPointInCell(G, i, j, x, y)
  type(ocean_grid_type), intent(in) :: G    !< Grid type
  integer,               intent(in) :: i, j !< i,j indices of cell to test
  real,                  intent(in) :: x, y !< x,y coordinates of point
! This is a crude calculation that assume a geographic coordinate system
  real :: xNE, xNW, xSE, xSW, yNE, yNW, ySE, ySW
  real :: p0, p1, p2, p3, l0, l1, l2, l3
  isPointInCell = .false.
  xNE = G%geoLonBu(i  ,j  ); yNE = G%geoLatBu(i  ,j  )
  xNW = G%geoLonBu(i-1,j  ); yNW = G%geoLatBu(i-1,j  )
  xSE = G%geoLonBu(i  ,j-1); ySE = G%geoLatBu(i  ,j-1)
  xSW = G%geoLonBu(i-1,j-1); ySW = G%geoLatBu(i-1,j-1)
  if (x<min(xNE,xNW,xSE,xSW) .or. x>max(xNE,xNW,xSE,xSW) .or. &
      y<min(yNE,yNW,ySE,ySW) .or. y>max(yNE,yNW,ySE,ySW) ) then
    return ! Avoid the more complicated calculation
  endif
  l0=(x-xSW)*(ySE-ySW)-(y-ySW)*(xSE-xSW)
  l1=(x-xSE)*(yNE-ySE)-(y-ySE)*(xNE-xSE)
  l2=(x-xNE)*(yNW-yNE)-(y-yNE)*(xNW-xNE)
  l3=(x-xNW)*(ySW-yNW)-(y-yNW)*(xSW-xNW)

  p0=sign(1., l0); if (l0.eq.0.) p0=0.
  p1=sign(1., l1); if (l1.eq.0.) p1=0.
  p2=sign(1., l2); if (l2.eq.0.) p2=0.
  p3=sign(1., l3); if (l3.eq.0.) p3=0.

  if ( (abs(p0)+abs(p2))+(abs(p1)+abs(p3)) .eq. abs((p0+p2)+(p1+p3)) ) then
    isPointInCell=.true.
  endif
end function isPointInCell

subroutine set_first_direction(G, y_first)
  type(ocean_grid_type), intent(inout) :: G
  integer,               intent(in) :: y_first

  G%first_direction = y_first
end subroutine set_first_direction

function get_thickness_units(G)
  character(len=48)                 :: get_thickness_units
  type(ocean_grid_type), intent(in) :: G
!   This subroutine returns the appropriate units for thicknesses,
! depending on whether the model is Boussinesq or not and the scaling for
! the vertical thickness.

! Arguments: G - The ocean's grid structure.
!  (ret)     get_thickness_units - The model's vertical thickness units.

  if (G%Boussinesq) then
    get_thickness_units = "meter"
  else
    get_thickness_units = "kilogram meter-2"
  endif
end function get_thickness_units

function get_flux_units(G)
  character(len=48)                 :: get_flux_units
  type(ocean_grid_type), intent(in) :: G
!   This subroutine returns the appropriate units for thickness fluxes,
! depending on whether the model is Boussinesq or not and the scaling for
! the vertical thickness.

! Arguments: G - The ocean's grid structure.
!  (ret)     get_flux_units - The model's thickness flux units.

  if (G%Boussinesq) then
    get_flux_units = "meter3 second-1"
  else
    get_flux_units = "kilogram second-1"
  endif
end function get_flux_units

function get_tr_flux_units(G, tr_units, tr_vol_conc_units,tr_mass_conc_units)
  character(len=48)                      :: get_tr_flux_units
  type(ocean_grid_type),      intent(in) :: G
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
    if (G%Boussinesq) then
      get_tr_flux_units = trim(tr_units)//" meter3 second-1"
    else
      get_tr_flux_units = trim(tr_units)//" kilogram second-1"
    endif
  endif
  if (present(tr_vol_conc_units)) then
    if (G%Boussinesq) then
      get_tr_flux_units = trim(tr_vol_conc_units)//" second-1"
    else
      get_tr_flux_units = trim(tr_vol_conc_units)//" m-3 kg s-1"
    endif
  endif
  if (present(tr_mass_conc_units)) then
    if (G%Boussinesq) then
      get_tr_flux_units = trim(tr_mass_conc_units)//" kg-1 m3 s-1"
    else
      get_tr_flux_units = trim(tr_mass_conc_units)//" second-1"
    endif
  endif

end function get_tr_flux_units


subroutine allocate_metrics(G)
  type(ocean_grid_type), intent(inout) :: G
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  ! This subroutine allocates the lateral elements of the ocean_grid_type that
  ! are always used and zeros them out.

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ALLOC_(G%dxT(isd:ied,jsd:jed))       ; G%dxT(:,:) = 0.0
  ALLOC_(G%dxCu(IsdB:IedB,jsd:jed))    ; G%dxCu(:,:) = 0.0
  ALLOC_(G%dxCv(isd:ied,JsdB:JedB))    ; G%dxCv(:,:) = 0.0
  ALLOC_(G%dxBu(IsdB:IedB,JsdB:JedB))  ; G%dxBu(:,:) = 0.0
  ALLOC_(G%IdxT(isd:ied,jsd:jed))      ; G%IdxT(:,:) = 0.0
  ALLOC_(G%IdxCu(IsdB:IedB,jsd:jed))   ; G%IdxCu(:,:) = 0.0
  ALLOC_(G%IdxCv(isd:ied,JsdB:JedB))   ; G%IdxCv(:,:) = 0.0
  ALLOC_(G%IdxBu(IsdB:IedB,JsdB:JedB)) ; G%IdxBu(:,:) = 0.0

  ALLOC_(G%dyT(isd:ied,jsd:jed))       ; G%dyT(:,:) = 0.0
  ALLOC_(G%dyCu(IsdB:IedB,jsd:jed))    ; G%dyCu(:,:) = 0.0
  ALLOC_(G%dyCv(isd:ied,JsdB:JedB))    ; G%dyCv(:,:) = 0.0
  ALLOC_(G%dyBu(IsdB:IedB,JsdB:JedB))  ; G%dyBu(:,:) = 0.0
  ALLOC_(G%IdyT(isd:ied,jsd:jed))      ; G%IdyT(:,:) = 0.0
  ALLOC_(G%IdyCu(IsdB:IedB,jsd:jed))   ; G%IdyCu(:,:) = 0.0
  ALLOC_(G%IdyCv(isd:ied,JsdB:JedB))   ; G%IdyCv(:,:) = 0.0
  ALLOC_(G%IdyBu(IsdB:IedB,JsdB:JedB)) ; G%IdyBu(:,:) = 0.0

  ALLOC_(G%areaT(isd:ied,jsd:jed))       ; G%areaT(:,:) = 0.0
  ALLOC_(G%IareaT(isd:ied,jsd:jed))      ; G%IareaT(:,:) = 0.0
  ALLOC_(G%areaBu(IsdB:IedB,JsdB:JedB))  ; G%areaBu(:,:) = 0.0
  ALLOC_(G%IareaBu(IsdB:IedB,JsdB:JedB)) ; G%IareaBu(:,:) = 0.0

  ALLOC_(G%mask2dT(isd:ied,jsd:jed))      ; G%mask2dT(:,:) = 0.0
  ALLOC_(G%mask2dCu(IsdB:IedB,jsd:jed))   ; G%mask2dCu(:,:) = 0.0
  ALLOC_(G%mask2dCv(isd:ied,JsdB:JedB))   ; G%mask2dCv(:,:) = 0.0
  ALLOC_(G%mask2dBu(IsdB:IedB,JsdB:JedB)) ; G%mask2dBu(:,:) = 0.0
  ALLOC_(G%geoLatT(isd:ied,jsd:jed))      ; G%geoLatT(:,:) = 0.0
  ALLOC_(G%geoLatCu(IsdB:IedB,jsd:jed))   ; G%geoLatCu(:,:) = 0.0
  ALLOC_(G%geoLatCv(isd:ied,JsdB:JedB))   ; G%geoLatCv(:,:) = 0.0
  ALLOC_(G%geoLatBu(IsdB:IedB,JsdB:JedB)) ; G%geoLatBu(:,:) = 0.0
  ALLOC_(G%geoLonT(isd:ied,jsd:jed))      ; G%geoLonT(:,:) = 0.0
  ALLOC_(G%geoLonCu(IsdB:IedB,jsd:jed))   ; G%geoLonCu(:,:) = 0.0
  ALLOC_(G%geoLonCv(isd:ied,JsdB:JedB))   ; G%geoLonCv(:,:) = 0.0
  ALLOC_(G%geoLonBu(IsdB:IedB,JsdB:JedB)) ; G%geoLonBu(:,:) = 0.0

  ALLOC_(G%dx_Cv(isd:ied,JsdB:JedB))     ; G%dx_Cv(:,:) = 0.0
  ALLOC_(G%dy_Cu(IsdB:IedB,jsd:jed))     ; G%dy_Cu(:,:) = 0.0
  ALLOC_(G%dx_Cv_obc(isd:ied,JsdB:JedB)) ; G%dx_Cv_obc(:,:) = 0.0
  ALLOC_(G%dy_Cu_obc(IsdB:IedB,jsd:jed)) ; G%dy_Cu_obc(:,:) = 0.0  

  ALLOC_(G%areaCu(IsdB:IedB,jsd:jed))  ; G%areaCu(:,:) = 0.0
  ALLOC_(G%areaCv(isd:ied,JsdB:JedB))  ; G%areaCv(:,:) = 0.0
  ALLOC_(G%IareaCu(IsdB:IedB,jsd:jed)) ; G%IareaCu(:,:) = 0.0
  ALLOC_(G%IareaCv(isd:ied,JsdB:JedB)) ; G%IareaCv(:,:) = 0.0

  ALLOC_(G%CoriolisBu(IsdB:IedB, JsdB:JedB)) ; G%CoriolisBu(:,:) = 0.0
  ALLOC_(G%dF_dx(isd:ied, jsd:jed)) ; G%dF_dx(:,:) = 0.0
  ALLOC_(G%dF_dy(isd:ied, jsd:jed)) ; G%dF_dy(:,:) = 0.0

end subroutine allocate_metrics

subroutine MOM_grid_end(G)
! Arguments: G - The ocean's grid structure.
  type(ocean_grid_type), intent(inout) :: G

  DEALLOC_(G%dxT)  ; DEALLOC_(G%dxCu)  ; DEALLOC_(G%dxCv)  ; DEALLOC_(G%dxBu)
  DEALLOC_(G%IdxT) ; DEALLOC_(G%IdxCu) ; DEALLOC_(G%IdxCv) ; DEALLOC_(G%IdxBu)

  DEALLOC_(G%dyT)  ; DEALLOC_(G%dyCu)  ; DEALLOC_(G%dyCv)  ; DEALLOC_(G%dyBu)
  DEALLOC_(G%IdyT) ; DEALLOC_(G%IdyCu) ; DEALLOC_(G%IdyCv) ; DEALLOC_(G%IdyBu)

  DEALLOC_(G%areaT)  ; DEALLOC_(G%IareaT)
  DEALLOC_(G%areaBu) ; DEALLOC_(G%IareaBu)
  DEALLOC_(G%areaCu) ; DEALLOC_(G%IareaCu)
  DEALLOC_(G%areaCv)  ; DEALLOC_(G%IareaCv)

  DEALLOC_(G%mask2dT)  ; DEALLOC_(G%mask2dCu)
  DEALLOC_(G%mask2dCv) ; DEALLOC_(G%mask2dBu)

  DEALLOC_(G%geoLatT)  ; DEALLOC_(G%geoLatCu)
  DEALLOC_(G%geoLatCv) ; DEALLOC_(G%geoLatBu)
  DEALLOC_(G%geoLonT)  ; DEALLOC_(G%geoLonCu)
  DEALLOC_(G%geoLonCv) ; DEALLOC_(G%geoLonBu)

  DEALLOC_(G%dx_Cv) ; DEALLOC_(G%dy_Cu)
  DEALLOC_(G%dx_Cv_obc) ; DEALLOC_(G%dy_Cu_obc)

  DEALLOC_(G%bathyT)  ; DEALLOC_(G%CoriolisBu)
  DEALLOC_(G%dF_dx)  ; DEALLOC_(G%dF_dy)
  DEALLOC_(G%g_prime) ; DEALLOC_(G%Rlay)
  deallocate(G%gridLonT) ; deallocate(G%gridLatT)
  deallocate(G%gridLonB) ; deallocate(G%gridLatB)

  call verticalGridEnd( G%GV )
end subroutine MOM_grid_end

end module MOM_grid
