!> Contains a shareable dynamic type for describing horizontal grids and metric data
!! and utilty routines that work on this type.
module MOM_dyn_horgrid

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_array_transform, only : rotate_array, rotate_array_pair
use MOM_domains,         only : MOM_domain_type, deallocate_MOM_domain
use MOM_error_handler,   only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_hor_index,       only : hor_index_type
use MOM_unit_scaling,    only : unit_scale_type

implicit none ; private

public create_dyn_horgrid, destroy_dyn_horgrid, set_derived_dyn_horgrid
public rescale_dyn_horgrid_bathymetry, rotate_dyn_horgrid

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Describes the horizontal ocean grid with only dynamic memory arrays
type, public :: dyn_horgrid_type
  type(MOM_domain_type), pointer :: Domain => NULL() !< Ocean model domain
  type(MOM_domain_type), pointer :: Domain_aux => NULL() !< A non-symmetric auxiliary domain type.
  type(hor_index_type) :: HI !< Horizontal index ranges

  integer :: isc !< The start i-index of cell centers within the computational domain
  integer :: iec !< The end i-index of cell centers within the computational domain
  integer :: jsc !< The start j-index of cell centers within the computational domain
  integer :: jec !< The end j-index of cell centers within the computational domain

  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain

  integer :: isg !< The start i-index of cell centers within the global domain
  integer :: ieg !< The end i-index of cell centers within the global domain
  integer :: jsg !< The start j-index of cell centers within the global domain
  integer :: jeg !< The end j-index of cell centers within the global domain

  integer :: IscB !< The start i-index of cell vertices within the computational domain
  integer :: IecB !< The end i-index of cell vertices within the computational domain
  integer :: JscB !< The start j-index of cell vertices within the computational domain
  integer :: JecB !< The end j-index of cell vertices within the computational domain

  integer :: IsdB !< The start i-index of cell vertices within the data domain
  integer :: IedB !< The end i-index of cell vertices within the data domain
  integer :: JsdB !< The start j-index of cell vertices within the data domain
  integer :: JedB !< The end j-index of cell vertices within the data domain

  integer :: IsgB !< The start i-index of cell vertices within the global domain
  integer :: IegB !< The end i-index of cell vertices within the global domain
  integer :: JsgB !< The start j-index of cell vertices within the global domain
  integer :: JegB !< The end j-index of cell vertices within the global domain

  integer :: isd_global !< The value of isd in the global index space (decompoistion invariant).
  integer :: jsd_global !< The value of isd in the global index space (decompoistion invariant).
  integer :: idg_offset !< The offset between the corresponding global and local i-indices.
  integer :: jdg_offset !< The offset between the corresponding global and local j-indices.
  logical :: symmetric  !< True if symmetric memory is used.

  logical :: nonblocking_updates  !< If true, non-blocking halo updates are
                                  !! allowed.  The default is .false. (for now).
  integer :: first_direction !< An integer that indicates which direction is to be updated first in
                             !! directionally split parts of the calculation.  This can be altered
                             !! during the course of the run via calls to set_first_direction.

  real, allocatable, dimension(:,:) :: &
    mask2dT, &   !< 0 for land points and 1 for ocean points on the h-grid [nondim].
    geoLatT, &   !< The geographic latitude at q points [degrees of latitude] or [m].
    geoLonT, &   !< The geographic longitude at q points [degrees of longitude] or [m].
    dxT, &       !< dxT is delta x at h points [L ~> m].
    IdxT, &      !< 1/dxT [L-1 ~> m-1].
    dyT, &       !< dyT is delta y at h points [L ~> m].
    IdyT, &      !< IdyT is 1/dyT [L-1 ~> m-1].
    areaT, &     !< The area of an h-cell [L2 ~> m2].
    IareaT       !< 1/areaT [L-2 ~> m-2].
  real, allocatable, dimension(:,:) :: sin_rot
                 !< The sine of the angular rotation between the local model grid's northward
                 !! and the true northward directions [nondim].
  real, allocatable, dimension(:,:) :: cos_rot
                 !< The cosine of the angular rotation between the local model grid's northward
                 !! and the true northward directions [nondim].

  real, allocatable, dimension(:,:) :: &
    mask2dCu, &  !< 0 for boundary points and 1 for ocean points on the u grid [nondim].
    OBCmaskCu, & !< 0 for boundary or OBC points and 1 for ocean points on the u grid [nondim].
    geoLatCu, &  !< The geographic latitude at u points [degrees of latitude] or [m].
    geoLonCu, &  !< The geographic longitude at u points [degrees of longitude] or [m].
    dxCu, &      !< dxCu is delta x at u points [L ~> m].
    IdxCu, &     !< 1/dxCu [L-1 ~> m-1].
    dyCu, &      !< dyCu is delta y at u points [L ~> m].
    IdyCu, &     !< 1/dyCu [L-1 ~> m-1].
    dy_Cu, &     !< The unblocked lengths of the u-faces of the h-cell [L ~> m].
    IareaCu, &   !< The masked inverse areas of u-grid cells [L-2 ~> m-2].
    areaCu       !< The areas of the u-grid cells [L2 ~> m2].

  real, allocatable, dimension(:,:) :: &
    mask2dCv, &  !< 0 for boundary points and 1 for ocean points on the v grid [nondim].
    OBCmaskCv, & !< 0 for boundary or OBC points and 1 for ocean points on the v grid [nondim].
    geoLatCv, &  !< The geographic latitude at v points [degrees of latitude] or [m].
    geoLonCv, &  !< The geographic longitude at v points [degrees of longitude] or [m].
    dxCv, &      !< dxCv is delta x at v points [L ~> m].
    IdxCv, &     !< 1/dxCv [L-1 ~> m-1].
    dyCv, &      !< dyCv is delta y at v points [L ~> m].
    IdyCv, &     !< 1/dyCv [L-1 ~> m-1].
    dx_Cv, &     !< The unblocked lengths of the v-faces of the h-cell [L ~> m].
    IareaCv, &   !< The masked inverse areas of v-grid cells [L-2 ~> m-2].
    areaCv       !< The areas of the v-grid cells [L2 ~> m2].

  real, allocatable, dimension(:,:) :: &
    porous_DminU, & !< minimum topographic height (deepest) of U-face [Z ~> m]
    porous_DmaxU, & !< maximum topographic height (shallowest) of U-face [Z ~> m]
    porous_DavgU    !< average topographic height of U-face [Z ~> m]

  real, allocatable, dimension(:,:) :: &
    porous_DminV, & !< minimum topographic height (deepest) of V-face [Z ~> m]
    porous_DmaxV, & !< maximum topographic height (shallowest) of V-face [Z ~> m]
    porous_DavgV    !< average topographic height of V-face [Z ~> m]

  real, allocatable, dimension(:,:) :: &
    mask2dBu, &  !< 0 for boundary points and 1 for ocean points on the q grid [nondim].
    geoLatBu, &  !< The geographic latitude at q points [degrees of latitude] or [m].
    geoLonBu, &  !< The geographic longitude at q points [degrees of longitude] or [m].
    dxBu, &      !< dxBu is delta x at q points [L ~> m].
    IdxBu, &     !< 1/dxBu [L-1 ~> m-1].
    dyBu, &      !< dyBu is delta y at q points [L ~> m].
    IdyBu, &     !< 1/dyBu [L-1 ~> m-1].
    areaBu, &    !< areaBu is the area of a q-cell [L ~> m]
    IareaBu      !< IareaBu = 1/areaBu [L-2 ~> m-2].

  real, pointer, dimension(:) :: gridLatT => NULL()
        !< The latitude of T points for the purpose of labeling the output axes,
        !! often in units of [degrees_N] or [km] or [m] or [gridpoints].
        !! On many grids this is the same as geoLatT.
  real, pointer, dimension(:) :: gridLatB => NULL()
        !< The latitude of B points for the purpose of labeling the output axes,
        !! often in units of [degrees_N] or [km] or [m] or [gridpoints].
        !! On many grids this is the same as geoLatBu.
  real, pointer, dimension(:) :: gridLonT => NULL()
        !< The longitude of T points for the purpose of labeling the output axes,
        !! often in units of [degrees_E] or [km] or [m] or [gridpoints].
        !! On many grids this is the same as geoLonT.
  real, pointer, dimension(:) :: gridLonB => NULL()
        !< The longitude of B points for the purpose of labeling the output axes,
        !! often in units of [degrees_E] or [km] or [m] or [gridpoints].
        !! On many grids this is the same as geoLonBu.
  character(len=40) :: &
    ! Except on a Cartesian grid, these are usually some variant of "degrees".
    x_axis_units, &     !< The units that are used in labeling the x coordinate axes.
    y_axis_units, &     !< The units that are used in labeling the y coordinate axes.
    ! These are internally generated names, including "m", "km", "deg_E" and "deg_N".
    x_ax_unit_short, &  !< A short description of the x-axis units for documenting parameter units
    y_ax_unit_short     !< A short description of the y-axis units for documenting parameter units

  real, allocatable, dimension(:,:) :: &
    bathyT        !< Ocean bottom depth at tracer points, in depth units [Z ~> m].

  logical :: bathymetry_at_vel  !< If true, there are separate values for the
                  !! basin depths at velocity points.  Otherwise the effects of
                  !! of topography are entirely determined from thickness points.
  real, allocatable, dimension(:,:) :: &
    Dblock_u, &   !< Topographic depths at u-points at which the flow is blocked [Z ~> m].
    Dopen_u       !< Topographic depths at u-points at which the flow is open at width dy_Cu [Z ~> m].
  real, allocatable, dimension(:,:) :: &
    Dblock_v, &   !< Topographic depths at v-points at which the flow is blocked [Z ~> m].
    Dopen_v       !< Topographic depths at v-points at which the flow is open at width dx_Cv [Z ~> m].
  real, allocatable, dimension(:,:) :: &
    CoriolisBu, & !< The Coriolis parameter at corner points [T-1 ~> s-1].
    Coriolis2Bu   !< The square of the Coriolis parameter at corner points [T-2 ~> s-2].
  real, allocatable, dimension(:,:) :: &
    df_dx, &      !< Derivative d/dx f (Coriolis parameter) at h-points [T-1 L-1 ~> s-1 m-1].
    df_dy         !< Derivative d/dy f (Coriolis parameter) at h-points [T-1 L-1 ~> s-1 m-1].

  ! These variables are global sums that are useful for 1-d diagnostics and should not be rescaled.
  real :: areaT_global  !< Global sum of h-cell area [m2]
  real :: IareaT_global !< Global sum of inverse h-cell area (1/areaT_global) [m-2]

  ! These parameters are run-time parameters that are used during some
  ! initialization routines (but not all)
  real :: south_lat     !< The latitude (or y-coordinate) of the first v-line [degrees_N] or [km] or [m]
  real :: west_lon      !< The longitude (or x-coordinate) of the first u-line [degrees_E] or [km] or [m]
  real :: len_lat       !< The latitudinal (or y-coord) extent of physical domain [degrees_N] or [km] or [m]
  real :: len_lon       !< The longitudinal (or x-coord) extent of physical domain [degrees_E] or [km] or [m]
  real :: Rad_Earth     !< The radius of the planet [m]
  real :: Rad_Earth_L   !< The radius of the planet in rescaled units [L ~> m]
  real :: max_depth     !< The maximum depth of the ocean [Z ~> m]
end type dyn_horgrid_type

contains

!---------------------------------------------------------------------
!> Allocate memory used by the dyn_horgrid_type and related structures.
subroutine create_dyn_horgrid(G, HI, bathymetry_at_vel)
  type(dyn_horgrid_type), pointer, intent(inout) :: G  !< A pointer to the dynamic horizontal grid type
  type(hor_index_type),   intent(in) :: HI !< A hor_index_type for array extents
  logical,        optional, intent(in) :: bathymetry_at_vel !< If true, there are
                             !! separate values for the basin depths at velocity
                             !! points.  Otherwise the effects of topography are
                             !! entirely determined from thickness points.
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isg, ieg, jsg, jeg

  ! This subroutine allocates the lateral elements of the dyn_horgrid_type that
  ! are always used and zeros them out.

  if (associated(G)) then
    call MOM_error(WARNING, "create_dyn_horgrid called with an associated horgrid_type.")
  else
    allocate(G)
  endif

  G%HI = HI

  G%isc = HI%isc ; G%iec = HI%iec ; G%jsc = HI%jsc ; G%jec = HI%jec
  G%isd = HI%isd ; G%ied = HI%ied ; G%jsd = HI%jsd ; G%jed = HI%jed
  G%isg = HI%isg ; G%ieg = HI%ieg ; G%jsg = HI%jsg ; G%jeg = HI%jeg

  G%IscB = HI%IscB ; G%IecB = HI%IecB ; G%JscB = HI%JscB ; G%JecB = HI%JecB
  G%IsdB = HI%IsdB ; G%IedB = HI%IedB ; G%JsdB = HI%JsdB ; G%JedB = HI%JedB
  G%IsgB = HI%IsgB ; G%IegB = HI%IegB ; G%JsgB = HI%JsgB ; G%JegB = HI%JegB

  G%idg_offset = HI%idg_offset ; G%jdg_offset = HI%jdg_offset
  G%isd_global = G%isd + HI%idg_offset ; G%jsd_global = G%jsd + HI%jdg_offset
  G%symmetric = HI%symmetric

  G%bathymetry_at_vel = .false.
  if (present(bathymetry_at_vel)) G%bathymetry_at_vel = bathymetry_at_vel

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  allocate(G%dxT(isd:ied,jsd:jed), source=0.0)
  allocate(G%dxCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%dxCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%dxBu(IsdB:IedB,JsdB:JedB), source=0.0)
  allocate(G%IdxT(isd:ied,jsd:jed), source=0.0)
  allocate(G%IdxCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%IdxCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%IdxBu(IsdB:IedB,JsdB:JedB), source=0.0)

  allocate(G%dyT(isd:ied,jsd:jed), source=0.0)
  allocate(G%dyCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%dyCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%dyBu(IsdB:IedB,JsdB:JedB), source=0.0)
  allocate(G%IdyT(isd:ied,jsd:jed), source=0.0)
  allocate(G%IdyCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%IdyCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%IdyBu(IsdB:IedB,JsdB:JedB), source=0.0)

  allocate(G%areaT(isd:ied,jsd:jed), source=0.0)
  allocate(G%IareaT(isd:ied,jsd:jed), source=0.0)
  allocate(G%areaBu(IsdB:IedB,JsdB:JedB), source=0.0)
  allocate(G%IareaBu(IsdB:IedB,JsdB:JedB), source=0.0)

  allocate(G%mask2dT(isd:ied,jsd:jed), source=0.0)
  allocate(G%mask2dCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%mask2dCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%mask2dBu(IsdB:IedB,JsdB:JedB), source=0.0)
  allocate(G%OBCmaskCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%OBCmaskCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%geoLatT(isd:ied,jsd:jed), source=0.0)
  allocate(G%geoLatCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%geoLatCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%geoLatBu(IsdB:IedB,JsdB:JedB), source=0.0)
  allocate(G%geoLonT(isd:ied,jsd:jed), source=0.0)
  allocate(G%geoLonCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%geoLonCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%geoLonBu(IsdB:IedB,JsdB:JedB), source=0.0)

  allocate(G%dx_Cv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%dy_Cu(IsdB:IedB,jsd:jed), source=0.0)

  allocate(G%areaCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%areaCv(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%IareaCu(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%IareaCv(isd:ied,JsdB:JedB), source=0.0)

  allocate(G%porous_DminU(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%porous_DmaxU(IsdB:IedB,jsd:jed), source=0.0)
  allocate(G%porous_DavgU(IsdB:IedB,jsd:jed), source=0.0)

  allocate(G%porous_DminV(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%porous_DmaxV(isd:ied,JsdB:JedB), source=0.0)
  allocate(G%porous_DavgV(isd:ied,JsdB:JedB), source=0.0)


  allocate(G%bathyT(isd:ied, jsd:jed), source=0.0)
  allocate(G%CoriolisBu(IsdB:IedB, JsdB:JedB), source=0.0)
  allocate(G%Coriolis2Bu(IsdB:IedB, JsdB:JedB), source=0.0)
  allocate(G%dF_dx(isd:ied, jsd:jed), source=0.0)
  allocate(G%dF_dy(isd:ied, jsd:jed), source=0.0)

  allocate(G%sin_rot(isd:ied,jsd:jed), source=0.0)
  allocate(G%cos_rot(isd:ied,jsd:jed), source=1.0)

  if (G%bathymetry_at_vel) then
    allocate(G%Dblock_u(IsdB:IedB, jsd:jed), source=0.0)
    allocate(G%Dopen_u(IsdB:IedB, jsd:jed), source=0.0)
    allocate(G%Dblock_v(isd:ied, JsdB:JedB), source=0.0)
    allocate(G%Dopen_v(isd:ied, JsdB:JedB), source=0.0)
  endif

  ! gridLonB and gridLatB are used as edge values in some cases, so they
  ! always need to use symmetric memory allcoations.
  allocate(G%gridLonT(isg:ieg), source=0.0)
  allocate(G%gridLonB(isg-1:ieg), source=0.0)
  allocate(G%gridLatT(jsg:jeg), source=0.0)
  allocate(G%gridLatB(jsg-1:jeg), source=0.0)

end subroutine create_dyn_horgrid


!> Copy the rotated contents of one horizontal grid type into another.  The input
!! and output grid type arguments can not use the same object.
subroutine rotate_dyn_horgrid(G_in, G, US, turns)
  type(dyn_horgrid_type), intent(in)    :: G_in   !< The input horizontal grid type
  type(dyn_horgrid_type), intent(inout) :: G      !< An output rotated horizontal grid type
                                                  !! that has already been allocated, but whose
                                                  !! contents are largely replaced here.
  type(unit_scale_type),  intent(in)    :: US     !< A dimensional unit scaling type
  integer, intent(in) :: turns                    !< Number of quarter turns

  ! Center point
  call rotate_array(G_in%geoLonT, turns, G%geoLonT)
  call rotate_array(G_in%geoLatT, turns, G%geoLatT)
  call rotate_array_pair(G_in%dxT, G_in%dyT, turns, G%dxT, G%dyT)
  call rotate_array(G_in%areaT, turns, G%areaT)
  call rotate_array(G_in%bathyT, turns, G%bathyT)

  call rotate_array_pair(G_in%df_dx, G_in%df_dy, turns, G%df_dx, G%df_dy)
  call rotate_array(G_in%sin_rot, turns, G%sin_rot)
  call rotate_array(G_in%cos_rot, turns, G%cos_rot)
  call rotate_array(G_in%mask2dT, turns, G%mask2dT)

  ! Face points
  call rotate_array_pair(G_in%geoLonCu, G_in%geoLonCv, turns, G%geoLonCu, G%geoLonCv)
  call rotate_array_pair(G_in%geoLatCu, G_in%geoLatCv, turns, G%geoLatCu, G%geoLatCv)
  call rotate_array_pair(G_in%dxCu, G_in%dyCv, turns, G%dxCu, G%dyCv)
  call rotate_array_pair(G_in%dxCv, G_in%dyCu, turns, G%dxCv, G%dyCu)
  call rotate_array_pair(G_in%dx_Cv, G_in%dy_Cu, turns, G%dx_Cv, G%dy_Cu)

  call rotate_array_pair(G_in%mask2dCu, G_in%mask2dCv, turns, G%mask2dCu, G%mask2dCv)
  call rotate_array_pair(G_in%OBCmaskCu, G_in%OBCmaskCv, turns, G%OBCmaskCu, G%OBCmaskCv)
  call rotate_array_pair(G_in%areaCu, G_in%areaCv, turns, G%areaCu, G%areaCv)
  call rotate_array_pair(G_in%IareaCu, G_in%IareaCv, turns, G%IareaCu, G%IareaCv)

  call rotate_array_pair(G_in%porous_DminU, G_in%porous_DminV, &
       turns, G%porous_DminU, G%porous_DminV)
  call rotate_array_pair(G_in%porous_DmaxU, G_in%porous_DmaxV, &
       turns, G%porous_DmaxU, G%porous_DmaxV)
  call rotate_array_pair(G_in%porous_DavgU, G_in%porous_DavgV, &
       turns, G%porous_DavgU, G%porous_DavgV)


  ! Vertex point
  call rotate_array(G_in%geoLonBu, turns, G%geoLonBu)
  call rotate_array(G_in%geoLatBu, turns, G%geoLatBu)
  call rotate_array_pair(G_in%dxBu, G_in%dyBu, turns, G%dxBu, G%dyBu)
  call rotate_array(G_in%areaBu, turns, G%areaBu)
  call rotate_array(G_in%CoriolisBu, turns, G%CoriolisBu)
  call rotate_array(G_in%Coriolis2Bu, turns, G%Coriolis2Bu)
  call rotate_array(G_in%mask2dBu, turns, G%mask2dBu)

  ! Topography at the cell faces
  G%bathymetry_at_vel = G_in%bathymetry_at_vel
  if (G%bathymetry_at_vel) then
    call rotate_array_pair(G_in%Dblock_u, G_in%Dblock_v, turns, G%Dblock_u, G%Dblock_v)
    call rotate_array_pair(G_in%Dopen_u, G_in%Dopen_v, turns, G%Dopen_u, G%Dopen_v)
  endif

  ! Nominal grid axes
  ! TODO: We should not assign lat values to the lon axis, and vice versa.
  !   We temporarily copy lat <-> lon since several components still expect
  !   lat and lon sizes to match the first and second dimension sizes.
  !   But we ought to instead leave them unchanged and adjust the references to
  !   these axes.
  if (modulo(turns, 2) /= 0) then
    G%gridLonT(:) = G_in%gridLatT(G_in%jeg:G_in%jsg:-1)
    G%gridLatT(:) = G_in%gridLonT(:)
    G%gridLonB(:) = G_in%gridLatB(G_in%jeg:(G_in%jsg-1):-1)
    G%gridLatB(:) = G_in%gridLonB(:)
  else
    G%gridLonT(:) = G_in%gridLonT(:)
    G%gridLatT(:) = G_in%gridLatT(:)
    G%gridLonB(:) = G_in%gridLonB(:)
    G%gridLatB(:) = G_in%gridLatB(:)
  endif

  G%x_axis_units = G_in%y_axis_units
  G%y_axis_units = G_in%x_axis_units
  G%x_ax_unit_short = G_in%y_ax_unit_short
  G%y_ax_unit_short = G_in%x_ax_unit_short
  G%south_lat = G_in%south_lat
  G%west_lon = G_in%west_lon
  G%len_lat = G_in%len_lat
  G%len_lon = G_in%len_lon

  ! Rotation-invariant fields
  G%areaT_global = G_in%areaT_global
  G%IareaT_global = G_in%IareaT_global
  G%Rad_Earth = G_in%Rad_Earth
  G%Rad_Earth_L = G_in%Rad_Earth_L
  G%max_depth = G_in%max_depth

  call set_derived_dyn_horgrid(G, US)
end subroutine rotate_dyn_horgrid


!> rescale_dyn_horgrid_bathymetry permits a change in the internal units for the bathymetry on the
!! grid, both rescaling the depths and recording the new internal depth units.
subroutine rescale_dyn_horgrid_bathymetry(G, m_in_new_units)
  type(dyn_horgrid_type), intent(inout) :: G !< The dynamic horizontal grid type
  real,                   intent(in)    :: m_in_new_units !< The new internal representation of 1 m depth [m Z-1 ~> 1]

  ! Local variables
  real :: rescale ! The inverse of m_in_new_units, used in rescaling bathymetry [Z m-1 ~> 1]
  integer :: i, j, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (m_in_new_units == 1.0) return
  if (m_in_new_units < 0.0) &
    call MOM_error(FATAL, "rescale_dyn_horgrid_bathymetry: Negative depth units are not permitted.")
  if (m_in_new_units == 0.0) &
    call MOM_error(FATAL, "rescale_dyn_horgrid_bathymetry: Zero depth units are not permitted.")

  rescale = 1.0 / m_in_new_units
  do j=jsd,jed ; do i=isd,ied
    G%bathyT(i,j) = rescale*G%bathyT(i,j)
  enddo ; enddo
  if (G%bathymetry_at_vel) then ; do j=jsd,jed ; do I=IsdB,IedB
    G%Dblock_u(I,j) = rescale*G%Dblock_u(I,j) ; G%Dopen_u(I,j) = rescale*G%Dopen_u(I,j)
  enddo ; enddo ; endif
  if (G%bathymetry_at_vel) then ; do J=JsdB,JedB ; do i=isd,ied
    G%Dblock_v(i,J) = rescale*G%Dblock_v(i,J) ; G%Dopen_v(i,J) = rescale*G%Dopen_v(i,J)
  enddo ; enddo ; endif
  G%max_depth = rescale*G%max_depth

end subroutine rescale_dyn_horgrid_bathymetry

!> set_derived_dyn_horgrid calculates metric terms that are derived from other metrics.
subroutine set_derived_dyn_horgrid(G, US)
  type(dyn_horgrid_type), intent(inout) :: G !< The dynamic horizontal grid type
  type(unit_scale_type), optional, intent(in) :: US !< A dimensional unit scaling type
!    Various inverse grid spacings and derived areas are calculated within this
!  subroutine.
  integer :: i, j, isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  do j=jsd,jed ; do i=isd,ied
    if (G%dxT(i,j) < 0.0) G%dxT(i,j) = 0.0
    if (G%dyT(i,j) < 0.0) G%dyT(i,j) = 0.0
    G%IdxT(i,j) = Adcroft_reciprocal(G%dxT(i,j))
    G%IdyT(i,j) = Adcroft_reciprocal(G%dyT(i,j))
    G%IareaT(i,j) = Adcroft_reciprocal(G%areaT(i,j))
  enddo ; enddo

  do j=jsd,jed ; do I=IsdB,IedB
    if (G%dxCu(I,j) < 0.0) G%dxCu(I,j) = 0.0
    if (G%dyCu(I,j) < 0.0) G%dyCu(I,j) = 0.0
    G%IdxCu(I,j) = Adcroft_reciprocal(G%dxCu(I,j))
    G%IdyCu(I,j) = Adcroft_reciprocal(G%dyCu(I,j))
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    if (G%dxCv(i,J) < 0.0) G%dxCv(i,J) = 0.0
    if (G%dyCv(i,J) < 0.0) G%dyCv(i,J) = 0.0
    G%IdxCv(i,J) = Adcroft_reciprocal(G%dxCv(i,J))
    G%IdyCv(i,J) = Adcroft_reciprocal(G%dyCv(i,J))
  enddo ; enddo

  do J=JsdB,JedB ; do I=IsdB,IedB
    if (G%dxBu(I,J) < 0.0) G%dxBu(I,J) = 0.0
    if (G%dyBu(I,J) < 0.0) G%dyBu(I,J) = 0.0

    G%IdxBu(I,J) = Adcroft_reciprocal(G%dxBu(I,J))
    G%IdyBu(I,J) = Adcroft_reciprocal(G%dyBu(I,J))
    ! areaBu has usually been set to a positive area elsewhere.
    if (G%areaBu(I,J) <= 0.0) G%areaBu(I,J) = G%dxBu(I,J) * G%dyBu(I,J)
    G%IareaBu(I,J) =  Adcroft_reciprocal(G%areaBu(I,J))
  enddo ; enddo

end subroutine set_derived_dyn_horgrid

!> Adcroft_reciprocal(x) = 1/x for |x|>0 or 0 for x=0.
function Adcroft_reciprocal(val) result(I_val)
  real, intent(in) :: val  !< The value being inverted in abitrary units [A ~> a]
  real :: I_val            !< The Adcroft reciprocal of val [A-1 ~> a-1].

  I_val = 0.0 ; if (val /= 0.0) I_val = 1.0/val
end function Adcroft_reciprocal

!---------------------------------------------------------------------
!> Release memory used by the dyn_horgrid_type and related structures.
subroutine destroy_dyn_horgrid(G)
  type(dyn_horgrid_type), pointer :: G !< The dynamic horizontal grid type

  if (.not.associated(G)) then
    call MOM_error(FATAL, "destroy_dyn_horgrid called with an unassociated horgrid_type.")
  endif

  deallocate(G%dxT)  ; deallocate(G%dxCu)  ; deallocate(G%dxCv)  ; deallocate(G%dxBu)
  deallocate(G%IdxT) ; deallocate(G%IdxCu) ; deallocate(G%IdxCv) ; deallocate(G%IdxBu)

  deallocate(G%dyT)  ; deallocate(G%dyCu)  ; deallocate(G%dyCv)  ; deallocate(G%dyBu)
  deallocate(G%IdyT) ; deallocate(G%IdyCu) ; deallocate(G%IdyCv) ; deallocate(G%IdyBu)

  deallocate(G%areaT)  ; deallocate(G%IareaT)
  deallocate(G%areaBu) ; deallocate(G%IareaBu)
  deallocate(G%areaCu) ; deallocate(G%IareaCu)
  deallocate(G%areaCv)  ; deallocate(G%IareaCv)

  deallocate(G%mask2dT)  ; deallocate(G%mask2dCu) ; deallocate(G%OBCmaskCu)
  deallocate(G%mask2dCv) ; deallocate(G%OBCmaskCv) ; deallocate(G%mask2dBu)

  deallocate(G%geoLatT)  ; deallocate(G%geoLatCu)
  deallocate(G%geoLatCv) ; deallocate(G%geoLatBu)
  deallocate(G%geoLonT)  ; deallocate(G%geoLonCu)
  deallocate(G%geoLonCv) ; deallocate(G%geoLonBu)

  deallocate(G%dx_Cv) ; deallocate(G%dy_Cu)

  deallocate(G%porous_DminU) ; deallocate(G%porous_DmaxU) ; deallocate(G%porous_DavgU)
  deallocate(G%porous_DminV) ; deallocate(G%porous_DmaxV) ; deallocate(G%porous_DavgV)

  deallocate(G%bathyT)  ; deallocate(G%CoriolisBu) ; deallocate(G%Coriolis2Bu)
  deallocate(G%dF_dx)   ; deallocate(G%dF_dy)
  deallocate(G%sin_rot) ; deallocate(G%cos_rot)

  if (allocated(G%Dblock_u)) deallocate(G%Dblock_u)
  if (allocated(G%Dopen_u)) deallocate(G%Dopen_u)
  if (allocated(G%Dblock_v)) deallocate(G%Dblock_v)
  if (allocated(G%Dopen_v)) deallocate(G%Dopen_v)

  deallocate(G%gridLonT) ; deallocate(G%gridLatT)
  deallocate(G%gridLonB) ; deallocate(G%gridLatB)

  ! CS%debug is required to validate Domain_aux, so use allocation test
  if (associated(G%Domain_aux)) call deallocate_MOM_domain(G%Domain_aux)

  call deallocate_MOM_domain(G%Domain)

  deallocate(G)

end subroutine destroy_dyn_horgrid

end module MOM_dyn_horgrid
