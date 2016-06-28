! This file is part of MOM6. See LICENSE.md for the license.
!> Controls where open boundary conditions are applied 
module MOM_open_boundary

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_domains, only : pass_var, pass_vector
use MOM_domains, only : To_All, SCALAR_PAIR, CGRID_NE
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type, log_param
use MOM_grid, only : ocean_grid_type
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_io, only : EAST_FACE, NORTH_FACE
use MOM_io, only : slasher, read_data
use MOM_tracer_registry, only : add_tracer_OBC_values, tracer_registry_type
use MOM_variables, only : thermo_var_ptrs

implicit none ; private

#include <MOM_memory.h>

public open_boundary_config
public open_boundary_init
public open_boundary_query
public open_boundary_end
public open_boundary_impose_normal_slope
public open_boundary_impose_land_mask
public Radiation_Open_Bdry_Conds
public set_Flather_positions
public set_Flather_data

integer, parameter, public :: OBC_NONE = 0, OBC_SIMPLE = 1, OBC_WALL = 2
integer, parameter, public :: OBC_FLATHER = 3
integer, parameter, public :: OBC_DIRECTION_N = 100 !< Indicates the boundary is an effective northern boundary
integer, parameter, public :: OBC_DIRECTION_S = 200 !< Indicates the boundary is an effective southern boundary
integer, parameter, public :: OBC_DIRECTION_E = 300 !< Indicates the boundary is an effective eastern boundary
integer, parameter, public :: OBC_DIRECTION_W = 400 !< Indicates the boundary is an effective western boundary

!> Open-boundary data
type, public :: ocean_OBC_type
  logical :: apply_OBC_u_flather_east = .false.  !< True if any zonal velocity points in the
                                                 !! local domain use east-facing Flather OBCs.
  logical :: apply_OBC_u_flather_west = .false.  !< True if any zonal velocity points in the
                                                 !! local domain use west-facing Flather OBCs.
  logical :: apply_OBC_v_flather_north = .false. !< True if any zonal velocity points in the
                                                 !! local domain use north-facing Flather OBCs.
  logical :: apply_OBC_v_flather_south = .false. !< True if any zonal velocity points in the
                                                 !! local domain use south-facing Flather OBCs.
  logical :: apply_OBC_u = .false.  !< True if any zonal velocity points in to local domain use OBCs.
  logical :: apply_OBC_v = .false.  !< True if any meridional velocity points in to local domain use OBCs.
  logical, pointer, dimension(:,:) :: &
    OBC_mask_u => NULL(), & !< True at zonal velocity points that have prescribed OBCs.
    OBC_mask_v => NULL()    !< True at meridional velocity points that have prescribed OBCs.
  ! These arrays indicate the kind of open boundary conditions that are to be applied at the u and v
  ! points, and can be OBC_NONE, OBC_SIMPLE, OBC_WALL, or OBC_FLATHER.  Generally these
  ! should be consistent with OBC_mask_[uv], with OBC_mask_[uv] .false. for OBC_kind_[uv] = NONE
  ! and true for all other values.
  integer, pointer, dimension(:,:) :: &
    OBC_kind_u => NULL(), & !< Type of OBC at u-points.
    OBC_kind_v => NULL()    !< Type of OBC at v-points.
  ! These arrays indicate the outward-pointing orientation of the open boundary and will be set to
  ! one of OBC_DIRECTION_N, OBC_DIRECTION_S, OBC_DIRECTION_E or OBC_DIRECTION_W.
  integer, pointer, dimension(:,:) :: &
    OBC_direction_u => NULL(), & !< Orientation of OBC at u-points.
    OBC_direction_v => NULL()    !< Orientation of OBC at v-points.
  ! The following apply at points with OBC_kind_[uv] = OBC_FLATHER.
  real, pointer, dimension(:,:,:) :: &
    rx_old_u => NULL(), &  !< The rx_old_u value for radiation coeff for u-velocity in x-direction
    ry_old_v => NULL(), &  !< The ry_old_v value for radiation coeff for v-velocity in y-direction
    rx_old_h => NULL(), &  !< The rx_old_h value for radiation coeff for layer thickness h in x-direction
    ry_old_h => NULL()     !< The ry_old_h value for radiation coeff for layer thickness h in y-direction

  !   The following can be used to specify the outer-domain values of the
  ! surface height and barotropic velocity.  If these are not allocated, the
  ! default with Flather boundary conditions is the same as if they were
  ! filled with zeros.  With simple OBCs, these should not be allocated.
  real, pointer, dimension(:,:) :: &
    ubt_outer => NULL(), &    !< The u-velocity in the outer domain, in m s-1.
    vbt_outer => NULL(), &    !< The v-velocity in the outer domain, in m s-1.
    eta_outer_u => NULL(), &  !< The SSH anomaly in the outer domain, in m or kg m-2.
    eta_outer_v => NULL()     !< The SSH anomaly in the outer domain, in m or kg m-2.

  ! The following apply at points with OBC_kind_[uv] = OBC_SIMPLE.
  real, pointer, dimension(:,:,:) :: &
    u => NULL(), &  !< The prescribed values of the zonal velocity (u) at OBC points.
    v => NULL(), &  !< The prescribed values of the meridional velocity (v) at OBC points.
    uh => NULL(), & !< The prescribed values of the zonal volume transport (uh) at OBC points.
    vh => NULL()    !< The prescribed values of the meridional volume transport (vh) at OBC points.

  ! The following parameters are used in the baroclinic radiation code:
  real :: gamma_uv !< The relative weighting for the baroclinic radiation
                   !! velocities (or speed of characteristics) at the
                   !! new time level (1) or the running mean (0) for velocities.
                   !! Valid values range from 0 to 1, with a default of 0.3.
  real :: gamma_h  !< The relative weighting for the baroclinic radiation
                   !! velocities (or speed of characteristics) at the
                   !! new time level (1) or the running mean (0) for thicknesses.
                   !! Valid values range from 0 to 1, with a default of 0.2.
  real :: rx_max   !< The maximum magnitude of the baroclinic radiation
                   !! velocity (or speed of characteristics), in m s-1.  The
                   !! default value is 10 m s-1.
end type ocean_OBC_type

integer :: id_clock_pass

character(len=40)  :: mod = "MOM_open_boundary" ! This module's name.
! This include declares and sets the variable "version".
#include "version_variable.h"

contains

!> Enables OBC module and reads configuration parameters
subroutine open_boundary_config(G, param_file, OBC)
  type(dyn_horgrid_type),  intent(in)    :: G !< Ocean grid structure
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handle
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary control structure
  ! Local variables
  logical :: flather_east, flather_west, flather_north, flather_south

  allocate(OBC)

  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "APPLY_OBC_U", OBC%apply_OBC_u, &
                 "If true, open boundary conditions may be set at some \n"//&
                 "u-points, with the configuration controlled by OBC_CONFIG", &
                 default=.false.)
  call get_param(param_file, mod, "APPLY_OBC_V", OBC%apply_OBC_v, &
                 "If true, open boundary conditions may be set at some \n"//&
                 "v-points, with the configuration controlled by OBC_CONFIG", &
                 default=.false.)
  call get_param(param_file, mod, "APPLY_OBC_U_FLATHER_EAST", OBC%apply_OBC_u_flather_east, &
                 "Apply a Flather open boundary condition on the eastern\n"//&
                 "side of the global domain", &
                 default=.false.)
  call get_param(param_file, mod, "APPLY_OBC_U_FLATHER_WEST", OBC%apply_OBC_u_flather_west, &
                 "Apply a Flather open boundary condition on the western\n"//&
                 "side of the global domain", &
                 default=.false.)
  call get_param(param_file, mod, "APPLY_OBC_V_FLATHER_NORTH", OBC%apply_OBC_v_flather_north, &
                 "Apply a Flather open boundary condition on the northern\n"//&
                 "side of the global domain", &
                 default=.false.)
  call get_param(param_file, mod, "APPLY_OBC_V_FLATHER_SOUTH", OBC%apply_OBC_v_flather_south, &
                 "Apply a Flather open boundary condition on the southern\n"//&
                 "side of the global domain", &
                 default=.false.)

  ! Safety check
  if ((OBC%apply_OBC_u_flather_west .or. OBC%apply_OBC_v_flather_south) .and. &
      .not.G%symmetric ) call MOM_error(FATAL, &
                 "MOM_open_boundary, open_boundary_config: "//&
                 "Symmetric memory must be used when APPLY_OBC_U_FLATHER_WEST "//&
                 "or APPLY_OBC_V_FLATHER_SOUTH is true.")

  if (.not.(OBC%apply_OBC_u .or. OBC%apply_OBC_v .or. &
            OBC%apply_OBC_v_flather_north .or. OBC%apply_OBC_v_flather_south .or. &
            OBC%apply_OBC_u_flather_east .or. OBC%apply_OBC_u_flather_west)) then
    ! No open boundaries have been requested
    call open_boundary_dealloc(OBC)
  endif

end subroutine open_boundary_config

!> Initialize open boundary control structure
subroutine open_boundary_init(G, param_file, OBC)
  type(ocean_grid_type), intent(in)    :: G !< Ocean grid structure
  type(param_file_type), intent(in)    :: param_file !< Parameter file handle
  type(ocean_OBC_type),  pointer       :: OBC !< Open boundary control structure
  ! Local variables

  if (.not.associated(OBC)) return

  if ( OBC%apply_OBC_v_flather_north .or. OBC%apply_OBC_v_flather_south .or. &
       OBC%apply_OBC_u_flather_east .or. OBC%apply_OBC_u_flather_west ) then
    call get_param(param_file, mod, "OBC_RADIATION_MAX", OBC%rx_max, &
                   "The maximum magnitude of the baroclinic radiation \n"//&
                   "velocity (or speed of characteristics).  This is only \n"//&
                   "used if one of the APPLY_OBC_[UV]_FLATHER_... is true.", &
                   units="m s-1", default=10.0)
    call get_param(param_file, mod, "OBC_RAD_VEL_WT", OBC%gamma_uv, &
                   "The relative weighting for the baroclinic radiation \n"//&
                   "velocities (or speed of characteristics) at the new \n"//&
                   "time level (1) or the running mean (0) for velocities. \n"//&
                   "Valid values range from 0 to 1. This is only used if \n"//&
                   "one of the APPLY_OBC_[UV]_FLATHER_...  is true.", &
                   units="nondim",  default=0.3)
    call get_param(param_file, mod, "OBC_RAD_THICK_WT", OBC%gamma_h, &
                   "The relative weighting for the baroclinic radiation \n"//&
                   "velocities (or speed of characteristics) at the new \n"//&
                   "time level (1) or the running mean (0) for thicknesses. \n"//&
                   "Valid values range from 0 to 1. This is only used if \n"//&
                   "one of the APPLY_OBC_[UV]_FLATHER_...  is true.", &
                   units="nondim",  default=0.2)
  endif

  id_clock_pass = cpu_clock_id('(Ocean OBC halo updates)', grain=CLOCK_ROUTINE)

end subroutine open_boundary_init

!> Query the state of open boundary module configuration
logical function open_boundary_query(OBC, apply_orig_OBCs, apply_orig_Flather)
  type(ocean_OBC_type), pointer     :: OBC !< Open boundary control structure
  logical, optional,    intent(in)  :: apply_orig_OBCs !< If present, returns True if APPLY_OBC_U/V was set
  logical, optional,    intent(in)  :: apply_orig_Flather !< If present, returns True if APPLY_OBC_*_FLATHER_* was set
  open_boundary_query = .false.
  if (.not. associated(OBC)) return
  if (present(apply_orig_OBCs)) open_boundary_query = OBC%apply_OBC_u .or. OBC%apply_OBC_v
  if (present(apply_orig_Flather)) open_boundary_query = OBC%apply_OBC_v_flather_north .or. &
                                                         OBC%apply_OBC_v_flather_south .or. &
                                                         OBC%apply_OBC_u_flather_east .or. &
                                                         OBC%apply_OBC_u_flather_west
end function open_boundary_query

!> Deallocate open boundary data
subroutine open_boundary_dealloc(OBC)
  type(ocean_OBC_type), pointer :: OBC !< Open boundary control structure
  if (.not. associated(OBC)) return
  if (associated(OBC%OBC_mask_u)) deallocate(OBC%OBC_mask_u)
  if (associated(OBC%OBC_mask_v)) deallocate(OBC%OBC_mask_v)
  if (associated(OBC%OBC_kind_u)) deallocate(OBC%OBC_kind_u)
  if (associated(OBC%OBC_kind_v)) deallocate(OBC%OBC_kind_v)
  if (associated(OBC%rx_old_u)) deallocate(OBC%rx_old_u)
  if (associated(OBC%ry_old_v)) deallocate(OBC%ry_old_v)
  if (associated(OBC%rx_old_h)) deallocate(OBC%rx_old_h)
  if (associated(OBC%ry_old_h)) deallocate(OBC%ry_old_h)
  if (associated(OBC%ubt_outer)) deallocate(OBC%ubt_outer)
  if (associated(OBC%vbt_outer)) deallocate(OBC%vbt_outer)
  if (associated(OBC%eta_outer_u)) deallocate(OBC%eta_outer_u)
  if (associated(OBC%eta_outer_v)) deallocate(OBC%eta_outer_v)
  if (associated(OBC%u)) deallocate(OBC%u)
  if (associated(OBC%v)) deallocate(OBC%v)
  if (associated(OBC%uh)) deallocate(OBC%uh)
  if (associated(OBC%vh)) deallocate(OBC%vh)
  deallocate(OBC)
end subroutine open_boundary_dealloc

!> Close open boundary data
subroutine open_boundary_end(OBC)
  type(ocean_OBC_type), pointer :: OBC !< Open boundary control structure
  call open_boundary_dealloc(OBC)
end subroutine open_boundary_end

!> Sets the slope of bathymetry normal to an open bounndary to zero.
subroutine open_boundary_impose_normal_slope(OBC, G, depth)
  type(ocean_OBC_type),             pointer       :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),           intent(in)    :: G !< Ocean grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: depth !< Bathymetry at h-points
  ! Local variables
  integer :: i, j

  if (.not.associated(OBC)) return

  if (associated(OBC%OBC_direction_u)) then
    do j=G%jsd,G%jed ; do I=G%isd,G%ied-1
      if (OBC%OBC_direction_u(I,j) == OBC_DIRECTION_E) depth(i+1,j) = depth(i,j)
      if (OBC%OBC_direction_u(I,j) == OBC_DIRECTION_W) depth(i,j) = depth(i+1,j)
    enddo ; enddo
  endif

  if (associated(OBC%OBC_kind_v)) then
    do J=G%jsd,G%jed-1 ; do i=G%isd,G%ied
      if (OBC%OBC_direction_v(i,J) == OBC_DIRECTION_N) depth(i,j+1) = depth(i,j)
      if (OBC%OBC_direction_v(i,J) == OBC_DIRECTION_S) depth(i,j) = depth(i,j+1)
    enddo ; enddo
  endif

end subroutine open_boundary_impose_normal_slope

!> Reconcile masks and open boundaries, deallocate OBC on PEs where it is not needed
subroutine open_boundary_impose_land_mask(OBC, G)
  type(ocean_OBC_type),              pointer       :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),            intent(in) :: G !< Ocean grid structure
  ! Local variables
  integer :: i, j
  logical :: any_U, any_V

  if (.not.associated(OBC)) return

  if (associated(OBC%OBC_kind_u)) then
    do j=G%jsd,G%jed ; do I=G%isd,G%ied-1
      if (G%mask2dCu(I,j) == 0) then
        OBC%OBC_kind_u(I,j) = OBC_NONE
        OBC%OBC_direction_u(I,j) = OBC_NONE
        OBC%OBC_mask_u(I,j) = .false.
      endif
    enddo ; enddo
  endif

  if (associated(OBC%OBC_kind_v)) then
    do J=G%jsd,G%jed-1 ; do i=G%isd,G%ied
      if (G%mask2dCv(i,J) == 0) then
        OBC%OBC_kind_v(i,J) = OBC_NONE
        OBC%OBC_direction_v(i,J) = OBC_NONE
        OBC%OBC_mask_v(i,J) = .false.
      endif
    enddo ; enddo
  endif

  any_U = .false.
  if (associated(OBC%OBC_mask_u)) then
    do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
      ! G%mask2du will be open wherever bathymetry allows it.
      ! Bathymetry outside of the open boundary was adjusted to match
      ! the bathymetry inside so these points will be open unless the
      ! bathymetry inside the boundary was do shallow and flagged as land.
      if (OBC%OBC_mask_u(I,j)) any_U = .true.
    enddo ; enddo
    if (.not. any_U) then
      deallocate(OBC%OBC_mask_u)
    endif
  endif

  any_V = .false.
  if (associated(OBC%OBC_mask_v)) then
    do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
      if (OBC%OBC_mask_v(i,J)) any_V = .true.
    enddo ; enddo
    if (.not. any_V) then
      deallocate(OBC%OBC_mask_v)
    endif
  endif

  if (.not.(any_U .or. any_V)) call open_boundary_dealloc(OBC)

end subroutine open_boundary_impose_land_mask

!> Diagnose radiation conditions at open boundaries
subroutine Radiation_Open_Bdry_Conds(OBC, u_new, u_old, v_new, v_old, &
                                     h_new, h_old, G)
  type(ocean_grid_type),                     intent(inout) :: G !< Ocean grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary control structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u_new
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u_old
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v_new
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v_old
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h_new
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h_old
  ! Local variables
  real :: dhdt, dhdx, gamma_u, gamma_h, gamma_v
  real :: rx_max, ry_max ! coefficients for radiation
  real :: rx_new, rx_avg ! coefficients for radiation

  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(OBC)) return
  if (.not.(OBC%apply_OBC_u_flather_east .or. OBC%apply_OBC_u_flather_west .or. &
            OBC%apply_OBC_v_flather_north .or. OBC%apply_OBC_v_flather_south)) &
    return

  gamma_u = OBC%gamma_uv ; gamma_v = OBC%gamma_uv ; gamma_h = OBC%gamma_h
  rx_max = OBC%rx_max ; ry_max = OBC%rx_max

  if (OBC%apply_OBC_u_flather_east .or. OBC%apply_OBC_u_flather_west) then
    do k=1,nz ; do j=js,je ; do I=is-1,ie ; if (OBC%OBC_mask_u(I,j)) then
      if (OBC%OBC_direction_u(I,j) == OBC_DIRECTION_E) then
        dhdt = u_old(I-1,j,k)-u_new(I-1,j,k) !old-new
        dhdx = u_new(I-1,j,k)-u_new(I-2,j,k) !in new time backward sasha for I-1
        rx_new = 0.0
        if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
        rx_avg = (1.0-gamma_u)*OBC%rx_old_u(I,j,k) + gamma_u*rx_new
        OBC%rx_old_u(I,j,k) = rx_avg
        u_new(I,j,k) = (u_old(I,j,k) + rx_avg*u_new(I-1,j,k)) / (1.0+rx_avg)

    !   dhdt = h_old(I,j,k)-h_new(I,j,k) !old-new
    !   dhdx = h_new(I,j,k)-h_new(I-1,j,k) !in new time
    !   rx_new = 0.0
    !   if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
    !   rx_avg = (1.0-gamma_h)*OBC%rx_old_h(I,j,k) + gamma_h*rx_new
    !   OBC%rx_old_h(I,j,k) = rx_avg
    !    h_new(I+1,j,k) = (h_old(I+1,j,k) + rx_avg*h_new(I,j,k)) / (1.0+rx_avg) !original
      endif
      if (OBC%OBC_direction_u(I,j) == OBC_DIRECTION_W) then
        dhdt = u_old(I+1,j,k)-u_new(I+1,j,k) !old-new
        dhdx = u_new(I+1,j,k)-u_new(I+2,j,k) !in new time backward sasha for I+1
        rx_new = 0.0
        if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
        rx_avg = (1.0-gamma_u)*OBC%rx_old_u(I,j,k) + gamma_u*rx_new
        OBC%rx_old_u(I,j,k) = rx_avg
        u_new(I,j,k) = (u_old(I,j,k) + rx_avg*u_new(I+1,j,k)) / (1.0+rx_avg)

    !   dhdt = h_old(I+1,j,k)-h_new(I+1,j,k) !old-new
    !   dhdx = h_new(I+1,j,k)-h_new(I+2,j,k) !in new time
    !   rx_new = 0.0
    !   if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
    !   rx_avg = (1.0-gamma_h)*OBC%rx_old_h(I,j,k) + gamma_h*rx_new
    !   OBC%rx_old_h(I,j,k) = rx_avg
    !   h_new(I,j,k) = (h_old(I,j,k) + rx_avg*h_new(I+1,j,k)) / (1.0+rx_avg) !original
      endif
    endif ; enddo ; enddo ; enddo
  endif

  if (OBC%apply_OBC_v_flather_north .or. OBC%apply_OBC_v_flather_south) then
    do k=1,nz ; do J=js-1,je ; do i=is,ie ; if (OBC%OBC_mask_v(i,J)) then
      if (OBC%OBC_direction_v(i,J) == OBC_DIRECTION_N) then
        dhdt = v_old(i,J-1,k)-v_new(i,J-1,k) !old-new
        dhdx = v_new(i,J-1,k)-v_new(i,J-2,k) !in new time backward sasha for J-1
        rx_new = 0.0
        if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
        rx_avg = (1.0-gamma_v)*OBC%ry_old_v(i,J,k) + gamma_v*rx_new
        OBC%ry_old_v(i,J,k) = rx_avg
        v_new(i,J,k) = (v_old(i,J,k) + rx_avg*v_new(i,J-1,k)) / (1.0+rx_avg)

    !   dhdt = h_old(i,J,k)-h_new(i,J,k) !old-new
    !   dhdx = h_new(i,J,k)-h_new(i,J-1,k) !in new time
    !   rx_new = 0.0
    !   if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
    !   rx_avg = (1.0-gamma_h)*OBC%ry_old_h(i,J,k) + gamma_h*rx_new
    !   OBC%ry_old_h(i,J,k) = rx_avg
    !   h_new(i,J+1,k) = (h_old(i,J+1,k) + rx_avg*h_new(i,J,k)) / (1.0+rx_avg) !original
      endif

      if (OBC%OBC_direction_v(i,J) == OBC_DIRECTION_S) then
        dhdt = v_old(i,J+1,k)-v_new(i,J+1,k) !old-new
        dhdx = v_new(i,J+1,k)-v_new(i,J+2,k) !in new time backward sasha for J+1
        rx_new = 0.0
        if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
        rx_avg = (1.0-gamma_v)*OBC%ry_old_v(i,J,k) + gamma_v*rx_new
        OBC%ry_old_v(i,J,k) = rx_avg
        v_new(i,J,k) = (v_old(i,J,k) + rx_avg*v_new(i,J+1,k)) / (1.0+rx_avg)

    !   dhdt = h_old(i,J+1,k)-h_new(i,J+1,k) !old-new
    !   dhdx = h_new(i,J+1,k)-h_new(i,J+2,k) !in new time
    !   rx_new = 0.0
    !   if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
    !   rx_avg = (1.0-gamma_h)*OBC%ry_old_h(i,J,k) + gamma_h*rx_new
    !   OBC%ry_old_h(i,J,k) = rx_avg
    !   h_new(i,J,k) = (h_old(i,J,k) + rx_avg*h_new(i,J+1,k)) / (1.0+rx_avg) !original
      endif

    endif ; enddo ; enddo ; enddo
  endif

  call cpu_clock_begin(id_clock_pass)
  call pass_vector(u_new, v_new, G%Domain)
  call pass_var(h_new, G%Domain)
  call cpu_clock_end(id_clock_pass)

end subroutine Radiation_Open_Bdry_Conds

!> Sets the domain boundaries as Flather open boundaries using the original
!! Flather run-time logicals
subroutine set_Flather_positions(G, OBC)
  type(dyn_horgrid_type),                 intent(inout) :: G
  type(ocean_OBC_type),                   pointer    :: OBC
  ! Local variables
  integer :: east_boundary, west_boundary, north_boundary, south_boundary
  integer :: i, j

  if (.not.associated(OBC%OBC_mask_u)) then
    allocate(OBC%OBC_mask_u(G%IsdB:G%IedB,G%jsd:G%jed)) ; OBC%OBC_mask_u(:,:) = .false.
  endif
  if (.not.associated(OBC%OBC_direction_u)) then
    allocate(OBC%OBC_direction_u(G%IsdB:G%IedB,G%jsd:G%jed)) ; OBC%OBC_direction_u(:,:) = OBC_NONE
  endif
  if (.not.associated(OBC%OBC_kind_u)) then
    allocate(OBC%OBC_kind_u(G%IsdB:G%IedB,G%jsd:G%jed)) ; OBC%OBC_kind_u(:,:) = OBC_NONE
  endif
  if (.not.associated(OBC%OBC_mask_v)) then
    allocate(OBC%OBC_mask_v(G%isd:G%ied,G%JsdB:G%JedB)) ; OBC%OBC_mask_v(:,:) = .false.
  endif
  if (.not.associated(OBC%OBC_direction_v)) then
    allocate(OBC%OBC_direction_v(G%isd:G%ied,G%JsdB:G%JedB)) ; OBC%OBC_direction_v(:,:) = OBC_NONE
  endif
  if (.not.associated(OBC%OBC_kind_v)) then
    allocate(OBC%OBC_kind_v(G%isd:G%ied,G%JsdB:G%JedB)) ; OBC%OBC_kind_v(:,:) = OBC_NONE
  endif

  ! This code should be modified to allow OBCs to be applied anywhere.

  if (G%symmetric) then
    east_boundary = G%ieg
    west_boundary = G%isg-1
    north_boundary = G%jeg
    south_boundary = G%jsg-1
  else
    ! I am not entirely sure that this works properly. -RWH
    east_boundary = G%ieg-1
    west_boundary = G%isg
    north_boundary = G%jeg-1
    south_boundary = G%jsg
  endif

  if (OBC%apply_OBC_u_flather_east) then
    ! Determine where u points are applied at east side
    do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
      if ((I+G%idg_offset) == east_boundary) then !eastern side
        OBC%OBC_mask_u(I,j) = .true.
        OBC%OBC_direction_u(I,j) = OBC_DIRECTION_E
        OBC%OBC_kind_u(I,j) = OBC_FLATHER
        OBC%OBC_mask_v(i+1,J) = .true.
        if (OBC%OBC_direction_v(i+1,J) == OBC_NONE) then
          OBC%OBC_direction_v(i+1,J) = OBC_DIRECTION_E
          OBC%OBC_kind_v(i+1,J) = OBC_FLATHER
        endif
        OBC%OBC_mask_v(i+1,J-1) = .true.
        if (OBC%OBC_direction_v(i+1,J-1) == OBC_NONE) then
          OBC%OBC_direction_v(i+1,J-1) = OBC_DIRECTION_E
          OBC%OBC_kind_v(i+1,J-1) = OBC_FLATHER
        endif
      endif
    enddo ; enddo
  endif

  if (OBC%apply_OBC_u_flather_west) then
    ! Determine where u points are applied at west side
    do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
      if ((I+G%idg_offset) == west_boundary) then !western side
        OBC%OBC_mask_u(I,j) = .true.
        OBC%OBC_direction_u(I,j) = OBC_DIRECTION_W
        OBC%OBC_kind_u(I,j) = OBC_FLATHER
        OBC%OBC_mask_v(i,J) = .true.
        if (OBC%OBC_direction_v(i,J) == OBC_NONE) then
          OBC%OBC_direction_v(i,J) = OBC_DIRECTION_W
          OBC%OBC_kind_v(i,J) = OBC_FLATHER
        endif
        OBC%OBC_mask_v(i,J-1) = .true.
        if (OBC%OBC_direction_v(i,J-1) == OBC_NONE) then
          OBC%OBC_direction_v(i,J-1) = OBC_DIRECTION_W
          OBC%OBC_kind_v(i,J-1) = OBC_FLATHER
        endif
      endif
    enddo ; enddo
  endif

  if (OBC%apply_OBC_v_flather_north) then
    ! Determine where v points are applied at north side
    do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
      if ((J+G%jdg_offset) == north_boundary) then         !northern side
        OBC%OBC_mask_v(i,J) = .true.
        OBC%OBC_direction_v(i,J) = OBC_DIRECTION_N
        OBC%OBC_kind_v(i,J) = OBC_FLATHER
        OBC%OBC_mask_u(I,j+1) = .true.
        if (OBC%OBC_direction_u(I,j+1) == OBC_NONE) then
          OBC%OBC_direction_u(I,j+1) = OBC_DIRECTION_N
          OBC%OBC_kind_u(I,j+1) = OBC_FLATHER
        endif
        OBC%OBC_mask_u(I-1,j+1) = .true.
        if (OBC%OBC_direction_u(I-1,j+1) == OBC_NONE) then
          OBC%OBC_direction_u(I-1,j+1) = OBC_DIRECTION_N
          OBC%OBC_kind_u(I-1,j+1) = OBC_FLATHER
        endif
      endif
    enddo ; enddo
  endif

  if (OBC%apply_OBC_v_flather_south) then
    ! Determine where v points are applied at south side
    do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
      if ((J+G%jdg_offset) == south_boundary) then         !southern side
        OBC%OBC_mask_v(i,J) = .true.
        OBC%OBC_direction_v(i,J) = OBC_DIRECTION_S
        OBC%OBC_kind_v(i,J) = OBC_FLATHER
        OBC%OBC_mask_u(I,j) = .true.
        if (OBC%OBC_direction_u(I,j) == OBC_NONE) then
          OBC%OBC_direction_u(I,j) = OBC_DIRECTION_S
          OBC%OBC_kind_u(I,j) = OBC_FLATHER
        endif
        OBC%OBC_mask_u(I-1,j) = .true.
        if (OBC%OBC_direction_u(I-1,j) == OBC_NONE) then
          OBC%OBC_direction_u(I-1,j) = OBC_DIRECTION_S
          OBC%OBC_kind_u(I-1,j) = OBC_FLATHER
        endif
      endif
    enddo ; enddo
  endif

  !   If there are no OBC points on this PE, there is no reason to keep the OBC
  ! type, and it could be deallocated.
end subroutine set_Flather_positions

!> Sets the initial definitions of the characteristic open boundary conditions.
!! \author Mehmet Ilicak
subroutine set_Flather_data(OBC, tv, h, G, PF, tracer_Reg)
  type(ocean_grid_type),                     intent(inout) :: G !< Ocean grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary structure
  type(thermo_var_ptrs),                     intent(inout) :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(inout) :: h !< Thickness
  type(param_file_type),                     intent(in)    :: PF !< Parameter file handle
  type(tracer_registry_type),                pointer       :: tracer_Reg !< Tracer registry
  ! Local variables
  logical :: read_OBC_eta = .false.
  logical :: read_OBC_uv = .false.
  logical :: read_OBC_TS = .false.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: isd_off, jsd_off
  integer :: IsdB, IedB, JsdB, JedB
  character(len=40)  :: mod = "set_Flather_Bdry_Conds" ! This subroutine's name.
  character(len=200) :: filename, OBC_file, inputdir ! Strings for file/path

  real :: temp_u(G%domain%niglobal+1,G%domain%njglobal)
  real :: temp_v(G%domain%niglobal,G%domain%njglobal+1)

  real, pointer, dimension(:,:,:) :: &
    OBC_T_u => NULL(), &    ! These arrays should be allocated and set to
    OBC_T_v => NULL(), &    ! specify the values of T and S that should come
    OBC_S_u => NULL(), & 
    OBC_S_v => NULL()     

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call get_param(PF, mod, "READ_OBC_UV", read_OBC_uv, &
                 "If true, read the values for the velocity open boundary \n"//&
                 "conditions from the file specified by OBC_FILE.", &
                 default=.false.)
  call get_param(PF, mod, "READ_OBC_ETA", read_OBC_eta, &
                 "If true, read the values for the sea surface height \n"//&
                 "open boundary conditions from the file specified by \n"//&
                 "OBC_FILE.", default=.false.)
  call get_param(PF, mod, "READ_OBC_TS", read_OBC_TS, &
                 "If true, read the values for the temperature and \n"//&
                 "salinity open boundary conditions from the file \n"//&
                 "specified by OBC_FILE.", default=.false.)
  if (read_OBC_uv .or. read_OBC_eta .or. read_OBC_TS) then
    call get_param(PF, mod, "OBC_FILE", OBC_file, &
                 "The file from which the appropriate open boundary \n"//&
                 "condition values are read.", default="MOM_OBC_FILE.nc")
    call get_param(PF, mod, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    filename = trim(inputdir)//trim(OBC_file)
    call log_param(PF, mod, "INPUTDIR/OBC_FILE", filename)
  endif

  if (.not.associated(OBC%vbt_outer)) then
    allocate(OBC%vbt_outer(isd:ied,JsdB:JedB)) ; OBC%vbt_outer(:,:) = 0.0
  endif

  if (.not.associated(OBC%ubt_outer)) then
    allocate(OBC%ubt_outer(IsdB:IedB,jsd:jed)) ; OBC%ubt_outer(:,:) = 0.0
  endif

  if (.not.associated(OBC%eta_outer_u)) then
    allocate(OBC%eta_outer_u(IsdB:IedB,jsd:jed)) ; OBC%eta_outer_u(:,:) = 0.0
  endif

  if (.not.associated(OBC%eta_outer_v)) then
    allocate(OBC%eta_outer_v(isd:ied,JsdB:JedB)) ; OBC%eta_outer_v(:,:) = 0.0
  endif
  
  if (read_OBC_uv) then
    call read_data(filename, 'ubt', OBC%ubt_outer, &
                   domain=G%Domain%mpp_domain, position=EAST_FACE)
    call read_data(filename, 'vbt', OBC%vbt_outer, &
                   domain=G%Domain%mpp_domain, position=NORTH_FACE)
  endif

  if (read_OBC_eta) then
    call read_data(filename, 'eta_outer_u', OBC%eta_outer_u, &
                   domain=G%Domain%mpp_domain, position=EAST_FACE)
    call read_data(filename, 'eta_outer_v', OBC%eta_outer_v, &
                   domain=G%Domain%mpp_domain, position=NORTH_FACE)
  endif

  call pass_vector(OBC%eta_outer_u,OBC%eta_outer_v,G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
  call pass_vector(OBC%ubt_outer,OBC%vbt_outer,G%Domain)

  ! Define radiation coefficients r[xy]_old_[uvh] as needed.  For now, there are
  ! no radiation conditions applied to the thicknesses, since the thicknesses
  ! might not be physically motivated.  Instead, sponges should be used to
  ! enforce the near-boundary layer structure.
  if (OBC%apply_OBC_u_flather_west .or. OBC%apply_OBC_u_flather_east) then
    allocate(OBC%rx_old_u(IsdB:IedB,jsd:jed,nz)) ; OBC%rx_old_u(:,:,:) = 0.0
 !   allocate(OBC%rx_old_h(Isd:Ied,jsd:jed,nz))   ; OBC%rx_old_h(:,:,:) = 0.0
  endif
  if (OBC%apply_OBC_v_flather_south .or. OBC%apply_OBC_v_flather_north) then
    allocate(OBC%ry_old_v(isd:ied,JsdB:JedB,nz)) ; OBC%ry_old_v(:,:,:) = 0.0
 !   allocate(OBC%ry_old_h(isd:ied,Jsd:Jed,nz))   ; OBC%ry_old_h(:,:,:) = 0.0
  endif


  if (associated(tv%T)) then
    allocate(OBC_T_u(IsdB:IedB,jsd:jed,nz)) ; OBC_T_u(:,:,:) = 0.0
    allocate(OBC_S_u(IsdB:IedB,jsd:jed,nz)) ; OBC_S_u(:,:,:) = 0.0
    allocate(OBC_T_v(isd:ied,JsdB:JedB,nz)) ; OBC_T_v(:,:,:) = 0.0
    allocate(OBC_S_v(isd:ied,JsdB:JedB,nz)) ; OBC_S_v(:,:,:) = 0.0

    if (read_OBC_TS) then
      call read_data(filename, 'OBC_T_u', OBC_T_u, &
                     domain=G%Domain%mpp_domain, position=EAST_FACE)
      call read_data(filename, 'OBC_S_u', OBC_S_u, &
                     domain=G%Domain%mpp_domain, position=EAST_FACE)

      call read_data(filename, 'OBC_T_v', OBC_T_v, &
                     domain=G%Domain%mpp_domain, position=NORTH_FACE)
      call read_data(filename, 'OBC_S_v', OBC_S_v, &
                     domain=G%Domain%mpp_domain, position=NORTH_FACE)
    else
      call pass_var(tv%T, G%Domain)
      call pass_var(tv%S, G%Domain)
      do k=1,nz ; do j=js,je ; do I=is-1,ie
        if (OBC%OBC_mask_u(I,j)) then
          if (OBC%OBC_direction_u(I,j) == OBC_DIRECTION_E) then
            OBC_T_u(I,j,k) = tv%T(i,j,k)
            OBC_S_u(I,j,k) = tv%S(i,j,k)
          elseif (OBC%OBC_direction_u(I,j) == OBC_DIRECTION_W) then
            OBC_T_u(I,j,k) = tv%T(i+1,j,k)
            OBC_S_u(I,j,k) = tv%S(i+1,j,k)
          elseif (G%mask2dT(i,j) + G%mask2dT(i+1,j) > 0) then
            OBC_T_u(I,j,k) = (G%mask2dT(i,j)*tv%T(i,j,k) + G%mask2dT(i+1,j)*tv%T(i+1,j,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i+1,j))
            OBC_S_u(I,j,k) = (G%mask2dT(i,j)*tv%S(i,j,k) + G%mask2dT(i+1,j)*tv%S(i+1,j,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i+1,j))
          else ! This probably shouldn't happen or maybe it doesn't matter?
            OBC_T_u(I,j,k) = 0.5*(tv%T(i,j,k)+tv%T(i+1,j,k))
            OBC_S_u(I,j,k) = 0.5*(tv%S(i,j,k)+tv%S(i+1,j,k))
          endif
        else
          OBC_T_u(I,j,k) = 0.5*(tv%T(i,j,k)+tv%T(i+1,j,k))
          OBC_S_u(I,j,k) = 0.5*(tv%S(i,j,k)+tv%S(i+1,j,k))
        endif
      enddo; enddo ; enddo

      do k=1,nz ; do J=js-1,je ; do i=is,ie
        if (OBC%OBC_mask_v(i,J)) then
          if (OBC%OBC_direction_v(i,J) == OBC_DIRECTION_N) then
            OBC_T_v(i,J,k) = tv%T(i,j,k)
            OBC_S_v(i,J,k) = tv%S(i,j,k)
          elseif (OBC%OBC_direction_v(i,J) == OBC_DIRECTION_S) then
            OBC_T_v(i,J,k) = tv%T(i,j+1,k)
            OBC_S_v(i,J,k) = tv%S(i,j+1,k)
          elseif (G%mask2dT(i,j) + G%mask2dT(i,j+1) > 0) then
            OBC_T_v(i,J,k) = (G%mask2dT(i,j)*tv%T(i,j,k) + G%mask2dT(i,j+1)*tv%T(i,j+1,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i,j+1))
            OBC_S_v(i,J,k) = (G%mask2dT(i,j)*tv%S(i,j,k) + G%mask2dT(i,j+1)*tv%S(i,j+1,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i,j+1))
          else ! This probably shouldn't happen or maybe it doesn't matter?
            OBC_T_v(i,J,k) = 0.5*(tv%T(i,j,k)+tv%T(i,j+1,k))
            OBC_S_v(i,J,k) = 0.5*(tv%S(i,j,k)+tv%S(i,j+1,k))
          endif
        else
          OBC_T_v(i,J,k) = 0.5*(tv%T(i,j,k)+tv%T(i,j+1,k))
          OBC_S_v(i,J,k) = 0.5*(tv%S(i,j,k)+tv%S(i,j+1,k))
        endif
      enddo; enddo ; enddo
    endif

    call pass_vector(OBC_T_u, OBC_T_v, G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
    call pass_vector(OBC_S_u, OBC_S_v, G%Domain, To_All+SCALAR_PAIR, CGRID_NE)

    call add_tracer_OBC_values("T", tracer_Reg, OBC_in_u=OBC_T_u, &
                                                OBC_in_v=OBC_T_v)
    call add_tracer_OBC_values("S", tracer_Reg, OBC_in_u=OBC_S_u, &
                                                OBC_in_v=OBC_S_v)
    do k=1,nz ; do j=jsd,jed ; do I=isd,ied-1
      if (OBC%OBC_direction_u(I,j) == OBC_DIRECTION_E) then
        tv%T(i+1,j,k) = tv%T(i,j,k) ; tv%S(i+1,j,k) = tv%S(i,j,k)
      elseif (OBC%OBC_direction_u(I,j) == OBC_DIRECTION_W) then
        tv%T(i,j,k) = tv%T(i+1,j,k) ; tv%S(i,j,k) = tv%S(i+1,j,k)
      endif
    enddo ; enddo ; enddo
    do k=1,nz ; do J=jsd,jed-1 ; do i=isd,ied
      if (OBC%OBC_direction_v(i,J) == OBC_DIRECTION_N) then
        tv%T(i,j+1,k) = tv%T(i,j,k) ; tv%S(i,j+1,k) = tv%S(i,j,k)
      elseif (OBC%OBC_direction_v(i,J) == OBC_DIRECTION_S) then
        tv%T(i,j,k) = tv%T(i,j+1,k) ; tv%S(i,j,k) = tv%S(i,j+1,k)
      endif
    enddo ; enddo ; enddo
  endif

  do k=1,nz ; do j=jsd,jed ; do I=isd,ied-1
    if (OBC%OBC_direction_u(I,j) == OBC_DIRECTION_E) h(i+1,j,k) = h(i,j,k)
    if (OBC%OBC_direction_u(I,j) == OBC_DIRECTION_W) h(i,j,k) = h(i+1,j,k)
  enddo ; enddo ; enddo
  do k=1,nz ; do J=jsd,jed-1 ; do i=isd,ied
    if (OBC%OBC_direction_v(i,J) == OBC_DIRECTION_N) h(i,j+1,k) = h(i,j,k)
    if (OBC%OBC_direction_v(i,J) == OBC_DIRECTION_S) h(i,j,k) = h(i,j+1,k)
  enddo ; enddo ; enddo

end subroutine set_Flather_data

!> \namespace mom_open_boundary
!! This module implements some aspects of internal open boundary
!! conditions in MOM.
!!
!! A small fragment of the grid is shown below:
!!
!!    j+1  x ^ x ^ x   At x:  q, CoriolisBu
!!    j+1  > o > o >   At ^:  v, tauy
!!    j    x ^ x ^ x   At >:  u, taux
!!    j    > o > o >   At o:  h, bathyT, buoy, tr, T, S, Rml, ustar
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!!
!! The boundaries always run through q grid points (x).

end module MOM_open_boundary
