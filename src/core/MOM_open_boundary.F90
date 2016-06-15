! This file is part of MOM6. See LICENSE.md for the license.
!> Controls where open boundary conditions are applied 
module MOM_open_boundary

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_domains, only : pass_var, pass_vector
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type

implicit none ; private

#include <MOM_memory.h>

public Radiation_Open_Bdry_Conds, open_boundary_init, open_boundary_end

!> The control structure for open-boundaries
type, public :: open_boundary_CS ; private
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
end type open_boundary_CS

integer, parameter, public :: OBC_NONE = 0, OBC_SIMPLE = 1, OBC_WALL = 2
integer, parameter, public :: OBC_FLATHER_E = 4, OBC_FLATHER_W = 5
integer, parameter, public :: OBC_FLATHER_N = 6, OBC_FLATHER_S = 7

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
  ! points, and can be OBC_NONE, OBC_SIMPLE, OBC_WALL, or one of OBC_FLATHER_[EWNS].  Generally these
  ! should be consistent with OBC_mask_[uv], with OBC_mask_[uv] .false. for OBC_kind_[uv] = NONE
  ! and true for all other values.
  integer, pointer, dimension(:,:) :: &
    OBC_kind_u => NULL(), & !< Type of OBC at u-points.
    OBC_kind_v => NULL()    !< Type of OBC at v-points.
  ! The following apply at points with OBC_kind_[uv] = OBC_FLATHER_x.
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
end type ocean_OBC_type

integer :: id_clock_pass

character(len=40)  :: mod = "MOM_open_boundary" ! This module's name.
! This include declares and sets the variable "version".
#include "version_variable.h"

contains

!> Diagnose radiation conditions at open boundaries
subroutine Radiation_Open_Bdry_Conds(OBC, u_new, u_old, v_new, v_old, &
                                     h_new, h_old, G, CS)
  type(ocean_grid_type),                     intent(inout) :: G !< Ocean grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary data
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u_new
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u_old
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v_new
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v_old
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h_new
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h_old
  type(open_boundary_CS),                    pointer       :: CS !< Open boundary control structure
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
  if (.not.associated(CS)) call MOM_error(FATAL, &
         "MOM_open_boundary: Module must be initialized before it is used.")

  gamma_u = CS%gamma_uv ; gamma_v = CS%gamma_uv ; gamma_h = CS%gamma_h
  rx_max = CS%rx_max ; ry_max = CS%rx_max

  if (OBC%apply_OBC_u_flather_east .or. OBC%apply_OBC_u_flather_west) then
    do k=1,nz ; do j=js,je ; do I=is-1,ie ; if (OBC%OBC_mask_u(I,j)) then
      if (OBC%OBC_kind_u(I,j) == OBC_FLATHER_E) then
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
      if (OBC%OBC_kind_u(I,j) == OBC_FLATHER_W) then
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
      if (OBC%OBC_kind_v(i,J) == OBC_FLATHER_N) then
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

      if (OBC%OBC_kind_v(i,J) == OBC_FLATHER_S) then
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

!> Initialize open boundary control structure
subroutine open_boundary_init(Time, G, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G !< Ocean grid structure
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handle
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(open_boundary_CS),  pointer       :: CS !< Open boundary control structure
  ! Local variables
  logical :: flather_east, flather_west, flather_north, flather_south

  if (associated(CS)) then
    call MOM_error(WARNING, "MOM_open_boundary: open_boundary_init called with associated control structure.")
    return
  endif

  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "APPLY_OBC_U_FLATHER_EAST", flather_east, &
                 "If true, some zonal velocity points use Flather open \n"//&
                 "boundary conditions on the east side of the ocean.", &
                 default=.false.)
  call get_param(param_file, mod, "APPLY_OBC_U_FLATHER_WEST", flather_west, &
                 "If true, some zonal velocity points use Flather open \n"//&
                 "boundary conditions on the west side of the ocean.", &
                 default=.false.)
  call get_param(param_file, mod, "APPLY_OBC_V_FLATHER_NORTH", flather_north, &
                 "If true, some meridional velocity points use Flather \n"//&
                 "open boundary conditions on the north side of the ocean.", &
                 default=.false.)
  call get_param(param_file, mod, "APPLY_OBC_V_FLATHER_SOUTH", flather_south, &
                 "If true, some meridional velocity points use Flather \n"//&
                 "open boundary conditions on the north side of the ocean.", &
                 default=.false.)
  if (.not.(flather_east .or. flather_west .or. flather_north .or. &
            flather_south)) return

  allocate(CS)
  call get_param(param_file, mod, "OBC_RADIATION_MAX", CS%rx_max, &
                 "The maximum magnitude of the baroclinic radiation \n"//&
                 "velocity (or speed of characteristics).  This is only \n"//&
                 "used if one of the APPLY_OBC_[UV]_FLATHER_... is true.", &
                 units="m s-1", default=10.0)
  call get_param(param_file, mod, "OBC_RAD_VEL_WT", CS%gamma_uv, &
                 "The relative weighting for the baroclinic radiation \n"//&
                 "velocities (or speed of characteristics) at the new \n"//&
                 "time level (1) or the running mean (0) for velocities. \n"//&
                 "Valid values range from 0 to 1. This is only used if \n"//&
                 "one of the APPLY_OBC_[UV]_FLATHER_...  is true.", &
                 units="nondim",  default=0.3)
  call get_param(param_file, mod, "OBC_RAD_THICK_WT", CS%gamma_h, &
                 "The relative weighting for the baroclinic radiation \n"//&
                 "velocities (or speed of characteristics) at the new \n"//&
                 "time level (1) or the running mean (0) for thicknesses. \n"//&
                 "Valid values range from 0 to 1. This is only used if \n"//&
                 "one of the APPLY_OBC_[UV]_FLATHER_...  is true.", &
                 units="nondim",  default=0.2)

  id_clock_pass = cpu_clock_id('(Ocean OBC halo updates)', grain=CLOCK_ROUTINE)

end subroutine open_boundary_init

!> Deallocate open boundary data
subroutine open_boundary_end(CS)
  type(open_boundary_CS), pointer :: CS !< Open boundary control structure
  deallocate(CS)
end subroutine open_boundary_end

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
