!> Solve the layer continuity equation using the PPM method for layer fluxes.
module MOM_continuity_PPM

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_segment_type, OBC_NONE
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : BT_cont_type, porous_barrier_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public continuity_PPM, continuity_PPM_init, continuity_PPM_stencil
public continuity_fluxes, continuity_adjust_vel
public zonal_mass_flux, meridional_mass_flux
public zonal_edge_thickness, meridional_edge_thickness
public continuity_zonal_convergence, continuity_merdional_convergence
public zonal_flux_thickness, meridional_flux_thickness
public zonal_BT_mass_flux, meridional_BT_mass_flux
public set_continuity_loop_bounds

!>@{ CPU time clock IDs
integer :: id_clock_reconstruct, id_clock_update, id_clock_correct
!>@}

!> Control structure for mom_continuity_ppm
type, public :: continuity_PPM_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  type(diag_ctrl), pointer :: diag !< Diagnostics control structure.
  logical :: upwind_1st      !< If true, use a first-order upwind scheme.
  logical :: monotonic       !< If true, use the Colella & Woodward monotonic
                             !! limiter; otherwise use a simple positive
                             !! definite limiter.
  logical :: simple_2nd      !< If true, use a simple second order (arithmetic
                             !! mean) interpolation of the edge values instead
                             !! of the higher order interpolation.
  real :: tol_eta            !< The tolerance for free-surface height
                             !! discrepancies between the barotropic solution and
                             !! the sum of the layer thicknesses [H ~> m or kg m-2].
  real :: tol_vel            !< The tolerance for barotropic velocity
                             !! discrepancies between the barotropic solution and
                             !! the sum of the layer thicknesses [L T-1 ~> m s-1].
  real :: CFL_limit_adjust   !< The maximum CFL of the adjusted velocities [nondim]
  logical :: aggress_adjust  !< If true, allow the adjusted velocities to have a
                             !! relative CFL change up to 0.5.  False by default.
  logical :: vol_CFL         !< If true, use the ratio of the open face lengths
                             !! to the tracer cell areas when estimating CFL
                             !! numbers.  Without aggress_adjust, the default is
                             !! false; it is always true with.
  logical :: better_iter     !< If true, stop corrective iterations using a
                             !! velocity-based criterion and only stop if the
                             !! iteration is better than all predecessors.
  logical :: use_visc_rem_max !< If true, use more appropriate limiting bounds
                             !! for corrections in strongly viscous columns.
  logical :: marginal_faces  !< If true, use the marginal face areas from the
                             !! continuity solver for use as the weights in the
                             !! barotropic solver.  Otherwise use the transport
                             !! averaged areas.
end type continuity_PPM_CS

!> A container for loop bounds
type, public :: cont_loop_bounds_type ; private
  !>@{ Loop bounds
  integer :: ish, ieh, jsh, jeh
  !>@}
end type cont_loop_bounds_type

!> Finds the thickness fluxes from the continuity solver or their vertical sum without
!! actually updating the layer thicknesses.
interface continuity_fluxes
  module procedure continuity_3d_fluxes, continuity_2d_fluxes
end interface continuity_fluxes

contains

!> Time steps the layer thicknesses, using a monotonically limit, directionally split PPM scheme,
!! based on Lin (1994).
subroutine continuity_PPM(u, v, hin, h, uh, vh, dt, G, GV, US, CS, OBC, pbv, uhbt, vhbt, &
                          visc_rem_u, visc_rem_v, u_cor, v_cor, BT_cont, du_cor, dv_cor)
  type(ocean_grid_type),   intent(in)    :: G   !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u   !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v   !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: hin !< Initial layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h   !< Final layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: uh  !< Zonal volume flux, u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out)   :: vh  !< Meridional volume flux, v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS  !< Module's control structure.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< pointers to porous barrier fractional cell metrics
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(in)    :: uhbt !< The summed volume flux through zonal faces
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(in)    :: vhbt !< The summed volume flux through meridional faces
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_u
                             !< The fraction of zonal momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_v
                             !< The fraction of meridional momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: u_cor
                             !< The zonal velocities that give uhbt as the depth-integrated transport [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(out)   :: v_cor
                             !< The meridional velocities that give vhbt as the depth-integrated
                             !! transport [L T-1 ~> m s-1].
  type(BT_cont_type), optional, pointer  :: BT_cont !< A structure with elements that describe
                             !!  the effective open face areas as a function of barotropic flow.
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(out)   :: du_cor !< The zonal velocity increments from u that give uhbt
                                                 !! as the depth-integrated transports [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(out)   :: dv_cor !< The meridional velocity increments from v that give vhbt
                                                 !! as the depth-integrated transports [L T-1 ~> m s-1].

  ! Local variables
  real :: h_W(SZI_(G),SZJ_(G),SZK_(GV)) ! West edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_E(SZI_(G),SZJ_(G),SZK_(GV)) ! East edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_S(SZI_(G),SZJ_(G),SZK_(GV)) ! South edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_N(SZI_(G),SZJ_(G),SZK_(GV)) ! North edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_min  ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  type(cont_loop_bounds_type) :: LB ! A type indicating the loop range for a phase of the updates
  logical :: x_first

  h_min = GV%Angstrom_H

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "MOM_continuity_PPM: Module must be initialized before it is used.")

  x_first = (MOD(G%first_direction,2) == 0)

  if (present(visc_rem_u) .neqv. present(visc_rem_v)) call MOM_error(FATAL, &
      "MOM_continuity_PPM: Either both visc_rem_u and visc_rem_v or neither"// &
      " one must be present in call to continuity_PPM.")

  ! update device visc_rem_u for zonal_mass_flux
  !$omp target update to(visc_rem_u, u)

  if (x_first) then
    !  First advect zonally, with loop bounds that accomodate the subsequent meridional advection.
    LB = set_continuity_loop_bounds(G, CS, i_stencil=.false., j_stencil=.true.)
    call zonal_edge_thickness(hin, h_W, h_E, G, GV, US, CS, OBC, LB)
    call zonal_mass_flux(u, hin, h_W, h_E, uh, dt, G, GV, US, CS, OBC, pbv%por_face_areaU, &
                         LB, uhbt, visc_rem_u, u_cor, BT_cont, du_cor)
    ! update device uh from zonal_mass_flux
    !$omp target update to(uh)
    call continuity_zonal_convergence(h, uh, dt, G, GV, LB, hin)

    ! update host h from continuity_zonal_convergence
    !$omp target update from(h)

    !  Now advect meridionally, using the updated thicknesses to determine the fluxes.
    LB = set_continuity_loop_bounds(G, CS, i_stencil=.false., j_stencil=.false.)
    call meridional_edge_thickness(h, h_S, h_N, G, GV, US, CS, OBC, LB)
    call meridional_mass_flux(v, h, h_S, h_N, vh, dt, G, GV, US, CS, OBC, pbv%por_face_areaV, &
                              LB, vhbt, visc_rem_v, v_cor, BT_cont, dv_cor)
    call continuity_merdional_convergence(h, vh, dt, G, GV, LB, hmin=h_min)

  else  ! .not. x_first
    !  First advect meridionally, with loop bounds that accomodate the subsequent zonal advection.
    LB = set_continuity_loop_bounds(G, CS, i_stencil=.true., j_stencil=.false.)
    call meridional_edge_thickness(hin, h_S, h_N, G, GV, US, CS, OBC, LB)
    call meridional_mass_flux(v, hin, h_S, h_N, vh, dt, G, GV, US, CS, OBC, pbv%por_face_areaV, &
                              LB, vhbt, visc_rem_v, v_cor, BT_cont, dv_cor)
    call continuity_merdional_convergence(h, vh, dt, G, GV, LB, hin)

    !  Now advect zonally, using the updated thicknesses to determine the fluxes.
    LB = set_continuity_loop_bounds(G, CS, i_stencil=.false., j_stencil=.false.)
    call zonal_edge_thickness(h, h_W, h_E, G, GV, US, CS, OBC, LB)
    call zonal_mass_flux(u, h, h_W, h_E, uh, dt, G, GV, US, CS, OBC, pbv%por_face_areaU, &
                         LB, uhbt, visc_rem_u, u_cor, BT_cont, du_cor)
    call continuity_zonal_convergence(h, uh, dt, G, GV, LB, hmin=h_min)
  endif

end subroutine continuity_PPM

!> Finds the thickness fluxes from the continuity solver without actually updating the
!! layer thicknesses.  Because the fluxes in the two directions are calculated based on the
!! input thicknesses, which are not updated between the direcitons, the fluxes returned here
!! are not the same as those that would be returned by a call to continuity.
subroutine continuity_3d_fluxes(u, v, h, uh, vh, dt, G, GV, US, CS, OBC, pbv)
  type(ocean_grid_type),   intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u   !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v   !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: uh  !< Thickness fluxes through zonal faces,
                                                !! u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out)   :: vh  !< Thickness fluxes through meridional faces,
                                                !! v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics

  ! Local variables
  real :: h_W(SZI_(G),SZJ_(G),SZK_(GV)) ! West edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_E(SZI_(G),SZJ_(G),SZK_(GV)) ! East edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_S(SZI_(G),SZJ_(G),SZK_(GV)) ! South edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_N(SZI_(G),SZJ_(G),SZK_(GV)) ! North edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]

  call zonal_edge_thickness(h, h_W, h_E, G, GV, US, CS, OBC)
  call zonal_mass_flux(u, h, h_W, h_E, uh, dt, G, GV, US, CS, OBC, pbv%por_face_areaU)

  call meridional_edge_thickness(h, h_S, h_N, G, GV, US, CS, OBC)
  call meridional_mass_flux(v, h, h_S, h_N, vh, dt, G, GV, US, CS, OBC, pbv%por_face_areaV)

end subroutine continuity_3d_fluxes

!> Find the vertical sum of the thickness fluxes from the continuity solver without actually
!! updating the layer thicknesses.  Because the fluxes in the two directions are calculated
!! based on the input thicknesses, which are not updated between the directions, the fluxes
!! returned here are not the same as those that would be returned by a call to continuity.
subroutine continuity_2d_fluxes(u, v, h, uhbt, vhbt, dt, G, GV, US, CS, OBC, pbv)
  type(ocean_grid_type),   intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u   !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v   !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(out)   :: uhbt !< Vertically summed thickness flux through
                                                !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(out)   :: vhbt !< Vertically summed thickness flux through
                                                !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics

  ! Local variables
  real :: h_W(SZI_(G),SZJ_(G),SZK_(GV)) ! West edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_E(SZI_(G),SZJ_(G),SZK_(GV)) ! East edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_S(SZI_(G),SZJ_(G),SZK_(GV)) ! South edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_N(SZI_(G),SZJ_(G),SZK_(GV)) ! North edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]

  call zonal_edge_thickness(h, h_W, h_E, G, GV, US, CS, OBC)
  call zonal_BT_mass_flux(u, h, h_W, h_E, uhbt, dt, G, GV, US, CS, OBC, pbv%por_face_areaU)

  call meridional_edge_thickness(h, h_S, h_N, G, GV, US, CS, OBC)
  call meridional_BT_mass_flux(v, h, h_S, h_N, vhbt, dt, G, GV, US, CS, OBC, pbv%por_face_areaV)

end subroutine continuity_2d_fluxes

!> Correct the velocities to give the specified depth-integrated transports by applying a
!! barotropic acceleration (subject to viscous drag) to the velocities.
subroutine continuity_adjust_vel(u, v, h, dt, G, GV, US, CS, OBC, pbv, uhbt, vhbt, visc_rem_u, visc_rem_v)
  type(ocean_grid_type),   intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u   !< Zonal velocity, which will be adjusted to
                                                !! give uhbt as the depth-integrated
                                                !! transport [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v   !< Meridional velocity, which will be adjusted
                                                !! to give vhbt as the depth-integrated
                                                !! transport [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: uhbt !< The vertically summed thickness flux through
                                                !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(in)    :: vhbt !< The vertically summed thickness flux through
                                                !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_u !< Both the fraction of the zonal momentum
                                                !! that remains after a time-step of viscosity, and
                                                !! the fraction of a time-step's worth of a barotropic
                                                !! acceleration that a layer experiences after viscosity
                                                !! is applied [nondim].  This goes between 0 (at the
                                                !! bottom) and 1 (far above the bottom).  When this
                                                !! column is under an ice shelf, this also goes to 0
                                                !! at the top due to the no-slip boundary condition there.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_v !< Both the fraction of the meridional momentum
                                                !! that remains after a time-step of viscosity, and
                                                !! the fraction of a time-step's worth of a barotropic
                                                !! acceleration that a layer experiences after viscosity
                                                !! is applied [nondim].  This goes between 0 (at the
                                                !! bottom) and 1 (far above the bottom).  When this
                                                !! column is under an ice shelf, this also goes to 0
                                                !! at the top due to the no-slip boundary condition there.

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: u_in  !< Input zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: v_in  !< Input meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uh  !< Volume flux through zonal faces =
                                                !! u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vh  !< Volume flux through meridional faces =
                                                !! v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: h_W(SZI_(G),SZJ_(G),SZK_(GV)) ! West edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_E(SZI_(G),SZJ_(G),SZK_(GV)) ! East edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_S(SZI_(G),SZJ_(G),SZK_(GV)) ! South edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_N(SZI_(G),SZJ_(G),SZK_(GV)) ! North edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]

  ! It might not be necessary to separate the input velocity array from the adjusted velocities,
  ! but it seems safer to do so, even if it might be less efficient.
  u_in(:,:,:) = u(:,:,:)
  v_in(:,:,:) = v(:,:,:)

  call zonal_edge_thickness(h, h_W, h_E, G, GV, US, CS, OBC)
  call zonal_mass_flux(u_in, h, h_W, h_E, uh, dt, G, GV, US, CS, OBC, pbv%por_face_areaU, &
                       uhbt=uhbt, visc_rem_u=visc_rem_u, u_cor=u)

  call meridional_edge_thickness(h, h_S, h_N, G, GV, US, CS, OBC)
  call meridional_mass_flux(v_in, h, h_S, h_N, vh, dt, G, GV, US, CS, OBC, pbv%por_face_areaV, &
                            vhbt=vhbt, visc_rem_v=visc_rem_v, v_cor=v)

end subroutine continuity_adjust_vel


!> Updates the thicknesses due to zonal thickness fluxes.
subroutine continuity_zonal_convergence(h, uh, dt, G, GV, LB, hin, hmin)
  type(ocean_grid_type),       intent(in)    :: G    !< Ocean's grid structure
  type(verticalGrid_type),     intent(in)    :: GV   !< Ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                               intent(inout) :: h    !< Final layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                               intent(in)    :: uh   !< Zonal thickness flux, u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,                        intent(in)    :: dt   !< Time increment [T ~> s]
  type(cont_loop_bounds_type), intent(in)    :: LB   !< Loop bounds structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                     optional, intent(in)    :: hin  !< Initial layer thickness [H ~> m or kg m-2].
                                                     !! If hin is absent, h is also the initial thickness.
  real,              optional, intent(in)    :: hmin !< The minimum layer thickness [H ~> m or kg m-2]

  real :: h_min  ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  integer :: i, j, k, ish, ieh, jsh, jeh, nz

  call cpu_clock_begin(id_clock_update)

  h_min = 0.0 ; if (present(hmin)) h_min = hmin

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  if (present(hin)) then
    !$omp target teams distribute parallel do collapse(3) &
    !$omp   map(to: hin(ish:ieh, :, :), G, G%IareaT(ish:ieh, jsh:jeh), uh(ish-1:ieh, :, :)) &
    !$omp   map(from: h(ish:ieh, :, :))
    do k=1,nz ; do j=jsh,jeh ; do i=ish, ieh
      h(i,j,k) = max( hin(i,j,k) - dt * G%IareaT(i,j) * (uh(I,j,k) - uh(I-1,j,k)), h_min )
    enddo ; enddo ; enddo
  else
    ! untested
    !$omp target teams distribute parallel do collapse(3) &
    !$omp   map(to: G, G%IareaT(ish:ieh, jsh:jeh), uh(ish-1:ieh, :, :)) &
    !$omp   map(tofrom: h(ish:ieh, :, :))
    do k=1,nz ; do j=jsh,jeh ; do i=ish, ieh
      h(i,j,k) = max( h(i,j,k) - dt * G%IareaT(i,j) * (uh(I,j,k) - uh(I-1,j,k)), h_min )
    enddo ; enddo ; enddo
  endif

  call cpu_clock_end(id_clock_update)

end subroutine continuity_zonal_convergence

!> Updates the thicknesses due to meridional thickness fluxes.
subroutine continuity_merdional_convergence(h, vh, dt, G, GV, LB, hin, hmin)
  type(ocean_grid_type),       intent(in)    :: G    !< Ocean's grid structure
  type(verticalGrid_type),     intent(in)    :: GV   !< Ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                               intent(inout) :: h    !< Final layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                               intent(in)    :: vh   !< Meridional thickness flux, v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,                        intent(in)    :: dt   !< Time increment [T ~> s]
  type(cont_loop_bounds_type), intent(in)    :: LB   !< Loop bounds structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                     optional, intent(in)    :: hin  !< Initial layer thickness [H ~> m or kg m-2].
                                                     !! If hin is absent, h is also the initial thickness.
  real,              optional, intent(in)    :: hmin !< The minimum layer thickness [H ~> m or kg m-2]

  real :: h_min  ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  integer :: i, j, k

  call cpu_clock_begin(id_clock_update)

  h_min = 0.0 ; if (present(hmin)) h_min = hmin

  if (present(hin)) then
    !$OMP parallel do default(shared)
    do k=1,GV%ke ; do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j,k) = max( hin(i,j,k) - dt * G%IareaT(i,j) * (vh(i,J,k) - vh(i,J-1,k)), h_min )
    enddo ; enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do k=1,GV%ke ; do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j,k) = max( h(i,j,k) - dt * G%IareaT(i,j) * (vh(i,J,k) - vh(i,J-1,k)), h_min )
    enddo ; enddo ; enddo
  endif

  call cpu_clock_end(id_clock_update)

end subroutine continuity_merdional_convergence


!> Set the reconstructed thicknesses at the eastern and western edges of tracer cells.
subroutine zonal_edge_thickness(h_in, h_W, h_E, G, GV, US, CS, OBC, LB_in)
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean's vertical grid structure.
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Tracer cell layer thickness [H ~> m or kg m-2].
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: h_W  !< Western edge layer thickness [H ~> m or kg m-2].
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: h_E  !< Eastern edge layer thickness [H ~> m or kg m-2].
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS   !< This module's control structure.
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundaries control structure.
  type(cont_loop_bounds_type), &
                 optional, intent(in)    :: LB_in !< Loop bounds structure.

  ! Local variables
  type(cont_loop_bounds_type) :: LB
  integer :: i, j, k, ish, ieh, jsh, jeh, nz

  call cpu_clock_begin(id_clock_reconstruct)

  if (present(LB_in)) then
    LB = LB_in
  else
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
  endif
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  if (CS%upwind_1st) then
    ! untested
    !$omp target teams distribute parallel do collapse(3) &
    !$omp   map(to: h_in(ish-1:ieh+1, :, :)) &
    !$omp   map(from: h_W(ish-1:ieh+1, :, :), h_E(ish-1:ieh+1, :, :))
    do k=1,nz ; do j=jsh,jeh ; do i=ish-1,ieh+1
      h_W(i,j,k) = h_in(i,j,k) ; h_E(i,j,k) = h_in(i,j,k)
    enddo ; enddo ; enddo
  else
    call PPM_reconstruction_x(h_in, h_W, h_E, G, GV, LB, &
                              2.0*GV%Angstrom_H, CS%monotonic, CS%simple_2nd, OBC)
  endif

  call cpu_clock_end(id_clock_reconstruct)

end subroutine zonal_edge_thickness


!> Set the reconstructed thicknesses at the eastern and western edges of tracer cells.
subroutine meridional_edge_thickness(h_in, h_S, h_N, G, GV, US, CS, OBC, LB_in)
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean's vertical grid structure.
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Tracer cell layer thickness [H ~> m or kg m-2].
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: h_S  !< Southern edge layer thickness [H ~> m or kg m-2].
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: h_N  !< Northern edge layer thickness [H ~> m or kg m-2].
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS   !< This module's control structure.
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundaries control structure.
  type(cont_loop_bounds_type), &
                 optional, intent(in)    :: LB_in !< Loop bounds structure.

  ! Local variables
  type(cont_loop_bounds_type) :: LB
  integer :: i, j, k, ish, ieh, jsh, jeh, nz

  call cpu_clock_begin(id_clock_reconstruct)

  if (present(LB_in)) then
    LB = LB_in
  else
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
  endif
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  if (CS%upwind_1st) then
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=jsh-1,jeh+1 ; do i=ish,ieh
      h_S(i,j,k) = h_in(i,j,k) ; h_N(i,j,k) = h_in(i,j,k)
    enddo ; enddo ; enddo
  else
    call PPM_reconstruction_y(h_in, h_S, h_N, G, GV, LB, &
                              2.0*GV%Angstrom_H, CS%monotonic, CS%simple_2nd, OBC)
  endif

  call cpu_clock_end(id_clock_reconstruct)

end subroutine meridional_edge_thickness


!> Calculates the mass or volume fluxes through the zonal faces, and other related quantities.
subroutine zonal_mass_flux(u, h_in, h_W, h_E, uh, dt, G, GV, US, CS, OBC, por_face_areaU, &
                           LB_in, uhbt, visc_rem_u, u_cor, BT_cont, du_cor)
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u    !< Zonal velocity [L T-1 ~> m s-1].
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Layer thickness used to calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_W !< Western edge thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_E !< Eastern edge thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: uh   !< Volume flux through zonal faces = u*h*dy
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS   !< This module's control structure.
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundaries control structure.
  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)), &
                           intent(in)    :: por_face_areaU !< fractional open area of U-faces [nondim]
  type(cont_loop_bounds_type), &
                 optional, intent(in)    :: LB_in !< Loop bounds structure.
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(in)    :: uhbt !< The summed volume flux through zonal faces
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_u
                     !< The fraction of zonal momentum originally in a layer that remains after a
                     !! time-step of viscosity, and the fraction of a time-step's worth of a barotropic
                     !! acceleration that a layer experiences after viscosity is applied [nondim].
                     !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: u_cor
                     !< The zonal velocities (u with a barotropic correction)
                     !! that give uhbt as the depth-integrated transport [L T-1 ~> m s-1]
  type(BT_cont_type), optional, pointer  :: BT_cont !< A structure with elements that describe the
                     !! effective open face areas as a function of barotropic flow.
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(out)   :: du_cor !< The zonal velocity increments from u that give uhbt
                                                 !! as the depth-integrated transports [L T-1 ~> m s-1].

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: duhdu ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    du, &         ! Corrective barotropic change in the velocity to give uhbt [L T-1 ~> m s-1].
    du_min_CFL, & ! Lower limit on du correction to avoid CFL violations [L T-1 ~> m s-1]
    du_max_CFL, & ! Upper limit on du correction to avoid CFL violations [L T-1 ~> m s-1]
    duhdu_tot_0, & ! Summed partial derivative of uh with u [H L ~> m2 or kg m-1].
    uh_tot_0, &   ! Summed transport with no barotropic correction [H L2 T-1 ~> m3 s-1 or kg s-1].
    visc_rem_max  ! The column maximum of visc_rem [nondim].
  logical, dimension(SZIB_(G), SZJ_(G)) :: do_I
  real, dimension(SZIB_(G),SZJ_(G), SZK_(GV)) :: &
    visc_rem_u_tmp      ! A 2-D copy of visc_rem_u or an array of 1's [nondim].
  real, dimension(SZIB_(G)) :: FAuI  ! A list of sums of zonal face areas [H L ~> m2 or kg m-1].
  real :: FA_u    ! A sum of zonal face areas [H L ~> m2 or kg m-1].
  real :: I_vrm   ! 1.0 / visc_rem_max [nondim]
  real :: CFL_dt  ! The maximum CFL ratio of the adjusted velocities divided by
                  ! the time step [T-1 ~> s-1].
  real :: I_dt    ! 1.0 / dt [T-1 ~> s-1].
  real :: du_lim  ! The velocity change that give a relative CFL of 1 [L T-1 ~> m s-1].
  real :: dx_E, dx_W ! Effective x-grid spacings to the east and west [L ~> m].
  type(cont_loop_bounds_type) :: LB
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz
  integer :: l_seg ! The OBC segment number
  logical :: use_visc_rem, set_BT_cont
  logical :: local_specified_BC, local_Flather_OBC, local_open_BC, any_simple_OBC  ! OBC-related logicals
  logical :: simple_OBC_pt(SZIB_(G))  ! Indicates points in a row with specified transport OBCs

  call cpu_clock_begin(id_clock_correct)

  use_visc_rem = present(visc_rem_u)

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))

  local_specified_BC = .false. ; local_Flather_OBC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_u_BCs_exist_globally
    local_Flather_OBC = OBC%Flather_u_BCs_exist_globally
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif ; endif

  if (present(du_cor)) du_cor(:,:) = 0.0

  if (present(LB_in)) then
    LB = LB_in
  else
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
  endif
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  CFL_dt = CS%CFL_limit_adjust / dt
  I_dt = 1.0 / dt
  if (CS%aggress_adjust) CFL_dt = I_dt

  ! a better solution is needed!
  if (.not.use_visc_rem) then
    !$omp target teams distribute parallel do collapse(3) &
    !$omp   map(from: visc_rem_u_tmp(ish-1:ieh, :, :))
    do k=1,nz ; do j=jsh,jeh ; do i=ish-1,ieh
      visc_rem_u_tmp(i,j,k) = 1.0
    enddo ; enddo ; enddo
  else
    !$omp target teams distribute parallel do collapse(3) &
    !$omp   map(to: visc_rem_u(ish-1:ieh, :, :)) &
    !$omp   map(from: visc_rem_u_tmp(ish-1:ieh, :, :))
    do k=1,nz ; do j=jsh,jeh ; do i=ish-1,ieh
      visc_rem_u_tmp(i,j,k) = visc_rem_u(i,j,k)
    enddo ; enddo ; enddo
  end if

  !$omp target teams distribute parallel do collapse(2) &
  !$omp   map(from: do_I(ish-1:ieh, jsh:jeh))
  do j=jsh,jeh ; do i=ish-1,ieh
    do_I(i, j) = .true.
  enddo ; enddo
  
  ! Set uh and duhdu.
  call zonal_flux_layer(u, h_in, h_W, h_E, &
                        uh, duhdu, visc_rem_u_tmp, &
                        dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaU, OBC)
  !$omp target update from(uh)
  
  ! untested!
  if (local_specified_BC) then
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh ; if (OBC%segnum_u(I,j) /= 0) then
      l_seg = abs(OBC%segnum_u(I,j))
      if (OBC%segment(l_seg)%specified) uh(I,j,k) = OBC%segment(l_seg)%normal_trans(I,j,k)
    endif ; enddo ; enddo ; enddo
  endif

  if (present(uhbt) .or. set_BT_cont) then
    if (use_visc_rem.and.CS%use_visc_rem_max) then
      visc_rem_max(:, :) = 0.0
      do k=1,nz ; do j=jsh,jeh ; do i=ish-1,ieh
        visc_rem_max(I,j) = max(visc_rem_max(I,j), visc_rem_u(I,j,k))
      enddo ; enddo ; enddo
    else
      visc_rem_max(:, :) = 1.0
    endif
    !   Set limits on du that will keep the CFL number between -1 and 1.
    ! This should be adequate to keep the root bracketed in all cases.
    do j=jsh,jeh ; do I=ish-1,ieh
      I_vrm = 0.0
      if (visc_rem_max(I,j) > 0.0) I_vrm = 1.0 / visc_rem_max(I,j)
      if (CS%vol_CFL) then
        dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
        dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
      else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif
      du_max_CFL(I,j) = 2.0* (CFL_dt * dx_W) * I_vrm
      du_min_CFL(I,j) = -2.0 * (CFL_dt * dx_E) * I_vrm
      uh_tot_0(I,j) = 0.0 ; duhdu_tot_0(I,j) = 0.0
    enddo ; enddo
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
      duhdu_tot_0(I,j) = duhdu_tot_0(I,j) + duhdu(I, j, k)
      uh_tot_0(I,j) = uh_tot_0(I,j) + uh(I,j,k)
    enddo ; enddo ; enddo
    if (use_visc_rem) then
      if (CS%aggress_adjust) then
        ! untested!
        do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
          if (CS%vol_CFL) then
            dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
            dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
          else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

          du_lim = 0.499*((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k)))
          if (du_max_CFL(I,j) * visc_rem_u_tmp(I,j,k) > du_lim) &
            du_max_CFL(I,j) = du_lim / visc_rem_u_tmp(I,j,k)

          du_lim = 0.499*((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k)))
          if (du_min_CFL(I,j) * visc_rem_u_tmp(I,j,k) < du_lim) &
            du_min_CFL(I,j) = du_lim / visc_rem_u_tmp(I,j,k)
        enddo ; enddo ; enddo
      else
        do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
          if (CS%vol_CFL) then
            dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
            dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
          else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

          if (du_max_CFL(I,j) * visc_rem_u_tmp(I,j,k) > dx_W*CFL_dt - u(I,j,k)*G%mask2dCu(I,j)) &
            du_max_CFL(I,j) = (dx_W*CFL_dt - u(I,j,k)) / visc_rem_u_tmp(I,j,k)
          if (du_min_CFL(I,j) * visc_rem_u_tmp(I,j,k) < -dx_E*CFL_dt - u(I,j,k)*G%mask2dCu(I,j)) &
            du_min_CFL(I,j) = -(dx_E*CFL_dt + u(I,j,k)) / visc_rem_u_tmp(I,j,k)
        enddo ; enddo ; enddo
      endif
    else
      ! untested!
      if (CS%aggress_adjust) then
        do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
          if (CS%vol_CFL) then
            dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
            dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
          else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

          du_max_CFL(I,j) = MIN(du_max_CFL(I,j), 0.499 * &
                      ((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k))) )
          du_min_CFL(I,j) = MAX(du_min_CFL(I,j), 0.499 * &
                      ((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k))) )
        enddo ; enddo ; enddo
      else
        do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
          if (CS%vol_CFL) then
            dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
            dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
          else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

          du_max_CFL(I,j) = MIN(du_max_CFL(I,j), dx_W*CFL_dt - u(I,j,k))
          du_min_CFL(I,j) = MAX(du_min_CFL(I,j), -(dx_E*CFL_dt + u(I,j,k)))
        enddo ; enddo ; enddo
      endif
    endif
    do j=jsh,jeh ; do I=ish-1,ieh
      du_max_CFL(I,j) = max(du_max_CFL(I,j),0.0)
      du_min_CFL(I,j) = min(du_min_CFL(I,j),0.0)
    enddo ; enddo

    any_simple_OBC = .false.
    if (present(uhbt) .or. set_BT_cont) then
      ! untested!
      if (local_specified_BC .or. local_Flather_OBC) then ; do j=jsh,jeh ; do I=ish-1,ieh
        l_seg = abs(OBC%segnum_u(I,j))

        ! Avoid reconciling barotropic/baroclinic transports if transport is specified
        simple_OBC_pt(I) = .false.
        if (l_seg /= OBC_NONE) simple_OBC_pt(I) = OBC%segment(l_seg)%specified
        do_I(I, j) = .not.simple_OBC_pt(I)
        any_simple_OBC = any_simple_OBC .or. simple_OBC_pt(I)
      enddo ; enddo ; else ; do j=jsh,jeh ; do I=ish-1,ieh
        do_I(I, j) = .true.
      enddo ; enddo ; endif
      if (present(uhbt)) then
        ! Find du and uh.
        call zonal_flux_adjust(u, h_in, h_W, h_E, uhbt, uh_tot_0, duhdu_tot_0, du, &
                              du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem_u_tmp, &
                              ish, ieh, jsh, jeh, do_I, por_face_areaU, uh, OBC=OBC)
        !$omp target update from(uh)

        if (present(u_cor)) then
          do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
            u_cor(i,j,k) = u(i,j,k) + du(i,j) * visc_rem_u_tmp(i,j,k)
          enddo ; enddo ; enddo
          if (any_simple_OBC) then 
            do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh ; if (simple_OBC_pt(I)) then
              u_cor(I,j,k) = OBC%segment(abs(OBC%segnum_u(I,j)))%normal_vel(I,j,k)
            endif ; enddo ; enddo ; enddo
          endif
        endif ! u-corrected

        if (present(du_cor)) then
          do j=jsh,jeh ; do I=ish-1,ieh ; du_cor(I,j) = du(I,j) ; enddo ; enddo
        endif
      endif
      if (set_BT_cont) then
        call set_zonal_BT_cont(u, h_in, h_W, h_E, BT_cont, uh_tot_0, duhdu_tot_0,&
                               du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem_u_tmp, &
                               visc_rem_max, ish, ieh, jsh, jeh, do_I, por_face_areaU)
      endif
    endif

  endif

  do j=jsh,jeh

    if (present(uhbt) .or. set_BT_cont) then

      if (set_BT_cont) then
        if (any_simple_OBC) then
          do I=ish-1,ieh
            if (simple_OBC_pt(I)) FAuI(I) = GV%H_subroundoff*G%dy_Cu(I,j)
          enddo
          ! NOTE: simple_OBC_pt(I) should prevent access to segment OBC_NONE
          do k=1,nz ; do I=ish-1,ieh ; if (simple_OBC_pt(I)) then
            l_seg = abs(OBC%segnum_u(I,j))
            if ((abs(OBC%segment(l_seg)%normal_vel(I,j,k)) > 0.0) .and. (OBC%segment(l_seg)%specified)) &
              FAuI(I) = FAuI(I) + OBC%segment(l_seg)%normal_trans(I,j,k) / OBC%segment(l_seg)%normal_vel(I,j,k)
          endif ; enddo ; enddo
          do I=ish-1,ieh ; if (simple_OBC_pt(I)) then
            BT_cont%FA_u_W0(I,j) = FAuI(I) ; BT_cont%FA_u_E0(I,j) = FAuI(I)
            BT_cont%FA_u_WW(I,j) = FAuI(I) ; BT_cont%FA_u_EE(I,j) = FAuI(I)
            BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
          endif ; enddo
        endif
      endif ! set_BT_cont

    endif ! present(uhbt) or set_BT_cont

  enddo ! j-loop

  if (local_open_BC .and. set_BT_cont) then
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          do j = OBC%segment(n)%HI%Jsd, OBC%segment(n)%HI%Jed
            FA_u = 0.0
            do k=1,nz ; FA_u = FA_u + h_in(i,j,k)*(G%dy_Cu(I,j)*por_face_areaU(I,j,k)) ; enddo
            BT_cont%FA_u_W0(I,j) = FA_u ; BT_cont%FA_u_E0(I,j) = FA_u
            BT_cont%FA_u_WW(I,j) = FA_u ; BT_cont%FA_u_EE(I,j) = FA_u
            BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
          enddo
        else
          do j = OBC%segment(n)%HI%Jsd, OBC%segment(n)%HI%Jed
            FA_u = 0.0
            do k=1,nz ; FA_u = FA_u + h_in(i+1,j,k)*(G%dy_Cu(I,j)*por_face_areaU(I,j,k)) ; enddo
            BT_cont%FA_u_W0(I,j) = FA_u ; BT_cont%FA_u_E0(I,j) = FA_u
            BT_cont%FA_u_WW(I,j) = FA_u ; BT_cont%FA_u_EE(I,j) = FA_u
            BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
          enddo
        endif
      endif
    enddo
  endif

  if  (set_BT_cont) then ; if (allocated(BT_cont%h_u)) then
    if (present(u_cor)) then
      call zonal_flux_thickness(u_cor, h_in, h_W, h_E, BT_cont%h_u, dt, G, GV, US, LB, &
                                CS%vol_CFL, CS%marginal_faces, OBC, por_face_areaU, visc_rem_u)
    else
      call zonal_flux_thickness(u, h_in, h_W, h_E, BT_cont%h_u, dt, G, GV, US, LB, &
                                CS%vol_CFL, CS%marginal_faces, OBC, por_face_areaU, visc_rem_u)
    endif
  endif ; endif

  call cpu_clock_end(id_clock_correct)

end subroutine zonal_mass_flux


!> Calculates the vertically integrated mass or volume fluxes through the zonal faces.
subroutine zonal_BT_mass_flux(u, h_in, h_W, h_E, uhbt, dt, G, GV, US, CS, OBC, por_face_areaU, LB_in)
  type(ocean_grid_type),                      intent(in)  :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: u    !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_in !< Layer thickness used to
                                                                  !! calculate fluxes [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_W  !< Western edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_E  !< Eastern edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)),          intent(out) :: uhbt !< The summed volume flux through zonal
                                                                  !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                                       intent(in)  :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                      intent(in)  :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),                    intent(in)  :: CS   !< This module's control structure.G
  type(ocean_OBC_type),                       pointer     :: OBC  !< Open boundary condition type
                                                                  !! specifies whether, where, and what
                                                                  !! open boundary conditions are used.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)),  intent(in)  :: por_face_areaU !< fractional open area of U-faces [nondim]
  type(cont_loop_bounds_type),      optional, intent(in)  :: LB_in !< Loop bounds structure.

  ! Local variables
  real :: uh(SZIB_(G),SZJ_(G),SZK_(GV))      ! Volume flux through zonal faces = u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: duhdu(SZIB_(G),SZJ_(G),SZK_(GV))   ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
  logical, dimension(SZIB_(G),SZJ_(G)) :: do_I
  real :: ones(SZIB_(G),SZJ_(G),SZK_(GV))    ! An array of 1's [nondim]
  integer :: i, j, k, ish, ieh, jsh, jeh, nz, l_seg
  logical :: local_specified_BC
  logical, dimension(SZJ_(G)) :: OBC_in_row

  call cpu_clock_begin(id_clock_correct)

  local_specified_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_v_BCs_exist_globally
  endif ; endif

  if (present(LB_in)) then
    ish = LB_in%ish ; ieh = LB_in%ieh ; jsh = LB_in%jsh ; jeh = LB_in%jeh ; nz = GV%ke
  else
    ish = G%isc ; ieh = G%iec ; jsh = G%jsc ; jeh = G%jec ; nz = GV%ke
  endif

  ones(:,:,:) = 1.0 ; do_I(:,:) = .true. ; OBC_in_row(:) = .false.

  uhbt(:,:) = 0.0

  ! Determining whether there are any OBC points outside of the k-loop should be more efficient.
  if (local_specified_BC) then
    do j=jsh,jeh ; do I=ish-1,ieh ; if (OBC%segnum_u(I,j) /= 0) then
      if (OBC%segment(abs(OBC%segnum_u(I,j)))%specified) OBC_in_row(j) = .true.
    endif ; enddo ; enddo
  endif

  ! This sets uh and duhdu.
  call zonal_flux_layer(u, h_in, h_W, h_E, uh, duhdu, ones, &
                        dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaU, OBC)
  
  do k=1,nz ; do j=jsh,jeh ; do i=ish-1,ieh
    if (OBC_in_row(j) .and. OBC%segnum_u(I,j) /= 0) then
      l_seg = abs(OBC%segnum_u(I,j))
      if (OBC%segment(l_seg)%specified) uh(I,j,k) = OBC%segment(l_seg)%normal_trans(I,j,k)
    endif
  enddo ; enddo ; enddo

  ! Accumulate the barotropic transport.
  do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
        uhbt(I,j) = uhbt(I,j) + uh(I,j,k)
  enddo ; enddo ; enddo ! j-loop

  call cpu_clock_end(id_clock_correct)

end subroutine zonal_BT_mass_flux

!> Evaluates the zonal mass or volume fluxes in a layer.
subroutine zonal_flux_layer(u, h, h_W, h_E, uh, duhdu, visc_rem, dt, G, GV, US, &
                            ish, ieh, jsh, jeh, nz, do_I, vol_CFL, por_face_areaU, OBC)
  type(ocean_grid_type),    intent(in)    :: G        !< Ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV       !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: u        !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: visc_rem !< Both the fraction of the
                        !! momentum originally in a layer that remains after a time-step
                        !! of viscosity, and the fraction of a time-step's worth of a barotropic
                        !! acceleration that a layer experiences after viscosity is applied [nondim].
                        !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h        !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h_W      !< West edge thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                            intent(in)    :: h_E      !< East edge thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(inout) :: uh       !< Zonal mass or volume
                                                      !! transport [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(inout) :: duhdu    !< Partial derivative of uh
                                                      !! with u [H L ~> m2 or kg m-1].
  real,                     intent(in)    :: dt       !< Time increment [T ~> s]
  type(unit_scale_type),    intent(in)    :: US       !< A dimensional unit scaling type.
  integer,                  intent(in)    :: ish      !< Start of i index range.
  integer,                  intent(in)    :: ieh      !< End of i index range.
  integer,                  intent(in)    :: jsh      !< Start of j index range.
  integer,                  intent(in)    :: jeh      !< End of j index range.
  integer,                  intent(in)    :: nz       !< Edn of k index range.
  logical, dimension(SZIB_(G),SZJ_(G)), &
                            intent(in)    :: do_I     !< Which i values to work on.
  logical,                  intent(in)    :: vol_CFL  !< If true, rescale the
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: por_face_areaU !< fractional open area of U-faces 
                                                            !! [nondim].
          !! ratio of face areas to the cell areas when estimating the CFL number.
  type(ocean_OBC_type), optional, pointer :: OBC !< Open boundaries control structure.
  ! Local variables
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length [H ~> m or kg m-2]
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  integer :: i, j, k
  integer :: l_seg
  logical :: local_open_BC

  local_open_BC = .false.
  if (present(OBC)) then ; if (associated(OBC)) then
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif ; endif

  !$omp target teams distribute parallel do collapse(3) &
  !$omp   private(CFL, curv_3, h_marg) &
  !$omp   map(to: do_I(ish-1:ieh, jsh:jeh), u(ish-1:ieh, :, :), G, G%dy_Cu(ish-1:ieh, jsh:jeh), &
  !$omp       G%IareaT(ish-1:ieh+1, jsh:jeh), G%IdxT(ish-1:ieh+1, jsh:jeh), h_W(ish-1:ieh+1, :, :), &
  !$omp       h_E(ish-1:ieh+1, :, :), h(ish-1:ieh+1, :, :), por_face_areaU(ish-1:ieh, :, :), &
  !$omp       visc_rem(ish-1:ieh, :, :)) &
  !$omp   map(tofrom: uh(ish-1:ieh, :, :), duhdu(ish-1:ieh, :, :)) ! tofrom so non-updated elems so we don't accidentally zero old values
  do k = 1, nz; do j = jsh, jeh; do I=ish-1,ieh ; if (do_I(I, j)) then
    ! Set new values of uh and duhdu.
    if (u(I, j, k) > 0.0) then
      if (vol_CFL) then ; CFL = (u(I, j, k) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
      else ; CFL = u(I, j, k) * dt * G%IdxT(i,j) ; endif
      curv_3 = (h_W(i, j, k) + h_E(i, j, k)) - 2.0*h(i, j, k)
      uh(i, j, k) = (G%dy_Cu(I,j) * por_face_areaU(I, j, k)) * u(I, j, k) * &
          (h_E(i, j, k) + CFL * (0.5*(h_W(i, j, k) - h_E(i, j, k)) + curv_3*(CFL - 1.5)))
      h_marg = h_E(i, j, k) + CFL * ((h_W(i, j, k) - h_E(i, j, k)) + 3.0*curv_3*(CFL - 1.0))
    elseif (u(I, j, k) < 0.0) then
      if (vol_CFL) then ; CFL = (-u(I, j, k) * dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
      else ; CFL = -u(I, j, k) * dt * G%IdxT(i+1,j) ; endif
      curv_3 = (h_W(i+1, j, k) + h_E(i+1, j, k)) - 2.0*h(i+1, j, k)
      uh(i, j, k) = (G%dy_Cu(I,j) * por_face_areaU(I, j, k)) * u(I, j, k) * &
          (h_W(i+1, j, k) + CFL * (0.5*(h_E(i+1, j, k)-h_W(i+1, j, k)) + curv_3*(CFL - 1.5)))
      h_marg = h_W(i+1, j, k) + CFL * ((h_E(i+1, j, k)-h_W(i+1, j, k)) + 3.0*curv_3*(CFL - 1.0))
    else
      uh(i, j, k) = 0.0
      h_marg = 0.5 * (h_W(i+1, j, k) + h_E(i, j, k))
    endif
    duhdu(I, j, k) = (G%dy_Cu(I,j) * por_face_areaU(I, j, k)) * h_marg * visc_rem(I, j, k)
  endif ; enddo; enddo; enddo

  if (local_open_BC) then
    ! untested
    !$omp target teams distribute parallel do collapse(3) &
    !$omp   map(to: do_I(ish-1:ieh, jsh:jeh), OBC, OBC%segnum_u(ish-1:ieh, jsh:jeh), &
    !$omp       OBC%segment(:), G, G%dy_Cu(ish-1:ieh, jsh:jeh), por_face_areaU(ish-1:ieh, :, :), &
    !$omp       u(ish-1:ieh, :, :), h(ish-1:ieh+1, :, :), visc_rem(ish-1:ieh, :, :)) &
    !$Omp   map(tofrom: uh(ish-1:ieh, :, :), duhdu(ish-1:ieh, :, :))
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh 
      if (do_I(I, j)) then ; if (OBC%segnum_u(I,j) /= 0) then
        if (OBC%segment(abs(OBC%segnum_u(I,j)))%open) then
          if (OBC%segnum_u(I,j) > 0) then !  OBC_DIRECTION_E
            uh(i, j, k) = (G%dy_Cu(I,j) * por_face_areaU(I, j, k)) * u(I, j, k) * h(i, j, k)
            duhdu(I, j, k) = (G%dy_Cu(I,j) * por_face_areaU(I, j, k)) * h(i, j, k) * visc_rem(I, j, k)
          else !  OBC_DIRECTION_W
            uh(i, j, k) = (G%dy_Cu(I,j) * por_face_areaU(I, j, k)) * u(I, j, k) * h(i+1, j, k)
            duhdu(I, j, k) = (G%dy_Cu(I,j)* por_face_areaU(I, j, k)) * h(i+1, j, k) * visc_rem(I, j, k)
          endif
        endif
      endif ; endif
    enddo ; enddo ; enddo
  endif

end subroutine zonal_flux_layer


!> Sets the effective interface thickness associated with the fluxes at each zonal velocity point,
!! optionally scaling back these thicknesses to account for viscosity and fractional open areas.
subroutine zonal_flux_thickness(u, h, h_W, h_E, h_u, dt, G, GV, US, LB, vol_CFL, &
                                marginal, OBC, por_face_areaU, visc_rem_u)
  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)   :: u    !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h    !< Layer thickness used to
                                                                   !! calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_W  !< West edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_E  !< East edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h_u !< Effective thickness at zonal faces,
                                                                   !! scaled down to account for the effects of
                                                                   !! viscosity and the fractional open area
                                                                   !! [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
  type(cont_loop_bounds_type),               intent(in)    :: LB   !< Loop bounds structure.
  logical,                                   intent(in)    :: vol_CFL !< If true, rescale the ratio
                          !! of face areas to the cell areas when estimating the CFL number.
  logical,                                   intent(in)    :: marginal !< If true, report the
                          !! marginal face thicknesses; otherwise report transport-averaged thicknesses.
  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)), &
                                   intent(in)    :: por_face_areaU !< fractional open area of U-faces [nondim]
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                   optional, intent(in)    :: visc_rem_u
                          !< Both the fraction of the momentum originally in a layer that remains after
                          !! a time-step of viscosity, and the fraction of a time-step's worth of a
                          !! barotropic acceleration that a layer experiences after viscosity is applied [nondim].
                          !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).

  ! Local variables
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length [H ~> m or kg m-2]
  real :: h_avg  ! The average thickness of a flux [H ~> m or kg m-2].
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  logical :: local_open_BC
  integer :: i, j, k, ish, ieh, jsh, jeh, nz, n
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  !$OMP parallel do default(shared) private(CFL,curv_3,h_marg,h_avg)
  do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
    if (u(I,j,k) > 0.0) then
      if (vol_CFL) then ; CFL = (u(I,j,k) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
      else ; CFL = u(I,j,k) * dt * G%IdxT(i,j) ; endif
      curv_3 = (h_W(i,j,k) + h_E(i,j,k)) - 2.0*h(i,j,k)
      h_avg = h_E(i,j,k) + CFL * (0.5*(h_W(i,j,k) - h_E(i,j,k)) + curv_3*(CFL - 1.5))
      h_marg = h_E(i,j,k) + CFL * ((h_W(i,j,k) - h_E(i,j,k)) + 3.0*curv_3*(CFL - 1.0))
    elseif (u(I,j,k) < 0.0) then
      if (vol_CFL) then ; CFL = (-u(I,j,k)*dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
      else ; CFL = -u(I,j,k) * dt * G%IdxT(i+1,j) ; endif
      curv_3 = (h_W(i+1,j,k) + h_E(i+1,j,k)) - 2.0*h(i+1,j,k)
      h_avg = h_W(i+1,j,k) + CFL * (0.5*(h_E(i+1,j,k)-h_W(i+1,j,k)) + curv_3*(CFL - 1.5))
      h_marg = h_W(i+1,j,k) + CFL * ((h_E(i+1,j,k)-h_W(i+1,j,k)) + &
                                    3.0*curv_3*(CFL - 1.0))
    else
      h_avg = 0.5 * (h_W(i+1,j,k) + h_E(i,j,k))
      !   The choice to use the arithmetic mean here is somewhat arbitrarily, but
      ! it should be noted that h_W(i+1,j,k) and h_E(i,j,k) are usually the same.
      h_marg = 0.5 * (h_W(i+1,j,k) + h_E(i,j,k))
 !    h_marg = (2.0 * h_W(i+1,j,k) * h_E(i,j,k)) / &
 !             (h_W(i+1,j,k) + h_E(i,j,k) + GV%H_subroundoff)
    endif

    if (marginal) then ; h_u(I,j,k) = h_marg
    else ; h_u(I,j,k) = h_avg ; endif
  enddo ; enddo ; enddo
  if (present(visc_rem_u)) then
    ! Scale back the thickness to account for the effects of viscosity and the fractional open
    ! thickness to give an appropriate non-normalized weight for each layer in determining the
    ! barotropic acceleration.
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
      h_u(I,j,k) = h_u(I,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
    enddo ; enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
      h_u(I,j,k) = h_u(I,j,k) * por_face_areaU(I,j,k)
    enddo ; enddo ; enddo
  endif

  local_open_BC = .false.
  if (associated(OBC)) local_open_BC = OBC%open_u_BCs_exist_globally
  if (local_open_BC) then
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          if (present(visc_rem_u)) then ; do k=1,nz
            do j = OBC%segment(n)%HI%jsd, OBC%segment(n)%HI%jed
              h_u(I,j,k) = h(i,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
            enddo
          enddo ; else ; do k=1,nz
            do j = OBC%segment(n)%HI%jsd, OBC%segment(n)%HI%jed
              h_u(I,j,k) = h(i,j,k) * por_face_areaU(I,j,k)
            enddo
          enddo ; endif
        else
          if (present(visc_rem_u)) then ; do k=1,nz
            do j = OBC%segment(n)%HI%jsd, OBC%segment(n)%HI%jed
              h_u(I,j,k) = h(i+1,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
            enddo
          enddo ; else ; do k=1,nz
            do j = OBC%segment(n)%HI%jsd, OBC%segment(n)%HI%jed
              h_u(I,j,k) = h(i+1,j,k) * por_face_areaU(I,j,k)
            enddo
          enddo ; endif
        endif
      endif
    enddo
  endif

end subroutine zonal_flux_thickness

!> Returns the barotropic velocity adjustment that gives the
!! desired barotropic (layer-summed) transport.
subroutine zonal_flux_adjust(u, h_in, h_W, h_E, uhbt, uh_tot_0, duhdu_tot_0, &
                             du, du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             ish, ieh, jsh, jeh, do_I_in, por_face_areaU, uh_3d, OBC)

  type(ocean_grid_type),                      intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                    intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u     !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h_in !< Layer thickness used to
                                                                    !! calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h_W  !< West edge thickness in the
                                                                    !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h_E  !< East edge thickness in the
                                                                    !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: visc_rem !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step of viscosity, and
                       !! the fraction of a time-step's worth of a barotropic acceleration that a layer
                       !! experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZIB_(G),SZJ_(G)),          intent(in)    :: uhbt !< The summed volume flux
                       !! through zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G)),          intent(in)    :: du_max_CFL  !< Maximum acceptable
                       !! value of du [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G)),          intent(in)    :: du_min_CFL  !< Minimum acceptable
                       !! value of du [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G)),          intent(in)    :: uh_tot_0    !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G)),          intent(in)    :: duhdu_tot_0 !< The partial derivative
                       !! of du_err with du at 0 adjustment [H L ~> m2 or kg m-1].
  real, dimension(SZIB_(G),SZJ_(G)),          intent(inout) :: du !<
                       !! The barotropic velocity adjustment [L T-1 ~> m s-1].
  real,                                       intent(in)    :: dt        !< Time increment [T ~> s].
  type(unit_scale_type),                      intent(in)    :: US        !< A dimensional unit scaling type.
  type(continuity_PPM_CS),                    intent(in)    :: CS        !< This module's control structure.
  integer,                                    intent(in)    :: ish, jsh  !< Start of index range.
  integer,                                    intent(in)    :: ieh, jeh  !< End of index range.
  logical, dimension(SZIB_(G),SZJ_(G)),       intent(in)    :: do_I_in   !< A logical flag
                       !! indicating which I values to work on.
  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)), intent(in)   :: por_face_areaU !< fractional open area 
                                                                              !! of U-faces [nondim].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                    optional, intent(inout) :: uh_3d !< Volume flux through zonal
                       !! faces = u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(ocean_OBC_type),             optional, pointer       :: OBC !< Open boundaries control structure.
  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    uh_aux, &      ! An auxiliary zonal volume flux [H L2 T-1 ~> m3 s-1 or kg s-1].
    duhdu          ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    uh_err, &      ! Difference between uhbt and the summed uh [H L2 T-1 ~> m3 s-1 or kg s-1].
    uh_err_best, & ! The smallest value of uh_err found so far [H L2 T-1 ~> m3 s-1 or kg s-1].
    duhdu_tot,&    ! Summed partial derivative of uh with u [H L ~> m2 or kg m-1].
    du_min, &      ! Lower limit on du correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
    du_max         ! Upper limit on du correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: u_new ! The velocity with the correction added [L T-1 ~> m s-1].
  real :: du_prev  ! The previous value of du [L T-1 ~> m s-1].
  real :: ddu      ! The change in du from the previous iteration [L T-1 ~> m s-1].
  real :: tol_eta  ! The tolerance for the current iteration [H ~> m or kg m-2].
  real :: tol_vel  ! The tolerance for velocity in the current iteration [L T-1 ~> m s-1].
  integer :: i, j, k, nz, itt, max_itts = 20
  logical :: domore, do_I(SZIB_(G),SZJ_(G))

  nz = GV%ke

  !$omp target enter data &
  !$omp   map(to: uh_3d(ish-1:ieh, :, :), do_I_in(ish-1:ieh, jsh:jeh), &
  !$omp       du_max_CFL(ish-1:ieh, jsh:jeh), du_min_CFL(ish-1:ieh, jsh:jeh), &
  !$omp       uh_tot_0(ish-1:ieh, jsh:jeh), uhbt(ish-1:ieh, jsh:jeh), &
  !$omp       duhdu_tot_0(ish-1:ieh, jsh:jeh), G, G%IareaT(ish-1:ieh+1, jsh:jeh), CS, &
  !$omp       u(ish-1:ieh, :, :), visc_rem(ish-1:ieh, :, :), h_W(ish-1:ieh+1, :, :), &
  !$omp       h_E(ish-1:ieh+1, :, :), h_in(ish-1:ieh+1, :, :), G%dy_Cu(ish-1:ieh, jsh:jeh), &
  !$omp       G%IdxT(ish-1:ieh+1, jsh:jeh), por_face_areaU(ish-1:ieh, :, :)) &
  !$omp   map(alloc: uh_aux(ish-1:ieh, :, :), duhdu(ish-1:ieh, :, :), du(ish-1:ieh, jsh:jeh), &
  !$omp       do_I(ish-1:ieh, jsh:jeh), du_max(ish-1:ieh, jsh:jeh), du_min(ish-1:ieh, jsh:jeh), &
  !$omp       duhdu_tot(ish-1:ieh, jsh:jeh), uh_err(ish-1:ieh, jsh:jeh), &
  !$omp       uh_err_best(ish-1:ieh, jsh:jeh), u_new(ish-1:ieh, :, :))

  !$omp target teams distribute parallel do collapse(3) &
  !$omp   map(from: uh_aux(ish-1:ieh, :, :), duhdu(ish-1:ieh, :, :))
  do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
    uh_aux(I,j,k) = 0.0 ; duhdu(I,j,k) = 0.0
  enddo ; enddo ; enddo

  if (present(uh_3d)) then
    !$omp target teams distribute parallel do collapse(3) &
    !$omp   map(to: uh_3d(ish-1:ieh, :, :), uh_aux(ish-1:ieh, :, :))
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
      uh_aux(i,j,k) = uh_3d(I,j,k)
    enddo ; enddo ; enddo
  endif

  !$omp target teams distribute parallel do collapse(2) &
  !$omp   map(to: do_I_in(ish-1:ieh, jsh:jeh), du_max_CFL(ish-1:ieh, jsh:jeh), &
  !$omp       du_min_CFL(ish-1:ieh, jsh:jeh), uh_tot_0(ish-1:ieh, jsh:jeh), &
  !$omp       uhbt(ish-1:ieh, jsh:jeh), duhdu_tot_0(ish-1:ieh, jsh:jeh)) &
  !$omp   map(from: du(ish-1:ieh, jsh:jeh), do_I(ish-1:ieh, jsh:jeh), du_max(ish-1:ieh, jsh:jeh), &
  !$omp       du_min(ish-1:ieh, jsh:jeh), uh_err(ish-1:ieh, jsh:jeh), duhdu_tot(ish-1:ieh, jsh:jeh), &
  !$omp       uh_err_best(ish-1:ieh, jsh:jeh))
  do j=jsh,jeh ; do I=ish-1,ieh
    du(I,j) = 0.0 ; do_I(I,j) = do_I_in(I,j)
    du_max(I,j) = du_max_CFL(I,j) ; du_min(I,j) = du_min_CFL(I,j)
    uh_err(I,j) = uh_tot_0(I,j) - uhbt(I,j) ; duhdu_tot(I,j) = duhdu_tot_0(I,j)
    uh_err_best(I,j) = abs(uh_err(I,j))
  enddo ; enddo

  do itt=1,max_itts
    select case (itt)
      case (:1) ; tol_eta = 1e-6 * CS%tol_eta
      case (2)  ; tol_eta = 1e-4 * CS%tol_eta
      case (3)  ; tol_eta = 1e-2 * CS%tol_eta
      case default ; tol_eta = CS%tol_eta
    end select
    tol_vel = CS%tol_vel

    !$omp target teams distribute parallel do collapse(2) &
    !$omp   map(to: uh_err(ish-1:ieh, jsh:jeh), du(ish-1:ieh, jsh:jeh)) &
    !$omp   map(tofrom: du_max(ish-1:ieh, jsh:jeh), du_min(ish-1:ieh, jsh:jeh), &
    !$omp       do_I(ish-1:ieh, jsh:jeh))
    do j=jsh, jeh ; do I=ish-1,ieh
      if (uh_err(I,j) > 0.0) then ; du_max(I,j) = du(I,j)
      elseif (uh_err(I,j) < 0.0) then ; du_min(I,j) = du(I,j)
      else ; do_I(I,j) = .false. ; endif
    enddo ; enddo
    domore = .false.
    !$omp target teams distribute parallel do collapse(2) &
    !$omp   private(ddu, du_prev) &
    !$omp   reduction(.or.:domore) &
    !$omp   map(to: do_I(ish-1:ieh, jsh:jeh), G, G%IareaT(ish-1:ieh+1, jsh:jeh), CS, &
    !$omp       uh_err(ish-1:ieh, jsh:jeh), duhdu_tot(ish-1:ieh, jsh:jeh), &
    !$omp       uh_err_best(ish-1:ieh, jsh:jeh), du_max(ish-1:ieh, jsh:jeh), &
    !$omp       du_min(ish-1:ieh, jsh:jeh)) &
    !$omp   map(tofrom: du(ish-1:ieh, jsh:jeh), do_I(ish-1:ieh, jsh:jeh), domore)
    do j=jsh,jeh ; do I=ish-1,ieh ; if (do_I(I,j)) then
      if ((dt * min(G%IareaT(i,j),G%IareaT(i+1,j))*abs(uh_err(I,j)) > tol_eta) .or. &
          (CS%better_iter .and. ((abs(uh_err(I,j)) > tol_vel * duhdu_tot(I,j)) .or. &
                                 (abs(uh_err(I,j)) > uh_err_best(I,j))) )) then
      !   Use Newton's method, provided it stays bounded.  Otherwise bisect
      ! the value with the appropriate bound.
        ddu = -uh_err(I,j) / duhdu_tot(I,j)
        du_prev = du(I,j)
        du(I,j) = du(I,j) + ddu
        if (abs(ddu) < 1.0e-15*abs(du(I,j))) then
          do_I(I,j) = .false. ! ddu is small enough to quit.
        elseif (ddu > 0.0) then
          if (du(I,j) >= du_max(I,j)) then
            du(I,j) = 0.5*(du_prev + du_max(I,j))
            if (du_max(I,j) - du_prev < 1.0e-15*abs(du(I,j))) do_I(I,j) = .false.
          endif
        else ! ddu < 0.0
          if (du(I,j) <= du_min(I,j)) then
            du(I,j) = 0.5*(du_prev + du_min(I,j))
            if (du_prev - du_min(I,j) < 1.0e-15*abs(du(I,j))) do_I(I,j) = .false.
          endif
        endif
        if (do_I(I,j)) domore = .true.
      else
        do_I(I,j) = .false.
      endif
    endif ; enddo ; enddo
    if (.not.domore) exit

    if ((itt < max_itts) .or. present(uh_3d)) then
      !$omp target teams distribute parallel do collapse(3) &
      !$omp   map(to: u(ish-1:ieh, :, :), du(ish-1:ieh, jsh:jeh), visc_rem(ish-1:ieh, :, :)) &
      !$omp   map(from: u_new(ish-1:ieh, :, :))
      do k=1,nz ; do j = jsh,jeh ; do I=ish-1,ieh
        u_new(I,j,k) = u(I,j,k) + du(I,j) * visc_rem(I,j,k)
      enddo ; enddo ; enddo
      call zonal_flux_layer(u_new, h_in, h_W, h_E, &
                            uh_aux, duhdu, visc_rem, &
                            dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaU, OBC)
    endif

    if (itt < max_itts) then
      !$omp target teams distribute parallel do collapse(2) &
      !$omp   map(to: uhbt(ish-1:ieh, jsh:jeh), uh_aux(ish-1:ieh, :, :), duhdu(ish-1:ieh, :, :)) &
      !$omp   map(from: duhdu_tot(ish-1:ieh, jsh:jeh), uh_err(ish-1:ieh, jsh:jeh)) &
      !$omp   map(tofrom: uh_err_best(ish-1:ieh, jsh:jeh))
      do j=jsh,jeh ; do I=ish-1,ieh
        uh_err(I,j) = -uhbt(I,j) ; duhdu_tot(i,j) = 0.0
        do k=1,nz
          uh_err(I,j) = uh_err(I,j) + uh_aux(I,j,k)
          duhdu_tot(I,j) = duhdu_tot(I,j) + duhdu(I,j,k)
        enddo
        uh_err_best(I,j) = min(uh_err_best(I,j), abs(uh_err(I,j)))
      enddo; enddo
    endif
  enddo ! itt-loop
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.

  if (present(uh_3d)) then
    !$omp target teams distribute parallel do collapse(3) &
    !$omp   map(to: uh_aux(ish-1:ieh, :, :)) &
    !$omp   map(from: uh_3d(ish-1:ieh, :, :))
    do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
      uh_3d(I,j,k) = uh_aux(I,j,k)
    enddo ; enddo ; enddo
  endif

  !$omp target exit data &
  !$omp   map(from: uh_3d(ish-1:ieh, :, :), du(ish-1:ieh, jsh:jeh)) &
  !$omp   map(release: uh_aux(ish-1:ieh, :, :), duhdu(ish-1:ieh, :, :), do_I_in(ish-1:ieh, jsh:jeh), &
  !$omp       du_max_CFL(ish-1:ieh, jsh:jeh), du_min_CFL(ish-1:ieh, jsh:jeh), &
  !$omp       uh_tot_0(ish-1:ieh, jsh:jeh), uhbt(ish-1:ieh, jsh:jeh), &
  !$omp       duhdu_tot_0(ish-1:ieh, jsh:jeh), do_I(ish-1:ieh, jsh:jeh), du_max(ish-1:ieh, jsh:jeh), &
  !$omp       du_min(ish-1:ieh, jsh:jeh), duhdu_tot(ish-1:ieh, jsh:jeh), uh_err(ish-1:ieh, jsh:jeh), &
  !$omp       uh_err_best(ish-1:ieh, jsh:jeh), G, G%IareaT(ish-1:ieh+1, jsh:jeh), CS, &
  !$omp       u(ish-1:ieh, :, :), visc_rem(ish-1:ieh, :, :), u_new(ish-1:ieh, :, :), &
  !$omp       h_W(ish-1:ieh+1, :, :), h_E(ish-1:ieh+1, :, :), h_in(ish-1:ieh+1, :, :), &
  !$omp       G%dy_Cu(ish-1:ieh, jsh:jeh), G%IdxT(ish-1:ieh+1, jsh:jeh), &
  !$omp       por_face_areaU(ish-1:ieh, :, :))

end subroutine zonal_flux_adjust


!> Sets a structure that describes the zonal barotropic volume or mass fluxes as a
!! function of barotropic flow to agree closely with the sum of the layer's transports.
subroutine set_zonal_BT_cont(u, h_in, h_W, h_E, BT_cont, uh_tot_0, duhdu_tot_0, &
                             du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             visc_rem_max, ish, ieh, jsh, jeh, do_I, por_face_areaU)
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u    !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Layer thickness used to calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_W  !< West edge thickness in the reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_E  !< East edge thickness in the reconstruction [H ~> m or kg m-2].
  type(BT_cont_type),      intent(inout) :: BT_cont !< A structure with elements
                       !! that describe the effective open face areas as a function of barotropic flow.
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: uh_tot_0    !< The summed transport with 0 adjustment 
                                                        !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: duhdu_tot_0 !< The partial derivative
                       !! of du_err with du at 0 adjustment [H L ~> m2 or kg m-1].
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: du_max_CFL  !< Maximum acceptable value of du [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: du_min_CFL  !< Minimum acceptable value of du [L T-1 ~> m s-1].
  real,                    intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS   !< This module's control structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: visc_rem !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step of viscosity, and
                       !! the fraction of a time-step's worth of a barotropic acceleration that a layer
                       !! experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: visc_rem_max !< Maximum allowable visc_rem [nondim].
  integer,                 intent(in)    :: ish      !< Start of i index range.
  integer,                 intent(in)    :: ieh      !< End of i index range.
  integer,                 intent(in)    :: jsh      !< Start of j index range.
  integer,                 intent(in)    :: jeh      !< End of j index range.
  logical, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: do_I     !< A logical flag indicating which I values to work on.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: por_face_areaU !< fractional open area of U-faces [nondim]
  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    du0, &            ! The barotropic velocity increment that gives 0 transport [L T-1 ~> m s-1].
    duL, duR, &       ! The barotropic velocity increments that give the westerly
                      ! (duL) and easterly (duR) test velocities [L T-1 ~> m s-1].
    zeros, &          ! An array of full of 0 transports [H L2 T-1 ~> m3 s-1 or kg s-1]
    du_CFL, &         ! The velocity increment that corresponds to CFL_min [L T-1 ~> m s-1].
    FAmt_L, FAmt_R, & ! The summed effective marginal face areas for the 3
    FAmt_0, &         ! test velocities [H L ~> m2 or kg m-1].
    uhtot_L, &        ! The summed transport with the westerly (uhtot_L) and
    uhtot_R           ! and easterly (uhtot_R) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    u_L, u_R, &   ! The westerly (u_L), easterly (u_R), and zero-barotropic
    u_0, &        ! transport (u_0) layer test velocities [L T-1 ~> m s-1].
    duhdu_L, &    ! The effective layer marginal face areas with the westerly
    duhdu_R, &    ! (_L), easterly (_R), and zero-barotropic (_0) test
    duhdu_0, &    ! velocities [H L ~> m2 or kg m-1].
    uh_L, uh_R, & ! The layer transports with the westerly (_L), easterly (_R),
    uh_0       ! and zero-barotropic (_0) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: FA_0    ! The effective face area with 0 barotropic transport [L H ~> m2 or kg m-1].
  real :: FA_avg  ! The average effective face area [L H ~> m2 or kg m-1], nominally given by
                  ! the realized transport divided by the barotropic velocity.
  real :: visc_rem_lim ! The larger of visc_rem and min_visc_rem [nondim]. This
                       ! limiting is necessary to keep the inverse of visc_rem
                       ! from leading to large CFL numbers.
  real :: min_visc_rem ! The smallest permitted value for visc_rem that is used
                       ! in finding the barotropic velocity that changes the
                       ! flow direction [nondim].  This is necessary to keep the inverse
                       ! of visc_rem from leading to large CFL numbers.
  real :: CFL_min ! A minimal increment in the CFL to try to ensure that the
                  ! flow is truly upwind [nondim]
  real :: Idt     ! The inverse of the time step [T-1 ~> s-1].
  logical :: domore
  integer :: i, j, k, nz

  nz = GV%ke ; Idt = 1.0 / dt
  min_visc_rem = 0.1 ; CFL_min = 1e-6

  !$omp target enter data &
  !$omp   map(to: u(ish-1:ieh, :, :), h_W(ish-1:ieh+1, :, :), h_E(ish-1:ieh+1, :, :), &
  !$omp       h_in(ish-1:ieh+1, :, :), uh_tot_0(ish-1:ieh, jsh:jeh), duhdu_tot_0(ish-1:ieh, jsh:jeh), &
  !$omp       G, G%IareaT(ish-1:ieh+1, jsh:jeh), G%dy_Cu(ish-1:ieh, jsh:jeh), &
  !$omp       G%IdxT(ish-1:ieh+1, jsh:jeh), G%dxCu(ish-1:ieh, jsh:jeh), BT_cont, &
  !$omp       du_max_CFL(ish-1:ieh, jsh:jeh), du_min_CFL(ish-1:ieh, jsh:jeh), CS, &
  !$omp       visc_rem(ish-1:ieh, :, :), visc_rem_max(ish-1:ieh, jsh:jeh), do_I(ish-1:ieh, jsh:jeh), &
  !$omp       por_face_areaU(ish-1:ieh, :, :)) &
  !$omp   map(alloc: zeros(ish-1:ieh, jsh:jeh), du0(ish-1:ieh, jsh:jeh), duL(ish-1:ieh, jsh:jeh), &
  !$omp       duR(ish-1:ieh, jsh:jeh), du_CFL(ish-1:ieh, jsh:jeh), FAmt_L(ish-1:ieh, jsh:jeh), &
  !$omp       FAmT_R(ish-1:ieh, jsh:jeh), FAmt_0(ish-1:ieh, jsh:jeh), uhtot_L(ish-1:ieh, jsh:jeh), &
  !$omp       uhtot_R(ish-1:ieh, jsh:jeh), u_L(ish-1:ieh, :, :), u_R(ish-1:ieh, :, :), &
  !$omp       u_0(ish-1:ieh, :, :), duhdu_L(ish-1:ieh, :, :), duhdu_R(ish-1:ieh, :, :), &
  !$omp       duhdu_0(ish-1:ieh, :, :), uh_L(ish-1:ieh, :, :), uh_R(ish-1:ieh, :, :), &
  !$omp       uh_0(ish-1:ieh, :, :), BT_cont%FA_u_W0(ish-1:ieh, jsh:jeh), &
  !$omp       BT_cont%FA_u_WW(ish-1:ieh, jsh:jeh), BT_cont%FA_u_E0(ish-1:ieh, jsh:jeh), &
  !$omp       BT_cont%FA_u_EE(ish-1:ieh, jsh:jeh), BT_cont%uBT_WW(ish-1:ieh, jsh:jeh), &
  !$omp       BT_cont%uBT_EE(ish-1:ieh, jsh:jeh))

  ! Diagnose the zero-transport correction, du0.
  !$omp target teams distribute parallel do collapse(2) &
  !$omp   map(from: zeros(ish-1:ieh, jsh:jeh))
  do j=jsh,jeh ; do I=ish-1,ieh ; zeros(I, j) = 0.0 ; enddo ; enddo
  call zonal_flux_adjust(u, h_in, h_W, h_E, zeros, uh_tot_0, duhdu_tot_0, du0, &
                         du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                         ish, ieh, jsh, jeh, do_I, por_face_areaU)

  ! Determine the westerly- and easterly- fluxes.  Choose a sufficiently
  ! negative velocity correction for the easterly-flux, and a sufficiently
  ! positive correction for the westerly-flux.
  domore = .false.
  !$omp target teams distribute parallel do collapse(2) &
  !$omp   reduction(.or.:domore) &
  !$omp   map(to: do_I(ish-1:ieh, jsh:jeh), G, G%dxCu(ish-1:ieh, jsh:jeh), du0(ish-1:ieh, jsh:jeh)) &
  !$omp   map(from: du_CFL(ish-1:ieh, jsh:jeh), duR(ish-1:ieh, jsh:jeh), duL(ish-1:ieh, jsh:jeh), &
  !$omp       FAmt_L(ish-1:ieh, jsh:jeh), FAmt_R(ish-1:ieh, jsh:jeh), FAmt_0(ish-1:ieh, jsh:jeh), &
  !$omp       uhtot_L(ish-1:ieh, jsh:jeh), uhtot_R(ish-1:ieh, jsh:jeh))
  do j=jsh,jeh ; do I=ish-1,ieh
    if (do_I(I,j)) domore = .true.
    du_CFL(I,j) = (CFL_min * Idt) * G%dxCu(I,j)
    duR(I,j) = min(0.0,du0(I,j) - du_CFL(I,j))
    duL(I,j) = max(0.0,du0(I,j) + du_CFL(I,j))
    FAmt_L(I,j) = 0.0 ; FAmt_R(I,j) = 0.0 ; FAmt_0(I,j) = 0.0
    uhtot_L(I,j) = 0.0 ; uhtot_R(I,j) = 0.0
  enddo ; enddo

  ! short circuit if none should be updated
  if (.not.domore) then
    !$omp target teams distribute parallel do collapse(2) &
    !$omp   map(alloc: BT_cont) &
    !$omp   map(from: BT_cont%FA_u_W0(ish-1:ieh, jsh:jeh), BT_cont%FA_u_WW(ish-1:ieh, jsh:jeh), &
    !$omp       BT_cont%FA_u_E0(ish-1:ieh, jsh:jeh), BT_cont%FA_u_EE(ish-1:ieh, jsh:jeh), &
    !$omp       BT_cont%uBT_WW(ish-1:ieh, jsh:jeh), BT_cont%uBT_EE(ish-1:ieh, jsh:jeh))
    do j=jsh,jeh ; do I=ish-1,ieh
      BT_cont%FA_u_W0(I,j) = 0.0 ; BT_cont%FA_u_WW(I,j) = 0.0
      BT_cont%FA_u_E0(I,j) = 0.0 ; BT_cont%FA_u_EE(I,j) = 0.0
      BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
    enddo ; enddo
    return
  endif

  !$omp target teams distribute parallel do collapse(2) &
  !$omp   private(visc_rem_lim) &
  !$omp   map(to: do_I(ish-1:ieh, jsh:jeh), visc_rem(ish-1:ieh, :, :), &
  !$omp       visc_rem_max(ish-1:ieh, jsh:jeh), u(ish-1:ieh, :, :), du_CFL(ish-1:ieh, jsh:jeh)) &
  !$omp   map(tofrom: duR(ish-1:ieh, jsh:jeh), duL(ish-1:ieh, jsh:jeh))
  do j=jsh,jeh ; do I=ish-1,ieh ; if (do_I(I,j)) then
    do k=1,nz ! k-loop is serialised
      visc_rem_lim = max(visc_rem(I,j,k), min_visc_rem*visc_rem_max(I,j))
      if (visc_rem_lim > 0.0) then ! This is almost always true for ocean points.
        if (u(I,j,k) + duR(I,j)*visc_rem_lim > -du_CFL(I,j)*visc_rem(I,j,k)) &
          duR(I,j) = -(u(I,j,k) + du_CFL(I,j)*visc_rem(I,j,k)) / visc_rem_lim
        if (u(I,j,k) + duL(I,j)*visc_rem_lim < du_CFL(I,j)*visc_rem(I,j,k)) &
          duL(I,j) = -(u(I,j,k) - du_CFL(I,j)*visc_rem(I,j,k)) / visc_rem_lim
      endif
    enddo
  endif ; enddo ; enddo

  !$omp target teams distribute parallel do collapse(3) &
  !$omp   map(to: do_I(ish-1:ieh, jsh:jeh), u(ish-1:ieh, :, :), duL(ish-1:ieh, jsh:jeh), &
  !$Omp       duR(ish-1:ieh, jsh:jeh), du0(ish-1:ieh, jsh:jeh), visc_rem(ish-1:ieh, :, :)) &
  !$Omp   map(from: u_L(ish-1:ieh, :, :), u_R(ish-1:ieh, :, :), u_0(ish-1:ieh, :, :))
  do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh ; if (do_I(I,j)) then
    u_L(I,j,k) = u(I,j,k) + duL(I,j) * visc_rem(I,j,k)
    u_R(I,j,k) = u(I,j,k) + duR(I,j) * visc_rem(I,j,k)
    u_0(I,j,k) = u(I,j,k) + du0(I,j) * visc_rem(I,j,k)
  endif ; enddo ; enddo ; enddo

  call zonal_flux_layer(u_0, h_in, h_W, h_E, uh_0, duhdu_0, &
                        visc_rem, dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaU)
  call zonal_flux_layer(u_L, h_in, h_W, h_E, uh_L, duhdu_L, &
                        visc_rem, dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaU)
  call zonal_flux_layer(u_R, h_in, h_W, h_E, uh_R, duhdu_R, &
                        visc_rem, dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaU)

  !$omp target teams distribute parallel do collapse(2) &
  !$omp   map(to: do_I(ish-1:ieh, jsh:jeh), duhdu_0(ish-1:ieh, :, :), duhdu_L(ish-1:ieh, :, :), &
  !$omp       duhdu_R(ish-1:ieh, :, :), uh_L(ish-1:ieh, :, :), uh_R(ish-1:ieh, :, :)) &
  !$omp   map(tofrom: FAmt_0(ish-1:ieh, jsh:jeh), FAmt_L(ish-1:ieh, jsh:jeh), &
  !$omp       FAmt_R(ish-1:ieh, jsh:jeh), uhtot_L(ish-1:ieh, jsh:jeh), uhtot_R(ish-1:ieh, jsh:jeh))
  do j=jsh,jeh ; do I=ish-1,ieh ; if (do_I(I,j)) then
    do k=1,nz
      FAmt_0(I,j) = FAmt_0(I,j) + duhdu_0(I,j,k)
      FAmt_L(I,j) = FAmt_L(I,j) + duhdu_L(I,j,k)
      FAmt_R(I,j) = FAmt_R(I,j) + duhdu_R(I,j,k)
      uhtot_L(I,j) = uhtot_L(I,j) + uh_L(I,j,k)
      uhtot_R(I,j) = uhtot_R(I,j) + uh_R(I,j,k)
    enddo
  endif ; enddo ; enddo

  !$omp target teams distribute parallel do collapse(2) &
  !$omp   private(FA_0, FA_avg) &
  !$omp   map(to: do_I(ish-1:ieh, jsh:jeh), FAmt_0(ish-1:ieh, jsh:jeh), FAmt_L(ish-1:ieh, jsh:jeh), &
  !$omp       FAmt_R(ish-1:ieh, jsh:jeh), du0(ish-1:ieh, jsh:jeh), duL(ish-1:ieh, jsh:jeh), &
  !$omp       duR(ish-1:ieh, jsh:jeh), uhtot_L(ish-1:ieh, jsh:jeh), uhtot_R(ish-1:ieh, jsh:jeh), &
  !$omp       BT_cont) &
  !$omp   map(from: BT_cont%FA_u_W0(ish-1:ieh, jsh:jeh), BT_cont%FA_u_WW(ish-1:ieh, jsh:jeh), &
  !$omp       BT_cont%FA_u_E0(ish-1:ieh, jsh:jeh), BT_cont%FA_u_EE(ish-1:ieh, jsh:jeh), &
  !$omp       BT_cont%uBT_WW(ish-1:ieh, jsh:jeh), BT_cont%uBT_EE(ish-1:ieh, jsh:jeh))
  do j=jsh,jeh ; do I=ish-1,ieh ; if (do_I(I,j)) then
    FA_0 = FAmt_0(I,j) ; FA_avg = FAmt_0(I,j)
    if ((duL(I,j) - du0(I,j)) /= 0.0) &
      FA_avg = uhtot_L(I,j) / (duL(I,j) - du0(I,j))
    if (FA_avg > max(FA_0, FAmt_L(I,j))) then ; FA_avg = max(FA_0, FAmt_L(I,j))
    elseif (FA_avg < min(FA_0, FAmt_L(I,j))) then ; FA_0 = FA_avg ; endif

    BT_cont%FA_u_W0(I,j) = FA_0 ; BT_cont%FA_u_WW(I,j) = FAmt_L(I,j)
    if (abs(FA_0-FAmt_L(I,j)) <= 1e-12*FA_0) then ; BT_cont%uBT_WW(I,j) = 0.0 ; else
      BT_cont%uBT_WW(I,j) = (1.5 * (duL(I,j) - du0(I,j))) * &
                            ((FAmt_L(I,j) - FA_avg) / (FAmt_L(I,j) - FA_0))
    endif

    FA_0 = FAmt_0(I,j) ; FA_avg = FAmt_0(I,j)
    if ((duR(I,j) - du0(I,j)) /= 0.0) &
      FA_avg = uhtot_R(I,j) / (duR(I,j) - du0(I,j))
    if (FA_avg > max(FA_0, FAmt_R(I,j))) then ; FA_avg = max(FA_0, FAmt_R(I,j))
    elseif (FA_avg < min(FA_0, FAmt_R(I,j))) then ; FA_0 = FA_avg ; endif

    BT_cont%FA_u_E0(I,j) = FA_0 ; BT_cont%FA_u_EE(I,j) = FAmt_R(I,j)
    if (abs(FAmt_R(I,j) - FA_0) <= 1e-12*FA_0) then ; BT_cont%uBT_EE(I,j) = 0.0 ; else
      BT_cont%uBT_EE(I,j) = (1.5 * (duR(I,j) - du0(I,j))) * &
                            ((FAmt_R(I,j) - FA_avg) / (FAmt_R(I,j) - FA_0))
    endif
  else
    BT_cont%FA_u_W0(I,j) = 0.0 ; BT_cont%FA_u_WW(I,j) = 0.0
    BT_cont%FA_u_E0(I,j) = 0.0 ; BT_cont%FA_u_EE(I,j) = 0.0
    BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
  endif ; enddo ; enddo

  !$omp target update from(BT_cont%FA_u_W0(ish-1:ieh, jsh:jeh), BT_cont%FA_u_WW(ish-1:ieh, jsh:jeh), &
  !$omp       BT_cont%FA_u_E0(ish-1:ieh, jsh:jeh), BT_cont%FA_u_EE(ish-1:ieh, jsh:jeh), &
  !$omp       BT_cont%uBT_WW(ish-1:ieh, jsh:jeh), BT_cont%uBT_EE(ish-1:ieh, jsh:jeh))
  !$omp target exit data &
  !$omp   map(from: BT_cont%FA_u_W0(ish-1:ieh, jsh:jeh), BT_cont%FA_u_WW(ish-1:ieh, jsh:jeh), &
  !$omp       BT_cont%FA_u_E0(ish-1:ieh, jsh:jeh), BT_cont%FA_u_EE(ish-1:ieh, jsh:jeh), &
  !$omp       BT_cont%uBT_WW(ish-1:ieh, jsh:jeh), BT_cont%uBT_EE(ish-1:ieh, jsh:jeh)) &
  !$omp   map(release: zeros(ish-1:ieh, jsh:jeh), u(ish-1:ieh, :, :), h_W(ish-1:ieh+1, :, :), &
  !$omp       h_E(ish-1:ieh+1, :, :), h_in(ish-1:ieh+1, :, :), uh_tot_0(ish-1:ieh, jsh:jeh), &
  !$omp       duhdu_tot_0(ish-1:ieh, jsh:jeh), du0(ish-1:ieh, jsh:jeh), duL(ish-1:ieh, jsh:jeh), &
  !$omp       duR(ish-1:ieh, jsh:jeh), du_CFL(ish-1:ieh, jsh:jeh), FAmt_L(ish-1:ieh, jsh:jeh), &
  !$omp       FAmT_R(ish-1:ieh, jsh:jeh), FAmt_0(ish-1:ieh, jsh:jeh), uhtot_L(ish-1:ieh, jsh:jeh), &
  !$omp       uhtot_R(ish-1:ieh, jsh:jeh), u_L(ish-1:ieh, :, :), u_R(ish-1:ieh, :, :), &
  !$omp       u_0(ish-1:ieh, :, :), duhdu_L(ish-1:ieh, :, :), duhdu_R(ish-1:ieh, :, :), &
  !$omp       duhdu_0(ish-1:ieh, :, :), uh_L(ish-1:ieh, :, :), uh_R(ish-1:ieh, :, :), &
  !$omp       uh_0(ish-1:ieh, :, :), G, G%IareaT(ish-1:ieh+1, jsh:jeh), G%dy_Cu(ish-1:ieh, jsh:jeh), &
  !$omp       G%IdxT(ish-1:ieh+1, jsh:jeh), G%dxCu(ish-1:ieh, jsh:jeh), BT_cont, &
  !$omp       du_max_CFL(ish-1:ieh, jsh:jeh), du_min_CFL(ish-1:ieh, jsh:jeh), CS, &
  !$omp       visc_rem(ish-1:ieh, :, :), visc_rem_max(ish-1:ieh, jsh:jeh), do_I(ish-1:ieh, jsh:jeh), &
  !$omp       por_face_areaU(ish-1:ieh, :, :))

end subroutine set_zonal_BT_cont

!> Calculates the mass or volume fluxes through the meridional faces, and other related quantities.
subroutine meridional_mass_flux(v, h_in, h_S, h_N, vh, dt, G, GV, US, CS, OBC, por_face_areaV, &
                                LB_in, vhbt, visc_rem_v, v_cor, BT_cont, dv_cor)
  type(ocean_grid_type),                      intent(in)  :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: v    !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_in !< Layer thickness used to
                                                                  !! calculate fluxes [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_S  !< South edge thickness in the
                                                                  !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_N  !< North edge thickness in the
                                                                  !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: vh   !< Volume flux through meridional
                                                                  !! faces = v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,                                       intent(in)  :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                      intent(in)  :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),                    intent(in)  :: CS   !< This module's control structure.G
  type(ocean_OBC_type),                       pointer     :: OBC  !< Open boundary condition type
                                                                  !! specifies whether, where, and what
                                                                  !! open boundary conditions are used.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)),  intent(in)  :: por_face_areaV !< fractional open area of V-faces [nondim]
  type(cont_loop_bounds_type),      optional, intent(in)  :: LB_in !< Loop bounds structure.
  real, dimension(SZI_(G),SZJB_(G)), optional, intent(in) :: vhbt !< The summed volume flux through meridional
                                                                  !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    optional, intent(in)  :: visc_rem_v !< Both the fraction of the momentum
                                   !! originally in a layer that remains after a time-step of viscosity,
                                   !! and the fraction of a time-step's worth of a barotropic acceleration
                                   !! that a layer experiences after viscosity is applied [nondim].
                                   !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    optional, intent(out) :: v_cor
                                   !< The meridional velocities (v with a barotropic correction)
                                   !! that give vhbt as the depth-integrated transport [L T-1 ~> m s-1].
  type(BT_cont_type),               optional, pointer     :: BT_cont !< A structure with elements that describe
                                   !! the effective open face areas as a function of barotropic flow.
  real, dimension(SZI_(G),SZJB_(G)), &
                                    optional, intent(out)   :: dv_cor !< The meridional velocity increments from v
                                                                  !! that give vhbt as the depth-integrated
                                                                  !! transports [L T-1 ~> m s-1].

  ! Local variables
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    dvhdv         ! Partial derivative of vh with v [H L ~> m2 or kg m-1].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    dv, &         ! Corrective barotropic change in the velocity to give vhbt [L T-1 ~> m s-1].
    dv_min_CFL, & ! Lower limit on dv correction to avoid CFL violations [L T-1 ~> m s-1]
    dv_max_CFL, & ! Upper limit on dv correction to avoid CFL violations [L T-1 ~> m s-1]
    dvhdv_tot_0, & ! Summed partial derivative of vh with v [H L ~> m2 or kg m-1].
    vh_tot_0, &   ! Summed transport with no barotropic correction [H L2 T-1 ~> m3 s-1 or kg s-1].
    visc_rem_max  ! The column maximum of visc_rem [nondim]
  logical, dimension(SZI_(G),SZJB_(G)) :: do_I
  real, dimension(SZI_(G),SZJB_(G)) :: FAvi  ! A list of sums of meridional face areas [H L ~> m2 or kg m-1].
  real :: FA_v    ! A sum of meridional face areas [H L ~> m2 or kg m-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: visc_rem_v_tmp ! A copy of visc_rem_v or an array of 1's [nondim]
  real :: I_vrm   ! 1.0 / visc_rem_max [nondim]
  real :: CFL_dt  ! The maximum CFL ratio of the adjusted velocities divided by
                  ! the time step [T-1 ~> s-1].
  real :: I_dt    ! 1.0 / dt [T-1 ~> s-1].
  real :: dv_lim  ! The velocity change that give a relative CFL of 1 [L T-1 ~> m s-1].
  real :: dy_N, dy_S ! Effective y-grid spacings to the north and south [L ~> m].
  type(cont_loop_bounds_type) :: LB
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz
  integer :: l_seg ! The OBC segment number
  logical :: use_visc_rem, set_BT_cont
  logical :: local_specified_BC, local_Flather_OBC, local_open_BC, any_simple_OBC  ! OBC-related logicals
  logical :: simple_OBC_pt(SZI_(G),SZJ_(G))  ! Indicates points in a row with specified transport OBCs
  type(OBC_segment_type), pointer :: segment => NULL()

  call cpu_clock_begin(id_clock_correct)

  use_visc_rem = present(visc_rem_v)

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))

  local_specified_BC = .false. ; local_Flather_OBC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_v_BCs_exist_globally
    local_Flather_OBC = OBC%Flather_v_BCs_exist_globally
    local_open_BC = OBC%open_v_BCs_exist_globally
  endif ; endif

  if (present(dv_cor)) dv_cor(:,:) = 0.0

  if (present(LB_in)) then
    LB = LB_in
  else
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
  endif
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  CFL_dt = CS%CFL_limit_adjust / dt
  I_dt = 1.0 / dt
  if (CS%aggress_adjust) CFL_dt = I_dt

  ! a better solution is needed
  if (.not.use_visc_rem) then
    visc_rem_v_tmp(:, :, :) = 1.0
  else
    visc_rem_v_tmp(:, :, :) = visc_rem_v(:, :, :)
  endif
  do j=jsh-1,jeh ; do i=ish,ieh ; do_I(i,j) = .true. ; enddo ; enddo
  ! This sets vh and dvhdv.
  call merid_flux_layer(v, h_in, h_S, h_N, &
                        vh, dvhdv, visc_rem_v_tmp, &
                        dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaV, OBC)
  ! untested
  if (local_specified_BC) then
    do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh ; if (OBC%segnum_v(i,J) /= 0) then
      l_seg = abs(OBC%segnum_v(i,J))
      if (OBC%segment(l_seg)%specified) vh(i,J,k) = OBC%segment(l_seg)%normal_trans(i,J,k)
    endif ; enddo ; enddo ; enddo
  endif
  if (present(vhbt) .or. set_BT_cont) then
    if (use_visc_rem .and. CS%use_visc_rem_max) then
      visc_rem_max(:, :) = 0.0
      do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh
        visc_rem_max(i, j) = max(visc_rem_max(i, j), visc_rem_v(i,j,k))
      enddo ; enddo ; enddo
    else
      visc_rem_max(:, :) = 1.0
    endif
    !   Set limits on dv that will keep the CFL number between -1 and 1.
    ! This should be adequate to keep the root bracketed in all cases.
    do j=jsh-1,jeh ; do i=ish,ieh
      I_vrm = 0.0
      if (visc_rem_max(i,j) > 0.0) I_vrm = 1.0 / visc_rem_max(i,j)
      if (CS%vol_CFL) then
        dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(i,J), 1000.0*G%dyT(i,j))
        dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(i,J), 1000.0*G%dyT(i,j+1))
      else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
      dv_max_CFL(i,j) = 2.0 * (CFL_dt * dy_S) * I_vrm
      dv_min_CFL(i,j) = -2.0 * (CFL_dt * dy_N) * I_vrm
      vh_tot_0(i,j) = 0.0 ; dvhdv_tot_0(i,j) = 0.0
    enddo ; enddo
    do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh
      dvhdv_tot_0(i,j) = dvhdv_tot_0(i,j) + dvhdv(i,j,k)
      vh_tot_0(i,j) = vh_tot_0(i,j) + vh(i,J,k)
    enddo ; enddo ; enddo
    if (use_visc_rem) then
      if (CS%aggress_adjust) then
        ! untested
        do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh
          if (CS%vol_CFL) then
            dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(I,j), 1000.0*G%dyT(i,j))
            dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(I,j), 1000.0*G%dyT(i,j+1))
          else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
          dv_lim = 0.499*((dy_S*I_dt - v(i,J,k)) + MIN(0.0,v(i,J-1,k)))
          if (dv_max_CFL(i,j) * visc_rem_v_tmp(I,j,k) > dv_lim) &
            dv_max_CFL(i,j) = dv_lim / visc_rem_v_tmp(I,j,k)

          dv_lim = 0.499*((-dy_N*CFL_dt - v(i,J,k)) + MAX(0.0,v(i,J+1,k)))
          if (dv_min_CFL(i,j) * visc_rem_v_tmp(I,j,k) < dv_lim) &
            dv_min_CFL(i,j) = dv_lim / visc_rem_v_tmp(I,j,k)
        enddo ; enddo ; enddo
      else
        do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh
          if (CS%vol_CFL) then
            dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(I,j), 1000.0*G%dyT(i,j))
            dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(I,j), 1000.0*G%dyT(i,j+1))
          else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
          if (dv_max_CFL(i,j) * visc_rem_v_tmp(I,j,k) > dy_S*CFL_dt - v(i,J,k)*G%mask2dCv(i,J)) &
            dv_max_CFL(i,j) = (dy_S*CFL_dt - v(i,J,k)) / visc_rem_v_tmp(I,j,k)
          if (dv_min_CFL(i,j) * visc_rem_v_tmp(I,j,k) < -dy_N*CFL_dt - v(i,J,k)*G%mask2dCv(i,J)) &
            dv_min_CFL(i,j) = -(dy_N*CFL_dt + v(i,J,k)) / visc_rem_v_tmp(I,j,k)
        enddo ; enddo ; enddo
      endif
    else
      if (CS%aggress_adjust) then
        ! untested
        do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh
          if (CS%vol_CFL) then
            dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(I,j), 1000.0*G%dyT(i,j))
            dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(I,j), 1000.0*G%dyT(i,j+1))
          else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
          dv_max_CFL(i,j) = min(dv_max_CFL(i,j), 0.499 * &
                      ((dy_S*I_dt - v(i,J,k)) + MIN(0.0,v(i,J-1,k))) )
          dv_min_CFL(i,j) = max(dv_min_CFL(i,j), 0.499 * &
                      ((-dy_N*I_dt - v(i,J,k)) + MAX(0.0,v(i,J+1,k))) )
        enddo ; enddo ; enddo
      else
        do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh
          if (CS%vol_CFL) then
            dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(I,j), 1000.0*G%dyT(i,j))
            dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(I,j), 1000.0*G%dyT(i,j+1))
          else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
          dv_max_CFL(i,j) = min(dv_max_CFL(i,j), dy_S*CFL_dt - v(i,J,k))
          dv_min_CFL(i,j) = max(dv_min_CFL(i,j), -(dy_N*CFL_dt + v(i,J,k)))
        enddo ; enddo ; enddo
      endif
    endif
    do j=jsh-1,jeh ; do i=ish,ieh
      dv_max_CFL(i,j) = max(dv_max_CFL(i,j),0.0)
      dv_min_CFL(i,j) = min(dv_min_CFL(i,j),0.0)
    enddo ; enddo

    ! untested
    any_simple_OBC = .false.
    if (present(vhbt) .or. set_BT_cont) then
      if (local_specified_BC .or. local_Flather_OBC) then ; do j=jsh-1,jeh ; do i=ish,ieh
        l_seg = abs(OBC%segnum_v(i,J))

        ! Avoid reconciling barotropic/baroclinic transports if transport is specified
        simple_OBC_pt(i,j) = .false.
        if (l_seg /= 0) simple_OBC_pt(i,j) = OBC%segment(l_seg)%specified
        do_I(i,j) = .not.simple_OBC_pt(i,j)
        any_simple_OBC = any_simple_OBC .or. simple_OBC_pt(i,j)
      enddo ; enddo ; else ; do j=jsh-1,jeh ; do i=ish,ieh
        do_I(i,j) = .true.
      enddo ; enddo ; endif ! local_specified_BC .or. local_Flather_OBC
    endif ! present(vhbt) .or. set_BT_cont (redundant?)
    if (present(vhbt)) then
      ! Find dv and vh.
      call meridional_flux_adjust(v, h_in, h_S, h_N, vhbt, vh_tot_0, dvhdv_tot_0, dv, &
                             dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem_v_tmp, &
                             ish, ieh, jsh, jeh, do_I, por_face_areaV, vh, OBC=OBC)

      if (present(v_cor)) then
        do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh
          v_cor(i,j,k) = v(i,j,k) + dv(i,j) * visc_rem_v_tmp(i,j,k)
        enddo ; enddo ; enddo
        if (any_simple_OBC) then
            ! untested
            do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh ; if (simple_OBC_pt(i,j)) then
              v_cor(i,J,k) = OBC%segment(abs(OBC%segnum_v(i,J)))%normal_vel(i,J,k)
            endif ; enddo ; enddo ; enddo
        endif
      endif ! v-corrected

      if (present(dv_cor)) then
        do j=jsh-1,jeh ; do i=ish,ieh ; dv_cor(i,J) = dv(i,j) ; enddo ; enddo
      endif ! dv-corrected
    endif
    
    if (set_BT_cont) then
      call set_merid_BT_cont(v, h_in, h_S, h_N, BT_cont, vh_tot_0, dvhdv_tot_0, &
                             dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem_v_tmp, &
                             visc_rem_max, ish, ieh, jsh, jeh, do_I, por_face_areaV)
        
      if (any_simple_OBC) then
        ! untested
        do j=jsh-1,jeh ; do i=ish,ieh
          if (simple_OBC_pt(i,j)) FAvi(i,j) = GV%H_subroundoff*G%dx_Cv(i,J)
        enddo ; enddo
        ! NOTE: simple_OBC_pt(i,j) should prevent access to segment OBC_NONE
        do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh ; if (simple_OBC_pt(i,j)) then
          segment => OBC%segment(abs(OBC%segnum_v(i,J)))
          if ((abs(segment%normal_vel(i,J,k)) > 0.0) .and. (segment%specified)) &
            FAvi(i,j) = FAvi(i,j) + segment%normal_trans(i,J,k) / segment%normal_vel(i,J,k)
        endif ; enddo ; enddo ; enddo
        do j=jsh-1,jeh ; do i=ish,ieh ; if (simple_OBC_pt(i,j)) then
          BT_cont%FA_v_S0(i,J) = FAvi(i,j) ; BT_cont%FA_v_N0(i,J) = FAvi(i,j)
          BT_cont%FA_v_SS(i,J) = FAvi(i,j) ; BT_cont%FA_v_NN(i,J) = FAvi(i,j)
          BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
        endif ; enddo ; enddo
      endif ! any_simple_OBC
    endif ! set_BT_cont
  endif ! present(vhbt) or set_BT_cont

  if (local_open_BC .and. set_BT_cont) then
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          do i = OBC%segment(n)%HI%Isd, OBC%segment(n)%HI%Ied
            FA_v = 0.0
            do k=1,nz ; FA_v = FA_v + h_in(i,j,k)*(G%dx_Cv(i,J)*por_face_areaV(i,J,k)) ; enddo
            BT_cont%FA_v_S0(i,J) = FA_v ; BT_cont%FA_v_N0(i,J) = FA_v
            BT_cont%FA_v_SS(i,J) = FA_v ; BT_cont%FA_v_NN(i,J) = FA_v
            BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
          enddo
        else
          do i = OBC%segment(n)%HI%Isd, OBC%segment(n)%HI%Ied
            FA_v = 0.0
            do k=1,nz ; FA_v = FA_v + h_in(i,j+1,k)*(G%dx_Cv(i,J)*por_face_areaV(i,J,k)) ; enddo
            BT_cont%FA_v_S0(i,J) = FA_v ; BT_cont%FA_v_N0(i,J) = FA_v
            BT_cont%FA_v_SS(i,J) = FA_v ; BT_cont%FA_v_NN(i,J) = FA_v
            BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
          enddo
        endif
      endif
    enddo
  endif

  if (set_BT_cont) then ; if (allocated(BT_cont%h_v)) then
    if (present(v_cor)) then
      call meridional_flux_thickness(v_cor, h_in, h_S, h_N, BT_cont%h_v, dt, G, GV, US, LB, &
                                    CS%vol_CFL, CS%marginal_faces, OBC, por_face_areaV, visc_rem_v)
    else
      call meridional_flux_thickness(v, h_in, h_S, h_N, BT_cont%h_v, dt, G, GV, US, LB, &
                                    CS%vol_CFL, CS%marginal_faces, OBC, por_face_areaV, visc_rem_v)
    endif
  endif ; endif

  call cpu_clock_end(id_clock_correct)

end subroutine meridional_mass_flux


!> Calculates the vertically integrated mass or volume fluxes through the meridional faces.
subroutine meridional_BT_mass_flux(v, h_in, h_S, h_N, vhbt, dt, G, GV, US, CS, OBC, por_face_areaV, LB_in)
  type(ocean_grid_type),                      intent(in)  :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: v    !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_in !< Layer thickness used to
                                                                  !! calculate fluxes [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_S  !< Southern edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_N  !< Northern edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G)),          intent(out) :: vhbt !< The summed volume flux through meridional
                                                                  !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                                       intent(in)  :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                      intent(in)  :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),                    intent(in)  :: CS   !< This module's control structure.G
  type(ocean_OBC_type),                       pointer     :: OBC  !< Open boundary condition type
                                                                  !! specifies whether, where, and what
                                                                  !! open boundary conditions are used.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)),  intent(in)  :: por_face_areaV !< fractional open area of V-faces [nondim]
  type(cont_loop_bounds_type),      optional, intent(in)  :: LB_in !< Loop bounds structure.

  ! Local variables
  real :: vh(SZI_(G),SZJB_(G),SZK_(GV))      ! Volume flux through meridional faces = v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
  real ::  dvhdv(SZI_(G),SZJB_(G),SZK_(GV))  ! Partial derivative of vh with v [H L ~> m2 or kg m-1].
  logical, dimension(SZI_(G),SZJB_(G)) :: do_I
  real :: ones(SZI_(G),SZJB_(G),SZK_(GV))    ! An array of 1's [nondim]
  integer :: i, j, k, ish, ieh, jsh, jeh, nz, l_seg
  logical :: local_specified_BC, OBC_in_row(SZJB_(G))

  call cpu_clock_begin(id_clock_correct)

  local_specified_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_v_BCs_exist_globally
  endif ; endif

  if (present(LB_in)) then
    ish = LB_in%ish ; ieh = LB_in%ieh ; jsh = LB_in%jsh ; jeh = LB_in%jeh ; nz = GV%ke
  else
    ish = G%isc ; ieh = G%iec ; jsh = G%jsc ; jeh = G%jec ; nz = GV%ke
  endif

  ones(:,:,:) = 1.0 ; do_I(:,:) = .true.

  vhbt(:,:) = 0.0

  ! Determining whether there are any OBC points outside of the k-loop should be more efficient.
  OBC_in_row(:) = .false.
  if (local_specified_BC) then
    do j=jsh-1,jeh ; do i=ish,ieh ; if (OBC%segnum_v(i,J) /= 0) then
      if (OBC%segment(abs(OBC%segnum_v(i,J)))%specified) OBC_in_row(j) = .true.
    endif ; enddo ; enddo
  endif
  
  ! This sets vh and dvhdv.
  call merid_flux_layer(v, h_in, h_S, h_N, vh, dvhdv, ones, &
                        dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaV, OBC)
  
  do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh
    if (OBC_in_row(j) .and. OBC%segnum_v(i,J) /= 0) then
      l_seg = abs(OBC%segnum_v(i,J))
      if (OBC%segment(l_seg)%specified) vh(i,j,k) = OBC%segment(l_seg)%normal_trans(i,J,k)
    endif
  enddo ; enddo ; enddo


  ! Accumulate the barotropic transport.
  do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
    vhbt(i,J) = vhbt(i,J) + vh(i,J,k)
  enddo ; enddo ; enddo

  call cpu_clock_end(id_clock_correct)

end subroutine meridional_BT_mass_flux


!> Evaluates the meridional mass or volume fluxes in a layer.
subroutine merid_flux_layer(v, h, h_S, h_N, vh, dvhdv, visc_rem, dt, G, GV, US, &
                            ish, ieh, jsh, jeh, nz, do_I, vol_CFL, por_face_areaV, OBC)
  type(ocean_grid_type),                      intent(in)    :: G        !< Ocean's grid structure.
  type(verticalGrid_type),                    intent(in)    :: GV       !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v        !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: visc_rem !< Both the fraction of the
         !! momentum originally in a layer that remains after a time-step
         !! of viscosity, and the fraction of a time-step's worth of a barotropic
         !! acceleration that a layer experiences after viscosity is applied [nondim].
         !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h       !< Layer thickness used to calculate fluxes,
                                                                       !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h_S     !< South edge thickness in the reconstruction
                                                                       !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h_N     !< North edge thickness in the reconstruction
                                                                       !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vh      !< Meridional mass or volume transport
                                                                       !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: dvhdv   !< Partial derivative of vh with v
                                                                       !! [H L ~> m2 or kg m-1].
  real,                                       intent(in)    :: dt      !< Time increment [T ~> s].
  type(unit_scale_type),                      intent(in)    :: US      !< A dimensional unit scaling type
  integer,                                    intent(in)    :: ish     !< Start of i index range.
  integer,                                    intent(in)    :: ieh     !< End of i index range.
  integer,                                    intent(in)    :: jsh     !< End of j index range.
  integer,                                    intent(in)    :: jeh     !< End of j index range.
  integer,                                    intent(in)    :: nz      !< End of k index range.
  logical, dimension(SZI_(G),SZJB_(G)),       intent(in)    :: do_I    !< Which i values to work on.
  logical,                                    intent(in)    :: vol_CFL !< If true, rescale the
         !! ratio of face areas to the cell areas when estimating the CFL number.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                              intent(in)    :: por_face_areaV !< fractional open area of V-faces [nondim]
  type(ocean_OBC_type),                   optional, pointer :: OBC !< Open boundaries control structure.
  ! Local variables
  real :: CFL ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h, i.e. [H ~> m or kg m-2].
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  integer :: i, j, k
  logical :: local_open_BC

  local_open_BC = .false.
  if (present(OBC)) then ; if (associated(OBC)) then
    local_open_BC = OBC%open_v_BCs_exist_globally
  endif ; endif

  do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh ; if (do_I(i,j)) then
    if (v(i,j,k) > 0.0) then
      if (vol_CFL) then ; CFL = (v(i,j,k) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
      else ; CFL = v(i,j,k) * dt * G%IdyT(i,j) ; endif
      curv_3 = (h_S(i,j,k) + h_N(i,j,k)) - 2.0*h(i,j,k)
      vh(i,j,k) = (G%dx_Cv(i,J)*por_face_areaV(i,J,k)) * v(i,j,k) * ( h_N(i,j,k) + CFL * &
          (0.5*(h_S(i,j,k) - h_N(i,j,k)) + curv_3*(CFL - 1.5)) )
      h_marg = h_N(i,j,k) + CFL * ((h_S(i,j,k) - h_N(i,j,k)) + &
                                  3.0*curv_3*(CFL - 1.0))
    elseif (v(i,j,k) < 0.0) then
      if (vol_CFL) then ; CFL = (-v(i,j,k) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
      else ; CFL = -v(i,j,k) * dt * G%IdyT(i,j+1) ; endif
      curv_3 = (h_S(i,j+1,k) + h_N(i,j+1,k)) - 2.0*h(i,j+1,k)
      vh(i,j,k) = (G%dx_Cv(i,J)*por_face_areaV(i,J,k)) * v(i,j,k) * ( h_S(i,j+1,k) + CFL * &
          (0.5*(h_N(i,j+1,k)-h_S(i,j+1,k)) + curv_3*(CFL - 1.5)) )
      h_marg = h_S(i,j+1,k) + CFL * ((h_N(i,j+1,k)-h_S(i,j+1,k)) + &
                                    3.0*curv_3*(CFL - 1.0))
    else
      vh(i,j,k) = 0.0
      h_marg = 0.5 * (h_S(i,j+1,k) + h_N(i,j,k))
    endif
    dvhdv(i,j,k) = (G%dx_Cv(i,J)*por_face_areaV(i,J,k)) * h_marg * visc_rem(i,j,k)
  endif ; enddo ; enddo; enddo

  if (local_open_BC) then
    do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh ; if (do_I(i,j)) then
      if (OBC%segnum_v(i,J) /= 0) then
        if (OBC%segment(abs(OBC%segnum_v(i,J)))%open) then
          if (OBC%segnum_v(i,J) > 0) then !  OBC_DIRECTION_N
            vh(i,j,k) = (G%dx_Cv(i,J)*por_face_areaV(i,J,k)) * v(i,j,k) * h(i,j,k)
            dvhdv(i,j,k) = (G%dx_Cv(i,J)*por_face_areaV(i,J,k)) * h(i,j,k) * visc_rem(i,j,k)
          else
            vh(i,j,k) = (G%dx_Cv(i,J)*por_face_areaV(i,J,k)) * v(i,j,k) * h(i,j+1,k)
            dvhdv(i,j,k) = (G%dx_Cv(i,J)*por_face_areaV(i,J,k)) * h(i,j+1,k) * visc_rem(i,j,k)
          endif
        endif
      endif
    endif ; enddo ; enddo ; enddo
  endif
end subroutine merid_flux_layer


!> Sets the effective interface thickness associated with the fluxes at each meridional velocity point,
!! optionally scaling back these thicknesses to account for viscosity and fractional open areas.
subroutine meridional_flux_thickness(v, h, h_S, h_N, h_v, dt, G, GV, US, LB, vol_CFL, &
                                     marginal, OBC, por_face_areaV, visc_rem_v)
  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)   :: v    !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h    !< Layer thickness used to calculate fluxes,
                                                                   !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_S  !< South edge thickness in the reconstruction,
                                                                   !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_N  !< North edge thickness in the reconstruction,
                                                                   !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: h_v !< Effective thickness at meridional faces,
                                                                   !! scaled down to account for the effects of
                                                                   !! viscosity and the fractional open area
                                                                   !! [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(cont_loop_bounds_type),               intent(in)    :: LB   !< Loop bounds structure.
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
  logical,                                   intent(in)    :: vol_CFL !< If true, rescale the ratio
                          !! of face areas to the cell areas when estimating the CFL number.
  logical,                                   intent(in)    :: marginal !< If true, report the marginal
                          !! face thicknesses; otherwise report transport-averaged thicknesses.
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                                     intent(in) :: por_face_areaV  !< fractional open area of V-faces [nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), optional, intent(in) :: visc_rem_v !< Both the fraction
                          !! of the momentum originally in a layer that remains after a time-step of
                          !! viscosity, and the fraction of a time-step's worth of a barotropic
                          !! acceleration that a layer experiences after viscosity is applied [nondim].
                          !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).

  ! Local variables
  real :: CFL ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h [H ~> m or kg m-2] .
  real :: h_avg  ! The average thickness of a flux [H ~> m or kg m-2].
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  logical :: local_open_BC
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  !$OMP parallel do default(shared) private(CFL,curv_3,h_marg,h_avg)
  do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
    if (v(i,J,k) > 0.0) then
      if (vol_CFL) then ; CFL = (v(i,J,k) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
      else ; CFL = v(i,J,k) * dt * G%IdyT(i,j) ; endif
      curv_3 = (h_S(i,j,k) + h_N(i,j,k)) - 2.0*h(i,j,k)
      h_avg = h_N(i,j,k) + CFL * (0.5*(h_S(i,j,k) - h_N(i,j,k)) + curv_3*(CFL - 1.5))
      h_marg = h_N(i,j,k) + CFL * ((h_S(i,j,k) - h_N(i,j,k)) + &
                                3.0*curv_3*(CFL - 1.0))
    elseif (v(i,J,k) < 0.0) then
      if (vol_CFL) then ; CFL = (-v(i,J,k)*dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
      else ; CFL = -v(i,J,k) * dt * G%IdyT(i,j+1) ; endif
      curv_3 = (h_S(i,j+1,k) + h_N(i,j+1,k)) - 2.0*h(i,j+1,k)
      h_avg = h_S(i,j+1,k) + CFL * (0.5*(h_N(i,j+1,k)-h_S(i,j+1,k)) + curv_3*(CFL - 1.5))
      h_marg = h_S(i,j+1,k) + CFL * ((h_N(i,j+1,k)-h_S(i,j+1,k)) + &
                                    3.0*curv_3*(CFL - 1.0))
    else
      h_avg = 0.5 * (h_S(i,j+1,k) + h_N(i,j,k))
      !   The choice to use the arithmetic mean here is somewhat arbitrarily, but
      ! it should be noted that h_S(i+1,j,k) and h_N(i,j,k) are usually the same.
      h_marg = 0.5 * (h_S(i,j+1,k) + h_N(i,j,k))
 !    h_marg = (2.0 * h_S(i,j+1,k) * h_N(i,j,k)) / &
 !             (h_S(i,j+1,k) + h_N(i,j,k) + GV%H_subroundoff)
    endif

    if (marginal) then ; h_v(i,J,k) = h_marg
    else ; h_v(i,J,k) = h_avg ; endif
  enddo ; enddo ; enddo

  if (present(visc_rem_v)) then
    ! Scale back the thickness to account for the effects of viscosity and the fractional open
    ! thickness to give an appropriate non-normalized weight for each layer in determining the
    ! barotropic acceleration.
    !$OMP parallel do default(shared)
    do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
      h_v(i,J,k) = h_v(i,J,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
    enddo ; enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
      h_v(i,J,k) = h_v(i,J,k) * por_face_areaV(i,J,k)
    enddo ; enddo ; enddo
  endif

  local_open_BC = .false.
  if (associated(OBC)) local_open_BC = OBC%open_v_BCs_exist_globally
  if (local_open_BC) then
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          if (present(visc_rem_v)) then ; do k=1,nz
            do i = OBC%segment(n)%HI%isd, OBC%segment(n)%HI%ied
              h_v(i,J,k) = h(i,j,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
            enddo
          enddo ; else ; do k=1,nz
            do i = OBC%segment(n)%HI%isd, OBC%segment(n)%HI%ied
              h_v(i,J,k) = h(i,j,k) * por_face_areaV(i,J,k)
            enddo
          enddo ; endif
        else
          if (present(visc_rem_v)) then ; do k=1,nz
            do i = OBC%segment(n)%HI%isd, OBC%segment(n)%HI%ied
              h_v(i,J,k) = h(i,j+1,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
            enddo
          enddo ; else ; do k=1,nz
            do i = OBC%segment(n)%HI%isd, OBC%segment(n)%HI%ied
              h_v(i,J,k) = h(i,j+1,k) * por_face_areaV(i,J,k)
            enddo
          enddo ; endif
        endif
      endif
    enddo
  endif

end subroutine meridional_flux_thickness


!> Returns the barotropic velocity adjustment that gives the desired barotropic (layer-summed) transport.
subroutine meridional_flux_adjust(v, h_in, h_S, h_N, vhbt, vh_tot_0, dvhdv_tot_0, &
                             dv, dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             ish, ieh, jsh, jeh, do_I_in, por_face_areaV, vh_3d, OBC)
  type(ocean_grid_type),             intent(in)  :: G    !< Ocean's grid structure.
  type(verticalGrid_type),           intent(in)  :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                     intent(in)  :: v    !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                     intent(in)  :: h_in !< Layer thickness used to calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),&
                                     intent(in)  :: h_S  !< South edge thickness in the reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                     intent(in)  :: h_N  !< North edge thickness in the reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                     intent(in)  :: visc_rem !< Both the fraction of the momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJB_(G)), intent(in)  :: vhbt !< The summed volume flux through meridional faces
                                                         !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)  :: dv_max_CFL !< Maximum acceptable value of dv [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)  :: dv_min_CFL !< Minimum acceptable value of dv [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)  :: vh_tot_0   !< The summed transport with 0 adjustment
                                                               !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)  :: dvhdv_tot_0 !< The partial derivative of dv_err with
                                                                !! dv at 0 adjustment [H L ~> m2 or kg m-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(out) :: dv         !< The barotropic velocity adjustment [L T-1 ~> m s-1].
  real,                              intent(in)  :: dt         !< Time increment [T ~> s].
  type(unit_scale_type),             intent(in)  :: US         !< A dimensional unit scaling type
  type(continuity_PPM_CS),           intent(in)  :: CS         !< This module's control structure.
  integer,                           intent(in)  :: ish        !< Start of i index range.
  integer,                           intent(in)  :: ieh        !< End of i index range.
  integer,                           intent(in)  :: jsh        !< Start of j index range.
  integer,                           intent(in)  :: jeh        !< End of j index range.
  logical, dimension(SZI_(G),SZJB_(G)), &
                                     intent(in)  :: do_I_in  !< A flag indicating which I values to work on.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                                     intent(in)  :: por_face_areaV !< fractional open area of V-faces [nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                         optional, intent(inout) :: vh_3d !< Volume flux through meridional
                             !! faces = v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(ocean_OBC_type), optional, pointer :: OBC !< Open boundaries control structure.
  ! Local variables
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    vh_aux, &  ! An auxiliary meridional volume flux [H L2 T-1 ~> m3 s-1 or kg s-1].
    dvhdv, &   ! Partial derivative of vh with v [H L ~> m2 or kg m-1].
    v_new      ! The velocity with the correction added [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    vh_err, &  ! Difference between vhbt and the summed vh [H L2 T-1 ~> m3 s-1 or kg s-1].
    vh_err_best, & ! The smallest value of vh_err found so far [H L2 T-1 ~> m3 s-1 or kg s-1].
    dvhdv_tot,&! Summed partial derivative of vh with u [H L ~> m2 or kg m-1].
    dv_min, &  ! Lower limit on dv correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
    dv_max     ! Upper limit on dv correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
  real :: dv_prev ! The previous value of dv [L T-1 ~> m s-1].
  real :: ddv     ! The change in dv from the previous iteration [L T-1 ~> m s-1].
  real :: tol_eta ! The tolerance for the current iteration [H ~> m or kg m-2].
  real :: tol_vel ! The tolerance for velocity in the current iteration [L T-1 ~> m s-1].
  integer :: i, j, k, nz, itt, max_itts = 20
  logical :: domore, do_I(SZI_(G),SZJB_(G))

  nz = GV%ke

  vh_aux(:,:,:) = 0.0 ; dvhdv(:,:,:) = 0.0
  
  if (present(vh_3d)) then ; do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh
    vh_aux(i,j,k) = vh_3d(i,J,k)
  enddo ; enddo ; enddo ; endif

  do j=jsh-1,jeh ; do i=ish,ieh
    dv(i,j) = 0.0 ; do_I(i,j) = do_I_in(i,j)
    dv_max(i,j) = dv_max_CFL(i,j) ; dv_min(i,j) = dv_min_CFL(i,j)
    vh_err(i,j) = vh_tot_0(i,j) - vhbt(i,j) ; dvhdv_tot(i,j) = dvhdv_tot_0(i,j)
    vh_err_best(i,j) = abs(vh_err(i,j))
  enddo ; enddo

  do itt=1,max_itts
    select case (itt)
      case (:1) ; tol_eta = 1e-6 * CS%tol_eta
      case (2)  ; tol_eta = 1e-4 * CS%tol_eta
      case (3)  ; tol_eta = 1e-2 * CS%tol_eta
      case default ; tol_eta = CS%tol_eta
    end select
    tol_vel = CS%tol_vel

    do j=jsh-1,jeh ; do i=ish,ieh
      if (vh_err(i,j) > 0.0) then ; dv_max(i,j) = dv(i,j)
      elseif (vh_err(i,j) < 0.0) then ; dv_min(i,j) = dv(i,j)
      else ; do_I(i,j) = .false. ; endif
    enddo ; enddo
    domore = .false.
    do j=jsh-1,jeh ; do i=ish,ieh ; if (do_I(i,j)) then
      if ((dt * min(G%IareaT(i,j),G%IareaT(i,j+1))*abs(vh_err(i,j)) > tol_eta) .or. &
          (CS%better_iter .and. ((abs(vh_err(i,j)) > tol_vel * dvhdv_tot(i,j)) .or. &
                                 (abs(vh_err(i,j)) > vh_err_best(i,j))) )) then
        !   Use Newton's method, provided it stays bounded.  Otherwise bisect
        ! the value with the appropriate bound.
        ddv = -vh_err(i,j) / dvhdv_tot(i,j)
        dv_prev = dv(i,j)
        dv(i,j) = dv(i,j) + ddv
        if (abs(ddv) < 1.0e-15*abs(dv(i,j))) then
          do_I(i,j) = .false. ! ddv is small enough to quit.
        elseif (ddv > 0.0) then
          if (dv(i,j) >= dv_max(i,j)) then
            dv(i,j) = 0.5*(dv_prev + dv_max(i,j))
            if (dv_max(i,j) - dv_prev < 1.0e-15*abs(dv(i,j))) do_I(i,j) = .false.
          endif
        else ! dvv(i) < 0.0
          if (dv(i,j) <= dv_min(i,j)) then
            dv(i,j) = 0.5*(dv_prev + dv_min(i,j))
            if (dv_prev - dv_min(i,j) < 1.0e-15*abs(dv(i,j))) do_I(i,j) = .false.
          endif
        endif
        if (do_I(i,j)) domore = .true.
      else
        do_I(i,j) = .false.
      endif
    endif ; enddo ; enddo
    if (.not.domore) exit

    if ((itt < max_itts) .or. present(vh_3d)) then ; do k=1,nz ; do j=jsh-1,jeh
      do i=ish,ieh ; v_new(i,j,k) = v(i,J,k) + dv(i,j) * visc_rem(i,j,k) ; enddo ; enddo ; enddo
      call merid_flux_layer(v_new, h_in, h_S, h_N, &
                            vh_aux, dvhdv, visc_rem, &
                            dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaV, OBC)
    endif

    if (itt < max_itts) then
      do j=jsh-1,jeh ; do i=ish,ieh
        vh_err(i,j) = -vhbt(i,j) ; dvhdv_tot(i,j) = 0.0
      enddo ; enddo
      do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh
        vh_err(i,j) = vh_err(i,j) + vh_aux(i,j,k)
        dvhdv_tot(i,j) = dvhdv_tot(i,j) + dvhdv(i,j,k)
      enddo ; enddo ; enddo
      do j=jsh-1,jeh ; do i=ish,ieh
        vh_err_best(i,j) = min(vh_err_best(i,j), abs(vh_err(i,j)))
      enddo ; enddo
    endif
  enddo ! itt-loop
  
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.

  if (present(vh_3d)) then ; do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh
    vh_3d(i,J,k) = vh_aux(i,j,k)
  enddo ; enddo ; enddo ; endif

end subroutine meridional_flux_adjust


!> Sets of a structure that describes the meridional barotropic volume or mass fluxes as a
!! function of barotropic flow to agree closely with the sum of the layer's transports.
subroutine set_merid_BT_cont(v, h_in, h_S, h_N, BT_cont, vh_tot_0, dvhdv_tot_0, &
                             dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             visc_rem_max, ish, ieh, jsh, jeh, do_I, por_face_areaV)
  type(ocean_grid_type),                      intent(in)    :: G    !< Ocean's grid structure.
  type(verticalGrid_type),                    intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v    !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h_in !< Layer thickness used to calculate fluxes,
                                                                    !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h_S  !< South edge thickness in the reconstruction,
                                                                    !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h_N  !< North edge thickness in the reconstruction,
                                                                    !! [H ~> m or kg m-2].
  type(BT_cont_type),                         intent(inout) :: BT_cont !< A structure with elements
                       !! that describe the effective open face areas as a function of barotropic flow.
  real, dimension(SZI_(G),SZJB_(G)),          intent(in)    :: vh_tot_0 !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)),          intent(in)    :: dvhdv_tot_0 !< The partial derivative
                       !! of du_err with dv at 0 adjustment [H L ~> m2 or kg m-1].
  real, dimension(SZI_(G),SZJB_(G)),          intent(in)    :: dv_max_CFL !< Maximum acceptable value
                                                                          !!  of dv [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)),          intent(in)    :: dv_min_CFL !< Minimum acceptable value
                                                                          !!  of dv [L T-1 ~> m s-1].
  real,                                       intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),                      intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS),                    intent(in)    :: CS   !< This module's control structure.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: visc_rem !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step
                       !! of viscosity, and the fraction of a time-step's worth of a barotropic
                       !! acceleration that a layer experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJB_(G)),          intent(in)    :: visc_rem_max !< Maximum allowable visc_rem [nondim]
  integer,                                    intent(in)    :: ish  !< Start of i index range.
  integer,                                    intent(in)    :: ieh  !< End of i index range.
  integer,                                    intent(in)    :: jsh  !< Start of j index range.
  integer,                                    intent(in)    :: jeh  !< End of j index range.
  logical, dimension(SZI_(G),SZJB_(G)),       intent(in)    :: do_I !< A logical flag indicating
                                                                    !! which I values to work on.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)),  intent(in)    :: por_face_areaV !< fractional open area of V-faces 
                                                                              !! [nondim]
  ! Local variables
  real, dimension(SZI_(G),SZJB_(G)) :: &
    dv0, &             ! The barotropic velocity increment that gives 0 transport [L T-1 ~> m s-1].
    dvL, dvR, &        ! The barotropic velocity increments that give the southerly
                       ! (dvL) and northerly (dvR) test velocities [L T-1 ~> m s-1].
    dv_CFL, &          ! The velocity increment that corresponds to CFL_min [L T-1 ~> m s-1].
    zeros, &           ! An array of full of 0 transports [H L2 T-1 ~> m3 s-1 or kg s-1]
    FAmt_L, FAmt_R, &  ! The summed effective marginal face areas for the 3
    FAmt_0, &          ! test velocities [H L ~> m2 or kg m-1].
    vhtot_L, &         ! The summed transport with the southerly (vhtot_L) and
    vhtot_R            ! and northerly (vhtot_R) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &    
    v_L, v_R, &        ! The southerly (v_L), northerly (v_R), and zero-barotropic
    v_0, &             ! transport (v_0) layer test velocities [L T-1 ~> m s-1].
    dvhdv_L, &         ! The effective layer marginal face areas with the southerly
    dvhdv_R, &         ! (_L), northerly (_R), and zero-barotropic (_0) test
    dvhdv_0, &         ! velocities [H L ~> m2 or kg m-1].
    vh_L, vh_R, &      ! The layer transports with the southerly (_L), northerly (_R)
    vh_0               ! and zero-barotropic (_0) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: FA_0         ! The effective face area with 0 barotropic transport [H L ~> m2 or kg m-1].
  real :: FA_avg       ! The average effective face area [H L ~> m2 or kg m-1], nominally given by
                       ! the realized transport divided by the barotropic velocity.
  real :: visc_rem_lim ! The larger of visc_rem and min_visc_rem [nondim]  This
                       ! limiting is necessary to keep the inverse of visc_rem
                       ! from leading to large CFL numbers.
  real :: min_visc_rem ! The smallest permitted value for visc_rem that is used
                       ! in finding the barotropic velocity that changes the
                       ! flow direction [nondim].  This is necessary to keep the inverse
                       ! of visc_rem from leading to large CFL numbers.
  real :: CFL_min      ! A minimal increment in the CFL to try to ensure that the
                       ! flow is truly upwind [nondim]
  real :: Idt          ! The inverse of the time step [T-1 ~> s-1].
  logical :: domore
  integer :: i, j, k, nz

  nz = GV%ke ; Idt = 1.0 / dt
  min_visc_rem = 0.1 ; CFL_min = 1e-6

 ! Diagnose the zero-transport correction, dv0.
  do j=jsh-1,jeh ; do i=ish,ieh ; zeros(i,j) = 0.0 ; enddo ; enddo
  call meridional_flux_adjust(v, h_in, h_S, h_N, zeros, vh_tot_0, dvhdv_tot_0, dv0, &
                         dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                         ish, ieh, jsh, jeh, do_I, por_face_areaV)

  ! Determine the southerly- and northerly- fluxes. Choose a sufficiently
  ! negative velocity correction for the northerly-flux, and a sufficiently
  ! positive correction for the southerly-flux.
  domore = .false.
  do j=jsh-1,jeh ; do i=ish,ieh ; if (do_I(i,j)) then
    domore = .true.
    dv_CFL(i,j) = (CFL_min * Idt) * G%dyCv(i,J)
    dvR(i,j) = min(0.0,dv0(i,j) - dv_CFL(i,j))
    dvL(i,j) = max(0.0,dv0(i,j) + dv_CFL(i,j))
    FAmt_L(i,j) = 0.0 ; FAmt_R(i,j) = 0.0 ; FAmt_0(i,j) = 0.0
    vhtot_L(i,j) = 0.0 ; vhtot_R(i,j) = 0.0
  endif ; enddo ; enddo

  if (.not.domore) then
    do j=jsh-1,jeh ; do i=ish,ieh
      BT_cont%FA_v_S0(i,J) = 0.0 ; BT_cont%FA_v_SS(i,J) = 0.0
      BT_cont%vBT_SS(i,J) = 0.0
      BT_cont%FA_v_N0(i,J) = 0.0 ; BT_cont%FA_v_NN(i,J) = 0.0
      BT_cont%vBT_NN(i,J) = 0.0
    enddo ; enddo
    return
  endif

  do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh ; if (do_I(i,j)) then
    visc_rem_lim = max(visc_rem(i,j,k), min_visc_rem*visc_rem_max(i,j))
    if (visc_rem_lim > 0.0) then ! This is almost always true for ocean points.
      if (v(i,J,k) + dvR(i,j)*visc_rem_lim > -dv_CFL(i,j)*visc_rem(i,j,k)) &
        dvR(i,j) = -(v(i,J,k) + dv_CFL(i,j)*visc_rem(i,j,k)) / visc_rem_lim
      if (v(i,J,k) + dvL(i,j)*visc_rem_lim < dv_CFL(i,j)*visc_rem(i,j,k)) &
        dvL(i,j) = -(v(i,J,k) - dv_CFL(i,j)*visc_rem(i,j,k)) / visc_rem_lim
    endif
  endif ; enddo ; enddo ; enddo

  do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh ; if (do_I(i,j)) then
    v_L(i,j,k) = v(I,j,k) + dvL(i,j) * visc_rem(i,j,k)
    v_R(i,j,k) = v(I,j,k) + dvR(i,j) * visc_rem(i,j,k)
    v_0(i,j,k) = v(I,j,k) + dv0(i,j) * visc_rem(i,j,k)
  endif ; enddo ; enddo ; enddo
  call merid_flux_layer(v_0, h_in, h_S, h_N, vh_0, dvhdv_0, &
                        visc_rem, dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaV)
  call merid_flux_layer(v_L, h_in, h_S, h_N, vh_L, dvhdv_L, &
                        visc_rem, dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaV)
  call merid_flux_layer(v_R, h_in, h_S, h_N, vh_R, dvhdv_R, &
                        visc_rem, dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaV)
  do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh ; if (do_I(i,j)) then
    FAmt_0(i,j) = FAmt_0(i,j) + dvhdv_0(i,j,k)
    FAmt_L(i,j) = FAmt_L(i,j) + dvhdv_L(i,j,k)
    FAmt_R(i,j) = FAmt_R(i,j) + dvhdv_R(i,j,k)
    vhtot_L(i,j) = vhtot_L(i,j) + vh_L(i,j,k)
    vhtot_R(i,j) = vhtot_R(i,j) + vh_R(i,j,k)
  endif ; enddo ; enddo ; enddo

  do j = jsh-1,jeh ; do i=ish,ieh ; if (do_I(i,j)) then
    FA_0 = FAmt_0(i,j) ; FA_avg = FAmt_0(i,j)
    if ((dvL(i,j) - dv0(i,j)) /= 0.0) &
      FA_avg = vhtot_L(i,j) / (dvL(i,j) - dv0(i,j))
    if (FA_avg > max(FA_0, FAmt_L(i,j))) then ; FA_avg = max(FA_0, FAmt_L(i,j))
    elseif (FA_avg < min(FA_0, FAmt_L(i,j))) then ; FA_0 = FA_avg ; endif
    BT_cont%FA_v_S0(i,J) = FA_0 ; BT_cont%FA_v_SS(i,J) = FAmt_L(i,j)
    if (abs(FA_0-FAmt_L(i,j)) <= 1e-12*FA_0) then ; BT_cont%vBT_SS(i,J) = 0.0 ; else
      BT_cont%vBT_SS(i,J) = (1.5 * (dvL(i,j) - dv0(i,j))) * &
                   ((FAmt_L(i,j) - FA_avg) / (FAmt_L(i,j) - FA_0))
    endif

    FA_0 = FAmt_0(i,j) ; FA_avg = FAmt_0(i,j)
    if ((dvR(i,j) - dv0(i,j)) /= 0.0) &
      FA_avg = vhtot_R(i,j) / (dvR(i,j) - dv0(i,j))
    if (FA_avg > max(FA_0, FAmt_R(i,j))) then ; FA_avg = max(FA_0, FAmt_R(i,j))
    elseif (FA_avg < min(FA_0, FAmt_R(i,j))) then ; FA_0 = FA_avg ; endif
    BT_cont%FA_v_N0(i,J) = FA_0 ; BT_cont%FA_v_NN(i,J) = FAmt_R(i,j)
    if (abs(FAmt_R(i,j) - FA_0) <= 1e-12*FA_0) then ; BT_cont%vBT_NN(i,J) = 0.0 ; else
      BT_cont%vBT_NN(i,J) = (1.5 * (dvR(i,j) - dv0(i,j))) * &
                   ((FAmt_R(i,j) - FA_avg) / (FAmt_R(i,j) - FA_0))
    endif
  else
    BT_cont%FA_v_S0(i,J) = 0.0 ; BT_cont%FA_v_SS(i,J) = 0.0
    BT_cont%FA_v_N0(i,J) = 0.0 ; BT_cont%FA_v_NN(i,J) = 0.0
    BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
  endif ; enddo ; enddo

end subroutine set_merid_BT_cont

!> Calculates left/right edge values for PPM reconstruction.
subroutine PPM_reconstruction_x(h_in, h_W, h_E, G, GV, LB, h_min, monotonic, simple_2nd, OBC)
  type(ocean_grid_type),             intent(in)  :: G    !< Ocean's grid structure.
  type(verticalGrid_type),           intent(in)  :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_in !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(out) :: h_W  !< West edge thickness in the reconstruction,
                                                         !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(out) :: h_E  !< East edge thickness in the reconstruction,
                                                         !! [H ~> m or kg m-2].
  type(cont_loop_bounds_type),       intent(in)  :: LB   !< Active loop bounds structure.
  real,                              intent(in)  :: h_min !< The minimum thickness
                    !! that can be obtained by a concave parabolic fit [H ~> m or kg m-2]
  logical,                           intent(in)  :: monotonic !< If true, use the
                    !! Colella & Woodward monotonic limiter.
                    !! Otherwise use a simple positive-definite limiter.
  logical,                           intent(in)  :: simple_2nd !< If true, use the
                    !! arithmetic mean thicknesses as the default edge values
                    !! for a simple 2nd order scheme.
  type(ocean_OBC_type),              pointer     :: OBC !< Open boundaries control structure.

  ! Local variables with useful mnemonic names.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: slp ! The slopes per grid point [H ~> m or kg m-2]
  real, parameter :: oneSixth = 1./6.  ! [nondim]
  real :: h_ip1, h_im1 ! Neighboring thicknesses or sensibly extrapolated values [H ~> m or kg m-2]
  real :: dMx, dMn     ! The difference between the local thickness and the maximum (dMx) or
                       ! minimum (dMn) of the surrounding values [H ~> m or kg m-2]
  character(len=256) :: mesg
  integer :: i, j, k, isl, iel, jsl, jel, nz, n, stencil
  logical :: local_open_BC
  type(OBC_segment_type), pointer :: segment => NULL()

  local_open_BC = .false.
  if (associated(OBC)) then
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif

  isl = LB%ish-1 ; iel = LB%ieh+1 ; jsl = LB%jsh ; jel = LB%jeh ; nz = GV%ke

  ! This is the stencil of the reconstruction, not the scheme overall.
  stencil = 2 ; if (simple_2nd) stencil = 1

  if ((isl-stencil < G%isd) .or. (iel+stencil > G%ied)) then
    write(mesg,'("In MOM_continuity_PPM, PPM_reconstruction_x called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               stencil + max(G%isd-isl,iel-G%ied)
    call MOM_error(FATAL,mesg)
  endif
  if ((jsl < G%jsd) .or. (jel > G%jed)) then
    write(mesg,'("In MOM_continuity_PPM, PPM_reconstruction_x called with a ", &
               & "y-halo that needs to be increased by ",i2,".")') &
               max(G%jsd-jsl,jel-G%jed)
    call MOM_error(FATAL,mesg)
  endif

  !$omp target enter data &
  !$omp   map(to: G, G%mask2dT(isl-2:iel+2, jsl:jel), h_in(isl-2:iel+2, :, :)) &
  !$omp   map(alloc: h_W, h_E, slp) ! whole arrays as unclear range of loops in OBC regions

  if (simple_2nd) then
    ! untested
    !$omp target teams distribute parallel do collapse(3) &
    !$omp   private(h_im1, h_ip1) &
    !$omp   map(to: G, G%mask2dT(isl-1:iel+1, jsl:jel), h_in(isl-1:iel+1, :, :)) &
    !$omp   map(from: h_W(isl:iel, :, :), h_E(isl:iel, :, :))
    do k =1,nz ; do j=jsl,jel ; do i=isl,iel
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j,k) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j,k)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j,k) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j,k)
      h_W(i,j,k) = 0.5*( h_im1 + h_in(i,j,k) )
      h_E(i,j,k) = 0.5*( h_ip1 + h_in(i,j,k) )
    enddo ; enddo ; enddo
  else
    !$omp target teams distribute parallel do collapse(3) &
    !$omp   private(dMx, dMn) &
    !$omp   map(to: G, G%mask2dT(isl-2:iel+2, jsl:jel), h_in(isl-2:iel+2, :, :)) &
    !$omp   map(from: slp(isl:iel, :, :))
    do k=1,nz ; do j=jsl,jel ; do i=isl-1,iel+1
      if ((G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j)) == 0.0) then
        slp(i,j,k) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j,k) = 0.5 * (h_in(i+1,j,k) - h_in(i-1,j,k))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i+1,j,k), h_in(i-1,j,k), h_in(i,j,k)) - h_in(i,j,k)
        dMn = h_in(i,j,k) - min(h_in(i+1,j,k), h_in(i-1,j,k), h_in(i,j,k))
        slp(i,j,k) = sign(1.,slp(i,j,k)) * min(abs(slp(i,j,k)), 2. * min(dMx, dMn))
                ! * (G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j))
      endif
    enddo ; enddo ; enddo

    if (local_open_BC) then
      ! untested
      do n=1, OBC%number_of_segments
        segment => OBC%segment(n)
        if (.not. segment%on_pe) cycle
        if (segment%is_E_or_W) then
          I=segment%HI%IsdB
          !$omp target teams distribute parallel do collapse(2) &
          !$omp   map(to: segment) &
          !$omp   map(tofrom: slp(i:i+1, :, :))
          do k=1, nz ; do j=segment%HI%jsd,segment%HI%jed
            slp(i+1,j,k) = 0.0
            slp(i,j,k) = 0.0
          enddo ; enddo
        endif
      enddo
    endif

    !$omp target teams distribute parallel do collapse(3) &
    !$omp   private(h_im1, h_ip1) &
    !$omp   map(to: G, G%mask2dT(isl-1:iel+1, jsl:jel), h_in(isl-1:iel+1, :, :), &
    !$omp       slp(isl-1:iel+1, :, :)) &
    !$omp   map(from: h_W(isl:iel, :, :), h_E(isl:iel, :, :))
    do k=1,nz ; do j=jsl,jel ; do i=isl,iel
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
    ! h_im1 = h_in(i-1,j,k) ; if (G%mask2dT(i-1,j) < 0.5) h_im1 = h_in(i,j)
    ! h_ip1 = h_in(i+1,j,k) ; if (G%mask2dT(i+1,j) < 0.5) h_ip1 = h_in(i,j)
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j,k) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j,k)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j,k) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j,k)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_W(i,j,k) = 0.5*( h_im1 + h_in(i,j,k) ) + oneSixth*( slp(i-1,j,k) - slp(i,j,k) )
      h_E(i,j,k) = 0.5*( h_ip1 + h_in(i,j,k) ) + oneSixth*( slp(i,j,k) - slp(i+1,j,k) )
    enddo ; enddo ; enddo
  endif

  if (local_open_BC) then
    ! untested
    do n=1, OBC%number_of_segments
      segment => OBC%segment(n)
      if (.not. segment%on_pe) cycle
      if (segment%direction == OBC_DIRECTION_E) then
        I=segment%HI%IsdB
        !$omp target teams distribute parallel do collapse(2) &
        !$omp   map(to: segment, h_in(i, :, :)) &
        !$omp   map(tofrom: h_W(i:i+1, :, :), h_E(i:i+1, :, :))
        do k=1,nz ; do j=segment%HI%jsd,segment%HI%jed
          h_W(i+1,j,k) = h_in(i,j,k)
          h_E(i+1,j,k) = h_in(i,j,k)
          h_W(i,j,k) = h_in(i,j,k)
          h_E(i,j,k) = h_in(i,j,k)
        enddo ; enddo
      elseif (segment%direction == OBC_DIRECTION_W) then
        I=segment%HI%IsdB
        !$omp target teams distribute parallel do collapse(2) &
        !$omp   map(to: segment, h_in(i, :, :)) &
        !$omp   map(tofrom: h_W(i:i+1, :, :), h_E(i:i+1, :, :))
        do k=1,nz ; do j=segment%HI%jsd,segment%HI%jed
          h_W(i,j,k) = h_in(i+1,j,k)
          h_E(i,j,k) = h_in(i+1,j,k)
          h_W(i+1,j,k) = h_in(i+1,j,k)
          h_E(i+1,j,k) = h_in(i+1,j,k)
        enddo ; enddo
      endif
    enddo
  endif

  if (monotonic) then
    ! untested
    call PPM_limit_CW84(h_in, h_W, h_E, G, GV, isl, iel, jsl, jel, nz)
  else
    call PPM_limit_pos(h_in, h_W, h_E, h_min, G, GV, isl, iel, jsl, jel, nz)
  endif

  !$omp target exit data &
  !$omp   map(from: h_W, h_E) & ! whole arrays as unclear range of loops in OBC regions
  !$omp   map(release: G, G%mask2dT(isl-2:iel+2, jsl:jel), h_in(isl-2:iel+2, :, :), slp)

  return
end subroutine PPM_reconstruction_x

!> Calculates left/right edge values for PPM reconstruction.
subroutine PPM_reconstruction_y(h_in, h_S, h_N, G, GV, LB, h_min, monotonic, simple_2nd, OBC)
  type(ocean_grid_type),             intent(in)  :: G    !< Ocean's grid structure.
  type(verticalGrid_type),           intent(in)  :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_in !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(out) :: h_S  !< South edge thickness in the reconstruction,
                                                         !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(out) :: h_N  !< North edge thickness in the reconstruction,
                                                         !! [H ~> m or kg m-2].
  type(cont_loop_bounds_type),       intent(in)  :: LB   !< Active loop bounds structure.
  real,                              intent(in)  :: h_min !< The minimum thickness
                    !! that can be obtained by a concave parabolic fit [H ~> m or kg m-2]
  logical,                           intent(in)  :: monotonic !< If true, use the
                    !! Colella & Woodward monotonic limiter.
                    !! Otherwise use a simple positive-definite limiter.
  logical,                           intent(in)  :: simple_2nd !< If true, use the
                    !! arithmetic mean thicknesses as the default edge values
                    !! for a simple 2nd order scheme.
  type(ocean_OBC_type),              pointer     :: OBC !< Open boundaries control structure.

  ! Local variables with useful mnemonic names.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: slp ! The slopes per grid point [H ~> m or kg m-2]
  real, parameter :: oneSixth = 1./6.      ! [nondim]
  real :: h_jp1, h_jm1 ! Neighboring thicknesses or sensibly extrapolated values [H ~> m or kg m-2]
  real :: dMx, dMn     ! The difference between the local thickness and the maximum (dMx) or
                       ! minimum (dMn) of the surrounding values [H ~> m or kg m-2]
  character(len=256) :: mesg
  integer :: i, j, k, isl, iel, jsl, jel, nz, n, stencil
  logical :: local_open_BC
  type(OBC_segment_type), pointer :: segment => NULL()

  local_open_BC = .false.
  if (associated(OBC)) then
    local_open_BC = OBC%open_v_BCs_exist_globally
  endif

  isl = LB%ish ; iel = LB%ieh ; jsl = LB%jsh-1 ; jel = LB%jeh+1 ; nz = G%ke

  ! This is the stencil of the reconstruction, not the scheme overall.
  stencil = 2 ; if (simple_2nd) stencil = 1

  if ((isl < G%isd) .or. (iel > G%ied)) then
    write(mesg,'("In MOM_continuity_PPM, PPM_reconstruction_y called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               max(G%isd-isl,iel-G%ied)
    call MOM_error(FATAL,mesg)
  endif
  if ((jsl-stencil < G%jsd) .or. (jel+stencil > G%jed)) then
    write(mesg,'("In MOM_continuity_PPM, PPM_reconstruction_y called with a ", &
                 & "y-halo that needs to be increased by ",i2,".")') &
                 stencil + max(G%jsd-jsl,jel-G%jed)
    call MOM_error(FATAL,mesg)
  endif

  if (simple_2nd) then
    ! untested
    do k=1,nz ; do j=jsl,jel ; do i=isl,iel
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1,k) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j,k)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1,k) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j,k)
      h_S(i,j,k) = 0.5*( h_jm1 + h_in(i,j,k) )
      h_N(i,j,k) = 0.5*( h_jp1 + h_in(i,j,k) )
    enddo ; enddo ; enddo
  else
    do k=1,nz ; do j=jsl-1,jel+1 ; do i=isl,iel
      if ((G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1)) == 0.0) then
        slp(i,j,k) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j,k) = 0.5 * (h_in(i,j+1,k) - h_in(i,j-1,k))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i,j+1,k), h_in(i,j-1,k), h_in(i,j,k)) - h_in(i,j,k)
        dMn = h_in(i,j,k) - min(h_in(i,j+1,k), h_in(i,j-1,k), h_in(i,j,k))
        slp(i,j,k) = sign(1.,slp(i,j,k)) * min(abs(slp(i,j,k)), 2. * min(dMx, dMn))
                ! * (G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1))
      endif
    enddo ; enddo ; enddo

    if (local_open_BC) then
      ! untested
      do n=1, OBC%number_of_segments
        segment => OBC%segment(n)
        if (.not. segment%on_pe) cycle
        if (segment%is_N_or_S) then
          J=segment%HI%JsdB
          do k=1,nz ; do i=segment%HI%isd,segment%HI%ied
            slp(i,j+1,k) = 0.0
            slp(i,j,k) = 0.0
          enddo ; enddo
        endif
      enddo
    endif

    do k=1,nz ; do j=jsl,jel ; do i=isl,iel
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1,k) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j,k)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1,k) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j,k)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_S(i,j,k) = 0.5*( h_jm1 + h_in(i,j,k) ) + oneSixth*( slp(i,j-1,k) - slp(i,j,k) )
      h_N(i,j,k) = 0.5*( h_jp1 + h_in(i,j,k) ) + oneSixth*( slp(i,j,k) - slp(i,j+1,k) )
    enddo ; enddo ; enddo
  endif

  if (local_open_BC) then
    ! untested
    do n=1, OBC%number_of_segments
      segment => OBC%segment(n)
      if (.not. segment%on_pe) cycle
      if (segment%direction == OBC_DIRECTION_N) then
        J=segment%HI%JsdB
        do k=1,nz ; do i=segment%HI%isd,segment%HI%ied
          h_S(i,j+1,k) = h_in(i,j,k)
          h_N(i,j+1,k) = h_in(i,j,k)
          h_S(i,j,k) = h_in(i,j,k)
          h_N(i,j,k) = h_in(i,j,k)
        enddo ; enddo
      elseif (segment%direction == OBC_DIRECTION_S) then
        J=segment%HI%JsdB
        do k=1,nz ; do i=segment%HI%isd,segment%HI%ied
          h_S(i,j,k) = h_in(i,j+1,k)
          h_N(i,j,k) = h_in(i,j+1,k)
          h_S(i,j+1,k) = h_in(i,j+1,k)
          h_N(i,j+1,k) = h_in(i,j+1,k)
        enddo ; enddo
      endif
    enddo
  endif

  if (monotonic) then
    ! untested
    call PPM_limit_CW84(h_in, h_S, h_N, G, GV, isl, iel, jsl, jel, nz)
  else
    call PPM_limit_pos(h_in, h_S, h_N, h_min, G, GV, isl, iel, jsl, jel, nz)
  endif

  return
end subroutine PPM_reconstruction_y

!> This subroutine limits the left/right edge values of the PPM reconstruction
!! to give a reconstruction that is positive-definite.  Here this is
!! reinterpreted as giving a constant thickness if the mean thickness is less
!! than h_min, with a minimum of h_min otherwise.
subroutine PPM_limit_pos(h_in, h_L, h_R, h_min, G, GV, iis, iie, jis, jie, nz)
  type(ocean_grid_type),             intent(in)  :: G    !< Ocean's grid structure.
  type(verticalGrid_type),           intent(in)  :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_in !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h_L !< Left thickness in the reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h_R !< Right thickness in the reconstruction [H ~> m or kg m-2].
  real,                              intent(in)  :: h_min !< The minimum thickness
                    !! that can be obtained by a concave parabolic fit [H ~> m or kg m-2]
  integer,                           intent(in)  :: iis      !< Start of i index range.
  integer,                           intent(in)  :: iie      !< End of i index range.
  integer,                           intent(in)  :: jis      !< Start of j index range.
  integer,                           intent(in)  :: jie      !< End of j index range.
  integer,                           intent(in)  :: nz       !< End of k index range.

! Local variables
  real    :: curv  ! The grid-normalized curvature of the three thicknesses  [H ~> m or kg m-2]
  real    :: dh    ! The difference between the edge thicknesses             [H ~> m or kg m-2]
  real    :: scale ! A scaling factor to reduce the curvature of the fit               [nondim]
  integer :: i,j,k

  !$omp target teams distribute parallel do collapse(3) &
  !$omp   private(curv, dh, scale) &
  !$omp   map(to: h_in(iis:iie, :, :)) &
  !$omp   map(tofrom: h_L(iis:iie, :, :), h_R(iis:iie, :, :))
  do k=1,nz ; do j=jis,jie ; do i=iis,iie
    ! This limiter prevents undershooting minima within the domain with
    ! values less than h_min.
    curv = 3.0*((h_L(i,j,k) + h_R(i,j,k)) - 2.0*h_in(i,j,k))
    if (curv > 0.0) then ! Only minima are limited.
      dh = h_R(i,j,k) - h_L(i,j,k)
      if (abs(dh) < curv) then ! The parabola's minimum is within the cell.
        if (h_in(i,j,k) <= h_min) then
          h_L(i,j,k) = h_in(i,j,k) ; h_R(i,j,k) = h_in(i,j,k)
        elseif (12.0*curv*(h_in(i,j,k) - h_min) < (curv**2 + 3.0*dh**2)) then
          ! The minimum value is h_in - (curv^2 + 3*dh^2)/(12*curv), and must
          ! be limited in this case.  0 < scale < 1.
          scale = 12.0*curv*(h_in(i,j,k) - h_min) / (curv**2 + 3.0*dh**2)
          h_L(i,j,k) = h_in(i,j,k) + scale*(h_L(i,j,k) - h_in(i,j,k))
          h_R(i,j,k) = h_in(i,j,k) + scale*(h_R(i,j,k) - h_in(i,j,k))
        endif
      endif
    endif
  enddo ; enddo ; enddo

end subroutine PPM_limit_pos

!> This subroutine limits the left/right edge values of the PPM reconstruction
!! according to the monotonic prescription of Colella and Woodward, 1984.
subroutine PPM_limit_CW84(h_in, h_L, h_R, G, GV, iis, iie, jis, jie, nz)
  type(ocean_grid_type),             intent(in)  :: G     !< Ocean's grid structure.
  type(verticalGrid_type),           intent(in)  :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h_in  !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h_L !< Left thickness in the reconstruction,
                                                          !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h_R !< Right thickness in the reconstruction,
                                                          !! [H ~> m or kg m-2].
  integer,                           intent(in)  :: iis   !< Start of i index range.
  integer,                           intent(in)  :: iie   !< End of i index range.
  integer,                           intent(in)  :: jis   !< Start of j index range.
  integer,                           intent(in)  :: jie   !< End of j index range.
  integer,                           intent(in)  :: nz    !< End of k index range.

  ! Local variables
  real    :: h_i      ! A copy of the cell-average layer thickness                [H ~> m or kg m-2]
  real    :: RLdiff   ! The difference between the input edge values              [H ~> m or kg m-2]
  real    :: RLdiff2  ! The squared difference between the input edge values   [H2 ~> m2 or kg2 m-4]
  real    :: RLmean   ! The average of the input edge thicknesses                 [H ~> m or kg m-2]
  real    :: FunFac   ! A curious product of the thickness slope and curvature [H2 ~> m2 or kg2 m-4]
  integer :: i, j, k

  ! untested
  !!$omp target teams distribute parallel do collapse(3) &
  !!$omp   private(h_i, RLdiff, RLdiff2, RLmean, FunFac) &
  !!$omp   map(to: h_in(iis:iie, :, :)) &
  !!$omp   map(tofrom: h_L(iis:iie, :, :), h_R(iis:iie, :, :))
  do k=1,nz ; do j=jis,jie ; do i=iis,iie
    ! This limiter monotonizes the parabola following
    ! Colella and Woodward, 1984, Eq. 1.10
    h_i = h_in(i,j,k)
    if ( ( h_R(i,j,k) - h_i ) * ( h_i - h_L(i,j,k) ) <= 0. ) then
      h_L(i,j,k) = h_i ; h_R(i,j,k) = h_i
    else
      RLdiff = h_R(i,j,k) - h_L(i,j,k)            ! Difference of edge values
      RLmean = 0.5 * ( h_R(i,j,k) + h_L(i,j,k) )  ! Mean of edge values
      FunFac = 6. * RLdiff * ( h_i - RLmean ) ! Some funny factor
      RLdiff2 = RLdiff * RLdiff               ! Square of difference
      if ( FunFac >  RLdiff2 ) h_L(i,j,k) = 3. * h_i - 2. * h_R(i,j,k)
      if ( FunFac < -RLdiff2 ) h_R(i,j,k) = 3. * h_i - 2. * h_L(i,j,k)
    endif
  enddo ; enddo ; enddo

end subroutine PPM_limit_CW84

!> Return the maximum ratio of a/b or maxrat.
function ratio_max(a, b, maxrat) result(ratio)
  real, intent(in) :: a       !< Numerator, in arbitrary units [A]
  real, intent(in) :: b       !< Denominator, in arbitrary units [B]
  real, intent(in) :: maxrat  !< Maximum value of ratio [A B-1]
  real :: ratio               !< Return value [A B-1]

  if (abs(a) > abs(maxrat*b)) then
    ratio = maxrat
  else
    ratio = a / b
  endif
end function ratio_max

!> Initializes continuity_ppm_cs
subroutine continuity_PPM_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< Vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure indicating
                  !! the open file to parse for model parameter values.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to
                  !! regulate diagnostic output.
  type(continuity_PPM_CS), intent(inout) :: CS   !< Module's control structure.

  !> This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_continuity_PPM" ! This module's name.

  CS%initialized = .true.

! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MONOTONIC_CONTINUITY", CS%monotonic, &
                 "If true, CONTINUITY_PPM uses the Colella and Woodward "//&
                 "monotonic limiter.  The default (false) is to use a "//&
                 "simple positive definite limiter.", default=.false.)
  call get_param(param_file, mdl, "SIMPLE_2ND_PPM_CONTINUITY", CS%simple_2nd, &
                 "If true, CONTINUITY_PPM uses a simple 2nd order "//&
                 "(arithmetic mean) interpolation of the edge values. "//&
                 "This may give better PV conservation properties. While "//&
                 "it formally reduces the accuracy of the continuity "//&
                 "solver itself in the strongly advective limit, it does "//&
                 "not reduce the overall order of accuracy of the dynamic "//&
                 "core.", default=.false.)
  call get_param(param_file, mdl, "UPWIND_1ST_CONTINUITY", CS%upwind_1st, &
                 "If true, CONTINUITY_PPM becomes a 1st-order upwind "//&
                 "continuity solver.  This scheme is highly diffusive "//&
                 "but may be useful for debugging or in single-column "//&
                 "mode where its minimal stencil is useful.", default=.false.)
  call get_param(param_file, mdl, "ETA_TOLERANCE", CS%tol_eta, &
                 "The tolerance for the differences between the "//&
                 "barotropic and baroclinic estimates of the sea surface "//&
                 "height due to the fluxes through each face.  The total "//&
                 "tolerance for SSH is 4 times this value.  The default "//&
                 "is 0.5*NK*ANGSTROM, and this should not be set less "//&
                 "than about 10^-15*MAXIMUM_DEPTH.", units="m", scale=GV%m_to_H, &
                 default=0.5*GV%ke*GV%Angstrom_m)

  call get_param(param_file, mdl, "VELOCITY_TOLERANCE", CS%tol_vel, &
                 "The tolerance for barotropic velocity discrepancies "//&
                 "between the barotropic solution and  the sum of the "//&
                 "layer thicknesses.", units="m s-1", default=3.0e8, scale=US%m_s_to_L_T)
                 ! The speed of light is the default.

  call get_param(param_file, mdl, "CONT_PPM_AGGRESS_ADJUST", CS%aggress_adjust,&
                 "If true, allow the adjusted velocities to have a "//&
                 "relative CFL change up to 0.5.", default=.false.)
  CS%vol_CFL = CS%aggress_adjust
  call get_param(param_file, mdl, "CONT_PPM_VOLUME_BASED_CFL", CS%vol_CFL, &
                 "If true, use the ratio of the open face lengths to the "//&
                 "tracer cell areas when estimating CFL numbers.  The "//&
                 "default is set by CONT_PPM_AGGRESS_ADJUST.", &
                 default=CS%aggress_adjust, do_not_read=CS%aggress_adjust)
  call get_param(param_file, mdl, "CONTINUITY_CFL_LIMIT", CS%CFL_limit_adjust, &
                 "The maximum CFL of the adjusted velocities.", units="nondim", &
                 default=0.5)
  call get_param(param_file, mdl, "CONT_PPM_BETTER_ITER", CS%better_iter, &
                 "If true, stop corrective iterations using a velocity "//&
                 "based criterion and only stop if the iteration is "//&
                 "better than all predecessors.", default=.true.)
  call get_param(param_file, mdl, "CONT_PPM_USE_VISC_REM_MAX", CS%use_visc_rem_max, &
                 "If true, use more appropriate limiting bounds for "//&
                 "corrections in strongly viscous columns.", default=.true.)
  call get_param(param_file, mdl, "CONT_PPM_MARGINAL_FACE_AREAS", CS%marginal_faces, &
                 "If true, use the marginal face areas from the continuity "//&
                 "solver for use as the weights in the barotropic solver. "//&
                 "Otherwise use the transport averaged areas.", default=.true.)
  CS%diag => diag

  id_clock_reconstruct = cpu_clock_id('(Ocean continuity reconstruction)', grain=CLOCK_ROUTINE)
  id_clock_update = cpu_clock_id('(Ocean continuity update)', grain=CLOCK_ROUTINE)
  id_clock_correct = cpu_clock_id('(Ocean continuity correction)', grain=CLOCK_ROUTINE)

end subroutine continuity_PPM_init

!> continuity_PPM_stencil returns the continuity solver stencil size
function continuity_PPM_stencil(CS) result(stencil)
  type(continuity_PPM_CS), intent(in) :: CS   !< Module's control structure.
  integer ::  stencil !< The continuity solver stencil size with the current settings.

  stencil = 3 ; if (CS%simple_2nd) stencil = 2 ; if (CS%upwind_1st) stencil = 1

end function continuity_PPM_stencil

!> Set up a structure that stores the sizes of the i- and j-loops to to work on in the continuity solver.
function set_continuity_loop_bounds(G, CS, i_stencil, j_stencil) result(LB)
  type(ocean_grid_type),   intent(in) :: G   !< The ocean's grid structure.
  type(continuity_PPM_CS), intent(in) :: CS  !< Module's control structure.
  logical,       optional, intent(in) :: i_stencil !< If present and true, extend the i-loop bounds
                                             !! by the stencil width of the continuity scheme.
  logical,       optional, intent(in) :: j_stencil !< If present and true, extend the j-loop bounds
                                             !! by the stencil width of the continuity scheme.
  type(cont_loop_bounds_type)         :: LB  !< A type storing the array sizes to work on in the continuity routines.

  ! Local variables
  logical :: add_i_stencil, add_j_stencil ! Local variables set based on i_stencil and j_stensil
  integer :: stencil    ! The continuity solver stencil size with the current continuity scheme.

  add_i_stencil = .false. ; if (present(i_stencil)) add_i_stencil = i_stencil
  add_j_stencil = .false. ; if (present(j_stencil)) add_j_stencil = j_stencil

  stencil = continuity_PPM_stencil(CS)

  if (add_i_stencil) then
    LB%ish = G%isc-stencil ; LB%ieh = G%iec+stencil
  else
    LB%ish = G%isc ; LB%ieh = G%iec
  endif

  if (add_j_stencil) then
    LB%jsh = G%jsc-stencil ; LB%jeh = G%jec+stencil
  else
    LB%jsh = G%jsc ; LB%jeh = G%jec
  endif

end function set_continuity_loop_bounds

!> \namespace mom_continuity_ppm
!!
!! This module contains the subroutines that advect layer
!! thickness.  The scheme here uses a Piecewise-Parabolic method with
!! a positive definite limiter.

end module MOM_continuity_PPM
