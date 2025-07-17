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
  !$omp target update to(visc_rem_u, u, visc_rem_v, v)
  ! problems with hin
  !$omp target enter data &
  !$omp   map(to: G, G%dy_Cu, G%IareaT, G%IdxT, G%areaT, G%dxT, G%mask2dCu, G%dxCu, G%IareaT, &
  !$omp       G%mask2dT, G%dx_Cv, G%dyCv, G%dyT, G%IdyT, G%mask2dCv, GV, u, v, h, CS, pbv, &
  !$omp       pbv%por_face_areaU, pbv%por_face_areaV, uhbt, vhbt, visc_rem_u, visc_rem_v, BT_cont) &
  !$omp   map(alloc: h_W, h_E, h_S, h_N, uh, vh, u_cor, v_cor, du_cor, dv_cor, BT_cont%FA_u_E0, &
  !$omp       BT_cont%FA_u_W0, BT_cont%FA_v_N0, BT_cont%FA_v_S0, BT_cont%FA_u_EE, BT_cont%FA_u_WW, &
  !$omp       BT_cont%FA_v_NN, BT_cont%FA_v_SS, BT_cont%uBT_EE, BT_cont%uBT_WW, BT_cont%vBT_NN, &
  !$omp       BT_cont%vBT_SS, BT_cont%h_u, BT_cont%h_V, LB)

  if (x_first) then
    !  First advect zonally, with loop bounds that accomodate the subsequent meridional advection.
    LB = set_continuity_loop_bounds(G, CS, i_stencil=.false., j_stencil=.true.)
    !$omp target update to(LB)
    call zonal_edge_thickness(hin, h_W, h_E, G, GV, US, CS, OBC, LB)
    call zonal_mass_flux(u, hin, h_W, h_E, uh, dt, G, GV, US, CS, OBC, pbv%por_face_areaU, &
                         LB, uhbt, visc_rem_u, u_cor, BT_cont, du_cor)
    call continuity_zonal_convergence(h, uh, dt, G, GV, LB, hin)

    ! update host h from continuity_zonal_convergence

    !  Now advect meridionally, using the updated thicknesses to determine the fluxes.
    LB = set_continuity_loop_bounds(G, CS, i_stencil=.false., j_stencil=.false.)
    !$omp target update to(LB)
    call meridional_edge_thickness(h, h_S, h_N, G, GV, US, CS, OBC, LB)
    call meridional_mass_flux(v, h, h_S, h_N, vh, dt, G, GV, US, CS, OBC, pbv%por_face_areaV, &
                              LB, vhbt, visc_rem_v, v_cor, BT_cont, dv_cor)
    call continuity_merdional_convergence(h, vh, dt, G, GV, LB, hmin=h_min)

  else  ! .not. x_first
    !  First advect meridionally, with loop bounds that accomodate the subsequent zonal advection.
    LB = set_continuity_loop_bounds(G, CS, i_stencil=.true., j_stencil=.false.)
    !$omp target update to(LB)
    call meridional_edge_thickness(hin, h_S, h_N, G, GV, US, CS, OBC, LB)
    call meridional_mass_flux(v, hin, h_S, h_N, vh, dt, G, GV, US, CS, OBC, pbv%por_face_areaV, &
                              LB, vhbt, visc_rem_v, v_cor, BT_cont, dv_cor)
    call continuity_merdional_convergence(h, vh, dt, G, GV, LB, hin)

    !  Now advect zonally, using the updated thicknesses to determine the fluxes.
    LB = set_continuity_loop_bounds(G, CS, i_stencil=.false., j_stencil=.false.)
    !$omp target update to(LB)
    call zonal_edge_thickness(h, h_W, h_E, G, GV, US, CS, OBC, LB)
    call zonal_mass_flux(u, h, h_W, h_E, uh, dt, G, GV, US, CS, OBC, pbv%por_face_areaU, &
                         LB, uhbt, visc_rem_u, u_cor, BT_cont, du_cor)
    call continuity_zonal_convergence(h, uh, dt, G, GV, LB, hmin=h_min)
  endif

  !$omp target update from(uh, u_cor, du_cor, vh, v_cor, dv_cor)
  !$omp target update from(BT_cont%FA_u_W0, BT_cont%FA_u_WW, BT_cont%FA_u_E0, BT_cont%FA_u_EE, &
  !$omp       BT_cont%uBT_WW, BT_cont%uBT_EE, BT_cont%h_u)
  !$omp target update from(BT_cont%FA_v_S0, BT_cont%FA_v_SS, BT_cont%FA_v_N0, BT_cont%FA_v_NN, &
  !$omp       BT_cont%vBT_SS, BT_cont%vBT_NN, BT_cont%h_v)
  !$omp target update from(h)

  !$omp target exit data &
  !$omp   map(from: uh, h, vh, u_cor, v_cor, du_cor, dv_cor, BT_cont%FA_u_E0, BT_cont%FA_u_W0, &
  !$omp       BT_cont%FA_v_N0, BT_cont%FA_v_S0, BT_cont%FA_u_EE, BT_cont%FA_u_WW, BT_cont%FA_v_NN, &
  !$omp       BT_cont%FA_v_SS, BT_cont%uBT_EE, BT_cont%uBT_WW, BT_cont%vBT_NN, BT_cont%vBT_SS, &
  !$omp       BT_cont%h_u, BT_cont%h_V) &
  !$omp   map(release: G, G%dy_Cu, G%IareaT, G%IdxT, G%areaT, G%dxT, G%mask2dCu, G%dxCu, G%IareaT, &
  !$omp       G%mask2dT, G%dyCv, G%dyT, G%IdyT, G%mask2dCv, GV, u, v, h_W, h_E, h_S, h_N, CS, pbv, &
  !$omp       pbv%por_face_areaU, pbv%por_face_areaV, uhbt, vhbt, visc_rem_u, visc_rem_v, LB)

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
    do concurrent (k=1:nz, j=jsh:jeh, i=ish:ieh)
      h(i,j,k) = max( hin(i,j,k) - dt * G%IareaT(i,j) * (uh(I,j,k) - uh(I-1,j,k)), h_min )
    enddo
  else
    ! untested
    do concurrent (k=1:nz, j=jsh:jeh, i=ish:ieh)
      h(i,j,k) = max( h(i,j,k) - dt * G%IareaT(i,j) * (uh(I,j,k) - uh(I-1,j,k)), h_min )
    enddo
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
    ! untested
    do concurrent (k=1:GV%ke, j=LB%jsh:LB%jeh, i=LB%ish:LB%ieh)
      h(i,j,k) = max( hin(i,j,k) - dt * G%IareaT(i,j) * (vh(i,J,k) - vh(i,J-1,k)), h_min )
    enddo
  else
    do concurrent (k=1:GV%ke, j=LB%jsh:LB%jeh, i=LB%ish:LB%ieh)
      h(i,j,k) = max( h(i,j,k) - dt * G%IareaT(i,j) * (vh(i,J,k) - vh(i,J-1,k)), h_min )
    enddo
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
    do concurrent (k=1:nz, j=jsh:jeh, i=ish-1:ieh+1)
      h_W(i,j,k) = h_in(i,j,k) ; h_E(i,j,k) = h_in(i,j,k)
    enddo
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
    ! untested
    do concurrent (k=1:nz, j=jsh-1:jeh+1, i=ish:ieh)
      h_S(i,j,k) = h_in(i,j,k) ; h_N(i,j,k) = h_in(i,j,k)
    enddo
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
  real, dimension(SZIB_(G),SZJ_(G), SZK_(GV)) :: &
    visc_rem_u_tmp      ! A 2-D copy of visc_rem_u or an array of 1's [nondim].
  real :: FAuI  ! A sum of zonal face areas [H L ~> m2 or kg m-1].
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

  call cpu_clock_begin(id_clock_correct)

  use_visc_rem = present(visc_rem_u)

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))

  local_specified_BC = .false. ; local_Flather_OBC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_u_BCs_exist_globally
    local_Flather_OBC = OBC%Flather_u_BCs_exist_globally
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif ; endif

  if (present(LB_in)) then
    LB = LB_in
  else
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
  endif
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  CFL_dt = CS%CFL_limit_adjust / dt
  I_dt = 1.0 / dt
  if (CS%aggress_adjust) CFL_dt = I_dt

  !$omp target enter data if (use_visc_rem) map(to: visc_rem_u)
  !$omp target enter data &
  !$omp   map(to: G, G%dy_Cu, G%IareaT, G%IdxT, G%areaT, G%dxT, G%mask2dCu, G%dxCu, u, h_in, &
  !$omp       h_W, h_E, CS, por_face_areaU, uhbt, BT_cont) &
  !$omp   map(alloc: visc_rem_u_tmp, uh, u_cor, BT_cont%FA_u_E0, BT_cont%FA_u_W0, &
  !$omp       BT_cont%FA_u_EE, BT_cont%FA_u_WW, BT_cont%uBT_EE, BT_cont%uBT_WW, BT_cont%h_u, &
  !$omp       du_cor, duhdu, du, du_min_CFL, du_max_CFL, duhdu_tot_0, uh_tot_0, visc_rem_max)

  do concurrent (j=jsh:jeh)

    if (present(du_cor)) then
      do concurrent (i=ish-1:ieh)
        du_cor(i,j) = 0.0
      enddo
    endif

    if (.not.use_visc_rem) then
      do concurrent (k=1:nz, i=ish-1:ieh)
        visc_rem_u_tmp(i,j,k) = 1.0
      enddo
    else
      ! this is expensive
      do concurrent (k=1:nz, i=ish-1:ieh)
        visc_rem_u_tmp(i,j,k) = visc_rem_u(i,j,k)
      enddo
    end if
  
    ! Set uh and duhdu.
    !!DIR$ FORCEINLINE
    do concurrent (k=1:nz , I=ish-1:ieh)
      call zonal_flux_layere(u(I,j,k), h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), h_E(I+1,j,k), &
                            uh(I,j,k), duhdu(I,j,k), visc_rem_u_tmp(I,j,k), G%dy_Cu(I,j), G%IareaT(I,j), G%IareaT(I+1,j), G%IdxT(I,j), G%IdxT(i+1,j), &
                            dt, G, GV, US, CS%vol_CFL, por_face_areaU(I,j,k))
    enddo
    if (local_open_BC) then
      !!DIR$ FORCEINLINE
      do concurrent (k=1:nz, I=ish-1:ieh)
        call zonal_flux_layere_OBC(u(I,j,k), h_in, uh(I,j,k), duhdu(I,j,k), visc_rem_u_tmp(I,j,k), &
                                  G, GV, I, j, k, por_face_areaU(I,j,k), OBC)
      enddo
    endif
  
    ! untested!
    if (local_specified_BC) then
      do concurrent (k=1:nz, I=ish-1:ieh, OBC%segnum_u(I,j) /= 0)
        l_seg = abs(OBC%segnum_u(I,j))
        if (OBC%segment(l_seg)%specified) uh(I,j,k) = OBC%segment(l_seg)%normal_trans(I,j,k)
      enddo
    endif

    if (present(uhbt) .or. set_BT_cont) then
      if (use_visc_rem.and.CS%use_visc_rem_max) then
        ! poor performance for nvfortran + do concurrent if k is inside loop
        do concurrent (I=ish-1:ieh)
          visc_rem_max(I,j) = visc_rem_u_tmp(I,j,1)
        enddo
        do k=2,nz ; do concurrent (I=ish-1:ieh)
          visc_rem_max(I,j) = max(visc_rem_max(I,j), visc_rem_u_tmp(I,j,k))
        enddo ; enddo
      else
        do concurrent (i=ish-1:ieh)
          visc_rem_max(i, j) = 1.0
        enddo
      endif
      !   Set limits on du that will keep the CFL number between -1 and 1.
      ! This should be adequate to keep the root bracketed in all cases.
      do concurrent (I=ish-1:ieh)
        I_vrm = 0.0
        if (visc_rem_max(I,j) > 0.0) I_vrm = 1.0 / visc_rem_max(I,j)
        if (CS%vol_CFL) then
          dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
          dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
        else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif
        du_max_CFL(I,j) = 2.0* (CFL_dt * dx_W) * I_vrm
        du_min_CFL(I,j) = -2.0 * (CFL_dt * dx_E) * I_vrm
        uh_tot_0(I,j) = 0.0 ; duhdu_tot_0(I,j) = 0.0
      enddo
      ! poor performance for nvfortran + do concurrent if k is inside loop
      do k=1,nz ; do concurrent (I=ish-1:ieh)
        duhdu_tot_0(I,j) = duhdu_tot_0(I,j) + duhdu(I, j, k)
        uh_tot_0(I,j) = uh_tot_0(I,j) + uh(I,j,k)
      enddo ; enddo

      if (use_visc_rem) then
        if (CS%aggress_adjust) then
          ! untested!
          do k=1,nz ; do concurrent (I=ish-1:ieh)
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
          enddo ; enddo
        else
          do k=1,nz ; do concurrent (I=ish-1:ieh)
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            if (du_max_CFL(I,j) * visc_rem_u_tmp(I,j,k) > dx_W*CFL_dt - u(I,j,k)*G%mask2dCu(I,j)) &
              du_max_CFL(I,j) = (dx_W*CFL_dt - u(I,j,k)) / visc_rem_u_tmp(I,j,k)
            if (du_min_CFL(I,j) * visc_rem_u_tmp(I,j,k) < -dx_E*CFL_dt - u(I,j,k)*G%mask2dCu(I,j)) &
              du_min_CFL(I,j) = -(dx_E*CFL_dt + u(I,j,k)) / visc_rem_u_tmp(I,j,k)
          enddo ; enddo
        endif
      else
        ! untested!
        if (CS%aggress_adjust) then
          do k=1,nz ; do concurrent (I=ish-1:ieh)
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            du_max_CFL(I,j) = MIN(du_max_CFL(I,j), 0.499 * &
                        ((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k))) )
            du_min_CFL(I,j) = MAX(du_min_CFL(I,j), 0.499 * &
                        ((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k))) )
          enddo ; enddo
        else
          do k=1,nz ; do concurrent (I=ish-1:ieh)
            if (CS%vol_CFL) then
              dx_W = ratio_max(G%areaT(i,j), G%dy_Cu(I,j), 1000.0*G%dxT(i,j))
              dx_E = ratio_max(G%areaT(i+1,j), G%dy_Cu(I,j), 1000.0*G%dxT(i+1,j))
            else ; dx_W = G%dxT(i,j) ; dx_E = G%dxT(i+1,j) ; endif

            du_max_CFL(I,j) = MIN(du_max_CFL(I,j), dx_W*CFL_dt - u(I,j,k))
            du_min_CFL(I,j) = MAX(du_min_CFL(I,j), -(dx_E*CFL_dt + u(I,j,k)))
          enddo ; enddo
        endif
      endif
      do concurrent (I=ish-1:ieh)
        du_max_CFL(I,j) = max(du_max_CFL(I,j),0.0)
        du_min_CFL(I,j) = min(du_min_CFL(I,j),0.0)
      enddo
    endif
  enddo

  call present_uhbt_or_set_BT_cont(u, h_in, h_W, h_E, uh_tot_0, duhdu_tot_0, du, du_max_CFL, du_min_CFL, visc_rem_u_tmp, &
                                   visc_rem_max, por_face_areaU, uhbt, uh, u_cor, du_cor, BT_cont, dt, set_BT_cont, &
                                   local_specified_BC, local_Flather_OBC, local_open_BC, ish, ieh, jsh, jeh, nz, G, GV, US, CS, &
                                   OBC, LB)

  !$omp target exit data &
  !$omp   map(from: uh, u_cor, BT_cont%FA_u_E0, BT_cont%FA_u_W0, BT_cont%FA_u_EE, BT_cont%FA_u_WW, &
  !$omp       BT_cont%uBT_EE, BT_cont%uBT_WW, BT_cont%h_u, du_cor) &
  !$omp   map(release: visc_rem_u_tmp, G, G%dy_Cu, G%IareaT, G%IdxT, G%areaT, G%dxT, &
  !$Omp       G%mask2dCu, G%dxCu, u, h_in, h_W, h_E, CS, por_face_areaU, uhbt, BT_cont, duhdu, du, &
  !$omp       du_min_CFL, du_max_CFL, duhdu_tot_0, uh_tot_0, visc_rem_max)
  !$omp target exit data if (use_visc_rem) map(release: visc_rem_u)

  call cpu_clock_end(id_clock_correct)

end subroutine zonal_mass_flux

subroutine present_uhbt_or_set_BT_cont(u, h_in, h_W, h_E, uh_tot_0, duhdu_tot_0, du, du_max_CFL, du_min_CFL, visc_rem_u, visc_rem_max, por_face_areaU, uhbt, uh, u_cor, du_cor, BT_cont, dt, set_BT_cont, local_specified_BC, local_Flather_OBC, local_open_BC, ish, ieh, jsh, jeh, nz, G, GV, US, CS, OBC, LB)
  type(ocean_grid_type), intent(in):: G
  type(verticalGrid_type), intent(in)    :: GV   !< Ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u    !< Zonal velocity [L T-1 ~> m s-1].
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Layer thickness used to calculate fluxes [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_W !< Western edge thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_E !< Eastern edge thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: uh_tot_0 !< Summed transport with no barotropic correction [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: duhdu_tot_0 !< Summed partial derivative of uh with u [H L ~> m2 or kg m-1].
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(inout) :: du !< Corrective barotropic change in the velocity to give uhbt [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: du_max_CFL !< Upper limit on du correction to avoid CFL violations [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: du_min_CFL !< Lower limit on du correction to avoid CFL violations [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: uh   !< Volume flux through zonal faces = u*h*dy
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: visc_rem_u
                     !< The fraction of zonal momentum originally in a layer that remains after a
                     !! time-step of viscosity, and the fraction of a time-step's worth of a barotropic
                     !! acceleration that a layer experiences after viscosity is applied [nondim].
                     !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: visc_rem_max
  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)), &
                           intent(in)    :: por_face_areaU !< fractional open area of U-faces [nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: u_cor
                     !< The zonal velocities (u with a barotropic correction)
                     !! that give uhbt as the depth-integrated transport [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(out)   :: du_cor !< The zonal velocity increments from u that give uhbt
                                                 !! as the depth-integrated transports [L T-1 ~> m s-1].
  type(BT_cont_type), optional, pointer  :: BT_cont !< A structure with elements that describe the
                     !! effective open face areas as a function of barotropic flow.
  logical, intent(in):: set_BT_cont
  logical, intent(in):: local_specified_BC
  logical, intent(in):: local_Flather_OBC
  logical, intent(in):: local_open_BC
  real, dimension(SZIB_(G),SZJ_(G)), optional, intent(in):: uhbt
  integer, intent(in):: ish
  integer, intent(in):: ieh
  integer, intent(in):: jsh
  integer, intent(in):: jeh
  integer, intent(in):: nz
  real,                    intent(in)    :: dt   !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(continuity_PPM_CS), intent(in)    :: CS   !< This module's control structure.
  type(ocean_OBC_type), pointer:: OBC
  type(cont_loop_bounds_type), intent(in) :: LB
  ! Local variables
  logical, dimension(SZIB_(G), SZJ_(G)) :: do_I
  logical, dimension(SZIB_(G), SZJ_(G)) :: simple_OBC_pt  ! Indicates points in a row with specified transport OBCs
  logical:: any_simple_OBC
  integer:: l_seg, i, j, k, n
  real :: FAuI, FA_u

  if (present(uhbt) .or. set_BT_cont) then
    !$omp target enter data map(alloc: do_I)
    !$omp target enter data map(alloc: simple_OBC_pt) if(local_specified_BC .or. local_Flather_OBC)
    any_simple_OBC = .false.
    if (local_specified_BC .or. local_Flather_OBC) then
      do concurrent (j=jsh:jeh, I=ish-1:ieh)
        l_seg = abs(OBC%segnum_u(I,j))
        ! Avoid reconciling barotropic/baroclinic transports if transport is specified
        simple_OBC_pt(I,j) = .false.
        if (l_seg /= OBC_NONE) simple_OBC_pt(I,j) = OBC%segment(l_seg)%specified
        do_I(I, j) = .not.simple_OBC_pt(I,j)
        any_simple_OBC = any_simple_OBC .or. simple_OBC_pt(I,j)
      enddo
    else
      do concurrent (j=jsh:jeh, I=ish-1:ieh)
        do_I(I, j) = .true.
      enddo
    endif

    if (present(uhbt)) then
      ! Find du and uh.
      call zonal_flux_adjust(u, h_in, h_W, h_E, uh_tot_0, duhdu_tot_0, du, &
                            du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem_u, &
                            ish, ieh, jsh, jeh, do_I, por_face_areaU, uhbt, uh, OBC=OBC)
      
      do concurrent (j=jsh:jeh)
        if (present(u_cor)) then
          do concurrent (k=1:nz, I=ish-1:ieh)
            u_cor(I,j,k) = u(I,j,k) + du(I,j) * visc_rem_u(I,j,k)
          enddo
          if (any_simple_OBC) then 
            do concurrent (k=1:nz, I=ish-1:ieh, simple_OBC_pt(I,j))
              u_cor(I,j,k) = OBC%segment(abs(OBC%segnum_u(I,j)))%normal_vel(I,j,k)
            enddo
          endif
        endif ! u-corrected

        if (present(du_cor)) then
          do concurrent (I=ish-1:ieh) ; du_cor(I,j) = du(I,j) ; enddo
        endif
      enddo
    endif
    if (set_BT_cont) then
      ! Diagnose the zero-transport correction, du0.
      call zonal_flux_adjust(u, h_in, h_W, h_E, uh_tot_0, duhdu_tot_0, du, &
                            du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem_u, &
                            ish, ieh, jsh, jeh, do_I, por_face_areaU)
      call set_zonal_BT_cont(u, h_in, h_W, h_E, BT_cont, du, uh_tot_0, duhdu_tot_0,&
                              du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem_u, &
                              visc_rem_max, ish, ieh, jsh, jeh, do_I, por_face_areaU)
      if (any_simple_OBC) then
        ! untested
        do concurrent (j=jsh:jeh, I=ish-1:ieh)
          ! NOTE: simple_OBC_pt(I, j) should prevent access to segment OBC_NONE
          if (simple_OBC_pt(I,j)) then
            FAuI = GV%H_subroundoff*G%dy_Cu(I,j)
            do k=1,nz
              l_seg = abs(OBC%segnum_u(I,j))
              if ((abs(OBC%segment(l_seg)%normal_vel(I,j,k)) > 0.0) .and. (OBC%segment(l_seg)%specified)) &
                FAuI = FAuI + OBC%segment(l_seg)%normal_trans(I,j,k) / OBC%segment(l_seg)%normal_vel(I,j,k)
            enddo
            BT_cont%FA_u_W0(I,j) = FAuI ; BT_cont%FA_u_E0(I,j) = FAuI
            BT_cont%FA_u_WW(I,j) = FAuI ; BT_cont%FA_u_EE(I,j) = FAuI
            BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
          endif
        enddo
      endif
    endif
    !$omp target exit data map(release: simple_OBC_pt) if(local_specified_BC .or. local_Flather_OBC)
    !$omp target exit data map(release: do_I)
  endif

  if (local_open_BC .and. set_BT_cont) then
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          do concurrent (j = OBC%segment(n)%HI%Jsd:OBC%segment(n)%HI%Jed)
            FA_u = 0.0
            do k=1,nz ; FA_u = FA_u + h_in(i,j,k)*(G%dy_Cu(I,j)*por_face_areaU(I,j,k)) ; enddo
            BT_cont%FA_u_W0(I,j) = FA_u ; BT_cont%FA_u_E0(I,j) = FA_u
            BT_cont%FA_u_WW(I,j) = FA_u ; BT_cont%FA_u_EE(I,j) = FA_u
            BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
          enddo
        else
          do concurrent (j = OBC%segment(n)%HI%Jsd:OBC%segment(n)%HI%Jed)
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

end subroutine present_uhbt_or_set_BT_cont


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
elemental subroutine zonal_flux_layere(u, h, h_p1, h_W, h_W_p1, h_E, h_E_p1, uh, duhdu, visc_rem, G_dy_Cu, G_IareaT, G_IareaT_p1, &
                                      G_IdxT, G_IdxT_p1, dt, G, GV, US, vol_CFL, por_face_areaU)
  type(ocean_grid_type),    intent(in)    :: G        !< Ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV       !< Ocean's vertical grid structure.
  real,                     intent(in)    :: u        !< Zonal velocity [L T-1 ~> m s-1].
  real,                     intent(in)    :: visc_rem !< Both the fraction of the
                        !! momentum originally in a layer that remains after a time-step
                        !! of viscosity, and the fraction of a time-step's worth of a barotropic
                        !! acceleration that a layer experiences after viscosity is applied [nondim].
                        !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real,                      intent(in)    :: h, h_p1        !< Layer thickness [H ~> m or kg m-2].
  real,                      intent(in)    :: h_W, h_W_p1      !< West edge thickness [H ~> m or kg m-2].
  real,                      intent(in)    :: h_E, h_E_p1      !< East edge thickness [H ~> m or kg m-2].
  real,                     intent(out) :: uh       !< Zonal mass or volume
                                                      !! transport [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                     intent(out) :: duhdu    !< Partial derivative of uh
                                                      !! with u [H L ~> m2 or kg m-1].
  real,                     intent(in)    :: dt       !< Time increment [T ~> s]
  type(unit_scale_type),    intent(in)    :: US       !< A dimensional unit scaling type.
  logical,                  intent(in)    :: vol_CFL  !< If true, rescale the
  real,                     intent(in)    :: por_face_areaU !< fractional open area of U-faces 
                                                            !! [nondim].
  real,                     intent(in)    :: G_dy_Cu, G_IareaT, G_IareaT_p1, G_IdxT, G_IdxT_p1
          !! ratio of face areas to the cell areas when estimating the CFL number.
  ! Local variables
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length [H ~> m or kg m-2]
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  real :: tmp, dh, c1, c2, c3

  !DIR$ ATTRIBUTES FORCEINLINE :: zonal_flux_layere

  ! Set new values of uh and duhdu.
  tmp = G_dy_Cu * por_face_areaU ! precalculate things
  CFL = abs(u*dt) ! increases inlining likelihood
  c1 = 0.5 * (h_W_p1 + h_E)
  c2 = 0.0
  c3 = 0.0
  if (u > 0.0) then
    CFL = CFL * merge(G_dy_Cu * G_IareaT, G_IdxT, vol_CFL)
    c1 = h_E
    c2 = h_W
    c3 = h
  elseif (u < 0.0) then
    CFL = CFL * merge(G_dy_Cu * G_IareaT_p1, G_IdxT_p1, vol_CFL)
    c1 = h_W_p1
    c2 = h_E_p1
    c3 = h_p1
  endif
  curv_3 = (c2 + c1) - 2.0*c3
  dh = (c2 - c1)
  uh = tmp * u * &
        (c1 + CFL * (0.5*dh + curv_3*(CFL - 1.5)))
  h_marg = c1 + CFL * (dh + 3.0*curv_3*(CFL-1.0))
  duhdu = tmp * h_marg * visc_rem

end subroutine zonal_flux_layere

pure subroutine zonal_flux_layere_OBC(u, h, uh, duhdu, visc_rem, G, GV, i, j, k, por_face_areaU, OBC)
  type(ocean_grid_type),    intent(in)    :: G        !< Ocean's grid structure.
  type(verticalGrid_type),  intent(in)    :: GV       !< Ocean's vertical grid structure.
  real,                     intent(in)    :: u        !< Zonal velocity [L T-1 ~> m s-1].
  real,                     intent(in)    :: visc_rem !< Both the fraction of the
                        !! momentum originally in a layer that remains after a time-step
                        !! of viscosity, and the fraction of a time-step's worth of a barotropic
                        !! acceleration that a layer experiences after viscosity is applied [nondim].
                        !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h        !< Layer thickness [H ~> m or kg m-2].
  real,                     intent(inout) :: uh       !< Zonal mass or volume
                                                      !! transport [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                     intent(inout) :: duhdu    !< Partial derivative of uh
                                                      !! with u [H L ~> m2 or kg m-1].
  integer, value,           intent(in)    :: i      !< Start of i index range.
  integer, value,           intent(in)    :: j      !< Start of j index range.
  integer, value,           intent(in)    :: k       !< Edn of k index range.
  real,                     intent(in)    :: por_face_areaU !< fractional open area of U-faces 
                                                            !! [nondim].
          !! ratio of face areas to the cell areas when estimating the CFL number.
  type(ocean_OBC_type),     intent(in)    :: OBC !< Open boundaries control structure.
  integer :: l_seg

  ! untested
  if (OBC%segnum_u(I,j) /= OBC_NONE) then
    l_seg = OBC%segnum_u(I,j)
    if (OBC%segment(l_seg)%open) then
      if (OBC%segment(l_seg)%direction == OBC_DIRECTION_E) then
        uh = (G%dy_Cu(I,j) * por_face_areaU) * u * h(i, j, k)
        duhdu = (G%dy_Cu(I,j) * por_face_areaU) * h(i, j, k) * visc_rem
      else
        uh = (G%dy_Cu(I,j) * por_face_areaU) * u * h(i+1, j, k)
        duhdu = (G%dy_Cu(I,j)* por_face_areaU) * h(i+1, j, k) * visc_rem
      endif
    endif
  endif

end subroutine zonal_flux_layere_OBC

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

  !$omp target enter data &
  !$omp   map(to: do_I, u, G, G%dy_Cu, G%IareaT, G%IdxT, h_W, h_E, h, por_face_areaU, visc_rem, &
  !$omp       uh, duhdu)

  do concurrent (k=1:nz, j=jsh:jeh, I=ish-1:ieh, do_I(I, j))
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
  enddo

  if (local_open_BC) then
    ! untested
    do concurrent (k=1:nz, j=jsh:jeh, I=ish-1:ieh )
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
    enddo
  endif
  !$omp target exit data &
  !$omp   map(from: uh, duhdu) &
  !$omp   map(release: do_I, u, G, G%dy_Cu, G%IareaT, G%IdxT, h_W, h_E, h, por_face_areaU, visc_rem)
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
  logical :: local_open_BC
  integer :: i, j, k, ish, ieh, jsh, jeh, nz, n
  real :: dh
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = GV%ke

  !$omp target enter data &
  !$omp   map(to: u, G, G%dy_Cu, G%IareaT, G%IdxT, h, h_W, h_E, por_face_areaU) &
  !$omp   map(alloc: h_u)

  do concurrent (k=1:nz, j=jsh:jeh, I=ish-1:ieh)
    if (u(I,j,k) > 0.0) then
      if (vol_CFL) then ; CFL = (u(I,j,k) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
      else ; CFL = u(I,j,k) * dt * G%IdxT(i,j) ; endif
      curv_3 = (h_W(i,j,k) + h_E(i,j,k)) - 2.0*h(i,j,k)
      dh = h_W(i,j,k) - h_E(i,j,k)
      if (marginal) then
        h_u(I,j,k) = h_E(i,j,k) + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
      else
        h_u(I,j,k) = h_E(i,j,k) + CFL * (0.5*dh + curv_3*(CFL - 1.5))
      endif
    elseif (u(I,j,k) < 0.0) then
      if (vol_CFL) then ; CFL = (-u(I,j,k)*dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
      else ; CFL = -u(I,j,k) * dt * G%IdxT(i+1,j) ; endif
      curv_3 = (h_W(i+1,j,k) + h_E(i+1,j,k)) - 2.0*h(i+1,j,k)
      dh = h_E(i+1,j,k)-h_W(i+1,j,k)
      if (marginal) then
        h_u(I,j,k) = h_W(i+1,j,k) + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
      else
        h_u(I,j,k) = h_W(i+1,j,k) + CFL * (0.5*dh + curv_3*(CFL - 1.5))
      endif
    else
      !   The choice to use the arithmetic mean here is somewhat arbitrarily, but
      ! it should be noted that h_W(i+1,j,k) and h_E(i,j,k) are usually the same.
      h_u(I,j,k) = 0.5 * (h_W(i+1,j,k) + h_E(i,j,k))
 !    h_marg = (2.0 * h_W(i+1,j,k) * h_E(i,j,k)) / &
 !             (h_W(i+1,j,k) + h_E(i,j,k) + GV%H_subroundoff)
    endif
    
    if (present(visc_rem_u)) then
      ! Scale back the thickness to account for the effects of viscosity and the fractional open
      ! thickness to give an appropriate non-normalized weight for each layer in determining the
      ! barotropic acceleration.
      h_u(I,j,k) = h_u(I,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
    else
      h_u(I,j,k) = h_u(I,j,k) * por_face_areaU(I,j,k)
    endif
  enddo

  local_open_BC = .false.
  if (associated(OBC)) local_open_BC = OBC%open_u_BCs_exist_globally
  if (local_open_BC) then
    ! untested
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          if (present(visc_rem_u)) then
            do concurrent (k=1:nz, j = OBC%segment(n)%HI%jsd:OBC%segment(n)%HI%jed)
              h_u(I,j,k) = h(i,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
            enddo
          else
            do concurrent (k=1:nz, j = OBC%segment(n)%HI%jsd:OBC%segment(n)%HI%jed)
              h_u(I,j,k) = h(i,j,k) * por_face_areaU(I,j,k)
            enddo
          endif
        else
          if (present(visc_rem_u)) then 
            do concurrent (k=1:nz, j = OBC%segment(n)%HI%jsd:OBC%segment(n)%HI%jed)
              h_u(I,j,k) = h(i+1,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
            enddo
          else
            do concurrent (k=1:nz, j = OBC%segment(n)%HI%jsd:OBC%segment(n)%HI%jed)
              h_u(I,j,k) = h(i+1,j,k) * por_face_areaU(I,j,k)
            enddo
          endif
        endif
      endif
    enddo
  endif

  !$omp target exit data &
  !$omp   map(from: h_u) &
  !$omp   map(release: u, G, G%dy_Cu, G%IareaT, G%IdxT, h, h_W, h_E, por_face_areaU)

end subroutine zonal_flux_thickness

!> Returns the barotropic velocity adjustment that gives the
!! desired barotropic (layer-summed) transport.
subroutine zonal_flux_adjust(u, h_in, h_W, h_E, uh_tot_0, duhdu_tot_0, &
                             du, du_max_CFL, du_min_CFL, dt, G, GV, US, CS, visc_rem, &
                             ish, ieh, jsh, jeh, do_I_in, por_face_areaU, uhbt, uh_3d, OBC)

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
  real, dimension(SZIB_(G),SZJ_(G)), optional, intent(in)    :: uhbt !< The summed volume flux
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
  real :: &
    duhdu          ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
  real, dimension(SZIB_(G),SZK_(GV)) :: &
    uh_aux      ! An auxiliary zonal volume flux [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G)) :: &
    uh_err, &      ! Difference between uhbt and the summed uh [H L2 T-1 ~> m3 s-1 or kg s-1].
    uh_err_best, & ! The smallest value of uh_err found so far [H L2 T-1 ~> m3 s-1 or kg s-1].
    duhdu_tot,&    ! Summed partial derivative of uh with u [H L ~> m2 or kg m-1].
    du_min, &      ! Lower limit on du correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
    du_max         ! Upper limit on du correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
  real :: u_new ! The velocity with the correction added [L T-1 ~> m s-1].
  real :: du_prev  ! The previous value of du [L T-1 ~> m s-1].
  real :: ddu      ! The change in du from the previous iteration [L T-1 ~> m s-1].
  real :: tol_eta  ! The tolerance for the current iteration [H ~> m or kg m-2].
  real :: tol_vel  ! The tolerance for velocity in the current iteration [L T-1 ~> m s-1].
  integer :: i, j, k, nz, itt
  logical :: domore, do_I(SZIB_(G)), local_OBC, use_uhbt
  integer, parameter:: max_itts = 20

  local_OBC = .false.
  if (present(OBC)) then
    if (associated(OBC)) then
      local_OBC = OBC%open_u_BCs_exist_globally
    endif
  endif

  use_uhbt = present(uhbt)

  nz = GV%ke

  !$omp target enter data &
  !$omp   map(to: uh_3d, do_I_in, du_max_CFL, du_min_CFL, uh_tot_0, uhbt, &
  !$omp       duhdu_tot_0, G, G%IareaT, CS, u, visc_rem, h_W, h_E, h_in, G%dy_Cu, &
  !$omp       G%IdxT, por_face_areaU) &
  !$omp   map(alloc: uh_aux, du, do_I, du_max, du_min, duhdu_tot, uh_err, uh_err_best, u_new)

  ! NVIDIA do concurrent doesn't work with private arrays (private scalars OK)
  !$omp target loop private(i, k, itt, uh_aux, uh_err, uh_err_best, duhdu_tot, du_min, du_max, do_I)
  do j=jsh,jeh
    if (present(uh_3d)) then
      do concurrent (k=1:nz, I=ish-1:ieh)
        uh_aux(I,k) = uh_3d(I,j,k)
      enddo
    endif
    do concurrent (I=ish-1:ieh)
      du(I,j) = 0.0 ; do_I(I) = do_I_in(I,j)
      du_max(I) = du_max_CFL(I,j) ; du_min(I) = du_min_CFL(I,j)
      uh_err(I) = uh_tot_0(I,j)
      if (use_uhbt) uh_err(I) = uh_err(I) - uhbt(I,j)
      duhdu_tot(I) = duhdu_tot_0(I,j)
      uh_err_best(I) = abs(uh_err(I))
    enddo
    do itt=1,max_itts
      select case (itt)
        case (:1) ; tol_eta = 1e-6 * CS%tol_eta
        case (2)  ; tol_eta = 1e-4 * CS%tol_eta
        case (3)  ; tol_eta = 1e-2 * CS%tol_eta
        case default ; tol_eta = CS%tol_eta
      end select
      tol_vel = CS%tol_vel

      domore = .false.
      do concurrent (I=ish-1:ieh, do_I(I))
        if (uh_err(I) > 0.0) then ; du_max(I) = du(I,j)
        elseif (uh_err(I) < 0.0) then ; du_min(I) = du(I,j)
        else ; do_I(I) = .false. ; endif
        if ((dt * min(G%IareaT(i,j),G%IareaT(i+1,j))*abs(uh_err(I)) > tol_eta) .or. &
            (CS%better_iter .and. ((abs(uh_err(I)) > tol_vel * duhdu_tot(I)) .or. &
                                  (abs(uh_err(I)) > uh_err_best(I))) )) then
        !   Use Newton's method, provided it stays bounded.  Otherwise bisect
        ! the value with the appropriate bound.
          ddu = -uh_err(I) / duhdu_tot(I)
          du_prev = du(I,j)
          du(I,j) = du(I,j) + ddu
          if (abs(ddu) < 1.0e-15*abs(du(I,j))) then
            do_I(I) = .false. ! ddu is small enough to quit.
          elseif (ddu > 0.0) then
            if (du(I,j) >= du_max(I)) then
              du(I,j) = 0.5*(du_prev + du_max(I))
              if (du_max(I) - du_prev < 1.0e-15*abs(du(I,j))) do_I(I) = .false.
            endif
          else ! ddu < 0.0
            if (du(I,j) <= du_min(I)) then
              du(I,j) = 0.5*(du_prev + du_min(I))
              if (du_prev - du_min(I) < 1.0e-15*abs(du(I,j))) do_I(I) = .false.
            endif
          endif
          if (do_I(I)) domore = .true.
        else
          do_I(I) = .false.
        endif
      enddo

      ! Below conditional compilation is to control whether early exit happens when compiled with
      ! OpenMP - compiling with OpenMP prevents early exit. Without OpenMP, enables early exit.
      ! Early exit saves time on CPU, but causes other loops to be serialized on GPU.
      !$ if (.false.) then
      if (.not.domore) exit
      !$ endif

      if ((itt < max_itts) .or. present(uh_3d)) then
        do concurrent (I=ish-1:ieh)
          uh_err(I) = 0.0 ; duhdu_tot(I) = 0.0
          if (use_uhbt) uh_err(I) = -uhbt(I,j)
        enddo
        do k=1,nz ; do concurrent (I=ish-1:ieh, do_I(I))
          u_new = u(I,j,k) + du(I,j) * visc_rem(I,j,k)
          call zonal_flux_layere(u_new, h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), h_E(I+1,j,k), &
                                uh_aux(I,k), duhdu, visc_rem(I,j,k), G%dy_Cu(I,j), G%IareaT(I,j), G%IareaT(I+1,j), G%IdxT(I,j), &
                                G%IdxT(i+1,j), dt, G, GV, US, CS%vol_CFL, por_face_areaU(I,j,k))
          ! Below if statement looks expensive in profiling results, but I believe it's
          ! masking the expensive update of uh_err beneath 
          if (local_OBC) call zonal_flux_layere_OBC(u_new, h_in, uh_aux(I,k), duhdu, visc_rem(I,j,k), &
                                G, GV, I, j, k, por_face_areaU(I,j,k), OBC)
          uh_err(I) = uh_err(I) + uh_aux(I,k)
          duhdu_tot(I) = duhdu_tot(I) + duhdu
        enddo ; enddo
        do concurrent (I=ish-1:ieh)
          uh_err_best(I) = min(uh_err_best(I), abs(uh_err(I)))
        enddo
      endif

    enddo ! itt-loop
    if (present(uh_3d)) then
      do concurrent (k=1:nz, I=ish-1:ieh)
        uh_3d(I,j,k) = uh_aux(I,k)
      enddo
    endif
  enddo ! j-loop
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.

  !$omp target exit data &
  !$omp   map(from: uh_3d, du) &
  !$omp   map(release: uh_aux, do_I_in, du_max_CFL, du_min_CFL, uh_tot_0, uhbt, &
  !$omp       duhdu_tot_0, do_I, du_max, du_min, duhdu_tot, uh_err, uh_err_best, G, G%IareaT, CS, &
  !$omp       u, visc_rem, u_new, h_W, h_E, h_in, G%dy_Cu, G%IdxT, por_face_areaU)

end subroutine zonal_flux_adjust


!> Sets a structure that describes the zonal barotropic volume or mass fluxes as a
!! function of barotropic flow to agree closely with the sum of the layer's transports.
subroutine set_zonal_BT_cont(u, h_in, h_W, h_E, BT_cont, du0, uh_tot_0, duhdu_tot_0, &
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
                           intent(in)    :: du0  !< The barotropic velocity increment that gives 0 transport [L T-1 ~> m s-1].
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
  real, dimension(SZIB_(G)) :: &
    duL, duR, &       ! The barotropic velocity increments that give the westerly
    du_CFL, &         ! The velocity increment that corresponds to CFL_min [L T-1 ~> m s-1].
                      ! (duL) and easterly (duR) test velocities [L T-1 ~> m s-1].
    FAmt_L, FAmt_R, & ! The summed effective marginal face areas for the 3
    FAmt_0, &         ! test velocities [H L ~> m2 or kg m-1].
    uhtot_L, &        ! The summed transport with the westerly (uhtot_L) and
    uhtot_R           ! and easterly (uhtot_R) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: &
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
  !$omp   map(to: u, h_W, h_E, h_in, uh_tot_0, duhdu_tot_0, G, G%IareaT, G%dy_Cu, G%IdxT, G%dxCu, &
  !$omp       BT_cont, du_max_CFL, du_min_CFL, CS, visc_rem, visc_rem_max, do_I, por_face_areaU, du0) &
  !$omp   map(alloc: duL, duR, du_CFL, FAmt_L, FAmT_R, FAmt_0, uhtot_L, uhtot_R, &
  !$omp       BT_cont%FA_u_W0, &
  !$omp       BT_cont%FA_u_WW, BT_cont%FA_u_E0, BT_cont%FA_u_EE, BT_cont%uBT_WW, BT_cont%uBT_EE)

  !$omp target loop private(I, k, duL, duR, du_CFL, FAmt_L, FAmt_R, FAmt_0, uhtot_L, uhtot_R)
  do j=jsh,jeh
    ! Determine the westerly- and easterly- fluxes.  Choose a sufficiently
    ! negative velocity correction for the easterly-flux, and a sufficiently
    ! positive correction for the westerly-flux.
    do concurrent (I=ish-1:ieh)
      du_CFL(I) = (CFL_min * Idt) * G%dxCu(I,j)
      duR(I) = min(0.0,du0(I,j) - du_CFL(I))
      duL(I) = max(0.0,du0(I,j) + du_CFL(I))
      FAmt_L(I) = 0.0 ; FAmt_R(I) = 0.0 ; FAmt_0(I) = 0.0
      uhtot_L(I) = 0.0 ; uhtot_R(I) = 0.0
    enddo

  ! nvfortran do concurrent bad performance if k is inside
    do k=1,nz ; do concurrent (I=ish-1:ieh, do_I(I,j))
      visc_rem_lim = max(visc_rem(I,j,k), min_visc_rem*visc_rem_max(I,j))
      if (visc_rem_lim > 0.0) then ! This is almost always true for ocean points.
        if (u(I,j,k) + duR(I)*visc_rem_lim > -du_CFL(I)*visc_rem(I,j,k)) &
          duR(I) = -(u(I,j,k) + du_CFL(I)*visc_rem(I,j,k)) / visc_rem_lim
        if (u(I,j,k) + duL(I)*visc_rem_lim < du_CFL(I)*visc_rem(I,j,k)) &
          duL(I) = -(u(I,j,k) - du_CFL(I)*visc_rem(I,j,k)) / visc_rem_lim
      endif
    enddo ; enddo

    do k=1,nz ; do concurrent (I=ish-1:ieh, do_I(I,j))
      u_L = u(I,j,k) + duL(I) * visc_rem(I,j,k)
      u_R = u(I,j,k) + duR(I) * visc_rem(I,j,k)
      u_0 = u(I,j,k) + du0(I,j) * visc_rem(I,j,k)
      call zonal_flux_layere(u_0, h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), h_E(I+1,j,k), uh_0, duhdu_0, &
                            visc_rem(I,j,k), G%dy_Cu(I,j), G%IareaT(I,j), G%IareaT(I+1,j), G%IdxT(I,j), G%IdxT(i+1,j), dt, G, GV, &
                            US, CS%vol_CFL, por_face_areaU(I,j,k))
      call zonal_flux_layere(u_L, h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), h_E(I+1,j,k), uh_L, duhdu_L, &
                            visc_rem(I,j,k), G%dy_Cu(I,j), G%IareaT(I,j), G%IareaT(I+1,j), G%IdxT(I,j), G%IdxT(i+1,j), dt, G, GV, &
                            US, CS%vol_CFL, por_face_areaU(I,j,k))
      call zonal_flux_layere(u_R, h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), h_E(I+1,j,k), uh_R, duhdu_R, &
                            visc_rem(I,j,k), G%dy_Cu(I,j), G%IareaT(I,j), G%IareaT(I+1,j), G%IdxT(I,j), G%IdxT(i+1,j), dt, G, GV, &
                            US, CS%vol_CFL, por_face_areaU(I,j,k))
      FAmt_0(I) = FAmt_0(I) + duhdu_0
      FAmt_L(I) = FAmt_L(I) + duhdu_L
      FAmt_R(I) = FAmt_R(I) + duhdu_R
      uhtot_L(I) = uhtot_L(I) + uh_L
      uhtot_R(I) = uhtot_R(I) + uh_R
    enddo ; enddo
    do concurrent (I=ish-1:ieh) ; if (do_I(I,j)) then
      FA_0 = FAmt_0(I) ; FA_avg = FAmt_0(I)
      if ((duL(I) - du0(I,j)) /= 0.0) &
        FA_avg = uhtot_L(I) / (duL(I) - du0(I,j))
      if (FA_avg > max(FA_0, FAmt_L(I))) then ; FA_avg = max(FA_0, FAmt_L(I))
      elseif (FA_avg < min(FA_0, FAmt_L(I))) then ; FA_0 = FA_avg ; endif

      BT_cont%FA_u_W0(I,j) = FA_0 ; BT_cont%FA_u_WW(I,j) = FAmt_L(I)
      if (abs(FA_0-FAmt_L(I)) <= 1e-12*FA_0) then ; BT_cont%uBT_WW(I,j) = 0.0 ; else
        BT_cont%uBT_WW(I,j) = (1.5 * (duL(I) - du0(I,j))) * &
                              ((FAmt_L(I) - FA_avg) / (FAmt_L(I) - FA_0))
      endif

      FA_0 = FAmt_0(I) ; FA_avg = FAmt_0(I)
      if ((duR(I) - du0(I,j)) /= 0.0) &
        FA_avg = uhtot_R(I) / (duR(I) - du0(I,j))
      if (FA_avg > max(FA_0, FAmt_R(I))) then ; FA_avg = max(FA_0, FAmt_R(I))
      elseif (FA_avg < min(FA_0, FAmt_R(I))) then ; FA_0 = FA_avg ; endif

      BT_cont%FA_u_E0(I,j) = FA_0 ; BT_cont%FA_u_EE(I,j) = FAmt_R(I)
      if (abs(FAmt_R(I) - FA_0) <= 1e-12*FA_0) then ; BT_cont%uBT_EE(I,j) = 0.0 ; else
        BT_cont%uBT_EE(I,j) = (1.5 * (duR(I) - du0(I,j))) * &
                              ((FAmt_R(I) - FA_avg) / (FAmt_R(I) - FA_0))
      endif
    else
      BT_cont%FA_u_W0(I,j) = 0.0 ; BT_cont%FA_u_WW(I,j) = 0.0
      BT_cont%FA_u_E0(I,j) = 0.0 ; BT_cont%FA_u_EE(I,j) = 0.0
      BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
    endif ; enddo
  enddo

  !$omp target exit data &
  !$omp   map(from: BT_cont%FA_u_W0, BT_cont%FA_u_WW, BT_cont%FA_u_E0, BT_cont%FA_u_EE, &
  !$omp       BT_cont%uBT_WW, BT_cont%uBT_EE) &
  !$omp   map(release: u, h_W, h_E, h_in, uh_tot_0, duhdu_tot_0, du0, duL, duR, du_CFL, &
  !$omp       FAmt_L, FAmT_R, FAmt_0, uhtot_L, uhtot_R, &
  !$omp       G, G%IareaT, G%dy_Cu, G%IdxT, G%dxCu, BT_cont, &
  !$omp       du_max_CFL, du_min_CFL, CS, visc_rem, visc_rem_max, do_I, por_face_areaU)

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

  !$omp target enter data &
  !$omp   map(to: G, G%dx_Cv, G%IdyT, G%dyT, G%dyCv, G%mask2dCv, G%areaT, G%IareaT, v, h_in, h_S, &
  !$omp       h_N, CS, por_face_areaV, vhbt, visc_rem_v, BT_cont) &
  !$omp   map(alloc: vh, v_cor, BT_cont%FA_v_S0, BT_cont%FA_v_SS, BT_cont%vBT_SS, BT_cont%FA_v_N0, &
  !$omp       BT_cont%FA_v_NN, BT_cont%vBT_NN, BT_cont%h_v, dv_cor, dvhdv, dv, dv_min_CFL, &
  !$omp       dv_max_CFL, dvhdv_tot_0, vh_tot_0, visc_rem_max, do_I, FAvi, visc_rem_v_tmp, &
  !$omp       simple_OBC_pt)

  ! a better solution is needed
  if (.not.use_visc_rem) then
    do concurrent (k=1:nz, j=jsh-1:jeh, i=ish:ieh)
      visc_rem_v_tmp(i,j,k) = 1.0
    enddo
  else
    do concurrent (k=1:nz, j=jsh-1:jeh, i=ish:ieh)
      visc_rem_v_tmp(i,j,k) = visc_rem_v(i,j,k)
    enddo
  endif
  do concurrent (j=jsh-1:jeh, i=ish:ieh)
    do_I(i,j) = .true.
  enddo
  ! This sets vh and dvhdv.
  call merid_flux_layer(v, h_in, h_S, h_N, &
                        vh, dvhdv, visc_rem_v_tmp, &
                        dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaV, OBC)
  ! untested
  if (local_specified_BC) then
    do concurrent (k=1:nz, j=jsh-1:jeh, i=ish:ieh, OBC%segnum_v(i,J) /= 0)
      l_seg = abs(OBC%segnum_v(i,J))
      if (OBC%segment(l_seg)%specified) vh(i,J,k) = OBC%segment(l_seg)%normal_trans(i,J,k)
    enddo
  endif
  if (present(vhbt) .or. set_BT_cont) then
    if (use_visc_rem .and. CS%use_visc_rem_max) then
      do concurrent (j=jsh-1:jeh, i=ish:ieh)
        visc_rem_max(i, j) = visc_rem_v(i,j,1)
      enddo
      do k=2,nz ; do concurrent (j=jsh-1:jeh, i=ish:ieh)
        visc_rem_max(i, j) = max(visc_rem_max(i, j), visc_rem_v(i,j,k))
      enddo ; enddo
    else
      do concurrent (j=jsh-1:jeh, i=ish:ieh)
        visc_rem_max(i, j) = 1.0
      enddo
    endif
    !   Set limits on dv that will keep the CFL number between -1 and 1.
    ! This should be adequate to keep the root bracketed in all cases.
    do concurrent (j=jsh-1:jeh, i=ish:ieh)
      I_vrm = 0.0
      if (visc_rem_max(i,j) > 0.0) I_vrm = 1.0 / visc_rem_max(i,j)
      if (CS%vol_CFL) then
        dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(i,J), 1000.0*G%dyT(i,j))
        dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(i,J), 1000.0*G%dyT(i,j+1))
      else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
      dv_max_CFL(i,j) = 2.0 * (CFL_dt * dy_S) * I_vrm
      dv_min_CFL(i,j) = -2.0 * (CFL_dt * dy_N) * I_vrm
      vh_tot_0(i,j) = 0.0 ; dvhdv_tot_0(i,j) = 0.0
    enddo
    do k=1,nz ; do concurrent (j=jsh-1:jeh, i=ish:ieh)
      dvhdv_tot_0(i,j) = dvhdv_tot_0(i,j) + dvhdv(i,j,k)
      vh_tot_0(i,j) = vh_tot_0(i,j) + vh(i,J,k)
    enddo ; enddo
    if (use_visc_rem) then
      if (CS%aggress_adjust) then
        ! untested
        do k=1,nz ; do concurrent (j=jsh-1:jeh, i=ish:ieh)
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
        enddo ; enddo
      else
        do k=1,nz ; do concurrent (j=jsh-1:jeh, i=ish:ieh)
          if (CS%vol_CFL) then
            dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(I,j), 1000.0*G%dyT(i,j))
            dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(I,j), 1000.0*G%dyT(i,j+1))
          else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
          if (dv_max_CFL(i,j) * visc_rem_v_tmp(I,j,k) > dy_S*CFL_dt - v(i,J,k)*G%mask2dCv(i,J)) &
            dv_max_CFL(i,j) = (dy_S*CFL_dt - v(i,J,k)) / visc_rem_v_tmp(I,j,k)
          if (dv_min_CFL(i,j) * visc_rem_v_tmp(I,j,k) < -dy_N*CFL_dt - v(i,J,k)*G%mask2dCv(i,J)) &
            dv_min_CFL(i,j) = -(dy_N*CFL_dt + v(i,J,k)) / visc_rem_v_tmp(I,j,k)
        enddo ; enddo
      endif
    else
      if (CS%aggress_adjust) then
        ! untested
        do k=1,nz ; do concurrent (j=jsh-1:jeh, i=ish:ieh)
          if (CS%vol_CFL) then
            dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(I,j), 1000.0*G%dyT(i,j))
            dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(I,j), 1000.0*G%dyT(i,j+1))
          else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
          dv_max_CFL(i,j) = min(dv_max_CFL(i,j), 0.499 * &
                      ((dy_S*I_dt - v(i,J,k)) + MIN(0.0,v(i,J-1,k))) )
          dv_min_CFL(i,j) = max(dv_min_CFL(i,j), 0.499 * &
                      ((-dy_N*I_dt - v(i,J,k)) + MAX(0.0,v(i,J+1,k))) )
        enddo ; enddo
      else
        do k=1,nz ; do concurrent (j=jsh-1:jeh, i=ish:ieh)
          if (CS%vol_CFL) then
            dy_S = ratio_max(G%areaT(i,j), G%dx_Cv(I,j), 1000.0*G%dyT(i,j))
            dy_N = ratio_max(G%areaT(i,j+1), G%dx_Cv(I,j), 1000.0*G%dyT(i,j+1))
          else ; dy_S = G%dyT(i,j) ; dy_N = G%dyT(i,j+1) ; endif
          dv_max_CFL(i,j) = min(dv_max_CFL(i,j), dy_S*CFL_dt - v(i,J,k))
          dv_min_CFL(i,j) = max(dv_min_CFL(i,j), -(dy_N*CFL_dt + v(i,J,k)))
        enddo ; enddo
      endif
    endif
    do concurrent (j=jsh-1:jeh, i=ish:ieh)
      dv_max_CFL(i,j) = max(dv_max_CFL(i,j),0.0)
      dv_min_CFL(i,j) = min(dv_min_CFL(i,j),0.0)
    enddo

    ! untested
    any_simple_OBC = .false.
    if (present(vhbt) .or. set_BT_cont) then
      if (local_specified_BC .or. local_Flather_OBC) then
        do concurrent (j=jsh-1:jeh, i=ish:ieh)
          l_seg = abs(OBC%segnum_v(i,J))

          ! Avoid reconciling barotropic/baroclinic transports if transport is specified
          simple_OBC_pt(i,j) = .false.
          if (l_seg /= 0) simple_OBC_pt(i,j) = OBC%segment(l_seg)%specified
          do_I(i,j) = .not.simple_OBC_pt(i,j)
          any_simple_OBC = any_simple_OBC .or. simple_OBC_pt(i,j)
        enddo
      else
        do concurrent (j=jsh-1:jeh, i=ish:ieh)
          do_I(i,j) = .true.
        enddo
      endif ! local_specified_BC .or. local_Flather_OBC
    endif ! present(vhbt) .or. set_BT_cont (redundant?)
    if (present(vhbt)) then
      ! Find dv and vh.
      call meridional_flux_adjust(v, h_in, h_S, h_N, vhbt, vh_tot_0, dvhdv_tot_0, dv, &
                             dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem_v_tmp, &
                             ish, ieh, jsh, jeh, do_I, por_face_areaV, vh, OBC=OBC)

      if (present(v_cor)) then
        do concurrent (k=1:nz, j=jsh-1:jeh, i=ish:ieh)
          v_cor(i,j,k) = v(i,j,k) + dv(i,j) * visc_rem_v_tmp(i,j,k)
        enddo
        if (any_simple_OBC) then
          ! untested
          do concurrent (k=1:nz, j=jsh-1:jeh, i=ish:ieh, simple_OBC_pt(i,j)) 
            v_cor(i,J,k) = OBC%segment(abs(OBC%segnum_v(i,J)))%normal_vel(i,J,k)
          enddo
        endif
      endif ! v-corrected

      if (present(dv_cor)) then
        do concurrent (j=jsh-1:jeh, i=ish:ieh) ; dv_cor(i,J) = dv(i,j) ; enddo
      endif ! dv-corrected
    endif
    
    if (set_BT_cont) then
      call set_merid_BT_cont(v, h_in, h_S, h_N, BT_cont, vh_tot_0, dvhdv_tot_0, &
                             dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem_v_tmp, &
                             visc_rem_max, ish, ieh, jsh, jeh, do_I, por_face_areaV)
        
      if (any_simple_OBC) then
        ! untested - these loops could probably be fused (need a test case to verify correctness)
        do concurrent (j=jsh-1:jeh, i=ish:ieh, simple_OBC_pt(i,j))
          FAvi(i,j) = GV%H_subroundoff*G%dx_Cv(i,J)
        enddo
        ! NOTE: simple_OBC_pt(i,j) should prevent access to segment OBC_NONE
        do concurrent (j=jsh-1:jeh, i=ish:ieh, simple_OBC_pt(i,j))
          segment => OBC%segment(abs(OBC%segnum_v(i,J)))
          do k=1,nz
          if ((abs(segment%normal_vel(i,J,k)) > 0.0) .and. (segment%specified)) &
            FAvi(i,j) = FAvi(i,j) + segment%normal_trans(i,J,k) / segment%normal_vel(i,J,k)
        enddo ; enddo
        do concurrent (j=jsh-1:jeh, i=ish:ieh, simple_OBC_pt(i,j))
          BT_cont%FA_v_S0(i,J) = FAvi(i,j) ; BT_cont%FA_v_N0(i,J) = FAvi(i,j)
          BT_cont%FA_v_SS(i,J) = FAvi(i,j) ; BT_cont%FA_v_NN(i,J) = FAvi(i,j)
          BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
        enddo
      endif ! any_simple_OBC
    endif ! set_BT_cont
  endif ! present(vhbt) or set_BT_cont

  ! untested - probably needs to be refactored to be performant on GPU
  if (local_open_BC .and. set_BT_cont) then
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          do concurrent (i = OBC%segment(n)%HI%Isd:OBC%segment(n)%HI%Ied)
            FA_v = 0.0
            do k=1,nz ; FA_v = FA_v + h_in(i,j,k)*(G%dx_Cv(i,J)*por_face_areaV(i,J,k)) ; enddo
            BT_cont%FA_v_S0(i,J) = FA_v ; BT_cont%FA_v_N0(i,J) = FA_v
            BT_cont%FA_v_SS(i,J) = FA_v ; BT_cont%FA_v_NN(i,J) = FA_v
            BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
          enddo
        else
          do concurrent (i = OBC%segment(n)%HI%Isd:OBC%segment(n)%HI%Ied)
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

  !$omp target exit data &
  !$omp   map(from: vh, v_cor, BT_cont%FA_v_S0, BT_cont%FA_v_SS, BT_cont%vBT_SS, BT_cont%FA_v_N0, &
  !$omp       BT_cont%FA_v_NN, BT_cont%vBT_NN, BT_cont%h_v, dv_cor) &
  !$omp   map(release: G, G%dx_Cv, G%IdyT, G%dyT, G%dyCv, G%mask2dCv, G%areaT, G%IareaT, v, h_in, &
  !$omp       h_S, h_N, CS, por_face_areaV, vhbt, visc_rem_v, dvhdv, dv, dv_min_CFL, dv_max_CFL, &
  !$omp       dvhdv_tot_0, vh_tot_0, visc_rem_max, do_I, FAvi, visc_rem_v_tmp, simple_OBC_pt)

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

  !$omp target enter data &
  !$omp   map(to: do_I, v, G, G%dx_Cv, G%IareaT, G%IdyT, h_S, h_N, h, por_face_areaV, visc_rem, &
  !$omp       vh, dvhdv)

  do concurrent (k=1:nz, j=jsh-1:jeh, i=ish:ieh, do_I(i,j))
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
  enddo

  if (local_open_BC) then
    do concurrent (k=1:nz, j=jsh-1:jeh, i=ish:ieh, do_I(i,j))
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
    enddo
  endif
  !$omp target exit data &
  !$omp   map(from: vh, dvhdv) &
  !$omp   map(release: do_I, v, G, G%dx_Cv, G%IareaT, G%IdyT, h_S, h_N, h, por_face_areaV, visc_rem)
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

  !$omp target enter data &
  !$omp   map(to: v, G, G%dx_Cv, G%IareaT, G%IdyT, h_S, h_N, h, visc_rem_v, por_face_areaV) &
  !$omp   map(alloc: h_v)

  do concurrent (k=1:nz, J=jsh-1:jeh, i=ish:ieh)
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
  enddo

  if (present(visc_rem_v)) then
    ! Scale back the thickness to account for the effects of viscosity and the fractional open
    ! thickness to give an appropriate non-normalized weight for each layer in determining the
    ! barotropic acceleration.
    do concurrent (k=1:nz, J=jsh-1:jeh, i=ish:ieh)
      h_v(i,J,k) = h_v(i,J,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
    enddo
  else
    do concurrent (k=1:nz, J=jsh-1:jeh, i=ish:ieh)
      h_v(i,J,k) = h_v(i,J,k) * por_face_areaV(i,J,k)
    enddo
  endif

  local_open_BC = .false.
  if (associated(OBC)) local_open_BC = OBC%open_v_BCs_exist_globally
  ! untested - will need to be refactored to be performant on GPUs
  if (local_open_BC) then
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          if (present(visc_rem_v)) then
            do concurrent (k=1:nz, i = OBC%segment(n)%HI%isd:OBC%segment(n)%HI%ied)
              h_v(i,J,k) = h(i,j,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
            enddo
          else
            do concurrent (k=1:nz, i = OBC%segment(n)%HI%isd:OBC%segment(n)%HI%ied)
              h_v(i,J,k) = h(i,j,k) * por_face_areaV(i,J,k)
            enddo
          endif
        else
          if (present(visc_rem_v)) then
            do concurrent (k=1:nz, i = OBC%segment(n)%HI%isd:OBC%segment(n)%HI%ied)
              h_v(i,J,k) = h(i,j+1,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
            enddo
          else
            do concurrent (k=1:nz, i = OBC%segment(n)%HI%isd:OBC%segment(n)%HI%ied)
              h_v(i,J,k) = h(i,j+1,k) * por_face_areaV(i,J,k)
            enddo
          endif
        endif
      endif
    enddo
  endif

  !$omp target exit data &
  !$omp   map(from: h_v) &
  !$omp   map(release: v, G, G%dx_Cv, G%IareaT, G%IdyT, h_S, h_N, h, visc_rem_v, por_face_areaV)

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

  !$omp target enter data &
  !$omp   map(to: G, G%IareaT, G%IdyT, G%dx_Cv, G%IdyT, v, h_in, h_S, h_N, visc_rem, vhbt, &
  !$omp       dv_max_CFL, dv_min_CFL, vh_tot_0, dvhdv_tot_0, CS, do_I_in, por_face_areaV, vh_3d) &
  !$omp   map(alloc: dv, vh_aux, dvhdv, v_new, vh_err, vh_err_best, dvhdv_tot, dv_min, dv_max, &
  !$omp       do_I)

  do concurrent (k=1:nz, j=jsh-1:jeh, i=ish:ieh)
    vh_aux(i,j,k) = 0.0 ; dvhdv(i,j,k) = 0.0
  enddo
  
  if (present(vh_3d)) then
    do concurrent (k=1:nz, j=jsh-1:jeh, i=ish:ieh)
      vh_aux(i,j,k) = vh_3d(i,J,k)
    enddo
  endif

  do concurrent (j=jsh-1:jeh, i=ish:ieh)
    dv(i,j) = 0.0 ; do_I(i,j) = do_I_in(i,j)
    dv_max(i,j) = dv_max_CFL(i,j) ; dv_min(i,j) = dv_min_CFL(i,j)
    vh_err(i,j) = vh_tot_0(i,j) - vhbt(i,j) ; dvhdv_tot(i,j) = dvhdv_tot_0(i,j)
    vh_err_best(i,j) = abs(vh_err(i,j))
  enddo

  do itt=1,max_itts
    select case (itt)
      case (:1) ; tol_eta = 1e-6 * CS%tol_eta
      case (2)  ; tol_eta = 1e-4 * CS%tol_eta
      case (3)  ; tol_eta = 1e-2 * CS%tol_eta
      case default ; tol_eta = CS%tol_eta
    end select
    tol_vel = CS%tol_vel

    do concurrent (j=jsh-1:jeh, i=ish:ieh)
      if (vh_err(i,j) > 0.0) then ; dv_max(i,j) = dv(i,j)
      elseif (vh_err(i,j) < 0.0) then ; dv_min(i,j) = dv(i,j)
      else ; do_I(i,j) = .false. ; endif
    enddo
    domore = .false.
    do concurrent (j=jsh-1:jeh, i=ish:ieh, do_I(i,j)) reduce(.or.:domore)
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
    enddo
    if (.not.domore) exit

    if ((itt < max_itts) .or. present(vh_3d)) then
      do concurrent (k=1:nz, j=jsh-1:jeh, i=ish:ieh)
        v_new(i,j,k) = v(i,J,k) + dv(i,j) * visc_rem(i,j,k)
      enddo
      call merid_flux_layer(v_new, h_in, h_S, h_N, &
                            vh_aux, dvhdv, visc_rem, &
                            dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaV, OBC)
    endif

    if (itt < max_itts) then
      do concurrent (j=jsh-1:jeh, i=ish:ieh)
        vh_err(i,j) = -vhbt(i,j) ; dvhdv_tot(i,j) = 0.0
      enddo
      do k = 1,nz ; do concurrent (j=jsh-1:jeh, i=ish:ieh)
        vh_err(i,j) = vh_err(i,j) + vh_aux(i,j,k)
        dvhdv_tot(i,j) = dvhdv_tot(i,j) + dvhdv(i,j,k)
      enddo ; enddo
      do concurrent (j=jsh-1:jeh, i=ish:ieh)
        vh_err_best(i,j) = min(vh_err_best(i,j), abs(vh_err(i,j)))
      enddo
    endif
  enddo ! itt-loop
  
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.

  if (present(vh_3d)) then
    do concurrent (k=1:nz, j=jsh-1:jeh, i=ish:ieh)
      vh_3d(i,J,k) = vh_aux(i,j,k)
    enddo
  endif

  !$omp target exit data &
  !$omp   map(from: dv, vh_3d) &
  !$omp   map(release: G, G%IareaT, G%IdyT, G%dx_Cv, G%IdyT, v, h_in, h_S, h_N, visc_rem, vhbt, &
  !$omp       dv_max_CFL, dv_min_CFL, vh_tot_0, dvhdv_tot_0, CS, do_I_in, por_face_areaV, vh_aux, &
  !$omp       dvhdv, v_new, vh_err, vh_err_best, dvhdv_tot, dv_min, dv_max, do_I)

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

  !$omp target enter data &
  !$omp   map(to: G, G%dx_Cv, G%dyCv, G%IdyT, v, h_in, h_S, h_N, BT_cont, vh_tot_0, dvhdv_tot_0, &
  !$omp       dv_max_CFL, dv_min_CFL, CS, visc_rem, visc_rem_max, do_I, por_face_areaV) &
  !$omp   map(alloc: BT_cont%FA_v_S0, BT_cont%FA_v_SS, BT_cont%vBT_SS, BT_cont%FA_v_N0, &
  !$omp       BT_cont%FA_v_NN, BT_cont%vBT_NN, dv0, dvL, dvR, dv_CFL, zeros, FAmt_L, FAmt_R, &
  !$omp       FAmt_0, vhtot_L, vhtot_R, v_L, v_R, v_0, dvhdv_L, dvhdv_R, dvhdv_0, vh_L, vh_R, vh_0)

 ! Diagnose the zero-transport correction, dv0.
  do concurrent (j=jsh-1:jeh, i=ish:ieh) ; zeros(i,j) = 0.0 ; enddo
  call meridional_flux_adjust(v, h_in, h_S, h_N, zeros, vh_tot_0, dvhdv_tot_0, dv0, &
                         dv_max_CFL, dv_min_CFL, dt, G, GV, US, CS, visc_rem, &
                         ish, ieh, jsh, jeh, do_I, por_face_areaV)

  ! Determine the southerly- and northerly- fluxes. Choose a sufficiently
  ! negative velocity correction for the northerly-flux, and a sufficiently
  ! positive correction for the southerly-flux.
  domore = .false.
  do concurrent (j=jsh-1:jeh, i=ish:ieh, do_I(i,j)) reduce(.or.:domore)
    domore = .true. ! might be better to do reduction on cpu to avoid reduce
    dv_CFL(i,j) = (CFL_min * Idt) * G%dyCv(i,J)
    dvR(i,j) = min(0.0,dv0(i,j) - dv_CFL(i,j))
    dvL(i,j) = max(0.0,dv0(i,j) + dv_CFL(i,j))
    FAmt_L(i,j) = 0.0 ; FAmt_R(i,j) = 0.0 ; FAmt_0(i,j) = 0.0
    vhtot_L(i,j) = 0.0 ; vhtot_R(i,j) = 0.0
  enddo

  if (.not.domore) then
    do concurrent (j=jsh-1:jeh, i=ish:ieh)
      BT_cont%FA_v_S0(i,J) = 0.0 ; BT_cont%FA_v_SS(i,J) = 0.0
      BT_cont%vBT_SS(i,J) = 0.0
      BT_cont%FA_v_N0(i,J) = 0.0 ; BT_cont%FA_v_NN(i,J) = 0.0
      BT_cont%vBT_NN(i,J) = 0.0
    enddo
    !$omp target exit data &
    !$omp   map(from: BT_cont%FA_v_S0, BT_cont%FA_v_SS, BT_cont%vBT_SS, BT_cont%FA_v_N0, &
    !$omp       BT_cont%FA_v_NN, BT_cont%vBT_NN) &
    !$omp   map(release: G, G%dx_Cv, G%dyCv, G%IdyT, v, h_in, h_S, h_N, BT_cont, vh_tot_0, &
    !$omp       dvhdv_tot_0, dv_max_CFL, dv_min_CFL, CS, visc_rem, visc_rem_max, do_I, &
    !$omp       por_face_areaV, dv0, dvL, dvR, dv_CFL, zeros, FAmt_L, FAmt_R, FAmt_0, vhtot_L, &
    !$omp       vhtot_R, v_L, v_R, v_0, dvhdv_L, dvhdv_R, dvhdv_0, vh_L, vh_R, vh_0)
    return
  endif

  ! not parallelized on k because of dvR/L are calculated per column
  ! nvfortran do concurrent poor performance when k is inside
  do k=1,nz ; do concurrent (j=jsh-1:jeh, i=ish:ieh, do_I(i,j))
    visc_rem_lim = max(visc_rem(i,j,k), min_visc_rem*visc_rem_max(i,j))
    if (visc_rem_lim > 0.0) then ! This is almost always true for ocean points.
      if (v(i,J,k) + dvR(i,j)*visc_rem_lim > -dv_CFL(i,j)*visc_rem(i,j,k)) &
        dvR(i,j) = -(v(i,J,k) + dv_CFL(i,j)*visc_rem(i,j,k)) / visc_rem_lim
      if (v(i,J,k) + dvL(i,j)*visc_rem_lim < dv_CFL(i,j)*visc_rem(i,j,k)) &
        dvL(i,j) = -(v(i,J,k) - dv_CFL(i,j)*visc_rem(i,j,k)) / visc_rem_lim
    endif
  enddo ; enddo

  do concurrent (k=1:nz, j=jsh-1:jeh, i=ish:ieh, do_I(i,j))
    v_L(i,j,k) = v(I,j,k) + dvL(i,j) * visc_rem(i,j,k)
    v_R(i,j,k) = v(I,j,k) + dvR(i,j) * visc_rem(i,j,k)
    v_0(i,j,k) = v(I,j,k) + dv0(i,j) * visc_rem(i,j,k)
  enddo
  call merid_flux_layer(v_0, h_in, h_S, h_N, vh_0, dvhdv_0, &
                        visc_rem, dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaV)
  call merid_flux_layer(v_L, h_in, h_S, h_N, vh_L, dvhdv_L, &
                        visc_rem, dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaV)
  call merid_flux_layer(v_R, h_in, h_S, h_N, vh_R, dvhdv_R, &
                        visc_rem, dt, G, GV, US, ish, ieh, jsh, jeh, nz, do_I, CS%vol_CFL, por_face_areaV)
  do k=1,nz ; do concurrent (j=jsh-1:jeh, i=ish:ieh, do_I(i,j))
    FAmt_0(i,j) = FAmt_0(i,j) + dvhdv_0(i,j,k)
    FAmt_L(i,j) = FAmt_L(i,j) + dvhdv_L(i,j,k)
    FAmt_R(i,j) = FAmt_R(i,j) + dvhdv_R(i,j,k)
    vhtot_L(i,j) = vhtot_L(i,j) + vh_L(i,j,k)
    vhtot_R(i,j) = vhtot_R(i,j) + vh_R(i,j,k)
  enddo ; enddo

  do concurrent (j=jsh-1:jeh, i=ish:ieh) ; if (do_I(i,j)) then
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
  endif ; enddo

  !$omp target exit data &
  !$omp   map(from: BT_cont%FA_v_S0, BT_cont%FA_v_SS, BT_cont%vBT_SS, BT_cont%FA_v_N0, &
  !$omp       BT_cont%FA_v_NN, BT_cont%vBT_NN) &
  !$omp   map(release: G, G%dx_Cv, G%dyCv, G%IdyT, v, h_in, h_S, h_N, BT_cont, vh_tot_0, &
  !$omp       dvhdv_tot_0, dv_max_CFL, dv_min_CFL, CS, visc_rem, visc_rem_max, do_I, &
  !$omp       por_face_areaV, dv0, dvL, dvR, dv_CFL, zeros, FAmt_L, FAmt_R, FAmt_0, vhtot_L, &
  !$omp       vhtot_R, v_L, v_R, v_0, dvhdv_L, dvhdv_R, dvhdv_0, vh_L, vh_R, vh_0)

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
  !$omp   map(to: G, G%mask2dT, h_in) &
  !$omp   map(alloc: h_W, h_E, slp)

  if (simple_2nd) then
    ! untested
    do concurrent (k =1:nz, j=jsl:jel, i=isl:iel)
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j,k) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j,k)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j,k) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j,k)
      h_W(i,j,k) = 0.5*( h_im1 + h_in(i,j,k) )
      h_E(i,j,k) = 0.5*( h_ip1 + h_in(i,j,k) )
    enddo
  else
    do concurrent (k=1:nz, j=jsl:jel, i=isl-1:iel+1)
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
    enddo

    if (local_open_BC) then
      ! untested
      do n=1, OBC%number_of_segments
        segment => OBC%segment(n)
        if (.not. segment%on_pe) cycle
        if (segment%is_E_or_W) then
          I=segment%HI%IsdB
          do concurrent (k=1:nz, j=segment%HI%jsd:segment%HI%jed)
            slp(i+1,j,k) = 0.0
            slp(i,j,k) = 0.0
          enddo
        endif
      enddo
    endif

    do concurrent (k=1:nz, j=jsl:jel, i=isl:iel)
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
    ! h_im1 = h_in(i-1,j,k) ; if (G%mask2dT(i-1,j) < 0.5) h_im1 = h_in(i,j)
    ! h_ip1 = h_in(i+1,j,k) ; if (G%mask2dT(i+1,j) < 0.5) h_ip1 = h_in(i,j)
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j,k) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j,k)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j,k) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j,k)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_W(i,j,k) = 0.5*( h_im1 + h_in(i,j,k) ) + oneSixth*( slp(i-1,j,k) - slp(i,j,k) )
      h_E(i,j,k) = 0.5*( h_ip1 + h_in(i,j,k) ) + oneSixth*( slp(i,j,k) - slp(i+1,j,k) )
    enddo
  endif

  if (local_open_BC) then
    ! untested
    do n=1, OBC%number_of_segments
      segment => OBC%segment(n)
      if (.not. segment%on_pe) cycle
      if (segment%direction == OBC_DIRECTION_E) then
        I=segment%HI%IsdB
        do concurrent (k=1:nz, j=segment%HI%jsd:segment%HI%jed)
          h_W(i+1,j,k) = h_in(i,j,k)
          h_E(i+1,j,k) = h_in(i,j,k)
          h_W(i,j,k) = h_in(i,j,k)
          h_E(i,j,k) = h_in(i,j,k)
        enddo
      elseif (segment%direction == OBC_DIRECTION_W) then
        I=segment%HI%IsdB
        do concurrent (k=1:nz, j=segment%HI%jsd:segment%HI%jed)
          h_W(i,j,k) = h_in(i+1,j,k)
          h_E(i,j,k) = h_in(i+1,j,k)
          h_W(i+1,j,k) = h_in(i+1,j,k)
          h_E(i+1,j,k) = h_in(i+1,j,k)
        enddo
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
  !$omp   map(from: h_W, h_E) &
  !$omp   map(release: G, G%mask2dT, h_in, slp)

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

  !$omp target enter data &
  !$omp   map(to: G, G%mask2dT, h_in) &
  !$omp   map(alloc: slp, h_S, h_N)

  if (simple_2nd) then
    ! untested
    do concurrent (k=1:nz, j=jsl:jel, i=isl:iel)
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1,k) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j,k)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1,k) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j,k)
      h_S(i,j,k) = 0.5*( h_jm1 + h_in(i,j,k) )
      h_N(i,j,k) = 0.5*( h_jp1 + h_in(i,j,k) )
    enddo
  else
    do concurrent (k=1:nz, j=jsl-1:jel+1, i=isl:iel)
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
    enddo

    if (local_open_BC) then
      ! untested
      do n=1, OBC%number_of_segments
        segment => OBC%segment(n)
        if (.not. segment%on_pe) cycle
        if (segment%is_N_or_S) then
          J=segment%HI%JsdB
          do concurrent (k=1:nz, i=segment%HI%isd:segment%HI%ied)
            slp(i,j+1,k) = 0.0
            slp(i,j,k) = 0.0
          enddo
        endif
      enddo
    endif

    do concurrent (k=1:nz, j=jsl:jel, i=isl:iel)
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1,k) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j,k)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1,k) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j,k)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_S(i,j,k) = 0.5*( h_jm1 + h_in(i,j,k) ) + oneSixth*( slp(i,j-1,k) - slp(i,j,k) )
      h_N(i,j,k) = 0.5*( h_jp1 + h_in(i,j,k) ) + oneSixth*( slp(i,j,k) - slp(i,j+1,k) )
    enddo
  endif

  if (local_open_BC) then
    ! untested
    do n=1, OBC%number_of_segments
      segment => OBC%segment(n)
      if (.not. segment%on_pe) cycle
      if (segment%direction == OBC_DIRECTION_N) then
        J=segment%HI%JsdB
        do concurrent (k=1:nz, i=segment%HI%isd:segment%HI%ied)
          h_S(i,j+1,k) = h_in(i,j,k)
          h_N(i,j+1,k) = h_in(i,j,k)
          h_S(i,j,k) = h_in(i,j,k)
          h_N(i,j,k) = h_in(i,j,k)
        enddo
      elseif (segment%direction == OBC_DIRECTION_S) then
        J=segment%HI%JsdB
        do concurrent (k=1:nz, i=segment%HI%isd:segment%HI%ied)
          h_S(i,j,k) = h_in(i,j+1,k)
          h_N(i,j,k) = h_in(i,j+1,k)
          h_S(i,j+1,k) = h_in(i,j+1,k)
          h_N(i,j+1,k) = h_in(i,j+1,k)
        enddo
      endif
    enddo
  endif

  if (monotonic) then
    ! untested
    call PPM_limit_CW84(h_in, h_S, h_N, G, GV, isl, iel, jsl, jel, nz)
  else
    call PPM_limit_pos(h_in, h_S, h_N, h_min, G, GV, isl, iel, jsl, jel, nz)
  endif

  !$omp target exit data &
  !$omp   map(release: G, G%mask2dT, h_in, slp) &
  !$omp   map(from: h_S, h_N)

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

  do concurrent (k=1:nz, j=jis:jie, i=iis:iie)
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
  enddo

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
  do concurrent (k=1:nz, j=jis:jie, i=iis:iie)
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
  enddo

end subroutine PPM_limit_CW84

!> Return the maximum ratio of a/b or maxrat.
pure function ratio_max(a, b, maxrat) result(ratio)
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
