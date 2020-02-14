!> Time steps the ocean dynamics with an unsplit quasi 2nd order Runge-Kutta scheme
module MOM_dynamics_unsplit_RK2

! This file is part of MOM6. See LICENSE.md for the license.

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Alistair Adcroft and Robert Hallberg, 2010-2012                 *
!*                                                                     *
!*    This file contains code that does the time-stepping of the       *
!*  adiabatic dynamic core, in this case with a pseudo-second order    *
!*  Runge-Kutta time stepping scheme for the momentum and a forward-   *
!*  backward coupling between the momentum and continuity equations,   *
!*  but without any splitting between the baroclinic and barotropic    *
!*  modes. Apart from the lack of splitting, this is closely analogous *
!*  to the split time stepping scheme, and efforts have been taken to  *
!*  ensure that for certain configurations (e.g., very short           *
!*  baroclinic time steps, a single barotropic step per baroclinic     *
!*  step, and particular choices about how to coupled the baroclinic   *
!*  and barotropic solves, the two solutions reproduce each other.     *
!*  Although this time stepping scheme is not very efficient with a    *
!*  large number of layers, it is valuable for verifying the proper    *
!*  behavior of the more complicated split time stepping scheme, and   *
!*  is not too inefficient for use with only a few layers.             *
!*                                                                     *
!*    The subroutine step_MOM_dyn_unsplit_RK2 actually does the time   *
!*  stepping, while register_restarts_dyn_unsplit_RK2 sets the fields  *
!*  that are found in a full restart file with this scheme, and        *
!*  initialize_dyn_unsplit_RK2 initializes the cpu clocks that are     *
!*  used in this module.  For largely historical reasons, this module  *
!*  does not have its own control structure, but shares the same       *
!*  control structure with MOM.F90 and the other MOM_dynamics_...      *
!*  modules.                                                           *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, CoriolisBu                            *
!*    j+1  > o > o >   At ^:  v, PFv, CAv, vh, diffv, tauy, vbt, vhtr  *
!*    j    x ^ x ^ x   At >:  u, PFu, CAu, uh, diffu, taux, ubt, uhtr  *
!*    j    > o > o >   At o:  h, bathyT, eta, T, S, tr                 *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1                                                  *
!*           i  i+1                                                    *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_variables, only : vertvisc_type, thermo_var_ptrs
use MOM_variables, only : ocean_internal_state, accel_diag_ptrs, cont_diag_ptrs
use MOM_forcing_type, only : mech_forcing
use MOM_checksum_packages, only : MOM_thermo_chksum, MOM_state_chksum, MOM_accel_chksum
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
use MOM_cpu_clock, only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : diag_mediator_init, enable_averages
use MOM_diag_mediator, only : disable_averaging, post_data, safe_alloc_ptr
use MOM_diag_mediator, only : register_diag_field, register_static_field
use MOM_diag_mediator, only : set_diag_mediator_grid, diag_ctrl
use MOM_domains, only : MOM_domains_init, pass_var, pass_vector
use MOM_domains, only : pass_var_start, pass_var_complete
use MOM_domains, only : pass_vector_start, pass_vector_complete
use MOM_domains, only : To_South, To_West, To_All, CGRID_NE, SCALAR_PAIR
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : MOM_set_verbosity
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_io, only : MOM_io_init
use MOM_restart, only : register_restart_field, query_initialized, save_restart
use MOM_restart, only : restart_init, MOM_restart_CS
use MOM_time_manager, only : time_type, time_type_to_real, operator(+)
use MOM_time_manager, only : operator(-), operator(>), operator(*), operator(/)

use MOM_ALE, only : ALE_CS
use MOM_boundary_update, only : update_OBC_data, update_OBC_CS
use MOM_barotropic, only : barotropic_CS
use MOM_continuity, only : continuity, continuity_init, continuity_CS, continuity_stencil
use MOM_CoriolisAdv, only : CorAdCalc, CoriolisAdv_init, CoriolisAdv_CS
use MOM_debugging, only : check_redundant
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_hor_visc, only : horizontal_viscosity, hor_visc_init, hor_visc_CS
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE_types, only : MEKE_type
use MOM_open_boundary, only : ocean_OBC_type
use MOM_open_boundary, only : radiation_open_bdry_conds
use MOM_open_boundary, only : open_boundary_zero_normal_flow
use MOM_PressureForce, only : PressureForce, PressureForce_init, PressureForce_CS
use MOM_set_visc, only : set_viscous_ML, set_visc_CS
use MOM_thickness_diffuse, only : thickness_diffuse_CS
use MOM_tidal_forcing, only : tidal_forcing_init, tidal_forcing_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_vert_friction, only : vertvisc, vertvisc_coef, vertvisc_init, vertvisc_CS
use MOM_verticalGrid, only : verticalGrid_type, get_thickness_units
use MOM_verticalGrid, only : get_flux_units, get_tr_flux_units

implicit none ; private

#include <MOM_memory.h>

!> MOM_dynamics_unsplit_RK2 module control structure
type, public :: MOM_dyn_unsplit_RK2_CS ; private
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    CAu, &    !< CAu = f*v - u.grad(u) [L T-2 ~> m s-2].
    PFu, &    !< PFu = -dM/dx [L T-2 ~> m s-2].
    diffu     !< Zonal acceleration due to convergence of the along-isopycnal stress tensor [L T-2 ~> m s-2].

  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    CAv, &    !< CAv = -f*u - u.grad(v) [L T-2 ~> m s-2].
    PFv, &    !< PFv = -dM/dy [L T-2 ~> m s-2].
    diffv     !< Meridional acceleration due to convergence of the along-isopycnal stress tensor [L T-2 ~> m s-2].

  real, pointer, dimension(:,:) :: taux_bot => NULL() !< frictional x-bottom stress from the ocean
                                                      !! to the seafloor [R L Z T-2 ~> Pa]
  real, pointer, dimension(:,:) :: tauy_bot => NULL() !< frictional y-bottom stress from the ocean
                                                      !! to the seafloor [R L Z T-2 ~> Pa]

  real    :: be      !< A nondimensional number from 0.5 to 1 that controls
                     !! the backward weighting of the time stepping scheme.
  real    :: begw    !< A nondimensional number from 0 to 1 that controls
                     !! the extent to which the treatment of gravity waves
                     !! is forward-backward (0) or simulated backward
                     !! Euler (1).  0 is almost always used.
  logical :: debug   !< If true, write verbose checksums for debugging purposes.

  logical :: module_is_initialized = .false. !< Record whether this mouled has been initialzed.

  !>@{ Diagnostic IDs
  integer :: id_uh = -1, id_vh = -1
  integer :: id_PFu = -1, id_PFv = -1, id_CAu = -1, id_CAv = -1
  !!@}

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(accel_diag_ptrs), pointer :: ADp => NULL() !< A structure pointing to the
                                   !! accelerations in the momentum equations,
                                   !! which can later be used to calculate
                                   !! derived diagnostics like energy budgets.
  type(cont_diag_ptrs), pointer :: CDp => NULL() !< A structure with pointers to
                                   !! various terms in the continuity equations,
                                   !! which can later be used to calculate
                                   !! derived diagnostics like energy budgets.

  ! The remainder of the structure points to child subroutines' control structures.
  !> A pointer to the horizontal viscosity control structure
  type(hor_visc_CS), pointer :: hor_visc_CSp => NULL()
  !> A pointer to the continuity control structure
  type(continuity_CS), pointer :: continuity_CSp => NULL()
  !> A pointer to the CoriolisAdv control structure
  type(CoriolisAdv_CS), pointer :: CoriolisAdv_CSp => NULL()
  !> A pointer to the PressureForce control structure
  type(PressureForce_CS), pointer :: PressureForce_CSp => NULL()
  !> A pointer to the vertvisc control structure
  type(vertvisc_CS), pointer :: vertvisc_CSp => NULL()
  !> A pointer to the set_visc control structure
  type(set_visc_CS), pointer :: set_visc_CSp => NULL()
  !> A pointer to the tidal forcing control structure
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL()
  !> A pointer to the ALE control structure.
  type(ALE_CS), pointer :: ALE_CSp => NULL()

  type(ocean_OBC_type), pointer :: OBC => NULL() !< A pointer to an open boundary
     !! condition type that specifies whether, where, and what open boundary
     !! conditions are used.  If no open BCs are used, this pointer stays
     !! nullified.  Flather OBCs use open boundary_CS as well.
  !> A pointer to the update_OBC control structure
  type(update_OBC_CS), pointer :: update_OBC_CSp => NULL()

end type MOM_dyn_unsplit_RK2_CS


public step_MOM_dyn_unsplit_RK2, register_restarts_dyn_unsplit_RK2
public initialize_dyn_unsplit_RK2, end_dyn_unsplit_RK2

!>@{ CPU time clock IDs
integer :: id_clock_Cor, id_clock_pres, id_clock_vertvisc
integer :: id_clock_horvisc, id_clock_continuity, id_clock_mom_update
integer :: id_clock_pass, id_clock_pass_init
!!@}

contains

! =============================================================================

!> Step the MOM6 dynamics using an unsplit quasi-2nd order Runge-Kutta scheme
subroutine step_MOM_dyn_unsplit_RK2(u_in, v_in, h_in, tv, visc, Time_local, dt, forces, &
                  p_surf_begin, p_surf_end, uh, vh, uhtr, vhtr, eta_av, G, GV, US, CS, &
                  VarMix, MEKE)
  type(ocean_grid_type),             intent(inout) :: G       !< The ocean's grid structure.
  type(verticalGrid_type),           intent(in)    :: GV      !< The ocean's vertical grid structure.
  type(unit_scale_type),             intent(in)    :: US      !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u_in !< The input and output zonal
                                                              !! velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v_in !< The input and output meridional
                                                              !! velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h_in !< The input and output layer thicknesses,
                                                              !! [H ~> m or kg m-2], depending on whether
                                                              !! the Boussinesq approximation is made.
  type(thermo_var_ptrs),             intent(in)    :: tv      !< A structure pointing to various
                                                              !! thermodynamic variables.
  type(vertvisc_type),               intent(inout) :: visc    !< A structure containing vertical
                                                              !! viscosities, bottom drag
                                                              !! viscosities, and related fields.
  type(time_type),                   intent(in)    :: Time_local   !< The model time at the end of
                                                              !! the time step.
  real,                              intent(in)    :: dt      !< The baroclinic dynamics time step [T ~> s].
  type(mech_forcing),                intent(in)    :: forces  !< A structure with the driving mechanical forces
  real, dimension(:,:),              pointer       :: p_surf_begin !< A pointer (perhaps NULL) to
                                                              !! the surface pressure at the beginning
                                                              !! of this dynamic step [Pa].
  real, dimension(:,:),              pointer       :: p_surf_end   !< A pointer (perhaps NULL) to
                                                              !! the surface pressure at the end of
                                                              !! this dynamic step [Pa].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uh !< The zonal volume or mass transport
                                                              !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vh  !< The meridional volume or mass
                                                              !! transport [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uhtr !< The accumulated zonal volume or
                                                              !! mass transport since the last
                                                              !! tracer advection [H L2 ~> m3 or kg].
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vhtr !< The accumulated meridional volume
                                                              !! or mass transport since the last
                                                              !! tracer advection [H L2 ~> m3 or kg].
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: eta_av  !< The time-mean free surface height
                                                              !! or column mass [H ~> m or kg m-2].
  type(MOM_dyn_unsplit_RK2_CS),      pointer       :: CS      !< The control structure set up by
                                                              !! initialize_dyn_unsplit_RK2.
  type(VarMix_CS),                   pointer       :: VarMix  !< A pointer to a structure with
                                                              !! fields that specify the spatially
                                                              !! variable viscosities.
  type(MEKE_type),                   pointer       :: MEKE    !< A pointer to a structure containing
                                                              !! fields related to the Mesoscale
                                                              !! Eddy Kinetic Energy.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_av, hp
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: up ! Predicted zonal velocities [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: vp ! Predicted meridional velocities [L T-1 ~> m s-1]
  real, dimension(:,:), pointer :: p_surf => NULL()
  real :: dt_pred   ! The time step for the predictor part of the baroclinic
                    ! time stepping [T ~> s].
  logical :: dyn_p_surf
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  dt_pred = dt * CS%BE

  h_av(:,:,:) = 0; hp(:,:,:) = 0
  up(:,:,:) = 0
  vp(:,:,:) = 0

  dyn_p_surf = associated(p_surf_begin) .and. associated(p_surf_end)
  if (dyn_p_surf) then
    call safe_alloc_ptr(p_surf,G%isd,G%ied,G%jsd,G%jed) ; p_surf(:,:) = 0.0
  else
    p_surf => forces%p_surf
  endif

! Runge-Kutta second order accurate two step scheme is used to step
! all of the fields except h.  h is stepped separately.

  if (CS%debug) then
    call MOM_state_chksum("Start Predictor ", u_in, v_in, h_in, uh, vh, G, GV, US)
  endif

! diffu = horizontal viscosity terms (u,h)
  call enable_averages(dt,Time_local, CS%diag)
  call cpu_clock_begin(id_clock_horvisc)
  call horizontal_viscosity(u_in, v_in, h_in, CS%diffu, CS%diffv, MEKE, VarMix, &
                            G, GV, US, CS%hor_visc_CSp)
  call cpu_clock_end(id_clock_horvisc)
  call disable_averaging(CS%diag)
  call pass_vector(CS%diffu, CS%diffv, G%Domain, clock=id_clock_pass)

! This continuity step is solely for the Coroilis terms, specifically in the
! denominator of PV and in the mass transport or PV.
! uh = u[n-1]*h[n-1/2]
! hp = h[n-1/2] + dt/2 div . uh
  call cpu_clock_begin(id_clock_continuity)
  ! This is a duplicate calculation of the last continuity from the previous step
  ! and could/should be optimized out. -AJA
  call continuity(u_in, v_in, h_in, hp, uh, vh, dt_pred, G, GV, US, &
                  CS%continuity_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_continuity)
  call pass_var(hp, G%Domain, clock=id_clock_pass)
  call pass_vector(uh, vh, G%Domain, clock=id_clock_pass)

! h_av = (h + hp)/2  (used in PV denominator)
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz
    do j=js-2,je+2 ; do i=is-2,ie+2
      h_av(i,j,k) = (h_in(i,j,k) + hp(i,j,k)) * 0.5
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_mom_update)

! CAu = -(f+zeta)/h_av vh + d/dx KE  (function of u[n-1] and uh[n-1])
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(u_in, v_in, h_av, uh, vh, CS%CAu, CS%CAv, CS%OBC, CS%ADp, &
                 G, GV, US, CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)

! PFu = d/dx M(h_av,T,S)  (function of h[n-1/2])
  call cpu_clock_begin(id_clock_pres)
  if (dyn_p_surf) then ; do j=js-2,je+2 ; do i=is-2,ie+2
    p_surf(i,j) = 0.5*p_surf_begin(i,j) + 0.5*p_surf_end(i,j)
  enddo ; enddo ; endif
  call PressureForce(h_in, tv, CS%PFu, CS%PFv, G, GV, US, &
                     CS%PressureForce_CSp, CS%ALE_CSp, p_surf)
  call cpu_clock_end(id_clock_pres)
  call pass_vector(CS%PFu, CS%PFv, G%Domain, clock=id_clock_pass)
  call pass_vector(CS%CAu, CS%CAv, G%Domain, clock=id_clock_pass)

  if (associated(CS%OBC)) then; if (CS%OBC%update_OBC) then
    call update_OBC_data(CS%OBC, G, GV, US, tv, h_in, CS%update_OBC_CSp, Time_local)
  endif; endif
  if (associated(CS%OBC)) then
    call open_boundary_zero_normal_flow(CS%OBC, G, CS%PFu, CS%PFv)
    call open_boundary_zero_normal_flow(CS%OBC, G, CS%CAu, CS%CAv)
    call open_boundary_zero_normal_flow(CS%OBC, G, CS%diffu, CS%diffv)
  endif

! up+[n-1/2] = u[n-1] + dt_pred * (PFu + CAu)
  call cpu_clock_begin(id_clock_mom_update)
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    up(I,j,k) = G%mask2dCu(I,j) * (u_in(I,j,k) + dt_pred * &
                   ((CS%PFu(I,j,k) + CS%CAu(I,j,k)) + CS%diffu(I,j,k)))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    vp(i,J,k) = G%mask2dCv(i,J) * (v_in(i,J,k) + dt_pred * &
                   ((CS%PFv(i,J,k) + CS%CAv(i,J,k)) + CS%diffv(i,J,k)))
  enddo ; enddo ; enddo
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) &
    call MOM_accel_chksum("Predictor 1 accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv,&
                          CS%diffu, CS%diffv, G, GV, US)

 ! up[n-1/2] <- up*[n-1/2] + dt/2 d/dz visc d/dz up[n-1/2]
  call cpu_clock_begin(id_clock_vertvisc)
  call enable_averages(dt, Time_local, CS%diag)
  call set_viscous_ML(up, vp, h_av, tv, forces, visc, dt_pred, G, GV, US, &
                      CS%set_visc_CSp)
  call disable_averaging(CS%diag)
  call vertvisc_coef(up, vp, h_av, forces, visc, dt_pred, G, GV, US, &
                     CS%vertvisc_CSp, CS%OBC)
  call vertvisc(up, vp, h_av, forces, visc, dt_pred, CS%OBC, CS%ADp, CS%CDp, &
                G, GV, US, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)
  call pass_vector(up, vp, G%Domain, clock=id_clock_pass)

! uh = up[n-1/2] * h[n-1/2]
! h_av = h + dt div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(up, vp, h_in, hp, uh, vh, dt, G, GV, US, CS%continuity_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_continuity)
  call pass_var(hp, G%Domain, clock=id_clock_pass)
  call pass_vector(uh, vh, G%Domain, clock=id_clock_pass)

! h_av <- (h + hp)/2   (centered at n-1/2)
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    h_av(i,j,k) = (h_in(i,j,k) + hp(i,j,k)) * 0.5
  enddo ; enddo ; enddo

  if (CS%debug) &
    call MOM_state_chksum("Predictor 1", up, vp, h_av, uh, vh, G, GV, US)

! CAu = -(f+zeta(up))/h_av vh + d/dx KE(up)  (function of up[n-1/2], h[n-1/2])
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(up, vp, h_av, uh, vh, CS%CAu, CS%CAv, CS%OBC, CS%ADp, &
                 G, GV, US, CS%CoriolisAdv_CSp)
  call cpu_clock_end(id_clock_Cor)
  if (associated(CS%OBC)) then
    call open_boundary_zero_normal_flow(CS%OBC, G, CS%CAu, CS%CAv)
  endif

! call enable_averages(dt, Time_local, CS%diag)  ?????????????????????/

! up* = u[n] + (1+gamma) * dt * ( PFu + CAu )  Extrapolated for damping
! u*[n+1] = u[n] + dt * ( PFu + CAu )
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    up(I,j,k) = G%mask2dCu(I,j) * (u_in(I,j,k) + dt * (1.+CS%begw) * &
            ((CS%PFu(I,j,k) + CS%CAu(I,j,k)) + CS%diffu(I,j,k)))
    u_in(I,j,k) = G%mask2dCu(I,j) * (u_in(I,j,k) +  dt * &
            ((CS%PFu(I,j,k) + CS%CAu(I,j,k)) + CS%diffu(I,j,k)))
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    vp(i,J,k) = G%mask2dCv(i,J) * (v_in(i,J,k) + dt * (1.+CS%begw) * &
            ((CS%PFv(i,J,k) + CS%CAv(i,J,k)) + CS%diffv(i,J,k)))
    v_in(i,J,k) = G%mask2dCv(i,J) * (v_in(i,J,k) + dt * &
            ((CS%PFv(i,J,k) + CS%CAv(i,J,k)) + CS%diffv(i,J,k)))
  enddo ; enddo ; enddo

! up[n] <- up* + dt d/dz visc d/dz up
! u[n] <- u*[n] + dt d/dz visc d/dz u[n]
  call cpu_clock_begin(id_clock_vertvisc)
  call vertvisc_coef(up, vp, h_av, forces, visc, dt, G, GV, US, &
                     CS%vertvisc_CSp, CS%OBC)
  call vertvisc(up, vp, h_av, forces, visc, dt, CS%OBC, CS%ADp, CS%CDp, &
                G, GV, US, CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot)
  call vertvisc_coef(u_in, v_in, h_av, forces, visc, dt, G, GV, US, &
                     CS%vertvisc_CSp, CS%OBC)
  call vertvisc(u_in, v_in, h_av, forces, visc, dt, CS%OBC, CS%ADp, CS%CDp,&
                G, GV, US, CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot)
  call cpu_clock_end(id_clock_vertvisc)
  call pass_vector(up, vp, G%Domain, clock=id_clock_pass)
  call pass_vector(u_in, v_in, G%Domain, clock=id_clock_pass)

! uh = up[n] * h[n]  (up[n] might be extrapolated to damp GWs)
! h[n+1] = h[n] + dt div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(up, vp, h_in, h_in, uh, vh, dt, G, GV, US, CS%continuity_CSp, OBC=CS%OBC)
  call cpu_clock_end(id_clock_continuity)
  call pass_var(h_in, G%Domain, clock=id_clock_pass)
  call pass_vector(uh, vh, G%Domain, clock=id_clock_pass)

! Accumulate mass flux for tracer transport
  do k=1,nz
    do j=js-2,je+2 ; do I=Isq-2,Ieq+2
      uhtr(I,j,k) = uhtr(I,j,k) + dt*uh(I,j,k)
    enddo ; enddo
    do J=Jsq-2,Jeq+2 ; do i=is-2,ie+2
      vhtr(i,J,k) = vhtr(i,J,k) + dt*vh(i,J,k)
    enddo ; enddo
  enddo

  if (CS%debug) then
    call MOM_state_chksum("Corrector", u_in, v_in, h_in, uh, vh, G, GV, US)
    call MOM_accel_chksum("Corrector accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv, &
                          CS%diffu, CS%diffv, G, GV, US)
  endif

  if (GV%Boussinesq) then
    do j=js,je ; do i=is,ie ; eta_av(i,j) = -GV%Z_to_H*G%bathyT(i,j) ; enddo ; enddo
  else
    do j=js,je ; do i=is,ie ; eta_av(i,j) = 0.0 ; enddo ; enddo
  endif
  do k=1,nz ; do j=js,je ; do i=is,ie
    eta_av(i,j) = eta_av(i,j) + h_av(i,j,k)
  enddo ; enddo ; enddo

  if (dyn_p_surf) deallocate(p_surf)

!   Here various terms used in to update the momentum equations are
! offered for averaging.
  if (CS%id_PFu > 0) call post_data(CS%id_PFu, CS%PFu, CS%diag)
  if (CS%id_PFv > 0) call post_data(CS%id_PFv, CS%PFv, CS%diag)
  if (CS%id_CAu > 0) call post_data(CS%id_CAu, CS%CAu, CS%diag)
  if (CS%id_CAv > 0) call post_data(CS%id_CAv, CS%CAv, CS%diag)

!   Here the thickness fluxes are offered for averaging.
  if (CS%id_uh > 0) call post_data(CS%id_uh, uh, CS%diag)
  if (CS%id_vh > 0) call post_data(CS%id_vh, vh, CS%diag)

end subroutine step_MOM_dyn_unsplit_RK2

! =============================================================================

!> Allocate the control structure for this module, allocates memory in it, and registers
!! any auxiliary restart variables that are specific to the unsplit RK2 time stepping scheme.
!!
!! All variables registered here should have the ability to be recreated if they are not present
!! in a restart file.
subroutine register_restarts_dyn_unsplit_RK2(HI, GV, param_file, CS, restart_CS)
  type(hor_index_type),         intent(in)    :: HI         !< A horizontal index type structure.
  type(verticalGrid_type),      intent(in)    :: GV         !< The ocean's vertical grid structure.
  type(param_file_type),        intent(in)    :: param_file !< A structure to parse for run-time
                                                            !! parameters.
  type(MOM_dyn_unsplit_RK2_CS), pointer       :: CS         !< The control structure set up by
                                                            !! initialize_dyn_unsplit_RK2.
  type(MOM_restart_CS),         pointer       :: restart_CS !< A pointer to the restart control
                                                            !! structure.
!   This subroutine sets up any auxiliary restart variables that are specific
! to the unsplit time stepping scheme.  All variables registered here should
! have the ability to be recreated if they are not present in a restart file.

  ! Local variables
  character(len=48) :: thickness_units, flux_units
  integer :: isd, ied, jsd, jed, nz, IsdB, IedB, JsdB, JedB
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke
  IsdB = HI%IsdB ; IedB = HI%IedB ; JsdB = HI%JsdB ; JedB = HI%JedB

! This is where a control structure that is specific to this module would be allocated.
  if (associated(CS)) then
    call MOM_error(WARNING, "register_restarts_dyn_unsplit_RK2 called with an associated "// &
                             "control structure.")
    return
  endif
  allocate(CS)

  ALLOC_(CS%diffu(IsdB:IedB,jsd:jed,nz)) ; CS%diffu(:,:,:) = 0.0
  ALLOC_(CS%diffv(isd:ied,JsdB:JedB,nz)) ; CS%diffv(:,:,:) = 0.0
  ALLOC_(CS%CAu(IsdB:IedB,jsd:jed,nz)) ; CS%CAu(:,:,:) = 0.0
  ALLOC_(CS%CAv(isd:ied,JsdB:JedB,nz)) ; CS%CAv(:,:,:) = 0.0
  ALLOC_(CS%PFu(IsdB:IedB,jsd:jed,nz)) ; CS%PFu(:,:,:) = 0.0
  ALLOC_(CS%PFv(isd:ied,JsdB:JedB,nz)) ; CS%PFv(:,:,:) = 0.0

  thickness_units = get_thickness_units(GV)
  flux_units = get_flux_units(GV)

!  No extra restart fields are needed with this time stepping scheme.

end subroutine register_restarts_dyn_unsplit_RK2

!> Initialize parameters and allocate memory associated with the unsplit RK2 dynamics module.
subroutine initialize_dyn_unsplit_RK2(u, v, h, Time, G, GV, US, param_file, diag, CS, &
                                      restart_CS, Accel_diag, Cont_diag, MIS, MEKE, &
                                      OBC, update_OBC_CSp, ALE_CSp, setVisc_CSp, &
                                      visc, dirs, ntrunc, cont_stencil)
  type(ocean_grid_type),                     intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(time_type),                   target, intent(in)    :: Time !< The current model time.
  type(param_file_type),                     intent(in)    :: param_file !< A structure to parse
                                                                         !! for run-time parameters.
  type(diag_ctrl),                   target, intent(inout) :: diag !< A structure that is used to
                                                                   !! regulate diagnostic output.
  type(MOM_dyn_unsplit_RK2_CS),              pointer       :: CS   !< The control structure set up
                                                                   !! by initialize_dyn_unsplit_RK2.
  type(MOM_restart_CS),                      pointer       :: restart_CS !< A pointer to the restart
                                                                         !! control structure.
  type(accel_diag_ptrs),             target, intent(inout) :: Accel_diag !< A set of pointers to the
                                      !! various accelerations in the momentum equations, which can
                                      !! be used for later derived diagnostics, like energy budgets.
  type(cont_diag_ptrs),              target, intent(inout) :: Cont_diag  !< A structure with pointers
                                                                         !! to various terms in the
                                                                         !! continuity equations.
  type(ocean_internal_state),                intent(inout) :: MIS  !< The "MOM6 Internal State"
                                                     !! structure, used to pass around pointers
                                                     !! to various arrays for diagnostic purposes.
  type(MEKE_type),                           pointer       :: MEKE !< MEKE data
  type(ocean_OBC_type),                      pointer       :: OBC  !< If open boundary conditions
                                                    !! are used, this points to the ocean_OBC_type
                                                    !! that was set up in MOM_initialization.
  type(update_OBC_CS),                       pointer       :: update_OBC_CSp !< If open boundary
                                                         !! condition updates are used, this points
                                                         !! to the appropriate control structure.
  type(ALE_CS),                              pointer       :: ALE_CSp     !< This points to the ALE
                                                                          !! control structure.
  type(set_visc_CS),                         pointer       :: setVisc_CSp !< This points to the
                                                                          !! set_visc control
                                                                          !! structure.
  type(vertvisc_type),                       intent(inout) :: visc !< A structure containing
                                                         !! vertical viscosities, bottom drag
                                                         !! viscosities, and related fields.
  type(directories),                         intent(in)    :: dirs !< A structure containing several
                                                                   !! relevant directory paths.
  integer, target,                           intent(inout) :: ntrunc !< A target for the variable
                                                       !! that records the number of times the
                                                       !! velocity is truncated (this should be 0).
  integer,                         optional, intent(out)   :: cont_stencil !< The stencil for
                                                       !! thickness from the continuity solver.

  !   This subroutine initializes all of the variables that are used by this
  ! dynamic core, including diagnostics and the cpu clocks.

  ! Local varaibles
  character(len=40) :: mdl = "MOM_dynamics_unsplit_RK2" ! This module's name.
  character(len=48) :: thickness_units, flux_units
  real :: H_convert
  logical :: use_tides
  integer :: isd, ied, jsd, jed, nz, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(CS)) call MOM_error(FATAL, &
      "initialize_dyn_unsplit_RK2 called with an unassociated control structure.")
  if (CS%module_is_initialized) then
    call MOM_error(WARNING, "initialize_dyn_unsplit_RK2 called with a control "// &
                            "structure that has already been initialized.")
    return
  endif
  CS%module_is_initialized = .true.

  CS%diag => diag

  call get_param(param_file, mdl, "BE", CS%be, &
                 "If SPLIT is true, BE determines the relative weighting "//&
                 "of a  2nd-order Runga-Kutta baroclinic time stepping "//&
                 "scheme (0.5) and a backward Euler scheme (1) that is "//&
                 "used for the Coriolis and inertial terms.  BE may be "//&
                 "from 0.5 to 1, but instability may occur near 0.5. "//&
                 "BE is also applicable if SPLIT is false and USE_RK2 "//&
                 "is true.", units="nondim", default=0.6)
  call get_param(param_file, mdl, "BEGW", CS%begw, &
                 "If SPLIT is true, BEGW is a number from 0 to 1 that "//&
                 "controls the extent to which the treatment of gravity "//&
                 "waves is forward-backward (0) or simulated backward "//&
                 "Euler (1).  0 is almost always used. "//&
                 "If SPLIT is false and USE_RK2 is true, BEGW can be "//&
                 "between 0 and 0.5 to damp gravity waves.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "TIDES", use_tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)

  allocate(CS%taux_bot(IsdB:IedB,jsd:jed)) ; CS%taux_bot(:,:) = 0.0
  allocate(CS%tauy_bot(isd:ied,JsdB:JedB)) ; CS%tauy_bot(:,:) = 0.0

  MIS%diffu => CS%diffu ; MIS%diffv => CS%diffv
  MIS%PFu => CS%PFu ; MIS%PFv => CS%PFv
  MIS%CAu => CS%CAu ; MIS%CAv => CS%CAv

  CS%ADp => Accel_diag ; CS%CDp => Cont_diag
  Accel_diag%diffu => CS%diffu ; Accel_diag%diffv => CS%diffv
  Accel_diag%PFu => CS%PFu ; Accel_diag%PFv => CS%PFv
  Accel_diag%CAu => CS%CAu ; Accel_diag%CAv => CS%CAv

  call continuity_init(Time, G, GV, US, param_file, diag, CS%continuity_CSp)
  if (present(cont_stencil)) cont_stencil = continuity_stencil(CS%continuity_CSp)
  call CoriolisAdv_init(Time, G, GV, US, param_file, diag, CS%ADp, CS%CoriolisAdv_CSp)
  if (use_tides) call tidal_forcing_init(Time, G, param_file, CS%tides_CSp)
  call PressureForce_init(Time, G, GV, US, param_file, diag, CS%PressureForce_CSp, &
                          CS%tides_CSp)
  call hor_visc_init(Time, G, US, param_file, diag, CS%hor_visc_CSp, MEKE)
  call vertvisc_init(MIS, Time, G, GV, US, param_file, diag, CS%ADp, dirs, &
                     ntrunc, CS%vertvisc_CSp)
  if (.not.associated(setVisc_CSp)) call MOM_error(FATAL, &
    "initialize_dyn_unsplit_RK2 called with setVisc_CSp unassociated.")
  CS%set_visc_CSp => setVisc_CSp

  if (associated(ALE_CSp)) CS%ALE_CSp => ALE_CSp
  if (associated(OBC)) CS%OBC => OBC

  flux_units = get_flux_units(GV)
  H_convert = GV%H_to_m ; if (.not.GV%Boussinesq) H_convert = GV%H_to_kg_m2
  CS%id_uh = register_diag_field('ocean_model', 'uh', diag%axesCuL, Time, &
      'Zonal Thickness Flux', flux_units, y_cell_method='sum', v_extensive=.true., &
      conversion=H_convert*US%L_to_m**2*US%s_to_T)
  CS%id_vh = register_diag_field('ocean_model', 'vh', diag%axesCvL, Time, &
      'Meridional Thickness Flux', flux_units, x_cell_method='sum', v_extensive=.true., &
      conversion=H_convert*US%L_to_m**2*US%s_to_T)
  CS%id_CAu = register_diag_field('ocean_model', 'CAu', diag%axesCuL, Time, &
      'Zonal Coriolis and Advective Acceleration', 'meter second-2', conversion=US%L_T2_to_m_s2)
  CS%id_CAv = register_diag_field('ocean_model', 'CAv', diag%axesCvL, Time, &
      'Meridional Coriolis and Advective Acceleration', 'meter second-2', conversion=US%L_T2_to_m_s2)
  CS%id_PFu = register_diag_field('ocean_model', 'PFu', diag%axesCuL, Time, &
      'Zonal Pressure Force Acceleration', 'meter second-2', conversion=US%L_T2_to_m_s2)
  CS%id_PFv = register_diag_field('ocean_model', 'PFv', diag%axesCvL, Time, &
      'Meridional Pressure Force Acceleration', 'meter second-2', conversion=US%L_T2_to_m_s2)

  id_clock_Cor = cpu_clock_id('(Ocean Coriolis & mom advection)', grain=CLOCK_MODULE)
  id_clock_continuity = cpu_clock_id('(Ocean continuity equation)', grain=CLOCK_MODULE)
  id_clock_pres = cpu_clock_id('(Ocean pressure force)', grain=CLOCK_MODULE)
  id_clock_vertvisc = cpu_clock_id('(Ocean vertical viscosity)', grain=CLOCK_MODULE)
  id_clock_horvisc = cpu_clock_id('(Ocean horizontal viscosity)', grain=CLOCK_MODULE)
  id_clock_mom_update = cpu_clock_id('(Ocean momentum increments)', grain=CLOCK_MODULE)
  id_clock_pass = cpu_clock_id('(Ocean message passing)', grain=CLOCK_MODULE)
  id_clock_pass_init = cpu_clock_id('(Ocean init message passing)', grain=CLOCK_ROUTINE)

end subroutine initialize_dyn_unsplit_RK2

!> Clean up and deallocate memory associated with the dyn_unsplit_RK2 module.
subroutine end_dyn_unsplit_RK2(CS)
  type(MOM_dyn_unsplit_RK2_CS), pointer :: CS !< dyn_unsplit_RK2 control structure that
                                              !! will be deallocated in this subroutine.

  DEALLOC_(CS%diffu) ; DEALLOC_(CS%diffv)
  DEALLOC_(CS%CAu)   ; DEALLOC_(CS%CAv)
  DEALLOC_(CS%PFu)   ; DEALLOC_(CS%PFv)

  deallocate(CS)
end subroutine end_dyn_unsplit_RK2

end module MOM_dynamics_unsplit_RK2
