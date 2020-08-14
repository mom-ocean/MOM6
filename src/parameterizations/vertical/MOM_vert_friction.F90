!> Implements vertical viscosity (vertvisc)
module MOM_vert_friction

! This file is part of MOM6. See LICENSE.md for the license.
use MOM_domains,       only : pass_var, To_All, Omit_corners
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl
use MOM_debugging, only : uvchksum, hchksum
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : mech_forcing
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_SIMPLE, OBC_NONE, OBC_DIRECTION_E
use MOM_open_boundary, only : OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_PointAccel, only : write_u_accel, write_v_accel, PointAccel_init
use MOM_PointAccel, only : PointAccel_CS
use MOM_time_manager, only : time_type, time_type_to_real, operator(-)
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs, vertvisc_type
use MOM_variables, only : cont_diag_ptrs, accel_diag_ptrs
use MOM_variables, only : ocean_internal_state
use MOM_verticalGrid, only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS
implicit none ; private

#include <MOM_memory.h>

public vertvisc, vertvisc_remnant, vertvisc_coef
public vertvisc_limit_vel, vertvisc_init, vertvisc_end
public updateCFLtruncationValue

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure with parameters and memory for the MOM_vert_friction module
type, public :: vertvisc_CS ; private
  real    :: Hmix            !< The mixed layer thickness in thickness units [H ~> m or kg m-2].
  real    :: Hmix_stress     !< The mixed layer thickness over which the wind
                             !! stress is applied with direct_stress [H ~> m or kg m-2].
  real    :: Kvml            !< The mixed layer vertical viscosity [Z2 T-1 ~> m2 s-1].
  real    :: Kv              !< The interior vertical viscosity [Z2 T-1 ~> m2 s-1].
  real    :: Hbbl            !< The static bottom boundary layer thickness [H ~> m or kg m-2].
  real    :: Kvbbl           !< The vertical viscosity in the bottom boundary
                             !! layer [Z2 T-1 ~> m2 s-1].

  real    :: maxvel          !< Velocity components greater than maxvel are truncated [L T-1 ~> m s-1].
  real    :: vel_underflow   !< Velocity components smaller than vel_underflow
                             !! are set to 0 [L T-1 ~> m s-1].
  logical :: CFL_based_trunc !< If true, base truncations on CFL numbers, not
                             !! absolute velocities.
  real    :: CFL_trunc       !< Velocity components will be truncated when they
                             !! are large enough that the corresponding CFL number
                             !! exceeds this value, nondim.
  real    :: CFL_report      !< The value of the CFL number that will cause the
                             !! accelerations to be reported, nondim.  CFL_report
                             !! will often equal CFL_trunc.
  real    :: truncRampTime   !< The time-scale over which to ramp up the value of
                             !! CFL_trunc from CFL_truncS to CFL_truncE
  real    :: CFL_truncS      !< The start value of CFL_trunc
  real    :: CFL_truncE      !< The end/target value of CFL_trunc
  logical :: CFLrampingIsActivated = .false. !< True if the ramping has been initialized
  type(time_type) :: rampStartTime !< The time at which the ramping of CFL_trunc starts

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NK_INTERFACE_) :: &
    a_u                !< The u-drag coefficient across an interface [Z T-1 ~> m s-1].
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    h_u                !< The effective layer thickness at u-points [H ~> m or kg m-2].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NK_INTERFACE_) :: &
    a_v                !< The v-drag coefficient across an interface [Z T-1 ~> m s-1].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    h_v                !< The effective layer thickness at v-points [H ~> m or kg m-2].
  real, pointer, dimension(:,:) :: a1_shelf_u => NULL() !< The u-momentum coupling coefficient under
                           !! ice shelves [Z T-1 ~> m s-1]. Retained to determine stress under shelves.
  real, pointer, dimension(:,:) :: a1_shelf_v => NULL() !< The v-momentum coupling coefficient under
                           !! ice shelves [Z T-1 ~> m s-1]. Retained to determine stress under shelves.

  logical :: split          !< If true, use the split time stepping scheme.
  logical :: bottomdraglaw  !< If true, the  bottom stress is calculated with a
                            !! drag law c_drag*|u|*u. The velocity magnitude
                            !! may be an assumed value or it may be based on the
                            !! actual velocity in the bottommost HBBL, depending
                            !! on whether linear_drag is true.
  logical :: Channel_drag   !< If true, the drag is exerted directly on each
                            !! layer according to what fraction of the bottom
                            !! they overlie.
  logical :: harmonic_visc  !< If true, the harmonic mean thicknesses are used
                            !! to calculate the viscous coupling between layers
                            !! except near the bottom.  Otherwise the arithmetic
                            !! mean thickness is used except near the bottom.
  real    :: harm_BL_val    !< A scale to determine when water is in the boundary
                            !! layers based solely on harmonic mean thicknesses
                            !! for the purpose of determining the extent to which
                            !! the thicknesses used in the viscosities are upwinded.
  logical :: direct_stress  !< If true, the wind stress is distributed over the
                            !! topmost Hmix_stress of fluid and KVML may be very small.
  logical :: dynamic_viscous_ML  !< If true, use the results from a dynamic
                            !! calculation, perhaps based on a bulk Richardson
                            !! number criterion, to determine the mixed layer
                            !! thickness for viscosity.
  logical :: answers_2018   !< If true, use the order of arithmetic and expressions that recover the
                            !! answers from the end of 2018.  Otherwise, use expressions that do not
                            !! use an arbitary and hard-coded maximum viscous coupling coefficient
                            !! between layers.
  logical :: debug          !< If true, write verbose checksums for debugging purposes.
  integer :: nkml           !< The number of layers in the mixed layer.
  integer, pointer :: ntrunc !< The number of times the velocity has been
                            !! truncated since the last call to write_energy.
  character(len=200) :: u_trunc_file  !< The complete path to a file in which a column of
                            !! u-accelerations are written if velocity truncations occur.
  character(len=200) :: v_trunc_file !< The complete path to a file in which a column of
                            !! v-accelerations are written if velocity truncations occur.
  logical :: StokesMixing   !< If true, do Stokes drift mixing via the Lagrangian current
                            !! (Eulerian plus Stokes drift).  False by default and set
                            !! via STOKES_MIXING_COMBINED.

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  !>@{ Diagnostic identifiers
  integer :: id_du_dt_visc = -1, id_dv_dt_visc = -1, id_au_vv = -1, id_av_vv = -1
  integer :: id_h_u = -1, id_h_v = -1, id_hML_u = -1 , id_hML_v = -1
  integer :: id_Ray_u = -1, id_Ray_v = -1, id_taux_bot = -1, id_tauy_bot = -1
  integer :: id_Kv_slow = -1, id_Kv_u = -1, id_Kv_v = -1
  ! integer :: id_hf_du_dt_visc    = -1, id_hf_dv_dt_visc    = -1
  integer :: id_hf_du_dt_visc_2d = -1, id_hf_dv_dt_visc_2d = -1
  !>@}

  type(PointAccel_CS), pointer :: PointAccel_CSp => NULL() !< A pointer to the control structure
                              !! for recording accelerations leading to velocity truncations

  ! real, pointer :: hf_du_dt_visc(:,:,:)  => NULL() ! Zonal friction accel. x fract. thickness [L T-2 ~> m s-2].
  ! real, pointer :: hf_dv_dt_visc(:,:,:)  => NULL() ! Merdional friction accel. x fract. thickness [L T-2 ~> m s-2].
  ! 3D diagnostics hf_du(dv)_dt_visc are commented because there is no clarity on proper remapping grid option.
  ! The code is retained for degugging purposes in the future.

end type vertvisc_CS

contains

!> Perform a fully implicit vertical diffusion
!! of momentum.  Stress top and bottom boundary conditions are used.
!!
!! This is solving the tridiagonal system
!! \f[ \left(h_k + a_{k + 1/2} + a_{k - 1/2} + r_k\right) u_k^{n+1}
!!     = h_k u_k^n + a_{k + 1/2} u_{k+1}^{n+1} + a_{k - 1/2} u_{k-1}^{n+1} \f]
!! where \f$a_{k + 1/2} = \Delta t \nu_{k + 1/2} / h_{k + 1/2}\f$
!! is the <em>interfacial coupling thickness per time step</em>,
!! encompassing background viscosity as well as contributions from
!! enhanced mixed and bottom layer viscosities.
!! $r_k$ is a Rayleight drag term due to channel drag.
!! There is an additional stress term on the right-hand side
!! if DIRECT_STRESS is true, applied to the surface layer.

subroutine vertvisc(u, v, h, forces, visc, dt, OBC, ADp, CDp, G, GV, US, CS, &
                    taux_bot, tauy_bot, Waves)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  type(mech_forcing),    intent(in)      :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),   intent(inout)   :: visc   !< Viscosities and bottom drag
  real,                  intent(in)      :: dt     !< Time increment [T ~> s]
  type(ocean_OBC_type),  pointer         :: OBC    !< Open boundary condition structure
  type(accel_diag_ptrs), intent(inout)   :: ADp    !< Accelerations in the momentum
                                                   !! equations for diagnostics
  type(cont_diag_ptrs),  intent(inout)   :: CDp    !< Continuity equation terms
  type(vertvisc_CS),     pointer         :: CS     !< Vertical viscosity control structure
  real, dimension(SZIB_(G),SZJ_(G)), &
                   optional, intent(out) :: taux_bot !< Zonal bottom stress from ocean to
                                                     !! rock [R L Z T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)), &
                   optional, intent(out) :: tauy_bot !< Meridional bottom stress from ocean to
                                                     !! rock [R L Z T-2 ~> Pa]
  type(wave_parameters_CS), &
                   optional, pointer     :: Waves !< Container for wave/Stokes information

  ! Fields from forces used in this subroutine:
  !   taux: Zonal wind stress [R L Z T-2 ~> Pa].
  !   tauy: Meridional wind stress [R L Z T-2 ~> Pa].

  ! Local variables

  real :: b1(SZIB_(G))          ! A variable used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
  real :: c1(SZIB_(G),SZK_(G))  ! A variable used by the tridiagonal solver [nondim].
  real :: d1(SZIB_(G))          ! d1=1-c1 is used by the tridiagonal solver [nondim].
  real :: Ray(SZIB_(G),SZK_(G)) ! Ray is the Rayleigh-drag velocity [Z T-1 ~> m s-1].
  real :: b_denom_1             ! The first term in the denominator of b1 [H ~> m or kg m-2].

  real :: Hmix             ! The mixed layer thickness over which stress
                           ! is applied with direct_stress [H ~> m or kg m-2].
  real :: I_Hmix           ! The inverse of Hmix [H-1 ~> m-1 or m2 kg-1].
  real :: Idt              ! The inverse of the time step [T-1 ~> s-1].
  real :: dt_Rho0          ! The time step divided by the mean density [T H Z-1 R-1 ~> s m3 kg-1 or s].
  real :: dt_Z_to_H        ! The time step times the conversion from Z to the
                           ! units of thickness - [T H Z-1 ~> s or s kg m-3].
  real :: h_neglect        ! A thickness that is so small it is usually lost
                           ! in roundoff and can be neglected [H ~> m or kg m-2].

  real :: stress           !   The surface stress times the time step, divided
                           ! by the density [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real :: zDS, hfr, h_a    ! Temporary variables used with direct_stress.
  real :: surface_stress(SZIB_(G))! The same as stress, unless the wind stress
                           ! stress is applied as a body force [H L T-1 ~> m2 s-1 or kg m-1 s-1].

  real, allocatable, dimension(:,:) :: hf_du_dt_visc_2d ! Depth sum of hf_du_dt_visc [L T-2 ~> m s-2]
  real, allocatable, dimension(:,:) :: hf_dv_dt_visc_2d ! Depth sum of hf_dv_dt_visc [L T-2 ~> m s-2]

  logical :: do_i(SZIB_(G))
  logical :: DoStokesMixing

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, n
  is = G%isc ; ie = G%iec; js = G%jsc; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = G%ke

  if (.not.associated(CS)) call MOM_error(FATAL,"MOM_vert_friction(visc): "// &
         "Module must be initialized before it is used.")

  if (CS%direct_stress) then
    Hmix = CS%Hmix_stress
    I_Hmix = 1.0 / Hmix
  endif
  dt_Rho0 = dt / GV%H_to_RZ
  dt_Z_to_H = dt*GV%Z_to_H
  h_neglect = GV%H_subroundoff
  Idt = 1.0 / dt

  !Check if Stokes mixing allowed if requested (present and associated)
  DoStokesMixing=.false.
  if (CS%StokesMixing) then
    if (present(Waves)) DoStokesMixing = associated(Waves)
    if (.not. DoStokesMixing) &
      call MOM_error(FATAL,"Stokes Mixing called without allocated"//&
                     "Waves Control Structure")
  endif

  do k=1,nz ; do i=Isq,Ieq ; Ray(i,k) = 0.0 ; enddo ; enddo

  !   Update the zonal velocity component using a modification of a standard
  ! tridagonal solver.

  !$OMP parallel do default(shared) firstprivate(Ray) &
  !$OMP                 private(do_i,surface_stress,zDS,stress,h_a,hfr, &
  !$OMP                         b_denom_1,b1,d1,c1)
  do j=G%jsc,G%jec
    do I=Isq,Ieq ; do_i(I) = (G%mask2dCu(I,j) > 0) ; enddo

    ! When mixing down Eulerian current + Stokes drift add before calling solver
    if (DoStokesMixing) then ; do k=1,nz ; do I=Isq,Ieq
      if (do_i(I)) u(I,j,k) = u(I,j,k) + US%m_s_to_L_T*Waves%Us_x(I,j,k)
    enddo ; enddo ; endif

    if (associated(ADp%du_dt_visc)) then ; do k=1,nz ; do I=Isq,Ieq
      ADp%du_dt_visc(I,j,k) = u(I,j,k)
    enddo ; enddo ; endif

!   One option is to have the wind stress applied as a body force
! over the topmost Hmix fluid.  If DIRECT_STRESS is not defined,
! the wind stress is applied as a stress boundary condition.
    if (CS%direct_stress) then
      do I=Isq,Ieq ; if (do_i(I)) then
        surface_stress(I) = 0.0
        zDS = 0.0
        stress = dt_Rho0 * forces%taux(I,j)
        do k=1,nz
          h_a = 0.5 * (h(I,j,k) + h(I+1,j,k)) + h_neglect
          hfr = 1.0 ; if ((zDS+h_a) > Hmix) hfr = (Hmix - zDS) / h_a
          u(I,j,k) = u(I,j,k) + I_Hmix * hfr * stress
          zDS = zDS + h_a ; if (zDS >= Hmix) exit
        enddo
      endif ; enddo ! end of i loop
    else ; do I=Isq,Ieq
      surface_stress(I) = dt_Rho0 * (G%mask2dCu(I,j)*forces%taux(I,j))
    enddo ; endif ! direct_stress

    if (CS%Channel_drag) then ; do k=1,nz ; do I=Isq,Ieq
      Ray(I,k) = visc%Ray_u(I,j,k)
    enddo ; enddo ; endif

    ! perform forward elimination on the tridiagonal system
    !
    ! denote the diagonal of the system as b_k, the subdiagonal as a_k
    ! and the superdiagonal as c_k. The right-hand side terms are d_k.
    !
    ! ignoring the rayleigh drag contribution,
    ! we have a_k = -dt_Z_to_H * a_u(k)
    !         b_k = h_u(k) + dt_Z_to_H * (a_u(k) + a_u(k+1))
    !         c_k = -dt_Z_to_H * a_u(k+1)
    !
    ! for forward elimination, we want to:
    ! calculate c'_k = - c_k                / (b_k + a_k c'_(k-1))
    ! and       d'_k = (d_k - a_k d'_(k-1)) / (b_k + a_k c'_(k-1))
    ! where c'_1 = c_1/b_1 and d'_1 = d_1/b_1
    ! (see Thomas' tridiagonal matrix algorithm)
    !
    ! b1 is the denominator term 1 / (b_k + a_k c'_(k-1))
    ! b_denom_1 is (b_k + a_k + c_k) - a_k(1 - c'_(k-1))
    !            = (b_k + c_k + c'_(k-1))
    ! this is done so that d1 = b1 * b_denom_1 = 1 - c'_(k-1)
    ! c1(k) is -c'_(k - 1)
    ! and the right-hand-side is destructively updated to be d'_k
    !
    do I=Isq,Ieq ; if (do_i(I)) then
      b_denom_1 = CS%h_u(I,j,1) + dt_Z_to_H * (Ray(I,1) + CS%a_u(I,j,1))
      b1(I) = 1.0 / (b_denom_1 + dt_Z_to_H*CS%a_u(I,j,2))
      d1(I) = b_denom_1 * b1(I)
      u(I,j,1) = b1(I) * (CS%h_u(I,j,1) * u(I,j,1) + surface_stress(I))
    endif ; enddo
    do k=2,nz ; do I=Isq,Ieq ; if (do_i(I)) then
      c1(I,k) = dt_Z_to_H * CS%a_u(I,j,K) * b1(I)
      b_denom_1 = CS%h_u(I,j,k) + dt_Z_to_H * (Ray(I,k) + CS%a_u(I,j,K)*d1(I))
      b1(I) = 1.0 / (b_denom_1 + dt_Z_to_H * CS%a_u(I,j,K+1))
      d1(I) = b_denom_1 * b1(I)
      u(I,j,k) = (CS%h_u(I,j,k) * u(I,j,k) + &
                  dt_Z_to_H * CS%a_u(I,j,K) * u(I,j,k-1)) * b1(I)
    endif ; enddo ; enddo

    ! back substitute to solve for the new velocities
    ! u_k = d'_k - c'_k x_(k+1)
    do k=nz-1,1,-1 ; do I=Isq,Ieq ; if (do_i(I)) then
      u(I,j,k) = u(I,j,k) + c1(I,k+1) * u(I,j,k+1)
    endif ; enddo ; enddo ! i and k loops

    if (associated(ADp%du_dt_visc)) then ; do k=1,nz ; do I=Isq,Ieq
      ADp%du_dt_visc(I,j,k) = (u(I,j,k) - ADp%du_dt_visc(I,j,k))*Idt
    enddo ; enddo ; endif

    if (associated(visc%taux_shelf)) then ; do I=Isq,Ieq
      visc%taux_shelf(I,j) = -GV%Rho0*CS%a1_shelf_u(I,j)*u(I,j,1) ! - u_shelf?
    enddo ; endif

    if (PRESENT(taux_bot)) then
      do I=Isq,Ieq
        taux_bot(I,j) = GV%Rho0 * (u(I,j,nz)*CS%a_u(I,j,nz+1))
      enddo
      if (CS%Channel_drag) then ; do k=1,nz ; do I=Isq,Ieq
        taux_bot(I,j) = taux_bot(I,j) + GV%Rho0 * (Ray(I,k)*u(I,j,k))
      enddo ; enddo ; endif
    endif

    ! When mixing down Eulerian current + Stokes drift subtract after calling solver
    if (DoStokesMixing) then ; do k=1,nz ; do I=Isq,Ieq
      if (do_i(I)) u(I,j,k) = u(I,j,k) - US%m_s_to_L_T*Waves%Us_x(I,j,k)
    enddo ; enddo ; endif

  enddo ! end u-component j loop

  ! Now work on the meridional velocity component.

  !$OMP parallel do default(shared) firstprivate(Ray) &
  !$OMP               private(do_i,surface_stress,zDS,stress,h_a,hfr, &
  !$OMP                       b_denom_1,b1,d1,c1)
  do J=Jsq,Jeq
    do i=is,ie ; do_i(i) = (G%mask2dCv(i,J) > 0) ; enddo

    ! When mixing down Eulerian current + Stokes drift add before calling solver
    if (DoStokesMixing) then ; do k=1,nz ; do i=is,ie
      if (do_i(i)) v(i,j,k) = v(i,j,k) + US%m_s_to_L_T*Waves%Us_y(i,j,k)
    enddo ; enddo ; endif

    if (associated(ADp%dv_dt_visc)) then ; do k=1,nz ; do i=is,ie
      ADp%dv_dt_visc(i,J,k) = v(i,J,k)
    enddo ; enddo ; endif

!   One option is to have the wind stress applied as a body force
! over the topmost Hmix fluid.  If DIRECT_STRESS is not defined,
! the wind stress is applied as a stress boundary condition.
    if (CS%direct_stress) then
      do i=is,ie ; if (do_i(i)) then
        surface_stress(i) = 0.0
        zDS = 0.0
        stress = dt_Rho0 * forces%tauy(i,J)
        do k=1,nz
          h_a = 0.5 * (h(i,J,k) + h(i,J+1,k)) + h_neglect
          hfr = 1.0 ; if ((zDS+h_a) > Hmix) hfr = (Hmix - zDS) / h_a
          v(i,J,k) = v(i,J,k) + I_Hmix * hfr * stress
          zDS = zDS + h_a ; if (zDS >= Hmix) exit
        enddo
      endif ; enddo ! end of i loop
    else ; do i=is,ie
      surface_stress(i) = dt_Rho0 * (G%mask2dCv(i,J)*forces%tauy(i,J))
    enddo ; endif ! direct_stress

    if (CS%Channel_drag) then ; do k=1,nz ; do i=is,ie
      Ray(i,k) = visc%Ray_v(i,J,k)
    enddo ; enddo ; endif

    do i=is,ie ; if (do_i(i)) then
      b_denom_1 = CS%h_v(i,J,1) + dt_Z_to_H * (Ray(i,1) + CS%a_v(i,J,1))
      b1(i) = 1.0 / (b_denom_1 + dt_Z_to_H*CS%a_v(i,J,2))
      d1(i) = b_denom_1 * b1(i)
      v(i,J,1) = b1(i) * (CS%h_v(i,J,1) * v(i,J,1) + surface_stress(i))
    endif ; enddo
    do k=2,nz ; do i=is,ie ; if (do_i(i)) then
      c1(i,k) = dt_Z_to_H * CS%a_v(i,J,K) * b1(i)
      b_denom_1 = CS%h_v(i,J,k) + dt_Z_to_H * (Ray(i,k) + CS%a_v(i,J,K)*d1(i))
      b1(i) = 1.0 / (b_denom_1 + dt_Z_to_H * CS%a_v(i,J,K+1))
      d1(i) = b_denom_1 * b1(i)
      v(i,J,k) = (CS%h_v(i,J,k) * v(i,J,k) + dt_Z_to_H * CS%a_v(i,J,K) * v(i,J,k-1)) * b1(i)
    endif ; enddo ; enddo
    do k=nz-1,1,-1 ; do i=is,ie ; if (do_i(i)) then
      v(i,J,k) = v(i,J,k) + c1(i,k+1) * v(i,J,k+1)
    endif ; enddo ; enddo ! i and k loops

    if (associated(ADp%dv_dt_visc)) then ; do k=1,nz ; do i=is,ie
      ADp%dv_dt_visc(i,J,k) = (v(i,J,k) - ADp%dv_dt_visc(i,J,k))*Idt
    enddo ; enddo ; endif

    if (associated(visc%tauy_shelf)) then ; do i=is,ie
      visc%tauy_shelf(i,J) = -GV%Rho0*CS%a1_shelf_v(i,J)*v(i,J,1) ! - v_shelf?
    enddo ; endif

    if (present(tauy_bot)) then
      do i=is,ie
        tauy_bot(i,J) = GV%Rho0 * (v(i,J,nz)*CS%a_v(i,J,nz+1))
      enddo
      if (CS%Channel_drag) then ; do k=1,nz ; do i=is,ie
        tauy_bot(i,J) = tauy_bot(i,J) + GV%Rho0 * (Ray(i,k)*v(i,J,k))
      enddo ; enddo ; endif
    endif

    ! When mixing down Eulerian current + Stokes drift subtract after calling solver
    if (DoStokesMixing) then ; do k=1,nz ; do i=is,ie
      if (do_i(i)) v(i,J,k) = v(i,J,k) - US%m_s_to_L_T*Waves%Us_y(i,J,k)
    enddo ; enddo ; endif

  enddo ! end of v-component J loop

  call vertvisc_limit_vel(u, v, h, ADp, CDp, forces, visc, dt, G, GV, US, CS)

  ! Here the velocities associated with open boundary conditions are applied.
  if (associated(OBC)) then
    do n=1,OBC%number_of_segments
      if (OBC%segment(n)%specified) then
        if (OBC%segment(n)%is_N_or_S) then
          J = OBC%segment(n)%HI%JsdB
          do k=1,nz ; do i=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
            v(i,J,k) = OBC%segment(n)%normal_vel(i,J,k)
          enddo ; enddo
        elseif (OBC%segment(n)%is_E_or_W) then
          I = OBC%segment(n)%HI%IsdB
          do k=1,nz ; do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
            u(I,j,k) = OBC%segment(n)%normal_vel(I,j,k)
          enddo ; enddo
        endif
      endif
    enddo
  endif

! Offer diagnostic fields for averaging.
  if (CS%id_du_dt_visc > 0) &
    call post_data(CS%id_du_dt_visc, ADp%du_dt_visc, CS%diag)
  if (CS%id_dv_dt_visc > 0) &
    call post_data(CS%id_dv_dt_visc, ADp%dv_dt_visc, CS%diag)
  if (present(taux_bot) .and. (CS%id_taux_bot > 0)) &
    call post_data(CS%id_taux_bot, taux_bot, CS%diag)
  if (present(tauy_bot) .and. (CS%id_tauy_bot > 0)) &
    call post_data(CS%id_tauy_bot, tauy_bot, CS%diag)

  ! Diagnostics for terms multiplied by fractional thicknesses

  ! 3D diagnostics hf_du(dv)_dt_visc are commented because there is no clarity on proper remapping grid option.
  ! The code is retained for degugging purposes in the future.
  !if (CS%id_hf_du_dt_visc > 0) then
  !  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
  !    CS%hf_du_dt_visc(I,j,k) = ADp%du_dt_visc(I,j,k) * ADp%diag_hfrac_u(I,j,k)
  !  enddo ; enddo ; enddo
  !  call post_data(CS%id_hf_du_dt_visc, CS%hf_du_dt_visc, CS%diag)
  !endif
  !if (CS%id_hf_dv_dt_visc > 0) then
  !  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
  !    CS%hf_dv_dt_visc(i,J,k) = ADp%dv_dt_visc(i,J,k) * ADp%diag_hfrac_v(i,J,k)
  !  enddo ; enddo ; enddo
  !  call post_data(CS%id_hf_dv_dt_visc, CS%hf_dv_dt_visc, CS%diag)
  !endif
  if (CS%id_hf_du_dt_visc_2d > 0) then
    allocate(hf_du_dt_visc_2d(G%IsdB:G%IedB,G%jsd:G%jed))
    hf_du_dt_visc_2d(:,:) = 0.0
    do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      hf_du_dt_visc_2d(I,j) = hf_du_dt_visc_2d(I,j) + ADp%du_dt_visc(I,j,k) * ADp%diag_hfrac_u(I,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_hf_du_dt_visc_2d, hf_du_dt_visc_2d, CS%diag)
    deallocate(hf_du_dt_visc_2d)
  endif
  if (CS%id_hf_dv_dt_visc_2d > 0) then
    allocate(hf_dv_dt_visc_2d(G%isd:G%ied,G%JsdB:G%JedB))
    hf_dv_dt_visc_2d(:,:) = 0.0
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      hf_dv_dt_visc_2d(i,J) = hf_dv_dt_visc_2d(i,J) + ADp%dv_dt_visc(i,J,k) * ADp%diag_hfrac_v(i,J,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_hf_dv_dt_visc_2d, hf_dv_dt_visc_2d, CS%diag)
    deallocate(hf_dv_dt_visc_2d)
  endif

end subroutine vertvisc

!> Calculate the fraction of momentum originally in a layer that remains
!! after a time-step of viscosity, and the fraction of a time-step's
!! worth of barotropic acceleration that a layer experiences after
!! viscosity is applied.
subroutine vertvisc_remnant(visc, visc_rem_u, visc_rem_v, dt, G, GV, US, CS)
  type(ocean_grid_type), intent(in)   :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV   !< Ocean vertical grid structure
  type(vertvisc_type),   intent(in)   :: visc !< Viscosities and bottom drag
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                         intent(inout) :: visc_rem_u !< Fraction of a time-step's worth of a
                                              !! barotopic acceleration that a layer experiences after
                                              !! viscosity is applied in the zonal direction [nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                         intent(inout) :: visc_rem_v !< Fraction of a time-step's worth of a
                                              !! barotopic acceleration that a layer experiences after
                                              !! viscosity is applied in the meridional direction [nondim]
  real,                  intent(in)    :: dt  !< Time increment [T ~> s]
  type(unit_scale_type), intent(in)    :: US  !< A dimensional unit scaling type
  type(vertvisc_CS),     pointer       :: CS  !< Vertical viscosity control structure

  ! Local variables

  real :: b1(SZIB_(G))          ! A variable used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
  real :: c1(SZIB_(G),SZK_(G))  ! A variable used by the tridiagonal solver [nondim].
  real :: d1(SZIB_(G))          ! d1=1-c1 is used by the tridiagonal solver [nondim].
  real :: Ray(SZIB_(G),SZK_(G)) ! Ray is the Rayleigh-drag velocity [Z T-1 ~> m s-1].
  real :: b_denom_1   ! The first term in the denominator of b1 [H ~> m or kg m-2].
  real :: dt_Z_to_H        ! The time step times the conversion from Z to the
                           ! units of thickness [T H Z-1 ~> s or s kg m-3].
  logical :: do_i(SZIB_(G))

  integer :: i, j, k, is, ie, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = G%ke

  if (.not.associated(CS)) call MOM_error(FATAL,"MOM_vert_friction(visc): "// &
         "Module must be initialized before it is used.")

  dt_Z_to_H = dt*GV%Z_to_H

  do k=1,nz ; do i=Isq,Ieq ; Ray(i,k) = 0.0 ; enddo ; enddo

  ! Find the zonal viscous using a modification of a standard tridagonal solver.
!$OMP parallel do default(none) shared(G,Isq,Ieq,CS,nz,visc,dt_Z_to_H,visc_rem_u) &
!$OMP                     firstprivate(Ray)                                       &
!$OMP                          private(do_i,b_denom_1,b1,d1,c1)
  do j=G%jsc,G%jec
    do I=Isq,Ieq ; do_i(I) = (G%mask2dCu(I,j) > 0) ; enddo

    if (CS%Channel_drag) then ; do k=1,nz ; do I=Isq,Ieq
      Ray(I,k) = visc%Ray_u(I,j,k)
    enddo ; enddo ; endif

    do I=Isq,Ieq ; if (do_i(I)) then
      b_denom_1 = CS%h_u(I,j,1) + dt_Z_to_H * (Ray(I,1) + CS%a_u(I,j,1))
      b1(I) = 1.0 / (b_denom_1 + dt_Z_to_H*CS%a_u(I,j,2))
      d1(I) = b_denom_1 * b1(I)
      visc_rem_u(I,j,1) = b1(I) * CS%h_u(I,j,1)
    endif ; enddo
    do k=2,nz ; do I=Isq,Ieq ; if (do_i(I)) then
      c1(I,k) = dt_Z_to_H * CS%a_u(I,j,K)*b1(I)
      b_denom_1 = CS%h_u(I,j,k) + dt_Z_to_H * (Ray(I,k) + CS%a_u(I,j,K)*d1(I))
      b1(I) = 1.0 / (b_denom_1 + dt_Z_to_H * CS%a_u(I,j,K+1))
      d1(I) = b_denom_1 * b1(I)
      visc_rem_u(I,j,k) = (CS%h_u(I,j,k) + dt_Z_to_H * CS%a_u(I,j,K) * visc_rem_u(I,j,k-1)) * b1(I)
    endif ; enddo ; enddo
    do k=nz-1,1,-1 ; do I=Isq,Ieq ; if (do_i(I)) then
      visc_rem_u(I,j,k) = visc_rem_u(I,j,k) + c1(I,k+1)*visc_rem_u(I,j,k+1)

    endif ; enddo ; enddo ! i and k loops

  enddo ! end u-component j loop

  ! Now find the meridional viscous using a modification.
!$OMP parallel do default(none) shared(Jsq,Jeq,is,ie,G,CS,visc,dt_Z_to_H,visc_rem_v,nz) &
!$OMP                     firstprivate(Ray)                                             &
!$OMP                          private(do_i,b_denom_1,b1,d1,c1)
  do J=Jsq,Jeq
    do i=is,ie ; do_i(i) = (G%mask2dCv(i,J) > 0) ; enddo

    if (CS%Channel_drag) then ; do k=1,nz ; do i=is,ie
      Ray(i,k) = visc%Ray_v(i,J,k)
    enddo ; enddo ; endif

    do i=is,ie ; if (do_i(i)) then
      b_denom_1 = CS%h_v(i,J,1) + dt_Z_to_H * (Ray(i,1) + CS%a_v(i,J,1))
      b1(i) = 1.0 / (b_denom_1 + dt_Z_to_H*CS%a_v(i,J,2))
      d1(i) = b_denom_1 * b1(i)
      visc_rem_v(i,J,1) = b1(i) * CS%h_v(i,J,1)
    endif ; enddo
    do k=2,nz ; do i=is,ie ; if (do_i(i)) then
      c1(i,k) = dt_Z_to_H * CS%a_v(i,J,K)*b1(i)
      b_denom_1 = CS%h_v(i,J,k) + dt_Z_to_H * (Ray(i,k) + CS%a_v(i,J,K)*d1(i))
      b1(i) = 1.0 / (b_denom_1 + dt_Z_to_H * CS%a_v(i,J,K+1))
      d1(i) = b_denom_1 * b1(i)
      visc_rem_v(i,J,k) = (CS%h_v(i,J,k) + dt_Z_to_H * CS%a_v(i,J,K) * visc_rem_v(i,J,k-1)) * b1(i)
    endif ; enddo ; enddo
    do k=nz-1,1,-1 ; do i=is,ie ; if (do_i(i)) then
      visc_rem_v(i,J,k) = visc_rem_v(i,J,k) + c1(i,k+1)*visc_rem_v(i,J,k+1)
    endif ; enddo ; enddo ! i and k loops
  enddo ! end of v-component J loop

  if (CS%debug) then
    call uvchksum("visc_rem_[uv]", visc_rem_u, visc_rem_v, G%HI, haloshift=0, &
                  scalar_pair=.true.)
  endif

end subroutine vertvisc_remnant


!> Calculate the coupling coefficients (CS%a_u and CS%a_v)
!! and effective layer thicknesses (CS%h_u and CS%h_v) for later use in the
!! applying the implicit vertical viscosity via vertvisc().
subroutine vertvisc_coef(u, v, h, forces, visc, dt, G, GV, US, CS, OBC)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  type(mech_forcing),      intent(in)    :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),     intent(in)    :: visc   !< Viscosities and bottom drag
  real,                    intent(in)    :: dt     !< Time increment [T ~> s]
  type(vertvisc_CS),       pointer       :: CS     !< Vertical viscosity control structure
  type(ocean_OBC_type),    pointer       :: OBC    !< Open boundary condition structure

  ! Field from forces used in this subroutine:
  !   ustar: the friction velocity [Z T-1 ~> m s-1], used here as the mixing
  !     velocity in the mixed layer if NKML > 1 in a bulk mixed layer.

  ! Local variables

  real, dimension(SZIB_(G),SZK_(G)) :: &
    h_harm, &   ! Harmonic mean of the thicknesses around a velocity grid point,
                ! given by 2*(h+ * h-)/(h+ + h-) [H ~> m or kg m-2].
    h_arith, &  ! The arithmetic mean thickness [H ~> m or kg m-2].
    h_delta, &  ! The lateral difference of thickness [H ~> m or kg m-2].
    hvel, &     ! hvel is the thickness used at a velocity grid point [H ~> m or kg m-2].
    hvel_shelf  ! The equivalent of hvel under shelves [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZK_(G)+1) :: &
    a_cpl, &    ! The drag coefficients across interfaces [Z T-1 ~> m s-1].  a_cpl times
                ! the velocity difference gives the stress across an interface.
    a_shelf, &  ! The drag coefficients across interfaces in water columns under
                ! ice shelves [Z T-1 ~> m s-1].
    z_i         ! An estimate of each interface's height above the bottom,
                ! normalized by the bottom boundary layer thickness, nondim.
  real, dimension(SZIB_(G)) :: &
    kv_bbl, &     ! The bottom boundary layer viscosity [Z2 T-1 ~> m2 s-1].
    bbl_thick, &  ! The bottom boundary layer thickness [H ~> m or kg m-2].
    I_Hbbl, &     ! The inverse of the bottom boundary layer thickness [H-1 ~> m-1 or m2 kg-1].
    I_Htbl, &     ! The inverse of the top boundary layer thickness [H-1 ~> m-1 or m2 kg-1].
    zcol1, &      ! The height of the interfaces to the north and south of a
    zcol2, &      ! v-point [H ~> m or kg m-2].
    Ztop_min, &   ! The deeper of the two adjacent surface heights [H ~> m or kg m-2].
    Dmin, &       ! The shallower of the two adjacent bottom depths converted to
                  ! thickness units [H ~> m or kg m-2].
    zh, &         ! An estimate of the interface's distance from the bottom
                  ! based on harmonic mean thicknesses [H ~> m or kg m-2].
    h_ml          ! The mixed layer depth [H ~> m or kg m-2].
  real, allocatable, dimension(:,:) :: hML_u ! Diagnostic of the mixed layer depth at u points [H ~> m or kg m-2].
  real, allocatable, dimension(:,:) :: hML_v ! Diagnostic of the mixed layer depth at v points [H ~> m or kg m-2].
  real, allocatable, dimension(:,:,:) :: Kv_u !< Total vertical viscosity at u-points [Z2 T-1 ~> m2 s-1].
  real, allocatable, dimension(:,:,:) :: Kv_v !< Total vertical viscosity at v-points [Z2 T-1 ~> m2 s-1].
  real :: zcol(SZI_(G)) ! The height of an interface at h-points [H ~> m or kg m-2].
  real :: botfn   ! A function which goes from 1 at the bottom to 0 much more
                  ! than Hbbl into the interior.
  real :: topfn   ! A function which goes from 1 at the top to 0 much more
                  ! than Htbl into the interior.
  real :: z2      ! The distance from the bottom, normalized by Hbbl, nondim.
  real :: z2_wt   ! A nondimensional (0-1) weight used when calculating z2.
  real :: z_clear ! The clearance of an interface above the surrounding topography [H ~> m or kg m-2].
  real :: a_cpl_max  ! The maximum drag doefficient across interfaces, set so that it will be
                     ! representable as a 32-bit float in MKS units  [Z T-1 ~> m s-1]
  real :: h_neglect  ! A thickness that is so small it is usually lost
                     ! in roundoff and can be neglected [H ~> m or kg m-2].

  real :: I_valBL ! The inverse of a scaling factor determining when water is
                  ! still within the boundary layer, as determined by the sum
                  ! of the harmonic mean thicknesses.
  logical, dimension(SZIB_(G)) :: do_i, do_i_shelf
  logical :: do_any_shelf
  integer, dimension(SZIB_(G)) :: &
    zi_dir   !  A trinary logical array indicating which thicknesses to use for
             !  finding z_clear.
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = G%ke

  if (.not.associated(CS)) call MOM_error(FATAL,"MOM_vert_friction(coef): "// &
         "Module must be initialized before it is used.")

  h_neglect = GV%H_subroundoff
  a_cpl_max = 1.0e37 * US%m_to_Z * US%T_to_s
  I_Hbbl(:) = 1.0 / (CS%Hbbl + h_neglect)
  I_valBL = 0.0 ; if (CS%harm_BL_val > 0.0) I_valBL = 1.0 / CS%harm_BL_val

  if (CS%id_Kv_u > 0) then
    allocate(Kv_u(G%IsdB:G%IedB,G%jsd:G%jed,G%ke)) ; Kv_u(:,:,:) = 0.0
  endif

  if (CS%id_Kv_v > 0) then
    allocate(Kv_v(G%isd:G%ied,G%JsdB:G%JedB,G%ke)) ; Kv_v(:,:,:) = 0.0
  endif

  if (CS%debug .or. (CS%id_hML_u > 0)) then
    allocate(hML_u(G%IsdB:G%IedB,G%jsd:G%jed)) ; hML_u(:,:) = 0.0
  endif
  if (CS%debug .or. (CS%id_hML_v > 0)) then
    allocate(hML_v(G%isd:G%ied,G%JsdB:G%JedB)) ; hML_v(:,:) = 0.0
  endif

  if ((associated(visc%taux_shelf) .or. associated(forces%frac_shelf_u)) .and. &
      .not.associated(CS%a1_shelf_u)) then
    allocate(CS%a1_shelf_u(G%IsdB:G%IedB,G%jsd:G%jed)) ; CS%a1_shelf_u(:,:)=0.0
  endif
  if ((associated(visc%tauy_shelf) .or. associated(forces%frac_shelf_v)) .and. &
      .not.associated(CS%a1_shelf_v)) then
    allocate(CS%a1_shelf_v(G%isd:G%ied,G%JsdB:G%JedB)) ; CS%a1_shelf_v(:,:)=0.0
  endif

  !$OMP parallel do default(private) shared(G,GV,US,CS,visc,Isq,Ieq,nz,u,h,forces,hML_u, &
  !$OMP                                     OBC,h_neglect,dt,I_valBL,Kv_u,a_cpl_max) &
  !$OMP                     firstprivate(i_hbbl)
  do j=G%Jsc,G%Jec
    do I=Isq,Ieq ; do_i(I) = (G%mask2dCu(I,j) > 0) ; enddo

    if (CS%bottomdraglaw) then ; do I=Isq,Ieq
      kv_bbl(I) = visc%Kv_bbl_u(I,j)
      bbl_thick(I) = visc%bbl_thick_u(I,j) * GV%Z_to_H + h_neglect
      if (do_i(I)) I_Hbbl(I) = 1.0 / bbl_thick(I)
    enddo ; endif

    do k=1,nz ; do I=Isq,Ieq ; if (do_i(I)) then
      h_harm(I,k) = 2.0*h(i,j,k)*h(i+1,j,k) / (h(i,j,k)+h(i+1,j,k)+h_neglect)
      h_arith(I,k) = 0.5*(h(i+1,j,k)+h(i,j,k))
      h_delta(I,k) = h(i+1,j,k) - h(i,j,k)
    endif ; enddo ; enddo
    do I=Isq,Ieq
      Dmin(I) = min(G%bathyT(i,j), G%bathyT(i+1,j)) * GV%Z_to_H
      zi_dir(I) = 0
    enddo

    ! Project thickness outward across OBCs using a zero-gradient condition.
    if (associated(OBC)) then ; if (OBC%number_of_segments > 0) then
      do I=Isq,Ieq ; if (do_i(I) .and. (OBC%segnum_u(I,j) /= OBC_NONE)) then
        if (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_E) then
          do k=1,nz ; h_harm(I,k) = h(i,j,k) ; h_arith(I,k) = h(i,j,k) ; h_delta(I,k) = 0. ; enddo
          Dmin(I) = G%bathyT(i,j) * GV%Z_to_H
          zi_dir(I) = -1
        elseif (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_W) then
          do k=1,nz ; h_harm(I,k) = h(i+1,j,k) ; h_arith(I,k) = h(i+1,j,k) ; h_delta(I,k) = 0. ; enddo
          Dmin(I) = G%bathyT(i+1,j) * GV%Z_to_H
          zi_dir(I) = 1
        endif
      endif ; enddo
    endif ; endif

!    The following block calculates the thicknesses at velocity
!  grid points for the vertical viscosity (hvel).  Near the
!  bottom an upwind biased thickness is used to control the effect
!  of spurious Montgomery potential gradients at the bottom where
!  nearly massless layers layers ride over the topography.
    if (CS%harmonic_visc) then
      do I=Isq,Ieq ; z_i(I,nz+1) = 0.0 ; enddo
      do k=nz,1,-1 ; do I=Isq,Ieq ; if (do_i(I)) then
        hvel(I,k) = h_harm(I,k)
        if (u(I,j,k) * h_delta(I,k) < 0) then
          z2 = z_i(I,k+1) ; botfn = 1.0 / (1.0 + 0.09*z2*z2*z2*z2*z2*z2)
          hvel(I,k) = (1.0-botfn)*h_harm(I,k) + botfn*h_arith(I,k)
        endif
        z_i(I,k) =  z_i(I,k+1) + h_harm(I,k)*I_Hbbl(I)
      endif ; enddo ; enddo ! i & k loops
    else ! Not harmonic_visc
      do I=Isq,Ieq ; zh(I) = 0.0 ; z_i(I,nz+1) = 0.0 ; enddo
      do i=Isq,Ieq+1 ; zcol(i) = -G%bathyT(i,j) * GV%Z_to_H ; enddo
      do k=nz,1,-1
        do i=Isq,Ieq+1 ; zcol(i) = zcol(i) + h(i,j,k) ; enddo
        do I=Isq,Ieq ; if (do_i(I)) then
          zh(I) = zh(I) + h_harm(I,k)

          z_clear = max(zcol(i),zcol(i+1)) + Dmin(I)
          if (zi_dir(I) < 0) z_clear = zcol(i) + Dmin(I)
          if (zi_dir(I) > 0) z_clear = zcol(i+1) + Dmin(I)

          z_i(I,k) = max(zh(I), z_clear) * I_Hbbl(I)

          hvel(I,k) = h_arith(I,k)
          if (u(I,j,k) * h_delta(I,k) > 0) then
            if (zh(I) * I_Hbbl(I) < CS%harm_BL_val) then
              hvel(I,k) = h_harm(I,k)
            else
              z2_wt = 1.0  ; if (zh(I) * I_Hbbl(I) < 2.0*CS%harm_BL_val) &
                z2_wt = max(0.0, min(1.0, zh(I) * I_Hbbl(I) * I_valBL - 1.0))
              z2 = z2_wt * (max(zh(I), z_clear) * I_Hbbl(I))
              botfn = 1.0 / (1.0 + 0.09*z2*z2*z2*z2*z2*z2)
              hvel(I,k) = (1.0-botfn)*h_arith(I,k) + botfn*h_harm(I,k)
            endif
          endif

        endif ; enddo ! i  loop
      enddo ! k loop
    endif

    call find_coupling_coef(a_cpl, hvel, do_i, h_harm, bbl_thick, kv_bbl, z_i, h_ml, &
                            dt, j, G, GV, US, CS, visc, forces, work_on_u=.true., OBC=OBC)
    if (allocated(hML_u)) then
      do i=isq,ieq ; if (do_i(i)) then ; hML_u(I,j) = h_ml(I) ; endif ; enddo
    endif

    do_any_shelf = .false.
    if (associated(forces%frac_shelf_u)) then
      do I=Isq,Ieq
        CS%a1_shelf_u(I,j) = 0.0
        do_i_shelf(I) = (do_i(I) .and. forces%frac_shelf_u(I,j) > 0.0)
        if (do_i_shelf(I)) do_any_shelf = .true.
      enddo
      if (do_any_shelf) then
        if (CS%harmonic_visc) then
          call find_coupling_coef(a_shelf, hvel, do_i_shelf, h_harm, bbl_thick, &
                                  kv_bbl, z_i, h_ml, dt, j, G, GV, US, CS, &
                                  visc, forces, work_on_u=.true., OBC=OBC, shelf=.true.)
        else  ! Find upwind-biased thickness near the surface.
          ! Perhaps this needs to be done more carefully, via find_eta.
          do I=Isq,Ieq ; if (do_i_shelf(I)) then
            zh(I) = 0.0 ; Ztop_min(I) = min(zcol(i), zcol(i+1))
            I_HTbl(I) = 1.0 / (visc%tbl_thick_shelf_u(I,j)*GV%Z_to_H + h_neglect)
          endif ; enddo
          do k=1,nz
            do i=Isq,Ieq+1 ; zcol(i) = zcol(i) - h(i,j,k) ; enddo
            do I=Isq,Ieq ; if (do_i_shelf(I)) then
              zh(I) = zh(I) + h_harm(I,k)

              hvel_shelf(I,k) = hvel(I,k)
              if (u(I,j,k) * h_delta(I,k) > 0) then
                if (zh(I) * I_HTbl(I) < CS%harm_BL_val) then
                  hvel_shelf(I,k) = min(hvel(I,k), h_harm(I,k))
                else
                  z2_wt = 1.0  ; if (zh(I) * I_HTbl(I) < 2.0*CS%harm_BL_val) &
                    z2_wt = max(0.0, min(1.0, zh(I) * I_HTbl(I) * I_valBL - 1.0))
                  z2 = z2_wt * (max(zh(I), Ztop_min(I) - min(zcol(i),zcol(i+1))) * I_HTbl(I))
                  topfn = 1.0 / (1.0 + 0.09*z2**6)
                  hvel_shelf(I,k) = min(hvel(I,k), (1.0-topfn)*h_arith(I,k) + topfn*h_harm(I,k))
                endif
              endif
            endif ; enddo
          enddo
          call find_coupling_coef(a_shelf, hvel_shelf, do_i_shelf, h_harm, &
                                  bbl_thick, kv_bbl, z_i, h_ml, dt, j, G, GV, US, CS, &
                                  visc, forces, work_on_u=.true., OBC=OBC, shelf=.true.)
        endif
        do I=Isq,Ieq ; if (do_i_shelf(I)) CS%a1_shelf_u(I,j) = a_shelf(I,1) ; enddo
      endif
    endif

    if (do_any_shelf) then
      do K=1,nz+1 ; do I=Isq,Ieq ; if (do_i_shelf(I)) then
        CS%a_u(I,j,K) = min(a_cpl_max, forces%frac_shelf_u(I,j)  * a_shelf(I,K) + &
                                       (1.0-forces%frac_shelf_u(I,j)) * a_cpl(I,K))
! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
!        CS%a_u(I,j,K) = min(a_cpl_max, forces%frac_shelf_u(I,j)  * max(a_shelf(I,K), a_cpl(I,K)) + &
!                                       (1.0-forces%frac_shelf_u(I,j)) * a_cpl(I,K))
      elseif (do_i(I)) then
        CS%a_u(I,j,K) = min(a_cpl_max, a_cpl(I,K))
      endif ; enddo ; enddo
      do k=1,nz ; do I=Isq,Ieq ; if (do_i_shelf(I)) then
        ! Should we instead take the inverse of the average of the inverses?
        CS%h_u(I,j,k) = forces%frac_shelf_u(I,j)  * hvel_shelf(I,k) + &
                   (1.0-forces%frac_shelf_u(I,j)) * hvel(I,k) + h_neglect
      elseif (do_i(I)) then
        CS%h_u(I,j,k) = hvel(I,k) + h_neglect
      endif ; enddo ; enddo
    else
      do K=1,nz+1 ; do I=Isq,Ieq ; if (do_i(I)) CS%a_u(I,j,K) = min(a_cpl_max, a_cpl(I,K)) ; enddo ; enddo
      do k=1,nz ; do I=Isq,Ieq ; if (do_i(I)) CS%h_u(I,j,k) = hvel(I,k) + h_neglect ; enddo ; enddo
    endif

    ! Diagnose total Kv at u-points
    if (CS%id_Kv_u > 0) then
      do k=1,nz ; do I=Isq,Ieq
        if (do_i(I)) Kv_u(I,j,k) = 0.5 * GV%H_to_Z*(CS%a_u(I,j,K)+CS%a_u(I,j,K+1)) * CS%h_u(I,j,k)
      enddo ; enddo
    endif

  enddo


  ! Now work on v-points.
  !$OMP parallel do default(private) shared(G,GV,CS,US,visc,is,ie,Jsq,Jeq,nz,v,h,forces,hML_v, &
  !$OMP                                  OBC,h_neglect,dt,I_valBL,Kv_v,a_cpl_max) &
  !$OMP                     firstprivate(i_hbbl)
  do J=Jsq,Jeq
    do i=is,ie ; do_i(i) = (G%mask2dCv(i,J) > 0) ; enddo

    if (CS%bottomdraglaw) then ; do i=is,ie
      kv_bbl(i) = visc%Kv_bbl_v(i,J)
      bbl_thick(i) = visc%bbl_thick_v(i,J) * GV%Z_to_H + h_neglect
      if (do_i(i)) I_Hbbl(i) = 1.0 / bbl_thick(i)
    enddo ; endif

    do k=1,nz ; do i=is,ie ; if (do_i(i)) then
      h_harm(i,k) = 2.0*h(i,j,k)*h(i,j+1,k) / (h(i,j,k)+h(i,j+1,k)+h_neglect)
      h_arith(i,k) = 0.5*(h(i,j+1,k)+h(i,j,k))
      h_delta(i,k) = h(i,j+1,k) - h(i,j,k)
    endif ; enddo ; enddo
    do i=is,ie
      Dmin(i) = min(G%bathyT(i,j), G%bathyT(i,j+1)) * GV%Z_to_H
      zi_dir(i) = 0
    enddo

    ! Project thickness outward across OBCs using a zero-gradient condition.
    if (associated(OBC)) then ; if (OBC%number_of_segments > 0) then
      do i=is,ie ; if (do_i(i) .and. (OBC%segnum_v(i,J) /= OBC_NONE)) then
        if (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_N) then
          do k=1,nz ; h_harm(I,k) = h(i,j,k) ; h_arith(I,k) = h(i,j,k) ; h_delta(i,k) = 0. ; enddo
          Dmin(I) = G%bathyT(i,j) * GV%Z_to_H
          zi_dir(I) = -1
        elseif (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_S) then
          do k=1,nz ; h_harm(i,k) = h(i,j+1,k) ; h_arith(i,k) = h(i,j+1,k) ; h_delta(i,k) = 0. ; enddo
          Dmin(i) = G%bathyT(i,j+1) * GV%Z_to_H
          zi_dir(i) = 1
        endif
      endif ; enddo
    endif ; endif

!    The following block calculates the thicknesses at velocity
!  grid points for the vertical viscosity (hvel).  Near the
!  bottom an upwind biased thickness is used to control the effect
!  of spurious Montgomery potential gradients at the bottom where
!  nearly massless layers layers ride over the topography.
    if (CS%harmonic_visc) then
      do i=is,ie ; z_i(i,nz+1) = 0.0 ; enddo

      do k=nz,1,-1 ; do i=is,ie ; if (do_i(i)) then
        hvel(i,k) = h_harm(i,k)
        if (v(i,J,k) * h_delta(i,k) < 0) then
          z2 = z_i(i,k+1) ; botfn = 1.0 / (1.0 + 0.09*z2*z2*z2*z2*z2*z2)
          hvel(i,k) = (1.0-botfn)*h_harm(i,k) + botfn*h_arith(i,k)
        endif
        z_i(i,k) = z_i(i,k+1)  + h_harm(i,k)*I_Hbbl(i)
      endif ; enddo ; enddo ! i & k loops
    else ! Not harmonic_visc
      do i=is,ie
        zh(i) = 0.0 ; z_i(i,nz+1) = 0.0
        zcol1(i) = -G%bathyT(i,j) * GV%Z_to_H
        zcol2(i) = -G%bathyT(i,j+1) * GV%Z_to_H
      enddo
      do k=nz,1,-1 ; do i=is,ie ; if (do_i(i)) then
        zh(i) = zh(i) + h_harm(i,k)
        zcol1(i) = zcol1(i) + h(i,j,k) ; zcol2(i) = zcol2(i) + h(i,j+1,k)

        z_clear = max(zcol1(i),zcol2(i)) + Dmin(i)
        if (zi_dir(i) < 0) z_clear = zcol1(i) + Dmin(I)
        if (zi_dir(i) > 0) z_clear = zcol2(i) + Dmin(I)

        z_i(I,k) = max(zh(i), z_clear) * I_Hbbl(i)

        hvel(i,k) = h_arith(i,k)
        if (v(i,J,k) * h_delta(i,k) > 0) then
          if (zh(i) * I_Hbbl(i) < CS%harm_BL_val) then
            hvel(i,k) = h_harm(i,k)
          else
            z2_wt = 1.0  ; if (zh(i) * I_Hbbl(i) < 2.0*CS%harm_BL_val) &
              z2_wt = max(0.0, min(1.0, zh(i) * I_Hbbl(i) * I_valBL - 1.0))
            z2 = z2_wt * (max(zh(i), max(zcol1(i),zcol2(i)) + Dmin(i)) * I_Hbbl(i))
            botfn = 1.0 / (1.0 + 0.09*z2*z2*z2*z2*z2*z2)
            hvel(i,k) = (1.0-botfn)*h_arith(i,k) + botfn*h_harm(i,k)
          endif
       endif

      endif ; enddo ; enddo ! i & k loops
    endif

    call find_coupling_coef(a_cpl, hvel, do_i, h_harm, bbl_thick, kv_bbl, z_i, h_ml, &
                            dt, j, G, GV, US, CS, visc, forces, work_on_u=.false., OBC=OBC)
    if ( allocated(hML_v)) then
       do i=is,ie ; if (do_i(i)) then ; hML_v(i,J) = h_ml(i) ; endif ; enddo
    endif
    do_any_shelf = .false.
    if (associated(forces%frac_shelf_v)) then
      do i=is,ie
        CS%a1_shelf_v(i,J) = 0.0
        do_i_shelf(i) = (do_i(i) .and. forces%frac_shelf_v(i,J) > 0.0)
        if (do_i_shelf(I)) do_any_shelf = .true.
      enddo
      if (do_any_shelf) then
        if (CS%harmonic_visc) then
          call find_coupling_coef(a_shelf, hvel, do_i_shelf, h_harm, bbl_thick, &
                                  kv_bbl, z_i, h_ml, dt, j, G, GV, US, CS, visc, &
                                  forces, work_on_u=.false., OBC=OBC, shelf=.true.)
        else  ! Find upwind-biased thickness near the surface.
          ! Perhaps this needs to be done more carefully, via find_eta.
          do i=is,ie ; if (do_i_shelf(i)) then
            zh(i) = 0.0 ; Ztop_min(I) = min(zcol1(i), zcol2(i))
            I_HTbl(i) = 1.0 / (visc%tbl_thick_shelf_v(i,J)*GV%Z_to_H + h_neglect)
          endif ; enddo
          do k=1,nz
            do i=is,ie ; if (do_i_shelf(i)) then
              zcol1(i) = zcol1(i) - h(i,j,k) ; zcol2(i) = zcol2(i) - h(i,j+1,k)
              zh(i) = zh(i) + h_harm(i,k)

              hvel_shelf(i,k) = hvel(i,k)
              if (v(i,J,k) * h_delta(i,k) > 0) then
                if (zh(i) * I_HTbl(i) < CS%harm_BL_val) then
                  hvel_shelf(i,k) = min(hvel(i,k), h_harm(i,k))
                else
                  z2_wt = 1.0  ; if (zh(i) * I_HTbl(i) < 2.0*CS%harm_BL_val) &
                    z2_wt = max(0.0, min(1.0, zh(i) * I_HTbl(i) * I_valBL - 1.0))
                  z2 = z2_wt * (max(zh(i), Ztop_min(i) - min(zcol1(i),zcol2(i))) * I_HTbl(i))
                  topfn = 1.0 / (1.0 + 0.09*z2**6)
                  hvel_shelf(i,k) = min(hvel(i,k), (1.0-topfn)*h_arith(i,k) + topfn*h_harm(i,k))
                endif
             endif
            endif ; enddo
          enddo
          call find_coupling_coef(a_shelf, hvel_shelf, do_i_shelf, h_harm, &
                                  bbl_thick, kv_bbl, z_i, h_ml, dt, j, G, GV, US, CS, &
                                  visc, forces, work_on_u=.false., OBC=OBC, shelf=.true.)
        endif
        do i=is,ie ; if (do_i_shelf(i)) CS%a1_shelf_v(i,J) = a_shelf(i,1) ; enddo
      endif
    endif

    if (do_any_shelf) then
      do K=1,nz+1 ; do i=is,ie ; if (do_i_shelf(i)) then
        CS%a_v(i,J,K) = min(a_cpl_max, forces%frac_shelf_v(i,J)  * a_shelf(i,k) + &
                                       (1.0-forces%frac_shelf_v(i,J)) * a_cpl(i,K))
! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
!        CS%a_v(i,J,K) = min(a_cpl_max, forces%frac_shelf_v(i,J)  * max(a_shelf(i,K), a_cpl(i,K)) + &
                    !                   (1.0-forces%frac_shelf_v(i,J)) * a_cpl(i,K))
      elseif (do_i(i)) then
        CS%a_v(i,J,K) = min(a_cpl_max, a_cpl(i,K))
      endif ; enddo ; enddo
      do k=1,nz ; do i=is,ie ; if (do_i_shelf(i)) then
        ! Should we instead take the inverse of the average of the inverses?
        CS%h_v(i,J,k) = forces%frac_shelf_v(i,J)  * hvel_shelf(i,k) + &
                   (1.0-forces%frac_shelf_v(i,J)) * hvel(i,k) + h_neglect
      elseif (do_i(i)) then
        CS%h_v(i,J,k) = hvel(i,k) + h_neglect
      endif ; enddo ; enddo
    else
      do K=1,nz+1 ; do i=is,ie ; if (do_i(i)) CS%a_v(i,J,K) = min(a_cpl_max, a_cpl(i,K)) ; enddo ; enddo
      do k=1,nz ; do i=is,ie ; if (do_i(i)) CS%h_v(i,J,k) = hvel(i,k) + h_neglect ; enddo ; enddo
    endif

    ! Diagnose total Kv at v-points
    if (CS%id_Kv_v > 0) then
      do k=1,nz ; do i=is,ie
        if (do_i(I)) Kv_v(i,J,k) = 0.5 * GV%H_to_Z*(CS%a_v(i,J,K)+CS%a_v(i,J,K+1)) * CS%h_v(i,J,k)
      enddo ; enddo
    endif

  enddo ! end of v-point j loop

  if (CS%debug) then
    call uvchksum("vertvisc_coef h_[uv]", CS%h_u, CS%h_v, G%HI, haloshift=0, &
                  scale=GV%H_to_m, scalar_pair=.true.)
    call uvchksum("vertvisc_coef a_[uv]", CS%a_u, CS%a_v, G%HI, haloshift=0, &
                  scale=US%Z_to_m*US%s_to_T, scalar_pair=.true.)
    if (allocated(hML_u) .and. allocated(hML_v)) &
      call uvchksum("vertvisc_coef hML_[uv]", hML_u, hML_v, G%HI, &
                    haloshift=0, scale=GV%H_to_m, scalar_pair=.true.)
  endif

! Offer diagnostic fields for averaging.
  if (associated(visc%Kv_slow) .and. (CS%id_Kv_slow > 0)) &
      call post_data(CS%id_Kv_slow, visc%Kv_slow, CS%diag)
  if (CS%id_Kv_u > 0) call post_data(CS%id_Kv_u, Kv_u, CS%diag)
  if (CS%id_Kv_v > 0) call post_data(CS%id_Kv_v, Kv_v, CS%diag)
  if (CS%id_au_vv > 0) call post_data(CS%id_au_vv, CS%a_u, CS%diag)
  if (CS%id_av_vv > 0) call post_data(CS%id_av_vv, CS%a_v, CS%diag)
  if (CS%id_h_u > 0) call post_data(CS%id_h_u, CS%h_u, CS%diag)
  if (CS%id_h_v > 0) call post_data(CS%id_h_v, CS%h_v, CS%diag)
  if (CS%id_hML_u > 0) call post_data(CS%id_hML_u, hML_u, CS%diag)
  if (CS%id_hML_v > 0) call post_data(CS%id_hML_v, hML_v, CS%diag)

  if (allocated(hML_u)) deallocate(hML_u)
  if (allocated(hML_v)) deallocate(hML_v)

end subroutine vertvisc_coef

!> Calculate the 'coupling coefficient' (a_cpl) at the interfaces.
!! If BOTTOMDRAGLAW is defined, the minimum of Hbbl and half the adjacent
!! layer thicknesses are used to calculate a_cpl near the bottom.
subroutine find_coupling_coef(a_cpl, hvel, do_i, h_harm, bbl_thick, kv_bbl, z_i, h_ml, &
                              dt, j, G, GV, US, CS, visc, forces, work_on_u, OBC, shelf)
  type(ocean_grid_type),     intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),   intent(in)  :: GV !< Ocean vertical grid structure
  type(unit_scale_type),     intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZK_(GV)+1), &
                             intent(out) :: a_cpl !< Coupling coefficient across interfaces [Z T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZK_(GV)), &
                             intent(in)  :: hvel !< Thickness at velocity points [H ~> m or kg m-2]
  logical, dimension(SZIB_(G)), &
                             intent(in)  :: do_i !< If true, determine coupling coefficient for a column
  real, dimension(SZIB_(G),SZK_(GV)), &
                             intent(in)  :: h_harm !< Harmonic mean of thicknesses around a velocity
                                                   !! grid point [H ~> m or kg m-2]
  real, dimension(SZIB_(G)), intent(in)  :: bbl_thick !< Bottom boundary layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G)), intent(in)  :: kv_bbl !< Bottom boundary layer viscosity [Z2 T-1 ~> m2 s-1].
  real, dimension(SZIB_(G),SZK_(GV)+1), &
                             intent(in)  :: z_i  !< Estimate of interface heights above the bottom,
                                                 !! normalized by the bottom boundary layer thickness
  real, dimension(SZIB_(G)), intent(out) :: h_ml !< Mixed layer depth [H ~> m or kg m-2]
  integer,                   intent(in)  :: j    !< j-index to find coupling coefficient for
  real,                      intent(in)  :: dt   !< Time increment [T ~> s]
  type(vertvisc_CS),         pointer     :: CS   !< Vertical viscosity control structure
  type(vertvisc_type),       intent(in)  :: visc !< Structure containing viscosities and bottom drag
  type(mech_forcing),        intent(in)  :: forces !< A structure with the driving mechanical forces
  logical,                   intent(in)  :: work_on_u !< If true, u-points are being calculated,
                                                  !! otherwise they are v-points
  type(ocean_OBC_type),      pointer     :: OBC   !< Open boundary condition structure
  logical,         optional, intent(in)  :: shelf !< If present and true, use a surface boundary
                                                  !! condition appropriate for an ice shelf.

  ! Local variables

  real, dimension(SZIB_(G)) :: &
    u_star, &   ! ustar at a velocity point [Z T-1 ~> m s-1].
    absf, &     ! The average of the neighboring absolute values of f [T-1 ~> s-1].
!      h_ml, &  ! The mixed layer depth [H ~> m or kg m-2].
    nk_visc, &  ! The (real) interface index of the base of mixed layer.
    z_t, &      ! The distance from the top, sometimes normalized
                ! by Hmix, [H ~> m or kg m-2] or [nondim].
    kv_TBL, &   ! The viscosity in a top boundary layer under ice [Z2 T-1 ~> m2 s-1].
    tbl_thick
  real, dimension(SZIB_(G),SZK_(GV)+1) :: &
    Kv_tot, &   ! The total viscosity at an interface [Z2 T-1 ~> m2 s-1].
    Kv_add      ! A viscosity to add [Z2 T-1 ~> m2 s-1].
  real :: h_shear ! The distance over which shears occur [H ~> m or kg m-2].
  real :: r       ! A thickness to compare with Hbbl [H ~> m or kg m-2].
  real :: visc_ml ! The mixed layer viscosity [Z2 T-1 ~> m2 s-1].
  real :: I_Hmix  ! The inverse of the mixed layer thickness [H-1 ~> m-1 or m2 kg-1].
  real :: a_ml    ! The layer coupling coefficient across an interface in
                  ! the mixed layer [Z T-1 ~> m s-1].
  real :: I_amax  ! The inverse of the maximum coupling coefficient [T Z-1 ~> s m-1].
  real :: temp1   ! A temporary variable [H Z ~> m2 or kg m-1]
  real :: h_neglect   ! A thickness that is so small it is usually lost
                      ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: z2      ! A copy of z_i [nondim]
  real :: botfn   ! A function that is 1 at the bottom and small far from it [nondim]
  real :: topfn   ! A function that is 1 at the top and small far from it [nondim]
  real :: kv_top  ! A viscosity associated with the top boundary layer [Z2 T-1 ~> m2 s-1]
  logical :: do_shelf, do_OBCs
  integer :: i, k, is, ie, max_nk
  integer :: nz

  a_cpl(:,:) = 0.0
  Kv_tot(:,:) = 0.0

  if (work_on_u) then ; is = G%IscB ; ie = G%IecB
  else ; is = G%isc ; ie = G%iec ; endif
  nz = G%ke
  h_neglect = GV%H_subroundoff

  if (CS%answers_2018) then
    !   The maximum coupling coefficent was originally introduced to avoid
    ! truncation error problems in the tridiagonal solver. Effectively, the 1e-10
    ! sets the maximum coupling coefficient increment to 1e10 m per timestep.
    I_amax = (1.0e-10*US%Z_to_m) * dt
  else
    I_amax = 0.0
  endif

  do_shelf = .false. ; if (present(shelf)) do_shelf = shelf
  do_OBCs = .false.
  if (associated(OBC)) then ; do_OBCS = (OBC%number_of_segments > 0) ; endif
  h_ml(:) = 0.0

!    The following loop calculates the vertical average velocity and
!  surface mixed layer contributions to the vertical viscosity.
  do i=is,ie ; Kv_tot(i,1) = 0.0 ; enddo
  if ((GV%nkml>0) .or. do_shelf) then ; do k=2,nz ; do i=is,ie
    if (do_i(i)) Kv_tot(i,K) = CS%Kv
  enddo ; enddo ; else
    I_Hmix = 1.0 / (CS%Hmix + h_neglect)
    do i=is,ie ; z_t(i) = h_neglect*I_Hmix ; enddo
    do K=2,nz ; do i=is,ie ; if (do_i(i)) then
      z_t(i) = z_t(i) + h_harm(i,k-1)*I_Hmix
      Kv_tot(i,K) = CS%Kv + CS%Kvml / ((z_t(i)*z_t(i)) *  &
               (1.0 + 0.09*z_t(i)*z_t(i)*z_t(i)*z_t(i)*z_t(i)*z_t(i)))
    endif ; enddo ; enddo
  endif

  do i=is,ie ; if (do_i(i)) then
    if (CS%bottomdraglaw) then
      r = hvel(i,nz)*0.5
      if (r < bbl_thick(i)) then
        a_cpl(i,nz+1) = kv_bbl(i) / (I_amax*kv_bbl(i) + (r+h_neglect)*GV%H_to_Z)
      else
        a_cpl(i,nz+1) = kv_bbl(i) / (I_amax*kv_bbl(i) + (bbl_thick(i)+h_neglect)*GV%H_to_Z)
      endif
    else
      a_cpl(i,nz+1) = CS%Kvbbl / ((0.5*hvel(i,nz)+h_neglect)*GV%H_to_Z + I_amax*CS%Kvbbl)
    endif
  endif ; enddo

  if (associated(visc%Kv_shear)) then
    ! The factor of 2 that used to be required in the viscosities is no longer needed.
    if (work_on_u) then
      do K=2,nz ; do i=is,ie ; if (do_i(i)) then
        Kv_add(i,K) = 0.5*(visc%Kv_shear(i,j,k) + visc%Kv_shear(i+1,j,k))
      endif ; enddo ; enddo
      if (do_OBCs) then
        do I=is,ie ; if (do_i(I) .and. (OBC%segnum_u(I,j) /= OBC_NONE)) then
          if (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_E) then
            do K=2,nz ; Kv_add(i,K) = visc%Kv_shear(i,j,k) ; enddo
          elseif (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_W) then
            do K=2,nz ; Kv_add(i,K) = visc%Kv_shear(i+1,j,k) ; enddo
          endif
        endif ; enddo
      endif
      do K=2,nz ; do i=is,ie ; if (do_i(i)) then
        Kv_tot(i,K) = Kv_tot(i,K) + Kv_add(i,K)
      endif ; enddo ; enddo
    else
      do K=2,nz ; do i=is,ie ; if (do_i(i)) then
        Kv_add(i,K) = 0.5*(visc%Kv_shear(i,j,k) + visc%Kv_shear(i,j+1,k))
      endif ; enddo ; enddo
      if (do_OBCs) then
        do i=is,ie ; if (do_i(i) .and. (OBC%segnum_v(i,J) /= OBC_NONE)) then
          if (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_N) then
            do K=2,nz ; Kv_add(i,K) = visc%Kv_shear(i,j,k) ; enddo
          elseif (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_S) then
            do K=2,nz ; Kv_add(i,K) = visc%Kv_shear(i,j+1,k) ; enddo
          endif
        endif ; enddo
      endif
      do K=2,nz ; do i=is,ie ; if (do_i(i)) then
        Kv_tot(i,K) = Kv_tot(i,K) + Kv_add(i,K)
      endif ; enddo ; enddo
    endif
  endif

  if (associated(visc%Kv_shear_Bu)) then
    if (work_on_u) then
      do K=2,nz ; do I=Is,Ie ; If (do_i(I)) then
        Kv_tot(I,K) = Kv_tot(I,K) + (0.5)*(visc%Kv_shear_Bu(I,J-1,k) + visc%Kv_shear_Bu(I,J,k))
      endif ; enddo ; enddo
    else
      do K=2,nz ; do i=is,ie ; if (do_i(i)) then
        Kv_tot(i,K) = Kv_tot(i,K) + (0.5)*(visc%Kv_shear_Bu(I-1,J,k) + visc%Kv_shear_Bu(I,J,k))
      endif ; enddo ; enddo
    endif
  endif

  do K=nz,2,-1 ; do i=is,ie ; if (do_i(i)) then
    !    botfn determines when a point is within the influence of the bottom
    !  boundary layer, going from 1 at the bottom to 0 in the interior.
    z2 = z_i(i,k)
    botfn = 1.0 / (1.0 + 0.09*z2*z2*z2*z2*z2*z2)

    if (CS%bottomdraglaw) then
      Kv_tot(i,K) = Kv_tot(i,K) + (kv_bbl(i) - CS%Kv)*botfn
      r = 0.5*(hvel(i,k) + hvel(i,k-1))
      if (r > bbl_thick(i)) then
        h_shear = ((1.0 - botfn) * r + botfn*bbl_thick(i)) + h_neglect
      else
        h_shear = r + h_neglect
      endif
    else
      Kv_tot(i,K) = Kv_tot(i,K) + (CS%Kvbbl-CS%Kv)*botfn
      h_shear = 0.5*(hvel(i,k) + hvel(i,k-1) + h_neglect)
    endif

    ! Calculate the coupling coefficients from the viscosities.
    a_cpl(i,K) = Kv_tot(i,K) / (h_shear*GV%H_to_Z + I_amax*Kv_tot(i,K))
  endif ; enddo ; enddo ! i & k loops

  if (do_shelf) then
    ! Set the coefficients to include the no-slip surface stress.
    do i=is,ie ; if (do_i(i)) then
      if (work_on_u) then
        kv_TBL(i) = visc%Kv_tbl_shelf_u(I,j)
        tbl_thick(i) = visc%tbl_thick_shelf_u(I,j) * GV%Z_to_H + h_neglect
      else
        kv_TBL(i) = visc%Kv_tbl_shelf_v(i,J)
        tbl_thick(i) = visc%tbl_thick_shelf_v(i,J) * GV%Z_to_H + h_neglect
      endif
      z_t(i) = 0.0

      ! If a_cpl(i,1) were not already 0, it would be added here.
      if (0.5*hvel(i,1) > tbl_thick(i)) then
        a_cpl(i,1) = kv_TBL(i) / (tbl_thick(i)*GV%H_to_Z + I_amax*kv_TBL(i))
      else
        a_cpl(i,1) = kv_TBL(i) / ((0.5*hvel(i,1)+h_neglect)*GV%H_to_Z + I_amax*kv_TBL(i))
      endif
    endif ; enddo

    do K=2,nz ; do i=is,ie ;  if (do_i(i)) then
      z_t(i) = z_t(i) + hvel(i,k-1) / tbl_thick(i)
      topfn = 1.0 / (1.0 + 0.09 * z_t(i)**6)

      r = 0.5*(hvel(i,k)+hvel(i,k-1))
      if (r > tbl_thick(i)) then
        h_shear = ((1.0 - topfn) * r + topfn*tbl_thick(i)) + h_neglect
      else
        h_shear = r + h_neglect
      endif

      kv_top = topfn * kv_TBL(i)
      a_cpl(i,K) = a_cpl(i,K) + kv_top / (h_shear*GV%H_to_Z + I_amax*kv_top)
    endif ; enddo ; enddo
  elseif (CS%dynamic_viscous_ML .or. (GV%nkml>0)) then
    max_nk = 0
    do i=is,ie ; if (do_i(i)) then
      if (GV%nkml>0) nk_visc(i) = real(GV%nkml+1)
      if (work_on_u) then
        u_star(I) = 0.5*(forces%ustar(i,j) + forces%ustar(i+1,j))
        absf(I) = 0.5*(abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I,J)))
        if (CS%dynamic_viscous_ML) nk_visc(I) = visc%nkml_visc_u(I,j) + 1
      else
        u_star(i) = 0.5*(forces%ustar(i,j) + forces%ustar(i,j+1))
        absf(i) = 0.5*(abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))
        if (CS%dynamic_viscous_ML) nk_visc(i) = visc%nkml_visc_v(i,J) + 1
      endif
      h_ml(i) = h_neglect ; z_t(i) = 0.0
      max_nk = max(max_nk,ceiling(nk_visc(i) - 1.0))
    endif ; enddo

    if (do_OBCS) then ; if (work_on_u) then
      do I=is,ie ; if (do_i(I) .and. (OBC%segnum_u(I,j) /= OBC_NONE)) then
        if (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_E) &
          u_star(I) = forces%ustar(i,j)
        if (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_W) &
          u_star(I) = forces%ustar(i+1,j)
      endif ; enddo
    else
      do i=is,ie ; if (do_i(i) .and. (OBC%segnum_v(i,J) /= OBC_NONE)) then
        if (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_N) &
          u_star(i) = forces%ustar(i,j)
        if (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_S) &
          u_star(i) = forces%ustar(i,j+1)
      endif ; enddo
    endif ; endif

    do k=1,max_nk ; do i=is,ie ; if (do_i(i)) then
      if (k+1 <= nk_visc(i)) then ! This layer is all in the ML.
        h_ml(i) = h_ml(i) + hvel(i,k)
      elseif (k < nk_visc(i)) then ! Part of this layer is in the ML.
        h_ml(i) = h_ml(i) + (nk_visc(i) - k) * hvel(i,k)
      endif
    endif ; enddo ; enddo

    do K=2,max_nk ; do i=is,ie ; if (do_i(i)) then ; if (k < nk_visc(i)) then
      ! Set the viscosity at the interfaces.
      z_t(i) = z_t(i) + hvel(i,k-1)
      temp1 = (z_t(i)*h_ml(i) - z_t(i)*z_t(i))*GV%H_to_Z
      !   This viscosity is set to go to 0 at the mixed layer top and bottom (in a log-layer)
      ! and be further limited by rotation to give the natural Ekman length.
      visc_ml = u_star(i) * 0.41 * (temp1*u_star(i)) / (absf(i)*temp1 + (h_ml(i)+h_neglect)*u_star(i))
      a_ml = visc_ml / (0.25*(hvel(i,k)+hvel(i,k-1) + h_neglect) * GV%H_to_Z + 0.5*I_amax*visc_ml)
      ! Choose the largest estimate of a.
      if (a_ml > a_cpl(i,K)) a_cpl(i,K) = a_ml
    endif ; endif ; enddo ; enddo
  endif

end subroutine find_coupling_coef

!> Velocity components which exceed a threshold for physically reasonable values
!! are truncated. Optionally, any column with excessive velocities may be sent
!! to a diagnostic reporting subroutine.
subroutine vertvisc_limit_vel(u, v, h, ADp, CDp, forces, visc, dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  type(accel_diag_ptrs),   intent(in)    :: ADp    !< Acceleration diagnostic pointers
  type(cont_diag_ptrs),    intent(in)    :: CDp    !< Continuity diagnostic pointers
  type(mech_forcing),      intent(in)    :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),     intent(in)    :: visc   !< Viscosities and bottom drag
  real,                    intent(in)    :: dt     !< Time increment [T ~> s]
  type(vertvisc_CS),       pointer       :: CS     !< Vertical viscosity control structure

  ! Local variables

  real :: maxvel           ! Velocities components greater than maxvel
  real :: truncvel         ! are truncated to truncvel, both [L T-1 ~> m s-1].
  real :: CFL              ! The local CFL number.
  real :: H_report         ! A thickness below which not to report truncations.
  real :: dt_Rho0          ! The timestep divided by the Boussinesq density [m2 T2 s-1 L-1 Z-1 R-1 ~> s m3 kg-1].
  real :: vel_report(SZIB_(G),SZJB_(G))   ! The velocity to report [L T-1 ~> m s-1]
  real :: u_old(SZIB_(G),SZJ_(G),SZK_(G)) ! The previous u-velocity [L T-1 ~> m s-1]
  real :: v_old(SZI_(G),SZJB_(G),SZK_(G)) ! The previous v-velocity [L T-1 ~> m s-1]
  logical :: trunc_any, dowrite(SZIB_(G),SZJB_(G))
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  maxvel = CS%maxvel
  truncvel = 0.9*maxvel
  H_report = 6.0 * GV%Angstrom_H
  dt_Rho0 = (US%L_T_to_m_s*US%Z_to_m) * dt / GV%Rho0

  if (len_trim(CS%u_trunc_file) > 0) then
    !$OMP parallel do default(shared) private(trunc_any,CFL)
    do j=js,je
      trunc_any = .false.
      do I=Isq,Ieq ; dowrite(I,j) = .false. ; enddo
      if (CS%CFL_based_trunc) then
        do I=Isq,Ieq ; vel_report(i,j) = 3.0e8*US%m_s_to_L_T ; enddo ! Speed of light default.
        do k=1,nz ; do I=Isq,Ieq
          if (abs(u(I,j,k)) < CS%vel_underflow) u(I,j,k) = 0.0
          if (u(I,j,k) < 0.0) then
            CFL = (-u(I,j,k) * dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
          else
            CFL = (u(I,j,k) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
          endif
          if (CFL > CS%CFL_trunc) trunc_any = .true.
          if (CFL > CS%CFL_report) then
            dowrite(I,j) = .true.
            vel_report(I,j) = MIN(vel_report(I,j), abs(u(I,j,k)))
          endif
        enddo ; enddo
      else
        do I=Isq,Ieq; vel_report(I,j) = maxvel; enddo
        do k=1,nz ; do I=Isq,Ieq
          if (abs(u(I,j,k)) < CS%vel_underflow) then ; u(I,j,k) = 0.0
          elseif (abs(u(I,j,k)) > maxvel) then
            dowrite(I,j) = .true. ; trunc_any = .true.
          endif
        enddo ; enddo
      endif

      do I=Isq,Ieq ; if (dowrite(I,j)) then
        u_old(I,j,:) = u(I,j,:)
      endif ; enddo

      if (trunc_any) then ; if (CS%CFL_based_trunc) then
        do k=1,nz ; do I=Isq,Ieq
          if ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i+1,j) < -CS%CFL_trunc) then
            u(I,j,k) = (-0.9*CS%CFL_trunc) * (G%areaT(i+1,j) / (dt * G%dy_Cu(I,j)))
            if (h(i,j,k) + h(i+1,j,k) > H_report) CS%ntrunc = CS%ntrunc + 1
          elseif ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i,j) > CS%CFL_trunc) then
            u(I,j,k) = (0.9*CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dy_Cu(I,j)))
            if (h(i,j,k) + h(i+1,j,k) > H_report) CS%ntrunc = CS%ntrunc + 1
          endif
        enddo ; enddo
      else
        do k=1,nz ; do I=Isq,Ieq ; if (abs(u(I,j,k)) > maxvel) then
          u(I,j,k) = SIGN(truncvel,u(I,j,k))
          if (h(i,j,k) + h(i+1,j,k) > H_report) CS%ntrunc = CS%ntrunc + 1
        endif ; enddo ; enddo
      endif ; endif
    enddo ! j-loop
  else  ! Do not report accelerations leading to large velocities.
    if (CS%CFL_based_trunc) then
!$OMP parallel do default(none) shared(nz,js,je,Isq,Ieq,u,dt,G,CS,h,H_report)
      do k=1,nz ; do j=js,je ; do I=Isq,Ieq
        if (abs(u(I,j,k)) < CS%vel_underflow) then ; u(I,j,k) = 0.0
        elseif ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i+1,j) < -CS%CFL_trunc) then
          u(I,j,k) = (-0.9*CS%CFL_trunc) * (G%areaT(i+1,j) / (dt * G%dy_Cu(I,j)))
          if (h(i,j,k) + h(i+1,j,k) > H_report) CS%ntrunc = CS%ntrunc + 1
        elseif ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i,j) > CS%CFL_trunc) then
          u(I,j,k) = (0.9*CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dy_Cu(I,j)))
          if (h(i,j,k) + h(i+1,j,k) > H_report) CS%ntrunc = CS%ntrunc + 1
        endif
      enddo ; enddo ; enddo
    else
!$OMP parallel do default(none) shared(nz,js,je,Isq,Ieq,u,G,CS,truncvel,maxvel,h,H_report)
      do k=1,nz ; do j=js,je ; do I=Isq,Ieq
        if (abs(u(I,j,k)) < CS%vel_underflow) then ; u(I,j,k) = 0.0
        elseif (abs(u(I,j,k)) > maxvel) then
          u(I,j,k) = SIGN(truncvel, u(I,j,k))
          if (h(i,j,k) + h(i+1,j,k) > H_report) CS%ntrunc = CS%ntrunc + 1
        endif
      enddo ; enddo ; enddo
    endif
  endif

  if (len_trim(CS%u_trunc_file) > 0) then
    do j=js,je; do I=Isq,Ieq ; if (dowrite(I,j)) then
!   Here the diagnostic reporting subroutines are called if
! unphysically large values were found.
      call write_u_accel(I, j, u_old, h, ADp, CDp, dt, G, GV, US, CS%PointAccel_CSp, &
               vel_report(I,j), forces%taux(I,j)*dt_Rho0, a=CS%a_u, hv=CS%h_u)
    endif ; enddo ; enddo
  endif

  if (len_trim(CS%v_trunc_file) > 0) then
    !$OMP parallel do default(shared) private(trunc_any,CFL)
    do J=Jsq,Jeq
      trunc_any = .false.
      do i=is,ie ; dowrite(i,J) = .false. ; enddo
      if (CS%CFL_based_trunc) then
        do i=is,ie ; vel_report(i,J) = 3.0e8*US%m_s_to_L_T ; enddo ! Speed of light default.
        do k=1,nz ; do i=is,ie
          if (abs(v(i,J,k)) < CS%vel_underflow) v(i,J,k) = 0.0
          if (v(i,J,k) < 0.0) then
            CFL = (-v(i,J,k) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
          else
            CFL = (v(i,J,k) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
          endif
          if (CFL > CS%CFL_trunc) trunc_any = .true.
          if (CFL > CS%CFL_report) then
            dowrite(i,J) = .true.
            vel_report(i,J) = MIN(vel_report(i,J), abs(v(i,J,k)))
          endif
        enddo ; enddo
      else
        do i=is,ie ; vel_report(i,J) = maxvel ; enddo
        do k=1,nz ; do i=is,ie
          if (abs(v(i,J,k)) < CS%vel_underflow) then ; v(i,J,k) = 0.0
          elseif (abs(v(i,J,k)) > maxvel) then
            dowrite(i,J) = .true. ; trunc_any = .true.
          endif
        enddo ; enddo
      endif

      do i=is,ie ; if (dowrite(i,J)) then
        v_old(i,J,:) = v(i,J,:)
      endif ; enddo

      if (trunc_any) then ; if (CS%CFL_based_trunc) then
        do k=1,nz; do i=is,ie
          if ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j+1) < -CS%CFL_trunc) then
            v(i,J,k) = (-0.9*CS%CFL_trunc) * (G%areaT(i,j+1) / (dt * G%dx_Cv(i,J)))
            if (h(i,j,k) + h(i,j+1,k) > H_report) CS%ntrunc = CS%ntrunc + 1
          elseif ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j) > CS%CFL_trunc) then
            v(i,J,k) = (0.9*CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dx_Cv(i,J)))
            if (h(i,j,k) + h(i,j+1,k) > H_report) CS%ntrunc = CS%ntrunc + 1
          endif
        enddo ; enddo
      else
        do k=1,nz ; do i=is,ie ; if (abs(v(i,J,k)) > maxvel) then
          v(i,J,k) = SIGN(truncvel,v(i,J,k))
          if (h(i,j,k) + h(i,j+1,k) > H_report) CS%ntrunc = CS%ntrunc + 1
        endif ; enddo ; enddo
      endif ; endif
    enddo ! J-loop
  else  ! Do not report accelerations leading to large velocities.
    if (CS%CFL_based_trunc) then
      !$OMP parallel do default(shared)
      do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
        if (abs(v(i,J,k)) < CS%vel_underflow) then ; v(i,J,k) = 0.0
        elseif ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j+1) < -CS%CFL_trunc) then
          v(i,J,k) = (-0.9*CS%CFL_trunc) * (G%areaT(i,j+1) / (dt * G%dx_Cv(i,J)))
          if (h(i,j,k) + h(i,j+1,k) > H_report) CS%ntrunc = CS%ntrunc + 1
        elseif ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j) > CS%CFL_trunc) then
          v(i,J,k) = (0.9*CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dx_Cv(i,J)))
          if (h(i,j,k) + h(i,j+1,k) > H_report) CS%ntrunc = CS%ntrunc + 1
        endif
      enddo ; enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
        if (abs(v(i,J,k)) < CS%vel_underflow) then ; v(i,J,k) = 0.0
        elseif (abs(v(i,J,k)) > maxvel) then
          v(i,J,k) = SIGN(truncvel, v(i,J,k))
          if (h(i,j,k) + h(i,j+1,k) > H_report) CS%ntrunc = CS%ntrunc + 1
        endif
      enddo ; enddo ; enddo
    endif
  endif

  if (len_trim(CS%v_trunc_file) > 0) then
    do J=Jsq,Jeq; do i=is,ie ; if (dowrite(i,J)) then
!   Here the diagnostic reporting subroutines are called if
! unphysically large values were found.
      call write_v_accel(i, J, v_old, h, ADp, CDp, dt, G, GV, US, CS%PointAccel_CSp, &
               vel_report(i,J), forces%tauy(i,J)*dt_Rho0, a=CS%a_v, hv=CS%h_v)
    endif ; enddo ; enddo
  endif

end subroutine vertvisc_limit_vel

!> Initialize the vertical friction module
subroutine vertvisc_init(MIS, Time, G, GV, US, param_file, diag, ADp, dirs, &
                         ntrunc, CS)
  type(ocean_internal_state), &
                   target, intent(in)    :: MIS    !< The "MOM Internal State", a set of pointers
                                                   !! to the fields and accelerations that make
                                                   !! up the ocean's physical state
  type(time_type), target, intent(in)    :: Time   !< Current model time
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< File to parse for parameters
  type(diag_ctrl), target, intent(inout) :: diag   !< Diagnostic control structure
  type(accel_diag_ptrs),   intent(inout) :: ADp    !< Acceleration diagnostic pointers
  type(directories),       intent(in)    :: dirs   !< Relevant directory paths
  integer, target,         intent(inout) :: ntrunc !< Number of velocity truncations
  type(vertvisc_CS),       pointer       :: CS     !< Vertical viscosity control structure

  ! Local variables

  real :: hmix_str_dflt
  real :: Kv_dflt ! A default viscosity [m2 s-1].
  real :: Hmix_m  ! A boundary layer thickness [m].
  logical :: default_2018_answers
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_vert_friction" ! This module's name.
  character(len=40)  :: thickness_units

  if (associated(CS)) then
    call MOM_error(WARNING, "vertvisc_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  if (GV%Boussinesq) then; thickness_units = "m"
  else; thickness_units = "kg m-2"; endif

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%diag => diag ; CS%ntrunc => ntrunc ; ntrunc = 0

! Default, read and log parameters
  call log_version(param_file, mdl, version, "", log_to_all=.true., debugging=.true.)
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  call get_param(param_file, mdl, "VERT_FRICTION_2018_ANSWERS", CS%answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the answers "//&
                 "from the end of 2018.  Otherwise, use expressions that do not use an arbitrary "//&
                 "hard-coded maximum viscous coupling coefficient between layers.", &
                 default=default_2018_answers)
  call get_param(param_file, mdl, "BOTTOMDRAGLAW", CS%bottomdraglaw, &
                 "If true, the bottom stress is calculated with a drag "//&
                 "law of the form c_drag*|u|*u. The velocity magnitude "//&
                 "may be an assumed value or it may be based on the "//&
                 "actual velocity in the bottommost HBBL, depending on "//&
                 "LINEAR_DRAG.", default=.true.)
  call get_param(param_file, mdl, "CHANNEL_DRAG", CS%Channel_drag, &
                 "If true, the bottom drag is exerted directly on each "//&
                 "layer proportional to the fraction of the bottom it "//&
                 "overlies.", default=.false.)
  call get_param(param_file, mdl, "DIRECT_STRESS", CS%direct_stress, &
                 "If true, the wind stress is distributed over the "//&
                 "topmost HMIX_STRESS of fluid (like in HYCOM), and KVML "//&
                 "may be set to a very small value.", default=.false.)
  call get_param(param_file, mdl, "DYNAMIC_VISCOUS_ML", CS%dynamic_viscous_ML, &
                 "If true, use a bulk Richardson number criterion to "//&
                 "determine the mixed layer thickness for viscosity.", &
                 default=.false.)
  call get_param(param_file, mdl, "U_TRUNC_FILE", CS%u_trunc_file, &
                 "The absolute path to a file into which the accelerations "//&
                 "leading to zonal velocity truncations are written. "//&
                 "Undefine this for efficiency if this diagnostic is not "//&
                 "needed.", default=" ", debuggingParam=.true.)
  call get_param(param_file, mdl, "V_TRUNC_FILE", CS%v_trunc_file, &
                 "The absolute path to a file into which the accelerations "//&
                 "leading to meridional velocity truncations are written. "//&
                 "Undefine this for efficiency if this diagnostic is not "//&
                 "needed.", default=" ", debuggingParam=.true.)
  call get_param(param_file, mdl, "HARMONIC_VISC", CS%harmonic_visc, &
                 "If true, use the harmonic mean thicknesses for "//&
                 "calculating the vertical viscosity.", default=.false.)
  call get_param(param_file, mdl, "HARMONIC_BL_SCALE", CS%harm_BL_val, &
                 "A scale to determine when water is in the boundary "//&
                 "layers based solely on harmonic mean thicknesses for "//&
                 "the purpose of determining the extent to which the "//&
                 "thicknesses used in the viscosities are upwinded.", &
                 default=0.0, units="nondim")
  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false.)

  if (GV%nkml < 1) &
    call get_param(param_file, mdl, "HMIX_FIXED", CS%Hmix, &
                 "The prescribed depth over which the near-surface "//&
                 "viscosity and diffusivity are elevated when the bulk "//&
                 "mixed layer is not used.", units="m", scale=GV%m_to_H, &
                 unscaled=Hmix_m, fail_if_missing=.true.)
  if (CS%direct_stress) then
    if (GV%nkml < 1) then
      call get_param(param_file, mdl, "HMIX_STRESS", CS%Hmix_stress, &
                 "The depth over which the wind stress is applied if "//&
                 "DIRECT_STRESS is true.", units="m", default=Hmix_m, scale=GV%m_to_H)
    else
      call get_param(param_file, mdl, "HMIX_STRESS", CS%Hmix_stress, &
                 "The depth over which the wind stress is applied if "//&
                 "DIRECT_STRESS is true.", units="m", fail_if_missing=.true., scale=GV%m_to_H)
    endif
    if (CS%Hmix_stress <= 0.0) call MOM_error(FATAL, "vertvisc_init: " // &
       "HMIX_STRESS must be set to a positive value if DIRECT_STRESS is true.")
  endif
  call get_param(param_file, mdl, "KV", CS%Kv, &
                 "The background kinematic viscosity in the interior. "//&
                 "The molecular value, ~1e-6 m2 s-1, may be used.", &
                 units="m2 s-1", fail_if_missing=.true., scale=US%m2_s_to_Z2_T, unscaled=Kv_dflt)

  if (GV%nkml < 1) call get_param(param_file, mdl, "KVML", CS%Kvml, &
                 "The kinematic viscosity in the mixed layer.  A typical "//&
                 "value is ~1e-2 m2 s-1. KVML is not used if "//&
                 "BULKMIXEDLAYER is true.  The default is set by KV.", &
                 units="m2 s-1", default=Kv_dflt, scale=US%m2_s_to_Z2_T)
  if (.not.CS%bottomdraglaw) call get_param(param_file, mdl, "KVBBL", CS%Kvbbl, &
                 "The kinematic viscosity in the benthic boundary layer. "//&
                 "A typical value is ~1e-2 m2 s-1. KVBBL is not used if "//&
                 "BOTTOMDRAGLAW is true.  The default is set by KV.", &
                 units="m2 s-1", default=Kv_dflt, scale=US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "HBBL", CS%Hbbl, &
                 "The thickness of a bottom boundary layer with a "//&
                 "viscosity of KVBBL if BOTTOMDRAGLAW is not defined, or "//&
                 "the thickness over which near-bottom velocities are "//&
                 "averaged for the drag law if BOTTOMDRAGLAW is defined "//&
                 "but LINEAR_DRAG is not.", units="m", fail_if_missing=.true., scale=GV%m_to_H)
  call get_param(param_file, mdl, "MAXVEL", CS%maxvel, &
                 "The maximum velocity allowed before the velocity "//&
                 "components are truncated.", units="m s-1", default=3.0e8, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "CFL_BASED_TRUNCATIONS", CS%CFL_based_trunc, &
                 "If true, base truncations on the CFL number, and not an "//&
                 "absolute speed.", default=.true.)
  call get_param(param_file, mdl, "CFL_TRUNCATE", CS%CFL_trunc, &
                 "The value of the CFL number that will cause velocity "//&
                 "components to be truncated; instability can occur past 0.5.", &
                 units="nondim", default=0.5)
  call get_param(param_file, mdl, "CFL_REPORT", CS%CFL_report, &
                 "The value of the CFL number that causes accelerations "//&
                 "to be reported; the default is CFL_TRUNCATE.", &
                 units="nondim", default=CS%CFL_trunc)
  call get_param(param_file, mdl, "CFL_TRUNCATE_RAMP_TIME", CS%truncRampTime, &
                 "The time over which the CFL truncation value is ramped "//&
                 "up at the beginning of the run.", &
                 units="s", default=0.)
  CS%CFL_truncE = CS%CFL_trunc
  call get_param(param_file, mdl, "CFL_TRUNCATE_START", CS%CFL_truncS, &
                 "The start value of the truncation CFL number used when "//&
                 "ramping up CFL_TRUNC.", &
                 units="nondim", default=0.)
  call get_param(param_file, mdl, "STOKES_MIXING_COMBINED", CS%StokesMixing, &
                 "Flag to use Stokes drift Mixing via the Lagrangian "//&
                 " current (Eulerian plus Stokes drift). "//&
                 " Still needs work and testing, so not recommended for use.",&
                 Default=.false.)
  !BGR 04/04/2018{
  ! StokesMixing is required for MOM6 for some Langmuir mixing parameterization.
  !   The code used here has not been developed for vanishing layers or in
  !   conjunction with any bottom friction.  Therefore, the following line is
  !   added so this functionality cannot be used without user intervention in
  !   the code.  This will prevent general use of this functionality until proper
  !   care is given to the previously mentioned issues.  Comment out the following
  !   MOM_error to use, but do so at your own risk and with these points in mind.
  !}
  if (CS%StokesMixing) then
    call MOM_error(FATAL, "Stokes mixing requires user intervention in the code.\n"//&
                          "  Model now exiting.  See MOM_vert_friction.F90 for \n"//&
                          "  details (search 'BGR 04/04/2018' to locate comment).")
  endif
  call get_param(param_file, mdl, "VEL_UNDERFLOW", CS%vel_underflow, &
                 "A negligibly small velocity magnitude below which velocity "//&
                 "components are set to 0.  A reasonable value might be "//&
                 "1e-30 m/s, which is less than an Angstrom divided by "//&
                 "the age of the universe.", units="m s-1", default=0.0, scale=US%m_s_to_L_T)

  ALLOC_(CS%a_u(IsdB:IedB,jsd:jed,nz+1)) ; CS%a_u(:,:,:) = 0.0
  ALLOC_(CS%h_u(IsdB:IedB,jsd:jed,nz))   ; CS%h_u(:,:,:) = 0.0
  ALLOC_(CS%a_v(isd:ied,JsdB:JedB,nz+1)) ; CS%a_v(:,:,:) = 0.0
  ALLOC_(CS%h_v(isd:ied,JsdB:JedB,nz))   ; CS%h_v(:,:,:) = 0.0

  CS%id_Kv_slow = register_diag_field('ocean_model', 'Kv_slow', diag%axesTi, Time, &
     'Slow varying vertical viscosity', 'm2 s-1', conversion=US%Z2_T_to_m2_s)

  CS%id_Kv_u = register_diag_field('ocean_model', 'Kv_u', diag%axesCuL, Time, &
     'Total vertical viscosity at u-points', 'm2 s-1', conversion=US%Z2_T_to_m2_s)

  CS%id_Kv_v = register_diag_field('ocean_model', 'Kv_v', diag%axesCvL, Time, &
     'Total vertical viscosity at v-points', 'm2 s-1', conversion=US%Z2_T_to_m2_s)

  CS%id_au_vv = register_diag_field('ocean_model', 'au_visc', diag%axesCui, Time, &
     'Zonal Viscous Vertical Coupling Coefficient', 'm s-1', conversion=US%Z_to_m*US%s_to_T)

  CS%id_av_vv = register_diag_field('ocean_model', 'av_visc', diag%axesCvi, Time, &
     'Meridional Viscous Vertical Coupling Coefficient', 'm s-1', conversion=US%Z_to_m*US%s_to_T)

  CS%id_h_u = register_diag_field('ocean_model', 'Hu_visc', diag%axesCuL, Time, &
     'Thickness at Zonal Velocity Points for Viscosity', thickness_units, &
     conversion=GV%H_to_m)

  CS%id_h_v = register_diag_field('ocean_model', 'Hv_visc', diag%axesCvL, Time, &
     'Thickness at Meridional Velocity Points for Viscosity', thickness_units, &
     conversion=GV%H_to_m)

  CS%id_hML_u = register_diag_field('ocean_model', 'HMLu_visc', diag%axesCu1, Time, &
     'Mixed Layer Thickness at Zonal Velocity Points for Viscosity', thickness_units, &
     conversion=GV%H_to_m)

  CS%id_hML_v = register_diag_field('ocean_model', 'HMLv_visc', diag%axesCv1, Time, &
     'Mixed Layer Thickness at Meridional Velocity Points for Viscosity', thickness_units, &
     conversion=GV%H_to_m)

  CS%id_du_dt_visc = register_diag_field('ocean_model', 'du_dt_visc', diag%axesCuL, &
     Time, 'Zonal Acceleration from Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_du_dt_visc > 0) call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
  CS%id_dv_dt_visc = register_diag_field('ocean_model', 'dv_dt_visc', diag%axesCvL, &
     Time, 'Meridional Acceleration from Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_dv_dt_visc > 0) call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)

  CS%id_taux_bot = register_diag_field('ocean_model', 'taux_bot', diag%axesCu1, &
     Time, 'Zonal Bottom Stress from Ocean to Earth', 'Pa', &
     conversion=US%RZ_to_kg_m2*US%L_T2_to_m_s2)
  CS%id_tauy_bot = register_diag_field('ocean_model', 'tauy_bot', diag%axesCv1, &
     Time, 'Meridional Bottom Stress from Ocean to Earth', 'Pa', &
     conversion=US%RZ_to_kg_m2*US%L_T2_to_m_s2)

  !CS%id_hf_du_dt_visc = register_diag_field('ocean_model', 'hf_du_dt_visc', diag%axesCuL, Time, &
  !    'Fractional Thickness-weighted Zonal Acceleration from Vertical Viscosity', 'm s-2', &
  !    v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_du_dt_visc > 0) then
  !  call safe_alloc_ptr(CS%hf_du_dt_visc,IsdB,IedB,jsd,jed,nz)
  !  call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
  !  call safe_alloc_ptr(ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)
  !endif

  !CS%id_hf_dv_dt_visc = register_diag_field('ocean_model', 'hf_dv_dt_visc', diag%axesCvL, Time, &
  !    'Fractional Thickness-weighted Meridional Acceleration from Vertical Viscosity', 'm s-2', &
  !    v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_dv_dt_visc > 0) then
  !  call safe_alloc_ptr(CS%hf_dv_dt_visc,isd,ied,JsdB,JedB,nz)
  !  call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  !  call safe_alloc_ptr(ADp%diag_hfrac_v,isd,ied,Jsd,JedB,nz)
  !endif

  CS%id_hf_du_dt_visc_2d = register_diag_field('ocean_model', 'hf_du_dt_visc_2d', diag%axesCu1, Time, &
      'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Vertical Viscosity', 'm s-2', &
      conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_du_dt_visc_2d > 0) then
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_hf_dv_dt_visc_2d = register_diag_field('ocean_model', 'hf_dv_dt_visc_2d', diag%axesCv1, Time, &
      'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Vertical Viscosity', 'm s-2', &
      conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_dv_dt_visc_2d > 0) then
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%diag_hfrac_v,isd,ied,Jsd,JedB,nz)
  endif

  if ((len_trim(CS%u_trunc_file) > 0) .or. (len_trim(CS%v_trunc_file) > 0)) &
    call PointAccel_init(MIS, Time, G, param_file, diag, dirs, CS%PointAccel_CSp)

end subroutine vertvisc_init

!> Update the CFL truncation value as a function of time.
!! If called with the optional argument activate=.true., record the
!! value of Time as the beginning of the ramp period.
subroutine updateCFLtruncationValue(Time, CS, activate)
  type(time_type), target, intent(in)    :: Time     !< Current model time
  type(vertvisc_CS),       pointer       :: CS       !< Vertical viscosity control structure
  logical, optional,       intent(in)    :: activate !< Specifiy whether to record the value of
                                                     !! Time as the beginning of the ramp period

  ! Local variables
  real :: deltaTime, wghtA
  character(len=12) :: msg

  if (CS%truncRampTime==0.) return ! This indicates to ramping is turned off

  ! We use the optional argument to indicate this Time should be recorded as the
  ! beginning of the ramp-up period.
  if (present(activate)) then
    if (activate) then
      CS%rampStartTime = Time ! Record the current time
      CS%CFLrampingIsActivated = .true.
    endif
  endif
  if (.not.CS%CFLrampingIsActivated) return
  deltaTime = max( 0., time_type_to_real( Time - CS%rampStartTime ) )
  if (deltaTime >= CS%truncRampTime) then
    CS%CFL_trunc = CS%CFL_truncE
    CS%truncRampTime = 0. ! This turns off ramping after this call
  else
    wghtA = min( 1., deltaTime / CS%truncRampTime ) ! Linear profile in time
    !wghtA = wghtA*wghtA ! Convert linear profile to parabolic profile in time
    !wghtA = wghtA*wghtA*(3. - 2.*wghtA) ! Convert linear profile to cosine profile
    wghtA = 1. - ( (1. - wghtA)**2 ) ! Convert linear profiel to nverted parabolic profile
    CS%CFL_trunc = CS%CFL_truncS + wghtA * ( CS%CFL_truncE - CS%CFL_truncS )
  endif
  write(msg(1:12),'(es12.3)') CS%CFL_trunc
  call MOM_error(NOTE, "MOM_vert_friction: updateCFLtruncationValue set CFL"// &
                       " limit to "//trim(msg))
end subroutine updateCFLtruncationValue

!> Clean up and deallocate the vertical friction module
subroutine vertvisc_end(CS)
  type(vertvisc_CS), pointer :: CS !< Vertical viscosity control structure that
                                   !! will be deallocated in this subroutine.

  DEALLOC_(CS%a_u) ; DEALLOC_(CS%h_u)
  DEALLOC_(CS%a_v) ; DEALLOC_(CS%h_v)
  if (associated(CS%a1_shelf_u)) deallocate(CS%a1_shelf_u)
  if (associated(CS%a1_shelf_v)) deallocate(CS%a1_shelf_v)
  deallocate(CS)
end subroutine vertvisc_end

!> \namespace mom_vert_friction
!! \author Robert Hallberg
!! \date April 1994 - October 2006
!!
!!  The vertical diffusion of momentum is fully implicit.  This is
!!  necessary to allow for vanishingly small layers.  The coupling
!!  is based on the distance between the centers of adjacent layers,
!!  except where a layer is close to the bottom compared with a
!!  bottom boundary layer thickness when a bottom drag law is used.
!!  A stress top b.c. and a no slip bottom  b.c. are used.  There
!!  is no limit on the time step for vertvisc.
!!
!!  Near the bottom, the horizontal thickness interpolation scheme
!!  changes to an upwind biased estimate to control the effect of
!!  spurious Montgomery potential gradients at the bottom where
!!  nearly massless layers layers ride over the topography.  Within a
!!  few boundary layer depths of the bottom, the harmonic mean
!!  thickness (i.e. (2 h+ h-) / (h+ + h-) ) is used if the velocity
!!  is from the thinner side and the arithmetic mean thickness
!!  (i.e. (h+ + h-)/2) is used if the velocity is from the thicker
!!  side.  Both of these thickness estimates are second order
!!  accurate.  Above this the arithmetic mean thickness is used.
!!
!!  In addition, vertvisc truncates any velocity component that
!!  exceeds maxvel to truncvel. This basically keeps instabilities
!!  spatially localized.  The number of times the velocity is
!!  truncated is reported each time the energies are saved, and if
!!  exceeds CS%Maxtrunc the model will stop itself and change the time
!!  to a large value.  This has proven very useful in (1) diagnosing
!!  model failures and (2) letting the model settle down to a
!!  meaningful integration from a poorly specified initial condition.
!!
!!  The same code is used for the two velocity components, by
!!  indirectly referencing the velocities and defining a handful of
!!  direction-specific defined variables.
!!
!!  Macros written all in capital letters are defined in MOM_memory.h.
!!
!!     A small fragment of the grid is shown below:
!! \verbatim
!!    j+1  x ^ x ^ x   At x:  q
!!    j+1  > o > o >   At ^:  v, frhatv, tauy
!!    j    x ^ x ^ x   At >:  u, frhatu, taux
!!    j    > o > o >   At o:  h
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!! \endverbatim
!!
!!  The boundaries always run through q grid points (x).
end module MOM_vert_friction
