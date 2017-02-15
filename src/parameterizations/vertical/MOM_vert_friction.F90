!> Implements vertical viscosity (vertvisc)

module MOM_vert_friction
! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl
use MOM_debugging, only : uchksum, vchksum, hchksum
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_SIMPLE
use MOM_PointAccel, only : write_u_accel, write_v_accel, PointAccel_init
use MOM_PointAccel, only : PointAccel_CS
use MOM_time_manager, only : time_type, time_type_to_real, operator(-)
use MOM_variables, only : thermo_var_ptrs, vertvisc_type
use MOM_variables, only : cont_diag_ptrs, accel_diag_ptrs
use MOM_variables, only : ocean_internal_state
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public vertvisc, vertvisc_remnant, vertvisc_coef
public vertvisc_limit_vel, vertvisc_init, vertvisc_end
public updateCFLtruncationValue

type, public :: vertvisc_CS ; private
  real    :: Hmix            !< The mixed layer thickness in m.
  real    :: Hmix_stress     !< The mixed layer thickness over which the wind
                             !! stress is applied with direct_stress, in m.
  real    :: Kvml            !< The mixed layer vertical viscosity in m2 s-1.
  real    :: Kv              !< The interior vertical viscosity in m2 s-1.
  real    :: Hbbl            !< The static bottom boundary layer thickness, in m.
  real    :: Kvbbl           !< The vertical viscosity in the bottom boundary
                             !! layer, in m2 s-1.

  real    :: maxvel          !< Velocity components greater than maxvel,
                             !! in m s-1, are truncated.
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
    a_u                !< The u-drag coefficient across an interface, in m s-1.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    h_u                !< The effective layer thickness at u-points, m or kg m-2.
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NK_INTERFACE_) :: &
    a_v                !< The v-drag coefficient across an interface, in m s-1.
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    h_v                !< The effective layer thickness at v-points, m or kg m-2.
  !>@{
  !! The surface coupling coefficient under ice shelves
  !! in m s-1. Retained to determine stress under shelves.
  real, pointer, dimension(:,:) :: &
    a1_shelf_u => NULL(), &
    a1_shelf_v => NULL()
  !>@}

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
  logical :: debug          !< If true, write verbose checksums for debugging purposes.
  integer :: nkml           !< The number of layers in the mixed layer.
  integer, pointer :: ntrunc  !< The number of times the velocity has been
                              !! truncated since the last call to write_energy.
  !>@{
  !! The complete path to files in which a column's worth of
  !! accelerations are written when velocity truncations occur.
  character(len=200) :: u_trunc_file
  character(len=200) :: v_trunc_file
  !>@}

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  !>@{
  !! Diagnostic identifiers
  integer :: id_du_dt_visc = -1, id_dv_dt_visc = -1, id_au_vv = -1, id_av_vv = -1
  integer :: id_h_u = -1, id_h_v = -1, id_hML_u = -1 , id_hML_v = -1
  integer :: id_Ray_u = -1, id_Ray_v = -1, id_taux_bot = -1, id_tauy_bot = -1
  !>@}

  type(PointAccel_CS), pointer :: PointAccel_CSp => NULL()
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

subroutine vertvisc(u, v, h, fluxes, visc, dt, OBC, ADp, CDp, G, GV, CS, &
                    taux_bot, tauy_bot)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  real, intent(inout), &
    dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: u      !< Zonal velocity in m s-1
  real, intent(inout), &
    dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: v      !< Meridional velocity in m s-1
  real, intent(in), &
    dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: h      !< Layer thickness in H
  type(forcing),         intent(in)      :: fluxes !< Structure containing forcing fields
  type(vertvisc_type),   intent(inout)   :: visc   !< Viscosities and bottom drag
  real,                  intent(in)      :: dt     !< Time increment in s
  type(ocean_OBC_type),  pointer         :: OBC    !< Open boundary condition structure
  type(accel_diag_ptrs), intent(inout)   :: ADp    !< Accelerations in the momentum
                                                   !! equations for diagnostics
  type(cont_diag_ptrs),  intent(inout)   :: CDp    !< Continuity equation terms
  type(vertvisc_CS),     pointer         :: CS     !< Vertical viscosity control structure
  !> Zonal bottom stress from ocean to rock in Pa
  real, optional, intent(out), dimension(SZIB_(G),SZJ_(G)) :: taux_bot
  !> Meridional bottom stress from ocean to rock in Pa
  real, optional, intent(out), dimension(SZI_(G),SZJB_(G)) :: tauy_bot

  ! Fields from fluxes used in this subroutine:
  !   taux: Zonal wind stress in Pa.
  !   tauy: Meridional wind stress in Pa.

  ! Local variables

  real :: b1(SZIB_(G))          ! b1 and c1 are variables used by the
  real :: c1(SZIB_(G),SZK_(G))  ! tridiagonal solver.  c1 is nondimensional,
                                ! while b1 has units of inverse thickness.
  real :: d1(SZIB_(G))          ! d1=1-c1 is used by the tridiagonal solver, ND.
  real :: Ray(SZIB_(G),SZK_(G)) ! Ray is the Rayleigh-drag velocity in m s-1
  real :: b_denom_1             ! The first term in the denominator of b1, in H.

  real :: Hmix             ! The mixed layer thickness over which stress
                           ! is applied with direct_stress, translated into
                           ! thickness units - either m or kg m-2.
  real :: I_Hmix           ! The inverse of Hmix, in m-1 or m2 kg-1.
  real :: Idt              ! The inverse of the time step, in s-1.
  real :: dt_Rho0          ! The time step divided by the mean
                           ! density, in s m3 kg-1.
  real :: Rho0             ! A density used to convert drag laws into stress in
                           ! Pa, in kg m-3.
  real :: dt_m_to_H        ! The time step times the conversion from m to the
                           ! units of thickness - either s or s m3 kg-1.
  real :: h_neglect        ! A thickness that is so small it is usually lost
                           ! in roundoff and can be neglected, in m or kg m-2.

  real :: stress           !   The surface stress times the time step, divided
                           ! by the density, in units of m2 s-1.
  real :: zDS, hfr, h_a    ! Temporary variables used with direct_stress.
  real :: surface_stress(SZIB_(G))! The same as stress, unless the wind
                           ! stress is applied as a body force, in
                           ! units of m2 s-1.

  logical :: do_i(SZIB_(G))

  integer :: i, j, k, is, ie, Isq, Ieq, Jsq, Jeq, nz, n
  is = G%isc ; ie = G%iec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = G%ke

  if (.not.associated(CS)) call MOM_error(FATAL,"MOM_vert_friction(visc): "// &
         "Module must be initialized before it is used.")

  if (CS%direct_stress) then
    Hmix = CS%Hmix_stress*GV%m_to_H
    I_Hmix = 1.0 / Hmix
  endif
  dt_Rho0 = dt/GV%H_to_kg_m2
  dt_m_to_H = dt*GV%m_to_H
  Rho0 = GV%Rho0
  h_neglect = GV%H_subroundoff
  Idt = 1.0 / dt

  do k=1,nz ; do i=Isq,Ieq ; Ray(i,k) = 0.0 ; enddo ; enddo

  !   Update the zonal velocity component using a modification of a standard
  ! tridagonal solver.
!$OMP parallel do default(none) shared(G,Isq,Ieq,ADp,nz,u,CS,dt_Rho0,fluxes,h, &
!$OMP                                  h_neglect,Hmix,I_Hmix,visc,dt_m_to_H,   &
!$OMP                                  Idt,taux_bot,Rho0)                      &
!$OMP                     firstprivate(Ray)                                    &
!$OMP                          private(do_i,surface_stress,zDS,stress,h_a,hfr, &
!$OMP                                     b_denom_1,b1,d1,c1)
  do j=G%jsc,G%jec
    do I=Isq,Ieq ; do_i(I) = (G%mask2dCu(I,j) > 0) ; enddo

    if (ASSOCIATED(ADp%du_dt_visc)) then ; do k=1,nz ; do I=Isq,Ieq
      ADp%du_dt_visc(I,j,k) = u(I,j,k)
    enddo ; enddo ; endif

!   One option is to have the wind stress applied as a body force
! over the topmost Hmix fluid.  If DIRECT_STRESS is not defined,
! the wind stress is applied as a stress boundary condition.
    if (CS%direct_stress) then
      do I=Isq,Ieq ; if (do_i(I)) then
        surface_stress(I) = 0.0
        zDS = 0.0
        stress = dt_Rho0 * fluxes%taux(I,j)
        do k=1,nz
          h_a = 0.5 * (h(I,j,k) + h(I+1,j,k)) + h_neglect
          hfr = 1.0 ; if ((zDS+h_a) > Hmix) hfr = (Hmix - zDS) / h_a
          u(I,j,k) = u(I,j,k) + I_Hmix * hfr * stress
          zDS = zDS + h_a ; if (zDS >= Hmix) exit
        enddo
      endif ; enddo ! end of i loop
    else ; do I=Isq,Ieq
      surface_stress(I) = dt_Rho0 * (G%mask2dCu(I,j)*fluxes%taux(I,j))
    enddo ; endif ! direct_stress

    if (CS%Channel_drag) then ; do k=1,nz ; do I=Isq,Ieq
      Ray(I,k) = visc%Ray_u(I,j,k)
    enddo ; enddo ; endif

    ! perform forward elimination on the tridiagonal system
    !
    ! denote the diagonal of the system as b_i, the subdiagonal as a_i
    ! and the superdiagonal as c_i. The right-hand side terms are d_i.
    !
    ! ignoring the rayleigh drag contribution,
    ! we have a_i = -dt_m_to_H * a_u(i)
    !         b_i = h_u(i) + dt_m_to_H * (a_u(i) + a_u(i+1))
    !         c_i = -dt_m_to_H * a_u(i+1)
    !
    ! for forward elimination, we want to:
    ! calculate c'_i = - c_i                / (b_i + a_i c'_(i-1))
    ! and       d'_i = (d_i - a_i d'_(i-1)) / (b_i + a_i c'_(i-1))
    ! where c'_1 = c_1/b_1 and d'_1 = d_1/b_1
    ! (see Thomas' tridiagonal matrix algorithm)
    !
    ! b1 is the denominator term 1 / (b_i + a_i c'_(i-1))
    ! b_denom_1 is (b_i + a_i + c_i) - a_i(1 - c'_(i-1))
    !            = (b_i + c_i + c'_(i-1))
    ! this is done so that d1 = b1 * b_denom_1 = 1 - c'_(i-1)
    ! c1(i) is -c'_(i - 1)
    ! and the right-hand-side is destructively updated to be d'_i
    !
    do I=Isq,Ieq ; if (do_i(I)) then
      b_denom_1 = CS%h_u(I,j,1) + dt_m_to_H * (Ray(I,1) + CS%a_u(I,j,1))
      b1(I) = 1.0 / (b_denom_1 + dt_m_to_H*CS%a_u(I,j,2))
      d1(I) = b_denom_1 * b1(I)
      u(I,j,1) = b1(I) * (CS%h_u(I,j,1) * u(I,j,1) + surface_stress(I))
    endif ; enddo
    do k=2,nz ; do I=Isq,Ieq ; if (do_i(I)) then
      c1(I,k) = dt_m_to_H * CS%a_u(I,j,K) * b1(I)
      b_denom_1 = CS%h_u(I,j,k) + dt_m_to_H * (Ray(I,k) + CS%a_u(I,j,K)*d1(I))
      b1(I) = 1.0 / (b_denom_1 + dt_m_to_H * CS%a_u(I,j,K+1))
      d1(I) = b_denom_1 * b1(I)
      u(I,j,k) = (CS%h_u(I,j,k) * u(I,j,k) + &
                  dt_m_to_H * CS%a_u(I,j,K) * u(I,j,k-1)) * b1(I)
    endif ; enddo ; enddo

    ! back substitute to solve for the new velocities
    ! u_i = d'_i - c'_i x_(i+1)
    do k=nz-1,1,-1 ; do I=Isq,Ieq ; if (do_i(I)) then
      u(I,j,k) = u(I,j,k) + c1(I,k+1) * u(I,j,k+1)
    endif ; enddo ; enddo ! i and k loops

    if (ASSOCIATED(ADp%du_dt_visc)) then ; do k=1,nz ; do I=Isq,Ieq
      ADp%du_dt_visc(I,j,k) = (u(I,j,k) - ADp%du_dt_visc(I,j,k))*Idt
    enddo ; enddo ; endif

    if (ASSOCIATED(visc%taux_shelf)) then ; do I=Isq,Ieq
      visc%taux_shelf(I,j) = -Rho0*CS%a1_shelf_u(I,j)*u(I,j,1) ! - u_shelf?
    enddo ; endif

    if (PRESENT(taux_bot)) then
      do I=Isq,Ieq
        taux_bot(I,j) = Rho0 * (u(I,j,nz)*CS%a_u(I,j,nz+1))
      enddo
      if (CS%Channel_drag) then ; do k=1,nz ; do I=Isq,Ieq
        taux_bot(I,j) = taux_bot(I,j) + Rho0 * (Ray(I,k)*u(I,j,k))
      enddo ; enddo ; endif
    endif
  enddo ! end u-component j loop

  ! Now work on the meridional velocity component.
!$OMP parallel do default(none) shared(G,Jsq,Jeq,ADp,nz,v,CS,dt_Rho0,fluxes,h, &
!$OMP                                  Hmix,I_Hmix,visc,dt_m_to_H,Idt,Rho0,    &
!$OMP                                  tauy_bot,is,ie)                         &
!$OMP                     firstprivate(Ray)                                    &
!$OMP                          private(do_i,surface_stress,zDS,stress,h_a,hfr, &
!$OMP                                  b_denom_1,b1,d1,c1)
  do J=Jsq,Jeq
    do i=is,ie ; do_i(i) = (G%mask2dCv(i,J) > 0) ; enddo

    if (ASSOCIATED(ADp%dv_dt_visc)) then ; do k=1,nz ; do i=is,ie
      ADp%dv_dt_visc(i,J,k) = v(i,J,k)
    enddo ; enddo ; endif

!   One option is to have the wind stress applied as a body force
! over the topmost Hmix fluid.  If DIRECT_STRESS is not defined,
! the wind stress is applied as a stress boundary condition.
    if (CS%direct_stress) then
      do i=is,ie ; if (do_i(i)) then
        surface_stress(i) = 0.0
        zDS = 0.0
        stress = dt_Rho0 * fluxes%tauy(i,J)
        do k=1,nz
          h_a = 0.5 * (h(i,J,k) + h(i,J+1,k))
          hfr = 1.0 ; if ((zDS+h_a) > Hmix) hfr = (Hmix - zDS) / h_a
          v(i,J,k) = v(i,J,k) + I_Hmix * hfr * stress
          zDS = zDS + h_a ; if (zDS >= Hmix) exit
        enddo
      endif ; enddo ! end of i loop
    else ; do i=is,ie
      surface_stress(i) = dt_Rho0 * (G%mask2dCv(i,J)*fluxes%tauy(i,J))
    enddo ; endif ! direct_stress

    if (CS%Channel_drag) then ; do k=1,nz ; do i=is,ie
      Ray(i,k) = visc%Ray_v(i,J,k)
    enddo ; enddo ; endif

    do i=is,ie ; if (do_i(i)) then
      b_denom_1 = CS%h_v(i,J,1) + dt_m_to_H * (Ray(i,1) + CS%a_v(i,J,1))
      b1(i) = 1.0 / (b_denom_1 + dt_m_to_H*CS%a_v(i,J,2))
      d1(i) = b_denom_1 * b1(i)
      v(i,J,1) = b1(i) * (CS%h_v(i,J,1) * v(i,J,1) + surface_stress(i))
    endif ; enddo
    do k=2,nz ; do i=is,ie ; if (do_i(i)) then
      c1(i,k) = dt_m_to_H * CS%a_v(i,J,K) * b1(i)
      b_denom_1 = CS%h_v(i,J,k) + dt_m_to_H *  (Ray(i,k) + CS%a_v(i,J,K)*d1(i))
      b1(i) = 1.0 / (b_denom_1 + dt_m_to_H * CS%a_v(i,J,K+1))
      d1(i) = b_denom_1 * b1(i)
      v(i,J,k) = (CS%h_v(i,J,k) * v(i,J,k) + dt_m_to_H * &
                  CS%a_v(i,J,K) * v(i,J,k-1)) * b1(i)
    endif ; enddo ; enddo
    do k=nz-1,1,-1 ; do i=is,ie ; if (do_i(i)) then
      v(i,J,k) = v(i,J,k) + c1(i,k+1) * v(i,J,k+1)
    endif ; enddo ; enddo ! i and k loops

    if (ASSOCIATED(ADp%dv_dt_visc)) then ; do k=1,nz ; do i=is,ie
      ADp%dv_dt_visc(i,J,k) = (v(i,J,k) - ADp%dv_dt_visc(i,J,k))*Idt
    enddo ; enddo ; endif

    if (ASSOCIATED(visc%tauy_shelf)) then ; do i=is,ie
      visc%tauy_shelf(i,J) = -Rho0*CS%a1_shelf_v(i,J)*v(i,J,1) ! - v_shelf?
    enddo ; endif

    if (present(tauy_bot)) then
      do i=is,ie
        tauy_bot(i,J) = Rho0 * (v(i,J,nz)*CS%a_v(i,J,nz+1))
      enddo
      if (CS%Channel_drag) then ; do k=1,nz ; do i=is,ie
        tauy_bot(i,J) = tauy_bot(i,J) + Rho0 * (Ray(i,k)*v(i,J,k))
      enddo ; enddo ; endif
    endif
  enddo ! end of v-component J loop

  call vertvisc_limit_vel(u, v, h, ADp, CDp, fluxes, visc, dt, G, GV, CS)

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

end subroutine vertvisc

!> Calculate the fraction of momentum originally in a layer that remains
!! after a time-step of viscosity, and the fraction of a time-step's
!! worth of barotropic acceleration that a layer experiences after
!! viscosity is applied.
subroutine vertvisc_remnant(visc, visc_rem_u, visc_rem_v, dt, G, GV, CS)
  type(ocean_grid_type), intent(in)   :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV   !< Ocean vertical grid structure
  type(vertvisc_type), intent(in)     :: visc !< Viscosities and bottom drag
  !> Fraction of a time-step's worth of a barotopic acceleration that
  !! a layer experiences after viscosity is applied in the zonal direction
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: visc_rem_u
  !> Fraction of a time-step's worth of a barotopic acceleration that
  !! a layer experiences after viscosity is applied in the meridional direction
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: visc_rem_v
  real, intent(in)           :: dt !< Time increment in s
  type(vertvisc_CS), pointer :: CS !< Vertical viscosity control structure

  ! Local variables

  real :: b1(SZIB_(G))          ! b1 and c1 are variables used by the
  real :: c1(SZIB_(G),SZK_(G))  ! tridiagonal solver.  c1 is nondimensional,
                                ! while b1 has units of inverse thickness.
  real :: d1(SZIB_(G))          ! d1=1-c1 is used by the tridiagonal solver, ND.
  real :: Ray(SZIB_(G),SZK_(G)) ! Ray is the Rayleigh-drag velocity times the
                                ! time step, in m.
  real :: b_denom_1   ! The first term in the denominator of b1, in m or kg m-2.
  real :: dt_m_to_H        ! The time step times the conversion from m to the
                           ! units of thickness - either s or s m3 kg-1.
  logical :: do_i(SZIB_(G))

  integer :: i, j, k, is, ie, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = G%ke

  if (.not.associated(CS)) call MOM_error(FATAL,"MOM_vert_friction(visc): "// &
         "Module must be initialized before it is used.")

  dt_m_to_H = dt*GV%m_to_H

  do k=1,nz ; do i=Isq,Ieq ; Ray(i,k) = 0.0 ; enddo ; enddo

  ! Find the zonal viscous using a modification of a standard tridagonal solver.
!$OMP parallel do default(none) shared(G,Isq,Ieq,CS,nz,visc,dt_m_to_H,visc_rem_u) &
!$OMP                     firstprivate(Ray)                                       &
!$OMP                          private(do_i,b_denom_1,b1,d1,c1)
  do j=G%jsc,G%jec
    do I=Isq,Ieq ; do_i(I) = (G%mask2dCu(I,j) > 0) ; enddo

    if (CS%Channel_drag) then ; do k=1,nz ; do I=Isq,Ieq
      Ray(I,k) = visc%Ray_u(I,j,k)
    enddo ; enddo ; endif

    do I=Isq,Ieq ; if (do_i(I)) then
      b_denom_1 = CS%h_u(I,j,1) + dt_m_to_H * (Ray(I,1) + CS%a_u(I,j,1))
      b1(I) = 1.0 / (b_denom_1 + dt_m_to_H*CS%a_u(I,j,2))
      d1(I) = b_denom_1 * b1(I)
      visc_rem_u(I,j,1) = b1(I) * CS%h_u(I,j,1)
    endif ; enddo
    do k=2,nz ; do I=Isq,Ieq ; if (do_i(I)) then
      c1(I,k) = dt_m_to_H * CS%a_u(I,j,K)*b1(I)
      b_denom_1 = CS%h_u(I,j,k) + dt_m_to_H * (Ray(I,k) + CS%a_u(I,j,K)*d1(I))
      b1(I) = 1.0 / (b_denom_1 + dt_m_to_H * CS%a_u(I,j,K+1))
      d1(I) = b_denom_1 * b1(I)
      visc_rem_u(I,j,k) = (CS%h_u(I,j,k) + dt_m_to_H * CS%a_u(I,j,K) * visc_rem_u(I,j,k-1)) * b1(I)
    endif ; enddo ; enddo
    do k=nz-1,1,-1 ; do I=Isq,Ieq ; if (do_i(I)) then
      visc_rem_u(I,j,k) = visc_rem_u(I,j,k) + c1(I,k+1)*visc_rem_u(I,j,k+1)

    endif ; enddo ; enddo ! i and k loops

  enddo ! end u-component j loop

  ! Now find the meridional viscous using a modification.
!$OMP parallel do default(none) shared(Jsq,Jeq,is,ie,G,CS,visc,dt_m_to_H,visc_rem_v,nz) &
!$OMP                     firstprivate(Ray)                                             &
!$OMP                          private(do_i,b_denom_1,b1,d1,c1)
  do J=Jsq,Jeq
    do i=is,ie ; do_i(i) = (G%mask2dCv(i,J) > 0) ; enddo

    if (CS%Channel_drag) then ; do k=1,nz ; do i=is,ie
      Ray(i,k) = visc%Ray_v(i,J,k)
    enddo ; enddo ; endif

    do i=is,ie ; if (do_i(i)) then
      b_denom_1 = CS%h_v(i,J,1) + dt_m_to_H * (Ray(i,1) + CS%a_v(i,J,1))
      b1(i) = 1.0 / (b_denom_1 + dt_m_to_H*CS%a_v(i,J,2))
      d1(i) = b_denom_1 * b1(i)
      visc_rem_v(i,J,1) = b1(i) * CS%h_v(i,J,1)
    endif ; enddo
    do k=2,nz ; do i=is,ie ; if (do_i(i)) then
      c1(i,k) = dt_m_to_H * CS%a_v(i,J,K)*b1(i)
      b_denom_1 = CS%h_v(i,J,k) + dt_m_to_H *  (Ray(i,k) + CS%a_v(i,J,K)*d1(i))
      b1(i) = 1.0 / (b_denom_1 + dt_m_to_H * CS%a_v(i,J,K+1))
      d1(i) = b_denom_1 * b1(i)
      visc_rem_v(i,J,k) = (CS%h_v(i,J,k) + dt_m_to_H * CS%a_v(i,J,K) * visc_rem_v(i,J,k-1)) * b1(i)
    endif ; enddo ; enddo
    do k=nz-1,1,-1 ; do i=is,ie ; if (do_i(i)) then
      visc_rem_v(i,J,k) = visc_rem_v(i,J,k) + c1(i,k+1)*visc_rem_v(i,J,k+1)
    endif ; enddo ; enddo ! i and k loops
  enddo ! end of v-component J loop

  if (CS%debug) then
    call uchksum(visc_rem_u,"visc_rem_u",G%HI,haloshift=0)
    call vchksum(visc_rem_v,"visc_rem_v",G%HI,haloshift=0)
  endif

end subroutine vertvisc_remnant


!> Calculate the coupling coefficients (CS%a_u and CS%a_v)
!! and effective layer thicknesses (CS%h_u and CS%h_v) for later use in the
!! applying the implicit vertical viscosity via vertvisc().
subroutine vertvisc_coef(u, v, h, fluxes, visc, dt, G, GV, CS)
  type(ocean_grid_type), intent(in)      :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  real, intent(in), &
    dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: u      !< Zonal velocity in m s-1
  real, intent(in), &
    dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: v      !< Meridional velocity in m s-1
  real, intent(in), &
    dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: h      !< Layer thickness in H
  type(forcing), intent(in)              :: fluxes !< Structure containing forcing fields
  type(vertvisc_type), intent(in)        :: visc   !< Viscosities and bottom drag
  real, intent(in)                       :: dt     !< Time increment in s
  type(vertvisc_CS), pointer             :: CS     !< Vertical viscosity control structure

  ! Field from fluxes used in this subroutine:
  !   ustar: the friction velocity in m s-1, used here as the mixing
  !     velocity in the mixed layer if NKML > 1 in a bulk mixed layer.

  ! Local variables

  real, dimension(SZIB_(G),SZK_(G)) :: &
    h_harm, &   ! Harmonic mean of the thicknesses around a velocity grid point,
                ! given by 2*(h+ * h-)/(h+ + h-), in m or kg m-2 (H for short).
    hvel, &     ! hvel is the thickness used at a velocity grid point, in H.
    hvel_shelf  ! The equivalent of hvel under shelves, in H.
  real, dimension(SZIB_(G),SZK_(G)+1) :: &
    a, &        ! The drag coefficients across interfaces, in m s-1.  a times
                ! the velocity difference gives the stress across an interface.
    a_shelf, &  ! The drag coefficients across interfaces in water columns under
                ! ice shelves, in m s-1.
    z_i         ! An estimate of each interface's height above the bottom,
                ! normalized by the bottom boundary layer thickness, nondim.
  real, dimension(SZIB_(G)) :: &
    kv_bbl, &     ! The bottom boundary layer viscosity in m2 s-1.
    bbl_thick, &  ! The bottom boundary layer thickness in m or kg m-2.
    I_Hbbl, &     ! The inverse of the bottom boundary layer thickness, in units
                  ! of H-1 (i.e., m-1 or m2 kg-1).
    I_Htbl, &     ! The inverse of the top boundary layer thickness, in units
                  ! of H-1 (i.e., m-1 or m2 kg-1).
    zcol1, &      ! The height of the interfaces to the north and south of a
    zcol2, &      ! v-point, in m or kg m-2.
    Ztop_min, &   ! The deeper of the two adjacent surface heights, in H.
    Dmin, &       ! The shallower of the two adjacent bottom depths converted to
                  ! thickness units, in m or kg m-2.
    zh, &         ! An estimate of the interface's distance from the bottom
                  ! based on harmonic mean thicknesses, in m or kg m-2.
    h_ml          ! The mixed layer depth, in m or kg m-2.
  real, allocatable, dimension(:,:) :: hML_u, hML_v
  real :: zcol(SZI_(G)) ! The height of an interface at h-points, in m or kg m-2.
  real :: botfn   ! A function which goes from 1 at the bottom to 0 much more
                  ! than Hbbl into the interior.
  real :: topfn   ! A function which goes from 1 at the top to 0 much more
                  ! than Htbl into the interior.
  real :: z2      ! The distance from the bottom, normalized by Hbbl, nondim.
  real :: z2_wt   ! A nondimensional (0-1) weight used when calculating z2.
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected, in H.
  real :: H_to_m, m_to_H ! Unit conversion factors.

  real :: h_arith ! The arithmetic mean thickness, in m or kg m-2.
  real :: I_valBL ! The inverse of a scaling factor determining when water is
                  ! still within the boundary layer, as determined by the sum
                  ! of the harmonic mean thicknesses.
  logical, dimension(SZIB_(G)) :: do_i, do_i_shelf
  logical :: do_any_shelf

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = G%ke

  if (.not.associated(CS)) call MOM_error(FATAL,"MOM_vert_friction(coef): "// &
         "Module must be initialized before it is used.")

  h_neglect = GV%H_subroundoff
  H_to_m = GV%H_to_m ; m_to_H = GV%m_to_H
  I_Hbbl(:) = 1.0 / (CS%Hbbl * GV%m_to_H + h_neglect)
  I_valBL = 0.0 ; if (CS%harm_BL_val > 0.0) I_valBL = 1.0 / CS%harm_BL_val

  if (CS%debug .or. (CS%id_hML_u > 0)) then
    allocate(hML_u(G%IsdB:G%IedB,G%jsd:G%jed)) ; hML_u(:,:) = 0.0
  endif
  if (CS%debug .or. (CS%id_hML_v > 0)) then
    allocate(hML_v(G%isd:G%ied,G%JsdB:G%JedB)) ; hML_v(:,:) = 0.0
  endif

  if ((associated(visc%taux_shelf) .or. associated(fluxes%frac_shelf_u)) .and. &
      .not.associated(CS%a1_shelf_u)) then
    allocate(CS%a1_shelf_u(G%IsdB:G%IedB,G%jsd:G%jed)) ; CS%a1_shelf_u(:,:)=0.0
  endif
  if ((associated(visc%tauy_shelf) .or. associated(fluxes%frac_shelf_u)) .and. &
      .not.associated(CS%a1_shelf_v)) then
    allocate(CS%a1_shelf_v(G%isd:G%ied,G%JsdB:G%JedB)) ; CS%a1_shelf_v(:,:)=0.0
  endif

!$OMP parallel do default(none) shared(G,GV,CS,visc,Isq,ieq,nz,u,h,fluxes,hML_u, &
!$OMP                                  h_neglect,dt,m_to_H,I_valBL) &
!$OMP                     firstprivate(i_hbbl)                                             &
!$OMP                          private(do_i,kv_bbl,bbl_thick,z_i,h_harm,h_arith,hvel,z2,   &
!$OMP                                  botfn,zh,Dmin,zcol,a,do_any_shelf,do_i_shelf,       &
!$OMP                                  a_shelf,Ztop_min,I_HTbl,hvel_shelf,topfn,h_ml,z2_wt)
  do j=G%Jsc,G%Jec
    do I=Isq,Ieq ; do_i(I) = (G%mask2dCu(I,j) > 0) ; enddo

    if (CS%bottomdraglaw) then ; do I=Isq,Ieq
      kv_bbl(I) = visc%kv_bbl_u(I,j)
      bbl_thick(I) = visc%bbl_thick_u(I,j) * m_to_H
      if (do_i(I)) I_Hbbl(I) = 1.0 / (bbl_thick(I) + h_neglect)
    enddo ; endif

!    The following block calculates the thicknesses at velocity
!  grid points for the vertical viscosity (hvel[k]).  Near the
!  bottom an upwind biased thickness is used to control the effect
!  of spurious Montgomery potential gradients at the bottom where
!  nearly massless layers layers ride over the topography.
    if (CS%harmonic_visc) then
      do I=Isq,Ieq ; z_i(I,nz+1) = 0.0 ; enddo
      do k=nz,1,-1 ; do I=Isq,Ieq ; if (do_i(I)) then
        h_harm(I,k) = 2.0*h(i,j,k)*h(i+1,j,k) / (h(i,j,k)+h(i+1,j,k)+h_neglect)
        h_arith = 0.5*(h(i+1,j,k)+h(i,j,k))

        hvel(I,k) = h_harm(I,k)
        if (u(I,j,k) * (h(i+1,j,k)-h(i,j,k)) < 0) then
          z2 = z_i(I,k+1) ; botfn = 1.0 / (1.0 + 0.09*z2*z2*z2*z2*z2*z2)
          hvel(I,k) = (1.0-botfn)*h_harm(I,k) + botfn*h_arith
        endif
        z_i(I,k) =  z_i(I,k+1) + h_harm(I,k)*I_Hbbl(I)
      endif ; enddo ; enddo ! i & k loops
    else ! Not harmonic_visc
      do I=Isq,Ieq
        zh(I) = 0.0 ; z_i(I,nz+1) = 0.0
        Dmin(I) = min(G%bathyT(i,j), G%bathyT(i+1,j)) * m_to_H
      enddo
      do i=Isq,Ieq+1 ; zcol(i) = -G%bathyT(i,j) * m_to_H ; enddo
      do k=nz,1,-1
        do i=Isq,Ieq+1 ; zcol(i) = zcol(i) + h(i,j,k) ; enddo
        do I=Isq,Ieq ; if (do_i(I)) then
          h_harm(I,k) = 2.0*h(i,j,k)*h(i+1,j,k) / (h(i,j,k)+h(i+1,j,k)+h_neglect)
          h_arith = 0.5*(h(i+1,j,k)+h(i,j,k))
          zh(I) = zh(I) + h_harm(I,k)
          z_i(I,k) = max(zh(I), max(zcol(i),zcol(i+1)) + Dmin(I)) * I_Hbbl(I)

          hvel(I,k) = h_arith
          if (u(I,j,k) * (h(i+1,j,k)-h(i,j,k)) > 0) then
            if (zh(I) * I_Hbbl(I) < CS%harm_BL_val) then
              hvel(I,k) = h_harm(I,k)
            else
              z2_wt = 1.0  ; if (zh(I) * I_Hbbl(I) < 2.0*CS%harm_BL_val) &
                z2_wt = max(0.0, min(1.0, zh(I) * I_Hbbl(I) * I_valBL - 1.0))
              z2 = z2_wt * (max(zh(I), max(zcol(i),zcol(i+1)) + Dmin(I)) * I_Hbbl(I))
              botfn = 1.0 / (1.0 + 0.09*z2*z2*z2*z2*z2*z2)
              hvel(I,k) = (1.0-botfn)*h_arith + botfn*h_harm(I,k)
            endif
          endif

        endif ; enddo ! i  loop
      enddo ! k loop
    endif

    call find_coupling_coef(a, hvel, do_i, h_harm, bbl_thick, kv_bbl, z_i, h_ml, &
                            dt, j, G, GV, CS, visc, fluxes, work_on_u=.true.)
    if (allocated(hML_u)) then
        do i=isq,ieq ; if (do_i(i)) then ; hML_u(I,j) = h_ml(I) ; endif ; enddo
    endif
    do_any_shelf = .false.
    if (associated(fluxes%frac_shelf_u)) then
      do I=Isq,Ieq
        CS%a1_shelf_u(I,j) = 0.0
        do_i_shelf(I) = (do_i(I) .and. fluxes%frac_shelf_u(I,j) > 0.0)
        if (do_i_shelf(I)) do_any_shelf = .true.
      enddo
      if (do_any_shelf) then
        if (CS%harmonic_visc) then
          call find_coupling_coef(a_shelf, hvel, do_i_shelf, h_harm, bbl_thick, &
                                  kv_bbl, z_i, h_ml, dt, j, G, GV, CS, &
                                  visc, fluxes, work_on_u=.true., shelf=.true.)
        else  ! Find upwind-biased thickness near the surface.
          ! Perhaps this needs to be done more carefully, via find_eta.
          do I=Isq,Ieq ; if (do_i_shelf(I)) then
            zh(I) = 0.0 ; Ztop_min(I) = min(zcol(i), zcol(i+1))
            I_HTbl(I) = 1.0 / (visc%tbl_thick_shelf_u(I,j)*m_to_H + h_neglect)
          endif ; enddo
          do k=1,nz
            do i=Isq,Ieq+1 ; zcol(i) = zcol(i) - h(i,j,k) ; enddo
            do I=Isq,Ieq ; if (do_i_shelf(I)) then
              h_arith = 0.5*(h(i+1,j,k)+h(i,j,k))
              zh(I) = zh(I) + h_harm(I,k)

              hvel_shelf(I,k) = hvel(I,k)
              if (u(I,j,k) * (h(i+1,j,k)-h(i,j,k)) > 0) then
                if (zh(I) * I_HTbl(I) < CS%harm_BL_val) then
                  hvel_shelf(I,k) = min(hvel(I,k), h_harm(I,k))
                else
                  z2_wt = 1.0  ; if (zh(I) * I_HTbl(I) < 2.0*CS%harm_BL_val) &
                    z2_wt = max(0.0, min(1.0, zh(I) * I_HTbl(I) * I_valBL - 1.0))
                  z2 = z2_wt * (max(zh(I), Ztop_min(I) - min(zcol(i),zcol(i+1))) * I_HTbl(I))
                  topfn = 1.0 / (1.0 + 0.09*z2**6)
                  hvel_shelf(I,k) = min(hvel(I,k), (1.0-topfn)*h_arith + topfn*h_harm(I,k))
                endif
              endif
            endif ; enddo
          enddo
          call find_coupling_coef(a_shelf, hvel_shelf, do_i_shelf, h_harm, &
                                  bbl_thick, kv_bbl, z_i, h_ml, dt, j, G, GV, CS, &
                                  visc, fluxes, work_on_u=.true., shelf=.true.)
        endif
        do I=Isq,Ieq ; if (do_i_shelf(I)) CS%a1_shelf_u(I,j) = a_shelf(I,1) ; enddo
      endif
    endif

    if (do_any_shelf) then
      do K=1,nz+1 ; do I=Isq,Ieq ; if (do_i_shelf(I)) then
        CS%a_u(I,j,K) = fluxes%frac_shelf_u(I,j)  * a_shelf(I,K) + &
                   (1.0-fluxes%frac_shelf_u(I,j)) * a(I,K)
! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
!        CS%a_u(I,j,K) = fluxes%frac_shelf_u(I,j)  * max(a_shelf(I,K), a(I,K)) + &
!                   (1.0-fluxes%frac_shelf_u(I,j)) * a(I,K)
      elseif (do_i(I)) then
        CS%a_u(I,j,K) = a(I,K)
      endif ; enddo ; enddo
      do k=1,nz ; do I=Isq,Ieq ; if (do_i_shelf(I)) then
        ! Should we instead take the inverse of the average of the inverses?
        CS%h_u(I,j,k) = fluxes%frac_shelf_u(I,j)  * hvel_shelf(I,k) + &
                   (1.0-fluxes%frac_shelf_u(I,j)) * hvel(I,k)
      elseif (do_i(I)) then
        CS%h_u(I,j,k) = hvel(I,k)
      endif ; enddo ; enddo
    else
      do K=1,nz+1 ; do I=Isq,Ieq ; if (do_i(I)) CS%a_u(I,j,K) = a(I,K) ; enddo ; enddo
      do k=1,nz ; do I=Isq,Ieq ; if (do_i(I)) CS%h_u(I,j,k) = hvel(I,k) ; enddo ; enddo
    endif

  enddo


  ! Now work on v-points.
!$OMP parallel do default(none) shared(G,GV,CS,visc,is,ie,Jsq,Jeq,nz,v,h,fluxes,hML_v, &
!$OMP                                  h_neglect,dt,m_to_H,I_valBL) &
!$OMP                     firstprivate(i_hbbl)                                                &
!$OMP                          private(do_i,kv_bbl,bbl_thick,z_i,h_harm,h_arith,hvel,z2,      &
!$OMP                                  botfn,zh,Dmin,zcol1,zcol2,a,do_any_shelf,do_i_shelf,  &
!$OMP                                  a_shelf,Ztop_min,I_HTbl,hvel_shelf,topfn,h_ml,z2_wt)
  do J=Jsq,Jeq
    do i=is,ie ; do_i(i) = (G%mask2dCv(i,J) > 0) ; enddo

    if (CS%bottomdraglaw) then ; do i=is,ie
      kv_bbl(i) = visc%kv_bbl_v(i,J)
      bbl_thick(i) = visc%bbl_thick_v(i,J) * m_to_H
      if (do_i(i)) I_Hbbl(i) = 1.0 / bbl_thick(i)
    enddo ; endif

!    The following block calculates the thicknesses at velocity
!  grid points for the vertical viscosity (hvel[k]).  Near the
!  bottom an upwind biased thickness is used to control the effect
!  of spurious Montgomery potential gradients at the bottom where
!  nearly massless layers layers ride over the topography.
    if (CS%harmonic_visc) then
      do i=is,ie ; z_i(i,nz+1) = 0.0 ; enddo

      do k=nz,1,-1 ; do i=is,ie ; if (do_i(i)) then
        h_harm(i,k) = 2.0*h(i,j,k)*h(i,j+1,k) / (h(i,j,k)+h(i,j+1,k)+h_neglect)
        h_arith = 0.5*(h(i,j+1,k)+h(i,j,k))

        hvel(i,k) = h_harm(i,k)
        if (v(i,J,k) * (h(i,j+1,k)-h(i,j,k)) < 0) then
          z2 = z_i(i,k+1) ; botfn = 1.0 / (1.0 + 0.09*z2*z2*z2*z2*z2*z2)
          hvel(i,k) = (1.0-botfn)*h_harm(i,k) + botfn*h_arith
        endif
        z_i(i,k) = z_i(i,k+1)  + h_harm(i,k)*I_Hbbl(i)
      endif ; enddo ; enddo ! i & k loops
    else ! Not harmonic_visc
      do i=is,ie
        zh(i) = 0.0 ; z_i(i,nz+1) = 0.0
        Dmin(i) = min(G%bathyT(i,j), G%bathyT(i,j+1)) * m_to_H
        zcol1(i) = -G%bathyT(i,j) * m_to_H
        zcol2(i) = -G%bathyT(i,j+1) * m_to_H
      enddo
      do k=nz,1,-1 ; do i=is,ie ; if (do_i(i)) then
        h_harm(i,k) = 2.0*h(i,j,k)*h(i,j+1,k) / (h(i,j,k)+h(i,j+1,k)+h_neglect)
        h_arith = 0.5*(h(i,j+1,k)+h(i,j,k))
        zh(i) = zh(i) + h_harm(i,k)
        zcol1(i) = zcol1(i) + h(i,j,k) ; zcol2(i) = zcol2(i) + h(i,j+1,k)
        z_i(I,k) = max(zh(i), max(zcol1(i),zcol2(i)) + Dmin(i)) * I_Hbbl(i)

        hvel(i,k) = h_arith
        if (v(i,J,k) * (h(i,j+1,k)-h(i,j,k)) > 0) then
          if (zh(i) * I_Hbbl(i) < CS%harm_BL_val) then
            hvel(i,k) = h_harm(i,k)
          else
            z2_wt = 1.0  ; if (zh(i) * I_Hbbl(i) < 2.0*CS%harm_BL_val) &
              z2_wt = max(0.0, min(1.0, zh(i) * I_Hbbl(i) * I_valBL - 1.0))
            z2 = z2_wt * (max(zh(i), max(zcol1(i),zcol2(i)) + Dmin(i)) * I_Hbbl(i))
            botfn = 1.0 / (1.0 + 0.09*z2*z2*z2*z2*z2*z2)
            hvel(i,k) = (1.0-botfn)*h_arith + botfn*h_harm(i,k)
          endif
       endif

      endif ; enddo ; enddo ! i & k loops
    endif

    call find_coupling_coef(a, hvel, do_i, h_harm, bbl_thick, kv_bbl, z_i, h_ml, &
                            dt, j, G, GV, CS, visc, fluxes, work_on_u=.false.)
    if ( allocated(hML_v)) then
       do i=is,ie ; if (do_i(i)) then ; hML_v(i,J) = h_ml(i) ; endif ; enddo
    endif
    do_any_shelf = .false.
    if (associated(fluxes%frac_shelf_v)) then
      do i=is,ie
        CS%a1_shelf_v(i,J) = 0.0
        do_i_shelf(i) = (do_i(i) .and. fluxes%frac_shelf_v(i,J) > 0.0)
        if (do_i_shelf(I)) do_any_shelf = .true.
      enddo
      if (do_any_shelf) then
        if (CS%harmonic_visc) then
          call find_coupling_coef(a_shelf, hvel, do_i_shelf, h_harm, bbl_thick, &
                                  kv_bbl, z_i, h_ml, dt, j, G, GV, CS, visc, &
                                  fluxes, work_on_u=.false., shelf=.true.)
        else  ! Find upwind-biased thickness near the surface.
          ! Perhaps this needs to be done more carefully, via find_eta.
          do i=is,ie ; if (do_i_shelf(i)) then
            zh(i) = 0.0 ; Ztop_min(I) = min(zcol1(i), zcol2(i))
            I_HTbl(i) = 1.0 / (visc%tbl_thick_shelf_v(i,J)*m_to_H + h_neglect)
          endif ; enddo
          do k=1,nz
            do i=is,ie ; if (do_i_shelf(i)) then
              zcol1(i) = zcol1(i) - h(i,j,k) ; zcol2(i) = zcol2(i) - h(i,j+1,k)
              h_arith = 0.5*(h(i,j+1,k)+h(i,j,k))
              zh(i) = zh(i) + h_harm(i,k)

              hvel_shelf(i,k) = hvel(i,k)
              if (v(i,j,k) * (h(i,j+1,k)-h(i,j,k)) > 0) then
                if (zh(i) * I_HTbl(i) < CS%harm_BL_val) then
                  hvel_shelf(i,k) = min(hvel(i,k), h_harm(i,k))
                else
                  z2_wt = 1.0  ; if (zh(i) * I_HTbl(i) < 2.0*CS%harm_BL_val) &
                    z2_wt = max(0.0, min(1.0, zh(i) * I_HTbl(i) * I_valBL - 1.0))
                  z2 = z2_wt * (max(zh(i), Ztop_min(i) - min(zcol1(i),zcol2(i))) * I_HTbl(i))
                  topfn = 1.0 / (1.0 + 0.09*z2**6)
                  hvel_shelf(i,k) = min(hvel(i,k), (1.0-topfn)*h_arith + topfn*h_harm(i,k))
                endif
             endif
            endif ; enddo
          enddo
          call find_coupling_coef(a_shelf, hvel_shelf, do_i_shelf, h_harm, &
                                  bbl_thick, kv_bbl, z_i, h_ml, dt, j, G, GV, CS, &
                                  visc, fluxes, work_on_u=.false., shelf=.true.)
        endif
        do i=is,ie ; if (do_i_shelf(i)) CS%a1_shelf_v(i,J) = a_shelf(i,1) ; enddo
      endif
    endif

    if (do_any_shelf) then
      do K=1,nz+1 ; do i=is,ie ; if (do_i_shelf(i)) then
        CS%a_v(i,J,K) = fluxes%frac_shelf_v(i,J)  * a_shelf(i,k) + &
                   (1.0-fluxes%frac_shelf_v(i,J)) * a(i,K)
! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
!        CS%a_v(i,J,K) = fluxes%frac_shelf_v(i,J)  * max(a_shelf(i,K), a(i,K)) + &
!                   (1.0-fluxes%frac_shelf_v(i,J)) * a(i,K)
      elseif (do_i(i)) then
        CS%a_v(i,J,K) = a(i,K)
      endif ; enddo ; enddo
      do k=1,nz ; do i=is,ie ; if (do_i_shelf(i)) then
        ! Should we instead take the inverse of the average of the inverses?
        CS%h_v(i,J,k) = fluxes%frac_shelf_v(i,J)  * hvel_shelf(i,k) + &
                   (1.0-fluxes%frac_shelf_v(i,J)) * hvel(i,k)
      elseif (do_i(i)) then
        CS%h_v(i,J,k) = hvel(i,k)
      endif ; enddo ; enddo
    else
      do K=1,nz+1 ; do i=is,ie ; if (do_i(i)) CS%a_v(i,J,K) = a(i,K) ; enddo ; enddo
      do k=1,nz ; do i=is,ie ; if (do_i(i)) CS%h_v(i,J,k) = hvel(i,k) ; enddo ; enddo
    endif
  enddo ! end of v-point j loop


  if (CS%debug) then
    call uchksum(CS%h_u*H_to_m,"vertvisc_coef h_u",G%HI,haloshift=0)
    call vchksum(CS%h_v*H_to_m,"vertvisc_coef h_v",G%HI,haloshift=0)
    call uchksum(CS%a_u,"vertvisc_coef a_u",G%HI,haloshift=0)
    call vchksum(CS%a_v,"vertvisc_coef a_v",G%HI,haloshift=0)
    if (allocated(hML_u)) call uchksum(hML_u*H_to_m,"vertvisc_coef hML_u",G%HI,haloshift=0)
    if (allocated(hML_v)) call vchksum(hML_v*H_to_m,"vertvisc_coef hML_v",G%HI,haloshift=0)
  endif

! Offer diagnostic fields for averaging.
  if (CS%id_au_vv > 0) call post_data(CS%id_au_vv, CS%a_u, CS%diag)
  if (CS%id_av_vv > 0) call post_data(CS%id_av_vv, CS%a_v, CS%diag)
  if (CS%id_h_u > 0) call post_data(CS%id_h_u, CS%h_u, CS%diag)
  if (CS%id_h_v > 0) call post_data(CS%id_h_v, CS%h_v, CS%diag)
  if (CS%id_hML_u > 0) call post_data(CS%id_hML_u, hML_u, CS%diag)
  if (CS%id_hML_v > 0) call post_data(CS%id_hML_v, hML_v, CS%diag)

  if (allocated(hML_u)) deallocate(hML_u)
  if (allocated(hML_v)) deallocate(hML_v)

end subroutine vertvisc_coef

!> Calculate the 'coupling coefficient' (a[k]) at the
!! interfaces. If BOTTOMDRAGLAW is defined, the minimum of Hbbl and half the
!! adjacent layer thicknesses are used to calculate a[k] near the bottom.
subroutine find_coupling_coef(a, hvel, do_i, h_harm, bbl_thick, kv_bbl, z_i, h_ml, &
                              dt, j, G, GV, CS, visc, fluxes, work_on_u, shelf)
  type(ocean_grid_type), intent(in)                 :: G  !< Ocean grid structure
  type(verticalGrid_type),               intent(in) :: GV !< Ocean vertical grid structure
  !> Coupling coefficient across interfaces, in m s-1
  real,    dimension(SZIB_(G),SZK_(GV)+1), intent(out) :: a
  !> Thickness at velocity points, in H
  real,    dimension(SZIB_(G),SZK_(GV)),   intent(in)  :: hvel
  !> If true, determine coupling coefficient for a column
  logical, dimension(SZIB_(G)),            intent(in)  :: do_i
  !> Harmonic mean of thicknesses around a velocity grid point, in H
  real,    dimension(SZIB_(G),SZK_(GV)),   intent(in)  :: h_harm
  !> Bottom boundary layer thickness, in H
  real,    dimension(SZIB_(G)),            intent(in)  :: bbl_thick
  !> Bottom boundary layer viscosity, in m2 s-1
  real,    dimension(SZIB_(G)),            intent(in)  :: kv_bbl
  !> Estimate of interface heights above the bottom,
  !! normalised by the bottom boundary layer thickness
  real,    dimension(SZIB_(G),SZK_(GV)+1), intent(in)  :: z_i
  !> Mixed layer depth, in H
  real,    dimension(SZIB_(G)),         intent(out) :: h_ml
  !> j-index to find coupling coefficient for
  integer,                              intent(in)  :: j
  !> Time increment, in s
  real,                                 intent(in)  :: dt
  !> Vertical viscosity control structure
  type(vertvisc_CS), pointer                        :: CS
  !> Structure containing viscosities and bottom drag
  type(vertvisc_type), intent(in)                   :: visc
  !> Structure containing forcing fields
  type(forcing), intent(in)                         :: fluxes
  !> If true, u-points are being calculated, otherwise v-points
  logical,                              intent(in)  :: work_on_u
  !> If present and true, use a surface boundary condition
  !! appropriate for an ice shelf.
  logical, optional,                    intent(in)  :: shelf

  ! Local variables

  real, dimension(SZIB_(G)) :: &
    u_star, &   ! ustar at a velocity point, in m s-1.
    absf, &     ! The average of the neighboring absolute values of f, in s-1.
!      h_ml, &     ! The mixed layer depth, in m or kg m-2.
    nk_visc, &  ! The (real) interface index of the base of mixed layer.
    z_t, &      ! The distance from the top, sometimes normalized
                ! by Hmix, in m or nondimensional.
    kv_tbl, &
    tbl_thick
  real :: h_shear ! The distance over which shears occur, m or kg m-2.
  real :: r       ! A thickness to compare with Hbbl, in m or kg m-2.
  real :: visc_ml ! The mixed layer viscosity, in m2 s-1.
  real :: I_Hmix  ! The inverse of the mixed layer thickness, in m-1 or m2 kg-1.
  real :: a_ml    ! The layer coupling coefficient across an interface in
                  ! the mixed layer, in m s-1.
  real :: temp1   ! A temporary variable in m2 s-1.
  real :: h_neglect   ! A thickness that is so small it is usually lost
                      ! in roundoff and can be neglected, in H.
  real :: dz_neglect  ! A thickness in m that is so small it is usually lost
                      ! in roundoff and can be neglected, in m.
  real :: z2      ! A copy of z_i, nondim.
  real :: H_to_m, m_to_H ! Unit conversion factors.
  real :: topfn
  real :: a_top
  logical :: do_shelf
  integer :: i, k, is, ie, max_nk
  integer :: nz
  real    :: botfn

  a(:,:) = 0.0

  if (work_on_u) then ; is = G%IscB ; ie = G%IecB
  else ; is = G%isc ; ie = G%iec ; endif
  nz = G%ke
  h_neglect = GV%H_subroundoff
  H_to_m = GV%H_to_m ; m_to_H = GV%m_to_H
  dz_neglect = GV%H_subroundoff*GV%H_to_m

  do_shelf = .false. ; if (present(shelf)) do_shelf = shelf
  h_ml(:) = 0.0

!    The following loop calculates the vertical average velocity and
!  surface mixed layer contributions to the vertical viscosity.
  do i=is,ie ; a(i,1) = 0.0 ; enddo
  if ((GV%nkml>0) .or. do_shelf) then ; do k=2,nz ; do i=is,ie
    if (do_i(i)) a(i,K) = 2.0*CS%Kv
  enddo ; enddo ; else
    I_Hmix = 1.0 / (CS%Hmix * m_to_H + h_neglect)
    do i=is,ie ; z_t(i) = h_neglect*I_Hmix ; enddo
    do K=2,nz ; do i=is,ie ; if (do_i(i)) then
      z_t(i) = z_t(i) + h_harm(i,k-1)*I_Hmix
      a(i,K) = 2.0*CS%Kv + 2.0*CS%Kvml / ((z_t(i)*z_t(i)) *  &
               (1.0 + 0.09*z_t(i)*z_t(i)*z_t(i)*z_t(i)*z_t(i)*z_t(i)))
    endif ; enddo ; enddo
  endif

  do i=is,ie ; if (do_i(i)) then
    if (CS%bottomdraglaw) then
      r = hvel(i,nz)*0.5
      if (r < bbl_thick(i)) then
        a(i,nz+1) = 1.0*kv_bbl(i) / (1e-10*dt*kv_bbl(i) + r*H_to_m)
      else
        a(i,nz+1) = 1.0*kv_bbl(i) / (1e-10*dt*kv_bbl(i) + bbl_thick(i)*H_to_m)
      endif
    else
      a(i,nz+1) = 2.0*CS%Kvbbl / (hvel(i,nz)*H_to_m + 2.0e-10*dt*CS%Kvbbl)
    endif
  endif ; enddo

  if (associated(visc%Kv_turb)) then
     ! BGR/ Add factor of 2. * the averaged Kv_turb.
     !      this is needed to reproduce the analytical solution to
     !      a simple diffusion problem, likely due to h_shear being
     !      equal to 2 x \delta z
    if (work_on_u) then
      do K=nz,2,-1 ; do i=is,ie ; if (do_i(i)) then
        a(i,K) = a(i,K) + (2.*0.5)*(visc%Kv_turb(i,j,k) + visc%Kv_turb(i+1,j,k))
      endif ; enddo ; enddo
    else
      do K=nz,2,-1 ; do i=is,ie ; if (do_i(i)) then
        a(i,K) = a(i,K) + (2.*0.5)*(visc%Kv_turb(i,j,k) + visc%Kv_turb(i,j+1,k))
      endif ; enddo ; enddo
    endif
  endif

  do K=nz,2,-1 ; do i=is,ie ; if (do_i(i)) then
    !    botfn determines when a point is within the influence of the bottom
    !  boundary layer, going from 1 at the bottom to 0 in the interior.
    z2 = z_i(i,k)
    botfn = 1.0 / (1.0 + 0.09*z2*z2*z2*z2*z2*z2)

    if (CS%bottomdraglaw) then
      a(i,K) = a(i,K) + 2.0*(kv_bbl(i)-CS%Kv)*botfn
      r = (hvel(i,k)+hvel(i,k-1))
      if (r > 2.0*bbl_thick(i)) then
        h_shear = ((1.0 - botfn) * r + botfn*2.0*bbl_thick(i))
      else
        h_shear = r
      endif
    else
      a(i,K) = a(i,K) + 2.0*(CS%Kvbbl-CS%Kv)*botfn
      h_shear = hvel(i,k) + hvel(i,k-1) + h_neglect
    endif

    !   Up to this point a has units of m2 s-1, but now is converted to m s-1.
    !   The term including 1e-10 in the denominators is here to avoid
    ! truncation error problems in the tridiagonal solver. Effectively, this
    ! sets the maximum coupling coefficient at 1e10 m.
    a(i,K) = a(i,K) / (h_shear*H_to_m + 1.0e-10*dt*a(i,K))
  endif ; enddo ; enddo ! i & k loops

  if (do_shelf) then
    ! Set the coefficients to include the no-slip surface stress.
    do i=is,ie ; if (do_i(i)) then
      if (work_on_u) then
        kv_tbl(i) = visc%kv_tbl_shelf_u(I,j)
        tbl_thick(i) = visc%tbl_thick_shelf_u(I,j) * m_to_H
      else
        kv_tbl(i) = visc%kv_tbl_shelf_v(i,J)
        tbl_thick(i) = visc%tbl_thick_shelf_v(i,J) * m_to_H
      endif
      z_t(i) = 0.0

      ! If a(i,1) were not already 0, it would be added here.
      if (0.5*hvel(i,1) > tbl_thick(i)) then
        a(i,1) = kv_tbl(i) / (tbl_thick(i) *H_to_m + (1.0e-10*dt)*kv_tbl(i))
      else
        a(i,1) = kv_tbl(i) / (0.5*hvel(i,1)*H_to_m + (1.0e-10*dt)*kv_tbl(i))
      endif
    endif ; enddo

    do K=2,nz ; do i=is,ie ;  if (do_i(i)) then
      z_t(i) = z_t(i) + hvel(i,k-1) / tbl_thick(i)
      topfn = 1.0 / (1.0 + 0.09 * z_t(i)**6)

      r = (hvel(i,k)+hvel(i,k-1))
      if (r > 2.0*tbl_thick(i)) then
        h_shear = ((1.0 - topfn) * r + topfn*2.0*tbl_thick(i))
      else
        h_shear = r
      endif
    !   The term including 1e-10 in the denominators is here to avoid
    ! truncation error problems in the tridiagonal solver. Effectively, this
    ! sets the maximum coupling coefficient increment to 1e10 m.
      a_top = 2.0 * topfn * kv_tbl(i)
      a(i,K) = a(i,K) + a_top / (h_shear*H_to_m + 1.0e-10*dt*a_top)
    endif ; enddo ; enddo
  elseif (CS%dynamic_viscous_ML .or. (GV%nkml>0)) then
    max_nk = 0
    do i=is,ie ; if (do_i(i)) then
      if (GV%nkml>0) nk_visc(i) = real(GV%nkml+1)
      if (work_on_u) then
        u_star(I) = 0.5*(fluxes%ustar(i,j) + fluxes%ustar(i+1,j))
        absf(I) = 0.5*(abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I,J)))
        if (CS%dynamic_viscous_ML) nk_visc(I) = visc%nkml_visc_u(I,j) + 1
      else
        u_star(i) = 0.5*(fluxes%ustar(i,j) + fluxes%ustar(i,j+1))
        absf(i) = 0.5*(abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))
        if (CS%dynamic_viscous_ML) nk_visc(i) = visc%nkml_visc_v(i,J) + 1
      endif
      h_ml(i) = h_neglect ; z_t(i) = 0.0
      max_nk = max(max_nk,ceiling(nk_visc(i) - 1.0))
    endif ; enddo

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
      temp1 = (z_t(i)*h_ml(i) - z_t(i)*z_t(i)) * H_to_m
      !   This viscosity is set to go to 0 at the mixed layer top and bottom
      ! (in a log-layer) and be further limited by rotation to give the
      ! natural Ekman length.
      visc_ml = u_star(i) * 0.41 * (temp1*u_star(i)) / &
                     (absf(i)*temp1 + h_ml(i)*u_star(i))
      a_ml = 4.0*visc_ml / ((hvel(i,k)+hvel(i,k-1) + h_neglect) * H_to_m + &
                            2.0e-10*dt*visc_ml)
      ! Choose the largest estimate of a.
      if (a_ml > a(i,K)) a(i,K) = a_ml
    endif ; endif ; enddo ; enddo
  endif

end subroutine find_coupling_coef

!> Velocity components which exceed a threshold for physically
!! reasonable values are truncated. Optionally, any column with excessive
!! velocities may be sent to a diagnostic reporting subroutine.
subroutine vertvisc_limit_vel(u, v, h, ADp, CDp, fluxes, visc, dt, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  real, intent(inout), &
    dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: u      !< Zonal velocity in m s-1
  real, intent(inout), &
    dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: v      !< Meridional velocity in m s-1
  real, intent(in), &
    dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: h      !< Layer thickness in H
  type(accel_diag_ptrs), intent(in)      :: ADp    !< Acceleration diagnostic pointers
  type(cont_diag_ptrs),  intent(in)      :: CDp    !< Continuity diagnostic pointers
  type(forcing),         intent(in)      :: fluxes !< Forcing fields
  type(vertvisc_type),   intent(in)      :: visc   !< Viscosities and bottom drag
  real,                  intent(in)      :: dt     !< Time increment in s
  type(vertvisc_CS),     pointer         :: CS     !< Vertical viscosity control structure

  ! Local variables

  real :: maxvel           ! Velocities components greater than maxvel
  real :: truncvel         ! are truncated to truncvel, both in m s-1.
  real :: CFL              ! The local CFL number.
  real :: H_report         ! A thickness below which not to report truncations.
  real :: dt_Rho0          ! The timestep divided by the Boussinesq density, in dt m3 kg-1.
  real :: vel_report(SZIB_(G),SZJB_(G))
  real :: u_old(SZIB_(G),SZJ_(G),SZK_(G))
  real :: v_old(SZI_(G),SZJB_(G),SZK_(G))
  logical :: trunc_any, dowrite(SZIB_(G),SZJB_(G))
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  maxvel = CS%maxvel
  truncvel = 0.9*maxvel
  H_report = 6.0 * GV%Angstrom
  dt_Rho0 = dt / GV%Rho0

  if (len_trim(CS%u_trunc_file) > 0) then
!$OMP parallel do default(none) shared(js,je,Isq,Ieq,nz,CS,G,fluxes,u,h,dt,maxvel,ADp,CDp,truncvel, &
!$OMP                                  u_old,vel_report,dowrite,H_report) &
!$OMP                          private(trunc_any,CFL)
    do j=js,je
      trunc_any = .false.
      do I=Isq,Ieq ; dowrite(I,j) = .false. ; enddo
      if (CS%CFL_based_trunc) then
        do I=Isq,Ieq ; vel_report(i,j) = 3.0e8 ; enddo ! Speed of light default.
        do k=1,nz ; do I=Isq,Ieq
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
        do k=1,nz ; do I=Isq,Ieq ; if (abs(u(I,j,k)) > maxvel) then
          dowrite(I,j) = .true. ; trunc_any = .true.
        endif ;enddo ; enddo
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
        endif ; enddo ;  enddo
      endif ; endif
    enddo ! j-loop
  else
    if (CS%CFL_based_trunc) then
!$OMP parallel do default(none) shared(nz,js,je,Isq,Ieq,u,dt,G,CS,h,H_report)
      do k=1,nz ; do j=js,je ; do I=Isq,Ieq
        if ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i+1,j) < -CS%CFL_trunc) then
          u(I,j,k) = (-0.9*CS%CFL_trunc) * (G%areaT(i+1,j) / (dt * G%dy_Cu(I,j)))
          if (h(i,j,k) + h(i+1,j,k) > H_report) CS%ntrunc = CS%ntrunc + 1
        elseif ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i,j) > CS%CFL_trunc) then
          u(I,j,k) = (0.9*CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dy_Cu(I,j)))
          if (h(i,j,k) + h(i+1,j,k) > H_report) CS%ntrunc = CS%ntrunc + 1
        endif
      enddo ; enddo ; enddo
    else
!$OMP parallel do default(none) shared(nz,js,je,Isq,Ieq,u,G,CS,truncvel,maxvel,h,H_report)
      do k=1,nz ; do j=js,je ; do I=Isq,Ieq ; if (abs(u(I,j,k)) > maxvel) then
        u(I,j,k) = SIGN(truncvel,u(I,j,k))
        if (h(i,j,k) + h(i+1,j,k) > H_report) CS%ntrunc = CS%ntrunc + 1
      endif ; enddo ; enddo ; enddo
    endif
  endif

  if (len_trim(CS%u_trunc_file) > 0) then
    do j=js,je; do I=Isq,Ieq ; if (dowrite(I,j)) then
!   Here the diagnostic reporting subroutines are called if
! unphysically large values were found.
        call write_u_accel(I, j, u_old, h, ADp, CDp, dt, G, GV, CS%PointAccel_CSp, &
               vel_report(I,j), -vel_report(I,j), fluxes%taux(I,j)*dt_Rho0, &
               a=CS%a_u(:,j,:), hv=CS%h_u(:,j,:))
    endif ; enddo; enddo
  endif


  if (len_trim(CS%v_trunc_file) > 0) then
!$OMP parallel do default(none) shared(Jsq,Jeq,is,ie,nz,CS,G,fluxes,v,h,dt,maxvel,ADp,CDp,truncvel, &
!$OMP                                  v_old,vel_report,dowrite,H_report)                           &
!$OMP                          private(trunc_any,CFL)
    do J=Jsq,Jeq
      trunc_any = .false.
      do i=is,ie ; dowrite(i,J) = .false. ; enddo
      if (CS%CFL_based_trunc) then
        do i=is,ie ; vel_report(i,J) = 3.0e8 ; enddo ! Speed of light default.
        do k=1,nz ; do i=is,ie
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
        do k=1,nz ; do i=is,ie ; if (abs(v(i,J,k)) > maxvel) then
          dowrite(i,J) = .true. ; trunc_any = .true.
        endif ; enddo ; enddo
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
        endif ; enddo ;  enddo
      endif ; endif
    enddo ! J-loop
  else
    if (CS%CFL_based_trunc) then
!$OMP parallel do default(none) shared(is,ie,Jsq,Jeq,nz,v,dt,G,CS,h,H_report)
      do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
        if ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j+1) < -CS%CFL_trunc) then
          v(i,J,k) = (-0.9*CS%CFL_trunc) * (G%areaT(i,j+1) / (dt * G%dx_Cv(i,J)))
          if (h(i,j,k) + h(i,j+1,k) > H_report) CS%ntrunc = CS%ntrunc + 1
        elseif ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j) > CS%CFL_trunc) then
          v(i,J,k) = (0.9*CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dx_Cv(i,J)))
          if (h(i,j,k) + h(i,j+1,k) > H_report) CS%ntrunc = CS%ntrunc + 1
        endif
      enddo ; enddo ; enddo
    else
!$OMP parallel do default(none) shared(is,ie,Jsq,Jeq,nz,v,G,CS,h,truncvel,maxvel,H_report)
      do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie ; if (abs(v(i,J,k)) > maxvel) then
        v(i,J,k) = SIGN(truncvel,v(i,J,k))
        if (h(i,j,k) + h(i,j+1,k) > H_report) CS%ntrunc = CS%ntrunc + 1
      endif ; enddo ; enddo ; enddo
    endif
  endif

  if (len_trim(CS%v_trunc_file) > 0) then
    do J=Jsq,Jeq; do i=is,ie ; if (dowrite(i,J)) then
!   Here the diagnostic reporting subroutines are called if
! unphysically large values were found.
        call write_v_accel(i, J, v_old, h, ADp, CDp, dt, G, GV, CS%PointAccel_CSp, &
               vel_report(i,J), -vel_report(i,J), fluxes%tauy(i,J)*dt_Rho0, &
               a=CS%a_v(:,J,:),hv=CS%h_v(:,J,:))
    endif ; enddo; enddo
  endif

end subroutine vertvisc_limit_vel

!> Initialise the vertical friction module
subroutine vertvisc_init(MIS, Time, G, GV, param_file, diag, ADp, dirs, &
                         ntrunc, CS)
  !> "MOM Internal State", a set of pointers to the fields and accelerations
  !! that make up the ocean's physical state
  type(ocean_internal_state), target, intent(in) :: MIS
  type(time_type), target, intent(in)    :: Time       !< Current model time
  type(ocean_grid_type),   intent(in)    :: G          !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV         !< Ocean vertical grid structure
  type(param_file_type),   intent(in)    :: param_file !< File to parse for parameters
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostic control structure
  type(accel_diag_ptrs),   intent(inout) :: ADp        !< Acceleration diagnostic pointers
  type(directories),       intent(in)    :: dirs       !< Relevant directory paths
  integer, target,         intent(inout) :: ntrunc     !< Number of velocity truncations
  type(vertvisc_CS),       pointer       :: CS         !< Vertical viscosity control structure

  ! Local variables

  real :: hmix_str_dflt
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_vert_friction" ! This module's name.
  character(len=40)  :: thickness_units = "meters or kg m-2"

  if (associated(CS)) then
    call MOM_error(WARNING, "vertvisc_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%diag => diag ; CS%ntrunc => ntrunc ; ntrunc = 0

! Default, read and log parameters
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "BOTTOMDRAGLAW", CS%bottomdraglaw, &
                 "If true, the bottom stress is calculated with a drag \n"//&
                 "law of the form c_drag*|u|*u. The velocity magnitude \n"//&
                 "may be an assumed value or it may be based on the \n"//&
                 "actual velocity in the bottommost HBBL, depending on \n"//&
                 "LINEAR_DRAG.", default=.true.)
  call get_param(param_file, mod, "CHANNEL_DRAG", CS%Channel_drag, &
                 "If true, the bottom drag is exerted directly on each \n"//&
                 "layer proportional to the fraction of the bottom it \n"//&
                 "overlies.", default=.false.)
  call get_param(param_file, mod, "DIRECT_STRESS", CS%direct_stress, &
                 "If true, the wind stress is distributed over the \n"//&
                 "topmost HMIX_STRESS of fluid (like in HYCOM), and KVML \n"//&
                 "may be set to a very small value.", default=.false.)
  call get_param(param_file, mod, "DYNAMIC_VISCOUS_ML", CS%dynamic_viscous_ML, &
                 "If true, use a bulk Richardson number criterion to \n"//&
                 "determine the mixed layer thickness for viscosity.", &
                 default=.false.)
  call get_param(param_file, mod, "U_TRUNC_FILE", CS%u_trunc_file, &
                 "The absolute path to a file into which the accelerations \n"//&
                 "leading to zonal velocity truncations are written.  \n"//&
                 "Undefine this for efficiency if this diagnostic is not \n"//&
                 "needed.", default=" ")
  call get_param(param_file, mod, "V_TRUNC_FILE", CS%v_trunc_file, &
                 "The absolute path to a file into which the accelerations \n"//&
                 "leading to meridional velocity truncations are written. \n"//&
                 "Undefine this for efficiency if this diagnostic is not \n"//&
                 "needed.", default=" ")
  call get_param(param_file, mod, "HARMONIC_VISC", CS%harmonic_visc, &
                 "If true, use the harmonic mean thicknesses for \n"//&
                 "calculating the vertical viscosity.", default=.false.)
  call get_param(param_file, mod, "HARMONIC_BL_SCALE", CS%harm_BL_val, &
                 "A scale to determine when water is in the boundary \n"//&
                 "layers based solely on harmonic mean thicknesses for \n"//&
                 "the purpose of determining the extent to which the \n"//&
                 "thicknesses used in the viscosities are upwinded.", &
                 default=0.0, units="nondim")
  call get_param(param_file, mod, "DEBUG", CS%debug, default=.false.)

  if (GV%nkml < 1) &
    call get_param(param_file, mod, "HMIX_FIXED", CS%Hmix, &
                 "The prescribed depth over which the near-surface \n"//&
                 "viscosity and diffusivity are elevated when the bulk \n"//&
                 "mixed layer is not used.", units="m", fail_if_missing=.true.)
  if (CS%direct_stress) then
    if (GV%nkml < 1) then
      call get_param(param_file, mod, "HMIX_STRESS", CS%Hmix_stress, &
                 "The depth over which the wind stress is applied if \n"//&
                 "DIRECT_STRESS is true.", units="m", default=CS%Hmix)
    else
      call get_param(param_file, mod, "HMIX_STRESS", CS%Hmix_stress, &
                 "The depth over which the wind stress is applied if \n"//&
                 "DIRECT_STRESS is true.", units="m", fail_if_missing=.true.)
    endif
    if (CS%Hmix_stress <= 0.0) call MOM_error(FATAL, "vertvisc_init: " // &
       "HMIX_STRESS must be set to a positive value if DIRECT_STRESS is true.")
  endif
  call get_param(param_file, mod, "KV", CS%Kv, &
                 "The background kinematic viscosity in the interior. \n"//&
                 "The molecular value, ~1e-6 m2 s-1, may be used.", &
                 units="m2 s-1", fail_if_missing=.true.)

! CS%Kvml = CS%Kv ; CS%Kvbbl = CS%Kv ! Needed? -AJA
  if (GV%nkml < 1) call get_param(param_file, mod, "KVML", CS%Kvml, &
                 "The kinematic viscosity in the mixed layer.  A typical \n"//&
                 "value is ~1e-2 m2 s-1. KVML is not used if \n"//&
                 "BULKMIXEDLAYER is true.  The default is set by KV.", &
                 units="m2 s-1", default=CS%Kv)
  if (.not.CS%bottomdraglaw) call get_param(param_file, mod, "KVBBL", CS%Kvbbl, &
                 "The kinematic viscosity in the benthic boundary layer. \n"//&
                 "A typical value is ~1e-2 m2 s-1. KVBBL is not used if \n"//&
                 "BOTTOMDRAGLAW is true.  The default is set by KV.", &
                 units="m2 s-1", default=CS%Kv)
  call get_param(param_file, mod, "HBBL", CS%Hbbl, &
                 "The thickness of a bottom boundary layer with a \n"//&
                 "viscosity of KVBBL if BOTTOMDRAGLAW is not defined, or \n"//&
                 "the thickness over which near-bottom velocities are \n"//&
                 "averaged for the drag law if BOTTOMDRAGLAW is defined \n"//&
                 "but LINEAR_DRAG is not.", units="m", fail_if_missing=.true.)
  call get_param(param_file, mod, "MAXVEL", CS%maxvel, &
                 "The maximum velocity allowed before the velocity \n"//&
                 "components are truncated.", units="m s-1", default=3.0e8)
  call get_param(param_file, mod, "CFL_BASED_TRUNCATIONS", CS%CFL_based_trunc, &
                 "If true, base truncations on the CFL number, and not an \n"//&
                 "absolute speed.", default=.true.)
  call get_param(param_file, mod, "CFL_TRUNCATE", CS%CFL_trunc, &
                 "The value of the CFL number that will cause velocity \n"//&
                 "components to be truncated; instability can occur past 0.5.", &
                 units="nondim", default=0.5)
  call get_param(param_file, mod, "CFL_REPORT", CS%CFL_report, &
                 "The value of the CFL number that causes accelerations \n"//&
                 "to be reported; the default is CFL_TRUNCATE.", &
                 units="nondim", default=CS%CFL_trunc)
  call get_param(param_file, mod, "CFL_TRUNCATE_RAMP_TIME", CS%truncRampTime, &
                 "The time over which the CFL trunction value is ramped\n"//&
                 "up at the beginning of the run.", &
                 units="s", default=0.)
  CS%CFL_truncE = CS%CFL_trunc
  call get_param(param_file, mod, "CFL_TRUNCATE_START", CS%CFL_truncS, &
                 "The start value of the truncation CFL number used when\n"//&
                 "ramping up CFL_TRUNC.", &
                 units="nondim", default=0.)

  ALLOC_(CS%a_u(IsdB:IedB,jsd:jed,nz+1)) ; CS%a_u(:,:,:) = 0.0
  ALLOC_(CS%h_u(IsdB:IedB,jsd:jed,nz))   ; CS%h_u(:,:,:) = 0.0
  ALLOC_(CS%a_v(isd:ied,JsdB:JedB,nz+1)) ; CS%a_v(:,:,:) = 0.0
  ALLOC_(CS%h_v(isd:ied,JsdB:JedB,nz))   ; CS%h_v(:,:,:) = 0.0

  CS%id_au_vv = register_diag_field('ocean_model', 'au_visc', diag%axesCui, Time, &
     'Zonal Viscous Vertical Coupling Coefficient', 'meter second-1')
  CS%id_av_vv = register_diag_field('ocean_model', 'av_visc', diag%axesCvi, Time, &
     'Meridional Viscous Vertical Coupling Coefficient', 'meter second-1')

  CS%id_h_u = register_diag_field('ocean_model', 'Hu_visc', diag%axesCuL, Time, &
     'Thickness at Zonal Velocity Points for Viscosity', thickness_units)
  CS%id_h_v = register_diag_field('ocean_model', 'Hv_visc', diag%axesCvL, Time, &
     'Thickness at Meridional Velocity Points for Viscosity', thickness_units)
  CS%id_hML_u = register_diag_field('ocean_model', 'HMLu_visc', diag%axesCu1, Time, &
     'Mixed Layer Thickness at Zonal Velocity Points for Viscosity', thickness_units)
  CS%id_hML_v = register_diag_field('ocean_model', 'HMLv_visc', diag%axesCv1, Time, &
     'Mixed Layer Thickness at Meridional Velocity Points for Viscosity', thickness_units)

  CS%id_du_dt_visc = register_diag_field('ocean_model', 'du_dt_visc', diag%axesCuL, &
     Time, 'Zonal Acceleration from Vertical Viscosity', 'meter second-2')
  if (CS%id_du_dt_visc > 0) call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
  CS%id_dv_dt_visc = register_diag_field('ocean_model', 'dv_dt_visc', diag%axesCvL, &
     Time, 'Meridional Acceleration from Vertical Viscosity', 'meter second-2')
  if (CS%id_dv_dt_visc > 0) call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)

  CS%id_taux_bot = register_diag_field('ocean_model', 'taux_bot', diag%axesCu1, &
     Time, 'Zonal Bottom Stress from Ocean to Earth', 'Pa')
  CS%id_tauy_bot = register_diag_field('ocean_model', 'tauy_bot', diag%axesCv1, &
     Time, 'Meridional Bottom Stress from Ocean to Earth', 'Pa')

  if ((len_trim(CS%u_trunc_file) > 0) .or. (len_trim(CS%v_trunc_file) > 0)) &
    call PointAccel_init(MIS, Time, G, param_file, diag, dirs, CS%PointAccel_CSp)

end subroutine vertvisc_init

!> Update the CFL truncation value as a function of time.
!! If called with the optional argument activate=.true., record the
!! value of Time as the beginning of the ramp period.
subroutine updateCFLtruncationValue(Time, CS, activate)
  type(time_type), target, intent(in)    :: Time     !< Current model time
  type(vertvisc_CS),       pointer       :: CS       !< Vertical viscosity control structure
  !> Whether to record the value of Time as the beginning of the ramp period
  logical, optional,       intent(in)    :: activate

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
  type(vertvisc_CS),   pointer       :: CS
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
